/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "libgretl.h"
#include "libset.h"
#include "missing_private.h"
#include "gretl_bfgs.h"
#include "gretl_normal.h"
#include "qr_estimate.h"
#include "matrix_extra.h"
#include "gretl_string_table.h"

#include <errno.h>

/**
 * SECTION:discrete
 * @short_description: models for limited dependent variables
 * and related cases.
 * @title: Limdep
 * @include: libgretl.h
 *
 * Covers logit (binary, ordered or multinomial), probit (binary
 * or ordered), logistic, tobit, interval regression, models for count data
 * and for duration data, and the heckit sample-selection model.
 * Plus a few utility functions.
 */

#define LPDEBUG 0

typedef struct op_container_ op_container;

/* structure for handling ordered probit or logit */

struct op_container_ {
    int ci;           /* model command index (PROBIT or LOGIT) */
    gretlopt opt;     /* option flags */
    int bootstrap;    /* state: doing bootstrap of normality test */
    int *y;           /* dependent variable */
    double **Z;       /* data */
    int *list;        /* dependent var plus regular regressors */
    int ymin;         /* minimum of original y values */
    int ymax;         /* max of (possibly normalized) y */
    int t1;           /* beginning of sample */
    int t2;           /* end of sample */
    int nobs;         /* number of observations */
    int nx;           /* number of explanatory variables */
    int k;            /* total number of parameters */
    double *theta;    /* real parameter estimates */
    double *ndx;      /* index variable */
    double *dP;       /* probabilities */
    MODEL *pmod;      /* model struct, initially containing OLS */
    gretl_matrix *G;  /* score matrix by observation */
    double *g;        /* total score vector */
    gretl_matrix *nty; /* dependent var for normality test */
    gretl_matrix *ntb; /* coefficients, normality test */
    double X20;       /* original value of normality test */
    int replics;      /* replications of normtest */
    int X2_ngt;       /* times X20 exceeded in bootstrap */
};

struct sorter {
    double x;
    int t;
};

static double lp_cdf (double x, int ci)
{
    switch (ci) {
    case PROBIT:
        return normal_cdf(x);
    case LOGIT:
        return 1.0 / (1.0 + exp(-x));
    default:
        return NADBL;
    }
}

static double lp_pdf (double x, int ci)
{
    double tmp;

    switch (ci) {
    case PROBIT:
        return normal_pdf(x);
    case LOGIT:
        tmp = 1.0 + exp(-x);
        return (tmp - 1.0) / (tmp * tmp);
    default:
        return NADBL;
    }
}

static void op_container_destroy (op_container *OC)
{
    free(OC->y);
    free(OC->ndx);
    free(OC->dP);
    free(OC->list);
    free(OC->g);
    free(OC->theta);

    gretl_matrix_free(OC->G);
    gretl_matrix_free(OC->nty);
    gretl_matrix_free(OC->ntb);

    free(OC);
}

static op_container *op_container_new (int ci, int ndum, int ymin,
                                       double **Z, MODEL *pmod,
                                       gretlopt opt)
{
    op_container *OC;
    int i, t, vy = pmod->list[1];
    int nobs = pmod->nobs;
    int err = 0;

    OC = malloc(sizeof *OC);
    if (OC == NULL) {
        return NULL;
    }

    OC->ci = ci;
    OC->opt = opt;
    OC->bootstrap = 0;

    OC->Z = Z;
    OC->pmod = pmod;
    OC->t1 = pmod->t1;
    OC->t2 = pmod->t2;
    OC->nobs = nobs;
    OC->k = pmod->ncoeff;
    OC->ymin = ymin;
    OC->ymax = ndum;
    OC->nx = OC->k - ndum;

    OC->y = NULL;
    OC->ndx = NULL;
    OC->dP = NULL;
    OC->list = NULL;
    OC->g = NULL;
    OC->theta = NULL;
    OC->G = NULL;
    OC->nty = NULL;
    OC->ntb = NULL;

    OC->y = malloc(nobs * sizeof *OC->y);
    OC->ndx = malloc(nobs * sizeof *OC->ndx);
    OC->dP = malloc(nobs * sizeof *OC->dP);

    OC->list = gretl_list_new(1 + OC->nx);
    OC->g = malloc(OC->k * sizeof *OC->g);
    OC->theta = malloc(OC->k * sizeof *OC->theta);

    if (OC->y == NULL || OC->ndx == NULL ||
        OC->dP == NULL || OC->list == NULL ||
        OC->g == NULL || OC->theta == NULL) {
        op_container_destroy(OC);
        return NULL;
    }

    if (ci == PROBIT) {
        /* include extra storage for normality test */
        OC->G = gretl_matrix_alloc(nobs, OC->k + 2);
        OC->nty = gretl_matrix_alloc(nobs, 1);
        OC->ntb = gretl_matrix_alloc(OC->k + 2, 1);
        if (OC->G == NULL || OC->nty == NULL || OC->ntb == NULL) {
            err = E_ALLOC;
        } else {
            /* "shrink" G to proper size for gradient */
            gretl_matrix_reuse(OC->G, -1, OC->k);
        }
    } else {
        OC->G = gretl_matrix_alloc(nobs, OC->k);
        if (OC->G == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        op_container_destroy(OC);
        return NULL;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (!na(pmod->uhat[t])) {
            OC->y[i++] = (int) Z[vy][t];
        }
    }

    OC->list[1] = vy;
    for (i=0; i<OC->nx; i++) {
        OC->list[i+2] = pmod->list[i+2];
    }

    /* for probit normality test bootstrap */
    OC->X20 = NADBL;
    OC->replics = 0;
    OC->X2_ngt = 0;

#if LPDEBUG
    fprintf(stderr, "nobs = %d\n", OC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", OC->t1, OC->t2);
    fprintf(stderr, "k = %d\n", OC->k);
    fprintf(stderr, "ndum = %d\n", ndum);
    fprintf(stderr, "nx = %d\n", OC->nx);
    fprintf(stderr, "Max(y) = M = %d\n", OC->ymax);
    printlist(OC->list, "list, in op_container_new");
#endif

    return OC;
}

static int op_compute_score (op_container *OC, int yt,
                             double ystar0, double ystar1,
                             double dP, int t, int s)
{
    double gsi, dm, mills0, mills1;
    int M = OC->ymax;
    int i, v;

    if (ystar1 < 6.0 || OC->ci == LOGIT) {
        mills0 = (yt == 0)? 0.0 : lp_pdf(ystar0, OC->ci) / dP;
        mills1 = (yt == M)? 0.0 : lp_pdf(ystar1, OC->ci) / dP;
    } else {
        /* L'Hopital-based approximation */
        mills0 = (yt == 0)? 0.0 : -ystar0;
        mills1 = (yt == M)? 0.0 : -ystar1;
    }

    dm = mills1 - mills0;

    for (i=0; i<OC->nx; i++) {
        v = OC->list[i+2];
        gsi = -dm * OC->Z[v][t];
        gretl_matrix_set(OC->G, s, i, gsi);
        OC->g[i] += gsi;
    }

    for (i=OC->nx; i<OC->k; i++) {
        gretl_matrix_set(OC->G, s, i, 0.0);
        if (i == OC->nx + yt - 1) {
            gsi = -mills0;
            gretl_matrix_set(OC->G, s, i, gsi);
            OC->g[i] += gsi;
        }
        if (i == OC->nx + yt) {
            gsi = mills1;
            gretl_matrix_set(OC->G, s, i, gsi);
            OC->g[i] += gsi;
        }
    }

    return 0;
}

#define dPMIN 1.0e-15

static int op_compute_probs (const double *theta, op_container *OC)
{
    double m0, m1, ystar0 = 0.0, ystar1 = 0.0;
    int M = OC->ymax;
    int nx = OC->nx;
    double P0, P1, h, adj, dP;
    int i, t, s, yt;

    /* initialize analytical score */
    for (i=0; i<OC->k; i++) {
        OC->g[i] = 0.0;
    }

    s = 0;
    for (t=OC->pmod->t1; t<=OC->pmod->t2; t++) {
        if (na(OC->pmod->uhat[t])) {
#if LPDEBUG > 1
            fprintf(stderr, "obs %4d excluded\n", t);
#endif
            continue;
        }
        yt = OC->y[s];
        if (yt == 0) {
            m0 = theta[nx];
            ystar1 = OC->ndx[s] + m0;
        } else {
            m0 = theta[nx + yt - 1];
            ystar0 = OC->ndx[s] + m0;
            if (yt < M) {
                m1 = theta[nx + yt];
                ystar1 = OC->ndx[s] + m1;
            }
        }
#if LPDEBUG > 1
        fprintf(stderr, "t:%4d/%d s=%d y=%d, ndx = %10.6f, ystar0 = %9.7f, ystar1 = %9.7f\n",
                t+1, OC->nobs, s, yt, OC->ndx[s], ystar0, ystar1);
#endif
        if (ystar0 < 6.0 || OC->ci == LOGIT) {
            P0 = (yt == 0)? 0.0 : lp_cdf(ystar0, OC->ci);
            P1 = (yt == M)? 1.0 : lp_cdf(ystar1, OC->ci);
            dP = P1 - P0;
        } else {
            /* Taylor-based 1st order approximation */
            h = ystar1 - ystar0;
            adj = lp_pdf(ystar1, OC->ci) + lp_pdf(ystar0, OC->ci);
            dP =  0.5 * h * adj;
        }
        if (dP > dPMIN) {
            OC->dP[s] = dP;
        } else {
#if LPDEBUG
            fprintf(stderr, "very small dP at obs %d; y=%d, ndx=%g, dP=%.15g\n",
                    t, yt, OC->ndx[s], dP);
#endif
            return 1;
        }
        op_compute_score(OC, yt, ystar0, ystar1, dP, t, s);
        s++;
    }

    return 0;
}

/* Below: method for getting around the "non-increasing cut point"
   issue in ordered models by construction: the 2nd and higher cut
   points are represented to the optimizer in the form of the
   log-difference from the previous cut point.
*/

static void op_transform_theta (op_container *OC, double *theta)
{
    int i;

    for (i=0; i<=OC->nx; i++) {
        theta[i] = OC->theta[i];
    }

    for (i=OC->nx+1; i<OC->k; i++) {
        /* convert cut point 2 and higher to log-difference form */
        theta[i] = log(OC->theta[i] - OC->theta[i-1]);
    }
}

/* Inverse operation for the transformation done by
   op_transform_theta() */

static void op_get_real_theta (op_container *OC, const double *theta)
{
    int i;

    for (i=0; i<=OC->nx; i++) {
        OC->theta[i] = theta[i];
    }

    for (i=OC->nx+1; i<OC->k; i++) {
        /* retrieve cut point 2 and higher from log-difference form */
        OC->theta[i] = exp(theta[i]) + OC->theta[i-1];
    }
}

static double op_loglik (const double *theta, void *ptr)
{
    op_container *OC = (op_container *) ptr;
    double x, ll = 0.0;
    int i, s, t, v;
    int err;

    if (theta != OC->theta) {
        op_get_real_theta(OC, theta);
    }

    s = 0;
    for (t=OC->t1; t<=OC->t2; t++) {
        if (na(OC->pmod->uhat[t])) {
            continue;
        }
        x = 0.0;
        for (i=0; i<OC->nx; i++) {
            /* the independent variables */
            v = OC->list[i+2];
            x -= OC->theta[i] * OC->Z[v][t];
        }
        OC->ndx[s++] = x;
#if LPDEBUG > 2
        fprintf(stderr, "t = %d, s = %d, x = %g\n", t, s, x);
#endif
    }

    err = op_compute_probs(OC->theta, OC);
    if (err) {
        ll = NADBL;
    } else {
        s = 0;
        for (t=OC->t1; t<=OC->t2; t++) {
            if (!na(OC->pmod->uhat[t])) {
                ll += log(OC->dP[s++]);
            }
        }
    }

#if LPDEBUG > 1
    fprintf(stderr, "ll = %16.10f\n", ll);
#endif

    return ll;
}

static int op_score (double *theta, double *s, int npar, BFGS_CRIT_FUNC ll,
                     void *ptr)
{
    op_container *OC = (op_container *) ptr;
    int i, j;

    for (i=0; i<npar; i++) {
        s[i] = OC->g[i];
    }

    for (i=OC->nx; i<npar; i++) {
        for (j=i+1; j<npar; j++) {
            /* add effects of changes in subsequent cut points */
            s[i] += OC->g[j];
        }
        if (i > OC->nx) {
            s[i] *= exp(theta[i]); /* apply jacobian */
        }
    }

    return 0;
}

static int ordered_hessian (op_container *OC, gretl_matrix *H)
{
    double smal = 1.0e-09;  /* "small" is some sort of macro on win32 */
    double dx, dx2;
    double ti, x, ll, *g0;
    int i, j, k = OC->k;
    int err = 0;

    g0 = malloc(k * sizeof *g0);
    if (g0 == NULL) {
        return E_ALLOC;
    }

    for (i=0; i<k; i++) {
        ti = OC->theta[i];
        dx = (fabs(ti) > 0.001) ? fabs(ti) * smal : smal;
        dx2 = 2.0 * dx;
        OC->theta[i] -= dx;
        ll = op_loglik(OC->theta, OC);
        if (na(ll)) {
            OC->theta[i] = ti;
            err = E_DATA;
            break;
        }
        for (j=0; j<k; j++) {
            g0[j] = OC->g[j];
        }
        OC->theta[i] += dx2;
        ll = op_loglik(OC->theta, OC);
        if (na(ll)) {
            OC->theta[i] = ti;
            err = E_DATA;
            break;
        }
        for (j=0; j<k; j++) {
            x = (OC->g[j] - g0[j]) / dx2;
            gretl_matrix_set(H, i, j, -x);
        }
        /* restore original theta */
        OC->theta[i] = ti;
    }

    free(g0);

    if (!err) {
        gretl_matrix_xtr_symmetric(H);
    }

    return err;
}

static gretl_matrix *ordered_hessian_inverse (op_container *OC,
                                              int *err)
{
    gretl_matrix *H = gretl_zero_matrix_new(OC->k, OC->k);

    if (H == NULL) {
        *err = E_ALLOC;
    } else {
        *err = ordered_hessian(OC, H);
    }

    if (!*err) {
        *err = gretl_invert_symmetric_matrix(H);
        if (*err) {
            fprintf(stderr, "ordered_hessian_inverse: inversion failed\n");
        }
    }

    if (H != NULL && *err) {
        gretl_matrix_free(H);
        H = NULL;
    }

    return H;
}

/**
 * ordered_model_prediction:
 * @pmod: model for ordered data, either logit or probit.
 * @Xb: X\beta, the value of the index function at a given
 * observation.
 * @ymin: the minimum value of the dependent variable.
 *
 * Returns: the "predicted value" of the (ordinal) dependent variable,
 * taken to be the value for which the estimated probability is greatest.
 */

double ordered_model_prediction (const MODEL *pmod, double Xb,
                                 int ymin)
{
    /* position of least cut point in coeff array */
    int k = gretl_model_get_int(pmod, "nx");
    int maxval = pmod->ncoeff - k;
    double prob, pmax, cut;
    double CDF, CDFbak;
    int i, pred = ymin;

    cut = pmod->coeff[k];
    pmax = CDFbak = lp_cdf(cut - Xb, pmod->ci);

    for (i=1; i<maxval; i++) {
        cut = pmod->coeff[++k];
        CDF = lp_cdf(cut - Xb, pmod->ci);
        prob = CDF - CDFbak;
        if (prob > pmax) {
            pmax = prob;
            pred = ymin + i;
        }
        CDFbak = CDF;
    }

    prob = 1 - CDFbak;
    if (prob > pmax) {
        pred = ymin + maxval;
    }

    return (double) pred;
}

gretl_matrix *ordered_probabilities (const MODEL *pmod,
                                     const double *zhat,
                                     int t1, int t2,
                                     const DATASET *dset,
                                     int *err)
{
    gretl_matrix *P;
    int k = gretl_model_get_int(pmod, "nx");
    const double *c = pmod->coeff + k;
    int ncut = pmod->ncoeff - k;
    int n = t2 - t1 + 1;
    int ci = pmod->ci;
    char **S = NULL;
    double zht, pij;
    int i, t, j;

    P = gretl_matrix_alloc(n, ncut+1);
    if (P == NULL) {
        *err = E_ALLOC;
        return NULL;
    }
    S = strings_array_new(n);

    for (t=t1, i=0; t<=t2; t++, i++) {
        zht = zhat[t];
        if (na(zht)) {
            for (j=0; j<=ncut; j++) {
                gretl_matrix_set(P, i, j, NADBL);
            }
        } else {
            pij = lp_cdf(c[0] - zht, ci);
            gretl_matrix_set(P, i, 0, pij);
            for (j=1; j<ncut; j++) {
                pij = lp_cdf(c[j] - zht, ci) - lp_cdf(c[j-1] - zht, ci);
                gretl_matrix_set(P, i, j, pij);
            }
            pij = 1.0 - lp_cdf(c[ncut-1] - zht, pmod->ci);
            gretl_matrix_set(P, i, ncut, pij);
        }
        if (S != NULL) {
            S[i] = retrieve_date_string(t+1, dset, err);
        }
    }

    gretl_matrix_set_t1(P, t1);
    gretl_matrix_set_t2(P, t2);
    if (S != NULL) {
        gretl_matrix_set_rownames(P, S);
    }

    return P;
}

/**
 * mn_logit_prediction:
 * @Xt: vector of regressors at observation t.
 * @b: array of coefficients.
 * @yvals: vector of dependent variable values.
  *
 * Returns: the predicted value of the dependent variable, that
 * is, the value for which the estimated probability is greatest.
 */

double mn_logit_prediction (const gretl_matrix *Xt,
                            const double *b,
                            const gretl_matrix *yvals)
{
    double *eXtb = NULL;
    double St, pj, pmax;
    int i, j, k, m, nx;
    int pidx = 0;

    nx = gretl_vector_get_length(Xt);
    m = gretl_vector_get_length(yvals);

    eXtb = malloc(m * sizeof *eXtb);
    if (eXtb == NULL) {
        return NADBL;
    }

    /* base case */
    eXtb[0] = St = 1.0;
    k = 0;

    /* loop across the other y-values */
    for (j=1; j<m; j++) {
        /* accumulate exp(X*beta) */
        eXtb[j] = 0.0;
        for (i=0; i<nx; i++) {
            eXtb[j] += Xt->val[i] * b[k++];
        }
        eXtb[j] = exp(eXtb[j]);
        St += eXtb[j];
    }

    pmax = 0.0;

    for (j=0; j<m; j++) {
        pj = eXtb[j] / St;
        if (pj > pmax) {
            pmax = pj;
            pidx = j;
        }
    }

    free(eXtb);

    return yvals->val[pidx];
}

/* compute generalized residual for ordered models */

static double op_gen_resid (op_container *OC, const double *theta, int t)
{
    double ndxt, m0, m1, ystar0, f0, f1;
    double ret, dP, ystar1 = 0.0;
    int M = OC->ymax;
    int nx = OC->nx;
    int yt;

    dP = OC->dP[t];
    yt = OC->y[t];
    ndxt = OC->ndx[t];

    if (yt == 0) {
        m0 = theta[nx];
        ystar1 = ndxt + m0;
    } else {
        m0 = theta[nx + yt - 1];
        ystar0 = ndxt + m0;
        if (yt < M) {
            m1 = theta[nx + yt];
            ystar1 = ndxt + m1;
        }
    }

    if (ystar1 < 6.0 || OC->ci == LOGIT || 1) {
        f0 = (yt == 0)? 0.0 : lp_pdf(ystar0, OC->ci) / dP;
        f1 = (yt == M)? 0.0 : lp_pdf(ystar1, OC->ci) / dP;
    } else {
        /* L'Hôpital-based approximation */
        f0 = (yt == 0)? 0.0 : -ystar0;
        f1 = (yt == M)? 0.0 : -ystar1;
    }

    ret = (f0 - f1);

    return ret;
}

/* Initialize the cut-points by counting the occurrences of each value
   of the (normalized) dependent variable, finding the sample
   proportion (cumulating as we go), and taking the inverse of the
   normal CDF.
*/

static void cut_points_init (op_container *OC,
                             const MODEL *pmod,
                             const double **Z)
{
    const double *y = Z[pmod->list[1]];
    double p = 0.0;
    int i, j, t, nj;

    for (i=OC->nx, j=0; i<OC->k; i++, j++) {
        nj = 0;
        for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t]) && y[t] == j) {
                nj++;
            }
        }
        p += (double) nj / pmod->nobs;
        OC->theta[i] = normal_cdf_inverse(p);
    }
}

static void add_pseudo_rsquared (MODEL *pmod, double L0,
                                 int k, int T, gretlopt opt)
{
    if (opt & OPT_S) {
        /* Estrella pseudo-R^2 */
        double expon = -2.0 * L0/T;

        pmod->rsq = 1.0 - pow(pmod->lnL/L0, expon);
        pmod->adjrsq = 1.0 - pow((pmod->lnL - k)/L0, expon);
        pmod->opt |= OPT_S;
    } else {
        /* McFadden pseudo-R^2 */
        pmod->rsq = 1.0 - pmod->lnL/L0;
        pmod->adjrsq = 1.0 - (pmod->lnL - k)/L0;
    }
}

static double binary_null_loglik (int *y, int T)
{
    double L0 = 0;
    int ones = 0;
    int zeros, t;

    for (t=0; t<T; t++) {
        ones += y[t];
    }
    zeros = T - ones;

    L0 = ones * log(ones / (double) T);
    L0 += zeros * log(zeros / (double) T);

    return L0;
}

static void op_LR_test (MODEL *pmod, op_container *OC,
                        const double **Z)
{
    int full_nx = OC->nx;
    int restore = 1;
    double L0;

    if (OC->ymax == 1) {
        restore = 0;
        L0 = binary_null_loglik(OC->y, OC->nobs);
    } else {
        OC->k -= OC->nx;
        OC->nx = 0;
        cut_points_init(OC, pmod, Z);
        L0 = op_loglik(OC->theta, OC);
    }

    if (!na(L0) && L0 <= pmod->lnL) {
        pmod->chisq = 2.0 * (pmod->lnL - L0);
        gretl_model_set_int(pmod, "lr_df", full_nx);
        add_pseudo_rsquared(pmod, L0, full_nx, OC->nobs, OC->opt);
    }

    if (restore) {
        /* restore original data on OC */
        OC->nx = full_nx;
        OC->k += OC->nx;
    }
}

static int oprobit_normtest (MODEL *pmod,
                             op_container *OC)
{
    gretl_matrix *ntX;
    double *theta = OC->theta;
    int k = OC->k;
    int nx = OC->nx;
    int t, s, yt;
    double u, v, a2v, b2u;
    double a = 0, b = 0;
    double e3, e4;
    int err = 0;

    /* augmented version of G matrix */
    ntX = gretl_matrix_reuse(OC->G, -1, k+2);

    s = 0;
    for (t=OC->pmod->t1; t<=OC->pmod->t2; t++) {
        if (na(OC->pmod->uhat[t])) {
            continue;
        }
        yt = OC->y[s];
        if (yt == 0) {
            b = OC->ndx[s] + theta[nx];
            u = gretl_matrix_get(OC->G, s, nx);
            b2u = b*b*u;
            v = a2v = 0;
        } else {
            a = OC->ndx[s] + theta[nx + yt - 1];
            v = - gretl_matrix_get(OC->G, s, nx + yt - 1);
            a2v = a*a*v;
            if (yt < OC->ymax) {
                b = OC->ndx[s] + theta[nx + yt];
                u = gretl_matrix_get(OC->G, s, nx + yt);
                b2u = b*b*u;
            } else {
                u = b2u = 0;
            }
        }
        e3 = 2*(v-u) + (a2v - b2u);
        e4 = 3*(a*v-b*u) + (a*a2v - b*b2u);
        gretl_matrix_set(ntX, s, k, e3);
        gretl_matrix_set(ntX, s, k+1, e4);
        s++;
    }

    /* dependent var should be all 1s */
    for (t=0; t<OC->nobs; t++) {
        OC->nty->val[t] = 1.0;
    }

    err = gretl_matrix_ols(OC->nty, ntX, OC->ntb, NULL, NULL, NULL);

    if (!err) {
        double X2 = OC->nobs;

        gretl_matrix_multiply(ntX, OC->ntb, OC->nty);
        for (t=0; t<OC->nobs; t++) {
            u = 1 - OC->nty->val[t];
            X2 -= u * u;
        }
#if 0
        fprintf(stderr, "normtest: X2 = %g\n", X2);
#endif
        if (X2 > 0) {
            if (OC->bootstrap) {
                OC->replics += 1;
                if (X2 > OC->X20) {
                    OC->X2_ngt += 1;
                }
            } else {
                ModelTest *test;

                OC->X20 = X2;
                gretl_model_add_normality_test(pmod, X2);
                test = gretl_model_get_test(pmod, GRETL_TEST_NORMAL);
                if (test != NULL) {
                    /* note asymptotic nature of test */
                    model_test_set_opt(test, OPT_A);
                }
            }
        }
    } else {
        fprintf(stderr, "oprobit_normtest: err = %d\n", err);
    }

    /* return G to correct size for gradient */
    gretl_matrix_reuse(OC->G, -1, k);

    return err;
}

static int fill_op_model (MODEL *pmod, const int *list,
                          const DATASET *dset,
                          op_container *OC,
                          int fncount, int grcount)
{
    gretl_matrix *H = NULL;
    int npar = OC->k;
    int nx = OC->nx;
    int correct = 0;
    double xti, Xb;
    int i, s, t, v, omp;
    int err = 0;

    H = ordered_hessian_inverse(OC, &err);
    if (err) {
        goto bailout;
    }

    if (OC->opt & OPT_R) {
        err = gretl_model_add_QML_vcv(pmod, OC->ci, H, OC->G,
                                      dset, OC->opt, NULL);
    } else {
        err = gretl_model_add_hessian_vcv(pmod, H);
    }

    gretl_matrix_free(H);

    if (err) {
        goto bailout;
    }

    pmod->ci = OC->ci;
    gretl_model_set_int(pmod, "ordered", 1);
    gretl_model_set_int(pmod, "nx", OC->nx);
    gretl_model_set_int(pmod, "ymin", OC->ymin);

    /* blank out invalid statistics */
    pmod->rsq = pmod->adjrsq = pmod->fstt = pmod->sigma = NADBL;
    gretl_model_destroy_data_item(pmod, "centered-R2");
    gretl_model_destroy_data_item(pmod, "uncentered");

    if (grcount > 0) {
        gretl_model_set_int(pmod, "fncount", fncount);
        gretl_model_set_int(pmod, "grcount", grcount);
    } else {
        gretl_model_set_int(pmod, "iters", fncount);
    }

    pmod->ncoeff = npar;
    for (i=0; i<npar; i++) {
        pmod->coeff[i] = OC->theta[i];
    }

    if (OC->ci == PROBIT) {
        oprobit_normtest(pmod, OC);
    }

    s = 0;
    for (t=OC->t1; t<=OC->t2; t++) {
        Xb = 0.0;
        for (i=0; i<OC->nx; i++) {
            v = OC->list[i+2];
            xti = OC->Z[v][t];
            if (na(xti)) {
                Xb = NADBL;
                break;
            } else {
                Xb += OC->theta[i] * xti;
            }
        }
        /* yhat = X\hat{beta} */
        pmod->yhat[t] = Xb;
        if (na(Xb) || na(pmod->uhat[t])) {
            continue;
        }
        omp = (int) ordered_model_prediction(pmod, Xb, 0);
        if (omp == OC->y[s]) {
            correct++;
        }
        /* compute generalized residual */
        pmod->uhat[t] = op_gen_resid(OC, OC->theta, s);
        s++;
    }

    gretl_model_set_int(pmod, "correct", correct);

    pmod->lnL = op_loglik(OC->theta, OC);
    mle_criteria(pmod, 0);
    pmod->rsq = pmod->adjrsq = NADBL;

    gretl_model_allocate_param_names(pmod, npar);

    if (pmod->errcode == 0) {
        char tmp[16];

        for (i=0; i<nx; i++) {
            v = OC->list[i+2];
            gretl_model_set_param_name(pmod, i, dset->varname[v]);
        }
        s = 1;
        for (i=nx; i<npar; i++) {
            sprintf(tmp, "cut%d", s++);
            gretl_model_set_param_name(pmod, i, tmp);
        }
    }

    /* trim the model list: remove references to the 'cut'
       dummy variables */
    for (i=pmod->list[0]; i>1; i--) {
        if (!in_gretl_list(list, pmod->list[i])) {
            gretl_list_delete_at_pos(pmod->list, i);
        }
    }

    if (nx > 0) {
        op_LR_test(pmod, OC, (const double **) dset->Z);
        gretl_model_set_coeff_separator(pmod, NULL, nx);
    }

 bailout:

    if (err && !pmod->errcode) {
        pmod->errcode = err;
    }

    return pmod->errcode;
}

static void record_bootstrap_pvalue (op_container *OC,
                                     MODEL *pmod)
{
    double pval = OC->X2_ngt / (double) OC->replics;
    ModelTest *test;

    test = gretl_model_get_test(pmod, GRETL_TEST_NORMAL);
    if (test != NULL) {
        model_test_set_pvalue(test, pval);
        model_test_set_opt(test, OPT_B);
    }
}

/* Prepare for a bootstrap iteration of the ordered probit
   normality test: create artificial y. Note: we're using
   pmod->coeff in creating y, thereby ensuring that we get
   the coefficients from the original estimation, which
   are saved onto @pmod before bootstrapping starts.
*/

static void op_boot_prep (op_container *OC,
                          const MODEL *pmod)
{
    double ystar, *cut = pmod->coeff + OC->nx;
    int y, ncut = OC->k - OC->nx;
    int i, v, t, s = 0;

    for (t=OC->t1; t<=OC->t2; t++) {
        if (na(OC->pmod->uhat[t])) {
            continue;
        }
        ystar = gretl_one_snormal();
        /* add regression effect */
        for (i=0; i<OC->nx; i++) {
            v = OC->list[i+2];
            ystar += pmod->coeff[i] * OC->Z[v][t];
        }
        y = 0;
        /* convert to observable using cut points */
        for (i=0; i<ncut; i++) {
            if (ystar > cut[i]) {
                y++;
            } else {
                break;
            }
        }
        OC->y[s++] = y;
    }
}

static void op_boot_init (op_container *OC,
                          int *bs_maxit)
{
    int K, err = 0;

    K = get_optval_int(PROBIT, OPT_B, &err);
    if (!err && K > 0) {
        *bs_maxit = K;
    }
    OC->bootstrap = 1;
}

/* Main ordered estimation function */

static int do_ordered (int ci, int ndum, int ymin,
                       DATASET *dset, MODEL *pmod,
                       const int *list,
                       gretlopt opt, PRN *prn)
{
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;
    op_container *OC;
    int i, npar;
    double *theta = NULL;
    double toler;
    int bs_iter = 0;
    int bs_maxit = 1000;
    int use_newton = 0;
    gretlopt maxopt;
    int err;

    OC = op_container_new(ci, ndum, ymin, dset->Z, pmod, opt);
    if (OC == NULL) {
        return E_ALLOC;
    }

    npar = OC->k;
    /* transformed theta to pass to optimizer */
    theta = malloc(npar * sizeof *theta);
    if (theta == NULL) {
        op_container_destroy(OC);
        return E_ALLOC;
    }

    if (libset_get_int(GRETL_OPTIM) == OPTIM_NEWTON) {
        use_newton = 1;
    }

    /* initialize slopes */
    for (i=0; i<OC->nx; i++) {
        OC->theta[i] = 0.0001;
    }

    /* initialize cut points */
    cut_points_init(OC, pmod, (const double **) dset->Z);

    /* transform theta to log-diff form */
    op_transform_theta(OC, theta);

#if LPDEBUG
    for (i=0; i<npar; i++) {
        fprintf(stderr, "theta[%d]: 'real' = %g, transformed = %g\n", i,
                OC->theta[i], theta[i]);
    }
    fprintf(stderr, "\ninitial loglikelihood = %.12g\n",
            op_loglik(theta, OC));
#endif

 reestimate:

    if (OC->bootstrap) {
        /* prepare for bootstrap iteration */
        fncount = grcount = 0;
        op_boot_prep(OC, pmod);
        bs_iter++;
    }

    maxopt = (prn != NULL)? (OPT_U | OPT_V) : OPT_U;

    if (use_newton) {
        double crittol = 1.0e-7;
        double gradtol = 1.0e-7;

        err = newton_raphson_max(theta, npar, maxit,
                                 crittol, gradtol, &fncount,
                                 C_LOGLIK, op_loglik,
                                 op_score, NULL, OC,
                                 maxopt, prn);
        fprintf(stderr, "use_newton: err = %d\n", err);
    } else {
        BFGS_defaults(&maxit, &toler, PROBIT);
        err = BFGS_max(theta, npar, maxit, toler,
                       &fncount, &grcount, op_loglik, C_LOGLIK,
                       op_score, OC, NULL, maxopt, prn);
    }

    if (!err && !OC->bootstrap) {
        /* transform back to 'real' theta and fill the model struct */
        op_get_real_theta(OC, theta);
        err = fill_op_model(pmod, list, dset, OC, fncount, grcount);
    }

    if ((err == E_NOCONV || err == E_NAN) && bs_iter > 0) {
        /* tolerate random numerical problems? */
        err = 0;
        bs_iter--;
        goto reestimate;
    }

    if (!err && (opt & OPT_B)) {
        /* bootstrapping the ordered probit normality test */
        if (bs_iter == 0) {
            /* start the procedure */
            op_boot_init(OC, &bs_maxit);
        } else {
            /* bootstrap in progress: run test */
            oprobit_normtest(NULL, OC);
        }
        if (bs_iter == bs_maxit) {
            record_bootstrap_pvalue(OC, pmod);
        } else {
            goto reestimate;
        }
    }

    free(theta);
    op_container_destroy(OC);

    return err;
}

/* We want to ensure that the values of the dependent variable
   actually used in the analysis (after dropping any bad
   observations) form a zero-based series of consecutive
   integers.
*/

static int maybe_fix_op_depvar (MODEL *pmod, DATASET *dset,
                                double **orig_y, int *ndum,
                                int *ymin)
{
    gretl_matrix *v = NULL;
    double *yvals = NULL;
    int dv = pmod->list[1];
    int i, t, n = 0;
    int fixit = 0;
    int nv = 0;
    int err = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (!na(pmod->uhat[t])) {
            n++;
        }
    }

    /* Transcribe the y values that were used in
       the initial OLS
    */

    yvals = malloc(n * sizeof *yvals);
    if (yvals == NULL) {
        return E_ALLOC;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (!na(pmod->uhat[t])) {
            yvals[i++] = dset->Z[dv][t];
        }
    }

    /* Make a sorted vector containing the distinct
       values of y
    */
    v = gretl_matrix_values(yvals, n, OPT_S, &err);

#if LPDEBUG
    gretl_matrix_print(v, "distinct y values");
#endif

    if (!err) {
        nv = gretl_vector_get_length(v);
        *ndum = nv - 1;
        if (v->val[0] != 0.0) {
            /* the minimum y-value is not zero */
            fixit = 1;
        } else {
            for (i=1; i<nv; i++) {
                if (v->val[i] != v->val[i-1] + 1) {
                    /* y values are not consecutive integers */
                    fixit = 1;
                    break;
                }
            }
        }
    }

    if (fixit) {
        double *normy = malloc(dset->n * sizeof *normy);

        if (normy == NULL) {
            err = E_ALLOC;
        } else {
            for (t=0; t<dset->n; t++) {
                normy[t] = NADBL;
                if (!na(pmod->uhat[t])) {
                    for (i=0; i<nv; i++) {
                        if (dset->Z[dv][t] == v->val[i]) {
                            normy[t] = i;
                            break;
                        }
                    }
                }
            }
            /* Back up the original y and replace it for
               the duration of the ordered analysis.
            */
            *orig_y = dset->Z[dv];
            dset->Z[dv] = normy;
            *ymin = v->val[0];
        }
    }

#if LPDEBUG
    if (!err && fixit) {
        fputs("ordered model: using normalized y\n", stderr);
    } else {
        fputs("ordered model: using original y\n", stderr);
    }
#endif

    free(yvals);
    gretl_matrix_free(v);

    return err;
}

static void restore_depvar (double **Z, double *y, int v)
{
    free(Z[v]);
    Z[v] = y;
}

static int *make_dummies_list (const int *list,
                               DATASET *dset,
                               int *err)
{
    int *dumlist = gretl_list_new(1);

    if (dumlist == NULL) {
        *err = E_ALLOC;
    } else {
        dumlist[1] = list[1];

        /* OPT_F -> drop first value */
        *err = list_dumgenr(&dumlist, dset, OPT_F);
        if (*err) {
            free(dumlist);
            dumlist = NULL;
        }
    }

    return dumlist;
}

/* make internal regression list for ordered model */

static int *make_op_list (const int *list, DATASET *dset,
                          int **pdumlist, int *err)
{
    int *dumlist;
    int *biglist;
    int i, k, nv;

    dumlist = make_dummies_list(list, dset, err);
    if (dumlist == NULL) {
        return NULL;
    }

    nv = list[0] + dumlist[0];

    biglist = gretl_list_new(nv);
    if (biglist == NULL) {
        free(dumlist);
        *err = E_ALLOC;
        return NULL;
    }

    k = 1;
    for (i=1; i<=list[0]; i++) {
        biglist[k++] = list[i];
    }
    for (i=1; i<=dumlist[0]; i++) {
        biglist[k++] = dumlist[i];
    }

    *pdumlist = dumlist;

    return biglist;
}

static int list_purge_const (int *list, DATASET *dset)
{
    MODEL tmpmod;
    int depvar = list[1];
    int i, j, n, ok = 0;
    int err = 0;

    /* first remove "const" (var 0) itself, if present */
    for (i=2; i<=list[0]; i++) {
        if (list[i] == 0) {
            gretl_list_delete_at_pos(list, i);
            break;
        }
    }

    /* drop other stuff possibly collinear with the constant
       (e.g. sets of dummies) */

    list[1] = 0;     /* substitute the constant as dependent */
    n = list[0] - 1; /* number of RHS terms */

    for (j=0; j<n && !ok; j++) {
        int vi, pos;

        tmpmod = lsq(list, dset, OLS, OPT_A);
        if (tmpmod.errcode) {
            err = tmpmod.errcode;
            break;
        }
        ok = (tmpmod.ess > 1.0e-6);
        if (!ok) {
            for (i=tmpmod.ncoeff-1; i>=0; i--) {
                if (fabs(tmpmod.coeff[i]) > 1.0e-06) {
                    /* tmpmod.list and list may not be identical */
                    vi = tmpmod.list[i+2];
                    pos = in_gretl_list(list, vi);
                    if (pos >= 2) {
                        gretl_list_delete_at_pos(list, pos);
                    }
                    break;
                }
            }
        }
        clear_model(&tmpmod);
    }

    /* reinstate the real dependent variable */
    list[1] = depvar;

    return err;
}

static int ordered_depvar_check (int v, const DATASET *dset)
{
    if (!series_is_discrete(dset, v) &&
        !gretl_is_oprobit_ok(dset->t1, dset->t2, dset->Z[v])) {
        gretl_errmsg_sprintf(_("The variable '%s' is not discrete"),
                             dset->varname[v]);
        return E_DATA;
    }

    return 0;
}

/* driver function for ordered logit/probit */

MODEL ordered_estimate (int ci, const int *list,
                        DATASET *dset, gretlopt opt,
                        PRN *prn)
{
    MODEL model;
    PRN *vprn;
    int orig_v = dset->v;
    double *orig_y = NULL;
    int *mylist = NULL;
    int *biglist = NULL;
    int *dumlist = NULL;
    int ymin = 0;
    int ndum = 0;

    vprn = (opt & OPT_V)? prn : NULL;

    mylist = gretl_list_copy(list);
    gretl_model_init(&model, dset);
    model.errcode = ordered_depvar_check(mylist[1], dset);

    if (!model.errcode) {
        /* remove the constant from the incoming list, if present */
        model.errcode = list_purge_const(mylist, dset);
    }

    if (!model.errcode) {
        /* construct augmented regression list, including dummies
           for the level of the dependent variable
        */
        biglist = make_op_list(mylist, dset, &dumlist, &model.errcode);
    }

    if (!model.errcode) {
        /* run initial OLS, with dummies added */
        model = lsq(biglist, dset, OLS, OPT_A);
        if (model.errcode) {
            fprintf(stderr, "ordered_estimate: initial OLS failed\n");
        }
    }

    if (model.errcode) {
        free(mylist);
        free(dumlist);
        free(biglist);
        return model;
    }

#if LPDEBUG
    PRN *dprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

    pputs(dprn, "ordered_estimate: initial OLS\n");
    printmodel(&model, dset, OPT_S, dprn);
    gretl_print_destroy(dprn);
#endif

    if (!model.errcode) {
        /* after accounting for any missing observations, normalize
           the dependent variable if necessary
        */
        model.errcode = maybe_fix_op_depvar(&model, dset, &orig_y,
                                            &ndum, &ymin);
    }

    /* do the actual ordered probit analysis */
    if (!model.errcode) {
        clear_model_xpx(&model);
        model.errcode = do_ordered(ci, ndum, ymin, dset, &model,
                                   list, opt, vprn);
    }

    free(dumlist);
    free(biglist);

    if (orig_y != NULL) {
        /* if we messed with the dependent var, put the original back */
        restore_depvar(dset->Z, orig_y, mylist[1]);
    }

    if (dset->v > orig_v) {
        /* clean up any automatically-added dummies */
        dataset_drop_last_variables(dset, dset->v - orig_v);
    }

    set_model_id(&model, opt);
    free(mylist);

    return model;
}

static double logit (double x)
{
    double l = 1.0 / (1.0 + exp(-x));

#if LPDEBUG
    if (x > 40 || x < -40) {
        fprintf(stderr, "x = %g, logit(x) = %.16f\n", x, l);
    }
#endif

    return l;
}

static double logit_pdf (double x)
{
    double l, z = exp(-x);

    l = z / ((1.0 + z) * (1.0 + z));

#if LPDEBUG
    if (x > 40 || x < -40) {
        fprintf(stderr, "x = %g, logit_pdf(x) = %g\n", x, l);
    }
#endif

    if (x < 0 && isnan(l)) {
#if LPDEBUG
        fprintf(stderr, "logit_pdf(): x = %g, forcing l to zero\n", x);
#endif
        l = 0;
    }

    return l;
}

/* Here we're checking for a dummy variable that acts as a
   "one-way perfect predictor" of the binary dependent variable;
   that is, an x such that Prob(y = A | x = B) = 1 for some
   assignment of values 0 or 1 to A and B. The MLE does not
   exist in the presence of such a regressor, so we'll remove
   it from the model (after alerting the user).

   The approach taken by Stata is not only to drop such a
   regressor but also to drop the observations that are thus
   perfectly predicted. As the outcome of disussions in
   November-December 2014 we decided not to follow Stata in
   this policy: we just drop the regressor. However, in case
   we want to revisit this point I'm leaving in place the
   apparatus required to implement the Stata policy.

   AC, 2014-12-04
*/

#define LIKE_STATA 0

static char *classifier_check (int *list, const DATASET *dset,
                               PRN *prn, int *ndropped, int *err)
{
    char *mask = NULL;
    int yno = list[1];
    const double *y = dset->Z[yno];
    int i, v, ni, t;

    *ndropped = 0;

    for (i=list[0]; i>=2; i--) {
        int getout = 0;

        v = list[i];
        if (v == 0) {
            continue;
        }

        ni = gretl_isdummy(dset->t1, dset->t2, dset->Z[v]);

        if (ni > 0) {
            const double *x = dset->Z[v];
            int xytab[4] = {0};
            int pp0 = 1, pp1 = 1;
            int maskval = -1;

            for (t=dset->t1; t<=dset->t2; t++) {
                xytab[0] += (x[t] == 0 && y[t] == 0);
                xytab[1] += (x[t] == 0 && y[t] == 1);
                xytab[2] += (x[t] == 1 && y[t] == 0);
                xytab[3] += (x[t] == 1 && y[t] == 1);
                if (xytab[1] && xytab[3]) {
                    /* x does not perfectly predict y == 0 */
                    pp0 = 0;
                }
                if (xytab[0] && xytab[2]) {
                    /* x does not perfectly predict y == 1 */
                    pp1 = 0;
                }
                if (!pp0 && !pp1) {
                    break;
                }
            }

#if LIKE_STATA
            if (pp0 && pp1) {
                pputc(prn, '\n');
                pprintf(prn, "Note: %s = %s%s at all observations\n",
                        dset->varname[yno], xytab[0] ? "" : "not-",
                        dset->varname[v]);
                *err = E_NOCONV;
                getout = 1;
            } else if (pp0) {
                maskval = xytab[1] ? 1 : 0;
                pputc(prn, '\n');
                pprintf(prn, "Note: Prob(%s = %d | %s = %d) = 1\n",
                        dset->varname[yno], 0, dset->varname[v],
                        maskval);
            } else if (pp1) {
                pputc(prn, '\n');
                maskval = xytab[0] ? 1 : 0;
                pprintf(prn, "Note: Prob(%s = %d | %s = %d) = 1\n",
                        dset->varname[yno], 1, dset->varname[v],
                        maskval);
            }

            if (maskval >= 0) {
                if (mask == NULL) {
                    mask = malloc(dset->n + 1);
                    if (mask == NULL) {
                        *err = E_ALLOC;
                    } else {
                        memset(mask, '0', dset->n);
                        mask[dset->n] = 0;
                    }
                }
                if (mask != NULL) {
                    for (t=dset->t1; t<=dset->t2; t++) {
                        if (dset->Z[v][t] == maskval) {
                            mask[t-dset->t1] = '1';
                            *ndropped += 1;
                        }
                    }
                    pprintf(prn, "%s dropped and %d observations not used\n",
                            dset->varname[v], *ndropped);
                }

                gretl_list_delete_at_pos(list, i);
                /* It'll get too confusing if we try doing
                   this for more than one regressor?
                */
                getout = 1;
            }
#else /* not like Stata */
            if (pp0 && pp1) {
                pputc(prn, '\n');
                pprintf(prn, _("Note: %s = %s%s at all observations\n"),
                        dset->varname[yno], xytab[0] ? "" : "!",
                        dset->varname[v]);
                *err = E_NOCONV;
                getout = 1;
            } else if (pp0) {
                maskval = xytab[1] ? 1 : 0;
                pputc(prn, '\n');
                pprintf(prn, _("Note: Prob(%s = %d | %s = %d) = 1\n"),
                        dset->varname[yno], 0, dset->varname[v],
                        maskval);
            } else if (pp1) {
                maskval = xytab[0] ? 1 : 0;
                pputc(prn, '\n');
                pprintf(prn, _("Note: Prob(%s = %d | %s = %d) = 1\n"),
                        dset->varname[yno], 1, dset->varname[v],
                        maskval);
            }

            if (maskval >= 0) {
                pprintf(prn, _("Dropping %s\n"), dset->varname[v]);
                gretl_list_delete_at_pos(list, i);
            }
#endif
        }
        if (getout) {
            break;
        }
    }

#if LPDEBUG
    fprintf(stderr, "classifier check: mask = %p\n", (void *) mask);
#endif

    return mask;
}

/* struct for holding multinomial logit info */

typedef struct mnl_info_ mnl_info;

struct mnl_info_ {
    int n;            /* number of categories (excluding base) */
    int k;            /* number of coeffs per category */
    int npar;         /* total number of parameters */
    int T;            /* number of observations */
    double *theta;    /* coeffs for Newton/BFGS */
    gretl_matrix_block *B;
    gretl_matrix *y;  /* dependent variable */
    gretl_matrix *X;  /* regressors */
    gretl_matrix *b;  /* coefficients, matrix form */
    gretl_matrix *Xb; /* coeffs times regressors */
    gretl_matrix *P;  /* probabilities */
};

static void mnl_info_destroy (mnl_info *mnl)
{
    if (mnl != NULL) {
        gretl_matrix_block_destroy(mnl->B);
        free(mnl->theta);
        free(mnl);
    }
}

static mnl_info *mnl_info_new (int n, int k, int T)
{
    mnl_info *mnl = malloc(sizeof *mnl);
    int i;

    if (mnl != NULL) {
        mnl->n = n;
        mnl->k = k;
        mnl->T = T;
        mnl->npar = k * n;
        mnl->theta = malloc(mnl->npar * sizeof *mnl->theta);
        if (mnl->theta == NULL) {
            free(mnl);
            return NULL;
        }
        mnl->B = gretl_matrix_block_new(&mnl->y, T, 1,
                                        &mnl->X, T, k,
                                        &mnl->b, k, n,
                                        &mnl->Xb, T, n,
                                        &mnl->P, T, n,
                                        NULL);
        if (mnl->B == NULL) {
            free(mnl->theta);
            free(mnl);
            mnl = NULL;
        } else {
            for (i=0; i<mnl->npar; i++) {
                mnl->theta[i] = 0.0;
            }
        }
    }

    return mnl;
}

/* compute loglikelihood for multinomial logit */

static double mn_logit_loglik (const double *theta, void *ptr)
{
    mnl_info *mnl = (mnl_info *) ptr;
    double x, xti, exti, ll = 0.0;
    int yt, i, t;

    for (i=0; i<mnl->npar; i++) {
        mnl->b->val[i] = theta[i];
    }

    gretl_matrix_multiply(mnl->X, mnl->b, mnl->Xb);

    for (t=0; t<mnl->T; t++) {
        x = 1.0;
        for (i=0; i<mnl->n; i++) {
            /* sum row i of exp(Xb) */
            xti = gretl_matrix_get(mnl->Xb, t, i);
            exti = exp(xti);
            x += exti;
        }
        ll -= log(x);
        yt = gretl_vector_get(mnl->y, t);
        if (yt > 0) {
            ll += gretl_matrix_get(mnl->Xb, t, yt-1);
        }
    }

    return ll;
}

static int mn_logit_score (double *theta, double *s, int npar,
                           BFGS_CRIT_FUNC ll, void *ptr)
{
    mnl_info *mnl = (mnl_info *) ptr;
    double x, xti, exti, pti, p, g;
    int i, j, k, t, yt;
    int err = 0;

    for (i=0; i<npar; i++) {
        s[i] = 0.0;
    }

    for (t=0; t<mnl->T && !errno; t++) {
        x = 1.0;
        for (i=0; i<mnl->n; i++) {
            /* sum row i of exp(Xb) */
            xti = gretl_matrix_get(mnl->Xb, t, i);
            exti = exp(xti);
            x += exti;
            gretl_matrix_set(mnl->P, t, i, exti);
        }
        yt = gretl_vector_get(mnl->y, t);
        k = 0;
        for (i=0; i<mnl->n; i++) {
            pti = gretl_matrix_get(mnl->P, t, i) / x;
            gretl_matrix_set(mnl->P, t, i, pti);
            p = (i == (yt-1)) - pti;
            for (j=0; j<mnl->k; j++) {
                g = p * gretl_matrix_get(mnl->X, t, j);
                s[k++] += g;
            }
        }
    }

    return err;
}

/* multinomial logit: form the negative of the analytical
   Hessian */

static int mnl_hessian (double *theta, gretl_matrix *H, void *data)
{
    mnl_info *mnl = data;
    gretl_matrix_block *B;
    gretl_matrix *x;
    gretl_matrix *xx;
    gretl_matrix *hjk;
    double xti, ptj, ptk;
    int r, c;
    int i, j, k, t;

    B = gretl_matrix_block_new(&x, 1, mnl->k,
                               &xx, mnl->k, mnl->k,
                               &hjk, mnl->k, mnl->k,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    r = c = 0;

    for (j=0; j<mnl->n; j++) {
        for (k=0; k<=j; k++) {
            gretl_matrix_zero(hjk);
            for (t=0; t<mnl->T; t++) {
                for (i=0; i<mnl->k; i++) {
                    xti = gretl_matrix_get(mnl->X, t, i);
                    gretl_vector_set(x, i, xti);
                }
                gretl_matrix_multiply_mod(x, GRETL_MOD_TRANSPOSE,
                                          x, GRETL_MOD_NONE,
                                          xx, GRETL_MOD_NONE);
                ptj = gretl_matrix_get(mnl->P, t, j);
                ptk = gretl_matrix_get(mnl->P, t, k);
                gretl_matrix_multiply_by_scalar(xx, ptj * ((j == k) - ptk));
                gretl_matrix_add_to(hjk, xx);
            }
            gretl_matrix_inscribe_matrix(H, hjk, r, c, GRETL_MOD_NONE);
            if (j != k) {
                gretl_matrix_inscribe_matrix(H, hjk, c, r, GRETL_MOD_NONE);
            }
            c += mnl->k;
        }
        r += mnl->k;
        c = 0;
    }

    gretl_matrix_block_destroy(B);

    return 0;
}

static gretl_matrix *mnl_hessian_inverse (mnl_info *mnl, int *err)
{
    gretl_matrix *H;

    H = gretl_zero_matrix_new(mnl->npar, mnl->npar);
    if (H == NULL) {
        *err = E_ALLOC;
    } else {
        *err = mnl_hessian(mnl->theta, H, mnl);
    }

    if (!*err) {
        *err = gretl_invert_symmetric_matrix(H);
    }

    return H;
}

static gretl_matrix *mnl_score_matrix (mnl_info *mnl, int *err)
{
    gretl_matrix *G;
    double p, g;
    int yt, i, j, k, t;

    G = gretl_matrix_alloc(mnl->T, mnl->npar);
    if (G == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (t=0; t<mnl->T; t++) {
        yt = gretl_vector_get(mnl->y, t);
        k = 0;
        for (i=0; i<mnl->n; i++) {
            p = (i == (yt-1)) - gretl_matrix_get(mnl->P, t, i);
            for (j=0; j<mnl->k; j++) {
                g = p * gretl_matrix_get(mnl->X, t, j);
                gretl_matrix_set(G, t, k++, g);
            }
        }
    }

    return G;
}

static int mnl_add_variance_matrix (MODEL *pmod, mnl_info *mnl,
                                    const DATASET *dset,
                                    gretlopt opt)
{
    gretl_matrix *H = NULL;
    gretl_matrix *G = NULL;
    int err = 0;

    H = mnl_hessian_inverse(mnl, &err);
    if (err) {
        return err;
    }

    if (opt & OPT_R) {
        G = mnl_score_matrix(mnl, &err);
    }

    if (!err) {
        if (opt & OPT_R) {
            err = gretl_model_add_QML_vcv(pmod, LOGIT, H, G,
                                          dset, opt, NULL);
        } else {
            err = gretl_model_add_hessian_vcv(pmod, H);
        }
    }

    gretl_matrix_free(H);
    gretl_matrix_free(G);

    return err;
}

/* Construct 'yhat' and 'uhat'.  Maybe this is too simple-minded; we
   just find, for each observation, the y-value for which the
   probability is maximized and set that as yhat[t].  We then set
   uhat[t] as a binary "hit" (residual = 0) or "miss" (residual = 1).
   Constructing a "quantitative" residual as y[t] - yhat[t] seems
   spurious, since there's no meaningful metric for the "distance"
   between y and yhat when y is an unordered response.

   Note: @yvals is non-NULL if and only if we had to transform the
   dependent variable, because it did not form a 0-based sequence of
   consecutive integers.
*/

static void mn_logit_yhat (MODEL *pmod, mnl_info *mnl,
                           const int *yvals)
{
    double p, pmax;
    int i, s, t, yidx;
    int ncorrect = 0;

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (na(pmod->yhat[t])) {
            continue;
        }
        pmax = 0.0;
        yidx = 0;
        for (i=0; i<mnl->n; i++) {
            p = gretl_matrix_get(mnl->Xb, s, i);
            if (p > pmax) {
                pmax = p;
                yidx = i + 1;
            }
        }
        if (yidx == (int) mnl->y->val[s]) {
            ncorrect++;
            pmod->uhat[t] = 0;
        } else {
            pmod->uhat[t] = 1;
        }
        if (yvals != NULL) {
            pmod->yhat[t] = yvals[yidx];
        } else {
            pmod->yhat[t] = yidx;
        }
        s++;
    }

    gretl_model_set_int(pmod, "correct", ncorrect);
}

/**
 * mn_logit_probabilities:
 * @pmod: pointer to multinomial logit model
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Computes the estimated probabilities of the outcomes
 * for a multinomial logit model. The returned matrix
 * is n x m, where n is the number of observations in
 * the sample range over which the model was estimated
 * and m is the number of distinct outcomes. Each element
 * represents the conditional probability of outcome j
 * given the values of the regressors at observation i.
 *
 * If any of the regressor values are missing at a given
 * observation the probability is set to NaN; provided the
 * regressor information is complete we compute the
 * outcome probabilities even if the actual outcome is
 * missing.
 *
 * Returns: allocated matrix or NULL on failure.
 */

gretl_matrix *mn_logit_probabilities (const MODEL *pmod,
                                      int t1, int t2,
                                      const DATASET *dset,
                                      int *err)
{
    gretl_matrix *P = NULL;
    const gretl_matrix *yvals = NULL;
    const double *b = NULL;
    double *eXbt = NULL;
    double St, ptj;
    char **S = NULL;
    int j, k, s, t, vi, ok;
    int i, m = 0;

    if (pmod == NULL || pmod->list == NULL || pmod->coeff == NULL) {
        *err = E_DATA;
        return NULL;
    }

    /* list of outcome values (including the base case) */
    yvals = gretl_model_get_data(pmod, "yvals");
    if (yvals == NULL) {
        *err = E_DATA;
    } else {
        m = gretl_vector_get_length(yvals);
    }

    if (!*err) {
        for (i=1; i<=pmod->list[0]; i++) {
            if (pmod->list[i] >= dset->v) {
                /* a regressor has disappeared */
                *err = E_DATA;
                break;
            }
        }
    }

    if (!*err) {
        int n = t2 - t1 + 1;

        P = gretl_matrix_alloc(n, m);
        if (P == NULL) {
            *err = E_ALLOC;
        } else {
            S = strings_array_new(n);
        }
    }

    if (!*err) {
        /* allocate required workspace */
        eXbt = malloc(m * sizeof *eXbt);
        if (eXbt == NULL) {
            *err = E_ALLOC;
        }
    }

    if (*err) {
        goto bailout;
    }

    b = pmod->coeff;

    for (t=t1, s=0; t<=t2; t++, s++) {
        ok = 1;
        for (i=2; i<=pmod->list[0]; i++) {
            vi = pmod->list[i];
            if (na(dset->Z[vi][t])) {
                ok = 0;
                break;
            }
        }
        if (!ok) {
            /* one or more regressors missing */
            for (j=0; j<m; j++) {
                gretl_matrix_set(P, s, j, NADBL);
            }
        } else {
            /* base case */
            eXbt[0] = St = 1.0;
            k = 0;
            /* loop across the other y-values */
            for (j=1; j<m; j++) {
                /* accumulate exp(X*beta) */
                eXbt[j] = 0.0;
                for (i=2; i<=pmod->list[0]; i++) {
                    vi = pmod->list[i];
                    eXbt[j] += dset->Z[vi][t] * b[k++];
                }
                eXbt[j] = exp(eXbt[j]);
                St += eXbt[j];
            }
            for (j=0; j<m; j++) {
                ptj = eXbt[j] / St;
                gretl_matrix_set(P, s, j, ptj);
            }
        }
        if (S != NULL) {
            S[s] = retrieve_date_string(t+1, dset, err);
        }
    }

    if (P != NULL) {
        gretl_matrix_set_t1(P, t1);
        gretl_matrix_set_t2(P, t2);
        if (S != NULL) {
            gretl_matrix_set_rownames(P, S);
        }
    }

 bailout:

    free(eXbt);

    return P;
}

/**
 * binary_logit_odds_ratios:
 * @pmod: pointer to binary logit model
 * @dset: dataset struct.
 *
 * Computes odds ratios, their standard errors, and 95 percent
 * confidence intervals for a binary logit model. Attaches to
 * @pmod a labeled matrix representation of the results.
 *
 * Returns: 0 on success, error code on error.
 */

static int binary_logit_odds_ratios (MODEL *pmod,
                                     const DATASET *dset)
{
    CoeffIntervals *cf;
    gretl_matrix *m;
    int err = 0;

    cf = gretl_model_get_coeff_intervals(pmod, dset, OPT_O | OPT_E);
    if (cf != NULL) {
        m = conf_intervals_matrix(cf);
        err = gretl_model_set_matrix_as_data(pmod, "oddsratios", m);
        free_coeff_intervals(cf);
    }

    return err;
}

/* In case the dependent variable is not in canonical form for
   multinomial logit, construct a transformed version */

static int make_canonical_depvar (MODEL *pmod, const double *y,
                                  gretl_matrix *yvec)
{
    struct sorter *s;
    double nexty, bady;
    int i, t, n = pmod->nobs;

    s = malloc(n * sizeof *s);
    if (s == NULL) {
        return E_ALLOC;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (!na(pmod->uhat[t])) {
            s[i].x = y[t];
            s[i].t = i;
            i++;
        }
    }

    qsort(s, n, sizeof *s, gretl_compare_doubles);

    /* normalize to a minimum of zero */
    if (s[0].x != 0) {
        double ymin = s[0].x;

        for (i=0; i<n && s[i].x == ymin; i++) {
            s[i].x = 0.0;
        }
    }

    /* ensure that the sorted values increase by steps of one */
    for (i=1; i<n; i++) {
        if (s[i].x != s[i-1].x) {
            nexty = s[i-1].x + 1;
            if (s[i].x != nexty) {
                bady = s[i].x;
                while (i < n && s[i].x == bady) {
                    s[i++].x = nexty;
                }
                i--; /* compensate for outer i++ */
            }
        }
    }

    /* write canonical version of y into yvec */
    for (i=0; i<n; i++) {
        yvec->val[s[i].t] = s[i].x;
    }

    free(s);

    return 0;
}

/* multinomial logit: count the distinct values taken on by the
   dependent variable.  If the variable does not take the form of a
   sequence of consecutive integers with a base of zero, flag this
   by returning via @yvals the vector of distinct values. To assist
   with later computations, we also return a frequency count for
   the distinct y values in @valcount.
*/

static int mn_value_count (const double *y, MODEL *pmod,
                           int **valcount, int **yvals)
{
    double *sy;
    int *vc = NULL;
    int *v = NULL;
    int want_yvals = 0;
    int s, t, n;

    sy = malloc(pmod->nobs * sizeof *sy);
    if (sy == NULL) {
        pmod->errcode = E_ALLOC;
        return 0;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (!na(pmod->uhat[t])) {
            sy[s++] = y[t];
        }
    }

    qsort(sy, pmod->nobs, sizeof *sy, gretl_compare_doubles);

    n = count_distinct_values(sy, pmod->nobs);

    vc = malloc(n * sizeof *vc);
    v = malloc(n * sizeof *v);

    if (vc == NULL || v == NULL) {
        pmod->errcode = E_ALLOC;
        free(sy);
        free(vc);
        return 0;
    }

    for (s=0; s<n; s++) {
        vc[s] = 0;
    }

    if (sy[0] != 0.0) {
        /* we'll need v later */
        want_yvals = 1;
    }

    v[0] = (int) sy[0];
    s = 0;

    for (t=0; t<pmod->nobs && !pmod->errcode; t++) {
        if (sy[t] != floor(sy[t])) {
            pmod->errcode = E_DATA;
            gretl_errmsg_set(_("logit: the dependent variable must form a sequence of "
                             "integers"));
        } else if (t > 0 && sy[t] != sy[t-1]) {
            v[++s] = (int) sy[t];
            if (sy[t] != sy[t-1] + 1) {
                /* not consecutive: need v later */
                want_yvals = 1;
            }
            vc[s] = 1;
        } else {
            vc[s] += 1;
        }
    }

    if (pmod->errcode) {
        free(vc);
        n = 0;
    } else {
        *valcount = vc;
        if (want_yvals) {
            *yvals = v;
            v = NULL;
        }
    }

    free(v);
    free(sy);

    return n;
}

/* transcribe multinomial logit results into @pmod, which was
   initialized via OLS, and add covariance matrix
*/

static void mnl_finish (mnl_info *mnl, MODEL *pmod,
                        const int *valcount,
                        const int *yvals,
                        const DATASET *dset,
                        gretlopt opt)
{
    int i;

    pmod->errcode = gretl_model_write_coeffs(pmod, mnl->theta, mnl->npar);

    if (!pmod->errcode) {
        pmod->errcode = mnl_add_variance_matrix(pmod, mnl, dset, opt);
    }

    if (!pmod->errcode) {
        pmod->ci = LOGIT;
        pmod->lnL = mn_logit_loglik(mnl->theta, mnl);
        mle_criteria(pmod, 0);

        gretl_model_set_int(pmod, "multinom", mnl->n);
        gretl_model_set_int(pmod, "cblock", mnl->k);

        pmod->ess = pmod->sigma = NADBL;
        pmod->fstt = pmod->rsq = pmod->adjrsq = NADBL;

        gretl_model_allocate_param_names(pmod, mnl->npar);
    }

    if (!pmod->errcode) {
        int j, vj, k = 0;

        for (i=0; i<mnl->n; i++) {
            for (j=0; j<mnl->k; j++) {
                vj = pmod->list[j+2];
                gretl_model_set_param_name(pmod, k++, dset->varname[vj]);
            }
        }
        mn_logit_yhat(pmod, mnl, yvals);
        set_model_id(pmod, opt);
    }

    if (!pmod->errcode) {
        /* add a record of the values of the dependent variable */
        gretl_matrix *yv = gretl_column_vector_alloc(mnl->n + 1);

        if (yv != NULL) {
            for (i=0; i<=mnl->n; i++) {
                if (yvals != NULL) {
                    yv->val[i] = yvals[i];
                } else {
                    yv->val[i] = i;
                }
            }
            gretl_model_set_matrix_as_data(pmod, "yvals", yv);
        }
    }

    if (!pmod->errcode) {
        /* add overall likelihood ratio test */
        int ni, df = pmod->ncoeff;
        double L0 = 0.0;

        if (pmod->ifc) {
            df -= mnl->n;
        }
        for (i=0; i<=mnl->n; i++) {
            ni = valcount[i];
            if (pmod->ifc) {
                L0 += ni * log((double) ni / mnl->T);
            } else {
                L0 += ni * log(1.0 / (mnl->n + 1));
            }
        }
        pmod->chisq = 2.0 * (pmod->lnL - L0);
        pmod->dfn = df;
    }

    clear_model_xpx(pmod);

    pmod->opt |= OPT_M;
}

/* multinomial logit */

static MODEL mnl_model (const int *list, DATASET *dset,
                        gretlopt opt, PRN *prn)
{
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;
    MODEL mod;
    mnl_info *mnl;
    int *valcount = NULL;
    int *yvals = NULL;
    int n, k = list[0] - 1;
    gretlopt maxopt;
    int use_bfgs = 0;
    int i, vi, t, s;

    /* we'll start with OLS to flush out data issues */
    mod = lsq(list, dset, OLS, OPT_A);
    if (mod.errcode) {
        return mod;
    }

    n = mn_value_count(dset->Z[list[1]], &mod, &valcount, &yvals);
    if (mod.errcode) {
        return mod;
    }

    if (opt & OPT_C) {
        /* cluster implies robust */
        opt |= OPT_R;
    }

    n--; /* exclude the first value */

    if (n * k > mod.nobs) {
        mod.errcode = E_DF;
        return mod;
    }

    mnl = mnl_info_new(n, k, mod.nobs);
    if (mnl == NULL) {
        mod.errcode = E_ALLOC;
        return mod;
    }

    if (yvals != NULL) {
        /* the dependent variable needs transforming */
        mod.errcode = make_canonical_depvar(&mod, dset->Z[list[1]], mnl->y);
        if (mod.errcode) {
            goto bailout;
        }
    } else {
        s = 0;
        for (t=mod.t1; t<=mod.t2; t++) {
            if (!na(mod.yhat[t])) {
                mnl->y->val[s++] = dset->Z[list[1]][t];
            }
        }
    }

    for (i=0; i<k; i++) {
        vi = list[i+2];
        s = 0;
        for (t=mod.t1; t<=mod.t2; t++) {
            if (!na(mod.yhat[t])) {
                gretl_matrix_set(mnl->X, s++, i, dset->Z[vi][t]);
            }
        }
    }

    if (mod.ifc) {
        double lf0 = log(valcount[0]);
        int j = 0;

        for (i=1; i<=mnl->n; i++) {
            mnl->theta[j] = log(valcount[i]) - lf0;
            j += mnl->k;
        }
    }

    if (libset_get_int(GRETL_OPTIM) == OPTIM_BFGS) {
        use_bfgs = 1;
    }

    maxopt = (prn != NULL)? (OPT_U | OPT_V) : OPT_U;

    if (use_bfgs) {
        mod.errcode = BFGS_max(mnl->theta, mnl->npar, maxit, 0.0,
                               &fncount, &grcount, mn_logit_loglik, C_LOGLIK,
                               mn_logit_score, mnl, NULL, maxopt, prn);
    } else {
        double crittol = 1.0e-8;
        double gradtol = 1.0e-7;

        maxit = 100;
        mod.errcode = newton_raphson_max(mnl->theta, mnl->npar, maxit,
                                         crittol, gradtol, &fncount,
                                         C_LOGLIK, mn_logit_loglik,
                                         mn_logit_score, mnl_hessian, mnl,
                                         maxopt, prn);
    }

    if (!mod.errcode) {
        mnl_finish(mnl, &mod, valcount, yvals, dset, opt);
    }

 bailout:

    mnl_info_destroy(mnl);
    free(valcount);
    free(yvals);

    return mod;
}

/**
 * biprobit_model:
 * @list: binary dependent variable 1, binary dependent variable 2,
 * list of regressors for y1. If @list ends here, it is assumed that the
 * explanatory variables for y2 are the same as y1. Otherwise, the list
 * must include a separator and the list of regressors for y2.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q for quiet operation, OPT_V for verbose
 * operation, OPT_R for robust covariance matrix, OPT_G for covariance
 * matrix based on Outer Product of Gradient.
 * @prn: printing struct.
 *
 * Computes estimates of the bivariate probit model specified by @list,
 * using maximum likelihood via Newton-Raphson.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL biprobit_model (const int *list, DATASET *dset,
                      gretlopt opt, PRN *prn)
{
    MODEL bpmod;
    MODEL (* biprobit_estimate) (const int *, DATASET *,
                                 gretlopt, PRN *);

    gretl_error_clear();

    biprobit_estimate = get_plugin_function("biprobit_estimate");

    if (biprobit_estimate == NULL) {
        gretl_model_init(&bpmod, dset);
        bpmod.errcode = E_FOPEN;
        return bpmod;
    }

    bpmod = (*biprobit_estimate) (list, dset, opt, prn);
    set_model_id(&bpmod, opt);

    return bpmod;
}

/* struct for holding binary probit/logit info */

typedef struct bin_info_ bin_info;

struct bin_info_ {
    int ci;           /* PROBIT or LOGIT */
    int k;            /* number of parameters */
    int T;            /* number of observations */
    int pp_err;       /* to record perfect-prediction error */
    double *theta;    /* coeffs for Newton-Raphson */
    int *y;           /* dependent variable */
    gretl_matrix *X;  /* regressors */
    gretl_matrix *Ri; /* inverse of R from decomp of X */
    gretl_matrix_block *B;
    gretl_matrix *pX; /* for use with Hessian */
    gretl_matrix *b;  /* coefficients in matrix form */
    gretl_matrix *Xb; /* index function values */
};

static void bin_info_destroy (bin_info *bin)
{
    if (bin != NULL) {
        gretl_matrix_block_destroy(bin->B);
        gretl_matrix_free(bin->X);
        gretl_matrix_free(bin->Ri);
        free(bin->theta);
        free(bin->y);
        free(bin);
    }
}

static bin_info *bin_info_new (int ci, int k, int T)
{
    bin_info *bin = malloc(sizeof *bin);
    int err = 0;

    if (bin == NULL) {
        return NULL;
    }

    bin->X = NULL;
    bin->Ri = NULL;
    bin->B = NULL;
    bin->ci = ci;
    bin->k = k;
    bin->T = T;
    bin->pp_err = 0;
    bin->theta = malloc(k * sizeof *bin->theta);
    bin->y = malloc(T * sizeof *bin->y);
    if (bin->theta == NULL || bin->y == NULL) {
        err = E_ALLOC;
    } else {
        bin->B = gretl_matrix_block_new(&bin->pX, T, k,
                                        &bin->b, k, 1,
                                        &bin->Xb, T, 1,
                                        NULL);
        if (bin->B == NULL) {
            err = E_ALLOC;
        }
    }
    if (err) {
        bin_info_destroy(bin);
        bin = NULL;
    }

    return bin;
}

/* compute loglikelihood for binary probit/logit */

static double binary_loglik (const double *theta, void *ptr)
{
    bin_info *bin = (bin_info *) ptr;
    double max0 = -1.0e200;
    double min1 = 1.0e200;
    double e, ndx, p;
    double ll = 0.0;
    int i, t;

    for (i=0; i<bin->k; i++) {
        bin->b->val[i] = theta[i];
    }
    gretl_matrix_multiply(bin->X, bin->b, bin->Xb);

    /* perfect prediction check */
    for (t=0; t<bin->T; t++) {
        ndx = gretl_vector_get(bin->Xb, t);
        if (bin->y[t] == 0 && ndx > max0) {
            max0 = ndx;
        } else if (bin->y[t] == 1 && ndx < min1) {
            min1 = ndx;
        }
    }
    if (min1 > max0) {
        bin->pp_err = 1;
        return NADBL;
    }

    /* compute loglikelihood */
    for (t=0; t<bin->T; t++) {
        ndx = gretl_vector_get(bin->Xb, t);
        if (bin->ci == PROBIT) {
            p = bin->y[t] ? normal_cdf(ndx) : normal_cdf(-ndx);
        } else {
            e = logit(ndx);
            p = bin->y[t] ? e : 1-e;
        }
        ll += log(p);
    }

    return ll;
}

static int binary_score (double *theta, double *s, int k,
                         BFGS_CRIT_FUNC ll, void *ptr)
{
    bin_info *bin = (bin_info *) ptr;
    double ndx, w;
    int t, j, yt;
    int err = 0;

    for (j=0; j<bin->k; j++) {
        s[j] = 0.0;
    }

    for (t=0; t<bin->T; t++) {
        yt = bin->y[t];
        ndx = gretl_vector_get(bin->Xb, t);
        if (bin->ci == PROBIT) {
            w = yt ? invmills(-ndx) : -invmills(ndx);
        } else {
            w = yt - logit(ndx);
        }
        for (j=0; j<bin->k; j++) {
            s[j] += w * gretl_matrix_get(bin->X, t, j);
        }
    }

    errno = 0;

    return err;
}

/* binary probit/logit: form the negative of the analytical
   Hessian */

static int binary_hessian (double *theta, gretl_matrix *H,
                           void *data)
{
    bin_info *bin = data;
    double w, p, ndx, xtj;
    int t, j;

    for (t=0; t<bin->T; t++) {
        ndx = gretl_vector_get(bin->Xb, t);
        if (bin->ci == PROBIT) {
            w = bin->y[t] ? invmills(-ndx) : -invmills(ndx);
            p = w * (ndx + w);
        } else {
            p = logit(ndx);
            p = p * (1-p);
        }
        for (j=0; j<bin->k; j++) {
            xtj = gretl_matrix_get(bin->X, t, j);
            gretl_matrix_set(bin->pX, t, j, p * xtj);
        }
    }

    gretl_matrix_multiply_mod(bin->pX, GRETL_MOD_TRANSPOSE,
                              bin->X, GRETL_MOD_NONE,
                              H, GRETL_MOD_NONE);

    return 0;
}

static gretl_matrix *binary_hessian_inverse (bin_info *bin, int *err)
{
    gretl_matrix *H;

    H = gretl_zero_matrix_new(bin->k, bin->k);
    if (H == NULL) {
        *err = E_ALLOC;
    } else {
        *err = binary_hessian(bin->theta, H, bin);
    }

    if (!*err) {
        *err = gretl_invert_symmetric_matrix(H);
        if (*err) {
            /* fallback: is this a good idea? (2025-11-04) */
            fprintf(stderr, "Warning: binary model Hessian is not p.d.\n");
            *err = gretl_invert_general_matrix(H);
        }
        if (*err) {
            gretl_errmsg_set("Failed to invert binary model Hessian");
        }
    }

    return H;
}

static gretl_matrix *binary_score_matrix (bin_info *bin, int *err)
{
    gretl_matrix *G;
    double w, ndx, xtj;
    int yt, j, t;

    G = gretl_matrix_alloc(bin->T, bin->k);

    if (G == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (t=0; t<bin->T && !errno; t++) {
        yt = bin->y[t];
        ndx = gretl_vector_get(bin->Xb, t);
        if (bin->ci == PROBIT) {
            w = yt ? invmills(-ndx) : -invmills(ndx);
        } else {
            w = yt - logit(ndx);
        }
        for (j=0; j<bin->k; j++) {
            xtj = gretl_matrix_get(bin->X, t, j);
            gretl_matrix_set(G, t, j, w * xtj);
        }
    }

    return G;
}

static int binary_variance_matrix (MODEL *pmod, bin_info *bin,
                                   const DATASET *dset,
                                   gretlopt opt)
{
    gretl_matrix *H = NULL;
    gretl_matrix *G = NULL;
    int err = 0;

    H = binary_hessian_inverse(bin, &err);

    if (!err && (opt & OPT_R)) {
        G = binary_score_matrix(bin, &err);
    }

    if (!err) {
        if (opt & OPT_R) {
            err = gretl_model_add_QML_vcv(pmod, bin->ci, H, G,
                                          dset, opt, NULL);
        } else {
            err = gretl_model_add_hessian_vcv(pmod, H);
        }
    }

    gretl_matrix_free(H);
    gretl_matrix_free(G);

    return err;
}

static void binary_model_chisq (bin_info *bin, MODEL *pmod,
                                gretlopt opt)
{
    double L0;

    if (pmod->ncoeff == 1 && pmod->ifc) {
        /* constant-only model */
        pmod->chisq = NADBL;
        pmod->rsq = 0;
        pmod->adjrsq = NADBL;
        return;
    }

    L0 = binary_null_loglik(bin->y, bin->T);

    if (!na(L0) && L0 <= pmod->lnL) {
        pmod->chisq = 2.0 * (pmod->lnL - L0);
        add_pseudo_rsquared(pmod, L0, bin->k, bin->T, opt);
    } else {
        pmod->rsq = pmod->adjrsq = pmod->chisq = NADBL;
    }
}

/* Binary probit normality test as in Bera, Jarque and Lee
   (International Economic Review, 1984), also quoted in Verbeek,
   chapter 7: we regress a column of 1s on the products of the
   generalized residual with X, (X\beta)^2 and (X\beta)^3.  The test
   statistic is T times the uncentered R-squared, and is distributed as
   chi-square(2). It can be shown that this test is numerically
   identical to the Chesher-Irish (87) test in the probit case
   (although C&I make no mention of this in their article).
*/

static int binary_probit_normtest (MODEL *pmod, bin_info *bin)
{
    gretl_matrix_block *B;
    gretl_matrix *X, *y, *b;
    double xti, Xb, et;
    int k = bin->k;
    int i, s, t;
    int err = 0;

    B = gretl_matrix_block_new(&X, bin->T, k+2,
                               &y, bin->T, 1,
                               &b, k+2, 1,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
        if (model_missing(pmod, t)) {
            continue;
        }
        et = pmod->uhat[t];
        Xb = gretl_vector_get(bin->Xb, s);
        for (i=0; i<k; i++) {
            xti = gretl_matrix_get(bin->X, s, i);
            gretl_matrix_set(X, s, i, et * xti);
        }
        gretl_matrix_set(X, s, k, et * Xb * Xb);
        gretl_matrix_set(X, s, k+1, et * Xb * Xb * Xb);
        gretl_vector_set(y, s, 1);
        s++;
    }

    err = gretl_matrix_ols(y, X, b, NULL, NULL, NULL);

    if (!err) {
        double X2 = bin->T;

        gretl_matrix_multiply(X, b, y);
        for (t=0; t<y->rows; t++) {
            X2 -= (1 - y->val[t]) * (1 - y->val[t]);
        }
        if (X2 > 0) {
            gretl_model_add_normality_test(pmod, X2);
        }
    }

    gretl_matrix_block_destroy(B);

    return err;
}

/* Special "slope" calculation for a dummy regressor in binary model:
   this is F(~Xb + b_j) - F(~Xb) where ~Xb denotes the sum of
   (coefficient times mean) for all regressors other than the dummy in
   question and b_j indicates the coefficient on the dummy.  That is,
   the calculation measures the effect on the probability of Y = 1 of
   the discrete change 0 to 1 in x_j.
*/

static double dumslope (MODEL *pmod, const double *xbar, int j)
{
    double s, Xb = 0.0;
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
        if (i != j) {
            Xb += pmod->coeff[i] * xbar[i];
        }
    }

    if (pmod->ci == LOGIT) {
        s = logit(Xb + pmod->coeff[j]) - logit(Xb);
    } else {
        s = normal_cdf(Xb + pmod->coeff[j]) - normal_cdf(Xb);
    }

    return s;
}

static int binary_model_add_slopes (MODEL *pmod, bin_info *bin,
                                    const DATASET *dset)
{
    double *xbar;
    double *slopes;
    double Xb, fXb;
    size_t ssize;
    int i, vi, t;
    int err = 0;

    xbar = malloc(bin->k * sizeof *xbar);
    ssize = bin->k * sizeof *slopes;
    slopes = malloc(ssize);

    if (slopes == NULL || xbar == NULL) {
        free(xbar);
        free(slopes);
        return E_ALLOC;
    }

    Xb = 0.0;
    for (i=0; i<bin->k; i++) {
        xbar[i] = 0.0;
        vi = pmod->list[i+2];
        for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t])) {
                xbar[i] += dset->Z[vi][t];
            }
        }
        xbar[i] /= bin->T;
        Xb += pmod->coeff[i] * xbar[i];
    }

    fXb = (bin->ci == LOGIT)? logit_pdf(Xb) : normal_pdf(Xb);
    gretl_model_set_double(pmod, "fXb", fXb);

    for (i=0; i<bin->k; i++) {
        vi = pmod->list[i+2];
        if (vi == 0) {
            slopes[i] = 0.0;
        } else if (gretl_isdummy(pmod->t1, pmod->t2, dset->Z[vi])) {
            slopes[i] = dumslope(pmod, xbar, i);
        } else {
            slopes[i] = pmod->coeff[i] * fXb;
        }
    }

    err = gretl_model_set_data(pmod, "slopes", slopes,
                               GRETL_TYPE_DOUBLE_ARRAY,
                               ssize);
    free(xbar);
    if (err) {
        free(slopes);
    }

    return err;
}

static double binary_model_fXb (MODEL *pmod, bin_info *bin,
                                const DATASET *dset)
{
    double xbar, Xb = 0.0;
    int i, vi, t;

    for (i=0; i<bin->k; i++) {
        xbar = 0.0;
        vi = pmod->list[i+2];
        for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t])) {
                xbar += dset->Z[vi][t];
            }
        }
        xbar /= bin->T;
        Xb += pmod->coeff[i] * xbar;
    }

    return bin->ci == LOGIT ? logit_pdf(Xb) : normal_pdf(Xb);
}

void binary_model_hatvars (MODEL *pmod,
                           const gretl_matrix *ndx,
                           const int *y,
                           gretlopt opt)
{
    int *act_pred;
    double *ll = NULL;
    double F;
    int n = pmod->full_n;
    int i, s, t;

    /* add space for actual/predicted matrix */
    act_pred = malloc(4 * sizeof *act_pred);
    if (act_pred != NULL) {
        for (i=0; i<4; i++) {
            act_pred[i] = 0;
        }
    }

    if (!(opt & OPT_E)) {
        /* OPT_E indicates random effects */
        ll = malloc(n * sizeof *ll);
        if (ll != NULL) {
            for (t=0; t<n; t++) {
                ll[t] = NADBL;
            }
        }
    }

    errno = 0;

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
        double ndxt, endxt;
        int yt;

        if (model_missing(pmod, t)) {
            continue;
        }

        ndxt = gretl_vector_get(ndx, s);
        yt = y[s++];

        if (act_pred != NULL) {
            i = 2 * yt + (ndxt > 0.0);
            act_pred[i] += 1;
        }

        if (pmod->ci == LOGIT) {
            endxt = exp(ndxt);
            if (errno == ERANGE) {
                F = (ndxt > 0)? 1 : 0;
                errno = 0;
            } else {
                F = endxt / (1.0 + endxt);
            }
            pmod->yhat[t] = F;
            pmod->uhat[t] = yt - pmod->yhat[t];
        } else {
            F = normal_cdf(ndxt);
            pmod->yhat[t] = F;
            pmod->uhat[t] = yt ? invmills(-ndxt) : -invmills(ndxt);
        }

        if (ll != NULL) {
            ll[t] = yt ? log(F) : log(1-F);
        }
    }

    if (act_pred != NULL) {
        gretl_model_set_data(pmod, "discrete_act_pred", act_pred,
                             GRETL_TYPE_INT_ARRAY,
                             4 * sizeof *act_pred);
    }
    if (ll != NULL) {
        gretl_model_set_data(pmod, "llt", ll,
                             GRETL_TYPE_DOUBLE_ARRAY,
                             n * sizeof *ll);
    }
}

/* After a binary model is estimated via QR, multiply the
   coefficients by R^{-1} to convert back to the statistics
   of the original data.
*/

static void binary_qr_finish_coeffs (MODEL *pmod, bin_info *bin)
{
    gretl_matrix bt = {0};
    gretl_matrix pc = {0};

    gretl_matrix_init_full(&bt, bin->k, 1, bin->theta);
    gretl_matrix_init_full(&pc, bin->k, 1, pmod->coeff);
    gretl_matrix_multiply(bin->Ri, &bt, &pc);
}

/* After a binary model is estimated via QR, correct its covariance
   matrix, V, which is relative to the decomposed data, to
   R^{-1} * V * R^{-1}'.
*/

static int binary_qr_finish_vcv (MODEL *pmod, bin_info *bin)
{
    gretl_matrix *V = NULL;
    gretl_matrix *RVR = NULL;
    int err = 0;

    RVR = gretl_matrix_alloc(bin->k, bin->k);
    if (RVR == NULL) {
        err = E_ALLOC;
    } else {
        V = gretl_vcv_matrix_from_model(pmod, NULL, &err);
    }

    if (!err) {
        gretl_matrix_qform(bin->Ri, GRETL_MOD_NONE,
                           V, RVR, GRETL_MOD_NONE);
        gretl_model_write_vcv(pmod, RVR);
    }

    gretl_matrix_free(RVR);
    gretl_matrix_free(V);

    return err;
}

static int binary_model_finish (bin_info *bin,
                                MODEL *pmod,
                                const DATASET *dset,
                                gretlopt opt)
{
    memcpy(pmod->coeff, bin->theta, bin->k * sizeof(double));
    pmod->ci = bin->ci;
    pmod->lnL = binary_loglik(pmod->coeff, bin);

    if (pmod->ci == PROBIT && (opt & OPT_X)) {
        /* this model is just the starting-point for
           random-effects probit estimation
        */
        pmod->opt |= OPT_P;
        return 0;
    }

    binary_model_chisq(bin, pmod, opt);
    pmod->errcode = binary_variance_matrix(pmod, bin, dset, opt);

    if (!pmod->errcode) {
        binary_qr_finish_coeffs(pmod, bin);
        if (opt & OPT_P) {
            /* showing p-values, not slopes */
            double fXb = binary_model_fXb(pmod, bin, dset);

            pmod->opt |= OPT_P;
            gretl_model_set_double(pmod, "fXb", fXb);
        } else {
            pmod->errcode = binary_model_add_slopes(pmod, bin, dset);
        }
    }

    if (!pmod->errcode) {
        binary_model_hatvars(pmod, bin->Xb, bin->y, opt);
        mle_criteria(pmod, 0);
        if (opt & OPT_A) {
            pmod->aux = AUX_AUX;
        } else {
            set_model_id(pmod, opt);
        }
        gretl_model_set_int(pmod, "binary", 1);
    }

    if (!pmod->errcode && bin->Ri != NULL) {
        pmod->errcode = binary_qr_finish_vcv(pmod, bin);
    }

    if (!pmod->errcode) {
        if (pmod->ci == PROBIT) {
            binary_probit_normtest(pmod, bin);
        } else {
            binary_logit_odds_ratios(pmod, dset);
        }
    }

    return pmod->errcode;
}

/* revise_lpm_coeffs: here we're revising the coefficients we got via
   initial OLS estimation (Linear Probability Model) when it turns out
   OLS used QR decomposition and therefore we're going to use QR in ML
   estimation of the binary model. In the non-QR case we just divide the
   coeffs by pmod->sigma, but in the QR case we first obtain revised
   coeffs as Q'y (where in context Q = bin->X).
*/

static int revise_lpm_coeffs (bin_info *bin,
                              MODEL *pmod)
{
    gretl_matrix *y = NULL;
    gretl_matrix b = {};
    int i, err = 0;

    y = gretl_matrix_alloc(bin->T, 1);
    for (i=0; i<pmod->nobs; i++) {
        y->val[i] = bin->y[i];
    }
    gretl_matrix_init_full(&b, bin->k, 1, bin->theta);
    gretl_matrix_multiply_mod(bin->X, GRETL_MOD_TRANSPOSE,
                              y, GRETL_MOD_NONE,
                              &b, GRETL_MOD_NONE);
    gretl_matrix_divide_by_scalar(&b, pmod->sigma);
    gretl_matrix_free(y);

    return err;
}

static int make_binary_y_and_X (bin_info *bin,
                                MODEL *pmod,
                                DATASET *dset)
{
    int v = pmod->list[1];
    int i, t, s;
    int err = 0;

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
        if (!na(pmod->yhat[t])) {
            bin->y[s++] = dset->Z[v][t] ? 1.0 : 0.0;
        }
    }

    if (gretl_model_get_int(pmod, "QR")) {
        fprintf(stderr, "make_binary_y_and_X: QR case\n");
        bin->X = gretl_model_steal_data(pmod, "Q");
        bin->Ri = gretl_model_steal_data(pmod, "R");
        err = revise_lpm_coeffs(bin, pmod);
    } else {
        bin->X = gretl_matrix_alloc(bin->T, bin->k);
        if (bin->X == NULL) {
            err = E_ALLOC;
        }
        for (i=0; i<bin->k && !err; i++) {
            bin->theta[i] = pmod->coeff[i] / pmod->sigma;
        }
        for (i=0; i<bin->k && !err; i++) {
            v = pmod->list[i+2];
            for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
                if (!na(pmod->yhat[t])) {
                    gretl_matrix_set(bin->X, s++, i, dset->Z[v][t]);
                }
            }
        }
    }

    return err;
}

MODEL binary_model (int ci, const int *list,
                    DATASET *dset, gretlopt opt,
                    PRN *prn)
{
    gretlopt maxopt = OPT_NONE;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *blist = NULL;
    char *mask = NULL;
    double crittol = 1.0e-8;
    double gradtol = 1.0e-7;
    int maxit = 100;
    int fncount = 0;
    int ndropped = 0;
    MODEL mod;
    bin_info *bin = NULL;
    PRN *vprn = NULL;
    int depvar;

    gretl_model_init(&mod, dset);
    depvar = list[1];

    if (!gretl_isdummy(dset->t1, dset->t2, dset->Z[depvar])) {
        gretl_errmsg_sprintf(_("'%s' is not a binary variable"),
                             dset->varname[depvar]);
        mod.errcode = E_DATA;
        return mod;
    }

    blist = gretl_list_copy(list);
    if (blist == NULL) {
        mod.errcode = E_ALLOC;
        goto bailout;
    }

    /* the cluster option implies robust */
    if (opt & OPT_C) {
        opt |= OPT_R;
    }

    list_adjust_sample(blist, &dset->t1, &dset->t2, dset, NULL);

    mask = classifier_check(blist, dset, prn, &ndropped, &mod.errcode);
    if (mod.errcode) {
        if (mod.errcode == E_NOCONV) {
            gretl_errmsg_set(_("Perfect prediction obtained: no MLE exists"));
        }
    } else if (mask != NULL) {
        mod.errcode = copy_to_reference_missmask(mask);
    }

    if (!mod.errcode) {
        gretlopt ols_opt = OPT_A | OPT_B;

        /* If we're doing logit/probit as an auxiliary regression,
           it might be safer to abort on perfect collinearity;
           otherwise we'll try automatic elimination.
        */
        if (opt & OPT_A) {
            ols_opt |= OPT_Z;
        }
        mod = lsq(blist, dset, OLS, ols_opt);
#if LPDEBUG
        printmodel(&mod, dset, OPT_NONE, prn);
#endif
    }

    if (!mod.errcode) {
        bin = bin_info_new(ci, mod.ncoeff, mod.nobs);
        if (bin == NULL) {
            mod.errcode = E_ALLOC;
        } else {
            mod.errcode = make_binary_y_and_X(bin, &mod, dset);
        }
    }
    if (mod.errcode) {
        goto bailout;
    }

    if (opt & OPT_V) {
        /* respect verbosity */
        maxopt = OPT_V;
        vprn = prn;
    }
    if (!(opt & OPT_X)) {
        /* not just an auxiliary estimator: respect user choices
           for Newton parameters
        */
        maxopt |= OPT_U;
    }

    mod.errcode = newton_raphson_max(bin->theta, bin->k, maxit,
                                     crittol, gradtol, &fncount,
                                     C_LOGLIK, binary_loglik,
                                     binary_score, binary_hessian, bin,
                                     maxopt, vprn);
    if (bin->pp_err) {
        /* trash any existing error message */
        gretl_error_clear();
        mod.errcode = E_NOCONV;
        gretl_errmsg_set(_("Perfect prediction obtained: no MLE exists"));
    }
    if (!mod.errcode) {
        binary_model_finish(bin, &mod, dset, opt);
        if (!mod.errcode && ndropped > 0) {
            gretl_model_set_int(&mod, "binary_obs_dropped", ndropped);
        }
    }

 bailout:

    bin_info_destroy(bin);
    free(blist);
    free(mask);

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return mod;
}

/**
 * binary_logit:
 * @list: binary dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: if includes OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix; if OPT_P arrange for
 * printing of p-values, not slopes at mean; if OPT_A treat as an
 * auxiliary regression.
 * @prn: printing struct in case additional information is
 * wanted (in which case add OPT_V to @opt).
 *
 * Computes estimates of the logit model specified by @list,
 * using Newton-Raphson.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL binary_logit (const int *list, DATASET *dset,
                    gretlopt opt, PRN *prn)
{
    return binary_model(LOGIT, list, dset, opt, prn);
}

/**
 * binary_probit:
 * @list: binary dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: if includes OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix; if OPT_P arrange for
 * printing of p-values, not slopes at mean; if OPT_A treat as an
 * auxiliary regression.
 * @prn: printing struct in case additional information is
 * wanted (in which case add OPT_V to @opt).
 *
 * Computes estimates of the probit model specified by @list,
 * using Newton-Raphson.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL binary_probit (const int *list, DATASET *dset,
                     gretlopt opt, PRN *prn)
{
    return binary_model(PROBIT, list, dset, opt, prn);
}

/**
 * ordered_logit:
 * @list: ordinal dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: if includes OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix; if includes OPT_V, produce verbose
 * output.
 * @prn: printing struct in case additional information is
 * wanted (OPT_V).
 *
 * Computes ML estimates of the ordered logit model specified by @list.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL ordered_logit (const int *list, DATASET *dset,
                     gretlopt opt, PRN *prn)
{
    PRN *vprn = (opt & OPT_V)? prn : NULL;

    return ordered_estimate(LOGIT, list, dset, opt, vprn);
}

/**
 * ordered_probit:
 * @list: ordinal dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: if includes OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix; if includes OPT_V, produce verbose
 * output.
 * @prn: printing struct in case additional information is
 * wanted (OPT_V).
 *
 * Computes ML estimates of the ordered probit model specified by @list.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL ordered_probit (const int *list, DATASET *dset,
                      gretlopt opt, PRN *prn)
{
    PRN *vprn = (opt & OPT_V)? prn : NULL;

    return ordered_estimate(PROBIT, list, dset, opt, vprn);
}

/**
 * multinomial_logit:
 * @list: discrete dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: if includes OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix; if includes OPT_V, produce verbose
 * output.
 * @prn: printing struct in case additional information is
 * wanted (OPT_V).
 *
 * Computes ML estimates of the multinomial model specified by @list.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL multinomial_logit (const int *list, DATASET *dset,
                         gretlopt opt, PRN *prn)
{
    PRN *vprn = (opt & OPT_V)? prn : NULL;

    return mnl_model(list, dset, opt, vprn);
}

/**
 * logistic_ymax_lmax:
 * @y: data series.
 * @dset: dataset information.
 * @ymax: location to receive max(y).
 * @lmax: location to receive a guess at a suitable
 * asymptote for a logistic curve fitted to @y.
 *
 * Checks that the non-missing values of @y are all positive,
 * and if so writes the maximum value of @y to @ymax. The
 * value written to @lmax is 1 if max(y) < 1, else 100
 * if max(y) < 100, else 1.1 * max(y).
 *
 * Returns: 0 on success, non-zero on error.
 */

int logistic_ymax_lmax (const double *y, const DATASET *dset,
                        double *ymax, double *lmax)
{
    int t;

    *ymax = 0.0;

    for (t=dset->t1; t<=dset->t2; t++) {
        if (na(y[t])) {
            continue;
        }
        if (y[t] <= 0.0) {
            gretl_errmsg_set(_("Illegal non-positive value of the "
                               "dependent variable"));
            return 1;
        }
        if (y[t] > *ymax) {
            *ymax = y[t];
        }
    }

    if (*ymax < 1.0) {
        *lmax = 1.0;
    } else if (*ymax < 100.0) {
        *lmax = 100.0;
    } else {
        /* admittedly arbitrary */
        *lmax = 1.1 * *ymax;
    }

    return 0;
}

static int make_logistic_depvar (DATASET *dset, int dv,
                                 double lmax)
{
    int t, v = dset->v;
    int err;

    err = dataset_add_series(dset, 1);

    if (!err) {
        for (t=0; t<dset->n; t++) {
            double p = dset->Z[dv][t];

            if (na(p)) {
                dset->Z[v][t] = NADBL;
            } else {
                dset->Z[v][t] = log(p / (lmax - p));
            }
        }
    }

    return err;
}

static int rewrite_logistic_stats (const DATASET *dset,
                                   MODEL *pmod, int dv,
                                   double lmax)
{
    double den, lam = M_PI/5.35;
    double x, sigma, ess, jll, jaic;
    int use_approx = 1;
    int t;

    if (pmod->depvar == NULL || *pmod->depvar == '\0') {
        if (pmod->depvar == NULL) {
            free(pmod->depvar);
        }
        pmod->depvar = gretl_strdup(dset->varname[dv]);
    }

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, dset->Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, dset->Z[dv]);

    if (pmod->vcv == NULL) {
        /* make the VCV matrix before messing with the model stats */
        makevcv(pmod, pmod->sigma);
    }

    if (use_approx) {
        double s2 = pmod->sigma * pmod->sigma;

        den = sqrt(1.0/(lam * lam) + s2);
    }

    ess = 0.0;
    jll = pmod->lnL;

    for (t=pmod->t1; t<=pmod->t2; t++) {
        x = pmod->yhat[t];
        if (!na(x)) {
            if (use_approx) {
                /* approximation via normal CDF */
                pmod->yhat[t] = lmax * normal_cdf(x/den);
            } else {
                /* naive version */
                pmod->yhat[t] = lmax / (1.0 + exp(-x));
            }
            pmod->uhat[t] = dset->Z[dv][t] - pmod->yhat[t];
            ess += pmod->uhat[t] * pmod->uhat[t];
            x = dset->Z[dv][t];
            jll -= log(x * (lmax-x) / lmax);
        }
    }

    sigma = sqrt(ess / pmod->dfd);
    jaic = -2.0 * jll + 2 * pmod->ncoeff;

    pmod->list[1] = dv;
    gretl_model_set_double(pmod, "lmax", lmax);
    gretl_model_set_double(pmod, "ess_orig", ess);
    gretl_model_set_double(pmod, "sigma_orig", sigma);
    gretl_model_set_double(pmod, "jll", jll);
    gretl_model_set_double(pmod, "jaic", jaic);
    pmod->ci = LOGISTIC;

    if (!(pmod->opt & OPT_F)) {
        ls_criteria(pmod);
    }

    return 0;
}

static double get_real_lmax (const double *y,
                             const DATASET *dset,
                             double user_lmax)
{
    double ymax, lmax;
    int err;

    err = logistic_ymax_lmax(y, dset, &ymax, &lmax);
    if (err) {
        return NADBL;
    }

    if (!na(user_lmax)) {
        if (user_lmax <= ymax) {
            gretl_errmsg_set(_("Invalid value for the maximum of the "
                               "dependent variable"));
            lmax = NADBL;
        } else {
            /* respect the user's choice */
            lmax = user_lmax;
        }
    }

    return lmax;
}

/**
 * logistic_model:
 * @list: dependent variable plus list of regressors.
 * @lmax: value for the asymptote of the logistic curve, or
 * %NADBL for automatic treatment of this.
 * @dset: pointer to dataset.
 * @opt: option flags.
 *
 * Estimate the model given in @list using the logistic transformation
 * of the dependent variable.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL logistic_model (const int *list, double lmax,
                      DATASET *dset, gretlopt opt)
{
    int *llist = NULL;
    int dv = list[1];
    double real_lmax;
    MODEL lmod;
    int err = 0;

    gretl_model_init(&lmod, dset);

    if ((opt & OPT_F) && !dataset_is_panel(dset)) {
        gretl_errmsg_set(_("This estimator requires panel data"));
        err = E_DATA;
    } else {
        llist = gretl_list_copy(list);
        if (llist == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        real_lmax = get_real_lmax(dset->Z[dv], dset, lmax);
        if (na(real_lmax)) {
            err = E_DATA;
        }
    }

    if (!err) {
        err = make_logistic_depvar(dset, dv, real_lmax);
    }

    if (err) {
        free(llist);
        lmod.errcode = err;
        return lmod;
    }

    /* replace with transformed dependent variable */
    llist[1] = dset->v - 1;

    if (opt & OPT_C) {
        set_cluster_vcv_ci(LOGISTIC);
    }

    if (opt & OPT_F) {
        lmod = panel_model(llist, dset, opt, NULL);
    } else {
        lmod = lsq(llist, dset, OLS, opt | OPT_A);
    }
    if (!lmod.errcode) {
        rewrite_logistic_stats(dset, &lmod, dv, real_lmax);
        if (!(opt & OPT_F)) {
            /* already done, for the fixed-effects case */
            set_model_id(&lmod, OPT_NONE);
        }
    }

    if (opt & OPT_C) {
        set_cluster_vcv_ci(0);
    }

    dataset_drop_last_variables(dset, 1);
    free(llist);

    return lmod;
}

static double choose (double n, double k, int *err)
{
    double c = 1.0;
    int i;

    for (i=0; i<k; i++) {
        c *= (n - i) / (k - i);
    }

    if (na(c)) {
        *err = 1;
    }

    return c;
}

static double table_prob (double a, double b, double c, double d,
                          double n, int *err)
{
    double p1, p2, p3, P = NADBL;

    p1 = choose(a+b, a, err);
    if (!*err) {
        p2 = choose(c+d, c, err);
    }
    if (!*err) {
        p3 = choose(n, a+c, err);
    }

    if (!*err) {
        P = p1 * p2 / p3;
        if (na(P)) {
            *err = 1;
            P = NADBL;
        }
    }

#if 0
    fprintf(stderr, "\ntable_prob: a=%g, b=%g, c=%g, d=%g, n=%g\n",
            a, b, c, d, n);
    fprintf(stderr, " p1=%g, p2=%g, p3=%g; P = %g\n",
            p1, p2, p3, P);
#endif

    return P;
}

/**
 * fishers_exact_test:
 * @tab: pointer to 2 x 2 cross-tabulation struct.
 * @prn: gretl printer.
 *
 * Computes and prints to @prn the p-value for Fisher's Exact Test for
 * association between the two variables represented in @tab.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int fishers_exact_test (const Xtab *tab, PRN *prn)
{
    double a, b, c, d, n, E0;
    double P0, Pi, PL, PR, P2;
    int err = 0;

    if (tab->rows != 2 || tab->cols != 2) {
        return E_DATA;
    }

    if (tab->n > 1000) {
        /* not worth trying? */
        return E_DATA;
    }

    a = tab->f[0][0];
    b = tab->f[0][1];
    c = tab->f[1][0];
    d = tab->f[1][1];
    n = tab->n;

    E0 = (tab->rtotal[0] * tab->ctotal[0]) / n;

    /* Probability of the observed table */
    PL = PR = P2 = P0 = table_prob(a, b, c, d, n, &err);

    if (!err) {
        while (a > 0 && d > 0) {
            a -= 1; d -= 1;
            c += 1; b += 1;
            Pi = table_prob(a, b, c, d, n, &err);
            if (err) {
                break;
            }
            if (Pi <= P0 || tab->f[0][0] > E0) {
                PL += Pi;
            }
            if (Pi <= P0) {
                P2 += Pi;
            }
        }
    }

    if (!err) {
        a = tab->f[0][0];
        b = tab->f[0][1];
        c = tab->f[1][0];
        d = tab->f[1][1];

        while (c > 0 && b > 0) {
            c -= 1; b -= 1;
            a += 1; d += 1;
            Pi = table_prob(a, b, c, d, n, &err);
            if (err) {
                break;
            }
            if (Pi <= P0 || tab->f[0][0] < E0) {
                PR += Pi;
            }
            if (Pi <= P0) {
                P2 += Pi;
            }
        }
    }

    if (!err) {
        pprintf(prn, "\n%s:\n", _("Fisher's Exact Test"));
        pprintf(prn, _("  Left:   P-value = %g\n"), PL);
        pprintf(prn, _("  Right:  P-value = %g\n"), PR);
        pprintf(prn, _("  2-Tail: P-value = %g\n"), P2);
        pputc(prn, '\n');
    }

    return err;
}

/**
 * interval_model:
 * @list: high/low (2 variables) plus list of regressors.
 * @dset: dataset struct.
 * @opt: if includes %OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix; if includes %OPT_V give verbose
 * operation.
 * @prn: printing struct in case additional information is
 * wanted (signalled via %OPT_V).
 *
 * Returns: a #MODEL struct, containing interval estimates of the
 * model specified by @list.
 */

MODEL interval_model (const int *list, DATASET *dset,
                      gretlopt opt, PRN *prn)
{
    MODEL intmod;
    MODEL (* interval_estimate) (const int *, DATASET *, gretlopt,
                                 PRN *);

    gretl_error_clear();

    interval_estimate = get_plugin_function("interval_estimate");

    if (interval_estimate == NULL) {
        gretl_model_init(&intmod, dset);
        intmod.errcode = E_FOPEN;
        return intmod;
    }

    intmod = (*interval_estimate) (list, dset, opt, prn);
    set_model_id(&intmod, opt);

    return intmod;
}

/**
 * tobit_model:
 * @list: dependent variable plus list of regressors.
 * @llim: left bound on dependent variable; use #NADBL for
 * no left-censoring.
 * @rlim: right bound on dependent variable; use #NADBL for
 * no right-censoring.
 * @dset: dataset struct.
 * @opt: may include OPT_V for verbose operation, OPT_R for
 * robust (QML) standard errors.
 * @prn: printing struct for iteration info (or NULL if this is not
 * wanted).
 *
 * Produce Tobit estimates of the model given in @list.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tobit_model (const int *list, double llim, double rlim,
                   DATASET *dset, gretlopt opt, PRN *prn)
{
    MODEL (* tobit_estimate) (const int *, double, double,
                              DATASET *, gretlopt, PRN *);
    MODEL tmod;

    gretl_error_clear();

    tobit_estimate = get_plugin_function("tobit_via_intreg");
    if (tobit_estimate == NULL) {
        gretl_model_init(&tmod, dset);
        tmod.errcode = E_FOPEN;
        return tmod;
    }

    tmod = (*tobit_estimate) (list, llim, rlim, dset, opt, prn);
    set_model_id(&tmod, opt);

    return tmod;
}

/* run several checks on the data supplied for a duration
   model, before invoking the duration plugin to complete
   the business
*/

static int duration_precheck (const int *list,
                              DATASET *dset, MODEL *pmod,
                              int *pcensvar)
{
    int *olslist = NULL;
    int seppos, censvar = 0;
    int l0 = list[0];
    int err = 0;

    /* such models must contain a constant */
    if (!gretl_list_const_pos(list, 2, dset)) {
        return E_NOCONST;
    }

    /* if there's a separator, it must be in second-last place */
    seppos = gretl_list_separator_position(list);
    if (seppos > 0 && seppos != l0 - 1) {
        return E_PARSE;
    }

    if (seppos) {
        /* the censoring variable, if present, must be a dummy */
        censvar = list[l0];
        if (!gretl_isdummy(dset->t1, dset->t2, dset->Z[censvar])) {
            gretl_errmsg_sprintf(_("The variable '%s' is not a 0/1 variable."),
                                 dset->varname[censvar]);
            err = E_DATA;
        } else {
            olslist = gretl_list_copy(list);
            if (olslist == NULL) {
                err = E_ALLOC;
            } else {
                /* include the censoring dummy. to ensure the
                   sample is right */
                olslist[l0 - 1] = censvar;
                olslist[0] -= 1;
            }
        }
    }

    if (!err) {
        /* run an initial OLS to "set the model up" and check for errors;
           the duration_estimate_driver function will overwrite the
           coefficients etc.
        */
        if (olslist != NULL) {
            *pmod = lsq(olslist, dset, OLS, OPT_A);
            if (!pmod->errcode) {
                /* remove reference to censoring var */
                pmod->list[0] -= 1;
                pmod->ncoeff -= 1;
                pmod->dfn -= 1;
                pmod->dfd += 1;
            }
            free(olslist);
        } else {
            *pmod = lsq(list, dset, OLS, OPT_A);
        }
        err = pmod->errcode;
        *pcensvar = censvar;
    }

    if (!err) {
        int t, yno = pmod->list[1];

        for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t]) && dset->Z[yno][t] <= 0) {
                gretl_errmsg_set(_("Durations must be positive"));
                err = E_DATA;
            }
        }
    }

    return err;
}

/**
 * duration_model:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: may include OPT_R for robust covariance matrix.
 * @prn: printing struct for iteration info (or NULL is this is not
 * wanted).
 *
 * Estimate the duration model given in @list using ML.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL duration_model (const int *list, DATASET *dset,
                      gretlopt opt, PRN *prn)
{
    MODEL dmod;
    int censvar = 0;
    int (* duration_estimate) (MODEL *, int, const DATASET *,
                               gretlopt, PRN *);

    gretl_error_clear();
    gretl_model_init(&dmod, dset);

    dmod.errcode = duration_precheck(list, dset, &dmod,
                                     &censvar);
    if (dmod.errcode) {
        return dmod;
    }

    duration_estimate = get_plugin_function("duration_estimate");

    if (duration_estimate == NULL) {
        dmod.errcode = E_FOPEN;
        return dmod;
    }

    (*duration_estimate) (&dmod, censvar, dset, opt, prn);
    set_model_id(&dmod, opt);

    return dmod;
}

static int get_trailing_var (int *list)
{
    int l0 = list[0];
    int ret = 0;

    if (list[l0 - 1] == LISTSEP) {
        ret = list[l0];
        list[0] -= 2;
    }

    return ret;
}

/**
 * count_model:
 * @list: dependent variable plus list of regressors.
 * @ci: either POISSON or NEGBIN.
 * @dset: dataset struct.
 * @opt: may include OPT_R for robust covariance matrix.
 * @prn: printing struct for iteration info (or NULL is this is not
 * wanted).
 *
 * Estimate the count data model given in @list using ML.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL count_model (const int *list, int ci, DATASET *dset,
                   gretlopt opt, PRN *prn)
{
    MODEL cmod;
    int *listcpy;
    int offvar;
    int (* count_data_estimate) (MODEL *, int, int,
                                 DATASET *, gretlopt,
                                 PRN *);

    gretl_error_clear();

    gretl_model_init(&cmod, dset);

    if (!gretl_iscount(dset->t1, dset->t2, dset->Z[list[1]])) {
        gretl_errmsg_sprintf(_("%s: the dependent variable must be count data"),
                             gretl_command_word(ci));
        cmod.errcode = E_DATA;
        return cmod;
    }

    listcpy = gretl_list_copy(list);
    if (listcpy == NULL) {
        cmod.errcode = E_ALLOC;
        return cmod;
    }

    offvar = get_trailing_var(listcpy);

    /* run an initial OLS to "set the model up" and check for errors.
       the count_data_estimate_driver function will overwrite the
       coefficients etc.
    */

    cmod = lsq(listcpy, dset, OLS, OPT_A);
    free(listcpy);

    if (cmod.errcode) {
        return cmod;
    }

    count_data_estimate = get_plugin_function("count_data_estimate");
    if (count_data_estimate == NULL) {
        cmod.errcode = E_FOPEN;
        return cmod;
    }

    (*count_data_estimate) (&cmod, ci, offvar, dset, opt, prn);
    set_model_id(&cmod, opt);

    return cmod;
}

/**
 * heckit_model:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: option flags (may include OPT_V for verbose output).
 * @prn: printing struct for iteration info (or NULL is this is not
 * wanted).
 *
 * Produce Heckit estimates of the model given in @list. The list must
 * include a separator to divide the main equation from the selection
 * equation.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL heckit_model (const int *list, DATASET *dset,
                    gretlopt opt, PRN *prn)
{
    MODEL model;
    MODEL (* heckit_estimate) (const int *, DATASET *,
                               gretlopt, PRN *);

    gretl_error_clear();

    heckit_estimate = get_plugin_function("heckit_estimate");
    if (heckit_estimate == NULL) {
        gretl_model_init(&model, dset);
        model.errcode = E_FOPEN;
        return model;
    }

    model = (*heckit_estimate) (list, dset, opt, prn);
    set_model_id(&model, opt);

    return model;
}

/**
 * reprobit_model:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: option flags (may include OPT_V for verbose output).
 * @prn: printing struct for iteration info (or NULL if this is not
 * wanted).
 *
 * Produce random-effects probit estimates of the model given in @list.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL reprobit_model (const int *list, DATASET *dset,
                      gretlopt opt, PRN *prn)
{
    MODEL model;
    MODEL (* reprobit_estimate) (const int *, DATASET *,
                                 gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();

    if (!dataset_is_panel(dset)) {
        err = E_PDWRONG;
    } else {
        reprobit_estimate = get_plugin_function("reprobit_estimate");
        if (reprobit_estimate == NULL) {
            err = E_FOPEN;
        }
    }

    if (err) {
        gretl_model_init(&model, dset);
        model.errcode = err;
        return model;
    }

    model = (*reprobit_estimate) (list, dset, opt, prn);
    set_model_id(&model, opt);

    return model;
}

/* apparatus for computing all permutations of @k 1s in @n bits */

/* The following clever algorithm to calculate the "next permutation"
   of bits in lexicographic order is attributed to Dario Sneidermanis
   at https://stackoverflow.com/questions/1851134/
   generate-all-binary-strings-of-length-n-with-k-bits-set
*/

static guint64 next_perm (guint64 v)
{
    guint64 w; /* next permutation of bits */

    /* t gets v's least significant 0 bits set to 1 */
    guint64 t = v | (v - 1LL);
    /* Next set to 1 the most significant bit to change, set to 0 the
       least significant ones, and add the necessary 1 bits.
    */
    w = (t + 1LL) | (((~t & -~t) - 1LL) >> (__builtin_ctzll(v) + 1LL));

    return w;
}

static void binrow (int i, int n, guint64 val, gretl_matrix *m)
{
    int k, j = 0;

    for (k=n-1; k>=0; k--) {
        if (val & (1LL << k)) {
            gretl_matrix_set(m, i, j, 1);
        }
        j++;
    }
}

static int nperms (int n, int k)
{
    double num = lngamma(n+1);
    double den1 = lngamma(k+1);
    double den2 = lngamma(n-k+1);

    return (int) nearbyint(exp(num - (den1 + den2)));
}

gretl_matrix *bit_permutations (int n, int k, int *err)
{
    gretl_matrix *ret;
    guint64 v = 0, vmax = 0;
    int i, np = 0;

    if (n > 64) {
        gretl_errmsg_set(_("binperms: n must be <= 64"));
        *err = E_INVARG;
        return NULL;
    }

    if (n < 0 || k < 0 || n < k) {
	gretl_errmsg_set(_("binperms: we need n >= k >= 0"));
	*err = E_INVARG;
	return NULL;
    } else if (n == 0 && k == 0) {
	ret = gretl_zero_matrix_new(1, 0);
    } else if (k == n) {
        ret = gretl_unit_matrix_new(1, n);
    } else if (k == 0) {
        ret = gretl_zero_matrix_new(1, n);
    } else if (k == 1) {
        ret = gretl_zero_matrix_new(n, n);
        if (ret != NULL) {
            for (i = n-1; i<n*n-1; i += n-1) {
                ret->val[i] = 1.0;
            }
        }
    } else if (k == n - 1) {
        ret = gretl_unit_matrix_new(n, n);
        if (ret != NULL) {
            for (i=0; i<n; i++) {
                gretl_matrix_set(ret, i, i, 0);
            }
        }
    } else {
        np = nperms(n, k);
        ret = gretl_zero_matrix_new(np, n);
    }

    if (ret == NULL) {
        *err = E_ALLOC;
        return NULL;
    } else if (np == 0) {
        /* result already settled */
        return ret;
    }

    /* set the min and max of @v */
    for (i=0; i<k; i++) {
        v |= (1LL << i);
        vmax |= (1LL << (n-1-i));
    }

    /* run the iteration */
    for (i=0; i<np; i++) {
        binrow(i, n, v, ret);
	if (i < np - 1) {
            v = next_perm(v);
        }
    }

    if (v != vmax) {
	/* "can't happen" */
	fprintf(stderr, "binperms: something went wrong!\n");
    }

    return ret;
}
