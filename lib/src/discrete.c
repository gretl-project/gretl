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

#define CHOL_TINY 1.0e-13

typedef struct op_container_ op_container;

/* structure for handling ordered probit or logit */

struct op_container_ {
    int ci;           /* model command index (PROBIT or LOGIT) */
    gretlopt opt;     /* option flags */
    int *y;           /* dependent variable */
    double **Z;       /* data */
    int *list;        /* dependent var plus regular regressors */
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
    gretl_matrix_free(OC->G);
    free(OC->g);
    free(OC->theta);

    free(OC);
}

static op_container *op_container_new (int ci, int ndum,
				       double **Z, MODEL *pmod,  
				       gretlopt opt)
{
    op_container *OC;
    int i, t, vy = pmod->list[1];
    int nobs = pmod->nobs;

    OC = malloc(sizeof *OC);
    if (OC == NULL) {
	return NULL;
    }

    OC->ci = ci;

    OC->Z = Z;
    OC->pmod = pmod;
    OC->t1 = pmod->t1;
    OC->t2 = pmod->t2;
    OC->nobs = nobs;
    OC->k = pmod->ncoeff;
    OC->ymax = ndum;
    OC->nx = OC->k - ndum;

    OC->opt = opt;

    OC->y = NULL;
    OC->ndx = NULL;
    OC->dP = NULL;
    OC->list = NULL;
    OC->G = NULL;
    OC->g = NULL;

    OC->y = malloc(nobs * sizeof *OC->y);
    OC->ndx = malloc(nobs * sizeof *OC->ndx);
    OC->dP = malloc(nobs * sizeof *OC->dP);

    OC->list = gretl_list_new(1 + OC->nx);
    OC->G = gretl_matrix_alloc(nobs, OC->k);
    OC->g = malloc(OC->k * sizeof *OC->g);
    OC->theta = malloc(OC->k * sizeof *OC->theta);

    if (OC->y == NULL || OC->ndx == NULL || 
	OC->dP == NULL || OC->list == NULL ||
	OC->G == NULL || OC->g == NULL ||
	OC->theta == NULL) {
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

#if LPDEBUG
    fprintf(stderr, "nobs = %d\n", OC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", OC->t1, OC->t2);
    fprintf(stderr, "k = %d\n", OC->k);
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
		t, OC->nobs, s, yt, OC->ndx[s], ystar0, ystar1);
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
	    fprintf(stderr, "very small dP at obs %d; y=%d, ndx=%g, dP=%g\n", 
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
   points are represented to BFGS in the form of the log-difference
   from the previous cut point.
*/

static void op_make_BFGS_theta (op_container *OC, double *theta)
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
   op_make_BFGS_theta() */

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

    return 1;
}

static int ordered_hessian (op_container *OC, gretl_matrix *H) 
{
    double smal = 1.0e-07;  /* "small" is some sort of macro on win32 */
    double smal2 = 2.0 * smal;
    double ti, x, ll, *g0;
    int i, j, k = OC->k;
    int err = 0;

    g0 = malloc(k * sizeof *g0);
    if (g0 == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	ti = OC->theta[i];
	OC->theta[i] -= smal;
	ll = op_loglik(OC->theta, OC);
	if (na(ll)) {
	    OC->theta[i] = ti;
	    err = E_DATA;
	    break;
	}
	for (j=0; j<k; j++) {
	    g0[j] = OC->g[j];
	}
	OC->theta[i] += smal2;
	ll = op_loglik(OC->theta, OC);
	if (na(ll)) {
	    OC->theta[i] = ti; 
	    err = E_DATA;
	    break;
	}	
	for (j=0; j<k; j++) {
	    x = (OC->g[j] - g0[j]) / smal2;
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

    return H;
}

/**
 * ordered_model_prediction:
 * @pmod: model for ordered data, either logit or probit.
 * @Xb: X\beta, the value of the index function at a given
 * observation.
 *
 * Returns: the predicted value of the (ordinal) dependent variable;
 * that is, the value for which the estimated probability is greatest.
 */

double ordered_model_prediction (const MODEL *pmod, double Xb)
{
    /* position of least cut point in coeff array */
    int k = gretl_model_get_int(pmod, "nx");
    int maxval = pmod->ncoeff - k;
    double prob, pmax, cut;
    double CDF, CDFbak;
    int i, pred = 0;

    cut = pmod->coeff[k];
    pmax = CDFbak = lp_cdf(cut - Xb, pmod->ci);

    for (i=1; i<maxval; i++) {
	cut = pmod->coeff[++k];
	CDF = lp_cdf(cut - Xb, pmod->ci);
	prob = CDF - CDFbak;
	if (prob > pmax) {
	    pmax = prob;
	    pred = i;
	}
	CDFbak = CDF;
    }

    prob = 1 - CDFbak;
    if (prob > pmax) {
	pred = maxval;
    }

    return (double) pred;
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
	/* L'Hopital-based approximation */
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

static void cut_points_init (op_container *OC, const MODEL *pmod, 
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

static void op_LR_test (MODEL *pmod, op_container *OC, 
			const double **Z)
{
    int nx = OC->nx;
    double L0;

    OC->k -= OC->nx;
    OC->nx = 0;

    cut_points_init(OC, pmod, Z);

    L0 = op_loglik(OC->theta, OC);

    if (!na(L0) && L0 <= pmod->lnL) {
	pmod->chisq = 2.0 * (pmod->lnL - L0);
	gretl_model_set_int(pmod, "lr_df", nx);
    }

    /* restore original data on OC */
    OC->nx = nx;
    OC->k += nx;
}

static int oprobit_normtest (MODEL *pmod, op_container *OC)
{
    int nobs = OC->nobs;
    int k = OC->k;
    int nx = OC->nx;
    int M = OC->ymax;
    double *theta = OC->theta;
    int t, s, i, yt, err = 0;
    gretl_matrix *CMtestmat;
    gretl_matrix *y;
    gretl_matrix *beta;
    double m0, m1, u, v, a2v, b2u;
    double gval, a = 0, b = 0;
    double e3, e4;

    CMtestmat = gretl_matrix_alloc(nobs, k+2);
    y = gretl_unit_matrix_new(nobs, 1);
    beta = gretl_matrix_alloc(k+2, 1);

    if (CMtestmat == NULL || y == NULL || beta == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    s = 0;
    for (t=OC->pmod->t1; t<=OC->pmod->t2; t++) {
	if (na(OC->pmod->uhat[t])) {
	    continue;
	}

	yt = OC->y[s];

	if (yt == 0) {
	    m0 = theta[nx];
	    b = OC->ndx[s] + m0;
	} else {
	    m0 = theta[nx + yt - 1];
	    a = OC->ndx[s] + m0;
	    if (yt < M) {
		m1 = theta[nx + yt];
		b = OC->ndx[s] + m1;
	    }
	} 

	if (yt == 0) {
	    u = gretl_matrix_get(OC->G, s, nx);
	    b2u = b*b*u;
	    v = a2v = 0;
	} else {
	    gval = gretl_matrix_get(OC->G, s, nx + yt - 1);
	    v = -gval;
	    a2v = a*a*v;
	    if (yt < M) {
		u = gretl_matrix_get(OC->G, s, nx + yt);
		b2u = b*b*u;
	    } else {
		b2u = u = 0;
	    }
	} 

	for (i=0; i<k; i++) {
	    gval = gretl_matrix_get(OC->G, s, i);
	    gretl_matrix_set(CMtestmat, s, i, gval);
	}

	e3 = 2*(v-u) + (a2v - b2u);
	e4 = 3*(a*v-b*u) + (a*a2v - b*b2u);

	gretl_matrix_set(CMtestmat, s, k, e3);
	gretl_matrix_set(CMtestmat, s, k+1, e4);

	s++;
    }

    err = gretl_matrix_ols(y, CMtestmat, beta, NULL, NULL, NULL);

    if (!err) {
	double X2 = nobs;

	gretl_matrix_multiply(CMtestmat, beta, y);
	for (t=0; t<y->rows; t++) {
	    X2 -= (1 - y->val[t]) * (1 - y->val[t]);
	}

	if (X2 > 0) {
	    gretl_model_add_normality_test(pmod, X2);
	}
    }

 bailout:

    gretl_matrix_free(CMtestmat);
    gretl_matrix_free(y);
    gretl_matrix_free(beta);

    return err;
}

static void fill_op_model (MODEL *pmod, const int *list,
			   const DATASET *dset, 
			   op_container *OC,
			   int fncount, int grcount)
{
    gretl_matrix *H = NULL;
    int npar = OC->k;
    int nx = OC->nx;
    int correct = 0;
    double Xb;
    int i, s, t, v;
    int err = 0;

    H = ordered_hessian_inverse(OC, &err);
    if (err) {
	goto bailout;
    }

    if (OC->opt & OPT_R) {
	err = gretl_model_add_QML_vcv(pmod, OC->ci, H, OC->G,
				      dset, OC->opt);
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
	if (na(OC->pmod->uhat[t])) {
	    continue;
	}

	Xb = 0.0;
	for (i=0; i<OC->nx; i++) {
	    v = OC->list[i+2];
	    Xb += OC->theta[i] * OC->Z[v][t];
	}

	/* yhat = X\hat{beta} */
	pmod->yhat[t] = Xb;
	if (ordered_model_prediction(pmod, Xb) == OC->y[s]) {
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
}

/* Main ordered estimation function */

static int do_ordered (int ci, int ndum, 
		       DATASET *dset, 
		       MODEL *pmod, const int *list, 
		       gretlopt opt, PRN *prn)
{
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;
    op_container *OC;
    int i, npar;
    double *theta = NULL;
    double toler;
    int use_newton = 0;
    int err;

    OC = op_container_new(ci, ndum, dset->Z, pmod, opt);
    if (OC == NULL) {
	return E_ALLOC;
    }

    npar = OC->k;

    /* transformed theta to pass to BFGS */
    theta = malloc(npar * sizeof *theta);
    if (theta == NULL) {
	op_container_destroy(OC);
	return E_ALLOC;
    }

    /* initialize slopes */
    for (i=0; i<OC->nx; i++) {
	OC->theta[i] = 0.0001;
    }

    /* initialize cut points */
    cut_points_init(OC, pmod, (const double **) dset->Z);

    /* transform theta to log-diff form */
    op_make_BFGS_theta(OC, theta);

#if LPDEBUG
    for (i=0; i<npar; i++) {
	fprintf(stderr, "theta[%d]: 'real' = %g, BFGS = %g\n", i, 
		OC->theta[i], theta[i]);
    }
    fprintf(stderr, "\ninitial loglikelihood = %.12g\n", 
	    op_loglik(theta, OC));
#endif

    if (libset_get_int(GRETL_OPTIM) == OPTIM_NEWTON) {
	use_newton = 1;
    }

    if (use_newton) {
	double crittol = 1.0e-7;
	double gradtol = 1.0e-7;

	err = newton_raphson_max(theta, npar, maxit, 
				 crittol, gradtol, &fncount, 
				 C_LOGLIK, op_loglik, 
				 op_score, NULL, OC, 
				 (prn != NULL)? OPT_V : OPT_NONE,
				 prn);
    } else {
	BFGS_defaults(&maxit, &toler, PROBIT);
	err = BFGS_max(theta, npar, maxit, toler, 
		       &fncount, &grcount, op_loglik, C_LOGLIK,
		       op_score, OC, NULL, 
		       (prn != NULL)? OPT_V : OPT_NONE, prn);
    }

    if (err) {
	goto bailout;
    } 

    /* transform back to 'real' theta */
    op_get_real_theta(OC, theta);

    if (!err) {
	fill_op_model(pmod, list, dset, OC, fncount, grcount);
    }

 bailout:

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
				double **orig_y, int *ndum)
{
    struct sorter *s = NULL;
    double *sy = NULL;
    double nexty, bady;
    int dv = pmod->list[1];
    int i, t, n = 0;
    int fixit = 0;
    int err = 0;

    /* make a copy of the observations on the dep. var. that
       were actually used in the initial OLS, recording
       the observation numbers, then sort
    */

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    n++;
	}
    }

    s = malloc(n * sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    s[i].x = dset->Z[dv][t];
	    s[i].t = t;
	    i++;
	}
    }

    qsort(s, n, sizeof *s, gretl_compare_doubles);

    /* normalize to minimum zero */
    if (s[0].x != 0) {
	double ymin = s[0].x;

	for (i=0; i<n && s[i].x == ymin; i++) {
	    s[i].x = 0.0;
	}
	fixit = 1;
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
		fixit = 1;
	    } 
	} 
    }

    /* the number of dummies actually used will equal the
       max of normalized y
    */
    *ndum = (int) s[n-1].x;

    if (fixit) {
	/* the dependent var needs transforming */
	sy = copyvec(dset->Z[dv], dset->n);
	if (sy == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		sy[s[i].t] = s[i].x;
	    }
	    /* back up the original dep. var and replace it for
	       the duration of the ordered analysis */
	    *orig_y = dset->Z[dv];
	    dset->Z[dv] = sy;
	}
    }

#if 0
    if (fixit) {
	fputs("oprobit: using normalized y\n", stderr);
    } else {
	fputs("oprobit: using original y\n", stderr);
    }
#endif

    free(s);

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

static void list_purge_const (int *list, DATASET *dset)
{
    MODEL tmpmod;
    int depvar = list[1];
    int i, ok = 0;

    /* first remove the constant itself, if present */

    for (i=2; i<=list[0]; i++) {
	if (list[i] == 0) {
	    gretl_list_delete_at_pos(list, i);
	    break;
	}
    }

    /* drop other stuff possibly collinear with the constant 
       (eg sets of dummies) */
    
    list[1] = 0; /* substitute the constant as dependent */

    while (!ok) {
	tmpmod = lsq(list, dset, OLS, OPT_A);
	ok = (tmpmod.ess > 1.0e-6);
	if (!ok) {
	    for (i=tmpmod.ncoeff-1; i>0; i--) {
		if (fabs(tmpmod.coeff[i]) > 1.0e-06) {
		    gretl_list_delete_at_pos(list, i+2);
		    break;
		}
	    }
	}
	clear_model(&tmpmod);
    }

    /* reinstate the real dependent variable */
    list[1] = depvar;
}

static int check_for_missing_dummies (MODEL *pmod, int *dlist)
{
    int i, dv;

    for (i=1; i<=dlist[0]; i++) {
	dv = dlist[i];
	if (!in_gretl_list(pmod->list, dv)) {
	    fprintf(stderr, "check for missing dummies: var %d is gone!\n", dv);
	    return E_SINGULAR;
	}
    }

    return 0;
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

/* driver function for ordered logit/probit: note, prn is non-NULL
   only if the verbose option has been selected
*/

static MODEL ordered_estimate (int *list, DATASET *dset, int ci,
			       gretlopt opt, PRN *prn) 
{
    MODEL model;
    int orig_v = dset->v;
    double *orig_y = NULL;
    int *biglist = NULL;
    int *dumlist = NULL;
    int ndum = 0;

    gretl_model_init(&model, dset);

    model.errcode = ordered_depvar_check(list[1], dset);

    if (!model.errcode) {
	/* remove the constant from the incoming list, if present */
	list_purge_const(list, dset);

	/* construct augmented regression list, including dummies 
	   for the level of the dependent variable
	*/
	biglist = make_op_list(list, dset, &dumlist, &model.errcode);
    }

    if (model.errcode) {
	return model;
    }

    /* run initial OLS, with dummies added */
    model = lsq(biglist, dset, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "ordered_estimate: initial OLS failed\n");
	free(dumlist);
	free(biglist);
	return model;
    }

#if LPDEBUG
    pprintf(prn, "ordered_estimate: initial OLS\n");
    printmodel(&model, dset, OPT_S, prn);
#endif

    if (model.list[0] < biglist[0]) {
	/* One or more regressors were dropped in OLS, most likely due
	   to collinearity.  We can accept the dropping of
	   user-selected regressors, but all the added dummies must be
	   retained, so we'd better investigate.
	*/
	model.errcode = check_for_missing_dummies(&model, dumlist);
    }
	
    if (!model.errcode) {
	/* after accounting for any missing observations, check that
	   the dependent variable is acceptable 
	*/
	model.errcode = maybe_fix_op_depvar(&model, dset,
					    &orig_y, &ndum);
    }

    /* do the actual ordered probit analysis */
    if (!model.errcode) {
	clear_model_xpx(&model);
	model.errcode = do_ordered(ci, ndum, dset, &model, list, 
				   opt, prn);
    }

    free(dumlist);
    free(biglist);

    if (orig_y != NULL) {
	/* if we messed with the dependent var, put the original back */
	restore_depvar(dset->Z, orig_y, list[1]);
    }

    if (dset->v > orig_v) {
	/* clean up any automatically-added dummies */
	dataset_drop_last_variables(dset, dset->v - orig_v);
    }

    set_model_id(&model);

    return model;
}

static double logit (double x)
{
    double l = 1.0 / (1.0 + exp(-x));

#if LPDEBUG
    if (x > 40 || x < -40) {
	fprintf(stderr, "x = %g, logit(x) = %g\n", x, l);
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

static char *classifier_check (int *list, const DATASET *dset,
			       PRN *prn, int *err)
{
    char *mask = NULL;
    const double *y = dset->Z[list[1]];
    int i, v, ni, t;

    for (i=2; i<=list[0]; i++) {
	v = list[i];
	if (v == 0) {
	    continue;
	}

	ni = gretl_isdummy(dset->t1, dset->t2, dset->Z[v]);

	if (ni > 0) {
	    int same[2] = {0};
	    int diff[2] = {0};
	    int maskval, pc = 1;

	    for (t=dset->t1; t<=dset->t2 && pc; t++) {
		if (dset->Z[v][t] > 0) {
		    if (y[t] > 0) {
			same[0] += 1;
		    } else {
			diff[0] += 1;
		    }
		} else {
		    if (y[t] == 0) {
			same[1] += 1;
		    } else {
			diff[1] += 1;
		    }
		}
		if (same[0] && diff[0] && same[1] && diff[1]) {
		    pc = 0;
		}
	    }

	    if (pc) {
		if (!(same[0] && diff[0])) {
		    maskval = 1;
		    pprintf(prn, "Note: %s != 0 predicts %s perfectly\n",
			    dset->varname[v], (same[0])? "success" : "failure");
		} else {
		    maskval = 0;
		    pprintf(prn, "Note: %s == 0 predicts %s perfectly\n",
			    dset->varname[v], (diff[1])? "success" : "failure");
		}

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
			}
		    }
		    pprintf(prn, "%d observations not used\n", ni);
		}
		gretl_list_delete_at_pos(list, i--);
	    }
	}
    }

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
    double x, xti, ll = 0.0;
    int yt, i, t;

    errno = 0;

    for (i=0; i<mnl->npar; i++) {
	mnl->b->val[i] = theta[i];
    }

    gretl_matrix_multiply(mnl->X, mnl->b, mnl->Xb);

    for (t=0; t<mnl->T && !errno; t++) {
	x = 1.0;
	for (i=0; i<mnl->n; i++) {
	    /* sum row i of exp(Xb) */
	    xti = exp(gretl_matrix_get(mnl->Xb, t, i));
	    x += xti;
	}
	ll -= log(x);
	yt = gretl_vector_get(mnl->y, t);
	if (yt > 0) {
	    ll += gretl_matrix_get(mnl->Xb, t, yt-1);
	}
    }

    if (errno != 0) {
	ll = NADBL;
    }

    return ll;
}

static int mn_logit_score (double *theta, double *s, int npar, 
			   BFGS_CRIT_FUNC ll, void *ptr)
{
    mnl_info *mnl = (mnl_info *) ptr;
    double x, xti, pti, p, g;
    int i, j, k, t, yt;
    int err = 0;

    errno = 0;

    for (i=0; i<npar; i++) {
	s[i] = 0.0;
    }

    for (t=0; t<mnl->T && !errno; t++) {
	x = 1.0;
	for (i=0; i<mnl->n; i++) {
	    /* sum row i of exp(Xb) */
	    xti = exp(gretl_matrix_get(mnl->Xb, t, i));
	    x += xti;
	    gretl_matrix_set(mnl->P, t, i, xti);
	}

	if (!errno) {
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
    }

    if (errno != 0) {
	err = E_NAN;
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
					  dset, opt);
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
 * observation the probabiity is set to NaN; provided the
 * regressor information is complete we compute the
 * outcome probabilities even if the actual outcome is
 * missing.
 *
 * Returns: allocated matrix or NULL on failure.
 */

gretl_matrix *mn_logit_probabilities (const MODEL *pmod,
				      const DATASET *dset,
				      int *err)
{
    gretl_matrix *P = NULL;
    const gretl_matrix *yvals = NULL;
    double *eXbt = NULL;
    int i, m = 0;

    if (pmod == NULL || pmod->list == NULL || pmod->coeff == NULL) {
	*err = E_DATA;
    }

    if (!*err) {
	/* list of outcome values (including the base case) */
	yvals = gretl_model_get_data(pmod, "yvals");
	if (yvals == NULL) {
	    *err = E_DATA;
	} else {
	    m = gretl_vector_get_length(yvals);
	}
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
	int n = pmod->t2 - pmod->t1 + 1;

	P = gretl_matrix_alloc(n, m);
	if (P == NULL) {
	    *err = E_ALLOC;
	} else {
	    /* allow for casting results to series */
	    gretl_matrix_set_t1(P, pmod->t1);
	    gretl_matrix_set_t2(P, pmod->t2);
	}
    }

    if (!*err) {
	/* allocate required workspace */
	eXbt = malloc(m * sizeof *eXbt);
	if (eXbt == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	const double *b = pmod->coeff;
	double St, ptj;
	int j, k, t, vi, ok;

	for (t=pmod->t1; t<=pmod->t2; t++) {
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
		    gretl_matrix_set(P, t, j, M_NA);
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
		    gretl_matrix_set(P, t, j, ptj);
		}
	    }
	}
    }

    free(eXbt);

    return P;
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
	    gretl_errmsg_set("logit: the dependent variable must form a sequence of "
			     "integers");
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
	set_model_id(pmod);
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

    if (use_bfgs) {
	mod.errcode = BFGS_max(mnl->theta, mnl->npar, maxit, 0.0, 
			       &fncount, &grcount, mn_logit_loglik, C_LOGLIK,
			       mn_logit_score, mnl, NULL,
			       (prn != NULL)? OPT_V : OPT_NONE, prn);
    } else {	
	double crittol = 1.0e-8;
	double gradtol = 1.0e-7;
	
	maxit = 100;
	mod.errcode = newton_raphson_max(mnl->theta, mnl->npar, maxit, 
					 crittol, gradtol, &fncount, 
					 C_LOGLIK, mn_logit_loglik, 
					 mn_logit_score, mnl_hessian, mnl,
					 (prn != NULL)? OPT_V : OPT_NONE,
					 prn);
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

MODEL biprobit_model (int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    MODEL bpmod;
    void *handle;
    MODEL (* biprobit_estimate) (const int *, DATASET *, 
				 gretlopt, PRN *);

    gretl_error_clear();

    biprobit_estimate = get_plugin_function("biprobit_estimate", &handle);

    if (biprobit_estimate == NULL) {
	gretl_model_init(&bpmod, dset);
	bpmod.errcode = E_FOPEN;
	return bpmod;
    }

    bpmod = (*biprobit_estimate) (list, dset, opt, prn);

    close_plugin(handle);

    set_model_id(&bpmod);

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
    gretl_matrix_block *B;
    gretl_matrix *X;  /* regressors */
    gretl_matrix *pX; /* for use with Hessian */
    gretl_matrix *b;  /* coefficients in matrix form */
    gretl_matrix *Xb; /* index function values */
};

static void bin_info_destroy (bin_info *bin)
{
    if (bin != NULL) {
	gretl_matrix_block_destroy(bin->B);
	free(bin->theta);
	free(bin->y);
	free(bin);
    }
}

static bin_info *bin_info_new (int ci, int k, int T)
{
    bin_info *bin = malloc(sizeof *bin);

    if (bin != NULL) {
	bin->ci = ci;
	bin->k = k;
	bin->T = T;
	bin->pp_err = 0;
	bin->theta = malloc(k * sizeof *bin->theta);
	if (bin->theta == NULL) {
	    free(bin);
	    return NULL;
	}
	bin->y = malloc(T * sizeof *bin->y);
	if (bin->y == NULL) {
	    free(bin->theta);
	    free(bin);
	    return NULL;
	}	
	bin->B = gretl_matrix_block_new(&bin->X, T, k,
					&bin->pX, T, k,
					&bin->b, k, 1,
					&bin->Xb, T, 1,
					NULL);
	if (bin->B == NULL) {
	    free(bin->theta);
	    free(bin->y);
	    free(bin);
	    bin = NULL;
	} 
    }

    return bin;
}

/*
  If min1 > max0, then there exists a separating hyperplane between
  all the zeros and all the ones; in this case no MLE exists.
*/

static int perfect_prediction_check (bin_info *bin)
{
    double max0 = -1.0e200;
    double min1 = 1.0e200;
    const double *ndx = bin->Xb->val;
    int t;

    for (t=0; t<bin->T; t++) {
	if (bin->y[t] == 0 && ndx[t] > max0) {
	    max0 = ndx[t];
	} else if (bin->y[t] == 1 && ndx[t] < min1) {
	    min1 = ndx[t];
	}
    }

    return (min1 > max0);
}

/* compute loglikelihood for binary probit/logit */

static double binary_loglik (const double *theta, void *ptr)
{
    bin_info *bin = (bin_info *) ptr;
    double e, ndx, p, ll = 0.0;
    int yt, i, t;

    errno = 0;

    for (i=0; i<bin->k; i++) {
	bin->b->val[i] = theta[i];
    }

    gretl_matrix_multiply(bin->X, bin->b, bin->Xb);

    if (perfect_prediction_check(bin)) {
	bin->pp_err = 1;
	return NADBL;
    }

    for (t=0; t<bin->T && !errno; t++) {
	yt = bin->y[t];
	ndx = gretl_vector_get(bin->Xb, t);
	if (bin->ci == PROBIT) {
	    p = yt ? normal_cdf(ndx) : normal_cdf(-ndx);
	} else {
	    e = logit(ndx);
	    p = yt ? e : 1-e;
	}
	ll += log(p);
    }

    if (errno != 0) {
	ll = NADBL;
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

    errno = 0;

    for (j=0; j<bin->k; j++) {
	s[j] = 0.0;
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
	    s[j] += w * gretl_matrix_get(bin->X, t, j);
	}
    }

    if (errno != 0) {
	err = E_NAN;
    }

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
    if (err) {
	return err;
    }

    if (opt & OPT_R) {
	G = binary_score_matrix(bin, &err);
    }

    if (!err) {
	if (opt & OPT_R) {
	    err = gretl_model_add_QML_vcv(pmod, bin->ci, H, G, 
					  dset, opt);
	} else {
	    err = gretl_model_add_hessian_vcv(pmod, H);
	}
    }

    gretl_matrix_free(H);
    gretl_matrix_free(G);

    return err;
}

static void binary_model_chisq (bin_info *bin, MODEL *pmod)
{
    int t, zeros, ones = 0;
    double Lr, chisq;

    if (pmod->ncoeff == 1 && pmod->ifc) {
	/* constant-only model */
	pmod->chisq = NADBL;
	pmod->rsq = 0;
	pmod->adjrsq = NADBL;
	return;
    }
    
    for (t=0; t<bin->T; t++) {
	ones += bin->y[t];
    }

    zeros = bin->T - ones;

    Lr = ones * log(ones / (double) bin->T);
    Lr += zeros * log(zeros / (double) bin->T);
    chisq = 2.0 * (pmod->lnL - Lr);

    if (chisq < 0) {
	pmod->rsq = pmod->adjrsq = pmod->chisq = NADBL;
    } else {
	pmod->chisq = chisq;
	/* McFadden pseudo-R^2 */
	pmod->rsq = 1.0 - pmod->lnL / Lr;
	pmod->adjrsq = 1.0 - (pmod->lnL - bin->k) / Lr;
    }
}

/* Binary probit normality test as in Jarque, Bera and Lee (1984),
   also quoted in Verbeek, chapter 7: we regress a column of 1s on the
   products of the generalized residual with X, (X\beta)^2 and
   (X\beta)^3.  The test statistic is T times the uncentered
   R-squared, and is distributed as chi-square(2). It can be shown
   that this test is numerically identical to the Chesher-Irish (87)
   test in the probit case (although C&I make no mention of this in
   their article).  
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

static double dumslope (MODEL *pmod, const gretl_matrix *xbar, int j)
{
    double s, Xb = 0.0;
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	if (i != j) {
	    Xb += pmod->coeff[i] * xbar->val[i];
	}
    } 

    if (pmod->ci == LOGIT) {
	s = logit(Xb + pmod->coeff[j]) - logit(Xb);
    } else {
	s = normal_cdf(Xb + pmod->coeff[j]) - normal_cdf(Xb);
    }

    return s;
}

static int binary_model_add_slopes (MODEL *pmod, bin_info *bin)
{
    gretl_matrix *xbar;
    const double *Xi;
    double *slopes;
    double Xb, fXb;
    size_t ssize;
    int i, err = 0;

    xbar = gretl_matrix_column_mean(bin->X, &err);
    if (err) {
	return err;
    }    

    ssize = pmod->ncoeff * sizeof *slopes;
    slopes = malloc(ssize);
    if (slopes == NULL) {
	gretl_matrix_free(xbar);
	return E_ALLOC;
    }

    Xb = 0.0;
    for (i=0; i<pmod->ncoeff; i++) {
	Xb += pmod->coeff[i] * xbar->val[i];
    }

    fXb = (bin->ci == LOGIT)? logit_pdf(Xb) : normal_pdf(Xb);
    gretl_model_set_double(pmod, "fXb", fXb);

    Xi = bin->X->val;

    for (i=0; i<bin->k; i++) {
	if (pmod->list[i+2] == 0) {
	    slopes[i] = 0.0;
	} else if (gretl_isdummy(0, bin->T-1, Xi)) {
	    slopes[i] = dumslope(pmod, xbar, i);
	} else {
	    slopes[i] = pmod->coeff[i] * fXb;
	}
	Xi += bin->T;
    }

    err = gretl_model_set_data(pmod, "slopes", slopes, 
			       GRETL_TYPE_DOUBLE_ARRAY,
			       ssize);

    gretl_matrix_free(xbar);
    if (err) {
	free(slopes);
    }

    return err;
}

static double binary_model_fXb (bin_info *bin)
{
    double xbar, Xb = 0.0;
    int i, t;

    for (i=0; i<bin->k; i++) {
	xbar = 0.0;
	for (t=0; t<bin->T; t++) {
	    xbar += gretl_matrix_get(bin->X, t, i);
	}
	xbar /= bin->T;
	Xb += bin->theta[i] * xbar;
    }

    return (bin->ci == LOGIT)? logit_pdf(Xb) : normal_pdf(Xb);    
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

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
	double ndxt;
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
	    F = exp(ndxt) / (1.0 + exp(ndxt));
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

static int binary_model_finish (bin_info *bin, MODEL *pmod,
				const DATASET *dset,
				gretlopt opt)
{
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = bin->theta[i];
    }

    pmod->ci = bin->ci;
    pmod->lnL = binary_loglik(pmod->coeff, bin);

    if (pmod->ci == PROBIT && (opt & OPT_X)) {
	/* this model is just the starting-point for
	   random-effects probit estimation */
	pmod->opt |= OPT_P;
	return 0;
    }

    binary_model_chisq(bin, pmod);
    pmod->errcode = binary_variance_matrix(pmod, bin, dset, opt);

    if (!pmod->errcode) {
	if (opt & OPT_P) {
	    /* showing p-values, not slopes */
	    double fXb = binary_model_fXb(bin);

	    pmod->opt |= OPT_P;
	    gretl_model_set_double(pmod, "fXb", fXb);
	} else {
	    pmod->errcode = binary_model_add_slopes(pmod, bin);
	}
    }

    if (!pmod->errcode) {
	binary_model_hatvars(pmod, bin->Xb, bin->y, opt);
	if (pmod->ci == PROBIT) {
	    binary_probit_normtest(pmod, bin);
	}
	mle_criteria(pmod, 0);
	if (opt & OPT_A) {
	    pmod->aux = AUX_AUX;
	} else {
	    pmod->ID = model_count_plus();
	}
    }

    return pmod->errcode;
}

static MODEL binary_model (int ci, const int *inlist, 
			   DATASET *dset, gretlopt opt, 
			   PRN *prn)
{
    gretlopt max_opt = OPT_NONE;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *list = NULL;
    char *mask = NULL;
    double crittol = 1.0e-8;
    double gradtol = 1.0e-7;
    int maxit = 100;
    int fncount = 0;
    MODEL mod;
    bin_info *bin = NULL;
    PRN *vprn = NULL;
    int depvar;
    int i, vi, t, s;

    gretl_model_init(&mod, dset);
    depvar = inlist[1];

    if (!gretl_isdummy(dset->t1, dset->t2, dset->Z[depvar])) {
	gretl_errmsg_sprintf(_("'%s' is not a binary variable"), 
			     dset->varname[depvar]);
	mod.errcode = E_DATA;
	return mod;
    }

    list = gretl_list_copy(inlist);
    if (list == NULL) {
	mod.errcode = E_ALLOC;
	goto bailout;
    }

    /* the cluster option implies robust */
    if (opt & OPT_C) {
	opt |= OPT_R;
    }

    list_adjust_sample(list, &dset->t1, &dset->t2, dset, NULL);

    mask = classifier_check(list, dset, prn, &mod.errcode);
    if (mod.errcode) {
	goto bailout;
    }

    if (mask != NULL) {
	mod.errcode = copy_to_reference_missmask(mask);
    }

    if (!mod.errcode) {
	gretlopt ols_opt = OPT_A;

	/* If we're doing logit/probit as an auxiliary regression,
	   it might be safer to abort on perfect collinearity;
	   otherwise we'll try automatic elimination.
	*/
	if (opt & OPT_A) {
	    ols_opt |= OPT_Z;
	}
	mod = lsq(list, dset, OLS, ols_opt);
#if LPDEBUG
	printmodel(&mod, dset, OPT_NONE, prn);
#endif
    }

    if (mod.errcode) {
	goto bailout;
    }

    bin = bin_info_new(ci, mod.ncoeff, mod.nobs);
    if (bin == NULL) {
	mod.errcode = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<bin->k; i++) {
	bin->theta[i] = mod.coeff[i] / mod.sigma;
    }

    s = 0;
    for (t=mod.t1; t<=mod.t2; t++) {
	if (!na(mod.yhat[t])) {
	    bin->y[s++] = (dset->Z[depvar][t] != 0);
	}
    }

    for (i=0; i<bin->k; i++) {
	vi = mod.list[i+2];
	s = 0;
	for (t=mod.t1; t<=mod.t2; t++) {
	    if (!na(mod.yhat[t])) {
		gretl_matrix_set(bin->X, s++, i, dset->Z[vi][t]);
	    }
	}
    } 

    if (opt & OPT_V) {
	max_opt = OPT_V;
	vprn = prn;
    }

    mod.errcode = newton_raphson_max(bin->theta, bin->k, maxit, 
				     crittol, gradtol, &fncount, 
				     C_LOGLIK, binary_loglik, 
				     binary_score, binary_hessian, bin,
				     max_opt, vprn);

    if (bin->pp_err) {
	/* trash any existing error message */
	gretl_error_clear();
	mod.errcode = E_NOCONV;
	gretl_errmsg_set(_("Perfect prediction obtained: no MLE exists"));
    }

    if (!mod.errcode) {
	binary_model_finish(bin, &mod, dset, opt);
    }  

 bailout:  

    bin_info_destroy(bin);

    free(list);
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

MODEL ordered_logit (int *list, DATASET *dset, 
		     gretlopt opt, PRN *prn)
{
    PRN *vprn = (opt & OPT_V)? prn : NULL;

    return ordered_estimate(list, dset, LOGIT, opt, vprn);
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

MODEL ordered_probit (int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    PRN *vprn = (opt & OPT_V)? prn : NULL;

    return ordered_estimate(list, dset, PROBIT, opt, vprn);
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

MODEL multinomial_logit (int *list, DATASET *dset, 
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
    double x, ess, sigma;
    int t;

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, dset->Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, dset->Z[dv]);

    /* make the VCV matrix before messing with the model stats */
    makevcv(pmod, pmod->sigma);

    ess = 0.0;

    for (t=0; t<dset->n; t++) {
	x = pmod->yhat[t];
	if (!na(x)) {
	    pmod->yhat[t] = lmax / (1.0 + exp(-x));
	    pmod->uhat[t] = dset->Z[dv][t] - pmod->yhat[t];
	    ess += pmod->uhat[t] * pmod->uhat[t];
	}
    }

    sigma = sqrt(ess / pmod->dfd);

    pmod->list[1] = dv;
    gretl_model_set_double(pmod, "lmax", lmax);
    gretl_model_set_double(pmod, "ess_orig", ess);
    gretl_model_set_double(pmod, "sigma_orig", sigma);
    pmod->ci = LOGISTIC;
    ls_criteria(pmod);

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
 * @dset: dataset struct.
 *
 * Estimate the model given in @list using the logistic transformation
 * of the dependent variable.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL logistic_model (const int *list, double lmax,
		      DATASET *dset)
{
    int *llist = NULL;
    int dv = list[1];
    double real_lmax;
    MODEL lmod;
    int err = 0;

    fprintf(stderr, "logistic model: lmax = %g\n", lmax);

    gretl_model_init(&lmod, dset); 

    llist = gretl_list_copy(list);
    if (llist == NULL) {
	err = E_ALLOC;
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

    lmod = lsq(llist, dset, OLS, OPT_A);
    if (!lmod.errcode) {
	rewrite_logistic_stats(dset, &lmod, dv, real_lmax);
	set_model_id(&lmod);
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

    if (xna(c)) {
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
	if (xna(P)) {
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
	pprintf(prn, "  Left:   P-value = %g\n", PL);
	pprintf(prn, "  Right:  P-value = %g\n", PR);
	pprintf(prn, "  2-Tail: P-value = %g\n", P2);
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

MODEL interval_model (int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    MODEL intmod;
    void *handle;
    MODEL (* interval_estimate) (int *, DATASET *, gretlopt, 
				 PRN *);

    gretl_error_clear();

    interval_estimate = get_plugin_function("interval_estimate", &handle);

    if (interval_estimate == NULL) {
	gretl_model_init(&intmod, dset);
	intmod.errcode = E_FOPEN;
	return intmod;
    }

    intmod = (*interval_estimate) (list, dset, opt, prn);

    close_plugin(handle);

    set_model_id(&intmod);

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
    void *handle;
    MODEL tmod;    

    gretl_error_clear();

    tobit_estimate = get_plugin_function("tobit_via_intreg", &handle);
    if (tobit_estimate == NULL) {
	gretl_model_init(&tmod, dset);
	tmod.errcode = E_FOPEN;
	return tmod;
    }

    tmod = (*tobit_estimate) (list, llim, rlim, dset, opt, prn);

    close_plugin(handle);
    set_model_id(&tmod);

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
    void *handle;
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

    duration_estimate = get_plugin_function("duration_estimate", 
					    &handle);

    if (duration_estimate == NULL) {
	dmod.errcode = E_FOPEN;
	return dmod;
    }

    (*duration_estimate) (&dmod, censvar, dset, opt, prn);

    close_plugin(handle);

    set_model_id(&dmod);

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
    void *handle;
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

    count_data_estimate = get_plugin_function("count_data_estimate", 
					      &handle);

    if (count_data_estimate == NULL) {
	cmod.errcode = E_FOPEN;
	return cmod;
    }

    (*count_data_estimate) (&cmod, ci, offvar, dset, opt, prn);

    close_plugin(handle);

    set_model_id(&cmod);

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
    void *handle;
    MODEL (* heckit_estimate) (const int *, DATASET *, 
			       gretlopt, PRN *);

    gretl_error_clear();

    heckit_estimate = get_plugin_function("heckit_estimate", &handle);
    if (heckit_estimate == NULL) {
	gretl_model_init(&model, dset);
	model.errcode = E_FOPEN;
	return model;
    }

    model = (*heckit_estimate) (list, dset, opt, prn);

    close_plugin(handle);

    set_model_id(&model);

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
    void *handle;
    MODEL (* reprobit_estimate) (const int *, DATASET *,
				 gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();

    if (!dataset_is_panel(dset)) {
	err = E_PDWRONG;
    } else {	
	reprobit_estimate = get_plugin_function("reprobit_estimate", &handle);
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

    close_plugin(handle);

    set_model_id(&model);

    return model;
}
