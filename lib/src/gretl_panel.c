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

/* panel estimation and dataset handling for gretl */

#include "libgretl.h"
#include "gretl_f2c.h"
#include "clapack_double.h"
#include "gretl_model.h"
#include "gretl_panel.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_string_table.h"
#include "matrix_extra.h" /* for testing */

/**
 * SECTION:gretl_panel
 * @short_description: estimation of panel data models
 * @title: Panel data
 * @include: libgretl.h
 *
 * Provides support for estimating panel data models
 * such as fixed and random effects.
 */

#define PDEBUG 0
#define CDEBUG 0

enum vcv_ops {
    VCV_INIT,
    VCV_SUBTRACT
};

typedef struct unit_ts_ unit_ts;

struct unit_ts_ {
    int *t0;  /* start of longest contiguous time series, per unit */
    int *T;   /* length of such sequence, per unit */
};

typedef struct panelmod_t_ panelmod_t;

struct panelmod_t_ {
    gretlopt opt;         /* option flags */
    int nunits;           /* total cross-sectional units in sample range */
    int effn;             /* effective (included) cross-section units */
    int T;                /* times-series length of panel */
    int Tmax;             /* maximum number of usable obs per unit */
    int Tmin;             /* minimum number of usable obs per unit > 0 */
    double Tbar;          /* harmonic mean of per-unit time-series lengths */
    int NT;               /* total observations used (based on pooled model) */
    int ntdum;            /* number of time dummies added */
    int *unit_obs;        /* array of number of observations per x-sect unit */
    char *varying;        /* array to record properties of pooled-model regressors */
    int *vlist;           /* list of time-varying variables from pooled model */
    int balanced;         /* 1 if the model dataset is balanced, else 0 */
    int dfH;              /* number of coeffs for Hausman test */
    int Fdfn;             /* numerator df, F for differing intercepts */
    double Fdfd;          /* denominator df, F for differing intercepts */
    double theta;         /* quasi-demeaning coefficient */
    double theta_bar;     /* mean of theta_i (unbalanced case) */
    double Ffe;           /* joint significance of fixed effects */
    double BP;            /* Breusch-Pagan test statistic */
    double H;             /* Hausman test statistic */
    gretl_matrix *bdiff;  /* array of coefficient differences */
    gretl_matrix *Sigma;  /* Hausman covariance matrix */
    double s2b;           /* residual variance, group means regression */
    double s2e;           /* residual variance, fixed-effects regression */
    double s2v;           /* estimate of between variance */
    double ubPub;         /* for use in unbalanced Swamy-Arora */
    double dw;            /* Durbin-Watson, for transfer to random effects */
    double rho;           /* for transfer to random effects */
    int *small2big;       /* data indexation array */
    int *big2small;       /* reverse data indexation array */
    MODEL *pooled;        /* reference model (pooled OLS) */
    MODEL *realmod;       /* fixed or random effects model */
    double *re_uhat;      /* "fixed" random-effects residuals */
};

struct {
    int n;      /* number of cross-sectional units */
    int T;      /* number of observations per unit */
    int offset; /* sampling offset into full dataset */
} panidx;

/* for testing purposes */
int stata_sa;
int IGLS;
char glsmat[MAXLEN];

static int varying_vars_list (const DATASET *dset, panelmod_t *pan);
static void calculate_Tbar (panelmod_t *pan);

/* translate from (i = unit, t = time period for that unit) to
   overall 0-based index into the data set */
#define panel_index(i,t) (i * panidx.T + t + panidx.offset)

#define panel_missing(p, t) (na(p->pooled->uhat[t]))


static void
panel_index_init (const DATASET *dset, int nunits, int T)
{
    panidx.n = nunits;
    panidx.T = T;
    panidx.offset = dset->t1;

#if PDEBUG
    fprintf(stderr, "panel_index_init: n=%d, T=%d, offset=%d\n",
            panidx.n, panidx.T, panidx.offset);
#endif
}

/* Allocate the indexation arrays that allow us to translate between
   the full-length data array and the possibly shorter array used for
   transformed data (de-meaned or quasi-demeaned).
*/

static int allocate_data_finders (panelmod_t *pan, int bign)
{
    int s;

    if (pan->small2big != NULL) {
        /* already done */
        return 0;
    }

    pan->small2big = malloc(pan->NT * sizeof *pan->small2big);
    pan->big2small = malloc(bign * sizeof *pan->big2small);

    if (pan->small2big == NULL || pan->big2small == NULL) {
        return E_ALLOC;
    }

    for (s=0; s<bign; s++) {
        pan->big2small[s] = -1;
    }

    return 0;
}

#define small_index(p,t) ((p->big2small == NULL)? t : p->big2small[t])
#define big_index(p,t)   ((p->small2big == NULL)? t : p->small2big[t])

static void panelmod_init (panelmod_t *pan)
{
    pan->nunits = 0;
    pan->effn = 0;
    pan->T = 0;
    pan->Tmax = 0;
    pan->Tmin = 0;
    pan->Tbar = 0;
    pan->NT = 0;
    pan->ntdum = 0;
    pan->unit_obs = NULL;
    pan->varying = NULL;
    pan->vlist = NULL;
    pan->opt = OPT_NONE;

    pan->balanced = 1;
    pan->theta = NADBL;
    pan->theta_bar = NADBL;
    pan->ubPub = NADBL;
    pan->dw = NADBL;
    pan->rho = NADBL;

    pan->Ffe = NADBL;
    pan->Fdfn = 0;
    pan->Fdfd = 0;

    pan->BP = NADBL;
    pan->H = NADBL;
    pan->dfH = 0;

    pan->bdiff = NULL;
    pan->Sigma = NULL;

    pan->small2big = NULL;
    pan->big2small = NULL;

    pan->pooled = NULL;
    pan->realmod = NULL;
    pan->re_uhat = NULL;
}

static void panelmod_free (panelmod_t *pan)
{
    free(pan->unit_obs);
    free(pan->varying);
    free(pan->vlist);

    gretl_matrix_free(pan->bdiff);
    gretl_matrix_free(pan->Sigma);

    free(pan->small2big);
    free(pan->big2small);
    free(pan->re_uhat);

    free(pan->realmod);
}

/* test variable number v against the (possibly reduced) regression
   list which contains only time-varying regressors
*/

static int var_is_varying (const panelmod_t *pan, int v)
{
    if (v != 0) {
        int i;

        for (i=2; i<=pan->vlist[0]; i++) {
            if (pan->vlist[i] == v) {
                return 1;
            }
        }
    }

    return 0;
}

/* retrieve X'X^{-1} from the pooled or fixed effects model */

static gretl_matrix *panel_model_xpxinv (MODEL *pmod, int *err)
{
    gretl_matrix *X;
    int k = pmod->ncoeff;
    double x;
    int i, j, m;

#if 0
    fprintf(stderr, "panel_model_xpxinv...\n");
#endif

    if (gretl_model_get_int(pmod, "vcv_xpx")) {
        /* already done */
        gretl_model_set_int(pmod, "vcv_xpx", 0);
    } else if (pmod->vcv != NULL) {
        double s2 = pmod->sigma * pmod->sigma;
        int nv = (k * k + k) / 2;

        for (i=0; i<nv; i++) {
            pmod->vcv[i] /= s2;
        }
    } else {
        *err = makevcv(pmod, 1.0);
        if (*err) {
            return NULL;
        }
        gretl_model_set_int(pmod, "vcv_xpx", 1);
    }

    X = gretl_matrix_alloc(k, k);
    if (X == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    m = 0;
    for (j=0; j<k; j++) {
        for (i=j; i<k; i++) {
            x = pmod->vcv[m++];
            gretl_matrix_set(X, i, j, x);
            gretl_matrix_set(X, j, i, x);
        }
    }

#if 0
    fprintf(stderr, "pmod->sigma = %g\n", pmod->sigma);
    printlist(pmod->list, "pmod->list");
    gretl_matrix_print(X, "panel (X'X)^{-1}");
#endif

    return X;
}

/* Beck and Katz, as outlined in Greene.  We offer the following
   for pooled OLS and FE.  Greene writes (in effect):

   Var(b) = A^{-1} W A^{-1}

   where A = \sum_{i=1}^n X'_i X_i (called "XX" below)
         W = \sum_{i=1}^n \sum_{j=1}^n \sigma_{ij} X'_{i} X_{j}

   and \sigma_{ij} is estimated as (1/T) \sum_{t=1}^T e_i e_j,
   with the e's being OLS (or FE) residuals.

   In computing this, we need to check that W is positive definite:
   this seems not to be guaranteed for unbalanced panels.
*/

static int
beck_katz_vcv (MODEL *pmod, panelmod_t *pan, const DATASET *dset,
               gretl_matrix *XX, gretl_matrix *W, gretl_matrix *V)
{
    gretl_matrix *Xi = NULL;
    gretl_matrix *Xj = NULL;
    int T = pan->T;
    int k = pmod->ncoeff;
    int s, si, sj;
    int i, j, p, v, t;
    int err = 0;

    Xi = gretl_matrix_alloc(pan->Tmax, k);
    Xj = gretl_matrix_alloc(pan->Tmax, k);

    if (Xi == NULL || Xj == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];
        double sii = 0.0;

        if (Ti == 0) {
            continue;
        }

        Xi = gretl_matrix_reuse(Xi, Ti, k);

        s = 0;
        for (t=0; t<T; t++) {
            si = panel_index(i, t);
            if (!na(pmod->uhat[si])) {
                sii += pmod->uhat[si] * pmod->uhat[si];
                si = small_index(pan, si);
                for (p=0; p<k; p++) {
                    v = pmod->list[p + 2];
                    gretl_matrix_set(Xi, s, p, dset->Z[v][si]);
                }
                s++;
            }
        }

        sii /= Ti;
        sii = sqrt(sii);

        /* "diagonal" component of W */
        gretl_matrix_multiply_by_scalar(Xi, sii);
        gretl_matrix_multiply_mod(Xi, GRETL_MOD_TRANSPOSE,
                                  Xi, GRETL_MOD_NONE,
                                  W, GRETL_MOD_CUMULATE);

        for (j=i+1; j<pan->nunits; j++) {
            int Tij = 0;
            double sij = 0.0;

            if (pan->unit_obs[j] == 0) {
                continue;
            }

            /* count matching observations */
            for (t=0; t<T; t++) {
                si = panel_index(i, t);
                sj = panel_index(j, t);
                if (!na(pmod->uhat[si]) && !na(pmod->uhat[sj])) {
                    sij += pmod->uhat[si] * pmod->uhat[sj];
                    Tij++;
                }
            }

            if (Tij == 0) {
                continue;
            }

            sij /= Tij;
            Xi = gretl_matrix_reuse(Xi, Tij, k);
            Xj = gretl_matrix_reuse(Xj, Tij, k);

            s = 0;
            for (t=0; t<T; t++) {
                si = panel_index(i, t);
                sj = panel_index(j, t);
                if (!na(pmod->uhat[si]) && !na(pmod->uhat[sj])) {
                    si = small_index(pan, si);
                    sj = small_index(pan, sj);
                    for (p=0; p<k; p++) {
                        v = pmod->list[p + 2];
                        gretl_matrix_set(Xi, s, p, dset->Z[v][si]);
                        gretl_matrix_set(Xj, s, p, dset->Z[v][sj]);
                    }
                    s++;
                }
            }

            /* cumulate s_ij * Xi'Xj into W (using V as storage) */
            gretl_matrix_multiply_by_scalar(Xi, sij);
            gretl_matrix_multiply_mod(Xi, GRETL_MOD_TRANSPOSE,
                                      Xj, GRETL_MOD_NONE,
                                      V, GRETL_MOD_NONE);
            gretl_matrix_add_to(W, V);
            /* and take in s_ij * Xj'Xi */
            gretl_matrix_add_transpose_to(W, V);
        }
    }

    /* check that the middle term is p.d. */
    gretl_matrix_copy_values(V, W);
    err = gretl_matrix_cholesky_decomp(V);
    if (err) {
        fprintf(stderr, "beck_katz_vcv: matrix W is not p.d.\n");
        gretl_model_set_int(pmod, "panel_bk_failed", 1);
        goto bailout;
    }

    /* form V = (Xi'Xi)^{-1} W (Xi'Xi)^{-1} */
    gretl_matrix_qform(XX, GRETL_MOD_NONE, W,
                       V, GRETL_MOD_NONE);

    gretl_model_set_vcv_info(pmod, VCV_PANEL, PANEL_BK);

 bailout:

    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);

    return err;
}

/* helper function for panel VCV clustered by period */

static int *get_panel_tobs (panelmod_t *pan)
{
    int *tobs = calloc(pan->T, sizeof *tobs);

    if (tobs != NULL) {
        int i, t;

        for (i=0; i<pan->nunits; i++) {
            for (t=0; t<pan->T; t++) {
                if (!panel_missing(pan, panel_index(i, t))) {
                    tobs[t] += 1;
                }
            }
        }
    }

    return tobs;
}

static void augment_clustered_W (const gretl_matrix *e,
                                 const gretl_matrix *X,
                                 gretl_matrix *eX,
                                 gretl_matrix *W)
{
    gretl_matrix_multiply_mod(e, GRETL_MOD_TRANSPOSE,
                              X, GRETL_MOD_NONE,
                              eX, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(eX, GRETL_MOD_TRANSPOSE,
                              eX, GRETL_MOD_NONE,
                              W, GRETL_MOD_CUMULATE);
}

static void finalize_clustered_vcv (MODEL *pmod,
                                    panelmod_t *pan,
                                    const gretl_matrix *XX,
                                    const gretl_matrix *W,
                                    gretl_matrix *V,
                                    int nc)
{
    double adj;

    /* form V(b) = (X'X)^{-1} W (X'X)^{-1} */
    gretl_matrix_qform(XX, GRETL_MOD_NONE, W,
                       V, GRETL_MOD_NONE);

    /* Now, what (if anything) are we going to do in terms of
       degrees-of-freedom adjustment?
    */

    if (pmod->opt & OPT_N) {
        /* --no-df-corr */
        return;
    }

    if (pmod->ci == IVREG) {
        /* compatible with gretl 2023b: do we want this? */
        return;
    }

    /* initial adjustment using the number of clusters */
    adj = nc / (nc - 1.0);

    if (pmod->ci == IVREG || pmod->opt & OPT_N) {
        /* Just apply the @adj calculated above? This is said to
           be the "asymtotic-like" approach in Stata; see
           https://www.stata.com/meeting/13uk/nichols_crse.pdf,
           under "Finite-Sample Adjustment". Right now this
           is never reached!
        */
        ;
    } else {
        /* Apply the full Cameron-Gelbach-Miller adjustment */
        adj *= (pmod->nobs - 1.0) / (pmod->nobs - pmod->ncoeff);
    }

    gretl_matrix_multiply_by_scalar(V, adj);
}

#define CI_TIME (-2)
#define CI_UNIT (-3)

typedef struct cluster_info_ {
    const double *cz;    /* series of all cluster-var values */
    gretl_matrix *cvals; /* sorted vector of unique cluster-var values */
    int dcid[2];         /* IDs of up to two @dset series */
    int pcid[3];         /* IDs of up to three @pset series */
    int nc[2];           /* number(s) of clusters */
    gint8 target;        /* 0 or 1 for first or second cluster var */
    gint8 pooled;        /* flag for the pooled OLS case */
} cluster_info;

#define by_time_and_unit(ci) (ci->dcid[0] < -1 && ci->dcid[1] < -1)
#define by_time_or_unit(ci) (ci->dcid[0] < -1 || ci->dcid[1] < -1)
#define by_time(ci) (ci->dcid[0] == CI_TIME)
#define by_unit(ci) (ci->dcid[0] == CI_UNIT)
#define generic(ci) (ci->dcid[0] > 0)
#define two_way(ci) (ci->dcid[0] != -1 && ci->dcid[1] != -1)
#define is_present(v) (v != -1)

static cluster_info *cluster_info_new (panelmod_t *pan)
{
    cluster_info *ci = calloc(1, sizeof *ci);

    ci->dcid[0] = ci->dcid[1] = -1;
    ci->target = 0;
    ci->pooled = (pan->opt & OPT_P)? 1 : 0;

    return ci;
}

static void cluster_info_free (cluster_info *ci)
{
    if (ci != NULL) {
        gretl_matrix_free(ci->cvals);
        free(ci);
    }
}

static const char *cluster_name (int v, const DATASET *dset)
{
    if (v > 0) {
	return dset->varname[v];
    } else if (v == CI_TIME) {
	return "time";
    } else {
	return "unit";
    }
}

static int
time_cluster_vcv (MODEL *pmod, panelmod_t *pan, const DATASET *dset,
                  const gretl_matrix *XX, gretl_matrix *W,
                  gretl_matrix *V, cluster_info *ci)
{
    gretl_vector *et = NULL;
    gretl_matrix *Xt = NULL;
    gretl_vector *eXt = NULL;
    int *tobs = NULL;
    int N = pan->effn;
    int k = pmod->ncoeff;
    int idx = ci->target;
    int i, j, v, s, t;
    int n_c = 0;
    int err = 0;

    tobs = get_panel_tobs(pan);
    if (tobs == NULL) {
        return E_ALLOC;
    }

    et  = gretl_column_vector_alloc(N);
    Xt  = gretl_matrix_alloc(N, k);
    eXt = gretl_vector_alloc(k);

    if (et == NULL || Xt == NULL || eXt == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    for (t=0; t<pan->T; t++) {
        int Nt = tobs[t]; /* # of units represented in period @t */
        int p = 0;

        if (Nt == 0) {
            continue;
        }

        n_c++;
        et = gretl_matrix_reuse(et, Nt, 1);
        Xt = gretl_matrix_reuse(Xt, Nt, k);

        for (i=0; i<pan->nunits; i++) {
            s = panel_index(i, t);
            if (na(pmod->uhat[s])) {
                continue;
            }
            gretl_vector_set(et, p, pmod->uhat[s]);
            s = small_index(pan, s);
            for (j=0; j<k; j++) {
                v = pmod->list[j+2];
                if (s < 0) {
                    gretl_matrix_set(Xt, p, j, 0);
                } else {
                    gretl_matrix_set(Xt, p, j, dset->Z[v][s]);
                }
            }
            p++;
            if (p == Nt) {
                /* we're done with this period */
                break;
            }
        }
        augment_clustered_W(et, Xt, eXt, W);
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, pan->Tmax);
    if (!two_way(ci)) {
	/* finish the job */
	gretl_model_set_vcv_info(pmod, VCV_PANEL, PANEL_TIME);
	gretl_model_set_int(pmod, "n_clusters", n_c);
    } else {
	ci->nc[idx] = n_c;
    }

 bailout:

    gretl_matrix_free(et);
    gretl_matrix_free(Xt);
    gretl_matrix_free(eXt);
    free(tobs);

    return err;
}

static int panel_cval_count (MODEL *pmod, panelmod_t *pan,
                             cluster_info *ci, int i,
			     const double *cvar)
{
    double cvi = ci->cvals->val[i];
    int s, t, cc = 0;

    if (cvar != NULL) {
	/* the pooled OLS case */
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t]) && cvar[t] == cvi) {
		cc++;
	    }
	}
    } else {
	/* other cases */
	for (s=0; s<pmod->nobs; s++) {
	    t = big_index(pan, s);
	    if (!na(pmod->uhat[t]) && ci->cz[s] == cvi) {
		cc++;
	    }
	}
    }

    return cc;
}

static int panel_cval_count_max (MODEL *pmod,
                                 panelmod_t *pan,
                                 cluster_info *ci,
				 const double *cvar)
{
    int n = gretl_vector_get_length(ci->cvals);
    int i, cc, cmax = 0;

    for (i=0; i<n; i++) {
        cc = panel_cval_count(pmod, pan, ci, i, cvar);
        if (cc > cmax) {
            cmax = cc;
        }
    }

    return cmax;
}

static int
generic_cluster_vcv (MODEL *pmod, panelmod_t *pan, const DATASET *dset,
                     const gretl_matrix *XX, gretl_matrix *W,
                     gretl_matrix *V, cluster_info *ci)
{
    const double *cvar = NULL;
    gretl_vector *ei = NULL;
    gretl_matrix *Xi = NULL;
    gretl_vector *eXi = NULL;
    int R, M = ci->cvals->rows;
    int k = pmod->ncoeff;
    int idx = ci->target;
    int i, j, v, t;
    int n_c = 0;
    int err = 0;

    if (ci->pooled) {
	cvar = dset->Z[ci->dcid[idx]];
    }

    R = panel_cval_count_max(pmod, pan, ci, cvar);
#if CDEBUG
    fprintf(stderr, "generic_cluster_vcv: M=%d, R=%d\n", M, R);
#endif
    ei  = gretl_column_vector_alloc(R);
    Xi  = gretl_matrix_alloc(R, k);
    eXi = gretl_vector_alloc(k);

    if (ei == NULL || Xi == NULL || eXi == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    for (i=0; i<M; i++) {
        double cvi = ci->cvals->val[i];
        int Ni = panel_cval_count(pmod, pan, ci, i, cvar);
        int s, p = 0;

#if CDEBUG > 1
        fprintf(stderr, "cvals[%d]=%g, count=%d\n", i, cvi, Ni);
#endif
        if (Ni == 0) {
            continue;
        }

        ei = gretl_matrix_reuse(ei, Ni, -1);
        Xi = gretl_matrix_reuse(Xi, Ni, -1);

	if (ci->pooled) {
	    /* working with the full dataset */
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!na(pmod->uhat[t]) && cvar[t] == cvi) {
		    gretl_vector_set(ei, p, pmod->uhat[t]);
		    for (j=0; j<k; j++) {
			v = pmod->list[j+2];
			gretl_matrix_set(Xi, p, j, dset->Z[v][t]);
		    }
		    p++;
		}
		if (p == Ni) {
		    break;
		}
	    }
	} else {
	    /* working with dataset from which NAs have been purged */
	    for (s=0; s<pmod->nobs; s++) {
		t = big_index(pan, s);
		if (!na(pmod->uhat[t]) && ci->cz[s] == cvi) {
		    gretl_vector_set(ei, p, pmod->uhat[t]);
		    for (j=0; j<k; j++) {
			v = pmod->list[j+2];
			gretl_matrix_set(Xi, p, j, dset->Z[v][s]);
		    }
		    p++;
		}
		if (p == Ni) {
		    break;
		}
	    }
	}
        augment_clustered_W(ei, Xi, eXi, W);
        n_c++;
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, n_c);

    if (!two_way(ci)) {
        /* not doing two-way clustering, so finish the job */
	const char *cname = dset->varname[ci->pcid[idx]];

        gretl_model_set_cluster_vcv_info(pmod, cname, NULL);
        gretl_model_set_int(pmod, "n_clusters", n_c);
    } else if (idx < 2) {
	ci->nc[idx] = n_c;
    }

 bailout:

    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

static int
panel_white_vcv (MODEL *pmod, panelmod_t *pan, const DATASET *dset,
		 const gretl_matrix *XX, gretl_matrix *W,
		 gretl_matrix *V, cluster_info *ci)
{
    gretl_vector *ei = NULL;
    gretl_matrix *Xi = NULL;
    gretl_vector *eXi = NULL;
    int k = pmod->ncoeff;
    int i, v, s, t;
    int n_c = 0;
    int err = 0;

    ei  = gretl_column_vector_alloc(1);
    Xi  = gretl_matrix_alloc(1, k);
    eXi = gretl_vector_alloc(k);

    if (ei == NULL || Xi == NULL || eXi == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    for (s=0; s<pmod->nobs; s++) {
	t = big_index(pan, s);
	if (!na(pmod->uhat[t])) {
	    ei->val[0] = pmod->uhat[t];
	    for (i=0; i<k; i++) {
		v = pmod->list[i+2];
		gretl_vector_set(Xi, i, dset->Z[v][s]);
	    }
	    augment_clustered_W(ei, Xi, eXi, W);
	    n_c++;
	}
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, n_c);
    if (!two_way(ci)) {
	/* should we offer "plain White" (HC0) as an option? */
	gretl_model_set_vcv_info(pmod, VCV_HC, 0);
    }

 bailout:

    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

/* HAC covariance matrix for pooled, fixed- or random-effects models,
   given "fixed T and large N".  In the case of "large T" a different
   form is needed for robustness in respect of autocorrelation.  See
   Arellano, "Panel Data Econometrics" (Oxford, 2003), pages 18-19.
*/

static int
unit_cluster_vcv (MODEL *pmod, panelmod_t *pan, const DATASET *dset,
		  const gretl_matrix *XX, gretl_matrix *W,
		  gretl_matrix *V, cluster_info *ci)
{
    gretl_vector *ei = NULL;
    gretl_matrix *Xi = NULL;
    gretl_vector *eXi = NULL;
    int T = pan->Tmax;
    int k = pmod->ncoeff;
    int idx = ci->target;
    int i, j, v, s, t;
    int err = 0;

    ei  = gretl_column_vector_alloc(T);
    Xi  = gretl_matrix_alloc(T, k);
    eXi = gretl_vector_alloc(k);

    if (ei == NULL || Xi == NULL || eXi == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];
        int p = 0;

        if (Ti == 0) {
            continue;
        }

        ei = gretl_matrix_reuse(ei, Ti, 1);
        Xi = gretl_matrix_reuse(Xi, Ti, k);

        for (t=0; t<pan->T; t++) {
            s = panel_index(i, t);
            if (na(pmod->uhat[s])) {
                continue;
            }
            gretl_vector_set(ei, p, pmod->uhat[s]);
            s = small_index(pan, s);
            for (j=0; j<k; j++) {
                v = pmod->list[j+2];
                if (s < 0) {
                    gretl_matrix_set(Xi, p, j, 0);
                } else {
                    gretl_matrix_set(Xi, p, j, dset->Z[v][s]);
                }
            }
            if (++p == Ti) {
                /* we're done with this unit */
                break;
            }
        }
        augment_clustered_W(ei, Xi, eXi, W);
    }

    finalize_clustered_vcv(pmod, pan, XX, W, V, pan->effn);

    if (!two_way(ci)) {
        /* finish the job if not doing two-way clustering */
	gretl_model_set_vcv_info(pmod, VCV_PANEL, PANEL_HAC);
	gretl_model_set_int(pmod, "n_clusters", pan->effn);
    } else {
	ci->nc[idx] = pan->effn;
    }

 bailout:

    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

/* In response to "$time" or "$unit" appearing in the context of the
   --cluster option for a basic panel-data estimator, we construct an
   on-the-fly series that represents the time period -- or if @srcv
   equals CI_UNIT, the cross-sectional unit. We need this only in
   case of two-way clustering with one of $time or $unit combined
   with a named series.
*/

static void make_unit_or_period (MODEL *pmod,
				 panelmod_t *pan,
				 double *y,
				 int srcv)
{
    int do_unit = (srcv == CI_UNIT);
    int i, t, s;

    for (i=0; i<pan->nunits; i++) {
        if (pan->unit_obs[i] == 0) {
	    continue;
	}
	for (t=0; t<pan->T; t++) {
	    s = panel_index(i, t);
	    if (na(pmod->uhat[s])) {
		continue;
	    }
	    s = small_index(pan, s);
	    y[s] = do_unit ? i+1 : t+1;
	}
    }
}

static int cluster_missval_error (DATASET *dset, int v, int t)
{
    char obs[OBSLEN];

    gretl_errmsg_sprintf("cluster variable '%s' missing at observation %s",
			 dset->varname[v], ntolabel(obs, t, dset));
    return E_MISSDATA;
}

/* The role of this function is to transcribe a clustering variable
   into the panel-estimation dataset, @pset, the destination series ID
   being given by ci->pcid. As for the source of the data, there are
   two cases:

   (1) a named series present in the main dataset, @dset; or

   (2) an automatically generated series, responding to the keyword
   "$time" or "$unit".

   The second case arises only in the case of "mixed" two-way
   clustering (where one of the variables comes from @dset and the
   other is automatic). If the cluster specification contains only
   automatic variables there's no need for transcription.

   On entry to this function ci->target should be set to either 0,
   indicating that we're working on the first (or only) cluster
   variable, or 1 to indicate that we're working on the second.
*/

static int transcribe_cluster_var (MODEL *pmod,
                                   panelmod_t *pan,
                                   DATASET *pset,
                                   DATASET *dset,
                                   cluster_info *ci)
{
    const double *src = NULL;
    double *dest = NULL;
    int idx = ci->target;
    int srcv, t, s;
    int err = 0;

    /* source and destination of transcription */
    srcv = ci->dcid[idx];
    dest = pset->Z[ci->pcid[idx]];

#if CDEBUG
    fprintf(stderr, "  transcribe: target=%d, pcid=%d, srcv=%d (%s)\n",
	    idx, ci->pcid[idx], srcv, srcv > 0 ? dset->varname[srcv] :
	    srcv == CI_TIME ? "TIME" : "UNIT");
#endif

    if (srcv == CI_TIME || srcv == CI_UNIT) {
	/* the "automatic" case */
	make_unit_or_period(pmod, pan, dest, srcv);
	return 0;
    }

    src = dset->Z[srcv];

    if (ci->pooled) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t]) && na(src[t])) {
		err = cluster_missval_error(dset, srcv, t);
		break;
            }
        }
    } else if (pset->n == dset->n) {
        /* balanced */
        for (t=0; t<dset->n; t++) {
            if (na(src[t])) {
                err = cluster_missval_error(dset, srcv, t);
                break;
            }
        }
        if (!err) {
            memcpy(dest, src, pset->n * sizeof **pset->Z);
        }
    } else {
	/* unbalanced */
        for (s=0; s<pmod->nobs && !err; s++) {
            t = big_index(pan, s);
            if (!na(pmod->uhat[t])) {
                if (na(src[t])) {
                    err = cluster_missval_error(dset, srcv, t);
                } else {
                    dest[s] = src[t];
                }
            }
        }
    }

    if (!err) {
	int nread = pmod->nobs;

        if (ci->cvals != NULL) {
            gretl_matrix_free(ci->cvals);
        }
	if (ci->pooled) {
	    dest += pmod->t1;
	    nread = pmod->t2 - pmod->t1 + 1;
	}
        ci->cvals = gretl_matrix_values(dest, nread, OPT_S, &err);
        ci->cz = dest;
#if CDEBUG > 1
        gretl_matrix_print(ci->cvals, "cluster var values");
#endif
    }

#if CDEBUG
    fprintf(stderr, "transcribe: returning %d\n", err);
#endif

    return err;
}

/* For a clustering variable named by @s, find out if it's a named
   series, or a keyword ("$time", "$unit"), or neither (an error). To
   signal the presence of one of the keywords we assign one of CI_TIME
   or CI_UNIT in place of a regular series ID number.
*/

static int get_id_or_special (const char *s,
			      const DATASET *dset,
			      cluster_info *ci)
{
    int *pid = &ci->dcid[ci->target];
    int id = current_series_index(dset, s);
    int err = 0;

    if (id > 0) {
	*pid = id;
    } else if (id == 0) {
	/* "const" will not do! */
	err = E_DATA;
    } else if (!strcmp(s, "$time")) {
	*pid = CI_TIME;
    } else if (!strcmp(s, "$unit")) {
	*pid = CI_UNIT;
    } else {
	err = E_UNKVAR;
    }

    return err;
}

static int process_cluster_string (const char *s,
				   const DATASET *dset,
				   cluster_info *ci)
{
    const char *p = strchr(s, ',');
    int err = 0;

#if CDEBUG
    fprintf(stderr, "process_cluster_string: '%s'\n", s);
#endif

    if (p == NULL) {
	/* no comma: should have a single name */
	ci->target = 0;
	err = get_id_or_special(s, dset, ci);
    } else {
	/* should have two names */
	gchar *s1 = NULL;
	gchar *s2 = NULL;

	s1 = g_strndup(s, p - s);
	ci->target = 0;
	err = get_id_or_special(s1, dset, ci);
	if (!err) {
	    s2 = g_strdup(p + 1);
	    ci->target = 1;
	    err = get_id_or_special(s2, dset, ci);
	    ci->target = 0;
	}
	g_free(s1);
	g_free(s2);
    }

    return err;
}

/* In this function we're looking for one or two identifiers supplied
   via the string parameter to the --cluster option.  We accept
   "$unit" and "$time" as automatic references to the cross-sectional
   or time dimension of the panel; otherwise names of series in @dset
   are needed.

   If this succeeds, we transcribe the first (or only) clustering
   series into the panel dataset, @pset.  If there's a second argument
   we record its ID as dcid[1] under the @ci struct; it will get
   transcribed later, in two_way_cluster_vcv().
*/

static int check_cluster_var (MODEL *pmod,
                              panelmod_t *pan,
                              DATASET *pset,
                              DATASET *dset,
                              cluster_info *ci)
{
    /* note: OPT_L means that we came here from an invocation of
       "ols", so the relevant command index for retrieval of an
       option parameter is OLS and not PANEL
    */
    int opt_ci = (pan->opt & OPT_L)? OLS : PANEL;
    const char *s = get_optval_string(opt_ci, OPT_C);
    int err = 0;

#if CDEBUG
    fprintf(stderr, "check_cluster_var: s = '%s'\n", s);
#endif

    if (s == NULL || *s == '\0') {
        err = E_DATA;
    } else {
	err = process_cluster_string(s, dset, ci);
    }

#if CDEBUG
    /* cluster series IDs, relative to full dataset */
    fprintf(stderr, "  dcid[0] %d, dcid[1] %d\n", ci->dcid[0], ci->dcid[1]);
#endif

    if (err) {
	return err;
    } else if (by_time_and_unit(ci)) {
	; /* nothing more to be done */
    } else if (by_time_or_unit(ci) && !two_way(ci)) {
	; /* nothing more to be done */
    } else {
	int adds[3] = {0};
        int v, n_add = 0;

	if (ci->pooled && ci->dcid[0] > 0) {
	    /* using an existing series for first var */
	    ci->pcid[0] = ci->dcid[0];
	} else {
	    adds[0] = 1;
	}
	if (is_present(ci->dcid[1])) {
	    if (pset == dset && ci->dcid[1] > 0) {
		/* we'll just need one extra series for the combination */
		ci->pcid[1] = ci->dcid[1];
		adds[2] = 1;
	    } else {
		/* we'll need more series plus the combination */
		adds[1] = adds[2] = 1;
	    }
	}

	n_add = adds[0] + adds[1] + adds[2];

	if (n_add > 0) {
	    err = dataset_add_series(pset, n_add);
	    if (err) {
		return err;
	    }
	    v = pset->v - n_add;
	    if (adds[0]) {
		ci->pcid[0] = v;
		strcpy(pset->varname[v++], cluster_name(ci->dcid[0], dset));
	    }
	    if (adds[1]) {
		ci->pcid[1] = v;
		strcpy(pset->varname[v++], cluster_name(ci->dcid[1], dset));
	    }
	    if (adds[2]) {
		ci->pcid[2] = v;
		strcpy(pset->varname[v], "combo");
	    }
	}

	/* transcribe (if needed) the first cluster var */
	err = transcribe_cluster_var(pmod, pan, pset, dset, ci);
    }

    return err;
}

/* Fix-up, à la Cameron, Gelbach and Miller, for non-positive
   semi-definite covariance matrix in case of non-nested
   two-way clustering: set any negative eigenvalues to zero and
   reconstitute @V.
*/

static int maybe_eigenfix (gretl_matrix *V)
{
    gretl_matrix *U = gretl_matrix_copy(V);
    gretl_matrix *lam = NULL;
    int n = V->rows;
    int i, fixit = 0;
    int err = 0;

    if (U == NULL) {
        err = E_ALLOC;
    } else {
        lam = gretl_symmetric_matrix_eigenvals(U, 1, &err);
    }

    for (i=0; i<n && !err; i++) {
        if (lam->val[i] < 0) {
            fixit = 1;
            break;
        }
    }

    if (fixit) {
        gretl_matrix *Tmp = gretl_matrix_copy(U);
        double uij, lvj;
        int j;

        if (Tmp == NULL) {
            err = E_ALLOC;
            goto bailout;
        }
        for (j=0; j<n; j++) {
            lvj = lam->val[j];
            for (i=0; i<n; i++) {
                if (lvj > 0) {
                    uij = gretl_matrix_get(Tmp, i, j);
                    gretl_matrix_set(Tmp, i, j, uij * lvj);
                } else {
                    gretl_matrix_set(Tmp, i, j, 0);
                }
            }
        }
        gretl_matrix_multiply_mod(Tmp, GRETL_MOD_NONE,
                                  U, GRETL_MOD_TRANSPOSE,
                                  V, GRETL_MOD_NONE);
        gretl_matrix_free(Tmp);
    }

 bailout:

    gretl_matrix_free(lam);
    gretl_matrix_free(U);

    return err;
}

static int two_way_cluster_vcv (MODEL *pmod,
				panelmod_t *pan,
				DATASET *pset,
				DATASET *dset,
				const gretl_matrix *XX,
				gretl_matrix *W,
				gretl_matrix *V,
				cluster_info *ci)
{
    gretl_matrix *Tmp;
    int k = V->rows;
    int *v = ci->pcid;
    int err = 0;

    Tmp = gretl_matrix_alloc(k, k);
    if (Tmp == NULL) {
        return E_ALLOC;
    }

    /* prepare for first step */
     ci->target = 0;

    /* Step 1: calculate @V for the first cluster var */
    if (ci->dcid[0] == CI_TIME) {
	err = time_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (ci->dcid[0] == CI_UNIT) {
	err = unit_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else {
	err = generic_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    }

    /* prepare for second step */
    ci->target = 1;

    if (!err && !by_time_and_unit(ci)) {
        /* Step 2a: transcribe the second cluster var if needed */
        err = transcribe_cluster_var(pmod, pan, pset, dset, ci);
    }

    if (!err) {
        /* Step 2b: calculate variance for the second cluster var */
        gretl_matrix_zero(W);
	if (ci->dcid[1] == CI_TIME) {
	    err = time_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
	} else if (ci->dcid[1] == CI_UNIT) {
	    err = unit_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
	} else {
	    err = generic_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
	}
	/* add second variance to first */
	gretl_matrix_add_to(V, Tmp);
    }

    if (!err && by_time_and_unit(ci)) {
	/* Step 3a: calculate variance for combination (= plain White) */
	gretl_matrix_zero(W);
	err = panel_white_vcv(pmod, pan, pset, XX, W, Tmp, ci);
	goto step4;
    }

    if (!err) {
        /* Step 3b: combine the two cluster vars */
        err = combine_categories(pset, v[0], v[1], v[2]);
    }

    if (!err) {
        /* Step 3c: update ci->cz and ci->cvals, then calculate
           variance for the combination
        */
        gretl_matrix_free(ci->cvals);
        ci->cz = pset->Z[ci->pcid[2]];
        ci->cvals = gretl_matrix_values(ci->cz, pmod->nobs, OPT_S, &err);
        gretl_matrix_zero(W);
	ci->target = 2;
        err = generic_cluster_vcv(pmod, pan, pset, XX, W, Tmp, ci);
	ci->target = 0;
    }

 step4:

    if (!err) {
        gretl_matrix_subtract_from(V, Tmp);
	err = maybe_eigenfix(V);
    }

    if (!err) {
        int nc = MIN(ci->nc[0], ci->nc[1]);

	if (by_time_and_unit(ci)) {
	    gretl_model_set_vcv_info(pmod, VCV_PANEL, PANEL_BOTH);
	    gretl_model_set_int(pmod, "n_clusters", nc);
	} else {
	    gretl_model_set_cluster_vcv_info(pmod,
					     pset->varname[v[0]],
					     pset->varname[v[1]]);
	    gretl_model_set_int(pmod, "n_clusters", nc);
	}
    }

    gretl_matrix_free(Tmp);

    return err;
}

/* The method of Driscoll and Kraay in "Consistent Covariance Matrix
   Estimation with Spatially Dependent Panel Data", REStat 80:4, 1998,
   pp. 549-560. See also Daniel Hoechle, "Robust Standard Errors for
   Panel Regressions with Cross-Sectional Dependence", Stata Journal
   7:3, 2007, for discussion of his implementation called xtscc.
*/

static int
driscoll_kraay_vcv (MODEL *pmod, panelmod_t *pan, const DATASET *dset,
                    const gretl_matrix *XX, gretl_matrix *W,
                    gretl_matrix *V)
{
    gretl_matrix_block *B;
    gretl_matrix *H = NULL;
    gretl_matrix *ht = NULL;
    gretl_matrix *htj = NULL;
    gretl_matrix *Wj = NULL;
    gretl_matrix *S = NULL;
    double bw; /* Bartlett weight */
    int *tobs = NULL;
    int k = pmod->ncoeff;
    int i, j, m, vj, s, t;
    int n_c = 0;
    int err = 0;

    tobs = get_panel_tobs(pan);
    if (tobs == NULL) {
        return E_ALLOC;
    }

    B = gretl_matrix_block_new(&H, pan->T, k,
			       &S, k, k,
			       &ht, k, 1,
			       &htj, k, 1,
			       &Wj, k, k,
			       NULL);
    if (B == NULL) {
	free(tobs);
	return E_ALLOC;
    }

    gretl_matrix_zero(H);
    gretl_matrix_zero(S);

    /* Maximum lag for Newey-West. Note: do "set hac_lag nw2"
       for agreement with Stata's xtscc */
    m = get_hac_lag(pan->Tmax);

    /* build the H matrix */
    for (t=0; t<pan->T; t++) {
        int Nt = tobs[t]; /* # of units represented in period @t */
        double eit, htj, hplus;
        int p = 0;

        if (Nt == 0) {
            continue;
        }
	n_c++;
        for (i=0; i<pan->nunits; i++) {
            s = panel_index(i, t);
            eit = pmod->uhat[s];
            if (na(eit)) {
                continue;
            }
            s = small_index(pan, s);
            if (s >= 0) {
                for (j=0; j<k; j++) {
                    /* cumulate xit*eit into row @t of @H */
                    vj = pmod->list[j+2];
                    hplus = dset->Z[vj][s] * eit;
                    htj = gretl_matrix_get(H, t, j);
                    gretl_matrix_set(H, t, j, htj + hplus);
                }
            }
            if (++p == Nt) {
                /* we're done with this period */
                break;
            }
        }
    }

    /* compute initial @S = Omega_0 */
    for (t=0; t<pan->T; t++) {
        for (i=0; i<k; i++) {
            ht->val[i] = gretl_matrix_get(H, t, i);
        }
        gretl_matrix_multiply_mod(ht, GRETL_MOD_NONE,
                                  ht, GRETL_MOD_TRANSPOSE,
                                  S, GRETL_MOD_CUMULATE);
    }

    /* cumulate the weighted cross-lag terms */
    for (j=1; j<=m; j++) {
        gretl_matrix_zero(Wj);
        for (t=j; t<pan->T; t++) {
            for (i=0; i<k; i++) {
                ht->val[i] = gretl_matrix_get(H, t, i);
                htj->val[i] = gretl_matrix_get(H, t-j, i);
            }
            gretl_matrix_multiply_mod(ht, GRETL_MOD_NONE,
                                      htj, GRETL_MOD_TRANSPOSE,
                                      Wj, GRETL_MOD_CUMULATE);
        }
        /* add Barlett weight * (Wj + Wj') to @S */
        bw = 1.0 - j / (m + 1.0);
        gretl_matrix_add_self_transpose(Wj);
        gretl_matrix_multiply_by_scalar(Wj, bw);
        gretl_matrix_add_to(S, Wj);
    }

    finalize_clustered_vcv(pmod, pan, XX, S, V, pan->Tmax);
    gretl_model_set_vcv_info(pmod, VCV_PANEL, PANEL_DK);
    gretl_model_set_hac_order(pmod, m);
    gretl_model_set_int(pmod, "DKT", n_c);

    gretl_matrix_block_destroy(B);
    free(tobs);

    return err;
}

/* Common setup for Arellano, Beck-Katz, Driscoll-Kraay and other
   clustered VCV estimators. For fixed and random effects, @pset is
   the special panel dataset comprising only usable observations, and
   @dset is the original dataset; for pooled OLS @pset and @dset are
   one and the same.
*/

static int
panel_robust_vcv (MODEL *pmod, panelmod_t *pan,
                  DATASET *pset, DATASET *dset)
{
    cluster_info *ci = NULL;
    gretl_matrix *W = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *XX = NULL;
    int k = pmod->ncoeff;
    int pan_robust = 0;
    int err = 0;

    XX = panel_model_xpxinv(pmod, &err);
    if (err) {
        fprintf(stderr, "panel_robust_vcv: failed at panel_model_xpxinv\n");
        goto bailout;
    }

    ci = cluster_info_new(pan);

#if CDEBUG
    fprintf(stderr, "start panel_robust_vcv (pooled = %d)\n", ci->pooled);
    gretl_matrix_print(XX, "X'X^{-1}");
#endif

    if (pan->opt & OPT_C) {
        err = check_cluster_var(pmod, pan, pset, dset, ci);
        if (err) {
            goto bailout;
        }
    } else {
        pan_robust = libset_get_int(PANEL_ROBUST);
    }

    W = gretl_zero_matrix_new(k, k);
    V = gretl_matrix_alloc(k, k);

    if (W == NULL || V == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    if (two_way(ci)) {
	err = two_way_cluster_vcv(pmod, pan, pset, dset, XX, W, V, ci);
    } else if (by_time(ci)) {
        err = time_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (by_unit(ci)) {
	err = unit_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (generic(ci)) {
        err = generic_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    } else if (pan_robust == BECK_KATZ) {
        err = beck_katz_vcv(pmod, pan, pset, XX, W, V);
    } else if (pan_robust == DRISCOLL_KRAAY) {
        driscoll_kraay_vcv(pmod, pan, pset, XX, W, V);
    } else {
	/* default: pan_robust == ARELLANO */
        err = unit_cluster_vcv(pmod, pan, pset, XX, W, V, ci);
    }

#if CDEBUG
    gretl_matrix_print(W, "W");
    gretl_matrix_print(V, "V");
#endif

    if (!err) {
        gretl_model_write_vcv(pmod, V);
    }

 bailout:

    gretl_matrix_free(W);
    gretl_matrix_free(V);
    gretl_matrix_free(XX);
    if (ci != NULL) {
        cluster_info_free(ci);
    }

    if (err && pmod->vcv != NULL) {
        free(pmod->vcv);
        pmod->vcv = NULL;
    }

    return err;
}

static void femod_regular_vcv (MODEL *pmod)
{
    if (pmod->vcv == NULL) {
        /* estimated via Cholesky: no vcv yet */
        makevcv(pmod, pmod->sigma);
    } else {
        /* estimated via QR: "vcv" = (X'X)^{-1} */
        int i, k = pmod->ncoeff;
        int n = k * (k + 1) / 2;
        double s2 = pmod->sigma * pmod->sigma;

        for (i=0; i<n; i++) {
            pmod->vcv[i] *= s2;
        }
    }
}

#if 0 /* not ready yet */

/* See Cameron and Trivedi, Microeconometrics, 21.3.4 */

static int panel_autocorr_1 (MODEL *pmod, const panelmod_t *pan)
{
    double *ubar;
    double *ctt, *css, *cst;
    double uit, uis, rho;
    int i, t, ti;
    int n = 0;

    ubar = malloc(pan->T * sizeof *ubar);
    ctt = malloc(pan->T * sizeof *ctt);
    css = malloc(pan->T * sizeof *css);
    cst = malloc(pan->T * sizeof *cst);

    if (ubar == NULL) {
        return E_ALLOC;
    }

    for (t=0; t<pan->T; t++) {
        ubar[t] = 0.0;
        ctt[t] = css[t] = cst[t] = 0.0;
    }

    for (i=0; i<pan->nunits; i++) {
        if (pan->unit_obs[i] > 0) {
            for (t=0; t<pan->T; t++) {
                ti = panel_index(i, t);
                uit = pmod->uhat[ti];
                if (!na(uit)) {
                    ubar[t] += uit;
                }
            }
            n++;
        }
    }

    for (t=0; t<pan->T; t++) {
        ubar[t] /= n;
    }

    for (i=0; i<pan->nunits; i++) {
        if (pan->unit_obs[i] > 0) {
            for (t=0; t<pan->T; t++) {
                ti = panel_index(i, t);
                uit = pmod->uhat[ti];
                if (t > 0) {
                    uis = pmod->uhat[ti-1];
                } else {
                    uis = NADBL;
                }
                if (!na(uit)) {
                    ctt[t] += (uit - ubar[t]) * (uit - ubar[t]);
                }
                if (!na(uis)) {
                    css[t] += (uis - ubar[t-1]) * (uis - ubar[t-1]);
                }
                if (!na(uis) && !na(uit)) {
                    cst[t] += (uis - ubar[t-1]) * (uit - ubar[t-1]);
                }
            }
        }
    }

    for (t=1; t<pan->T; t++) {
        ctt[t] /= (n - 1);
        css[t] /= (n - 1);
        cst[t] /= (n - 1);
        rho = cst[t] / sqrt(css[t] * ctt[t]);
        fprintf(stderr, "rho(%d) = %g\n", t, rho);
    }

    free(ubar);
    free(ctt);
    free(cst);
    free(css);

    return 0;
}

#endif

#define DWPVAL_TESTING 1

int panel_DW_pval_ok (const MODEL *pmod)
{
    if (na(pmod->dw)) {
        return 0;
    } else {
#if DWPVAL_TESTING
        return 1;
#else
        int Tmax = gretl_model_get_int(pmod, "Tmax");
        int Tmin = gretl_model_get_int(pmod, "Tmin");

        /* Too restrictive? Or not restrictive enough? */
        return Tmax == Tmin;
#endif
    }
}

/* See Bhargava, Franzini and Narendranathan, "Serial Correlation and
   the Fixed Effects Model", Review of Economic Studies 49, 1982,
   page 536. Strictly speaking what's calculated here is the marginal
   significance level of the DW stat when considered as a d_L value.
   It's unclear if this method can really be extended to unbalanced
   panels.
*/

double BFN_panel_DW_pvalue (MODEL *pmod, const DATASET *dset, int *err)
{
    gretl_matrix *lam = NULL;
    double r, pv, lamq, sinarg, pi_2T;
    int T = gretl_model_get_int(pmod, "Tmax");
    int N = gretl_model_get_int(pmod, "n_included_units");
    int nlam, k = pmod->ncoeff;
    int i, q;

    if (pmod->ifc) {
        k--; /* don't include the constant */
    }

    nlam = pmod->nobs - N - k;
    lam = gretl_column_vector_alloc(nlam);
    if (lam == NULL) {
        *err = E_ALLOC;
        return NADBL;
    }

    pi_2T = M_PI / (2.0*T);
    sinarg = sin(pi_2T);
    lamq = 4 * sinarg * sinarg;
    r = pmod->dw;

    q = 1;
    for (i=0; i<nlam; i++) {
        lam->val[i] = lamq - r;
        if ((i+1) % N == 0) {
            q++;
            sinarg = sin(q*pi_2T);
            lamq = 4 * sinarg * sinarg;
        }
    }

    pv = imhof(lam, 0, err);
    if (!*err) {
        if (pv < 0) {
            pv = 0;
        }
        gretl_model_set_double(pmod, "dw_pval", pv);
    }

#if 0
    fprintf(stderr, "DW: T=%d, Tmin=%d, N=%d, nlam=%d, DW=%g, pv=%g\n",
            T, gretl_model_get_int(pmod, "Tmin"), N, nlam, pmod->dw, pv);
#endif

    gretl_matrix_free(lam);

    return pv;
}

/* Durbin-Watson statistic for the pooled or fixed effects model.

   See Bhargava, Franzini and Narendranathan, "Serial Correlation and
   the Fixed Effects Model", Review of Economic Studies 49, 1982,
   pp. 533-549, and also Baltagi and Wu, "Unequally Spaced Panel Data
   Regressions With AR(1) Disturbances", Econometric Theory, 15, 1999,
   pp. 814-823 for discussion of the unbalanced case.
*/

static void panel_dwstat (MODEL *pmod, panelmod_t *pan)
{
    double dwnum = 0.0;
    double rnum = 0.0;
    double rden = 0.0;
    double ut, u1;
    int in_bounds = 1;
    int i, t, s;

    pmod->dw = pmod->rho = NADBL;

    if (pmod->ess <= 0.0) {
        return;
    }

#if 0
    fprintf(stderr, "DW: pmod=%p, pooled=%p, T=%d\n",
            (void *) pmod, (void *) pan->pooled, pan->T);
    fprintf(stderr, " t1=%d, t2=%d, nobs=%d, full_n=%d\n", pmod->t1,
            pmod->t2, pmod->nobs, pmod->full_n);
#endif

    for (i=0; i<pan->nunits && in_bounds; i++) {
        int started = 0;

        if (pan->unit_obs[i] == 0) {
            continue;
        }
        for (t=0; t<pan->T; t++) {
            s = panel_index(i, t);
            if (s >= pmod->t2) {
                in_bounds = 0;
                break;
            }
            ut = pmod->uhat[s];
            if (!na(ut)) {
                if (t == 0) {
                    started = 1;
                    continue;
                }
                u1 = pmod->uhat[s-1];
                if (na(u1)) {
                    /* implicitly take u1 as 0 */
                    if (started) {
                        dwnum += ut * ut;
                    }
                } else {
                    dwnum += (ut - u1) * (ut - u1);
                    rnum += ut * u1;
                    rden += u1 * u1;
                }
                started = 1;
            }
        }
    }

    if (dwnum > 0.0) {
        pmod->dw = dwnum / pmod->ess;
    }
    if (na(pmod->rho) && rden > 0.0 && !na(rden)) {
        pmod->rho = rnum / rden;
    }
    if (!na(pmod->rho) && (pmod->rho <= -1.0 || pmod->rho >= 1.0)) {
        pmod->rho = NADBL;
    }

    if (pan->opt & OPT_U) {
        /* make these stats available for random effects reporting */
        pan->dw = pmod->dw;
        pan->rho = pmod->rho;
    }
}

/* Allocate the arrays needed to perform the Hausman test,
   in its matrix formulation.
*/

static int matrix_hausman_allocate (panelmod_t *pan)
{
    int k = pan->vlist[0] - 2;
    int err = 0;

    /* array to hold differences between coefficient estimates */
    pan->bdiff = gretl_vector_alloc(k);
    if (pan->bdiff == NULL) {
        err = E_ALLOC;
    } else {
        /* array to hold covariance matrix */
        pan->Sigma = gretl_matrix_alloc(k, k);
        if (pan->Sigma == NULL) {
            gretl_matrix_free(pan->bdiff);
            pan->bdiff = NULL;
            err = E_ALLOC;
        }
    }

    if (!err) {
        pan->dfH = k;
    }

    return err;
}

static int *real_varying_list (panelmod_t *pan)
{
    int *vlist = gretl_list_copy(pan->vlist);
    int i;

    if (vlist != NULL) {
        for (i=2; i<=vlist[0]; i++) {
            if (vlist[i] == 0) {
                gretl_list_delete_at_pos(vlist, i);
                break;
            }
        }
    }

    return vlist;
}

/* Construct a version of the dataset from which the group means are
   subtracted, for the "within" regression.  Nota bene: this auxiliary
   dataset is not necessarily of full length: missing observations
   are skipped.
*/

static DATASET *within_groups_dataset (const DATASET *dset,
                                       panelmod_t *pan)
{
    DATASET *wset = NULL;
    int *vlist = NULL;
    int i, j, vj, nv;
    int s, t, bigt;
    int err = 0;

    pan->balanced = 1;

    for (i=0; i<pan->nunits; i++) {
        if (pan->unit_obs[i] > 0) {
            if (pan->unit_obs[i] != pan->Tmax) {
                pan->balanced = 0;
            }
        }
    }

    if (pan->NT < dset->n) {
	/* allow for panel obs book-keeping */
        err = allocate_data_finders(pan, dset->n);
        if (err) {
            return NULL;
        }
    }

    vlist = real_varying_list(pan);
    if (vlist == NULL) {
        return NULL;
    }

    nv = pan->vlist[0];

#if PDEBUG
    fprintf(stderr, "within_groups dataset: nvars=%d, nobs=%d\n",
            pan->vlist[0], pan->NT);
#endif

    wset = create_auxiliary_dataset(nv, pan->NT, 0);
    if (wset == NULL) {
        free(vlist);
        return NULL;
    }

    for (j=1; j<=vlist[0]; j++) {
        double xbar, gxbar = 0.0;
        int allzero = 1;

        vj = vlist[j];
        gxbar = 0.0;
        s = 0;

#if PDEBUG
        strcpy(wset->varname[j], dset->varname[vj]);
        fprintf(stderr, "de-meaning: working on list[%d], %s\n",
                j, dset->varname[vj]);
#endif

        for (i=0; i<pan->nunits; i++) {
            int Ti = pan->unit_obs[i];
            int grpzero = 1;
            int got = 0;

            if (Ti == 0) {
                continue;
            }

            /* first pass: find the group mean */
            xbar = 0.0;
            for (t=0; t<pan->T; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    xbar += dset->Z[vj][bigt];
                }
            }

            gxbar += xbar;
            xbar /= Ti;

            /* second pass: calculate de-meaned values */
            got = 0;
            for (t=0; t<pan->T && got<Ti; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    wset->Z[j][s] = dset->Z[vj][bigt] - xbar;
                    if (wset->Z[j][s] != 0.0) {
                        grpzero = 0;
                    }
                    got++;
                    if (pan->small2big != NULL) {
                        pan->small2big[s] = bigt;
                        pan->big2small[bigt] = s;
                    }
                    s++;
                }
            }

            if (!grpzero) {
                allzero = 0;
            }
        } /* end loop over units */

        gxbar /= pan->NT;

        if (j == 1 && allzero) {
            /* the dependent variable is not time-varying */
            fprintf(stderr, "fixed effects: dependent var is time-invariant\n");
        }

        /* wZ = data - group mean + grand mean */
        for (s=0; s<pan->NT; s++) {
            wset->Z[j][s] += gxbar;
        }
    }

    free(vlist);

    return wset;
}

/* Construct a quasi-demeaned version of the dataset so we can apply
   least squares to estimate the random effects model.  This dataset
   is not necessarily of full length.  If we're implementing the
   Hausman test using the regression approach this dataset will also
   include "straight" de-meaned counterparts of the time-varying
   variables.
*/

static DATASET *
random_effects_dataset (const DATASET *dset,
                        const DATASET *gset,
                        int *relist,
                        int *hlist,
                        panelmod_t *pan)
{
    DATASET *rset;
    double xbar, theta_i;
    int hreg = (hlist != NULL);
    int v1 = relist[0];
    int v2 = 0;
    int i, j, k, k2, t;
    int vj, s, bigt, u;
    int err = 0;

    if (hreg) {
        /* Apparatus for regression version of Hausman test:
           note that we shouldn't include time dummies here
           since the de-meaned versions would be perfectly
           collinear with the quasi-demeaned ones.
        */
        int hmax = pan->vlist[0] - pan->ntdum;

        for (i=1; i<hmax; i++) {
            if (pan->vlist[i+1] != 0) {
                hlist[0] = ++v2;
                hlist[i-1] = v1 + i - 2;
            }
        }
    }

    if (pan->NT < dset->n) {
        err = allocate_data_finders(pan, dset->n);
        if (err) {
            return NULL;
        }
    }

#if PDEBUG
    fprintf(stderr, "random_effects_dataset: nvars=%d, nobs=%d\n",
            v1 + v2, pan->NT);
#endif

    rset = create_auxiliary_dataset(v1 + v2, pan->NT, 0);
    if (rset == NULL) {
        return NULL;
    }

    /* build GLS regression list, and process varnames for regression
       version of Hausman test if wanted
    */

    k = 0;
    k2 = v1 - 1;
    for (j=1; j<=v1; j++) {
        vj = pan->pooled->list[j];
        if (vj == 0) {
            relist[j] = 0;
        } else {
            relist[j] = ++k;
            strcpy(rset->varname[k], dset->varname[vj]);
            if (hreg && k2 < rset->v - 1 && j > 1 && var_is_varying(pan, vj)) {
                /* hausman-related term */
                k2++;
                strcpy(rset->varname[k2], "_");
                strncat(rset->varname[k2], dset->varname[vj],
                        VNAMELEN - 2);
            }
        }
    }

    /* Now create the transformed variables: original data minus theta
       times the appropriate group mean (and in addition, original
       data minus the group mean for time-varying vars, if doing the
       Hausman test by the regression method).
    */

    s = u = 0;
    pan->theta_bar = 0.0;

    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];
        int got = 0;

        if (Ti == 0) {
            continue;
        }

        if (Ti != pan->Tmax) {
            theta_i = 1.0 - sqrt(pan->s2e / (Ti * pan->s2v + pan->s2e));
        } else {
            theta_i = pan->theta;
        }

        pan->theta_bar += theta_i;

        for (t=0; t<pan->T && got<Ti; t++) {
            bigt = panel_index(i, t);
            k = 0;
            k2 = v1 - 1;
            if (!panel_missing(pan, bigt)) {
                for (j=0; j<v1; j++) {
                    vj = pan->pooled->list[j+1];
                    if (vj == 0) {
                        rset->Z[0][s] -= theta_i;
                    } else {
                        k++;
                        xbar = (k < gset->v)? gset->Z[k][u] : 1.0 / pan->Tmax;
                        rset->Z[k][s] = dset->Z[vj][bigt] - theta_i * xbar;
                        if (hreg && k2 < rset->v - 1 && var_is_varying(pan, vj)) {
                            /* hausman-related term */
                            rset->Z[++k2][s] = dset->Z[vj][bigt] - xbar;
                        }
                    }
                }
                got++;
                if (pan->small2big != NULL) {
                    pan->small2big[s] = bigt;
                    pan->big2small[bigt] = s;
                }
                s++;
            }
        }
        u++;
    }

    return rset;
}

#define BETWEEN_DEBUG 0

/* Construct a "mini-dataset" containing the group means, in
   order to run the group-means or "between" regression.
*/

static DATASET *group_means_dataset (panelmod_t *pan,
                                     const DATASET *dset)
{
    DATASET *gset;
    double x;
    int gn = pan->effn;
    int gv = pan->pooled->list[0];
    int i, j, k;
    int s, t, bigt;

    if (pan->balanced && pan->ntdum > 0) {
        gv -= pan->ntdum;
    }

#if PDEBUG
    fprintf(stderr, "group_means_dataset: nvars=%d, nobs=%d\n",
            gv, gn);
#endif

    gset = create_auxiliary_dataset(gv, gn, 0);
    if (gset == NULL) {
        return NULL;
    }

    k = 1;
    for (j=1; j<=gv; j++) {
        int vj = pan->pooled->list[j];

        if (vj == 0) {
            continue;
        }

#if BETWEEN_DEBUG
        strcpy(ginfo->varname[k], dset->varname[vj]);
#else
        if (pan->opt & OPT_B) {
            /* will save the "between" model: so name the variables */
            strcpy(gset->varname[k], dset->varname[vj]);
        }
#endif

        s = 0;
        for (i=0; i<pan->nunits; i++) {
            int Ti = pan->unit_obs[i];

            if (Ti == 0) {
                continue;
            }

            x = 0.0;
            for (t=0; t<pan->T; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    x += dset->Z[vj][bigt];
                }
            }
            gset->Z[k][s++] = x / Ti;
        }
        k++;
    }

#if BETWEEN_DEBUG
    if (1) {
        PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

        printdata(NULL, NULL, gset, OPT_O, prn);
        gretl_print_destroy(prn);
    }
#endif

    return gset;
}

/* spruce up the between model and attach it to pan */

static int save_between_model (MODEL *pmod, const int *blist,
                               DATASET *gset, panelmod_t *pan)
{
    gretl_matrix *uh, *yh;
    int i, err = 0;

    pmod->ci = PANEL;
    pmod->opt |= OPT_B;
    pmod->dw = NADBL;
    gretl_model_add_panel_varnames(pmod, gset, NULL);

    uh = gretl_column_vector_alloc(pmod->nobs);
    yh = gretl_column_vector_alloc(pmod->nobs);

    if (uh == NULL || yh == NULL) {
        err = E_ALLOC;
    } else {
        for (i=0; i<pmod->nobs; i++) {
            uh->val[i] = pmod->uhat[i];
            yh->val[i] = pmod->yhat[i];
        }
        gretl_model_set_matrix_as_data(pmod, "uhat", uh);
        gretl_model_set_matrix_as_data(pmod, "yhat", yh);
    }

#if 0
    /* this is risky at present: too many functions want
       to read pmod->uhat directly */
    free(pmod->uhat); free(pmod->yhat);
    pmod->uhat = pmod->yhat = NULL;
#endif

    *pan->realmod = *pmod;

    return err;
}

/* Compute @ubPub as the Ti-weighted sum of the squared
   residuals from the Between model, as per Stata (but
   in disagreement with Baltagi and Chang, 1994).
*/

static int alt_compute_ubPub (panelmod_t *pan, MODEL *bmod)
{
    int i, Ti, t = 0;

    pan->ubPub = 0.0;

    for (i=0; i<pan->nunits; i++) {
        Ti = pan->unit_obs[i];
        if (Ti > 0) {
            pan->ubPub += Ti * bmod->uhat[t] * bmod->uhat[t];
            t++;
        }
    }

    return 0;
}

static void adjust_gset_data (panelmod_t *pan, DATASET *gset,
                              int step)
{
    int i, j, Ti, t = 0;
    double adj;

    for (i=0; i<pan->nunits; i++) {
        Ti = pan->unit_obs[i];
        if (Ti > 0) {
            adj = step == 0 ? sqrt(Ti) : 1.0/sqrt(Ti);
            for (j=0; j<gset->v; j++) {
                gset->Z[j][t] *= adj;
            }
            t++;
        }
    }
}

/* Compute @ubPub as the sum of squared residuals from a
   Ti-weighted Between regression, as per Baltagi and Chang,
   1994, and also Baltagi, 2013.
*/

static int compute_ubPub (panelmod_t *pan, MODEL *bmod,
                          int *blist, DATASET *gset)
{
    int err;

    /* multiply all data by sqrt(Ti) */
    adjust_gset_data(pan, gset, 0);
    clear_model(bmod);
    *bmod = lsq(blist, gset, OLS, OPT_A);
    err = bmod->errcode;
    if (!err) {
        pan->ubPub = bmod->ess;
    }
    /* put the original data back */
    adjust_gset_data(pan, gset, 1);

    return err;
}

/* calculate the group means or "between" regression and its error
   variance */

static int between_variance (panelmod_t *pan, DATASET *gset)
{
    gretlopt bopt;
    MODEL bmod;
    int *blist;
    int i, j, k;
    int err = 0;

    blist = gretl_list_new(gset->v);
    if (blist == NULL) {
        return E_ALLOC;
    }

    j = k = 1;
    for (i=1; i<=gset->v; i++) {
        if (pan->pooled->list[i] == 0) {
            blist[k++] = 0;
        } else {
            blist[k++] = j++;
        }
    }

    bopt = (pan->opt & OPT_B)? OPT_NONE : OPT_A;
    bmod = lsq(blist, gset, OLS, bopt);

    if (bmod.errcode == 0) {
        pan->s2b = bmod.ess / (bmod.nobs - bmod.ncoeff);
    } else {
        err = bmod.errcode;
    }

#if PDEBUG
    if (err) {
        fprintf(stderr, "error %d in between_variance\n", err);
    } else {
        fprintf(stderr, "pan->s2b = %g\n", pan->s2b);
    }
#endif

    if (!err && (pan->opt & OPT_B)) {
        err = save_between_model(&bmod, blist, gset, pan);
    } else {
        if (!err && !pan->balanced && (pan->opt & OPT_U) &&
            (pan->opt & OPT_X) && !(pan->opt & OPT_E)) {
            /* Prepare for the Baltagi-Chang take on Swamy-Arora
               in the case of an unbalanced panel
            */
            if (stata_sa) {
                err = alt_compute_ubPub(pan, &bmod);
            } else {
                err = compute_ubPub(pan, &bmod, blist, gset);
            }
        }
        clear_model(&bmod);
    }

    free(blist);

    return err;
}

/* op is VCV_SUBTRACT for the random effects model.  With the fixed
   effects model we don't have to worry about excluding vcv entries
   for time-invariant variables, since none of these are included.
*/

static int
vcv_skip (const MODEL *pmod, int i, const panelmod_t *pan, int op)
{
    int skip = 0;

    if (pmod->list[i+2] == 0) {
        /* always skip the constant */
        skip = 1;
    } else if (op == VCV_SUBTRACT && !pan->varying[i]) {
        /* random effects, time-invariant var */
        skip = 1;
    }

    return skip;
}

/* Fill out the covariance matrix for use with the Hausman test:
   the entries that get transcribed here are only those for
   slopes with respect to time-varying variables.
*/

static void
vcv_slopes (panelmod_t *pan, const MODEL *pmod, int op)
{
    int idx, i, j;
    int mj, mi = 0;
    int sj, si = 0;
    double x;

    for (i=0; i<pan->dfH; i++) {
        if (vcv_skip(pmod, mi, pan, op)) {
            i--;
            mi++;
            continue;
        }
        mj = mi;
        sj = si;
        for (j=i; j<pan->dfH; j++) {
            if (vcv_skip(pmod, mj, pan, op)) {
                j--;
                mj++;
                continue;
            }

            idx = ijton(mi, mj, pmod->ncoeff);

            if (op == VCV_SUBTRACT) {
                x = gretl_matrix_get(pan->Sigma, si, sj);
                x -= pmod->vcv[idx];
            } else {
                x = pmod->vcv[idx];
            }

            gretl_matrix_set(pan->Sigma, si, sj, x);
            if (si != sj) {
                gretl_matrix_set(pan->Sigma, sj, si, x);
            }

            mj++;
            sj++;
        }
        mi++;
        si++;
    }
}

/* calculate Hausman test statistic, matrix diff style */

static int bXb (panelmod_t *pan)
{
    int err;

    err = gretl_invert_symmetric_matrix(pan->Sigma);

    if (!err) {
        pan->H = gretl_scalar_qform(pan->bdiff, pan->Sigma, &err);
    }

    if (err) {
        pan->H = NADBL;
    }

    return err;
}

static void fixed_effects_df_correction (MODEL *pmod, int k)
{
    double dfcorr = sqrt((double) pmod->dfd / (pmod->dfd - k));
    int i;

    pmod->dfd -= k;
    pmod->dfn += k;

    pmod->sigma *= dfcorr;

    for (i=0; i<pmod->ncoeff; i++) {
        pmod->sderr[i] *= dfcorr;
    }
}

/* used for printing fixed- or random-effects estimates
   in the context of the "panel diagnostics" routine
*/

static void simple_print_panel_model (MODEL *pmod,
                                      DATASET *dset,
                                      const int *xlist,
                                      PRN *prn)
{
    int *savelist = gretl_list_copy(pmod->list);
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
        pmod->list[i+2] = xlist[i+2];
    }
    printmodel(pmod, dset, OPT_P, prn);
    free(pmod->list);
    pmod->list = savelist;
}

static void print_re_results (panelmod_t *pan,
                              MODEL *pmod,
                              DATASET *dset,
                              PRN *prn)
{
    pputs(prn, _("Variance estimators:"));
    pputc(prn, '\n');
    pprintf(prn, _(" between = %g"), pan->s2v);
    pputc(prn, '\n');
    pprintf(prn, _(" within = %g"), pan->s2e);
    pputc(prn, '\n');

    if (pan->balanced || pan->s2v == 0) {
        pprintf(prn, _("theta used for quasi-demeaning = %g"), pan->theta);
        pputc(prn, '\n');
    } else {
        pputs(prn, _("Panel is unbalanced: theta varies across units"));
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    pputs(prn, _("Random effects estimator\n"
                 "allows for a unit-specific component to the error term\n"));
    pputc(prn, '\n');

    simple_print_panel_model(pmod, dset, pan->pooled->list, prn);
}

static int print_fe_results (panelmod_t *pan,
                             MODEL *pmod,
                             DATASET *dset,
                             PRN *prn)
{
    int dfn = pan->effn - 1;

    pputs(prn, _("Fixed effects estimator\n"
                 "allows for differing intercepts by cross-sectional unit\n"));
    pputc(prn, '\n');

    simple_print_panel_model(pmod, dset, pan->vlist, prn);

    pprintf(prn, _("Residual variance: %g/(%d - %d) = %g\n"),
            pmod->ess, pmod->nobs, pan->vlist[0] - 1 + dfn, pan->s2e);
    pputc(prn, '\n');

    if (!na(pan->Ffe)) {
        pprintf(prn, _("Joint significance of differing group means:\n"));
        pprintf(prn, " F(%d, %g) = %g %s %g\n", pan->Fdfn, pan->Fdfd, pan->Ffe,
                _("with p-value"), snedecor_cdf_comp(pan->Fdfn, pan->Fdfd, pan->Ffe));
        pputs(prn, _("(A low p-value counts against the null hypothesis that "
                     "the pooled OLS model\nis adequate, in favor of the fixed "
                     "effects alternative.)\n\n"));
    } else {
        pputc(prn, '\n');
    }

    return 0;
}

static int time_dummies_wald_test (panelmod_t *pan, MODEL *pmod)
{
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    double x;
    int i, j, k, bigk;
    int di, dj;
    int err;

    if (pan->ntdum == 0) {
        return 0;
    }

    k = pan->ntdum;
    bigk = pmod->ncoeff;

    if (pmod->vcv == NULL) {
        err = makevcv(pmod, pmod->sigma);
        if (err) {
            return err;
        }
    }

    b = gretl_column_vector_alloc(k);
    vcv = gretl_matrix_alloc(k, k);
    if (b == NULL || vcv == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    di = bigk - k;
    for (i=0; i<k; i++) {
        b->val[i] = pmod->coeff[di++];
    }

    di = bigk - k;
    for (i=0; i<k; i++) {
        dj = bigk - k;
        for (j=0; j<=i; j++) {
            x = pmod->vcv[ijton(di, dj++, bigk)];
            gretl_matrix_set(vcv, i, j, x);
            gretl_matrix_set(vcv, j, i, x);
        }
        di++;
    }

    err = gretl_invert_symmetric_matrix(vcv);
    if (err) {
        goto bailout;
    }

    x = gretl_scalar_qform(b, vcv, &err);
    if (err) {
        fprintf(stderr, "Failed to compute test statistic\n");
        goto bailout;
    }

    if (!err) {
        ModelTest *test = model_test_new(GRETL_TEST_PANEL_TIMEDUM);

        if (test != NULL) {
            model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
            model_test_set_dfn(test, k);
            model_test_set_value(test, x);
            model_test_set_pvalue(test, chisq_cdf_comp(k, x));
            maybe_add_test_to_model(pmod, test);
        }
    }

 bailout:

    gretl_matrix_free(vcv);
    gretl_vector_free(b);

    return err;
}

static void save_fixed_effects_F (panelmod_t *pan, MODEL *wmod)
{
    int robust = (pan->opt & OPT_R);
    ModelTest *test;

    if (na(pan->Ffe)) {
        return;
    }

    test = model_test_new(robust ? GRETL_TEST_PANEL_WELCH :
                          GRETL_TEST_PANEL_F);

    if (test != NULL) {
        model_test_set_teststat(test, robust ? GRETL_STAT_WF : GRETL_STAT_F);
        model_test_set_dfn(test, pan->Fdfn);
        model_test_set_dfd(test, pan->Fdfd);
        model_test_set_value(test, pan->Ffe);
        model_test_set_pvalue(test, snedecor_cdf_comp(pan->Fdfn, pan->Fdfd, pan->Ffe));
        maybe_add_test_to_model(wmod, test);
    }
}

/* "regular" = not robust: sums-of-squares based joint test on
   the fixed effects */

static void regular_fixed_effects_F (panelmod_t *pan, MODEL *wmod)
{
    int k_pooled = pan->pooled->list[0];
    int k_fe = pan->vlist[0];

    pan->Fdfn = pan->effn - 1;
    pan->Fdfd = wmod->dfd;

    if (k_pooled > k_fe) {
        pan->Fdfn -= k_pooled - k_fe;
        if (pan->Fdfn <= 0) {
            pan->Ffe = NADBL;
            return;
        }
    }

    pan->Ffe = (pan->pooled->ess - wmod->ess) * pan->Fdfd /
        (wmod->ess * pan->Fdfn);
    if (na(pan->Ffe)) {
        pan->Ffe = NADBL;
    } else if (pan->Ffe < 0) {
        pan->Ffe = 0;
    }
}

static int fix_panelmod_list (MODEL *targ, panelmod_t *pan)
{
    int i;

#if PDEBUG
    printlist(targ->list, "targ->list");
    printlist(pan->pooled->list, "pan->pooled->list");
    printlist(pan->vlist, "pan->vlist");
#endif

    free(targ->list);
    targ->list = gretl_list_copy(pan->pooled->list);
    if (targ->list == NULL) {
        return E_ALLOC;
    }

    if (pan->opt & OPT_F) {
        /* fixed effects: remove any non-varying variables */
        for (i=2; i<=targ->list[0]; i++) {
            if (!in_gretl_list(pan->vlist, targ->list[i])) {
                gretl_list_delete_at_pos(targ->list, i--);
            }
        }
    }

#if PDEBUG
    printlist(targ->list, "new targ->list");
#endif

    return 0;
}

/* Compute F-test or chi-square test for the regular
   regressors, skipping @nskip trailing coefficients
   (which can be used to skip time dummies)
*/

static double panel_overall_test (MODEL *pmod, panelmod_t *pan,
                                  int nskip, gretlopt opt)
{
    double test = NADBL;
    int *omitlist = NULL;

    if (pmod->ncoeff == 1) {
        return test;
    }

    if (nskip > 0) {
        int i, k = pmod->list[0] - nskip - 2;

        omitlist = gretl_list_new(k);
        for (i=1; i<=k; i++) {
            omitlist[i] = pmod->list[i+2];
        }
    }

    if (opt & OPT_X) {
        test = wald_omit_chisq(omitlist, pmod);
    } else {
        test = wald_omit_F(omitlist, pmod);
    }

    free(omitlist);

    return test;
}

/* Correct various model statistics, in the case where we estimated
   the fixed effects or "within" model on an auxiliary dataset
   from which the group means were subtracted.
*/

static int fix_within_stats (MODEL *fmod, panelmod_t *pan)
{
    double wrsq, wfstt = NADBL;
    int wdfn, nc = fmod->ncoeff;
    int err = 0;

    err = fix_panelmod_list(fmod, pan);
    if (err) {
        return err;
    }

    fmod->ybar = pan->pooled->ybar;
    fmod->sdy = pan->pooled->sdy;
    fmod->tss = pan->pooled->tss;
    fmod->ifc = 1;

    wrsq = fmod->rsq;
    if (wrsq < 0.0) {
        wrsq = 0.0;
    } else if (na(wrsq)) {
        wrsq = NADBL;
    }

    /* Should we differentiate "regular" regressors from
       time dummies, if included? For now, yes.
    */

    if (pan->ntdum > 0) {
        wdfn = fmod->ncoeff - 1 - pan->ntdum;
        if (wdfn > 0) {
            wfstt = panel_overall_test(fmod, pan, pan->ntdum,
                                       OPT_NONE);
        }
    } else {
        wdfn = fmod->ncoeff - 1;
        if (wdfn > 0) {
            if (pan->opt & OPT_R) {
                wfstt = panel_overall_test(fmod, pan, 0, OPT_NONE);
            } else {
                wfstt = (wrsq / (1.0 - wrsq)) * ((double) fmod->dfd / wdfn);
            }
        }
    }

    if (!na(wfstt) && wfstt >= 0.0) {
        ModelTest *test = model_test_new(GRETL_TEST_WITHIN_F);
        int nc, wdfd = fmod->dfd;

        if (pan->opt & OPT_C) {
            nc = gretl_model_get_int(fmod, "n_clusters");
            fmod->dfd = wdfd = nc - 1;
        } else if (pan->opt & OPT_R) {
	    nc = gretl_model_get_int(fmod, "DKT");
	    if (nc > 0) {
		/* Driscoll-Kraay */
		fmod->dfd = wdfd = nc - 1;
	    } else {
		fmod->dfd = wdfd = pan->effn - 1;
	    }
        }

        if (test != NULL) {
            model_test_set_teststat(test, GRETL_STAT_F);
            model_test_set_dfn(test, wdfn);
            model_test_set_dfd(test, wdfd);
            model_test_set_value(test, wfstt);
            model_test_set_pvalue(test, snedecor_cdf_comp(wdfn, wdfd, wfstt));
            maybe_add_test_to_model(fmod, test);
        }
    }

    /* note: this member is being borrowed for the "Within R-squared" */
    fmod->adjrsq = wrsq;

    /* LSDV-based statistics (FIXME: can we do this in
       the --robust case?)
    */
    if (pan->opt & OPT_R) {
        fmod->rsq = 1.0 - (fmod->ess / fmod->tss);
        fmod->fstt = NADBL;
    } else {
        fmod->rsq = 1.0 - (fmod->ess / fmod->tss);
        if (fmod->rsq < 0.0) {
            fmod->rsq = 0.0;
        } else {
            fmod->fstt = (fmod->rsq / (1.0 - fmod->rsq)) *
                ((double) fmod->dfd / fmod->dfn);
        }
    }

    fmod->ncoeff = fmod->dfn + 1; /* number of params estimated */
    ls_criteria(fmod);
    fmod->ncoeff = nc;

    return err;
}

/* Fixed-effects model: add the per-unit intercept estimates
   to the model in case the user wants to retrieve them.
   Random-effects model: add estimate of the individual effects.

   By this point the model -- even if it been estimated on a short
   dataset -- should have a full-length residual series.
*/

static int panel_model_add_ahat (MODEL *pmod, const DATASET *dset,
                                 panelmod_t *pan)
{
    double *ahat = NULL;
    const double *x;
    int i, j, t, bigt;
    int n, err = 0;

    n = pmod->full_n;

    ahat = malloc(n * sizeof *ahat);
    if (ahat == NULL) {
        return E_ALLOC;
    }

    for (t=0; t<n; t++) {
        ahat[t] = NADBL;
    }

    if (pan->opt & OPT_F) {
        /* fixed effects */
        for (i=0; i<pan->nunits; i++) {
            if (pan->unit_obs[i] > 0) {
                double a = 0.0;

                /* a = y - Xb, where the 'b' is based on de-meaned data */

                for (t=0; t<pan->T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod->uhat[bigt])) {
                        a += dset->Z[pmod->list[1]][bigt];
                        for (j=1; j<pmod->ncoeff; j++) {
                            x = dset->Z[pmod->list[j+2]];
                            a -= pmod->coeff[j] * x[bigt];
                        }
                    }
                }

                a /= pan->unit_obs[i];

                for (t=0; t<pan->T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod->uhat[bigt])) {
                        ahat[bigt] = a;
                    }
                }
            }
        }
    } else {
        /* random effects */
        double uhbar, frac = 0;

        if (pan->balanced) {
            frac = 1.0 - pan->theta;
            frac = 1.0 - frac * frac;
        }

        for (i=0; i<pan->nunits; i++) {
            if (pan->unit_obs[i] > 0) {
                /* get mean residual */
                uhbar = 0.0;
                for (t=0; t<pan->T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod->uhat[bigt])) {
                        uhbar += pmod->uhat[bigt];
                    }
                }
                uhbar /= pan->unit_obs[i];
                /* ahat = frac * uhbar */
                if (!pan->balanced) {
                    frac = pan->s2v / (pan->s2v + pan->s2e / pan->unit_obs[i]);
                }
                for (t=0; t<pan->T; t++) {
                    bigt = panel_index(i, t);
                    if (!na(pmod->uhat[bigt])) {
                        ahat[bigt] = frac * uhbar;
                    }
                }
            }
        }
    }

    err = gretl_model_set_data(pmod, "ahat", ahat,
                               GRETL_TYPE_DOUBLE_ARRAY,
                               n * sizeof *ahat);

    return err;
}

/* If we're estimating the fixed effects model with the --robust
   flag, we should do a robust version of the joint test on
   the fixed effects. Here we use the algorithm of B. L. Welch,
   "On the Comparison of Several Mean Values: An Alternative
   Approach" (Biometrika 38, 1951, pp. 330-336). The variable
   we're testing for difference of means (by individual) is the
   residual from pooled OLS.
*/

static int robust_fixed_effects_F (panelmod_t *pan)
{
    MODEL *pmod = pan->pooled;
    double *u, *w, *h, *xbar;
    double muhat = 0.0;
    double x, s2, W = 0.0;
    double A, B, sum_h;
    int k = pan->effn;
    int i, j, t, s;
    int Ti, bigt;
    int err = 0;

    u = malloc(pan->Tmax * sizeof *u);
    w = malloc(3 * k * sizeof *w);

    if (u == NULL || w == NULL) {
        free(u);
        free(w);
        return E_ALLOC;
    }

    h = w + k;
    xbar = h + k;

    j = 0;
    for (i=0; i<pan->nunits; i++) {
        Ti = pan->unit_obs[i];
        if (Ti > 1) {
            s = 0;
            for (t=0; t<pan->T; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    u[s++] = pmod->uhat[bigt];
                }
            }
            xbar[j] = gretl_mean(0, s-1, u);
            s2 = gretl_variance(0, s-1, u);
            w[j] = Ti / s2;
            W += w[j];
            muhat += w[j] * xbar[j];
            j++;
        }
    }

    muhat /= W;
    A = sum_h = 0.0;

    j = 0;
    for (i=0; i<pan->nunits; i++) {
        Ti = pan->unit_obs[i];
        if (Ti > 1) {
            x = 1 - w[j]/W;
            h[j] = (x * x) / (Ti - 1);
            sum_h += h[j];
            x = xbar[j] - muhat;
            A += w[j] * x * x;
            j++;
        }
    }

    A /= (k - 1);
    B = 1.0 + (2.0*(k-2.0)/(k * k - 1.0)) * sum_h;

    pan->Ffe = A / B;
    pan->Fdfn = k - 1;
    pan->Fdfd = (k * k - 1.0)/(3 * sum_h);

    free(u);
    free(w);

    return err;
}

static int *real_FE_list (panelmod_t *pan)
{
    int *list = gretl_list_copy(pan->pooled->list);
    int i;

    if (list != NULL) {
        /* purge any non-time-varying variables */
        for (i=2; i<=list[0]; i++) {
            if (!in_gretl_list(pan->vlist, list[i])) {
                gretl_list_delete_at_pos(list, i--);
            }
        }
    }

    return list;
}

/* For testing purposes: retrieve user-specified values for
   s2v and s2e. These would typically be known good values
   in the context of a simulation.
*/

static int read_true_variances (panelmod_t *pan)
{
    gretl_matrix *m;
    int err = 0;

    m = gretl_matrix_read_from_file(glsmat, 0, &err);
    if (m == NULL) {
        fprintf(stderr, "read_true_variances: no matrix!\n");
        if (!err) {
            err = E_DATA;
        }
    } else {
        pan->s2v = m->val[0];
        pan->s2e = m->val[1];
    }

    return err;
}

/* computation of $\hat{\sigma}^2_v$ a la Nerlove, if wanted */

static int nerlove_s2v (MODEL *pmod, const DATASET *dset,
                        panelmod_t *pan)
{
    double a, amean, *ahat;
    double wmean, *wi = NULL;
    const double *x;
    int i, j, t, k, bigt;
    int *list = NULL;
    int err = 0;

    ahat = malloc(pan->effn * sizeof *ahat);
    if (ahat == NULL) {
        return E_ALLOC;
    }

    if (pan->opt & OPT_X) {
        /* --unbalanced */
        wi = malloc(pan->effn * sizeof *wi);
        if (wi == NULL) {
            free(ahat);
            return E_ALLOC;
        }
    }

    list = real_FE_list(pan);
    if (list == NULL) {
        free(ahat);
        free(wi);
        return E_ALLOC;
    }

    wmean = amean = 0.0;
    k = 0;

    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];

        if (Ti > 0) {
            a = 0.0;
            for (t=0; t<pan->T; t++) {
                bigt = panel_index(i, t);
                if (!na(pan->pooled->uhat[bigt])) {
                    a += dset->Z[list[1]][bigt];
                    for (j=1; j<pmod->ncoeff; j++) {
                        x = dset->Z[list[j+2]];
                        a -= pmod->coeff[j] * x[bigt];
                    }
                }
            }
            a /= Ti;
            ahat[k] = a;
            amean += a;
            if (wi != NULL) {
                wi[k] = Ti / (double) pmod->nobs;
                wmean += wi[k] * a;
            }
            k++;
        }
    }

    pan->s2v = 0.0;

    if (wi != NULL) {
        for (i=0; i<pan->effn; i++) {
            pan->s2v += wi[i] * (ahat[i] - wmean) * (ahat[i] - wmean);
        }
        pan->s2v /= (pan->effn - 1.0) / (double) pan->effn;
        free(wi);
    } else {
        amean /= pan->effn;
        for (i=0; i<pan->effn; i++) {
            pan->s2v += (ahat[i] - amean) * (ahat[i] - amean);
        }
        pan->s2v /= pan->effn - 1;
    }

    free(ahat);
    free(list);

    return err;
}

/* robust RE: we need an extra residuals array */

static int re_hatvars_prep (panelmod_t *pan)
{
    int t, n = pan->pooled->full_n;

    pan->re_uhat = malloc(n * sizeof *pan->re_uhat);

    if (pan->re_uhat == NULL) {
        return E_ALLOC;
    } else {
        for (t=0; t<n; t++) {
            pan->re_uhat[t] = NADBL;
        }
        return 0;
    }
}

/* When calculating DW using the fixed-effects residuals,
   in the context of estimation of the random-effects
   specification, we need to expand up the residual series
   temporarily.
*/

static int expand_fe_uhat (MODEL *femod, panelmod_t *pan)
{
    double *uhat = NULL;
    int NT = pan->pooled->full_n;
    int i, s, t;

    uhat = malloc(NT * sizeof *uhat);
    if (uhat == NULL) {
        return E_ALLOC;
    }

    for (t=0; t<NT; t++) {
        uhat[t] = NADBL;
    }

    s = 0;
    for (i=0; i<pan->nunits; i++) {
        int ti, Ti = pan->unit_obs[i];

        if (Ti == 0) {
            continue;
        }
        for (ti=0; ti<Ti; ti++) {
            t = big_index(pan, s);
            uhat[t] = femod->uhat[s];
            if (s == 0) {
                femod->t1 = t;
            } else if (s == femod->nobs - 1) {
                femod->t2 = t;
            }
            s++;
        }
    }

    /* replace the uhat array on @femod */
    free(femod->uhat);
    femod->uhat = uhat;

    return 0;
}

/* Fix uhat and yhat in two cases.

   (a) When we estimated fixed effects using a de-meaned dataset we
   need to ensure that the uhat and yhat values get written to the
   right observation slots in relation to the full dataset, and also
   that the yhat values get corrected, putting the means back in.

   (b) When estimating the random effects model we need to compute
   residuals based on the untransformed data (and again, place them
   correctly in relation to the full dataset). However, if we're
   going to produce robust standard errors we also need to preserve
   the GLS residuals.

   The placement issue arises because the special datasets used in
   these cases are not necessarily of full length, since they are
   purged of missing observations.
*/

static int
fix_panel_hatvars (MODEL *pmod, panelmod_t *pan, const double **Z)
{
    const double *y = NULL;
    double *yhat = pan->pooled->yhat;
    double *uhat = NULL;
    int n = pan->pooled->full_n;
    double yht, SSR = 0.0;
    int i, j, s, t;
    int err = 0;

    if (yhat == NULL) {
        fprintf(stderr, "fix_panel_hatvars: pan->pooled->yhat is NULL\n");
        return E_DATA;
    }

    y = Z[pan->pooled->list[1]];

    uhat = malloc(n * sizeof *uhat);
    if (uhat == NULL) {
        return E_ALLOC;
    }

    if (pan->opt & OPT_U) {
        /* random effects */
        if (pan->opt & OPT_R) {
            err = re_hatvars_prep(pan);
        }
        if (err) {
            return err;
        }
    }

    for (t=0; t<n; t++) {
        uhat[t] = NADBL;
    }

    s = 0;
    for (i=0; i<pan->nunits; i++) {
        int ti, Ti = pan->unit_obs[i];

        if (Ti == 0) {
            continue;
        }

        for (ti=0; ti<Ti; ti++) {
            t = big_index(pan, s);
            if (pan->opt & OPT_U) {
                /* random effects */
                yht = 0.0;
                for (j=0; j<pmod->ncoeff; j++) {
                    yht += pmod->coeff[j] * Z[pan->pooled->list[j+2]][t];
                }
                yhat[t] = yht;
                uhat[t] = y[t] - yht;
                SSR += uhat[t] * uhat[t];
                if (pan->re_uhat != NULL) {
                    /* store both the "fixed" and the GLS residuals */
                    pan->re_uhat[t] = uhat[t];
                    uhat[t] = pmod->uhat[s];
                }
            } else {
                /* fixed effects */
                uhat[t] = pmod->uhat[s];
                yhat[t] = y[t] - uhat[t];
            }
            if (s == 0) {
                pmod->t1 = t;
            } else if (s == pmod->nobs - 1) {
                pmod->t2 = t;
            }
            s++;
        }
    }

    if (pan->opt & OPT_U) {
        double r;

        /* Defer rewriting of pmod->ess and pmod->sigma to avoid
           screwing up the covariance matrix estimator and Hausman
           test.
        */
        gretl_model_set_double(pmod, "fixed_SSR", SSR);
        r = gretl_corr(pmod->t1, pmod->t2, y, yhat, NULL);
        if (!na(r)) {
            gretl_model_set_double(pmod, "corr-rsq", r * r);
        }
    }

    pmod->full_n = n;

    /* replace the uhat and yhat arrays on @pmod */
    free(pmod->uhat);
    pmod->uhat = uhat;
    free(pmod->yhat);
    pmod->yhat = yhat;

    /* NULLify stolen pointer */
    pan->pooled->yhat = NULL;

    return err;
}

static int
hausman_move_uhat (MODEL *pmod, panelmod_t *pan)
{
    int n = pan->pooled->full_n;
    double *uhat;
    int i, s, t, ti;

    uhat = malloc(n * sizeof *uhat);
    if (uhat == NULL) {
        return E_ALLOC;
    }

    for (t=0; t<n; t++) {
        uhat[t] = NADBL;
    }

    s = 0;
    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];

        if (Ti > 0) {
            for (ti=0; ti<Ti; ti++) {
                t = big_index(pan, s);
                uhat[t] = pmod->uhat[s];
                s++;
            }
        }
    }

    free(pmod->uhat);
    pmod->uhat = uhat;

    return 0;
}

#if PDEBUG > 1

static void verbose_femod_print (MODEL *femod, DATASET *wset,
                                 PRN *prn)
{
    pprintf(prn, "*** initial FE model (on within data)\n");
    printmodel(femod, wset, OPT_O, prn);

# if PDEBUG > 2
    int i, j;

    fprintf(stderr, "femod: data series length = %d\n", wset->n);
    for (i=0; i<wset->n; i++) {
        fprintf(stderr, "femod.uhat[%d] = %g, ", i, femod->uhat[i]);
        fprintf(stderr, "data: ");
        for (j=0; j<wset->v; j++) {
            fprintf(stderr, "%g ", wset->Z[j][i]);
        }
        fputc('\n', stderr);
    }
# endif
}

#endif

/* Estimate the fixed-effects model using a parallel dataset
   with the group means subtracted from all variables.
*/

static MODEL
fixed_effects_model (panelmod_t *pan, DATASET *dset, PRN *prn)
{
    MODEL femod;
    gretlopt lsqopt = OPT_A | OPT_X;
    DATASET *wset = NULL;
    int *felist = NULL;
    int save_qr = 0;
    int i;

#if PDEBUG
    fprintf(stderr, "fixed_effects: using de-meaned data\n");
#endif

    gretl_model_init(&femod, dset);

    felist = gretl_list_new(pan->vlist[0]);
    if (felist == NULL) {
        femod.errcode = E_ALLOC;
        return femod;
    }

    wset = within_groups_dataset(dset, pan);
    if (wset == NULL) {
        free(felist);
        femod.errcode = E_ALLOC;
        return femod;
    }

    felist[1] = 1;
    felist[2] = 0;
    for (i=3; i<=felist[0]; i++) {
        felist[i] = i - 1;
    }

    if (pan->opt & OPT_F) {
        /* suppress auto-removal of collinear terms */
        lsqopt |= OPT_Z;
    }
    if (pan->opt & OPT_N) {
        /* suppress df correction */
        lsqopt |= OPT_N;
    }

    save_qr = libset_get_int(USE_QR);
    libset_set_int(USE_QR, 1);

    femod = lsq(felist, wset, OLS, lsqopt);

    libset_set_int(USE_QR, save_qr);

    if (femod.errcode) {
        fprintf(stderr, "femod.errcode = %d\n", femod.errcode);
    } else if ((pan->opt & OPT_F) && femod.list[0] < felist[0]) {
        femod.errcode = E_SINGULAR;
    } else {
        if (!(pan->opt & OPT_N)) {
            /* we estimated a bunch of group means, and have to
               subtract degrees of freedom */
            fixed_effects_df_correction(&femod, pan->effn - 1);
        }
#if PDEBUG > 1
        verbose_femod_print(&femod, wset, prn);
#endif
        if (pan->opt & OPT_F) {
            /* estimating the FE model in its own right */
            if (pan->opt & OPT_R) {
                /* we have to do this before the pooled residual
                   array is "stolen" for the fixed-effects model
                */
                robust_fixed_effects_F(pan);
            }
            fix_panel_hatvars(&femod, pan, (const double **) dset->Z);
            if (pan->opt & OPT_R) {
                femod.errcode = panel_robust_vcv(&femod, pan, wset, dset);
            } else {
                femod_regular_vcv(&femod);
            }
        } else if (pan->opt & OPT_E) {
            if (IGLS) {
                read_true_variances(pan);
            } else {
                femod.errcode = nerlove_s2v(&femod, dset, pan);
            }
        }
    }

    destroy_dataset(wset);
    free(felist);

    return femod;
}

/* Construct a gretl list containing the index numbers of the
   cross-sectional units included in the fixed-effects
   regression. This is for the purpose of naming the per-unit
   intercepts.
*/

static int *fe_units_list (const panelmod_t *pan)
{
    int *ulist = NULL;
    int i, j, n = 0;

    for (i=0; i<pan->nunits; i++) {
        if (pan->unit_obs[i] > 0) {
            n++;
        }
    }

    ulist = gretl_list_new(n);

    if (ulist != NULL) {
        j = 1;
        for (i=0; i<pan->nunits; i++) {
            if (pan->unit_obs[i] > 0) {
                ulist[j++] = i + 1;
            }
        }
    }

    return ulist;
}

/* Compose a list referencing all variables that were dropped from the
   final panel model relative to the incoming regression
   specification.  This may include some variables that were dropped
   at the stage of running the baseline pooled model (presumably
   because of perfect collinearity).

   In the case of fixed effects it may include additional variables
   dropped due to the fact that they are time-invariant.
*/

static int compose_panel_droplist (MODEL *pmod, panelmod_t *pan)
{
    int fixed_effects = (pmod->opt & OPT_F);
    const int *pooldrop;
    int *dlist;
    int i, ndrop = 0;

    /* regressors dropped at the stage of estimating the
       initial pooled model */
    pooldrop = gretl_model_get_list(pan->pooled, "droplist");
    if (pooldrop != NULL) {
        ndrop += pooldrop[0];
    }

    if (fixed_effects) {
        /* regressors dropped because time-invariant */
        ndrop += pan->pooled->list[0] - pan->vlist[0];
    }

    if (ndrop == 0) {
        /* nothing to be done */
        return 0;
    }

    dlist = gretl_list_new(ndrop);
    if (dlist == NULL) {
        return E_ALLOC;
    }

    i = 1;

    if (pooldrop != NULL) {
        for (i=1; i<=pooldrop[0]; i++) {
            dlist[i] = pooldrop[i];
        }
    }

    if (fixed_effects) {
        if (pan->vlist[0] < pan->pooled->list[0]) {
            int j, vj;

            for (j=2; j<=pan->pooled->list[0]; j++) {
                vj = pan->pooled->list[j];
                if (!in_gretl_list(pan->vlist, vj)) {
                    dlist[i++] = vj;
                }
            }
        }
    }

    return gretl_model_set_list_as_data(pmod, "droplist", dlist);
}

static void fix_gls_stats (MODEL *pmod, panelmod_t *pan)
{
    double chisq = NADBL;
    int nc, df = 0;

    fix_panelmod_list(pmod, pan);

    pmod->ybar = pan->pooled->ybar;
    pmod->sdy = pan->pooled->sdy;
    pmod->tss = pan->pooled->tss;

    pmod->rsq = NADBL;
    pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;
    pmod->rho = pan->rho;
    pmod->dw = pan->dw;

    /* add joint chi-square test on regressors */

    if (pan->ntdum > 0) {
        df = pmod->ncoeff - 1 - pan->ntdum;
        if (df > 0) {
            chisq = panel_overall_test(pmod, pan, pan->ntdum,
                                       OPT_X);
        }
    } else {
        df = pmod->ncoeff - 1;
        if (df > 0) {
            chisq = panel_overall_test(pmod, pan, 0, OPT_X);
        }
    }

    /* deferred rewriting of ess and sigma, etc. */
    pmod->ess = gretl_model_get_double(pmod, "fixed_SSR");
    pmod->sigma = sqrt(pmod->ess / (pmod->nobs - (pmod->ncoeff - 1)));
    gretl_model_destroy_data_item(pmod, "fixed_SSR");
    nc = pmod->ncoeff;
    pmod->ncoeff = pmod->dfn + 1;
    ls_criteria(pmod);
    pmod->ncoeff = nc;

    if (!na(chisq) && chisq >= 0.0) {
        ModelTest *test = model_test_new(GRETL_TEST_RE_WALD);

        if (test != NULL) {
            model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
            model_test_set_dfn(test, df);
            model_test_set_value(test, chisq);
            model_test_set_pvalue(test, chisq_cdf_comp(df, chisq));
            maybe_add_test_to_model(pmod, test);
        }
    }
}

static void add_panel_obs_info (MODEL *pmod, panelmod_t *pan)
{
    gretl_model_set_int(pmod, "n_included_units", pan->effn);
    gretl_model_set_int(pmod, "panel_T", pan->T);
    gretl_model_set_int(pmod, "Tmin", pan->Tmin);
    gretl_model_set_int(pmod, "Tmax", pan->Tmax);
    if (pan->Tmax > pan->Tmin && pan->Tbar == 0) {
        calculate_Tbar(pan);
    }
    if (pan->Tbar > 0) {
        gretl_model_set_double(pmod, "Tbar", pan->Tbar);
    }
}

static void replace_re_residuals (MODEL *pmod, panelmod_t *pan)
{
    int t;

    /* replace the GLS residuals with the "fixed up" ones */

    for (t=0; t<pmod->full_n; t++) {
        pmod->uhat[t] = pan->re_uhat[t];
    }
}

/* We use this to "finalize" models estimated via fixed effects
   and random effects */

static int save_panel_model (MODEL *pmod, panelmod_t *pan,
                             const double **Z,
                             const DATASET *dset)
{
    int err = 0;

    pmod->ci = PANEL;

    if (pan->opt & OPT_F) {
        /* fixed effects */
        int *ulist;

        pmod->opt |= OPT_F;
        err = fix_within_stats(pmod, pan);
        if (err) {
            return err;
        }
        ulist = fe_units_list(pan);
        gretl_model_add_panel_varnames(pmod, dset, ulist);
        free(ulist);
        panel_model_add_ahat(pmod, dset, pan);
        save_fixed_effects_F(pan, pmod);
    } else {
        /* random effects */
        pmod->opt |= OPT_U;
        gretl_model_set_double(pmod, "s2v", pan->s2v);
        gretl_model_set_double(pmod, "s2e", pan->s2e);
        if (pan->balanced) {
            gretl_model_set_double(pmod, "theta", pan->theta);
        } else if (!na(pan->theta_bar)) {
            pan->theta_bar /= pan->effn;
            gretl_model_set_double(pmod, "theta_bar", pan->theta_bar);
        }
        gretl_model_add_panel_varnames(pmod, dset, NULL);
        if (pan->re_uhat != NULL) {
            replace_re_residuals(pmod, pan);
        }
        fix_gls_stats(pmod, pan);
        panel_model_add_ahat(pmod, dset, pan);
        if (pan->opt & OPT_E) {
            /* record use of Nerlove transformation */
            pmod->opt |= OPT_E;
        }
        if ((pan->opt & OPT_X) && !IGLS) {
            /* record use of special unbalanced ANOVA */
            pmod->opt |= OPT_X;
            if (!(pan->opt & OPT_E)) {
                const char *meth = stata_sa ? "stata" : "bc";

                gretl_model_set_string_as_data(pmod, "anova_method",
                                               gretl_strdup(meth));
            }
        }
    }

    /* compose list of dropped variables, if any */
    compose_panel_droplist(pmod, pan);

    if (!(pan->opt & OPT_A)) {
        set_model_id(pmod, pan->opt);
    }

    if (pan->opt & OPT_F) {
        panel_dwstat(pmod, pan);
    }

    if (pan->opt & OPT_R) {
        pmod->opt |= OPT_R;
    }

    time_dummies_wald_test(pan, pmod);
    add_panel_obs_info(pmod, pan);
    *pan->realmod = *pmod;

    return err;
}

/* drive the calculation of the fixed effects regression, print the
   results (if wanted), and compute the "within" error variance */

static int within_variance (panelmod_t *pan,
                            DATASET *dset,
                            PRN *prn)
{
    MODEL femod;
    int i, err = 0;

    femod = fixed_effects_model(pan, dset, prn);

    if (femod.errcode) {
        pputs(prn, _("Error estimating fixed effects model\n"));
        errmsg(femod.errcode, prn);
        err = femod.errcode;
        clear_model(&femod);
    } else {
        int den;

        if (pan->opt & OPT_E) {
            /* Nerlove */
            den = femod.nobs;
        } else {
            /* as per Greene: nT - n - K */
            den = femod.nobs - pan->effn - (pan->vlist[0] - 2);
        }

        if (den == 0) {
            gretl_errmsg_set(_("Inadequate data for panel estimation"));
            return E_DF;
        } else {
            pan->s2e = femod.ess / den;
        }

#if PDEBUG
        fprintf(stderr, "nT = %d, n = %d, K = %d\n", femod.nobs,
                pan->effn, pan->vlist[0] - 2);
        fprintf(stderr, "pan->s2e = %g / %d = %g\n", femod.ess,
                den, pan->s2e);
        fprintf(stderr, "sqrt(pan->s2e) = %g\n", sqrt(pan->s2e));
#endif

        if (!(pan->opt & OPT_R)) {
            regular_fixed_effects_F(pan, &femod);
        }
        if (pan->opt & OPT_V) {
            print_fe_results(pan, &femod, dset, prn);
        }
        if (pan->bdiff != NULL && pan->Sigma != NULL) {
            for (i=1; i<femod.ncoeff; i++) {
                pan->bdiff->val[i-1] = femod.coeff[i];
            }
            femod_regular_vcv(&femod);
            vcv_slopes(pan, &femod, VCV_INIT);
        }
        if (pan->opt & OPT_F) {
            err = save_panel_model(&femod, pan, (const double **) dset->Z,
                                   dset);
            if (err) {
                clear_model(&femod);
            }
        } else {
            if (pan->opt & OPT_U) {
                if (expand_fe_uhat(&femod, pan) == 0) {
                    panel_dwstat(&femod, pan);
                }
            }
            clear_model(&femod);
        }
    }

    return err;
}

static void print_hausman_result (panelmod_t *pan, PRN *prn)
{
    if (na(pan->H)) {
        pputs(prn, _("Hausman test matrix is not positive definite!\n"));
    } else {
        pprintf(prn, _("Hausman test statistic:\n"
                       " H = %g with p-value = prob(chi-square(%d) > %g) = %g\n"),
                pan->H, pan->dfH, pan->H, chisq_cdf_comp(pan->dfH, pan->H));
        pputs(prn, _("(A low p-value counts against the null hypothesis that "
                     "the random effects\nmodel is consistent, in favor of the fixed "
                     "effects model.)\n"));
    }
}

static void save_hausman_result (panelmod_t *pan)
{
    ModelTest *test;

    if (pan->realmod == NULL || na(pan->H) || pan->dfH == 0) {
        return;
    }

    test = model_test_new(GRETL_TEST_PANEL_HAUSMAN);

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
        model_test_set_dfn(test, pan->dfH);
        model_test_set_value(test, pan->H);
        if (na(pan->H)) {
            model_test_set_pvalue(test, NADBL);
        } else {
            model_test_set_pvalue(test, chisq_cdf_comp(pan->dfH, pan->H));
        }
        maybe_add_test_to_model(pan->realmod, test);
    }
}

/* Handle the case where collinear terms were dropped
   when doing the Hausman test via the regression
   method in robust mode.
*/

static double robust_hausman_fixup (const int *hlist,
                                    MODEL *pmod)
{
    int *wlist = gretl_list_copy(hlist);
    double H = NADBL;

    if (wlist != NULL) {
        int i;

        for (i=wlist[0]; i>0; i--) {
            if (!in_gretl_list(pmod->list, wlist[i])) {
                gretl_list_delete_at_pos(wlist, i);
            }
        }
        H = wald_omit_chisq(wlist, pmod);
        free(wlist);
    }

    return H;
}

/* Estimate the augmented GLS model for the Hausman test;
   return either the error sum of squares or, in the
   robust case, a Wald chi-square statistic.
*/

static double hausman_regression_result (panelmod_t *pan,
                                         const int *relist,
                                         const int *hlist,
                                         DATASET *rset,
                                         DATASET *dset,
                                         PRN *prn)
{
    double ret = NADBL;
    int *biglist = NULL;
    int err = 0;

    biglist = gretl_list_add(relist, hlist, &err);

    if (biglist != NULL) {
        MODEL hmod;

        gretl_model_init(&hmod, NULL);
        hmod = lsq(biglist, rset, OLS, OPT_A);
#if PDEBUG > 1
        pputs(prn, "Hausman test regression\n");
        printmodel(&hmod, rset, OPT_NONE, prn);
#endif
        if (hmod.errcode == 0) {
            /* Find the number of additional regressors actually used,
               relative to @relist, allowing for the possibility
               that one or more elements of @biglist were dropped
               in estimation of @hmod due to excessive collinearity.
            */
            pan->dfH = hmod.list[0] - relist[0];
            if (pan->dfH > 0) {
                if (pan->opt & OPT_R) {
                    /* do robust Wald test */
                    if (hmod.full_n < pan->pooled->full_n) {
                        hausman_move_uhat(&hmod, pan);
                    }
                    panel_robust_vcv(&hmod, pan, rset, dset);
                    if (hmod.vcv != NULL) {
                        if (pan->dfH < hlist[0]) {
                            ret = robust_hausman_fixup(hlist, &hmod);
                        } else {
                            ret = wald_omit_chisq(hlist, &hmod);
                        }
                    }
                } else {
                    /* just record the unrestricted SSR */
                    ret = hmod.ess;
                }
            }
        }
        clear_model(&hmod);
    }

    free(biglist);

    return ret;
}

/* Computation of s2_v in the manner of Swamy and Arora, for an
   unbalanced panel. See Baltagi, Econometric Analysis of Panel
   Data, 3e, section 9.2.1. */

static int unbalanced_SA_s2v (panelmod_t *pan,
                              DATASET *dset)
{
    gretl_matrix_block *B;
    gretl_matrix *ZmZ = NULL;
    gretl_matrix *PZ = NULL;
    gretl_matrix *Z = NULL;
    gretl_matrix *ZPZ = NULL;
    gretl_matrix *D2 = NULL;
    gretl_matrix *trmat = NULL;
    int k = pan->pooled->ncoeff;
    double zjt, z, tr;
    int i, j, s, t, p;
    int bigt, got;
    int err = 0;

    B = gretl_matrix_block_new(&ZmZ, pan->effn, k,
                               &PZ,  pan->NT, k,
                               &Z,   pan->NT, k,
                               &ZPZ, k, k,
                               &D2,  k, k, NULL);

    if (B == NULL) {
        return E_ALLOC;
    }

    s = p = 0;
    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];

        if (Ti == 0) {
            continue;
        }

        for (j=0; j<k; j++) {
            zjt = 0.0;
            got = 0;
            /* cumulate sum of observations for unit */
            for (t=0; t<pan->T && got<Ti; t++) {
                bigt = panel_index(i, t);
                if (!panel_missing(pan, bigt)) {
                    z = dset->Z[pan->pooled->list[j+2]][bigt];
                    zjt += z;
                    gretl_matrix_set(Z, p + got, j, z);
                    got++;
                }
            }
            gretl_matrix_set(ZmZ, s, j, zjt);
            for (t=0; t<Ti; t++) {
                gretl_matrix_set(PZ, p + t, j, zjt / Ti);
            }
        }
        s++;
        p += Ti;
    }

    gretl_matrix_multiply_mod(Z, GRETL_MOD_TRANSPOSE,
                              PZ, GRETL_MOD_NONE,
                              ZPZ, GRETL_MOD_NONE);
    gretl_invert_symmetric_matrix(ZPZ);

    gretl_matrix_multiply_mod(ZmZ, GRETL_MOD_TRANSPOSE,
                              ZmZ, GRETL_MOD_NONE,
                              D2, GRETL_MOD_NONE);

    trmat = gretl_matrix_reuse(PZ, k, k);
    gretl_matrix_multiply(ZPZ, D2, trmat);
    tr = gretl_matrix_trace(trmat);

#if 0
    gretl_matrix_print(trmat, "trmat");
    fprintf(stderr, "S-A: effn=%d, k=%d, NT=%d, tr=%g\n", pan->effn,
            k, pan->NT, tr);
#endif

    pan->s2v = (pan->ubPub - (pan->effn - k) * pan->s2e) / (pan->NT - tr);

#if PDEBUG
    fprintf(stderr, "S-A: ubPub=%#.8g, tr=%#.8g, s2v=%#.8g, sv=%#.8g\n",
            pan->ubPub, tr, pan->s2v, sqrt(pan->s2v));
#endif

    gretl_matrix_block_destroy(B);

    return err;
}

/* Calculate the random effects regression.  Print the results
   here if we're doing the "panel diagnostics" test, otherwise
   save the results.
*/

static int random_effects (panelmod_t *pan,
                           DATASET *dset,
                           DATASET *gset,
                           PRN *prn)
{
    DATASET *rset;
    MODEL remod;
    gretlopt lsqopt = OPT_A | OPT_Z;
    int *relist = NULL;
    int *hlist = NULL;
    double hres = NADBL;
    int i, err = 0;

    gretl_model_init(&remod, dset);

    /* FGLS regression list */
    relist = gretl_list_new(pan->pooled->list[0]);
    if (relist == NULL) {
        return E_ALLOC;
    }

    if (!(pan->opt & OPT_M)) {
        /* extra regressors for Hausman test, regression approach */
        hlist = gretl_list_new(pan->vlist[0] - 1);
        if (hlist == NULL) {
            free(relist);
            return E_ALLOC;
        }
    }

    /* If OPT_E (--nerlove) was given, we've already calculated
       pan->s2v as the variance of the fixed effects. Otherwise
       we're using the Swamy and Arora method, and pan->s2v still
       needs to be computed.

       Note: for unbalanced panels, theta will vary across the
       units in the final calculation.
    */
    if (!(pan->opt & OPT_E)) {
        if (IGLS) {
            /* get user-specified values */
            err = read_true_variances(pan);
        } else if (!pan->balanced && !na(pan->ubPub)) {
            err = unbalanced_SA_s2v(pan, dset);
        } else {
            pan->s2v = pan->s2b - pan->s2e / pan->Tbar;
#if PDEBUG
            fprintf(stderr, "Swamy-Arora: initial s2v = %g - %g = %g\n",
                    pan->s2b, pan->s2e / pan->Tbar, pan->s2v);
#endif
        }
        if (pan->s2v < 0) {
            pan->s2v = 0.0;
        }
    }

    /* theta, the quasi-demeaning coefficient */
    pan->theta = 1.0 - sqrt(pan->s2e / (pan->s2e + pan->Tmax * pan->s2v));

#if PDEBUG
    fprintf(stderr, "pan->s2v = %.8g, pan->s2e = %.8g\n", pan->s2v, pan->s2e);
    fprintf(stderr, "random_effects theta = %g\n", pan->theta);
#endif

    /* make special transformed dataset, and regression list */
    rset = random_effects_dataset(dset, gset, relist, hlist, pan);
    if (rset == NULL) {
        free(relist);
        free(hlist);
        return E_ALLOC;
    }

    if (hlist != NULL) {
        hres = hausman_regression_result(pan, relist, hlist,
                                         rset, dset, prn);
    }

    /* regular random-effects model */
#if PDEBUG
    fprintf(stderr, "estimate regular GLS model\n");
#endif
    remod = lsq(relist, rset, OLS, lsqopt);

    if ((err = remod.errcode)) {
        pputs(prn, _("Error estimating random effects model\n"));
        errmsg(err, prn);
    } else {
#if PDEBUG > 1
        for (i=0; i<20 && i<remod.nobs; i++) {
            fprintf(stderr, "remod uhat[%d] = %g\n", i, remod.uhat[i]);
        }
#endif
        if (pan->bdiff != NULL) {
            /* matrix-diff Hausman variant, if wanted */
            int vi, k = 0;

            for (i=0; i<remod.ncoeff; i++) {
                vi = pan->pooled->list[i+2];
                if (var_is_varying(pan, vi)) {
                    pan->bdiff->val[k++] -= remod.coeff[i];
                }
            }
        }

        /* note: this call moved to here 2017-10-18 so we
           get computation of robust standard errors right
        */
        fix_panel_hatvars(&remod, pan, (const double **) dset->Z);

        if (pan->opt & OPT_R) {
            err = panel_robust_vcv(&remod, pan, rset, dset);
            if (err) {
                goto bailout;
            }
        } else {
            double sigma = remod.sigma; /* or... ? */
#if PDEBUG
            fprintf(stderr, "GLS sigma = %g\n", sigma);
            fprintf(stderr, "GLS SSR   = %.7g\n", remod.ess);
#endif
            makevcv(&remod, sigma);
        }

        if (pan->opt & OPT_V) {
            print_re_results(pan, &remod, dset, prn);
        }

        if (!na(hres)) {
            if (pan->opt & OPT_R) {
                /* @hres is already the (robust) Hausman test statistic */
                pan->H = hres;
            } else {
                /* @hres is the unrestricted ess */
                pan->H = (remod.ess / hres - 1.0) * (remod.nobs);
            }
        } else if (pan->Sigma != NULL) {
            vcv_slopes(pan, &remod, VCV_SUBTRACT);
        }
    }

#if PDEBUG > 1
    if (remod.errcode == 0) {
        printmodel(&remod, rset, OPT_NONE, prn);
    }
#endif

 bailout:

    if (!err && (pan->opt & OPT_U)) {
        save_panel_model(&remod, pan, (const double **) dset->Z, rset);
    } else {
        clear_model(&remod);
    }

    destroy_dataset(rset);

    free(relist);
    free(hlist);

    return err;
}

static void save_breusch_pagan_result (panelmod_t *pan)
{
    ModelTest *test;

    if (pan->realmod == NULL || na(pan->BP)) {
        return;
    }

    test = model_test_new(GRETL_TEST_PANEL_BP);

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
        model_test_set_dfn(test, 1);
        model_test_set_value(test, pan->BP);
        model_test_set_pvalue(test, chisq_cdf_comp(1, pan->BP));
        maybe_add_test_to_model(pan->realmod, test);
    }
}

/* do the panel Breusch-Pagan test and either print or save
   the results */

static int breusch_pagan_LM (panelmod_t *pan, PRN *prn)
{
    double A = 0.0;
    int n = pan->pooled->nobs;
    int print_means = 0;
    int i, t, M = 0;

    if ((pan->opt & OPT_V) && pan->effn <= 10) {
        print_means = 1;
    }

    if (print_means) {
        pputs(prn, _("Means of pooled OLS residuals for cross-sectional units:"));
        pputs(prn, "\n\n");
    }

    for (i=0; i<pan->nunits; i++) {
        int Ti = pan->unit_obs[i];

        if (Ti > 0) {
            double u, usum = 0.0;

            for (t=0; t<pan->T; t++) {
                u = pan->pooled->uhat[panel_index(i, t)];
                if (!na(u)) {
                    usum += u;
                }
            }
            if (print_means) {
                pprintf(prn, _(" unit %2d: %13.5g\n"), i + 1,
                        usum / (double) Ti);
            }
            A += usum * usum;
            M += Ti * Ti;
        }
    }

    pan->BP = ((double) n * n /(2.0 * (M - n))) * pow((A / pan->pooled->ess) - 1.0, 2);

    if (pan->opt & OPT_V) {
        pputc(prn, '\n');
        pprintf(prn, _("Breusch-Pagan test statistic:\n"
                       " LM = %g with p-value = prob(chi-square(1) > %g) = %g\n"),
                pan->BP, pan->BP, chisq_cdf_comp(1, pan->BP));

        pputs(prn, _("(A low p-value counts against the null hypothesis that "
                     "the pooled OLS model\nis adequate, in favor of the random "
                     "effects alternative.)\n\n"));
    }

    return 0;
}

static void save_panspec_result (panelmod_t *pan, PRN *prn)
{
    gretl_matrix *tests = gretl_column_vector_alloc(3);
    gretl_matrix *pvals = gretl_column_vector_alloc(3);
    char **S = strings_array_new(3);

    /* fixed-effects poolability F-test */
    tests->val[0] = pan->Ffe;
    pvals->val[0] = snedecor_cdf_comp(pan->Fdfn, pan->Fdfd, pan->Ffe);
    /* random-effects poolability test */
    tests->val[1] = pan->BP;
    pvals->val[1] = chisq_cdf_comp(1, pan->BP);
    /* Hausman test */
    tests->val[2] = pan->H;
    pvals->val[2] = chisq_cdf_comp(pan->dfH, pan->H);

    if (S != NULL) {
        /* add row names */
        char **S2;

        S[0] = gretl_strdup("Poolability (Wald)");
        S[1] = gretl_strdup("Poolability (B-P)");
        S[2] = gretl_strdup("Hausman");
        S2 = strings_array_dup(S, 3);
        gretl_matrix_set_rownames(tests, S);
        if (S2 != NULL) {
            gretl_matrix_set_rownames(pvals, S2);
        }
    }

    if (pan->opt & OPT_V) {
        print_hausman_result(pan, prn);
    }
    record_matrix_test_result(tests, pvals);
}

static int finalize_hausman_test (panelmod_t *pan, int ci, PRN *prn)
{
    int mdiff = 0, err = 0;

    if (pan->bdiff != NULL && pan->Sigma != NULL) {
        /* matrix approach */
        mdiff = 1;
        err = bXb(pan);
    } else if (na(pan->H)) {
        /* regression approach bombed somehow? */
        err = E_DATA;
    }

    if (ci == PANSPEC) {
        /* the context is the "panspec" command */
        if (!err || (mdiff && err == E_NOTPD)) {
            save_panspec_result(pan, prn);
        }
    } else {
        /* the context is random-effects estimation */
        if (mdiff && err == E_NOTPD) {
            pputs(prn, _("Hausman test matrix is not positive definite"));
            pputc(prn, '\n');
        }
        save_hausman_result(pan);
    }

    return err;
}

/* Based on the residuals from pooled OLS, do some accounting to see
   (a) how many cross-sectional units were actually included (after
   omitting any missing values); (b) how many time-series observations
   were included for each unit; and (c) what was the maximum number
   of time-series observations used.  Return 0 if all goes OK,
   non-zero otherwise.
*/

static int panel_obs_accounts (panelmod_t *pan)
{
    int *uobs;
    int i, t, bigt;

    uobs = malloc(pan->nunits * sizeof *uobs);
    if (uobs == NULL) {
        return E_ALLOC;
    }

    pan->NT = 0;        /* to hold total complete observations */
    pan->effn = 0;      /* to hold number of included units */
    pan->Tmax = 0;      /* to hold max time-series length */
    pan->Tmin = pan->T; /* to hold min time-series length > 0 */

    for (i=0; i<pan->nunits; i++) {
        uobs[i] = 0;
        for (t=0; t<pan->T; t++) {
            bigt = panel_index(i, t);
#if PDEBUG > 2
            fprintf(stderr, "unit %d, bigt=%d, pmod->uhat[%d]: %s\n", i, t,
                    bigt, (panel_missing(pan, bigt))? "NA" : "OK");
#endif
            if (!panel_missing(pan, bigt)) {
                uobs[i] += 1;
            }
        }
        if (uobs[i] > 0) {
            pan->effn += 1;
            if (uobs[i] > pan->Tmax) {
                pan->Tmax = uobs[i];
            }
            if (uobs[i] < pan->Tmin) {
                pan->Tmin = uobs[i];
            }
            pan->NT += uobs[i];
        }
    }

    for (i=0; i<pan->nunits; i++) {
        if (uobs[i] > 0 && uobs[i] != pan->Tmax) {
            pan->balanced = 0;
            break;
        }
    }

    pan->unit_obs = uobs;

    return 0;
}

/* Construct an array of char to record which parameters in the full
   model correspond to time-varying regressors
*/

static int panel_set_varying (panelmod_t *pan, const MODEL *pmod)
{
    int i;

    pan->varying = calloc(pmod->ncoeff, 1);
    if (pan->varying == NULL) {
        return E_ALLOC;
    }

    for (i=0; i<pmod->ncoeff; i++) {
        if (var_is_varying(pan, pmod->list[i+2])) {
            pan->varying[i] = 1;
        }
    }

    return 0;
}

/* find harmonic mean of the number of time-series observations
   per included group */

static void calculate_Tbar (panelmod_t *pan)
{
    int i;

    if (pan->balanced) {
        pan->Tbar = pan->Tmax;
    } else {
        double den = 0.0;

        for (i=0; i<pan->nunits; i++) {
            if (pan->unit_obs[i] > 0) {
                den += 1.0 / pan->unit_obs[i];
            }
        }
        pan->Tbar = pan->effn / den;
    }
}

static int
panelmod_setup (panelmod_t *pan, MODEL *pmod, const DATASET *dset,
                int ntdum, gretlopt opt)
{
    int nobs, dsetT = sample_size(dset);
    int err = 0;

    pan->opt = opt;
    pan->pooled = pmod;

    /* assumes (possibly padded) balanced panel dataset */
    pan->nunits = dsetT / dset->pd;
    pan->T = dset->pd;
    nobs = pan->nunits * pan->T;

    if (nobs < dsetT) {
        fprintf(stderr, "dataset nobs = %d, but panel nobs = %d\n",
                dsetT, nobs);
    }

    panel_index_init(dset, pan->nunits, pan->T);
    pan->ntdum = ntdum;

    err = panel_obs_accounts(pan);

    if (!err && (pan->opt & (OPT_U | OPT_F | OPT_B))) {
        pan->realmod = malloc(sizeof *pan->realmod);
        if (pan->realmod == NULL) {
            err = E_ALLOC;
        }
    }

    IGLS = stata_sa = 0;

    if (!err && (opt & OPT_X)) {
        if (!(pan->opt & OPT_U)) {
            /* --unbalanced --random-effects */
            err = E_BADOPT;
        } else {
            /* probe optional argument to --unbalanced */
            const char *s = get_optval_string(PANEL, OPT_X);

            if (s != NULL) {
                if (!strcmp(s, "stata")) {
                    stata_sa = 1;
                } else if (!strcmp(s, "bc")) {
                    stata_sa = 0;
                } else if (has_suffix(s, ".mat")) {
                    /* hidden testing option */
                    strcpy(glsmat, s);
                    IGLS = 1;
                } else {
                    gretl_errmsg_sprintf(_("%s: invalid option argument"), s);
                    err = E_INVARG;
                }
            }
        }
    }

    if (err && pan->unit_obs != NULL) {
        free(pan->unit_obs);
        pan->unit_obs = NULL;
    }

    return err;
}

/* Called in relation to a model estimated by pooled OLS: test for
   both fixed and random effects. Implements the "panspec" command.
*/

int panel_diagnostics (MODEL *pmod, DATASET *dset,
                       gretlopt opt, PRN *prn)
{
    int nerlove = 0;
    int quiet = (opt & OPT_Q);
    gretlopt psopt = opt;
    panelmod_t pan;
    int xdf, err = 0;

#if PDEBUG
    fputs("\n*** Starting panel_diagnostics ***\n", stderr);
    fprintf(stderr, " OPT_M? %s\n", (opt & OPT_M)? "yes" : "no");
#endif

    if (pmod->ifc == 0) {
        /* at many points we assume the base regression has an
           intercept included */
        return E_NOCONST;
    }

#if PDEBUG /* not really ready yet! */
    if (pmod->opt & OPT_R) {
        /* forward the robust option */
        opt |= OPT_R;
    }
#endif

    if (opt & OPT_E) {
        nerlove = 1;
    }

    /* Add OPT_V to make the fixed and random effects functions verbose,
       unless we've been passed the --quiet option.
    */
    if (!quiet) {
        psopt |= OPT_V;
    }
    panelmod_init(&pan);
    err = panelmod_setup(&pan, pmod, dset, 0, psopt);
    if (err) {
        goto bailout;
    }

    if (pan.effn < pan.nunits) {
        fprintf(stderr, "number of units included = %d\n", pan.effn);
        if (pan.effn <= 0) {
            return E_DATA;
        }
    }

    /* figure out which of the original regressors are time-varying */
    err = varying_vars_list(dset, &pan);
    if (err) {
        goto bailout;
    }

    err = panel_set_varying(&pan, pmod);
    if (err) {
        goto bailout;
    }

    calculate_Tbar(&pan);

    /* degrees of freedom relative to # of x-sectional units */
    xdf = pan.effn - pmod->ncoeff;

#if PDEBUG
    fprintf(stderr, "nunits=%d, T=%d, effn=%d, Tmax=%d, xdf=%d\n",
            pan.nunits, pan.T, pan.effn, pan.Tmax, xdf);
#endif

    if (xdf <= 0 && !nerlove) {
        ; /* can't do Random Effects */
    } else if (pan.opt & OPT_M) {
        /* matrix version of Hausman test */
        err = matrix_hausman_allocate(&pan);
        if (err) {
            goto bailout;
        }
    }

    if (!quiet) {
        /* header */
        pputc(prn, '\n');
        pprintf(prn, _("Diagnostics: using n = %d cross-sectional units\n"),
                pan.effn);
        pputc(prn, '\n');
    }

    err = within_variance(&pan, dset, prn);
    if (err) {
        goto bailout;
    }

    if (xdf <= 0) {
        pprintf(prn, "Omitting group means regression: "
                "insufficient degrees of freedom\n");
        goto bailout;
    }

    if (xdf > 0 && !na(pan.s2e)) {
        DATASET *gset;

        gset = group_means_dataset(&pan, dset);
        if (gset == NULL) {
            err = E_ALLOC;
        } else {
            err = between_variance(&pan, gset);
        }

        if (err) {
            pputs(prn, _("Couldn't estimate group means regression\n"));
            if (err == E_SINGULAR) {
                /* treat this as non-fatal */
                err = 0;
            }
        } else {
            random_effects(&pan, dset, gset, prn);
            if (pan.theta > 0) {
                /* this test hardly makes sense if GLS = OLS */
                breusch_pagan_LM(&pan, prn);
            }
            finalize_hausman_test(&pan, PANSPEC, prn);
        }

        if (gset != NULL) {
            destroy_dataset(gset);
        }
    }

 bailout:

    panelmod_free(&pan);

    return err;
}

static int between_model (panelmod_t *pan, const DATASET *dset)
{
    DATASET *gset;
    int err = 0;

    gset = group_means_dataset(pan, dset);

    if (gset == NULL) {
        err = E_ALLOC;
    } else {
        err = between_variance(pan, gset);
        if (err) {
            destroy_dataset(gset);
        } else {
            pan->realmod->dataset = gset;
        }
    }

    return err;
}

/* When we automatically add dummy variables to the dataset,
   in estimating a panel model with the --time-dummies option,
   should we delete these when we're done, or leave them
   in the dataset? 2013-11-18: we'll delete them if the dataset
   is subsampled, otherwise leave them.
*/

static int
process_time_dummies (MODEL *pmod, const DATASET *dset, int v)
{
    int i, vi, n = 0;

    for (i=pmod->list[0]; i>1; i--) {
        if (pmod->list[i] >= v) {
            n++;
        }
    }

    if (n > 0) {
        gretl_model_allocate_param_names(pmod, pmod->ncoeff);

        if (pmod->errcode) {
            return pmod->errcode;
        }

        for (i=2; i<=pmod->list[0]; i++) {
            vi = pmod->list[i];
            gretl_model_set_param_name(pmod, i-2, dset->varname[vi]);
        }

        for (i=pmod->list[0]; i>1; i--) {
            if (pmod->list[i] >= v) {
                gretl_list_delete_at_pos(pmod->list, i);
            }
        }
    }

    /* Record the fact that the model used time dummies, in case
       the user wants to add/omit variables or something.
    */
    pmod->opt |= OPT_D;

    return 0;
}

static int add_time_dummies_to_list (const int *list,
                                     DATASET *dset,
                                     int **plist)
{
    char dname[VNAMELEN];
    int *biglist = NULL;
    int ntdum = dset->pd - 1;
    int i, j, v;
    int err = 0;

    biglist = gretl_list_new(list[0] + ntdum);
    if (biglist == NULL) {
        return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
        biglist[i] = list[i];
    }

    j = list[0] + 1;

    for (i=2; i<=dset->pd; i++) {
        sprintf(dname, "dt_%d", i);
        v = series_greatest_index(dset, dname);
        if (v == dset->v) {
            err = E_DATA;
            break;
        } else if (in_gretl_list(list, v)) {
            /* already present */
            biglist[0] -= 1;
        } else {
            biglist[j++] = v;
        }
    }

    if (err) {
        free(biglist);
    } else {
        *plist = biglist;
    }

    return err;
}

static int panel_check_for_const (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
        if (list[i] == 0) {
            return 0; /* OK */
        }
    }

    return E_NOCONST;
}

static int get_ntdum (const int *orig, const int *new)
{
    int i, n = 0;

    for (i=2; i<=new[0]; i++) {
        if (!in_gretl_list(orig, new[i])) {
            n++;
        }
    }

    return n;
}

static int save_pooled_model (MODEL *pmod, panelmod_t *pan,
			      DATASET *dset)
{
    int err = 0;

    if (pan->opt & OPT_R) {
        err = panel_robust_vcv(pmod, pan, dset, dset);
	if (!err) {
	    pmod->opt |= OPT_R;
	    pmod->fstt = panel_overall_test(pmod, pan, 0, OPT_NONE);
	    pmod->dfd = pan->effn - 1;
	}
    }

    if (!err) {
	gretl_model_set_int(pmod, "pooled", 1);
	add_panel_obs_info(pmod, pan);
	if (!(pan->opt & OPT_A)) {
	    set_model_id(pmod, pan->opt);
	}
	panel_dwstat(pmod, pan);
    }

    return err;
}

/* fixed effects | random effects | between | pooled */

#define estimator_specified(o) (o & (OPT_F|OPT_U|OPT_B|OPT_P))

/**
 * real_panel_model:
 * @list: list containing model specification.
 * @dset: dataset struct.
 * @opt: may include %OPT_U for the random effects model;
 * %OPT_R for robust standard errors; %OPT_M to use the
 * matrix-difference variant of the Hausman test (random
 * effects only); %OPT_B for the "between" model; %OPT_P for
 * pooled OLS; and %OPT_D to include time dummies.
 * %OPT_C for clustered standard errors is also accepted.
 * If %OPT_U is given, either of the mutually incompatible options
 * %OPT_E and %OPT_X may be given to inflect the calculation of the
 * variance of the individual effects: %OPT_E means use Nerlove's
 * method, %OPT_X means use the special "exact" method of Swamy
 * and Arora in the case of an unbalanced panel. The default is
 * to use the "generic" Swamy-Arora method.
 * @prn: printing struct.
 *
 * Estimates a panel model, by default using the fixed effects
 * estimator.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL real_panel_model (const int *list, DATASET *dset,
                        gretlopt opt, PRN *prn)
{
    int orig_v = dset->v;
    MODEL mod;
    panelmod_t pan;
    gretlopt pan_opt = opt;
    gretlopt ols_opt = OPT_A;
    int *olslist = NULL;
    int nerlove = 0;
    int ntdum = 0;
    int err = 0;

    gretl_model_init(&mod, dset);
    panelmod_init(&pan);

    if (!estimator_specified(opt)) {
        /* default: add OPT_F to save the fixed effects model */
        pan_opt |= OPT_F;
    }

    if (opt & OPT_P) {
        /* doing pooled OLS */
        if (opt & OPT_N) {
            /* no-df-corr */
            ols_opt |= OPT_N;
        }
    } else {
        /* all other variants */
        mod.errcode = panel_check_for_const(list);
        if (mod.errcode) {
            return mod;
        }
    }

    if (opt & OPT_D) {
        /* time dummies wanted */
        err = gen_panel_dummies(dset, OPT_T, prn);
        if (!err) {
            err = add_time_dummies_to_list(list, dset, &olslist);
        }
    } else {
        olslist = gretl_list_copy(list);
        if (olslist == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        goto bailout;
    }

    /* baseline: estimate via pooled OLS */
    mod = lsq(olslist, dset, OLS, ols_opt);
    if (mod.errcode) {
        err = mod.errcode;
        fprintf(stderr, "real_panel_model: error %d in initial OLS\n",
                mod.errcode);
        goto bailout;
    }

#if PDEBUG > 1
    pprintf(prn, "*** initial baseline OLS\n");
    printlist(olslist, "olslist");
    printmodel(&mod, dset, OPT_NONE, prn);
#endif

    free(olslist);

#if PDEBUG
    if (pan_opt & OPT_F) {
        fprintf(stderr, "\n*** Doing fixed effects\n");
    } else if (pan_opt & OPT_U) {
        fprintf(stderr, "\n*** Doing random effects\n");
    } else if (pan_opt & OPT_P) {
	fprintf(stderr, "\n*** Doing pooled OLS\n");
    }
#endif

    if (pan_opt & OPT_C) {
        /* cluster implies robust */
        pan_opt |= OPT_R;
    }
    if (pan_opt & OPT_N) {
        /* no-df-corr (not working!)  */
        pan_opt |= OPT_N;
    }

    if ((opt & OPT_D) && !(opt & OPT_B)) {
        ntdum = get_ntdum(list, mod.list);
    }

    err = panelmod_setup(&pan, &mod, dset, ntdum, pan_opt);
    if (err) {
        goto bailout;
    }

    if (opt & OPT_P) {
        err = save_pooled_model(&mod, &pan, dset);
        goto bailout;
    }

    if (opt & OPT_B) {
        err = between_model(&pan, dset);
        goto bailout;
    }

    /* figure out which of the original regressors are time-varying,
       or unit-varying as the case may be
    */
    err = varying_vars_list(dset, &pan);
    if (err) {
        goto bailout;
    }

    err = panel_set_varying(&pan, &mod);
    if (err) {
        goto bailout;
    }

    if (opt & OPT_U) {
        nerlove = (opt & OPT_E)? 1 : 0;
    }

    calculate_Tbar(&pan);

    if (opt & OPT_U) {
        /* trying to do random effects */
        int xdf = pan.effn - (mod.ncoeff - pan.ntdum);

        if (xdf <= 0 && !nerlove) {
            gretl_errmsg_set(_("Couldn't estimate group means regression"));
            err = mod.errcode = E_DF;
            goto bailout;
        } else if (pan.opt & OPT_M) {
            /* matrix version of Hausman test */
            err = matrix_hausman_allocate(&pan);
            if (err) {
                goto bailout;
            }
        }
    }

    err = within_variance(&pan, dset, prn);
    if (err) {
        goto bailout;
    }

    if ((opt & OPT_U) && !na(pan.s2e)) {
        DATASET *gset;

        breusch_pagan_LM(&pan, prn);

        gset = group_means_dataset(&pan, dset);
        if (gset == NULL) {
            err = E_ALLOC;
        } else {
            err = between_variance(&pan, gset);
        }

        if (err && !nerlove) {
            pputs(prn, _("Couldn't estimate group means regression\n"));
        } else {
            err = random_effects(&pan, dset, gset, prn);
            if (!err) {
                save_breusch_pagan_result(&pan);
                finalize_hausman_test(&pan, PANEL, prn);
            }
        }

        if (gset != NULL) {
            destroy_dataset(gset);
        }
    }

 bailout:

    if (!err) {
        if ((pan_opt & (OPT_F | OPT_U)) && mod.missmask != NULL) {
            /* preserve missing obs mask, if any */
            char *mask = mod.missmask;

            mod.missmask = NULL;
            clear_model(&mod);
            mod = *pan.realmod;
            mod.missmask = mask;
        } else if (!(opt & OPT_P)) {
            /* not doing pooled OLS, in which case @mod would
               already hold the correct model
            */
            clear_model(&mod);
            mod = *pan.realmod;
        }
        gretl_model_smpl_init(&mod, dset);
        if (opt & OPT_D) {
            if (pan.ntdum > 0) {
                gretl_model_set_int(&mod, "ntdum", pan.ntdum);
                maybe_suppress_time_dummies(&mod, pan.ntdum);
            }
            if (complex_subsampled()) {
                process_time_dummies(&mod, dset, orig_v);
            }
        }
    }

    panelmod_free(&pan);

    if (err && mod.errcode == 0) {
        mod.errcode = err;
    }

    if (complex_subsampled()) {
        dataset_drop_last_variables(dset, dset->v - orig_v);
    }

#if PDEBUG
    fprintf(stderr, "real_panel_model return: mod.errcode = %d\n",
	    mod.errcode);
#endif

    return mod;
}

/* Called from qr_estimate.c in case robust VCV estimation is called
   for in the context of TSLS estimation on panel data.
*/

int panel_tsls_robust_vcv (MODEL *pmod, DATASET *dset)
{
    panelmod_t pan;
    int err = 0;

    panelmod_init(&pan);

    err = panelmod_setup(&pan, pmod, dset, 0, OPT_NONE);
    if (!err) {
        err = panel_robust_vcv(pmod, &pan, dset, dset);
    }

    panelmod_free(&pan);

    return err;
}

/* write weights for groupwise weighted least squares into the
   last variable in the dataset */

static int
write_weights_to_dataset (double *uvar, int nunits, int T,
                          DATASET *dset)
{
    int w = dset->v - 1;
    double wi;
    int i, t;

    for (i=0; i<nunits; i++) {
        if (uvar[i] <= 0.0 || na(uvar[i])) {
            wi = 0.0;
        } else {
            wi = 1.0 / uvar[i]; /* sqrt? */
        }
        for (t=0; t<T; t++) {
            dset->Z[w][panel_index(i, t)] = wi;
        }
    }

    return 0;
}

static int allocate_weight_var (DATASET *dset)
{
    if (dataset_add_series(dset, 1)) {
        return E_ALLOC;
    }

    strcpy(dset->varname[dset->v - 1], "unit_wt");

    return 0;
}

static double max_coeff_diff (const MODEL *pmod, const double *bvec)
{
    double diff, maxdiff = 0.0;
    int i;

    for (i=0; i<pmod->ncoeff; i++) {
        diff = fabs(pmod->coeff[i] - bvec[i]);
        if (diff > maxdiff) {
            maxdiff = diff;
        }
    }

    return maxdiff;
}

#define S2MINOBS 1

/* Wald test for groupwise heteroskedasticity, without assuming
   normality of the errors: see Greene, 4e, p. 598.  Note that the
   computation involves a factor of (1/(T_i - 1)) so we cannot
   include groups that have only one observation.
*/

static double
wald_hetero_test (const MODEL *pmod, double s2,
                  const double *uvar, panelmod_t *pan,
                  int *df)
{
    double x, W = 0.0;
    int i, t, Ti;

    *df = 0;

    for (i=0; i<pan->nunits; i++) {
        double fii = 0.0;

        Ti = pan->unit_obs[i];
        if (Ti < 2) {
            continue;
        }

        for (t=0; t<pan->T; t++) {
            x = pmod->uhat[panel_index(i, t)];
            if (!na(x)) {
                x = x * x - uvar[i];
                fii += x * x;
            }
        }

        if (fii <= 0) {
            W = NADBL;
            break;
        }

        fii *= (1.0 / Ti) * (1.0 / (Ti - 1.0));
        x = uvar[i] - s2;
        W += x * x / fii;
        *df += 1;
    }

    if (*df < 2) {
        W = NADBL;
    }

    return W;
}

/* Likelihood-ratio test for groupwise heteroskedasticity: see Greene,
   4e, pp. 597 and 599.  This test requires that we're able to
   calculate the ML estimates (using iterated WLS).
*/

static int
ml_hetero_test (MODEL *pmod, double s2, const double *uvar,
                int nunits, const int *unit_obs)
{
    ModelTest *test;
    double x2, s2h = 0.0;
    int i, Ti, df = 0;
    int err = 0;

    for (i=0; i<nunits; i++) {
        Ti = unit_obs[i];
        if (Ti > 0) {
            s2h += Ti * log(uvar[i]);
            df++;
        }
    }

    x2 = pmod->nobs * log(s2) - s2h;
    df--;

    test = model_test_new(GRETL_TEST_GROUPWISE);
    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_LR);
        model_test_set_dfn(test, df);
        model_test_set_value(test, x2);
        model_test_set_pvalue(test, chisq_cdf_comp(df, x2));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = 1;
    }

    return err;
}

static double pooled_ll (const MODEL *pmod)
{
    double n = pmod->nobs;

    return -(n / 2.0) * (1.0 + LN_2_PI - log(n) + log(pmod->ess));
}

/* calculate log-likelihood for panel MLE with groupwise
   heteroskedasticity
*/

static void panel_ML_ll (MODEL *pmod, const double *uvar,
                         int nunits, const int *unit_obs)
{
    double ll = -(pmod->nobs / 2.0) * LN_2_PI;
    int i, Ti;

    for (i=0; i<nunits; i++) {
        Ti = unit_obs[i];
        if (Ti > 0) {
            ll -= (Ti / 2.0) * (1.0 + log(uvar[i]));
        }
    }

    pmod->lnL = ll;
    mle_criteria(pmod, 0);
}

/* we can't estimate a group-specific variance based on just one
   observation */

static int singleton_check (const int *unit_obs, int nunits)
{
    int i;

    for (i=0; i<nunits; i++) {
        if (unit_obs[i] == 1) {
            return 1;
        }
    }

    return 0;
}

/* compute per-unit error variances */

static void unit_error_variances (double *uvar, const MODEL *pmod,
                                  panelmod_t *pan, int *df)
{
    int i, t;
    double uit;

    *df = 0;

    for (i=0; i<pan->nunits; i++) {
        if (pan->unit_obs[i] < S2MINOBS) {
            uvar[i] = NADBL;
        } else {
            uvar[i] = 0.0;
            for (t=0; t<pan->T; t++) {
                uit = pmod->uhat[panel_index(i, t)];
                if (!na(uit)) {
                    uvar[i] += uit * uit;
                }
            }
            uvar[i] /= pan->unit_obs[i];
            *df += 1;
        }
    }
}

static void print_unit_variances (panelmod_t *pan, double *uvar,
                                  int wald_test, PRN *prn)
{
    int i, Ti;

    pputs(prn, " unit    variance\n");

    for (i=0; i<pan->nunits; i++) {
        Ti = pan->unit_obs[i];
        if (na(uvar[i]) || (wald_test && Ti < 2)) {
            pprintf(prn, "%5d%12s (T = %d)\n", i+1, "NA", Ti);
        } else {
            pprintf(prn, "%5d%#12g (T = %d)\n", i+1, uvar[i], Ti);
        }
    }
}

#define SMALLDIFF 0.0001
#define WLS_MAX   30

/**
 * panel_wls_by_unit:
 * @list: regression list.
 * @dset: dataset struct.
 * @opt: may include %OPT_I (iterate to ML solution),
 * %OPT_V for verbose operation.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Performs weighted least squares (possibly iterated) on
 * a panel model, allowing for groupwise heteroskedasticity.
 *
 * Returns: the estimates, in a %MODEL.
 */

MODEL panel_wls_by_unit (const int *list, DATASET *dset,
                         gretlopt opt, PRN *prn)
{
    MODEL mdl;
    panelmod_t pan;
    gretlopt wlsopt = OPT_A | OPT_Z;
    double *uvar = NULL;
    double *bvec = NULL;
    double s2, diff = 1.0;
    int *wlist = NULL;
    int orig_v = dset->v;
    int i, iter = 0;

    gretl_error_clear();
    panelmod_init(&pan);

    if (opt & OPT_I) {
        /* iterating: no degrees-of-freedom correction in WLS */
        wlsopt |= OPT_N;
    }

    /* baseline pooled model */
    mdl = lsq(list, dset, OLS, OPT_A);
    if (mdl.errcode) {
        goto bailout;
    }

    mdl.errcode = panelmod_setup(&pan, &mdl, dset, 0, OPT_NONE);
    if (mdl.errcode) {
        goto bailout;
    }

    uvar = malloc(pan.nunits * sizeof *uvar);
    if (uvar == NULL) {
        free(pan.unit_obs);
        mdl.errcode = E_ALLOC;
        return mdl;
    }

    if (opt & OPT_I) {
	/* the --iterate option */
        if (singleton_check(pan.unit_obs, pan.nunits)) {
            pprintf(prn, _("Can't produce ML estimates: "
                           "some units have only one observation"));
            pputc(prn, '\n');
            opt ^= OPT_I;
        }
    }

    s2 = mdl.ess / mdl.nobs;

    if ((opt & OPT_V) && (opt & OPT_I)) {
        pprintf(prn, _("\nOLS error variance = %g\n"), s2);
        pprintf(prn, "log-likelihood = %g\n", pooled_ll(&mdl));
    }

    if (allocate_weight_var(dset)) {
        mdl.errcode = E_ALLOC;
        goto bailout;
    }

    if (opt & OPT_I) {
        bvec = malloc(mdl.ncoeff * sizeof *bvec);
        if (bvec == NULL) {
            mdl.errcode = E_ALLOC;
            goto bailout;
        }
    }

    /* allocate and construct WLS regression list */
    wlist = gretl_list_new(mdl.list[0] + 1);
    if (wlist == NULL) {
        mdl.errcode = E_ALLOC;
        goto bailout;
    }
    wlist[1] = dset->v - 1; /* weight variable: the last var added */
    for (i=2; i<=wlist[0]; i++) {
        wlist[i] = mdl.list[i-1];
    }

    /* If wanted (and possible) iterate to ML solution; otherwise just do
       one-step FGLS estimation.
    */
    while (diff > SMALLDIFF) {
        int df = 0;

        iter++;
        unit_error_variances(uvar, &mdl, &pan, &df);

        if (opt & OPT_V) {
            if (opt & OPT_I) {
                pprintf(prn, "\n*** %s %d ***\n", _("iteration"),
                        iter);
            } else {
                pputc(prn, '\n');
            }
            print_unit_variances(&pan, uvar, 0, prn);
        }

        write_weights_to_dataset(uvar, pan.nunits, pan.T, dset);

        if (opt & OPT_I) {
            /* save coefficients for comparison */
            for (i=0; i<mdl.ncoeff; i++) {
                bvec[i] = mdl.coeff[i];
            }
        }

        clear_model(&mdl);
        mdl = lsq(wlist, dset, WLS, wlsopt);
        if (mdl.errcode) {
            break;
        } else if (iter > WLS_MAX) {
            mdl.errcode = E_NOCONV;
            break;
        }

        if (opt & OPT_I) {
            diff = max_coeff_diff(&mdl, bvec);
            if ((opt & OPT_V) && iter == 1) {
                pprintf(prn, _("\nFGLS pooled error variance = %g\n"),
                        mdl.ess / mdl.nobs);
            }
        } else {
            /* one-step FGLS */
            break;
        }
    }

    if (!mdl.errcode) {
        mdl.ci = PANEL;
        if (!(opt & OPT_A)) {
            set_model_id(&mdl, pan.opt);
        }
        gretl_model_set_int(&mdl, "n_included_units", pan.effn);
        mdl.opt |= OPT_H; /* --unit-weights */
        mdl.nwt = 0;

        if (opt & OPT_I) {
            int df = 0;

            gretl_model_set_int(&mdl, "iters", iter);
            ml_hetero_test(&mdl, s2, uvar, pan.nunits, pan.unit_obs);
            unit_error_variances(uvar, &mdl, &pan, &df);
            panel_ML_ll(&mdl, uvar, pan.nunits, pan.unit_obs);
            if (opt & OPT_V) {
                pputc(prn, '\n');
            }
        }
    }

 bailout:

    free(pan.unit_obs);
    free(uvar);
    free(wlist);
    free(bvec);

    if (!(opt & OPT_V)) {
        dataset_drop_last_variables(dset, dset->v - orig_v);
    }

    return mdl;
}

static void print_wald_test (double W, int df, double pval,
                             panelmod_t *pan, double *uvar,
                             double s2, PRN *prn)
{
    pprintf(prn, "\n%s:\n",
            _("Distribution free Wald test for heteroskedasticity"));
    pprintf(prn, " %s(%d) = %g, ",  _("Chi-square"), df, W);
    pprintf(prn, "%s = %g\n\n", _("with p-value"), pval);

    if (pan->nunits <= 30) {
        pprintf(prn, "%s = %g\n\n", _("Pooled error variance"), s2);
        print_unit_variances(pan, uvar, 1, prn);
    }
}

/**
 * groupwise_hetero_test:
 * @pmod: panel model to be tested.
 * @dset: information on the (panel) data set.
 * @opt: may contain OPT_S to attach the test result
 * to @pmod, OPT_I for silent operation.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Performs a Wald test for the null hypothesis that the
 * error variance is uniform across the units.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int groupwise_hetero_test (MODEL *pmod, DATASET *dset,
                           gretlopt opt, PRN *prn)
{
    panelmod_t pan;
    double *uvar = NULL;
    double s2, W = NADBL;
    int df, err = 0;

    if (pmod->ci == OLS || (pmod->ci == PANEL && (pmod->opt & OPT_F))) {
        ; /* OK for pooled or fixed effects */
    } else if (0 && pmod->ci == PANEL && (pmod->opt & OPT_U)) {
        ; /* just testing! */
    } else {
        return E_NOTIMP;
    }

    gretl_error_clear();
    panelmod_init(&pan);

    err = panelmod_setup(&pan, pmod, dset, 0, OPT_NONE);
    if (err) {
        return err;
    }

    uvar = malloc(pan.nunits * sizeof *uvar);
    if (uvar == NULL) {
        free(pan.unit_obs);
        return E_ALLOC;
    }

    s2 = pmod->ess / pmod->nobs;

    unit_error_variances(uvar, pmod, &pan, &df);

    if (df >= 2) {
        W = wald_hetero_test(pmod, s2, uvar, &pan, &df);
    }

    if (na(W)) {
        err = E_DATA;
    } else {
        double pval = chisq_cdf_comp(df, W);

        if (!(opt & OPT_I)) {
            print_wald_test(W, df, pval, &pan, uvar, s2, prn);
        }

        if (opt & OPT_S) {
            ModelTest *test = model_test_new(GRETL_TEST_GROUPWISE);

            if (test != NULL) {
                model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
                model_test_set_dfn(test, df);
                model_test_set_value(test, W);
                model_test_set_pvalue(test, pval);
                maybe_add_test_to_model(pmod, test);
            }
        }

        record_test_result(W, pval);
    }

    free(pan.unit_obs);
    free(uvar);

    return err;
}

static int print_ar_aux_model (MODEL *pmod, DATASET *dset,
                               int j, PRN *prn)
{
    const char *heads[] = {
        N_("First differenced equation (dependent, d_y)"),
        N_("Autoregression of residuals (dependent, uhat)"),
        N_("Auxiliary regression including lagged residual"),
    };
    gretl_matrix *cse;
    gretl_array *S;
    int i, vi, err = 0;

    cse = gretl_matrix_alloc(pmod->ncoeff, 2);
    S = gretl_array_new(GRETL_TYPE_STRINGS, pmod->ncoeff, &err);

    if (cse == NULL || S == NULL) {
        return E_ALLOC;
    }

    for (i=0; i<pmod->ncoeff; i++) {
        gretl_matrix_set(cse, i, 0, pmod->coeff[i]);
        gretl_matrix_set(cse, i, 1, pmod->sderr[i]);
        vi = pmod->list[i+2];
        gretl_array_set_string(S, i, dset->varname[vi], 0);
    }

    if (j >= 0) {
        pprintf(prn, "%s:\n", _(heads[j]));
    }
    print_model_from_matrices(cse, NULL, S, pmod->dfd, OPT_NONE, prn);
    pprintf(prn, "  n = %d, R-squared = %.4f\n\n", pmod->nobs, pmod->rsq);

    gretl_matrix_free(cse);
    gretl_array_nullify_elements(S);
    gretl_array_destroy(S);

    return 0;
}

static void finalize_autocorr_test (MODEL *pmod, ModelTest *test,
                                    gretlopt opt, PRN *prn)
{
    if (!(opt & OPT_I)) {
        if (opt & OPT_Q) {
            pputc(prn, '\n');
        }
        gretl_model_test_print_direct(test, 1, prn);
    }
    if (opt & OPT_S) {
        maybe_add_test_to_model(pmod, test);
    } else {
        free(test);
    }
}

/* See Wooldridge's Econometric Analysis of Cross Section
   and Panel Data (MIT Press, 2002), pages 176-7. Test
   for first order autocorrelation in the context of pooled OLS.
*/

static int pooled_autocorr_test (MODEL *pmod, DATASET *dset,
                                 gretlopt opt, PRN *prn)
{
    int quiet = (opt & OPT_Q);
    int *ulist = NULL;
    int i, vi = 0;
    int err = 0;

    if (pmod->full_n != dset->n) {
        return E_DATA;
    }

    err = dataset_add_NA_series(dset, 1);
    if (!err) {
        vi = dset->v - 1;
        strcpy(dset->varname[vi], "uhat");
        for (i=0; i<dset->n; i++) {
            dset->Z[vi][i] = pmod->uhat[i];
        }
    }

    if (!err) {
        ulist = gretl_list_new(pmod->ncoeff + 2);
        if (ulist == NULL) {
            err = E_ALLOC;
        } else {
            for (i=1; i<=pmod->list[0]; i++) {
                ulist[i] = pmod->list[i];
            }
            /* append lagged residual to list */
            ulist[i] = laggenr(vi, 1, dset);
            if (ulist[i] < 0) {
                err = E_DATA;
            } else {
                strcpy(dset->varname[ulist[i]], "uhat(-1)");
            }
        }
    }

    if (!err) {
        gretlopt auxopt = OPT_P | OPT_R | OPT_Q;
        MODEL tmp;

        /* estimate model including lagged residual */
        tmp = real_panel_model(ulist, dset, auxopt, NULL);
        err = tmp.errcode;
        if (!err && !quiet) {
            if (gretl_echo_on()) {
                pputc(prn, '\n');
            }
            print_ar_aux_model(&tmp, dset, 2, prn);
        }
        if (!err) {
            double tstat, pval;
            int df = tmp.dfd;

            i = tmp.ncoeff - 1;
            tstat = tmp.coeff[i] / tmp.sderr[i];
            pval = student_pvalue_2(df, tstat);
            record_test_result(tstat, pval);

            if ((opt & OPT_S) || !(opt & OPT_I)) {
                /* saving the test onto @pmod, or not silent */
                ModelTest *test = model_test_new(GRETL_TEST_PANEL_AR);

                if (test != NULL) {
                    model_test_set_teststat(test, GRETL_STAT_STUDENT);
                    model_test_set_value(test, tstat);
                    model_test_set_dfn(test, df);
                    model_test_set_pvalue(test, pval);
                    finalize_autocorr_test(pmod, test, opt, prn);
                }
            }
        }
        clear_model(&tmp);
    }

    free(ulist);

    return err;
}

/* See Wooldridge's Econometric Analysis of Cross Section
   and Panel Data (MIT Press, 2002), pages 282-3. Test
   for first order autocorrelation based on the residuals
   from the first-differenced model.
*/

static int wooldridge_autocorr_test (MODEL *pmod, DATASET *dset,
                                     gretlopt opt, PRN *prn)
{
    MODEL tmp;
    gretlopt tmp_opt;
    int quiet = (opt & OPT_Q);
    int *dlist = NULL;
    int i, j, vi;
    int clearit = 0;
    int err = 0;

    dlist = gretl_list_new(pmod->ncoeff + 1);
    if (dlist == NULL) {
        return E_ALLOC;
    }

    dlist[0] = 0;
    for (i=1, j=1; i<=pmod->list[0] && !err; i++) {
        vi = pmod->list[i];
        if (vi != 0 && !gretl_isconst(dset->t1, dset->t2, dset->Z[vi])) {
            dlist[j] = diffgenr(vi, DIFF, dset);
            if (dlist[j] < 0) {
                err = E_DATA;
            } else {
                j++;
                dlist[0] += 1;
            }
        } else if (i == 1) {
            gretl_errmsg_set(_("The dependent variable is constant"));
            err = E_DATA;
        }
    }

    /* aux panel models: pooled, robust, quiet */
    tmp_opt = OPT_P | OPT_R | OPT_Q;

    if (!err) {
        /* estimate model in first-differenced form */
        tmp = real_panel_model(dlist, dset, tmp_opt, NULL);
        err = tmp.errcode;
        if (!err && !quiet) {
            if (gretl_echo_on()) {
                pputc(prn, '\n');
            }
            print_ar_aux_model(&tmp, dset, 0, prn);
        }
    }

    if (!err) {
        err = dataset_add_allocated_series(dset, tmp.uhat);
        tmp.uhat = NULL;
        clear_model(&tmp);
    }

    if (!err) {
        vi = dset->v - 1; /* the last series added */
        strcpy(dset->varname[vi], "uhat");
        dlist[0] = 2;
        dlist[1] = vi;
        dlist[2] = laggenr(vi, 1, dset);
        if (dlist[2] < 0) {
            err = E_DATA;
        } else {
            /* regress residual on its first lag */
            tmp = real_panel_model(dlist, dset, tmp_opt, NULL);
            err = tmp.errcode;
            if (!err && !quiet) {
                strcpy(dset->varname[dlist[2]], "uhat(-1)");
                print_ar_aux_model(&tmp, dset, 1, prn);
            }
            clearit = 1;
        }
    }

    if (!err) {
        double c, s, F, pval;

        c = tmp.coeff[0] + 0.5;
        s = tmp.sderr[0];
        F = c * c / (s * s);
        pval = snedecor_cdf_comp(1, tmp.dfd, F);
        record_test_result(F, pval);

        if ((opt & OPT_S) || !(opt & OPT_I)) {
            /* saving the test onto @pmod, or not silent */
            ModelTest *test = model_test_new(GRETL_TEST_PANEL_AR);

            if (test != NULL) {
                model_test_set_teststat(test, GRETL_STAT_F);
                model_test_set_value(test, F);
                model_test_set_dfn(test, 1);
                model_test_set_dfd(test, tmp.dfd);
                model_test_set_pvalue(test, pval);
                finalize_autocorr_test(pmod, test, opt, prn);
            }
        }
    }

    if (clearit) {
        clear_model(&tmp);
    }
    free(dlist);

    return err;
}

int panel_autocorr_test (MODEL *pmod, DATASET *dset,
                         gretlopt opt, PRN *prn)
{
    int save_robust = libset_get_int(PANEL_ROBUST);
    int orig_v = dset->v;
    int err;

    libset_set_int(PANEL_ROBUST, ARELLANO); /* ? */

    if (pmod->ci == OLS || (pmod->opt & OPT_P)) {
        err = pooled_autocorr_test(pmod, dset, opt, prn);
    } else {
        err = wooldridge_autocorr_test(pmod, dset, opt, prn);
    }

    libset_set_int(PANEL_ROBUST, save_robust);
    dataset_drop_last_variables(dset, dset->v - orig_v);

    return err;
}

/* test for cross-sectional dependence */

int panel_xdepend_test (MODEL *pmod, DATASET *dset,
                        gretlopt opt, PRN *prn)
{
    const double *u;
    double rij, rsum = 0.0;
    double arsum = 0.0;
    double ssx, ssy, sxy;
    double xbar, ybar;
    int N1, N2, T, Tij;
    int i, j, t, si, sj;
    int N = 0, Nr = 0;
    int err = 0;

    if (dset->structure != STACKED_TIME_SERIES) {
        return E_PDWRONG;
    } else if (pmod->uhat == NULL) {
        return E_DATA;
    }

    T = dset->pd;
    N1 = pmod->t1 / T;
    N2 = pmod->t2 / T;
    u = pmod->uhat;

    for (i=N1; i<N2; i++) {
        int Nj = 0;

        for (j=i+1; j<=N2; j++) {
            xbar = ybar = 0.0;
            Tij = 0;
            for (t=0; t<T; t++) {
                si = i * T + t;
                sj = j * T + t;
                if (!na(u[si]) && !na(u[sj])) {
                    Tij++;
                    xbar += u[si];
                    ybar += u[sj];
                }
            }
            if (Tij >= 2) {
                ssx = ssy = sxy = 0.0;
                xbar /= Tij;
                ybar /= Tij;
                for (t=0; t<T; t++) {
                    si = i * T + t;
                    sj = j * T + t;
                    if (!na(u[si]) && !na(u[sj])) {
                        ssx += (u[si] - xbar) * (u[si] - xbar);
                        ssy += (u[sj] - ybar) * (u[sj] - ybar);
                        sxy += (u[si] - xbar) * (u[sj] - ybar);
                    }
                }
                rij = sxy / sqrt(ssx * ssy);
                rsum += sqrt(Tij) * rij;
                arsum += fabs(rij);
                Nr++;
                Nj++;
            }
        }
        if (Nj > 0) {
            N++;
        }
    }

    if (N == 0) {
        err = E_TOOFEW;
    }

    if (!err) {
        double CD, pval;

        N = N + 1;
        CD = sqrt(2.0 / (N * (N - 1.0))) * rsum;
        pval = normal_pvalue_2(CD);

        if (!(opt & OPT_I)) {
            pputs(prn, _("Pesaran CD test for cross-sectional dependence"));
            pprintf(prn, "\n%s: z = %f,\n", _("Test statistic"), CD);
            pprintf(prn, "%s = P(|z| > %g) = %.3g\n", _("with p-value"),
                    CD, pval);
            pprintf(prn, _("Average absolute correlation = %.3f"), arsum / Nr);
            pputc(prn, '\n');
        }

        if (opt & OPT_S) {
            ModelTest *test = model_test_new(GRETL_TEST_XDEPEND);

            if (test != NULL) {
                model_test_set_teststat(test, GRETL_STAT_Z);
                model_test_set_value(test, CD);
                model_test_set_pvalue(test, pval);
                maybe_add_test_to_model(pmod, test);
            }
        }

        record_test_result(CD, pval);
    }

    return err;
}

/**
 * switch_panel_orientation:
 * @dset: pointer to dataset struct.
 *
 * Reorganizes @dset, transforming from stacked cross sections to
 * stacked time series or vice versa.  If the transformation is from
 * stacked time series to stacked cross section, the dataset will
 * no longer be acceptable as a panel for gretl's purposes; it
 * may be useful for export purposes, though.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int switch_panel_orientation (DATASET *dset)
{
    int oldmode = dset->structure;
    double *tmp;
    char **markers = NULL;
    double pdx;
    int T, n;
    int i, j, s, t;

    if (dset->Z == NULL) {
        return E_NODATA;
    } else if (oldmode != STACKED_TIME_SERIES &&
        oldmode != STACKED_CROSS_SECTION) {
        return E_DATA;
    }

    tmp = malloc(dset->n * sizeof *tmp);
    if (tmp == NULL) {
        return E_ALLOC;
    }

    if (oldmode == STACKED_CROSS_SECTION) {
        n = dset->pd;
        T = dset->n / n;
    } else {
        T = dset->pd;
        n = dset->n / T;
    }

    /* copy the data series across in transformed order */
    for (i=1; i<dset->v; i++) {
        for (t=0; t<dset->n; t++) {
            /* transcribe to tmp in original order */
            tmp[t] = dset->Z[i][t];
        }
        s = 0;
        if (oldmode == STACKED_CROSS_SECTION) {
            /* convert to stacked time-series */
            for (j=0; j<n; j++) {
                for (t=0; t<T; t++) {
                    dset->Z[i][s++] = tmp[t * n + j];
                }
            }
        } else {
            /* convert to stacked cross-sections */
            for (t=0; t<T; t++) {
                for (j=0; j<n; j++) {
                    dset->Z[i][s++] = tmp[j * T + t];
                }
            }
        }
    }

    /* rearrange observations markers if relevant */
    if (dset->S != NULL) {
        markers = strings_array_new_with_length(dset->n, OBSLEN);
        if (markers != NULL) {
            for (t=0; t<dset->n; t++) {
                strcpy(markers[t], dset->S[t]);
            }
            s = 0;
            if (oldmode == STACKED_CROSS_SECTION) {
                for (j=0; j<n; j++) {
                    for (t=0; t<T; t++) {
                        strcpy(dset->S[s++], markers[t * n + j]);
                    }
                }
            } else {
                for (t=0; t<T; t++) {
                    for (j=0; j<n; j++) {
                        strcpy(dset->S[s++], markers[j * T + t]);
                    }
                }
            }
            strings_array_free(markers, dset->n);
        } else {
            /* should we flag an error? */
            dataset_destroy_obs_markers(dset);
        }
    }

    dset->sd0 = 1.0;
    pdx = 0.1;

    /* change the datainfo setup */
    if (oldmode == STACKED_CROSS_SECTION) {
        dset->structure = STACKED_TIME_SERIES;
        dset->pd = T;
        while (T /= 10) {
            pdx *= 0.1;
        }
    } else {
        dset->structure = STACKED_CROSS_SECTION;
        dset->pd = n;
        while (n /= 10) {
            pdx *= 0.1;
        }
    }

    dset->sd0 += pdx;
    ntolabel(dset->stobs, 0, dset);
    ntolabel(dset->endobs, dset->n - 1, dset);

    free(tmp);

    return 0;
}

/**
 * panel_isconst:
 * @t1: starting observation.
 * @t2: ending observation.
 * @pd: panel time-series length.
 * @x: data series to examine.
 * @bygroup: use 1 to check for constancy across groups,
 * 0 for constancy across time.
 *
 * Check whether series @x is constant (either over time or
 * across groups) in a panel dataset. Missing values are
 * ignored.
 *
 * Returns: 1 if the variable is constant, otherwise 0.
 */

int panel_isconst (int t1, int t2, int pd, const double *x,
                   int bygroup)
{
    double x0 = NADBL;
    int t, ret = 1;

    if (bygroup) {
        /* check for variation across groups */
        int tref;

        for (t=t1; t<=t2 && ret; t++) {
            if (na(x[t])) {
                continue;
            }
            tref = t - pd;
            while (tref >= t1) {
                /* check last obs for the same period */
                x0 = x[t-pd];
                if (na(x0)) {
                    tref -= pd;
                } else {
                    if (x[t] != x0) {
                        ret = 0;
                    }
                    break;
                }
            }
        }
    } else {
        /* check for time-variation */
        int u, ubak = -1;

        for (t=t1; t<=t2 && ret; t++) {
            if (na(x[t])) {
                continue;
            }
            u = t / pd;
            if (u == ubak) {
                /* same group */
                if (!na(x0) && x[t] != x0) {
                    ret = 0;
                }
            } else {
                /* starting a new group */
                x0 = x[t];
                ubak = u;
            }
        }
    }

    return ret;
}

static int varying_vars_list (const DATASET *dset, panelmod_t *pan)
{
    int i, j, k, t, bigt;
    int err = 0;

    pan->vlist = gretl_list_new(pan->pooled->list[0]);
    if (pan->vlist == NULL) {
        return E_ALLOC;
    }

    pan->vlist[0] = 1;
    pan->vlist[1] = pan->pooled->list[1];

    k = 2;

    for (j=2; j<=pan->pooled->list[0]; j++) {
        int vj = pan->pooled->list[j];
        int varies = 0;
        double x0;

        if (vj == 0) {
            /* const */
            pan->vlist[k++] = 0;
            pan->vlist[0] += 1;
            continue;
        }

        for (i=0; i<pan->nunits && !varies; i++) {
            if (pan->unit_obs[i] == 0) {
                continue;
            }
            x0 = NADBL;
            for (t=0; t<pan->T && !varies; t++) {
                bigt = panel_index(i, t);
                if (panel_missing(pan, bigt)) {
                    continue;
                }
                if (na(x0)) {
                    x0 = dset->Z[vj][bigt];
                } else if (dset->Z[vj][bigt] != x0) {
                    varies = 1;
                }
            }
        }

        if (varies) {
            pan->vlist[k++] = vj;
            pan->vlist[0] += 1;
        } else {
            fprintf(stderr, "Variable %d '%s' is time-invariant\n",
                    vj, dset->varname[vj]);
        }
    }

#if PDEBUG
    printlist(pan->pooled->list, "original regressors");
    printlist(pan->vlist, "time-varying regressors");
#endif

    return err;
}

/* apparatus for setting panel structure by reading two variables,
   indexing the cross-sectional unit and the time period
   respectively
*/

typedef struct s_point_t s_point;
typedef struct sorter_t sorter;

struct s_point_t {
    int obsnum;
    double val1;
    double val2;
};

struct sorter_t {
    int sortvar1;
    int sortvar2;
    s_point *points;
};

static void sorter_fill_points (sorter *s, const double **Z,
                                int n)
{
    int t;

    for (t=0; t<n; t++) {
        s->points[t].obsnum = t;
        s->points[t].val1 = Z[s->sortvar1][t];
        s->points[t].val2 = Z[s->sortvar2][t];
    }
}

static int compare_obs (const void *a, const void *b)
{
    s_point *pa = (s_point *) a;
    s_point *pb = (s_point *) b;
    int ret;

    ret = (pa->val1 > pb->val1) - (pa->val1 < pb->val1);
    if (ret == 0) {
        ret = (pa->val2 > pb->val2) - (pa->val2 < pb->val2);
    }

    if (ret == 0) {
        /* mark an error: got a duplicated value */
        pa->obsnum = pb->obsnum = -1;
    }

    return ret;
}

static int check_indices (sorter *s, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (s->points[i].obsnum < 0) {
            gretl_errmsg_sprintf(_("Error: unit %g, period %g: duplicated observation"),
                                 s->points[i].val1, s->points[i].val2);
            return E_DATA;
        }
    }

    return 0;
}

static int panel_is_sorted (DATASET *dset, int uv, int tv)
{
    const double *unit = dset->Z[uv];
    const double *period = dset->Z[tv];
    int t, ret = 1;

    for (t=1; t<dset->n && ret; t++) {
        if (unit[t] < unit[t-1]) {
            /* the unit variable must be non-decreasing */
            ret = 0;
        } else if (unit[t] == unit[t-1] && period[t] <= period[t-1]) {
            /* the period variable must be increasing within the
               observations for each unit */
            ret = 0;
        }
    }

    return ret;
}

static int panel_data_sort_by (DATASET *dset, int uv, int tv, int *ustrs)
{
    int n = dset->n;
    char **S = NULL;
    double *tmp = NULL;
    int i, t;
    sorter s;
    int err = 0;

    tmp = malloc(n * sizeof *tmp);
    if (tmp == NULL) {
        return E_ALLOC;
    }

    s.points = malloc(n * sizeof *s.points);
    if (s.points == NULL) {
        free(tmp);
        return E_ALLOC;
    }

    if (dset->S != NULL) {
        S = strings_array_new_with_length(n, OBSLEN);
        if (S == NULL) {
            free(tmp);
            free(s.points);
            return E_ALLOC;
        }
    }

    s.sortvar1 = uv;
    s.sortvar2 = tv;

    sorter_fill_points(&s, (const double **) dset->Z, n);

    qsort((void *) s.points, (size_t) n, sizeof s.points[0],
          compare_obs);

    err = check_indices(&s, n);
    if (err) {
        goto bailout;
    }

    for (i=1; i<dset->v; i++) {
        for (t=0; t<n; t++) {
            tmp[t] = dset->Z[i][t];
        }
        for (t=0; t<n; t++) {
            dset->Z[i][t] = tmp[s.points[t].obsnum];
        }
    }

    if (S != NULL) {
        *ustrs = 1;
        for (t=0; t<n; t++) {
            strcpy(S[t], dset->S[t]);
            if (S[t][0] == '\0') {
                *ustrs = 0;
            }
        }
        for (t=0; t<n; t++) {
            strcpy(dset->S[t], S[s.points[t].obsnum]);
        }
        for (t=1; t<n && *ustrs; t++) {
            if (dset->Z[uv][t] == dset->Z[uv][t-1] &&
                strcmp(dset->S[t], dset->S[t-1])) {
                *ustrs = 0;
            }
        }
        strings_array_free(S, n);
    }

 bailout:

    free(s.points);
    free(tmp);

    return err;
}

/* Given the variables coding for panel unit and time period,
   construct parallel arrays of int that code for the same
   information, but are zero-based and consecutive.
*/

static int normalize_uid_tid (const double *tid, int T,
                              const DATASET *dset,
                              int uv, int tv,
                              int **pnuid, int **pntid)
{
    int *nuid = NULL;
    int *ntid = NULL;
    int ui, tmin;
    int i, t;

    nuid = malloc(dset->n * sizeof *nuid);
    ntid = malloc(dset->n * sizeof *ntid);
    if (nuid == NULL || ntid == NULL) {
        free(nuid);
        free(ntid);
        return E_ALLOC;
    }

    ui = tmin = 0;

    for (i=0; i<dset->n; i++) {
        if (i > 0 && dset->Z[uv][i] > dset->Z[uv][i-1]) {
            tmin = 0;
            ui++;
        } else if (i > 0) {
            tmin++;
        }
        nuid[i] = ui;
        for (t=tmin; t<T; t++) {
            if (dset->Z[tv][i] == tid[t]) {
                ntid[i] = t;
                break;
            }
        }
    }

    *pnuid = nuid;
    *pntid = ntid;

    return 0;
}

/* Handle the case where a sub-sampled panel dataset has been padded
   with rows containing NAs to recreate a (nominally) balanced panel
   structure. We have to get rid of the padding before we try
   restoring the full data range. Note that while this function
   "shrinks" @dset in the sense of reducing the series length as
   recorded in dset->n, it does not "physically" shrink the array
   members of dset->Z.
*/

int undo_panel_padding (DATASET *dset)
{
    char *mask = dset->padmask;
    double *Zi;
    char **S = NULL;
    int padded_n = dset->n; /* current series length */
    int real_n = dset->n;   /* target series length */
    int i, t;
    int err = 0;

    for (t=0; t<dset->n; t++) {
        if (mask[t]) {
            real_n--;
        }
    }

    fprintf(stderr, "undo_panel_padding: padded n*T = %d, original dset->n = %d\n",
            padded_n, real_n);

    if (real_n == padded_n) {
        fprintf(stderr, "strange, couldn't find any padding!\n");
        return E_DATA;
    }

    /* temporary holder for shorter series */
    Zi = malloc(real_n * sizeof *Zi);

    if (Zi == NULL) {
        return E_ALLOC;
    }

    if (dset->S != NULL) {
        /* holder for shorter obs labels array */
        S = strings_array_new_with_length(real_n, OBSLEN);
    }

    if (!err) {
        /* write non-padding rows from dset->Z (and dset->S, if present)
           into the right places in the reduced-size arrays */
        int s;

        for (i=0; i<dset->v; i++) {
            s = 0;
            for (t=0; t<padded_n; t++) {
                if (!mask[t]) {
                   Zi[s] = dset->Z[i][t];
                   if (i == 0 && S != NULL) {
                       strcpy(S[s], dset->S[t]);
                   }
                   s++;
                }
            }
            /* copy data back to dset->Z */
            memcpy(dset->Z[i], Zi, real_n * sizeof *Zi);
        }

        if (dset->S != NULL && S != NULL) {
            strings_array_free(dset->S, padded_n);
            dset->S = S;
        }
    }

    free(Zi);

    dset->n = real_n;
    dset->t2 = dset->n - 1;
    free(dset->padmask);
    dset->padmask = NULL;

    return err;
}

static int pad_panel_dataset (const double *uid, int uv, int nunits,
                              const double *tid, int tv, int nperiods,
                              DATASET *dset, int ustrs, char *mask)
{
    double **bigZ = NULL;
    char **S = NULL;
    int *nuid = NULL;
    int *ntid = NULL;
    int big_n, n_orig = dset->n;
    int i, s, t;
    int err = 0;

    big_n = nunits * nperiods;

    fprintf(stderr, "pad_panel_dataset: n*T = %d*%d = %d but dset->n = %d\n",
            nunits, nperiods, big_n, n_orig);

    err = normalize_uid_tid(tid, nperiods, dset, uv, tv,
                            &nuid, &ntid);
    if (err) {
        return err;
    }

    bigZ = doubles_array_new(dset->v, big_n);
    if (bigZ == NULL) {
        err = E_ALLOC;
    } else {
        dset->n = big_n;
        dset->t2 = dset->n - 1;
        for (i=0; i<dset->v; i++) {
            for (t=0; t<big_n; t++) {
                bigZ[i][t] = (i == 0)? 1.0 : NADBL;
            }
        }
    }

    if (!err && dset->S != NULL && ustrs) {
        S = strings_array_new_with_length(dset->n, OBSLEN);
    }

    if (!err) {
        int j, buv = 0, btv = 0;

        if (mask != NULL) {
            for (t=0; t<big_n; t++) {
                mask[t] = 1;
            }
        }

        /* write rows from original Z into the right places in bigZ */
        for (t=0; t<n_orig; t++) {
            s = nuid[t] * nperiods + ntid[t];
            for (i=1; i<dset->v; i++) {
                bigZ[i][s] = dset->Z[i][t];
            }
            if (S != NULL) {
                strcpy(S[s], dset->S[t]);
            }
            if (mask != NULL) {
                /* record (non-)padding in @mask */
                mask[s] = 0;
            }
        }

        /* where are uv and tv in bigZ? */
        j = 1;
        for (i=1; i<dset->v; i++) {
            if (i == uv) {
                buv = j;
            } else if (i == tv) {
                btv = j;
            }
            j++;
        }

        /* complete the index info in slots uv and tv */
        i = j = 0;
        for (t=0; t<dset->n; t++) {
            if (t > 0 && t % nperiods == 0) {
                i++;
                j = 0;
            }
            bigZ[buv][t] = uid[i];
            bigZ[btv][t] = tid[j++];
        }

        /* swap the padded arrays into Z */
        for (i=0; i<dset->v; i++) {
            free(dset->Z[i]);
            dset->Z[i] = bigZ[i];
        }

        free(bigZ);
    }

    free(nuid);
    free(ntid);

    if (!err && dset->S != NULL) {
        /* expand the obs (unit) marker strings appropriately */
        if (S == NULL) {
            strings_array_free(dset->S, n_orig);
            dset->S = NULL;
            dset->markers = NO_MARKERS;
        } else {
            char si[OBSLEN];
            int j;

            for (i=0; i<nunits; i++) {
                t = i * nperiods;
                for (j=0; j<nperiods; j++) {
                    if (S[t][0] != '\0') {
                        strcpy(si, S[t]);
                        break;
                    }
                    t++;
                }
                t = i * nperiods;
                for (j=0; j<nperiods; j++) {
                    strcpy(S[t++], si);
                }
            }
            strings_array_free(dset->S, n_orig);
            dset->S = S;
        }
    }

    return err;
}

static int check_index_values (const double *x, int n)
{
    int i;

    for (i=1; i<n; i++) {
        if (x[i] < 0) {
            return E_DATA;
        } else if (na(x[i])) {
            return E_MISSDATA;
        }
    }

    return 0;
}

static int uv_tv_from_varnames (const char *uname, const char *tname,
                                const DATASET *dset, int *uv, int *tv)
{
    *uv = current_series_index(dset, uname);
    if (*uv <= 0) {
        return E_DATA;
    }

    *tv = current_series_index(dset, tname);
    if (*tv <= 0) {
        return E_DATA;
    }

    return 0;
}

#if PDEBUG
static void print_unit_var (int uv, double **Z, int n, int after)
{
    int i, imax = 20;

    fprintf(stderr, "Z[uv], %s sorting (first %d obs):\n",
            (after)? "after" : "before", imax);
    for (i=0; i<n && i<imax; i++) {
        fprintf(stderr, " Z[%d][%d] = %g\n", uv, i, Z[uv][i]);
    }
}
#endif

static void finalize_panel_datainfo (DATASET *dset, int nperiods)
{
    int pdp = nperiods;
    int den = 10.0;

    while ((pdp = pdp / 10)) {
        den *= 10;
    }
    dset->structure = STACKED_TIME_SERIES;
    dset->pd = nperiods;
    dset->sd0 = 1.0 + 1.0 / den;
    ntolabel(dset->stobs, 0, dset);
    ntolabel(dset->endobs, dset->n - 1, dset);
}

static int check_full_dataset (void)
{
    DATASET *fset = fetch_full_dataset();

    if (!dataset_is_panel(fset)) {
        const char *msg =
            N_("You cannot use the --panel-vars option with the setobs command when\n"
            "\n"
            "* the dataset is currently sub-sampled and\n"
            "* the full dataset is not a panel.\n"
            "\n"
            "If you first structure your full dataset as a panel, you can then\n"
            "do what you are trying to do.");

        gretl_errmsg_set(_(msg));
        return E_DATA;
    }

    return 0;
}

/**
 * set_panel_structure_from_vars:
 * @uv: index of series uniquely identifying units/groups.
 * @tv: index of series uniquely identifying time periods.
 * @dset: pointer to dataset.
 *
 * Sets the panel structure of @dset based on the information
 * in the series with indices @uv and @tv, if possible.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int set_panel_structure_from_vars (int uv, int tv, DATASET *dset)
{
    int subsampled, presorted;
    double *uid = NULL;
    double *tid = NULL;
    char *mask = NULL;
    int n = dset->n;
    int totmiss = 0;
    int fulln = 0;
    int nunits = 0;
    int nperiods = 0;
    int ustrs = 0;
    int err = 0;

    subsampled = complex_subsampled();

    if (subsampled) {
        err = check_full_dataset();
        if (err) {
            return err;
        }
    }

#if PDEBUG
    fprintf(stderr, "\n*** set_panel_structure_from_vars:\n "
            "using var %d ('%s') for unit, var %d ('%s') for time\n",
            uv, dset->varname[uv], tv, dset->varname[tv]);
    print_unit_var(uv, dset->Z, n, 0);
#endif

    /* start by making sorted copies of the unit and time
       variables and counting the unique values of each
    */

    uid = copyvec(dset->Z[uv], n);
    tid = copyvec(dset->Z[tv], n);

    if (uid == NULL || tid == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    qsort(uid, n, sizeof *uid, gretl_compare_doubles);
    nunits = count_distinct_values(uid, dset->n);

    qsort(tid, n, sizeof *tid, gretl_compare_doubles);
    nperiods = count_distinct_values(tid, dset->n);

    /* check that we have a possible panel structure */

    if (nunits == 1 || nperiods == 1 ||
        nunits == n || nperiods == n ||
        n > nunits * nperiods) {
        fprintf(stderr, "Dataset does not have a panel structure\n");
        fprintf(stderr, " (nunits=%d, nperiods=%d, n=%d)\n", nunits, nperiods, n);
        err = E_DATA;
        goto bailout;
    }

    /* figure how many panel-data rows are implied, and how many
       "padding" rows are needed to complete the structure
    */

    fulln = nunits * nperiods;
    totmiss = fulln - dset->n;

#if PDEBUG
    fprintf(stderr, "Found %d units, %d periods\n", nunits, nperiods);
    fprintf(stderr, "Units: min %g, max %g\n", uid[0], uid[n - 1]);
    fprintf(stderr, "Periods: min %g, max %g\n", tid[0], tid[n - 1]);
    fprintf(stderr, "Required rows = %d * %d = %d\n", nunits, nperiods, fulln);
    fprintf(stderr, "Missing rows = %d - %d = %d\n", fulln, dset->n, totmiss);
#endif

    /* determine if the data rows are already in sort order by unit,
       and by period for each unit */

    presorted = panel_is_sorted(dset, uv, tv);

#if PDEBUG
    fprintf(stderr, "panel_is_sorted? %s\n", presorted ? "yes" : "no");
#endif

    if (!presorted) {
        if (subsampled) {
            /* We can't re-order the rows in a subsampled dataset, or the
               data will get scrambled when restoration of the full dataset
               is attempted.
            */
            gretl_errmsg_set(_("Sorry, can't do this with a sub-sampled dataset"));
            err = E_DATA;
        } else {
            /* sort full dataset by unit and period */
            err = panel_data_sort_by(dset, uv, tv, &ustrs);
        }
#if PDEBUG
        if (!err) print_unit_var(uv, dset->Z, n, 1);
#endif
    }

    if (!err && totmiss > 0) {
        /* establish a nominally balanced panel */
        mask = malloc(fulln);
        if (mask == NULL) {
            err = E_ALLOC;
        } else {
            rearrange_id_array(uid, nunits, n);
            rearrange_id_array(tid, nperiods, n);
            err = pad_panel_dataset(uid, uv, nunits,
                                    tid, tv, nperiods,
                                    dset, ustrs,
                                    mask);
        }
    }

    if (!err) {
        finalize_panel_datainfo(dset, nperiods);
    }

 bailout:

    free(uid);
    free(tid);

    if (!err && complex_subsampled()) {
        dset->padmask = mask;
    } else {
        free(mask);
    }

    return err;
}

int
set_panel_structure_from_varnames (const char *uname, const char *tname,
                                   DATASET *dset)
{
    int n = dset->n;
    int uv = 0, tv = 0;
    int err = 0;

    err = uv_tv_from_varnames(uname, tname, dset, &uv, &tv);

    if (!err) {
        err = check_index_values(dset->Z[uv], n);
    }

    if (!err) {
        err = check_index_values(dset->Z[tv], n);
    }

    if (!err) {
        err = set_panel_structure_from_vars(uv, tv, dset);
    }

    return err;
}

static int group_uniqueness_check (char **S, int n)
{
    int i, j;

    for (i=0; i<n-1; i++) {
        for (j=i+1; j<n; j++) {
            if (!strcmp(S[i], S[j])) {
                gretl_errmsg_sprintf(_("The string '%s' is given for "
                                     "two or more groups"), S[i]);
                return E_DATA;
            }
        }
    }

    return 0;
}

static char **group_names_from_array (const char *aname,
                                      int n_groups,
                                      int *err)
{
    gretl_array *A = get_array_by_name(aname);
    char **S = NULL;
    int ns = 0;

    if (A == NULL) {
        *err = E_DATA;
    } else {
        char **AS = gretl_array_get_strings(A, &ns);

        if (ns != n_groups) {
            *err = E_DATA;
        } else {
            S = strings_array_dup(AS, ns);
            if (S == NULL) {
                *err = E_ALLOC;
            }
        }
    }

    return S;
}

/* usable_groups_series: to be usable as the basis for panel
   groups strings, a series must have values that are (a)
   constant within-group and (b) unique across groups. Note
   that we only get to this test if the series in question
   has already been found to be string-valued.
*/

static int usable_groups_series (DATASET *dset, int v,
                                 const char *vname)
{
    const double *x = dset->Z[v];
    int i, j, g = 0;
    int ok = 1;

    if (na(x[0])) {
        return 0;
    }

    for (i=1; i<dset->n && ok; i++) {
        if (na(x[i])) {
            ok = 0;
        } else if (i % dset->pd == 0) {
            /* starting a new group: x-value must not repeat a
               prior group's value */
            for (j=0; j<=g; j++) {
                if (x[i] == x[j*dset->pd]) {
                    gretl_errmsg_sprintf(_("%s: values are not unique per group"),
                                         vname);
                    ok = 0;
                }
            }
            g++;
        } else if (x[i] != x[i-1]) {
            gretl_errmsg_sprintf(_("%s: is not constant within group"), vname);
            ok = 0;
        }
    }

    return ok;
}

static int maybe_use_strval_series (DATASET *dset,
                                    const char *vname,
                                    int ng)
{
    int v = current_series_index(dset, vname);
    series_table *st = NULL;
    int ns = 0;
    int err = 0;

    if (v <= 0) {
        return E_INVARG;
    }

    st = series_get_string_table(dset, v);
    if (st == NULL) {
        gretl_errmsg_sprintf(_("The series %s is not string-valued"), vname);
        return E_INVARG;
    }

    series_table_get_strings(st, &ns);

#if 0
    fprintf(stderr, "maybe_use_strval_series (%s, ID %d, ns=%d)\n", vname, v, ns);
#endif

    if (ns < ng) {
        gretl_errmsg_sprintf(_("The series %s holds %d strings but %d "
                             "are needed"), vname, ns, ng);
        err = E_INVARG;
    } else {
        /* note: don't mess with numerical values, just check them */
        if (usable_groups_series(dset, v, vname)) {
            set_panel_groups_name(dset, vname);
        } else {
            err = E_INVARG;
        }
    }

    return err;
}

/* Enable construction of a string-valued series holding a name for
   each panel group/unit. The name of this series is set on the
   'pangrps' member of @dset, and the names are then used in panel
   graphs where appropriate.

   @vname should contain the name of a series, either an existing one
   (which may be overwritten) or a new one to create; and @grpnames
   should be (1) a string literal, or (2) the name of a string variable
   (which in either case should hold N space-separated strings, where N
   is the number of panel groups), or (3) the name of an array variable
   holding N strings.
*/

int set_panel_group_strings (const char *vname,
                             const char *grpnames,
                             DATASET *dset)
{
    const char *namestr = NULL;
    char **S = NULL;
    int ng = dset->n / dset->pd;
    int v, orig_v = dset->v;
    int err = 0;

    if (vname == NULL || *vname == '\0') {
        return E_DATA;
    }

    if (grpnames == NULL || *grpnames == '\0') {
        /* We just got the name of a series? That may be
           OK if it's string-valued and has a suitable
           set of values.
        */
        return maybe_use_strval_series(dset, vname, ng);
    }

    if (strchr(grpnames, ' ')) {
        /* group names as string literal */
        namestr = grpnames;
    } else if (gretl_is_string(grpnames)) {
        /* group names in a single string variable */
        namestr = get_string_by_name(grpnames);
    } else {
        /* try for array of strings */
        S = group_names_from_array(grpnames, ng, &err);
    }

    if (namestr != NULL) {
        /* we must obtain @S by splitting */
        int ns = 0;

        if (strchr(namestr, '"') != NULL) {
            S = gretl_string_split_quoted(namestr, &ns, NULL, &err);
        } else {
            S = gretl_string_split(namestr, &ns, " \n\t");
        }
        if (!err && S == NULL) {
            err = E_ALLOC;
        }
        if (!err && ns != ng) {
            /* FIXME subsampled case? */
            fprintf(stderr, "Got %d strings but there are %d groups\n",
                    ns, ng);
            strings_array_free(S, ns);
            S = NULL;
            err = E_DATA;
        }
    }

    if (!err && S != NULL) {
        err = group_uniqueness_check(S, ng);
    }

    if (!err) {
        v = current_series_index(dset, vname);
        if (v < 0) {
            /* we need to add a series */
            char *gen = gretl_strdup_printf("%s", vname);

            err = generate(gen, dset, GRETL_TYPE_SERIES, OPT_Q, NULL);
            if (!err) {
                v = dset->v - 1;
            }
            free(gen);
        }
    }

    if (!err) {
        series_table *st = series_table_new(S, ng, &err);

        if (!err) {
            int i, g = 0;

            series_attach_string_table(dset, v, st);
            for (i=0; i<dset->n; i++) {
                if (i % dset->pd == 0) {
                    g++;
                }
                dset->Z[v][i] = g;
            }
        }
    }

    if (err) {
        if (S != NULL) {
            strings_array_free(S, ng);
        }
        if (dset->v > orig_v) {
            dataset_drop_last_variables(dset, dset->v - orig_v);
        }
    } else {
        set_panel_groups_name(dset, vname);
    }

    return err;
}

/* utility functions */

static int find_time_var (const DATASET *dset)
{
    const char *tnames[] = {
        "year",
        "Year",
        "period",
        "Period",
	"time",
        NULL
    };
    int i, v;

    for (i=0; tnames[i] != NULL; i++) {
        v = series_index(dset, tnames[i]);
        if (v < dset->v) {
            return v;
        }
    }

    return 0;
}

int guess_panel_structure (double **Z, DATASET *dset)
{
    int ret, v = find_time_var(dset);

    if (v == 0) {
        ret = 0; /* can't guess */
    } else if (floateq(Z[v][0], Z[v][1])) {
        /* "year" or whatever is same for first two obs */
        dset->structure = STACKED_CROSS_SECTION;
        ret = STACKED_CROSS_SECTION;
    } else {
        dset->structure = STACKED_TIME_SERIES;
        ret = STACKED_TIME_SERIES;
    }

    return ret;
}

int balanced_panel (const DATASET *dset)
{
    int n = dset->t2 - dset->t1 + 1;
    int ret = 0;

    if (n % dset->pd == 0) {
        char unit[OBSLEN], period[OBSLEN];

        if (sscanf(dset->endobs, "%[^:]:%s", unit, period) == 2) {
            if (atoi(period) == dset->pd) {
                ret = 1;
            }
        }
    }

    return ret;
}

static int may_be_time_name (const char *vname)
{
    char test[VNAMELEN];

    strcpy(test, vname);
    gretl_lower(test);

    return strcmp(test, "year") == 0 ||
        strcmp(test, "time") == 0 ||
	strcmp(test, "period") == 0;
}

/* See if the panel dataset contains a variable that plausibly
   represents the year of the observations. We look for series labeled
   "year", "time" or "period" (in a case insensitive comparison) and,
   if found, check that the series has the same sequence of values for
   each individual, and the same increment between successive
   time-series observations.  (The increment does not have to be 1, to
   accommodate, e.g., quinquennial or decennial observations.)

   Return the ID number of a suitable variable, or 0 if there's nothing
   suitable.

   This may be used in setting the x-axis tic marks in a panel time-
   series plot.
*/

int plausible_panel_time_var (const DATASET *dset)
{
    int i, t, ret = 0;

    for (i=1; i<dset->v && ret==0; i++) {
        if (may_be_time_name(dset->varname[i])) {
            const double *x = dset->Z[i];
            int val0 = x[0];
            int incr0 = x[1] - x[0];
            int ok = 1;

            for (t=0; t<dset->n && ok; t++) {
                if (na(x[t]) || x[t] < 0) {
                    ok = 0;
                } else if (t > 0 && t % dset->pd == 0) {
                    if (x[t] != val0) {
                        ok = 0;
                    }
                } else if (t > 1 && x[t] - x[t-1] != incr0) {
                    ok = 0;
                }
            }
            if (ok) {
                ret = i;
            }
        }
    }

    return ret;
}

/* Tests whether the series with index number @v codes for the
   time dimension of a panel dataset. If so, returns 1 and
   records the minimum value and (constant) increment of the
   index; otherwise returns 0.
*/

int is_panel_time_var (const DATASET *dset, int v,
                       int tmax, int *minval,
                       int *incr)
{
    int t, ret = 0;

    *minval = 0;
    *incr = 0;

    if (may_be_time_name(dset->varname[v])) {
        const double *x = dset->Z[v];
        int val0 = (int) x[0];
        int incr0 = (int) x[1] - (int) x[0];
        int ok = 1;

        for (t=0; t<tmax && ok; t++) {
            if (na(x[t]) || x[t] < 0 || x[t] != floor(x[t])) {
                ok = 0;
            } else if (t > 0 && t % dset->pd == 0) {
                if (x[t] != val0) {
                    ok = 0;
                }
            } else if (t > 1 && x[t] - x[t-1] != incr0) {
                ok = 0;
            }
        }
        if (ok) {
            *minval = val0;
            *incr = incr0;
            ret = 1;
        }
    }

    return ret;
}

int is_panel_unit_var (const DATASET *dset, int v, int tmax)
{
    const double *x = dset->Z[v];
    int t, ret = 1;

    if (x[0] != 1) {
        return 0;
    }

    for (t=1; t<tmax && ret; t++) {
        if (na(x[t]) || x[t] < 0 || x[t] != floor(x[t])) {
            ret = 0;
        } else if (t % dset->pd == 0) {
            if (x[t] != x[t-1] + 1) {
                ret = 0;
            }
        } else if (x[t] != x[t-1]) {
            ret = 0;
        }
    }

    return ret;
}

/* FIXME: this does not yet handle the dropping of instruments */

static int *dpanel_list_omit (const MODEL *orig, const int *drop, int *err)
{
    const int *old = orig->list;
    int *new = gretl_list_copy(old);
    int sep = 0;
    int i, j;

    if (new == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (i=2; i<=new[0]; i++) {
        if (new[i] == LISTSEP) {
            sep++;
        }
        if (sep == 1) {
            for (j=1; j<=drop[0]; j++) {
                if (drop[j] == new[i]) {
                    gretl_list_delete_at_pos(new, i--);
                }
            }
        }
    }

#if 0
    printlist(old, "old");
    printlist(drop, "drop");
    printlist(new, "new");
#endif

    return new;
}

/**
 * panel_list_omit:
 * @orig: list specifying original panel model.
 * @drop: list of variables to be omitted.
 * @err: pointer to receive error code.
 *
 * Creates a new panel regression list, by first reconstructing
 * the regression specification from @orig then deleting from
 * the reconstruction the variables found in @drop.
 *
 * Returns: the new, reduced list or %NULL on error.
 */

int *panel_list_omit (const MODEL *orig, const int *drop, int *err)
{
    int *newlist = NULL;
    int i;

    if (orig->ci == DPANEL) {
        return dpanel_list_omit(orig, drop, err);
    }

    /* sorry, can't drop the constant */
    if (drop != NULL) {
        int cpos = in_gretl_list(drop, 0);

        if (cpos >= 2) {
            gretl_errmsg_set(_("Panel models must include an intercept"));
            *err = E_DATA;
            return NULL;
        }
    }

    if (orig->opt & OPT_F) {
        int *panlist;

        /* fixed-effects lists have the constant removed,
           so we need to put it back first
        */
        panlist = gretl_list_new(orig->list[0] + 1);
        if (panlist != NULL) {
            panlist[1] = orig->list[1];
            panlist[2] = 0;
            for (i=3; i<=panlist[0]; i++) {
                panlist[i] = orig->list[i-1];
            }
            if (drop == NULL) {
                newlist = gretl_list_omit_last(panlist, err);
            } else {
                newlist = gretl_list_omit(panlist, drop, 2, err);
            }
            free(panlist);
        }
    } else if (drop == NULL) {
        newlist = gretl_list_omit_last(orig->list, err);
    } else {
        newlist = gretl_list_omit(orig->list, drop, 2, err);
    }

    return newlist;
}

/* FIXME doesn't handle adding instruments */

static int *dpanel_list_add (const MODEL *orig, const int *add, int *err)
{
    const int *old = orig->list;
    int *new = gretl_list_copy(old);
    int sep = 0, pos = old[0] + 1;
    int i;

    if (new == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (i=2; i<=old[0]; i++) {
        if (old[i] == LISTSEP) {
            sep++;
            if (sep == 2) {
                pos = i - 1;
            }
        }
    }

    gretl_list_insert_list(&new, add, pos);
    if (new == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

#if 0
    printlist(old, "old");
    printlist(add, "add");
    printlist(new, "new");
#endif

    return new;
}

/**
 * panel_list_add:
 * @orig: list specifying original panel model.
 * @add: list of variables to be added.
 * @err: pointer to receive error code.
 *
 * Creates a new panel regression list, by first reconstructing
 * the regression specification from @orig then adding to
 * the reconstruction the variables found in @add.
 *
 * Returns: the new, augmented list or %NULL on error.
 */

int *panel_list_add (const MODEL *orig, const int *add, int *err)
{
    if (orig->ci == DPANEL) {
        return dpanel_list_add(orig, add, err);
    } else {
        return gretl_list_add(orig->list, add, err);
    }
}

/* calculate the within and between standard deviations for a given
   variable in a panel data set */

int panel_variance_info (const double *x, const DATASET *dset,
                         double xbar, double *psw, double *psb)
{
    double sw = 0.0, sb = 0.0;
    double xibar, d;
    int effn, effnT;
    int n, T, nT, Ti;
    int i, t, s;

    if (!dataset_is_panel(dset)) {
        return E_PDWRONG;
    }

    nT = dset->t2 - dset->t1 + 1;
    T = dset->pd;
    n = nT / T;

    effn = 0;
    effnT = 0;

    for (i=0; i<n; i++) {
        Ti = 0;
        xibar = 0.0;
        for (t=0; t<T; t++) {
            s = dset->t1 + i * T + t;
            if (!na(x[s])) {
                Ti++;
                xibar += x[s];
            }
        }
        if (Ti > 1) {
            xibar /= Ti;
            for (t=0; t<T; t++) {
                s = dset->t1 + i * T + t;
                if (!na(x[s])) {
                    d = x[s] - xibar;
                    sw += d * d;
                }
            }
        }
        if (Ti > 0) {
            /* is this right for singleton observations? */
            d = xibar - xbar;
            sb += d * d;
            effn++;
            effnT += Ti;
        }
    }

    if (effn > 1) {
        /* between variance: \sum_{i=1}^N (\bar{x}_i - \bar{x})^2
           \times (N-1)^{-1}
        */
        sb /= (effn - 1);
        sb = sqrt(sb);
    } else {
        sb = NADBL;
    }

    if (effnT - effn > 0) {
        /* within variance: \sum_{i=1}^N \sum_{t=1}^T (x_{it} - \bar{x}_i)^2
           \times (\sum_{i=1}^N T_i - N)^{-1}
        */
        sw /= (effnT - effn);
        sw = sqrt(sw);
    } else {
        sw = NADBL;
    }

    *psw = sw;
    *psb = sb;

    return 0;
}

int series_is_group_invariant (const DATASET *dset, int v)
{
    const double *x = dset->Z[v];
    int T = dset->pd;
    int N = dset->n / T;
    int i, t, s;

    for (i=1; i<N; i++) {
        for (t=0; t<T; t++) {
            s = t + i * T;
            if (x[s] != x[t]) {
                return 0;
            }
        }
    }

    return 1;
}

int panel_padding_rows (const DATASET *dset)
{
    int missrow, nmiss = 0;
    int i, t;

    for (t=dset->t1; t<=dset->t2; t++) {
        missrow = 1;
        for (i=1; i<dset->v; i++) {
            if (!na(dset->Z[i][t])) {
                missrow = 0;
                break;
            }
        }
        if (missrow) {
            nmiss++;
        }
    }

    return nmiss;
}

/* Given auxiliary time-series info attached to the panel
   dataset @pdset, transcribe it to the time-series dataset
   @tsset.
*/

int time_series_from_panel (DATASET *tset, const DATASET *pset)
{
    int err = 0;

    if (pset->panel_pd == 0) {
        /* the panel time dimension is not set */
        tset->structure = SPECIAL_TIME_SERIES;
        tset->pd = 1;
        return 0;
    }

    tset->structure = TIME_SERIES;
    tset->pd = pset->panel_pd;
    tset->sd0 = pset->panel_sd0;

    if (tset->pd == 1) {
        sprintf(tset->stobs, "%d", (int) tset->sd0);
    } else if (tset->pd == 4 || tset->pd == 12) {
        double dyr = floor(tset->sd0);
        double dp = tset->sd0 - dyr;
        int yr, p;

        yr = (int) dyr;
        dp *= (tset->pd == 4)? 10 : 100;
        p = nearbyint(dp);
        if (yr > 0 && yr < 9999 && p > 0) {
            if (tset->pd == 4 && p < 4) {
                sprintf(tset->stobs, "%d:%d", yr, p);
            } else if (p < 12) {
                sprintf(tset->stobs, "%d:%02d", yr, p);
            }
        } else {
            err = 1;
        }
    } else if (calendar_data(tset)) {
        calendar_date_string(tset->stobs, 0, tset);
    }

    return err;
}

/* Given an annual, quarterly or monthly date string @s in gretl's
   standard format (YYYY, YYYY:Q or YYYY:MM), convert to a zero-based
   observation index relative to a series of frequency @pd and length
   @n starting at @t0.

   @t0 is assumed to be given in a form that makes the calculation
   easy, namely the number of periods since the start of the year 0.
   For annual data this is just the year; for quarterly its 4 * year
   plus initial quarter and for monthly 12 * year plus initial month.

   If the resulting index is out of bounds, -1 is returned.
*/

int obs_index_from_aqm (const char *s, int pd, int t0, int n)
{
    const char *digits = "0123456789";
    int ok_len = (pd == 1)? 4 : (pd == 4)? 6 : 7;
    int len = strlen(s);
    int t = 0;

    if (len != ok_len || strspn(s, digits) != 4) {
        return -1;
    }

    if (pd == 1) {
        t = atoi(s);
    } else if (s[4] != ':' || strspn(s + 5, digits) != len - 5) {
        t = -1;
    } else {
        int sub = atoi(s + 5);

        if (sub > pd) {
            t = -1;
        } else {
            t = pd * atoi(s) + sub;
        }
    }

    if (t >= 0) {
        t -= t0;
        if (t < 0 || t >= n) {
            t = -1;
        }
    }

    return t;
}

/* Scan @dset for a string-valued series that's suitable for defining
   panel group names: its number of string values must equal the
   number of groups, and its value must be constant within each
   group's time series. In addition, if @maxlen > 0 the strings must
   be no longer than @maxlen UTF-8 characters (we're looking to use
   these strings for plotting and we don't want them too long).

   We return the ID number of the first series to match these
   criteria, or 0 if we can't find any such series.
*/

int usable_group_names_series_id (const DATASET *dset, int maxlen)
{
    char **S;
    int N = dset->n / dset->pd;
    int j, t, s, ns, ok;
    int i, ret = 0;

    for (i=1; i<dset->v; i++) {
        if (!is_string_valued(dset, i)) {
            continue;
        }
        S = series_get_string_vals(dset, i, &ns, 0);
        if (ns < N) {
            continue;
        }
        ok = 1;
        for (j=0, s=0; j<N && ok; j++) {
            if (maxlen > 0 && g_utf8_strlen(S[j], -1) > maxlen) {
                /* too long to be usable for plots? */
                ok = 0;
            }
            for (t=0; t<dset->pd && ok; t++) {
                if (t > 0 && dset->Z[i][s] != dset->Z[i][s-1]) {
                    ok = 0;
                } else if (j > 0 && t == 0 &&
                           dset->Z[i][s] == dset->Z[i][s-1]) {
                    ok = 0;
                }
                s++;
            }
        }
        if (ok) {
            ret = i;
            break;
        }
    }

    return ret;
}
