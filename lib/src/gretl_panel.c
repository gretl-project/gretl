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
    int Tmax;             /* effective times-series length (max usable obs per unit) */
    int Tmin;             /* shortest usable times-series */
    double Tbar;          /* harmonic mean of per-units time-series lengths */
    int NT;               /* total observations used (based on pooled model) */
    int ntdum;            /* number of time dummies added */
    int *unit_obs;        /* array of number of observations per x-sect unit */
    char *varying;        /* array to record properties of pooled-model regressors */
    int *vlist;           /* list of time-varying variables from pooled model */
    int balanced;         /* 1 if the model dataset is balanced, else 0 */
    int nbeta;            /* number of slope coeffs for Hausman test */
    int Fdfn;             /* numerator df, F for differing intercepts */
    int Fdfd;             /* denominator df, F for differing intercepts */
    double theta;         /* quasi-demeaning coefficient */
    double F;             /* joint significance of differing unit intercepts */
    double BP;            /* Breusch-Pagan test statistic */
    double H;             /* Hausman test statistic */
    gretl_matrix *bdiff;  /* array of coefficient differences */
    gretl_matrix *Sigma;  /* Hausman covariance matrix */
    double s2b;           /* "between" error variance */
    double s2e;           /* \hat{sigma}^2_e, from fixed-effects estimator */
    double s2v;           /* \hat{sigma}^2_v measure */
    int *small2big;       /* data indexation array */
    int *big2small;       /* reverse data indexation array */
    MODEL *pooled;        /* reference model (pooled OLS) */
    MODEL *realmod;       /* fixed or random effects model */
    unit_ts *tsinfo;      /* per unit time-series info */
};

struct {
    int n;      /* number of cross-sectional units */
    int T;      /* number of observations per unit */
    int offset; /* sampling offset into full dataset */
} panidx;

static int varying_vars_list (const DATASET *dset, panelmod_t *pan);

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
    pan->nbeta = 0;
    pan->theta = NADBL;

    pan->F = NADBL;
    pan->Fdfn = 0;
    pan->Fdfd = 0;

    pan->BP = NADBL;
    pan->H = NADBL;

    pan->bdiff = NULL;
    pan->Sigma = NULL;

    pan->small2big = NULL;
    pan->big2small = NULL;
    
    pan->pooled = NULL;
    pan->realmod = NULL;
    pan->tsinfo = NULL;
}

static void panelmod_free_tsinfo (panelmod_t *pan)
{
    if (pan->tsinfo != NULL) {
	free(pan->tsinfo->t0);
	free(pan->tsinfo->T);
	free(pan->tsinfo);
	pan->tsinfo = NULL;
    }
}

static int panelmod_add_tsinfo (panelmod_t *pan)
{
    int err = 0;

    pan->tsinfo = malloc(sizeof *pan->tsinfo);
    if (pan->tsinfo == NULL) {
	return E_ALLOC;
    }

    pan->tsinfo->t0 = malloc(pan->nunits * sizeof *pan->tsinfo->t0);
    pan->tsinfo->T = malloc(pan->nunits * sizeof *pan->tsinfo->T);

    if (pan->tsinfo->t0 == NULL || pan->tsinfo->T == NULL) {
	panelmod_free_tsinfo(pan);
	err = E_ALLOC;
    }

    return err;
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

    free(pan->realmod);

    panelmod_free_tsinfo(pan);
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
beck_katz_vcv (MODEL *pmod, panelmod_t *pan, const double **Z,
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
		    gretl_matrix_set(Xi, s, p, Z[v][si]);
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
			gretl_matrix_set(Xi, s, p, Z[v][si]);
			gretl_matrix_set(Xj, s, p, Z[v][sj]);
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

/* HAC covariance matrix for the pooled or fixed-effects model, given
   "fixed T and large N".  In the case of "large T" a different form
   is needed for robustness in respect of autocorrelation.  See
   Arellano, "Panel Data Econometrics" (Oxford, 2003), pages 18-19.

   In the case of fixed effects there should be no missing values,
   since we created a special "within" dataset purged of same.
   But with the pooled model we need to index into the full dataset,
   and work around any missing values.
*/

static int 
arellano_vcv (MODEL *pmod, panelmod_t *pan, const double **Z,
	      gretl_matrix *XX, gretl_matrix *W, gretl_matrix *V)
{
    gretl_vector *e = NULL;
    gretl_matrix *Xi = NULL;
    gretl_vector *eXi = NULL;
    int T = pan->Tmax;
    int k = pmod->ncoeff;
    int i, j, v, s, t;
    int err = 0;

    e   = gretl_column_vector_alloc(T);
    Xi  = gretl_matrix_alloc(T, k);
    eXi = gretl_vector_alloc(k);

    if (e == NULL || Xi == NULL || eXi == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    s = 0;

    for (i=0; i<pan->nunits; i++) {
	int Ti = pan->unit_obs[i];
	int p = 0;

	if (Ti == 0) {
	    continue;
	}

	e = gretl_matrix_reuse(e, Ti, 1);
	Xi = gretl_matrix_reuse(Xi, Ti, k);

	for (t=0; t<pan->T; t++) {
	    s = panel_index(i, t);
	    if (na(pmod->uhat[s])) {
		continue;
	    }
	    gretl_vector_set(e, p, pmod->uhat[s]);
	    s = small_index(pan, s);
	    for (j=0; j<k; j++) {
		v = pmod->list[j+2];
		gretl_matrix_set(Xi, p, j, Z[v][s]);
	    }
	    p++;
	}

	gretl_matrix_multiply_mod(e, GRETL_MOD_TRANSPOSE,
				  Xi, GRETL_MOD_NONE,
				  eXi, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(eXi, GRETL_MOD_TRANSPOSE,
				  eXi, GRETL_MOD_NONE,
				  W, GRETL_MOD_CUMULATE);
    }

    /* form V(b_W) = (X'X)^{-1} W (X'X)^{-1} */
    gretl_matrix_qform(XX, GRETL_MOD_NONE, W,
		       V, GRETL_MOD_NONE);

#if 0
    gretl_matrix_print(XX, "X'X^{-1}");
    gretl_matrix_print(W, "W");
    gretl_matrix_print(V, "V");
#endif

    gretl_model_set_vcv_info(pmod, VCV_PANEL, PANEL_HAC);

 bailout:

    gretl_matrix_free(e);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

/* common setup for Arellano and Beck-Katz VCV estimators
   (note that pooled OLS with the --cluster option is
   handled separately)
*/

static int 
panel_robust_vcv (MODEL *pmod, panelmod_t *pan, const double **Z)
{
    gretl_matrix *W = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *XX = NULL;
    int k = pmod->ncoeff;
    int err = 0;

    W = gretl_zero_matrix_new(k, k);
    V = gretl_matrix_alloc(k, k);

    if (W == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    XX = panel_model_xpxinv(pmod, &err);
    if (err) {
	fprintf(stderr, "panel_robust_vcv: failed at panel_model_xpxinv\n");
	goto bailout;
    }

    /* call the appropriate function */
    if (libset_get_bool(PCSE)) {
	err = beck_katz_vcv(pmod, pan, Z, XX, W, V);
    } else {
	err = arellano_vcv(pmod, pan, Z, XX, W, V);
    }

    if (!err) {
	int i, j, s = 0;

	for (i=0; i<k; i++) {
	    for (j=i; j<k; j++) {
		pmod->vcv[s++] = gretl_matrix_get(V, i, j);
	    }
	    pmod->sderr[i] = sqrt(gretl_matrix_get(V, i, i));
	}
    }

 bailout:

    gretl_matrix_free(W);
    gretl_matrix_free(V);
    gretl_matrix_free(XX);

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

/* Determine and record the starting point and length of the longest
   unbroken sequence of time-series observations for each panel unit.
   We use this in computing rho and Durbin-Watson.
*/

static int 
panel_unit_AR_info (const MODEL *pmod, panelmod_t *pan)
{
    int i, t, s, Ti, Timax, t0, t0max, rem;
    int err, okunits = 0;

    err = panelmod_add_tsinfo(pan);
    if (err) {
	return err;
    }

    for (i=0; i<pan->nunits; i++) {
	pan->tsinfo->t0[i] = 0;
	pan->tsinfo->T[i] = 0;

	if (pan->unit_obs[i] < 2) {
	    continue;
	}

	Timax = t0max = t0 = 0;

	while (1) {
	    Ti = 0;
	    for (t=t0; t<pan->T; t++) {
		s = panel_index(i, t);
		if (na(pmod->uhat[s])) {
		    break;
		}
		Ti++; /* length of current non-NA sequence */
	    }
	    if (Ti > Timax) {
		Timax = Ti;
		t0max = t0;
	    }
	    rem = pan->T - t;
	    if (rem <= Timax) {
		break;
	    } else {
		t0 = t + 1;
	    }
	}

	if (Timax >= 2) {
	    pan->tsinfo->t0[i] = t0max;
	    pan->tsinfo->T[i] = Timax;
	    okunits++;
	} 
    }

    if (okunits == 0) {
	/* nothing doing */
	panelmod_free_tsinfo(pan);
	err = E_DATA;
    }

    return err;
}

static int panel_DW_pvalue (MODEL *pmod, const panelmod_t *pan,
			    const double **Z)
{
    gretl_matrix *X = NULL;
    gretl_matrix *XX = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *A = NULL;
    gretl_matrix *MA = NULL;
    gretl_matrix *Ai = NULL;
    gretl_matrix *E = NULL;
    double pv, sz;
    int T = 0;
    int k = pmod->ncoeff;
    int jmin = 0, effj = 0;
    int bigt = 0, ins = 0;
    int i, j, jj, t, s, v;
    int err = 0;

    /* determine the total observations to be used */
    for (i=0; i<pan->nunits; i++) {
	T += pan->tsinfo->T[i];
    }

    if (pan->opt & OPT_F) {
	/* fixed effects */
	k += pan->effn - 1; /* allow for unit dummies */
	jmin = 1;           /* offset to avoid average constant */
    }

    /* allocation size? */
    sz = T * k + 3 * T * T + k * k;
    sz *= 8.0 / (1024 * 1024);

    if (sz > 50) {
	/* more than 50 MB of storage: we'll not try this */
	fprintf(stderr, "panel_DW_pvalue: size = %g MB\n", sz);
    }

    X = gretl_matrix_alloc(T, k);
    M = gretl_identity_matrix_new(T);
    A = gretl_zero_matrix_new(T, T);
    MA = gretl_zero_matrix_new(T, T);
    XX = gretl_matrix_alloc(k, k);

    if (X == NULL || M == NULL || A == NULL || MA == NULL ||
	XX == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* note: we use only contiguous time-series observations,
       as recorded in pan->tsinfo */

    for (i=0; i<pan->nunits; i++) {
	int Ti = pan->tsinfo->T[i];
	int tmin, tmax;

	if (Ti == 0) {
	    continue;
	}

	Ai = gretl_DW_matrix_new(Ti);
	if (Ai == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	err = gretl_matrix_inscribe_matrix(A, Ai, ins, ins,
					   GRETL_MOD_NONE);
	if (err) {
	    break;
	}

	ins += Ti;
	gretl_matrix_free(Ai);

	tmin = pan->tsinfo->t0[i];
	tmax = tmin + Ti - 1;

	for (t=tmin; t<=tmax; t++) {
	    s = panel_index(i, t);
	    for (j=jmin; j<pmod->ncoeff; j++) {
		/* write regressors at t into X */
		v = pmod->list[j+2];
		gretl_matrix_set(X, bigt, j - jmin, Z[v][s]);
	    }
	    if (jmin > 0) {
		/* fixed effects: append unit dummies */
		for (j=0; j<pan->effn; j++) {
		    jj = j + pmod->ncoeff - 1;
		    gretl_matrix_set(X, bigt, jj, (effj == j)? 1.0 : 0.0);
		}
	    }
	    bigt++;
	}
	effj++;
    }

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
			      X, GRETL_MOD_NONE,
			      XX, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(XX);

    if (!err) {
	/* M = I - X(X'X)^{-1}X' */
	err = gretl_matrix_qform(X, GRETL_MOD_NONE,
				 XX, M, GRETL_MOD_DECREMENT);
    }

    if (!err) {
	err = gretl_matrix_multiply(M, A, MA);
    }

    if (!err) {
	E = gretl_general_matrix_eigenvals(MA, 0, &err);
    }

    if (!err) {
	k = T - k; /* number of non-zero eigenvalues */
	for (i=0; i<k; i++) {
	    E->val[i] -= pmod->dw;
	}
	gretl_matrix_reuse(E, k, 1);
	pv = imhof(E, 0.0, &err);
	if (!err) {
	    gretl_model_set_double(pmod, "dw_pval", pv);
	}
    }   

 bailout:

    gretl_matrix_free(XX);
    gretl_matrix_free(X);
    gretl_matrix_free(M);
    gretl_matrix_free(A);
    gretl_matrix_free(MA);
    gretl_matrix_free(E);

    return err;
}

/* Durbin-Watson statistic for the pooled or fixed effects model.  We
   only use units that have at least two consecutive time-series
   observations, and we use only consecutive observations.
   
   See Bhargava, Franzini and Narendranathan, "Serial Correlation and
   the Fixed Effects Model", Review of Economic Studies 49, 1982,
   pp. 533-549.
*/

static void panel_dwstat (MODEL *pmod, panelmod_t *pan)
{
    double ut, u1;
    double dwnum = 0.0, dwden = 0.0;
    double rnum = 0.0, rden = 0.0;
    int i, t, ti;

    pmod->dw = pmod->rho = NADBL;

    if (pmod->ess <= 0.0) {
	return;
    }

    if (panel_unit_AR_info(pmod, pan) != 0) {
	/* can't do this */
	return;
    }

    for (i=0; i<pan->nunits; i++) {
	int tmin, tmax;

	if (pan->tsinfo->T[i] == 2) {
	    continue;
	}

	tmin = pan->tsinfo->t0[i] + 1;
	tmax = pan->tsinfo->t0[i] + pan->tsinfo->T[i] - 1;

	for (ti=tmin; ti<=tmax; ti++) {
	    t = panel_index(i, ti);
	    ut = pmod->uhat[t];
	    u1 = pmod->uhat[t-1];
	    dwnum += (ut - u1) * (ut - u1);
	    dwden += ut * ut;
	    rnum += ut * u1;
	    rden += u1 * u1;
	    if (ti == tmin) {
		/* include first \hat{u}_t squared */
		dwden += u1 * u1;
	    }
	}
    }

    if (dwden > 0.0 && !na(dwden)) {
	pmod->dw = dwnum / dwden;
    }

    if (rden > 0.0 && !na(rden)) {
	pmod->rho = rnum / rden;
	if (pmod->rho <= -1.0 || pmod->rho >= 1.0) {
	    pmod->rho = NADBL;
	}
    } 
}

/* Allocate the arrays needed to perform the Hausman test,
   in its matrix formulation.
*/

static int hausman_allocate (panelmod_t *pan)
{
    int k = pan->vlist[0] - 2;

    pan->nbeta = k;

    if (pan->opt & OPT_M) {
	/* array to hold differences between coefficient estimates */
	pan->bdiff = gretl_vector_alloc(k);
	if (pan->bdiff == NULL) {
	    return E_ALLOC;
	}
	/* array to hold covariance matrix */
	pan->Sigma = gretl_matrix_alloc(k, k);
	if (pan->Sigma == NULL) {
	    gretl_matrix_free(pan->bdiff);
	    pan->bdiff = NULL;
	    return E_ALLOC; 
	}
    }

    return 0;
}   

/* printing routines for the "panel diagnostics" test */

static void print_panel_coeff (const MODEL *pmod,
 			       const char *vname,
 			       int i, int maxlen,
			       PRN *prn)
{
    double tstat = pmod->coeff[i] / pmod->sderr[i];
    char errstr[18];
    char pvstr[18];
 
    sprintf(errstr, "(%.5g)", pmod->sderr[i]);
    sprintf(pvstr, "[%.5f]", student_pvalue_2(pmod->dfd, tstat));
    pprintf(prn, "%*s: %14.5g %15s %15s\n", maxlen, vname,
 	    pmod->coeff[i], errstr, pvstr);
}

static void print_re_model_top (const panelmod_t *pan, PRN *prn)
{
    pputs(prn, "Variance estimators:\n");
    pprintf(prn, " between = %g\n", pan->s2b);
    pprintf(prn, " within = %g\n", pan->s2e);

    if (pan->balanced) {
	pprintf(prn, "theta used for quasi-demeaning = %g\n", pan->theta);
    } else {
	pputs(prn, "Panel is unbalanced: theta varies across units\n");
    }
    pputc(prn, '\n');

    pputs(prn,
	  _("                         Random effects estimator\n"
	    "           allows for a unit-specific component to the "
	    "error term\n"
	    "           (standard errors in parentheses, p-values in brackets)\n\n"));

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

	vj = vlist[j];
	gxbar = 0.0;
	s = 0;

#if PDEBUG
	strcpy(wset->varname[j], dset->varname[vj]);
#endif

	for (i=0; i<pan->nunits; i++) {
	    int Ti = pan->unit_obs[i];
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
		    got++;
		    if (pan->small2big != NULL) {
			pan->small2big[s] = bigt;
			pan->big2small[bigt] = s;
		    }
		    s++;
		}
	    }
	}

	/* wZ = data - group mean + grand mean */
	gxbar /= pan->NT;
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
random_effects_dataset (const DATASET *dset, const DATASET *gset,
			int *relist, int *hlist, 
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
	/* apparatus for regression version of Hausman test */
	for (i=1; i<pan->vlist[0]; i++) {
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
       version of Hausman test if wanted */

    k = 0;
    k2 = v1 - 1;
    for (j=1; j<=v1; j++) {
	vj = pan->pooled->list[j];
	if (vj == 0) {
	    relist[j] = 0;
	} else {
	    relist[j] = ++k;
	    strcpy(rset->varname[k], dset->varname[vj]);
	    if (hreg && j > 1 && var_is_varying(pan, vj)) {
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
			if (hreg && var_is_varying(pan, vj)) {
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

/* Construct a mini-dataset containing the group means, in
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
    int *droplist;
    int i, j, dpos;
    int err = 0;

    pmod->ci = PANEL;
    pmod->opt |= OPT_B;
    pmod->dw = NADBL;

    gretl_model_add_panel_varnames(pmod, gset, NULL);

    droplist = gretl_model_get_data(pmod, "droplist");

    /* replace both the model's regression list and its list of
       dropped variables, if any, with the ID numbers of the
       corresponding variables in the main dataset 
    */
    j = 1;
    for (i=1; i<=pmod->list[0]; i++) {
	if (i > 1 && droplist != NULL) {
	    while ((dpos = in_gretl_list(droplist, blist[j]))) {
		droplist[dpos] = pan->pooled->list[j];
		j++;
	    }
	}
	pmod->list[i] = pan->pooled->list[j++];
    }

    *pan->realmod = *pmod;

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
	pan->s2b = bmod.sigma * bmod.sigma;
    } else {
	err = bmod.errcode;
#if PDEBUG
	fprintf(stderr, "error %d in between_variance\n", err);
#endif
    }

    if (!err && (pan->opt & OPT_B)) {
	err = save_between_model(&bmod, blist, gset, pan);
    } else {
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

    for (i=0; i<pan->nbeta; i++) {
	if (vcv_skip(pmod, mi, pan, op)) {
	    i--;
	    mi++;
	    continue;
	}
	mj = mi;
	sj = si;
	for (j=i; j<pan->nbeta; j++) {
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

static void panel_df_correction (MODEL *pmod, int k)
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

/* The function below is used in the context of the "panel diagnostics"
   test (only): print the fixed-effects or "within" model in simple
   form.
*/

static int print_fe_results (panelmod_t *pan, 
			     MODEL *pmod,
			     DATASET *dset,
			     PRN *prn)
{
    int n, maxlen = 0;
    int dfn, i, vi;

    pputs(prn, 
	  _("Fixed effects estimator\n"
	    "allows for differing intercepts by cross-sectional unit\n"
	    "slope standard errors in parentheses, p-values in brackets\n"));
    pputc(prn, '\n');

    for (i=0; i<pmod->ncoeff; i++) {
	vi = pan->vlist[i+2];
	n = strlen(dset->varname[vi]);
	if (n > maxlen) {
	    maxlen = n;
	}
    }

    maxlen = maxlen < 12 ? 12 : maxlen;

    /* print the slope coefficients, for varying regressors */
    for (i=0; i<pmod->ncoeff; i++) {
	vi = pan->vlist[i+2];
	print_panel_coeff(pmod, dset->varname[vi], i, maxlen, prn);
    }
    pputc(prn, '\n');   

    pprintf(prn, _("%d group means were subtracted from the data"), pan->effn);
    pputc(prn, '\n');

    dfn = pan->effn - 1;
    pprintf(prn, _("\nResidual variance: %g/(%d - %d) = %g\n"), 
	    pmod->ess, pmod->nobs, pan->vlist[0] - 1 + dfn, pan->s2e);

    pprintf(prn, _("Joint significance of differing group means:\n"));
    pprintf(prn, " F(%d, %d) = %g %s %g\n", pan->Fdfn, pan->Fdfd, pan->F, 
	    _("with p-value"), snedecor_cdf_comp(pan->Fdfn, pan->Fdfd, pan->F));

    pputs(prn, _("(A low p-value counts against the null hypothesis that "
		 "the pooled OLS model\nis adequate, in favor of the fixed "
		 "effects alternative.)\n\n"));

    return 0;
}

static int time_dummies_wald_test (panelmod_t *pan, MODEL *wmod)
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
    bigk = wmod->ncoeff;

    err = makevcv(wmod, wmod->sigma);
    if (err) {
	return err;
    }

    b = gretl_column_vector_alloc(k);
    vcv = gretl_matrix_alloc(k, k);
    if (b == NULL || vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    di = bigk - k;
    for (i=0; i<k; i++) {
	b->val[i] = wmod->coeff[di++];
    }

    di = bigk - k;
    for (i=0; i<k; i++) {
	dj = bigk - k;
	for (j=0; j<=i; j++) {
	    x = wmod->vcv[ijton(di, dj++, bigk)];
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
	fprintf(stderr, _("Failed to compute test statistic\n"));
	goto bailout;
    }

    if (!err) {
	ModelTest *test = model_test_new(GRETL_TEST_PANEL_TIMEDUM);

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
	    model_test_set_dfn(test, k);
	    model_test_set_value(test, x);
	    model_test_set_pvalue(test, chisq_cdf_comp(k, x));
	    maybe_add_test_to_model(wmod, test);
	}
    }		

 bailout:

    gretl_matrix_free(vcv);
    gretl_vector_free(b);
    
    return err;
}

static void save_fixed_effects_F (panelmod_t *pan, MODEL *wmod)
{
    ModelTest *test = model_test_new(GRETL_TEST_PANEL_F);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_F);
	model_test_set_dfn(test, pan->Fdfn);
	model_test_set_dfd(test, pan->Fdfd);
	model_test_set_value(test, pan->F);
	model_test_set_pvalue(test, snedecor_cdf_comp(pan->Fdfn, pan->Fdfd, pan->F));
	maybe_add_test_to_model(wmod, test);
    }	    
}

static void fixed_effects_F (panelmod_t *pan, MODEL *wmod)
{
    pan->Fdfn = pan->effn - 1;
    pan->Fdfd = wmod->dfd;

    pan->F = (pan->pooled->ess - wmod->ess) * pan->Fdfd / 
	(wmod->ess * pan->Fdfn);

    if (pan->F < 0) {
	pan->F = 0;
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

/* Correct various model statistics, in the case where we estimated
   the fixed effects or "within" model on an auxiliary dataset
   from which the group means were subtracted.
*/

static int fix_within_stats (MODEL *fmod, panelmod_t *pan)
{
    double wfstt, wrsq;
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
    }

    wdfn = fmod->ncoeff - 1;

    wfstt = (wrsq / (1.0 - wrsq)) * ((double) fmod->dfd / wdfn);
    if (wfstt >= 0.0) {
	ModelTest *test = model_test_new(GRETL_TEST_WITHIN_F);

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_F);
	    model_test_set_dfn(test, wdfn);
	    model_test_set_dfd(test, fmod->dfd);
	    model_test_set_value(test, wfstt);
	    model_test_set_pvalue(test, snedecor_cdf_comp(wdfn, fmod->dfd, wfstt));
	    maybe_add_test_to_model(fmod, test);
	}
    }

    /* note: this member is being borrowed for the "Within R-squared" */
    fmod->adjrsq = wrsq;

    /* LSDV-based statistics */
    fmod->rsq = 1.0 - (fmod->ess / fmod->tss);
    if (fmod->rsq < 0.0) {
	fmod->rsq = 0.0;
    } else {
	fmod->fstt = (fmod->rsq / (1.0 - fmod->rsq)) * 
	    ((double) fmod->dfd / fmod->dfn);
    }

    fmod->ncoeff = fmod->dfn + 1; /* number of params estimated */
    ls_criteria(fmod);
    fmod->ncoeff = nc;

    return err;
}

/* Fixed-effects model: add the per-unit intercept estimates ("ahat")
   to the model in case the user wants to retrieve them.  By this
   point the model -- even if it been estimated on a short dataset --
   should have a full-length residual series.
*/

static int fe_model_add_ahat (MODEL *pmod, const DATASET *dset,
			      panelmod_t *pan)
{
    double *ahat = NULL;
    const double *x;
    int i, j, t, bigt;
    int err = 0;

    ahat = malloc(dset->n * sizeof *ahat);
    if (ahat == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<dset->n; t++) {
	ahat[t] = NADBL;
    }

    for (i=0; i<pan->nunits; i++) {
	if (pan->unit_obs[i] > 0) {
	    double a = 0.0;

	    /* a = y - Xb, where the 'b' is based on de-meaned data */

	    a = 0.0;
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

    err = gretl_model_set_data(pmod, "ahat", ahat, 
			       GRETL_TYPE_DOUBLE_ARRAY, 
			       dset->n * sizeof *ahat);

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

/* computation of $\hat{\sigma}^2_u$ a la Nerlove, if wanted */

static int nerlove_s2v (MODEL *pmod, const DATASET *dset,
			panelmod_t *pan)
{
    double amean, *ahat;
    const double *x;
    int i, j, t, k, bigt;
    int *list;
    int err = 0;

    ahat = malloc(pan->effn * sizeof *ahat);
    if (ahat == NULL) {
	return E_ALLOC;
    }

    list = real_FE_list(pan);
    if (list == NULL) {
	free(ahat);
	return E_ALLOC;
    }    

    amean = 0.0;
    k = 0;

    for (i=0; i<pan->nunits; i++) {
	if (pan->unit_obs[i] > 0) {
	    double a = 0.0;

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

	    a /= pan->unit_obs[i];
	    ahat[k++] = a;
	    amean += a;
	}
    }

    amean /= pan->effn;
    pan->s2v = 0.0;

    for (i=0; i<pan->effn; i++) {
	pan->s2v += (ahat[i] - amean) * (ahat[i] - amean);
    }

    pan->s2v /= pan->effn - 1;

    free(ahat);
    free(list);

    return err;
}

/* Fix uhat and yhat in two cases: 

   (a) when we estimated fixed effects using a de-meaned dataset we
   need to ensure that the uhat and yhat values get written to the
   right observation slots in relation to the full dataset, and also
   that the yhat values get corrected, putting the means back in.

   (b) when estimating the random effects model we need to compute
   residuals based on the untransformed data (and again, place them
   correctly in relation to the full dataset).

   The placement issue arises because the special datasets used in
   these cases are not necessarily of full length, since they are
   purged of missing observations.
*/

static void 
fix_panel_hatvars (MODEL *pmod, panelmod_t *pan, const double **Z)
{
    const double *y = NULL;
    double *uhat = pan->pooled->uhat;
    double *yhat = pan->pooled->yhat;
    int n = pan->pooled->full_n;
    int re_n = 0;
    double yht;
    int ti;
    int i, j, s, t;

    y = Z[pan->pooled->list[1]];

    if (pan->opt & OPT_U) {
	/* random effects model */
	pmod->ess = 0.0;
    } 

    s = 0;

    for (i=0; i<pan->nunits; i++) {
	int Ti = pan->unit_obs[i];

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
		re_n++;
		uhat[t] = y[t] - yht;
		pmod->ess += uhat[t] * uhat[t];
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
	pmod->sigma = sqrt(pmod->ess / (re_n - (pmod->ncoeff - 1)));
    }

    pmod->full_n = n;

    free(pmod->uhat);
    pmod->uhat = uhat;
    
    free(pmod->yhat);
    pmod->yhat = yhat;

    /* we've stolen these */
    pan->pooled->uhat = NULL;
    pan->pooled->yhat = NULL;
}

#if PDEBUG > 1

static void verbose_femod_print (MODEL *femod, DATASET *wset, 
				 PRN *prn)
{
    int i, j;

    pprintf(prn, "*** initial FE model (on within data)\n");
    printmodel(femod, wset, OPT_O, prn);

    fprintf(stderr, "femod: data series length = %d\n", wset->n);
    for (i=0; i<wset->n; i++) {
	fprintf(stderr, "femod.uhat[%d] = %g, ", i, femod->uhat[i]);
	fprintf(stderr, "data: ");
	for (j=0; j<wset->v; j++) {
	    fprintf(stderr, "%g ", wset->Z[j][i]);
	}
	fputc('\n', stderr);
    }    
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

    femod = lsq(felist, wset, OLS, lsqopt);

    if (femod.errcode) {
	fprintf(stderr, "femod.errcode = %d\n", femod.errcode);
    } else if ((pan->opt & OPT_F) && femod.list[0] < felist[0]) {
	femod.errcode = E_SINGULAR;
    } else {
	/* we estimated a bunch of group means, and have to
	   subtract degrees of freedom */
	panel_df_correction(&femod, pan->effn - 1);
#if PDEBUG > 1
	verbose_femod_print(&femod, wset, prn);
#endif
	if (pan->opt & OPT_F) {
	    /* estimating the FE model in its own right */
	    fix_panel_hatvars(&femod, pan, (const double **) dset->Z);
	    if (pan->opt & OPT_R) {
		panel_robust_vcv(&femod, pan, (const double **) wset->Z);
	    } else {
		femod_regular_vcv(&femod);
	    }
	} else if (pan->opt & OPT_N) {
	    femod.errcode = nerlove_s2v(&femod, dset, pan);
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

    /* regressors dropped at the stage of estmating the 
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
    int nc;

    fix_panelmod_list(pmod, pan);

    pmod->ybar = pan->pooled->ybar;
    pmod->sdy = pan->pooled->sdy;
    pmod->tss = pan->pooled->tss;

    pmod->rsq = NADBL;
    pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;
    pmod->rho = NADBL;
    pmod->dw = NADBL;

    nc = pmod->ncoeff;
    pmod->ncoeff = pmod->dfn + 1;
    ls_criteria(pmod);
    pmod->ncoeff = nc;
}

static void add_panel_obs_info (MODEL *pmod, panelmod_t *pan)
{
    gretl_model_set_int(pmod, "n_included_units", pan->effn);
    gretl_model_set_int(pmod, "Tmin", pan->Tmin);
    gretl_model_set_int(pmod, "Tmax", pan->Tmax);
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
	fe_model_add_ahat(pmod, dset, pan);
	save_fixed_effects_F(pan, pmod);
    } else {
	/* random effects */
	pmod->opt |= OPT_U;
	gretl_model_set_double(pmod, "within-variance", pan->s2e);
	gretl_model_set_double(pmod, "between-variance", pan->s2b);
	if (pan->balanced) {
	    gretl_model_set_double(pmod, "gls-theta", pan->theta);
	}
	gretl_model_add_panel_varnames(pmod, dset, NULL);
	fix_panel_hatvars(pmod, pan, Z);
	fix_gls_stats(pmod, pan);
	if (pan->opt & OPT_N) {
	    /* record use of Nerlove transformation */
	    pmod->opt |= OPT_N;
	}
    }

    /* compose list of dropped variables, if any */
    compose_panel_droplist(pmod, pan);

    if (!(pan->opt & OPT_A)) {
	set_model_id(pmod);
    }  

    if (pan->opt & OPT_F) {
	panel_dwstat(pmod, pan);
	if (!na(pmod->dw) && (pan->opt & OPT_I)) {
	    panel_DW_pvalue(pmod, pan, Z);
	}
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

	if (pan->opt & OPT_N) {
	    /* Nerlove */
	    den = femod.nobs;
	} else {
	    /* as per Greene: nT - n - K */
	    den = femod.nobs - pan->effn - (pan->vlist[0] - 2);
	}

	pan->s2e = femod.ess / den;
#if PDEBUG
	fprintf(stderr, "nT = %d, n = %d, K = %d\n", femod.nobs,
		pan->effn, pan->vlist[0] - 2);
	fprintf(stderr, "pan->s2e = %g / %d = %g\n", femod.ess,
	       den, pan->s2e);
#endif
	fixed_effects_F(pan, &femod);

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
	    clear_model(&femod);
	}
    }

    return err;
}

static void print_hausman_result (panelmod_t *pan, PRN *prn)
{
    if (na(pan->H)) {
	pputs(prn, _("\nHausman test matrix is not positive definite (this "
		     "result may be treated as\n\"fail to reject\" the random effects "
		     "specification).\n"));
    } else {
	pprintf(prn, _("\nHausman test statistic:\n"
		       " H = %g with p-value = prob(chi-square(%d) > %g) = %g\n"),
		pan->H, pan->nbeta, pan->H, 
		chisq_cdf_comp(pan->nbeta, pan->H));
	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the random effects\nmodel is consistent, in favor of the fixed "
		     "effects model.)\n"));
    }
}

static void save_hausman_result (panelmod_t *pan)
{
    ModelTest *test;

    if (pan->realmod == NULL) {
	return;
    }

    test = model_test_new(GRETL_TEST_PANEL_HAUSMAN);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
	model_test_set_dfn(test, pan->nbeta);
	model_test_set_value(test, pan->H);
	if (na(pan->H)) {
	    model_test_set_pvalue(test, NADBL);
	} else {
	    model_test_set_pvalue(test, chisq_cdf_comp(pan->nbeta, pan->H));
	}
	maybe_add_test_to_model(pan->realmod, test);
    }	    
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
    double URSS = NADBL;
    int i, k, err = 0;

    gretl_model_init(&remod, dset);

    /* GLS regression list */
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

    /* If OPT_N (--nerlove) was given, we've already calculated
       pan->s2v as the variance of the fixed effects. Otherwise
       we're using the Swamy and Arora method, and pan->s2v still 
       needs to be computed.

       Note: for unbalanced panels, theta will actually vary across 
       the units in the final calculation.
    */
    if (!(pan->opt & OPT_N)) {
	pan->s2v = pan->s2b - pan->s2e / pan->Tbar;
	if (pan->s2v < 0) {
	    pan->s2v = 0.0;
	}
    }

    /* theta, the quasi-demeaning coefficient */
    pan->theta = 1.0 - sqrt(pan->s2e / (pan->s2e + pan->Tmax * pan->s2v));

#if PDEBUG
    fprintf(stderr, "s_v = %.8g, s_e = %.8g\n", sqrt(pan->s2v), sqrt(pan->s2e));
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
	/* estimate the augmented random-effects model (GLS) */
	int *biglist = gretl_list_add(relist, hlist, &err);

	if (!err) {
	    remod = lsq(biglist, rset, OLS, OPT_A);
	    if (remod.errcode == 0) {
		/* record unrestricted SSR */
		URSS = remod.ess;
	    }
	    clear_model(&remod);
	}
	free(biglist);
    }	

    /* regular random-effects model */
    remod = lsq(relist, rset, OLS, lsqopt);

    if ((err = remod.errcode)) {
	pputs(prn, _("Error estimating random effects model\n"));
	errmsg(err, prn);
    } else {
	if (pan->opt & OPT_V) {
	    print_re_model_top(pan, prn);
	}

	k = 0;
	for (i=0; i<remod.ncoeff; i++) {
	    int vi = pan->pooled->list[i+2];

	    if (pan->opt & OPT_V) {
		print_panel_coeff(&remod, dset->varname[vi], i, 15, prn);
	    }
	    if (pan->bdiff != NULL && var_is_varying(pan, vi)) {
		pan->bdiff->val[k++] -= remod.coeff[i];
	    }
	}

	makevcv(&remod, remod.sigma);

	if (!na(URSS)) {
	    /* it appears that Baltagi uses T-k instead of T here. Why? */
	    pan->H = (remod.ess / URSS - 1.0) * (remod.nobs);
	} else if (pan->Sigma != NULL) {
	    vcv_slopes(pan, &remod, VCV_SUBTRACT);
	}
    }

#if PDEBUG > 1
    if (remod.errcode == 0) {
	printmodel(&remod, rset, OPT_NONE, prn);
    }
#endif	

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
	pputs(prn, _("\nMeans of pooled OLS residuals for cross-sectional "
		     "units:\n\n"));
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
	pprintf(prn, _("\nBreusch-Pagan test statistic:\n"
		       " LM = %g with p-value = prob(chi-square(1) > %g) = %g\n"), 
		pan->BP, pan->BP, chisq_cdf_comp(1, pan->BP));

	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the pooled OLS model\nis adequate, in favor of the random "
		     "effects alternative.)\n\n"));
    } 

    return 0;
}

static int finalize_hausman_test (panelmod_t *pan, PRN *prn)
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

    if (pan->opt & OPT_V) {
	if (!err || (mdiff && err == E_NOTPD)) {
	    print_hausman_result(pan, prn);
	} 
    } else {
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

    pan->NT = 0;
    pan->effn = 0;
    pan->Tmax = 0;
    pan->Tmin = pan->T;

    for (i=0; i<pan->nunits; i++) {
	uobs[i] = 0;
	for (t=0; t<pan->T; t++) {
	    bigt = panel_index(i, t);
#if PDEBUG > 1
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
    int err = 0;

    pan->opt = opt;
    pan->pooled = pmod;

    /* assumes (possibly padded) balanced panel dataset */
    pan->nunits = (dset->t2 - dset->t1 + 1) / dset->pd;
    pan->T = dset->pd;

    panel_index_init(dset, pan->nunits, pan->T);
    pan->ntdum = ntdum;
    
    err = panel_obs_accounts(pan);

    if (!err && (pan->opt & (OPT_U | OPT_F | OPT_B))) {
	pan->realmod = malloc(sizeof *pan->realmod);
	if (pan->realmod == NULL) {
	    err = E_ALLOC;
	}
    }
    
    if (err && pan->unit_obs != NULL) {
	free(pan->unit_obs);
	pan->unit_obs = NULL;
    }

    return err;
}

/* Called in relation to a model estimated by pooled OLS: test for
   both fixed and random effects.
*/

int panel_diagnostics (MODEL *pmod, DATASET *dset, 
		       gretlopt opt, PRN *prn)
{
    panelmod_t pan;
    int xdf, err = 0;

#if PDEBUG
    fputs("\n*** Starting panel_diagnostics ***\n", stderr);
#endif

    if (pmod->ifc == 0) {
	/* at many points we assume the base regression has an
	   intercept included */
	return E_NOCONST;
    }

    /* add OPT_V to make the fixed and random effects functions verbose */
    panelmod_init(&pan);
    err = panelmod_setup(&pan, pmod, dset, 0, opt | OPT_V);
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

    /* can we do the Hausman test or not? */
    if (xdf > 0) {
	err = hausman_allocate(&pan);
	if (err) {
	    goto bailout;
	}
    } 

    if (pan.balanced) {
	pprintf(prn, _("      Diagnostics: assuming a balanced panel with %d "
		       "cross-sectional units\n "
		       "                        observed over %d periods\n\n"), 
		pan.effn, pan.Tmax);
    }

    err = within_variance(&pan, dset, prn);
    if (err) {
	goto bailout;
    }

    breusch_pagan_LM(&pan, prn);

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
	    finalize_hausman_test(&pan, prn);
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
    }

    if (gset != NULL) {
	destroy_dataset(gset);
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

static int
add_dummies_to_list (const int *list, DATASET *dset, 
		     int **plist)
{
    char dname[VNAMELEN];
    int *biglist = NULL;
    int ndum = dset->pd - 1;
    int i, j, v;
    int err = 0;

    biglist = gretl_list_new(list[0] + ndum);
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

static void save_pooled_model (MODEL *pmod, panelmod_t *pan,
			       const double **Z)
{
    gretl_model_set_int(pmod, "pooled", 1);
    add_panel_obs_info(pmod, pan);
    
    if (!(pan->opt & OPT_A)) {
	set_model_id(pmod);
    }

    if (pan->opt & OPT_R) {
	panel_robust_vcv(pmod, pan, Z);
    }

    panel_dwstat(pmod, pan);

    if (!na(pmod->dw) && (pan->opt & OPT_I)) {
	panel_DW_pvalue(pmod, pan, Z);
    }
}

/* fixed effects | random effects | between | pooled */

#define estimator_specified(o) (o & (OPT_F|OPT_U|OPT_B|OPT_P))

/* real_panel_model:
 * @list: list containing model specification.
 * @dset: dataset struct.
 * @opt: may include %OPT_U for the random effects model;
 * %OPT_R for robust standard errors (fixed effects model
 * and pooled OLS only); %OPT_M to use the matrix-difference 
 * version of the Hausman test (random effects only); %OPT_B for 
 * the "between" model; %OPT_P for pooled OLS; and %OPT_D to 
 * include time dummies. If and only if %OPT_P is given, %OPT_C
 * (clustered standard errors) is accepted.
 * @prn: printing struct.
 *
 * Estimates a panel model, by default the fixed effects model.
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
    int ntdum = 0;
    int err = 0;

    gretl_model_init(&mod, dset);
    panelmod_init(&pan);

    if (opt & OPT_P) {
	/* doing pooled OLS */
	if (opt & OPT_C) {
	    /* clustered */
	    ols_opt |= (OPT_C | OPT_R);
	}
    } else {
	/* not just pooled OLS */
	mod.errcode = panel_check_for_const(list);
	if (mod.errcode) {
	    return mod;
	}
    }

    /* add time dummies to list? */
    if (opt & OPT_D) {
	err = panel_dummies(dset, OPT_T, prn);
	if (!err) {
	    err = add_dummies_to_list(list, dset, &olslist);
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
	fprintf(stderr, "real_panel_model: error %d in intial OLS\n", 
		mod.errcode);
	goto bailout;
    } 

    free(olslist);

#if PDEBUG
    pprintf(prn, "*** initial baseline OLS\n");
    printmodel(&mod, dset, OPT_NONE, prn);
#endif

    if (!estimator_specified(opt)) {
	/* default: add OPT_F to save the fixed effects model */
	pan_opt |= OPT_F;
    }

#if PDEBUG
    if (pan_opt & OPT_F) {
	fprintf(stderr, "\n*** Doing fixed effects\n");
    } else if (pan_opt & OPT_U) {
	fprintf(stderr, "\n*** Doing random effects\n");
    }
#endif

    if ((opt & OPT_D) && !(opt & OPT_B)) {
	ntdum = get_ntdum(list, mod.list);
    }

    err = panelmod_setup(&pan, &mod, dset, ntdum, pan_opt);
    if (err) {
	goto bailout;
    }   

    if (opt & OPT_P) {
	save_pooled_model(&mod, &pan, (const double **) dset->Z);
	goto bailout;
    }

    if (opt & OPT_B) {
	err = between_model(&pan, dset);
	goto bailout;
    }	

    calculate_Tbar(&pan);

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
	/* trying to do random effects */
	int xdf = pan.effn - mod.ncoeff;

	if (xdf <= 0) {
#if PDEBUG
	    fprintf(stderr, "xdf = %d - %d = %d\n", pan.effn, mod.ncoeff, xdf);
#endif
	    err = mod.errcode = E_DF;
	    goto bailout;
	} else {
	    err = hausman_allocate(&pan);
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

	if (err) { 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	} else {
	    err = random_effects(&pan, dset, gset, prn);
	    if (!err) {
		save_breusch_pagan_result(&pan);
		finalize_hausman_test(&pan, prn);
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
	    clear_model(&mod);
	    mod = *pan.realmod;
	}
	gretl_model_smpl_init(&mod, dset);
	if ((opt & OPT_D) && pan.ntdum > 0) {
	    maybe_suppress_time_dummies(&mod, pan.ntdum);
	}
	if (complex_subsampled() && (opt & OPT_D)) {
	    process_time_dummies(&mod, dset, orig_v);
	}
    }

    panelmod_free(&pan);

    if (err && mod.errcode == 0) {
	mod.errcode = err;
    }

    if (complex_subsampled()) {
	dataset_drop_last_variables(dset, dset->v - orig_v);
    }

    return mod;    
}

/* Called from qr_estimate.c in case robust VCV estimation is called
   for in the context of TSLS estimation on panel data.
*/

int panel_tsls_robust_vcv (MODEL *pmod, const DATASET *dset)
{
    panelmod_t pan;
    int err = 0;

    panelmod_init(&pan);

    err = panelmod_setup(&pan, pmod, dset, 0, OPT_NONE);
    if (!err) {
	err = panel_robust_vcv(pmod, &pan, (const double **) dset->Z);
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
    gretlopt wlsopt = OPT_A;
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
	if (singleton_check(pan.unit_obs, pan.nunits)) {
	    pprintf(prn, _("Can't produce ML estimates: "
			   "some units have only one observation"));
	    pputc(prn, '\n');
	    opt ^= OPT_I;
	}
    }

    s2 = mdl.ess / mdl.nobs;

    if ((opt & OPT_V) && (opt & OPT_I)) {
	pprintf(prn, "\nOLS error variance = %g\n", s2);
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
	}

	if (iter > WLS_MAX) {
	    mdl.errcode = E_NOCONV;
	    break;
	}

	if (opt & OPT_I) {
	    diff = max_coeff_diff(&mdl, bvec);
	    if ((opt & OPT_V) && iter == 1) {
		pprintf(prn, "\nFGLS pooled error variance = %g\n",
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
	    set_model_id(&mdl);
	}
	gretl_model_set_int(&mdl, "n_included_units", pan.effn);
	mdl.opt |= OPT_W;
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
 * to @pmod.
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

	print_wald_test(W, df, pval, &pan, uvar, s2, prn);

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

	record_test_result(W, pval, _("groupwise heteroskedasticity"));
    }

    free(pan.unit_obs);
    free(uvar);

    return err;
}

static void panel_copy_var (DATASET *targ, int targv,
			    double *x, DATASET *src, int srcv,
			    int order)
{
    int t, j = 0;

    for (t=src->t1; t<=src->t2; t++) {
	if (t % src->pd >= order) {
	    targ->Z[targv][j++] = x[t];
	} 
    }

    if (srcv == -1) {
	strcpy(targ->varname[targv], "uhat");
	series_set_label(targ, targv, _("residual"));
    } else {
	strcpy(targ->varname[targv], src->varname[srcv]);
	series_set_label(targ, targv, series_get_label(src, srcv));
    }
}

static void make_reduced_data_info (DATASET *targ, DATASET *src, 
				    int order)
{
    targ->pd = src->pd - order;
    ntodate(targ->stobs, src->t1 + order, src);
    targ->sd0 = obs_str_to_double(targ->stobs); 
    targ->structure = src->structure;
}

static void panel_lag (DATASET *tmpset, 
		       double *src, DATASET *srcset, 
		       int v, int order, int lag)
{
    int t, j = 0;

    for (t=srcset->t1; t<=srcset->t2; t++) {
	if (t % srcset->pd >= order) {
	    tmpset->Z[v][j++] = src[t - lag];
	}
    }

    sprintf(tmpset->varname[v], "uhat_%d", lag);
    series_set_label(tmpset, v, "");
}

/* - do some sanity checks
   - create a local copy of the required portion of the data set,
     skipping the obs that will be missing
   - copy in the lags of uhat
   - estimate the aux model
   - destroy the temporary data set
*/

int panel_autocorr_test (MODEL *pmod, int order, DATASET *dset, 
			 gretlopt opt, PRN *prn)
{
    int *aclist;
    DATASET *tmpset;
    MODEL aux;
    double trsq, LMF;
    int i, nv, nunits, nobs, err = 0;
    int sn = dset->t2 - dset->t1 + 1;

    if (pmod->ci != OLS) {
	return E_NOTIMP;
    }

    if (pmod->missmask != NULL) {
	return E_DATA;
    }

    /* basic checks */
    if (order <= 0) order = 1;
    if (order > dset->pd - 1) return E_DF;
    if (pmod->ncoeff + order >= sn) return E_DF;

    /* get number of cross-sectional units */
    nunits = sn / dset->pd;

    /* we lose "order" observations for each unit */
    nobs = sn - nunits * order;

    /* the required number of variables */
    nv = pmod->list[0] + order;

    /* create temporary reduced dataset */
    tmpset = create_auxiliary_dataset(nv, nobs, 0);
    if (tmpset == NULL) {
	return E_ALLOC;
    }

    make_reduced_data_info(tmpset, dset, order);

#if PDEBUG
    fprintf(stderr, "Created data set, n=%d, pd=%d, vars=%d, stobs='%s'\n", 
	    tmpset->n, tmpset->pd, tmpset->v, tmpset->stobs);
#endif

    /* allocate the auxiliary regression list */
    aclist = malloc((nv + 1) * sizeof *aclist);
    if (aclist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	int k, v;

	aclist[0] = pmod->list[0] + order;
	/* copy model uhat to position 1 in temp data set */
	aclist[1] = 1;
	panel_copy_var(tmpset, 1,
		       &pmod->uhat[0], dset, -1,
		       order);
	/* copy across the original indep vars, making
	   the new regression list while we're at it */
	k = 2;
	for (i=2; i<=pmod->list[0]; i++) {
	    v = pmod->list[i];
	    if (v == 0) { /* the constant */
		aclist[i] = 0;
	    } else {
		aclist[i] = k;
		panel_copy_var(tmpset, k, 
			       dset->Z[v], dset, v, 
			       order);
		k++;
	    }
	}
    }

    if (!err) {
	int v = pmod->list[0] - 1;

	/* add lags of uhat to temp data set */
	for (i=1; i<=order; i++) {
	    panel_lag(tmpset, pmod->uhat, 
		      dset, v + i, order, i);
	    aclist[v + i + 1] = v + i;
	}
    }

    if (!err) {
	aux = lsq(aclist, tmpset, OLS, OPT_A);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    } 

    if (!err) {
	int dfd = aux.nobs - pmod->ncoeff - order;
	double pval;

	aux.aux = AUX_AR;
	gretl_model_set_int(&aux, "BG_order", order);
	printmodel(&aux, tmpset, OPT_NONE, prn);
	trsq = aux.rsq * aux.nobs;
	/* FIXME LMF */
	LMF = (aux.rsq / (1.0 - aux.rsq)) * dfd / order; 
	pval = snedecor_cdf_comp(order, dfd, LMF);

	pprintf(prn, "\n%s: LMF = %f,\n", _("Test statistic"), LMF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, dfd, LMF, pval);

	pprintf(prn, "\n%s: TR^2 = %f,\n", 
		_("Alternative statistic"), trsq);
	pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", 	_("with p-value"), 
		_("Chi-square"), order, trsq, chisq_cdf_comp(order, trsq));

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_AUTOCORR);

	    if (test != NULL) {
		int dfd = aux.nobs - pmod->ncoeff - order;

		model_test_set_teststat(test, GRETL_STAT_LMF);
		model_test_set_order(test, order);
		model_test_set_dfn(test, order);
		model_test_set_dfd(test, dfd);
		model_test_set_value(test, LMF);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	    
	}
    }

    free(aclist);
    clear_model(&aux); 

    destroy_dataset(tmpset);

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

    if (oldmode != STACKED_TIME_SERIES &&
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
    ntodate(dset->stobs, 0, dset);
    ntodate(dset->endobs, dset->n - 1, dset);

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
    int i, j, k, t;
    int bigt;

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

	if (vj == 0) {
	    pan->vlist[k++] = 0;
	    pan->vlist[0] += 1;
	    continue;
	}

	for (i=0; i<pan->nunits; i++) { 
	    int started = 0;
	    double xval = NADBL;

	    if (pan->unit_obs[i] == 0) {
		continue;
	    }

	    for (t=0; t<pan->T; t++) {
		bigt = panel_index(i, t);
		if (panel_missing(pan, bigt)) {
		    continue;
		}
		if (!started) {
		    xval = dset->Z[vj][bigt];
		    started = 1;
		} else if (dset->Z[vj][bigt] != xval) {
		    varies = 1;
		    break;
		}
	    }
	    if (varies) break;
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

    return 0;
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

static int uv_tv_from_line (const char *line, const DATASET *dset,
			    int *uv, int *tv)
{
    char uvname[VNAMELEN];
    char tvname[VNAMELEN];
    char fmt[12];
    int err = 0;

    sprintf(fmt, "%%%ds %%%ds", VNAMELEN-1, VNAMELEN-1);
    
    if (sscanf(line, fmt, uvname, tvname) != 2) {
	return E_PARSE;
    }

    *uv = series_index(dset, uvname);
    if (*uv == dset->v) {
	/* FIXME "not a series" */
	gretl_errmsg_sprintf(_("Unknown variable '%s'"), uvname);
	err = E_UNKVAR;
    } 

    if (err) {
	return err;
    }

    *tv = series_index(dset, tvname);
    if (*tv == dset->v) {
	gretl_errmsg_sprintf(_("Unknown variable '%s'"), tvname);
	err = E_UNKVAR;
    } 

    return err;
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
    ntodate(dset->stobs, 0, dset); 
    ntodate(dset->endobs, dset->n - 1, dset);
}

static int check_full_dataset (void)
{
    DATASET *fset = fetch_full_dataset();

    if (!dataset_is_panel(fset)) {
	const char *msg =
	    "You cannot use the --panel-vars option with the setobs command when\n"
	    "\n"
	    "* the dataset is currently sub-sampled and\n"
	    "* the full dataset is not a panel.\n"
	    "\n"
	    "If you first structure your full dataset as a panel, you can then\n"
	    "do what you are trying to do.";

	gretl_errmsg_set(msg);
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

int set_panel_structure_from_line (const char *line, DATASET *dset)
{
    int n = dset->n;
    int uv = 0, tv = 0;
    int err = 0;

    if (!strncmp(line, "setobs", 6)) {
	line += 7;
    }

    err = uv_tv_from_line(line, dset, &uv, &tv);

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
		gretl_errmsg_sprintf("The string '%s' is given for "
				     "two or more groups", S[i]);
		return E_DATA;
	    }
	}
    }

    return 0;
}

/* 2013-09-21: enable construction of a string-valued series
   holding a name for each panel group/unit. The name of this
   series is set on the 'pangrps' member of @dset, and the
   names are then used in panel graphs where appropriate.

   @line should contain: (a) the name of a series, either an
   existing one (which will be overwritten) or a new one to
   create; and (b) either a string literal or the name of
   a string variable, which in either case should hold N
   space-separated strings, where N is the number of panel
   groups.
*/

int set_panel_group_strings (const char *line, DATASET *dset)
{
    char vname[VNAMELEN];
    char strname[VNAMELEN];
    char *namestr = NULL;
    char **S = NULL;
    int ng = dset->n / dset->pd;
    int v, orig_v = dset->v;
    int freeit = 0;
    int err = 0;

    if (!strncmp(line, "setobs", 6)) {
	line += 7;
	line += strspn(line, " ");
    }

    if (sscanf(line, "%31s", vname) != 1) {
	err = E_PARSE;
    } else {
	/* get the string containing the names */
	line += strcspn(line, " ");
	line += strspn(line, " ");
	if (*line == '"') {
	    /* got a string literal */
	    const char *s = strchr(line + 1, '"');
	    
	    if (s == NULL) {
		err = E_PARSE;
	    } else {
		namestr = gretl_strndup(line + 1, s - line - 1);
		if (namestr == NULL) {
		    err = E_ALLOC;
		} else {
		    freeit = 1;
		}
	    }
	} else {
	    /* should be the name of a string variable */
	    if (sscanf(line, "%31s", strname) != 1) {
		err = E_PARSE;
	    } else if (!gretl_is_string(strname)) {
		err = E_DATA;
	    } else {
		namestr = get_string_by_name(strname);
	    }
	}
    }
    
    if (namestr != NULL) {
	int ngtest = 0;

	if (strchr(namestr, '"') != NULL) {
	    S = gretl_string_split_quoted(namestr, &ngtest, NULL, &err);
	} else {
	    S = gretl_string_split(namestr, &ngtest, " \n\t");
	}
	if (!err && S == NULL) {
	    err = E_ALLOC;
	}
	if (!err && ngtest != ng) {
	    fprintf(stderr, "Got %d strings but there are %d groups\n",
		    ngtest, ng);
	    err = E_DATA;
	}
	if (!err) {
	    err = group_uniqueness_check(S, ng);
	}
	if (freeit) {
	    free(namestr);
	}
    }

    if (!err) {
	v = current_series_index(dset, vname);
	if (v < 0) {
	    /* we need to add a series */
	    char *gen = gretl_strdup_printf("series %s", vname);

	    err = generate(gen, dset, OPT_Q, NULL);
	    if (!err) {
		v = dset->v - 1;
	    }
	    free(gen);
	}
    }

    if (!err) {
	series_table *st = series_table_new(S, ng);

	if (st == NULL) {
	    err = E_ALLOC;
	} else {
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
	strcmp(test, "period") == 0;
}

/* See if the panel dataset contains a variable that plausibly represents
   the year of the observations. We look for series labeled "year" or
   "period" (in a case insensitive comparison) and, if found, check
   that the series has the same sequence of values for each individual,
   and the same increment between successive time-series observations.
   (The increment does not have to be 1, to accommodate, e.g.,
   quinquennial or decennial observations.)

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

/* FIXME: this does not yet handle the dropping of instruments */

static int *arbond_list_omit (const MODEL *orig, const int *drop, int *err)
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

    if (orig->ci == ARBOND) {
	return arbond_list_omit(orig, drop, err);
    }

    /* sorry, can't drop the constant */
    if (drop != NULL) {
	int cpos = in_gretl_list(drop, 0);

	if (cpos >= 2) {
	    gretl_errmsg_set("Panel models must include an intercept");
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

static int *arbond_list_add (const MODEL *orig, const int *add, int *err)
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
    int *newlist = NULL;
    int i;

    if (orig->ci == ARBOND) {
	return arbond_list_add(orig, add, err);
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
	    newlist = gretl_list_add(panlist, add, err);
	    free(panlist);
	}
    } else {
	newlist = gretl_list_add(orig->list, add, err);
    }
	
    return newlist;
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
