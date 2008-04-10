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

#define PDEBUG 0

enum vcv_ops {
    VCV_INIT,
    VCV_SUBTRACT
};

typedef struct panelmod_t_ panelmod_t;

struct panelmod_t_ {
    gretlopt opt;         /* option flags */
    int nunits;           /* total cross-sectional units in sample range */
    int effn;             /* effective (included) cross-section units */
    int N_fe;             /* number of units included in fixed effects model */
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
    double sigma_e;       /* fixed-effects standard error */
    double theta;         /* quasi-demeaning coefficient */
    double F;             /* joint significance of differing unit intercepts */
    double BP;            /* Breusch-Pagan test statistic */
    double H;             /* Hausman test statistic */
    double *bdiff;        /* array of coefficient differences */
    double *sigma;        /* Hausman covariance matrix */
    double between_s2;    /* group-means or "between" error variance */
    double s2e;           /* \hat{sigma}^2_e, from fixed-effects estimator */
    double s2u;           /* \hat{sigma}^2_u measure */
    int *small2big;       /* data indexation array */
    int *big2small;       /* reverse data indexation array */
    MODEL *pooled;        /* reference model (pooled OLS) */
    MODEL *realmod;       /* fixed or random effects model */
};

struct {
    int n;      /* number of cross-sectional units */
    int T;      /* number of observations per unit */
    int offset; /* sampling offset into full dataset */
} panidx;

static int 
varying_vars_list (const double **Z, const DATAINFO *pdinfo,
		   panelmod_t *pan);

/* translate from (i = unit, t = time period for that unit) to
   overall 0-based index into the data set */
#define panel_index(i,t) (i * panidx.T + t + panidx.offset)

#define panel_missing(p, t) (na(p->pooled->uhat[t]))


static void 
panel_index_init (const DATAINFO *pdinfo, int nunits, int T)
{
    panidx.n = nunits;
    panidx.T = T;
    panidx.offset = pdinfo->t1;

#if PDEBUG
    fprintf(stderr, "panel_index_init: n=%d, T=%d, offset=%d\n", 
	    panidx.n, panidx.T, panidx.offset);
#endif
}

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
#define big_index(p,t) ((p->small2big == NULL)? t : p->small2big[t])

static void panelmod_init (panelmod_t *pan)
{
    pan->nunits = 0;
    pan->effn = 0;
    pan->N_fe = 0;
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
    pan->sigma_e = NADBL;
    pan->theta = NADBL;

    pan->F = NADBL;
    pan->Fdfn = 0;
    pan->Fdfd = 0;

    pan->BP = NADBL;
    pan->H = NADBL;

    pan->bdiff = NULL;
    pan->sigma = NULL;

    pan->small2big = NULL;
    pan->big2small = NULL;
    
    pan->pooled = NULL;
    pan->realmod = NULL;
}

static void panelmod_free (panelmod_t *pan)
{
    free(pan->unit_obs);
    free(pan->varying);
    free(pan->vlist);

    free(pan->bdiff);
    free(pan->sigma);
    
    free(pan->small2big);
    free(pan->big2small);

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

static gretl_matrix *panel_model_xpxinv (MODEL *pmod, 
					 const panelmod_t *pan,
					 int *err)
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
	*err = makevcv(pmod, 1.0); /* invert X'X into pmod->vcv */
	if (*err) {
	    return NULL;
	}
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
	gretl_model_set_int(pmod, "panel_bk_failed", 1);
	goto bailout;
    }

    /* form V = (Xi'Xi)^{-1} W (Xi'Xi)^{-1} */
    gretl_matrix_qform(XX, GRETL_MOD_NONE, W,
		       V, GRETL_MOD_NONE);

    gretl_model_set_int(pmod, "panel_bk", 1);

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

    e   = gretl_vector_alloc(T);
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

    gretl_model_set_int(pmod, "panel_hac", 1);

 bailout:

    gretl_matrix_free(e);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    return err;
}

/* common setup for Arellano and Beck-Katz VCV estimators */

static int 
panel_robust_vcv (MODEL *pmod, panelmod_t *pan, const double **Z)
{
    gretl_matrix *W = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *XX = NULL;
    int k = pmod->ncoeff;
    int err = 0;

    W  = gretl_zero_matrix_new(k, k);
    V  = gretl_matrix_alloc(k, k);

    if (W == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    XX = panel_model_xpxinv(pmod, pan, &err);
    if (err) {
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

/* Durbin-Watson statistic for the fixed effects model.  We only
   use units that have at least two time-series observations.
*/

static void panel_dwstat (MODEL *pmod, const panelmod_t *pan)
{
    double ut, u1;
    double num = 0.0;
    double den = 0.0;
    int i, t, ti, Ti;
    int started;

    pmod->dw = NADBL;
    pmod->rho = NADBL;

    if (pmod->ess <= 0.0) {
	return;
    }

    for (i=0; i<pan->nunits; i++) {
	Ti = pan->unit_obs[i];

	if (Ti < 2) {
	    continue;
	}

	started = 0;
	for (ti=1; ti<pan->T; ti++) {
	    t = panel_index(i, ti);
	    ut = pmod->uhat[t];
	    u1 = pmod->uhat[t-1];
	    if (!na(ut) && !na(u1)) {
		num += (ut - u1) * (ut - u1);
		den += ut * ut;
		if (!started) {
		    den += u1 * u1;
		    started = 1;
		}
	    }
	}
    }

    if (den > 0.0) {
	pmod->dw = num / den;
    }
}

/* Allocate the arrays needed to perform the Hausman test,
   in its matrix formulation.
*/

static int hausman_allocate (panelmod_t *pan)
{
    int ns, k = pan->vlist[0] - 2;

    pan->nbeta = k;

    if (pan->opt & OPT_H) {
	/* taking the regression approach to Hausman: we don't need
	   the allocations below */
	return 0;
    }

    /* array to hold differences between coefficient estimates */
    pan->bdiff = malloc(k * sizeof *pan->bdiff);
    if (pan->bdiff == NULL) {
	return E_ALLOC;
    }

    ns = k * (k + 1) / 2;

    /* array to hold covariance matrix */
    pan->sigma = malloc(ns * sizeof *pan->sigma);
    if (pan->sigma == NULL) {
	free(pan->bdiff);
	pan->bdiff = NULL;
	return E_ALLOC; 
    }

    return 0;
}   

/* printing routines for the "panel diagnostics" test */

#if 0

static void print_panel_coeff (const MODEL *pmod,
 			       const char *vname,
 			       int i, PRN *prn)
{
    model_coeff mc;

    model_coeff_init(&mc);
    mc.b = pmod->coeff[i];
    mc.se = pmod->sderr[i];
    mc.tval = mc.b / mc.se;
    mc.pval = student_pvalue_2(pmod->dfd, mc.tval);
    strcpy(mc.name, vname);
    print_coeff(&mc, prn);
}

#else

static void print_panel_coeff (const MODEL *pmod,
 			       const char *vname,
 			       int i, PRN *prn)
{
    double tstat = pmod->coeff[i] / pmod->sderr[i];
    char errstr[18];
    char pvstr[18];
 
    sprintf(errstr, "(%.5g)", pmod->sderr[i]);
    sprintf(pvstr, "[%.5f]", student_pvalue_2(pmod->dfd, tstat));
    pprintf(prn, "%*s: %14.5g %15s %15s\n", VNAMELEN, vname,
 	    pmod->coeff[i], errstr, pvstr);
}

#endif

static void print_re_model_top (const panelmod_t *pan, PRN *prn)
{
    pputs(prn, "Variance estimators:\n");
    pprintf(prn, " between = %g\n", pan->between_s2);
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

/* With FE_NEW we follow stata's approach for fixed effects: we
   subtract the groups means but add back in the grand means, for
   each variable. That way, we can estimate an "average" constant
   and get its standard error.
*/

#define FE_NEW 1

/* Construct a version of the dataset from which the group means are
   subtracted, for the "within" regression.  Nota bene: this auxiliary
   dataset is not necessarily of full length: missing observations
   are skipped.
*/

static DATAINFO *
within_groups_dataset (const double **Z, const DATAINFO *pdinfo,
		       double ***wZ, panelmod_t *pan)
{
    DATAINFO *winfo = NULL;
    double *xbar = NULL;
    int *vlist = NULL;
    int i, j, vj;
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

    if (pan->NT < pdinfo->n) {
	err = allocate_data_finders(pan, pdinfo->n);
	if (err) {
	    return NULL;
	}
    }

    vlist = real_varying_list(pan);
    if (vlist == NULL) {
	goto bailout;
    }

    xbar = malloc(vlist[0] * sizeof *xbar);
    if (xbar == NULL) {
	goto bailout;
    }

#if PDEBUG
    fprintf(stderr, "within_groups dataset: nvars=%d, nobs=%d\n", 
	    pan->vlist[0], pan->NT);
#endif

    winfo = create_auxiliary_dataset(wZ, pan->vlist[0], pan->NT);
    if (winfo == NULL) {
	goto bailout;
    }

#if FE_NEW

    for (j=1; j<=vlist[0]; j++) {
	double xbar, gxbar = 0.0;

	vj = vlist[j];
	gxbar = 0.0;
	s = 0;

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
		    xbar += Z[vj][bigt];
		}
	    }

	    gxbar += xbar;
	    xbar /= Ti;
	    
	    /* second pass: calculate de-meaned values */
	    got = 0;
	    for (t=0; t<pan->T && got<Ti; t++) { 
		bigt = panel_index(i, t);
		if (!panel_missing(pan, bigt)) {
		    (*wZ)[j][s] = Z[vj][bigt] - xbar;
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
	    (*wZ)[j][s] += gxbar;
	}
    }

#else

    s = 0;

    for (i=0; i<pan->nunits; i++) { 
	int Ti = pan->unit_obs[i];
	int got = 0;

	if (Ti == 0) {
	    continue;
	}

	/* first pass: find the group means */
	for (j=0; j<vlist[0]; j++) {
	    xbar[j] = 0.0;
	}
	for (t=0; t<pan->T; t++) {
	    bigt = panel_index(i, t);
	    if (!panel_missing(pan, bigt)) {
		for (j=0; j<vlist[0]; j++) {
		    vj = vlist[j+1];
		    xbar[j] += Z[vj][bigt];
		}
	    }
	}
	for (j=0; j<vlist[0]; j++) {
	    xbar[j] /= Ti;
	}	

	/* second pass: calculate de-meaned values */
	for (t=0; t<pan->T && got<Ti; t++) { 
	    bigt = panel_index(i, t);
	    if (!panel_missing(pan, bigt)) {
		for (j=0; j<vlist[0]; j++) {
		    vj = vlist[j+1];
		    (*wZ)[j+1][s] = Z[vj][bigt] - xbar[j];
		}
		got++;
		if (pan->small2big != NULL) {
		    pan->small2big[s] = bigt;
		    pan->big2small[bigt] = s;
		}
		s++;
	    }
	}	    
    }

#endif

 bailout:

    free(vlist);
    free(xbar);

    return winfo;
}

/* Construct a quasi-demeaned version of the dataset so we can apply
   least squares to estimate the random effects model.  This dataset
   is not necessarily of full length.  If we're implementing the
   Hausman test using the regression approach this dataset will also
   include "straight" de-meaned counterparts of the time-varying
   variables.
*/

static DATAINFO *
random_effects_dataset (const double **Z, const DATAINFO *pdinfo,
			const double **gZ, const DATAINFO *ginfo,
			double ***reZ, int *relist, int *hlist, 
			panelmod_t *pan)
{
    DATAINFO *reinfo;
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

    if (pan->NT < pdinfo->n) {
	err = allocate_data_finders(pan, pdinfo->n);
	if (err) {
	    return NULL;
	}
    }

#if PDEBUG
    fprintf(stderr, "random_effects_dataset: nvars=%d, nobs=%d\n",
	    v1 + v2, pan->NT);
#endif

    reinfo = create_auxiliary_dataset(reZ, v1 + v2, pan->NT);
    if (reinfo == NULL) {
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
	    strcpy(reinfo->varname[k], pdinfo->varname[vj]);
	    if (hreg && j > 1 && var_is_varying(pan, vj)) {
		k2++;
		strcpy(reinfo->varname[k2], "_");
		strncat(reinfo->varname[k2], pdinfo->varname[vj],
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
	    theta_i = 1.0 - sqrt(pan->s2e / (Ti * pan->s2u + pan->s2e));
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
			(*reZ)[0][s] -= theta_i;
		    } else {
			k++;
			xbar = (k < ginfo->v)? gZ[k][u] : 1.0 / pan->Tmax;
			(*reZ)[k][s] = Z[vj][bigt] - theta_i * xbar;
			if (hreg && var_is_varying(pan, vj)) {
			    (*reZ)[++k2][s] = Z[vj][bigt] - xbar;
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

    return reinfo;
}

#define BETWEEN_DEBUG 0

/* Construct a mini-dataset containing the group means, in
   order to run the group-means or "between" regression.
*/

static DATAINFO *
group_means_dataset (panelmod_t *pan,
		     const double **Z, const DATAINFO *pdinfo,
		     double ***gZ)
{
    DATAINFO *ginfo;
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

    ginfo = create_auxiliary_dataset(gZ, gv, gn);
    if (ginfo == NULL) {
	return NULL;
    }

    k = 1;
    for (j=1; j<=gv; j++) { 
	int vj = pan->pooled->list[j];

	if (vj == 0) {
	    continue;
	}

#if BETWEEN_DEBUG
	strcpy(ginfo->varname[k], pdinfo->varname[vj]);
#else
	if (pan->opt & OPT_B) {
	    /* will save the "between" model: so name the variables */
	    strcpy(ginfo->varname[k], pdinfo->varname[vj]);
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
		    x += Z[vj][bigt];
		}
	    }
	    (*gZ)[k][s++] = x / Ti;
	}
	k++;
    }

#if BETWEEN_DEBUG
    if (1) {
	PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	printdata(NULL, NULL, (const double **) *gZ, ginfo, OPT_O, prn);
	gretl_print_destroy(prn);
    }
#endif

    return ginfo;
}

/* spruce up the between model and attach it to pan */

static int save_between_model (MODEL *pmod, const int *blist,
			       DATAINFO *ginfo, panelmod_t *pan)
{
    int *droplist;
    int i, j, dpos;
    int err = 0;

    pmod->ci = PANEL;
    gretl_model_set_int(pmod, "between", 1);
    pmod->dw = NADBL;

    gretl_model_add_panel_varnames(pmod, ginfo, NULL);

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

static int
between_variance (panelmod_t *pan, double ***gZ, DATAINFO *ginfo)
{
    gretlopt bopt;
    MODEL bmod;
    int *blist;
    int i, j, k;
    int err = 0;

    blist = gretl_list_new(ginfo->v);
    if (blist == NULL) {
	return E_ALLOC;
    }

    j = k = 1;
    for (i=1; i<=ginfo->v; i++) { 
	if (pan->pooled->list[i] == 0) {
	    blist[k++] = 0;
	} else {
	    blist[k++] = j++;
	} 
    }

    bopt = (pan->opt & OPT_B)? OPT_NONE : OPT_A;
    bmod = lsq(blist, gZ, ginfo, OLS, bopt);

    if (bmod.errcode == 0) {
	pan->between_s2 = bmod.sigma * bmod.sigma;
    } else {
	err = bmod.errcode;
#if PDEBUG
	fprintf(stderr, "error %d in between_variance\n", err);
#endif
    }

    if (!err && (pan->opt & OPT_B)) {
	err = save_between_model(&bmod, blist, ginfo, pan);
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
    int k = 0;

    for (i=0; i<pan->nbeta; i++) {
	if (vcv_skip(pmod, mi, pan, op)) {
	    i--;
	    mi++;
	    continue;
	}
	mj = mi;
	for (j=i; j<pan->nbeta; j++) {
	    if (vcv_skip(pmod, mj, pan, op)) {
		j--;
		mj++;
		continue;
	    }

	    idx = ijton(mi, mj, pmod->ncoeff);
	    if (op == VCV_SUBTRACT) {
		pan->sigma[k++] -= pmod->vcv[idx];
	    } else {
		pan->sigma[k++] = pmod->vcv[idx];
	    }
	    mj++;
	}
	mi++;
    }
}

/* calculate Hausman test statistic, matrix diff style */

static int bXb (panelmod_t *pan)
{
    char uplo = 'L'; 
    integer nrhs = 1;
    integer n = pan->nbeta;
    integer ldb = n;
    integer info = 0;
    integer *ipiv;
    double *x;
    int i;

    pan->H = NADBL;

    x = copyvec(pan->bdiff, pan->nbeta);
    if (x == NULL) {
	return E_ALLOC;
    }

    ipiv = malloc(pan->nbeta * sizeof *ipiv);
    if (ipiv == NULL) {
	free(x);
	return E_ALLOC;
    }

    /* solve for X-inverse * b */
    dspsv_(&uplo, &n, &nrhs, pan->sigma, ipiv, x, &ldb, &info);

    if (info > 0) {
	fprintf(stderr, "Hausman sigma matrix is singular\n");
    } else if (info < 0) {
	fprintf(stderr, "Illegal entry in Hausman sigma matrix\n");
    } else {
	pan->H = 0.0;
	for (i=0; i<pan->nbeta; i++) {
	    pan->H += x[i] * pan->bdiff[i];
	}
	if (pan->H < 0.0) {
	    pan->H = NADBL;
	}
    }

    free(x);
    free(ipiv);

    return info != 0;
}

static void panel_df_correction (MODEL *pmod, int k)
{
    double dfcorr = sqrt((double) pmod->dfd / (pmod->dfd - k));
    int i;

    pmod->dfd -= k;
#if FE_NEW
    pmod->dfn += k;
#else
    pmod->dfn += k - 1; 
#endif

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
			     DATAINFO *pdinfo,
			     PRN *prn)
{
    int dfn, i;

    pputs(prn, 
	  _("Fixed effects estimator\n"
	    "allows for differing intercepts by cross-sectional unit\n"
	    "slope standard errors in parentheses, p-values in brackets\n"));
    pputc(prn, '\n');

    /* print the slope coefficients, for varying regressors */
#if FE_NEW
    for (i=0; i<pmod->ncoeff; i++) {
	int vi = pan->vlist[i+2];

	print_panel_coeff(pmod, pdinfo->varname[vi], i, prn);
    }
#else
    for (i=0; i<pmod->ncoeff; i++) {
	int vi = pan->vlist[i+3];

	print_panel_coeff(pmod, pdinfo->varname[vi], i, prn);
    }
#endif
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
    pan->Fdfn = pan->N_fe - 1;
    pan->Fdfd = wmod->dfd;

    pan->F = (pan->pooled->ess - wmod->ess) * pan->Fdfd / 
	(wmod->ess * pan->Fdfn);
}

static int
fix_panelmod_list (MODEL *targ, panelmod_t *pan)
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
#if !FE_NEW
	/* and remove the const */
	gretl_list_delete_at_pos(targ->list, 2);
#endif
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
    int nc = fmod->ncoeff;
    int err = 0;

    err = fix_panelmod_list(fmod, pan);
    if (err) {
	return err;
    }

    fmod->ybar = pan->pooled->ybar;
    fmod->sdy = pan->pooled->sdy;
    fmod->tss = pan->pooled->tss;
    fmod->ifc = 1;

    gretl_model_set_double(fmod, "rsq_within", fmod->rsq);
    fmod->fstt = (fmod->rsq / (1.0 - fmod->rsq)) * 
	((double) fmod->dfd / (fmod->ncoeff - 1));
    gretl_model_set_double(fmod, "F_variables", fmod->fstt); /* ?? */

    /* should we modify R^2 in this way? */
    fmod->rsq = 1.0 - (fmod->ess / fmod->tss);

    if (fmod->dfd > 0) {
	double den = fmod->tss * fmod->dfd;

	fmod->adjrsq = 1 - (fmod->ess * (fmod->nobs - 1) / den);
    }

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

static int fe_model_add_ahat (MODEL *pmod, const double **Z, 
			      const DATAINFO *pdinfo,
			      panelmod_t *pan)
{
    double *ahat = NULL;
    double ahi;
    int i, j, t, bigt;
    int err = 0;

    ahat = malloc(pdinfo->n * sizeof *ahat);
    if (ahat == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	ahat[t] = NADBL;
    }

    for (i=0; i<pan->nunits; i++) {
	int Ti = pan->unit_obs[i];

	if (Ti == 0) {
	    continue;
	}

	/* a = y - Xb, where the 'b' is based on de-meaned data */

	ahi = 0.0;
	for (t=0; t<pan->T; t++) {
	    bigt = panel_index(i, t);
	    if (!na(pmod->uhat[bigt])) {
		ahi += Z[pmod->list[1]][bigt];
#if FE_NEW
		for (j=1; j<pmod->ncoeff; j++) {
		    ahi -= pmod->coeff[j] * Z[pmod->list[j+2]][bigt];
		}
#else
		for (j=0; j<pmod->ncoeff; j++) {
		    ahi -= pmod->coeff[j] * Z[pmod->list[j+2]][bigt];
		}
#endif
	    }
	}

	ahi /= Ti;

	for (t=0; t<pan->T; t++) {
	    bigt = panel_index(i, t);
	    if (!na(pmod->uhat[bigt])) {
		ahat[bigt] = ahi;
	    }
	}
    }

    err = gretl_model_set_data(pmod, "ahat", ahat, 
			       GRETL_TYPE_DOUBLE_ARRAY, 
			       pdinfo->n * sizeof *ahat);

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
fix_panel_hatvars (MODEL *pmod, const DATAINFO *dinfo, 
		   panelmod_t *pan, const double **Z)
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

static void verbose_femod_print (MODEL *femod, double **wZ,
				 DATAINFO *winfo, PRN *prn)
{
    int i, j;

    printmodel(femod, winfo, OPT_O, prn);

    fprintf(stderr, "femod: data series length = %d\n", winfo->n);
    for (i=0; i<winfo->n; i++) {
	fprintf(stderr, "femod.uhat[%d] = %g, ", i, femod->uhat[i]);
	fprintf(stderr, "data: ");
	for (j=0; j<winfo->v; j++) {
	    fprintf(stderr, "%g ", wZ[j][i]);
	}
	fputc('\n', stderr);
    }    
}

#endif

/* Estimate the fixed-effects model using a parallel dataset
   with the group means subtracted from all variables.
*/

static MODEL 
fixed_effects_model (panelmod_t *pan, const double **Z, 
		     DATAINFO *pdinfo, PRN *prn)
{
    MODEL femod;
    gretlopt lsqopt = OPT_A | OPT_Z | OPT_X;
    double **wZ = NULL;
    DATAINFO *winfo = NULL;
    int *felist = NULL;
    int i;

#if PDEBUG
    fprintf(stderr, "fixed_effects: using de-meaned data\n");
#endif

    gretl_model_init(&femod);

#if FE_NEW
    felist = gretl_list_new(pan->vlist[0]); 
#else
    felist = gretl_list_new(pan->vlist[0] - 1);
#endif
    if (felist == NULL) {
	femod.errcode = E_ALLOC;
	return femod;
    }	

    winfo = within_groups_dataset(Z, pdinfo, &wZ, pan);
    if (winfo == NULL) {
	free(felist);
	femod.errcode = E_ALLOC;
	return femod;
    } 

#if FE_NEW
    felist[1] = 1;
    felist[2] = 0;
    for (i=3; i<=felist[0]; i++) {
	felist[i] = i - 1;
    }
#else
    for (i=1; i<=felist[0]; i++) {
	felist[i] = i;
    }
#endif

    femod = lsq(felist, &wZ, winfo, OLS, lsqopt);

    if (femod.errcode) {
	; /* pass on */
    } else if (femod.list[0] < felist[0]) {
	/* one or more variables were dropped, because they were
	   all zero -- this is a symptom of collinearity */
	femod.errcode = E_SINGULAR;
    } else {
	/* we estimated a bunch of group means, and have to
	   subtract degrees of freedom */
#if FE_NEW
	panel_df_correction(&femod, pan->N_fe - 1);
#else
	panel_df_correction(&femod, pan->N_fe);
#endif
#if PDEBUG > 1
	verbose_femod_print(&femod, wZ, winfo, prn);
#endif
	if (pan->opt & OPT_F) {
	    /* estimating the FE model in its own right */
	    fix_panel_hatvars(&femod, winfo, pan, Z);
	    if (pan->opt & OPT_R) {
		panel_robust_vcv(&femod, pan, (const double **) wZ);
	    } else {
		femod_regular_vcv(&femod);
	    }
	}
    }

    destroy_dataset(wZ, winfo);
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
   at the stage of running the baseline pooled model (e.g. because of
   perfect collinearity).  In the case of fixed effects it may include
   additional variables dropped due to the fact that they are
   time-invariant.  We want to be able to show this list to the user
   when printing the model.
*/

static int compose_panel_droplist (MODEL *pmod, panelmod_t *pan)
{
    const int *pooldrop;
    int *dlist;
    int ndrop = 0;
    int j, vj, i = 1;

    if (gretl_model_get_int(pmod, "fixed-effects")) {
	ndrop = pan->pooled->list[0] - pan->vlist[0];
    }

    pooldrop = gretl_model_get_list(pan->pooled, "droplist");
    if (pooldrop != NULL) {
	ndrop += pooldrop[0];
    }

    if (ndrop == 0) {
	/* nothing to be done */
	return 0;
    }    

    dlist = gretl_list_new(ndrop);
    if (dlist == NULL) {
	return E_ALLOC;
    }

    if (pooldrop != NULL) {
	for (i=1; i<=pooldrop[0]; i++) {
	    dlist[i] = pooldrop[i];
	}
    } 

    if (pan->vlist[0] < pan->pooled->list[0]) {
	for (j=2; j<=pan->pooled->list[0]; j++) {
	    vj = pan->pooled->list[j];
	    if (!in_gretl_list(pan->vlist, vj)) {
		dlist[i++] = vj;
	    }
	}
    }

    return gretl_model_set_list_as_data(pmod, "droplist", dlist);
}

static void add_panel_obs_info (MODEL *pmod, panelmod_t *pan)
{
    gretl_model_set_int(pmod, "n_included_units", pan->effn);
    gretl_model_set_int(pmod, "Tmin", pan->Tmin);
    gretl_model_set_int(pmod, "Tmax", pan->Tmax);
}

/* spruce up the fixed-effects model and attach it to pan */

static int save_fixed_effects_model (MODEL *pmod, panelmod_t *pan,
				     const double **Z,
				     DATAINFO *pdinfo)
{
    int *ulist;
    int err = 0;

    pmod->ci = PANEL;
    gretl_model_set_int(pmod, "fixed-effects", 1);

    err = fix_within_stats(pmod, pan);
    if (err) {
	return err;
    }

    /* compose list of dropped variables, if any */
    compose_panel_droplist(pmod, pan);

    ulist = fe_units_list(pan);
    gretl_model_add_panel_varnames(pmod, pdinfo, ulist);
    free(ulist);

    fe_model_add_ahat(pmod, Z, pdinfo, pan);
    set_model_id(pmod);
    panel_dwstat(pmod, pan);
    save_fixed_effects_F(pan, pmod);
    time_dummies_wald_test(pan, pmod);

    add_panel_obs_info(pmod, pan);

    *pan->realmod = *pmod;

    return err;
}

/* drive the calculation of the fixed effects regression, print the
   results (if wanted), and compute the "within" error variance */

static int within_variance (panelmod_t *pan,
			    double ***pZ, DATAINFO *pdinfo,
			    PRN *prn)
{
    MODEL femod;
    int i, err = 0;

    femod = fixed_effects_model(pan, (const double **) *pZ, pdinfo, prn);

    if (femod.errcode) {
	pputs(prn, _("Error estimating fixed effects model\n"));
	errmsg(femod.errcode, prn);
	err = femod.errcode;
	clear_model(&femod);
    } else {
	/* Greene: nT - n - K */
	int den = femod.nobs - pan->effn - (pan->vlist[0] - 2);

	pan->s2e = femod.ess / den;
#if PDEBUG
	fprintf(stderr, "nT = %d, n = %d, K = %d\n", femod.nobs,
		pan->effn, pan->vlist[0] - 2);
	fprintf(stderr, "pan->s2e = %g / %d = %g\n", femod.ess,
	       den, pan->s2e);
#endif
	fixed_effects_F(pan, &femod);

	if (pan->opt & OPT_V) {
	    print_fe_results(pan, &femod, pdinfo, prn);
	}

	if (pan->bdiff != NULL && pan->sigma != NULL) {
#if FE_NEW
	    for (i=1; i<femod.ncoeff; i++) {
		pan->bdiff[i-1] = femod.coeff[i];
	    }
#else
	    for (i=0; i<femod.ncoeff; i++) {
		pan->bdiff[i] = femod.coeff[i];
	    }
#endif
	    femod_regular_vcv(&femod);
	    pan->sigma_e = femod.sigma;
	    vcv_slopes(pan, &femod, VCV_INIT);
	}

	if (pan->opt & OPT_F) {
	    err = save_fixed_effects_model(&femod, pan, 
					   (const double **) *pZ,
					   pdinfo);
	    if (err) {
		clear_model(&femod);
	    }
	} else {
	    clear_model(&femod);
	}
    }

    return err;
}

static void fix_gls_stats (MODEL *gmod, panelmod_t *pan)
{
    int nc;

    fix_panelmod_list(gmod, pan);

    gmod->ybar = pan->pooled->ybar;
    gmod->sdy = pan->pooled->sdy;
    gmod->tss = pan->pooled->tss;

    gmod->rsq = NADBL;
    gmod->adjrsq = NADBL;
    gmod->fstt = NADBL;
    gmod->rho = NADBL;
    gmod->dw = NADBL;

    nc = gmod->ncoeff;
    gmod->ncoeff = gmod->dfn + 1;
    ls_criteria(gmod);
    gmod->ncoeff = nc;
}

/* spruce up GLS model and attach it to pan */

static void save_random_effects_model (MODEL *pmod, panelmod_t *pan,
				       const double **Z,
				       const DATAINFO *reinfo)
{
    pmod->ci = PANEL;

    gretl_model_set_int(pmod, "random-effects", 1);
    gretl_model_set_double(pmod, "within-variance", pan->s2e);
    gretl_model_set_double(pmod, "between-variance", pan->between_s2);

    if (pan->balanced) {
	gretl_model_set_double(pmod, "gls-theta", pan->theta);
    }

    add_panel_obs_info(pmod, pan);

    /* compose list of dropped variables, if any */
    compose_panel_droplist(pmod, pan);

    gretl_model_add_panel_varnames(pmod, reinfo, NULL);
    fix_panel_hatvars(pmod, reinfo, pan, Z);
    fix_gls_stats(pmod, pan);
    set_model_id(pmod);

    *pan->realmod = *pmod;
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

    if (pan->realmod == NULL || na(pan->H)) {
	return;
    }

    test = model_test_new(GRETL_TEST_PANEL_HAUSMAN);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
	model_test_set_dfn(test, pan->nbeta);
	model_test_set_value(test, pan->H);
	model_test_set_pvalue(test, chisq_cdf_comp(pan->nbeta, pan->H));
	maybe_add_test_to_model(pan->realmod, test);
    }	    
}

/* Calculate the random effects regression.  Print the results
   here if we're doing the "panel diagnostics" test, otherwise
   save the results.
*/

static int random_effects (panelmod_t *pan, 
			   const double **Z, DATAINFO *pdinfo, 
			   const double **gZ, DATAINFO *ginfo, 
			   PRN *prn)
{
    double **reZ;
    DATAINFO *reinfo;
    MODEL remod;
    gretlopt lsqopt = OPT_A | OPT_Z;
    int *relist = NULL;
    int *hlist = NULL;
    double URSS = NADBL;
    int i, k, err = 0;

    gretl_model_init(&remod);

    /* GLS regression list */
    relist = gretl_list_new(pan->pooled->list[0]);
    if (relist == NULL) {
	return E_ALLOC;
    }

    if (pan->opt & OPT_H) {
	/* extra regressors for Hausman test, reg. approach */
	hlist = gretl_list_new(pan->vlist[0] - 1);
	if (hlist == NULL) {
	    free(relist);
	    return E_ALLOC;
	}
    }    

    /* Calculate the quasi-demeaning coefficient, theta, using the
       Swamy and Arora method.  Note: for unbalanced panels, theta
       will actually vary across the units in the final calculation.
    */
    pan->s2u = pan->between_s2 - pan->s2e / pan->Tbar;
    if (pan->s2u < 0) {
	pan->s2u = 0.0;
    }
    pan->theta = 1.0 - sqrt(pan->s2e / (pan->Tmax * pan->s2u + pan->s2e));

#if PDEBUG
    fprintf(stderr, "s_u = %.8g, s_e = %.8g\n", sqrt(pan->s2u), sqrt(pan->s2e));
    fprintf(stderr, "random_effects theta = %g\n", pan->theta);
#endif

    /* make special transformed dataset, and regression list */
    reinfo = random_effects_dataset(Z, pdinfo, gZ, ginfo, &reZ, relist, 
				    hlist, pan);
    if (reinfo == NULL) {
	free(relist);
	free(hlist);
	return E_ALLOC;
    }

    if (hlist != NULL) {
	/* estimate the augmented random-effects model (GLS) */
	int *biglist = gretl_list_add(relist, hlist, &err);

	if (!err) {
	    remod = lsq(biglist, &reZ, reinfo, OLS, OPT_A);
	    if (remod.errcode == 0) {
		URSS = remod.ess;
	    }
	    clear_model(&remod);
	}
	free(biglist);
    }	

    /* regular random-effects model */
    remod = lsq(relist, &reZ, reinfo, OLS, lsqopt);

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
		print_panel_coeff(&remod, pdinfo->varname[vi], i, prn);
	    }
	    if (pan->bdiff != NULL && var_is_varying(pan, vi)) {
		pan->bdiff[k++] -= remod.coeff[i];
	    }
	}

	if (!na(URSS)) {
	    /* it appears that Baltagi uses T-k instead of T here. Why? */
	    pan->H = (remod.ess / URSS - 1.0) * (remod.nobs);
	} else if (pan->sigma != NULL) {
	    makevcv(&remod, remod.sigma);
	    vcv_slopes(pan, &remod, VCV_SUBTRACT);
	}
    }

#if PDEBUG > 1
    if (remod.errcode == 0) {
	printmodel(&remod, reinfo, OPT_NONE, prn);
    }
#endif	

    if (!err && (pan->opt & OPT_U)) {
	save_random_effects_model(&remod, pan, Z, reinfo);
    } else {
	clear_model(&remod);
    }

    destroy_dataset(reZ, reinfo);

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

static int 
breusch_pagan_LM (panelmod_t *pan, const DATAINFO *pdinfo, PRN *prn)
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

    pan->BP = (n * n /(2.0 * (M - n))) * pow((A / pan->pooled->ess) - 1.0, 2);

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
    int err = 0;

    if (pan->bdiff != NULL && pan->sigma != NULL) {
	/* matrix approach */
	err = bXb(pan);
    } else if (na(pan->H)) {
	/* regression approach bombed somehow? */
	err = 1;
    }

    if (pan->opt & OPT_V) {
	if (err) {
	    pputs(prn, _("Error attempting to invert vcv difference matrix\n"));
	} else {
	    print_hausman_result(pan, prn);
	}
    } else if (!err) {
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
    pan->N_fe = 0;
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
	if (uobs[i] > 0) {
	    pan->N_fe += 1;
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

/* find harmomic mean of the number of time-series observations
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
panelmod_setup (panelmod_t *pan, MODEL *pmod, const DATAINFO *pdinfo, 
		int ntdum, gretlopt opt)
{
    int err = 0;

    pan->opt = opt;
    pan->pooled = pmod;

    /* assumes (possibly padded) balanced panel dataset */
    pan->nunits = (pdinfo->t2 - pdinfo->t1 + 1) / pdinfo->pd;
    pan->T = pdinfo->pd;

    panel_index_init(pdinfo, pan->nunits, pan->T);
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

int panel_diagnostics (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
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
    err = panelmod_setup(&pan, pmod, pdinfo, 0, opt | OPT_V);
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
    err = varying_vars_list((const double **) *pZ, pdinfo, &pan);
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

    err = within_variance(&pan, pZ, pdinfo, prn);
    if (err) {
	goto bailout;
    }

    breusch_pagan_LM(&pan, pdinfo, prn);

    if (xdf <= 0) {
	pprintf(prn, "Omitting group means regression: "
		"insufficient degrees of freedom\n");
	goto bailout;
    }
    
    if (xdf > 0 && !na(pan.s2e)) {
	double **gZ = NULL;
	DATAINFO *ginfo;

	ginfo = group_means_dataset(&pan, (const double **) *pZ,
				    pdinfo, &gZ);
	if (ginfo == NULL) {
	    err = E_ALLOC;
	} else {
	    err = between_variance(&pan, &gZ, ginfo);
	}

	if (err) { 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	    if (err == E_SINGULAR) {
		/* treat this as non-fatal */
		err = 0;
	    }
	} else {
	    random_effects(&pan, (const double **) *pZ, pdinfo, 
			   (const double **) gZ, ginfo, prn);
	    finalize_hausman_test(&pan, prn);
	}

	if (ginfo != NULL) {
	    destroy_dataset(gZ, ginfo);
	}
    }

 bailout:

    panelmod_free(&pan);

    return err;
}

static int between_model (panelmod_t *pan, const double **Z,
			  const DATAINFO *pdinfo)
{
    double **gZ = NULL;
    DATAINFO *ginfo;
    int err = 0;

    ginfo = group_means_dataset(pan, Z, pdinfo, &gZ);
    if (ginfo == NULL) {
	err = E_ALLOC;
    } else {
	err = between_variance(pan, &gZ, ginfo);
    }

    if (ginfo != NULL) {
	destroy_dataset(gZ, ginfo);
    }

    return err;
}

static int
add_dummies_to_list (const int *list, DATAINFO *pdinfo, int **pbiglist)
{
    char dname[VNAMELEN];
    int *biglist = NULL;
    int i, j, v;
    int err = 0;

    biglist = gretl_list_new(list[0] + pdinfo->pd - 1);
    if (biglist == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
	biglist[i] = list[i];
    }

    j = list[0] + 1;
    for (i=2; i<=pdinfo->pd; i++) {
	sprintf(dname, "dt_%d", i);
	v = varindex(pdinfo, dname);
	if (v == pdinfo->v) {
	    err = E_DATA;
	    break;
	} else {
	    biglist[j++] = v;
	}
    }

    if (err) {
	free(biglist);
    } else {
	*pbiglist = biglist;
    }

    return err;
}

static int panel_check_for_const (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
	if (list[i] == 0) {
	    return 0;
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
    set_model_id(pmod);

    if (pan->opt & OPT_R) {
	panel_robust_vcv(pmod, pan, Z);
    }
}

/* fixed effects | random effects | between | pooled */

#define estimator_specified(o) (o & (OPT_F|OPT_U|OPT_B|OPT_P))

/* real_panel_model:
 * @list: list containing model specification.
 * @pZ: pointer to data array.  
 * @pdinfo: data info pointer. 
 * @opt: may include %OPT_U for the random effects model;
 * %OPT_R for robust standard errors (fixed effects model
 * and pooled OLS only); %OPT_H to use the regression approach 
 * to the Hausman test (random effects only); %OPT_B for 
 * the "between" model; %OPT_P for pooled OLS; and %OPT_D to 
 * include time dummies.
 * @prn: printing struct.
 *
 * Estimates a panel model, by default the fixed effects model.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL real_panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn)
{
    MODEL mod;
    panelmod_t pan;
    gretlopt pan_opt = opt;
    int *olslist = NULL;
    int orig_v = pdinfo->v;
    int ntdum = 0;
    int err = 0;

    gretl_model_init(&mod);
    panelmod_init(&pan);

    if (!(opt & OPT_P)) {
	mod.errcode = panel_check_for_const(list);
	if (mod.errcode) {
	    return mod;
	}
    }

    /* add time dummies to list? */
    if (opt & OPT_D) {
	err = panel_dummies(pZ, pdinfo, OPT_T);
	if (!err) {
	    err = add_dummies_to_list(list, pdinfo, &olslist);
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
    mod = lsq(olslist, pZ, pdinfo, OLS, OPT_A);
    if (mod.errcode) {
	err = mod.errcode;
	fprintf(stderr, "real_panel_model: error %d in intial OLS\n", mod.errcode);
	goto bailout;
    }

    free(olslist);

#if PDEBUG
    printmodel(&mod, pdinfo, OPT_NONE, prn);
#endif

    if (!estimator_specified(opt)) {
	/* default: add OPT_F to save the fixed effects model */
	pan_opt |= OPT_F;
    }

    if ((opt & OPT_D) && !(opt & OPT_B)) {
	ntdum = get_ntdum(list, mod.list);
    }

    err = panelmod_setup(&pan, &mod, pdinfo, ntdum, pan_opt);
    if (err) {
	goto bailout;
    }   

    if (opt & OPT_P) {
	save_pooled_model(&mod, &pan, (const double **) *pZ);
	goto bailout;
    }

    if (opt & OPT_B) {
	err = between_model(&pan, (const double **) *pZ, pdinfo);
	goto bailout;
    }	

    calculate_Tbar(&pan);

    /* figure out which of the original regressors are time-varying,
       or unit-varying as the case may be 
    */
    err = varying_vars_list((const double **) *pZ, pdinfo, &pan);
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
	    mod.errcode = E_DF;
	    goto bailout;
	} else {
	    err = hausman_allocate(&pan);
	    if (err) {
		goto bailout;
	    }
	}
    } 

    err = within_variance(&pan, pZ, pdinfo, prn);
    if (err) {
	goto bailout;
    }

    if ((opt & OPT_U) && !na(pan.s2e)) {
	double **gZ = NULL;
	DATAINFO *ginfo;

	breusch_pagan_LM(&pan, pdinfo, prn);

	ginfo = group_means_dataset(&pan, (const double **) *pZ,
				    pdinfo, &gZ);
	if (ginfo == NULL) {
	    err = E_ALLOC;
	} else {
	    err = between_variance(&pan, &gZ, ginfo);
	}

	if (err) { 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	} else {
	    random_effects(&pan, (const double **) *pZ, pdinfo, 
			   (const double **) gZ, ginfo, prn);
	    save_breusch_pagan_result(&pan);
	    finalize_hausman_test(&pan, prn);
	}

	if (ginfo != NULL) {
	    destroy_dataset(gZ, ginfo);
	}
    }

 bailout:

    if (!err) {
	if (!(opt & OPT_P)) {
	    clear_model(&mod);
	    mod = *pan.realmod;
	}
	gretl_model_smpl_init(&mod, pdinfo);
    }

    panelmod_free(&pan);

    if (err && mod.errcode == 0) {
	mod.errcode = err;
    }

#if 1
    dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
#endif

    return mod;    
}

int panel_tsls_robust_vcv (MODEL *pmod, const double **Z, 
			   const DATAINFO *pdinfo)
{
    panelmod_t pan;
    int err = 0;

    panelmod_init(&pan);

    err = panelmod_setup(&pan, pmod, pdinfo, 0, OPT_NONE);
    if (!err) {
	err = panel_robust_vcv(pmod, &pan, Z);
    }

    panelmod_free(&pan); 

    return err;
}

/* write weights for groupwise weighted least squares into the
   last variable in the dataset */

static int 
write_weights_to_dataset (double *uvar, int nunits, int T,
			  double **Z, DATAINFO *pdinfo)
{
    int w = pdinfo->v - 1;
    double wi;
    int i, t;

    for (i=0; i<nunits; i++) {
	if (uvar[i] <= 0.0 || na(uvar[i])) {
	    wi = 0.0;
	} else {
	    wi = 1.0 / uvar[i];
	}
	for (t=0; t<T; t++) {
	    Z[w][panel_index(i, t)] = wi;
	}
    }

    return 0;
}

static int
allocate_weight_var (double ***pZ, DATAINFO *pdinfo)
{
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    strcpy(pdinfo->varname[pdinfo->v - 1], "unit_wt");

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

#define S2MINOBS 2

static void
print_wald_test (double W, int nunits, const int *unit_obs, PRN *prn)
{
    int i, df = 0;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] >= S2MINOBS) df++;
    }

    pprintf(prn, "\n%s\n%s:\n",
	    _("Distribution free Wald test for heteroskedasticity"),
	    _("based on the FGLS residuals"));
    pprintf(prn, "%s(%d) = %g, ",  _("Chi-square"), df, W);
    pprintf(prn, _("with p-value = %g\n\n"), chisq_cdf_comp(df, W));
}

/* Wald test for groupwise heteroskedasticity, without assuming
   normality of the errors: see Greene, 4e, p. 598.  Note that the
   computation involves a factor of (1/(T_i - 1)) so this is not
   usable when some of the groups have only one observation.
*/

static double 
wald_hetero_test (const MODEL *pmod, const DATAINFO *pdinfo, 
		  double s2, const double *uvar,
		  panelmod_t *pan)
{
    double x, W = 0.0;
    int i, t, Ti;

    for (i=0; i<pan->nunits; i++) {
	double Vi = 0.0;

	Ti = pan->unit_obs[i];
	if (Ti == 0) {
	    continue;
	}

	if (Ti < S2MINOBS) {
	    W = NADBL;
	    break;
	}

	for (t=0; t<pan->T; t++) {
	    x = pmod->uhat[panel_index(i, t)];
	    if (!na(x)) {
		x = x * x - uvar[i];
		Vi += x * x;
	    }
	}

	if (Vi <= 0) {
	    W = NADBL;
	    break;
	}
	    
	Vi *= (1.0 / Ti) * (1.0 / (Ti - 1.0));
	x = uvar[i] - s2;
	W += x * x / Vi;
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
   observation? */

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

static void 
unit_error_variances (double *uvar, const MODEL *pmod, 
		      const DATAINFO *pdinfo,
		      panelmod_t *pan)
{
    int i, t;
    double uit;

    for (i=0; i<pan->nunits; i++) {
	if (pan->unit_obs[i] == 0) {
	    uvar[i] = NADBL;
	    continue;
	}
	uvar[i] = 0.0;
	for (t=0; t<pan->T; t++) {
	    uit = pmod->uhat[panel_index(i, t)];
	    if (!na(uit)) {
		uvar[i] += uit * uit;
	    }
	}
	uvar[i] /= pan->unit_obs[i]; 
    }
}

#define SMALLDIFF 0.0001
#define WLS_MAX   30

MODEL panel_wls_by_unit (const int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn)
{
    MODEL mdl;
    panelmod_t pan;
    gretlopt wlsopt = OPT_A;
    double *uvar = NULL;
    double *bvec = NULL;
    double s2, diff = 1.0;
    double W = NADBL;
    int *wlist = NULL;
    int orig_v = pdinfo->v;
    int i, iter = 0;

    gretl_error_clear();
    panelmod_init(&pan);

    if (opt & OPT_T) {
	/* iterating: no degrees-of-freedom correction */
	wlsopt |= OPT_N; 
    } 

    if (opt & OPT_A) {
	wlsopt |= OPT_A;
    }

    /* baseline pooled model */
    mdl = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (mdl.errcode) {
	goto bailout;
    }

    mdl.errcode = panelmod_setup(&pan, &mdl, pdinfo, 0, OPT_NONE);
    if (mdl.errcode) {
	goto bailout;
    }

    uvar = malloc(pan.nunits * sizeof *uvar);
    if (uvar == NULL) {
	free(pan.unit_obs);
	mdl.errcode = E_ALLOC;
	return mdl;
    }  

    if (opt & OPT_T) {
	if (singleton_check(pan.unit_obs, pan.nunits)) {
	    pprintf(prn, _("Can't produce ML estimates: "
			   "some units have only one observation"));
	    pputc(prn, '\n');
	    opt ^= OPT_T;
	}
    }

    s2 = mdl.ess / mdl.nobs;

    if ((opt & OPT_V) && (opt & OPT_T)) {
	pprintf(prn, "\nOLS error variance = %g\n", s2);
	pprintf(prn, "log-likelihood = %g\n", pooled_ll(&mdl));
    }

    if (allocate_weight_var(pZ, pdinfo)) {
	mdl.errcode = E_ALLOC;
	goto bailout;
    }

    if (opt & OPT_T) {
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

    wlist[1] = pdinfo->v - 1; /* weight variable: the last var added */
    for (i=2; i<=wlist[0]; i++) {
	wlist[i] = mdl.list[i-1];
    }

    /* if wanted (and possible) iterate to ML solution; otherwise just do
       one-step FGLS estimation 
    */

    while (diff > SMALLDIFF) {

	iter++;

	unit_error_variances(uvar, &mdl, pdinfo, &pan);

	if (opt & OPT_V) {
	    if (opt & OPT_T) {
		pprintf(prn, "\n*** %s %d ***\n", _("iteration"), 
			iter);
	    } else {
		pputc(prn, '\n');
	    }
	    pputs(prn, " unit    variance\n");
	    for (i=0; i<pan.nunits; i++) {
		if (pan.unit_obs[i] > 0) {
		    pprintf(prn, "%5d%12g (T = %d)\n", i + 1, uvar[i], pan.unit_obs[i]);
		} else {
		    pprintf(prn, "%5d%12s (T = %d)\n", i + 1, "NA", pan.unit_obs[i]);
		}
	    }
	}

	if ((opt & OPT_T) && iter == 2) {
	    W = wald_hetero_test(&mdl, pdinfo, s2, uvar, &pan);
	}

	write_weights_to_dataset(uvar, pan.nunits, pan.T, *pZ, pdinfo);

	if (opt & OPT_T) {
	    /* save coefficients for comparison */
	    for (i=0; i<mdl.ncoeff; i++) {
		bvec[i] = mdl.coeff[i];
	    }
	}

	clear_model(&mdl);
	mdl = lsq(wlist, pZ, pdinfo, WLS, wlsopt);
	if (mdl.errcode) {
	    break;
	}

	if (iter > WLS_MAX) {
	    mdl.errcode = E_NOCONV;
	    break;
	}

	if (opt & OPT_T) {
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
	gretl_model_set_int(&mdl, "unit-weights", 1);
	mdl.nwt = 0;

	if (opt & OPT_T) {
	    gretl_model_set_int(&mdl, "iters", iter);
	    ml_hetero_test(&mdl, s2, uvar, pan.nunits, pan.unit_obs);
	    unit_error_variances(uvar, &mdl, pdinfo, &pan);
	    panel_ML_ll(&mdl, uvar, pan.nunits, pan.unit_obs);
	    if (opt & OPT_V) {
		pputc(prn, '\n');
	    }
	} else {
	    unit_error_variances(uvar, &mdl, pdinfo, &pan);
	    W = wald_hetero_test(&mdl, pdinfo, s2, uvar, &pan);
	}

	if (!na(W)) {
	    print_wald_test(W, pan.nunits, pan.unit_obs, prn);
	}
    }    

 bailout:

    free(pan.unit_obs);
    free(uvar);
    free(wlist);
    free(bvec);

    if (!(opt & OPT_V)) {
	dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
    }
    
    return mdl;
}

static void panel_copy_var (double **targZ, DATAINFO *targinfo, int targv,
			    double *src, DATAINFO *srcinfo, int srcv,
			    int order)
{
    int t, j = 0;

    for (t=srcinfo->t1; t<=srcinfo->t2; t++) {
	if (t % srcinfo->pd >= order) {
	    targZ[targv][j++] = src[t];
	} 
    }

    if (srcv == -1) {
	strcpy(targinfo->varname[targv], "uhat");
	strcpy(VARLABEL(targinfo, targv), _("residual"));
    } else {
	strcpy(targinfo->varname[targv], srcinfo->varname[srcv]);
	strcpy(VARLABEL(targinfo, targv), VARLABEL(srcinfo, srcv));
    }
}

static void make_reduced_data_info (DATAINFO *targ, DATAINFO *src, int order)
{
    targ->pd = src->pd - order;
    ntodate(targ->stobs, src->t1 + order, src);
    targ->sd0 = obs_str_to_double(targ->stobs); 
    targ->structure = src->structure;
}

static void panel_lag (double **tmpZ, DATAINFO *tmpinfo, 
		       double *src, DATAINFO *srcinfo, 
		       int v, int order, int lag)
{
    int t, j = 0;

    for (t=srcinfo->t1; t<=srcinfo->t2; t++) {
	if (t % srcinfo->pd >= order) {
	    tmpZ[v][j++] = src[t - lag];
	}
    }

    sprintf(tmpinfo->varname[v], "uhat_%d", lag);
    *VARLABEL(tmpinfo, v) = 0;
}

/* - do some sanity checks
   - create a local copy of the required portion of the data set,
     skipping the obs that will be missing
   - copy in the lags of uhat
   - estimate the aux model
   - destroy the temporary data set
*/

int panel_autocorr_test (MODEL *pmod, int order, 
			 double **Z, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn)
{
    int *aclist;
    double **tmpZ;
    DATAINFO *tmpinfo;
    MODEL aux;
    double trsq, LMF;
    int i, nv, nunits, nobs, err = 0;
    int sn = pdinfo->t2 - pdinfo->t1 + 1;

    if (pmod->ci != OLS) {
	return E_NOTIMP;
    }

    if (pmod->missmask != NULL) {
	return E_DATA;
    }

    /* basic checks */
    if (order <= 0) order = 1;
    if (order > pdinfo->pd - 1) return E_DF;
    if (pmod->ncoeff + order >= sn) return E_DF;

    /* get number of cross-sectional units */
    nunits = sn / pdinfo->pd;

    /* we lose "order" observations for each unit */
    nobs = sn - nunits * order;

    /* the required number of variables */
    nv = pmod->list[0] + order;

    /* create temporary reduced dataset */
    tmpinfo = create_auxiliary_dataset(&tmpZ, nv, nobs);
    if (tmpinfo == NULL) {
	return E_ALLOC;
    }

    make_reduced_data_info(tmpinfo, pdinfo, order);

#if PDEBUG
    fprintf(stderr, "Created data set, n=%d, pd=%d, vars=%d, stobs='%s'\n", 
	    tmpinfo->n, tmpinfo->pd, tmpinfo->v, tmpinfo->stobs);
#endif

    /* allocate the auxiliary regression list */
    aclist = malloc((nv + 1) * sizeof *aclist);
    if (aclist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	int k;

	aclist[0] = pmod->list[0] + order;
	/* copy model uhat to position 1 in temp data set */
	aclist[1] = 1;
	panel_copy_var(tmpZ, tmpinfo, 1,
		       &pmod->uhat[0], pdinfo, -1,
		       order);
	/* copy across the original indep vars, making
	   the new regression list while we're at it */
	k = 2;
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == 0) { /* the constant */
		aclist[i] = 0;
	    } else {
		aclist[i] = k;
		panel_copy_var(tmpZ, tmpinfo, k, 
			       &Z[pmod->list[i]][0], pdinfo, pmod->list[i], 
			       order);
		k++;
	    }
	}
    }

    if (!err) {
	int v = pmod->list[0] - 1;

	/* add lags of uhat to temp data set */
	for (i=1; i<=order; i++) {
	    panel_lag(tmpZ, tmpinfo, &pmod->uhat[0], 
		      pdinfo, v + i, order, i);
	    aclist[v + i + 1] = v + i;
	}
    }

    if (!err) {
	aux = lsq(aclist, &tmpZ, tmpinfo, OLS, OPT_A);
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
	printmodel(&aux, tmpinfo, OPT_NONE, prn);
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

    destroy_dataset(tmpZ, tmpinfo);

    return err;
}

/**
 * switch_panel_orientation:
 * @Z: data array.
 * @pdinfo: dataset information struct.
 * 
 * Reorganizes the data array @Z and rewrites the dataset information
 * @pdinfo, transforming from stacked cross sections to stacked
 * time series or vice versa.  If the transformation is from
 * stacked time series to stacked cross section, the dataset will
 * no longer be acceptable as a panel for gretl's purposes; it
 * may be useful for export purposes, though.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int switch_panel_orientation (double **Z, DATAINFO *pdinfo)
{
    int oldmode = pdinfo->structure;
    double *tmp;
    char **markers = NULL;
    double pdx;
    int T, n;
    int i, j, s, t;

    if (oldmode != STACKED_TIME_SERIES &&
	oldmode != STACKED_CROSS_SECTION) {
	return E_DATA;
    }

    tmp = malloc(pdinfo->n * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    if (oldmode == STACKED_CROSS_SECTION) {
	n = pdinfo->pd;
	T = pdinfo->n / n;
    } else {
	T = pdinfo->pd;
	n = pdinfo->n / T;
    }

    /* copy the data series across in transformed order */
    for (i=1; i<pdinfo->v; i++) {
	if (var_is_scalar(pdinfo, i)) {
	    continue;
	}
	for (t=0; t<pdinfo->n; t++) {
	    /* transcribe to tmp in original order */
	    tmp[t] = Z[i][t];
	}
	s = 0;
	if (oldmode == STACKED_CROSS_SECTION) {
	    /* convert to stacked time-series */
	    for (j=0; j<n; j++) {
		for (t=0; t<T; t++) {
		    Z[i][s++] = tmp[t * n + j];
		}
	    }
	} else {
	    /* convert to stacked cross-sections */
	    for (t=0; t<T; t++) {
		for (j=0; j<n; j++) {
		    Z[i][s++] = tmp[j * T + t];
		}
	    } 
	}  
    }

    /* rearrange observations markers if relevant */
    if (pdinfo->S != NULL) {
	markers = strings_array_new_with_length(pdinfo->n, OBSLEN);
	if (markers != NULL) {
	    for (t=0; t<pdinfo->n; t++) {
		strcpy(markers[t], pdinfo->S[t]);
	    }
	    s = 0;
	    if (oldmode == STACKED_CROSS_SECTION) {
		for (j=0; j<n; j++) {
		    for (t=0; t<T; t++) {
			strcpy(pdinfo->S[s++], markers[t * n + j]);
		    }
		}
	    } else {
		for (t=0; t<T; t++) {
		    for (j=0; j<n; j++) {
			strcpy(pdinfo->S[s++], markers[j * T + t]);
		    }
		}
	    }
	    free_strings_array(markers, pdinfo->n);
	} else {
	    /* should we flag an error? */
	    dataset_destroy_obs_markers(pdinfo); 
	}
    }

    pdinfo->sd0 = 1.0;
    pdx = 0.1;

    /* change the datainfo setup */
    if (oldmode == STACKED_CROSS_SECTION) {
	pdinfo->structure = STACKED_TIME_SERIES;
	pdinfo->pd = T;
	while (T /= 10) {
	    pdx *= 0.1;
	}	
    } else {
	pdinfo->structure = STACKED_CROSS_SECTION;
	pdinfo->pd = n;
	while (n /= 10) {
	    pdx *= 0.1;
	}
    }
	
    pdinfo->sd0 += pdx;
    ntodate(pdinfo->stobs, 0, pdinfo);
    ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);

    free(tmp);

    return 0;
}

static int 
varying_vars_list (const double **Z, const DATAINFO *pdinfo,
		   panelmod_t *pan)
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
		    xval = Z[vj][bigt];
		    started = 1;
		} else if (Z[vj][bigt] != xval) {
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
		    vj, pdinfo->varname[vj]);
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
    const s_point *pa = (const s_point *) a;
    const s_point *pb = (const s_point *) b;
    int ret;

    ret = (pa->val1 > pb->val1) - (pa->val1 < pb->val1);
    if (ret == 0) {
	ret = (pa->val2 > pb->val2) - (pa->val2 < pb->val2);
    }

    return ret;
}

static int compose_panel_indices (DATAINFO *pdinfo, int T,
				  char *mask)
{
    int i, u, t;
    int err;

    err = dataset_allocate_panel_info(pdinfo);
    if (err) {
	return err;
    }
    
    pdinfo->paninfo->padmask = mask;

    u = t = 0;
    for (i=0; i<pdinfo->n; i++) {
	if (i > 0 && i % T == 0) {
	    t = 0;
	    u++;
	} 
	pdinfo->paninfo->unit[i] = u;
	pdinfo->paninfo->period[i] = t++;
    }

#if PDEBUG
    fprintf(stderr, "transcribe_panel_indices:\n");
    for (i=0; i<pdinfo->n; i++) {
	fprintf(stderr, " i=%d, unit=%d, period=%d\n", i, 
		pdinfo->paninfo->unit[i],
		pdinfo->paninfo->period[i]);
    }
#endif

    err = dataset_finalize_panel_indices(pdinfo);

    return err;
}

static int panel_data_sort_by (double **Z, DATAINFO *pdinfo,
			       int uv, int tv, int *ustrs)
{
    int n = pdinfo->n;
    char **S = NULL;
    double *tmp = NULL;
    int i, t;
    sorter s;

    tmp = malloc(n * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    s.points = malloc(n * sizeof *s.points);
    if (s.points == NULL) {
	free(tmp);
	return E_ALLOC;
    } 

    if (pdinfo->S != NULL) {
	S = strings_array_new_with_length(n, OBSLEN);
	if (S == NULL) {
	    free(tmp);
	    free(s.points);
	    return E_ALLOC;
	}
    }

    s.sortvar1 = uv;
    s.sortvar2 = tv;
    
    sorter_fill_points(&s, (const double **) Z, n);

    qsort((void *) s.points, (size_t) n, sizeof s.points[0], 
	  compare_obs);

    for (i=1; i<pdinfo->v; i++) {
	if (var_is_series(pdinfo, i)) {
	    for (t=0; t<n; t++) {
		tmp[t] = Z[i][t];
	    }
	    for (t=0; t<n; t++) {
		Z[i][t] = tmp[s.points[t].obsnum];
	    }
	}	
    }

    if (S != NULL) {
	*ustrs = 1;
	for (t=0; t<n; t++) {
	    strcpy(S[t], pdinfo->S[t]);
	    if (S[t][0] == '\0') {
		*ustrs = 0;
	    }
	}
	for (t=0; t<n; t++) {
	    strcpy(pdinfo->S[t], S[s.points[t].obsnum]);
	}
	for (t=1; t<n && *ustrs; t++) {
	    if (Z[uv][t] == Z[uv][t-1] &&
		strcmp(pdinfo->S[t], pdinfo->S[t-1])) {
		*ustrs = 0;
	    }
	}
	free_strings_array(S, n);
    }	

    free(s.points);
    free(tmp);

    return 0;
}

/* Given the variables coding for panel unit and time period,
   construct parallel arrays of int that code for the same
   information, but are zero-based and consecutive.
*/

static int normalize_uid_tid (const double *tid, int T,
			      const double **Z, const DATAINFO *pdinfo,
			      int uv, int tv, int **pnuid, int **pntid)
{
    int *nuid = NULL;
    int *ntid = NULL;
    int ui, tmin;
    int i, t;

    nuid = malloc(pdinfo->n * sizeof *nuid);
    ntid = malloc(pdinfo->n * sizeof *ntid);
    if (nuid == NULL || ntid == NULL) {
	free(nuid);
	free(ntid);
	return E_ALLOC;
    }

    ui = tmin = 0;

    for (i=0; i<pdinfo->n; i++) {
	if (i > 0 && Z[uv][i] > Z[uv][i-1]) {
	    tmin = 0;
	    ui++;
	} else if (i > 0) {
	    tmin++;
	}
	nuid[i] = ui;
	for (t=tmin; t<T; t++) {
	    if (Z[tv][i] == tid[t]) {
		ntid[i] = t;
		break;
	    }
	}
    }

    *pnuid = nuid;
    *pntid = ntid;

    return 0;
}

static int pad_panel_dataset (const double *uid, int uv, int nunits,
			      const double *tid, int tv, int nperiods, 
			      double **Z, DATAINFO *pdinfo,
			      int ustrs, char *mask)
{
    double **bigZ = NULL;
    char **S = NULL;
    int *nuid = NULL;
    int *ntid = NULL;
    int n_scalars = 0;
    int n_orig = pdinfo->n;
    int t2_orig = pdinfo->t2;
    int i, j, s, t;
    int err = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (var_is_scalar(pdinfo, i)) {
	    n_scalars++;
	}
    }

    err = normalize_uid_tid(tid, nperiods, (const double **) Z,
			    pdinfo, uv, tv, &nuid, &ntid);

    if (!err) {
	/* allocate temporary storage, skipping any scalars */
	pdinfo->v -= n_scalars;
	pdinfo->n = nunits * nperiods;
	pdinfo->t2 = pdinfo->n - 1;
	err = allocate_Z(&bigZ, pdinfo);
	pdinfo->v += n_scalars;
    }

    if (!err && pdinfo->S != NULL && ustrs) {
	S = strings_array_new_with_length(pdinfo->n, OBSLEN);
    }

    if (err) {
	pdinfo->n = n_orig;
	pdinfo->t2 = t2_orig;
    } else {
	int buv = 0, btv = 0;
	int tref = 0;

	/* write rows from original Z into the right places in bigZ */
	for (t=0; t<n_orig; t++) {
	    j = 1;
	    s = nuid[t] * nperiods + ntid[t];
	    for (i=1; i<pdinfo->v; i++) {
		if (var_is_series(pdinfo, i)) {
		    bigZ[j++][s] = Z[i][t];
		}
	    }
	    if (S != NULL) {
		strcpy(S[s], pdinfo->S[t]);
	    }
	    if (mask != NULL && s > tref) {
		/* recording the padding in "mask" */
		for (j=tref; j<s; j++) {
		    mask[j] = 1;
		}
		tref = s;
	    }
	    tref++;
	}

	/* where are uv and tv in bigZ? */
	j = 1;
	for (i=1; i<pdinfo->v; i++) {
	    if (i == uv) {
		buv = j;
	    } else if (i == tv) {
		btv = j;
	    }
	    if (var_is_series(pdinfo, i)) {
		j++;
	    }
	}

	/* complete the index info in slots uv and tv */
	i = j = 0;
	for (t=0; t<pdinfo->n; t++) {
	    if (t > 0 && t % nperiods == 0) {
		i++;
		j = 0;
	    }
	    bigZ[buv][t] = uid[i];
	    bigZ[btv][t] = tid[j++];
	}

	/* swap the padded arrays into Z */
	for (i=0, j=0; i<pdinfo->v; i++) {
	    if (var_is_series(pdinfo, i)) {
		free(Z[i]);
		Z[i] = bigZ[j++];
	    }
	}

	free(bigZ);
    }

    free(nuid);
    free(ntid);

    if (pdinfo->S != NULL) {
	/* expand the obs (unit) marker strings appropriately */
	if (S == NULL) {
	    dataset_destroy_obs_markers(pdinfo);
	} else {
	    char si[OBSLEN];

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
	    free_strings_array(pdinfo->S, n_orig);
	    pdinfo->S = S;
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

static int uv_tv_from_line (const char *line, const DATAINFO *pdinfo,
			    int *uv, int *tv)
{
    char uvname[VNAMELEN];
    char tvname[VNAMELEN];
    int err = 0;
    
    if (sscanf(line, "%15s %15s", uvname, tvname) != 2) {
	return E_PARSE;
    }

    *uv = varindex(pdinfo, uvname);
    if (*uv == pdinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), uvname);
	err = E_UNKVAR;
    } else if (!var_is_series(pdinfo, *uv)) {
	err = E_DATATYPE;
    }

    if (err) {
	return err;
    }

    *tv = varindex(pdinfo, tvname);
    if (*tv == pdinfo->v) {
	sprintf(gretl_errmsg, _("Unknown variable '%s'"), tvname);
	err = E_UNKVAR;
    } else if (!var_is_series(pdinfo, *tv)) {
	err = E_DATATYPE;
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

static void finalize_panel_datainfo (DATAINFO *pdinfo, int nperiods)
{
    int pdp = nperiods;
    int den = 10.0;

    while ((pdp = pdp / 10)) {
	den *= 10;
    }
    pdinfo->structure = STACKED_TIME_SERIES;
    pdinfo->pd = nperiods;
    pdinfo->sd0 = 1.0 + 1.0 / den;
    ntodate_full(pdinfo->stobs, 0, pdinfo); 
    ntodate_full(pdinfo->endobs, pdinfo->n - 1, pdinfo);
}

int set_panel_structure_from_vars (int uv, int tv, 
				   double **Z, 
				   DATAINFO *pdinfo)
{
    double *uid = NULL;
    double *tid = NULL;
    char *mask = NULL;
    int n = pdinfo->n;
    int totmiss = 0;
    int fulln = 0;
    int nunits = 0;
    int nperiods = 0;
    int ustrs = 0;
    int err = 0;

    /* FIXME sub-sampled dataset (needs to be disallowed?) */

#if PDEBUG
    fprintf(stderr, "set_panel_structure_from_vars:\n "
	    "using var %d ('%s') for unit, var %d ('%s') for time\n",
	    uv, pdinfo->varname[uv], tv, pdinfo->varname[tv]);
    print_unit_var(uv, Z, n, 0);
#endif

    uid = copyvec(Z[uv], n);
    tid = copyvec(Z[tv], n);

    if (uid == NULL || tid == NULL) {
	err = E_ALLOC;
	goto bailout;
    }  

    qsort(uid, n, sizeof *uid, gretl_compare_doubles);
    nunits = count_distinct_values(uid, pdinfo->n);

    qsort(tid, n, sizeof *tid, gretl_compare_doubles);
    nperiods = count_distinct_values(tid, pdinfo->n);

    if (nunits == 1 || nperiods == 1 || 
	nunits == n || nperiods == n ||
	n > nunits * nperiods) {
	fprintf(stderr, "Dataset does not have a panel structure\n");
	err = E_DATA;
	goto bailout;
    }

    fulln = nunits * nperiods;
    totmiss = fulln - pdinfo->n;

#if PDEBUG
    fprintf(stderr, "Found %d units, %d periods\n", nunits, nperiods);
    fprintf(stderr, "Units: min %g, max %g\n", uid[0], uid[n - 1]);
    fprintf(stderr, "Periods: min %g, max %g\n", tid[0], tid[n - 1]);
    fprintf(stderr, "Required rows = %d * %d = %d\n", nunits, nperiods, fulln);
    fprintf(stderr, "Missing rows = %d - %d = %d\n", fulln, pdinfo->n, totmiss);
#endif

    /* sort full dataset by unit and period */
    err = panel_data_sort_by(Z, pdinfo, uv, tv, &ustrs);

#if PDEBUG
    print_unit_var(uv, Z, n, 1);
#endif

    if (!err && totmiss > 0) {
	/* do we want this? */
	mask = calloc(fulln, 1);
	rearrange_id_array(uid, nunits, n);
	rearrange_id_array(tid, nperiods, n);
	err = pad_panel_dataset(uid, uv, nunits, 
				tid, tv, nperiods, 
				Z, pdinfo, ustrs, 
				mask);
    }

    if (!err) {
	err = compose_panel_indices(pdinfo, nperiods, mask);
    }

    if (!err) {
	finalize_panel_datainfo(pdinfo, nperiods);
    }

 bailout:

    free(uid);
    free(tid);

    return err;
}

int set_panel_structure_from_line (const char *line, 
				   double **Z, 
				   DATAINFO *pdinfo)
{
    int n = pdinfo->n;
    int uv = 0, tv = 0;
    int err = 0;

    if (!strncmp(line, "setobs", 6)) {
	line += 7;
    }

    err = uv_tv_from_line(line, pdinfo, &uv, &tv);

    if (!err) {
	err = check_index_values(Z[uv], n);
    }

    if (!err) {
	err = check_index_values(Z[tv], n);
    }

    if (!err) {
	err = set_panel_structure_from_vars(uv, tv, Z, pdinfo);
    }

    return err;
}

/* utility functions */

int guess_panel_structure (double **Z, DATAINFO *pdinfo)
{
    int ret, v = varindex(pdinfo, "year");

    if (v == pdinfo->v) {
	v = varindex(pdinfo, "Year");
    }

    if (v == pdinfo->v) {
	ret = 0; /* can't guess */
    } else if (floateq(Z[v][0], Z[v][1])) { 
	/* "year" is same for first two obs */
	pdinfo->structure = STACKED_CROSS_SECTION; 
	ret = STACKED_CROSS_SECTION;
    } else {
	pdinfo->structure = STACKED_TIME_SERIES; 
	ret = STACKED_TIME_SERIES;
    }

    return ret;
}

int balanced_panel (const DATAINFO *pdinfo)
{
    int n = pdinfo->t2 - pdinfo->t1 + 1;
    int ret = 0;

    if (n % pdinfo->pd == 0) {
	char unit[OBSLEN], period[OBSLEN];

	if (sscanf(pdinfo->endobs, "%[^:]:%s", unit, period) == 2) {
	    if (atoi(period) == pdinfo->pd) {
		ret = 1;
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
	    strcpy(gretl_errmsg, "Panel models must include an intercept");
	    *err = E_DATA;
	    return NULL;
	}
    }

    if (gretl_model_get_int(orig, "fixed-effects")) {
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

    if (gretl_model_get_int(orig, "fixed-effects")) {
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

#if 0 /* not ready */

static int panel_obs_freq (FreqDist *fr, int *x, int n)
{
    int *ifreq = NULL, *ivals = NULL;
    int *sorted = NULL;
    int i, t, last;
    int err = 0;

    sorted = malloc(n * sizeof *sorted);
    if (sorted == NULL) {
	return E_ALLOC;
    }
	
    for (t=0; t<n; t++) {
	sorted[t] = x[t];
    }
	
    qsort(sorted, n, sizeof *sorted, gretl_compare_ints); 
    nbins = count_distinct_int_values(sorted, n);

    ifreq = malloc(nbins * sizeof *ifreq);
    ivals = malloc(nbins * sizeof *ivals);
    if (ifreq == NULL || ivals == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    ivals[0] = last = sorted[0];
    ifreq[0] = i = 1;

    for (t=1; t<n; t++) {
	if (sorted[t] != last) {
	    last = sorted[t];
	    ifreq[i] = 1;
	    ivals[i++] = last;
	} else {
	    ifreq[i-1] += 1;
	}
    }

    if (freq_add_arrays(freq, nbins)) {
	*err = E_ALLOC;
    } else {
	for (k=0; k<nbins; k++) {
	    freq->endpt[k] = freq->midpt[k] = ivals[k];
	    freq->f[k] = ifreq[k];
	}
	freq->endpt[nbins] = xmax;
    }

 bailout:

    free(sorted);
    free(ivals);
    free(ifreq);

    return err;
}

#endif

/* calculate the within and between standard deviations for a given
   variable in a panel data set */

int panel_variance_info (const double *x, const DATAINFO *pdinfo,
			 double xbar, double *psw, double *psb)
{
    double sw = 0.0, sb = 0.0;
    double xibar, d;
    int effn, effnT;
    int n, T, nT, Ti;
    int i, t, s;
    
    if (pdinfo->paninfo == NULL) {
	return E_PDWRONG;
    }

    nT = pdinfo->t2 - pdinfo->t1 + 1;
    T = pdinfo->pd;
    n = nT / T;

    effn = 0;
    effnT = 0;

    for (i=0; i<n; i++) {
	Ti = 0;
	xibar = 0.0;
	for (t=0; t<T; t++) {
	    s = pdinfo->t1 + i * T + t;
	    if (!na(x[s])) {
		Ti++;
		xibar += x[s];
	    }
	}
	if (Ti > 1) {
	    xibar /= Ti;
	    for (t=0; t<T; t++) {
		s = pdinfo->t1 + i * T + t;
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
	sb /= (effn - 1);
	sb = sqrt(sb);
    } else {
	sb = NADBL;
    }

    if (effnT - effn > 0) {
	sw /= (effnT - effn);
	sw = sqrt(sw);
    } else {
	sw = NADBL;
    }

    *psw = sw;
    *psb = sb;

    return 0;
}   

int panel_obs_info (const int *list, const double **Z, const DATAINFO *pdinfo,
		    PRN *prn)
{
    int *uobs = NULL;
    const int *unit;
    int minTi, maxTi;
    int jmax, vj;
    int n, Ti, ok;
    int i, j, t;
    
    if (pdinfo->paninfo == NULL) {
	return E_PDWRONG;
    }

    n = (pdinfo->t2 - pdinfo->t1 + 1) / pdinfo->pd;

    uobs = malloc(n * sizeof *uobs);
    if (uobs == NULL) {
	return E_ALLOC;
    }

    unit = pdinfo->paninfo->unit;

    jmax = (list != NULL)? list[0] : pdinfo->v - 1;

    maxTi = 0;
    minTi = pdinfo->pd;
    Ti = 0;
    i = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (t > pdinfo->t1 && unit[t] != unit[t-1]) {
	    if (Ti < minTi) {
		minTi = Ti;
	    } else if (Ti > maxTi) {
		maxTi = Ti;
	    }
	    uobs[i] = Ti;
	    Ti = 0;
	    i++;
	}
	ok = 1;
	for (j=1; j<=jmax; j++) {
	    vj = (list != NULL)? list[j] : j;
	    if (na(Z[vj][t])) {
		ok = 0;
		break;
	    }
	}
	Ti += ok;
	if (t == pdinfo->t2) {
	    if (Ti < minTi) {
		minTi = Ti;
	    } else if (Ti > maxTi) {
		maxTi = Ti;
	    }
	    uobs[i] = Ti;
	}
    }

    pprintf(prn, "Panel observations info\n");

    if (minTi == maxTi) {
	pprintf(prn, "%d units, each with %d observations\n", n, maxTi);
    } else {
	for (i=0; i<n; i++) {
	    pprintf(prn, "unit %d: %d observations\n", i+1, uobs[i]);
	}
	/* do a frequency distribution */
    }

    free(uobs);

    return 0;
}   
