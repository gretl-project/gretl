/*
 *   Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* panel estimation and dataset handling for gretl */

#include "libgretl.h"
#include "f2c.h"
#include "clapack_double.h"
#include "gretl_model.h"
#include "gretl_panel.h"

#define PDEBUG 0

/* The minimum number of observations we'll accept for a given
   cross-sectional unit, to include that unit in the fixed-effects
   regression.  This was previously set at 2.
*/

#define FE_MINOBS 1

enum vcv_ops {
    VCV_INIT,
    VCV_SUBTRACT
};

typedef struct panelmod_t_ panelmod_t;

struct panelmod_t_ {
    gretlopt opt;         /* option flags */
    int nunits;           /* total cross-sectional units */
    int effn;             /* effective (included) cross-section units */
    int T;                /* times-series length of panel */
    int effT;             /* effective times-series length (max usable obs per unit) */
    int ndum;             /* number of dummy variables added for FE model */
    int *unit_obs;        /* array of number of observations per x-sect unit */
    char *varying;        /* array to record properties of pooled-model regressors */
    int *vlist;           /* list of time-varying variables from pooled model */
    int balanced;         /* 1 if the model dataset is balanced, else 0 */
    int nbeta;            /* number of slope coeffs for Hausman test */
    int Fdfn;             /* numerator df, F for differing intercepts */
    int Fdfd;             /* denominator df, F for differing intercepts */
    double ybar;          /* mean of dependent variable */
    double sdy;           /* standard deviation of dependent variable */
    double tss;           /* total sum of squares, dependent variable */
    double sigma_e;       /* fixed-effects standard error */
    double theta;         /* quasi-demeaning coefficient */
    double F;             /* joint significance of differing unit intercepts */
    double BP;            /* Breusch-Pagan test statistic */
    double H;             /* Hausman test statistic */
    double *bdiff;        /* array of coefficient differences */
    double *sigma;        /* Hausman covariance matrix */
    double within_s2;     /* fixed-effects or "within" error variance */
    double between_s2;    /* group-means or "between" error variance */
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
   overall index into the data set */
#define panel_index(i,t) (i * panidx.T + t + panidx.offset)

static void 
panel_index_init (const DATAINFO *pdinfo, int nunits, int T)
{
    panidx.n = nunits;
    panidx.T = T;
    panidx.offset = pdinfo->t1;
}

static void 
panelmod_init (panelmod_t *pan, MODEL *pmod, gretlopt opt)
{
    pan->nunits = 0;
    pan->effn = 0;
    pan->T = 0;
    pan->effT = 0;
    pan->ndum = 0;
    pan->unit_obs = NULL;
    pan->varying = NULL;
    pan->vlist = NULL;
    pan->opt = opt;

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
    
    pan->pooled = pmod;
    pan->realmod = NULL;
}

static void panelmod_free (panelmod_t *pan)
{
    free(pan->unit_obs);
    free(pan->varying);
    free(pan->vlist);

    free(pan->bdiff);
    free(pan->sigma);

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

/* retrieve X'X^{-1} from the fixed effects model */

static gretl_matrix *fe_model_xpx (MODEL *pmod)
{
    gretl_matrix *X;
    int k = pmod->ncoeff;
    double x;
    int i, j, m;

    makevcv(pmod, 1.0); /* invert X'X into pmod->vcv */

    X = gretl_matrix_alloc(k, k);
    if (X == NULL) {
	return NULL;
    }  

    m = 0;
    for (j=0; j<k; j++) {
	for (i=j; i<k; i++) {
	    x = pmod->vcv[m++];
	    X->val[mdx(X, i, j)] = X->val[mdx(X, j, i)] = x;
	}
    }

    return X;
}

/* HAC covariance matrix for the fixed-effects model */

static int 
fe_robust_vcv (MODEL *pmod, panelmod_t *pan, const double **Z)
{
    gretl_vector *e = NULL;
    gretl_matrix *Xi = NULL;
    gretl_vector *eXi = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *XX = NULL;
    gretl_matrix *W = NULL;
    int Ti, T = pan->effT;
    int k = pmod->ncoeff;
    int i, j, v, s, t;
    int err = 0;

    e = gretl_vector_alloc(T);
    Xi = gretl_matrix_alloc(T, k);
    eXi = gretl_vector_alloc(k);
    tmp = gretl_matrix_alloc(k, k);
    XX = fe_model_xpx(pmod);
    W = gretl_matrix_alloc(k, k);

    if (e == NULL || Xi == NULL || eXi == NULL ||
	tmp == NULL || XX == NULL || W == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_zero(W);

    s = 0;
    for (i=0; i<pan->nunits; i++) {
	Ti = pan->unit_obs[i];
	if (Ti < FE_MINOBS) {
	    continue;
	}
	e = gretl_matrix_reuse(e, Ti, 1);
	Xi = gretl_matrix_reuse(Xi, Ti, k);
	for (t=0; t<Ti; t++) {
	    gretl_vector_set(e, t, pmod->uhat[s]);
	    for (j=0; j<k; j++) {
		v = pmod->list[j + 2];
		gretl_matrix_set(Xi, t, j, Z[v][s]);
	    }
	    s++;
	}
	gretl_matrix_multiply_mod(e, GRETL_MOD_TRANSPOSE,
				  Xi, GRETL_MOD_NONE,
				  eXi);
	gretl_matrix_multiply_mod(eXi, GRETL_MOD_TRANSPOSE,
				  eXi, GRETL_MOD_NONE,
				  tmp);
	gretl_matrix_add_to(W, tmp);
    }

#if 0
    gretl_matrix_print(W, "sum_eXi");
#endif

    /* now form V(b_W) = (X'X)^{-1} W (X'X)^{-1} */

    gretl_matrix_multiply(XX, W, tmp);
    gretl_matrix_multiply(tmp, XX, W);

#if 0
    gretl_matrix_print(W, "V(b_W)");
#endif

    s = 0;
    for (i=0; i<k; i++) {
	for (j=i; j<k; j++) {
	    pmod->vcv[s++] = W->val[mdx(W, i, j)];
	}
	pmod->sderr[i] = sqrt(W->val[mdx(W, i, i)]);
    }

    gretl_model_set_int(pmod, "panel_hac", 1);

 bailout:

    gretl_matrix_free(e);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);
    gretl_matrix_free(tmp);
    gretl_matrix_free(XX);
    gretl_matrix_free(W);

    return err;
}

/* Durbin-Watson statistic for pooled or fixed effects model on a 
   balanced panel */

void panel_dwstat (MODEL *pmod, const DATAINFO *pdinfo)
{
    int T = pdinfo->pd;
    double ut, u1;
    double num = 0.0;
    int i, t, s, n;

    pmod->dw = NADBL;
    pmod->rho = NADBL;

    if (pmod->ess <= 0.0) {
	return;
    }

    if (pmod->nobs % T != 0) {
#if PDEBUG
	fprintf(stderr, "panel_dwstat: skipping: pmod->nobs = %d, T = %d\n", 
		pmod->nobs, T);
#endif
	return;
    }

    n = pmod->nobs / T;

    s = pmod->t1;
    for (i=0; i<n; i++) {
#if PDEBUG
	fprintf(stderr, "panel_dwstat: skipping obs %d\n", s);
#endif
	s++;
	for (t=1; t<T; t++) {
	    ut = pmod->uhat[s];
	    u1 = pmod->uhat[s-1];
	    if (!na(ut) && !na(u1)) {
		num += (ut - u1) * (ut - u1);
	    }
	    s++;
	}
    }

    pmod->dw = num / pmod->ess;
}

static void private_panel_dwstat (MODEL *pmod, DATAINFO *pdinfo, 
				  const panelmod_t *pan)
{
    if (pan->ndum > 0) {
	panel_dwstat(pmod, pdinfo);
    } else if (pan->balanced) {
	int pd = pdinfo->pd;

	pdinfo->pd = pan->effT;
	panel_dwstat(pmod, pdinfo);
	pdinfo->pd = pd;
    } else {
	pmod->rho = pmod->dw = NADBL;
    }
}

/* Allocate the arrays needed to perform the Hausman test,
   in its matrix formulation.
*/

static int hausman_allocate (panelmod_t *pan)
{
    int k = pan->vlist[0] - 2;
    int nsigma;

    pan->nbeta = k;

    if (pan->opt & OPT_H) {
	/* regression approach to Hausman: we don't need
	   the allocations below */
	return 0;
    }

    nsigma = k * (k + 1) / 2;

    /* array to hold differences between coefficient estimates */
    pan->bdiff = malloc(k * sizeof *pan->bdiff);
    if (pan->bdiff == NULL) {
	return E_ALLOC;
    }

    /* array to hold covariance matrix */
    pan->sigma = malloc(nsigma * sizeof *pan->sigma);
    if (pan->sigma == NULL) {
	free(pan->bdiff);
	pan->bdiff = NULL;
	return E_ALLOC; 
    }

    return 0;
}   

/* printing routines for the old "panel diagnostics" test */

static void print_panel_coeff (const MODEL *pmod, 
			       const char *vname,
			       int i, PRN *prn)
{
    double tstat = pmod->coeff[i] / pmod->sderr[i];
    char errstr[18];
    char pvstr[18];

    sprintf(errstr, "(%.5g)", pmod->sderr[i]);
    sprintf(pvstr, "[%.5f]", t_pvalue_2(tstat, pmod->dfd));
    pprintf(prn, "%*s: %14.5g %15s %15s\n", VNAMELEN, vname,
	    pmod->coeff[i], errstr, pvstr);
}

static void print_re_model_top (const panelmod_t *pan, PRN *prn)
{
    pputs(prn, "Variance estimators:\n");
    pprintf(prn, " between = %g\n", pan->between_s2);
    pprintf(prn, " within = %g\n", pan->within_s2);

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

/* Fix for statistics on the dependent variable, for the case
   where we do the "within" regression on an auxiliary de-meaned
   data set.
*/

static void within_depvarstats (panelmod_t *pan, double *y, int n)
{
    double x = 0.0;
    int t;

    for (t=0; t<n; t++) {
	x += y[t];
    }

    pan->ybar = x / n;
    pan->tss = 0.0;

    for (t=0; t<n; t++) {
	x = y[t] - pan->ybar;
	pan->tss += x * x;
    }

    pan->sdy = sqrt(pan->tss / (n - 1));
}

static DATAINFO *dataset_with_panel_info (double ***pZ, int k, int n)
{
    DATAINFO *dinfo;

    dinfo = create_new_dataset(pZ, k, n, 0);
    if (dinfo == NULL) {
	return NULL;
    }

    if (dataset_allocate_panel_info(dinfo)) {
	destroy_dataset(*pZ, dinfo);
	return NULL;
    }

    return dinfo;
}

/* Construct a version of the dataset from which the group means
   are subtracted, for the "within" regression when not using
   the dummy variables approach.  Nota bene: this auxiliary
   dataset is not necessarily of full length.
*/

static DATAINFO *
within_groups_dataset (const double **Z, double ***wZ, panelmod_t *pan)
{
    DATAINFO *winfo;
    double *wy = NULL;
    int wnobs = 0;
    double xbar;
    int i, j, k;
    int s, t, bigt;

    /* real number of included units */
    pan->effn = 0;
    pan->balanced = 1;

    for (i=0; i<pan->nunits; i++) { 
	if (pan->unit_obs[i] >= FE_MINOBS) {
	    wnobs += pan->unit_obs[i];
	    pan->effn += 1;
	} else if (pan->unit_obs[i] != pan->effT) {
	    pan->balanced = 0;
	}
    }

    wy = malloc(wnobs * sizeof *wy);
    if (wy == NULL) {
	return NULL;
    }

#if PDEBUG
    fprintf(stderr, "within_groups dataset: nvars=%d, nobs=%d\n", 
	    pan->vlist[0], wnobs);
#endif

    winfo = dataset_with_panel_info(wZ, pan->vlist[0], wnobs);
    if (winfo == NULL) {
	free(wy);
	return NULL;
    }

    k = 0;
    for (j=1; j<=pan->vlist[0]; j++) { 
	if (pan->vlist[j] == 0) {
	    continue;
	} 

	s = 0;
	k++;

#if PDEBUG
	fprintf(stderr, "working on list[%d] (original var %d, var %d in wZ)\n",
		j, pan->vlist[j], k);
#endif
	for (i=0; i<pan->nunits; i++) { 
	    int Ti = pan->unit_obs[i];
	    int got = 0;

#if PDEBUG
	    fprintf(stderr, "looking at x-sect unit %d: Ti = %d\n", i, Ti);
#endif
	    if (Ti < FE_MINOBS) {
		continue;
	    }

	    xbar = 0.0;
	    for (t=0; t<pan->T; t++) {
		bigt = panel_index(i, t);
		if (!na(pan->pooled->uhat[bigt])) {
		    xbar += Z[pan->vlist[j]][bigt];
		}
	    }
	    xbar /= (double) Ti;

#if PDEBUG > 1
	    fprintf(stderr, "xbar for var %d, unit %d = %g\n", 
		    pan->vlist[j], i, xbar);
#endif
	    for (t=0; t<pan->T && got < Ti; t++) {
		if (s >= winfo->n) {
		    fprintf(stderr, "*** Error: overflow of wZ at unit %d:\n" 
			    "  pan->T = %d, winfo->n = %d, hit s = %d at t = %d\n",  
			    i, pan->T, winfo->n, s, t);
		    break;
		}
		bigt = panel_index(i, t);
		if (!na(pan->pooled->uhat[bigt])) {
		    (*wZ)[k][s] = Z[pan->vlist[j]][bigt] - xbar;
#if PDEBUG > 1
		    fprintf(stderr, "Set wZ[%d][%d] = %g\n", k, s, (*wZ)[k][s]);
#endif
		    if (j == 1) {
			wy[s] = Z[pan->vlist[j]][bigt];
			winfo->paninfo->unit[s] = i;
			winfo->paninfo->period[s] = t;
		    }
		    got++;
		    s++;
		}
	    }
	}
    }

    within_depvarstats(pan, wy, wnobs);
    free(wy);

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
			const double **gZ, double ***reZ, 
			int *relist, int *hlist, panelmod_t *pan)
{
    DATAINFO *reinfo;
    int hreg = (hlist != NULL);
    int v1 = relist[0];
    int v2 = 0;
    int re_n = 0;
    int i, j, k, k2, t;

    for (i=0; i<pan->nunits; i++) { 
	if (pan->unit_obs[i] > 0) {
	    re_n += pan->unit_obs[i];
	}
    }

    if (hreg) {
	/* apparatus for regression version of Hausman test */
	for (i=1; i<pan->vlist[0]; i++) {
	    if (pan->vlist[i+1] != 0) {
		hlist[0] = ++v2;
		hlist[i-1] = v1 + i - 2;
	    }
	}
    }  

#if PDEBUG
    fprintf(stderr, "random_effects_dataset: nvars=%d, nobs=%d\n",
	    v1 + v2, re_n);
#endif

    reinfo = dataset_with_panel_info(reZ, v1 + v2, re_n);
    if (reinfo == NULL) {
	return NULL;
    }

    pan->balanced = 1;

    /* Now create the transformed variables: original data minus theta
       times the appropriate group mean.
    */    
    k = 0;
    k2 = v1 - 1;
    for (j=1; j<=v1; j++) {
	int pj = pan->pooled->list[j];
	const double *xj = Z[pj];
	const double *gm = NULL;
	int vvar = var_is_varying(pan, pj);
	double theta_i;
	int bigt, s, u = 0;

	if (pj == 0) {
	    relist[j] = 0;
	} else {
	    k++;
	    relist[j] = k; /* build GLS regression list */
	    gm = gZ[k];
	    strcpy(reinfo->varname[k], pdinfo->varname[pj]);
	    if (hreg && vvar && j > 1) {
		k2++;
		strcpy(reinfo->varname[k2], "_");
		strncat(reinfo->varname[k2], pdinfo->varname[pj],
			VNAMELEN - 2);
	    }
	}

	s = 0;
	for (i=0; i<pan->nunits; i++) {
	    int Ti = pan->unit_obs[i];

	    if (Ti == 0) {
		continue;
	    }

	    if (Ti != pan->effT) {
		theta_i = 1.0 - sqrt(pan->within_s2 / (Ti * pan->between_s2));
#if PDEBUG
		fprintf(stderr, "unit %d: T_i = %d, theta_i = %g\n",
			i, Ti, theta_i);
#endif
		pan->balanced = 0;
	    } else {
		theta_i = pan->theta;
	    }

	    for (t=0; t<pan->T; t++) {
		bigt = panel_index(i, t);
		if (!na(pan->pooled->uhat[bigt])) {
		    if (relist[j] == 0) {
			/* the intercept */
			(*reZ)[0][s] -= theta_i;
		    } else {
			(*reZ)[k][s] = xj[bigt] - theta_i * gm[u];
#if PDEBUG > 1
			fprintf(stderr, "Set reZ[%d][%d] = %g\n", k, s, (*reZ)[k][s]);
#endif
			if (hreg && vvar) {
			    (*reZ)[k2][s] = xj[bigt] - gm[u];
			}
			if (k == 1) {
			    reinfo->paninfo->unit[s] = i;
			    reinfo->paninfo->period[s] = t;
			}
		    }
		    s++;
		}
	    }
	    u++;
	}
    }

    return reinfo;
}

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
    int gn = 0;
    int i, j, k;
    int s, t, bigt;

    for (i=0; i<pan->nunits; i++) { 
	if (pan->unit_obs[i] > 0) {
	    gn++;
	}
    }

#if PDEBUG
    fprintf(stderr, "group_means_dataset: nvars=%d, nobs=%d\n", 
	    pan->pooled->list[0], gn);
#endif

    ginfo = create_new_dataset(gZ, pan->pooled->list[0], gn, 0);
    if (ginfo == NULL) {
	return NULL;
    }

    k = 1;
    for (j=1; j<=pan->pooled->list[0]; j++) { 
	int pj = pan->pooled->list[j];

	if (pj == 0) {
	    continue;
	}
	s = 0;
	for (i=0; i<pan->nunits; i++) { 
	    int Ti = pan->unit_obs[i];

	    if (Ti > 0) {
		x = 0.0;
		for (t=0; t<pan->T; t++) {
		    bigt = panel_index(i, t);
		    if (!na(pan->pooled->uhat[bigt])) {
			x += Z[pj][bigt];
		    }
		}
		(*gZ)[k][s] = x / (double) Ti;
#if PDEBUG > 1
		fprintf(stderr, "Set gZ[%d][%d] = %g\n", k, i, (*gZ)[k][s]);
#endif
		s++;
	    }
	}
	k++;
    }

    return ginfo;
}

/* calculate the group means or "between" regression and its error
   variance */

static int
between_variance (panelmod_t *pan, double ***gZ, DATAINFO *ginfo)
{
    MODEL bmod;
    int *blist;
    int i, j;
    int err = 0;

    blist = gretl_list_new(pan->pooled->list[0]);
    if (blist == NULL) {
	return E_ALLOC;
    }

    j = 1;
    for (i=1; i<=blist[0]; i++) { 
	if (pan->pooled->list[i] != 0) {
	    blist[i] = j++;
	}
    }

    bmod = lsq(blist, gZ, ginfo, OLS, OPT_A | OPT_Z);

    if (bmod.errcode == 0) {
	pan->between_s2 = bmod.sigma * bmod.sigma;
    } else {
	err = bmod.errcode;
#if PDEBUG
	fprintf(stderr, "error %d in between_variance\n", err);
#endif
    }

    clear_model(&bmod);
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
#if PDEBUG
	    fprintf(stderr, "setting sigma[%d] using vcv[%d] (%d,%d) = %g\n",
		    k, idx, mi, mj, pmod->vcv[idx]);
#endif
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

static void panel_df_correction (MODEL *pmod, int ndf)
{
    double dfcorr = sqrt((double) pmod->dfd / (pmod->dfd - ndf));
    int i;

    pmod->dfd -= ndf;
    pmod->dfn += ndf - 1; 

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
			     MODEL *wmod,
			     DATAINFO *pdinfo,
			     PRN *prn)
{
    int nslopes = wmod->ncoeff - pan->ndum;
    int dfn, i, j;

#if PDEBUG
    fprintf(stderr, "f.e. variance, wmod.ncoeff=%d, var=%g, "
	    "list[0]=%d\n", wmod->ncoeff, pan->within_s2, pan->vlist[0]);
#endif

    pputs(prn, 
	  _("Fixed effects estimator\n"
	    "allows for differing intercepts by cross-sectional unit\n"
	    "slope standard errors in parentheses, p-values in brackets\n"));

    if (pan->ndum > 0) {
	pputs(prn, _("a_i = intercepts"));
	pputs(prn, "\n\n");
    } else {
	pputc(prn, '\n');
    }

    /* print the slope coefficients, for varying regressors */
    for (i=0; i<nslopes; i++) {
	int vi = pan->vlist[i+3];

	print_panel_coeff(wmod, pdinfo->varname[vi], i, prn);
    }
    pputc(prn, '\n');   

    /* if we used LSDV approach, print the per-unit intercept 
       estimates */
    if (pan->ndum > 0) {
	for (i=0, j=0; i<pan->nunits; i++) {
	    if (pan->unit_obs[i] >= FE_MINOBS) {
		char dumstr[VNAMELEN];
		double b;

		b = wmod->coeff[j + nslopes];
		sprintf(dumstr, "a_%d", i + 1);
		pprintf(prn, "%*s: %14.4g\n", VNAMELEN, dumstr, b);
		j++;
	    }
	}
    } else {
	pprintf(prn, _("%d group means were subtracted from the data"), pan->effn);
	pputc(prn, '\n');
    }

    dfn = pan->effn - 1;

    pprintf(prn, _("\nResidual variance: %g/(%d - %d) = %g\n"), 
	    wmod->ess, wmod->nobs, pan->vlist[0] - 1 + dfn, pan->within_s2);

    if (pan->ndum > 0) {
	pputs(prn, _("Joint significance of unit dummy variables:\n"));
    } else {
	pprintf(prn, _("Joint significance of differing group means:\n"));
    }

    pprintf(prn, " F(%d, %d) = %g %s %g\n", pan->Fdfn, pan->Fdfd, pan->F, 
	    _("with p-value"), f_cdf_comp(pan->F, pan->Fdfn, pan->Fdfd));

    pputs(prn, _("(A low p-value counts against the null hypothesis that "
		 "the pooled OLS model\nis adequate, in favor of the fixed "
		 "effects alternative.)\n\n"));

    return 0;
}

static void save_fixed_effects_F (panelmod_t *pan, MODEL *wmod)
{
    ModelTest *test = model_test_new(GRETL_TEST_PANEL_F);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_F);
	model_test_set_dfn(test, pan->Fdfn);
	model_test_set_dfd(test, pan->Fdfd);
	model_test_set_value(test, pan->F);
	model_test_set_pvalue(test, f_cdf_comp(pan->F, pan->Fdfn, pan->Fdfd));
	maybe_add_test_to_model(wmod, test);
    }	    
}

static void fixed_effects_F (panelmod_t *pan, MODEL *wmod)
{
    pan->Fdfn = pan->effn - 1;
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

    /* remove any non-varying variables */
    for (i=2; i<=targ->list[0]; i++) {
	if (!in_gretl_list(pan->vlist, targ->list[i])) {
	    gretl_list_delete_at_pos(targ->list, i--);
	}
    }

    if (pan->opt & OPT_F) {
	/* fixed effects: remove the const */
	gretl_list_delete_at_pos(targ->list, 2);
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

static int fix_within_stats (MODEL *targ, panelmod_t *pan)
{
    int nc = targ->ncoeff;
    int err = 0;

    err = fix_panelmod_list(targ, pan);
    if (err) {
	return err;
    }

    targ->ybar = pan->ybar;
    targ->sdy = pan->sdy;
    targ->ifc = 1;

    targ->rsq = 1.0 - (targ->ess / pan->tss);

    if (targ->dfd > 0) {
	double den = pan->tss * targ->dfd;

	targ->adjrsq = 1 - (targ->ess * (targ->nobs - 1) / den);
    }

    if (targ->rsq < 0.0) {
	targ->rsq = 0.0;
    } else {
	targ->fstt = (targ->rsq / (1.0 - targ->rsq)) * 
	    ((double) targ->dfd / targ->dfn);
    }

    targ->ncoeff = targ->dfn + 1;
    ls_criteria(targ);
    targ->ncoeff = nc;

    return err;
}

/* Fixed-effects model: add the per-unit intercept estimates
   ("ahat") to the model in case the user wants to retrieve
   them.
*/

static int fe_model_add_ahat (MODEL *pmod, const double **Z, 
			      const DATAINFO *pdinfo,
			      panelmod_t *pan)
{
    double *ahat = NULL;
    size_t sz = pdinfo->n * sizeof *ahat;
    int i, j, s, t;
    int err = 0;

    ahat = malloc(sz);
    if (ahat == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	ahat[t] = NADBL;
    }

    if (pan->ndum > 0) {
	/* used LSDV approach */
	j = pmod->ncoeff - pan->ndum;
	for (i=0; i<pan->nunits; i++) {
	    if (pan->unit_obs[i] > 0) {
		for (t=0; t<pan->T; t++) {
		    s = panel_index(i, t);
		    if (!na(pan->pooled->uhat[s])) {
			ahat[s] = pmod->coeff[j];
		    }
		}
		j++;
	    }
	}
    } else {
	/* used de-meaned data: FIXME indexing? */
	for (i=0; i<pan->nunits; i++) {
	    int Ti = pan->unit_obs[i];

	    if (Ti > 1) {
		double ahi = 0.0;

		for (t=0; t<pan->T; t++) {
		    s = panel_index(i, t);
		    if (!na(pmod->uhat[s])) {
			ahi += Z[pmod->list[1]][s];
			for (j=0; j<pmod->ncoeff; j++) {
			    ahi -= pmod->coeff[j] * Z[pmod->list[j+2]][s];
			}
		    }
		}
		ahi /= Ti;
		for (t=0; t<pan->T; t++) {
		    s = panel_index(i, t);
		    if (!na(pmod->uhat[s])) {
			ahat[s] = ahi;
		    }
		}
	    }
	}
    }

    err = gretl_model_set_data(pmod, "ahat", ahat, 
			       MODEL_DATA_DOUBLE_ARRAY, 
			       sz);

    return err;
}

/* Fix uhat and yhat in two cases: (a) when using a de-meaned dataset,
   we need to ensure that the uhat, yhat values get written to the
   right observation slots in the full dataset; (b) when estimating
   the random effects model we need to compute residuals based on
   the untransformed data (and again, place them correctly in
   relation to the full dataset).
*/

static void 
fix_panel_hatvars (MODEL *pmod, const DATAINFO *dinfo, 
		   panelmod_t *pan, const double **Z)
{
    const double *y = NULL;
    double *uhat = pan->pooled->uhat;
    double *yhat = pan->pooled->yhat;
    int n = pan->pooled->full_n;
    double yht;
    int i, j, t, u, p;

    if (Z != NULL) {
	y = Z[pan->pooled->list[0]];
	pmod->ess = 0.0;
    }

    pan->pooled->uhat = NULL;
    pan->pooled->yhat = NULL;

    for (i=0; i<n; i++) {
	uhat[i] = yhat[i] = NADBL;
    }

    for (i=0; i<pmod->full_n; i++) {
	u = dinfo->paninfo->unit[i];
	p = dinfo->paninfo->period[i];
	t = panel_index(u, p);
	if (Z != NULL) {
	    yht = 0.0;
	    for (j=0; j<pmod->ncoeff; j++) {
		yht += pmod->coeff[j] * Z[pan->pooled->list[j+2]][t];
	    }
	    yhat[t] = yht;
	    uhat[t] = y[t] - yht;
	    pmod->ess += uhat[t] * uhat[t];
	} else {
	    uhat[t] = pmod->uhat[i];
	    yhat[t] = pmod->yhat[i];
	}
	if (i == 0) {
	    pmod->t1 = t;
	} else if (i == pmod->full_n - 1) {
	    pmod->t2 = t;
	}
    }

    pmod->full_n = n;

    free(pmod->uhat);
    pmod->uhat = uhat;
    
    free(pmod->yhat);
    pmod->yhat = yhat;
}

/* Estimate the fixed-effects model using the approach of creating a
   parallel dataset with the group means subtracted from all
   variables.
*/

static MODEL 
fixed_effects_by_demeaning (panelmod_t *pan, const double **Z, 
			    DATAINFO *pdinfo, PRN *prn)
{
    MODEL femod;
    gretlopt lsqopt = OPT_A | OPT_Z | OPT_C;
    double **wZ = NULL;
    DATAINFO *winfo = NULL;
    int *felist = NULL;
    int i;

#if PDEBUG
    fprintf(stderr, "fixed_effects: using de-meaned data\n");
#endif

    gretl_model_init(&femod);

    felist = gretl_list_new(pan->vlist[0] - 1);
    if (felist == NULL) {
	femod.errcode = E_ALLOC;
	return femod;
    }	

    winfo = within_groups_dataset(Z, &wZ, pan);
    if (winfo == NULL) {
	free(felist);
	femod.errcode = E_ALLOC;
	return femod;
    } 

    for (i=1; i<=felist[0]; i++) {
	felist[i] = i;
    }

    femod = lsq(felist, &wZ, winfo, OLS, lsqopt);

    if (femod.errcode) {
	; /* pass on */
    } else if (femod.list[0] < felist[0]) {
	/* one or more variables were dropped, because they were
	   all zero -- this is a symptom of collinearity */
	femod.errcode = E_SINGULAR;
    } else {
	/* we estimated a bunch of group means, and have to
	   subtract that many degrees of freedom */
	panel_df_correction(&femod, pan->effn);
#if PDEBUG
	printmodel(&femod, winfo, OPT_O, prn);
#endif
	if (pan->opt & OPT_F) {
	    /* saving the FE model */
	    fix_panel_hatvars(&femod, winfo, pan, NULL);
	    if (pan->opt & OPT_R) {
		fe_robust_vcv(&femod, pan, (const double **) wZ);
	    } 
	}
    }

    destroy_dataset(wZ, winfo);
    free(felist);

    return femod;
}

/* Estimate the fixed-effects model using the Least Squares Dummy
   Variables approach */

static MODEL 
fixed_effects_by_LSDV (panelmod_t *pan, double ***pZ, DATAINFO *pdinfo,
		       PRN *prn)
{
    MODEL femod;
    gretlopt lsqopt = OPT_A | OPT_Z | OPT_C;
    int *felist;
    int oldv = pdinfo->v;
    int i, j;

#if PDEBUG
    fprintf(stderr, "fixed_effects: using LSDV\n");
#endif

    gretl_model_init(&femod);

    felist = gretl_list_new(pan->vlist[0] + pan->ndum - 1);
    if (felist == NULL) {
	femod.errcode = E_ALLOC;
	return femod;
    }

    if (dataset_add_series(pan->ndum, pZ, pdinfo)) {
	free(felist);
	femod.errcode = E_ALLOC;
	return femod;
    }

    /* create the per-unit dummies for all included units */

    for (i=0, j=0; i<pan->nunits && j<pan->ndum; i++) {
	int t, dv = oldv + j;

	if (pan->unit_obs[i] < FE_MINOBS) {
	    continue;
	}

	sprintf(pdinfo->varname[dv], "unit_%d", i + 1);

	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[dv][t] = 0.0;
	}

	for (t=0; t<pan->T; t++) {
	    (*pZ)[dv][panel_index(i, t)] = 1.0;
	}

	j++;
    }

    /* regression list: copy pooled list, dropping the constant and
       adding the unit dummies
    */
    j = 1;
    for (i=1; i<=pan->vlist[0]; i++) {
	if (pan->vlist[i] != 0) {
	    felist[j++] = pan->vlist[i];
	}
    }
    for (i=0; i<pan->ndum; i++) {
	felist[j++] = oldv + i;
    }

#if PDEBUG
    printlist(felist, "felist");
#endif

    femod = lsq(felist, pZ, pdinfo, OLS, lsqopt);

#if PDEBUG
    if (!femod.errcode) {
	printmodel(&femod, pdinfo, OPT_NONE, prn);
    }
#endif

    dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo);
    free(felist);

    return femod;
}

/* Determine whether we should estimate the fixed effects model
   using the LSDV approach, or by de-meaning the data.  We actually
   return the number of dummy variables to be used, which is
   zero in the case of demeaning.
*/

static int get_fixed_effects_method (panelmod_t *pan)
{
    int i, nd = 0;

#if PDEBUG
    /* allow override via an environment variable */
    if (getenv("NODUMMIES") != NULL) {
	return 0;
    }
#endif

#if 1
    if ((pan->opt & OPT_R) && (pan->opt & OPT_F)) {
	fprintf(stderr, "Robust: will use demeaning\n");
	return 0;
    }
#endif

    for (i=0; i<pan->nunits; i++) {
#if PDEBUG
	fprintf(stderr, "pan unit %d: obs (based on OLS) = %d\n", 
		i, pan->unit_obs[i]);
#endif
	if (pan->unit_obs[i] >= FE_MINOBS) {
	    nd++;
	} else if (pan->unit_obs[i] > 0) {
	    /* a unit we need to exclude: hard to handle with
	       the dummies approach */
	    nd = 0;
	    break;
	}
    }

    if (nd <= 20 || nd < pan->vlist[0]) {
	pan->ndum = nd;
    }

    return pan->ndum;
}

/* Construct a gretl list containing the index numbers of the
   cross-sectional units included in the fixed-effects regression
   FIXME: we should actually read those index numbers from 
   datainfo.  This is for the purpose of naming the per-unit
   intercepts.
*/

int *fe_units_list (const panelmod_t *pan)
{
    int *ulist = NULL;
    int i, j, n = 0;

    for (i=0; i<pan->nunits; i++) { 
	if (pan->unit_obs[i] >= FE_MINOBS) {
	    n++;
	}
    }

    ulist = gretl_list_new(n);

    if (ulist != NULL) {
	j = 1;
	for (i=0; i<pan->nunits; i++) { 
	    if (pan->unit_obs[i] >= FE_MINOBS) {
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
    int j, pj, i = 1;

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
	    pj = pan->pooled->list[j];
	    if (!in_gretl_list(pan->vlist, pj)) {
		dlist[i++] = pj;
	    }
	}
    }

    return gretl_model_set_list_as_data(pmod, "droplist", dlist);
}

/* spruce up femod and attach it to pan */

static int save_fixed_effects_model (MODEL *femod, panelmod_t *pan,
				     const double **Z,
				     DATAINFO *pdinfo)
{
    int *ulist;
    int err = 0;

    femod->ci = PANEL;
    gretl_model_set_int(femod, "fixed-effects", 1);

    if (pan->ndum == 0) {
	err = fix_within_stats(femod, pan);
	if (err) {
	    return err;
	}
    } 

    /* compose list of dropped variables, if any */
    compose_panel_droplist(femod, pan);

    ulist = fe_units_list(pan);
    gretl_model_add_panel_varnames(femod, pdinfo, ulist);
    free(ulist);

    fe_model_add_ahat(femod, Z, pdinfo, pan);
    femod->list[0] -= pan->ndum;
    set_model_id(femod);
    private_panel_dwstat(femod, pdinfo, pan);
    save_fixed_effects_F(pan, femod);

    *pan->realmod = *femod;

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

    get_fixed_effects_method(pan);

    if (pan->ndum > 0) {
	femod = fixed_effects_by_LSDV(pan, pZ, pdinfo, prn);
    } else {
	femod = fixed_effects_by_demeaning(pan, (const double **) *pZ, 
					   pdinfo, prn);
    }

    if (femod.errcode) {
	pputs(prn, _("Error estimating fixed effects model\n"));
	errmsg(femod.errcode, prn);
	err = femod.errcode;
	clear_model(&femod);
    } else {
	pan->within_s2 = femod.sigma * femod.sigma;

	fixed_effects_F(pan, &femod);

	if (pan->opt & OPT_V) {
	    print_fe_results(pan, &femod, pdinfo, prn);
	}

	if (pan->bdiff != NULL && pan->sigma != NULL) {
	    int nslopes = femod.ncoeff - pan->ndum;

	    for (i=0; i<nslopes; i++) {
		pan->bdiff[i] = femod.coeff[i];
	    }
	    makevcv(&femod, femod.sigma);
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

/* spruce up remod and attach it to pan */

static void save_random_effects_model (MODEL *remod, panelmod_t *pan,
				       const double **Z,
				       const DATAINFO *reinfo)
{
    remod->ci = PANEL;

    gretl_model_set_int(remod, "random-effects", 1);
    gretl_model_set_double(remod, "within-variance", pan->within_s2);
    gretl_model_set_double(remod, "between-variance", pan->between_s2);

    if (pan->balanced) {
	gretl_model_set_double(remod, "gls-theta", pan->theta);
    }

    /* compose list of dropped variables, if any */
    compose_panel_droplist(remod, pan);

    gretl_model_add_panel_varnames(remod, reinfo, NULL);
    fix_panel_hatvars(remod, reinfo, pan, Z);
    fix_gls_stats(remod, pan);
    set_model_id(remod);

    *pan->realmod = *remod;
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
		chisq_cdf_comp(pan->H, pan->nbeta));
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
	model_test_set_pvalue(test, chisq_cdf_comp(pan->H, pan->nbeta));
	maybe_add_test_to_model(pan->realmod, test);
    }	    
}

/* Calculate the random effects regression.  Print the results
   here if we're doing the "panel diagnostics" test, otherwise
   save the results.
*/

static int random_effects (panelmod_t *pan, 
			   const double **Z, DATAINFO *pdinfo, 
			   const double **gZ, PRN *prn)
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
       will actually vary across the units.
    */
    pan->theta = 1.0 - sqrt(pan->within_s2 / (pan->effT * pan->between_s2));
#if PDEBUG
    fprintf(stderr, "random_effects: theta = %g\n", pan->theta);
#endif

    /* make special transformed dataset, and regression list */
    reinfo = random_effects_dataset(Z, pdinfo, gZ, &reZ, relist, hlist, pan);
    if (reinfo == NULL) {
	free(relist);
	free(hlist);
	return E_ALLOC;
    }

    if (hlist != NULL) {
	/* estimate the augmented random-effects model (GLS) */
	int *biglist = gretl_list_add(relist, hlist, &err);

	if (!err) {
	    remod = lsq(biglist, &reZ, reinfo, OLS, lsqopt);
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

#if PDEBUG
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

static void 
unit_error_variances (double *uvar, const MODEL *pmod, const DATAINFO *pdinfo,
		      int nunits, int T, const int *unit_obs)
{
    int i, t;
    double x;

    for (i=0; i<nunits; i++) {
	uvar[i] = 0.0;
	for (t=0; t<T; t++) {
	    x = pmod->uhat[panel_index(i, t)];
	    if (!na(x)) {
		uvar[i] += x * x;
	    }
	}
	if (unit_obs[i] > 1) {
	    uvar[i] /= (double) unit_obs[i]; 
	}	
    }
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
	model_test_set_pvalue(test, chisq_cdf_comp(pan->BP, 1));
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
		pan->BP, pan->BP, chisq_cdf_comp(pan->BP, 1));

	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the pooled OLS model\nis adequate, in favor of the random "
		     "effects alternative.)\n\n"));
    } 

    return 0;
}

#if PDEBUG
static void print_hausman_details (panelmod_t *pan)
{
    int i, k = pan->nbeta;

    for (i=0; i<k; i++) {
 	fprintf(stderr, "b%d_FE - beta%d_RE = %g\n", i, i, pan->bdiff[i]);
    }
    fputc('\n', stderr);

    k = k * (k + 1) / 2;

    for (i=0; i<k; i++) {
 	fprintf(stderr, "vcv_diff[%d] = %g\n", i, pan->sigma[i]);
    }
}
#endif

static int complete_hausman_test (panelmod_t *pan, PRN *prn)
{
    int err = 0;

    if (pan->bdiff != NULL && pan->sigma != NULL) {
	/* matrix approach */
#if PDEBUG
	print_hausman_details(pan);
#endif
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

static int panel_obs_accounts (const MODEL *pmod, int nunits, int T, 
			       int *effn, int *effT, 
			       int **unit_obs, int *bal)
{
    int *uobs;
    int ninc = 0;
    int obsmax = 0;
    int i, t;
    
    *effn = *effT = 0;

    uobs = malloc(nunits * sizeof *uobs);
    if (uobs == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nunits; i++) {
	uobs[i] = 0;
	for (t=0; t<T; t++) {
	    if (!na(pmod->uhat[panel_index(i, t)])) {
		uobs[i] += 1;
	    }
	}
	if (uobs[i] > 0) {
	    ninc++;
	    if (uobs[i] > obsmax) {
		obsmax = uobs[i];
	    }
	}
    }

    if (bal != NULL) {
	*bal = 1;
	for (i=0; i<nunits; i++) {
	    if (uobs[i] > 0 && uobs[i] != obsmax) {
		*bal = 0;
		break;
	    }
	}
    }

    *effn = ninc;
    *effT = obsmax;
    *unit_obs = uobs;

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

static int 
panelmod_setup (panelmod_t *pan, MODEL *pmod,
		const DATAINFO *pdinfo, gretlopt opt)
{
    int err = 0;

    panelmod_init(pan, pmod, opt);

    /* assumes (possibly padded) balanced panel dataset */
    pan->nunits = pdinfo->n / pdinfo->pd;
    pan->T = pdinfo->pd;

    panel_index_init(pdinfo, pan->nunits, pan->T);
    
    err = panel_obs_accounts(pmod, pan->nunits, pan->T,
			     &pan->effn, &pan->effT, 
			     &pan->unit_obs,
			     &pan->balanced);

    if (!err && (pan->opt & (OPT_U | OPT_F))) {
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
	/* at many points in these functions we assume the base
	   regression has an intercept included */
	return 1;
    }

    /* add OPT_V to make the fixed and random effects functions verbose */
    err = panelmod_setup(&pan, pmod, pdinfo, opt | OPT_V);
    if (err) {
	goto bailout;
    }   

    if (pan.effn < pan.nunits) {
	fprintf(stderr, "number of units included = %d\n", pan.effn);
	if (pan.effn <= 0) {
	    return E_DATA;
	}
    }

    /* figure out which of the original regressors are time-varying,
       or unit-varying as the case may be 
    */
    err = varying_vars_list((const double **) *pZ, pdinfo, &pan);
    if (err) {
	goto bailout;
    }

    err = panel_set_varying(&pan, pmod);
    if (err) {
	goto bailout;
    }    

    /* degrees of freedom relative to # of x-sectional units */
    xdf = pan.effn - pmod->ncoeff;

#if PDEBUG
    fprintf(stderr, "nunits=%d, T=%d, effn=%d, effT=%d, xdf=%d\n",
	    pan.nunits, pan.T, pan.effn, pan.effT, xdf);
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
		pan.effn, pan.effT);
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
    
    if (xdf > 0 && !na(pan.within_s2)) {
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
	} else {
	    random_effects(&pan, (const double **) *pZ, pdinfo, 
			   (const double **) gZ, prn);
	    complete_hausman_test(&pan, prn);
	}

	if (ginfo != NULL) {
	    destroy_dataset(gZ, ginfo);
	}
    }

 bailout:

    panelmod_free(&pan);

    return err;
}

/* real_panel_model:
 * @list: list containing model specification.
 * @pZ: pointer to data array.  
 * @pdinfo: data info pointer. 
 * @opt: may include %OPT_U for the random effects model,
 * %OPT_R for robust standard errors (fixed effects model
 * only), %OPT_H to use the regression approach to the Hausman
 * test (random effects only)...
 * @prn: printing struct.
 *
 * Estimates a panel model, either fixed effects or random
 * effects if %OPT_U is included.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL real_panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn)
{
    MODEL mod;
    panelmod_t pan;
    gretlopt pan_opt = opt;
    int err = 0;

    /* baseline: estimated via pooled OLS */
    mod = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (mod.errcode) {
	return mod;
    }

    if (mod.ifc == 0) {
	/* at many points in these functions we assume the base
	   regression has an intercept included */
	strcpy(gretl_errmsg, "Panel models must include an intercept");
	mod.errcode = E_DATA;
	return mod;
    }

    if (!(opt & OPT_U)) {
	/* default: add OPT_F to save the fixed effects model */
	pan_opt |= OPT_F;
    }

    err = panelmod_setup(&pan, &mod, pdinfo, pan_opt);
    if (err) {
	goto bailout;
    }   

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

    if ((opt & OPT_U) && !na(pan.within_s2)) {
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
			   (const double **) gZ, prn);
	    save_breusch_pagan_result(&pan);
	    complete_hausman_test(&pan, prn);
	}

	if (ginfo != NULL) {
	    destroy_dataset(gZ, ginfo);
	}
    }

 bailout:

    if (!err) {
	clear_model(&mod);
	mod = *pan.realmod;
    }

    panelmod_free(&pan);

    if (err && mod.errcode == 0) {
	mod.errcode = err;
    }

    return mod;    
}

static int 
write_uvar_to_dataset (double *uvar, int nunits, int T,
		       double **Z, DATAINFO *pdinfo)
{
    int i, t;
    int uv = pdinfo->v - 1;

    for (i=0; i<nunits; i++) {
	for (t=0; t<T; t++) {
	    if (uvar[i] <= 0.0) {
		Z[uv][panel_index(i, t)] = 0.0;
	    } else {
		Z[uv][panel_index(i, t)] = 1.0 / uvar[i];
	    }
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

static void
print_wald_test (double W, int nunits, const int *unit_obs, PRN *prn)
{
    int i, df = 0;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] > 1) df++;
    }

    pprintf(prn, "\n%s\n%s:\n",
	    _("Distribution free Wald test for heteroskedasticity"),
	    _("based on the FGLS residuals"));
    pprintf(prn, "%s(%d) = %g, ",  _("Chi-square"), df, W);
    pprintf(prn, _("with p-value = %g\n\n"), chisq_cdf_comp(W, df));
}

static double 
wald_hetero_test (const MODEL *pmod, const DATAINFO *pdinfo, 
		  double s2, const double *uvar,
		  int T, int nunits, const int *unit_obs)
{
    double x, W = 0.0;
    int i, t, Ti;

    for (i=0; i<nunits; i++) {
	double fii = 0.0;

	Ti = unit_obs[i];
	if (Ti == 1) {
	    W = NADBL;
	    break;
	}
	for (t=0; t<T; t++) {
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
    }

    return W;
}

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
	model_test_set_pvalue(test, chisq_cdf_comp(x2, df));
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

static double real_ll (const MODEL *pmod, const double *uvar, 
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

    return ll;
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

#define SMALLDIFF 0.0001
#define WLS_MAX   20

MODEL panel_wls_by_unit (const int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn)
{
    MODEL mdl;
    gretlopt wlsopt = OPT_A;
    double *uvar = NULL;
    double *bvec = NULL;
    double s2, diff = 1.0;
    double W = NADBL;
    int *unit_obs = NULL;
    int *wlist = NULL;
    int nunits, effn, T, effT;
    int orig_v = pdinfo->v;
    int i, iter = 0;

    gretl_errmsg_clear();

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

    nunits = pdinfo->n / pdinfo->pd;
    T = pdinfo->pd;

    panel_index_init(pdinfo, nunits, T);

    mdl.errcode = panel_obs_accounts(&mdl, nunits, T,
				     &effn, &effT,
				     &unit_obs,
				     NULL);
    if (mdl.errcode) {
	goto bailout;
    }

    uvar = malloc(nunits * sizeof *uvar);
    if (uvar == NULL) {
	free(unit_obs);
	mdl.errcode = E_ALLOC;
	return mdl;
    }  

    if (opt & OPT_T) {
	if (singleton_check(unit_obs, nunits)) {
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

    /* if wanted, iterate to ML solution; otherwise just do
       one-step FGLS estimation 
    */

    while (diff > SMALLDIFF) {

	iter++;

	unit_error_variances(uvar, &mdl, pdinfo, nunits, T, unit_obs);

	if (opt & OPT_V) {
	    if (opt & OPT_T) {
		pprintf(prn, "\n*** %s %d ***\n", _("iteration"), 
			iter);
	    } else {
		pputc(prn, '\n');
	    }
	    pputs(prn, " unit    variance\n");
	    for (i=0; i<nunits; i++) {
		if (unit_obs[i] > 0) {
		    pprintf(prn, "%5d%12g\n", i + 1, uvar[i]);
		}
	    }
	}

	if ((opt & OPT_T) && iter == 2) {
	    W = wald_hetero_test(&mdl, pdinfo, s2, uvar, T, 
				 nunits, unit_obs);
	}

	write_uvar_to_dataset(uvar, nunits, T, *pZ, pdinfo);

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
	if (!(opt & OPT_A)) {
	    set_model_id(&mdl);
	}
	gretl_model_set_int(&mdl, "n_included_units", effn);
	gretl_model_set_int(&mdl, "unit_weights", 1);
	mdl.nwt = 0;

	if (opt & OPT_T) {
	    gretl_model_set_int(&mdl, "iters", iter);
	    ml_hetero_test(&mdl, s2, uvar, nunits, unit_obs);
	    unit_error_variances(uvar, &mdl, pdinfo, nunits, T, unit_obs);
	    mdl.lnL = real_ll(&mdl, uvar, nunits, unit_obs);
	    if (opt & OPT_V) {
		pputc(prn, '\n');
	    }
	} else {
	    unit_error_variances(uvar, &mdl, pdinfo, nunits, T, unit_obs);
	    W = wald_hetero_test(&mdl, pdinfo, s2, uvar, T, 
				 nunits, unit_obs);
	}

	if (!na(W)) {
	    print_wald_test(W, nunits, unit_obs, prn);
	}
    }    

 bailout:

    free(unit_obs);
    free(uvar);
    free(wlist);
    free(bvec);

    dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
    
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

    /* basic checks */
    if (order <= 0) order = 1;
    if (order > pdinfo->pd - 1) return E_DF;
    if (pmod->ncoeff + order >= sn) return E_DF;

    if (!balanced_panel(pdinfo)) { 
	/* FIXME: check should be on balance of model? */
        return E_DATA;
    }

    if (pmod->missmask != NULL) {
	return E_MISSDATA;
    }

    /* get number of cross-sectional units */
    nunits = sn / pdinfo->pd;

    /* we lose "order" observations for each unit */
    nobs = sn - nunits * order;

    /* the required number of variables */
    nv = pmod->list[0] + order;

    /* create temporary reduced dataset */
    tmpinfo = create_new_dataset(&tmpZ, nv, nobs, 0);
    if (tmpinfo == NULL) return E_ALLOC;

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
	LMF = (aux.rsq / (1.0 - aux.rsq)) * dfd / order; 
	pval = f_cdf_comp(LMF, order, dfd);

	pprintf(prn, "\n%s: LMF = %f,\n", _("Test statistic"), LMF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, dfd, LMF, pval);

	pprintf(prn, "\n%s: TR^2 = %f,\n", 
		_("Alternative statistic"), trsq);
	pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", 	_("with p-value"), 
		_("Chi-square"), order, trsq, chisq_cdf_comp(trsq, order));

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
    int sts = (pdinfo->structure == STACKED_TIME_SERIES);
    double *tmp;
    char **markers = NULL;
    double pdx;
    int pd = pdinfo->pd;
    int nblocks = pdinfo->n / pd;
    int i, j, t;

    if (pdinfo->structure != STACKED_TIME_SERIES &&
	pdinfo->structure != STACKED_CROSS_SECTION) {
	return E_DATA;
    }

    tmp = malloc(pdinfo->n * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    /* copy the data series across in transformed order */
    for (i=1; i<pdinfo->v; i++) {
	if (var_is_scalar(pdinfo, i)) {
	    continue;
	}
	for (t=0; t<pdinfo->n; t++) {
	    tmp[t] = Z[i][t];
	}
	for (j=0; j<pd; j++) {
	    for (t=0; t<nblocks; t++) {
		Z[i][j * nblocks + t] = tmp[j + pd * t];
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
	    for (j=0; j<pd; j++) {
		for (t=0; t<nblocks; t++) {
		    strcpy(pdinfo->S[j * nblocks + t], markers[j + pd * t]);
		}
	    }
	    free_strings_array(markers, pdinfo->n);
	} else {
	    /* should we flag an error? */
	    dataset_destroy_obs_markers(pdinfo); 
	}
    }

    /* change the datainfo setup */
    pdinfo->structure = (sts)? STACKED_CROSS_SECTION : STACKED_TIME_SERIES;
    pdinfo->pd = nblocks;

    pdinfo->sd0 = 1.0;
    pdx = 0.1;
    while (nblocks /= 10) {
	pdx *= 0.1;
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
		if (na(pan->pooled->uhat[bigt])) {
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
    printlist(pan->vlist, "varying regressors");
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

static int dataset_sort_by (double **Z, DATAINFO *pdinfo,
			    int uv, int tv)
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
	for (t=0; t<n; t++) {
	    tmp[t] = Z[i][t];
	}
	for (t=0; t<n; t++) {
	    Z[i][t] = tmp[s.points[t].obsnum];
	}	
    }

    if (S != NULL) {
	for (t=0; t<n; t++) {
	    strcpy(S[t], pdinfo->S[t]);
	}
	for (t=0; t<n; t++) {
	    strcpy(pdinfo->S[t], S[s.points[t].obsnum]);
	}
	free_strings_array(S, n);
    }	

    free(s.points);
    free(tmp);

    return 0;
}

static int get_mask_skip (const char *mask, int n, int t)
{
    int s, skip = 1;

    for (s=t+1; s<n; s++) {
	if (mask[s]) skip++;
	else break;
    }

    return skip;
}

/* For when we've subsampled a panel dataset, and have padded out
   the subsample for balance: before restoring the full dataset
   we need to get rid of the added padding.
*/

int unpad_panel_dataset (double ***pZ, DATAINFO *pdinfo)
{
    PANINFO *pan = pdinfo->paninfo;
    int n = pdinfo->n;
    int i, t, s, drop;
    int skip, nrem;
    size_t srem;
    int err = 0;

    if (pan == NULL) {
	return E_DATA;
    }

    if (pan->padmask == NULL) {
	return 0;
    }

    s = 0;
    for (t=0; t<n; t++) {
	if (pan->padmask[s]) {
	    skip = get_mask_skip(pan->padmask, pdinfo->n, s);
	    nrem = n - skip - t;
	    srem = (nrem * sizeof *pan->unit);
	    memmove(pan->unit + t, pan->unit + t + skip, srem);
	    memmove(pan->period + t, pan->period + t + skip, srem);
	    srem = (nrem * sizeof ***pZ);
	    for (i=1; i<pdinfo->v; i++) {
		if (var_is_series(pdinfo, i)) {
		    memmove((*pZ)[i] + t, (*pZ)[i] + t + skip, srem);
		}
	    }
	    s += skip;
	    n -= skip;
	} 
	s++;
    }

    drop = pdinfo->n - n;

#if PDEBUG
    fprintf(stderr, "unpad: dropping %d observations\n", drop);
#endif

    if (drop > 0) {
	err = dataset_drop_observations(drop, pZ, pdinfo);
    }

    return err;
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
			      char *mask)
{
    double **bigZ = NULL;
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

    return err;
}

static int is_positive (const double *x, int n)
{
    int i;

    for (i=1; i<n; i++) {
	if (x[i] <= 0.0) {
	    return 0;
	}
    }    

    return 1;
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
    err = dataset_sort_by(Z, pdinfo, uv, tv);

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
				Z, pdinfo, mask);
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
    int n = pdinfo->n; /* ? */
    int uv, tv;
    int err = 0;

    if (!strncmp(line, "setobs", 6)) {
	line += 7;
    }

    err = uv_tv_from_line(line, pdinfo, &uv, &tv);
    if (err) {
	return err;
    }

    if (!is_positive(Z[uv], n) || !is_positive(Z[tv], n)) {
	return E_DATA;
    }

    return set_panel_structure_from_vars(uv, tv, Z, pdinfo);
}

/* used for readjusting panel data structure after sub-sampling */

int set_panel_structure_from_indices (double ***pZ, DATAINFO *pdinfo)
{
    int err;

    err = dataset_finalize_panel_indices(pdinfo);

    if (!err) {
	int t, v = pdinfo->v;

	err = dataset_add_series(2, pZ, pdinfo);
	if (!err) {
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[v][t] = pdinfo->paninfo->unit[t];
		(*pZ)[v+1][t] = pdinfo->paninfo->period[t];
	    }
	    err = set_panel_structure_from_vars(v, v+1, *pZ, pdinfo);
	    if (!err) {
		err = dataset_drop_last_variables(2, pZ, pdinfo);
	    }
	}
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
    } else if (floateq(Z[v][0], Z[v][1])) { /* "year" is same for first two obs */
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





