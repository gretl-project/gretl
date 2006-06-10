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

enum vcv_ops {
    VCV_INIT,
    VCV_SUBTRACT
};

typedef struct panelmod_t_ panelmod_t;

struct panelmod_t_ {
    int nunits;           /* total cross-sectional units */
    int effn;             /* effective (included) cross-section units */
    int T;                /* times-series length of panel */
    int effT;             /* effective times-series length */
    int ndum;             /* number of dummy variables added for FE model */
    int *unit_obs;        /* array of number of observations per x-sect unit */
    char *varying;        /* array to record properties of pooled-model regressors */
    int *vlist;           /* list of time-varying variables from pooled model */
    gretlopt opt;         /* option flags */
    int nbeta;            /* number of slope coeffs for Hausman test */
    double sigma_e;       /* fixed-effects standard error */
    double H;             /* Hausman test statistic */
    double *bdiff;        /* array of coefficient differences */
    double *sigma;        /* Hausman covariance matrix */
    double fe_var;        /* fixed-effects error variance */
    double gm_var;        /* group means error variance */
    MODEL *pooled;        /* reference model (pooled OLS) */
};

struct {
    int n;
    int T;
    int offset;
} panidx;

static int 
varying_vars_list (const double **Z, const DATAINFO *pdinfo,
		   panelmod_t *pan);

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

    pan->nbeta = 0;
    pan->sigma_e = NADBL;
    pan->H = NADBL;
    pan->bdiff = NULL;
    pan->sigma = NULL;
    
    pan->pooled = pmod;
}

static void panelmod_free (panelmod_t *pan)
{
    free(pan->unit_obs);
    free(pan->varying);
    free(pan->vlist);

    free(pan->bdiff);
    free(pan->sigma);
    
    pan->pooled = NULL;
}

/* test variable number v against the (possibly reduced) regression
   list which contains only time-varying regressors (or perhaps
   unit-varying regressors if we're doing time effects) */

static int var_is_varying (const int *list, int v)
{
    int ret = 0;

    if (v != 0) {
	int i;

	for (i=2; i<=list[0]; i++) {
	    if (list[i] == v) {
		ret = 1;
		break;
	    }
	}
    }

    return ret;
}

/* Durbin-Watson statistic for fixed effects model on a 
   balanced panel */

static void panel_dwstat (MODEL *pmod, const DATAINFO *pdinfo)
{
    int T = pdinfo->pd;
    double ut, u1;
    double num = 0.0;
    int i, t, s, n;

    if (pmod->ess <= 0.0) {
	return;
    }

    if (pmod->nobs % T != 0) {
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

static int hausman_allocate (panelmod_t *pan)
{
    int nbeta = pan->vlist[0] - 2;
    int nsigma = (nbeta * nbeta + nbeta) / 2;

    pan->nbeta = nbeta;

    /* array to hold differences between coefficient estimates */
    pan->bdiff = malloc(nbeta * sizeof *pan->bdiff);
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

#if PDEBUG
    fprintf(stderr, "hausman_allocate: alloc'd %d terms for bdiff, "
	    "%d for sigma\n", nbeta, nsigma);
#endif

    return 0;
}   

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

static void print_theta (const panelmod_t *pan, double theta,
			 int balanced, PRN *prn)
{
    pputs(prn, "Variance estimators:\n");
    pprintf(prn, "sigma_e = %14.5g, sigma_a = %14.5g \n", pan->gm_var, 
	    pan->fe_var);

    if (balanced) {
	pprintf(prn, "theta used for pseudo-differencing = %14.5g\n", theta);
    } else {
	pputs(prn, "Panel is unbalanced: theta varies across units\n");
    }
    pputc(prn, '\n');
}

/* construct a version of the dataset from which the group means
   are subtracted */

static DATAINFO *
within_groups_dataset (const double **Z, double ***wZ, panelmod_t *pan)
{
    DATAINFO *winfo;
    int wnobs = 0;
    int i, j, k;
    int t, bigt;

#if PDEBUG
    fprintf(stderr, "within_groups: nvar=%d, pooled model nobs=%d, *wZ=%p\n", 
	    pan->vlist[0], pan->pooled->nobs, (void *) *wZ);
#endif

    for (i=0; i<pan->nunits; i++) { 
	if (pan->unit_obs[i] > 1) {
	    /* we need more than one observation to compute any
	       within group variation */
	    wnobs += pan->unit_obs[i];
	}
    }

    winfo = create_new_dataset(wZ, pan->vlist[0], wnobs, 0);
    if (winfo == NULL) {
	return NULL;
    }

#if PDEBUG
    fprintf(stderr, "within_groups_dataset: now *wZ=%p (wnobs = %d)\n", 
	    (void *) *wZ, wnobs);
#endif

    k = 0;
    for (j=1; j<=pan->vlist[0]; j++) { 
	int s = 0;

	if (pan->vlist[j] == 0) {
	    continue;
	} else {
	    k++;
	}

#if PDEBUG
	fprintf(stderr, "working on list[%d] (original var %d, var %d in wZ)\n",
		j, pan->vlist[j], k);
#endif

	for (i=0; i<pan->nunits; i++) { 
	    int Ti = pan->unit_obs[i];
	    double xbar = 0.0;

#if PDEBUG
	    fprintf(stderr, "looking at x-sect unit %d\n", i);
#endif

	    if (Ti <= 1) {
#if PDEBUG
		fprintf(stderr, " skipping because Ti = %d\n", Ti);
#endif
		continue;
	    }

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
	    for (t=0; t<pan->T; t++) {
		if (s >= winfo->n) {
		    fprintf(stderr, "*** Error: overflow of wZ at s = %d!\n", s);
		    break;
		}
		bigt = panel_index(i, t);
		if (!na(pan->pooled->uhat[bigt])) {
		    (*wZ)[k][s] = Z[pan->vlist[j]][bigt] - xbar;
#if PDEBUG > 1
		    fprintf(stderr, "Set wZ[%d][%d] = %g\n", k, s, (*wZ)[k][s]);
#endif
		    s++;
		}
	    }
	}
    }

    return winfo;
}

/* construct a mini-dataset containing the group means */

static DATAINFO *
group_means_dataset (panelmod_t *pan,
		     const double **Z, const DATAINFO *pdinfo,
		     double ***gZ)
{
    DATAINFO *ginfo;
    int i, j, k;
    int t, bigt;

#if PDEBUG
    fprintf(stderr, "group_means_dataset: nvar=%d, nobs=%d, *gZ=%p\n", 
	    pan->pooled->list[0], pan->effn, (void *) *gZ);
#endif

    ginfo = create_new_dataset(gZ, pan->pooled->list[0], pan->effn, 0);
    if (ginfo == NULL) {
	return NULL;
    }

    k = 1;
    for (j=1; j<=pan->pooled->list[0]; j++) { 
	int s = 0;

	if (pan->pooled->list[j] == 0) {
	    continue;
	}

	for (i=0; i<pan->nunits; i++) { 
	    int Ti = pan->unit_obs[i];
	    double xx = 0.0;

	    if (Ti == 0) {
		continue;
	    }

	    for (t=0; t<pan->T; t++) {
		bigt = panel_index(i, t);
		if (!na(pan->pooled->uhat[bigt])) {
		    xx += Z[pan->pooled->list[j]][bigt];
		}
	    }

	    (*gZ)[k][s] = xx / (double) Ti;
#if PDEBUG > 1
	    fprintf(stderr, "Set gZ[%d][%d] = %g\n", k, i, (*gZ)[k][s]);
#endif
	    s++;
	}
	k++;
    }

    return ginfo;
}

/* calculate the group means regression and its error variance */

static int
group_means_variance (panelmod_t *pan, double ***gZ, DATAINFO *ginfo)
{
    MODEL gmmod;
    int *gmlist;
    int i, j;
    int err = 0;

    gmlist = gretl_list_new(pan->pooled->list[0]);
    if (gmlist == NULL) {
	return E_ALLOC;
    }

    j = 1;
    for (i=1; i<=gmlist[0]; i++) { 
	if (pan->pooled->list[i] == 0) {
	    gmlist[i] = 0;
	} else {
	    gmlist[i] = j++;
	}
    }

    gmmod = lsq(gmlist, gZ, ginfo, OLS, OPT_A | OPT_Z);

    if (gmmod.errcode == 0) {
	pan->gm_var = gmmod.sigma * gmmod.sigma;
    } else {
	err = gmmod.errcode;
#if PDEBUG
	fprintf(stderr, "error %d estimating group_means_variance\n", err);
#endif
    }

    clear_model(&gmmod);
    free(gmlist);

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

/* fill out the covariance matrix for use with the Hausman test:
   the entries that get transcribed here are only those for
   slopes with respect to time-varying variables */

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

/* calculate Hausman test statistic */

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
#if PDEBUG
	    fprintf(stderr, "added %g * %g to pan->H: now = %g\n",
		    x[i], pan->bdiff[i], pan->H);
#endif

	}
	if (pan->H < 0.0) {
	    pan->H = NADBL;
	}
    }

    free(x);
    free(ipiv);

    return info;
}

static void apply_panel_df_correction (MODEL *pmod, int ndf)
{
    double dfcorr = sqrt((double) pmod->dfd / (pmod->dfd - ndf));
    int i, j, idx;

    pmod->dfd -= ndf;
    pmod->dfn += ndf - 1; 

    pmod->sigma *= dfcorr;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->sderr[i] *= dfcorr;
    }

    dfcorr *= dfcorr;

    if (pmod->vcv != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    for (j=i; j<pmod->ncoeff; j++) {
		idx = ijton(i, j, pmod->ncoeff);
		pmod->vcv[idx] *= dfcorr;
	    }
	}
    }
}

#define MINOBS 1

/* calculate the fixed effects regression, either using dummy variables
   or by subtracting the group means from the data
*/

static MODEL
fixed_effects_model (panelmod_t *pan, double ***pZ, DATAINFO *pdinfo)
{
    int i, t, oldv = pdinfo->v;
    int dvlen, ndum = 0;
    int *felist = NULL;
    MODEL lsdv;
#if PDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
#endif

    gretl_model_init(&lsdv);

    for (i=0; i<pan->nunits; i++) {
	if (pan->unit_obs[i] >= MINOBS) {
	    ndum++;
	}
    }

    dvlen = pan->vlist[0] + ndum;
    ndum--;

#if PDEBUG
    fprintf(stderr, "ndum = %d, dvlen = pan->vlist[0] + ndum = %d + %d = %d\n", 
	    ndum + 1, pan->vlist[0], ndum, dvlen);
#endif

    /* should we use a set of dummy variables, or subtract
       the group means? */
    if (getenv("NODUMMIES") != NULL) {
	pan->ndum = 0;
    } else if (ndum <= 20 || ndum < pan->vlist[0]) {
	pan->ndum = ndum;
    }

#if PDEBUG
    if (pan->ndum > 0) {
	fprintf(stderr, " using dummy variables approach\n");
    } else {
	fprintf(stderr, " subtracting group means\n");
    }
#endif

    if (pan->ndum == 0) {
	double **wZ = NULL;
	DATAINFO *winfo;

	felist = gretl_list_new(pan->vlist[0] - 1);
	if (felist == NULL) {
	    lsdv.errcode = E_ALLOC;
	    return lsdv;
	}	

	winfo = within_groups_dataset((const double **) *pZ, &wZ, pan);
	if (winfo == NULL) {
	    free(felist);
	    lsdv.errcode = E_ALLOC;
	    return lsdv;
	} 

	for (i=1; i<=felist[0]; i++) {
	    felist[i] = i;
	}

	/* OPT_Z: do _not_ automatically eliminate perfectly
	   collinear variables */
	lsdv = lsq(felist, &wZ, winfo, OLS, OPT_A | OPT_Z);

	if (lsdv.errcode) {
	    ; /* pass on */
	} else if (lsdv.list[0] < felist[0]) {
	    /* one or more variables were dropped, because they were
	       all zero -- this is a symptom of collinearity */
	    lsdv.errcode = E_SINGULAR;
	} else {
	    /* we estimated a bunch of group means, and have to
	       subtract that many degrees of freedom */
	    apply_panel_df_correction(&lsdv, pan->nunits);
#if PDEBUG
	    printmodel(&lsdv, winfo, OPT_O, prn);
#endif
	}

	destroy_dataset(wZ, winfo);
	    
    } else {
	int j = 0;

	/* We can be assured there's an intercept in the original
	   regression */

	felist = gretl_list_new(dvlen - 1);
	if (felist == NULL) {
	    lsdv.errcode = E_ALLOC;
	    return lsdv;
	}

	if (dataset_add_series(ndum, pZ, pdinfo)) {
	    free(felist);
	    lsdv.errcode = E_ALLOC;
	    return lsdv;
	}

	/* create the per-unit dummy variables */

	for (i=0; i<pan->nunits; i++) {
	    int dv = oldv + j;

	    if (pan->unit_obs[i] < MINOBS) {
		continue;
	    }

	    sprintf(pdinfo->varname[dv], "unit_%d", i + 1);

	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[dv][t] = 0.0;
	    }

	    for (t=0; t<pan->T; t++) {
		(*pZ)[dv][panel_index(i, t)] = 1.0;
	    }

	    if (++j == ndum) {
		break;
	    }
	}

	/* construct the regression list */

	for (i=1; i<=pan->vlist[0]; i++) {
	    felist[i] = pan->vlist[i];
	}

	for (i=0; i<ndum; i++) {
	    felist[pan->vlist[0] + i + 1] = oldv + i;
	}

	lsdv = lsq(felist, pZ, pdinfo, OLS, OPT_A | OPT_Z);

#if PDEBUG
	if (!lsdv.errcode) {
	    printmodel(&lsdv, pdinfo, OPT_NONE, prn);
	}
#endif

	dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo);

    }

#if PDEBUG
    gretl_print_destroy(prn);
#endif

    free(felist);

    return lsdv;
}

static int print_fe_results (panelmod_t *pan, 
			     MODEL *lsdv,
			     DATAINFO *pdinfo,
			     PRN *prn)
{
    int i, j, k;
    int ndum;
    double F;

#if PDEBUG
    fprintf(stderr, "f.e. variance, lsdv.ncoeff=%d, var=%g, "
	    "list[0]=%d\n", lsdv->ncoeff, pan->fe_var, pan->vlist[0]);
    fprintf(stderr, "haus->bdiff = %p, haus->sigma = %p\n",
	    pan->bdiff, pan->sigma);
#endif

    ndum = lsdv->list[0] - pan->vlist[0];

    pputs(prn, 
	  _("Fixed effects estimator\n"
	    "allows for differing intercepts by cross-sectional unit\n"
	    "slope standard errors in parentheses, p-values in brackets\n"));

    if (ndum > 0) {
	pputs(prn, _("a_i = intercepts"));
	pputs(prn, "\n\n");
    } else {
	pputc(prn, '\n');
    }

    /* print the slope coefficients, for varying regressors */
    k = 0;
    for (i=1; i<pan->vlist[0] - 1; i++) {
	int vi = pan->vlist[i+2];

	j = (pan->ndum > 0)? i : i - 1;
	print_panel_coeff(lsdv, pdinfo->varname[vi], j, prn);
    }

    pputc(prn, '\n');   

    /* if we used dummy variables, print the per-unit intercept 
       estimates */

    if (ndum > 0) {
	j = 0;
	for (i=0; i<pan->nunits; i++) {
	    char dumstr[VNAMELEN];
	    double b;

	    if (pan->unit_obs[i] < MINOBS) {
		continue;
	    }

	    if (j == ndum) {
		b = lsdv->coeff[0];
	    } else {
		b = lsdv->coeff[j + pan->vlist[0] - 1] + lsdv->coeff[0];
	    }
	    sprintf(dumstr, "a_%d", i + 1);
	    pprintf(prn, "%*s: %14.4g\n", VNAMELEN, dumstr, b);
	    j++;
	}
    } else {
	pprintf(prn, _("%d group means were subtracted from the data"), pan->effn);
	pputc(prn, '\n');
	ndum = pan->effn - 1;
    }

    pprintf(prn, _("\nResidual variance: %g/(%d - %d) = %g\n"), 
	    lsdv->ess, lsdv->nobs, pan->vlist[0] - 1 + ndum, pan->fe_var);

    if (pan->ndum > 0) {
	pputs(prn, _("Joint significance of unit dummy variables:\n"));
    } else {
	pprintf(prn, _("Joint significance of differing group means:\n"));
    }

    F = (pan->pooled->ess - lsdv->ess) * lsdv->dfd / (lsdv->ess * ndum);

    pprintf(prn, " F(%d, %d) = %g %s %g\n", ndum, lsdv->dfd, F, 
	    _("with p-value"), fdist(F, ndum, lsdv->dfd));

    pputs(prn, _("(A low p-value counts against the null hypothesis that "
		 "the pooled OLS model\nis adequate, in favor of the fixed "
		 "effects alternative.)\n\n"));

    return 0;
}

static void fix_panelmod_list (MODEL *targ, MODEL *src, panelmod_t *pan)
{
    int i;

#if PDEBUG
    printlist(targ->list, "targ->list");
    printlist(src->list, "src->list");
    printlist(pan->vlist, "pan->vlist");
#endif

    free(targ->list);
    targ->list = src->list;
    src->list = NULL;

    /* remove any non-varying variables */

    for (i=2; i<=targ->list[0]; i++) {
	if (!in_gretl_list(pan->vlist, targ->list[i])) {
	    gretl_list_delete_at_pos(targ->list, i--);
	}
    }

#if PDEBUG
    printlist(targ->list, "new targ->list");
#endif

}

static void fix_within_stats (MODEL *targ, MODEL *src, panelmod_t *pan)
{
    int nc = targ->ncoeff;

    /* FIXME reporting of const */

    fix_panelmod_list(targ, src, pan);

    targ->ybar = src->ybar;
    targ->sdy = src->sdy;
    targ->ifc = 1;

    targ->rsq = 1.0 - (targ->ess / src->tss);

    if (targ->dfd > 0) {
	double den = targ->tss * targ->dfd;

	targ->adjrsq = 1 - (targ->ess * (targ->nobs - 1) / den);
    }

    if (targ->rsq < 0.0) {
	targ->rsq = 0.0;
    } else {
	targ->fstt = (targ->rsq / (1.0 - targ->rsq)) * ((double) targ->dfd / targ->dfn);
    }

    targ->ncoeff = targ->dfn + 1;
    ls_criteria(targ, NULL, NULL);
    targ->ncoeff = nc;
}

/* drive the calculation of the fixed effects regression, print the
   results (if wanted), and compute the error variance */

static int
fixed_effects_variance (panelmod_t *pan,
			double ***pZ, DATAINFO *pdinfo,
			PRN *prn)
{
    MODEL lsdv;
    int err = 0;

    lsdv = fixed_effects_model(pan, pZ, pdinfo);

    if (lsdv.errcode) {
	pputs(prn, _("Error estimating fixed effects model\n"));
	errmsg(lsdv.errcode, prn);
	err = lsdv.errcode;
    } else {
	int i, j, k;

	pan->fe_var = lsdv.sigma * lsdv.sigma;

	if (pan->opt & OPT_V) {
	    print_fe_results(pan, &lsdv, pdinfo, prn);
	}

	if (pan->bdiff != NULL) {
	    /* record the slope coefficients for varying regressors */
	    k = 0;
	    for (i=1; i<pan->vlist[0] - 1; i++) {
		j = (pan->ndum > 0)? i : i - 1;
		pan->bdiff[k++] = lsdv.coeff[j];
	    }
	}

	if (pan->sigma != NULL) {
	    makevcv(&lsdv);
	    pan->sigma_e = lsdv.sigma;
	    vcv_slopes(pan, &lsdv, VCV_INIT);
	}

	if (pan->opt & OPT_F) {
	    /* Swap content of fixed effects model into the MODEL
	       struct attached to pan. 
	    */
	    MODEL tmp;
		
	    tmp = *pan->pooled;
	    *pan->pooled = lsdv;
	    lsdv = tmp;
	    pan->pooled->ci = PANEL;
	    gretl_model_set_int(pan->pooled, "fixed-effects", 1);
	    if (pan->ndum == 0) {
		fix_within_stats(pan->pooled, &lsdv, pan);
	    } 
	    gretl_model_add_panel_varnames(pan->pooled, pdinfo);
	    if (pan->ndum > 0) {
		pan->pooled->list[0] -= pan->ndum;
	    }
	    set_model_id(pan->pooled);
	    panel_dwstat(pan->pooled, pdinfo);
	}
    }

    clear_model(&lsdv);

    return err;
}

/* "src" is the orginal pooled OLS model; "targ" is the GLS model
   which we'll print out when we're done.
*/

static void fix_gls_stats (MODEL *targ, MODEL *src, panelmod_t *pan,
			   const double **Z)
{
    const double *y;
    int i, t, nc;
    double yht, xit;

    fix_panelmod_list(targ, src, pan);

    targ->ybar = src->ybar;
    targ->sdy = src->sdy;

    y = Z[targ->list[1]];

#if 1 /* FIXME! The following is broken (isn't it?) */
    targ->ess = 0.0;
    for (t=0; t<src->full_n; t++) {
	if (na(y[t])) {
	    targ->uhat[t] = NADBL;
	    targ->yhat[t] = NADBL;
	    continue;
	}
	yht = 0.0;
	for (i=0; i<targ->ncoeff; i++) {
	    xit = Z[targ->list[i+2]][t]; /* vlist?? */
	    if (na(xit)) {
		yht = NADBL;
		break;
	    } else {
		yht += targ->coeff[i] * xit;
	    }
	}
	targ->yhat[t] = yht;
	if (!na(yht)) {
	    targ->uhat[t] = y[t] - yht;
	    targ->ess += targ->uhat[t] * targ->uhat[t];
	}
    } 
#endif

    targ->rsq = 1.0 - (targ->ess / src->tss);

    if (targ->dfd > 0) {
	double den = targ->tss * targ->dfd;

	targ->adjrsq = 1 - (targ->ess * (targ->nobs - 1) / den);
    }

    if (targ->rsq < 0.0) {
	targ->rsq = 0.0;
    } else {
	targ->fstt = (targ->rsq / (1.0 - targ->rsq)) * ((double) targ->dfd / targ->dfn);
    }

    nc = targ->ncoeff;
    targ->ncoeff = targ->dfn + 1;
    ls_criteria(targ, NULL, NULL);
    targ->ncoeff = nc;
}

/* calculate the random effects regression and print the results */

static int random_effects (panelmod_t *pan, 
			   const double **Z, DATAINFO *pdinfo, 
			   const double **gZ, PRN *prn)
{
    double **reZ;
    DATAINFO *reinfo;
    MODEL remod;
    double theta;
    int *relist;
    int re_n = pan->T * pan->effn;
    int i, j, k, t;
    int err = 0;

    /* regression list */
    relist = gretl_list_new(pan->pooled->list[0]);
    if (relist == NULL) {
	return E_ALLOC;
    }

    /* special transformed dataset */
    reinfo = create_new_dataset(&reZ, pan->pooled->list[0], re_n, 0);
    if (reinfo == NULL) {
	free(relist);
	return E_ALLOC;
    }

#if PDEBUG
    fprintf(stderr, "reZ: series length = T * effn = %d * %d = %d\n",
	    pan->T, pan->effn, pan->T * pan->effn);
#endif

    /* Create transformed variables: original data minus theta
       times the appropriate group mean.  Note: for unbalanced
       panels, theta varies across the units.
    */

    theta = 1.0 - sqrt(pan->fe_var / (pan->effT * pan->gm_var));

    k = 0;
    for (j=1; j<=relist[0]; j++) {
	const double *xj = Z[pan->pooled->list[j]];
	const double *gm = NULL;
	double theta_i;
	int bigt, rt, u = 0;

	if (pan->pooled->list[j] == 0) {
	    relist[j] = 0;
	} else {
	    k++;
	    relist[j] = k;
	    gm = gZ[k];
	    strcpy(reinfo->varname[k], pdinfo->varname[pan->pooled->list[j]]);
	}

	for (i=0; i<pan->nunits; i++) {
	    int Ti = pan->unit_obs[i];

	    if (Ti == 0) {
		continue;
	    }

	    if (Ti != pan->effT) {
		theta_i = 1.0 - sqrt(pan->fe_var / (Ti * pan->gm_var));
	    } else {
		theta_i = theta;
	    }

	    for (t=0; t<pan->T; t++) {
		bigt = panel_index(i, t);
		rt = u * pan->T + t;
		if (relist[j] == 0) {
		    /* the intercept */
		    reZ[0][rt] -= theta_i;
		} else if (na(pan->pooled->uhat[bigt])) {
		    reZ[k][rt] = NADBL;
		} else {
		    reZ[k][rt] = xj[bigt] - theta_i * gm[u];
		}
	    }
	    u++;
	}
    }

    remod = lsq(relist, &reZ, reinfo, OLS, OPT_A | OPT_Z);

    if ((err = remod.errcode)) {
	pputs(prn, _("Error estimating random effects model\n"));
	errmsg(err, prn);
    } else {
	if (pan->opt & OPT_V) {
	    pputs(prn,
		  _("                         Random effects estimator\n"
		    "           allows for a unit-specific component to the "
		    "error term\n"
		    "           (standard errors in parentheses, p-values in brackets)\n\n"));

	    print_theta(pan, theta, balanced_panel(pdinfo), prn);
	}

	k = 0;
	for (i=0; i<remod.ncoeff; i++) {
	    int vi = pan->pooled->list[i+2];

	    if (pan->opt & OPT_V) {
		print_panel_coeff(&remod, pdinfo->varname[vi], i, prn);
	    }
	    if (pan->bdiff != NULL && var_is_varying(pan->vlist, vi)) {
		pan->bdiff[k++] -= remod.coeff[i];
	    }
	}
	makevcv(&remod);
	vcv_slopes(pan, &remod, VCV_SUBTRACT);
    }

    if (pan->opt & OPT_R) {
	/* Swap content of random effects model into the MODEL
	   struct attached to pan.  
	*/
	MODEL tmp;

#if 0
	printmodel(&remod, reinfo, OPT_NONE, prn);
#endif		
	tmp = *pan->pooled;
	*pan->pooled = remod;
	remod = tmp;
	pan->pooled->ci = PANEL;
	gretl_model_add_panel_varnames(pan->pooled, reinfo);
	fix_gls_stats(pan->pooled, &remod, pan, Z);
	set_model_id(pan->pooled);
    }

    clear_model(&remod);
    destroy_dataset(reZ, reinfo);
    free(relist);    

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

static int 
breusch_pagan_LM (panelmod_t *pan, const DATAINFO *pdinfo, PRN *prn)
{
    double LM, A = 0.0;
    int n = pan->pooled->nobs;
    int i, t, M = 0;

    pputs(prn, _("\nMeans of pooled OLS residuals for cross-sectional "
		 "units:\n\n"));

    for (i=0; i<pan->nunits; i++) {
	int Ti = pan->unit_obs[i];
	double u, usum = 0.0;

	if (Ti == 0) {
	    continue;
	}

	for (t=0; t<pan->T; t++) {
	    u = pan->pooled->uhat[panel_index(i, t)];
	    if (!na(u)) {
		usum += u;
	    }
	}

	pprintf(prn, _(" unit %2d: %13.5g\n"), i + 1, usum / (double) Ti);
	A += usum * usum;
	M += Ti * Ti;
    }

    LM = (n * n /(2.0 * (M - n))) * pow((A / pan->pooled->ess) - 1.0, 2);

    pprintf(prn, _("\nBreusch-Pagan test statistic:\n"
		   " LM = %g with p-value = prob(chi-square(1) > %g) = %g\n"), 
	    LM, LM, chisq(LM, 1));

    pputs(prn, _("(A low p-value counts against the null hypothesis that "
		 "the pooled OLS model\nis adequate.)\n\n"));

    return 0;
}

static int do_hausman_test (panelmod_t *pan, PRN *prn)
{
#if PDEBUG
    int i, ns = pan->nbeta;
    int nterms = (ns * ns + ns) / 2;

    for (i=0; i<ns; i++) {
 	fprintf(stderr, "b%d_FE - beta%d_RE = %g\n", i, i, pan->bdiff[i]);
    }
    fputc('\n', stderr);

    for (i=0; i<nterms; i++) {
 	fprintf(stderr, "vcv_diff[%d] = %g\n", i, pan->sigma[i]);
    }
#endif

    if (bXb(pan)) { 
	pputs(prn, _("Error attempting to invert vcv difference matrix\n"));
	return 1;
    }

    if (na(pan->H)) {
	pputs(prn, _("\nHausman test matrix is not positive definite (this "
		     "result may be treated as\n\"fail to reject\" the random effects "
		     "specification).\n"));
    } else {
	pprintf(prn, _("\nHausman test statistic:\n"
		       " H = %g with p-value = prob(chi-square(%d) > %g) = %g\n"),
		pan->H, pan->nbeta, pan->H, chisq(pan->H, pan->nbeta));
	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the random effects\nmodel is consistent, in favor of the fixed "
		     "effects model.)\n"));
    }

    return 0;
}

static int get_maj_min (const DATAINFO *pdinfo, int *maj, int *min)
{
    int startmaj, startmin;
    int endmaj, endmin;

    if (sscanf(pdinfo->stobs, "%d:%d", &startmaj, &startmin) != 2) {
	return 1;
    }

    if (sscanf(pdinfo->endobs, "%d:%d", &endmaj, &endmin) != 2) {
	return 1;
    } 

    *maj = endmaj - startmaj + 1;
    *min = endmin - startmin + 1;

    return 0;
}

/* Based on the residuals from pooled OLS, determine how many
   cross-sectional units were actually included in the regression
   (after omitting any missing values); in addition, determine the
   number of observations included for each unit.
*/

static int n_included_units (const MODEL *pmod, const DATAINFO *pdinfo,
			     int *unit_obs)
{
    int nunits, T;
    int i, t;
    int ninc = 0;

    if (get_maj_min(pdinfo, &nunits, &T)) {
	return -1;
    }

    for (i=0; i<nunits; i++) {
	unit_obs[i] = 0;
	for (t=0; t<T; t++) {
	    if (!na(pmod->uhat[panel_index(i, t)])) {
		unit_obs[i] += 1;
	    }
	}
	if (unit_obs[i] > 0) {
	    ninc++;
	}
    }

    return ninc;
}

/* Determine the maximum number of observations included for any
   cross-sectional unit.  (This may vary across units due to omission
   of missing values.)
*/

static int effective_T (const int *unit_obs, int nunits)
{
    int i, effT = 0;

    for (i=0; i<nunits; i++) {
	if (unit_obs[i] > effT) {
	    effT = unit_obs[i];
	}
    }

    return effT;
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
	if (var_is_varying(pan->vlist, pmod->list[i+2])) {
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

    pan->nunits = pdinfo->n / pdinfo->pd;
    pan->T = pdinfo->pd;

    panel_index_init(pdinfo, pan->nunits, pan->T);
    
    pan->unit_obs = malloc(pan->nunits * sizeof *pan->unit_obs);
    if (pan->unit_obs == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	pan->effn = n_included_units(pmod, pdinfo, pan->unit_obs);
	pan->effT = effective_T(pan->unit_obs, pan->nunits);
    }

    return err;
}

int panel_diagnostics (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    int unbal = gretl_model_get_int(pmod, "unbalanced");
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
    fprintf(stderr, "nunits=%d, T=%d, effn=%d, effT=%d, unbal=%d, xdf=%d\n",
	    pan.nunits, pan.T, pan.effn, pan.effT, unbal, xdf);
#endif

    /* can we do the Hausman test or not? */
    if (xdf > 0) {
	err = hausman_allocate(&pan);
	if (err) {
	    goto bailout;
	}
    } 

    if (!unbal) {
	pprintf(prn, _("      Diagnostics: assuming a balanced panel with %d "
		       "cross-sectional units\n "
		       "                        observed over %d periods\n\n"), 
		pan.effn, pan.effT);
    }

    err = fixed_effects_variance(&pan, pZ, pdinfo, prn);
    if (err) {
	goto bailout;
    }

    breusch_pagan_LM(&pan, pdinfo, prn);

    if (xdf <= 0) {
	pprintf(prn, "Omitting group means regression: "
		"insufficient degrees of freedom\n");
	goto bailout;
    }
    
    if (xdf > 0 && !na(pan.fe_var)) {
	double **gZ = NULL;
	DATAINFO *ginfo;

	ginfo = group_means_dataset(&pan, (const double **) *pZ,
				    pdinfo, &gZ);

	if (ginfo != NULL) {
	    err = group_means_variance(&pan, &gZ, ginfo);
	}

	if (err) { 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	} else {
	    pprintf(prn, _("Residual variance for group means "
			   "regression: %g\n\n"), pan.gm_var);    
	    random_effects(&pan, (const double **) *pZ, pdinfo, 
			   (const double **) gZ, prn);
	    do_hausman_test(&pan, prn);
	}

	if (ginfo != NULL) {
	    destroy_dataset(gZ, ginfo);
	}
    }

 bailout:

    panelmod_free(&pan);

    return err;
}

MODEL real_panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn)
{
    MODEL mod;
    panelmod_t pan;
    gretlopt pan_opt = opt;
    int unbal;
    int err = 0;

    mod = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (mod.errcode) {
	return mod;
    }

    unbal = gretl_model_get_int(&mod, "unbalanced");

    if (mod.ifc == 0) {
	/* at many points in these functions we assume the base
	   regression has an intercept included */
	mod.errcode = E_DATA; /* FIXME error reporting */
	return mod;
    }

    if (!(opt & OPT_R) && !(opt & OPT_W)) {
	/* OPT_F: save the fixed effects model */
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

    if (opt & OPT_R) {
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

    err = fixed_effects_variance(&pan, pZ, pdinfo, prn);
    if (err) {
	goto bailout;
    }

    if ((opt & OPT_R) && !na(pan.fe_var)) {
	double **gZ = NULL;
	DATAINFO *ginfo;

	ginfo = group_means_dataset(&pan, (const double **) *pZ,
				    pdinfo, &gZ);

	if (ginfo != NULL) {
	    err = group_means_variance(&pan, &gZ, ginfo);
	}

	if (err) { 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	} else {
	    pprintf(prn, _("Residual variance for group means "
			   "regression: %g\n\n"), pan.gm_var);    
	    random_effects(&pan, (const double **) *pZ, pdinfo, 
			   (const double **) gZ, prn);
#if 0 /* should be printed after the model */
	    do_hausman_test(&pan, prn);
#endif
	}

	if (ginfo != NULL) {
	    destroy_dataset(gZ, ginfo);
	}
    }

 bailout:

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
    pprintf(prn, _("with p-value = %g\n\n"), chisq(W, df));
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
	model_test_set_pvalue(test, chisq(x2, df));
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

    gretl_model_init(&mdl);

    nunits = pdinfo->n / pdinfo->pd;
    T = pdinfo->pd;

    panel_index_init(pdinfo, nunits, T);

    unit_obs = malloc(nunits * sizeof *unit_obs);
    if (unit_obs == NULL) {
	mdl.errcode = E_ALLOC;
	return mdl;
    }

    uvar = malloc(nunits * sizeof *uvar);
    if (unit_obs == NULL) {
	free(unit_obs);
	mdl.errcode = E_ALLOC;
	return mdl;
    }    
    
    mdl = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (mdl.errcode) {
	goto bailout;
    }

    effn = n_included_units(&mdl, pdinfo, unit_obs);
    effT = effective_T(unit_obs, nunits);  

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
	pval = fdist(LMF, order, dfd);

	pprintf(prn, "\n%s: LMF = %f,\n", _("Test statistic"), LMF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, dfd, LMF, pval);

	pprintf(prn, "\n%s: TR^2 = %f,\n", 
		_("Alternative statistic"), trsq);
	pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", 	_("with p-value"), 
		_("Chi-square"), order, trsq, chisq(trsq, order));

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

	if (pan->opt & OPT_T) {
	    /* doing time effects */
	    for (t=0; t<pan->T; t++) { 
		int started = 0;
		double xval = NADBL;

		for (i=0; i<pan->nunits; i++) {
		    if (pan->unit_obs[i] == 0) {
			continue;
		    }
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
	} else {
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
	}

	if (varies) {
	    pan->vlist[k++] = vj;
	    pan->vlist[0] += 1;
	} else {
	    fprintf(stderr, "Variable %d '%s' is %s-invariant\n",
		    vj, pdinfo->varname[vj], (pan->opt & OPT_T)? 
		    "unit" : "time");
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

static void shift_data_forward (double **Z, DATAINFO *pdinfo,
				int t, int unit, const double *tid,
				int tmiss, int uv, int tv,
				int shift)
{
    int i, s;

    for (i=1; i<pdinfo->v; i++) {
	for (s=pdinfo->n-1; s>=t+shift; s--) {
	    Z[i][s] = Z[i][s - shift];
	}
	for (s=t; s<t+shift; s++) {
	    if (i == uv) {
		Z[i][s] = unit;
	    } else if (i == tv) {
		Z[i][s] = tid[tmiss++];
	    } else {
		Z[i][s] = NADBL;
	    }
	}
    }
}

static int really_pad_dataset (const double *uid, int uv, int nunits,
			       const double *tid, int tv, int nperiods, 
			       const int *uobs, double **Z, 
			       DATAINFO *pdinfo)
{
    int ni, shift;
    int i, j, t = 0;

    for (i=0; i<nunits; i++) {
	if (uobs[i] == nperiods) {
	    /* no missing obs for this unit */
	    t += nperiods;
	    continue;
	}
	ni = uobs[i];
	for (j=0; j<nperiods; ) {
	    if (j > ni || Z[tv][t] != tid[j]) {
		if (j > ni) {
		    shift = nperiods - ni;
		} else {
		    shift = 1; /* this is lame */
		}
		shift_data_forward(Z, pdinfo, t, uid[i], 
				   tid, j, uv, tv, shift);
		ni += shift;
		j += shift;
		t += shift;
	    } else {
		j++;
		t++;
	    }
	}
    }

    return 0;
}

static int maybe_pad_dataset (const double *uid, int uv, int nunits,
			      const double *tid, int tv, int nperiods, 
			      double ***pZ, DATAINFO *pdinfo)
{
    int *uobs;
    int totmiss = 0;
    int i, t;
    int err = 0;

    uobs = malloc(nunits * sizeof *uobs);
    if (uobs == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nunits; i++) {
	uobs[i] = 0;
	for (t=0; t<pdinfo->n; t++) {
	    if ((*pZ)[uv][t] == uid[i]) {
		uobs[i] += 1;
	    }
	}
#if PDEBUG
	fprintf(stderr, "unit %d: found %d obs\n", i, uobs[i]);
#endif
	totmiss += nperiods - uobs[i];
    }

#if PDEBUG
    fprintf(stderr, "calculated %d total missing obs\n", totmiss);
#endif

    if (totmiss > 0) {
	err = dataset_add_observations(totmiss, pZ, pdinfo, OPT_NONE);
	if (!err) {
	    really_pad_dataset(uid, uv, nunits,
			       tid, tv, nperiods,
			       uobs, *pZ, pdinfo);
	}
    }
		
    free(uobs);

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

int set_panel_structure_from_vars (int uv, int tv, 
				   double ***pZ, 
				   DATAINFO *pdinfo)
{
    double *uid = NULL;
    double *tid = NULL;
    int n = pdinfo->n;
    int nunits = 0;
    int nperiods = 0;
    int err = 0;

    /* FIXME sub-sampled dataset?? */

#if PDEBUG
    fprintf(stderr, "set_panel_structure_from_vars:\n "
	    "using var %d ('%s') for unit, var %d ('%s') for time\n",
	    uv, pdinfo->varname[uv], tv, pdinfo->varname[tv]);
#endif

    uid = copyvec((*pZ)[uv], n);
    tid = copyvec((*pZ)[tv], n);

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

#if PDEBUG
    fprintf(stderr, "Found %d units, %d periods\n", nunits, nperiods);
    fprintf(stderr, "Units: min %g, max %g\n", uid[0], uid[n - 1]);
    fprintf(stderr, "Periods: min %g, max %g\n", tid[0], tid[n - 1]);
#endif

    /* sort full dataset by unit and period */
    err = dataset_sort_by(*pZ, pdinfo, uv, tv);

    if (!err) {
	/* in case of unbalanced panel, pad out with missing
	   values */
	rearrange_id_array(uid, nunits, n);
	rearrange_id_array(tid, nperiods, n);
	err = maybe_pad_dataset(uid, uv, nunits, 
				tid, tv, nperiods, 
				pZ, pdinfo);
    }

    if (!err) {
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

 bailout:

    free(uid);
    free(tid);

    return err;
}

int set_panel_structure_from_line (const char *line, 
				   double ***pZ, 
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

    if (!is_positive((*pZ)[uv], n) || !is_positive((*pZ)[tv], n)) {
	return E_DATA;
    }

    return set_panel_structure_from_vars(uv, tv, pZ, pdinfo);
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




