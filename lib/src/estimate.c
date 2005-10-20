/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* estimate.c - gretl estimation procedures */

#include "libgretl.h"
#include "qr_estimate.h"
#include "libset.h"
#include "compat.h"
#include "missing_private.h"

/* There's a balancing act with 'TINY' here.  It's the minimum value
   for test that libgretl will accept before rejecting a
   data matrix as too highly collinear.  If you set it too high,
   data sets for which gretl could produce reasonable estimates will
   be rejected.  If you set it too low (and even 100 * DBL_EPSILON
   is definitely too low), gretl will produce more or less worthless
   coefficient estimates when given highly collinear data.  If you're
   tempted to change the value of TINY, check how gretl does on the
   NIST reference data sets for linear regression and ensure you're
   not getting any garbage results.  The setting of 2.1e-09 enables
   me to get decent results on the NIST nonlinear regression test
   suite, but it could be a bit too low for some contexts.
*/

#define TINY      2.1e-09 /* was 4.75e-09 (last changed 2004/07/16) */
#define SMALL     1.0e-08 /* threshold for printing a warning for collinearity */
#define YBARZERO  0.5e-14 /* threshold for treating mean of dependent
			     variable as effectively zero */
#define ESSZERO   1.0e-22 /* threshold for considering a tiny error-sum-of-
			     squares value to be effectively zero */

/* define for lots of debugging info */
#undef XPX_DEBUG

/* private function prototypes */
static int form_xpxxpy (const int *list, int t1, int t2, 
			double **Z, int nwt, double rho, int pwe,
			double *xpx, double *xpy, const char *mask);
static void regress (MODEL *pmod, double *xpy, double **Z, 
		     int n, double rho);
static int cholbeta (double *xpx, double *xpy, double *coeff, double *rss,
		     int nv);
static void diaginv (double *xpx, double *xpy, double *diag, int nv);

static double dwstat (int order, MODEL *pmod, double **Z);
static double rhohat (int order, int t1, int t2, const double *uhat);
static int hatvar (MODEL *pmod, int n, double **Z);
static void omitzero (MODEL *pmod, const double **Z, const DATAINFO *pdinfo);
static void tsls_omitzero (int *list, const double **Z, int t1, int t2);
static int depvar_zero (int t1, int t2, int yno, int nwt, 
			const double **Z);
static int lagdepvar (const int *list, const double **Z, const DATAINFO *pdinfo); 
static int jackknife_vcv (MODEL *pmod, const double **Z);
/* end private protos */


/* compute statistics for the dependent variable in a model */

static void model_depvar_stats (MODEL *pmod, const double **Z)
{
    double xx, sum = 0.0;
    int yno = pmod->list[1];
    int t, dwt = 0;

    if (pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy")) {
	dwt = pmod->nwt;
    }

    pmod->ybar = pmod->sdy = NADBL;

    if (pmod->nobs <= 0) {
	return;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}
	if (!model_missing(pmod, t)) {
	    sum += Z[yno][t];
	}
    }

    pmod->ybar = sum / pmod->nobs;

    sum = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}	
	if (!model_missing(pmod, t)) {
	    sum += (Z[yno][t] - pmod->ybar); 
	}
    }

    pmod->ybar = pmod->ybar + sum / pmod->nobs;

    if (fabs(pmod->ybar) < YBARZERO) {
	pmod->ybar = 0.0;
    }

    sum = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}
	if (!model_missing(pmod, t)) {
	    xx = Z[yno][t] - pmod->ybar;
	    sum += xx * xx;
	}
    }

    sum = (pmod->nobs > 1)? sum / (pmod->nobs - 1) : 0.0;

    pmod->sdy = (sum >= 0)? sqrt(sum) : NADBL;
}

/* initial setup for structure to hold info on autoregressive
   coefficients 
*/

static int ar_info_init (MODEL *pmod, int nterms)
{
    int i;

    pmod->arinfo = malloc(sizeof *pmod->arinfo);
    if (pmod->arinfo == NULL) {
	return 1;
    }

    pmod->arinfo->arlist = gretl_list_new(nterms);
    if (pmod->arinfo->arlist == NULL) {
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    pmod->arinfo->rho = malloc(nterms * sizeof *pmod->arinfo->rho);
    if (pmod->arinfo->rho == NULL) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    pmod->arinfo->sderr = malloc(nterms * sizeof *pmod->arinfo->sderr);
    if (pmod->arinfo->sderr == NULL) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo->rho);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    for (i=0; i<nterms; i++) {
	pmod->arinfo->sderr[i] = pmod->arinfo->rho[i] = NADBL;
    }

    return 0;
}

/* determine the degrees of freedom for a model */

static int get_model_df (MODEL *pmod)
{
    int err = 0;

    pmod->ncoeff = pmod->list[0] - 1;

    pmod->dfd = pmod->nobs - pmod->ncoeff;
    if (pmod->dfd < 0) {
	pmod->errcode = E_DF;
        sprintf(gretl_errmsg, _("No. of obs (%d) is less than no. "
		"of parameters (%d)"), pmod->nobs, pmod->ncoeff);
	err = 1;
    } else {
	pmod->dfn = pmod->ncoeff - pmod->ifc;
    }

    return err;
}

#define LDDEBUG 0

static int
transcribe_ld_vcv (MODEL *targ, MODEL *src)
{
    int nv = targ->ncoeff;
    int nxpx = (nv * nv + nv) / 2;
    int i, j;

    if (makevcv(src)) {
	return 1;
    }

    if (targ->vcv == NULL) {
	targ->vcv = malloc(nxpx * sizeof *targ->vcv);
	if (targ->vcv == NULL) {
	    return 1;
	}
    }

    for (i=0; i<nv; i++) {
	for (j=i; j<nv; j++) {
	    targ->vcv[ijton(i, j, nv)] = 
		src->vcv[ijton(i, j, src->ncoeff)];
	}
    }

    return 0;
}

/* Calculate consistent standard errors (and VCV matrix) when doing
   AR(1) estimation of a model with lagged dependent variable.  See
   Ramanathan, Introductory Econometrics, 5e, p. 450.
*/

static int 
ldepvar_std_errors (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    MODEL emod;
    const double *x;

    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;

    double rho = gretl_model_get_double(pmod, "rho_in");
    int origv = pdinfo->v;
    int vnew = pmod->list[0] + 1 - pmod->ifc;

    int *list;
    int vi, vm;
    int i, t;
    int err = 0;

#if LDDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
    printlist(pmod->list, "pmod->list");
    printf("vnew = %d\n", vnew);
    printf("rho = %g\n", rho);
#endif

    list = gretl_list_new(vnew + pmod->ifc);
    if (list == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    err = dataset_add_series(vnew, pZ, pdinfo);
    if (err) {
	free(list);
	pmod->errcode = E_ALLOC;
	return 1;
    }

    vi = origv;

    /* dependent var is residual from original model */
    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[vi][t] = pmod->uhat[t];
    }    
    strcpy(pdinfo->varname[vi], "eps");
    list[1] = vi++;

    /* indep vars are rho-differenced vars from original model */
    for (i=2; i<=pmod->list[0]; i++) {
	vm = pmod->list[i];
	if (vm == 0) {
	    list[i] = 0;
	    continue;
	}
	sprintf(pdinfo->varname[vi], "%.6s_r", pdinfo->varname[vm]);
	x = (*pZ)[vm];
	for (t=0; t<pdinfo->n; t++) {
	    if (t == 0 || na(x[t]) || na(x[t-1])) {
		(*pZ)[vi][t] = NADBL;
	    } else {
		(*pZ)[vi][t] = x[t] - rho * x[t-1];
	    }
	}
	list[i] = vi++;
    }

    /* last indep var is lagged u-hat */
    for (t=0; t<pdinfo->n; t++) {
	if (t == 0) {
	    (*pZ)[vi][t] = NADBL;
	} else { 
	    (*pZ)[vi][t] = (*pZ)[pmod->list[1]][t-1];
	}
	if (na((*pZ)[vi][t])) {
	    continue;
	}
	for (i=0; i<pmod->ncoeff; i++) {
	    x = (*pZ)[pmod->list[i+2]];
	    if (na(x[t-1])) {
		(*pZ)[vi][t] = NADBL;
		break;
	    } else {
		(*pZ)[vi][t] -= pmod->coeff[i] * x[t-1];
	    }
	}
    }

    list[list[0]] = vi;
    strcpy(pdinfo->varname[vi], "uhat_1");

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    emod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (emod.errcode) {
	err = emod.errcode;
    } else {
#if LDDEBUG
	printmodel(&emod, pdinfo, OPT_NONE, prn);
	gretl_print_destroy(prn);
#endif
	for (i=0; i<pmod->ncoeff; i++) {
	    pmod->sderr[i] = emod.sderr[i];
	}

	err = transcribe_ld_vcv(pmod, &emod);
    }
    
    clear_model(&emod);
    
    free(list);
    dataset_drop_last_variables(vnew, pZ, pdinfo);

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    if (err) {
	pmod->errcode = err;
    }

    return err;
}

/* special computation of statistics for autoregressive models */

static int compute_ar_stats (MODEL *pmod, const double **Z, double rho)
{
    int i, t, yno = pmod->list[1];
    double x, pw1 = 0.0;

    if (ar_info_init(pmod, 1)) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    if (pmod->ci == PWE) {
	pw1 = sqrt(1.0 - rho * rho);
    }

    pmod->arinfo->arlist[0] = pmod->arinfo->arlist[1] = 1;
    pmod->arinfo->rho[0] = rho;
    gretl_model_set_double(pmod, "rho_in", rho);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (t == pmod->t1 && pmod->ci == PWE) {
	    x = pw1 * Z[yno][t];
	    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
		x -= pmod->coeff[i] * pw1 * Z[pmod->list[i+2]][t];
	    }
	    if (pmod->ifc) {
		x -= pw1 * pmod->coeff[0];
	    }
	} else {
	    x = Z[yno][t] - rho * Z[yno][t-1];
	    for (i=0; i<pmod->ncoeff; i++) {
		x -= pmod->coeff[i] * 
		    (Z[pmod->list[i+2]][t] - 
		     rho * Z[pmod->list[i+2]][t-1]);
	    }
	}
	pmod->uhat[t] = x;
	pmod->yhat[t] = Z[yno][t] - x;
    }

    pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, Z[yno], pmod->yhat);

    pmod->adjrsq = 
	1.0 - ((1.0 - pmod->rsq) * (pmod->t2 - pmod->t1) / 
	       (double) pmod->dfd);

    return 0;
}

/* calculation of WLS stats in agreement with GNU R */

static void get_wls_stats (MODEL *pmod, const double **Z)
{
    int t, wobs = pmod->nobs, yno = pmod->list[1];
    double x, dy, w2, wmean = 0.0, wsum = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) { 
	if (model_missing(pmod, t)) {
	    continue;
	}
	if (Z[pmod->nwt][t] == 0.0) {
	    wobs--;
	    pmod->dfd -= 1;
	} else {
	    w2 = Z[pmod->nwt][t] * Z[pmod->nwt][t];
	    wmean += w2 * Z[yno][t];
	    wsum += w2;
	}
    }

    wmean /= wsum;
    x = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t) || Z[pmod->nwt][t] == 0.0) {
	    continue;
	}	
	w2 = Z[pmod->nwt][t] * Z[pmod->nwt][t]; 
	dy = Z[yno][t] - wmean;
	x += w2 * dy * dy;
    }

    pmod->fstt = ((x - pmod->ess) * pmod->dfd)/(pmod->dfn * pmod->ess);
    pmod->rsq = (1 - (pmod->ess / x));
    pmod->adjrsq = 1 - ((1 - pmod->rsq) * (pmod->nobs - 1)/pmod->dfd);
}

static void fix_wls_values (MODEL *pmod, double **Z)
{
    int t;

    if (gretl_model_get_int(pmod, "wt_dummy")) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    }
	}
    } else {
	pmod->ess_wt = pmod->ess;
	pmod->sigma_wt = pmod->sigma;
	pmod->ess = 0.0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
		pmod->nobs -= 1;
	    } else {
		pmod->yhat[t] /= Z[pmod->nwt][t];
		pmod->uhat[t] /= Z[pmod->nwt][t];
		pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    }
	}
	pmod->sigma = sqrt(pmod->ess / pmod->dfd);
    }
}

/* drop the weight var from the list of regressors (WLS) */

static void dropwt (int *list)
{
    int i;

    list[0] -= 1;
    for (i=1; i<=list[0]; i++) {
	list[i] = list[i+1];
    }
}

static void model_stats_init (MODEL *pmod)
{
    pmod->ess = pmod->ess_wt = NADBL;
    pmod->sigma = pmod->sigma_wt = NADBL;
    pmod->fstt = pmod->lnL = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
}

#undef SMPL_DEBUG 

static int 
lsq_check_for_missing_obs (MODEL *pmod, gretlopt opts,
			   DATAINFO *pdinfo, const double **Z, 
			   int *misst)
{
    int missv = 0;
    int reject_missing = 0;

    if (reference_missmask_present()) {
	int err = apply_reference_missmask(pmod);

#if SMPL_DEBUG
	fprintf(stderr, "missmask found, applied with err = %d\n",
		err);
#endif
	/* If there was a reference mask present, it was put there
	   as part of a hypothesis test on some original model, and
	   it has to be respected in estimation of this model */

	if (err) {
	    pmod->errcode = E_ALLOC;
	    return 1;
	} else {
	    return 0;
	}
    }

    /* can't do HAC VCV with missing obs in middle */
    if ((opts & OPT_R) && dataset_is_time_series(pdinfo) &&
	!get_force_hc()) {
	reject_missing = 1;
    } 

    if (opts & OPT_M) {
	reject_missing = 1;
    }

    if (reject_missing) {
	/* reject missing obs within adjusted sample */
	missv = adjust_t1t2(pmod, pmod->list, &pmod->t1, &pmod->t2,
			    pdinfo->n, Z, misst);
    } else {
	/* we'll try to compensate for missing obs */
	missv = adjust_t1t2(pmod, pmod->list, &pmod->t1, &pmod->t2,
			    pdinfo->n, Z, NULL);
	if (pmod->ci == POOLED && pmod->missmask != NULL) {
	    if (!model_mask_leaves_balanced_panel(pmod, pdinfo)) {
		gretl_model_set_int(pmod, "unbalanced", 1);
	    } 
	}
    }

#if SMPL_DEBUG
    if (1) {
	char t1s[OBSLEN], t2s[OBSLEN];
	int misscount = model_missval_count(pmod);

	ntodate(t1s, pmod->t1, pdinfo);
	ntodate(t2s, pmod->t2, pdinfo);
	fprintf(stderr, "*** after adjustment, t1=%d (%s), t2=%d (%s)\n", 
		pmod->t1, t1s, pmod->t2, t2s);
	fprintf(stderr, "Valid observations in range = %d\n", 
		pmod->t2 - pmod->t1 + 1 - misscount);
    }
#endif

    return missv;
}

static int const_pos (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
        if (list[i] == 0) {
	    return i;
	}
    }

    return 0;
}

static void 
lagged_depvar_check (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    int ldv = lagdepvar(pmod->list, Z, pdinfo);

    if (ldv) {
	gretl_model_set_int(pmod, "ldepvar", ldv);
    } else if (gretl_model_get_int(pmod, "ldepvar")) {
	gretl_model_set_int(pmod, "ldepvar", 0);
    }
}

#define COLL_DEBUG 0

int redundant_var (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, int trim)
{
    MODEL cmod;
    int targ, ml0 = pmod->list[0];
    int *list;
    int err = E_SINGULAR;
    int i, ret = 0;

    if (ml0 < 3) {
	/* shouldn't happen */
	return 0;
    }

    /* can't handle compound lists */
    for (i=1; i<=ml0; i++) {
	if (pmod->list[i] == LISTSEP) {
	    return 0;
	}
    }

    list = malloc(ml0 * sizeof *list);
    if (list == NULL) {
	return 0;
    }

#if COLL_DEBUG
    fprintf(stderr, "\n*** redundant_var called (trim = %d) ***\n", trim);
    printlist(pmod->list, "original model list");
#endif

    while (err == E_SINGULAR && ml0 > 3) {

	list[0] = ml0 - 1;

	for (targ=ml0; targ>2; targ--) {
	    double ess = 1.0, rsq = 0.0;
	    int j = 2;

	    list[1] = pmod->list[targ];

	    for (i=2; i<=ml0; i++) {
		if (i != targ) {
		    list[j++] = pmod->list[i];
		}
	    }

#if COLL_DEBUG
	    fprintf(stderr, "target list position = %d\n", targ);
	    printlist(list, "temp list for redundancy check");
#endif

	    cmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_Z, 0.0);

	    err = cmod.errcode;
	    if (err == 0) {
		ess = cmod.ess;
		rsq = cmod.rsq;
	    } 

	    clear_model(&cmod);

	    if (err && err != E_SINGULAR) {
		break;
	    } else if (ess == 0.0 || rsq == 1.0) {
		ret = 1;
		break;
	    }
	}

	if (ret) break;

	ml0--;
    }

    if (ret == 1) {
	static char msg[ERRLEN];
	int v = pmod->list[targ];

	/* remove var from list and reduce number of coeffs */
	gretl_list_delete_at_pos(pmod->list, targ);
	pmod->ncoeff -= 1;

	/* compose a message */
	if (trim == 0) {
	    strcpy(msg, _("Omitted due to exact collinearity:"));
	}
	if (pdinfo->varname[v][0] != 0) {
	    strcat(msg, " ");
	    strcat(msg, pdinfo->varname[v]);
	}

	strcpy(gretl_msg, msg);

	/* if there's a lagged dep var, it may have moved */
	if (gretl_model_get_int(pmod, "ldepvar")) {
	    lagged_depvar_check(pmod, (const double **) *pZ, pdinfo);
	}
    }

    free(list);

    return ret;
}

/**
 * lsq:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: one of the command indices in #LSQ_MODEL.
 * @opt: option flags: zero or more of the following --
 *   %OPT_R compute robust standard errors;
 *   %OPT_C force use of Cholesky decomp;
 *   %OPT_A treat as auxiliary regression (don't bother checking
 *     for presence of lagged dependent var, don't augment model count);
 *   %OPT_P use Prais-Winsten for first obs.
 *   %OPT_N don't use degrees of freedom correction for standard
 *      error of regression
 *   %OPT_M reject missing observations within sample range
 *   %OPT_Z (internal use) suppress the automatic elimination of 
 *      perfectly collinear variables
 * @rho: coefficient for rho-differencing the data (0.0 for no
 * differencing)
 *
 * Computes least squares estimates of the model specified by @list,
 * using an estimator determined by the value of @ci.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	   GretlCmdIndex ci, gretlopt opt, double rho)
{
    int l0, yno, i;
    int effobs = 0;
    int missv = 0, misst = 0;
    int jackknife = 0;
    int use_qr = get_use_qr();
    int pwe = (ci == PWE || (opt & OPT_P));
    double *xpy;
    MODEL mdl;

    *gretl_errmsg = '\0';

    if (list == NULL || pZ == NULL || pdinfo == NULL) {
	fprintf(stderr, "E_DATA: lsq: list = %p, pZ = %p, pdinfo = %p\n",
		(void *) list, (void *) pZ, (void *) pdinfo);
	mdl.errcode = E_DATA;
        return mdl;
    }

    if (ci == HSK) {
	return hsk_func(list, pZ, pdinfo);
    } 

    gretl_model_init(&mdl);
    gretl_model_smpl_init(&mdl, pdinfo);
    model_stats_init(&mdl);

    if (pwe) {
	gretl_model_set_int(&mdl, "pwe", 1);
    }

    if (list[0] == 1 || pdinfo->v == 1) {
	fprintf(stderr, "E_DATA: lsq: list[0] = %d, pdinfo->v = %d\n",
		list[0], pdinfo->v);
	mdl.errcode = E_DATA;
        return mdl;
    }

    /* preserve a copy of the list supplied, for future reference */
    mdl.list = gretl_list_copy(list);
    if (mdl.list == NULL) {
        mdl.errcode = E_ALLOC;
        return mdl;
    }

    mdl.t1 = pdinfo->t1;
    mdl.t2 = pdinfo->t2;
    mdl.ci = ci;

    /* Doing weighted least squares? */
    if (ci == WLS) { 
	mdl.nwt = mdl.list[1];
	if (gretl_iszero(mdl.t1, mdl.t2, (*pZ)[mdl.nwt])) {
	    mdl.errcode = E_WTZERO;
	    return mdl;
	}
	effobs = gretl_isdummy(mdl.t1, mdl.t2, (*pZ)[mdl.nwt]);
	if (effobs) {
	    /* the weight var is a dummy, with effobs 1s */
	    gretl_model_set_int(&mdl, "wt_dummy", 1);
	}
    } else {
	mdl.nwt = 0;
    }

    /* sanity check */
    if (mdl.t1 < 0 || mdl.t2 > pdinfo->n - 1) {
        mdl.errcode = E_NODATA;
        goto lsq_abort;
    }

    /* adjust sample range and check for missing obs: this
       may set the model errcode */
    missv = lsq_check_for_missing_obs(&mdl, opt, pdinfo,
				      (const double **) *pZ,
				      &misst);
    if (mdl.errcode) {
        goto lsq_abort;
    }

    /* react to presence of unhandled missing obs */
    if (missv) {
	if (dated_daily_data(pdinfo)) {
	    if (repack_missing_daily_obs(&mdl, *pZ, pdinfo)) {
		return mdl;
	    }
	} else {
	    sprintf(gretl_errmsg, _("Missing value encountered for "
		    "variable %d, obs %d"), missv, misst);
	    mdl.errcode = E_DATA;
	    return mdl;
	} 
    }

    if (ci == WLS) {
	dropwt(mdl.list);
    }

    yno = mdl.list[1];
    
    /* check for unknown vars in list */
    for (i=1; i<=mdl.list[0]; i++) {
        if (mdl.list[i] > pdinfo->v - 1) {
            mdl.errcode = E_UNKVAR;
            goto lsq_abort;
        }
    }       

    /* check for zero dependent var */
    if (depvar_zero(mdl.t1, mdl.t2, yno, mdl.nwt, (const double **) *pZ)) {  
        mdl.errcode = E_ZERO;
        goto lsq_abort; 
    } 

    /* drop any vars that are all zero and repack the list */
    omitzero(&mdl, (const double **) *pZ, pdinfo);

    /* if regressor list contains a constant, place it first */
    i = const_pos(mdl.list);
    mdl.ifc = (i > 1);
    if (i > 2) {
	rearrange_list(mdl.list);
    }

    /* Check for presence of lagged dependent variable? 
       (Don't bother if this is an auxiliary regression.) */
    if (!(opt & OPT_A)) {
	lagged_depvar_check(&mdl, (const double **) *pZ, pdinfo);
    }

    /* AR1: advance the starting observation by one? */
    if (rho != 0.0 && !pwe) {
	mdl.t1 += 1;
    }

    l0 = mdl.list[0];  /* holds 1 + number of coeffs */
    mdl.ncoeff = l0 - 1; 
    if (effobs) {
	mdl.nobs = effobs; /* FIXME? */
    } else {
	mdl.nobs = mdl.t2 - mdl.t1 + 1;
	if (has_missing_obs(&mdl)) {
	    mdl.nobs -= model_missval_count(&mdl);
	}
    }

    /* check degrees of freedom */
    if (get_model_df(&mdl)) {
        goto lsq_abort; 
    }

    /* if df correction is not wanted, record this fact */
    if (opt & OPT_N) {
	gretl_model_set_int(&mdl, "no-df-corr", 1);
    }

    if (dataset_is_time_series(pdinfo)) {
	opt |= OPT_T;
    }

    if (mdl.ci == HCCM || ((opt & OPT_R) && get_hc_version() == 4)) {
	jackknife = 1;
    }

    if (!jackknife && ((opt & OPT_R) || (use_qr && !(opt & OPT_C)))) { 
	mdl.rho = rho;
	gretl_qr_regress(&mdl, pZ, pdinfo, opt);
    } else {
	int trim = 0;
	int l, nxpx;

    trim_var:

	if (trim) {
	    l0 = mdl.list[0];
	    free(mdl.xpx);
	    free(mdl.coeff);
	    free(mdl.sderr);
	    mdl.errcode = 0;
	}
 
	l = l0 - 1;
	nxpx = l * (l + 1) / 2;

	xpy = malloc((l0 + 1) * sizeof *xpy);
	mdl.xpx = malloc(nxpx * sizeof *mdl.xpx);
	mdl.coeff = malloc(mdl.ncoeff * sizeof *mdl.coeff);
	mdl.sderr = malloc(mdl.ncoeff * sizeof *mdl.sderr);

	if (mdl.yhat == NULL) {
	    mdl.yhat = malloc(pdinfo->n * sizeof *mdl.yhat);
	} 
	if (mdl.uhat == NULL) {
	    mdl.uhat = malloc(pdinfo->n * sizeof *mdl.uhat);
	}

	if (xpy == NULL || mdl.xpx == NULL || mdl.coeff == NULL ||
	    mdl.sderr == NULL || mdl.yhat == NULL || mdl.uhat == NULL) {
	    mdl.errcode = E_ALLOC;
	    return mdl;
	}

	for (i=0; i<=l0; i++) {
	    xpy[i] = 0.0;
	}
	for (i=0; i<nxpx; i++) {
	    mdl.xpx[i] = 0.0;
	}

	/* calculate regression results, Cholesky style */
	form_xpxxpy(mdl.list, mdl.t1, mdl.t2, *pZ, mdl.nwt, rho, 
		    pwe, mdl.xpx, xpy, mdl.missmask);

#ifdef XPX_DEBUG
	for (i=0; i<=l0; i++) {
	    fprintf(stderr, "xpy[%d] = %g\n", i, xpy[i]);
	}
	for (i=0; i<nxpx; i++) {
	    fprintf(stderr, "xpx[%d] = %g\n", i, mdl.xpx[i]);
	}
	fputc('\n', stderr);
#endif

	regress(&mdl, xpy, *pZ, pdinfo->n, rho);
	free(xpy);

	if (mdl.errcode == E_SINGULAR && !(opt & OPT_Z) &&
	    redundant_var(&mdl, pZ, pdinfo, trim++)) {
	    goto trim_var;
	}
    }

    if (mdl.errcode) {
	goto lsq_abort;
    }

    /* get the mean and sd of dep. var. and make available */
    model_depvar_stats(&mdl, (const double **) *pZ);

    /* Doing an autoregressive procedure? */
    if (ci == CORC || ci == HILU || ci == PWE) {
	if (compute_ar_stats(&mdl, (const double **) *pZ, rho)) { 
	    goto lsq_abort;
	}
	if (gretl_model_get_int(&mdl, "ldepvar")) {
	    if (ldepvar_std_errors(&mdl, pZ, pdinfo)) {
		goto lsq_abort;
	    }
	}
	if (ci == HILU && (opt & OPT_B)) {
	    gretl_model_set_int(&mdl, "no-corc", 1);
	}
    }

    /* weighted least squares: fix fitted values, ESS, sigma */
    if (ci == WLS) {
	get_wls_stats(&mdl, (const double **) *pZ);
	fix_wls_values(&mdl, *pZ);
    }

    if ((opt & OPT_T) && mdl.missmask == NULL) {
	mdl.rho = rhohat(0, mdl.t1, mdl.t2, mdl.uhat);
	mdl.dw = dwstat(0, &mdl, *pZ);
    } else {
	mdl.rho = mdl.dw = NADBL;
    }

    /* weird special case: degenerate model */
    if (mdl.ncoeff == 1 && mdl.ifc) {
	mdl.rsq = mdl.adjrsq = 0.0;
	mdl.fstt = NADBL;
    }

    /* Generate model selection statistics */
    gretl_calculate_criteria(mdl.criterion, 
			     (mdl.ci == WLS)? mdl.ess_wt : mdl.ess, 
			     mdl.nobs, mdl.ncoeff);

    /* hccm command or HC3a */
    if (jackknife) {
	mdl.errcode = jackknife_vcv(&mdl, (const double **) *pZ);
    }

 lsq_abort:

    /* If we reshuffled any missing observations, put them
       back in their right places now */
    if (gretl_model_get_int(&mdl, "daily_repack")) {
	undo_daily_repack(&mdl, *pZ, pdinfo);
    }

    if (!(opt & OPT_A)) {
	/* if it's not an auxiliary regression, set an ID number
	   on the model */
	set_model_id(&mdl);
    }

    return mdl;
}

/*
  form_xpxxpy: form the X'X matrix and X'y vector

  - if rho is non-zero, quasi-difference the data first
  - if nwt is non-zero, use that variable as weight
  - if pwe is non-zero (as well as rho) construct the
    first observation as per Prais-Winsten

    Z[v][t] = observation t on variable v
    n = number of obs in data set
    t1, t2 = starting and ending observations
    rho = first order serial correlation coefficent
    nwt = ID number of variable used as weight

    xpx = X'X matrix as a lower triangle
          stacked by columns
    xpy = X'y vector
    xpy[0] = sum of y
    xpy[list[0]] = y'y
*/

static int form_xpxxpy (const int *list, int t1, int t2, 
			double **Z, int nwt, double rho, int pwe,
			double *xpx, double *xpy, const char *mask)
{
    int i, j, t;
    int li, lj, m;
    int l0 = list[0], yno = list[1];
    double x, z1, pw1;
    int qdiff = (rho != 0.0);

    /* Prais-Winsten term */
    if (qdiff && pwe) {
	pw1 = sqrt(1.0 - rho * rho);
    } else {
	pwe = 0;
	pw1 = 0.0;
    }

    xpy[0] = xpy[l0] = 0.0;

    for (t=t1; t<=t2; t++) {
	if (missing_masked(mask, t)) {
	    continue;
	}
	x = Z[yno][t]; 
	if (qdiff) {
	    if (pwe && t == t1) {
		x = pw1 * Z[yno][t];
	    } else {
		x -= rho * Z[yno][t-1];
	    }
	} else if (nwt) {
	    x *= Z[nwt][t];
	}
        xpy[0] += x;
        xpy[l0] += x * x;
    }

    if (xpy[l0] <= 0.0) {
         return yno; 
    }    

    m = 0;

    if (qdiff) {
	/* quasi-difference the data */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (pwe && t == t1) {
			x += pw1 * Z[li][t1] * pw1 * Z[lj][t];
		    } else {
			x += (Z[li][t] - rho * Z[li][t-1]) * 
			    (Z[lj][t] - rho * Z[lj][t-1]);
		    }
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		if (pwe && t == t1) {
		    x += pw1 * Z[yno][t] * pw1 * Z[li][t];
		} else {
		    x += (Z[yno][t] - rho * Z[yno][t-1]) *
			(Z[li][t] - rho * Z[li][t-1]);
		}
	    }
	    xpy[i-1] = x;
	}
    } else if (nwt) {
	/* weight the data */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!missing_masked(mask, t)) {
			z1 = Z[nwt][t];
			x += z1 * z1 * Z[li][t] * Z[lj][t];
		    }
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}   
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		if (!missing_masked(mask, t)) {
		    z1 = Z[nwt][t];
		    x += z1 * z1 * Z[yno][t] * Z[li][t];
		}
	    }
	    xpy[i-1] = x;
	}
    } else {
	/* no quasi-differencing or weighting wanted */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!missing_masked(mask, t)) {
			x += Z[li][t] * Z[lj][t];
		    }
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		if (!missing_masked(mask, t)) {
		    x += Z[yno][t] * Z[li][t];
		}
	    }
	    xpy[i-1] = x;
	}
    }

    return 0; 
}

static int make_ess (MODEL *pmod, double **Z)
{
    int i, t, yno = pmod->list[1], l0 = pmod->list[0];
    int nwt = pmod->nwt;
    double yhat, resid;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (nwt && Z[nwt][t] == 0.0) {
	    continue;
	}
	if (model_missing(pmod, t)) {
	    continue;
	}
	yhat = 0.0;
	for (i=2; i<=l0; i++) {
	    yhat += pmod->coeff[i-2] * Z[pmod->list[i]][t];
	}
	resid = Z[yno][t] - yhat;
	if (nwt) {
	    resid *= Z[nwt][t];
	}
	pmod->ess += resid * resid;
    }

    return 0;
}

static void compute_r_squared (MODEL *pmod, double *y)
{
    pmod->rsq = 1.0 - (pmod->ess / pmod->tss);

    if (pmod->dfd > 0) {
	double den = pmod->tss * pmod->dfd;

	if (pmod->ifc) {
	    pmod->adjrsq = 1 - (pmod->ess * (pmod->nobs - 1) / den);
	} else {
	    pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, y, pmod->yhat);
	    pmod->adjrsq = 
		1.0 - ((1.0 - pmod->rsq) * (pmod->nobs - 1.0) / pmod->dfd);
	} 
    }

    if (pmod->rsq < 0.0) {
	pmod->rsq = 0.0;
    }
}

/*
  regress: takes xpx, the X'X matrix produced by form_xpxxpy(), and
  xpy (X'y), and computes ols estimates and associated statistics.

  n = number of observations per series in data set
  ifc = 1 if constant is present, else = 0

  ess = error sum of squares
  sigma = standard error of regression
  fstt = F-statistic
  coeff = vector of regression coefficients
  sderr = vector of standard errors of regression coefficients
*/

static void regress (MODEL *pmod, double *xpy, double **Z, 
		     int n, double rho)
{
    int v, yno = pmod->list[1];
    double ysum, ypy, zz, rss = 0.0;
    double sgmasq = 0.0;
    double *diag = NULL;
    int i, err = 0;

    for (i=0; i<n; i++) {
	pmod->yhat[i] = pmod->yhat[i] = NADBL;
    }    

    ysum = xpy[0];
    ypy = xpy[pmod->ncoeff + 1];

#ifdef NO_LHS_ZERO
    if (floateq(ypy, 0.0)) { 
        pmod->errcode = E_YPY;
        return; 
    }
#endif

    zz = ysum * ysum / pmod->nobs;
    pmod->tss = ypy - zz;

    /*  Cholesky-decompose X'X and find the coefficients */
    err = cholbeta(pmod->xpx, xpy, pmod->coeff, &rss, pmod->ncoeff);
    if (err) {
        pmod->errcode = err;
        return;
    }   
    
    if (rho != 0.0) {
	pmod->ess = ypy - rss;
    } else {
	make_ess(pmod, Z);
	rss = ypy - pmod->ess;
    }

    if (pmod->ess < ESSZERO && pmod->ess > (-ESSZERO)) {
	pmod->ess = 0.0;
    } else if (pmod->ess < 0.0) { 
	sprintf(gretl_errmsg, _("Error sum of squares (%g) is not > 0"),
		pmod->ess);
        return; 
    }

    if (pmod->dfd == 0) {
	pmod->sigma = 0.0;
	pmod->adjrsq = NADBL;
    } else {
	if (gretl_model_get_int(pmod, "no-df-corr")) {
	    sgmasq = pmod->ess / pmod->nobs;
	} else {
	    sgmasq = pmod->ess / pmod->dfd;
	}
	pmod->sigma = sqrt(sgmasq);
    }

    if (floatlt(pmod->tss, 0.0) || floateq(pmod->tss, 0.0)) {
	pmod->rsq = pmod->adjrsq = NADBL;
    } 

    hatvar(pmod, n, Z); 
    if (pmod->errcode) return;

    if (pmod->tss > 0.0) {
	compute_r_squared(pmod, Z[yno]);
    }

#if 0
    if (pmod->ifc && pmod->ncoeff == 1) {
        zz = 0.0;
        pmod->dfn = 1;
    }
#endif

    if (sgmasq <= 0.0 || pmod->dfd == 0 || pmod->dfn == 0) {
	pmod->fstt = NADBL;
    } else {
	pmod->fstt = (rss - zz * pmod->ifc) / (sgmasq * pmod->dfn);
	if (pmod->fstt < 0.0) {
	    pmod->fstt = 0.0;
	}
    }

    diag = malloc(pmod->ncoeff * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    diaginv(pmod->xpx, xpy, diag, pmod->ncoeff);

    for (v=0; v<pmod->ncoeff; v++) {
	if (diag[v] >= 0.0) {
	    pmod->sderr[v] = pmod->sigma * sqrt(diag[v]);
	} else {
	    pmod->sderr[v] = 0.0;
	}
    }

    free(diag); 
    
    return;  
}

/*
  cholbeta: does an in-place Choleski decomposition of xpx (lower
  triangular matrix stacked in columns) and solves the normal
  equations for coeff.  

  xpx = X'X on input and Choleski decomposition on output
  xpy = the X'y vector on input and Choleski-transformed t
        vector on output 
  coeff = array of estimated coefficients 
  nv = number of regression coefficients including the constant

  The number of floating-point operations is basically 3.5 * nv^2
  plus (nv^3) / 3.
*/

static int 
cholbeta (double *xpx, double *xpy, double *coeff, double *rss, int nv)
{
    int i, j, k, kk, l, jm1;
    double e, d, d1, d2, test, xx;

    if (xpx[0] <= 0.0) {
	return E_NAN;
    }

    e = 1.0 / sqrt(xpx[0]);
    xpx[0] = e;
    xpy[1] *= e;
    for (i=1; i<nv; i++) {
	xpx[i] *= e;
    }

    kk = nv;

    for (j=2; j<=nv; j++) {
	/* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d1 += xx * xpy[l];
            d += xx * xx;
            k += nv-l;
        }
        d2 = xpx[kk] - d;
	test = d2 / xpx[kk];
        if (test < TINY) {
#if 0
	    fprintf(stderr, "cholbeta: test = %g\n", test);
#endif
	    if (rss != NULL) *rss = -1.0;
	    return E_SINGULAR;
        }
	if (test < SMALL) {
	    strcpy(gretl_msg, _("Warning: data matrix close to singularity!"));
	}
        e = 1 / sqrt(d2);
        xpx[kk] = e;
        xpy[j] = (xpy[j] - d1) * e;
        for (i=j+1; i<=nv; i++) {
	    /* off-diagonal elements */
            kk++;
            d = 0.0;
            k = j - 1;
            for (l=1; l<=jm1; l++) {
                d += xpx[k] * xpx[k-j+i];
                k += nv - l;
            }
            xpx[kk] = (xpx[kk] - d) * e;
        }
        kk++;
    }

    kk--;

    /* find regression sum of squares */
    if (rss != NULL) {
	d = 0.0;
	for (j=1; j<=nv; j++) {
	    d += xpy[j] * xpy[j];
	}
	*rss = d;
    }

    /* solve for the regression coefficients */
    if (coeff != NULL) {
	for (j=0; j<nv-1; j++) {
	    coeff[j] = 0.0;
	}
	coeff[nv-1] = xpy[nv] * xpx[kk];
	for (j=nv-1; j>=1; j--) {
	    d = xpy[j];
	    for (i=nv-1; i>=j; i--) {
		kk--;
		d -= coeff[i] * xpx[kk];
	    }
	    kk--;
	    coeff[j-1] = d * xpx[kk];
	}
	for (j=0; j<nv; j++) {
	    if (isnan(coeff[j])) {
		return E_NAN;
	    }
	}	
    } 

    return 0; 
}

/*
  diaginv: solves for the diagonal elements of the X'X inverse matrix.

  xpx = Cholesky-decomposed X'X matrix (input)
  xpy = X'y vector (input) used as work array
  diag = diagonal elements of X'X (output)
  nv = number of regression coefficients
*/

static void diaginv (double *xpx, double *xpy, double *diag, int nv)
{
    int kk, l, m, k, i, j;
    const int nxpx = nv * (nv + 1) / 2;
    double d, e;

    kk = 0;

    for (l=1; l<=nv-1; l++) {
        d = xpx[kk];
        xpy[l] = d;
        e = d * d;
        m = 0;
        if (l > 1) {
	    for (j=1; j<=l-1; j++) m += nv - j;
	}
        for (i=l+1; i<=nv; i++) {
            d = 0.0;
            k = i + m - 1;
            for (j=l; j<=i-1; j++) {
                d += xpy[j] * xpx[k];
                k += nv - j;
            }
            d = -d * xpx[k];
            xpy[i] = d;
            e += d * d;
        }
        kk += nv + 1 - l;
        diag[l-1] = e;
    }

    diag[nv-1] = xpx[nxpx-1] * xpx[nxpx-1];
}

/**
 * makevcv:
 * @pmod: pointer to model.
 *
 * Inverts the Choleski-decomposed X'X and computes the 
 * coefficient covariance matrix.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int makevcv (MODEL *pmod)
{
    int dec, mst, kk, i, j, kj, icnt, m, k, l = 0;
    const int nv = pmod->ncoeff;
    const int nxpx = (nv * nv + nv) / 2; 
    double d;

    if (pmod->vcv != NULL) {
	return 0;
    }

    if (pmod->xpx == NULL) {
	fprintf(stderr, "makevcv: pmod->xpx = NULL\n");
	return 1;
    }

    mst = nxpx;
    kk = nxpx - 1;

    pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nv; i++) {
	mst -= i;
	/* find diagonal element */
	d = pmod->xpx[kk];
	if (i > 0) {
	    for (j=kk+1; j<=kk+i; j++) {
		d -= pmod->xpx[j] * pmod->vcv[j];
	    }
	}
	pmod->vcv[kk] = d * pmod->xpx[kk];
	/* find off-diagonal elements indexed by kj */
	kj = kk;
	kk = kk - i - 2;
	if (i > nv - 2) {
	    continue;
	}
	for (j=i+1; j<nv; j++) {
	    icnt = i+1;
	    kj -= j;
	    d = 0.0;
	    m = mst + 1;
	    for (k=0; k<=j-1; k++) {
		if (icnt > 0) {
		    dec = 1;
		    icnt--;
		} else {
		    dec = k;
		}
		m -= dec;
		l = kj + i - k;
		d += pmod->vcv[m-1] * pmod->xpx[l];
	    }
	    pmod->vcv[kj] = -d * pmod->xpx[l-1];
	}
    }

    if (pmod->ci == CUSUM) {
	return 0;
    }

    /* some estimators need special treatment */

    if (pmod->ci != HCCM && pmod->ci != LOGIT && pmod->ci != PROBIT) {
	double sigma = pmod->sigma;

	if ((pmod->ci == WLS && !(gretl_model_get_int(pmod, "wt_dummy"))) || 
	    pmod->ci == ARCH || pmod->ci == HSK) {
	    sigma = pmod->sigma_wt;
	} 
	for (k=0; k<nxpx; k++) {
	    pmod->vcv[k] *= sigma * sigma;
	}
    }

    return 0;
}

/*  dwstat: computes Durbin-Watson statistic
    order is the order of autoregression, 0 for OLS.
*/

static double dwstat (int order, MODEL *pmod, double **Z)
{
    double diff, ut, ut1;
    double diffsq = 0.0;
    int t;

    if (order) order--;

    if (pmod->ess <= 0.0) {
	return NADBL;
    }

    for (t=pmod->t1+1+order; t<=pmod->t2; t++)  {
        ut = pmod->uhat[t];
        ut1 = pmod->uhat[t-1];
        if (na(ut) || na(ut1) ||
	    (pmod->nwt && (floateq(Z[pmod->nwt][t], 0.0) || 
			   floateq(Z[pmod->nwt][t-1], 0.0)))) { 
	    continue;
	}
        diff = ut - ut1;
        diffsq += diff * diff;
    }

    return diffsq / pmod->ess;
}

/* altrho: alternative calculation of rho */

static double altrho (int order, int t1, int t2, const double *uhat)
{
    double *ut, *ut1;    
    int t, n, len = t2 - (t1 + order) + 1;
    double uh, uh1, rho;

    ut = malloc(len * sizeof *ut);
    if (ut == NULL) {
	return NADBL;
    }

    ut1 = malloc(len * sizeof *ut1);
    if (ut1 == NULL) {
	free(ut);
	return NADBL;
    }

    n = 0;

    for (t=t1+order; t<=t2; t++) { 
        uh = uhat[t];
	uh1 = (t > 0)? uhat[t-1] : NADBL;
        if (!na(uh) && !na(uh1)) {
	    ut[n] = uh;
	    ut1[n] = uh1;
	    n++;
	}
    }

    rho = gretl_corr(0, n - 1, ut, ut1, NULL);

    free(ut);
    free(ut1);

    return rho;
}

/*  rhohat: computes first order serial correlation coefficient
    order is the order of autoregression, 0 for OLS.
*/

static double rhohat (int order, int t1, int t2, const double *uhat)
{
    double ut, ut1, uu = 0.0, xx = 0.0;
    double rho;
    int t;

    if (order) {
	order--;
    }

    for (t=t1+order+1; t<=t2; t++) { 
        ut = uhat[t];
        ut1 = uhat[t-1];
        if (na(ut) || na(ut1)) {
	    continue;
	}
        uu += ut * ut1;
        xx += ut1 * ut1;
    }

    if (floateq(xx, 0.0)) {
	return NADBL;
    }

    rho = uu / xx;

    if (rho > 1.0 || rho < -1.0) {
	rho = altrho(order, t1, t2, uhat);
    }

    return rho;
}

/* compute fitted values and residuals */

static int hatvar (MODEL *pmod, int n, double **Z)
{
    int xno, i, t;
    int yno = pmod->list[1];
    double x;

    for (t=0; t<n; t++) {
	pmod->yhat[t] = pmod->uhat[t] = NADBL;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}
	pmod->yhat[t] = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            xno = pmod->list[i+2];
	    x = Z[xno][t];
	    if (pmod->nwt) {
		x *= Z[pmod->nwt][t];
	    }
            pmod->yhat[t] += pmod->coeff[i] * x;
        }
	x = Z[yno][t];
	if (pmod->nwt) {
	    x *= Z[pmod->nwt][t];
	}
        pmod->uhat[t] = x - pmod->yhat[t];                
    }

    return 0;
}

static int hilu_plot (double *ssr, double *rho, int n)
{
    FILE *fp;
    int i;

    if (gnuplot_init(PLOT_REGULAR, &fp)) {
	return E_FOPEN; 
    }

    fputs("# hildreth-lu\n", fp);
    fputs("set xlabel 'rho'\n", fp);

    fprintf(fp, "set ylabel '%s'\n", I_("ESS"));

    fputs("set nokey\n", fp);
    fputs("set xrange [-1.0:1.0]\n", fp);
    fputs("plot '-' using 1:2 w impulses\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    for (i=0; i<n; i++) {
	fprintf(fp, "%g %g\n", rho[i], ssr[i]);
    }
    fputs("e\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    gnuplot_make_graph();

    return 0;
}

#undef USE_DW /* base rho-hat on the D-W statistic? */

static double autores (MODEL *pmod, const double **Z, int ci)
{
    int t, v, t1 = pmod->t1;
    double x, num = 0.0, den = 0.0;
    double rhohat;

    if (ci == CORC || ci == HILU) {
	t1--;
    }

    for (t=t1; t<=pmod->t2; t++) {
	x = Z[pmod->list[1]][t];
	for (v=0; v<pmod->ncoeff; v++) {
	    x -= pmod->coeff[v] * Z[pmod->list[v+2]][t];
	}
	pmod->uhat[t] = x;
	if (t > t1) {
#ifdef USE_DW
	    x = pmod->uhat[t] - pmod->uhat[t-1];
	    num += x * x;
#else
	    num += pmod->uhat[t] * pmod->uhat[t-1];
#endif
	    den += pmod->uhat[t-1] * pmod->uhat[t-1];
	}
    } 

#ifdef USE_DW
    den += pmod->uhat[pmod->t2] * pmod->uhat[pmod->t2];
    rhohat = 1.0 - num / (den * 2.0);
#else
    rhohat = num / den;
#endif

    return rhohat;
}

#define AR_DEBUG 0

/**
 * estimate_rho:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: %CORC for Cochrane-Orcutt, %HILU for Hildreth-Lu,
 *      %PWE for Prais-Winsten estimator.
 * @err: pointer for error code.
 * @opt: option flags: may include %OPT_B to suppress Cochrane-Orcutt
 *       fine-tuning of Hildreth-Lu results, %OPT_P to generate
 *       a gnuplot graph of the search in case @ci = %HILU.
 * @prn: gretl printing struct.
 *
 * Estimate the quasi-differencing coefficient for use with the
 * Cochrane-Orcutt, Hildreth-Lu or Prais-Winsten procedures for
 * handling first-order serial correlation.  Print a trace of the
 * search for rho.
 * 
 * Returns: rho estimate on successful completion, %NADBL on error.
 */

double estimate_rho (const int *list, double ***pZ, DATAINFO *pdinfo,
		     GretlCmdIndex ci, int *err, gretlopt opt, PRN *prn)
{
    double rho = 0.0, rho0 = 0.0, diff;
    double finalrho = 0.0, essmin = 1.0e8;
    double ess, ssr[199], rh[199]; 
    int iter, nn = 0;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int missv = 0, misst = 0;
    gretlopt lsqopt = OPT_A;
    int ascii = !(opt & OPT_P);
    MODEL corc_model;

    *gretl_errmsg = '\0';
    *err = 0;

    missv = adjust_t1t2(NULL, list, &pdinfo->t1, &pdinfo->t2, 
			pdinfo->n, (const double **) *pZ, &misst);
    if (missv) {
	sprintf(gretl_errmsg, _("Missing value encountered for "
				"variable %d, obs %d"), missv, misst);
	*err = E_DATA;
	goto bailout;
    }

    gretl_model_init(&corc_model);

    if (ci == PWE) {
	/* Prais-Winsten treatment of first observation */
	lsqopt |= OPT_P;
    } 

    if (ci == HILU) { 
	/* Do Hildreth-Lu first */
	for (rho = -0.990, iter = 0; rho < 1.0; rho += .01, iter++) {
	    clear_model(&corc_model);
	    corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
	    if ((*err = corc_model.errcode)) {
		clear_model(&corc_model);
		goto bailout;
	    }
	    ess = corc_model.ess;
	    if (ascii) {
		char num[16];
		int chk;
		
		if (iter == 0) {
		    pprintf(prn, "\n RHO       %s      RHO       %s      "
			    "RHO       %s      RHO       %s     \n",
			    _("ESS"), _("ESS"), _("ESS"), _("ESS"));
		}
		sprintf(num, "%f", 100 * fabs(rho));
		chk = atoi(num);
		if (chk == 99 || chk % 10 == 0) {
		    ssr[nn] = ess;
		    rh[nn++] = rho;
		    pprintf(prn, "%5.2f %10.4g", rho, ess);
		    if (nn % 4 == 0) {
			pputc(prn, '\n');
		    } else {
			bufspace(3, prn);
		    }
		} 
	    } else {
		ssr[nn] = ess;
		rh[nn++] = rho;
	    }
	    if (iter == 0 || ess < essmin) {
		essmin = ess;
		finalrho = rho;
	    }
	} /* end of basic iteration */
	
	if (finalrho > 0.989) {
	    /* try exploring this funny region? */
	    for (rho = 0.99; rho <= 0.999; rho += .001) {
		clear_model(&corc_model);
		corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
		if ((*err = corc_model.errcode)) {
		    clear_model(&corc_model);
		    goto bailout;
		}
		ess = corc_model.ess;
		if (ess < essmin) {
		    essmin = ess;
		    finalrho = rho;
		}
	    }
	}

	if (finalrho > 0.9989) {
	    /* this even funnier one? */
	    for (rho = 0.9991; rho <= 0.9999; rho += .0001) {
		clear_model(&corc_model);
		corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
		if ((*err = corc_model.errcode)) {
		    clear_model(&corc_model);
		    goto bailout;
		}
		ess = corc_model.ess;
		if (ess < essmin) {
		    essmin = ess;
		    finalrho = rho;
		}
	    }
	}

	rho0 = rho = finalrho;
	pprintf(prn, _("\n\nESS is minimum for rho = %g\n\n"), rho);
	if (ascii) {
	    graphyzx(NULL, ssr, NULL, rh, nn, "ESS", "RHO", NULL, 0, prn); 
	    pputs(prn, "\n");
	} else {
	    hilu_plot(ssr, rh, nn);
	}
    } else { 
	/* Go straight to Cochrane-Orcutt (or Prais-Winsten) */
	corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (!corc_model.errcode && corc_model.dfd == 0) {
	    corc_model.errcode = E_DF;
	}
	if ((*err = corc_model.errcode)) {
	    clear_model(&corc_model);
	    goto bailout;
	}
#ifdef USE_DW
	rho0 = rho = 1.0 - (corc_model.dw / 2.0);
#else
	rho0 = rho = corc_model.rho;
#endif
    }

    if (na(rho)) {
	*err = E_NOCONV;
	clear_model(&corc_model);
	goto bailout;
    }

    if (ci != HILU || !(opt & OPT_B)) {

	if (ci == HILU) {
	    pputs(prn, _("\nFine-tune rho using the CORC procedure...\n\n"));
	} else {
	    pputs(prn, _("\nPerforming iterative calculation of rho...\n\n"));
	}

	pputs(prn, _("                 ITER       RHO        ESS"));
	pputc(prn, '\n');

	iter = 0;
	diff = 1.0;

	while (diff > 0.001) {
	    pprintf(prn, "          %10d %12.5f", ++iter, rho);
	    clear_model(&corc_model);
	    corc_model = lsq(list, pZ, pdinfo, OLS, lsqopt, rho);
#if AR_DEBUG
	    fprintf(stderr, "corc_model: t1=%d, first two uhats: %g, %g\n",
		    corc_model.t1, 
		    corc_model.uhat[corc_model.t1],
		    corc_model.uhat[corc_model.t1+1]);
#endif
	    if ((*err = corc_model.errcode)) {
		clear_model(&corc_model);
		goto bailout;
	    }
	    pprintf(prn, "   %g\n", corc_model.ess);

	    rho = autores(&corc_model, (const double **) *pZ, ci);

#if AR_DEBUG
	    pputs(prn, "CORC model (using rho-transformed data)\n");
	    printmodel(&corc_model, pdinfo, OPT_NONE, prn);
	    pprintf(prn, "autores gives rho = %g\n", rho);
#endif

	    if (rho > .99999 || rho < -.99999) {
		*err = E_NOCONV;
		clear_model(&corc_model);
		goto bailout;
	    }

	    diff = (rho > rho0) ? rho - rho0 : rho0 - rho;
	    rho0 = rho;
	    if (iter == 30) break;
	}

	pprintf(prn, _("                final %11.5f\n\n"), rho);
    }

    clear_model(&corc_model);

 bailout:

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    if (*err) {
	rho = NADBL;
    }

    return rho;
}

/* .......................................................... */

static int get_pos (const int *list)
{
    int i, ret = -1;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    ret = i;
	    break;
	}
    }

    return ret;
}

void tsls_free_data (const MODEL *pmod)
{
    const char *endog = gretl_model_get_data(pmod, "endog");
    double **X = gretl_model_get_data(pmod, "tslsX");
    int i, m = 0;

    if (endog != NULL && X != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    if (endog[i]) m++;
	}
	for (i=0; i<m; i++) {
	    free(X[i]);
	}
    }
}

static int tsls_save_data (MODEL *pmod, const int *list, 
			   double **Z, DATAINFO *pdinfo)
{
    double **X = NULL;
    char *endog = NULL;
    int addvars = list[0];
    int i, j, k, pos, m, err = 0;
    size_t esize, Xsize = list[0] * sizeof *X;
    size_t es_old, xs_old;
    int recycle_X = 0;
    int recycle_e = 0;

    pos = get_pos(pmod->list);
    m = pos - 2;
    esize = m * sizeof *endog;

    /* re-use old pointers if applicable */

    X = gretl_model_get_data_and_size(pmod, "tslsX", &xs_old);
    if (X == NULL) {
	X = malloc(Xsize);
    } else if (Xsize != xs_old) {
	tsls_free_data(pmod);
	gretl_model_detach_data_item(pmod, "tslsX");
	free(X);
	X = malloc(Xsize);
    } else {
	recycle_X = 1;
    }

    endog = gretl_model_get_data_and_size(pmod, "endog", &es_old);
    if (endog == NULL) {
	endog = malloc(esize);
    } else if (esize != es_old) {
	gretl_model_detach_data_item(pmod, "endog");
	free(endog);
	endog = malloc(esize);
    } else {
	recycle_e = 1;
    }

    if (X == NULL || endog == NULL) {
	free(X);
	free(endog);
	return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
	k = pdinfo->v - 1 + i - addvars;
	X[i-1] = Z[k];
	Z[k] = NULL;
    }

    for (i=0; i<m; i++) {
	k = pmod->list[i+2];
	endog[i] = 0;
	for (j=1; j<=list[0]; j++) {
	    if (list[j] == k) {
		endog[i] = 1;
		break;
	    }
	}
    }

    /* now attach X and endog to the model */
    if (!recycle_X) {
	gretl_model_set_data(pmod, "tslsX", X, Xsize);
    }
    if (!recycle_e) {
	gretl_model_set_data(pmod, "endog", endog, esize);
    }

    return err;
}

double *tsls_get_Xi (const MODEL *pmod, double **Z, int i)
{
    const char *endog;
    double **X;
    double *ret;

    endog = gretl_model_get_data(pmod, "endog");
    X = gretl_model_get_data(pmod, "tslsX");

    if (endog == NULL || X == NULL) return NULL;

    if (!endog[i]) {
	ret = Z[pmod->list[i+2]];
    } else {
	int j, k = 0;

	for (j=0; j<i; j++) {
	    if (endog[j]) k++;
	}
	ret = X[k];
    }

    return ret;
}

/*
  tsls_make_replist: determines which variables in reglist, when
  checked against the predetermined and exogenous vars in instlist,
  need to be instrumented, and populates replist accordingly
*/

static int 
tsls_make_replist (const int *reglist, int *instlist, int *replist)
{
    int i, j, k = 0;
    int endog, fixup = 0;

    for (i=2; i<=reglist[0]; i++) {  
	endog = 1;
	for (j=1; j<=instlist[0]; j++) {
	    if (instlist[j] == reglist[i]) {
		endog = 0;
		break;
	    }
	}
	if (reglist[i] == 0 && endog) {
	    /* const is in reglist but not instlist: needs fixing */
	    fixup = 1;
	} else if (endog) {
	    replist[++k] = reglist[i];
	} 
    }

    replist[0] = k;

    if (fixup) {
	/* add const to list of instruments */
	instlist[0] += 1;
	for (i=instlist[0]; i>=2; i--) {
	    instlist[i] = instlist[i - 1];
	}
	instlist[1] = 0;
    }

    return 0;
}

/**
 * tsls_func:
 * @list: dependent variable plus list of regressors.
 * @pos_in: position in the list for the separator between list
 *   of variables and list of instruments (or 0 if unknown).
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain %OPT_R for robust VCV, %OPT_S to save second-
 * stage regressors (%OPT_S is used in context of three-stage least 
 * squares), %OPT_N (no df correction).
 *
 * Estimate the model given in @list by means of Two-Stage Least
 * Squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tsls_func (const int *list, int pos_in, double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt)
{
    int i, t;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int *reglist = NULL, *instlist = NULL, *replist = NULL;
    int *s1list = NULL, *s2list = NULL;
    int pos, nelem, orig_nvar = pdinfo->v;
    MODEL tsls;

    if (pos_in > 0) {
	pos = pos_in;
    } else {
	pos = get_pos(list);
    }

    /* initialize model (in case we bail out early on error) */
    gretl_model_init(&tsls);

    *gretl_errmsg = '\0';

    /* adjust sample range for missing observations */
    varlist_adjust_sample(list, &pdinfo->t1, &pdinfo->t2, 
			  (const double **) *pZ); 

    /* 
       reglist: dependent var plus list of regressors
       instlist: list of instruments
       s1list: regression list for first-stage regressions
       s2list: regression list for second-stage regression
       replist: list of vars to be replaced by fitted values
                for the second-stage regression
    */
    reglist = malloc(pos * sizeof *reglist);
    instlist = malloc((list[0] - pos + 2) * sizeof *instlist);
    s1list = malloc((list[0] - pos + 3) * sizeof *s1list);
    s2list = malloc(pos * sizeof *s2list);
    replist = malloc(pos * sizeof *replist);

    if (reglist == NULL || instlist == NULL || s1list == NULL ||
	s2list == NULL || replist == NULL) {
	tsls.errcode = E_ALLOC;
	goto tsls_bailout;
    }

    /* dep. var. plus regressors: first portion of input list */
    reglist[0] = pos - 1;
    for (i=1; i<pos; i++) {
	reglist[i] = list[i];
    }    

    /* set up list of instruments: second portion of input list */
    instlist[0] = list[0] - pos;
    for (i=1; i<=instlist[0]; i++) {
	instlist[i] = list[i + pos];
    }	

    /* drop any vars that are all zero, and reshuffle the constant
       into first position among the independent vars 
    */
    tsls_omitzero(reglist, (const double **) *pZ, pdinfo->t1, pdinfo->t2);
    rearrange_list(reglist);
    tsls_omitzero(instlist, (const double **) *pZ, pdinfo->t1, pdinfo->t2);

    /* initial composition of second-stage regression list (will be
       modified below) 
    */
    for (i=0; i<pos; i++) {
	s2list[i] = reglist[i];
    }

    /* determine the list of variables (replist) for which we need to
       obtain fitted values in the first stage 
    */
    tsls_make_replist(reglist, instlist, replist);

    /* check for order condition */
    if (instlist[0] < reglist[0] - 1) {
        sprintf(gretl_errmsg, 
		_("Order condition for identification is not satisfied.\n"
		"varlist 2 needs at least %d more variable(s) not in "
		"varlist1."), reglist[0] - 1 - instlist[0]);
	tsls.errcode = E_UNSPEC; 
	goto tsls_bailout;
    }

    /* replist[0] holds the number of fitted vars to create */
    if (dataset_add_series(replist[0], pZ, pdinfo)) {
	tsls.errcode = E_ALLOC;
	goto tsls_bailout;
    } 

    /* common setup for first-stage regressions: regressors are copied
       from instlist */
    s1list[0] = instlist[0] + 1;
    for (i=2; i<=s1list[0]; i++) {
	s1list[i] = instlist[i-1];
    }    

    /* 
       deal with the variables for which instruments are needed: cycle
       through the list of variables to be instrumented (replist), run
       the first stage regression, and add the fitted values into the
       data matrix Z
    */
    for (i=1; i<=replist[0]; i++) { 
	int newv = orig_nvar + i - 1;
	int j;

	/* select the dependent variable */
        s1list[1] = replist[i];

	/* run the first-stage regression */
	tsls = lsq(s1list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (tsls.errcode) {
	    goto tsls_bailout;
	}

	/* write the fitted values into data matrix, Z */
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[newv][t] = tsls.yhat[t];
	}

	clear_model(&tsls);

	/* give the fitted series the same name as the original */
	strcpy(pdinfo->varname[newv], pdinfo->varname[replist[i]]);

	/* substitute the newly created var into the right place in the
	   second-stage regression list */
	for (j=2; j<=s2list[0]; j++) {
	    if (s2list[j] == replist[i]) {
		s2list[j] = newv;
		break;
	    }
	}
    } 

    /* second-stage regression */
    tsls = lsq(s2list, pZ, pdinfo, OLS, OPT_NONE, 0.0);
    if (tsls.errcode) {
	goto tsls_bailout;
    }

    /* special: we need to use the original RHS vars to compute
       residuals and associated statistics */

    tsls.ess = 0.0;
    for (t=tsls.t1; t<=tsls.t2; t++) {
	double yh = 0.0;

	if (model_missing(&tsls, t)) {
	    continue;
	}
	for (i=0; i<tsls.ncoeff; i++) {
	    yh += tsls.coeff[i] * (*pZ)[reglist[i+2]][t];
	}
	tsls.yhat[t] = yh; 
	tsls.uhat[t] = (*pZ)[tsls.list[1]][t] - yh;
	tsls.ess += tsls.uhat[t] * tsls.uhat[t];
    }

    if (tsls.ess <= 0.0) {
	tsls.sigma = 0.0;
    } else {
	tsls.sigma = sqrt(tsls.ess / ((opt & OPT_N)? tsls.nobs : tsls.dfd));
    }

    /* computation of covariance matrix of parameter estimates */

    if ((opt & OPT_R) || get_use_qr()) {
	/* QR decomp in force, or robust standard errors called for */
	qr_tsls_vcv(&tsls, (const double **) *pZ, opt);
    } else {
	double *xpx = NULL, *xpy = NULL, *diag = NULL;
	int nxpx = tsls.ncoeff * (tsls.ncoeff + 1) / 2;

	xpx = malloc(nxpx * sizeof *xpx);
	xpy = malloc((s2list[0] + 1) * sizeof *xpy);
	diag = malloc(tsls.ncoeff * sizeof *diag);

	if (xpy == NULL || xpx == NULL || diag == NULL) {
	    free(xpx);
	    free(xpy);
	    free(diag);
	    tsls.errcode = E_ALLOC;
	    goto tsls_bailout;
	}

	form_xpxxpy(s2list, tsls.t1, tsls.t2, *pZ, 0, 0.0, 0,
		    xpx, xpy, tsls.missmask);

	cholbeta(xpx, xpy, NULL, NULL, tsls.ncoeff);    
	diaginv(xpx, xpy, diag, tsls.ncoeff);

	for (i=0; i<tsls.ncoeff; i++) {
	    tsls.sderr[i] = tsls.sigma * sqrt(diag[i]); 
	}

	free(diag); 
	free(xpx);
	free(xpy);
    }

    /* computation of additional statistics (R^2, F, etc.) */

    tsls.rsq = gretl_corr_rsq(tsls.t1, tsls.t2, (*pZ)[tsls.list[1]],
			      tsls.yhat);
    tsls.adjrsq = 
	1.0 - ((1.0 - tsls.rsq) * (tsls.nobs - 1.0) / tsls.dfd);
    tsls.fstt = tsls.rsq * tsls.dfd / (tsls.dfn * (1.0 - tsls.rsq));
    ls_aic_bic(&tsls);

    if (tsls.missmask == NULL) {
	/* no missing obs within sample range */
	tsls.rho = rhohat(0, tsls.t1, tsls.t2, tsls.uhat);
	tsls.dw = dwstat(0, &tsls, *pZ);
    } else {
	tsls.rho = tsls.dw = NADBL;
    }

    /* set command code on the model */
    tsls.ci = TSLS;

    /* write the full tsls list (dep. var. and regressors, followed by
       instruments, possibly purged of zero elements) into the model
       for future reference
    */
    nelem = reglist[0] + instlist[0];
    tsls.list = realloc(tsls.list, (nelem + 2) * sizeof *tsls.list);
    if (tsls.list == NULL) {
	tsls.errcode = E_ALLOC;
    } else {
	int j = 1;

	tsls.list[0] = nelem + 1;
	for (i=1; i<=reglist[0]; i++) {
	    tsls.list[j++] = reglist[i];
	}
	tsls.list[j++] = LISTSEP;
	for (i=1; i<=instlist[0]; i++) {
	    tsls.list[j++] = instlist[i];
	}
    }

 tsls_bailout:

    free(reglist); 
    free(instlist);
    free(s1list); 
    free(s2list);

    /* save first-stage fitted values, if wanted */
    if ((opt & OPT_S) && tsls.errcode == 0) {
	tsls_save_data(&tsls, replist, *pZ, pdinfo);
    }

    free(replist);

    /* delete first-stage fitted values from dataset */
    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    if (tsls.errcode) {
	model_count_minus();
    }

    /* restore original sample range */
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    return tsls;
}

/* Given an original regression list, augment it by adding the squares
   (and cross-products, if wanted) or logs of the original independent
   variables.  Generate these variables if need be.  Return the
   augmented list, or NULL on failure.
*/

int *augment_regression_list (const int *orig, int aux, 
			      double ***pZ, DATAINFO *pdinfo)
{
    int *list;
    int listlen;
    int i, k;

    if (aux == AUX_WHITE) {
	int trv = orig[0] - 1 - gretl_list_has_const(orig);
	int nt = (trv * trv + trv) / 2;

	listlen = orig[0] + nt + 1;
    } else {
	listlen = 2 * orig[0];
    }

    list = malloc(listlen * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    /* transcribe original list */
    for (i=0; i<=orig[0]; i++) {
	list[i] = orig[i];
    }

    /* add squares, cross-products or logs of independent vars */
    k = list[0];
    for (i=2; i<=orig[0]; i++) {
	int vnew, vi = orig[i];

	if (vi == 0) {
	    continue;
	}

	if (aux == AUX_SQ || aux == AUX_WHITE) {
	    vnew = xpxgenr(vi, vi, pZ, pdinfo);
	    if (vnew > 0) list[++k] = vnew;
	    if (aux == AUX_WHITE) {
		int j, vj;

		for (j=i+1; j<=orig[0]; j++) {
		    vj = orig[j];
		    if (vj == 0) {
			continue;
		    }
		    vnew = xpxgenr(vi, vj, pZ, pdinfo);
		    if (vnew > 0) list[++k] = vnew;
		}
	    }
	} else if (aux == AUX_LOG) {
	    vnew = loggenr(vi, pZ, pdinfo);
	    if (vnew > 0) list[++k] = vnew;
	}
    }

    list[0] = k;

    return list;
}

/* get_hsk_weights: take the residuals from the model pmod, square them
   and take logs; find the fitted values for this series using an
   auxiliary regression including the original independent variables
   and their squares; transform the fitted values by exponentiating
   and taking the square root; and add the resulting series to the
   data set
*/

static int get_hsk_weights (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    int oldv = pdinfo->v;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int *list = NULL;
    int err = 0, shrink = 0;
    double xx;
    MODEL aux;

    /* allocate space for an additional variable */
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    /* add transformed pmod residuals to data set */
    for (t=0; t<pdinfo->n; t++) {
	if (na(pmod->uhat[t])) {
	    (*pZ)[oldv][t] = NADBL;
	} else {
	    xx = pmod->uhat[t];
	    (*pZ)[oldv][t] = log(xx * xx);
	}
    }

    /* build regression list, adding the squares of the original
       independent vars */
    list = augment_regression_list(pmod->list, AUX_SQ, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    list[1] = oldv; /* the newly added uhat-squared */

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    aux = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.);
    err = aux.errcode;
    if (err) {
	shrink = pdinfo->v - oldv;
    } else {
	/* write into the data set the required weights */
	for (t=aux.t1; t<=aux.t2; t++) {
	    if (na(aux.yhat[t])) {
		(*pZ)[oldv][t] = NADBL;
	    } else {
		xx = aux.yhat[t];
		(*pZ)[oldv][t] = 1.0 / sqrt(exp(xx));
	    }
	}
	shrink = pdinfo->v - oldv - 1;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    clear_model(&aux);

    if (shrink > 0) {
	dataset_drop_last_variables(shrink, pZ, pdinfo);
    }

    free(list);

    return err;
}

/**
 * hsk_func:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using a correction for
 * heteroskedasticity.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL hsk_func (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int i, err;
    int orig_nvar = pdinfo->v;
    int *hsklist;
    MODEL hsk;

    *gretl_errmsg = '\0';

    /* run initial OLS */
    hsk = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (hsk.errcode) {
	return hsk;
    }

    /* use the residuals from the initial OLS to form weights */
    err = get_hsk_weights(&hsk, pZ, pdinfo);
    if (err) {
	hsk.errcode = err;
	return hsk;
    }

    /* allocate regression list for weighted least squares */
    hsklist = gretl_list_new(list[0] + 1);
    if (hsklist == NULL) {
	hsk.errcode = E_ALLOC;
	return hsk;
    }

    /* the last variable in the dataset will be the weight var */
    hsklist[1] = pdinfo->v - 1;

    /* put the original dependent variable in at position 2 */
    hsklist[2] = list[1];

    /* add the original independent vars into the WLS list */
    for (i=3; i<=hsklist[0]; i++) {
	hsklist[i] = list[i-1];
    }

    clear_model(&hsk);
    hsk = lsq(hsklist, pZ, pdinfo, WLS, OPT_NONE, 0.0);
    hsk.ci = HSK;

    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    free(hsklist);

    return hsk;
}

static double **allocate_hccm_p (int k, int n)
{
    double **p = malloc(k * sizeof *p);
    int i;

    if (p == NULL) return NULL;

    for (i=0; i<k; i++) {
	p[i] = malloc(n * sizeof **p);
	if (p[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) {
		free(p[j]);
	    }
	    free(p);
	    p = NULL;
	    break;
	}
    }

    return p;
}

static void free_hccm_p (double **p, int m)
{
    if (p != NULL) {
	int i;

	for (i=0; i<m; i++) {
	    free(p[i]);
	}
	free(p);
    }
}

static int jackknife_vcv (MODEL *pmod, const double **Z)
{
    double *st = NULL, *ustar = NULL;
    double **p = NULL;
    int nobs, tp, nc = 0;
    int i, j, k, t;
    int t1, t2;
    double xx;
    int err = 0;

    *gretl_errmsg = '\0';

    t1 = pmod->t1;
    t2 = pmod->t2;
    nobs = pmod->nobs;
    nc = pmod->ncoeff;

    st = malloc(nc * sizeof *st);
    ustar = malloc(nobs * sizeof *ustar);
    p = allocate_hccm_p(nc, nobs);

    if (st == NULL || p == NULL || ustar == NULL) {
	err = E_ALLOC;
	goto bailout;
    }  

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
	pmod->vcv = NULL;
    }

    pmod->ci = HCCM;

    if (makevcv(pmod)) {
	err = E_ALLOC;
	goto bailout;
    }

    /* form elements of (X'X)^{-1}X' */

    for (i=0; i<nc; i++) {
	tp = 0;
	for (t=t1; t<=t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    xx = 0.0;
	    for (j=0; j<nc; j++) {
		if (i <= j) {
		    k = ijton(i, j, nc);
		} else {
		    k = ijton(j, i, nc);
		}
		xx += pmod->vcv[k] * Z[pmod->list[j+2]][t];
	    }
	    p[i][tp++] = xx;
	}
    }

    tp = 0;
    for (t=t1; t<=t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}	
	xx = 0.0;
	for (i=0; i<nc; i++) {
	    xx += Z[pmod->list[i+2]][t] * p[i][tp];
	}
	if (floateq(xx, 1.0)) {
	    xx = 0.0;
	}
	ustar[tp++] = pmod->uhat[t] / (1.0 - xx);
    }

    for (i=0; i<nc; i++) {
	xx = 0.0;
	for (t=0; t<nobs; t++) {
	    xx += p[i][t] * ustar[t];
	}
	st[i] = xx;
    }

    for (t=0; t<nobs; t++) {
	for (i=0; i<nc; i++) {
	    p[i][t] *= ustar[t];
	}
    }

    /* MacKinnon and White, 1985, equation (13) */

    k = 0;
    for (i=0; i<nc; i++) {
	for (j=i; j<nc; j++) {
	    xx = 0.0;
	    for (t=0; t<nobs; t++) {
		xx += p[i][t] * p[j][t];
	    }
	    xx -= st[i] * st[j] / nobs;
	    /* MacKinnon and White: "It is tempting to omit the factor
	       (n - 1)/n from HC3" (1985, p. 309).  Here we leave it in
	       place, as in their simulations.
	    */
	    xx *= (nobs - 1.0) / nobs;
	    if (i == j) {
		pmod->sderr[i] = sqrt(xx);
	    }
	    pmod->vcv[k++] = xx;
	}
    }

    /* substitute robust F stat */
    if (pmod->dfd > 0 && pmod->dfn > 1) {
	pmod->fstt = robust_omit_F(NULL, pmod);
    }

    gretl_model_set_int(pmod, "robust", 1);
    gretl_model_set_int(pmod, "hc", 1);
    gretl_model_set_int(pmod, "hc_version", 4);

 bailout:

    pmod->ci = OLS;

    free(st);
    free(ustar);
    free_hccm_p(p, nc);

    return err;
}

/**
 * whites_test:
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save results to model.
 * @prn: gretl printing struct.
 *
 * Runs White's test for heteroskedasticity on the given model,
 * putting the results into @test.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int whites_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    int lo, ncoeff, yno, t;
    int v = pdinfo->v;
    int *list = NULL;
    double zz;
    MODEL white;
    int err = 0;

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == LOGISTIC) { 
	return E_NOTIMP;
    }

    if ((err = list_members_replaced(pmod->list, pdinfo, pmod->ID))) {
	return err;
    }

    gretl_model_init(&white);

    lo = pmod->list[0];
    yno = pmod->list[1];
    ncoeff = pmod->list[0] - 1;

    /* make space in data set */
    if (dataset_add_series(1, pZ, pdinfo)) {
	err = E_ALLOC;
    }

    if (!err) {
	/* get residuals, square and add to data set */
	for (t=0; t<pdinfo->n; t++) {
	    zz = pmod->uhat[t];
	    if (na(zz)) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = zz * zz;
	    }
	}
	strcpy(pdinfo->varname[v], "uhatsq");
    }

    if (!err) {
	/* build aux regression list, adding squares and
	   cross-products of the original independent vars */
	list = augment_regression_list(pmod->list, AUX_WHITE, pZ, pdinfo);
	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    list[1] = v; /* the newly added uhat-squared */
	}
    }

    if (!err) {
	/* run auxiliary regression */
	white = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.);
	err = white.errcode;
    }

    if (!err) {
	double TR2, pval;

	white.aux = AUX_WHITE;
	printmodel(&white, pdinfo, OPT_NONE, prn);

	TR2 = white.rsq * white.nobs;
	pval = chisq(TR2, white.ncoeff - 1);

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_WHITES);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_TR2);
		model_test_set_dfn(test, white.ncoeff - 1);
		model_test_set_value(test, TR2);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	  
	}

	record_test_result(TR2, pval, "White's");
    }

    clear_model(&white);

    dataset_drop_last_variables(pdinfo->v - v, pZ, pdinfo);

    free(list);

    return err;
}

static int ar_list_max (const int *list) 
{
    int i, lmax = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] > lmax) {
	    lmax = list[i];
	}
    }

    return lmax;
}

/**
 * ar_func:
 * @list: list of lags plus dependent variable and list of regressors.
 * @pos: position in list of separator (dummy element) between lag list
 * and specification of dependent and independent variables.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list using the generalized 
 * Cochrane-Orcutt procedure for autoregressive errors.
 * 
 * Returns: #MODEL struct containing the results.
 */

MODEL ar_func (const int *list, int pos, double ***pZ, 
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    double diff, ess, tss, xx;
    int i, j, t, t1, t2, vc, yno, ryno, iter;
    int err, lag, maxlag, v = pdinfo->v;
    int *arlist = NULL, *rholist = NULL;
    int *reglist = NULL, *reglist2 = NULL;
    MODEL ar, rhomod;

    *gretl_errmsg = '\0';

    gretl_model_init(&ar);
    gretl_model_init(&rhomod);

    arlist = malloc(pos * sizeof *arlist);
    reglist = malloc((list[0] - pos + 2) * sizeof *reglist);
    reglist2 = malloc((list[0] - pos + 2) * sizeof *reglist2);
    rholist = malloc((pos + 2) * sizeof *rholist);

    if (arlist == NULL || reglist == NULL || reglist2 == NULL ||
	rholist == NULL) {
	ar.errcode = E_ALLOC;
	goto bailout;
    }
    
    arlist[0] = pos - 1;
    for (i=1; i<pos; i++) {
	arlist[i] = list[i];
    }

    reglist2[0] = reglist[0] = list[0] - pos;
    for (i=1; i<=reglist[0]; i++) {
	reglist[i] = list[i + pos];
    }

    rholist[0] = arlist[0] + 1;
    maxlag = ar_list_max(arlist);

    if (const_pos(reglist) > 2) {
	rearrange_list(reglist);
    }

    /* special case: ar 1 ; ... => use CORC */
    if (arlist[0] == 1 && arlist[1] == 1) {
	xx = estimate_rho(reglist, pZ, pdinfo, CORC, &err, OPT_NONE, prn);
	if (err) {
	    ar.errcode = err;
	} else {
	    ar = lsq(reglist, pZ, pdinfo, CORC, OPT_NONE, xx);
	    printmodel(&ar, pdinfo, opt, prn); 
	}
	goto bailout;
    }

    /* first pass: estimate model via OLS: use OPT_M to generate an
       error in case of missing values within sample range 
    */
    ar = lsq(reglist, pZ, pdinfo, OLS, OPT_A | OPT_M, 0.0);
    if (ar.errcode) {
	goto bailout;
    }

    /* allocate space for the uhat terms and transformed data */
    if (dataset_add_series(arlist[0] + 1 + reglist[0], pZ, pdinfo)) {
	ar.errcode = E_ALLOC;
	goto bailout;
    }

    yno = reglist[1];
    t1 = ar.t1; t2 = ar.t2;
    rholist[1] = v;

    pprintf(prn, "%s\n\n", _("Generalized Cochrane-Orcutt estimation"));
    bufspace(17, prn);
    /* xgettext:no-c-format */
    pputs(prn, _("ITER             ESS           % CHANGE"));
    pputs(prn, "\n\n");

    /* now loop while ess is changing */
    diff = 1.0e6;
    ess = 0.0;
    for (iter = 1; iter <= 20 && diff > 0.005; iter++) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t > t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		/* special computation of uhat */
		xx = (*pZ)[yno][t];
		for (j=0; j<reglist[0]-1; j++) {
		    xx -= ar.coeff[j] * (*pZ)[reglist[j+2]][t];
		}
		(*pZ)[v][t] = xx;
	    }
	}		
	for (i=1; i<=arlist[0]; i++) {
	    lag = arlist[i];
	    rholist[1+i] = v + i;
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 + lag || t > t2) {
		    (*pZ)[v+i][t] = NADBL;
		} else {
		    (*pZ)[v+i][t] = (*pZ)[v][t-lag];
		}
	    }
	}

	/* now estimate the rho terms */
	if (iter > 1) {
	    clear_model(&rhomod);
	}
	rhomod = lsq(rholist, pZ, pdinfo, OLS, OPT_A, 0.0);

	/* and rho-transform the data */
	ryno = vc = v + i;
	for (i=1; i<=reglist[0]; i++) {
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 + maxlag || t > t2) {
		    (*pZ)[vc][t] = NADBL;
		} else {
		    xx = (*pZ)[reglist[i]][t];
		    for (j=1; j<=arlist[0]; j++) {
			lag = arlist[j];
			xx -= rhomod.coeff[j-1] * (*pZ)[reglist[i]][t-lag];
		    }
		    (*pZ)[vc][t] = xx;
		}
	    }
	    reglist2[i] = vc++;
	}

	/* estimate the transformed model */
	clear_model(&ar);
	ar = lsq(reglist2, pZ, pdinfo, OLS, OPT_A, 0.0);

        if (iter > 1) {
	    diff = 100 * (ar.ess - ess) / ess;
	}

        if (diff < 0.0) {
	    diff = -diff;
	}

	ess = ar.ess;
	pprintf(prn, "%16c%3d %20f ", ' ', iter, ess);

	if (iter > 1) {
	    pprintf(prn, "%13.3f\n", diff);
	} else {
	    pprintf(prn, "%*s\n", UTF_WIDTH(_("undefined"), 15), 
		     _("undefined")); 
	}
    } /* end "ess changing" loop */

    for (i=0; i<=reglist[0]; i++) {
	ar.list[i] = reglist[i];
    }
    if (gretl_list_has_const(reglist)) {
	ar.ifc = 1;
    }
    if (ar.ifc) {
	ar.dfn -= 1;
    }
    ar.ci = AR;

    /* special computation of fitted values */
    for (t=t1; t<=t2; t++) {
	xx = 0.0;
	for (j=2; j<=reglist[0]; j++) { 
	    xx += ar.coeff[j-2] * (*pZ)[reglist[j]][t];
	}
	ar.uhat[t] = (*pZ)[yno][t] - xx;
	for (j=1; j<=arlist[0]; j++) {
	    if (t - t1 >= arlist[j]) {
		xx += rhomod.coeff[j-1] * ar.uhat[t - arlist[j]];
	    }
	}
	ar.yhat[t] = xx;
    }

    for (t=t1; t<=t2; t++) { 
	ar.uhat[t] = (*pZ)[yno][t] - ar.yhat[t];
    }

    ar.rsq = gretl_corr_rsq(ar.t1, ar.t2, (*pZ)[reglist[1]], ar.yhat);
    ar.adjrsq = 1.0 - ((1.0 - ar.rsq) * (ar.nobs - 1.0) / ar.dfd);

    /* special computation of TSS */
    xx = gretl_mean(ar.t1, ar.t2, (*pZ)[ryno]);
    tss = 0.0;
    for (t=ar.t1; t<=ar.t2; t++) {
	tss += ((*pZ)[ryno][t] - xx) * ((*pZ)[ryno][t] - xx);
    }
    ar.fstt = ar.dfd * (tss - ar.ess) / (ar.dfn * ar.ess);
    ls_aic_bic(&ar);
    ar.dw = dwstat(maxlag, &ar, *pZ);
    ar.rho = rhohat(maxlag, ar.t1, ar.t2, ar.uhat);

    dataset_drop_last_variables(arlist[0] + 1 + reglist[0], pZ, pdinfo);

    if (ar_info_init(&ar, maxlag)) {
	ar.errcode = E_ALLOC;
    } else {
	for (i=0; i<=arlist[0]; i++) { 
	    ar.arinfo->arlist[i] = arlist[i];
	    if (i >= 1) {
		ar.arinfo->rho[i-1] = rhomod.coeff[i-1];
		ar.arinfo->sderr[i-1] = rhomod.sderr[i-1];
	    }
	}
    }
    clear_model(&rhomod);

    if (!ar.errcode) {
	set_model_id(&ar);
	printmodel(&ar, pdinfo, opt, prn);  
    }  

 bailout:

    free(reglist);
    free(reglist2);
    free(rholist);
    free(arlist);

    return ar;
}

/* From 2 to end of list, omits variables with all zero observations
   and re-packs the rest of them */

static void omitzero (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    int v, lv, offset, dropmsg = 0;
    double xx = 0.0;
    char vnamebit[20];

    offset = (pmod->ci == WLS)? 3 : 2;

    for (v=offset; v<=pmod->list[0]; v++) {
        lv = pmod->list[v];
        if (gretl_iszero(pmod->t1, pmod->t2, Z[lv])) {
	    gretl_list_delete_at_pos(pmod->list, v);
	    if (pdinfo->varname[lv][0] != 0) {
		sprintf(vnamebit, "%s ", pdinfo->varname[lv]);
		strcat(gretl_msg, vnamebit);
		dropmsg = 1;
		v--;
	    }
	}
    }

    if (pmod->nwt) {
	int t, wtzero;

	for (v=offset; v<=pmod->list[0]; v++) {
	    lv = pmod->list[v];
	    wtzero = 1;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		xx = Z[lv][t] * Z[pmod->nwt][t];
		if (floatneq(xx, 0.0)) {
		    wtzero = 0;
		    break;
		}
	    }
	    if (wtzero) {
		gretl_list_delete_at_pos(pmod->list, v);
		sprintf(vnamebit, "%s ", pdinfo->varname[lv]);
		strcat(gretl_msg, vnamebit);
		dropmsg = 1;
		v--;
	    }
	}
    }

    if (dropmsg) {
	strcat(gretl_msg, _("omitted because all obs are zero."));
    }
}

static void tsls_omitzero (int *list, const double **Z, int t1, int t2)
{
    int i, v;

    for (i=2; i<=list[0]; i++) {
        v = list[i];
        if (gretl_iszero(t1, t2, Z[v])) {
	    gretl_list_delete_at_pos(list, i);
	    i--;
	}
    }
}

/* ...........................................................*/

static int depvar_zero (int t1, int t2, int yno, int nwt,
			const double **Z)
{
    double y;
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
	y = Z[yno][t];
	if (na(y)) {
	    continue;
	}
	if (nwt) {
	    y *= Z[nwt][t];
	}
	if (y != 0.0) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/* lagdepvar: attempt to detect presence of a lagged dependent
   variable among the regressors -- if found, return the position of
   this lagged var in the list; otherwise return 0
*/

static int 
lagdepvar (const int *list, const double **Z, const DATAINFO *pdinfo) 
{
    char depvar[VNAMELEN], othervar[VNAMELEN];
    char *p;
    int i, t, ret = 0;

    strcpy(depvar, pdinfo->varname[list[1]]);

    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    break;
	}
	strcpy(othervar, pdinfo->varname[list[i]]);
	p = strrchr(othervar, '_');
	if (p != NULL && isdigit(*(p + 1))) {
	    /* looks like a lag */
	    size_t len = strlen(othervar) - strlen(p);

	    if (!strncmp(depvar, othervar, len)) {
		int gotlag = 1;

		/* strong candidate for lagged depvar, but make sure */
		for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) {
		    if (Z[list[1]][t-1] != Z[list[i]][t]) {
			gotlag = 0;
			break;
		    }
		}
		if (gotlag) {
		    ret = i;
		    break;
		}
	    }
	}
    } 

    return ret;
}

/* if 'full' is non-zero, do the whole thing (print the ARCH test
   model, re-estimate if p-value is < .10); otherwise just do
   the ARCH test itself and print the test result
*/

static MODEL 
real_arch_test (MODEL *pmod, int order, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn, int full)
{
    MODEL archmod;
    int *wlist = NULL, *arlist = NULL;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int oldv = pdinfo->v;
    int i, t, nwt, nv, n = pdinfo->n;
    double LM, xx;
    int err = 0;

    *gretl_errmsg = '\0';

    gretl_model_init(&archmod);

    /* assess the lag order */
    if (order < 1 || order > T - pmod->list[0]) {
	archmod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, _("Invalid lag order for arch (%d)"), order);
	err = 1;
    }

    if (!err) {
	/* allocate workspace */
	if (dataset_add_series(order + 1, pZ, pdinfo) || 
	    (arlist = malloc((order + 3) * sizeof *arlist)) == NULL) {
	    err = archmod.errcode = E_ALLOC;
	}
    }

    if (!err) {
	/* start list for aux regression */
	arlist[0] = 2 + order;
	arlist[1] = pdinfo->v - order - 1;
	arlist[2] = 0;

	/* run OLS and get squared residuals */
	archmod = lsq(pmod->list, pZ, pdinfo, OLS, OPT_A | OPT_M, 0.0);
	err = archmod.errcode;
    }

    if (!err) {
	nv = pdinfo->v - order - 1;
	strcpy(pdinfo->varname[nv], "utsq");
	for (t=0; t<n; t++) {
	    (*pZ)[nv][t] = NADBL;
	}
	for (t=archmod.t1; t<=archmod.t2; t++) {
	    xx = archmod.uhat[t];
	    (*pZ)[nv][t] = xx * xx;
	}
	/* also lags of squared resids */
	for (i=1; i<=order; i++) {
	    nv =  pdinfo->v - order + i - 1;
	    arlist[i+2] = nv;
	    sprintf(pdinfo->varname[nv], "utsq_%d", i);
	    for (t=0; t<n; t++) {
		(*pZ)[nv][t] = NADBL;
	    }
	    for (t=archmod.t1+i; t<=archmod.t2; t++) {
		(*pZ)[nv][t] = (*pZ)[arlist[1]][t-i];
	    }
	}

	/* run aux. regression */
	clear_model(&archmod);
	archmod = lsq(arlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	err = archmod.errcode;
    }

    if (!err) {
	archmod.aux = AUX_ARCH;
	archmod.order = order;
	LM = archmod.nobs * archmod.rsq;
	xx = chisq(LM, order);

	if (full) {
	    printmodel(&archmod, pdinfo, OPT_NONE, prn);
	    pprintf(prn, _("No of obs. = %d, unadjusted R^2 = %f\n"),
		    archmod.nobs, archmod.rsq);
	}

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_ARCH);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_TR2);
		model_test_set_order(test, order);
		model_test_set_dfn(test, order);
		model_test_set_value(test, LM);
		model_test_set_pvalue(test, xx);
		maybe_add_test_to_model(pmod, test);
	    }	    
	}

	if (!full) {
	    goto arch_test_exit;
	}

	record_test_result(LM, xx, "ARCH");

	pprintf(prn, _("LM test statistic (%f) is distributed as Chi-square "
		"(%d)\nArea to the right of LM = %f  "), LM, order, xx);

	if (xx > 0.1) {
	    pprintf(prn, "\n%s.\n%s.\n",
		    _("ARCH effect is insignificant at the 10 percent level"),
		    _("Weighted estimation not done"));
	} else {
	    pprintf(prn, "\n%s.\n",
		    _("ARCH effect is significant at the 10 percent level"));
	    /* do weighted estimation */
	    wlist = gretl_list_new(pmod->list[0] + 1);
	    if (wlist == NULL) {
		archmod.errcode = E_ALLOC;
	    } else {
		nwt = wlist[1] = pdinfo->v - 1; /* weight var */
		for (i=2; i<=wlist[0]; i++) {
		    wlist[i] = pmod->list[i-1];
		}
		nv = pdinfo->v - order - 1;
		for (t=archmod.t1; t<=archmod.t2; t++) {
		    xx = archmod.yhat[t];
		    if (xx <= 0.0) {
			xx = (*pZ)[nv][t];
		    }
		    (*pZ)[nwt][t] = 1.0 / sqrt(xx);
		}

		strcpy(pdinfo->varname[nwt], "1/sigma");

		clear_model(&archmod);
		archmod = lsq(wlist, pZ, pdinfo, WLS, OPT_NONE, 0.0);

		archmod.ci = ARCH;
		archmod.order = order;
		printmodel(&archmod, pdinfo, opt, prn);
	    }
	}
    }

 arch_test_exit:

    if (arlist != NULL) free(arlist);
    if (wlist != NULL) free(wlist);

    dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo); 

    return archmod;
}

/**
 * arch_test:
 * @pmod: model to be tested.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain %OPT_O to print covariance matrix, %OPT_S
 *       to save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests @pmod for Auto-Regressive Conditional Heteroskedasticity.  
 * If this effect is significant, re-restimates the model using 
 * weighted least squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch_test (MODEL *pmod, int order, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    return real_arch_test(pmod, order, pZ, pdinfo, opt, prn, 1);
}

/**
 * arch_test_simple:
 * @pmod: model to be tested.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Tests @pmod for Auto-Regressive Conditional Heteroskedasticity.  
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int arch_test_simple (MODEL *pmod, int order, double ***pZ, DATAINFO *pdinfo, 
		      PRN *prn)
{
    MODEL amod;
    int err;

    amod = real_arch_test(pmod, order, pZ, pdinfo, OPT_S, prn, 0);
    err = amod.errcode;
    clear_model(&amod);

    return err;
}

/**
 * arch_model:
 * @list: dependent variable plus list of regressors.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_O to print covariance matrix. (?)
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list via OLS, and test for Auto-
 * Regressive Conditional Heteroskedasticity.  If the latter is
 * significant, re-restimate the model using weighted least
 * squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch_model (const int *list, int order, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    MODEL lmod, amod;

    gretl_model_init(&lmod);
    lmod.list = gretl_list_copy(list);
    if (lmod.list == NULL) {
	lmod.errcode = E_ALLOC;
	return lmod;
    } 

    /* FIXME: vcv option? */
    amod = real_arch_test(&lmod, order, pZ, pdinfo, opt, prn, 1);

    free(lmod.list);

    return amod;
}

/**
 * lad:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using the method of Least
 * Absolute Deviation (LAD).
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lad (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    MODEL lad_model;
    void *handle;
    int (*lad_driver) (MODEL *, double **, DATAINFO *);

    /* run an initial OLS to "set the model up" and check for errors.
       the lad_driver function will overwrite the coefficients etc.
    */

    lad_model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);

    if (lad_model.errcode) {
        return lad_model;
    }

    lad_driver = get_plugin_function("lad_driver", &handle);

    if (lad_driver == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	lad_model.errcode = E_FOPEN;
	return lad_model;
    }

    (*lad_driver) (&lad_model, *pZ, pdinfo);
    close_plugin(handle);

    set_model_id(&lad_model);

    return lad_model;
}

/**
 * arma:
 * @list: dependent variable, AR and MA orders, and any exogenous
 * regressors.
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @opt: options: may include %OPT_S to suppress intercept, %OPT_V
 * for verbose results, %OPT_X to use X-12-ARIMA.
 * @PRN: for printing details of iterations (or %NULL). 
 *
 * Calculate ARMA estimates, using either native gretl code or
 * by invoking X-12-ARIMA.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arma (const int *list, const double **Z, const DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn)
{
    MODEL armod;
    void *handle;
    MODEL (*arma_func) (const int *, const double **, const DATAINFO *, 
			 gretlopt, PRN *);

    *gretl_errmsg = '\0';

    if (opt & OPT_X) {
	arma_func = get_plugin_function("arma_x12_model", &handle);
    } else {
	arma_func = get_plugin_function("arma_model", &handle);
    }

    if (arma_func == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	gretl_model_init(&armod);
	armod.errcode = E_FOPEN;
	return armod;
    }

    armod = (*arma_func) (list, Z, pdinfo, opt, prn);

    close_plugin(handle);
    set_model_id(&armod);

    return armod;
} 

/**
 * logistic_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @param: 
 *
 * Estimate the model given in @list using the logistic transformation
 * of the dependent variable.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL logistic_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		      const char *param)
{
    MODEL lmod;
    void *handle;
    MODEL (*logistic_estimate) (const int *, double ***, DATAINFO *, const char *);

    *gretl_errmsg = '\0';

    logistic_estimate = get_plugin_function("logistic_estimate", &handle);
    if (logistic_estimate == NULL) {
	gretl_model_init(&lmod);
	lmod.errcode = E_FOPEN;
	return lmod;
    }

    lmod = (*logistic_estimate) (list, pZ, pdinfo, param);

    close_plugin(handle);

    set_model_id(&lmod);

    return lmod;
}

/**
 * tobit_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: printing struct for iteration info (or NULL is this is not
 * wanted).
 *
 * Estimate the model given in @list using Tobit.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tobit_model (const int *list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL tmod;
    void *handle;
    MODEL (* tobit_estimate) (const int *, double ***, DATAINFO *, PRN *);

    *gretl_errmsg = '\0';

    tobit_estimate = get_plugin_function("tobit_estimate", &handle);
    if (tobit_estimate == NULL) {
	gretl_model_init(&tmod);
	tmod.errcode = E_FOPEN;
	return tmod;
    }

    tmod = (*tobit_estimate) (list, pZ, pdinfo, prn);

    close_plugin(handle);

    set_model_id(&tmod);

    return tmod;
}

static int get_offset_var (int *list)
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
 * poisson_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: printing struct for iteration info (or NULL is this is not
 * wanted).
 *
 * Estimate the Poisson regression model given in @list using ML.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL poisson_model (const int *list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL pmodel;
    void *handle;
    int *listcpy;
    int offvar;
    int (* poisson_estimate) (MODEL *, int, double ***, DATAINFO *, PRN *);

    *gretl_errmsg = '\0';

    gretl_model_init(&pmodel);

    listcpy = gretl_list_copy(list);
    if (listcpy == NULL) {
	pmodel.errcode = E_ALLOC;
        return pmodel;
    }

    offvar = get_offset_var(listcpy);

    /* run an initial OLS to "set the model up" and check for errors.
       the poisson_estimate_driver function will overwrite the
       coefficients etc.
    */

    pmodel = lsq(listcpy, pZ, pdinfo, OLS, OPT_A, 0.0);
    free(listcpy);

    if (pmodel.errcode) {
        return pmodel;
    }

    poisson_estimate = get_plugin_function("poisson_estimate", &handle);

    if (poisson_estimate == NULL) {
	pmodel.errcode = E_FOPEN;
	return pmodel;
    }

    (*poisson_estimate) (&pmodel, offvar, pZ, pdinfo, prn);

    close_plugin(handle);

    set_model_id(&pmodel);

    return pmodel;
}

/**
 * garch:
 * @list: dependent variable plus arch and garch orders.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: can specify robust standard errors and VCV.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Calculate GARCH estimates.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL garch (const int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn)
{
    MODEL gmod;
    void *handle;
    PRN *myprn;
    MODEL (*garch_model) (const int *, double ***, DATAINFO *, PRN *,
			  gretlopt);

    *gretl_errmsg = '\0';

    garch_model = get_plugin_function("garch_model", &handle);

    if (garch_model == NULL) {
	gretl_model_init(&gmod);
	gmod.errcode = E_FOPEN;
	return gmod;
    }

    if (opt & OPT_V) {
	myprn = prn;
    } else {
	myprn = NULL;
    }

    gmod = (*garch_model) (list, pZ, pdinfo, myprn, opt);

    close_plugin(handle);

    set_model_id(&gmod);

    return gmod;
} 

MODEL pooled (const int *list, double ***pZ, DATAINFO *pdinfo,
	      gretlopt opt, PRN *prn)
{
    MODEL wmod;

    *gretl_errmsg = '\0';

    if (opt & OPT_W) {
	void *handle;
	MODEL (*panel_wls_by_unit) (const int *, double ***, DATAINFO *,
				    gretlopt, PRN *);

	panel_wls_by_unit = get_plugin_function("panel_wls_by_unit", &handle);

	if (panel_wls_by_unit == NULL) {
	    gretl_model_init(&wmod);
	    wmod.errcode = E_FOPEN;
	    return wmod;
	}

	wmod = (*panel_wls_by_unit) (list, pZ, pdinfo, opt, prn);

	close_plugin(handle);
    } else {
	wmod = lsq(list, pZ, pdinfo, POOLED, opt, 0.0);
    }

    return wmod;
}

int groupwise_hetero_test (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			   PRN *prn)
{
    MODEL wmod;
    int err;

    if (!dataset_is_panel(pdinfo)) {
	strcpy(gretl_errmsg, _("This test is only available for panel data"));
	return 1;
    }

    wmod = pooled(pmod->list, pZ, pdinfo, OPT_W | OPT_T | OPT_A, prn);
    err = wmod.errcode;

    if (!err) {
	gretl_model_set_auxiliary(&wmod, AUX_GROUPWISE);
	printmodel(&wmod, pdinfo, OPT_NONE, prn);
    }

    clear_model(&wmod);

    return err;
}

