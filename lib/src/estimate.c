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
#include "gretl_private.h"
#include "gretl_list.h"

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
#define STATZERO  0.5e-14
#define ESSZERO   1.0e-22

extern void _print_rho (int *arlist, const MODEL *pmod, 
			int c, PRN *prn);

/* private function prototypes */
static int form_xpxxpy (const int *list, int t1, int t2, 
			double **Z, int nwt, double rho, int pwe,
			double *xpx, double *xpy);
static void regress (MODEL *pmod, double *xpy, double **Z, 
		     int n, double rho);
static int cholbeta (double *xpx, double *xpy, double *coeff, double *rss,
		     int nv);
static void diaginv (double *xpx, double *xpy, double *diag, int nv);

static double dwstat (int order, MODEL *pmod, double **Z);
static double rhohat (int order, int t1, int t2, const double *uhat);
static double altrho (int order, int t1, int t2, const double *uhat);
static int hatvar (MODEL *pmod, int n, double **Z);
static void dropwt (int *list);
static int get_aux_uhat (MODEL *pmod, double *uhat1, double ***pZ, 
			 DATAINFO *pdinfo);
static void omitzero (MODEL *pmod, const DATAINFO *pdinfo, double **Z);
static void tsls_omitzero (int *list, double **Z, int t1, int t2);
static int zerror (int t1, int t2, int yno, int nwt, double ***pZ);
static int lagdepvar (const int *list, const DATAINFO *pdinfo, 
		       double ***pZ);
static int tsls_match (const int *list1, const int *list2, int *newlist);
static double wt_dummy_mean (const MODEL *pmod, double **Z); 
static double wt_dummy_stddev (const MODEL *pmod, double **Z);
/* end private protos */

static int reorganize_uhat_yhat (MODEL *pmod) 
{
    int t, g;
    MISSOBS *mobs = (MISSOBS *) pmod->data;
    double *tmp;

    tmp = malloc(pmod->nobs * sizeof *tmp);
    if (tmp == NULL) return 1;

    /* first do uhat */
    for (t=0; t<pmod->nobs; t++)
	tmp[t] = pmod->uhat[t];

    g = 0;
    for (t=pmod->t1; t<=pmod->t2 + mobs->misscount; t++) {
	if (mobs->missvec[t - pmod->t1]) pmod->uhat[t] = NADBL;
	else pmod->uhat[t] = tmp[g++];
    }

    /* then yhat */
    for (t=0; t<pmod->nobs; t++)
	tmp[t] = pmod->yhat[t];

    g = 0;
    for (t=pmod->t1; t<=pmod->t2 + mobs->misscount; t++) {
	if (mobs->missvec[t - pmod->t1]) pmod->yhat[t] = NADBL;
	else pmod->yhat[t] = tmp[g++];
    }

    free(tmp);
    return 0;
}

/* with daily data, try eliminating the missing obs? */

static int handle_missing_daily_obs (MODEL *pmod, double **Z, 
				     DATAINFO *pdinfo)
{
    int misscount;
    char *missvec = missobs_vector(Z, pdinfo, &misscount);
    MISSOBS *mobs = NULL;

    if (missvec == NULL || (mobs = malloc(sizeof *mobs)) == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    } else {
	repack_missing(Z, pdinfo, missvec, misscount);
	pmod->t2 -= misscount;
	mobs->misscount = misscount;
	mobs->missvec = missvec;
	pmod->data = mobs;
    }

    return 0;
}

static int ar_info_init (MODEL *pmod, int nterms)
{
    int i;

    pmod->arinfo = malloc(sizeof *pmod->arinfo);
    if (pmod->arinfo == NULL) return 1;

    pmod->arinfo->arlist = malloc(nterms * sizeof *pmod->arinfo->arlist);
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
	pmod->arinfo->arlist[i] = 0;
	pmod->arinfo->sderr[i] = pmod->arinfo->rho[i] = NADBL;
    }

    return 0;
}

static int get_model_df (MODEL *pmod)
{
    pmod->ncoeff = pmod->list[0] - 1;

    pmod->dfd = pmod->nobs - pmod->ncoeff;
    if (pmod->dfd < 0) {
	pmod->errcode = E_DF;
        sprintf(gretl_errmsg, _("No. of obs (%d) is less than no. "
		"of parameters (%d)"), pmod->nobs, pmod->ncoeff);
	return 1;
    }

    pmod->dfn = pmod->ncoeff - pmod->ifc;

    return 0;
}

static int compute_ar_stats (MODEL *pmod, const double **Z, double rho)
{
    int i, t, yno = pmod->list[1];
    double x, pw1 = 0.0;

    if (ar_info_init(pmod, 2)) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    if (pmod->ci == PWE) {
	pw1 = sqrt(1.0 - rho * rho);
    }

    pmod->arinfo->arlist[0] = pmod->arinfo->arlist[1] = 1;
    pmod->arinfo->rho[1] = rho;
    gretl_model_set_double(pmod, "rho_in", rho);

    t = pmod->t1;
    if (pmod->ci == PWE) {
	x = pw1 * Z[yno][t];
	for (i=pmod->ifc; i<pmod->ncoeff; i++) {
	    x -= pmod->coeff[i] * pw1 * Z[pmod->list[i+2]][t];
	}
	if (pmod->ifc) {
	    x -= pw1 * pmod->coeff[0];
	}
	pmod->uhat[t] = x;
	pmod->yhat[t] = Z[yno][t] - x;
    } else {
	pmod->uhat[t] = NADBL;
	pmod->yhat[t] = NADBL;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	x = Z[yno][t] - rho * Z[yno][t-1];
	for (i=0; i<pmod->ncoeff; i++) {
	    x -= pmod->coeff[i] * 
		(Z[pmod->list[i+2]][t] - 
		 rho * Z[pmod->list[i+2]][t-1]);
	}
	pmod->uhat[t] = x;
	pmod->yhat[t] = Z[yno][t] - x;
    }

    pmod->rsq = 
	corrrsq(pmod->t2 - pmod->t1 + 1, &Z[yno][pmod->t1], 
		pmod->yhat + pmod->t1);
    pmod->adjrsq = 
	1.0 - ((1.0 - pmod->rsq) * (pmod->t2 - pmod->t1) / 
	       (double) pmod->dfd);

    return 0;
}

/* Calculation of WLS stats in agreement with GNU R */

static void get_wls_stats (MODEL *pmod, const double **Z)
{
    int t, wobs = pmod->nobs, yno = pmod->list[1];
    double x, w2, wmean = 0.0, wsum = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) { 
	if (Z[pmod->nwt][t] == 0) {
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
	if (Z[pmod->nwt][t] == 0.0) continue;
	w2 = Z[pmod->nwt][t] * Z[pmod->nwt][t]; 
	x += w2 * (Z[yno][t] - wmean) * (Z[yno][t] - wmean);
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
	    if (floateq(Z[pmod->nwt][t], 0.0)) 
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	}
    } else {
	double x;

	pmod->ess_wt = pmod->ess;
	pmod->sigma_wt = pmod->sigma;
	pmod->ess = 0.0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
		pmod->nobs -= 1;
	    } else {
		pmod->yhat[t] /= Z[pmod->nwt][t];
		x = pmod->uhat[t] /= Z[pmod->nwt][t];
		pmod->ess += x * x;
	    }
	}
	pmod->sigma = sqrt(pmod->ess / pmod->dfd);
    }
}

static void model_stats_init (MODEL *pmod)
{
    pmod->ess = pmod->ess_wt = NADBL;
    pmod->sigma = pmod->sigma_wt = NADBL;
    pmod->fstt = pmod->lnL = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
}

/**
 * lsq:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: command index (see gretl_commands.h)
 * @opts: option flags: 
 *   if & OPT_R compute robust standard errors;
 *   if & OPT_C force use of Cholesky decomp;
 *   if & OPT_A treat as auxiliary regression (don't bother checking
 *     for presence of lagged dependent var, don't augment model count);
 *   if & OPT_P use Prais-Winsten for first obs.
 * @rho: coefficient for rho-differencing the data (0.0 for no
 * differencing)
 *
 * Computes least squares estimates of the model specified by @list,
 * using an estimator determined by the value of @ci.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lsq (LIST list, double ***pZ, DATAINFO *pdinfo, 
	   int ci, gretlopt opts, double rho)
{
    int l0, yno, i;
    int effobs = 0;
    int missv = 0, misst = 0;
    int ldepvar = 0;
    int use_qr = get_use_qr();
    int pwe = (ci == PWE || (opts & OPT_P));
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
    if (ci == HCCM) {
	return hccm_func(list, pZ, pdinfo);
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
    mdl.list = copylist(list);
    if (mdl.list == NULL) {
        mdl.errcode = E_ALLOC;
        return mdl;
    }

    mdl.t1 = pdinfo->t1;
    mdl.t2 = pdinfo->t2;
    mdl.ci = ci;

    /* Doing weighted least squares? */
    mdl.nwt = 0;
    if (ci == WLS) { 
	mdl.nwt = mdl.list[1];
	if (gretl_iszero(mdl.t1, mdl.t2, (*pZ)[mdl.nwt])) {
	    mdl.errcode = E_WTZERO;
	    return mdl;
	}
	effobs = isdummy((*pZ)[mdl.nwt], mdl.t1, mdl.t2);
	if (effobs) {
	    gretl_model_set_int(&mdl, "wt_dummy", 1);
	}
    }

    /* check for missing obs in sample */
    if ((missv = adjust_t1t2(&mdl, mdl.list, &mdl.t1, &mdl.t2, 
			     (const double **) *pZ, &misst))) {
	if (!dated_daily_data(pdinfo)) {
	    sprintf(gretl_errmsg, _("Missing value encountered for "
		    "variable %d, obs %d"), missv, misst);
	    mdl.errcode = E_DATA;
	    return mdl;
	} else if (handle_missing_daily_obs(&mdl, *pZ, pdinfo)) {
	    return mdl;
	}
    }

    if (ci == WLS) dropwt(mdl.list);
    yno = mdl.list[1];
    
    /* check for availability of data */
    if (mdl.t1 < 0 || mdl.t2 > pdinfo->n - 1) {
        mdl.errcode = E_NODATA;
        goto lsq_abort;
    }                   

    /* check for unknown vars in list */
    for (i=1; i<=mdl.list[0]; i++) {
        if (mdl.list[i] > pdinfo->v - 1) {
            mdl.errcode = E_UNKVAR;
            goto lsq_abort;
        }
    }       

    /* check for zero dependent var */
    if (zerror(mdl.t1, mdl.t2, yno, mdl.nwt, pZ)) {  
        mdl.errcode = E_ZERO;
        goto lsq_abort; 
    } 

    /* drop any vars that are all zero and repack the list */
    omitzero(&mdl, pdinfo, *pZ);

    /* if regressor list contains a constant, place it first */
    i = gretl_hasconst(mdl.list);
    mdl.ifc = (i > 1);
    if (i > 2) rearrange_list(mdl.list);

    /* Check for presence of lagged dependent variable? 
       Don't bother if this is an auxiliary regression. */
    if (!(opts & OPT_A)) {
	ldepvar = lagdepvar(mdl.list, pdinfo, pZ);
	if (ldepvar) {
	    gretl_model_set_int(&mdl, "ldepvar", ldepvar);
	}
    }

    /* AR1: advance the starting observation by one? */
    if (rho != 0.0 && !pwe) {
	mdl.t1 += 1;
    }

    l0 = mdl.list[0];  /* holds 1 + number of coeffs */
    mdl.ncoeff = l0 - 1; 
    if (effobs) mdl.nobs = effobs;
    else mdl.nobs = mdl.t2 - mdl.t1 + 1;

    /* check degrees of freedom */
    if (get_model_df(&mdl)) {
        goto lsq_abort; 
    }

    if (dataset_is_time_series(pdinfo)) {
	opts |= OPT_T;
    }

    if ((opts & OPT_R) || (use_qr && !(opts & OPT_C))) { 
	mdl.rho = rho;
	gretl_qr_regress(&mdl, (const double **) *pZ, pdinfo->n, opts);
    } else {
	int l = l0 - 1;
	int nxpx = l * (l + 1) / 2;

	xpy = malloc((l0 + 1) * sizeof *xpy);
	mdl.xpx = malloc(nxpx * sizeof *mdl.xpx);
	mdl.coeff = malloc(mdl.ncoeff * sizeof *mdl.coeff);
	if (xpy == NULL || mdl.xpx == NULL || mdl.coeff == NULL) {
	    mdl.errcode = E_ALLOC;
	    return mdl;
	}

	for (i=0; i<=l0; i++) xpy[i] = 0.0;
	for (i=0; i<nxpx; i++) mdl.xpx[i] = 0.0;

	/* calculate regression results, Cholesky style */
	form_xpxxpy(mdl.list, mdl.t1, mdl.t2, *pZ, mdl.nwt, rho, 
		    pwe, mdl.xpx, xpy);

	regress(&mdl, xpy, *pZ, pdinfo->n, rho);
	free(xpy);
    }

    if (mdl.errcode) {
	goto lsq_abort;
    }

    /* get the mean and sd of depvar and make available */
    if (mdl.ci == WLS && gretl_model_get_int(&mdl, "wt_dummy")) {
	mdl.ybar = wt_dummy_mean(&mdl, *pZ);
	mdl.sdy = wt_dummy_stddev(&mdl, *pZ);
    } else {
	mdl.ybar = gretl_mean(mdl.t1, mdl.t2, (*pZ)[yno]);
	mdl.sdy = gretl_stddev(mdl.t1, mdl.t2, (*pZ)[yno]);
	if (fabs(mdl.ybar) < STATZERO) mdl.ybar = 0.0;
    }

    /* Doing an autoregressive procedure? */
    if (ci == CORC || ci == HILU || ci == PWE) {
	if (compute_ar_stats(&mdl, (const double **) *pZ, rho)) 
	    goto lsq_abort;
    }

    /* weighted least squares: fix fitted values, ESS, sigma */
    if (ci == WLS) {
	get_wls_stats(&mdl, (const double **) *pZ);
	fix_wls_values(&mdl, *pZ);
    }

    if (opts & OPT_T) {
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
    gretl_aic_etc(&mdl);

 lsq_abort:

    /* If we eliminated any missing observations, restore
       them now */
    if (mdl.data != NULL) {
	MISSOBS *mobs = (MISSOBS *) mdl.data;

	undo_repack_missing(*pZ, pdinfo, mobs->missvec,
			    mobs->misscount);
	reorganize_uhat_yhat(&mdl);
    }

    if (!(opts & OPT_A)) {
	set_model_id(&mdl);
    }

    return mdl;
}

/* .......................................................... */

static int form_xpxxpy (const int *list, int t1, int t2, 
			double **Z, int nwt, double rho, int pwe,
			double *xpx, double *xpy)
/*
        This function forms the X'X matrix and X'y vector

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
        xpy[0] = sum of y's
        xpy[list[0]] = y'y
*/
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
    }
    else if (nwt) {
	/* weight the data */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    z1 = Z[nwt][t];
		    x += z1 * z1 * Z[li][t] * Z[lj][t];
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}   
		xpx[m++] = x;
	    }
	    x = 0;
	    for (t=t1; t<=t2; t++) {
		z1 = Z[nwt][t];
		x += z1 * z1 * Z[yno][t] * Z[li][t];
	    }
	    xpy[i-1] = x;
	}
    }
    else {
	/* no quasi-differencing or weighting wanted */
	for (i=2; i<=l0; i++) {
	    li = list[i];
	    for (j=i; j<=l0; j++) {
		lj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    x += Z[li][t] * Z[lj][t];
		}
		if (floateq(x, 0.0) && li == lj)  {
		    return li;
		}
		xpx[m++] = x;
	    }
	    x = 0.0;
	    for (t=t1; t<=t2; t++) {
		x += Z[yno][t] * Z[li][t];
	    }
	    xpy[i-1] = x;
	}
    }

    return 0; 
}

/* .......................................................... */

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

/* .......................................................... */

static void compute_r_squared (MODEL *pmod, double *y)
{
    pmod->rsq = 1 - (pmod->ess / pmod->tss);

    if (pmod->dfd > 0) {
	double den = pmod->tss * pmod->dfd;

	if (pmod->ifc) {
	    pmod->adjrsq = 1 - (pmod->ess * (pmod->nobs - 1) / den);
	} else {
	    pmod->rsq = corrrsq(pmod->nobs, y, pmod->yhat + pmod->t1);
	    pmod->adjrsq = 
		1 - ((1 - pmod->rsq) * (pmod->nobs - 1) / pmod->dfd);
	} 
    }
}

/* .......................................................... */

static void regress (MODEL *pmod, double *xpy, double **Z, 
		     int n, double rho)
/*
        This function takes xpx, the X'X matrix output
        by form_xpxxpy(), and xpy (X'y), and
        computes ols estimates and associated statistics.

        n = no. of observations per series in data set
        ifc = 1 if constant is present else = 0

        ess = error sum of squares
        sigma = standard error of regression
        fstt = F-statistic
        coeff = vector of regression coefficients
        sderr = vector of standard errors of regression coefficients
*/
{
    int v, yno = pmod->list[1];
    double ysum, ypy, zz, rss = 0.0;
    double sgmasq = 0.0;
    double *diag = NULL;
    int i, err = 0;

    pmod->sderr = malloc(pmod->ncoeff * sizeof *pmod->sderr);
    pmod->yhat = malloc(n * sizeof *pmod->yhat);
    pmod->uhat = malloc(n * sizeof *pmod->uhat);

    if (pmod->sderr == NULL || pmod->yhat == NULL || pmod->uhat == NULL) {
        pmod->errcode = E_ALLOC;
        return;
    }

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

    /* the idea below was broken */
    if (pmod->ess < ESSZERO && pmod->ess > (-ESSZERO)) {
	pmod->ess = 0.0;
    } else if (pmod->ess < 0.0) { 
        /*  pmod->errcode = E_ESS; */ 
	sprintf(gretl_errmsg, _("Error sum of squares (%g) is not > 0"),
		pmod->ess);
        return; 
    }

    if (pmod->dfd == 0) {
	pmod->sigma = 0.0;
	pmod->adjrsq = NADBL;
    } else {
	sgmasq = pmod->ess / pmod->dfd;
	pmod->sigma = sqrt(sgmasq);
    }

    if (floatlt(pmod->tss, 0.0) || floateq(pmod->tss, 0.0)) {
       pmod->rsq = pmod->adjrsq = NADBL;
    } 

    hatvar(pmod, n, Z); 
    if (pmod->errcode) return;

    if (pmod->tss > 0.0) {
	compute_r_squared(pmod, &Z[yno][pmod->t1]);
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
	if (pmod->fstt < 0.0) pmod->fstt = 0.0;
    }

    diag = malloc(pmod->ncoeff * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    diaginv(pmod->xpx, xpy, diag, pmod->ncoeff);

    for (v=0; v<pmod->ncoeff; v++) { 
       pmod->sderr[v] = pmod->sigma * sqrt(diag[v]); 
    }

    free(diag); 
    
    return;  
}

/* .......................................................... */

static int cholbeta (double *xpx, double *xpy, double *coeff, double *rss,
		     int nv)
/*

  This function does an in-place Choleski decomposition of xpx
  (lower triangular matrix stacked in columns) and then
  solves the normal equations for coeff.  xpx is X'X
  on input and Choleski decomposition on output; xpy is
  the X'y vector on input and Choleski-transformed t vector
  on output. coeff is the vector of estimated coefficients; 
  nv is the number of regression coefficients including the 
  constant. 

  The number of floating-point operations is basically 3.5 * nv^2
  plus (nv^3) / 3.

*/

{
    int i, j, k, kk, l, jm1;
    double e, d, d1, d2, test, xx;

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
	    fprintf(stderr, "cholbeta: test = %g\n", test);
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
	for (j=0; j<nv-1; j++) coeff[j] = 0.0;
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
    } 

    return 0; 
}

/* ...............................................................    */

static void diaginv (double *xpx, double *xpy, double *diag, int nv)
/*
        Solves for the diagonal elements of the X'X inverse matrix.

        xpx = Cholesky-decomposed X'X matrix (input)
        xpy = X'y vector (input) used as work array
        diag = diagonal elements of X'X (output)
        nv = number of regression coefficients
*/
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
 * @pmod: gretl MODEL.
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

    if (pmod->vcv != NULL) return 0;
    if (pmod->xpx == NULL) return 1;

    mst = nxpx;
    kk = nxpx - 1;

    pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) return E_ALLOC;

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
	if (i > nv - 2) continue;
	for (j=i+1; j<nv; j++) {
	    icnt = i+1;
	    kj -= j;
	    d = 0.0;
	    m = mst + 1;
	    for (k=0; k<=j-1; k++) {
		if(icnt > 0) {
		    dec = 1;
		    icnt--;
		}
		else dec = k;
		m -= dec;
		l = kj + i - k;
		d += pmod->vcv[m-1] * pmod->xpx[l];
	    }
	    pmod->vcv[kj] = (-1.0) * d * pmod->xpx[l-1];
	}
    }

    if (pmod->ci == CUSUM) return 0;

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

/* ............................................................... */

/**
 * get_vcv:
 * @pmod: pointer to model.
 * 
 * Save the variance-covariance matrix for the parameter
 * estimates in @pmod.
 *
 * Returns: VCV struct or NULL on error.
 */

VCV *get_vcv (MODEL *pmod)
{
    int i, nv = pmod->ncoeff;
    VCV *vcv;

    vcv = malloc(sizeof *vcv);
    if (vcv == NULL) return NULL;

    vcv->list = malloc((nv + 1) * sizeof *vcv->list);
    if (vcv->list == NULL) {
	free(vcv);
	return NULL;
    }

    vcv->list[0] = nv;
    for (i=1; i<=nv; i++) {
	vcv->list[i] = pmod->list[i+1];
    }

    if (pmod->vcv == NULL && makevcv(pmod)) {
	free(vcv->list);
	free(vcv);
	return NULL;
    }

    /* calculate number of elements in vcv */
    nv = (nv * nv + nv) / 2;

    /* copy vcv */
    vcv->vec = copyvec(pmod->vcv, nv + 1);
    if (vcv->vec == NULL) {
	free(vcv->list);
	free(vcv);
	return NULL;
    }

    vcv->ci = pmod->ci;
    
    return vcv;
}

/* ............................................................... */

void free_vcv (VCV *vcv)
{
    free(vcv->vec);
    free(vcv->list);
    free(vcv);
}

/* ............................................................... */

static double dwstat (int order, MODEL *pmod, double **Z)
/*  computes durbin-watson statistic
    order is the order of autoregression, 0 for OLS.
*/
{
    double diff, ut, ut1;
    double diffsq = 0.0;
    int t;

    if (order) order--;

    if (pmod->ess <= 0.0) return NADBL;

    for (t=pmod->t1+1+order; t<=pmod->t2; t++)  {
        ut = pmod->uhat[t];
        ut1 = pmod->uhat[t-1];
        if (na(ut) || na(ut1) ||
	    (pmod->nwt && (floateq(Z[pmod->nwt][t], 0.0) || 
			   floateq(Z[pmod->nwt][t-1], 0.0)))) 
	    continue;
        diff = ut - ut1;
        diffsq += diff * diff;
    }

    return diffsq / pmod->ess;
}

/* ......................................................  */

static double rhohat (int order, int t1, int t2, const double *uhat)
/*  computes first order serial correlation coefficient
    order is the order of autoregression, 0 for OLS.
*/
{
    double ut, ut1, uu = 0.0, xx = 0.0;
    double rho;
    int t;

    if (order) order--;

#if 1 /* the original */
    for (t=t1+order+1; t<=t2; t++) { 
        ut = uhat[t];
        ut1 = uhat[t-1];
        if (na(ut) || na(ut1)) continue;
        uu += ut * ut1;
        xx += ut1 * ut1;
    }
#else
    for (t=t1; t<=t2; t++) { 
        ut = uhat[t];
	if (na(ut)) continue;
	if (t > t1) {
	    ut1 = uhat[t-1];
	    if (!na(ut1)) {
		uu += ut * ut1;
	    }
	}
        xx += ut * ut;
    }
#endif

    if (floateq(xx, 0.0)) return NADBL;

    rho = uu / xx;
    if (rho > 1.0 || rho < -1.0) {
	rho = altrho(order, t1, t2, uhat);
    }

    return rho;
}

/* .........................................................   */

static double altrho (int order, int t1, int t2, const double *uhat)
/* alternative calculation of rho */
{
    double *ut, *ut1;    
    int t, n, len = t2 - (t1 + order) + 1;
    double uh, uh1, rho;

    ut = malloc(len * sizeof *ut);
    if (ut == NULL) return NADBL;

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

    rho = gretl_corr(n, ut, ut1);

    free(ut);
    free(ut1);

    return rho;
}

/* ........................................................... */

double corrrsq (int nobs, const double *y, const double *yhat)
/* finds alternative R^2 value when there's no intercept */
{
    double x;

    x = gretl_corr(nobs, y, yhat);
    if (na(x)) return NADBL;
    else return x * x;
}

/* ........................................................... */

/* compute fitted values and residuals */

#undef OUT_OF_SAMPLE_OK

static int hatvar (MODEL *pmod, int n, double **Z)
{
    int xno, i, t;
    int yno = pmod->list[1];
    double x;

#ifdef OUT_OF_SAMPLE_OK
    double y;

    for (t=0; t<n; t++) {
	int missing = 0;

	y = Z[yno][t];
	if (na(y)) {
	    pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    continue;
	}
	pmod->yhat[t] = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            xno = pmod->list[i+2];
	    x = Z[xno][t];
	    if (na(x)) {
		missing = 1;
		break;
	    }
	    if (pmod->nwt) {
		if (na(Z[pmod->nwt][t])) {
		    missing = 1;
		    break;
		}
		x *= Z[pmod->nwt][t];
	    }
            pmod->yhat[t] += pmod->coeff[i] * x;
        }
	if (missing) {
	    pmod->yhat[t] = pmod->uhat[t] = NADBL;
	} else {
	    if (pmod->nwt) { 
		y *= Z[pmod->nwt][t];
	    }
	    pmod->uhat[t] = y - pmod->yhat[t]; 
	}               
    }
#else
    for (t=0; t<n; t++) {
	pmod->yhat[t] = pmod->uhat[t] = NADBL;
    }
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->yhat[t] = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            xno = pmod->list[i+2];
	    x = Z[xno][t];
	    if (pmod->nwt) x *= Z[pmod->nwt][t];
            pmod->yhat[t] += pmod->coeff[i] * x;
        }
	x = Z[yno][t];
	if (pmod->nwt) x *= Z[pmod->nwt][t];
        pmod->uhat[t] = x - pmod->yhat[t];                
    }
#endif

    return 0;
}

/* ........................................................... */

static void dropwt (int *list)
/* drop the weight var from the list of regressors (WLS) */
{
    int i;

    list[0] -= 1;
    for (i=1; i<=list[0]; i++) {
	list[i] = list[i+1];
    }
}

static int hilu_plot (double *ssr, double *rho, int n, 
		      PATHS *ppaths)
{
    FILE *fp;
    int i;

    if (ppaths == NULL) return 1;

    if (gnuplot_init(ppaths, PLOT_REGULAR, &fp)) return E_FOPEN; 

    fputs("# hildreth-lu\n", fp);
    fputs("set xlabel 'rho'\n", fp);
    fprintf(fp, "set ylabel '%s'\n", _("ESS"));
    fputs("set nokey\n", fp);
    fputs("set xrange [-1.0:1.0]\n", fp);
    fprintf(fp, "plot '-' using 1:2 w impulses\n");

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
    gnuplot_display(ppaths);

    return 0;
}

static double autores (MODEL *pmod, double rho, const double **Z,
		       int opt)
{
    int t, v, t1 = pmod->t1;
    double x, num = 0.0, den = 0.0;

    if (opt == CORC || opt == HILU) t1--;

    for (t=t1; t<=pmod->t2; t++) {
	x = Z[pmod->list[1]][t];
	for (v=0; v<pmod->ncoeff; v++) {
	    x -= pmod->coeff[v] * Z[pmod->list[v+2]][t];
	}
	pmod->uhat[t] = x;
	if (t > t1) {
	    num += pmod->uhat[t] * pmod->uhat[t-1];
	    den += pmod->uhat[t-1] * pmod->uhat[t-1];
	}
    }

    return num / den;
}

/**
 * hilu_corc:
 * @toprho: pointer to receive final rho value.
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ppaths: pointer to gretl paths info struct
 * @batch: = 1 if in batch mode
 * @opt: option flag: CORC for Cochrane-Orcutt, HILU for Hildreth-Lu,
 *                    PWE for Prais-Winsten estimator
 * @prn: gretl printing struct
 *
 * Estimate the model given in @list using either the Cochrane-Orcutt
 * procedure or Hildreth-Lu (for first-order serial correlation).
 * Print a trace of the search for the appropriate quasi-differencing
 * coefficient, rho.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int hilu_corc (double *toprho, LIST list, double ***pZ, DATAINFO *pdinfo,
	       PATHS *ppaths, int batch, int opt, PRN *prn)
{
    double rho = 0.0, rho0 = 0.0, diff;
    double finalrho = 0.0, essmin = 1.0e8;
    double ess, ssr[199], rh[199]; 
    int iter, nn = 0, err = 0;
    gretlopt lsqopt = OPT_NONE;
    MODEL corc_model;

    *gretl_errmsg = '\0';

    gretl_model_init(&corc_model);

    if (opt == PWE) lsqopt |= OPT_P;

    if (opt == HILU) { /* Do Hildreth-Lu first */
	for (rho = -0.990, iter = 0; rho < 1.0; rho += .01, iter++) {
	    clear_model(&corc_model);
	    corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
	    if ((err = corc_model.errcode)) {
		clear_model(&corc_model);
		return err;
	    }
	    ess = corc_model.ess;
	    if (batch) {
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
		    if (nn % 4 == 0) pputc(prn, '\n');
		    else bufspace(3, prn);
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
		if ((err = corc_model.errcode)) {
		    clear_model(&corc_model);
		    break;
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
		if ((err = corc_model.errcode)) {
		    clear_model(&corc_model);
		    break;
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
	if (batch) {
	    graphyzx(NULL, ssr, NULL, rh, nn, "ESS", "RHO", NULL, 0, prn); 
	    pputs(prn, "\n");
	} else {
	    hilu_plot(ssr, rh, nn, ppaths);
	}
	pputs(prn, _("\nFine-tune rho using the CORC procedure...\n\n")); 
    } else { 
	/* Go straight to Cochrane-Orcutt (or Prais-Winsten) */
	corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (!corc_model.errcode && corc_model.dfd == 0) {
	    corc_model.errcode = E_DF;
	}
	if ((err = corc_model.errcode)) {
	    clear_model(&corc_model);
	    return err;
	}
	rho0 = rho = corc_model.rho;
	pputs(prn, _("\nPerforming iterative calculation of rho...\n\n"));
    }

    pprintf(prn, "                 %s       RHO        %s\n",
	    _("ITER"), _("ESS"));

    iter = 0;
    diff = 1.0;
    while (diff > 0.001) {
	pprintf(prn, "          %10d %12.5f", ++iter, rho);
	clear_model(&corc_model);
	corc_model = lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
#ifdef AR_DEBUG
	fprintf(stderr, "corc_model: t1=%d, first two uhats: %g, %g\n",
		corc_model.t1, 
		corc_model.uhat[corc_model.t1],
		corc_model.uhat[corc_model.t1+1]);
#endif
	if ((err = corc_model.errcode)) {
	    clear_model(&corc_model);
	    return err;
	}
	pprintf(prn, "   %f\n", corc_model.ess);
#ifdef AR_DEBUG
	printmodel(&corc_model, pdinfo, OPT_NONE, prn);
#endif
	rho = autores(&corc_model, rho, (const double **) *pZ, opt);
	diff = (rho > rho0) ? rho - rho0 : rho0 - rho;
	rho0 = rho;
	if (iter == 20) break;
    }

    pprintf(prn, _("                final %11.5f\n\n"), rho);
    clear_model(&corc_model);

    *toprho = rho;

    return 0;
}

/* .......................................................... */

static int get_pos (const int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) return i;
    }

    return -1;
}

/* .......................................................... */

static int tsls_save_data (MODEL *pmod, const int *list, 
			   double **Z, DATAINFO *pdinfo)
{
    double **X = NULL;
    char *endog = NULL;
    int addvars = list[0];
    int i, j, k, pos, m, err = 0;

    pos = get_pos(pmod->list);
    m = pos - 2;

    X = malloc(list[0] * sizeof *X);
    
    endog = malloc(m * sizeof *endog);
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
    gretl_model_set_data(pmod, "tslsX", X, sizeof X);
    gretl_model_set_data(pmod, "endog", endog, m);

    return err;
}

const double *tsls_get_Xi (const MODEL *pmod, const double **Z, int i)
{
    const char *endog;
    double **X;
    const double *ret;

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

/**
 * tsls_func:
 * @list: dependent variable plus list of regressors.
 * @pos_in: position in the list for the separator between list
 *   of variables and list of instruments.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_R for robust VCV, OPT_S to save second-
 * stage regressors.
 *
 * Estimate the model given in @list by means of Two-Stage Least
 * Squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tsls_func (LIST list, int pos_in, double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt)
{
    int i, j, t, v, ncoeff;
    int *list1 = NULL, *list2 = NULL, *newlist = NULL;
    int *s1list = NULL, *s2list = NULL;
    int yno, n = pdinfo->n, orig_nvar = pdinfo->v;
    int nv, nxpx, pos, addvars = 0;
    MODEL tsls;
    double xx;
    double *yhat = NULL;

    if (pos_in > 0) {
	pos = pos_in;
    } else {
	pos = get_pos(list);
    }

    *gretl_errmsg = '\0';

    gretl_model_init(&tsls);

    list1 = malloc(pos * sizeof *list1);
    list2 = malloc((list[0] - pos + 1) * sizeof *list2);
    s1list = malloc((list[0] - pos + 2) * sizeof *s1list);
    s2list = malloc(pos * sizeof *s2list);
    newlist = malloc(pos * sizeof *newlist);

    if (list1 == NULL || list2 == NULL || s1list == NULL ||
	s2list == NULL || newlist == NULL) {
	tsls.errcode = E_ALLOC;
	goto tsls_bailout;
    }	

    list1[0] = pos - 1;
    for (i=1; i<pos; i++) {
	list1[i] = list[i];
    }
    tsls_omitzero(list1, *pZ, pdinfo->t1, pdinfo->t2);
    rearrange_list(list1);

    for (i=0; i<pos; i++) {
	s2list[i] = list1[i];
    }

    list2[0] = list[0] - pos;
    for (i=1; i<=list2[0]; i++) {
	list2[i] = list[i + pos];
    }
    tsls_omitzero(list2, *pZ, pdinfo->t1, pdinfo->t2);

    ncoeff = list2[0];
    if (ncoeff < list1[0] - 1) {
        sprintf(gretl_errmsg, 
		_("Order condition for identification is not satisfied.\n"
		"varlist 2 needs at least %d more variable(s) not in "
		"varlist1."), list1[0] - 1 - ncoeff);
	tsls.errcode = E_UNSPEC; 
	goto tsls_bailout;
    }

    /* now determine which fitted vals to obtain */
    if (tsls_match(list1, list2, newlist)) {
	strcpy(gretl_errmsg, 
	       _("Constant term is in varlist1 but not in varlist2"));
	tsls.errcode = E_UNSPEC;
	goto tsls_bailout;
    }

    /* newlist[0] holds the number of new vars to create */
    if (dataset_add_vars(newlist[0], pZ, pdinfo)) {
	tsls.errcode = E_ALLOC;
	goto tsls_bailout;
    } else {
	addvars = newlist[0];
    }

    /* deal with the variables for which instruments are needed */
    for (i=1; i<=newlist[0]; i++) { 
	yno = newlist[i];
        s1list[0] = ncoeff + 1;
        s1list[1] = yno;

        for (v=2; v<=s1list[0]; v++) {
	    s1list[v] = list2[v-1];
	}

	/* run first-stage regression */
	clear_model(&tsls);
	tsls = lsq(s1list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (tsls.errcode) {
	    goto tsls_bailout;
	}

        /* grab fitted values and stick into Z */
	for (j=2; j<=list1[0]; j++) {
	    if (list1[j] == newlist[i]) {
		s2list[j] = orig_nvar + i - 1;
		break;
	    }
	}

	for (t=0; t<n; t++) {
	    if (t >= tsls.t1 && t <= tsls.t2) {
		(*pZ)[orig_nvar+i-1][t] = tsls.yhat[t];
	    } else {
		(*pZ)[orig_nvar+i-1][t] = NADBL;
	    }
	}

	strcpy(pdinfo->varname[orig_nvar+i-1], pdinfo->varname[newlist[i]]);
    } 

    /* second-stage regression */
    clear_model(&tsls);
    tsls = lsq(s2list, pZ, pdinfo, OLS, OPT_NONE, 0.0);
    if (tsls.errcode) {
	goto tsls_bailout;
    }

    /* special: need to use the original RHS vars to compute residuals 
       and associated statistics */
    yhat = malloc(n * sizeof *yhat);
    if (yhat == NULL) {
	tsls.errcode = E_ALLOC;
	goto tsls_bailout;
    }

    tsls.ess = 0.0;
    for (t=tsls.t1; t<=tsls.t2; t++) {
	xx = 0.0;
	for (i=0; i<tsls.ncoeff; i++) {
	    xx += tsls.coeff[i] * (*pZ)[list1[i+2]][t];
	}
	yhat[t] = xx; 
	tsls.uhat[t] = (*pZ)[tsls.list[1]][t] - xx;
	tsls.ess += tsls.uhat[t] * tsls.uhat[t];
    }
    tsls.sigma = (tsls.ess >= 0.0) ? sqrt(tsls.ess / tsls.dfd) : 0.0;

    if (opt & OPT_R) {
	qr_tsls_vcv(&tsls, (const double **) *pZ, opt);
    } else {
	double *xpx = NULL, *xpy = NULL, *diag = NULL;

	nv = s2list[0] - 1;
	nxpx = nv * (nv + 1) / 2;
	xpx = malloc(nxpx * sizeof *xpx);
	xpy = malloc((s2list[0] + 1) * sizeof *xpy);
	diag = malloc((s2list[0] - 1) * sizeof *diag);

	if (xpy == NULL || xpx == NULL || diag == NULL) {
	    free(xpx);
	    free(xpy);
	    free(diag);
	    tsls.errcode = E_ALLOC;
	    goto tsls_bailout;
	}

	form_xpxxpy(s2list, tsls.t1, tsls.t2, *pZ, 0, 0.0, 0,
		    xpx, xpy);
	cholbeta(xpx, xpy, NULL, NULL, nv);    
	diaginv(xpx, xpy, diag, nv);
	for (i=0; i<tsls.ncoeff; i++) {
	    tsls.sderr[i] = tsls.sigma * sqrt(diag[i]); 
	}
	free(diag); 
	free(xpx);
	free(xpy);
    }

    tsls.rsq = corrrsq(tsls.nobs, &(*pZ)[tsls.list[1]][tsls.t1], 
		       yhat + tsls.t1);
    tsls.adjrsq = 
	1 - ((1 - tsls.rsq) * (tsls.nobs - 1) / tsls.dfd);
    tsls.fstt = tsls.rsq * tsls.dfd / (tsls.dfn * (1-tsls.rsq));
    gretl_aic_etc(&tsls);
    tsls.rho = rhohat(0, tsls.t1, tsls.t2, tsls.uhat);
    tsls.dw = dwstat(0, &tsls, *pZ);

    tsls.ci = TSLS;

    /* put the full list (possibly purged of zero elements) back in place */
    tsls.list = realloc(tsls.list, 
			(list1[0] + list2[0] + 2) * sizeof *tsls.list);
    if (tsls.list == NULL) {
	tsls.errcode = E_ALLOC;
    } else {
	tsls.list[0] = list1[0] + list2[0] + 1;
	j = 1;
	for (i=1; i<=list1[0]; i++) {
	    tsls.list[j++] = list1[i];
	}
	tsls.list[j++] = LISTSEP;
	for (i=1; i<=list2[0]; i++) {
	    tsls.list[j++] = list2[i];
	}
    }

    /* put the yhats into the model */
    for (t=tsls.t1; t<=tsls.t2; t++) {
	tsls.yhat[t] = yhat[t];
    }

 tsls_bailout:

    free(list1); 
    free(list2);
    free(s1list); 
    free(s2list);
    free(yhat); 

    if (opt & OPT_S) {
	tsls_save_data(&tsls, newlist, *pZ, pdinfo);
    } 
	
    dataset_drop_vars(addvars, pZ, pdinfo);
    free(newlist);

    if (tsls.errcode) {
	model_count_minus();
    }

    return tsls;
}

/* ........................................................ */

static int get_aux_uhat (MODEL *pmod, double *uhat1, double ***pZ, 
			 DATAINFO *pdinfo)
     /* feed in a uhat series -- this func finds the fitted values
	for the series using an aux. regression with squares, and
	adds that series to the data set */
{
    int *tmplist, *list, t, v = pdinfo->v;
    int i, j, l0 = pmod->list[0], check, shrink;
    MODEL aux;

    gretl_model_init(&aux);

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    /* add uhat1 to data set temporarily */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	(*pZ)[v][t] = uhat1[t];
    }

    /* tmplist is first used to construct a list of
       vars to be squared */
    tmplist = gretl_list_new((pmod->ifc)? l0 - 2 : l0 - 1);
    if (tmplist == NULL) return E_ALLOC;

    j = 1;
    for (i=2; i<=pmod->list[0]; i++) {
	if (pmod->list[i] != 0) {
	    tmplist[j++] = pmod->list[i];
	}
    }

    /* now generate the required squares */
    check = xpxgenr(tmplist, pZ, pdinfo, 0, 0);
    if (check < 1) {
	printf(_("generation of squares failed\n"));
	free(tmplist);
	return E_SQUARES;
    }

    tmplist = realloc(tmplist, (check + 2) * sizeof *tmplist);
    if (tmplist == NULL) return E_ALLOC;

    tmplist[0] = pdinfo->v - v - 1;
    for (i=1; i<=tmplist[0]; i++) {
	tmplist[i] = i + v;
    }

    list = gretl_list_add(pmod->list, tmplist, &check);
    if (check) {
	free(tmplist);
	return check;
    }

    list[1] = v;

    aux = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.);
    check = aux.errcode;
    if (check) shrink = pdinfo->v - v;
    else {
	for (t=aux.t1; t<=aux.t2; t++)
	    (*pZ)[v][t] = aux.yhat[t]; 
	shrink = pdinfo->v - v - 1;
    }

    if (shrink > 0) dataset_drop_vars(shrink, pZ, pdinfo);

    clear_model(&aux);

    free(tmplist);
    free(list);

    return check;
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

MODEL hsk_func (LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int err, lo, ncoeff, yno, t, nwt, v, n = pdinfo->n;
    int shrink, orig_nvar = pdinfo->v;
    double *uhat1, zz;
    int *hsklist;
    MODEL hsk;

    *gretl_errmsg = '\0';

    lo = list[0];         /* number of vars in original list */
    yno = list[1];        /* ID number of original dependent variable */
    ncoeff = list[0] - 1;
    rearrange_list(list);

    hsk = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (hsk.errcode) return hsk;

    uhat1 = malloc(n * sizeof *uhat1);
    if (uhat1 == NULL) {
	hsk.errcode = E_ALLOC;
	return hsk;
    }

    /* get residuals, square and log */
    for (t=hsk.t1; t<=hsk.t2; t++) {
	zz = hsk.uhat[t];
	uhat1[t] = log(zz * zz);
    }

    /* run auxiliary regression */
    err = get_aux_uhat(&hsk, uhat1, pZ, pdinfo);
    if (err) {
	hsk.errcode = err;
	free(uhat1);
	return hsk;
    }

    /* get fitted value from last regression and process */
    for (t=hsk.t1; t<=hsk.t2; t++) {
	zz = (*pZ)[pdinfo->v - 1][t];
	(*pZ)[pdinfo->v - 1][t] = 1.0 / sqrt(exp(zz));
    }    

    /* prepare to run weighted least squares */
    hsklist = malloc((lo + 2) * sizeof *hsklist);
    if (hsklist == NULL) {
	hsk.errcode = E_ALLOC;
	free(uhat1);
	return hsk;
    }

    /* "hsklist" will be one variable longer than the original
       regression list, because it includes a weight variable */
    hsklist[0] = lo + 1;

    /* the variable last added to the dataset will be the weight var */
    nwt = hsklist[1] = pdinfo->v - 1;

    /* put the original dependent variable in at position 2 */
    hsklist[2] = yno;

    /* put the original indep vars into the WLS list */
    for (v=lo+1; v>=3; v--) hsklist[v] = list[v-1];

    clear_model(&hsk);
    hsk = lsq(hsklist, pZ, pdinfo, WLS, OPT_NONE, 0.0);
    hsk.ci = HSK;

    shrink = pdinfo->v - orig_nvar;
    if (shrink > 0) dataset_drop_vars(shrink, pZ, pdinfo);

    free(hsklist);
    free(uhat1);

    return hsk;
}

static double **allocate_hccm_p (int lo, int n)
{
    int i;
    double **p = malloc(lo * sizeof *p);

    if (p == NULL) return NULL;

    for (i=0; i<lo; i++) {
	p[i] = malloc(n * sizeof **p);
	if (p[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) free(p[j]);
	    free(p);
	    p = NULL;
	}
    }

    return p;
}

/**
 * hccm_func:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using OLS, compute
 * heteroskedasticity-consistent covariance matrix using the
 * McKinnon-White procedure, and report standard errors using this matrix.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL hccm_func (LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int nobs, m3, lo, index, ncoeff, i, j, n, t, t1, t2;
    double xx, *st = NULL, *uhat1 = NULL, **p = NULL;
    MODEL hccm;

    *gretl_errmsg = '\0';

    gretl_model_init(&hccm);

    n = pdinfo->n;
    t1 = pdinfo->t1; t2 = pdinfo->t2;
    lo = list[0];
    m3 = 1 + (lo - 1) * (lo - 1);

    /* first try allocating memory */
    st = malloc(lo * sizeof *st);
    uhat1 = malloc(pdinfo->n * sizeof *uhat1);
    p = allocate_hccm_p(lo, t2 - t1 + 1);

    if (st == NULL || p == NULL || uhat1 == NULL) {
	hccm.errcode = E_ALLOC;
	goto bailout;
    }

    ncoeff = list[0] - 1;
    rearrange_list(list);

    /* run a regular OLS */
    hccm = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (hccm.errcode) {
	goto bailout;
    }

    hccm.ci = HCCM;
    nobs = hccm.nobs;

    if (get_use_qr()) {
	/* vcv is already computed */
	int nt = (ncoeff * ncoeff + ncoeff) / 2; 
	double sgmasq = hccm.sigma * hccm.sigma;

	for (i=0; i<nt; i++) {
	    hccm.vcv[i] /= sgmasq;
	}
    } else if (makevcv(&hccm)) {
	hccm.errcode = E_ALLOC;
	goto bailout;
    } 

    for (i=1; i<=ncoeff; i++) {
	for (t=t1; t<=t2; t++) {
	    xx = 0.0;
	    for (j=1; j<=ncoeff; j++) {
		if (i <= j) {
		    index = ijton(i-1, j-1, ncoeff);
		} else {
		    index = ijton(j-1, i-1, ncoeff);
		}
		xx += hccm.vcv[index] * (*pZ)[list[j+1]][t];
	    }
	    p[i][t] = xx;
	}
    }

    for (t=t1; t<=t2; t++) {
	xx = 0.0;
	for (i=1; i<=ncoeff; i++) {
	    xx += (*pZ)[list[i+1]][t] * p[i][t];
	}
	if (floateq(xx, 1.0)) xx = 0.0;
	uhat1[t] = hccm.uhat[t] / (1 - xx);
    }

    for (i=1; i<=ncoeff; i++) {
	xx = 0.0;
	for (t=t1; t<=t2; t++) {
	    xx += p[i][t] * uhat1[t];
	}
	st[i] = xx;
    }

    for (t=t1; t<=t2; t++) {
	for (i=1; i<=ncoeff; i++) {
	    p[i][t] *= uhat1[t];
	}
    }

    index = 0;
    for (i=1; i<=ncoeff; i++) {
	for (j=i; j<=ncoeff; j++) {
	    xx = 0.0;
	    for (t=t1; t<=t2; t++) {
		xx += p[i][t] * p[j][t];
	    }
	    xx = xx * (nobs - 1) / nobs -
		(nobs - 1) * st[i] * st[j] / (nobs * nobs);
	    if (i == j) {
		hccm.sderr[i-1] = sqrt(xx);
	    }
	    hccm.vcv[index++] = xx;
	}
    }

    /* substitute robust F stat */
    if (hccm.dfd > 0 && hccm.dfn > 1) {
	hccm.fstt = robust_omit_F(NULL, &hccm);
    }

 bailout:

    free(st);
    free(uhat1);

    if (p != NULL) {
	for (i=0; i<lo; i++) {
	    free(p[i]);
	}
	free(p);
    }

    set_model_id(&hccm);

    return hccm;
}

static int var_already_there (MODEL *pmod, double **Z, int testv)
{
    int i;
    int n = pmod->t2 - pmod->t1 + 1;
    int ret = 0;

    for (i=2; i<=pmod->list[0]; i++) {
	if (pmod->list[i] == 0) continue;
	if (vars_identical(Z[pmod->list[i]], Z[testv], n)) {
	    ret = 1;
	    break;
	}
    }

    return ret;
}

/**
 * whites_test:
 * @pmod: #MODEL struct.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @test: hypothesis test results struct.
 *
 * Runs White's test for heteroskedasticity on the given model,
 * putting the results into @test.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int whites_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		 PRN *prn, GRETLTEST *test)
{
    int lo, ncoeff, yno, i, t, check = 0;
    int shrink, v = pdinfo->v;
    int *tmplist = NULL, *list = NULL;
    double zz;
    MODEL white;
    int err = 0;

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == LOGISTIC) 
	return E_NOTIMP;

    if ((err = list_members_replaced(pmod->list, pdinfo, pmod->ID)))
	return err;

    gretl_model_init(&white);

    lo = pmod->list[0];
    yno = pmod->list[1];
    ncoeff = pmod->list[0] - 1;

    /* make space in data set */
    if (dataset_add_vars(1, pZ, pdinfo)) err = E_ALLOC;

    if (!err) {
	/* get residuals, square and add to data set */
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    zz = pmod->uhat[t];
	    (*pZ)[v][t] = zz * zz;
	}
	strcpy(pdinfo->varname[v], "uhatsq");

	tmplist = gretl_list_new((pmod->ifc)? lo - 2 : lo - 1);
	if (tmplist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=1; i<=tmplist[0]; i++) {
		tmplist[i] = pmod->list[i + 1 + pmod->ifc];
	    }
	}
    }

    if (!err) {
	/* now add squares */
	check = xpxgenr(tmplist, pZ, pdinfo, 0, 0);
	if (check < 1) {
	    fprintf(stderr, I_("generation of squares failed\n"));
	    free(tmplist);
	    err = E_SQUARES;
	}
    }

    if (!err) {
	tmplist = realloc(tmplist, (check + 2) * sizeof *tmplist);
	if (tmplist == NULL) err = E_ALLOC;
	else {
	    int k = 1;

	    tmplist[0] = pdinfo->v - v - 1; 
	    for (i=1; i<=tmplist[0]; i++) {
		/* check here for identical vars? */
		if (!var_already_there(pmod, *pZ, i + v)) { 
		    tmplist[k++] = i + v;
		}
	    }
	    tmplist[0] = k - 1;
	}
    }

    if (!err) {
	list = gretl_list_add(pmod->list, tmplist, &err);
	if (err) {
	    fprintf(stderr, I_("didn't add to list\n"));
	}
    }

    if (!err) {
	list[1] = v; 
	/* run auxiliary regression and print results */
	white = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.);
	err = white.errcode;
    }

    if (!err) {
	white.aux = AUX_WHITE;
	printmodel(&white, pdinfo, OPT_NONE, prn);

	if (test != NULL) {
	    gretl_test_init(test);
	    strcpy(test->type, N_("White's test for heteroskedasticity"));
	    strcpy(test->h_0, N_("heteroskedasticity not present"));
	    test->teststat = GRETL_TEST_TR2;
	    test->dfn = white.ncoeff - 1;
	    test->value = white.rsq * white.nobs;
	    test->pvalue = chisq(test->value, test->dfn);
	}
    }

    clear_model(&white);

    shrink = pdinfo->v - v;
    if (shrink > 0) dataset_drop_vars(shrink, pZ, pdinfo);

    free(tmplist);
    free(list);

    return err;
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

MODEL ar_func (LIST list, int pos, double ***pZ, 
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    double diff = 100.0, ess = 0, tss = 0, xx;
    int i, j, t, t1, t2, p, vc, yno, ryno = 0, iter = 0;
    int err, lag, n = pdinfo->n, v = pdinfo->v;
    int *arlist = NULL, *reglist = NULL, *reglist2 = NULL, *rholist = NULL;
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
	free(arlist);
	free(reglist);
	free(reglist2);
	free(rholist);
	ar.errcode = E_ALLOC;
	return ar;
    }
    
    arlist[0] = pos - 1;
    for (i=1; i<pos; i++) arlist[i] = list[i];

    reglist2[0] = reglist[0] = list[0] - pos;
    for (i=1; i<=reglist[0]; i++) reglist[i] = list[i + pos];

    rholist[0] = arlist[0] + 1;
    p = arlist[arlist[0]];
#if 0
    printf("arlist:\n"); printlist(arlist); 
    printf("reglist:\n"); printlist(reglist); 
#endif

    if (gretl_hasconst(reglist)) rearrange_list(reglist);

    /* special case: ar 1 ; ... => use CORC */
    if (arlist[0] == 1 && arlist[1] == 1) {
	err = hilu_corc(&xx, reglist, pZ, pdinfo, NULL, 1, CORC, prn);
	if (err) {
	    ar.errcode = err;
	} else {
	    ar = lsq(reglist, pZ, pdinfo, CORC, OPT_NONE, xx);
	}
	printmodel(&ar, pdinfo, opt, prn); 
	free(arlist);
	free(reglist);
	free(reglist2);
	free(rholist);
	return ar; 
    }

    /* make room for the uhat terms and transformed data */
    if (dataset_add_vars(arlist[0] + 1 + reglist[0], pZ, pdinfo)) {
	free(reglist);
	ar.errcode = E_ALLOC;
	return ar;
    }

    /*  rearrange_list(reglist); */ 
    yno = reglist[1];

    /* first pass: estimate model via OLS */
    ar = lsq(reglist, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (ar.errcode) {
	free(reglist);	
	return ar;
    }

    t1 = ar.t1; t2 = ar.t2;
    rholist[1] = v;

    pprintf(prn, "%s\n\n", _("Generalized Cochrane-Orcutt estimation"));
    bufspace(17, prn);
    /* xgettext:no-c-format */
    pputs(prn, _("ITER             ESS           % CHANGE"));
    pputs(prn, "\n\n");

    /* now loop while ess is changing */
    while (diff > 0.005) {
	iter++;
	for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
	/* special computation of uhat */
	for (t=t1; t<=t2; t++) {
	    xx = (*pZ)[yno][t];
	    for (j=0; j<reglist[0]-1; j++) {
		xx -= ar.coeff[j] * (*pZ)[reglist[j+2]][t];
	    }
	    (*pZ)[v][t] = xx;
	}
	for (i=1; i<=arlist[0]; i++) {
	    lag = arlist[i];
	    rholist[1+i] = v + i;
	    for (t=0; t<t1+lag; t++) (*pZ)[v+i][t] = NADBL;
	    for (t=t1+lag; t<=t2; t++)
		(*pZ)[v+i][t] = (*pZ)[v][t-lag];
	}
	ryno = vc = v + i;

	/* now estimate the rho terms */
	if (iter > 1) clear_model(&rhomod);
	rhomod = lsq(rholist, pZ, pdinfo, OLS, OPT_A, 0.0);

	/* and rho-transform the data */
	for (i=1; i<=reglist[0]; i++) {
	    for (t=0; t<n; t++) (*pZ)[vc][t] = NADBL;
	    for (t=t1+p; t<=t2; t++) {
		xx = (*pZ)[reglist[i]][t];
		for (j=1; j<=arlist[0]; j++) {
		    lag = arlist[j];
		    xx -= rhomod.coeff[j-1] * (*pZ)[reglist[i]][t-lag];
		}
		(*pZ)[vc][t] = xx;
	    }
	    reglist2[i] = vc;
	    vc++;
	}

	/* estimate the transformed model */
	clear_model(&ar);
	ar = lsq(reglist2, pZ, pdinfo, OLS, OPT_A, 0.0);

        if (iter > 1) {
	    diff = 100 * (ar.ess - ess)/ess;
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

	if (iter == 20) break;
    } /* end "ess changing" loop */

    for (i=0; i<=reglist[0]; i++) {
	ar.list[i] = reglist[i];
    }
    i = gretl_hasconst(reglist);
    if (i > 1) ar.ifc = 1;
    if (ar.ifc) ar.dfn -= 1;
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

    ar.rsq = corrrsq(ar.nobs, &(*pZ)[reglist[1]][ar.t1], ar.yhat + ar.t1);
    ar.adjrsq = 1 - ((1 - ar.rsq)*(ar.nobs - 1)/ar.dfd);

    /* special computation of TSS */
    xx = gretl_mean(ar.t1, ar.t2, (*pZ)[ryno]);
    for (t=ar.t1; t<=ar.t2; t++) {
	tss += ((*pZ)[ryno][t] - xx) * ((*pZ)[ryno][t] - xx);
    }
    ar.fstt = ar.dfd * (tss - ar.ess) / (ar.dfn * ar.ess);
    gretl_aic_etc(&ar);
    ar.dw = dwstat(p, &ar, *pZ);
    ar.rho = rhohat(p, ar.t1, ar.t2, ar.uhat);

    dataset_drop_vars(arlist[0] + 1 + reglist[0], pZ, pdinfo);
    free(reglist);
    free(reglist2);
    free(rholist);

    if (ar_info_init (&ar, 1 + p)) {
	ar.errcode = E_ALLOC;
    } else {
	for (i=0; i<=arlist[0]; i++) { 
	    ar.arinfo->arlist[i] = arlist[i];
	    if (i >= 1) {
		ar.arinfo->rho[i] = rhomod.coeff[i-1];
		ar.arinfo->sderr[i] = rhomod.sderr[i-1];
	    }
	}
    }

    set_model_id(&ar);
    printmodel(&ar, pdinfo, opt, prn);    

    free(arlist);
    clear_model(&rhomod);

    return ar;
}

/* ..........................................................  */

static void omitzero (MODEL *pmod, const DATAINFO *pdinfo, double **Z)
/* From 2 to end of list, omits variables with all zero observations
   and packs the rest of them */
{
    int t = 0, v, lv, offset, wtzero = 0, drop = 0;
    double xx = 0.;
    char vnamebit[20];

    offset = (pmod->ci == WLS)? 3 : 2;

    for (v=offset; v<=pmod->list[0]; v++) {
        lv = pmod->list[v];
        if (gretl_iszero(pmod->t1, pmod->t2, Z[lv])) {
	    list_exclude(v, pmod->list);
	    if (pdinfo->varname[lv][0] != 0) {
		sprintf(vnamebit, "%s ", pdinfo->varname[lv]);
	    } else {
		sprintf(vnamebit, "%s %d ", _("variable"), lv);
	    }
	    strcat(gretl_msg, vnamebit);
	    drop = 1;
	}
    }

    if (pmod->nwt) {
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
		list_exclude(v, pmod->list);
		sprintf(vnamebit, _("weighted %s "), pdinfo->varname[lv]);
		strcat(gretl_msg, vnamebit);
		drop = 1;
	    }
	}
    }

    if (drop) strcat(gretl_msg, _("omitted because all obs are zero."));
}

/* .........................................................   */

static void tsls_omitzero (int *list, double **Z, int t1, int t2)
{
    int i, v;

    for (i=2; i<=list[0]; i++) {
        v = list[i];
        if (gretl_iszero(t1, t2, Z[v])) {
	    list_exclude(i, list);
	    i--;
	}
    }
}

/* ...........................................................*/

static int zerror (int t1, int t2, int yno, int nwt, double ***pZ)
{
    double xx, yy;
    int t;

    xx = gretl_mean(t1, t2, (*pZ)[yno]);
    yy = gretl_stddev(t1, t2, (*pZ)[yno]);

    if (floateq(xx, 0.0) && floateq(yy, 0.0)) return 1;

    if (nwt) {
	xx = 0.0;
	for (t=t1; t<=t2; t++) {
	    xx = (*pZ)[nwt][t] * (*pZ)[yno][t];
	    if (floatneq(xx, 0.0)) return 0;
	}
	return 1;
    }

    return 0;
}

/* .......................................................... */

static int lagdepvar (const int *list, const DATAINFO *pdinfo, 
		      double ***pZ)
/* attempt to detect presence of a lagged dependent variable
   among the regressors -- if found, return the position of this
   lagged var in the list; otherwise return 0 */
{
    int i, t;
    char depvar[VNAMELEN], othervar[VNAMELEN];
    char *p;

    strcpy(depvar, pdinfo->varname[list[1]]);

    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) break;
	strcpy(othervar, pdinfo->varname[list[i]]);
	p = strrchr(othervar, '_');
	if (p != NULL && isdigit(*(p + 1))) {
	    /* looks like a lag */
	    size_t len = strlen(othervar) - strlen(p);

	    if (!strncmp(depvar, othervar, len)) {
		int gotlag = 1;

		/* strong candidate for lagged depvar, but make sure */
		for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) {
		    if ((*pZ)[list[1]][t-1] != (*pZ)[list[i]][t]) {
			gotlag = 0;
			break;
		    }
		}
		if (gotlag) return i;
	    }
	}
    } 

    return 0;
}

/* ............................................................ */

static int tsls_match (const int *list1, const int *list2, int *newlist)
/*
  Determines which variables in list1, when compared to all
  predetermined and exogenous variables in list2, need to have a
  reduced form ols regression run on them.  Returns the newlist of 
  dependent variables so that a reduced form ols regression
  can be run on each of them.  */
{
    int i, j, m, index = 0;
    int lo = list1[0], l2o = list2[0];

    for (i=2; i<=lo; i++) {     
	m = 0;
	for (j=1; j<=l2o; j++) {
	    if (list1[i] == list2[j]) j = l2o + 1;
	    else m++;
	    if (m == l2o) {
		if (list1[i] == 0) return 1;
		newlist[++index] = list1[i];
	    }
	}
    }  
    newlist[0] = index;

    return 0;
}

/* ............................................................  */

static double wt_dummy_mean (const MODEL *pmod, double **Z) 
/* returns mean of dependent variable in WLS model w. dummy weight */
{
    int m = pmod->nobs, yno = pmod->list[1];
    register int t;
    double sum = 0.0;

    if (m <= 0) return NADBL;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (floateq(Z[pmod->nwt][t], 0.0)) continue;
	else {
	    if (na(Z[yno][t])) {
		m--;
		continue;
	    } else {
		sum += Z[yno][t]; 
	    }
	}
    }

    return sum / m;
}

/* .............................................................  */

static double wt_dummy_stddev (const MODEL *pmod, double **Z) 
/*  returns standard deviation of dep. var. in WLS model
    with a dummy variable for weight.
*/
{
    int m = pmod->nobs, yno = pmod->list[1];
    register int t;
    double sumsq, xx, xbar;

    if (m == 0) return NADBL;
    xbar = wt_dummy_mean(pmod, Z);
    if (na(xbar)) return NADBL;
    sumsq = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
        xx = Z[yno][t] - xbar;
        if (floatneq(Z[pmod->nwt][t], 0.0) && !na(Z[yno][t]))
	    sumsq += xx * xx;
    }

    sumsq = (m > 1)? sumsq / (m-1) : 0.0;

    if (sumsq >= 0) return sqrt(sumsq);
    else return NADBL;
}

/**
 * arch:
 * @order: lag order for ARCH process.
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @test: hypothesis test results struct.
 * @opt: mat contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list via OLS, and test for Auto-
 * Regressive Conditional Heteroskedasticity.  If the latter is
 * significant, re-restimate the model using weighted least
 * squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch (int order, LIST list, double ***pZ, DATAINFO *pdinfo, 
	    GRETLTEST *test, gretlopt opt, PRN *prn)
{
    MODEL archmod;
    int *wlist = NULL, *arlist = NULL;
    int i, t, nwt, nv, n = pdinfo->n;
    double LM, xx;
    int err = 0;

    *gretl_errmsg = '\0';

    gretl_model_init(&archmod);

    /* assess the lag order */
    if (order < 1) {
	archmod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, _("Invalid lag order for arch (%d)"), order);
	err = 1;
    }

    if (!err) {
	/* allocate workspace */
	if (dataset_add_vars(order + 1, pZ, pdinfo) || 
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
	archmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	err = archmod.errcode;
    }

    if (!err) {
	nv = pdinfo->v - order - 1;
	strcpy(pdinfo->varname[nv], "utsq");
	for (t=0; t<n; t++) (*pZ)[nv][t] = NADBL;
	for (t=archmod.t1; t<=archmod.t2; t++) {
	    xx = archmod.uhat[t];
	    (*pZ)[nv][t] = xx * xx;
	}
	/* also lags of squared resids */
	for (i=1; i<=order; i++) {
	    nv =  pdinfo->v - order + i - 1;
	    arlist[i+2] = nv;
	    sprintf(pdinfo->varname[nv], "utsq_%d", i);
	    for (t=0; t<n; t++) (*pZ)[nv][t] = NADBL;
	    for (t=archmod.t1+i; t<=archmod.t2; t++) 
		(*pZ)[nv][t] = (*pZ)[arlist[1]][t-i];
	}

	/* run aux. regression */
	clear_model(&archmod);
	archmod = lsq(arlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	err = archmod.errcode;
    }

    if (!err) {
	/* print results */
	archmod.aux = AUX_ARCH;
	archmod.order = order;
	printmodel(&archmod, pdinfo, OPT_NONE, prn);
	pprintf(prn, _("No of obs. = %d, unadjusted R^2 = %f\n"),
		archmod.nobs, archmod.rsq);
	LM = archmod.nobs * archmod.rsq;
	xx = chisq(LM, order);

	if (test != NULL) {
	    gretl_test_init(test);
	    strcpy(test->type, N_("Test for ARCH of order %s"));
	    sprintf(test->param, "%d", order);
	    strcpy(test->h_0, N_("no ARCH effect is present"));
	    test->teststat = GRETL_TEST_TR2;
	    test->dfn = order;
	    test->value = LM;
	    test->pvalue = xx;
	}

	pprintf(prn, _("LM test statistic (%f) is distributed as Chi-square "
		"(%d)\nArea to the right of LM = %f  "), LM, order, xx);
	if (xx > 0.1) 
	    pprintf(prn, "\n%s.\n%s.\n",
		    _("ARCH effect is insignificant at the 10 percent level"),
		    _("Weighted estimation not done"));
	else {
	    pprintf(prn, "\n%s.\n",
		    _("ARCH effect is significant at the 10 percent level"));
	    /* weighted estimation */
	    wlist = malloc((list[0] + 2) * sizeof *wlist);
	    if (wlist == NULL) {
		archmod.errcode = E_ALLOC;
	    } else {
		wlist[0] = list[0] + 1;
		nwt = wlist[1] = pdinfo->v - 1; /* weight var */
		for (i=2; i<=wlist[0]; i++) wlist[i] = list[i-1];
		nv = pdinfo->v - order - 1;
		for (t=archmod.t1; t<=archmod.t2; t++) {
		    xx = archmod.yhat[t];
		    if (xx <= 0.0) xx = (*pZ)[nv][t];
		    (*pZ)[nwt][t] = 1/sqrt(xx);
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

    if (arlist != NULL) free(arlist);
    if (wlist != NULL) free(wlist);

    dataset_drop_vars(order + 1, pZ, pdinfo); 

    return archmod;
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

MODEL lad (LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int err = 0;
    MODEL lad_model;
    void *handle;
    int (*lad_driver)(MODEL *, double **, DATAINFO *);

    /* run an initial OLS to "set the model up" and check for errors.
       the lad function will overwrite the coefficients etc.
    */
    lad_model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);

    if ((err = lad_model.errcode)) {
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
 * @list: dependent variable plus AR and MA orders
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @PRN: for printing details of iterations (or NULL) 
 *
 * Calculate ARMA estimates.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arma (int *list, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    MODEL armod;
    void *handle;
    MODEL (*arma_model) (int *, const double **, DATAINFO *, PRN *);

    *gretl_errmsg = '\0';

    arma_model = get_plugin_function("arma_model", &handle);

    if (arma_model == NULL) {
	gretl_model_init(&armod);
	armod.errcode = E_FOPEN;
	return armod;
    }

    armod = (*arma_model) (list, Z, pdinfo, prn);

    close_plugin(handle);

    set_model_id(&armod);

    return armod;
} 

#ifdef HAVE_X12A

/**
 * arma_x12a:
 * @list: dependent variable plus AR and MA orders
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @PRN: for printing details of iterations (or NULL).
 * @ppaths: gretl path info struct (so we can find X-12-ARIMA)
 *
 * Calculate ARMA estimates, via a call to X-12-ARIMA.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arma_x12 (int *list, const double **Z, DATAINFO *pdinfo, PRN *prn,
		const PATHS *ppaths)
{
    MODEL armod;
    void *handle;
    MODEL (*arma_x12_model) (int *, const double **, DATAINFO *, PRN *, 
			     const char *, const char *, int);
    int gui = gretl_using_gui(ppaths);

    *gretl_errmsg = '\0';

    arma_x12_model = get_plugin_function("arma_x12_model", &handle);
    if (arma_x12_model == NULL) {
	gretl_model_init(&armod);
	armod.errcode = E_FOPEN;
	return armod;
    }

    armod = (*arma_x12_model) (list, Z, pdinfo, prn, ppaths->x12a, 
			       ppaths->x12adir, gui);

    close_plugin(handle);

    set_model_id(&armod);

    return armod;
}  

#endif /* HAVE_X12A */

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

MODEL logistic_model (int *list, double ***pZ, DATAINFO *pdinfo,
		      const char *param)
{
    MODEL lmod;
    void *handle;
    MODEL (*logistic_estimate) (int *, double ***, DATAINFO *, const char *);

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

MODEL tobit_model (LIST list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL tmod;
    void *handle;
    MODEL (* tobit_estimate) (int *, double ***, DATAINFO *, PRN *);

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

/**
 * garch:
 * @list: dependent variable plus arch and garch orders
 * @pZ: pointer to data matrix
 * @pdinfo: information on the data set
 * @opt: can specify robust standard errors and VCV
 * @prn: for printing details of iterations (or NULL)
 *
 * Calculate GARCH estimates.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL garch (int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn)
{
    MODEL gmod;
    void *handle;
    PRN *myprn;
    MODEL (*garch_model) (int *, double ***, DATAINFO *, PRN *,
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
