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
#include "internal.h"

/* There's a balancing act with 'TINY' here.  It's the minimum value
   for determinant that libgretl will accept before rejecting a
   data matrix as too highly collinear.  If you set it too high,
   data sets for which gretl could produce reasonable estimates will
   be rejected.  If you set it too low (and even 100 * DBL_EPSILON
   is definitely too low), gretl will produce more or less worthless
   coefficient estimates when given highly collinear data.  If you're
   tempted to change the value of TINY, check how gretl does on the
   NIST reference data sets for linear regression and ensure you're
   not getting any garbage results.
*/

#define TINY 0.5e-8 /* was 1.0e-13, produced poor results on NIST Filip
		       test */
#define SMALL 1.0e-8

typedef struct {
    double *xpx;
    double *xpy;
    int ivalue;
    int nv;
    int errcode;
} XPXXPY;	

typedef struct {
    XPXXPY xpxxpy;
    double *coeff;
    double rss;
    int errcode;
} CHOLBETA;

extern void _print_rho (int *arlist, const MODEL *pmod, 
			const int c, PRN *prn);
extern int _addtolist (const int *oldlist, const int *addvars, 
		       int **pnewlist, const DATAINFO *pdinfo, 
		       const int model_count);

/* private function prototypes */
static XPXXPY xpxxpy_func (const int *list, int t1, int t2, 
			   double **Z, int nwt, double rho);
static void regress (MODEL *pmod, XPXXPY xpxxpy, double **Z, 
		     int n, const double rho);
static CHOLBETA cholbeta (XPXXPY xpxxpy);
static void diaginv (XPXXPY xpxxpy, double *diag);
static double dwstat (int order, MODEL *pmod, double **Z);
static double corrrsq (int nobs, const double *y, const double *yhat);
static double rhohat (int order, int t1, int t2, const double *uhat);
static double altrho (int order, int t1, int t2, const double *uhat);
static int hatvar (MODEL *pmod, double **Z);
static void dropwt (int *list);
static void autores (int i, double **Z, const MODEL *pmod, double *uhat);
static int get_aux_uhat (MODEL *pmod, double *uhat1, double ***pZ, 
			 DATAINFO *pdinfo);
static void omitzero (MODEL *pmod, const DATAINFO *pdinfo, double **Z);
static void tsls_omitzero (int *list, double **Z, int t1, int t2);
static void rearrange (int *list);
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

/**
 * lsq:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: command index (see commands.h)
 * @opt: option flag: If = 1, then residuals, dw stat and rhohat are obtained.
 * @rho: coefficient for rho-differencing the data.
 *
 * Computes least squares estimates of the model specified by @list,
 * using an estimator determined by the value of @ci.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lsq (LIST list, double ***pZ, DATAINFO *pdinfo, 
	   const int ci, const int opt, const double rho)
{
    int l0, ifc, nwt, yno, i, n;
    int t, t1, t2, v, order, effobs = 0;
    int missv = 0, misst = 0;
    double xx;
    XPXXPY xpxxpy;
    MODEL model;

    if (list == NULL || pZ == NULL || pdinfo == NULL) {
	model.errcode = E_DATA;
        return model;
    }

    if (ci == HSK)
	return hsk_func(list, pZ, pdinfo);
    if (ci == HCCM)
	return hccm_func(list, pZ, pdinfo);

    _init_model(&model, pdinfo);

    if (list[0] == 1 || pdinfo->v == 1) {
	model.errcode = E_DATA;
        return model;
    }

    /* preserve a copy of the list supplied, for future reference */
    copylist(&(model.list), list);
    if (model.list == NULL) {
        model.errcode = E_ALLOC;
        return model;
    }

    n = pdinfo->n;
    model.t1 = pdinfo->t1;
    model.t2 = pdinfo->t2;
    model.ci = ci;

    /* Doing weighted least squares? */
    model.wt_dummy = 0;
    model.nwt = nwt = 0;
    if (ci == WLS) { 
	model.nwt = nwt = model.list[1];
	if (_iszero(model.t1, model.t2, (*pZ)[nwt])) {
	    model.errcode = E_WTZERO;
	    return model;
	}
	effobs = isdummy(nwt, model.t1, model.t2, *pZ);
	if (effobs) model.wt_dummy = 1;
    }

    /* check for missing obs in sample */
    if ((missv = _adjust_t1t2(&model, model.list, &model.t1, &model.t2, 
			      *pZ, &misst))) {
	if (!dated_daily_data(pdinfo)) {
	    sprintf(gretl_errmsg, _("Missing value encountered for "
		    "variable %d, obs %d"), missv, misst);
	    model.errcode = E_DATA;
	    return model;
	} else {
	    /* with daily data, try eliminating the missing obs? */
	    int misscount;
	    char *missvec = missobs_vector(*pZ, pdinfo, &misscount);
	    MISSOBS *mobs = NULL;

	    if (missvec == NULL || (mobs = malloc(sizeof *mobs)) == NULL) {
		model.errcode = E_ALLOC;
		return model;
	    } else {
		repack_missing(*pZ, pdinfo, missvec, misscount);
		model.t2 -= misscount;
		mobs->misscount = misscount;
		mobs->missvec = missvec;
		model.data = mobs;
	    }
	}
    }    
    t1 = model.t1; 
    t2 = model.t2; 

    if (ci == WLS) dropwt(model.list);
    yno = model.list[1];
    
    /* check for availability of data */
    if (t1 < 0 || t2 > n - 1) {
        model.errcode = E_NODATA;
        return model;
    }                   
    for (i=1; i<=model.list[0]; i++) {
        if (model.list[i] > pdinfo->v - 1) {
            model.errcode = E_UNKVAR;
            return model;
        }
    }       

    /* check for zero dependent var */
    if (zerror(t1, t2, yno, nwt, pZ)) {  
        model.errcode = E_ZERO;
        return model; 
    } 

    /* drop any vars that are all zero and repack the list */
    omitzero(&model, pdinfo, *pZ);

    /* see if the regressor list contains a constant (ID 0) */
    model.ifc = ifc = _hasconst(model.list);
    /* if so, move it to the last place */
    if (ifc) rearrange(model.list);

    /* check for presence of lagged dependent variable */
    model.ldepvar = lagdepvar(model.list, pdinfo, pZ);

    l0 = model.list[0];  /* holds 1 + number of coeffs */
    model.ncoeff = l0 - 1; 
    model.nobs = t2 - t1 + 1;
    if (effobs) model.nobs = effobs;

    /* check degrees of freedom */
    if (model.nobs < model.ncoeff) { 
	model.errcode = E_DF;
        sprintf(gretl_errmsg, _("No. of obs (%d) is less than no. "
		"of parameters (%d)"), model.nobs, model.ncoeff);
        return model; 
    }

    /* calculate regression results */
    xpxxpy = xpxxpy_func(model.list, t1, t2, *pZ, nwt, rho);
    model.tss = xpxxpy.xpy[l0];

    regress(&model, xpxxpy, *pZ, n, rho);
    free(xpxxpy.xpy);
    if (model.errcode) return model;

    /* get the mean and sd of depvar and make available */
    if (model.ci == WLS && model.wt_dummy) {
	model.ybar = wt_dummy_mean(&model, *pZ);
	model.sdy = wt_dummy_stddev(&model, *pZ);
    } else {
	model.ybar = _esl_mean(t1, t2, (*pZ)[yno]);
	model.sdy = _esl_stddev(t1, t2, (*pZ)[yno]);
    }

    /* Doing an autoregressive procedure? */
    if (ci == CORC || ci == HILU) {
	model.arlist = malloc(2 * sizeof(int));
	model.rhot = malloc(2 * sizeof(double));
	if (model.arlist == NULL || model.rhot == NULL) {
	    model.errcode = E_ALLOC;
	    return model;
	}
	model.arlist[0] = model.arlist[1] = 1;
	model.rhot[1] = model.rho_in = rho;
	if (model.ifc) {
	    model.coeff[model.ncoeff] /= (1.0 - rho);
	    model.sderr[model.ncoeff] /= (1.0 - rho);
	}
	model.uhat[t1] = NADBL;
	model.yhat[t1] = NADBL;
	for (t=t1+1; t<=t2; t++) {
	    xx = (*pZ)[yno][t] - rho * (*pZ)[yno][t-1];
	    for (v=1; v<=model.ncoeff-model.ifc; v++)
		xx -= model.coeff[v] * 
		    ((*pZ)[model.list[v+1]][t] - 
		    rho * (*pZ)[model.list[v+1]][t-1]);
	    if (model.ifc) xx -= (1 - rho) * model.coeff[model.ncoeff];
	    model.uhat[t] = xx;
	    model.yhat[t] = (*pZ)[yno][t] - xx;
	}
	model.rsq = 
	    corrrsq(t2-t1, &(*pZ)[yno][t1+1], model.yhat + t1+1);
    	model.adjrsq = 
           1 - ((1 - model.rsq)*(t2 - t1 - 1)/model.dfd);
    }

    /* weighted least squares: fix fitted values, ESS, sigma */
    if (ci == WLS && !(model.wt_dummy)) {
	model.ess_wt = model.ess;
	model.sigma_wt = model.sigma;
	model.ess = 0.0;
	for (t=t1; t<=t2; t++) {
	    if ((*pZ)[nwt][t] == 0.0) {
		model.yhat[t] = model.uhat[t] = NADBL;
		model.nobs -= 1;
	    } else {
		model.yhat[t] /= (*pZ)[nwt][t];
		xx = model.uhat[t] /= (*pZ)[nwt][t];
		model.ess += xx * xx;
	    }
	}
	model.sigma = sqrt(model.ess/model.dfd);
    }
    if (ci == WLS && model.wt_dummy) {
	for (t=t1; t<=t2; t++) {
	    if (floateq((*pZ)[nwt][t], 0.0)) 
		model.yhat[t] = model.uhat[t] = NADBL;
	}
    }

    /* if opt = 1, compute residuals and rhohat */
    if (opt) {
	order = (ci == CORC || ci == HILU)? 1 : 0;
	model.rho = rhohat(order, t1, t2, model.uhat);
	model.dw = dwstat(order, &model, *pZ);
    }

    /* Generate model selection statistics */
    _aicetc(&model);

    /* If we eliminated any missing observations, restore
       them now */
    if (model.data != NULL) {
	MISSOBS *mobs = (MISSOBS *) model.data;

	undo_repack_missing(*pZ, pdinfo, mobs->missvec,
			    mobs->misscount);
	reorganize_uhat_yhat(&model);
    }

    return model;
}

/* .......................................................... */

static XPXXPY xpxxpy_func (const int *list, int t1, int t2, 
			   double **Z, int nwt, double rho)
/*
        This function forms the X'X matrix and X'y vector
        - if rho is non-zero transforms data first
        - if nwt is non-zero, uses that variable as weight
        Z[v][t] = t-th observation for the v-th variable
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
    int i, j, li, lj, m, l0 = list[0], yno = list[1], t;
    double xx, z1, z2;
    XPXXPY xpxxpy;

    i = l0 - 1;
    m = i * (i + 1) / 2;

    if ((xpxxpy.xpy = malloc((l0 + 1) * sizeof(double))) == NULL ||
	(xpxxpy.xpx = malloc((m + 1) * sizeof(double))) == NULL) {
        xpxxpy.errcode = E_ALLOC;
        return xpxxpy;
    }
    xpxxpy.xpy[0] = xpxxpy.xpy[l0] = 0;

    xpxxpy.nv = l0 - 1;
    if (rho) for (t=t1+1; t<=t2; t++) {
        xx = Z[yno][t] - rho * Z[yno][t-1];
        xpxxpy.xpy[0] += xx;
        xpxxpy.xpy[l0] += xx * xx;
    }
    else if (nwt) for (t=t1; t<=t2; t++) {
        xx = Z[yno][t] * Z[nwt][t];       
        xpxxpy.xpy[0] += xx;
        xpxxpy.xpy[l0] += xx * xx;
    }
    else for (t=t1; t<=t2; t++) {
        xx = Z[yno][t]; 
        xpxxpy.xpy[0] += xx;
        xpxxpy.xpy[l0] += xx * xx;
    }
    if (floateq(xpxxpy.xpy[l0], 0.0)) {
         xpxxpy.ivalue = yno; 
         return xpxxpy; 
    }    
    m = 0;

    if (rho) for (i=2; i<=l0; i++) {
        li = list[i];
        z1 = li? rho: 0.0;
        for (j=i; j<=l0; j++) {
            lj = list[j];
            z2 = lj? rho: 0.0;
            xx = 0.0;
            for (t=t1+1; t<=t2; t++)
                xx += (Z[li][t] - z1 * Z[li][t-1]) * 
		    (Z[lj][t] - z2 * Z[lj][t-1]);
                if (floateq(xx, 0.0) && li == lj)  {
                    xpxxpy.ivalue = li;
                    return xpxxpy; 
                }
                xpxxpy.xpx[++m] = xx;
        }
        xx = 0;
        for (t=t1+1; t<=t2; t++)
            xx = xx + (Z[yno][t] - rho * Z[yno][t-1]) *
		(Z[li][t] - z1 * Z[li][t-1]);
        xpxxpy.xpy[i-1] = xx;
    }
    else if (nwt) for (i=2; i<=l0; i++) {
        li = list[i];
        for (j=i; j<=l0; j++) {
            lj = list[j];
            xx = 0.0;
            for (t=t1; t<=t2; t++) {
                z1 = Z[nwt][t];
                xx += z1 * z1 * Z[li][t] * Z[lj][t];
            }
            if (floateq(xx, 0.0) && li == lj)  {
                xpxxpy.ivalue = li;
                return xpxxpy;
            }   
            xpxxpy.xpx[++m] = xx;
        }
        xx = 0;
        for(t=t1; t<=t2; t++) {
            z1 = Z[nwt][t];
            xx += z1 * z1 * Z[yno][t] * Z[li][t];
        }
        xpxxpy.xpy[i-1] = xx;
    }
    else for (i=2; i<=l0; i++) {
        li = list[i];
        for (j=i; j<=l0; j++) {
            lj = list[j];
            xx = 0.0;
            for (t=t1; t<=t2; t++) xx += Z[li][t] * Z[lj][t];
            if (floateq(xx, 0.0) && li == lj)  {
                xpxxpy.ivalue = li;
                return xpxxpy;  
            }
            xpxxpy.xpx[++m] = xx;
        }
        xx = 0;
        for (t=t1; t<=t2; t++) xx += Z[yno][t] * Z[li][t];
        xpxxpy.xpy[i-1] = xx;
    }
    xpxxpy.ivalue = 0;
/*      for (i=1; i<=m; i++)  */
/*  	printf("xpx[%d] = %10.18f\n", i, xpxxpy.xpx[i]); */
    return xpxxpy; 
}

/* .......................................................... */

static void regress (MODEL *pmod, XPXXPY xpxxpy, double **Z, 
		     int n, const double rho)
/*
        This function takes xpx, the X'X matrix ouput
        by xpxxpy_func(), and xpy, which is X'y, and
        computes ols estimates and associated statistics.

        n = no. of observations per series in data set
        ifc = 1 if constant is present else = 0

        ess = error sum of squares
        sigma = standard error of regression
        fstt = f-statistics
        coeff = vector of regression coefficients
        sderr = vector of standard errors of regression coefficients
*/
{
    int t, v, nobs, nv, yno, nwt = pmod->nwt;
    int t1 = pmod->t1, t2 = pmod->t2;
    double *diag, ysum, ypy, zz, ess, rss, tss;
    double den = 0.0, sgmasq = 0.0;
    CHOLBETA cb;

    nv = xpxxpy.nv;
    yno = pmod->list[1];

    if ((pmod->sderr = calloc(nv+1, sizeof(double))) == NULL ||
	(pmod->yhat = calloc(n, sizeof(double))) == NULL ||
	(pmod->uhat = calloc(n, sizeof(double))) == NULL) {
        pmod->errcode = E_ALLOC;
        return;
    }
    nobs = pmod->nobs;
    if (rho) pmod->nobs = nobs = t2-t1;
    pmod->ncoeff = nv;
    pmod->dfd = nobs - nv;
    if (pmod->dfd < 0) { 
       pmod->errcode = E_DF; 
       return; 
    }
    pmod->dfn = nv - pmod->ifc;
    ysum = xpxxpy.xpy[0];
    ypy = xpxxpy.xpy[nv+1];
    if (floateq(ypy, 0.0)) { 
        pmod->errcode = E_YPY;
        return; 
    }
    zz = ysum * ysum/nobs;
    tss = ypy - zz;
    if (floatlt(tss, 0.0)) { 
        pmod->errcode = E_TSS; 
        return; 
    }

    /*  Choleski-decompose X'X and find the coefficients */
    cb = cholbeta(xpxxpy);
    pmod->coeff = cb.coeff;
    pmod->xpx = cb.xpxxpy.xpx;

    if (cb.errcode) {
        pmod->errcode = E_ALLOC;
        return;
    }   
    rss = cb.rss;
    if (rss == -1.0) { 
        pmod->errcode = E_SINGULAR;
        return; 
    }

    pmod->ess = ess = ypy - rss;
    if (ess < SMALL && ess > (-SMALL)) pmod->ess = ess = 0.0;
    else if (ess < 0.0) { 
        /*  pmod->errcode = E_ESS; */ 
	sprintf(gretl_errmsg, _("Error sum of squares (%g) is not > 0"),
		ess);
        return; 
    }

    if (pmod->dfd == 0) {
	pmod->sigma = 0.0;
	pmod->adjrsq = NADBL;
    } else {
	sgmasq = ess/pmod->dfd;
	pmod->sigma = sqrt(sgmasq);
	den = tss * pmod->dfd;
    }

    if (floatlt(tss, 0.0) || floateq(tss, 0.0)) {
       pmod->rsq = pmod->adjrsq = NADBL;
       pmod->errcode = E_TSS;
       return;
    }       

    hatvar(pmod, Z); 
    if (pmod->errcode) return;

    pmod->rsq = 1 - (ess/tss);
    if (pmod->dfd > 0) {
	pmod->adjrsq = 1 - (ess * (nobs-1)/den);
	if (!pmod->ifc) {  
	    pmod->rsq = corrrsq(nobs, &Z[yno][t1], pmod->yhat + t1);
	    pmod->adjrsq = 
		1 - ((1 - pmod->rsq)*(nobs - 1)/pmod->dfd);
	}
    }

    if (pmod->ifc && nv == 1) {
        zz = 0.0;
        pmod->dfn = 1;
    }
    if (sgmasq <= 0.0 || pmod->dfd == 0) pmod->fstt = NADBL;
    else pmod->fstt = (rss - zz * pmod->ifc)/(sgmasq * pmod->dfn);

    /* special treatment for weighted least squares (Ramanathan) */
#ifdef RAMANATHAN
    if (nwt && !(pmod->wt_dummy)) {
	int wobs = nobs;
	double *altyhat;

	altyhat = malloc(n * sizeof(double));
	zz = 0.0;
	for (t=pmod->t1; t<=pmod->t2; t++) { 
	    if (Z[nwt][t] == 0) {
		altyhat[t] = NADBL;
		wobs--;
		pmod->dfd -= 1;
	    } else {
		zz += pmod->yhat[t] * pmod->yhat[t];
		altyhat[t] = pmod->yhat[t] / Z[nwt][t];
	    }
	}
	pmod->dfn += 1;
	pmod->fstt = (zz * pmod->dfd)/(pmod->dfn * ess);
	pmod->rsq = corrrsq(nobs, &Z[yno][pmod->t1], altyhat + pmod->t1);
	pmod->adjrsq = 
           1 - ((1 - pmod->rsq)*(wobs - 1)/pmod->dfd);
	free(altyhat);
    }
#else
    /* alternative calculation of WLS stats, in agreement with GNU R */
    if (nwt && !(pmod->wt_dummy)) {
	int wobs = nobs;
	double w2, wmean = 0.0, wsum = 0.0;

	for (t=pmod->t1; t<=pmod->t2; t++) { 
	    if (Z[nwt][t] == 0) {
		wobs--;
		pmod->dfd -= 1;
	    } else {
		w2 = Z[nwt][t] * Z[nwt][t];
		wmean += w2 * Z[yno][t];
		wsum += w2;
	    }
	}
	wmean /= wsum;
	zz = 0.0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (Z[nwt][t] == 0) continue;
	    w2 = Z[nwt][t] * Z[nwt][t]; 
	    zz += w2 * (Z[yno][t] - wmean) * (Z[yno][t] - wmean);
	}
	pmod->fstt = ((zz - ess) * pmod->dfd)/(pmod->dfn * ess);
	pmod->rsq = (1 - (ess/zz));
	pmod->adjrsq = 1 - ((1 - pmod->rsq)*(nobs - 1)/pmod->dfd);
    }
#endif /* RAMANATHAN */

    diag = malloc((nv+1) * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }
    diaginv(xpxxpy, diag);
    for (v=1; v <= nv; v++) { 
       pmod->sderr[v] = pmod->sigma * sqrt(diag[v]); 
    }
    free(diag); 
    
    return;  
}

/* .......................................................... */

static CHOLBETA cholbeta (XPXXPY xpxxpy)
/*
  This function does an inplace Choleski decomposition of xpx
  (lower triangular matrix stacked in columns) and then
  solves the normal equations for coeff.  xpx is the X'X
  on input and Choleski decomposition on output; xpy is
  the X'y vector on input and chol. transformed t vector
  on output. coeff is the vector of estimated coefficients; 
  nv is the number of regression coefficients including the 
  constant.  */
{
    int nm1, i, j, k, kk, l, jm1, nv;
    double e, d, d1, test, xx;
    CHOLBETA cb;

    nv = xpxxpy.nv; 
    cb.errcode = 0; 

    if ((cb.coeff = malloc((nv + 1) * sizeof(double))) == NULL) {
        cb.errcode = E_ALLOC;
        return cb;
    }
    for (j=0; j<nv+1; j++) cb.coeff[j] = 0.0;

    cb.xpxxpy = xpxxpy;

    nm1 = nv - 1;
    e = 1/sqrt(xpxxpy.xpx[1]);
    xpxxpy.xpx[1] = e;
    xpxxpy.xpy[1] *= e;
    for (i=2; i<=nv; i++) xpxxpy.xpx[i] *= e;
    kk = nv + 1;

    for (j=2; j<=nv; j++) {
    /* diagonal elements */
        d = d1 = 0.0;
        k = j;
        jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
            xx = xpxxpy.xpx[k];
            d1 += xx * xpxxpy.xpy[l];
            d += xx * xx;
            k += nv-l;
        }
        test = xpxxpy.xpx[kk] - d;
        if (test <= TINY) {
           cb.rss = -1.0; 
           return cb;
        }
	if (test < SMALL) {
	    strcpy(gretl_msg, _("Warning: data matrix close to singularity!"));
	}
        e = 1/sqrt(test);
        xpxxpy.xpx[kk] = e;
        xpxxpy.xpy[j] = (xpxxpy.xpy[j] - d1) * e;
        /* off-diagonal elements */
        for (i=j+1; i<=nv; i++) {
            kk++;
            d = 0.0;
            k = j;
            for (l=1; l<=jm1; l++) {
                d += xpxxpy.xpx[k] * xpxxpy.xpx[k-j+i];
                k += nv - l;
            }
            xpxxpy.xpx[kk] = (xpxxpy.xpx[kk] - d) * e;
        }
        kk++;
    }
    kk--;
    d = 0.0;
    for(j=1; j<=nv; j++) {
        d += xpxxpy.xpy[j] * xpxxpy.xpy[j];
    }
    cb.rss = d;
    cb.coeff[nv] = xpxxpy.xpy[nv] * xpxxpy.xpx[kk];
    for(j=nm1; j>=1; j--) {
        d = xpxxpy.xpy[j];
        for (i=nv; i>=j+1; i--) {
            kk--;
            d = d - cb.coeff[i] * xpxxpy.xpx[kk];
        }
        kk--;
        cb.coeff[j] = d * xpxxpy.xpx[kk];
    }    
    return cb; 
}

/* ...............................................................    */

static void diaginv (XPXXPY xpxxpy, double *diag)
/*
        Solves for the diagonal elements of the X'X inverse matrix.

        xpx = Choleski-decomposed X'X matrix (input)
        xpy = X'y vector (input) used as work array
        diag = diagonal elements of X'X (output)
        nv = no. of regression coefficients
*/
{
    int kk = 1, l, m, nstop, k, i, j, nv;
    double d, e;

    nv = xpxxpy.nv;
    nstop = nv * (nv+1)/2;
    for (l=1; l<=nv-1; l++) {
        d = xpxxpy.xpx[kk];
        xpxxpy.xpy[l] = d;
        e = d * d;
        m = 0;
        if (l > 1) 
	    for (j=1; j<=l-1; j++) m += nv - j;
        for (i=l+1; i<=nv; i++) {
            d = 0.0;
            k = i + m;
            for (j=l; j<=i-1; j++) {
                d += xpxxpy.xpy[j] * xpxxpy.xpx[k];
                k += nv - j;
            }
            d = (-1.0) * d * xpxxpy.xpx[k];
            xpxxpy.xpy[i] = d;
            e += d * d;
        }
        kk += nv + 1 - l;
        diag[l] = e;
    }
    diag[nv] = xpxxpy.xpx[nstop] * xpxxpy.xpx[nstop];
}

/**
 * makevcv:
 * @pmod: gretl MODEL.
 *
 * Inverts the Choleski-decomposed stacked vector xpx and
 * computes the coefficient variance-covariance matrix.
 * 
 * Returns: 0 on successful completion, error code on error
 */

int makevcv (MODEL *pmod)
{
    int nv, dec, nm1, mst, kk, i, j, kj, icnt, m, k, l = 0;
    int idxpx;
    double sigma, d;

    nv = pmod->ncoeff;
    nm1 = nv - 1;
    mst = kk = idxpx = (nv * nv + nv)/2;

    pmod->vcv = malloc((mst + 1) * sizeof(double));
    if (pmod->vcv == NULL) return E_ALLOC;

    for (i=0; i<=nm1; i++) {
	mst -= i;
	/* find diagonal element */
	d = pmod->xpx[kk];
	if (i > 0) {
	    for (j=kk+1; j<=kk+i; j++) 
		d = d - pmod->xpx[j] * pmod->vcv[j];
	}
	pmod->vcv[kk] = d * pmod->xpx[kk];
	/* find off-diagonal elements indexed by kj */
	kj = kk;
	kk = kk - i - 2;
	if (i > nv-2) continue;
	for (j=i+1; j<=nm1; j++) {
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
		d += pmod->vcv[m] * pmod->xpx[l];
	    }
	    pmod->vcv[kj] = (-1.0) * d * pmod->xpx[l-1];
	}
    }

    if (pmod->ci == CUSUM) return 0;
    if ((pmod->ci == WLS && !(pmod->wt_dummy)) || 
	pmod->ci == ARCH || pmod->ci == HSK) sigma = pmod->sigma_wt;
    else sigma = pmod->sigma; 
    /* some estimators need special treatment */
    if (pmod->ci != HCCM && pmod->ci != LOGIT && pmod->ci != PROBIT) 
	for (k=0; k<idxpx; k++) pmod->vcv[k] = 
				    pmod->vcv[k+1] * sigma * sigma;
    if (pmod->ci == LOGIT || pmod->ci == PROBIT) 
	for (k=0; k<idxpx; k++) pmod->vcv[k] = pmod->vcv[k+1];

    if ((pmod->ci == CORC || pmod->ci == HILU) && pmod->ifc) {
	d = 1.0/(1.0 - pmod->rho_in);
	kk = -1;
	for (i=1; i<=nv; i++) {
	    for (j=1; j<=nv; j++) {
		if (j < i) continue;
		kk++;
		if (j == nv) {
		    pmod->vcv[kk] *= d;
		    if (j == i) pmod->vcv[kk] *= d;
		}
	    }
	}
    }

    return 0;
}

/* ............................................................... */

static double dwstat (int order, MODEL *pmod, double **Z)
/*  computes durbin-watson statistic
    opt is the order of autoregression, 0 for OLS and WLS
*/
{
    double diff, ut, ut1;
    double diffsq = 0.0;
    int t;

    if (order) order -= 1;

    if (pmod->ess <= 0.0) return NADBL;

    for (t=pmod->t1+1+order; t<=pmod->t2; t++)  {
        ut = pmod->uhat[t];
        ut1 = pmod->uhat[t-1];
        if (na(ut) || na(ut1) ||
	    (pmod->nwt && (floateq(Z[pmod->nwt][t], 0.0) || 
	    floateq(Z[pmod->nwt][t-1], 0.0)))) continue;
        diff = ut - ut1;
        diffsq += diff * diff;
    }

    return diffsq/pmod->ess;
}

/* ......................................................  */

static double rhohat (int order, int t1, int t2, const double *uhat)
/*  computes first order serial correlation coefficient
    order is the order of autoregression, 0 for OLS.
*/
{
    double ut, ut1, uu, xx, rho;
    int t;

    xx = uu = 0.0;
    if (order) order -= 1;
    for (t=t1+order+1; t<=t2; t++) { 
        ut = uhat[t];
        ut1 = uhat[t-1];
        if (na(ut) || na(ut1)) continue;
        uu += ut * ut1;
        xx += ut1 * ut1;
    }
    if (floateq(xx, 0.0)) return NADBL;
    rho = uu/xx;
    if (rho > 1.0 || rho < -1.0) {
	rho = altrho(order, t1, t2, uhat);
    }
    return rho;
}

/* .........................................................   */

static double altrho (const int order, int t1, int t2, const double *uhat)
/* alternative calculation of rho */
{
    double *ut, *ut1;    
    int t, n = 0;
    double uh, uh1, rho;

    if ((ut = calloc(t2-(t1+order)+1, sizeof *ut)) == NULL) 
	return NADBL;
    if ((ut1 = calloc(t2-(t1+order)+1, sizeof *ut1)) == NULL) 
	return NADBL;

    for (t=t1+order; t<=t2; t++) { 
        uh = uhat[t];
        if (t > 0) uh1 = uhat[t-1];
	else uh1 = NADBL;
        if (na(uh) || na(uh1)) continue;
        ut[n] = uh;
        ut1[n] = uh1;
	n++;
    }
    rho = _corr(n, ut, ut1);
    free(ut);
    free(ut1);
    return rho;
}

/* ........................................................... */

static double corrrsq (int nobs, const double *y, const double *yhat)
/* finds alternative R^2 value when there's no intercept */
{
    double xx;

    xx = _corr(nobs, y, yhat);
    return xx * xx;
}

/* ........................................................... */

static int hatvar (MODEL *pmod, double **Z)
/* finds fitted values and residuals */
{
    int yno, xno, i, t, nwt = pmod->nwt;
    double xx;

    yno = pmod->list[1];
    for (t=pmod->t1; t<=pmod->t2; t++) {
        for (i=1; i<=pmod->ncoeff; i++) {
            xno = pmod->list[i+1];
	    xx = Z[xno][t];
	    if (nwt) xx = xx * Z[nwt][t];
            pmod->yhat[t] += pmod->coeff[i] * xx;
        }
	xx = Z[yno][t];
	if (nwt) xx = xx * Z[nwt][t];
        pmod->uhat[t] = xx - pmod->yhat[t];                
    }
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

/**
 * hilu_corc:
 * @toprho: pointer to receive final rho value.
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: option flag: CORC for Cochrane-Orcutt, HILU for Hildreth-Lu.
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
	       const int opt, PRN *prn)
{
    double rho = 0.0, rho0 = 0.0, diff = 1.0, *uhat;
    double finalrho = 0, ess = 0, essmin = 0, ssr[22], rh[22]; 
    int step, iter = 0, nn = 0, err = 0;
    MODEL corc_model;

    _init_model(&corc_model, pdinfo);

    uhat = malloc(pdinfo->n * sizeof *uhat);
    if (uhat == NULL) return E_ALLOC;

    if (opt == HILU) { /* Do Hildreth-Lu first */
	rho = -1.0;
	diff = 0.1;
	for (step=1; step<=21; step++) {
	    if (rho < -.995) rho = -0.99;
	    else if (rho > 0.995) rho = 0.99;
	    if (step == 2) rho = -.90;
	    clear_model(&corc_model, pdinfo);
	    corc_model = lsq(list, pZ, pdinfo, OLS, 1, rho);
	    if ((err = corc_model.errcode)) {
		free(uhat);
		clear_model(&corc_model, pdinfo);
		return err;
	    }
	    if (step == 1) {
		pprintf(prn, "\n RHO       ESS      RHO       ESS      "
			"RHO       ESS      RHO       ESS     \n");
	    }
	    ssr[nn] = ess = corc_model.ess;
	    rh[nn] = rho;
	    nn++;
	    pprintf(prn, "%5.2f %10.4g", rho, ess);
	    if (step%4 == 0) pprintf(prn, "\n");
	    else _bufspace(3, prn);
	    if (step == 1) essmin = ess;
	    essmin = (ess < essmin)? ess : essmin;
	    if (ess-essmin > -SMALL && ess-essmin < SMALL)
		finalrho = rho;
	    if (rho > 0.989) break;
	    rho = rho + diff;
	}					
	rho0 = rho = finalrho;
	pprintf(prn, _("\n\nESS is minimum for rho = %.2f\n\n"), rho);
	_graphyzx(NULL, ssr, NULL, rh, nn, "ESS", "RHO", NULL, 0, prn); 
	pprintf(prn, _("\n\nFine-tune rho using the CORC procedure...\n\n")); 
    } else { /* Go straight to Cochrane-Orcutt */
	corc_model = lsq(list, pZ, pdinfo, OLS, 1, rho);
	if ((err = corc_model.errcode)) {
	    free(uhat);
	    clear_model(&corc_model, pdinfo);
	    return err;
	}
	rho0 = rho = corc_model.rho;
	pprintf(prn, _("\nPerforming iterative calculation of rho...\n\n"));
    }

    pprintf(prn, "                 ITER       RHO        ESS\n");
    while (diff > 0.001) {
	iter++;
	pprintf(prn, "          %10d %12.5f", iter, rho);
	clear_model(&corc_model, pdinfo);
	corc_model = lsq(list, pZ, pdinfo, OLS, 1, rho);
	if ((err = corc_model.errcode)) {
	    free(uhat);
	    clear_model(&corc_model, pdinfo);
	    return err;
	}
	pprintf(prn, "   %f\n", corc_model.ess);
	corc_model.dw = 1 - rho;
	autores(corc_model.list[1], *pZ, &corc_model, uhat);
	rho = rhohat(1, corc_model.t1, corc_model.t2, uhat);
	diff = (rho > rho0) ? rho - rho0 : rho0 - rho;
	rho0 = rho;
	if (iter == 20) break;
    }
    pprintf(prn, _("                final %11.5f\n\n"), rho);
    free(uhat);
    clear_model(&corc_model, pdinfo);

    *toprho = rho;
    return 0;
}

/* .......................................................... */

static void autores (int i, double **Z, const MODEL *pmod, double *uhat)
{
    int t, v;
    double xx;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	xx = Z[i][t];
	for (v=1; v<=pmod->ncoeff - pmod->ifc; v++)
	    xx -= pmod->coeff[v] * Z[pmod->list[v+1]][t];
	if (pmod->ifc) xx -= pmod->coeff[pmod->ncoeff] / pmod->dw;
	uhat[t] = xx;
    }
}

/**
 * tsls_func:
 * @list: dependent variable plus list of regressors.
 * @pos: position in the list for the separator between list
 *   of variables and list of instruments.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list by means of Two-Stage Least
 * Squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tsls_func (LIST list, const int pos, double ***pZ, 
		 DATAINFO *pdinfo)
{
    int i, j, t, v, ncoeff, *list1, *list2, *newlist, *s1list, *s2list;
    int yno, n = pdinfo->n, orig_nvar = pdinfo->v;
    MODEL tsls;
    XPXXPY xpxxpy;
    CHOLBETA cb;
    double xx, *diag, *yhat;

    /*  printf("tsls list:\n"); */
    /*  printlist(list); */
    /*  printf("pos = %d\n", pos); */

    _init_model(&tsls, pdinfo);

    if ((newlist = malloc(pos * sizeof(int))) == NULL ||
	(list1 = malloc(pos * sizeof(int))) == NULL ||
	(s2list = malloc(pos * sizeof(int))) == NULL ||
	(list2 = malloc((list[0] - pos + 1) * sizeof(int))) == NULL ||
	(s1list = malloc((list[0] - pos + 2) * sizeof(int))) == NULL) {
	tsls.errcode = E_ALLOC;
	return tsls;
    }	

    list1[0] = pos - 1;
    for (i=1; i<pos; i++) list1[i] = list[i];
    tsls_omitzero(list1, *pZ, pdinfo->t1, pdinfo->t2);
    rearrange(list1);
    for (i=0; i<pos; i++) s2list[i] = list1[i];
    list2[0] = list[0] - pos;
    for (i=1; i<=list2[0]; i++) list2[i] = list[i+pos];
/*      printf("list1:\n"); */
/*      printlist(list1); */
/*      printf("list2:\n"); */
/*      printlist(list2); */

    ncoeff = list2[0];
    if (ncoeff < list1[0]-1) {
        sprintf(gretl_errmsg, 
		_("Order condition for identification is not satisfied.\n"
		"varlist 2 needs at least %d more variable(s) not in "
		"varlist1."), list1[0] - 1 - ncoeff);
	free(list1); free(list2);
	free(s1list); free(s2list);
	free(newlist);
	tsls.errcode = E_UNSPEC; 
	return tsls;
    }

    /* now determine which fitted vals to obtain */
    if (tsls_match(list1, list2, newlist)) {
	free(list1); free(list2);
	free(s1list); free(s2list);
	free(newlist);
	strcpy(gretl_errmsg, 
	       _("Constant term is in varlist1 but not in varlist2"));
	tsls.errcode = E_UNSPEC;
	return tsls;
    }

/*      printf("newlist:\n"); */
/*      printlist(newlist);  */
    /* newlist[0] holds the number of new vars to create */
    if (dataset_add_vars(newlist[0], pZ, pdinfo)) {
	free(list1); free(list2);
	free(s1list); free(s2list);
	free(newlist);
	tsls.errcode = E_ALLOC;
	return tsls;	
    }

    /* deal with the variables for which instruments are needed */
    for (i=1; i<=newlist[0]; i++) { 
	yno = newlist[i];
        s1list[0] = ncoeff + 1;
        s1list[1] = yno;
        for (v=2; v<=s1list[0]; v++) s1list[v] = list2[v-1];
	/* run first-stage regression */
/*  	printf("running 1st stage:\n"); */
/*  	printlist(s1list); */
	clear_model(&tsls, pdinfo);
	tsls = lsq(s1list, pZ, pdinfo, OLS, 0, 0.0);
	if (tsls.errcode) {
	    free(list1); free(list2);
	    free(s1list); free(s2list);
	    dataset_drop_vars(newlist[0], pZ, pdinfo);
	    free(newlist);
	    return tsls;
	}
        /* grab fitted values and stick into Z */
	for (j=2; j<=list1[0]; j++) {
	    if (list1[j] == newlist[i]) {
		s2list[j] = orig_nvar + i - 1;
		break;
	    }
	}
	for (t=0; t<n; t++) (*pZ)[orig_nvar+i-1][t] = NADBL;
	for (t=tsls.t1; t<=tsls.t2; t++)
	    (*pZ)[orig_nvar+i-1][t] = tsls.yhat[t];
	strcpy(pdinfo->varname[orig_nvar+i-1], pdinfo->varname[newlist[i]]);
    } 

    /* second-stage regression */
    clear_model(&tsls, pdinfo);
    tsls = lsq(s2list, pZ, pdinfo, OLS, 1, 0.0);
/*      printf("second stage\n"); */
/*      printlist(s2list); */
    if (tsls.errcode) {
	free(list1); free(list2);
	free(s1list); free(s2list);
	dataset_drop_vars(newlist[0], pZ, pdinfo);
	free(newlist);
	return tsls;
    }

    /* special: need to use the original RHS vars to compute residuals 
       and associated statistics */
    yhat = malloc(n * sizeof(double));
    if (yhat == NULL) {
	free(list1); free(list2);
	free(s1list); free(s2list);
	dataset_drop_vars(newlist[0], pZ, pdinfo);
	free(newlist);
	tsls.errcode = E_ALLOC;
	return tsls;
    }
    tsls.ess = 0.0;
    for (t=tsls.t1; t<=tsls.t2; t++) {
	xx = 0.0;
	for (i=1; i<=tsls.ncoeff; i++) {
/*  	    printf("coeff[%d] = %f, Z(%d, %d) = %f\n",  */
/*  		   i, tsls.coeff[i], list1[i+1], t,  */
/*  		   (*pZ)[list1[i+1]*n + t]);  */
	    xx += tsls.coeff[i] * (*pZ)[list1[i+1]][t];
	}
	yhat[t] = xx; 
	tsls.uhat[t] = (*pZ)[tsls.list[1]][t] - xx;
	tsls.ess += tsls.uhat[t] * tsls.uhat[t];
    }
    tsls.sigma = (tsls.ess >= 0.0) ? sqrt(tsls.ess/tsls.dfd) : 0.0;

    xpxxpy = xpxxpy_func(s2list, tsls.t1, tsls.t2, *pZ, 0, 0.0);
    diag = malloc((xpxxpy.nv + 1) * sizeof(double));
    if (diag == NULL) {
	free(list1); free(list2);
	free(s1list); free(s2list);
	clear_model(&tsls, pdinfo);
	dataset_drop_vars(newlist[0], pZ, pdinfo);
	free(newlist);
	free(yhat);
	tsls.errcode = E_ALLOC;
	return tsls;
    }
    cb = cholbeta(xpxxpy);    
    diaginv(cb.xpxxpy, diag);
    for (i=1; i<=tsls.ncoeff; i++) 
	tsls.sderr[i] = tsls.sigma * sqrt(diag[i]); 
    if (diag != NULL) free(diag); 
    if (cb.xpxxpy.xpx != NULL) free(cb.xpxxpy.xpx);
    if (cb.xpxxpy.xpy != NULL) free(cb.xpxxpy.xpy);
    if (cb.coeff != NULL) free(cb.coeff);

    tsls.rsq = corrrsq(tsls.nobs, &(*pZ)[tsls.list[1]][tsls.t1], 
			yhat + tsls.t1);
    tsls.adjrsq = 
	1 - ((1 - tsls.rsq)*(tsls.nobs - 1)/tsls.dfd);
    tsls.fstt = tsls.rsq*tsls.dfd/(tsls.dfn*(1-tsls.rsq));
    _aicetc(&tsls);
    tsls.rho = rhohat(0, tsls.t1, tsls.t2, tsls.uhat);
    tsls.dw = dwstat(0, &tsls, *pZ);

    tsls.ci = TSLS;
    /* put the original list back */
    for (i=2; i<=list1[0]; i++) tsls.list[i] = list1[i];
    /* put the yhats into the model */
    for (t=tsls.t1; t<=tsls.t2; t++) tsls.yhat[t] = yhat[t];

    free(list1); free(list2);
    free(s1list); free(s2list);
    free(yhat); 
    dataset_drop_vars(newlist[0], pZ, pdinfo);
    free(newlist);
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
    int i, l0 = pmod->list[0], listlen, check, shrink;
    MODEL aux;

    _init_model(&aux, pdinfo);

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    /* add uhat1 to data set temporarily */
    for (t=pmod->t1; t<=pmod->t2; t++)
	(*pZ)[v][t] = uhat1[t];

    listlen = (pmod->ifc)? l0 - 1 : l0;
    tmplist = malloc(listlen * sizeof(int));
    if (tmplist == NULL) return E_ALLOC;
    tmplist[0] = listlen - 1;
    for (i=1; i<=tmplist[0]; i++) tmplist[i] = pmod->list[i+1];
    /*  printlist(tmplist); */

    /* now add squares */
    check = xpxgenr(tmplist, pZ, pdinfo, 0, 0);
    if (check < 1) {
	printf(_("generation of squares failed\n"));
	free(tmplist);
	return E_SQUARES;
    }

    free(tmplist);
    tmplist = malloc((check + 2) * sizeof(int));
    if (tmplist == NULL) return E_ALLOC;
    tmplist[0] = pdinfo->v - v - 1;
    for (i=1; i<=tmplist[0]; i++) 
	tmplist[i] = i + v;
    check = _addtolist(pmod->list, tmplist, &list, pdinfo, 999);
    if (check && check != E_VARCHANGE) {
	free(tmplist);
	return check;
    }
    list[1] = v;

    aux = lsq(list, pZ, pdinfo, OLS, 0, 0.);
    check = aux.errcode;
    if (check) shrink = pdinfo->v - v;
    else {
	for (t=aux.t1; t<=aux.t2; t++)
	    (*pZ)[v][t] = aux.yhat[t]; 
	shrink = pdinfo->v - v - 1;
    }
    if (shrink > 0) dataset_drop_vars(shrink, pZ, pdinfo);

    clear_model(&aux, pdinfo);
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

    _init_model(&hsk, pdinfo);

    lo = list[0];         /* number of vars in original list */
    yno = list[1];        /* ID number of original dependent variable */
    ncoeff = list[0] - 1;
    rearrange(list);

    hsk = lsq(list, pZ, pdinfo, OLS, 1, 0.0);
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
	(*pZ)[pdinfo->v - 1][t] = 1.0/sqrt(exp(zz));
    }    

    /* prepare to run weighted least squares */
    hsklist = malloc((lo + 2) * sizeof(int));
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
    for (v=lo+1; v>=3; v--) hsklist[v] = list[v-1];
    /* put the original dependent variable in a position 2 */
    hsklist[2] = yno;

    clear_model(&hsk, pdinfo);
    hsk = lsq(hsklist, pZ, pdinfo, WLS, 1, 0.0);
    hsk.ci = HSK;

    shrink = pdinfo->v - orig_nvar;
    if (shrink > 0) dataset_drop_vars(shrink, pZ, pdinfo);
    free(hsklist);
    free(uhat1);
    return hsk;
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
    double xx, *st, *uhat1, **p;
    MODEL hccm;

    _init_model(&hccm, pdinfo);

    n = pdinfo->n;
    t1 = pdinfo->t1; t2 = pdinfo->t2;
    lo = list[0];
    m3 = 1 + (lo - 1) * (lo - 1);

    /* first try allocating memory */
    if ((st = malloc(list[0] * sizeof *st)) == NULL ||
	(p = malloc(lo * sizeof *p)) == NULL) {
	hccm.errcode = E_ALLOC;
	return hccm;
    }
    for (i=0; i<lo; i++) {
	p[i] = malloc ((t2 - t1 + 1) * sizeof **p);
	if (p[i] == NULL) {
	    free(st);
	    hccm.errcode = E_ALLOC;
	    return hccm;
	}
    }
    uhat1 = malloc(pdinfo->n * sizeof *uhat1);
    if (uhat1 == NULL) {
	free(st);
	hccm.errcode = E_ALLOC;
	return hccm;
    }

    ncoeff = list[0] - 1;
    rearrange(list);

    /* run a regular OLS */
    hccm = lsq(list, pZ, pdinfo, OLS, 1, 0.0);
    if (hccm.errcode) {
	free(uhat1);
	free(st);
	for (i=0; i<lo; i++) free(p[i]);
	free(p);
	return hccm;
    }
    hccm.ci = HCCM;
    nobs = hccm.nobs;

    if (makevcv(&hccm)) {
	hccm.errcode = E_ALLOC;
	free(uhat1);
	free(st);
	for (i=0; i<lo; i++) free(p[i]);
	free(p);
	return hccm;
    } 

    for (i=1; i<=ncoeff; i++) {
	for (t=t1; t<=t2; t++) {
	    xx = 0.0;
	    for (j=1; j<=ncoeff; j++) {
		if (i <= j) index = 1 + ijton(i, j, ncoeff);
		else index = 1 + ijton(j, i, ncoeff);
		xx += hccm.vcv[index] * (*pZ)[list[j+1]][t];
	    }
	    p[i][t] = xx;
	}
    }
    for (t=t1; t<=t2; t++) {
	xx = 0.0;
	for (i=1; i<=ncoeff; i++) xx += (*pZ)[list[i+1]][t] * p[i][t];
	if (floateq(xx, 1.0)) xx = 0.0;
	uhat1[t] = hccm.uhat[t] / (1 - xx);
    }
    for (i=1; i<=ncoeff; i++) {
	xx = 0.0;
	for (t=t1; t<=t2; t++) xx += p[i][t] * uhat1[t];
	st[i] = xx;
    }
    for (t=t1; t<=t2; t++) 
	for (i=1; i<=ncoeff; i++) p[i][t] *= uhat1[t];

    index = 1;
    for (i=1; i<=ncoeff; i++) {
	for (j=i; j<=ncoeff; j++) {
	    xx = 0.0;
	    for (t=t1; t<=t2; t++) xx += p[i][t] * p[j][t];
	    xx = xx * (nobs-1)/nobs -
		(nobs-1) * st[i] * st[j]/(nobs*nobs);
	    if (i == j) hccm.sderr[i] = sqrt(xx);
	    hccm.vcv[index++] = xx;
	}
    }

    free(st);
    free(uhat1);
    for (i=0; i<lo; i++) free(p[i]);
    free(p);

    return hccm;
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
    int lo, ncoeff, yno, i, t, check;
    int shrink, v = pdinfo->v, listlen;
    int *tmplist = NULL, *list = NULL;
    double zz;
    MODEL white;
    int err = 0;

    _init_model(&white, pdinfo);

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

	listlen = (pmod->ifc)? lo - 1 : lo;

	tmplist = malloc(listlen * sizeof *tmplist);
	if (tmplist == NULL) err = E_ALLOC;
	else {
	    tmplist[0] = listlen - 1;
	    for (i=1; i<=tmplist[0]; i++) 
		tmplist[i] = pmod->list[i+1];
	}
    }

    if (!err) {
	/* now add squares */
	check = xpxgenr(tmplist, pZ, pdinfo, 0, 0);
	if (check < 1) {
	    fprintf(stderr, _("generation of squares failed\n"));
	    free(tmplist);
	    err = E_SQUARES;
	}
    }

    if (!err) {
	tmplist = realloc(tmplist, (check + 2) * sizeof *tmplist);
	if (tmplist == NULL) err = E_ALLOC;
	else {
	    tmplist[0] = pdinfo->v - v - 1; 
	    for (i=1; i<=tmplist[0]; i++) 
		tmplist[i] = i + v;
	}
    }

    if (!err) {
	err = _addtolist(pmod->list, tmplist, &list, pdinfo, 999);
	if (err) {
	    if (err != E_VARCHANGE) 
		fprintf(stderr, _("didn't add to list\n"));
	    else {
		err = 0;
	    }
	}
    }

    if (!err) {
	list[1] = v; 
	/* run auxiliary regression and print results */
	white = lsq(list, pZ, pdinfo, OLS, 0, 0.);
	err = white.errcode;
    }

    if (!err) {
	white.aux = AUX_WHITE;
	printmodel(&white, pdinfo, prn);

	if (test != NULL) {
	    strcpy(test->type, _("White's test for heteroskedasticity"));
	    strcpy(test->h_0, _("heteroskedasticity not present"));
	    sprintf(test->teststat, "TR^2 = %f", white.rsq * white.nobs);
	    sprintf(test->pvalue, _("prob(Chi-square(%d) > %f) = %f"), 
		    white.ncoeff - 1, white.rsq * white.nobs, 
		    chisq(white.rsq * white.nobs, white.ncoeff - 1));
	}
    }

    clear_model(&white, pdinfo);
    shrink = pdinfo->v - v;
    if (shrink > 0) dataset_drop_vars(shrink, pZ, pdinfo);

    free(tmplist);
    free(list);

    return err;
}

/**
 * ar_func:
 * @list: dependent variable plus list of regressors and list of lags.
 * @pos: position in list of separator between independent variables and
 * list of lags.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @model_count: count of models estimated so far.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list using the generalized 
 * Cochrane-Orcutt procedure for autoregressive errors.
 * 
 * Returns: #MODEL struct containing the results.
 */

MODEL ar_func (LIST list, const int pos, double ***pZ, 
	       DATAINFO *pdinfo, int *model_count, PRN *prn)
{
    double diff = 100.0, ess = 0, tss = 0, xx;
    int i, j, t, t1, t2, p, vc, yno, ryno = 0, iter = 0;
    int err, lag, n = pdinfo->n, v = pdinfo->v;
    int *arlist = NULL, *reglist = NULL, *reglist2 = NULL, *rholist = NULL;
    MODEL ar, rhomod;

    _init_model(&ar, pdinfo);
    _init_model(&rhomod, pdinfo);

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
    for (i=1; i<=reglist[0]; i++) reglist[i] = list[i+pos];
    rholist[0] = arlist[0] + 1;
    p = arlist[arlist[0]];
    /*  printf("arlist:\n"); printlist(arlist); */
    /*  printf("reglist:\n"); printlist(reglist); */

    if (_hasconst(reglist)) rearrange(reglist);

    /* special case: ar 1 ; ... => use CORC */
    if (arlist[0] == 1 && arlist[1] == 1) {
	err = hilu_corc(&xx, reglist, pZ, pdinfo, CORC, prn);
	if (err) ar.errcode = err;
	else ar = lsq(reglist, pZ, pdinfo, CORC, 1, xx);
	*model_count += 1;
	ar.ID = *model_count;
	printmodel(&ar, pdinfo, prn); 
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

    /*  rearrange(reglist); */ 
    yno = reglist[1];
    /* first pass: estimate model via OLS */
    ar = lsq(reglist, pZ, pdinfo, OLS, 0, 0.0);
    if (ar.errcode) {
	free(reglist);	
	return ar;
    }
    t1 = ar.t1; t2 = ar.t2;
    rholist[1] = v;

    pprintf(prn, _("Generalized Cochrane-Orcutt estimation\n\n"));
    _bufspace(17, prn);
    pprintf(prn, "ITER             ESS           %% CHANGE\n\n");

    /* now loop while ess is changing */
    while (diff > 0.005) {
	iter++;
	for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
	/* special computation of uhat */
	for (t=t1; t<=t2; t++) {
	    xx = (*pZ)[yno][t];
	    for (j=1; j<reglist[0]; j++) {
		xx -= ar.coeff[j] * (*pZ)[reglist[j+1]][t];
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
	if (iter > 1) clear_model(&rhomod, pdinfo);
	rhomod = lsq(rholist, pZ, pdinfo, OLS, 0, 0.0);

	/* and rho-transform the data */
	for (i=1; i<=reglist[0]; i++) {
/*  	    printf("i = %d, vc = %d, t1 = %d\n" */
/*  		   "using var %d, setting var %d\n", */
/*  		   i, vc, t1, reglist[i], vc); */
	    for (t=0; t<n; t++) (*pZ)[vc][t] = NADBL;
	    for (t=t1+p; t<=t2; t++) {
		xx = (*pZ)[reglist[i]][t];
		for (j=1; j<=arlist[0]; j++) {
		    lag = arlist[j];
		    xx -= rhomod.coeff[j] * (*pZ)[reglist[i]][t-lag];
		}
		(*pZ)[vc][t] = xx;
	    }
	    reglist2[i] = vc;
	    vc++;
	}

	/* estimate the transformed model */
	clear_model(&ar, pdinfo);
	ar = lsq(reglist2, pZ, pdinfo, OLS, 0, 0.0);

        if (iter > 1) diff = 100 * (ar.ess - ess)/ess;
        if (diff < 0.0) diff = -diff;
	ess = ar.ess;
	pprintf(prn, "%16c%3d %20f ", ' ', iter, ess);
	if (iter > 1) pprintf(prn, "%13.3f\n", diff);
	else pprintf(prn, _("      undefined\n")); 
	if (iter == 20) break;
    } /* end loop */

    for (i=0; i<=reglist[0]; i++) ar.list[i] = reglist[i];
    ar.ifc = _hasconst(reglist);
    if (ar.ifc) ar.dfn -= 1;
    ar.ci = AR;
    *model_count += 1;
    ar.ID = *model_count;
    printmodel(&ar, pdinfo, prn);

    pprintf(prn, _("Estimates of the AR coefficients:\n\n"));
    xx = 0.0;
    for (i=1; i<=arlist[0]; i++) {
	_print_rho(arlist, &rhomod, i, prn);
	xx += rhomod.coeff[i];
    }
    pprintf(prn, _("\nSum of AR coefficients = %f\n\n"), xx);
    ar.rho_in = xx;

    /* special computation of fitted values */
    for (t=t1; t<=t2; t++) {
	xx = 0.0;
	for (j=2; j<=reglist[0]; j++) 
	    xx += ar.coeff[j-1] * (*pZ)[reglist[j]][t];
	ar.uhat[t] = (*pZ)[yno][t] - xx;
	for (j=1; j<=arlist[0]; j++)
	    if (t - t1 >= arlist[j]) 
		xx += rhomod.coeff[j] * ar.uhat[t - arlist[j]];
	ar.yhat[t] = xx;
	/*  printf("yhat[%d] = %f\n", t, ar.yhat[t]); */
    }

    for (t=t1; t<=t2; t++) 
	ar.uhat[t] = (*pZ)[yno][t] - ar.yhat[t];
    ar.rsq = corrrsq(ar.nobs, &(*pZ)[reglist[1]][ar.t1], ar.yhat + ar.t1);
    ar.adjrsq = 
	1 - ((1 - ar.rsq)*(ar.nobs - 1)/ar.dfd);
    /*  ar.fstt = ar.rsq*ar.dfd/(ar.dfn*(1 - ar.rsq)); */
    /* special computation of TSS */
    xx = _esl_mean(ar.t1, ar.t2, (*pZ)[ryno]);
    for (t=ar.t1; t<=ar.t2; t++)
	tss += ((*pZ)[ryno][t] - xx) * ((*pZ)[ryno][t] - xx);
    ar.fstt = ar.dfd * (tss - ar.ess) / (ar.dfn * ar.ess);
    _aicetc(&ar);
    ar.dw = dwstat(p, &ar, *pZ);
    ar.rho = rhohat(p, ar.t1, ar.t2, ar.uhat);

    _print_ar(&ar, prn);

    dataset_drop_vars(arlist[0] + 1 + reglist[0], pZ, pdinfo);
    free(reglist);
    free(reglist2);
    free(rholist);

    ar.rhot = malloc((1 + p) * sizeof(double));
    if (ar.rhot == NULL) ar.errcode = E_ALLOC;
    else {
	for (i=1; i<=arlist[0]; i++) 
	    ar.rhot[i] = rhomod.coeff[i];
    }
    ar.arlist = NULL;
    copylist(&(ar.arlist), arlist);
    free(arlist);
    clear_model(&rhomod, pdinfo);
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
        if (_iszero(pmod->t1, pmod->t2, Z[lv])) {
	    list_exclude(v, pmod->list);
	    sprintf(vnamebit, "%s ", pdinfo->varname[lv]);
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
    int v, lv;

    for (v=2; v<=list[0]; v++) {
        lv = list[v];
        if (_iszero(t1, t2, Z[lv])) 
	    list_exclude(v, list);
    }
}

/* .........................................................   */

static void rearrange (int *list)
/* checks a list for a constant term (ID # 0), and if present, 
   move it to the last position
*/
{
    int lo = list[0], v;

    for (v=2; v<=lo; v++) {
        if (list[v] == 0)  {
            list_exclude(v, list);
            list[0] = lo;
            list[lo] = 0;
            return;
        }
    }
}

/* ...........................................................*/

static int zerror (int t1, int t2, int yno, int nwt, double ***pZ)
{
    double xx, yy;
    int t;

    xx = _esl_mean(t1, t2, (*pZ)[yno]);
    yy = _esl_stddev(t1, t2, (*pZ)[yno]);
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
    int i, c, t;
    char depvar[9], othervar[9];

    /* this may be an auxiliary regression */
    if (pdinfo->extra) return 0;

    strcpy(depvar, pdinfo->varname[list[1]]);

    for (i=2; i<=list[0]; i++) {
	strcpy(othervar, pdinfo->varname[list[i]]);
	c = haschar('_', othervar);
	if (c > 0 && isdigit(othervar[c+1]) 
	    && strncmp(depvar, othervar, c-1) == 0) {
	    /* strong candidate for lagged depvar, but make sure */
	    for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) 
		if ((*pZ)[list[1]][t-1] 
		    != (*pZ)[list[i]][t]) return 0;
	    return i; 
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
	    } else sum += Z[yno][t]; 
	}
    }
    return sum/m;
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
	    sumsq += xx*xx;
    }
    sumsq = (m > 1)? sumsq/(m-1) : 0.0;
    if (sumsq >= 0) return sqrt(sumsq);
    else return NADBL;
}

/**
 * arch:
 * @order: lag order for ARCH process.
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @model_count: count of models estimated so far.
 * @prn: gretl printing struct.
 * @test: hypothesis test results struct.
 *
 * Estimate the model given in @list via OLS, and test for Auto-
 * Regressive Conditional Heteroskedasticity.  If the latter is
 * significant, re-restimate the model using weighted least
 * squares.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch (int order, LIST list, double ***pZ, DATAINFO *pdinfo, 
	    int *model_count, PRN *prn, GRETLTEST *test)
{
    MODEL archmod;
    int *wlist = NULL, *arlist = NULL;
    int i, t, nwt, nv, n = pdinfo->n;
    double LM, xx;
    int err = 0;

    _init_model(&archmod, pdinfo);

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
	archmod = lsq(list, pZ, pdinfo, OLS, 0, 0.0);
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
	clear_model(&archmod, pdinfo);
	archmod = lsq(arlist, pZ, pdinfo, OLS, 1, 0.0);
	err = archmod.errcode;
    }

    if (!err) {
	/* print results */
	archmod.aux = AUX_ARCH;
	archmod.order = order;
	printmodel(&archmod, pdinfo, prn);
	pprintf(prn, _("No of obs. = %d, unadjusted R^2 = %f\n"),
		archmod.nobs, archmod.rsq);
	LM = archmod.nobs * archmod.rsq;
	xx = chisq(LM, order);

	if (test != NULL) {
	    sprintf(test->type, _("Test for ARCH of order %d"), order);
	    strcpy(test->h_0, _("no ARCH effect is present"));
	    sprintf(test->teststat, "TR^2 = %f", LM);
	    sprintf(test->pvalue, _("prob(Chi-square(%d) > %f) = %f"), 
		    order, LM, xx);
	}

	pprintf(prn, _("LM test statistic (%f) is distributed as Chi-square "
		"(%d)\nArea to the right of LM = %f  "), LM, order, xx);
	if (xx > 0.1) 
	    pprintf(prn, _("\nARCH effect is insignificant at the 10 "
		    "percent level.\nWeighted estimation not done.\n"));
	else {
	    pprintf(prn, _("\nARCH effect is significant at the 10 "
		    "percent level.\n"));
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
		clear_model(&archmod, pdinfo);
		archmod = lsq(wlist, pZ, pdinfo, WLS, 1, 0.0);
		if (model_count != NULL) {
		    *model_count += 1;
		    archmod.ID = *model_count;
		} else 
		    archmod.ID = -1;
		archmod.ci = ARCH;
		archmod.order = order;
		printmodel(&archmod, pdinfo, prn);
	    }
	}
    }

    if (arlist != NULL) free(arlist);
    if (wlist != NULL) free(wlist);
    dataset_drop_vars(order + 1, pZ, pdinfo); 
    return archmod;
}

