/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* The beginings of a GARCH plugin using the Fiorentini, Calzolari and
   Panattoni (fcp) fortran code.
*/

#include "libgretl.h"
#include "internal.h"

#include "fcp.h"

static void add_garch_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				const int *list)
{
    int p = list[1];
    int q = list[2];
    int r = list[0] - 4;
    int i, j, np = 3 + p + q + r;

    free(pmod->list);
    pmod->list = NULL;
    copylist(&pmod->list, list);

    pmod->params = malloc(np * sizeof pmod->params);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    pmod->nparams = np;

    for (i=0; i<np; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    pmod->errcode = E_ALLOC;
	    return;
	}
    }

    strcpy(pmod->params[0], pdinfo->varname[pmod->list[4]]);
    strcpy(pmod->params[1], pdinfo->varname[0]);

    j = 2;
    for (i=0; i<r; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[pmod->list[5+i]]);
    }

    strcpy(pmod->params[j++], "alpha(0)");

    for (i=0; i<p; i++) {
	sprintf(pmod->params[j++], "alpha(%d)", i + 1);
    }
    for (i=0; i<q; i++) {
	sprintf(pmod->params[j++], "beta(%d)", i + 1);
    }
}

static int write_garch_stats (MODEL *pmod, const double **Z,
			      const DATAINFO *pdinfo,
			      const int *list, const double *theta, 
			      int nparam, const double *res)
{
    int err = 0;
    double *coeff, *sderr;
    int maxlag, i, ynum = list[4];

    if (list[2] > list[1]) maxlag = list[2];
    else maxlag = list[1];

    coeff = realloc(pmod->coeff, nparam * sizeof *pmod->coeff);
    sderr = realloc(pmod->sderr, nparam * sizeof *pmod->sderr);

    if (coeff == NULL || sderr == NULL) return 1;

    for (i=0; i<nparam; i++) {
	coeff[i] = theta[i+1];
	sderr[i] = theta[i+nparam+1];
    }
    
    pmod->coeff = coeff;
    pmod->sderr = sderr;
   
    pmod->ncoeff = nparam;

    pmod->ess = 0.0;
    for (i=pmod->t1; i<=pmod->t2; i++) {
	pmod->uhat[i] = res[i+maxlag];
	pmod->ess += pmod->uhat[i] * pmod->uhat[i];
	pmod->yhat[i] =  Z[ynum][i] - pmod->uhat[i];
    }

    pmod->sigma = NADBL;
    
    pmod->ci = GARCH;
    
    add_garch_varnames(pmod, pdinfo, list);

    return err;
}

int do_fcp (const int *list, const double **Z, 
	    const DATAINFO *pdinfo, MODEL *pmod,
	    PRN *prn)
{
    int t1, t2;
    double **X;
    int nx, nobs;
    double *yhat;
    int ncoeff;
    double *coeff;
    double *vc, *res, *res2;
    double *ystoc;
    double *amax;
    double *b;
    int err = 0, iters = 0;
    int maxlag;
    int nobsmod;
    int nparam;
    int i, p, q, ynum;

    t1 = pmod->t1;
    t2 = pmod->t2;
    ncoeff = pmod->ncoeff;
    nx = ncoeff - 1;

    p = list[1];
    q = list[2];
    ynum = list[4];
    maxlag = (p > q)? p : q; 

    nparam = ncoeff + p + q + 1;

    nobs = t2 + 1;
    nobsmod = nobs + maxlag;

    yhat = malloc(nobsmod * sizeof *yhat);
    ystoc = malloc(nobsmod * sizeof *ystoc);

    res2 = malloc(nobsmod * sizeof *res2);
    for (i=0; i<nobsmod; i++) {
	res2[i] = 0.0;
    }

    res = malloc(nobsmod * sizeof *res);
    for (i=0; i<nobsmod; i++) {
	res[i] = 0.0;
    }    

    amax = malloc(nobsmod * sizeof *amax);
    for (i=0; i<nobsmod; i++) {
	amax[i] = 0.0;
    }

    coeff = malloc(ncoeff * sizeof *coeff);
    b = malloc(ncoeff * sizeof *b);
    for (i=0; i<ncoeff; i++) {
	coeff[i] = b[i] = 0.0;
    }    

    vc = malloc((ncoeff * ncoeff) * sizeof *vc);
    for (i=0; i<ncoeff; i++) {
	vc[i] = 0.0;
    } 

    if (nx > 0) {
	X = malloc(nx * sizeof *X);
	/* FIXME */
    } else {
	X = NULL;
    }

    for (i=0; i<maxlag; i++) {
	ystoc[i] = yhat[i] = 0.0;
    }    
    for (i=maxlag; i<nobsmod; i++) {
	ystoc[i] = yhat[i] = Z[ynum][i-maxlag];
    }

    /* initial coefficients from OLS */
    for (i=0; i<ncoeff; i++) {
	coeff[i] = pmod->coeff[i];
    }

    /* initialize elements of alpha, beta such that 
       alpha_0/(1 - alpha_1 - beta_1) = unconditional
       variance of y (?)
    */
    amax[0] = pmod->sigma * pmod->sigma;
    amax[1] = p;
    amax[2] = q; 
    for (i=0; i<p+q; i++) {
	/* initial alpha, beta values */
	amax[3+i] = 0.1;
    }

    /* Need to set starting point high enough to allow for lags? */

    err = garch_estimate(t1 + maxlag, t2 + maxlag, 
			 nobsmod, 
			 (const double **) X, nx, 
			 yhat, 
			 coeff, ncoeff, 
			 vc, 
			 res2, 
			 res, 
			 ystoc, 
			 amax, b, &iters, prn);

    if (err != 0) {
	fprintf(stderr, "garch_estimate returned %d\n", err);
    }

    if (err == 0) {
	int nparam = ncoeff + p + q + 1;

	for (i=1; i<=nparam; i++) {
	    pprintf(prn, "theta[%d]: %#14.6g (%#.6g)\n", i-1, amax[i], 
		    amax[i+nparam]);
	}
    }

    write_garch_stats(pmod, Z, pdinfo, list, amax, nparam, res);
    gretl_model_set_int(pmod, "iters", iters);
    pmod->lnL = amax[0];

    free(X); /* FIXME */
    free(yhat);
    free(coeff);
    free(vc);
    free(res2);
    free(res);
    free(ystoc);
    free(amax);
    free(b);

    return 0;
}

static int *make_ols_list (const int *list)
{
    int *olist;
    int i, ifc = 0;

    olist = malloc((list[0] - 1) * sizeof *olist);
    if (olist == NULL) return NULL;

    olist[0] = list[0] - 3;
    for (i=4; i<=list[0]; i++) {
	olist[i-3] = list[i];
	if (list[i] == 0) ifc = 1;
    }

    /* add constant, if absent */
    if (!ifc) {
	olist[0] += 1;
	olist[olist[0]] = 0;
    }

    return olist;
}

/* the driver function for the plugin */

MODEL garch_model (int *list, double ***pZ, DATAINFO *pdinfo,
		   PRN *prn) 
{
    MODEL model;
    int *ols_list;

    gretl_model_init(&model, NULL);

    ols_list = make_ols_list(list);
    if (ols_list == NULL) {
	model.errcode = E_ALLOC;
	return model;
    }

    /* run initial OLS */
    model = lsq(ols_list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (model.errcode) {
        return model;
    } 

    do_fcp(list, (const double **) *pZ, pdinfo, &model, prn); 

    return model;
}
