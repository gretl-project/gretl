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

#include "libgretl.h"
#include "internal.h"

enum {
    LMAX_ONE,
    LMAX_HUNDRED,
    LMAX_BAD
};

static int get_lmax (const double *y, const DATAINFO *pdinfo)
{
    int t, ret = LMAX_ONE;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (y[t] <= 0.0 || y[t] >= 100.0) {
	    ret = LMAX_BAD;
	    break;
	}
	if (y[t] >= 1.0) {
	    ret = LMAX_HUNDRED;
	}
    }
	    
    return ret;
}

static int make_logistic_depvar (double ***pZ, DATAINFO *pdinfo, 
				 int dv, int lmax)
{
    int t, v = pdinfo->v;
    int err;

    err = dataset_add_vars(1, pZ, pdinfo);
    if (err) return 1;

    for (t=0; t<pdinfo->n; t++) {
	double p = (*pZ)[dv][t];

	if (na(p)) continue;
	if (lmax == LMAX_HUNDRED) p /= 100.0;
	(*pZ)[v][t] = log(p / (1.0 - p));
    }

    return 0;
}

static int rewrite_logistic_stats (const double **Z, const DATAINFO *pdinfo,
				   MODEL *pmod, int dv, int lmax)
{
    int t;
    double x;

    pmod->ybar = _esl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = _esl_stddev(pmod->t1, pmod->t2, Z[dv]);

    pmod->ess = 0.0;
    for (t=0; t<pdinfo->n; t++) {
	x = pmod->yhat[t];
	if (na(x)) continue;
	pmod->yhat[t] = 1.0 / (1.0 + exp(-x));
	if (lmax == LMAX_HUNDRED) {
	    pmod->yhat[t] *= 100.0;
	} 
	pmod->uhat[t] = Z[dv][t] - pmod->yhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->tss += (Z[dv][t] - pmod->ybar) * (Z[dv][t] - pmod->ybar);
    }

    pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    pmod->list[1] = dv;

    if (lmax == LMAX_HUNDRED) {
	gretl_model_set_int(pmod, "logistic_percent", 1);
    }

    pmod->ci = LOGISTIC;

    gretl_aic_etc(pmod);

    return 0;
}

MODEL logistic_estimate (int *list, double ***pZ, DATAINFO *pdinfo) 
{
    int lmax;
    int dv = list[1];
    MODEL lmod;

    gretl_model_init(&lmod, pdinfo); 

    lmax = get_lmax((*pZ)[dv], pdinfo);
 
    if (lmax == LMAX_BAD) {
	lmod.errcode = E_DATA;
	return lmod;
    }

    if (make_logistic_depvar(pZ, pdinfo, dv, lmax)) {
	lmod.errcode = E_ALLOC;	
	return lmod;
    }

    list[1] = pdinfo->v - 1;

    lmod = lsq(list, pZ, pdinfo, OLS, 1, 0.0);
    if (!lmod.errcode) {
	rewrite_logistic_stats((const double **) *pZ, pdinfo, &lmod,
			       dv, lmax);
    }

    dataset_drop_vars(1, pZ, pdinfo);
    
    return lmod;
}



    

    

    
    
