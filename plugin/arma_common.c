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
 */

#include "usermat.h"

#define MAX_ARMA_ORDER 128
#define MAX_ARIMA_DIFF 2

static void 
real_arima_difference_series (double *dx, const double *x,
			      int t1, int t2, int *delta, 
			      int k);

static void 
arma_info_init (arma_info *ainfo, gretlopt opt, 
		const int *pqspec, const DATASET *dset)
{
    ainfo->yno = 0;
    ainfo->flags = 0;
    ainfo->pflags = 0;
    ainfo->alist = NULL;

    if (opt & OPT_X) {
	/* we got --x-12-arima */
	ainfo->flags |= ARMA_X12A;
    }    

    if (!(opt & OPT_C)) {
	/* we didn't get --conditional */
	ainfo->flags |= ARMA_EXACT;
    }

    ainfo->ll = NADBL;

    ainfo->pqspec = pqspec;
    ainfo->pmask = NULL;
    ainfo->qmask = NULL;

    ainfo->p = 0;
    ainfo->d = 0;
    ainfo->q = 0;
    ainfo->P = 0;
    ainfo->D = 0;
    ainfo->Q = 0; 
    
    ainfo->np = 0;
    ainfo->nq = 0;

    ainfo->maxlag = 0;
    ainfo->ifc = 0;
    ainfo->nexo = 0;
    ainfo->nc = 0;

    ainfo->t1 = dset->t1;
    ainfo->t2 = dset->t2;
    ainfo->pd = dset->pd;
    ainfo->T = 0;
    ainfo->r0 = 0;

    ainfo->y = NULL;
    ainfo->e = NULL;
    ainfo->Z = NULL;
    ainfo->yscale = 1.0;

    ainfo->xlist = NULL;
    ainfo->misslist = NULL;
    ainfo->dX = NULL;
    ainfo->G = NULL;
    ainfo->V = NULL;

    ainfo->n_aux = 0;
    ainfo->aux = NULL;

    ainfo->prn = NULL;
}

static void arma_info_cleanup (arma_info *ainfo)
{
    free(ainfo->alist);
    free(ainfo->pmask);
    free(ainfo->qmask);
    free(ainfo->e);
    free(ainfo->Z);
    free(ainfo->xlist);
    free(ainfo->misslist);

    gretl_matrix_free(ainfo->dX);
    gretl_matrix_free(ainfo->G);
    gretl_matrix_free(ainfo->V);

    if (arima_ydiff(ainfo)) {
	free(ainfo->y);
    }

    doubles_array_free(ainfo->aux, ainfo->n_aux);
}

enum {
    AR_MASK,
    MA_MASK
};

/* Create a mask for skipping certain intermediate lags, 
   AR or MA.  This function also sets ainfo->np and ainfo->nq,
   which record the actual number of non-seasonal AR and MA
   lags used.
*/

static char *mask_from_list (const int *list, 
			     arma_info *ainfo,
			     int m, int *err)
{
    int mlen = (m == AR_MASK)? ainfo->p : ainfo->q;
    int nv = 0, nmax = 0;
    char *mask;
    int i, k;

    mask = malloc(mlen + 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<mlen; i++) {
	mask[i] = '0';
    }
    mask[mlen] = '\0';

    for (i=1; i<=list[0]; i++) {
	k = list[i];
	if (k > 0) {
	    mask[k-1] = '1'; 
	    nv++;
	    if (k > nmax) {
		nmax = k;
	    }
	}
    }

    if (m == AR_MASK) {
	ainfo->p = nmax;
	ainfo->np = nv;
    } else {
	ainfo->q = nmax;
	ainfo->nq = nv;
    }

    if (nv == 0) {
	free(mask);
	mask = NULL;
    }

    return mask;
}

static int arma_make_masks (arma_info *ainfo)
{
    int *plist = NULL, *qlist = NULL;
    int err = 0;

    if (ainfo->pqspec != NULL) {
	if (gretl_list_has_separator(ainfo->pqspec)) {
	    gretl_list_split_on_separator(ainfo->pqspec, &plist, &qlist);
	} else {
	    plist = gretl_list_copy(ainfo->pqspec);
	}
    }

    if (ainfo->p > 0) {
	ainfo->np = ainfo->p;
	if (plist != NULL && plist[0] > 0) {
	    ainfo->pmask = mask_from_list(plist, ainfo, AR_MASK, &err);
	}
    }

    if (ainfo->q > 0 && !err) {
	ainfo->nq = ainfo->q;
	if (qlist != NULL && qlist[0] > 0) {
	    ainfo->qmask = mask_from_list(qlist, ainfo, MA_MASK, &err);
	}
    }

    free(plist);
    free(qlist);

    return err;
}

int arma_list_y_position (arma_info *ainfo)
{
    int ypos;

    if (arma_is_arima(ainfo)) {
	ypos = arma_has_seasonal(ainfo) ? 9 : 5;
    } else {
	ypos = arma_has_seasonal(ainfo) ? 7 : 4;
    }

    return ypos;
}

static int arima_integrate (double *dx, const double *x,
			    int t1, int t2, int d, int D, int s)
{
    double *ix;
    int *c;
    int k = d + s * D;
    int i, t;

    ix = malloc((t2 + 1) * sizeof *ix);
    if (ix == NULL) {
	return E_ALLOC;
    }

    c = arima_delta_coeffs(d, D, s);
    if (c == NULL) {
	free(ix);
	return E_ALLOC;
    }    

    for (t=0; t<t1; t++) {
	ix[t] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	ix[t] = dx[t];
	for (i=0; i<k; i++) {
	    if (c[i] != 0) {
		ix[t] += c[i] * x[t-i-1];
	    }
	}
    }

    /* transcribe integrated result back into "dx" */
    for (t=0; t<=t2; t++) {
	if (t < t1) {
	    dx[t] = NADBL;
	} else {
	    dx[t] = ix[t];
	}
    }

    free(ix);
    free(c);

    return 0;
}

static void ainfo_data_to_model (arma_info *ainfo, MODEL *pmod)
{
    pmod->ifc = ainfo->ifc;
    pmod->dfn = ainfo->nc - pmod->ifc;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = ainfo->nc;

    if (arma_has_seasonal(ainfo)) {
	gretl_model_set_int(pmod, "arma_P", ainfo->P);
	gretl_model_set_int(pmod, "arma_Q", ainfo->Q);
	gretl_model_set_int(pmod, "arma_pd", ainfo->pd);	
    }

    if (ainfo->d > 0 || ainfo->D > 0) {
	gretl_model_set_int(pmod, "arima_d", ainfo->d);
	gretl_model_set_int(pmod, "arima_D", ainfo->D);
    }

    if (ainfo->nexo > 0) {
	gretl_model_set_int(pmod, "armax", 1);
    }

    if (ainfo->pmask != NULL) {
	gretl_model_set_string_as_data(pmod, "pmask", 
				       gretl_strdup(ainfo->pmask));
    }

    if (ainfo->qmask != NULL) {
	gretl_model_set_string_as_data(pmod, "qmask", 
				       gretl_strdup(ainfo->qmask));
    }
}

static void arma_depvar_stats (MODEL *pmod, arma_info *ainfo,
			       const DATASET *dset)
{
    if (arma_is_arima(ainfo) && !arima_ydiff(ainfo)) {
	/* calculate differenced y for stats */
	int d = ainfo->d, D = ainfo->D;
	int T = pmod->t2 - pmod->t1 + 1;
	double *dy = malloc(T * sizeof *dy);
	int *delta = arima_delta_coeffs(d, D, ainfo->pd);

	if (dy != NULL && delta != NULL) {
	    int k = d + ainfo->pd * D;

	    real_arima_difference_series(dy, dset->Z[ainfo->yno], 
					 pmod->t1, pmod->t2, delta, k);
	    pmod->ybar = gretl_mean(0, T - 1, dy);
	    pmod->sdy = gretl_stddev(0, T - 1, dy);
	}
	free(dy);
	free(delta);
    } else {
	pmod->ybar = gretl_mean(pmod->t1, pmod->t2, ainfo->y);
	pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, ainfo->y);
    }
}

#define USE_ARIMA_INTEGRATE 1

/* write the various statistics from ARMA estimation into
   a gretl MODEL struct */

void write_arma_model_stats (MODEL *pmod, arma_info *ainfo,
			     const DATASET *dset)
{
    double mean_error;
    int do_criteria = 1;
    int t;

    pmod->ci = ARMA;

    ainfo_data_to_model(ainfo, pmod);

    free(pmod->list);
    pmod->list = gretl_list_copy(ainfo->alist);

    if (!arma_least_squares(ainfo)) {
	arma_depvar_stats(pmod, ainfo, dset);
    }

    mean_error = pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(ainfo->y[t]) && !na(pmod->uhat[t])) {
#if USE_ARIMA_INTEGRATE == 0
	    if (arma_is_arima(ainfo) && arima_ydiff(ainfo)) {
		pmod->yhat[t] = Z[ainfo->yno][t] - pmod->uhat[t];
	    }
#else
	    pmod->yhat[t] = ainfo->y[t] - pmod->uhat[t];
#endif
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    mean_error += pmod->uhat[t];
	} 
    }

#if USE_ARIMA_INTEGRATE
    if (arma_is_arima(ainfo) && arima_ydiff(ainfo)) {
	arima_integrate(pmod->yhat, dset->Z[ainfo->yno], 
			pmod->t1, pmod->t2, 
			ainfo->d, ainfo->D, ainfo->pd);
    }
#endif

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    if (na(pmod->sigma)) {
	/* in x12a or native exact cases this is already done */
	pmod->sigma = sqrt(pmod->ess / pmod->nobs);
    } 

    pmod->rsq = pmod->adjrsq = pmod->fstt = pmod->chisq = NADBL;
    pmod->tss = NADBL;

    if (arma_least_squares(ainfo)) {
	/* not applicable */
	do_criteria = 0;
    } else if (arma_by_x12a(ainfo) && !na(pmod->criterion[C_AIC])) {
	/* already given by x12a */
	do_criteria = 0;
    }

    if (do_criteria) {
	mle_criteria(pmod, 1);
    }

    gretl_model_add_arma_varnames(pmod, dset, ainfo->yno,
				  ainfo->p, ainfo->q, 
				  ainfo->pmask, ainfo->qmask,
				  ainfo->P, ainfo->Q,
				  ainfo->nexo);
}

static void calc_max_lag (arma_info *ainfo)
{
    if (arma_exact_ml(ainfo)) {
	ainfo->maxlag = ainfo->d + ainfo->D * ainfo->pd;
    } else {
	/* conditional ML */
	int pmax = ainfo->p + ainfo->P * ainfo->pd;
	int dmax = ainfo->d + ainfo->D * ainfo->pd;

	ainfo->maxlag = pmax + dmax;
    }

#if ARMA_DEBUG
    fprintf(stderr, "calc_max_lag: ainfo->maxlag = %d\n", ainfo->maxlag);
#endif
}

static int arma_adjust_sample (arma_info *ainfo, 
			       const DATASET *dset,
			       int *missv, int *misst)
{
    int *list = ainfo->alist;
    int ypos = arma_list_y_position(ainfo);
    int t0, t1 = dset->t1, t2 = dset->t2;
    int i, vi, vlmax, k, t;
    int missing;
    int err = 0;

#if ARMA_DEBUG
    fprintf(stderr, "arma_adjust_sample: at start, t1=%d, t2=%d, maxlag = %d\n",
	    t1, t2, ainfo->maxlag);
#endif

    t0 = t1 - ainfo->maxlag;
    if (t0 < 0) {
	t1 -= t0;
    }

    /* list position of last var to check for lags */
    if (arma_xdiff(ainfo)) {
	vlmax = list[0];
    } else {
	vlmax = ypos;
    }

    /* advance the starting point if need be */

    for (t=t1; t<=t2; t++) {
	missing = 0;
	for (i=ypos; i<=list[0] && !missing; i++) {
	    vi = list[i];
	    if (na(dset->Z[vi][t])) {
		/* current value missing */
		missing = 1;
	    }
	    if (i <= vlmax) {
		for (k=1; k<=ainfo->maxlag && !missing; k++) {
		    if (na(dset->Z[vi][t-k])) {
			/* lagged value missing */
			missing = 1;
		    }
		}
	    }
	}
	if (missing) {
	    t1++;
	} else {
	    break;
	}
    }

    /* retard the ending point if need be */

    for (t=t2; t>=t1; t--) {
	missing = 0;
	for (i=ypos; i<=list[0] && !missing; i++) {
	    vi = list[i];
	    if (na(dset->Z[vi][t])) {
		missing = 1;
	    }
	}
	if (missing) {
	    t2--;
	} else {
	    break;
	}
    }

    missing = 0;

    /* check for missing obs within the adjusted sample range */
    for (t=t1; t<t2; t++) {
	int tmiss = 0;

	for (i=ypos; i<=list[0]; i++) {
	    vi = list[i];
	    if (na(dset->Z[vi][t])) {
		if (missv != NULL && misst != NULL && *missv == 0) {
		    /* record info on first missing obs */
		    *missv = vi;
		    *misst = t + 1;
		}
		tmiss = 1;
	    }
	}
	if (tmiss) {
	    missing++;
	}
    }

    if (missing > 0 && !arma_na_ok(ainfo)) {
	err = E_MISSDATA;
    }

    if (!err) {
	ainfo->fullT = t2 - t1 + 1;
	ainfo->T = ainfo->fullT - missing;
	if (ainfo->T <= ainfo->nc) {
	    /* insufficient observations */
	    err = E_DF; 
	}
    }

    if (!err) {
#if ARMA_DEBUG
	fprintf(stderr, "arma_adjust_sample: at end, t1=%d, t2=%d\n",
		t1, t2);
#endif
	ainfo->t1 = t1;
	ainfo->t2 = t2;
    }

    return err;
}

/* remove the intercept from list of regressors */

static int arma_remove_const (arma_info *ainfo, 
			      const DATASET *dset)
{
    int *list = ainfo->alist;
    int seasonal = arma_has_seasonal(ainfo);
    int diffs = arma_is_arima(ainfo);
    int xstart, ret = 0;
    int i, j;

    if (diffs) {
	xstart = (seasonal)? 10 : 6;
    } else {
	xstart = (seasonal)? 8 : 5;
    }

    for (i=xstart; i<=list[0]; i++) {
	if (list[i] == 0 || true_const(list[i], dset)) {
	    for (j=i; j<list[0]; j++) {
		list[j] = list[j+1];
	    }
	    list[0] -= 1;
	    ret = 1;
	    break;
	}
    }

    return ret;
}

static int check_arma_sep (arma_info *ainfo, int sep1)
{
    int *list = ainfo->alist;
    int sep2 = (sep1 == 3)? 6 : 8;
    int i, err = 0;

    for (i=sep1+1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    if (i == sep2) {
		/* there's a second list separator in the right place:
		   we've got a seasonal specification */
		set_arma_has_seasonal(ainfo);
	    } else {
		err = 1;
	    }
	}
    }

    if (!err && sep1 == 4) {
	/* check for apparent but not "real" arima spec */
	if (arma_has_seasonal(ainfo)) {
	    if (list[2] == 0 && list[6] == 0) {
		gretl_list_delete_at_pos(list, 2);
		gretl_list_delete_at_pos(list, 5);
		unset_arma_is_arima(ainfo);
	    }
	} else {
	    if (list[2] == 0) {
		gretl_list_delete_at_pos(list, 2);
		unset_arma_is_arima(ainfo);
	    }
	}
    }

    return err;
}

static int arma_add_xlist (arma_info *ainfo, int ypos)
{
    int i, err = 0;

    ainfo->xlist = gretl_list_new(ainfo->nexo);

    if (ainfo->xlist == NULL) {
	err = E_ALLOC;
    } else {
	for (i=1; i<=ainfo->nexo; i++) {
	    ainfo->xlist[i] = ainfo->alist[ypos + i];
	}
    }

    return err;
}

#define count_arma_coeffs(a) (a->ifc + a->np + a->nq + a->P + a->Q + a->nexo)

static int check_arma_list (arma_info *ainfo, 
			    const DATASET *dset,
			    gretlopt opt)
{
    int *list = ainfo->alist;
    int ypos = arma_has_seasonal(ainfo) ? 7 : 4;
    int armax = (list[0] > ypos);
    int hadconst = 0;
    int err = 0;

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARMA_ORDER) {
	err = 1;
    } 

    if (!err) {
	ainfo->p = list[1];
	ainfo->q = list[2];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 7) {
	    err = 1;
	} else if (list[4] < 0 || list[4] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} 
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo->P = list[4];
	ainfo->Q = list[5];
    }

    /* now that we have p and q we can check for masked lags */

    if (!err) {
	err = arma_make_masks(ainfo);
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(ainfo, dset);
	}
	if ((opt & OPT_N) || (armax && !hadconst)) {
	    ; /* no constant present */
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->nexo = list[0] - ypos;
	ainfo->nc = count_arma_coeffs(ainfo);
	ainfo->yno = list[ypos];
	if (ainfo->nexo > 0) {
	    err = arma_add_xlist(ainfo, ypos);
	}
    }

    return err;
}

static int check_arima_list (arma_info *ainfo,
			     const DATASET *dset,
			     gretlopt opt)
{
    int *list = ainfo->alist;
    int ypos = arma_has_seasonal(ainfo) ? 9 : 5;
    int armax = (list[0] > ypos);
    int hadconst = 0;
    int err = 0;

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARIMA_DIFF) {
	err = 1;
    } else if (list[3] < 0 || list[3] > MAX_ARMA_ORDER) {
	err = 1;
    } 

    if (!err) {
	ainfo->p = list[1];
	ainfo->d = list[2];
	ainfo->q = list[3];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 9) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[6] < 0 || list[6] > MAX_ARIMA_DIFF) {
	    err = 1;
	} else if (list[7] < 0 || list[7] > MAX_ARMA_ORDER) {
	    err = 1;
	} 
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo->P = list[5];
	ainfo->D = list[6];
	ainfo->Q = list[7];
    }

    /* now that we have p and q we can check for masked lags */

    if (!err) {
	err = arma_make_masks(ainfo);
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(ainfo, dset);
	}
	if ((opt & OPT_N) || (armax && !hadconst)) {
	    ;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->nexo = list[0] - ypos;
	ainfo->nc = count_arma_coeffs(ainfo);
	ainfo->yno = list[ypos];
	if (ainfo->nexo > 0) {
	    err = arma_add_xlist(ainfo, ypos);
	}
    }

    return err;
}

static int arma_check_list (arma_info *ainfo, 
			    const DATASET *dset,
			    gretlopt opt)
{
    int *list = ainfo->alist;
    int sep1 = gretl_list_separator_position(list);
    int err = 0;

    if (sep1 == 3) {
	if (list[0] < 4) {
	    err = E_PARSE;
	}
    } else if (sep1 == 4) {
	if (list[0] < 5) {
	    err = E_PARSE;
	} else {
	    set_arma_is_arima(ainfo);
	}
    } else {
	err = E_PARSE;
    }

    if (!err) {
	err = check_arma_sep(ainfo, sep1);
    }

    if (!err) {
	if (arma_is_arima(ainfo)) {
	    /* check for arima spec */
	    err = check_arima_list(ainfo, dset, opt);
	} else {	    
	    /* check for simple arma spec */
	    err = check_arma_list(ainfo, dset, opt);
	} 
    }

    /* catch null model */
    if (ainfo->nc == 0) {
	err = E_ARGS;
    }

    return err;
}

static void 
real_arima_difference_series (double *dx, const double *x,
			      int t1, int t2, int *delta, 
			      int k)
{
    int i, p, t, s = 0;
    
    for (t=t1; t<=t2; t++) {
	dx[s] = x[t];
	for (i=0; i<k && !na(dx[s]); i++) {
	    if (delta[i] != 0) {
		p = t - i - 1;
		if (p < 0 || na(x[p])) {
		    dx[s]  = NADBL;
		} else {
		    dx[s] -= delta[i] * x[p];
		}
	    }
	}
	s++;
    }
}

#ifndef X12A_CODE

/* Add to the ainfo struct a full-length series y holding 
   the differenced version of the dependent variable.
   If the "xdiff" flag is set on ainfo, in addition 
   create a matrix dX holding the differenced regressors;
   in that case the time-series length of dX depends on
   the @fullX flag -- if fullX = 0, this equals
   ainfo->T but if fullX = 0 it equals ainfo->t2 + 1.
*/

int arima_difference (arma_info *ainfo, const DATASET *dset, 
		      int fullX)
{
    const double *y = dset->Z[ainfo->yno];
    double *dy = NULL;
    int *delta = NULL;
    int s = ainfo->pd;
    int k, t, t1 = 0;
    int err = 0;

#if ARMA_DEBUG
    fprintf(stderr, "doing arima_difference: d = %d, D = %d\n",
	    ainfo->d, ainfo->D);
    fprintf(stderr, "ainfo->t1 = %d, ainfo->t2 = %d\n", ainfo->t1,
	    ainfo->t2);
#endif

    /* note: dy is a full length series (dset->n) */

    dy = malloc(dset->n * sizeof *dy);
    if (dy == NULL) {
	return E_ALLOC;
    }

    delta = arima_delta_coeffs(ainfo->d, ainfo->D, s);
    if (delta == NULL) {
	free(dy);
	return E_ALLOC;
    }    

    for (t=0; t<dset->n; t++) {
	dy[t] = NADBL;
    }

    for (t=0; t<dset->n; t++) {
	if (na(y[t])) {
	    t1++;
	} else {
	    break;
	}
    }

    t1 += ainfo->d + ainfo->D * s;
    k = ainfo->d + s * ainfo->D;

    real_arima_difference_series(dy + t1, y, t1, ainfo->t2, delta, k);

#if ARMA_DEBUG > 1
    for (t=0; t<dset->n; t++) {
	fprintf(stderr, "dy[%d] = % 12.7g\n", t, dy[t]);
    }
#endif    

    ainfo->y = dy;
    set_arima_ydiff(ainfo);

    if (arma_xdiff(ainfo)) {
	/* also difference the ARIMAX regressors */
	int xt1 = ainfo->t1, xT = ainfo->T;

	if (fullX) {
	    xt1 = 0;
	    xT = ainfo->t2 + 1;
	} 

	ainfo->dX = gretl_matrix_alloc(xT, ainfo->nexo);

	if (ainfo->dX == NULL) {
	    err = E_ALLOC;
	} else {
	    double *val = ainfo->dX->val;
	    int i, vi;

	    for (i=0; i<ainfo->nexo; i++) {
		vi = ainfo->xlist[i+1];
		real_arima_difference_series(val, dset->Z[vi], xt1, 
					     ainfo->t2, delta, k);
		val += xT;
	    }
	}
    }

    free(delta);

    return err;
}

void arima_difference_undo (arma_info *ainfo, const DATASET *dset)
{
    free(ainfo->y);
    ainfo->y = (double *) dset->Z[ainfo->yno];

    if (ainfo->dX != NULL) {
	gretl_matrix_free(ainfo->dX);
	ainfo->dX = NULL;
    }

    unset_arima_ydiff(ainfo);
}

#endif /* X12A_CODE not defined */
