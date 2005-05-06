/*
 * Copyright (C) 1999-2005 Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
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

int gretl_forecast (int t1, int t2, int nv, 
		    const MODEL *pmod, double ***pZ)
{
    double xx, zz, zr;
    int i, k, maxlag = 0, yno;
    int v, t, miss;
    const int *arlist = NULL;
    int ar = AR_MODEL(pmod->ci);

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == GARCH) {
	for (t=t1; t<=t2; t++) {
	    (*pZ)[nv][t] = pmod->yhat[t];
	}
	return 0;
    }

    yno = pmod->list[1];

    if (ar) {
	arlist = pmod->arinfo->arlist;
	maxlag = arlist[arlist[0]];
	if (t1 < maxlag) {
	    t1 = maxlag; 
	}
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	zz = 0.0;
	if (ar) 
	    for (k=1; k<=arlist[0]; k++) {
	    xx = (*pZ)[yno][t - arlist[k]];
	    zr = pmod->arinfo->rho[k];
	    if (na(xx)) {
		if (zr == 0) {
		    continue;
		}
		xx = (*pZ)[nv][t - arlist[k]];
		if (na(xx)) {
		    (*pZ)[nv][t] = NADBL;
		    miss = 1;
		}
	    }
	    zz = zz + xx * zr;
	} /* end if ar */
	for (v=0; !miss && v<pmod->ncoeff; v++) {
	    k = pmod->list[v+2];
	    xx = (*pZ)[k][t];
	    if (na(xx)) {
		zz = NADBL;
		miss = 1;
	    }
	    if (!miss && ar) {
		xx = (*pZ)[k][t];
		for (i=1; i<=arlist[0]; i++) {
		    xx -= pmod->arinfo->rho[i] * (*pZ)[k][t - arlist[i]];
		}
	    }
	    if (!miss) {
		zz = zz + xx * pmod->coeff[v];
	    }
	}

	if (pmod->ci == LOGISTIC) {
	    double lmax = gretl_model_get_double(pmod, "lmax");

	    zz = lmax / (1.0 + exp(-zz));
	}

	(*pZ)[nv][t] = zz;
    }

    return 0;
}

static int allocate_fit_resid_arrays (FITRESID *fr, int n, int errs)
{
    fr->actual = malloc(n * sizeof *fr->actual);
    if (fr->actual == NULL) {
	return 1;
    }

    fr->fitted = malloc(n * sizeof *fr->fitted);
    if (fr->fitted == NULL) {
	free(fr->actual);
	fr->actual = NULL;
	return 1;
    }

    if (errs) {
	fr->sderr = malloc(n * sizeof *fr->sderr);
	if (fr->sderr == NULL) {
	    free(fr->actual);
	    fr->actual = NULL;
	    free(fr->fitted);
	    fr->fitted = NULL;
	    return 1;
	}
    } else {
	fr->sderr = NULL;
    }

    return 0;
}

FITRESID *fit_resid_new (int n, int errs)
{
    FITRESID *fr;

    fr = malloc(sizeof *fr);
    if (fr == NULL) {
	return NULL;
    }

    fr->err = 0;
    fr->t1 = 0;
    fr->t2 = 0;
    fr->nobs = 0;
    fr->real_nobs = 0;

    if (n == 0) {
	fr->actual = NULL;
	fr->fitted = NULL;
	fr->sderr = NULL;
	return fr;
    }

    if (allocate_fit_resid_arrays(fr, n, errs)) {
	free(fr);
	return NULL;
    }

    fr->nobs = n;
    
    return fr;
}

FITRESID *get_fit_resid (const MODEL *pmod, double ***pZ, 
			 DATAINFO *pdinfo)
{
    int depvar, t, ft;
    FITRESID *fr;

    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	depvar = pmod->list[4];
    } else {
	depvar = pmod->list[1];
    }

    fr = fit_resid_new(pmod->t2 - pmod->t1 + 1, 0);
    if (fr == NULL) {
	return NULL;
    }

    fr->sigma = pmod->sigma;

    fr->t1 = pmod->t1;
    fr->t2 = pmod->t2;
    fr->real_nobs = pmod->nobs;

    for (t=fr->t1; t<=fr->t2; t++) {
	ft = t - fr->t1;
	fr->actual[ft] = (*pZ)[depvar][t];
	fr->fitted[ft] = pmod->yhat[t];
    }

    if (gretl_isdummy(0, fr->nobs, fr->actual) > 0) {
	fr->pmax = get_precision(fr->fitted, fr->nobs, 8);
    } else {
	fr->pmax = get_precision(fr->actual, fr->nobs, 8);
    }
    
    strcpy(fr->depvar, pdinfo->varname[depvar]);
    
    return fr;
}

static void fcast_adjust_t1_t2 (const MODEL *pmod,
				const double **Z,
				int *t1, int *t2)
{
    int i, t;
    int my_t1 = *t1, my_t2 = *t2;
    int imin = (pmod->ifc)? 3 : 2;
    int miss;

    for (t=*t1; t<=*t2; t++) {
	miss = 0;
	for (i=imin; i<=pmod->list[0]; i++) {
	    if (na(Z[pmod->list[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) my_t1++;
	else break;
    }

    for (t=*t2; t>0; t--) {
	miss = 0;
	for (i=imin; i<=pmod->list[0]; i++) {
	    if (na(Z[pmod->list[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) my_t2--;
	else break;
    }

    *t1 = my_t1;
    *t2 = my_t2;
}

static int 
fcast_x_missing (const int *list, const double **Z, int t)
{
    int i, ret = 0;

    for (i=2; i<=list[0]; i++) {
	if (na(Z[list[i]][t])) {
	    ret = 1;
	    break;
	}
    }

    return ret;
}

/* 
   Below: the method for generating forecasts and prediction errors
   that is presented in Wooldridge's Introductory Econometrics,
   Chapter 6.
*/

FITRESID *get_fcast_with_errs (const char *str, const MODEL *pmod, 
			       double ***pZ, DATAINFO *pdinfo, 
			       PRN *prn)
{
    double **fZ = NULL;
    DATAINFO *finfo = NULL;
    MODEL fmod; 
    FITRESID *fr;
    int *flist = NULL;
    int i, k, t, n_est, nv;
    int yno = pmod->list[1];
    char t1str[OBSLEN], t2str[OBSLEN];

    fr = fit_resid_new(0, 1); 
    if (fr == NULL) {
	return NULL;
    }

    if (pmod->ci != OLS) {
	fr->err = E_OLSONLY;
	return fr;
    }

    /* bodge (reject in case of subsampled data) */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    /* parse dates */
    if (sscanf(str, "%*s %10s %10s", t1str, t2str) != 2) {
	fr->err = E_OBS;
	return fr;
    }

    fr->t1 = dateton(t1str, pdinfo);
    fr->t2 = dateton(t2str, pdinfo);

    if (fr->t1 < 0 || fr->t2 < 0 || fr->t2 <= fr->t1) {
	fr->err = E_OBS;
	return fr;
    }

    /* move endpoints if there are missing vals for the
       independent variables */
    fcast_adjust_t1_t2(pmod, (const double **) *pZ, &fr->t1, &fr->t2);

    /* number of obs for which forecasts will be generated */
    fr->nobs = fr->t2 - fr->t1 + 1;

    if (allocate_fit_resid_arrays(fr, fr->nobs, 1)) {
	fr->err = E_ALLOC;
	return fr;
    }

    nv = pmod->list[0];
    if (!pmod->ifc) nv++;

    n_est = pmod->t2 - pmod->t1 + 1;

    finfo = create_new_dataset(&fZ, nv, n_est, 0);
    if (finfo == NULL) {
	fr->err = E_ALLOC;
	return fr;
    }

    /* insert depvar at position 1 */
    for (t=0; t<finfo->n; t++) {
	fZ[1][t] = (*pZ)[yno][t + pmod->t1];
    }

    /* create new list */
    flist = malloc((finfo->v + 1) * sizeof *flist);
    if (flist == NULL) {
	fr->err = E_ALLOC;
	goto fcast_bailout;
    }

    flist[0] = finfo->v;
    flist[1] = 1;
    flist[2] = 0;
    for (i=3; i<=flist[0]; i++) {
	flist[i] = i - 1;
    }

    gretl_model_init(&fmod);

    /* loop across the observations for which we want forecasts
       and standard errors */

#ifdef FCAST_DEBUG
    printf("get_fcast_with_errs: ft1=%d, ft2=%d, pmod->t1=%d, pmod->t2=%d\n",
	   fr->t1, fr->t2, pmod->t1, pmod->t2);
#endif

    for (k=0; k<fr->nobs; k++) {
	int tk = k + fr->t1;

	fr->actual[k] = (*pZ)[yno][tk];

	if (fcast_x_missing(pmod->list, (const double **) *pZ, tk)) {
	    fr->sderr[k] = fr->fitted[k] = NADBL;
	    continue;
	}

	/* form modified indep vars: original data minus the values
	   to be used for the forecast 
	*/
	for (i=3; i<=flist[0]; i++) {
	    int v = (pmod->ifc)? pmod->list[i] : pmod->list[i-1];
	    const double *xv = (*pZ)[v];

	    for (t=0; t<finfo->n; t++) {
		int tp = t + pmod->t1;

		if (na(xv[tp])) {
		    fZ[i-1][t] = NADBL;
		} else {
		    fZ[i-1][t] = xv[tp] - xv[tk];
		}
	    }
	}

	fmod = lsq(flist, &fZ, finfo, OLS, OPT_A, 0.0);

	if (fmod.errcode) {
	    fr->err = fmod.errcode;
	    clear_model(&fmod);
	    goto fcast_bailout;
	}

	fr->fitted[k] = fmod.coeff[0];

	/* what exactly do we want here? */
#ifdef GIVE_SDERR_OF_EXPECTED_Y
	fr->sderr[k] = fmod.sderr[0];
#else
	fr->sderr[k] = sqrt(fmod.sderr[0] * fmod.sderr[0] + 
			    fmod.sigma * fmod.sigma);
#endif
	clear_model(&fmod);
    }

    fr->tval = tcrit95(pmod->dfd);
    strcpy(fr->depvar, pdinfo->varname[yno]);
    fr->df = pmod->dfd;

 fcast_bailout:

    free_Z(fZ, finfo);
    free(flist);
    clear_datainfo(finfo, CLEAR_FULL);
    free(finfo);

    return fr;
}

void free_fit_resid (FITRESID *fr)
{
    free(fr->actual);
    free(fr->fitted);
    free(fr->sderr);
    free(fr);
}

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn,
		     int plot)
{
    FITRESID *fr;
    int err;

    fr = get_fcast_with_errs(str, pmod, pZ, pdinfo, prn);

    if (fr == NULL) {
	return E_ALLOC;
    }

    if ((err = fr->err) == 0) {
	err = text_print_fcast_with_errs(fr, pZ, pdinfo, prn, plot);
    }

    free_fit_resid(fr);
    
    return err;
}

