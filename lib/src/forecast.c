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

static int 
allocate_fit_resid_arrays (FITRESID *fr, int n, int errs)
{
    fr->actual = malloc(n * sizeof *fr->actual);
    if (fr->actual == NULL) {
	return E_ALLOC;
    }

    fr->fitted = malloc(n * sizeof *fr->fitted);
    if (fr->fitted == NULL) {
	free(fr->actual);
	fr->actual = NULL;
	return E_ALLOC;
    }

    if (errs) {
	fr->sderr = malloc(n * sizeof *fr->sderr);
	if (fr->sderr == NULL) {
	    free(fr->actual);
	    fr->actual = NULL;
	    free(fr->fitted);
	    fr->fitted = NULL;
	    return E_ALLOC;
	}
    } else {
	fr->sderr = NULL;
    }

    return 0;
}

FITRESID *fit_resid_new (int n, int errs)
{
    FITRESID *fr = malloc(sizeof *fr);

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
    } else {
	if (allocate_fit_resid_arrays(fr, n, errs)) {
	    free(fr);
	    fr = NULL;
	} else {
	    fr->nobs = n;
	}
    }
    
    return fr;
}

void free_fit_resid (FITRESID *fr)
{
    free(fr->actual);
    free(fr->fitted);
    free(fr->sderr);
    free(fr);
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
				int *pt1, int *pt2)
{
    int i, t;
    int t1 = *pt1, t2 = *pt2;
    int imin = (pmod->ifc)? 3 : 2;
    int miss;

    for (t=t1; t<=t2; t++) {
	miss = 0;
	for (i=imin; i<=pmod->list[0]; i++) {
	    if (na(Z[pmod->list[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) t1++;
	else break;
    }

    for (t=t2; t>t1; t--) {
	miss = 0;
	for (i=imin; i<=pmod->list[0]; i++) {
	    if (na(Z[pmod->list[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) t2--;
	else break;
    }

    *pt1 = t1;
    *pt2 = t2;
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

/* initialize FITRESID struct based on string giving starting
   and ending observations */

static int 
fit_resid_init (const char *line, const MODEL *pmod, 
		const double **Z, const DATAINFO *pdinfo,
		FITRESID *fr)
{
    char t1str[OBSLEN], t2str[OBSLEN];

    /* parse the date strings giving the limits of the forecast
       range */
    if (sscanf(line, "%*s %10s %10s", t1str, t2str) != 2) {
	fr->err = E_OBS;
    }

    if (!fr->err) {
	fr->t1 = dateton(t1str, pdinfo);
	fr->t2 = dateton(t2str, pdinfo);

	if (fr->t1 < 0 || fr->t2 < 0 || fr->t2 <= fr->t1) {
	    fr->err = E_OBS;
	}
    }

    /* move endpoints if there are missing values for the
       independent variables */
    if (!fr->err) {
	fcast_adjust_t1_t2(pmod, Z, &fr->t1, &fr->t2);
	fr->nobs = fr->t2 - fr->t1 + 1;
	if (fr->nobs == 0) {
	    fr->err = E_OBS;
	}
    }

    if (!fr->err) {
	fr->err = allocate_fit_resid_arrays(fr, fr->nobs, 1);
    }

    return fr->err;
}

/* The following is based on Davidson and MacKinnon, Econometric
   Theory and Methods, chapter 3 (p. 104), which shows how the
   variance of forecast errors can be computed given the covariance
   matrix of the parameter estimates, provided the error term is
   assumed to be serially uncorrelated.
*/

FITRESID *get_fcast_with_errs (const char *str, MODEL *pmod, 
			       const double **Z, DATAINFO *pdinfo, 
			       PRN *prn)
{
    gretl_matrix *V = NULL;
    gretl_vector *Xs = NULL;
    gretl_vector *b = NULL;
    FITRESID *fr = NULL;

    double s2 = pmod->sigma * pmod->sigma;
    int yno = pmod->list[1];
    int k = pmod->ncoeff;
    int i, vi, s, t;

    fr = fit_resid_new(0, 1); 
    if (fr == NULL) {
	return NULL;
    }

    if (AR_MODEL(pmod->ci)) {
	/* need to test this this works for "exotic" non-AR models */
	fr->err = E_NOTIMP;
	return fr;
    }

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    fit_resid_init(str, pmod, Z, pdinfo, fr);
    if (fr->err) {
	return fr;
    }    
    
    V = gretl_vcv_matrix_from_model(pmod, NULL);
    if (V == NULL) {
	fr->err = E_ALLOC;
	goto bailout;
    }

    Xs = gretl_vector_alloc(k);
    if (Xs == NULL) {
	fr->err = E_ALLOC;
	goto bailout;
    }

    b = gretl_coeff_vector_from_model(pmod, NULL);
    if (b == NULL) {
	fr->err = E_ALLOC;
	goto bailout;
    }

    for (t=fr->t1; t<=fr->t2 && !fr->err; t++) {
	int missing = 0;
	double vyh;

	s = t - fr->t1;

	fr->actual[s] = Z[yno][t];

	/* skip if we can't compute forecast */
	if (t >= pmod->t1 && t <= pmod->t2) {
	    missing = na(pmod->yhat[t]);
	} else {
	    missing = fcast_x_missing(pmod->list, Z, t);
	}
	if (missing) {
	    fr->sderr[s] = fr->fitted[s] = NADBL;
	    continue;
	}	    

	/* populate Xs vector for observation */
	for (i=0; i<k; i++) {
	    vi = pmod->list[i + 2];
	    gretl_vector_set(Xs, i, Z[vi][t]);
	}

	/* forecast value */
	fr->fitted[s] = gretl_matrix_dot_product(Xs, GRETL_MOD_NONE,
						 b, GRETL_MOD_NONE,
						 NULL);

	/* forecast variance */
	vyh = gretl_scalar_b_X_b_prime(Xs, V, NULL);
	if (na(vyh)) {
	    fr->err = 1;
	} else {
	    vyh += s2;
	    if (vyh >= 0.0) {
		fr->sderr[s] = sqrt(vyh);
	    } else {
		fr->err = 1;
	    }
	}
    }

    fr->tval = tcrit95(pmod->dfd);
    strcpy(fr->depvar, pdinfo->varname[yno]);
    fr->df = pmod->dfd;

 bailout:

    gretl_matrix_free(V);
    gretl_vector_free(Xs);
    gretl_vector_free(b);

    return fr;
}

int fcast_with_errs (const char *str, MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    FITRESID *fr;
    int err;

    fr = get_fcast_with_errs(str, pmod, (const double **) *pZ, 
			     pdinfo, prn);

    if (fr == NULL) {
	return E_ALLOC;
    }

    if ((err = fr->err) == 0) {
	err = text_print_fcast_with_errs(fr, pZ, pdinfo, opt, prn);
    }

    free_fit_resid(fr);
    
    return err;
}

/* Compute forecasts for various sorts of gretl models, including
   those with autoregressive errors.

   The "ar" code below generates one-step ahead forecasts that
   incorporate the predictable portion of an AR error term:

       u_t = r1 u_{t-1} + r2 u_{t-1} + ... + e_t

   where e_t is white noise.  The forecasts are based on the
   representation of a model with such an error term as

       (1 - r(L)) y_t = (1 - r(L) X_t b + e_t

   where r(L) is a polynomial in the lag operator.  In effect, we
   generate a forecast that incorporates the error process by using
   rho-differenced X on the RHS, and applying a corresponding
   compensation on the LHS so that the forecast is of y_t, not (1 -
   r(L)) y_t.
*/

static int gretl_forecast (int t1, int t2, int nv, 
			   const MODEL *pmod, double ***pZ)
{
    double xval, yh;
    int i, k, maxlag = 0, yno;
    int v, t, miss;
    const int *arlist = NULL;
    double *yhat = (*pZ)[nv];
    int ar = SIMPLE_AR_MODEL(pmod->ci);

    /* bodge: for now we're not going to forecast out of sample
       for these estimators. TODO */
    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == GARCH) {
	for (t=t1; t<=t2; t++) {
	    yhat[t] = pmod->yhat[t];
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
	yh = 0.0;

	if (pmod->ci == PWE && t == pmod->t1) {
	    /* PWE first obs is special */
	    yhat[t] = pmod->yhat[t];
	    continue;
	}

	if (ar) { 
	    /* LHS adjustment */
	    double rk, ylag;

	    for (k=1; k<=arlist[0]; k++) {
		rk = pmod->arinfo->rho[k];
		/* use actual lagged y by preference */
		ylag = (*pZ)[yno][t - arlist[k]];
		if (na(ylag)) {
		    /* use forecast of lagged y */
		    ylag = yhat[t - arlist[k]];
		}
		if (na(ylag)) {
		    yhat[t] = NADBL;
		    miss = 1;
		} else {
		    yh += rk * ylag;
		}
	    }
	} /* end if ar */

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    v = pmod->list[i+2];
	    xval = (*pZ)[v][t];
	    if (na(xval)) {
		miss = 1;
	    } else {
		if (ar) {
		    /* use rho-differenced X on RHS */
		    double rk, xlag;

		    for (k=1; k<=arlist[0]; k++) {
			rk = pmod->arinfo->rho[k];
			xlag = (*pZ)[v][t - arlist[k]];
			if (!na(xlag)) {
			    xval -= rk * xlag;
			}
		    }
		}
		yh += xval * pmod->coeff[i];
	    }
	}

	if (!miss && pmod->ci == LOGISTIC) {
	    double lmax = gretl_model_get_double(pmod, "lmax");

	    yh = lmax / (1.0 + exp(-yh));
	}

	if (miss) {
	    yhat[t] = NADBL;
	} else { 
	    yhat[t] = yh;
	}
#if FCAST_DEBUG
	fprintf(stderr, "gretl_forecast: set yhat[%d] = %g\n", t, yhat[t]);
#endif
    }

    return 0;
}

/**
 * fcast:
 * @line: command line, giving a starting observation, ending
 * observation, and variable name to use for the forecast values
 * (the starting and ending observations may be omitted).
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: pointer to data information struct.
 *
 * Adds to the dataset a new variable containing predicted values for the
 * dependent variable in @pmod over the specified range of observations,
 * or, by default, over the sample range currently defined in @pdinfo.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int fcast (const char *line, const MODEL *pmod, double ***pZ,
	   DATAINFO *pdinfo)
{
    int t, t1, t2, vi;
    char t1str[OBSLEN], t2str[OBSLEN], varname[VNAMELEN];
    int err = 0;

    *t1str = '\0'; *t2str = '\0';

    /* the varname should either be in the 2nd or 4th position */
    if (sscanf(line, "%*s %8s %8s %8s", t1str, t2str, varname) != 3) {
	if (sscanf(line, "%*s" "%8s", varname) != 1) {
	    return E_PARSE;
	}
    }

    if (*t1str && *t2str) {
	t1 = dateton(t1str, pdinfo);
	t2 = dateton(t2str, pdinfo);
	if (t1 < 0 || t2 < 0 || t2 < t1) {
	    return E_DATA;
	}
    } else {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    }

    if (check_varname(varname)) {
	return 1;
    }

    vi = varindex(pdinfo, varname);

    if (vi == pdinfo->v) {
	err = dataset_add_series(1, pZ, pdinfo);
    }

    if (!err) {
	strcpy(pdinfo->varname[vi], varname);
	strcpy(VARLABEL(pdinfo, vi), _("predicted values"));

	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[vi][t] = NADBL;
	}

	gretl_forecast(t1, t2, vi, pmod, pZ);
    }

    return err;
}


