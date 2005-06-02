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
#include "forecast.h"

/* estimators where a simple X*b does _not_ give the
   predicted value of the dependent variable */

#define FCAST_SPECIAL(c) (c == LOGIT || \
                          c == LOGISTIC || \
                          c == NLS || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == TOBIT)

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

/**
 * fit_resid_new:
 * @n: the number of observations to allow for, or 0 if this
 * information will be added later.
 * @errs: set = 1 to allocate space in the structure for
 * forecast standard errors, otherwise 0. 
 *
 * Allocates a #FITRESID struct for holding fitted values and
 * residuals from a model (or out-of-sample forecasts).  If
 * @n is greater than 0 the arrays required for that number
 * of observations will be allocated.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

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

/**
 * free_fit_resid:
 * @fr: the pointer to be freed.
 *
 * Frees all resources associated with @fr, then frees the pointer
 * itself.
 */

void free_fit_resid (FITRESID *fr)
{
    free(fr->actual);
    free(fr->fitted);
    free(fr->sderr);
    free(fr);
}

/**
 * get_fit_resid:
 * @pmod: the model for which actual and fitted values
 * are wanted.
 * @Z: data array using which @pmod was estimated.
 * @pdinfo: dataset information.
 *
 * Allocates a #FITRESID structure and fills it out with
 * the actual and predicted values of the dependent variable
 * in @pmod.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FITRESID *get_fit_resid (const MODEL *pmod, const double **Z, 
			 const DATAINFO *pdinfo)
{
    int depvar, t, ft;
    FITRESID *fr;

    depvar = gretl_model_get_depvar(pmod);

    fr = fit_resid_new(pmod->t2 - pmod->t1 + 1, 0);
    if (fr == NULL) {
	return NULL;
    }

    if (LIMDEP(pmod->ci)) {
	fr->sigma = NADBL;
    } else {
	fr->sigma = pmod->sigma;
    }

    fr->t1 = pmod->t1;
    fr->t2 = pmod->t2;
    fr->real_nobs = pmod->nobs;

    for (t=fr->t1; t<=fr->t2; t++) {
	ft = t - fr->t1;
	fr->actual[ft] = Z[depvar][t];
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

static int fcast_adjust_t1_t2 (const MODEL *pmod,
			       const double **Z,
			       int *pt1, int *pt2)
{
    int *xlist = NULL;
    int t1 = *pt1, t2 = *pt2;
    int imin, miss;
    int i, t;

    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	xlist = gretl_arma_model_get_x_list(pmod);
	imin = (pmod->ifc)? 2 : 1;
    } else {
	xlist = gretl_list_copy(pmod->list);
	imin = (pmod->ifc)? 3 : 2;
    }

    if (xlist == NULL) {
	return 1;
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	for (i=imin; i<=xlist[0]; i++) {
	    if (na(Z[xlist[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) t1++;
	else break;
    }

    for (t=t2; t>t1; t--) {
	miss = 0;
	for (i=imin; i<=xlist[0]; i++) {
	    if (na(Z[xlist[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) t2--;
	else break;
    }

    *pt1 = t1;
    *pt2 = t2;

    free(xlist);

    return 0;
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
		FITRESID *fr, int errs)
{
    char t1str[OBSLEN], t2str[OBSLEN];

    if (!strncmp(line, "fcasterr", 8)) {
	line += 9;
    }

    /* parse the date strings giving the limits of the forecast
       range */
    if (sscanf(line, "%10s %10s", t1str, t2str) != 2) {
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
	fr->err = allocate_fit_resid_arrays(fr, fr->nobs, errs);
    }

    return fr->err;
}

/* Get forecasts plus standard errors for same, for models without
   autoregressive errors and without "special requirements"
   (e.g. nonlinearity).  The forecast standard errors include both
   uncertainty over the error process and parameter uncertainty
   (Davidson and MacKinnon method).
*/

static int
get_static_fcast_with_errs (FITRESID *fr, MODEL *pmod, 
			    const double **Z, const DATAINFO *pdinfo) 
{
    gretl_matrix *V = NULL;
    gretl_vector *Xs = NULL;
    gretl_vector *b = NULL;

    double s2 = pmod->sigma * pmod->sigma;
    int yno = gretl_model_get_depvar(pmod);
    int k = pmod->ncoeff;
    int i, vi, s, t;
    int err = 0;

    V = gretl_vcv_matrix_from_model(pmod, NULL);
    if (V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    Xs = gretl_vector_alloc(k);
    if (Xs == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    b = gretl_coeff_vector_from_model(pmod, NULL);
    if (b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (t=fr->t1; t<=fr->t2 && !err; t++) {
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
	    err = 1;
	} else {
	    vyh += s2;
	    if (vyh >= 0.0) {
		fr->sderr[s] = sqrt(vyh);
	    } else {
		fr->sderr[s] = NADBL;
		err = 1;
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

    return err;
}

/* generate forecasts from nonlinear least squares model, using
   the string specification of the regression function that
   was saved as data on the model (see nls.c) 
*/

static int nls_fcast (int t1, int t2, double *fcast, const MODEL *pmod, 
		      double ***pZ, DATAINFO *pdinfo)
{
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    const char *nlfunc;    
    char formula[MAXLINE];
    int t, oldv = pdinfo->v;
    int err = 0;

    nlfunc = gretl_model_get_data(pmod, "nl_regfunc");
    if (nlfunc == NULL) {
	err = E_DATA;
    }

    if (!err) {
	pdinfo->t1 = t1;
	pdinfo->t2 = t2;
	sprintf(formula, "$nl_y = %s", nlfunc);
	err = generate(formula, pZ, pdinfo, NULL, OPT_P);
    }

    if (!err) {
	/* transcribe values from last generated var to target */
	for (t=t1; t<=t2; t++) {
	    fcast[t] = (*pZ)[oldv][t];
	}
	err = dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo);
    }

    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    return err;
}

#define ARF_DEBUG 0

#if ARF_DEBUG
# include <stdarg.h>
static void dprintf (const char *format, ...)
{
   va_list args;

   va_start(args, format);
   vfprintf(stderr, format, args);
   va_end(args);

   return;
}
# define DPRINTF(x) dprintf x
#else 
# define DPRINTF(x)
#endif 

static int 
garch_fcast (int t1, int t2, double *yhat, const MODEL *pmod, 
	     const double **Z, const DATAINFO *pdinfo)
{
    double xval;
    int xvars, yno;
    int *xlist = NULL;
    int i, t;

    /* use pre-calculated fitted values over model estimation range */
    for (t=t1; t<=pmod->t2; t++) {
	yhat[t] = pmod->yhat[t];
    }

    if (t2 <= pmod->t2) {
	/* no real forecasts called for */
	return 0;
    }

    xlist = gretl_arma_model_get_x_list(pmod);
    yno = gretl_model_get_depvar(pmod);

    if (xlist != NULL) {
	xvars = xlist[0];
    } else {
	xvars = 0;
    }

    /* FIXME below: we could do something interesting for
       the forecast error variance */

    /* do real forecast */
    for (t=pmod->t2 + 1; t<=t2; t++) {
	int miss = 0;
	double yh = 0.0;

	for (i=1; i<=xvars; i++) {
	    xval = Z[xlist[i]][t];
	    if (na(xval)) {
		miss = 1;
	    } else {
		yh += pmod->coeff[i-1] * xval;
	    }
	}

	if (miss) {
	    yhat[t] = NADBL;
	} else {
	    yhat[t] = yh;
	}
    }

    free(xlist);

    return 0;
}

/* generate forecasts for ARMA (or ARMAX) models */

static int 
arma_fcast (int t1, int t2, double *yhat, const MODEL *pmod, 
	    const double **Z, const DATAINFO *pdinfo)
{
    const double *ar_coeff;
    const double *ma_coeff;

    double xval, yval;
    int xvars, yno;
    int *xlist = NULL;
    int p, q;
    int i, s, t;

    /* use pre-calculated fitted values over model estimation range */
    for (t=t1; t<=pmod->t2; t++) {
	yhat[t] = pmod->yhat[t];
    }

    if (t2 <= pmod->t2) {
	/* no real forecasts called for */
	return 0;
    }

    p = gretl_arma_model_get_AR_order(pmod);
    q = gretl_arma_model_get_MA_order(pmod);
    xlist = gretl_arma_model_get_x_list(pmod);
    yno = gretl_model_get_depvar(pmod);

    DPRINTF(("forecasting variable %d (%s)\n", yno, pdinfo->varname[yno]));

    if (xlist != NULL) {
	xvars = xlist[0];
    } else {
	xvars = 0;
    }

    ar_coeff = pmod->coeff + xvars;
    ma_coeff = ar_coeff + p;

    /* do real forecast */
    for (t=pmod->t2 + 1; t<=t2; t++) {
	int miss = 0;
	double yh = 0.0;

	DPRINTF(("\n *** Doing forecast for obs %d\n", t));

	/* contribution of independent variables */
	for (i=1; i<=xvars; i++) {
	    xval = Z[xlist[i]][t];
	    if (na(xval)) {
		miss = 1;
	    } else {
		yh += pmod->coeff[i-1] * xval;
	    }
	}

	DPRINTF((" x contribution = %g\n", yh));

	/* AR contribution */
	for (i=0; i<p && !miss; i++) {
	    s = t - i - 1;
	    if (s < 0) {
		yval = NADBL;
	    } else if (s <= pmod->t2) {
		yval = Z[yno][s];
		DPRINTF(("  AR: lag %d, y[%d] = %g\n", i+1, s, yval));
	    } else {
		yval = yhat[s];
		DPRINTF(("  AR: lag %d, yhat[%d] = %g\n", i+1, s, yval));
	    }
	    if (na(yval)) {
		miss = 1;
	    } else {
		DPRINTF(("  AR: lag %d, s=%d, using coeff %g\n", i+1, s, ar_coeff[i]));
		yh += ar_coeff[i] * yval;
	    }
	}

	DPRINTF((" with AR contribution: %g\n", yh));

	/* MA contribution */
	for (i=0; i<q && !miss; i++) {
	    s = t - i - 1;
	    if (s >= pmod->t1 && s <= pmod->t2) {
		DPRINTF(("  MA: lag %d, e[%d] = %g, coeff %g\n", i+1, s, 
			 pmod->uhat[s], ma_coeff[i]));
		yh += ma_coeff[i] * pmod->uhat[s];
	    }
	}

	DPRINTF((" with MA contribution: %g\n", yh));

	if (miss) {
	    yhat[t] = NADBL;
	} else {
	    yhat[t] = yh;
	}
    }

    free(xlist);

    return 0;
}

/* 
   The code below generates one-step ahead forecasts that
   incorporate the predictable portion of an AR error term:

       u_t = r1 u_{t-1} + r2 u_{t-1} + ... + e_t

   where e_t is white noise.  The forecasts are based on the
   representation of a model with such an error term as

       (1 - r(L)) y_t = (1 - r(L) X_t b + e_t

   where r(L) is a polynomial in the lag operator.  In effect, we
   generate a forecast that incorporates the error process by using
   rho-differenced X on the RHS, and applying a corresponding
   compensation on the LHS so that the forecast is of y_t, not 
   (1 - r(L)) y_t.
*/

static int 
ar_fcast (int t1, int t2, double *yhat, const MODEL *pmod, 
	  const double **Z, const DATAINFO *pdinfo)
{
    const int *arlist;
    double xval, yh;
    double rk, ylag, xlag;
    int miss, maxlag, yno;
    int i, k, v, t;

    yno = pmod->list[1];
    arlist = pmod->arinfo->arlist;
    maxlag = arlist[arlist[0]];

    if (t1 < maxlag) {
	t1 = maxlag; 
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	yh = 0.0;

	if (pmod->ci == PWE && t == pmod->t1) {
	    /* PWE first obs is special */
	    yhat[t] = pmod->yhat[t];
	    continue;
	}

	/* LHS adjustment */
	for (k=1; k<=arlist[0]; k++) {
	    rk = pmod->arinfo->rho[k];
	    /* use actual lagged y by preference */
	    ylag = Z[yno][t - arlist[k]];
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

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    v = pmod->list[i+2];
	    xval = Z[v][t];
	    if (na(xval)) {
		miss = 1;
	    } else {
		/* use rho-differenced X on RHS */
		for (k=1; k<=arlist[0]; k++) {
		    rk = pmod->arinfo->rho[k];
		    xlag = Z[v][t - arlist[k]];
		    if (!na(xlag)) {
			xval -= rk * xlag;
		    }
		}
		yh += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    yhat[t] = NADBL;
	} else {
	    yhat[t] = yh;
	}
    }

    return 0;
}

/* Calculates the transformation required to get from xb (= X*b) to
   the actual prediction for the dependent variable, for models of
   type LOGISTIC, LOGIT, PROBIT, TOBIT and POISSON.
 */

static double fcast_transform (double xb, int ci, int t, 
			       const double *offset,
			       double lmax)
{
    double yf = xb;

    if (ci == TOBIT) {
	if (xb < 0.0) {
	    yf = 0.0;
	}
    } else if (ci == LOGIT) {
	yf = exp(xb) / (1.0 + exp(xb));
    } else if (ci == PROBIT) {
	yf = normal_cdf(xb);
    } else if (ci == LOGISTIC) {
	if (na(lmax)) {
	    yf = 1.0 / (1.0 + exp(-xb));
	} else {
	    yf = lmax / (1.0 + exp(-xb));
	}
    } else if (ci == POISSON) {
	if (offset != NULL) {
	    if (na(offset[t])) {
		yf = NADBL;
	    } else {
		yf = exp(xb + log(offset[t]));
	    }
	} else {
	    yf = exp(xb);
	}
    }

    return yf;
}

/* compute forecasts for static linear models */

static int 
linear_fcast (int t1, int t2, double *yhat, const MODEL *pmod, 
	      const double **Z, const DATAINFO *pdinfo)
{
    const double *offset = NULL;
    double lmax = NADBL;
    double xval;
    int i, t;

    if (pmod->ci == POISSON) {
	/* special for poisson "offset" variable */
	int offnum = gretl_model_get_int(pmod, "offset_var");

	if (offnum > 0) {
	    offset = Z[offnum];
	}
    } else if (pmod->ci == LOGISTIC) {
	lmax = gretl_model_get_double(pmod, "lmax");
    }

    for (t=t1; t<=t2; t++) {
	int miss = 0;
	double yh = 0.0;

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    xval = Z[pmod->list[i+2]][t];
	    if (na(xval)) {
		miss = 1;
	    } else {
		yh += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    yhat[t] = NADBL;
	} else if (FCAST_SPECIAL(pmod->ci)) {
	    /* special handling for LOGIT and others */
	    yhat[t] = fcast_transform(yh, pmod->ci, t, offset, lmax);
	} else {
	    yhat[t] = yh;
	}
    }

    return 0;
}

/* driver for various functions that compute forecasts
   for different sorts of models */

static int
get_fcast_without_errs (FITRESID *fr, MODEL *pmod, 
			double ***pZ, DATAINFO *pdinfo) 
{
    int yno = gretl_model_get_depvar(pmod);
    double *yhat;
    int s, t;
    int err = 0;

    yhat = malloc(pdinfo->n * sizeof *yhat);
    if (yhat == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	yhat[t] = NADBL;
    }

    if (pmod->ci == NLS) {
	err = nls_fcast(fr->t1, fr->t2, yhat, pmod, pZ, pdinfo);
    } else if (SIMPLE_AR_MODEL(pmod->ci)) {
	err = ar_fcast(fr->t1, fr->t2, yhat, pmod, 
		       (const double **) *pZ, pdinfo);
    } else if (pmod->ci == ARMA) {
	err = arma_fcast(fr->t1, fr->t2, yhat, pmod, 
			 (const double **) *pZ, pdinfo);
    } else if (pmod->ci == GARCH) {
	err = garch_fcast(fr->t1, fr->t2, yhat, pmod, 
			  (const double **) *pZ, pdinfo);
    } else {
	err = linear_fcast(fr->t1, fr->t2, yhat, pmod, 
			   (const double **) *pZ, pdinfo);
    }

    if (err) {
	free(yhat);
	return err;
    }

    for (t=fr->t1; t<=fr->t2; t++) {
	s = t - fr->t1;
	fr->actual[s] = (*pZ)[yno][t];
	fr->fitted[s] = yhat[t];
    }

    free(yhat);

    fr->tval = tcrit95(pmod->dfd);
    strcpy(fr->depvar, pdinfo->varname[yno]);
    fr->df = pmod->dfd;

    return err;
}

/* public forecast-related functions follow */

/**
 * get_forecast:
 * @str: string giving starting and ending observations, separated
 * by a space.
 * @pmod: the model from which forecasts are wanted.
 * @pZ: pointer to data array using which @pmod was estimated.
 * @pdinfo: dataset information.
 *
 * Allocates a #FITRESID structure and fills it out with forecasts
 * based on @pmod, over the specified range of observations.  
 * For some sorts of models forecast standard errors are also
 * computed (these appear in the %sderr member of the structure
 * to which a pointer is returned; otherwise the %sderr member is
 * %NULL).
 * 
 * The calculation of forecast errors, where applicable, is based 
 * on Davidson and MacKinnon, Econometric Theory and Methods, 
 * chapter 3 (p. 104), which shows how the variance of forecast errors 
 * can be computed given the covariance matrix of the parameter 
 * estimates, provided the error term may be assumed to be serially 
 * uncorrelated.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 * The %err member of the returned object should be checked:
 * a non-zero value indicates an error condition.
 */

FITRESID *get_forecast (const char *str, MODEL *pmod, 
			double ***pZ, DATAINFO *pdinfo) 
{
    FITRESID *fr;
    int do_errs = 1;

    if (AR_MODEL(pmod->ci) || FCAST_SPECIAL(pmod->ci)) {
	do_errs = 0;
    }

    fr = fit_resid_new(0, do_errs); 
    if (fr == NULL) {
	return NULL;
    }

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    fit_resid_init(str, pmod, (const double **) *pZ, pdinfo, fr, do_errs);
    if (fr->err) {
	return fr;
    }     
    
    if (do_errs) {
	fr->err = get_static_fcast_with_errs(fr, pmod, (const double **) *pZ, 
					     pdinfo);
    } else {
	fr->err = get_fcast_without_errs(fr, pmod, pZ, pdinfo);
    }

    return fr;
}

/**
 * add_forecast:
 * @str: command string, giving a starting observation, ending
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
 * In the case of "simple" models with an autoregressive error term 
 * (%AR, %CORC, %HILU, %PWE) the predicted values incorporate
 * the forecastable portion of the error.  
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int add_forecast (const char *str, const MODEL *pmod, double ***pZ,
		  DATAINFO *pdinfo)
{
    int t, t1, t2, vi;
    char t1str[OBSLEN], t2str[OBSLEN], varname[VNAMELEN];
    int err = 0;

    *t1str = '\0'; *t2str = '\0';

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	return E_DATA;
    }

    /* the varname should either be in the 2nd or 4th position */
    if (sscanf(str, "%*s %8s %8s %8s", t1str, t2str, varname) != 3) {
	if (sscanf(str, "%*s" "%8s", varname) != 1) {
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

	/* write forecast values into the newly added variable: note
	   that in this case we're not interested in computing
	   forecast standard errors
	*/
	if (pmod->ci == NLS) {
	    nls_fcast(t1, t2, (*pZ)[vi], pmod, pZ, pdinfo);
	} else if (SIMPLE_AR_MODEL(pmod->ci)) {
	    ar_fcast(t1, t2, (*pZ)[vi], pmod, (const double **) *pZ, 
		     pdinfo);
	} else if (pmod->ci == ARMA) {
	    arma_fcast(t1, t2, (*pZ)[vi], pmod, (const double **) *pZ, 
		       pdinfo);
	} else if (pmod->ci == GARCH) {
	    garch_fcast(t1, t2, (*pZ)[vi], pmod, (const double **) *pZ, 
			pdinfo);
	} else {
	    linear_fcast(t1, t2, (*pZ)[vi], pmod, (const double **) *pZ, 
			 pdinfo);
	}
    }

    return err;
}

/**
 * display_forecast:
 * @str: string giving starting and ending observations, separated
 * by a space.
 * @pmod: the model from which forecasts are wanted.
 * @pZ: pointer to data array using which @pmod was estimated.
 * @pdinfo: dataset information.
 * @opt: if includes %OPT_P, make a plot of the forecasts.
 * @prn: printing structure.
 *
 * Computes forecasts based on @pmod, over the range of observations
 * given in @str.  Forecast standard errors are also computed
 * if possible.  The results are printed to @prn, and are also
 * plotted if %OPT_P is given.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int display_forecast (const char *str, MODEL *pmod, 
		      double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    FITRESID *fr;
    int err;

    fr = get_forecast(str, pmod, pZ, pdinfo);

    if (fr == NULL) {
	return E_ALLOC;
    }

    err = fr->err;

    if (!err) {
	err = text_print_forecast(fr, pZ, pdinfo, opt, prn);
    }

    free_fit_resid(fr);
    
    return err;
}

