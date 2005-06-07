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

#define CHECK_LAGGED_DEPVAR(c) (c != NLS && c != ARMA)

static int 
allocate_fit_resid_arrays (FITRESID *fr, int errs)
{
    fr->actual = malloc(fr->nobs * sizeof *fr->actual);
    if (fr->actual == NULL) {
	return E_ALLOC;
    }

    fr->fitted = malloc(fr->nobs * sizeof *fr->fitted);
    if (fr->fitted == NULL) {
	free(fr->actual);
	fr->actual = NULL;
	return E_ALLOC;
    }

    if (errs) {
	fr->sderr = malloc(fr->nobs * sizeof *fr->sderr);
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
 * forecast standard errors, otherwise 0 (relevant only if
 * @n is greater than zero).
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

    fr->model_ID = 0;
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
	fr->nobs = n;
	if (allocate_fit_resid_arrays(fr, errs)) {
	    free(fr);
	    fr = NULL;
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

/* Make a "double list" to keep track of any "independent variables"
   that are really lags of the dependent variable.  The first element
   of the list is the number of _pairs_ of integers that follow; in
   each following pair the leading number is the ID number of a
   variable in the regression list and the second is the lag order.
   Thus for example the pair (13, 2) records the fact that the
   variable with ID number 13 is actually lag 2 of the dependent
   variable.
*/

static int process_lagged_depvar (const MODEL *pmod, 
				  const DATAINFO *pdinfo,
				  int **depvar_lags)
{
    int *dvlags = NULL;
    const char *label;
    const char *yname;
    char xname[9], tmp[9];
    int lag, nlags = 0;
    int vi, i;
    int err = 0;

    if (pmod->ci == NLS || pmod->ci == ARMA) {
	/* we won't do this for these models */
	return 0;
    }

    yname = pdinfo->varname[pmod->list[1]];

    for (i=2; i<=pmod->list[0]; i++) {
	vi = pmod->list[i];
	if (vi == LISTSEP) {
	    /* two-stage least squares */
	    break;
	}
	label = VARLABEL(pdinfo, vi);
	if ((sscanf(label, "= %8[^(](t - %d)", xname, &lag) == 2 ||
	     sscanf(label, "%8[^=]=%8[^(](-%d)", tmp, xname, &lag) == 3) &&
	    !strcmp(xname, yname)) {
	    nlags++;
	}
    }

    if (nlags > 0) {
	dvlags = malloc((nlags * 2 + 1) * sizeof *dvlags);
	if (dvlags == NULL) {
	    err = 1;
	} else {
	    dvlags[0] = nlags;
	}
    }

    if (nlags > 0 && !err) {
	int j = 1;

	for (i=2; i<=pmod->list[0]; i++) {
	    vi = pmod->list[i];
	    if (vi == LISTSEP) {
		break;
	    }
	    label = VARLABEL(pdinfo, vi);
	    if ((sscanf(label, "= %8[^(](t - %d)", xname, &lag) == 2 ||
		 sscanf(label, "%8[^=]=%8[^(](-%d)", tmp, xname, &lag) == 3) &&
		!strcmp(xname, yname)) {
		dvlags[j++] = vi;
		dvlags[j++] = lag;
	    }
	}
    } 

    *depvar_lags = dvlags;

    return err;
}

/* initialize FITRESID struct based on model and string giving
   starting and ending observations */

static int 
fit_resid_init (const char *line, const MODEL *pmod, 
		const double **Z, const DATAINFO *pdinfo,
		FITRESID *fr, int errs, gretlopt opt)
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

    if (!fr->err) {
	fr->nobs = fr->t2 - fr->t1 + 1;
	if (pmod->ci == ARMA && fr->t2 > pmod->t2) {
	    errs = 1;
	} 
	if (SIMPLE_AR_MODEL(pmod->ci) && !(opt & OPT_D)) {
	    errs = 1;
	}
	fr->err = allocate_fit_resid_arrays(fr, errs);
    }

    fr->model_ID = pmod->ID;

    return fr->err;
}

/* if an "independent" variable needed in generating a forecast is
   missing, check whether it might in fact be a lagged value of the
   dependent variable, for which we have a prior forecast available 
*/

static double 
maybe_get_yhat_lag (int vi, const int *dvlags, const double *yhat, int s,
		    int *order)
{
    double yhlag = NADBL;
    int i, lag = 0;

    for (i=1; i<=dvlags[0]; i++) {
	if (dvlags[2*i-1] == vi) {
	    lag = dvlags[2*i];
	    if (s - lag >= 0) {
		yhlag = yhat[s - lag];
	    }
	    break;
	}
    }

    if (order != NULL) {
	*order = lag;
    }

    return yhlag;
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
	}   

	/* populate Xs vector for observation */
	for (i=0; i<k && !missing; i++) {
	    double xval;

	    vi = pmod->list[i + 2];
	    xval = Z[vi][t];
	    if (na(xval)) {
		fr->sderr[s] = fr->fitted[s] = NADBL;
		missing = 1;
	    } else {
		gretl_vector_set(Xs, i, xval);
	    }
	}

	if (missing) {
	    fr->sderr[s] = fr->fitted[s] = NADBL;
	    continue;
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

/* Generate forecasts from nonlinear least squares model, using
   the string specification of the regression function that
   was saved as data on the model (see nls.c).  For now we don't
   attempt to calculate forecast error variance.
*/

static int nls_fcast (int t1, int t2, double *fcast, int offset,
		      const MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
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
	    fcast[t - offset] = (*pZ)[oldv][t];
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

/* forecasts for GARCH models -- seems as if we ought to be able to
   do something interesting with forecast error variance, but right
   now we don't do anything */

static int 
garch_fcast (int t1, int t2, double *yhat, int offset, const MODEL *pmod, 
	     const int *dvlags, const double **Z, const DATAINFO *pdinfo)
{
    double xval;
    int xvars, yno;
    int *xlist = NULL;
    int i, v, s, t;

    xlist = gretl_arma_model_get_x_list(pmod);
    yno = gretl_model_get_depvar(pmod);

    if (xlist != NULL) {
	xvars = xlist[0];
    } else {
	xvars = 0;
    }

    for (t=t1; t<=t2; t++) {
	int miss = 0;
	double yh = 0.0;

	s = t - offset;

	for (i=1; i<=xvars; i++) {
	    v = xlist[i];
	    xval = Z[v][t];
	    if (na(xval) && dvlags != NULL) {
		xval = maybe_get_yhat_lag(v, dvlags, yhat, s, NULL);
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		yh += pmod->coeff[i-1] * xval;
	    }
	}

	if (miss) {
	    yhat[s] = NADBL;
	} else {
	    yhat[s] = yh;
	}
    }

    free(xlist);

    return 0;
}

/* Compute ARMA forecast error variance (ignoring parameter
   uncertainty, as is common), via recursion. 
*/

static double 
arma_variance_machine (const double *phi, int p, 
		       const double *theta, int q,
		       double *psi, int s, double *ss_psi)
{
    int i, h = s - 1;

    /* s = number of steps ahead of forecast baseline
       h = index into "psi" array of infinite-order MA
           coefficients
    */

    if (h > p) {
	/* we don't retain more that (AR order + 1) 
	   psi values (the last one being used as
           workspace) */
	h = p;
    }

    if (s == 1) {
	/* one step ahead: initialize to unity */
	psi[h] = 1.0;
    } else {
	/* init to zero prior to adding components */
	psi[h] = 0.0;
    }

    /* add AR-derived psi[h] components */
    for (i=1; i<=p && i<s; i++) {
	psi[h] += phi[i-1] * psi[h-i];
    }

    /* add MA-derived psi[h] components */
    if (s > 1 && s <= q+1) {
	psi[h] += theta[s-2];
    }

    /* increment running sum of psi squared terms */
    *ss_psi += psi[h] * psi[h];
    
    if (s > p) {
	/* drop the oldest psi and make space for a new one */
	for (i=0; i<p; i++) {
	    psi[i] = psi[i+1];
	}
    }

    return *ss_psi;
}

/* generate forecasts for ARMA (or ARMAX) models, including
   forecast standard errors if we're doing out-of-sample
   forecasting
*/

static int 
arma_fcast (int t1, int t2, double *yhat, double *sderr,
	    int offset, const MODEL *pmod, 
	    const double **Z, const DATAINFO *pdinfo)
{
    const double *phi;
    const double *theta;
    const double *beta;

    double *psi = NULL;
    double ss_psi = 0.0;
    double xval, yval;
    int xvars, yno;
    int *xlist = NULL;
    int p, q;
    int i, s, t, tt;

    /* use pre-calculated fitted values over model estimation range,
       and don't calculate forecast error variance */
    for (t=t1; t<=pmod->t2; t++) {
	tt = t - offset;
	yhat[tt] = pmod->yhat[t];
	if (sderr != NULL) {
	    sderr[tt] = NADBL;
	}
    }

    if (t2 <= pmod->t2) {
	/* no "real" forecasts were called for, we're done */
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

    phi = pmod->coeff + pmod->ifc; /* AR coeffs */
    theta = phi + p;               /* MA coeffs */
    beta = theta + q;              /* coeffs on indep vars */

    /* setup for forecast error variance */
    if (sderr != NULL) {
	psi = malloc((p + 1) * sizeof *psi);
    }

    /* do real forecast */
    for (t=pmod->t2 + 1; t<=t2; t++) {
	int miss = 0;
	double yh = 0.0;

	tt = t - offset;

	DPRINTF(("\n *** Doing forecast for obs %d\n", t));

	/* contribution of independent variables */
	for (i=1; i<=xvars; i++) {
	    int j = 0;

	    if (xlist[i] == 0) {
		yh += pmod->coeff[0];
	    } else {
		xval = Z[xlist[i]][t];
		if (na(xval)) {
		    miss = 1;
		} else {
		    yh += beta[j++] * xval;
		}
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
		yval = yhat[s - offset];
		DPRINTF(("  AR: lag %d, yhat[%d] = %g\n", i+1, s-offset, yval));
	    }
	    if (na(yval)) {
		miss = 1;
	    } else {
		DPRINTF(("  AR: lag %d, s=%d, using coeff %g\n", i+1, s, phi[i]));
		yh += phi[i] * yval;
	    }
	}

	DPRINTF((" with AR contribution: %g\n", yh));

	/* MA contribution */
	for (i=0; i<q && !miss; i++) {
	    s = t - i - 1;
	    if (s >= pmod->t1 && s <= pmod->t2) {
		DPRINTF(("  MA: lag %d, e[%d] = %g, coeff %g\n", i+1, s, 
			 pmod->uhat[s], theta[i]));
		yh += theta[i] * pmod->uhat[s];
	    }
	}

	DPRINTF((" with MA contribution: %g\n", yh));

	if (miss) {
	    yhat[tt] = NADBL;
	} else {
	    yhat[tt] = yh;
	}

	/* forecast error variance */
	if (psi != NULL) {
	    arma_variance_machine(phi, p, theta, q,
				  psi, t - pmod->t2, 
				  &ss_psi);
	    sderr[tt] = pmod->sigma * sqrt(ss_psi);
	}
    }

    free(xlist);

    if (psi != NULL) {
	free(psi);
    }

    return 0;
}

/* construct the "phi" array of AR coefficients, based on the
   ARINFO that was added to the model at estimation time.
   The latter's "rho" member may be a compacted array, with
   zero elements omitted, but here we need a full-length
   array with zeros inserted as required */

static double *make_phi_from_arinfo (const ARINFO *arinfo, int pmax)
{
    const int *arlist = arinfo->arlist;
    int maxlag = arlist[arlist[0]];
    double *phi;
    int i;

    if (pmax > maxlag) {
	maxlag = pmax;
    }

    phi = malloc(maxlag * sizeof *phi);

    if (phi != NULL) {
	int lag;

	for (i=0; i<maxlag; i++) {
	    phi[i] = 0.0;
	}

	for (i=1; i<=arlist[0]; i++) {
	    lag = arlist[i];
	    phi[lag-1] = arinfo->rho[i-1];
	}
    }

    return phi;
}

/* determine the highest lag order of lagged dependent
   variable */

static int get_max_dv_lag (const int *dvlags)
{
    int i, pmax = 0;

    for (i=1; i<=dvlags[0]; i++) {
	if (dvlags[2*i] > pmax) {
	    pmax = dvlags[2*i];
	}
    }

    return pmax;
}

/* Set things up for computing forecast error variance for ar models
   (AR, CORC, HILU, PWE).  This is complicated by the fact that there
   may be a lagged dependent variable in the picture.  If there is,
   the effective AR coefficients have to be incremented, for the
   purposes of calculating forecast variance.  But I'm not sure this
   is quite right yet.
*/

static void set_up_ar_fcast_variance (const MODEL *pmod, int pmax,
				      double **phi, double **psi,
				      double **errphi)
{
    *errphi = make_phi_from_arinfo(pmod->arinfo, pmax);

    if (*errphi != NULL) {
	*psi = malloc((pmax + 1) * sizeof **psi);
	if (*psi == NULL) {
	    free(*errphi);
	    *errphi = NULL;
	} else {
	    *phi = malloc(pmax * sizeof **phi);
	    if (*phi == NULL) {
		free(*errphi);
		*errphi = NULL;
		free(*psi);
		*psi = NULL;
	    }
	}
    } 
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

   We also attempt to calculate forecast error variance for
   out-of-sample forecasts.  These calculations, like those for
   ARMA, do not take into account parameter uncertainty.

   This code is used for AR, CORC, HILU and PWE models.
*/

static int 
ar_fcast (int t1, int t2, double *yhat, double *sderr, 
	  int offset, const MODEL *pmod, const int *dvlags,
	  const double **Z, const DATAINFO *pdinfo,
	  gretlopt opt)
{
    const int *arlist;
    double *phi = NULL;
    double *psi = NULL;
    double *errphi = NULL;
    double ss_psi = 0.0;
    double xval, yh;
    double rk, ylag, xlag;
    int miss, p, pmax, yno;
    int i, k, v, s, t;

    yno = pmod->list[1];
    arlist = pmod->arinfo->arlist;
    p = arlist[arlist[0]]; /* AR order of error term */

    if (dvlags != NULL) {
	/* highest order of lagged dependent var */
	pmax = get_max_dv_lag(dvlags);
    } else {
	pmax = p;
    }

    if (t2 > pmod->t2 && sderr != NULL) {
	set_up_ar_fcast_variance(pmod, pmax, &phi, &psi, &errphi);
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	yh = 0.0;
	s = t - offset;

	if (t < p) {
	    yhat[s] = NADBL;
	    if (sderr != NULL) {
		sderr[s] = NADBL;
	    }
	    continue;
	}

        if (pmod->ci == PWE && t == pmod->t1) {
            /* PWE first obs is special */
            yhat[s] = pmod->yhat[t];
            continue;
        }

	if (phi != NULL) {
	    /* initialize the phi's based on the AR error process
	       alone */
	    for (i=0; i<pmax; i++) {
		phi[i] = errphi[i];
	    }
	}

	/* LHS adjustment */
	for (k=1; k<=arlist[0]; k++) {
	    rk = pmod->arinfo->rho[k-1];
	    /* FIXME static versus dynamic forecasts */
	    ylag = Z[yno][t - arlist[k]];
	    if (na(ylag) && s - arlist[k] >= 0) {
		/* use prior forecast of lagged y */
		ylag = yhat[s - arlist[k]];
	    }
	    if (na(ylag)) {
		miss = 1;
	    } else {
		yh += rk * ylag;
	    }
	}

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    v = pmod->list[i+2];
	    xval = Z[v][t];
	    if (na(xval) && dvlags != NULL) {
		int order;

		/* lagged dependent variable? */
		xval = maybe_get_yhat_lag(v, dvlags, yhat, s, &order);
		if (!na(xval)) {
		    /* augment the effective AR coefficient for the
		       purpose of computing variance */
		    phi[order - 1] += pmod->coeff[i];
		}
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		/* use rho-differenced X on RHS */
		for (k=1; k<=arlist[0]; k++) {
		    rk = pmod->arinfo->rho[k-1];
		    xlag = Z[v][t - arlist[k]];
		    if (na(xlag) && dvlags != NULL) {
			/* lagged dependent variable? */
			xlag = maybe_get_yhat_lag(v, dvlags, yhat, s - arlist[k], 
						  NULL);
		    }
		    if (!na(xlag)) {
			xval -= rk * xlag;
		    }
		}
		yh += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    yhat[s] = NADBL;
	} else {
	    yhat[s] = yh;
	}

	/* forecast error variance */
	if (phi != NULL) {
	    if (t > pmod->t2) {
		arma_variance_machine(phi, pmax, NULL, 0,
				      psi, t - pmod->t2, 
				      &ss_psi);
		sderr[s] = pmod->sigma * sqrt(ss_psi);
	    } else {
		sderr[s] = NADBL;
	    }
	}
    }

    if (psi != NULL) {
	free(phi);
	free(psi);
	free(errphi);
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

/* compute forecasts for linear models without autoregressive errors,
   version without computation of forecast standard errors
*/

static int 
linear_fcast (int t1, int t2, double *yhat, int offset, const MODEL *pmod, 
	      const int *dvlags, const double **Z, const DATAINFO *pdinfo)
{
    const double *offvar = NULL;
    double lmax = NADBL;
    double xval;
    int i, vi, t, s;

    if (pmod->ci == POISSON) {
	/* special for poisson "offset" variable */
	int offnum = gretl_model_get_int(pmod, "offset_var");

	if (offnum > 0) {
	    offvar = Z[offnum];
	}
    } else if (pmod->ci == LOGISTIC) {
	lmax = gretl_model_get_double(pmod, "lmax");
    }

    for (t=t1; t<=t2; t++) {
	int miss = 0;
	double yh = 0.0;

	s = t - offset;

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    vi = pmod->list[i+2];
	    xval = Z[vi][t];
	    if (na(xval) && dvlags != NULL) {
		/* lagged dependent variable? */
		xval = maybe_get_yhat_lag(vi, dvlags, yhat, s, NULL);
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		yh += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    yhat[s] = NADBL;
	} else if (FCAST_SPECIAL(pmod->ci)) {
	    /* special handling for LOGIT and others */
	    yhat[s] = fcast_transform(yh, pmod->ci, t, offvar, lmax);
	} else {
	    yhat[s] = yh;
	}
    }

    return 0;
}

/* driver for various functions that compute forecasts
   for different sorts of models */

static int get_fcast (FITRESID *fr, MODEL *pmod, const int *dvlags,
		      double ***pZ, DATAINFO *pdinfo,
		      gretlopt opt) 
{
    int yno = gretl_model_get_depvar(pmod);
    int s, t;
    int err = 0;

    if (pmod->ci == NLS) {
	err = nls_fcast(fr->t1, fr->t2, fr->fitted, fr->t1, pmod, pZ, pdinfo);
    } else if (SIMPLE_AR_MODEL(pmod->ci)) {
	err = ar_fcast(fr->t1, fr->t2, fr->fitted, fr->sderr, fr->t1, pmod, 
		       dvlags, (const double **) *pZ, pdinfo, opt);
    } else if (pmod->ci == ARMA) {
	err = arma_fcast(fr->t1, fr->t2, fr->fitted, fr->sderr, fr->t1, pmod, 
			 (const double **) *pZ, pdinfo);
    } else if (pmod->ci == GARCH) {
	err = garch_fcast(fr->t1, fr->t2, fr->fitted, fr->t1, pmod, 
			  dvlags, (const double **) *pZ, pdinfo);
    } else {
	err = linear_fcast(fr->t1, fr->t2, fr->fitted, fr->t1, pmod, 
			   dvlags, (const double **) *pZ, pdinfo);
    }

    for (t=fr->t1; t<=fr->t2; t++) {
	s = t - fr->t1;
	fr->actual[s] = (*pZ)[yno][t];
    }

    if (pmod->ci == ARMA) {
	/* asymptotic normal */
	fr->tval = 1.96;
    } else {
	fr->tval = tcrit95(pmod->dfd);
    }

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
 * @opt: if OPT_D, force a dynamic forecast; if OPT_S, force
 * a static forecast.
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
			double ***pZ, DATAINFO *pdinfo,
			gretlopt opt) 
{
    FITRESID *fr;
    int *dvlags = NULL;
    int full_errs = 1;

    if (AR_MODEL(pmod->ci) || FCAST_SPECIAL(pmod->ci)) {
	full_errs = 0;
    }

    fr = fit_resid_new(0, 0); 
    if (fr == NULL) {
	return NULL;
    }

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    if (dataset_is_time_series(pdinfo) && !(opt & OPT_D)) {
	/* "dynamic" (one-step ahead) forecasts need the actual value
	   of any lagged dependent variable */
	process_lagged_depvar(pmod, pdinfo, &dvlags);
    }

    fit_resid_init(str, pmod, (const double **) *pZ, pdinfo, fr, 
		   full_errs, opt);
    if (fr->err) {
	return fr;
    }   

    if (full_errs) {
	fr->err = get_static_fcast_with_errs(fr, pmod, (const double **) *pZ, 
					     pdinfo);
    } else {
	fr->err = get_fcast(fr, pmod, dvlags, pZ, pdinfo, opt);
    }

    if (dvlags != NULL) {
	free(dvlags);
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
 * @opt: if OPT_D, force a dynamic forecast; if OPT_S, force a
 * static forecast.
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
		  DATAINFO *pdinfo, gretlopt opt)
{
    int t, t1, t2, vi;
    char t1str[OBSLEN], t2str[OBSLEN], varname[VNAMELEN];
    int *dvlags = NULL;
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

    if (dataset_is_time_series(pdinfo) && !(opt & OPT_D)) {
	process_lagged_depvar(pmod, pdinfo, &dvlags);
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
	    nls_fcast(t1, t2, (*pZ)[vi], 0, pmod, pZ, pdinfo);
	} else if (SIMPLE_AR_MODEL(pmod->ci)) {
	    ar_fcast(t1, t2, (*pZ)[vi], NULL, 0, pmod, dvlags,
		     (const double **) *pZ, pdinfo, opt);
	} else if (pmod->ci == ARMA) {
	    arma_fcast(t1, t2, (*pZ)[vi], NULL, 0, pmod, (const double **) *pZ, 
		       pdinfo);
	} else if (pmod->ci == GARCH) {
	    garch_fcast(t1, t2, (*pZ)[vi], 0, pmod, dvlags,
			(const double **) *pZ, pdinfo);
	} else {
	    linear_fcast(t1, t2, (*pZ)[vi], 0, pmod, dvlags,
			 (const double **) *pZ, pdinfo);
	}
    }

    if (dvlags != NULL) {
	free(dvlags);
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
 * @opt: if includes %OPT_P, make a plot of the forecasts; if
 * includes OPT_S, force a static forecast; if includes OPT_D,
 * force a dynamic forecast.
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

    fr = get_forecast(str, pmod, pZ, pdinfo, opt);

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

