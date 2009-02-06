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
 * 
 */

#include "libgretl.h"
#include "matrix_extra.h"
#include "gretl_scalar.h"
#include "forecast.h"
#include "var.h"
#include "system.h"
#include "libset.h"

#define ARF_DEBUG 0

#ifdef max
# undef max
#endif
#define max(x,y) (((x) > (y))? (x) : (y))

/* estimators where a simple X*b does _not_ give the
   predicted value of the dependent variable */

#define FCAST_SPECIAL(c) (c == LOGIT || \
                          c == LOGISTIC || \
                          c == NLS || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == TOBIT)

#define CHECK_LAGGED_DEPVAR(c) (c != NLS && c != ARMA)

typedef struct Forecast_ Forecast;

struct Forecast_ {
    int method;       /* static, dynamic or auto */
    double *yhat;     /* array of forecast values */
    double *sderr;    /* array of forecast standard errors */
    double *eps;      /* array of estimated forecast errors */
    int *dvlags;      /* info on lagged dependent variable */
    int t1;           /* start of forecast range */
    int t2;           /* end of forecast range */
    int model_t2;     /* end of period over which model was estimated */
};

/* create an empty, dummy AR info structure for use with models
   that don't have an explicit AR error process, but that do
   have a lagged dependent variable that in effect produces
   an AR error, for forecasting purposes */

static int dummy_ar_info_init (MODEL *pmod)
{
    pmod->arinfo = malloc(sizeof *pmod->arinfo);
    if (pmod->arinfo == NULL) {
	return 1;
    }

    pmod->arinfo->arlist = gretl_null_list();
    if (pmod->arinfo->arlist == NULL) {
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    pmod->arinfo->rho = NULL;
    pmod->arinfo->sderr = NULL;

    return 0;
}

static int fit_resid_add_sderr (FITRESID *fr)
{
    int t, err = 0;

    fr->sderr = malloc(fr->nobs * sizeof *fr->sderr);

    if (fr->sderr == NULL) {
	err = E_ALLOC;
    } else {
	for (t=0; t<fr->nobs; t++) {
	    fr->sderr[t] = NADBL;
	}
    }

    return err;
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
    if (fr != NULL) {
	free(fr->actual);
	free(fr->fitted);
	free(fr->resid);
	free(fr->sderr);
	free(fr);
    }
}

/**
 * fit_resid_new_with_length:
 * @n: the number of observations to allow for, or 0 if this
 * information will be added later.
 *
 * Allocates a #FITRESID struct for holding fitted values and
 * residuals from a model (or out-of-sample forecasts).  If
 * @n is greater than 0 the arrays required for that number
 * of observations will be allocated.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

static FITRESID *fit_resid_new_with_length (int n, int add_errs)
{
    FITRESID *f = malloc(sizeof *f);

    if (f == NULL) {
	return NULL;
    }

    f->method = 0;
    f->model_ID = 0;
    f->asymp = 0;
    f->model_t1 = 0;
    f->t0 = 0;
    f->t1 = 0;
    f->t2 = 0;
    f->df = 0;
    f->k = 0;
    f->nobs = 0;
    f->pmax = PMAX_NOT_AVAILABLE;

    f->sigma = NADBL;
    f->alpha = 0.05;

    f->actual = NULL;
    f->fitted = NULL;
    f->resid = NULL;
    f->sderr = NULL;

    *f->depvar = '\0';

    f->actual = malloc(n * sizeof *f->actual);
    f->fitted = malloc(n * sizeof *f->fitted);
    f->resid = malloc(n * sizeof *f->resid);
    if (add_errs) {
	f->sderr = malloc(n * sizeof *f->sderr);
    }

    if (f->actual == NULL || f->fitted == NULL || f->resid == NULL ||
	(add_errs && f->sderr == NULL)) {
	free_fit_resid(f);
	f = NULL;
    } else {
	int t;

	for (t=0; t<n; t++) {
	    f->actual[t] = f->fitted[t] = f->resid[t] = NADBL;
	    if (f->sderr != NULL) {
		f->sderr[t] = NADBL;
	    }
	}
	f->nobs = n;
    }

    return f;
}

static FITRESID *fit_resid_new_for_model (const MODEL *pmod, 
					  const DATAINFO *pdinfo,
					  int t1, int t2, int pre_n,
					  int *err)
{
    FITRESID *fr;

    if (t1 < 0 || t2 < 0 || t2 < t1) {
	*err = E_OBS;
	return NULL;
    }

    fr = fit_resid_new_with_length(pdinfo->n, 0);

    if (fr == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    fr->t1 = t1;
    fr->t2 = t2;

    if (pre_n > 0) {
	fr->t0 = fr->t1 - pre_n;
    } else {
	fr->t0 = t1;
    }

    fr->model_ID = pmod->ID;
    fr->model_t1 = pmod->t1;

    fr->asymp = ASYMPTOTIC_MODEL(pmod->ci);

    return fr;
}

static FITRESID *fit_resid_new_for_system (int asy, 
					   const DATAINFO *pdinfo,
					   int t1, int t2, int pre_n,
					   int *err)
{
    FITRESID *fr;

    if (t1 < 0 || t2 < 0 || t2 < t1) {
	*err = E_OBS;
	return NULL;
    }

    fr = fit_resid_new_with_length(pdinfo->n, 1);

    if (fr == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    fr->t1 = t1;
    fr->t2 = t2;

    if (pre_n > 0) {
	fr->t0 = fr->t1 - pre_n;
    } else {
	fr->t0 = t1;
    }

    fr->asymp = asy;

    return fr;
}

static void fit_resid_set_dec_places (FITRESID *fr)
{
    int n = fr->t2 - fr->t0 + 1;

    if (gretl_isdummy(fr->t0, fr->t2, fr->actual) > 0) {
	fr->pmax = get_precision(fr->fitted + fr->t0, n, 6);
    } else {
	fr->pmax = get_precision(fr->actual + fr->t0, n, 6);
    }

    if (fr->pmax == 0) {
	fr->pmax = 2;
    }
}

/**
 * get_fit_resid:
 * @pmod: the model for which actual and fitted values
 * are wanted.
 * @Z: data array using which @pmod was estimated.
 * @pdinfo: dataset information.
 * @err: location to receive error code.
 *
 * Allocates a #FITRESID structure and fills it out with
 * the actual and predicted values of the dependent variable
 * in @pmod.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FITRESID *get_fit_resid (const MODEL *pmod, const double **Z, 
			 const DATAINFO *pdinfo, int *err)
{
    FITRESID *fr;
    int t, dv = -1;

    if (pmod->ci != INTREG) {
	dv = gretl_model_get_depvar(pmod);
	if (dv < 0 || dv >= pdinfo->v) {
	    *err = E_DATA;
	    return NULL;
	}
    }

    fr = fit_resid_new_for_model(pmod, pdinfo, pmod->t1, pmod->t2,
				 0, err);
    if (*err) {
	return NULL;
    }

    if (LIMDEP(pmod->ci)) {
	fr->sigma = NADBL;
    } else {
	fr->sigma = gretl_model_get_double(pmod, "sigma_orig");
	if (na(fr->sigma)) {
	    fr->sigma = pmod->sigma;
	}
    }

    for (t=0; t<fr->nobs; t++) {
	if (dv < 0) {
	    if (na(pmod->yhat[t]) || na(pmod->uhat[t])) {
		fr->actual[t] = NADBL;
	    } else {
		fr->actual[t] = pmod->yhat[t] + pmod->uhat[t];
	    }
	} else {
	    fr->actual[t] = Z[dv][t];
	}
	fr->fitted[t] = pmod->yhat[t];
	fr->resid[t] = pmod->uhat[t];
    }

    fit_resid_set_dec_places(fr);
    
    if (dv < 0) {
	strcpy(fr->depvar, "implicit y");
    } else {
	strcpy(fr->depvar, pdinfo->varname[dv]);
    }
    
    return fr;
}

/* local shortcut to get a model's list of regressors */

static const int *model_xlist (MODEL *pmod)
{
    int *xlist = (int *) gretl_model_get_data(pmod, "xlist");

    if (xlist == NULL) {
	xlist = gretl_model_get_x_list(pmod);
	if (xlist != NULL) {
	    gretl_model_set_list_as_data(pmod, "xlist", xlist);
	} 
    }

    return xlist;
}

/* try to figure out if a model has any lags of the dependent
   variable among the regressors: return 1 if so, 0 if not
*/

static int 
has_depvar_lags (MODEL *pmod, const DATAINFO *pdinfo)
{
    const char *yname;
    const int *xlist;
    int i, vi;

    xlist = model_xlist(pmod);
    if (xlist == NULL) {
	return 0;
    }

    yname = gretl_model_get_depvar_name(pmod, pdinfo);

    for (i=1; i<=xlist[0]; i++) {
        vi = xlist[i];
	if (pdinfo->varinfo[vi]->lag > 0 &&
	    !strcmp(yname, pdinfo->varinfo[vi]->parent)) {
	    return 1;
	}
    }

    return 0;
}

/* Makes a list to keep track of any "independent variables" that are
   really lags of the dependent variable.  The list has as many
   elements as the model has independent variables, and in each place
   we either write a zero (if the coefficient does not correspond to a
   lag of the dependent variable) or a positive integer corresponding
   to the lag order.  However, In case the list of independent vars
   contains no lagged dependent var, *depvar_lags is set to NULL.
   Returns 1 on error, 0 otherwise.
*/

static int process_lagged_depvar (MODEL *pmod, 
				  const DATAINFO *pdinfo,
				  int **depvar_lags)
{
    const int *xlist = NULL;
    int *dvlags = NULL;
    int anylags;
    int err = 0;

    anylags = has_depvar_lags(pmod, pdinfo);

    if (!anylags) {
	*depvar_lags = NULL;
	return 0;
    }

    xlist = model_xlist(pmod);

    dvlags = malloc(xlist[0] * sizeof *dvlags);

    if (dvlags == NULL) {
	err = E_ALLOC;
    } else {
	const char *yname;
	int i, vi;

	yname = gretl_model_get_depvar_name(pmod, pdinfo);

	for (i=1; i<=xlist[0]; i++) {
            vi = xlist[i];
	    if (pdinfo->varinfo[vi]->lag > 0 &&
		!strcmp(yname, pdinfo->varinfo[vi]->parent)) {
		dvlags[i-1] = pdinfo->varinfo[vi]->lag;
	    } else {
		dvlags[i-1] = 0;
	    }
	}
    } 

    *depvar_lags = dvlags;

    return err;
}

/* Tries to determine if a model has any "real" exogenous regressors:
   we discount a simple time trend and periodic dummy variables, since
   these can be extended automatically.  If dvlags is non-NULL we use
   it to screen out "independent vars" that are really lags of the
   dependent variable.
*/

static int 
has_real_exog_regressors (MODEL *pmod, const int *dvlags,
			  const double **Z, const DATAINFO *pdinfo)
{
    const int *xlist = model_xlist(pmod);
    int i, xi, ret = 0;

    if (xlist != NULL) {
	for (i=0; i<xlist[0]; i++) {
	    xi = xlist[i + 1];
	    if (xi != 0 && (dvlags == NULL || dvlags[i] == 0)) {
		if (is_trend_variable(Z[xi], pdinfo->n)) {
		    continue;
		} else if (is_periodic_dummy(Z[xi], pdinfo)) {
		    continue;
		} else {
		    ret = 1;
		    break;
		}
	    }
	}
    } 

    return ret;
}

#define AR_DEBUG 0

/* Get a value for a lag of the dependent variable.  If method is
   dynamic we prefer lagged prediction to lagged actual.  If method is
   static, we never want the lagged prediction, only the actual.  If
   method is "auto", which value we prefer depends on whether we're in
   or out of sample (actual within, lagged prediction without).
*/

static double fcast_get_ldv (Forecast *fc, int i, int t, int lag,
			     const double **Z)
{
    double ldv;

    /* initialize to actual lagged value, if available */
    if (t - lag < 0) {
	ldv = NADBL;
    } else {
	ldv = Z[i][t-lag];
    }

#if AR_DEBUG
    fprintf(stderr, "fcast_get_ldv: i=%d, t=%d, lag=%d; "
	    "initial ldv = Z[%d][%d] = %g\n", i, t, lag, i, t-lag, ldv);
#endif

    if (fc->method != FC_STATIC) {
	int yht = t - lag;

#if AR_DEBUG
	fprintf(stderr, "fcast_get_ldv (non-static): yht = %d\n", yht);
#endif
	if (fc->method == FC_DYNAMIC && yht >= 0) {
	    if (!na(fc->yhat[yht])) {
		ldv = fc->yhat[yht];
	    }
	} else if (fc->method == FC_AUTO && yht >= 0) {
	    if (t > fc->model_t2 + lag || na(ldv)) {
		ldv = fc->yhat[yht];
#if AR_DEBUG
		fprintf(stderr, "fcast_get_ldv: reset ldv = yhat[%d] = %g\n",
			yht, ldv);
#endif
	    } 
	}
    }

    return ldv;
}

/* Get forecasts, plus standard errors for same, for models without
   autoregressive errors and without "special requirements"
   (e.g. nonlinearity).  The forecast standard errors include both
   uncertainty over the error process and parameter uncertainty
   (Davidson and MacKinnon method).
*/

static int
static_fcast_with_errs (Forecast *fc, MODEL *pmod, 
			const double **Z, const DATAINFO *pdinfo) 
{
    gretl_matrix *V = NULL;
    gretl_vector *Xs = NULL;
    gretl_vector *b = NULL;
    double s2 = pmod->sigma * pmod->sigma;
    double vyh, xval;
    int k = pmod->ncoeff;
    int i, vi, t;
    int err = 0;

    V = gretl_vcv_matrix_from_model(pmod, NULL, &err);
    if (err) {
	goto bailout;
    }

    b = gretl_coeff_vector_from_model(pmod, NULL, &err);
    if (err) {
	goto bailout;
    }

    Xs = gretl_vector_alloc(k);
    if (Xs == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (t=fc->t1; t<=fc->t2 && !err; t++) {
	int missing = 0;

	/* skip if we can't compute forecast */
	if (t >= pmod->t1 && t <= pmod->t2) {
	    missing = na(pmod->yhat[t]);
	}   

	/* populate Xs vector for observation */
	for (i=0; i<k && !missing; i++) {
	    vi = pmod->list[i + 2];
	    xval = Z[vi][t];
	    if (na(xval)) {
		fc->sderr[t] = fc->yhat[t] = NADBL;
		missing = 1;
	    } else {
		gretl_vector_set(Xs, i, xval);
	    }
	}

	if (missing) {
	    fc->sderr[t] = fc->yhat[t] = NADBL;
	    continue;
	}

	/* forecast value */
	fc->yhat[t] = gretl_vector_dot_product(Xs, b, NULL);

	/* forecast variance */
	vyh = gretl_scalar_qform(Xs, V, &err);
	if (na(vyh)) {
	    err = 1;
	} else {
	    vyh += s2;
	    if (vyh >= 0.0) {
		fc->sderr[t] = sqrt(vyh);
	    } else {
		fc->sderr[t] = NADBL;
		err = 1;
	    }
	}
    }

 bailout:

    gretl_matrix_free(V);
    gretl_vector_free(Xs);
    gretl_vector_free(b);

    return err;
}

/* Generate forecasts from nonlinear least squares model, using the
   string specification of the regression function that was saved as
   data on the model (see nls.c).  If the NLS formula is dynamic and
   the user has not requested a static forecast, we do an
   autoregressive genr out of sample.  If we're doing a static
   forecast we add a simple-minded forecast error, namely the standard
   error of the NLS regression.
*/

static int nls_fcast (Forecast *fc, const MODEL *pmod, 
		      double ***pZ, DATAINFO *pdinfo)
{
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    int fcv = pdinfo->v;
    int yno = 0;
    const char *nlfunc;   
    double *y = NULL;
    char formula[MAXLINE];
    int t, err = 0;

    nlfunc = gretl_model_get_data(pmod, "nl_regfunc");
    if (nlfunc == NULL) {
	err = E_DATA;
    }

    if (fc->method == FC_AUTO) {
	yno = pmod->list[1];
	y = copyvec((*pZ)[yno], pdinfo->n);
	if (y == NULL) {
	    err = E_ALLOC;
	} 
    }

    if (!err) {
	/* upper limit of static forecast */
	int t2 = (fc->method == FC_STATIC)? fc->t2 : pmod->t2;

	if (t2 >= fc->t1) {
	    /* non-null static range */
	    pdinfo->t1 = fc->t1;
	    pdinfo->t2 = t2;
	    sprintf(formula, "$nl_y = %s", nlfunc);
	    err = generate(formula, pZ, pdinfo, OPT_P, NULL);
	    if (!err) {
		for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		    fc->yhat[t] = (*pZ)[fcv][t];
		}
	    }
	}

	if (!err && fc->method == FC_AUTO && fc->t2 > pmod->t2) {
	    /* dynamic forecast out of sample */
	    pdinfo->t1 = pmod->t2 + 1;
	    pdinfo->t2 = fc->t2;
	    strcpy(formula, pmod->depvar);
	    err = generate(formula, pZ, pdinfo, OPT_P, NULL);
	    if (!err) {
		for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		    fc->yhat[t] = (*pZ)[yno][t];
		}
	    }
	}
    }

    if (pdinfo->v > fcv) {
	err = dataset_drop_last_variables(pdinfo->v - fcv, pZ, pdinfo);
    }

    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    if (y != NULL) {
	/* restore original dependent variable */
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[yno][t] = y[t];
	}
	free(y);
    }

    if (!err && fc->method == FC_STATIC && fc->sderr != NULL) {
#if 0   /* not quite yet */
	nls_boot_calc(pmod, pZ, pdinfo, fc->t1, fc->t2, fc->sderr);
#else
	/* by request, but is it a good idea? */
	for (t=fc->t1; t<=fc->t2; t++) {
	    fc->sderr[t] = pmod->sigma;
	}
#endif
    }

    return err;
}

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

#define depvar_lag(f,i) ((f->dvlags != NULL)? f->dvlags[i] : 0)

#define GFC_DEBUG 0

/*   
   GARCH error variance process:

   h(t) = a(0) + sum(i=1 to q) a(i) * u(t-i)^2 
                   + sum(j=1 to p) b(j) * h(t-j)

   This is then complexified if the model includes lags of the
   dependent variable among the regressors, but in the following
   function we just calculate the successive h's.
*/

static double *garch_h_hat (const MODEL *pmod, int t1, int t2,
			    int xvars)
{
    const double *alpha;
    const double *beta;
    double hlag;
    double *h, *mh;
    int nf, q, p, lmax;
    int i, s, t, ti;
    int err = 0;

    mh = gretl_model_get_data(pmod, "garch_h");
    if (mh == NULL) {
	return NULL;
    }

    nf = t2 - t1 + 1;
    if (nf <= 0) {
	return NULL;
    }

    h = malloc(nf * sizeof *h);
    if (h == NULL) {
	return NULL;
    }

    q = pmod->list[1];
    p = pmod->list[2];
    alpha = pmod->coeff + xvars;
    beta = alpha + q + 1;

    lmax = max(p, q);

    for (t=t1; t<=t2 && !err; t++) {
	s = t - t1;

	h[s] = alpha[0];
#if GFC_DEBUG
	fprintf(stderr, "h[%d or %d] init'd to %g\n", s, t, h[s]);
#endif
	for (i=1; i<=lmax; i++) {
	    ti = t - i;
	    if (ti < 0) {
		break;
	    }
	    hlag = 0.0;
	    if (ti <= pmod->t2) {
		if (!na(mh[ti])) {
		    hlag = mh[ti];
		}
	    } else {
		hlag = h[s-i];
	    }
	    if (i <= q) {
		if (ti <= pmod->t2 && !na(pmod->uhat[ti])) {
		    h[s] += alpha[i] * pmod->uhat[ti] * pmod->uhat[ti];
#if GFC_DEBUG
		    fprintf(stderr, " added alpha[%d] * uhat[%d]^2 = %g * %g = %g\n",
			    i, ti, alpha[i], pmod->uhat[ti] * pmod->uhat[ti],
			    alpha[i] * pmod->uhat[ti] * pmod->uhat[ti]);
#endif
		} else {
		    /* lagged uhat^2 not available: substitute its
		       expectation, h at the same date */
		    h[s] += alpha[i] * hlag;
		}
	    }
	    if (i <= p) {
		h[s] += beta[i-1] * hlag;
#if GFC_DEBUG
		fprintf(stderr, " added beta[%d] * h[-%d] = %g * %g = %g\n",
			i-1, i, beta[i-1], hlag, beta[i-1] * hlag);
#endif
	    }
	}
	if (h[s] < 0.0) {
	    err = 1;
	}
#if GFC_DEBUG
	fprintf(stderr, "h[%d or %d] = %g (%g)\n", s, t, h[s], sqrt(h[s]));
#endif
    }

    if (err) {
	free(h);
	h = NULL;
    }

    return h;
}

/* Infinite MA representation, in case GARCH model has lags of the
   dependent var among the regressors (the autoregressive coeffs are
   recorded in the array phi).
*/

static double *garch_psi (const double *phi, int p, int nf)
{
    double *psi;
    int i, s;

    if (phi == NULL) {
	return NULL;
    }

    psi = malloc(nf * sizeof *psi);
    if (psi == NULL) {
	return NULL;
    }

    psi[0] = 1.0;

    /* Are we OK with regard to signs below?  Do we need to keep a
       record of the squares of the psi elements too? (Or instead?)
    */

    for (s=1; s<nf; s++) {
	psi[s] = 0.0;

	/* add psi[s] components derived from dep var lags */
	for (i=0; i<p && i<s; i++) {
	    psi[s] += phi[i] * psi[s-i-1];
#if GFC_DEBUG
	    fprintf(stderr, "psi[%d]: adding phi[%d] * psi[%d] = %g\n",
		    s, i, s-i-1, phi[i] * psi[s-i-1]);
#endif
	}
    }

#if GFC_DEBUG
    for (s=0; s<nf; s++) {
	fprintf(stderr, "psi[%d] = %g\n", s, psi[s]);
    }    
#endif

    return psi;
}

/* Construct an array, "phi", of autoregressive coefficients, if a
   GARCH model has any lags of the dependent variable among the
   regressors; the array has dimension equal to the highest-order lag
   of the dependent variable.
*/

static double *garch_ldv_phi (const MODEL *pmod, const int *xlist,
			      const int *dvlags, int *ppmax)
{
    double *phi = NULL;
    int i, pmax = 0;

    if (dvlags != NULL && xlist != NULL) {
	for (i=0; i<xlist[0]; i++) {
	    if (dvlags[i] > pmax) {
		pmax = dvlags[i];
	    }
	}

	phi = malloc(pmax * sizeof *phi);
	if (phi != NULL) {
	    for (i=0; i<pmax; i++) {
		phi[i] = 0.0;
	    }
	    for (i=0; i<xlist[0]; i++) {
		if (dvlags[i] > 0) {
		    phi[dvlags[i] - 1] = pmod->coeff[i];
		}
	    }
	}
    }

    if (phi != NULL) {
#if GFC_DEBUG
	fprintf(stderr, "garch_ldv_phi: pmax = %d\n", pmax);
	for (i=0; i<pmax; i++) {
	    fprintf(stderr, "phi[%d] = %g\n", i, phi[i]);
	}
#endif
	*ppmax = pmax;
    }

    return phi;
}

/* Compute forecast standard errors for a GARCH model that includes
   lags of the dependent variable among the regressors.  This could do
   with some scrutiny (or the calculation, above, of the "psi" array
   could do with checking).
*/

static double garch_ldv_sderr (const double *h,
			       const double *psi,
			       int s)
{
    double ss, vs = h[s];
    int i;

    for (i=1; i<=s; i++) {
	vs += h[s-i] * psi[i] * psi[i];
    }

    if (vs >= 0.0) {
	ss = sqrt(vs);
    } else {
	ss = NADBL;
    }

    return ss;
}

static int garch_fcast (Forecast *fc, MODEL *pmod, 
			const double **Z, const DATAINFO *pdinfo)
{
    double xval;
    int xvars, yno;
    const int *xlist = NULL;
    double *h = NULL;
    double *mh = NULL;
    double *phi = NULL;
    double *psi = NULL;
    int i, v, t;

    xlist = model_xlist(pmod);
    yno = gretl_model_get_depvar(pmod);

    if (xlist != NULL) {
	xvars = xlist[0];
    } else {
	xvars = 0;
    }

    if (fc->sderr != NULL) {
	h = garch_h_hat(pmod, pmod->t2 + 1, fc->t2, xvars);
	mh = gretl_model_get_data(pmod, "garch_h");
    }

    if (fc->dvlags != NULL) {
	int ns = fc->t2 - pmod->t2;
	int pmax = 0;
	
	if (ns > 0) {
	    phi = garch_ldv_phi(pmod, xlist, fc->dvlags, &pmax);
	    psi = garch_psi(phi, pmax, ns);
	    free(phi);
	}
    }

    for (t=fc->t1; t<=fc->t2; t++) {
	int lag, miss = 0;
	double yh = 0.0;

	if (fc->method != FC_DYNAMIC && 
	    t >= pmod->t1 && t <= pmod->t2) {
	    fc->yhat[t] = pmod->yhat[t];
	    if (fc->sderr != NULL && mh != NULL) {
		fc->sderr[t] = mh[t];
	    }
	    continue;
	}

	for (i=0; i<xvars; i++) {
	    v = xlist[i+1];
	    if ((lag = depvar_lag(fc, i))) {
		xval = fcast_get_ldv(fc, yno, t, lag, Z);
	    } else {
		xval = Z[v][t];
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		yh += pmod->coeff[i] * xval;
	    }
	}

	if (miss) {
	    fc->yhat[t] = NADBL;
	} else {
	    fc->yhat[t] = yh;
	}

	if (h != NULL) {
	    if (t > pmod->t2) {
		if (psi != NULL) {
		    /* build in effect of lagged dependent var */
		    fc->sderr[t] = garch_ldv_sderr(h, psi, t - pmod->t2 - 1);
		} else {
		    /* no lagged dependent variable */
		    fc->sderr[t] = sqrt(h[t - pmod->t2 - 1]);
		}
	    } else {
		fc->sderr[t] = NADBL;
	    }
	}
    }

    if (h != NULL) {
	free(h);
    }
    if (psi != NULL) {
	free(psi);
    }

    return 0;
}

/* Compute ARMA forecast error variance (ignoring parameter
   uncertainty, as is common), via recursion.  Cf. Box and Jenkins,
   1976, p. 508, "Program 4", V(l) algorithm (with the sign of theta
   changed).
*/

static double arma_variance (const double *phi, int p, 
			     const double *theta, int q,
			     double *psi, int npsi, int l)
{
    static double sspsi;
    int i, j;

    if (l == 1) {
	sspsi = 0.0;
	psi[0] = 1.0;
	for (j=1; j<npsi; j++) {
	    psi[j] = 0.0;
	    for (i=1; i<=j; i++) {
		if (i <= p) {
		    psi[j] += phi[i] * psi[j-i];
		}
		if (j <= q && theta != NULL) {
		    psi[j] += theta[j];
		}
	    }
	}
    } 

    sspsi += psi[l-1] * psi[l-1];

    return sspsi;
}

/* When forecasting based on an armax model estimated using X12A,
   or via the Kalman filter, we need to form the series X\beta
   so that we can subtract X\beta_{t-i} from y_{t-i} in
   computing the AR portion of the forecast. "beta" below
   is the array of ARMAX coefficients, not including the
   constant.
*/

static double *create_Xb_series (Forecast *fc, const MODEL *pmod,
				 const double *beta, const int *xlist, 
				 const double **Z)
{
    double *Xb;
    double x;
    int miss;
    int i, j, t;

    Xb = malloc((fc->t2 + 1) * sizeof *Xb);
    if (Xb == NULL) {
	return NULL;
    }

    for (t=0; t<=fc->t2; t++) {
	Xb[t] = 0.0;
	miss = 0;
	j = 0;
	for (i=1; i<=xlist[0] && !miss; i++) {
	    if (xlist[i] == 0) {
		Xb[t] += pmod->coeff[0];
	    } else {
		x = Z[xlist[i]][t];
		if (na(x)) {
		    Xb[t] = NADBL;
		    miss = 1;
		} else {
		    Xb[t] += beta[j++] * x;
		}
	    }
	}
    }

    return Xb;
}

static int want_x_beta_prep (const MODEL *pmod, const int *xlist)
{
    int ret = 0;

    if (xlist != NULL) {
	int aflags = gretl_model_get_int(pmod, "arma_flags");

	if ((aflags & ARMA_EXACT) || (aflags & ARMA_X12A)) {
	    ret = 1;
	}
    }

    return ret;
}

/* generate forecasts for AR(I)MA (or ARMAX) models, including
   forecast standard errors if we're doing out-of-sample
   forecasting
*/

static int arma_fcast (Forecast *fc, MODEL *pmod, 
		       const double **Z, const DATAINFO *pdinfo)
{
    double *psi = NULL;
    double *phi = NULL;
    double *theta = NULL;
    double *phi0 = NULL;
    double *Xb = NULL;
    const double *beta;
    const double *y;

    double xval, yval, vl;
    double mu = NADBL;
    int xvars, yno;
    const int *xlist = NULL;
    int p, q, px = 0, npsi = 0;
    int fcstart = fc->t1;
    int ar_smax, ma_smax;
    int regarma = 0;
    int i, s, t;
    int err = 0;

    DPRINTF(("\n\n*** arma_fcast: METHOD = %d\n", fc->method));

    if (fc->method != FC_DYNAMIC) {
	/* use pre-calculated fitted values over model estimation range,
	   and don't bother calculating forecast error variance */
	for (t=fc->t1; t<=pmod->t2; t++) {
	    fc->yhat[t] = pmod->yhat[t];
	    if (fc->sderr != NULL) {
		fc->sderr[t] = NADBL;
	    }
	}
	if (fc->t2 <= pmod->t2) {
	    /* no "real" forecasts were called for, we're done */
	    return 0;
	} 
	fcstart = pmod->t2 + 1;
    }

    p = arma_model_max_AR_lag(pmod);
    q = arma_model_max_MA_lag(pmod);

    xlist = model_xlist(pmod);
    yno = gretl_model_get_depvar(pmod);

    DPRINTF(("forecasting variable %d (%s), obs %d to %d, with p=%d, q=%d\n", yno, 
	     pdinfo->varname[yno], fc->t1, fc->t2, p, q));

    xvars = (xlist != NULL)? xlist[0] : 0;

    err = arma_model_integrated_AR_MA_coeffs(pmod, &phi, &theta);
    if (err) {
	goto bailout;
    }

    beta = arma_model_get_x_coeffs(pmod);

    if (want_x_beta_prep(pmod, xlist)) {
	if (gretl_is_arima_model(pmod)) {
	    regarma = 1;
	    err = regarma_model_AR_coeffs(pmod, &phi0, &px);
	}
	if (xlist[0] == 1 && xlist[1] == 0) {
	    /* just a const, no ARMAX */
	    mu = pmod->coeff[0];
	    if (phi0 != NULL) {
		mu = 1.0;
		for (i=1; i<=px; i++) {
		    mu -= phi0[i];
		}
		mu *= pmod->coeff[0];
	    } else {
		mu = pmod->coeff[0];
	    }
	} else {
	    /* we have ARMAX terms */
	    Xb = create_Xb_series(fc, pmod, beta, xlist, Z);
	    if (Xb == NULL) {
		err = E_ALLOC;
	    }
	}
	if (err) {
	    goto bailout;
	}
    }

    /* setup for forecast error variance */
    if (fc->sderr != NULL) {
	npsi = fc->t2 - fcstart + 1;
	psi = malloc(npsi * sizeof *psi);
    }

    /* cut-off points for using actual rather than forecast
       values of y in generating further forecasts (FIXME??) */
    if (fc->method == FC_STATIC) {
	ar_smax = fc->t2;
	ma_smax = pmod->t2;
    } else if (fc->method == FC_DYNAMIC) {
	ar_smax = max(fc->t1 - 1, p);
	ma_smax = max(fc->t1 - 1, q);
    } else {
	ar_smax = pmod->t2;
	ma_smax = pmod->t2;
    }

    DPRINTF(("ar_smax = %d, ma_smax = %d\n", ar_smax, ma_smax));

    /* dependent variable */
    y = Z[yno];

    /* do real forecast */
    for (t=fcstart; t<=fc->t2 && !err; t++) {
	double yh = 0.0;
	int miss = 0;

	DPRINTF(("\n *** Doing forecast for obs %d\n", t));

	/* contribution of const and/or independent variables */

	if (!na(mu)) {
	    yh = mu;
	} else if (Xb != NULL) {
	    /* X\beta series is pre-computed */
	    if (na(Xb[t])) {
		miss = 1;
	    } else { 
		yh = Xb[t];
	    }
	} else {
	    int j = 0;

	    for (i=1; i<=xvars; i++) {
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
	}

	if (phi0 != NULL && Xb != NULL) {
	    /* the regarma case */
	    for (i=1; i<=px; i++) {
		s = t - i;
		if (s < 0) {
		    miss = 1;
		} else if (na(Xb[s])) {
		    miss = 1;
		} else {
		    yh -= phi0[i] * Xb[s];
		} 
	    }
	}

	DPRINTF((" Xb contribution = %g\n", yh));

	/* AR contribution (incorporating any differencing) */

	for (i=1; i<=p && !miss; i++) {
	    if (phi[i] == 0.0) {
		continue;
	    }
	    s = t - i;
	    if (s < 0) {
		yval = NADBL;
	    } else if (s <= ar_smax) {
		yval = y[s];
		DPRINTF(("  AR: lag %d, y[%d] = %g\n", i, s, yval));
	    } else {
		yval = fc->yhat[s];
		DPRINTF(("  AR: lag %d, yhat[%d] = %g\n", i, s, yval));
	    }
	    if (na(yval)) {
		DPRINTF(("  AR: lag %d, missing value\n", i));
		miss = 1;
	    } else {
		DPRINTF(("  AR: lag %d, using coeff %#.8g\n", i, phi[i]));
		if (!regarma && Xb != NULL) {
		    if (na(Xb[s])) {
			miss = 1;
		    } else {
			yh += phi[i] * (yval - Xb[s]);
		    }
		} else if (!regarma && !na(mu)) {
		    yh += phi[i] * (yval - mu);
		} else {
		    yh += phi[i] * yval;
		}
	    }
	}

	DPRINTF((" with AR contribution: %g\n", yh));

	/* MA contribution */

	for (i=1; i<=q && !miss; i++) {
	    if (theta[i] == 0.0) {
		continue;
	    }
	    s = t - i;
	    if (s >= pmod->t1 && s <= ma_smax) {
		DPRINTF(("  MA: lag %d, e[%d] = %g, theta[%d] = %g\n", i, s, 
			 pmod->uhat[s], i, theta[i]));
		yh += theta[i] * pmod->uhat[s];
	    } else if (fc->eps != NULL) {
		DPRINTF(("  MA: lag %d, ehat[%d] = %g, theta[%d] = %g\n", i, s, 
			 fc->eps[s], i, theta[i]));
		yh += theta[i] * fc->eps[s];
	    }
	}

	DPRINTF((" with MA contribution: %g\n", yh));

	if (miss) {
	    fc->yhat[t] = NADBL;
	} else {
	    fc->yhat[t] = yh;
	    if (fc->eps != NULL) {
		/* form estimated error: we do this only in
		   the case of a "static" forecast */
		fc->eps[t] = y[t] - yh;
		DPRINTF(("\n setting ehat[%d] = %g - %g = %g\n", t, 
			 y[t], yh, fc->eps[t]));
	    }
	}

	/* forecast error variance */
	if (psi != NULL) {
	    vl = arma_variance(phi, p, theta, q, psi, npsi, t - fcstart + 1);
	    fc->sderr[t] = pmod->sigma * sqrt(vl);
	}
    }

 bailout:

    free(psi);
    free(phi);
    free(theta);
    free(Xb);
    free(phi0);

    return err;
}

/* construct the "phi" array of AR coefficients, based on the ARINFO
   that was added to the model at estimation time.  The latter's rho
   member may be a compacted array, with zero elements omitted, but
   here we need a full-length array with zeros inserted as required.
*/

static double *make_phi_from_arinfo (const ARINFO *arinfo, int pmax)
{
    double *phi = malloc((pmax + 1) * sizeof *phi);

    if (phi != NULL) {
	int i, lag;

	for (i=0; i<=pmax; i++) {
	    phi[i] = 0.0;
	}

	for (i=1; i<=arinfo->arlist[0]; i++) {
	    lag = arinfo->arlist[i];
	    phi[lag] = arinfo->rho[i-1];
	}
    }

    return phi;
}

/* Determine the greatest lag order for a model, either via explicit
   AR error process or via inclusion of lagged dependent var as
   regressor.
*/

static int max_ar_lag (Forecast *fc, const MODEL *pmod, int p)
{
    int i, pmax = p; /* explicit AR order */

    if (fc->dvlags != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    if (fc->dvlags[i] > pmax) {
		pmax = fc->dvlags[i];
	    }
	}
    }

    return pmax;
}

/* Set things up for computing forecast error variance for ar models
   (AR, AR1).  This is complicated by the fact that there
   may be a lagged dependent variable in the picture.  If there is,
   the effective AR coefficients have to be incremented, for the
   purposes of calculating forecast variance.  But I'm not sure this
   is quite right yet.
*/

static int
set_up_ar_fcast_variance (Forecast *fc, const MODEL *pmod, 
			  int pmax, int npsi, double **pphi, double **ppsi,
			  double **perrphi)
{
    double *errphi = NULL;
    double *psi = NULL;
    double *phi = NULL;
    int err = 0;
    
    errphi = make_phi_from_arinfo(pmod->arinfo, pmax);
    psi = malloc(npsi * sizeof *psi);
    phi = malloc((pmax + 1) * sizeof *phi);

    if (errphi == NULL || psi == NULL || phi == NULL) {
	free(errphi);
	free(psi);
	free(phi);
	err = E_ALLOC;
    } else {
	*perrphi = errphi;
	*pphi = phi;
	*ppsi = psi;
    }

    return err;
}

/* 
   The code below generates forecasts that incorporate the 
   predictable portion of an AR error term:

       u_t = r1 u_{t-1} + r2 u_{t-1} + ... + e_t

   where e_t is white noise.  The forecasts are based on the
   representation of a model with such an error term as

       (1 - r(L)) y_t = (1 - r(L)) X_t b + e_t

   or

       y_t = r(L) y_t + (1 - r(L)) X_t b + e_t

   where r(L) is a polynomial in the lag operator.

   We also attempt to calculate forecast error variance for
   out-of-sample forecasts.  These calculations, like those for
   ARMA, do not take into account parameter uncertainty.

   This code is used for AR and AR1 models; it is also used for
   dynamic forecasting with models that do not have an explicit AR
   error process but that have one or more lagged values of the
   dependent variable as regressors.
*/

static int ar_fcast (Forecast *fc, MODEL *pmod, 
		     const double **Z, const DATAINFO *pdinfo)
{
    const int *arlist;
    double *phi = NULL;
    double *psi = NULL;
    double *errphi = NULL;
    double xval, yh, vl;
    double rk, ylag, xlag;
    int miss, yno;
    int i, k, v, t, tk;
    int p, dvlag, pmax = 0;
    int pwe, npsi = 0;
    int err = 0;

#if AR_DEBUG
    fprintf(stderr, "\n*** ar_fcast, method = %d\n\n", fc->method);
#endif

    yno = pmod->list[1];
    arlist = pmod->arinfo->arlist;
    p = arlist[arlist[0]]; /* AR order of error term */

    if (fc->t2 > pmod->t2 && fc->sderr != NULL) {
	/* we'll compute variance only if we're forecasting out of
	   sample */
	pmax = max_ar_lag(fc, pmod, p);
	npsi = fc->t2 - fc->t1 + 1;
	set_up_ar_fcast_variance(fc, pmod, pmax, npsi, &phi, &psi, &errphi);
    }

    pwe = (pmod->opt & OPT_P);

    for (t=fc->t1; t<=fc->t2; t++) {
	miss = 0;
	yh = 0.0;

	if (t < p) {
	    fc->yhat[t] = NADBL;
	    if (fc->sderr != NULL) {
		fc->sderr[t] = NADBL;
	    }
	    continue;
	}

        if (pwe && t == pmod->t1) {
            /* PWE first obs is special */
            fc->yhat[t] = pmod->yhat[t];
            continue;
        }

	if (phi != NULL) {
	    /* initialize the phi's based on the AR error process
	       alone */
	    for (i=1; i<=pmax; i++) {
		phi[i] = errphi[i];
	    }
	}

	/* r(L) y_t */
	for (k=1; k<=arlist[0]; k++) {
	    rk = pmod->arinfo->rho[k-1];
	    tk = t - arlist[k];
	    ylag = fcast_get_ldv(fc, yno, tk, 0, Z);
	    if (na(ylag)) {
		miss = 1;
	    } else {
		yh += rk * ylag;
	    }
	}

	/* (1 - r(L)) X_t b */
	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    v = pmod->list[i+2];
	    if ((dvlag = depvar_lag(fc, i))) {
		xval = fcast_get_ldv(fc, yno, t, dvlag, Z);
	    } else {
		xval = Z[v][t];
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		if (dvlag > 0 && phi != NULL) {
		    /* augment phi for computation of variance */
		    phi[dvlag] += pmod->coeff[i];
		}
		for (k=1; k<=arlist[0]; k++) {
		    rk = pmod->arinfo->rho[k-1];
		    tk = t - arlist[k];
		    if (dvlag > 0) {
			xlag = fcast_get_ldv(fc, yno, tk, dvlag, Z);
		    } else {
			xlag = Z[v][tk];
		    }
		    if (!na(xlag)) {
			xval -= rk * xlag;
		    }
		}
		yh += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    fc->yhat[t] = NADBL;
	} else {
	    fc->yhat[t] = yh;
	}

	/* forecast error variance */
	if (phi != NULL && pmod->ci != GARCH) {
	    if (t > pmod->t2) {
		vl = arma_variance(phi, pmax, NULL, 0, psi, npsi, t - pmod->t2); 
		fc->sderr[t] = pmod->sigma * sqrt(vl);
	    } else {
		fc->sderr[t] = NADBL;
	    }
	}
    }

    if (psi != NULL) {
	free(phi);
	free(psi);
	free(errphi);
    }

    return err;
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

static int linear_fcast (Forecast *fc, const MODEL *pmod, int yno,
			 const double **Z, const DATAINFO *pdinfo,
			 gretlopt opt)
{
    const double *offvar = NULL;
    double lmax = NADBL;
    double xval, yht;
    int i, vi, t;

    if (pmod->ci == POISSON) {
	/* special for poisson "offset" variable */
	int offnum = gretl_model_get_int(pmod, "offset_var");

	if (offnum > 0) {
	    offvar = Z[offnum];
	}
    } else if (pmod->ci == LOGISTIC) {
	lmax = gretl_model_get_double(pmod, "lmax");
    }

    if (opt & OPT_I) {
	yht = Z[yno][fc->t1 - 1];
    }

    for (t=fc->t1; t<=fc->t2; t++) {
	int miss = 0;

	if (!(opt & OPT_I)) {
	    /* not integrating forecast */
	    yht = 0.0;
	}

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    int lag;

	    vi = pmod->list[i+2];
	    if ((lag = depvar_lag(fc, i))) {
		xval = fcast_get_ldv(fc, yno, t, lag, Z);
	    } else {
		xval = Z[vi][t];
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		yht += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    if (opt & OPT_I) {
		/* integrating: all subsequent values are NA */
		int s;

		for (s=t; s<=fc->t2; s++) {
		    fc->yhat[s] = NADBL;
		}
		break;
	    } else {
		fc->yhat[t] = NADBL;
	    }
	} else if (FCAST_SPECIAL(pmod->ci)) {
	    /* special handling for LOGIT and others */
	    fc->yhat[t] = fcast_transform(yht, pmod->ci, t, offvar, lmax);
	} else {
	    fc->yhat[t] = yht;
	}
    }

    return 0;
}

#define dynamic_nls(m) (m->ci == NLS && gretl_model_get_int(m, "dynamic"))

static int get_forecast_method (Forecast *fc,
				MODEL *pmod, 
				const DATAINFO *pdinfo,
				gretlopt opt)
{
    int dyn_ok = 0;
    int dyn_errs_ok = 0;
    int err;

    err = incompatible_options(opt, OPT_D | OPT_S);
    if (err) {
	return err;
    }

    fc->dvlags = NULL;
    fc->method = FC_STATIC;

    /* do setup for possible lags of the dependent variable,
       unless OPT_S for "static" has been given */
    if (dataset_is_time_series(pdinfo) && !(opt & OPT_S) &&
	pmod->ci != ARMA) {
	process_lagged_depvar(pmod, pdinfo, &fc->dvlags);
    }

    if (!(opt & OPT_S)) {
	/* user didn't give the "static" option */
	if (pmod->ci == ARMA || fc->dvlags != NULL || dynamic_nls(pmod)) {
	    /* dynamic forecast is possible */
	    dyn_ok = 1;
	}
	if (SIMPLE_AR_MODEL(pmod->ci) || pmod->ci == GARCH) {
	    dyn_errs_ok = 1;
	}
    } 

    if (!dyn_ok && (opt & OPT_D)) {
	/* "dynamic" option given, but can't be honored */
	fputs("requested dynamic option, but it is not applicable\n", stderr);
	return E_NOTIMP;
    }

    /* NLS: we can only do dynamic out of sample (fc->t1 > pmod->t2) */
    if (pmod->ci == NLS && dyn_ok && (opt & OPT_D)) {
	if (fc->t1 <= pmod->t2) {
	    return E_NOTIMP;
	}
    }

    if (opt & OPT_D) {
	/* user requested dynamic forecast and it seems OK */
	fc->method = FC_DYNAMIC;
    } else if ((dyn_ok || dyn_errs_ok) && fc->t2 > pmod->t2) {
	/* do dynamic f'cast out of sample */
	fc->method = FC_AUTO;
    }

    return err;
}

static void forecast_init (Forecast *fc)
{
    fc->method = FC_AUTO;

    fc->t1 = 0;
    fc->t2 = 0;
    fc->model_t2 = 0;

    fc->yhat = NULL;
    fc->sderr = NULL;
    fc->eps = NULL;
    fc->dvlags = NULL;
}

static void forecast_free (Forecast *fc)
{
    if (fc->dvlags != NULL) {
	free(fc->dvlags);
    }

    if (fc->eps != NULL) {
	free(fc->eps);
    }
} 

static int fc_add_eps (Forecast *fc, int n)
{
    int t, err = 0;

    fc->eps = malloc(n * sizeof *fc->eps);

    if (fc->eps == NULL) {
	err = E_ALLOC;
    } else {
	for (t=0; t<n; t++) {
	    fc->eps[t] = 0.0;
	}
    }

    return err;
}

/* driver for various functions that compute forecasts
   for different sorts of models */

static int real_get_fcast (FITRESID *fr, MODEL *pmod, 
			   double ***pZ, DATAINFO *pdinfo,
			   gretlopt opt) 
{
    Forecast fc;
    const double **Z = (const double **) *pZ;
    int yno = gretl_model_get_depvar(pmod);
    int dummy_AR = 0;
    int DM_errs = 0;
    int dyn_errs = 0;
    int asy_errs = 0;
    int nf = 0;
    int t, err;

    forecast_init(&fc);

    fc.t1 = fr->t1;
    fc.t2 = fr->t2;
    fc.model_t2 = pmod->t2;

    err = get_forecast_method(&fc, pmod, pdinfo, opt);
    if (err) {
	return err;
    }

    if (opt & OPT_I) {
	/* integrated forecast wanted */
	int d, parent;

	d = is_standard_diff(yno, pdinfo, &parent);
	if (d == 0) {
	    return E_DATA;
	} 
	yno = parent;
    }

    if (pmod->ci == NLS && fc.method == FC_STATIC) {
	asy_errs = 1;
    }

    if (!FCAST_SPECIAL(pmod->ci)) {
	if (!AR_MODEL(pmod->ci) && fc.dvlags == NULL) {
	    /* we'll do Davidson-MacKinnon error variance */
	    DM_errs = 1;
	} else if (fc.method == FC_DYNAMIC) {
	    /* we'll do dynamic forecast errors throughout */
	    dyn_errs = 1;
	} else if (fc.method == FC_AUTO && fc.t2 > pmod->t2) {
	    /* do dynamic forecast errors out of sample */
	    dyn_errs = 1;
	}
    }

    /* FIXME temporary */
    if (DM_errs && (opt & OPT_I)) {
	DM_errs = 0;
    }

    if (dyn_errs && !AR_MODEL(pmod->ci) && fc.dvlags != NULL) {
	/* create dummy AR info structure for model with lagged
	   dependent variable */
	int err = dummy_ar_info_init(pmod);

	if (err) {
	    dyn_errs = 0;
	} else {
	    dummy_AR = 1;
	}
    }

    if (DM_errs || dyn_errs || asy_errs) {
	err = fit_resid_add_sderr(fr);
    }

    if (err) {
	return err;
    }

    fc.yhat = fr->fitted;
    fc.sderr = fr->sderr;

    if (pmod->ci == ARMA && fc.method == FC_STATIC) {
	fc_add_eps(&fc, pdinfo->n);
    }

    /* compute the actual forecast */
    if (DM_errs) {
	err = static_fcast_with_errs(&fc, pmod, Z, pdinfo);
    } else if (pmod->ci == NLS) {
	err = nls_fcast(&fc, pmod, pZ, pdinfo);
    } else if (SIMPLE_AR_MODEL(pmod->ci) || dummy_AR) {
	err = ar_fcast(&fc, pmod, Z, pdinfo);
    } else if (pmod->ci == ARMA) {
	err = arma_fcast(&fc, pmod, Z, pdinfo);
    } else if (pmod->ci == GARCH) {
	err = garch_fcast(&fc, pmod, Z, pdinfo);
    } else {
	err = linear_fcast(&fc, pmod, yno, Z, pdinfo, opt);
    }

    /* free any auxiliary info */
    forecast_free(&fc);

    if (dummy_AR) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
    }

    if (opt & OPT_I) {
	t = fr->t0 - 1;
	fprintf(stderr, "t-prior = %d\n", t);
	fr->fitted[t] = (*pZ)[yno][t];
	fr->resid[t] = 0.0;
    }

    for (t=0; t<fr->nobs; t++) {
	if (t >= fr->t1 && t <= fr->t2) {
	    if (!na(fr->fitted[t])) {
		nf++;
	    }	    
	} else if (t >= fr->t0 && t >= pmod->t1 && t <= pmod->t2) {
	    if (opt & OPT_I) {
		fr->fitted[t] = fr->fitted[t-1] + pmod->yhat[t];
		fr->resid[t] = fr->resid[t-1] + pmod->uhat[t];
	    } else {
		fr->fitted[t] = pmod->yhat[t];
		fr->resid[t] = pmod->uhat[t];
	    }
	}	    
	fr->actual[t] = (*pZ)[yno][t];
    }

    if (nf == 0) {
	err = E_MISSDATA;
    } else {
	fit_resid_set_dec_places(fr);
	strcpy(fr->depvar, pdinfo->varname[yno]);
	fr->df = pmod->dfd;
    }

    return err;
}

static int fcast_get_limit (const char *s, double ***pZ, DATAINFO *pdinfo)
{
    double x;
    int t, err = 0;

    if (gretl_is_scalar(s)) {
	x = gretl_scalar_get_value(s);
    } else {
	x = generate_scalar(s, pZ, pdinfo, &err);
    }

    if (x < 1 || x > pdinfo->n) {
	gretl_errmsg_set(_("Observation number out of bounds"));
	t = -1;
    } else {
	t = x - 1;
    }

    return t;
}

static int parse_forecast_string (const char *s, 
				  gretlopt opt, 
				  int t2est,
				  double ***pZ,
				  DATAINFO *pdinfo,
				  int *pt1, int *pt2,
				  int *pk, char *vname)
{
    char f[4][32];
    char *t1str = NULL, *t2str = NULL;
    char *kstr = NULL, *vstr = NULL;
    int nmax, nmin = (opt & OPT_R)? 1 : 0;
    int t1 = 0, t2 = 0;
    int nf, err = 0;
    
    /* "static" and "rolling" can't be combined */
    err = incompatible_options(opt, OPT_S | OPT_R);
    if (err) {
	return err;
    }

    if (!strncmp(s, "fcasterr", 8)) {
        s += 8;
    } else if (!strncmp(s, "fcast", 5)) {
	s += 5;
    }

    *vname = '\0';

    /* How many fields should we be looking for in the user input?
       If OPT_R ("rolling") is given, the max is 4: 

       t1, t2, k (steps-ahead), vname

       Otherwise the max is 3 (no k-value):

       t1, t2, vname

       However, if OPT_O ("out-of-sample") is given, we should not 
       expect t1 and t2, so the max becomes 2 (if OPT_R), else 1.

       Also note: t1 and t2 need not be given (even in the absence
       if OPT_O) since these default to the current sample range.
    */

    if (opt & OPT_O) {
	nmax = (opt & OPT_R)? 2 : 1;
    } else {
	nmax = (opt & OPT_R)? 4 : 3;
    }    

    nf = sscanf(s, "%31s %31s %31s %31s", f[0], f[1], f[2], f[3]);

    if (nf < nmin || nf > nmax) {
	fprintf(stderr, "fcast: expected %d to %d fields in input, got %d\n",
		nmin, nmax, nf);
	return E_DATA;
    }

    if (opt & OPT_R) {
	if (nf == 4) {
	    /* t1, t2, k, vname */
	    t1str = f[0];
	    t2str = f[1];
	    kstr = f[2];
	    vstr = f[3];
	} else if (nf == 3) {
	    /* t1, t2, k */
	    t1str = f[0];
	    t2str = f[1];
	    kstr = f[2];
	} else if (nf == 2) {
	    /* k, vname */
	    kstr = f[0];
	    vstr = f[1];
	} else if (nf == 1) {
	    /* just k */
	    kstr = f[0];
	}
    } else {
	if (nf == 3) {
	    /* t1, t2, vname */
	    t1str = f[0];
	    t2str = f[1];
	    vstr = f[2];
	} else if (nf == 2) {
	    /* t1, t2 */
	    t1str = f[0];
	    t2str = f[1];
	} else if (nf == 1) {
	    /* just vname */
	    vstr = f[0];
	}
    }

    if ((opt & OPT_O) && (t1str != NULL || t2str != NULL)) {
	fprintf(stderr, "fcast: got unexpected t1 and/or t2 field in input\n");
	return E_DATA;
    }	

    if (kstr != NULL && pk != NULL) {
	*pk = positive_int_from_string(kstr);
    }

    if (vstr != NULL) {
	strncat(vname, vstr, VNAMELEN - 1);
    }

    if (t1str != NULL && t2str != NULL) {
	t1 = dateton(t1str, pdinfo);
	if (t1 < 0) {
	    t1 = fcast_get_limit(t1str, pZ, pdinfo);
	}
	t2 = dateton(t2str, pdinfo);
	if (t2 < 0) {
	    t2 = fcast_get_limit(t2str, pZ, pdinfo);
	}	
	if (t1 < 0 || t2 < 0 || t2 < t1) {
	    err = E_DATA;
	}
    } else if (opt & OPT_O) {
	/* out of sample, if possible */
	if (pdinfo->n - t2est - 1 > 0) {
	    t1 = t2est + 1;
	    t2 = pdinfo->n - 1;
	} else {
	    err = E_OBS;
	}
    } else {
	/* default: sample range */
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    } 

    if (!err) {
	*pt1 = t1;
	*pt2 = t2;
    }

    return err;
}

/**
 * get_forecast:
 * @pmod: the model from which forecasts are wanted.
 * @t1: start of forecast range.
 * @t2: end of forecast range.
 * @pre_n: number of pre-forecast observations to include.
 * @pZ: pointer to data array using which @pmod was estimated.
 * @pdinfo: dataset information.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where a dynamic forecast is meaningful).
 * @err: location to receive error code.
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
 * Returns: pointer to allocated structure, or %NULL on failure,
 * in which case an error code is assigned via @err.
 */

FITRESID *get_forecast (MODEL *pmod, int t1, int t2, int pre_n,
			double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, int *err) 
{
    FITRESID *fr;

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	*err = E_DATA;
	return NULL;
    }

    fr = fit_resid_new_for_model(pmod, pdinfo, t1, t2, pre_n, err);

    if (!*err) {
	*err = real_get_fcast(fr, pmod, pZ, pdinfo, opt);
	if (*err) {
	    free_fit_resid(fr);
	    fr = NULL;
	}
    }

    return fr;
}

static gretl_matrix *fcast_matrix;
static gretl_matrix *fcerr_matrix;

void forecast_matrix_cleanup (void)
{
    if (fcast_matrix != NULL) {
	gretl_matrix_free(fcast_matrix);
	fcast_matrix = NULL;
    }

    if (fcerr_matrix != NULL) {
	gretl_matrix_free(fcerr_matrix);
	fcerr_matrix = NULL;
    }
}    

gretl_matrix *get_forecast_matrix (int idx, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *M = NULL;

    if (idx == M_FCAST) {
	M = fcast_matrix;
    } else if (idx == M_FCERR) {
	M = fcerr_matrix;
    } 

    if (M == NULL) {
	*err = E_BADSTAT;
    } else {
	ret = gretl_matrix_copy(M);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

static void fr_adjust_sample (FITRESID *fr, int *ft1, int *ft2,
			      int *et1, int *et2)
{
    int t;

    t = fr->t0;
    while (t<=fr->t2 && na(fr->fitted[t])) {
	*ft1 = ++t;
    }

    t = fr->t2;
    while (t>=*ft1 && na(fr->fitted[t])) {
	*ft2 = --t;
    }

    if (fr->sderr != NULL) {
	t = fr->t0;
	while (t<=fr->t2 && na(fr->sderr[t])) {
	    *et1 = ++t;
	}

	t = fr->t2;
	while (t>=*et1 && na(fr->sderr[t])) {
	    *et2 = --t;
	}
    }	
}

static int set_forecast_matrices_from_fr (FITRESID *fr)
{
    gretl_matrix *f = NULL;
    gretl_matrix *e = NULL;
    int ft1 = fr->t0;
    int ft2 = fr->t2;
    int et1 = fr->t0;
    int et2 = fr->t2;
    int T, t, s;
    int err = 0;

    fr_adjust_sample(fr, &ft1, &ft2, &et1, &et2);

    T = ft2 - ft1 + 1;
    if (T <= 0) {
	return 0;
    }

    f = gretl_matrix_alloc(T, 1);
    if (f == NULL) {
	return E_ALLOC;
    }

    f->t1 = ft1;
    f->t2 = ft2;

    if (fr->sderr != NULL) {
	T = et2 - et1 + 1;
	if (T > 0) {
	    e = gretl_matrix_alloc(T, 1);
	    if (e == NULL) {
		err = E_ALLOC;
	    } else {
		e->t1 = et1;
		e->t2 = et2;
	    }
	}	    
    }

    s = 0;
    for (t=ft1; t<=ft2; t++) {
	f->val[s++] = fr->fitted[t];
    }

    if (e != NULL) {
	s = 0;
	for (t=et1; t<=et2; t++) {
	    e->val[s++] = fr->sderr[t];
	}
    }

    gretl_matrix_free(fcast_matrix);
    gretl_matrix_free(fcerr_matrix);

    fcast_matrix = f;
    fcerr_matrix = e;

    return err;
}

static int matrix_first_ok_row (const gretl_matrix *a, int jmin, int cols)
{
    int jmax = jmin + cols;
    double x;
    int i, j, miss;

    for (i=0; i<a->rows; i++) {
	miss = 0;
	for (j=jmin; j<jmax && !miss; j++) {
	    x = gretl_matrix_get(a, i, j);
	    if (na(x)) {
		miss = 1;
	    }
	}
	if (!miss) {
	    break;
	}
    }

    return i;
}

static int matrix_last_ok_row (const gretl_matrix *a, int jmin, int cols)
{
    int jmax = jmin + cols;
    double x;
    int i, j, miss;

    for (i=a->rows-1; i>=0; i--) {
	miss = 0;
	for (j=jmin; j<jmax && !miss; j++) {
	    x = gretl_matrix_get(a, i, j);
	    if (na(x)) {
		miss = 1;
	    }
	}
	if (!miss) {
	    break;
	}
    }

    return i;
}

static void F_matrix_adjust_sample (const gretl_matrix *F, 
				    int *f0, int *fn, 
				    int *e0, int *en)
{
    int k = F->cols / 2;

    *f0 = matrix_first_ok_row(F, 0, k);
    *fn = matrix_last_ok_row(F, 0, k);

    *e0 = matrix_first_ok_row(F, k, k);
    *en = matrix_last_ok_row(F, k, k);
}

static int set_forecast_matrices_from_F (const gretl_matrix *F,
					 int imin, int imax)
{
    gretl_matrix *F1 = NULL;
    gretl_matrix *f = NULL;
    gretl_matrix *e = NULL;
    int n = F->cols;
    int k = n / 2;
    int fT = F->rows;
    int eT = F->rows;
    int f0 = 0, fn = fT;
    int e0 = 0, en = eT;
    double x;
    int i, j;
    int err = 0;

    if (imin == imax) {
	/* extract one pair of columns */
	F1 = gretl_matrix_alloc(fT, 2);
	if (F1 == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<fT; i++) {
	    x = gretl_matrix_get(F, i, imin);
	    gretl_matrix_set(F1, i, 0, x);
	    x = gretl_matrix_get(F, i, imin + k);
	    gretl_matrix_set(F1, i, 1, x);
	}
	F = F1;
	k = 1;
    }

    F_matrix_adjust_sample(F, &f0, &fn, &e0, &en);
    
    fT = fn - f0 + 1;
    if (fT <= 0) {
	return E_MISSDATA;
    }

    f = gretl_matrix_alloc(fT, k);
    if (f == NULL) {
	return E_ALLOC;
    }

    f->t1 = F->t1 + f0;
    f->t2 = f->t1 + fT - 1;

    eT = en - e0 + 1;

    if (eT > 0) {
	e = gretl_matrix_alloc(eT, k);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    e->t1 = F->t1 + e0;
	    e->t2 = e->t1 + eT - 1;
	}	    
    }

    for (j=0; j<k; j++) {
	for (i=0; i<fT; i++) {
	    x = gretl_matrix_get(F, i + f0, j);
	    gretl_matrix_set(f, i, j, x);
	}
    }

    if (e != NULL) {
	for (j=0; j<k; j++) {
	    for (i=0; i<eT; i++) {
		x = gretl_matrix_get(F, i + e0, k + j);
		gretl_matrix_set(e, i, j, x);
	    }
	}
    }    

    gretl_matrix_free(fcast_matrix);
    gretl_matrix_free(fcerr_matrix);

    fcast_matrix = f;
    fcerr_matrix = e;

    if (F1 != NULL) {
	gretl_matrix_free(F1);
    }

    return err;
}

/* compatibility with behaviour of old fcast command */

static int add_fcast_to_dataset (FITRESID *fr, const char *vname,
				 double ***pZ, DATAINFO *pdinfo,
				 PRN *prn)
{
    int oldv = pdinfo->v;
    int v, err = 0;

    v = series_index(pdinfo, vname);

    if (v == pdinfo->v) {
	/* new variable */
	if (check_varname(vname)) {
	    err = E_DATA;
	} else {
	    err = dataset_add_series(1, pZ, pdinfo);
	    if (!err) {
		strcpy(pdinfo->varname[v], vname);
	    }
	}
    }

    if (!err) {
	int t;

	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[v][t] = fr->fitted[t];
	}

	strcpy(VARLABEL(pdinfo, v), _("predicted values"));

	if (gretl_messages_on()) {
	    if (v < oldv) {
		pprintf(prn, _("Replaced series %s (ID %d)"),
			vname, v);
	    } else {
		pprintf(prn, _("Generated series %s (ID %d)"),
			vname, v);
	    }
	    pputc(prn, '\n');
	}
    }

    return err;
}

static int model_do_forecast (const char *str, MODEL *pmod, 
			      double ***pZ, DATAINFO *pdinfo, 
			      gretlopt opt, PRN *prn)
{
    char vname[VNAMELEN];
    FITRESID *fr;
    int t1, t2, k = -1;
    int err;

    if (pmod->ci == ARBOND || pmod->ci == HECKIT) {
	/* FIXME */
	return E_NOTIMP;
    }

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	return E_DATA;
    }

    err = parse_forecast_string(str, opt, pmod->t2, pZ, pdinfo, 
				&t1, &t2, &k, vname);
    if (err) {
	return err;
    }

    if (*vname != '\0') {
	/* saving to named series */
	opt |= (OPT_A | OPT_Q);
    }

    if (opt & OPT_R) {
	fr = rolling_OLS_k_step_fcast(pmod, pZ, pdinfo, t1, t2, 
				      k, 0, &err);
    } else {
	fr = get_forecast(pmod, t1, t2, 0, pZ, pdinfo, opt, &err);
    }

    if (!err && !(opt & OPT_Q)) {
	err = text_print_forecast(fr, pdinfo, opt, prn);
    }

    if (!err && (opt & OPT_A)) {
	/* add forecast directly */
	err = add_fcast_to_dataset(fr, vname, pZ, pdinfo, prn);
    }

    if (!err) {
	set_forecast_matrices_from_fr(fr);
    }

    free_fit_resid(fr);

    return err;
}

static int get_sys_fcast_var (const int *ylist, const char *vname,
			      DATAINFO *pdinfo)
{
    int i, vi;

    for (i=0; i<ylist[0]; i++) {
	vi = ylist[i+1];
	if (!strcmp(vname, pdinfo->varname[vi])) {
	    return i;
	}
    }

    return -1;
}

/* grab the relevant range from system forecast matrix and write it
   into a FITRESID struct for printing or graphing
*/

static int fill_system_forecast (FITRESID *fr, int i, int yno,
				 GRETL_VAR *var, equation_system *sys,
				 const gretl_matrix *F,
				 const double **Z, DATAINFO *pdinfo,
				 gretlopt opt)
{
    int m = F->cols / 2;
    int s, t, nf;
    int err = 0;

    strcpy(fr->depvar, pdinfo->varname[yno]);

    /* "pre-forecast" observations */
    for (t=fr->t0; t<fr->t1; t++) {
	fr->actual[t] = Z[yno][t];
	if (sys != NULL) {
	    if (i < sys->neqns && sys->lists[i][1] == yno) {
		/* Note: right now we can only handle endogenous
		   variables that appear on the LHS of stochastic
		   equations, because only such variables have an
		   associated column in sys->yhat.
		*/
		if (t >= sys->t1 && t <= sys->t2) {
		    s = t - sys->t1;
		    fr->fitted[t] = gretl_matrix_get(sys->yhat, s, i);
		}
	    }
	} else if (var != NULL) {
	    if (t >= var->t1 && t <= var->t2) {
		s = t - var->t1;
		fr->fitted[t] = fr->actual[t] - gretl_matrix_get(var->E, s, i);
	    }
	}
    }

    /* actual forecasts */
    nf = 0;
    for (t=fr->t1, s=0; t<=fr->t2; t++, s++) {
	fr->actual[t] = Z[yno][t];
	fr->fitted[t] = gretl_matrix_get(F, s, i);
	if (!na(fr->fitted[t])) {
	    nf++;
	}
	if (fr->sderr != NULL) {
	    fr->sderr[t] = gretl_matrix_get(F, s, i + m);
	}
    }

    if (nf == 0) {
	err = E_MISSDATA;
    } else {
	fit_resid_set_dec_places(fr);
    }

    return err;
}

static int system_do_forecast (const char *str, void *ptr, int type,
			       double ***pZ, DATAINFO *pdinfo, 
			       gretlopt opt, PRN *prn)
{
    char vname[VNAMELEN];
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const double **Z;
    const int *ylist = NULL;
    const gretl_matrix *F = NULL;
    int t1 = 0, t2 = 0;
    int t2est, ci;
    int imax, imin = 0;
    int df = 0;
    int err = 0;

    if (type == GRETL_OBJ_VAR) {
	var = (GRETL_VAR *) ptr;
	imax = var->neqns - 1;
	t2est = var->t2;
	ci = var->ci;
	df = var->df;
	ylist = var->ylist;
    } else {
	sys = (equation_system *) ptr;
	imax = sys->neqns + sys->nidents - 1;
	t2est = sys->t2;
	ci = SYSTEM;
	df = sys->df;
	ylist = sys->ylist;
    }

    err = parse_forecast_string(str, opt, t2est, pZ, pdinfo, &t1, &t2, NULL, vname);

    Z = (const double **) *pZ;

    if (!err) {
	if (var != NULL) {
	    F = gretl_VAR_get_forecast_matrix(var, t1, t2, Z, pdinfo, 
					      opt, &err);
	} else {
	    F = system_get_forecast_matrix(sys, t1, t2, Z, pdinfo, 
					   opt, &err);
	} 
    }
    
    if (!err && *vname != '\0') {
	imin = get_sys_fcast_var(ylist, vname, pdinfo);
	if (imin < 0) {
	    err = E_DATA;
	} else {
	    imax = imin;
	}
    }

    if (!err) {
	/* arrange to save forecast and errors */
	err = set_forecast_matrices_from_F(F, imin, imax);
    }	

    if (err) {
	return err;
    }

    if (!(opt & OPT_Q)) {
	/* assemble and print per-equation forecasts */
	gretlopt printopt = OPT_NONE;
	FITRESID *fr = NULL;
	int i, asy;

	asy = (ci == VECM);

	fr = fit_resid_new_for_system(asy, pdinfo, t1, t2, 0, &err);
	if (err) {
	    return err;
	}

    	if (asy) {
	    /* asymptotic normal */
	    fr->df = var->T;
	} else {
	    fr->df = df;
	}
    
	for (i=imin; i<=imax && !err; i++) {
	    err = fill_system_forecast(fr, i, ylist[i+1], var, sys, 
				       F, Z, pdinfo, opt);
	    if (!err) {
		err = text_print_forecast(fr, pdinfo, printopt, prn);
	    }
	    printopt |= OPT_Q;
	}

	free_fit_resid(fr);
    }

    return err;
}

/**
 * do_forecast:
 * @str: command string, which may include a starting 
 * observation and ending observation, and/or the name of a 
 * variable for saving the forecast values.
 * @pZ: pointer to data matrix.
 * @pdinfo: pointer to data information struct.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where this distinction is meaningful).
 * %OPT_R: do rolling/recursive forecast.
 * %OPT_Q: suppress printing of the forecast;
 * %OPT_P: ensure that the values are printed.
 *
 * In the case of "simple" models with an autoregressive error term 
 * (%AR, %AR1) the predicted values incorporate the forecastable portion 
 * of the error.  
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int do_forecast (const char *str, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    void *ptr;
    GretlObjType type;
    int err;

    ptr = get_last_model(&type);
    if (ptr == NULL) {
	return E_BADSTAT;
    }

    if ((opt & OPT_R) && type != GRETL_OBJ_EQN) {
	/* "rolling": single equations only */
	err = E_NOTIMP;
    } else if (type == GRETL_OBJ_EQN) {
	err = model_do_forecast(str, ptr, pZ, pdinfo, opt, prn);
    } else if (type == GRETL_OBJ_SYS || type == GRETL_OBJ_VAR) {
	err = system_do_forecast(str, ptr, type, pZ, pdinfo, opt, prn);
    } else {
	err = E_DATA;
    }

    return err;
}

/* try to determine in advance how far we can go with a forecast,
   either dynamic or static (ftype) */

static int 
fcast_get_t2max (const int *list, const int *dvlags, const MODEL *pmod,
		 const double **Z, const DATAINFO *pdinfo, int ftype)
{
    const double *ay = NULL;
    int i, vi, t;

    if (pmod->ci == ARMA && ftype == FC_STATIC) {
	int yno = gretl_model_get_depvar(pmod);

	ay = Z[yno];
    }

    for (t=pmod->t2; t<pdinfo->n; t++) {
	int all_ok = 1;

	if (ay != NULL && na(ay[t-1])) {
	    /* FIXME? */
	    break;
	}

	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == 0) {
		continue;
	    } else if (dvlags != NULL && dvlags[i-1] != 0) {
		continue;
	    } else if (is_trend_variable(Z[vi], pdinfo->n)) {
		continue;
	    } else if (is_periodic_dummy(Z[vi], pdinfo)) {
		continue;
	    }
	    if (na(Z[vi][t])) {
		all_ok = 0;
		break;
	    }
	}

	if (!all_ok) {
	    t--;
	    break;
	} else if (t == pdinfo->n - 1) {
	    break;
	}
    }

    return t;
}

/**
 * get_system_forecast:
 * @p: pointer to the VAR or equation system from which 
 * forecasts are wanted.
 * @ci: command index for system (%VAR, %VECM or %SYSTEM)
 * @i: 0-based index for the variable to forecast, within
 * the equation system.
 * @t1: start of forecast range.
 * @t2: end of forecast range.
 * @pre_n: number of pre-forecast observations to include.
 * @Z: data array on which system was estimated.
 * @pdinfo: dataset information.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample.
 * @err: location to receive error code.
 *
 * Allocates a #FITRESID structure and fills it out with forecasts
 * based on the system at location @p, over the specified range of 
 * observations.
 * 
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FITRESID *get_system_forecast (void *p, int ci, int i, 
			       int t1, int t2, int pre_n,
			       const double **Z, DATAINFO *pdinfo,
			       gretlopt opt, int *err) 
{
    FITRESID *fr;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *F = NULL;
    int nf = t2 - t1 + 1;
    int asy, df = 0, yno = 0;

    if (nf <= 0) {
	*err = E_DATA;
	return NULL;
    }

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) p;
	yno = var->ylist[i+1];
	df = var->df;
	F = gretl_VAR_get_forecast_matrix(var, t1, t2, Z, pdinfo, opt, err);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) p;
	yno = sys->ylist[i+1];
	df = sys->df;
	F = system_get_forecast_matrix(sys, t1, t2, Z, pdinfo, opt, err);
    } else {
	*err = E_DATA;
    }

    if (*err) {
	fprintf(stderr, "get_system_forecast: matrix F is NULL\n");
	return NULL;
    }

    asy = (ci == VECM);

    fr = fit_resid_new_for_system(asy, pdinfo, t1, t2, pre_n, err);
    if (*err) {
	return NULL;
    }
    
    if (asy) {
	/* asymptotic normal */
	fr->df = var->T;
    } else {
	fr->df = df;
    }

    *err = fill_system_forecast(fr, i, yno, var, sys, 
				F, Z, pdinfo, opt);

    if (*err) {
	free_fit_resid(fr);
	fr = NULL;
    } 

    return fr;
}

/**
 * forecast_options_for_model:
 * @pmod: the model from which forecasts are wanted.
 * @Z: data array.
 * @pdinfo: dataset information.
 * @flags: location to receive flags from among #ForecastFlags.
 * @dt2max: location to receive the last observation that can
 * be supported for a dynamic forecast.
 * @st2max: location to receive the last observation that can
 * be supported for a static forecast.
 *
 * Examines @pmod and determines which forecasting options are
 * applicable.
 */

void forecast_options_for_model (MODEL *pmod, const double **Z,
				 const DATAINFO *pdinfo, int *flags, 
				 int *dt2max, int *st2max)
{
    int *dvlags = NULL;
    int dv, exo = 1;

    *flags = 0;

    dv = gretl_model_get_depvar(pmod);

    if (is_standard_diff(dv, pdinfo, NULL)) {
	*flags |= FC_INTEGRATE_OK;
    }

    if (pmod->ci == NLS) {
	/* we'll try winging it! */
	if (gretl_model_get_int(pmod, "dynamic") && pmod->t2 < pdinfo->n - 1) {
	    *flags |= FC_AUTO_OK;
	}
	return;
    }

    *dt2max = pmod->t2;
    *st2max = pmod->t2;

    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	*flags |= FC_DYNAMIC_OK;
    } else if (AR_MODEL(pmod->ci)) {
	*flags |= FC_DYNAMIC_OK;
    } else if (dataset_is_time_series(pdinfo) &&
	       has_depvar_lags(pmod, pdinfo)) {
	*flags |= FC_DYNAMIC_OK;
    }

    if (*flags & FC_DYNAMIC_OK) {
	int err = process_lagged_depvar(pmod, pdinfo, &dvlags);

	if (!err) {
	    exo = has_real_exog_regressors(pmod, dvlags, Z, pdinfo);
	}
	if (!exo) {
	    *flags |= FC_ADDOBS_OK;
	    *dt2max = pdinfo->n - 1;
	}
    } 

    if (exo) {
	const int *xlist = model_xlist(pmod);

	if (xlist != NULL) {
	    *dt2max = fcast_get_t2max(xlist, dvlags, pmod, Z, pdinfo, FC_DYNAMIC);
	    *st2max = fcast_get_t2max(xlist, dvlags, pmod, Z, pdinfo, FC_STATIC);
	}
    }

    if (dvlags != NULL) {
	free(dvlags);
    }
}

static int y_lag (int v, int parent, const DATAINFO *pdinfo)
{
    if (pdinfo->varinfo[v]->transform == LAGS) {
	int pv = series_index(pdinfo, pdinfo->varinfo[v]->parent);

	if (pv == parent) {
	    return pdinfo->varinfo[v]->lag;
	}
    }

    return 0;
}

static int k_step_init (MODEL *pmod, const DATAINFO *pdinfo, 
			double **Z, int k, 
			double **py, int **pllist)
{
    double *y = NULL;
    int *llist = NULL;
    int vy = pmod->list[1];
    int i, nl = 0;

    for (i=2; i<=pmod->list[0]; i++) {
	if (y_lag(pmod->list[i], vy, pdinfo)) {
	    nl++;
	}
    }

    if (nl == 0) {
	return 0;
    }

    y = malloc(pdinfo->n * sizeof *y);
    llist = gretl_list_new(pmod->list[0] - 1);

    if (y == NULL || llist == NULL) {
	free(y);
	free(llist);
	return E_ALLOC;
    }

    for (i=0; i<pdinfo->n; i++) {
	y[i] = Z[vy][i];
    }

    for (i=2; i<=pmod->list[0]; i++) {
	llist[i-1] = y_lag(pmod->list[i], vy, pdinfo);
    }    

    *py = y;
    *pllist = llist;

    return 0;
}

static int rolling_fcast_adjust_obs (MODEL *pmod, int *t1, int t2, int k)
{
    /* the earliest possible forecast start */
    int t1min = pmod->t1 + pmod->ncoeff + k - 1;
    int err = 0;

    if (*t1 < t1min) {
	*t1 = t1min;
    }

    /* minimal forecast range = 1, otherwise error */
    if (t2 < *t1) {
	err = E_OBS;
    }

    return err;
}

/* recursive k-step ahead forecasts, for models estimated via OLS */

FITRESID * 
rolling_OLS_k_step_fcast (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			  int t1, int t2, int k, int pre_n, int *err)
{
    FITRESID *fr;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    double *y = NULL;
    int *llist = NULL;
    double xit, yf = NADBL;
    MODEL mod;
    int j, p, vi, nf;
    int i, s, t;

    if (pmod->ci != OLS) {
	*err = E_OLSONLY;
	return NULL;
    }

    if (k < 1) {
	gretl_errmsg_set("recursive forecast: steps-ahead must be >= 1");
	*err = E_DATA;
	return NULL;
    }    

    if (gretl_model_get_int(pmod, "daily_repack")) {
	*err = E_DATA;
	return NULL;
    }

    /* check feasibility of forecast range */
    *err = rolling_fcast_adjust_obs(pmod, &t1, t2, k);
    if (*err) {
	return NULL;
    }

    if (k > 1) {
	*err = k_step_init(pmod, pdinfo, *pZ, k, &y, &llist);
	if (*err) {
	    return NULL;
	}
    }

    fr = fit_resid_new_for_model(pmod, pdinfo, t1, t2, pre_n, err); 
    if (*err) {
	free(y);
	free(llist);
	return NULL;
    }

    fr->method = FC_KSTEP;
    fr->k = k;

    /* sample range for initial estimation */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = t1 - k; /* start of fcast range minus k */

    /* number of forecasts */
    nf = t2 - t1 + 1;

    fprintf(stderr, "rolling fcast: pdinfo->t1=%d, pdinfo->t2=%d, t1=%d, t2=%d, k=%d, nf=%d\n",
	    pdinfo->t1, pdinfo->t2, t1, t2, k, nf);

    for (t=0; t<pdinfo->n; t++) {
	fr->actual[t] = (*pZ)[pmod->list[1]][t];
    }

    for (s=0; s<nf; s++) {
	mod = lsq(pmod->list, pZ, pdinfo, OLS, OPT_A | OPT_Z);
	if (mod.errcode) {
	    *err = mod.errcode;
	    clear_model(&mod);
	    break;
	}

	/* the first obs following the estimation sample */
	t = pdinfo->t2 + 1;

	for (j=0; j<k; j++) {
	    yf = 0.0;
	    /* steps ahead */
	    for (i=0; i<mod.ncoeff; i++) {
		vi = mod.list[i+2];
		p = (llist != NULL)? llist[i+1] : 0;
		if (p > 0 && p <= j) {
		    xit = y[t-p];
		} else {
		    xit = (*pZ)[vi][t];
		}
		if (na(xit)) {
		    yf = NADBL;
		    break;
		} else {
		    yf += mod.coeff[i] * xit;
		}
	    }
	    if (y != NULL && j < k - 1) {
		y[t++] = yf;
	    }
	}

	fr->fitted[t] = yf;

	if (!na(fr->actual[t]) && !na(fr->fitted[t])) {
	    fr->resid[t] = fr->actual[t] - fr->fitted[t];
	}

	clear_model(&mod);
	pdinfo->t2 += 1;
    }

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    if (*err) {
	free_fit_resid(fr);
	fr = NULL;
    } else {
	fit_resid_set_dec_places(fr);
	strcpy(fr->depvar, pdinfo->varname[pmod->list[1]]);
    }

    free(y);
    free(llist);
	
    return fr;
}
