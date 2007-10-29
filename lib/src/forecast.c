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
#include "forecast.h"
#include "var.h"

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

static int 
allocate_basic_fit_resid_arrays (FITRESID *fr)
{
    int t;

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

    fr->resid = malloc(fr->nobs * sizeof *fr->resid);
    if (fr->resid == NULL) {
	free(fr->actual);
	free(fr->fitted);
	fr->actual = NULL;
	fr->fitted = NULL;
	return E_ALLOC;
    }

    for (t=0; t<fr->nobs; t++) {
	fr->actual[t] = fr->fitted[t] = fr->resid[t] = NADBL;
    }

    fr->sderr = NULL;

    return 0;
}

static int fit_resid_add_sderr (FITRESID *fr)
{
    int t;

    fr->sderr = malloc(fr->nobs * sizeof *fr->sderr);

    if (fr->sderr == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<fr->nobs; t++) {
	fr->sderr[t] = NADBL;
    }

    return 0;
}

/**
 * fit_resid_new:
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

static FITRESID *fit_resid_new (int n)
{
    FITRESID *fr = malloc(sizeof *fr);

    if (fr == NULL) {
	return NULL;
    }

    fr->model_ID = 0;
    fr->model_ci = 0;
    fr->model_t1 = 0;
    fr->err = 0;
    fr->t0 = 0;
    fr->t1 = 0;
    fr->t2 = 0;
    fr->nobs = 0;

    if (n > 0) {
	fr->nobs = n;
	if (allocate_basic_fit_resid_arrays(fr)) {
	    free(fr);
	    fr = NULL;
	} 
    } else {
	fr->actual = NULL;
	fr->fitted = NULL;
	fr->resid = NULL;
	fr->sderr = NULL;
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
    free(fr->resid);
    free(fr->sderr);
    free(fr);
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
    int depvar, t;
    FITRESID *fr;

    fr = fit_resid_new(pdinfo->n);
    if (fr == NULL) {
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

    depvar = gretl_model_get_depvar(pmod);

    fr->t0 = pmod->t1;
    fr->t1 = pmod->t1;
    fr->t2 = pmod->t2;

    for (t=0; t<fr->nobs; t++) {
	fr->actual[t] = Z[depvar][t];
	fr->fitted[t] = pmod->yhat[t];
	fr->resid[t] = pmod->uhat[t];
    }

    fit_resid_set_dec_places(fr);

    strcpy(fr->depvar, pdinfo->varname[depvar]);
    
    return fr;
}

/* local shortcut to get a model's list of regressors */

static int *model_xlist (MODEL *pmod)
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
    int *xlist;
    char xname[VNAMELEN], tmp[VNAMELEN];  
    const char *label;
    const char *yname;
    int i, lag;
    int ret = 0;

    xlist = model_xlist(pmod);
    if (xlist == NULL) {
	return 0;
    }

    yname = pdinfo->varname[gretl_model_get_depvar(pmod)];

    for (i=1; i<=xlist[0]; i++) {
	label = VARLABEL(pdinfo, xlist[i]);
	if ((sscanf(label, "= %15[^(](t - %d)", xname, &lag) == 2 ||
	     sscanf(label, "%15[^=]=%15[^(](-%d)", tmp, xname, &lag) == 3) &&
	    !strcmp(xname, yname)) {
	    ret = 1;
	    break;
	}
    }

    return ret;
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
    int *xlist = NULL;
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
	err = 1;
    } else {
	char xname[VNAMELEN], tmp[VNAMELEN];
	const char *label;
	const char *yname;
	int i, lag;

	yname = pdinfo->varname[gretl_model_get_depvar(pmod)];

	for (i=0; i<xlist[0]; i++) {
	    label = VARLABEL(pdinfo, xlist[i+1]);
	    if ((sscanf(label, "= %15[^(](t - %d)", xname, &lag) == 2 ||
		 sscanf(label, "%15[^=]=%15[^(](-%d)", tmp, xname, &lag) == 3) &&
		!strcmp(xname, yname)) {
		dvlags[i] = lag;
	    } else {
		dvlags[i] = 0;
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
    int *xlist = model_xlist(pmod);
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

static int 
fit_resid_init (int t0, int t1, int t2, const MODEL *pmod, 
		const DATAINFO *pdinfo, FITRESID *fr)
{
    if (t0 > t1) { 
	t0 = t1;
    }

    fr->t0 = t0;
    fr->t1 = t1;
    fr->t2 = t2;

    if (fr->t1 < 0 || fr->t2 < 0 || fr->t2 < fr->t1) {
	fr->err = E_OBS;
    }

    if (!fr->err) {
	fr->nobs = pdinfo->n;
	fr->err = allocate_basic_fit_resid_arrays(fr);
    }

    fr->model_ID = pmod->ID;
    fr->model_ci = pmod->ci;
    fr->model_t1 = pmod->t1;

    fr->pmax = PMAX_NOT_AVAILABLE;

    return fr->err;
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
    int k = pmod->ncoeff;
    int i, vi, t;
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

    for (t=fc->t1; t<=fc->t2 && !err; t++) {
	int missing = 0;
	double vyh;

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
	fc->yhat[t] = gretl_matrix_dot_product(Xs, GRETL_MOD_NONE,
					       b, GRETL_MOD_NONE,
					       NULL);

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

/* Generate forecasts from nonlinear least squares model, using
   the string specification of the regression function that
   was saved as data on the model (see nls.c).  For now we don't
   attempt to calculate forecast error variance.
*/

static int nls_fcast (Forecast *fc, const MODEL *pmod, 
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
	pdinfo->t1 = fc->t1;
	pdinfo->t2 = fc->t2;
	sprintf(formula, "$nl_y = %s", nlfunc);
	err = generate(formula, pZ, pdinfo, OPT_P, NULL);
    }

    if (!err) {
	/* transcribe values from last generated var to target */
	for (t=fc->t1; t<=fc->t2; t++) {
	    fc->yhat[t] = (*pZ)[oldv][t];
	}
	err = dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo);
    }

    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

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
    int *xlist = NULL;
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
   computing the AR portion of the forecast.
*/

static double *create_Xb_series (Forecast *fc, const MODEL *pmod,
				 const double *beta, int *xlist, 
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
    int *xlist = NULL;
    int p, q, px = 0, npsi = 0;
    int t1 = fc->t1;
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
	t1 = pmod->t2 + 1;
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
	npsi = fc->t2 - t1 + 1;
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
    for (t=t1; t<=fc->t2 && !err; t++) {
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
	    vl = arma_variance(phi, p, theta, q, psi, npsi, t - t1 + 1);
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
   (AR, CORC, HILU, PWE).  This is complicated by the fact that there
   may be a lagged dependent variable in the picture.  If there is,
   the effective AR coefficients have to be incremented, for the
   purposes of calculating forecast variance.  But I'm not sure this
   is quite right yet.
*/

static int
set_up_ar_fcast_variance (Forecast *fc, const MODEL *pmod, 
			  int pmax, int npsi, double **phi, double **psi,
			  double **errphi)
{
    int err = 0;
    
    *errphi = NULL;
    *psi = NULL;
    *phi = NULL;

    *errphi = make_phi_from_arinfo(pmod->arinfo, pmax);
    if (*errphi == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	*psi = malloc(npsi * sizeof **psi);
	if (*psi == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	*phi = malloc((pmax + 1) * sizeof **phi);
	if (*phi == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	free(*errphi);
	*errphi = NULL;
	free(*psi);
	*psi = NULL;
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

   This code is used for AR, CORC, HILU and PWE models; it
   is also used for dynamic forecasting with models that do
   not have an explicit AR error process but that have one
   or more lagged values of the dependent variable as
   regressors.  
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
    int npsi = 0;
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

        if (pmod->ci == PWE && t == pmod->t1) {
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

static int linear_fcast (Forecast *fc, const MODEL *pmod, 
			 const double **Z, const DATAINFO *pdinfo)
{
    const double *offvar = NULL;
    int yno = pmod->list[1];
    double lmax = NADBL;
    double xval;
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

    for (t=fc->t1; t<=fc->t2; t++) {
	int miss = 0;
	double yh = 0.0;

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
		yh += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    fc->yhat[t] = NADBL;
	} else if (FCAST_SPECIAL(pmod->ci)) {
	    /* special handling for LOGIT and others */
	    fc->yhat[t] = fcast_transform(yh, pmod->ci, t, offvar, lmax);
	} else {
	    fc->yhat[t] = yh;
	}
    }

    return 0;
}

static int get_forecast_method (Forecast *fc,
				MODEL *pmod, 
				const DATAINFO *pdinfo,
				gretlopt opt)
{
    int dyn_ok = 0;
    int dyn_errs_ok = 0;
    int err = 0;

    fc->dvlags = NULL;
    fc->method = FC_STATIC;

    if ((opt & OPT_D) && (opt & OPT_S)) {
	/* conflicting options: remove them */
	fputs("got conflicting options, static and dynamic\n", stderr);
	opt &= ~OPT_D;
	opt &= ~OPT_S;
	err = 1;
    }

    /* do setup for possible lags of the dependent variable,
       unless OPT_S for "static" has been given */
    if (dataset_is_time_series(pdinfo) && !(opt & OPT_S) &&
	pmod->ci != ARMA) {
	process_lagged_depvar(pmod, pdinfo, &fc->dvlags);
    }

    if (!(opt & OPT_S)) {
	if (pmod->ci == ARMA || fc->dvlags != NULL) {
	    /* dynamic forecast is possible, and not ruled out by 
	       "static" option */
	    dyn_ok = 1;
	}
	if (SIMPLE_AR_MODEL(pmod->ci) || pmod->ci == GARCH) {
	    dyn_errs_ok = 1;
	}
    }    

    if (!dyn_ok && (opt & OPT_D)) {
	/* "dynamic" option given, but can't be honored */
	fputs("requested dynamic option, but it is not applicable\n", stderr);
	opt &= ~OPT_D;
	err = 1;
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
    int nf = 0;
    int t, err = 0;

    forecast_init(&fc);

    fc.t1 = fr->t1;
    fc.t2 = fr->t2;
    fc.model_t2 = pmod->t2;

    get_forecast_method(&fc, pmod, pdinfo, opt);

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

    if (DM_errs || dyn_errs) {
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
	err = linear_fcast(&fc, pmod, Z, pdinfo);
    }

    forecast_free(&fc);

    if (dummy_AR) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
    }

    for (t=0; t<fr->nobs; t++) {
	if (t < fr->t1) {
	    if (t >= pmod->t1 && t <= pmod->t2) {
		fr->fitted[t] = pmod->yhat[t];
		fr->resid[t] = pmod->uhat[t];
	    } else {
		fr->fitted[t] = fr->resid[t] = NADBL;
	    }
	    if (fr->sderr != NULL) {
		fr->sderr[t] = NADBL;
	    }
	} else if (!na(fr->fitted[t])) {
	    nf++;
	}
	fr->actual[t] = (*pZ)[yno][t];
    }

    if (nf == 0) {
	err = E_MISSDATA;
    } else {
	if (pmod->ci == ARMA) {
	    /* asymptotic normal */
	    fr->tval = 1.96;
	} else {
	    fr->tval = tcrit95(pmod->dfd);
	}

	fit_resid_set_dec_places(fr);

	strcpy(fr->depvar, pdinfo->varname[yno]);
	fr->df = pmod->dfd;
    }

    return err;
}

static int parse_forecast_string (const char *s, const MODEL *pmod,
				  const DATAINFO *pdinfo,
				  int *t1, int *t2)
{
    char t1str[OBSLEN], t2str[OBSLEN];
    int err = 0;

    if (!strncmp(s, "fcasterr", 8)) {
        s += 9;
    }

    if (sscanf(s, "%10s %10s", t1str, t2str) == 2) {
	*t1 = dateton(t1str, pdinfo);
	*t2 = dateton(t2str, pdinfo);
    } else if (pmod != NULL && pdinfo->n - pmod->t2 - 1 > 0) {
	/* default, if it works: out-of-sample range */
	*t1 = pmod->t2 + 1;
	*t2 = pdinfo->n - 1;
    } else {
        err = E_OBS;
    } 

    return err;
}

/* public forecast-related functions follow */

/**
 * get_forecast:
 * @pmod: the model from which forecasts are wanted.
 * @t0: start of data display range.
 * @t1: start of forecast range.
 * @t2: end of forecast range.
 * @pZ: pointer to data array using which @pmod was estimated.
 * @pdinfo: dataset information.
 * @opt: if OPT_D, force a dynamic forecast; if OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where a dynamic forecast is meaningful).
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

FITRESID *get_forecast (MODEL *pmod, int t0, int t1, int t2,
			double ***pZ, DATAINFO *pdinfo,
			gretlopt opt) 
{
    FITRESID *fr;

    fr = fit_resid_new(0); 
    if (fr == NULL) {
	return NULL;
    }

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    fit_resid_init(t0, t1, t2, pmod, pdinfo, fr); 
    if (fr->err) {
	return fr;
    }

    fr->err = real_get_fcast(fr, pmod, pZ, pdinfo, opt);

    return fr;
}

static int organize_fcast_vars (const char *yhname, const char *sdname,
				int *vi, int *vj, double ***pZ,
				DATAINFO *pdinfo)
{
    int nadd = 0;
    int err = 0;

    *vi = *vj = 0;

    if (check_varname(yhname)) {
	return 1;
    }

    if (*sdname && check_varname(sdname)) {
	return 1;
    }

    *vi = varindex(pdinfo, yhname);
    if (*vi == pdinfo->v) {
	nadd++;
    }

    if (*sdname) {
	*vj = varindex(pdinfo, sdname);
	if (*vj == pdinfo->v) {
	    if (nadd > 0) {
		*vj += 1;
	    }
	    nadd++;
	}
    }

    if (nadd > 0) {
	err = dataset_add_series(nadd, pZ, pdinfo);
    }

    return err;
}

/**
 * add_forecast:
 * @str: command string, giving a starting observation, ending
 * observation, and variable name to use for the forecast values
 * (the starting and ending observations may be omitted).
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: pointer to data information struct.
 * @opt: if OPT_D, force a dynamic forecast; if OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where this distinction is meaningful).
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

int add_forecast (const char *str, MODEL *pmod, double ***pZ,
		  DATAINFO *pdinfo, gretlopt opt)
{
    int oldv = pdinfo->v;
    int t, t1, t2, vi, vj;
    char t1str[OBSLEN], t2str[OBSLEN];
    char yhname[VNAMELEN];
    char sdname[VNAMELEN];
    int nf, err = 0;

    if (pmod->ci == ARBOND) {
	return E_NOTIMP; /* FIXME */
    }

    *t1str = '\0'; *t2str = '\0';

    /* Reject in case model was estimated using repacked daily
       data: this case should be handled more elegantly */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	return E_DATA;
    }

    *yhname = *sdname = 0;

    /* the name for forecast values should either be in the 2nd or 4th
       position; it may be followed by a name for the standard errors
    */
    nf = sscanf(str, "%*s %10s %10s %15s %15s", t1str, t2str, yhname, sdname);
    if (nf < 3) {
	if (sscanf(str, "%*s" "%15s %15s", yhname, sdname) < 1) {
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

    err = organize_fcast_vars(yhname, sdname, &vi, &vj, pZ, pdinfo);
    if (err) {
	return err;
    }

    nf = 0;

    if (!err) {
	const double **Z = (const double **) *pZ;
	Forecast fc;

	forecast_init(&fc);

	strcpy(pdinfo->varname[vi], yhname);
	strcpy(VARLABEL(pdinfo, vi), _("predicted values"));

	if (vj > 0) {
	    strcpy(pdinfo->varname[vj], sdname);
	    strcpy(VARLABEL(pdinfo, vj), _("forecast standard errors"));
	}

	fc.yhat = (*pZ)[vi];
	if (vj > 0) {
	    fc.sderr = (*pZ)[vj];
	} 

	for (t=0; t<pdinfo->n; t++) {
	    fc.yhat[t] = NADBL;
	    if (fc.sderr != NULL) {
		fc.sderr[t] = NADBL;
	    }	
	}

	fc.t1 = t1;
	fc.t2 = t2;
	fc.model_t2 = pmod->t2;

	get_forecast_method(&fc, pmod, pdinfo, opt);

	if (pmod->ci == ARMA && fc.method == FC_STATIC) {
	    fc_add_eps(&fc, pdinfo->n);
	}

	/* write forecast values into the newly added variable(s) */
	if (pmod->ci == NLS) {
	    nls_fcast(&fc, pmod, pZ, pdinfo);
	} else if (SIMPLE_AR_MODEL(pmod->ci)) {
	    ar_fcast(&fc, pmod, Z, pdinfo);
	} else if (pmod->ci == ARMA) {
	    arma_fcast(&fc, pmod, Z, pdinfo);
	} else if (pmod->ci == GARCH) {
	    garch_fcast(&fc, pmod, Z, pdinfo);
	} else {
	    linear_fcast(&fc, pmod, Z, pdinfo);
	}

	forecast_free(&fc);

	for (t=0; t<pdinfo->n; t++) {
	    if (!na(fc.yhat[t])) {
		nf++;
	    }
	}
    }

    if (nf == 0) {
	dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo);
	err = E_DATA;
    }

    return err;
}

static int model_display_forecast (const char *str, MODEL *pmod, 
				   double ***pZ, DATAINFO *pdinfo, 
				   gretlopt opt, PRN *prn)
{
    FITRESID *fr;
    int t1, t2;
    int err;

    if (pmod->ci == ARBOND) {
	return E_NOTIMP;
    }

    err = parse_forecast_string(str, pmod, pdinfo, &t1, &t2);
    if (err) {
	return err;
    }

    fr = get_forecast(pmod, t1, t1, t2, pZ, pdinfo, opt);

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

static int VAR_display_forecast (const char *str, GRETL_VAR *var, 
				 double ***pZ, DATAINFO *pdinfo, 
				 gretlopt opt, PRN *prn)
{
    FITRESID *fr;
    const MODEL *pmod;
    int t1, t2;
    int i, n, err;

    pmod = gretl_VAR_get_model(var, 0);

    err = parse_forecast_string(str, pmod, pdinfo, &t1, &t2);
    if (err) {
	return err;
    }

    n = gretl_VAR_get_n_equations(var);

    for (i=0; i<n && !err; i++) {

	if (i > 0) {
	    opt |= OPT_Q;
	}

	fr = get_VAR_forecast(var, i, t1, t1, t2, (const double **) *pZ, 
			      pdinfo, opt);

	if (fr == NULL) {
	    return E_ALLOC;
	}

	err = fr->err;

	if (!err) {
	    err = text_print_forecast(fr, pZ, pdinfo, opt, prn);
	}

	free_fit_resid(fr);
    }

    return err;
}

/**
 * display_forecast:
 * @str: string giving starting and ending observations, separated
 * by a space.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: if OPT_D, force a dynamic forecast; if OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where this distinction is meaningful).
 * @prn: printing structure.
 *
 * Computes forecasts based on the last model estimated, over the 
 * range of observations given in @str.  Forecast standard errors 
 * are also computed if possible.  The results are printed to @prn, 
 * and are also plotted if %OPT_P is given.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int display_forecast (const char *str, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    GretlObjType type;
    void *ptr;
    int err;

    ptr = get_last_model(&type);  
    if (ptr == NULL) {
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN) {
	err = model_display_forecast(str, ptr, pZ, pdinfo, opt, prn);
    } else if (type == GRETL_OBJ_SYS) {
	err = E_NOTIMP;
    } else if (type == GRETL_OBJ_VAR) {
	err = VAR_display_forecast(str, ptr, pZ, pdinfo, opt, prn);
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
 * get_VAR_forecast:
 * @var: the VAR system from which forecasts are wanted.
 * @i: 0-based index for the variable to forecast, within
 * the VAR system (the dependent variable in the ith equation,
 * counting from 0).
 * @t0: start of data display range.
 * @t1: start of forecast range.
 * @t2: end of forecast range.
 * @Z: data array using which @var was estimated.
 * @pdinfo: dataset information.
 * @opt: if OPT_D, force a dynamic forecast; if OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample.
 *
 * Allocates a #FITRESID structure and fills it out with forecasts
 * based on @var, over the specified range of observations.  
 * If @t0 is less than @t1, the observations prior to @t1 will
 * not contain forecasts.
 * 
 * Returns: pointer to allocated structure, or %NULL on failure.
 * The %err member of the returned object should be checked:
 * a non-zero value indicates an error condition.
 */

FITRESID *get_VAR_forecast (GRETL_VAR *var, int i, int t0, int t1, int t2,
			    const double **Z, DATAINFO *pdinfo,
			    gretlopt opt)
{
    FITRESID *fr;
    const gretl_matrix *F;
    const MODEL *pmod = NULL;
    int nf = t2 - t1 + 1;
    int yno, m, s, t;

    if (nf <= 0) {
	return NULL;
    }

    pmod = gretl_VAR_get_model(var, i);
    if (pmod == NULL) {
	return NULL;
    }

    F = gretl_VAR_get_forecast_matrix(var, t0, t1, t2, Z, pdinfo, opt);
    if (F == NULL) {
	fprintf(stderr, "gretl_VAR_get_forecast_matrix() gave NULL\n");
	return NULL;
    }

    fr = fit_resid_new(pdinfo->n);
    if (fr == NULL) {
	return NULL;
    }
    
    if (!(opt & OPT_S)) {
	if (fit_resid_add_sderr(fr)) {
	    free_fit_resid(fr);
	    return NULL;
	}
    }

    fr->model_ci = var->ci;

    fr->t0 = t0;
    fr->t1 = t1;
    fr->t2 = t2;

    yno = pmod->list[1];

    strcpy(fr->depvar, gretl_model_get_depvar_name(pmod, pdinfo));

    m = var->neqns;

    nf = 0;

    for (t=fr->t0, s=0; t<=fr->t2; t++, s++) {
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
	fr->err = E_MISSDATA;
    } else {
	if (var->ci == VECM) {
	    /* asymptotic normal */
	    fr->df = var->T;
	    fr->tval = 1.96;
	} else {
	    fr->df = pmod->dfd;
	    fr->tval = tcrit95(fr->df);
	}
	fit_resid_set_dec_places(fr);
	strcpy(fr->depvar, pdinfo->varname[yno]);
    }

    return fr;
}

/**
 * forecast_options_for_model:
 * @pmod: the model from which forecasts are wanted.
 * @Z: data array.
 * @pdinfo: dataset information.
 * @dyn_ok: location to receive 1 if the "dynamic" option is
 * applicable, 0 otherwise.
 * @add_obs_ok: location to receive 1 if it looks as if we can
 * extend the forecast by adding blank observations, 0 otherwise.
 * @dt2max: location to receive the last observation that can
 * be supported for a dynamic forecast.
 * @st2max: location to receive the last observation that can
 * be supported for a static forecast.
 *
 * Examines @pmod and determines which forecasting options are
 * applicable.
 */

void forecast_options_for_model (MODEL *pmod, const double **Z,
				 const DATAINFO *pdinfo,
				 int *dyn_ok, int *add_obs_ok,
				 int *dt2max, int *st2max)
{
    int *dvlags = NULL;
    int exo = 1;

    *dyn_ok = 0;
    *add_obs_ok = 0;
    *dt2max = pmod->t2;
    *st2max = pmod->t2;

    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	*dyn_ok = 1;
    } else if (AR_MODEL(pmod->ci)) {
	*dyn_ok = 1;
    } else if (dataset_is_time_series(pdinfo) &&
	       has_depvar_lags(pmod, pdinfo)) {
	*dyn_ok = 1;
    }

    if (*dyn_ok) {
	int err = process_lagged_depvar(pmod, pdinfo, &dvlags);

	if (!err) {
	    exo = has_real_exog_regressors(pmod, dvlags, Z, pdinfo);
	}
	if (!exo) {
	    *add_obs_ok = 1;
	    *dt2max = pdinfo->n - 1;
	}
    } 

    if (exo) {
	int *xlist = model_xlist(pmod);

	if (xlist != NULL) {
	    *dt2max = fcast_get_t2max(xlist, dvlags, pmod, Z, pdinfo, FC_DYNAMIC);
	    *st2max = fcast_get_t2max(xlist, dvlags, pmod, Z, pdinfo, FC_STATIC);
	}
    }

    if (dvlags != NULL) {
	free(dvlags);
    }
}

/* recursive one-step ahead forecasts and forecast errors, for models
   estimated via OLS */

FITRESID * 
rolling_OLS_one_step_fcast (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			    int t0, int t1, int t2)
{
    FITRESID *fr;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    double xit, yf;
    int nf = t2 - t1 + 1;
    MODEL mod;
    int i, s, t;
    int err = 0;

    fr = fit_resid_new(0); 
    if (fr == NULL) {
	return NULL;
    }

    if (pmod->ci != OLS) {
	fr->err = E_OLSONLY;
	return fr;
    }

    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    if (t1 - t0 < pmod->ncoeff || t2 < t1) {
	fr->err = E_DATA;
	return fr;
    }

    fit_resid_init(t1, t1, t2, pmod, pdinfo, fr); 
    if (fr->err) {
	return fr;
    }

    fr->method = FC_ONESTEP;

    pdinfo->t1 = t0;
    pdinfo->t2 = t1 - 1;

    for (t=0; t<pdinfo->n; t++) {
	fr->actual[t] = (*pZ)[pmod->list[1]][t];
    }

    for (s=0; s<nf; s++) {
	mod = lsq(pmod->list, pZ, pdinfo, OLS, OPT_A);
	if (mod.errcode) {
	    err = mod.errcode;
	    clear_model(&mod);
	    break;
	}

	t = pdinfo->t2 + 1;
	yf = 0.0;
	for (i=0; i<mod.ncoeff; i++) {
	    xit = (*pZ)[mod.list[i+2]][t];
	    if (na(xit)) {
		yf = NADBL;
		break;
	    } else {
		yf += mod.coeff[i] * xit;
	    }
	}

	fr->fitted[t] = yf;
	
	clear_model(&mod);
	pdinfo->t2 += 1;
    }

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    if (!fr->err) {
	fit_resid_set_dec_places(fr);
	strcpy(fr->depvar, pdinfo->varname[pmod->list[1]]);
    }
	
    return fr;
}
