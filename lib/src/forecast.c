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
#include "uservar.h"
#include "gretl_midas.h"
#include "forecast.h"
#include "var.h"
#include "system.h"
#include "libset.h"

#define AR_DEBUG 0

/**
 * SECTION:forecast
 * @short_description: Forecasting based on models
 * @title: Forecasting
 * @include: gretl/libgretl.h, gretl/forecast.h
 *
 * Write-up needed here.
 *
 */

/* estimators where a simple X*b does _not_ give the
   predicted value of the dependent variable */

#define FCAST_SPECIAL(c) (c == LOGIT || \
                          c == LOGISTIC || \
			  c == MIDASREG || \
			  c == NEGBIN || \
                          c == NLS || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == TOBIT)

typedef struct Forecast_ Forecast;

struct Forecast_ {
    int method;       /* static, dynamic or auto */
    double *yhat;     /* array of forecast values */
    double *sderr;    /* array of forecast standard errors */
    double *eps;      /* array of estimated forecast errors */
    int *dvlags;      /* info on lags of dependent variable */
    int t1;           /* start of forecast range */
    int t2;           /* end of forecast range */
    int model_t2;     /* end of period over which model was estimated */
};

static gretl_matrix *fcast_matrix;
static gretl_matrix *fcerr_matrix;

static void set_fcast_matrices (gretl_matrix *fc,
				gretl_matrix *se)
{
    gretl_matrix_free(fcast_matrix);
    gretl_matrix_free(fcerr_matrix);

    fcast_matrix = fc;
    fcerr_matrix = se;
}

/* Create an empty, dummy AR info structure for use with models
   that don't have an explicit AR error process, but that do
   have a lagged dependent variable that in effect produces
   an AR error, for forecasting purposes. Or add a single-
   element AR info structure for the AR1 case.
*/

static int dummy_ar_info_init (MODEL *pmod)
{
    int err = 0;

    pmod->arinfo = malloc(sizeof *pmod->arinfo);

    if (pmod->arinfo == NULL) {
	err = E_ALLOC;
    } else {
	pmod->arinfo->rho = NULL;
	pmod->arinfo->sderr = NULL;
    }

    if (!err && pmod->ci == AR1) {
	pmod->arinfo->arlist = gretl_list_new(1);
	pmod->arinfo->rho = malloc(sizeof(double));
	pmod->arinfo->arlist[1] = 1;
	pmod->arinfo->rho[0] = gretl_model_get_double(pmod, "rho_gls");
    } else if (!err) {
	pmod->arinfo->arlist = gretl_null_list();
    }

    if (err) {
	free(pmod->arinfo);
	pmod->arinfo = NULL;
    }

    return err;
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
    f->std = 0;
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
					  const DATASET *dset,
					  int t1, int t2, int pre_n,
					  int *err)
{
    FITRESID *fr;

    if (t1 < 0 || t2 < 0 || t2 < t1) {
	*err = E_OBS;
	return NULL;
    }

    if (dset != NULL) {
        fr = fit_resid_new_with_length(dset->n, 0);
    } else {
        fr = fit_resid_new_with_length(t2 - t1 + 1, 0);
    }

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
					   const DATASET *dset,
					   int t1, int t2, int pre_n,
					   int *err)
{
    FITRESID *fr;

    if (t1 < 0 || t2 < 0 || t2 < t1) {
	*err = E_OBS;
	return NULL;
    }

    fr = fit_resid_new_with_length(dset->n, 1);

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

static int special_depvar (const MODEL *pmod)
{
    int ret = 0;

    if (pmod->ci == INTREG) {
	ret = 1;
    } else if (pmod->ci == PANEL) {
	if (pmod->opt & OPT_B) {
	    /* panel "between" model */
	    ret = 1;
	}
    }

    return ret;
}

/**
 * get_fit_resid:
 * @pmod: the model for which actual and fitted values
 * are wanted.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Allocates a #FITRESID structure and fills it out with
 * the actual and predicted values of the dependent variable
 * in @pmod.
 *
 * Returns: pointer to allocated structure, or %NULL on failure.
 */

FITRESID *get_fit_resid (const MODEL *pmod, const DATASET *dset,
			 int *err)
{
    double *uhat = pmod->uhat;
    double *yhat = pmod->yhat;
    FITRESID *fr;
    int t, dv = -1;

    if (!special_depvar(pmod)) {
	dv = gretl_model_get_depvar(pmod);
	if (dv < 0 || dv >= dset->v) {
	    *err = E_DATA;
	    return NULL;
	}
    }

    if (uhat == NULL || yhat == NULL) {
	gretl_matrix *uh = gretl_model_get_data(pmod, "uhat");
	gretl_matrix *yh = gretl_model_get_data(pmod, "yhat");

	if (uh == NULL || yh == NULL) {
	    *err = E_DATA;
	    return NULL;
	} else {
	    uhat = uh->val;
	    yhat = yh->val;
	}
    }

    fr = fit_resid_new_for_model(pmod, dset, pmod->t1, pmod->t2,
				 0, err);
    if (*err) {
	return NULL;
    }

    if (LIMDEP(pmod->ci)) {
	fr->sigma = NADBL;
    } else if (pmod->ci == GARCH && (pmod->opt & OPT_Z)) {
	fr->sigma = 1.0;
	fr->std = 1;
    } else {
	fr->sigma = gretl_model_get_double(pmod, "sigma_orig");
	if (na(fr->sigma)) {
	    fr->sigma = pmod->sigma;
	}
    }

    for (t=0; t<fr->nobs; t++) {
	if (dv < 0) {
	    if (na(yhat[t]) || na(uhat[t])) {
		fr->actual[t] = NADBL;
	    } else {
		fr->actual[t] = yhat[t] + uhat[t];
	    }
	} else {
	    fr->actual[t] = dset->Z[dv][t];
	}
	fr->fitted[t] = yhat[t];
	fr->resid[t] = uhat[t];
    }

    if (dv < 0) {
	fr->pmax = PMAX_NOT_AVAILABLE;
    } else {
	fit_resid_set_dec_places(fr);
    }

    if (dv < 0) {
	const char *dvname = pmod->depvar;

	if (dvname != NULL) {
	    strcpy(fr->depvar, dvname);
	} else {
	    strcpy(fr->depvar, "implicit y");
	}
    } else {
	strcpy(fr->depvar, dset->varname[dv]);
    }

    return fr;
}

/* local shortcut to get a model's list of regressors */

static const int *model_xlist (MODEL *pmod)
{
    int *xlist = gretl_model_get_list(pmod, "xlist");

    if (xlist == NULL) {
	xlist = gretl_model_get_x_list(pmod);
	if (xlist != NULL) {
	    gretl_model_set_list_as_data(pmod, "xlist", xlist);
	}
    }

    return xlist;
}

static inline int xy_equal(double x, double y)
{
    return x == y || (na(x) && na(y));
}

static int detect_lag (int vy, int vx, int t1, int t2,
                       int pmax, const DATASET *dset)
{
    const double *y = dset->Z[vy];
    const double *x = dset->Z[vx];
    int p, t, ok;
    int ret = 0;

    for (p=1; p<=pmax; p++) {
        ok = 1;
        for (t=t1+p; t<=t2; t++) {
            if (!xy_equal(x[t], y[t-p])) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            ret = p;
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
   to the lag order.  However, in case the list of independent vars
   contains no lagged dependent var, *depvar_lags is set to NULL.
   Returns 1 on error, 0 otherwise.
*/

static int *process_lagged_depvar (MODEL *pmod, const DATASET *dset)
{
    const int *xlist = NULL;
    int *dvlags = NULL;
    int nl = 0;

    xlist = model_xlist(pmod);
    if (xlist != NULL) {
        dvlags = malloc(xlist[0] * sizeof *dvlags);
    }

    if (dvlags != NULL) {
        int yv = gretl_model_get_depvar(pmod);
	int i, vi, lag;

	for (i=1; i<=xlist[0]; i++) {
            lag = 0;
            vi = xlist[i];
            if (series_get_transform(dset, vi) == LAGS) {
                if (series_get_parent_id(dset, vi) == yv) {
                    lag = series_get_lag(dset, vi);
                }
            } else if (vi > 0) {
                int pmax = dset->pd < 4 ? 4 : dset->pd;

                lag = detect_lag(yv, vi, dset->t1, dset->t2,
                                 pmax, dset);
            }
            dvlags[i-1] = lag;
            if (lag > 0) {
                nl++;
            }
        }
    }

    if (nl == 0) {
        /* didn't find any */
        free(dvlags);
        dvlags = NULL;
    }

    return dvlags;
}

/* Get a value for a lag of the dependent variable.  If method is
   dynamic we prefer lagged prediction to lagged actual.  If method is
   static, we never want the lagged prediction, only the actual.  If
   method is "auto", which value we prefer depends on whether we're in
   or out of sample (actual within, lagged prediction without).
*/

static double fcast_get_ldv (Forecast *fc, int i, int t, int lag,
			     const DATASET *dset)
{
    int p = t - lag;
    double ldv;

    /* initialize to actual lagged value, if available */
    if (p < 0) {
	ldv = NADBL;
    } else {
	ldv = dset->Z[i][p];
    }

#if AR_DEBUG
    fprintf(stderr, "fcast_get_ldv: i=%d, t=%d, lag=%d; "
	    "initial ldv = Z[%d][%d] = %g\n", i, t, lag, i, p, ldv);
#endif

    if (fc->method != FC_STATIC && p >= 0) {
	double yhat = fc->yhat[p];
	int fc_ok = !na(yhat);

	if (fc->method == FC_DYNAMIC) {
	    if (fc_ok) {
		ldv = yhat;
	    }
	} else if (fc->method == FC_AUTO) {
	    if (t > fc->model_t2 + lag || na(ldv)) {
		if (fc_ok) {
		    ldv = yhat;
		}
	    }
	}
#if AR_DEBUG
	fprintf(stderr, "fcast_get_ldv (non-static): using %g\n", ldv);
#endif
    }

    return ldv;
}

/* Get forecasts, plus standard errors for same, for models without
   autoregressive errors and without "special requirements"
   (e.g. nonlinearity).

   The forecast standard errors include both uncertainty over the
   error process and parameter uncertainty (Davidson and MacKinnon
   method) -- unless @opt contains OPT_M, in which case the standard
   errors pertain to prediction of mean Y.
*/

static int
static_fcast_with_errs (Forecast *fc, MODEL *pmod,
			const DATASET *dset,
			gretlopt opt)
{
    gretl_matrix *V = NULL;
    gretl_vector *Xt = NULL;
    gretl_vector *b = NULL;
    double s2 = pmod->sigma * pmod->sigma;
    double vyh, xval;
    int k = pmod->ncoeff;
    int same_data = 0;
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

    Xt = gretl_vector_alloc(k);
    if (Xt == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    same_data = same_dataset(pmod, dset);

    for (t=fc->t1; t<=fc->t2 && !err; t++) {
	int missing = 0;

	/* skip if we can't compute forecast */
	if (same_data && t >= pmod->t1 && t <= pmod->t2) {
	    missing = na(pmod->yhat[t]);
	}

	/* populate Xt vector for observation */
	for (i=0; i<k && !missing; i++) {
	    vi = pmod->list[i + 2];
	    xval = dset->Z[vi][t];
	    if (na(xval)) {
		fc->sderr[t] = fc->yhat[t] = NADBL;
		missing = 1;
	    } else {
		gretl_vector_set(Xt, i, xval);
	    }
	}

	if (missing) {
	    fc->sderr[t] = fc->yhat[t] = NADBL;
	    continue;
	}

	/* forecast value */
	fc->yhat[t] = gretl_vector_dot_product(Xt, b, NULL);

	/* forecast variance */
	vyh = gretl_scalar_qform(Xt, V, &err);
	if (na(vyh)) {
	    err = 1;
	} else {
	    if (!(opt & OPT_M)) {
		vyh += s2;
	    }
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
    gretl_vector_free(Xt);
    gretl_vector_free(b);

    return err;
}

#define individual_effects_model(m) (m != NULL && m->ci == PANEL && \
				     (m->opt & (OPT_F | OPT_U)))

static void special_set_fcast_matrix (gretl_matrix *fc)
{
    size_t sz;

    /* copy forecast from second column to first */
    sz = fc->rows * sizeof(double);
    memcpy(fc->val, fc->val + fc->rows, sz);

    /* and release some storage */
    gretl_matrix_realloc(fc, fc->rows, 1);

    set_fcast_matrices(fc, NULL);
}

static void special_print_fc_stats (gretl_matrix *fc,
				    PRN *prn)
{
    const double *y = fc->val;
    const double *f = y + fc->rows;
    gretl_matrix *fcs;
    int n_used;
    int err = 0;

    fcs = forecast_stats(y, f, 0, fc->rows - 1, &n_used,
			 OPT_O, &err);
    if (fcs != NULL) {
	print_fcast_stats_matrix(fcs, n_used, OPT_NONE, prn);
	gretl_matrix_free(fcs);
    }
}

static void transcribe_to_matrix (gretl_matrix *fc,
				  double yt, double fct,
				  char **S,
				  const DATASET *dset,
				  int t, int j)
{
    if (S != NULL) {
	char obsstr[OBSLEN];

	ntolabel(obsstr, t, dset);
	S[j] = gretl_strdup(obsstr);
    }
    if (fc->cols == 2) {
	gretl_matrix_set(fc, j, 0, yt);
	gretl_matrix_set(fc, j, 1, fct);
    } else {
	gretl_matrix_set(fc, j, 0, fct);
    }
}

/* specific to "dpanel": handles both within-sample
   and out of sample forecasts
*/

static int real_dpanel_fcast (double *yhat,
			      gretl_matrix *fc,
			      char **rlabels,
			      MODEL *pmod,
			      const char *mask,
			      const DATASET *dset,
			      gretlopt opt)
{
    const double *y, *xi, *b;
    double *yi = NULL;
    int *xlist, *ylags = NULL;
    int yno, free_ylags = 0;
    int fc_static = (opt & OPT_S);
    int dpdstyle = (pmod->opt & OPT_X);
    int t1min, ntdum, ifc, tdt = 0;
    int i, vi, s, t, k, u;
    int T, ti, row;

    ylags = gretl_model_get_data(pmod, "ylags");
    if (ylags == NULL) {
	ylags = gretl_consecutive_list_new(1, pmod->list[1]);
	free_ylags = 1;
    }

    /* per unit time series */
    T = dset->pd;
    yi = malloc(T * sizeof *yi);

    if (ylags == NULL || yi == NULL) {
	return E_ALLOC;
    }

#if 0
    fprintf(stderr, "real_dpanel_fcast: yhat=%p, fc=%p, mask=%p\n",
	    (void *) yhat, (void *) fc, (void *) mask);
    fprintf(stderr, " dset: n=%d, pd=%d, t1=%d, t2=%d\n",
	    dset->n, dset->pd, dset->t1, dset->t2);
#endif

    if (mask == NULL && !(opt & OPT_D)) {
	/* within-sample: default to static */
	fc_static = 1;
    }

    /* first usable time-series observation, zero based */
    t1min = gretl_model_get_int(pmod, "t1min");
    if (t1min > 0) {
	t1min--;
    } else {
	/* fallback? problematic if we have time dummies */
	t1min = 1 + ylags[ylags[0]];
    }

    /* is a constant included? */
    ifc = gretl_model_get_int(pmod, "ifc");

    /* do we have time dummies? */
    ntdum = gretl_model_get_int(pmod, "ntdum");
    if (ntdum > 0) {
	/* period to which first time dummy refers */
	tdt = t1min + ifc;
    }

#if 0
    fprintf(stderr, "ntdum=%d, t1min=%d, ifc=%d, tdt=%d\n",
	    ntdum, t1min, ifc, tdt);
#endif

    yno = gretl_model_get_depvar(pmod);
    xlist = gretl_model_get_x_list(pmod);
    y = dset->Z[yno];
    b = pmod->coeff;
    t = u = row = 0;

    for (s=0; s<=dset->t2; s++) {
	double fct = NADBL;
	int missing = 0;
	int fc_skip = 0;
	int j = 0;

	if (t == 0) {
	    /* fill dependent var vector for unit */
	    for (ti=0; ti<T; ti++) {
		yi[ti] = y[s+ti];
	    }
	}

	if (s < dset->t1) {
	    /* we haven't got to the requested sample yet */
	    goto next_t;
	}

	if (mask != NULL && mask[s] == 1) {
	    /* we're doing out-of-sample forecast and this is
	       an in-sample observation
	    */
	    fc_skip = 1;
	    goto transcribe;
	}

	if (t < t1min) {
	    /* required lags not available */
	    goto transcribe;
	}

	/* get lagged level for integration */
	if (na(yi[t-1])) {
	    missing = 1;
	} else {
	    fct = yi[t-1];
	}

	/* handle autoregressive terms */
	for (i=1; i<=ylags[0] && !missing; i++) {
	    k = ylags[i];
	    if (t - k - 1 < 0 || na(yi[t-k]) || na(yi[t-k-1])) {
		missing = 1;
	    } else {
		fct += b[j++] * (yi[t-k] - yi[t-k-1]);
	    }
	}

	/* then regular regressors, differenced */
	for (i=1; i<=xlist[0] && !missing; i++) {
	    vi = xlist[i];
	    if (vi == 0) {
		fct += b[j++];
	    } else {
		xi = dset->Z[vi];
		if (t - 1 < 0 || na(xi[s]) || na(xi[s-1])) {
		    missing = 1;
		} else {
		    fct += b[j++] * (xi[s] - xi[s-1]);
		}
	    }
	}

	/* then a time effect, if applicable */
	if (!missing && ntdum > 0 && t >= tdt) {
	    fct += b[j + t - tdt];
	    if (!dpdstyle) {
		/* time dummies differenced */
		if (t - tdt > 0) {
		    fct -= b[j + t - tdt - 1];
		}
	    }
	}

	if (missing) {
	    fct = NADBL;
	} else if (!fc_static) {
	    /* dynamic: replace observed with forecast */
	    yi[t] = fct;
	}

    transcribe:

	/* transcribe result? */
	if (yhat != NULL) {
	    yhat[s] = fct;
	} else if (fc_skip) {
	    ; /* skip */
	} else if (row >= fc->rows) {
	    fprintf(stderr, "dpanel_fcast out of bounds! s=%d, row=%d\n", s, row);
	} else {
	    transcribe_to_matrix(fc, y[s], fct, rlabels, dset, s, row++);
	}

    next_t:

	if (t == dset->pd - 1) {
	    t = 0;
	    u++;
	} else {
	    t++;
	}
    }

    free(xlist);
    free(yi);
    if (free_ylags) {
	free(ylags);
    }

    return 0;
}

/* for pooled OLS, fixed effects and random effects:
   handles both within-sample and out of sample
*/

static int real_panel_fcast (double *yhat,
			     gretl_matrix *fc,
			     char **rlabels,
			     MODEL *pmod,
			     const char *mask,
			     const DATASET *dset,
			     gretlopt opt)
{
    const double *y, *b;
    gretl_vector *a_vec = NULL;
    int i, vi, s, t;
    int T, row, err = 0;

#if 0
    fprintf(stderr, "real_panel_fcast: yhat=%p, fc=%p, mask=%p\n",
	    (void *) yhat, (void *) fc, (void *) mask);
    fprintf(stderr, " dset: n=%d, pd=%d, t1=%d, t2=%d\n",
	    dset->n, dset->pd, dset->t1, dset->t2);
#endif

    if (individual_effects_model(pmod)) {
	a_vec = gretl_model_ahat_vec(pmod, &err);
	if (err) {
	    fprintf(stderr, "real_panel_fcast: ahat vector is missing\n");
	    return err;
	}
    }

    y = dset->Z[pmod->list[1]];
    b = pmod->coeff;
    s = 0; /* panel unit index */
    T = dset->pd;
    row = 0;

    for (t=0; t<=dset->t2; t++) {
	double xit, fct = NADBL;
	int missing = 0;
	int fc_skip = 0;
	int j = 0;

	if (t > 0 && t % T == 0) {
	    /* advance unit index */
	    s++;
	}

	if (t < dset->t1) {
	    continue;
	}

	if (mask != NULL && mask[t] == 1) {
	    /* we're doing out-of-sample forecast and this is
	       an in-sample observation
	    */
	    fc_skip = 1;
	    goto transcribe;
	}

	if (a_vec != NULL) {
	    /* individual effect */
	    if (na(a_vec->val[s])) {
		missing = 1;
	    } else {
		fct = a_vec->val[s];
	    }
	} else {
	    fct = 0.0;
	}

	/* regular regressors */
	for (i=2, j=0; i<=pmod->list[0] && !missing; i++, j++) {
	    vi = pmod->list[i];
	    if (vi == 0) {
		if (!(pmod->opt & OPT_F)) {
		    fct += b[j];
		}
	    } else {
		xit = dset->Z[vi][t];
		if (na(xit)) {
		    missing = 1;
		} else {
		    fct += b[j] * xit;
		}
	    }
	}

	/* FIXME handling of time dummies? */

	if (missing) {
	    fct = NADBL;
	}

    transcribe:

	/* transcribe result? */
	if (yhat != NULL) {
	    yhat[t] = fct;
	} else if (fc_skip) {
	    ; /* skip */
	} else if (row >= fc->rows) {
	    fprintf(stderr, "panel_fcast out of bounds! t=%d, row=%d\n", t, row);
	} else {
	    transcribe_to_matrix(fc, y[t], fct, rlabels, dset, t, row++);
	}
    }

    gretl_matrix_free(a_vec);

    return 0;
}

static int panel_fcast (Forecast *fc, MODEL *pmod,
			DATASET *dset, gretlopt opt)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int err;

    dset->t1 = fc->t1;
    dset->t2 = fc->t2;
    if (pmod->ci == DPANEL) {
	err = real_dpanel_fcast(fc->yhat, NULL, NULL, pmod,
				NULL, dset, opt);
    } else {
	err = real_panel_fcast(fc->yhat, NULL, NULL, pmod,
			       NULL, dset, opt);
    }
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

/* Out-of-sample (in the time dimension) forecast routine
   for panel data, supporting pooled OLS, fixed effects,
   random effects and dpanel. This is the top-level
   driver.
*/

static int panel_os_special (MODEL *pmod, DATASET *dset,
			     const char *vname,
			     gretlopt opt,
			     PRN *prn)
{
    DATASET *fset = fetch_full_dataset();
    int named_save = *vname != '\0';
    gretl_matrix *fc = NULL;
    double *yhat = NULL;
    char *mask = dset->submask;
    char **Sr = NULL;
    int vi, sn;
    int err = 0;

    sn = get_dataset_submask_size(dset);
    if (sn != fset->n) {
	fprintf(stderr, "fullset->n = %d but submask size = %d\n",
		fset->n, sn);
	return E_DATA;
    }

    if (named_save) {
	yhat = malloc(fset->n * sizeof *yhat);
	if (yhat == NULL) {
	    return E_ALLOC;
	}
    } else {
	int nfc = fset->n - dset->n;
	int cols = opt & OPT_Q ? 1 : 2;

	fc = gretl_matrix_alloc(nfc, cols);
	if (fc == NULL) {
	    return E_ALLOC;
	}
	Sr = strings_array_new(nfc);
    }

    if (pmod->ci == DPANEL) {
	err = real_dpanel_fcast(yhat, fc, Sr, pmod, mask, fset, opt);
    } else {
	err = real_panel_fcast(yhat, fc, Sr, pmod, mask, fset, opt);
    }

    if (!err && yhat != NULL) {
	err = dataset_add_NA_series(dset, 1);
	if (!err) {
            vi = dset->v - 1;
	    err = rename_series(dset, vi, vname, OPT_NONE);
	    series_set_label(dset, vi, _("predicted values"));
	}
	if (!err) {
	    /* keep number of series in sync */
	    err = dataset_add_allocated_series(fset, yhat);
	}
	yhat = NULL;
    } else if (!err) {
	gretl_matrix_set_rownames(fc, Sr);
	if (opt & OPT_Q) {
	    /* @fc should be suitable to purpose */
	    set_fcast_matrices(fc, NULL);
	} else {
	    if (!(opt & OPT_T)) {
		/* not --stats-only */
		char **Sc = strings_array_new(2);
		int yno = gretl_model_get_depvar(pmod);

		Sc[0] = gretl_strdup(fset->varname[yno]);
		Sc[1] = gretl_strdup(_("prediction"));
		gretl_matrix_set_colnames(fc, Sc);
		gretl_matrix_print_to_prn(fc, NULL, prn);
	    }
	    if (!(opt & OPT_N)) {
		/* not --no-stats */
		special_print_fc_stats(fc, prn);
	    }
	    /* @fc needs some adjustment */
	    special_set_fcast_matrix(fc);
	}
	fc = NULL;
    }

    gretl_matrix_free(fc);
    free(yhat);

    return err;
}

#define NLS_DEBUG 0

/* Generate forecasts from nonlinear least squares model, using the
   string specification of the regression function that was saved as
   data on the model (see nls.c).  If the NLS formula is dynamic and
   the user has not requested a static forecast, we do an
   autoregressive genr out of sample.  If we're doing a static
   forecast we add a simple-minded forecast error, namely the standard
   error of the NLS regression.
*/

static int nls_fcast (Forecast *fc, const MODEL *pmod,
		      DATASET *dset)
{
    int oldt1 = dset->t1;
    int oldt2 = dset->t2;
    int origv = dset->v;
    int yno = 0;
    const char *nlfunc;
    double *y = NULL;
    char formula[MAXLINE];
    int t, err = 0;

    nlfunc = gretl_model_get_data(pmod, "nl_regfunc");
    if (nlfunc == NULL) {
	err = E_DATA;
    }

#if NLS_DEBUG
    fprintf(stderr, "nls_fcast: method=%d, nlfunc='%s'\n",
	    fc->method, nlfunc);
#endif

    if (fc->method == FC_AUTO) {
	yno = pmod->list[1];
	y = copyvec(dset->Z[yno], dset->n);
	if (y == NULL) {
	    err = E_ALLOC;
	}
    }

    /* The computation of $nl_y below may depend on computing
       some auxiliary quantities used within the nls block, so
       try that first below: nl_model_run_aux_genrs().
    */

    if (!err) {
	/* upper limit of static forecast */
	int t2 = (fc->method == FC_STATIC)? fc->t2 : pmod->t2;

	if (t2 >= fc->t1) {
	    /* non-null static range */
	    int fcv;

	    dset->t1 = fc->t1;
	    dset->t2 = t2;
	    err = nl_model_run_aux_genrs(pmod, dset);
	    fcv = dset->v;
#if NLS_DEBUG
	    fprintf(stderr, " static range %d to %d, fcv=%d\n",
		    dset->t1, dset->t2, fcv);
#endif
	    if (!err) {
		sprintf(formula, "$nl_y = %s", nlfunc);
		err = generate(formula, dset, GRETL_TYPE_SERIES,
			       OPT_P, NULL);
	    }
	    if (!err) {
		for (t=dset->t1; t<=dset->t2; t++) {
		    fc->yhat[t] = dset->Z[fcv][t];
		}
	    }
	}

	if (!err && fc->method == FC_AUTO && fc->t2 > pmod->t2) {
	    /* dynamic forecast out of sample: in this context
	       pmod->depvar is the expression to generate
	       the series in question under its own name
	    */
	    dset->t1 = pmod->t2 + 1;
	    dset->t2 = fc->t2;
#if NLS_DEBUG
	    fprintf(stderr, " dynamic range %d to %d\n", dset->t1, dset->t2);
#endif
	    err = nl_model_run_aux_genrs(pmod, dset);
	    if (!err) {
		strcpy(formula, pmod->depvar);
		err = generate(formula, dset, GRETL_TYPE_SERIES,
			       OPT_P, NULL);
	    }
	    if (!err) {
		for (t=dset->t1; t<=dset->t2; t++) {
		    fc->yhat[t] = dset->Z[yno][t];
		}
	    }
	}
    }

    /* restore dataset state */
    if (dset->v > origv) {
	err = dataset_drop_last_variables(dset, dset->v - origv);
    }
    dset->t1 = oldt1;
    dset->t2 = oldt2;

    if (y != NULL) {
	/* restore original dependent variable */
	memcpy(dset->Z[yno], y, dset->n * sizeof *y);
	free(y);
    }

    if (!err && fc->method == FC_STATIC && fc->sderr != NULL) {
#if 0   /* Not much tested: FIXME make this an option? */
	err = nls_boot_calc(pmod, dset, fc->t1, fc->t2, fc->sderr);
	if (!err) {
	    double et, s2 = pmod->sigma * pmod->sigma;

	    for (t=fc->t1; t<=fc->t2; t++) {
		if (!na(fc->sderr[t])) {
		    et = fc->sderr[t];
		    fc->sderr[t] = sqrt(et * et + s2);
		}
	    }
	}
#else
	/* kinda simple-minded */
	for (t=fc->t1; t<=fc->t2; t++) {
	    fc->sderr[t] = pmod->sigma;
	}
#endif
    }

    return err;
}

#define MIDAS_DEBUG 0

static int midas_fcast (Forecast *fc, const MODEL *pmod,
			DATASET *dset)
{
    int oldt1 = dset->t1;
    int oldt2 = dset->t2;
    int origv = dset->v;
    char formula[MAXLINE];
    char *mdsfunc = NULL;
    double *y = NULL;
    int yno = pmod->list[1];
    int t, err;

    err = midas_forecast_setup(pmod, dset, fc->method, &mdsfunc);

#if MIDAS_DEBUG
    fprintf(stderr, "midas_fcast: method %d, err = %d\n", fc->method, err);
#endif

    if (!err && fc->method == FC_AUTO) {
	y = copyvec(dset->Z[yno], dset->n);
	if (y == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* upper limit of static forecast */
	int t2 = (fc->method == FC_STATIC)? fc->t2 : pmod->t2;

	if (t2 >= fc->t1) {
	    /* non-null static range */
	    int fcv = dset->v;

	    dset->t1 = fc->t1;
	    dset->t2 = t2;
	    sprintf(formula, "$midas_y=%s", mdsfunc);
	    err = generate(formula, dset, GRETL_TYPE_SERIES,
			   OPT_P, NULL);
#if MIDAS_DEBUG
	    fprintf(stderr, " static range %d to %d, fcv=%d, err=%d\n",
		    dset->t1, dset->t2, fcv, err);
#endif
	    if (!err) {
		for (t=dset->t1; t<=dset->t2; t++) {
		    fc->yhat[t] = dset->Z[fcv][t];
		}
	    }
	}

	if (!err && fc->method == FC_AUTO && fc->t2 > pmod->t2) {
	    /* dynamic forecast out of sample */
	    dset->t1 = pmod->t2 + 1;
	    dset->t2 = fc->t2;
	    sprintf(formula, "%s=%s", dset->varname[yno], mdsfunc);
#if MIDAS_DEBUG
	    fprintf(stderr, " dynamic range %d to %d\n %s\n", dset->t1,
		    dset->t2, formula);
#endif
	    err = generate(formula, dset, GRETL_TYPE_SERIES,
			   OPT_P, NULL);
	    if (!err) {
		for (t=dset->t1; t<=dset->t2; t++) {
		    fc->yhat[t] = dset->Z[yno][t];
		}
	    }
	}
    }

    /* restore dataset state */
    if (dset->v > origv) {
	err = dataset_drop_last_variables(dset, dset->v - origv);
    }
    dset->t1 = oldt1;
    dset->t2 = oldt2;

    if (y != NULL) {
	/* restore original dependent variable */
	memcpy(dset->Z[yno], y, dset->n * sizeof *y);
	free(y);
    }

    if (!err && fc->method == FC_STATIC && fc->sderr != NULL) {
	for (t=fc->t1; t<=fc->t2; t++) {
	    fc->sderr[t] = pmod->sigma;
	}
    }

    /* clean up any uservars the MIDAS mechanism created */
    destroy_private_uvars();

    free(mdsfunc);

    return err;
}

#if AR_DEBUG
# include <stdarg.h>
static void my_dprintf (const char *format, ...)
{
   va_list args;

   va_start(args, format);
   vfprintf(stderr, format, args);
   va_end(args);

   return;
}
# define DPRINTF(x) my_dprintf x
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

    lmax = MAX(p, q);

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
			const DATASET *dset)
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
		xval = fcast_get_ldv(fc, yno, t, lag, dset);
	    } else {
		xval = dset->Z[v][t];
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
   uncertainty, as is common), via recursion. Cf. Box and Jenkins,
   1976, p. 508, "Program 4", V(l) algorithm (with the sign of theta
   changed).
*/

static double arma_variance (const double *phi, int p,
			     const double *theta, int q,
			     double *psi, int npsi, int l)
{
    /* the sum of squared psi's */
    static double sspsi;

    DPRINTF(("arma_variance: p=%d, q=%d, npsi=%d, l=%d\n", p, q, npsi, l));

    if (l == 1) {
	int i, j;

	psi[0] = 1.0;
	for (j=1; j<npsi; j++) {
	    psi[j] = 0.0;
	    for (i=1; i<=j; i++) {
		if (i <= p) {
		    psi[j] += phi[i] * psi[j-i];
		}
	    }
	    if (theta != NULL && j <= q) {
		psi[j] += theta[j];
	    }
	    DPRINTF(("psi[%d] = %g\n", j, psi[j]));
	}
	sspsi = 1.0;
    } else {
	sspsi += psi[l-1] * psi[l-1];
    }

    DPRINTF(("augmented 'sspsi' using psi(%d) = %g (sspsi = %g)\n",
	     l-1, psi[l-1], sspsi));

    return sspsi;
}

static double arima_difference_obs (const double *x, int t,
				    int *c, int k)
{
    double dx = x[t];

    if (!na(dx)) {
	int i, s;

	for (i=0; i<k; i++) {
	    if (c[i] != 0) {
		s = t - i - 1;
		if (s < 0 || na(x[s])) {
		    dx = NADBL;
		    break;
		} else {
		    dx -= c[i] * x[s];
		}
	    }
	}
    }

    return dx;
}

/* When forecasting based on an armax model estimated using X12A,
   or via the Kalman filter, we need to form the series X\beta
   so that we can subtract X\beta_{t-i} from y_{t-i} in
   computing the AR portion of the forecast. "beta" below
   is the array of ARMAX coefficients, not including the
   constant.

   If the model is ARIMAX and was NOT estimated with the
   --y-diff-only option (OPT_Y) then the regressors X have
   to be differenced in this context.
*/

static double *create_Xb_series (Forecast *fc, const MODEL *pmod,
				 const double *beta, const int *xlist,
				 const double **Z)
{
    double *Xb;
    int *delta = NULL;
    double x;
    int miss, diff = 0;
    int d, D, s = 0, k = 0;
    int i, j, t, vi;

    Xb = malloc((fc->t2 + 1) * sizeof *Xb);
    if (Xb == NULL) {
	return NULL;
    }

    d = gretl_model_get_int(pmod, "arima_d");
    D = gretl_model_get_int(pmod, "arima_D");

    if ((d > 0 || D > 0) && !(pmod->opt & OPT_Y)) {
	s = gretl_model_get_int(pmod, "arma_pd");
	diff = 1;
    }

    if (diff) {
	delta = arima_delta_coeffs(d, D, s);
	if (delta == NULL) {
	    free(Xb);
	    return NULL;
	}
	k = d + s * D;
    }

    for (t=0; t<=fc->t2; t++) {
	Xb[t] = 0.0;
	miss = 0;
	j = 0;
	for (i=1; i<=xlist[0] && !miss; i++) {
	    vi = xlist[i];
	    if (vi == 0) {
		Xb[t] += pmod->coeff[0];
	    } else {
		if (diff) {
		    x = arima_difference_obs(Z[vi], t, delta, k);
		} else {
		    x = Z[vi][t];
		}
		if (na(x)) {
		    Xb[t] = NADBL;
		    miss = 1;
		} else {
		    Xb[t] += beta[j++] * x;
		}
	    }
	}
    }

    if (delta != NULL) {
	free(delta);
    }

    return Xb;
}

static int want_x_beta_prep (const MODEL *pmod, const int *xlist)
{
    int ret = 0;

    if (xlist != NULL) {
	int aflags = gretl_model_get_int(pmod, "arma_flags");

	if (aflags & (ARMA_EXACT | ARMA_X12A)) {
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
		       const DATASET *dset)
{
    double *psi = NULL;
    gretl_vector *phi = NULL;
    gretl_vector *theta = NULL;
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
	     dset->varname[yno], fc->t1, fc->t2, p, q));

    xvars = (xlist != NULL)? xlist[0] : 0;

    err = arma_model_AR_MA_coeffs(pmod, &phi, &theta, OPT_I);
    if (err) {
	goto bailout;
    }

    beta = arma_model_get_x_coeffs(pmod);

#if AR_DEBUG
    fprintf(stderr, "beta = %p\n", (void *) beta);
    printlist(xlist, "xlist");
#endif

    if (want_x_beta_prep(pmod, xlist)) {
	if (gretl_is_arima_model(pmod)) {
	    regarma = 1;
	    err = regarma_model_AR_coeffs(pmod, &phi0, &px);
	}
	if (xlist[0] == 1 && xlist[1] == 0) {
	    /* just a const, no ARMAX */
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
	    Xb = create_Xb_series(fc, pmod, beta, xlist,
				  (const double **) dset->Z);
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
	ar_smax = MAX(fc->t1 - 1, p - 1);
	ma_smax = MAX(fc->t1 - 1, q);
    } else {
	ar_smax = pmod->t2;
	ma_smax = pmod->t2;
    }

    DPRINTF(("ar_smax = %d, ma_smax = %d\n", ar_smax, ma_smax));

    /* dependent variable */
    y = dset->Z[yno];

    /* do real forecast */
    for (t=fcstart; t<=fc->t2 && !err; t++) {
	double yh = 0.0;
	int miss = 0;

#if AR_DEBUG
	char obsstr[OBSLEN];
	ntolabel(obsstr, t, dset);
#endif
	DPRINTF(("\n *** Doing forecast for obs %d (%s)\n", t, obsstr));

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
		    xval = dset->Z[xlist[i]][t];
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

	DPRINTF((" Xb contribution: yh = %g\n", yh));

	/* AR contribution (incorporating any differencing) */

	for (i=1; i<=p && !miss; i++) {
	    if (phi->val[i] == 0.0) {
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
		DPRINTF(("  AR: lag %d, using coeff %#.8g\n", i, phi->val[i]));
		if (!regarma && Xb != NULL) {
		    if (na(Xb[s])) {
			miss = 1;
		    } else {
			yh += phi->val[i] * (yval - Xb[s]);
		    }
		} else if (!regarma && !na(mu)) {
		    DPRINTF(("    !regarma && !na(mu): yh += %g * (%g - %g)\n",
			     phi->val[i], yval, mu));
		    yh += phi->val[i] * (yval - mu);
		} else {
		    yh += phi->val[i] * yval;
		}
	    }
	}

	DPRINTF((" with AR contribution: %g\n", yh));

	/* MA contribution */

	for (i=1; i<=q && !miss; i++) {
	    if (theta->val[i] == 0.0) {
		continue;
	    }
	    s = t - i;
	    if (s >= pmod->t1 && s <= ma_smax) {
		DPRINTF(("  MA: lag %d, e[%d] = %g, theta[%d] = %g\n", i, s,
			 pmod->uhat[s], i, theta->val[i]));
		if (na(pmod->uhat[s])) {
		    char obsstr[OBSLEN];

		    ntolabel(obsstr, s, dset);
		    fprintf(stderr, "Uh oh, pmod->uhat[%d] is missing! (%s)\n",
			    s, obsstr);
		} else {
		    yh += theta->val[i] * pmod->uhat[s];
		}
	    } else if (fc->eps != NULL) {
		DPRINTF(("  MA: lag %d, ehat[%d] = %g, theta[%d] = %g\n", i, s,
			 fc->eps[s], i, theta->val[i]));
		yh += theta->val[i] * fc->eps[s];
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
	    const double *arvec = (phi != NULL)? phi->val : NULL;
	    const double *mavec = (theta != NULL)? theta->val : NULL;

	    vl = arma_variance(arvec, p, mavec, q,
			       psi, npsi, t - fcstart + 1);
	    fc->sderr[t] = pmod->sigma * sqrt(vl);
	}
    }

 bailout:

    free(psi);
    free(Xb);
    free(phi0);

    gretl_vector_free(phi);
    gretl_vector_free(theta);

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
set_up_ar_fcast_variance (const MODEL *pmod, int pmax, int npsi,
			  double **pphi, double **ppsi,
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
		     const DATASET *dset)
{
    const int *arlist;
    double *phi = NULL;
    double *psi = NULL;
    double *errphi = NULL;
    double *rho;
    double xval, yh, vl;
    double ylag, xlag;
    int miss, yno;
    int i, k, v, t, tk;
    int p, dvlag, pmax = 0;
    int pwe, npsi = 0;
    int err = 0;

#if AR_DEBUG
    fprintf(stderr, "\n*** ar_fcast, method = %d\n\n", fc->method);
#endif

    yno = pmod->list[1];

    if (pmod->ci == AR1 && pmod->arinfo == NULL) {
	dummy_ar_info_init(pmod);
    }

    arlist = pmod->arinfo->arlist;
    p = arlist[arlist[0]]; /* AR order of error term */
    rho = pmod->arinfo->rho;

    if (fc->t2 > pmod->t2 && fc->sderr != NULL) {
	/* we'll compute variance only if we're forecasting out of
	   sample */
	pmax = max_ar_lag(fc, pmod, p);
	npsi = fc->t2 - fc->t1 + 1;
	/* npsi = fc->t2 - pmod->t2; */
	DPRINTF(("pmax = %d, npsi = %d\n", pmax, npsi));
	set_up_ar_fcast_variance(pmod, pmax, npsi, &phi, &psi, &errphi);
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
	    tk = t - arlist[k];
	    ylag = fcast_get_ldv(fc, yno, tk, 0, dset);
	    if (na(ylag)) {
		miss = 1;
	    } else {
		yh += rho[k-1] * ylag;
	    }
	}

	/* (1 - r(L)) X_t b */
	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    v = pmod->list[i+2];
	    if ((dvlag = depvar_lag(fc, i))) {
		xval = fcast_get_ldv(fc, yno, t, dvlag, dset);
	    } else {
		xval = dset->Z[v][t];
	    }
	    if (na(xval)) {
		miss = 1;
	    } else {
		if (dvlag > 0 && phi != NULL) {
		    /* augment phi for computation of variance */
		    phi[dvlag] += pmod->coeff[i];
		}
		for (k=1; k<=arlist[0]; k++) {
		    tk = t - arlist[k];
		    if (dvlag > 0) {
			xlag = fcast_get_ldv(fc, yno, tk, dvlag, dset);
		    } else {
			xlag = dset->Z[v][tk];
		    }
		    if (!na(xlag)) {
			xval -= rho[k-1] * xlag;
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
		DPRINTF(("sderr[%d] = %g * sqrt(%g) = %g\n", t, pmod->sigma,
			 vl, fc->sderr[t]));
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
   type LOGISTIC, LOGIT, PROBIT, TOBIT, POISSON and NEGBIN.
 */

static double fcast_transform (double xb, const MODEL *pmod,
			       int t, const double *offset,
			       double lmax)
{
    int ymin, ci = pmod->ci;
    double yf = xb;

    if (ci == TOBIT) {
	if (xb < 0.0) {
	    yf = 0.0;
	}
    } else if (ci == LOGIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    ymin = gretl_model_get_int(pmod, "ymin");
	    yf = ordered_model_prediction(pmod, xb, ymin);
	} else {
	    yf = exp(xb) / (1.0 + exp(xb));
	}
    } else if (ci == PROBIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    ymin = gretl_model_get_int(pmod, "ymin");
	    yf = ordered_model_prediction(pmod, xb, ymin);
	} else {
	    yf = normal_cdf(xb);
	}
    } else if (ci == LOGISTIC) {
	if (na(lmax)) {
	    yf = 1.0 / (1.0 + exp(-xb));
	} else {
	    yf = lmax / (1.0 + exp(-xb));
	}
    } else if (ci == POISSON || ci == NEGBIN) {
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

static int
integrated_fcast (Forecast *fc, const MODEL *pmod, int yno,
		  const DATASET *dset)
{
    double s2 = NADBL;
    int dyn_start; /* start of dynamic forecast */
    double xval, yht = 0.0;
    int i, vi, t;

    if (fc->method == FC_STATIC) {
	dyn_start = dset->n; /* i.e., never */
    } else if (fc->method == FC_DYNAMIC) {
	dyn_start = fc->t1;
    } else {
	/* FC_AUTO: dynamic out of sample */
	dyn_start = pmod->t2 + 1;
    }

    if (fc->sderr != NULL) {
	s2 = pmod->sigma * pmod->sigma;
	fc->sderr[dyn_start - 1] = 0.0;
    }

    for (t=fc->t1; t<=fc->t2; t++) {
	int miss = 0;

	if (t <= dyn_start) {
	    yht = dset->Z[yno][t-1];
	    if (na(yht)) {
		miss = 1;
	    }
	}

	for (i=0; i<pmod->ncoeff && !miss; i++) {
	    vi = pmod->list[i+2];
	    xval = dset->Z[vi][t];
	    if (na(xval)) {
		miss = 1;
	    } else {
		yht += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    if (t >= dyn_start) {
		/* all subsequent values are NA */
		int s;

		for (s=t; s<=fc->t2; s++) {
		    fc->yhat[s] = NADBL;
		}
		break;
	    } else {
		fc->yhat[t] = NADBL;
	    }
	} else {
	    fc->yhat[t] = yht;
	    if (fc->sderr != NULL && t >= dyn_start) {
		fc->sderr[t] = fc->sderr[t-1] + s2;
	    }
	}
    }

    if (fc->sderr != NULL) {
	fc->sderr[dyn_start - 1] = NADBL;
	for (t=dyn_start; t<=fc->t2; t++) {
	    if (!na(fc->sderr[t])) {
		fc->sderr[t] = sqrt(fc->sderr[t]);
	    }
	}
    }

    return 0;
}

static int mlogit_allprobs (Forecast *fc, const MODEL *pmod,
			    const DATASET *dset)
{
    gretl_matrix *P;
    int err = 0;

    P = mn_logit_probabilities(pmod, fc->t1, fc->t2, dset, &err);
    if (P != NULL) {
	set_fcast_matrices(P, NULL);
    }

    return err;
}

static int mlogit_fcast (Forecast *fc, const MODEL *pmod,
			 const DATASET *dset)
{
    const gretl_matrix *yvals;
    gretl_matrix *Xt;
    int k, i, vi, t;

    yvals = gretl_model_get_data(pmod, "yvals");
    if (yvals == NULL) {
	return E_DATA;
    }

    k = pmod->list[0] - 1;
    Xt = gretl_matrix_alloc(1, k);
    if (Xt == NULL) {
	return E_ALLOC;
    }

    for (t=fc->t1; t<=fc->t2; t++) {
	double xti, yht;
	int miss = 0;

	for (i=0; i<k && !miss; i++) {
	    vi = pmod->list[i+2];
	    xti = dset->Z[vi][t];
	    if (na(xti)) {
		miss = 1;
	    } else {
		Xt->val[i] = xti;
	    }
	}

	if (miss) {
	    fc->yhat[t] = NADBL;
	} else {
	    yht = mn_logit_prediction(Xt, pmod->coeff, yvals);
	    fc->yhat[t] = yht;
	}
    }

    gretl_matrix_free(Xt);

    return 0;
}

/* Compute forecasts for linear models without autoregressive errors.
   We don't compute forecast standard errors, unless we're integrating
   the forecast from a static model estimated via OLS.
*/

static int linear_fcast (Forecast *fc,
                         const MODEL *pmod,
                         int yno,
			 const DATASET *dset,
                         const gretl_matrix *X,
                         gretlopt opt)
{
    const double *offvar = NULL;
    double xval, lmax = NADBL;
    int all_probs = 0;
    int k = pmod->ncoeff;
    int i, vi, t;
    int err = 0;

    if (COUNT_MODEL(pmod->ci) && dset != NULL) {
	/* special for "offset" variable */
	int offnum = gretl_model_get_int(pmod, "offset_var");

	if (offnum > 0) {
	    offvar = dset->Z[offnum];
	}
    } else if (pmod->ci == LOGISTIC) {
	lmax = gretl_model_get_double(pmod, "lmax");
    } else if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    /* we need to know how many coeffs are not
	       just estimated cut-points */
	    k = gretl_model_get_int(pmod, "nx");
	    if (k <= 0) {
		return E_MISSDATA;
	    } else if (opt & OPT_L) {
		all_probs = 1;
	    }
	}
    }

    for (t=fc->t1; t<=fc->t2; t++) {
	double yht = 0.0;
	int miss = 0;

	for (i=0; i<k && !miss; i++) {
	    int lag = depvar_lag(fc, i);

	    if (lag != 0) {
		xval = fcast_get_ldv(fc, yno, t, lag, dset);
	    } else if (dset != NULL) {
                vi = pmod->list[i+2];
		xval = dset->Z[vi][t];
	    } else {
                xval = X->val[i * X->rows + t];
            }
	    if (na(xval)) {
		miss = 1;
	    } else {
		yht += xval * pmod->coeff[i];
	    }
	}

	if (miss) {
	    fc->yhat[t] = NADBL;
	} else if (FCAST_SPECIAL(pmod->ci) && !all_probs) {
	    /* special handling for LOGIT and others */
	    fc->yhat[t] = fcast_transform(yht, pmod, t, offvar, lmax);
	} else {
	    fc->yhat[t] = yht;
	}
    }

    if (all_probs && dset != NULL) {
	gretl_matrix *P;

	P = ordered_probabilities(pmod, fc->yhat, fc->t1, fc->t2,
				  dset, &err);
	set_fcast_matrices(P, NULL);
    }

    return err;
}

#define dynamic_nls(m) ((m->ci == NLS || m->ci == MIDASREG) \
			&& gretl_model_get_int(m, "dynamic"))

static int get_forecast_method (Forecast *fc,
				MODEL *pmod,
				const DATASET *dset,
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
       unless OPT_S for "static" has been given
    */
    if (dataset_is_time_series(dset) && !(opt & OPT_S) && pmod->ci != ARMA) {
	fc->dvlags = process_lagged_depvar(pmod, dset);
    }

    if (!(opt & OPT_S)) {
	/* user didn't give the "static" option */
	if (pmod->ci == ARMA || fc->dvlags != NULL ||
	    dynamic_nls(pmod) || (opt & OPT_I)) {
	    /* dynamic forecast is possible */
	    dyn_ok = 1;
	}
	if (SIMPLE_AR_MODEL(pmod->ci) || pmod->ci == GARCH) {
	    dyn_errs_ok = 1;
	}
    }

    if (!dyn_ok && (opt & OPT_D)) {
	/* "dynamic" option given, but can't be honored */
	return inapplicable_option_error(FCAST, OPT_D);
    }

    /* NLS: we can only do dynamic out of sample (fc->t1 > pmod->t2) */
    if (pmod->ci == NLS && dyn_ok && (opt & OPT_D)) {
	if (fc->t1 <= pmod->t2) {
	    return inapplicable_option_error(FCAST, OPT_D);
	}
    }

    if (opt & OPT_D) {
	/* user requested dynamic forecast and it seems OK */
	fc->method = FC_DYNAMIC;
    } else if ((dyn_ok || dyn_errs_ok) && fc->t2 > pmod->t2) {
	/* do dynamic f'cast out of sample */
	fc->method = FC_AUTO;
    }

    if (fc->method == FC_DYNAMIC || fc->method == FC_AUTO) {
	if (fc->t1 > fc->model_t2 + 1) {
	    /* we'll probably need intervening forecasts */
	    fc->t1 = fc->model_t2 + 1;
	}
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

/* find the earlist starting point at which we have a previous
   valid observation on the level of Z[yno] with which to
   initialize an integrated forecast
*/

static int revise_fr_start (FITRESID *fr, int yno, const double **Z,
			    const DATASET *dset)
{
    int t, t0;

    if (fr->t0 == 0) {
	/* certainly can't start before obs 1 */
	fr->t0 = 1;
    }

    /* previous obs */
    t0 = fr->t0 - 1;

    /* find first non-missing value */
    for (t=t0; t<dset->n; t++) {
	if (!na(Z[yno][t])) {
	    break;
	}
    }

    t0 = t;

    if (t0 >= dset->n - 2) {
	return E_MISSDATA;
    }

    fr->t0 = t0 + 1;

    if (fr->t1 < fr->t0) {
	fr->t1 = fr->t0;
    }

    if (fr->t2 < fr->t1) {
	fr->t2 = fr->t1;
    }

    fr->fitted[t0] = Z[yno][t0];
    fr->resid[t0] = 0.0;

    return 0;
}

/* check whether we really can do an integrated forecast */

static int check_integrated_forecast_option (MODEL *pmod,
					     DATASET *dset,
					     int *pyno)
{
    int err = 0;

    if (pmod->ci != OLS) {
	err = inapplicable_option_error(FCAST, OPT_I);
    } else {
	int yno = gretl_model_get_depvar(pmod);
	int d, parent;

	d = is_standard_diff(yno, dset, &parent);
	if (!d) {
	    err = inapplicable_option_error(FCAST, OPT_I);
	} else if (pyno != NULL) {
	    /* in caller, make yno = ID of level variable */
	    *pyno = parent;
	}
    }

    return err;
}

static int check_all_probs_option (MODEL *pmod)
{
    int err = E_BADOPT;

    if (pmod->ci == LOGIT) {
	if (gretl_model_get_int(pmod, "multinom") ||
	    gretl_model_get_int(pmod, "ordered")) {
	    err = 0;
	}
    } else if (pmod->ci == PROBIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    err = 0;
	}
    }

    return err;
}

/* driver for various functions that compute forecasts
   for different sorts of models */

static int real_get_fcast (FITRESID *fr, MODEL *pmod,
			   DATASET *dset, gretlopt opt)
{
    Forecast fc;
    const double **Z = (const double **) dset->Z;
    int yno = gretl_model_get_depvar(pmod);
    int integrate = (opt & OPT_I);
    int dummy_AR = 0;
    int DM_errs = 0;
    int dyn_errs = 0;
    int asy_errs = 0;
    int int_errs = 0;
    int same_data = 0;
    int nf = 0;
    int t, err;

    if (opt & OPT_L) {
	err = check_all_probs_option(pmod);
	if (err) {
	    return err;
	}
    }

    forecast_init(&fc);

    if (integrate) {
	err = check_integrated_forecast_option(pmod, dset, &yno);
	if (!err) {
	    err = revise_fr_start(fr, yno, Z, dset);
	}
	if (err) {
	    return err;
	}
    }

    fc.t1 = fr->t1;
    fc.t2 = fr->t2;
    fc.model_t2 = pmod->t2;

    err = get_forecast_method(&fc, pmod, dset, opt);
    if (err) {
	return err;
    }

    if (pmod->ci == NLS && fc.method == FC_STATIC) {
	asy_errs = 1;
    }

    if (pmod->ci == PANEL || pmod->ci == DPANEL) {
	; /* don't do the things below */
    } else if (!FCAST_SPECIAL(pmod->ci) && !integrate) {
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

    if (integrate && pmod->ci == OLS &&
	fc.dvlags == NULL && fc.method != FC_STATIC) {
	int_errs = 1;
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

    if (DM_errs || dyn_errs || asy_errs || int_errs) {
	err = fit_resid_add_sderr(fr);
    }

    if (err) {
	return err;
    }

    fc.yhat = fr->fitted;
    fc.sderr = fr->sderr;

    if (pmod->ci == ARMA && fc.method == FC_STATIC) {
	fc_add_eps(&fc, dset->n);
    }

    /* compute the actual forecast */
    if (pmod->ci == DPANEL || individual_effects_model(pmod)) {
	err = panel_fcast(&fc, pmod, dset, opt);
    } else if (DM_errs) {
	err = static_fcast_with_errs(&fc, pmod, dset, opt);
    } else if (pmod->ci == NLS) {
	err = nls_fcast(&fc, pmod, dset);
    } else if (pmod->ci == MIDASREG) {
	err = midas_fcast(&fc, pmod, dset);
    } else if (SIMPLE_AR_MODEL(pmod->ci) || dummy_AR) {
	err = ar_fcast(&fc, pmod, dset);
    } else if (pmod->ci == ARMA) {
	err = arma_fcast(&fc, pmod, dset);
    } else if (pmod->ci == GARCH) {
	err = garch_fcast(&fc, pmod, dset);
    } else if (integrate) {
	err = integrated_fcast(&fc, pmod, yno, dset);
    } else if (pmod->ci == LOGIT && gretl_model_get_int(pmod, "multinom")) {
	if (opt & OPT_L) {
	    /* --all-probs */
	    err = mlogit_allprobs(&fc, pmod, dset);
	} else {
	    err = mlogit_fcast(&fc, pmod, dset);
	}
    } else {
	err = linear_fcast(&fc, pmod, yno, dset, NULL, opt);
    }

    /* free any auxiliary info */
    forecast_free(&fc);

    if (opt & OPT_L) {
        /* --all-probs: we're done */
	return err;
    }

    if (dummy_AR) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
    }

    same_data = same_dataset(pmod, dset);

    for (t=0; t<fr->nobs; t++) {
	if (t >= fr->t1 && t <= fr->t2) {
	    if (!na(fr->fitted[t])) {
		nf++;
	    }
	} else if (same_data && t >= fr->t0 && t <= fr->t2 &&
		   t >= pmod->t1 && t <= pmod->t2) {
	    if (integrate) {
		fr->fitted[t] = fr->fitted[t-1] + pmod->yhat[t];
		fr->resid[t] = fr->resid[t-1] + pmod->uhat[t];
	    } else {
		fr->fitted[t] = pmod->yhat[t];
		fr->resid[t] = pmod->uhat[t];
	    }
	}
	fr->actual[t] = dset->Z[yno][t];
    }

    if (nf == 0) {
	err = E_MISSDATA;
    } else {
	fit_resid_set_dec_places(fr);
	strcpy(fr->depvar, dset->varname[yno]);
	fr->df = pmod->dfd;
    }

    return err;
}

static int matrix_forecast_supported (MODEL *pmod)
{
    const int *xlist;

    if (pmod->errcode) {
        return 0;
    } else if (pmod->ncoeff < 1) {
        return 0;
    } else if (gretl_model_get_int(pmod, "ldepvar")) {
        return 0;
    } else if (pmod->ci != OLS &&
               pmod->ci != LOGIT &&
               pmod->ci != PROBIT) {
        /* FIXME relax this! */
        return 0;
    } else if ((xlist = model_xlist(pmod)) == NULL) {
        return 0;
    }

    if (pmod->ncoeff != xlist[0]) {
        return 0;
    } else {
        fprintf(stderr, "HERE ncoeff = nx = %d, OK?\n", pmod->ncoeff);
    }

    return 1;
}

/* April 2025: simple forecast interface: X-data given in matrix
   form and prediction returned in matrix form.
*/

gretl_matrix *matrix_forecast (MODEL *pmod,
                               const gretl_matrix *X,
                               int *err)
{
    gretl_matrix *F = NULL;
    FITRESID *fr = NULL;
    Forecast fc;
    int n = 0;

    if (!matrix_forecast_supported(pmod)) {
        *err = E_DATA;
    } else if (X == NULL || X->rows < 1) {
        *err = E_DATA;
    } else {
        n = X->rows;
        fr = fit_resid_new_for_model(pmod, NULL, 0, n, 0, err);
    }

    if (!*err) {
        int yno = pmod->list[1]; /* FIXME */

        forecast_init(&fc);
        fc.method = FC_STATIC;
        fc.yhat = fr->fitted;
        fc.sderr = fr->sderr;
        fc.t1 = 0;
        fc.t2 = n - 1;
        *err = linear_fcast(&fc, pmod, yno, NULL, X, OPT_NONE);
    }

    if (!*err) {
        int cols = fc.sderr == NULL ? 1 : 2;
        int i;

        F = gretl_zero_matrix_new(n, cols);
        if (F == NULL) {
            *err = E_ALLOC;
        } else {
            for (i=0; i<n; i++) {
                gretl_matrix_set(F, i, 0, fc.yhat[i]);
                if (cols > 1) {
                    gretl_matrix_set(F, i, 1, fc.sderr[i]);
                }
            }
        }
    }

    free_fit_resid(fr);

    return F;
}

static int fcast_get_limit (const char *s, DATASET *dset)
{
    double x;
    int t, err = 0;

    if (gretl_is_scalar(s)) {
	x = gretl_scalar_get_value(s, NULL);
    } else {
	x = generate_scalar(s, dset, &err);
    }

    if (x < 1 || x > dset->n) {
	gretl_error_clear();
	gretl_errmsg_set(_("Observation number out of bounds"));
	t = -1;
    } else {
	t = x - 1;
    }

    return t;
}

enum {
    OS_OK = 0,
    OS_ERR,
    OS_PANEL
};

static int out_of_sample_check (MODEL *pmod, DATASET *dset)
{
    int panel, ret;

    if (pmod->ci == PANEL || pmod->ci == DPANEL ||
	(pmod->ci == OLS && dataset_is_panel(dset))) {
	panel = 1;
	ret = OS_ERR;
	if (gretl_model_get_int(pmod, "ntdum") > 0) {
	    gretl_errmsg_set(_("Specification includes time dummies: cannot "
			       "forecast out of sample"));
	    return ret;
	}
    } else {
	panel = 0;
	ret = OS_OK;
    }

    if (pmod->smpl.t1 == 0 && pmod->smpl.t2 == dset->n - 1) {
	/* no out-of-sample obs reachable via t1, t2 */
	ret = OS_ERR;
	if (!panel) {
	    if (pmod->t2 < dset->n - 1) {
		/* 2019-11-15: may be OK for dynamic model */
		ret = OS_OK;
	    }
	} else {
	    /* the panel case */
	    DATASET *fullset = fetch_full_dataset();

	    if (dataset_is_panel(fullset)) {
		int Tfull = fullset->pd;
		int Tcurr = dset->pd;
		int Nfull = fullset->n / Tfull;
		int Ncurr = dset->n / Tcurr;

		if (Ncurr == Nfull && Tcurr < Tfull &&
		    dset->v == fullset->v) {
		    /* sub-sampled in the time dimension only, and
		       no series added or deleted (we trust)
		    */
		    ret = OS_PANEL;
		}
	    }
	}
	if (ret == OS_ERR) {
	    gretl_errmsg_set(_("No out-of-sample observations are available"));
	}
    } else if (panel) {
	gretl_errmsg_set(_("The --out-of-sample option is only supported in the "
			   "time dimension"));
    }

    return ret;
}

static int parse_forecast_string (const char *s,
				  gretlopt opt,
				  MODEL *pmod,
				  int t2est,
				  DATASET *dset,
				  int *pt1, int *pt2,
				  int *pk, char *vname,
				  int *os_case)
{
    char f[4][32] = {0};
    char *t1str = NULL, *t2str = NULL;
    char *kstr = NULL, *vstr = NULL;
    int nmax, nmin = (opt & OPT_R)? 1 : 0;
    int t1 = 0, t2 = 0;
    int nf = 0;
    int err = 0;

    /* "static" and "recursive" can't be combined */
    err = incompatible_options(opt, OPT_S | OPT_R);
    if (err) {
	return err;
    }

    *vname = '\0';

    /* How many fields should we be looking for in the user input?
       If OPT_R ("recursive") is given, the max is 4:

       t1, t2, k (steps-ahead), vname

       Otherwise the max is 3 (no k-value):

       t1, t2, vname

       However, if OPT_O ("out-of-sample") is given, we should not
       expect t1 and t2, so the max becomes 2 (if OPT_R), else 1.

       Also note: t1 and t2 need not be given (even in the absence
       of OPT_O) since these default to the current sample range.
    */

    if (opt & OPT_O) {
	nmax = (opt & OPT_R)? 2 : 1;
    } else {
	nmax = (opt & OPT_R)? 4 : 3;
    }

    if (*s == '\0') {
	nf = 0;
    } else {
	nf = sscanf(s, "%31s %31s %31s %31s", f[0], f[1], f[2], f[3]);
	if (nf > nmax) {
	    /* try for parenthesized t1, t2 terms? */
	    nf = sscanf(s, "(%31[^)]) (%31[^)]) %31s %31s", f[0], f[1], f[2], f[3]);
	}
    }

    if (nf < nmin || nf > nmax) {
	gretl_errmsg_sprintf("fcast: expected %d to %d fields in input, got %d",
			     nmin, nmax, nf);
	return E_PARSE;
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

    if (t1str != NULL || t2str != NULL) {
	if (opt & OPT_O) {
	    /* t1, t2 should not be given with --out-of-sample */
	    gretl_errmsg_set("fcast: unexpected t1 and/or t2 field in input");
	    return E_DATA;
	}
    }

    if (kstr != NULL && pk != NULL) {
	*pk = gretl_int_from_string(kstr, &err);
	if (err) {
	    return err;
	}
    }

    if (vstr != NULL) {
	strncat(vname, vstr, VNAMELEN - 1);
    }

    if (t1str != NULL && t2str != NULL) {
	t1 = dateton(t1str, dset);
	if (t1 < 0) {
	    t1 = fcast_get_limit(t1str, dset);
	}
	t2 = dateton(t2str, dset);
	if (t2 < 0) {
	    t2 = fcast_get_limit(t2str, dset);
	}
	if (t1 < 0 || t2 < 0 || t2 < t1) {
	    err = E_DATA;
	}
    } else if (opt & OPT_O) {
	/* out of sample, if possible */
	if (pmod != NULL) {
	    *os_case = out_of_sample_check(pmod, dset);
	}
	if (*os_case == OS_ERR) {
	    err = E_DATA;
	} else if (*os_case == OS_OK) {
	    if (dset->n - t2est - 1 > 0) {
		t1 = t2est + 1;
		t2 = dset->n - 1;
	    } else {
		err = E_OBS;
	    }
	}
    } else {
	/* default: current sample range */
	t1 = dset->t1;
	t2 = dset->t2;
    }

    if (!err) {
	/* in case we hit any "temporary" errors above */
	gretl_error_clear();
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
 * @dset: dataset struct.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where a dynamic forecast is meaningful).
 * If @opt includes %OPT_I, integrate the forecast (only relevant
 * if the dependent variable in the model in question is recognized
 * as the first difference of another variable).
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
			DATASET *dset, gretlopt opt, int *err)
{
    FITRESID *fr;

    if (pmod->errcode) {
	*err = E_DATA;
	return NULL;
    }

    fr = fit_resid_new_for_model(pmod, dset, t1, t2, pre_n, err);

    if (!*err) {
	*err = real_get_fcast(fr, pmod, dset, opt);
	if (*err) {
	    free_fit_resid(fr);
	    fr = NULL;
	}
    }

    return fr;
}

void forecast_matrix_cleanup (void)
{
    set_fcast_matrices(NULL, NULL);
}

gretl_matrix *get_forecast_matrix (int idx, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *M = NULL;

    if (idx == M_FCAST) {
	M = fcast_matrix;
    } else if (idx == M_FCSE) {
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

    gretl_matrix_set_t1(f, ft1);
    gretl_matrix_set_t2(f, ft2);

    if (fr->sderr != NULL) {
	T = et2 - et1 + 1;
	if (T > 0) {
	    e = gretl_matrix_alloc(T, 1);
	    if (e == NULL) {
		err = E_ALLOC;
	    } else {
		gretl_matrix_set_t1(e, et1);
		gretl_matrix_set_t2(e, et2);
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

    set_fcast_matrices(f, e);

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
    int i, j, Ft1, mt1;
    int err = 0;

    Ft1 = gretl_matrix_get_t1(F);

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

    mt1 = Ft1 + f0;
    gretl_matrix_set_t1(f, mt1);
    gretl_matrix_set_t2(f, mt1 + fT - 1);

    eT = en - e0 + 1;

    if (eT > 0) {
	e = gretl_matrix_alloc(eT, k);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    mt1 = Ft1 + e0;
	    gretl_matrix_set_t1(e, mt1);
	    gretl_matrix_set_t2(e, mt1 + eT - 1);
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

    set_fcast_matrices(f, e);

    if (F1 != NULL) {
	gretl_matrix_free(F1);
    }

    return err;
}

/* compatibility with behaviour of old fcast command */

static int add_fcast_to_dataset (FITRESID *fr, const char *vname,
				 DATASET *dset, PRN *prn)
{
    int oldv = dset->v;
    int v, err = 0;

    v = series_index(dset, vname);

    if (v == dset->v) {
	/* new variable */
	if (check_varname(vname)) {
	    err = E_DATA;
	} else {
	    err = dataset_add_series(dset, 1);
	    if (!err) {
		strcpy(dset->varname[v], vname);
	    }
	}
    }

    if (!err) {
	int t;

	for (t=0; t<dset->n; t++) {
	    dset->Z[v][t] = fr->fitted[t];
	}

	series_set_label(dset, v, _("predicted values"));

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

static void print_probs_matrix (PRN *prn)
{
    if (fcast_matrix != NULL) {
        if (gretl_echo_on()) {
            pputc(prn, '\n');
        }
        pprintf(prn, "%s\n\n", _("Estimated outcome probabilities"));
        gretl_matrix_print_to_prn(fcast_matrix, NULL, prn);
    }
}

static int add_fcast_matrix (const char *name,
			     DATASET *dset,
			     PRN *prn)
{
    int err = 0;

    if (gretl_reserved_word(name)) {
        err = E_TYPES;
    } else if (current_series_index(dset, name) >= 0) {
	err = E_TYPES;
    } else {
	err = user_var_add_or_replace(name, GRETL_TYPE_MATRIX,
				      fcast_matrix);
	fcast_matrix = NULL;
    }

    return err;
}

static int model_do_forecast (const char *str, MODEL *pmod,
			      DATASET *dset, gretlopt opt,
			      PRN *prn)
{
    char vname[VNAMELEN];
    FITRESID *fr = NULL;
    int t1, t2, k = -1;
    int all_probs = (opt & OPT_L);
    int os_case = 0;
    int quiet = 0;
    int err;

    /* FIXME relax some of these limitations? */
    if (pmod->ci == HECKIT || pmod->ci == DURATION) {
	return E_NOTIMP;
    } else if (pmod->ci == PANEL && !(pmod->opt & OPT_P)) {
	if (opt & OPT_R) {
	    return E_NOTIMP;
	}
    } else if (pmod->ci == PROBIT && (pmod->opt & OPT_E)) {
	return E_NOTIMP;
    }

    /* can't use --plot with --all-probs */
    if (all_probs && (opt & OPT_U)) {
        return E_BADOPT;
    }

    /* OPT_I for integrate: reject for non-OLS, or if the dependent
       variable is not recognized as a first difference */
    if (opt & OPT_I) {
	err = check_integrated_forecast_option(pmod, dset, NULL);
	if (err) {
	    return err;
	}
    }

    err = parse_forecast_string(str, opt, pmod, pmod->t2, dset,
				&t1, &t2, &k, vname, &os_case);
    if (err) {
	return err;
    }

    if (*vname != '\0') {
	/* saving to named series or matrix: be quiet */
	opt |= (OPT_A | OPT_Q);
    }

    quiet = (opt & OPT_Q);

    if (os_case == OS_PANEL) {
	return panel_os_special(pmod, dset, vname, opt, prn);
    } else if (opt & OPT_R) {
	fr = recursive_OLS_k_step_fcast(pmod, dset, t1, t2,
					k, 0, &err);
    } else {
	fr = get_forecast(pmod, t1, t2, 0, dset, opt, &err);
    }

    if (err) {
        goto bailout;
    }

    if (!err && (opt & OPT_A)) {
        /* named output */
	if (opt & OPT_L) {
            /* add forecast matrix directly */
	    err = add_fcast_matrix(vname, dset, prn);
	} else {
            /* add forecast series directly */
	    err = add_fcast_to_dataset(fr, vname, dset, prn);
	}
    }

    if (!quiet && (opt & OPT_L)) {
        print_probs_matrix(prn);
    } else if (!quiet) {
        gretlopt printopt = opt;

        if (opt & OPT_U) {
            /* do graph (from command line) */
            printopt |= OPT_P;
        }
        err = text_print_forecast(fr, dset, printopt, prn);
    }

    if (!err && !(opt & OPT_L)) {
	set_forecast_matrices_from_fr(fr);
    }

 bailout:

    free_fit_resid(fr);

    return err;
}

static int get_sys_fcast_var (const int *ylist, const char *vname,
			      DATASET *dset)
{
    int i, vi;

    for (i=0; i<ylist[0]; i++) {
	vi = ylist[i+1];
	if (!strcmp(vname, dset->varname[vi])) {
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
				 DATASET *dset)
{
    int m = F->cols / 2;
    int s, t, nf;
    int err = 0;

    strcpy(fr->depvar, dset->varname[yno]);

    /* "pre-forecast" observations */
    for (t=fr->t0; t<fr->t1; t++) {
	fr->actual[t] = dset->Z[yno][t];
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
	fr->actual[t] = dset->Z[yno][t];
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
			       DATASET *dset, gretlopt opt,
			       PRN *prn)
{
    char vname[VNAMELEN];
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const int *ylist = NULL;
    const gretl_matrix *F = NULL;
    int t1 = 0, t2 = 0;
    int t2est, ci;
    int imax, imin = 0;
    int os_case = 0;
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

    err = parse_forecast_string(str, opt, NULL, t2est, dset,
				&t1, &t2, NULL, vname, &os_case);

    if (!err) {
	if (var != NULL) {
	    F = gretl_VAR_get_forecast_matrix(var, t1, t2, dset,
					      opt, &err);
	} else {
	    F = system_get_forecast_matrix(sys, t1, t2, dset,
					   opt, &err);
	}
    }

    if (!err && *vname != '\0') {
	imin = get_sys_fcast_var(ylist, vname, dset);
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
	gretlopt printopt = (opt & OPT_N)? OPT_N : OPT_NONE;
	FITRESID *fr = NULL;
	int i, asy;

	asy = (ci == VECM);

	fr = fit_resid_new_for_system(asy, dset, t1, t2, 0, &err);
	if (err) {
	    return err;
	}

    	if (asy) {
	    /* asymptotic normal */
	    fr->df = var->T;
	} else {
	    fr->df = df;
	}

	if ((opt & OPT_U) && *vname != '\0') {
	    /* produce a plot */
	    printopt |= OPT_P;
	}

	for (i=imin; i<=imax && !err; i++) {
	    err = fill_system_forecast(fr, i, ylist[i+1], var, sys,
				       F, dset);
	    if (!err) {
		err = text_print_forecast(fr, dset, printopt, prn);
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
 * @dset: dataset struct.
 * @opt: if %OPT_D, force a dynamic forecast; if %OPT_S, force
 * a static forecast.  By default, the forecast is static within
 * the data range over which the model was estimated, and dynamic
 * out of sample (in cases where this distinction is meaningful).
 * %OPT_R: do recursive forecast.
 * %OPT_Q: suppress printing of the forecast;
 * %OPT_P: ensure that the values are printed.
 * %OPT_U: produce gnuplot plot.
 * @prn: gretl printing struct.
 *
 * In the case of "simple" models with an autoregressive error term
 * (%AR, %AR1) the forecast values incorporate the predictable
 * component of the error.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int do_forecast (const char *str, DATASET *dset,
		 gretlopt opt, PRN *prn)
{
    void *ptr;
    GretlObjType type;
    int err;

    ptr = get_last_model(&type);
    if (ptr == NULL) {
	return E_BADSTAT;
    }

    if ((opt & OPT_U) && dataset_is_panel(dset)) {
	gretl_errmsg_set(_("Forecast plot not implemented for panel data"));
	err = E_NOTIMP;
    } else if ((opt & (OPT_R | OPT_I)) && type != GRETL_OBJ_EQN) {
	/* "recursive", "integrate": single equations only */
	err = E_BADOPT;
    } else if (type == GRETL_OBJ_EQN) {
	err = model_do_forecast(str, ptr, dset, opt, prn);
    } else if (type == GRETL_OBJ_SYS || type == GRETL_OBJ_VAR) {
	err = system_do_forecast(str, ptr, type, dset, opt, prn);
    } else {
	err = E_DATA;
    }

    return err;
}

/* try to determine in advance how far we can go with a forecast,
   either dynamic or static (@ftype) */

static int fcast_get_t2max (const int *list, const int *dvlags,
			    const MODEL *pmod, const DATASET *dset,
			    int ftype)
{
    const double *ay = NULL;
    int i, vi, t, l0;

    if (pmod->ci == ARMA && ftype == FC_STATIC) {
	int yno = gretl_model_get_depvar(pmod);

	ay = dset->Z[yno];
    }

    l0 = list == NULL ? 0 : list[0];

    for (t=pmod->t2; t<dset->n; t++) {
	int p, vj, all_ok = 1;

	if (ay != NULL && na(ay[t-1])) {
	    /* FIXME? */
	    break;
	}

	for (i=1; i<=l0; i++) {
	    vi = list[i];
	    if (vi == 0) {
		continue;
	    } else if (dvlags != NULL && dvlags[i-1] != 0) {
		continue;
	    } else if (is_trend_variable(dset->Z[vi], dset->n)) {
		continue;
	    } else if (is_periodic_dummy(dset->Z[vi], dset)) {
		continue;
	    } else if ((p = series_get_lag(dset, vi)) > 0) {
		vj = series_get_parent_id(dset, vi);
		if (vj > 0 && !na(dset->Z[vj][t-p])) {
		    continue;
		}
	    }
	    if (na(dset->Z[vi][t])) {
		all_ok = 0;
		break;
	    }
	}

	if (!all_ok) {
	    t--;
	    break;
	} else if (t == dset->n - 1) {
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
 * @dset: dataset struct.
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
			       DATASET *dset, gretlopt opt,
			       int *err)
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
	F = gretl_VAR_get_forecast_matrix(var, t1, t2, dset, opt, err);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) p;
	yno = sys->ylist[i+1];
	df = sys->df;
	F = system_get_forecast_matrix(sys, t1, t2, dset, opt, err);
    } else {
	*err = E_DATA;
    }

    if (*err) {
	fprintf(stderr, "get_system_forecast: matrix F is NULL\n");
	return NULL;
    }

    asy = (ci == VECM);

    fr = fit_resid_new_for_system(asy, dset, t1, t2, pre_n, err);
    if (*err) {
	return NULL;
    }

    if (asy) {
	/* asymptotic normal */
	fr->df = var->T;
    } else {
	fr->df = df;
    }

    *err = fill_system_forecast(fr, i, yno, var, sys, F, dset);

    if (*err) {
	free_fit_resid(fr);
	fr = NULL;
    }

    return fr;
}

/* Try to determine whether adding observations to the dataset
   (without actually adding more "real" data) can serve to extend the
   range of out-of-sample prediction. The answer will in general be No
   if the specification includes exogenous regressors other than
   deterministic terms that can be extended automatically.  However,
   if an exogenous regressor is a lag series of order p we can extend
   it automatically for p periods.

   The @dvlags argument will be non-NULL only if the specification
   includes at least one lag of the dependent variable.
*/

static int addobs_can_help (MODEL *pmod, const int *dvlags,
			    const DATASET *dset)
{
    const int *xlist = model_xlist(pmod);
    int i, xi, ret = 1;

    if (xlist != NULL) {
	for (i=0; i<xlist[0]; i++) {
	    xi = xlist[i + 1];
	    if (xi != 0 && (dvlags == NULL || dvlags[i] == 0)) {
		if (is_trend_variable(dset->Z[xi], dset->n)) {
		    continue;
		} else if (is_periodic_dummy(dset->Z[xi], dset)) {
		    continue;
		} else if (series_get_lag(dset, xi) &&
			   series_get_parent_id(dset, xi)) {
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

/**
 * forecast_options_for_model:
 * @pmod: the model from which forecasts are wanted.
 * @dset: dataset struct.
 * @flags: location to receive flags from among #FcastFlags.
 * @dt2max: location to receive the last observation that can
 * be supported for a dynamic forecast.
 * @st2max: location to receive the last observation that can
 * be supported for a static forecast.
 *
 * Examines @pmod and determines which forecasting options are
 * applicable.
 */

void forecast_options_for_model (MODEL *pmod, const DATASET *dset,
				 FcastFlags *flags, int *dt2max,
				 int *st2max)
{
    const int *xlist;
    int *dvlags = NULL;
    int dvcheck = 0;
    int dv;

    *flags = 0;
    dv = gretl_model_get_depvar(pmod);

    if (pmod->ci == OLS) {
	if (is_standard_diff(dv, dset, NULL)) {
	    *flags |= FC_INTEGRATE_OK;
	} else {
	    *flags |= FC_MEAN_OK;
	}
    } else if (pmod->ci == NLS) {
	/* we'll try winging it! */
	if (gretl_model_get_int(pmod, "dynamic") && pmod->t2 < dset->n - 1) {
	    *flags |= FC_AUTO_OK;
	}
	return;
    }

    *dt2max = *st2max = pmod->t2;

    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	*flags |= FC_DYNAMIC_OK;
    } else if (AR_MODEL(pmod->ci)) {
	*flags |= FC_DYNAMIC_OK;
    } else if (dataset_is_time_series(dset)) {
        /* assume OK until shown otherwise */
	*flags |= FC_DYNAMIC_OK;
        dvcheck = 1;
    }

    if (*flags & FC_DYNAMIC_OK) {
	dvlags = process_lagged_depvar(pmod, dset);
        if (dvcheck && dvlags == NULL) {
            /* plain time-series, no dvlags found */
            *flags ^= FC_DYNAMIC_OK;
        }
    }

    if (addobs_can_help(pmod, dvlags, dset)) {
	*flags |= FC_ADDOBS_OK;
    }

    xlist = model_xlist(pmod);
    *dt2max = fcast_get_t2max(xlist, dvlags, pmod, dset, FC_DYNAMIC);
    *st2max = fcast_get_t2max(xlist, dvlags, pmod, dset, FC_STATIC);

    if (dvlags != NULL) {
	free(dvlags);
    }
}

static int ols_y_lag (int v, int parent_id, const DATASET *dset)
{
    int lag = 0;

    if (series_get_transform(dset, v) == LAGS) {
	int pv = series_get_parent_id(dset, v);

	if (pv == parent_id) {
	    lag = series_get_lag(dset, v);
	}
    } else {
        int pmax = dset->pd < 4 ? 4 : dset->pd;

        lag = detect_lag(parent_id, v, dset->t1, dset->t2,
                         pmax, dset);
    }

    return lag;
}

#define KSTEP_DEBUG 0

/* initialize for recursive k-step-ahead forecast */

static int k_step_init (MODEL *pmod, const DATASET *dset,
			int k, double **py, int **pllist)
{
    int *llist = NULL;
    int vy = pmod->list[1];
    int i, p;
    int nl = 0;
    int err = 0;

    llist = gretl_list_new(pmod->list[0] - 1);
    if (llist == NULL) {
        return E_ALLOC;
    }

    for (i=2; i<=pmod->list[0]; i++) {
        if (pmod->list[i] > 0) {
            p = ols_y_lag(pmod->list[i], vy, dset);
            if (p > 0) {
                /* vi is lagged dependent variable */
                llist[i-1] = p;
                nl++;
            }
        }
    }

#if KSTEP_DEBUG
    fprintf(stderr, "k_step_init: found %d y-lag terms\n", nl);
#endif

    if (nl > 0) {
        double *y = malloc(dset->n * sizeof *y);

        if (y == NULL) {
            err = E_ALLOC;
        } else {
            for (i=0; i<dset->n; i++) {
                y[i] = dset->Z[vy][i];
            }
            *py = y;
            *pllist = llist;
            llist = NULL;
        }
    }

    free(llist);

    return err;
}

static int recursive_fcast_adjust_obs (MODEL *pmod, int *t1, int t2, int k)
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

#if KSTEP_DEBUG
    fprintf(stderr, "k-step: set t1 = %d\n", *t1);
#endif

    return err;
}

/* recursive k-step ahead forecasts, for models estimated via OLS */

FITRESID *
recursive_OLS_k_step_fcast (MODEL *pmod, DATASET *dset,
			    int t1, int t2, int k,
			    int pre_n, int *err)
{
    FITRESID *fr;
    int orig_t1 = dset->t1;
    int orig_t2 = dset->t2;
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

    /* check feasibility of forecast range */
    *err = recursive_fcast_adjust_obs(pmod, &t1, t2, k);
    if (*err) {
	return NULL;
    }

    if (k > 1) {
        /* check the regressors for y-lags */
        *err = k_step_init(pmod, dset, k, &y, &llist);
        if (*err) {
            return NULL;
        }
    }

    fr = fit_resid_new_for_model(pmod, dset, t1, t2, pre_n, err);
    if (*err) {
	free(y);
	free(llist);
	return NULL;
    }

    fr->method = FC_KSTEP;
    fr->k = k;

    /* sample range for initial estimation */
    dset->t1 = pmod->t1;
    dset->t2 = t1 - k; /* start of fcast range minus k */

    /* number of forecasts */
    nf = t2 - t1 + 1;

#if KSTEP_DEBUG
    fprintf(stderr, "recursive fcast:\n"
            "  initial sample range: %d to %d\n"
            "  forecast range: %d to %d\n"
            "  steps-ahead %d, number of forecasts %d\n",
            dset->t1, dset->t2, t1, t2, k, nf);
#endif

    for (t=0; t<dset->n; t++) {
	fr->actual[t] = dset->Z[pmod->list[1]][t];
    }

    /* estimate over successive sample ranges */
    for (s=0; s<nf; s++) {
	mod = lsq(pmod->list, dset, OLS, OPT_A | OPT_Z);
	if (mod.errcode) {
	    *err = mod.errcode;
	    clear_model(&mod);
	    break;
	}
	/* the first obs following the estimation sample */
	t = dset->t2 + 1;
        /* iterate over steps */
	for (j=0; j<k; j++) {
	    yf = 0.0;
	    for (i=0; i<mod.ncoeff; i++) {
		vi = mod.list[i+2];
		p = (llist != NULL)? llist[i+1] : 0;
		if (p > 0 && p <= j) {
		    xit = y[t-p];
		} else {
		    xit = dset->Z[vi][t];
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
	dset->t2 += 1;
    }

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    if (*err) {
	free_fit_resid(fr);
	fr = NULL;
    } else {
	fit_resid_set_dec_places(fr);
	strcpy(fr->depvar, dset->varname[pmod->list[1]]);
    }

    free(y);
    free(llist);

    return fr;
}

void fcast_get_continuous_range (const FITRESID *fr, int *pt1, int *pt2)
{
    int t, t1 = fr->t1, t2 = fr->t2;

    for (t=t1; t<=t2; t++) {
	if (na(fr->actual[t]) || na(fr->fitted[t])) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=t2; t>=t1; t--) {
	if (na(fr->actual[t]) || na(fr->fitted[t])) {
	    t2--;
	} else {
	    break;
	}
    }

    *pt1 = t1;
    *pt2 = t2;
}
