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

#include "libgretl.h"
#include "version.h"
#include "libset.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"
#include "arma_priv.h"

#include "../cephes/libprob.h"

#define ARMA_DEBUG 0
#define ARMA_MDEBUG 0
#define SHOW_INIT 0

#include "arma_common.c"

static const double *as197_llt_callback (const double *b,
					 int i, void *data);

static const double *as154_llt_callback (const double *b,
					 int i, void *data);

static const double *kalman_arma_llt_callback (const double *b,
					       int i,
                                               void *data);

int maybe_correct_MA (arma_info *ainfo,
		      double *theta,
		      double *Theta)
{
    int err = 0;

    if (ainfo->q > 0) {
	err = flip_poly(theta, ainfo, 0, 0);
    }
    if (!err && ainfo->Q > 0) {
	err = flip_poly(Theta, ainfo, 0, 1);
    }

    return err;
}

/*
  Given an ARMA process $A(L)B(L) y_t = C(L)D(L) \epsilon_t$, finds the
  roots of the four polynomials -- or just two polynomials if seasonal
  AR and MA effects, B(L) and D(L) are not present -- and attaches
  this information to the ARMA model.

  pmod: MODEL pointer to which the roots info should be attached.

  ainfo: gives various pieces of information on the ARMA model,
  including seasonal and non-seasonal AR and MA orders.

  coeff: ifc + p + q + P + Q vector of coefficients (if an intercept
  is present it is element 0 and is ignored)

  returns: zero on success, non-zero on failure
*/

int arma_model_add_roots (MODEL *pmod, arma_info *ainfo,
			  const double *coeff)
{
    const double *phi =   coeff + ainfo->ifc;
    const double *Phi =     phi + ainfo->np;
    const double *theta =   Phi + ainfo->P;
    const double *Theta = theta + ainfo->nq;
    int nr = ainfo->p + ainfo->P + ainfo->q + ainfo->Q;
    int pmax, qmax, lmax;
    double *temp = NULL, *tmp2 = NULL;
    cmplx *rptr, *roots = NULL;
    int i, k, cerr = 0;

    pmax = (ainfo->p > ainfo->P)? ainfo->p : ainfo->P;
    qmax = (ainfo->q > ainfo->Q)? ainfo->q : ainfo->Q;
    lmax = (pmax > qmax)? pmax : qmax;

    if (pmax == 0 && qmax == 0) {
	return 0;
    }

    temp = malloc((lmax + 1) * sizeof *temp);
    tmp2 = malloc((lmax + 1) * sizeof *tmp2);
    roots = malloc(nr * sizeof *roots);

    if (temp == NULL || tmp2 == NULL || roots == NULL) {
	free(temp);
	free(tmp2);
	free(roots);
	return E_ALLOC;
    }

    temp[0] = 1.0;
    rptr = roots;

    if (ainfo->p > 0) {
	/* A(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		temp[i+1] = -phi[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}
	cerr = polrt(temp, tmp2, ainfo->p, rptr);
	rptr += ainfo->p;
    }

    if (!cerr && ainfo->P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo->P; i++) {
	    temp[i+1] = -Phi[i];
	}
	cerr = polrt(temp, tmp2, ainfo->P, rptr);
	rptr += ainfo->P;
    }

    if (!cerr && ainfo->q > 0) {
	/* C(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo->q; i++) {
	    if (MA_included(ainfo, i)) {
		temp[i+1] = theta[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}
	cerr = polrt(temp, tmp2, ainfo->q, rptr);
	rptr += ainfo->q;
    }

    if (!cerr && ainfo->Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo->Q; i++) {
	    temp[i+1] = Theta[i];
	}
	cerr = polrt(temp, tmp2, ainfo->Q, rptr);
    }

    free(temp);
    free(tmp2);

    if (cerr) {
	free(roots);
    } else {
	gretl_model_set_data(pmod, "roots", roots, GRETL_TYPE_CMPLX_ARRAY,
			     nr * sizeof *roots);
    }

    return 0;
}

/* add covariance matrix and standard errors based on Outer Product of
   Gradient
*/

static int arma_OPG_vcv (MODEL *pmod, void *data, int algo,
			 double *b, double s2,
			 int k, int T,
			 PRN *prn)
{
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    int err = 0;

    if (algo == 154) {
	G = numerical_score_matrix(b, T, k, as154_llt_callback,
				   data, &err);
    } else if (algo == 197) {
	G = numerical_score_matrix(b, T, k, as197_llt_callback,
				   data, &err);
    } else {
	G = numerical_score_matrix(b, T, k, kalman_arma_llt_callback,
				   data, &err);
    }

    if (!err) {
	V = gretl_matrix_XTX_new(G);
	if (V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	double rcond = gretl_symmetric_matrix_rcond(V, &err);

	if (!err && rcond < 1.0E-10) {
	    pprintf(prn, _("OPG: rcond = %g; will try Hessian\n"), rcond);
	    err = 1;
	}
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(V, s2);
	err = gretl_model_write_vcv(pmod, V);
    }

    gretl_matrix_free(G);
    gretl_matrix_free(V);

    return err;
}

static int arma_QML_vcv (MODEL *pmod, gretl_matrix *H,
			 void *data, int algo,
			 double *b, double s2, int k, int T,
			 PRN *prn)
{
    gretl_matrix *G;
    int err = 0;

    if (algo == 154) {
	G = numerical_score_matrix(b, T, k, as154_llt_callback,
				   data, &err);
    } else if (algo == 197) {
	G = numerical_score_matrix(b, T, k, as197_llt_callback,
				   data, &err);
    } else {
	G = numerical_score_matrix(b, T, k, kalman_arma_llt_callback,
				   data, &err);
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(G, sqrt(s2));
	err = gretl_model_add_QML_vcv(pmod, ARMA, H, G,
				      NULL, OPT_NONE, NULL);
    }

    gretl_matrix_free(G);

    return err;
}

static int arima_ydiff_only (arma_info *ainfo)
{
    if ((ainfo->d > 0 || ainfo->D > 0) &&
	ainfo->nexo > 0 && !arma_xdiff(ainfo)) {
	return 1;
    } else {
	return 0;
    }
}

static int arma_use_opg (gretlopt opt)
{
    int ret = 0; /* use of the Hessian is the default */

    if (opt & OPT_G) {
	ret = 1;
    } else if (libset_get_int(ARMA_VCV) == ML_OP) {
	ret = 1;
    }

    return ret;
}

static void kalman_rescale_y (gretl_vector *y, arma_info *ainfo)
{
    int i;

#if ARMA_DEBUG
    fprintf(stderr, "kalman_rescale_y: multiplying by %g\n",
	    ainfo->yscale);
#endif

    for (i=0; i<y->rows; i++) {
	if (!isnan(y->val[i])) {
	    y->val[i] -= ainfo->yshift;
	    y->val[i] *= ainfo->yscale;
	}
    }
}

/* for Kalman: convert from full-length y series to
   y vector of length ainfo->T */

static gretl_matrix *form_arma_y_vector (arma_info *ainfo,
					 int *err)
{
    gretl_matrix *yvec;

    yvec = gretl_vector_from_series(ainfo->y, ainfo->t1, ainfo->t2);

    if (yvec == NULL) {
	*err = E_ALLOC;
    } else {
	if (ainfo->yscale != 1.0) {
	    kalman_rescale_y(yvec, ainfo);
	}
#if ARMA_DEBUG
	gretl_matrix_print(yvec, "arma y vector");
#endif
    }

    return yvec;
}

static gretl_matrix *form_arma_X_matrix (arma_info *ainfo,
					 const DATASET *dset,
					 int *err)
{
    gretl_matrix *X;
    int missop;

#if ARMA_DEBUG
    printlist(ainfo->xlist, "ainfo->xlist (exog vars)");
#endif

    if (arma_na_ok(ainfo)) {
	missop = M_MISSING_OK;
    } else {
	missop = M_MISSING_ERROR;
    }

    X = gretl_matrix_data_subset(ainfo->xlist, dset,
				 ainfo->t1, ainfo->t2,
				 missop, err);

#if ARMA_DEBUG
    gretl_matrix_print(X, "X");
#endif

    return X;
}

#define y_missing(y) (na(y) || isnan(y))

/* support for AS 154 and AS 197 */

# include "as197.c"
# include "as154.c"
# include "as_driver.c"

/* support for native Kalman filter */

#include "kalman_arma.c"

static void arma_init_message (arma_info *ainfo)
{
    pprintf(ainfo->prn, "\n%s: ", _("ARMA initialization"));

    if (ainfo->init == INI_USER) {
	pprintf(ainfo->prn, "%s\n\n", _("user-specified values"));
    } else if (ainfo->init == INI_HR) {
	pprintf(ainfo->prn, "%s\n\n", _("Hannan-Rissanen method"));
    } else if (ainfo->init == INI_SMALL) {
	pprintf(ainfo->prn, "%s\n\n", _("small MA values"));
    } else if (ainfo->init == INI_NLS) {
	pprintf(ainfo->prn, "%s\n\n", _("using nonlinear AR model"));
    } else if (ainfo->init == INI_OLS) {
	pprintf(ainfo->prn, "%s\n\n", _("using linear AR model"));
    }
}

static int user_arma_init (double *coeff, arma_info *ainfo)
{
    int i, nc = n_initvals();

    if (nc == 0) {
	return 0;
    } else if (nc < ainfo->nc) {
	pprintf(ainfo->prn, _("ARMA initialization: need %d coeffs but got %d\n"),
		ainfo->nc, nc);
	return E_DATA;
    }

    if (arma_exact_ml(ainfo)) {
	/* user-specified initializer is handled within BFGSmax */
	for (i=0; i<ainfo->nc; i++) {
	    coeff[i] = 0.0;
	}
    } else {
	gretl_matrix *m = get_initvals();

	for (i=0; i<ainfo->nc; i++) {
	    coeff[i] = gretl_vector_get(m, i);
	}
	gretl_matrix_free(m);
    }

    ainfo->init = INI_USER;

    return 0;
}

/* Should we try Hannan-Rissanen initialization of ARMA
   coefficients? */

static int prefer_hr_init (arma_info *ainfo)
{
    int ret = 0;

    if (ainfo->q > 1 || ainfo->Q > 0) {
	ret = 1;
	if (arma_xdiff(ainfo)) {
	    /* don't use for ARIMAX (yet?) */
	    ret = 0;
	} else if (ainfo->T < 100) {
	    /* unlikely to work well with small sample */
	    ret = 0;
	} else if (ainfo->p > 0 && ainfo->P > 0) {
	    /* not sure about this: HR catches the MA terms, but NLS
	       handles the seasonal/non-seasonal AR interactions
	       better?
	    */
	    ret = 0;
	} else if ((ainfo->P > 0 && ainfo->p >= ainfo->pd) ||
		   (ainfo->Q > 0 && ainfo->q >= ainfo->pd)) {
	    /* overlapping seasonal/non-seasonal orders screw things up */
	    ret = 0;
	} else if (ret && arma_exact_ml(ainfo)) {
	    /* screen for cases where we'll use NLS */
	    if (ainfo->P > 0) {
		ret = 0;
	    } else if (ainfo->p + ainfo->P > 0 && ainfo->nexo > 0) {
		ret = 0;
	    } else if (ainfo->Q > 0 && arma_missvals(ainfo)) {
		ret = 0;
	    }
	}
    }

#if ARMA_DEBUG
    fprintf(stderr, "prefer_hr_init? %s\n", ret? "yes" : "no");
#endif

    return ret;
}

/* estimate an ARIMA (0,d,0) x (0,D,0) model via OLS */

static int arima_by_ls (const DATASET *dset, arma_info *ainfo,
			MODEL *pmod)
{
    gretl_matrix *X;
    gretl_matrix *b, *u, *V;
    double x, s2;
    int i, t, k = ainfo->dX->cols;
    int err = 0;

    if (ainfo->ifc) {
	/* the constant will not have been included in ainfo->dX */
	X = gretl_matrix_alloc(ainfo->T, k + 1);
	if (X == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<=k; i++) {
	    for (t=0; t<ainfo->T; t++) {
		if (i == 0) {
		    gretl_matrix_set(X, t, i, 1.0);
		} else {
		    x = gretl_matrix_get(ainfo->dX, t, i-1);
		    gretl_matrix_set(X, t, i, x);
		}
	    }
	}
	k++;
    } else {
	X = ainfo->dX;
    }

    b = gretl_column_vector_alloc(k);
    u = gretl_column_vector_alloc(ainfo->T);
    V = gretl_matrix_alloc(k, k);

    if (b == NULL || u == NULL || V == NULL) {
	err = E_ALLOC;
    } else {
	gretl_vector y;

	gretl_matrix_init(&y);
	y.rows = ainfo->T;
	y.cols = 1;
	y.val = ainfo->y + ainfo->t1;
	gretl_matrix_set_t1(&y, ainfo->t1);
	gretl_matrix_set_t2(&y, ainfo->t2);

	err = gretl_matrix_ols(&y, X, b, V, u, &s2);
    }

    if (!err) {
	pmod->ncoeff = k;
	pmod->full_n = dset->n;
	err = gretl_model_allocate_storage(pmod);
    }

    if (!err) {
	for (i=0; i<k; i++) {
	    pmod->coeff[i] = b->val[i];
	}
	for (t=0; t<ainfo->T; t++) {
	    pmod->uhat[t + ainfo->t1] = u->val[t];
	}
	err = gretl_model_write_vcv(pmod, V);
    }

    if (!err) {
    	pmod->ybar = gretl_mean(ainfo->t1, ainfo->t2, ainfo->y);
	pmod->sdy = gretl_stddev(ainfo->t1, ainfo->t2, ainfo->y);
	pmod->nobs = ainfo->T;
    }

    gretl_matrix_free(b);
    gretl_matrix_free(u);
    gretl_matrix_free(V);

    if (X != ainfo->dX) {
	gretl_matrix_free(X);
    }

    return err;
}

/* calculate info criteria for compatibility with ML? */
#define ML_COMPAT 1 /* 2017-03-23 */

static int arma_via_OLS (arma_info *ainfo, const double *coeff,
			 const DATASET *dset, MODEL *pmod)
{
    int err = 0;

    ainfo->flags |= ARMA_LS;

    if (arma_xdiff(ainfo)) {
	err = arima_by_ls(dset, ainfo, pmod);
    } else {
	err = arma_by_ls(coeff, dset, ainfo, pmod);
    }

    if (!err) {
	ArmaFlags f = arma_exact_ml(ainfo) ? ARMA_OLS : ARMA_LS;

	pmod->t1 = ainfo->t1;
	pmod->t2 = ainfo->t2;
	pmod->full_n = dset->n;
	write_arma_model_stats(pmod, ainfo, dset);
	if (arma_exact_ml(ainfo)) {
#if ML_COMPAT
	    /* In the case of ainfo->nc == 0 (no coefficients
	       actually estimated), pmod->ncoeff will be 1, since
	       we add a dummy constant with value 0. That "1"
	       will account for the variance estimate so that
	       addk should be zero. Otherwise we add 1 for the
	       variance estimate.
	    */
	    int addk = ainfo->nc == 0 ? 0 : 1;

	    mle_criteria(pmod, addk);
#else
	    ls_criteria(pmod);
#endif
	} else {
	    arma_model_add_roots(pmod, ainfo, pmod->coeff);
	}
	gretl_model_set_int(pmod, "arma_flags", f);
    }

    if (!err && pmod->errcode) {
	err = pmod->errcode;
    }

    return err;
}

/* Set flag to indicate differencing of exogenous regressors, in the
   case of an ARIMAX model using native exact ML -- unless this is
   forbidden by OPT_Y (--y-diff-only).  Note that we don't do this
   when we're using conditional ML (BHHH).
*/

static void maybe_set_xdiff_flag (arma_info *ainfo, gretlopt opt)
{
    if (arma_exact_ml(ainfo) &&
	(ainfo->d > 0 || ainfo->D > 0) &&
	ainfo->nexo > 0 && !(opt & OPT_Y)) {
	ainfo->pflags |= ARMA_XDIFF;
    }
}

/* Respond to OPT_S (--stdx): standardize exogenous regressors */

static int arma_standardize_x (arma_info *ainfo,
			       DATASET *dset)
{
    int orig_v = dset->v;
    int err = 0;

    ainfo->xstats = gretl_matrix_alloc(ainfo->nexo, 2);
    if (ainfo->xstats == NULL) {
	return E_ALLOC;
    }

    err = dataset_add_series(dset, ainfo->nexo);

    if (!err) {
	double xbar, sdx;
	int i, vi, vj, t;

	for (i=0; i<ainfo->nexo && !err; i++) {
	    vi = ainfo->xlist[i+1];
	    err = gretl_moments(ainfo->t1, ainfo->t2, dset->Z[vi],
				NULL, &xbar, &sdx, NULL, NULL, 1);
	    if (!err) {
		vj = orig_v + i;
		for (t=0; t<dset->n; t++) {
		    dset->Z[vj][t] = (dset->Z[vi][t] - xbar) / sdx;
		}
		/* replace x-ref with standardized version */
		ainfo->xlist[i+1] = vj;
		/* and record the stats used */
		gretl_matrix_set(ainfo->xstats, i, 0, xbar);
		gretl_matrix_set(ainfo->xstats, i, 1, sdx);
	    }
	}
    }

    if (!err) {
	set_arma_stdx(ainfo);
    }

    return err;
}

/* Set flag to allow NAs within the sample range for an
   ARMA model using native exact ML, or for one that can
   be estimated via OLS (with no MA term).
*/

static void maybe_allow_missvals (arma_info *ainfo)
{
    if (arma_exact_ml(ainfo)) {
	ainfo->pflags |= ARMA_NAOK;
    } else if (ainfo->q == 0 && ainfo->Q == 0) {
	ainfo->pflags |= ARMA_NAOK;
    }
}

static int check_arma_options (gretlopt opt)
{
    int err;

    /* can't specify LBFGS or --robust with conditional ML */
    err = options_incompatible_with(opt, OPT_C, OPT_L | OPT_R);

    if (!err) {
	/* nor more than one of AS 154, CML, Kalman */
	err = incompatible_options(opt, OPT_A | OPT_C | OPT_K);
    }

    return err;
}

typedef struct dim_info_ dim_info;

struct dim_info_ {
    int p; /* max non-seasonal AR order */
    int d; /* non-seasonal difference */
    int q; /* max non-seasonal MA order */
    int P; /* seasonal AR order */
    int D; /* seasonal difference */
    int Q; /* seasonal MA order */
};

static int arma_precheck (const int *list,
			  const int *pqspec,
			  gretlopt opt,
			  DATASET *dset,
			  arma_info *painfo,
			  dim_info *dinfo)
{
    arma_info ainfo = {0};
    int err;

    arma_info_init(&ainfo, opt, pqspec, dset);
    err = check_arma_options(opt);
    if (!err) {
	ainfo.alist = gretl_list_copy(list);
	if (ainfo.alist == NULL) {
	    err = E_ALLOC;
	}
    }
    if (!err) {
	err = arma_check_list(&ainfo, dset, opt);
    }
    if (!err) {
	if (painfo != NULL) {
	    *painfo = ainfo;
	} else if (dinfo != NULL) {
	    /* just transcribe the basics */
	    dinfo->p = ainfo.p;
	    dinfo->d = ainfo.d;
	    dinfo->q = ainfo.q;
	    dinfo->P = ainfo.P;
	    dinfo->D = ainfo.D;
	    dinfo->Q = ainfo.Q;
	}
	if (painfo == NULL) {
	    free(ainfo.alist);
	}
    }

    return err;
}

MODEL arma_model (const int *list, const int *pqspec,
		  DATASET *dset, gretlopt opt, PRN *prn)
{
    double *coeff = NULL;
    MODEL armod;
    arma_info ainfo_s, *ainfo;
    int missv = 0, misst = 0;
    int orig_v = dset->v;
    int err;

    gretl_model_init(&armod, dset);

    err = arma_precheck(list, pqspec, opt, dset, &ainfo_s, NULL);
    if (err) {
	armod.errcode = err;
	return armod;
    }

    ainfo = &ainfo_s;
    if (opt & OPT_V) {
	ainfo->prn = prn;
    }

    if (!err) {
	/* calculate maximum lag */
	maybe_set_xdiff_flag(ainfo, opt);
	calc_max_lag(ainfo);
    }

    if (!err) {
	/* adjust sample range if need be */
	maybe_allow_missvals(ainfo);
	err = arma_adjust_sample(ainfo, dset, &missv, &misst);
	if (err) {
	    if (missv > 0 && misst > 0) {
		gretl_errmsg_sprintf(_("Missing value encountered for "
				       "variable %d, obs %d"), missv, misst);
	    }
	} else if (missv > 0) {
	    set_arma_missvals(ainfo);
	}
    }

    if (!err && ainfo->nexo > 0 && arma_exact_ml(ainfo)) {
	/* FIXME check conditionality more rigorously */
	if ((opt & OPT_S) && !arma_xdiff(ainfo)) {
	    /* --stdx */
	    err = arma_standardize_x(ainfo, dset);
	}
    }

    if (!err) {
	/* allocate initial coefficient vector */
	coeff = malloc(ainfo->nc * sizeof *coeff);
	if (coeff == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* organize the dependent variable */
	ainfo->y = (double *) dset->Z[ainfo->yno];
	if (ainfo->d > 0 || ainfo->D > 0) {
	    if (arma_missvals(ainfo)) {
		/* for now: insist on native Kalman, since only it
		   handles the levels formulation of ARIMA
		*/
		opt &= ~OPT_A;
		opt |= OPT_K;
		set_arima_levels(ainfo);
	    } else {
		/* note: this replaces ainfo->y */
		err = arima_difference(ainfo, dset, 0);
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    if (ainfo->p == 0 && ainfo->P == 0 &&
	ainfo->q == 0 && ainfo->Q == 0 &&
	arma_exact_ml(ainfo) && !arma_missvals(ainfo)) {
	/* pure "I" model, no NAs: OLS provides the MLE */
	err = arma_via_OLS(ainfo, NULL, dset, &armod);
	goto bailout; /* estimation handled */
    }

    /* start initialization of the coefficients */

    /* first see if the user specified some values */
    err = user_arma_init(coeff, ainfo);
    if (err) {
	goto bailout;
    }

    if (!arma_exact_ml(ainfo) && ainfo->q == 0 && ainfo->Q == 0) {
	/* for a pure AR model, the conditional MLE is least
	   squares (OLS or NLS); in the NLS case a user-specified
	   initializer may be useful, if present
	*/
	const double *b = ainfo->init ? coeff : NULL;

	err = arma_via_OLS(ainfo, b, dset, &armod);
	goto bailout; /* estimation handled */
    }

    /* see if it may be helpful to scale the dependent
       variable, if we're doing exact ML */
    if (arma_exact_ml(ainfo) && ainfo->ifc && ainfo->init != INI_USER) {
	maybe_set_yscale(ainfo);
#if SHOW_INIT
	fprintf(stderr, "yscale = %g\n", ainfo->yscale);
#endif
    }

    /* try Hannan-Rissanen init, if suitable */
    if (!ainfo->init && prefer_hr_init(ainfo)) {
	hr_arma_init(coeff, dset, ainfo);
#if SHOW_INIT
	fprintf(stderr, "HR init (%d %d): %s\n", ainfo->p, ainfo->q,
		ainfo->init ? "success" : "fail");
#endif
    }

    /* initialize via AR model by OLS or NLS, adding minimal
       MA coefficients if needed: this is the fallback if
       Hannan-Rissanen fails, but also the default if the
       conditions of applicability of H-R are not met
    */
    if (!err && !ainfo->init) {
	err = ar_arma_init(coeff, dset, ainfo, &armod, opt);
#if SHOW_INIT
	fprintf(stderr, "AR init: err = %d\n", err);
#endif
    }

    if (ainfo->prn != NULL && ainfo->init) {
	arma_init_message(ainfo);
    }

    if (!err) {
	clear_model_xpx(&armod);
	if (arma_exact_ml(ainfo)) {
	    if (opt & OPT_K) {
		err = kalman_arma(coeff, dset, ainfo, &armod, opt);
	    } else {
		err = as_arma(coeff, dset, ainfo, &armod, opt);
	    }
	} else {
	    err = bhhh_arma(coeff, dset, ainfo, &armod, opt);
	}
    }

 bailout:

    if (err && !armod.errcode) {
	armod.errcode = err;
    }

    if (!err) {
	transcribe_extra_info(ainfo, &armod);
    }

    if (!armod.errcode) {
	gretl_model_smpl_init(&armod, dset);
    }

    free(coeff);
    arma_info_cleanup(ainfo);

    if (dset->v > orig_v) {
	dataset_drop_last_variables(dset, dset->v - orig_v);
    }

    return armod;
}

static gretl_matrix *arma_select_matrix (int pmax, int qmax)
{
    gretl_matrix *m;

    m = gretl_matrix_alloc((pmax+1) * (qmax+1), 5);

    if (m != NULL) {
        char **S = strings_array_new(5);

        if (S != NULL) {
            S[0] = gretl_strdup("p");
            S[1] = gretl_strdup("q");
            S[2] = gretl_strdup("AIC");
            S[3] = gretl_strdup("BIC");
            S[4] = gretl_strdup("HQC");
            gretl_matrix_set_colnames(m, S);
        }
    }

    return m;
}

static void arma_select_header (dim_info *dinfo, int w, PRN *prn)
{
    char word[8], s1[32], s2[32] = {0};

    if (dinfo->d) {
	sprintf(s1, "(p %d q)", dinfo->d);
	strcpy(word, "ARIMA");
    } else {
	strcpy(s1, "(p q)");
	strcpy(word, "ARMA");
    }
    if (dinfo->D) {
	sprintf(s2, "(%d %d %d)", dinfo->P, dinfo->D, dinfo->Q);
	strcpy(word, "SARIMA");
    } else if (dinfo->P || dinfo->Q) {
	sprintf(s2, "(%d %d)", dinfo->P, dinfo->Q);
	strcpy(word, "SARMA");
    }

    pprintf(prn, "\nInformation Criteria for %s%s%s specifications\n",
	    word, s1, s2);
    pputs(prn, "-----------------------------------------------\n");
    pprintf(prn, " p, q %*s %*s %*s\n", w, "AIC", w+1, "BIC", w+1, "HQC");
    pputs(prn, "-----------------------------------------------\n");
}

/* A simple start on providing built-in means of selecting the AR and
   MA orders of an AR(I)MA model via Information Criteria.  As things
   stand we accept a general (p d q)(P D Q) specification but we only
   actually search over p and q, the other paramaters being clamped at
   their input values.
*/

int arma_select (const int *list, const int *pqspec,
                 DATASET *dset, gretlopt opt, PRN *prn)
{
    dim_info dinfo = {0};
    MODEL amod;
    int print = !(opt & OPT_Q);
    gretlopt aopt = opt | OPT_Q;
    gretl_matrix *m = NULL;
    double cij, mincrit[3];
    int minidx[3] = {0};
    int *mylist = NULL;
    int p, q, pmax, qmax;
    int qpos = 2;
    int i, j;
    int err = 0;

    if (pqspec != NULL) {
        err = E_DATA;
    } else {
	err = arma_precheck(list, NULL, opt, dset, NULL, &dinfo);
    }
    if (err) {
	return err;
    }

    mylist = gretl_list_copy(list);
    if (gretl_list_separator_position(mylist) == 4) {
	qpos = 3;
    }
    pmax = mylist[1];
    qmax = mylist[qpos];
    if (pmax < 0 || qmax < 0) {
	err = E_INVARG;
    }
    if (!err) {
        m = arma_select_matrix(pmax, qmax);
        if (m == NULL) {
            err = E_ALLOC;
        }
    }
    if (err) {
	free(mylist);
        return err;
    }

    /* fill the results matrix */
    i = 0;
    for (p=0; p<=pmax; p++) {
        for (q=0; q<=qmax; q++) {
            gretl_matrix_set(m, i, 0, p);
            gretl_matrix_set(m, i, 1, q);
            mylist[1] = p;
	    mylist[qpos] = q;
            amod = arma_model(mylist, NULL, dset, aopt, NULL);
            if (amod.errcode) {
                for (j=0; j<3; j++) {
                    gretl_matrix_set(m, i, j+2, NADBL);
                }
            } else {
                for (j=0; j<3; j++) {
                    cij = amod.criterion[j];
                    gretl_matrix_set(m, i, j+2, cij);
                    if (i == 0) {
                        mincrit[j] = cij;
                    } else if (cij < mincrit[j]) {
                        mincrit[j] = cij;
                        minidx[j] = i;
                    }
                }
            }
            clear_model(&amod);
            i++;
        }
    }

    if (print) {
        /* print as table */
        int width = 12, ndec = 4;

        arma_select_header(&dinfo, width, prn);
        i = 0;
        for (p=0; p<=pmax; p++) {
            for (q=0; q<=qmax; q++) {
                pprintf(prn, "%2d,%2d", p, q);
                for (j=0; j<3; j++) {
                    cij = gretl_matrix_get(m, i, j+2);
                    if (na(cij)) {
                        pprintf(prn, " %*s", width, "NA");
                    } else {
                        pprintf(prn, " %*.*f%", width, ndec, cij);
                        pputc(prn, i == minidx[j] ? '*' : ' ');
                    }
                }
                pputc(prn, '\n');
                i++;
            }
        }
    }

    free(mylist);
    if (m != NULL) {
        record_matrix_test_result(m, NULL);
    }

    return err;
}
