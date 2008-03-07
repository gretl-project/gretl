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
#include "libset.h"
#include "matrix_extra.h"
#include "bhhh_max.h"

#define BHHH_DEBUG 0

struct _model_info {

    /* members that may be set by caller of bhhh_max, via 
       initialization or accessor functions */

    int k;              /* number of parameters */
    int t1, t2;         /* starting and ending point of sample */
    int bign;           /* total observations in current real dataset */
    int n_series;       /* number of additional series needed in the
                           likelihood and/or score calculations */
    double tol;         /* tolerance for convergence */
    unsigned char opts; /* options from among BHHH_opts */

    void *extra_info;   /* pointer to additional info which may be
			   used in likelihood callback */

    /* members set within bhhh_max */

    int n;            /* length of series */
    int iters;        /* number of iterations taken */
    double ll;        /* log-likelihood */
    double ll2;       /* test log-likelihood value */
    double *theta;    /* vector of parameters */
    double **series;  /* additional series */

    /* full VCV matrix from OPG regression, if (opts & FULL_VCV_MATRIX) */
    gretl_matrix *VCV; 

    /* pointer to OPG model, if (opts & PRESERVE_OPG_MODEL) */
    MODEL *pmod;
};

/**
 * model_info_free:
 * @minfo: model info pointer.
 *
 * Frees the dynamically allocated members of @minfo, then
 * frees @minfo itself.
 */

void model_info_free (model_info *minfo)
{
    int i;

    if (minfo == NULL) return;

    free(minfo->theta);

    if (minfo->series != NULL) {
	for (i=0; i<minfo->n_series; i++) {
	    free(minfo->series[i]);
	}
	free(minfo->series);
    }

    if (minfo->VCV != NULL) {
	gretl_matrix_free(minfo->VCV);
    }

    free(minfo);
}

/**
 * model_info_capture_OPG_model:
 * @minfo: model info pointer.
 *
 * Returns: pointer to a gretl #MODEL, which contains the
 * results of the OPG (Outer Product of the Gradient) regression
 * associated with @minfo, if the preservation of this model
 * has been flagged by use of the option %PRESERVE_OPG_MODEL:
 * see see model_info_set_opts(). Otherwise returns %NULL.
 */

MODEL *model_info_capture_OPG_model (model_info *minfo)
{
    MODEL *pmod = minfo->pmod;

    minfo->pmod = NULL;
    return pmod;
}

/**
 * model_info_get_VCV:
 * @minfo: model info pointer.
 *
 * Returns: pointer to the covariance matrix of @minfo, if
 * the creation of this matrix has been flagged by use of the
 * option %FULL_VCV_MATRIX: see model_info_set_opts().
 * Otherwise returns %NULL.
 */

gretl_matrix *model_info_get_VCV (model_info *minfo)
{
    return minfo->VCV;
}

/**
 * model_info_get_theta:
 * @minfo: model info pointer.
 *
 * Returns: the parameter vector @theta from @minfo. 
 * This array is freed when @minfo is freed, so if the
 * results are wanted by the caller, they should be
 * copied.
 */

double *model_info_get_theta (model_info *minfo)
{
    return minfo->theta;
}

/**
 * model_info_get_t1:
 * @minfo: model info pointer.
 *
 * Returns: the (zero-based) start of the sample range for 
 * @minfo.
 */

int model_info_get_t1 (const model_info *minfo)
{
    return minfo->t1;
}

/**
 * model_info_get_t2:
 * @minfo: model info pointer.
 *
 * Returns: the end of the sample range for @minfo.
 */

int model_info_get_t2 (const model_info *minfo)
{
    return minfo->t2;
}

/**
 * model_info_get_n:
 * @minfo: model info pointer.
 *
 * Returns: the length of the data series used in @minfo.
 */

int model_info_get_n (const model_info *minfo)
{
    return minfo->n;
} 

/**
 * model_info_get_iters:
 * @minfo: model info pointer.
 *
 * Returns: the number of iterations taken in estimating
 * the model in @minfo.
 */  

int model_info_get_iters (const model_info *minfo)
{
    return minfo->iters;
}

/**
 * model_info_get_series:
 * @minfo: model info pointer.
 *
 * Returns: the allocated two-dimensional array set on
 * @minfo using model_info_set_n_series().  This array
 * is freed when @minfo is freed, so the returned pointer
 * itself should not be modified.  
 */  

double **model_info_get_series (const model_info *minfo)
{
    return minfo->series;
}

/**
 * model_info_get_ll:
 * @minfo: model info pointer.
 *
 * Returns: the log-likelihood for @minfo.
 */ 

double model_info_get_ll (const model_info *minfo)
{
    return minfo->ll;
}

/**
 * model_info_set_ll:
 * @minfo: model info pointer.
 * @ll: log-likehood.
 * @do_score:
 *
 * Sets the log-likelihood for @minfo.  If @do_score is non-zero,
 * sets the primary ll value, otherwise sets the secondary value,
 * which is used for comparison (to see if we have succeeded in
 * increasing the likelihood by following the estimated gradient).
 */ 

void model_info_set_ll (model_info *minfo, double ll, int do_score)
{
    if (do_score) {
	minfo->ll = ll;
    } else {
	minfo->ll2 = ll;
    }
}

/**
 * model_info_set_opts:
 * @minfo: model info pointer.
 * @opts: option flags to set, from #BHHH_opts.
 *
 * Sets the option flags for @minfo.  If the option %PRESERVE_OPG_MODEL
 * is set, a pointer to the #MODEL struct used internally by bhhh_max() 
 * is available on exit via model_info_capture_OPG_model(); otherwise 
 * this is freed on exit.  If %FULL_VCV_MATRIX is set, a pointer to
 * the full covariance matrix from the last iteration of the OPG
 * model is available on exit using model_info_get_VCV().
 */ 

void model_info_set_opts (model_info *minfo, unsigned char opts)
{
    minfo->opts = opts;
}

/**
 * model_info_get_k:
 * @minfo: model info pointer.
 *
 * Returns: the number of regressors in @minfo.
 */ 

int model_info_get_k (model_info *minfo)
{
    return minfo->k;
}

/**
 * model_info_set_n_series:
 * @minfo: model info pointer.
 * @n: number of auxiliary data series.
 *
 * Sets the number of auxiliary data series needed for the
 * calculations associated with @minfo.  These series are
 * allocated within the bhhh_max() function, and are freed
 * when @minfo is freed; see model_info_free().  The length
 * of each series is given by the %n member of @minfo, which
 * is set implicitly when the start and end of the sample
 * range are set, using model_info_new(), and which
 * can be retrieved using model_info_get_n().  This length
 * is in fact t2 (zero-based ending observation of sample range) 
 * plus one.  
 * 
 * These series (if any) are not actually used within bhhh_max() 
 * itself: they are intended for use within the log-likelihood 
 * callback (serving as persistent storage, and saving the 
 * callback function from having to allocate and deallocate 
 * storage on each call).  The series may be written to in that 
 * context, but the pointers themselves should not be altered 
 * in any way.
 */ 

void model_info_set_n_series (model_info *minfo, int n)
{
    minfo->n_series = n;
}

/**
 * model_info_set_extra_info:
 * @minfo: model info pointer.
 * @extra: pointer to set on @minfo.
 *
 * Set the content of the "extra" pointer member of @minfo.
 * This pointer is not used within bhhh_max(), it is intended
 * for use within the log-likelihood callback.  Setting this
 * pointer is therefore optional.
 */

void model_info_set_extra_info (model_info *minfo, void *extra)
{
    minfo->extra_info = extra;
}

/**
 * model_info_get_extra_info:
 * @minfo: model info pointer.
 *
 * Retrieves the content of the "extra" pointer member of @minfo.
 *
 * Returns: the pointer that was set with model_info_set_extra_info(),
 * or %NULL if none was set.
 */

void *model_info_get_extra_info (model_info *minfo)
{
    return minfo->extra_info;
}

/* Construct the regression list for the OPG regression, with the
   appropriate indices into the temporary artificial dataset.
*/

static int *make_opg_list (int k)
{
    int *list;
    int i;

    list = gretl_list_new(k + 1);
    if (list == NULL) return NULL;

    list[1] = 0;  /* dep var is the constant */
    for (i=0; i<k; i++) {
	list[i+2] = i + 1; 
    }

#if BHHH_DEBUG
    printlist(list, "OPG regression list");
#endif

    return list;
}

static int model_info_init (model_info *minfo, const double *init_coeff)
{
    int i, t, err = 0;
    int n_series = minfo->n_series;

    minfo->theta = malloc(minfo->k * sizeof *minfo->theta);

    if (minfo->theta == NULL) {
	model_info_free(minfo);
	return 1;
    }    

    /* allocate series */
    if (n_series > 0) {
	minfo->series = malloc(n_series * sizeof *minfo->series);
	if (minfo->series == NULL) {
	    model_info_free(minfo);
	    return 1;
	}
	minfo->n_series = 0;
	for (i=0; i<n_series; i++) {
	    minfo->series[i] = malloc(minfo->n * sizeof **minfo->series);
	    if (minfo->series[i] == NULL) {
		model_info_free(minfo);
		return 1;
	    }
	    for (t=0; t<minfo->n; t++) {
		minfo->series[i][t] = 0.0;
	    }
	    minfo->n_series += 1;
	}
    }

    /* initialize parameters */
    for (i=0; i<minfo->k; i++) {
	minfo->theta[i] = init_coeff[i];
    }

    return err;
}

/**
 * model_info_new:
 * @k: number of parameters to be estimated.
 * @t1: starting observation of estimation range.
 * @t2: ending observation of estimation range.
 * @bign: total obserations in current dataset.
 * @tol: tolerance for assessing convergence.
 *
 * Creates the basic information structure required by
 * bhhh_max().  Within that function, iteration is terminated 
 * when the difference in the log-likelihood from one round to 
 * the next falls below @tol.  Note that @t1 and @t2 are
 * zero-based.  
 *
 * Returns: pointer to newly allocated model_info, or %NULL 
 * on failure.
 */
 
model_info *model_info_new (int k, int t1, int t2, int bign, double tol)
{
    model_info *mi;

    mi = malloc(sizeof *mi);

    if (mi != NULL) {
	mi->theta = NULL;
	mi->series = NULL;
	mi->pmod = NULL;
	mi->VCV = NULL;
	mi->n_series = 0;

	mi->k = k;
	mi->t1 = t1;
	mi->t2 = t2;
	mi->n = t2 + 1;
	mi->bign = bign;
	mi->tol = tol;
    }

    return mi;
}

static int expand_model_series (MODEL *pmod, model_info *minfo)
{
    double *x;
    int t;

    x = realloc(pmod->uhat, minfo->bign * sizeof *x);
    if (x == NULL) {
	pmod->errcode = E_ALLOC;
    } else {
	pmod->uhat = x;
	for (t=minfo->t2; t<minfo->bign; t++) {
	    pmod->uhat[t] = NADBL;
	}
    }

    x = realloc(pmod->yhat, minfo->bign * sizeof *x);
    if (x == NULL) {
	pmod->errcode = E_ALLOC;
    } else {
	pmod->yhat = x;
	for (t=minfo->t2; t<minfo->bign; t++) {
	    pmod->yhat[t] = NADBL;
	}
    }

    if (pmod->errcode == 0) {
	pmod->full_n = minfo->bign;
    }

    return 0;
}

/**
 * bhhh_max:
 * @loglik: pointer to function for calculating log-likelihood and
 * score matrix; see below.
 * @X: primary data set (not used within this function, but passed 
 * on for use by @loglik).
 * @init_coeff: starting values for coefficients.
 * @minfo: model info struct; see below.
 * @opt: can include %OPT_V for verbose output.
 * @prn: printing struct for iteration info (or %NULL).
 *
 * Maximize likelihood using the BHHH conditional ML method,
 * implemented via iteration of the Outer Product of the Gradient 
 * (OPG) regression with line search.
 *
 * The @minfo pointer is obtained using model_info_new().
 * Optional settings may be applied using model_info_set_opts(), 
 * model_info_set_n_series() and model_info_set_extra_info().
 *
 * @loglik is called to calculate the log-likelihood for the model
 * in question.  The parameters passed to this function are:
 * (1) the array of estimated coefficients; (2) the primary data 
 * array @X; (3) a further data array %Z; (4) @minfo; and (5) an 
 * integer that is 1 if the score matrix should be calculated in
 * %Z, otherwise 0.  The %Z array is allocated within %bhhh_max
 * and is freed on exit; it is a two-dimensional array with %k + 1
 * members, each of length %n (where %k and %n may be obtained
 * using model_info_get_k() and model_info_get_n()).  The first
 * member, %Z[0], represents the constant or intercept, and should
 * not be modified.  
 *
 * For an example of the use of this function, see arma.c in the
 * %plugin directory of the gretl source.
 * 
 * Returns: 0 on successful completion, non-zero error code otherwise.
 */

int bhhh_max (LL_FUNC loglik,
	      const double **X, 
	      const double *init_coeff,
	      model_info *minfo, 
	      gretlopt opt,
	      PRN *prn)
{
    /* OPG model */
    MODEL *bmod;
    int *blist;

    /* temporary artificial dataset */
    double **tZ = NULL;
    DATAINFO *tinfo = NULL;

    int iters, itermax;
    double minstep = 1.0e-06;
    double crit = 1.0;
    double stepsize = 0.25;

    double *delta = NULL, *ctemp = NULL;

    int i, t, err, k;

    err = model_info_init(minfo, init_coeff);
    if (err) {
	return E_ALLOC;
    }

    k = minfo->k;
    blist = make_opg_list(k);
    if (blist == NULL) {
	return E_ALLOC;
    }

    delta = malloc(k * sizeof *delta);
    ctemp = malloc(k * sizeof *ctemp);
    if (delta == NULL || ctemp == NULL) {
	free(delta);
	free(ctemp);
	free(blist);
	return E_ALLOC;
    }

    /* create temp dataset for OPG regression: k vars plus constant */    
    tinfo = create_auxiliary_dataset(&tZ, k + 1, minfo->n);
    if (tinfo == NULL) {
	free(delta);
	free(ctemp);
	free(blist);
	return E_ALLOC;
    } 

    /* respect the incoming sample range */
    tinfo->t1 = minfo->t1;
    tinfo->t2 = minfo->t2;

    /* zero the dataset */
    for (i=1; i<=k; i++) {
#if BHHH_DEBUG
	sprintf(tinfo->varname[i], "tZ%d", i);
#endif
	for (t=0; t<minfo->n; t++) {
	    tZ[i][t] = 0.0;
	}
    }

    /* initialize OPG model */
    bmod = gretl_model_new();

    iters = 0;
    itermax = libset_get_int(BHHH_MAXITER);

    while (crit > minfo->tol && iters++ < itermax) {

	/* compute loglikelihood and score matrix */
	err = loglik(minfo->theta, X, tZ, minfo, 1); 
	if (err) {
	    pputs(prn, "Error calculating log-likelihood\n");
	    err = E_NOCONV;
	    break;
	}

#if BHHH_DEBUG
	pprintf(prn, "Top of loop: ll = %g\n", minfo->ll);
	pprintf(prn, "Dataset for OPG regression:\n");
	printdata(blist, (const double **) tZ, tinfo, OPT_O, prn);
#endif

	/* BHHH via OPG regression */
	*bmod = lsq(blist, &tZ, tinfo, OLS, OPT_A);
	if (bmod->errcode) {
	    fprintf(stderr, "BHHH model error code = %d\n", bmod->errcode);
	    err = E_NOCONV;
	    break;
	}

#if BHHH_DEBUG
	pprintf(prn, "OLS model, OPG regression:\n");
	printmodel(bmod, tinfo, 0, prn);
#endif
	for (i=0; i<k; i++) {
	    delta[i] = bmod->coeff[i] * stepsize;
	    ctemp[i] = minfo->theta[i] + delta[i];
	} 
	
	clear_model(bmod);

	/* see if we've gone up... (0 means "don't compute score") */
	err = loglik(ctemp, X, tZ, minfo, 0); 

#if BHHH_DEBUG
	pprintf(prn, "bhhh loop: initial ll2 = %g\n", minfo->ll2);
#endif

	while (minfo->ll2 < minfo->ll || err) { 
	    /* ... if not, halve the steplength */
	    stepsize *= 0.5;
	    if (stepsize < minstep) {
		err = E_NOCONV;
		break;
	    }
	    for (i=0; i<k; i++) {
		delta[i] *= 0.5;
		ctemp[i] = minfo->theta[i] + delta[i];
	    }
	    err = loglik(ctemp, X, tZ, minfo, 0);
#if BHHH_DEBUG
	    pprintf(prn, "bhhh loop: modified ll2 = %g\n", minfo->ll2);
#endif
	}

	if (err) break;

	/* actually update parameter estimates */
	for (i=0; i<k; i++) {
	    minfo->theta[i] = ctemp[i];
	}	

	/* double the steplength? (was < 4.0 below) */
	if (stepsize < 1.0) {
	    stepsize *= 2.0;
	}

	/* print interation info, if wanted */
	if (opt & OPT_V) {
	    print_iter_info(iters, minfo->ll, C_LOGLIK, k, minfo->theta, delta, 
			    stepsize, prn);
	}

	crit = minfo->ll2 - minfo->ll;  
    }

    if (opt & OPT_V) {
	print_iter_info(-1, minfo->ll, C_LOGLIK, k, minfo->theta, delta, 
			stepsize, prn);
    }

    if (crit > minfo->tol || err != 0) {
	fprintf(stderr, "bhhh_max: iters = %d, crit = %g, tol = %g, err = %d\n",
		iters, crit, minfo->tol, err);
	err = E_NOCONV;
    }

    free(delta);
    free(ctemp);

    if (!err) {
	if (minfo->opts & FULL_VCV_MATRIX) {
	    gretl_matrix *G, *VCV;

	    G = gretl_matrix_from_2d_array((const double **) tZ + 1, 
					   minfo->n, minfo->k);
	    VCV = gretl_matrix_vcv(G);
	    minfo->VCV = VCV;
	    gretl_matrix_free(G);
	}
	if (minfo->opts & PRESERVE_OPG_MODEL) {
	    int qr_bak = libset_get_bool(USE_QR);

	    /* run OPG once more using QR, to get packed VCV */
	    libset_set_bool(USE_QR, 1);
	    *bmod = lsq(blist, &tZ, tinfo, OLS, OPT_A);
	    libset_set_bool(USE_QR, qr_bak);
	    minfo->pmod = bmod;
	    gretl_model_set_int(bmod, "iters", iters);
	    if (minfo->bign > minfo->t2) {
		expand_model_series(bmod, minfo);
	    }
	} 
	minfo->iters = iters;
    }

    /* free all remaining temp stuff */
    destroy_dataset(tZ, tinfo);
    free(blist);

    if (bmod != minfo->pmod) {
	free(bmod);
    }

    return err;
}
