/* 
 * Copyright (C) 2004 Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
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
#include "libset.h"
#include "bhhh_max.h"

#undef BHHH_DEBUG

struct _model_info {

    /* members that may be set by caller of bhhh_max, via accessor
       functions */

    int k;              /* number of parameters */
    int t1, t2;         /* starting and ending point of sample */
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

#define DEFAULT_MAX_ITER 1000

static int get_maxiter (void)
{
    char *mistr = getenv("GRETL_MAX_ITER");
    int mi = DEFAULT_MAX_ITER;

    if (mistr != NULL) {
	if (!sscanf(mistr, "%d", &mi)) {
	    mi = DEFAULT_MAX_ITER;
	}
    }

    return mi;
}

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
 * associated with @minfo.
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
 * Returns: pointer to the covariance matrix of @minfo.
 */

gretl_matrix *model_info_get_VCV (model_info *minfo)
{
    return minfo->VCV;
}

/**
 * model_info_get_theta:
 * @minfo: model info pointer.
 *
 * Returns: parameter vector @theta from @minfo.
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
 * Returns: the number of observations used in @minfo.
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
 * Returns: FIXME.
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
 * Sets the option flags for @minfo. 
 */ 

void model_info_set_opts (model_info *minfo, unsigned char opts)
{
    minfo->opts = opts;
}

/**
 * model_info_set_tol:
 * @minfo: model info pointer.
 * @tol: tolerance for convergence of estimates.
 *
 * Sets the convergence tolerance for @minfo. 
 */ 

void model_info_set_tol (model_info *minfo, double tol)
{
    minfo->tol = tol;
}

/**
 * model_info_set_k:
 * @minfo: model info pointer.
 * @k: number of regressors in minfo.
 *
 * Sets the total number of regressors.
 */ 

void model_info_set_k (model_info *minfo, int k)
{
    minfo->k = k;
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
 * @n: number of data series.
 *
 * Sets the number of auxiliary data series needed for the
 * OPG regression associated with @minfo.
 */ 

void model_info_set_n_series (model_info *minfo, int n)
{
    minfo->n_series = n;
}

/**
 * model_info_set_t1_t2:
 * @minfo: model info pointer.
 * @t1: starting observation number (zero-based).
 * @t2: ending observation number.
 *
 * Sets the sample range for @minfo.
 */ 

void model_info_set_t1_t2 (model_info *minfo, int t1, int t2)
{
    minfo->t1 = t1;
    minfo->t2 = t2;
    minfo->n = t2 + 1;
}

/**
 * model_info_set_extra_info:
 * @minfo: model info pointer.
 * @extra: pointer to set on @minfo.
 *
 * Set the content of the "extra" pointer member of @minfo.
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
 * or NULL if none was set.
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

    list = malloc((k + 2) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = k + 1;
    list[1] = 0;  /* dep var is the constant */
    for (i=0; i<k; i++) {
	list[i+2] = i + 1; 
    }

#ifdef BHHH_DEBUG
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
 *
 * Returns: pointer to newly allocated model_info, or
 * %NULL on failure.
 */
 
model_info *model_info_new (void)
{
    model_info *mi;

    mi = malloc(sizeof *mi);

    if (mi != NULL) {
	mi->theta = NULL;
	mi->series = NULL;
	mi->pmod = NULL;
	mi->VCV = NULL;
	mi->n_series = 0;
    }

    return mi;
}

static int bhhh_iter_info (int iter, double *theta, int m, double ll,
			   double steplength, PRN *prn)
{
    int i;

    pprintf(prn, "\n*** %s %d: theta, ll ***\n", _("iteration"), iter);

    for (i=0; i<m; i++) {
	if (i && i % 5 == 0) pputc(prn, '\n');
	if (na(theta[i]) || isnan(theta[i])) {
	    pprintf(prn, "Invalid value for theta[%d]\n", i);
	    return 1;
	}
	pprintf(prn, "%#12.5g ", theta[i]);
    }

    pprintf(prn, "\n    %s = %g, ll = %g\n", _("step length"),
	    steplength, ll);

    return 0;
}

/**
 * bhhh_max:
 * @loglik: pointer to function for calculating log-likelihood and
 * score matrix.
 * @X: original data set (passed on for use by @loglik).
 * @init_coeff: starting values for coefficients.
 * @minfo: model info struct, with some initialization carried out by
 * the caller.
 * @prn: printing struct for iteration info (or %NULL).
 *
 * Maximize likelihood using the BHHH conditional ML method,
 * implemented via iteration of the Outer Product of the Gradient 
 * (OPG) regression.
 * 
 * Returns: 0 on successful completion, non-zero error code otherwise.
 */

int bhhh_max (LL_FUNC loglik,
	      const double **X, 
	      const double *init_coeff,
	      model_info *minfo, 
	      PRN *prn)
{
    /* OPG model */
    MODEL *bmod;
    int *blist;

    /* temporary artificial dataset */
    double **tZ = NULL;
    DATAINFO *tinfo = NULL;

    int iters, itermax;
    double minstep = 1.0e-06; /* in arma, was 1.0e-08 */
    double crit = 1.0;
    double stepsize = 0.25;   /* in arma, was 0.125 */

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
    tinfo = create_new_dataset(&tZ, k + 1, minfo->n, 0);
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
#ifdef BHHH_DEBUG
	sprintf(tinfo->varname[i], "tZ%d", i);
#endif
	for (t=0; t<minfo->n; t++) {
	    tZ[i][t] = 0.0;
	}
    }

    /* initialize OPG model */
    bmod = gretl_model_new();

    iters = 0;
    itermax = get_maxiter();

    while (crit > minfo->tol && iters++ < itermax) {

	/* compute loglikelihood and score matrix */
	err = loglik(minfo->theta, X, tZ, minfo, 1); 
	if (err) {
	    pputs(prn, "Error calculating log-likelihood\n");
	    err = E_NOCONV;
	    break;
	}

#ifdef BHHH_DEBUG
	pprintf(prn, "Top of loop: ll = %g\n", minfo->ll);
#endif

	/* BHHH via OPG regression */
	*bmod = lsq(blist, &tZ, tinfo, OLS, OPT_A, 0.0);
	if (bmod->errcode) {
	    err = E_NOCONV;
	    break;
	}

#ifdef BHHH_DEBUG
	pprintf(prn, "Dataset for OPG regression:\n");
	printdata(blist, (const double **) tZ, tinfo, OPT_O, prn);
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

#ifdef BHHH_DEBUG
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
#ifdef BHHH_DEBUG
	    pprintf(prn, "bhhh loop: modified ll2 = %g\n", minfo->ll2);
#endif
	}

	if (err) break;

	/* actually update parameter estimates */
	for (i=0; i<k; i++) {
	    minfo->theta[i] = ctemp[i];
	}	

	/* double the steplength? */
	if (stepsize < 4.0) {
	    stepsize *= 2.0;
	}

	/* print interation info, if wanted */
	bhhh_iter_info(iters, minfo->theta, k, minfo->ll, stepsize, prn);

	crit = minfo->ll2 - minfo->ll;  
    }

    if (crit > minfo->tol || err != 0) {
	fprintf(stderr, "bhhh_max: crit = %g, tol = %g, err = %d\n",
		crit, minfo->tol, err);
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
	    int qr_bak = get_use_qr();

	    /* run OPG once more using QR, to get packed VCV */
	    set_use_qr(1);
	    *bmod = lsq(blist, &tZ, tinfo, OLS, OPT_A, 0.0);
	    set_use_qr(qr_bak);
	    minfo->pmod = bmod;
	    gretl_model_set_int(bmod, "iters", iters);
	} 
	minfo->iters = iters;
    }

    /* free all remaining temp stuff */
    free_Z(tZ, tinfo);
    free_datainfo(tinfo); 
    free(blist);

    if (bmod != minfo->pmod) {
	free(bmod);
    }

    return err;
}
