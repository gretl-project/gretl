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
#include "gretl_private.h"

#include "bhhh_max.h"

#undef BHHH_DEBUG

struct _model_info {

    /* members that may be set by caller of bhhh_max, via accessor
       functions: */

    int k;              /* number of parameters */
    int p, q, r;        /* for use with ARMA: AR, MA orders and number of
                           other regressors.  Otherwise unused. */
    int t1, t2;         /* starting and ending point of sample */
    int n_series;       /* number of additional series needed in the
                           likelihood and/or score calculations */
    double tol;         /* tolerance for convergence */
    unsigned char opts; /* options from among bhhh_opts */

    /* members set within bhhh_max: */

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
 * @model: model info pointer.
 *
 * Frees the dynamically allocated members of @model, then
 * frees @model itself.
 * 
 */

void model_info_free (model_info *model)
{
    int i;

    free(model->theta);

    if (model->series != NULL) {
	for (i=0; i<model->n_series; i++) {
	    free(model->series[i]);
	}
	free(model->series);
    }

    if (model->VCV != NULL) {
	gretl_matrix_free(model->VCV);
    }

    free(model);
}

/**
 * model_info_capture_OPG_model:
 * @model: model info pointer.
 *
 * Returns: pointer to a gretl #MODEL, which contains the
 * results of the OPG (Outer Product of the Gradient) regression
 * associated with @model.
 */

MODEL *model_info_capture_OPG_model (model_info *model)
{
    MODEL *pmod = model->pmod;

    model->pmod = NULL;
    return pmod;
}

/**
 * model_info_get_VCV:
 * @model: model info pointer.
 *
 * Returns: pointer to the covariance matrix of @model.
 */

gretl_matrix *model_info_get_VCV (model_info *model)
{
    return model->VCV;
}

/**
 * model_info_get_theta:
 * @model: model info pointer.
 *
 * Returns: parameter vector @theta from @model.
 */

double *model_info_get_theta (model_info *model)
{
    return model->theta;
}

/**
 * model_info_get_t1:
 * @model: model info pointer.
 *
 * Returns: the (zero-based) start of the sample range for 
 * @model.
 */

int model_info_get_t1 (const model_info *model)
{
    return model->t1;
}

/**
 * model_info_get_t2:
 * @model: model info pointer.
 *
 * Returns: the end of the sample range for @model.
 */

int model_info_get_t2 (const model_info *model)
{
    return model->t2;
}

/**
 * model_info_get_n:
 * @model: model info pointer.
 *
 * Returns: the number of observations used in @model.
 */

int model_info_get_n (const model_info *model)
{
    return model->n;
} 

/**
 * model_info_get_iters:
 * @model: model info pointer.
 *
 * Returns: the number of iterations taken in estimating
 * @model.
 */  

int model_info_get_iters (const model_info *model)
{
    return model->iters;
}

/**
 * model_info_get_pqr:
 * @model: model info pointer.
 * @p: pointer to receive the AR order of @model.
 * @q: pointer to receive the MA order of @model.
 * @r: pointer to receive the number of regular regressors in 
 * ARMA(X) @model.
 *
 */  
 
void model_info_get_pqr (const model_info *model, 
			 int *p, int *q, int *r)
{
    *p = model->p;
    *q = model->q;
    *r = model->r;
}

/**
 * model_info_get_series:
 * @model: model info pointer.
 *
 * Returns: FIXME.
 */  

double **model_info_get_series (const model_info *model)
{
    return model->series;
}

/**
 * model_info_get_ll:
 * @model: model info pointer.
 *
 * Returns: the log-likelihood for @model.
 */ 

double model_info_get_ll (const model_info *model)
{
    return model->ll;
}

/**
 * model_info_set_ll:
 * @model: model info pointer.
 * @ll: log-likehood.
 * @do_score:
 *
 * Sets the log-likelihood for @model.  If @do_score is non-zero,
 * sets the primary ll value, otherwise sets the secondary value.
 * FIXME: explain this.
 */ 

void model_info_set_ll (model_info *model, double ll, int do_score)
{
    if (do_score) {
	model->ll = ll;
    } else {
	model->ll2 = ll;
    }
}

/**
 * model_info_set_opts:
 * @model: model info pointer.
 * @opts: option flags to set.
 *
 * Sets the option flags for @model. 
 * FIXME: explain this.
 */ 

void model_info_set_opts (model_info *model, unsigned char opts)
{
    model->opts = opts;
}

/**
 * model_info_set_tol:
 * @model: model info pointer.
 * @tol: tolerance for convergence of estimates.
 *
 * Sets the convergence tolerance for @model. 
 * FIXME: explain this.
 */ 

void model_info_set_tol (model_info *model, double tol)
{
    model->tol = tol;
}

/**
 * model_info_set_pqr:
 * @model: model info pointer.
 * @p: AR order for @model.
 * @q: MA order for @model.
 * @r: Number of regular regressors in ARMA(X) @model.
 *
 * Sets the specified parameters for @model.
 */ 

void model_info_set_pqr (model_info *model, int p, int q, int r)
{
    model->p = p;
    model->q = q;
    model->r = r;
    model->k = p + q + r + 1;
}

/**
 * model_info_set_k:
 * @model: model info pointer.
 * @k: number of regressors in (non-ARMA) model.
 *
 * Sets the number of regressors.
 */ 

void model_info_set_k (model_info *model, int k)
{
    model->k = k;
    model->p = model->q = model->r = 0;
}

/**
 * model_info_get_k:
 * @model: model info pointer.
 *
 * Returns: the number of regressors in (non-ARMA) @model.
 */ 

int model_info_get_k (model_info *model)
{
    return model->k;
}

/**
 * model_info_set_n_series:
 * @model: model info pointer.
 * @n: number of data series.
 *
 * Sets the number of auxiliary data series needed for @model.
 */ 

void model_info_set_n_series (model_info *model, int n)
{
    model->n_series = n;
}

/**
 * model_info_set_t1_t2:
 * @model: model info pointer.
 * @t1: starting observation number (zero-based).
 * @t2: ending observation number.
 *
 * Sets the sample range for @model.
 */ 

void model_info_set_t1_t2 (model_info *model, int t1, int t2)
{
    model->t1 = t1;
    model->t2 = t2;
    model->n = t2 + 1;
}

/* Below: construct the regression list for the OPG regression, with
   the appropriate indices into the temporary artificial dataset.
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

static int model_info_init (model_info *model, const double *init_coeff)
{
    int i, t, k, err = 0;
    int n_series = model->n_series;

    k = model->k;  

    model->theta = malloc(k * sizeof *model->theta);

    if (model->theta == NULL) {
	model_info_free(model);
	return 1;
    }    

    /* allocate series */
    if (n_series > 0) {
	model->series = malloc(n_series * sizeof *model->series);
	if (model->series == NULL) {
	    model_info_free(model);
	    return 1;
	}
	model->n_series = 0;
	for (i=0; i<n_series; i++) {
	    model->series[i] = malloc(model->n * sizeof **model->series);
	    if (model->series[i] == NULL) {
		model_info_free(model);
		return 1;
	    }
	    for (t=0; t<model->n; t++) {
		model->series[i][t] = 0.0;
	    }
	    model->n_series += 1;
	}
    }

    /* initialize parameters */
#ifdef BHHH_DEBUG
    fprintf(stderr, "model_info_init: initializing theta\n");
#endif
    for (i=0; i<k; i++) {
#ifdef BHHH_DEBUG
	fprintf(stderr, "theta[%d] = %g\n", i, init_coeff[i]);
#endif
	model->theta[i] = init_coeff[i];
    }

    return err;
}

/**
 * model_info_new:
 *
 * Returns: pointer to newly allocated model_info.
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

    pprintf(prn, "\n*** %s %d: theta, ll ***\n", ("iteration"), iter);
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
 * @X: data set (used by @loglik).
 * @init_coeff: starting values for coefficients.
 * @model: model info struct, with some initialization carried out by
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
	      model_info *model, 
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

    err = model_info_init(model, init_coeff);
    if (err) return E_ALLOC;

    k = model->k;
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
    tinfo = create_new_dataset(&tZ, k + 1, model->n, 0);
    if (tinfo == NULL) {
	free(delta);
	free(ctemp);
	free(blist);
	return E_ALLOC;
    } 

    /* respect the incoming sample range */
    tinfo->t1 = model->t1;
    tinfo->t2 = model->t2;

    /* zero the dataset */
    for (i=1; i<=k; i++) {
	for (t=0; t<model->n; t++) {
	    tZ[i][t] = 0.0;
	}
    }

    /* initialize OPG model */
    bmod = gretl_model_new();

    iters = 0;
    itermax = get_maxiter();

    while (crit > model->tol && iters++ < itermax) {

	/* compute loglikelihood and score matrix */
	err = loglik(model->theta, X, tZ, model, 1); 
	if (err) {
	    pputs(prn, "Error calculating log-likelihood\n");
	    err = E_NOCONV;
	    break;
	}

	/* BHHH via OPG regression */
	*bmod = lsq(blist, &tZ, tinfo, OLS, OPT_A, 0.0);
	if (bmod->errcode) {
	    err = E_NOCONV;
	    break;
	}

#ifdef BHHH_DEBUG
	pprintf(prn, "Dataset for OPG regression\n");
	printdata(blist, (const double **) tZ, tinfo, OPT_O, prn);
	pprintf(prn, "OLS model, OPG regression\n");
	printmodel(bmod, tinfo, 0, prn);
#endif

	for (i=0; i<k; i++) {
	    delta[i] = bmod->coeff[i] * stepsize;
	    ctemp[i] = model->theta[i] + delta[i];
	} 
	
	clear_model(bmod);

	/* see if we've gone up... (0 = "don't compute score") */
	err = loglik(ctemp, X, tZ, model, 0); 

	while (model->ll2 < model->ll || err) { 
	    /* ... if not, halve the steplength */
	    stepsize *= 0.5;
	    if (stepsize < minstep) {
		err = E_NOCONV;
		break;
	    }
	    for (i=0; i<k; i++) {
		delta[i] *= 0.5;
		ctemp[i] = model->theta[i] + delta[i];
	    }
	    err = loglik(ctemp, X, tZ, model, 0);
	}

	if (err) break;

	/* actually update parameter estimates */
	for (i=0; i<k; i++) {
	    model->theta[i] = ctemp[i];
	}	

	/* double the steplength? */
	if (stepsize < 4.0) stepsize *= 2.0;

	/* print interation info, if wanted */
	bhhh_iter_info(iters, model->theta, k, model->ll, stepsize, prn);

	crit = model->ll2 - model->ll;  
    }

    if (crit > model->tol || err != 0) {
	fprintf(stderr, "bhhh_max: crit = %g, tol = %g, err = %d\n",
		crit, model->tol, err);
	err = E_NOCONV;
    }

    free(delta);
    free(ctemp);

    if (!err) {
	if (model->opts & FULL_VCV_MATRIX) {
	    gretl_matrix *G, *VCV;

	    G = gretl_matrix_from_2d_array((const double **) tZ + 1, 
					   model->n, model->k);
	    VCV = gretl_matrix_vcv(G);
	    model->VCV = VCV;
	    gretl_matrix_free(G);
	}
	if (model->opts & PRESERVE_OPG_MODEL) {
	    int qr_bak = get_use_qr();

	    /* run OPG once more using QR, to get packed VCV */
	    set_use_qr(1);
	    *bmod = lsq(blist, &tZ, tinfo, OLS, OPT_A, 0.0);
	    set_use_qr(qr_bak);
	    model->pmod = bmod;
	    gretl_model_set_int(bmod, "iters", iters);
	} 
	model->iters = iters;
    }

    /* free all remaining temp stuff */
    free_Z(tZ, tinfo);
    free_datainfo(tinfo); 
    free(blist);

    if (bmod != model->pmod) {
	free(bmod);
    }

    return err;
}
