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
#include "internal.h"

#include "bhhh_max.h"

/* #define DEBUG */

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
 * Frees the dynamically allocated members of @model.
 * 
 */

void model_info_free (model_info *model)
{
    int i;

    free(model->theta);
    free(model->list);

    if (model->series != NULL) {
	for (i=0; i<model->n_series; i++) {
	    free(model->series[i]);
	}
	free(model->series);
    }

    if (model->VCV != NULL) {
	gretl_matrix_free(model->VCV);
    }
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

#ifdef DEBUG
    printlist(list, "OPG regression list");
#endif

    return list;
}

static int model_info_init (model_info *model, const double *init_coeff)
{
    int i, t, k, err = 0;
    int n_series = model->n_series;

    model->n = model->t2 + 1;
    k = model->k;  

    model->theta = NULL;
    model->series = NULL;
    model->list = NULL;
    model->pmod = NULL;
    model->VCV = NULL;

    model->list = make_opg_list(k);
    model->theta = malloc(k * sizeof *model->theta);

    if (model->list == NULL || model->theta == NULL) {
	model_info_free(model);
	return 1;
    }    

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
    for (i=0; i<k; i++) {
	model->theta[i] = init_coeff[i];
    }

    model->s2 = 1.0;

    return err;
}

static void bhhh_iter_info (int iter, double *theta, int m, double ll,
			    double steplength, PRN *prn)
{
    int i;

    pprintf(prn, "\n*** iteration %d: theta and ll ***\n", iter);
    for (i=0; i<m; i++) {
	if (i && i % 5 == 0) pputc(prn, '\n');
	pprintf(prn, "%#12.5g ", theta[i]);
    }
    pprintf(prn, "\n    steplength = %g, ll = %g\n", steplength, ll);
}

/**
 * bhhh_max:
 * @loglik: pointer to function for calculating log-likelihood and
 * score matrix.
 * @X: data set (used by @loglik).
 * @init_coeff: starting values for coefficients.
 * @model: model info struct, with some initialization carried out by
 * the caller.
 * @prn: printing struct for iteration info (or NULL).
 *
 * Maximize likelihood using the BHHH conditional ML method,
 * implemented via iteration of the Outer Product of the Gradient 
 * (OPG) regression.
 * 
 * Returns: 0 on successful completion, non-zero error code otherwise.
 */

int bhhh_max (int (*loglik) (double *, const double **, double **,
			     model_info *, int), 
	      const double **X, const double *init_coeff,
	      model_info *model, PRN *prn)
{
    /* OPG model */
    MODEL *bmod;

    /* temporary artificial dataset */
    double **tZ = NULL;
    DATAINFO *tinfo = NULL;

    int iters, itermax;
    double minstep = 1.0e-08; /* in Tobit, was 1.0e-06 */
    double crit = 1.0;
    double stepsize = 0.125;  /* in Tobit, was 0.25 */

    double *delta = NULL, *ctemp = NULL;

    int i, t, err, k;

    err = model_info_init(model, init_coeff);
    if (err) return E_ALLOC;

    k = model->k;

    delta = malloc(k * sizeof *delta);
    ctemp = malloc(k * sizeof *ctemp);
    if (delta == NULL || ctemp == NULL) {
	free(delta);
	free(ctemp);
	return E_ALLOC;
    }

    /* create temp dataset for OPG regression: k vars plus constant */    
    tinfo = create_new_dataset(&tZ, k + 1, model->n, 0);
    if (tinfo == NULL) {
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
    bmod = gretl_model_new(NULL);

    iters = 0;
    itermax = get_maxiter();

    while (crit > model->tol && iters++ < itermax && !err) {

	/* compute loglikelihood and score matrix */
	loglik(model->theta, X, tZ, model, 1); 

	/* BHHH via OPG regression */
	*bmod = lsq(model->list, &tZ, tinfo, OLS, OPT_A, 0.0);

	for (i=0; i<k; i++) {
	    delta[i] = bmod->coeff[i] * stepsize;
	    ctemp[i] = model->theta[i] + delta[i];
	} 
	
	clear_model(bmod, NULL);

	/* see if we've gone up... (0 = "don't compute score") */
	err = loglik(ctemp, X, tZ, model, 0); 

	while (model->ll2 < model->ll || err) { 
	    /* ... if not, halve steplength, as with ARMA models */
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

	/* double the steplength? */
	if (stepsize < 4.0) stepsize *= 2.0;

	/* actually update parameter estimates */
	for (i=0; i<k; i++) {
	    model->theta[i] += delta[i];
	}

	bhhh_iter_info(iters, model->theta, k, model->ll, stepsize, prn);

	if (model->s2 < 0.0) {
	    /* if the variance is negative here, we're stuck */
	    err = E_NOCONV;
	    break;
	}

	crit = model->ll2 - model->ll;  
    }

    if (crit > model->tol || err != 0) {
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
	    *bmod = lsq(model->list, &tZ, tinfo, OLS, OPT_A, 0.0);
	    set_use_qr(qr_bak);
	    model->pmod = bmod;
	    gretl_model_set_int(bmod, "iters", iters);
	} 
    }

    /* free all temp stuff */
    free_Z(tZ, tinfo);
    free_datainfo(tinfo);  

    if (bmod != model->pmod) {
	free(bmod);
    }

    return err;
}
