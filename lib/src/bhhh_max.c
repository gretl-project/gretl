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

#ifdef STANDALONE
# include <gretl/libgretl.h>
#else
# include "libgretl.h"
# include "internal.h"
#endif

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

void model_info_free (model_info *model)
{
    int i;

    free(model->theta);
    free(model->delta);
    free(model->deltmp);
    free(model->list);

    if (model->series != NULL) {
	for (i=0; i<model->n_series; i++) {
	    free(model->series[i]);
	}
	free(model->series);
    }

    gretl_matrix_free(model->VCV);
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

static int model_info_init (model_info *model, const double *init_coeff,
			    int n_init_coeff)
{
    int i, t, k, err = 0;
    int n_series = model->n_series;

    model->n = model->t2 + 1;
    k = model->k = n_init_coeff;  

    model->theta = model->delta = model->deltmp = NULL;
    model->series = NULL;
    model->list = NULL;

    model->list = make_opg_list(k + 1);
    model->theta = malloc((k + 1) * sizeof *model->theta);
    model->delta = malloc((k + 1) * sizeof *model->delta);
    model->deltmp = malloc((k + 1) * sizeof *model->deltmp);

    if (model->list == NULL || model->theta == NULL || 
	model->delta == NULL || model->deltmp == NULL) {
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

    /* initialize the coefficients */
    for (i=0; i<k; i++) {
	model->theta[i] = init_coeff[i];
    }

    /* initialize variance */
    model->theta[k] = 1.0;

    model->VCV = NULL;

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

int bhhh_max (int (*loglik) (double *, const double **, double **,
			     model_info *, int), 
	      const double **X, const double *init_coeff,
	      int n_init_coeff, model_info *model,
	      PRN *prn)
{
    /* OPG model */
    MODEL bmod;

    /* temporary artificial dataset */
    double **tZ = NULL;
    DATAINFO *tinfo = NULL;

    int iters, itermax;
    double tol = 1.0e-09;
    double smallstep = 1.0e-06;
    double crit = 1.0e20;
    double stepsize = 0.25;

    int i, t, err, k;

    err = model_info_init(model, init_coeff, n_init_coeff);
    if (err) return E_ALLOC;

    k = model->k;

    fprintf(stderr, "bhhh: k=%d, dataset will have %d vars\n", k, k + 2);

    /* create temp dataset for OPG regression: k+1 vars plus constant */    
    tinfo = create_new_dataset(&tZ, k + 2, model->n, 0);
    if (tinfo == NULL) {
	return E_ALLOC;
    } 

    /* zero the dataset */
    for (i=1; i<k+2; i++) {
	for (t=0; t<model->n; t++) {
	    tZ[i][t] = 0.0;
	}
    }

    iters = 0;
    itermax = get_maxiter();

    while (crit > tol && iters++ < itermax && !err) {

	/* compute loglikelihood and score matrix */
	loglik(model->theta, X, tZ, model, 1); 

	/* BHHH via OPG regression */
	bmod = lsq(model->list, &tZ, tinfo, OLS, OPT_A, 0.0);

	for (i=0; i<=k; i++) {
	    model->delta[i] = bmod.coeff[i] * stepsize;
	    model->deltmp[i] = model->theta[i] + model->delta[i];
	} 
	
	clear_model(&bmod, NULL);

	/* see if we've gone up... (0 = "don't compute score") */
	err = loglik(model->deltmp, X, tZ, model, 0); 

	while ((model->ll2 < model->ll || err) && stepsize > smallstep) { 
	    /* ... if not, halve steplength, as with ARMA models */
	    stepsize *= 0.5;
	    for (i=0; i<=k; i++) {
		model->delta[i] *= 0.5;
		model->deltmp[i] = model->theta[i] + model->delta[i];
	    }
	    err = loglik(model->deltmp, X, tZ, model, 0);
	}

	/* double the steplength? */
	if (stepsize < 4.0) stepsize *= 2.0;

	/* actually update parameter estimates */
	for (i=0; i<=k; i++) {
	    model->theta[i] += model->delta[i];
	}

	bhhh_iter_info(iters, model->theta, k+1, model->ll, stepsize, prn);

	if (model->theta[k] < 0.0) {
	    /* if the variance is negative here, we're stuck */
	    err = E_NOCONV;
	    break;
	}

	crit = model->ll2 - model->ll;  
    }

    if (crit > tol || err != 0) {
	err = E_NOCONV;
    }

    if (!err) {
	gretl_matrix *G, *VCV;

	/* fill out VCV */
	G = gretl_matrix_from_2d_array((const double **) tZ + 1, 
				       model->n, model->k + 1);
	VCV = gretl_matrix_vcv(G);
	model->VCV = VCV;
	gretl_matrix_free(G);

	/* free all temp stuff */
	free_Z(tZ, tinfo);
	free_datainfo(tinfo);
    }

    return err;
}
