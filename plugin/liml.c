/*
 *  Copyright (c) 2004 by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"
#include "system.h"

#define LDEBUG 1

static int on_exo_list (const int *exlist, int v)
{
    int i;

    for (i=1; i<=exlist[0]; i++) {
	if (exlist[i] == v) return 1;
    }

    return 0;
}

/* compose E0 or E1 as in Greene, 4e, p. 686 */

static int resids_to_E (gretl_matrix *E, MODEL *lmod, int *reglist,
			const int *exlist, const int *list, int T,
			double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int i, j, t;
    int t1 = pdinfo->t1;
    int err = 0;

    /* loop across the k endog vars in list */
    j = 0;
    for (i=1; i<=list[0]; i++) {
	if (!on_exo_list(exlist, list[i])) {
	    reglist[1] = list[i];

#if LDEBUG
	    fprintf(stderr, "resids_to_E, aux reg: dependent var %s\n",
		    pdinfo->varname[reglist[1]]);
#endif
	    *lmod = lsq(reglist, pZ, pdinfo, OLS, OPT_NONE, 0.0);
	    if ((err = lmod->errcode)) {
		clear_model(lmod);
		break;
	    }

#if LDEBUG > 1
	    printmodel(lmod, pdinfo, OPT_NONE, prn);
#endif
	    /* put resids into appropriate column of E */
	    for (t=0; t<T; t++) {
		gretl_matrix_set(E, t, j, lmod->uhat[t + t1]);
	    }
	    j++;
	    clear_model(lmod);
	}
    }

    return err;
}

/* find the least characteristic root */

static double lambda_min (const double *lambda, int k)
{
    double lmin = 1.0;
    int i;

    for (i=0; i<k; i++) {
	if (i == 0) {
	    lmin = lambda[i];
	} else if (lambda[i] < lmin) {
	    lmin = lambda[i];
	}
    }

    return lmin;
}

static int *
liml_make_reglist (const int *exlist, const int *list,
		   int *k)
{
    int nexo = exlist[0];
    int *reglist;
    int i, j;

    reglist = malloc((nexo + 2) * sizeof *reglist);
    if (reglist == NULL) {
	return NULL;
    }

#if LDEBUG > 1
    fprintf(stderr, "Found %d exog vars; reglist allocated at length %d\n",
	    nexo, nexo + 2);
#endif

    /* at first, put all included exog vars in reglist */
    *k = 1;
    reglist[0] = 1;
    reglist[1] = 0;
    j = 2;
    for (i=2; i<=list[0]; i++) {
	if (on_exo_list(exlist, list[i])) {
	    reglist[0] += 1;
	    reglist[j++] = list[i];
	} else {
	    *k += 1;
	}
    }

    return reglist;
}

/* implements Greene's equation (16-28), 4e, p. 686 */

static gretl_matrix * 
liml_get_gamma_hat (const gretl_matrix *W0, const gretl_matrix *W1, 
		    double lmin)
{
    gretl_matrix *Wdiff = NULL;
    gretl_matrix *g = NULL;
    int k = W0->rows;
    double w0, w1;
    int i, j;
    int err = 0;

    Wdiff = gretl_matrix_alloc(k - 1, k - 1);
    if (Wdiff == NULL) {
	return NULL;
    }

    g = gretl_matrix_alloc(k - 1, 1);
    if (g == NULL) {
	gretl_matrix_free(Wdiff);
	return NULL;
    }

    for (i=0; i<Wdiff->rows; i++) {
	for (j=0; j<Wdiff->cols; j++) {
	    w0 = gretl_matrix_get(W0, i+1, j+1);
	    w1 = gretl_matrix_get(W1, i+1, j+1);
	    gretl_matrix_set(Wdiff, i, j, w0 - lmin * w1);
	}
    }

#if LDEBUG
    gretl_matrix_print(Wdiff, "W0 - lmin * W1", NULL);
#endif

    for (i=0; i<g->rows; i++) {
	w0 = gretl_matrix_get(W0, i+1, 0);
	w1 = gretl_matrix_get(W1, i+1, 0);
	gretl_vector_set(g, i, w0 - lmin * w1);
    }

#if LDEBUG
    gretl_matrix_print(g, "w0 - lmin * w1", NULL);
#endif

    err = gretl_LU_solve(Wdiff, g);
    if (err) {
	gretl_matrix_free(g);
	g = NULL;
    }

    gretl_matrix_free(Wdiff);

    return g;
}

/* implements the regression on p. 687 of Greene, 4e,
   following equation (16-28)
*/

static int 
liml_get_betahat (MODEL *olsmod, int *list, const int *exlist,
		  const gretl_matrix *g, double ***pZ, 
		  DATAINFO *pdinfo)
{
    MODEL bmod;
    double *y;
    double gj;
    int depvar = olsmod->list[1];
    int oldv = pdinfo->v;
    int i, j, t;
    int err = 0;

    /* create special dependent variable: the original y minus fitted
      y based on the RHS endogenous variables
    */

    y = malloc(pdinfo->n * sizeof *y);
    if (y == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	if (t < olsmod->t1 || t > olsmod->t2) {
	    y[t] = NADBL;
	} else {
	    y[t] = (*pZ)[depvar][t];
	    j = 0;
	    for (i=2; i<=olsmod->list[0]; i++) {
		if (!on_exo_list(exlist, olsmod->list[i])) {
		    gj = gretl_vector_get(g, j++);
		    y[t] -= (*pZ)[olsmod->list[i]][t] * gj;
		}
	    }
	}
    }

    if (dataset_add_allocated_var(y, pZ, pdinfo)) {
	free(y);
	return E_ALLOC;
    }
    
    list[0] = 1;
    list[1] = oldv; /* dep var is the newly added var */
    j = 2;
    /* put all included exog vars in list */
    for (i=2; i<=olsmod->list[0]; i++) {
	if (on_exo_list(exlist, olsmod->list[i])) {
	    list[0] += 1;
	    list[j++] = olsmod->list[i];
	} 
    }

    bmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    err = bmod.errcode;

    if (!err) {
	double bi;
	int k;

	k = j = 0;
	for (i=0; i<olsmod->ncoeff; i++) {
	    if (on_exo_list(exlist, olsmod->list[i+2])) {
		bi = bmod.coeff[j++];
	    } else {
		bi = gretl_vector_get(g, k++);
	    }
	    olsmod->coeff[i] = bi;
#if LDEBUG
	    fprintf(stderr, "beta_hat[%d] = %g\n", i, bi);
#endif
	}
    }

    clear_model(&bmod);
    dataset_drop_vars(pdinfo->v - oldv, pZ, pdinfo);

    return err;
}

static int liml_do_equation (gretl_equation_system *sys, int eq, double ***pZ,
			     DATAINFO *pdinfo, PRN *prn)
{
    const int *exlist = system_get_instr_vars(sys);
    const int *list = system_get_list(sys, eq);

    gretl_matrix *E = NULL;
    gretl_matrix *W0 = NULL;
    gretl_matrix *W1 = NULL;
    gretl_matrix *W2 = NULL;
    gretl_matrix *Inv = NULL;
    gretl_matrix *g = NULL;

    double *lambda, lmin = 1.0;
    MODEL *olsmod;
    MODEL lmod;
    int *reglist;
    int T = system_n_obs(sys);
    int i, k;
    int err = 0;

#if LDEBUG
    fprintf(stderr, "\nWorking on equation for %s\n", pdinfo->varname[list[1]]);
#endif

    olsmod = system_get_model(sys, eq);

#if LDEBUG
    printlist(olsmod->list, "OLS model list");
#endif

    reglist = liml_make_reglist(exlist, olsmod->list, &k);
    if (reglist == NULL) {
	return E_ALLOC;
    }

#if LDEBUG
    printf("number of endogenous vars in equation: k = %d\n", k);
#endif

    /* allocate matrices */
    E = gretl_matrix_alloc(T, k);
    W0 = gretl_matrix_alloc(k, k);
    W1 = gretl_matrix_alloc(k, k);    
    W2 = gretl_matrix_alloc(k, k);
    Inv = gretl_matrix_alloc(k, k);

    if (E == NULL || W0 == NULL || W1 == NULL || W2 == NULL ||
	Inv == NULL) {
	goto bailout;
    }

    err = resids_to_E(E, &lmod, reglist, exlist, list, 
		      T, pZ, pdinfo, prn);
    if (err) goto bailout;

    err = gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
				    E, GRETL_MOD_NONE,
				    W0);
    if (err) goto bailout;

#if LDEBUG
    gretl_matrix_print(W0, "W0", NULL);
#endif

    /* now re-estimate using all exogenous vars */
    reglist[0] = 1 + exlist[0];
    for (i=2; i<=reglist[0]; i++) {
	reglist[i] = exlist[i-1];
    }

    err = resids_to_E(E, &lmod, reglist, exlist, list, 
		      T, pZ, pdinfo, prn);
    if (err) goto bailout;
    
    err = gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
				    E, GRETL_MOD_NONE,
				    W1);
    if (err) goto bailout;

#if LDEBUG
    gretl_matrix_print(W1, "W1", NULL);
#endif

    gretl_matrix_copy_values(Inv, W1);

    err = gretl_invert_symmetric_matrix(Inv);
    if (err) goto bailout;

    err = gretl_matrix_multiply(Inv, W0, W2);
    if (err) goto bailout;
    
    lambda = gretl_general_matrix_eigenvals(W2, NULL);
    if (lambda == NULL) {
	err = 1;
	goto bailout;
    }

    lmin = lambda_min(lambda, k);
    printf("smallest eigenvalue: %g\n\n", lmin);
    free(lambda);

    g = liml_get_gamma_hat(W0, W1, lmin);
    if (g == NULL) {
	err = 1;
	goto bailout;
    }

#if LDEBUG
    for (i=0; i<k-1; i++) {
	printf("gamma_hat[%d] = %g\n", i, gretl_matrix_get(g, i, 0));
	    
    }
#endif

    err = liml_get_betahat(olsmod, reglist, exlist, g, pZ, pdinfo);
    if (err) {
	goto bailout;
    }

 bailout:

    free(reglist);
    gretl_matrix_free(E);
    gretl_matrix_free(W0);
    gretl_matrix_free(W1);
    gretl_matrix_free(W2);
    gretl_matrix_free(Inv);
    gretl_matrix_free(g);

    return err;
}

/* Driver function for LIML */

int liml_driver (gretl_equation_system *sys, double ***pZ, 
		 gretl_matrix *sigma, DATAINFO *pdinfo, 
		 PRN *prn)
{
    int g = system_n_equations(sys);
    int i, err = 0;

    pputs(prn, "\n*** LIML: experimental, work in progress ***\n\n");

    for (i=0; i<g; i++) {
	liml_do_equation(sys, i, pZ, pdinfo, prn);
    }

    return err;
}



    
	

    
    
    

    
    
