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

#define LDEBUG 0

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
			double ***pZ, DATAINFO *pdinfo)
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
liml_make_reglist (const gretl_equation_system *sys, const int *list, 
		   int *k)
{
    const int *exlist = system_get_instr_vars(sys);
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
	    /* an endogenous var */
	    *k += 1;
	}
    }

    return reglist;
}

/* set the special LIML k-class data on the model: these data will be
   retrieved when calculating the LIML coefficients and their
   covariance matrix
*/

static int 
liml_set_model_data (MODEL *pmod, const gretl_matrix *E, 
		     const int *exlist, const int *list, 
		     int T, int t1, int n, double lmin,
		     const double **Z)
{
    double *Xi = NULL;
    double *ymod = NULL;
    int ysize = n * sizeof *ymod;
    double yt, xit, eit;
    int m = list[0] - 1;
    int i, j, t;
    int err = 0;

    ymod = malloc(ysize);
    if (ymod == NULL) {
	return 1;
    }

    for (t=0; t<n; t++) {
	ymod[t] = NADBL;
    }

    for (t=0; t<T; t++) {
	yt = Z[list[1]][t + t1];
	eit = gretl_matrix_get(E, t, 0);
	ymod[t + t1] = yt - lmin * eit;
	j = 1;
	for (i=0; i<m; i++) {
	    if (on_exo_list(exlist, list[i+2])) {
		continue;
	    }
	    Xi = tsls_get_Xi(pmod, Z, i);
	    if (Xi == NULL) {
		err = 1;
		break;
	    }
	    xit = Z[list[i+2]][t + t1];
	    eit = gretl_matrix_get(E, t, j++);
	    Xi[t + t1] = xit - lmin * eit;
	}
	if (err) break;
    }

    if (err) {
	free(ymod);
    } else {
	gretl_model_set_data(pmod, "liml_y", ymod, ysize);
    }

    return err;
}

static int liml_do_equation (gretl_equation_system *sys, int eq, 
			     double ***pZ, DATAINFO *pdinfo, 
			     PRN *prn)
{
    const int *exlist = system_get_instr_vars(sys);
    const int *list = system_get_list(sys, eq);

    gretl_matrix *E = NULL;
    gretl_matrix *W0 = NULL;
    gretl_matrix *W1 = NULL;
    gretl_matrix *W2 = NULL;
    gretl_matrix *Inv = NULL;
    gretl_matrix *g = NULL;
    gretl_matrix *V = NULL;

    double *lambda, lmin = 1.0;
    double ll = 0.0;
    MODEL *pmod;
    MODEL lmod;
    int *reglist;
    int idf;
    int T = system_n_obs(sys);
    int i, k;
    int err = 0;

#if LDEBUG
    fprintf(stderr, "\nWorking on equation for %s\n", pdinfo->varname[list[1]]);
    printlist(list, "original equation list");
#endif

    /* degrees of freedom for over-identification test:
       total exog vars minus number of params in equation.
       FIXME: handle case of idf = 0.
    */
    idf = exlist[0] - (list[0] - 1);
    
    /* pointer to corresponding TSLS model */
    pmod = system_get_model(sys, eq);

    /* make regression list using only included exogenous vars */
    reglist = liml_make_reglist(sys, list, &k);
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
		      T, pZ, pdinfo);
    if (err) goto bailout;

    err = gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
				    E, GRETL_MOD_NONE,
				    W0);
    if (err) goto bailout;

#if LDEBUG
    gretl_matrix_print(W0, "W0", NULL);
#endif

    /* re-set the regression list using all exogenous vars */
    reglist[0] = 1 + exlist[0];
    for (i=2; i<=reglist[0]; i++) {
	reglist[i] = exlist[i-1];
    }

    err = resids_to_E(E, &lmod, reglist, exlist, list, 
		      T, pZ, pdinfo);
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
    free(lambda);
    gretl_model_set_double(pmod, "lmin", lmin);
    gretl_model_set_int(pmod, "idf", idf);

#if LDEBUG
    printf("smallest eigenvalue: %g\n\n", lmin);
#endif

    err = liml_set_model_data(pmod, E, exlist, list, T,
			      pdinfo->t1, pdinfo->n, lmin,
			      (const double **) *pZ);
    if (err) {
	fprintf(stderr, "error in liml_set_model_data\n");
    }

    /* compute and set log-likelihood, etc */
    ll = system_n_equations(sys) * LN_2_PI;
    ll += log(lmin);
    ll += gretl_matrix_log_determinant(W1);
    ll *= -(T / 2.0);
    pmod->lnL = ll;

#if LDEBUG
    printf("log-likelihood = %g\n", ll);
#endif

    mle_aic_bic(pmod, 0); /* check the "0" (additional params) here */

 bailout:

    free(reglist);
    gretl_matrix_free(E);
    gretl_matrix_free(W0);
    gretl_matrix_free(W1);
    gretl_matrix_free(W2);
    gretl_matrix_free(Inv);
    gretl_matrix_free(g);
    gretl_matrix_free(V);

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
#if LDEBUG > 1
	printmodel(system_get_model(sys, i), pdinfo, OPT_NONE, prn);
#endif
	liml_do_equation(sys, i, pZ, pdinfo, prn);
    }

    return err;
}



    
	

    
    
    

    
    
