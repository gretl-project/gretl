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

/* Implements Greene's equation (16-28), 4e, p. 686. */

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

/* Implements the regression on p. 687 of Greene, 4e,
   following equation (16-28).
*/

static int 
liml_get_betahat (MODEL *olsmod, int *reglist, const int *list,
		  const int *exlist, const gretl_matrix *g, 
		  double ***pZ, DATAINFO *pdinfo)
{
    MODEL bmod;
    double *y;
    double gj;
    int depvar = list[1];
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
	    for (i=2; i<=list[0]; i++) {
		if (!on_exo_list(exlist, list[i])) {
		    gj = gretl_vector_get(g, j++);
		    y[t] -= (*pZ)[list[i]][t] * gj;
		}
	    }
	}
    }

    if (dataset_add_allocated_var(y, pZ, pdinfo)) {
	free(y);
	return E_ALLOC;
    }
    
    reglist[0] = 1;
    reglist[1] = oldv; /* dep var is the newly added var */
    j = 2;
    /* put all included exog vars on the RHS */
    for (i=2; i<=list[0]; i++) {
	if (on_exo_list(exlist, list[i])) {
	    reglist[0] += 1;
	    reglist[j++] = list[i];
	} 
    }

    bmod = lsq(reglist, pZ, pdinfo, OLS, OPT_A, 0.0);
    err = bmod.errcode;

    if (!err) {
	double bi;
	int k;

	k = j = 0;
	for (i=0; i<olsmod->ncoeff; i++) {
	    if (on_exo_list(exlist, list[i+2])) {
		bi = bmod.coeff[j++];
	    } else {
		bi = gretl_vector_get(g, k++);
	    }
	    olsmod->coeff[i] = bi;
	}
    }

    clear_model(&bmod);
    dataset_drop_vars(pdinfo->v - oldv, pZ, pdinfo);

    return err;
}

static void 
liml_set_residuals_etc (MODEL *pmod, const int *list, 
			const gretl_matrix *XTX, 
			const double **Z)
{
    double yhat;
    int depvar = list[1];
    int i, t;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	yhat = 0.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    yhat += pmod->coeff[i] * Z[list[i+2]][t];
	}
	pmod->yhat[t] = yhat;
	pmod->uhat[t] = Z[depvar][t] - yhat;
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    /* could apply a df correction here (TSP seems to) */
    pmod->sigma = sqrt(pmod->ess / pmod->nobs); 

    for (i=0; i<pmod->ncoeff; i++) {
	double xii = gretl_matrix_get(XTX, i, i);

	pmod->sderr[i] = pmod->sigma * sqrt(xii);
    } 
}

/* Compute and invert the special counterpart of "X-transpose X"
   that forms the basis for var(betahat) in LIML.  See Davidson
   and MacKinnon, Econometric Theory and Methods, equation (12.96).
*/

static gretl_matrix *
liml_get_XTX_inverse (const gretl_matrix *E, const int *exlist, 
		      const int *list, int T, int t1, double lmin, 
		      const double **Z)
{
    gretl_matrix *X = NULL;
    gretl_matrix *Xmod = NULL;
    gretl_matrix *XTX = NULL;
    double xit, eit;
    int m = list[0] - 1;
    int i, j, t;
    int err = 0;

    X = gretl_matrix_alloc(T, m);
    Xmod = gretl_matrix_alloc(T, m);
    XTX = gretl_matrix_alloc(m, m);

    if (X == NULL || Xmod == NULL || XTX == NULL) {
	err = 1;
	goto bailout;
    }

    for (t=0; t<T; t++) {
	j = 1;
	for (i=2; i<=list[0]; i++) {
	    xit = Z[list[i]][t + t1];
	    gretl_matrix_set(X, t, i-2, xit);
	    if (on_exo_list(exlist, list[i])) {
		gretl_matrix_set(Xmod, t, i-2, xit);
	    } else {
		eit = gretl_matrix_get(E, t, j++);
		gretl_matrix_set(Xmod, t, i-2, xit - lmin * eit);
	    } 
	}
    }

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
			      Xmod, GRETL_MOD_NONE,
			      XTX);
    err = gretl_invert_general_matrix(XTX);

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(Xmod);

    if (err) {
	gretl_matrix_free(XTX);
	XTX = NULL;
    }

    return XTX;
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
    gretl_matrix *XTX = NULL;

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
    
    /* pointer to TSLS model */
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

    /* estimates of coeffs on included endogenous vars */
    g = liml_get_gamma_hat(W0, W1, lmin);
    if (g == NULL) {
	err = 1;
	goto bailout;
    }

    /* estimates of coeffs on included exogenous vars */
    err = liml_get_betahat(pmod, reglist, list, exlist, g, 
			   pZ, pdinfo);
    if (err) {
	goto bailout;
    }

    /* for covariance matrix of coefficient estimates */
    XTX = liml_get_XTX_inverse(E, exlist, list, T, pdinfo->t1, lmin,
			       (const double **) *pZ);
    if (XTX == NULL) {
	err = 1;
	goto bailout;
    }

    /* correct yhat, uhat, ESS, sigma, standard errors */
    liml_set_residuals_etc(pmod, list, XTX, (const double **) *pZ);

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
    gretl_matrix_free(XTX);

    return err;
}

/* re-write the cross-equation VCV using the LIML residuals */

static void 
liml_rewrite_sigma (gretl_matrix *sigma, gretl_equation_system *sys,
		    int t1)
{
    const MODEL *imod, *jmod;
    int g = system_n_equations(sys);
    int T = system_n_obs(sys);
    int i, j, t;
    double xx;

    for (i=0; i<g; i++) {
	imod = system_get_model(sys, i);
	for (j=i; j<g; j++) {
	    jmod = system_get_model(sys, j);
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += imod->uhat[t + t1] * jmod->uhat[t + t1];
	    }
	    xx /= T;
	    gretl_matrix_set(sigma, i, j, xx);
	    if (j != i) {
		gretl_matrix_set(sigma, j, i, xx);
	    }
	}
    }
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

    liml_rewrite_sigma(sigma, sys, pdinfo->t1);

    return err;
}



    
	

    
    
    

    
    
