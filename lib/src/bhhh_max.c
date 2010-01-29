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
#include "bhhh_max.h"

#define BHHH_DEBUG 0

/**
 * bhhh_max:
 * @theta: array of adjustable coefficients.
 * @k: number of elements in @theta.
 * @G: T x @k matrix to hold the gradient.
 * @loglik: pointer to function for calculating log-likelihood and
 * score.
 * @toler: tolerance for convergence.
 * @itcount: location to receive count of iterations.
 * @data: pointer to be passed to @loglik.
 * @V: matrix to receive covariance, or %NULL.
 * @opt: can include %OPT_V for verbose output.
 * @prn: printing struct for iteration info (or %NULL).
 *
 * Maximize likelihood using the BHHH conditional ML method,
 * implemented via iteration of the Outer Product of the Gradient 
 * (OPG) regression with line search.
 *
 * @loglik is called to calculate the log-likelihood for the model
 * in question.  The parameters passed to this function are:
 * (1) the current array of estimated coefficients; (2) @G; 
 * (3) the @data pointer; and (4) an integer that is 1 if 
 * the gradient should be calculated in @G, otherwise 0.
 *
 * Note that @G does nto have to initialized on entry.
 *
 * For an example of the use of such a function, see arma.c in the
 * %plugin directory of the gretl source.
 * 
 * Returns: 0 on successful completion, non-zero error code otherwise.
 */

int bhhh_max (double *theta, int k, 
	      gretl_matrix *G,
	      BHHH_FUNC loglik,
	      double toler, int *itcount,
	      void *data, 
	      gretl_matrix *V,
	      gretlopt opt,
	      PRN *prn)
{
    gretl_matrix *c = NULL;
    gretl_matrix *gcoeff = NULL;
    double *delta = NULL, *ctemp = NULL;
    int iters, itermax;
    double minstep = 1.0e-06;
    double crit = 1.0;
    double stepsize = 0.25;
    double ll2, ll = 0.0;
    int i, T, err = 0;

    if (gretl_matrix_cols(G) != k) {
	return E_NONCONF;
    }

    T = gretl_matrix_rows(G);
    c = gretl_unit_matrix_new(T, 1);
    gcoeff = gretl_column_vector_alloc(k);

    if (c == NULL || gcoeff == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    delta = malloc(k * sizeof *delta);
    ctemp = malloc(k * sizeof *ctemp);

    if (delta == NULL || ctemp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    iters = 0;
    itermax = libset_get_int(BHHH_MAXITER);

    while (crit > toler && iters++ < itermax) {

	/* compute loglikelihood and score */
	ll = loglik(theta, G, data, 1, &err); 
	if (err) {
	    pputs(prn, "Error calculating log-likelihood\n");
	    break;
	}

#if BHHH_DEBUG
	pprintf(prn, "Top of loop: ll = %g\n", ll);
	gretl_matrix_print(S, "RHS in OPG regression");
#endif

	/* BHHH via OPG regression */
	err = gretl_matrix_ols(c, G, gcoeff, NULL, NULL, NULL);

	if (err) {
	    fprintf(stderr, "BHHH OLS error code = %d\n", err);
	    break;
	} 

	for (i=0; i<k; i++) {
	    delta[i] = gcoeff->val[i] * stepsize;
	    ctemp[i] = theta[i] + delta[i];
	} 
	
	/* see if we've gone up...  (0 means don't compute score) */
	ll2 = loglik(ctemp, G, data, 0, &err); 

#if BHHH_DEBUG
	pprintf(prn, "bhhh loop: initial ll2 = %g\n", ll2);
#endif

	while (err || ll2 < ll) { 
	    /* ... if not, halve the steplength */
	    stepsize *= 0.5;
	    if (stepsize < minstep) {
		err = E_NOCONV;
		break;
	    }
	    for (i=0; i<k; i++) {
		delta[i] *= 0.5;
		ctemp[i] = theta[i] + delta[i];
	    }
	    ll2 = loglik(ctemp, G, data, 0, &err);
#if BHHH_DEBUG
	    pprintf(prn, "bhhh loop: modified ll2 = %g\n", ll2);
#endif
	}

	if (err) break;

	/* actually update parameter estimates */
	for (i=0; i<k; i++) {
	    theta[i] = ctemp[i];
	}	

	/* double the steplength? (was < 4.0 below) */
	if (stepsize < 1.0) {
	    stepsize *= 2.0;
	}

	/* print iteration info, if wanted */
	if (opt & OPT_V) {
	    print_iter_info(iters, ll, C_LOGLIK, k, theta, delta, 
			    stepsize, prn);
	}

	crit = ll2 - ll;  
    }

    if (opt & OPT_V) {
	print_iter_info(-1, ll, C_LOGLIK, k, theta, delta, 
			stepsize, prn);
    }

    if (crit > toler && !err) {
	err = E_NOCONV;
    }

    if (err) {
	fprintf(stderr, "bhhh_max: iters = %d, crit = %g, tol = %g, err = %d\n",
		iters, crit, toler, err);
    } else {
	if (V != NULL) {
	    /* run OPG once more, to get VCV */
	    double s2 = 0.0;

	    err = gretl_matrix_ols(c, G, gcoeff, V, NULL, &s2);
	} 
	*itcount = iters;
    }

 bailout:

    gretl_matrix_free(c);
    gretl_matrix_free(gcoeff);

    free(delta);
    free(ctemp);

    return err;
}
