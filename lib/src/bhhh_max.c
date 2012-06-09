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

/* just experimental at this point: call @func to compute the
   per-observation contributions to the log-likelihood in @L,
   and use that info to compute @G numerically.
*/

static double bhhh_numeric (double *b, int k, 
			    gretl_matrix *L,
			    gretl_matrix *G,
			    BHHH_FUNC func, void *data,
			    int *err)
{
    const double h = 1.0e-8;
    double gti, ret;
    int i, t, T = G->rows;

    /* loglik at current evaluation point */
    ret = func(b, L, data, 0, err); 

    for (i=0; i<k && !*err; i++) {
	double ll, bi0 = b[i];

	b[i] = bi0 - h;
	ll = func(b, L, data, 0, err);
	
	if (na(ll)) {
	    ret = NADBL;
	    break;
	}

	for (t=0; t<T; t++) {
	    gretl_matrix_set(G, t, i, L->val[t]);
	}

	b[i] = bi0 + h;
	ll = func(b, L, data, 0, err); 

	b[i] = bi0;

	if (na(ll)) {
	    ret = NADBL;
	    break;
	}

	for (t=0; t<T; t++) {
	    gti = gretl_matrix_get(G, t, i);
	    gti = (L->val[t] - gti) / (2.0 * h);
	    gretl_matrix_set(G, t, i, gti);
	}
    }

    return ret;
}

/* in case we want to compute the score numerically */

static gretl_matrix *make_score_matrix (gretl_matrix *L, int k, int *err)
{
    int T = gretl_matrix_rows(L);
    gretl_matrix *G = NULL;

    if (T == 0) {
	*err = E_DATA;
    } else if (T <= k) {
	*err = E_DF;
    } else {
	G = gretl_zero_matrix_new(T, k);
	if (G == NULL) {
	    *err = E_ALLOC;
	}
    }

    return G;
}

/**
 * bhhh_max:
 * @theta: array of adjustable coefficients.
 * @k: number of elements in @theta.
 * @M: T x @k matrix to hold the gradient.
 * @callback: pointer to function for calculating log-likelihood and
 * score.
 * @toler: tolerance for convergence.
 * @fncount: location to receive count of function evaluations.
 * @grcount: location to receive count of gradient evaluations.
 * @data: pointer to be passed to @loglik.
 * @V: matrix in which to store covariance, or %NULL.
 * @opt: can include %OPT_V for verbose output.
 * @prn: printing struct for iteration info (or %NULL).
 *
 * Maximize likelihood using the BHHH method, implemented via 
 * iteration of the Outer Product of the Gradient (OPG) regression 
 * with line search.
 *
 * @callback is called to calculate the loglikelihood for the model
 * in question.  The parameters passed to this function are:
 * (1) the current array of estimated coefficients; (2) @G; 
 * (3) the @data pointer; (4) an integer indicator that is 1 if 
 * the gradient should be calculated in @G, 0 if only the
 * loglikelihood is needed; and (5) an int pointer to receive
 * an error code. The return value from @callback should be the
 * loglikelihood (or #NADBL on error).
 *
 * For an example of the use of such a function, see arma.c in the
 * %plugin directory of the gretl source.
 * 
 * Returns: 0 on successful completion, non-zero error code otherwise.
 */

int bhhh_max (double *theta, int k, 
	      gretl_matrix *M,
	      BHHH_FUNC callback,
	      double toler, 
	      int *fncount, int *grcount,
	      void *data, 
	      gretl_matrix *V,
	      gretlopt opt,
	      PRN *prn)
{
    gretl_matrix *c = NULL;
    gretl_matrix *g = NULL;
    gretl_matrix *G = NULL;
    double *delta = NULL, *ctemp = NULL;
    double *grad;
    int itermax, iter = 0;
    int fcount = 0;
    int gcount = 0;
    double minstep = 1.0e-06;
    double crit = 1.0;
    double stepsize = 0.25;
    double ll2, ll = 0.0;
    int numeric = 0;
    int i, T, err = 0;

    if (opt & OPT_N) {
	/* numerical score (not ready) */
	G = make_score_matrix(M, k, &err);
	if (err) {
	    return err;
	}
	numeric = 1;
    } else {	
	/* analytical score */
	if (gretl_matrix_cols(M) != k) {
	    return E_NONCONF;
	} else {
	    G = M;
	}
    } 

    T = gretl_matrix_rows(G);
    c = gretl_unit_matrix_new(T, 1);
    g = gretl_column_vector_alloc(k);

    if (c == NULL || g == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    delta = malloc(k * sizeof *delta);
    ctemp = malloc(k * sizeof *ctemp);

    if (delta == NULL || ctemp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    itermax = libset_get_int(BHHH_MAXITER);
    grad = g->val;

    while (crit > toler && iter++ < itermax) {

#if BHHH_DEBUG
	fprintf(stderr, "BHHH: iter = %d\n", iter);
#endif

	/* compute loglikelihood and score */
	if (numeric) {
	    ll = bhhh_numeric(theta, k, M, G, callback, data, &err);
	} else {
	    ll = callback(theta, G, data, 1, &err); 
	}

	fcount++;
	gcount++;

	if (err) {
	    pputs(prn, "Error calculating log-likelihood\n");
	    break;
	}

#if BHHH_DEBUG
	fprintf(stderr, "Top of loop: ll = %g\n", ll);
#endif
#if BHHH_DEBUG > 1
	gretl_matrix_print(G, "RHS in OPG regression");
#endif

	/* BHHH via OPG regression */
	err = gretl_matrix_ols(c, G, g, NULL, NULL, NULL);

	if (err) {
	    fprintf(stderr, "BHHH OLS error code = %d\n", err);
	    break;
	} 

	for (i=0; i<k; i++) {
	    delta[i] = grad[i] * stepsize;
	    ctemp[i] = theta[i] + delta[i];
	} 
	
	/* see if we've gone up...  (0 means don't compute score) */
	ll2 = callback(ctemp, G, data, 0, &err); 
	fcount++;

#if BHHH_DEBUG
	fprintf(stderr, "bhhh loop: initial ll2 = %#.14g\n", ll2);
#endif

	while (err || ll2 < ll) { 
	    /* ... or if not, halve the steplength */
	    stepsize *= 0.5;
	    if (stepsize < minstep) {
		fprintf(stderr, "BHHH: hit minimum step size %g\n", minstep);
		err = E_NOCONV;
		break;
	    }
	    for (i=0; i<k; i++) {
		delta[i] *= 0.5;
		ctemp[i] = theta[i] + delta[i];
	    }
	    ll2 = callback(ctemp, G, data, 0, &err);
	    fcount++;
#if BHHH_DEBUG
	    fprintf(stderr, "bhhh loop: modified ll2 = %#.14g\n", ll2);
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
	    print_iter_info(iter, ll, C_LOGLIK, k, theta, grad, 
			    stepsize, prn);
	}

	crit = ll2 - ll;  
    }

    *fncount = fcount;
    *grcount = gcount;

    if (opt & OPT_V) {
	print_iter_info(-1, ll, C_LOGLIK, k, theta, grad, 
			stepsize, prn);
    }

    if (crit > toler && !err) {
	err = E_NOCONV;
    }

    if (err) {
	fprintf(stderr, "bhhh_max: iters = %d, crit = %g, tol = %g, err = %d\n",
		iter, crit, toler, err);
    } else if (V != NULL) {
	/* run OPG once more, to get VCV */
	double s2 = 0.0;

	err = gretl_matrix_ols(c, G, g, V, NULL, &s2);
    }

 bailout:

    gretl_matrix_free(c);
    gretl_matrix_free(g);

    if (G != M) {
	gretl_matrix_free(G);
    }

    free(delta);
    free(ctemp);

    return err;
}
