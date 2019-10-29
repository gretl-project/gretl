/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2017 Allin Cottrell and Riccardo "Jack" Lucchetti
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

/* Code to use ADMM to solve the Lasso problem. Based on Boyd et al,
   "Distributed Optimization and Statistical Learning via the
   Alternating Direction Method of Multipliers", Foundations and
   Trends in Machine Learning Vol. 3, No. 1 (2010) 1­122.
*/

#include "libgretl.h"
#include "matrix_extra.h"
#include "version.h"

static gretl_matrix *gretl_vector_calloc (int n)
{
    return gretl_zero_matrix_new(n, 1);
}

static double abs_sum (const gretl_vector *z)
{
    int i, n = gretl_vector_get_length(z);
    double ret = 0;

    for (i=0; i<n; i++) {
	ret += fabs(z->val[i]);
    }

    return ret;
}

static double objective (const gretl_matrix *A,
			 const gretl_vector *b,
			 const gretl_vector *z,
			 double lambda,
			 gretl_vector *Azb)
{
    double Azb_nrm2;
    double obj = 0;

    gretl_matrix_multiply(A, z, Azb);
    gretl_matrix_subtract_from(Azb, b);
    Azb_nrm2 = gretl_vector_dot_product(Azb, Azb, NULL);
    obj = 0.5 * Azb_nrm2 + lambda * abs_sum(z);

    return obj / A->rows;
}

static void soft_threshold (gretl_vector *v, double k)
{
    double vi;
    int i;

    for (i = 0; i < v->rows; i++) {
	vi = v->val[i];
	if (vi > k)       { v->val[i] = vi - k; }
	else if (vi < -k) { v->val[i] = vi + k; }
	else              { v->val[i] = 0; }
    }
}

int admm_lasso (const gretl_matrix *A,
		const gretl_matrix *b,
		gretl_bundle *bun)
{
    const int MAX_ITER  = 20000; // was 50 in Boyd
    const double RELTOL = 1e-3;  // 1e-2 in Boyd
    const double ABSTOL = 1e-5;  // 1e-4 in Boyd
    double critmin = 1e200;
    gretl_matrix *L = NULL;
    gretl_matrix *B = NULL;
    double lmax, d, nrm2;
    int verbose = 0;
    int skinny, nlam;
    int m, n, i, j;
    int jbest = 0;
    int err = 0;

    gretl_matrix *lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);
    double rho = gretl_bundle_get_scalar(bun, "rho", &err);
    int stdize = gretl_bundle_get_scalar(bun, "stdize", &err);

    if (err) {
	return err;
    }

    m = A->rows;
    n = A->cols;
    skinny = (m >= n);

    gretl_vector *x      = gretl_vector_calloc(n);
    gretl_vector *u      = gretl_vector_calloc(n);
    gretl_vector *z      = gretl_vector_calloc(n);
    gretl_vector *y      = gretl_vector_calloc(n);
    gretl_vector *r      = gretl_vector_calloc(n);
    gretl_vector *zprev  = gretl_vector_calloc(n);
    gretl_vector *zdiff  = gretl_vector_calloc(n);

    gretl_vector *q      = gretl_vector_calloc(n);
    gretl_vector *w      = gretl_vector_calloc(n);
    gretl_vector *p      = gretl_vector_calloc(m);

    gretl_vector *Atb    = gretl_vector_calloc(n);
    gretl_vector *Azb    = gretl_vector_calloc(m);

    double nxstack  = 0;
    double nystack  = 0;
    double prires   = 0;
    double dualres  = 0;
    double eps_pri  = 0;
    double eps_dual = 0;

    /* Precompute and cache factorizations */

    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);

    lmax = gretl_matrix_infinity_norm(Atb);
    fprintf(stderr, "lambda-max = %g\n", lmax);

    nlam = gretl_vector_get_length(lfrac);
    printf("using lambda-fraction sequence of length %d, starting at %g\n",
	   nlam, lfrac->val[0]);

    /* Use the matrix inversion lemma for efficiency */
    if (skinny) {
	/* L = chol(AtA + rho*I) */
	L = gretl_matrix_alloc(n, n);
	gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
				  A, GRETL_MOD_NONE,
				  L, GRETL_MOD_NONE);
	for (i=0; i<n; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + rho);
	}
	gretl_matrix_cholesky_decomp(L);
    } else {
	/* L = chol(I + 1/rho*AAt) */
	L = gretl_matrix_alloc(m, m);
	gretl_matrix_multiply_mod(A, GRETL_MOD_NONE,
				  A, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_NONE);
	gretl_matrix_multiply_by_scalar(L, 1/rho);
	for (i=0; i<m; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + 1.0);
	}
	err = gretl_matrix_cholesky_decomp(L);
    }

    B = gretl_zero_matrix_new(n + stdize, nlam);

    if (nlam > 0) {
	printf("     lambda  coeffs    criterion\n");
    }

    /* loop across lambda values */

    for (j=0; j<nlam; j++) {
	double lambda, crit;
	int iter = 0;
	int conv = 0;
	int nnz = 0;

	/* pull current lambda fraction and scale it */
	lambda = lfrac->val[j] * lmax;

	/* ADMM solver loop for given lambda */

	while (iter < MAX_ITER && !err) {
	    /* u-update: u = u + x - z */
	    gretl_matrix_subtract_from(x, z);
	    gretl_matrix_add_to(u, x);

	    /* x-update: x = (A^T A + rho I) \ (A^T b + rho z - y) */
	    gretl_matrix_copy_values(q, z);
	    gretl_matrix_subtract_from(q, u);
	    gretl_matrix_multiply_by_scalar(q, rho);
	    gretl_matrix_add_to(q, Atb);   // q = A^T b + rho*(z - u)

	    if (skinny) {
		/* x = U \ (L \ q) */
		gretl_cholesky_solve(L, q);
		gretl_matrix_copy_values(x, q);
	    } else {
		/* x = q/rho - 1/rho^2 * A^T * (U \ (L \ (A*q))) */
		gretl_matrix_multiply(A, q, p);
		err = gretl_cholesky_solve(L, p);
		gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
					  p, GRETL_MOD_NONE,
					  x, GRETL_MOD_NONE);
		gretl_matrix_multiply_by_scalar(x, -1/(rho*rho));
		gretl_matrix_multiply_by_scalar(q, 1/rho);
		gretl_matrix_add_to(x, q);
	    }

	    gretl_matrix_copy_values(w, x);
	    gretl_matrix_add_to(w, u);   /* w = x + u */

	    /* sqrt(sum ||r_i||_2^2) */
	    prires  = sqrt(gretl_vector_dot_product(r, r, NULL));
	    /* sqrt(sum ||r_i||_2^2) */
	    nxstack = sqrt(gretl_vector_dot_product(x, x, NULL));
	    /* sqrt(sum ||y_i||_2^2) */
	    nystack = gretl_vector_dot_product(u, u, NULL) / pow(rho, 2);
	    nystack = sqrt(nystack);

	    gretl_matrix_copy_values(zprev, z);
	    gretl_matrix_copy_values(z, w);
	    soft_threshold(z, lambda/rho);

	    /* Termination checks */

	    /* dual residual */
	    gretl_matrix_copy_values(zdiff, z);
	    gretl_matrix_subtract_from(zdiff, zprev);
	    /* ||s^k||_2^2 = N rho^2 ||z - zprev||_2^2 */
	    nrm2 = sqrt(gretl_vector_dot_product(zdiff, zdiff, NULL));
	    dualres = rho * nrm2;

	    /* compute primal and dual feasibility tolerances */
	    nrm2 = sqrt(gretl_vector_dot_product(z, z, NULL));
	    eps_pri  = sqrt(n)*ABSTOL + RELTOL * fmax(nxstack, nrm2);
	    eps_dual = sqrt(n)*ABSTOL + RELTOL * nystack;

	    if (iter > 1 && prires <= eps_pri && dualres <= eps_dual) {
		if (verbose) {
		    printf("breaking at iter %d: prires = %g, dualres = %g\n",
			   iter, prires, dualres);
		}
		conv = iter + 1;
		break;
	    }

	    /* Compute residual: r = x - z */
	    gretl_matrix_copy_values(r, x);
	    gretl_matrix_subtract_from(r, z);

	    iter++;
	} /* end ADMM solve */

	for (i=0; i<n; i++) {
	    if (z->val[i] != 0.0) {
		nnz++;
	    }
	    gretl_matrix_set(B, i+stdize, j, z->val[i]);
	}
	crit = objective(A, b, z, lambda, Azb);
	printf("%#12.6g  %5d    %#g (%d iters%s)\n", lambda/m, nnz, crit,
	       conv ? conv : MAX_ITER, conv ? ", full convergence" : "");
	if (crit < critmin) {
	    critmin = crit;
	    jbest = j;
	}
    } /* end lambda values */

    gretl_bundle_set_scalar(bun, "lmax", lmax);
    gretl_bundle_set_scalar(bun, "jbest", jbest + 1);
    gretl_bundle_set_scalar(bun, "lfbest", lfrac->val[jbest]);
    gretl_bundle_set_scalar(bun, "crit", critmin);
    if (nlam > 1) {
	gretl_bundle_donate_data(bun, "B", B, GRETL_TYPE_MATRIX, 0);
    } else {
	gretl_bundle_donate_data(bun, "b", B, GRETL_TYPE_MATRIX, 0);
	gretl_bundle_set_scalar(bun, "lambda", lfrac->val[0] * lmax);
    }

    /* cleanup */
    gretl_matrix_free(L);
    gretl_matrix_free(x);
    gretl_matrix_free(u);
    gretl_matrix_free(z);
    gretl_matrix_free(y);
    gretl_matrix_free(r);
    gretl_matrix_free(w);
    gretl_matrix_free(zprev);
    gretl_matrix_free(zdiff);
    gretl_matrix_free(q);
    gretl_matrix_free(Atb);
    gretl_matrix_free(Azb);
    gretl_matrix_free(p);

    return err;
}
