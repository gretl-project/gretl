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
   Trends in Machine Learning, Vol. 3, No. 1 (2010) 1­122.
*/

#include "libgretl.h"
#include "matrix_extra.h"
#include "version.h"

#define MAX_ITER 20000 // was 50 in Boyd
#define RELTOL 1e-3    // 1e-2 in Boyd
#define ABSTOL 1e-5    // 1e-4 in Boyd

static void vector_copy_values (gretl_vector *targ,
				const gretl_vector *src,
				int n)
{
    memcpy(targ->val, src->val, n * sizeof *targ->val);
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

static double xv_score (const gretl_matrix *A,
			const gretl_vector *b,
			const gretl_vector *z,
			gretl_vector *Azb)
{
    double SSR;

    gretl_matrix_multiply(A, z, Azb);
    gretl_matrix_subtract_from(Azb, b);
    SSR = gretl_vector_dot_product(Azb, Azb, NULL);

    return SSR / A->rows;
}

static void soft_threshold (gretl_vector *v, double k)
{
    double vi;
    int i;

    for (i=0; i<v->rows; i++) {
	vi = v->val[i];
	if (vi > k)       { v->val[i] = vi - k; }
	else if (vi < -k) { v->val[i] = vi + k; }
	else              { v->val[i] = 0; }
    }
}

static int get_cholesky_factor (const gretl_matrix *A,
				gretl_matrix *L,
				double rho)
{
    double d;
    int i, err = 0;

    if (A->rows >= A->cols) {
	/* "skinny": L = chol(AtA + rho*I) */
	gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
				  A, GRETL_MOD_NONE,
				  L, GRETL_MOD_NONE);
	for (i=0; i<A->cols; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + rho);
	}
	gretl_matrix_cholesky_decomp(L);
    } else {
	/* "fat": L = chol(I + 1/rho*AAt) */
	gretl_matrix_multiply_mod(A, GRETL_MOD_NONE,
				  A, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_NONE);
	gretl_matrix_multiply_by_scalar(L, 1/rho);
	for (i=0; i<A->rows; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + 1.0);
	}
	err = gretl_matrix_cholesky_decomp(L);
    }

    return err;
}

static int real_admm_lasso (const gretl_matrix *A,
			    const gretl_matrix *b,
			    gretl_bundle *bun)
{
    gretl_matrix_block *MB;
    double critmin = 1e200;
    double abstol;
    gretl_matrix *B = NULL;
    gretl_matrix *MSE = NULL;
    gretl_matrix *lfrac;
    double lmax, nrm2;
    double nxstack, nystack;
    double prires, dualres;
    double eps_pri, eps_dual;
    int verbose = 0;
    int ldim, nlam;
    int m, n, i, j;
    int jbest = 0;
    int err = 0;

    gretl_vector *x, *u, *z, *y, *r, *zprev, *zdiff;
    gretl_vector *q, *p, *w, *Atb, *Azb;
    gretl_matrix *L;

    double rho = gretl_bundle_get_scalar(bun, "rho", &err);
    int stdize = gretl_bundle_get_scalar(bun, "stdize", &err);
    double rho2 = rho * rho;

    if (gretl_bundle_has_key(bun, "lxv")) {
	lfrac = gretl_bundle_get_matrix(bun, "lxv", &err);
    } else {
	lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);
    }

    if (err) {
	return err;
    }

    /* dimensions */
    nlam = gretl_vector_get_length(lfrac);
    m = A->rows;
    n = A->cols;
    ldim = m >= n ? n : m;
    abstol = sqrt(n) * ABSTOL;

    MB = gretl_matrix_block_new(&x, n, 1, &u, n, 1,
				&z, n, 1, &y, n, 1,
				&r, n, 1, &zprev, n, 1,
				&zdiff, n, 1, &q, n, 1,
				&q, n, 1, &w, n, 1,
				&p, m, 1, &Atb, n, 1,
				&Azb, m, 1, &L, ldim, ldim,
				NULL);
    if (MB == NULL) {
	return E_ALLOC;
    }
    gretl_matrix_block_zero(MB);

    /* Precompute and cache factorizations */

    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);

    lmax = gretl_matrix_infinity_norm(Atb);

    fprintf(stderr, "lambda-max = %g\n", lmax);
    if (nlam > 1) {
	printf("using lambda-fraction sequence of length %d, starting at %g\n",
	       nlam, lfrac->val[0]);
    } else {
	printf("using lambda-fraction %g\n", lfrac->val[0]);
    }

    get_cholesky_factor(A, L, rho);

    B = gretl_zero_matrix_new(n + stdize, nlam);
    if (nlam > 0) {
	printf("     lambda  coeffs    criterion\n");
    }

    /* loop across lambda values */

    for (j=0; j<nlam && !err; j++) {
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
	    vector_copy_values(q, z, n);
	    gretl_matrix_subtract_from(q, u);
	    gretl_matrix_multiply_by_scalar(q, rho);
	    gretl_matrix_add_to(q, Atb);   // q = A^T b + rho*(z - u)

	    if (m >= n) {
		/* x = U \ (L \ q) */
		gretl_cholesky_solve(L, q);
		vector_copy_values(x, q, n);
	    } else {
		/* x = q/rho - 1/rho^2 * A^T * (U \ (L \ (A*q))) */
		gretl_matrix_multiply(A, q, p);
		err = gretl_cholesky_solve(L, p);
		gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
					  p, GRETL_MOD_NONE,
					  x, GRETL_MOD_NONE);
		gretl_matrix_multiply_by_scalar(x, -1/rho2);
		gretl_matrix_multiply_by_scalar(q, 1/rho);
		gretl_matrix_add_to(x, q);
	    }

	    vector_copy_values(w, x, n);
	    gretl_matrix_add_to(w, u);   /* w = x + u */

	    /* sqrt(sum ||r_i||_2^2) */
	    prires  = sqrt(gretl_vector_dot_product(r, r, NULL));
	    /* sqrt(sum ||r_i||_2^2) */
	    nxstack = sqrt(gretl_vector_dot_product(x, x, NULL));
	    /* sqrt(sum ||y_i||_2^2) */
	    nystack = gretl_vector_dot_product(u, u, NULL) / rho2;
	    nystack = sqrt(nystack);

	    vector_copy_values(zprev, z, n);
	    vector_copy_values(z, w, n);
	    soft_threshold(z, lambda/rho);

	    /* Termination checks */

	    /* dual residual */
	    vector_copy_values(zdiff, z, n);
	    gretl_matrix_subtract_from(zdiff, zprev);
	    /* ||s^k||_2^2 = N rho^2 ||z - zprev||_2^2 */
	    nrm2 = sqrt(gretl_vector_dot_product(zdiff, zdiff, NULL));
	    dualres = rho * nrm2;

	    /* compute primal and dual feasibility tolerances */
	    nrm2 = sqrt(gretl_vector_dot_product(z, z, NULL));
	    eps_pri  = abstol + RELTOL * fmax(nxstack, nrm2);
	    eps_dual = abstol + RELTOL * nystack;

	    if (iter > 1 && prires <= eps_pri && dualres <= eps_dual) {
		if (verbose) {
		    printf("breaking at iter %d: prires = %g, dualres = %g\n",
			   iter, prires, dualres);
		}
		conv = iter + 1;
		break;
	    }

	    /* Compute residual: r = x - z */
	    vector_copy_values(r, x, n);
	    gretl_matrix_subtract_from(r, z);

	    iter++;
	} /* end ADMM solve loop */

	if (!err) {
	    for (i=0; i<n; i++) {
		if (z->val[i] != 0.0) {
		    nnz++;
		}
		gretl_matrix_set(B, i+stdize, j, z->val[i]);
	    }
	    crit = objective(A, b, z, lambda, Azb);
	    printf("%#12.6g  %5d    %#.8g (%d iters%s)\n", lambda/m, nnz, crit,
		   conv ? conv : MAX_ITER, conv ? ", full convergence" : "");
	    if (crit < critmin) {
		critmin = crit;
		jbest = j;
	    }
	}
    } /* end loop over lambda values */

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
    gretl_matrix_block_destroy(MB);

    return err;
}

static int lasso_xv_round (const gretl_matrix *A,
			   const gretl_matrix *b,
			   const gretl_matrix *A_out,
			   const gretl_matrix *b_out,
			   gretl_bundle *bun)
{
    static gretl_matrix_block *MB;
    gretl_matrix *MSE;
    double abstol;
    double lmax, nrm2;
    double nxstack, nystack;
    double prires, dualres;
    double eps_pri, eps_dual;
    int verbose = 0;
    int ldim, nlam;
    int m, n, i, j;
    int err = 0;

    if (A == NULL) {
	/* cleanup signal */
	gretl_matrix_block_destroy(MB);
	MB = NULL;
	return 0;
    }

    static gretl_vector *x, *u, *z, *y, *r, *zprev, *zdiff;
    static gretl_vector *q, *p, *w, *Atb, *Azb;
    static gretl_matrix *L;

    double rho = gretl_bundle_get_scalar(bun, "rho", &err);
    int stdize = gretl_bundle_get_scalar(bun, "stdize", &err);
    gretl_matrix *lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);
    double rho2 = rho * rho;

    if (err) {
	return err;
    }

    /* dimensions */
    nlam = gretl_vector_get_length(lfrac);
    m = A->rows;
    n = A->cols;
    ldim = m >= n ? n : m;
    abstol = sqrt(n) * ABSTOL;

    if (MB == NULL) {
	MB = gretl_matrix_block_new(&x, n, 1, &u, n, 1,
				    &z, n, 1, &y, n, 1,
				    &r, n, 1, &zprev, n, 1,
				    &zdiff, n, 1, &q, n, 1,
				    &q, n, 1, &w, n, 1,
				    &p, m, 1, &Atb, n, 1,
				    &Azb, m, 1, &L, ldim, ldim,
				    NULL);
	if (MB == NULL) {
	    return E_ALLOC;
	}
	gretl_matrix_block_zero(MB);
    }
    // gretl_matrix_block_zero(MB);

    /* Precompute and cache factorizations */

    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);

    lmax = gretl_matrix_infinity_norm(Atb);

    get_cholesky_factor(A, L, rho);

    MSE = gretl_bundle_get_matrix(bun, "MSE", &err);

    /* loop across lambda values */

    for (j=0; j<nlam && !err; j++) {
	double lambda, score;
	int iter = 0;

	/* pull current lambda fraction and scale it */
	lambda = lfrac->val[j] * lmax;

	/* ADMM solver loop for given lambda */

	while (iter < MAX_ITER && !err) {
	    /* u-update: u = u + x - z */
	    gretl_matrix_subtract_from(x, z);
	    gretl_matrix_add_to(u, x);

	    /* x-update: x = (A^T A + rho I) \ (A^T b + rho z - y) */
	    vector_copy_values(q, z, n);
	    gretl_matrix_subtract_from(q, u);
	    gretl_matrix_multiply_by_scalar(q, rho);
	    gretl_matrix_add_to(q, Atb);   // q = A^T b + rho*(z - u)

	    if (m >= n) {
		/* x = U \ (L \ q) */
		gretl_cholesky_solve(L, q);
		vector_copy_values(x, q, n);
	    } else {
		/* x = q/rho - 1/rho^2 * A^T * (U \ (L \ (A*q))) */
		gretl_matrix_multiply(A, q, p);
		err = gretl_cholesky_solve(L, p);
		gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
					  p, GRETL_MOD_NONE,
					  x, GRETL_MOD_NONE);
		gretl_matrix_multiply_by_scalar(x, -1/rho2);
		gretl_matrix_multiply_by_scalar(q, 1/rho);
		gretl_matrix_add_to(x, q);
	    }

	    vector_copy_values(w, x, n);
	    gretl_matrix_add_to(w, u);   /* w = x + u */

	    /* sqrt(sum ||r_i||_2^2) */
	    prires  = sqrt(gretl_vector_dot_product(r, r, NULL));
	    /* sqrt(sum ||r_i||_2^2) */
	    nxstack = sqrt(gretl_vector_dot_product(x, x, NULL));
	    /* sqrt(sum ||y_i||_2^2) */
	    nystack = sqrt(gretl_vector_dot_product(u, u, NULL) / rho2);

	    vector_copy_values(zprev, z, n);
	    vector_copy_values(z, w, n);
	    soft_threshold(z, lambda/rho);

	    /* Termination checks */

	    /* dual residual */
	    vector_copy_values(zdiff, z, n);
	    gretl_matrix_subtract_from(zdiff, zprev);
	    /* ||s^k||_2^2 = N rho^2 ||z - zprev||_2^2 */
	    nrm2 = sqrt(gretl_vector_dot_product(zdiff, zdiff, NULL));
	    dualres = rho * nrm2;

	    /* compute primal and dual feasibility tolerances */
	    nrm2 = sqrt(gretl_vector_dot_product(z, z, NULL));
	    eps_pri  = abstol + RELTOL * fmax(nxstack, nrm2);
	    eps_dual = abstol + RELTOL * nystack;

	    if (iter > 1 && prires <= eps_pri && dualres <= eps_dual) {
		break;
	    }

	    /* Compute residual: r = x - z */
	    vector_copy_values(r, x, n);
	    gretl_matrix_subtract_from(r, z);

	    iter++;
	} /* end ADMM solve loop */

	/* cumulate out-of-sample MSE */
	gretl_matrix_reuse(Azb, A_out->rows, 1);
	score = xv_score(A_out, b_out, z, Azb);
	MSE->val[j] += score;
#if 0
	printf("  lambda-frac=%g, score=%g, cum=%g\n", lfrac->val[j],
	       score, MSE->val[j]);
#endif
	gretl_matrix_reuse(Azb, m, 1);
    } /* end loop over lambda values */

    /* cleanup */
    // gretl_matrix_block_destroy(MB);

    return err;
}

static int admm_lasso_xv (const gretl_matrix *A,
			  const gretl_matrix *b,
			  gretl_bundle *bun)
{
    gretl_matrix *A_est = NULL;
    gretl_matrix *b_est = NULL;
    gretl_matrix *A_out = NULL;
    gretl_matrix *b_out = NULL;
    gretl_matrix *lfrac = NULL;
    gretl_matrix *MSE = NULL;
    int i, j, f, nf, fsize;
    int nlam;
    int err = 0;

    nf = gretl_bundle_get_scalar(bun, "nfolds", &err);
    fsize = A->rows / nf;

    printf("admm_lasso_xv: nf=%d, fsize=%d\n", nf, fsize);

    A_est = gretl_matrix_alloc((nf-1)*fsize, A->cols);
    A_out = gretl_matrix_alloc(fsize, A->cols);
    b_est = gretl_matrix_alloc((nf-1)*fsize, 1);
    b_out = gretl_matrix_alloc(fsize, 1);

    lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);
    nlam = gretl_vector_get_length(lfrac);

    MSE = gretl_zero_matrix_new(nlam, 1);
    gretl_bundle_donate_data(bun, "MSE", MSE, GRETL_TYPE_MATRIX, 0);

    for (f=0; f<nf && !err; f++) {
	double aij;
	int ke, ko;

	/* subset the data */
	for (j=0; j<A->cols; j++) {
	    ke = ko = 0;
	    for (i=0; i<A->rows; i++) {
		aij = gretl_matrix_get(A, i, j);
		if (i/fsize == f) {
		    /* "out of sample" range */
		    gretl_matrix_set(A_out, ko, j, aij);
		    if (j == 0) {
			b_out->val[ko] = b->val[i];
		    }
		    ko++;
		} else {
		    /* estimation sample */
		    gretl_matrix_set(A_est, ke, j, aij);
		    if (j == 0) {
			b_est->val[ke] = b->val[i];
		    }
		    ke++;
		}
	    }
	}
	err = lasso_xv_round(A_est, b_est, A_out, b_out, bun);
    }

    /* send cleanup signal */
    lasso_xv_round(NULL, NULL, NULL, NULL, NULL);

    if (!err) {
	gretl_matrix *lxv = gretl_matrix_alloc(1, 1);
	double minMSE = MSE->val[0];
	int jbest = 0;

	for (j=0; j<nlam; j++) {
	    printf("s = %#g -> MSE %#g\n", lfrac->val[j], MSE->val[j]);
	    if (MSE->val[j] < minMSE) {
		jbest = j;
		minMSE = MSE->val[j];
	    }
	}
	printf("\nOut-of-sample MSE minimized at %g for s=%g\n",
	       minMSE, lfrac->val[jbest]);
	/* now determine coefficient vector on full training set */
	lxv->val[0] = lfrac->val[jbest];
	gretl_bundle_donate_data(bun, "lxv", lxv, GRETL_TYPE_MATRIX, 0);
	err = real_admm_lasso(A, b, bun);
    }

    gretl_matrix_free(A_est);
    gretl_matrix_free(b_est);
    gretl_matrix_free(A_out);
    gretl_matrix_free(b_out);

    return err;
}

int admm_lasso (const gretl_matrix *A,
		const gretl_matrix *b,
		gretl_bundle *bun)
{
    int xv, err = 0;

    xv = gretl_bundle_get_int(bun, "xvalidate", &err);

    if (xv) {
	return admm_lasso_xv(A, b, bun);
    } else {
	return real_admm_lasso(A, b, bun);
    }
}
