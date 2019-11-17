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

#if defined(USE_AVX)
# define USE_SIMD
# if defined(HAVE_IMMINTRIN_H)
#  include <immintrin.h>
# else
#  include <mmintrin.h>
#  include <xmmintrin.h>
#  include <emmintrin.h>
# endif
#endif

#define MAX_ITER 20000
#define USE_YBAR 1

double reltol = 1.0e-4;
double abstol = 1.0e-6;
#if USE_YBAR
double ybar = 0.0;
#endif

enum {
    CRIT_MSE,
    CRIT_MAE,
    CRIT_PCC
};

static const char *crit_string (int crit)
{
    if (crit == CRIT_MSE) {
	return "MSE";
    } else if (crit == CRIT_MAE) {
	return "MAE";
    } else {
	return "pc correct";
    }
}

static int randomize_rows (gretl_matrix *A, gretl_matrix *b)
{
    gretl_vector *vp;
    double x, tmp;
    int i, j, src;

    vp = gretl_matrix_alloc(A->rows, 1);
    if (vp == NULL) {
	return E_ALLOC;
    }

    fill_permutation_vector(vp, A->rows);

    for (i=0; i<A->rows; i++) {
	src = vp->val[i] - 1;
	if (src == i) {
	    continue;
	}
	for (j=0; j<A->cols; j++) {
	    tmp = gretl_matrix_get(A, i, j);
	    x = gretl_matrix_get(A, src, j);
	    gretl_matrix_set(A, i, j, x);
	    gretl_matrix_set(A, src, j, tmp);
	}
	tmp = b->val[i];
	b->val[i] = b->val[src];
	b->val[src] = tmp;
    }

    gretl_matrix_free(vp);

    return 0;
}

static void vector_copy_values (gretl_vector *targ,
				const gretl_vector *src,
				int n)
{
    memcpy(targ->val, src->val, n * sizeof *targ->val);
}

#if defined(USE_SIMD)

static void vector_add_into (const gretl_vector *a,
			     const gretl_vector *b,
			     gretl_vector *c, int n)
{
    const double *ax = a->val;
    const double *bx = b->val;
    double *cx = c->val;
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, c256;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(ax);
	b256 = _mm256_loadu_pd(bx);
	c256 = _mm256_add_pd(a256, b256);
	_mm256_storeu_pd(cx, c256);
	ax += 4;
	bx += 4;
	cx += 4;
    }
    for (i=0; i<rem; i++) {
	cx[i] = ax[i] + bx[i];
    }
}

static void vector_add_to (gretl_vector *a,
			   const gretl_vector *b,
			   int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, sum;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(ax);
	b256 = _mm256_loadu_pd(bx);
	sum = _mm256_add_pd(a256, b256);
	_mm256_storeu_pd(ax, sum);
	ax += 4;
	bx += 4;
    }
    for (i=0; i<rem; i++) {
	ax[i] += bx[i];
    }
}

/* a = a - b */

static void vector_subtract_from (gretl_vector *a,
				  const gretl_vector *b,
				  int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, dif;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(ax);
	b256 = _mm256_loadu_pd(bx);
	dif = _mm256_sub_pd(a256, b256);
	_mm256_storeu_pd(ax, dif);
	ax += 4;
	bx += 4;
    }
    for (i=0; i<rem; i++) {
	ax[i] -= bx[i];
    }
}

/* c = a - b */

static void vector_subtract_into (const gretl_vector *a,
				  const gretl_vector *b,
				  gretl_vector *c, int n,
				  int cumulate)
{
    const double *ax = a->val;
    const double *bx = b->val;
    double *cx = c->val;
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, c256;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(ax);
	b256 = _mm256_loadu_pd(bx);
	if (cumulate) {
	    __m256d d256 = _mm256_sub_pd(a256, b256);

	    c256 = _mm256_loadu_pd(cx);
	    d256 = _mm256_add_pd(c256, d256);
	    _mm256_storeu_pd(cx, d256);
	} else {
	    c256 = _mm256_sub_pd(a256, b256);
	    _mm256_storeu_pd(cx, c256);
	}
	ax += 4;
	bx += 4;
	cx += 4;
    }
    for (i=0; i<rem; i++) {
	if (cumulate) {
	    cx[i] += ax[i] - bx[i];
	} else {
	    cx[i] = ax[i] - bx[i];
	}
    }
}

/* compute q = rho * (z - u) + A'b */

static inline void compute_q (gretl_vector *q,
			      const gretl_vector *z,
			      const gretl_vector *u,
			      const gretl_vector *a, /* A'b */
			      double rho, int n)
{
    __m256d z256, u256, a256;
    __m256d r256, tmp;
    const double *zx = z->val;
    const double *ux = u->val;
    const double *ax = a->val;
    double *qx = q->val;
    const int mul = rho != 1.0;
    int imax = n / 4;
    int rem = n % 4;
    int i;

    if (mul) {
	/* broadcast rho */
	r256 = _mm256_broadcast_sd(&rho);
    }

    /* FIXME check for _mm256_fmadd_pd() and use it
       if available? */

    for (i=0; i<imax; i++) {
	z256 = _mm256_loadu_pd(zx);
	u256 = _mm256_loadu_pd(ux);
	a256 = _mm256_loadu_pd(ax);
	/* subtract u from z */
	tmp = _mm256_sub_pd(z256, u256);
	if (mul) {
	    /* multiply by rho */
	    tmp = _mm256_mul_pd(tmp, r256);
	}
	/* add a */
	tmp = _mm256_add_pd(tmp, a256);
	/* write result into q */
	_mm256_storeu_pd(qx, tmp);
	zx += 4;
	ux += 4;
	ax += 4;
	qx += 4;
    }

    for (i=0; i<rem; i++) {
	if (mul) {
	    qx[i] = rho * (zx[i] - ux[i]) + ax[i];
	} else {
	    qx[i] = zx[i] - ux[i] + ax[i];
	}
    }
}

static void vector_add_scalar (gretl_vector *v,
			       double x, int n)
{
    double *vx = v->val;
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d x256, v256, sum;

    /* broadcast x */
    x256 = _mm256_broadcast_sd(&x);

    for (i=0; i<imax; i++) {
	v256 = _mm256_loadu_pd(vx);
	sum = _mm256_add_pd(v256, x256);
	_mm256_storeu_pd(vx, sum);
	vx += 4;
    }
    for (i=0; i<rem; i++) {
	vx[i] += x;
    }
}

#else

static void vector_add_into (const gretl_vector *a,
			     const gretl_vector *b,
			     gretl_vector *c, int n)
{
    int i;

    for (i=0; i<n; i++) {
	c->val[i] = a->val[i] + b->val[i];
    }
}

static void vector_add_to (gretl_vector *a,
			   const gretl_vector *b,
			   int n)
{
    int i;

    for (i=0; i<n; i++) {
	a->val[i] += b->val[i];
    }
}

static void vector_subtract_from (gretl_vector *a,
				  const gretl_vector *b,
				  int n)
{
    int i;

    for (i=0; i<n; i++) {
	a->val[i] -= b->val[i];
    }
}

static void vector_subtract_into (const gretl_vector *a,
				  const gretl_vector *b,
				  gretl_vector *c, int n,
				  int cumulate)
{
    int i;

    for (i=0; i<n; i++) {
	if (cumulate) {
	    c->val[i] += a->val[i] - b->val[i];
	} else {
	    c->val[i] = a->val[i] - b->val[i];
	}
    }
}

static inline void compute_q (gretl_vector *q,
			      const gretl_vector *z,
			      const gretl_vector *u,
			      const gretl_vector *Atb,
			      double rho, int n)
{
    const int mul = rho != 1.0;
    int i;

    for (i=0; i<n; i++) {
	if (mul) {
	    q->val[i] = rho * (z->val[i] - u->val[i]) + Atb->val[i];
	} else {
	    q->val[i] = z->val[i] - u->val[i] + Atb->val[i];
	}
    }
}

static void vector_add_scalar (gretl_vector *v,
			       double x, int n)
{
    int i;

    for (i=0; i<n; i++) {
	v->val[i] += x;
    }
}

#endif /* AVX or not */

static double abs_sum (const gretl_vector *z)
{
    const int n = gretl_vector_get_length(z);
    double ret = 0;
    int i;

    for (i=0; i<n; i++) {
	ret += fabs(z->val[i]);
    }

    return ret;
}

/* calculate the lasso criterion */

static double objective (const gretl_matrix *A,
			 const gretl_vector *b,
			 const gretl_vector *z,
			 double lambda,
			 gretl_vector *Azb)
{
    double SSR;
    double obj = 0;

    gretl_matrix_multiply(A, z, Azb);
#if USE_YBAR
    if (ybar != 0) {
	vector_add_scalar(Azb, ybar, A->rows);
    }
#endif
    vector_subtract_from(Azb, b, A->rows);
    SSR = gretl_vector_dot_product(Azb, Azb, NULL);
    obj = 0.5 * SSR + lambda * abs_sum(z);

    return obj / A->rows;
}

/* calculate the cross validation criterion */

static double xv_score (const gretl_matrix *A,
			const gretl_vector *b,
			const gretl_vector *z,
			gretl_vector *Azb,
			int crit_type)
{
    double sum = 0;

    /* get fitted values */
    gretl_matrix_multiply(A, z, Azb);
#if USE_YBAR
    if (ybar != 0) {
	vector_add_scalar(Azb, ybar, A->rows);
    }
#endif
    if (crit_type == CRIT_PCC) {
	/* count incorrect classifications */
	double yhat;
	int i, icc = 0;

	for (i=0; i<A->rows; i++) {
	    yhat = gretl_round(Azb->val[i]);
	    icc += yhat != b->val[i];
	}
	sum = 100 * icc;
    } else {
	/* compute and process residuals */
	vector_subtract_from(Azb, b, A->rows);
	if (crit_type == CRIT_MSE) {
	    sum = gretl_vector_dot_product(Azb, Azb, NULL);
	} else {
	    sum = abs_sum(Azb);
	}
    }

    return sum / A->rows;
}

static void soft_threshold (gretl_vector *v, double lambda,
			    double rho)
{
    double vi, k;
    int i;

    k = rho == 1.0 ? lambda : lambda / rho;

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
    int i;

    if (A->rows >= A->cols) {
	/* "skinny": L = chol(A'A + rho*I) */
	gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
				  A, GRETL_MOD_NONE,
				  L, GRETL_MOD_NONE);
	for (i=0; i<A->cols; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + rho);
	}
    } else {
	/* "fat": L = chol(I + 1/rho*AA') */
	gretl_matrix_multiply_mod(A, GRETL_MOD_NONE,
				  A, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_NONE);
	if (rho != 1.0) {
	    gretl_matrix_multiply_by_scalar(L, 1/rho);
	}
	for (i=0; i<A->rows; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + 1.0);
	}
    }

    return gretl_matrix_cholesky_decomp(L);
}

static int admm_iteration (const gretl_matrix *A,
			   const gretl_matrix *L,
			   const gretl_vector *Atb,
			   gretl_vector *x, gretl_vector *z,
			   gretl_vector *u, gretl_vector *q,
			   gretl_vector *p, gretl_vector *r,
			   gretl_vector *zprev, gretl_vector *zdiff,
			   double lambda, double rho,
			   int *iters)
{
    double nxstack, nystack;
    double prires, dualres;
    double eps_pri, eps_dual;
    double nrm2, rho2 = rho*rho;
    int mul = rho != 1.0;
    int n = A->cols;
    int iter = 0;
    int err = 0;

    while (iter < MAX_ITER && !err) {
	/* u-update: u = u + r */
	vector_add_to(u, r, n);

	/* x-update: x = (A^T A + rho I) \ (A^T b + rho z - y) */

	compute_q(q, z, u, Atb, rho, n);
	if (A->rows >= A->cols) {
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
	    if (mul) {
		gretl_matrix_multiply_by_scalar(x, -1/rho2);
		gretl_matrix_multiply_by_scalar(q, 1/rho);
	    } else {
		gretl_matrix_multiply_by_scalar(x, -1);
	    }
	    vector_add_to(x, q, n);
	}

	/* sqrt(sum ||r_i||_2^2) */
	prires  = sqrt(gretl_vector_dot_product(r, r, NULL));
	/* sqrt(sum ||r_i||_2^2) */
	nxstack = sqrt(gretl_vector_dot_product(x, x, NULL));
	/* sqrt(sum ||y_i||_2^2) */
	if (mul) {
	    nystack = gretl_vector_dot_product(u, u, NULL) / rho2;
	    nystack = sqrt(nystack);
	} else {
	    nystack = sqrt(gretl_vector_dot_product(u, u, NULL));
	}

	vector_copy_values(zprev, z, n);
	vector_add_into(x, u, z, n);
	soft_threshold(z, lambda, rho);

	/* Termination checks */

	/* dual residual */
	vector_subtract_into(z, zprev, zdiff, n, 0); /* zdiff = z - zprev */
	/* ||s^k||_2^2 = N rho^2 ||z - zprev||_2^2 */
	nrm2 = sqrt(gretl_vector_dot_product(zdiff, zdiff, NULL));
	dualres = mul ? rho * nrm2 : nrm2;

	/* compute primal and dual feasibility tolerances */
	nrm2 = sqrt(gretl_vector_dot_product(z, z, NULL));
	eps_pri  = abstol + reltol * fmax(nxstack, nrm2);
	eps_dual = abstol + reltol * nystack;

	if (iter > 1 && prires <= eps_pri && dualres <= eps_dual) {
	    break;
	}

	/* Compute residual: r = x - z */
	vector_subtract_into(x, z, r, n, 0);

	iter++;
    }

    *iters = iter;

    return err;
}

static int real_admm_lasso (const gretl_matrix *A,
			    const gretl_matrix *b,
			    gretl_bundle *bun,
			    double rho)
{
    gretl_matrix_block *MB;
    double critmin = 1e200;
    gretl_matrix *B = NULL;
    gretl_matrix *lfrac;
    double lmax;
    int ldim, nlam;
    int m, n, i, j;
    int jbest = 0;
    int err = 0;

    gretl_vector *x, *u, *z, *y, *r, *zprev, *zdiff;
    gretl_vector *q, *Atb, *m1;
    gretl_matrix *L;

    int stdize = gretl_bundle_get_int(bun, "stdize", &err);
    int xval = gretl_bundle_get_int(bun, "xvalidate", &err);
    int verbo = gretl_bundle_get_int(bun, "verbosity", &err);

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

    MB = gretl_matrix_block_new(&x, n, 1, &u, n, 1,
				&z, n, 1, &y, n, 1,
				&r, n, 1, &zprev, n, 1,
				&zdiff, n, 1, &q, n, 1,
				&m1, m, 1, &Atb, n, 1,
				&L, ldim, ldim, NULL);
    if (MB == NULL) {
	return E_ALLOC;
    }
    gretl_matrix_block_zero(MB);

    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);

    lmax = gretl_matrix_infinity_norm(Atb);

    if (verbo > 0) {
	if (nlam > 1) {
	    printf("using lambda-fraction sequence of length %d, starting at %g\n",
		   nlam, lfrac->val[0]);
	} else {
	    printf("using lambda-fraction %g\n", lfrac->val[0]);
	}
    }

    get_cholesky_factor(A, L, rho);

    B = gretl_zero_matrix_new(n + stdize, nlam);
    if (verbo > 0 && nlam > 1) {
	printf("     lambda  coeffs    criterion\n");
    }

    for (j=0; j<nlam && !err; j++) {
	/* loop across lambda values */
	double crit, lambda = lfrac->val[j] * lmax;
	int iters = 0;
	int nnz = 0;

	err = admm_iteration(A, L, Atb, x, z, u, q, m1, r, zprev, zdiff,
			     lambda, rho, &iters);

	if (!err) {
	    for (i=0; i<n; i++) {
		if (z->val[i] != 0.0) {
		    nnz++;
		}
		gretl_matrix_set(B, i+stdize, j, z->val[i]);
	    }
	    crit = objective(A, b, z, lambda, m1);
	    if (verbo > 0 && nlam > 1) {
		printf("%#12.6g  %5d    %#.8g (%d iters)\n",
		       lambda/m, nnz, crit, iters);
	    }
	    if (crit < critmin) {
		critmin = crit;
		jbest = j;
	    }
	}
    }

    gretl_bundle_set_scalar(bun, "lmax", lmax);
    if (nlam > 1 || xval) {
	gretl_bundle_set_scalar(bun, "best_idx", jbest + 1);
	gretl_bundle_set_scalar(bun, "lfbest", lfrac->val[jbest]);
    }
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
			   const gretl_matrix *lfrac,
			   gretl_matrix *XVC,
			   double lmax, double rho,
			   int fold, int crit_type)
{
    static gretl_vector *x, *u, *z, *y;
    static gretl_vector *r, *zprev, *zdiff;
    static gretl_vector *q, *Atb, *m1, *L;
    static gretl_matrix_block *MB;
    int ldim, nlam;
    int m, n, j;
    int err = 0;

    if (A == NULL) {
	/* cleanup signal */
	gretl_matrix_block_destroy(MB);
	MB = NULL;
	return 0;
    }

    /* dimensions */
    nlam = gretl_vector_get_length(lfrac);
    m = A->rows;
    n = A->cols;
    ldim = m >= n ? n : m;

    if (MB == NULL) {
	MB = gretl_matrix_block_new(&x, n, 1, &u, n, 1,
				    &z, n, 1, &y, n, 1,
				    &r, n, 1, &zprev, n, 1,
				    &zdiff, n, 1, &q, n, 1,
				    &m1, m, 1, &Atb, n, 1,
				    &L, ldim, ldim, NULL);
	if (MB == NULL) {
	    return E_ALLOC;
	}
	gretl_matrix_block_zero(MB);
    }

    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);
#if 0
    lmax = gretl_matrix_infinity_norm(Atb);
#endif

    get_cholesky_factor(A, L, rho);

    for (j=0; j<nlam && !err; j++) {
	/* loop across lambda values */
	double score, lambda = lfrac->val[j] * lmax;
	int iters = 0;

	err = admm_iteration(A, L, Atb, x, z, u, q, m1, r, zprev, zdiff,
			     lambda, rho, &iters);

	if (!err) {
	    /* record out-of-sample criterion */
	    gretl_matrix_reuse(m1, A_out->rows, 1);
	    score = xv_score(A_out, b_out, z, m1, crit_type);
	    gretl_matrix_reuse(m1, m, 1);
	    gretl_matrix_set(XVC, j, fold, score);
	}
    }

    return err;
}

static void prepare_xv_data (const gretl_matrix *X,
			     const gretl_matrix *y,
			     gretl_matrix *Ae,
			     gretl_matrix *be,
			     gretl_matrix *Af,
			     gretl_matrix *bf,
			     int f)
{
    int i, j, ke, ko;
    double xij;

    for (j=0; j<X->cols; j++) {
	ke = ko = 0;
	for (i=0; i<X->rows; i++) {
	    xij = gretl_matrix_get(X, i, j);
	    if (i/Af->rows == f) {
		/* "out of sample" range */
		gretl_matrix_set(Af, ko, j, xij);
		if (j == 0) {
		    bf->val[ko] = y->val[i];
		}
		ko++;
	    } else {
		/* estimation sample */
		gretl_matrix_set(Ae, ke, j, xij);
		if (j == 0) {
		    be->val[ke] = y->val[i];
		}
		ke++;
	    }
	}
    }
}

/* Given @C holding criterion values per lambda (rows) and
   per fold (columns), write into the last two columns
   the means and standard errors.
*/

static void process_xv_criterion (gretl_matrix *XVC,
				  gretl_matrix *lfrac,
				  int *ibest, int *i1se,
				  int nf, int crit_type)
{
    double avg, d, v, se, se1, avgmin = 1e200;
    int i, j, ialt, imin = 0;

    for (i=0; i<XVC->rows; i++) {
	v = avg = 0;
	for (j=0; j<nf; j++) {
	    avg += gretl_matrix_get(XVC, i, j);
	}
	avg /= nf;
	if (i == 0) {
	    avgmin = avg;
	} else if (avg < avgmin) {
	    avgmin = avg;
	    imin = i;
	}
	gretl_matrix_set(XVC, i, nf, avg);
	if (crit_type == CRIT_PCC) {
	    printf("s = %#g -> %s %#g\n", lfrac->val[i],
		   crit_string(crit_type), 100 - avg);
	    continue;
	}
	for (j=0; j<nf; j++) {
	    d = gretl_matrix_get(XVC, i, j) - avg;
	    v += d * d;
	}
	v /= (nf - 1);
	se = sqrt(v/nf);
	gretl_matrix_set(XVC, i, nf+1, se);
	printf("s = %#g -> %s %#g (%#g)\n", lfrac->val[i],
	       crit_string(crit_type), avg, se);
    }

    *ibest = ialt = imin;

    if (crit_type == CRIT_PCC) {
	return;
    }

    /* estd. standard error of minimum average XVC */
    se1 = gretl_matrix_get(XVC, imin, nf+1);

    /* Find the index of the largest lamba that gives
       an average XVC within one standard error of the
       minimum (glmnet's "$lambda.1se").
    */
    for (i=imin-1; i>=0; i--) {
	avg = gretl_matrix_get(XVC, i, nf);
	if (avg - avgmin < se1) {
	    ialt = i;
	} else {
	    break;
	}
    }

    *i1se = ialt;
}

/* Shrink the cross validation criterion matrix down to what the user
   might be interested in: the per-lambda criteria and their estimated
   std errors.
*/

static void shrink_crit (gretl_matrix **pXVC)
{
    gretl_matrix *XVC = *pXVC;
    gretl_matrix *tmp;

    tmp = gretl_matrix_alloc(XVC->rows, 2);

    if (tmp != NULL) {
	int n = XVC->rows;
	int offset = n * (XVC->cols - 2);

	memcpy(tmp->val, XVC->val + offset, 2*n * sizeof *tmp->val);
	gretl_matrix_free(XVC);
	*pXVC = tmp;
    }
}

static int get_crit_type (gretl_bundle *bun)
{
    const char *s = gretl_bundle_get_string(bun, "xvcrit", NULL);
    int ret = 0;

    if (s != NULL) {
	if (g_ascii_strcasecmp(s, "mse") == 0) {
	    ret = CRIT_MSE;
	} else if (g_ascii_strcasecmp(s, "mae") == 0) {
	    ret = CRIT_MAE;
	} else if (g_ascii_strcasecmp(s, "rank") == 0) {
	    ret = CRIT_PCC;
	} else {
	    gretl_errmsg_sprintf("'%s' invalid criterion", s);
	    ret = -1;
	}
    }

    return ret;
}

static int admm_lasso_xv (gretl_matrix *A,
			  gretl_matrix *b,
			  gretl_bundle *bun,
			  double rho)
{
    gretl_matrix_block *AB;
    gretl_matrix *Ae, *Af;
    gretl_matrix *be, *bf;
    gretl_matrix *lfrac;
    gretl_matrix *Atb, *XVC;
    double lmax;
    int nlam, fsize, esize;
    int randfolds = 0;
    int crit_type = 0;
    int f, nf;
    int err = 0;

    nf = gretl_bundle_get_int(bun, "nfolds", &err);
    randfolds = gretl_bundle_get_int(bun, "randfolds", &err);
    lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);
    if (err) {
	return err;
    }

    crit_type = get_crit_type(bun);
    if (crit_type < 0) {
	return E_INVARG;
    }

    fsize = A->rows / nf;
    esize = (nf - 1) * fsize;

    printf("admm_lasso_xv: nf=%d, fsize=%d, randfolds=%d, crit=%s\n",
	   nf, fsize, randfolds, crit_string(crit_type));

    AB = gretl_matrix_block_new(&Ae, esize, A->cols,
				&Af, fsize, A->cols,
				&be, esize, 1,
				&bf, fsize, 1, NULL);
    if (AB == NULL) {
	return E_ALLOC;
    }

    nlam = gretl_vector_get_length(lfrac);

    /* determine the infnorm for all training data */
    Atb = gretl_matrix_alloc(A->cols, 1);
    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);
    lmax = gretl_matrix_infinity_norm(Atb);
    /* and scale it down for the folds */
    lmax *= esize / (double) A->rows;
    fprintf(stderr, "cross validation lambda max = %#g\n", lmax);
    gretl_matrix_free(Atb);

    if (randfolds) {
	/* scramble the row order of A and b */
	randomize_rows(A, b);
    }

    XVC = gretl_zero_matrix_new(nlam, nf + 2);

    for (f=0; f<nf && !err; f++) {
	prepare_xv_data(A, b, Ae, be, Af, bf, f);
	err = lasso_xv_round(Ae, be, Af, bf, lfrac, XVC, lmax, rho,
			     f, crit_type);
    }

    /* send cleanup signal */
    lasso_xv_round(NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0);

    if (!err) {
	gretl_matrix *lxv = gretl_matrix_alloc(1, 1);
	int ibest = 0, i1se = 0;

	process_xv_criterion(XVC, lfrac, &ibest, &i1se, nf, crit_type);
	printf("\nAverage out-of-sample %s minimized at %g for s=%g\n",
	       crit_string(crit_type), gretl_matrix_get(XVC, ibest, nf),
	       lfrac->val[ibest]);
	printf("Largest s within one s.e. of best criterion: %g\n",
	       lfrac->val[i1se]);
	/* now determine coefficient vector on full training set */
	lxv->val[0] = lfrac->val[ibest];
	gretl_bundle_donate_data(bun, "lxv", lxv, GRETL_TYPE_MATRIX, 0);
	err = real_admm_lasso(A, b, bun, rho);
    }

    if (!err) {
	shrink_crit(&XVC);
	gretl_bundle_donate_data(bun, "XVC", XVC, GRETL_TYPE_MATRIX, 0);
    } else {
	gretl_matrix_free(XVC);
    }

    gretl_matrix_block_destroy(AB);

    return err;
}

int admm_lasso (gretl_matrix *A,
		gretl_matrix *b,
		gretl_bundle *bun)
{
    gretl_matrix *ctrl;
    double rho = 8.0; /* was 1.0 */
    int xv, err = 0;

    ctrl = gretl_bundle_get_matrix(bun, "admmctrl", NULL);
    if (ctrl != NULL) {
	if (ctrl->val[0] > 0) {
	    rho = ctrl->val[0];
	}
	if (ctrl->val[1] > 0) {
	    reltol = ctrl->val[1];
	}
	if (ctrl->val[2] > 0) {
	    abstol = ctrl->val[2];
	}
    }

#if USE_YBAR
    if (gretl_bundle_get_int(bun, "stdize_y", NULL) == 0) {
	/* we'll need to add mean(y) */
	ybar = gretl_mean(0, b->rows-1, b->val);
    }
#endif

    /* scale the absolute tolerance */
    abstol = sqrt(A->cols) * abstol;

    xv = gretl_bundle_get_int(bun, "xvalidate", &err);

    if (xv) {
	return admm_lasso_xv(A, b, bun, rho);
    } else {
	return real_admm_lasso(A, b, bun, rho);
    }
}
