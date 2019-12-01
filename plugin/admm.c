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

#ifdef HAVE_MPI
# include "gretl_mpi.h"
# include "gretl_foreign.h"
#endif

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

#define RELTOL_DEFAULT 1.0e-4
#define ABSTOL_DEFAULT 1.0e-6

double reltol;
double abstol;
double ybar;

enum {
    CRIT_MSE,
    CRIT_MAE,
    CRIT_PCC
};

#ifdef HAVE_MPI
static int mpi_parent_action (gretl_matrix *A,
			      gretl_matrix *b,
			      gretl_bundle *bun,
			      double rho,
			      PRN *prn);
#endif

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
    if (ybar != 0) {
	vector_add_scalar(Azb, ybar, A->rows);
    }
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
    if (ybar != 0) {
	vector_add_scalar(Azb, ybar, A->rows);
    }
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

#define VAR_RHO 1 /* seems to work well */
#define RHO_DEBUG 0

static int admm_iteration (const gretl_matrix *A,
			   const gretl_vector *Atb,
			   gretl_matrix *L,
			   gretl_vector *x, gretl_vector *z,
			   gretl_vector *u, gretl_vector *q,
			   gretl_vector *p, gretl_vector *r,
			   gretl_vector *zprev, gretl_vector *zdiff,
			   double lambda, double *prho,
			   int tune_rho, int *iters)
{
    double nxstack, nystack;
    double prires, dualres;
    double eps_pri, eps_dual;
    double rho = *prho;
    double nrm2, rho2 = rho*rho;
    int itermin = 1;
    int n = A->cols;
    int iter = 0;
    int err = 0;

#if RHO_DEBUG
    fprintf(stderr, "*** admm: lambda %g, rho %g ***\n", lambda, rho);
#endif

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
	    gretl_matrix_multiply_by_scalar(x, -1/rho2);
	    gretl_matrix_multiply_by_scalar(q, 1/rho);
	    vector_add_to(x, q, n);
	}

	/* sqrt(sum ||r_i||_2^2) */
	prires  = sqrt(gretl_vector_dot_product(r, r, NULL));
	/* sqrt(sum ||r_i||_2^2) */
	nxstack = sqrt(gretl_vector_dot_product(x, x, NULL));
	/* sqrt(sum ||y_i||_2^2) */
	nystack = gretl_vector_dot_product(u, u, NULL) / rho2;
	nystack = sqrt(nystack);

	vector_copy_values(zprev, z, n);
	vector_add_into(x, u, z, n);
	soft_threshold(z, lambda, rho);

	/* Termination checks */

	/* dual residual */
	vector_subtract_into(z, zprev, zdiff, n, 0); /* zdiff = z - zprev */
	/* ||s^k||_2^2 = N rho^2 ||z - zprev||_2^2 */
	nrm2 = sqrt(gretl_vector_dot_product(zdiff, zdiff, NULL));
	dualres = rho * nrm2;

	/* compute primal and dual feasibility tolerances */
	nrm2 = sqrt(gretl_vector_dot_product(z, z, NULL));
	eps_pri  = abstol + reltol * fmax(nxstack, nrm2);
	eps_dual = abstol + reltol * nystack;

	if (iter >= itermin && prires <= eps_pri && dualres <= eps_dual) {
	    break;
	}

	/* Compute residual: r = x - z */
	vector_subtract_into(x, z, r, n, 0);

	if (tune_rho && iter > 0 && (iter == 32 || iter % 200 == 0)) {
	    double mult = 10;
	    double adj = 0.0;

	    if (prires > mult * dualres) {
		adj = 2.0;
	    } else if (dualres > mult * prires) {
		adj = 0.5;
	    }
	    if (adj > 0) {
		rho *= adj;
# if RHO_DEBUG
		fprintf(stderr, "  iter %d: rho *= %g (now %g)\n",
			iter, adj, rho);
# endif
		rho2 = rho * rho;
		gretl_matrix_multiply_by_scalar(u, 1.0/adj);
		gretl_matrix_multiply_by_scalar(r, 1.0/adj);
		get_cholesky_factor(A, L, rho);
		/* ensure a fair number of subsequent iterations */
		itermin = iter + 100;
	    }
	}

	iter++;
    }

    *iters = iter;
    *prho = rho;

    return err;
}

static gretl_matrix *make_coeff_matrix (gretl_bundle *bun, int xvalid,
					int rows, int nlam,
					int *jmin, int *jmax)
{
    gretl_matrix *B = NULL;
    int xv_single_b = 0;

    if (xvalid) {
	/* do we want just the "best" coeff vector? */
	xv_single_b = gretl_bundle_get_bool(bun, "single_b", 0);
    }

    if (xv_single_b) {
	int use_1se = gretl_bundle_get_bool(bun, "use_1se", 0);
	const char *ikey = use_1se ? "idx1se" : "idxmin";
	int idx;

	idx = gretl_bundle_get_int(bun, ikey, NULL);
	B = gretl_zero_matrix_new(rows, 1);
	*jmin = idx - 1; /* zero-based */
	*jmax = *jmin + 1;
    } else {
	B = gretl_zero_matrix_new(rows, nlam);
	*jmin = 0;
	*jmax = nlam;
    }

    if (B != NULL) {
	const char *bkey = B->cols == 1 ? "b" : "B";

	gretl_bundle_donate_data(bun, bkey, B, GRETL_TYPE_MATRIX, 0);
    }

    return B;
}

/* This function is executed when we want to obtain a set
   of coefficients using the full training data, with either
   a single value of lambda or a vector of lambdas. We come
   here straight away if the user has not requested cross
   validation; we also come here after cross validation.
*/

static int real_admm_lasso (const gretl_matrix *A,
			    const gretl_matrix *b,
			    gretl_bundle *bun,
			    double rho0, PRN *prn)
{
    gretl_matrix_block *MB;
    double lcritmin = 1e200;
    gretl_matrix *B = NULL;
    gretl_matrix *lfrac;
    double lmax, rho = rho0;
    int ldim, nlam = 1;
    int m, n, i, j;
    int jmin, jmax;
    int idxmin = 0;
    int err = 0;

    gretl_vector *x, *u, *z, *y, *r, *zprev, *zdiff;
    gretl_vector *q, *Atb, *m1;
    gretl_matrix *L;

    int stdize = gretl_bundle_get_int(bun, "stdize", &err);
    int xvalid = gretl_bundle_get_int(bun, "xvalidate", &err);
    int verbose = gretl_bundle_get_bool(bun, "verbosity", 1);

    lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);

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

    if (!xvalid && verbose > 0) {
	if (nlam > 1) {
	    pprintf(prn, "using lambda-fraction sequence of length %d, starting at %g\n",
		    nlam, lfrac->val[0]);
	} else {
	    pprintf(prn, "using lambda-fraction %g\n", lfrac->val[0]);
	}
    }

    get_cholesky_factor(A, L, rho);

    B = make_coeff_matrix(bun, xvalid, n + stdize, nlam, &jmin, &jmax);
    if (B == NULL) {
	gretl_matrix_block_destroy(MB);
	return E_ALLOC;
    }

    if (!xvalid && verbose > 0 && nlam > 1) {
	pputc(prn, '\n');
	pprintf(prn, "      lambda     df     criterion   iters\n");
    }

    for (j=jmin; j<jmax && !err; j++) {
	/* loop across lambda values */
	double lcrit, lambda = lfrac->val[j] * lmax;
	int tune_rho = 0;
	int iters = 0;
	int nnz = 0;

#if VAR_RHO
	tune_rho = 1;
#endif

	err = admm_iteration(A, Atb, L, x, z, u, q, m1, r, zprev, zdiff,
			     lambda, &rho, tune_rho, &iters);

	if (!err) {
	    for (i=0; i<n; i++) {
		if (z->val[i] != 0.0) {
		    nnz++;
		}
		if (B->cols == 1) {
		    gretl_matrix_set(B, i+stdize, 0, z->val[i]);
		} else {
		    gretl_matrix_set(B, i+stdize, j, z->val[i]);
		}
	    }
	    if (!xvalid) {
		lcrit = objective(A, b, z, lambda, m1);
		if (verbose > 0 && nlam > 1) {
		    pprintf(prn, "%#12.6g  %5d    %#.8g   %5d\n",
			    lambda/m, nnz, lcrit, iters);
		}
		if (lcrit < lcritmin) {
		    lcritmin = lcrit;
		    idxmin = j;
		}
	    }
	}
    }

    gretl_bundle_set_scalar(bun, "lmax", lmax);
    if (!xvalid) {
	if (nlam > 1) {
	    gretl_bundle_set_scalar(bun, "idxmin", idxmin + 1);
	    gretl_bundle_set_scalar(bun, "lfmin", lfrac->val[idxmin]);
	}
	gretl_bundle_set_scalar(bun, "lcrit", lcritmin);
    }
    if (nlam == 1) {
	gretl_bundle_set_scalar(bun, "lambda", lfrac->val[0] * lmax);
    }

    gretl_bundle_delete_data(bun, "verbosity");

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
			   double lmax, double rho0,
			   int fold, int crit_type)
{
    static gretl_vector *x, *u, *z, *y;
    static gretl_vector *r, *zprev, *zdiff;
    static gretl_vector *q, *Atb, *m1, *L;
    static gretl_matrix_block *MB;
    double rho = rho0;
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
#if 0 /* ?? */
    lmax = gretl_matrix_infinity_norm(Atb);
#endif

    get_cholesky_factor(A, L, rho);

    for (j=0; j<nlam && !err; j++) {
	/* loop across lambda values */
	double score, lambda = lfrac->val[j] * lmax;
	int tune_rho = 0;
	int iters = 0;

#if VAR_RHO
	tune_rho = 1;
#endif
	err = admm_iteration(A, Atb, L, x, z, u, q, m1, r, zprev, zdiff,
			     lambda, &rho, tune_rho, &iters);

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
			     gretl_matrix *Xe,
			     gretl_matrix *ye,
			     gretl_matrix *Xf,
			     gretl_matrix *yf,
			     int f)
{
    int i, j, ke, ko;
    double xij;

    for (j=0; j<X->cols; j++) {
	ke = ko = 0;
	for (i=0; i<X->rows; i++) {
	    xij = gretl_matrix_get(X, i, j);
	    if (i/Xf->rows == f) {
		/* "out of sample" range */
		gretl_matrix_set(Xf, ko, j, xij);
		if (j == 0) {
		    yf->val[ko] = y->val[i];
		}
		ko++;
	    } else {
		/* estimation sample */
		gretl_matrix_set(Xe, ke, j, xij);
		if (j == 0) {
		    ye->val[ke] = y->val[i];
		}
		ke++;
	    }
	}
    }
}

/* Given @C holding criterion values per lambda (rows) and
   per fold (columns), compose a matrix holding the means
   and possibly standard errors.
*/

static gretl_matrix *process_xv_criterion (gretl_matrix *XVC,
					   gretl_matrix *lfrac,
					   int *imin, int *i1se,
					   int crit_type,
					   PRN *prn)
{
    gretl_matrix *metrics;
    double avg, d, v, se, se1, avgmin = 1e200;
    int mcols = (crit_type == CRIT_PCC)? 1 : 2;
    int nf = XVC->cols;
    int i, j;

    metrics = gretl_zero_matrix_new(XVC->rows, mcols);
    if (metrics == NULL) {
	return NULL;
    }

    *imin = 0;

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
	    *imin = i;
	}
	gretl_matrix_set(metrics, i, 0, avg);
	if (crit_type == CRIT_PCC && prn != NULL) {
	    pprintf(prn, "s = %#g -> %s %#g\n", lfrac->val[i],
		    crit_string(crit_type), 100 - avg);
	    continue;
	}
	for (j=0; j<nf; j++) {
	    d = gretl_matrix_get(XVC, i, j) - avg;
	    v += d * d;
	}
	v /= (nf - 1);
	se = sqrt(v/nf);
	gretl_matrix_set(metrics, i, 1, se);
	if (prn != NULL) {
	    pprintf(prn, "s = %#g -> %s %#g (%#g)\n", lfrac->val[i],
		    crit_string(crit_type), avg, se);
	}
    }

    *i1se = *imin;

    if (crit_type != CRIT_PCC) {
	/* estd. standard error of minimum average XVC */
	se1 = gretl_matrix_get(metrics, *imin, 1);

	/* Find the index of the largest lamba that gives
	   an average XVC within one standard error of the
	   minimum (glmnet's "$lambda.1se").
	*/
	for (i=*imin-1; i>=0; i--) {
	    avg = gretl_matrix_get(metrics, i, 0);
	    if (avg - avgmin < se1) {
		*i1se = i;
	    } else {
		break;
	    }
	}
    }

    return metrics;
}

/* Analyse results after cross-validation. Return the
   optimal lambda value or NADBL on failure.
*/

static int post_xvalidation_task (gretl_matrix *XVC,
				  gretl_matrix *lfrac,
				  int crit_type,
				  gretl_bundle *b,
				  PRN *prn)
{
    gretl_matrix *metrics;
    int imin = 0, i1se = 0;

    metrics = process_xv_criterion(XVC, lfrac, &imin, &i1se,
				   crit_type, prn);
    if (metrics == NULL) {
	return E_ALLOC;
    }

    if (prn != NULL) {
	pprintf(prn, "\nAverage out-of-sample %s minimized at %g for s=%g\n",
		crit_string(crit_type), gretl_matrix_get(metrics, imin, 0),
		lfrac->val[imin]);
	pprintf(prn, "Largest s within one s.e. of best criterion: %g\n",
		lfrac->val[i1se]);
    }

    gretl_bundle_donate_data(b, "XVC", metrics, GRETL_TYPE_MATRIX, 0);
    gretl_bundle_set_int(b, "idxmin", imin + 1);
    gretl_bundle_set_int(b, "idx1se", i1se + 1);
    gretl_bundle_set_scalar(b, "lfmin", lfrac->val[imin]);
    gretl_bundle_set_scalar(b, "lf1se", lfrac->val[i1se]);

    return 0;
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

static double get_xvalidation_lmax (gretl_matrix *A,
				    gretl_matrix *b,
				    int esize)
{
    gretl_matrix *Atb;
    double lmax;

    /* determine the infnorm for all training data */
    Atb = gretl_matrix_alloc(A->cols, 1);
    gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
			      b, GRETL_MOD_NONE,
			      Atb, GRETL_MOD_NONE);
    lmax = gretl_matrix_infinity_norm(Atb);

    /* and scale it down for the folds */
    lmax *= esize / (double) A->rows;

    gretl_matrix_free(Atb);

    return lmax;
}

static int get_xvalidation_details (gretl_bundle *bun,
				    int *nf, int *randfolds,
				    gretl_matrix **lfrac,
				    int *crit_type)
{
    int err = 0;

    *nf = gretl_bundle_get_int(bun, "nfolds", &err);
    *randfolds = gretl_bundle_get_int(bun, "randfolds", &err);
    *lfrac = gretl_bundle_get_matrix(bun, "lfrac", &err);

    if (!err && *nf < 2) {
	err = E_INVARG;
    }
    if (!err) {
	*crit_type = get_crit_type(bun);
	if (*crit_type < 0) {
	    err = E_INVARG;
	}
    }

    return err;
}

static int admm_lasso_xv (gretl_matrix *A,
			  gretl_matrix *b,
			  gretl_bundle *bun,
			  double rho,
			  PRN *prn)
{
    gretl_matrix_block *AB;
    gretl_matrix *Ae, *Af;
    gretl_matrix *be, *bf;
    gretl_matrix *lfrac;
    gretl_matrix *XVC;
    double lmax;
    int nlam, fsize, esize;
    int randfolds = 0;
    int crit_type = 0;
    int verbose;
    int f, nf;
    int err;

    err = get_xvalidation_details(bun, &nf, &randfolds,
				  &lfrac, &crit_type);
    if (err) {
	return err;
    }

    verbose = gretl_bundle_get_bool(bun, "verbosity", 1);

    fsize = A->rows / nf;
    esize = (nf - 1) * fsize;

    if (verbose) {
	pprintf(prn, "admm_lasso_xv: nf=%d, fsize=%d, randfolds=%d, crit=%s\n\n",
		nf, fsize, randfolds, crit_string(crit_type));
    }
    manufacture_gui_callback(FLUSH);

    AB = gretl_matrix_block_new(&Ae, esize, A->cols,
				&Af, fsize, A->cols,
				&be, esize, 1,
				&bf, fsize, 1, NULL);
    if (AB == NULL) {
	return E_ALLOC;
    }

    nlam = gretl_vector_get_length(lfrac);
    lmax = get_xvalidation_lmax(A, b, esize);

    if (randfolds) {
	/* scramble the row order of A and b */
	randomize_rows(A, b);
    }

    XVC = gretl_zero_matrix_new(nlam, nf);

    for (f=0; f<nf && !err; f++) {
	prepare_xv_data(A, b, Ae, be, Af, bf, f);
	err = lasso_xv_round(Ae, be, Af, bf, lfrac, XVC, lmax, rho,
			     f, crit_type);
    }

    /* send cleanup signal */
    lasso_xv_round(NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0);

    if (!err) {
	PRN *myprn = verbose ? prn : NULL;

	err = post_xvalidation_task(XVC, lfrac, crit_type, bun, myprn);
	if (!err) {
	    /* determine coefficient vector on full training set */
	    err = real_admm_lasso(A, b, bun, rho, myprn);
	}
    }

    gretl_matrix_free(XVC);
    gretl_matrix_block_destroy(AB);

    return err;
}

#ifdef HAVE_MPI

static int mpi_admm_lasso_xv (gretl_matrix *A,
			      gretl_matrix *b,
			      gretl_bundle *bun,
			      double rho,
			      PRN *prn)
{
    gretl_matrix_block *AB = NULL;
    gretl_matrix *XVC = NULL;
    gretl_matrix *Ae = NULL;
    gretl_matrix *Af = NULL;
    gretl_matrix *be = NULL;
    gretl_matrix *bf = NULL;
    gretl_matrix *lfrac;
    double lmax;
    int fsize, esize;
    int folds_per;
    int folds_rem;
    int randfolds;
    int nlam, rank;
    int crit_type = 0;
    int np, rankmax = 0;
    int verbose;
    int f, nf, r;
    int my_f = 0;
    int err = 0;

    rank = gretl_mpi_rank();
    np = gretl_mpi_n_processes();
    rankmax = np - 1;

    err = get_xvalidation_details(bun, &nf, &randfolds,
				  &lfrac, &crit_type);
    if (err) {
	return err;
    }

    verbose = gretl_bundle_get_int_deflt(bun, "verbosity", 1);

    nlam = gretl_vector_get_length(lfrac);
    fsize = A->rows / nf;
    esize = (nf - 1) * fsize;
    folds_per = nf / np;
    folds_rem = nf % np;

    /* matrix-space for per-fold data */
    AB = gretl_matrix_block_new(&Ae, esize, A->cols,
				&Af, fsize, A->cols,
				&be, esize, 1,
				&bf, fsize, 1, NULL);
    if (AB == NULL) {
	return E_ALLOC;
    }

    if (rank == 0) {
	lmax = get_xvalidation_lmax(A, b, esize);
    }

    if (randfolds) {
	/* generate the same random folds in all processes */
	unsigned seed;

	if (rank == 0) {
	    if (gretl_bundle_has_key(bun, "seed")) {
		seed = gretl_bundle_get_unsigned(bun, "seed", NULL);
	    } else {
		seed = gretl_rand_get_seed();
	    }
	}
	gretl_mpi_bcast(&seed, GRETL_TYPE_UNSIGNED, 0);
	gretl_rand_set_seed(seed);
	randomize_rows(A, b);
    }

    if (rank < folds_rem) {
	XVC = gretl_zero_matrix_new(nlam, folds_per + 1);
    } else {
	XVC = gretl_zero_matrix_new(nlam, folds_per);
    }

    /* send @lmax to workers */
    gretl_mpi_bcast(&lmax, GRETL_TYPE_DOUBLE, 0);

    if (rank == 0) {
	if (verbose) {
	    pprintf(prn, "admm_lasso_xv: nf=%d, fsize=%d, randfolds=%d, crit=%s\n\n",
		    nf, fsize, randfolds, crit_string(crit_type));
	}
	manufacture_gui_callback(FLUSH);
    }

    /* process all folds */
    r = 0;
    for (f=0; f<nf && !err; f++) {
	if (rank == r) {
	    prepare_xv_data(A, b, Ae, be, Af, bf, f);
	    if (verbose > 1) {
		pprintf(prn, "rank %d: taking fold %d\n", rank, f+1);
	    }
	    err = lasso_xv_round(Ae, be, Af, bf, lfrac, XVC, lmax, rho,
				 my_f++, crit_type);
	}
	if (r == rankmax) {
	    r = 0;
	} else {
	    r++;
	}
    }

    /* reduce @XVC to root by column concatenation */
    gretl_matrix_mpi_reduce(XVC, &XVC, GRETL_MPI_HCAT, 0, OPT_NONE);

    /* send cleanup signal, all processes */
    lasso_xv_round(NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0);

    if (rank == 0 && !err) {
	PRN *myprn = verbose ? prn : NULL;

	err = post_xvalidation_task(XVC, lfrac, crit_type, bun, myprn);
	if (!err) {
	    /* determine coefficient vector on full training set */
	    err = real_admm_lasso(A, b, bun, rho, myprn);
	}
    }

    gretl_matrix_free(XVC);
    gretl_matrix_block_destroy(AB);

    return err;
}

#endif /* HAVE_MPI */

static void prepare_admm_params (gretl_matrix *A,
				 gretl_matrix *b,
				 gretl_bundle *bun,
				 double *rho)
{
    gretl_matrix *ctrl;
    int len;

    /* set defaults */
    reltol = RELTOL_DEFAULT;
    abstol = ABSTOL_DEFAULT;

    ctrl = gretl_bundle_get_matrix(bun, "admmctrl", NULL);
    len = gretl_vector_get_length(ctrl);

    if (len > 0 && ctrl->val[0] > 0) {
	*rho = ctrl->val[0];
    }
    if (len > 1 && ctrl->val[1] > 0) {
	reltol = ctrl->val[1];
    }
    if (len > 2 && ctrl->val[2] > 0) {
	abstol = ctrl->val[2];
    }

    if (gretl_bundle_get_bool(bun, "stdize_y", 0)) {
	/* we need to add mean(y) */
	ybar = gretl_mean(0, b->rows-1, b->val);
    } else {
	ybar = 0.0;
    }

    /* scale the absolute tolerance */
    abstol = sqrt(A->cols) * abstol;
}

int admm_lasso (gretl_matrix *A,
		gretl_matrix *b,
		gretl_bundle *bun,
		PRN *prn)
{
    double rho = 8.0;
    int xv;

    prepare_admm_params(A, b, bun, &rho);

    xv = gretl_bundle_get_bool(bun, "xvalidate", 0);

    if (xv) {
#ifdef HAVE_MPI
	int no_mpi = gretl_bundle_get_bool(bun, "no_mpi", 0);

	if (!no_mpi) {
	    if (gretl_mpi_n_processes() > 1) {
		return mpi_admm_lasso_xv(A, b, bun, rho, prn);
	    } else if (auto_mpi_ok()) {
		return mpi_parent_action(A, b, bun, rho, prn);
	    }
	}
#endif
	return admm_lasso_xv(A, b, bun, rho, prn);
    } else {
	return real_admm_lasso(A, b, bun, rho, prn);
    }
}

#ifdef HAVE_MPI

/* We come here if a parent process has called our
   automatic local MPI routine for cross validation:
   this function will be executed by all gretlmpi
   instances.
*/

/* Using shared memory to transfer matrices is perhaps
   a little faster, but it's not yet tested on Windows.
*/
#define MPI_USE_SHM 0

int admm_xv_mpi (PRN *prn)
{
    gretl_bundle *bun = NULL;
    gretl_matrix *A;
    gretl_matrix *b;
    double rho = 8.0;
    int rank, err = 0;

    rank = gretl_mpi_rank();

    /* read matrices deposited by parent process */
#if MPI_USE_SHM
    A = shm_read_matrix("lasso_A.shm", 0, &err);
    b = shm_read_matrix("lasso_b.shm", 0, &err);
#else
    A = gretl_matrix_read_from_file("lasso_A.bin", 1, &err);
    b = gretl_matrix_read_from_file("lasso_b.bin", 1, &err);
#endif

    if (!err) {
	bun = gretl_bundle_read_from_file("lasso_bun.xml", 1, &err);
    }

    if (!err) {
	prepare_admm_params(A, b, bun, &rho);
    }

    if (!err) {
	err = mpi_admm_lasso_xv(A, b, bun, rho, prn);
	if (!err && rank == 0) {
	    /* write results, to be picked up by parent */
	    gretl_bundle_write_to_file(bun, "lasso_XV_result.xml", 1);
	}
    }

#if MPI_USE_SHM
    if (rank == 0) {
	shm_finalize_matrix("lasso_A.shm");
	shm_finalize_matrix("lasso_b.shm");
    }
#endif

    gretl_matrix_free(A);
    gretl_matrix_free(b);
    gretl_bundle_destroy(bun);

    return err;
}

static int mpi_parent_action (gretl_matrix *A,
			      gretl_matrix *b,
			      gretl_bundle *bun,
			      double rho,
			      PRN *prn)
{
    int err;

#if MPI_USE_SHM
    err = shm_write_matrix(A, "lasso_A.shm");
    if (!err) {
	err = shm_write_matrix(b, "lasso_b.shm");
    }
#else
    err = gretl_matrix_write_to_file(A, "lasso_A.bin", 1);
    if (!err) {
	err = gretl_matrix_write_to_file(b, "lasso_b.bin", 1);
    }
#endif

    if (!err) {
	err = gretl_bundle_write_to_file(bun, "lasso_bun.xml", 1);
    }

    if (!err) {
	/* compose and execute MPI script */
	err = foreign_start(MPI, NULL, OPT_NONE, prn);
	if (!err) {
	    foreign_append("_admm_lasso()", MPI);
	    err = foreign_execute(NULL, OPT_L | OPT_S | OPT_Q, prn);
	    if (err) {
		fprintf(stderr, "mpi_parent: foreign exec error %d\n", err);
	    }
	}
    }

    if (!err) {
	/* retrieve results bundle written by gretlmpi */
	gretl_bundle *res;

	res = gretl_bundle_read_from_file("lasso_XV_result.xml", 1, &err);
	if (!err) {
	    gretl_bundles_swap_content(bun, res);
	    gretl_bundle_destroy(res);
	    gretl_bundle_delete_data(bun, "verbosity");
	}
    }

    return err;
}

#endif /* HAVE_MPI */
