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

/* Code for regularized least squares. Includes these methods:

   ADMM: Based on Boyd et al, "Distributed Optimization and
   Statistical Learning via the Alternating Direction Method of
   Multipliers", Foundations and Trends in Machine Learning, Vol. 3,
   No. 1 (2010) 1-122.

   CCD (Cyclical Coordinate Descent): Based on the Fortran code
   employed by R's glmnet for the Gaussian case and the "covariance"
   algorithm.

   SVD: for Ridge.
*/

#include "libgretl.h"
#include "matrix_extra.h"
#include "version.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
# include "gretl_foreign.h"
#endif

#if defined(HAVE_IMMINTRIN_H) && defined(USE_AVX)
# define USE_SIMD
# include <immintrin.h>
#endif

#define ADMM_MAX_ITER 20000
#define ADMM_RELTOL_DEFAULT 1.0e-4
#define ADMM_ABSTOL_DEFAULT 1.0e-6

double admm_reltol;
double admm_abstol;

#define CCD_MAX_ITER 100000
#define CCD_TOLER_DEFAULT 1.0e-7
#define BIG_LAMBDA 9.9e35

#define RHO_DEBUG 0
#define LAMBDA_DEBUG 0
#define RIDGE_DEBUG 0

double ccd_toler;

enum {
    LAMSCALE_NONE,
    LAMSCALE_GLMNET,
    LAMSCALE_FROB
};

typedef struct regls_info_ {
    gretl_bundle *b;
    gretl_matrix *X;
    gretl_matrix *y;
    gretl_matrix *lfrac;
    gretl_matrix *Xty;
    gretl_matrix *R2;
    gretl_matrix *crit;
    gretl_matrix *BIC;
    gretl_matrix *edf;
    double rho;
    double infnorm;
    double alpha;
    int nlam;
    int n;
    int k;
    int nf;
    gint8 ccd;
    gint8 ridge;
    gint8 stdize;
    gint8 xvalid;
    gint8 verbose;
    gint8 lamscale;
    gint8 randfolds;
    gint8 use_1se;
    gint8 free_lf;
    PRN *prn;
} regls_info;

typedef struct ccd_info_ {
    gretl_matrix_block *MB;
    gretl_matrix *Xty;
    gretl_matrix *xv;
    gretl_matrix *B;
    gretl_matrix *lam;
    double lmax;
} ccd_info;

#ifdef HAVE_MPI
static int mpi_parent_action (regls_info *ri);
#endif

static void prepare_ccd_param (regls_info *ri)
{
    double tol;

    tol = gretl_bundle_get_scalar(ri->b, "ccd_toler", NULL);

    if (!na(tol) && tol > 0.0 && tol < 1.0) {
	ccd_toler = tol;
    } else {
	ccd_toler = CCD_TOLER_DEFAULT;
    }
}

static void maybe_set_lambda_scale (regls_info *ri)
{
    if (gretl_bundle_has_key(ri->b, "lambda_scale")) {
	ri->lamscale = gretl_bundle_get_int(ri->b, "lambda_scale", NULL);
    }
}

static void prepare_admm_params (regls_info *ri)
{
    gretl_matrix *ctrl;
    int len;

    /* set defaults */
    admm_reltol = ADMM_RELTOL_DEFAULT;
    admm_abstol = ADMM_ABSTOL_DEFAULT;

    ctrl = gretl_bundle_get_matrix(ri->b, "admmctrl", NULL);
    len = gretl_vector_get_length(ctrl);

    if (len > 0 && ctrl->val[0] > 0) {
	ri->rho = ctrl->val[0];
    }
    if (len > 1 && ctrl->val[1] > 0) {
	admm_reltol = ctrl->val[1];
    }
    if (len > 2 && ctrl->val[2] > 0) {
	admm_abstol = ctrl->val[2];
    }

    /* scale the absolute tolerance */
    admm_abstol *= sqrt(ri->X->cols);
}

static int get_xvalidation_details (regls_info *ri)
{
    int err = 0;

    ri->nf = gretl_bundle_get_int(ri->b, "nfolds", &err);
    ri->randfolds = gretl_bundle_get_bool(ri->b, "randfolds", 0);

    if (!err && (ri->nf < 2 || ri->nf > ri->n)) {
	gretl_errmsg_set("Invalid number of folds for cross validation");
	err = E_INVARG;
    }

    return err;
}

static int get_bundled_lfrac (regls_info *ri,
                              gretl_bundle *b)
{
    GretlType t = gretl_bundle_get_member_type(b, "lfrac", NULL);
    int err = 0;

    if (t == GRETL_TYPE_MATRIX) {
        ri->lfrac = gretl_bundle_get_matrix(b, "lfrac", &err);
    } else if (t == GRETL_TYPE_DOUBLE) {
        double x = gretl_bundle_get_scalar(b, "lfrac", &err);

        if (!err) {
            ri->lfrac = gretl_matrix_from_scalar(x);
            ri->free_lf = 1;
        }
    }

    return err;
}

static regls_info *regls_info_new (gretl_matrix *X,
				   gretl_matrix *y,
				   gretl_bundle *b,
				   PRN *prn, int *err)
{
    regls_info *ri = calloc(1, sizeof *ri);

    if (ri == NULL) {
	*err = E_ALLOC;
    } else {
        *err = get_bundled_lfrac(ri, b);
    }

    if (!*err) {
	ri->b = b;
	ri->X = X;
	ri->y = y;
	ri->stdize =  gretl_bundle_get_int(b, "stdize", err);
	ri->xvalid =  gretl_bundle_get_int(b, "xvalidate", err);
	ri->verbose = gretl_bundle_get_bool(b, "verbosity", 1);
	if (gretl_bundle_has_key(b, "alpha")) {
	    ri->alpha = gretl_bundle_get_scalar(b, "alpha", NULL);
	    if (ri->alpha == 0) {
		ri->ridge = 1;
	    }
	} else {
	    ri->ridge = gretl_bundle_get_bool(b, "ridge", 0);
	    ri->alpha = ri->ridge ? 0 : 1;
	}
	if (ri->alpha > 0 && ri->alpha < 1) {
	    ri->ccd = 1;
	} else {
	    ri->ccd = gretl_bundle_get_bool(b, "ccd", 0);
	}
    }

    if (*err) {
	free(ri);
	ri = NULL;
    } else {
	ri->prn = prn;
	ri->R2 = ri->crit = ri->BIC = ri->edf = NULL;
	ri->n = ri->X->rows;
	ri->k = ri->X->cols;
	ri->nlam = gretl_vector_get_length(ri->lfrac);
	ri->rho = 8.0;
	ri->infnorm = 0.0;
	ri->lamscale = LAMSCALE_GLMNET;
	ri->Xty = NULL;
	if (ri->ccd) {
	    prepare_ccd_param(ri);
	} else if (!ri->ridge && !ri->ccd) {
	    prepare_admm_params(ri);
	}
	if (ri->alpha < 1) {
	    maybe_set_lambda_scale(ri);
	    ri->edf = gretl_matrix_alloc(ri->nlam, 1);
	    if (ri->edf == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (!*err && ri->xvalid) {
	    *err = get_xvalidation_details(ri);
	} else if (!*err) {
	    ri->nf = ri->randfolds = ri->use_1se = 0;
	    ri->crit = gretl_zero_matrix_new(ri->nlam, 1);
	    ri->R2 = gretl_zero_matrix_new(ri->nlam, 1);
	    ri->BIC = gretl_zero_matrix_new(ri->nlam, 1);
	    if (ri->R2 == NULL || ri->crit == NULL || ri->BIC == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    return ri;
}

static void regls_info_free (regls_info *ri)
{
    if (ri != NULL) {
	gretl_matrix_free(ri->Xty);
	gretl_matrix_free(ri->R2);
	gretl_matrix_free(ri->crit);
	gretl_matrix_free(ri->BIC);
        if (ri->free_lf) {
            gretl_matrix_free(ri->lfrac);
        }
	free(ri);
    }
}

/* For when we're not doing cross validation: push various
   statistics into the output bundle
*/

static void regls_set_crit_data (regls_info *ri)
{
    if (ri->nlam > 1) {
	gretl_bundle_donate_data(ri->b, "crit", ri->crit, GRETL_TYPE_MATRIX, 0);
	if (ri->BIC != NULL) {
	    gretl_bundle_donate_data(ri->b, "BIC", ri->BIC, GRETL_TYPE_MATRIX, 0);
	}
	if (ri->R2 != NULL) {
	    gretl_bundle_donate_data(ri->b, "R2", ri->R2, GRETL_TYPE_MATRIX, 0);
	}
	if (ri->edf != NULL) {
	    gretl_bundle_donate_data(ri->b, "edf", ri->edf, GRETL_TYPE_MATRIX, 0);
	}
	ri->crit = ri->BIC = ri->R2 = ri->edf = NULL;
    } else {
	gretl_bundle_set_scalar(ri->b, "crit", ri->crit->val[0]);
	if (ri->BIC != NULL) {
	    gretl_bundle_set_scalar(ri->b, "BIC", ri->BIC->val[0]);
	}
	if (ri->R2 != NULL) {
	    gretl_bundle_set_scalar(ri->b, "R2", ri->R2->val[0]);
	}
	if (ri->edf != NULL) {
	    gretl_bundle_set_scalar(ri->b, "edf", ri->edf->val[0]);
	}
    }
}

static double vector_infnorm (const gretl_vector *z)
{
    const int n = gretl_vector_get_length(z);
    double azi, ret = 0;
    int i;

    for (i=0; i<n; i++) {
	azi = fabs(z->val[i]);
	if (azi > ret) ret = azi;
    }

    return ret;
}

/* compute X'y and its infinity-norm for all training data */

static int regls_set_Xty (regls_info *ri)
{
    int err = 0;

    if (ri->Xty == NULL) {
	ri->Xty = gretl_matrix_alloc(ri->X->cols, 1);
	if (ri->Xty == NULL) {
	    err = E_ALLOC;
	}
        ri->Xty = gretl_matrix_alloc(ri->X->cols, 1);
        if (ri->Xty == NULL) {
            err = E_ALLOC;
        }
    }
    if (!err) {
	gretl_matrix_multiply_mod(ri->X, GRETL_MOD_TRANSPOSE,
				  ri->y, GRETL_MOD_NONE,
				  ri->Xty, GRETL_MOD_NONE);
	ri->infnorm = vector_infnorm(ri->Xty);
	if (ri->ccd || ri->ridge) {
	    ri->infnorm /= ri->n;
	}
#if LAMBDA_DEBUG
	fprintf(stderr, "regls_set_Xty: infnorm = %g\n", ri->infnorm);
#endif
    }

    return err;
}

static int randomize_rows (gretl_matrix *X, gretl_matrix *y)
{
    gretl_vector *vp;
    double x, tmp;
    int i, j, src;

    vp = gretl_matrix_alloc(X->rows, 1);
    if (vp == NULL) {
	return E_ALLOC;
    }

    fill_permutation_vector(vp, X->rows);

    for (i=0; i<X->rows; i++) {
	src = vp->val[i] - 1;
	if (src == i) {
	    continue;
	}
	for (j=0; j<X->cols; j++) {
	    tmp = gretl_matrix_get(X, i, j);
	    x = gretl_matrix_get(X, src, j);
	    gretl_matrix_set(X, i, j, x);
	    gretl_matrix_set(X, src, j, tmp);
	}
	tmp = y->val[i];
	y->val[i] = y->val[src];
	y->val[src] = tmp;
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

static inline double max (double x, double y)
{
    return x >= y ? x : y;
}

#if defined(USE_SIMD)

static inline double hsum_double_avx (__m256d v)
{
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    __m128d high64;

    vlow   = _mm_add_pd(vlow, vhigh);
    high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}

static void vector_add_into (const double * __restrict a,
			     const double * __restrict b,
			     double *c, int n)
{
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, c256;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(a);
	b256 = _mm256_loadu_pd(b);
	c256 = _mm256_add_pd(a256, b256);
	_mm256_storeu_pd(c, c256);
	a += 4;
	b += 4;
	c += 4;
    }
    for (i=0; i<rem; i++) {
	c[i] = a[i] + b[i];
    }
}

static void vector_add_to (double *a,
			   const double * __restrict b,
			   int n)
{
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, sum;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(a);
	b256 = _mm256_loadu_pd(b);
	sum = _mm256_add_pd(a256, b256);
	_mm256_storeu_pd(a, sum);
	a += 4;
	b += 4;
    }
    for (i=0; i<rem; i++) {
	a[i] += b[i];
    }
}

/* a = a - b */

static void vector_subtract_from (double *a,
				  const double * __restrict b,
				  int n)
{
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, dif;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(a);
	b256 = _mm256_loadu_pd(b);
	dif = _mm256_sub_pd(a256, b256);
	_mm256_storeu_pd(a, dif);
	a += 4;
	b += 4;
    }
    for (i=0; i<rem; i++) {
	a[i] -= b[i];
    }
}

/* c = a - b */

static void vector_subtract_into (const double * __restrict a,
				  const double * __restrict b,
				  double *c, int n,
				  int cumulate)
{
    int imax = n / 4;
    int rem = n % 4;
    int i;

    __m256d a256, b256, c256;

    for (i=0; i<imax; i++) {
	a256 = _mm256_loadu_pd(a);
	b256 = _mm256_loadu_pd(b);
	if (cumulate) {
	    __m256d d256 = _mm256_sub_pd(a256, b256);

	    c256 = _mm256_loadu_pd(c);
	    d256 = _mm256_add_pd(c256, d256);
	    _mm256_storeu_pd(c, d256);
	} else {
	    c256 = _mm256_sub_pd(a256, b256);
	    _mm256_storeu_pd(c, c256);
	}
	a += 4;
	b += 4;
	c += 4;
    }
    for (i=0; i<rem; i++) {
	if (cumulate) {
	    c[i] += a[i] - b[i];
	} else {
	    c[i] = a[i] - b[i];
	}
    }
}

/* compute q = rho * (b - u) + X'y */

static inline void compute_q (gretl_vector *q,
			      const gretl_vector *b,
			      const gretl_vector *u,
			      const gretl_vector *a, /* X'y */
			      double rho, int n)
{
    __m256d b256, u256, a256;
    __m256d r256, tmp;
    const double *bx = b->val;
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

    for (i=0; i<imax; i++) {
	b256 = _mm256_loadu_pd(bx);
	u256 = _mm256_loadu_pd(ux);
	a256 = _mm256_loadu_pd(ax);
	/* subtract u from b */
	tmp = _mm256_sub_pd(b256, u256);
	if (mul) {
	    /* multiply by rho */
	    tmp = _mm256_mul_pd(tmp, r256);
	}
	/* add a */
	tmp = _mm256_add_pd(tmp, a256);
	/* write result into q */
	_mm256_storeu_pd(qx, tmp);
	bx += 4;
	ux += 4;
	ax += 4;
	qx += 4;
    }

    for (i=0; i<rem; i++) {
	if (mul) {
	    qx[i] = rho * (bx[i] - ux[i]) + ax[i];
	} else {
	    qx[i] = bx[i] - ux[i] + ax[i];
	}
    }
}

static double dot_product (const double * __restrict x,
                           const double * __restrict y, int n)
{
    double ret = 0.0;
    int i, imax = n / 4;
    int rem = n % 4;

    __m256d x256, y256, tmp;

    for (i=0; i<imax; i++) {
	x256 = _mm256_loadu_pd(x);
	y256 = _mm256_loadu_pd(y);
	tmp = _mm256_mul_pd(x256, y256);
	ret += hsum_double_avx(tmp);
	x += 4;
	y += 4;
    }

    for (i=0; i<rem; i++) {
	ret += x[i] * y[i];
    }

    return ret;
}

#else

static void vector_add_into (const double * __restrict a,
			     const double * __restrict b,
			     double *c, int n)
{
    int i;

    for (i=0; i<n; i++) {
	c[i] = a[i] + b[i];
    }
}

static void vector_add_to (double *a,
                           const double * __restrict b,
			   int n)
{
    int i;

    for (i=0; i<n; i++) {
	a[i] += b[i];
    }
}

static void vector_subtract_from (double *a,
				  const double * __restrict b,
				  int n)
{
    int i;

    for (i=0; i<n; i++) {
	a[i] -= b[i];
    }
}

static void vector_subtract_into (const double * __restrict a,
				  const double * __restrict b,
				  double *c, int n,
				  int cumulate)
{
    int i;

    for (i=0; i<n; i++) {
	if (cumulate) {
	    c[i] += a[i] - b[i];
	} else {
	    c[i] = a[i] - b[i];
	}
    }
}

static double dot_product (const double * __restrict x,
                           const double * __restrict y,
                           int n)
{
    double ret = 0.0;
    int i;

    for (i=0; i<n; i++) {
	ret += x[i] * y[i];
    }
    return ret;
}

static inline void compute_q (gretl_vector *q,
			      const gretl_vector *b,
			      const gretl_vector *u,
			      const gretl_vector *Xty,
			      double rho, int n)
{
    const int mul = rho != 1.0;
    int i;

    for (i=0; i<n; i++) {
	if (mul) {
	    q->val[i] = rho * (b->val[i] - u->val[i]) + Xty->val[i];
	} else {
	    q->val[i] = b->val[i] - u->val[i] + Xty->val[i];
	}
    }
}

#endif /* AVX or not */

static double own_dot_product (const gretl_vector *x)
{
    int i, n = gretl_vector_get_length(x);
    double ret = 0.0;

    for (i=0; i<n; i++) {
        ret += x->val[i] * x->val[i];
    }

    return ret;
}

/* fortran: dot_product(X(:,j), X(:,k)) for @X with @n rows */

static double dot_prod_jk (const gretl_matrix *X, int j, int k, int n)
{
    const double *xj = X->val + n * j;
    const double *xk = X->val + n * k;

    return dot_product(xj, xk, n);
}

/* fortran: dot_product(v(1:n), m(j,1:n)) */

static double dot_prod_vm (const double *v,
			   const gretl_matrix *m,
			   int j, int n)
{
    double ret = 0;
    int i;

    for (i=0; i<n; i++) {
	ret += v[i] * gretl_matrix_get(m, j, i);
    }

    return ret;
}

/* implement these fortran lines:

   x(1:n) = y(idx(1:n)) - x(1:n) !! sub = 1
   x(1:n) = y(idx(1:n))          !! sub = 0
*/

static void range_set_sub (double *x, const double *y,
			   const int *idx, int n, int sub)
{
    int i;

    if (sub) {
	for (i=0; i<n; i++) {
	    x[i] = y[idx[i]] - x[i];
	}
    } else {
	for (i=0; i<n; i++) {
	    x[i] = y[idx[i]];
	}
    }
}

/* fortran: B(1:n,j) = a(idx(1:n)) */

static void fill_coeff_column (gretl_matrix *B, int nx, int j,
			       const double *a, const int *idx,
			       int n)
{
    int i, offset = B->rows > nx;
    double *b = B->val + j * B->rows;

    for (i=0; i<n; i++) {
	b[i+offset] = a[idx[i]];
    }
}

/* sign(x,y): gives "the value of x with the sign of y",
   but in context @x will always be positive.
*/

static inline double sign (double x, double y)
{
    return y >= 0 ? x : -x;
}

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

/* Cyclical Coordinate Descent (CCD) auxiliary functions */

static int ccd_scale (gretl_matrix *x, double *y,
		      double *xty, double *xv)
{
    int i, j, n = x->rows;
    double *xj, v = sqrt(1.0/n);

    for (i=0; i<n; i++) {
	y[i] *= v;
    }
    for (j=0; j<x->cols; j++) {
	xj = x->val + j * n;
	for (i=0; i<n; i++) {
	    xj[i] *= v;
	}
	if (xv != NULL) {
	    xv[j] = dot_product(xj, xj, n);
	}
	if (xty != NULL) {
	    xty[j] = dot_product(y, xj, n);
	}
    }

    return 0;
}

static void finalize_ccd_coeffs (gretl_matrix *B,
				 double *a, int nx,
				 int *ia)
{
    int offset = B->rows > nx;
    size_t asize = nx * sizeof *a;
    double *bj;
    int i, j;

    for (j=0; j<B->cols; j++) {
	bj = B->val + j*B->rows + offset;
	memcpy(a, bj, asize);
	for (i=0; i<nx; i++) {
	    bj[i] = 0.0;
	}
	for (i=0; i<nx; i++) {
	    if (a[i] != 0) {
		bj[ia[i]] = a[i];
	    }
	}
    }
}

static int ccd_iteration (double alpha, const gretl_matrix *X, double *g,
			  int nlam, const double *ulam, double thr,
			  int maxit, const double *xv, int *lmu,
			  gretl_matrix *B, int *ia, int *kin,
			  double *Rsq, int *pnlp)
{
    gretl_matrix *C;
    double alm, u, v, rsq = 0;
    double ak, del, dlx, cij;
    double omb, dem, ab;
    double *a, *da;
    int *mm, nin, jz, iz = 0;
    int j, k, l, m, nlp = 0;
    int nx = X->cols;
    int bad_R2 = 0;
    int err = 0;

#if 0
    fputs("ccd_iteration args:\n", stderr);
    fprintf(stderr, "alpha = %g\n", alpha);
    gretl_matrix_print(X, "X");
    fprintf(stderr, "g = %g\n", g[0]);
    fprintf(stderr, "nlam = %d\n", nlam);
    fprintf(stderr, "ulam = %g\n", ulam[0]);
    fprintf(stderr, "thr = %g\n", thr);
    fprintf(stderr, "maxit = %d\n", maxit);
    fprintf(stderr, "xv = %g\n", xv[0]);
    fprintf(stderr, "lmu = %d\n", lmu[0]);
    gretl_matrix_print(B, "b");
    fprintf(stderr, "ia = %d\n", ia[0]);
    fprintf(stderr, "kin = %d\n", kin[0]);
    fprintf(stderr, "nlp = %d\n", *pnlp);
#endif

    C = gretl_matrix_alloc(nx, nx);
    a = malloc(nx * sizeof *a);
    da = malloc(nx * sizeof *da);
    mm = malloc(nx * sizeof *mm);
    if (C == NULL || a == NULL || da == NULL || mm == NULL) {
	fprintf(stderr, "ccd: allocation failure (nx = %d)\n", nx);
	return E_ALLOC;
    }
    /* "zero" @a and @mm */
    for (j=0; j<nx; j++) {
	a[j] = 0.0;
	mm[j] = -1;
    }
    nin = nlp = *pnlp = 0;
    omb = 1.0 - alpha; /* = 0 for lasso */

    for (m=0; m<nlam; m++) {
	alm = ulam[m];
	dem = alm*omb;
	ab = alm*alpha;
	jz = 1;
    maybe_restart:
	if (iz * jz == 0) {
            nlp++;
            dlx = 0.0;
	    for (k=0; k<nx; k++) {
		ak = a[k];
		u = g[k] + ak*xv[k];
		v = fabs(u) - ab;
		a[k] = v > 0.0 ? sign(v,u) / (xv[k]+dem) : 0.0;
		if (a[k] != ak) {
		    if (mm[k] < 0) {
			if (nin >= nx) goto check_conv;
			for (j=0; j<nx; j++) {
			    if (mm[j] >= 0) {
				cij = gretl_matrix_get(C, k, mm[j]);
				gretl_matrix_set(C, j, nin, cij);
			    } else if (j != k) {
				cij = dot_prod_jk(X, j, k, X->rows);
				gretl_matrix_set(C, j, nin, cij);
			    } else {
				gretl_matrix_set(C, j, nin, xv[j]);
			    }
			}
			mm[k] = nin;
			ia[nin] = k;
			nin++;
		    }
		    del = a[k] - ak;
		    rsq += del * (2*g[k] - del*xv[k]);
		    dlx = max(xv[k]*del*del, dlx);
		    for (j=0; j<nx; j++) {
			cij = gretl_matrix_get(C, j, mm[k]);
			g[j] -= cij*del;
		    }
		}
            }
	check_conv:
	    if (dlx < thr || nin > nx) {
		goto m_finish;
	    } else if (nlp > maxit) {
		fputs("ccd: max iters reached\n", stderr);
		err = E_NOCONV;
		goto getout;
            }
	}
	iz = 1;
	range_set_sub(da, a, ia, nin, 0);
    nlp_plus:
	nlp++;
	dlx = 0.0;
	for (l=0; l<nin; l++) {
	    k = ia[l];
	    ak = a[k];
	    u = g[k] + ak*xv[k];
	    v = fabs(u) - ab;
	    a[k] = v > 0.0 ? sign(v,u) / (xv[k]+dem) : 0.0;
	    if (a[k] != ak) {
		del = a[k] - ak;
		rsq += del * (2*g[k] - del*xv[k]);
		dlx = max(xv[k]*del*del, dlx);
		for (j=0; j<nin; j++) {
		    cij = gretl_matrix_get(C, ia[j], mm[k]);
		    g[ia[j]] -= cij*del;
		}
	    }
	}
	if (dlx < thr) {
	    range_set_sub(da, a, ia, nin, 1);
	    for (j=0; j<nx; j++) {
		if (mm[j] < 0) {
		    g[j] -= dot_prod_vm(da, C, j, nin);
		}
	    }
	    jz = 0;
	    goto maybe_restart;
	} else if (nlp <= maxit) {
	    /* try another iteration */
	    goto nlp_plus;
	} else {
	    /* reached max iterations */
	    err = E_NOCONV;
	    fprintf(stderr, "ccd: nlp = %d, maxit %d\n", nlp, maxit);
	    goto getout;
	}
    m_finish:
	if (nin <= nx) {
	    if (nin > 0) {
		fill_coeff_column(B, nx, m, a, ia, nin);
	    }
	    kin[m] = nin;
	    if (Rsq != NULL) {
		if (rsq > 1) {
		    bad_R2 = 1;
		}
		Rsq[m] = rsq;
	    }
	    *lmu = m + 1;
	} else {
            err = E_NOCONV;
	    fprintf(stderr, "ccd: error at foot of loop\n");
            goto getout;
	}
    } /* end loop over lambda values */

 getout:

    if (!err) {
	finalize_ccd_coeffs(B, a, nx, ia);
	if (bad_R2) {
	    Rsq[0] = NADBL;
	}
    }

    *pnlp = nlp;
    free(a);
    free(mm);
    free(da);
    gretl_matrix_free(C);

    return err;
}

static int ccd_get_crit (const gretl_matrix *B,
		         const gretl_matrix *lam,
			 regls_info *ri)
{
    double *bj, l1, l2, SSR;
    double ll, llc, edf = 0;
    double lambda, nulldev = 1.0;
    double alpha = ri->alpha;
    double penalty;
    int imin = B->rows > ri->k;
    int dfj, n = ri->n;
    int i, j;

    if (!ri->stdize) {
	/* in case @y is in fact non-standard */
	const double *y = ri->y->val;

	nulldev = 0.0;
	for (i=0; i<n; i++) {
	    nulldev += y[i] * y[i];
	}
    }

    llc = -0.5 * n * (1 + LN_2_PI - log(n));

    for (j=0; j<B->cols; j++) {
	lambda = lam->val[j];
	l1 = l2 = 0;
	dfj = 0;
	bj = B->val + j*B->rows;
	for (i=imin; i<B->rows; i++) {
	    if (alpha == 1) {
		/* lasso */
		l1 += fabs(bj[i]);
		dfj += bj[i] != 0;
	    } else if (alpha == 0) {
		/* ridge */
		l2 += bj[i] * bj[i];
	    } else {
		l1 += alpha * fabs(bj[i]);
		l2 += bj[i] * bj[i];
		dfj += bj[i] != 0;
	    }
	}
	SSR = nulldev * (1.0 - ri->R2->val[j]);
	/* with CCD, y and X are scaled by 1/sqrt(n) */
	ll = llc - 0.5 * n * log(n*SSR);
	if (alpha == 1) {
	    /* lasso */
	    gretl_vector_set(ri->crit, j, 0.5 * SSR + lambda * l1);
	    gretl_vector_set(ri->BIC, j, -2 * ll + dfj * log(n));
	} else if (alpha == 0) {
	    /* ridge */
	    edf = ri->edf->val[j];
	    gretl_vector_set(ri->crit, j, SSR + lambda * l2);
	    gretl_vector_set(ri->BIC, j, -2 * ll + edf * log(n));
	} else {
	    /* elnet */
	    edf = ri->edf->val[j];
	    penalty = 0.5 * (1 - alpha) * l2 + alpha * l1;
	    gretl_vector_set(ri->crit, j, 0.5 * SSR + lambda * penalty);
	    gretl_vector_set(ri->BIC, j, -2 * ll + edf * log(n));
	}
    }

    return 0;
}

/* We call ridge_effective_df() only if we're doing ridge
   via CCD -- otherwise the effective df gets computed as
   part of the larger SVD calculation.
*/

static int ridge_effective_df (const gretl_matrix *lam,
			       regls_info *ri)
{
    gretl_matrix *s = NULL;
    int err;

    err = gretl_matrix_SVD(ri->X, NULL, &s, NULL, 0);

    if (!err) {
	int i, j, k = gretl_vector_get_length(s);
	double sv2, edfj;

	for (i=0; i<k; i++) {
	    sv2 = s->val[i] * s->val[i];
	    s->val[i] = sv2;
	}
	for (j=0; j<ri->nlam; j++) {
	    edfj = 0;
	    for (i=0; i<k; i++) {
		edfj += s->val[i] / (s->val[i] + lam->val[j]);
	    }
	    ri->edf->val[j] = edfj;
	}
	gretl_matrix_free(s);
    }

    return err;
}

static int elnet_effective_df (const gretl_matrix *lam,
			       const gretl_matrix *B,
			       regls_info *ri)
{
    gretl_matrix *XTX = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *xi = NULL;
    double *dest, *src;
    double xii, dfj;
    size_t csize;
    int i, j, t;
    int err = 0;

    X = gretl_matrix_copy(ri->X);
    XTX = gretl_matrix_alloc(ri->k, ri->k);
    xi = gretl_matrix_alloc(1, ri->k);

    if (X == NULL || XTX == NULL || xi == NULL) {
	return E_ALLOC;
    }

    csize = ri->n * sizeof(double);

    for (j=0; j<ri->nlam; j++) {
	double lam2 = lam->val[j] * (1 - ri->alpha)/2;
	int inv_err, kj = 0;

	dest = X->val;
	dfj = 0;

	for (i=0; i<ri->k; i++) {
	    if (gretl_matrix_get(B, i, j) != 0) {
		src = ri->X->val + i * ri->n;
		memcpy(dest, src, csize);
		dest += ri->n;
		kj++;
	    }
	}
	if (kj == 0) {
	    ri->edf->val[j] = 0;
	    continue;
	}
	gretl_matrix_reuse(X, -1, kj);
	gretl_matrix_reuse(XTX, kj, kj);
	gretl_matrix_reuse(xi, 1, kj);

	gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
				  X, GRETL_MOD_NONE,
				  XTX, GRETL_MOD_NONE);
	for (i=0; i<kj; i++) {
	    xii = gretl_matrix_get(XTX, i, i);
	    gretl_matrix_set(XTX, i, i, xii + lam2);
	}
	inv_err = gretl_invert_symmetric_matrix(XTX);
	if (inv_err) {
	    fprintf(stderr, "elnet df: inversion failed for j=%d\n", j);
	    ri->edf->val[j] = NADBL;
	} else {
	    for (t=0; t<ri->n; t++) {
		for (i=0; i<kj; i++) {
		    xi->val[i] = gretl_matrix_get(X, t, i);
		}
		dfj += gretl_scalar_qform(xi, XTX, &err);
	    }
	    ri->edf->val[j] = dfj;
	}
    }

    gretl_matrix_free(X);
    gretl_matrix_free(XTX);
    gretl_matrix_free(xi);

    return err;
}

static gchar *crit_print_format (const gretl_matrix *crit,
				 int ridge)
{
    gchar *fmt = NULL;

    if (ridge) {
	fmt = g_strdup_printf("%%12f  %%6.2f   %%.4f   %%#g\n");
    } else {
	fmt = g_strdup_printf("%%12f  %%5d    %%f   %%.4f  %%#g\n");
    }

    return fmt;
}

static void lambda_sequence_header (PRN *prn)
{
    pputc(prn, '\n');
    pputs(prn, "    lambda/n     df   criterion      R^2      BIC\n");
}

static void ccd_print (const gretl_matrix *B,
		       const gretl_matrix *lam,
		       regls_info *ri)
{
    gchar *cfmt = NULL;
    double *bj;
    int k = B->rows;
    int nlam = B->cols;
    int i, j, dfj;

    if (ri->crit != NULL) {
	/* header for output showing penalized criterion */
	lambda_sequence_header(ri->prn);
    } else {
	/* as per R, more or less */
	pputc(ri->prn, '\n');
	pputs(ri->prn, "    df     R^2  lambda    BIC\n");
    }

    cfmt = crit_print_format(ri->crit, 0);

    for (j=0; j<nlam; j++) {
	bj = B->val + j*k;
	dfj = 0;
	for (i=0; i<k; i++) {
	    dfj += fabs(bj[i]) > 0;
	}
	if (ri->crit != NULL) {
	    pprintf(ri->prn, cfmt, lam->val[j], dfj, ri->crit->val[j],
		    ri->R2->val[j], ri->BIC->val[j]);
	} else {
	    pprintf(ri->prn, "%-2d  %2d  %.4f  %.4f  %#g\n", j+1, dfj,
		    ri->R2->val[j], lam->val[j], ri->BIC->val[j]);
	}
    }

    g_free(cfmt);
}

/* This also serves for printing elastic net results */

static void ridge_print (const gretl_matrix *lam,
			 regls_info *ri)
{
    gchar *cfmt = NULL;
    int j;

    pprintf(ri->prn, "\n  %s\n\n", _("df = effective number of free parameters"));
    pputs(ri->prn, "      lambda      df      R^2       BIC\n");

    cfmt = crit_print_format(ri->crit, 1);

    for (j=0; j<ri->nlam; j++) {
	pprintf(ri->prn, cfmt, lam->val[j], ri->edf->val[j],
		ri->R2->val[j], ri->BIC->val[j]);
    }
    g_free(cfmt);
}

static void xv_ridge_print (const gretl_matrix *lam,
			    regls_info *ri)
{
    int j;

    pputc(ri->prn, '\n');
    pputs(ri->prn, "      lambda     df\n");
    for (j=0; j<ri->nlam; j++) {
	pprintf(ri->prn, "%12f  %.3f\n", lam->val[j], ri->edf->val[j]);
    }
}

/* end functions specific to CCD */

/* calculate the lasso objective function */

static double lasso_objective (const gretl_matrix *X,
			       const gretl_vector *y,
			       const gretl_vector *b,
			       double lambda,
			       gretl_vector *u,
			       double *pSSR,
			       double *pR2)
{
    double TSS, SSR, obj;

    TSS = own_dot_product(y);
    gretl_matrix_multiply(X, b, u);
    vector_subtract_from(u->val, y->val, y->rows);
    SSR = own_dot_product(u);
    obj = 0.5 * SSR + lambda * abs_sum(b);
    *pR2 = 1.0 - SSR/TSS;
    if (pSSR != NULL) {
	*pSSR = SSR;
    }

    return obj / y->rows;
}

/* calculate the cross validation criterion */

static double xv_score (const gretl_matrix *X,
			const gretl_vector *y,
			const gretl_vector *b,
			gretl_vector *Xb)
{
    double sum = 0;

    /* get fitted values */
    gretl_matrix_multiply(X, b, Xb);
    /* compute and process residuals */
    vector_subtract_from(Xb->val, y->val, X->rows);
    sum = own_dot_product(Xb);

    return sum / X->rows;
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

static int get_cholesky_factor (const gretl_matrix *X,
				gretl_matrix *L,
				double rho)
{
    double d;
    int i;

    if (X->rows >= X->cols) {
	/* "skinny": L = chol(X'X + rho*I) */
	gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
				  X, GRETL_MOD_NONE,
				  L, GRETL_MOD_NONE);
	for (i=0; i<X->cols; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + rho);
	}
    } else {
	/* "fat": L = chol(I + 1/rho*XX') */
	gretl_matrix_multiply_mod(X, GRETL_MOD_NONE,
				  X, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_NONE);
	if (rho != 1.0) {
	    gretl_matrix_multiply_by_scalar(L, 1/rho);
	}
	for (i=0; i<X->rows; i++) {
	    d = gretl_matrix_get(L, i, i);
	    gretl_matrix_set(L, i, i, d + 1.0);
	}
    }

    return gretl_matrix_cholesky_decomp(L);
}

static int admm_iteration (const gretl_matrix *X,
			   const gretl_vector *Xty,
			   gretl_matrix *L,
			   gretl_vector *v, gretl_vector *b,
			   gretl_vector *u, gretl_vector *q,
			   gretl_vector *p, gretl_vector *r,
			   gretl_vector *bprev, gretl_vector *bdiff,
			   double lambda, double *prho,
			   int tune_rho, int *iters)
{
    double nxstack, nystack;
    double prires, dualres;
    double eps_pri, eps_dual;
    double rho = *prho;
    double nrm2, rho2 = rho*rho;
    int itermin = 1;
    int n = X->cols;
    int iter = 0;
    int err = 0;

#if RHO_DEBUG
    fprintf(stderr, "*** admm: lambda %g, rho %g ***\n", lambda, rho);
#endif

    while (iter < ADMM_MAX_ITER && !err) {
	/* u-update: u = u + r */
	vector_add_to(u->val, r->val, n);

	/* v-update: v = (X^T X + rho I) \ (X^T y + rho b - u) */

	compute_q(q, b, u, Xty, rho, n);
	if (X->rows >= X->cols) {
	    /* v = U \ (L \ q) */
	    gretl_cholesky_solve(L, q);
	    vector_copy_values(v, q, n);
	} else {
	    /* v = q/rho - 1/rho^2 * X^T * (U \ (L \ (X*q))) */
	    gretl_matrix_multiply(X, q, p);
	    err = gretl_cholesky_solve(L, p);
	    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
				      p, GRETL_MOD_NONE,
				      v, GRETL_MOD_NONE);
	    gretl_matrix_multiply_by_scalar(v, -1/rho2);
	    gretl_matrix_multiply_by_scalar(q, 1/rho);
	    vector_add_to(v->val, q->val, n);
	}

	/* sqrt(sum ||r_i||_2^2) */
	prires  = sqrt(own_dot_product(r));
	/* sqrt(sum ||v_i||_2^2) */
	nxstack = sqrt(own_dot_product(v));
	/* sqrt(sum ||u_i||_2^2) */
	nystack = own_dot_product(u) / rho2;
	nystack = sqrt(nystack);

	vector_copy_values(bprev, b, n);
	vector_add_into(v->val, u->val, b->val, n);
	soft_threshold(b, lambda, rho);

	/* Termination checks */

	/* dual residual */
	vector_subtract_into(b->val, bprev->val, bdiff->val, n, 0); /* bdiff = b - bprev */
	/* ||s^k||_2^2 = N rho^2 ||b - bprev||_2^2 */
	nrm2 = sqrt(own_dot_product(bdiff));
	dualres = rho * nrm2;

	/* compute primal and dual feasibility tolerances */
	nrm2 = sqrt(own_dot_product(b));
	eps_pri  = admm_abstol + admm_reltol * fmax(nxstack, nrm2);
	eps_dual = admm_abstol + admm_reltol * nystack;

	if (iter >= itermin && prires <= eps_pri && dualres <= eps_dual) {
	    /* converged */
	    break;
	}

	/* Compute residual: r = v - b */
	vector_subtract_into(v->val, b->val, r->val, n, 0);

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
		get_cholesky_factor(X, L, rho);
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

static gretl_matrix *make_coeff_matrix (regls_info *ri)
{
    gretl_matrix *B;
    int rows = ri->k + ri->stdize;

    B = gretl_zero_matrix_new(rows, ri->nlam);

    if (B != NULL) {
	gretl_bundle_donate_data(ri->b, "B", B, GRETL_TYPE_MATRIX, 0);
    }

    return B;
}

static void ccd_make_lambda (regls_info *ri,
			     gretl_matrix *lam,
			     double *lmax)
{
    int i;

#if LAMBDA_DEBUG
    fprintf(stderr, "ccd_make_lambda: lmax = %g\n", *lmax);
    fprintf(stderr, "ri->lamscale = %d, ri->n = %d\n", ri->lamscale, ri->n);
#endif

    gretl_matrix_copy_values(lam, ri->lfrac);
    if (ri->lamscale == LAMSCALE_NONE) {
	for (i=0; i<ri->nlam; i++) {
	    lam->val[i] /= ri->n;
	}
	return;
    }
    if (ri->alpha < 1.0) {
	*lmax /= max(ri->alpha, 1.0e-3);
#if LAMBDA_DEBUG
	fprintf(stderr, "revised lmax = %g\n", *lmax);
#endif
    }
    for (i=0; i<ri->nlam; i++) {
	lam->val[i] *= *lmax;
    }
    if (ri->alpha < 1.0 && ri->nlam > 1) {
	lam->val[0] = BIG_LAMBDA;
    }
}

static void lasso_lambda_report (regls_info *ri)
{
    pprintf(ri->prn, "lambda-max = %g\n", ri->infnorm);
}

/* Remedial R^2 calculation for CCD: it seems that we end up
   coming here only when standardization is turned off but
   the dependent variable is substantially non-standard,
   in which case ccd_iteration() can produce R^2 > 1.
*/

static int ccd_alt_R2 (regls_info *ri, gretl_matrix *B)
{
    gretl_matrix *bj, *yh;
    int n = ri->y->rows;
    int k = ri->X->cols;
    int err = 0;

    bj = gretl_matrix_alloc(k, 1);
    yh = gretl_matrix_alloc(n, 1);

    if (bj == NULL || yh == NULL) {
	err = E_ALLOC;
    } else {
	const double *y = ri->y->val;
	size_t sz = k * sizeof(double);
	double ui, SSR, TSS = 0;
	int i, j;

	for (i=0; i<n; i++) {
	    TSS += y[i] * y[i];
	}
	for (j=0; j<ri->nlam; j++) {
	    memcpy(bj->val, B->val + j*B->rows, sz);
	    gretl_matrix_multiply(ri->X, bj, yh);
	    SSR = 0;
	    for (i=0; i<n; i++) {
		ui = y[i] - yh->val[i];
		SSR += ui * ui;
	    }
	    ri->R2->val[j] = 1.0 - SSR/TSS;
	}
    }

    gretl_matrix_free(yh);
    gretl_matrix_free(bj);

    return err;
}

static int ccd_prep (regls_info *ri, ccd_info *ci)
{
    int nlam = ri->nlam;
    int k = ri->k;

#if 0
    fprintf(stderr, "*** ccd_prep ***\n");
#endif

    ci->MB = gretl_matrix_block_new(&ci->xv, k, 1,
				    &ci->Xty, k, 1,
				    &ci->lam, nlam, 1,
				    NULL);
    ci->B = gretl_zero_matrix_new(k + ri->stdize, nlam);

    if (ci->MB == NULL || ci->B == NULL) {
	return E_ALLOC;
    }

    /* scale data by sqrt(1/n) */
    ccd_scale(ri->X, ri->y->val, ci->Xty->val, ci->xv->val);

    /* and compute lambda sequence */
    ci->lmax = vector_infnorm(ci->Xty);
    ccd_make_lambda(ri, ci->lam, &ci->lmax);

    return 0;
}

/* Cyclical Coordinate Descent driver: we come here either
   to get coefficient estimates right away, or after
   cross validation. Handles both LASSO and Ridge.
*/

static int ccd_regls (regls_info *ri)
{
    ccd_info ci = {0};
    double *Rsq = NULL;
    int maxit = CCD_MAX_ITER;
    int *ia, *nnz;
    int nlp = 0, lmu = 0;
    int nlam = ri->nlam;
    int k = ri->k;
    int err;

    err = ccd_prep(ri, &ci);
    if (err) {
	return err;
    }

    /* integer workspace */
    ia = calloc(k + nlam, sizeof *ia);
    if (ia == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    nnz = ia + k;

    if (!ri->xvalid) {
	Rsq = ri->R2->val;
    }

    if (ri->edf != NULL && ri->alpha == 0) {
	/* we'll want this but it's not calculated by CCD */
	err = ridge_effective_df(ci.lam, ri);
    }

    if (ri->alpha == 1 && !ri->xvalid && ri->verbose) {
	lasso_lambda_report(ri);
    }

#if LAMBDA_DEBUG
    fprintf(stderr, "ccd_regls: ci.lmax = %g\n", ci.lmax);
    gretl_matrix_print(ci.lam, "lam in ccd_regls");
#endif
    err = ccd_iteration(ri->alpha, ri->X, ci.Xty->val, nlam, ci.lam->val,
			ccd_toler, maxit, ci.xv->val, &lmu, ci.B, ia,
			nnz, Rsq, &nlp);
    if (err) {
	goto bailout;
    }

    if (ri->edf != NULL && ri->alpha > 0 && ri->alpha < 1) {
	/* elastic net */
	elnet_effective_df(ci.lam, ci.B, ri);
    }

    if (Rsq != NULL && na(Rsq[0])) {
	/* remedy spurious R^2 > 1 */
	err = ccd_alt_R2(ri, ci.B);
	if (err) {
	    goto bailout;
	}
    }

    if (ri->lamscale == LAMSCALE_NONE) {
	gretl_matrix_multiply_by_scalar(ri->y, sqrt(ri->n));
	gretl_matrix_copy_values(ci.lam, ri->lfrac);
    } else if (ri->alpha < 1.0) {
	/* not entirely truthful! */
	ci.lam->val[0] = ri->lfrac->val[0] * ci.lmax;
    }

    if (ri->xvalid && ri->verbose > 1 && ri->ridge && nlam > 1) {
	xv_ridge_print(ci.lam, ri);
    }

    if (!ri->xvalid) {
	ccd_get_crit(ci.B, ci.lam, ri);
	if (ri->verbose) {
	    if (ri->alpha < 1) {
		ridge_print(ci.lam, ri);
	    } else {
		ccd_print(ci.B, ci.lam, ri);
	    }
	}
	if (nlam > 1) {
	    double BICmin = 1e200;
	    int j, idxmin = 0;

	    for (j=0; j<nlam; j++) {
		if (ri->BIC->val[j] < BICmin) {
		    BICmin = ri->BIC->val[j];
		    idxmin = j;
		}
	    }
	    gretl_bundle_set_scalar(ri->b, "idxmin", idxmin + 1);
	    gretl_bundle_set_scalar(ri->b, "lfmin", ri->lfrac->val[idxmin]);
	}
	regls_set_crit_data(ri);
    }

    if (!err) {
	gretl_bundle_donate_data(ri->b, "B", ci.B, GRETL_TYPE_MATRIX, 0);
	ci.B = NULL;
	if (ri->lamscale != LAMSCALE_NONE) {
	    gretl_bundle_set_scalar(ri->b, "lmax", ci.lmax * ri->n);
	}
	if (nlam == 1) {
	    double lambda = ri->lfrac->val[0];

	    if (ri->lamscale != LAMSCALE_NONE) {
		/* show a value comparable with ADMM (??) */
		lambda *= ci.lmax * ri->n;
	    }
	    gretl_bundle_set_scalar(ri->b, "lambda", lambda);
	}
    }

 bailout:

    gretl_matrix_free(ci.B);
    gretl_matrix_block_destroy(ci.MB);
    free(ia);

    return err;
}

/* Variant of SVD ridge that computes the covariance
   matrix as well as the parameter vector */

static int svd_ridge_vcv (regls_info *ri,
			  double lam,
			  gretl_matrix *B,
			  gretl_matrix **pV)
{
    gretl_matrix_block *MB = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *Vt = NULL;
    gretl_matrix *sv = NULL;
    gretl_matrix *sve = NULL;
    gretl_matrix *RI = NULL;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *Ve = NULL;
    gretl_matrix *u = NULL;
    gretl_matrix *b = NULL;
    double vij, SSR, s2;
    int n = ri->X->rows;
    int k = ri->X->cols;
    int offset = 0;
    int i, j, err = 0;

#if RIDGE_DEBUG
    fprintf(stderr, "*** svd_ridge_vcv, lam = %g ***\n", lam);
#endif

    err = gretl_matrix_SVD(ri->X, NULL, &sv, &Vt, 0);

    if (!err) {
	MB = gretl_matrix_block_new(&sve, 1, k,
				    &u, n, 1,
				    &RI, k, k,
				    &Ve, k, k,
				    &Tmp, k, k,
				    &b, k, 1, NULL);
	if (MB == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (!err) {
	V = gretl_matrix_alloc(k, k);
	if (V == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (ri->edf != NULL) {
	ri->edf->val[0] = 0.0;
    }

    /* sve = 1 / (sv.^2 + lambda) */
    for (i=0; i<k; i++) {
	sve->val[i] = 1.0 / (sv->val[i] * sv->val[i] + lam);
	if (ri->edf != NULL) {
	    ri->edf->val[0] += sv->val[i] * sv->val[i] * sve->val[i];
	}
    }

    /* Ve = Vt' .* sve */
    for (j=0; j<k; j++) {
	for (i=0; i<k; i++) {
	    vij = gretl_matrix_get(Vt, j, i);
	    gretl_matrix_set(Ve, i, j, vij * sve->val[j]);
	}
    }

    /* RI = Ve * Vt */
    gretl_matrix_multiply(Ve, Vt, RI);

    /* b = RI * Xty */
    gretl_matrix_multiply(RI, ri->Xty, b);

    /* transcribe b to @B */
    offset = B->rows > k ? 1 : 0;
    memcpy(B->val + offset, b->val, k * sizeof *b->val);

    /* u = X*b - y */
    gretl_matrix_multiply(ri->X, b, u);
    gretl_matrix_subtract_from(u, ri->y);

    /* residual variance */
    SSR = own_dot_product(u);
    s2 = SSR / n;

    /* V = s2 * ridgeI * X'X * ridgeI */
    gretl_matrix_multiply_mod(ri->X, GRETL_MOD_TRANSPOSE,
			      ri->X, GRETL_MOD_NONE,
			      Ve, GRETL_MOD_NONE);
    gretl_matrix_multiply(RI, Ve, Tmp);
    gretl_matrix_multiply(Tmp, RI, V);
    gretl_matrix_multiply_by_scalar(V, s2);

    if (ri->R2 != NULL) {
	double TSS = own_dot_product(ri->y);

	ri->R2->val[0] = 1.0 - SSR/TSS;
    }

 bailout:

    if (!err) {
	*pV = V;
    } else {
	gretl_matrix_free(V);
    }

    gretl_matrix_block_destroy(MB);

    return err;
}

/* Variant of SVD ridge that just computes the
   parameter vector */

static int svd_ridge_bhat (double *lam, int nlam, gretl_matrix *X,
			   gretl_matrix *y, gretl_matrix *B,
			   gretl_matrix *R2, gretl_matrix *edf)
{
    gretl_matrix_block *MB = NULL;
    gretl_matrix *U = NULL;
    gretl_matrix *Vt = NULL;
    gretl_matrix *sv = NULL;
    gretl_matrix *sve = NULL;
    gretl_matrix *Uty = NULL;
    gretl_matrix *L = NULL;
    gretl_matrix *yh = NULL;
    gretl_matrix *b = NULL;
    double vij, ui, SSR;
    double edfl, TSS = 0;
    double *targ;
    int offset = 0;
    int n = X->rows;
    int k = X->cols;
    int i, j, l;
    int err;

#if RIDGE_DEBUG
    fprintf(stderr, "*** svd_ridge_bhat ***\n");
#endif

    err = gretl_matrix_SVD(X, &U, &sv, &Vt, 0);

    if (!err) {
	MB = gretl_matrix_block_new(&sve, 1, sv->cols,
				    &Uty, U->cols, 1,
				    &L, Vt->cols, Vt->rows,
				    &b, k, 1,
				    &yh, n, 1, NULL);
	if (MB == NULL) {
	    err = E_ALLOC;
	}
    }
    if (err) {
	goto bailout;
    }

    if (R2 != NULL) {
	for (i=0; i<n; i++) {
	    TSS += y->val[i] * y->val[i];
	}
    }

    offset = B->rows > k ? 1 : 0;

    gretl_matrix_multiply_mod(U, GRETL_MOD_TRANSPOSE,
			      y, GRETL_MOD_NONE,
			      Uty, GRETL_MOD_NONE);

    for (l=0; l<nlam; l++) {
	edfl = 0;
	for (j=0; j<sv->cols; j++) {
	    sve->val[j] = sv->val[j] / (sv->val[j] * sv->val[j] + lam[l]);
	    if (edf != NULL) {
		edfl += sv->val[j] * sve->val[j];
	    }
	}
	if (edf != NULL) {
	    edf->val[l] = edfl;
	}
	/* L = Vt' .* sve */
	for (j=0; j<L->cols; j++) {
	    for (i=0; i<L->rows; i++) {
		vij = gretl_matrix_get(Vt, j, i);
		gretl_matrix_set(L, i, j, vij * sve->val[j]);
	    }
	}
	gretl_matrix_multiply(L, Uty, b);
	gretl_matrix_multiply(X, b, yh);
	if (R2 != NULL) {
	    SSR = 0.0;
	    for (i=0; i<n; i++) {
		ui = y->val[i] - yh->val[i];
		SSR += ui * ui;
	    }
	    R2->val[l] = 1.0 - SSR/TSS;
	}
	targ = B->val + l * B->rows + offset;
	memcpy(targ, b->val, k * sizeof *targ);
    }

 bailout:

    gretl_matrix_block_destroy(MB);
    gretl_matrix_free(U);
    gretl_matrix_free(sv);
    gretl_matrix_free(Vt);

    return err;
}

/* called only from svd_ridge() */

static double ridge_scale (regls_info *ri,
			   gretl_matrix *lam)
{
    double lmax = NADBL;
    int i;

#if RIDGE_DEBUG
    fprintf(stderr, "*** ridge_scale, lamscale = %d ***\n",
	    ri->lamscale);
#endif

    if (ri->lamscale == LAMSCALE_GLMNET) {
	gretl_matrix *Xty = gretl_matrix_alloc(ri->X->cols, 1);

	if (Xty == NULL) {
	    return lmax;
	} else if (ri->nlam == 1) {
	    gretl_matrix_multiply_mod(ri->X, GRETL_MOD_TRANSPOSE,
				      ri->y, GRETL_MOD_NONE,
				      Xty, GRETL_MOD_NONE);
	    lmax = 1000 * vector_infnorm(Xty);
	} else {
	    /* as per glmnet, scale data by sqrt(1/n) */
	    ccd_scale(ri->X, ri->y->val, Xty->val, NULL);
	    lmax = 1000 * vector_infnorm(Xty);
	    for (i=0; i<ri->nlam; i++) {
		lam->val[i] *= lmax;
	    }
	    if (ri->nlam > 1) {
		lam->val[0] = BIG_LAMBDA;
	    }
	    gretl_matrix_free(Xty);
	}
    } else {
	/* max = squared Frobenius norm = X->cols */
	lmax = ri->X->cols;
	for (i=0; i<ri->nlam; i++) {
	    lam->val[i] *= lmax;
	}
    }

    return lmax;
}

static int svd_ridge (regls_info *ri)
{
    gretl_matrix *B = NULL;
    gretl_matrix *lam = NULL;
    gretl_matrix *V = NULL;
    double lmax = 1.0;
    double lam0 = 0.0;
    int err = 0;

#if RIDGE_DEBUG
    fprintf(stderr, "\n*** svd_ridge ***\n");
#endif

    lam = gretl_matrix_copy(ri->lfrac);
    B = gretl_zero_matrix_new(ri->k + ri->stdize, ri->nlam);
    if (lam == NULL || B == NULL) {
	return E_ALLOC;
    }

    if (ri->lamscale != LAMSCALE_NONE) {
	lmax = ridge_scale(ri, lam);
    }

#if RIDGE_DEBUG
    fprintf(stderr, "lfrac[0] = %g, lmax = %g\n", ri->lfrac->val[0], lmax);
#endif

    if (ri->nlam == 1) {
	/* calculate the covariance matrix */
	lam0 = ri->lfrac->val[0] * lmax;
	err = svd_ridge_vcv(ri, lam0, B, &V);
    } else {
	/* just calculate the parameters */
	err = svd_ridge_bhat(lam->val, ri->nlam, ri->X, ri->y, B,
			     ri->R2, ri->edf);
    }
    if (err) {
	goto bailout;
    }

    if (ri->lamscale == LAMSCALE_GLMNET) {
	/* not entirely truthful! */
	lam->val[0] = ri->lfrac->val[0] * lmax;
	if (ri->nlam == 1) {
	    lam->val[0] /= ri->n;
	}
    }

    if (!ri->xvalid) {
	ccd_get_crit(B, lam, ri);
	if (ri->verbose) {
	    ridge_print(lam, ri);
	}
	if (ri->nlam > 1) {
	    double BICmin = 1e200;
	    int j, idxmin = 0;

	    for (j=0; j<ri->nlam; j++) {
		if (ri->BIC->val[j] < BICmin) {
		    BICmin = ri->BIC->val[j];
		    idxmin = j;
		}
	    }
	    gretl_bundle_set_scalar(ri->b, "idxmin", idxmin + 1);
	    gretl_bundle_set_scalar(ri->b, "lfmin", ri->lfrac->val[idxmin]);
	}
	regls_set_crit_data(ri);
    }

    if (!err) {
	gretl_bundle_donate_data(ri->b, "B", B, GRETL_TYPE_MATRIX, 0);
	B = NULL;
	if (ri->lamscale != LAMSCALE_NONE) {
	    gretl_bundle_set_scalar(ri->b, "lmax", lmax * ri->n);
	}
	if (ri->nlam == 1) {
	    gretl_bundle_set_scalar(ri->b, "lambda", lam0);
	    // ri->lam->val[0] /= ri->n;
	    if (V != NULL) {
		gretl_bundle_donate_data(ri->b, "vcv", V,
					 GRETL_TYPE_MATRIX, 0);
	    }
	}
    }

 bailout:

    gretl_matrix_free(B);
    gretl_matrix_free(lam);

    return err;
}

/* This function is executed when we want to obtain a set
   of coefficients using the full training data, with either
   a single value of lambda or a vector of lambdas. We come
   here straight away if the user has not requested cross
   validation; we also come here after cross validation.
*/

static int admm_lasso (regls_info *ri)
{
    gretl_matrix_block *MB;
    double BICmin = 1e200;
    gretl_matrix *B = NULL;
    double lmax, rho = ri->rho;
    double llc = 0;
    int k = ri->k;
    int n = ri->n;
    int i, j, ldim;
    int idxmin = 0;
    int err = 0;

    gretl_vector *v, *u, *b, *r, *bprev, *bdiff;
    gretl_vector *q, *n1;
    gretl_matrix *L;

    ldim = n >= k ? k : n;
    MB = gretl_matrix_block_new(&v, k, 1,
				&u, k, 1,
				&b, k, 1,
				&r, k, 1,
				&bprev, k, 1,
				&bdiff, k, 1,
				&q, k, 1,
				&n1, n, 1,
				&L, ldim, ldim,
				NULL);
    if (MB == NULL) {
	return E_ALLOC;
    }
    gretl_matrix_block_zero(MB);

    lmax = ri->infnorm;

    if (!ri->xvalid && ri->verbose > 0) {
	lasso_lambda_report(ri);
    }

    if (!err) {
	get_cholesky_factor(ri->X, L, rho);
	B = make_coeff_matrix(ri);
	if (B == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	gretl_matrix_block_destroy(MB);
	return err;
    }

    if (!ri->xvalid && ri->verbose > 0) {
	lambda_sequence_header(ri->prn);
	llc = -0.5 * n * (1 + LN_2_PI - log(n));
    }

    for (j=0; j<ri->nlam && !err; j++) {
	/* loop across lambda values */
	double critj, lambda = ri->lfrac->val[j] * lmax;
	int tune_rho = 1;
	int iters = 0;
	int nnz = 0;

	err = admm_iteration(ri->X, ri->Xty, L, v, b, u, q, n1, r,
			     bprev, bdiff, lambda, &rho, tune_rho,
			     &iters);

	if (!err) {
	    for (i=0; i<k; i++) {
		if (b->val[i] != 0.0) {
		    nnz++;
		}
		if (B->cols == 1) {
		    gretl_matrix_set(B, i + ri->stdize, 0, b->val[i]);
		} else {
		    gretl_matrix_set(B, i + ri->stdize, j, b->val[i]);
		}
	    }
	    if (!ri->xvalid) {
		double R2, SSR, ll;

		critj = lasso_objective(ri->X, ri->y, b, lambda, n1, &SSR, &R2);
		ll = llc - 0.5 * n * log(SSR);
		ri->BIC->val[j] = -2 * ll + nnz * log(n);
		if (ri->verbose > 0) {
		    pprintf(ri->prn, "%12f  %5d    %f   %.4f  %#g\n",
			    lambda/n, nnz, critj, R2, ri->BIC->val[j]);
		}
		if (ri->BIC->val[j] < BICmin) {
		    BICmin = ri->BIC->val[j];
		    idxmin = j;
		}
		ri->crit->val[j] = critj;
		ri->R2->val[j] = R2;
	    }
	}
    }

    gretl_bundle_set_scalar(ri->b, "lmax", lmax);

    if (!ri->xvalid) {
	if (ri->nlam > 1) {
	    gretl_bundle_set_scalar(ri->b, "idxmin", idxmin + 1);
	    gretl_bundle_set_scalar(ri->b, "lfmin", ri->lfrac->val[idxmin]);
	}
	regls_set_crit_data(ri);
    }
    if (ri->nlam == 1) {
	gretl_bundle_set_scalar(ri->b, "lambda", ri->lfrac->val[0] * lmax);
    }

    gretl_matrix_block_destroy(MB);

    return err;
}

static int admm_do_fold (const gretl_matrix *X,
			 const gretl_matrix *y,
			 const gretl_matrix *X_out,
			 const gretl_matrix *y_out,
			 const gretl_matrix *lfrac,
			 gretl_matrix *XVC,
			 double lmax, double rho0,
			 int fold)
{
    static gretl_vector *v, *u, *b;
    static gretl_vector *r, *bprev, *bdiff;
    static gretl_vector *q, *Xty, *n1, *L;
    static gretl_matrix_block *MB;
    double rho = rho0;
    int ldim, nlam;
    int n, k, j;
    int err = 0;

    if (X == NULL) {
	/* cleanup signal */
	gretl_matrix_block_destroy(MB);
	MB = NULL;
	return 0;
    }

    nlam = gretl_vector_get_length(lfrac);
    n = X->rows;
    k = X->cols;
    ldim = n >= k ? k : n;

    if (MB == NULL) {
	MB = gretl_matrix_block_new(&v, k, 1, &u, k, 1,
				    &b, k, 1, &r, k, 1,
				    &bprev, k, 1, &bdiff, k, 1,
				    &q, k, 1, &n1, n, 1,
				    &Xty, k, 1, &L, ldim, ldim,
				    NULL);
	if (MB == NULL) {
	    return E_ALLOC;
	}
	gretl_matrix_block_zero(MB);
    }

    /* compute X'y for the estimation sample */
    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
			      y, GRETL_MOD_NONE,
			      Xty, GRETL_MOD_NONE);

    get_cholesky_factor(X, L, rho);

    for (j=0; j<nlam && !err; j++) {
	/* loop across lambda values */
	double score, lambda = lfrac->val[j] * lmax;
	int tune_rho = 1;
	int iters = 0;

	err = admm_iteration(X, Xty, L, v, b, u, q, n1, r, bprev, bdiff,
			     lambda, &rho, tune_rho, &iters);

	if (!err) {
	    /* record out-of-sample criterion */
	    gretl_matrix_reuse(n1, X_out->rows, 1);
	    score = xv_score(X_out, y_out, b, n1);
	    gretl_matrix_reuse(n1, n, 1);
	    gretl_matrix_set(XVC, j, fold, score);
	}
    }

    return err;
}

static int ccd_do_fold (gretl_matrix *X,
			gretl_matrix *y,
			gretl_matrix *X_out,
			gretl_matrix *y_out,
			const gretl_matrix *lam,
			gretl_matrix *XVC,
			int fold,
			double alpha)
{
    static gretl_matrix_block *MB;
    static gretl_matrix *Xty, *xv;
    static gretl_matrix *B;
    static gretl_matrix *u;
    static gretl_matrix *b;
    static int *ia, *nnz;
    int maxit = CCD_MAX_ITER;
    int nlp = 0, lmu = 0;
    int nlam, nout;
    int k, j;
    int err = 0;

    if (X == NULL) {
	/* cleanup signal */
	gretl_matrix_block_destroy(MB);
	MB = NULL;
	free(ia);
	ia = NULL;
	return 0;
    }

    /* dimensions */
    nlam = gretl_vector_get_length(lam);
    nout = X_out->rows;
    k = X->cols;

    if (MB == NULL) {
	MB = gretl_matrix_block_new(&xv, k, 1, &Xty, k, 1,
				    &B, k, nlam, &u, nout, 1,
				    &b, k, 1, NULL);
	ia = calloc(k + nlam, sizeof *ia);
	if (MB == NULL || ia == NULL) {
	    return E_ALLOC;
	}
	nnz = ia + k;
    }
    gretl_matrix_zero(B);

#if LAMBDA_DEBUG
    gretl_matrix_print(lam, "lam, in ccd_do_fold");
#endif

    /* scale the estimation subset by sqrt(1/n) */
    ccd_scale(X, y->val, Xty->val, xv->val);

    err = ccd_iteration(alpha, X, Xty->val, nlam, lam->val,
			ccd_toler, maxit, xv->val, &lmu, B,
			ia, nnz, NULL, &nlp);

    if (err) {
	fprintf(stderr, "ccd_do_fold: ccd_iteration returned %d\n", err);
    } else {
	/* record out-of-sample criteria */
	size_t bsize = k * sizeof(double);
	double score;

	for (j=0; j<nlam; j++) {
	    memcpy(b->val, B->val + j*k, bsize);
	    score = xv_score(X_out, y_out, b, u);
	    gretl_matrix_set(XVC, j, fold, score);
	}
    }

    return err;
}

static int svd_do_fold (gretl_matrix *X,
			gretl_matrix *y,
			gretl_matrix *X_out,
			gretl_matrix *y_out,
			const gretl_matrix *lam,
			gretl_matrix *XVC,
			int fold,
			gint8 lamscale)
{
    static gretl_matrix_block *MB;
    static gretl_matrix *B;
    static gretl_matrix *u;
    static gretl_matrix *b;
    int nlam, nout;
    int k, j;
    int err = 0;

    if (X == NULL) {
	/* cleanup signal */
	gretl_matrix_block_destroy(MB);
	MB = NULL;
	return 0;
    }

    nlam = gretl_vector_get_length(lam);
    nout = X_out->rows;
    k = X->cols;

    if (MB == NULL) {
	MB = gretl_matrix_block_new(&B, k, nlam, &u, nout, 1,
				    &b, k, 1, NULL);
	if (MB == NULL) {
	    return E_ALLOC;
	}
    }
    gretl_matrix_zero(B);

    if (lamscale == LAMSCALE_GLMNET) {
	/* scale the estimation sample by sqrt(1/n) */
	ccd_scale(X, y->val, NULL, NULL);
    }

    err = svd_ridge_bhat(lam->val, nlam, X, y, B, NULL, NULL);

#if 0
    fprintf(stderr, "svd: err=%d, nlp=%d, lmu=%d\n", err, nlp, lmu);
#endif

    if (!err) {
	/* record out-of-sample criteria */
	size_t bsize = k * sizeof(double);
	double score;

	for (j=0; j<nlam; j++) {
	    memcpy(b->val, B->val + j*k, bsize);
	    score = xv_score(X_out, y_out, b, u);
	    gretl_matrix_set(XVC, j, fold, score);
	}
    }

    return err;
}

/* Note: @X and @y are the full data matrices. @Xe and @ye will hold
   the estimation sample, and @Xf and @yf will hold the data for
   which prediction is to be performed. We need to be careful not to
   write values out of bounds in case the two disjoint sub-samples do
   not exhaust the full data (i.e. the full number of observations is
   not divisible by the fold size without remainder).
*/

static void prepare_xv_data (const gretl_matrix *X,
			     const gretl_matrix *y,
			     gretl_matrix *Xe,
			     gretl_matrix *ye,
			     gretl_matrix *Xf,
			     gretl_matrix *yf,
			     int f)
{
    int i, j, re, rf;
    double xij;

    for (j=0; j<X->cols; j++) {
	re = rf = 0;
	for (i=0; i<X->rows; i++) {
	    xij = gretl_matrix_get(X, i, j);
	    if (i/Xf->rows == f) {
		/* "out of sample" range */
		if (rf < Xf->rows) {
		    gretl_matrix_set(Xf, rf, j, xij);
		    if (j == 0) {
			yf->val[rf] = y->val[i];
		    }
		}
		rf++;
	    } else {
		/* estimation sample */
		if (re < Xe->rows) {
		    gretl_matrix_set(Xe, re, j, xij);
		    if (j == 0) {
			ye->val[re] = y->val[i];
		    }
		}
		re++;
	    }
	}
    }
}

/* Given @XVC holding criterion values per lambda (rows) and per fold
   (columns), compose a matrix holding the means, plus standard errors
   if wanted.
*/

static gretl_matrix *process_xv_criterion (gretl_matrix *XVC,
					   gretl_matrix *lfrac,
					   int *imin, int *i1se,
					   PRN *prn)
{
    gretl_matrix *metrics;
    double avg, d, v, se, se1, avgmin = 1e200;
    int mcols = 2;
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
	for (j=0; j<nf; j++) {
	    d = gretl_matrix_get(XVC, i, j) - avg;
	    v += d * d;
	}
	v /= (nf - 1);
	se = sqrt(v/nf);
	gretl_matrix_set(metrics, i, 1, se);
    }

    *i1se = *imin;

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

    if (prn != NULL) {
	int common = (*i1se == *imin);

	pprintf(prn, "          s        %s         se\n", "MSE");

	for (i=0; i<XVC->rows; i++) {
	    avg = gretl_matrix_get(metrics, i, 0);
	    se = gretl_matrix_get(metrics, i, 1);
	    pprintf(prn, "%11f %10f %10f", lfrac->val[i], avg, se);
	    if (i == *imin && common) {
		pputs(prn, " *+");
	    } else if (i == *imin) {
		pputs(prn, " *");
	    } else if (i == *i1se) {
		pputs(prn, " +");
	    }
	    pputc(prn, '\n');
	}
    }

    return metrics;
}

/* Analyse and record results after cross-validation */

static int post_xvalidation_task (regls_info *ri,
				  gretl_matrix *XVC,
				  PRN *prn)
{
    gretl_matrix *metrics;
    int imin = 0, i1se = 0;
    char **S = NULL;

    metrics = process_xv_criterion(XVC, ri->lfrac, &imin, &i1se, prn);
    if (metrics == NULL) {
	return E_ALLOC;
    }

    if (prn != NULL) {
	pputs(prn, _("\nNote: s = lambda/lambda-max\n"));
	pprintf(prn, _("Average out-of-sample %s minimized at %#g for s=%#g (\"*\")\n"),
		"MSE", gretl_matrix_get(metrics, imin, 0), ri->lfrac->val[imin]);
	pprintf(prn, _("Largest s within one s.e. of best criterion: %#g (\"+\")\n"),
		ri->lfrac->val[i1se]);
    }

    S = strings_array_new(2);
    S[0] = gretl_strdup("mean_MSE");
    S[1] = gretl_strdup("se_MSE");
    gretl_matrix_set_colnames(metrics, S);

    gretl_bundle_donate_data(ri->b, "crit", metrics, GRETL_TYPE_MATRIX, 0);
    gretl_bundle_set_int(ri->b, "idxmin", imin + 1);
    gretl_bundle_set_int(ri->b, "idx1se", i1se + 1);
    gretl_bundle_set_scalar(ri->b, "lfmin", ri->lfrac->val[imin]);
    gretl_bundle_set_scalar(ri->b, "lf1se", ri->lfrac->val[i1se]);

    return 0;
}

/* called by the cross-validation driver functions, regls_xv()
   and real_regls_xv_mpi(), for all algorithms: ADMM, CCD,
   SVD.
*/

static double get_xvalidation_lmax (regls_info *ri, int esize)
{
    double lmax = ri->infnorm;

#if LAMBDA_DEBUG
    fprintf(stderr, "get_xvalidation_lmax: ri->infnorm = %g\n", lmax);
#endif

    if (ri->ccd) {
	if (ri->alpha < 1.0) {
	    lmax /= max(ri->alpha, 1.0e-3);
#if LAMBDA_DEBUG
	    fprintf(stderr, "revised lmax = %g\n", lmax);
#endif
	}
    } else if (ri->ridge && ri->lamscale == LAMSCALE_GLMNET) {
	if (ri->alpha < 1.0) {
	    lmax /= max(ri->alpha, 1.0e-3);
	}
    } else if (ri->ridge && ri->lamscale == LAMSCALE_FROB) {
	lmax = ri->X->cols; /* ?? */
    }

    return lmax;
}

static void xv_cleanup (regls_info *ri)
{
    if (ri->ccd) {
	ccd_do_fold(NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
    } else if (ri->ridge) {
	svd_do_fold(NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
    } else {
	admm_do_fold(NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
    }
}

static gretl_matrix *make_xv_lambda (regls_info *ri,
				     double lmax,
				     int *err)
{
    gretl_matrix *lam;
    int i;

    lam = gretl_matrix_copy(ri->lfrac);
    if (lam == NULL) {
	*err = E_ALLOC;
    } else if (ri->lamscale != LAMSCALE_NONE) {
	for (i=0; i<ri->nlam; i++) {
	    lam->val[i] *= lmax;
	}
	if (ri->alpha < 1 && ri->lamscale == LAMSCALE_GLMNET) {
	    lam->val[0] = BIG_LAMBDA;
	}
    }

    return lam;
}

/* unified cross validation function, employed when we're
   not doing MPI
*/

static int regls_xv (regls_info *ri)
{
    PRN *prn = ri->prn;
    gretl_matrix_block *XY;
    gretl_matrix *Xe, *Xf;
    gretl_matrix *ye, *yf;
    gretl_matrix *lam = NULL;
    gretl_matrix *XVC = NULL;
    double lmax;
    int f, fsize, csize;
    int err = 0;

    /* the size of each fold */
    fsize = ri->n / ri->nf;
    /* the size of the complement of each fold */
    csize = (ri->nf - 1) * fsize;

    if (ri->verbose) {
	pprintf(prn, "regls_xv: nf=%d, fsize=%d, randfolds=%d, "
		"ridge=%d, ccd=%d\n", ri->nf, fsize, ri->randfolds,
		ri->ridge, ri->ccd);
	gretl_flush(prn);
    }

    XY = gretl_matrix_block_new(&Xe, csize, ri->k,
				&Xf, fsize, ri->k,
				&ye, csize, 1,
				&yf, fsize, 1, NULL);
    if (XY == NULL) {
	return E_ALLOC;
    }

    lmax = get_xvalidation_lmax(ri, csize);
    if (ri->verbose) {
	pprintf(prn, "cross-validation lmax = %g\n\n", lmax);
	gretl_flush(prn);
    }

    if (ri->ccd || ri->ridge) {
	lam = make_xv_lambda(ri, lmax, &err);
    }

    if (!err && ri->randfolds) {
	/* scramble the row order of X and y */
	randomize_rows(ri->X, ri->y);
    }

    if (!err) {
	XVC = gretl_zero_matrix_new(ri->nlam, ri->nf);
	if (XVC == NULL) {
	    err = E_ALLOC;
	}
    }

    for (f=0; f<ri->nf && !err; f++) {
	prepare_xv_data(ri->X, ri->y, Xe, ye, Xf, yf, f);
	if (ri->ccd) {
	    err = ccd_do_fold(Xe, ye, Xf, yf, lam, XVC, f, ri->alpha);

	} else if (ri->ridge) {
	    err = svd_do_fold(Xe, ye, Xf, yf, lam, XVC, f,
			      ri->lamscale);
	} else {
	    err = admm_do_fold(Xe, ye, Xf, yf, ri->lfrac, XVC,
			       lmax, ri->rho, f);
	}
    }

    /* send deallocation signal */
    xv_cleanup(ri);

    if (!err) {
	PRN *myprn = ri->verbose ? prn : NULL;

	err = post_xvalidation_task(ri, XVC, myprn);
	if (!err) {
	    /* determine coefficient vector(s) on full training set */
	    if (ri->ccd) {
		err = ccd_regls(ri);
	    } else if (ri->ridge) {
		err = svd_ridge(ri);
	    } else {
		err = admm_lasso(ri);
	    }
	}
    }

    gretl_matrix_free(lam);
    gretl_matrix_free(XVC);
    gretl_matrix_block_destroy(XY);

    return err;
}

#ifdef HAVE_MPI

static int real_regls_xv_mpi (regls_info *ri)
{
    gretl_matrix_block *XY = NULL;
    gretl_matrix *XVC = NULL;
    gretl_matrix *Xe = NULL;
    gretl_matrix *Xf = NULL;
    gretl_matrix *ye = NULL;
    gretl_matrix *yf = NULL;
    gretl_matrix *lam = NULL;
    double lmax;
    int fsize, csize;
    int folds_per;
    int folds_rem;
    int rank;
    int np, rankmax = 0;
    int f, r;
    int my_f = 0;
    int err = 0;
    PRN *prn = ri->prn;

    rank = gretl_mpi_rank();
    np = gretl_mpi_n_processes();
    rankmax = np - 1;

    fsize = ri->n / ri->nf;
    csize = (ri->nf - 1) * fsize;
    folds_per = ri->nf / np;
    folds_rem = ri->nf % np;

    /* matrix-space for per-fold data */
    XY = gretl_matrix_block_new(&Xe, csize, ri->k,
				&Xf, fsize, ri->k,
				&ye, csize, 1,
				&yf, fsize, 1, NULL);
    if (XY == NULL) {
	return E_ALLOC;
    }

    if (rank == 0) {
	lmax = get_xvalidation_lmax(ri, csize);
    }
    gretl_mpi_bcast(&lmax, GRETL_TYPE_DOUBLE, 0);

    if (ri->randfolds) {
	/* generate the same random folds in all processes */
	unsigned seed;

	if (rank == 0) {
	    if (gretl_bundle_has_key(ri->b, "seed")) {
		seed = gretl_bundle_get_unsigned(ri->b, "seed", NULL);
	    } else {
		seed = gretl_rand_get_seed();
	    }
	}
	gretl_mpi_bcast(&seed, GRETL_TYPE_UNSIGNED, 0);
	gretl_rand_set_seed(seed);
	randomize_rows(ri->X, ri->y);
    }

    if (rank < folds_rem) {
	XVC = gretl_zero_matrix_new(ri->nlam, folds_per + 1);
    } else {
	XVC = gretl_zero_matrix_new(ri->nlam, folds_per);
    }
    if (XVC == NULL) {
	err = E_ALLOC;
    }

    if (ri->ccd || ri->ridge) {
	lam = make_xv_lambda(ri, lmax, &err);
    }

    if (rank == 0) {
	if (ri->verbose) {
	    pprintf(prn, "regls_xv_mpi: nf=%d, fsize=%d, randfolds=%d\n\n",
		    ri->nf, fsize, ri->randfolds);
	    gretl_flush(prn);
	}
    }

    /* process all folds */
    r = 0;
    for (f=0; f<ri->nf && !err; f++) {
	if (rank == r) {
	    prepare_xv_data(ri->X, ri->y, Xe, ye, Xf, yf, f);
	    if (ri->verbose > 1) {
		pprintf(ri->prn, "rank %d: taking fold %d\n", rank, f+1);
	    }
	    if (ri->ccd) {
		err = ccd_do_fold(Xe, ye, Xf, yf, lam, XVC, my_f++,
				  ri->alpha);
	    } else if (ri->ridge) {
		err = svd_do_fold(Xe, ye, Xf, yf, lam, XVC, my_f++,
				  ri->lamscale);
	    } else {
		err = admm_do_fold(Xe, ye, Xf, yf, ri->lfrac, XVC, lmax,
				   ri->rho, my_f++);
	    }
	}
	if (r == rankmax) {
	    r = 0;
	} else {
	    r++;
	}
    }

    /* reduce @XVC to root by column concatenation */
    gretl_matrix_mpi_reduce(XVC, &XVC, GRETL_MPI_HCAT, 0, OPT_NONE);

    /* send deallocation signal, all processes */
    xv_cleanup(ri);

    if (rank == 0 && !err) {
	PRN *myprn = ri->verbose ? prn : NULL;

	err = post_xvalidation_task(ri, XVC, myprn);
	if (!err) {
	    /* determine coefficient vector on full training set */
	    if (ri->ccd) {
		err = ccd_regls(ri);
	    } else if (ri->ridge) {
		err = svd_ridge(ri);
	    } else {
		err = admm_lasso(ri);
	    }
	}
    }

    gretl_matrix_free(lam);
    gretl_matrix_free(XVC);
    gretl_matrix_block_destroy(XY);

    return err;
}

static int xv_use_mpi (regls_info *ri)
{
    int no_mpi = gretl_bundle_get_bool(ri->b, "no_mpi", 0);
    int ret = (no_mpi == 0);

    /* It's not yet clear whether MPI is useful for the ridge case, or
       for ccd in general. This may depend on the size of the data;
       experimentation is needed!
    */
    if (ret && (ri->ccd || ri->ridge)) {
	ret = 0;
    }

    return ret;
}

#endif /* HAVE_MPI or not */

int gretl_regls (gretl_matrix *X,
		 gretl_matrix *y,
		 gretl_bundle *bun,
		 PRN *prn)
{
    int (*regfunc) (regls_info *) = NULL;
    regls_info *ri;
    int err = 0;

    ri = regls_info_new(X, y, bun, prn, &err);
    if (err) {
	fprintf(stderr, "err %d from regls_info_new\n", err);
	return err;
    }

    if (ri->xvalid) {
#ifdef HAVE_MPI
	if (xv_use_mpi(ri)) {
	    if (gretl_mpi_n_processes() > 1) {
		regfunc = real_regls_xv_mpi;
	    } else if (auto_mpi_ok()) {
		regfunc = mpi_parent_action;
	    }
	}
#endif
	if (regfunc == NULL) {
	    regfunc = regls_xv;
	}
    } else if (ri->ccd) {
	regfunc = ccd_regls;
    } else if (ri->ridge) {
	regfunc = svd_ridge;
    } else {
	regfunc = admm_lasso;
    }

#ifdef HAVE_MPI
    if (regfunc != mpi_parent_action) {
	err = regls_set_Xty(ri);
    }
#else
    err = regls_set_Xty(ri);
#endif

    if (!err) {
	err = regfunc(ri);
    }

    regls_info_free(ri);

    return err;
}

#ifdef HAVE_MPI

/* We come here if a parent process has called our automatic local MPI
   routine for cross validation: this function will be executed by all
   gretlmpi instances.
*/

int regls_xv_mpi (PRN *prn)
{
    regls_info *ri = NULL;
    gretl_bundle *bun = NULL;
    gretl_matrix *X;
    gretl_matrix *y;
    int err = 0;

    /* read matrices deposited by parent process */
    X = gretl_matrix_read_from_file("regls_X.bin", 1, &err);
    y = gretl_matrix_read_from_file("regls_y.bin", 1, &err);

    if (!err) {
	bun = gretl_bundle_read_from_file("regls_bun.xml", 1, &err);
    }

    if (!err) {
	ri = regls_info_new(X, y, bun, prn, &err);
    }
    if (!err) {
	err = regls_set_Xty(ri);
    }

    if (!err) {
	err = real_regls_xv_mpi(ri);
	if (!err && gretl_mpi_rank() == 0) {
	    /* write results, to be picked up by parent */
	    gretl_bundle_write_to_file(bun, "regls_XV_result.xml", 1);
	}
    }

    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_bundle_destroy(bun);
    regls_info_free(ri);

    return err;
}

static int mpi_parent_action (regls_info *ri)
{
    int err;

    err = gretl_matrix_write_to_file(ri->X, "regls_X.bin", 1);
    if (!err) {
	err = gretl_matrix_write_to_file(ri->y, "regls_y.bin", 1);
    }
    if (!err) {
	err = gretl_bundle_write_to_file(ri->b, "regls_bun.xml", 1);
    }

    if (!err) {
	/* compose and execute MPI script */
	err = foreign_start(MPI, NULL, OPT_NONE, ri->prn);
	if (!err) {
	    int np = gretl_bundle_get_int(ri->b, "mpi_np", NULL);
	    int mpi_local = gretl_bundle_get_int(ri->b, "mpi_local", NULL);
	    gretlopt mpi_opt = OPT_S | OPT_Q;

	    if (np > 0) {
		/* user-specified number of processes */
		mpi_opt |= OPT_N;
		set_optval_int(MPI, OPT_N, np);
	    }
	    if (mpi_local) {
		/* local machine only */
		mpi_opt |= OPT_L;
	    }
	    if (ri->verbose) {
		pputs(ri->prn, _("Invoking MPI...\n\n"));
		gretl_flush(ri->prn);
	    } else {
		fprintf(stderr, "doing MPI\n");
	    }
	    foreign_append("_regls()", MPI);
	    err = foreign_execute(NULL, mpi_opt, ri->prn);
	    if (err) {
		fprintf(stderr, "mpi_parent: foreign exec error %d\n", err);
	    }
	}
    }

    if (!err) {
	/* retrieve results bundle written by gretlmpi */
	gretl_bundle *res;

	res = gretl_bundle_read_from_file("regls_XV_result.xml", 1, &err);
	if (!err) {
	    gretl_bundles_swap_content(ri->b, res);
	    gretl_bundle_destroy(res);
	}
    }

    return err;
}

#endif /* HAVE_MPI */

/* apparatus to support glasso ("graphical lasso") */

static int any_missing (const gretl_matrix *X)
{
    int i, n = X->rows * X->cols;

    for (i=0; i<n; i++) {
        if (na(X->val[i])) {
            return 1;
        }
    }

    return 0;
}

static void Mj_matrix (gretl_matrix *M_,
                       const gretl_matrix *M,
                       int cut)
{
    int i, j, k = 0;

    for (j=0; j<M->cols; j++) {
        for (i=0; i<M->rows; i++) {
            if (i != cut && j != cut) {
                M_->val[k++] = gretl_matrix_get(M, i, j);
            }
        }
    }
}

static void Mj_vector (gretl_matrix *v,
                       const gretl_matrix *M,
                       int j)
{
    int i, k = 0;

    for (i=0; i<M->rows; i++) {
        if (i != j) {
            v->val[k++] = gretl_matrix_get(M, i, j);
        }
    }
}

static void Wj_revise (gretl_matrix *W, int j,
                       const gretl_matrix *Wb)
{
    double wbij;
    int i, k = 0;

    for (i=0; i<W->rows; i++) {
        if (i != j) {
            wbij = Wb->val[k++];
            gretl_matrix_set(W, i, j, wbij);
            gretl_matrix_set(W, j, i, wbij);
        }
    }
}

static void vsqrt (gretl_matrix *v)
{
    int i;

    for (i=0; i<v->rows; i++) {
        v->val[i] = sqrt(v->val[i]);
    }
}

/* compute X = E * (v .* E') */

static void ev_dotmul1 (gretl_matrix *X,
                        const gretl_matrix *v,
                        const gretl_matrix *E,
                        gretl_matrix *Tmp)
{
    double xij;
    int i, j;

    gretl_matrix_copy_values(Tmp, E);
    for (j=0; j<Tmp->cols; j++) {
        for (i=0; i<Tmp->rows; i++) {
            xij = gretl_matrix_get(Tmp, i, j);
            gretl_matrix_set(Tmp, i, j, xij * v->val[i]);
        }
    }

    gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
                              Tmp, GRETL_MOD_NONE,
                              X, GRETL_MOD_NONE);
}

/* compute Y = [E * (1 ./ v) * E'] * R */

static void ev_dotmul2 (gretl_matrix *Y,
                        const gretl_matrix *v,
                        const gretl_matrix *E,
                        const gretl_matrix *R,
                        gretl_matrix *Tmp1,
                        gretl_matrix *Tmp2)
{
    double xij;
    int i, j;

    gretl_matrix_copy_values(Tmp1, E);
    for (j=0; j<Tmp1->cols; j++) {
        for (i=0; i<Tmp1->rows; i++) {
            xij = gretl_matrix_get(Tmp1, i, j);
            gretl_matrix_set(Tmp1, i, j, xij / v->val[i]);
        }
    }

    gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
                              Tmp1, GRETL_MOD_NONE,
                              Tmp2, GRETL_MOD_NONE);
    gretl_matrix_multiply(Tmp2, R, Y);
}

struct glasso_info_ {
    gretl_matrix_block *B;
    gretl_matrix *W;
    gretl_matrix *W0;
    gretl_matrix *Wd;
    gretl_matrix *Sj;
    gretl_matrix *Wj;
    gretl_matrix *Ev;
    gretl_matrix *X;
    gretl_matrix *y;
    gretl_matrix *Xty;
    gretl_matrix *xv;
    gretl_matrix *b;
    gretl_matrix *n1;
    gretl_matrix *Tmp1;
    gretl_matrix *Tmp2;
    double *lam;
    int *ia;
    int *nnz;
};

typedef struct glasso_info_ glasso_info;

static glasso_info *glasso_info_new (int p, int n, double rho)
{
    glasso_info *gi = malloc(sizeof *gi);

    if (gi == NULL) {
	return NULL;
    }

    gi->B = gretl_matrix_block_new(&gi->W0, p, p,
				   &gi->Wd, p, p,
				   &gi->Sj, n, 1,
				   &gi->Wj, n, n,
				   &gi->Ev, n, n,
				   &gi->X, n, n,
				   &gi->y, n, 1,
				   &gi->Xty, n, 1,
				   &gi->xv, n, 1,
				   &gi->b, n, 1,
				   &gi->n1, n, 1,
				   &gi->Tmp1, n, n,
				   &gi->Tmp2, n, n,
				   NULL);
    if (gi->B == NULL) {
	free(gi);
	return NULL;
    }

    gi->lam = malloc(sizeof *gi->lam);
    gi->lam[0] = rho;

    gi->ia = calloc(n + 1, sizeof *gi->ia);
    if (gi->ia != NULL) {
        gi->nnz = gi->ia + n;
    }

    return gi;
}

static void glasso_info_free (glasso_info *gi)
{
    if (gi != NULL) {
	gretl_matrix_block_destroy(gi->B);
	free(gi->lam);
	free(gi->ia);
	free(gi);
    }
}

static int ccd_glasso (glasso_info *gi, double tol)
{
    int maxit = CCD_MAX_ITER;
    int nlp = 0, lmu = 0;
    int err;

    ccd_toler = tol;
    ccd_scale(gi->X, gi->y->val, gi->Xty->val, gi->xv->val);
    gretl_matrix_zero(gi->b);

#if LAMBDA_DEBUG
    fprintf(stderr, "ccd_glasso: lambda = %g\n", gi->lam[0]);
#endif

    err = ccd_iteration(1.0, gi->X, gi->Xty->val, 1, gi->lam,
			ccd_toler, maxit, gi->xv->val, &lmu, gi->b, gi->ia,
			gi->nnz, NULL, &nlp);

    return err;
}

static int glasso_converged (const gretl_matrix *W0,
			     const gretl_matrix *W1,
			     double tol)
{
    const double *x = W0->val;
    const double *y = W1->val;
    double csum = 0;
    int i, j;

    for (j=0; j<W0->cols; j++) {
        csum = 0.0;
        for (i=0; i<W0->rows; i++) {
	    csum = fabs(x[i] - y[i]);
        }
        if (csum > tol) {
	    return 0;
        }
	x += W0->rows;
	y += W0->rows;
    }

#if 0
    fprintf(stderr, "glasso_converged: tol %g, csum %g\n", tol, csum);
#endif

    return 1;
}

/* The following function may be redundant if we're getting option
   values that have already been vetted by graphlasso.gfn. The @tol
   and @maxit defaults below are as in the glasso R package.
*/

static int handle_glasso_options (gretl_bundle *b,
				  double *prho,
				  double *ptol,
				  int *pmaxit)
{
    double rho, tol = 1.0e-4;
    int maxit = 1.0e4;
    int err = 0;

    if (gretl_bundle_has_key(b, "rho")) {
	rho = gretl_bundle_get_scalar(b, "rho", &err);
	if (!err && rho < 0) {
	    err = E_INVARG;
	}
    } else {
	err = E_ARGS;
    }

    if (!err && gretl_bundle_has_key(b, "tol")) {
	tol = gretl_bundle_get_scalar(b, "tol", &err);
	if (!err && tol <= 0) {
	    err = E_INVARG;
	}
    }

    if (!err && gretl_bundle_has_key(b, "maxit")) {
	maxit = gretl_bundle_get_int(b, "maxit", &err);
	if (!err && maxit <= 0) {
	    err = E_INVARG;
	}
    }

    if (!err) {
	*prho = rho;
	*ptol = tol;
	*pmaxit = maxit;
    }

    return err;
}

gretl_matrix *gretl_glasso (const gretl_matrix *S,
                            gretl_bundle *b,
                            PRN *prn, int *err)
{
    glasso_info *gi = NULL;
    gretl_matrix *W = NULL;
    gretl_matrix *dsqrt;
    double rho, tol;
    int maxit;
    double wij;
    int p = S->rows;
    int n = p-1;
    int i, j;

    if (p < 2 || S->cols != p) {
        *err = E_NONCONF;
    } else if (any_missing(S)) {
        *err = E_MISSDATA;
    } else {
	*err = handle_glasso_options(b, &rho, &tol, &maxit);
    }
    if (*err) {
        return NULL;
    }

    /* this will become the return value */
    W = gretl_matrix_copy(S);
    if (W == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    gi = glasso_info_new(p, n, rho/n);
    if (gi == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(W);
	return NULL;
    }

    for (i=0; i<p; i++) {
        wij = gretl_matrix_get(W, i, i);
        gretl_matrix_set(W, i, i, wij + rho);
    }
    gretl_matrix_copy_values(gi->W0, W);

    for (i=0; i<maxit; i++) {
        for (j=p-1; j>=0 && !*err; j--) {
            Mj_matrix(gi->Wj, W, j);
            Mj_vector(gi->Sj, S, j);
            gretl_matrix_copy_values(gi->Ev, gi->Wj);
            dsqrt = gretl_symmetric_matrix_eigenvals(gi->Ev, 1, err);
            vsqrt(dsqrt);
            gretl_square_matrix_transpose(gi->Ev);
            ev_dotmul1(gi->X, dsqrt, gi->Ev, gi->Tmp1);
            ev_dotmul2(gi->y, dsqrt, gi->Ev, gi->Sj, gi->Tmp1, gi->Tmp2);
            *err = ccd_glasso(gi, tol);
            if (!*err) {
                gretl_matrix_multiply(gi->Wj, gi->b, gi->n1);
                Wj_revise(W, j, gi->n1);
            } else {
		fprintf(stderr, "gretl_glasso: err = %d from ccd_glasso\n", *err);
	    }
	    gretl_matrix_free(dsqrt);
        }
	if (*err || glasso_converged(W, gi->W0, tol)) {
            break;
        }
        gretl_matrix_copy_values(gi->W0, W);
    }

    glasso_info_free(gi);

    if (*err) {
        gretl_matrix_free(W);
        W = NULL;
    } else {
	pprintf(prn, "glasso iterations: %d\n", i + 1);
    }

    return W;
}
