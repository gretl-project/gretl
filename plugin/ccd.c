#include "libgretl.h"

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

static double dot_product (const double *x, const double *y, int n)
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

static double dot_product (const double *x, const double *y, int n)
{
    double ret = 0.0;
    int i;

    for (i=0; i<n; i++) {
	ret += x[i] * y[i];
    }
    return ret;
}

#endif

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

static void fill_coeff_column (gretl_matrix *B, int j,
			       const double *a, const int *idx,
			       int n)
{
    double *b = B->val + j * B->rows;
    int i;

    for (i=0; i<n; i++) {
	b[i] = a[idx[i]];
    }
}

/* sign(x,y): gives "the value of x with the sign of y",
   but in context @x will always be positive.
*/

static inline double sign (double x, double y)
{
    return y >= 0 ? x : -x;
}

static int ccd_scale (int n, int k, gretl_matrix *x, double *y,
		      double *g, double *xv)
{
    double *xj, v = sqrt(1.0/n);
    int i, j;

    for (i=0; i<n; i++) {
	y[i] *= v;
    }
    for (j=0; j<k; j++) {
	xj = x->val + j * n;
	for (i=0; i<n; i++) {
	    xj[i] *= v;
	}
	xv[j] = dot_product(xj, xj, n);
	g[j] = dot_product(y, xj, n);
    }

    return 0;
}

static int ccd_iteration (const gretl_matrix *X, double *g,
			  int nlam, const double *ulam, double thr,
			  int maxit, const double *xv, int *lmu,
			  gretl_matrix *B, int *ia, int *kin,
			  double *Rsq, int *pnlp)
{
    gretl_matrix *C;
    double alm, u, v, rsq = 0;
    double ak, del, dlx, cij;
    double *a, *da;
    int *mm, nin, jz, iz = 0;
    int j, k, l, m, nlp = 0;
    int nx = X->cols;
    int err = 0;

    C = gretl_matrix_alloc(nx, nx);
    a = malloc(nx * sizeof *a);
    da = malloc(nx * sizeof *da);
    mm = malloc(nx * sizeof *mm);
    if (C == NULL || a == NULL || da == NULL || mm == NULL) {
	return E_ALLOC;
    }
    /* "zero" @a and @mm */
    for (j=0; j<nx; j++) {
	a[j] = 0.0;
	mm[j] = -1;
    }
    nin = nlp = *pnlp = 0;

    for (m=0; m<nlam; m++) {
	alm = ulam[m];
	jz = 1;
    maybe_restart:
	if (iz*jz == 0) {
            nlp++;
            dlx = 0.0;
	    for (k=0; k<nx; k++) {
		ak = a[k];
		u = g[k] + ak*xv[k];
		v = fabs(u) - alm;
		a[k] = v > 0.0 ? sign(v,u) / xv[k] : 0.0;
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
		fprintf(stderr, "ccd: max iters reached\n");
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
	    v = fabs(u) - alm;
	    a[k] = v > 0.0 ? sign(v,u)/xv[k] : 0.0;
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
	    goto getout;
	}
    m_finish:
	if (nin <= nx) {
	    if (nin > 0) {
		fill_coeff_column(B, m, a, ia, nin);
	    }
	    kin[m] = nin;
	    Rsq[m] = rsq;
	    *lmu = m + 1;
	} else {
            err = E_NOCONV;
	    fprintf(stderr, "ccd: error at foot of loop\n");
            goto getout;
	}
    } /* end loop over lambda values */

 getout:

    *pnlp = nlp;
    free(a);
    free(mm);
    free(da);
    gretl_matrix_free(C);

    return err;
}

static gretl_matrix *ccd_c_unpack (gretl_matrix *B,
				   gretl_matrix *lam,
				   const gretl_matrix *R2,
				   int *nin, int *ia,
				   int lmu, int stdize,
				   gretl_matrix *crit,
				   double nulldev,
				   int *err)
{
    gretl_matrix *BB;
    double *bj, bij;
    double bsum, SSR;
    int k = B->rows;
    int nlam = B->cols;
    int i, j;

    BB = gretl_zero_matrix_new(k + stdize, nlam);
    if (BB == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* figure lasso criterion, if wanted */
    if (crit != NULL) {
	for (j=0; j<lmu; j++) {
	    bsum = 0;
	    for (i=0; i<nin[j]; i++) {
		bij = gretl_matrix_get(B, i, j);
		bsum += fabs(bij);
	    }
	    SSR = nulldev * (1.0 - R2->val[j]);
	    bsum *= lam->val[j];
	    gretl_vector_set(crit, j, 0.5 * SSR + bsum);
	}
    }

    /* "unpack" @B into @BB */
    for (j=0; j<nlam; j++) {
	bj = B->val + j*k;
	for (i=0; i<k; i++) {
	    if (fabs(bj[i]) > 0) {
		gretl_matrix_set(BB, ia[i] + stdize, j, bj[i]);
	    }
	}
    }

    return BB;
}

static void ccd_print (const gretl_matrix *B,
		       const gretl_matrix *R2,
		       const gretl_matrix *lam,
		       const gretl_matrix *crit,
		       PRN *prn)
{
    double *bj;
    int k = B->rows;
    int nlam = B->cols;
    int i, j, dfj;

    pputc(prn, '\n');
    if (crit != NULL) {
	/* header for output showing lasso criterion */
	pputs(prn, "      lambda     df   criterion      R^2\n");
    } else {
	/* as per R, more or less */
	pputs(prn, "    df     R^2  lambda\n");
    }
    for (j=0; j<nlam; j++) {
	bj = B->val + j*k;
	dfj = 0;
	for (i=0; i<k; i++) {
	    dfj += fabs(bj[i]) > 0;
	}
	if (crit != NULL) {
	    pprintf(prn, "%12f  %5d    %f   %.4f\n",
		    lam->val[j], dfj, crit->val[j], R2->val[j]);
	} else {
	    pprintf(prn, "%-2d  %2d  %.4f  %.4f\n", j+1, dfj, R2->val[j],
		    lam->val[j]);
	}
    }
}

int ccd_driver (gretl_matrix *X, gretl_matrix *y,
		gretl_matrix *lam, int stdize,
		gretl_matrix **pB,
		gretl_matrix **pR2,
		gretl_matrix **pcrit,
		PRN *prn)
{
    gretl_matrix_block *MB;
    gretl_matrix *g, *xv;
    gretl_matrix *B;
    gretl_matrix *R2;
    gretl_matrix *crit = NULL;
    double thresh = 1.0e-7;
    double lmax;
    int maxit = 100000;
    int *nnz, *ia;
    int nlp = 0, lmu = 0;
    int nlam;
    int n, k, i;
    int err = 0;

    /* note: set @thresh to 1.0e-9 to get results as accurate as
       those from ADMM with its default tolerances
    */

    n = X->rows;
    k = X->cols;
    nlam = lam->rows;

    MB = gretl_matrix_block_new(&xv, k, 1, &g, k, 1,
				&B, k, nlam, NULL);
    if (MB == NULL) {
	return E_ALLOC;
    }
    gretl_matrix_zero(B);

    /* scale data by sqrt(1/n) */
    ccd_scale(n, k, X, y->val, g->val, xv->val);

    lmax = gretl_matrix_infinity_norm(g);
    for (i=0; i<nlam; i++) {
	lam->val[i] *= lmax;
    }

    R2 = gretl_matrix_alloc(nlam, 1);
    nnz = malloc(nlam * sizeof *nnz);
    ia = malloc(k * sizeof *ia);
    if (pcrit != NULL) {
	crit = gretl_matrix_alloc(nlam, 1);
    }

    if (R2 == NULL || nnz == NULL || ia == NULL ||
	(pcrit != NULL && crit == NULL)) {
	err = E_ALLOC;
	goto bailout;
    }

    err = ccd_iteration(X, g->val, nlam, lam->val,
			thresh, maxit, xv->val, &lmu, B,
			ia, nnz, R2->val, &nlp);
    fprintf(stderr, "ccd: err=%d, nlp=%d, lmu=%d\n", err, nlp, lmu);

    if (!err) {
	gretl_matrix *BB;
	double nulldev = 1.0;

	if (pcrit != NULL && !stdize) {
	    for (i=0; i<n; i++) {
		nulldev += y->val[i] * y->val[i];
	    }
	}
	BB = ccd_c_unpack(B, lam, R2, nnz, ia, lmu, stdize,
			  crit, nulldev, &err);
	if (!err) {
	    ccd_print(B, R2, lam, crit, prn);
	    *pB = BB;
	    *pR2 = R2;
	    if (pcrit != NULL) {
		*pcrit = crit;
	    }
	}
    }

 bailout:

    if (err) {
	gretl_matrix_free(R2);
	gretl_matrix_free(crit);
    }
    gretl_matrix_block_destroy(MB);
    free(nnz);
    free(ia);

    return 0;
}
