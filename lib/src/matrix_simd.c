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

/* All code here is conditional on AVX (128-bit SSE is not really
   worth the bother when working with doubles)
*/

#define SHOW_SIMD 0

static int gretl_matrix_simd_add_to (gretl_matrix *a,
				     const gretl_matrix *b,
				     int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int i, imax = n / 4;
    int rem = n % 4;

#if SHOW_SIMD
    fprintf(stderr, "AVX: gretl_matrix_simd_add_to (%d x %d)\n",
	    a->rows, a->cols);
#endif

    for (i=0; i<imax; i++) {
	/* add 4 doubles in parallel */
	__m256d Ymm_A = _mm256_loadu_pd(ax);
	__m256d Ymm_B = _mm256_loadu_pd(bx);
	__m256d Ymm_C = _mm256_add_pd(Ymm_A, Ymm_B);

	_mm256_storeu_pd(ax, Ymm_C);
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	ax[i] += bx[i];
    }

    return 0;
}

static int gretl_matrix_simd_subt_from (gretl_matrix *a,
					const gretl_matrix *b,
					int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int i, imax = n / 4;
    int rem = n % 4;

#if SHOW_SIMD
    fprintf(stderr, "AVX: gretl_matrix_simd_subt_from (%d x %d, rem=%d)\n",
	    a->rows, a->cols, rem);
#endif

    for (i=0; i<imax; i++) {
	/* subtract 4 doubles in parallel */
	__m256d Ymm_A = _mm256_loadu_pd(ax);
	__m256d Ymm_B = _mm256_loadu_pd(bx);
	__m256d Ymm_C = _mm256_sub_pd(Ymm_A, Ymm_B);

	_mm256_storeu_pd(ax, Ymm_C);
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	ax[i] -= bx[i];
    }

    return 0;
}

static int gretl_matrix_simd_add (const double *ax,
				  const double *bx,
				  double *cx,
				  int n)
{
    int i, imax = n / 4;
    int rem = n % 4;

#if SHOW_SIMD
    fprintf(stderr, "AVX: gretl_matrix_simd_add (n = %d\n", n);
#endif

    for (i=0; i<imax; i++) {
	/* process 4 doubles in parallel */
	__m256d Ymm_A = _mm256_loadu_pd(ax);
	__m256d Ymm_B = _mm256_loadu_pd(bx);
	__m256d Ymm_C = _mm256_add_pd(Ymm_A, Ymm_B);

	_mm256_storeu_pd(cx, Ymm_C);
	ax += 4;
	bx += 4;
	cx += 4;
    }

    for (i=0; i<rem; i++) {
	cx[i] = ax[i] + bx[i];
    }

    return 0;
}

static int gretl_matrix_simd_subtract (const double *ax,
				       const double *bx,
				       double *cx,
				       int n)
{
    int i, imax = n / 4;
    int rem = n % 4;

#if SHOW_SIMD
    fprintf(stderr, "AVX: gretl_matrix_simd_subtract (n = %d, rem = %d)\n",
	    n, rem);
#endif

    for (i=0; i<imax; i++) {
	/* process 4 doubles in parallel */
	__m256d Ymm_A = _mm256_loadu_pd(ax);
	__m256d Ymm_B = _mm256_loadu_pd(bx);
	__m256d Ymm_C = _mm256_sub_pd(Ymm_A, Ymm_B);

	_mm256_storeu_pd(cx, Ymm_C);
	cx += 4;
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	cx[i] = ax[i] - bx[i];
    }

    return 0;
}

/* very fast but restrictive: both A and B must be 4 x 4 */

static int gretl_matrix_avx_mul4 (const double *aval,
				  const double *bval,
				  double *cval)
{
    __m256d b1, b2, b3, b4;
    __m256d mul, col;
    int j;

    /* load the columns of A */
    __m256d a1 = _mm256_loadu_pd(aval);
    __m256d a2 = _mm256_loadu_pd(aval + 4);
    __m256d a3 = _mm256_loadu_pd(aval + 8);
    __m256d a4 = _mm256_loadu_pd(aval + 12);

    for (j=0; j<4; j++) {
	/* broadcast the elements of col j of B */
        b1 = _mm256_broadcast_sd(&bval[4*j]);
        b2 = _mm256_broadcast_sd(&bval[4*j + 1]);
        b3 = _mm256_broadcast_sd(&bval[4*j + 2]);
        b4 = _mm256_broadcast_sd(&bval[4*j + 3]);

	mul = _mm256_mul_pd(b1, a1);
	b1  = _mm256_mul_pd(b2, a2);
	col = _mm256_add_pd(mul, b1);
	mul = _mm256_mul_pd(b3, a3);
        b2  = _mm256_mul_pd(b4, a4);
	b3  = _mm256_add_pd(mul, b2);
	col = _mm256_add_pd(col, b3);

	_mm256_storeu_pd(&cval[4*j], col);
    }

    return 0;
}

/* very fast but restrictive: both A and B must be 8 x 8 */

static int gretl_matrix_avx_mul8 (const double *aval,
				  const double *bval,
				  double *cval)
{
    __m256d a1, a2, a3, a4, a5, a6, a7, a8;
    __m256d b1, b2, b3, b4, b5, b6, b7, b8;
    __m256d mul, col;
    int i, j;

    for (i=0; i<2; i++) {
	/* half-load of columns of A */
	a1 = _mm256_loadu_pd(aval);
	a2 = _mm256_loadu_pd(aval + 8);
	a3 = _mm256_loadu_pd(aval + 16);
	a4 = _mm256_loadu_pd(aval + 24);
	a5 = _mm256_loadu_pd(aval + 32);
	a6 = _mm256_loadu_pd(aval + 40);
	a7 = _mm256_loadu_pd(aval + 48);
	a8 = _mm256_loadu_pd(aval + 56);

	for (j=0; j<8; j++) {
	    /* broadcast the elements of col i of B */
	    b1 = _mm256_broadcast_sd(&bval[8*j]);
	    b2 = _mm256_broadcast_sd(&bval[8*j + 1]);
	    b3 = _mm256_broadcast_sd(&bval[8*j + 2]);
	    b4 = _mm256_broadcast_sd(&bval[8*j + 3]);
	    b5 = _mm256_broadcast_sd(&bval[8*j + 4]);
	    b6 = _mm256_broadcast_sd(&bval[8*j + 5]);
	    b7 = _mm256_broadcast_sd(&bval[8*j + 6]);
	    b8 = _mm256_broadcast_sd(&bval[8*j + 7]);

	    mul = _mm256_mul_pd(b1, a1);
	    b1  = _mm256_mul_pd(b2, a2);
	    col = _mm256_add_pd(mul, b1);
	    mul = _mm256_mul_pd(b3, a3);
	    b2  = _mm256_mul_pd(b4, a4);
	    b3  = _mm256_add_pd(mul, b2);
	    col = _mm256_add_pd(col, b3);
	    mul = _mm256_mul_pd(b5, a5);
	    b4  = _mm256_mul_pd(b6, a6);
	    b5  = _mm256_add_pd(mul, b4);
	    col = _mm256_add_pd(col, b5);
	    mul = _mm256_mul_pd(b7, a7);
	    b6  = _mm256_mul_pd(b8, a8);
	    b7  = _mm256_add_pd(mul, b6);
	    col = _mm256_add_pd(col, b7);

	    _mm256_storeu_pd(&cval[8*j], col);
	}
	aval += 4;
	cval += 4;
    }

    return 0;
}

/* Note: this is probably usable only for k <= 8 (shortage of AVX
   registers). It is much faster if m is a multiple of 4; n is
   unconstrained.
*/

static int gretl_matrix_simd_mul (const gretl_matrix *A,
				  const gretl_matrix *B,
				  gretl_matrix *C)
{
    int m = A->rows;
    int n = B->cols;
    int k = A->cols;
    const double *aval = A->val;
    const double *bval = B->val;
    double *cval = C->val;
    int hmax, hrem, i, j;

#if SHOW_SIMD
    fprintf(stderr, "AVX: gretl_matrix_simd_mul (MNK = %d, %d, %d)\n",
	    m, n, k);
#endif

    if (m == 4 && n == 4 && k == 4) {
	return gretl_matrix_avx_mul4(aval, bval, cval);
    }

    if (m == 8 && n == 8 && k == 8) {
	return gretl_matrix_avx_mul8(aval, bval, cval);
    }

    hmax = m / 4;
    hrem = m % 4;

    if (m >= 4) {
	__m256d a[k], b[k];
	__m256d mult, ccol;
	int h;

	for (h=0; h<hmax; h++) {
	    /* loop across the k columns of A, loading
	       4 elements from each */
	    for (j=0; j<k; j++) {
		a[j] = _mm256_loadu_pd(aval + j*m);
	    }
	    /* loop across the n columns of B */
	    for (j=0; j<n; j++) {
		/* broadcast the k elements of col j of B */
		for (i=0; i<k; i++) {
		    b[i] = _mm256_broadcast_sd(&bval[j*k + i]);
		}
		/* cumulate the products */
		ccol = _mm256_setzero_pd();
		for (i=0; i<k; i++) {
		    mult = _mm256_mul_pd(b[i], a[i]);
		    ccol = _mm256_add_pd(ccol, mult);
		}
		/* write 4 rows of each column of C */
		_mm256_storeu_pd(&cval[m*j], ccol);
	    }
	    /* advance by the number of doubles we can store in
	       one go */
	    aval += 4;
	    cval += 4;
	}
    }

    if (hrem >= 2) {
	/* do a single 128-bit run */
	__m128d a[k], b[k];
	__m128d mult, ccol;

	for (j=0; j<k; j++) {
	    a[j] = _mm_loadu_pd(aval + j*m);
	}
	for (j=0; j<n; j++) {
	    for (i=0; i<k; i++) {
		b[i] = _mm_set1_pd(bval[j*k + i]);
	    }
	    ccol = _mm_setzero_pd();
	    for (i=0; i<k; i++) {
		mult = _mm_mul_pd(b[i], a[i]);
		ccol = _mm_add_pd(ccol, mult);
	    }
	    _mm_storeu_pd(&cval[m*j], ccol);
	}
	hrem -= 2;
	aval += 2;
	cval += 2;
    }

    if (hrem) {
	/* odd-valued m: compute the last row */
	double ccol, a[k];

	for (j=0; j<k; j++) {
	    a[j] = aval[j*m];
	}
	for (j=0; j<n; j++) {
	    ccol = 0.0;
	    for (i=0; i<k; i++) {
		ccol += bval[j*k + i] * a[i];
	    }
	    cval[m*j] = ccol;
	}
    }

    return 0;
}

/* See https://stackoverflow.com/questions/49941645,
   Peter Cordes's answer on how efficiently to sum the
   contents of an __m256d into a single double.
*/

static inline double hsum_double_avx (__m256d v)
{
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    __m128d high64;

    vlow   = _mm_add_pd(vlow, vhigh);
    high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}

double gretl_vector_simd_dot_product (const gretl_vector *a,
				      const gretl_vector *b)
{
    __m256d Ymm_A, Ymm_B, Ymm_C;
    const double *ax = a->val;
    const double *bx = b->val;
    int n = gretl_vector_get_length(a);
    int i, imax = n / 4;
    int rem = n % 4;
    double ret = 0.0;

    for (i=0; i<imax; i++) {
	/* multiply 4 doubles in parallel */
	Ymm_A = _mm256_loadu_pd(ax);
	Ymm_B = _mm256_loadu_pd(bx);
	Ymm_C = _mm256_mul_pd(Ymm_A, Ymm_B);
	ret += hsum_double_avx(Ymm_C);
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	ret += ax[i] * bx[i];
    }

    return ret;
}

void gretl_matrix_simd_scalar_mul (double *mx, double x, int n)
{
    __m256d mxi, mul, res;
    int i, imax = n / 4;
    int rem = n % 4;

    mul = _mm256_broadcast_sd(&x);

    for (i=0; i<imax; i++) {
	/* multiply 4 doubles in parallel */
	mxi = _mm256_loadu_pd(mx);
	res = _mm256_mul_pd(mul, mxi);
	_mm256_storeu_pd(mx, res);
	mx += 4;
    }

    for (i=0; i<rem; i++) {
	mx[i] *= x;
    }
}
