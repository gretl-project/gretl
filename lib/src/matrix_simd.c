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

/* All code here is conditional on either AVX or SSE */

static double *mval_realloc (gretl_matrix *m, int n)
{
    double *newval = mval_malloc(n * sizeof(double));

    if (newval != NULL && m->val != NULL) {
	int old_n = m->rows * m->cols;

	n = n < old_n ? n : old_n;
	memcpy(newval, m->val, n * sizeof(double));
	mval_free(m->val);
    }

    return newval;
}

#ifdef USE_AVX

static void gretl_matrix_simd_add_to (gretl_matrix *a,
				      const gretl_matrix *b,
				      int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int i, imax = n / 4;
    int rem = n % 4;

    for (i=0; i<imax; i++) {
	/* add 4 doubles in parallel */
	__m256d Ymm_A = _mm256_load_pd(ax);
	__m256d Ymm_B = _mm256_load_pd(bx);
	__m256d Ymm_C = _mm256_add_pd(Ymm_A, Ymm_B);
	_mm256_store_pd(ax, Ymm_C);
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	ax[i] += bx[i];
    }
}

static void gretl_matrix_simd_subt_from (gretl_matrix *a,
					 const gretl_matrix *b,
					 int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int i, imax = n / 4;
    int rem = n % 4;

    for (i=0; i<imax; i++) {
	/* subtract 4 doubles in parallel */
	__m256d Ymm_A = _mm256_load_pd(ax);
	__m256d Ymm_B = _mm256_load_pd(bx);
	__m256d Ymm_C = _mm256_sub_pd(Ymm_A, Ymm_B);
	_mm256_store_pd(ax, Ymm_C);
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	ax[i] -= bx[i];
    }
}

#else /* SSE */

static void gretl_matrix_simd_add_to (gretl_matrix *a,
				      const gretl_matrix *b,
				      int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int i, imax = n / 2;

    for (i=0; i<imax; i++) {
	/* add 2 doubles in parallel */
	__m128d Xmm_A = _mm_load_pd(ax);
	__m128d Xmm_B = _mm_load_pd(bx);
	__m128d Xmm_C = _mm_add_pd(Xmm_A, Xmm_B);
	_mm_store_pd(ax, Xmm_C);
	ax += 2;
	bx += 2;
    }

    if (n % 2) {
	ax[0] += bx[0];
    }
}

static void gretl_matrix_simd_subt_from (gretl_matrix *a,
					 const gretl_matrix *b,
					 int n)
{
    double *ax = a->val;
    const double *bx = b->val;
    int i, imax = n / 2;

    for (i=0; i<imax; i++) {
	/* subtract 2 doubles in parallel */
	__m128d Xmm_A = _mm_load_pd(ax);
	__m128d Xmm_B = _mm_load_pd(bx);
	__m128d Xmm_C = _mm_sub_pd(Xmm_A, Xmm_B);
	_mm_store_pd(ax, Xmm_C);
	ax += 2;
	bx += 2;
    }

    if (n % 2) {
	ax[0] -= bx[0];
    }
}

#endif /* AVX vs SSE */

#ifdef USE_AVX

static void gretl_matrix_simd_add (const gretl_matrix *a,
				   const gretl_matrix *b,
				   gretl_matrix *c)
{
    const double *ax = a->val;
    const double *bx = b->val;
    double *cx = c->val;
    int i, n = a->rows * a->cols;
    int imax = n / 4;
    int rem = n % 4;

    for (i=0; i<imax; i++) {
	/* process 4 doubles in parallel */
	__m256d Ymm_A = _mm256_load_pd(ax);
	__m256d Ymm_B = _mm256_load_pd(bx);
	__m256d Ymm_C = _mm256_add_pd(Ymm_A, Ymm_B);
	_mm256_store_pd(cx, Ymm_C);
	cx += 4;
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	cx[i] = ax[i] + bx[i];
    }
}

static void gretl_matrix_simd_subtract (const gretl_matrix *a,
					const gretl_matrix *b,
					gretl_matrix *c)
{
    const double *ax = a->val;
    const double *bx = b->val;
    double *cx = c->val;
    int i, n = a->rows * a->cols;
    int imax = n / 4;
    int rem = n % 4;

    for (i=0; i<imax; i++) {
	/* process 4 doubles in parallel */
	__m256d Ymm_A = _mm256_load_pd(ax);
	__m256d Ymm_B = _mm256_load_pd(bx);
	__m256d Ymm_C = _mm256_sub_pd(Ymm_A, Ymm_B);
	_mm256_store_pd(cx, Ymm_C);
	cx += 4;
	ax += 4;
	bx += 4;
    }

    for (i=0; i<rem; i++) {
	cx[i] = ax[i] - bx[i];
    }
}

#else /* SSE */

static void gretl_matrix_simd_add (const gretl_matrix *a,
				   const gretl_matrix *b,
				   gretl_matrix *c)
{
    const double *ax = a->val;
    const double *bx = b->val;
    double *cx = c->val;
    int i, n = a->rows * a->cols;
    int imax = n / 2;
			
    for (i=0; i<imax; i++) {
	/* process 2 doubles in parallel */
	__m128d Xmm_A = _mm_load_pd(ax);
	__m128d Xmm_B = _mm_load_pd(bx);
	__m128d Xmm_C = _mm_add_pd(Xmm_A, Xmm_B);
	_mm_store_pd(cx, Xmm_C);
	cx += 2;
	ax += 2;
	bx += 2;
    }

    if (n % 2) {
	cx[0] = ax[0] + bx[0];
    }
}

static void gretl_matrix_simd_subtract (const gretl_matrix *a,
					const gretl_matrix *b,
					gretl_matrix *c)
{
    const double *ax = a->val;
    const double *bx = b->val;
    double *cx = c->val;
    int i, n = a->rows * a->cols;
    int imax = n / 2;
			
    for (i=0; i<imax; i++) {
	/* process 2 doubles in parallel */
	__m128d Xmm_A = _mm_load_pd(ax);
	__m128d Xmm_B = _mm_load_pd(bx);
	__m128d Xmm_C = _mm_sub_pd(Xmm_A, Xmm_B);
	_mm_store_pd(cx, Xmm_C);
	cx += 2;
	ax += 2;
	bx += 2;
    }

    if (n % 2) {
	cx[0] = ax[0] - bx[0];
    }
}

#endif /* AVX vs SSE */

#if 0 /* not yet */
#ifdef USE_AVX

/* Note: this is probably usable only for n <= 8 (shortage of AVX
   registers). It is very much faster if m is a multiple of 4; k is
   unconstrained.
*/

static void gretl_matrix_simd_mul (const gretl_matrix *A,
				   const gretl_matrix *B,
				   gretl_matrix *C)
{
    int m = A->rows;
    int n = A->cols;
    int k = B->cols;
    const double *aval = A->val;
    const double *bval = B->val;
    double *cval = C->val;
    int hmax = m / 4;
    int hrem = m % 4;
    int i, j;

    if (m >= 4) {
	__m256d a[n], b[n];
	__m256d mult, ccol;
	int h;

	for (h=0; h<hmax; h++) {
	    /* loop across the n columns of A, loading
	       4 elements from each
	     */
	    if (hrem) {
		/* unaligned */
		for (j=0; j<n; j++) {
		    a[j] = _mm256_loadu_pd(aval + j*m);
		}
	    } else {
		/* aligned */
		for (j=0; j<n; j++) {
		    a[j] = _mm256_load_pd(aval + j*m);
		}
	    }

	    /* loop across the k columns of B */
	    for (j=0; j<k; j++) {
		/* broadcast the n elements of col j of B */
		for (i=0; i<n; i++) {
		    b[i] = _mm256_broadcast_sd(&bval[j*n + i]);
		}

		/* cumulate the products */
		ccol = _mm256_setzero_pd();
		for (i=0; i<n; i++) {
		    mult = _mm256_mul_pd(b[i], a[i]);
		    ccol = _mm256_add_pd(ccol, mult);
		}

		/* write 4 rows of each column of C */
		if (hrem) {
		    /* unaligned */
		    _mm256_storeu_pd(&cval[m*j], ccol);
		} else {
		    /* aligned */
		    _mm256_store_pd(&cval[m*j], ccol);
		}
		    
	    }

	    /* advance by the number of doubles we can store in
	       one go */
	    aval += 4;
	    cval += 4;
	}
    }
    
    if (hrem >= 2) {
	/* do a single 128-bit run */
	__m128d a[n], b[n];
	__m128d mult, ccol;
	int realign = m % 2;

	if (realign) {
	    for (j=0; j<n; j++) {
		a[j] = _mm_loadu_pd(aval + j*m);
	    }
	} else {
	    for (j=0; j<n; j++) {
		a[j] = _mm_load_pd(aval + j*m);
	    }
	}

	for (j=0; j<k; j++) {
	    for (i=0; i<n; i++) {
		b[i] = _mm_set1_pd(bval[j*n + i]);
	    }

	    ccol = _mm_setzero_pd();
	    for (i=0; i<n; i++) {
		mult = _mm_mul_pd(b[i], a[i]);
		ccol = _mm_add_pd(ccol, mult);
	    }

	    if (realign) {
		_mm_storeu_pd(&cval[m*j], ccol);
	    } else {
		_mm_store_pd(&cval[m*j], ccol);
	    }
	}
	hrem -= 2;
	aval += 2;
	cval += 2;
    }
    
    if (hrem) {
	/* odd-valued m: compute the last row */
	double ccol, a[n];

	for (j=0; j<n; j++) {
	    a[j] = aval[j*m];
	}

	for (j=0; j<k; j++) {
	    ccol = 0.0; 
	    for (i=0; i<n; i++) {
		ccol += bval[j*n + i] * a[i];
	    }
	    cval[m*j] = ccol;
	}
    }
}

#else /* SSE */

/* Note: this is probably usable only for n <= 8 (shortage of SSE
   registers). It is faster if m is a multiple of 2; but k is
   unconstrained.
*/

static void gretl_matrix_simd_mul (const gretl_matrix *A,
				   const gretl_matrix *B,
				   gretl_matrix *C)
{
    int m = A->rows;
    int n = A->cols;
    int k = B->cols;
    const double *aval = A->val;
    const double *bval = B->val;
    double *cval = C->val;
    int hmax = m / 2;
    int hrem = m % 2;
    int i, j;

    if (m >= 2) {
	__m128d a[n], b[n];
	__m128d mult, ccol;
	int h;

	for (h=0; h<hmax; h++) {
	    /* loop across the n columns of A, loading
	       2 elements from each
	    */
	    if (hrem) {
		for (j=0; j<n; j++) {
		    a[j] = _mm_loadu_pd(aval + j*m);
		}
	    } else {
		for (j=0; j<n; j++) {
		    a[j] = _mm_load_pd(aval + j*m);
		}
	    }

	    /* loop across the k columns of B */
	    for (j=0; j<k; j++) {
		/* broadcast the n elements of col j of B */
		for (i=0; i<n; i++) {
		    b[i] = _mm_set1_pd(bval[j*n + i]);
		}

		/* cumulate the products */
		ccol = _mm_setzero_pd();
		for (i=0; i<n; i++) {
		    mult = _mm_mul_pd(b[i], a[i]);
		    ccol = _mm_add_pd(ccol, mult);
		}

		/* write 2 rows of each column of C */
		if (hrem) {
		    _mm_storeu_pd(&cval[m*j], ccol);
		} else {
		    _mm_store_pd(&cval[m*j], ccol);
		}
	    }

	    /* advance by the number of doubles we can store in
	       one go */
	    aval += 2;
	    cval += 2;
	}
    }

    if (hrem) {
	/* odd-valued m: compute the last row */
	double ccol, a[n];

	for (j=0; j<n; j++) {
	    a[j] = aval[j*m];
	}

	for (j=0; j<k; j++) {
	    ccol = 0.0; 
	    for (i=0; i<n; i++) {
		ccol += bval[j*n + i] * a[i];
	    }
	    cval[m*j] = ccol;
	}
    }
}

#endif /* AVX vs SIMD */
#endif /* not yet */

