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

/* complex matrices */

#include "clapack_complex.h"
#include "gretl_cmatrix.h"

#include <fftw3.h>

/* Get two matrices out of array @A, checking that they
   of the same dimensions, and if @square is non-zero
   that they are square.
*/

static int get_two_matrices (gretl_array *A,
			     gretl_matrix **pmr,
			     gretl_matrix **pmi,
			     int square)
{
    gretl_matrix *mr = NULL, *mi = NULL;
    int err = 0;
    
    mr = gretl_array_get_element(A, 0, NULL, &err);
    if (!err) {
	mi = gretl_array_get_element(A, 1, NULL, &err);
    }

    if (!err && (gretl_is_null_matrix(mr) || gretl_is_null_matrix(mi))) {
	err = E_INVARG;
    }

    if (!err) {
	int r = mr->rows;
	int c = mr->cols;

	if (mi->rows != r || mi->cols != c) {
	    err = E_NONCONF;
	}
	if (!err && square && r != c) {
	    err = E_NONCONF;
	}
    }

    if (!err) {
	*pmr = mr;
	*pmi = mi;
    }

    return err;
}

/* Put two matrices (real part, imaginary part) into array @A
   given dimensions @r and @c and complex source @cx.
*/

static int complex_mat_into_array (cmplx *cx, int r, int c,
				   gretl_array *A)
{
    gretl_matrix *mr = gretl_matrix_alloc(r, c);
    gretl_matrix *mi = gretl_matrix_alloc(r, c);
    int i, j, k;
    int err = 0;

    if (mr == NULL || mi == NULL) {
	return E_ALLOC;
    }
	
    for (j=0; j<c; j++) {
	k = j * r;
	for (i=0; i<r; i++) {
	    gretl_matrix_set(mr, i, j, cx[k].r);
	    gretl_matrix_set(mi, i, j, cx[k].i);
	    k++;
	}
    }

    /* note: array takes ownership of these matrices */
    gretl_array_set_matrix(A, 0, mr, 0);
    gretl_array_set_matrix(A, 1, mi, 0);    

    return err;
}

/* Write the content of matrices @mr (real part) and @mi
   (imaginary part) into complex array @cx.
*/

static void matrices_into_complex (const gretl_matrix *mr,
				   const gretl_matrix *mi,
				   cmplx *cx)
{
    int r = mr->rows;
    int c = mr->cols;
    int i, j, k;

    for (j=0; j<c; j++) {
	k = j * r;
	for (i=0; i<r; i++) {
	    cx[k].r = gretl_matrix_get(mr, i, j);
	    cx[k].i = gretl_matrix_get(mi, i, j);
	    k++;
	}
    }
}

/* Compute eigenvalues (and optionally eigenvectors, if @V is non-NULL)
   of a Hermitian matrix, using Lapack.
*/

gretl_matrix *gretl_zheev (gretl_array *A, gretl_array *V, int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *mr, *mi;
    integer n, info, lwork;
    double *w = NULL;
    double *rwork = NULL;
    cmplx *a = NULL;
    cmplx *work = NULL;
    cmplx wsz;
    char jobz = V != NULL ? 'V' : 'N';
    char uplo = 'U';
    int i, j, k;

    *err = get_two_matrices(A, &mr, &mi, 1);

    if (!*err) {
	n = mr->rows;
	ret = gretl_matrix_alloc(n, 1);
	a = malloc(n * n * sizeof *a);
	if (ret == NULL || a == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	goto bailout;
    }

    w = ret->val;

    /* write upper triangle of complex matrix into @a */
    for (j=0; j<n; j++) {
	k = j * n;
	for (i=0; i<n; i++) { /* FIXME triangle */
	    a[k].r = gretl_matrix_get(mr, i, j);
	    a[k].i = gretl_matrix_get(mi, i, j);
	    k++;
	}
    }

    /* get optimal workspace size */
    lwork = -1;
    zheev_(&jobz, &uplo, &n, a, &n, w, &wsz, &lwork, rwork, &info);

    lwork = (integer) wsz.r;
    work = malloc(lwork * sizeof *work);
    rwork = malloc((3 * n - 2) * sizeof *rwork);
    if (work == NULL || rwork == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* do the actual eigen decomposition */
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
	fprintf(stderr, "zheev: info = %d\n", info);
	*err = E_DATA;
    } else if (V != NULL) {
	*err = complex_mat_into_array(a, n, n, V);
    }

 bailout:

    free(rwork);
    free(work);
    free(a);

    if (*err) {
	gretl_matrix_free(ret);
	ret = NULL;
    }

    return ret;
}

/* Compute the inverse of a complex matrix represented by the array
   @A via LU decomposition; uses the Lapack functions zgetrf() and
   zgetri().
*/

gretl_array *gretl_zgetri (gretl_array *A, int *err)
{
    gretl_array *Ainv = NULL;
    gretl_matrix *mr = NULL;
    gretl_matrix *mi = NULL;
    integer lwork = -1;
    integer *ipiv;
    cmplx *a, *work = NULL;
    integer n, info;

    *err = get_two_matrices(A, &mr, &mi, 1); /* square? */
    if (*err) {
	return NULL;
    }

    n = mr->rows;

    a = malloc(n * n *sizeof *a);
    ipiv = malloc(2 * n * sizeof *ipiv);
    if (a == NULL || ipiv == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    matrices_into_complex(mr, mi, a);

    zgetrf_(&n, &n, a, &n, ipiv, &info);
    if (info != 0) {
	printf("zgetrf: info = %d\n", info);
	*err = E_DATA;
    }

    if (!*err) {
	/* workspace size query */
	cmplx wsz;

	zgetri_(&n, a, &n, ipiv, &wsz, &lwork, &info);
	lwork = (integer) wsz.r;
	work = malloc(lwork * sizeof *work);
	if (work == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* actual computation */
	zgetri_(&n, a, &n, ipiv, work, &lwork, &info);
	if (info != 0) {
	    printf("zgetri: info = %d\n", info);
	    *err = E_DATA;
	}
    }

    if (!*err) {
	Ainv = gretl_array_new(GRETL_TYPE_MATRICES, 2, err);
	if (!*err) {
	    *err = complex_mat_into_array(a, n, n, Ainv);
	}
    }

 bailout:
    
    free(work);
    free(ipiv);
    free(a);

    if (*err && Ainv != NULL) {
	gretl_array_destroy(Ainv);
	Ainv = NULL;
    }

    return Ainv;
}

/* Multiplication of complex matrices via Lapack's zgemm(). */

gretl_array *gretl_zgemm (gretl_array *A, gretl_array *B, int *err)
{
    gretl_array *C = NULL;
    gretl_matrix *ar = NULL;
    gretl_matrix *ai = NULL;
    gretl_matrix *br = NULL;
    gretl_matrix *bi = NULL;
    cmplx *a = NULL, *b = NULL;
    cmplx *c = NULL;
    cmplx alpha = {1, 0};
    cmplx beta = {0, 0};
    char transa = 'N';
    char transb = 'N';
    integer m, n, k;

    *err = get_two_matrices(A, &ar, &ai, 0);
    if (*err) {
	return NULL;
    }

    *err = get_two_matrices(B, &br, &bi, 0);
    if (*err) {
	return NULL;
    }

    /* FIXME allow for transposition? */

    m = ar->rows;
    k = ar->cols;
    n = br->cols;

    a = malloc(m * k * sizeof *a);
    b = malloc(k * n * sizeof *b);
    c = malloc(m * n * sizeof *c);

    if (a == NULL || b == NULL || c == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	matrices_into_complex(ar, ai, a);
	matrices_into_complex(br, bi, b);
    }

    if (!*err) {
	integer lda = m;
	integer ldb = k;
	integer ldc = m;

	zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda,
	       b, &ldb, &beta, c, &ldc);
    }

    if (!*err) {
	C = gretl_array_new(GRETL_TYPE_MATRICES, 2, err);
	if (!*err) {
	    *err = complex_mat_into_array(c, m, n, C);
	}
    }

    free(a);
    free(b);
    free(c);

    if (*err && C != NULL) {
	gretl_array_destroy(C);
	C = NULL;
    }

    return C;
}

int complex_matrix_print (gretl_array *A, PRN *prn)
{
    gretl_matrix *mr = NULL;
    gretl_matrix *mi = NULL;
    double re, im;
    char s[4] = "   ";
    int r, c, i, j;
    int err = 0;

    err = get_two_matrices(A, &mr, &mi, 0);
    if (err) {
	return err;
    }

    r = mr->rows;
    c = mr->cols;

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    re = gretl_matrix_get(mr, i, j);
	    im = gretl_matrix_get(mi, i, j);
	    s[1] = (im >= 0) ? '+' : '-';
	    pprintf(prn, "%7.4f%s%6.4fi", re, s, fabs(im));
	    if (j < c - 1) {
		pputs(prn, "  ");
	    }
        }
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    return 0;
}

gretl_array *gretl_complex_fft (gretl_array *A, int inverse, int *err)
{
    gretl_array *B = NULL;
    gretl_matrix *mr, *mi;
    fftw_complex *tmp, *ptr;
    fftw_plan p;
    int sign;
    int r, c, j;

    *err = get_two_matrices(A, &mr, &mi, 0);
    if (*err) {
	return NULL;
    }

    r = mr->rows;
    c = mr->cols;

    tmp = fftw_malloc(r * c * sizeof *tmp);
    if (tmp == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    matrices_into_complex(mr, mi, (cmplx *) tmp);
    sign = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    ptr = tmp;
    for (j=0; j<c; j++) {
	p = fftw_plan_dft_1d(r, ptr, ptr, sign, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	/* advance pointer to next column */
	ptr += r;
    }

    if (inverse) {
	/* "FFTW computes an unnormalized transform: computing a
	    forward followed by a backward transform (or vice versa)
	    will result in the original data multiplied by the size of
	    the transform (the product of the dimensions)."
	*/
	for (j=0; j<r*c; j++) {
	    tmp[j][0] = tmp[j][0] / r;
	    tmp[j][1] = tmp[j][1] / r;
	}
    }

    B = gretl_array_new(GRETL_TYPE_MATRICES, 2, err);
    if (!*err) {
	*err = complex_mat_into_array((cmplx *) tmp, r, c, B);
    }

    fftw_free(tmp);

    return B;
}
