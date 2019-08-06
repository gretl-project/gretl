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
#include "clapack_complex.h"
#include "gretl_cmatrix.h"

/* Note: since we include gretl_cmatrix.h (which in turn includes
   C99's complex.h) before fftw3.h, FFTW's fftw_complex will be
   typedef'd to C99's "double complex" and can be manipulated as
   such.
*/

#include <fftw3.h>
#include <errno.h>

#define gretl_cmatrix_get(z,r,i,j) (z[(j)*r+(i)])
#define gretl_cmatrix_set(z,r,i,j,v) (z[(j)*r+(i)]=v)

#define cscalar(m) (m->rows == 2 && m->cols == 1)

/* helper function for fftw-based real FFT functions */

static int fft_allocate (double **px, gretl_matrix **pm,
			 fftw_complex **pc, int r, int c)
{
    *pm = gretl_matrix_alloc(r, c);
    if (*pm == NULL) {
	return E_ALLOC;
    }

    *px = fftw_malloc(r * sizeof **px);
    if (*px == NULL) {
	gretl_matrix_free(*pm);
	return E_ALLOC;
    }

    *pc = fftw_malloc((r/2 + 1 + r % 2) * sizeof **pc);
    if (*pc == NULL) {
	gretl_matrix_free(*pm);
	free(*px);
	return E_ALLOC;
    }

    return 0;
}

/* start fftw-based real FFT functions */

static gretl_matrix *
real_gretl_matrix_fft (const gretl_matrix *y, int inverse, int *err)
{
    gretl_matrix *ft = NULL;
    fftw_plan p = NULL;
    double *ffx = NULL;
    double xr, xi;
    fftw_complex *ffz;
    int r, c, m, odd, cr, ci;
    int i, j, cdim;

    if (y->rows < 2) {
	*err = E_DATA;
	return NULL;
    }

    r = y->rows;
    m = r / 2;
    odd = r % 2;
    c = inverse ? y->cols / 2 : y->cols;

    if (c == 0) {
	*err = E_NONCONF;
	return NULL;
    }

    cdim = inverse ? c : 2 * c;
    *err = fft_allocate(&ffx, &ft, &ffz, r, cdim);
    if (*err) {
	return NULL;
    }

    cr = 0;
    ci = 1;
    for (j=0; j<c; j++) {
	/* load the data */
	if (inverse) {
	    for (i=0; i<=m+odd; i++) {
		xr = gretl_matrix_get(y, i, cr);
		xi = gretl_matrix_get(y, i, ci);
		ffz[i] = xr + xi * I;
	    }
	} else {
	    for (i=0; i<r; i++) {
		ffx[i] = gretl_matrix_get(y, i, j);
	    }
	}

	if (j == 0) {
	    /* make the plan just once */
	    if (inverse) {
		p = fftw_plan_dft_c2r_1d(r, ffz, ffx, FFTW_ESTIMATE);
	    } else {
		p = fftw_plan_dft_r2c_1d(r, ffx, ffz, FFTW_ESTIMATE);
	    }
	}

	/* run the transform */
	fftw_execute(p);

	/* transcribe the result */
	if (inverse) {
	    for (i=0; i<r; i++) {
		gretl_matrix_set(ft, i, j, ffx[i] / r);
	    }
	} else {
	    for (i=0; i<=m+odd; i++) {
		gretl_matrix_set(ft, i, cr, creal(ffz[i]));
		gretl_matrix_set(ft, i, ci, cimag(ffz[i]));
	    }
	    for (i=m; i>0; i--) {
		gretl_matrix_set(ft, r-i, cr,  creal(ffz[i]));
		gretl_matrix_set(ft, r-i, ci, -cimag(ffz[i]));
	    }
	}
	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(ffz);
    fftw_free(ffx);

    return ft;
}

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err)
{
    return real_gretl_matrix_fft(y, 0, err);
}

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err)
{
    return real_gretl_matrix_fft(y, 1, err);
}

/* end fftw-based real FFT functions */

static int cmatrix_validate (const gretl_matrix *m, int square)
{
    if (gretl_is_null_matrix(m)) {
	return 0;
    } else if (m->rows % 2 != 0) {
	return 0;
    } else if (square && m->rows != 2 * m->cols) {
	return 0;
    } else {
	return 1;
    }
}

/* Construct a complex matrix, with all-zero imaginary part,
   from real matrix @A.
*/

static gretl_matrix *complex_from_real (const gretl_matrix *A,
					int *err)
{
    gretl_matrix *C = NULL;

    if (gretl_is_null_matrix(A)) {
	*err = E_DATA;
	return NULL;
    }

    C = gretl_zero_matrix_new(A->rows * 2, A->cols);

    if (C == NULL) {
	*err = E_ALLOC;
    } else {
	double complex *cz = (double complex *) C->val;
	int i, n = A->rows * A->cols;

	for (i=0; i<n; i++) {
	    cz[i] = A->val[i];
	}
	C->is_complex = 1;
    }

    return C;
}

/* Multiplication of complex matrices via BLAS zgemm().
   We allow for the possibility that the conjugate
   transpose of @A should be used. This could be
   generalized to allow [conjugate] transposition of
   B too.
*/

static gretl_matrix *gretl_zgemm (const gretl_matrix *A,
				  GretlMatrixMod amod,
				  const gretl_matrix *B,
				  int *err)
{
    gretl_matrix *C;
    cmplx *a, *b, *c;
    cmplx alpha = {1, 0};
    cmplx beta = {0, 0};
    char transa = 'N';
    char transb = 'N';
    integer lda, ldb, ldc;
    integer m, n, k;

    if (!cmatrix_validate(A, 0) || !cmatrix_validate(B, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    if (amod == GRETL_MOD_CTRANSP) {
	transa = 'C';
	m = A->cols;     /* rows of op(A) */
	k = A->rows / 2; /* cols of op(A) */
    } else {
	m = A->rows / 2; /* complex rows of A */
	k = A->cols;     /* cols of A */
    }
    n = B->cols;

    if (k != B->rows / 2) {
	*err = E_NONCONF;
	return NULL;
    }

    C = gretl_matrix_alloc(m*2, n);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    a = (cmplx *) A->val;
    b = (cmplx *) B->val;
    c = (cmplx *) C->val;

    lda = A->rows / 2;
    ldb = B->rows / 2;
    ldc = C->rows / 2;

    zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda,
	   b, &ldb, &beta, c, &ldc);

    C->is_complex = 1;

    return C;
}

/* Multiplication of @A times @B, where we know that at
   least one of them is complex; the other, if it's not
   complex, must be converted to a complex matrix with
   a zero imaginary part first.
*/

gretl_matrix *gretl_cmatrix_multiply (const gretl_matrix *A,
				      const gretl_matrix *B,
				      int *err)
{
    gretl_matrix *T = NULL;
    gretl_matrix *C = NULL;

    if (A->is_complex && B->is_complex) {
	if (cscalar(A) || cscalar(B)) {
	    return gretl_cmatrix_dot_op(A, B, '*', err);
	} else {
	    C = gretl_zgemm(A, 0, B, err);
	}
    } else if (A->is_complex) {
	/* case of real B */
	T = complex_from_real(B, err);
	if (T != NULL) {
	    C = gretl_zgemm(A, 0, T, err);
	}
    } else if (B->is_complex) {
	/* case of real A */
	T = complex_from_real(A, err);
	if (T != NULL) {
	    C = gretl_zgemm(T, 0, B, err);
	}
    } else {
	*err = E_TYPES;
    }

    gretl_matrix_free(T);

    return C;
}

/* Returns (conjugate transpose of A, or A^H) times B,
   allowing for the possibility that either A or B (but
   not both!) may be a real matrix on input.
*/

gretl_matrix *gretl_cmatrix_AHB (const gretl_matrix *A,
				 const gretl_matrix *B,
				 int *err)
{
    gretl_matrix *C = NULL;

    if (A->is_complex && B->is_complex) {
	C = gretl_zgemm(A, GRETL_MOD_CTRANSP, B, err);
    } else if (A->is_complex) {
	/* case of real B */
	gretl_matrix *CB = complex_from_real(B, err);

	if (CB != NULL) {
	    C = gretl_zgemm(A, GRETL_MOD_CTRANSP, CB, err);
	    gretl_matrix_free(CB);
	}
    } else {
	/* case of real A */
	gretl_matrix *CA = complex_from_real(A, err);

	if (CA != NULL) {
	    C = gretl_zgemm(CA, GRETL_MOD_CTRANSP, B, err);
	    gretl_matrix_free(CA);
	}
    }

    return C;
}

/* Eigen decomposition of complex (Hermitian) matrix using
   LAPACK's zheev() */

gretl_matrix *
gretl_zheev (const gretl_matrix *A, int eigenvecs, int *err)
{
    gretl_matrix *evals = NULL;
    integer n, info, lwork;
    double *w = NULL;
    double *rwork = NULL;
    cmplx *a = NULL;
    cmplx *work = NULL;
    cmplx wsz;
    char jobz = eigenvecs ? 'V' : 'N';
    char uplo = 'U';

    if (!cmatrix_validate(A, 1)) {
	*err = E_INVARG;
	return NULL;
    }

    n = A->rows / 2;

    evals = gretl_matrix_alloc(n, 1);
    if (evals == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    w = evals->val;
    a = (cmplx *) A->val;

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

    /* do the actual decomposition */
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
	fprintf(stderr, "zheev: info = %d\n", info);
	*err = E_DATA;
    }

 bailout:

    free(rwork);
    free(work);

    if (*err) {
	gretl_matrix_free(evals);
	evals = NULL;
    }

    return evals;
}

/* Eigen decomposition of complex (non-Hermitian) matrix using
   LAPACK's zgeev() */

gretl_matrix *gretl_zgeev (const gretl_matrix *A,
			   gretl_matrix *VL,
			   gretl_matrix *VR,
			   int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *Acpy = NULL;
    gretl_matrix *Ltmp = NULL;
    gretl_matrix *Rtmp = NULL;
    integer n, info, lwork;
    integer ldvl, ldvr;
    double *w = NULL;
    double *rwork = NULL;
    cmplx *a = NULL;
    cmplx *work = NULL;
    cmplx *vl = NULL, *vr = NULL;
    cmplx wsz;
    char jobvl = VL != NULL ? 'V' : 'N';
    char jobvr = VR != NULL ? 'V' : 'N';

    if (!cmatrix_validate(A, 1)) {
	*err = E_INVARG;
	return NULL;
    }

    n = A->rows / 2;
    ldvl = VL != NULL ? n : 1;
    ldvr = VR != NULL ? n : 1;

    /* we need a copy of @A, which gets overwritten */
    Acpy = gretl_matrix_copy(A);
    if (Acpy == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    a = (cmplx *) Acpy->val;

    if (VL != NULL) {
	/* left eigenvectors wanted */
	if (VL->rows * VL->cols == A->rows * A->cols) {
	    /* VL is useable as is */
	    VL->rows = A->rows;
	    VL->cols = A->cols;
	    vl = (cmplx *) VL->val;
	} else {
	    /* we need to allocate storage */
	    Ltmp = gretl_zero_matrix_new(A->rows, A->cols);
	    if (Ltmp == NULL) {
		*err = E_ALLOC;
		goto bailout;
	    }
	    vl = (cmplx *) Ltmp->val;
	}
    }

    if (VR != NULL) {
	/* right eigenvectors wanted */
	if (VR->rows * VR->cols == A->rows * A->cols) {
	    /* VR is useable as is */
	    VR->rows = A->rows;
	    VR->cols = A->cols;
	    vr = (cmplx *) VR->val;
	} else {
	    /* we need to allocate storage */
	    Rtmp = gretl_zero_matrix_new(A->rows, A->cols);
	    if (Rtmp == NULL) {
		*err = E_ALLOC;
		goto bailout;
	    }
	    vr = (cmplx *) Rtmp->val;
	}
    }

    rwork = malloc(2 * n * sizeof *rwork);
    ret = gretl_matrix_alloc(2 * n, 1);
    if (rwork == NULL || ret == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    w = ret->val;
    ret->is_complex = 1;

    /* get optimal workspace size */
    lwork = -1;
    zgeev_(&jobvl, &jobvr, &n, a, &n, w, vl, &ldvl, vr, &ldvr,
	   &wsz, &lwork, rwork, &info);
    lwork = (integer) wsz.r;
    work = malloc(lwork * sizeof *work);
    if (work == NULL || rwork == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* do the actual decomposition */
    zgeev_(&jobvl, &jobvr, &n, a, &n, w, vl, &ldvl, vr, &ldvr,
	   work, &lwork, rwork, &info);
    if (info != 0) {
	fprintf(stderr, "zgeev: info = %d\n", info);
	*err = E_DATA;
    } else {
	if (Ltmp != NULL) {
	    gretl_matrix_replace_content(VL, Ltmp);
	    VL->is_complex = 1;
	}
	if (Rtmp != NULL) {
	    gretl_matrix_replace_content(VR, Rtmp);
	    VR->is_complex = 1;
	}
    }

 bailout:

    free(rwork);
    free(work);
    gretl_matrix_free(Acpy);
    gretl_matrix_free(Ltmp);
    gretl_matrix_free(Rtmp);

    if (*err) {
	gretl_matrix_free(ret);
	ret = NULL;
    }

    return ret;
}

/* Inverse of a complex matrix via LU decomposition using the
   LAPACK functions zgetrf() and zgetri()
*/

gretl_matrix *gretl_zgetri (const gretl_matrix *A, int *err)
{
    gretl_matrix *Ainv = NULL;
    integer lwork = -1;
    integer *ipiv;
    cmplx *a, *work = NULL;
    integer n, info;

    if (!cmatrix_validate(A, 1)) {
	*err = E_INVARG;
	return NULL;
    }

    Ainv = gretl_matrix_copy(A);
    if (Ainv == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    n = A->cols;

    ipiv = malloc(2 * n * sizeof *ipiv);
    if (ipiv == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    a = (cmplx *) Ainv->val;

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
	/* actual inversion */
	zgetri_(&n, a, &n, ipiv, work, &lwork, &info);
	if (info != 0) {
	    printf("zgetri: info = %d\n", info);
	    *err = E_DATA;
	}
    }

 bailout:

    free(work);
    free(ipiv);

    if (*err && Ainv != NULL) {
	gretl_matrix_free(Ainv);
	Ainv = NULL;
    }

    return Ainv;
}

enum {
    SVD_THIN,
    SVD_FULL
};

/* SVD of a complex matrix via the LAPACK function zgesvd() */

int gretl_cmatrix_SVD (const gretl_matrix *x, gretl_matrix **pu,
		       gretl_vector **ps, gretl_matrix **pvt,
		       int smod)
{
    integer m, n, lda;
    integer ldu = 1, ldvt = 1;
    integer lwork = -1L;
    double *rwork = NULL;
    integer info;
    gretl_matrix *a = NULL;
    gretl_matrix *s = NULL;
    gretl_matrix *u = NULL;
    gretl_matrix *vt = NULL;
    cmplx *az;
    char jobu = 'N', jobvt = 'N';
    cmplx zu, zvt;
    cmplx *uval = &zu;
    cmplx *vtval = &zvt;
    cmplx *work = NULL;
    int xsize, k, err = 0;

    if (pu == NULL && ps == NULL && pvt == NULL) {
	/* no-op */
	return 0;
    }

    if (gretl_is_null_matrix(x)) {
	return E_DATA;
    }

    lda = m = x->rows / 2;
    n = x->cols;
    xsize = lda * n;

    if (smod == SVD_THIN && m < n) {
	fprintf(stderr, "real_gretl_matrix_SVD: a is %d x %d, should be 'thin'\n",
		a->rows, a->cols);
	return E_NONCONF;
    }

    az = malloc(xsize * sizeof *az);
    if (az == NULL) {
	return E_ALLOC;
    }
    memcpy(az, x->val, xsize * sizeof *az);

    k = (m < n)? m : n;

    s = gretl_vector_alloc(k);
    if (s == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (pu != NULL) {
	ldu = m;
	if (smod == SVD_FULL) {
	    u = gretl_matrix_alloc(2*ldu, m);
	} else {
	    u = gretl_matrix_alloc(2*ldu, n);
	}
	if (u == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
	    u->is_complex = 1;
	    uval = (cmplx *) u->val;
	    jobu = (smod == SVD_FULL)? 'A' : 'S';
	}
    }

    if (pvt != NULL) {
	ldvt = n;
	vt = gretl_matrix_alloc(2*ldvt, n);
	if (vt == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
	    vt->is_complex = 1;
	    vtval = (cmplx *) vt->val;
	    jobvt = 'A';
	}
    }

    work = malloc(sizeof *work);
    rwork = malloc((5 * MIN(m,n)) * sizeof *rwork);
    if (work == NULL || rwork == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* workspace query */
    lwork = -1;
    zgesvd_(&jobu, &jobvt, &m, &n, az, &lda, s->val, uval, &ldu,
	    vtval, &ldvt, work, &lwork, rwork, &info);

    if (info != 0 || work[0].r <= 0.0) {
	fprintf(stderr, "zgesvd: workspace query failed\n");
	err = E_DATA;
	goto bailout;
    }

    lwork = (integer) work[0].r;
    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* actual computation */
    zgesvd_(&jobu, &jobvt, &m, &n, az, &lda, s->val, uval, &ldu,
	    vtval, &ldvt, work, &lwork, rwork, &info);

    if (info != 0) {
	fprintf(stderr, "gretl_cmatrix_SVD: info = %d\n", (int) info);
	err = E_DATA;
	goto bailout;
    }

    if (ps != NULL) {
	*ps = s;
	s = NULL;
    }
    if (pu != NULL) {
	*pu = u;
	u = NULL;
    }
    if (pvt != NULL) {
	*pvt = vt;
	vt = NULL;
    }

 bailout:

    free(az);
    free(work);
    free(rwork);
    gretl_matrix_free(s);
    gretl_matrix_free(u);
    gretl_matrix_free(vt);

    return err;
}

/* Complex FFT (or inverse) via fftw */

gretl_matrix *gretl_complex_fft (const gretl_matrix *A,
				 int inverse, int *err)
{
    gretl_matrix *B = NULL;
    fftw_complex *tmp, *ptr;
    fftw_plan p;
    int sign;
    int r, c, j;

    if (!cmatrix_validate(A, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    B = gretl_matrix_copy(A);
    if (B == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    r = A->rows / 2;
    c = A->cols;

    tmp = (fftw_complex *) B->val;
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
	    So should we do the following?
	*/
	for (j=0; j<r*c; j++) {
	    tmp[j] /= r;
	}
    }

    return B;
}

#define cmatrix_get_re(m,i,j) (m->val[(j)*m->rows+(i)*2])
#define cmatrix_get_im(m,i,j) (m->val[(j)*m->rows+(i)*2+1])

int complex_matrix_print (const gretl_matrix *A,
			  const char *name,
			  PRN *prn)
{
    double complex *az;
    double complex aij;
    double re, im, xmax;
    char s[4] = "   ";
    int all_ints = 1;
    int zwidth = 0;
    int r, c, i, j;

    if (!cmatrix_validate(A, 0)) {
	return E_INVARG;
    }

    r = A->rows / 2;
    c = A->cols;

    az = (double complex *) A->val;
    xmax = 0;

    for (j=0; j<c && all_ints; j++) {
	for (i=0; i<r && all_ints; i++) {
	    aij = gretl_cmatrix_get(az, r, i, j);
	    re = creal(aij);
	    im = cimag(aij);
	    if (floor(re) != re || floor(im) != im) {
		all_ints = 0;
	    } else {
		re = MAX(fabs(re), fabs(im));
		if (re > xmax) {
		    xmax = re;
		}
	    }
	}
    }

    if (all_ints && xmax > 0) {
	/* try for a more compact format */
	xmax = log10(xmax);
	if (xmax > 0 && xmax < 3) {
	    zwidth = floor(xmax) + 2;
	}
    }

    if (name != NULL && *name != '\0') {
	pprintf(prn, "%s (%d x %d)", name, r, c);
	pputs(prn, "\n\n");
    }

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    re = cmatrix_get_re(A, i, j);
	    im = cmatrix_get_im(A, i, j);
	    s[1] = (im >= 0) ? '+' : '-';
	    if (zwidth > 0) {
		pprintf(prn, "%*g%s%*gi", zwidth, re, s, zwidth-1, fabs(im));
	    } else {
		pprintf(prn, "%7.4f%s%6.4fi", re, s, fabs(im));
	    }
	    if (j < c - 1) {
		pputs(prn, "  ");
	    }
        }
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    return 0;
}

int complex_matrix_printf (const gretl_matrix *A,
			   const char *fmt,
			   PRN *prn)
{
    double re, im;
    gchar *fmtstr = NULL;
    char s[3] = "  ";
    int r, c, i, j;

    if (!cmatrix_validate(A, 0)) {
	return E_INVARG;
    }

    if (fmt == NULL) {
	fmt = "%7.4f";
    } else {
	/* we only accept floating-point formats */
	char c = fmt[strlen(fmt) - 1];

	if (c != 'f' && c != 'g') {
	    return E_INVARG;
	}
    }

    fmtstr = g_strdup_printf("%s%%s%si", fmt, fmt);

    r = A->rows / 2;
    c = A->cols;

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    re = cmatrix_get_re(A, i, j);
	    im = cmatrix_get_im(A, i, j);
	    s[1] = (im >= 0) ? '+' : '-';
	    pprintf(prn, fmtstr, re, s, fabs(im));
	    if (j < c - 1) {
		pputs(prn, "  ");
	    }
        }
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    g_free(fmtstr);

    return 0;
}

/* Compose a complex matrix from its real and imaginary
   components. If @Im is NULL the matrix will have a
   constant imaginary part given by @ival; otherwise
   the matrices @Re and @Im must be of the same
   dimensions.
*/

gretl_matrix *gretl_cmatrix (const gretl_matrix *Re,
			     const gretl_matrix *Im,
			     double ival, int *err)
{
    gretl_matrix *C = NULL;
    int i, n;

    if (gretl_is_null_matrix(Re)) {
	*err = E_INVARG;
    } else if (Im != NULL) {
	if (Im->rows != Re->rows || Im->cols != Re->cols) {
	    *err = E_NONCONF;
	}
    }

    if (!*err) {
	n = Re->rows * Re->cols;
	C = gretl_matrix_alloc(2 * Re->rows, Re->cols);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	double complex *cz = (double complex *) C->val;

	for (i=0; i<n; i++) {
	    if (Im != NULL) {
		ival = Im->val[i];
	    }
	    cz[i] = Re->val[i] + ival * I;
	}
	C->is_complex = 1;
    }

    return C;
}

/* Extract the real part of @A if @im = 0 or the imaginary part
   if @im is non-zero.
*/

gretl_matrix *gretl_cxtract (const gretl_matrix *A, int im,
			     int *err)
{
    gretl_matrix *B = NULL;
    int i, r, n;

    if (!cmatrix_validate(A, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    r = A->rows / 2;
    n = r * A->cols;

    B = gretl_matrix_alloc(r, A->cols);

    if (B == NULL) {
	*err = E_ALLOC;
    } else {
	double complex *az = (double complex *) A->val;

	for (i=0; i<n; i++) {
	    B->val[i] = im ? cimag(az[i]) : creal(az[i]);
	}
    }

    return B;
}

/* Return [conjugate] transpose of complex matrix @A */

gretl_matrix *gretl_ctrans (const gretl_matrix *A,
			    int conjugate, int *err)
{
    gretl_matrix *C = NULL;
    double complex *az, *cz;
    int i, j, ra, rc;

    if (!cmatrix_validate(A, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    C = gretl_matrix_alloc(2 * A->cols, A->rows / 2);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    az = (double complex *) A->val;
    cz = (double complex *) C->val;
    ra = A->rows / 2;
    rc = C->rows / 2;

    for (j=0; j<A->cols; j++) {
	for (i=0; i<ra; i++) {
	    cz[i*rc+j] = conjugate ? conj(az[j*ra+i]) : az[j*ra+i];
	}
    }

    C->is_complex = 1;

    return C;
}

/* Convert complex matrix @A to its conjugate transpose */

int gretl_ctrans_in_place (gretl_matrix *A)
{
    gretl_matrix *C = NULL;
    double complex *a, *c;
    int i, j, ra, rc;
    int err = 0;

    if (!cmatrix_validate(A, 0)) {
	return E_INVARG;
    }

    C = gretl_matrix_alloc(2 * A->cols, A->rows / 2);
    if (C == NULL) {
	return E_ALLOC;
    }

    a = (double complex *) A->val;
    c = (double complex *) C->val;
    ra = A->rows / 2;
    rc = C->rows / 2;

    for (j=0; j<A->cols; j++) {
	for (i=0; i<ra; i++) {
	    c[i*rc+j] = conj(a[j*ra+i]);
	}
    }

    A->rows = C->rows;
    A->cols = C->cols;

    memcpy(A->val, C->val, C->rows * C->cols * sizeof(double));
    gretl_matrix_destroy_info(A);

    gretl_matrix_free(C);

    return err;
}

/* Addition or subtraction of matrices: handle the case
   where one operand is complex and the other is real.
   In addition, if both matrices are complex, handle the
   case where one of them is a complex scalar.
*/

gretl_matrix *cmatrix_add_sub (const gretl_matrix *A,
			       const gretl_matrix *B,
			       int sgn, int *err)
{
    gretl_matrix *C = NULL;
    int cr = A->rows;
    int cc = A->cols;
    int a_scalar = 0;
    int b_scalar = 0;

    if (A->is_complex && cscalar(B)) {
	b_scalar = 1;
	goto allocate;
    } else if (cscalar(A) && B->is_complex) {
	cr = B->rows;
	cc = B->cols;
	a_scalar = 1;
	goto allocate;
    }

    if (B->cols != A->cols) {
	*err = E_NONCONF;
	return NULL;
    }

    if (A->is_complex && B->is_complex) {
	/* both complex */
	if (B->rows != cr) {
	    *err = E_NONCONF;
	}
    } else if (A->is_complex) {
	/* A complex, B real */
	if (B->rows != cr/2) {
	    *err = E_NONCONF;
	}
    } else {
	/* A real, B complex */
	cr = B->rows;
	if (A->rows != cr/2) {
	    *err = E_NONCONF;
	}
    }

 allocate:

    if (!*err) {
	C = gretl_matrix_alloc(cr, cc);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	double complex *cz = (double complex *) C->val;
	double complex *az = (double complex *) A->val;
	double complex *bz = (double complex *) B->val;
	int i, n = cc * cr / 2;

	if (b_scalar) {
	    for (i=0; i<n; i++) {
		cz[i] = sgn < 0 ? az[i] - bz[0] : az[i] + bz[0];
	    }
	} else if (a_scalar) {
	    for (i=0; i<n; i++) {
		cz[i] = sgn < 0 ? az[0] - bz[i] : az[0] + bz[i];
	    }
	} else if (A->is_complex && B->is_complex) {
	    for (i=0; i<n; i++) {
		cz[i] = sgn < 0 ? az[i] - bz[i] : az[i] + bz[i];
	    }
	} else if (A->is_complex) {
	    for (i=0; i<n; i++) {
		cz[i] = sgn < 0 ? az[i] - B->val[i] : az[i] + B->val[i];
	    }
	} else {
	    for (i=0; i<n; i++) {
		cz[i] = sgn < 0 ? A->val[i] - bz[i] : A->val[i] + bz[i];
	    }
	}
	C->is_complex = 1;
    }

    return C;
}

/* When adding a real scalar to a complex matrix, only
   the real elements of the matrix should be changed.
*/

int cmatrix_add_scalar (gretl_matrix *targ,
			const gretl_matrix *A,
			double x, int Asign)
{
    double complex *tz = (double complex *) targ->val;
    double complex *az = (double complex *) A->val;
    int i, n = (A->rows / 2) * A->cols;

    for (i=0; i<n; i++) {
	tz[i] = Asign < 0 ? x - az[i] : x + az[i];
    }

    return 0;
}

int apply_cmatrix_dfunc (gretl_matrix *targ,
			const gretl_matrix *src,
			double (*dfunc) (double complex))
{
    double complex *csrc = (double complex *) src->val;
    int n = src->cols * src->rows / 2;
    int i, err = 0;

    errno = 0;

    if (dfunc != NULL) {
	for (i=0; i<n; i++) {
	    targ->val[i] = dfunc(csrc[i]);
	}
    }

    if (errno) {
	gretl_errmsg_set_from_errno(NULL, errno);
	err = E_DATA;
    }

    return err;
}

int apply_cmatrix_cfunc (gretl_matrix *targ,
			const gretl_matrix *src,
			double complex (*cfunc) (double complex))
{
    double complex *csrc = (double complex *) src->val;
    double complex *ctarg = (double complex *) targ->val;
    int n = src->cols * src->rows / 2;
    int i, err = 0;

    errno = 0;

    targ->is_complex = 1;
    for (i=0; i<n; i++) {
	ctarg[i] = cfunc(csrc[i]);
    }

    if (errno) {
	gretl_errmsg_set_from_errno(NULL, errno);
	err = E_DATA;
    }

    return err;
}

int apply_cmatrix_unary_op (gretl_matrix *targ,
			    const gretl_matrix *src,
			    int op)
{
    double complex *csrc = (double complex *) src->val;
    double complex *ctarg = (double complex *) targ->val;
    int i, n = src->cols * src->rows / 2;
    int err = 0;

    targ->is_complex = 1;
    for (i=0; i<n && !err; i++) {
	if (op == 1) {
	    /* U_NEG */
	    ctarg[i] = -csrc[i];
	} else if (op == 2) {
	    /* U_POS */
	    ctarg[i] = csrc[i];
	} else if (op == 3) {
	    /* U_NOT */
	    ctarg[i] = (csrc[i] == 0);
	} else {
	    err = E_INVARG;
	}
    }

    return err;
}

static gretl_matrix *complex_scalar_to_mat (double complex z,
					    int *err)
{
    gretl_matrix *ret = gretl_matrix_alloc(2, 1);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	ret->val[0] = creal(z);
	ret->val[1] = cimag(z);
	ret->is_complex = 1;
    }

    return ret;
}

gretl_matrix *gretl_cmatrix_determinant (const gretl_matrix *X,
					 int *err)
{
    gretl_matrix *E = NULL;
    gretl_matrix *ret = NULL;

    if (!cmatrix_validate(X, 1)) {
	*err = E_INVARG;
	return ret;
    }

    E = gretl_zgeev(X, NULL, NULL, err);
    if (E != NULL) {
	E->is_complex = 1;
    }

    if (E != NULL) {
	double complex cret = 1;
	double complex csl;
	int i, n = E->rows / 2;
	int k = 0;

	for (i=0; i<n; i++) {
	    csl = E->val[k] + E->val[k+1] * I;
	    cret *= csl;
	    k += 2;
	}
	gretl_matrix_free(E);
	ret = complex_scalar_to_mat(cret, err);
    }

    return ret;
}

gretl_matrix *gretl_cmatrix_trace (const gretl_matrix *X,
				   int *err)
{
    gretl_matrix *ret = NULL;

    if (!cmatrix_validate(X, 1)) {
	*err = E_INVARG;
    } else {
	double complex *zx = (double complex *) X->val;
	double complex tr = 0;
	int i, r = X->rows / 2;

	for (i=0; i<r; i++) {
	    tr += gretl_cmatrix_get(zx, r, i, i);
	}
	ret = complex_scalar_to_mat(tr, err);
    }

    return ret;
}

/* Retrieve the diagonal of complex matrix @X in the form
   of a complex column vector.
*/

gretl_matrix *gretl_cmatrix_get_diagonal (const gretl_matrix *X,
					  int *err)
{
    gretl_matrix *ret = NULL;

    if (!cmatrix_validate(X, 0)) {
	*err = E_INVARG;
    } else {
	int r = X->rows / 2;
	int d = MIN(r, X->cols);

	ret = gretl_matrix_alloc(2*d, 1);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    double complex *zx = (double complex *) X->val;
	    double complex *zd = (double complex *) ret->val;
	    int i;

	    for (i=0; i<d; i++) {
		zd[i] = gretl_cmatrix_get(zx, r, i, i);
	    }
	    ret->is_complex = 1;
	}
    }

    return ret;
}

/* Set the diagonal of complex matrix @targ using either
   @src (if not NULL) or @x. In the first case @src can
   be either a complex vector of the right length, or a
   real vector, or a complex scalar.
*/

int gretl_cmatrix_set_diagonal (gretl_matrix *targ,
				const gretl_matrix *src,
				double x)
{
    double complex *zsrc = NULL;
    double complex *ztarg = NULL;
    double complex zi = 0;
    int r, d, i;
    int match = 0;
    int err = 0;

    if (!cmatrix_validate(targ, 0)) {
	return E_INVARG;
    }

    r = targ->rows / 2;
    d = MIN(r, targ->cols);

    if (src != NULL) {
	if (src->is_complex) {
	    if ((src->rows == 2 && src->cols == d) ||
		(src->cols == 1 && src->rows == 2*d)) {
		/* conformable complex vector */
		zsrc = (double complex *) src->val;
		match = 1;
	    } else if (src->rows == 2 && src->cols == 1) {
		/* complex scalar */
		zi = *(double complex *) src->val;
		match = 2;
	    }
	} else if (gretl_vector_get_length(src) == d) {
	    /* conformable real vector */
	    match = 3;
	}
    } else {
	/* use real scalar, @x */
	zi = x;
	match = 4;
    }

    if (match == 0) {
	return E_NONCONF;
    }

    ztarg = (double complex *) targ->val;

    for (i=0; i<d; i++) {
	if (match == 1) {
	    gretl_cmatrix_set(ztarg, r, i, i, zsrc[i]);
	} else if (match == 3) {
	    gretl_cmatrix_set(ztarg, r, i, i, src->val[i]);
	} else {
	    gretl_cmatrix_set(ztarg, r, i, i, zi);
	}
    }

    return err;
}

gretl_matrix *gretl_cmatrix_vech (const gretl_matrix *X,
				  int *err)
{
    gretl_matrix *ret = NULL;

    if (!cmatrix_validate(X, 1)) {
	*err = E_INVARG;
    } else {
	int r = X->rows / 2;
	int m = r * (r+1) / 2;

	ret = gretl_matrix_alloc(2*m, 1);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    double complex *zx = (double complex *) X->val;
	    double complex *zv = (double complex *) ret->val;
	    int i, j;

	    m = 0;
	    for (i=0; i<r; i++) {
		for (j=i; j<r; j++) {
		    zv[m++] = gretl_cmatrix_get(zx, r, i, j);
		}
	    }
	    ret->is_complex = 1;
	}
    }

    return ret;
}

gretl_matrix *gretl_cmatrix_unvech (const gretl_matrix *X,
				    int *err)
{
    gretl_matrix *ret = NULL;

    if (!cmatrix_validate(X, 0) || X->cols != 1) {
	*err = E_INVARG;
    } else {
	int r = X->rows / 2;
	int n = (int) ((sqrt(1.0 + 8.0 * r) - 1.0) / 2.0);

	ret = gretl_matrix_alloc(2*n, n);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    double complex *zx = (double complex *) X->val;
	    double complex *zu = (double complex *) ret->val;
	    double complex zk;
	    int i, j, k = 0;

	    for (j=0; j<n; j++) {
		for (i=j; i<n; i++) {
		    zk = zx[k++];
		    gretl_cmatrix_set(zu, n, i, j, zk);
		    gretl_cmatrix_set(zu, n, j, i, zk);
		}
	    }
	    ret->is_complex = 1;
	}
    }

    return ret;
}

/*
  Copies the values from row @is of @src into row @id
  of @dest, provided @src and @dest have the same number
  of columns.
*/

static void cmatrix_copy_row (gretl_matrix *dest, int id,
			      const gretl_matrix *src, int is)
{
    double complex *zd = (double complex *) dest->val;
    double complex *zs = (double complex *) src->val;
    double complex zj;
    int j, r = src->rows / 2;

    for (j=0; j<src->cols; j++) {
	zj = gretl_cmatrix_get(zs, r, is, j);
	gretl_cmatrix_set(zd, r, id, j, zj);
    }
}

gretl_matrix *gretl_cmatrix_reverse_rows (const gretl_matrix *X,
					  int *err)
{
    gretl_matrix *ret;
    int i, r, c;

    if (!cmatrix_validate(X, 0)) {
	*err = E_INVARG;
	return NULL;
    } else if (gretl_is_null_matrix(X)) {
	return gretl_null_matrix_new();
    }

    r = X->rows;
    c = X->cols;
    ret = gretl_matrix_alloc(r, c);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	r /= 2;
	for (i=0; i<r; i++) {
	    cmatrix_copy_row(ret, i, X, r-i-1);
	}
	ret->is_complex = 1;
    }

    return ret;
}

gretl_matrix *gretl_cmatrix_shape (const gretl_matrix *A,
				   int r, int c, int *err)
{
    gretl_matrix *B;
    double complex *az;
    double complex *bz;
    int i, k, nA, nB;

    if (gretl_is_null_matrix(A) || r < 0 || c < 0) {
	*err = E_INVARG;
	return NULL;
    }

    if (r == 0 && c == 0) {
	return gretl_null_matrix_new();
    }

    B = gretl_matrix_alloc(2*r, c);
    if (B == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    az = (double complex *) A->val;
    bz = (double complex *) B->val;
    nA = A->cols * A->rows/2;
    nB = r * c;
    B->is_complex = 1;

    k = 0;
    for (i=0; i<nB; i++) {
	bz[i] = az[k++];
	if (k == nA) {
	    k = 0;
	}
    }

    return B;
}

int gretl_cmatrix_zero_triangle (gretl_matrix *m, char t)
{
    double complex *zm;
    double complex z0;
    int i, j, r;

    if (!cmatrix_validate(m, 1)) {
	return E_INVARG;
    }

    zm = (double complex *) m->val;
    z0 = 0;
    r = m->rows / 2;

    if (t == 'U') {
	for (i=0; i<r-1; i++) {
	    for (j=i+1; j<m->cols; j++) {
		gretl_cmatrix_set(zm, r, i, j, z0);
	    }
	}
    } else {
	for (i=1; i<r; i++) {
	    for (j=0; j<i; j++) {
		gretl_cmatrix_set(zm, r, i, j, z0);
	    }
	}
    }

    return 0;
}

/* switch between "legacy" and new representations of a
   complex matrix */

gretl_matrix *gretl_cmatrix_switch (const gretl_matrix *m,
				    int to_new, int *err)
{
    gretl_matrix *ret = NULL;
    int r, c, i, j, k, jj;
    double x;

    if (gretl_is_null_matrix(m)) {
	return gretl_null_matrix_new();
    }

    if (to_new) {
	r = m->rows * 2;
	c = m->cols / 2;
    } else {
	r = m->rows / 2;
	c = m->cols * 2;
    }

    ret = gretl_matrix_alloc(r, c);
    if (ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    jj = 0;

    if (to_new) {
	/* make re and im components contiguous */
	for (j=0; j<ret->cols; j++) {
	    k = 0;
	    for (i=0; i<m->rows; i++) {
		x = gretl_matrix_get(m, i, jj);
		gretl_matrix_set(ret, k++, j, x);
		x = gretl_matrix_get(m, i, jj+1);
		gretl_matrix_set(ret, k++, j, x);
	    }
	    jj += 2;
	}
	ret->is_complex = 1;
    } else {
	/* put re and im components in adjacent columns */
	for (j=0; j<m->cols; j++) {
	    k = 0;
	    for (i=0; i<ret->rows; i++) {
		x = gretl_matrix_get(m, k++, j);
		gretl_matrix_set(ret, i, jj, x);
		x = gretl_matrix_get(m, k++, j);
		gretl_matrix_set(ret, i, jj+1, x);
	    }
	    jj += 2;
	}
    }

    return ret;
}

gretl_matrix *gretl_cmatrix_vector_stat (const gretl_matrix *m,
					 GretlVecStat vs, int rowwise,
					 int *err)
{
    double complex *rz, *mz;
    double complex z;
    gretl_matrix *ret;
    int r, c, i, j;
    int mr, rr;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    /* note: 8-byte rows */
    r = rowwise ? m->rows : 2;
    c = rowwise ? 1 : m->cols;

    ret = gretl_matrix_alloc(r, c);
    if (ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    ret->is_complex = 1;

    mz = (double complex *) m->val;
    rz = (double complex *) ret->val;

    mr = m->rows / 2;
    rr = r / 2;

    if (rowwise) {
	/* by rows */
	int jmin = vs == V_PROD ? 1 : 0;

	for (i=0; i<mr; i++) {
	    z = vs == V_PROD ? mz[0] : 0;
	    for (j=jmin; j<m->cols; j++) {
		if (vs == V_PROD) {
		    z *= gretl_cmatrix_get(mz, mr, i, j);
		} else {
		    z += gretl_cmatrix_get(mz, mr, i, j);
		}
	    }
	    if (vs == V_MEAN) {
		z /= m->cols;
	    }
	    gretl_cmatrix_set(rz, rr, i, 0, z);
	}
    } else {
	/* by columns */
	int imin = vs == V_PROD ? 1 : 0;

	for (j=0; j<m->cols; j++) {
	    z = vs == V_PROD ? mz[0] : 0;
	    for (i=imin; i<mr; i++) {
		if (vs == V_PROD) {
		    z *= gretl_cmatrix_get(mz, mr, i, j);
		} else {
		    z += gretl_cmatrix_get(mz, mr, i, j);
		}
	    }
	    if (vs == V_MEAN) {
		z /= mr;
	    }
	    gretl_cmatrix_set(rz, rr, 0, j, z);
	}
    }

    return ret;
}

int gretl_cmatrix_fill (gretl_matrix *m, double complex z)
{
    double complex *mz = (double complex *) m->val;
    int i, n = m->cols * m->rows / 2;

    for (i=0; i<n; i++) {
	mz[i] = z;
    }

    return 0;
}

gretl_matrix *scalar_to_complex (double x, int *err)
{
    gretl_matrix *m = gretl_matrix_alloc(2, 1);

    if (m != NULL) {
	m->val[0] = x;
	m->val[1] = 0;
	m->is_complex = 1;
    } else {
	*err = E_ALLOC;
    }

    return m;
}

static void vec_x_op_vec_y (double complex *z,
			    const double complex *x,
			    const double complex *y,
			    int n, int op)
{
    int i;

    switch (op) {
    case '*':
	for (i=0; i<n; i++) {
	    z[i] = x[i] * y[i];
	}
	break;
    case '/':
	for (i=0; i<n; i++) {
	    z[i] = x[i] / y[i];
	}
	break;
    case '+':
	for (i=0; i<n; i++) {
	    z[i] = x[i] + y[i];
	}
	break;
    case '-':
	for (i=0; i<n; i++) {
	    z[i] = x[i] - y[i];
	}
	break;
    case '^':
	for (i=0; i<n; i++) {
	    z[i] = cpow(x[i], y[i]);
	}
	break;
    case '=':
	for (i=0; i<n; i++) {
	    z[i] = x[i] == y[i];
	}
	break;
    case '!':
	for (i=0; i<n; i++) {
	    z[i] = x[i] != y[i];
	}
	break;
    default:
	break;
    }
}

static void vec_x_op_y (double complex *z,
			const double complex *x,
			double complex y,
			int n, int op)
{
    int i;

    switch (op) {
    case '*':
	for (i=0; i<n; i++) {
	    z[i] = x[i] * y;
	}
	break;
    case '/':
	for (i=0; i<n; i++) {
	    z[i] = x[i] / y;
	}
	break;
    case '+':
	for (i=0; i<n; i++) {
	    z[i] = x[i] + y;
	}
	break;
    case '-':
	for (i=0; i<n; i++) {
	    z[i] = x[i] - y;
	}
	break;
    case '^':
	for (i=0; i<n; i++) {
	    z[i] = cpow(x[i], y);
	}
	break;
    case '=':
	for (i=0; i<n; i++) {
	    z[i] = x[i] == y;
	}
	break;
    case '!':
	for (i=0; i<n; i++) {
	    z[i] = x[i] != y;
	}
	break;
    default:
	break;
    }
}

static void x_op_vec_y (double complex *z,
			double complex x,
			const double complex *y,
			int n, int op)
{
    int i;

    switch (op) {
    case '*':
	for (i=0; i<n; i++) {
	    z[i] = x * y[i];
	}
	break;
    case '/':
	for (i=0; i<n; i++) {
	    z[i] = x / y[i];
	}
	break;
    case '+':
	for (i=0; i<n; i++) {
	    z[i] = x + y[i];
	}
	break;
    case '-':
	for (i=0; i<n; i++) {
	    z[i] = x - y[i];
	}
	break;
    case '^':
	for (i=0; i<n; i++) {
	    z[i] = cpow(x, y[i]);
	}
	break;
    case '=':
	for (i=0; i<n; i++) {
	    z[i] = x == y[i];
	}
	break;
    case '!':
	for (i=0; i<n; i++) {
	    z[i] = x != y[i];
	}
	break;
    default:
	break;
    }
}

static double complex x_op_y (double complex x,
			      double complex y,
			      int op)
{
    switch (op) {
    case '*':
	return x * y;
    case '/':
	return x / y;
    case '+':
	return x + y;
    case '-':
	return x - y;
    case '^':
	return cpow(x, y);
    case '=':
	return x == y;
    case '!':
	return x != y;
    default:
	return 0;
    }
}

static gretl_matrix *cmatrix_dot_op (const gretl_matrix *a,
				     const gretl_matrix *b,
				     int op, int *err)
{
    gretl_matrix *c = NULL;
    double complex *az, *bz, *cz;
    double complex x, y;
    int ar, br, cr;
    int nr, nc;
    int conftype;
    int i, j, off;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	*err = E_DATA;
	return NULL;
    }

    conftype = dot_operator_conf(a, b, &nr, &nc);

    if (conftype == CONF_NONE) {
	fputs("gretl_cmatrix_dot_op: matrices not conformable\n", stderr);
	fprintf(stderr, " op = '%c', A is %d x %d, B is %d x %d\n",
		(char) op, a->rows/2, a->cols, b->rows/2, b->cols);
	*err = E_NONCONF;
	return NULL;
    }

    c = gretl_matrix_alloc(nr * 2, nc);
    if (c == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    c->is_complex = 1;

    az = (double complex *) a->val;
    bz = (double complex *) b->val;
    cz = (double complex *) c->val;

    ar = a->rows / 2;
    br = b->rows / 2;
    cr = c->rows / 2;

#if 0 /* maybe expose this in gretl_matrix.h? */
    math_err_init();
#endif

    switch (conftype) {
    case CONF_ELEMENTS:
	vec_x_op_vec_y(cz, az, bz, nr*nc, op);
	break;
    case CONF_A_COLVEC:
	for (i=0; i<nr; i++) {
	    x = az[i];
	    for (j=0; j<nc; j++) {
		y = gretl_cmatrix_get(bz, br, i, j);
		y = x_op_y(x, y, op);
		gretl_cmatrix_set(cz, cr, i, j, y);
	    }
	}
	break;
    case CONF_B_COLVEC:
	for (i=0; i<nr; i++) {
	    y = bz[i];
	    for (j=0; j<nc; j++) {
		x = gretl_cmatrix_get(az, ar, i, j);
		x = x_op_y(x, y, op);
		gretl_cmatrix_set(cz, cr, i, j, x);
	    }
	}
	break;
    case CONF_A_ROWVEC:
	off = 0;
	for (j=0; j<nc; j++) {
	    x = az[j];
	    x_op_vec_y(cz + off, x, bz + off, nr, op);
	    off += nr;
	}
	break;
    case CONF_B_ROWVEC:
	off = 0;
	for (j=0; j<nc; j++) {
	    y = bz[j];
	    vec_x_op_y(cz + off, az + off, y, nr, op);
	    off += nr;
	}
	break;
    case CONF_A_SCALAR:
	x_op_vec_y(cz, az[0], bz, nr*nc, op);
	break;
    case CONF_B_SCALAR:
	vec_x_op_y(cz, az, bz[0], nr*nc, op);
	break;
    case CONF_AC_BR:
	for (i=0; i<nr; i++) {
	    x = az[i];
	    for (j=0; j<nc; j++) {
		y = bz[j];
		y = x_op_y(x, y, op);
		gretl_cmatrix_set(cz, cr, i, j, y);
	    }
	}
	break;
    case CONF_AR_BC:
	for (j=0; j<nc; j++) {
	    x = az[j];
	    for (i=0; i<nr; i++) {
		y = bz[i];
		y = x_op_y(x, y, op);
		gretl_cmatrix_set(cz, cr, i, j, y);
	    }
	}
	break;
    default: /* hush a warning */
	break;
    }

    if (errno) {
#if 0 /* not yet */
	*err = math_err_check("gretl_matrix_dot_op", errno);
#endif
	if (*err) {
	    gretl_matrix_free(c);
	    c = NULL;
	}
    }

    return c;
}

gretl_matrix *gretl_cmatrix_dot_op (const gretl_matrix *A,
				    const gretl_matrix *B,
				    int op, int *err)
{
    gretl_matrix *T = NULL;
    gretl_matrix *C = NULL;

    if (A->is_complex && B->is_complex) {
	C = cmatrix_dot_op(A, B, op, err);
    } else if (A->is_complex) {
	/* case of real B */
	T = complex_from_real(B, err);
	if (T != NULL) {
	    C = cmatrix_dot_op(A, T, op, err);
	}
    } else if (B->is_complex) {
	/* case of real A */
	T = complex_from_real(A, err);
	if (T != NULL) {
	    C = cmatrix_dot_op(T, B, op, err);
	}
    } else {
	*err = E_TYPES;
    }

    gretl_matrix_free(T);

    return C;
}
