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

#define cscalar(m) (m->rows == 1 && m->cols == 1)

/* helper function for fftw-based real FFT functions */

static int fft_allocate (double **ffx, double complex **ffz,
			 gretl_matrix **ret, int r, int c,
			 int inverse, int newstyle)
{
    /* real workspace, per series */
    *ffx = fftw_malloc(r * sizeof **ffx);
    if (*ffx == NULL) {
	return E_ALLOC;
    }

    /* complex workspace, per series */
    *ffz = fftw_malloc((r/2 + 1 + r % 2) * sizeof **ffz);
    if (*ffz == NULL) {
	free(*ffx);
	return E_ALLOC;
    }

    /* matrix to hold output */
    if (newstyle && !inverse) {
	*ret = gretl_cmatrix_new(r, c);
    } else {
	*ret = gretl_matrix_alloc(r, c);
    }
    if (*ret == NULL) {
	free(*ffx);
	free(*ffz);
	return E_ALLOC;
    }

    return 0;
}

/* FFT for real input -> complex output and
   FFTI for Hermetian input -> real output.
   Both old and new-style complex formats
   are supported.
*/

static gretl_matrix *
real_matrix_fft (const gretl_matrix *y, int inverse,
		 int newstyle, int *err)
{
    gretl_matrix *ret = NULL;
    fftw_plan p = NULL;
    double *ffx = NULL;
    double complex *ffz = NULL;
    double xr, xi;
    int r, c, m, odd, cr, ci;
    int incols, outcols;
    int i, j;

    if (y->rows < 2) {
	*err = E_DATA;
	return NULL;
    }

    r = y->rows;
    m = r / 2;
    odd = r % 2;

    incols = y->cols;
    if (newstyle) {
	/* the number of columns is invariant wrt real vs complex */
	outcols = incols;
    } else {
	/* the number of columns is double for complex what it is
	   for real */
	outcols = inverse ? incols / 2 : incols * 2;
    }

    *err = fft_allocate(&ffx, &ffz, &ret, r, outcols,
			inverse, newstyle);
    if (*err) {
	return NULL;
    }

    c = MIN(incols, outcols);
    cr = 0;
    ci = 1;

    for (j=0; j<c; j++) {
	/* load the data */
	if (newstyle && inverse) {
	    for (i=0; i<=m+odd; i++) {
		ffz[i] = gretl_cmatrix_get(y, i, j);
	    }
	} else if (inverse) {
	    for (i=0; i<=m+odd; i++) {
		xr = gretl_matrix_get(y, i, cr);
		xi = gretl_matrix_get(y, i, ci);
		ffz[i] = xr + xi * I;
	    }
	} else {
	    /* going in the real -> complex direction */
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
		gretl_matrix_set(ret, i, j, ffx[i] / r);
	    }
	} else if (newstyle) {
	    for (i=0; i<=m+odd; i++) {
		gretl_cmatrix_set(ret, i, j, ffz[i]);
	    }
	    for (i=m; i>0; i--) {
		gretl_cmatrix_set(ret, r-i, j, conj(ffz[i]));
	    }
	} else {
	    for (i=0; i<=m+odd; i++) {
		gretl_matrix_set(ret, i, cr, creal(ffz[i]));
		gretl_matrix_set(ret, i, ci, cimag(ffz[i]));
	    }
	    for (i=m; i>0; i--) {
		gretl_matrix_set(ret, r-i, cr,  creal(ffz[i]));
		gretl_matrix_set(ret, r-i, ci, -cimag(ffz[i]));
	    }
	}
	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(ffz);
    fftw_free(ffx);

    return ret;
}

static int row_is_real (const gretl_matrix *y, int i)
{
    double complex z;
    int j;

    for (j=0; j<y->cols; j++) {
	z = gretl_cmatrix_get(y, i, j);
	if (cimag(z) != 0) {
	    return 0;
	}
    }

    return 1;
}

static int rows_are_conjugate (const gretl_matrix *y, int i, int k)
{
    double complex z1, z2;
    int j;

    for (j=0; j<y->cols; j++) {
	z1 = gretl_cmatrix_get(y, i, j);
	z2 = gretl_cmatrix_get(y, k, j);
	if (z1 != conj(z2)) {
	    return 0;
	}
    }

    return 1;
}

static int fft_is_hermitian (const gretl_matrix *y)
{
    if (!row_is_real(y, 0)) {
	return 0;
    } else {
	int i, k, n2 = y->rows / 2;

	if (y->rows % 2 == 0 && !row_is_real(y, n2)) {
	    /* In the Hermitian case the middle row of the
	       odd-numbered remainder, on excluding the first
	       row, must be real.
	    */
	    return 0;
	}
	for (i=1; i<n2; i++) {
	    k = y->rows - i;
	    if (!rows_are_conjugate(y, i, k)) {
		return 0;
	    }
	}
    }

    return 1;
}

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int cmat, int *err)
{
    return real_matrix_fft(y, 0, cmat, err);
}

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err)
{
    if (y->is_complex) {
	/* new-style */
	if (fft_is_hermitian(y)) {
	    /* its inverse should be real */
	    return real_matrix_fft(y, 1, 1, err);
	} else {
	    return gretl_cmatrix_fft(y, 1, err);
	}
    } else {
	/* old-style */
	return real_matrix_fft(y, 1, 0, err);
    }
}

static int cmatrix_validate (const gretl_matrix *m, int square)
{
    int ret = 1;

    if (gretl_is_null_matrix(m)) {
	ret = 0;
    } else if (!m->is_complex || m->z == NULL) {
	ret = 0;
    } else if (square && m->rows != m->cols) {
	ret = 0;
    }

    if (!ret) {
	fprintf(stderr, "cmatrix_validate: failed\n");
    }

    return ret;
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

    C = gretl_cmatrix_new0(A->rows, A->cols);

    if (C == NULL) {
	*err = E_ALLOC;
    } else {
	int i, n = A->rows * A->cols;

	for (i=0; i<n; i++) {
	    C->z[i] = A->val[i];
	}
    }

    return C;
}

/* Multiplication of complex matrices via BLAS zgemm(),
   allowing for conjugate transposition of @A or @B.
*/

static gretl_matrix *gretl_zgemm (const gretl_matrix *A,
				  GretlMatrixMod amod,
				  const gretl_matrix *B,
				  GretlMatrixMod bmod,
				  int *err)
{
    gretl_matrix *C;
    cmplx *a, *b, *c;
    cmplx alpha = {1, 0};
    cmplx beta = {0, 0};
    char transa = 'N';
    char transb = 'N';
    integer lda, ldb, ldc;
    integer m, n, k, kb;

    if (!cmatrix_validate(A, 0) || !cmatrix_validate(B, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    if (amod == GRETL_MOD_CTRANSP) {
	transa = 'C';
	m = A->cols; /* rows of op(A) */
	k = A->rows; /* cols of op(A) */
    } else {
	m = A->rows; /* complex rows of A */
	k = A->cols; /* cols of A */
    }

    if (bmod == GRETL_MOD_CTRANSP) {
	transb = 'C';
	n = B->rows;  /* columns of op(B) */
	kb = B->cols; /* rows of op(B) */
    } else {
	n = B->cols;
	kb = B->rows;
    }

    if (k != kb) {
	*err = E_NONCONF;
	return NULL;
    }

    C = gretl_cmatrix_new(m, n);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    a = (cmplx *) A->val;
    b = (cmplx *) B->val;
    c = (cmplx *) C->val;

    lda = A->rows;
    ldb = B->rows;
    ldc = C->rows;

    zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda,
	   b, &ldb, &beta, c, &ldc);

    return C;
}

static gretl_matrix *real_cmatrix_multiply (const gretl_matrix *A,
					    const gretl_matrix *B,
					    GretlMatrixMod amod,
					    int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *T = NULL;

    if (A->is_complex && B->is_complex) {
	if (amod == 0 && (cscalar(A) || cscalar(B))) {
	    return gretl_cmatrix_dot_op(A, B, '*', err);
	} else {
	    C = gretl_zgemm(A, amod, B, 0, err);
	}
    } else if (A->is_complex) {
	/* case of real B */
	T = complex_from_real(B, err);
	if (T != NULL) {
	    C = gretl_zgemm(A, amod, T, 0, err);
	}
    } else if (B->is_complex) {
	/* case of real A */
	T = complex_from_real(A, err);
	if (T != NULL) {
	    C = gretl_zgemm(T, amod, B, 0, err);
	}
    } else {
	*err = E_TYPES;
    }

    gretl_matrix_free(T);

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
    return real_cmatrix_multiply(A, B, 0, err);
}

/* Returns (conjugate transpose of A, or A^H) times B,
   allowing for the possibility that either A or B (but
   not both!) may be a real matrix on input.
*/

gretl_matrix *gretl_cmatrix_AHB (const gretl_matrix *A,
				 const gretl_matrix *B,
				 int *err)
{
    return real_cmatrix_multiply(A, B, GRETL_MOD_CTRANSP, err);
}

/* Eigen decomposition of complex (Hermitian) matrix using
   LAPACK's zheev() */

gretl_matrix *gretl_zheev (gretl_matrix *A, int eigenvecs,
			   int *err)
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

    n = A->rows;
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

static int zgeev_eigvecs_alloc (gretl_matrix *m,
				gretl_matrix **pev,
				cmplx **evz,
				int n)
{
    gretl_matrix *ev = NULL;
    int mrc = m->rows * m->cols;
    int dim = n * n;

    /* We need an n x n complex matrix for output:
       is @m usable or do we need to allocate a
       new matrix?
    */
    if (m->is_complex && mrc == dim) {
	m->rows = m->cols = n;
	*evz = (cmplx *) m->val;
    } else if (!m->is_complex && mrc == 2*dim) {
	m->rows = 2*n;
	m->cols = n;
	gretl_matrix_set_complex_full(m, 1);
	*evz = (cmplx *) m->val;
    } else {
	/* have to allocate */
	ev = gretl_cmatrix_new0(n, n);
	if (ev == NULL) {
	    return E_ALLOC;
	} else {
	    *pev = ev;
	    *evz = (cmplx *) ev->val;
	}
    }

    return 0;
}

/* Eigen decomposition of complex (non-Hermitian) matrix using
   LAPACK's zgeev() */

gretl_matrix *gretl_zgeev (const gretl_matrix *A,
			   gretl_matrix *VR,
			   gretl_matrix *VL,
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

    n = A->rows;
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
	*err = zgeev_eigvecs_alloc(VL, &Ltmp, &vl, n);
	if (*err) {
	    goto bailout;
	}
    }

    if (VR != NULL) {
	/* right eigenvectors wanted */
	*err = zgeev_eigvecs_alloc(VR, &Rtmp, &vr, n);
	if (*err) {
	    goto bailout;
	}
    }

    rwork = malloc(2 * n * sizeof *rwork);
    ret = gretl_cmatrix_new(n, 1);
    if (rwork == NULL || ret == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    w = ret->val;

    /* get optimal workspace size */
    lwork = -1;
    zgeev_(&jobvl, &jobvr, &n, a, &n, w, vl, &ldvl, vr, &ldvr,
	   &wsz, &lwork, rwork, &info);
    lwork = (integer) wsz.r;
    work = malloc(lwork * sizeof *work);
    if (work == NULL) {
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
	}
	if (Rtmp != NULL) {
	    gretl_matrix_replace_content(VR, Rtmp);
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

gretl_matrix *gretl_cmatrix_inverse (const gretl_matrix *A, int *err)
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

    if (!cmatrix_validate(x, 0)) {
	return E_INVARG;
    }

    lda = m = x->rows;
    n = x->cols;
    xsize = lda * n;

    if (smod == SVD_THIN && m < n) {
	fprintf(stderr, "gretl_cmatrix_SVD: X is %d x %d, should be 'thin'\n",
		x->rows, x->cols);
	return E_NONCONF;
    }

    az = malloc(xsize * sizeof *az);
    if (az == NULL) {
	return E_ALLOC;
    }
    memcpy(az, x->val, xsize * sizeof *az);

    k = MIN(m, n);

    s = gretl_vector_alloc(k);
    if (s == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (pu != NULL) {
	ldu = m;
	if (smod == SVD_FULL) {
	    u = gretl_cmatrix_new(ldu, m);
	} else {
	    u = gretl_cmatrix_new(ldu, n);
	}
	if (u == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
	    uval = (cmplx *) u->val;
	    jobu = (smod == SVD_FULL)? 'A' : 'S';
	}
    }

    if (pvt != NULL) {
	ldvt = n;
	vt = gretl_cmatrix_new(ldvt, n);
	if (vt == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
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

static double csvd_smin (const gretl_matrix *a, double smax)
{
    const double macheps = 2.20e-16;
    int dmax = (a->rows > a->cols)? a->rows : a->cols;

    /* as per numpy, Matlab (2015-09-28) */
    return dmax * macheps * smax;
}

int gretl_cmatrix_rank (const gretl_matrix *A, int *err)
{
    gretl_matrix *S = NULL;
    int i, k, rank = 0;

    if (!cmatrix_validate(A, 0)) {
	*err = E_INVARG;
	return 0;
    }

    k = MIN(A->rows, A->cols);

    if (A->rows > 4 * k || A->cols > 4 * k) {
	GretlMatrixMod mod1, mod2;
	gretl_matrix *B;

	mod1 = A->rows > k ? GRETL_MOD_CTRANSP : 0;
	mod2 = A->cols > k ? GRETL_MOD_CTRANSP : 0;
	B = gretl_zgemm(A, mod1, A, mod2, err);
	if (!*err) {
	    *err = gretl_cmatrix_SVD(B, NULL, &S, NULL, 1);
	}
	gretl_matrix_free(B);
    } else {
	*err = gretl_cmatrix_SVD(A, NULL, &S, NULL, 1);
    }

    if (!*err) {
	double smin = csvd_smin(A, S->val[0]);

	for (i=0; i<k; i++) {
	    if (S->val[i] > smin) {
		rank++;
	    }
	}
    }

    gretl_matrix_free(S);

    return rank;
}

/* generalized inverse of a complex matrix via its SVD */

gretl_matrix *gretl_cmatrix_ginv (const gretl_matrix *A, int *err)
{
    gretl_matrix *U = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *s = NULL;
    gretl_matrix *Vt = NULL;
    gretl_matrix *ret = NULL;

    *err = gretl_cmatrix_SVD(A, &U, &s, &V, 1);

    if (!*err) {
	double complex vij;
	int i, j, h = 0;

	for (i=0; i<s->cols; i++) {
	    h += s->val[i] > 1.0e-13;
	}
	Vt = gretl_ctrans(V, 1, err);
	if (!*err) {
	    for (j=0; j<h; j++) {
		for (i=0; i<Vt->rows; i++) {
		    vij = gretl_cmatrix_get(Vt, i, j);
		    gretl_cmatrix_set(Vt, i, j, vij / s->val[j]);
		}
	    }
	}
	if (!*err) {
	    Vt->cols = U->cols = h;
	    ret = gretl_zgemm(Vt, 0, U, GRETL_MOD_CTRANSP, err);
	}
    }

    gretl_matrix_free(U);
    gretl_matrix_free(V);
    gretl_matrix_free(s);
    gretl_matrix_free(Vt);

    return ret;
}

/* Horizontal direct product of complex @A and @B */

static gretl_matrix *real_cmatrix_hdp (const gretl_matrix *A,
				       const gretl_matrix *B,
				       int *err)
{
    gretl_matrix *C = NULL;
    int r, p, q;

    if (!cmatrix_validate(A,0) || !cmatrix_validate(B,0)) {
	*err = E_INVARG;
	return NULL;
    }

    if (B->rows != A->rows) {
	*err = E_NONCONF;
	return NULL;
    }

    r = A->rows;
    p = A->cols;
    q = B->cols;

    C = gretl_cmatrix_new0(r, p*q);

    if (C == NULL) {
	*err = E_ALLOC;
    } else {
	double complex aij, bik;
	int i, j, k, joff;

	for (i=0; i<r; i++) {
	    for (j=0; j<p; j++) {
		aij = gretl_cmatrix_get(A, i, j);
		if (aij != 0.0) {
		    joff = j * q;
		    for (k=0; k<q; k++) {
			bik = gretl_cmatrix_get(B, i, k);
			gretl_cmatrix_set(C, i, joff + k, aij*bik);
		    }
		}
	    }
	}
    }

    return C;
}

/* Kronecker product of complex @A and @B */

static gretl_matrix *real_cmatrix_kron (const gretl_matrix *A,
					const gretl_matrix *B,
					int *err)
{
    gretl_matrix *K = NULL;
    int p, q, r, s;

    if (!cmatrix_validate(A,0) || !cmatrix_validate(B,0)) {
	*err = E_INVARG;
	return NULL;
    }

    p = A->rows;
    q = A->cols;
    r = B->rows;
    s = B->cols;

    K = gretl_cmatrix_new0(p*r, q*s);

    if (K == NULL) {
	*err = E_ALLOC;
    } else {
	double complex x, aij, bkl;
	int i, j, k, l;
	int ioff, joff;
	int Ki, Kj;

	for (i=0; i<p; i++) {
	    ioff = i * r;
	    for (j=0; j<q; j++) {
		/* block ij is an r * s matrix, a_{ij} * B */
		aij = gretl_cmatrix_get(A, i, j);
		joff = j * s;
		for (k=0; k<r; k++) {
		    Ki = ioff + k;
		    for (l=0; l<s; l++) {
			bkl = gretl_cmatrix_get(B, k, l);
			Kj = joff + l;
			x = aij * bkl;
			gretl_cmatrix_set(K, Ki, Kj, x);
		    }
		}
	    }
	}
    }

    return K;
}

gretl_matrix *gretl_cmatrix_kronlike (const gretl_matrix *A,
				      const gretl_matrix *B,
				      int hdp, int *err)
{
    gretl_matrix *L = (gretl_matrix *) A;
    gretl_matrix *R = (gretl_matrix *) B;
    gretl_matrix *C = NULL;

    if (A->is_complex && B->is_complex) {
	; /* OK */
    } else if (A->is_complex) {
	R = complex_from_real(B, err);
    } else if (B->is_complex) {
	L = complex_from_real(A, err);
    } else {
	*err = E_TYPES;
    }

    if (!*err) {
	if (hdp) {
	    C = real_cmatrix_hdp(L, R, err);
	} else {
	    C = real_cmatrix_kron(L, R, err);
	}
    }

    if (L != A) gretl_matrix_free(L);
    if (R != B) gretl_matrix_free(R);

    return C;
}

gretl_matrix *gretl_cmatrix_hdprod (const gretl_matrix *A,
				    const gretl_matrix *B,
				    int *err)
{
    return gretl_cmatrix_kronlike(A, B, 1, err);
}

gretl_matrix *gretl_cmatrix_kronecker (const gretl_matrix *A,
				       const gretl_matrix *B,
				       int *err)
{
    return gretl_cmatrix_kronlike(A, B, 0, err);
}

/* Complex FFT (or inverse) via fftw */

gretl_matrix *gretl_cmatrix_fft (const gretl_matrix *A,
				 int inverse, int *err)
{
    gretl_matrix *B = NULL;
    double complex *tmp, *ptr;
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

    r = A->rows;
    c = A->cols;

    tmp = (double complex *) B->val;
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

#define CDEC5 1 /* five decimal places in default format */

int gretl_cmatrix_print_range (const gretl_matrix *A,
			       const char *name,
			       int rmin, int rmax,
			       PRN *prn)
{
    double complex aij;
    double re, im, ai, xmax;
    int alt_default = 0;
    int all_ints = 1;
    int intmax = 1000;
#if CDEC5
    int altmin = 1000;
#else
    int altmin = 100;
#endif
    int zwidth = 0;
    int rdecs, idecs;
    char pm;
    int r, c, i, j;

    if (!cmatrix_validate(A, 0)) {
	return E_INVARG;
    }

    r = A->rows;
    c = A->cols;

    if (rmin < 0) rmin = 0;
    if (rmax < 0) rmax = r;

    xmax = 0;

    for (j=0; j<c; j++) {
	for (i=rmin; i<rmax; i++) {
	    aij = gretl_cmatrix_get(A, i, j);
	    re = creal(aij);
	    im = cimag(aij);
	    if (all_ints && (floor(re) != re || floor(im) != im)) {
		all_ints = 0;
	    }
	    re = MAX(fabs(re), fabs(im));
	    if (re > xmax) {
		xmax = re;
	    }
	    if (all_ints && xmax >= intmax) {
		all_ints = 0;
	    }
	    if (!all_ints && xmax >= altmin) {
		break;
	    }
	}
	if (!all_ints && xmax >= altmin) {
	    break;
	}
    }

    if (all_ints) {
	/* apply a more compact format */
	zwidth = 2 + (xmax >= 10) + (xmax >= 100);
    } else if (xmax >= altmin) {
	/* we'll want a different default format */
	alt_default = 1;
    }

    if (name != NULL && *name != '\0') {
	pprintf(prn, "%s (%d x %d)", name, r, c);
	pputs(prn, "\n\n");
    }

    for (i=rmin; i<rmax; i++) {
	for (j=0; j<c; j++) {
	    aij = gretl_cmatrix_get(A, i, j);
	    re = creal(aij);
	    im = cimag(aij);
	    pm = (im >= -0) ? '+' : '-';
	    ai = fabs(im);
	    if (zwidth > 0) {
		pprintf(prn, "%*g %c %*gi", zwidth, re, pm, zwidth-1, ai);
	    } else if (alt_default) {
		pprintf(prn, "%# 9.4e %c %#8.4ei", re, pm, ai);
	    } else {
#if CDEC5
		rdecs = re <= -10 ? 4 : 5;
		idecs = ai >= 10 ? 4 : 5;
		pprintf(prn, "%8.*f %c %7.*fi", rdecs, re, pm, idecs, ai);
#else
		rdecs = re <= -10 ? 3 : 4;
		idecs = ai >= 10 ? 3 : 4;
		pprintf(prn, "%7.*f %c %6.*fi", rdecs, re, pm, idecs, ai);
#endif
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

int gretl_cmatrix_print (const gretl_matrix *A,
			 const char *name,
			 PRN *prn)
{
    return gretl_cmatrix_print_range(A, name, -1, -1, prn);
}

int gretl_cmatrix_printf (const gretl_matrix *A,
			  const char *fmt,
			  PRN *prn)
{
    double complex aij;
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

    r = A->rows;
    c = A->cols;

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    aij = gretl_cmatrix_get(A, i, j);
	    re = creal(aij);
	    im = cimag(aij);
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

gretl_matrix *gretl_cmatrix_build (const gretl_matrix *Re,
				   const gretl_matrix *Im,
				   double ival, int *err)
{
    gretl_matrix *C = NULL;
    int i, n;

    if (gretl_is_null_matrix(Re) || Re->is_complex) {
	*err = E_INVARG;
    } else if (Im != NULL) {
	if (Im->rows != Re->rows || Im->cols != Re->cols) {
	    *err = E_NONCONF;
	} else if (Im->is_complex) {
	    *err = E_INVARG;
	}
    }

    if (!*err) {
	n = Re->rows * Re->cols;
	C = gretl_cmatrix_new(Re->rows, Re->cols);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	for (i=0; i<n; i++) {
	    if (Im != NULL) {
		ival = Im->val[i];
	    }
	    C->z[i] = Re->val[i] + ival * I;
	}
    }

    return C;
}

/* Extract the real part of @A if @im = 0 or the imaginary part
   if @im is non-zero.
*/

gretl_matrix *gretl_cmatrix_extract (const gretl_matrix *A,
				     int im, int *err)
{
    gretl_matrix *B = NULL;
    int i, n;

    if (!cmatrix_validate(A, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    B = gretl_matrix_alloc(A->rows, A->cols);

    if (B == NULL) {
	*err = E_ALLOC;
    } else {
	n = A->rows * A->cols;
	for (i=0; i<n; i++) {
	    B->val[i] = im ? cimag(A->z[i]) : creal(A->z[i]);
	}
    }

    return B;
}

/* Return [conjugate] transpose of complex matrix @A */

gretl_matrix *gretl_ctrans (const gretl_matrix *A,
			    int conjugate, int *err)
{
    gretl_matrix *C = NULL;
    double complex aij;
    int i, j;

    if (!cmatrix_validate(A, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    C = gretl_cmatrix_new(A->cols, A->rows);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (j=0; j<A->cols; j++) {
	for (i=0; i<A->rows; i++) {
	    aij = gretl_cmatrix_get(A, i, j);
	    if (conjugate) {
		aij = conj(aij);
	    }
	    gretl_cmatrix_set(C, j, i, aij);
	}
    }

    return C;
}

/* Convert complex matrix @A to its conjugate transpose */

int gretl_ctrans_in_place (gretl_matrix *A)
{
    gretl_matrix *C = NULL;
    double complex aij;
    int i, j, n;
    int err = 0;

    if (!cmatrix_validate(A, 0)) {
	return E_INVARG;
    }

    /* temporary matrix */
    C = gretl_cmatrix_new(A->cols, A->rows);
    if (C == NULL) {
	return E_ALLOC;
    }

    for (j=0; j<A->cols; j++) {
	for (i=0; i<A->rows; i++) {
	    aij = gretl_cmatrix_get(A, i, j);
	    gretl_cmatrix_set(C, j, i, conj(aij));
	}
    }

    /* now rejig @A */
    A->rows = C->rows;
    A->cols = C->cols;
    n = C->rows * C->cols;
    memcpy(A->z, C->z, n * sizeof *A->z);
    gretl_matrix_destroy_info(A);

    gretl_matrix_free(C);

    return err;
}

/* Addition or subtraction of matrices: handle the case
   where one operand is complex and the other is real.
   In addition handle the case where one of them is a
   scalar matrix, real or complex.
*/

gretl_matrix *gretl_cmatrix_add_sub (const gretl_matrix *A,
				     const gretl_matrix *B,
				     int sgn, int *err)
{
    gretl_matrix *C = NULL;
    double complex aval = 0;
    double complex bval = 0;
    int cr = A->rows;
    int cc = A->cols;
    int a_scalar = 0;
    int b_scalar = 0;

    if (A->is_complex && B->rows == 1 && B->cols == 1) {
	b_scalar = 1;
	bval = B->is_complex ? B->z[0] : B->val[0];
	goto allocate;
    } else if (B->is_complex && A->rows == 1 && A->cols == 1) {
	cr = B->rows;
	cc = B->cols;
	a_scalar = 1;
	aval = A->is_complex ? A->z[0] : A->val[0];
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
	if (B->rows != cr) {
	    *err = E_NONCONF;
	}
    } else {
	/* A real, B complex */
	cr = B->rows;
	if (A->rows != cr) {
	    *err = E_NONCONF;
	}
    }

 allocate:

    if (!*err) {
	C = gretl_cmatrix_new(cr, cc);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	int i, n = cc * cr;

	if (b_scalar) {
	    for (i=0; i<n; i++) {
		C->z[i] = sgn < 0 ? A->z[i] - bval : A->z[i] + bval;
	    }
	} else if (a_scalar) {
	    for (i=0; i<n; i++) {
		C->z[i] = sgn < 0 ? aval - B->z[i] : aval + B->z[i];
	    }
	} else if (A->is_complex && B->is_complex) {
	    for (i=0; i<n; i++) {
		C->z[i] = sgn < 0 ? A->z[i] - B->z[i] : A->z[i] + B->z[i];
	    }
	} else if (A->is_complex) {
	    for (i=0; i<n; i++) {
		C->z[i] = sgn < 0 ? A->z[i] - B->val[i] : A->z[i] + B->val[i];
	    }
	} else {
	    for (i=0; i<n; i++) {
		C->z[i] = sgn < 0 ? A->val[i] - B->z[i] : A->val[i] + B->z[i];
	    }
	}
    }

    return C;
}

/* Apply a function which maps from complex to real:
   creal, cimag, carg, cmod.
*/

int apply_cmatrix_dfunc (gretl_matrix *targ,
			const gretl_matrix *src,
			double (*dfunc) (double complex))
{
    int n = src->cols * src->rows;
    int i, err = 0;

    if (!cmatrix_validate(src, 0) || targ->is_complex) {
	return E_INVARG;
    }

    errno = 0;

    for (i=0; i<n; i++) {
	targ->val[i] = dfunc(src->z[i]);
    }

    if (errno) {
	gretl_errmsg_set_from_errno(NULL, errno);
	err = E_DATA;
    }

    return err;
}

/* Apply a function which maps from complex to complex;
   includes trigonometric functions, log, exp.
*/

int apply_cmatrix_cfunc (gretl_matrix *targ,
			 const gretl_matrix *src,
			 double complex (*cfunc) (double complex))
{
    int n = src->cols * src->rows;
    int i, err = 0;

    if (!cmatrix_validate(src, 0) || !cmatrix_validate(targ, 0)) {
	return E_INVARG;
    }

    errno = 0;

    for (i=0; i<n; i++) {
	targ->z[i] = cfunc(src->z[i]);
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
    int i, n = src->cols * src->rows;
    int err;

    if (!cmatrix_validate(src, 0) || !cmatrix_validate(targ, 0)) {
	return E_INVARG;
    }

    for (i=0; i<n && !err; i++) {
	if (op == 1) {
	    /* U_NEG */
	    targ->z[i] = -src->z[i];
	} else if (op == 2) {
	    /* U_POS */
	    targ->z[i] = src->z[i];
	} else if (op == 3) {
	    /* U_NOT */
	    targ->z[i] = (src->z[i] == 0);
	} else {
	    err = E_INVARG;
	}
    }

    return err;
}

gretl_matrix *gretl_cmatrix_from_scalar (double complex z, int *err)
{
    gretl_matrix *ret = gretl_cmatrix_new(1, 1);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	ret->z[0] = z;
    }

    return ret;
}

/* Determinant via eigenvalues */

gretl_matrix *gretl_cmatrix_determinant (const gretl_matrix *X,
					 int log, int *err)
{
    gretl_matrix *E = NULL;
    gretl_matrix *ret = NULL;

    if (!cmatrix_validate(X, 1)) {
	*err = E_INVARG;
	return ret;
    }

    E = gretl_zgeev(X, NULL, NULL, err);

    if (E != NULL) {
	double complex cret = 1;
	int i;

	for (i=0; i<E->rows; i++) {
	    cret *= E->z[i];
	}
	gretl_matrix_free(E);
	if (log) {
	    cret = clog(cret);
	}
	ret = gretl_cmatrix_from_scalar(cret, err);
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
	double complex tr = 0;
	int i;

	for (i=0; i<X->rows; i++) {
	    tr += gretl_cmatrix_get(X, i, i);
	}
	ret = gretl_cmatrix_from_scalar(tr, err);
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
    double complex zi = 0;
    int d, i;
    int match = 0;
    int err = 0;

    if (!cmatrix_validate(targ, 0)) {
	return E_INVARG;
    }

    d = MIN(targ->rows, targ->cols);

    if (src != NULL) {
	if (src->is_complex) {
	    if (gretl_vector_get_length(src) == d) {
		/* conformable complex vector */
		match = 1;
	    } else if (cscalar(src)) {
		/* complex scalar */
		zi = src->z[0];
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

    for (i=0; i<d; i++) {
	if (match == 1) {
	    gretl_cmatrix_set(targ, i, i, src->z[i]);
	} else if (match == 3) {
	    gretl_cmatrix_set(targ, i, i, src->val[i]);
	} else {
	    gretl_cmatrix_set(targ, i, i, zi);
	}
    }

    return err;
}

/* Set the lower or upper part of square complex matrix
   @targ using either @src (if not NULL) or @x. In the first
   case @src can be either a complex vector of the right length,
   or a real vector, or a complex scalar.
*/

int gretl_cmatrix_set_triangle (gretl_matrix *targ,
				const gretl_matrix *src,
				double x, int upper)
{
    double complex zi = 0;
    int r, c, p, i, j, n;
    int match = 0;
    int err = 0;

    if (!cmatrix_validate(targ, 0)) {
	return E_INVARG;
    }

    r = targ->rows;
    c = targ->cols;

    if ((c == 1 && upper) || (r == 1 && !upper)) {
	/* no such part */
	return E_INVARG;
    }

    p = MIN(r, c);
    n = (p * (p-1)) / 2;

    if (r > c && !upper) {
	n += (r - c) * c;
    } else if (c > r && upper) {
	n += (c - r) * r;
    }

    if (src != NULL) {
	if (src->is_complex) {
	    if (gretl_vector_get_length(src) == n) {
		/* conformable complex vector */
		match = 1;
	    } else if (cscalar(src)) {
		/* complex scalar */
		zi = src->z[0];
		match = 2;
	    }
	} else if (gretl_vector_get_length(src) == n) {
	    /* conformable real vector */
	    match = 3;
	}
    } else {
	/* use real scalar, @x */
	zi = x;
	match = 2;
    }

    if (match == 0) {
	return E_NONCONF;
    } else {
	int jmin = upper ? 1 : 0;
	int jmax = upper ? c : r;
	int imin = upper ? 0 : 1;
	int imax = upper ? 1 : r;
	int k = 0;

	for (j=jmin; j<jmax; j++) {
	    for (i=imin; i<imax; i++) {
		if (match == 1) {
		    gretl_cmatrix_set(targ, i, j, src->z[k++]);
		} else if (match == 3) {
		    gretl_cmatrix_set(targ, i, j, src->val[k++]);
		} else {
		    gretl_cmatrix_set(targ, i, j, zi);
		}
	    }
	    if (!upper) {
		imin++;
	    } else if (imax < r) {
		imax++;
	    }
	}
    }

    return err;
}

/* switch between "legacy" and new representations of a
   complex matrix */

gretl_matrix *gretl_cmatrix_switch (const gretl_matrix *m,
				    int to_new, int *err)
{
    gretl_matrix *ret = NULL;
    int r, c, i, j, k;

    if (gretl_is_null_matrix(m)) {
	return gretl_null_matrix_new();
    }

    if ((to_new && m->is_complex) ||
	(!to_new && !cmatrix_validate(m, 0))) {
	*err = E_INVARG;
	return NULL;
    }

    r = m->rows;
    c = to_new ? m->cols / 2 : m->cols * 2;

    if (to_new) {
	ret = gretl_cmatrix_new(r, c);
    } else {
	ret = gretl_matrix_alloc(r, c);
    }
    if (ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    k = 0;

    if (to_new) {
	/* make re and im components contiguous */
	double xr, xi;

	for (j=0; j<ret->cols; j++) {
	    for (i=0; i<m->rows; i++) {
		xr = gretl_matrix_get(m, i, k);
		xi = gretl_matrix_get(m, i, k+1);
		gretl_cmatrix_set(ret, i, j, xr + xi * I);
	    }
	    k += 2;
	}
    } else {
	/* put re and im components in adjacent columns */
	double complex mij;

	for (j=0; j<m->cols; j++) {
	    for (i=0; i<ret->rows; i++) {
		mij = gretl_cmatrix_get(m, i, j);
		gretl_matrix_set(ret, i, k, creal(mij));
		gretl_matrix_set(ret, i, k+1, cimag(mij));
	    }
	    k += 2;
	}
    }

    return ret;
}

/* Compute column or row sums, means or products */

gretl_matrix *gretl_cmatrix_vector_stat (const gretl_matrix *m,
					 GretlVecStat vs, int rowwise,
					 int *err)
{
    double complex z;
    gretl_matrix *ret;
    int r, c, i, j;

    if (!cmatrix_validate(m, 0)) {
	*err = E_INVARG;
	return NULL;
    }

    r = rowwise ? m->rows : 1;
    c = rowwise ? 1 : m->cols;

    ret = gretl_cmatrix_new(r, c);
    if (ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (rowwise) {
	/* by rows */
	int jmin = vs == V_PROD ? 1 : 0;

	for (i=0; i<m->rows; i++) {
	    z = vs == V_PROD ? m->z[0] : 0;
	    for (j=jmin; j<m->cols; j++) {
		if (vs == V_PROD) {
		    z *= gretl_cmatrix_get(m, i, j);
		} else {
		    z += gretl_cmatrix_get(m, i, j);
		}
	    }
	    if (vs == V_MEAN) {
		z /= m->cols;
	    }
	    gretl_cmatrix_set(ret, i, 0, z);
	}
    } else {
	/* by columns */
	int imin = vs == V_PROD ? 1 : 0;

	for (j=0; j<m->cols; j++) {
	    z = vs == V_PROD ? m->z[0] : 0;
	    for (i=imin; i<m->rows; i++) {
		if (vs == V_PROD) {
		    z *= gretl_cmatrix_get(m, i, j);
		} else {
		    z += gretl_cmatrix_get(m, i, j);
		}
	    }
	    if (vs == V_MEAN) {
		z /= m->rows;
	    }
	    gretl_cmatrix_set(ret, 0, j, z);
	}
    }

    return ret;
}

int gretl_cmatrix_fill (gretl_matrix *m, double complex z)
{
    if (!cmatrix_validate(m, 0)) {
	return E_TYPES;
    } else {
	int i, n = m->cols * m->rows;

	for (i=0; i<n; i++) {
	    m->z[i] = z;
	}
	return 0;
    }
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

static gretl_matrix *cmatrix_dot_op (const gretl_matrix *A,
				     const gretl_matrix *B,
				     int op, int *err)
{
    gretl_matrix *C = NULL;
    double complex x, y;
    int nr, nc;
    int conftype;
    int i, j, off;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    conftype = dot_operator_conf(A, B, &nr, &nc);

    if (conftype == CONF_NONE) {
	fputs("gretl_cmatrix_dot_op: matrices not conformable\n", stderr);
	fprintf(stderr, " op = '%c', A is %d x %d, B is %d x %d\n",
		(char) op, A->rows/2, A->cols, B->rows/2, B->cols);
	*err = E_NONCONF;
	return NULL;
    }

    C = gretl_cmatrix_new(nr, nc);
    if (C == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

#if 0 /* maybe expose this in gretl_matrix.h? */
    math_err_init();
#endif

    switch (conftype) {
    case CONF_ELEMENTS:
	vec_x_op_vec_y(C->z, A->z, B->z, nr*nc, op);
	break;
    case CONF_A_COLVEC:
	for (i=0; i<nr; i++) {
	    x = A->z[i];
	    for (j=0; j<nc; j++) {
		y = gretl_cmatrix_get(B, i, j);
		y = x_op_y(x, y, op);
		gretl_cmatrix_set(C, i, j, y);
	    }
	}
	break;
    case CONF_B_COLVEC:
	for (i=0; i<nr; i++) {
	    y = B->z[i];
	    for (j=0; j<nc; j++) {
		x = gretl_cmatrix_get(A, i, j);
		x = x_op_y(x, y, op);
		gretl_cmatrix_set(C, i, j, x);
	    }
	}
	break;
    case CONF_A_ROWVEC:
	off = 0;
	for (j=0; j<nc; j++) {
	    x = A->z[j];
	    x_op_vec_y(C->z + off, x, B->z + off, nr, op);
	    off += nr;
	}
	break;
    case CONF_B_ROWVEC:
	off = 0;
	for (j=0; j<nc; j++) {
	    y = B->z[j];
	    vec_x_op_y(C->z + off, A->z + off, y, nr, op);
	    off += nr;
	}
	break;
    case CONF_A_SCALAR:
	x_op_vec_y(C->z, A->z[0], B->z, nr*nc, op);
	break;
    case CONF_B_SCALAR:
	vec_x_op_y(C->z, A->z, B->z[0], nr*nc, op);
	break;
    case CONF_AC_BR:
	for (i=0; i<nr; i++) {
	    x = A->z[i];
	    for (j=0; j<nc; j++) {
		y = B->z[j];
		y = x_op_y(x, y, op);
		gretl_cmatrix_set(C, i, j, y);
	    }
	}
	break;
    case CONF_AR_BC:
	for (j=0; j<nc; j++) {
	    x = A->z[j];
	    for (i=0; i<nr; i++) {
		y = B->z[i];
		y = x_op_y(x, y, op);
		gretl_cmatrix_set(C, i, j, y);
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
	    gretl_matrix_free(C);
	    C = NULL;
	}
    }

    return C;
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

gretl_matrix *gretl_cmatrix_divide (const gretl_matrix *A,
				    const gretl_matrix *B,
				    GretlMatrixMod mod,
				    int *err)
{
    gretl_matrix *T = NULL;
    gretl_matrix *C = NULL;

    if (A->is_complex && B->is_complex) {
	C = gretl_matrix_divide(A, B, mod, err);
    } else if (A->is_complex) {
	/* case of real B */
	T = complex_from_real(B, err);
	if (T != NULL) {
	    C = gretl_matrix_divide(A, T, mod, err);
	}
    } else if (B->is_complex) {
	/* case of real A */
	T = complex_from_real(A, err);
	if (T != NULL) {
	    C = gretl_matrix_divide(T, B, mod, err);
	}
    } else {
	*err = E_TYPES;
    }

    gretl_matrix_free(T);

    return C;
}

gretl_matrix *cmatrix_get_element (const gretl_matrix *M,
				   int i, int *err)
{
    gretl_matrix *ret = NULL;

    if (M == NULL) {
	*err = E_DATA;
    } else if (i < 0 || i >= M->rows * M->cols) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	*err = E_INVARG;
    } else {
	ret = gretl_cmatrix_from_scalar(M->z[i], err);
    }

    return ret;
}

/* Set either the real or the imaginary part of @targ,
   using @src if non-NULL or @x otherwise.
*/

int gretl_cmatrix_set_part (gretl_matrix *targ,
			    const gretl_matrix *src,
			    double x, int im)
{
    double complex z;
    int i, j;

    if (!cmatrix_validate(targ, 0)) {
	return E_TYPES;
    } else if (src != NULL) {
	if (gretl_is_null_matrix(src) ||
	    src->rows != targ->rows ||
	    src->cols != targ->cols ||
	    src->is_complex) {
	    return E_INVARG;
	}
    }

    for (j=0; j<targ->cols; j++) {
	for (i=0; i<targ->rows; i++) {
	    if (src != NULL) {
		x = gretl_matrix_get(src, i, j);
	    }
	    z = gretl_cmatrix_get(targ, i, j);
	    if (im) {
		z = creal(z) + x * I;
	    } else {
		z = x + cimag(z) * I;
	    }
	    gretl_cmatrix_set(targ, i, j, z);
	}
    }

    return 0;
}
