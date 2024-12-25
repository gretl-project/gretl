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
#include "gretl_f2c.h"
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

#define cna(z) (na(creal(z)) || na(cimag(z)))

/* helper function for fftw-based real FFT functions */

static int fft_allocate (double **ffx, double complex **ffz,
                         gretl_matrix **ret, int r, int c,
                         int inverse)
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
    if (!inverse) {
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

static gretl_matrix *real_matrix_fft (const gretl_matrix *y,
                                      int inverse, int *err)
{
    gretl_matrix *ret = NULL;
    fftw_plan p = NULL;
    double *ffx = NULL;
    double complex *ffz = NULL;
    int r, c, m, odd;
    int i, j;

    if (y->rows < 2) {
        *err = E_DATA;
        return NULL;
    }

    r = y->rows;
    c = y->cols;
    m = r / 2;
    odd = r % 2;

    *err = fft_allocate(&ffx, &ffz, &ret, r, c, inverse);
    if (*err) {
        return NULL;
    }

    for (j=0; j<c; j++) {
        /* load the data */
        if (inverse) {
            for (i=0; i<=m+odd; i++) {
                ffz[i] = gretl_cmatrix_get(y, i, j);
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
        } else {
            for (i=0; i<=m+odd; i++) {
                gretl_cmatrix_set(ret, i, j, ffz[i]);
            }
            for (i=m; i>0; i--) {
                gretl_cmatrix_set(ret, r-i, j, conj(ffz[i]));
            }
        }
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

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err)
{
    return real_matrix_fft(y, 0, err);
}

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err)
{
    if (!y->is_complex) {
        *err = E_TYPES;
        return NULL;
    } else if (fft_is_hermitian(y)) {
        /* its inverse should be real */
        return real_matrix_fft(y, 1, err);
    } else {
        return gretl_cmatrix_fft(y, 1, err);
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

#if 0 /* currently unused but could be useful! */

/* Multiplication of complex matrices via BLAS zgemm(),
   allowing for conjugate transposition of @A or @B.
*/

static int gretl_zgemm_full (cmplx alpha,
                             const gretl_matrix *A,
                             char transa,
                             const gretl_matrix *B,
                             char transb,
                             cmplx beta,
                             gretl_matrix *C)
{
    integer lda, ldb, ldc;
    integer m, n, k, kb;

    if (transa == 'C') {
        m = A->cols; /* rows of op(A) */
        k = A->rows; /* cols of op(A) */
    } else {
        m = A->rows; /* complex rows of A */
        k = A->cols; /* cols of A */
    }

    if (transb == 'C') {
        n = B->rows;  /* columns of op(B) */
        kb = B->cols; /* rows of op(B) */
    } else {
        n = B->cols;
        kb = B->rows;
    }

    if (k != kb || C->rows != m || C->cols != n) {
        return E_NONCONF;
    }

    lda = A->rows;
    ldb = B->rows;
    ldc = C->rows;

    zgemm_(&transa, &transb, &m, &n, &k, &alpha,
           (cmplx *) A->val, &lda, (cmplx *) B->val, &ldb,
           &beta, (cmplx *) C->val, &ldc);

    return 0;
}

#endif /* currently unused */

/* Variant of zgemm: simplified version of gretl_zgemm_full
   which allocates the product matrix, C.
*/

static gretl_matrix *gretl_zgemm (const gretl_matrix *A,
                                  char transa,
                                  const gretl_matrix *B,
                                  char transb,
                                  int *err)
{
    gretl_matrix *C;
    cmplx alpha = {1, 0};
    cmplx beta = {0, 0};
    integer lda, ldb, ldc;
    integer m, n, k, kb;

    if (!cmatrix_validate(A, 0) || !cmatrix_validate(B, 0)) {
        *err = E_INVARG;
        return NULL;
    }

    if (transa == 'C') {
        m = A->cols; /* rows of op(A) */
        k = A->rows; /* cols of op(A) */
    } else {
        m = A->rows; /* complex rows of A */
        k = A->cols; /* cols of A */
    }

    if (transb == 'C') {
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

    lda = A->rows;
    ldb = B->rows;
    ldc = C->rows;

    zgemm_(&transa, &transb, &m, &n, &k, &alpha,
           (cmplx *) A->val, &lda, (cmplx *) B->val, &ldb,
           &beta, (cmplx *) C->val, &ldc);

    return C;
}

static gretl_matrix *real_cmatrix_multiply (const gretl_matrix *A,
                                            const gretl_matrix *B,
                                            char transa, int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *T = NULL;

    if (A->is_complex && B->is_complex) {
        if (transa == 'N' && (cscalar(A) || cscalar(B))) {
            return gretl_cmatrix_dot_op(A, B, '*', err);
        } else {
            C = gretl_zgemm(A, transa, B, 'N', err);
        }
    } else if (A->is_complex) {
        /* case of real B */
        T = complex_from_real(B, err);
        if (T != NULL) {
            C = gretl_zgemm(A, transa, T, 'N', err);
        }
    } else if (B->is_complex) {
        /* case of real A */
        T = complex_from_real(A, err);
        if (T != NULL) {
            C = gretl_zgemm(T, transa, B, 'N', err);
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
    return real_cmatrix_multiply(A, B, 'N', err);
}

/* Returns (conjugate transpose of A, or A^H) times B,
   allowing for the possibility that either A or B (but
   not both!) may be a real matrix on input.
*/

gretl_matrix *gretl_cmatrix_AHB (const gretl_matrix *A,
                                 const gretl_matrix *B,
                                 int *err)
{
    return real_cmatrix_multiply(A, B, 'C', err);
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

static int ensure_aux_cmatrix (gretl_matrix *m,
                               gretl_matrix **paux,
                               cmplx **ppz,
                               int r, int c)
{
    gretl_matrix *aux = NULL;
    int mrc, dim = r * c;

    mrc = (m == NULL)? 0 : m->rows * m->cols;

    /* We need an r x c complex matrix for output:
       can we just use @m?
    */
    if (mrc == dim && m->is_complex) {
        /* OK, reusable complex matrix */
        m->rows = r;
        m->cols = c;
        if (ppz != NULL) {
            *ppz = (cmplx *) m->val;
        }
    } else if (mrc == 2*dim && !m->is_complex) {
        /* OK, reusable real matrix */
        m->rows = 2*r;
        m->cols = c;
        gretl_matrix_set_complex_full(m, 1);
        if (ppz != NULL) {
            *ppz = (cmplx *) m->val;
        }
    } else {
        /* nope, have to allocate anew */
        aux = gretl_cmatrix_new0(r, c);
        if (aux == NULL) {
            return E_ALLOC;
        } else {
            *paux = aux;
            if (ppz != NULL) {
                *ppz = (cmplx *) aux->val;
            }
        }
    }

    return 0;
}

static int ensure_aux_matrix (gretl_matrix *m,
                              gretl_matrix **paux,
                              int r, int c)
{
    gretl_matrix *aux = NULL;
    int mrc, dim = r * c;

    mrc = (m == NULL)? 0 : m->rows * m->cols;

    if (mrc == dim && !m->is_complex) {
        /* OK, reusable real matrix */
        m->rows = r;
        m->cols = c;
    } else {
        /* have to allocate anew */
        aux = gretl_matrix_alloc(r, c);
        if (aux == NULL) {
            return E_ALLOC;
        } else {
            *paux = aux;
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
        *err = ensure_aux_cmatrix(VL, &Ltmp, &vl, n, n);
        if (*err) {
            goto bailout;
        }
    }

    if (VR != NULL) {
        /* right eigenvectors wanted */
        *err = ensure_aux_cmatrix(VR, &Rtmp, &vr, n, n);
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

/* Schur factorization of @A, with optional assignment of the
   matrix of Schur vectors to @Z and/or the eigenvalues of @A
   to @W, via LAPACK zgees().
*/

gretl_matrix *gretl_zgees (const gretl_matrix *A,
                           gretl_matrix *Z,
                           gretl_matrix *W,
                           int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *Ztmp = NULL;
    gretl_matrix *Wtmp = NULL;
    integer n, info, lwork;
    integer ldvs;
    double *rwork = NULL;
    cmplx *a = NULL;
    cmplx *work = NULL;
    cmplx *vs = NULL;
    cmplx *w = NULL;
    char jobvs = Z != NULL ? 'V' : 'N';
    char srt = 'N';
    integer sdim = 0;

    if (!cmatrix_validate(A, 1)) {
        *err = E_INVARG;
        return NULL;
    }

    n = A->rows;
    ldvs = Z != NULL ? n : 1;

    /* we need a copy of @A, which gets overwritten */
    ret = gretl_matrix_copy(A);
    if (ret == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    a = (cmplx *) ret->z;

    if (Z != NULL) {
        /* Schur vectors wanted */
        *err = ensure_aux_cmatrix(Z, &Ztmp, &vs, n, n);
        if (*err) {
            goto bailout;
        }
    }

    /* Eigenvalues vector: seems like we need this,
       regardless of whether @W is passed?
    */
    *err = ensure_aux_cmatrix(W, &Wtmp, &w, n, 1);
    if (*err) {
        goto bailout;
    }

    work = malloc(sizeof *work);
    rwork = malloc(n * sizeof *rwork);
    if (work == NULL || rwork == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    /* get optimal workspace size */
    lwork = -1;
    zgees_(&jobvs, &srt, NULL, &n, a, &n, &sdim, w, vs, &ldvs,
           work, &lwork, rwork, NULL, &info);
    lwork = (integer) work[0].r;
    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    /* do the actual factorization */
    zgees_(&jobvs, &srt, NULL, &n, a, &n, &sdim, w, vs, &ldvs,
           work, &lwork, rwork, NULL, &info);
    if (info != 0) {
        *err = E_DATA;
    } else {
        if (Ztmp != NULL) {
            gretl_matrix_replace_content(Z, Ztmp);
        }
        if (W != NULL && Wtmp != NULL) {
            gretl_matrix_replace_content(W, Wtmp);
        }
    }

 bailout:

    free(rwork);
    free(work);
    gretl_matrix_free(Ztmp);
    gretl_matrix_free(Wtmp);

    if (*err) {
        gretl_matrix_free(ret);
        ret = NULL;
    }

    return ret;
}

/* SVD of a complex matrix via the LAPACK function zgesvd() */

int gretl_cmatrix_SVD (const gretl_matrix *x, gretl_matrix **pu,
                       gretl_vector **ps, gretl_matrix **pvt,
                       int full)
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
        if (full) {
            u = gretl_cmatrix_new(ldu, m);
        } else {
            u = gretl_cmatrix_new(ldu, k);
        }
        if (u == NULL) {
            err = E_ALLOC;
            goto bailout;
        } else {
            uval = (cmplx *) u->val;
            jobu = full ? 'A' : 'S';
        }
    }

    if (pvt != NULL) {
        ldvt = full ? n : k;
        vt = gretl_cmatrix_new(ldvt, n);
        if (vt == NULL) {
            err = E_ALLOC;
            goto bailout;
        } else {
            vtval = (cmplx *) vt->val;
            jobvt = full ? 'A' : 'S';
        }
    }

    work = malloc(sizeof *work);
    rwork = malloc(5 * k * sizeof *rwork);
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
        char trans1 = A->rows > k ? 'C' : 'N';
        char trans2 = A->cols > k ? 'C' : 'N';
        gretl_matrix *B;

        B = gretl_zgemm(A, trans1, A, trans2, err);
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

/* Inverse of a complex matrix via its SVD. If the  @generalized
   flag is non-zero the generalized inverse is produced in case
   of rank deficiency, othewise an error is flagged if @A is not
   of full rank.
*/

static gretl_matrix *cmatrix_SVD_inverse (const gretl_matrix *A,
                                          int generalized,
                                          int *err)
{
    gretl_matrix *U = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *s = NULL;
    gretl_matrix *Vt = NULL;
    gretl_matrix *ret = NULL;

    *err = gretl_cmatrix_SVD(A, &U, &s, &V, 1);

    if (!*err) {
        double smin = csvd_smin(A, s->val[0]);
        double complex vij;
        int i, j, h = 0;

        for (i=0; i<s->cols; i++) {
            h += s->val[i] > smin;
        }
        if (!generalized && h < s->cols) {
            gretl_errmsg_set(_("Matrix is singular"));
            *err = E_SINGULAR;
            goto bailout;
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
            ret = gretl_zgemm(Vt, 'N', U, 'C', err);
        }
    }

 bailout:

    gretl_matrix_free(U);
    gretl_matrix_free(V);
    gretl_matrix_free(s);
    gretl_matrix_free(Vt);

    return ret;
}

gretl_matrix *gretl_cmatrix_ginv (const gretl_matrix *A, int *err)
{
    return cmatrix_SVD_inverse(A, 1, err);
}

gretl_matrix *gretl_cmatrix_inverse (const gretl_matrix *A, int *err)
{
    return cmatrix_SVD_inverse(A, 0, err);
}

/* Horizontal direct product of complex @A and @B, or of
   @A with itself if B is NULL.
*/

static gretl_matrix *real_cmatrix_hdp (const gretl_matrix *A,
                                       const gretl_matrix *B,
                                       int *err)
{
    gretl_matrix *C = NULL;
    double complex aij, bik;
    int do_symmetric = 0;
    int i, j, k, ndx;
    int r, p, q, ccols;

    if (!cmatrix_validate(A,0)) {
        *err = E_INVARG;
    } else if (B != NULL) {
        if (!cmatrix_validate(B,0)) {
            *err = E_INVARG;
        } else if (B->rows != A->rows) {
            *err = E_NONCONF;
        }
    } else {
        do_symmetric = 1;
    }

    if (*err) {
        return NULL;
    }

    r = A->rows;
    p = A->cols;

    if (do_symmetric) {
        q = p;
        ccols = p * (p+1) / 2;
    } else {
        q = B->cols;
        ccols = p * q;
    }

    C = gretl_cmatrix_new0(r, ccols);

    if (C == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (i=0; i<r; i++) {
        ndx = 0;
        for (j=0; j<p; j++) {
            aij = gretl_cmatrix_get(A, i, j);
            if (do_symmetric) {
                for (k=j; k<q; k++) {
                    bik = gretl_cmatrix_get(A, i, k);
                    gretl_cmatrix_set(C, i, ndx++, aij*conj(bik));
                }
            } else if (aij != 0.0) {
                ndx = j * q;
                for (k=0; k<q; k++) {
                    bik = gretl_cmatrix_get(B, i, k);
                    gretl_cmatrix_set(C, i, ndx + k, aij*bik);
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

    if (A->is_complex && hdp && B == NULL) {
        ; /* OK */
    } else if (A->is_complex && B->is_complex) {
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
   components. If @Re is NULL the result will have a
   constant real part given by @x; and if @Im is NULL it
   will have a constant imaginary part given by @y. If
   both matrix arguments are non-NULL they must be of
   the same dimensions.
*/

gretl_matrix *gretl_cmatrix_build (const gretl_matrix *Re,
                                   const gretl_matrix *Im,
                                   double x, double y,
                                   int *err)
{
    gretl_matrix *C = NULL;
    int r = 1, c = 1;
    int i, n;

    if (Re != NULL) {
        if (Re->is_complex) {
            *err = E_INVARG;
        } else {
            r = Re->rows;
            c = Re->cols;
        }
    }
    if (!*err && Im != NULL) {
        if (Im->is_complex) {
            *err = E_INVARG;
        } else if (Re != NULL) {
            if (Im->rows != r || Im->cols != c) {
                *err = E_NONCONF;
            }
        } else {
            r = Im->rows;
            c = Im->cols;
        }
    }

    if (!*err) {
        n = r * c;
        C = gretl_cmatrix_new(r, c);
        if (C == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        for (i=0; i<n; i++) {
            if (Re != NULL) {
                x = Re->val[i];
            }
            if (Im != NULL) {
                y = Im->val[i];
            }
            C->z[i] = x + y * I;
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

double gretl_cquad (double complex z)
{
    double complex zr = creal(z);
    double complex zi = cimag(z);

    return zr*zr + zi*zi;
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
    int err = 0;

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
                                         int ldet, int *err)
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
        if (ldet) {
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
    int lower = !upper;
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

    if (r > c && lower) {
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
            if (lower) {
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
                                         int skip_na, int *err)
{
    double complex z;
    double complex mij;
    gretl_matrix *ret;
    int product;
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

    product = (vs == V_PROD);

    if (rowwise) {
        /* by rows */
        double complex z0 = product ? 1 : 0;
        int ncols;

        for (i=0; i<m->rows; i++) {
            ncols = m->cols;
            z = z0;
            for (j=0; j<m->cols; j++) {
                mij = gretl_cmatrix_get(m, i, j);
                if (cna(mij) && skip_na) {
                    ncols--;
                    continue;
                } else if (product) {
                    z *= mij;
                } else {
                    z += mij;
                }
            }
            if (ncols == 0) {
                z = NADBL + 0*I;
            } else if (vs == V_MEAN) {
                z /= ncols;
            }
            gretl_cmatrix_set(ret, i, 0, z);
        }
    } else {
        /* by columns */
        double complex z0 = product ? 1 : 0;
        int nrows;

        for (j=0; j<m->cols; j++) {
            nrows = m->rows;
            z = z0;
            for (i=0; i<m->rows; i++) {
                mij = gretl_cmatrix_get(m, i, j);
                if (cna(mij) && skip_na) {
                    nrows--;
                    continue;
                } else if (product) {
                    z *= mij;
                } else {
                    z += mij;
                }
            }
            if (nrows == 0) {
                z = NADBL + 0*I;
            } else if (vs == V_MEAN) {
                z /= nrows;
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

static int imag_part_zero (const gretl_matrix *A)
{
    int i, n = A->rows * A->cols;
    double tol = 1.0e-15; /* ?? */

    for (i=0; i<n; i++) {
        if (fabs(cimag(A->z[i])) > tol) {
            fprintf(stderr, "imag_part_zero? no, got %g\n",
                    cimag(A->z[i]));
            return 0;
        }
    }
    return 1;
}

/* Matrix logarithm, for a diagonalizable matrix */

gretl_matrix *gretl_matrix_log (const gretl_matrix *A, int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *Vi = NULL;
    gretl_matrix *w = NULL;
    gretl_matrix *R = NULL;
    int i;

    if (gretl_is_null_matrix(A) || A->rows != A->cols) {
        /* we require a square matrix, real or complex */
        *err = E_INVARG;
        return NULL;
    }

    if (A->is_complex) {
        C = gretl_matrix_copy(A);
        if (C == NULL) {
            *err = E_ALLOC;
        }
    } else {
        C = gretl_cmatrix_build(A, NULL, 0, 0, err);
    }

    if (!*err) {
        V = gretl_cmatrix_new(A->rows, A->rows);
        if (V == NULL) {
            *err = E_ALLOC;
        } else {
            w = gretl_zgeev(C, V, NULL, err);
        }
    }

    if (!*err) {
        /* The following step will fail if
           @A is not diagonalizable
        */
        Vi = gretl_cmatrix_inverse(V, err);
    }

    if (!*err) {
        for (i=0; i<A->rows; i++) {
            w->z[i] = clog(w->z[i]);
        }
        R = cmatrix_dot_op(w, Vi, '*', err);
        if (!*err) {
            gretl_matrix_free(Vi);
            Vi = gretl_cmatrix_multiply(V, R, err);
        }
    }

    gretl_matrix_free(C);

    if (!*err) {
        if (imag_part_zero(Vi)) {
            /* should we do this? */
            C = gretl_cmatrix_extract(Vi, 0, err);
        } else {
            C = Vi;
            Vi = NULL;
        }
    } else {
        C = NULL;
    }

    gretl_matrix_free(V);
    gretl_matrix_free(Vi);
    gretl_matrix_free(w);
    gretl_matrix_free(R);

    return C;
}

/* Matrix exponential for a diagonalizable complex matrix */

gretl_matrix *gretl_cmatrix_exp (const gretl_matrix *A, int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *Vi = NULL;
    gretl_matrix *w = NULL;
    gretl_matrix *R = NULL;
    int i;

    if (!cmatrix_validate(A, 1)) {
        *err = E_INVARG;
        return NULL;
    }

    C = gretl_matrix_copy(A);
    if (C == NULL) {
        *err = E_ALLOC;
    }

    if (!*err) {
        V = gretl_cmatrix_new(A->rows, A->rows);
        if (V == NULL) {
            *err = E_ALLOC;
        } else {
            w = gretl_zgeev(C, V, NULL, err);
        }
    }

    if (!*err) {
        /* The following step will fail if
           @A is not diagonalizable
        */
        Vi = gretl_cmatrix_inverse(V, err);
    }

    gretl_matrix_free(C);
    C = NULL;

    if (!*err) {
        for (i=0; i<A->rows; i++) {
            w->z[i] = cexp(w->z[i]);
        }
        R = cmatrix_dot_op(w, Vi, '*', err);
        if (!*err) {
            C = gretl_cmatrix_multiply(V, R, err);
        }
    }

    gretl_matrix_free(V);
    gretl_matrix_free(Vi);
    gretl_matrix_free(w);
    gretl_matrix_free(R);

    return C;
}

static int matrix_is_hermitian (const gretl_matrix *z)
{
    double complex zij, zji;
    int i, j, n = z->rows;

    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            zij = gretl_cmatrix_get(z, i, j);
            zji = gretl_cmatrix_get(z, j, i);
            if (zji != conj(zij)) {
                return 0;
            }
        }
    }
    return 1;
}

/* Cholesky decomposition of Hermitian matrix using LAPACK
   function zpotrf()
*/

gretl_matrix *gretl_cmatrix_cholesky (const gretl_matrix *A,
                                      int *err)
{
    gretl_matrix *C = NULL;
    char uplo = 'L';
    integer n, info;

    if (!cmatrix_validate(A, 1)) {
        *err = E_INVARG;
    } else if (!matrix_is_hermitian(A)) {
        gretl_errmsg_set(_("Matrix is not Hermitian"));
        *err = E_INVARG;
    } else {
        C = gretl_matrix_copy(A);
        if (C == NULL) {
            *err = E_ALLOC;
        }
    }

    if (*err) {
        return NULL;
    }

    n = A->rows;
    zpotrf_(&uplo, &n, (cmplx *) C->val, &n, &info);
    if (info != 0) {
        fprintf(stderr, "gretl_cmatrix_cholesky: "
                "zpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
        *err = (info > 0)? E_NOTPD : E_DATA;
    }

    if (*err) {
        gretl_matrix_free(C);
        C = NULL;
    } else {
        gretl_matrix_zero_upper(C);
    }

    return C;
}

gretl_matrix *gretl_cmatrix_QR_pivot_decomp (const gretl_matrix *A,
                                             gretl_matrix *R,
                                             gretl_matrix *P,
                                             int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *Rtmp = NULL;
    gretl_matrix *Ptmp = NULL;
    integer m, n, lda;
    integer info = 0;
    integer lwork = -1;
    integer *jpvt = NULL;
    cmplx *tau = NULL;
    cmplx *work = NULL;
    double *rwork = NULL;
    int i, j;

    if (!cmatrix_validate(A, 0)) {
        *err = E_INVARG;
        return NULL;
    }

    lda = m = A->rows;
    n = A->cols;

    if (n > m) {
        gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns"));
        *err = E_NONCONF;
        return NULL;
    }

    Q = gretl_matrix_copy(A);
    if (Q == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    if (R != NULL) {
        *err = ensure_aux_cmatrix(R, &Rtmp, NULL, n, n);
        if (*err) {
            goto bailout;
        } else if (Rtmp != NULL) {
            gretl_matrix_replace_content(R, Rtmp);
        }
    }
    if (P != NULL) {
        *err = ensure_aux_matrix(P, &Ptmp, 1, n);
        if (*err) {
            goto bailout;
        } else if (Ptmp != NULL) {
            gretl_matrix_replace_content(P, Ptmp);
        }
    }

    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = malloc(sizeof *work);
    rwork = malloc(2 * n * sizeof *rwork);
    if (tau == NULL || work == NULL || rwork == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    /* pivot array */
    jpvt = malloc(n * sizeof *jpvt);
    if (jpvt == NULL) {
        *err = E_ALLOC;
        goto bailout;
    } else {
        for (i=0; i<n; i++) {
            jpvt[i] = 0;
        }
    }

    /* workspace size query */
    zgeqp3_(&m, &n, (cmplx *) Q->val, &lda, jpvt, tau, work, &lwork, rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqp3: info = %d\n", (int) info);
        *err = E_DATA;
    } else {
        /* optimally sized work array */
        lwork = (integer) work[0].r;
        work = realloc(work, (size_t) lwork * sizeof *work);
        if (work == NULL) {
            *err = E_ALLOC;
        }
    }

    if (*err) {
        goto bailout;
    }

    /* run actual QR factorization */
    zgeqp3_(&m, &n, (cmplx *) Q->val, &lda, jpvt, tau, work, &lwork, rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqp3: info = %d\n", (int) info);
        *err = E_DATA;
        goto bailout;
    }

    if (R != NULL) {
        /* copy the upper triangular R out of Q */
        double complex z;

        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (i <= j) {
                    z = gretl_cmatrix_get(Q, i, j);
                    gretl_cmatrix_set(R, i, j, z);
                } else {
                    gretl_cmatrix_set(R, i, j, 0.0);
                }
            }
        }
    }

    /* turn Q into "the real" Q, with the help of tau */
    zungqr_(&m, &n, &n, (cmplx *) Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zungqr: info = %d\n", (int) info);
        *err = E_DATA;
    }

    if (P != NULL) {
        for (i=0; i<n; i++) {
            P->val[i] = jpvt[i];
        }
    }

 bailout:

    free(tau);
    free(work);
    free(rwork);
    free(jpvt);
    gretl_matrix_free(Rtmp);
    gretl_matrix_free(Ptmp);

    if (*err) {
        gretl_matrix_free(Q);
        Q = NULL;
    }

    return Q;
}

/* QR decomposition of complex matrix using LAPACK functions
   zgeqrf() and zungqr().
*/

gretl_matrix *gretl_cmatrix_QR_decomp (const gretl_matrix *A,
                                       gretl_matrix *R,
                                       int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *Rtmp = NULL;
    integer m, n, lda;
    integer info = 0;
    integer lwork = -1;
    cmplx *tau = NULL;
    cmplx *work = NULL;
    int i, j;

    if (!cmatrix_validate(A, 0)) {
        *err = E_INVARG;
        return NULL;
    }

    lda = m = A->rows;
    n = A->cols;

    if (n > m) {
        gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns"));
        *err = E_NONCONF;
        return NULL;
    }

    Q = gretl_matrix_copy(A);
    if (Q == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    if (R != NULL) {
        *err = ensure_aux_cmatrix(R, &Rtmp, NULL, n, n);
        if (*err) {
            goto bailout;
        } else if (Rtmp != NULL) {
            gretl_matrix_replace_content(R, Rtmp);
        }
    }

    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = malloc(sizeof *work);
    if (tau == NULL || work == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    /* workspace size query */
    zgeqrf_(&m, &n, (cmplx *) Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqrf: info = %d\n", (int) info);
        *err = E_DATA;
    } else {
        /* optimally sized work array */
        lwork = (integer) work[0].r;
        work = realloc(work, (size_t) lwork * sizeof *work);
        if (work == NULL) {
            *err = E_ALLOC;
        }
    }

    if (*err) {
        goto bailout;
    }

    /* run actual QR factorization */
    zgeqrf_(&m, &n, (cmplx *) Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqrf: info = %d\n", (int) info);
        *err = E_DATA;
        goto bailout;
    }

    if (R != NULL) {
        /* copy the upper triangular R out of Q */
        double complex z;

        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (i <= j) {
                    z = gretl_cmatrix_get(Q, i, j);
                    gretl_cmatrix_set(R, i, j, z);
                } else {
                    gretl_cmatrix_set(R, i, j, 0.0);
                }
            }
        }
    }

    /* turn Q into "the real" Q, with the help of tau */
    zungqr_(&m, &n, &n, (cmplx *) Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zungqr: info = %d\n", (int) info);
        *err = E_DATA;
    }

 bailout:

    free(tau);
    free(work);
    gretl_matrix_free(Rtmp);

    if (*err) {
        gretl_matrix_free(Q);
        Q = NULL;
    }

    return Q;
}

/* Copy data from real matrix @src to complex matrix *targ,
   starting at offset (@r0, @c0) into @targ. @src is treated
   as a complex matrix with zero imaginary part.
*/

void real_to_complex_fill (gretl_matrix *targ,
                           const gretl_matrix *src,
                           int r0, int c0)
{
    double complex z;
    int i, j, r, c;

    for (j=0, c=c0; j<src->cols && c<targ->cols; j++, c++) {
        for (i=0, r=r0; i<src->rows && r<targ->rows; i++, r++) {
            z = gretl_matrix_get(src, i, j);
            gretl_cmatrix_set(targ, r, c, z);
        }
    }
}

/* Tests whether a matrix has the is_complex flag set, and also
   whether it's "really complex" (has a non-zero imaginary part).
*/

int matrix_is_complex (const gretl_matrix *M)
{
    if (M == NULL || !M->is_complex) {
        return 0;
    } else {
        int i, n = M->rows * M->cols;
        int ret = 1;

        for (i=0; i<n; i++) {
            if (cimag(M->z[i]) != 0) {
                ret = 2;
                break;
            }
        }
        return ret;
    }
}
