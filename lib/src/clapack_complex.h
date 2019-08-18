#ifndef CLAPACK_COMPLEX_H
#define CLAPACK_COMPLEX_H

/* LAPACK subroutines: double-precision complex versions only */

#ifndef integer
# include <stdint.h>
typedef int32_t integer;
#endif

void zheev_ (const char *jobz, const char *uplo, integer *n,
	     cmplx *a, integer *lda, double *w, cmplx *work,
	     integer *lwork, double *rwork, integer *info);

void zgeev_ (const char *jobvl, const char *jobvr, integer *n,
	     cmplx *a, integer *lda, double *w, cmplx *vl,
	     integer *ldvl, cmplx *vr, integer *ldvr, cmplx *work,
	     integer *lwork, double *rwork, integer *info);

void zgelsy_(integer *m, integer *n, integer *nrhs, cmplx *a, integer *lda,
	     cmplx *b, integer *ldb, integer *jpvt, double *rcond, integer *rank,
	     cmplx *work, integer *lwork, double *rwork, integer *info);

void zgetrs_ (char *trans, integer *n, integer *nrhs, cmplx *a, integer *lda,
	      integer *ipiv, cmplx *b, integer *ldb, integer *info);

void zgetrf_ (integer *m, integer *n, cmplx *a, integer *lda,
	      integer *ipiv, integer *info);

void zgetri_ (integer *n, cmplx *a, integer *lda, integer *ipiv,
	      cmplx *work, integer *lwork, integer *info);

void zgemm_ (const char *transa, const char *transb,
	     integer *m, integer *n, integer *k,
	     cmplx *alpha, cmplx *a, integer *lda,
	     cmplx *b, integer *ldb, cmplx *beta,
	     cmplx *c, integer *ldc);

void zsyrk_ (const char *uplo, const char *trans, integer *n,
	     integer *k, cmplx *alpha, cmplx *a, integer *lda,
	     cmplx *beta, cmplx *c, integer *ldc);

void zgesvd_ (const char *jobu, const char *jobvt,
	      integer *m, integer *n, cmplx *a, integer *lda,
	      double *s, cmplx *u, integer *ldu, cmplx *vt, integer *ldvt,
	      cmplx *work, integer *lwork, double *rwork, integer *info);

#endif /* CLAPACK_COMPLEX_H */
