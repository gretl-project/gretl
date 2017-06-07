#ifndef CLAPACK_COMPLEX_H
#define CLAPACK_COMPLEX_H

/* LAPACK subroutines: double-precision complex versions only */

void zheev_ (const char *jobz, const char *uplo, integer *n,
	     cmplx *a, integer *lda, double *w, cmplx *work,
	     integer *lwork, double *rwork, integer *info);

void zgetrf_ (integer *m, integer *n, cmplx *a, integer *lda,
	      integer *ipiv, integer *info);

void zgetri_ (integer *n, cmplx *a, integer *lda, integer *ipiv,
	      cmplx *work, integer *lwork, integer *info);

#endif /* CLAPACK_COMPLEX_H */
