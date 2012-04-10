#include "minpack.h"
#include <math.h>

int qrfac_(int m, int n, double *a, int lda, 
	   logical *pivot, int *ipvt, int lipvt, 
	   double *rdiag, double *acnorm, double *wa)
{
    /* Initialized data */
    const double p05 = .05;

    /* System generated locals */
    int a_dim1, a_offset, i2, i3;
    double d1, d2, d3;

    /* Local variables */
    int i, j, k, jp1;
    double sum;
    int kmax;
    double temp;
    int minmn;
    double epsmch;
    double ajnorm;

/*     subroutine qrfac */

/*     this subroutine uses householder transformations with column */
/*     pivoting (optional) to compute a qr factorization of the */
/*     m by n matrix a. that is, qrfac determines an orthogonal */
/*     matrix q, a permutation matrix p, and an upper trapezoidal */
/*     matrix r with diagonal elements of nonincreasing magnitude, */
/*     such that a*p = q*r. the householder transformation for */
/*     column k, k = 1,2,...,min(m,n), is of the form */

/*                           t */
/*           i - (1/u(k))*u*u */

/*     where u has zeros in the first k-1 positions. the form of */
/*     this transformation and the method of pivoting first */
/*     appeared in the corresponding linpack subroutine. */

/*     the subroutine statement is */

/*       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of a. */

/*       n is a positive integer input variable set to the number */
/*         of columns of a. */

/*       a is an m by n array. on input a contains the matrix for */
/*         which the qr factorization is to be computed. on output */
/*         the strict upper trapezoidal part of a contains the strict */
/*         upper trapezoidal part of r, and the lower trapezoidal */
/*         part of a contains a factored form of q (the non-trivial */
/*         elements of the u vectors described above). */

/*       lda is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array a. */

/*       pivot is a logical input variable. if pivot is set true, */
/*         then column pivoting is enforced. if pivot is set false, */
/*         then no column pivoting is done. */

/*       ipvt is an integer output array of length lipvt. ipvt */
/*         defines the permutation matrix p such that a*p = q*r. */
/*         column j of p is column ipvt(j) of the identity matrix. */
/*         if pivot is false, ipvt is not referenced. */

/*       lipvt is a positive integer input variable. if pivot is false, */
/*         then lipvt may be as small as 1. if pivot is true, then */
/*         lipvt must be at least n. */

/*       rdiag is an output array of length n which contains the */
/*         diagonal elements of r. */

/*       acnorm is an output array of length n which contains the */
/*         norms of the corresponding columns of the input matrix a. */
/*         if this information is not needed, then acnorm can coincide */
/*         with rdiag. */

/*       wa is a work array of length n. if pivot is false, then wa */
/*         can coincide with rdiag. */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

    /* Parameter adjustments */
    --wa;
    --acnorm;
    --rdiag;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */

    /* epsmch is the machine precision. */
    epsmch = dpmpar_(1);

    /* compute the initial column norms and initialize several arrays. */

    for (j = 1; j <= n; ++j) {
	acnorm[j] = enorm_(m, &a[j * a_dim1 + 1]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (*pivot) {
	    ipvt[j] = j;
	}
    }

    /* reduce a to r with householder transformations. */

    minmn = min(m,n);
    for (j = 1; j <= minmn; ++j) {
	if (! (*pivot)) {
	    goto L40;
	}

	/* bring the column of largest norm into the pivot position. */

	kmax = j;
	for (k = j; k <= n; ++k) {
	    if (rdiag[k] > rdiag[kmax]) {
		kmax = k;
	    }
	}
	if (kmax == j) {
	    goto L40;
	}
	for (i = 1; i <= m; ++i) {
	    temp = a[i + j * a_dim1];
	    a[i + j * a_dim1] = a[i + kmax * a_dim1];
	    a[i + kmax * a_dim1] = temp;
	}
	rdiag[kmax] = rdiag[j];
	wa[kmax] = wa[j];
	k = ipvt[j];
	ipvt[j] = ipvt[kmax];
	ipvt[kmax] = k;
L40:

/*        compute the householder transformation to reduce the */
/*        j-th column of a to a multiple of the j-th unit vector. */

	i2 = m - j + 1;
	ajnorm = enorm_(i2, &a[j + j * a_dim1]);
	if (ajnorm == 0.0) {
	    goto L100;
	}
	if (a[j + j * a_dim1] < 0.0) {
	    ajnorm = -ajnorm;
	}
	for (i = j; i <= m; ++i) {
	    a[i + j * a_dim1] /= ajnorm;
	}
	a[j + j * a_dim1] += 1.0;

/*        apply the transformation to the remaining columns */
/*        and update the norms. */

	jp1 = j + 1;
	if (n < jp1) {
	    goto L100;
	}
	for (k = jp1; k <= n; ++k) {
	    sum = 0.0;
	    for (i = j; i <= m; ++i) {
		sum += a[i + j * a_dim1] * a[i + k * a_dim1];
	    }
	    temp = sum / a[j + j * a_dim1];
	    for (i = j; i <= m; ++i) {
		a[i + k * a_dim1] -= temp * a[i + j * a_dim1];
	    }
	    if (!(*pivot) || rdiag[k] == 0.0) {
		goto L80;
	    }
	    temp = a[j + k * a_dim1] / rdiag[k];
	    d3 = temp;
	    d1 = 0.0, d2 = 1.0 - d3 * d3;
	    rdiag[k] *= sqrt((max(d1,d2)));
	    d1 = rdiag[k] / wa[k];
	    if (p05 * (d1 * d1) > epsmch) {
		goto L80;
	    }
	    i3 = m - j;
	    rdiag[k] = enorm_(i3, &a[jp1 + k * a_dim1]);
	    wa[k] = rdiag[k];
L80:
	    ;
	}
L100:
	rdiag[j] = -ajnorm;
    }

    return 0;
}

