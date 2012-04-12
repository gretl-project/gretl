/* 
   This source file based on Minpack: initially converted from 
   fortran using f2c, then rendered into relatively idiomatic
   C with zero-based indexing throughout and pass-by-value for
   parameters that do not function as pointers. We also rely
   on <float.h> for the machine precision rather than Minpack's
   dpmpar().

   See README in this directory for the Minpack Copyright.

   Allin Cottrell, Wake Forest University, April 2012
*/

#include "minpack.h"
#include <math.h>

/*
c     qrsolv:
c
c     given an m by n matrix a, an n by n diagonal matrix d,
c     and an m-vector b, the problem is to determine an x which
c     solves the system
c
c           a*x = b ,     d*x = 0 ,
c
c     in the least squares sense.
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, if
c     a*p = q*r, where p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then qrsolv expects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q transpose)*b. the system
c     a*x = b, d*x = 0, is then equivalent to
c
c                  t       t
c           r*z = q *b ,  p *d*p*z = 0 ,
c
c     where x = p*z. if this system does not have full rank,
c     then a least squares solution is obtained. on output qrsolv
c     also provides an upper triangular matrix s such that
c
c            t   t               t
c           p *(a *a + d*d)*p = s *s .
c
c     s is computed within qrsolv and may be of separate interest.
c
c     the subroutine statement is
c
c       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an n by n array. on input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         on output the full upper triangle is unaltered, and the
c         strict lower triangle contains the strict upper triangle
c         (transposed) of the upper triangular matrix s.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       ipvt is an integer input array of length n which defines the
c         permutation matrix p such that a*p = q*r. column j of p
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q transpose)*b.
c
c       x is an output array of length n which contains the least
c         squares solution of the system a*x = b, d*x = 0.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix s.
c
c       wa is a work array of length n.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int qrsolv_(int n, double *r, int ldr, int *ipvt, 
	    double *diag, double *qtb, double *x, 
	    double *sdiag, double *wa)
{
    const double p5 = .5;
    const double p25 = .25;
    double tanx, cosx, sinx, cotan;
    double sum, temp, qtbpj;
    int i, j, k, l;
    int nsing;

    /* copy r and Q'*b to preserve input and initialize s;
       in particular, save the diagonal elements of r in x 
    */

    for (j = 0; j < n; ++j) {
	for (i = j; i < n; ++i) {
	    r[i + j * ldr] = r[j + i * ldr];
	}
	x[j] = r[j + j * ldr];
	wa[j] = qtb[j];
    }

    /* eliminate the diagonal matrix d using a Givens rotation */

    for (j = 0; j < n; ++j) {
	/* prepare the row of d to be eliminated, locating the
	   diagonal element using p from the QR factorization 
	*/
	l = ipvt[j];
	if (diag[l] == 0.0) {
	    goto store;
	}
	for (k = j; k < n; ++k) {
	    sdiag[k] = 0.0;
	}
	sdiag[j] = diag[l];

	/* the transformations to eliminate the row of d
	   modify only a single element of Q'*b beyond the 
	   first n, which is initially zero 
	*/

	qtbpj = 0.0;
	for (k = j; k < n; ++k) {
	    /* determine a Givens rotation which eliminates the
	       appropriate element in the current row of d
	    */
	    if (sdiag[k] == 0.0) {
		continue;
	    }
	    if (fabs(r[k + k * ldr]) < fabs(sdiag[k])) {
		cotan = r[k + k * ldr] / sdiag[k];
		sinx = p5 / sqrt(p25 + p25 * (cotan * cotan));
		cosx = sinx * cotan;
	    } else {
		tanx = sdiag[k] / r[k + k * ldr];
		cosx = p5 / sqrt(p25 + p25 * (tanx * tanx));
		sinx = cosx * tanx;
	    }

	    /* compute the modified diagonal element of r and
	       the modified element of (Q'*b,0).
	    */
	    r[k + k * ldr] = cosx * r[k + k * ldr] + sinx * sdiag[k];
	    temp = cosx * wa[k] + sinx * qtbpj;
	    qtbpj = -sinx * wa[k] + cosx * qtbpj;
	    wa[k] = temp;

	    /* accumulate the tranformation in the row of s */
	    for (i = k+1; i < n; ++i) {
		temp = cosx * r[i + k * ldr] + sinx * sdiag[i];
		sdiag[i] = -sinx * r[i + k * ldr] + cosx * sdiag[i];
		r[i + k * ldr] = temp;
	    }
	}

    store:
	/* store the diagonal element of s and restore
	   the corresponding diagonal element of r
	*/
	sdiag[j] = r[j + j * ldr];
	r[j + j * ldr] = x[j];
    }

    /* solve the triangular system for z; if the system is
       singular, then obtain a least squares solution 
    */

    nsing = n;
    for (j = 0; j < n; ++j) {
	if (sdiag[j] == 0.0 && nsing == n) {
	    nsing = j;
	}
	if (nsing < n) {
	    wa[j] = 0.0;
	}
    }
    if (nsing > 0) {
	for (k = 0; k < nsing; ++k) {
	    j = nsing - k - 1;
	    sum = 0.0;
	    for (i = j+1; i < nsing; ++i) {
		sum += r[i + j * ldr] * wa[i];
	    }
	    wa[j] = (wa[j] - sum) / sdiag[j];
	}
    }

    /* permute the components of z back to components of x */

    for (j = 0; j < n; ++j) {
	l = ipvt[j];
	x[l] = wa[j];
    }

    return 0;
}

