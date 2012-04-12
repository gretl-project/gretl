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
#include <float.h>

/*
c     lmpar:
c
c     given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta,
c     the problem is to determine a value for the parameter
c     par such that if x solves the system
c
c           a*x = b ,     sqrt(par)*d*x = 0 ,
c
c     in the least squares sense, and dxnorm is the euclidean
c     norm of d*x, then either par is zero and
c
c           (dxnorm-delta) .le. 0.1*delta ,
c
c     or par is positive and
c
c           abs(dxnorm-delta) .le. 0.1*delta .
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, if
c     a*p = q*r, where p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, then lmpar expects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q transpose)*b. on output
c     lmpar also provides an upper triangular matrix s such that
c
c            t   t                   t
c           p *(a *a + par*d*d)*p = s *s .
c
c     s is employed within lmpar and may be of separate interest.
c
c     only a few iterations are generally needed for convergence
c     of the algorithm. if, however, the limit of 10 iterations
c     is reached, then the output par will contain the best
c     value obtained so far.
c
c     the subroutine statement is
c
c       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
c                        wa1,wa2)
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
c       delta is a positive input variable which specifies an upper
c         bound on the euclidean norm of d*x.
c
c       par is a nonnegative variable. on input par contains an
c         initial estimate of the levenberg-marquardt parameter.
c         on output par contains the final estimate.
c
c       x is an output array of length n which contains the least
c         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
c         for the output par.
c
c       sdiag is an output array of length n which contains the
c         diagonal elements of the upper triangular matrix s.
c
c       wa1 and wa2 are work arrays of length n.
c
c     subprograms called
c
c       minpack-supplied ... enorm,qrsolv
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int lmpar_(int n, double *r, int ldr, int *ipvt, double *diag, 
	   double *qtb, double delta, double *par, double *x, 
	   double *sdiag, double *wa1, double *wa2)
{
    const double p1 = .1;
    const double p001 = .001;
    double fp, sum, parc, parl;
    double temp, paru, dwarf;
    double d, gnorm, dxnorm;
    int iter, nsing;
    int i, j, k, l;
    
    /* dwarf is the smallest positive magnitude */
    dwarf = DBL_MIN;

    /* compute and store in x the Gauss-Newton direction: if the
       jacobian is rank-deficient, obtain a least squares solution 
    */

    nsing = n;
    for (j = 0; j < n; ++j) {
	wa1[j] = qtb[j];
	if (r[j + j * ldr] == 0.0 && nsing == n) {
	    nsing = j;
	}
	if (nsing < n) {
	    wa1[j] = 0.0;
	}
    }

    if (nsing > 0) {
	for (k = 0; k < nsing; ++k) {
	    j = nsing - k - 1;
	    wa1[j] /= r[j + j * ldr];
	    temp = wa1[j];
	    for (i = 0; i < j; ++i) {
		wa1[i] -= r[i + j * ldr] * temp;
	    }
	}
    }

    for (j = 0; j < n; ++j) {
	l = ipvt[j];
	x[l] = wa1[j];
    }

    /* initialize the iteration counter.
       evaluate the function at the origin, and test
       for acceptance of the Gauss-Newton direction 
    */

    iter = 0;
    for (j = 0; j < n; ++j) {
	wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm_(n, wa2);
    fp = dxnorm - delta;
    if (fp <= p1 * delta) {
	*par = 0.0;
	return 0;
    }

    /* if the jacobian is not rank deficient, the Newton
       step provides a lower bound, parl, for the zero of
       the function. otherwise set this bound to zero 
    */

    parl = 0.0;
    if (nsing >= n) {
	for (j = 0; j < n; ++j) {
	    l = ipvt[j];
	    wa1[j] = diag[l] * (wa2[l] / dxnorm);
	}
	for (j = 0; j < n; ++j) {
	    sum = 0.0;
	    for (i = 0; i < j; ++i) {
		sum += r[i + j * ldr] * wa1[i];
	    }
	    wa1[j] = (wa1[j] - sum) / r[j + j * ldr];
	}
	temp = enorm_(n, wa1);
	parl = fp / delta / temp / temp;
    }

    /* calculate an upper bound, paru, for the zero of the function */

    for (j = 0; j < n; ++j) {
	sum = 0.0;
	for (i = 0; i <= j; ++i) {
	    sum += r[i + j * ldr] * qtb[i];
	}
	l = ipvt[j];
	wa1[j] = sum / diag[l];
    }
    gnorm = enorm_(n, wa1);
    paru = gnorm / delta;
    if (paru == 0.0) {
	paru = dwarf / min(delta, p1);
    }

    /* if the input par lies outside of the interval (parl, paru),
       set par to the closer endpoint */

    *par = max(*par, parl);
    *par = min(*par, paru);
    if (*par == 0.0) {
	*par = gnorm / dxnorm;
    }

    /* beginning of iteration */

    while (1) {
	++iter;

	/* evaluate the function at the current value of par */

	if (*par == 0.0) {
	    d = p001 * paru;
	    *par = max(dwarf, d);
	}
	temp = sqrt(*par);
	for (j = 0; j < n; ++j) {
	    wa1[j] = temp * diag[j];
	}

	qrsolv_(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);

	for (j = 0; j < n; ++j) {
	    wa2[j] = diag[j] * x[j];
	}

	dxnorm = enorm_(n, wa2);
	temp = fp;
	fp = dxnorm - delta;

	/* If the function is small enough, accept the current value
	   of par. Also test for the exceptional cases where parl
	   is zero or the number of iterations has reached 10. 
	*/

	if (fabs(fp) <= p1 * delta || 
	    (parl == 0.0 && fp <= temp && temp < 0.0) ||
	    iter == 10) {
	    break;
	}

	/* compute the Newton correction */

	for (j = 0; j < n; ++j) {
	    l = ipvt[j];
	    wa1[j] = diag[l] * (wa2[l] / dxnorm);
	}
	for (j = 0; j < n; ++j) {
	    wa1[j] /= sdiag[j];
	    temp = wa1[j];
	    for (i = j+1; i < n; ++i) {
		wa1[i] -= r[i + j * ldr] * temp;
	    }
	}

	temp = enorm_(n, wa1);
	parc = fp / delta / temp / temp;

	/* depending on the sign of the function, update parl or paru */
	if (fp > 0.0) {
	    parl = max(parl, *par);
	} else if (fp < 0.0) {
	    paru = min(paru, *par);
	}

	/* compute an improved estimate for par */
	d = *par + parc;
	*par = max(parl, d);
    }

    if (iter == 0) {
	*par = 0.0;
    }

    return 0;
}

