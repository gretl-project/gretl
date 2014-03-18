/* 
   This source file based on Minpack: initially converted from 
   fortran using f2c, then rendered into relatively idiomatic
   C with zero-based indexing throughout and pass-by-value for
   parameters that do not function as pointers. We also rely
   on <float.h> for the machine precision rather than Minpack's
   dpmpar().

   See README in this directory for the Minpack Copyright.

   Allin Cottrell, Wake Forest University, April 2012

   Modified in March 2014 by Jack Lucchetti to add a "quality"
   switch (see below).

*/

#include "minpack.h"
#include <math.h>
#include <float.h>

/*
c     fdjac2:
c
c     this subroutine computes a forward-difference approximation
c     to the m by n jacobian matrix associated with a specified
c     problem of m functions in n variables.
c
c     the subroutine statement is
c
c       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       quality is a [0-2] int: 0 gives old-style forward-difference,
c         1 bilateral difference and 2 gives 6th order Richardson
c         extrapolation; higher quality means lower speed.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int fdjac2_(S_fp fcn, int m, int n, int quality, 
	    double *x, double *fvec, 
	    double *fjac, int ldfjac, int *iflag, 
	    double epsfcn, double *wa, void *p)
{
    double eps = sqrt(max(epsfcn, DBL_EPSILON));
    double h, temp;
    int i, j;

    if (quality == 0) {
	/* plain old minpack fdjac2 */
	for (j = 0; j < n; j++) {
	    temp = x[j];
	    h = eps * fabs(temp);
	    if (h == 0.0) {
		h = eps;
	    }
	    x[j] = temp + h;
	    (*fcn)(m, n, x, wa, iflag, p);
	    if (*iflag < 0) {
		return 0;
	    }
	    x[j] = temp;
	    for (i = 0; i < m; i++) {
		fjac[i + j * ldfjac] = (wa[i] - fvec[i]) / h;
	    }
	}
    } else {
	int k, dim, d;
	int a[4], b[4];
	double y[m];

	if (quality == 1) {
	    eps *= 2;
	    dim = 2;
	    d = 2;
	    a[0] = b[0] = -1;
	    a[1] = b[1] = 1;
	} else if (quality == 2) {
	    dim = 4;
	    d = 12;
	    a[0] = -2; a[1] = -1; a[2] = 1; a[3] =  2;
	    b[0] =  1; b[1] = -8; b[2] = 8; b[3] = -1;
	} else {
	    *iflag = -1;
	    return 0;
	}

	for (j = 0; j < n; j++) {
	    temp = x[j];
	    h = eps * fabs(temp);
	    if (h == 0.0) {
		h = eps;
	    }

	    for (i = 0; i < m; i++) {
		y[i] = 0.0;
	    }
	    
	    for (k = 0; k < dim; k++) {
		x[j] = temp + a[k]*h;
		(*fcn)(m, n, x, wa, iflag, p);
		if (*iflag < 0) {
		    return 0;
		}
		for (i = 0; i < m; i++) {
		    y[i] += b[k]*wa[i];
		}
	    }

	    for (i = 0; i < m; i++) {
		fjac[i + j * ldfjac] = y[i] / (d*h);
	    }
	    
	    x[j] = temp;
	}
    }
    
    return 0;
}

