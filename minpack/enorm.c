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
c     enorm:
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c       double precision function enorm(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

double enorm_(int n, double *x)
{
    const double rdwarf = 3.834e-20;
    const double rgiant = 1.304e19;
    double agiant = rgiant / n;
    double d, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    double xabs, x1max = 0.0, x3max = 0.0;
    double ret_val = 0.0;
    int i;

    for (i = 0; i < n; i++) {
	xabs = fabs(x[i]);
	if (xabs > rdwarf && xabs < agiant) {
	    /* sum for intermediate components */
	    s2 += xabs * xabs;
	} else if (xabs > rdwarf) {
	    /* sum for large components */
	    if (xabs > x1max) {
		d = x1max / xabs;
		s1 = 1.0 + s1 * (d * d);
		x1max = xabs;
	    } else {
		d = xabs / x1max;
		s1 += d * d;
	    }
	} else {
	    /* sum for small components */
	    if (xabs > x3max) {
		d = x3max / xabs;
		s3 = 1.0 + s3 * (d * d);
		x3max = xabs;
	    } else if (xabs != 0.0) {
		d = xabs / x3max;
		s3 += d * d;
	    }
	}
    }

    /* calculation of norm */
    if (s1 != 0.0) {
	ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
    } else if (s2 != 0.0) {
	if (s2 >= x3max) {
	    ret_val = sqrt(s2 * (1.0 + x3max / s2 * (x3max * s3)));
	} else {
	    ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
	}
    } else {
	ret_val = x3max * sqrt(s3);
    }

    return ret_val;
}

