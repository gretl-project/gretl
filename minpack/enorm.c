#include "minpack.h"
#include <math.h>

/* given an n-vector x, this function calculates the
   euclidean norm of x.

   n is a positive integer input variable.

   x is an input array of length n.

   the euclidean norm is computed by accumulating the sum of
   squares in three different sums. the sums of squares for the
   small and large components are scaled so that no overflows
   occur. non-destructive underflows are permitted. underflows
   and overflows do not occur in the computation of the unscaled
   sum of squares for the intermediate components.
   the definitions of small, intermediate and large components
   depend on two constants, rdwarf and rgiant. the main
   restrictions on these constants are that rdwarf**2 not
   underflow and rgiant**2 not overflow. the constants
   given here are suitable for every known computer.

   argonne national laboratory. minpack project. march 1980.
   burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

double enorm_(int n, double *x)
{
    const double rdwarf = 3.834e-20;
    const double rgiant = 1.304e19;
    double agiant = rgiant / n;
    double d1, ret_val = 0.0;
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    double xabs, x1max = 0.0, x3max = 0.0;
    int i;

    /* Parameter adjustments */
    --x;

    for (i = 1; i <= n; ++i) {
	xabs = fabs(x[i]);
	if (xabs > rdwarf && xabs < agiant) {
	    /* sum for intermediate components */
	    s2 += xabs * xabs;
	} else if (xabs > rdwarf) {
	    /* sum for large components */
	    if (xabs > x1max) {
		d1 = x1max / xabs;
		s1 = 1.0 + s1 * (d1 * d1);
		x1max = xabs;
	    } else {
		d1 = xabs / x1max;
		s1 += d1 * d1;
	    }
	} else {
	    /* sum for small components */
	    if (xabs > x3max) {
		d1 = x3max / xabs;
		s3 = 1.0 + s3 * (d1 * d1);
		x3max = xabs;
	    } else if (xabs != 0.0) {
		d1 = xabs / x3max;
		s3 += d1 * d1;
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

