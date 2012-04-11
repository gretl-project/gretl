#include "minpack.h"
#include <math.h>
#include <float.h>

#define log10e 0.43429448190325182765

static double d_log10 (double x)
{
    return log10e * log(x);
}

/*
c     chkder:
c
c     this subroutine checks the gradients of m nonlinear functions
c     in n variables, evaluated at a point x, for consistency with
c     the functions themselves. the user must call chkder twice,
c     first with mode = 1 and then with mode = 2.
c
c     mode = 1. on input, x must contain the point of evaluation.
c               on output, xp is set to a neighboring point.
c
c     mode = 2. on input, fvec must contain the functions and the
c                         rows of fjac must contain the gradients
c                         of the respective functions each evaluated
c                         at x, and fvecp must contain the functions
c                         evaluated at xp.
c               on output, err contains measures of correctness of
c                          the respective gradients.
c
c     the subroutine does not perform reliably if cancellation or
c     rounding errors cause a severe loss of significance in the
c     evaluation of a function. therefore, none of the components
c     of x should be unusually small (in particular, zero) or any
c     other value which may cause loss of significance.
c
c     the subroutine statement is
c
c       subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables.
c
c       x is an input array of length n.
c
c       fvec is an array of length m. on input when mode = 2,
c         fvec must contain the functions evaluated at x.
c
c       fjac is an m by n array. on input when mode = 2,
c         the rows of fjac must contain the gradients of
c         the respective functions evaluated at x.
c
c       ldfjac is a positive integer input parameter not less than m
c         which specifies the leading dimension of the array fjac.
c
c       xp is an array of length n. on output when mode = 1,
c         xp is set to a neighboring point of x.
c
c       fvecp is an array of length m. on input when mode = 2,
c         fvecp must contain the functions evaluated at xp.
c
c       mode is an integer input variable set to 1 on the first call
c         and 2 on the second. other values of mode are equivalent
c         to mode = 1.
c
c       err is an array of length m. on output when mode = 2,
c         err contains measures of correctness of the respective
c         gradients. if there is no severe loss of significance,
c         then if err(i) is 1.0 the i-th gradient is correct,
c         while if err(i) is 0.0 the i-th gradient is incorrect.
c         for values of err between 0.0 and 1.0, the categorization
c         is less certain. in general, a value of err(i) greater
c         than 0.5 indicates that the i-th gradient is probably
c         correct, while a value of err(i) less than 0.5 indicates
c         that the i-th gradient is probably incorrect.
c
c     subprograms called
c
c       minpack supplied ... dpmpar
c
c       fortran supplied ... dabs,dlog10,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
*/

int chkder_(int m, int n, double *x, 
	    double *fvec, double *fjac, int ldfjac, 
	    double *xp, double *fvecp, int mode, 
	    double *err)
{
    const double factor = 100.;

    int fjac_offset;
    double d;

    int i, j;
    double eps, temp, epsmch;

    --err;
    --fvecp;
    --fvec;
    --xp;
    --x;
    fjac_offset = 1 + ldfjac;
    fjac -= fjac_offset;

    epsmch = DBL_EPSILON;
    eps = sqrt(epsmch);

    if (mode == 1) {
	for (j = 1; j <= n; ++j) {
	    temp = eps * fabs(x[j]);
	    if (temp == 0.0) {
		temp = eps;
	    }
	    xp[j] = x[j] + temp;
	}
    } else {
	/* mode = 2 */
	double epsf = factor * epsmch;
	double epslog = d_log10(eps);

	for (i = 1; i <= m; ++i) {
	    err[i] = 0.0;
	}
	for (j = 1; j <= n; ++j) {
	    temp = fabs(x[j]);
	    if (temp == 0.0) {
		temp = 1.0;
	    }
	    for (i = 1; i <= m; ++i) {
		err[i] += temp * fjac[i + j * ldfjac];
	    }
	}
	for (i = 1; i <= m; ++i) {
	    temp = 1.0;
	    if (fvec[i] != 0.0 && fvecp[i] != 0.0 && 
		fabs(fvecp[i] - fvec[i]) >= epsf * fabs(fvec[i])) {
		d = fabs((fvecp[i] - fvec[i]) / eps - err[i]);
		temp = eps * d / (fabs(fvec[i]) + fabs(fvecp[i]));
	    }
	    err[i] = 1.0;
	    if (temp > epsmch && temp < eps) {
		err[i] = (d_log10(temp) - epslog) / epslog;
	    }
	    if (temp >= eps) {
		err[i] = 0.0;
	    }
	}
    }

    return 0;
}

