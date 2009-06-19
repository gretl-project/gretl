/* Complement to the cephes Bessel code.  Derived
   from GNU R, and in turn derived from netlib.
   Modified for gretl by Allin Cottrell, June 2009
*/

/* From http://www.netlib.org/specfun/rkbesl	
   Fortran translated by f2c,...
   Martin Maechler, ETH Zurich for GNU R
*/

#include "mconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define xmax_BESS_K   705.342 /* maximal x for UNscaled answer */
#define sqxmin_BESS_K 1.49e-154
#define M_SQRT_2dPI   0.797884560802865355879892119869 /* sqrt(2/pi) */

static int imin2 (int x, int y)
{
    return (x < y)? x : y;
}

static int imax2 (int x, int y)
{
    return (x > y)? x : y;
}

static double ftrunc(double x)
{
    return (x >= 0)? floor(x) : ceil(x);
}

/* Allin Cottrell mods below: pass x, alpha, nb and ize by value since
   they are not modified; also make K_bessel return the value of
   interest directly (in context, bk[nb-1]).
*/

static double K_bessel (double x, double alpha, int nb,
			int ize, double *bk, int *ncalc)
{
/*
  This routine calculates modified Bessel functions
  of the third kind, K_(N+ALPHA) (X), for non-negative
  argument X, and non-negative order N+ALPHA, with or without
  exponential scaling.

  Explanation of variables in the calling sequence

 X     - Non-negative argument for which
	 K's or exponentially scaled K's (K*EXP(X))
	 are to be calculated.	If K's are to be calculated,
	 X must not be greater than XMAX_BESS_K.
 ALPHA - Fractional part of order for which
	 K's or exponentially scaled K's (K*EXP(X)) are
	 to be calculated.  0 <= ALPHA < 1.0.
 NB    - Number of functions to be calculated, NB > 0.
	 The first function calculated is of order ALPHA, and the
	 last is of order (NB - 1 + ALPHA).
 IZE   - Type.	IZE = 1 if unscaled K's are to be calculated,
		    = 2 if exponentially scaled K's are to be calculated.
 BK    - Output vector of length NB.	If the
	 routine terminates normally (NCALC=NB), the vector BK
	 contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
	 or the corresponding exponentially scaled functions.
	 If (0 < NCALC < NB), BK(I) contains correct function
	 values for I <= NCALC, and contains the ratios
	 K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
 NCALC - Output variable indicating possible errors.
	 Before using the vector BK, the user should check that
	 NCALC=NB, i.e., all orders have been calculated to
	 the desired accuracy.	See error returns below.

 *******************************************************************

 Error returns

  In case of an error, NCALC != NB, and not all K's are
  calculated to the desired accuracy.

  NCALC < -1:  An argument is out of range. For example,
	NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= XMAX_BESS_K.
	In this case, the B-vector is not calculated,
	and NCALC is set to MIN0(NB,0)-2	 so that NCALC != NB.
  NCALC = -1:  Either  K(ALPHA,X) >= XINF  or
	K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) >= XINF.	 In this case,
	the B-vector is not calculated.	Note that again
	NCALC != NB.

  0 < NCALC < NB: Not all requested function values could
	be calculated accurately.  BK(I) contains correct function
	values for I <= NCALC, and contains the ratios
	K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.


 Intrinsic functions required are:

     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT


 Acknowledgement

	This program is based on a program written by J. B. Campbell
	(2) that computes values of the Bessel functions K of float
	argument and float order.  Modifications include the addition
	of non-scaled functions, parameterization of machine
	dependencies, and the use of more accurate approximations
	for SINH and SIN.

 References: "On Temme's Algorithm for the Modified Bessel
	      Functions of the Third Kind," Campbell, J. B.,
	      TOMS 6(4), Dec. 1980, pp. 581-586.

	     "A FORTRAN IV Subroutine for the Modified Bessel
	      Functions of the Third Kind of Real Order and Real
	      Argument," Campbell, J. B., Report NRC/ERB-925,
	      National Research Council, Canada.

  Latest modification: May 30, 1989

  Modified by: W. J. Cody and L. Stoltz
	       Applied Mathematics Division
	       Argonne National Laboratory
	       Argonne, IL  60439

*/
    /*---------------------------------------------------------------------
     * Mathematical constants
     *	A = LOG(2) - Euler's constant
     *	D = SQRT(2/PI)
     ---------------------------------------------------------------------*/
    static const double a = .11593151565841244881;

    /*---------------------------------------------------------------------
      P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA + Euler's constant
      Coefficients converted from hex to decimal and modified
      by W. J. Cody, 2/26/82 */
    static const double p[8] = { .805629875690432845,20.4045500205365151,
	    157.705605106676174,536.671116469207504,900.382759291288778,
	    730.923886650660393,229.299301509425145,.822467033424113231 };
    static const double q[7] = { 29.4601986247850434,277.577868510221208,
	    1206.70325591027438,2762.91444159791519,3443.74050506564618,
	    2210.63190113378647,572.267338359892221 };
    /* R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA) */
    static const double r[5] = { -.48672575865218401848,13.079485869097804016,
	    -101.96490580880537526,347.65409106507813131,
	    3.495898124521934782e-4 };
    static const double s[4] = { -25.579105509976461286,212.57260432226544008,
	    -610.69018684944109624,422.69668805777760407 };
    /* T    - Approximation for SINH(Y)/Y */
    static const double t[6] = { 1.6125990452916363814e-10,
	    2.5051878502858255354e-8,2.7557319615147964774e-6,
	    1.9841269840928373686e-4,.0083333333333334751799,
	    .16666666666666666446 };
    /*---------------------------------------------------------------------*/
    static const double estm[6] = { 
	52.0583,5.7607,2.7782,14.4303,185.3004, 9.3715 
    };
    static const double estf[7] = { 
	41.8341,7.1075,6.4306,42.511,1.35633,84.5096,20.
    };

    /* Local variables */
    int iend, i, j, k, m, mplus1, ii = 0;
    double x2by4, twox, c, blpha, ratio, wminf;
    double d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, twonu;
    double dm, bk1, bk2, nu, ret;

    nu = alpha;
    *ncalc = imin2(nb, 0) - 2;

    /* check for invalid parameter settings */
    if (nb <= 0 || (nu < 0 || nu >= 1) || (ize != 1 && ize != 2)) {
	fprintf(stderr, "internal error: bad arg for K_bessel\n");
	return NAN;
    }

    /* check the argument, x */
    if (x <= 0 || (ize == 1 && x > xmax_BESS_K)) {
	if (x < 0.0) {
	    mtherr("bessel_k", CEPHES_DOMAIN);
	    return NAN;
	} else if (x == 0.0) {
	    ret = INFINITY;
	} else {
	    /* would only have underflow */
	    ret = 0.0;
	}
	*ncalc = nb;
	return ret;
    }

    k = 0;
    if (nu < sqxmin_BESS_K) {
	nu = 0.;
    } else if (nu > .5) {
	k = 1;
	nu -= 1.;
    }

    twonu = nu + nu;
    iend = nb + k - 1;
    c = nu * nu;
    d3 = -c;

    if (x <= 1.0) {
	/* ------------------------------------------------------------
	   Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
	   Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
	   ------------------------------------------------------------ */
	d1 = 0.; d2 = p[0];
	t1 = 1.; t2 = q[0];
	for (i = 2; i <= 7; i += 2) {
	    d1 = c * d1 + p[i - 1];
	    d2 = c * d2 + p[i];
	    t1 = c * t1 + q[i - 1];
	    t2 = c * t2 + q[i];
	}
	d1 = nu * d1;
	t1 = nu * t1;
	f1 = log(x);
	f0 = a + nu * (p[7] - nu * (d1 + d2) / (t1 + t2)) - f1;
	q0 = exp(-nu * (a - nu * (p[7] + nu * (d1-d2) / (t1-t2)) - f1));
	f1 = nu * f0;
	p0 = exp(f1);
	/* -----------------------------------------------------------
	   Calculation of F0 =
	   ----------------------------------------------------------- */
	d1 = r[4];
	t1 = 1.;
	for (i = 0; i < 4; ++i) {
	    d1 = c * d1 + r[i];
	    t1 = c * t1 + s[i];
	}
	/* d2 := sinh(f1)/ nu = sinh(f1)/(f1/f0)
	 *	   = f0 * sinh(f1)/f1 */
	if (fabs(f1) <= .5) {
	    f1 *= f1;
	    d2 = 0.;
	    for (i = 0; i < 6; ++i) {
		d2 = f1 * d2 + t[i];
	    }
	    d2 = f0 + f0 * f1 * d2;
	} else {
	    d2 = sinh(f1) / nu;
	}
	f0 = d2 - nu * d1 / (t1 * p0);
	if (x <= 1e-10) {
	    /* ---------------------------------------------------------
	       X <= 1.0E-10
	       Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
	       --------------------------------------------------------- */
	    bk[0] = f0 + x * f0;
	    if (ize == 1) {
		bk[0] -= x * bk[0];
	    }
	    ratio = p0 / f0;
	    c = x * DBL_MAX;
	    if (k != 0) {
		/* ---------------------------------------------------
		   Calculation of K(ALPHA,X)
		   and  X*K(ALPHA+1,X)/K(ALPHA,X),	ALPHA >= 1/2
		   --------------------------------------------------- */
		*ncalc = -1;
		if (bk[0] >= c / ratio) {
		    goto finish;
		}
		bk[0] = ratio * bk[0] / x;
		twonu += 2.;
		ratio = twonu;
	    }
	    *ncalc = 1;
	    if (nb == 1) {
		goto finish;
	    }

	    /* -----------------------------------------------------
	       Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),
	       L = 1, 2, ... , NB-1
	       ----------------------------------------------------- */
	    *ncalc = -1;
	    for (i = 1; i < nb; ++i) {
		if (ratio >= c) {
		    goto finish;
		}
		bk[i] = ratio / x;
		twonu += 2.;
		ratio = twonu;
	    }
	    *ncalc = 1;
	    goto set_all_bk;
	} else {
	    /* ------------------------------------------------------
	       10^-10 < X <= 1.0
	       ------------------------------------------------------ */
	    c = 1.;
	    x2by4 = x * x / 4.;
	    p0 = .5 * p0;
	    q0 = .5 * q0;
	    d1 = -1.;
	    d2 = 0.;
	    bk1 = 0.;
	    bk2 = 0.;
	    f1 = f0;
	    f2 = p0;
	    do {
		d1 += 2.;
		d2 += 1.;
		d3 = d1 + d3;
		c = x2by4 * c / d2;
		f0 = (d2 * f0 + p0 + q0) / d3;
		p0 /= d2 - nu;
		q0 /= d2 + nu;
		t1 = c * f0;
		t2 = c * (p0 - d2 * f0);
		bk1 += t1;
		bk2 += t2;
	    } while (fabs(t1 / (f1 + bk1)) > DBL_EPSILON ||
		     fabs(t2 / (f2 + bk2)) > DBL_EPSILON);
	    bk1 = f1 + bk1;
	    bk2 = 2. * (f2 + bk2) / x;
	    if (ize == 2) {
		d1 = exp(x);
		bk1 *= d1;
		bk2 *= d1;
	    }
	    wminf = estf[0] * x + estf[1];
	}
    } else if (DBL_EPSILON * x > 1.0) {
	/* -------------------------------------------------
	   X > 1/EPS
	   ------------------------------------------------- */
	*ncalc = nb;
	bk1 = 1.0 / (M_SQRT_2dPI * sqrt(x));
	bk[nb-1] = bk1;
	goto finish;

    } else {
	/* -------------------------------------------------------
	   X > 1.0
	   ------------------------------------------------------- */
	twox = x + x;
	blpha = 0.0;
	ratio = 0.0;
	if (x <= 4.0) {
	    /* ----------------------------------------------------------
	       Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0
	       ----------------------------------------------------------*/
	    d2 = ftrunc(estm[0] / x + estm[1]);
	    m = (int) d2;
	    d1 = d2 + d2;
	    d2 -= .5;
	    d2 *= d2;
	    for (i = 2; i <= m; ++i) {
		d1 -= 2.0;
		d2 -= d1;
		ratio = (d3 + d2) / (twox + d1 - ratio);
	    }
	    /* -----------------------------------------------------------
	       Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
	       recurrence and K(ALPHA,X) from the wronskian
	       -----------------------------------------------------------*/
	    d2 = ftrunc(estm[2] * x + estm[3]);
	    m = (int) d2;
	    c = fabs(nu);
	    d3 = c + c;
	    d1 = d3 - 1.;
	    f1 = DBL_MIN;
	    f0 = (2. * (c + d2) / x + .5 * x / (c + d2 + 1.)) * DBL_MIN;
	    for (i = 3; i <= m; ++i) {
		d2 -= 1.;
		f2 = (d3 + d2 + d2) * f0;
		blpha = (1. + d1 / d2) * (f2 + blpha);
		f2 = f2 / x + f1;
		f1 = f0;
		f0 = f2;
	    }
	    f1 = (d3 + 2.) * f0 / x + f1;
	    d1 = 0.;
	    t1 = 1.;
	    for (i = 1; i <= 7; ++i) {
		d1 = c * d1 + p[i - 1];
		t1 = c * t1 + q[i - 1];
	    }
	    p0 = exp(c * (a + c * (p[7] - c * d1 / t1) - log(x))) / x;
	    f2 = (c + .5 - ratio) * f1 / x;
	    bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0;
	    if (ize == 1) {
		bk1 *= exp(-x);
	    }
	    wminf = estf[2] * x + estf[3];
	} else {
	    /* ---------------------------------------------------------
	       Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by
	       backward recurrence, for  X > 4.0
	       ----------------------------------------------------------*/
	    dm = ftrunc(estm[4] / x + estm[5]);
	    m = (int) dm;
	    d2 = dm - .5;
	    d2 *= d2;
	    d1 = dm + dm;
	    for (i = 2; i <= m; ++i) {
		dm -= 1.;
		d1 -= 2.;
		d2 -= d1;
		ratio = (d3 + d2) / (twox + d1 - ratio);
		blpha = (ratio + ratio * blpha) / dm;
	    }
	    bk1 = 1. / ((M_SQRT_2dPI + M_SQRT_2dPI * blpha) * sqrt(x));
	    if (ize == 1)
		bk1 *= exp(-x);
	    wminf = estf[4] * (x - fabs(x - estf[6])) + estf[5];
	}
	/* ---------------------------------------------------------
	   Calculation of K(ALPHA+1,X)
	   from K(ALPHA,X) and  K(ALPHA+1,X)/K(ALPHA,X)
	   --------------------------------------------------------- */
	bk2 = bk1 + bk1 * (nu + .5 - ratio) / x;
    }

    /*--------------------------------------------------------------------
      Calculation of 'NCALC', K(ALPHA+I,X),	I  =  0, 1, ... , NCALC-1,
      &	  K(ALPHA+I,X)/K(ALPHA+I-1,X),	I = NCALC, NCALC+1, ... , NB-1
      -------------------------------------------------------------------*/

    *ncalc = nb;
    bk[0] = bk1;
    if (iend == 0) {
	goto finish;
    }

    j = 1 - k;
    if (j >= 0) {
	bk[j] = bk2;
    }

    if (iend == 1) {
	goto finish;
    }

    m = imin2((int) (wminf - nu),iend);
    for (i = 2; i <= m; ++i) {
	t1 = bk1;
	bk1 = bk2;
	twonu += 2.0;
	if (x < 1.) {
	    if (bk1 >= DBL_MAX / twonu * x) {
		break;
	    }
	} else if (bk1 / x >= DBL_MAX / twonu) {
	    break;
	}
	bk2 = twonu / x * bk1 + t1;
	ii = i;
	if (++j >= 0) {
	    bk[j] = bk2;
	}
    }

    m = ii;
    if (m == iend) {
	goto finish;
    }
    ratio = bk2 / bk1;
    mplus1 = m + 1;
    *ncalc = -1;
    for (i = mplus1; i <= iend; ++i) {
	twonu += 2.0;
	ratio = twonu / x + 1./ratio;
	if (++j >= 1) {
	    bk[j] = ratio;
	} else {
	    if (bk2 >= DBL_MAX / ratio) {
		goto finish;
	    }
	    bk2 *= ratio;
	}
    }
    *ncalc = imax2(1, mplus1 - k);
    if (*ncalc == 1) {
	bk[0] = bk2;
    }
    if (nb == 1) {
	goto finish;
    }

 set_all_bk:

    /* fill out values not directly computed: kicks in
       only if ncalc < nb 
    */
    for (i = *ncalc; i < nb; ++i) {
	bk[i] *= bk[i-1];
	(*ncalc)++;
    }

 finish:

    return bk[nb-1];
}

static double *zero_dbl_array (int n)
{
    double *x = malloc(n * sizeof *x);

    if (x == NULL) {
	mtherr ("allocation", CEPHES_UNKNOWN);
    } else {
	int i;

	for (i=0; i<n; i++) {
	    x[i] = 0.0;
	}
    }

    return x;
}

/* ize == 1 -> unscaled; 
   ize == 2 -> exponentially scaled 
*/

double netlib_bessel_K (double v, double x, int ize)
{
    int nb, ncalc;
    double *bk;

    if (isnan(x) || isnan(v)) {
	return x + v;
    }

    if (x < 0) {
	mtherr("netlib_bessel_K", CEPHES_DOMAIN);
	return NAN;
    }

    if (v < 0) {
	v = -v;
    }

    nb = 1 + (int) floor(v); /* nb-1 <= |v| < nb */
    v -= (nb - 1);

    bk = zero_dbl_array(nb);
    if (bk == NULL) {
	return NAN;
    }

    x = K_bessel(x, v, nb, ize, bk, &ncalc);

    if (ncalc != nb) {
	if (ncalc < 0) {
	    mtherr("netlib_bessel_K", CEPHES_TLOSS);
	    x = NAN;
	} else {
	    mtherr("netlib_bessel_K", CEPHES_PLOSS);
	}
    }

    free(bk);

    return x;
}


