#ifndef LIBPROB_H
#define LIBPROB_H

/* area under the left hand tail (from 0 to x)
 * of the Chi square probability density function with
 * v degrees of freedom.
 */
double chdtr (double v, double x);

/*
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 */
double chdtri (double df, double y);

/*
 * Returns the area from x to infinity under the F density
 * function (also known as Snedcor's density or the
 * variance ratio density).
 */
double fdtrc (int ia, int ib, double x);

/*
 * Finds the F density argument x such that the integral
 * from x to infinity of the F density is equal to the
 * given probability p.
 */
double fdtri (int ia, int ib, double y);

/*
 * Computes the integral from minus infinity to t of the Student
 * t distribution with integer k > 0 degrees of freedom.
 */
double stdtr (int k, double t);

/*
 * Given probability p, finds the argument t such that stdtr(k,t)
 * is equal to p.
 */
double stdtri (int k, double p);

/*
 * Returns the area under the Gaussian probability density
 * function, integrated from minus infinity to x.
 */
double ndtr (double a);

/*
 * Returns the argument, x, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to x) is equal to y.
 */
double ndtri (double y0);

/*
 * Returns gamma function of the argument.  The result is
 * correctly signed.
*/

double cephes_gamma (double x); /* alias for gamma() */

#endif /* LIBPROB_H */





