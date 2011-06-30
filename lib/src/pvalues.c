/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/*  pvalues.c - routines relating to computation of pvalues
    of sample statistics
*/  

#include "libgretl.h" 
#include "../../cephes/libprob.h"

#include <errno.h>

/**
 * SECTION:pvalues
 * @short_description: probability values for test statistics and
 * related functionality
 * @title: P-values
 * @include: libgretl.h
 *
 * Libgretl uses the cephes library, developed by Stephen
 * Moshier, as the basic engine for much of the functionality
 * herein. We add some extra distributions, and wrap the cephes
 * functions for ease of use with libgretl (e.g. on failure
 * they return the libgretl missing value code, %NADBL).
 */

/**
 * gamma_function:
 * @x: argument.
 *
 * Returns: the gamma function of @x, or #NADBL on failure.
 */

double gamma_function (double x)
{
    double ret = cephes_gamma(x);

    if (get_cephes_errno()) {
	ret = NADBL;
    }

    return ret;
}

/**
 * ln_gamma:
 * @x: argument.
 *
 * Returns: the log gamma function of @x, or #NADBL on failure.
 */

double ln_gamma (double x)
{
    double ret = cephes_lgamma(x);

    if (get_cephes_errno()) {
	ret = NADBL;
    }

    return ret;
}

/**
 * digamma:
 * @x: argument.
 *
 * Returns: the digamma (or Psi) function of @x, or #NADBL on failure.
 */

double digamma (double x)
{
    double ret = psi(x);

    if (get_cephes_errno()) {
	ret = NADBL;
    }

    return ret;
}
 
/**
 * binomial_cdf:
 * @p: probability of success on each trial.
 * @n: number of trials.
 * @k: maximum number of successes. 
 *
 * Returns: the probability of @k or less successes on
 * @n trials given binomial probability @p, or
 * #NADBL on failure.
 */

double binomial_cdf (double p, int n, int k)
{
    double x = NADBL;

    if (p >= 0 && n >= 0 && k >= 0) {
	x = bdtr(k, n, p);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

/**
 * binomial_cdf_comp:
 * @p: probability of success on each trial.
 * @n: number of trials.
 * @k: maximum number of successes.
 *
 * Returns: the probability of @k + 1 or more successes on
 * @n trials given binomial probability @p, or
 * #NADBL on failure.
 */

double binomial_cdf_comp (double p, int n, int k)
{
    double x = NADBL;

    if (p >= 0 && n >= 0 && k >= 0) {
	x = bdtrc(k, n, p);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

/*  
    gives the event probability p such that the sum of the terms 0
    through k of the Binomial probability density, for n trials, is
    equal to the given cumulative probability a.
*/

static double binomial_cdf_inverse (int n, int k, double a)
{
    double p = NADBL;

    if (a >= 0 && n >= 0 && k >= 0) {
	p = bdtri(k, n, a);
	if (get_cephes_errno()) {
	    p = NADBL;
	}
    }

    return p;
}

/* The following is no doubt horribly inefficient */

static double binomial_critval (double p, int n, double a)
{
    double pk, ac = 1 - a;
    int k;

    if (n <= 0 || p <= 0 || p >= 1 || a <= 0 || a >= 1) {
	return NADBL;
    }

    for (k=n; k>0; k--) {
	pk = binomial_cdf(p, n, k);
	if (pk < ac) {
	    break;
	}
    }

    return (double) (k + 1);
}

/**
 * binomial_pmf:
 * @p: success probability.
 * @n: number of trials.
 * @k: number of successes.
 * 
 * Returns: the probability mass for @k successes in @n
 * binomial trials with success probability @p.
 */

double binomial_pmf (double p, int n, int k)
{
    double pm = binomial_cdf(p, n, k);

    if (k > 0 && !na(pm)) {
	pm -= binomial_cdf(p, n, k - 1);
    }

    return pm;
}

static int binomial_pmf_array (double p, int n, double *x, int T)
{
    double pm;
    int i, k;

    for (i=0; i<T; i++) {
	k = x[i];
	pm = binomial_cdf(p, n, k);
	if (k > 0 && !na(pm)) {
	    pm -= binomial_cdf(p, n, k - 1);
	}
	x[i] = pm;
    }

    return 0;
}

/**
 * x_factorial:
 * @x: input value.
 * 
 * Returns: the factorial of int(x), cast to a double, or
 * NADBL on failure.
 */

double x_factorial (double x)
{
    double fact;
    int n = x;

    if (x < 0) {
	fact = NADBL;
    } else if (x > 12) {
	fact = cephes_gamma(1 + x);
	if (get_cephes_errno()) {
	    fact = NADBL;
	}
    } else if (n == 0) {
	fact = 1;
    } else {
	fact = n;
	while (--n > 1) {
	    fact *= n;
	}
    }

    return fact;
}

/**
 * log_x_factorial:
 * @x: input value.
 * 
 * Returns: the log of the factorial of int(x), cast to a double, or
 * NADBL on failure.
 */

double log_x_factorial (double x)
{
    double lfact;
    int n = x;

    if (x < 0) {
	lfact = NADBL;
    } else if (x > 12) {
	lfact = cephes_lgamma(1 + x);
	if (get_cephes_errno()) {
	    lfact = NADBL;
	}
    } else if (n == 0) {
	lfact = 0;
    } else {
	lfact = n;
	while (--n > 1) {
	    lfact *= n;
	}
	lfact = log(lfact);
    }

    return lfact;
}

/**
 * tcrit95:
 * @df: degrees of freedom.
 * 
 * Returns: the two-sided 95 percent critical value for the t 
 * distribution with @df degrees of freedom, or #NADBL on
 * failure.
 */

double tcrit95 (int df)
{
    double x = NADBL;

    if (df > 0) {
	x = stdtri(df, 0.975);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

/**
 * rhocrit95:
 * @n: sample size.
 *
 * Computes the two-sided 5 percent critical value of the sample
 * correlation coefficient for a sample of size @n. This is based 
 * on the inverse of the function which maps from the correlation 
 * coefficient, r, to a student t statistic, namely
 *
 *            t = r / sqrt[(1—r^2) / (n-2)]
 *
 * The inverse is r = sqrt(t^2 / (t^2 + n - 2)).
 * 
 * Returns: the critical value, or #NADBL on failure.
 */

double rhocrit95 (int n)
{
    double rc = NADBL;

    if (n - 2 > 0) {
	double tc = stdtri(n - 2, 0.975);

	if (get_cephes_errno() == 0) {
	    double tc2 = tc * tc;

	    rc = sqrt(tc2 / (tc2 + n - 2));
	}
    }
    
    return rc;
}

/**
 * normal_pvalue_2:
 * @x: double-precision value.
 *
 * Calculates the two-sided p-value for @x in relation to the
 * standard normal distribution.
 *
 * Returns: 2 times (1 minus the value of the standard normal
 * CDF evaluated at abs(@x)), or 0 on underflow.
 */

double normal_pvalue_2 (double x)
{
    double p = (x < 0)? ndtr(x) : ndtr(-x);

    return 2 * p;
}

/**
 * normal_pvalue_1:
 * @x: double-precision value.
 * 
 * Calculates the one-sided p-value for @x in relation to the
 * standard normal distribution (that is, the probability that a 
 * random variable distributed as N(0, 1) is greater than @x).
 *
 * Returns: 1 minus the value of the standard normal CDF
 * evaluated at @x.
 */

double normal_pvalue_1 (double x)
{
    double p = ndtr(x);

    if (get_cephes_errno()) {
	p = NADBL;
    } else {
	p = 1 - p;
    }

    return p;
}

/**
 * student_cdf:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral from minus infinity to @x of
 * the Student's t distribution with @df degrees of freedom, or
 * #NADBL on failure.
 */

double student_cdf (double df, double x)
{
    double p = NADBL;

    if (df > 0) {
	p = stdtr(df, x);
	if (get_cephes_errno()) {
	    p = NADBL;
	}
    }

    return p;
}

/**
 * student_cdf_comp:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral from @x to infinity of
 * the t distribution with @df degrees of freedom, or
 * #NADBL on failure.
 */

static double student_cdf_comp (double df, double x)
{
    double p = NADBL;

    if (df > 0) {
	if (x > 0) {
	    p = stdtr(df, -x);
	    if (get_cephes_errno()) {
		p = NADBL;
	    }
	} else {
	    p = stdtr(df, x);
	    if (get_cephes_errno()) {
		p = NADBL;
	    } else {
		p = 1 - p;
	    }
	}
    }

    return p;
}

/**
 * normal_cdf_comp:
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral from @x to infinity of the standard 
 * normal distribution, or #NADBL on failure.
 */

double normal_cdf_comp (double x)
{
    double p;

    if (x > 0) {
	p = ndtr(-x);
	if (get_cephes_errno()) {
	    p = NADBL;
	}	
    } else {
	p = ndtr(x);
	if (get_cephes_errno()) {
	    p = NADBL;
	} else {
	    p = 1 - p;
	}
    } 

    return p;
}

/**
 * student_pvalue_1:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the probability that t(@df) is greater than @x,
 * or #NADBL on failure.
 */

double student_pvalue_1 (double df, double x)
{
    double p = NADBL;

    if (df > 0) {
	p = stdtr(df, x);
	if (get_cephes_errno()) {
	    p = NADBL;
	} else {
	    p = 1 - p;
	}
    }

    return p;
}

/**
 * student_pvalue_2:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the probability that t(@df) is greater than @x
 * (two-sided, using the absolute value of @x), or
 * #NADBL on failure.
 */

double student_pvalue_2 (double df, double x)
{
    double p = NADBL;

    if (df > 0) {
	if (x < 0.0) {
	    p = stdtr(df, x);
	} else {
	    p = stdtr(df, -x);
	}
	if (get_cephes_errno()) {
	    p = NADBL;
	} else {
	    p *= 2;
	}
    }

    return p;
}

/**
 * student_critval:
 * @df: degrees of freedom.
 * @a: right-tail probability.
 *
 * Returns: the argument x such that the integral from x to 
 * infinity of the t(@df) density is equal to the given
 * probability @a, or #NADBL on failure.
 */

double student_critval (double df, double a)
{
    double x;

    if (df < 0) {
	return NADBL;
    }    

    if (a > .10) {
	x = stdtri(df, 1 - a);
    } else {
	x = -stdtri(df, a);
    }

    if (get_cephes_errno()) {
	x = NADBL;
    } 

    return x;
}

/**
 * student_cdf_inverse:
 * @df: degrees of freedom.
 * @a: probability.
 *
 * Returns: the argument x such that the integral from 
 * minus infinity to @x of the t(@df) density is equal to 
 * the given probability @a, or #NADBL on failure.
 */

double student_cdf_inverse (double df, double a)
{
    double x;

    if (df < 0) {
	return NADBL;
    }

    x = stdtri(df, a);

    if (get_cephes_errno()) {
	x = NADBL;
    } 

    return x;
}

/**
 * chisq_cdf:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral from 0 to @x of the chi-square
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

double chisq_cdf (int df, double x)
{
    double p = NADBL;

    if (df > 0 && x >= 0) {
	p = chdtr(df, x);
	if (get_cephes_errno()) {
	    p = NADBL;
	}
    }

    return p;
}

/**
 * chisq_cdf_comp:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral from @x to infinity of the chi-square
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

double chisq_cdf_comp (int df, double x)
{
    double p = NADBL;

    if (df > 0 && x >= 0) {
	p = chdtrc(df, x);
	if (get_cephes_errno()) {
	    p = NADBL;
	}
    }

    return p;
}

/**
 * chisq_critval:
 * @df: degrees of freedom.
 * @a: right-tail probability.
 * 
 * Returns: the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given probability @a, or #NADBL on failure.
 */

static double chisq_critval (int df, double a)
{
    double x = NADBL;

    if (df > 0 && a >= 0) {
	x = chdtri(df, a);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

static double chisq_cdf_inverse (int df, double a)
{
    double x = NADBL;

    if (df > 0 && a >= 0) {
	x = chdtri(df, 1 - a);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

/**
 * snedecor_cdf:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral of the F distribution with @dfn and
 * @dfd degrees of freedom, from 0 to @x, or #NADBL on failure.
 */

static double snedecor_cdf (int dfn, int dfd, double x)
{
    double p = NADBL;

    if (dfn > 0 && dfd > 0 && x >= 0) {
	p = fdtr(dfn, dfd, x);
	if (get_cephes_errno()) {
	    p = NADBL;
	}
    }

    return p;
}

/**
 * snedecor_cdf_comp:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the integral of the F distribution with @dfn and
 * @dfd degrees of freedom, from @x to infinity, or #NADBL 
 * on failure.
 */

double snedecor_cdf_comp (int dfn, int dfd, double x)
{
    double p = NADBL;

    if (dfn > 0 && dfd > 0 && x >= 0) {
	p = fdtrc(dfn, dfd, x);
	if (get_cephes_errno()) {
	    p = NADBL;
	}
    }

    return p;
}

/**
 * snedecor_critval:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @a: right-tail probability.
 * 
 * Returns: the F argument x such that the integral
 * from x to infinity of the F density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double snedecor_critval (int dfn, int dfd, double a)
{
    double x = NADBL;

    if (dfn > 0 && dfd > 0 && a >= 0) {
	x = fdtri(dfn, dfd, a);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

static double snedecor_cdf_inverse (int dfn, int dfd, double a)
{
    double x = NADBL;

    if (dfn > 0 && dfd > 0 && a >= 0) {
	x = fdtri(dfn, dfd, 1 - a);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

/* PDFs */

static double Binv (double p, double q)
{
    double f = NADBL;

    errno = 0;

    if (p > 0 && q > 0) {
	double x1 = ln_gamma(p + q);
	double x2 = ln_gamma(p);
	double x3 = ln_gamma(q);

	if (!na(x1) && !na(x2) && !na(x3)) {
	    f = exp(x1 - x2 - x3);
	    if (errno) {
		f = NADBL;
	    }
	}
    }

    return f;
}    

static int snedecor_pdf_array (int v1, int v2, double *x, int n)
{
    int i, err = 0;

    errno = 0;

    if (v1 > 0 && v2 > 0) {
	double xm = v1, xn = v2;
	double vr = xm / xn;
	double x1 = Binv(0.5*xm, 0.5*xn);
	double x2 = pow(vr, 0.5*xm);
	double x3, x4;

	if (errno) {
	    err = E_NAN;
	} else {
	    for (i=0; i<n; i++) {
		if (!na(x[i]) && x[i] > 0) {
		    errno = 0;
		    x3 = pow(x[i], 0.5*xm - 1.0);
		    x4 = pow(1.0 + vr * x[i], 0.5 * (xm+xn));
		    x[i] = x1 * x2 * x3 / x4;
		    if (errno) {
			x[i] = NADBL;
		    }
		} else {
		    x[i] = NADBL;
		}
	    }
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    errno = 0;

    return err;
}

static double snedecor_pdf (int m, int n, double x)
{
    snedecor_pdf_array(m, n, &x, 1);

    return x;
}

static int chisq_pdf_array (int m, double *x, int n)
{
    int i, err = 0;

    errno = 0;

    if (m > 0) {
	double m2 = m / 2.0;
	double x1 = pow(.5, m2);
	double x2 = gamma_function(m2);
	double x3, x4;

	if (errno) {
	    err = E_NAN;
	} else {
	    for (i=0; i<n; i++) {
		if (!na(x[i]) && x[i] > 0) {
		    errno = 0;
		    x3 = pow(x[i], m2 - 1.0);
		    x4 = exp(-x[i] / 2.0);
		    x[i] = (x1/x2) * x3 * x4;
		    if (errno) {
			x[i] = NADBL;
		    }
		} else {
		    x[i] = NADBL;
		}
	    }
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    errno = 0;

    return err;
}

static double chisq_pdf (int m, double x)
{
    chisq_pdf_array(m, &x, 1);

    return x;
}

static int student_pdf_array (double m, double *x, int n)
{
    int i, err = 0;

    errno = 0;

    if (m > 0 && !na(m)) {
	double x1 = Binv(0.5 * m, 0.5) / sqrt(m);
	double x3 = 0.5 * (m + 1.0);
	double x2;

	if (errno) {
	    err = E_NAN;
	} else {
	    for (i=0; i<n; i++) {
		if (!na(x[i])) {
		    errno = 0;
		    x2 = m / (m + x[i] * x[i]);
		    x[i] = x1 * pow(x2, x3);
		    if (errno) {
			x[i] = NADBL;
		    }
		}
	    }
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    errno = 0;

    return err;
}

static double student_pdf (double m, double x)
{
    student_pdf_array(m, &x, 1);

    return x;
}

static int weibull_pdf_array (double k, double l, 
			      double *x, int n)
{
    int i, err = 0;

    errno = 0;

    if (!na(k) && k > 0 && !na(l) && l > 0) {
	double x1 = k / l;
	double x2, x3, x4;

	if (errno) {
	    err = E_NAN;
	} else {
	    for (i=0; i<n; i++) {
		if (!na(x[i]) && x[i] >= 0) {
		    errno = 0;
		    x3 = x[i] / l;
		    x2 = pow(x3, k - 1.0);
		    x3 *= x2;
		    x4 = exp(-x3);
		    x[i] = x1 * x2 * x4;
		    if (errno) {
			x[i] = NADBL;
		    }
		} else {
		    x[i] = NADBL;
		}
	    }
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    errno = 0;

    return err;
}

static double weibull_pdf (double k, double l, double x)
{
    weibull_pdf_array(k, l, &x, 1);

    return x;
}

/**
 * normal_cdf:
 * @x: double-precision value.
 * 
 * Returns: the value of the standard normal CDF evaluated
 * at @x, or #NADBL on failure.
 */

double normal_cdf (double x)
{
    double y = ndtr(x);

    if (get_cephes_errno()) {
	y = NADBL;
    } else if (y == 1.0) {
	/* do we want to do this? */
	y = 0.9999999999999999;
    }

    return y;
}

/**
 * normal_cdf_inverse:
 * @x: double-precision value.
 * 
 * Returns: the argument, y, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to y) is equal to x, or #NADBL on failure.
 */

double normal_cdf_inverse (double x)
{
    double y = ndtri(x);

    if (get_cephes_errno()) {
	y = NADBL;
    }

    return y;
}

/**
 * normal_critval:
 * @a: right-tail probability.
 *
 * Returns: the argument z such that the integral from z to 
 * infinity of the standard normal density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double normal_critval (double a)
{
    double z;

    if (a > 0.10) {
	z = ndtri(1.0 - a);
    } else {
	z = -ndtri(a);
    }

    if (get_cephes_errno()) {
	z = NADBL;
    } 

    return z;
}

static int normal_pdf_array (double *x, int n)
{
    double s = 1.0 / SQRT_2_PI;
    int i;

    for (i=0; i<n; i++) {
	if (!na(x[i])) {
	    errno = 0;
	    x[i] = s * exp(-0.5 * x[i] * x[i]);
	    if (errno) {
		x[i] = NADBL;
	    }
	}
    }

    errno = 0;

    return 0;
}

/**
 * normal_pdf:
 * @x: double-precision value.
 * 
 * Returns: the value of the standard normal PDF evaluated
 * at @x.
 */

double normal_pdf (double x)
{
    return exp(-0.5 * x * x) / SQRT_2_PI;
}

/**
 * log_normal_pdf:
 * @x: double-precision value.
 * 
 * Returns: the value of the log-normal PDF evaluated
 * at @x.
 */

double log_normal_pdf (double x)
{
    return -0.5 * x * x - LN_SQRT_2_PI;
}

/**
 * invmills:
 * @x: double-precision value.
 *
 * Adapted by putting together code from gsl and TDA (Univ. Bochum). 
 * The latter is, in turn, based on 
 * A. V. Swan, The Reciprocal of Mills's Ratio, Algorithm AS 17,
 * Journal of the Royal Statistical Society. Series C (Applied Statistics), 
 * Vol. 18, No. 1 (1969), 115-116.
 *
 * Returns: the inverse Mills ratio, that is the ratio between the
 * normal density function and the complement of the distribution 
 * function, both evaluated at @x.
 */

#define SQRT_HALF_PI 1.2533141373155002512078826424
#define MILLS_BOTTOM -22.9
#define MILLS_TOP 25
#define MILLS_EPS 1.0e-09

double invmills (double x)
{
    double a, a0, a1, a2;
    double b, b0, b1, b2;
    double r, s, t, d;
    double ret;

    if (x == 0.0) {
        return 1.0 / SQRT_HALF_PI;
    }

    if (x < MILLS_BOTTOM) {
	return 0;
    }

    if (x > MILLS_TOP) {
	a0 = 1.0/(x * x);
	a1 = 1.0 - 9.0 * a0 * (1.0 - 11.0 * a0);
	a2 = 1.0 - 5.0 * a0 * (1.0 - 7.0 * a0 * a1);
	d = 1.0 - a0 * (1.0 - 3.0 * a0 * a2);
	return x / d;
    }

    d = (x < 0.0)? -1.0 : 1.0;
    x = fabs(x);

    if (x <= 2.0) {
	s = 0.0;
	a = 1.0;
	r = t = x;
	b = x * x;
	while (fabs(s-t) > MILLS_EPS) {
	    a += 2.0;
	    s = t;
	    r *= b / a;
	    t += r;
	}
	ret = 1.0 / (SQRT_HALF_PI * exp(0.5 * b) - d * t);
    } else {
	a = 2.0;
	r = s = b1 = x;
	a1 = x * x + 1.0;
	a2 = x * (a1 + 2.0);
	b2 = a1 + 1.0;
	t  = a2 / b2;
	while (fabs(r-t) > MILLS_EPS && fabs(s-t) > MILLS_EPS) {
	    a += 1.0;
	    a0 = a1;
	    a1 = a2;
	    a2 = x * a1 + a * a0;
	    b0 = b1;
	    b1 = b2;
	    b2 = x * b1 + a * b0;
	    r  = s;
	    s  = t;
	    t  = a2 / b2;
	}
	ret = t;
	if (d < 0.0) {
	    ret /= (2.0 * SQRT_HALF_PI * exp(0.5 * x * x) * t - 1.0);
	}
    }

    return ret;
}

/**
 * bvnorm_cdf:
 * @rho: correlation coefficient.
 * @a: abscissa value, first Gaussian r.v.
 * @b: abscissa value, second Gaussian r.v.
 *
 * Ripped and adapted from Gnumeric, with a bug corrected for the case
 * (a * b < 0) && (rho < 0).
 *
 * Returns: for (x, y) a bivariate standard Normal rv with correlation
 * coefficient @rho, the joint probability that (x < @a) and (y < @b), or
 * #NADBL on failure.
 */

double bvnorm_cdf (double rho, double a, double b)
{
    static const double x[] = {0.24840615, 0.39233107, 0.21141819, 
			       0.03324666, 0.00082485334};
    static const double y[] = {0.10024215, 0.48281397, 1.0609498, 
			       1.7797294, 2.6697604};
    double ret = NADBL;
    double a1, b1, den;
    int i, j;

    if (fabs(rho) > 1) {
	return NADBL;
    }	

    if (rho == 0.0) {
	/* joint prob is just the product of the marginals */
	return normal_cdf(a) * normal_cdf(b);
    }

    if (rho == 1.0) {
	/* the two variables are in fact the same */
	return normal_cdf(a<b ? a : b);
    }
    
    if (rho == -1.0) {
	/* the two variables are perfectly negatively correlated: 
	   P(x<a, y<b) = P((x<a) && (x>b)) = P(x \in (b,a))
	*/
	ret = (a<=b) ? 0 : normal_cdf(a)-normal_cdf(b);
	return ret;
    }
    
    den = sqrt(2.0 * (1 - rho * rho));

    a1 = a / den;
    b1 = b / den;

    if (a <= 0 && b <= 0 && rho < 0) {
	/* standard case */
	double sum = 0.0;

	for (i=0; i<5; i++) {
	    for (j=0; j<5; j++) {
		sum += x[i] * x[j] * 
		    exp (a1 * (2 * y[i] - a1) + 
			 b1 * (2 * y[j] - b1) + 
			 2 * rho * (y[i] - a1) * (y[j] - b1));
	    }
	}
	ret = (sqrt(1 - (rho * rho)) / M_PI * sum);
    } else if (a <= 0 && b >= 0 && rho > 0) {
	ret = normal_cdf(a) - bvnorm_cdf(-rho, a, -b);
    } else if (a >= 0 && b <= 0 && rho > 0) {
	ret = normal_cdf(b) - bvnorm_cdf(-rho, -a, b);
    } else if (a >= 0 && b >= 0 && rho < 0) {
	ret = normal_cdf(a) + normal_cdf(b) - 1 + bvnorm_cdf(rho, -a, -b);
    } else if ((a * b * rho) > 0) {
	int sgna = (a < 0)? -1 : 1;
	int sgnb = (b < 0)? -1 : 1;
	double rho1, rho2, tmp, delta;

	tmp = sqrt((a * a) - 2 * rho * a * b + (b * b));
	rho1 = (rho * a - b) * sgna / tmp;
	rho2 = (rho * b - a) * sgnb / tmp;
	delta = (sgna * sgnb && (rho > 0))? 0 : 0.5;

	ret = (bvnorm_cdf(rho1, a, 0) + bvnorm_cdf(rho2, b, 0) - delta);
    }    

    return ret;
}

/**
 * gamma_cdf:
 * @s1: first parameter.
 * @s2: second parameter.
 * @x: reference value.
 * @control: see below.
 *
 * Calculates the value of the CDF of the gamma distribution
 * at @x.  If @control equals 1, then it is assumed that the
 * parameters @s1 and @s2 represent the shape and scale,
 * respectively, otherwise it is assumed they give mean and
 * variance.

 * Returns: the calculated probability, or #NADBL on failure.
 */

double gamma_cdf (double s1, double s2, double x, int control)
{
    double scale, shape, p;

    if (control == 1) {
	shape = s1; 
	scale = s2; 
    } else {
	scale = s2 / s1; 
	shape = s1 / scale; 
    }

    /* for the cephes functions, the parameterization is
       inverse-scale (or "rate"), shape */

    p = gdtr(1.0 / scale, shape, x);
    if (get_cephes_errno()) {
	p = NADBL;
    }

    return p;
}

/**
 * gamma_cdf_comp:
 * @s1: first parameter.
 * @s2: second parameter.
 * @x: reference value.
 * @control: see below.
 *
 * Calculates the complement of the CDF of the gamma distribution
 * at @x.  If @control equals 1, then it is assumed that the
 * parameters @s1 and @s2 represent the shape and scale,
 * respectively, otherwise it is assumed they give mean and
 * variance.

 * Returns: the calculated probability, or #NADBL on failure.
 */

double gamma_cdf_comp (double s1, double s2, double x, int control)
{
    double shape, scale, p;

    if (control == 1) {
	shape = s1; 
	scale = s2; 
    } else {
	scale = s2 / s1;    /* variance / mean */
	shape = s1 / scale; /* mean / scale */
    }

    p = gdtrc(1.0 / scale, shape, x);
    if (get_cephes_errno()) {
	p = NADBL;
    }

    return p;
}

/**
 * gamma_cdf_inverse:
 * @shape: shape.
 * @scale: scale.
 * @p: probability.
 *
 * Returns: the argument x such that the integral from minus infinity
 * to @x of the gamma density with given scale and shape parameters is
 * equal to the given probability @p, or #NADBL on failure. Note that
 * the alternate parametrization (mean, variance) is not supported.
 */

double gamma_cdf_inverse (double shape, double scale, double p)
{
    double x = NADBL;

    if (p==0) {
	return 0;
    }

    if (shape > 0 && scale > 0 && p > 0 && p < 1) {
	x = igami(shape, 1-p) * scale;
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

static int gamma_pdf_array (double shape, double scale, 
			    double *x, int n)
{
    int i, err = 0;

    errno = 0;

    if (!na(shape) && shape > 0 && !na(scale) && scale > 0) {
	double x1, x2;
	double x3 = pow(scale, shape);
	double x4 = gamma_function(shape);

	if (errno || na(x4)) {
	    err = E_NAN;
	} else {
	    for (i=0; i<n; i++) {
		if (!na(x[i]) && x[i] > 0) {
		    errno = 0;
		    x1 = pow(x[i], shape - 1.0);
		    x2 = exp(-x[i]/scale);
		    x[i] = x1 * x2 / (x3 * x4);
		    if (errno) {
			x[i] = NADBL;
		    }
		} else {
		    x[i] = NADBL;
		}
	    }
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    errno = 0;

    return err;
}

/**
 * gamma_pdf:
 * @shape: shape parameter.
 * @scale: scale parameter.
 * @x: reference value.
 *
 * Returns: the value of the gamma pdf at @x, or #NADBL on failure.
 */

double gamma_pdf (double shape, double scale, double x)
{
    gamma_pdf_array(shape, scale, &x, 1);

    return x;
}

/**
 * poisson_cdf:
 * @lambda: mean (also variance).
 * @k: test value.
 *
 * Returns: the probability of X <= @k, for X an r.v. that follows
 * the Poisson distribution with parameter @lambda.
 */

static double poisson_cdf (double lambda, int k)
{
    double x = NADBL;

    if (lambda >= 0 && k >= 0) {
	x = pdtr(k, lambda);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

/**
 * poisson_cdf_comp:
 * @lambda: mean (also variance).
 * @k: test value.
 *
 * Returns: the probability of X > @k, for X an r.v. that follows
 * the Poisson distribution with parameter @lambda.
 */

static double poisson_cdf_comp (double lambda, int k)
{
    double x = NADBL;

    if (lambda >= 0 && k >= 0) {
	x = pdtrc(k, lambda);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

static int poisson_pmf_array (double lambda, double *x, int n)
{
    double den, l0, p;
    int i, j, k;

    if (lambda <= 0) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
	return 0;
    }

    l0 = exp(-lambda);

    for (i=0; i<n; i++) {
	k = x[i];
	den = x_factorial((double) k);
	if (xna(den)) {
	    p = NADBL;
	} else {
	    p = l0 * pow(lambda, (double) k) / den;
	}
	if (xna(p)) {
	    p = l0;
	    for (j=1; j<=k; j++) {
		p *= lambda / j;
	    }
	}
	x[i] = p;
    }	    

    return 0;
}

/**
 * poisson_pmf:
 * @lambda: mean (also variance).
 * @k: test value.
 *
 * Returns: the probability mass at @k, for an r.v. that follows
 * the Poisson distribution with parameter @lambda.
 */

double poisson_pmf (double lambda, int k)
{
    double den, l0, p;

    if (lambda <= 0 || k < 0) {
	return NADBL;
    }

    den = x_factorial((double) k);
    l0 = exp(-lambda);

    if (xna(den)) {
	p = NADBL;
    } else {
	p = l0 * pow(lambda, (double) k) / den;
    }

    if (xna(p)) {
	int i;

	p = l0;
	for (i=1; i<=k; i++) {
	    p *= lambda / i;
	}
    } 

    return p;
}

/* The following is probably horribly inefficient --
   investigate the possibility of using Chebyshev,
   or perhaps binary search.
*/

static double poisson_critval (double mu, double a)
{
    double pk = 0.0;
    double ac = 1 - a;
    int k, k0 = 0;

    if (mu <= 0 || a <= 0 || a >= 1) {
	return NADBL;
    }

    if (mu >= 10 && a < 0.5) {
	k0 = mu - 1;
	pk = poisson_cdf(mu, k0++);
    }     

    for (k=k0; ; k++) {
	pk = poisson_cdf(mu, k);
	if (pk >= ac) {
	    break;
	}
    }  

    return (double) k;
}

/**
 * poisson_cdf_inverse:
 * @k: test value.
 * @p: cumulative probability.
 *
 * Returns: the Poisson parameter such that the integral
 * from 0 to @k of the Poisson density is equal to the
 * given probability @p.
 */

static double poisson_cdf_inverse (int k, double p)
{
    double x = NADBL;

    if (k >= 0 && p >= 0 && p <= 1) {
	x = pdtri(k, p);
	if (get_cephes_errno()) {
	    x = NADBL;
	}
    }

    return x;
}

static double 
weibull_critval (double shape, double scale, double rtail)
{
    double ret = NADBL;

    if (shape > 0 && scale > 0 && rtail > 0 && rtail < 1) {
	ret = scale * pow(-log(rtail), 1.0 / shape);
    }

    return ret;
}

/**
 * weibull_cdf:
 * @shape: shape parameter > 0.
 * @scale: scale parameter > 0.
 * @x: test value.
 *
 * Returns: the probability of X <= @x, for X an r.v. that follows
 * the Weibull distribution with parameters @shape and @scale.
 */

double weibull_cdf (double shape, double scale, double x)
{
    double ret = NADBL;

    if (shape > 0 && scale > 0) {
	if (x == 0.0) {
	    ret = 0.0;
	} else if (x > 0.0) {
	    ret = 1.0 - exp(-pow(x/scale, shape));
	}
    }

    return ret;
}

static double weibull_cdf_comp (double shape, double scale, double x)
{
    double ret = NADBL;

    if (shape > 0 && scale > 0) {
	if (x == 0.0) {
	    ret = 1.0;
	} else if (x > 0.0) {
	    ret = exp(-pow(x/scale, shape));
	}
    }

    return ret;
}

/**
 * GED_pdf:
 * @nu: shape parameter.
 * @x: double.
 * 
 * Returns: the density function of the Generalized Error distribution
 * with shape parameter @nu at @x, or #NADBL on failure.
 */

double GED_pdf (double nu, double x)
{
    if (nu > 0) {
	double lg1 = ln_gamma(1/nu);
	double lg3 = ln_gamma(3/nu);
	double lC  = 0.5*(lg3 - 3*lg1);
	double k   = pow(0.5, 1/nu) * exp(0.5*(lg1 - lg3));
	double znu = pow(fabs(x/k), nu);

	return (0.5 * nu) * exp(lC - 0.5 * znu);
    } else {
	return NADBL;
    }
}

static int GED_pdf_array (double nu, double *x, int n)
{
    int i, err = 0;

    if (nu > 0) {
	double lg1 = ln_gamma(1/nu);
	double lg3 = ln_gamma(3/nu);
	double lC  = 0.5*(lg3 - 3*lg1);
	double k   = pow(0.5, 1/nu) * exp(0.5*(lg1 - lg3));
	double znu; 

	for (i=0; i<n; i++) {
	    if (!na(x[i])) {
		znu = pow(fabs(x[i]/k), nu);
		x[i] = (0.5 * nu) * exp(lC - 0.5 * znu);
	    } else {
		x[i] = NADBL;
	    }
	}
    } else {
	err = E_DATA;
    }

    if (err) {
	for (i=0; i<n; i++) {
	    x[i] = NADBL;
	}
    }

    return err;
}

/**
 * GED_cdf:
 * @nu: shape parameter.
 * @x: reference value.
 *
 * Calculates the value of the CDF of the Generalized Error
 * distribution with parameter @nu at @x. We exploit the fact that
 * if x ~ GED(n), then |x/k|^n is a Gamma rv.
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double GED_cdf (double nu, double x)
{
    if (nu > 0) {
	int sgn    = (x > 0)? 1 : -1;
	double p   = 1/nu;
	double lg1 = ln_gamma(p);
	double lg3 = ln_gamma(3*p);
	double k   = pow(0.5, p) * exp(0.5*(lg1 - lg3));
	double znu = pow(fabs(x/k), nu);
	double P   = gamma_cdf(p, 2, znu, 1);

	P = 0.5 * (1 + sgn*P);
	return P;
    } else {
	return NADBL;
    }
}

/**
 * GED_cdf_comp:
 * @nu: shape parameter.
 * @x: reference value.
 *
 * Calculates the complement of the CDF of the Generalized Error
 * distribution with parameter @nu at @x. We exploit the fact that
 * if x ~ GED(n), then |x/k|^n is a Gamma rv.
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double GED_cdf_comp (double nu, double x)
{
    if (nu > 0) {
	int sgn    = (x > 0)? 1 : -1;
	double p   = 1/nu;
	double lg1 = ln_gamma(p);
	double lg3 = ln_gamma(3*p);
	double k   = pow(0.5, p) * exp(0.5*(lg1 - lg3));
	double znu = pow(fabs(x/k), nu);
	double P   = gamma_cdf_comp(p, 2, znu, 1);

	P = (sgn == 1) ? 0.5 * P : 1 - 0.5 * P;
	return P;
    } else {
	return NADBL;
    }
}

/**
 * johansen_trace_pval:
 * @N: the number of potentially cointegrated variables
 * minus the cointegration rank under H0.
 * @det: index of deterministic setup of model (0 to 4).
 * @T: the sample size, or 0 for an asymptotic result.
 * @tr: the trace statistic.
 *
 * Returns: the p-value of the trace statistic, computed
 * via Doornik's gamma approximation, or #NADBL on failure.
 */

double johansen_trace_pval (int N, int det, int T, double tr)
{
    double (*pvfunc) (double, int, int, int);
    void *handle;
    double pv = NADBL;

    pvfunc = get_plugin_function("trace_pvalue", &handle);

    if (pvfunc == NULL) {
	fputs(I_("Couldn't load plugin function\n"), stderr);
    } else {
	pv = (*pvfunc) (tr, N, det, T);
	close_plugin(handle);
    }

    return pv;
}

/**
 * print_critval:
 * @st: distribution code.
 * @parm: array holding 1 or 2 parameter values.
 * @a: alpha.
 * @c: the critical value.
 * @prn: gretl printer.
 *
 * Prints the critical value information in a consistent manner.
 */

void print_critval (char st, const double *parm, double a, double c, PRN *prn)
{
    switch (st) {
    case 'z':
	pprintf(prn, "%s", _("Standard normal distribution"));
	break;
    case 't':
	pprintf(prn, "t(%g)", parm[0]);
	break;
    case 'X':
	pprintf(prn, "%s(%g)", _("Chi-square"), parm[0]);
	break;
    case 'F':
	pprintf(prn, "F(%g, %g)", parm[0], parm[1]);
	break;
    case 'B':
	pprintf(prn, "Binomial (P = %g, %g trials)", parm[0], parm[1]);
	break;
    case 'P':
	pprintf(prn, "Poisson (mean = %g)", parm[0]);
	break;
    case 'W':
	pprintf(prn, "Weibull (shape = %g, scale = %g)", parm[0], parm[1]);
	break;
    }

    pputs(prn, "\n ");
    pprintf(prn, _("right-tail probability = %g"), a);
    pputs(prn, "\n ");
    pprintf(prn, _("complementary probability = %g"), 1.0 - a);
    if (a < 0.5 && (st == 'z' || st == 't')) {
	pputs(prn, "\n ");
	pprintf(prn, _("two-tailed probability = %g"), 2.0 * a);
    }
    pputs(prn, "\n\n ");
    pprintf(prn, _("Critical value = %g"), c);
    pputc(prn, '\n');    
}

/* This apparatus is for use with the "batch p-value" 
   routine: it remembers the parameters (p) and
   argument (x) from the last internal p-value 
   assessment.
*/

static double pvargs[3];

static void remember_pvalue_args (const double *p, double x)
{
    pvargs[0] = p[0];
    pvargs[1] = p[1];
    pvargs[2] = x;
}

/* end remember parameters */

static int pdist_check_input (char st, const double *parm,
			      double x)
{
    int i, np = 1;

    if (na(x)) {
	return E_MISSDATA;
    }

    if (st == 'z') {
	np = 0;
    } else if (st == 'F' || st == 'G' ||
	       st == 'B' || st == 'W') {
	np = 2;
    } else if (st == 'J') {
	np = 3;
    }

    for (i=0; i<np; i++) {
	if (na(parm[i])) {
	    return E_MISSDATA;
	}
    }

    return 0;
}

/**
 * gretl_get_cdf_inverse:
 * @st: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @a: probability value.
 *
 * Returns: the argument, y, for which the area under the PDF
 * specified by @st and @parm, integrated from its minimum to y,
 * is equal to @a, or #NADBL on failure.
 */

double gretl_get_cdf_inverse (char st, const double *parm, 
			      double a)
{
    double y = NADBL;

    if (pdist_check_input(st, parm, a) == E_MISSDATA) {
	return y;
    }

    if (st == 'z') {
	y = normal_cdf_inverse(a);
    } else if (st == 't') {
	y = student_cdf_inverse(parm[0], a);
    } else if (st == 'X') {
	y = chisq_cdf_inverse((int) parm[0], a);
    } else if (st == 'G') {
	y = gamma_cdf_inverse(parm[0], parm[1], a);
    } else if (st == 'F') {
	y = snedecor_cdf_inverse((int) parm[0], (int) parm[1], a);
    } else if (st == 'B') {
	y = binomial_cdf_inverse((int) parm[0], (int) parm[1], a);
    } else if (st == 'P') {
	y = poisson_cdf_inverse((int) parm[0], a);
    } 

    return y;
}

/**
 * gretl_get_critval:
 * @st: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @a: right-tail probability.
 *
 * Returns: the abcsissa value x for the distribution specified
 * by @st and @parm, such that P(X >= x) = @a, or #NADBL on 
 * failure.
 */

double gretl_get_critval (char st, const double *parm, double a)
{
    double x = NADBL;

    if (pdist_check_input(st, parm, a) == E_MISSDATA) {
	return x;
    }

    if (st == 'z') {
	x = normal_critval(a);
    } else if (st == 't') {
	x = student_critval(parm[0], a);
    } else if (st == 'X') {	
	x = chisq_critval((int) parm[0], a);
    } else if (st == 'F') {
	x = snedecor_critval((int) parm[0], (int) parm[1], a);
    } else if (st == 'B') {
	x = binomial_critval(parm[0], (int) parm[1], a);
    } else if (st == 'P') {
	x = poisson_critval(parm[0], a);
    } else if (st == 'W') {
	x = weibull_critval(parm[0], parm[1], a);
    } 

    return x;
}

/**
 * gretl_get_cdf:
 * @st: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Evaluates the CDF for the distribution specified by
 * @st and @parm applicable at @x. 
 *
 * Returns: the CDF value, or #NADBL on error.
 */

double gretl_get_cdf (char st, const double *parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(st, parm, x) == E_MISSDATA) {
	return y;
    }    

    if (st == 'z') {
	y = normal_cdf(x);
    } else if (st == 't') {
	y = student_cdf(parm[0], x);
    } else if (st == 'X') {
	y = chisq_cdf((int) parm[0], x);
    } else if (st == 'F') {
	y = snedecor_cdf((int) parm[0], (int) parm[1], x);
    } else if (st == 'G') {
	y = gamma_cdf(parm[0], parm[1], x, 1);
    } else if (st == 'B') {
	y = binomial_cdf(parm[0], (int) parm[1], (int) x);
    } else if (st == 'P') {
	y = poisson_cdf(parm[0], (int) x);
    } else if (st == 'W') {
	y = weibull_cdf(parm[0], parm[1], x);
    } else if (st == 'E') {
	y = GED_cdf(parm[0], x);
    }

    return y;
}

/**
 * gretl_get_pdf:
 * @st: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Evaluates the PDF for the distribution specified by
 * @st and @parm at @x. 
 *
 * Returns: the PDF value, or #NADBL on error.
 */

double gretl_get_pdf (char st, const double *parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(st, parm, x) == E_MISSDATA) {
	return y;
    } 

    if (st == 'z') {
	y = normal_pdf(x);
    } else if (st == 't') {
	y = student_pdf(parm[0], x);
    } else if (st == 'X') {
	y = chisq_pdf((int) parm[0], x);
    } else if (st == 'F') {
	y = snedecor_pdf((int) parm[0], (int) parm[1], x);
    } else if (st == 'G') {
	y = gamma_pdf(parm[0], parm[1], x);
    } else if (st == 'B') {
	y = binomial_pmf(parm[0], parm[1], x);
    } else if (st == 'P') {
	y = poisson_pmf(parm[0], x);
    } else if (st == 'W') {
	y = weibull_pdf(parm[0], parm[1], x);
    } else if (st == 'E') {
	y = GED_pdf(parm[0], x);
    }

    return y;
}

/**
 * gretl_fill_pdf_array:
 * @st: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @x: see below.
 * @n: number of elements in @x.
 *
 * On input, @x contains an array of abscissae at which the
 * PDF specified by @st and @parm should be evaluated.  On
 * output it contains the corresponding PDF values.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int gretl_fill_pdf_array (char st, const double *parm, 
			  double *x, int n)
{
    int err = E_DATA;

    if (pdist_check_input(st, parm, 0) == E_MISSDATA) {
	return E_MISSDATA;
    } 

    if (st == 'z') {
	err = normal_pdf_array(x, n);
    } else if (st == 't') {
	err = student_pdf_array(parm[0], x, n);
    } else if (st == 'X') {
	err = chisq_pdf_array(parm[0], x, n);
    } else if (st == 'F') {
	err = snedecor_pdf_array((int) parm[0], (int) parm[1], x, n);
    } else if (st == 'G') {
	err = gamma_pdf_array(parm[0], parm[1], x, n); 
    } else if (st == 'B') {
	err = binomial_pmf_array(parm[0], parm[1], x, n);
    } else if (st == 'P') {
	err = poisson_pmf_array(parm[0], x, n);
    } else if (st == 'W') {
	err = weibull_pdf_array(parm[0], parm[1], x, n);
    } else if (st == 'E') {
	err = GED_pdf_array(parm[0], x, n);
    }

    return err;
}

/**
 * gretl_get_pvalue:
 * @st: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Returns: the integral of the PDF specified by @st and
 * @parm from @x to infinity, or #NADBL on error.
 */

double gretl_get_pvalue (char st, const double *parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(st, parm, x) == E_MISSDATA) {
	return y;
    } 

    if (st == 'z') {
	y = normal_cdf_comp(x);
    } else if (st == 't') {
	y = student_cdf_comp(parm[0], x);
    } else if (st == 'X') {
	y = chisq_cdf_comp((int) parm[0], x);
    } else if (st == 'F') {
	y = snedecor_cdf_comp((int) parm[0], (int) parm[1], x);
    } else if (st == 'G') {
	y = gamma_cdf_comp(parm[0], parm[1], x, 1);
    } else if (st == 'B') {
	y = binomial_cdf_comp(parm[0], (int) parm[1], x);
    } else if (st == 'P') {
	y = poisson_cdf_comp(parm[0], x);
    } else if (st == 'W') {
	y = weibull_cdf_comp(parm[0], parm[1], x);
    } else if (st == 'E') {
	y = GED_cdf_comp(parm[0], x);
    } else if (st == 'J') {
	y = johansen_trace_pval((int) parm[0], (int) parm[1], 
				(int) parm[2], x);
    }

    remember_pvalue_args(parm, x);

    return y;
}

/**
 * gretl_get_random_series:
 * @st: distribution code.
 * @parm: array holding either one or two scalar 
 * parameter values, depending on the distribution.
 * @serp1: series containing values for first param,
 * or %NULL.
 * @serp2: series containing values for second param,
 * or %NULL.
 * @pdinfo: dataset information.
 * @err: location to receive error code.
 *
 * Produces a random series conforming to the distribution
 * given by @st, which may require specification of either
 * one or two parameters. These parameters are either
 * given as (scalar) elements of @parm, or they may vary
 * by observation, in which case they are given as the
 * elements of @serp1 or @serp2.
 *
 * Returns: the array of pseudo-random values, or %NULL
 * on error.
 */

double *gretl_get_random_series (char st, const double *parm,
				 const double *serp1, 
				 const double *serp2,
				 const DATAINFO *pdinfo,
				 int *err)
{
    double *x = malloc(pdinfo->n * sizeof *x);
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    int t;

    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (t=0; t<pdinfo->n; t++) {
	x[t] = NADBL;
    }

    if (st == 'u') {
	/* uniform */
	double min = parm[0], max = parm[1];

	if (serp1 != NULL || serp2 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) min = serp1[t];
		if (serp2 != NULL) max = serp2[t];
		*err = gretl_rand_uniform_minmax(x, t, t, min, max);
	    }
	} else {
	    *err = gretl_rand_uniform_minmax(x, t1, t2, min, max);
	}
    } else if (st == 'z') {
	/* normal */
	double mu = parm[0], sd = parm[1];

	if (serp1 != NULL || serp2 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) mu = serp1[t];
		if (serp2 != NULL) sd = serp2[t];
		*err = gretl_rand_normal_full(x, t, t, mu, sd);
	    }
	} else {
	    *err = gretl_rand_normal_full(x, t1, t2, mu, sd);
	}
    } else if (st == 't') {
	/* Student's t */
	int v = parm[0];

	if (serp1 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		v = serp1[t];
		*err = gretl_rand_student(x, t, t, v);
	    }
	} else {	
	    *err = gretl_rand_student(x, t1, t2, v);
	}
    } else if (st == 'X') {
	/* chi-square */
	int v = parm[0];

	if (serp1 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		v = serp1[t];
		*err = gretl_rand_chisq(x, t, t, v);
	    }
	} else {
	    *err = gretl_rand_chisq(x, t1, t2, v);
	}
    } else if (st == 'F') {
	/* Snedecor F */
	int v1 = parm[0], v2 = parm[1];

	if (serp1 != NULL || serp2 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) v1 = serp1[t];
		if (serp2 != NULL) v2 = serp2[t];
		*err = gretl_rand_F(x, t, t, v1, v2);
	    }
	} else {
	    *err = gretl_rand_F(x, t1, t2, v1, v2);
	}
    } else if (st == 'G') {
	/* gamma */
	double shape = parm[0], scale = parm[1];

	if (serp1 != NULL || serp2 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) shape = serp1[t];
		if (serp2 != NULL) scale = serp2[t];
		*err = gretl_rand_gamma(x, t, t, shape, scale);
	    }
	} else {
	    *err = gretl_rand_gamma(x, t1, t2, shape, scale);
	}
    } else if (st == 'B') {
	/* binomial */
	double pr = parm[0];
	int n = parm[1];

	if (serp1 != NULL || serp2 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) pr = serp1[t];
		if (serp2 != NULL) n = serp2[1];
		*err = gretl_rand_binomial(x, t, t, n, pr);
	    }
	} else {
	    *err = gretl_rand_binomial(x, t1, t2, n, pr);
	}
    } else if (st == 'P') {
	/* Poisson */
	double m = parm[0];

	if (serp1 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		m = serp1[t];
		*err = gretl_rand_poisson(x, t, t, &m, 0);
	    }
	} else {
	    *err = gretl_rand_poisson(x, t1, t2, &m, 0);
	}
    } else if (st == 'W') {
	/* Weibull */ 
	double shape = parm[0], scale = parm[1];

	if (serp1 != NULL || serp2 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) shape = serp1[t];
		if (serp2 != NULL) scale = serp2[t];
		*err = gretl_rand_weibull(x, t, t, shape, scale);
	    }
	} else {
	    *err = gretl_rand_weibull(x, t1, t2, shape, scale);
	}
    } else if (st == 'E') {
	/* GED */ 
	double nu = parm[0];

	if (serp1 != NULL) {
	    for (t=t1; t<=t2 && !*err; t++) {
		if (serp1 != NULL) nu = serp1[t];
		*err = gretl_rand_GED(x, t, t, nu);
	    }
	} else {
	    *err = gretl_rand_GED(x, t1, t2, nu);
	}
    }	

    return x;
}

static int 
print_pv_string (double x, double p, PRN *prn)
{
    char numstr[32];

    if (na(p)) {
	pprintf(prn, _("area to the right of %g: NA\n"), x);
	return 1;
    }

    sprintf(numstr, "%g", p);

    if (!strcmp(numstr, "1") || !strcmp(numstr, "0")) {
	pprintf(prn, _("area to the right of %g =~ %g\n"), x, p);
    } else {
	pprintf(prn, _("area to the right of %g = %g\n"), x, p);
    }

    return 0;
}

/**
 * print_pvalue:
 * @st: distribution code.
 * @parm: array holding 1 or 2 parameter values.
 * @x: the value in the distribution.
 * @pv: the p-value.
 * @prn: gretl printer.
 *
 * Prints the p-value information in a consistent manner.
 */

void print_pvalue (char st, const double *parm, double x,
		   double pv, PRN *prn)
{
    double pc;
    int err;

    switch (st) {

    case 'z':
    case 'n':
    case 'N':
    case '1':
	pprintf(prn, "%s: ", _("Standard normal"));
	err = print_pv_string(x, pv, prn);
	if (err) return;
	if (pv < 0.5) {
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pv, 1 - 2 * pv);
	} else {
	    pc = normal_cdf(x);
	    pprintf(prn, _("(to the left: %g)\n"), pc);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pc, 1 - 2 * pc);
	}
	break;

    case 't':
    case '2':
	pprintf(prn, "t(%d): ", (int) parm[0]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	if (pv < 0.5) {
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pv, 1 - 2 * pv);
	} else {
	    pc = student_cdf(parm[0], x);
	    pprintf(prn, _("(to the left: %g)\n"), pc);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pc, 1 - 2 * pc);
	}
	break;

    case 'X':
    case 'x':
    case 'c':
    case '3':
	pprintf(prn, "%s(%d): ", _("Chi-square"), (int) parm[0]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = chisq_cdf(parm[0], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'F':
    case 'f':
    case '4':
	pprintf(prn, "F(%d, %d): ", (int) parm[0], (int) parm[1]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = snedecor_cdf((int) parm[0], (int) parm[1], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'G':
    case 'g':
    case '5':
	pprintf(prn, _("Gamma (shape %g, scale %g, mean %g, variance %g):"
		       "\n area to the right of %g = %g\n"), 
		parm[0], parm[1], parm[0] * parm[1], 
		parm[0] * parm[1] * parm[1],
		x, pv);
	break;

    case 'B':
    case 'b':
    case '6':
	pprintf(prn, _("Binomial (p = %g, n = %d):"
		       "\n Prob(x > %d) = %g\n"), 
		parm[0], (int) parm[1], (int) x, pv);
	pc = binomial_cdf(parm[0], parm[1], x);
	if (x > 0) {
	    pprintf(prn, _(" Prob(x <= %d) = %g\n"), (int) x, pc);
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) x,
		    pc - binomial_cdf(parm[0], parm[1], x - 1));
	} else {
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) x, pc);
	}
	break;

    case 'p':
    case 'P':
    case '8':
	pprintf(prn, _("Poisson (mean = %g): "), parm[0]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = poisson_cdf(parm[0], (int) x);
	if (x > 0) {
	    pprintf(prn, _(" Prob(x <= %d) = %g\n"), (int) x, pc);
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) x,
		    poisson_pmf(parm[0], x));
	} else {
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) x, pc);
	}
	break;	

    case 'w':
    case 'W':
    case '9':
	pprintf(prn, _("Weibull (shape = %g, scale = %g): "), 
		parm[0], parm[1]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = weibull_cdf(parm[0], parm[1], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'E':
	pprintf(prn, _("GED (shape = %g): "), parm[0]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = GED_cdf(parm[0], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    default:
	break;
    }
}

/**
 * batch_pvalue:
 * @str: the command line, which should be of one of the following forms:
 * pvalue z x (Normal distribution);
 * pvalue t df x (t-distribution);
 * pvalue X df x (Chi-square);
 * pvalue F dfn dfd x (F-distribution); or
 * pvalue G mean variance x (Gamma distribution).
 * pvalue B prob n x (Binomial distribution).
 * pvalue P mean k (Poisson distribution).
 * pvalue W shape scale x (Weibull distribution).
 * @pZ: pointer to the data array.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Calculates and prints the probability that a random variable 
 * distributed as specified in the command line @str exceeds the 
 * value indicated in @str.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int batch_pvalue (const char *str, 
		  double ***pZ, DATAINFO *pdinfo, 
		  PRN *prn)
{
    double pv = NADBL;
    char line[MAXLEN];
    char **S;
    char st;
    int i, n, m;
    int err = 0;
    
    if (!strncmp(str, "pvalue ", 7)) {
	str += 7;
    }

    S = gretl_string_split(str, &n);
    if (S == NULL) {
	return E_ALLOC;
    }

    st = S[0][0];

    strcpy(line, "pvalue(");
    m = 8;
    for (i=0; i<n && !err; i++) {
	m += strlen(S[i]) + 1;
	if (m > MAXLEN) {
	    err = E_DATA;
	} else {
	    strcat(line, S[i]);
	    strcat(line, (i == n - 1)? ")" : ",");
	}
    }

    free_strings_array(S, n);

    if (!err) {
	pv = generate_scalar(line, pZ, pdinfo, &err);
    }

    if (!err) {
	print_pvalue(st, pvargs, pvargs[2], pv, prn);
    }

    return err;
}

/**
 * gretl_get_DW:
 * @n: number of observations
 * @k: number of regressors excluding the constant.
 * @err: location to receive error code.
 *
 * Consults a table of Durbin-Watson critical values and
 * returns the results in a gretl_matrix.  
 *
 * Returns: on success, a 4-vector containing the lower
 * and upper Durbin-Watson values, dl and du, along with
 * the values actually used for @n and @k (which may differ
 * from those given on input if the exact values are not
 * found in the table and have to be approximated).  
 * On error, returns %NULL.
 */

gretl_matrix *gretl_get_DW (int n, int k, int *err)
{
    void *handle;
    int (*dw_lookup) (int, int, gretl_matrix **);
    gretl_matrix *m = NULL;

    dw_lookup = get_plugin_function("dw_lookup", &handle);

    if (dw_lookup == NULL) {
	*err = E_FOPEN;
	return NULL;
    }

    *err = (*dw_lookup) (n, k, &m);
    close_plugin(handle);

    return m;
}

