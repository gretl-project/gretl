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
#include "libset.h"
#include "../../cephes/libprob.h"

#include <errno.h>

#define NORM_CDF_MAX 0.9999999999999999

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

double chisq_cdf (double df, double x)
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

double chisq_cdf_comp (double df, double x)
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
	y = NORM_CDF_MAX;
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
 * Returns: the argument x such that the integral from zero to @x of
 * the gamma density with given scale and shape parameters is equal to
 * the given probability @p, or #NADBL on failure. Note that the
 * alternate parametrization (mean, variance) is not supported.
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
 * GED_cdf_inverse:
 * @nu: shape parameter.
 * @a: probability.
 *
 * Returns: the argument x such that the integral from minus infinity
 * to @x of the standardized GED density with shape parameter @nu is
 * equal to the given probability @a, or #NADBL on failure. We exploit
 * the well-known relationship between the standardized GED and the
 * Gamma variates.
 */

double GED_cdf_inverse (double nu, double a)
{
    if (nu > 0 && a < 1 && a > 0) {
	double a2, p, lg1, lg3, sd, x;
	int sgn;

	if (a > 0.5) {
	    a2 = 2*a - 1;
	    sgn = 1;
	} else {
	    a2 = 1 - 2*a;
	    sgn = -1;
	}

	p   = 1/nu;
	lg1 = ln_gamma(p);
	lg3 = ln_gamma(3*p);
	sd  = pow(2.0, p) * exp(0.5*(lg3 - lg1));
	x   = gamma_cdf_inverse(p, 2.0, a2);

	return sgn * pow(x, p) / sd;
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

struct distmap {
    int code;
    char *s;
};

int dist_code_from_string (const char *s)
{
    struct distmap dmap[] = {
	{ D_UNIFORM,  "u" },
	{ D_UDISCRT,  "i" },
	{ D_NORMAL,   "z" },
	{ D_STUDENT,  "t" },
	{ D_CHISQ,    "x" },
	{ D_SNEDECOR, "f" },
	{ D_BINOMIAL, "b" },
	{ D_POISSON,  "p" },
	{ D_WEIBULL,  "w" },
	{ D_GAMMA,    "g" },
	{ D_GED,      "e" },
	{ D_BETA,     "beta" },
	{ D_DW,       "d" },
	{ D_BINORM,   "D" },
	{ D_JOHANSEN, "J" },
	{ D_BETABIN,  "bb" },
	{ D_NONE,     NULL }
    };
    char test[8];
    int i;

    if (!strcmp(s, "D")) {
	/* special: case counts for bivariate normal */
	return D_BINORM;
    }

    /* otherwise we'll ignore case */
    for (i=0; i<8 && s[i]; i++) {
	test[i] = tolower(s[i]);
    }
    test[i] = '\0';

    for (i=0; dmap[i].code; i++) {
	if (!strcmp(test, dmap[i].s)) {
	    return dmap[i].code;
	}
    }

    /* backward compatibility (do we need this?) */
    if (!strcmp(test, "n")) {
	return D_NORMAL;
    } else if (!strcmp(test, "c")) {
	return D_CHISQ;
    }

    return D_NONE;
}

/**
 * print_critval:
 * @dist: distribution code.
 * @parm: array holding 0 to 2 parameter values.
 * @a: alpha.
 * @c: the critical value.
 * @prn: gretl printer.
 *
 * Prints the critical value information in a consistent manner.
 */

void print_critval (int dist, const double *parm, double a, double c, PRN *prn)
{
    switch (dist) {
    case D_NORMAL:
	pprintf(prn, "%s", _("Standard normal distribution"));
	break;
    case D_STUDENT:
	pprintf(prn, "t(%g)", parm[0]);
	break;
    case D_CHISQ:
	pprintf(prn, "%s(%g)", _("Chi-square"), parm[0]);
	break;
    case D_SNEDECOR:
	pprintf(prn, "F(%g, %g)", parm[0], parm[1]);
	break;
    case D_BINOMIAL:
	pprintf(prn, "Binomial (P = %g, %g trials)", parm[0], parm[1]);
	break;
    case D_POISSON:
	pprintf(prn, "Poisson (mean = %g)", parm[0]);
	break;
    case D_WEIBULL:
	pprintf(prn, "Weibull (shape = %g, scale = %g)", parm[0], parm[1]);
	break;
    }

    pputs(prn, "\n ");
    pprintf(prn, _("right-tail probability = %g"), a);
    pputs(prn, "\n ");
    pprintf(prn, _("complementary probability = %g"), 1.0 - a);
    if (a < 0.5 && (dist == D_NORMAL || dist == D_STUDENT)) {
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

static int pdist_check_input (int dist, const double *parm,
			      double x)
{
    int i, np = 1;

    if (na(x)) {
	return E_MISSDATA;
    }

    if (dist == D_NORMAL) {
	np = 0;
    } else if (dist == D_SNEDECOR || dist == D_GAMMA ||
	       dist == D_BINOMIAL || dist == D_WEIBULL) {
	np = 2;
    } else if (dist == D_JOHANSEN || dist == D_BETABIN) {
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
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @a: probability value.
 *
 * Returns: the argument, y, for which the area under the PDF
 * specified by @dist and @parm, integrated from its minimum to y,
 * is equal to @a, or #NADBL on failure.
 */

double gretl_get_cdf_inverse (int dist, const double *parm, 
			      double a)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, a) == E_MISSDATA) {
	return y;
    }

    if (dist == D_NORMAL) {
	y = normal_cdf_inverse(a);
    } else if (dist == D_STUDENT) {
	y = student_cdf_inverse(parm[0], a);
    } else if (dist == D_CHISQ) {
	y = chisq_cdf_inverse((int) parm[0], a);
    } else if (dist == D_GAMMA) {
	y = gamma_cdf_inverse(parm[0], parm[1], a);
    } else if (dist == D_SNEDECOR) {
	y = snedecor_cdf_inverse((int) parm[0], (int) parm[1], a);
    } else if (dist == D_BINOMIAL) {
	y = binomial_cdf_inverse((int) parm[0], (int) parm[1], a);
    } else if (dist == D_POISSON) {
	y = poisson_cdf_inverse((int) parm[0], a);
    } else if (dist == D_GED) {
	y = GED_cdf_inverse(parm[0], a);
    } 

    return y;
}

/**
 * gretl_get_critval:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @a: right-tail probability.
 *
 * Returns: the abcsissa value x for the distribution specified
 * by @dist and @parm, such that P(X >= x) = @a, or #NADBL on 
 * failure.
 */

double gretl_get_critval (int dist, const double *parm, double a)
{
    double x = NADBL;

    if (pdist_check_input(dist, parm, a) == E_MISSDATA) {
	return x;
    }

    if (dist == D_NORMAL) {
	x = normal_critval(a);
    } else if (dist == D_STUDENT) {
	x = student_critval(parm[0], a);
    } else if (dist == D_CHISQ) {	
	x = chisq_critval((int) parm[0], a);
    } else if (dist == D_SNEDECOR) {
	x = snedecor_critval((int) parm[0], (int) parm[1], a);
    } else if (dist == D_BINOMIAL) {
	x = binomial_critval(parm[0], (int) parm[1], a);
    } else if (dist == D_POISSON) {
	x = poisson_critval(parm[0], a);
    } else if (dist == D_WEIBULL) {
	x = weibull_critval(parm[0], parm[1], a);
    } 

    return x;
}

/**
 * gretl_get_cdf:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Evaluates the CDF for the distribution specified by
 * @dist and @parm applicable at @x. 
 *
 * Returns: the CDF value, or #NADBL on error.
 */

double gretl_get_cdf (int dist, const double *parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, x) == E_MISSDATA) {
	return y;
    }    

    if (dist == D_NORMAL) {
	y = normal_cdf(x);
    } else if (dist == D_STUDENT) {
	y = student_cdf(parm[0], x);
    } else if (dist == D_CHISQ) {
	y = chisq_cdf((int) parm[0], x);
    } else if (dist == D_SNEDECOR) {
	y = snedecor_cdf((int) parm[0], (int) parm[1], x);
    } else if (dist == D_GAMMA) {
	y = gamma_cdf(parm[0], parm[1], x, 1);
    } else if (dist == D_BINOMIAL) {
	y = binomial_cdf(parm[0], (int) parm[1], (int) x);
    } else if (dist == D_POISSON) {
	y = poisson_cdf(parm[0], (int) x);
    } else if (dist == D_WEIBULL) {
	y = weibull_cdf(parm[0], parm[1], x);
    } else if (dist == D_GED) {
	y = GED_cdf(parm[0], x);
    }

    return y;
}

/**
 * gretl_get_pdf:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Evaluates the PDF for the distribution specified by
 * @dist and @parm at @x. 
 *
 * Returns: the PDF value, or #NADBL on error.
 */

double gretl_get_pdf (int dist, const double *parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, x) == E_MISSDATA) {
	return y;
    } 

    if (dist == D_NORMAL) {
	y = normal_pdf(x);
    } else if (dist == D_STUDENT) {
	y = student_pdf(parm[0], x);
    } else if (dist == D_CHISQ) {
	y = chisq_pdf((int) parm[0], x);
    } else if (dist == D_SNEDECOR) {
	y = snedecor_pdf((int) parm[0], (int) parm[1], x);
    } else if (dist == D_GAMMA) {
	y = gamma_pdf(parm[0], parm[1], x);
    } else if (dist == D_BINOMIAL) {
	y = binomial_pmf(parm[0], parm[1], x);
    } else if (dist == D_POISSON) {
	y = poisson_pmf(parm[0], x);
    } else if (dist == D_WEIBULL) {
	y = weibull_pdf(parm[0], parm[1], x);
    } else if (dist == D_GED) {
	y = GED_pdf(parm[0], x);
    }

    return y;
}

/**
 * gretl_fill_pdf_array:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @x: see below.
 * @n: number of elements in @x.
 *
 * On input, @x contains an array of abscissae at which the
 * PDF specified by @dist and @parm should be evaluated.  On
 * output it contains the corresponding PDF values.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int gretl_fill_pdf_array (int dist, const double *parm, 
			  double *x, int n)
{
    int err = E_DATA;

    if (pdist_check_input(dist, parm, 0) == E_MISSDATA) {
	return E_MISSDATA;
    } 

    if (dist == D_NORMAL) {
	err = normal_pdf_array(x, n);
    } else if (dist == D_STUDENT) {
	err = student_pdf_array(parm[0], x, n);
    } else if (dist == D_CHISQ) {
	err = chisq_pdf_array(parm[0], x, n);
    } else if (dist == D_SNEDECOR) {
	err = snedecor_pdf_array((int) parm[0], (int) parm[1], x, n);
    } else if (dist == D_GAMMA) {
	err = gamma_pdf_array(parm[0], parm[1], x, n); 
    } else if (dist == D_BINOMIAL) {
	err = binomial_pmf_array(parm[0], parm[1], x, n);
    } else if (dist == D_POISSON) {
	err = poisson_pmf_array(parm[0], x, n);
    } else if (dist == D_WEIBULL) {
	err = weibull_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == D_GED) {
	err = GED_pdf_array(parm[0], x, n);
    }

    return err;
}

/**
 * gretl_get_pvalue:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values, 
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Returns: the integral of the PDF specified by @dist and
 * @parm from @x to infinity, or #NADBL on error.
 */

double gretl_get_pvalue (int dist, const double *parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, x) == E_MISSDATA) {
	return y;
    } 

    if (dist == D_NORMAL) {
	y = normal_cdf_comp(x);
    } else if (dist == D_STUDENT) {
	y = student_cdf_comp(parm[0], x);
    } else if (dist == D_CHISQ) {
	y = chisq_cdf_comp((int) parm[0], x);
    } else if (dist == D_SNEDECOR) {
	y = snedecor_cdf_comp((int) parm[0], (int) parm[1], x);
    } else if (dist == D_GAMMA) {
	y = gamma_cdf_comp(parm[0], parm[1], x, 1);
    } else if (dist == D_BINOMIAL) {
	y = binomial_cdf_comp(parm[0], (int) parm[1], x);
    } else if (dist == D_POISSON) {
	y = poisson_cdf_comp(parm[0], x);
    } else if (dist == D_WEIBULL) {
	y = weibull_cdf_comp(parm[0], parm[1], x);
    } else if (dist == D_GED) {
	y = GED_cdf_comp(parm[0], x);
    } else if (dist == D_JOHANSEN) {
	y = johansen_trace_pval((int) parm[0], (int) parm[1], 
				(int) parm[2], x);
    }

    remember_pvalue_args(parm, x);

    return y;
}

static int gretl_fill_random_array (double *x, int t1, int t2,
				    int dist, const double *parm,
				    const double *vecp1, 
				    const double *vecp2)
{
    int t, err = 0;

    if (dist == D_UNIFORM) {
	/* uniform, continuous */
	double min = parm[0], max = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) min = vecp1[t];
		if (vecp2 != NULL) max = vecp2[t];
		err = gretl_rand_uniform_minmax(x, t, t, min, max);
	    }
	} else {
	    err = gretl_rand_uniform_minmax(x, t1, t2, min, max);
	}
    } else if (dist == D_UDISCRT) {
	/* uniform, discrete */
	int min = (int) parm[0], max = (int) parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) min = (int) vecp1[t];
		if (vecp2 != NULL) max = (int) vecp2[t];
		err = gretl_rand_uniform_int_minmax(x, t, t, min, max,
						    OPT_NONE);
	    }
	} else {
	    err = gretl_rand_uniform_int_minmax(x, t1, t2, min, max,
						OPT_NONE);
	}	
    } else if (dist == D_NORMAL) {
	double mu = parm[0], sd = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) mu = vecp1[t];
		if (vecp2 != NULL) sd = vecp2[t];
		err = gretl_rand_normal_full(x, t, t, mu, sd);
	    }
	} else {
	    err = gretl_rand_normal_full(x, t1, t2, mu, sd);
	}
    } else if (dist == D_STUDENT) {
	/* Student's t */
	double v = parm[0];

	if (vecp1 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		v = vecp1[t];
		err = gretl_rand_student(x, t, t, v);
	    }
	} else {	
	    err = gretl_rand_student(x, t1, t2, v);
	}
    } else if (dist == D_CHISQ) {
	/* chi-square */
	int v = parm[0];

	if (vecp1 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		v = vecp1[t];
		err = gretl_rand_chisq(x, t, t, v);
	    }
	} else {
	    err = gretl_rand_chisq(x, t1, t2, v);
	}
    } else if (dist == D_SNEDECOR) {
	int v1 = parm[0], v2 = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) v1 = vecp1[t];
		if (vecp2 != NULL) v2 = vecp2[t];
		err = gretl_rand_F(x, t, t, v1, v2);
	    }
	} else {
	    err = gretl_rand_F(x, t1, t2, v1, v2);
	}
    } else if (dist == D_GAMMA) {
	double shape = parm[0], scale = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) shape = vecp1[t];
		if (vecp2 != NULL) scale = vecp2[t];
		err = gretl_rand_gamma(x, t, t, shape, scale);
	    }
	} else {
	    err = gretl_rand_gamma(x, t1, t2, shape, scale);
	}
    } else if (dist == D_BINOMIAL) {
	double pr = parm[0];
	int n = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) pr = vecp1[t];
		if (vecp2 != NULL) n = vecp2[1];
		err = gretl_rand_binomial(x, t, t, n, pr);
	    }
	} else {
	    err = gretl_rand_binomial(x, t1, t2, n, pr);
	}
    } else if (dist == D_POISSON) {
	double m = parm[0];

	if (vecp1 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		m = vecp1[t];
		err = gretl_rand_poisson(x, t, t, &m, 0);
	    }
	} else {
	    err = gretl_rand_poisson(x, t1, t2, &m, 0);
	}
    } else if (dist == D_WEIBULL) {
	double shape = parm[0], scale = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) shape = vecp1[t];
		if (vecp2 != NULL) scale = vecp2[t];
		err = gretl_rand_weibull(x, t, t, shape, scale);
	    }
	} else {
	    err = gretl_rand_weibull(x, t1, t2, shape, scale);
	}
    } else if (dist == D_GED) {
	double nu = parm[0];

	if (vecp1 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) nu = vecp1[t];
		err = gretl_rand_GED(x, t, t, nu);
	    }
	} else {
	    err = gretl_rand_GED(x, t1, t2, nu);
	}
    } else if (dist == D_BETA) {
	double shape1 = parm[0], shape2 = parm[1];

	if (vecp1 != NULL || vecp2 != NULL) {
	    for (t=t1; t<=t2 && !err; t++) {
		if (vecp1 != NULL) shape1 = vecp1[t];
		if (vecp2 != NULL) shape2 = vecp2[t];
		err = gretl_rand_beta(x, t, t, shape1, shape2);
	    }
	} else {
	    err = gretl_rand_beta(x, t1, t2, shape1, shape2);
	}
    } else if (dist == D_BETABIN) {
	double shape1 = parm[1], shape2 = parm[2];
	int n = parm[0];

	err = gretl_rand_beta_binomial(x, t1, t2, n, shape1, shape2);
    }

    return err;
}

/**
 * gretl_get_random_series:
 * @dist: distribution code.
 * @parm: array holding either one or two scalar 
 * parameter values, depending on the distribution.
 * @vecp1: series containing values for first param,
 * or %NULL.
 * @vecp2: series containing values for second param,
 * or %NULL.
 * @dset: dataset information.
 * @err: location to receive error code.
 *
 * Produces a random series conforming to the distribution
 * given by @dist, which may require specification of either
 * one or two parameters. These parameters are either
 * given as (scalar) elements of @parm, or they may vary
 * by observation, in which case they are given as the
 * elements of @serp1 or @serp2.
 *
 * Returns: the array of pseudo-random values, or %NULL
 * on error.
 */

double *gretl_get_random_series (int dist, const double *parm,
				 const double *vecp1, 
				 const double *vecp2,
				 const DATASET *dset,
				 int *err)
{
    double *x = malloc(dset->n * sizeof *x);

    if (x == NULL) {
	*err = E_ALLOC;
    } else {
	int t;

	for (t=0; t<dset->n; t++) {
	    x[t] = NADBL;
	}

	*err = gretl_fill_random_array(x, dset->t1, dset->t2, dist,
				       parm, vecp1, vecp2);
    }

    return x;
}

gretl_matrix *gretl_get_random_matrix (int dist, const double *parm,
				       int rows, int cols,
				       int *err)
{
    gretl_matrix *m = NULL;
    int n = rows * cols;

    if (n <= 0) {
	*err = E_INVARG;
    } else {
	m = gretl_matrix_alloc(rows, cols);
	if (m == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	} else {
	    *err = gretl_fill_random_array(m->val, 0, n - 1, dist,
					   parm, NULL, NULL);
	}
    }

    return m;
}

double gretl_get_random_scalar (int dist, const double *parm,
				int *err)
{
    double x;

    *err = gretl_fill_random_array(&x, 0, 0, dist,
				   parm, NULL, NULL);

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
 * @dist: distribution code.
 * @parm: array holding 1 or 2 parameter values.
 * @x: the value in the distribution.
 * @pv: the p-value.
 * @prn: gretl printer.
 *
 * Prints the p-value information in a consistent manner.
 */

void print_pvalue (int dist, const double *parm, double x,
		   double pv, PRN *prn)
{
    double pc;
    int err;

    switch (dist) {

    case D_NORMAL:
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

    case D_STUDENT:
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

    case D_CHISQ:
	pprintf(prn, "%s(%d): ", _("Chi-square"), (int) parm[0]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = chisq_cdf(parm[0], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case D_SNEDECOR:
	pprintf(prn, "F(%d, %d): ", (int) parm[0], (int) parm[1]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = snedecor_cdf((int) parm[0], (int) parm[1], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case D_GAMMA:
	pprintf(prn, _("Gamma (shape %g, scale %g, mean %g, variance %g):"
		       "\n area to the right of %g = %g\n"), 
		parm[0], parm[1], parm[0] * parm[1], 
		parm[0] * parm[1] * parm[1],
		x, pv);
	break;

    case D_BINOMIAL:
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

    case D_POISSON:
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

    case D_WEIBULL:
	pprintf(prn, _("Weibull (shape = %g, scale = %g): "), 
		parm[0], parm[1]);
	err = print_pv_string(x, pv, prn);
	if (err) return;
	pc = weibull_cdf(parm[0], parm[1], x);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case D_GED:
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
 * @dset: dataset struct.
 * @prn: gretl printing struct.
 *
 * Calculates and prints the probability that a random variable 
 * distributed as specified in the command line @str exceeds the 
 * value indicated in @str.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int batch_pvalue (const char *str, DATASET *dset, PRN *prn)
{
    double pv = NADBL;
    char line[MAXLEN];
    char **S;
    int dist;
    int i, n, m;
    int err = 0;
    
    if (!strncmp(str, "pvalue ", 7)) {
	str += 7;
    }

    S = gretl_string_split(str, &n, NULL);
    if (S == NULL) {
	return E_ALLOC;
    }

    dist = dist_code_from_string(S[0]);

    if (dist == D_NONE) {
	err = E_INVARG;
    } else {
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
    }

    strings_array_free(S, n);

    if (!err) {
	pv = generate_scalar(line, dset, &err);
    }

    if (!err) {
	print_pvalue(dist, pvargs, pvargs[2], pv, prn);
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

