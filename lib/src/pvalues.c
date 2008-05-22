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

double gamma_function (double x)
{
    double ret = cephes_gamma(x);

    if (get_cephes_errno()) {
	fprintf(stderr, "gamma_function: arg = %g\n", x);
	ret = NADBL;
    }

    return ret;
}

double log_gamma_function (double x)
{
    double ret = cephes_lgamma(x);

    if (get_cephes_errno()) {
	fprintf(stderr, "log_gamma_function: arg = %g\n", x);
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
    equal to the given cumulative probability y.
*/

static double binomial_cdf_inverse (double y, int n, int k)
{
    double p = NADBL;

    if (y >= 0 && n >= 0 && k >= 0) {
	p = bdtri(k, n, y);
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
 * Returns: the two-sided 95 percent critical value for the sample
 * correlation coefficient, sample size @n, or #NADBL on
 * failure.
 */

double rhocrit95 (int n)
{
    double x = NADBL;

    if (n - 2 > 0) {
	x = stdtri(n - 2, 0.975);
	if (get_cephes_errno()) {
	    x = NADBL;
	} else {
	    double x2 = x * x;

	    x = sqrt(x2 / (x2 - 2 + n));
	}
    }
    
    return x;
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

double student_cdf (int df, double x)
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

static double student_cdf_comp (int df, double x)
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
 * student_pvalue_2:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the probability that t(@df) is greater than @x
 * (two-sided, using the absolute value of @x), or
 * #NADBL on failure.
 */

double student_pvalue_2 (int df, double x)
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

#define df_ok(d) (floor(d) == d && d < (double) INT_MAX)

/**
 * student_critval:
 * @df: degrees of freedom.
 * @a: right-tail probability.
 *
 * Returns: the argument x such that the integral from x to 
 * infinity of the t(@df) density is equal to the given
 * probability @a, or #NADBL on failure.
 */

static double student_critval (double df, double a)
{
    double x;

    if (df < 1) {
	return NADBL;
    }    

    if (df_ok(df)) {
	if (a > .10) {
	    x = stdtri((int) df, 1 - a);
	} else {
	    x = -stdtri((int) df, a);
	}
    } else {
	if (a > .10) {
	    x = ndtri(1 - a);
	} else {
	    x = -ndtri(a);
	}
    }

    if (get_cephes_errno()) {
	x = NADBL;
    } 

    return x;
}

double student_cdf_inverse (double df, double a)
{
    double x;

    if (df < 1) {
	return NADBL;
    }

    if (df_ok(df)) {
	x = stdtri((int) df, a);
    } else {
	x = ndtri(a);
    }    

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
 * Returns: the integral from 0 to @x of the chi-squared
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
 * Returns: the integral from @x to infinity of the chi-squared
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
	double x1 = log_gamma_function(p + q);
	double x2 = log_gamma_function(p);
	double x3 = log_gamma_function(q);

	if (!na(x1) && !na(x2) && !na(x3)) {
	    f = exp(x1 - x2 - x3);
	    if (errno) {
		f = NADBL;
	    }
	}
    }

    return f;
}    

static double snedecor_pdf (int m, int n, double x)
{
    double fx = NADBL;

    errno = 0;

    if (m > 0 && n > 0 && x > 0) {
	double xm = m, xn = n;
	double x1 = Binv(0.5*xm, 0.5*xn);
	double x2 = pow(xm / xn, 0.5*xm);
	double x3 = pow(x, 0.5*xm - 1.0);
	double x4 = pow(1.0 + xm/xn * x, 0.5 * (xm+xn));

	if (!errno && !na(x1)) {
	    fx = x1 * x2 * x3 / x4;
	}
    }

    return fx;
}

static double chisq_pdf (int m, double x)
{
    double fx = NADBL;

    errno = 0;

    if (m > 0 && x >= 0) {
	double m2 = m / 2.0;
	double x1 = pow(.5, m2);
	double x2 = gamma_function(m2);
	double x3 = pow(x, m2 - 1.0);
	double x4 = exp(-x / 2.0);

	if (!errno && !na(x2)) {
	    fx = (x1/x2) * x3 * x4;
	}
    }

    return fx;
}

static double student_pdf (double m, double x)
{
    double fx = NADBL;

    errno = 0;

    if (m > 0) {
	double x1 = Binv(0.5 * m, 0.5) / sqrt(m);
	double x2 = m / (m + x * x);
	double x3 = 0.5 * (m + 1.0);

	if (!errno && !na(x1)) {
	    fx = x1 * pow(x2, x3);
	    if (errno) {
		fx = NADBL;
	    }
	}
    }

    return fx;
}

static double weibull_pdf (double k, double l, double x)
{
    double fx = NADBL;

    errno = 0;

    if (k > 0 && l > 0 && x >= 0) {
	double x1 = k / l;
	double x2 = pow(x / l, k - 1.0);
	double x3 = pow(x / l, k);
	double x4 = exp(-x3);

	if (!errno) {
	    fx = x1 * x2 * x4;
	}
    }

    return fx;
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
    }

    return y;
}

/**
 * normal_cdf_inverse:
 * @x: double-precision value.
 * 
 * Returns the argument, y, for which the area under the
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

/**
 * normal_pdf:
 * @x: double-precision value.
 * 
 * Returns: the value of the standard normal PDF evaluated
 * at @x.
 */

double normal_pdf (double x)
{
    return (1 / sqrt(M_2PI)) * exp(-0.5 * x * x);
}

/**
 * general_normal_pdf:
 * @mu: mean.
 * @sigma: standard deviation.
 * @x: abscissa value.
 * 
 * Returns: the value of the normal PDF with mean @mu
 * and standard deviation @sigma evaluated at @x.
 */

double general_normal_pdf (double mu, double sigma, double x)
{
    if (na(x) || na(mu) || na(sigma) || sigma <= 0.0) {
	return NADBL;
    } else {
	double d = x - mu;

	return (1 / (sigma * sqrt(M_2PI))) * 
	    exp(-0.5 * d * d / (sigma * sigma));
    }
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
    return (x * x) / 2 - 0.91893853320467274178;
}

/**
 * bvnorm_cdf:
 * @a: input value.
 * @b: input value.
 * @rho: input value, correlation coefficient.
 *
 * Ripped and adapted from Gnumeric, with a bug corrected for the case
 * (a * b < 0) && (rho < 0).
 *
 * Returns: for (x,y) a bivariate standard Normal rv with correlation
 * coefficient rho, the joint probability that (x < a) && (y < b), or
 * #NADBL on failure.
 */

double bvnorm_cdf (double a, double b, double rho)
{
    static const double x[] = {0.24840615, 0.39233107, 0.21141819, 
			       0.03324666, 0.00082485334};
    static const double y[] = {0.10024215, 0.48281397, 1.0609498, 
			       1.7797294, 2.6697604};

    double den, a1, b1;
    double ret = NADBL;
    int i, j;

    if (fabs(rho) > 1) {
	return NADBL;
    }	

    if (rho == 0.0) {
	/* joint prob is just the product of the marginals */
	return normal_cdf(a) * normal_cdf(b);
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
	ret = normal_cdf(a) - bvnorm_cdf(a, -b, -rho);
    } else if (a >= 0 && b <= 0 && rho > 0) {
	ret = normal_cdf(b) - bvnorm_cdf(-a, b, -rho);
    } else if (a >= 0 && b >= 0 && rho < 0) {
	ret = normal_cdf(a) + normal_cdf(b) - 1 + bvnorm_cdf(-a, -b, rho);
    } else if ((a * b * rho) > 0) {
	int sgna = (a < 0)? -1 : 1;
	int sgnb = (b < 0)? -1 : 1;
	double rho1, rho2, tmp, delta;

	tmp = sqrt((a * a) - 2 * rho * a * b + (b * b));
	rho1 = (rho * a - b) * sgna / tmp;
	rho2 = (rho * b - a) * sgnb / tmp;
	delta = (sgna * sgnb && (rho > 0))? 0 : 0.5;

	ret = (bvnorm_cdf(a, 0, rho1) + bvnorm_cdf(b, 0, rho2) - delta);
    }    

    return ret;
}

/* Control 1 : s1, s2 = shape, scale
           2 : s1, s2 = mean, variance
   Returns NADBL on error.
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

    /* cephes expects inverse-scale, shape, I think >8-} */

    p = gdtr(1.0 / scale, shape, x);
    if (get_cephes_errno()) {
	p = NADBL;
    }

    return p;
}

/* Control 1 : s1, s2 = shape, scale
           2 : s1, s2 = mean, variance
   Returns NADBL on error.
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

double gamma_pdf (double shape, double scale, double x)
{
    double p = NADBL;

    errno = 0;

    if (shape > 0 && scale > 0 && x > 0) {
	double x1 = pow(x, shape - 1.0);
	double x2 = exp(-x/scale);
	double x3 = pow(scale, shape);
	double x4 = gamma_function(shape);

	if (!errno && !na(x4)) {
	    p = x1 * x2 / (x3 * x4);
	}
    }

    return p;
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

static double poisson_pmf (double lambda, int k)
{
    double den, l0, p;

    if (lambda <= 0 || k < 0) {
	return NADBL;
    }

    den = x_factorial((double) k);
    l0 = exp(-lambda);

    if (na(den) || isinf(den) || isnan(den)) {
	p = NADBL;
    } else {
	p = l0 * pow(lambda, (double) k) / den;
    }

    if (na(p) || isinf(p) || isnan(p)) {
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

static double weibull_cdf (double shape, double scale, double x)
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

/* order in x: [params], alpha, critval */

void print_critval (char st, double *x, double c, PRN *prn)
{
    int apos = 0;
    double a;

    switch (st) {
    case 'z':
	pprintf(prn, "%s", _("Standard normal distribution"));
	break;
    case 't':
	pprintf(prn, "t(%g)", x[0]);
	apos = 1;
	break;
    case 'X':
	pprintf(prn, _("Chi-square(%g)"), x[0]);
	apos = 1;
	break;
    case 'F':
	pprintf(prn, "F(%g, %g)", x[0], x[1]);
	apos = 2;
	break;
    case 'B':
	pprintf(prn, "Binomial (P = %g, %g trials)", x[0], x[1]);
	apos = 2;
	break;
    case 'P':
	pprintf(prn, "Poisson (mean = %g)", x[0]);
	apos = 1;
	break;
    case 'W':
	pprintf(prn, "Weibull (shape = %g, scale = %g)", x[0], x[1]);
	apos = 2;
	break;
    }

    a = x[apos];

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

static double dparm[3];

static void dparm_set (const double *p)
{
    int i;

    for (i=0; i<3; i++) {
	dparm[i] = p[i];
    }
}

double gretl_get_cdf_inverse (char st, double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = normal_cdf_inverse(p[0]);
    } else if (st == 't') {
	x = student_cdf_inverse(p[0], p[1]);
    } else if (st == 'X') {
	x = chisq_cdf_inverse((int) p[0], p[1]);
    } else if (st == 'F') {
	x = snedecor_cdf_inverse((int) p[0], (int) p[1], p[2]);
    } else if (st == 'B') {
	x = binomial_cdf_inverse((int) p[0], (int) p[1], p[2]);
    } else if (st == 'P') {
	x = poisson_cdf_inverse((int) p[0], p[1]);
    } 

    return x;
}

double gretl_get_critval (char st, double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = normal_critval(p[0]);
    } else if (st == 't') {
	x = student_critval(p[0], p[1]);
    } else if (st == 'X') {	
	x = chisq_critval((int) p[0], p[1]);
    } else if (st == 'F') {
	x = snedecor_critval((int) p[0], (int) p[1], p[2]);
    } else if (st == 'B') {
	x = binomial_critval(p[0], (int) p[1], p[2]);
    } else if (st == 'P') {
	x = poisson_critval(p[0], p[1]);
    } else if (st == 'W') {
	x = weibull_critval(p[0], p[1], p[2]);
    }

    return x;
}

double gretl_get_cdf (char st, double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = normal_cdf(p[0]);
    } else if (st == 't') {
	x = student_cdf((int) p[0], p[1]);
    } else if (st == 'X') {
	x = chisq_cdf((int) p[0], p[1]);
    } else if (st == 'F') {
	x = snedecor_cdf((int) p[0], (int) p[1], p[2]);
    } else if (st == 'G') {
	x = gamma_cdf(p[0], p[1], p[2], 1);
    } else if (st == 'B') {
	x = binomial_cdf(p[0], (int) p[1], (int) p[2]);
    } else if (st == 'D') {
	x = bvnorm_cdf(p[1], p[2], p[0]);
    } else if (st == 'P') {
	x = poisson_cdf(p[0], (int) p[1]);
    } else if (st == 'W') {
	x = weibull_cdf(p[0], p[1], p[2]);
    }

    return x;
}

double gretl_get_pdf (char st, double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = normal_pdf(p[0]);
    } else if (st == 't') {
	x = student_pdf(p[0], p[1]);
    } else if (st == 'X') {
	x = chisq_pdf((int) p[0], p[1]);
    } else if (st == 'F') {
	x = snedecor_pdf((int) p[0], (int) p[1], p[2]);
    } else if (st == 'G') {
	x = gamma_pdf(p[0], p[1], p[2]);
    } else if (st == 'W') {
	x = weibull_pdf(p[0], p[1], p[2]);
    }

    return x;
}

double gretl_get_pvalue (char st, const double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = normal_cdf_comp(p[0]);
    } else if (st == 't') {
	x = student_cdf_comp((int) p[0], p[1]);
    } else if (st == 'X') {
	x = chisq_cdf_comp((int) p[0], p[1]);
    } else if (st == 'F') {
	x = snedecor_cdf_comp((int) p[0], (int) p[1], p[2]);
    } else if (st == 'G') {
	x = gamma_cdf_comp(p[0], p[1], p[2], 1);
    } else if (st == 'B') {
	x = binomial_cdf_comp(p[0], (int) p[1], (int) p[2]);
    } else if (st == 'P') {
	x = poisson_cdf_comp(p[0], (int) p[1]);
    } else if (st == 'W') {
	x = weibull_cdf_comp(p[0], p[1], p[2]);
    }

    dparm_set(p);

    return x;
}

double *gretl_get_random_series (char st, const double *p,
				 const DATAINFO *pdinfo,
				 int *err)
{
    double *x = malloc(pdinfo->n * sizeof *x);
    int t;

    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (t=0; t<pdinfo->n; t++) {
	x[t] = NADBL;
    }

    if (st == 'u') {
	double min = p[0], max = p[1];

	*err = gretl_rand_uniform_minmax(x, pdinfo->t1, pdinfo->t2, min, max);
    } else if (st == 'z') {
	double mu = p[0], sd = p[1];

	*err = gretl_rand_normal_full(x, pdinfo->t1, pdinfo->t2, mu, sd);
    } else if (st == 't') {
	int v = p[0];

	*err = gretl_rand_student(x, pdinfo->t1, pdinfo->t2, v);
    } else if (st == 'X') {
	int v = p[0];

	*err = gretl_rand_chisq(x, pdinfo->t1, pdinfo->t2, v);
    } else if (st == 'F') {
	int v1 = p[0], v2 = p[1];

	*err = gretl_rand_F(x, pdinfo->t1, pdinfo->t2, v1, v2);
    } else if (st == 'G') {
	double shape = p[0], scale = p[1];

	*err = gretl_rand_gamma(x, pdinfo->t1, pdinfo->t2, shape, scale);
    } else if (st == 'B') {
	double pr = p[0];
	int n = p[1];

	*err = gretl_rand_binomial(x, pdinfo->t1, pdinfo->t2, n, pr);
    } else if (st == 'P') {
	double m = p[0]; /* FIXME? */

	gretl_rand_poisson(x, pdinfo->t1, pdinfo->t2, &m, 0);
    } else if (st == 'W') {
	double shape = p[0], scale = p[1];

	*err = gretl_rand_weibull(x, pdinfo->t1, pdinfo->t2, shape, scale);
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

void print_pvalue (char st, double *p, double pv, PRN *prn)
{
    double pc;
    int err;

    switch (st) {

    case 'z':
    case 'n':
    case 'N':
    case '1':
	pprintf(prn, "\n%s: ", _("Standard normal"));
	err = print_pv_string(p[0], pv, prn);
	if (err) return;
	if (pv < 0.5) {
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pv, 1 - 2 * pv);
	} else {
	    pc = normal_cdf(p[0]);
	    pprintf(prn, _("(to the left: %g)\n"), pc);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pc, 1 - 2 * pc);
	}
	break;

    case 't':
    case '2':
	pprintf(prn, "\nt(%d): ", (int) p[0]);
	err = print_pv_string(p[1], pv, prn);
	if (err) return;
	if (pv < 0.5) {
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pv, 1 - 2 * pv);
	} else {
	    pc = student_cdf((int) p[0], p[1]);
	    pprintf(prn, _("(to the left: %g)\n"), pc);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pc, 1 - 2 * pc);
	}
	break;

    case 'X':
    case 'x':
    case 'c':
    case '3':
	pprintf(prn, "\n%s(%d): ", _("Chi-square"), (int) p[0]);
	err = print_pv_string(p[1], pv, prn);
	if (err) return;
	pc = chisq_cdf(p[0], p[1]);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'F':
    case 'f':
    case '4':
	pprintf(prn, "\nF(%d, %d): ", (int) p[0], (int) p[1]);
	err = print_pv_string(p[2], pv, prn);
	if (err) return;
	pc = snedecor_cdf((int) p[0], (int) p[1], p[2]);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'G':
    case 'g':
    case '5':
	pprintf(prn, _("\nGamma (shape %g, scale %g, mean %g, variance %g):"
		       "\n area to the right of %g = %g\n"), 
		p[0], p[1], p[0] * p[1], p[0] * p[1] * p[1],
		p[2], pv);
	break;

    case 'B':
    case 'b':
    case '6':
	pprintf(prn, _("\nBinomial (p = %g, n = %d):"
		       "\n Prob(x > %d) = %g\n"), 
		p[0], (int) p[1], (int) p[2], pv);
	pc = binomial_cdf(p[0], p[1], p[2]);
	if (p[2] > 0) {
	    pprintf(prn, _(" Prob(x <= %d) = %g\n"), (int) p[2], pc);
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) p[2],
		    pc - binomial_cdf(p[0], p[1], p[2] - 1));
	} else {
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) p[2], pc);
	}
	break;

    case 'p':
    case 'P':
    case '8':
	pprintf(prn, _("\nPoisson (mean = %g): "), p[0]);
	err = print_pv_string(p[1], pv, prn);
	if (err) return;
	pc = poisson_cdf(p[0], (int) p[1]);
	if (p[1] > 0) {
	    pprintf(prn, _(" Prob(x <= %d) = %g\n"), (int) p[1], pc);
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) p[1],
		    poisson_pmf(p[0], p[1]));
	} else {
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) p[1], pc);
	}
	break;	

    case 'w':
    case 'W':
    case '9':
	pprintf(prn, _("\nWeibull (shape = %g, scale = %g): "), 
		p[0], p[1]);
	err = print_pv_string(p[2], pv, prn);
	if (err) return;
	pc = weibull_cdf(p[0], p[1], p[2]);
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
    double x = NADBL;
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
	x = generate_scalar(line, pZ, pdinfo, &err);
    }

    if (!err) {
	print_pvalue(st, dparm, x, prn);
    }

    return err;
}

