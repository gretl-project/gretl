/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/*  pvalues.c - routines relating to computation of pvalues
    of sample statistics
*/  

#include "libgretl.h" 
#include "../../cephes/libprob.h"

#include <errno.h>
 
/**
 * binomial_cdf:
 * @k: maximum number of successes.
 * @n: number of trials.
 * @p: probability of success on each trial.
 *
 * Returns: the probability of @k or less successes on
 * @n trials given binomial probability @p, or
 * #NADBL on failure.
 */

double binomial_cdf (int k, int n, double p)
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
 * @k: maximum number of successes.
 * @n: number of trials.
 * @p: probability of success on each trial.
 *
 * Returns: the probability of @k + 1 or more successes on
 * @n trials given binomial probability @p, or
 * #NADBL on failure.
 */

double binomial_cdf_comp (int k, int n, double p)
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
 * CDF evaluated at abs(@x)).
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
    return 1 - ndtr(x);
}

/**
 * t_cdf:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the integral from minus infinity to @x of
 * the t distribution with @df degrees of freedom, or
 * #NADBL on failure.
 */

double t_cdf (double x, int df)
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
 * t_cdf_comp:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the integral from @x to infinity of
 * the t distribution with @df degrees of freedom, or
 * #NADBL on failure.
 */

double t_cdf_comp (double x, int df)
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
 * t_pvalue_2:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the probability that t(@df) is greater than @x
 * (two-sided, using the absolute value of @x), or
 * #NADBL on failure.
 */

double t_pvalue_2 (double x, int df)
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
 * t_critval:
 * @a: right-tail probability.
 * @df: degrees of freedom.
 *
 * Returns: the argument x such that the integral from x to 
 * infinity of the t(@df) density is equal to the given
 * probability @a, or #NADBL on failure.
 */

double t_critval (double a, int df)
{
    double x = stdtri(df, 1 - a);

    if (get_cephes_errno()) {
	x = NADBL;
    } 

    return x;
}

/**
 * chisq_cdf:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the integral from 0 to @x of the chi-squared
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

double chisq_cdf (double x, int df)
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
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the integral from @x to infinity of the chi-squared
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

double chisq_cdf_comp (double x, int df)
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
 * @a: right-tail probability.
 * @df: degrees of freedom.
 * 
 * Returns: the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double chisq_critval (double a, int df)
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

/**
 * f_cdf:
 * @x: the cutoff point in the distribution.
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * 
 * Returns: the integral of the F distribution with @dfn and
 * @dfd degrees of freedom, from 0 to @x, or #NADBL on failure.
 */

double f_cdf (double x, int dfn, int dfd)
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
 * f_cdf_comp:
 * @x: the cutoff point in the distribution.
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * 
 * Returns: the integral of the F distribution with @dfn and
 * @dfd degrees of freedom, from @x to infinity, or #NADBL 
 * on failure.
 */

double f_cdf_comp (double x, int dfn, int dfd)
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
 * f_critval:
 * @a: right-tail probability.
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * 
 * Returns: the F argument x such that the integral
 * from x to infinity of the F density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double f_critval (double a, int dfn, int dfd)
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
 * normal_critval:
 * @a: right-tail probability.
 *
 * Returns: the argument z such that the integral from z to 
 * infinity of the standard normal density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double normal_critval (double a)
{
    double z = ndtri(1.0 - a);

    if (get_cephes_errno()) {
	z = NADBL;
    } 

    return z;
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

static double poisson_pmf (double lambda, int k)
{
    double den = x_factorial((double) k);
    double l0 = exp(-lambda);
    double p;

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
    double x, l0 = exp(-lambda);
    int i;

    if (k > 133) {
	double p;

	x = p = l0;
	for (i=1; i<=k; i++) {
	    p *= lambda / i;
	    x += p;
	}	    
    } else {
	double den, l1 = 1.0;

	x = 0.0;
	for (i=0; i<=k; i++) {
	    den = x_factorial((double) i);
	    x += l0 * l1 / den;
	    l1 *= lambda;
	}
    }

    return x;
}

#if 0 /* not ready */

static double bincoeff (int n, int r)
{
    double c;
    int i, j;

    if (r > n) {
	return 0;
    } else if (r == n) {
	return 1;
    }

    c = (n - r > r)? n - r : r;
    if (c < n) c += 1;

    for (i=c+1, j=2; i<=n; i+=1, j++) {
	c *= i;
	c /= j;
    }

    return c;
}

static double hypergeom_pmf (int k, int N, int m, int n)
{
    double c1, c2, c3;

    c1 = (n + m - N > 0)? n + m - N : 0;
    c2 = (m < n)? m : n;

    if (c1 > c2) {
	c3 = c1;
	c1 = c2;
	c2 = c3;
    }

    if (k < c1 || k > c2) {
	return 0;
    }

    c1 = bincoeff(m, k);
    c2 = bincoeff(N - m, n - k);
    c3 = bincoeff(N, n);

    return c1 * c2 / c3;
}

#endif

static double dparm[3];

static void dparm_set (const double *p)
{
    int i;

    for (i=0; i<3; i++) {
	dparm[i] = p[i];
    }
}

double gretl_get_critval (char st, double *p)
{
    double x = NADBL;

    if (st == 'z') {
	if (p[0] > 0.5) {
	    x = ndtri(1.0 - p[0]);
	} else {
	    x = -ndtri(p[0]);
	}
    } else if (st == 't') {
	if (p[1] > 0.5) {
	    x = stdtri((int) p[0], 1 - p[1]);
	} else {
	    x = -stdtri((int) p[0], p[1]);
	}
    } else if (st == 'X') {	
	x = chisq_critval(p[1], (int) p[0]);
    } else if (st == 'F') {
	x = f_critval(p[2], (int) p[0], (int) p[1]);
    }

    return x;
}

double gretl_get_cdf (char st, double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = normal_cdf(p[0]);
    } else if (st == 't') {
	x = t_cdf(p[1], (int) p[0]);
    } else if (st == 'X') {
	x = chisq_cdf(p[1], (int) p[0]);
    } else if (st == 'F') {
	x = f_cdf(p[2], (int) p[0], (int) p[1]);
    } else if (st == 'G') {
	x = 1.0 - gamma_cdf_comp(p[0], p[1], p[2], 1);
    } else if (st == 'B') {
	x = binomial_cdf((int) p[2], (int) p[1], p[0]);
    } else if (st == 'D') {
	x = bvnorm_cdf(p[0], p[1], p[2]);
    } else if (st == 'P') {
	x = poisson_cdf(p[0], (int) p[1]);
    }

    return x;
}

double gretl_get_pvalue (char st, const double *p)
{
    double x = NADBL;

    if (st == 'z') {
	x = 1.0 - normal_cdf(p[0]);
    } else if (st == 't') {
	x = t_cdf_comp(p[1], (int) p[0]);
    } else if (st == 'X') {
	x = chisq_cdf_comp(p[1], (int) p[0]);
    } else if (st == 'F') {
	x = f_cdf_comp(p[2], (int) p[0], (int) p[1]);
    } else if (st == 'G') {
	x = gamma_cdf_comp(p[0], p[1], p[2], 1);
    } else if (st == 'B') {
	x = binomial_cdf_comp((int) p[2], (int) p[1], p[0]);
    } else if (st == 'P') {
	x = 1.0 - poisson_cdf(p[0], (int) p[1]);
    }

    if (!na(x)) {
	dparm_set(p);
    }

    return x;
}

static void 
print_pv_string (double x, double p, PRN *prn)
{
    char numstr[32];

    sprintf(numstr, "%g", p);

    if (!strcmp(numstr, "1")) {
	pprintf(prn, _("area to the right of %g =~ %g\n"), x, p);
    } else {
	pprintf(prn, _("area to the right of %g = %g\n"), x, p);
    }
}

void print_pvalue (char st, double *p, double pv, PRN *prn)
{
    double pc;

    switch (st) {

    case 'z':
    case 'n':
    case 'N':
    case '1':
	pprintf(prn, "\n%s: ", _("Standard normal"));
	print_pv_string(p[0], pv, prn);
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
	print_pv_string(p[1], pv, prn);
	if (pv < 0.5) {
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2 * pv, 1 - 2 * pv);
	} else {
	    pc = t_cdf(p[1], (int) p[0]);
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
	print_pv_string(p[1], pv, prn);
	pc = chisq_cdf(p[1], p[0]);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'F':
    case 'f':
    case '4':
	pprintf(prn, "\nF(%d, %d): ", (int) p[0], (int) p[1]);
	print_pv_string(p[2], pv, prn);
	pc = f_cdf(p[2], (int) p[0], (int) p[1]);
	pprintf(prn, _("(to the left: %g)\n"), pc);
	break;

    case 'G':
    case 'g':
    case '5':
	pprintf(prn, _("\nGamma (shape %g, scale %g, mean %g, variance %g):"
		       "\n area to the right of %g = %g\n"), 
		p[0], p[1], p[0] * p[1], p[0] * p[1] * p[1],
		p[2], 1 - pv);
	break;

    case 'B':
    case 'b':
    case '6':
	pprintf(prn, _("\nBinomial (p = %g, n = %d):"
		       "\n Prob(x > %d) = %g\n"), 
		p[0], (int) p[1], (int) p[2], pv);
	pc = binomial_cdf(p[2], p[1], p[0]);
	pprintf(prn, _(" Prob(x <= %d) = %g\n"), (int) p[2], pc);
	if (p[2] > 0) {
	    pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) p[2],
		    pc - binomial_cdf(p[2] - 1, p[1], p[0]));
	}		
	break;

    case 'p':
    case 'P':
    case '8':
	pprintf(prn, _("\nPoisson (mean = %g): "), p[0]);
	print_pv_string(p[1], pv, prn);
	pprintf(prn, _(" Prob(x = %d) = %g\n"), (int) p[1],
		poisson_pmf(p[0], (int) p[1]));
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

/* Functions relating to the gamma distribution.
   Allin Cottrell (cottrell@wfu.edu), October 2000.
   Draws upon the pascal code specialf.pas by Bent Nielsen.
*/

static const double gamma_tol = 1e-7;

/* Gamma distribution function.
   See Johnson, Kotz and Balakrishnan: 
   Continuous Univariate Distributions
   vol 1, 2nd ed, Wiley 1994
*/

static double gammadist_wilson_hilferty (double shape, double scale, double x)
{
    double df = 2 * shape;
    double xscaled = x * 2 / scale;
    double xx;

    xx = exp(log(xscaled/df)/3) - 1 + (double)(2) / 9 / df;
    xx *= sqrt(9 * df / 2);

    return normal_cdf(xx);
} 

/* Expansion of Gamma Integral int_0^x t^lambda-1 exp(-t)
   Abramowitz and Stegun p. 262
   Note that the series is alternating.
*/

static double gamma_integral_expansion (double lambda, double x)
{
    double g, x1 = 1;
    double x2, x3 = 1 / lambda;
    int i = 0;

    do {
	i++;
	x1 *= (-x)/i;
	x2 = x1 / (lambda + i);
	x3 += x2;
    } while (fabs(x2) >= gamma_tol && i <= 100);

    if (i == 100) {
	g = NADBL;
    } else {
	g = x3 * exp(lambda * log(x));
    }

    return g;
}

/* Continued Fraction Expansion for Gamma Integral
   int_0^x t^lambda-1 exp(-t) dx
   Abramowitz and Stegun p. 263
   Implemented in Fortran by
   B. L. Shea (1988): Chi-squared and incomplete gamma integral,
   Applied Statistics, vol 37, pp. 466-473.
   See also Schwartz p. 120.      
*/

static double gamma_integral_fraction (double lambda, double x)
{
    double a = 1 - lambda;
    double b = a + x + 1;
    double p1 = 1, p2 = 1 + x;
    double q1 = x, q2 = b * x;
    double r2 = p2 / q2;
    double d, p0, q0, r1, xx;
    int c = 0;

    do {
	p0 = p1; p1 = p2; q0 = q1; q1 = q2; r1 = r2;
	a = a + 1; b = b + 2; c = c + 1; d = a * c;
	p2 = b * p1 - d * p0;
	q2 = b * q1 - d * q0;
	if (fabs(q2) > 0) {
	    r2 = p2 / q2;
	}
	xx = fabs(r2 - r1);
    } while (!((xx < gamma_tol ) || (xx < gamma_tol * r2) || (c == 100)));

    if (c == 100) {
	xx = NADBL;
    } else {
	xx = cephes_gamma(lambda);
	xx -= exp(-x + lambda * log(x)) * r2;
    }

    return xx;
}

static double gamma_integral (double lambda, double x)
{
    double g;

    if (x < 0)  { 
	g = NADBL;
    } else if (x < gamma_tol) {
	g = 0;
    } else if (x <= 1 || x < 0.9 * lambda) {
	g = gamma_integral_expansion(lambda, x);
    } else {
	g = gamma_integral_fraction(lambda, x);
    }

    return g;
}

/* Control 1 : s1, s2 = shape, scale
           2 : s1, s2 = expectation, variance
   Returns NADBL on error 
*/

double gamma_cdf_comp (double s1, double s2, double x, int control)
{
    double shape, scale, xx;

    if (control == 1) {
	shape = s1; 
	scale = s2; 
    } else {
	scale = s2 / s1; 
	shape = s1 / scale; 
    }	

    if ((shape > 20) && (x / scale < 0.9 * shape) && (x > 1)) {
	xx = gammadist_wilson_hilferty(shape, scale, x);
    } else {
	xx = gamma_integral(shape, x / scale);
	if (!na(xx)) {
	    xx /= cephes_gamma(shape);
	}
    }

    return xx;
}
