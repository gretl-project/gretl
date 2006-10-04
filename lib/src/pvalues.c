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

    if (p >= 0.0 && n >= 0 && k >= 0) {
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

    if (p >= 0.0 && n >= 0 && k >= 0) {
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

    if (x < 0.0) {
	fact = NADBL;
    } else if (x > 12.0) {
	fact = cephes_gamma(1.0 + x);
	if (get_cephes_errno()) {
	    fact = NADBL;
	}
    } else if (n == 0) {
	fact = 1.0;
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

    if (x < 0.0) {
	lfact = NADBL;
    } else if (x > 12.0) {
	lfact = cephes_lgamma(1.0 + x);
	if (get_cephes_errno()) {
	    lfact = NADBL;
	}
    } else if (n == 0) {
	lfact = 0.0;
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
    double p = (x < 0.0)? ndtr(x) : ndtr(-x);

    return 2.0 * p;
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
    return 1.0 - ndtr(x);
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
	    p = 1.0 - p;
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
	    p *= 2.0;
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
    double x = stdtri(df, 1.0 - a);

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

    if (df > 0 && x >= 0.0) {
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

    if (df > 0 && x >= 0.0) {
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

    if (df > 0 && a >= 0.0) {
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

    if (dfn > 0 && dfd > 0 && x >= 0.0) {
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

    if (dfn > 0 && dfd > 0 && x >= 0.0) {
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

    if (dfn > 0 && dfd > 0 && a >= 0.0) {
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
    return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
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
    return (x * x) / 2.0 - 0.91893853320467274178;
}

static double get_number_or_val (const char *s, 
				 const double **Z,
				 const DATAINFO *pdinfo)
{
    if (numeric_string(s)) {
	return dot_atof(s);
    } else {
	int v = varindex(pdinfo, s);

	if (v > 0 && v < pdinfo->v && var_is_scalar(pdinfo, v)) {
	    return Z[v][0];
	}
    } 

    return NADBL;
}

static int val_from_string (char *s, const double **Z,
			    const DATAINFO *pdinfo, 
			    double *px, int *pn)
{
    double x = NADBL;
    int n = 0, err = 0;

    if (isalpha((unsigned char) *s)) {
	int v = varindex(pdinfo, s);

	if (v < pdinfo->v) {
	    if (var_is_series(pdinfo, v)) {
		x = Z[v][pdinfo->t1];
	    } else {
		x = Z[v][0];
	    }
	    if (na(x)) {
		strcpy(gretl_errmsg, _("Missing values encountered"));
		err = 1;
	    } else {
		n = (int) x;
	    }
	} else {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), s);
	    err = 1;
	}
    } else if (*s != '\0') {
	if (check_atof(s)) {
	    err = 1;
	} else {
	    n = atoi(s);
	    x = atof(s);
	}
    }	

    *px = x;
    *pn = n;

    return err;
}

static char normalize_stat (char c)
{
    switch (c) {
    case '1':
    case 'z':
    case 'n':
    case 'N':
	return 'z'; /* Normal */
    case '2':
    case 't':
	return 't'; /* Student's t */
    case '3':
    case 'c':
    case 'x':
    case 'X':
	return 'X'; /* Chi-square */
    case '4':
    case 'f':
    case 'F':
	return 'F'; /* F */
    case '5':
    case 'g':
    case 'G':
	return 'G'; /* Gamma */
    case '6':
    case 'b':
    case 'B':
	return 'B'; /* Binomial */
    default:
	break;
    }

    return 0;
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
	    x = stdtri((int) p[0], 1.0 - p[1]);
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
	x = 1.0 - gamma_cdf_comp(p[0], p[1], p[2], 2);
    } else if (st == 'B') {
	x = binomial_cdf(p[2], p[1], p[0]);
    }

    return x;
}

static double find_cdf (char st, int n[3], double x[3])
{
    double p = NADBL;

    if (st == 'z') {
	p = normal_cdf(x[0]);
    } else if (st == 't') {
	p = t_cdf(x[1], n[0]);
    } else if (st == 'X') {
	p = chisq_cdf(x[1], n[0]);
    } else if (st == 'F') {
	p = f_cdf(x[2], n[0], n[1]);
    } else if (st == 'G') {
	p = 1.0 - gamma_cdf_comp(x[0], x[1], x[2], 2);
    } else if (st == 'B') {
	p = binomial_cdf(n[2], n[1], x[0]);
    }

    return p;
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
	x = gamma_cdf_comp(p[0], p[1], p[2], 2);
    } else if (st == 'B') {
	x = binomial_cdf_comp(p[2], (int) p[1], (int) p[0]);
    }

    return x;
}

static double find_pvalue (char st, int n[3], double x[3])
{
    double p = NADBL;

    if (st == 'z') {
	p = 1.0 - normal_cdf(x[0]);
    } else if (st == 't') {
	p = t_cdf_comp(x[1], n[0]);
    } else if (st == 'X') {
	p = chisq_cdf_comp(x[1], n[0]);
    } else if (st == 'F') {
	p = f_cdf_comp(x[2], n[0], n[1]);
    } else if (st == 'G') {
	p = gamma_cdf_comp(x[0], x[1], x[2], 2);
    } else if (st == 'B') {
	p = binomial_cdf_comp(n[2], n[1], x[0]);
    }

    return p;
}

static void 
print_pv_string (double x, double p, int left, PRN *prn)
{
    if (left) {
	if (p == 1.0) {
	    pprintf(prn, _("area to the left of %g =~ %g\n"), x, p);
	} else {
	    pprintf(prn, _("area to the left of %g = %g\n"), x, p);
	}
    } else {
	if (p == 1.0) {
	    pprintf(prn, _("area to the right of %g =~ %g\n"), x, p);
	} else {
	    pprintf(prn, _("area to the right of %g = %g\n"), x, p);
	}
    }
}

static double 
find_and_print_pvalue (char st, int n[3], double x[3], PRN *prn)
{
    double pv = NADBL;

    switch (st) {

    case 'z':
	if (x[0] < 0.0) {
	    pv = normal_cdf(x[0]);
	} else {
	    pv = normal_cdf(-x[0]);
	}
	if (!na(pv)) {	
	    pprintf(prn, "\n%s: ", _("Standard normal"));
	    print_pv_string(x[0], pv, (x[0] < 0), prn);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2.0 * pv, 1.0 - 2.0 * pv);
	}
	break;

    case 't':
	if (x[1] < 0.0) {
	    pv = t_cdf(x[1], n[0]);
	} else {
	    pv = t_cdf_comp(x[1], n[0]);
	}
	if (!na(pv)) {
	    pprintf(prn, "\nt(%d): ", n[0]);
	    print_pv_string(x[1], pv, (x[1] < 0), prn);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2.0 * pv, 1.0 - 2.0 * pv);
	}
	break;

    case 'X':
	pv = chisq_cdf_comp(x[1], n[0]);
	if (!na(pv)) {
	    pprintf(prn, "\n%s(%d): ", _("Chi-square"), n[0]);
	    print_pv_string(x[1], pv, 0, prn);
	    pprintf(prn, _("(to the left: %g)\n"), chisq_cdf(x[1], n[0]));
	}
	break;

    case 'F':
	pv = f_cdf_comp(x[2], n[0], n[1]);
	if (!na(pv)) {
	    pprintf(prn, "\nF(%d, %d): ", n[0], n[1]);
	    print_pv_string(x[2], pv, 0, prn);
	    pprintf(prn, _("(to the left: %g)\n"), f_cdf(x[2], n[0], n[1]));
	}
	break;

    case 'G':
	pv = gamma_cdf_comp(x[0], x[1], x[2], 2);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {
	    pprintf(prn, _("\nGamma (mean %g, variance %g, shape %g, scale %g):"
			   "\n area to the right of %g = %g\n"), 
		    x[0], x[1], x[0] * x[0] / x[1], x[1] / x[0],
		    x[2], 1.0 - pv);
	}
	break;

    case 'B':
	pv = binomial_cdf_comp(n[2], n[1], x[0]);
	if (!na(pv)) {
	    double pc = binomial_cdf(n[2], n[1], x[0]);

	    pprintf(prn, _("\nBinomial (p = %g, n = %d):"
			   "\n Prob(x > %d) = %g\n"), 
		    x[0], n[1], n[2], pv);
	    pprintf(prn, _(" Prob(x <= %d) = %g\n"), n[2], pc);
	    if (n[2] > 0) {
		pprintf(prn, _(" Prob(x = %d) = %g\n"), n[2],
			pc - binomial_cdf(n[2] - 1, n[1], x[0]));
	    }		
	}
	break;

    default:
	break;
    }

    return pv;
}

#if 0
double new_batch_pvalue (const char *str, 
			 double **pZ, DATAINFO *pdinfo, 
			 PRN *prn, int *err)
{
    double x = NADBL;
    char line[MAXLEN];
    char **S;
    int n, m;
    
    if (!strncmp(str, "pvalue ", 7)) {
	str += 7;
    }

    while (*str == ' ') str++;

    S = gretl_string_split(s, &n);
    if (S == NULL) {
	*err = E_ALLOC;
	return x;
    }

    strcpy(line, "pvalue(");
    m = 8;
    for (i=0; i<n && !*err; i++) {
	m += strlen(S[i]) + 1;
	if (m > MAXLEN) {
	    *err = E_DATA;
	} else {
	    strcat(line, S[i]);
	    strcat(line, (i == n - 1)? ")" : ",");
	}
    }

    if (*err) {
	x = generate_scalar(line, pZ, pdinfo, err);
    }

    free_strings_array(S, n);

    return x;
}
#endif

/**
 * batch_pvalue:
 * @str: the command line, which should be of one of the following forms:
 * pvalue z x (Normal distribution);
 * pvalue t df x (t-distribution);
 * pvalue X df x (Chi-square);
 * pvalue F dfn dfd x (F-distribution); or
 * pvalue G mean variance x (Gamma distribution).
 * pvalue B prob n x (Binomial distribution).
 * @pZ: pointer to the data array.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 * @opt: (for internal use) %OPT_G forces uses of '.' as
 * decimal separator; %OPT_C actually evaluate cdf instead of
 * p-value.
 * 
 * Returns: the probability that a random variable distributed as
 * specified in the command line @str exceeds the value indicated
 * in @str, or #NADBL in case of failure.
 */

double batch_pvalue (const char *str, 
		     const double **Z, const DATAINFO *pdinfo, 
                     PRN *prn, gretlopt opt)
{
    int n[3] = {0};
    double x[3] = {0.0};
    char st = 0;
    double pv = NADBL;
    char s1[32] = {0};
    char s2[32] = {0};
    char s3[32] = {0};
    int commas = 1;
    int nf = 0;
    int err = 0;

    if (!strncmp(str, "pvalue ", 7)) {
	str += 7;
	commas = 0;
    }

    while (*str == ' ') str++;

#if 0 /* not yet */
    S = gretl_string_split(&n);
    if (S == NULL) {
	*err = E_ALLOC;
    } else {
	strcpy(line, "pvalue(");
	for (i=0; i<n; i++) {
	    strcat(line, S[i]);
	    strcat(line, (i == n - 1)? ")" : ",");
	}
	pv = generate_scalar(line, pZ, pdinfo, err);
	free_strings_array(S, n);
    }
#endif    
    
    if (!sscanf(str, "%c", &st) || 
	(st = normalize_stat(st)) == 0) {
	pputs(prn, _("\nunrecognized pvalue code\n"));
	return NADBL;
    }

    str++;
    while (*str == ' ' || *str == ',') str++;

    if (st == 'z') {
	nf = sscanf(str, "%31s", s1);
	if (nf != 1) err = 1;
    } else if (st == 't' || st == 'X') {
	if (commas) {
	    nf = sscanf(str, "%31[^,],%31s", s1, s2);
	} else {
	    nf = sscanf(str, "%31s %31s", s1, s2);
	}
	if (nf != 2) err = 1;
    } else {
	if (commas) {
	    nf = sscanf(str, "%31[^,],%31[^,],%31s", s1, s2, s3);
	} else {
	    nf = sscanf(str, "%31s %31s %31s", s1, s2, s3);
	}
	if (nf != 3) err = 1;
    }	

    if (err) {
	pputs(prn, _("\npvalue: missing parameter\n"));
	return NADBL;
    }

    if (opt & OPT_G) {
	gretl_push_c_numeric_locale();
    }

    err = val_from_string(s1, Z, pdinfo, &x[0], &n[0]);
    if (!err && nf > 1) {
	err = val_from_string(s2, Z, pdinfo, &x[1], &n[1]);
    } 
    if (!err && nf > 2) {
	err = val_from_string(s3, Z, pdinfo, &x[2], &n[2]);
    }	

    if (opt & OPT_G) {
	gretl_pop_c_numeric_locale();
    }

    if (err) {
	print_gretl_errmsg(prn);
	return NADBL;
    }

    if (opt & OPT_C) {
	pv = find_cdf(st, n, x);
    } else if (opt & OPT_G) {
	pv = find_pvalue(st, n, x);
    } else {
	pv = find_and_print_pvalue(st, n, x, prn);
    }

    if (na(pv)) {
	pputs(prn, _("\nError computing pvalue\n"));
    }

    return pv;
}

static int parse_genr_critical_input (const char *str, 
				      const double **Z, 
				      const DATAINFO *pdinfo,
				      int *st, int *dfn, int *dfd, 
				      double *a)
{
    char dfnstr[VNAMELEN], dfdstr[VNAMELEN];
    char astr[32];
    double val;
    int err = 0;

    dfnstr[0] = dfdstr[0] = astr[0] = '\0';

    gretl_push_c_numeric_locale();

    if (sscanf(str, "F,%8[^,],%8[^,],%24s", dfnstr, dfdstr, astr) == 3) {
	*st = 'F';
    } else if (sscanf(str, "X,%8[^,],%24s", dfnstr, astr) == 2) {
	*st = 'X';
    } else if (sscanf(str, "t,%8[^,],%24s", dfnstr, astr) == 2) {
	*st = 't';
    } else if (sscanf(str, "N,%24s", astr) || sscanf(str, "z,%24s", astr)) {
	*st = 'z';
	*dfn = 500;
    } else {
	err = 1;
    }

    gretl_pop_c_numeric_locale();

    if (err) return err;

    if (*dfnstr != '\0') {
	val = get_number_or_val(dfnstr, Z, pdinfo);
	if (na(val)) {
	    err = 1;
	} else {
	    *dfn = val;
	}
    }

    if (*dfdstr != '\0') {
	val = get_number_or_val(dfdstr, Z, pdinfo);
	if (na(val)) {
	    err = 1;
	} else {
	    *dfd = val;
	}
    }

    if (*astr != '\0') {
	*a = get_number_or_val(astr, Z, pdinfo);
	if (na(*a) || *a < 0.0) {
	    err = 1;
	}
    }

    return err;
}

double genr_get_critical (const char *line, const double **Z, 
			  const DATAINFO *pdinfo)
{
    double alpha = 0.0, ret = NADBL;
    int st = 0, dfn = -1, dfd = -1;

    if (parse_genr_critical_input(line, Z, pdinfo,
				  &st, &dfn, &dfd, &alpha)) {
	return NADBL;
    }

    if ((st == 't' || st == 'X' || st == 'F') && dfn <= 0) {
	strcpy(gretl_errmsg, _("Invalid degrees of freedom\n"));
	return NADBL;
    } else if (st == 'F' && dfd <= 0) {
	strcpy(gretl_errmsg, _("Invalid degrees of freedom\n"));
	return NADBL;
    }	

    if (st == 'F') {
	ret = f_critval(alpha, dfn, dfd);
    } else if (st == 'X') {
	ret = chisq_critval(alpha, dfn);
    } else if (st == 't') {
	if (alpha > 0.5) {
	    ret = stdtri(dfn, 1.0 - alpha);
	} else {
	    ret = -stdtri(dfn, alpha);
	}
    } else {
	/* normal */
	if (alpha > 0.5) {
	    ret = ndtri(1.0 - alpha);
	} else {
	    ret = -ndtri(alpha);
	}
    } 

    return ret;
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
    double df = 2.0 * shape;
    double xscaled = x * 2.0 / scale;
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
    double g, x1 = 1.0;
    double x2, x3 = 1.0 / lambda;
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

    if (x < 0.0)  { 
	g = NADBL;
    } else if (x < gamma_tol) {
	g = 0.0;
    } else if (x <= 1.0 || x < 0.9 * lambda) {
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














