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
 
static void putxx (double xx);
static void pnormal (void);
static void ptvalue (void);
static void pchisq (void);
static void pfvalue (void);
static void pgamma (void);
static double getx (void);
static void getdf (const char *str);

const char negval[] = N_("\nEnter x value (value < 0 will exit menu): "); 

/**
 * binomial_cdf:
 * @k: maximum number of successes.
 * @n: number of trials.
 * @p: probability of success on each trial.
 *
 * Returns: the probability of @k or less successes on
 * @n trials given binomial probability @p.
 */

double binomial_cdf (int k, int n, double p)
{
    double x = bdtr(k, n, p);

    return (isnan(x))? NADBL : x;
}

/**
 * binomial_pvalue:
 * @k: number of successes.
 * @n: number of trials.
 * @p: probability of success on each trial.
 *
 * Returns: the probability of @k + 1 or more successes on
 * @n trials given binomial probability @p.
 */

double binomial_pvalue (int k, int n, double p)
{
    double x = bdtrc(k, n, p);

    return (isnan(x))? NADBL : x;
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
 * tcrit95:
 * @df: degrees of freedom.
 * 
 * Returns: the two-sided 95 percent critical value for the t 
 * distribution with @df degrees of freedom.
 */

double tcrit95 (int df)
{
    return stdtri(df, 0.975);
}

/**
 * rhocrit95:
 * @n: sample size.
 * 
 * Returns: the two-sided 95 percent critical value for the sample
 * correlation coefficient, sample size @n.
 */

double rhocrit95 (int n)
{
    double x = stdtri(n - 2, 0.975);
    
    return sqrt(x*x / (x*x - 2 + n));
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
    return 2.0 * (1.0 - ndtr(fabs(x)));
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
 * t_pvalue_2:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the probability that t(@df) is greater than @x
 * (two-sided, using the absolute value of @x).
 */

double t_pvalue_2 (double x, int df)
{
    double ret = -1.0;

    if (df > 0) {
	double s = stdtr(df, fabs(x));
	
	ret = 2.0 * (1.0 - s);
	if (ret < 0.0) {
	    ret = 0.0;
	}
    }

    return ret;
}

/**
 * fdist:
 * @x: the cutoff point in the distribution.
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * 
 * Returns: the probability a random variable distributed as
 * F(@dfn, @dfd) is greater than @x, or -1 if either @dfn or @dfd is
 * negative.
 */

double fdist (double x, int dfn, int dfd)
{
    double ret = 1.0;

    if (dfn < 1 || dfd < 1) {
        ret = -1.0;
    } else if (x > 0.0) {
	if (0 && dfd > 300) {
	    /* FIXME igamc underflow ?? */
	    ret = chisq(x * dfn, dfn);
	} else {
	    ret = fdtrc(dfn, dfd, x);
	}
    }

    return ret;
}

/**
 * chisq:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the probability that a random variable distributed as
 * Chi-squared(@df) is greater than @x.
 */

double chisq (double x, int df)
{
    double ret = 1.0;

    if (df < 0) {
	ret = -1.0;
    } else if (x > 0.0) {
	ret = chdtrc(df, x);
    }

    return ret;
}

/**
 * normal_cdf:
 * @x: double-precision value.
 * 
 * Returns: the value of the standard normal CDF evaluated
 * at @x.
 */

double normal_cdf (double x)
{
    return ndtr(x);
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

	if (v > 0 && v < pdinfo->v && !pdinfo->vector[v]) {
	    return Z[v][0];
	}
    } 

    return NADBL;
}

static int val_from_string (const char *s, const double **Z,
			    const DATAINFO *pdinfo, 
			    double *px, int *pn)
{
    double x = NADBL;
    int n = 0, err = 0;

    if (isalpha((unsigned char) *s)) {
	int v = varindex(pdinfo, s);

	if (v < pdinfo->v) {
	    if (pdinfo->vector[v]) {
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

/**
 * batch_pvalue:
 * @str: the command line, which should be of one of the following forms:
 * pvalue 1 x (Normal distribution);
 * pvalue 2 df x (t-distribution);
 * pvalue 3 df x (Chi-square);
 * pvalue 4 dfn dfd x (F-distribution); or
 * pvalue 5 mean variance x (Gamma distribution).
 * pvalue 6 prob n x (Binomial distribution).
 * @Z: the data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 * 
 * Returns: the probability that a random variable distributed as
 * specified in the command line @str exceeds the value indicated
 * in @str, or #NADBL in case of failure.
 */

double batch_pvalue (const char *str, 
		     const double **Z, const DATAINFO *pdinfo, 
                     PRN *prn)
{
    int n1 = 0, n2 = 0, n3 = 0;
    double x1 = 0, x2 = 0, x3 = 0;
    char st = 0;
    double tmp, pv = NADBL;
    char s1[9] = {0};
    char s2[9] = {0};
    char s3[9] = {0};
    char cmd[7];
    int err = 0;

    for (;;) {
	if (sscanf(str, "%c,%[^,],%[^,],%s", &st, s1, s2, s3) == 4) {
	    break;
	}
	*s1 = *s2 = *s3 = '\0';
	if (sscanf(str, "%c,%[^,],%s", &st, s1, s3) == 3) {
	    break;
	} 
	*s1 = *s2 = *s3 = '\0';
	if (sscanf(str, "%c,%s", &st, s3) == 2) {
	    break;
	} 
	*s1 = *s2 = *s3 = '\0';
	if (sscanf(str, "%s %c %s %s %s", cmd, &st, s1, s2, s3) == 5) {
	    break;
	} 
	*s1 = *s2 = *s3 = '\0';
	if (sscanf(str, "%s %c %s %s", cmd, &st, s1, s3) == 4) {
	    break;
	} 
	*s1 = *s2 = *s3 = '\0';
	if (sscanf(str, "%s %c %s", cmd, &st, s3) == 3) {
	    break;
	} 
	*s1 = *s2 = *s3 = '\0';
	break;
    }

    st = normalize_stat(st);

    if (st == 0) {
	pputs(prn, _("\nunrecognized pvalue code\n"));
	return NADBL;
    }

    err = val_from_string(s1, Z, pdinfo, &x1, &n1);

    if (!err) {
	err = val_from_string(s2, Z, pdinfo, &x2, &n2);
    }

    if (!err) {
	err = val_from_string(s3, Z, pdinfo, &x3, &n3);
    }

    if (err) {
	print_gretl_errmsg(prn);
	return NADBL;
    }

    /* check for missing params */
    if (st == 'z' && !*s3) {
	err = 1;
    } else if ((st == 't' || st == 'X') && (!*s1 || !*s3)) {
	err = 1;
    } else if ((st == 'F' || st == 'G' || st == 'B') &&
	(!*s1 || !*s2 || !*s3)) {
	err = 1;
    }

    if (err) {
	pputs(prn, _("\npvalue: missing parameter\n"));
	return NADBL;
    }	

    switch (st) {

    case 'z':
	tmp = x3;
	if (x3 > 0.0) tmp = -tmp;
	pv = normal_cdf(tmp);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {	
	    pprintf(prn, _("\nStandard normal: area to the %s "
			   "of %g = %g\n"), (x3 > 0)? _("right"): _("left"), 
		    x3, pv);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2.0 * pv, 1.0 - 2.0 * pv);
	}
	break;

    case 't':
	pv = t_pvalue_2(x3, n1);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {
	    pv *= 0.5;
	    pprintf(prn, _("\nt(%d): area to the %s of %g = %g\n"), 
		    n1, (x3 > 0)? _("right"): _("left"),
		    x3, pv);
	    pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		    2.0 * pv, 1.0 - 2.0 * pv);
	}
	break;

    case 'X':
	pv = chisq(x3, n1);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {
	    pprintf(prn, _("\nChi-square(%d): area to the right of %g = %g\n"), 
		    n1, x3, pv);
	    pprintf(prn, _("(to the left: %g)\n"), 1.0 - pv);
	}
	break;

    case 'F':
	pv = fdist(x3, n1, n2);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {
	    pprintf(prn, _("\nF(%d, %d): area to the right of %g = %g\n"), 
		    n1, n2, x3, pv);
	    pprintf(prn, _("(to the left: %g)\n"), 1.0 - pv);
	}
	break;

    case 'G':
	pv = gamma_dist(x1, x2, x3, 2);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {
	    pprintf(prn, _("\nGamma (mean %g, variance %g, shape %g, scale %g):"
			   "\n area to the right of %g = %g\n"), 
		    x1, x2, x1*x1/x2, x2/x1,
		    x3, 1.0 - pv);
	}
	break;

    case 'B':
	pv = binomial_pvalue(n3, n2, x1);
	if (pv < 0) {
	    pv = NADBL;
	} else if (!na(pv)) {
	    pprintf(prn, _("\nBinomial (p = %g, n = %d):"
			   "\n Prob(x > %d) = %g\n"), 
		    x1, n2, n3, pv);
	}
	break;

    default:
	break;
    }

    if (na(pv)) {
	pputs(prn, _("\nError computing pvalue\n"));
    }

    return pv;
}

static void putxx (double xx)
{
    if (xx < 0.0001) {
	puts("< 0.0001");
    } else {
	printf("%g\n", xx);
    }
}

/**
 * interact_pvalue:
 * 
 * P-value finder function for interactive use at the command prompt.
 */

void interact_pvalue (void)
{
    int choice, v;
    char ans[3];

    do {
	printf(_("\n\nChoose one of the following distributions: "
	       "\n\n\t1) Standard normal\t\t2) Student's t\n\t3) "
	       "Chi-square\t\t\t4) F\n"
	       "\t5) Gamma\n\n"
	       "Enter your choice (a number < 0 to exit gretl, 0 to quit "
	       "menu, or\n1, 2, 3, 4, or 5): "));

	fflush(stdout);
	v = fscanf(stdin, "%d", &choice);

	if (v == EOF || v == 0) return;
	if (choice < 0) exit(0);
	printf("%d ", choice);

	switch (choice) {
	case 0:
	    putchar('\n');
	    return;
	case 1:		
	    pnormal();
	    break;
	case 2:		
	    ptvalue();
	    break;
	case 3:		
	    pchisq();
	    break;
	case 4:		
	    pfvalue();
	    break;
	case 5:
	    pgamma();
	    break;
	default:	
	    puts(_("\ninvalid choice"));
	    break;
	}

	printf(_("\nDo you want to continue with more pvalues (y or n)? "));
	fflush(stdout);
	fscanf(stdin, "%s", ans);

    } while (ans[0] == 'Y' || ans[0] == 'y');
}

static void pnormal (void)
{
    double xx, zx; 

    printf("%s", _(negval));
    zx = getx();
    if (zx < 0.0) return;
    xx = normal_pvalue_1(zx);
    printf(_("\nFor the standard normal, area (one-tail) to the "
	   "right of %g is "), zx);
    putxx(xx);
}

static void ptvalue (void)
{
    int n;
    double xx, zx, xsq; 

    getdf(" ");
    n = (int) getx();
    if (n <= 0) return;
    printf("%s", _(negval));
    zx = getx();
    if(zx < 0.0) return;
    xsq = zx * zx;
    xx = fdist(xsq, 1, n)/2.0;
    printf(_("\nFor Student's t(%d), area (one-tail) to the "
	   "right of %g is "), n, zx);
    putxx(xx);
}

static void pchisq (void)
{
    int n;
    double xx, zx; 

    getdf(" ");
    n = (int) getx();
    if (n <= 0) return;
    printf("%s", _(negval));
    zx = getx();
    if(zx < 0.0) return;
    xx = chisq(zx, n);
    printf(_("\nFor Chi-square(%d), area to the right of %g is "), 
	   n, zx);
    putxx(xx);
}

static void pfvalue (void)
{
    int m, n;
    double xx, zx; 

    getdf(_(" for the numerator "));
    m = (int) getx();
    if (m <= 0) return;
    getdf(_(" for the denominator "));
    n = (int) getx();
    if (n <= 0) return;
    printf("%s", _(negval));
    zx = getx();
    if (zx < 0.0) return;
    xx = fdist(zx, m, n);
    printf(_("\nFor F(%d, %d), area to the right of %g is "),
	   m, n, zx);
    putxx(xx);
}

static void pgamma (void) 
{
    double mean, variance;
    double xx, zx; 

    printf(_("\nEnter the mean: "));
    mean = getx();
    if (mean <= 0) return;
    printf(_("\nEnter the variance: "));
    variance = getx();
    if (variance <= 0) return;
    printf("%s", _(negval));
    zx = getx();
    if (zx < 0.0) return;
    xx = 1.0 - gamma_dist(mean, variance, zx, 2);
    printf(_("\nFor Gamma (mean %g, variance %g), area to the right of %g is "),
	   mean, variance, zx);
    putxx(xx);
}

static double getx (void)
{
    double aa;

    if ((fscanf(stdin, "%lf", &aa)) == 1) {
	return aa;
    }
    return -1;
}

static void getdf (const char *str)
{
    printf(_("\nEnter d.f.%s(value <= 0 will exit menu): "), str);
}

/**
 * f_crit_a:
 * @a: significance level.
 * @df1: numerator degrees of freedom.
 * @df2: denominator degrees of freedom.
 *
 * Returns: the one-sided critical value for F(@df1, @df2, @a).
 */

double f_crit_a (double a, int df1, int df2)
{
    if (df1 < 1 || df2 < 1 || a < 0.0) {
	return NADBL;
    } else {
	return fdtri(df1, df2, a);
    }
}

static double chi_crit_a (double a, int df)
{
    if (df < 1 || a < 0.0) {
	return NADBL;
    } else {
	return chdtri(df, a);
    }
}

static int 
parse_critical_input (const char *s, int *df, int *n)
{
    int ret = -1;

    if (sscanf(s, "F %d %d", df, n) == 2) {
	ret ='F';
    } else if (sscanf(s, "X %d", df)) {
	ret = 'X';
    } else if (sscanf(s, "t %d", df)) {
	ret = 't';
    } else if (sscanf(s, "d %d", n)) {
	ret = 'd';
    } else if (*s == 'N' || *s == 'z') {
	ret = 'z';
    } 

    return ret;
}

/**
 * print_critical:
 * @line: the command line, which should be of one of the following forms:
 * critical z (Normal); or
 * critical N (Normal); or
 * critical t df (student's t); or
 * critical X df (chi-square); or
 * critical F dfn dfd (F distribution); or
 * critical d n (Durbin-Watson, sample size n).
 * @prn: gretl printing struct.
 *
 * Prints critical values for the specified distribution at the
 * commonly used significance levels.
 *
 * Returns: 0 if successful, 1 on error.
 */

int print_critical (const char *line, PRN *prn)
{
    void *handle = NULL;
    void *funp = NULL;
    void (*norm_table)(PRN *, int) = NULL;
    void (*dw)(int, PRN *) = NULL;
    void (*tcrit)(int, PRN *, int) = NULL;
    void (*chicrit)(int, PRN *, int) = NULL;
    int st, n = -1, df = -1;
    int err = 0;

    st = parse_critical_input(line + 9, &df, &n);

    if (st < 0) {
	pputs(prn, _("Invalid input\n"));
	err = 1;
    } else if ((st == 't' || st == 'X' || st == 'F') && df <= 0) {
	pputs(prn, _("Invalid degrees of freedom\n"));
    } else if (st == 'F' && n <= 0) {
	pputs(prn, _("Invalid degrees of freedom\n"));
	err = 1;
    } else if (st == 'd' && n <= 0) {
	pputs(prn, _("Invalid sample size\n"));
	err = 1;
    }    

    if (err) return 1;

    switch (st) {
    case 'z': /* normal */
	funp = norm_table = get_plugin_function("norm_lookup", &handle);
	break;
    case 't': /* t */
	funp = tcrit = get_plugin_function("t_lookup", &handle);
	break;
    case 'X': /* chi-square */
	funp = chicrit = get_plugin_function("chisq_lookup", &handle);
	break;
    case 'F': /* F */
	break;
    case 'd': /* DW */
	funp = dw = get_plugin_function("dw_lookup", &handle);
	break;
    default:
	break;
    }

    if (st != 'F' && funp == NULL)  {
	pputs(prn, _("Couldn't load plugin function\n"));
	return 1;
    }
    
    switch (st) {
    case 'z':
	(*norm_table)(prn, 0);
	break;
    case 't':
	(*tcrit)(df, prn, 0);
	break;
    case 'X':
	(*chicrit)(df, prn, 0);
	break;	
    case 'F':
	pprintf(prn, _("Approximate critical values of F(%d, %d)\n\n"),
		df, n);
	pprintf(prn, _(" 10%% in right tail %.2f\n"), f_crit_a(.10, df, n));
	pprintf(prn, "  5%%               %.2f\n", f_crit_a(.05, df, n));	
	pprintf(prn, "  1%%               %.2f\n", f_crit_a(.01, df, n));
	break;
    case 'd':
	(*dw)(n, prn);
	break;
    default:
	break;
    }

    if (handle != NULL) {
	close_plugin(handle);
    }

    return 0;
}

static int parse_genr_critical_input (const char *str, 
				      const double **Z, 
				      const DATAINFO *pdinfo,
				      int *st, int *dfn, int *dfd, 
				      double *a)
{
    char dfnstr[VNAMELEN], dfdstr[VNAMELEN];
    char astr[VNAMELEN];
    double val;
    int err = 0;

    dfnstr[0] = dfdstr[0] = astr[0] = '\0';

    if (sscanf(str, "F,%8[^,],%8[^,],%8s", dfnstr, dfdstr, astr) == 3) {
	*st = 'F';
    } else if (sscanf(str, "X,%8[^,],%8s", dfnstr, astr) == 2) {
	*st = 'X';
    } else if (sscanf(str, "t,%8[^,],%8s", dfnstr, astr) == 2) {
	*st = 't';
    } else if (sscanf(str, "N,%8s", astr) || sscanf(str, "z,%8s", astr)) {
	*st = 'z';
	*dfn = 500;
    } else {
	err = 1;
    }

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
    int st, dfn = -1, dfd = -1;

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
	ret = f_crit_a(alpha, dfn, dfd);
    } else if (st == 'X') {
	ret = chi_crit_a(alpha, dfn);
    } else {
	/* normal or t */
	ret = sqrt(f_crit_a(2.0 * alpha, 1, dfn));
    } 

    if (ret < 0) {
	ret = NADBL;
    }

    return ret;
}

/* Functions relating to the gamma distribution.
   Allin Cottrell (cottrell@wfu.edu), October 2000.
   Draws upon the pascal code specialf.pas by Bent Nielsen.
*/

static const double gamma_tol = 1e-7;

/* internal functions */

static double gamma_integral (double lambda, double x);
static double gamma_integral_expansion (double lambda, double x);
static double gamma_integral_fraction (double lambda, double x);
static double gammadist_wilson_hilferty (double shape, double scale, double x);

/* Control 1 : s1, s2 = shape, scale
           2 : s1, s2 = expectation, variance
   Returns NADBL on error 
*/

double gamma_dist (double s1, double s2, double x, int control)
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

/* end exported functions */

static double gamma_integral_expansion (double lambda, double x)
     /* Expansion of Gamma Integral int_0^x t^lambda-1 exp(-t)
	Abramowitz and Stegun p. 262
	Note that the series is alternating.
     */
{
    double xx, x1 = 1.0;
    double x2, x3 = 1.0 / lambda;
    int i = 0;

    do {
	i++;
	x1 *= (-x)/i;
	x2 = x1 / (lambda + i);
	x3 += x2;
    } while (fabs(x2) >= gamma_tol && i <= 100);

    if (i == 100) {
	xx = NADBL;
    } else {
	xx = x3 * exp(lambda * log(x));
    }

    return xx;
}

static double gamma_integral_fraction (double lambda, double x)
     /* Continued Fraction Expansion for Gamma Integral
	int_0^x t^lambda-1 exp(-t) dx
	Abramowitz and Stegun p. 263
	Implemented in Fortran by
	B. L. Shea (1988): Chi-squared and incomplete gamma integral,
	Applied Statistics, vol 37, pp. 466-473.
	See also Schwartz p. 120.      
     */
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
    double xx;

    if (x < 0.0)  { 
	xx = NADBL;
    } else if (x < gamma_tol) {
	xx = 0.0;
    } else if (x <= 1.0 || x < 0.9 * lambda) {
	xx = gamma_integral_expansion(lambda, x);
    } else {
	xx = gamma_integral_fraction(lambda, x);
    }

    return xx;
}

static double gammadist_wilson_hilferty (double shape, double scale, double x)
     /* Gamma distribution function.
	See Johnson, Kotz and Balakrishnan: 
	Continuous Univariate Distributions
        vol 1, 2nd ed, Wiley 1994
     */
{
    double df = 2.0 * shape;
    double xscaled = x * 2.0 / scale;
    double xx;

    xx = exp(log(xscaled/df)/3) - 1 + (double)(2) / 9 / df;
    xx *= sqrt(9 * df / 2);

    return normal_cdf(xx);
} 











