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
#include "gretl_private.h"
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
 * x_factorial:
 * @x: input value.
 * 
 * Returns: the factorial of int(x), cast to a double, or
 * NADBL on failure.
 *
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
 * Returns: the 95 percent critical value for the t distribution
 * with @df degrees of freedom (two-sided)
 *
 */

double tcrit95 (int df)
{
    return stdtri(df, 0.975);
}

/**
 * rhocrit95:
 * @n: sample size.
 * 
 * Returns: the 95 percent critical value for the sample correlation
 * coefficient, sample size @n.
 *
 */

double rhocrit95 (int n)
{
    double x = stdtri(n - 2, 0.975);
    
    return sqrt(x*x / (x*x - 2 + n));
}

/**
 * gaussprob:
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the probability that z is greater than @x
 * (two-sided, using absolute values).
 */

double gaussprob (double x)
{
    return 2.0 * (1.0 - ndtr(fabs(x)));
}

/**
 * tprob:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the probability that t(@df) is greater than @x
 * (two-sided, using absolute values).
 */

double tprob (double x, int df)
{
    if (df <= 0) {
	return -1.0;
    } else {
#if 0
	return fdist(x * x, 1, df);
#else
	double s = stdtr(df, fabs(x));
	double ret = 2.0 * (1.0 - s);

	if (ret < 0.0) ret = 0.0;
	return ret;
#endif
    }
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
 * normal:
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the probability that a random variable distributed as
 * N(0, 1) is greater than @x.
 *
 */

double normal (double x)
{
    return 1.0 - ndtr(x);
}

double normal_cdf (double x)
{
    return ndtr(x);
}

double normal_pdf (double x)
{
    return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
}

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

/**
 * batch_pvalue:
 * @str: the command line, which should be of one of the following forms:
 * pvalue 1 x (Normal distribution);
 * pvalue 2 df x (t-distribution);
 * pvalue 3 df x (Chi-square);
 * pvalue 4 dfn dfd x (F-distribution); or
 * pvalue 5 mean variance x (Gamma distribution).
 * @Z: the data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 * 
 * Returns: the probability that a random variable distributed as
 * specified in the command line @str exceeds the value indicated
 * in @str, or a negative number in case of failure.
 *
 */

double batch_pvalue (const char *str, 
		     const double **Z, const DATAINFO *pdinfo, 
                     PRN *prn)
{
    int i, df1 = 0, df2 = 0;
    char stat = 0;
    double xx = NADBL, mean = 0, variance = 0, xval = 0, tmp;
    char cmd[7], df1str[9], df2str[9], fstr[9]; 
    int gotvar, err = 0;

    for (;;) {
	if (sscanf(str, "%c,%[^,],%[^,],%s", &stat, df1str, df2str, fstr) == 4)
	    break;
	else *df1str = *df2str = *fstr = '\0';
	if (sscanf(str, "%c,%[^,],%s", &stat, df1str, fstr) == 3)
	    break;
	else *df1str = *df2str = *fstr = '\0';
	if (sscanf(str, "%c,%s", &stat, fstr) == 2)
	    break;
	else *df1str = *df2str = *fstr = '\0';
	if (sscanf(str, "%s %c %s %s %s", cmd, &stat, df1str, df2str, fstr) == 5)
	    break;
	else *df1str = *df2str = *fstr = '\0';
	if (sscanf(str, "%s %c %s %s", cmd, &stat, df1str, fstr) == 4)
	    break;
	else *df1str = *df2str = *fstr = '\0';
	if (sscanf(str, "%s %c %s", cmd, &stat, fstr) == 3)
	    break;
	else *df1str = *df2str = *fstr = '\0';
	break;
    }

    if (isalpha((unsigned char) *df1str)) {
	gotvar = 0;
	for (i=0; i<pdinfo->v; i++) {
	    if (strcmp(df1str, pdinfo->varname[i]) == 0) {
		gotvar = 1;
		df1 = (int) Z[i][0];
		mean = Z[i][0];
		break;
	    }
	}
	if (!gotvar) {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), df1str);
	    err = 1;
	}
    } else {
	if (*df1str && check_atof(df1str)) {
	    err = 1;
	} else {
	    df1 = atoi(df1str);
	    mean = atof(df1str);
	}
    }

    if (isalpha((unsigned char) *df2str)) {
	gotvar = 0;
	for (i=0; i<pdinfo->v; i++) {
	    if (strcmp(df2str, pdinfo->varname[i]) == 0) {
		gotvar = 1;
		df2 = (int) Z[i][0];
		variance = Z[i][0];
		break;
	    }
	}
	if (!gotvar) {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), df2str);
	    err = 1;
	}
    } else {
	if (*df2str && check_atof(df2str)) {
	    err = 1;
	} else {
	    df2 = atoi(df2str);
	    variance = atof(df2str);
	}
    }

    if (isalpha((unsigned char) *fstr)) {
	gotvar = 0;
	for (i=0; i<pdinfo->v; i++) {
	    if (strcmp(fstr, pdinfo->varname[i]) == 0) {
		gotvar = 1;
		xval = get_xvalue(i, Z, pdinfo);
		if (na(xval)) {
		    pputs(prn, _("\nstatistic has missing value code\n"));
		    return NADBL;
		}		
		break;
	    }
	}
	if (!gotvar) {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), fstr);
	    err = 1;
	}
    } else {
	if (*fstr && check_atof(fstr)) {
	    err = 1;
	} else {
	    xval = atof(fstr);
	}
    }

    if (err) {
	print_gretl_errmsg(prn);
	return NADBL;
    }

    switch (stat) {

    case '1':
    case 'z':
    case 'n':
	tmp = xval;
	if (xval > 0.0) tmp = -tmp;
	xx = normal_cdf(tmp);
	if (xx < 0) {
	    pputs(prn, _("\np-value calculation failed\n"));
	    return -1;
	}	
	pprintf(prn, _("\nStandard normal: area to the %s "
		"of %g = %g\n"), (xval > 0)? _("right"): _("left"), 
		xval, xx);
	pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		2.0 * xx, 1.0 - 2.0 * xx);
	return xx;

    case '2':
    case 't':
	if (!*fstr || !*df1str) {
	    pputs(prn, _("\npvalue for t: missing parameter\n"));
	    return -1;
	}
	xx = tprob(xval, df1); /* this is two-tailed */
	if (xx < 0) {
	    pputs(prn, _("\np-value calculation failed\n"));
	    return -1;
	}
	xx *= 0.5;
	pprintf(prn, _("\nt(%d): area to the %s of %g = %g\n"), 
		df1, (xval > 0)? _("right"): _("left"),
		xval, xx);
	pprintf(prn, _("(two-tailed value = %g; complement = %g)\n"), 
		2.0 * xx, 1.0 - 2.0 * xx);
	return xx;

    case '3':
    case 'c':
    case 'x':
    case 'X':
	if (!*fstr || !*df1str) {
	    pputs(prn, _("\npvalue for chi-square: missing parameter\n"));
	    return -1;
	}
	xx = chisq(xval, df1);
	if (xx < 0) {
	    pputs(prn, _("\np-value calculation failed\n"));
	    return -1;
	}
	pprintf(prn, _("\nChi-square(%d): area to the right of %g = %g\n"), 
		df1, xval, xx);
	pprintf(prn, _("(to the left: %g)\n"), 1.0 - xx);
	return xx;

    case '4':
    case 'f':
    case 'F':
	if (!*fstr || !*df1str || !*df2str) {
	    pputs(prn, _("\npvalue for F: missing parameter\n"));
	    return -1;
	}
	xx = fdist(xval, df1, df2);
	if (xx < 0) {
	    pputs(prn, _("\np-value calculation failed\n"));
	    return -1;
	}
	pprintf(prn, _("\nF(%d, %d): area to the right of %g = %g\n"), 
		df1, df2, xval, xx);
	pprintf(prn, _("(to the left: %g)\n"), 1.0 - xx);
	return xx;

    case '5':
    case 'g':
    case 'G':
	xx = gamma_dist(mean, variance, xval, 2);
	if (na(xx))
	    pputs(prn, _("\nError computing gamma distribution\n"));
	else
	    pprintf(prn, _("\nGamma (mean %g, variance %g, shape %g, scale %g):"
		    "\n area to the right of %g = %g\n"), 
		    mean, variance, mean*mean/variance, variance/mean,
		    xval, 1.0 - xx);
	return xx;

    default:
	pputs(prn, _("\nunrecognized pvalue code\n"));
	return NADBL;
    }
}

static void putxx (double xx)
{
    if (xx < 0.0001) puts("< 0.0001");
    else printf("%g\n", xx);
}

/**
 * interact_pvalue:
 * 
 * P-value finder function for interactive use at the command prompt.
 *
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

/* ........................................................ */

static void pnormal (void)
{
    double xx, zx; 

    printf("%s", _(negval));
    zx = getx();
    if (zx < 0.0) return;
    xx = normal(zx);
    printf(_("\nFor the standard normal, area (one-tail) to the "
	   "right of %g is "), zx);
    putxx(xx);
}

/* ........................................................ */

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

/* ........................................................ */

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

/* ........................................................ */

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

/* ........................................................ */

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

/* ........................................................ */

static double getx (void)
{
    double aa;

    if ((fscanf(stdin, "%lf", &aa)) == 1) {
	return aa;
    }
    return -1;
}

/* ........................................................ */

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
 * Returns the critical value for F(df1, df2, a).
 *
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

/* ........................................................ */

static int parse_critical_input (const char *str, int *i, 
				 int *df, int *n)
{
    *i = -1;

    if (sscanf(str, "critical F %d %d", df, n) == 2) *i = 3;
    else if (sscanf(str, "critical X %d", df) == 1) *i = 2;
    else if (sscanf(str, "critical t %d", df) == 1) *i = 1;
    else if (sscanf(str, "critical d %d", n) == 1) *i = 4;
    else if (!strncmp(str, "critical N", 10)) *i = 0;

    return (*i == -1);
}

/**
 * print_critical:
 * @line: the command line, which should be of one of the following forms:
 * critical t df (student's t)
 * critical X df (chi-square)
 * critical F dfn dfd (F distribution)
 * @prn: gretl printing struct.
 *
 * Prints critical values for the specified distribution at the
 * commonly used significance levels.
 *
 * Returns: 0 if successful, 1 on error.
 *
 */

int print_critical (const char *line, PRN *prn)
{
    void *handle = NULL;
    void *funp = NULL;
    void (*norm_table)(PRN *, int) = NULL;
    void (*dw)(int, PRN *) = NULL;
    void (*tcrit)(int, PRN *, int) = NULL;
    void (*chicrit)(int, PRN *, int) = NULL;
    int i, n = -1, df = -1, err = 0;

    if (parse_critical_input(line, &i, &df, &n)) {
	pputs(prn, _("Invalid input\n"));
	err = 1;
    }

    if ((0 < i && i < 4 && df <= 0) || (i == 3 && n <= 0)) {
	pputs(prn, _("Invalid degrees of freedom\n"));
	err = 1;
    }
    else if (i == 4 && n <= 0) {
	pputs(prn, _("Invalid sample size\n"));
	err = 1;
    }    

    if (err) return 1;

    switch (i) {
    case 0: /* normal */
	funp = norm_table = get_plugin_function("norm_lookup", &handle);
	break;
    case 1: /* t */
	funp = tcrit = get_plugin_function("t_lookup", &handle);
	break;
    case 2: /* chi-square */
	funp = chicrit = get_plugin_function("chisq_lookup", &handle);
	break;
    case 3: /* F */
	break;
    case 4: /* DW */
	funp = dw = get_plugin_function("dw_lookup", &handle);
	break;
    default:
	break;
    }

    if (i != 3 && funp == NULL)  {
	pputs(prn, _("Couldn't load plugin function\n"));
	return 1;
    }
    
    switch (i) {
    case 0:
	(*norm_table)(prn, 0);
	break;
    case 1:
	(*tcrit)(df, prn, 0);
	break;
    case 2:
	(*chicrit)(df, prn, 0);
	break;	
    case 3:
	pprintf(prn, _("Approximate critical values of F(%d, %d)\n\n"),
		df, n);
	pprintf(prn, _(" 10%% in right tail %.2f\n"), f_crit_a(.10, df, n));
	pprintf(prn, "  5%%               %.2f\n", f_crit_a(.05, df, n));	
	pprintf(prn, "  1%%               %.2f\n", f_crit_a(.01, df, n));
	break;
    case 4:
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
				      int *i, int *dfn, int *dfd, 
				      double *a)
{
    char dfnstr[VNAMELEN], dfdstr[VNAMELEN];
    char astr[VNAMELEN];
    double val;
    *i = -1;

    dfnstr[0] = dfdstr[0] = astr[0] = '\0';

    if (sscanf(str, "F,%8[^,],%8[^,],%8s", dfnstr, dfdstr, astr) == 3) {
	*i = 3;
    } else if (sscanf(str, "X,%8[^,],%8s", dfnstr, astr) == 2) {
	*i = 2;
    } else if (sscanf(str, "t,%8[^,],%8s", dfnstr, astr) == 2) {
	*i = 1;
    } else if (sscanf(str, "N,%8s", astr) == 1) {
	*i = 0;
	*dfn = 500;
    }

    if (*i == -1) return 1;

    if (dfnstr[0] != '\0') {
	val = get_number_or_val(dfnstr, Z, pdinfo);
	if (na(val)) {
	    *i = -1;
	} else {
	    *dfn = val;
	}
    }
    if (dfdstr[0] != '\0') {
	val = get_number_or_val(dfdstr, Z, pdinfo);
	if (na(val)) {
	    *i = -1;
	} else {
	    *dfd = val;
	}
    }
    if (astr[0] != '\0') {
	*a = get_number_or_val(astr, Z, pdinfo);
	if (na(*a) || *a < 0.0) {
	    *i = -1;
	}
    }

    return (*i == -1);
}

double genr_get_critical (const char *line, const double **Z, 
			  const DATAINFO *pdinfo)
{
    double alpha, ret = NADBL;
    int st, dfn = -1, dfd = -1;

    if (parse_genr_critical_input(line, Z, pdinfo,
				  &st, &dfn, &dfd, &alpha)) {
	return NADBL;
    }

    if ((0 < st && st < 4 && dfn <= 0) || (st == 3 && dfd <= 0)) {
	strcpy(gretl_errmsg, _("Invalid degrees of freedom\n"));
	return NADBL;
    }

    if (st == 3) {
	ret = f_crit_a(alpha, dfn, dfd);
    } else if (st == 2) {
	ret = chi_crit_a(alpha, dfn);
    } else {
	ret = sqrt(f_crit_a(2.0 * alpha, 1, dfn));
    } 

    if (ret < 0) {
	ret = NADBL;
    }

    return ret;
}

/* ........................................................ */

/* Functions relating to the gamma distribution.
   Allin Cottrell (cottrell@wfu.edu), October 2000.
   Draws upon the pascal code specialf.pas by Bent Nielsen.
*/

static const double tolerance = 1e-7;

/* internal functions */

static double gamma_integral (double lambda, double x);
static double gamma_integral_expansion (double lambda, double x);
static double gamma_integral_fraction (double lambda, double x);
static double gammadist_wilson_hilferty (double shape, double scale, double x);

/* exported functions */

double gamma_dist (double s1, double s2, double x, int control)
     /* Control 1 : s1, s2 = shape, scale
                2 : s1, s2 = expectation, variance
        Returns NADBL (-999.0) on error 
     */
{
    double shape = 0, scale = 0, xx;

    switch (control) {
    case 1: 
	shape = s1; 
	scale = s2; 
	break;
    case 2: 
	scale = s2/s1; 
	shape = s1/scale; 
	break;
    }
    if ((shape > 20) && (x/scale < 0.9*shape) && (x > 1))
	xx = gammadist_wilson_hilferty(shape, scale, x);
    else {
	xx = gamma_integral(shape, x/scale);
	if (na(xx)) return xx;
	xx /= cephes_gamma(shape);
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
    double x1, x2, x3, xx;
    int i = 0;

    x1 = 1;
    x3 = 1/lambda;
    do {
	i++;
	x1 *= (-x)/i;
	x2 = x1/(lambda + i);
	x3 += x2;
    } while (fabs(x2) >= tolerance && i <= 100);
    xx = x3 * exp(lambda*log(x));
    if (i == 100) return NADBL;
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
    double a = 1-lambda;
    double b = a+x+1;
    double p1 = 1, p2 = 1+x;
    double q1 = x, q2 = b*x;
    double d, p0, q0, r1, r2, xx;
    int c = 0;

    r2 = p2/q2;
    do {
	p0 = p1; p1 = p2; q0 = q1; q1 = q2; r1 = r2;
	a = a+1; b = b+2; c = c+1; d = a*c;
	p2 = b*p1 - d*p0;
	q2 = b*q1 - d*q0;
	if (fabs(q2) > 0)  r2 = p2/q2;
    } while (!((fabs(r2-r1) < tolerance ) || (fabs(r2-r1) < tolerance * r2)
	       || (c == 100)));
    xx = cephes_gamma(lambda);
    xx -= exp(-x + lambda * log(x)) * r2;
    if (c == 100) return NADBL;
    return xx;
}

static double gamma_integral (double lambda, double x)
{
    double xx;

    if (x < 0)  { 
	return NADBL;
    }
    if (x < tolerance) {
	xx = 0;
    } else if ((x <= 1) || (x < 0.9*lambda)) {
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
    double xx, df, xscaled;

    df = 2 * shape;
    xscaled = x * 2/scale;
    xx = exp(log(xscaled/df)/3) - 1 + (double)(2)/9/df;
    xx *= sqrt(9*df/2);
    return 1.0 - normal(xx);
} 











