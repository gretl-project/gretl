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
 
extern void _putxx (const double xx);

static void _pnormal (void);
static void _ptvalue (void);
static void _pchisq (void);
static void _pfvalue (void);
static void _pgamma (void);
static double _getvalue (void);
static void _enterdf (const char *str);

double _gammadist (double s1, double s2, double x, int control);

const char negval[] = "\nEnter x value (value < 0 will exit menu): "; 

/**
 * tcrit95:
 * @df: degrees of freedom.
 * 
 * Returns: the 95 percent critical value for the t distribution
 * with @df degrees of freedom.
 *
 */

double tcrit95 (const int df)
{
    double x = 1.960;
    
    while (tprob(x, df) > .05) x += .001;
    return x;
}

/**
 * rhocrit95:
 * @n: sample size.
 * 
 * Returns: the 95 percent critical value for the sample correlation
 * coefficient, sample size @n.
 *
 */

double rhocrit95 (const int n)
{
    double x = 1.960;
    
    while (tprob(x, n - 2) > .05) x += .001;
    return sqrt(x*x / (x*x - 2 + n));
}

/**
 * tprob:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the probability that t(@df) is greater than @x.
 *
 */

double tprob (const double x, const int df)
{
    double xx;

    if (df <= 0) return -1.0;
    xx = x*x;
    return fdist(xx, 1, df);
}

/**
 * fdist:
 * @x: the cutoff point in the distribution.
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * 
 * Returns: the probability a random variable disstributed as
 * F(@dfn, @dfd) is greater than @x, or -1 if either @dfn or @dfd is
 * negative.
 */

double fdist (const double x, const int dfn, const int dfd)
{
    int ia, ib, i, j;
    double xi, fia, fi1, fi2, w, zz, p, zy, d, zk;

    if (dfn <= 0 || dfd <= 0) return -1;
    if (x <= 0) return 1.0;
    ia = 2*(dfn/2)-dfn+2;
    ib = 2*(dfd/2)-dfd+2;
    w = x*dfn/dfd;
    zz = 1./(1.+w);
    if (ia == 1) {
	if(ib == 1) {
	    p = sqrt(w);
	    zy = .3183098862;
	    d = zy*zz/p;
	    p = 2.0*zy*atan(p);
	} else {
	    p = sqrt(w*zz);
	    d = .5*p*zz/w;
	}
    } else {
	if(ib == 1) {
	    p = sqrt(zz);
	    d = .5*zz*p;
	    p = 1.0-p;
	} else {
	    d = zz*zz;
	    p = w*zz;
	}
    }
    zy = 2.*w/zz;
    if (ia == 1) {
	fia = (double) ia;
	for (i=ib+2; i<=dfd; i=i+2) {
	    fi1 = (double) (i-1);
	    fi2 = (double) (i-2);
	    d = (1.0+fia/fi2)*d*zz;
	    p = p+d*zy/fi1;
	}
    } else {
	xi = (double) ((dfd - 1)/2);
	zk = pow(zz,xi);
	d = d*zk*dfd/ib;
	p = p*zk+w*zz*(zk-1.0)/(zz-1.0);
    }
    zy = w*zz;
    zz = 2.0/zz;
    ib = dfd-2;
    for (i=ia+2; i<=dfn; i=i+2) {
	j = i+ib;
	d = zy*d*j/(i - 2);
	p = p - zz*d/j;
    }
    if (p <= 0.0) p = 0.0;
    if (p >= 1.0) p = 1.0;
    return (1.0 - p);
}

/**
 * chisq:
 * @x: the cutoff point in the distribution.
 * @df: degrees of freedom.
 * 
 * Returns: the probability that a random variable distributed as
 * Chi-squared(@df) is greater than @x.
 *
 */

double chisq (const double x, const int df)
{
    double aa, bb, absx, p, zy, zz, h, d, sum, xs, xh, xi, xx, x2;
    int i, m, iy;

    if (x <= 0.0 || df <= 0) return 1.0;
    if (df == 1) {
	zy = sqrt(x);
	p = 2.0 * normal(zy);
	return p;
    }
    h = (double) df;
    if (df <= 30) {
	x2 = -x/2.0;
	zy = h/2.0;
	iy = (int) zy;
	zy = zy - iy;
	if (zy > 0.01) {
	    aa = (h - .999)/2.0;
	    m = (int) aa;
	    d = 1.0;
	    zy = sqrt(x);
	    sum = 0.0;
	    for (i=1; i<=m; i++) {
		xi = (double) i;
		d = d*(2.0*xi-1.0);
		sum += (pow(x, xi) / (zy*d));
	    }
	    xs = normal(zy);
	    p = 2.*xs + 0.7978845612587234 * sum * exp(x2);
	} else {
	    aa = (h - 1.999)/2.0;
	    m = (int) aa; 
	    d = 1.0;
	    sum = 0.0;
	    for (i=1; i<=m; i++) {
		xi = (double) i;
		d *= 2.0*i;
		sum += pow(x, xi) / d;
	    }
	    if (x <= 175.0) p = (1.0 + sum) * exp(x2);
	    else p = 0.0;
	}
	return p;
    }
    xs = h - 1.0;
    d = x - h + 2.0/3.0 - .08/h;
    xx = 2.0 * xs;
    if (floateq(x, xs)) zz = - (1.0/3.0 + .08/h) / sqrt(xx);
    else {
	xh = x - xs;
	absx = (xh>0.0)? xh : -xh;
	aa = xs / x;
	bb = xs * log(aa) + x - xs;
	zz = d * sqrt(bb)/absx;
    }
    p = normal(zz);
    return p;
}

/**
 * normal:
 * @x: the cutoff point in the distribution.
 * 
 * Returns: the probability that a random variable distributed as
 * N(0, 1) is greater than @x.
 *
 */

double normal (const double x)
{
    const double a1 = .0705230784;
    const double a2 = .0422820123;
    const double a3 = .0092705272;
    const double a4 = .0001520143;
    const double a5 = .0002765672;
    const double a6 = .0000430638;
    double absx, xx, zz, p;

    absx = (x > 0.0)? x : -x;
    if (absx <= 14.14) zz = .7071067812*absx;
    else zz = 10.0;
    xx = a6*zz + a5;
    xx = xx*zz + a4;
    xx = xx*zz + a3;
    xx = xx*zz + a2;
    xx = xx*zz + a1;
    xx = xx*zz + 1.0;
    p = 0.5 * pow(xx, -16.0);
    if (x > 0.0) p = 1.0 - p;
    return (1.0 - p);
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
		     const double *Z, const DATAINFO *pdinfo, 
                     print_t *prn)
{
    int i, df1 = 0, df2 = 0;
    char stat;
    double xx = NADBL, mean = 0, variance = 0, xval = 0;
    char cmd[7], df1str[9], df2str[9], fstr[9]; 

    for (;;) {
    if (sscanf(str, "%s %c %s %s %s", cmd, &stat, df1str, df2str, fstr) == 5)
	break;
    if (sscanf(str, "%s %c %s %s", cmd, &stat, df1str, fstr) == 4)
	break;
    if (sscanf(str, "%s %c %s", cmd, &stat, fstr) == 3)
	break;
    else break;
    }

    if (isalpha((unsigned char) df1str[0])) {
	for (i=0; i<pdinfo->v; i++) {
	    if (strcmp(df1str, pdinfo->varname[i]) == 0) {
		df1 = (int) Z[pdinfo->n * i + 1];
		mean = Z[pdinfo->n * i + 1];
		break;
	    }
	}
    } else {
	df1 = atoi(df1str);
	mean = atof(df1str);
    }

    if (isalpha((unsigned char) df2str[0])) {
	for (i=0; i<pdinfo->v; i++) {
	    if (strcmp(df2str, pdinfo->varname[i]) == 0) {
		df2 = (int) Z[pdinfo->n * i + 1];
		variance = Z[pdinfo->n * i + 1];
		break;
	    }
	}
    } else {
	df2 = atoi(df2str);
	variance = atof(df2str);
    }

    if (isalpha((unsigned char) fstr[0])) {
	for (i=0; i<pdinfo->v; i++) {
	    if (strcmp(fstr, pdinfo->varname[i]) == 0) {
		xval = Z[pdinfo->n * i + pdinfo->t1];
		if (na(xval)) {
		    pprintf(prn, "\nstatistic has missing value code\n");
		    return NADBL;
		}		
		break;
	    }
	}
    } else xval = atof(fstr);

    switch (stat) {

    case '1':
    case 'z':
    case 'n':
	xx = normal(xval);
	if (xx < 0) {
	    pprintf(prn, "\np-value calculation failed\n");
	    return -1;
	}	
	pprintf(prn, "\nStandard normal: area to the %s "
		"of %f = %.4g\n", (xval > 0)? "right": "left", 
		xval, xx);
	pprintf(prn, "(two-tailed value = %.4g; complement = %.4g)\n", 
		2.0 * xx, 1.0 - 2.0 * xx);
	return xx;

    case '2':
    case 't':
	xx = tprob(xval, df1);
	if (xx < 0) {
	    pprintf(prn, "\np-value calculation failed\n");
	    return -1;
	}
	pprintf(prn, "\nt(%d): area to the %s of %f = %.4g\n", 
		df1, (xval > 0)? "right": "left",
		xval, 0.5 * xx);
	pprintf(prn, "(two-tailed value = %.4g; complement = %.4g)\n", 
		xx, 1.0 - xx);
	return xx;
    case '3':
    case 'c':
    case 'x':
    case 'X':
	xx = chisq(xval, df1);
	if (xx < 0) {
	    pprintf(prn, "\np-value calculation failed\n");
	    return -1;
	}
	pprintf(prn, "\nChi-square(%d): area to the right of %f = %.4g\n", 
		df1, xval, xx);
	pprintf(prn, "(to the left: %.4g)\n", 1.0 - xx);
	return xx;

    case '4':
    case 'f':
    case 'F':
	xx = fdist(xval, df1, df2);
	if (xx < 0) {
	    pprintf(prn, "\np-value calculation failed\n");
	    return -1;
	}
	pprintf(prn, "\nF(%d, %d): area to the right of %f = %.4g\n", 
		df1, df2, xval, xx);
	pprintf(prn, "(to the left: %.4g)\n", 1.0 - xx);
	return xx;

    case '5':
    case 'g':
    case 'G':
	xx = _gammadist(mean, variance, xval, 2);
	if (na(xx))
	    pprintf(prn, "\nError computing gamma distribution\n");
	else
	    pprintf(prn, "\nGamma (mean %g, variance %g, shape %g, scale %g):"
		    "\n area to the right of %f = %.4g\n", 
		    mean, variance, mean*mean/variance, variance/mean,
		    xval, 1.0 - xx);
	return xx;

    default:
	pprintf(prn, "\nunrecognized pvalue code\n");
	return NADBL;
    }
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
	printf("\n\nChoose one of the following distributions: "
	       "\n\n\t1) Standard normal\t\t2) Student's t\n\t3) "
	       "Chi-square\t\t\t4) F\n"
	       "\t5) Gamma\n\n"
	       "Enter your choice (a number < 0 to exit gretl, 0 to quit "
	       "menu, or\n1, 2, 3, 4, or 5): ");

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
	    _pnormal();
	    break;
	case 2:		
	    _ptvalue();
	    break;
	case 3:		
	    _pchisq();
	    break;
	case 4:		
	    _pfvalue();
	    break;
	case 5:
	    _pgamma();
	    break;
	default:	
	    puts("\ninvalid choice");
	    break;
	}

	printf("\nDo you want to continue with more pvalues (y or n)? ");
	fflush(stdout);
	fscanf(stdin, "%s", ans);

    } while (ans[0] == 'Y' || ans[0] == 'y');
}

/* ........................................................ */

static void _pnormal (void)
{
    double xx, zx; 

    printf("%s", negval);
    zx = _getvalue();
    if(zx < 0.0) return;
    xx = normal(zx);
    printf("\nFor the standard normal, area (one-tail) to the "
	   "right of %g is ", zx);
    _putxx(xx);
}

/* ........................................................ */

static void _ptvalue (void)
{
    int n;
    double xx, zx, xsq; 

    _enterdf(" ");
    n = (int) _getvalue();
    if (n <= 0) return;
    printf("%s", negval);
    zx = _getvalue();
    if(zx < 0.0) return;
    xsq = zx * zx;
    xx = fdist(xsq, 1, n)/2.0;
    printf("\nFor Student's t(%d), area (one-tail) to the "
	   "right of %g is ", n, zx);
    _putxx(xx);
}

/* ........................................................ */

static void _pchisq (void)
{
    int n;
    double xx, zx; 

    _enterdf(" ");
    n = (int) _getvalue();
    if (n <= 0) return;
    printf("%s", negval);
    zx = _getvalue();
    if(zx < 0.0) return;
    xx = chisq(zx, n);
    printf("\nFor Chi-square(%d), area to the right of %g is ", 
	   n, zx);
    _putxx(xx);
}

/* ........................................................ */

static void _pfvalue (void)
{
    int m, n;
    double xx, zx; 

    _enterdf(" for the numerator ");
    m = (int) _getvalue();
    if (m <= 0) return;
    _enterdf(" for the denominator ");
    n = (int) _getvalue();
    if (n <= 0) return;
    printf("%s", negval);
    zx = _getvalue();
    if (zx < 0.0) return;
    xx = fdist(zx, m, n);
    printf("\nFor F(%d, %d), area to the right of %g is ",
	   m, n, zx);
    _putxx(xx);
}

/* ........................................................ */

static void _pgamma (void) 
{
    double mean, variance;
    double xx, zx; 

    printf("\nEnter the mean: ");
    mean = _getvalue();
    if (mean <= 0) return;
    printf("\nEnter the variance: ");
    variance = _getvalue();
    if (variance <= 0) return;
    printf("%s", negval);
    zx = _getvalue();
    if (zx < 0.0) return;
    xx = 1.0 - _gammadist(mean, variance, zx, 2);
    printf("\nFor Gamma (mean %g, variance %g), area to the right of %g is ",
	   mean, variance, zx);
    _putxx(xx);
}

/* ........................................................ */

static double _getvalue (void)
{
    double aa;

    if ((fscanf(stdin, "%lf", &aa)) == 1) {
	return aa;
    }
    return -1;
}

/* ........................................................ */

static void _enterdf (const char *str)
{
    printf("\nEnter d.f.%s(value <= 0 will exit menu): ", str);
}

/* ........................................................ */

/* Functions relating to the gamma distribution.
   Allin Cottrell (cottrell@wfu.edu), October 2000.
   Draws upon the pascal code specialf.pas by Bent Nielsen.
*/

static const double tolerance = 1e-7;

/* internal functions */

static double gamma_stirling (double x);
static double gamma_12 (double x);
static double gamma_integral (double lambda, double x);
static double gamma_integral_expansion (double lambda, double x);
static double gamma_integral_fraction (double lambda, double x);
static double gammadist_wilson_hilferty (double shape, double scale, double x);

/* exported functions */

double _gamma_func (double x)
{
    double xx;

    if (x > 171) { 
	return NADBL; 
    }
    if (x >= 6)
	xx = gamma_stirling(x);
    else if (x > 2)
	xx = (x-1) * _gamma_func(x-1);
    else if (x < 1)
	xx = _gamma_func(x+1)/x;
    else 
	xx = gamma_12(x);
    return xx;
} 

/* ........................................................ */

double _gammadist (double s1, double s2, double x, int control)
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
	xx /= _gamma_func(shape);
    }
    return xx;
}

/* end exported functions */

static double gamma_stirling (double x)
     /* Stirling's series, only to be applied for large x.
	Abramowitz and Stegun, p. 257 
     */
{
    double x2, xx;

    x2 = x*x;
    xx = 1 + (double)(1)/12/x + (double)(1)/288/x2
	- (double)(139)/51840/x/x2 
	- (double)(571)/2488320/(x2*x2);
    xx *= exp(-x + (x-(double)(1)/2) *log(x))*sqrt(2*M_PI);
    return xx;
}

static double gamma_12 (double x)
     /* Polynomial approximation, 6.1.36, 
	Abramowitz and Stegun, p. 257
	error < 3e-7. 
     */
{
    const double b1 = -0.577191652;
    const double b2 = 0.988205891;
    const double b3 = -0.897056937;
    const double b4 = 0.918206857;
    const double b5 = -0.756704078;
    const double b6 = 0.482199394;
    const double b7 = -0.193527818;
    const double b8 = 0.035868343;
    double xx, x1, x2, x4;

    x1 = x-1;
    x2 = x1*x1;
    x4 = x2*x2;
    xx = 1 + b1*x1 + b2*x2 + b3*x1*x2 + b4*x4 + b5*x4*x1;
    xx += b6*x4*x2 + b7*x4*x2*x1 + b8*x4*x4;
    return xx;
}

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
    xx = _gamma_func(lambda);
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











