/* Driver plugin for James MacKinnons's "urcval" function, to calculate
   p-values for Dickey-Fuller tests.

   See James G. MacKinnon, "Numerical Distribution Functions for Unit
   Root and Cointegration Tests", Journal of Applied Econometrics,
   Vol. 11, No. 6, 1996, pp. 601-618.

   The calculation code here is Copyright (c) James G. MacKinnon.
   The "wrapper" is written by Allin Cottrell.

   The plugin requires MacKinnpn's datafiles, probs.tab and urc-1.tab,
   available in tabs-unix.zip, at

   http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

   The datafiles should be placed in a subdirectory named "data" of
   the gretl plugin directory.

*/

#include "libgretl.h"
#include "f2c.h"

enum urc_errs {
    URC_OK,
    URC_BAD_PARAM,
    URC_NOT_FOUND,
    URC_SMALL_SAMPLE
};

#define MAXVARS 8
#define URCLEN 221

/* The following code translated from James MacKinnon's FORTRAN file
   urcrouts.f, with the help of f2c (version 20030306).
*/

/* private FORTRAN-derived functions */

static double fpval (double *beta, double *cnorm, double *wght, 
		     double *probs, double stat, 
		     double precrt, int nobs, int model, 
		     int nreg, int np);

static double eval (double *beta, int model, int nreg, int nobs);

static int gls (double *xmat, double *yvect, double *omega, 
		double *beta, double *xomx, double *fits, 
		double *resid, double *ssr, double *ssrt, int nobs, 
		int nvar, int nomax, int nvmax, int ivrt);

static int cholx (double *amat, int m, int n, int *kxx);

static double ddnor (double ystar);

/* Copyright (c) James G. MacKinnon, 1996 (corrected 2003-5-5) */

/* urcrouts.f: This is a set of subroutines to estimate critical values
   and P values for unit root and cointegration tests. It is written in
   Fortran 77. Simply call urcval, specifying the first seven arguments.
   The result comes back in the last argument.

   These routines and the associated data files may be used freely for
   non-commercial purposes, provided that proper attribution is made.
   Please cite the paper

   James G. MacKinnon, "Numerical distribution functions for unit root
   and cointegration tests," Journal of Applied Econometrics, 11,
   1996, 601-618.

   The routines and data files may not be incorporated into any book
   or computer program without the express, written consent of the author.

   The routines must have access to the files probs.tab and urc-#.tab
   for # = 1, 2, 3 ..., 12.
*/

int urcval (int niv, int itv, int nobs, double arg, 
	    double *val, const char *path)
{
    FILE *fp;
    char line[80];

    /* return value */
    int urc_ret = URC_OK;

    int i, j;
    int nvar;
    int iskip;
    char datfile[FILENAME_MAX];

    struct {
	double probs[URCLEN], cnorm[URCLEN], beta[884], wght[URCLEN];
	int nz, nreg, model, minsize;
    } urc;

    int urc_offsets[] = {
	39, 
	60685,
	121331,
	178662,
	239308,
	303269,
	360600,
	427876,
	481892
    };

    /* 
       niv = # of integrated variables
       itv = 1, 2, 3, 4 for nc, c, ct, ctt
       nobs = sample size (0 for asymptotic)
       arg = test statistic
       val = P value (returned by routine)
    */

    /* Check that parameters are valid. */
    if (niv < 1 || niv > MAXVARS ||
	itv < 1 || itv > 4) {
	return URC_BAD_PARAM;
    }

    /* Open data file */
    sprintf(datfile, "%sdata%curcdata", path, SLASH);
    fp = fopen(datfile, "r");
    if (fp == NULL) {
	return URC_NOT_FOUND;
    }

    /* skip to appropriate location in file */
    fseek(fp, (long) urc_offsets[niv - 1], SEEK_SET);

    fgets(line, sizeof line, fp);
    sscanf(line, "%*s %d %d %d %d", &urc.nz, &urc.nreg,
	   &urc.model, &urc.minsize);

#if 0
    printf("nz=%d, nreg=%d, model=%d, minsize=%d\n",
	   urc.nz, urc.nreg, urc.model, urc.minsize);
#endif

    if (urc.model == 2 || urc.model == 4) {
	nvar = 3;
    } else {
	nvar = 4;
    }

    for (i = 1; i <= URCLEN; i++) {
	for (j = 1; j <= nvar; j++) {
	    fscanf(fp, "%lf", &urc.beta[j + (i << 2) - 5]);
	}
	fscanf(fp, "%lf", &urc.wght[i - 1]);
    }

    /* read from embedded "probs.tab" */
    fseek(fp, (long) urc_offsets[MAXVARS], SEEK_SET);
    for (i = 0; i < URCLEN; i++) {
	fscanf(fp, "%lf", &urc.probs[i]);
	fscanf(fp, "%lf", &urc.cnorm[i]);
    }

    if (nobs > 0 && nobs < urc.minsize) {
	urc_ret = URC_SMALL_SAMPLE;
    }

    *val = fpval(urc.beta, urc.cnorm, urc.wght, urc.probs,
		 arg, 2.0, nobs, urc.model, urc.nreg, 9);

    fclose(fp);

    return urc_ret;
}

static double fpval (double *beta, double *cnorm, double *wght, 
		     double *probs, double stat, 
		     double precrt, int nobs, int model, 
		     int nreg, int np)
{
    /* System generated locals */
    double d1, d2;

    /* Local variables */
    static int i, j, ic, jc;
    static double sd4;
    static int np1;
    static double bot;
    static int nph;
    static double top, ssr, diff;
    static int imin;
    static double fits[20], xmat[80], xomx[16],	
	ssrt, gamma[4], diffm, omega[400],
	resid[20], crfit;
    static double crits[URCLEN], yvect[20];
    static int nptop;
    static double ttest;

    double pval = 0.0;

    /* Copyright (c) James G. MacKinnon, 1995.
       Routine to find P value for any specified test statistic. 
    */

    /* first, compute all the estimated critical values */

    /* Parameter adjustments */
    --probs;
    --wght;
    --cnorm;
    beta -= 5;

    /* Function Body */
    for (i = 1; i <= URCLEN; ++i) {
	crits[i - 1] = eval(&beta[(i << 2) + 1], model, nreg, nobs);
    }

    /* find critical value closest to test statistic */
    diffm = 1.0e3;
    imin = 0;
    for (i = 1; i <= URCLEN; ++i) {
	diff = (d1 = stat - crits[i - 1], abs(d1));
	if (diff < diffm) {
	    diffm = diff;
	    imin = i;
	}
    }

    nph = np / 2;
    nptop = URCLEN - nph;
    if (imin > nph && imin < nptop) {

	/* imin is not too close to the end. Use np points around stat. */
	for (i = 1; i <= np; ++i) {
	    ic = imin - nph - 1 + i;
	    yvect[i - 1] = cnorm[ic];
	    xmat[i - 1] = 1.;
	    xmat[i + 19] = crits[ic - 1];
	    xmat[i + 39] = xmat[i + 19] * crits[ic - 1];
	    xmat[i + 59] = xmat[i + 39] * crits[ic - 1];
	}

	/* form omega matrix */
	for (i = 1; i <= np; ++i) {
	    for (j = i; j <= np; ++j) {
		ic = imin - nph - 1 + i;
		jc = imin - nph - 1 + j;
		top = probs[ic] * (1. - probs[jc]);
		bot = probs[jc] * (1. - probs[ic]);
		omega[i + j * 20 - 21] = 
		    wght[ic] * wght[jc] * sqrt(top / bot);
	    }
	}
	for (i = 1; i <= np; ++i) {
	    for (j = i; j <= np; ++j) {
		omega[j + i * 20 - 21] = omega[i + j * 20 - 21];
	    }
	}

	gls(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, np, 
	    4, 20, 4, 0);

	/* check to see if gamma(4) is needed */
	sd4 = sqrt(ssrt / (np - 4) * xomx[15]);
	ttest = abs(gamma[3]) / sd4;
	if (ttest > precrt) {
	    /* Computing 2nd power */
	    d1 = stat;
	    /* Computing 3rd power */
	    d2 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1) + 
		gamma[3] * (d2 * (d2 * d2));
	    return ddnor(crfit);
	} else {
	    gls(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, 
		np, 3, 20, 4, 1);
	    /* Computing 2nd power */
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1);
	    return ddnor(crfit);
	}
    } else {

	/* imin is close to one of the ends. Use points from imin +/- nph to end. */
	if (imin < np) {
	    np1 = imin + nph;
	    if (np1 < 5) {
		np1 = 5;
	    }
	    for (i = 1; i <= np1; ++i) {
		yvect[i - 1] = cnorm[i];
		xmat[i - 1] = 1.;
		xmat[i + 19] = crits[i - 1];
		xmat[i + 39] = xmat[i + 19] * crits[i - 1];
		xmat[i + 59] = xmat[i + 39] * crits[i - 1];
	    }
	} else {
	    np1 = (URCLEN + 1) - imin + nph;
	    if (np1 < 5) {
		np1 = 5;
	    }
	    for (i = 1; i <= np1; ++i) {
		ic = (URCLEN + 1) - i;
		yvect[i - 1] = cnorm[ic];
		xmat[i - 1] = 1.;
		xmat[i + 19] = crits[ic - 1];
		xmat[i + 39] = xmat[i + 19] * crits[ic - 1];
		xmat[i + 59] = xmat[i + 39] * crits[ic - 1];
	    }
	}

	/* form omega matrix */
	for (i = 1; i <= np1; ++i) {
	    for (j = i; j <= np1; ++j) {
		if (imin < np) {
		    top = probs[i] * (1. - probs[j]);
		    bot = probs[j] * (1. - probs[i]);
		    omega[i + j * 20 - 21] = wght[i] * wght[j] * sqrt(top / bot);
		} else {
		    /* This is to avoid numerical singularities at the upper end */
		    omega[i + j * 20 - 21] = 0.;
		    if (i == j) {
			omega[i + i * 20 - 21] = 1.;
		    }
		}
	    }
	}
	for (i = 1; i <= np1; ++i) {
	    for (j = i; j <= np1; ++j) {
		omega[j + i * 20 - 21] = omega[i + j * 20 - 21];
	    }
	}

	gls(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, np1, 
	    4, 20, 4, 0);

	/* check to see if gamma(4) is needed */
	sd4 = sqrt(ssrt / (np1 - 4) * xomx[15]);
	ttest = abs(gamma[3]) / sd4;
	if (ttest > precrt) {
	    /* Computing 2nd power */
	    d1 = stat;
	    /* Computing 3rd power */
	    d2 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1) + 
		gamma[3] * (d2 * (d2 * d2));
	    pval = ddnor(crfit);
	} else {
	    gls(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt,
		np1, 3, 20, 4, 1);
	    /* Computing 2nd power */
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1);
	    pval = ddnor(crfit);
	}

	/* check that nothing crazy has happened at the ends */
	if (imin == 1 && pval > probs[1]) {
	    pval = probs[1];
	}
	if (imin == URCLEN && pval < probs[URCLEN]) {
	    pval = probs[URCLEN];
	}
	return pval;
    }

    return pval;
}

static double eval (double *beta, int model, int nreg, int nobs)
{
    double d1, d2;
    double onobs;
    double cval = 0.0;

    /* Copyright (c) James G. MacKinnon, 1995. Routine to evaluate
       response surface for specified betas and sample size. 
    */

    /* Parameter adjustments */
    --beta;

    /* Function Body */
    if (nobs == 0) {
	cval = beta[1];
    } else if (model == 2) {
	onobs = 1. / nobs;
	/* Computing 2nd power */
	d1 = onobs;
	cval = beta[1] + beta[2] * onobs + beta[3] * (d1 * d1);
    } else if (model == 3) {
	onobs = 1. / nobs;
	/* Computing 2nd power */
	d1 = onobs;
	/* Computing 3rd power */
	d2 = onobs;
	cval = beta[1] + beta[2] * onobs + beta[3] * (d1 * d1) 
	    + beta[4] * (d2 * (d2 * d2));
    } else if (model == 4) {
	onobs = 1. / (nobs - nreg);
	/* Computing 2nd power */
	d1 = onobs;
	cval = beta[1] + beta[2] * onobs + beta[3] * (d1 * d1);
    } else if (model == 5) {
	onobs = 1. / (nobs - nreg);
	/* Computing 2nd power */
	d1 = onobs;
	/* Computing 3rd power */
	d2 = onobs;
	cval = beta[1] + beta[2] * onobs + 
	    beta[3] * (d1 * d1) + beta[4] * (d2 * (d2 * d2));
    } else {
	fputs("*** Warning! Error in input file. ***", stderr);
    }

    return cval;
}

static int gls (double *xmat, double *yvect, double *omega, 
		double *beta, double *xomx, double *fits, 
		double *resid, double *ssr, double *ssrt, int nobs, 
		int nvar, int nomax, int nvmax, int ivrt)
{
    double d1;
    int omega_offset = 1 + nomax;
    int xomx_offset = 1 + nvmax;
    int xmat_offset = 1 + nomax;
    int i, j, k, l, kxx;
    static double xomy[50];

    /* Copyright (c) James G. MacKinnon, 1995.  Subroutine to do GLS
       estimation the obvious way.  Use only when sample size is small
       (nobs <= 50). 1995-1-3 
    */

    /* xomx is covariance matrix of parameter estimates if omega is
       truly known.  First, invert omega matrix if ivrt=0. Original one
       gets replaced. */

    /* Parameter adjustments */
    --resid;
    --fits;
    omega -= omega_offset;
    --yvect;
    xomx -= xomx_offset;
    --beta;
    xmat -= xmat_offset;

    if (ivrt == 0) {
	cholx(&omega[omega_offset], nomax, nobs, &kxx);
    }

    /* form xomx matrix and xomy vector */

    for (j = 1; j <= nvar; ++j) {
	xomy[j - 1] = 0.;
	for (l = j; l <= nvar; ++l) {
	    xomx[j + l * nvmax] = 0.;
	}
    }

    for (i = 1; i <= nobs; ++i) {
	for (k = 1; k <= nobs; ++k) {
	    for (j = 1; j <= nvar; ++j) {
		xomy[j - 1] += xmat[i + j * nomax] * 
		    omega[k + i * nomax] * yvect[k];
		for (l = j; l <= nvar; ++l) {
		    xomx[j + l * nvmax] += xmat[i + j * nomax] * 
			omega[k + i * nomax] * 
			xmat[k + l * nomax];
		}
	    }
	}
    }

    for (j = 1; j <= nvar; ++j) {
	for (l = j; l <= nvar; ++l) {
	    xomx[l + j * nvmax] = xomx[j + l * nvmax];
	}
    }

    /* invert xomx matrix */
    cholx(&xomx[xomx_offset], nvmax, nvar, &kxx);

    /*  now form estimates of beta. */
    for (i = 1; i <= nvar; ++i) {
	beta[i] = 0.;
	for (j = 1; j <= nvar; ++j) {
	    beta[i] += xomx[i + j * nvmax] * xomy[j - 1];
	}
    }

    /* find ssr, fitted values, and residuals */
    *ssr = 0.;
    for (i = 1; i <= nobs; ++i) {
	fits[i] = 0.;
	for (j = 1; j <= nvar; ++j) {
	    fits[i] += xmat[i + j * nomax] * beta[j];
	}
	resid[i] = yvect[i] - fits[i];
	/* Computing 2nd power */
	d1 = resid[i];
	*ssr += d1 * d1;
    }

    /* find ssr from transformed regression */
    *ssrt = 0.;
    for (i = 1; i <= nobs; ++i) {
	for (k = 1; k <= nobs; ++k) {
	    *ssrt += resid[i] * omega[k + i * nomax] * resid[k];
	}
    }

    return 0;
}

static int cholx (double *amat, int m, int n, int *kxx)
{
    int i, j, k;
    double t;
    int kl;
    double ooa = 0.0;

    /* Copyright (c) James G. MacKinnon, 1993.  This routine uses the
       cholesky decomposition to invert a real symmetric matrix.
    */

    /* Parameter adjustment */
    amat -= 1 + m;

    *kxx = 0;
    for (i = 1; i <= n; ++i) {
	kl = i - 1;
	for (j = i; j <= n; ++j) {
	    if (i > 1) {
		for (k = 1; k <= kl; ++k) {
		    amat[i + j * m] -= 
			amat[k + i * m] * amat[k + j * m];
		}
	    } else {
		if (amat[i + i * m] <= 0.0) {
		    *kxx = i;
		    goto L20;
		}
	    }
	    if (i == j) {
		amat[i + i * m] = sqrt(amat[i + i * m]);
	    } else {
		if (j == i + 1) {
		    ooa = 1. / amat[i + i * m];
		}
		amat[i + j * m] *= ooa;
	    }
	}
    }

    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    ooa = 1. / amat[j + j * m];
	    if (i >= j) {
		t = 1.;
		goto L12;
	    }
	    kl = j - 1;
	    t = 0.;
	    for (k = i; k <= kl; ++k) {
		t -= amat[i + k * m] * amat[k + j * m];
	    }
L12:
	    amat[i + j * m] = t * ooa;
	}
    }

    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    t = 0.;
	    for (k = j; k <= n; ++k) {
		t += amat[i + k * m] * amat[j + k * m];
	    }
	    amat[i + j * m] = t;
	    amat[j + i * m] = t;
	}
    }

L20:
    return 0;
}

static double ddnor (double ystar)
{
    double c[5] = { 3209.377589138469472562,
		    377.4852376853020208137,
		    113.8641541510501556495,
		    3.161123743870565596947,
		    .185777706184603152673 
    };

    double d[4] = { 2844.236833439170622273,
		    1282.616526077372275645,
		    244.0246379344441733056,
		    23.60129095234412093499 
    };

    double orpi = .5641895835477562869483;

    double root2 = .70710678118654752440083;

    double p[6] = { -6.58749161529837803157e-4,
		    -.0160837851487422766278,
		    -.125781726111229246204,
		    -.360344899949804439429,
		    -.305326634961232344035,
		    -.0163153871373020978498 
    };

    double q[5] = { .00233520497626869185443,
		    .0605183413124413191178,
		    .527905102951428412248,
		    1.87295284992346047209,
		    2.56852019228982242072 
    };

    double a[9] = { 1230.33935479799725272,
		    2051.07837782607146532,
		    1712.04761263407058314,
		    881.952221241769090411,
		    298.635138197400131132,
		    66.1191906371416294775,
		    8.88314979438837594118,
		    .56418849698867008918,
		    2.15311535474403846343e-8 
    };

    double b[8] = { 1230.33935480374942043,
		    3439.36767414372163696,
		    4362.6190901432471582,
		    3290.79923573345962678,
		    1621.38957456669018874,
		    537.181101862009857509,
		    117.693950891312499305,
		    15.7449261107098347253 
    };

    double x, y, x2, x3, x4, x5, x6, x7, x8, xm2, xm4, xm6, xm8;
    double erf, bot, xm10;
    double top, erfc, crap;
    int isw;

    /* Copyright (c) James G. MacKinnon, 1993.  Routine to evaluate
       cumulative normal distribution Written originally in late
       1970's -- modified 1993 to avoid changing the argument.

       This subroutine uses Cody's method to evaluate the cumulative
       normal distribution. It is probably accurate to 19 or 20
       significant digits. It was written in 1977, based on the Cody
       article referred to in the documentation for IMSL subroutine
       mdnor.
    */

    isw = 1;
    y = ystar;
    if (ystar < -16.) {
	y = -16.;
    }
    if (ystar > 16.) {
	y = 16.;
    }
    x = -y * root2;
    if (x > 0.) {
	goto L1;
    }
    if (x < 0.) {
	goto L2;
    }

    return .5;

 L2:
    x = -x;
    isw = -1;

 L1:
    if (x < .477) {
	goto L10;
    }
    if (x <= 4.) {
	goto L20;
    }

    /* evaluate erfc for x.gt.4.0 */
    x2 = x * x;
    xm2 = 1. / x2;
    xm4 = xm2 * xm2;
    xm6 = xm4 * xm2;
    xm8 = xm4 * xm4;
    xm10 = xm6 * xm4;
    top = p[0] + p[1] * xm2 + p[2] * xm4 + p[3] * xm6 + p[4] * xm8 + p[5] * 
	xm10;
    bot = q[0] + q[1] * xm2 + q[2] * xm4 + q[3] * xm6 + q[4] * xm8 + xm10;
    crap = orpi + top / (bot * x2);
    erfc = exp(-x2) * crap / x;

    if (isw == -1) {
	erfc = 2. - erfc;
    }

    return erfc * .5;
 L20:

    /* evaluate erfc for .477.lt.x.le.4.0 */
    x2 = x * x;
    x3 = x2 * x;
    x4 = x2 * x2;
    x5 = x3 * x2;
    x6 = x3 * x3;
    x7 = x3 * x4;
    x8 = x4 * x4;
    top = a[0] + a[1] * x + a[2] * x2 + a[3] * x3 + a[4] * 
	x4 + a[5] * x5 + a[6] * x6 + a[7] * x7 + a[8] * x8;
    bot = b[0] + b[1] * x + b[2] * x2 + b[3] * x3 + 
	b[4] * x4 + b[5] * x5 + b[6] * x6 + b[7] * x7 + x8;
    erfc = exp(-x2) * top / bot;

    if (isw == -1) {
	erfc = 2. - erfc;
    }

    return erfc * .5;
 L10:

    /* evaluate erf for x.lt..477 */
    x2 = x * x;
    x4 = x2 * x2;
    x6 = x4 * x2;
    x8 = x4 * x4;
    top = c[0] + c[1] * x2 + c[2] * x4 + c[3] * x6 + c[4] * x8;
    bot = d[0] + d[1] * x2 + d[2] * x4 + d[3] * x6 + x8;
    erf = x * top / bot;

    erf *= isw;
    erfc = 1. - erf;

    return erfc * .5;
}

double mackinnon_pvalue (double tval, int n, char *path)
{
    int niv = 1; /* number of variables */
    int itv = 2; /* model "c", with constant */
    double val;
    int check;

    check = urcval(niv, itv, n, tval, &val, path);

    if (check == URC_NOT_FOUND) {
	*path = '\0';
    }

    if (check != URC_OK && check != URC_SMALL_SAMPLE) {
	val = NADBL;
    }

    return val;
}

