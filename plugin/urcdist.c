/* Driver for James MacKinnons's "urcval" function, to calculate
   p-values for unit root tests.

   See James G. MacKinnon, "Numerical Distribution Functions for Unit
   Root and Cointegration Tests", Journal of Applied Econometrics,
   Vol. 11, No. 6, 1996, pp. 601-618, and also 

   http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

   The calculation code here is Copyright (c) James G. MacKinnon, 
   1996 (corrected 2003-5-5).

   The "wrapper" is written by Allin Cottrell, 2004.
*/

#include "libgretl.h"
#include "var.h"

#include <zlib.h>

#undef URDEBUG

#ifdef URDEBUG
FILE *fdb;
#endif

enum urc_errs {
    URC_OK,
    URC_BAD_PARAM,
    URC_NOT_FOUND,
    URC_SMALL_SAMPLE
};

#define MAXVARS 8
#define URCLEN 221
#define BIGLEN 884

static double fpval (double *beta, double *cnorm, double *wght, 
		     double *probs, double stat, int nobs, 
		     int model, int nreg);

static double eval (double *beta, int model, int nreg, int nobs);

static int gls (double *xmat, double *yvect, double *omega, 
		double *beta, double *xomx, double *fits, 
		double *resid, double *ssr, double *ssrt, int nobs, 
		int nvar, int nomax, int nvmax, int ivrt);

static int cholx (double *a, int m, int n);

static double ddnor (double ystar);

static char *read_double_and_advance (double *val, char *s)
{
    char valstr[16];

    while (isspace((unsigned char) *s)) s++;

    sscanf(s, "%s", valstr);
    *val = atof(valstr);
    s += strlen(valstr);

    return s;
}

#if 0 /* work in progress here */

/* npr = number of restrictions (p - r) 
   itt = 1 or 2 for lambda-max test or trace test
   itv = 0, 1, 2, 3, or 4 for cases 0, 1, 2, 1*, or 2*.
   arg = test statistic
   val = P value (returned by routine) 
*/

static int johval (int npr, int itt, int itv, 
		   double arg, double *val)
{
    gzFile fz;
    int i1, i2[2];

    int i, ii, np, nx;
    int junk;
    double size;
    int iskip;
    double precrt;

    char datfile[FILENAME_MAX];

    struct {
	double crits[URCLEN];
	double cnorm[URCLEN];
	double beta[BIGLEN];
	double wght[URCLEN];
    } joh;

    /* byte offsets into data : FIXME */
    int joh_offsets[] = {
	39,     /* joh-1 */
	60685,  /* joh-2 */
	121331, /* joh-3 */
	178662, /* joh-4 */
	239308, /* joh-5 */
	303269, /* joh-6 */
	360600, /* joh-7 */
	427876, /* joh-8 */
	481892  /* probs */
    };

    int urc_ret = URC_OK;

    if (npr < 1 || npr > 12) {
	fprintf(stderr, "Number of restrictions must be between 1 and 8\n");
	return URC_BAD_PARAM;
    }

    if (itt < 1 || itt > 2) {
	fprintf(stderr, "The valid options for itt are 1 and 2.\n");
	return URC_BAD_PARAM;
    }

    if (itv < 0 || itv > 4) {
	fprintf(stderr, "The valid options for itv are 0, 1, 2, 3 and 4\n");
	return URC_BAD_PARAM;
    }

    /* Open data file */
    sprintf(datfile, "%sdata%cjohdata.gz", path, SLASH);
    fz = gzopen(datfile, "rb");
    if (fz == NULL) {
	return URC_NOT_FOUND;
    }

    /* skip to appropriate location in data file */
    gzseek(fz, (z_off_t) joh_offsets[npr - 1], SEEK_SET);

    if (itt != 1) {
	iskip = 1110;
    }    
    iskip += itv * (URCLEN + 1);
    for (i = 0; i < iskip; i++) {
	gzgets(fz, line, sizeof line);
    }

    /* johN table */
    for (i = 1; i <= URCLEN; i++) {
	char *s = gzgets(fz, line, sizeof line);

	/* ii = 222 - i; */
	read_double_and_advance(&joh.crits[i - 1], s); /* &crits[ii - 1] */
	read_double_and_advance(&joh.wght[i - 1], s);
    }

    /* skip to probs now */
    for (i = 1; i <= URCLEN; i++) {
	char *s = gzgets(fz, line, sizeof line);

	read_double_and_advance(&joh.probs[i - 1], s);
	read_double_and_advance(&joh.cnorm[i - 1], s);
    }

    np = 11;
    precrt = 2.; /* check args to fpval */
    *val = fpval(joh.crits, joh.cnorm, joh.wght, joh.probs, 
		 arg, 2.0, 11, nx);

    gzclose(fz);

    return urc_ret;
}

#endif /* work in progress */

/* 
   niv = # of integrated variables
   itv = appropriate ur_code for nc, c, ct, ctt models
   nobs = sample size (0 for asymptotic)
   arg = test statistic
   val = P-value (returned by routine)
*/

static int urcval (int niv, int itv, int nobs, double arg, 
		   double *val, const char *path)
{
    gzFile fz;
    char line[80], code[8];

    int urc_ret = URC_OK;

    int i, j, iskip, nvar;
    char datfile[FILENAME_MAX];

    struct {
	double probs[URCLEN];
	double cnorm[URCLEN];
	double beta[BIGLEN];
	double wght[URCLEN];
	int nz, nreg, model, minsize;
    } urc;

    /* byte offsets into data */
    int urc_offsets[] = {
	39,     /* urc-1 */
	60685,  /* urc-2 */
	121331, /* urc-3 */
	178662, /* urc-4 */
	239308, /* urc-5 */
	303269, /* urc-6 */
	360600, /* urc-7 */
	427876, /* urc-8 */
	481892  /* probs */
    };

    /* Check that parameters are valid */
    if (niv < 1 || niv > MAXVARS) {
	return URC_BAD_PARAM;
    }
    if (itv < UR_NO_CONST || itv > UR_TREND_SQUARED) {
	return URC_BAD_PARAM;
    }

    /* Open data file */
    sprintf(datfile, "%sdata%curcdata.gz", path, SLASH);
    fz = gzopen(datfile, "rb");
    if (fz == NULL) {
	return URC_NOT_FOUND;
    }

    /* skip to appropriate location in data file */
    gzseek(fz, (z_off_t) urc_offsets[niv - 1], SEEK_SET);

    iskip = (itv - 1) * (URCLEN + 1);
    for (i = 0; i < iskip; i++) {
	gzgets(fz, line, sizeof line);
    }

    gzgets(fz, line, sizeof line);
    sscanf(line, "%s %d %d %d %d", code, &urc.nz, &urc.nreg,
	   &urc.model, &urc.minsize);

#ifdef URDEBUG
    fprintf(fdb, "code=%s nz=%d, nreg=%d, model=%d, minsize=%d\n",
	    code, urc.nz, urc.nreg, urc.model, urc.minsize);
    fflush(fdb);
#endif

    if (urc.model == 2 || urc.model == 4) {
	nvar = 3;
    } else {
	nvar = 4;
    }

    for (i = 1; i <= URCLEN; i++) {
	char *s = gzgets(fz, line, sizeof line);

	for (j = 1; j <= nvar; j++) {
	    s = read_double_and_advance(&urc.beta[j + (i << 2) - 5], s);
	}
	read_double_and_advance(&urc.wght[i - 1], s);
    }

    /* read from embedded "probs.tab" */
    gzseek(fz, (z_off_t) urc_offsets[MAXVARS], SEEK_SET);
    for (i = 0; i < URCLEN; i++) {
	gzgets(fz, line, sizeof line);
	sscanf(line, "%lf %lf", &urc.probs[i], &urc.cnorm[i]);
    }

    if (nobs > 0 && nobs < urc.minsize) {
	urc_ret = URC_SMALL_SAMPLE;
    }

    *val = fpval(urc.beta, urc.cnorm, urc.wght, urc.probs,
		 arg, nobs, urc.model, urc.nreg);

    gzclose(fz);

    return urc_ret;
}

static double fpval (double *beta, double *cnorm, double *wght, 
		     double *prob, double stat, int nobs, 
		     int model, int nreg)
{
    double d1, precrt = 2.0;
    int i, j, ic, jc, imin = 0;
    int np1, nph, nptop, np = 9;
    double bot, top, ssr, diff;
    double sd4, ttest, crfit;
    double ssrt, diffm = 1000.0; 
    double yvec[20], fits[20], resid[20];
    double xmat[80], xomx[16], gamma[4], omega[400];
    double crits[URCLEN];
    double pval = 0.0;

    /* Copyright (c) James G. MacKinnon, 1995.
       Routine to find P value for any specified test statistic. 
    */

    /* Parameter adjustments */
    --prob;
    --wght;
    --cnorm;
    beta -= 5;

    /* first, compute all the estimated critical values */

    for (i = 1; i <= URCLEN; ++i) {
	crits[i - 1] = eval(&beta[(i << 2) + 1], model, nreg, nobs);
    }

    /* find critical value closest to test statistic */
    for (i = 0; i < URCLEN; i++) {
	diff = fabs(stat - crits[i]);
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
	    yvec[i - 1] = cnorm[ic];
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
		top = prob[ic] * (1. - prob[jc]);
		bot = prob[jc] * (1. - prob[ic]);
		omega[i + j * 20 - 21] = wght[ic] * wght[jc] * sqrt(top / bot);
	    }
	}
	for (i = 1; i <= np; ++i) {
	    for (j = i; j <= np; ++j) {
		omega[j + i * 20 - 21] = omega[i + j * 20 - 21];
	    }
	}

	gls(xmat, yvec, omega, gamma, xomx, fits, resid, &ssr, &ssrt, np, 
	    4, 20, 4, 0);

	/* check to see if gamma(4) is needed */
	sd4 = sqrt(ssrt / (np - 4) * xomx[15]);
	ttest = fabs(gamma[3]) / sd4;
	if (ttest > precrt) {
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1) + 
		gamma[3] * (d1 * d1 * d1);
	    return ddnor(crfit);
	} else {
	    gls(xmat, yvec, omega, gamma, xomx, fits, resid, &ssr, &ssrt, 
		np, 3, 20, 4, 1);
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
		yvec[i - 1] = cnorm[i];
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
		yvec[i - 1] = cnorm[ic];
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
		    top = prob[i] * (1. - prob[j]);
		    bot = prob[j] * (1. - prob[i]);
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

	gls(xmat, yvec, omega, gamma, xomx, fits, resid, &ssr, &ssrt, np1, 
	    4, 20, 4, 0);

	/* check to see if gamma(4) is needed */
	sd4 = sqrt(ssrt / (np1 - 4) * xomx[15]);
	ttest = fabs(gamma[3]) / sd4;
	if (ttest > precrt) {
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1) + 
		gamma[3] * (d1 * d1 * d1);
	    pval = ddnor(crfit);
	} else {
	    gls(xmat, yvec, omega, gamma, xomx, fits, resid, &ssr, &ssrt,
		np1, 3, 20, 4, 1);
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1);
	    pval = ddnor(crfit);
	}

	/* check that nothing crazy has happened at the ends */
	if (imin == 1 && pval > prob[1]) {
	    pval = prob[1];
	}
	if (imin == URCLEN && pval < prob[URCLEN]) {
	    pval = prob[URCLEN];
	}
	return pval;
    }

    return pval;
}

static double eval (double *beta, int model, int nreg, int nobs)
{
    double d, cval = 0.0;

    /* Copyright (c) James G. MacKinnon, 1995. Routine to evaluate
       response surface for specified betas and sample size. 
    */

    /* Parameter adjustments */
    --beta;

    /* Function Body */
    if (nobs == 0) {
	cval = beta[1];
    } else if (model == 2) {
	d = 1. / nobs;
	cval = beta[1] + beta[2] * d + beta[3] * (d * d);
    } else if (model == 3) {
	d = 1. / nobs;
	cval = beta[1] + beta[2] * d + beta[3] * (d * d) 
	    + beta[4] * (d * d * d);
    } else if (model == 4) {
	d = 1. / (nobs - nreg);
	cval = beta[1] + beta[2] * d + beta[3] * (d * d);
    } else if (model == 5) {
	d = 1. / (nobs - nreg);
	cval = beta[1] + beta[2] * d + beta[3] * (d * d) 
	    + beta[4] * (d * d * d);
    } else {
	fputs("*** Warning! Error in input file. ***", stderr);
    }

    return cval;
}

static int gls (double *xmat, double *yvec, double *omega, 
		double *beta, double *xomx, double *fits, 
		double *resid, double *ssr, double *ssrt, int nobs, 
		int nvar, int nomax, int nvmax, int ivrt)
{
    int omega_offset = 1 + nomax;
    int xomx_offset = 1 + nvmax;
    int xmat_offset = 1 + nomax;
    int i, j, k, l;
    double xomy[50];

    /* Copyright (c) James G. MacKinnon, 1995.  Subroutine to do GLS
       estimation the obvious way.  Use only when sample size is small
       (nobs <= 50). 1995-1-3 
    */

    /* xomx is covariance matrix of parameter estimates if omega is
       truly known.  First, invert omega matrix if ivrt=0. Original one
       gets replaced. */

    /* Parameter adjustments */
    omega -= omega_offset;
    xomx -= xomx_offset;
    xmat -= xmat_offset;
    --resid;
    --fits;
    --yvec;
    --beta;

    if (ivrt == 0) {
	cholx(&omega[omega_offset], nomax, nobs);
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
		    omega[k + i * nomax] * yvec[k];
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
    cholx(&xomx[xomx_offset], nvmax, nvar);

    /* form estimates of beta */
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
	resid[i] = yvec[i] - fits[i];
	*ssr += resid[i] * resid[i];
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

static int cholx (double *a, int m, int n)
{
    int i, j, k, kl;
    double t, ooa = 0.0;
    int err = 0;

    /* Copyright (c) James G. MacKinnon, 1993.  This routine uses the
       cholesky decomposition to invert a real symmetric matrix.
    */

    /* Parameter adjustment */
    a -= 1 + m;

    for (i = 1; i <= n; ++i) {
	kl = i - 1;
	for (j = i; j <= n; ++j) {
	    if (i > 1) {
		for (k = 1; k <= kl; ++k) {
		    a[i + j * m] -= a[k + i * m] * a[k + j * m];
		}
	    } else {
		if (a[i + i * m] <= 0.0) {
		    err = i;
		    goto cholx_exit;
		}
	    }
	    if (i == j) {
		a[i + i * m] = sqrt(a[i + i * m]);
	    } else {
		if (j == i + 1) {
		    ooa = 1. / a[i + i * m];
		}
		a[i + j * m] *= ooa;
	    }
	}
    }

    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    ooa = 1. / a[j + j * m];
	    if (i >= j) {
		t = 1.;
		goto cholx_jump;
	    }
	    kl = j - 1;
	    t = 0.;
	    for (k = i; k <= kl; ++k) {
		t -= a[i + k * m] * a[k + j * m];
	    }
	cholx_jump:
	    a[i + j * m] = t * ooa;
	}
    }

    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    t = 0.;
	    for (k = j; k <= n; ++k) {
		t += a[i + k * m] * a[j + k * m];
	    }
	    a[i + j * m] = t;
	    a[j + i * m] = t;
	}
    }

 cholx_exit:

    return err;
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

    double x, x2, x3, x4, x5, x6, x7, x8, xm2, xm4, xm6, xm8;
    double erf, bot, xm10;
    double top, erfc, crap;
    double y = ystar;
    int isw = 1;

    /* Copyright (c) James G. MacKinnon, 1993.  Routine to evaluate
       cumulative normal distribution Written originally in late
       1970's -- modified 1993 to avoid changing the argument.

       This subroutine uses Cody's method to evaluate the cumulative
       normal distribution. It is probably accurate to 19 or 20
       significant digits. It was written in 1977, based on the Cody
       article referred to in the documentation for IMSL subroutine
       mdnor.
    */

    if (ystar < -16.) {
	y = -16.;
    } else if (ystar > 16.) {
	y = 16.;
    }

    x = -y * root2;

    if (x == 0.0) {
	return .5;
    } else if (x < 0.0) {
	x = -x;
	isw = -1;
    }

    /* evaluate erfc for x > 4.0 */
    if (x > 4.0) {
	x2 = x * x;
	xm2 = 1.0 / x2;
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
    } else if (x <= 4.0 && x > .477) {
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
    } else {
	/* evaluate erf for x < .477 */
	x2 = x * x;
	x4 = x2 * x2;
	x6 = x4 * x2;
	x8 = x4 * x4;
	top = c[0] + c[1] * x2 + c[2] * x4 + c[3] * x6 + c[4] * x8;
	bot = d[0] + d[1] * x2 + d[2] * x4 + d[3] * x6 + x8;
	erf = x * top / bot;

	erf *= isw;
	erfc = 1.0 - erf;
    }

    return erfc * .5;
}

double mackinnon_pvalue (double tval, int n, int niv, int itv, char *path)
{
    double val;
    int check;

#ifdef URDEBUG
    fdb = fopen("debug.txt", "w");
    fprintf(fdb, "mackinnon_pvalue: tval=%g, n=%d, niv=%d, itv=%d\n",
	    tval, n, niv, itv);
    fprintf(fdb, "mackinnon_pvalue: path='%s'\n", path);
    fflush(fdb);
#endif

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    check = urcval(niv, itv, n, tval, &val, path);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#ifdef URDEBUG
    fclose(fdb);
#endif

    if (check == URC_NOT_FOUND) {
	path[0] = '\0';
    }

    if (check != URC_OK && check != URC_SMALL_SAMPLE) {
	val = NADBL;
    }

    return val;
}

