/* Driver for James MacKinnons's "urcval" function, to calculate
   p-values for unit root tests.

   See James G. MacKinnon, "Numerical Distribution Functions for Unit
   Root and Cointegration Tests", Journal of Applied Econometrics,
   Vol. 11, No. 6, 1996, pp. 601-618, and also 

   http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

   The calculation code here is Copyright (c) James G. MacKinnon, 
   1996 (corrected 2003-5-5).

   This "wrapper" written by Allin Cottrell, 2004.
*/

#include "libgretl.h"
#include "version.h"

#define URDEBUG 0

#if URDEBUG
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

/* Copyright (c) James G. MacKinnon, 1995. Routine to evaluate
   response surface for specified betas and sample size. 
*/

static double eval_crit (double *b, int model, int nreg, int nobs)
{
    double d, cval = 0.0;

    if (nobs == 0) {
	cval = b[0];
    } else if (model == 2) {
	d = 1.0 / nobs;
	cval = b[0] + b[1] * d + b[2] * (d * d);
    } else if (model == 3) {
	d = 1.0 / nobs;
	cval = b[0] + b[1] * d + b[2] * (d * d) + b[3] * (d * d * d);
    } else if (model == 4) {
	d = 1.0 / (nobs - nreg);
	cval = b[0] + b[1] * d + b[2] * (d * d);
    } else if (model == 5) {
	d = 1.0 / (nobs - nreg);
	cval = b[0] + b[1] * d + b[2] * (d * d) + b[3] * (d * d * d);
    } else {
	fputs("*** Warning! Error in input file. ***", stderr);
    }

    return cval;
}

/* Copyright (c) James G. MacKinnon, 1993.  This routine uses the
   Cholesky decomposition to invert a real symmetric matrix.
*/

static int cholx (double *a, int m, int n)
{
    int i, j, k, kl;
    double t, ooa = 0.0;
    int err = 0;

    /* Parameter adjustment */
    a -= 1 + m;

    for (i = 1; i <= n; ++i) {
	kl = i - 1;
	for (j = i; j <= n; ++j) {
	    if (i > 1) {
		for (k = 1; k <= kl; ++k) {
		    a[i + j * m] -= a[k + i * m] * a[k + j * m];
		}
	    } else if (a[i + i * m] <= 0.0) {
		/* error: get out */
		err = i;
		goto cholx_exit;
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
	    ooa = 1.0 / a[j + j * m];
	    if (i >= j) {
		t = 1.0;
	    } else {
		kl = j - 1;
		t = 0.0;
		for (k = i; k <= kl; ++k) {
		    t -= a[i + k * m] * a[k + j * m];
		}
	    }
	    a[i + j * m] = t * ooa;
	}
    }

    for (i = 1; i <= n; ++i) {
	for (j = i; j <= n; ++j) {
	    t = 0.0;
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

/* Copyright (c) James G. MacKinnon, 1995.  Subroutine to do GLS
   estimation the obvious way.  Use only when sample size is small
   (nobs <= 50). 1995-1-3 
*/

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

/* Copyright (c) James G. MacKinnon, 1995.
   Routine to find P value for any specified test statistic. 
*/

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

    /* Parameter adjustments */
    --prob;
    --wght;
    --cnorm;
    beta -= 5;

    /* first, compute all the estimated critical values */

    for (i = 1; i <= URCLEN; ++i) {
	crits[i - 1] = eval_crit(&beta[(i << 2) + 1], model, nreg, nobs);
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
	/* imin is not too close to the end. 
	   Use np points around stat. 
	*/
	for (i=1; i<=np; i++) {
	    ic = imin - nph - 1 + i;
	    yvec[i - 1] = cnorm[ic];
	    xmat[i - 1] = 1.0;
	    xmat[i + 19] = crits[ic - 1];
	    xmat[i + 39] = xmat[i + 19] * crits[ic - 1];
	    xmat[i + 59] = xmat[i + 39] * crits[ic - 1];
	}

	/* form omega matrix */
	for (i=1; i<=np; i++) {
	    for (j=i; j<=np; j++) {
		ic = imin - nph - 1 + i;
		jc = imin - nph - 1 + j;
		top = prob[ic] * (1. - prob[jc]);
		bot = prob[jc] * (1. - prob[ic]);
		omega[i + j * 20 - 21] = wght[ic] * wght[jc] * sqrt(top / bot);
	    }
	}
	for (i=1; i<=np; i++) {
	    for (j=i; j<=np; j++) {
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
	} else {
	    gls(xmat, yvec, omega, gamma, xomx, fits, resid, &ssr, &ssrt, 
		np, 3, 20, 4, 1);
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1);
	}
	pval = normal_cdf(crfit);
    } else {
	/* imin is close to one of the ends. Use points from 
	   imin +/- nph to end. 
	*/
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
		    /* Avoid numerical singularities at the upper end */
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
	    pval = normal_cdf(crfit);
	} else {
	    gls(xmat, yvec, omega, gamma, xomx, fits, resid, &ssr, &ssrt,
		np1, 3, 20, 4, 1);
	    d1 = stat;
	    crfit = gamma[0] + gamma[1] * stat + gamma[2] * (d1 * d1);
	    pval = normal_cdf(crfit);
	}

	/* check that nothing crazy has happened at the ends */
	if (imin == 1 && pval > prob[1]) {
	    pval = prob[1];
	}
	if (imin == URCLEN && pval < prob[URCLEN]) {
	    pval = prob[URCLEN];
	}
    }

    return pval;
}

static char *read_double_and_advance (double *val, char *s)
{
    char valstr[16];

    while (isspace((unsigned char) *s)) s++;

    sscanf(s, "%s", valstr);
    *val = atof(valstr);
    s += strlen(valstr);

    return s;
}

/* 
   niv = # of integrated variables
   itv = appropriate ur_code for nc, c, ct, ctt models
   nobs = sample size (0 for asymptotic)
   arg = test statistic
   pval = location to receive P-value
*/

static int urcval (int niv, int itv, int nobs, double arg, 
		   const char *path, double *pval)
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

    if (itv < 1 || itv > 4) {
	/* these limits correspond to UR_NO_CONST and UR_QUAD_TREND
	   in lib/src/adf_kpss.c */
	return URC_BAD_PARAM;
    }

    /* Open data file */
    sprintf(datfile, "%sdata%curcdata.gz", path, SLASH);
    fz = gretl_gzopen(datfile, "rb");
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

#if URDEBUG
    fprintf(fdb, "code=%s nz=%d, nreg=%d, model=%d, minsize=%d\n",
	    code, urc.nz, urc.nreg, urc.model, urc.minsize);
    fflush(fdb);
#endif

    nvar = (urc.model == 2 || urc.model == 4)? 3 : 4;

    for (i=1; i<=URCLEN; i++) {
	char *s = gzgets(fz, line, sizeof line);

	for (j=1; j<=nvar; j++) {
	    s = read_double_and_advance(&urc.beta[j + (i << 2) - 5], s);
	}
	read_double_and_advance(&urc.wght[i - 1], s);
    }

    /* read from embedded "probs.tab" */
    gzseek(fz, (z_off_t) urc_offsets[MAXVARS], SEEK_SET);
    for (i=0; i<URCLEN; i++) {
	gzgets(fz, line, sizeof line);
	sscanf(line, "%lf %lf", &urc.probs[i], &urc.cnorm[i]);
    }

    if (nobs > 0 && nobs < urc.minsize) {
	urc_ret = URC_SMALL_SAMPLE;
    }

    *pval = fpval(urc.beta, urc.cnorm, urc.wght, urc.probs,
		  arg, nobs, urc.model, urc.nreg);

    gzclose(fz);

    return urc_ret;
}

/* 
   tau = test statistic
   n = sample size (or 0 for asymptotic)
   niv = # of integrated variables
   itv = 1, 2, 3, or 4 for nc, c, ct, ctt models.
   path = path to urc data files
   
   returns: the computed P-value
*/

double mackinnon_pvalue (double tau, int n, int niv, int itv, char *path)
{
    double pval = NADBL;
    int err;

#if URDEBUG
    fdb = fopen("debug.txt", "w");
    fprintf(fdb, "mackinnon_pvalue: tau=%g, n=%d, niv=%d, itv=%d\n",
	    tau, n, niv, itv);
    fprintf(fdb, "mackinnon_pvalue: path='%s'\n", path);
    fflush(fdb);
#endif

    gretl_push_c_numeric_locale();
    err = urcval(niv, itv, n, tau, path, &pval);
    gretl_pop_c_numeric_locale();

#if URDEBUG
    fclose(fdb);
#endif

    if (err == URC_NOT_FOUND) {
	path[0] = '\0';
    }

    if (err != URC_OK && err != URC_SMALL_SAMPLE) {
	pval = NADBL;
    }

    return pval;
}

