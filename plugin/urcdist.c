/* Driver for James MacKinnons's "urcval" function, to calculate
   p-values for unit root tests.

   See James G. MacKinnon, "Numerical Distribution Functions for Unit
   Root and Cointegration Tests", Journal of Applied Econometrics,
   Vol. 11, No. 6, 1996, pp. 601-618, and also

   http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

   The calculation code here is Copyright (c) James G. MacKinnon,
   1996 (corrected 2003-5-5).

   This "wrapper" written by Allin Cottrell, 2004; revised to
   load MacKinnon's tables in binary format, 2015.
*/

#include "libgretl.h"
#include "version.h"
#include "swap_bytes.h"
#include "matrix_extra.h"

#define URDEBUG 0

#if URDEBUG
FILE *fdb;
#endif

#define NIVMAX 8
#define URCLEN 221
#define BIGLEN 884

/* Based on Fortran code copyright (c) James G. MacKinnon,
   1995. Routine to evaluate response surface for specified betas and
   sample size.
*/

static void eval_all_crit (double tau, double *b, int nb,
			   int T, double *crit, int *imin)
{
    double d1 = 0, d2 = 0, d3 = 0;
    double diff, diffm = 1000;
    int i;

    if (T > 0) {
	d1 = 1.0 / T;
	d2 = d1 * d1;
	d3 = d1 * d2;
    }

    for (i=0; i<URCLEN; i++) {
	if (T == 0) {
	    crit[i] = b[0];
	} else if (nb == 3) {
	    crit[i] = b[0] + b[1]*d1 + b[2]*d2;
	} else if (nb == 4) {
	    crit[i] = b[0] + b[1]*d1 + b[2]*d2 + b[3]*d3;
	}
	diff = fabs(tau - crit[i]);
	if (diff < diffm) {
	    diffm = diff;
	    *imin = i + 1; /* has to be 1-based */
	}
	b += nb;
    }
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
		fprintf(stderr, "cholx: failed at i = %d\n", i);
		err = E_NOTPD;
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
		double *resid, double *ssr, double *ssrt,
		int T, int ivrt)
{
    int nomax = 20;
    int nvmax = 4;
    int nvar = 4 - ivrt;
    int omega_offset = 1 + nomax;
    int xomx_offset = 1 + nvmax;
    int xmat_offset = 1 + nomax;
    int i, j, k, l;
    double xomy[50];
    int err = 0;

    /* xomx is covariance matrix of parameter estimates if omega is
       truly known. First, invert omega matrix if ivrt=0. Original one
       gets replaced.
    */

    /* parameter adjustments */
    omega -= omega_offset;
    xomx -= xomx_offset;
    xmat -= xmat_offset;
    --resid;
    --fits;
    --yvec;
    --beta;

    if (ivrt == 0) {
	err = cholx(&omega[omega_offset], nomax, T);
	if (err) {
	    return err;
	}
    }

    /* form xomx matrix and xomy vector */
    for (j = 1; j <= nvar; ++j) {
	xomy[j - 1] = 0.;
	for (l = j; l <= nvar; ++l) {
	    xomx[j + l * nvmax] = 0.;
	}
    }
    for (i = 1; i <= T; ++i) {
	for (k = 1; k <= T; ++k) {
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
    err = cholx(&xomx[xomx_offset], nvmax, nvar);
    if (err) {
	return err;
    }

    /* form estimates of beta */
    for (i = 1; i <= nvar; ++i) {
	beta[i] = 0.0;
	for (j = 1; j <= nvar; ++j) {
	    beta[i] += xomx[i + j * nvmax] * xomy[j - 1];
	}
    }

    /* find ssr, fitted values, and residuals */
    *ssr = 0.0;
    for (i = 1; i <= T; ++i) {
	fits[i] = 0.0;
	for (j = 1; j <= nvar; ++j) {
	    fits[i] += xmat[i + j * nomax] * beta[j];
	}
	resid[i] = yvec[i] - fits[i];
	*ssr += resid[i] * resid[i];
    }

    /* find ssr from transformed regression */
    *ssrt = 0.0;
    for (i = 1; i <= T; ++i) {
	for (k = 1; k <= T; ++k) {
	    *ssrt += resid[i] * omega[k + i * nomax] * resid[k];
	}
    }

    return err;
}

/* Based on Fortran code copyright (c) James G. MacKinnon, 1995.
   Routine to find P-value for any specified test statistic.
*/

static double fpval (double *beta, int nbeta, double *wght,
		     double *prob, double *cnorm,
		     double tau, int T, int *err)
{
    double d1, precrt = 2.0;
    int i, j, ic, jc, imin = 0;
    int np1, nph, nptop, np = 9;
    double bot, top, ssr, ssrt;
    double se3, ttest, crfit;
    double yvec[20], fits[20], resid[20];
    double xmat[80], xomx[16], gamma[4], omega[400];
    double crits[URCLEN];
    double pval = 0.0;

    /* first compute all the estimated critical values,
       and find the one closest to the test statistic,
       indexed by @imin.
    */
    eval_all_crit(tau, beta, nbeta, T, crits, &imin);

    nph = np / 2;
    nptop = URCLEN - nph;

    if (imin > nph && imin < nptop) {
	/* imin is not too close to the end.
	   Use np points around tau.
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
		top = prob[ic] * (1 - prob[jc]);
		bot = prob[jc] * (1 - prob[ic]);
		omega[i + j * 20 - 21] = wght[ic] * wght[jc] * sqrt(top / bot);
	    }
	}
	for (i=1; i<=np; i++) {
	    for (j=i; j<=np; j++) {
		omega[j + i * 20 - 21] = omega[i + j * 20 - 21];
	    }
	}

	*err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
		   &ssr, &ssrt, np, 0);
	if (*err) {
	    goto bailout;
	}

	/* check: is a cubic term actually needed? */
	se3 = sqrt(ssrt / (np - 4) * xomx[15]);
	ttest = fabs(gamma[3]) / se3;
	d1 = tau;

	if (ttest > precrt) {
	    crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1) +
		gamma[3] * (d1 * d1 * d1);
	} else {
	    *err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
		       &ssr, &ssrt, np, 1);
	    if (*err) {
		goto bailout;
	    }
	    crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1);
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
		xmat[i - 1] = 1.0;
		xmat[i + 19] = crits[i - 1];
		xmat[i + 39] = xmat[i + 19] * crits[i-1];
		xmat[i + 59] = xmat[i + 39] * crits[i-1];
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
		xmat[i + 39] = xmat[i + 19] * crits[ic-1];
		xmat[i + 59] = xmat[i + 39] * crits[ic-1];
	    }
	}

	/* form omega matrix */
	for (i = 1; i <= np1; ++i) {
	    for (j = i; j <= np1; ++j) {
		if (imin < np) {
		    top = prob[i] * (1.0 - prob[j]);
		    bot = prob[j] * (1.0 - prob[i]);
		    omega[i + j * 20 - 21] = wght[i] * wght[j] * sqrt(top / bot);
		} else {
		    /* avoid numerical singularities at the upper end */
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

	*err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
		   &ssr, &ssrt, np1, 0);
	if (*err) {
	    goto bailout;
	}

	/* is gamma[3] needed? */
	se3 = sqrt(ssrt / (np1 - 4) * xomx[15]);
	ttest = fabs(gamma[3]) / se3;
	d1 = tau;

	if (ttest > precrt) {
	    crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1) +
		gamma[3] * (d1 * d1 * d1);
	} else {
	    *err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
		       &ssr, &ssrt, np1, 1);
	    if (*err) {
		goto bailout;
	    }
	    crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1);
	}
	pval = normal_cdf(crfit);

	/* check that nothing crazy has happened at the ends */
	if (imin == 1 && pval > prob[1]) {
	    pval = prob[1];
	}
	if (imin == URCLEN && pval < prob[URCLEN]) {
	    pval = prob[URCLEN];
	}
    }

 bailout:

    if (*err) {
	pval = NADBL;
    }

    return pval;
}

#if G_BYTE_ORDER == G_BIG_ENDIAN

static void urc_swap_endianness (double *beta,
				 int nbeta,
				 double *wght,
				 double *prob,
				 double *cnorm)
{
    int i, n = nbeta * URCLEN;

    for (i=0; i<n; i++) {
	reverse_double(beta[i]);
    }

    /* the following arrays are padded by one */

    for (i=1; i<=URCLEN; i++) {
	reverse_double(wght[i]);
    }
    for (i=1; i<=URCLEN; i++) {
	reverse_double(prob[i]);
    }
    for (i=1; i<=URCLEN; i++) {
	reverse_double(cnorm[i]);
    }
}

#endif

struct urcinfo {
    int nz;
    int nreg;
    int model;
    int Tmin;
    int pos;
};

/*
   niv = # of integrated variables
   itv = appropriate ur_code for nc, c, ct, ctt models
   T = sample size (0 for asymptotic)
   tau = test statistic
   err = location to receive error code
*/

static double urcval (int niv, int itv, int T, double tau,
		      int *err)
{
    FILE *fp;
    gchar *datapath = NULL;
    double beta[BIGLEN];
    double wght[URCLEN+1];
    double prob[URCLEN+1];
    double cnorm[URCLEN+1];
    struct urcinfo uis[] = {
	{0,  1, 2, 20,      0}, /* dfnc: table 1 */
	{0,  2, 2, 20,   7072}, /* dfc */
	{0,  3, 3, 20,  14144}, /* dfct */
	{0,  4, 3, 20,  22984}, /* dfctt */
	{1,  2, 2, 20,  31824}, /* conc: table 2 */
	{1,  3, 2, 20,  38896}, /* coc */
	{1,  4, 3, 25,  45968}, /* coct */
	{1,  5, 3, 20,  54808}, /* coctt */
	{2,  3, 2, 25,  63648}, /* conc: table 3 */
	{2,  4, 2, 20,  70720}, /* coc */
	{2,  5, 2, 20,  77792}, /* coct */
	{2,  6, 3, 20,  84864}, /* coctt */
	{3,  4, 3, 20,  93704}, /* conc: table 4 */
	{3,  5, 2, 25, 102544}, /* coc */
	{3,  6, 3, 20, 109616}, /* coct */
	{3,  7, 2, 30, 118456}, /* coctt */
	{4,  5, 2, 25, 125528}, /* conc: table 5 */
	{4,  6, 3, 20, 132600}, /* coc */
	{4,  7, 3, 20, 141440}, /* coct */
	{4,  8, 3, 20, 150280}, /* coctt */
	{5,  6, 2, 30, 159120}, /* conc: table 6 */
	{5,  7, 2, 30, 166192}, /* coc */
	{5,  8, 2, 30, 173264}, /* coct */
	{5,  9, 3, 25, 180336}, /* coctt */
	{6,  7, 3, 25, 189176}, /* conc: table 7 */
	{6,  8, 3, 25, 198016}, /* coc */
	{6,  9, 3, 30, 206856}, /* coct */
	{6, 10, 3, 30, 215696}, /* coctt */
	{7,  8, 2, 40, 224536}, /* conc: table 8 */
	{7,  9, 2, 35, 231608}, /* coc */
	{7, 10, 2, 40, 238680}, /* coct */
	{7, 11, 2, 40, 245752}, /* coctt */
	{0,  0, 0,  0, 252824}, /* prob */
	{0,  0, 0, -1, 254592}  /* cnorm */
    };
    struct urcinfo *ui;
    size_t nr1, nr2;
    int i, nbeta;
    double pval = NADBL;

    /* Check that parameters are valid */
    if (niv < 1 || niv > NIVMAX) {
	*err = E_DATA;
	return pval;
    }

    if (itv < 1 || itv > 4) {
	/* these limits correspond to UR_NO_CONST and UR_QUAD_TREND
	   in lib/src/adf_kpss.c */
	*err = E_DATA;
	return pval;
    }

    /* Open data file */
    datapath = g_strdup_printf("%sdata%curcdata.bin",
			       gretl_plugin_path(), SLASH);
    fp = gretl_fopen(datapath, "rb");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", datapath);
	*err = E_FOPEN;
	g_free(datapath);
	return pval;
    }

    g_free(datapath);

    /* skip to appropriate location in data file */
    i = (niv-1) * 4 + (itv - 1);
    ui = &uis[i];
    fseek(fp, ui->pos, SEEK_SET);

    /* the number of coefficients in the critical values
       equations */
    nbeta = ui->model == 2 ? 3 : 4;

#if URDEBUG
    fprintf(fdb, "nz=%d, nreg=%d, model=%d, Tmin=%d, offset=%d\n",
	    ui->nz, ui->nreg, ui->model, ui->Tmin, ui->pos);
    fflush(fdb);
#endif

    /* these arrays are padded by one for fortran */
    wght[0] = prob[0] = cnorm[0] = 0.0;

    /* read coefficients and weights */
    nr1 = fread(beta, sizeof(double), nbeta * URCLEN, fp);
    nr2 = fread(wght + 1, sizeof(double), URCLEN, fp);
    if (nr1 != nbeta * URCLEN || nr2 != URCLEN) {
	fprintf(stderr, "error reading urcdata\n");
	*err = E_DATA;
    }

    if (!*err) {
	/* read from embedded "probs.tab" data */
	fseek(fp, uis[32].pos, SEEK_SET);
	nr1 = fread(prob + 1, sizeof(double), URCLEN, fp);
	nr2 = fread(cnorm + 1, sizeof(double), URCLEN, fp);
	if (nr1 != URCLEN || nr2 != URCLEN) {
	    fprintf(stderr, "error reading urcdata\n");
	    *err = E_DATA;
	}
    }

    fclose(fp);

#if G_BYTE_ORDER == G_BIG_ENDIAN
    if (!*err) {
	urc_swap_endianness(beta, nbeta, wght, prob, cnorm);
    }
#endif

    if (!*err && T > 0 && T < ui->Tmin) {
	/* error, or warning? */
	fprintf(stderr, "Warning, too few observations!\n");
	/* *err = E_TOOFEW; */
    }

    if (!*err) {
	pval = fpval(beta, nbeta, wght, prob, cnorm,
		     tau, T, err);
    }

    return pval;
}

/*
   tau = test statistic
   T = sample size (or 0 for asymptotic)
   niv = # of integrated variables
   itv = 1, 2, 3, or 4 for nc, c, ct, ctt models.

   returns: the computed P-value
*/

double mackinnon_pvalue (double tau, int T, int niv, int itv)
{
    double pval = NADBL;
    int err = 0;

#if URDEBUG
    fdb = fopen("debug.txt", "w");
    if (fdb != NULL) {
	fprintf(fdb, "mackinnon_pvalue: tau=%g, T=%d, niv=%d, itv=%d\n",
		tau, T, niv, itv);
	fflush(fdb);
    }
#endif

    pval = urcval(niv, itv, T, tau, &err);

#if URDEBUG
    if (fdb != NULL) {
	fclose(fdb);
    }
#endif

    return pval;
}

/* Extra: code to obtain approximate finite-sample p-values for
   DF-GLS tests a la Elliott-Rothenberg-Stock, using Sephton's
   response surfaces for the test-down cases, otherwise Cottrell
   and Komashko surfaces.
*/

#define N_ALPHA_C 25
#define N_ALPHA_T 26

static const double alpha_c[N_ALPHA_C] = {
    0.001, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
    0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.15, 0.2, 0.3,
    0.4, 0.5, 0.6, 0.7, 0.9, 0.99
};

static const double alpha_t[N_ALPHA_T] = {
    0.001, 0.0025, 0.005, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06,
    0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.15, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.9, 0.99
};

static const char *bnames[] = {
    "dfgls-beta-c.bin", "dfgls-beta-t.bin"
};

static const char *snames[] = {
    "npc.bin", "npt.bin", "pqc.bin", "pqt.bin"
};

static gretl_matrix *get_data_matrix (const char *base, int *err)
{
    gretl_matrix *m;
    gchar *fname;

    fname = g_strdup_printf("%sdata%c%s", gretl_plugin_path(),
			    SLASH, base);
    m = gretl_matrix_read_from_file(fname, 0, err);
    if (*err) {
	fprintf(stderr, "Couldn't open %s\n", fname);
    }
    g_free(fname);

    return m;
}

#define PDEBUG 0

double dfgls_pvalue (double tau, int T, int trend,
		     int kmax, int PQ, int *err)
{
    const char *fbase;
    gretl_matrix_block *B;
    gretl_matrix *alpha = NULL;
    gretl_matrix *beta = NULL;
    gretl_matrix *C, *V, *g;
    gretl_matrix *r, *y, *X;
    const double *a, *b;
    double d, dmin = 1.0e6;
    double ci, s2, pval = NADBL;
    int i1, i2, npoints, np2;
    int sephton = (kmax > 0);
    int nreg, ncoeff, nalpha;
    int i, j, k, imin = 0;

    fbase = sephton ? snames[trend + 2*PQ] : bnames[trend];
    beta = get_data_matrix(fbase, err);

    if (!*err && sephton) {
	alpha = get_data_matrix("s_alpha.bin", err);
	if (!*err) {
	    a = alpha->val;
	    nalpha = alpha->rows;
	    npoints = 13;
	}
    } else if (!*err) {
	a = trend ? alpha_t : alpha_c;
	nalpha = trend ? N_ALPHA_T : N_ALPHA_C;
	npoints = 5;
    }

    if (*err) {
	gretl_matrix_free(beta);
	return pval;
    }

    ncoeff = beta->rows;
    b = beta->val;
    np2 = npoints / 2;

#if PDEBUG
    fprintf(stderr, "dfgls: tau %.8g, trend %d, kmax %d, PQ %d\n",
	    tau, trend, kmax, PQ);
    fprintf(stderr, "dfgls: ncoeff %d, nalpha %d\n", ncoeff, nalpha);
#endif

    /* max # of second-stage regressors */
    nreg = 4;

    /* allocate all storage */
    B = gretl_matrix_block_new(&C, nalpha, 1,
			       &y, npoints, 1,
			       &X, npoints, nreg,
			       &g, nreg, 1,
			       &V, nreg, nreg,
			       &r, ncoeff, 1, NULL);
    if (B == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (T == 0) {
	/* asymptotic value */
	r->val[0] = 1.0;
	for (j=1; j<ncoeff; j++) {
	    r->val[j] = 0;
	}
    } else {
	/* finite sample */
	int npow = ncoeff - 2 - sephton;
	double Tr = 1.0 / T;

	r->val[0] = 1.0;
	r->val[1] = Tr;
	for (j=0; j<npow; j++) {
	    r->val[j+2] = pow(Tr, j+2);
	}
	if (sephton) {
	    r->val[ncoeff-1] = kmax * Tr;
	}
    }

    /* compute all @nalpha critical values and determine
       which is closest to @tau.
    */
    k = 0;
    for (i=0; i<nalpha; i++) {
	ci = 0;
	for (j=0; j<ncoeff; j++) {
	    ci += b[k+j] * r->val[j];
	}
	gretl_matrix_set(C, i, 0, ci);
	if (i == 0) {
	    if (tau < 1.25 * C->val[0]) {
		/* extrapolation unlikely to work well for big
		   negative tau: employ conservative fudge
		*/
		pval = 1.0e-5;
		goto bailout;
	    }
	}
	d = fabs(ci - tau);
	if (d < dmin) {
	    dmin = d;
	    imin = i;
	}
	k += ncoeff;
    }

#if PDEBUG
    if (sephton) {
	/* these value agree with Octave on Sephton's programs */
	fprintf(stderr, "sephton: compute critical values for T = %d\n", T);
	fprintf(stderr, "cv 0.01: %g\n", C->val[12]);
	fprintf(stderr, "cv 0.05: %g\n", C->val[20]);
	fprintf(stderr, "cv 0.10: %g\n", C->val[30]);
	fprintf(stderr, "sephton: dmin = %g, imin = %d\n", dmin, imin);
    }
#endif

    /* select starting point @i1 for critvals range */
    i1 = imin - np2; i2 = imin + np2;
    if (i1 < 0) {
	i1 = 0;
    } else if (i2 >= nalpha) {
	i2 = nalpha - 1;
	i1 = i2 - npoints + 1;
    }

    /* fill out y */
    for (i=0; i<npoints; i++) {
	y->val[i] = normal_cdf_inverse(a[i1 + i]);
    }

#if PDEBUG
    fprintf(stderr, " i1 = %d, i2 = %d\n", i1, i2);
    for (i=0; i<npoints; i++) {
	fprintf(stderr, " %d: a=%g -> y=%g\n", i, a[i1 + i], y->val[i]);
    }
#endif

    /* fill out X */
    for (j=0; j<nreg; j++) {
	for (i=0; i<npoints; i++) {
	    if (j == 0) {
		gretl_matrix_set(X, i, j, 1.0);
	    } else {
		ci = C->val[i1 + i];
		if (j == 1) {
		    gretl_matrix_set(X, i, j, ci);
		} else {
		    gretl_matrix_set(X, i, j, pow(ci, j));
		}
	    }
	}
    }

    /* we use SVD here just to avoid spamming stderr with
       reports of near-perfect collinearity
    */
    *err = gretl_matrix_SVD_ols(y, X, g, V, NULL, &s2);

    if (!*err) {
	double t3 = g->val[3] / sqrt(gretl_matrix_get(V, 3, 3));

	if (fabs(t3) < 2) {
	    /* drop the cubic term */
	    X->cols = 3;
	    g->rows = 3;
	    *err = gretl_matrix_SVD_ols(y, X, g, NULL, NULL, NULL);
	}
	if (!*err) {
	    ci = g->val[0] + g->val[1] * tau;
	    for (i=2; i<g->rows; i++) {
		ci += pow(tau, i) * g->val[i];
	    }
	    pval = normal_cdf(ci);
	}
    }

    if (!*err && (pval < 0 || pval > 1)) {
	/* allow for a little inaccuracy at the extremes */
	pval = pval < 0 ? 0 : 1;
    }

 bailout:

    gretl_matrix_free(beta);
    gretl_matrix_free(alpha);
    gretl_matrix_block_destroy(B);

    return pval;
}
