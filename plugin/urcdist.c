/*
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

/* urcrouts.f -- translated by f2c (version 20030306) */

struct {
    doublereal probs[221], cnorm[221], beta[884]	/* was [4][221] */, 
	    wght[221];
    integer nz, nreg, model, minsize;
} saveall_;

#define saveall_1 saveall_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__20 = 20;
static integer c__4 = 4;
static integer c__0 = 0;

/* Copyright (c) James G. MacKinnon, 1996 (corrected 2003-5-5) */

/* urcrouts.f: This is a set of subroutines to estimate critical values */
/* and P values for unit root and cointegration tests. It is written in */
/* Fortran 77. Simply call urcval, specifying the first seven arguments. */
/* The result comes back in the last argument. */

/* A standalone program called urcdist.f is also available. */

/* These routines and the associated data files may be used freely for */
/* non-commercial purposes, provided that proper attribution is made. */
/* Please cite the paper */

/*  James G. MacKinnon, "Numerical distribution functions for unit root */
/*  and cointegration tests," Journal of Applied Econometrics, 11, */
/*  1996, 601-618. */

/* The routines and data files may not be incorporated into any book */
/* or computer program without the express, written consent of the author. */

/* The routines must have access to the files probs.tab and urc-#.tab */
/* for # = 1, 2, 3 ..., 12. As currently written, these files must be */
/* in the current directory or in the directory /usr/local/urcdist. Make */
/* sure that these files are in the proper format for your computer (i.e. */
/* that lines are terminated by CR/LF for DOS, Windows, and OS/2 systems */
/* and by LF alone for Unix systems). */

/* Subroutine */ int urcval_(integer *niv, integer *itt, integer *itv, 
	integer *nobs, doublereal *arg, doublereal *val, const char *path)
{
    /* Format strings */
    static char fmt_301[] = "(1a1)";
    static char fmt_302[] = "(8x,3i3,i5)";
    static char fmt_303[] = "(5d15.8)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* return value */
    int urc_ret = URC_OK;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_rsle(cilist *), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, np, nx;
    static doublereal pval, stat;
    static integer nvar, junk;
    extern /* Subroutine */ int fcrit_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *), fpval_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer iskip;
    static doublereal precrt;
    static char usrfile[FILENAME_MAX];

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 11, 0, fmt_301, 0 };
    static cilist io___11 = { 0, 11, 0, fmt_301, 0 };
    static cilist io___12 = { 0, 11, 0, fmt_302, 0 };
    static cilist io___14 = { 0, 11, 0, fmt_303, 0 };
    static cilist io___16 = { 0, 12, 0, 0, 0 };


/* niv = # of integrated variables */
/* itt = 1 or 2 for tau or z test */
/* itv = 1, 2, 3, 4 for nc, c, ct, ctt */
/* nobs = sample size (0 for asymptotic) */
/* arg = level of test (between .0001 and .9999) if nc = 1 */
/*     = test statistic if nc = 2 */
/* val = critical value if nc = 1 (returned by routine) */
/*     = P value if nc = 2 (returned by routine) */
/* ** Note that arg and val are double precision (real*8) variables ** */


/* common block saveall should not appear in any other routine */

/* Check that parameters are valid. */

    if (*niv < 1 || *niv > 12) {
	return URC_BAD_PARAM;
    }
    if (*itt < 1 || *itt > 2) {
	return URC_BAD_PARAM;
    }
    if (*itv < 1 || *itv > 4) {
	return URC_BAD_PARAM;
    }

    /* Open data files */

    sprintf(usrfile, "%sdata%cprobs.tab", path, SLASH);
    o__1.oerr = 1;
    o__1.ounit = 12;
    o__1.ofnmlen = strlen(usrfile);
    o__1.ofnm = usrfile;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	return URC_NOT_FOUND;
    }

    sprintf(usrfile, "%sdata%curc-%d.tab", path, SLASH, (int) *niv);
    o__1.oerr = 1;
    o__1.ounit = 11;
    o__1.ofnmlen = strlen(usrfile);
    o__1.ofnm = usrfile;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	cl__1.cerr = 0;
	cl__1.cunit = 12;
	cl__1.csta = 0;
	f_clos(&cl__1);
	return URC_NOT_FOUND;
    }

/* read data from unit 11. */

/* skip copyright line */

    s_rsfe(&io___7);
    do_fio(&c__1, (char *)&junk, (ftnlen)sizeof(integer));
    e_rsfe();
    iskip = 0;

/* skip groups of 222 lines as necessary. */

    if (*itt != 1) {
	iskip = 888;
    }
    iskip += (*itv - 1) * 222;
    if (iskip > 0) {
	i__1 = iskip;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_rsfe(&io___11);
	    do_fio(&c__1, (char *)&junk, (ftnlen)sizeof(integer));
	    e_rsfe();
	}
    }

    s_rsfe(&io___12);
    do_fio(&c__1, (char *)&saveall_1.nz, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&saveall_1.nreg, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&saveall_1.model, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&saveall_1.minsize, (ftnlen)sizeof(integer));
    e_rsfe();
    if (saveall_1.model == 2 || saveall_1.model == 4) {
	nvar = 3;
    } else {
	nvar = 4;
    }

    for (i__ = 1; i__ <= 221; ++i__) {
	s_rsfe(&io___14);
	i__1 = nvar;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&saveall_1.beta[j + (i__ << 2) - 5], (
		    ftnlen)sizeof(doublereal));
	}
	do_fio(&c__1, (char *)&saveall_1.wght[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsfe();
    }
    for (i__ = 1; i__ <= 221; ++i__) {
	s_rsle(&io___16);
	do_lio(&c__5, &c__1, (char *)&saveall_1.probs[i__ - 1], (ftnlen)
		sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&saveall_1.cnorm[i__ - 1], (ftnlen)
		sizeof(doublereal));
	e_rsle();
    }

    stat = *arg;

    if (*nobs > 0 && *nobs < saveall_1.minsize) {
	urc_ret = URC_SMALL_SAMPLE;
    }

    np = 9;
    precrt = 2.;
    fpval_(saveall_1.beta, saveall_1.cnorm, saveall_1.wght, 
	   saveall_1.probs, &pval, &stat, &precrt, nobs, &
	   saveall_1.model, &saveall_1.nreg, &np, &nx);
    *val = pval;

    cl__1.cerr = 0;
    cl__1.cunit = 11;
    cl__1.csta = 0;
    f_clos(&cl__1);
    cl__1.cerr = 0;
    cl__1.cunit = 12;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return urc_ret;
} /* urcval_ */

/* Subroutine */ int fcrit_(doublereal *probs, doublereal *cnorm, doublereal *
	beta, doublereal *wght, doublereal *cval, doublereal *size, 
	doublereal *precrt, integer *nobs, integer *model, integer *nreg, 
	integer *np, integer *nx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, ic, jc;
    static doublereal sd4;
    static integer np1;
    static doublereal bot;
    static integer nph;
    extern /* Subroutine */ int gls_(doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *);
    static doublereal top, ssr, diff;
    extern /* Subroutine */ int eval_(doublereal *, doublereal *, integer *, 
	    integer *, integer *);
    static integer imin;
    static doublereal fits[20], xmat[80]	/* was [20][4] */, xomx[16]	
	    /* was [4][4] */, ssrt, gamma[4], diffm, omega[400]	/* was [20][
	    20] */, resid[20], anorm, crits[221], yvect[20];
    static integer nptop;
    static doublereal ttest;
    extern /* Subroutine */ int innorz_(doublereal *, doublereal *);


/* Copyright (c) James G. MacKinnon, 1995 */
/* Routine to find a critical value for any specified test size. */
/* Uses GLS to estimate approximating regression. */

    /* Parameter adjustments */
    --wght;
    beta -= 5;
    --cnorm;
    --probs;

    /* Function Body */
    diffm = 1e3;
    imin = 0;
    for (i__ = 1; i__ <= 221; ++i__) {
	diff = (d__1 = *size - probs[i__], abs(d__1));
	if (diff < diffm) {
	    diffm = diff;
	    imin = i__;
	    if (diffm < 1e-6) {
		goto L100;
	    }
	}
    }
L100:

    nph = *np / 2;
    nptop = 221 - nph;
    if (imin > nph && imin < nptop) {

/* imin is not too close to the end. Use np points around stat. */

	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ic = imin - nph - 1 + i__;
	    eval_(&beta[(ic << 2) + 1], &crits[ic - 1], model, nreg, nobs);
	    yvect[i__ - 1] = crits[ic - 1];
	    xmat[i__ - 1] = 1.;
	    xmat[i__ + 19] = cnorm[ic];
	    xmat[i__ + 39] = xmat[i__ + 19] * cnorm[ic];
	    xmat[i__ + 59] = xmat[i__ + 39] * cnorm[ic];
	}

/* form omega matrix */

	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *np;
	    for (j = i__; j <= i__2; ++j) {
		ic = imin - nph - 1 + i__;
		jc = imin - nph - 1 + j;
		top = probs[ic] * (1. - probs[jc]);
		bot = probs[jc] * (1. - probs[ic]);
		omega[i__ + j * 20 - 21] = wght[ic] * wght[jc] * sqrt(top / 
			bot);
	    }
	}
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *np;
	    for (j = i__; j <= i__2; ++j) {
		omega[j + i__ * 20 - 21] = omega[i__ + j * 20 - 21];
	    }
	}

	*nx = 4;
	gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, np, 
		nx, &c__20, &c__4, &c__0);

/* check to see if gamma(4) is needed */

	sd4 = sqrt(ssrt / (*np - *nx) * xomx[15]);
	ttest = abs(gamma[3]) / sd4;
	if (ttest > *precrt) {
	    innorz_(size, &anorm);
/* Computing 2nd power */
	    d__1 = anorm;
/* Computing 3rd power */
	    d__2 = anorm;
	    *cval = gamma[0] + gamma[1] * anorm + gamma[2] * (d__1 * d__1) + 
		    gamma[3] * (d__2 * (d__2 * d__2));
	    return 0;
	} else {
	    *nx = 3;
	    gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, 
		    np, nx, &c__20, &c__4, &c__1);
	    innorz_(size, &anorm);
/* Computing 2nd power */
	    d__1 = anorm;
	    *cval = gamma[0] + gamma[1] * anorm + gamma[2] * (d__1 * d__1);
	    return 0;
	}

/* imin is close to one of the ends. Use points from imin +/- nph to end. */

    } else {
	if (imin < *np) {
	    np1 = imin + nph;
	    if (np1 < 5) {
		np1 = 5;
	    }
	    i__1 = np1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		eval_(&beta[(i__ << 2) + 1], &crits[i__ - 1], model, nreg, 
			nobs);
		yvect[i__ - 1] = crits[i__ - 1];
		xmat[i__ - 1] = 1.;
		xmat[i__ + 19] = cnorm[i__];
		xmat[i__ + 39] = xmat[i__ + 19] * cnorm[i__];
		xmat[i__ + 59] = xmat[i__ + 39] * cnorm[i__];
	    }
	} else {
	    np1 = 222 - imin + nph;
	    if (np1 < 5) {
		np1 = 5;
	    }
	    i__1 = np1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		eval_(&beta[(222 - i__ << 2) + 1], &crits[222 - i__ - 1], 
			model, nreg, nobs);
		ic = 222 - i__;
		yvect[i__ - 1] = crits[ic - 1];
		xmat[i__ - 1] = 1.;
		xmat[i__ + 19] = cnorm[ic];
		xmat[i__ + 39] = xmat[i__ + 19] * cnorm[ic];
		xmat[i__ + 59] = xmat[i__ + 39] * cnorm[ic];
	    }
	}

/* form omega matrix */

	i__1 = np1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = np1;
	    for (j = i__; j <= i__2; ++j) {
		if (imin < *np) {
		    top = probs[i__] * (1. - probs[j]);
		    bot = probs[j] * (1. - probs[i__]);
		    omega[i__ + j * 20 - 21] = wght[i__] * wght[j] * sqrt(top 
			    / bot);
		} else {

/* This is to avoid numerical singularities at the upper end */

		    omega[i__ + j * 20 - 21] = 0.;
		    if (i__ == j) {
			omega[i__ + i__ * 20 - 21] = 1.;
		    }
		}
	    }
	}
	i__1 = np1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = np1;
	    for (j = i__; j <= i__2; ++j) {
		omega[j + i__ * 20 - 21] = omega[i__ + j * 20 - 21];
	    }
	}

	*nx = 4;
	gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, &np1, 
		nx, &c__20, &c__4, &c__0);

/* check to see if gamma(4) is needed */

	sd4 = sqrt(ssrt / (np1 - *nx) * xomx[15]);
	ttest = (d__1 = gamma[3] / sd4, abs(d__1));
	if (ttest > *precrt) {
	    innorz_(size, &anorm);
/* Computing 2nd power */
	    d__1 = anorm;
/* Computing 3rd power */
	    d__2 = anorm;
	    *cval = gamma[0] + gamma[1] * anorm + gamma[2] * (d__1 * d__1) + 
		    gamma[3] * (d__2 * (d__2 * d__2));
	    return 0;
	} else {
	    *nx = 3;
	    gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, &
		    np1, nx, &c__20, &c__4, &c__1);
	    innorz_(size, &anorm);
/* Computing 2nd power */
	    d__1 = anorm;
	    *cval = gamma[0] + gamma[1] * anorm + gamma[2] * (d__1 * d__1);
	    return 0;
	}

    }
    return 0;
} /* fcrit_ */

/* Subroutine */ int fpval_(doublereal *beta, doublereal *cnorm, doublereal *
	wght, doublereal *probs, doublereal *pval, doublereal *stat, 
	doublereal *precrt, integer *nobs, integer *model, integer *nreg, 
	integer *np, integer *nx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, ic, jc;
    static doublereal sd4;
    static integer np1;
    static doublereal bot;
    static integer nph;
    extern /* Subroutine */ int gls_(doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *);
    static doublereal top, ssr, diff;
    extern /* Subroutine */ int eval_(doublereal *, doublereal *, integer *, 
	    integer *, integer *);
    static integer imin;
    static doublereal fits[20], xmat[80]	/* was [20][4] */, xomx[16]	
	    /* was [4][4] */, ssrt, gamma[4], diffm, omega[400]	/* was [20][
	    20] */, resid[20], crfit;
    extern /* Subroutine */ int ddnor_(doublereal *, doublereal *);
    static doublereal crits[221], yvect[20];
    static integer nptop;
    static doublereal ttest;


/* Copyright (c) James G. MacKinnon, 1995 */
/* Routine to find P value for any specified test statistic. */


/* first, compute all the estimated critical values */

    /* Parameter adjustments */
    --probs;
    --wght;
    --cnorm;
    beta -= 5;

    /* Function Body */
    for (i__ = 1; i__ <= 221; ++i__) {
	eval_(&beta[(i__ << 2) + 1], &crits[i__ - 1], model, nreg, nobs);
    }

/* find critical value closest to test statistic */

    diffm = 1e3;
    imin = 0;
    for (i__ = 1; i__ <= 221; ++i__) {
	diff = (d__1 = *stat - crits[i__ - 1], abs(d__1));
	if (diff < diffm) {
	    diffm = diff;
	    imin = i__;
	}
    }

    nph = *np / 2;
    nptop = 221 - nph;
    if (imin > nph && imin < nptop) {

/* imin is not too close to the end. Use np points around stat. */

	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ic = imin - nph - 1 + i__;
	    yvect[i__ - 1] = cnorm[ic];
	    xmat[i__ - 1] = 1.;
	    xmat[i__ + 19] = crits[ic - 1];
	    xmat[i__ + 39] = xmat[i__ + 19] * crits[ic - 1];
	    xmat[i__ + 59] = xmat[i__ + 39] * crits[ic - 1];
	}

/* form omega matrix */

	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *np;
	    for (j = i__; j <= i__2; ++j) {
		ic = imin - nph - 1 + i__;
		jc = imin - nph - 1 + j;
		top = probs[ic] * (1. - probs[jc]);
		bot = probs[jc] * (1. - probs[ic]);
		omega[i__ + j * 20 - 21] = wght[ic] * wght[jc] * sqrt(top / 
			bot);
	    }
	}
	i__1 = *np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *np;
	    for (j = i__; j <= i__2; ++j) {
		omega[j + i__ * 20 - 21] = omega[i__ + j * 20 - 21];
	    }
	}

	*nx = 4;
	gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, np, 
		nx, &c__20, &c__4, &c__0);

/* check to see if gamma(4) is needed */

	sd4 = sqrt(ssrt / (*np - *nx) * xomx[15]);
	ttest = abs(gamma[3]) / sd4;
	if (ttest > *precrt) {
/* Computing 2nd power */
	    d__1 = *stat;
/* Computing 3rd power */
	    d__2 = *stat;
	    crfit = gamma[0] + gamma[1] * *stat + gamma[2] * (d__1 * d__1) + 
		    gamma[3] * (d__2 * (d__2 * d__2));
	    ddnor_(&crfit, pval);
	    return 0;
	} else {
	    *nx = 3;
	    gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, 
		    np, nx, &c__20, &c__4, &c__1);
/* Computing 2nd power */
	    d__1 = *stat;
	    crfit = gamma[0] + gamma[1] * *stat + gamma[2] * (d__1 * d__1);
	    ddnor_(&crfit, pval);
	    return 0;
	}
    } else {

/* imin is close to one of the ends. Use points from imin +/- nph to end. */

	if (imin < *np) {
	    np1 = imin + nph;
	    if (np1 < 5) {
		np1 = 5;
	    }
	    i__1 = np1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		yvect[i__ - 1] = cnorm[i__];
		xmat[i__ - 1] = 1.;
		xmat[i__ + 19] = crits[i__ - 1];
		xmat[i__ + 39] = xmat[i__ + 19] * crits[i__ - 1];
		xmat[i__ + 59] = xmat[i__ + 39] * crits[i__ - 1];
	    }
	} else {
	    np1 = 222 - imin + nph;
	    if (np1 < 5) {
		np1 = 5;
	    }
	    i__1 = np1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ic = 222 - i__;
		yvect[i__ - 1] = cnorm[ic];
		xmat[i__ - 1] = 1.;
		xmat[i__ + 19] = crits[ic - 1];
		xmat[i__ + 39] = xmat[i__ + 19] * crits[ic - 1];
		xmat[i__ + 59] = xmat[i__ + 39] * crits[ic - 1];
	    }
	}

/* form omega matrix */

	i__1 = np1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = np1;
	    for (j = i__; j <= i__2; ++j) {
		if (imin < *np) {
		    top = probs[i__] * (1. - probs[j]);
		    bot = probs[j] * (1. - probs[i__]);
		    omega[i__ + j * 20 - 21] = wght[i__] * wght[j] * sqrt(top 
			    / bot);
		} else {

/* This is to avoid numerical singularities at the upper end */

		    omega[i__ + j * 20 - 21] = 0.;
		    if (i__ == j) {
			omega[i__ + i__ * 20 - 21] = 1.;
		    }
		}
	    }
	}
	i__1 = np1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = np1;
	    for (j = i__; j <= i__2; ++j) {
		omega[j + i__ * 20 - 21] = omega[i__ + j * 20 - 21];
	    }
	}

	*nx = 4;
	gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, &np1, 
		nx, &c__20, &c__4, &c__0);

/* check to see if gamma(4) is needed */

	sd4 = sqrt(ssrt / (np1 - *nx) * xomx[15]);
	ttest = abs(gamma[3]) / sd4;
	if (ttest > *precrt) {
/* Computing 2nd power */
	    d__1 = *stat;
/* Computing 3rd power */
	    d__2 = *stat;
	    crfit = gamma[0] + gamma[1] * *stat + gamma[2] * (d__1 * d__1) + 
		    gamma[3] * (d__2 * (d__2 * d__2));
	    ddnor_(&crfit, pval);
	} else {
	    *nx = 3;
	    gls_(xmat, yvect, omega, gamma, xomx, fits, resid, &ssr, &ssrt, &
		    np1, nx, &c__20, &c__4, &c__1);
/* Computing 2nd power */
	    d__1 = *stat;
	    crfit = gamma[0] + gamma[1] * *stat + gamma[2] * (d__1 * d__1);
	    ddnor_(&crfit, pval);
	}

/* check that nothing crazy has happened at the ends */

	if (imin == 1 && *pval > probs[1]) {
	    *pval = probs[1];
	}
	if (imin == 221 && *pval < probs[221]) {
	    *pval = probs[221];
	}
	return 0;
    }
    return 0;
} /* fpval_ */

/* Subroutine */ int eval_(doublereal *beta, doublereal *cval, integer *model,
	 integer *nreg, integer *nobs)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal onobs;

    /* Fortran I/O blocks */
    static cilist io___83 = { 0, 6, 0, 0, 0 };



/* Copyright (c) James G. MacKinnon, 1995 */
/* Routine to evaluate response surface for specified betas and sample size. */

    /* Parameter adjustments */
    --beta;

    /* Function Body */
    if (*nobs == 0) {
	*cval = beta[1];
	return 0;
    }
    if (*model == 2) {
	onobs = 1. / *nobs;
/* Computing 2nd power */
	d__1 = onobs;
	*cval = beta[1] + beta[2] * onobs + beta[3] * (d__1 * d__1);
	return 0;
    }
    if (*model == 3) {
	onobs = 1. / *nobs;
/* Computing 2nd power */
	d__1 = onobs;
/* Computing 3rd power */
	d__2 = onobs;
	*cval = beta[1] + beta[2] * onobs + beta[3] * (d__1 * d__1) + beta[4] 
		* (d__2 * (d__2 * d__2));
	return 0;
    }
    if (*model == 4) {
	onobs = 1. / (*nobs - *nreg);
/* Computing 2nd power */
	d__1 = onobs;
	*cval = beta[1] + beta[2] * onobs + beta[3] * (d__1 * d__1);
	return 0;
    }
    if (*model == 5) {
	onobs = 1. / (*nobs - *nreg);
/* Computing 2nd power */
	d__1 = onobs;
/* Computing 3rd power */
	d__2 = onobs;
	*cval = beta[1] + beta[2] * onobs + beta[3] * (d__1 * d__1) + beta[4] 
		* (d__2 * (d__2 * d__2));
	return 0;
    }
    s_wsle(&io___83);
    do_lio(&c__9, &c__1, "*** Warning! Error in input file. ***", (ftnlen)37);
    e_wsle();
    return 0;
} /* eval_ */

/* Subroutine */ int gls_(doublereal *xmat, doublereal *yvect, doublereal *
	omega, doublereal *beta, doublereal *xomx, doublereal *fits, 
	doublereal *resid, doublereal *ssr, doublereal *ssrt, integer *nobs, 
	integer *nvar, integer *nomax, integer *nvmax, integer *ivrt)
{
    /* System generated locals */
    integer xmat_dim1, xmat_offset, omega_dim1, omega_offset, xomx_dim1, 
	    xomx_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, kxx;
    static doublereal xomy[50];
    extern /* Subroutine */ int cholx_(doublereal *, integer *, integer *, 
	    integer *);


/* Copyright (c) James G. MacKinnon, 1995 */
/* Subroutine to do GLS estimation the obvious way */
/* Use only when sample size is small (nobs <= 50) */
/* 1995-1-3 */


/* xomx is covariance matrix of parameter estimates if omega is truly known */
/* First, invert omega matrix if ivrt=0. Original one gets replaced. */

    /* Parameter adjustments */
    --resid;
    --fits;
    omega_dim1 = *nomax;
    omega_offset = 1 + omega_dim1;
    omega -= omega_offset;
    --yvect;
    xomx_dim1 = *nvmax;
    xomx_offset = 1 + xomx_dim1;
    xomx -= xomx_offset;
    --beta;
    xmat_dim1 = *nomax;
    xmat_offset = 1 + xmat_dim1;
    xmat -= xmat_offset;

    /* Function Body */
    if (*ivrt == 0) {
	cholx_(&omega[omega_offset], nomax, nobs, &kxx);
    }

/* form xomx matrix and xomy vector */

    i__1 = *nvar;
    for (j = 1; j <= i__1; ++j) {
	xomy[j - 1] = 0.;
	i__2 = *nvar;
	for (l = j; l <= i__2; ++l) {
	    xomx[j + l * xomx_dim1] = 0.;
	}
    }

    i__1 = *nobs;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nobs;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *nvar;
	    for (j = 1; j <= i__3; ++j) {
		xomy[j - 1] += xmat[i__ + j * xmat_dim1] * omega[k + i__ * 
			omega_dim1] * yvect[k];
		i__4 = *nvar;
		for (l = j; l <= i__4; ++l) {
		    xomx[j + l * xomx_dim1] += xmat[i__ + j * xmat_dim1] * 
			    omega[k + i__ * omega_dim1] * xmat[k + l * 
			    xmat_dim1];
/* L24: */
		}
	    }
/* L21: */
	}
    }

    i__2 = *nvar;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *nvar;
	for (l = j; l <= i__1; ++l) {
	    xomx[l + j * xomx_dim1] = xomx[j + l * xomx_dim1];
	}
    }

/* invert xomx matrix */

    cholx_(&xomx[xomx_offset], nvmax, nvar, &kxx);

/*  now form estimates of beta. */

    i__2 = *nvar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	beta[i__] = 0.;
	i__1 = *nvar;
	for (j = 1; j <= i__1; ++j) {
	    beta[i__] += xomx[i__ + j * xomx_dim1] * xomy[j - 1];
/* L5: */
	}
    }

/* find ssr, fitted values, and residuals */

    *ssr = 0.;
    i__1 = *nobs;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fits[i__] = 0.;
	i__2 = *nvar;
	for (j = 1; j <= i__2; ++j) {
	    fits[i__] += xmat[i__ + j * xmat_dim1] * beta[j];
	}
	resid[i__] = yvect[i__] - fits[i__];
/* Computing 2nd power */
	d__1 = resid[i__];
	*ssr += d__1 * d__1;
    }

/* find ssr from transformed regression */

    *ssrt = 0.;
    i__1 = *nobs;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nobs;
	for (k = 1; k <= i__2; ++k) {
	    *ssrt += resid[i__] * omega[k + i__ * omega_dim1] * resid[k];
	}
    }

    return 0;
} /* gls_ */

/* Subroutine */ int cholx_(doublereal *amat, integer *m, integer *n, integer 
	*kxx)
{
    /* System generated locals */
    integer amat_dim1, amat_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer kl;
    static doublereal ooa;


/* Copyright (c) James G. MacKinnon, 1993 */
/* This routine uses the cholesky decomposition to invert a real */
/* symmetric matrix. */

    /* Parameter adjustments */
    amat_dim1 = *m;
    amat_offset = 1 + amat_dim1;
    amat -= amat_offset;

    /* Function Body */
    *kxx = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kl = i__ - 1;
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    if (i__ > 1) {
		i__3 = kl;
		for (k = 1; k <= i__3; ++k) {
/* L3: */
		    amat[i__ + j * amat_dim1] -= amat[k + i__ * amat_dim1] * 
			    amat[k + j * amat_dim1];
		}
	    } else {
		if (amat[i__ + i__ * amat_dim1] <= 0.) {
		    *kxx = i__;
		    goto L20;
		}
	    }
	    if (i__ == j) {
		amat[i__ + i__ * amat_dim1] = sqrt(amat[i__ + i__ * amat_dim1]
			);
	    } else {
		if (j == i__ + 1) {
		    ooa = 1. / amat[i__ + i__ * amat_dim1];
		}
		amat[i__ + j * amat_dim1] *= ooa;
	    }
/* L7: */
	}
/* L8: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    ooa = 1. / amat[j + j * amat_dim1];
	    if (i__ >= j) {
		t = 1.;
		goto L12;
	    }
	    kl = j - 1;
	    t = 0.;
	    i__3 = kl;
	    for (k = i__; k <= i__3; ++k) {
/* L11: */
		t -= amat[i__ + k * amat_dim1] * amat[k + j * amat_dim1];
	    }
L12:
	    amat[i__ + j * amat_dim1] = t * ooa;
	}
/* L13: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    t = 0.;
	    i__3 = *n;
	    for (k = j; k <= i__3; ++k) {
/* L14: */
		t += amat[i__ + k * amat_dim1] * amat[j + k * amat_dim1];
	    }
	    amat[i__ + j * amat_dim1] = t;
/* L19: */
	    amat[j + i__ * amat_dim1] = t;
/* L15: */
	}
/* L16: */
    }
L20:
    return 0;
} /* cholx_ */

/* Subroutine */ int ddnor_(doublereal *ystar, doublereal *gauss)
{
    /* Initialized data */

    static doublereal c__[5] = { 3209.377589138469472562,
	    377.4852376853020208137,113.8641541510501556495,
	    3.161123743870565596947,.185777706184603152673 };
    static doublereal d__[4] = { 2844.236833439170622273,
	    1282.616526077372275645,244.0246379344441733056,
	    23.60129095234412093499 };
    static doublereal orpi = .5641895835477562869483;
    static doublereal root2 = .70710678118654752440083;
    static doublereal p[6] = { -6.58749161529837803157e-4,
	    -.0160837851487422766278,-.125781726111229246204,
	    -.360344899949804439429,-.305326634961232344035,
	    -.0163153871373020978498 };
    static doublereal q[5] = { .00233520497626869185443,
	    .0605183413124413191178,.527905102951428412248,
	    1.87295284992346047209,2.56852019228982242072 };
    static doublereal a[9] = { 1230.33935479799725272,2051.07837782607146532,
	    1712.04761263407058314,881.952221241769090411,
	    298.635138197400131132,66.1191906371416294775,
	    8.88314979438837594118,.56418849698867008918,
	    2.15311535474403846343e-8 };
    static doublereal b[8] = { 1230.33935480374942043,3439.36767414372163696,
	    4362.6190901432471582,3290.79923573345962678,
	    1621.38957456669018874,537.181101862009857509,
	    117.693950891312499305,15.7449261107098347253 };

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal x, y, x2, x3, x4, x5, x6, x7, x8, xm2, xm4, xm6, xm8, 
	    erf, bot, xm10;
    static integer isw;
    static doublereal top, erfc, crap;


/* Copyright (c) James G. MacKinnon, 1993 */
/* Routine to evaluate cumulative normal distribution */
/* Written originally in late 1970's */
/* Modified 1993 to avoid changing the argument */

/* This subroutine uses Cody's method to evaluate the cumulative */
/* normal distribution. It is probably accurate to 19 or 20 */
/* significant digits. It was written in 1977, based on the Cody */
/* article referred to in the documentation for IMSL subroutine mdnor. */


    isw = 1;
    y = *ystar;
    if (*ystar < -16.) {
	y = -16.;
    }
    if (*ystar > 16.) {
	y = 16.;
    }
    x = -y * root2;
    if (x > 0.) {
	goto L1;
    }
    if (x < 0.) {
	goto L2;
    }
    *gauss = .5;
    return 0;
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

/*  evaluate erfc for x.gt.4.0 */

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
    *gauss = erfc * .5;
    return 0;
L20:

/*  evaluate erfc for .477.lt.x.le.4.0 */

    x2 = x * x;
    x3 = x2 * x;
    x4 = x2 * x2;
    x5 = x3 * x2;
    x6 = x3 * x3;
    x7 = x3 * x4;
    x8 = x4 * x4;
    top = a[0] + a[1] * x + a[2] * x2 + a[3] * x3 + a[4] * x4 + a[5] * x5 + a[
	    6] * x6 + a[7] * x7 + a[8] * x8;
    bot = b[0] + b[1] * x + b[2] * x2 + b[3] * x3 + b[4] * x4 + b[5] * x5 + b[
	    6] * x6 + b[7] * x7 + x8;
    erfc = exp(-x2) * top / bot;

    if (isw == -1) {
	erfc = 2. - erfc;
    }
    *gauss = erfc * .5;
    return 0;
L10:

/*  evaluate erf for x.lt..477 */

    x2 = x * x;
    x4 = x2 * x2;
    x6 = x4 * x2;
    x8 = x4 * x4;
    top = c__[0] + c__[1] * x2 + c__[2] * x4 + c__[3] * x6 + c__[4] * x8;
    bot = d__[0] + d__[1] * x2 + d__[2] * x4 + d__[3] * x6 + x8;
    erf = x * top / bot;

    erf *= isw;
    erfc = 1. - erf;
    *gauss = erfc * .5;
    return 0;
} /* ddnor_ */

/* Subroutine */ int innorz_(doublereal *prob, doublereal *anorm)
{
    /* Initialized data */

    static doublereal c0 = 2.515517;
    static doublereal d1 = 1.432788;
    static doublereal c1 = .802853;
    static doublereal c2 = .010328;
    static doublereal d3 = .001308;
    static doublereal d2 = .189269;
    static doublereal const__ = .398942280401432678;

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double log(doublereal), sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal t, pr, pr2, arg, dens, prob2;
    extern /* Subroutine */ int ddnor_(doublereal *, doublereal *);
    static doublereal error, anorm2;

    /* Fortran I/O blocks */
    static cilist io___131 = { 0, 6, 0, 0, 0 };

/* Copyright (c) James G. MacKinnon, 1995 */
/* Inverse normal routine that adjusts crude result twice. */
/* It seems to be accurate to about 14 digits. */
/* Crude result is taken from Abramowitz & Stegun (1968) */
/* It should have abs. error < 4.5 * 10^-4 */

    if (*prob < 0. || *prob > 1.) {
	s_wsle(&io___131);
	do_lio(&c__9, &c__1, "Attempt to find inverse normal of ", (ftnlen)34)
		;
	do_lio(&c__5, &c__1, (char *)&(*prob), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    pr = *prob;
    if (*prob > .5) {
	pr = 1. - *prob;
    }
/* Computing 2nd power */
    d__1 = pr;
    arg = 1 / (d__1 * d__1);
    t = sqrt(log(arg));
/* Computing 2nd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
/* Computing 3rd power */
    d__3 = t;
    *anorm = t - (c0 + c1 * t + c2 * (d__1 * d__1)) / (d1 * t + 1 + d2 * (
	    d__2 * d__2) + d3 * (d__3 * (d__3 * d__3)));

/* now correct crude result by direct method */

    ddnor_(anorm, &prob2);
    pr2 = 1. - prob2;
/* Computing 2nd power */
    d__1 = pr2;
    arg = 1 / (d__1 * d__1);
    t = sqrt(log(arg));
/* Computing 2nd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
/* Computing 3rd power */
    d__3 = t;
    anorm2 = t - (c0 + c1 * t + c2 * (d__1 * d__1)) / (d1 * t + 1 + d2 * (
	    d__2 * d__2) + d3 * (d__3 * (d__3 * d__3)));
    *anorm = *anorm + *anorm - anorm2;
    if (*prob < .5) {
	*anorm = -(*anorm);
    }

/* now correct better result, using Taylor series approximation */

    ddnor_(anorm, &prob2);
    error = prob2 - *prob;
/* Computing 2nd power */
    d__1 = *anorm;
    dens = const__ * exp(d__1 * d__1 * -.5);
    *anorm -= error / dens;
    return 0;
} /* innorz_ */

double mackinnon_pvalue (double tval, int n, const char *path)
{
    integer niv = 1; /* number of variables */
    integer itt = 1; /* tau test (2 for z test) */
    integer itv = 2; /* model "c", with constant */
    integer nobs = n;
    
    doublereal arg = tval;
    doublereal val;
    
    int check;

    check = urcval_(&niv, &itt, &itv, &nobs, &arg, &val, path);

    if (check != URC_OK && check != URC_SMALL_SAMPLE) {
	val = NADBL;
    }

    return (double) val;
}

