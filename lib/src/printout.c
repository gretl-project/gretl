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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/*  printout.c - simple text print routines for some gretl structs */ 

#include "libgretl.h"
#include "version.h"
#include "libset.h"
#include "forecast.h"

#include <time.h>

#undef PRN_DEBUG

void bufspace (int n, PRN *prn)
{
    while (n-- > 0) {
	pputc(prn, ' ');
    }
}

/**
 * printxx:
 * @xx: number to print.
 * @str: buffer into which to print.
 * @ci: command index (PRINT or SUMMARY).
 *
 * Print a string representation of the double-precision value @xx
 * to the buffer @str, in a format that depends on @ci.
 */

static void printxx (const double xx, char *str, int ci)
{
    int d = (ci == PRINT)? 8 : 6;

    sprintf(str, "%#*.*g", d, GRETL_DIGITS, xx);
}

static void covhdr (PRN *prn)
{
    pprintf(prn, "\n%s\n\n", _("COVARIANCE MATRIX OF REGRESSION COEFFICIENTS"));
}

/**
 * session_time:
 * @prn: where to print.
 *
 * Print the current time to the specified printing object,
 * or to %stdout if @prn is %NULL.
 */

void session_time (PRN *prn)
{
    time_t runtime = time(NULL);
    PRN *myprn = NULL;

    if (prn == NULL) {
	myprn = gretl_print_new(GRETL_PRINT_STDOUT);
	prn = myprn;
    }

    pprintf(prn, "%s: %s\n", _("Current session"), print_time(&runtime));
    
    if (myprn != NULL) {
	gretl_print_destroy(myprn);
    }    
}

/**
 * logo:
 *
 * Print to stdout gretl version information.
 */

void logo (void)
{
    printf(_("gretl version %s\n"), version_string);
    puts(_("Copyright Ramu Ramanathan and Allin Cottrell"));
    puts(_("This is free software with ABSOLUTELY NO WARRANTY"));
}

/**
 * gui_logo:
 * @prn: where to print.
 *
 * Print gretl GUI version information to the specified printing
 * object, or to %stdout if @prn is %NULL.
 */

void gui_logo (PRN *prn)
{
    PRN *myprn = NULL;

    if (prn == NULL) {
	myprn = gretl_print_new(GRETL_PRINT_STDOUT);
	prn = myprn;
    }
	
    pprintf(prn, _("gretl: gui client for gretl version %s,\n"), version_string);
    pputs(prn, _("copyright Allin Cottrell.\n"));
    pputs(prn, _("This is free software with ABSOLUTELY NO WARRANTY.\n"));

    if (myprn != NULL) {
	gretl_print_destroy(myprn);
    }
}

/**
 * lib_logo:
 *
 * Print gretl library version information to stdout.
 */

void lib_logo (void)
{
    printf("\nLibgretl-1.0, revision %d\n", LIBGRETL_REVISION);
}

/**
 * gui_script_logo:
 * @prn: gretl printing struct.
 *
 * Print to @prn a header for script output in gui program.
 */

void gui_script_logo (PRN *prn)
{
    time_t runtime = time(NULL);

    pprintf(prn, _("gretl version %s\n"), version_string);
    pprintf(prn, "%s: %s\n", _("Current session"), print_time(&runtime));
}

/* -------------------------------------------------------- */

static void 
print_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    pprintf(prn, " %8s ", cf->names[i]);

    bufspace(3, prn);

    if (isnan(cf->coeff[i])) {
	pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 16), _("undefined"));
    } else {
	gretl_print_value(cf->coeff[i], prn);
    }

    bufspace(2, prn);

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 10), _("undefined"));
    } else {
	pprintf(prn, " (%#.*g, %#.*g)", 
		GRETL_DIGITS, cf->coeff[i] - cf->maxerr[i],
		GRETL_DIGITS, cf->coeff[i] + cf->maxerr[i]);
    }

    pputc(prn, '\n');
}

/**
 * text_print_model_confints:
 * @cf: pointer to confidence intervals.
 * @prn: gretl printing struct.
 *
 * Print to @prn the 95 percent confidence intervals for parameter
 * estimates contained in @cf.
 */

void text_print_model_confints (const CoeffIntervals *cf, PRN *prn)
{
    int i;

    pprintf(prn, "t(%d, .025) = %.3f\n\n", cf->df, tcrit95(cf->df));
    /* xgettext:no-c-format */
    pputs(prn, _("      VARIABLE      COEFFICIENT      95% CONFIDENCE "
	    "INTERVAL\n\n"));      

    for (i=0; i<cf->ncoeff; i++) {
	print_coeff_interval(cf, i, prn);
    }

    pputc(prn, '\n');
}

/**
 * printcorr:
 * @corrmat: gretl correlation matrix struct.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print correlation matrix to @prn in a simple columnar format.
 * 
 */

void printcorr (const CorrMat *corrmat, const DATAINFO *pdinfo, 
		PRN *prn)
{
    int i = 1, j, k = 0, m, ncoeffs;
    char corrstring[32];

    m = corrmat->list[0];
    ncoeffs = (m * (m + 1)) / 2;

    pputs(prn, _("\nPairwise correlation coefficients:\n\n"));

    while (k < ncoeffs) {
        for (i=1; i<=m; i++) {
	    k++;
	    for (j=i+1; j<=m; j++) {
		sprintf(corrstring, "corr(%s, %s)", 
			pdinfo->varname[corrmat->list[i]], 
			pdinfo->varname[corrmat->list[j]]);
		if (na(corrmat->xpx[k])) {
		    pprintf(prn, "  %-24s    %s\n", 
			    corrstring, _("undefined"));
		} else if (corrmat->xpx[k] < 0.) {
		    pprintf(prn, "  %-24s = %.4f\n", corrstring, 
			    corrmat->xpx[k]);
		} else {
		    pprintf(prn, "  %-24s =  %.4f\n", corrstring, 
			    corrmat->xpx[k]);
		}
		k++;
	    }
        }
    }

    pputc(prn, '\n');
}

static void print_freq_test (const FREQDIST *freq, PRN *prn)
{
    double pval = NADBL;

    if (freq->dist == DIST_NORMAL) {
	pval = chisq(freq->test, 2);
	pprintf(prn, "\n%s:\n", 
		_("Test for null hypothesis of normal distribution"));
	pprintf(prn, "%s(2) = %.3f %s %.5f\n", 
		_("Chi-square"), freq->test, 
		_("with p-value"), pval);
    } else if (freq->dist == DIST_GAMMA) {
	pval = normal_pvalue_2(freq->test);
	pprintf(prn, "\n%s:\n", 
		_("Test for null hypothesis of gamma distribution"));
	pprintf(prn, "z = %.3f %s %.5f\n", freq->test, 
		_("with p-value"), pval);
    }	

    if (!na(pval)) {
	record_test_result(freq->test, pval, 
			   (freq->dist == DIST_NORMAL)? 
			   "normality" : "gamma");
    }
}

/**
 * print_freq:
 * @freq: gretl frequency distribution struct.
 * @prn: gretl printing struct.
 *
 * Print frequency distribution to @prn.
 * 
 */

void print_freq (const FREQDIST *freq, PRN *prn)
{
    int i, k, nlw, K = freq->numbins - 1;
    char word[64];

    if (freq == NULL) {
	return;
    }

    pprintf(prn, _("\nFrequency distribution for %s, obs %d-%d "
		   "(%d valid observations)\n"),
	    freq->varname, freq->t1 + 1, freq->t2 + 1, freq->n);

    if (freq->numbins == 0) {
	if (!na(freq->test)) {
	    print_freq_test(freq, prn);
	}
	return;
    } 

    pprintf(prn, _("number of bins = %d, mean = %g, sd = %g\n"), 
	    freq->numbins, freq->xbar, freq->sdx);

    pputs(prn, _("\n       interval          midpt      frequency\n\n"));

    for (k=0; k<=K; k++) {
	*word = '\0';
	if (k == 0) {
	    pputs(prn, "          <  ");
	} else if (k == K) {
	    pputs(prn, "          >= ");
	} else {
	    pprintf(prn, "%10g - ", freq->endpt[k]);
	}
	if (k == K) {
	    sprintf(word, "%g", freq->endpt[k]);
	} else {
	    sprintf(word, "%g", freq->endpt[k+1]);
	}
	pprintf(prn, "%s", word);

	nlw = 10 - strlen(word);
	bufspace(nlw, prn);

	sprintf(word, " %g", freq->midpt[k]);
	pputs(prn, word);

	nlw = 10 - strlen(word);
	bufspace(nlw, prn);

	pprintf(prn, "%6d  ", freq->f[k]);
	i = 36.0 * freq->f[k]/freq->n;
	while (i--) {
	    pputc(prn, '*');
	}
	pputc(prn, '\n');
    }

    if (!na(freq->test)) {
	print_freq_test(freq, prn);
    }
}

/**
 * print_smpl:
 * @pdinfo: data information struct
 * @fulln: full length of data series.
 * @prn: gretl printing struct.
 *
 * Print current sample information to @prn.
 * 
 */

void print_smpl (const DATAINFO *pdinfo, int fulln, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];

    if (fulln && !dataset_is_panel(pdinfo)) {
	pprintf(prn, _("Full data set: %d observations\n"),
		fulln);
	pprintf(prn, _("Current sample: %d observations\n"),
		pdinfo->n);
	return;
    }

    ntodate_full(date1, pdinfo->t1, pdinfo);
    ntodate_full(date2, pdinfo->t2, pdinfo);

    if (fulln) {
	pprintf(prn, _("Full data set: %d observations\n"), fulln);
    } else {
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Full data range"), 
		pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    }

    pprintf(prn, "%s:  %s - %s", _("Current sample"), date1, date2);

    if (pdinfo->t1 > 0 || pdinfo->t2 < pdinfo->n - 1 ||
	(fulln && dataset_is_panel(pdinfo))) {
	pprintf(prn, " (n = %d)\n", pdinfo->t2 - pdinfo->t1 + 1);
    } else {
	pputc(prn, '\n');
    } 
}

/**
 * gretl_fix_exponent:
 * @s: string representation of floating-point number.
 * 
 * Some C libraries (e.g. MS) print an "extra" zero in the exponent
 * when using scientific notation, e.g. "1.45E-002".  This function
 * checks for this and cuts it out if need be.
 *
 * Returns: the fixed numeric string.
 */

char *gretl_fix_exponent (char *s)
{
    char *p;

    if ((p = strstr(s, "+00")) || (p = strstr(s, "-00"))) {
	memmove(p+1, p+2, strlen(p+1));
    }

    return s;
}

/* For some reason sprintf using "%#G" seems to stick an extra
   zero on the end of some numbers -- i.e. when using a precision
   of 6 you can get a result of "1.000000", with 6 trailing
   zeros.  The following function checks for this and lops it
   off if need be. */

static void cut_extra_zero (char *numstr, int digits)
{
    if (strchr(numstr, 'E') == NULL && strchr(numstr, 'e') == NULL) {
	int s = strspn(numstr, "-.,0");
	int p = (strchr(numstr + s, '.') || strchr(numstr + s, ','));

	numstr[s + p + digits] = '\0';
    }
}

/* The following function formats a double in such a way that the
   decimal point will be printed in the same position for all
   numbers printed this way.  The total width of the number
   string (including possible padding on left or right) is 
   2*P + 5 characters, where P denotes the precision ("digits"). 
*/

void gretl_print_fullwidth_double (double x, int digits, PRN *prn)
{
    char numstr[36], final[36];
    char *p;
    int i, tmp, forept = 0;
    char decpoint = '.';

#ifdef ENABLE_NLS
    decpoint = get_local_decpoint();
#endif

    /* let's not print non-zero values for numbers smaller than
       machine zero */
    x = screen_zero(x);

    sprintf(numstr, "%#.*G", digits, x);

    gretl_fix_exponent(numstr);

    p = strchr(numstr, decpoint);
    if (p != NULL) {
	forept = p - numstr;
    }
    tmp = digits + 1 - forept;
    *final = 0;
    for (i=0; i<tmp; i++) {
	strcat(final, " ");
    }

    tmp = strlen(numstr) - 1;
    if (numstr[tmp] == decpoint) {
	numstr[tmp] = 0;
    }

    cut_extra_zero(numstr, digits);

    strcat(final, numstr);

    tmp = 2 * digits + 5 - strlen(final);
    for (i=0; i<tmp; i++) {
	strcat(final, " ");
    }

    pputs(prn, final);
}

/* ......................................................... */ 

void gretl_print_value (double x, PRN *prn)
{
    gretl_print_fullwidth_double(x, GRETL_DIGITS, prn);  
}

/* ......................................................... */ 


/**
 * outcovmx:
 * @pmod: pointer to model.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 * 
 * Print to @prn the variance-covariance matrix for the parameter
 * estimates in @pmod.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int outcovmx (MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int k, nbeta = 0;
    int *tmplist = NULL;

    if (pmod->ci == TSLS) {
	k = 2;
	while (pmod->list[k++] != LISTSEP) {
	    nbeta++;
	}
    } else if (pmod->ci == ARMA || pmod->ci == GARCH) {
	nbeta = pmod->ifc + pmod->list[1] + pmod->list[2] + pmod->list[0] - 4;
    } else {
	nbeta = pmod->list[0] - 1;
    }

    tmplist = gretl_list_new(nbeta);
    if (tmplist == NULL) {
	return E_ALLOC;
    }

    for (k=1; k<=tmplist[0]; k++) {
	tmplist[k] = pmod->list[k+1];
    }

    if (pmod->vcv == NULL) {
	if (makevcv(pmod)) {
	    return E_ALLOC;
	}
    } 

    text_print_matrix(pmod->vcv, tmplist, pmod, pdinfo, prn);  

    free(tmplist);

    return 0;
}

/* ......................................................... */ 

static void outxx (const double xx, int ci, PRN *prn)
{
    if (isnan(xx) || na(xx)) { 
	if (ci == CORR) {
	    pprintf(prn, " %*s", UTF_WIDTH(_("undefined"), 13), 
		    _("undefined"));
	} else {
	    pputs(prn, "              ");
	}
    }
	
    else if (ci == CORR) {
	pprintf(prn, " %13.4f", xx);
    }

    else {
	char numstr[18];

	if (xx > -0.001 && xx < 0.001) {
	    sprintf(numstr, "%.5e", xx);
	} else {
	    sprintf(numstr, "%g", xx);
	}
	gretl_fix_exponent(numstr);
	pprintf(prn, "%14s", numstr);
    }
}

/* ......................................................... */ 

static int takenotes (int quit_opt)
{
    char resp[4];

    if (quit_opt) {
	puts(_("\nTake notes then press return key to continue (or q to quit)"));
    } else {
	puts(_("\nTake notes then press return key to continue"));
    }

    fflush(stdout);

    fgets(resp, sizeof resp, stdin);

    if (quit_opt && *resp == 'q') {
	return 1;
    }

    return 0;
}

/**
 * scroll_pause:
 * 
 * Pause after a "page" of text at the console.
 */

void scroll_pause (void)
{
    takenotes(0);
}

/**
 * scroll_pause_or_quit:
 * 
 * Pause after a "page" of text, and give the user the option of
 * breaking out of the printing routine.
 * 
 * Returns: 1 if the user chose to quit, otherwise 0.
 */

int scroll_pause_or_quit (void)
{
    return takenotes(1);
}

/*  Given a one dimensional array which represents a symmetric
    matrix, prints out an upper triangular matrix of any size.

    Due to screen and printer column limitations the program breaks up
    a large upper triangular matrix into 5 variables at a time. For
    example, if there were 10 variables the program would first print
    an upper triangular matrix of the first 5 rows and columns, then
    it would print a rectangular matrix of the first 5 rows but now
    columns 6 - 10, and finally an upper triangular matrix of rows 6 -
    10 and columns 6 - 10
*/

void text_print_matrix (const double *rr, const int *list, 
			MODEL *pmod, const DATAINFO *pdinfo, 
			PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, m, index, ij2, lineno = 0;
    int nls = (pmod != NULL && pmod->ci == NLS);
    int arma = (pmod != NULL && pmod->ci == ARMA);
    int garch = (pmod != NULL && pmod->ci == GARCH);
    int pause = gretl_get_text_pause();
    char s[16];
    enum { FIELDS = 5 };

    if (pmod != NULL) covhdr(prn);

    m = 1;
    lo = list[0];

    for (i=0; i<=lo/FIELDS; i++) {
	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	if (pause && i > 0) {
	    takenotes(0);
	}

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    if (nls || arma || garch) {
		ljnf = j + nf;
		strcpy(s, pmod->params[ljnf]);
	    } else {
		ljnf = list[j + nf];
		strcpy(s, pdinfo->varname[ljnf]);
	    }
	    bufspace(9 - strlen(s), prn);
	    pprintf(prn, "%3d) %s", ljnf, s);
	}
	pputc(prn, '\n');
	lineno += 2;

	/* print rectangular part, if any, of matrix */
	lineno = 1;
	for (j=0; j<nf; j++) {
	    if (pause && (lineno % PAGELINES == 0)) {
		takenotes(0);
		lineno = 1;
	    }
	    for (k=0; k<p; k++) {
		index = ijton(j, nf+k, lo);
		outxx(rr[index], (pmod == NULL)? CORR : 0, prn);
	    }
	    pprintf(prn, "   (%d\n", (nls || arma || garch)? 
		    j+1 : list[j+1]);
	    lineno++;
	}

	/* print upper triangular part of matrix */
	lineno = 1;
	for (j=0; j<p; ++j) {
	    if (pause && (lineno % PAGELINES == 0)) {
		takenotes(0);
		lineno = 1;
	    }
	    ij2 = nf + j;
	    bufspace(14 * j, prn);
	    for (k=j; k<p; k++) {
		index = ijton(ij2, nf+k, lo);
		outxx(rr[index], (pmod == NULL)? CORR : 0, prn);
	    }
	    pprintf(prn, "   (%d\n", (nls || arma || garch)? 
		    ij2+1 : list[ij2+1]);
	    lineno++;
	}
	pputc(prn, '\n');
    }
}

/* ........................................................... */

static void fit_resid_head (const FITRESID *fr, 
			    const DATAINFO *pdinfo, 
			    PRN *prn)
{
    int i;
    char label[16];
    char mdate1[OBSLEN], mdate2[OBSLEN];

    ntodate(mdate1, fr->t1, pdinfo);
    ntodate(mdate2, fr->t2, pdinfo);

    pprintf(prn, _("Model estimation range: %s - %s"), mdate1, mdate2);
    pprintf(prn, " (n = %d)\n", fr->real_nobs);

    if (!na(fr->sigma)) {
	pprintf(prn, _("Standard error of residuals = %f\n"), fr->sigma);
    }
    
    pprintf(prn, "\n     %s ", _("Obs"));

    for (i=1; i<4; i++) {
	if (i == 1) strcpy(label, fr->depvar);
	if (i == 2) strcpy(label, _("fitted"));
	if (i == 3) strcpy(label, _("residuals"));
	pprintf(prn, "%*s", UTF_WIDTH(label, 13), label); 
    }

    pputs(prn, "\n\n");
}

/*  skips to new page and prints names of variables
    from v1 to v2 */

static void varheading (int v1, int v2, 
			const DATAINFO *pdinfo, const int *list,
			PRN *prn)
{
    int i;
        
    pputs(prn, "\n     Obs ");

    for (i=v1; i<=v2; i++) { 
	pprintf(prn, "%13s", pdinfo->varname[list[i]]);
    }

    pputs(prn, "\n\n");
}

/**
 * gretl_printxn:
 * @x: number to print.
 * @n: controls width of output.
 * @prn: gretl printing struct.
 *
 * Print a string representation of the double-precision value @x
 * in a format that depends on @n.
 */

void gretl_printxn (double x, int n, PRN *prn)
{
    char s[32];
    int ls;

    if (na(x)) {
	*s = '\0';
    } else {
	printxx(x, s, PRINT);
    }

    ls = strlen(s);

    pputc(prn, ' ');
    bufspace(n - 3 - ls, prn);
    pputs(prn, s);
}

/* ........................................................... */

static void printstr_ten (PRN *prn, double xx, int *ls)
{
    int lwrd;    
    char str[32];

    if (na(xx)) {
	strcpy(str, "NA");
    } else {
	sprintf(str, "%.10g", xx);
    }
    strcat(str, "  ");
    lwrd = strlen(str);
    if (*ls + lwrd > 78) {
	*ls = 0;
	pputc(prn, '\n');
    }
    pputs(prn, str);
    *ls += lwrd;
}

/* ........................................................ */

static void printstr (PRN *prn, double xx, int *ls)
{
    int lwrd;
    char str[32];

    if (na(xx)) {
	strcpy(str, "NA");
    } else {
	printxx(xx, str, 0);
    }
    strcat(str, "  ");
    lwrd = strlen(str);
    if (*ls + lwrd > 78) {
	*ls = 0;
	pputc(prn, '\n');
    }
    pputs(prn, str);
    *ls += lwrd;
}

/* prints series z from current sample t1 to t2 */

static void printz (const double *z, const DATAINFO *pdinfo, 
		    PRN *prn, gretlopt opt)
{
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2, ls = 0;
    double xx;

    if (gretl_isconst(t1, t2, z)) {
	if (opt & OPT_T) {
	    printstr_ten(prn, z[t1], &ls);
	} else {
	    printstr(prn, z[t1], &ls);
	}
    }
    else for (t=t1; t<=t2; t++) {
	xx = z[t];
	if (opt & OPT_T) {
	    printstr_ten(prn, xx, &ls);
	} else {
	    printstr(prn, xx, &ls);
	}
    }

    pputc(prn, '\n');
}

#define SMAX 7            /* stipulated max. significant digits */
#define TEST_PLACES 12    /* # of decimal places to use in test string */

static int get_signif (const double *x, int n)
     /* return either (a) the number of significant digits in
	a data series (+), or (b) the number of decimal places to
	use when printing the series (-) */
{
    static char numstr[48];
    int i, j, s, smax = 0; 
    int lead, leadmax = 0, leadmin = 99;
    int gotdec, trail, trailmax = 0;
    double xx;
    int allfrac = 1;
    char decpoint = '.';

#ifdef ENABLE_NLS
    decpoint = get_local_decpoint();
#endif

    for (i=0; i<n; i++) {

	if (na(x[i])) continue;

	xx = fabs(x[i]);
	if (xx >= 1.0) allfrac = 0;

	sprintf(numstr, "%.*f", TEST_PLACES, xx);
	s = strlen(numstr) - 1;
	trail = TEST_PLACES;
	gotdec = 0;

	for (j=s; j>0; j--) {
	    if (numstr[j] == '0') {
		s--;
		if (!gotdec) trail--;
	    } else if (numstr[j] == decpoint) {
		gotdec = 1;
		if (xx < 10000) break;
		else continue;
	    } else {
		break;
	    }
	}

	if (trail > trailmax) {
	    trailmax = trail;
	}

	if (xx < 1.0) s--; /* don't count leading zero */

	if (s > smax) {
	    smax = s;
	}

#ifdef PRN_DEBUG
	fprintf(stderr, "get_signif: set smax = %d\n", smax);
#endif

	lead = 0;
	for (j=0; j<=s; j++) {
	    if (xx >= 1.0 && numstr[j] != decpoint) lead++;
	    else break;
	}

	if (lead > leadmax) leadmax = lead;
	if (lead < leadmin) leadmin = lead;
    } 

    if (smax > SMAX) smax = SMAX;

    if (trailmax > 0 && (leadmax + trailmax <= SMAX)) {
	smax = -trailmax;
    } else if ((leadmin < leadmax) && (leadmax < smax)) {
#ifdef PRN_DEBUG
	fprintf(stderr, "get_signif: setting smax = -(%d - %d)\n", 
		smax, leadmax);
#endif	
	smax = -1 * (smax - leadmax); /* # of decimal places */
    } else if (leadmax == smax) {
	smax = 0;
    } else if (leadmax == 0 && !allfrac) {
#ifdef PRN_DEBUG
	fprintf(stderr, "get_signif: setting smax = -(%d - 1)\n", smax);
#endif
	smax = -1 * (smax - 1);
    } 

    return smax;
}

/* ........................................................... */

static int g_too_long (double x, int signif)
{
    char n1[32], n2[32];

    sprintf(n1, "%.*G", signif, x);
    sprintf(n2, "%.0f", x);
    
    return (strlen(n1) > strlen(n2));
}

/* ........................................................... */

static int bufprintnum (char *buf, double x, int signif, int width)
{
    static char numstr[32];
    int i, l;

    /* guard against monster numbers that will smash the stack */
    if (fabs(x) > 1.0e20) {
	sprintf(numstr, "%g", x);
	goto finish;
    }

    if (signif < 0) {
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for signif: "
		    "printing with %%.%df\n", signif, -signif);
#endif
	sprintf(numstr, "%.*f", -signif, x);
    } else if (signif == 0) {
#ifdef PRN_DEBUG
	    fprintf(stderr, "got 0 for signif: "
		    "printing with %%.0f\n");
#endif
	sprintf(numstr, "%.0f", x);
    } else {
	double z = fabs(x);

	if (z < 1) l = 0;
	else if (z < 10) l = 1;
	else if (z < 100) l = 2;
	else if (z < 1000) l = 3;
	else if (z < 10000) l = 4;
	else if (z < 100000) l = 5;
	else if (z < 1000000) l = 6;
	else l = 7;

	if (l == 6 && signif < 6) {
	   sprintf(numstr, "%.0f", x); 
	} else if (l >= signif) { 
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%dG\n", l, signif, signif);
#endif
	    if (g_too_long(x, signif)) {
		sprintf(numstr, "%.0f", x);
	    } else {
		sprintf(numstr, "%.*G", signif, x);
	    }
	} else if (z >= .10) {
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%df\n", l, signif, signif-l);
#endif
	    sprintf(numstr, "%.*f", signif - l, x);
	} else {
	    if (signif > 4) signif = 4;
#ifdef PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%#.%dG\n", l, signif, signif);
#endif
	    sprintf(numstr, "%#.*G", signif, x); /* # wanted? */
	}
    }

 finish:

    l = width - strlen(numstr);
    for (i=0; i<l; i++) {
	strcat(buf, " ");
    }
    strcat(buf, numstr);

    return 0;
}

/**
 * print_obs_marker:
 * @t: observation number.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print a string (label, date or obs number) representing the given @t.
 *
 */

void print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn)
{
    if (pdinfo->markers) { 
	pprintf(prn, "%8s ", pdinfo->S[t]); 
    } else {
	char tmp[OBSLEN]; 

	ntodate(tmp, t, pdinfo);
	pprintf(prn, "%8s ", tmp);
    }
}

/* See if there is a variable that is an outcome of sorting,
   and has sorted case-markers attached to it.  If so,
   we'll arrange to print the case-markers in the correct
   sequence, provided the variable in question is being
   printed by itself, or as the last in a short list of
   variables.
*/

static int 
check_for_sorted_var (int *list, const DATAINFO *pdinfo)
{
    int i, v, ret = 0;
    int l0 = list[0];
    int pos = 0;

    if (l0 < 5 && !complex_subsampled()) {
	for (i=1; i<=l0; i++) {
	    v = list[i];
	    if (pdinfo->varinfo[v]->sorted_markers != NULL) {
		if (ret == 0) {
		    ret = v;
		    pos = i;
		} else {
		    ret = 0;
		    pos = 0;
		    break;
		}
	    }
	}
    }

    if (ret && pos != list[0]) {
	/* sorted var should be last in list */
	int tmp = list[l0];

	list[l0] = list[pos];
	list[pos] = tmp;
    }

    return ret;
}

static void print_scalar (double x, const char *vname, 
			  gretlopt opt, int allconst,
			  PRN *prn)
{
    if (!allconst) {
	pputc(prn, '\n');
    }

    if (na(x)) {
	pprintf(prn, "%8s = NA", vname);
    } else if (opt & OPT_T) {
	pprintf(prn, "%8s = %.10g", vname, x);
    } else {
	pprintf(prn, "%8s = %10g", vname, x);
    }

    if (allconst) {
	pputc(prn, '\n');
    }
}

/**
 * printdata:
 * @list: list of variables to print.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: if = OPT_O, print the data by observation (series in columns);
 *       if = OPT_T, print the data to 10 significant digits.
 * @prn: gretl printing struct.
 *
 * Print the data for the variables in @list, from observations t1 to
 * t2.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int printdata (const int *list, const double **Z, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    int j, v, v1, v2, j5, nvj5, lineno, ncol;
    int allconst, scalars = 0;
    int nvars = 0, sortvar = 0;
    int *plist = NULL;
    int *pmax = NULL; 
    int t, nsamp;
    char line[128];
    int err = 0;

    int pause = gretl_get_text_pause();

    if (list == NULL) {
	plist = full_var_list(pdinfo, &nvars);
    } else {
	nvars = list[0];
	if (nvars > 0) {
	    plist = gretl_list_copy(list);
	}
    }

    if (plist == NULL) {
	if (nvars == 0) {
	    pputs(prn, "No data\n");
	    goto endprint;
	} else {
	    return E_ALLOC;
	}
    }

    lineno = 1;

    /* screen out any scalars and print them first */
    for (j=1; j<=plist[0]; j++) {
	if (!pdinfo->vector[plist[j]]) {
	    print_scalar(Z[plist[j]][0], pdinfo->varname[plist[j]],
			 opt, 0, prn);
	    scalars = 1;
	    gretl_list_delete_at_pos(plist, j);
	    j--;
	} 
    }

    if (scalars) {
	pputc(prn, '\n');
    }

    /* special case: all vars have constant value over sample */
    allconst = 1;
    for (j=1; j<=plist[0]; j++) {
	double xx = Z[plist[j]][pdinfo->t1];

	for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) {
	    if (floatneq(Z[plist[j]][t], xx)) {
		allconst = 0;
		break;
	    }
	}
	if (!allconst) break;
    }

    if (allconst) {
	for (j=1; j<=plist[0]; j++) {
	    print_scalar(Z[plist[j]][pdinfo->t1], pdinfo->varname[plist[j]],
			 opt, 1, prn);
	}
	goto endprint;
    }

    if (!(opt & OPT_O)) { /* not by observations, but by variable */
	if (plist[0] > 0) {
	    pputc(prn, '\n');
	}
	for (j=1; j<=plist[0]; j++) {
	    pprintf(prn, _("Varname: %s\n"), pdinfo->varname[plist[j]]);
	    print_smpl(pdinfo, 0, prn);
	    pputc(prn, '\n');
	    printz(Z[plist[j]], pdinfo, prn, opt);
	    pputc(prn, '\n');
	}
	goto endprint;
    }

    pmax = malloc(plist[0] * sizeof *pmax);
    if (pmax == NULL) {
	err = E_ALLOC;
	goto endprint;
    }

    nsamp = pdinfo->t2 - pdinfo->t1 + 1;
    for (j=1; j<=plist[0]; j++) {
	/* this runs fairly quickly, even for large dataset */
	pmax[j-1] = get_signif(Z[plist[j]] + pdinfo->t1, nsamp);
    }

    sortvar = check_for_sorted_var(plist, pdinfo);

    /* print data by observations */
    ncol = 5;
    for (j=0; j<=plist[0]/ncol; j++) {
	char obs_string[OBSLEN];

	j5 = j * ncol;
	nvj5 = plist[0] - j5;
	v1 = j5 +1;
	if (nvj5) {
	    /* starting a new block of variables */
	    v2 = (ncol > nvj5)? nvj5 : ncol;
	    v2 += j5;
	    varheading(v1, v2, pdinfo, plist, prn);

	    if (pause && j > 0 && takenotes(1)) {
		goto endprint;
	    }

	    lineno = 1;

	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {

		if (sortvar && plist[0] == 1) {
		    strcpy(obs_string, SORTED_MARKER(pdinfo, sortvar, t));
		} else {
		    get_obs_string(obs_string, t, pdinfo);
		}
		
		sprintf(line, "%8s ", obs_string);
		
		for (v=v1; v<=v2; v++) {
		    double xx = Z[plist[v]][t];

		    if (na(xx)) {
			strcat(line, "             ");
		    } else { 
			bufprintnum(line, xx, pmax[v-1], 13);
		    }
		}

		if (sortvar && plist[0] > 1) {
		    sprintf(obs_string, "%8s", SORTED_MARKER(pdinfo, sortvar, t));
		    strcat(line, obs_string);
		}

		strcat(line, "\n");

		if (pputs(prn, line) < 0) {
		    err = E_ALLOC;
		    goto endprint;
		}

		if (pause && (lineno % PAGELINES == 0)) {
		    if (takenotes(1)) {
			goto endprint;
		    }
		    lineno = 1;
		}

		lineno++;
	    } /* end of printing obs (t) loop */
	} /* end if nvj5 */
    } /* end for j loop */

    pputc(prn, '\n');

 endprint:

    free(plist);
    free(pmax);

    return err;
}

int
text_print_fit_resid (const FITRESID *fr, const DATAINFO *pdinfo, PRN *prn)
{
    int t, anyast = 0;
    double xx;

    fit_resid_head(fr, pdinfo, prn); 

    for (t=0; t<fr->nobs; t++) {
	print_obs_marker(t + fr->t1, pdinfo, prn);

	if (na(fr->actual[t])) {
	    pputc(prn, '\n');
	} else if (na(fr->fitted[t])) {
	    pprintf(prn, "%13g\n", fr->actual[t]);
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) {
		anyast = 1;
	    }
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13.*f%13.*f%13.*f%s\n", 
			fr->pmax, fr->actual[t],
			fr->pmax, fr->fitted[t], fr->pmax, xx,
			(ast)? " *" : "");
	    } else {
		pprintf(prn, "%13g%13g%13g%s\n", 
			fr->actual[t],
			fr->fitted[t], xx,
			(ast)? " *" : "");
	    }
	}
    }

    pputc(prn, '\n');

    if (anyast) {
	pputs(prn, _("Note: * denotes a residual in excess of "
		     "2.5 standard errors\n"));
    }

    return 0;
}

/**
 * text_print_forecast:
 * @fr: pointer to structure containing forecasts.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: if includes %OPT_P, make a plot of the forecasts.
 * @prn: printing structure.
 *
 * Print the forecasts in @fr to @prn, and also plot the
 * forecasts if %OPT_P is given.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int text_print_forecast (const FITRESID *fr, 
			 double ***pZ, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn)
{
    int t, pv, err = 0;
    int do_errs = (fr->sderr != NULL);
    double *maxerr = NULL;
    int time_series = (pdinfo->structure == TIME_SERIES);
    int plot = (opt & OPT_P);

    if (do_errs) {
	maxerr = malloc(fr->nobs * sizeof *maxerr);
	if (maxerr == NULL) {
	    return E_ALLOC;
	}
    }

    if (do_errs) {
	/* FIXME arma */
	pprintf(prn, _(" For 95%% confidence intervals, t(%d, .025) = %.3f\n"), 
		fr->df, fr->tval);
    }

    pputs(prn, "\n     Obs ");
    pprintf(prn, "%12s", fr->depvar);
    pprintf(prn, "%*s", UTF_WIDTH(_("prediction"), 14), _("prediction"));

    if (do_errs) {
	pprintf(prn, "%*s", UTF_WIDTH(_(" std. error"), 14), _(" std. error"));
	pprintf(prn, _("   95%% confidence interval\n"));
    } else {
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    for (t=0; t<fr->nobs; t++) {
	print_obs_marker(t + fr->t1, pdinfo, prn);
	gretl_printxn(fr->actual[t], 15, prn);

	if (na(fr->fitted[t])) {
	    pputc(prn, '\n');
	    continue;
	}
	gretl_printxn(fr->fitted[t], 15, prn);

	if (do_errs) {
	    if (na(fr->sderr[t])) {
		maxerr[t] = NADBL;
	    } else {
		gretl_printxn(fr->sderr[t], 15, prn);
		maxerr[t] = fr->tval * fr->sderr[t];
		gretl_printxn(fr->fitted[t] - maxerr[t], 15, prn);
		pputs(prn, " -");
		gretl_printxn(fr->fitted[t] + maxerr[t], 10, prn);
	    }
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    /* do we really want a plot for non-time series? */

    if (plot && fr->nobs > 0) {
	if (time_series) {
	    switch (pdinfo->pd) {
	    case 1:
		pv = plotvar(pZ, pdinfo, "annual");
		break;
	    case 4:
		pv = plotvar(pZ, pdinfo, "qtrs");
		break;
	    case 12:
		pv = plotvar(pZ, pdinfo, "months");
		break;
	    case 24:
		pv = plotvar(pZ, pdinfo, "hrs");
		break;
	    case 10:
		pv = plotvar(pZ, pdinfo, "decdate");
		break;
	    default:
		pv = plotvar(pZ, pdinfo, "time");
	    }
	} else {
	    pv = plotvar(pZ, pdinfo, "index");
	}

	if (pv < 0) {
	    err = 1;
	} else {
	    err = plot_fcast_errs(fr->nobs, &(*pZ)[pv][fr->t1], 
				  fr->actual, fr->fitted, maxerr, 
				  fr->depvar, 
				  (time_series)? pdinfo->pd : 0);
	}
    }

    if (maxerr != NULL) {
	free(maxerr);
    }

    return err;
}

/**
 * print_fit_resid:
 * @pmod: pointer to gretl model.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print to @prn the fitted values and residuals from @pmod.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int print_fit_resid (const MODEL *pmod, const double **Z, 
		     const DATAINFO *pdinfo, PRN *prn)
{
    FITRESID *fr;

    fr = get_fit_resid(pmod, Z, pdinfo);
    if (fr == NULL) {
	return 1;
    }

    text_print_fit_resid(fr, pdinfo, prn);
    free_fit_resid(fr);

    return 0;
}

/* apparatus for user-defined printf statements */

#undef PRINTF_DEBUG

#define is_format_char(c) (c == 'f' || c == 'g' || c == 'd' || c == 's')

static int print_arg (const char **pfmt, double val, 
		      const char *str, PRN *prn)
{
    char fmt[16];
    int fc = *(*pfmt + 1);
    size_t n = 0;
    int err = 0;

    fmt[0] = '%';

    if (is_format_char(fc)) {
	fmt[1] = fc;
	fmt[2] = '\0';
	n = 2;
    } else {
	sscanf(*pfmt + 1, "%14[^gfsd]", fmt + 1);
	n = strlen(fmt);
	fc = *(*pfmt + n);
	fmt[n] = fc;
	fmt[++n] = '\0';
    }

    if (n == 0 || !is_format_char(fc)) {
	err = 1;
    } else if (fc == 's') {
	if (str == NULL) {
	    fprintf(stderr, "NULL string in printf\n");
	    err = 1;
	} else {
#ifdef PRINTF_DEBUG
	    fprintf(stderr, "fmt='%s', val = '%s'\n", fmt, str);
#endif
	    pprintf(prn, fmt, str);
	    *pfmt += n;	
	}
    } else {
#ifdef PRINTF_DEBUG
	if (fc == 'd') {
	    fprintf(stderr, "fmt='%s', val = %d\n", fmt, (int) val);
	} else {
	    fprintf(stderr, "fmt='%s', val = %g\n", fmt, val);
	}
#endif
	if (fc == 'd') {
	    pprintf(prn, fmt, (int) val);
	} else {
	    pprintf(prn, fmt, val);
	}
	*pfmt += n;
    }
    
    return err;
}

static int handle_escape (int c, PRN *prn)
{
    int err = 0;

    switch (c) {
    case 'n':
	pputc(prn, '\n');
	break;
    case 't':
	pputc(prn, '\t');
	break;
    case 'v':
	pputc(prn, '\v');
	break;
    case '\\':
	pputc(prn, '\\');
	break;
    default:
	err = 1;
    }

    return err;
}

static int output_format_only (const char *s, PRN *prn)
{
    int err = 0;

    while (*s) {
	if (*s == '\\') {
	    err = handle_escape(*(s+1), prn);
	    s++;
	} else {
	    pputc(prn, *s);
	}
	s++;
    }

    return err;
}

static char *get_arg (char *line)
{
    int inparen = 0;
    static char *p = NULL;
    char *q, *ret = NULL;

    if (line != NULL) p = line;

    q = p;
    while (*p && ret == NULL) {
	if (*p == '(') inparen++;
	else if (*p == ')') inparen--;
	if (!inparen && *p == ',') {
	    *p = '\0';
	    ret = q;
	}
	p++;
    }

    if (*p == '\0') ret = q;

    return ret;
}

static int get_marker_offset (const char *s)
{
    int off = 0;

    if (sscanf(s, "marker+%d", &off)) {
	if (off < 0) {
	    off = 0;
	}
    }

    return off;
}

static int real_do_printf (const char *line, double ***pZ, 
			   DATAINFO *pdinfo, MODEL *pmod,
			   PRN *prn, int t)
{
    const char *p;
    char format[128];
    char *argv, *str = NULL;
    double *vals = NULL;
    int argc = 0, cnvc = 0, inparen = 0;
    int markerpos = -1;
    int markeroffset = 0;
    int i, err = 0;

    if (t < 0) {
	t = pdinfo->t1;
    }

#ifdef PRINTF_DEBUG
    fprintf(stderr, "do_printf: line='%s'\n", line);
#endif

    *gretl_errmsg = '\0';

    if (!strncmp(line, "printf ", 7)) {
	line += 7;
    }

    if (sscanf(line, "\"%127[^\"]\"", format) != 1) {
	return 1;
    }

#ifdef PRINTF_DEBUG
    fprintf(stderr, "do_printf: format='%s'\n", format);
#endif

    p = format;
    while (*p) {
	if (*p == '%') {
	    if (*(p+1) == '%') p++;
	    else cnvc++;
	}
	p++;
    }

    line += strlen(format) + 2;
    if (*line != ',') {
	err = output_format_only(format, prn);
	return err;
    }

    line++;
    p = line;
    while (*p) {
	if (*p == '(') inparen++;
	else if (*p == ')') inparen--;
	if (!inparen && *p == ',') argc++;
	p++;
    }

    argc++;
    if (argc != cnvc) {
	fprintf(stderr, "do_printf: argc = %d but conversions = %d\n",
		argc, cnvc);
	err = 1;
	goto printf_bailout;
    }

    vals = malloc(argc * sizeof *vals);
    str = malloc(strlen(line) + 1);
    if (vals == NULL || str == NULL) {
	err = E_ALLOC;
	goto printf_bailout;
    }

    strcpy(str, line);

    for (i=0; i<argc; i++) {
	argv = get_arg((i > 0)? NULL : str);
	chopstr(argv);

#ifdef PRINTF_DEBUG
	fprintf(stderr, "do_printf: processing argv[%d] '%s'\n", i, argv);	
#endif

	if (numeric_string(argv)) {
	    vals[i] = atof(argv);
	} else if (!strncmp(argv, "marker", 6)) {
	    if (markerpos >= 0 || pdinfo->S == NULL) {
		err = 1;
	    } else {
		markerpos = i;
		vals[i] = 0.0;
		markeroffset = get_marker_offset(argv);
	    }
	} else {
	    int v = varindex(pdinfo, argv);

	    if (v < pdinfo->v) {
		/* simple existent varname */
		if (pdinfo->vector[v]) {
		    vals[i] = (*pZ)[v][t];
		} else {
		    vals[i] = (*pZ)[v][0];
		}
	    } else {
		err = get_generated_value(argv, &vals[i], pZ, pdinfo, 
					  pmod, t);
	    }
	}

#ifdef PRINTF_DEBUG
	fprintf(stderr, " after processing arg, vals[%d] = %g, err = %d\n", 
		i, vals[i], err);	
#endif

	if (err) {
	    goto printf_bailout;
	}
    }    

    p = format;
    i = 0;
    while (*p && !err) {
	const char *marker = NULL;

	if (*p == '%') {
	    if (*(p + 1) == '%') {
		pputc(prn, '%');
		p += 2;
	    } else {
		if (i == markerpos) {
		    marker = pdinfo->S[t];
		    if (markeroffset > 0 && 
			markeroffset < strlen(pdinfo->S[t])) {
			marker += markeroffset;
		    }
		} 
		err = print_arg(&p, vals[i++], marker, prn);
	    }
	} else if (*p == '\\') {
	    err = handle_escape(*(p + 1), prn);
	    p += 2;
	} else {
	    pputc(prn, *p);
	    p++;
	}
    }

    if (err) {
	pputc(prn, '\n');
    }

 printf_bailout:
    free(vals);
    free(str);

    return err;
}

int do_printf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, MODEL *pmod,
	       PRN *prn)
{
    return real_do_printf(line, pZ, pdinfo, pmod, prn, -1);
}

/* originating command is of form:

     genr markers=f1,f2,f3,...

   we're assuming that we're just getting the f* part here
*/

int generate_obs_markers (double ***pZ, DATAINFO *pdinfo, char *s)
{
    PRN *prn;
    int t, err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER);
    if (prn == NULL) {
	err = E_ALLOC;
    }

    if (!err && pdinfo->S == NULL) {
	err = dataset_allocate_obs_markers(pdinfo);
    }

    if (!err) {
	const char *buf;

	for (t=0; t<pdinfo->n && !err; t++) {
	    gretl_print_reset_buffer(prn);
	    err = real_do_printf(s, pZ, pdinfo, NULL, prn, t);
	    if (!err) {
		buf = gretl_print_get_buffer(prn);
		pdinfo->S[t][0] = '\0';
		strncat(pdinfo->S[t], buf, OBSLEN - 1);
	    }
	}
    }

    gretl_print_destroy(prn);
	
    return err;
}

int in_usa (void)
{
    static int ustime = -1;

    if (ustime < 0) {
	char test[12];
	struct tm t = {0};

	t.tm_year = 100;
	t.tm_mon = 0;
	t.tm_mday = 31;

	strftime(test, sizeof test, "%x", &t);

	if (!strncmp(test, "01/31", 5)) {
	    ustime = 1;
	} else {
	    ustime = 0;
	}
    }

    return ustime;
}

/**
 * bufgets:
 * @s: target string (must be pre-allocated)
 * @size: max number of characters to print
 * @buf: source buffer
 *
 * This function (which works rather line fgets) must be initialized 
 * via the call: bufgets(NULL, 0, buf);
 * 
 * Returns: s (NULL if nothing more can be read from buf).
 */

char *bufgets (char *s, size_t size, const char *buf)
{
    enum {
	END_OF_BUF,
	GOT_LF,
	GOT_CR,
	GOT_CRLF
    };
    int i, status = END_OF_BUF;
    static const char *p;

    /* mechanism for resetting p */
    if (s == NULL || size == 0) {
	p = NULL;
	return 0;
    }

    /* start at beginning of buffer */
    if (p == NULL) p = buf;

    /* signal that we've reached the end of the buffer */
    if (p && *p == 0) return NULL;

    *s = 0;
    /* advance to line-end, end of buffer, or maximum size,
       whichever comes first */
    for (i=0; ; i++) {
	s[i] = p[i];
	if (p[i] == 0) {
	    break;
	}
	if (p[i] == '\r') {
	    s[i] = 0;
	    if (p[i+1] == '\n') {
		status = GOT_CRLF;
	    } else {
		status = GOT_CR;
	    }
	    break;
	}
	if (p[i] == '\n') {
	    s[i] = 0;
	    status = GOT_LF;
	    break;
	}
	if (i == size - 1) {
	    fprintf(stderr, "bufgets: line too long: max %d characters\n", 
		    (int) size);
	    s[i] = '\0';
	    break;
	}
    }

    /* advance the buffer pointer */
    p += i;
    if (status == GOT_CR || status == GOT_LF) p++;
    else if (status == GOT_CRLF) p += 2;

    return s;
}
