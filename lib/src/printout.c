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
#include "gretl_func.h"

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
    pprintf(prn, "\n  %s\n\n", 
	    _("Covariance matrix of regression coefficients"));
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
    printf(_("gretl version %s\n"), GRETL_VERSION);
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
	
    pprintf(prn, _("gretl: gui client for gretl version %s,\n"), GRETL_VERSION);
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

    pprintf(prn, _("gretl version %s\n"), GRETL_VERSION);
    pprintf(prn, "%s: %s\n", _("Current session"), print_time(&runtime));
}

static void 
print_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    pprintf(prn, " %*s ", VNAMELEN - 1, cf->names[i]);

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
 * print_centered:
 * @s: string to print.
 * @width: width of field.
 * @prn: gretl printing struct.
 *
 * If the string @s is shorter than width, print it centered
 * in a field of the given width (otherwise just print it
 * straight).
 */

void print_centered (const char *s, int width, PRN *prn)
{
    int rem = width - strlen(s);

    if (rem <= 1) {
	pprintf(prn, "%s", s);
    } else {
	int i, off = rem / 2;

	for (i=0; i<off; i++) {
	    pputs(prn, " ");
	}
	pprintf(prn, "%-*s", width - off, s);
    }
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
    pputs(prn, _("      VARIABLE         COEFFICIENT      95% CONFIDENCE "
	    "INTERVAL\n\n"));      

    for (i=0; i<cf->ncoeff; i++) {
	print_coeff_interval(cf, i, prn);
    }

    pputc(prn, '\n');
}

/**
 * printcorr:
 * @corrmat: gretl correlation matrix struct.
 * @prn: gretl printing struct.
 *
 * Print correlation matrix to @prn in a simple columnar format.
 */

void printcorr (const VMatrix *corrmat, PRN *prn)
{
    int i, j, k = 0;
    int m, nterms;
    char corrstring[32];

    m = corrmat->dim;
    nterms = (m * (m + 1)) / 2;

    pputs(prn, _("\nPairwise correlation coefficients:\n\n"));

    while (k < nterms) {
        for (i=1; i<=m; i++) {
	    k++;
	    for (j=i+1; j<=m; j++) {
		sprintf(corrstring, "corr(%s, %s)", 
			corrmat->names[i-1], corrmat->names[j-1]);
		if (na(corrmat->vec[k])) {
		    pprintf(prn, "  %-24s    %s\n", 
			    corrstring, _("undefined"));
		} else if (corrmat->vec[k] < 0.0) {
		    pprintf(prn, "  %-24s = %.4f\n", corrstring, 
			    corrmat->vec[k]);
		} else {
		    pprintf(prn, "  %-24s =  %.4f\n", corrstring, 
			    corrmat->vec[k]);
		}
		k++;
	    }
        }
    }

    pputc(prn, '\n');
}

static void print_freq_test (const FreqDist *freq, PRN *prn)
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
 */

void print_freq (const FreqDist *freq, PRN *prn)
{
    int i, k, nlw, K;
    int total, valid, missing;
    char word[64];
    double f, cumf = 0;

    if (freq == NULL) {
	return;
    }

    K = freq->numbins - 1;
    valid = freq->n;
    total = freq->t2 - freq->t1 + 1;

    pprintf(prn, _("\nFrequency distribution for %s, obs %d-%d\n"),
	    freq->varname, freq->t1 + 1, freq->t2 + 1);

    if (freq->numbins == 0) {
	if (!na(freq->test)) {
	    print_freq_test(freq, prn);
	}
	return;
    } 

    if (!freq->discrete) {
	pprintf(prn, _("number of bins = %d, mean = %g, sd = %g\n"), 
		freq->numbins, freq->xbar, freq->sdx);

	pputs(prn, 
	      _("\n       interval          midpt   frequency    rel.     cum.\n\n"));
    } else {
	pputs(prn, _("\n          frequency    rel.     cum.\n\n"));
    }

    for (k=0; k<=K; k++) {
	*word = '\0';
	if (freq->discrete) {
	    sprintf(word, "%4.0f", freq->midpt[k]);
	} else {
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
	}

	pputs(prn, word);
	nlw = 10 - strlen(word);
	bufspace(nlw, prn);

	pprintf(prn, "%6d  ", freq->f[k]);
	f = 100.0 * freq->f[k] / valid;
	cumf += f;
	pprintf(prn, "  %6.2f%% %7.2f%% ", f, cumf);
	i = 0.36 * f;
	while (i--) {
	    pputc(prn, '*');
	}
	pputc(prn, '\n');
    }

    missing = total - valid;

    if (missing > 0) {
	pputc(prn, '\n');
	pprintf(prn, _("Invalid cases = %d (%5.2f%%)"), missing, 
		100 * (double) missing / total);
	pputc(prn, '\n');
    }

    if (!na(freq->test)) {
	print_freq_test(freq, prn);
    }
}

/**
 * print_xtab:
 * @tab: gretl cross-tabulation struct.
 * @prn: gretl printing struct.
 *
 * Print crosstab to @prn.
 */

void print_xtab (const Xtab *tab, gretlopt opt, PRN *prn)
{
    int i, j, r, c;
    r = tab->rows;
    c = tab->cols;
    double x, y;
    int n5 = 0;
    double pearson = 0.0;

    pprintf(prn,"\n       ");
    for (j=0; j<c; j++) {
	pprintf(prn, "[%4d]", (tab->cval)[j]);
    } 

    pprintf(prn,"  TOT.\n\n");

    for (i=0; i<r; i++) {

	if (tab->rtotal[i]) {
	    pprintf(prn,"[%4d] ",(tab->rval)[i]);

	    for (j=0; j<c; j++) {
		if (tab->ctotal[j]) {
		    if ((tab->f)[i][j] || (opt & OPT_Z)) {
			if (opt & (OPT_C | OPT_R)) {
			    if (opt & OPT_C) {
				x = 100.0 * (tab->f)[i][j] / tab->ctotal[j];
			    } else {
				x = 100.0 * (tab->f)[i][j] / tab->rtotal[i];
			    }
			    pprintf(prn, "%5.1f%%", x);
			} else {
			    pprintf(prn, "%5d ", (tab->f)[i][j]);
			}
		    } else {
			pputs(prn,"      ");
		    }
		} 

		if (!na(pearson)) {
		    y = (double) (tab->rtotal[i] * tab->ctotal[j]) / tab->n;
		    x = (double) (tab->f)[i][j] - y;
		    if (y < 1.0) {
			pearson = NADBL;
		    } else {
			pearson += x * x / y;
			if (y >= 5) n5++;
		    }
		}
	    }

	    if (opt & OPT_C) {
		x = 100.0 * tab->rtotal[i] / tab->n;
		pprintf(prn, "%5.1f%%\n", x);
	    } else {
		pprintf(prn, "%6d\n", tab->rtotal[i]);
	    }
	}
    }

    pputs(prn,_("\nTOTAL  "));

    for (j=0; j<c; j++) {
	if (opt & OPT_R) {
	    x = 100.0 * tab->ctotal[j] / tab->n;
	    pprintf(prn, "%5.1f%%", x);
	} else {
	    pprintf(prn, "%5d ", tab->ctotal[j]);
	}
    }
    
    pprintf(prn, "%6d\n", tab->n);

    if (tab->missing) {
	pprintf(prn, "\n%d missing values\n", tab->missing);
    }

    if (!na(pearson)) {
	double n5p = (double) n5 / (r * c);
	int df;

	if (n5p >= .80) {
	    df = (r - 1) * (c - 1);
	    pprintf(prn, "\nPearson chi-square test = %g (%d df, p-value = %g)\n", 
		    pearson, df, chisq(pearson, df));
	}
    }
}

/**
 * print_smpl:
 * @pdinfo: data information struct
 * @fulln: full length of data series.
 * @prn: gretl printing struct.
 *
 * Prints the current sample information to @prn.
 */

void print_smpl (const DATAINFO *pdinfo, int fulln, PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];

    if (!gretl_messages_on()) {
	return;
    }

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
    int n;

    if ((p = strstr(s, "+00")) || (p = strstr(s, "-00"))) {
	memmove(p+1, p+2, strlen(p+1));
    }

    n = strlen(s);
    if (s[n-1] == '.') {
	s[n-1] = 0;
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
   2 * P + 5 characters, where P denotes the precision ("digits"). 
*/

void gretl_print_fullwidth_double (double x, int digits, PRN *prn)
{
    char numstr[36], final[36];
    int totlen = 2 * digits + 5; /* try changing this? */
    int i, tmp, forept = 0;
    char decpoint = '.';
    char *p;

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
    } else {
	/* handle case of no decimal point, added Sept 2, 2005 */
	forept = strlen(numstr);
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

    tmp = totlen - strlen(final);
    for (i=0; i<tmp; i++) {
	strcat(final, " ");
    }

    pputs(prn, final);
}

void gretl_print_value (double x, PRN *prn)
{
    gretl_print_fullwidth_double(x, GRETL_DIGITS, prn);  
}

/**
 * print_contemporaneous_covariance_matrix:
 * @m: covariance matrix.
 * @ldet: log-determinant of @m.
 * @prn: gretl printing struct.
 * 
 * Print to @prn the covariance matrix @m, with correlations
 * above the diagonal, and followed by the log determinant.
 */

void
print_contemp_covariance_matrix (const gretl_matrix *m, 
				 double ldet, PRN *prn)
{
    int rows = gretl_matrix_rows(m);
    int cols = gretl_matrix_cols(m);
    int jmax = 1;
    char numstr[16];
    double x;
    int i, j;

    pprintf(prn, "%s\n(%s)\n\n",
	    _("Cross-equation VCV for residuals"),
	    _("correlations above the diagonal"));

    for (i=0; i<rows; i++) {
	for (j=0; j<jmax; j++) {
	    pprintf(prn, "%#13.5g", gretl_matrix_get(m, i, j));
	}
	for (j=jmax; j<cols; j++) {
	    x = gretl_matrix_get(m, i, i) * gretl_matrix_get(m, j, j);
	    x = sqrt(x);
	    x = gretl_matrix_get(m, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    pprintf(prn, "%13s", numstr);
	}
	pputc(prn, '\n');
	if (jmax < cols) {
	    jmax++;
	}
    }

    if (!na(ldet)) {
	pprintf(prn, "\n%s = %g\n", _("log determinant"), ldet);
    }
}

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
    VMatrix *vmat;
    int err = 0;

    vmat = gretl_model_get_vcv(pmod, pdinfo);

    if (vmat == NULL) {
	err = E_ALLOC;
    } else {
	text_print_vmatrix(vmat, prn);
	free_vmatrix(vmat);
    }  

    return err;
}

static void outxx (const double xx, int ci, int wid, PRN *prn)
{
    if (isnan(xx) || na(xx)) { 
	if (ci == CORR) {
	    pprintf(prn, " %*s", UTF_WIDTH(_("undefined"), wid), 
		    _("undefined"));
	} else {
	    bufspace(wid, prn);
	}
    } else if (ci == CORR) {
	pprintf(prn, " %*.4f", wid - 1, xx);
    } else {
	char numstr[18];

	if (xx > -0.001 && xx < 0.001) {
	    sprintf(numstr, "%.5e", xx);
	} else {
	    sprintf(numstr, "%g", xx);
	}
	gretl_fix_exponent(numstr);
	pprintf(prn, "%*s", wid, numstr);
    }
}

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

static int vmat_maxlen (VMatrix *vmat)
{
    int len, maxlen = 0;
    int i;

    for (i=0; i<vmat->dim; i++) {
	len = strlen(vmat->names[i]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    return maxlen;
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

void text_print_vmatrix (VMatrix *vmat, PRN *prn)
{
    register int i, j;
    int nf, li2, p, k, m, idx, ij2, lineno = 0;
    int pause = gretl_get_text_pause();
    int maxlen = 0;
    int fwidth = 14;
    int fields = 5;
    const char *s;

    if (vmat->ci != CORR) {
	covhdr(prn);
    }

    maxlen = vmat_maxlen(vmat);
    if (maxlen > 10) {
	fields = 4;
	fwidth = 16;
    }

    m = 1;

    for (i=0; i<=vmat->dim/fields; i++) {
	nf = i * fields;
	li2 = vmat->dim - nf;
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	if (pause && i > 0) {
	    takenotes(0);
	}

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    s = vmat->names[j + nf - 1];
	    bufspace(fwidth - strlen(s), prn);
	    pputs(prn, s); 
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
		idx = ijton(j, nf+k, vmat->dim);
		outxx(vmat->vec[idx], vmat->ci, fwidth, prn);
	    }
	    if (fwidth < 15) pputc(prn, ' ');
	    pprintf(prn, " %s\n", vmat->names[j]);
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
	    bufspace(fwidth * j, prn);
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, vmat->dim);
		outxx(vmat->vec[idx], vmat->ci, fwidth, prn);
	    }
	    if (fwidth < 15) pputc(prn, ' ');
	    pprintf(prn, " %s\n", vmat->names[ij2]);
	    lineno++;
	}
	pputc(prn, '\n');
    }
}

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
    pputc(prn, '\n');

    if (!na(fr->sigma)) {
	pprintf(prn, _("Standard error of residuals = %g\n"), fr->sigma);
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

/* prints names of variables in @list, positions v1 to v2 */

static void varheading (const int *list, int v1, int v2, int wid,
			const DATAINFO *pdinfo, PRN *prn)
{
    int i;

    if (csv_format(prn)) {
	pputs(prn, "Obs");
	pputc(prn, pdinfo->delim);
	for (i=v1; i<=v2; i++) { 
	    pprintf(prn, "%s", pdinfo->varname[list[i]]);
	    if (i < v2) {
		pputc(prn, pdinfo->delim);
	    } 
	}
	pputc(prn, '\n');
    } else {
	pputs(prn, "\n     Obs ");
	for (i=v1; i<=v2; i++) { 
	    pprintf(prn, "%*s", wid, pdinfo->varname[list[i]]);
	}
	pputs(prn, "\n\n");
    }
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

static void fcast_print_x (double x, int n, int pmax, PRN *prn)
{
    if (pmax != PMAX_NOT_AVAILABLE && !na(x)) {
	pprintf(prn, "%*.*f", n - 2, pmax, x);
    } else {
	gretl_printxn(x, n, prn);
    }
}

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
    } else for (t=t1; t<=t2; t++) {
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

/**
 * get_signif:
 * @x: array to examine
 * @n: length of the array
 * 
 * Examines array @x from the point of view of printing the
 * data.  Tries to determine the most economical yet faithful
 * string representation of the data.
 *
 * Returns: if successful, either a positive integer representing
 * the number of significant digits to use when printing the
 * series (e.g. when using the %%g conversion in printf), or a
 * negative integer representing the number of decimal places
 * to use (e.g. with the %%f conversion).  If unsuccessful,
 * returns #PMAX_NOT_AVAILABLE.
 */

static int get_signif (const double *x, int n)
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

	if (na(x[i])) {
	    continue;
	}

	xx = fabs(x[i]);

	if (xx > 0 && (xx < 1.0e-6 || xx > 1.0e+8)) {
	    return PMAX_NOT_AVAILABLE;
	}	

	if (xx >= 1.0) {
	    allfrac = 0;
	}

	sprintf(numstr, "%.*f", TEST_PLACES, xx);
	s = strlen(numstr) - 1;
	trail = TEST_PLACES;
	gotdec = 0;

	for (j=s; j>0; j--) {
	    if (numstr[j] == '0') {
		s--;
		if (!gotdec) {
		    trail--;
		}
	    } else if (numstr[j] == decpoint) {
		gotdec = 1;
		if (xx < 10000) {
		    break;
		} else {
		    continue;
		}
	    } else {
		break;
	    }
	}

	if (trail > trailmax) {
	    trailmax = trail;
	}

	if (xx < 1.0) {
	    s--; /* don't count leading zero */
	}

	if (s > smax) {
	    smax = s;
	}

#ifdef PRN_DEBUG
	fprintf(stderr, "get_signif: set smax = %d\n", smax);
#endif

	lead = 0;
	for (j=0; j<=s; j++) {
	    if (xx >= 1.0 && numstr[j] != decpoint) {
		lead++;
	    } else {
		break;
	    }
	}

	if (lead > leadmax) {
	    leadmax = lead;
	}
	if (lead < leadmin) {
	    leadmin = lead;
	}
    } 

    if (smax > SMAX) {
	smax = SMAX;
    }

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

static int g_too_long (double x, int signif)
{
    char n1[32], n2[32];

    sprintf(n1, "%.*G", signif, x);
    sprintf(n2, "%.0f", x);
    
    return (strlen(n1) > strlen(n2));
}

static int bufprintnum (char *buf, double x, int signif, int width)
{
    static char numstr[32];
    int i, l;

    /* guard against monster numbers that will smash the stack */
    if (fabs(x) > 1.0e20 || signif == PMAX_NOT_AVAILABLE) {
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
 * obs_marker_init:
 * @pdinfo: data information struct.
 *
 * Check the length to which observation markers should
 * be printed, in a tabular context.  (We don't want to
 * truncate 10-character date strings by chopping off the day.)
 */

static int oprintlen = 8;

void obs_marker_init (const DATAINFO *pdinfo)
{
    int t, datestrs = 0;

    if (pdinfo->markers) {
	for (t=0; t<pdinfo->n; t++) {
	    if (strlen(pdinfo->S[t]) == 10 && 
		isdigit(pdinfo->S[t][0]) &&
		strchr(pdinfo->S[t], '/')) {
		datestrs = 1;
		break;
	    }
	}
    } 

    if (datestrs) {
	oprintlen = 10;
    } else {
	oprintlen = 8;
    }
}

/**
 * print_obs_marker:
 * @t: observation number.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print a string (label, date or obs number) representing the given @t.
 */

void print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn)
{
    char tmp[OBSLEN] = {0};

    if (pdinfo->markers) { 
	strncat(tmp, pdinfo->S[t], oprintlen);
	pprintf(prn, "%*s ", oprintlen, tmp); 
    } else {
	ntodate(tmp, t, pdinfo);
	pprintf(prn, "%8s ", tmp);
    }
}

/**
 * varlist:
 * @pdinfo: data information struct.
 * @prn: gretl printing struct
 *
 * Prints a list of the names of the variables currently defined.
 */

void varlist (const DATAINFO *pdinfo, PRN *prn)
{
    int level = gretl_function_stack_depth();
    int len, maxlen = 0;
    int nv = 4;
    int i, j, n = 0;

    for (i=0; i<pdinfo->v; i++) {
	if (STACK_LEVEL(pdinfo, i) == level) {
	    len = strlen(pdinfo->varname[i]);
	    if (len > maxlen) {
		maxlen = len;
	    }
	    n++;
	}
    }

    if (maxlen < 9) {
	nv = 5;
    } else if (maxlen > 13) {
	nv = 3;
    }

    pprintf(prn, _("Listing %d variables:\n"), n);

    j = 1;
    for (i=0; i<pdinfo->v; i++) {
	if (level > 0 && STACK_LEVEL(pdinfo, i) != level) {
	    continue;
	}
	pprintf(prn, "%3d) %-*s", i, maxlen + 2, pdinfo->varname[i]);
	if (j % nv == 0) {
	    pputc(prn, '\n');
	}
	j++;
    }

    if (n % nv) {
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/**
 * maybe_list_vars:
 * @pdinfo: data information struct.
 * @prn: gretl printing struct
 *
 * Prints a list of the names of the variables currently defined,
 * unless gretl messaging is turned off.
 */

void maybe_list_vars (const DATAINFO *pdinfo, PRN *prn)
{
    if (gretl_messages_on()) {
	varlist(pdinfo, prn);
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

    pprintf(prn, "%15s = ", vname);

    if (na(x)) {
	pputs(prn, "NA");
    } else {
	if (x >= 0.0) {
	    pputc(prn, ' ');
	}
	if (opt & OPT_T) {
	    pprintf(prn, "%.10e", x);
	} else {
	    pprintf(prn, "%#.6g", x);
	}
    }

    if (allconst) {
	pputc(prn, '\n');
    }
}

static int printdata_blocks;

int get_printdata_blocks (void)
{
    return printdata_blocks;
}

/**
 * printdata:
 * @list: list of variables to print.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: if OPT_O, print the data by observation (series in columns);
 *       if OPT_N, use simple obs numbers, not dates; 
 *       if OPT_T, print the data to 10 significant digits.
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
    int j, v, v1, v2, jc, nvjc, lineno, ncol;
    int allconst, scalars = 0;
    int nvars = 0, sortvar = 0;
    int maxlen = 0, bplen = 13;
    int *plist = NULL;
    int *pmax = NULL; 
    int t, nsamp;
    char line[128];
    int err = 0;

    int pause = gretl_get_text_pause();

    printdata_blocks = 0;

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
	int len, v = plist[j];

	if (!pdinfo->vector[v]) {
	    print_scalar(Z[v][0], pdinfo->varname[v], opt, 0, prn);
	    scalars = 1;
	    gretl_list_delete_at_pos(plist, j);
	    j--;
	} else {
	    len = strlen(pdinfo->varname[v]);
	    if (len > maxlen) {
		maxlen = len;
	    }
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

    if (!(opt & OPT_O)) { 
	/* not by observations, but by variable */
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

    if (maxlen > 13) {
	ncol = 4;
	bplen = 16;
    } else {
	ncol = 5;
    }

    /* print data by observations */
    for (j=0; j<=plist[0]/ncol; j++) {
	char obs_string[OBSLEN];

	jc = j * ncol;
	nvjc = plist[0] - jc;
	v1 = jc +1;
	if (nvjc) {
	    /* starting a new block of variables */
	    v2 = (ncol > nvjc)? nvjc : ncol;
	    v2 += jc;
	    varheading(plist, v1, v2, bplen, pdinfo, prn);
	    printdata_blocks++;

	    if (pause && j > 0 && takenotes(1)) {
		goto endprint;
	    }

	    lineno = 1;

	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {

		if (sortvar && plist[0] == 1) {
		    strcpy(obs_string, SORTED_MARKER(pdinfo, sortvar, t));
		} else if (opt & OPT_N) {
		    sprintf(obs_string, "%d", t + 1);
		} else {
		    get_obs_string(obs_string, t, pdinfo);
		}
		
		sprintf(line, "%8s ", obs_string);
		
		for (v=v1; v<=v2; v++) {
		    double xx = Z[plist[v]][t];

		    if (na(xx)) {
			strcat(line, "             ");
			if (bplen == 16) {
			    strcat(line, "   ");
			}
		    } else { 
			bufprintnum(line, xx, pmax[v-1], bplen);
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
	} /* end if nvjc */
    } /* end for j loop */

    pputc(prn, '\n');

 endprint:

    free(plist);
    free(pmax);

    return err;
}

/**
 * print_data_sorted:
 * @list: list of variables to print.
 * @obsvec: list of observation numbers.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print the data for the variables in @list, using the sort order 
 * given in @obsvec.  The first element of @obsvec must contain the
 * number of observations that follow.  By default, printing is plain 
 * text, formatted in columns using space characters, but if the @prn 
 * format is set to %GRETL_FORMAT_CSV then printing respects the user's 
 * choice of column delimiter.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int print_data_sorted (const int *list, const int *obsvec, 
		       const double **Z, const DATAINFO *pdinfo, 
		       PRN *prn)
{
    char sdelim[2] = {0};
    int csv = csv_format(prn);
    int *pmax = NULL; 
    double xx;
    char obs_string[OBSLEN];
    char line[128];
    int bplen = 16;
    int T = obsvec[0];
    int i, s, t;

    /* must have a list of up to 4 variables... */
    if (list == NULL || list[0] > 4) {
	return E_DATA;
    }

    /* ...with no scalars or bad variable numbers */
    for (i=1; i<=list[0]; i++) {
	if (list[i] >= pdinfo->v || !pdinfo->vector[list[i]]) {
	    return E_DATA;
	}
    }

    /* and T must be in bounds */
    if (T > pdinfo->n - pdinfo->t1) {
	return E_DATA;
    }

    pmax = malloc(list[0] * sizeof *pmax);
    if (pmax == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
	pmax[i-1] = get_signif(Z[list[i]] + pdinfo->t1, T);
    }

    varheading(list, 1, list[0], bplen, pdinfo, prn);

    if (csv) {
	sdelim[0] = pdinfo->delim;
    }

    /* print data by observations */
    for (s=0; s<T; s++) {
	t = obsvec[s+1];
	if (t >= pdinfo->n) {
	    continue;
	}
	get_obs_string(obs_string, t, pdinfo);
	if (csv) {
	    sprintf(line, "%s", obs_string);
	    strcat(line, sdelim);
	} else {
	    sprintf(line, "%8s ", obs_string);
	}
	for (i=1; i<=list[0]; i++) {
	    xx = Z[list[i]][t];
	    if (na(xx)) {
		if (csv) {
		    strcat(line, "NA");
		} else {
		    strcat(line, "                ");
		}
	    } else { 
		if (csv) {
		    bufprintnum(line, xx, pmax[i-1], 0);
		} else {
		    bufprintnum(line, xx, pmax[i-1], bplen);
		}
	    }
	    if (csv && i < list[0]) {
		strcat(line, sdelim);
	    }
	}
	pputs(prn, line);
	pputc(prn, '\n');
    } 

    pputc(prn, '\n');
    free(pmax);

    return 0;
}

int
text_print_fit_resid (const FITRESID *fr, const DATAINFO *pdinfo, PRN *prn)
{
    int t, anyast = 0;
    double xx;

    fit_resid_head(fr, pdinfo, prn); 

    obs_marker_init(pdinfo);

    for (t=fr->t1; t<=fr->t2; t++) {
	print_obs_marker(t, pdinfo, prn);

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
    int pmax = fr->pmax;
    int errpmax = fr->pmax;
    double *maxerr = NULL;
    int time_series = (pdinfo->structure == TIME_SERIES);
    int plot = (opt & OPT_P);

    if (do_errs) {
	maxerr = malloc(fr->nobs * sizeof *maxerr);
	if (maxerr == NULL) {
	    return E_ALLOC;
	}
    }

    pputc(prn, '\n');

    if (do_errs) {
	if (fr->model_ci == ARMA || fr->model_ci == VECM) {
	    pprintf(prn, _(" For 95%% confidence intervals, z(.025) = %.2f\n"), 
		    1.96);
	} else {
	    pprintf(prn, _(" For 95%% confidence intervals, t(%d, .025) = %.3f\n"), 
		    fr->df, fr->tval);
	}
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

    if (do_errs) {
	for (t=0; t<fr->t1; t++) {
	    maxerr[t] = NADBL;
	}
	if (pmax < 4) {
	    errpmax = pmax + 1;
	}
    }

    obs_marker_init(pdinfo);

    for (t=fr->t0; t<=fr->t2; t++) {
	print_obs_marker(t, pdinfo, prn);
	fcast_print_x(fr->actual[t], 15, pmax, prn);

	if (na(fr->fitted[t])) {
	    pputc(prn, '\n');
	    continue;
	}

	fcast_print_x(fr->fitted[t], 15, pmax, prn);

	if (do_errs) {
	    if (na(fr->sderr[t])) {
		maxerr[t] = NADBL;
	    } else {
		fcast_print_x(fr->sderr[t], 15, errpmax, prn);
		maxerr[t] = fr->tval * fr->sderr[t];
		fcast_print_x(fr->fitted[t] - maxerr[t], 15, pmax, prn);
		pputs(prn, " - ");
		fcast_print_x(fr->fitted[t] + maxerr[t], 10, pmax, prn);
	    }
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    /* do we really want a plot for non-time series? */

    if (plot && fr->nobs > 0) {
	pv = plotvar(pZ, pdinfo);
	if (pv < 0) {
	    err = 1;
	} else {
	    err = plot_fcast_errs(fr->t0, fr->t2, (*pZ)[pv], 
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

#define PRINTF_DEBUG 0

#define is_format_char(c) (c == 'f' || c == 'g' || c == 'd' || c == 's')
#define numeric_conv(c)   (c == 'f' || c == 'g' || c == 'd')

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
    } else if (fc != 's') {
	/* numeric */
	if (na(val)) {
	    fmt[n-1] = 's';
	    pprintf(prn, fmt, "NA");
	} else if (fc == 'd') {
	    pprintf(prn, fmt, (int) val);
	} else {
	    pprintf(prn, fmt, val);
	}
	*pfmt += n;
    } else {
	/* string */
	if (str == NULL) {
	    fprintf(stderr, "NULL string in printf\n");
	    err = 1;
	} else {
	    pprintf(prn, fmt, str);
	    *pfmt += n;	
	}
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

static char *date_string (const char *s, double **Z, const DATAINFO *pdinfo)
{
    char darg[16];
    char *ret = NULL;
    int v, t = -1;

    if (!strncmp(s, "date(", 5) && sscanf(s + 5, "%15[^)])", darg)) {
	if (isdigit(*darg)) {
	    t = atoi(darg);
	} else {
	    v = varindex(pdinfo, darg);
	    if (v < pdinfo->v) {
		t = Z[v][0];
	    }
	}
	if (t > 0 && t <= pdinfo->n) {
	    char tstr[OBSLEN];

	    ntodate(tstr, t - 1, pdinfo);
	    ret = gretl_strdup(tstr);
	}
    }

    return ret;
}

static char *varname_string (const char *s, const DATAINFO *pdinfo)
{
    char *vname = NULL;
    int v;

    if (sscanf(s, "varname(%d)", &v)) {
	if (v >= 0 && v < pdinfo->v) {
	    vname = gretl_strdup(pdinfo->varname[v]);
	}
    }

    return vname;
}

static char *literal_string (const char *s)
{
    char *ret = NULL;
    int n = strlen(s);

    if (*s == '"' && s[n - 1] == '"') {
	ret = gretl_strndup(s + 1, n - 2);
    }

    return ret;
}

static int get_conversion (const char *s, int *skip)
{
    *skip = strspn(s, "#0123456789.");
    return (*(s + *skip));
}

static int real_do_printf (const char *line, double ***pZ, 
			   DATAINFO *pdinfo, PRN *prn, int t)
{
    const char *p;
    char format[128];
    char *argv, *str = NULL;
    double *xvals = NULL;
    char **svals = NULL;
    int argc = 0, inparen = 0;
    int xcnv = 0, scnv = 0;
    int markerpos = -1;
    int markeroffset = 0;
    int i, err = 0;

    if (t < 0) {
	t = pdinfo->t1;
    }

#if PRINTF_DEBUG
    fprintf(stderr, "do_printf: line='%s'\n", line);
#endif

    *gretl_errmsg = '\0';

    if (!strncmp(line, "printf ", 7)) {
	line += 7;
    }

    if (sscanf(line, "\"%127[^\"]\"", format) != 1) {
	return 1;
    }

#if PRINTF_DEBUG
    fprintf(stderr, "do_printf: format='%s'\n", format);
#endif

    p = format;
    while (*p) {
	int c, skip;

	if (*p == '%') {
	    c = get_conversion(p + 1, &skip);
	    if (c == '%') {
		p++;
	    } else if (numeric_conv(c)) {
		xcnv++;
		p += skip;
	    } else if (c == 's') {
		scnv++;
		p += skip;
	    } else {
		err = 1;
		break;
	    }
	}
	p++;
    }

    if (err) {
	return err;
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
    if (argc != xcnv + scnv) {
	fprintf(stderr, "do_printf: argc = %d but conversions = %d\n",
		argc, xcnv + scnv);
	err = 1;
	goto printf_bailout;
    }

    /* play safe with sizes here */
    xvals = malloc(argc * sizeof *xvals);
    str = malloc(strlen(line) + 1);
    svals = create_strings_array(argc);

    if (xvals == NULL || svals == NULL || str == NULL) {
	err = E_ALLOC;
	goto printf_bailout;
    }

    strcpy(str, line);

    for (i=0; i<argc; i++) {
	char *special;

	argv = get_arg((i > 0)? NULL : str);
	chopstr(argv);

	xvals[i] = NADBL;
	svals[i] = NULL;

#if PRINTF_DEBUG
	fprintf(stderr, "do_printf: processing argv[%d] '%s'\n", i, argv);	
#endif
	if (numeric_string(argv)) {
	    xvals[i] = atof(argv);
	} else if (!strncmp(argv, "marker", 6)) {
	    if (markerpos >= 0 || pdinfo->S == NULL) {
		err = 1;
	    } else {
		markerpos = i;
		xvals[i] = 0.0;
		markeroffset = get_marker_offset(argv);
	    }
	} else if ((special = varname_string(argv, pdinfo)) != NULL) {
	    svals[i] = special;
	} else if ((special = date_string(argv, *pZ, pdinfo)) != NULL) {
	    svals[i] = special;
	} else if ((special = literal_string(argv)) != NULL) {
	    svals[i] = special;
	} else {
	    int v = varindex(pdinfo, argv);

	    if (v < pdinfo->v) {
#if PRINTF_DEBUG
		fprintf(stderr, "'%s' is variable #%d (vector = %d)\n",
			argv, v, pdinfo->vector[v]);
#endif
		/* simple existent varname */
		if (pdinfo->vector[v]) {
		    xvals[i] = (*pZ)[v][t];
		} else {
		    xvals[i] = (*pZ)[v][0];
		}
	    } else {
#if PRINTF_DEBUG
		fprintf(stderr, "looked up '%s' as variable: not found\n", argv);
#endif
		err = get_generated_value(argv, &xvals[i], pZ, pdinfo, t);
	    }
	}

#if PRINTF_DEBUG
	fprintf(stderr, " after processing arg, xvals[%d] = %g, err = %d\n", 
		i, xvals[i], err);	
#endif
	if (err) {
	    goto printf_bailout;
	}
    }    

    p = format;
    i = 0;
    while (*p && !err) {
	const char *s = NULL;

	if (*p == '%') {
	    if (*(p + 1) == '%') {
		pputc(prn, '%');
		p += 2;
	    } else {
		if (i == markerpos) {
		    s = pdinfo->S[t];
		    if (markeroffset > 0 && 
			markeroffset < strlen(pdinfo->S[t])) {
			s += markeroffset;
		    }
		} else if (svals[i] != NULL) {
		    s = svals[i];
		}
		err = print_arg(&p, xvals[i], s, prn);
		i++;
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

    free(xvals);
    free(str);
    free_strings_array(svals, argc);

    return err;
}

int do_printf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, PRN *prn)
{
    return real_do_printf(line, pZ, pdinfo, prn, -1);
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
	    err = real_do_printf(s, pZ, pdinfo, prn, t);
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
 * This function (which works rather like fgets) must be initialized 
 * via the call: bufgets_init(buf);
 * 
 * Returns: @s (or %NULL if nothing more can be read from @buf).
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
	return NULL;
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

/**
 * bufgets_init:
 * @buf: source buffer.
 *
 * Initializes a text buffer for use with bufgets(). 
 */

void bufgets_init (const char *buf)
{
    bufgets(NULL, 0, buf);
}
