/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/*  printout.c - simple text print routines for some gretl structs */ 

#include "libgretl.h"
#include "version.h"
#include "libset.h"
#include "forecast.h"
#include "gretl_func.h"
#include "matrix_extra.h"
#include "gretl_scalar.h"

#include <time.h>

#define PRN_DEBUG 0

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
	myprn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
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
    puts(_("Copyright Ramu Ramanathan, Allin Cottrell and Riccardo \"Jack\" Lucchetti"));
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
	myprn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
	prn = myprn;
    }
	
    pprintf(prn, _("gretl: gui client for gretl version %s,\n"), GRETL_VERSION);
    pputs(prn, _("Copyright Allin Cottrell and Riccardo \"Jack\" Lucchetti"));
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
	pprintf(prn, " %#12.6g %#12.6g", 
		cf->coeff[i] - cf->maxerr[i],
		cf->coeff[i] + cf->maxerr[i]);
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
    double tail = cf->alpha / 2;
    int i;

    if (cf->asy) {
	pprintf(prn, "z(%g) = %.4f\n\n", tail, cf->t);
    } else {
	pprintf(prn, "t(%d, %g) = %.3f\n\n", cf->df, tail, cf->t);
    }

    /* xgettext:no-c-format */
    pprintf(prn, _("      VARIABLE         COEFFICIENT      %g%% CONFIDENCE "
		   "INTERVAL\n\n"), 100 * (1 - cf->alpha));      

    for (i=0; i<cf->ncoeff; i++) {
	print_coeff_interval(cf, i, prn);
    }

    pputc(prn, '\n');
}

static void print_freq_test (const FreqDist *freq, PRN *prn)
{
    double pval = NADBL;

    if (freq->dist == D_NORMAL) {
	pval = chisq_cdf_comp(2, freq->test);
	pprintf(prn, "\n%s:\n", 
		_("Test for null hypothesis of normal distribution"));
	pprintf(prn, "%s(2) = %.3f %s %.5f\n", 
		_("Chi-square"), freq->test, 
		_("with p-value"), pval);
    } else if (freq->dist == D_GAMMA) {
	pval = normal_pvalue_2(freq->test);
	pprintf(prn, "\n%s:\n", 
		_("Test for null hypothesis of gamma distribution"));
	pprintf(prn, "z = %.3f %s %.5f\n", freq->test, 
		_("with p-value"), pval);
    }	

    if (!na(pval)) {
	record_test_result(freq->test, pval, 
			   (freq->dist == D_NORMAL)? 
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

    if (K == 0) {
	pprintf(prn, "\n constant: value = %g\n", freq->midpt[0]);
	return;
    }

    if (freq->discrete) {
	pputs(prn, _("\n          frequency    rel.     cum.\n\n"));

	for (k=0; k<=K; k++) {
	    sprintf(word, "%4g", freq->midpt[k]);
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
    } else {
	int digits = 5;
	int someneg = 0, somemneg = 0;
	int len, xlen, mxlen;
	double x;

	pprintf(prn, _("number of bins = %d, mean = %g, sd = %g\n"), 
		freq->numbins, freq->xbar, freq->sdx);
	pputs(prn, 
	      _("\n       interval          midpt   frequency    rel.     cum.\n\n"));

    tryagain:

	xlen = mxlen = 0;

	for (k=0; k<=K; k++) {
	    x = freq->endpt[k];
	    if (x < 0) {
		someneg = 1;
	    }
	    sprintf(word, "%#.*g", digits, x);
	    len = strlen(word);
	    if (len > xlen) {
		xlen = len;
	    }
	    x = freq->midpt[k];
	    if (x < 0) {
		somemneg = 1;
	    }
	    sprintf(word, "%#.*g", digits, x);
	    len = strlen(word);
	    if (len > mxlen) {
		mxlen = len;
	    }
	}

	if (xlen > 10 && digits == 5) {
	    digits--;
	    goto tryagain;
	}

	xlen++;
	xlen = (xlen > 10)? xlen : 10;

	mxlen++;
	mxlen = (mxlen > 10)? mxlen : 10;
	
	for (k=0; k<=K; k++) {
	    *word = '\0';
	    if (k == 0) {
		pprintf(prn, "%*s", xlen + 3, " < ");
	    } else if (k == K) {
		pprintf(prn, "%*s", xlen + 3, ">= ");
	    } else {
		sprintf(word, "%#.*g", digits, freq->endpt[k]);
		pprintf(prn, "%*s", xlen, word);
		pputs(prn, " - ");
	    }

	    x = (k == K)? freq->endpt[k] : freq->endpt[k+1];
	    if (x > 0 && someneg) {
		sprintf(word, " %#.*g", digits, x);
	    } else {
		sprintf(word, "%#.*g", digits, x);
	    }
	    pprintf(prn, "%-*s", xlen, word);

	    x = freq->midpt[k];
	    if (x > 0 && somemneg) {
		sprintf(word, " %#.*g", digits, x);
	    } else {
		sprintf(word, "%#.*g", digits, x);
	    }
	    pprintf(prn, "%-*s", mxlen, word);

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
    }

    missing = total - valid;

    if (missing > 0) {
	pprintf(prn, "\n%s = %d (%5.2f%%)\n", _("Missing observations"), 
		missing, 100 * (double) missing / total);
    }

    if (!na(freq->test)) {
	print_freq_test(freq, prn);
    } else {
	pputc(prn, '\n');
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
    int r = tab->rows;
    int c = tab->cols;
    double x, y;
    int n5 = 0;
    double ymin = 1.0e-7;
    double pearson = 0.0;
    int i, j;

    pputc(prn, '\n');
    pprintf(prn, _("Cross-tabulation of %s (rows) against %s (columns)"),
	    tab->rvarname, tab->cvarname);

    pputs(prn, "\n\n       ");
    for (j=0; j<c; j++) {
	pprintf(prn, "[%4g]", tab->cval[j]);
    } 

    pprintf(prn,"  %s\n  \n", _("TOT."));

    for (i=0; i<r; i++) {

	if (tab->rtotal[i] > 0) {
	    pprintf(prn, "[%4g] ", tab->rval[i]);

	    for (j=0; j<c; j++) {
		if (tab->ctotal[j]) {
		    if (tab->f[i][j] || (opt & OPT_Z)) {
			if (opt & (OPT_C | OPT_R)) {
			    if (opt & OPT_C) {
				x = 100.0 * tab->f[i][j] / tab->ctotal[j];
			    } else {
				x = 100.0 * tab->f[i][j] / tab->rtotal[i];
			    }
			    pprintf(prn, "%5.1f%%", x);
			} else {
			    pprintf(prn, "%5d ", tab->f[i][j]);
			}
		    } else {
			pputs(prn, "      ");
		    }
		} 

		if (!na(pearson)) {
		    y = (double) (tab->rtotal[i] * tab->ctotal[j]) / tab->n;
		    x = (double) tab->f[i][j] - y;
		    if (y < ymin) {
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

    pputc(prn, '\n');
    pputs(prn, _("TOTAL  "));

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
	pputc(prn, '\n');
	pprintf(prn, _("%d missing values"), tab->missing);
	pputc(prn, '\n');
    }

    if (na(pearson)) {
	pprintf(prn, _("Pearson chi-square test not computed: some "
		       "expected frequencies were less\n"
		       "than %g\n"), ymin);
    } else {
	double n5p = (double) n5 / (r * c);
	int df = (r - 1) * (c - 1);

	pputc(prn, '\n');
	pprintf(prn, _("Pearson chi-square test = %g (%d df, p-value = %g)"), 
		pearson, df, chisq_cdf_comp(df, pearson));
	pputc(prn, '\n');
	if (n5p < 0.80) {
	    pputs(prn, "Warning: Less than of 80% of cells had expected "
		  "values of 5 or greater.\n");
	}
    }

    if (r == 2 && c == 2) {
	fishers_exact_test(tab, prn);
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
    if (!gretl_messages_on() || pdinfo->v == 0 || gretl_looping_quietly()) {
	return;
    }

    if (fulln && !dataset_is_panel(pdinfo)) {
	pprintf(prn, _("Full data set: %d observations\n"),
		fulln);
	pprintf(prn, _("Current sample: %d observations\n"),
		pdinfo->n);
	return;
    }

    if (fulln) {
	pprintf(prn, _("Full data set: %d observations\n"), fulln);
    } else {
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Full data range"), 
		pdinfo->stobs, pdinfo->endobs, pdinfo->n);
    }

    if (pdinfo->t1 > 0 || pdinfo->t2 < pdinfo->n - 1 ||
	(fulln && dataset_is_panel(pdinfo))) {
	char d1[OBSLEN], d2[OBSLEN];
	ntodate_full(d1, pdinfo->t1, pdinfo);
	ntodate_full(d2, pdinfo->t2, pdinfo);

	pprintf(prn, "%s:  %s - %s", _("Current sample"), d1, d2);
	pprintf(prn, " (n = %d)\n", pdinfo->t2 - pdinfo->t1 + 1);
    }

    pputc(prn, '\n');
}

static void print_var_smpl (int v, const double **Z, 
			    const DATAINFO *pdinfo, 
			    PRN *prn)
{
    int t, n = 0;

    if (pdinfo->t1 > 0 || pdinfo->t2 < pdinfo->n - 1) {
	char d1[OBSLEN], d2[OBSLEN];
	ntodate_full(d1, pdinfo->t1, pdinfo);
	ntodate_full(d2, pdinfo->t2, pdinfo);

	pprintf(prn, "%s:  %s - %s", _("Current sample"), d1, d2);
    } else {
	pprintf(prn, "%s: %s - %s", _("Full data range"), 
		pdinfo->stobs, pdinfo->endobs);
    }

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (!na(Z[v][t])) {
	    n++;
	}
    }

    pprintf(prn, " (n = %d)\n", n);
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
	if (*(p + 3)) {
	    memmove(p+1, p+2, strlen(p+1));
	}
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
    int tex = tex_format(prn);
    int rows = gretl_matrix_rows(m);
    int cols = gretl_matrix_cols(m);
    int jmax = 1;
    char numstr[16];
    double x;
    int i, j;

    if (tex) {
	pputs(prn, "\\begin{center}\n");
	pprintf(prn, "%s \\\\\n", I_("Cross-equation VCV for residuals"));
	pprintf(prn, "(%s)\n\n", I_("correlations above the diagonal"));
	pputs(prn, "\\[\n\\begin{array}{");
	for (j=0; j<cols; j++) {
	    pputc(prn, 'c');
	}
	pputs(prn, "}\n");
    } else {
	pprintf(prn, "%s\n", _("Cross-equation VCV for residuals"));
	pprintf(prn, "(%s)\n\n", _("correlations above the diagonal"));
    }

    for (i=0; i<rows; i++) {
	for (j=0; j<jmax; j++) {
	    pprintf(prn, "%#13.5g", gretl_matrix_get(m, i, j));
	    if (tex && j < cols - 1) {
		pputs(prn, " & ");
	    }
	}
	for (j=jmax; j<cols; j++) {
	    x = gretl_matrix_get(m, i, i) * gretl_matrix_get(m, j, j);
	    x = sqrt(x);
	    x = gretl_matrix_get(m, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    pprintf(prn, "%13s", numstr);
	    if (tex && j < cols - 1) {
		pputs(prn, " & ");
	    }	    
	}
	if (tex) {
	    pputs(prn, "\\\\\n");
	} else {
	    pputc(prn, '\n');
	}
	if (jmax < cols) {
	    jmax++;
	}
    }

    if (tex) {
	pputs(prn, "\\end{array}\n\\]\n");
    }    

    if (!na(ldet)) {
	if (tex) {
	    if (ldet < 0) {
		pprintf(prn, "\n%s = ", I_("log determinant"));
		pprintf(prn, "$-$%g\n", -ldet);
	    } else {
		pprintf(prn, "\n%s = %g\n", I_("log determinant"), ldet);
	    }
	} else {
	    pprintf(prn, "\n%s = %g\n", _("log determinant"), ldet);
	}
    }

    if (tex) {
	pputs(prn, "\n\\end{center}\n");
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

static void outxx (double x, int ci, int wid, PRN *prn)
{
    if (isnan(x) || na(x)) { 
	if (ci == CORR) {
	    pprintf(prn, "%*s", UTF_WIDTH(_("NA"), wid), 
		    _("NA"));
	} else {
	    bufspace(wid, prn);
	}
    } else if (ci == CORR) {
	pprintf(prn, " %*.4f", wid - 1, x);
    } else {
	char numstr[18];

	if (x == -0) x = 0.0;

	if (x != 0 && x > -0.001 && x < 0.001) {
	    sprintf(numstr, "%.5e", x);
	} else {
	    sprintf(numstr, "%g", x);
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
    char label[16];
    char obs1[OBSLEN], obs2[OBSLEN];
    int onestep = fr->method == FC_ONESTEP;
    int i;

    if (onestep) {
	ntodate(obs1, fr->model_t1, pdinfo);   
	pputs(prn, _("Recursive one-step ahead forecasts"));
	pputs(prn, "\n\n");
	pprintf(prn, _("The forecast for time t is based on (a) coefficients obtained by\n"
		       "estimating the model over the sample %s to t-1, and (b) the\n"
		       "regressors evaluated at time t."), obs1);
	pputs(prn, "\n\n");
	pputs(prn, _("This is truly a forecast only if all the stochastic regressors\n"
		     "are in fact lagged values."));
	pputs(prn, "\n\n");
    } else {
	ntodate(obs1, fr->t1, pdinfo);
	ntodate(obs2, fr->t2, pdinfo);
	pprintf(prn, _("Model estimation range: %s - %s"), obs1, obs2);
	pputc(prn, '\n');

	if (!na(fr->sigma)) {
	    pprintf(prn, _("Standard error of residuals = %g\n"), fr->sigma);
	}
    }
    
    pprintf(prn, "\n     %s ", _("Obs"));

    for (i=1; i<4; i++) {
	if (i == 1) strcpy(label, fr->depvar);
	if (i == 2) strcpy(label, (onestep)? _("forecast") : _("fitted"));
	if (i == 3) strcpy(label, (onestep)? _("error") : _("residual"));
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

static void printstr_long (PRN *prn, double xx, int d, int *ls)
{
    char str[64];
    int n;

    if (na(xx)) {
	strcpy(str, "NA  ");
	n = 4;
    } else {
	sprintf(str, "%#.*E  ", d, xx);
	n = strlen(str);
    }
    if (*ls + n > 78) {
	pputc(prn, '\n');
	*ls = 0;
    }
    pputs(prn, str);
    *ls += n;
}

/* prints series z from current sample t1 to t2 */

static void printz (const double *z, const DATAINFO *pdinfo, 
		    PRN *prn, gretlopt opt)
{
    int t, dig = 10, ls = 0;
    double xx;

    if (opt & OPT_L) {
	dig = libset_get_int(LONGDIGITS);
	opt = OPT_T;
    }

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xx = z[t];
	if (opt & OPT_T) {
	    printstr_long(prn, xx, dig, &ls);
	} else {
	    char str[32];
	    int n;

	    if (na(xx)) {
		strcpy(str, "NA  ");
		n = 4;
	    } else {
		sprintf(str, "%#.*g  ", GRETL_DIGITS, xx);
		n = strlen(str);
	    }
	    if (ls + n > 78) {
		pputc(prn, '\n');
		ls = 0;
	    }
	    pputs(prn, str);
	    ls += n;
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

#if PRN_DEBUG
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
#if PRN_DEBUG
	fprintf(stderr, "get_signif: setting smax = -(%d - %d)\n", 
		smax, leadmax);
#endif	
	smax = -1 * (smax - leadmax); /* # of decimal places */
    } else if (leadmax == smax) {
	smax = 0;
    } else if (leadmax == 0 && !allfrac) {
#if PRN_DEBUG
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
#if PRN_DEBUG
	    fprintf(stderr, "got %d for signif: "
		    "printing with %%.%df\n", signif, -signif);
#endif
	sprintf(numstr, "%.*f", -signif, x);
    } else if (signif == 0) {
#if PRN_DEBUG
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
#if PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%dG\n", l, signif, signif);
#endif
	    if (g_too_long(x, signif)) {
		sprintf(numstr, "%.0f", x);
	    } else {
		sprintf(numstr, "%.*G", signif, x);
	    }
	} else if (z >= .10) {
#if PRN_DEBUG
	    fprintf(stderr, "got %d for leftvals, %d for signif: "
		    "printing with %%.%df\n", l, signif, signif-l);
#endif
	    sprintf(numstr, "%.*f", signif - l, x);
	} else {
	    if (signif > 4) signif = 4;
#if PRN_DEBUG
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
    int level = gretl_function_depth();
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

static void print_varlist (const char *name, const int *list, 
			   const DATAINFO *pdinfo, PRN *prn)
{
    int i, v, len = 0;

    if (list[0] == 0) {
	pprintf(prn, " %s\n", _("list is empty"));
    } else {
	len += pprintf(prn, " %s: ", name);
	for (i=1; i<=list[0]; i++) {
	    v = list[i];
	    if (v == LISTSEP) {
		len += pputs(prn, "; ");
	    } else if (v >= 0 && v < pdinfo->v) {
		len += pprintf(prn, "%s ", pdinfo->varname[v]);
	    } else {
		len += pprintf(prn, "%d ", v);
	    }
	    if (i < list[0] && len > 68) {
		pputs(prn, " \\\n ");
		len = 0;
	    }
	}
	pputc(prn, '\n');
    }
}

static void print_scalar (double x, const char *vname, 
			  gretlopt opt, PRN *prn)
{
    pputc(prn, '\n');
    pprintf(prn, "%15s = ", vname);

    if (na(x)) {
	pputs(prn, "NA");
    } else {
	if (x >= 0.0) {
	    pputc(prn, ' ');
	}
	if (opt & OPT_L) {
	    pprintf(prn, "%#.*E", libset_get_int(LONGDIGITS), x);
	} else if (opt & OPT_T) {
	    pprintf(prn, "%#.10E", x);
	} else {
	    pprintf(prn, "%#.6g", x);
	}
    }

    pputc(prn, '\n');
}

static void 
print_listed_objects (const char *s, const DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    const gretl_matrix *m;
    const int *list;
    char *name;

    while ((name = gretl_word_strdup(s, &s)) != NULL) {
	if (gretl_is_scalar(name)) {
	    print_scalar(gretl_scalar_get_value(name), 
			 name, opt, prn);
	} else if ((m = get_matrix_by_name(name)) != NULL) {
	    gretl_matrix_print_to_prn(m, name, prn);
	} else if ((list = get_list_by_name(name)) != NULL) {
	    print_varlist(name, list, pdinfo, prn);
	}
	free(name);
    }
}

static int printdata_blocks;

int get_printdata_blocks (void)
{
    return printdata_blocks;
}

static int adjust_print_list (int *list, int *screenvar,
			      gretlopt opt)
{
    int pos;

    if (!(opt & OPT_O)) {
	return E_PARSE;
    }

    pos = gretl_list_separator_position(list);

    if (list[0] < 3 || pos != list[0] - 1) {
	return E_PARSE;
    } else {
	*screenvar = list[list[0]];
	list[0] = pos - 1;
    }

    return 0;
}

/**
 * printdata:
 * @list: list of variables to print.
 * @mstr: optional string holding names of matrices to print.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * @opt: if %OPT_O, print the data by observation (series in columns);
 * if %OPT_N, use simple obs numbers, not dates; if %OPT_T, print the 
 * data to 10 significant digits; if %OPT_L, print the data to a number 
 * of digits set by "set longdigits" (default 10).
 * @prn: gretl printing struct.
 *
 * Print the data for the variables in @list, from observations t1 to
 * t2.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int printdata (const int *list, const char *mstr, 
	       const double **Z, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    int pause = gretl_get_text_pause();
    int j, v, v1, v2, jc, nvjc, lineno, ncol;
    int screenvar = 0;
    int sortvar = 0;
    int maxlen = 0, bplen = 13;
    int *plist = NULL;
    int *pmax = NULL; 
    int t, nsamp;
    char line[128];
    int err = 0;

    printdata_blocks = 0;

    if (list != NULL && list[0] == 0) {
	/* explicitly empty list given */
	if (mstr == NULL) {
	    return 0; /* no-op */
	} else {
	    goto endprint;
	}
    } else if (list == NULL) {
	/* no list given */
	if (mstr == NULL) {
	    int nvars = 0;

	    plist = full_var_list(pdinfo, &nvars);
	    if (nvars == 0) {
		/* no-op */
		return 0;
	    }
	} else {
	    goto endprint;
	} 
    } else {
	plist = gretl_list_copy(list);
    }

    /* at this point plist should have something in it */
    if (plist == NULL) {
	return E_ALLOC;
    } 

    if (gretl_list_has_separator(plist)) {
	err = adjust_print_list(plist, &screenvar, opt);
	if (err) {
	    free(plist);
	    return err;
	}
    }

    lineno = 1;

    if (plist[0] == 0) {
	/* no series */
	pputc(prn, '\n');
	goto endprint;
    }

    if (gretl_print_has_buffer(prn)) {
	/* how big a job do we have? */
	int T = pdinfo->t2 - pdinfo->t1 + 1;
	int nx = plist[0] * T;

	if (nx > 1000) {
	    err = gretl_print_alloc(prn, nx * 12);
	}
    }

    if (!(opt & OPT_O)) { 
	/* not by observations, but by variable */
	pputc(prn, '\n');
	for (j=1; j<=plist[0]; j++) {
	    int vj = plist[j];

	    if (plist[0] > 1) {
		pprintf(prn, _("Varname: %s\n"), pdinfo->varname[vj]);
	    }
	    print_var_smpl(vj, Z, pdinfo, prn);
	    pputc(prn, '\n');
	    printz(Z[vj], pdinfo, prn, opt);
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

    /* main block: print data by observations */

    for (j=0; j<=plist[0]/ncol; j++) {
	char obs_string[OBSLEN];

	jc = j * ncol;
	nvjc = plist[0] - jc;
	v1 = jc + 1;

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

		if (screenvar && Z[screenvar][t] == 0.0) {
		    /* screened out by boolean */
		    continue;
		}

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

    if (mstr != NULL) {
	print_listed_objects(mstr, pdinfo, opt, prn);
    }

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

    /* ...with no bad variable numbers */
    for (i=1; i<=list[0]; i++) {
	if (list[i] >= pdinfo->v) {
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

#define SIGMA_MIN 1.0e-18

int
text_print_fit_resid (const FITRESID *fr, const DATAINFO *pdinfo, 
		      PRN *prn)
{
    int onestep = fr->method == FC_ONESTEP;
    int t, anyast = 0;
    double yt, yf, et;
    double MSE = 0.0;
    double AE = 0.0;
    int effn = 0;
    int err = 0;

    fit_resid_head(fr, pdinfo, prn); 

    obs_marker_init(pdinfo);

    for (t=fr->t1; t<=fr->t2; t++) {
	print_obs_marker(t, pdinfo, prn);

	yt = fr->actual[t];
	yf = fr->fitted[t];
	et = fr->resid[t];

	if (na(yt)) {
	    pputc(prn, '\n');
	} else if (na(yf)) {
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13.*f\n", fr->pmax, yt);
	    } else {
		pprintf(prn, "%13g\n", yt);
	    }
	} else if (na(et)) {
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13.*f%13.*f\n", fr->pmax, yt, yf);
	    } else {
		pprintf(prn, "%13g%13g\n", yt, yf);
	    }
	} else {
	    int ast = 0;

	    if (onestep) {
		MSE += et * et;
		AE += fabs(et);
		effn++;
	    } else if (fr->sigma > SIGMA_MIN) {
		ast = (fabs(et) > 2.5 * fr->sigma);
		if (ast) {
		    anyast = 1;
		}
	    }

	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13.*f%13.*f%13.*f%s\n", 
			fr->pmax, yt, fr->pmax, yf, fr->pmax, et,
			(ast)? " *" : "");
	    } else {
		pprintf(prn, "%13g%13g%13g%s\n", 
			yt, yf, et,
			(ast)? " *" : "");
	    }
	}
    }

    pputc(prn, '\n');

    if (anyast) {
	pputs(prn, _("Note: * denotes a residual in excess of "
		     "2.5 standard errors\n"));
    }

    if (effn > 0) {
	MSE /= effn;
	pprintf(prn, "%s = %g\n", _("Mean Squared Error"), MSE);
	pprintf(prn, "%s = %g\n", _("Root Mean Squared Error"), sqrt(MSE));
	pprintf(prn, "%s = %g\n", _("Mean Absolute Error"), AE / effn);
    }

    if (onestep && fr->nobs > 0 && gretl_in_gui_mode()) {
	const double *obs = gretl_plotx(pdinfo);
	int ts = dataset_is_time_series(pdinfo);
	int t0 = (fr->t0 >= 0)? fr->t0 : 0;

	if (obs == NULL) {
	    err = 1;
	} else {
	    err = plot_fcast_errs(t0, fr->t2, obs, 
				  fr->actual, fr->fitted, NULL, 
				  fr->depvar, (ts)? pdinfo->pd : 0,
				  OPT_NONE);
	}
    }


    return err;
}

/**
 * text_print_forecast:
 * @fr: pointer to structure containing forecasts.
 * @pdinfo: dataset information.
 * @opt: if includes %OPT_P, make a plot of the forecasts.
 * @prn: printing structure.
 *
 * Print the forecasts in @fr to @prn, and also plot the
 * forecasts if %OPT_P is given.  If %OPT_F is given in
 * addition to %OPT_P, use fill style for the confidence
 * bands when plotting a forecast with standard errors.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int text_print_forecast (const FITRESID *fr, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn)
{
    int do_errs = (fr->sderr != NULL);
    int pmax = fr->pmax;
    int errpmax = fr->pmax;
    double *maxerr = NULL;
    int t, err = 0;

    if (do_errs) {
	maxerr = malloc(fr->nobs * sizeof *maxerr);
	if (maxerr == NULL) {
	    return E_ALLOC;
	}
    }

    if (!(opt & OPT_Q)) {
	pputc(prn, '\n');
    }

    if (do_errs && !(opt & OPT_Q)) {
	if (fr->asymp) {
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

    if ((opt & OPT_P) && fr->nobs > 0) {
	const double *obs = gretl_plotx(pdinfo);
	int ts = dataset_is_time_series(pdinfo);

	/* FIXME daily data x-axis */

	if (obs == NULL) {
	    err = 1;
	} else {
	    err = plot_fcast_errs(fr->t0, fr->t2, obs, 
				  fr->actual, fr->fitted, maxerr, 
				  fr->depvar, (ts)? pdinfo->pd : 0,
				  opt);
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
    int err = 0;

    fr = get_fit_resid(pmod, Z, pdinfo, &err);

    if (!err) {
	text_print_fit_resid(fr, pdinfo, prn);
	free_fit_resid(fr);
    }

    return err;
}

static void print_iter_val (double x, int i, int k, PRN *prn)
{
    if (na(x)) {
	pprintf(prn, "%-12s", "NA");
    } else {
	pprintf(prn, "%#12.5g", x);
    }
    if (i && i % 6 == 5 && i < k-1) {
	pprintf(prn, "\n%12s", " ");
    }
}

/**
 * print_iter_info:
 * @iter: iteration number.
 * @crit: criterion (e.g. log-likelihood).
 * @type: type of criterion (%C_LOGLIK or %C_OTHER)
 * @k: number of parameters.
 * @b: parameter array.
 * @g: gradient array.
 * @sl: step length.
 * @prn: gretl printing struct.
 *
 * Print to @prn information pertaining to step @iter of an 
 * iterative estimation process.
 */

void 
print_iter_info (int iter, double crit, int type, int k, 
		 const double *b, const double *g, 
		 double sl, PRN *prn)
{
    const char *cstrs[] = {
	N_("Log-likelihood"),
	N_("GMM criterion"),
	N_("Criterion"),
    };
    const char *cstr = cstrs[type];
    int i;

    if (type == C_GMM) {
	crit = -crit;
    }

    if (iter < 0) {
	pputs(prn, _("\n--- FINAL VALUES: \n"));
    } else {
	pprintf(prn, "%s %d: ", _("Iteration"), iter);
    }

    if (na(crit) || na(-crit)) {
	pprintf(prn, "%s = NA", _(cstr));
    } else {
	pprintf(prn, "%s = %#.12g", _(cstr), crit);
    }

    if (sl > 0.0 && !na(sl)) {
	pprintf(prn, _(" (steplength = %g)"), sl);
    }	

    pputc(prn, '\n');
	
    pputs(prn, _("Parameters: "));
    for (i=0; i<k; i++) {
	print_iter_val(b[i], i, k, prn);
    }
    pputc(prn, '\n');

    pputs(prn, _("Gradients:  "));
    for (i=0; i<k; i++) {
	print_iter_val(g[i], i, k, prn);
    }
    pputs(prn, "\n\n");

    if (iter < 0 || (iter % 20 == 0)) {
	/* experimental */
	iter_print_callback((iter < 0)? 0 : iter, prn);
    }
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

struct readbuf {
    const char *start;
    const char *point;
};

static struct readbuf *rbuf;
static int n_bufs;

static int rbuf_push (const char *s)
{
    struct readbuf *tmp = NULL;
    int i;

    for (i=0; i<n_bufs; i++) {
	if (rbuf[i].start == NULL) {
	    /* re-use existing slot */
	    rbuf[i].start = s;
	    rbuf[i].point = s;
	    return 0;
	}
    }    

    tmp = realloc(rbuf, (n_bufs + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    rbuf = tmp;
    rbuf[n_bufs].start = s;
    rbuf[n_bufs].point = s;
    n_bufs++;

    return 0;
}

static const char *rbuf_get_point (const char *s)
{
    int i;

    for (i=0; i<n_bufs; i++) {
	if (rbuf[i].start == s) {
	    return rbuf[i].point;
	}
    }

    return NULL;
}

static void rbuf_set_point (const char *s, const char *p)
{
    int i;

    for (i=0; i<n_bufs; i++) {
	if (rbuf[i].start == s) {
	    rbuf[i].point = p;
	    break;
	}
    }
}

static void rbuf_finalize (const char *s)
{
    int i;

    for (i=0; i<n_bufs; i++) {
	if (rbuf[i].start == s) {
	    rbuf[i].start = NULL;
	    rbuf[i].point = NULL;
	    break;
	}
    }    
}

/**
 * bufgets:
 * @s: target string (must be pre-allocated).
 * @size: maximum number of characters to read.
 * @buf: source buffer.
 *
 * This function works much like fgets, reading successive lines 
 * from a buffer rather than a file.  
 * Important note: use of bufgets() on a particular buffer must be 
 * preceded by a call to bufgets_init() on the same buffer, and must be
 * followed by a call to bufgets_finalize(), again on the same
 * buffer.
 * 
 * Returns: @s (or %NULL if nothing more can be read from @buf).
 */

char *bufgets (char *s, size_t size, const char *buf)
{
    enum {
	GOT_LF = 1,
	GOT_CR,
	GOT_CRLF
    };
    int i, status = 0;
    const char *p;

    if (s == NULL && size == 1) {
	/* signal for end-of-read */
	rbuf_finalize(buf);
	return NULL;
    }

    if (s == NULL || size == 0) {
	/* signal for initialization */
	rbuf_push(buf);
	return NULL;
    }

    p = rbuf_get_point(buf);
    if (p == NULL) {
	return NULL;
    }

    if (*p == '\0') {
	/* reached the end of the buffer */
	return NULL;
    }

    *s = 0;
    /* advance to line-end, end of buffer, or maximum size,
       whichever comes first */
    for (i=0; ; i++) {
	s[i] = p[i];
	if (p[i] == '\0') {
	    break;
	}
	if (p[i] == '\r') {
	    s[i] = '\0';
	    if (p[i+1] == '\n') {
		status = GOT_CRLF;
	    } else {
		status = GOT_CR;
	    }
	    break;
	}
	if (p[i] == '\n') {
	    s[i] = '\0';
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
    if (status == GOT_CR || status == GOT_LF) {
	p++;
    } else if (status == GOT_CRLF) {
	p += 2;
    }

    /* replace newline */
    if (status && i < size - 1) {
	strcat(s, "\n");
    }

    rbuf_set_point(buf, p);

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

/**
 * bufgets_finalize:
 * @buf: source buffer.
 *
 * Signals that we are done reading from @buf.
 */

void bufgets_finalize (const char *buf)
{
    bufgets(NULL, 1, buf);
}
