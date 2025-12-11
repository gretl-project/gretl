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
#include "uservar.h"
#include "gretl_string_table.h"
#include "gretl_midas.h"
#include "matrix_extra.h"

#include <time.h>

#define PDEBUG 0

static int gretl_digits = 6;

#define have_range_option(o) (o & (OPT_R | OPT_H | OPT_T))

static int get_print_range (gretlopt opt, int len,
                            int *start, int *stop);

void bufspace (int n, PRN *prn)
{
    while (n-- > 0) {
	pputc(prn, ' ');
    }
}

void dashline (int n, PRN *prn)
{
    while (n-- > 0) {
	pputc(prn, '-');
    }
    pputc(prn, '\n');
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

    sprintf(str, "%#*.*g", d, gretl_digits, xx);
}

static void covhdr (PRN *prn)
{
    pprintf(prn, "%s:\n\n",
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
    char timestr[48];
    PRN *myprn = NULL;

    if (prn == NULL) {
	myprn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
	prn = myprn;
    }

    print_time(timestr);
    pprintf(prn, "%s: %s\n", _("Current session"), timestr);

    if (myprn != NULL) {
	gretl_print_destroy(myprn);
    }
}

/**
 * logo:
 * @quiet: if non-zero, just print version ID, else print
 * copyright info also.
 *
 * Print to stdout gretl version information.
 */

void logo (int quiet)
{
    printf(_("gretl version %s\n"), GRETL_VERSION);

    if (!quiet) {
	puts(_("Copyright Ramu Ramanathan, Allin Cottrell and Riccardo \"Jack\" Lucchetti"));
	puts(_("This is free software with ABSOLUTELY NO WARRANTY"));
    }
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
    pputc(prn, '\n');
    pputs(prn, _("This is free software with ABSOLUTELY NO WARRANTY"));
    pputc(prn, '\n');

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
    printf("\nLibgretl-1.0, version %d\n", LIBGRETL_CURRENT);
}

/**
 * gui_script_logo:
 * @prn: gretl printing struct.
 *
 * Print to @prn a header for script output in gui program.
 */

void gui_script_logo (PRN *prn)
{
    char timestr[48];

    pprintf(prn, _("gretl version %s\n"), GRETL_VERSION);
    print_time(timestr);
    pprintf(prn, "%s: %s\n", _("Current session"), timestr);
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
 * max_obs_marker_length:
 * @dset: dataset information.
 *
 * Returns: the length of the longest observation marker
 * within the current sample range.
 */

int max_obs_marker_length (const DATASET *dset)
{
    char s[OBSLEN];
    int t, n, nmax = 0;

    if (dset->S != NULL) {
	/* we have specific observation strings */
	for (t=dset->t1; t<=dset->t2; t++) {
	    get_obs_string(s, t, dset);
	    n = g_utf8_strlen(s, -1);
	    if (n > nmax) {
		nmax = n;
	    }
	    if (nmax == OBSLEN - 1) {
		break;
	    }
	}
    } else if (dated_daily_data(dset)) {
	get_obs_string(s, dset->t2, dset);
	nmax = strlen(s);
    } else if (dataset_is_time_series(dset)) {
	switch (dset->pd) {
	case 1:   /* annual: YYYY */
	case 10:  /* decennial: YYYY */
	    nmax = 4;
	    break;
	case 4:   /* quarterly: YYYY:Q */
	    nmax = 6;
	    break;
	case 12:  /* monthly: YYYY:MM */
	    nmax = 7;
	    break;
	default:
	    break;
	}
	if (nmax == 0) {
	    get_obs_string(s, dset->t2, dset);
	    nmax = strlen(s);
	}
    } else {
	int T = dset->t2 - dset->t1 + 1;
	int incr = (T < 120)? 1 : (T / 100.0);

	for (t=dset->t1; t<=dset->t2; t+=incr) {
	    get_obs_string(s, t, dset);
	    n = strlen(s);
	    if (n > nmax) {
		nmax = n;
	    }
	}
    }

    return nmax;
}

void print_freq_test (const FreqDist *freq, PRN *prn)
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

    pputc(prn, '\n');

    if (!na(pval)) {
	record_test_result(freq->test, pval);
    }
}

/**
 * print_freq:
 * @freq: gretl frequency distribution struct.
 * @varno: ID number of the series in question.
 * @dset: pointer to dataset.
 * @prn: gretl printing struct.
 *
 * Print frequency distribution to @prn.
 */

void print_freq (const FreqDist *freq, int varno, const DATASET *dset,
		 PRN *prn)
{
    int i, k, nlw, K;
    int total, valid;
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

    if (freq->strvals) {
	int len, maxlen = 0;

	for (i=0; i<freq->numbins; i++) {
	    /* FIXME UTF-8? */
	    len = g_utf8_strlen(freq->S[i], -1);
	    if (len > maxlen) {
		maxlen = len;
	    }
	}
	len = maxlen > 31 ? 31 : maxlen < 10 ? 10 : maxlen;

	pputc(prn, '\n');
	bufspace(len, prn);
	pputs(prn, _("frequency   percent\n\n"));
	for (k=0; k<=K; k++) {
	    *word = '\0';
	    gretl_utf8_strncat(word, freq->S[k], len-2);
	    pputs(prn, word);
	    nlw = len - g_utf8_strlen(word, -1);
	    bufspace(nlw, prn);
	    pprintf(prn, "%6d   ", freq->f[k]);
	    f = 100.0 * freq->f[k] / valid;
	    pprintf(prn, "  %6.2f%% ", f);
	    if (f < 100) {
		i = 0.36 * f;
		while (i--) {
		    pputc(prn, '*');
		}
	    }
	    pputc(prn, '\n');
	}
    } else if (freq->discrete) {
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
	    if (f < 100) {
		i = 0.36 * f;
		while (i--) {
		    pputc(prn, '*');
		}
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

	    x = (k == K && K > 0)? freq->endpt[k] : freq->endpt[k+1];
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
	    if (K > 1) {
		while (i--) {
		    pputc(prn, '*');
		}
	    }
	    pputc(prn, '\n');
	}
    }

    if (valid < total) {
	int missing = total - valid;

	pprintf(prn, "\n%s = %d (%5.2f%%)\n", _("Missing observations"),
		missing, 100 * (double) missing / total);
    }

    if (!na(freq->test)) {
	print_freq_test(freq, prn);
    } else {
	pputc(prn, '\n');
    }
}

static void tex_xtab_heading (const Xtab *tab, PRN *prn)
{
    char s1[9 + VNAMELEN * 2] = {0};
    char s2[9 + VNAMELEN * 2] = {0};
    const char *src;
    char *targ;
    int i, j, k;

    strcpy(s1, "\\texttt{");
    strcpy(s2, "\\texttt{");

    /* fix any underscores in series names */
    for (k=0; k<2; k++) {
	targ = k == 0 ? s1 : s2;
	src = k == 0 ? tab->rvarname : tab->cvarname;
	for (i=0, j=8; src[i]; i++, j++) {
	    if (src[i] == '_') {
		strcat(targ, "\\_");
		j++;
	    } else {
		targ[j] = src[i];
	    }
	}
	strcat(targ, "}");
    }

    pprintf(prn, _("Cross-tabulation of %s (rows) against %s (columns)"),
	    s1, s2);
    pputs(prn, "\n\n\\vspace{1em}\n\n");
}

static int row_strlen (const Xtab *tab)
{
    int i, n, nmax = 0;

    for (i=0; i<tab->rows; i++) {
	n = g_utf8_strlen(tab->Sr[i], -1);
	if (n > nmax) {
	    nmax = n;
	}
    }

    return nmax;
}

static int col_strlen (const Xtab *tab)
{
    int j, n, nmax = 0;

    for (j=0; j<tab->cols; j++) {
	n = g_utf8_strlen(tab->Sc[j], -1);
	if (n > nmax) {
	    nmax = n;
	}
    }

    return nmax;
}

/* allow for wider columns when cell counts are large numbers */

static void get_xtab_col_widths (const Xtab *tab, int stdwidth,
				 int *cwidth, int *twidth)
{
    int i, di, dmax = 0;

    /* column width */
    for (i=0; i<tab->cols; i++) {
	di = (int) (floor(log10(tab->ctotal[i])));
	dmax = di > dmax ? di : dmax;
    }
    *cwidth = (dmax >= stdwidth - 2)? dmax + 3 : stdwidth;

    /* totals width */
    dmax = 0;
    for (i=0; i<tab->rows; i++) {
	di = (int) (floor(log10(tab->rtotal[i])));
	dmax = di > dmax ? di : dmax;
    }
    *twidth = (dmax >= stdwidth - 2)? dmax + 3 : stdwidth;
}

static void real_print_xtab (const Xtab *tab, const DATASET *dset,
			     gretlopt opt, PRN *prn)
{
    double x, y, cj, ri;
    int n5 = 0;
    double ymin = 1.0e-7;
    double pearson = 0.0;
    double pval = NADBL;
    char lbl[64];
    int dashlen = 0;
    int rlen = 0;
    int clen = 0;
    int cw, tw;
    int stdw = 6;
    int totals = 1;
    int tex = 0;
    int bold = 0;
    int i, j;

    if (opt & OPT_T) {
	/* LaTeX output */
	tex = 1;
	if ((opt & OPT_E) && !(tab->rstrs || tab->cstrs)) {
	    /* bold-face equal values */
	    bold = 1;
	}
	/* don't bother with chi-square */
	pearson = NADBL;
    }

    if (opt & OPT_N) {
	totals = 0;
    }

    if (*tab->rvarname != '\0' && *tab->cvarname != '\0') {
	pputc(prn, '\n');
	if (tex) {
	    tex_xtab_heading(tab, prn);
	} else {
	    pprintf(prn, _("Cross-tabulation of %s (rows) against %s (columns)"),
		    tab->rvarname, tab->cvarname);
	    pputs(prn, "\n\n");
	}
    } else {
	pputc(prn, '\n');
    }

    if (tex) {
	int nc = tab->cols + 1 + totals;

	pputs(prn, "\\begin{tabular}{");
	pputs(prn, "r|");
	for (j=1; j<nc; j++) {
	    pputc(prn, 'r');
	}
	pputs(prn, "}\n");
    }

    if (tex) {
	pputs(prn, "     & ");
    }

    if (tab->Sr != NULL) {
	rlen = 2 + row_strlen(tab);
	if (rlen > 16) {
	    rlen = 16;
	} else if (rlen < 7) {
	    rlen = 7;
	}
    } else {
	rlen = 7;
    }

    /* allow for BIG integers */
    get_xtab_col_widths(tab, stdw, &cw, &tw);

    if (tab->Sc != NULL) {
	clen = 2 + col_strlen(tab);
	if (clen > 12) {
	    clen = 12;
	} else if (clen < cw) {
	    clen = 6;
	}
    } else {
	clen = cw;
    }

    if (totals) {
        dashlen = clen + 1 + tab->cols*cw + tw;
    }

    bufspace(rlen, prn);

    /* header row: column labels */

    for (j=0; j<tab->cols; j++) {
	if (tab->Sc != NULL) {
	    *lbl = '\0';
	    strncat(lbl, tab->Sc[j], 63);
	    gretl_utf8_truncate(lbl, clen-2);
	    if (tex) {
		pputs(prn, lbl);
	    } else {
		bufspace(clen - g_utf8_strlen(lbl, -1), prn);
		pputs(prn, lbl);
	    }
	} else {
	    cj = tab->cval[j];
	    if (tex) {
		pprintf(prn, "%4g", cj);
	    } else {
		pprintf(prn, "[%*g]", cw-2, cj);
	    }
	}
	if (tex) {
	    if (j < tab->cols - 1 || totals) {
		pputs(prn, " & ");
	    } else {
		pputs(prn, "\\\\ \\hline\n");
	    }
	}
    }

    if (totals) {
	if (tex) {
	    pprintf(prn,"$\\Sigma$\\\\ \\hline\n");
	} else {
	    bufspace(2 + (tw - stdw), prn);
	    pprintf(prn,"%s\n", _("TOT."));
            dashline(dashlen, prn);
	}
    } else if (!tex) {
	pputc(prn, '\n');
    }

    /* body of table */

    for (i=0; i<tab->rows; i++) {
	if (tab->rtotal[i] == 0) {
	    continue;
	}
	if (tab->Sr != NULL) {
	    *lbl = '\0';
	    strncat(lbl, tab->Sr[i], 63);
	    gretl_utf8_truncate(lbl, rlen-2);
	    pputs(prn, lbl);
	    if (!tex) {
		bufspace(rlen - g_utf8_strlen(lbl, -1), prn);
	    }
	} else {
	    ri = tab->rval[i];
	    if (tex) {
		pprintf(prn, "%4g", ri);
	    } else {
		pprintf(prn, "[%4g] ", ri);
	    }
	}
	if (tex) {
	    pputs(prn, " & ");
	}
	/* row counts */
	for (j=0; j<tab->cols; j++) {
	    if (tab->ctotal[j] > 0) {
		if (tab->f[i][j] || (opt & OPT_Z)) {
		    if (opt & (OPT_C | OPT_R)) {
			if (opt & OPT_C) {
			    x = 100.0 * tab->f[i][j] / tab->ctotal[j];
			} else {
			    x = 100.0 * tab->f[i][j] / tab->rtotal[i];
			}
			if (tex) {
			    /* FIXME strvals! */
			    if (bold && tab->cval[j] == tab->rval[i]) {
				pprintf(prn, "\\textbf{%.1f%%%%}", x);
			    } else {
				pprintf(prn, "%*.1f%%%%", clen, x);
			    }
			} else {
			    pprintf(prn, "%*.1f%%", clen-1, x);
			}
		    } else {
			if (bold && tab->cval[j] == tab->rval[i]) {
			    pprintf(prn, "\\textbf{%d} ", tab->f[i][j]);
			} else {
			    pprintf(prn, "%*d", clen, tab->f[i][j]);
			}
		    }
		} else if (!tex) {
		    bufspace(clen, prn);
		}
		if (tex && (totals || j < tab->cols-1)) {
		    pputs(prn, "& ");
		}
	    }
	    if (!na(pearson)) {
		/* cumulate chi-square */
		y = ((double) tab->rtotal[i] * tab->ctotal[j]) / tab->n;
		if (y < ymin) {
		    pearson = NADBL;
		} else {
		    x = (double) tab->f[i][j] - y;
		    pearson += x * x / y;
		    if (y >= 5) {
			n5++;
		    }
		}
	    }
	}
	if (totals) {
	    /* row totals */
	    if (opt & OPT_C) {
		x = 100.0 * tab->rtotal[i] / tab->n;
		if (tex) {
		    pprintf(prn, "%5.1f%%%%", x);
		} else {
		    pprintf(prn, "%5.1f%%", x);
		}
	    } else {
		pprintf(prn, "%*d", tw, tab->rtotal[i]);
	    }
	}
	/* terminate row */
	if (tex) {
	    if (totals && i == tab->rows-1) {
		pputs(prn, "\\\\ [2pt]\n");
	    } else {
		pputs(prn, "\\\\\n");
	    }
	} else {
	    pputc(prn, '\n');
	}
    }

    /* footer row */

    if (totals) {
	/* column totals */
	if (tex) {
	    pputs(prn, "$\\Sigma$ & ");
	} else {
            if (dashlen > 0) {
                dashline(dashlen, prn);
            } else {
                pputc(prn, '\n');
            }
	    pprintf(prn, "%-*s", rlen, _("TOTAL"));
	}
	for (j=0; j<tab->cols; j++) {
	    if (opt & OPT_R) {
		x = 100.0 * tab->ctotal[j] / tab->n;
		if (tex) {
		    pprintf(prn, "%*.1f%%%%", clen, x);
		} else {
		    pprintf(prn, "%*.1f%%", clen-1, x);
		}
	    } else {
		pprintf(prn, "%*d", clen, tab->ctotal[j]);
	    }
	    if (tex) {
		pputs(prn, "& ");
	    }
	}
	pprintf(prn, "%*d\n", tw, tab->n);
    }

    if (tex) {
	pputs(prn, "\\end{tabular}\n");
	return;
    }

    /* additional information, if applicable */

    if (tab->missing) {
	pputc(prn, '\n');
	pprintf(prn, _("%d missing values"), tab->missing);
	pputc(prn, '\n');
    }

    if (na(pearson)) {
	pputc(prn, '\n');
	pprintf(prn, _("Pearson chi-square test not computed: some "
		       "expected frequencies were less\n"
		       "than %g\n"), ymin);
    } else {
	double n5p = (double) n5 / (tab->rows * tab->cols);
	int df = (tab->rows - 1) * (tab->cols - 1);

	pval = chisq_cdf_comp(df, pearson);
	if (na(pval)) {
	    pearson = NADBL;
	} else {
	    pputc(prn, '\n');
	    pprintf(prn, _("Pearson chi-square test = %g (%d df, p-value = %g)"),
		    pearson, df, pval);
	    pputc(prn, '\n');
	    if (!tex && n5p < 0.80) {
		/* xgettext:no-c-format */
		pputs(prn, _("Warning: Less than of 80% of cells had expected "
			     "values of 5 or greater.\n"));
	    }
	}
    }

    if (opt & OPT_S) {
	/* saving Pearson test result */
	record_test_result(pearson, pval);
    }

    if (!(opt & OPT_F) && !tex && tab->rows == 2 && tab->cols == 2) {
	fishers_exact_test(tab, prn);
    }
}

/**
 * print_xtab:
 * @tab: gretl cross-tabulation struct.
 * @dset: pointer to dataset, or NULL.
 * @opt: may contain %OPT_R to print row percentages, %OPT_C
 * to print column percentages, %OPT_Z to display zero entries,
 * %OPT_T to print as TeX (LaTeX), %OPT_N to omit marginal
 * totals, %OPT_B to bold-face counts where the row and
 * column values are equal (TeX only), %OPT_S to record
 * Pearson test result.
 * @prn: gretl printing struct.
 *
 * Print crosstab to @prn.
 */

void print_xtab (const Xtab *tab, const DATASET *dset,
		 gretlopt opt, PRN *prn)
{
    int done = 0;

    if (opt & OPT_T) {
	/* --tex : are we printing to file? */
	const char *fname = get_optval_string(XTAB, OPT_T);

	if (fname != NULL && *fname != '\0') {
	    int err = 0;
	    PRN *xprn;

	    gretl_maybe_switch_dir(fname);
	    xprn = gretl_print_new_with_filename(fname, &err);
	    if (err) {
		pprintf(prn, _("Couldn't write to %s"), fname);
		pputc(prn, '\n');
	    } else {
		real_print_xtab(tab, dset, opt, xprn);
		gretl_print_destroy(xprn);
	    }
	    done = 1;
	}
    }

    if (!done) {
	real_print_xtab(tab, dset, opt, prn);
    }
}

/**
 * print_smpl:
 * @dset: data information struct
 * @fulln: full length of data series, if dataset is
 * subsampled, or 0 if not applicable/known.
 * @opt: may include OPT_F to force showing the sample
 * range even when "messages" are turned off.
 * @prn: gretl printing struct.
 *
 * Prints the current sample information to @prn.
 */

void print_smpl (const DATASET *dset, int fulln,
		 gretlopt opt, PRN *prn)
{
    if (dset == NULL || dset->v == 0 || prn == NULL) {
	return;
    }

    if (!(opt & OPT_F) && !gretl_messages_on()) {
	/* hush */
	return;
    }

    if (fulln && !dataset_is_panel(dset)) {
	pprintf(prn, _("Full data set: %d observations\n"), fulln);
	if (sample_size(dset) < dset->n ||
	    (dataset_is_time_series(dset) && dset->sd0 > 1)) {
	    print_sample_obs(dset, prn);
	} else {
	    pprintf(prn, _("Current sample: %d observations\n"),
		    sample_size(dset));
	}
	return;
    }

    if (fulln) {
	pprintf(prn, _("Full data set: %d observations\n"), fulln);
    } else {
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Full data range"),
		dset->stobs, dset->endobs, dset->n);
    }

    if (dset->t1 > 0 || dset->t2 < dset->n - 1 ||
	(fulln && dataset_is_panel(dset))) {
	print_sample_obs(dset, prn);
    }

    pputc(prn, '\n');
}

static void print_var_smpl (int v, const DATASET *dset, PRN *prn)
{
    int t, n = 0;

    if (dset->t1 > 0 || dset->t2 < dset->n - 1) {
	char d1[OBSLEN], d2[OBSLEN];
	ntolabel(d1, dset->t1, dset);
	ntolabel(d2, dset->t2, dset);

	pprintf(prn, "%s:  %s - %s", _("Current sample"), d1, d2);
    } else {
	pprintf(prn, "%s: %s - %s", _("Full data range"),
		dset->stobs, dset->endobs);
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	if (!na(dset->Z[v][t])) {
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
 * Returns: the corrected numeric string.
 */

char *gretl_fix_exponent (char *s)
{
    char *p;
    int n;

    if ((p = strstr(s, "+00")) || (p = strstr(s, "-00"))) {
	if (*(p+3)) {
	    memmove(p+1, p+2, strlen(p+1));
	}
    }

    n = strlen(s);
    if (s[n-1] == '.' || s[n-1] == ',') {
	/* delete trailing junk */
	s[n-1] = '\0';
    }

    return s;
}

/* determine the max number of characters that will be output when
   printing Z[v] over the current sample range using format %.*f
   or %#.*g, with precision 'digits'
*/

static int max_number_length (int v, const DATASET *dset,
			      char fmt, int digits)
{
    double a, x, amax = 0.0, amin = 1.0e300;
    int t, n, maxsgn = 0, minsgn = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
	x = dset->Z[v][t];
	if (!na(x)) {
	    a = fabs(x);
	    if (a > amax) {
		amax = a;
		maxsgn = (x < 0);
	    }
	    if (fmt == 'g' && a < amin) {
		amin = a;
		minsgn = (x < 0);
	    }
	}
    }

    if (fmt == 'f') {
	if (amax <= 1.0) {
	    n = 1;
	} else {
	    n = ceil(log10(amax)) + (fmod(amax, 10) == 0);
	}
	n += digits + 1 + maxsgn;
    } else {
	double l10 = log10(amax);
	int amaxn = digits + 1, aminn = digits + 1;

	if (l10 >= digits) {
	    amaxn += 5 + maxsgn;
	}
	l10 = log10(amin);
	if (l10 < -4) {
	    aminn += 5 + minsgn;
	} else if (l10 < 0) {
	    aminn += (int) ceil(-l10) + minsgn;
	}
	n = (amaxn > aminn)? amaxn : aminn;
#if 0
	fprintf(stderr, "var %d, amax=%g, amin=%g, n=%d\n",
		v, amax, amin, n);
#endif
    }

    return n;
}

static int series_column_width (int v, const DATASET *dset,
				char fmt, int digits)
{
    int namelen = strlen(dset->varname[v]);
    int numlen = max_number_length(v, dset, fmt, digits);

    return (namelen > numlen)? namelen : numlen;
}


/* For some reason sprintf using "%#G" seems to stick an extra
   zero on the end of some numbers -- i.e. when using a precision
   of 6 you can get a result of "1.000000", with 6 trailing
   zeros.  The following function checks for this and lops it
   off if need be.
*/

static char *cut_extra_zero (char *s, int digits)
{
    if (strchr(s, 'E') == NULL && strchr(s, 'e') == NULL) {
	int n = strspn(s, "-.,0");
	int m = (strchr(s + n, '.') || strchr(s + n, ','));

	s[n + m + digits] = '\0';
    }

    return s;
}

static char *cut_trailing_point (char *s)
{
    int n = strlen(s);

    if (s[n-1] == '.' || s[n-1] == ',') {
	s[n-1] = '\0';
    }

    return s;
}

/* below: targ should be 36 bytes long */

void gretl_sprint_fullwidth_double (double x, int digits, char *targ,
				    PRN *prn)
{
    char decpoint;
    int n;

    *targ = '\0';

    if (na(x)) {
	strcpy(targ, "NA");
	return;
    }

    decpoint = get_local_decpoint();

    if (digits == -4) {
	if (x < .0001 && x > 0.0) {
	    sprintf(targ, "%#.3g", x);
	    digits = 3;
	} else {
	    sprintf(targ, "%.4f", x);
	    return;
	}
    } else {
	/* let's not print non-zero values for numbers smaller than
	   machine zero */
	x = screen_zero(x);
	sprintf(targ, "%#.*g", digits, x);
    }

    gretl_fix_exponent(targ);

    n = strlen(targ) - 1;
    if (targ[n] == decpoint) {
	targ[n] = '\0';
    }

    cut_extra_zero(targ, digits);

    if (*targ == '-' && gretl_print_has_minus(prn)) {
	char tmp[36];

	strcpy(tmp, targ + 1);
	*targ = '\0';
	strcat(targ, "−"); /* U+2212: minus */
	strcat(targ, tmp);
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
    char decpoint;
    char *p;

    decpoint = get_local_decpoint();

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

int get_gretl_digits (void)
{
    return gretl_digits;
}

int set_gretl_digits (int d)
{
    if (d >= 3 && d <= 6) {
	gretl_digits = d;
	return 0;
    } else {
	/* In the contexts where we're using @gretl_digits
	   we don't want it to be less than 3 or greater
	   than 6.
	*/
	return E_INVARG;
    }
}

void gretl_print_value (double x, PRN *prn)
{
    gretl_print_fullwidth_double(x, gretl_digits, prn);
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
	pprintf(prn, "%s \\\\\n", _("Cross-equation VCV for residuals"));
	pprintf(prn, "(%s)\n\n", _("correlations above the diagonal"));
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
		pprintf(prn, "\n%s = ", _("log determinant"));
		pprintf(prn, "$-$%g\n", -ldet);
	    } else {
		pprintf(prn, "\n%s = %g\n", _("log determinant"), ldet);
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
 * @dset: data information struct.
 * @prn: gretl printing struct.
 *
 * Print to @prn the variance-covariance matrix for the parameter
 * estimates in @pmod.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int outcovmx (MODEL *pmod, const DATASET *dset, PRN *prn)
{
    VMatrix *vmat;
    int err = 0;

    vmat = gretl_model_get_vcv(pmod, dset);

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

static int vmat_maxlen (VMatrix *vmat)
{
    int i, len, maxlen = 0;

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
    int n, nf, li2, p, k, idx, ij2;
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

    for (i=0; i<=vmat->dim/fields; i++) {
	nf = i * fields;
	li2 = vmat->dim - nf;
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    s = vmat->names[j + nf - 1];
	    n = strlen(s);
	    if (n > fwidth - 1) {
		pprintf(prn, " %.*s~", fwidth - 2, s);
	    } else {
		bufspace(fwidth - n, prn);
		pputs(prn, s);
	    }
	}
	pputc(prn, '\n');

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    for (k=0; k<p; k++) {
		idx = ijton(j, nf+k, vmat->dim);
		outxx(vmat->vec[idx], vmat->ci, fwidth, prn);
	    }
	    if (fwidth < 15) pputc(prn, ' ');
	    n = strlen(vmat->names[j]);
	    if (n > fwidth - 1) {
		pprintf(prn, " %.*s~\n", fwidth - 2, vmat->names[j]);
	    } else {
		pprintf(prn, " %s\n", vmat->names[j]);
	    }
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    ij2 = nf + j;
	    bufspace(fwidth * j, prn);
	    for (k=j; k<p; k++) {
		idx = ijton(ij2, nf+k, vmat->dim);
		outxx(vmat->vec[idx], vmat->ci, fwidth, prn);
	    }
	    if (fwidth < 15) pputc(prn, ' ');
	    n = strlen(vmat->names[ij2]);
	    if (n > fwidth - 1) {
		pprintf(prn, " %.*s~\n", fwidth - 2, vmat->names[ij2]);
	    } else {
		pprintf(prn, " %s\n", vmat->names[ij2]);
	    }
	}
	pputc(prn, '\n');
    }
}

static int fit_resid_head (const FITRESID *fr,
			   const DATASET *dset,
			   int obslen,
			   PRN *prn)
{
    char label[32];
    char obs1[OBSLEN], obs2[OBSLEN];
    int kstep = fr->method == FC_KSTEP;
    int ywidth;

    if (kstep) {
	ntolabel(obs1, fr->model_t1, dset);
	pprintf(prn, _("Recursive %d-step ahead forecasts"), fr->k);
	pputs(prn, "\n\n");
	pprintf(prn, _("The forecast for time t is based on (a) coefficients obtained by\n"
		       "estimating the model over the sample %s to t-%d, and (b) the\n"
		       "regressors evaluated at time t."), obs1, fr->k);
	pputs(prn, "\n\n");
	pputs(prn, _("This is truly a forecast only if all the stochastic regressors\n"
		     "are in fact lagged values."));
	pputs(prn, "\n\n");
    } else {
	ntolabel(obs1, fr->t1, dset);
	ntolabel(obs2, fr->t2, dset);
	pprintf(prn, _("Model estimation range: %s - %s"), obs1, obs2);
	pputc(prn, '\n');

	if (fr->std) {
	    pprintf(prn, "%s\n", _("The residuals are standardized"));
	} else if (!na(fr->sigma)) {
	    pprintf(prn, "%s = %.*g\n", _("Standard error of the regression"),
		    gretl_digits, fr->sigma);
	}
    }

    pputc(prn, '\n');
    bufspace(obslen, prn);

    /* column 1 */
    maybe_trim_varname(label, fr->depvar);
    ywidth = strlen(label) + 1;
    if (ywidth < 13) {
	ywidth = 13;
    }
    pprintf(prn, "%*s", ywidth, label);

    /* column 2 */
    strcpy(label, (kstep)? _("forecast") : _("fitted"));
    pprintf(prn, "%*s", UTF_WIDTH(label, 13), label);

    /* column 3 */
    strcpy(label, (kstep)? _("error") : _("residual"));
    pprintf(prn, "%*s", UTF_WIDTH(label, 13), label);

    pputs(prn, "\n\n");

    return ywidth;
}

static const char *get_series_name (const DATASET *dset, int i,
				    gchar **altname, int debug)
{
    gchar *showname = NULL;
    const char *sname = dset->varname[i];
    const char *lname;
    const char *ret;

    if (debug) {
	ret = sname;
    } else if (series_is_listarg(dset, i, &lname)) {
	if (lname != NULL) {
	    showname = g_strdup_printf("%s.%s", lname, sname);
	    *altname = showname;
	    ret = showname;
	} else {
	    ret = "[masked]";
	}
    } else {
	ret = sname;
    }

    return ret;
}

/* prints a heading with the names of the variables in @list */

static void varheading (const int *list, int leader,
			const int *wid, const DATASET *dset,
			char delim, PRN *prn)
{
    const char *name;
    gchar *tmp;
    int i, vi;

    if (delim) {
	if (leader >= 0) {
	    pprintf(prn, "obs%c", delim);
	}
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    pputs(prn, dset->varname[vi]);
	    if (i < list[0]) {
		pputc(prn, delim);
	    }
	}
	pputc(prn, '\n');
    } else if (rtf_format(prn)) {
	pputs(prn, "{obs\\cell ");
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    pprintf(prn, "%s\\cell ", dset->varname[vi]);
	}
	pputs(prn, "}\n");
    } else {
	pputc(prn, '\n');
	bufspace(leader, prn);
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    tmp = NULL;
	    name = get_series_name(dset, vi, &tmp, 0);
	    pprintf(prn, "%*s", wid[i], name);
	    g_free(tmp);
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

static void print_stringvals_for_var (const DATASET *dset, int v, PRN *prn)
{
    const char *s;
    int ls = 0;
    int t, n;

    for (t=dset->t1; t<=dset->t2; t++) {
	s = series_get_string_for_obs(dset, v, t);
	if (s == NULL) {
	    s = "(null)";
	}
	n = strlen(s);
	if (ls + n > 78) {
	    pputc(prn, '\n');
	    ls = 0;
	}
	pprintf(prn, "%s ", s);
	ls += n;
    }

    pputc(prn, '\n');
}

/* prints series @v from current sample t1 to t2 */

static void print_series_by_var (const DATASET *dset, int v, PRN *prn)
{
    const double *z = dset->Z[v];
    char format[12];
    int t, ls = 0;
    int anyneg = 0;
    double x;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (z[t] < 0) {
	    anyneg = 1;
	    break;
	}
    }

    if (anyneg) {
	sprintf(format, "%% #.%dg  ", gretl_digits);
    } else {
	sprintf(format, "%%#.%dg  ", gretl_digits);
    }

    for (t=dset->t1; t<=dset->t2; t++) {
	char str[32];
	int n;

	x = z[t];

	if (na(x)) {
	    sprintf(str, "%*s  ", gretl_digits + 1 + anyneg, "NA");
	} else if (isnan(x)) {
	    sprintf(str, "%*s  ", gretl_digits + 1 + anyneg, "NaN");
	} else if (isinf(x)) {
	    sprintf(str, "%*s  ", gretl_digits + 1 + anyneg,
		    (x < 0)? "-inf" : "inf");
	} else {
	    sprintf(str, format, x);
	}

	n = strlen(str);
	if (ls + n > 78) {
	    pputc(prn, '\n');
	    ls = 0;
	}

	pputs(prn, str);
	ls += n;
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
    char decpoint;

    decpoint = get_local_decpoint();

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
#if PDEBUG
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
#if PDEBUG
	fprintf(stderr, "get_signif: setting smax = -(%d - %d)\n",
		smax, leadmax);
#endif
	smax = -1 * (smax - leadmax); /* # of decimal places */
    } else if (leadmax == smax) {
	smax = 0;
    } else if (leadmax == 0 && !allfrac) {
#if PDEBUG
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

static char *bufprintnum (char *buf, double x, int signif,
			  int gprec, int width)
{
    static char numstr[32];
    int i, l;

    *buf = '\0';

    if (isnan(x)) {
	strcpy(numstr, "NaN");
	goto finish;
    } else if (isinf(x)) {
	strcpy(numstr, (x < 0)? "-inf" : "inf");
	goto finish;
    }

    /* guard against monster numbers that will smash the stack */
    if (fabs(x) > 1.0e20 || signif == PMAX_NOT_AVAILABLE) {
	sprintf(numstr, "%.*g", gprec, x);
	goto finish;
    }

    if (signif < 0) {
	sprintf(numstr, "%.*f", -signif, x);
    } else if (signif == 0) {
	sprintf(numstr, "%.0f", x);
    } else {
	double a = fabs(x);

	if (a < 1) l = 0;
	else if (a < 10) l = 1;
	else if (a < 100) l = 2;
	else if (a < 1000) l = 3;
	else if (a < 10000) l = 4;
	else if (a < 100000) l = 5;
	else if (a < 1000000) l = 6;
	else l = 7;

#if PDEBUG
	fprintf(stderr, "%g: got %d for leftvals, %d for signif\n",
		x, l, signif);
#endif

	if (l == 6 && signif < 6) {
	   sprintf(numstr, "%.0f", x);
	} else if (l >= signif) {
#if PDEBUG
	    fprintf(stderr, " printing with '%%.%dG'\n", signif);
#endif
	    if (g_too_long(x, signif)) {
		sprintf(numstr, "%.0f", x);
	    } else {
		sprintf(numstr, "%.*G", signif, x);
	    }
	} else if (a >= .10) {
#if PDEBUG
	    fprintf(stderr, " printing with '%%.%df'\n", signif-l);
#endif
	    sprintf(numstr, "%.*f", signif - l, x);
	} else {
	    if (signif > 4) {
		signif = 4;
	    }
#if PDEBUG
	    fprintf(stderr, " printing with '%%#.%dG'\n", signif);
#endif
	    sprintf(numstr, "%#.*G", signif, x); /* # wanted? */
	}
    }

 finish:

    if (width > 0) {
	/* pad on left as needed */
	l = width - strlen(numstr);
	for (i=0; i<l; i++) {
	    strcat(buf, " ");
	}
    }
    strcat(buf, numstr);

    return buf;
}

static void real_print_obs_marker (int t, const DATASET *dset,
				   int len, int pad, PRN *prn)
{
    char tmp[OBSLEN] = {0};
    int thislen = len;

    if (dataset_has_markers(dset)) {
	strcpy(tmp, dset->S[t]);
	thislen = get_utf_width(tmp, len);
    } else {
	ntolabel(tmp, t, dset);
    }

    if (pad) {
	pprintf(prn, "%*s ", thislen, tmp);
    } else {
	pprintf(prn, "%*s", thislen, tmp);
    }
}

/**
 * print_obs_marker:
 * @t: observation number.
 * @dset: data information struct.
 * @len: width to which to print.
 * @prn: gretl printing struct.
 *
 * Print a string (label, date or obs number) representing the given @t.
 */

void print_obs_marker (int t, const DATASET *dset, int len, PRN *prn)
{
    real_print_obs_marker(t, dset, len, 1, prn);
}

char *maybe_trim_varname (char *targ, const char *src)
{
    int srclen = strlen(src);

    if (srclen < NAMETRUNC) {
	strcpy(targ, src);
    } else {
	const char *p = strrchr(src, '_');

	*targ = '\0';

	if (p != NULL && isdigit(*(p+1)) && strlen(p) < 4) {
	    /* preserve lag identifier? */
	    int snip = srclen - NAMETRUNC + 2;
	    int fore = p - src;

	    strncat(targ, src, fore - snip);
	    strcat(targ, "~");
	    strcat(targ, p);
	} else {
	    strncat(targ, src, NAMETRUNC - 2);
	    strcat(targ, "~");
	}
    }

    return targ;
}

int max_namelen_in_list (const int *list, const DATASET *dset)
{
    int i, vi, ni, n = 0;

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (vi >= 0 && vi < dset->v) {
	    ni = strlen(dset->varname[list[i]]);
	    if (ni > n) {
		n = ni;
	    }
	}
    }

    if (n >= NAMETRUNC) {
	n = NAMETRUNC - 1;
    }

    return n;
}

static int show_series (int i, int fd, const DATASET *dset,
			int debug)
{
    if (debug) {
	/* show all */
	return 1;
    } else {
	return (i == 0 || fd == 0 || series_get_stack_level(dset, i) == fd);
    }
}

/**
 * list_series:
 * @dset: data information struct.
 * @opt: may contain OPT_D for debugging.
 * @prn: gretl printing struct
 *
 * Prints a list of the names of the series currently defined.
 */

void list_series (const DATASET *dset, gretlopt opt, PRN *prn)
{
    int fd = gretl_function_depth();
    int debug = (opt & OPT_D);
    const char *name;
    gchar *tmp;
    int len, maxlen = 0;
    int nv = 4, n = 0;
    int i, j;

    if (dset->v == 0) {
	pprintf(prn, _("No series are defined\n"));
	return;
    }

    for (i=0; i<dset->v; i++) {
	if (show_series(i, fd, dset, debug)) {
	    tmp = NULL;
	    name = get_series_name(dset, i, &tmp, debug);
	    len = strlen(name);
	    if (len > maxlen) {
		maxlen = len;
	    }
	    g_free(tmp);
	    n++;
	}
    }

    if (maxlen < 9) {
	nv = 5;
    } else if (maxlen > 20) {
	nv = 1;
    } else if (maxlen > 13) {
	nv = 3;
    }

    pprintf(prn, _("Listing %d variables:\n"), n);

    j = 1;
    for (i=0; i<dset->v; i++) {
	if (show_series(i, fd, dset, debug)) {
	    tmp = NULL;
	    name = get_series_name(dset, i, &tmp, debug);
	    if (debug) {
		pprintf(prn, "%3d) %-*s level %d\n", i, maxlen + 2,
			name, series_get_stack_level(dset, i));
	    } else {
		pprintf(prn, "%3d) %-*s", i, maxlen + 2, name);
		if (j % nv == 0) {
		    pputc(prn, '\n');
		}
	    }
	    g_free(tmp);
	    j++;
	}
    }

    if (n % nv) {
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/**
 * maybe_list_series:
 * @dset: data information struct.
 * @prn: gretl printing struct
 *
 * Prints a list of the names of the series currently defined,
 * unless gretl messaging is turned off.
 */

void maybe_list_series (const DATASET *dset, PRN *prn)
{
    if (gretl_messages_on() && !gretl_looping() && dset->v < 51) {
	list_series(dset, OPT_NONE, prn);
    }
}

static int line_count (const char *s)
{
    int n = 1;

    while (*s) {
	if (*s == '\n') n++;
	s++;
    }

    return n;
}

static void string_print_range (const char *s, int lmin, int lmax,
				PRN *prn)
{
    int len, l = 0;

    while (*s && l < lmax) {
	len = strcspn(s, "\r\n");
	if (l >= lmin) {
	    pprintf(prn, "%.*s\n", len, s);
	}
	s += len;
	if (*s == '\r') s++;
	if (*s == '\n') s++;
	l++;
    }
}

static int print_evaled_object_range (const char *s,
                                      DATASET *dset,
                                      gretlopt opt,
                                      PRN *prn)
{
    int start, stop;
    int err;

    err = get_print_range(opt, -1, &start, &stop);

    if (!err) {
        gchar *tmp;

        if (stop < 0) {
            tmp = g_strdup_printf("%s[%d:]", s, start);
        } else {
            tmp = g_strdup_printf("%s[%d:%d]", s, start, stop);
        }
        err = generate(tmp, dset, GRETL_TYPE_NONE, OPT_P, prn);
        g_free(tmp);
    }

    return err;
}

static int print_listed_objects (const char *s,
				 const DATASET *dset,
				 gretlopt opt,
				 PRN *prn)
{
    const char *syms = "=+-/*<>?|~^!%&.,:;\\'[({$";
    char *name = NULL;
    int err = 0;

    if (strcspn(s, syms) < strlen(s)) {
	/* try treating as expression to be evaluated */
        if (have_range_option(opt)) {
            return print_evaled_object_range(s, (DATASET *) dset, opt, prn);
        } else {
            return generate(s, (DATASET *) dset, GRETL_TYPE_NONE, OPT_P, prn);
        }
    }

    while ((name = gretl_word_strdup(s, &s, OPT_S | OPT_U, &err)) != NULL) {
	user_var *uv = get_user_var_by_name(name);

	if (uv == NULL) {
	    err = E_UNKVAR;
	    break;
	} else if (have_range_option(opt)) {
	    GretlType t = user_var_get_type(uv);
	    int start, stop;

	    if (t == GRETL_TYPE_ARRAY) {
		gretl_array *a = user_var_get_value(uv);
		int len = gretl_array_get_length(a);

		err = get_print_range(opt, len, &start, &stop);
		if (!err) {
		    err = gretl_array_print_range(a, start, stop+1, prn);
		}
	    } else if (t == GRETL_TYPE_MATRIX) {
		gretl_matrix *m = user_var_get_value(uv);
		int len = gretl_matrix_rows(m);

		err = get_print_range(opt, len, &start, &stop);
		if (!err) {
		    gretl_matrix_print_range(m, name, start, stop+1, prn);
		}
	    } else if (t == GRETL_TYPE_STRING) {
		const char *s = user_var_get_value(uv);
		int len = line_count(s);

		err = get_print_range(opt, len, &start, &stop);
		if (!err) {
		    string_print_range(s, start, stop+1, prn);
		}
	    } else {
		err = print_user_var_by_name(name, dset, opt, prn);
	    }
	} else {
	    err = print_user_var_by_name(name, dset, opt, prn);
	}
	free(name);
	if (err) {
	    break;
	}
    }

    return err;
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

static int obslen_from_t (int t)
{
    char s[OBSLEN];

    sprintf(s, "%d", t + 1);
    return strlen(s);
}

/* in case we're printing a lot of data to a PRN that uses a
   buffer, pre-allocate a relatively big chunk of memory
*/

static int check_prn_size (const int *list, const DATASET *dset,
			   PRN *prn)
{
    int nx = list[0] * (dset->t2 - dset->t1 + 1);
    int err = 0;

    if (nx > 1000) {
	err = gretl_print_alloc(prn, nx * 12);
    }

    return err;
}

static int *get_pmax_array (const int *list, const DATASET *dset)
{
    int *pmax = malloc(list[0] * sizeof *pmax);
    int i, vi, T = sample_size(dset);

    if (pmax == NULL) {
	return NULL;
    }

    /* this runs fairly quickly, even for large dataset */

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	pmax[i-1] = get_signif(dset->Z[vi] + dset->t1, T);
    }

    return pmax;
}

static void bufprint_string (char *buf, const char *s,
			     int width, PRN *prn)
{
    int n = width - g_utf8_strlen(s, -1);

    *buf = '\0';

    if (n > 0) {
	bufspace(n, prn);
	pputs(prn, s);
    } else {
	gretl_utf8_strncat(buf, s, width - 3);
	n = width - g_utf8_strlen(buf, -1);
	bufspace(n-1, prn);
	pputs(prn, buf);
	pputc(prn, '~');
    }
}

static int get_range_params (int len, int *k1, int *k2)
{
    const char *s = get_optval_string(PRINT, OPT_R);
    int nf = 0;
    int err = 0;

    if (s == NULL || *s == '\0' || strchr(s, ':') == NULL) {
	err = E_PARSE;
    } else {
	char **S = gretl_string_split(s, &nf, ":");

	if (S == NULL || nf != 2) {
	    err = E_PARSE;
	} else {
	    if (S[0][0] == '\0') {
		*k1 = 1;
	    } else {
		*k1 = gretl_int_from_string(S[0], &err);
	    }
	    if (!err) {
		if (S[1][0] == '\0') {
		    *k2 = len;
		} else {
		    *k2 = gretl_int_from_string(S[1], &err);
		}
	    }
	}
        strings_array_free(S, nf);
    }

    return err;
}

static int get_head_tail_params (gretlopt opt, int len,
                                 int *k1, int *k2)
{
    int head = (opt & OPT_H)? 1 : 0;
    const char *s;
    int err = 0;

    s = get_optval_string(PRINT, head ? OPT_H : OPT_T);
    if (s == NULL) {
        *k1 = head ? 1 : -5;
        *k2 = head ? 5 : len;
    } else {
        int htk = positive_int_from_string(s);

        if (htk < 0) {
            err = E_INVARG;
        } else if (head) {
            *k1 = 1;
            *k2 = htk;
        } else if (!err) {
            *k1 = -htk;
            *k2 = len;
        }
    }

    return err;
}

static int get_print_range (gretlopt opt, int len,
                            int *start, int *stop)
{
    int k1 = 0, k2 = 0;
    int err = 0;

    if (opt & OPT_R) {
        err = get_range_params(len, &k1, &k2);
    } else {
        err = get_head_tail_params(opt, len, &k1, &k2);
    }

    if (err) {
        return err;
    }

    if (len == -1) {
        *start = k1;
        *stop = k2;
    } else if (k1 == 0 || k2 == 0) {
        fprintf(stderr, "get_print_range: got a zero value\n");
        err = E_INVARG;
    } else if (k1 < 0 || k2 < 0) {
        if (k1 < 0) {
            k1 = len + k1 + 1;
        }
        if (k2 < 0) {
            k2 = len + k2 + 1;
        }
    }
    if (!err && k2 < k1) {
        fprintf(stderr, "get_print_range: got empty range\n");
    }
    if (!err && (k1 < 1 || k1 > len || k2 < 1 || k2 > len)) {
        fprintf(stderr, "get_print_range: out of bounds\n");
        err = E_INVARG;
    }

    if (!err) {
        /* convert to zero-based */
        *start = k1 - 1;
        *stop = k2 - 1;
    }

    return err;
}

static int *column_widths_from_list (const int *list,
				     const DATASET *dset)
{
    int *ret = gretl_list_new(list[0]);
    int i, w, maxw = 0;

    for (i=1; i<=list[0]; i++) {
	w = 1 + series_get_string_width(dset, list[i]);
	if (w > maxw) {
	    if (list[0] > 1 && w > 32) {
		w = 32;
	    }
	    maxw = w;
	}
	/* check this minimum? */
	ret[i] = w < 13 ? 13 : w;
    }

    if (maxw <= 13) {
	/* standardize on moderate width */
	for (i=1; i<=list[0]; i++) {
	    ret[i] = 13;
	}
    }

    return ret;
}

static void print_plain_numbers (int *list, const DATASET *dset,
				 PRN *prn)
{
    int i, vi, t;

    for (t=dset->t1; t<=dset->t2; t++) {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (na(dset->Z[vi][t])) {
		pputs(prn, "NA");
	    } else {
		pprintf(prn, "%.8g", dset->Z[vi][t]);
	    }
	    if (i < list[0]) {
		pputc(prn, ' ');
	    } else {
		pputc(prn, '\n');
	    }
	}
    }
}

/* print the series referenced in @list by observation */

static int print_by_obs (int *list, const DATASET *dset,
			 gretlopt opt, int screenvar,
			 PRN *prn)
{
    int BMAX = libset_get_int(DATACOLS);
    int i, j, j0, k, t, nrem;
    int *colwidths = NULL;
    int obslen = 0;
    int *pmax = NULL;
    char buf[128];
    int blist[BMAX+1];
    int wlist[BMAX+1];
    int gprec = 6;
    int vi, wi;
    double x;
    int err = 0;

    if (!(opt & OPT_S)) {
	pmax = get_pmax_array(list, dset);
	if (pmax == NULL) {
	    return E_ALLOC;
	}
    }

    if (opt & OPT_D) {
	obslen = obslen_from_t(dset->t2);
    } else {
	obslen = max_obs_marker_length(dset);
    }

    colwidths = column_widths_from_list(list, dset);

    nrem = list[0];
    k = 1;

    while (nrem > 0) {
	/* fill the "block" list */
	j0 = k;
	wlist[0] = blist[0] = 0;
	for (i=1; i<=BMAX && nrem>0; i++) {
	    blist[i] = list[k];
	    blist[0] += 1;
	    wlist[i] = colwidths[k];
	    wlist[0] += 1;
	    k++;
	    nrem--;
	}

	varheading(blist, obslen, wlist, dset, 0, prn);

	for (t=dset->t1; t<=dset->t2; t++) {
	    if (screenvar && dset->Z[screenvar][t] == 0.0) {
		/* screened out by boolean */
		continue;
	    }
	    if (opt & OPT_D) {
		pprintf(prn, "%*d", obslen, t + 1);
	    } else {
		real_print_obs_marker(t, dset, obslen, 0, prn);
	    }
	    for (i=1, j=j0; i<=blist[0]; i++, j++) {
		vi = blist[i];
		wi = wlist[i];
		if (!(opt & OPT_U) && is_string_valued(dset, vi)) {
		    const char *s = series_get_string_for_obs(dset, vi, t);

		    if (s == NULL || *s == '\0') {
			bufspace(wi, prn);
		    } else {
			bufprint_string(buf, s, wi, prn);
		    }
		} else {
		    x = dset->Z[vi][t];
		    if (na(x)) {
			bufspace(wi, prn);
		    } else {
			bufprintnum(buf, x, pmax[j-1], gprec, wi);
			pputs(prn, buf);
		    }
		}
	    }
	    pputc(prn, '\n');
	}
    }

    pputc(prn, '\n');

    free(pmax);
    free(colwidths);

    return err;
}

static int print_by_var (const int *list, const DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    const char *name;
    gchar *tmp;
    int i, vi;

    pputc(prn, '\n');

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (vi < 0 || vi > dset->v) {
	    continue;
	}
	if (list[0] > 1) {
	    tmp = NULL;
	    name = get_series_name(dset, vi, &tmp, 0);
	    pprintf(prn, "%s:\n", name);
	    g_free(tmp);
	}
	print_var_smpl(vi, dset, prn);
	pputc(prn, '\n');
	if (!(opt & OPT_U) && is_string_valued(dset, vi)) {
	    print_stringvals_for_var(dset, vi, prn);
	} else {
	    print_series_by_var(dset, vi, prn);
	}
	pputc(prn, '\n');
    }

    return 0;
}

static int midas_print_list (const int *list,
			     const DATASET *dset,
			     PRN *prn)
{
    DATASET *tmpset = NULL;
    int err = 0;

    tmpset = midas_aux_dataset(list, dset, &err);

    if (!err) {
	int mlist[2] = {1, 0};

	err = print_by_obs(mlist, tmpset, OPT_NONE, 0, prn);
	destroy_dataset(tmpset);
    }

    return err;
}

/**
 * printdata:
 * @list: list of variables to print.
 * @ostr: optional string holding names of non-series objects to print.
 * @dset: dataset struct.
 * @opt: if OPT_O, print the data by observation (series in columns);
 * if OPT_D, use simple obs numbers, not dates; if OPT_M, print midas
 * list in original time-series order; if OPT_R print specified range
 * of object; if OPT_X (relevant only for series), print the data
 * by observation without any header or observation info.
 * @prn: gretl printing struct.
 *
 * Print the data for the variables in @list over the currently
 * defined sample range, or the objects named in @ostr.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int printdata (const int *list, const char *ostr,
	       DATASET *dset, gretlopt opt,
	       PRN *prn)
{
    int screenvar = 0;
    int *plist = NULL;
    int err = 0;

    err = incompatible_options(opt, OPT_R | OPT_H | OPT_T);
    if (err) {
        return err;
    }

    if (list != NULL && list[0] == 0) {
	/* an explicitly empty list was given */
	if (ostr == NULL) {
	    return 0; /* no-op */
	} else {
	    goto endprint;
	}
    } else if (list == NULL) {
	/* no list was given */
	if (ostr == NULL && dset != NULL) {
	    int nvars = 0;

	    plist = full_var_list(dset, &nvars);
	    if (nvars == 0) {
		/* no-op */
		if (gretl_messages_on()) {
		    pputs(prn, "no data\n");
		}
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

    if (opt & OPT_M) {
	err = midas_print_list(plist, dset, prn);
	free(plist);
	return err;
    }

    if (gretl_list_has_separator(plist)) {
	err = adjust_print_list(plist, &screenvar, opt);
	if (err) {
	    free(plist);
	    return err;
	}
    }

    if (plist[0] == 0) {
	/* no series */
	pputc(prn, '\n');
	goto endprint;
    }

    if (gretl_print_has_buffer(prn)) {
	err = check_prn_size(plist, dset, prn);
	if (err) {
	    goto endprint;
	}
    }

    if (have_range_option(opt)) {
	/* --range, --head or --tail */
	int save_t1 = dset->t1;
	int save_t2 = dset->t2;
	int start = 0, stop = 0;

        /* --byobs is implicit */
        opt |= OPT_O;

	err = get_print_range(opt, sample_size(dset), &start, &stop);
	if (err) {
	    return err;
	} else if (stop < start) {
	    goto endprint;
	}
	dset->t1 = save_t1 + start;
	dset->t2 = save_t1 + stop;
	if (opt & OPT_X) {
	    print_plain_numbers(plist, dset, prn);
	} else if (opt & OPT_O) {
	    err = print_by_obs(plist, dset, opt, screenvar, prn);
	} else {
	    err = print_by_var(plist, dset, opt, prn);
	}
	dset->t1 = save_t1;
	dset->t2 = save_t2;
    } else {
	if (opt & OPT_X) {
	    print_plain_numbers(plist, dset, prn);
	} else if (opt & OPT_O) {
	    err = print_by_obs(plist, dset, opt, screenvar, prn);
	} else {
	    err = print_by_var(plist, dset, opt, prn);
	}
    }

 endprint:

    if (!err && ostr != NULL) {
	err = print_listed_objects(ostr, dset, opt, prn);
    }

    free(plist);

    return err;
}

int print_series_with_format (const int *list,
			      const DATASET *dset,
			      char fmt, int digits,
			      PRN *prn)
{
    int BMAX = libset_get_int(DATACOLS);
    int i, j, j0, v, t, k, nrem = 0;
    int *colwidths, blist[BMAX+1];
    char obslabel[OBSLEN];
    char format[16];
    char *buf = NULL;
    int buflen, obslen;
    double x;
    int err = 0;

    if (list == NULL || list[0] == 0) {
	return 0;
    }

    for (i=1; i<=list[0]; i++) {
	if (list[i] >= dset->v) {
	    return E_DATA;
	}
    }

    colwidths = gretl_list_new(list[0]);
    if (colwidths == NULL) {
	return E_ALLOC;
    }

    nrem = list[0];

    buflen = 0;
    for (i=1; i<=list[0]; i++) {
	colwidths[i] = series_column_width(list[i], dset, fmt, digits);
	colwidths[i] += 3;
	if (colwidths[i] > buflen) {
	    buflen = colwidths[i];
	}
    }

    buf = malloc(buflen);
    if (buf == NULL) {
	free(colwidths);
	return E_ALLOC;
    }

    if (gretl_print_has_buffer(prn)) {
	err = check_prn_size(list, dset, prn);
	if (err) {
	    goto bailout;
	}
    }

    if (fmt == 'f') {
	sprintf(format, "%%.%df", digits);
    } else {
	sprintf(format, "%%#.%dg", digits);
    }

    obslen = max_obs_marker_length(dset);

    k = 1;

    while (nrem > 0) {
	/* fill the "block" list */
	j0 = k;
	blist[0] = 0;
	for (i=1; i<=BMAX && nrem>0; i++) {
	    blist[i] = list[k++];
	    blist[0] += 1;
	    nrem--;
	}

	/* print block heading */
	bufspace(obslen, prn);
	for (i=1, j=j0; i<=blist[0]; i++, j++) {
	    v = blist[i];
	    pprintf(prn, "%*s", colwidths[j], dset->varname[v]);
	}
	pputs(prn, "\n\n");

	/* print block observations */
	for (t=dset->t1; t<=dset->t2; t++) {
	    get_obs_string(obslabel, t, dset);
	    pprintf(prn, "%*s", obslen, obslabel);
	    for (i=1, j=j0; i<=blist[0]; i++, j++) {
		v = blist[i];
		x = dset->Z[v][t];
		if (na(x)) {
		    bufspace(colwidths[j], prn);
		} else {
		    sprintf(buf, format, x);
		    if (fmt == 'g') {
			/* post-process ugliness */
			cut_trailing_point(cut_extra_zero(buf, digits));
		    }
		    pprintf(prn, "%*s", colwidths[j], buf);
		}
	    }
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    }

 bailout:

    free(colwidths);
    free(buf);

    return err;
}

enum {
    RTF_HEADER,
    RTF_TRAILER
};

static void rtf_print_row_spec (int ncols, int type, PRN *prn)
{
    int j;

    if (type == RTF_TRAILER) {
	pputc(prn, '{');
    }

    pputs(prn, "\\trowd\\trautofit1\n\\intbl\n");
    for (j=1; j<=ncols; j++) {
	pprintf(prn, "\\cellx%d\n", j);
    }

    if (type == RTF_TRAILER) {
	pputs(prn, "\\row }\n");
    }
}

/**
 * print_data_in_columns:
 * @list: list of variables to print.
 * @obsvec: list of observation numbers (or %NULL)
 * @dset: dataset struct.
 * @opt: may include OPT_X to exclude the observations
 * column that is usually printed first.
 * @prn: gretl printing struct.
 *
 * Print the data for the variables in @list.  If @obsvec is not %NULL,
 * it should specify a sort order; the first element of @obsvec must
 * contain the number of observations that follow.  By default, printing is
 * plain text, formatted in columns using space characters, but if the @prn
 * format is set to %GRETL_FORMAT_CSV then printing respects the user's
 * choice of column delimiter, and if the format is set to %GRETL_FORMAT_RTF
 * the data are printed as an RTF table.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int print_data_in_columns (const int *list, const int *obsvec,
			   const DATASET *dset, gretlopt opt,
			   PRN *prn)
{
    int csv = csv_format(prn);
    int rtf = rtf_format(prn);
    const char *csv_na = "";
    int print_obs = 1;
    char delim = 0;
    int *pmax = NULL;
    int *colwidths = NULL;
    double xx;
    char obs_string[OBSLEN];
    char buf[128];
    int ncols = 0, obslen = 0;
    int gprec = 6;
    int i, s, t, T;

    if (obsvec != NULL) {
	T = obsvec[0];
    } else {
	T = sample_size(dset);
    }

    /* we must have a non-empty list of variables */
    if (list == NULL || list[0] < 1) {
	return E_DATA;
    }

    /* ...with no bad variable numbers */
    for (i=1; i<=list[0]; i++) {
	if (list[i] < 0 || list[i] >= dset->v) {
	    return E_DATA;
	}
    }

    /* and T must be in bounds */
    if (T > dset->n - dset->t1) {
	return E_DATA;
    }

    pmax = get_pmax_array(list, dset);
    if (pmax == NULL) {
	return E_ALLOC;
    }

    if (csv) {
	/* columns delimited by some character */
	gprec = 15;
	delim = get_data_export_delimiter();
	if (get_local_decpoint() == ',' && delim == ',') {
	    delim = ';';
	}
	csv_na = get_csv_na_write_string();
	if (opt & OPT_X) {
	    print_obs = 0;
	    obslen = -1;
	}
    } else if (rtf) {
	ncols = list[0] + 1;
    } else {
	colwidths = column_widths_from_list(list, dset);
	obslen = max_obs_marker_length(dset);
    }

    if (rtf) {
	pputs(prn, "{\\rtf1\n");
	rtf_print_row_spec(ncols, RTF_HEADER, prn);
    }

    varheading(list, obslen, colwidths, dset, delim, prn);

    if (rtf) {
	rtf_print_row_spec(ncols, RTF_TRAILER, prn);
    }

    /* print data by observations */
    for (s=0; s<T; s++) {
	t = (obsvec != NULL)? obsvec[s+1] : (dset->t1 + s);
	if (t >= dset->n) {
	    continue;
	}

	if (rtf) {
	    rtf_print_row_spec(ncols, RTF_HEADER, prn);
	    pputc(prn, '{');
	}

	if (print_obs) {
	    get_obs_string(obs_string, t, dset);
	    if (csv) {
		pprintf(prn, "%s%c", obs_string, delim);
	    } else if (rtf) {
		pprintf(prn, "%s\\cell ", obs_string);
	    } else {
		pprintf(prn, "%*s", obslen, obs_string);
	    }
	}

	for (i=1; i<=list[0]; i++) {
	    const char *strval = NULL;
	    int wi = 0;

	    if (colwidths != NULL) {
		wi = colwidths[i];
	    }

	    if (is_string_valued(dset, list[i])) {
		strval = series_get_string_for_obs(dset, list[i], t);
	    }

	    if (strval != NULL) {
		/* display string value */
		if (*strval == '\0') {
		    if (csv) {
			pputs(prn, "\"\"");
		    } else if (rtf) {
			pputs(prn, "\\qr \\cell ");
		    } else {
			bufspace(wi, prn);
		    }
		} else if (csv) {
		    pprintf(prn, "\"%s\"", strval);
		} else if (rtf) {
		    pprintf(prn, "\\qr %s\\cell ", strval);
		} else {
		    bufprint_string(buf, strval, wi, prn);
		}
	    } else {
		/* numerical value */
		xx = dset->Z[list[i]][t];
		if (na(xx)) {
		    if (csv) {
			pputs(prn, csv_na);
		    } else if (rtf) {
			pputs(prn, "\\qr NA\\cell ");
		    } else {
			bufspace(wi, prn);
		    }
		} else {
		    if (rtf) {
			bufprintnum(buf, xx, pmax[i-1], gprec, 0);
			pprintf(prn, "\\qr %s\\cell ", buf);
		    } else {
			bufprintnum(buf, xx, pmax[i-1], gprec, wi);
			pputs(prn, buf);
		    }
		}
	    }
	    if (csv && i < list[0]) {
		pputc(prn, delim);
	    }
	}
	if (rtf) {
	    pputs(prn, "}\n");
	    rtf_print_row_spec(ncols, RTF_TRAILER, prn);
	} else {
	    pputc(prn, '\n');
	}
    }

    if (rtf) {
	pputs(prn, "}\n");
    } else if (!csv) {
	pputc(prn, '\n');
    }

    free(pmax);
    free(colwidths);

    return 0;
}

int print_fcast_stats_matrix (const gretl_matrix *m,
			      int T, gretlopt opt,
			      PRN *prn)
{
    const char *strs[] = {
	N_("Mean Error"),
	N_("Root Mean Squared Error"),
	N_("Mean Absolute Error"),
	N_("Mean Percentage Error"),
	N_("Mean Absolute Percentage Error"),
	N_("Theil's U1"),
	N_("Bias proportion, UM"),
	N_("Regression proportion, UR"),
	N_("Disturbance proportion, UD")
    };
    const char *U2_str = N_("Theil's U2");
    double x;
    int i, n, nmax = 0;
    int len, err = 0;

    len = gretl_vector_get_length(m);

    for (i=0; i<len; i++) {
	x = gretl_vector_get(m, i);
	if (!isnan(x)) {
	    n = g_utf8_strlen(_(strs[i]), -1);
	    if (n > nmax) {
		nmax = n;
	    }
	}
    }

    nmax += 2;

    pputs(prn, "  ");
    pputs(prn, _("Forecast evaluation statistics"));
    pputc(prn, ' ');
    pprintf(prn, _("using %d observations"), T);
    pputs(prn, "\n\n");

    for (i=0; i<len; i++) {
	const char *si;

	x = gretl_vector_get(m, i);
	if (!isnan(x)) {
	    si = (i == 5 && (opt & OPT_T))? U2_str : strs[i];
	    pprintf(prn, "  %-*s % .5g\n", UTF_WIDTH(_(si), nmax), _(si), x);
	}
    }
    pputc(prn, '\n');

    return err;
}

static int fr_print_fc_stats (const FITRESID *fr, gretlopt opt,
			      PRN *prn)
{
    gretl_matrix *m;
    int t1, t2, n_used;
    int err = 0;

    fcast_get_continuous_range(fr, &t1, &t2);

    if (t2 - t1 + 1 <= 0) {
	return E_MISSDATA;
    }

    m = forecast_stats(fr->actual, fr->fitted, t1, t2, &n_used,
		       opt, &err);

    if (!err) {
	err = print_fcast_stats_matrix(m, n_used, opt, prn);
    }

    gretl_matrix_free(m);

    return err;
}

#define SIGMA_MIN 1.0e-18

int text_print_fit_resid (const FITRESID *fr,
			  const DATASET *dset,
			  PRN *prn)
{
    gretlopt fc_opt = OPT_NONE;
    int kstep = fr->method == FC_KSTEP;
    int t, anyast = 0;
    double yt, yf, et;
    int ywidth;
    int obslen;
    int err = 0;

    obslen = max_obs_marker_length(dset);
    ywidth = fit_resid_head(fr, dset, obslen, prn);

    for (t=fr->t1; t<=fr->t2; t++) {
	real_print_obs_marker(t, dset, obslen, 0, prn);

	yt = fr->actual[t];
	yf = fr->fitted[t];
	et = fr->resid[t];

	if (na(yt)) {
	    pputc(prn, '\n');
	} else if (na(yf)) {
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%*.*f\n", ywidth, fr->pmax, yt);
	    } else {
		pprintf(prn, "%#*g\n", ywidth, yt);
	    }
	} else if (na(et)) {
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%*.*f", ywidth, fr->pmax, yt);
		pprintf(prn, "%13.*f", fr->pmax, yf);
	    } else {
		pprintf(prn, "%#*g\n", ywidth, yt);
		pprintf(prn, "%#13g\n", yf);
	    }
	} else {
	    int ast = 0;

	    if (!kstep && fr->sigma > SIGMA_MIN) {
		ast = (fabs(et) > 2.5 * fr->sigma);
		if (ast) {
		    anyast = 1;
		}
	    }
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%*.*f%13.*f%13.*f%s\n", ywidth,
			fr->pmax, yt, fr->pmax, yf, fr->pmax, et,
			(ast)? " *" : "");
	    } else {
		pprintf(prn, "%#*g%#13g%#13g%s\n",
			ywidth, yt, yf, et,
			(ast)? " *" : "");
	    }
	}
    }

    pputc(prn, '\n');

    if (anyast) {
	pputs(prn, _("Note: * denotes a residual in excess of "
		     "2.5 standard errors\n"));
    }

    if (dataset_is_time_series(dset)) {
	fc_opt |= OPT_T;
    }
    fr_print_fc_stats(fr, fc_opt, prn);

    if (kstep && fr->nobs > 0 && gretl_in_gui_mode()) {
	err = plot_fcast_errs(fr, NULL, dset, OPT_NONE);
    }

    return err;
}

/**
 * text_print_forecast:
 * @fr: pointer to structure containing forecasts.
 * @dset: dataset information.
 * @opt: if includes %OPT_P, make a plot of the forecasts;
 * if includes %OPT_N, skip printing of the forecast
 * evaluation statistics.
 * @prn: printing struct.
 *
 * Prints the forecasts in @fr to @prn, and also plots the
 * forecasts if %OPT_P is given. If a plot is requested and
 * @fr includes forecast standard errors, then the options
 * %OPT_F or %OPT_L may be given to use "fill" style or lines,
 * respectively, for the confidence bands (the default style
 * being vertical bars per observation).
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int text_print_forecast (const FITRESID *fr, DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    int do_errs = (fr->sderr != NULL);
    int obslen, pmax = fr->pmax;
    int errpmax = fr->pmax;
    int quiet = (opt & OPT_Q);
    int ywidth;
    double *maxerr = NULL;
    double conf = 100 * (1 - fr->alpha);
    double tval = 0;
    char label[32];
    int t, err = 0;

    if (opt & OPT_T) {
	/* --stats-only */
	return fr_print_fc_stats(fr, OPT_D, prn);
    }

    if (do_errs) {
	maxerr = malloc(fr->nobs * sizeof *maxerr);
	if (maxerr == NULL) {
	    return E_ALLOC;
	}
    }

    if (!quiet) {
	pputc(prn, '\n');
    }

    if (do_errs) {
	double a2 = fr->alpha / 2;

	tval = (fr->asymp)? normal_critval(a2) : student_critval(fr->df, a2);

	if (!quiet) {
	    if (fr->asymp) {
		pprintf(prn, _(" For %g%% confidence intervals, z(%g) = %.2f\n"),
			conf, a2, tval);
	    } else {
		pprintf(prn, _(" For %g%% confidence intervals, t(%d, %g) = %.3f\n"),
			conf, fr->df, a2, tval);
	    }
	}
    }

    obslen = max_obs_marker_length(dset);

    if (!quiet) {
	pputc(prn, '\n');
    }

    bufspace(obslen + 1, prn);

    maybe_trim_varname(label, fr->depvar);
    ywidth = strlen(label) + 1;
    if (ywidth < 12) {
	ywidth = 12;
    }
    pprintf(prn, "%*s", ywidth, label);

    pprintf(prn, "%*s", UTF_WIDTH(_("prediction"), 14), _("prediction"));

    if (do_errs) {
	pprintf(prn, "%*s", UTF_WIDTH(_(" std. error"), 14), _(" std. error"));
	pprintf(prn, _("        %g%% interval\n"), conf);
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

    for (t=fr->t0; t<=fr->t2; t++) {
	print_obs_marker(t, dset, obslen, prn);
	fcast_print_x(fr->actual[t], ywidth + 2, pmax, prn);

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
		maxerr[t] = tval * fr->sderr[t];
		fcast_print_x(fr->fitted[t] - maxerr[t], 15, pmax, prn);
		pputs(prn, " - ");
		fcast_print_x(fr->fitted[t] + maxerr[t], 10, pmax, prn);
	    }
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    if (!(opt & OPT_N)) {
	gretlopt fc_opt = OPT_D;

	if (dataset_is_time_series(dset)) {
	    fc_opt |= OPT_T;
	}
	fr_print_fc_stats(fr, fc_opt, prn);
    }

    /* do we really want a plot for non-time series? */

    if ((opt & OPT_P) && fr->nobs > 0) {
	err = plot_fcast_errs(fr, maxerr, dset, opt);
    }

    if (maxerr != NULL) {
	free(maxerr);
    }

    return err;
}

/**
 * print_fit_resid:
 * @pmod: pointer to gretl model.
 * @dset: dataset struct.
 * @prn: gretl printing struct.
 *
 * Print to @prn the fitted values and residuals from @pmod.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int print_fit_resid (const MODEL *pmod, const DATASET *dset,
		     PRN *prn)
{
    FITRESID *fr;
    int err = 0;

    fr = get_fit_resid(pmod, dset, &err);

    if (!err) {
	text_print_fit_resid(fr, dset, prn);
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
	N_("loglikelihood"),
	N_("GMM criterion"),
	N_("SSR"),
	N_("Criterion"),
    };
    const char *cstr = cstrs[type];
    double x;
    int details;
    int sldone = 0;
    int i;

    details = libset_get_int(MAX_VERBOSE) == 2;

    if (type == C_GMM) {
	crit = -crit;
    } else if (type == C_SSR && crit < 0) {
	crit = -crit;
    }

    if (iter < 0) {
	pprintf(prn, "--- %s:\n", _("FINAL VALUES"));
    } else if (details) {
	pprintf(prn, "%s %4d: ", _("Iteration"), iter);
    } else {
	if (iter == 1) {
	    gretl_print_ensure_vspace(prn);
	    pprintf(prn, "--- %s\n", _("Iteration"));
	}
	pprintf(prn, "%4d: ", iter);
    }

    if (na(crit) || na(-crit)) {
	pprintf(prn, "%s = NA", _(cstr));
    } else if (details) {
	pprintf(prn, "%s = %#.12g", _(cstr), crit);
    } else {
	pprintf(prn, "%s %#.12g", _(cstr), crit);
    }

    if (sl > 0.0 && !na(sl)) {
	if (details || g == NULL) {
	    pprintf(prn, _(" (steplength = %g)"), sl);
	} else {
	    pprintf(prn, _(" (step %g"), sl);
	    sldone = 1;
	}
    }

    if (details) {
	pputc(prn, '\n');
	if (b != NULL) {
	    pputs(prn, _("Parameters: "));
	    for (i=0; i<k; i++) {
		print_iter_val(b[i], i, k, prn);
	    }
	    pputc(prn, '\n');
	}
    }

    if (g != NULL) {
	if (details) {
	    pputs(prn, _("Gradients:  "));
	}
	x = 0.0;
	for (i=0; i<k; i++) {
	    x += fabs(b[i] * g[i]);
	    if (details) {
		print_iter_val(g[i], i, k, prn);
	    }
	}
	if (details) {
	    pprintf(prn, " (%s %.2e)", _("norm"), sqrt(x/k));
	} else if (sldone) {
	    pprintf(prn, ", %s %.2e)", _("norm"), sqrt(x/k));
	} else {
	    pprintf(prn, " (%s %.2e)", _("norm"), sqrt(x/k));
	}
	if (details || iter == -1) {
	    pputc(prn, '\n');
	}
    }

    if (iter >= 0) {
	pputc(prn, '\n');
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

typedef struct readbuf_ readbuf;

struct readbuf_ {
    const char *start;
    const char *point;
};

static readbuf *rbuf;
static int n_bufs;

static readbuf *matching_buffer (const char *s)
{
    int i;

    for (i=0; i<n_bufs; i++) {
	if (rbuf[i].start == s) {
	    return &rbuf[i];
	}
    }

    return NULL;
}

/**
 * bufgets_init:
 * @buf: source buffer.
 *
 * Initializes a text buffer for use with bufgets().
 *
 * Returns: 0 on success, non-zero on error.
 */

int bufgets_init (const char *buf)
{
    readbuf *tmp = matching_buffer(buf);
    int i, err = 0;

    if (buf == NULL) {
	fprintf(stderr, "bufgets_init: got NULL argument\n");
	return 1;
    } else if (tmp != NULL) {
	fprintf(stderr, "GRETL ERROR: buffer at %p is already "
		"initialized\n", (void *) buf);
	return 1;
    }

    for (i=0; i<n_bufs; i++) {
	if (rbuf[i].start == NULL) {
	    /* OK, re-use an existing slot */
	    rbuf[i].start = rbuf[i].point = buf;
	    return 0;
	}
    }

    tmp = realloc(rbuf, (n_bufs + 1) * sizeof *tmp);

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	rbuf = tmp;
	rbuf[n_bufs].start = rbuf[n_bufs].point = buf;
	n_bufs++;
    }

    return err;
}

int query_bufgets_init (const char *buf)
{
    if (matching_buffer(buf) != NULL) {
	return 0; /* OK */
    } else {
	return bufgets_init(buf);
    }
}

static const char *rbuf_get_point (const char *s)
{
    readbuf *rbuf = matching_buffer(s);

    return (rbuf == NULL)? NULL : rbuf->point;
}

static void rbuf_set_point (const char *s, const char *p)
{
    readbuf *rbuf = matching_buffer(s);

    if (rbuf != NULL) {
	rbuf->point = p;
    }
}

/**
 * bufgets_finalize:
 * @buf: source buffer.
 *
 * Signals that we are done reading from @buf.
 */

void bufgets_finalize (const char *buf)
{
    readbuf *rbuf = matching_buffer(buf);

    if (rbuf != NULL) {
	rbuf->start = rbuf->point = NULL;
    }
}

enum {
    GOT_LF = 1,
    GOT_CR,
    GOT_CRLF
};

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
    int i, status = 0;
    const char *p;

    p = rbuf_get_point(buf);
    if (p == NULL) {
	return NULL;
    }

    if (*p == '\0') {
	/* reached the end of the buffer */
	return NULL;
    }

    *s = '\0';

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
	    fprintf(stderr, "*** bufgets: line too long: max %d characters\n",
		    (int) size);
	    s[i] = '\0';
	    fprintf(stderr, " '%.32s...'\n", s);
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

char *safe_bufgets (char **pdest, size_t *psize, const char *buf)
{
    char *s = *pdest;
    size_t sz = *psize;
    int i, status = 0;
    const char *p;

    p = rbuf_get_point(buf);
    if (p == NULL) {
	return NULL;
    }

    if (*p == '\0') {
	/* reached the end of the buffer */
	return NULL;
    }

    *s = '\0';

    /* advance to line-end, end of buffer, or maximum size,
       whichever comes first */
    for (i=0; ; i++) {
	if (p[i] == '\0') {
	    break;
	}
	if (p[i] == '\r') {
	    if (p[i+1] == '\n') {
		status = GOT_CRLF;
	    } else {
		status = GOT_CR;
	    }
	    break;
	}
	if (p[i] == '\n') {
	    status = GOT_LF;
	    break;
	}
    }

    if (i < sz - 1) {
	/* dest is big enough */
	memcpy(s, p, i);
	s[i] = '\0';
    } else {
	/* dest is too small */
	size_t new_sz = 2*sz;
	char *tmp;

	if (i < new_sz - 1) {
	    tmp = realloc(*pdest, new_sz);
	} else {
	    new_sz = i + 1;
	    tmp = realloc(*pdest, new_sz);
	}
	if (tmp == NULL) {
	    fputs("*** bufgets: line too long\n", stderr);
	    fprintf(stderr, " '%.32s...'\n", s);
	    return NULL;
	} else {
	    s = *pdest = tmp;
	    *psize = new_sz;
	    memcpy(s, p, i);
	    s[i] = '\0';
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
    if (status) {
	strcat(s, "\n");
    }

    rbuf_set_point(buf, p);

    return s;
}

size_t bufgets_peek_line_length (const char *buf)
{
    const char *p = rbuf_get_point(buf);
    size_t len = 0;

    if (p == NULL || *p == '\0') {
	return 0;
    }

    while (*p) {
	if (*p == '\r' || *p == '\n') {
	    break;
	} else {
	    len++;
	    p++;
	}
    }

    return len + 1;
}

/**
 * bufseek:
 * @buf: char buffer.
 * @offset: offset from start of @buf.
 *
 * Buffer equivalent of fseek() with SEEK_SET.  Note that @buf
 * must first be initialized via bufgets_init().
 *
 * Returns: 0 on success, 1 on error.
 */

int bufseek (const char *buf, long int offset)
{
    readbuf *rbuf = matching_buffer(buf);

    if (rbuf != NULL) {
	rbuf->point = rbuf->start + offset;
	return 0;
    }

    return 1;
}

/**
 * buf_rewind:
 * @buf: char buffer.
 *
 * Buffer equivalent of rewind().  Note that @buf
 * must first be initialized using bufgets_init().
 */

void buf_rewind (const char *buf)
{
    bufseek(buf, 0);
}

/**
 * buf_back_lines:
 * @buf: char buffer.
 * @n: number of lines.
 *
 * Applies only when @buf has been initialized with bufgets_init().
 * Moves the reading point of @buf for bufgets() back by @n lines,
 * if possible, or to the start of the buffer if @n is greater than
 * the number of previous lines.
 *
 * Returns: 0 on success, 1 on error.
 */

int buf_back_lines (const char *buf, int n)
{
    readbuf *rbuf = matching_buffer(buf);

    if (rbuf != NULL) {
	const char *p = rbuf->point;
	int i, len = p - rbuf->start;
	int count = 0;

	for (i=0; i<len; i++) {
	    p--;
	    if (*p == '\n') {
		count++;
	    }
	    if (count == n + 1) {
		p++;
		break;
	    }
	}
	rbuf->point = p;
	return 0;
    }

    return 1;
}

/**
 * buftell:
 * @buf: char buffer.
 *
 * Buffer equivalent of ftell.  Note that @buf
 * must first be initialized via bufgets_init().
 *
 * Returns: offset from start of buffer.
 */

long buftell (const char *buf)
{
    readbuf *rbuf = matching_buffer(buf);

    return (rbuf == NULL)? 0 : rbuf->point - rbuf->start;
}

/* for internal use */

void bufgets_cleanup (void)
{
    if (n_bufs > 0) {
	free(rbuf);
	rbuf = NULL;
	n_bufs = 0;
    }
}
