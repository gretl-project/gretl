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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include "monte_carlo.h"
#include "gretl_func.h"
#include "objstack.h"
#include "cmd_private.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_panel.h"
#include "loop_private.h"
#include "texprint.h"
#include "gretl_string_table.h"

#include <errno.h>

#ifndef WIN32
# include <signal.h>
# ifdef USE_GLIB2
#  include <glib.h>
#  define USE_GSPAWN
# endif
#endif

/**
 * date:
 * @nt: observation number (zero-based).
 * @pd: data periodicity or frequency.
 * @sd0: floating point representation of starting date.
 *
 * Returns: the date corresponding to @nt, as a double-precision number.
 */

double date (int nt, int pd, const double sd0)
{
    int ysd = (int) sd0, yy, pp, yp;
    int p10 = 10;

    if (pd == 1) {
	return (double) (ysd + nt);  
    } 

    pp = pd;
    while ((pp = pp / 10)) {
	p10 *= 10;
    }

    pp = nt % pd + p10 * (sd0 - ysd) + .5;
    if (pp != pd)  {
        yy = ysd + nt/pd  + pp/pd + .5;
        yp = pp % pd;
    }  else {
        yy = ysd + nt/pd + .5;
        yp = pp;
    }

    return yy + (double) yp / p10;
}

/**
 * ijton:
 * @i: row number (0-based)
 * @j: column number (0-based)
 * @nrows: number of rows (and columns) in symmetric matrix.
 *
 * Given a (row, column) reference into a symmetric 2-dimensional 
 * matrix A, finds the index into a 1-dimensional array x
 * composed of the non-redundant (lower) elements of A.
 *
 * E.g. for the 3 x 3 case with 6 non-redundant elements, 0 to 5,
 *
 *    A(0,0) = x[0]  A(0,1) = x[1]  A(0,2) = x[2]
 *    A(1,0) = x[1]  A(1,1) = x[3]  A(1,2) = x[4]
 *    A(2,0) = x[2]  A(2,1) = x[4]  A(2,1) = x[5]
 *
 * Returns: 0-based index into flat array.
 */

int ijton (int i, int j, int nrows)
{
    if (i > j) {
	int tmp = i;

	i = j;
	j = tmp;
    }

    return nrows * i + j - i - ((i - 1) * i / 2);
}

/**
 * ztox:
 * @i: index number of variable to extract.
 * @px: array into which to write the series.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 * 
 * Pull one series from data matrix and put it into @px.
 *
 * Returns: the number of valid observations put into @px.
 */

int ztox (int i, double *px, const double **Z, const DATAINFO *pdinfo) 
{
    int t, m = 0;
    double xx;

    if (var_is_scalar(pdinfo, i)) {
	px[0] = Z[i][0];
	return 1;
    }
    
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xx = Z[i][t];
	if (na(xx)) {
	    continue;
	}
	else px[m++] = xx;
    }

    if (m == 0) {
	fprintf(stderr, "\nztox: No valid observations for variable %s\n", 
		pdinfo->varname[i]);
    } 

    return m;
}

/**
 * gretl_isdummy:
 * @t1: starting observation.
 * @t2: ending observation. 
 * @x: data series to examine.
 * 
 * Check whether variable @x has only 0 or 1 values over the
 * given sample range (or possibly missing values).
 *
 * Returns: 0 if the variable is not a 0/1 dummy, otherwise the
 * number of 1s in the series.
 */

int gretl_isdummy (int t1, int t2, const double *x)
{
    int t, m = 0, goodobs = 0;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	if (x[t] != 0.0 && x[t] != 1.0) {
	    return 0;
	}
	if (x[t] == 1.0) {
	    m++;
	}
	goodobs++;
    }

    if (m < goodobs) {
	return m;
    }

    return 0;
} 

/**
 * gretl_iszero:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation. 
 * 
 * Check whether variable @x has only zero values over the
 * given sample range (or possibly missing values).
 *
 * Returns: 1 if the variable is all zeros, otherwise 0.
 */

int gretl_iszero (int t1, int t2, const double *x)
{
    double sum = 0.0;
    int t;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    sum += x[t] * x[t];
	}
    }

    return floateq(sum, 0.0);
}

/**
 * gretl_isconst:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation. 
 * 
 * Check whether variable @x is constant over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if the variable is constant, otherwise 0.
 */

int gretl_isconst (int t1, int t2, const double *x)
{
    int t, ret = 1;

    while (na(x[t1]) && t1 <= t2) {
	t1++;
    }

    for (t=t1+1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	if (floatneq(x[t], x[t1])) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * gretl_isunits:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation. 
 * 
 * Check whether variable @x equals 1 over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_isunits (int t1, int t2, const double *x)
{
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t]) && x[t] != 1.0) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

#define FEWVALS 8

static int few_vals (int t1, int t2, const double *x)
{
    double test[FEWVALS];
    int match;
    int i, t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    match = 0;
	    for (i=0; i<n; i++) {
		if (x[t] == test[i]) {
		    match = 1;
		    break;
		}
	    }
	    if (!match) {
		if (n == FEWVALS) {
		    n++;
		    break;
		}
		test[n++] = x[t];
	    }
	}
    }

    return n;
}

 /**
 * gretl_isdiscrete:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation. 
 *
 * Checks the variable @x over the range @t1 to @t2 for discreteness.
 * This is a heuristic whose components are (a) whether the values
 * are "fairly round" (multiples of 0.25) or not, and, if test (a) is
 * passed, (b) whether the variable takes on only "few" distinct
 * values.
 * 
 * Returns: 0 if test (a) is not passed or the number of distinct values
 * is > 8; else 1 if the number of distinct values is <= 8; else 2 if 
 * the number of distinct values is <= 4.  A return of 1 is supposed
 * to indicate that it's "reasonable" to treat @x as discrete, while
 * a return of 2 indicates that it's probably ureasonable _not_ to
 * treat @x as discrete, for the purpose of drawing up a frequency
 * distribution.
 */

int gretl_isdiscrete (int t1, int t2, const double *x)
{
    int t, n = 0, d = 1;
    double r;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	n++;
	if (!ok_int(x[t])) {
	    d = 0;
	    break;
	}
	r = x[t] - floor(x[t]);
	if (r != 0.0 && r != 0.25 && r != 0.5 && r != 0.75) {
	    d = 0;
	    break;
	}	    
    }

    if (n == 0) {
	d = 0;
    }

    if (d) {
	n = few_vals(t1, t2, x);
	if (n > FEWVALS) {
	    d = 0;
	} else if (n < 5) {
	    d = 2;
	}
    }

    return d;
}

/**
 * true_const:
 * @v: index number of variable to test.
 * @Z: data array.
 * @pdinfo: dataset information. 
 * 
 * Check whether variable Z[v] equals 1 over the sample
 * range given in @pdinfo, (aside from any missing values).
 *
 * Returns: 1 if so, otherwise 0.
 */

int true_const (int v, const double **Z, const DATAINFO *pdinfo)
{
    if (v < 0 || v >= pdinfo->v) {
	return 0;
    }

    return gretl_isunits(pdinfo->t1, pdinfo->t2, Z[v]);
}

/**
 * gretl_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with %qsort.  Sorts doubles in
 * ascending order.
 * 
 * Returns: appropriate value for %qsort.
 */

int gretl_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
}

/**
 * gretl_inverse_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 * 
 * Comparison function for use with %qsort.  Sorts doubles in
 * descending order.
 * 
 * Returns: appropriate value for %qsort.
 */

int gretl_inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da < *db) - (*da > *db);
}

/**
 * gretl_compare_ints:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with %qsort.  Sorts integers in
 * ascending order.
 * 
 * Returns: appropriate value for %qsort.
 */

int gretl_compare_ints (const void *a, const void *b)
{
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
     
    return *ia - *ib;
}

/**
 * count_distinct_values:
 * @x: sorted array of doubles.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_values (const double *x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
	if (x[i] != x[i-1]) {
	    c++;
	}
    }

    return c;
}

/**
 * count_distinct_int_values:
 * @x: sorted array of ints.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_int_values (const int *x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
	if (x[i] != x[i-1]) c++;
    }

    return c;
}

/**
 * rearrange_id_array:
 * @x: sorted array of doubles.
 * @m: number of distinct values in array.
 * @n: number of elements in array.
 *
 * Rearranges the sorted array @x such that the first @m 
 * elements contain the @m distinct values in sorted order.
 *
 * Returns: 0 on success, 1 on error (in case @m is greater
 * than @n).
 */

int rearrange_id_array (double *x, int m, int n)
{
    int i, k = 1;

    if (m >= n || m == 1) {
	return 1;
    }

    for (i=1; i<n && k<m; i++) {
	if (x[i] != x[i-1]) {
	    x[k++] = x[i];
	}
    }

    return 0;
}

/**
 * printlist:
 * @list: array of integers.
 * @msg: message to print along with @list (or NULL).
 * 
 * Prints to stderr the given @list of integers along with a message.
 */

void printlist (const int *list, const char *msg)
{
    int i;

    if (msg) {
	fprintf(stderr, "%s:\n", msg);
    } else {
	fprintf(stderr, "list: ");
    }

    if (list == NULL) {
	fputs( "list is NULL", stderr);
    } else {
	fprintf(stderr, "%d : ", list[0]);
	for (i=1; i<=list[0]; i++) {
	    fprintf(stderr, "%d ", list[i]);
	}
    }

    fputc('\n', stderr);
}

/**
 * gretl_calculate_criteria:
 * @ess: error sum of squares.
 * @n: number of observations.
 * @k: number of parameters estimated.
 * @ll: pointer to recieve loglikelihood.
 * @aic: pointer to recieve Akaike criterion.
 * @bic: pointer to recieve Schwartz Bayesian criterion.
 * @hqc: pointer to recieve Hannan-Quinn criterion.
 *
 * Calculates model selection criteria based on @ess, @nobs and
 * @k, for a model estimated via least squares.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_calculate_criteria (double ess, int n, int k,
			      double *ll, double *aic, double *bic,
			      double *hqc)
{
    double lnl, c[3];
    int err = 0;

    if (na(ess) || ess <= 0.0 || k < 1 || n <= k) {
	err = 1;
    } else {
	const double ln2pi1 = 2.837877066409345; /* log(2*pi) + 1 */

	errno = 0;

	lnl = -.5 * n * log(ess);

	if (errno == EDOM || errno == ERANGE) {
	    err = 1;
	} else {
	    lnl += -.5 * n * (ln2pi1 - log((double) n));
	    c[0] = -2.0 * lnl + 2 * k;
	    c[1] = -2.0 * lnl + k * log(n);
	    c[2] = -2.0 * lnl + 2 * k * log(log(n));
	}
    }

    if (err) {
	*ll = NADBL;
	*aic = NADBL;
	*bic = NADBL;
	*hqc = NADBL;
    } else {
	*ll = lnl;
	*aic = c[0];
	*bic = c[1];
	*hqc = c[2];
    }	

    return err;
}

/**
 * ls_criteria:
 * @pmod: pointer to gretl model structure.
 *
 * Fills out the model selection criteria members of @pmod, using
 * gretl_calculate_criteria().
 *
 * Returns: 0 on success, non-zero on error.
 */

int ls_criteria (MODEL *pmod)
{
    double ll, aic, bic, hqc;
    int err;

    err = gretl_calculate_criteria(pmod->ess, pmod->nobs, pmod->ncoeff,
				   &ll, &aic, &bic, &hqc);

    pmod->lnL = ll;
    pmod->criterion[C_AIC] = aic;
    pmod->criterion[C_BIC] = bic;
    pmod->criterion[C_HQC] = hqc;

    return err;
}

/**
 * gretl_print_criteria:
 * @ess: error sum of squares.
 * @nobs: number of observations.
 * @ncoeff: number of parameters estimated.
 * @prn: printing struct.
 *
 * Prints the values of the model selection criteria AIC, 
 * BIC and HQC for the given parameters.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_print_criteria (double ess, int nobs, int ncoeff, PRN *prn)
{
    double ll, aic, bic, hqc;
    int err;

    err = gretl_calculate_criteria(ess, nobs, ncoeff, &ll, &aic, &bic, &hqc);

    if (err) {
	pputs(prn, _("Error calculating model selection criteria\n"));
    } else {
	pprintf(prn, _("Using ess = %g, %d observations, %d coefficients\n"), 
		ess, nobs, ncoeff);
	pprintf(prn, "\nAIC = %g\nBIC = %g\nHQC = %g\n\n", aic, bic, hqc);
    }

    return err;
}

static char *
real_format_obs (char *obs, int maj, int min, int pd, char sep)
{
    if (pd >= 10) {
	int pdp = pd / 10, minlen = 2;
	char fmt[16];

	while ((pdp = pdp / 10)) minlen++;
	sprintf(fmt, "%%d%c%%0%dd", sep, minlen);
	sprintf(obs, fmt, maj, min);
    } else {
	sprintf(obs, "%d%c%d", maj, sep, min);
    }

    return obs;
}

/**
 * format_obs:
 * @obs: target string (should be of length #OBSLEN).
 * @maj: major period (e.g. year).
 * @min: minor period (e.g. quarter, month).
 * @pd: data frequency.
 *
 * Prints to @obs the gretl-type date string representing 
 * the observation given by @maj, @min and @pd.
 *
 * Returns: @obs.
 */

char *format_obs (char *obs, int maj, int min, int pd)
{
    return real_format_obs(obs, maj, min, pd, ':');
}

static int get_stobs_maj_min (char *stobs, int *maj, int *min)
{
    int dotc = 0;
    char *p = stobs;
    int err = 0;

    while (*p) {
	if (*p == ':') {
	    *p = '.';
	    dotc++;
	} else if (*p == '.') {
	    dotc++;
	} else if (!isdigit((unsigned char) *p)) {
	    err = 1;
	    break;
	}
	p++;
    }

    if (!err) {
	if (dotc > 1 || *stobs == '.' || 
	    stobs[strlen(stobs) - 1] == '.') {
	    err = 1;
	}
    }

    if (!err) {
	if (dotc > 0) {
	    sscanf(stobs, "%d.%d", maj, min);
	    if (*maj <= 0 || *min <= 0) {
		err = 1;
	    }
	} else {
	    sscanf(stobs, "%d", maj);
	    if (*maj <= 0) {
		err = 1;
	    }
	}
    }

    return err;
}

static int 
catch_setobs_errors (const char *stobs, int pd, int n, int min, gretlopt opt)
{
    int err = 0;

    if (pd == 1) {
	if (min > 0) {
	    strcpy(gretl_errmsg, _("no ':' allowed in starting obs with "
				   "frequency 1"));
	    err = 1;
	} else if (opt == OPT_S || opt == OPT_C) {
	    strcpy(gretl_errmsg, _("panel data must have frequency > 1"));
	    err = 1;
	}
    } else {
	if (min == 0) {
	    strcpy(gretl_errmsg, _("starting obs must contain a ':' with "
				   "frequency > 1"));
	    err = 1;
	} else if (min > pd) {
	    sprintf(gretl_errmsg, 
		    _("starting obs '%s' is incompatible with frequency"), 
		    stobs);
	    err = 1;
	} else if (opt == OPT_X) {
	    strcpy(gretl_errmsg, _("cross-sectional data: frequency must be 1"));
	    err = 1;
	} else if (n % pd != 0) {
	    if (opt == OPT_S || opt == OPT_C) {
		sprintf(gretl_errmsg, _("Panel datasets must be balanced.\n"
					"The number of observations (%d) is not a multiple\n"
					"of the number of %s (%d)."), 
			n, ((opt == OPT_S)? _("periods") : _("units")), pd);
		err = 1;
	    }
	}
    }

    return err;
}

#define recognized_ts_frequency(f) (f == 4 || f == 12 || f == 24)

/**
 * set_obs:
 * @line: command line.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: %OPT_S for stacked time-series, %OPT_C for stacked cross-section,
 * %OPT_T for time series, %OPT_X for cross section, %OPT_P to set
 * panel structure via two variables representing unit and period
 * respectively.
 * 
 * Set the frequency and initial observation for a dataset.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int set_obs (const char *line, double **Z, DATAINFO *pdinfo, 
	     gretlopt opt)
{
    char stobs[OBSLEN];
    int structure = STRUCTURE_UNKNOWN;
    int pd, dated = 0;
    int err = 0;

    gretl_error_clear();

    if (opt & OPT_R) {
	/* restructure panel: "hidden" option */
	return switch_panel_orientation(Z, pdinfo);
    }    

    if (!strcmp(line, "setobs")) {
	/* we'll just print current obs info */
	return 0;
    }

    if (opt & OPT_P) {
	return set_panel_structure_from_line(line, Z, pdinfo);
    }

    /* now we get down to business */

    if (sscanf(line, "%*s %d %10s", &pd, stobs) != 2) {
	strcpy(gretl_errmsg, _("Failed to parse line as frequency, startobs"));
	return 1;
    }

    /* truncate stobs if not a calendar date */
    if (strchr(stobs, '/') != NULL) {
	dated = 1;
    } else {
	stobs[8] = '\0';
    }

    /* does frequency make sense? */
    if (pd < 1 || (pdinfo->n > 0 && pd > pdinfo->n && opt != OPT_T)) {
	sprintf(gretl_errmsg, 
		_("frequency (%d) does not make seem to make sense"), pd);
	return 1;
    }

    /* if an explicit structure option was passed in, respect it */
    if (opt == OPT_X) {
	structure = CROSS_SECTION;
    } else if (opt == OPT_T) {
	structure = TIME_SERIES;
    } else if (opt == OPT_S) {
	structure = STACKED_TIME_SERIES;
    } else if (opt == OPT_C) {
	structure = STACKED_CROSS_SECTION;
    } else if (opt == OPT_N) {
	structure = SPECIAL_TIME_SERIES;
    }

    if (dated) {
	if (opt == OPT_X || opt == OPT_S || opt == OPT_C) {
	    sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
	    return 1;
	}

	if (pd == 5 || pd == 6 || pd == 7 || pd == 52) {
	    /* calendar-dated data, daily or weekly */
	    double ed0 = get_epoch_day(stobs);

	    if (ed0 < 0) {
		sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
		return 1;
	    }

	    pdinfo->sd0 = ed0;
	    structure = TIME_SERIES;

	    /* replace any existing markers with date strings */
	    dataset_destroy_obs_markers(pdinfo);
	} else {
	    sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
	    return 1;
	}
    } else if (structure == TIME_SERIES && pd == 10) {
	/* decennial data */
	pdinfo->sd0 = (double) atoi(stobs);
    } else {
	int maj = 0, min = 0;

	if (get_stobs_maj_min(stobs, &maj, &min)) {
	    sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
	    return 1;
	}

	/* catch undated daily or weekly data */
	if ((pd == 5 || pd == 6 || pd == 7 || pd == 52)  
	    && min == 0 && opt != OPT_X && opt != OPT_S && opt != OPT_C) {
	    pdinfo->structure = TIME_SERIES;
	} else {
	    if (catch_setobs_errors(stobs, pd, pdinfo->n, min, opt)) {
		return 1;
	    }
	    if (pd == 1) {
		sprintf(stobs, "%d", maj);
		if (structure == STRUCTURE_UNKNOWN) {
		    if (maj > 1) {
			structure = TIME_SERIES; /* annual? */
		    } else {
			structure = CROSS_SECTION;
		    }
		}
	    } else {
		real_format_obs(stobs, maj, min, pd, '.');
		if (structure == STRUCTURE_UNKNOWN && 
		    recognized_ts_frequency(pd)) {
		    structure = TIME_SERIES;
		}
	    }
	}

	/* for non-calendar data */
	pdinfo->sd0 = dot_atof(stobs);
    }

    pdinfo->pd = pd;
    pdinfo->structure = structure;

    if (pdinfo->structure != STACKED_TIME_SERIES &&
	pdinfo->structure != STACKED_CROSS_SECTION &&
	pdinfo->paninfo != NULL) {
	/* This could be a problem in some cases, like
	   if the dataset is subsampled? */
	dataset_destroy_panel_info(pdinfo);
    }

    ntodate_full(pdinfo->stobs, 0, pdinfo); 
    ntodate_full(pdinfo->endobs, pdinfo->n - 1, pdinfo);

    /* pre-process stacked cross-sectional panels: put into canonical
       stacked time series form
    */
    if (pdinfo->structure == STACKED_CROSS_SECTION) {
	err = switch_panel_orientation(Z, pdinfo);
    }	

    if (!err && pdinfo->structure == STACKED_TIME_SERIES) {
	err = dataset_add_default_panel_indices(pdinfo);
    }

    return err;
}

/**
 * gretl_int_from_string:
 * @s: string to examine.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @err: location to receive error code.
 * 
 * If @s is a valid string representation of an integer,
 * return that integer, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * otherwise set the content of @err to a non-zero value.
 *
 * Returns: integer value.
 */

int gretl_int_from_string (const char *s, const double **Z, 
			   const DATAINFO *pdinfo, int *err)
{
    char *test;
    int n = 0;

    if (s == NULL || *s == 0) {
	*err = E_DATA;
	return 0;
    }

    errno = 0;

    n = strtol(s, &test, 10);

    if (errno == ERANGE) {
	*err = E_DATA;
	errno = 0;
	return 0;
    }

    if (*test == '\0') {
	return n;
    } else if (Z == NULL || pdinfo == NULL) {
	*err = E_DATA;
    } else if (test[1] == '\0' && is_active_index_loop_char(*test)) {
	n = loop_scalar_read(*s);
    } else {
	int v = varindex(pdinfo, s);
	double x;

	if (v >= pdinfo->v) {
	    *err = E_UNKVAR;
	} else if (var_is_series(pdinfo, v)) {
	    *err = E_TYPES;
	} else {
	    x = Z[v][0];
	    if (na(x)) {
		*err = E_MISSDATA;
	    } else {
		n = (int) x;
	    }
	} 
    }

    return n;    
}

/**
 * positive_int_from_string:
 * @s: string to examine.
 * 
 * If @s is a valid string representation of a positive integer,
 * return that integer, otherwise return -1.
 *
 * Returns: integer value.
 */

int positive_int_from_string (const char *s)
{
    int ret = -1;

    if (s != NULL && *s != '\0') {
	char *test;

	errno = 0;

	ret = strtol(s, &test, 10);
	if (*test != '\0' || !strcmp(s, test) || errno == ERANGE) {
	    ret = -1;
	} 
    }

    return ret;
}

/**
 * varnum_from_string:
 * @str: string representation of an integer ID number.
 * @pdinfo: dataset information.
 * 
 * Returns: integer ID number, or -1 on failure.
 */

int varnum_from_string (const char *str, DATAINFO *pdinfo)
{
    int v = positive_int_from_string(str);

    if (v <= 0 || v >= pdinfo->v) {
	v = -1;
    } 
    
    return v;
}

/**
 * rename_var_by_id:
 * @idstr: string representation of the ID number of the
 * variable to be renamed.
 * @vname: new name to give the variable.
 * @pdinfo: dataset information.
 * 
 * Returns: 0 on sucess, %E_DATA on error.
 */

int rename_var_by_id (const char *idstr, const char *vname, 
		      DATAINFO *pdinfo)
{
    int v = varnum_from_string(idstr, pdinfo);

    if (v < 0) {
	return E_DATA;
    }

    if (var_is_const(pdinfo, v)) {
	return overwrite_err(pdinfo, v);
    }

    /* should be pre-checked for validity of varname and
       non-duplication (see interact.c under RENAME)
    */

    strcpy(pdinfo->varname[v], vname);

    return 0;
}

/**
 * copyvec:
 * @src: array of doubles.
 * @n: number of elements to copy.
 * 
 * Returns: an allocated copy of the first @n elements of
 * array @src, or %NULL on failure.
 */

double *copyvec (const double *src, int n)
{
    double *targ;
    int i;

    if (n == 0 || src == NULL) {
	return NULL;
    }

    targ = malloc(n * sizeof *targ);
    if (targ == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	targ[i] = src[i];
    }

    return targ;
}

/**
 * doubles_array_free:
 * @X: 2-dimensional array of doubles.
 * @m: number of sub-arrays.
 *
 * Frees a 2-dimensional array of doubles, first freeing
 * each sub-array.
 */

void doubles_array_free (double **X, int m)
{
    int i;

    if (X != NULL) {
	for (i=0; i<m; i++) {
	    free(X[i]);
	}
	free(X);
    }
}

/**
 * doubles_array_new:
 * @m: number of sub-arrays.
 * @n: length of each sub-array.
 *
 * Allocates a 2-dimensional array of doubles, that is,
 * @m arrays each containing @n elements.
 * 
 * Returns: the allocated array, or %NULL on failure.
 */

double **doubles_array_new (int m, int n)
{
    double **X;
    int i;

    X = malloc(m * sizeof *X);

    if (X == NULL) {
	return X;
    }

    for (i=0; i<m; i++) {
	X[i] = NULL;
    }

    for (i=0; i<m; i++) {
	X[i] = malloc(n * sizeof **X);
	if (X[i] == NULL) {
	    doubles_array_free(X, m);
	    X = NULL;
	    break;
	}
    }

    return X;
}

/**
 * data_array_from_model:
 * @pmod: reference model.
 * @Z: main data array.
 * @missv: should equal 1 if there are missing values to be
 * skipped, else 0.
 *
 * Constructs a dataset containing all the variables referenced in
 * @pmod.  The arrays start at the correct sample offset for @pmod,
 * and are contiguous.  If @missvals equals 0, this is done by creating
 * a set of pointers into the main dataset, but if there are missing
 * values to be handled, the sub-arrays are newly allocated and purged
 * of NAs.
 *
 * Returns: two-dimensional array, or %NULL on failure.
 */

double **data_array_from_model (const MODEL *pmod, double **Z, int missv)
{
    double **X;
    int nv = pmod->list[0];
    int offset = pmod->t1;
    int v, i;

    if (missv) {
	X = doubles_array_new(nv, pmod->nobs);
    } else {
	X = malloc(nv * sizeof *X);
    }

    if (X == NULL) return NULL;

    if (missv) {
	int t, s;

	for (t=0; t<pmod->nobs; t++) {
	    X[0][t] = 1.0;
	}

	for (i=1; i<nv; i++) {
	    v = (i == 1)? pmod->list[1] : pmod->list[i + 1];
	    s = 0;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!na(pmod->uhat[t])) {
		    X[i][s++] = Z[v][t];
		}
	    }
	}
    } else {
	/* constant in slot 0 */
	X[0] = Z[0] + offset;

	/* dependent var in slot 1 */
	X[1] = Z[pmod->list[1]] + offset;

	/* independent vars in slots 2, 3, ... */
	for (i=2; i<nv; i++) {
	    v = pmod->list[i + 1];
	    X[i] = Z[v] + offset;
	}
    }

    return X;
}

/**
 * re_estimate:
 * @model_spec: estimation command line.
 * @tmpmod: pointer to model to recieve results.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 *
 * Returns: 0 on success, non-zero on failure.
 */

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo) 
{
    CMD cmd;
    double rho = 0.0;
    int err = 0;

    if (gretl_cmd_init(&cmd)) {
	return 1;
    }

    err = parse_command_line(model_spec, &cmd, pZ, pdinfo);
    if (err) {
	gretl_cmd_free(&cmd);
	return err;
    }

    gretl_model_init(tmpmod);

    switch (cmd.ci) {
    case AR:
	*tmpmod = ar_func(cmd.list, pZ, pdinfo, OPT_NONE, NULL);
	break;
    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, pZ, pdinfo, cmd.ci, 
			   &err, cmd.opt, NULL);
	if (!err) {
	    *tmpmod = ar1_lsq(cmd.list, pZ, pdinfo, cmd.ci, 0, rho);
	}
	break;
    case HSK:
	*tmpmod = hsk_func(cmd.list, pZ, pdinfo);
	break;
    case LOGIT:
    case PROBIT:
	*tmpmod = logit_probit(cmd.list, pZ, pdinfo, cmd.ci, cmd.opt, 
			       NULL);
	break;
    case TOBIT:
	*tmpmod = tobit_model(cmd.list, pZ, pdinfo, NULL);
	break;
    case POISSON:
	*tmpmod = poisson_model(cmd.list, pZ, pdinfo, NULL);
	break;
    case OLS:
    case WLS:
    case HCCM:
	*tmpmod = lsq(cmd.list, pZ, pdinfo, cmd.ci, cmd.opt);
	break;
    case TSLS:
	break;
    default:
	break;
    }

    if (tmpmod->errcode) {
	err = 1;
	clear_model(tmpmod);
    }

    gretl_cmd_free(&cmd);

    return err;
}

double get_xvalue (int i, const double **Z, const DATAINFO *pdinfo)
{
    if (var_is_series(pdinfo, i)) {
	return Z[i][pdinfo->t1];
    } else {
	return Z[i][0];
    }	
}

#ifndef WIN32
# ifdef USE_GSPAWN

static int font_not_found (const char *s)
{
    /* "Could not find/open font when opening font X, using default" 
       or "gnuplot_x11: Some character sets not available" */

    if (strstr(s, "using default") || strstr(s, "character sets not available")) {
	return 1;
    } else {
	return 0;
    }
}

int gretl_spawn (char *cmdline)
{
    GError *error = NULL;
    gchar *errout = NULL, *sout = NULL;
    int ok, status;
    int ret = 0;

    gretl_error_clear();

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_command_line_sync (cmdline,
				    &sout,   /* standard output */
				    &errout, /* standard error */
				    &status, /* exit status */
				    &error);

    if (!ok) {
	strcpy(gretl_errmsg, error->message);
	fprintf(stderr, "gretl_spawn: '%s'\n", error->message);
	g_error_free(error);
	ret = 1;
    } else if (errout && *errout) {
	fprintf(stderr, "stderr: '%s'\n", errout);
	if (!font_not_found(errout)) {
	    strcpy(gretl_errmsg, errout);
	    fprintf(stderr, "gretl_errmsg: '%s'\n", gretl_errmsg);
	    ret = 1;
	}
    } else if (status != 0) {
	if (sout != NULL) {
	    sprintf(gretl_errmsg, "%s\n%s", 
		    _("Command failed"),
		    sout);
	    fprintf(stderr, "status=%d: '%s'\n", status, sout);
	} else {
	    strcpy(gretl_errmsg, _("Command failed"));
	    fprintf(stderr, "status=%d\n", status);
	}
	ret = 1;
    }

    if (errout != NULL) g_free(errout);
    if (sout != NULL) g_free(sout);

    if (ret) {
	fprintf(stderr, "Failed command: '%s'\n", cmdline);
    } 

    return ret;
}

# else /* now non-glib2 version */

int gretl_spawn (char *cmdline)
{
    int err;

    errno = 0;

    signal(SIGCHLD, SIG_DFL);

    err = system(cmdline);
    if (err) {
	fprintf(stderr, "Failed command: '%s'\n", cmdline);
	perror(NULL);
    }

    return err;
}

# endif
#endif /* !WIN32 */

/* file copying */

int gretl_copy_file (const char *src, const char *dest) 
{
    FILE *srcfd, *destfd;
    char buf[8192];
    size_t n;

    if (!strcmp(src, dest)) {
	return 1;
    }
   
    if ((srcfd = gretl_fopen(src, "rb")) == NULL) {
	sprintf(gretl_errmsg, _("Couldn't open %s"), src);
	return 1; 
    }

    if ((destfd = gretl_fopen(dest, "wb")) == NULL) {
	sprintf(gretl_errmsg, _("Couldn't write to %s"), dest);
	fclose(srcfd);
	return 1;
    }

    while ((n = fread(buf, 1, sizeof buf, srcfd)) > 0) {
	fwrite(buf, 1, n, destfd);
    }

    fclose(srcfd);
    fclose(destfd);

    return 0;
}    

/* library init and cleanup functions */

void libgretl_init (void)
{
    libset_init();
    gretl_rand_init();
    set_gretl_tex_preamble(); 
}

void libgretl_session_cleanup (void)
{
    gretl_saved_objects_cleanup();
    gretl_transforms_cleanup();
    gretl_lists_cleanup();
    destroy_user_matrices();
    gretl_plotx(NULL);
}

void libgretl_cleanup (void)
{
    const char *p;

    libgretl_session_cleanup();

    gretl_rand_free();
    gretl_functions_cleanup();
    libset_cleanup();
    gretl_command_hash_cleanup();
    saved_strings_cleanup();
    lapack_mem_free();

    p = strstr(gretl_plotfile(), "gpttmp");
    if (p != NULL) {
	int pnum;

	if (!sscanf(p, "gpttmp%d.plt", &pnum)) {
	    remove(gretl_plotfile());
	}
    }
}

/* record and retrieve hypothesis test results */

enum {
    SET_TEST_STAT,
    GET_TEST_STAT,
    GET_TEST_PVAL
};

static double
record_or_get_test_result (double teststat, double pval, char *instr,
			   int code)
{
    static char savestr[MAXLABEL] = {0};
    static double val = NADBL;
    static double pv = NADBL;

    double ret = NADBL;

    if (code == SET_TEST_STAT) {
	val = teststat;
	pv = pval;
	*savestr = '\0';
	if (instr != NULL) {
	    strncat(savestr, instr, MAXLABEL - 1);
	} 
    } else if (code == GET_TEST_STAT || code == GET_TEST_PVAL) {
	if (instr != NULL) {
	    if (code == GET_TEST_STAT) {
		sprintf(instr, _("%s test"), savestr);
	    } else {
		/* GET_TEST_PVAL */
		sprintf(instr, _("p-value for %s test"), savestr);
	    }
	}
	ret = (code == GET_TEST_STAT)? val : pv;
    } 
	
    return ret;
}

void record_test_result (double teststat, double pval, char *blurb)
{
    record_or_get_test_result(teststat, pval, blurb, SET_TEST_STAT);
}

double get_last_test_statistic (char *blurb)
{
    return record_or_get_test_result(0, 0, blurb, GET_TEST_STAT);
}

double get_last_pvalue (char *blurb)
{
    return record_or_get_test_result(0, 0, blurb, GET_TEST_PVAL);
}
