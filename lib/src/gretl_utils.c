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

/* utils.c for gretl  */

#include "libgretl.h"
#include "gretl_private.h"
#include "system.h"

#include <errno.h>

#ifndef WIN32
# include <glib.h>
# include <signal.h>
# if GLIB_CHECK_VERSION(2,0,0)
#  define GRETL_GLIB 2
# endif /* GLIB_CHECK_VERSION */
#endif /* ! WIN32 */

static int allocate_fit_resid_arrays (FITRESID *fr, int n, int errs);

/*
  returns the simple correlation coefficient between the the arrays zx
  and zy, for the n observations 0 to n-1.  returns NADBL if square
  root argument is invalid or no of observations is zero
*/

double gretl_corr (int n, const double *x, const double *y)
{
    int i, nn;
    double sx, sy, sxx, syy, sxy, den, xbar, ybar;
    double cval = 0.0;

    if (n == 0) {
	return NADBL;
    }

    if (gretl_isconst(0, n-1, x) || gretl_isconst(0, n-1, y)) {
	return NADBL;
    }

    nn = n;
    sx = sy = 0.0;
    for (i=0; i<n; ++i) {
        if (na(x[i]) || na(y[i])) {
            nn--;
            continue;
        }
        sx += x[i];
        sy += y[i];
    }

    if (nn == 0) {
	return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxx = syy = sxy = 0.0;

    for (i=0; i<n; ++i) {
        if (na(x[i]) || na(y[i])) {
	    continue;
	}
        sx = x[i] - xbar;
        sy = y[i] - ybar;
	sxx += sx * sx;
	syy += sy * sy;
	sxy += sx * sy;
    }

    if (sxy != 0.0) {
        den = sxx * syy;
        if (den > 0.0) {
	    cval = sxy / sqrt(den);
        } else {
	    cval = NADBL;
	}
    }

    return cval;
}

/* .......................................................  */

double gretl_covar (int n, const double *x, const double *y)
{
    int i, nn;
    double sx, sy, sxy, xi, yi, xbar, ybar;

    if (n == 0) {
	return NADBL;
    }

    nn = n;
    sx = sy = 0.0;

    for (i=0; i<n; ++i) {
        xi = x[i];
        yi = y[i];
        if (na(xi) || na(yi)) {
            nn--;
            continue;
        }
        sx += xi;
        sy += yi;
    }

    if (nn == 0) {
	return NADBL;
    }

    xbar = sx / nn;
    ybar = sy / nn;
    sxy = 0.0;

    for (i=0; i<n; i++) {
        xi = x[i];
        yi = y[i];
        if (na(xi) || na(yi)) {
	    continue;
	}
        sx = xi - xbar;
        sy = yi - ybar;
        sxy = sxy + (sx * sy);
    }

    return sxy / (nn - 1);
}

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

#if 0
    fprintf(stderr, "ztox: working on %s\n", pdinfo->varname[i]);
#endif

    if (!pdinfo->vector[i]) {
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

#if 0
    else if (m < pdinfo->t2 - pdinfo->t1 + 1) {
	fprintf(stderr, "\nztox: Dropped missing obs for var %s\n",
		pdinfo->varname[i]);
    }
#endif

    return m;
}

/**
 * isdummy:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation. 
 * 
 * Check whether variable @x has only 0 or 1 values over the
 * given sample range (or possibly missing values).
 *
 * Returns: 0 if the variable is not a 0/1 dummy, otherwise the
 * number of 1s in the series.
 */

int isdummy (const double *x, int t1, int t2)
{
    int t, m = 0;

    for (t=t1; t<=t2; t++) {
#if 1
	if (floatneq(x[t], 0.0) && floatneq(x[t], 1.0) && !na(x[t])) {
	    return 0;
	}
#else
	if (floatneq(x[t], 0.0) && floatneq(x[t], 1.0)) {
	    return 0;
	}
#endif
	if (floateq(x[t], 1.0)) m++;
    }

    if (m < t2 - t1 + 1) {
	return m;
    }

    return 0;
} 

/* ........................................................  */

int gretl_iszero (int t1, int t2, const double *x)
/*  checks whether all obs are zero for variable x from t1 to t2 */
{
    double sum = 0.0;
    int t;

    for (t=t1; t<=t2; t++) {
        sum += x[t] * x[t];
    }

    return floateq(sum, 0.0);
}

/**
 * list_exclude:
 * @n: position of element to be removed (zero-based). 
 * @list: array of integers.
 * 
 * Removes the element at position @n within @list.
 */

void list_exclude (int n, int *list)
{
    int i;

    for (i=n; i<list[0]; i++) {
	list[i] = list[i+1];
    }

    list[0] = list[0] - 1;
}

/* ........................................................  */

int gretl_isconst (int t1, int t2, const double *x)
{
    
    int t;

    while (na(x[t1]) && t1 <= t2) {
	t1++;
    }

    for (t=t1+1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	if (floatneq(x[t], x[t1])) {
	    return 0;
	}
    }

    return 1;
}

/* returns mean of array x from obs t1 through t2 */

double gretl_mean (int t1, int t2, const double *x)
{
    int n;
    register int t;
    double xbar, sum = 0.0;

    n = t2 - t1 + 1;
    if (n <= 0) {
	return NADBL;
    }

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += x[t];
	} else {
	    n--;
	}
    }

    xbar = sum / n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    sum += (x[t] - xbar); 
	}
    }

    return xbar + sum / n;
}

/*  returns min and max of array zx for sample t1 through t2  */

void gretl_minmax (int t1, int t2, const double *x, 
		   double *min, double *max)
{
    int t;

    while (na(x[t1]) && t1 <= t2) {
	t1++;
    }

    if (t1 >= t2) {
        *min = *max = NADBL;
        return;
    }

    *min = x[t1];
    *max = x[t1];

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) {
	    *max = x[t] > *max ? x[t] : *max;
	    *min = x[t] < *min ? x[t] : *min;
	}
    }
}

/* check if a var list contains a constant (variable with ID
   number 0) */

int gretl_hasconst (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
        if (list[i] == 0) {
	    return 1;
	}
    }

    return 0;
}

/* ...................................................... */

int gretl_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
}

/*  returns standard deviation of array x from t1 through t2
    return NADBL if square root argument is invalid
    or there are no observations
*/

double gretl_stddev (int t1, int t2, const double *x)
{
    double xx;

    xx = gretl_variance(t1, t2, x);

    return (na(xx))? xx : sqrt(xx);
}

double gretl_variance (int t1, int t2, const double *x)
{
    int i, n;
    double sumsq, xx, xbar;

    n = t2 - t1 + 1;

    if (n == 0) {
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (i=t1; i<=t2; i++) {
	if (!na(x[i])) {
	    xx = x[i] - xbar;
	    sumsq += xx*xx;
	} else {
	    n--;
	}
    }

    sumsq = (n > 1)? sumsq / (n - 1) : 0.0;

    return (sumsq >= 0)? sumsq : NADBL;
}

double gretl_sst (int t1, int t2, const double *x)
{
    int i;
    double sumsq, xx, xbar;

    if (t2 - t1 + 1 == 0) {
	return NADBL;
    }

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) {
	return NADBL;
    }

    sumsq = 0.0;

    for (i=t1; i<=t2; i++) {
	if (!na(x[i])) {
	    xx = x[i] - xbar;
	    sumsq += xx * xx;
	} 
    }

    return sumsq;
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
	for (i=0; i<=list[0]; i++) {
	    fprintf(stderr, "%d ", list[i]);
	}
    }

    fputc('\n', stderr);
}

/**
 * print_list_to_buffer:
 * @list: array of integers.
 * @buf: buffer to which list should be printed.
 * @len: length of the buffer
 * 
 * Prints to @buf the given @list of integers.
 *
 * Returns: 0 on success, 1 if the buffer is too small to contain
 * the printed list.
 */

int print_list_to_buffer (const int *list, char *buf, size_t len)
{
    int i;
    char numstr[16];
    size_t test = 0;

    *buf = '\0';

    for (i=1; i<=list[0]; i++) {
	sprintf(numstr, "%d ", list[i]);
	test += strlen(numstr);
	if (test >= len) {
	    *buf = '\0';
	    return 1;
	}
	strcat(buf, numstr);
    }

    return 0;
}

/* ....................................................... */

int calculate_criteria (double *x, double ess, int nobs, int ncoeff)
{
    if (na(ess) || ess <= 0.0 || ncoeff < 1 || nobs <= ncoeff) {
	x[C_AIC] = NADBL;
	x[C_BIC] = NADBL;

	return 1;
    } else {
	const double ln2pi1 = 2.837877066409345;
	double ll;

	errno = 0;
	ll = -.5 * nobs * log(ess);

	if (errno == EDOM || errno == ERANGE) {
	    x[C_AIC] = NADBL;
	    x[C_BIC] = NADBL;
	} else {
	    ll += -.5 * nobs * (ln2pi1 - log((double) nobs));
	    x[C_AIC] = -2.0 * ll + 2 * ncoeff;
	    x[C_BIC] = -2.0 * ll + ncoeff * log(nobs);
	}

	return 0;
    }
}

/* Compute model selection criteria */

int ls_aic_bic (MODEL *pmod)
{
    return calculate_criteria(pmod->criterion, 
			      pmod->ess, pmod->nobs,
			      pmod->ncoeff);
}

int gretl_criteria (double ess, int nobs, int ncoeff, PRN *prn)
{
    double x[2];
    int err;

    err = calculate_criteria(x, ess, nobs, ncoeff);

    if (err) {
	pputs(prn, _("Error calculating model selection criteria\n"));
    } else {
	pprintf(prn, _("Using ess = %g, %d observations, %d coefficients\n"), 
		ess, nobs, ncoeff);
	pprintf(prn, "\nAIC = %g\nBIC = %g\n\n", x[C_AIC], x[C_BIC]);
    }

    return err;
}

char *real_format_obs (char *obs, int maj, int min, int pd, char sep)
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

#define recognized_ts_frequency(f) (f == 4 || f == 12 || f == 24)

/**
 * set_obs:
 * @line: command line.
 * @pdinfo: data information struct.
 * @opt: OPT_S for stacked time-series, OPT_C for stacked cross-section,
 * OPT_T for time series, OPT_X for cross section.
 * 
 * Set the frequency and initial observation for a dataset.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int set_obs (const char *line, DATAINFO *pdinfo, gretlopt opt)
{
    char stobs[OBSLEN];
    int structure = STRUCTURE_UNKNOWN;
    int pd, dated = 0;

    *gretl_errmsg = '\0';

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
	    destroy_dataset_markers(pdinfo);
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
	    int balanced = 1;
	    int err = 0;

	    /* various pathologies */

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
		} else if (pdinfo->n % pd != 0) {
		    balanced = 0;
		    if (opt == OPT_S || opt == OPT_C) {
			sprintf(gretl_errmsg, _("Panel datasets must be balanced.\n"
						"The number of observations (%d) is not a multiple\n"
						"of the number of %s (%d)."), 
				pdinfo->n, ((opt == OPT_S)? _("periods") : _("units")), pd);
			err = 1;
		    }
		}
	    }

	    if (err) {
		return 1;
	    }

	    /* OK? */
		
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
		if (structure == STRUCTURE_UNKNOWN) {
		    if (maj > 1500 && recognized_ts_frequency(pd)) {
			structure = TIME_SERIES;
		    } else if (balanced) {
			structure = STACKED_TIME_SERIES; /* panel? */
		    } else {
			structure = TIME_SERIES; /* ?? */
		    }
		}
	    }
	}

	/* for non-calendar data */
	pdinfo->sd0 = dot_atof(stobs);
    }

    pdinfo->pd = pd;
    pdinfo->structure = structure;

    ntodate_full(pdinfo->stobs, 0, pdinfo); 
    ntodate_full(pdinfo->endobs, pdinfo->n - 1, pdinfo);

    return 0;
}

/* ........................................................... */

int is_model_cmd (const char *s)
{
    if (s == NULL || *s == '\0') return 0;

    if (!strncmp(s, "ols", 3)  ||
	!strncmp(s, "corc", 4) ||
	!strncmp(s, "hilu", 4) ||
	!strncmp(s, "wls", 3)  ||
	!strncmp(s, "pwe", 3)  ||
	!strncmp(s, "pooled", 6)  ||
	!strncmp(s, "hccm", 4) ||
	!strncmp(s, "hsk", 3)  ||
	!strncmp(s, "add", 3)  ||
	!strncmp(s, "lad", 3)  ||
	!strncmp(s, "omit", 4) ||
	!strncmp(s, "tsls", 4) ||
	!strncmp(s, "logit", 5)  ||
	!strncmp(s, "probit", 6) ||
	!strncmp(s, "tobit", 5) ||
	!strncmp(s, "poisson", 7) ||
	!strncmp(s, "garch", 5) ||
	!strncmp(s, "logistic", 8) ||
	!strncmp(s, "end nls", 7) ||
	!strncmp(s, "arma", 4) ||
	!strncmp(s, "ar ", 3) ||
	!strcmp(s, "ar")) {
	return 1;
    }

    return 0;
}

int is_model_ref_cmd (int ci)
{
    if (ci == ADD ||
	ci == OMIT ||
	ci == ARCH ||
	ci == CHOW ||
	ci == CUSUM ||
	ci == LMTEST ||
	ci == LEVERAGE ||
	ci == VIF ||
	ci == RESTRICT ||
	ci == FCAST ||
	ci == FCASTERR ||
	ci == FIT) {
	return 1;
    }

    return 0;
}

int is_quiet_model_test (int ci, gretlopt opt)
{
    if ((opt & OPT_Q) && (ci == OMIT || ci == ADD ||
			  ci == OMITFROM || ci == ADDTO)) {
	return 1;
    }

    return 0;
}

/* .......................................................... */

int list_dups (const int *list, int ci)
{
    int i, j, start = 2;

    if (ci == ARCH) start = 3;

    if (ci == TSLS || ci == AR || ci == ARMA || 
	ci == SCATTERS || ci == MPOLS || ci == GARCH) {
	for (i=2; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
    }
    
    for (i=start; i<list[0]; i++) {
	for (j=start+1; j<=list[0]; j++) {
	    if (i != j && list[i] == list[j]) {
		return list[i];
	    }
	}
    }

    return 0;
}

/* ........................................................... */

int *copylist (const int *src)
{
    int *targ;
    int i;

    if (src == NULL) {
	return NULL;
    }

    targ = malloc((src[0] + 1) * sizeof *targ);
    if (targ == NULL) {
	return NULL;
    }

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }

    return targ;
}

/* ......................................................  */

static int reallocate_markers (DATAINFO *pdinfo, int n)
{
    char **S;
    int t;

    S = realloc(pdinfo->S, n * sizeof *S);
    if (S == NULL) {
	return 1;
    }

    for (t=pdinfo->n; t<n; t++) {
	S[t] = malloc(OBSLEN);
	if (S[t] == NULL) {
	    int j;

	    for (j=pdinfo->n; j<t; j++) {
		free(S[j]);
	    }
	    free(S);
	    return 1;
	}
	S[t][0] = '\0';	    
    }

    pdinfo->S = S;

    return 0;
}

int grow_nobs (int newobs, double ***pZ, DATAINFO *pdinfo)
{
    double *x;
    int i, t, bign;

    if (newobs <= 0) return 0;

    bign = pdinfo->n + newobs;

    for (i=0; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
	    x = realloc((*pZ)[i], bign * sizeof *x);
	    if (x == NULL) {
		return E_ALLOC;
	    }
	    (*pZ)[i] = x;
	    for (t=pdinfo->n; t<bign; t++) {
		(*pZ)[i][t] = (i == 0)? 1.0 : NADBL;
	    }	    
	}
    }
    
    if (pdinfo->markers && pdinfo->S != NULL) {
	if (reallocate_markers(pdinfo, bign)) {
	    return E_ALLOC;
	}
    }
    
    if (pdinfo->t2 == pdinfo->n - 1) {
	pdinfo->t2 = bign - 1;
    }

    pdinfo->n = bign;

    /* does daily data need special handling? */
    ntodate(pdinfo->endobs, bign - 1, pdinfo);

    return 0;
}

/* ......................................................  */

static int real_dataset_add_vars (int newvars, double *x,
				  double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char *vector;
    VARINFO **varinfo;
    int i, n = pdinfo->n, v = pdinfo->v;    

    newZ = realloc(*pZ, (v + newvars) * sizeof *newZ);  
    if (newZ == NULL) return E_ALLOC;

    if (newvars == 1 && x != NULL) {
	/* new var is pre-allocated */
	newZ[v] = x;
    } else {
	for (i=0; i<newvars; i++) {
	    newZ[v+i] = malloc(n * sizeof **newZ);
	    if (newZ[v+i] == NULL) return E_ALLOC;
	}
    }

    *pZ = newZ;

    varname = realloc(pdinfo->varname, (v + newvars) * sizeof *varname);
    if (varname == NULL) {
	return E_ALLOC;
    }

    pdinfo->varname = varname;

    for (i=0; i<newvars; i++) {
	pdinfo->varname[v+i] = malloc(VNAMELEN);
	if (pdinfo->varname[v+i] == NULL) {
	    return E_ALLOC;
	}
	pdinfo->varname[v+i][0] = '\0';
    }

    if (pdinfo->varinfo != NULL) {
	varinfo = realloc(pdinfo->varinfo, (v + newvars) * sizeof *varinfo);
	if (varinfo == NULL) {
	    return E_ALLOC;
	} else {
	    pdinfo->varinfo = varinfo;
	}
	for (i=0; i<newvars; i++) {
	    pdinfo->varinfo[v+i] = malloc(sizeof **varinfo);
	    if (pdinfo->varinfo[v+i] == NULL) {
		return E_ALLOC;
	    }
	    gretl_varinfo_init(pdinfo->varinfo[v+i]);
	}
    }

    vector = realloc(pdinfo->vector, (v + newvars));
    if (vector == NULL) {
	return E_ALLOC;
    }

    pdinfo->vector = vector;

    for (i=0; i<newvars; i++) {
	pdinfo->vector[v+i] = 1;
    }

    pdinfo->v += newvars;

    return 0;
}

/* ......................................................  */

int dataset_add_vars (int newvars, double ***pZ, DATAINFO *pdinfo)
{
    return real_dataset_add_vars(newvars, NULL, pZ, pdinfo);
}

int dataset_add_allocated_var (double *x, double ***pZ, DATAINFO *pdinfo)
{
    return real_dataset_add_vars(1, x, pZ, pdinfo);
}

/* ......................................................  */

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char *vector;
    VARINFO **varinfo;
    int n = pdinfo->n, v = pdinfo->v;    

    newZ = realloc(*pZ, (v + 1) * sizeof *newZ);  
    if (newZ == NULL) return E_ALLOC;

    newZ[v] = malloc(n * sizeof **newZ);
    if (newZ[v] == NULL) return E_ALLOC;

    *pZ = newZ;

    varname = realloc(pdinfo->varname, (v + 1) * sizeof *varname);
    if (varname == NULL) return E_ALLOC;
    pdinfo->varname = varname;

    pdinfo->varname[v] = malloc(VNAMELEN);
    if (pdinfo->varname[v] == NULL) return E_ALLOC;

    pdinfo->varname[v][0] = '\0';

    if (pdinfo->varinfo != NULL) {
	varinfo = realloc(pdinfo->varinfo, (v + 1) * sizeof *varinfo);
	if (varinfo == NULL) {
	    return E_ALLOC;
	}
	pdinfo->varinfo = varinfo;
	pdinfo->varinfo[v] = malloc(sizeof **varinfo);
	if (pdinfo->varinfo[v] == NULL) {
	    return E_ALLOC;
	}
	gretl_varinfo_init(pdinfo->varinfo[v]);
    }

    vector = realloc(pdinfo->vector, (v + 1));
    if (vector == NULL) return E_ALLOC;
    pdinfo->vector = vector;

    pdinfo->vector[v] = 0;

    pdinfo->v += 1;

    return 0;
}

/* ......................................................  */

int positive_int_from_string (const char *s)
{
    int ret;
    char *test;

    errno = 0;

    ret = strtol(s, &test, 10);

    if (*test != '\0' || !strcmp(s, test) || errno == ERANGE) {
        ret = -1;
    } 

    return ret;
}

/* ......................................................  */

int varnum_from_string (const char *str, DATAINFO *pdinfo)
{
    int varno = positive_int_from_string(str);

    if (varno <= 0 || varno >= pdinfo->v) {
	varno = -1;
    } 
    
    return varno;
}

/* ......................................................  */

static int shrink_dataset_to_size (double ***pZ,
				   DATAINFO *pdinfo,
				   int nv)
{
    char **varname;
    char *vector;
    VARINFO **varinfo;
    double **newZ;
    
    varname = realloc(pdinfo->varname, nv * sizeof *varname);
    if (varname == NULL) return E_ALLOC;
    pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, nv * sizeof *vector);
    if (vector == NULL) return E_ALLOC;
    pdinfo->vector = vector;

    varinfo = realloc(pdinfo->varinfo, nv * sizeof *varinfo);
    if (varinfo == NULL) return E_ALLOC;
    pdinfo->varinfo = varinfo;

    newZ = realloc(*pZ, nv * sizeof *newZ); 
    if (newZ == NULL) return E_ALLOC;
    *pZ = newZ;

    pdinfo->v = nv;

    return 0;
}

#undef DROPDBG

int dataset_drop_listed_vars (const int *list, double ***pZ, 
			      DATAINFO *pdinfo, int *renumber)
{
    int oldv = pdinfo->v, vmax = pdinfo->v;
    int i, v, ndel = 0; 

    if (renumber != NULL) {
	*renumber = 0;
    }

#if DROPDBG
    printlist(list, "vars to be deleted");
#endif

    /* free and set to NULL all the vars to be deleted */

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0 && v < oldv) {
	    free((*pZ)[v]);
	    (*pZ)[v] = NULL;
	    free(pdinfo->varname[v]);
	    if (pdinfo->varinfo[v] != NULL) {
		free(pdinfo->varinfo[v]);
	    }
	    ndel++;
	}
    }

    /* rearrange pointers if necessary */

    for (v=1; v<vmax; v++) {
	if ((*pZ)[v] == NULL) {
	    int gap = 1;

	    for (i=v+1; i<vmax; i++) {
		if ((*pZ)[i] == NULL) {
		    gap++;
		} else {
		    break;
		}
	    }

	    if (i < vmax) {
		vmax -= gap;
		for (i=v; i<vmax; i++) {
		    if (renumber != NULL && !hidden_var(i + gap, pdinfo)) {
			*renumber = 1;
		    }
		    pdinfo->varname[i] = pdinfo->varname[i + gap];
		    pdinfo->varinfo[i] = pdinfo->varinfo[i + gap];
		    pdinfo->vector[i] = pdinfo->vector[i + gap];
		    (*pZ)[i] = (*pZ)[i + gap];
		}		    
	    } else {
		/* deleting all subsequent vars */
		break;
	    }
	}
    }

    return shrink_dataset_to_size(pZ, pdinfo, oldv - ndel);
}

/* drop specified number of variables at the end of the dataset */

int dataset_drop_vars (int delvars, double ***pZ, DATAINFO *pdinfo)
{
    int i, v = pdinfo->v;   

    if (delvars <= 0) {
	return 0;
    }

    if (pdinfo->v <= 1) {
	return E_DATA;
    }

    for (i=v-delvars; i<v; i++) {
	if (pdinfo->varname[i] != NULL) {
	    free(pdinfo->varname[i]);
	}
	if (pdinfo->varinfo[i] != NULL) {
	    free(pdinfo->varinfo[i]);
	}
	if ((*pZ)[i] != NULL) {
	    free((*pZ)[i]);
	}
    }

    return shrink_dataset_to_size(pZ, pdinfo, v - delvars);
}

/* ........................................................... */

static void make_stack_label (char *label, char *s)
{
    char *p = strstr(s, "--");
    int len = strlen(s);

    if (p == NULL) {
	if (len > MAXLABEL - 1) {
	    strncat(label, s, MAXLABEL - 4);
	    strcat(label, "...");
	} else {
	    strcat(label, s);
	}
    } else {
	int llen = strlen(p);
	char *q = strstr(p + 2, "--");
	int sp = 1 + (q != NULL);

	len++;
	*p = '\0';

	if (len + sp > MAXLABEL - 1) {
	    strncat(label, s, MAXLABEL - 4 - (llen + sp));
	    strcat(label, "...");
	} else {
	    strcat(label, s);
	}
	strcat(label, " -");
	if (q == NULL) {
	    strcat(label, p + 1);
	} else {
	    strncat(label, p + 1, q - p - 1);
	    strcat(label, " ");
	    strcat(label, q);
	}
    }
}

/* ........................................................... */

static int get_stack_param_val (const char *s, const double **Z,
				const DATAINFO *pdinfo)
{
    int val = -1;

    if (isdigit(*s)) {
	val = atoi(s);
    } else {
	char vname[VNAMELEN];
	int i, len = strcspn(s, " -");

	if (len > VNAMELEN - 1) len = VNAMELEN - 1;
	*vname = '\0';
	strncat(vname, s, len);
	i = varindex(pdinfo, vname);
	if (i < pdinfo->v) {
	    val = (int) Z[i][0];
	}
    }

    return val;
}

static int get_optional_offset (const char *s, const double **Z,
				const DATAINFO *pdinfo, int *err)
{
    const char *p = strstr(s, "--o");
    int off = 0;

    if (p != NULL) {
	if (strncmp(p, "--offset=", 9)) {
	    *err = E_SYNTAX;
	} else {
	    off = get_stack_param_val(p + 9, Z, pdinfo);
	    if (off < 0 || off > pdinfo->n - 1) {
		*err = E_DATA;
	    }
	}
    }

    return off;
}

static int get_optional_length (const char *s, const double **Z,
				const DATAINFO *pdinfo, int *err)
{
    const char *p = strstr(s, "--l");
    int len = 0;

    if (p != NULL) {
	if (strncmp(p, "--length=", 9)) {
	    *err = E_SYNTAX;
	} else {
	    len = get_stack_param_val(p + 9, Z, pdinfo);
	    if (len < 0 || len > pdinfo->n) {
		*err = E_DATA;
	    }
	}
    }

    return len;
}

/* ........................................................... */

static int missing_tail (const double *x, int n)
{
    int i, nmiss = 0;

    for (i=n-1; i>=0; i--) {
	if (na(x[i])) nmiss++;
	else break;
    }

    return nmiss;
}

/* ........................................................... */

int dataset_stack_vars (double ***pZ, DATAINFO *pdinfo, 
			char *newvar, char *s)
{
    char vn1[VNAMELEN], vn2[VNAMELEN];
    char format[16];
    char *p, *scpy;
    int *vnum = NULL;
    double *bigx = NULL;
    int i, v1 = 0, v2 = 0, nv = 0;
    int maxok, offset;
    int oldn, bign, genv;
    int err = 0;

    scpy = gretl_strdup(s);
    if (scpy == NULL) return E_ALLOC;

    genv = varindex(pdinfo, newvar);

    s += 6;
    if (*s == ',') return E_SYNTAX;

    p = strrchr(s, ')');
    if (p == NULL) return E_SYNTAX;
    *p = '\0';

    /* do we have a range of vars? */
    sprintf(format, "%%%d[^.]..%%%ds", VNAMELEN-1, VNAMELEN-1);
    if (sscanf(s, format, vn1, vn2) == 2) {
	if (isdigit(*vn1) && isdigit(*vn2)) {
	    v1 = atoi(vn1);
	    v2 = atoi(vn2);
	} else {
	    v1 = varindex(pdinfo, vn1);
	    v2 = varindex(pdinfo, vn2);
	}
	if (v1 >= 0 && v2 > v1 && v2 < pdinfo->v) {
	    nv = v2 - v1 + 1;
	} else {
	    fputs("stack vars: range is invalid\n", stderr);
	    err = E_DATA;
	}
    } else {
	/* or do we have a comma separated list of vars? */
	char *p = s;

	while (*p) {
	    if (*p == ',') nv++;
	    p++;
	}
	nv++;

	if (nv < 2) return E_SYNTAX;

	vnum = malloc(nv * sizeof *vnum);
	if (vnum == NULL) {
	    err = E_ALLOC;
	}

	for (i=0; i<nv && !err; i++) {
	    p = strtok((i == 0)? s : NULL, ",");
	    if (isdigit(*p)) {
		v1 = atoi(p);
	    } else {
		v1 = varindex(pdinfo, p);
	    }
	    if (v1 < 0 || v1 >= pdinfo->v) {
		err = E_UNKVAR;
	    } else {
		vnum[i] = v1;
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    /* get offset specified by user? */
    offset = get_optional_offset(scpy, (const double **) *pZ, 
				 pdinfo, &err);
    if (err) {
	goto bailout;
    }

    /* get length specified by user? */
    maxok = get_optional_length(scpy, (const double **) *pZ, 
				pdinfo, &err);
    if (err) {
	goto bailout;
    }

    if (offset + maxok > pdinfo->n) {
	err = E_DATA;
	goto bailout;
    }

    if (maxok > 0) {
	bign = nv * maxok;
	if (bign < pdinfo->n) {
	    bign = pdinfo->n;
	}
    } else {
	/* calculate required series length */	
	maxok = 0;
	for (i=0; i<nv; i++) {
	    int j, ok;

	    j = (vnum == NULL)? i + v1 : vnum[i];

	    if (pdinfo->vector[j]) {
		ok = pdinfo->n - missing_tail((*pZ)[j], pdinfo->n);
	    } else {
		ok = 1;
	    }
	    if (ok > maxok) maxok = ok;
	}

	if (maxok * nv <= pdinfo->n && pdinfo->n % maxok == 0) {
	    /* suggests that at least one var has already been stacked */
	    bign = pdinfo->n;
	    maxok -= offset;
	} else {
	    /* no stacking done: need to expand series length */
	    bign = nv * (pdinfo->n - offset);
	    maxok = 0;
	}
    }

    /* allocate stacked series */
    bigx = malloc(bign * sizeof *bigx);
    if (bigx == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* extend length of all series? */
    oldn = pdinfo->n;
    if (bign > oldn) {
	err = grow_nobs(bign - oldn, pZ, pdinfo);
	if (err) {
	    free(bigx);
	    goto bailout;
	}
    }    

    /* construct stacked series */
    for (i=0; i<nv; i++) {
	int j, t, bigt, tmax;

	j = (vnum == NULL)? i + v1 : vnum[i];

	if (maxok > 0) {
	    bigt = maxok * i;
	    tmax = offset + maxok;
	} else {
	    bigt = oldn * i;
	    tmax = oldn;
	}

	for (t=offset; t<tmax; t++) {
	    if (pdinfo->vector[j]) {
		bigx[bigt] = (*pZ)[j][t];
	    } else {
		bigx[bigt] = (*pZ)[j][0];
	    }
	    if (pdinfo->S != NULL && bigt != t && 
		pdinfo->S[bigt][0] == '\0') {
		strcpy(pdinfo->S[bigt], pdinfo->S[t]);
	    }
	    bigt++;
	}

	if (i == nv - 1) {
	    for (t=bigt; t<bign; t++) {
		bigx[bigt++] = NADBL;
	    }	
	}    
    }

    /* add stacked series to dataset */
    if (genv == pdinfo->v) {
	/* add as new variable */
	err = dataset_add_allocated_var(bigx, pZ, pdinfo);
	if (err) {
	    free(bigx);
	    goto bailout;
	}
    } else {
	/* replace existing variable of same name */
	free((*pZ)[genv]);
	(*pZ)[genv] = bigx;
	gretl_varinfo_init(pdinfo->varinfo[genv]);
    }
    
    /* complete the details */
    if (!err) {
	strcpy(pdinfo->varname[genv], newvar);
	make_stack_label(VARLABEL(pdinfo, genv), scpy);
	sprintf(gretl_msg, "%s %s %s (ID %d)", 
		(genv == pdinfo->v - 1)? _("Generated") : _("Replaced"),
		_("vector"), newvar, genv);
    }

 bailout:

    free(vnum);

    return err;
}

/* ........................................................... */

int rename_var_by_id (const char *str, const char *vname, 
		      DATAINFO *pdinfo)
{
    int varno = varnum_from_string(str, pdinfo);

    if (varno < 0) return E_DATA;

    /* should be pre-checked for validity of varname and
       non-duplication (see interact.c under RENAME)
    */

    strcpy(pdinfo->varname[varno], vname);

    return 0;
}

/* ........................................................... */

int hidden_var (int i, const DATAINFO *pdinfo)
{
    if (strcmp(pdinfo->varname[i], "subdum") == 0 ||
	strcmp(pdinfo->varname[i], "annual") == 0 ||
	strcmp(pdinfo->varname[i], "qtrs") == 0 ||
	strcmp(pdinfo->varname[i], "months") == 0 ||
	strcmp(pdinfo->varname[i], "hrs") == 0 ||
	strcmp(pdinfo->varname[i], "decdate") == 0) {
	return 1;
    } else {
	return 0;
    }
}

/* ........................................................... */

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

/* ........................................................... */

int gretl_forecast (int t1, int t2, int nv, 
		    const MODEL *pmod, double ***pZ)
{
    double xx, zz, zr;
    int i, k, maxlag = 0, yno;
    int v, t, miss;
    const int *arlist = NULL;
    int ar = AR_MODEL(pmod->ci);

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == GARCH) {
	for (t=t1; t<=t2; t++) {
	    (*pZ)[nv][t] = pmod->yhat[t];
	}
	return 0;
    }

    yno = pmod->list[1];

    if (ar) {
	arlist = pmod->arinfo->arlist;
	maxlag = arlist[arlist[0]];
	if (t1 < maxlag) {
	    t1 = maxlag; 
	}
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	zz = 0.0;
	if (ar) 
	    for (k=1; k<=arlist[0]; k++) {
	    xx = (*pZ)[yno][t - arlist[k]];
	    zr = pmod->arinfo->rho[k];
	    if (na(xx)) {
		if (zr == 0) {
		    continue;
		}
		xx = (*pZ)[nv][t - arlist[k]];
		if (na(xx)) {
		    (*pZ)[nv][t] = NADBL;
		    miss = 1;
		}
	    }
	    zz = zz + xx * zr;
	} /* end if ar */
	for (v=0; !miss && v<pmod->ncoeff; v++) {
	    k = pmod->list[v+2];
	    xx = (*pZ)[k][t];
	    if (na(xx)) {
		zz = NADBL;
		miss = 1;
	    }
	    if (!miss && ar) {
		xx = (*pZ)[k][t];
		for (i=1; i<=arlist[0]; i++) {
		    xx -= pmod->arinfo->rho[i] * (*pZ)[k][t - arlist[i]];
		}
	    }
	    if (!miss) {
		zz = zz + xx * pmod->coeff[v];
	    }
	}

	if (pmod->ci == LOGISTIC) {
	    double lmax = gretl_model_get_double(pmod, "lmax");

	    zz = lmax / (1.0 + exp(-zz));
	}

	(*pZ)[nv][t] = zz;
    }

    return 0;
}

/* ........................................................... */

FITRESID *get_fit_resid (const MODEL *pmod, double ***pZ, 
			 DATAINFO *pdinfo)
{
    int depvar, t, ft;
    FITRESID *fr;

    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	depvar = pmod->list[4];
    } else {
	depvar = pmod->list[1];
    }

    fr = fit_resid_new(pmod->t2 - pmod->t1 + 1, 0);
    if (fr == NULL) {
	return NULL;
    }

    fr->sigma = pmod->sigma;

    fr->t1 = pmod->t1;
    fr->t2 = pmod->t2;
    fr->real_nobs = pmod->nobs;

    for (t=fr->t1; t<=fr->t2; t++) {
	ft = t - fr->t1;
	fr->actual[ft] = (*pZ)[depvar][t];
	fr->fitted[ft] = pmod->yhat[t];
    }

    if (isdummy(fr->actual, 0, fr->nobs) > 0) {
	fr->pmax = get_precision(fr->fitted, fr->nobs, 8);
    } else {
	fr->pmax = get_precision(fr->actual, fr->nobs, 8);
    }
    
    strcpy(fr->depvar, pdinfo->varname[depvar]);
    
    return fr;
}

/* ........................................................... */

static void fcast_adjust_t1_t2 (const MODEL *pmod,
				const double **Z,
				int *t1, int *t2)
{
    int i, t;
    int my_t1 = *t1, my_t2 = *t2;
    int imin = (pmod->ifc)? 3 : 2;
    int miss;

    for (t=*t1; t<=*t2; t++) {
	miss = 0;
	for (i=imin; i<=pmod->list[0]; i++) {
	    if (na(Z[pmod->list[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) my_t1++;
	else break;
    }

    for (t=*t2; t>0; t--) {
	miss = 0;
	for (i=imin; i<=pmod->list[0]; i++) {
	    if (na(Z[pmod->list[i]][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) my_t2--;
	else break;
    }

    *t1 = my_t1;
    *t2 = my_t2;
}

static int 
fcast_x_missing (const int *list, const double **Z, int t)
{
    int i;

    for (i=2; i<=list[0]; i++) {
	if (na(Z[list[i]][t])) return 1;
    }

    return 0;
}

/* 
   Below: the method for generating forecasts and prediction errors
   that is presented in Wooldridge's Introductory Econometrics,
   Chapter 6.
*/

FITRESID *get_fcast_with_errs (const char *str, const MODEL *pmod, 
			       double ***pZ, DATAINFO *pdinfo, 
			       PRN *prn)
{
    double **fZ = NULL;
    DATAINFO *finfo = NULL;
    MODEL fmod; 
    FITRESID *fr;
    int *flist = NULL;
    int i, k, t, n_est, nv;
    int yno = pmod->list[1];
    char t1str[OBSLEN], t2str[OBSLEN];

    fr = fit_resid_new(0, 1); 
    if (fr == NULL) {
	return NULL;
    }

    if (pmod->ci != OLS) {
	fr->err = E_OLSONLY;
	return fr;
    }

    /* bodge (reject in case of subsampled data) */
    if (gretl_model_get_int(pmod, "daily_repack")) {
	fr->err = E_DATA;
	return fr;
    }

    /* parse dates */
    if (sscanf(str, "%*s %10s %10s", t1str, t2str) != 2) {
	fr->err = E_OBS;
	return fr;
    }

    fr->t1 = dateton(t1str, pdinfo);
    fr->t2 = dateton(t2str, pdinfo);

    if (fr->t1 < 0 || fr->t2 < 0 || fr->t2 <= fr->t1) {
	fr->err = E_OBS;
	return fr;
    }

    /* move endpoints if there are missing vals for the
       independent variables */
    fcast_adjust_t1_t2(pmod, (const double **) *pZ, &fr->t1, &fr->t2);

    /* number of obs for which forecasts will be generated */
    fr->nobs = fr->t2 - fr->t1 + 1;

    if (allocate_fit_resid_arrays(fr, fr->nobs, 1)) {
	fr->err = E_ALLOC;
	return fr;
    }

    nv = pmod->list[0];
    if (!pmod->ifc) nv++;

    n_est = pmod->t2 - pmod->t1 + 1;

    finfo = create_new_dataset(&fZ, nv, n_est, 0);
    if (finfo == NULL) {
	fr->err = E_ALLOC;
	return fr;
    }

    /* insert depvar at position 1 */
    for (t=0; t<finfo->n; t++) {
	fZ[1][t] = (*pZ)[yno][t + pmod->t1];
    }

    /* create new list */
    flist = malloc((finfo->v + 1) * sizeof *flist);
    if (flist == NULL) {
	fr->err = E_ALLOC;
	goto fcast_bailout;
    }

    flist[0] = finfo->v;
    flist[1] = 1;
    flist[2] = 0;
    for (i=3; i<=flist[0]; i++) {
	flist[i] = i - 1;
    }

    gretl_model_init(&fmod);

    /* loop across the observations for which we want forecasts
       and standard errors */

#ifdef FCAST_DEBUG
    printf("get_fcast_with_errs: ft1=%d, ft2=%d, pmod->t1=%d, pmod->t2=%d\n",
	   fr->t1, fr->t2, pmod->t1, pmod->t2);
#endif

    for (k=0; k<fr->nobs; k++) {
	int tk = k + fr->t1;

	fr->actual[k] = (*pZ)[yno][tk];

	if (fcast_x_missing(pmod->list, (const double **) *pZ, tk)) {
	    fr->sderr[k] = fr->fitted[k] = NADBL;
	    continue;
	}

	/* form modified indep vars: original data minus the values
	   to be used for the forecast 
	*/
	for (i=3; i<=flist[0]; i++) {
	    int v = (pmod->ifc)? pmod->list[i] : pmod->list[i-1];
	    const double *xv = (*pZ)[v];

	    for (t=0; t<finfo->n; t++) {
		int tp = t + pmod->t1;

		if (na(xv[tp])) {
		    fZ[i-1][t] = NADBL;
		} else {
		    fZ[i-1][t] = xv[tp] - xv[tk];
		}
	    }
	}

	fmod = lsq(flist, &fZ, finfo, OLS, OPT_A, 0.0);

	if (fmod.errcode) {
	    fr->err = fmod.errcode;
	    clear_model(&fmod);
	    goto fcast_bailout;
	}

	fr->fitted[k] = fmod.coeff[0];

	/* what exactly do we want here? */
#ifdef GIVE_SDERR_OF_EXPECTED_Y
	fr->sderr[k] = fmod.sderr[0];
#else
	fr->sderr[k] = sqrt(fmod.sderr[0] * fmod.sderr[0] + 
			    fmod.sigma * fmod.sigma);
#endif
	clear_model(&fmod);
    }

    fr->tval = tcrit95(pmod->dfd);
    strcpy(fr->depvar, pdinfo->varname[yno]);
    fr->df = pmod->dfd;

 fcast_bailout:

    free_Z(fZ, finfo);
    free(flist);
    clear_datainfo(finfo, CLEAR_FULL);
    free(finfo);

    return fr;
}

/* ........................................................... */

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn,
		     int plot)
{
    FITRESID *fr;
    int err;

    fr = get_fcast_with_errs(str, pmod, pZ, pdinfo, prn);

    if (fr == NULL) {
	return E_ALLOC;
    }

    if ((err = fr->err) == 0) {
	err = text_print_fcast_with_errs(fr, pZ, pdinfo, prn, plot);
    }

    free_fit_resid(fr);
    
    return err;
}

/* ........................................................... */

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo) 
{
    CMD cmd;
    int err = 0, ignore = 0;
    double rho = 0;
    PRN prn;

    prn.fp = NULL;
    prn.buf = NULL;

    cmd.list = malloc(sizeof *cmd.list);
    cmd.param = malloc(1);

    if (cmd.list == NULL || cmd.param == NULL) {
	return 1;
    }

    getcmd(model_spec, pdinfo, &cmd, &ignore, pZ, NULL);

    gretl_model_init(tmpmod);

    switch(cmd.ci) {
    case AR:
	*tmpmod = ar_func(cmd.list, atoi(cmd.param), pZ, 
			  pdinfo, OPT_NONE, &prn);
	break;
    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, pZ, pdinfo, 1, cmd.ci, 
			   &err, &prn);
	if (!err) {
	    *tmpmod = lsq(cmd.list, pZ, pdinfo, cmd.ci, 0, rho);
	}
	break;
    case HSK:
	*tmpmod = hsk_func(cmd.list, pZ, pdinfo);
	break;
    case LOGIT:
    case PROBIT:
	*tmpmod = logit_probit(cmd.list, pZ, pdinfo, cmd.ci);
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
    case POOLED:
	*tmpmod = lsq(cmd.list, pZ, pdinfo, cmd.ci, cmd.opt, 0.0);
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

    free(cmd.list);
    free(cmd.param);

    return err;
}

/* ........................................................... */

int guess_panel_structure (double **Z, DATAINFO *pdinfo)
{
    int v, panel;

    v = varindex(pdinfo, "year");

    if (v == pdinfo->v) {
	v = varindex(pdinfo, "Year");
    }

    if (v == pdinfo->v) {
	panel = 0; /* can't guess */
    } else {
	if (floateq(Z[v][0], Z[v][1])) { /* "year" is same for first two obs */
	    pdinfo->structure = STACKED_CROSS_SECTION; 
	    panel = STACKED_CROSS_SECTION;
	} else {
	    pdinfo->structure = STACKED_TIME_SERIES; 
	    panel = STACKED_TIME_SERIES;
	}
    }

    return panel;
}

/* 
   nunits = number of cross-sectional units
   T = number pf time-periods per cross-sectional unit
*/

int get_panel_structure (const DATAINFO *pdinfo, int *nunits, int *T)
{
    int err = 0;

    if (pdinfo->structure == STACKED_TIME_SERIES) {
        *nunits = pdinfo->n / pdinfo->pd;
        *T = pdinfo->pd;
    } else if (pdinfo->structure == STACKED_CROSS_SECTION) {
	*nunits = pdinfo->pd;
	*T = pdinfo->n / pdinfo->pd;
    } else {
	err = 1;
    }

    return err;
}

/* ........................................................... */

int set_panel_structure (gretlopt opt, DATAINFO *pdinfo, PRN *prn)
{
    int nunits, T;
    int old_ts = pdinfo->structure;

    if (pdinfo->pd == 1) {
	pputs(prn, _("The current data frequency, 1, is not "
		"compatible with panel data.\nPlease see the 'setobs' "
		"command.\n"));
	return 1;
    }

    if (opt == OPT_C) {
	pdinfo->structure = STACKED_CROSS_SECTION;
    } else {
	pdinfo->structure = STACKED_TIME_SERIES;
    }

    if (get_panel_structure(pdinfo, &nunits, &T)) {
	pputs(prn, _("Failed to set panel structure\n"));
	pdinfo->structure = old_ts;
	return 1;
    } else {
	pprintf(prn, _("Panel structure set to %s\n"),
		(pdinfo->structure == STACKED_CROSS_SECTION)? 
		_("stacked cross sections") : _("stacked time series"));
	pprintf(prn, _("(%d units observed in each of %d periods)\n"),
		nunits, T);
    }

    return 0;
}

/* ........................................................... */

int balanced_panel (const DATAINFO *pdinfo)
{
    char unit[OBSLEN], period[OBSLEN];

    if ((pdinfo->t2 - pdinfo->t1 + 1) % pdinfo->pd) {
        return 0;
    }

    if (sscanf(pdinfo->endobs, "%[^:]:%s", unit, period) == 2) {
        if (atoi(period) != pdinfo->pd) return 0;
    } else {
        return 0;
    }

    return 1;
}

/* ........................................................... */

double get_xvalue (int i, const double **Z, const DATAINFO *pdinfo)
{
    if (pdinfo->vector[i]) {
	return Z[i][pdinfo->t1];
    } else {
	return Z[i][0];
    }	
}

/* ........................................................... */

static void free_mp_varnames (mp_results *mpvals)
{
    int i, n = mpvals->ncoeff + 1;

    if (mpvals->varnames != NULL) {
	for (i=0; i<n; i++) {
	    free(mpvals->varnames[i]);
	}
	free(mpvals->varnames);
    }
}

/* ........................................................... */

int allocate_mp_varnames (mp_results *mpvals)
{
    int i, j, n = mpvals->ncoeff + 1;

    mpvals->varnames = malloc(n * sizeof *mpvals->varnames);

    if (mpvals->varnames == NULL) {
	return 1;
    }

    for (i=0; i<n; i++) {
	mpvals->varnames[i] = malloc(12);
	if (mpvals->varnames[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(mpvals->varnames[j]);
	    }
	    free(mpvals->varnames);
	    return 1;
	}
	mpvals->varnames[i][0] = 0;
    }

    return 0;
}

/* ........................................................... */

void free_gretl_mp_results (mp_results *mpvals)
{
    if (mpvals != NULL) {
	free(mpvals->coeff);
	free(mpvals->sderr);
	if (mpvals->varnames != NULL) {
	    free_mp_varnames(mpvals);
	}
	if (mpvals->varlist != NULL) {
	    free(mpvals->varlist);
	}
	free(mpvals);
    }
}

/* ........................................................... */

mp_results *gretl_mp_results_new (int nc)
{
    mp_results *mpvals;
    int i;

    mpvals = malloc(sizeof *mpvals);
    if (mpvals == NULL) return NULL;

    mpvals->ncoeff = nc;

    mpvals->coeff = malloc(nc * sizeof *mpvals->coeff);
    mpvals->sderr = malloc(nc * sizeof *mpvals->sderr);
    mpvals->varnames = NULL;
    mpvals->varlist = NULL;

    if (mpvals->coeff == NULL || 
	mpvals->sderr == NULL) {
	free_gretl_mp_results(mpvals);
	return NULL;
    }

    for (i=0; i<nc; i++) mpvals->coeff[i] = NADBL;
    for (i=0; i<nc; i++) mpvals->sderr[i] = NADBL;

    mpvals->sigma = mpvals->ess = NADBL;
    mpvals->rsq = mpvals->fstt = NADBL;
    mpvals->adjrsq = NADBL;

    mpvals->t1 = mpvals->t2 = mpvals->ifc = 0;
    mpvals->dfn = mpvals->dfd = 0;

    return mpvals;
}

/* ........................................................... */

static int allocate_fit_resid_arrays (FITRESID *fr, int n, int errs)
{
    fr->actual = malloc(n * sizeof *fr->actual);
    if (fr->actual == NULL) {
	return 1;
    }

    fr->fitted = malloc(n * sizeof *fr->fitted);
    if (fr->fitted == NULL) {
	free(fr->actual);
	fr->actual = NULL;
	return 1;
    }

    if (errs) {
	fr->sderr = malloc(n * sizeof *fr->sderr);
	if (fr->sderr == NULL) {
	    free(fr->actual);
	    fr->actual = NULL;
	    free(fr->fitted);
	    fr->fitted = NULL;
	    return 1;
	}
    } else {
	fr->sderr = NULL;
    }

    return 0;
}

/* ........................................................... */

FITRESID *fit_resid_new (int n, int errs)
{
    FITRESID *fr;

    fr = malloc(sizeof *fr);
    if (fr == NULL) {
	return NULL;
    }

    fr->err = 0;
    fr->t1 = 0;
    fr->t2 = 0;
    fr->nobs = 0;
    fr->real_nobs = 0;

    if (n == 0) {
	fr->actual = NULL;
	fr->fitted = NULL;
	fr->sderr = NULL;
	return fr;
    }

    if (allocate_fit_resid_arrays(fr, n, errs)) {
	free(fr);
	return NULL;
    }

    fr->nobs = n;
    
    return fr;
}

/* ........................................................... */

void free_fit_resid (FITRESID *fr)
{
    free(fr->actual);
    free(fr->fitted);
    free(fr->sderr);
    free(fr);
}

static int copy_main_list (int **targ, const int *src)
{
    int i, n = 0;

    if (src == NULL) return 1;

    for (i=1; i<=src[0] && src[i]!=LISTSEP; i++) n++;

    if (*targ != NULL) free(*targ);

    *targ = malloc((n + 2) * sizeof *targ);
    if (*targ == NULL) {
	return 1;
    }
    
    (*targ)[0] = n;
    for (i=1; i<=n; i++) {
	(*targ)[i] = src[i];
    }

    return 0;
}

/* ........................................................... */

/**
 * get_model_confints:
 * @pmod: pointer to gretl model.
 *
 * Save the 95 percent confidence intervals for the parameter
 * estimates in @pmod.
 * 
 * Returns: pointer to #CONFINT struct containing the results.
 */

CONFINT *get_model_confints (const MODEL *pmod)
{
    int i;
    double t = tcrit95(pmod->dfd);
    CONFINT *cf;

    cf = malloc(sizeof *cf);
    if (cf == NULL) return NULL;

    cf->coeff = malloc(pmod->ncoeff * sizeof *cf->coeff);
    if (cf->coeff == NULL) {
	free(cf);
	return NULL;
    }

    cf->maxerr = malloc(pmod->ncoeff * sizeof *cf->maxerr);
    if (cf->maxerr == NULL) {
	free(cf);
	free(cf->coeff);
	return NULL;
    }

    cf->list = NULL;
    if (copy_main_list(&cf->list, pmod->list)) {
	free(cf);
	free(cf->coeff);
	free(cf->maxerr);
	return NULL;
    }

    for (i=0; i<pmod->ncoeff; i++) { 
	cf->coeff[i] = pmod->coeff[i];
	cf->maxerr[i] = (pmod->sderr[i] > 0)? t * pmod->sderr[i] : 0;
    }

    cf->df = pmod->dfd;
    cf->ifc = pmod->ifc;

    return cf;
}

/* ........................................................... */

void free_confint (CONFINT *cf)
{
    free(cf->coeff);
    free(cf->maxerr);
    free(cf->list);
    free(cf);
}

/* ........................................................... */

#if GRETL_GLIB

static int font_not_found (const char *s)
{
    /* "Could not find/open font when opening font X, using default" */

    if (strstr(s, "using default")) {
	return 1;
    } else {
	return 0;
    }
}

int gretl_spawn (const char *cmdline)
{
    GError *error = NULL;
    gchar *errout = NULL, *sout = NULL;
    int ok, status;
    int ret = 0;

    *gretl_errmsg = '\0';

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
	strcpy(gretl_errmsg, errout);
	fprintf(stderr, "stderr: '%s'\n", errout);
	if (!font_not_found(errout)) {
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

#endif

#if !defined(WIN32) && !defined(GRETL_GLIB)

int gretl_spawn (const char *cmdline)
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

#endif

/* library init and cleanup functions */

void libgretl_init (CMD *cmd)
{
    if (cmd != NULL && gretl_cmd_init(cmd)) {
	exit(EXIT_FAILURE);
    }

    gretl_rand_init();

    set_gretl_tex_preamble(); 
}

void libgretl_cleanup (CMD *cmd)
{
    const char *p;

    if (cmd != NULL) {
	gretl_cmd_free(cmd);
    }

    gretl_rand_free();
    gretl_functions_cleanup();
    gretl_equation_systems_cleanup();
    testvec(0);

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
record_or_get_test_result (double teststat, double pval, char *blurb,
			   int code)
{
    static double val = NADBL;
    static double pv = NADBL;
    static char info[128] = {0};

    double ret = NADBL;

    if (code == SET_TEST_STAT) {
	val = teststat;
	pv = pval;
	if (blurb != NULL) {
	    strcpy(info, blurb);
	} else {
	    *info = '\0';
	}
    } else if (code == GET_TEST_STAT || code == GET_TEST_PVAL) {
	if (blurb != NULL) {
	    strcpy(blurb, info);
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
