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

#include <errno.h>

#ifndef WIN32
# include <glib.h>
# include <signal.h>
# if GLIB_CHECK_VERSION(2,0,0)
#  define GRETL_GLIB 2
# endif /* GLIB_CHECK_VERSION */
#endif /* ! WIN32 */

static int allocate_fit_resid_arrays (FITRESID *fr, int n, int errs);

/* .......................................................  */

double gretl_corr (int n, const double *zx, const double *zy)
     /*
       returns the simple correlation coefficient between the the
       arrays zx and zy, for the n observations 0 to n-1.  returns
       NADBL if square root argument is invalid or no of observations
       is zero 
     */
{
    int i, nn;
    double sx, sy, sxx, syy, sxy, den, zxbar, zybar;
    double cval = 0.0;

    if (n == 0) return NADBL;

    if (gretl_isconst(0, n-1, zx) || gretl_isconst(0, n-1, zy)) 
	return NADBL;

    nn = n;
    sx = sy = 0.0;
    for (i=0; i<n; ++i) {
        if (na(zx[i]) || na(zy[i])) {
            nn--;
            continue;
        }
        sx += zx[i];
        sy += zy[i];
    }

    if (nn == 0) return NADBL;

    zxbar = sx / nn;
    zybar = sy / nn;
    sxx = syy = sxy = 0.0;

    for (i=0; i<n; ++i) {
        if (na(zx[i]) || na(zy[i])) continue;
        sx = zx[i] - zxbar;
        sy = zy[i] - zybar;
	sxx += sx * sx;
	syy += sy * sy;
	sxy += sx * sy;
    }

    if (sxy != 0.0) {
        den = sxx * syy;
        if (den > 0.0) cval = sxy / sqrt(den);
        else cval = NADBL;
    }

    return cval;
}

/* .......................................................  */

double gretl_covar (int n, const double *zx, const double *zy)
{
    register int i;
    int nn;
    double sx, sy, sxy, zxi, zyi, zxbar, zybar;

    if (n == 0) return NADBL;
    nn = n;
    sx = sy = 0.0;

    for (i=0; i<n; ++i) {
        zxi = zx[i];
        zyi = zy[i];
        if (na(zxi) || na(zyi)) {
            nn--;
            continue;
        }
        sx += zxi;
        sy += zyi;
    }

    if (nn == 0) return NADBL;
    zxbar = sx/nn;
    zybar = sy/nn;
    sxy = 0.0;

    for (i=0; i<n; i++) {
        zxi = zx[i];
        zyi = zy[i];
        if (na(zxi) || na(zyi)) continue;
        sx = zxi - zxbar;
        sy = zyi - zybar;
        sxy = sxy + (sx*sy);
    }

    return sxy/(nn - 1);
}

#define PDTON(p) (((p) == 1)? 1 : ((p) < 10)? 10 : 100)

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
    double dd;

    if (pd == 1) {
	return ((double) (ysd + nt));  
    }

    pp = nt % pd + PDTON(pd) * (sd0 - ysd) + .5;
    if (pp != pd)  {
        yy = ysd + nt/pd  + pp/pd + .5;
        yp = pp % pd;
    }  else {
        yy = ysd + nt/pd + .5;
        yp = pp;
    }

    dd = (pd < 10)? 0.1 : 0.01;

    return (yy + yp * dd);
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
    int idx;

    if (i > j) {
	int tmp = i;

	i = j;
	j = tmp;
    }

    idx = nrows * i + j - i - ((i - 1) * i / 2);
    return idx;
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

int ztox (int i, double *px, double **Z, const DATAINFO *pdinfo) 
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
	if (na(xx)) continue;
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
 * given sample range. 
 *
 * Returns: 0 if the variable is not a 0/1 dummy, otherwise the
 * number of 1s in the series.
 */

int isdummy (const double *x, int t1, int t2)
{
    int t, m = 0;

    for (t=t1; t<=t2; t++) {
	if (floatneq(x[t], 0.0) && floatneq(x[t], 1.0)) {
	    return 0;
	}
	if (floateq(x[t], 1.0)) m++;
    }

    if (m < t2 - t1 + 1) return m;

    return 0;
} 

/* ........................................................  */

int gretl_iszero (int t1, int t2, const double *x)
/*  checks whether all obs are zero for variable x from t1 to t2 */
{
    int t;
    double xx, sum = 0.0;

    for (t=t1; t<=t2; t++) {
        xx = x[t];
        sum = sum + xx * xx;
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

    for (i=n; i<list[0]; i++) list[i] = list[i+1];
    list[0] = list[0] - 1;
}

/* ........................................................  */

int gretl_isconst (int t1, int t2, const double *x)
{
    int t;
    double xx = x[t1];

    for (t=t1+1; t<=t2; t++) {
	if (floatneq(x[t], xx)) return 0;
    }
    return 1;
}

/* ............................................................  */

double gretl_mean (int t1, int t2, const double *x)
/* returns mean of array x from obs t1 through t2 */
{
    int n;
    register int t;
    double xbar, sum = 0.0;

    n = t2 - t1 + 1;
    if (n <= 0) return NADBL;

    for (t=t1; t<=t2; t++) {
	if (!(na(x[t]))) 
	    sum += x[t];
	else 
	    n--;
    }

    xbar = sum/n;
    sum = 0.0;

    for (t=t1; t<=t2; t++) 
	if (!(na(x[t]))) sum += (x[t] - xbar); 

    return xbar + sum / n;
}

/* ......................................................  */

void gretl_minmax (int t1, int t2, const double zx[], 
		   double *min, double *max)
/*  returns min and max of array zx for sample t1 through t2  */
{
    register int t;
    double xt;

    *min = zx[t1];
    *max = zx[t1];

    if (t2-t1+1 == 0) {
        *min = *max = NADBL;
        return;
    }

    for (t=t1; t<=t2; t++) {
        xt = zx[t];
	if (!(na(xt))) {
	    *max = xt > *max ? xt : *max;
	    *min = xt < *min ? xt : *min;
	}
    }
}

/* ..........................................................  */

int gretl_hasconst (const int *list)
/* check if a var list contains a constant (variable with ID
   number 0) */
{
    int i;

    for (i=2; i<=list[0]; i++) 
        if (list[i] == 0) return i;

    return 0;
}

/* ...................................................... */

int gretl_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
}

/* .............................................................  */

double gretl_stddev (int t1, int t2, const double *x)
/*  returns standard deviation of array x from t1 through t2
    return NADBL if square root argument is invalid
    or there are no observations
*/
{
    double xx;

    xx = gretl_variance(t1, t2, x);

    return (na(xx))? xx : sqrt(xx);
}

/* .............................................................  */

double gretl_variance (int t1, int t2, const double *x)
{
    int n;
    register int i;
    double sumsq, xx, xbar;

    n = t2 - t1 + 1;
    if (n == 0) return NADBL;

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) return NADBL;

    sumsq = 0.0;
    for (i=t1; i<=t2; i++) {
	if (!na(x[i])) {
	    xx = x[i] - xbar;
	    sumsq += xx*xx;
	} else n--;
    }

    sumsq = (n > 1)? sumsq/(n - 1) : 0.0;

    return (sumsq >= 0)? sumsq : NADBL;
}

/* .............................................................  */

double gretl_sst (int t1, int t2, const double *x)
{
    register int i;
    double sumsq, xx, xbar;

    if (t2 - t1 + 1 == 0) return NADBL;

    xbar = gretl_mean(t1, t2, x);
    if (na(xbar)) return NADBL;

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

    if (msg) fprintf(stderr, "%s:\n", msg);
    else fprintf(stderr, "list: ");
    for (i=0; i<=list[0]; i++) fprintf(stderr, "%d ", list[i]);
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

static void calculate_criteria (double *criterion, double ess, 
				int nobs, int ncoeff)
{
    double zz, zx, ersq, zn;

    zz = (double) (nobs - ncoeff);
    criterion[C_SGMASQ]  = ess / zz;
    ersq = ess / nobs;
    criterion[C_FPE]     = ersq * (nobs + ncoeff) / zz;
    zz = 2.0 * ncoeff / nobs;
    criterion[C_AIC]     = ersq * exp(zz);
    criterion[C_SHIBATA] = ersq * (1.0 + zz);
    criterion[C_RICE]    = ((1-zz) > 0.0)? ersq / (1 - zz) : NADBL;
    zn = (double) nobs;
    zx = log(zn);
    criterion[C_HQ]      = ersq * pow(zx, zz);
    zz = (double) ncoeff;
    zz /= zn;
    criterion[C_BIC]     = ersq * pow(zn, zz);
    zz = 1.0 - zz;
    criterion[C_GCV]     = ersq / (zz * zz);
}

/* Compute model selection criteria */

void gretl_aic_etc (MODEL *pmod)
{
    calculate_criteria(pmod->criterion, pmod->ess, pmod->nobs,
		       pmod->ncoeff);
}

void gretl_criteria (const double ess, int nobs, int ncoeff, 
		     PRN *prn)
{
    double criterion[8];

    calculate_criteria(criterion, ess, nobs, ncoeff);
    
    pprintf(prn, _("Using ess = %f, %d observations, %d coefficients\n"), 
	    ess, nobs, ncoeff);
    pputs(prn, _("\nMODEL SELECTION STATISTICS\n\n"));	
    pprintf(prn, "SGMASQ    %13g     AIC       %13g     FPE       %12g\n"
	    "HQ        %13g     SCHWARZ   %13g     SHIBATA   %12g\n"
	    "GCV       %13g",
	    criterion[C_SGMASQ], criterion[C_AIC], 
	    criterion[C_FPE], criterion[C_HQ], 
	    criterion[C_BIC], criterion[C_SHIBATA], criterion[C_GCV]);
    if (criterion[C_RICE] > 0.0) {
	pprintf(prn, "     RICE      %13g\n", criterion[C_RICE]);
    } else {
	pputs(prn, "     RICE          undefined\n");
    }
    pputc(prn, '\n');
}

/* ....................................................... */

int adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		 const double **Z, int *misst)
     /* drop first/last observations from sample if missing obs 
	encountered -- also check for missing vals within the
        remaining sample */
{
    int i, t, dwt = 0, t1min = *t1, t2max = *t2;
    double xx;

    if (pmod != NULL && gretl_model_get_int(pmod, "wt_dummy")) 
	dwt = pmod->nwt;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) continue;
	for (t=t1min; t<t2max; t++) {
	    xx = Z[list[i]][t];
	    if (dwt) xx *= Z[dwt][t];
	    if (na(xx)) t1min += 1;
	    else break;
	}
    }
    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) continue;
	for (t=t2max; t>t1min; t--) {
	    xx = Z[list[i]][t];
	    if (dwt) xx *= Z[dwt][t];
	    if (na(xx)) t2max -= 1;
	    else break;
	}
    } 
    if (misst != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) continue;
	    for (t=t1min; t<=t2max; t++) {
		xx = Z[list[i]][t];
		if (dwt) xx *= Z[dwt][t];
		if (na(xx)) {
		    *misst = t + 1;
		    return list[i];
		}
	    }
	}     
    }

    *t1 = t1min; *t2 = t2max;

    return 0;
}

/* ........................................................... */

static int real_setmiss (double missval, int varno, 
			 double **Z, DATAINFO *pdinfo) 
{
    int i, t, count = 0;
    int start = 1, end = pdinfo->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	for (t=0; t<pdinfo->n; t++) {
	    if (Z[i][t] == missval) {
		Z[i][t] = NADBL;
		count++;
	    }
	}	
    }

    return count;
}

/**
 * set_miss:
 * @LIST: list of variables to process.
 * @param: string with specification of value to treat as missing.
 * @Z: data matrix.
 * @pdinfo: pointer to data information struct.
 * @PRN: pointer to printing struct.
 * 
 * Set to "missing" each observation of each series in list that
 * has the specified value, as in @param.
 *
 */

void set_miss (LIST list, const char *param, double **Z,
	       DATAINFO *pdinfo, PRN *prn)
{
    double missval;
    int i, count;

    missval = atof(param);

    if (list[0] == 0) {
	count = real_setmiss(missval, 0, Z, pdinfo);
	if (count) { 
	    pprintf(prn, _("Set %d values to \"missing\"\n"), count);
	} else {
	    pputs(prn, _("Didn't find any matching observations\n"));
	}
	return;
    }

    for (i=1; i<=list[0]; i++) {
	if (!pdinfo->vector[list[i]]) {
	    pprintf(prn, _("The variable %s is a scalar\n"), 
		    pdinfo->varname[list[i]]);
	    continue;
	}
	count = real_setmiss(missval, list[i], Z, pdinfo);
	if (count) { 
	    pprintf(prn, _("%s: set %d observations to \"missing\"\n"), 
		    pdinfo->varname[list[i]], count);
	} else { 
	    pprintf(prn, _("%s: Didn't find any matching observations\n"),
		    pdinfo->varname[list[i]]);
	}
    }
}

/**
 * set_obs:
 * @line: command line.
 * @pdinfo: data information struct.
 * @opt: OPT_S for stacked time-series, OPT_C for stacked cross-section.
 * 
 * Impose a time-series or panel interpretation on a data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int set_obs (const char *line, DATAINFO *pdinfo, gretlopt opt)
{
    int pd, i, len, bad = 0;
    char stobs[OBSLEN], endobs[OBSLEN], endbit[7], *p;
    long ed0 = 0L;

    *gretl_errmsg = '\0';

    if (sscanf(line, "%*s %d %8s", &pd, stobs) != 2) {
	strcpy(gretl_errmsg, _("Failed to parse line as frequency, startobs"));
	return 1;
    }

    /* does frequency make sense? */
    if (pd < 1 || (pdinfo->n > 0 && pd > pdinfo->n)) {
	sprintf(gretl_errmsg, 
		_("frequency (%d) does not make seem to make sense"), pd);
	return 1;
    }

    p = stobs;
    while (*p) {
	if (*p == ':') *p = '.';
	p++;
    }

    /* special case: daily data (dated or undated) */
    if ((pd == 5 || pd == 7) && (strstr(stobs, "/") || !strcmp(stobs, "1"))) {
	if (strcmp(stobs, "1")) {
	    /* dated */
	    ed0 = get_epoch_day(stobs);
	    if (ed0 < 0) {
		sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
		return 1;
	    }
	    else pdinfo->sd0 = (double) ed0;
	} else {
	    pdinfo->sd0 = 1.0;
	}
    } else {
	int maj = 0, min = 0, dc = 0, pos;

	len = pos = strlen(stobs);
	for (i=0; i<len; i++) {
	    if (stobs[i] != '.' && !isdigit((unsigned char) stobs[i])) {
		bad = 1;
		break;
	    }
	    if (stobs[i] == '.') {
		if (dc == 0) pos = i;
		dc++;
	    }
	}
	if (bad || dc > 1) {
	    sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
	    return 1;
	}
	if (pd > 1 && pos == len) {
	    strcpy(gretl_errmsg, _("starting obs must contain a '.' with "
		   "frequency > 1"));
	    return 1;
	}
	if (pd == 1 && pos < len) {
	    strcpy(gretl_errmsg, _("no '.' allowed in starting obs with "
		   "frequency 1"));
	    return 1;
	}  

	if (pd > 1) {
	    maj = atoi(stobs);
	    strcpy(endbit, stobs + pos + 1);
	    min = atoi(endbit);
	    if (min < 0 || min > pd) {
		sprintf(gretl_errmsg, 
			_("starting obs '%s' is incompatible with frequency"), 
			stobs);
		return 1;
	    }	    
	    if (pd > 10) {
		int pdp = pd / 10, minlen = 2;
		char fmt[12];

		while ((pdp = pdp / 10)) minlen++;
		sprintf(fmt, "%%d.%%0%dd", minlen);
		sprintf(stobs, fmt, maj, min);
	    } else {
		sprintf(stobs, "%d.%d", maj, min);
	    }
	}
    }

    /* adjust data info struct */
    pdinfo->pd = pd;

    if (ed0 == 0L) {
	pdinfo->sd0 = dot_atof(stobs);
    } else {
	pdinfo->time_series = TIME_SERIES;
    }

    ntodate(pdinfo->stobs, 0, pdinfo);
    ntodate(endobs, pdinfo->n - 1, pdinfo);
    strcpy(pdinfo->endobs, endobs);

    if (opt == OPT_S) {
	pdinfo->time_series = STACKED_TIME_SERIES;
    } else if (opt == OPT_C) {
	pdinfo->time_series = STACKED_CROSS_SECTION;
    } else if (pdinfo->sd0 >= 1.0) {
        pdinfo->time_series = TIME_SERIES; /* but might be panel? */
    } else {
	pdinfo->time_series = 0;
    }

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
	    if (i != j && list[i] == list[j]) return list[i];
	}
    }

    return 0;
}

/* ........................................................... */

int *copylist (const int *src)
{
    int *targ;
    int i, n;

    if (src == NULL) return NULL;
    n = src[0];

    targ = malloc((n + 1) * sizeof *targ);
    if (targ == NULL) return NULL;

    for (i=0; i<=n; i++) {
	targ[i] = src[i];
    }

    return targ;
}

/* ......................................................  */

int grow_nobs (int newobs, double ***pZ, DATAINFO *pdinfo)
{
    double *x;
    int i, t, n = pdinfo->n, v = pdinfo->v;
    char endobs[12];

    if (newobs <= 0) return 0;

    for (i=0; i<v; i++) {
	x = realloc((*pZ)[i], (n + newobs) * sizeof *x);
	if (x == NULL) return E_ALLOC;
	else (*pZ)[i] = x;
    }
    
    if (pdinfo->markers && pdinfo->S != NULL) {
	char **S;

	if (allocate_case_markers(&S, n + newobs)) return E_ALLOC;
	else pdinfo->S = S;
    }

    pdinfo->n += newobs;
    pdinfo->t2 = pdinfo->n - 1;
    ntodate(endobs, pdinfo->t2, pdinfo);
    strcpy(pdinfo->endobs, endobs);

    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[0][t] = 1.0;
    }

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
	/* need to allocate for new vars(s) */
	for (i=0; i<newvars; i++) {
	    newZ[v+i] = malloc(n * sizeof **newZ);
	    if (newZ[v+i] == NULL) return E_ALLOC;
	}
    }

    *pZ = newZ;

    varname = realloc(pdinfo->varname, (v + newvars) * sizeof *varname);
    if (varname == NULL) return E_ALLOC;

    pdinfo->varname = varname;
    for (i=0; i<newvars; i++) {
	pdinfo->varname[v+i] = malloc(VNAMELEN);
	if (pdinfo->varname[v+i] == NULL) return E_ALLOC;
	pdinfo->varname[v+i][0] = '\0';
    }

    if (pdinfo->varinfo != NULL) {
	varinfo = realloc(pdinfo->varinfo, (v + newvars) * sizeof *varinfo);
	if (varinfo == NULL) return E_ALLOC;
	else pdinfo->varinfo = varinfo;
	for (i=0; i<newvars; i++) {
	    pdinfo->varinfo[v+i] = malloc(sizeof **varinfo);
	    if (pdinfo->varinfo[v+i] == NULL) return E_ALLOC;
	    gretl_varinfo_init(pdinfo->varinfo[v+i]);
	}
    }

    vector = realloc(pdinfo->vector, (v + newvars));
    if (vector == NULL) return E_ALLOC;

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
    else pdinfo->varname = varname;

    pdinfo->varname[v] = malloc(VNAMELEN);
    if (pdinfo->varname[v] == NULL) return E_ALLOC;
    pdinfo->varname[v][0] = '\0';

    if (pdinfo->varinfo != NULL) {
	varinfo = realloc(pdinfo->varinfo, (v + 1) * sizeof *varinfo);
	if (varinfo == NULL) return E_ALLOC;
	else pdinfo->varinfo = varinfo;

	pdinfo->varinfo[v] = malloc(sizeof **varinfo);
	if (pdinfo->varinfo[v] == NULL) return E_ALLOC;
	gretl_varinfo_init(pdinfo->varinfo[v]);
    }

    vector = realloc(pdinfo->vector, (v + 1));
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;

    pdinfo->vector[v] = 0;

    pdinfo->v += 1;

    return 0;
}

/* ......................................................  */

int varnum_from_string (const char *str, DATAINFO *pdinfo)
{
    int varno;
    char *test;

    errno = 0;

    strtol(str, &test, 10);
    if (*test != '\0' || !strcmp(str, test) || errno == ERANGE) {
        return -1;
    }

    varno = atoi(str);

    if (varno <= 0 || varno >= pdinfo->v) {
	return -1;
    } 
    
    return varno;
}

/* ......................................................  */

int dataset_drop_listed_vars (const int *list, double ***pZ, 
			      DATAINFO *pdinfo, int *renumber)
{
    double **newZ;
    char **varname;
    char *vector;
    VARINFO **varinfo;
    int oldv = pdinfo->v, vmax = pdinfo->v;
    int i, v, ndel = 0; 

    *renumber = 0;

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

    for (v=1; v<vmax; v++) {
	if ((*pZ)[v] == NULL) {
	    int gap = 1;

	    for (i=v+1; i<vmax; i++) {
		if ((*pZ)[i] == NULL) gap++;
		else break;
	    }
	    if (i < vmax) {
		*renumber = 1;
		vmax -= gap;
		for (i=v; i<vmax; i++) {
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

    varname = realloc(pdinfo->varname, (oldv - ndel) * sizeof *varname);
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, (oldv - ndel) * sizeof *vector);
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;

    varinfo = realloc(pdinfo->varinfo, (oldv - ndel) * sizeof *varinfo);
    if (varinfo == NULL) return E_ALLOC;
    else pdinfo->varinfo = varinfo;

    newZ = realloc(*pZ, (oldv - ndel) * sizeof *newZ); 
    if (newZ == NULL) return E_ALLOC;
    else *pZ = newZ;

    pdinfo->v -= ndel;

    return 0;
}

/* drop specified number of variables at the end of the dataset */

int dataset_drop_vars (int delvars, double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char *vector;
    VARINFO **varinfo;
    int i, v = pdinfo->v;   

    if (delvars <= 0) return 0;

    if (pdinfo->v <= 1) return E_DATA;

    for (i=v-delvars; i<v; i++) {
	if (pdinfo->varname[i] != NULL) free(pdinfo->varname[i]);
	if (pdinfo->varinfo[i] != NULL) free(pdinfo->varinfo[i]);
	if ((*pZ)[i] != NULL) free((*pZ)[i]);
    }

    newZ = realloc(*pZ, (v - delvars) * sizeof *newZ); 
    if (newZ == NULL) return E_ALLOC;
    else *pZ = newZ;
        
    varname = realloc(pdinfo->varname, (v - delvars) * sizeof *varname);
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, (v - delvars) * sizeof *vector);
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;

    varinfo = realloc(pdinfo->varinfo, (v - delvars) * sizeof *varinfo);
    if (varinfo == NULL) return E_ALLOC;
    else pdinfo->varinfo = varinfo;

    pdinfo->v -= delvars;

    return 0;
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
    int i;
    double *xx;

    if (n == 0 || src == NULL) return NULL;

    xx = malloc(n * sizeof *xx);
    if (xx == NULL) return NULL;

    for (i=0; i<n; i++) xx[i] = src[i];

    return xx;
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
	if (t1 < maxlag) t1 = maxlag; 
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	zz = 0.0;
	if (ar) 
	    for (k=1; k<=arlist[0]; k++) {
	    xx = (*pZ)[yno][t - arlist[k]];
	    zr = pmod->arinfo->rho[k];
	    if (na(xx)) {
		if (zr == 0) continue;
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
    int depvar, t;
    int t1 = pmod->t1, t2 = pmod->t2, n = pdinfo->n;
    FITRESID *fr;

    if (pmod->ci == ARMA) {
	depvar = pmod->list[4];
    } else {
	depvar = pmod->list[1];
    }

    if (pmod->data != NULL) {
	t2 += get_misscount(pmod);
    }

    fr = fit_resid_new(n, 0);
    if (fr == NULL) return NULL;

    fr->sigma = pmod->sigma;

    for (t=0; t<n; t++) {
	fr->actual[t] = (*pZ)[depvar][t];
	fr->fitted[t] = pmod->yhat[t];
    }

    if (isdummy(fr->actual, 0, n) > 0) {
	fr->pmax = get_precision(fr->fitted, n, 8);
    } else {
	fr->pmax = get_precision(fr->actual, n, 8);
    }
    
    strcpy(fr->depvar, pdinfo->varname[depvar]);
    
    fr->t1 = t1;
    fr->t2 = t2;
    fr->nobs = pmod->nobs;

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
    int *list = NULL;
    int ft1, ft2;
    int i, j, k, t, nfcast, n_est, nv;
    int yno = pmod->list[1];
    char t1str[OBSLEN], t2str[OBSLEN];

    fr = fit_resid_new(0, 1); 
    if (fr == NULL) return NULL;

    if (pmod->ci != OLS) {
	fr->err = E_OLSONLY;
	return fr;
    }

    /* bodge (reject in case of subsampled data) */
    if (pmod->data != NULL) {
	fr->err = E_DATA;
	return fr;
    }

    /* parse dates */
    if (sscanf(str, "%*s %8s %8s", t1str, t2str) != 2) {
	fr->err = E_OBS;
	return fr;
    }

    ft1 = dateton(t1str, pdinfo);
    ft2 = dateton(t2str, pdinfo);
    if (ft1 < 0 || ft2 < 0 || ft2 <= ft1) {
	fr->err = E_OBS;
	return fr;
    }

    /* move endpoints if there are missing vals for the
       independent variables */
    fcast_adjust_t1_t2(pmod, (const double **) *pZ, &ft1, &ft2);

    /* number of obs for which forecasts will be generated */
    nfcast = ft2 - ft1 + 1;

    if (allocate_fit_resid_arrays(fr, nfcast, 1)) {
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

    list = malloc((finfo->v + 1) * sizeof *list);
    if (list == NULL) {
	fr->err = E_ALLOC;
	goto fcast_bailout;
    }

    list[0] = finfo->v;
    list[1] = 1;
    list[2] = 0;
    for (i=3; i<=list[0]; i++) {
	list[i] = i - 1;
    }

    gretl_model_init(&fmod);

    /* loop across the observations for which we want forecasts
       and standard errors */

#ifdef FCAST_DEBUG
    printf("get_fcast_with_errs: ft1=%d, ft2=%d, pmod->t1=%d, pmod->t2=%d\n",
	   ft1, ft2, pmod->t1, pmod->t2);
#endif

    for (k=0; k<nfcast; k++) {
	/* form modified indep vars: original data minus the values
	   to be used for the forecast */
	for (i=3; i<=list[0]; i++) {
	    j = (pmod->ifc)? pmod->list[i] : pmod->list[i-1];
	    for (t=0; t<finfo->n; t++) {
		fZ[i-1][t] = (*pZ)[j][t + pmod->t1] 
		    - (*pZ)[j][k + ft1]; 
	    }
	}
	clear_model(&fmod);
	fmod = lsq(list, &fZ, finfo, OLS, OPT_A, 0.0);
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
	fr->actual[k] = (*pZ)[pmod->list[1]][k + ft1];
    }

    clear_model(&fmod);

    fr->tval = tcrit95(pmod->dfd);
    strcpy(fr->depvar, pdinfo->varname[yno]);

    fr->t1 = ft1;
    fr->t2 = ft2;
    fr->nobs = ft2 - ft1 + 1;
    fr->df = pmod->dfd;

 fcast_bailout:
    free_Z(fZ, finfo);
    free(list);
    clear_datainfo(finfo, CLEAR_FULL);
    free(finfo);

    return fr;
}

/* ........................................................... */

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn,
		     PATHS *ppaths, int plot)
{
    FITRESID *fr;
    int err;

    fr = get_fcast_with_errs(str, pmod, pZ, pdinfo, prn);
    if (fr == NULL) {
	return E_ALLOC;
    } 
    else if (fr->err) {
	err = fr->err;
	free_fit_resid(fr);
	return err;
    }

    err = text_print_fcast_with_errs (fr, pZ, pdinfo, prn,
				      ppaths, plot);

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
    if (cmd.list == NULL || cmd.param == NULL) 
	return 1;

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
	err = hilu_corc(&rho, cmd.list, pZ, pdinfo, 
			NULL, 1, cmd.ci, &prn);
	if (!err) {
	    *tmpmod = lsq(cmd.list, pZ, pdinfo, cmd.ci, 0, rho);
	}
	break;
    case HCCM:
	*tmpmod = hccm_func(cmd.list, pZ, pdinfo);
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
    case OLS:
    case WLS:
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

    if (cmd.list) free(cmd.list);
    if (cmd.param) free(cmd.param);

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
	    pdinfo->time_series = STACKED_CROSS_SECTION; 
	    panel = STACKED_CROSS_SECTION;
	} else {
	    pdinfo->time_series = STACKED_TIME_SERIES; 
	    panel = STACKED_TIME_SERIES;
	}
    }

    return panel;
}

/* ........................................................... */

int get_panel_structure (DATAINFO *pdinfo, int *nunits, int *T)
{
    int err = 0;

    if (pdinfo->time_series == STACKED_TIME_SERIES) {
        *nunits = pdinfo->n / pdinfo->pd;
        *T = pdinfo->pd;
    } 
    else if (pdinfo->time_series == STACKED_CROSS_SECTION) {
        char Tstr[8];

        if (sscanf(pdinfo->endobs, "%[^:]:%d", Tstr, nunits) != 2) {
            err = 1;
        } else { 
            *T = atoi(Tstr);
	}
    } 
    else err = 1;

    return err;
}

/* ........................................................... */

int set_panel_structure (unsigned char flag, DATAINFO *pdinfo, PRN *prn)
{
    int nunits, T;
    int old_ts = pdinfo->time_series;

    if (pdinfo->pd == 1) {
	pputs(prn, _("The current data frequency, 1, is not "
		"compatible with panel data.\nPlease see the 'setobs' "
		"command.\n"));
	return 1;
    }

    if (flag == 'c') {
	pdinfo->time_series = STACKED_CROSS_SECTION;
    } else {
	pdinfo->time_series = STACKED_TIME_SERIES;
    }

    if (get_panel_structure(pdinfo, &nunits, &T)) {
	pputs(prn, _("Failed to set panel structure\n"));
	pdinfo->time_series = old_ts;
	return 1;
    } else {
	pprintf(prn, _("Panel structure set to %s\n"),
		(pdinfo->time_series == STACKED_CROSS_SECTION)? 
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

    if ((pdinfo->t2 - pdinfo->t1 + 1) % pdinfo->pd)
        return 0;

    if (sscanf(pdinfo->endobs, "%[^:]:%s", unit, period) == 2) {
        if (atoi(period) != pdinfo->pd) return 0;
    } else {
        return 0;
    }

    return 1;
}

/* ........................................................... */

double get_xvalue (int i, double **Z, const DATAINFO *pdinfo)
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
    int i, n = mpvals->ncoeff + 1;

    mpvals->varnames = malloc(n * sizeof *mpvals->varnames);
    if (mpvals->varnames == NULL) return 1;

    for (i=0; i<n; i++) {
	mpvals->varnames[i] = malloc(12);
	if (mpvals->varnames[i] == NULL) {
	    free_mp_varnames(mpvals);
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
    if (fr == NULL) return NULL;

    fr->err = 0;
    fr->t1 = 0;
    fr->t2 = 0;
    fr->nobs = 0;

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
    if (*targ == NULL) return 1;
    
    (*targ)[0] = n;
    for (i=1; i<=n; i++) (*targ)[i] = src[i];

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

static int real_gretl_spawn (const char *cmdline, int verbose)
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
	if (verbose) {
	    fprintf(stderr, "gretl_spawn: '%s'\n", error->message);
	}
	g_error_free(error);
	ret = 1;
    } else if (errout && *errout) {
	strcpy(gretl_errmsg, errout);
	if (verbose) {
	    fprintf(stderr, "stderr: '%s'\n", errout);
	}
	if (!font_not_found(errout)) {
	    ret = 1;
	}
    } else if (status != 0) {
	sprintf(gretl_errmsg, "%s\n%s", 
		_("Command failed"),
		sout);
	if (verbose) {
	    fprintf(stderr, "status=%d: '%s'\n", status, sout);
	}
	ret = 1;
    }

    if (errout != NULL) g_free(errout);
    if (sout != NULL) g_free(sout);

    if (ret && verbose) {
	fprintf(stderr, "Failed command: '%s'\n", cmdline);
    } 

    return ret;
}

#endif

#if !defined(WIN32) && !defined(GRETL_GLIB)

static int real_gretl_spawn (const char *cmdline, int verbose)
{
    int err;

    errno = 0;

    signal(SIGCHLD, SIG_DFL);

    err = system(cmdline);
    if (err && verbose) {
	fprintf(stderr, "Failed command: '%s'\n", cmdline);
	perror(NULL);
    }

    return err;
}

#endif

#ifndef WIN32

int gretl_spawn (const char *cmdline)
{
    return real_gretl_spawn(cmdline, 1);
}

int gretl_spawn_quiet (const char *cmdline)
{
    return real_gretl_spawn(cmdline, 0);
}

#endif
