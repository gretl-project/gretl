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
#include "internal.h"

#include <dirent.h>
#include <unistd.h>
#include <errno.h>

#ifndef WIN32
# include <glib.h>
# if GLIB_CHECK_VERSION(2,0,0)
#  define GLIB2
#  include <signal.h>
# endif /* GLIB_CHECK_VERSION */
#endif /* ! WIN32 */

static DIR *dir;

static int _pdton (int pd);
static int allocate_fit_resid_arrays (FITRESID *fr, int n, int errs);

/* .......................................................... */

static int path_append (char *file, const char *path)
{
    char temp[MAXLEN];
    int n, pathlen = strlen(file) + strlen(path) + 1;

    if (pathlen > MAXLEN) return 1;
    strcpy(temp, path);
    n = strlen(temp);
    if (temp[n - 1] != SLASH && n < MAXLEN - 1) {
	temp[n] = SLASH;
	temp[n + 1] = '\0';
    }
    strcat(temp, file);
    strcpy(file, temp);
    return 0;
}

/* .......................................................  */

double _corr (int n, const double *zx, const double *zy)
/*
        returns the simple correlation coefficient between the the
        arrays zx and zy, for the n observations 0 to n-1.  returns
        -999 if square root argument is invalid or no of observations
        is zero 
*/
{
    register int i;
    int nn;
    double sx, sy, sxx, syy, sxy, den, zxi, zyi, zxbar, zybar;

    if (n == 0) return NADBL;
    if (_isconst(0, n-1, zx) || _isconst(0, n-1, zy)) return NADBL;

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
    sxx = syy = sxy = 0.0;

    for (i=0; i<n; ++i) {
        zxi = zx[i];
        zyi = zy[i];
        if (na(zxi) || na(zyi)) continue;
        sx = zxi - zxbar;
        sy = zyi - zybar;
        sxx = sxx + (sx*sx);
        syy = syy + (sy*sy); 
        sxy = sxy + (sx*sy);
    }

    if (sxy != 0.0) {
        den = sxx * syy;
        if (den > 0.0) return sxy/sqrt(den);
        else return NADBL;
    }
    else return 0.0;
}

/* .......................................................  */

double _covar (int n, const double *zx, const double *zy)
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

    if (pd == 1) 
	return ((double) (ysd + nt));  
    pp = (nt) % pd + _pdton(pd) * (sd0 - ysd) + .5;
    if (pp != pd)  {
        yy = ysd + (nt)/pd  + (pp/pd) + .5;
        yp = pp % pd;
    }  else {
        yy = ysd + (nt)/pd + .5;
        yp = pp;
    }
    dd = (pd < 10)? 0.1: 0.01;
    return (yy + yp * dd);
}

/**
 * ijton:
 * @i: row number (1-based)
 * @j: column number (1-based)
 * @nrows: number of rows (and columns) in symmetric matrix.
 *
 * Given a (row, column) reference into a symmetric 2-dimensional 
 * matrix A, finds the 0-based index into a 1-dimensional array 
 * x composed of the non-redundant elements of A.
 *
 * E.g. for the 3 x 3 case with 6 non-redundant elements, 0 to 5,
 *
 *    A(1,1) = x[0]  A(1,2) = x[1]  A(1,3) = x[2]
 *    A(2,1) = x[1]  A(2,2) = x[3]  A(2,3) = x[4]
 *    A(3,1) = x[2]  A(3,2) = x[4]  A(3,3) = x[5]
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

    idx = nrows * (i - 1) + j - i - ((i - 2) * (i - 1)/2);
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

#ifdef notdef
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
    if (m == 0)
	fprintf(stderr, "\nztox: No valid observations for variable %s\n", 
		pdinfo->varname[i]);
    else if (m < pdinfo->t2 - pdinfo->t1 + 1)
	fprintf(stderr, "\nztox: Dropped missing obs for var %s\n",
		pdinfo->varname[i]);
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

int isdummy (double *x, int t1, int t2)
{
    int t, m = 0;
    double xx;

    for (t=t1; t<=t2; t++) {
	xx = x[t];
	if (floatneq(xx, 0.0) && floatneq(xx, 1.0)) {
	    return 0;
	}
	if (floateq(xx, 1.0)) m++;
    }

    if (m < t2 - t1 + 1) return m;

    return 0;
} 

/* ........................................................  */

int _iszero (int t1, int t2, const double *x)
/*  checks whether all obs are zero for variable x from t1 to t2 */
{
    int t;
    double xx, sum = 0.0;

    for (t=t1; t<=t2; t++) {
        xx = x[t];
        sum = sum + xx*xx;
    }
    if (floateq(sum, 0.0)) return 1;
    else return 0;
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

int _isconst (int t1, int t2, const double *x)
{
    int t;
    double xx = x[t1];

    for (t=t1+1; t<=t2; t++) if (floatneq(x[t], xx)) return 0;
    return 1;
}

/* ............................................................  */

double _esl_mean (int t1, int t2, const double *x)
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

void _minmax (int t1, int t2, const double zx[], 
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

/* .......................................................     */

static int _pdton (int pd)
{
    if (pd == 1) return 1;
    if (pd < 10) return 10;
    return 100;
}

/* ..........................................................  */

int _hasconst (const int *list)
/* check if a var list contains a constant (variable with ID
   number 0) */
{
    int i;

    for (i=2; i<=list[0]; i++) 
        if (list[i] == 0) return i;

    return 0;
}

/* ...................................................... */

int _compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
}

/* .............................................................  */

double _esl_stddev (int t1, int t2, const double *x)
/*  returns standard deviation of array x from t1 through t2
    return -999 if square root argument is invalid
    or there are no observations
*/
{
    double xx;

    xx = _esl_variance(t1, t2, x);
    if (na(xx)) return xx;
    return sqrt(xx);
}

/* .............................................................  */

double _esl_variance (int t1, int t2, const double *x)
{
    int n;
    register int i;
    double sumsq, xx, xbar;

    n = t2 - t1 + 1;
    if (n == 0) return NADBL;
    xbar = _esl_mean(t1, t2, x);
    if (na(xbar)) return NADBL;
    sumsq = 0.0;
    for (i=t1; i<=t2; i++) {
	if (!na(x[i])) {
	    xx = x[i] - xbar;
	    sumsq += xx*xx;
	} else n--;
    }
    sumsq = (n > 1)? sumsq/(n - 1) : 0.0;
    if (sumsq >= 0) return sumsq;
    else return NADBL;
}

/* .............................................................  */

double _esl_sst (int t1, int t2, const double *x)
{
    register int i;
    double sumsq, xx, xbar;

    if (t2 - t1 + 1 == 0) return NADBL;
    xbar = _esl_mean(t1, t2, x);
    if (na(xbar)) return NADBL;
    sumsq = 0.0;
    for (i=t1; i<=t2; i++) {
	if (!na(x[i])) {
	    xx = x[i] - xbar;
	    sumsq += xx*xx;
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

/* ....................................................... */

void gretl_aic_etc (MODEL *pmod)
/*
    Compute model selection criteria -- needs nobs, ncoeff
    and ess from model.
*/
{
    double zz, zx, ersq, zn;
    double ess = pmod->ess;
    int nobs = pmod->nobs, ncoeff = pmod->ncoeff;

    zz = (double) (nobs - ncoeff);
    pmod->criterion[0] = ess/zz;
    ersq = ess/nobs;
    pmod->criterion[2] = ersq * (nobs + ncoeff) / zz;
    zz = 2.0 * ncoeff / nobs;
    pmod->criterion[1] = ersq * exp(zz);
    pmod->criterion[5] = ersq * (1.0 + zz);
    pmod->criterion[7] = ((1-zz) > 0.0)? ersq/(1-zz) : NADBL;
    zn = (double) nobs;
    zx = log(zn);
    pmod->criterion[3] = ersq * pow(zx, zz);
    zz = (double) ncoeff;
    zz = zz/zn;
    pmod->criterion[4] = ersq * pow(zn, zz);
    zz = 1.0 - zz;
    pmod->criterion[6] = ersq/(zz*zz);
}

/* ....................................................... */

void _criteria (const double ess, int nobs, int ncoeff, 
		PRN *prn)
{
    double zz, zx, ersq, zn;
    double criterion[8];
    
    zz = (double) (nobs - ncoeff);
    criterion[0] = ess/zz;
    ersq = ess/nobs;
    criterion[2] = ersq * (nobs + ncoeff) / zz;
    zz = 2.0 * ncoeff/nobs;
    criterion[1] = ersq * exp(zz);
    criterion[5] = ersq * (1.0 + zz);
    criterion[7] = ((1-zz) > 0.0)? ersq/(1-zz) : NADBL;
    zn = (double) nobs;
    zx = log(zn);
    criterion[3] = ersq * pow(zx, zz);
    zz = (double) ncoeff;
    zz = zz/zn;
    criterion[4] = ersq * pow(zn, zz);
    zz = 1.0 - zz;
    criterion[6] = ersq/(zz*zz);

    pprintf(prn, _("Using ess = %f, %d observations, %d coefficients\n"), 
	   ess, nobs, ncoeff);
    pputs(prn, _("\nMODEL SELECTION STATISTICS\n\n"));	
    pprintf(prn, "SGMASQ    %13g     AIC       %13g     FPE       %12g\n"
	    "HQ        %13g     SCHWARZ   %13g     SHIBATA   %12g\n"
	    "GCV       %13g",
	    criterion[0], criterion[1], 
	    criterion[2], criterion[3], 
	    criterion[4], criterion[5], criterion[6]);
    if (criterion[7] > 0.0) pprintf(prn, "     RICE      %13g\n", 
					  criterion[7]);
    else pputs(prn, "     RICE          undefined\n");
    pputs(prn, "\n");
}

/* ....................................................... */

int _adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		  double **Z, int *misst)
     /* drop first/last observations from sample if missing obs 
	encountered -- also check for missing vals within the
        remaining sample */
{
    int i, t, dwt = 0, t1min = *t1, t2max = *t2;
    double xx;

    if (pmod != NULL && pmod->wt_dummy) dwt = pmod->nwt;

    for (i=1; i<=list[0]; i++) {
	for (t=t1min; t<t2max; t++) {
	    xx = Z[list[i]][t];
	    if (dwt) xx *= Z[dwt][t];
	    if (na(xx)) t1min += 1;
	    else break;
	}
    }
    for (i=1; i<=list[0]; i++) {
	for (t=t2max; t>t1min; t--) {
	    xx = Z[list[i]][t];
	    if (dwt) xx *= Z[dwt][t];
	    if (na(xx)) t2max -= 1;
	    else break;
	}
    } 
    if (misst != NULL) {
	for (i=1; i<=list[0]; i++) {
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
	if (count) 
	    pprintf(prn, _("%s: set %d observations to \"missing\"\n"), 
		    pdinfo->varname[list[i]], count);
	else 
	    pprintf(prn, _("%s: Didn't find any matching observations\n"),
		    pdinfo->varname[list[i]]);
    }
}

/**
 * set_obs:
 * @line: command line.
 * @pdinfo: data information struct.
 * @opt: 's' for stacked time-series, 'c' for stacked cross-section.
 * 
 * Impose a time-series or panel interpretation on a data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int set_obs (char *line, DATAINFO *pdinfo, unsigned char opt)
{
    int pd, pos, i, len, dc = 0, bad = 0;
    char stobs[9], endobs[9], endbit[7], *p;
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

    /* special case: dated daily data */
    if ((pd == 5 || pd == 7) && strstr(stobs, "/")) {
	ed0 = get_epoch_day(stobs);
	if (ed0 < 0) {
	    sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
	    return 1;
	}
	else pdinfo->sd0 = (double) ed0;
    } else {
	/* is stobs acceptable? */
	len = strlen(stobs);
	for (i=0; i<len; i++) {
	    if (stobs[i] != '.' && !isdigit((unsigned char) stobs[i])) {
		bad = 1;
		break;
	    }
	    if (stobs[i] == '.') dc++;
	}
	if (bad || dc > 1) {
	    sprintf(gretl_errmsg, _("starting obs '%s' is invalid"), stobs);
	    return 1;
	}
	pos = dotpos(stobs);
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
	if ((pd > 1 && pd < 10 && strlen(stobs + pos) != 2) ||
	    (pd >= 10 && pd < 100 && strlen(stobs + pos) != 3)) {
	    sprintf(gretl_errmsg, _("starting obs '%s' is incompatible with "
		    "frequency"), stobs);
	    return 1;
	}
	if (pd > 1) {
	    strcpy(endbit, stobs + pos + 1);
	    dc = atoi(endbit);
	    if (dc < 0 || dc > pd) {
		sprintf(gretl_errmsg, 
			_("starting obs '%s' is incompatible with frequency"), 
			stobs);
		return 1;
	    }	    
	}
    }

    /* adjust data info struct */
    pdinfo->pd = pd;

    if (ed0 == 0L) {
	pdinfo->sd0 = atof(stobs);
    } else {
	pdinfo->time_series = TIME_SERIES;
    }

    ntodate(pdinfo->stobs, 0, pdinfo);
    ntodate(endobs, pdinfo->n - 1, pdinfo);
    strcpy(pdinfo->endobs, endobs);

    if (opt == 's') pdinfo->time_series = STACKED_TIME_SERIES;
    else if (opt == 'c') pdinfo->time_series = STACKED_CROSS_SECTION;
    else if (pdinfo->sd0 >= 1.0) 
        pdinfo->time_series = TIME_SERIES; /* but might be panel? */
    else pdinfo->time_series = 0;

    /* and report? */
#if 0
    fprintf(stderr, I_("setting data frequency = %d\n"), pd);
    if (pdinfo->n > 0) {
	fprintf(stderr, I_("data range: %s - %s\n"), pdinfo->stobs, pdinfo->endobs);
    }
#endif

    return 0;
}

#ifdef WIN32
char *unslash (const char *src)
{
    size_t n = strlen(src);
    char *dest = malloc(n);

    if (dest != NULL) strncpy(dest, src, n-1);
    return dest;
}
#endif

/* .......................................................... */

static int get_subdir (const char *topdir, int first, char *fname)
{
    DIR *try;
    struct dirent *dirent;
#ifdef WIN32
    char *tmp = unslash(topdir);

    if (tmp == NULL) return -1;
#endif

    if (first) {
#ifdef WIN32
	if ((dir = opendir(tmp)) == NULL) {
	    free(tmp);
	    return -1;
	}
	free(tmp);
#else
	if ((dir = opendir(topdir)) == NULL) return -1;
#endif
    } else {
	if ((dirent = readdir(dir)) == NULL) {
	    closedir(dir);
	    dir = NULL;
	    return -1;
	} else {
	    if (strcmp(dirent->d_name, ".") == 0 ||
		strcmp(dirent->d_name, "..") == 0) return 0;
	    strcpy(fname, topdir);
	    strcat(fname, dirent->d_name);
	    if ((try = opendir(fname)) != NULL) {
		closedir(try);
		return 1;
	    } else return 0;
	}
    }
    return 1;
}

/* .......................................................... */

static char *search_dir (char *filename, const char *topdir, 
			 int recurse)
{
    FILE *test;
    int got = 0;
    char origfile[MAXLEN], trypath[MAXLEN];

    strcpy(origfile, filename);

    if (path_append(filename, topdir) == 0) {
#ifdef PATH_DEBUG
	fprintf(stderr, I_("Trying %s\n"), filename);
#endif
	test = fopen(filename, "r");
	if (test != NULL) {
	    fclose(test);
	    return filename;
	}
	if (!recurse) return NULL;
	if (get_subdir(topdir, 1, trypath) > 0) {
	    while ((got = get_subdir(topdir, 0, trypath)) >= 0) {
		strcpy(filename, origfile);
		if (got && path_append(filename, trypath) == 0) {
#ifdef PATH_DEBUG
		    fprintf(stderr, I_("Trying %s\n"), filename);
#endif
		    test = fopen(filename, "r");
		    if (test != NULL) {
			fclose(test);
			return filename;
		    }		    
		}
	    }
	}
    }
    return NULL;
}

/**
 * addpath:
 * @fname: initially given file name.
 * @ppaths: path information struct.
 * @script: if non-zero, suppose the file is a command script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.  Usually called by getopenfile().
 *
 * Returns: the full name of the file that was found, or NULL if no
 * file could be found.
 */

char *addpath (char *fname, PATHS *ppaths, int script)
{
    char orig[MAXLEN];
    char *thisdir, *tmp = fname;
    FILE *test;

    strcpy(orig, fname);

    /* try opening filename as given */
    test = fopen(fname, "r");
    if (test != NULL) { 
	fclose(test); 
	/* if a relative path was given, convert it to absolute */
#ifdef WIN32
	if (fname[1] == ':') return fname;
#endif
	if (fname[0] != SLASH && (thisdir = malloc(MAXLEN)) != NULL) {
	    int i = 0;

	    if (getcwd(thisdir, MAXLEN-1) != NULL) {
#ifdef WIN32		
		lower(thisdir); /* hmmm */
#endif
		if (strstr(fname, thisdir) == NULL) {
		    strcpy(fname, thisdir);
		    strcat(fname, SLASHSTR);
		    if (orig[0] == '.' && orig[1] == SLASH &&
			strlen(orig) > 2) i = 2;
		    strcat(fname, orig + i);
		}
	    }
	    free(thisdir);
	} /* end conversion to absolute path */
	return fname;
    } else {  
	/* not able to open file as given */
	if (fname[0] == '.' || fname[0] == SLASH)
	    return NULL;
    }

    /* try looking where script was found */
    if (ppaths->currdir[0]) {
	if ((fname = search_dir(fname, ppaths->currdir, 0))) 
	    return fname;
    }

    fname = tmp;
    strcpy(fname, orig);
    if (!script) {
	/* if it's a data file we want, try system data dir */
	if ((fname = search_dir(fname, ppaths->datadir, 1))) 
	    return fname;
    } else {
	/* for a script, try system script dir */
	if ((fname = search_dir(fname, ppaths->scriptdir, 1))) 
	    return fname;
    }

    /* try looking in user's dir */
    fname = tmp;
    strcpy(fname, orig);
    if ((fname = search_dir(fname, ppaths->userdir, 1))) 
	return fname;

    fname = tmp;
    strcpy(fname, orig);

    return NULL;
}

static int get_quoted_filename (const char *line, char *fname)
{
    char *p;
    int quote = '"';

    p = strchr(line, quote);
    if (p == NULL) {
	quote = '\'';
	p = strchr(line, quote);
	if (p == NULL) return 0;
    }

    if (p != NULL) {
	char *q = strrchr(line, quote);

	if (q == NULL) return 0;
	else {
	    size_t len = q - p;

	    if (len > 0) {
		*fname = 0;
		strncat(fname, p+1, len-1);
	    } else
		return 0;
	}
    }
    return 1;
}


/**
 * getopenfile:
 * @line: command line (e.g. "open foo").
 * @fname: filename to be filled out.
 * @ppaths: path information struct.
 * @setpath: if non-zero, set @ppaths->currdir based on the file
 * that is found (if any).
 * @script: if non-zero, suppose the file is a command script.
 * 
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened.
 *
 * Returns: 0 on successful parsing of @line, 1 on error.
 */

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 int setpath, int script)
{
    int spos, n;

    /* get the initial filename off the command line */
    if (get_quoted_filename(line, fname)) return 0; 
    if (sscanf(line, "%*s %s", fname) != 1) return 1;
    /* try a basic path search on this filename */
    addpath(fname, ppaths, script);
    if (addpath != NULL && setpath) {
	ppaths->currdir[0] = '.';
	ppaths->currdir[1] = SLASH;
	ppaths->currdir[2] = '\0';
	spos = slashpos(fname);
	if (spos) {
	    strncpy(ppaths->currdir, fname, (size_t) spos);
	    n = strlen(ppaths->currdir);
	    ppaths->currdir[n] = SLASH;
	    ppaths->currdir[n+1] = '\0';
	}
    }
    if (dir != NULL) {  /* dir is static, declared outside of funcs */
	closedir(dir);     
	dir = NULL;
    }
    return 0;
}

/* .......................................................... */

static unsigned char isflag (unsigned char c)
{
    if (c == 'o' || c == 'c' || c == 'i' || c == 'm' ||
	c == 'r' || c == 's' || c == 't' || c == 'l' ||
	c == 'a' || c == 'z' || c == 'w' || c == 'q') {
	return c;
    }
    return 0;
}

/* .......................................................... */

int catchflag (char *line, unsigned char *oflag)
     /* check for "-<char>" in line: if found, chop it out and set
	oflag value accordingly.  
	Strip trailing semicolon while we're at it. 
     */
{
    int i, opt, n = strlen(line);

    *oflag = 0;

    if (n < 2) return 0;

    /* to enable reading of trad. esl input files */
    if (line[n-2] == ';' && isspace(line[n-1])) {
	line[n-2] = '\0';
	n = strlen(line);
    } else if (line[n-1] == ';') {
	line[n-1] = '\0';
	n = strlen(line);
    }

    for (i=4; i<n-1; i++) {
	if (line[i] == '-' && 
	    isspace((unsigned char) line[i-1]) && 
	    (opt = isflag(line[i+1])) &&
	    (i+2 == n || isspace((unsigned char) line[i+2]))) {
		*oflag = opt;
		_delete(line, i, 2);
		return 1;
	    }
    }

    return 0;
}

/* .......................................................... */

int _list_dups (const int *list, int ci)
{
    int i, j, start = 2;

    if (ci == ARCH) start = 3;

    if (ci == TSLS || ci == AR || ci == SCATTERS || ci == MPOLS) {
	for (i=2; i<list[0]; i++) {
	    if (list[i] == 999) {
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

/* .......................................................... */

void _init_model (MODEL *pmod, const DATAINFO *pdinfo)
{
    int i;

    if (pmod == NULL) return;
    pmod->list = NULL;
    pmod->subdum = NULL;
    if (pdinfo != NULL) {
	pmod->smpl.t1 = pdinfo->t1;
	pmod->smpl.t2 = pdinfo->t2;
    }
    pmod->coeff = NULL;
    pmod->sderr = NULL;
    pmod->yhat = NULL;
    pmod->uhat = NULL;
    pmod->xpx = NULL;
    pmod->vcv = NULL;
    pmod->arinfo = NULL;
    pmod->slope = NULL;
    pmod->name = NULL;
    pmod->params = NULL;
    pmod->ntests = 0;
    pmod->tests = NULL;
    pmod->data = NULL;
    pmod->dataset = NULL;
    pmod->errcode = 0;
    pmod->ci = 0;
    pmod->ifc = 0;
    pmod->rho_in = 0.0;
    pmod->aux = AUX_NONE;
    *gretl_msg = '\0';
    
    for (i=0; i<8; i++) {
	pmod->criterion[i] = 0.0;
    }
}

/**
 * gretl_model_new:
 * @pdinfo: pointer to data information struct.
 * 
 * Allocates memory for a gretl MODEL struct and initializes the struct.
 *
 * Returns: pointer to #MODEL (or NULL if allocation fails).
 */

MODEL *gretl_model_new (const DATAINFO *pdinfo)
{
    MODEL *pmod = malloc(sizeof *pmod);

    _init_model(pmod, pdinfo);
    return pmod;
}

/* ........................................................... */

void exchange_smpl (MODEL *pmod, DATAINFO *pdinfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;

    pdinfo->t1 = pmod->smpl.t1;
    pdinfo->t2 = pmod->smpl.t2;
    pmod->smpl.t1 = t1;
    pmod->smpl.t2 = t2;
}

/* .......................................................... */

static void clear_ar_info (MODEL *pmod)
{
    if (pmod->arinfo->arlist) free(pmod->arinfo->arlist);
    if (pmod->arinfo->rho) free(pmod->arinfo->rho);
    if (pmod->arinfo->sderr) free(pmod->arinfo->sderr);

    free(pmod->arinfo);
    pmod->arinfo = NULL;
}

void debug_print_model_info (const MODEL *pmod, const char *msg)
{
    fprintf(stderr, "%s:\n"
	    " pmod = %p\n"
	    " pmod->list = %p\n"
	    " pmod->subdum = %p\n"
	    " pmod->coeff = %p\n"
	    " pmod->sderr = %p\n"
	    " pmod->yhat = %p\n"
	    " pmod->uhat = %p\n"
	    " pmod->xpx = %p\n"
	    " pmod->vcv = %p\n"
	    " pmod->name = %p\n"
	    " pmod->params = %p\n"
	    " pmod->arinfo = %p\n"
	    " pmod->slope = %p\n"
	    " pmod->tests = %p\n"
	    " pmod->data = %p\n", msg,
	    pmod, pmod->list, pmod->subdum, pmod->coeff,
	    pmod->sderr, pmod->yhat, pmod->uhat, pmod->xpx,
	    pmod->vcv, pmod->name, pmod->params, pmod->arinfo, 
	    pmod->slope, pmod->tests, pmod->data);
}

/* .......................................................... */

void clear_model (MODEL *pmod, const DATAINFO *pdinfo)
{
    if (pmod != NULL) {
#if 0
	debug_print_model_info(pmod, "Doing clear_model");
#endif
	if (pmod->list) free(pmod->list);
	if (pmod->subdum) free(pmod->subdum);
	if (pmod->coeff) free(pmod->coeff);
	if (pmod->sderr) free(pmod->sderr);
	if (pmod->yhat) free(pmod->yhat);
	if (pmod->uhat) free(pmod->uhat);
	if (pmod->xpx) free(pmod->xpx);
	if (pmod->vcv) free(pmod->vcv);
	if (pmod->name) free(pmod->name);
	if (pmod->arinfo) {
	    clear_ar_info(pmod);
	}
	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    if (pmod->slope) free(pmod->slope);
	}
	if (pmod->ntests) free(pmod->tests);
	if (pmod->params) {
	    int i;

	    for (i=0; i<=pmod->ncoeff; i++) {
		free(pmod->params[i]);
	    }
	    free(pmod->params);
	}
	if (pmod->data) {
	    MISSOBS *mobs = (MISSOBS *) pmod->data;
	    free(mobs->missvec);
	    free(pmod->data);
	}
	if (pmod->dataset) {
	    free_model_dataset(pmod);
	}
    }

    _init_model(pmod, pdinfo);
}

/* .......................................................... */

static char *internal_path_stuff (int code, const char *path)
{
    static char gretl_lib_path[MAXLEN];

    if (code == 1) {
#ifdef WIN32
	strcpy(gretl_lib_path, path);
#else
# ifdef GLIB2 
	const char *sfx = "-gtk2/";
# else
	const char *sfx = "-gtk1/";
# endif
	char *p = strstr(path, "/share");

	if (p) {
	    size_t len = p - path;

	    *gretl_lib_path = 0;
	    strncat(gretl_lib_path, path, len);
	    strcat(gretl_lib_path, "/lib/gretl");
	    strcat(gretl_lib_path, sfx);
	} else {
	    sprintf(gretl_lib_path, "%s/lib/gretl%s", path, sfx);
	}
#endif
	return NULL;
    } 
    else if (code == 0) {
	return gretl_lib_path;
    }
    return NULL;
}

const char *fetch_gretl_lib_path (void)
{
    return internal_path_stuff (0, NULL);
}
    
/* .......................................................... */

void show_paths (PATHS *ppaths)
{
    printf(_("gretl: using these basic search paths:\n"));
    printf("gretldir: %s\n", ppaths->gretldir);
    printf("userdir: %s\n", ppaths->userdir);
    printf("datadir: %s\n", ppaths->datadir);
    printf("scriptdir: %s\n", ppaths->scriptdir);
    printf("gnuplot: %s\n", ppaths->gnuplot);
}

/* .......................................................... */

#ifdef WIN32

int set_paths (PATHS *ppaths, int defaults, int gui)
{
    char envstr[MAXLEN];

    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL)
	    strcpy(ppaths->gretldir, home);
	else
	    strcpy(ppaths->gretldir, "c:\\userdata\\gretl");

	sprintf(ppaths->binbase, "%s\\db\\", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "f:\\"); 

	if (gui) {
	    strcpy(ppaths->dbhost_ip, "152.17.150.2");
	} else {
	    ppaths->dbhost_ip[0] = '\0';
	}
	ppaths->currdir[0] = '\0';

	strcpy(ppaths->pngfont, "verdana 8");
    }

    sprintf(ppaths->datadir, "%s\\data\\", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%s\\scripts\\", ppaths->gretldir);

    if (gui) {
	sprintf(ppaths->helpfile, "%s\\%s", ppaths->gretldir,
		_("gretl_hlp.txt"));
	sprintf(ppaths->cmd_helpfile, "%s\\%s", ppaths->gretldir,
		_("gretlcli_hlp.txt"));
    } else { 
	sprintf(ppaths->helpfile, "%s\\%s", ppaths->gretldir,
		_("gretlcli_hlp.txt"));
    }

    if (ppaths->userdir[strlen(ppaths->userdir) - 1] != SLASH)
	strcat(ppaths->userdir, "\\");

    *ppaths->plotfile = '\0';

    sprintf(envstr, "GTKSOURCEVIEW_LANGUAGE_DIR=%s\\share\\gtksourceview-1.0"
	    "\\language-specs", ppaths->gretldir);
    putenv(envstr);

    internal_path_stuff (1, ppaths->gretldir);

    return 0;
}

#else

int set_paths (PATHS *ppaths, int defaults, int gui)
{
    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL) {
	    strcpy(ppaths->gretldir, home);
	} else {
	    strcpy(ppaths->gretldir, GRETL_PREFIX);
	    strcat(ppaths->gretldir, "/share/gretl/");
	} 

	sprintf(ppaths->binbase, "%sdb/", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "/mnt/dosc/userdata/rats/oecd/");

	if (gui) {
	    strcpy(ppaths->dbhost_ip, "152.17.150.2");
	} else {
	    ppaths->dbhost_ip[0] = '\0';
	}

	strcpy(ppaths->gnuplot, "gnuplot");
	*ppaths->pngfont = 0;

	ppaths->currdir[0] = '\0';	

	/* try to set a default userdir */
	home = getenv("HOME");
	if (home != NULL) {
	    strcpy(ppaths->userdir, home);
	    strcat(ppaths->userdir, "/gretl/");
	} else {
	    strcpy(ppaths->userdir, "");
	}
    } 

    sprintf(ppaths->datadir, "%sdata/", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts/", ppaths->gretldir);
    
    if (gui) {
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretl.hlp"));
	sprintf(ppaths->cmd_helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli.hlp"));
    } else
	sprintf(ppaths->helpfile, "%s%s", ppaths->gretldir,
		_("gretlcli.hlp"));

    *ppaths->plotfile = '\0';

    internal_path_stuff (1, ppaths->gretldir);

    return 0;
}

#endif /* win32 versus unix */

/* ........................................................... */

int copylist (int **target, const int *src)
{
    int i, n;

    if (src == NULL) return 1;
    n = src[0];

    if (*target != NULL) free(*target);
    *target = malloc((n + 2) * sizeof *target);
    if (*target == NULL) return 1;

    for (i=0; i<=n; i++) (*target)[i] = src[i];

    return 0;
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
    for (t=0; t<pdinfo->n; t++) (*pZ)[0][t] = 1.0;
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
	pdinfo->varname[v+i] = malloc(9);
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
    pdinfo->varname[v] = malloc(9);
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

static int varnum_from_string (const char *str, DATAINFO *pdinfo)
{
    int varno;
    char *test;
    extern int errno;

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

int 
dataset_drop_var_wrapper (const char *str, double ***pZ, DATAINFO *pdinfo)
{
    int varno = varnum_from_string(str, pdinfo);

    if (varno < 0) return E_DATA;
    else return dataset_drop_var(varno, pZ, pdinfo);
}

/* ......................................................  */

int dataset_drop_var (int varno, double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char *vector;
    VARINFO **varinfo;
    int i, v = pdinfo->v; 

    if (varno <= 0 || varno >= v) {
	return E_DATA;
    }

    free(pdinfo->varname[varno]);
    if (pdinfo->varinfo[varno] != NULL) {
	free(pdinfo->varinfo[varno]);
    }
    free((*pZ)[varno]);

    for (i=varno; i<v-1; i++) {
	pdinfo->varname[i] = pdinfo->varname[i+1];
	pdinfo->varinfo[i] = pdinfo->varinfo[i+1];
	(*pZ)[i] = (*pZ)[i+1];
    }

    varname = realloc(pdinfo->varname, (v-1) * sizeof *varname);
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, (v-1) * sizeof *vector);
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;

    varinfo = realloc(pdinfo->varinfo, (v-1) * sizeof *varinfo);
    if (varinfo == NULL) return E_ALLOC;
    else pdinfo->varinfo = varinfo;

    newZ = realloc(*pZ, (v-1) * sizeof *newZ); 
    if (newZ == NULL) return E_ALLOC;
    else *pZ = newZ;

    pdinfo->v -= 1;

    return 0;
}

/* ......................................................  */

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
	strcmp(pdinfo->varname[i], "hrs") == 0)
	return 1;
    return 0;
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

static ARINFO *copy_ar_info (const ARINFO *src)
{
    ARINFO *targ;

    targ = malloc(sizeof *targ);
    if (targ == NULL) return NULL;

    if (src->arlist != NULL) {
	int i, m = src->arlist[0];

      	if ((targ->rho = copyvec(src->rho, m)) == NULL) {
	    free(targ);
	    targ = NULL;
	    return NULL; 
	}

      	if ((targ->sderr = copyvec(src->sderr, m)) == NULL) { 
	    free(targ->rho);
	    free(targ);
	    targ = NULL;
	    return NULL; 
	}

	targ->arlist = malloc((m + 1) * sizeof(int));
	if (targ->arlist == NULL) {
	    free(targ->rho);
	    free(targ->sderr);
	    free(targ);
	    targ = NULL;
	    return NULL; 
	}
	for (i=0; i<=m; i++) {
	    targ->arlist[i] = src->arlist[i];
	}
    }    

    return targ;
}

/* ........................................................... */

int copy_model (MODEL *targ, const MODEL *src, const DATAINFO *pdinfo)
{
    int i = src->list[0] - 1;
    int m = i * (i + 1) / 2;

    /* monolithic copy of structure */
    *targ = *src;

    /* now work on pointer members */
    _init_model(targ, pdinfo);
    if ((targ->coeff = copyvec(src->coeff, src->ncoeff)) == NULL)
	return 1;
    if ((targ->sderr = copyvec(src->sderr, src->ncoeff))  == NULL)  
	return 1;
    if ((targ->uhat = copyvec(src->uhat, pdinfo->n)) == NULL) return 1;
    if ((targ->yhat = copyvec(src->yhat, pdinfo->n)) == NULL) return 1;
    if (src->subdum != NULL && 
	(targ->subdum = copyvec(src->subdum, pdinfo->n)) == NULL) return 1;

    if (src->xpx != NULL &&
	(targ->xpx = copyvec(src->xpx, m)) == NULL) return 1;
    if (src->vcv != NULL && 
	(targ->vcv = copyvec(src->vcv, m)) == NULL) return 1;

    if (src->arinfo != NULL) {
	targ->arinfo = copy_ar_info(src->arinfo);
	if (targ->arinfo == NULL) 
	    return 1; 
    }

    if (src->slope != NULL &&
	(targ->slope = copyvec(src->slope, src->ncoeff)) == NULL)
	    return 1;

    m = src->list[0];
    targ->list = malloc((m + 1) * sizeof *targ->list);
    if (targ->list == NULL) return 1;
    for (i=0; i<=m; i++) targ->list[i] = src->list[i];    

    return 0;
}

/* ........................................................... */

int swap_models (MODEL **targ, MODEL **src)
{
    MODEL *tmp = *targ;

    *targ = *src;
    *src = tmp;
    return 0;
}

/* ........................................................... */

int _forecast (int t1, int t2, int nv, 
	       const MODEL *pmod, double ***pZ)
{
    double xx, zz, zr;
    int i, k, maxlag = 0, yno = pmod->list[1], ARMODEL;
    int v, t, miss;
    const int *arlist = NULL;

    ARMODEL = (pmod->ci == AR || pmod->ci == CORC || 
	       pmod->ci == HILU)? 1 : 0;

    if (ARMODEL) {
	arlist = pmod->arinfo->arlist;
	maxlag = arlist[arlist[0]];
	if (t1 < maxlag) t1 = maxlag; 
    }

    for (t=t1; t<=t2; t++) {
	miss = 0;
	zz = 0.0;
	if (ARMODEL) 
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
	} /* end if ARMODEL */
	for (v=0; !miss && v<pmod->ncoeff; v++) {
	    k = pmod->list[v+2];
	    xx = (*pZ)[k][t];
	    if (na(xx)) {
		zz = NADBL;
		miss = 1;
	    }
	    if (!miss && ARMODEL) {
		xx = (*pZ)[k][t];
		for (i=1; i<=arlist[0]; i++) {
		    xx -= pmod->arinfo->rho[i] * (*pZ)[k][t - arlist[i]];
		}
	    }
	    if (!miss) {
		zz = zz + xx * pmod->coeff[v];
	    }
	}
	(*pZ)[nv][t] = zz;
    }
    return 0;
}

/* ........................................................... */

int _full_model_list (MODEL *pmod, int **plist)
/* reconstitute full varlist for WLS and AR models */
{
    int i, pos = 0, len, *mylist = NULL;

    if (pmod->ci != WLS && pmod->ci != AR) return 0;

    if (pmod->ci == WLS) { 
	mylist = malloc(((*plist)[0] + 2) * sizeof *mylist);
	if (mylist == NULL) return -1;
	for (i=1; i<=(*plist)[0]; i++) 
	    mylist[i+1] = (*plist)[i];
	mylist[0] = (*plist)[0] + 1;
	mylist[1] = pmod->nwt;
    }
    else if (pmod->ci == AR) {
	pos = pmod->arinfo->arlist[0] + 1;
	len = pos + (*plist)[0] + 2;
	mylist = malloc(len * sizeof *mylist);
	if (mylist == NULL) return -1;
	mylist[0] = len - 2;
	for (i=1; i<pos; i++) mylist[i] = pmod->arinfo->arlist[i];
	mylist[pos] = 999;
	for (i=1; i<=(*plist)[0]; i++) {
	    mylist[pos+i] = (*plist)[i];
	}
    }

    copylist(plist, mylist);
    free(mylist);

    return pos;
}

/* ........................................................... */

FITRESID *get_fit_resid (const MODEL *pmod, double ***pZ, 
			 DATAINFO *pdinfo)
{
    int depvar, t, nfit = 0;
    int t1 = pmod->t1, t2 = pmod->t2, n = pdinfo->n;
    char fcastline[32];
    FITRESID *fr;

    depvar = pmod->list[1];

    if (pmod->data != NULL) {
	t2 += get_misscount(pmod);
    }

    if (pmod->ci != NLS) {
	sprintf(fcastline, "fcast %s %s fitted", pdinfo->stobs, 
		pdinfo->endobs);
	nfit = fcast(fcastline, pmod, pdinfo, pZ); 
	if (nfit < 0) return NULL; 
    }

    fr = fit_resid_new(n, 0);
    if (fr == NULL) return NULL;

    fr->sigma = pmod->sigma;
    for (t=0; t<n; t++) {
	fr->actual[t] = (*pZ)[depvar][t];
	if (pmod->ci == NLS) fr->fitted[t] = pmod->yhat[t];
	else fr->fitted[t] = (*pZ)[nfit][t];
    }

    if (isdummy(fr->actual, 0, n) > 0) {
	fr->pmax = get_precision(fr->fitted, n);
    } else {
	fr->pmax = get_precision(fr->actual, n);
    }

    strcpy(fr->depvar, pdinfo->varname[depvar]);
    
    fr->t1 = t1;
    fr->t2 = t2;
    fr->nobs = pmod->nobs;

    /* should we delete the fitted value from *pZ? */

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
    char t1str[9], t2str[9];

    fr = fit_resid_new(0, 1); 
    if (fr == NULL) return NULL;

    if (pmod->ci != OLS) {
	fr->err = E_OLSONLY;
	return fr;
    }

    /* bodge (rejected in case of subsampled data) */
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
    finfo->extra = 1;

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

    _init_model(&fmod, finfo);

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
	clear_model(&fmod, finfo);
	fmod = lsq(list, &fZ, finfo, OLS, 1, 0.0);
	if (fmod.errcode) {
	    fr->err = fmod.errcode;
	    clear_model(&fmod, finfo);
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

    clear_model(&fmod, finfo);

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

    fr = get_fcast_with_errs (str, pmod, pZ, pdinfo, prn);
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

int is_model_cmd (const char *line)
{
    if (!strncmp(line, "ols", 3)  ||
	!strncmp(line, "corc", 4) ||
	!strncmp(line, "hilu", 4) ||
	!strncmp(line, "wls", 3)  ||
	!strncmp(line, "pooled", 6)  ||
	!strncmp(line, "hccm", 4) ||
	!strncmp(line, "hsk", 3)  ||
	!strncmp(line, "add", 3)  ||
	!strncmp(line, "lad", 3)  ||
	!strncmp(line, "omit", 4) ||
	!strncmp(line, "tsls", 4) ||
	!strncmp(line, "logit", 5)  ||
	!strncmp(line, "probit", 6) ||
	!strncmp(line, "end nls", 7) ||
	!strncmp(line, "ar", 2))
	return 1;
    return 0;
}

/* ........................................................... */

int is_model_ref_cmd (int ci)
{
    if (ci == ADD ||
	ci == OMIT ||
	ci == ARCH ||
	ci == CHOW ||
	ci == CUSUM ||
	ci == LMTEST ||
	ci == FCAST ||
	ci == FCASTERR ||
	ci == FIT)
	return 1;
    return 0;
}

/* ........................................................... */

static void store_list (int *list, char *buf)
{
    int i;
    char numstr[5];

    for (i=1; i<=list[0]; i++) {
	sprintf(numstr, "%d ", list[i]);
	strcat(buf, numstr);
    }
}

/* ........................................................... */

int save_model_spec (MODEL *pmod, MODELSPEC *spec, DATAINFO *fullinfo)
{
    if (pmod->list == NULL) return 1;

    sprintf(spec->cmd, "%s ", gretl_commands[pmod->ci]);
    
    if (pmod->ci == AR) {
	store_list(pmod->arinfo->arlist, spec->cmd);
	strcat(spec->cmd, "; ");
    }

    store_list(pmod->list, spec->cmd);

    if (pmod->subdum != NULL) {
	int t;

	spec->subdum = malloc(fullinfo->n * sizeof *spec->subdum);
	if (spec->subdum == NULL) return 1;
	for (t=0; t<fullinfo->n; t++)
	    spec->subdum[t] = pmod->subdum[t];
    }

    return 0;
}

/* ........................................................... */

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo) 
{
    CMD command;
    int err = 0, ignore = 0, model_count = 0;
    double rho = 0;
    PRN prn;

    prn.fp = NULL;
    prn.buf = NULL;

    command.list = malloc(sizeof *command.list);
    command.param = malloc(1);
    if (command.list == NULL || command.param == NULL) 
	return 1;

    getcmd(model_spec, pdinfo, &command, &ignore, pZ, NULL);

    _init_model(tmpmod, pdinfo);

    switch(command.ci) {
    case AR:
	*tmpmod = ar_func(command.list, atoi(command.param), pZ, 
			  pdinfo, &model_count, &prn);
	break;
    case CORC:
    case HILU:
	err = hilu_corc(&rho, command.list, pZ, pdinfo, 
			NULL, 1, command.ci, &prn);
	if (!err) {
	    *tmpmod = lsq(command.list, pZ, pdinfo, command.ci, 0, rho);
	}
	break;
    case HCCM:
	*tmpmod = hccm_func(command.list, pZ, pdinfo);
	break;
    case HSK:
	*tmpmod = hsk_func(command.list, pZ, pdinfo);
	break;
    case LOGIT:
    case PROBIT:
	*tmpmod = logit_probit(command.list, pZ, pdinfo, command.ci);
	break;
    case OLS:
    case WLS:
    case POOLED:
	*tmpmod = lsq(command.list, pZ, pdinfo, command.ci, 0, 0.0);
	break;
    case TSLS:
	break;
    default:
	break;
    }

    if (tmpmod->errcode) {
	err = 1;
	clear_model(tmpmod, pdinfo);
    }
    if (command.list) free(command.list);
    if (command.param) free(command.param);

    return err;
}

/* ........................................................... */

int guess_panel_structure (double **Z, DATAINFO *pdinfo)
{
    int v, panel;

    v = varindex(pdinfo, "year");
    if (v == pdinfo->v)
	v = varindex(pdinfo, "Year");
    if (v == pdinfo->v)
	panel = 0; /* can't guess */
    else {
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
    char unit[9], period[9];

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
    int i, ncoeff = pmod->list[0];
    double t = tcrit95(pmod->dfd);
    CONFINT *cf;

    cf = malloc(sizeof *cf);
    if (cf == NULL) return NULL;

    cf->coeff = malloc(ncoeff * sizeof *cf->coeff);
    if (cf->coeff == NULL) {
	free(cf);
	return NULL;
    }

    cf->maxerr = malloc(ncoeff * sizeof *cf->maxerr);
    if (cf->maxerr == NULL) {
	free(cf);
	free(cf->coeff);
	return NULL;
    }

    cf->list = NULL;
    if (copylist(&cf->list, pmod->list)) {
	free(cf);
	free(cf->coeff);
	free(cf->maxerr);
	return NULL;
    }

    for (i=0; i<ncoeff; i++) { /* FIXME? */
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

#ifdef GLIB2

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
	sprintf(gretl_errmsg, "%s\n%s", 
		_("Command failed"),
		sout);
	fprintf(stderr, "status=%d: '%s'\n", status, sout);
	ret = 1;
    }

    if (errout != NULL) g_free(errout);
    if (sout != NULL) g_free(sout);

    if (ret) {
	fprintf(stderr, "Failed command: '%s'\n", cmdline);
    } 

    return ret;
}

#elif !defined(WIN32)

#include <signal.h>

int gretl_spawn (const char *cmdline)
{
    int err;
    extern int errno;

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
