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

#ifdef OS_WIN32
# include "../../winconfig.h"
#else
# include "../../config.h"
#endif
#include "libgretl.h"
#include "internal.h"
#include <dirent.h>
#include <unistd.h>

static DIR *dir;

static int _pdton (const int pd);

extern int plot_fcast_errs (const int n, const double *obs, 
			    const double *depvar, const double *yhat, 
			    const double *maxerr, const char *varname, 
			    const PATHS *ppaths);

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

double _corr (const int n, const double *zx, const double *zy)
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
    for (i = 0; i < n; ++i) {
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

double _covar (const int n, const double *zx, const double *zy)
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

double date (const int nt, const int pd, const double sd0)
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
 * @i: row number.
 * @j: column number.
 * @lo: total number of elements in array.
 *
 * Given references i (row) and j (column) of a 2-dimensional array,
 * finds the corresponding position in the 1-dimensional array of
 * the same elements.
 *
 * Returns: position in 1-dimensional array.
 */

int ijton (const int i, const int j, const int lo)
{
    int n;

    n = lo * (i - 1) + j - i - ((i - 2) * (i - 1)/2);
    return n;
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

int ztox (const int i, double *px, double **Z, const DATAINFO *pdinfo) 
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
 * @varnum: index number of variable to examine.
 * @t1: starting observation.
 * @t2: ending observation. 
 * @Z: data matrix.
 * @n: full length of data series in @Z.
 * 
 * Check whether variable @varnum has only 0 or 1 values over the
 * given sample range. 
 *
 * Returns: 0 if the variable is not a 0/1 dummy, otherwise the
 * number of 1s in the series.
 */

int isdummy (const int varnum, const int t1, const int t2, 
	     double **Z, const int n)
{
    int t, m = 0;
    double xx;

    for (t=t1; t<=t2; t++) {
	xx = Z[varnum][t];
	if (floatneq(xx, 0.0) && floatneq(xx, 1.0)) return 0;
	if (floateq(xx, 1.0)) m++;
    }
    if (m < t2 - t1 + 1) return m;
    return 0;
} 

/* ........................................................  */

int _iszero (const int t1, const int t2, const double *x)
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

void list_exclude (const int n, int *list)
{
    int i;

    for (i=n; i<list[0]; i++) list[i] = list[i+1];
    list[0] = list[0] - 1;
}

/* ........................................................  */

int _isconst (const int t1, const int t2, const double *x)
{
    int t;
    double xx = x[t1];

    for (t=t1+1; t<=t2; t++) if (floatneq(x[t], xx)) return 0;
    return 1;
}

/* ............................................................  */

double _esl_mean (const int t1, const int t2, const double *x)
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
    xbar = xbar + sum/n;
    return xbar;
}

/* ......................................................  */

void _minmax (const int t1, const int t2, const double zx[], 
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

static int _pdton (const int pd)
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
        if (list[i] == 0) return 1;

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

double _esl_stddev (const int t1, const int t2, const double *x)
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

double _esl_variance (const int t1, const int t2, const double *x)
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

void _aicetc (MODEL *pmod)
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
    pmod->criterion[2] = ersq * (nobs+ncoeff)/zz;
    zz = 2.0 * ncoeff/nobs;
    pmod->criterion[1] = ersq * exp(zz);
    pmod->criterion[5] = ersq *(1.0+zz);
    pmod->criterion[7] = (1-zz)>0.0? ersq/(1-zz): NADBL;
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

void _criteria (const double ess, const int nobs, const int ncoeff, 
		PRN *prn)
{
    double zz, zx, ersq, zn;
    double criterion[8];
    
    zz = (double) (nobs - ncoeff);
    criterion[0] = ess/zz;
    ersq = ess/nobs;
    criterion[2] = ersq * (nobs+ncoeff)/zz;
    zz = 2.0 * ncoeff/nobs;
    criterion[1] = ersq * exp(zz);
    criterion[5] = ersq *(1.0+zz);
    criterion[7] = (1-zz)>0.0? ersq/(1-zz): NADBL;
    zn = (double) nobs;
    zx = log(zn);
    criterion[3] = ersq * pow(zx, zz);
    zz = (double) ncoeff;
    zz = zz/zn;
    criterion[4] = ersq * pow(zn, zz);
    zz = 1.0 - zz;
    criterion[6] = ersq/(zz*zz);

    pprintf(prn, "Using ess = %f, %d observations, %d coefficients\n", 
	   ess, nobs, ncoeff);
    pprintf(prn, "\nMODEL SELECTION STATISTICS\n\n");	
    pprintf(prn, "SGMASQ    %13g     AIC       %13g     FPE       %12g\n"
	    "HQ        %13g     SCHWARZ   %13g     SHIBATA   %12g\n"
	    "GCV       %13g",
	    criterion[0], criterion[1], 
	    criterion[2], criterion[3], 
	    criterion[4], criterion[5], criterion[6]);
    if (criterion[7] > 0.0) pprintf(prn, "     RICE      %13g\n", 
					  criterion[7]);
    else pprintf(prn, "     RICE          undefined\n");
    pprintf(prn, "\n");
}

/* ....................................................... */

int _adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		  double **Z, const int n, int *misst)
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

int set_obs (char *line, DATAINFO *pdinfo, int opt)
{
    int pd, pos, i, len, dc = 0, bad = 0;
    char stobs[8], endobs[8], endbit[7];

    gretl_errmsg[0] = '\0';

    if (sscanf(line, "%*s %d %7s", &pd, stobs) != 2) {
	strcpy(gretl_errmsg, "Failed to parse line as frequency, startobs");
	return 1;
    }

    /* does frequency make sense? */
    if (pd < 1 || pd > pdinfo->n) {
	sprintf(gretl_errmsg, 
		"frequency (%d) does not make seem to make sense", pd);
	return 1;
    }
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
	sprintf(gretl_errmsg, "starting obs '%s' is invalid", stobs);
	return 1;
    }
    pos = dotpos(stobs);
    if (pd > 1 && pos == len) {
	strcpy(gretl_errmsg, "starting obs must contain a '.' with frequency > 1");
	return 1;
    }
    if (pd == 1 && pos < len) {
	strcpy(gretl_errmsg, "no '.' allowed in starting obs with frequency 1");
	return 1;
    }    
    if ((pd > 1 && pd < 10 && strlen(stobs+pos) != 2) ||
	(pd >= 10 && pd < 100 && strlen(stobs+pos) != 3)) {
	sprintf(gretl_errmsg, "starting obs '%s' is incompatible with frequency", 
		stobs);
	return 1;
    }
    if (pd > 1) {
	strcpy(endbit, stobs+pos+1);
	dc = atoi(endbit);
	if (dc < 0 || dc > pd) {
	    sprintf(gretl_errmsg, 
		    "starting obs '%s' is incompatible with frequency", 
		    stobs);
	    return 1;
	}	    
    }

    /* adjust data info struct */
    pdinfo->pd = pd;
    strcpy(pdinfo->stobs, stobs);
    pdinfo->sd0 = atof(stobs);
    ntodate(endobs, pdinfo->n - 1, pdinfo);
    strcpy(pdinfo->endobs, endobs);

    if (opt == OPT_S) pdinfo->time_series = STACKED_TIME_SERIES;
    else if (opt == OPT_C) pdinfo->time_series = STACKED_CROSS_SECTION;
    else if (pdinfo->sd0 >= 2.0) 
        pdinfo->time_series = TIME_SERIES; /* actual time series? */
    else if (pdinfo->sd0 > 1.0) 
	pdinfo->time_series = STACKED_TIME_SERIES; /* panel data? */
    else pdinfo->time_series = 0;

    /* and report */
    fprintf(stderr, "setting data frequency = %d\n", pd);
    fprintf(stderr, "data range: %s - %s\n", stobs, endobs);
    return 0;
}

#ifdef OS_WIN32
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
#ifdef OS_WIN32
    char *tmp = unslash(topdir);

    if (tmp == NULL) return -1;
#endif

    if (first) {
#ifdef OS_WIN32
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
			 const int recurse)
{
    FILE *test;
    int got = 0;
    char origfile[MAXLEN], trypath[MAXLEN];

    strcpy(origfile, filename);

    if (path_append(filename, topdir) == 0) {
	fprintf(stderr, "Trying %s\n", filename);
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
		    fprintf(stderr, "Trying %s\n", filename);
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
    if (test != NULL) { /* found it */
	fclose(test); 
	/* if a relative path was given, convert it to absolute */
#ifdef OS_WIN32
	if (fname[1] == ':') return fname;
#endif
	if (fname[0] != SLASH && (thisdir = malloc(MAXLEN)) != NULL) {
	    int i = 0;

	    if (getcwd(thisdir, MAXLEN-1) != NULL) {
#ifdef OS_WIN32		
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
    } else {  /* not able to open file as given */
	if (fname[0] == '.' || fname[0] == SLASH)
	    return NULL;
    }

    /* try looking where script was found */
    if (ppaths->currdir[0]) {
	if ((fname = search_dir(fname, ppaths->currdir, 0))) 
	    return fname;
    }

    /* try looking in user's dir */
    fname = tmp;
    strcpy(fname, orig);
    if ((fname = search_dir(fname, ppaths->userdir, 1))) 
	return fname;

    /* if it's a data file we want, try system data dir */
    fname = tmp;
    strcpy(fname, orig);
    if (!script) {
	if ((fname = search_dir(fname, ppaths->datadir, 1))) 
	    return fname;
    } else {
	/* for a script, try system script dir */
	if ((fname = search_dir(fname, ppaths->scriptdir, 1))) 
	    return fname;
    }
    fname = tmp;
    strcpy(fname, orig);
    return NULL;
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

static int isflag (int c)
{
    switch (c) {
    case 'o': return OPT_O;
    case 'c': return OPT_C;
    case 'm': return OPT_M;
    case 'r': return OPT_R;
    case 's': return OPT_S;
    case 't': return OPT_T;
    case 'l': return OPT_L;
    case 'z': return OPT_Z;
    default: return 0;
    }
}

/* ........................................................... */

char getflag (int opt)
{
    switch (opt) {
    case OPT_R: return 'r';
    case OPT_S: return 's';
    case OPT_T: return 't';
    case OPT_O: return 'o';
    case OPT_C: return 'c';
    case OPT_M: return 'm';
    case OPT_L: return 'l';
    case OPT_Z: return 'z';
    default: return 0;
    }
}

/* .......................................................... */

int catchflag (char *line, int *oflag)
     /* check for "-<char>" in line: if found, chop it out and set
	oflag value accordingly.  
	Strip trailing semicolon while we're at it. 
     */
{
    int i, opt, n = strlen(line);

    *oflag = 0;

    /* to enable reading of trad. esl input files */
    if (line[n-2] == ';') {
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

    if (ci == TSLS || ci == AR || ci == SCATTERS) {
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
    pmod->arlist = NULL;
    pmod->rhot = NULL;
    pmod->slope = NULL;
    pmod->infomsg[0] = '\0';
    pmod->name = NULL;
    pmod->ntests = 0;
    pmod->tests = NULL;
    pmod->errcode = 0;
    pmod->aux = 0;
    gretl_errmsg[0] = '\0';
}

/* ........................................................... */

MODEL *gretl_model_new (DATAINFO *pdinfo)
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

/* #define SESSION_DEBUG */

/* ........................................................... */

int silent_remember (MODEL **ppmod, SESSION *psession, SESSIONBUILD *rebuild,
		     DATAINFO *pdinfo)
{
    int i = psession->nmodels;
    MODEL *pmod = *ppmod;
    MODEL *tmp;

#ifdef SESSION_DEBUG
    fprintf(stderr, "psession->nmodels = %d\n", psession->nmodels);
#endif

    if ((pmod->name = malloc(64)) == NULL) return 1;
    strncpy(pmod->name, rebuild->model_name[i], 63);

    if (psession->nmodels)
	psession->models = realloc(psession->models, 
				   (i + 1) * sizeof(MODEL *));
    else
	psession->models = malloc(sizeof(MODEL *));
    if (psession->models == NULL) return 1;
    psession->nmodels += 1;
    psession->models[i] = pmod;
    tmp = malloc(sizeof *tmp);
    if (tmp == NULL) return 1;
    *ppmod = tmp;
    _init_model(tmp, pdinfo);

#ifdef SESSION_DEBUG
    fprintf(stderr, "copied '%s' to psession->models[%d]\n" 
	    " nmodels = %d\n", rebuild->model_name[i], i, psession->nmodels); 
#endif
    return 0;
}
    
/* .......................................................... */

int clear_model (void *ptr, SESSION *psession, SESSIONBUILD *rebuild,
		 DATAINFO *pdinfo)
{
    int i;
    static int save;
    MODEL *pmod;
    MODEL **ppmod;

    if (rebuild && psession) {
	ppmod = (MODEL **) ptr;
	pmod = *ppmod;

#ifdef SESSION_DEBUG
	fprintf(stderr, "clear_model: rebuild & psession non-NULL, save = %d\n", 
		save);
#endif
	if (save) {
	    for (i=0; i<rebuild->nmodels; i++) {
#ifdef SESSION_DEBUG
		fprintf(stderr, "i=%d, pmod->ID=%d, rebuild->model_ID[%d] = %d\n",
			i, pmod->ID, i, rebuild->model_ID[i]);
#endif		
		if (pmod->ID == rebuild->model_ID[i]) {
#ifdef SESSION_DEBUG
  		    fprintf(stderr, "Rebuilding saved model %d (%s)\n",  
  			   pmod->ID, rebuild->model_name[i]);
#endif	 
		    return silent_remember(ppmod, psession, rebuild, pdinfo);
		}
	    }
	}
	save = 1;
    } else 
	pmod = (MODEL *) ptr;

    if (pmod != NULL) {
	if (pmod->list) free(pmod->list);
	if (pmod->subdum) free(pmod->subdum);
	if (pmod->coeff) free(pmod->coeff);
	if (pmod->sderr) free(pmod->sderr);
	if (pmod->yhat) free(pmod->yhat);
	if (pmod->uhat) free(pmod->uhat);
	if (pmod->xpx) free(pmod->xpx);
	if (pmod->vcv) free(pmod->vcv);
	if (pmod->name) free(pmod->name);
	if (pmod->ci == AR || pmod->ci == CORC || pmod->ci == HILU) {
	    if (pmod->arlist) free(pmod->arlist);
	    if (pmod->rhot) free(pmod->rhot);
	}
	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    if (pmod->slope) free(pmod->slope);
	}
	if (pmod->ntests) free(pmod->tests);
    }
    _init_model(pmod, pdinfo);

    return 0;
}

/* .......................................................... */

void show_paths (PATHS *ppaths)
{
    printf("gretl: using these basic search paths:\n");
    printf("gretldir: %s\n", ppaths->gretldir);
    printf("userdir: %s\n", ppaths->userdir);
    printf("datadir: %s\n", ppaths->datadir);
    printf("scriptdir: %s\n", ppaths->scriptdir);
    printf("gnuplot: %s\n", ppaths->gnuplot);
}

/* .......................................................... */

#ifdef OS_WIN32

int set_paths (PATHS *ppaths, const int defaults, const int gui)
{
    if (defaults) {
	char *home;

	home = getenv("GRETL_HOME");
	if (home != NULL)
	    strcpy(ppaths->gretldir, home);
	else
	    strcpy(ppaths->gretldir, "c:\\userdata\\gretl"); 	
	sprintf(ppaths->binbase, "%s\\db\\", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "f:\\"); 
	strcpy(ppaths->dbhost_ip, "152.17.150.2");
	ppaths->hdrfile[0] = '\0';
	ppaths->currdir[0] = '\0';
    }

    sprintf(ppaths->datadir, "%s\\data\\", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%s\\scripts\\", ppaths->gretldir);

    if (gui) {
	sprintf(ppaths->helpfile, "%s\\gretl_hlp.txt", ppaths->gretldir);
	sprintf(ppaths->cmd_helpfile, "%s\\gretlcli_hlp.txt", ppaths->gretldir);
    } else 
	sprintf(ppaths->helpfile, "%s\\gretlcli_hlp.txt", ppaths->gretldir);

    if (ppaths->userdir[strlen(ppaths->userdir) - 1] != SLASH)
	strcat(ppaths->userdir, "\\");

    sprintf(ppaths->plotfile, "%sgpttmp.plt", ppaths->userdir);

    if (get_base(ppaths->pgnuplot, ppaths->gnuplot, SLASH))
	strcat(ppaths->pgnuplot, "pgnuplot.exe");
    else
	strcpy(ppaths->pgnuplot, "pgnuplot.exe");	

    return 0;
}

#else

int set_paths (PATHS *ppaths, const int defaults, const int gui)
{
    if (defaults) {
	char *home;
	DIR *try = NULL;

	home = getenv("GRETL_HOME");
	if (home != NULL)
	    strcpy(ppaths->gretldir, home);
	else
	    strcpy(ppaths->gretldir, GRETL); 
	sprintf(ppaths->binbase, "%sdb/", ppaths->gretldir);
	strcpy(ppaths->ratsbase, "/mnt/dosc/userdata/rats/oecd/");
	strcpy(ppaths->dbhost_ip, "152.17.150.2");
	strcpy(ppaths->gnuplot, "gnuplot");
	ppaths->hdrfile[0] = '\0';
	ppaths->currdir[0] = '\0';

	/* figure out user's home gretl directory */
	home = getenv("HOME");
	if (home != NULL) {
	    strcpy(ppaths->userdir, home);
	    strcat(ppaths->userdir, "/gretl/");
	    if ((try = opendir(ppaths->userdir)) == NULL) {
		/* see if there's an old-style user dir */
		strcpy(ppaths->userdir, home);
		strcat(ppaths->userdir, "/esl/");
		if ((try = opendir(ppaths->userdir)) == NULL) {
		    strcpy(ppaths->userdir, home);
		    strcat(ppaths->userdir, "/gretl/");
		} 
	    } 
	} else 
	    strcpy(ppaths->userdir, "");
	if (try != NULL) closedir(try);
    } 

    sprintf(ppaths->datadir, "%sdata/", ppaths->gretldir);
    sprintf(ppaths->scriptdir, "%sscripts/", ppaths->gretldir);
    
    if (gui) {
	sprintf(ppaths->helpfile, "%sgretl.hlp", ppaths->gretldir);
	sprintf(ppaths->cmd_helpfile, "%sgretlcli.hlp", ppaths->gretldir);
    } else
	sprintf(ppaths->helpfile, "%sgretlcli.hlp", ppaths->gretldir);

    sprintf(ppaths->plotfile, "%sgpttmp.plt", ppaths->userdir);

    return 0;
}

#endif

/* ........................................................... */

int copylist (int **target, const int *src)
{
    int i, n;

    if (src == NULL) return 1;
    n = src[0];
    if (*target != NULL) free(*target);
    *target = malloc((n + 2) * sizeof(int));
    if (*target == NULL) return 1;
    for (i=0; i<=n; i++) (*target)[i] = src[i];
    return 0;
}

/* ......................................................  */

int grow_nobs (const int newobs, double ***pZ, DATAINFO *pdinfo)
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

	if (_allocate_case_markers(&S, n + newobs)) return E_ALLOC;
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

int dataset_add_vars (const int newvars, double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char **label;
    char *vector;
    int i, n = pdinfo->n, v = pdinfo->v;    

/*      printf("grow_Z: n = %d, v = %d, newvars = %d\n", n, v, newvars);  */
/*      printf("Z size wanted: %d\n", (v + newvars) * n * sizeof(double));  */

    newZ = realloc(*pZ, (v + newvars) * sizeof *newZ);  
    if (newZ == NULL) return E_ALLOC;
    for (i=0; i<newvars; i++) {
	newZ[v+i] = malloc(n * sizeof **newZ);
	if (newZ[v+i] == NULL) return E_ALLOC;
    }
    *pZ = newZ;

    varname = realloc(pdinfo->varname, (v + newvars) * sizeof(char *));
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;
    for (i=0; i<newvars; i++) {
	pdinfo->varname[v+i] = malloc(9);
	if (pdinfo->varname[v+i] == NULL) return E_ALLOC;
	strcpy(pdinfo->varname[v+i], "");
    }

    if (pdinfo->label != NULL) {
	label = realloc(pdinfo->label, (v + newvars) * sizeof(char *));
	if (label == NULL) return E_ALLOC;
	else pdinfo->label = label;
	for (i = 0; i<newvars; i++) {
	    pdinfo->label[v+i] = malloc(MAXLABEL);
	    if (pdinfo->label[v+i] == NULL) return E_ALLOC;
	    strcpy(pdinfo->label[v+i], "");
	}
    }

    vector = realloc(pdinfo->vector, (v + newvars));
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;
    for (i=0; i<newvars; i++) 
	pdinfo->vector[v+i] = 1;

    pdinfo->v += newvars;
    return 0;
}

/* ......................................................  */

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char **label;
    char *vector;
    int n = pdinfo->n, v = pdinfo->v;    

    newZ = realloc(*pZ, (v + 1) * sizeof *newZ);  
    if (newZ == NULL) return E_ALLOC;
    newZ[v] = malloc(n * sizeof **newZ);
    if (newZ[v] == NULL) return E_ALLOC;
    *pZ = newZ;

    varname = realloc(pdinfo->varname, (v + 1) * sizeof(char *));
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;
    pdinfo->varname[v] = malloc(9);
    if (pdinfo->varname[v] == NULL) return E_ALLOC;
    pdinfo->varname[v][0] = '\0';

    if (pdinfo->label != NULL) {
	label = realloc(pdinfo->label, (v + 1) * sizeof(char *));
	if (label == NULL) return E_ALLOC;
	else pdinfo->label = label;
	pdinfo->label[v] = malloc(MAXLABEL);
	if (pdinfo->label[v] == NULL) return E_ALLOC;
	pdinfo->label[v][0] = '\0';
    }

    vector = realloc(pdinfo->vector, (v + 1));
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;
    pdinfo->vector[v] = 0;

    pdinfo->v += 1;
    return 0;
}

/* ......................................................  */

int dataset_drop_var (int varno, double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char **label;
    char *vector;
    int i, v = pdinfo->v; 

    free(pdinfo->varname[varno]);
    if (pdinfo->label[varno] != NULL) 
	free(pdinfo->label[varno]);
    free((*pZ)[varno]);

    for (i=varno; i<v-1; i++) {
	pdinfo->varname[i] = pdinfo->varname[i+1];
	pdinfo->label[i] = pdinfo->label[i+1];
	(*pZ)[i] = (*pZ)[i+1];
    }

    varname = realloc(pdinfo->varname, (v-1) * sizeof(char *));
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, (v-1));
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;

    label = realloc(pdinfo->label, (v-1) * sizeof(char *));
    if (label == NULL) return E_ALLOC;
    else pdinfo->label = label;

    newZ = realloc(*pZ, (v-1) * sizeof *newZ); 
    if (newZ == NULL) return E_ALLOC;
    else *pZ = newZ;

    pdinfo->v -= 1;

    return 0;
}

/* ......................................................  */

int dataset_drop_vars (const int delvars, double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    char **varname;
    char **label;
    char *vector;
    int i, v = pdinfo->v;   

    if (delvars <= 0) return 0;

    for (i=v-delvars; i<v; i++) {
	if (pdinfo->varname[i] != NULL) free(pdinfo->varname[i]);
	if (pdinfo->label[i] != NULL) free(pdinfo->label[i]);
	if ((*pZ)[i] != NULL) free((*pZ)[i]);
    }

    newZ = realloc(*pZ, (v - delvars) * sizeof *newZ); 
    if (newZ == NULL) return E_ALLOC;
    else *pZ = newZ;
        
    varname = realloc(pdinfo->varname, (v - delvars) * sizeof(char *));
    if (varname == NULL) return E_ALLOC;
    else pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, (v - delvars));
    if (vector == NULL) return E_ALLOC;
    else pdinfo->vector = vector;

    label = realloc(pdinfo->label, (v - delvars) * sizeof(char *));
    if (label == NULL) return E_ALLOC;
    else pdinfo->label = label;

    pdinfo->v -= delvars;
    return 0;
}

/* ........................................................... */

int hidden_var (const int i, const DATAINFO *pdinfo)
{
    if (strcmp(pdinfo->varname[i], "subdum") == 0 ||
	strcmp(pdinfo->varname[i], "annual") == 0 ||
	strcmp(pdinfo->varname[i], "qtrs") == 0 ||
	strcmp(pdinfo->varname[i], "months") == 0 ||
	strcmp(pdinfo->varname[i], "hours") == 0)
	return 1;
    return 0;
}

/* ........................................................... */

double *copyvec (const double *src, const int n)
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

int copy_model (MODEL *targ, const MODEL *src, const DATAINFO *pdinfo)
{
    int i = src->list[0] - 1;
    int m = i * (i + 1) / 2;

    /* monolithic copy of structure */
    *targ = *src;

    /* now work on pointer members */
    _init_model(targ, pdinfo);
    if ((targ->coeff = copyvec(src->coeff, src->ncoeff + 1)) == NULL)
	return 1;
    if ((targ->sderr = copyvec(src->sderr, src->ncoeff + 1))  == NULL)  
	return 1;
    if ((targ->uhat = copyvec(src->uhat, pdinfo->n)) == NULL) return 1;
    if ((targ->yhat = copyvec(src->yhat, pdinfo->n)) == NULL) return 1;
    if ((targ->xpx = copyvec(src->xpx, m + 1)) == NULL) return 1;
    if (src->subdum != NULL && 
	(targ->subdum = copyvec(src->subdum, pdinfo->n)) == NULL) return 1;
    if (src->vcv != NULL && 
	(targ->vcv = copyvec(src->vcv, m + 1)) == NULL) return 1;
    if (src->arlist != NULL) {
      	if ((targ->rhot = copyvec(src->rhot, src->arlist[0])) == NULL) 
	    return 1; 
	m = src->arlist[0];
	targ->arlist = malloc((m + 1) * sizeof(int));
	if (targ->arlist == NULL) return 1;
	for (i=0; i<=m; i++) targ->arlist[i] = src->arlist[i];
    }
    if (src->slope != NULL &&
	(targ->slope = copyvec(src->slope, src->ncoeff + 1)) == NULL)
	    return 1;

    m = src->list[0];
    targ->list = malloc((m + 1) * sizeof(int));
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

int _forecast (int t1, const int t2, const int nv, 
	       const MODEL *pmod, DATAINFO *pdinfo, double ***pZ)
{
    double xx, zz, zr;
    int i, k, maxlag = 0, yno = pmod->list[1], ARMODEL;
    int v, t;

    ARMODEL = (pmod->ci == AR || pmod->ci == CORC || 
	       pmod->ci == HILU)? 1: 0;
    if (ARMODEL) {
	maxlag = pmod->arlist[pmod->arlist[0]];
	if (t1 < maxlag) t1 = maxlag; 
    }
    for (t=t1; t<=t2; t++) {
	zz = 0.0;
	if (ARMODEL) for (k=1; k<=pmod->arlist[0]; k++) {
	    xx = (*pZ)[yno][t-pmod->arlist[k]];
	    zr = pmod->rhot[k];
	    if (na(xx)) {
		if (zr == 0) continue;
		xx = (*pZ)[nv][t-pmod->arlist[k]];
		if (na(xx)) {
		    (*pZ)[nv][t] = NADBL;
		    goto ENDIT;
		}
	    }
	    zz = zz + xx * zr;
	}
	for (v=1; v<=pmod->ncoeff; v++) {
	    k = pmod->list[v+1];
	    xx = (*pZ)[k][t];
	    if (na(xx)) {
		zz = NADBL;
		break;
	    }
	    if (ARMODEL) {
		xx = (*pZ)[k][t];
		for (i=1; i<=pmod->arlist[0]; i++) 
		    xx -= pmod->rhot[i] * (*pZ)[k][t-pmod->arlist[i]];
	    }
	    zz = zz + xx * pmod->coeff[v];
	}
	(*pZ)[nv][t] = zz;
    ENDIT:  ;
    }
    return 0;
}

/* ........................................................... */

int _full_model_list (MODEL *pmod, int **plist)
/* reconstitute full varlist for WLS and AR models */
{
    int i, pos = 0, len, *mylist;

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
	pos = pmod->arlist[0] + 1;
	len = pos + (*plist)[0] + 2;
	mylist = malloc(len * sizeof *mylist);
	if (mylist == NULL) return -1;
	mylist[0] = len - 2;
	for (i=1; i<pos; i++) mylist[i] = pmod->arlist[i];
	mylist[pos] = 999;
	for (i=1; i<=(*plist)[0]; i++) 
	    mylist[pos+i] = (*plist)[i];
    }
    copylist(plist, mylist);
    free(mylist);
    return pos;
}

/* ........................................................... */

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn,
		     const PATHS *ppaths, const int plot)
     /* use Salkever's method to generate forecasts plus forecast
	variances -- FIXME ifc = 0, and methods other than OLS */
{
    double **fZ;
    DATAINFO fdatainfo;
    MODEL fmod; 
    int *list, orig_v, idate, ft1, ft2, v1, err = 0;
    int i, j, k, t, nfcast, fn, fv;
    double xdate, tval, maxerr, *yhat, *sderr, *depvar;
    char t1str[8], t2str[8];

    if (pmod->ci != OLS || !pmod->ifc) return E_OLSONLY;

    /* parse dates */
    if (sscanf(str, "%*s %7s %7s", t1str, t2str) != 2) 
	return E_OBS; 
    ft1 = dateton(t1str, pdinfo->pd, pdinfo->stobs);
    ft2 = dateton(t2str, pdinfo->pd, pdinfo->stobs);
    if (ft1 < 0 || ft2 < 0 || ft2 < ft1) return E_OBS;

    orig_v = pmod->list[0];
    v1 = pmod->list[1];

    /* FIXME -- need to handle case where orig. regression was
       run on a sub-sample */

    nfcast = ft2 - ft1 + 1;
    fn = fdatainfo.n = nfcast + pdinfo->t2 + 1;
    fv = fdatainfo.v = nfcast + orig_v;

    fZ = malloc(fv * sizeof *fZ);
    if (fZ != NULL) 
	for (i=0; i<fv; i++) 
	    fZ[i] = malloc(fn * sizeof **fZ);
    list = malloc((fv + 1) * sizeof *list);
    yhat = malloc(nfcast * sizeof *yhat);
    sderr = malloc(nfcast * sizeof *sderr);
    depvar = malloc(nfcast * sizeof *depvar);
    if (fZ == NULL || list == NULL || yhat == NULL 
	|| sderr == NULL || depvar == NULL) return E_ALLOC;

    strcpy(fdatainfo.stobs, pdinfo->stobs);
    fdatainfo.t1 = pdinfo->t1;
    fdatainfo.t2 = fn - 1;
    fdatainfo.varname = NULL;
    fdatainfo.label = NULL;
    fdatainfo.S = NULL;

    /* start new list */
    list[0] = fv;
    for (i=1; i<=list[0]; i++) list[i] = i;
    if (pmod->ifc) list[list[0]] = 0;
    /*  printlist(list); */

    /* set entire new data matrix to zero */
    for (i=0; i<fv; i++)
	for (t=0; t<fn; t++) fZ[i][t] = 0.0;
    /* insert const at pos. 0 */
    for (t=0; t<fn; t++) fZ[0][t] = 1.0;
    /* insert orig model vars into fZ */
    k = pmod->ifc? orig_v-1: orig_v;
    for (i=1; i<=k; i++) {
	for (t=0; t<=pdinfo->t2; t++) 
	    fZ[i][t] = (*pZ)[pmod->list[i]][t];
	if (i == 1) continue;
	for (t=pdinfo->t2+1; t<fn; t++)
	    fZ[i][t] = (*pZ)[pmod->list[i]][t - (pdinfo->t2+1) + ft1];
    }
    /* insert -I section */
    for (i=orig_v; i<fv; i++) {
	k = orig_v - i;
	for (t=pdinfo->t2+1; t<fn; t++) {
	    j = pdinfo->t2 + 1 - t;
	    if (k == j) fZ[i][t] = -1.0;
	}
    }

#ifdef notdef
    /* check: print matrix */
    for (t=0; t<fn; t++) {
 	for (i=0; i<fv; i++)
 	    printf("%.2f ", fZ[i][t]);
 	putc('\n', stdout);
    }
#endif
    
    _init_model(&fmod, &fdatainfo);
    fdatainfo.extra = 1;
    fmod = lsq(list, &fZ, &fdatainfo, OLS, 1, 0.0);
    if (fmod.errcode) {
	err = fmod.errcode;
	clear_model(&fmod, NULL, NULL, &fdatainfo);
	free_Z(fZ, &fdatainfo);
	free(list);
	free(yhat);
	free(sderr);
	free(depvar);
	return err;
    }

    /* find the fitted values */
    t = 0;
    for (i=orig_v-1; i<fv-1; i++) {
	yhat[t] = fmod.coeff[i];
	t++;
    }    

    /* and the variances */
    if (makevcv(&fmod)) return E_ALLOC;

#ifdef notdef
    nv = fv - 1;
    k = (nv * nv + nv)/2;
    for (i=0; i<k; i++)
 	printf("vcv[%d] = %f\n", i, fmod.vcv[i]);
#endif

    k = -1;
    t = 0;
    for (i=1; i<fv; i++) {
	for (j=1; j<fv; j++) {
	    if (j < i) continue;
	    k++;
	    if (j == i && i < fv - 1 && i > orig_v - 2) {
		sderr[t] = sqrt(fmod.vcv[k]);
		t++;
	    }
	}
    } 

    /* print results */
    for (t=0; t<nfcast; t++) 
	depvar[t] = (*pZ)[v1][ft1 + t];
    tval = _tcrit95(pmod->dfd);
    pprintf(prn, " For 95%% confidence intervals, t(%d, .025) = %.3f\n", 
	    pmod->dfd, tval);
    if (pdinfo->pd == 1) pprintf(prn, "\n Obs ");
    else pprintf(prn, "\n\n     Obs ");
    pprintf(prn, "%13s", pdinfo->varname[v1]);
    pprintf(prn, "%13s", "prediction");
    pprintf(prn, "%14s", " std. error");
    pprintf(prn, "   95%% confidence interval\n");
    pprintf(prn, "\n");

    for (t=0; t<nfcast; t++) {
	if (pdinfo->markers) { 
	    pprintf(prn, "%8s ", pdinfo->S[t]); 
	} else {
	    xdate = date(t + ft1, pdinfo->pd, pdinfo->sd0);
	    idate = (int) xdate;
	    if (pdinfo->pd == 1) pprintf(prn, "%4d ", idate);
	    else if (pdinfo->pd < 10) pprintf(prn, "%8.1f ", xdate);
	    else pprintf(prn, "%8.2f ", xdate);
	}
	_printxs(depvar[t], 15, PRINT, prn);
	_printxs(yhat[t], 15, PRINT, prn);
	_printxs(sderr[t], 15, PRINT, prn);
	maxerr = tval * sderr[t];
	_printxs(yhat[t] - maxerr, 15, PRINT, prn);
	pprintf(prn, " -");
	_printxs(yhat[t] + maxerr, 10, PRINT, prn);
	pprintf(prn, "\n");
	sderr[t] = maxerr;
    }

    if (plot) {
	if (pdinfo->time_series == TIME_SERIES) {
	    switch (pdinfo->pd) {
	    case 1:
		plotvar(pZ, pdinfo, "annual");
		break;
	    case 4:
		plotvar(pZ, pdinfo, "qtrs");
		break;
	    case 12:
		plotvar(pZ, pdinfo, "months");
		break;
	    case 24:
		plotvar(pZ, pdinfo, "hours");
		break;
	    default:
		plotvar(pZ, pdinfo, "time");
	    }
	} else plotvar(pZ, pdinfo, "index");
	err = plot_fcast_errs(nfcast, &(*pZ)[pdinfo->v - 1][ft1], 
			      depvar, yhat, sderr, pdinfo->varname[v1], ppaths);
    }

    clear_model(&fmod, NULL, NULL, &fdatainfo);
    free_Z(fZ, &fdatainfo);
    free(list);
    free(yhat);
    free(sderr);
    free(depvar);
    clear_datainfo(&fdatainfo, 0);

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
	!strncmp(line, "omit", 4) ||
	!strncmp(line, "tsls", 4) ||
	!strncmp(line, "logit", 5)    ||
	!strncmp(line, "probit", 6)   ||
	!strncmp(line, "ar", 2))
	return 1;
    return 0;
}

/* ........................................................... */

int is_model_ref_cmd (const int ci)
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
    sprintf(spec->cmd, "%s ", commands[pmod->ci]);
    
    if (pmod->ci == AR) {
	store_list(pmod->arlist, spec->cmd);
	strcat(spec->cmd, "; ");
    }
    store_list(pmod->list, spec->cmd);

    if (pmod->subdum != NULL) {
	int t;

	spec->subdum = malloc(fullinfo->n * sizeof(double));
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

    command.list = malloc(sizeof(int));
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
			command.ci, &prn);
	if (!err)
	    *tmpmod = lsq(command.list, pZ, pdinfo, command.ci, 0, rho);
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
	clear_model(tmpmod, NULL, NULL, pdinfo);
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

int set_panel_structure (int flag, DATAINFO *pdinfo, PRN *prn)
{
    if (pdinfo->pd == 1) {
	pprintf(prn, "The current data frequency, 1, is not "
		"compatible with panel data.\nPlease see the 'setobs' "
		"command.\n");
	return 1;
    }

    if (flag) {
	pprintf(prn, "Panel structure set to stacked cross sections\n");
	pdinfo->time_series = STACKED_CROSS_SECTION;
    } else {
	pprintf(prn, "Panel structure set to stacked time series\n");
	pdinfo->time_series = STACKED_TIME_SERIES;
    }

    return 0;
}

/* ........................................................... */

int balanced_panel (const DATAINFO *pdinfo)
{
    char unit[9], period[9];

    if ((pdinfo->t2 - pdinfo->t1 + 1) % pdinfo->pd)
        return 0;

    if (sscanf(pdinfo->endobs, "%[^.].%s", unit, period) == 2) {
        if (atoi(period) != pdinfo->pd)
            return 0;
    } else 
        return 0;

    return 1;
}




