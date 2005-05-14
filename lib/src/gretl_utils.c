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
#include "gretl_func.h"
#include "system.h"

#include <errno.h>

#ifndef WIN32
# include <glib.h>
# include <signal.h>
# if GLIB_CHECK_VERSION(2,0,0)
#  define GRETL_GLIB 2
# endif /* GLIB_CHECK_VERSION */
#endif /* ! WIN32 */

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

int gretl_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
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

/* Compute model selection criteria */

int gretl_calculate_criteria (double *x, double ess, int nobs, int ncoeff)
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

int ls_aic_bic (MODEL *pmod)
{
    return gretl_calculate_criteria(pmod->criterion, 
				    pmod->ess, pmod->nobs,
				    pmod->ncoeff);
}

int gretl_print_criteria (double ess, int nobs, int ncoeff, PRN *prn)
{
    double x[2];
    int err;

    err = gretl_calculate_criteria(x, ess, nobs, ncoeff);

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

    if (newZ == NULL) {
	return E_ALLOC;
    }

    if (newvars == 1 && x != NULL) {
	/* new var is pre-allocated */
	newZ[v] = x;
    } else {
	for (i=0; i<newvars; i++) {
	    newZ[v+i] = malloc(n * sizeof **newZ);
	    if (newZ[v+i] == NULL) {
		return E_ALLOC;
	    }
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

    if (newZ == NULL) {
	return E_ALLOC;
    }

    newZ[v] = malloc(n * sizeof **newZ);

    if (newZ[v] == NULL) {
	return E_ALLOC;
    }

    *pZ = newZ;

    varname = realloc(pdinfo->varname, (v + 1) * sizeof *varname);

    if (varname == NULL) {
	return E_ALLOC;
    }
    pdinfo->varname = varname;

    pdinfo->varname[v] = malloc(VNAMELEN);
    if (pdinfo->varname[v] == NULL) {
	return E_ALLOC;
    }

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
    if (vector == NULL) {
	return E_ALLOC;
    }
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

int varnum_from_string (const char *str, DATAINFO *pdinfo)
{
    int varno = positive_int_from_string(str);

    if (varno <= 0 || varno >= pdinfo->v) {
	varno = -1;
    } 
    
    return varno;
}

/* ......................................................  */

static int 
shrink_dataset_to_size (double ***pZ, DATAINFO *pdinfo, int nv)
{
    char **varname;
    char *vector;
    VARINFO **varinfo;
    double **newZ;
    
    varname = realloc(pdinfo->varname, nv * sizeof *varname);
    if (varname == NULL) {
	return E_ALLOC;
    }
    pdinfo->varname = varname;

    vector = realloc(pdinfo->vector, nv * sizeof *vector);
    if (vector == NULL) {
	return E_ALLOC;
    }
    pdinfo->vector = vector;

    varinfo = realloc(pdinfo->varinfo, nv * sizeof *varinfo);
    if (varinfo == NULL) {
	return E_ALLOC;
    }
    pdinfo->varinfo = varinfo;

    newZ = realloc(*pZ, nv * sizeof *newZ); 
    if (newZ == NULL) {
	return E_ALLOC;
    }
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

/* drop all hidden (automatic) variables from the dataset */

int dataset_destroy_hidden_vars (double ***pZ, DATAINFO *pdinfo)
{
    int i, nhid = 0;
    int err = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (hidden_var(i, pdinfo)) {
	    nhid++;
	}
    }

    if (nhid > 0) {
	int *hidlist = gretl_list_new(nhid);

	if (hidlist == NULL) {
	    err = 1;
	} else {
	    int j = 1;

	    for (i=1; i<pdinfo->v; i++) {
		if (hidden_var(i, pdinfo)) {
		    hidlist[j++] = i;
		}
	    }	    
	    err = dataset_drop_listed_vars(hidlist, pZ, pdinfo, NULL);
	    free(hidlist);
	}
    }

    return err;
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
	    free_varinfo(pdinfo, i);
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
	if (na(x[i])) {
	    nmiss++;
	} else {
	    break;
	}
    }

    return nmiss;
}

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

    if (varno < 0) {
	return E_DATA;
    }

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

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo) 
{
    CMD cmd;
    int err = 0, ignore = 0;
    double rho = 0;

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
			  pdinfo, OPT_NONE, NULL);
	break;
    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, pZ, pdinfo, 1, cmd.ci, 
			   &err, cmd.opt, NULL);
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
    gretl_transforms_cleanup();

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
