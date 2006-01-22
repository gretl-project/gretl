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
#include "objstack.h"
#include "cmd_private.h"
#include "libset.h"
#include "usermat.h"

#include <errno.h>

#ifndef WIN32
# include <signal.h>
# ifdef USE_GTK2
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
	fprintf(stderr, "%d : ", list[0]);
	for (i=1; i<=list[0]; i++) {
	    fprintf(stderr, "%d ", list[i]);
	}
    }

    fputc('\n', stderr);
}

/* Compute model selection criteria */

int gretl_calculate_criteria (double ess, int nobs, int ncoeff,
			      double *ll, double *aic, double *bic)
{
    int err = 0;

    if (na(ess) || ess <= 0.0 || ncoeff < 1 || nobs <= ncoeff) {
	*ll = NADBL;
	*aic = NADBL;
	*bic = NADBL;
	err = 1;
    } else {
	const double ln2pi1 = 2.837877066409345;

	errno = 0;

	*ll = -.5 * nobs * log(ess);

	if (errno == EDOM || errno == ERANGE) {
	    *ll = NADBL;
	    *aic = NADBL;
	    *bic = NADBL;
	} else {
	    *ll += -.5 * nobs * (ln2pi1 - log((double) nobs));
	    *aic = -2.0 * *ll + 2 * ncoeff;
	    *bic = -2.0 * *ll + ncoeff * log(nobs);
	}
    }

    return err;
}

int ls_aic_bic (MODEL *pmod)
{
    double ll, aic, bic;
    int err;

    err = gretl_calculate_criteria(pmod->ess, pmod->nobs, pmod->ncoeff,
				   &ll, &aic, &bic);

    pmod->lnL = ll;
    pmod->criterion[C_AIC] = aic;
    pmod->criterion[C_BIC] = bic;

    return err;
}

int gretl_print_criteria (double ess, int nobs, int ncoeff, PRN *prn)
{
    double ll, aic, bic;
    int err;

    err = gretl_calculate_criteria(ess, nobs, ncoeff, &ll, &aic, &bic);

    if (err) {
	pputs(prn, _("Error calculating model selection criteria\n"));
    } else {
	pprintf(prn, _("Using ess = %g, %d observations, %d coefficients\n"), 
		ess, nobs, ncoeff);
	pprintf(prn, "\nAIC = %g\nBIC = %g\n\n", aic, bic);
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

/**
 * gretl_integer_from_string:
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

int gretl_integer_from_string (const char *s, const double **Z, 
			       const DATAINFO *pdinfo, int *err)
{
    char *test;
    int n = 0;

    n = strtol(s, &test, 10);
    if (*test == '\0') {
	return n;
    } else {
	int v = varindex(pdinfo, s);
	double x;

	if (v < pdinfo->v && !pdinfo->vector[v]) {
	    x = Z[v][0];
	    if (na(x)) {
		*err = 1;
	    } else {
		n = (int) x;
	    }
	} else {
	    *err = 1;
	}
    }

    return n;    
}

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

int varnum_from_string (const char *str, DATAINFO *pdinfo)
{
    int varno = positive_int_from_string(str);

    if (varno <= 0 || varno >= pdinfo->v) {
	varno = -1;
    } 
    
    return varno;
}

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
    double rho = 0.0;
    int err;

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
	*tmpmod = ar_func(cmd.list, atoi(cmd.param), pZ, 
			  pdinfo, OPT_NONE, NULL);
	break;
    case CORC:
    case HILU:
    case PWE:
	rho = estimate_rho(cmd.list, pZ, pdinfo, cmd.ci, 
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
	*tmpmod = logit_probit(cmd.list, pZ, pdinfo, cmd.ci, cmd.opt);
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

    gretl_cmd_free(&cmd);

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

double get_xvalue (int i, const double **Z, const DATAINFO *pdinfo)
{
    if (pdinfo->vector[i]) {
	return Z[i][pdinfo->t1];
    } else {
	return Z[i][0];
    }	
}

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
	mpvals->varnames[i] = malloc(VNAMELEN + 8);
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
    if (mpvals == NULL) {
	return NULL;
    }

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

    for (i=0; i<nc; i++) {
	mpvals->coeff[i] = NADBL;
	mpvals->sderr[i] = NADBL;
    }

    mpvals->sigma = mpvals->ess = NADBL;
    mpvals->rsq = mpvals->fstt = NADBL;
    mpvals->adjrsq = NADBL;

    mpvals->t1 = mpvals->t2 = mpvals->ifc = 0;
    mpvals->dfn = mpvals->dfd = 0;

    return mpvals;
}

#ifndef WIN32
# ifdef USE_GSPAWN

static int font_not_found (const char *s)
{
    /* "Could not find/open font when opening font X, using default" */

    if (strstr(s, "using default")) {
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

void libgretl_cleanup (void)
{
    const char *p;

    gretl_rand_free();
    gretl_functions_cleanup();
    gretl_saved_objects_cleanup();
    gretl_transforms_cleanup();
    libset_cleanup();
    gretl_lists_cleanup();
    gretl_command_hash_cleanup();
    destroy_user_matrices();

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
