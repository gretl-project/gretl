/*
 *  Copyright (c) 2004 by Allin Cottrell
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
#include "gretl_private.h"

#undef TRDEBUG

/* newlag is a library global, used when auto-generating a lag in the
   context of reading a regression list (interact.c)
*/
int newlag;

enum {
    VARS_IDENTICAL,
    X_HAS_MISSING,
    Y_HAS_MISSING,
    VARS_DIFFER
} varcomp_codes;

static int check_vals (const double *x, const double *y, int n)
{
    int ret = VARS_IDENTICAL;
    int t;

    for (t=0; t<n; t++) {
	if (floatneq(x[t], y[t])) {
	    if (na(x[t]) && !na(y[t])) {
		if (ret == VARS_IDENTICAL || ret == X_HAS_MISSING) {
		    ret = X_HAS_MISSING;
		} else {
		    ret = VARS_DIFFER;
		}
	    } else if (!na(x[t]) && na(y[t])) {
		if (ret == VARS_IDENTICAL || ret == Y_HAS_MISSING) {
		    ret = Y_HAS_MISSING;
		} else {
		    ret = VARS_DIFFER;
		}
	    } else {
		ret = VARS_DIFFER;
	    }
	    
	}
	if (ret == VARS_DIFFER) {
	    break;
	}
    }

    return ret;
}

static int 
make_transform_varname (char *vname, const char *orig, int cmd, 
			int lag, int len)
{
    *vname = '\0';

    if (cmd == DIFF) {
	strcpy(vname, "d_");
	strncat(vname, orig, len - 2);
    } else if (cmd == LDIFF) {
	strcpy(vname, "ld_");
	strncat(vname, orig, len - 3);
    } else if (cmd == LOGS) {
	strcpy(vname, "l_");
	strncat(vname, orig, len - 2);
    } else if (cmd == SQUARE) {
	strcpy(vname, "sq_");
	strncat(vname, orig, len - 3);
    } else if (cmd == LAGS) {
	char ext[6];

	sprintf(ext, "_%d", lag);
	strncat(vname, orig, len - strlen(ext));
	strcat(vname, ext);
    }

    return 0;
}

static int
make_transform_label (char *label, const char *parent,
		      int cmd, int lag)
{
    int err = 0;

    if (cmd == DIFF) {
	sprintf(label, _("= first difference of %s"), parent);
    } else if (cmd == LDIFF) {
	sprintf(label, _("= log difference of %s"), parent);
    } else if (cmd == LOGS) {
	sprintf(label, _("= log of %s"), parent);
    } else if (cmd == SQUARE) {
	sprintf(label, _("= %s squared"), parent);
    } else if (cmd == LAGS) {
	sprintf(label, "= %s(t - %d)", parent, lag);
    } else {
	err = 1;
    }

    return err;
}

static void make_xp_varname (char *vname, const char *v1, 
			     const char *v2, int len)
{
    int v2len = len / 2;
    int v1len = (len % 2)? v2len : v2len - 1;

    *vname = '\0';
    strncat(vname, v1, v1len);
    strcat(vname, "_");
    strncat(vname, v2, v2len);
}

/* Array into which to write a generated variable, prior to testing
   whether or not the same var already exists.  Declared as non-static
   so that it can be called in context of library cleanup (with n =
   0).
*/

double *testvec (int n)
{
    static double *x;
    static int nx;
    int t;

    if (n == 0) {
	/* clean up */
	free(x);
	x = NULL;
	nx = 0;
	return NULL;
    }	

    if (n > nx) {
	double *y = realloc(x, n * sizeof *x);
	
	if (y == NULL) {
	    free(x);
	    x = NULL;
	    nx = 0;
	    return NULL;
	}

	nx = n;
	x = y;
    }

    for (t=0; t<n; t++) {
	x[t] = NADBL;
    }    

    return x;
}

/* write lagged values of variable v into lagvec */

static int get_lag (int v, int lag, double *lagvec, 
		    const double **Z, const DATAINFO *pdinfo)
{
    int t, t1, lt;

    t1 = (lag > pdinfo->t1)? lag : pdinfo->t1;

    for (t=0; t<pdinfo->n; t++) {
	lagvec[t] = NADBL;
    }

    if (pdinfo->structure == STACKED_CROSS_SECTION) {
	/* needs rather special handling */
	for (t=t1; t<=pdinfo->t2; t++) { 
	    lt = t - lag * pdinfo->pd;
	    if (lt < 0 || lt >= pdinfo->n) {
		continue;
	    }
	    lagvec[t] = Z[v][lt];
	}
    } else if (dated_daily_data(pdinfo)) {
	for (t=t1; t<=pdinfo->t2; t++) {
	    lt = t - lag;
	    while (lt >= 0 && na(Z[v][lt])) {
		lt--;
	    }
	    lagvec[t] = Z[v][lt];
	}
    } else { 
	/* the "standard" time-series case */
	for (t=t1; t<=pdinfo->t2; t++) {
	    lt = t - lag;
	    if (lt < 0 || lt >= pdinfo->n) {
		continue;
	    }
	    lagvec[t] = Z[v][lt];
	}
    }

    /* post-process missing panel values */
    if (pdinfo->structure == STACKED_TIME_SERIES) {
	char *p, obs[OBSLEN];
	int j;

	for (t=t1; t<=pdinfo->t2; t++) {
	    ntodate(obs, t, pdinfo);
	    p = strchr(obs, ':');
	    j = atoi(p + 1);
	    if (j <= lag) {
		lagvec[t] = NADBL;
	    }
	}
    }

    return 0;
}

/* write log of variable v into logvec */

static int get_log (int v, double *logvec, const double **Z, 
		    const DATAINFO *pdinfo)
{
    double xx;
    int t, err = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2 && !err; t++) {
	xx = (pdinfo->vector[v])? Z[v][t] : Z[v][0];
	if (na(xx)) {
	    logvec[t] = NADBL;
	} else if (xx <= 0.0) {
	    sprintf(gretl_errmsg, 
		    _("Log error: Variable '%s', obs %d,"
		      " value = %g\n"), pdinfo->varname[v],
		    t+1, xx);
	    err = 1;
	} else {
	    logvec[t] = log(xx); 
	}
    }

    return err;
}

/* write first diff of variable v into diffvec */

static int get_diff (int v, double *diffvec, int ci,
		     const double **Z, const DATAINFO *pdinfo)
{
    double x0, x1;
    int t, t1;

    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;

    for (t=t1; t<=pdinfo->t2; t++) {
	if (pdinfo->structure == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}
	x0 = Z[v][t];
	if (pdinfo->structure == STACKED_CROSS_SECTION) {
	    x1 = (t - pdinfo->pd >= 0)? Z[v][t-pdinfo->pd] : NADBL;
	} else {
	    x1 = Z[v][t-1];
	}
	if (ci == LDIFF) {
	    if (!na(x0) && !na(x1) && x0 > 0 && x1 > 0) {
		diffvec[t] = log(x0) - log(x1);
	    }
	} else {
	    if (!na(x0) && !na(x1)) {
		diffvec[t] = x0 - x1;
	    }
	}
    }

    return 0;
}

/* write fractional difference of variable v into diffvec */

static int get_fracdiff (int v, double *diffvec, double d,
			 const double **Z, const DATAINFO *pdinfo)
{
    int dd, t, t1, T;
    const double TOL = 1.0E-07;
    double phi = -d;

    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    T = pdinfo->t2 - t1 + 1;

    for (t=t1; t<=pdinfo->t2; t++) {
	diffvec[t] = Z[v][t];
    }

    dd = 1;

    while ((dd < T) && fabs(phi) > TOL) {
        for (t = dd; t < T; t++) {
	    diffvec[t] += phi * Z[v][t - dd];
	}
        phi *= (dd - d) / (dd + 1);
	dd++;
    }


    return 0;
}

/* write square or cross-product into xvec */

static int get_xpx (int vi, int vj, double *xvec, const double **Z, 
		    const DATAINFO *pdinfo)
{
    double xit, xjt;
    int t;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xit = (pdinfo->vector[vi])? Z[vi][t] : Z[vi][0];
	xjt = (pdinfo->vector[vj])? Z[vj][t] : Z[vj][0];
	if (na(xit) || na(xjt)) {
	    xvec[t] = NADBL;
	} else {
	    xvec[t] = xit * xjt;
	}
    }

    return 0;
}

enum transform_results {
    VAR_ADDED_OK,
    VAR_EXISTS_OK,
    VAR_ADD_FAILED,
    VARNAME_DUPLICATE
};

static int 
check_add_transform (int vnum, const double *x,
		     const char *vname, const char *label,
		     DATAINFO *pdinfo, double ***pZ)
{
    int t, ret = VAR_ADDED_OK;

    if (vnum < pdinfo->v) {
	/* a variable of this name already exists */
	int chk = check_vals(x, (*pZ)[vnum], pdinfo->n);

	/* heuristic: we'll assume that if the variables have the same
	   name, and differ only in respect of one being more complete
	   than the other (having valid values where the other has
	   missing values), then we should not create a new variable,
	   but rather use the more complete series.

	   Watch out for cases where this is not the desired behavior!
	*/

	if (chk == VARS_IDENTICAL) {
	    ret = VAR_EXISTS_OK;
	} else if (chk == X_HAS_MISSING) { /* is this right? */
	    ret = VAR_EXISTS_OK;
	} else if (chk == Y_HAS_MISSING) {
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[vnum][t] = x[t];
	    }
	    ret = VAR_EXISTS_OK;
	} else {
	    if (!strcmp(label, VARLABEL(pdinfo, vnum))) {
		/* labels match: update the values */
		for (t=0; t<pdinfo->n; t++) {
		    (*pZ)[vnum][t] = x[t];
		}
		ret = VAR_EXISTS_OK;
	    } else {
		/* labels do not match: problem */
		ret = VARNAME_DUPLICATE;
	    }
	}
    } else {
	/* no var of this name, working from scratch */
	if (dataset_add_vars(1, pZ, pdinfo)) {
	    ret = VAR_ADD_FAILED;
	} else {
	    strcpy(pdinfo->varname[vnum], vname);
	    strcpy(VARLABEL(pdinfo, vnum), label);
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[vnum][t] = x[t];
	    }
	}
    }

    return ret;
}

/* get_transform: create specified transformation of variable v if
   this variable does not already exist.

   Return the ID number of the transformed var, or -1 on error.
*/

static int get_transform (int ci, int v, int aux, 
			  double ***pZ, DATAINFO *pdinfo,
			  int startlen)
{
    char vname[VNAMELEN];
    char label[MAXLABEL];
    int len, vno = -1, err = 0;
    double *vx;

    vx = testvec(pdinfo->n);
    if (vx == NULL) {
	return -1;
    }

    if (ci == LAGS) {
	/* "aux" = lag */
	err = get_lag(v, aux, vx, (const double **) *pZ, pdinfo);
    } else if (ci == LOGS) {
	err = get_log(v, vx, (const double **) *pZ, pdinfo);
    } else if (ci == DIFF || ci == LDIFF) {
	err = get_diff(v, vx, ci, (const double **) *pZ, pdinfo);
    } else if (ci == SQUARE) {
	/* "aux" = second variable number */
	err = get_xpx(v, aux, vx, (const double **) *pZ, pdinfo);
    }

#if 0
    } else if (ci == FRACDIFF) {
	err = get_fracdiff(v, vx, d, (const double **) *pZ, pdinfo);
    }
#endif

    if (err) {
	return -1;
    }

    if (ci == SQUARE && v != aux) {
	sprintf(label, _("= %s times %s"), pdinfo->varname[v], 
		pdinfo->varname[aux]);
    } else {
	make_transform_label(label, pdinfo->varname[v], ci, aux);
    }
    
    for (len=startlen; len<=VNAMELEN; len++) {

	if (len == VNAMELEN) {
	    /* last resort: hack the name */
	    make_varname_unique(vname, 0, pdinfo);
	} else if (ci == SQUARE && v != aux) {
	    make_xp_varname(vname, pdinfo->varname[v],
			    pdinfo->varname[aux], len);
	} else {
	    make_transform_varname(vname, pdinfo->varname[v], ci, 
				   aux, len);
	}
	vno = varindex(pdinfo, vname);

	err = check_add_transform(vno, vx, vname, label, pdinfo, pZ);
	if (err != VAR_EXISTS_OK && err != VAR_ADDED_OK) {
	    vno = -1;
	}
	if (err != VARNAME_DUPLICATE) {
	    /* either we're OK, or there was a fatal error */
	    break;
	}
    }

    return vno;
}

/* laggenr: create specified lag of variable v if this variable does
   not already exist.

   Return the ID number of the lag var, or -1 on error.
*/

int laggenr (int v, int lag, double ***pZ, DATAINFO *pdinfo)
{
    int lno, oldv = pdinfo->v;

    if (!pdinfo->vector[v]) {
	return -1;
    }

    /* sanity check on lag length */
    if (lag > pdinfo->n) {
	return -1;
    }

    newlag = 1;

    lno = get_transform(LAGS, v, lag, pZ, pdinfo, 8);

    if (lno < oldv) {
	newlag = 0;
    }

    return lno;
}

/* loggenr: create log of variable v if this variable does not
   already exist.  

   Return the ID number of the log var, or -1 on error.
*/

int loggenr (int v, double ***pZ, DATAINFO *pdinfo)
{
    return get_transform(LOGS, v, 0, pZ, pdinfo, 8);
}

/* diffgenr: create first difference (or log difference) of variable v,
   if this variable does not already exist.

   Return the ID number of the differenced var, or -1 on error.
*/

int diffgenr (int v, double ***pZ, DATAINFO *pdinfo, int ldiff)
{
    if (!pdinfo->vector[v]) {
	return -1;
    }

    return get_transform((ldiff)? LDIFF : DIFF, v, 0, pZ, pdinfo, 8);
}

/* xpxgenr: create square or cross-product, if the target variable
   does not already exist.

   Return the ID number of the squared or cross var, or -1 on error.
*/

int xpxgenr (int vi, int vj, double ***pZ, DATAINFO *pdinfo)
{
    if (vi == vj) {
	if (isdummy((*pZ)[vi], pdinfo->t1, pdinfo->t2)) {
	    return -1;
	}
    }

    return get_transform(SQUARE, vi, vj, pZ, pdinfo, 8);
}

static int 
get_starting_length (const int *list, DATAINFO *pdinfo, int trim)
{
    int width = 8 - trim;
    int len, maxlen = 0;
    int i, j;

    for (i=1; i<=list[0]; i++) {
	len = strlen(pdinfo->varname[list[i]]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    if (maxlen <= width) {
	/* no problem: generated names will fit in 8 chars */
	return 8;
    }

    for (len=width; len<=maxlen; len++) {
	int conflict = 0;

	for (i=1; i<=list[0]; i++) {
	    for (j=i+1; j<=list[0]; j++) {
		if (!strncmp(pdinfo->varname[list[i]],
			     pdinfo->varname[list[j]],
			     len)) {
		    conflict = 1;
		    break;
		}
	    }
	    if (conflict) {
		break;
	    }
	}
	if (!conflict) {
	    break;
	}
    }

    len += trim;

    if (len < 8) {
	len = 8;
    } else if (len > VNAMELEN) {
	len = VNAMELEN;
    }

    return len;
}

/**
 * list_loggenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set the natural logs of the
 * variables given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_loggenr (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int lognum, i, v;
    int startlen;
    int n_ok = 0;

    startlen = get_starting_length(list, pdinfo, 2);

    for (i=1; i<=list[0]; i++) {
	v = list[i];

	if (v == 0 || !pdinfo->vector[v]) {
	    continue; 
	}
	if (isdummy((*pZ)[v], pdinfo->t1, pdinfo->t2)) {
	    continue;
	}

	lognum = get_transform(LOGS, v, 0, pZ, pdinfo, startlen);
	if (lognum > 0) {
	    n_ok++;
	}
    }

    return (n_ok > 0)? 0 : E_LOGS;
}

int 
real_list_laggenr (const int *list, double ***pZ, DATAINFO *pdinfo,
		   int maxlag, int **lagnums)
{
    int lagnum, l, i, v;
    int startlen;
    int *record = NULL;

    startlen = get_starting_length(list, pdinfo, (maxlag > 9)? 3 : 2);
    
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v == 0 || !pdinfo->vector[v]) {
	    continue;
	}
	if (lagnums != NULL) {
	    record = lagnums[i-1];
	}
	for (l=1; l<=maxlag; l++) {
	    lagnum = get_transform(LAGS, v, l, pZ, pdinfo, startlen);
#if TRDEBUG
	    fprintf(stderr, "base var '%s', lag %d: lagnum = %d\n",
		    pdinfo->varname[v], l, lagnum);
#endif
	    if (lagnum < 0) {
		return 1;
	    }
#if TRDEBUG
	    fprintf(stderr, "lag var name '%s', label '%s'\n",
		    pdinfo->varname[lagnum], VARLABEL(pdinfo, lagnum));
#endif
	    if (record != NULL) {
		record[l] = lagnum;
	    }
	}
    }

    return 0;
}

/**
 * list_laggenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set lagged values of the 
 * variables given in @list (up to the frequency of the data).
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int list_laggenr (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int maxlag;

    if (pdinfo->pd < 52) {
	maxlag = pdinfo->pd;
    } else {
	maxlag = 14;
    } 

    /* play safe with panel data */
    if (dataset_is_panel(pdinfo)) {
	maxlag = 1;
    }

    return real_list_laggenr(list, pZ, pdinfo, maxlag, NULL);  
}

/**
 * list_diffgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate first-differences of the variables in @list, and add them
 * to the data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int list_diffgenr (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int i, v, startlen;

    startlen = get_starting_length(list, pdinfo, 2);
    
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (get_transform(DIFF, v, 0, pZ, pdinfo, startlen) < 0) {
	    return 1;
	}
    }

    return 0;
}

/**
 * list_ldiffgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate log-differences of the variables in @list, and add them
 * to the data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int list_ldiffgenr (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int i, v, startlen;

    startlen = get_starting_length(list, pdinfo, 3);
    
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (get_transform(LDIFF, v, 0, pZ, pdinfo, startlen) < 0) {
	    return 1;
	}
    }

    return 0;
}

/**
 * list_xpxgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: If = 0, only squares are generated, if OPT_O, both
 * squares and cross-products are generated.
 *
 * Generates and adds to the data set squares and (if @opt is non-zero) 
 * cross-products of the variables given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_xpxgenr (const int *list, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt)
{
    int xnum, i, j, vi, vj;
    int startlen;
    int n_ok = 0;

    startlen = get_starting_length(list, pdinfo, 3);

    for (i=1; i<=list[0]; i++) {
	vi = list[i];

	if (vi == 0 || !pdinfo->vector[vi]) {
	    continue; 
	}
	if (isdummy((*pZ)[vi], pdinfo->t1, pdinfo->t2)) {
	    continue;
	}

	xnum = get_transform(SQUARE, vi, vi, pZ, pdinfo, startlen);
	if (xnum > 0) {
	    n_ok++;
	}

	if (opt & OPT_O) {
	    for (j=i+1; j<=list[0]; j++) {
		vj = list[j];
		xnum = xpxgenr(vi, vj, pZ, pdinfo);
		if (xnum > 0) {
		    n_ok++;
		}
	    }
	}
    }

    return (n_ok > 0)? 0 : E_SQUARES;
}

int lagvarnum (int v, int l, const DATAINFO *pdinfo)
{
    char label[MAXLABEL];
    int i;

    make_transform_label(label, pdinfo->varname[v], LAGS, l);

#if TRDEBUG
    fprintf(stderr, "Looking for lag var with label '%s'...\n", label);
#endif

    for (i=1; i<pdinfo->v; i++) {
	if (!strcmp(label, VARLABEL(pdinfo, i))) {
	    return i;
	}
    }

    /* second attempt? */
    sprintf(label, "= %s(-%d)", pdinfo->varname[v], l);
    for (i=1; i<pdinfo->v; i++) {
	if (strstr(VARLABEL(pdinfo, i), label)) {
	    return i;
	}
    }

#if TRDEBUG
    fputs("*** not found\n", stderr);
#endif

    return -1;
}
