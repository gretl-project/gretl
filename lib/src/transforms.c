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

#define TRDEBUG 0

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
    } else if (cmd == SDIFF) {
	strcpy(vname, "sd_");
	strncat(vname, orig, len - 3);
    } else if (cmd == LOGS) {
	strcpy(vname, "l_");
	strncat(vname, orig, len - 2);
    } else if (cmd == SQUARE) {
	strcpy(vname, "sq_");
	strncat(vname, orig, len - 3);
    } else if (cmd == LAGS) {
	char ext[6];

	if (lag >= 0) {
	    /* an actual lag */
	    sprintf(ext, "_%d", lag);
	} else {
	    /* in fact a lead */
	    sprintf(ext, "%d", -lag);
	}
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
    } else if (cmd == SDIFF) {
	sprintf(label, _("= seasonal difference of %s"), parent);
    } else if (cmd == LOGS) {
	sprintf(label, _("= log of %s"), parent);
    } else if (cmd == SQUARE) {
	sprintf(label, _("= %s squared"), parent);
    } else if (cmd == LAGS) {
	if (lag >= 0) {
	    sprintf(label, "= %s(t - %d)", parent, lag);
	} else {
	    sprintf(label, "= %s(t + %d)", parent, -lag);
	}
    } else {
	err = 1;
    }

    return err;
}

int is_standard_lag (int v, const DATAINFO *pdinfo)
{
    const char *test = VARLABEL(pdinfo, v);
    char pm, vname[VNAMELEN];
    int lag, ret = 0;

    if (sscanf(test, "= %15[^(](t %c %d)", vname, &pm, &lag) == 3) {
	ret = 1;
    }

    return ret;
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
   whether or not the same var already exists.  
*/

static double *testvec (int n)
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

void gretl_transforms_cleanup (void)
{
    testvec(0);
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

/* write some sort of difference of variable v into diffvec */

static int get_diff (int v, double *diffvec, int ci,
		     const double **Z, const DATAINFO *pdinfo)
{
    double x0, x1;
    int t, t0, t1;

    t0 = (ci == SDIFF)? pdinfo->pd : 1;
    t1 = (pdinfo->t1 > t0)? pdinfo->t1 : t0;

    for (t=t1; t<=pdinfo->t2; t++) {
	if (pdinfo->structure == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}
	x0 = Z[v][t];
	if (pdinfo->structure == STACKED_CROSS_SECTION) {
	    x1 = (t - pdinfo->pd >= 0)? Z[v][t-pdinfo->pd] : NADBL;
	} else {
	    x1 = Z[v][t - t0];
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
#if TRDEBUG
		fprintf(stderr, "check_add_transform: updating var %d (%s)\n",
			vnum, vname);
#endif
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
	if (dataset_add_series(1, pZ, pdinfo)) {
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
    } else if (ci == DIFF || ci == LDIFF || ci == SDIFF) {
	err = get_diff(v, vx, ci, (const double **) *pZ, pdinfo);
    } else if (ci == SQUARE) {
	/* "aux" = second variable number */
	err = get_xpx(v, aux, vx, (const double **) *pZ, pdinfo);
    }

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

/**
 * laggenr: 
 * @v: ID number in dataset of source variable.
 * @lag: the order of the lag to create.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 *
 * Creates the specified lag of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the lagged variable, or -1 on error.
 */

int laggenr (int v, int lag, double ***pZ, DATAINFO *pdinfo)
{
    int lno;

    if (!pdinfo->vector[v] || lag > pdinfo->n) {
	lno = -1;
    } else if (lag == 0) {
	lno = v;
    } else {
	lno = get_transform(LAGS, v, lag, pZ, pdinfo, VNAMELEN - 3);
    }

    return lno;
}

/**
 * loggenr: 
 * @v: ID number in dataset of source variable.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 *
 * Creates the natural log of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the log variable, or -1 on error.
 */

int loggenr (int v, double ***pZ, DATAINFO *pdinfo)
{
    return get_transform(LOGS, v, 0, pZ, pdinfo, VNAMELEN - 3);
}

/**
 * diffgenr: 
 * @v: ID number in dataset of source variable.
 * @ci: %DIFF (first difference), %LDIFF (log difference) or %SDIFF
 * (seasonal difference).
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 *
 * Creates the first difference (or log- or seasonal difference, 
 * depending on the value of @ci) of variable @v, if the
 * differenced variable does not already exist.
 *
 * Returns: the ID number of the differenced variable, or -1 on error.
 */

int diffgenr (int v, int ci, double ***pZ, DATAINFO *pdinfo)
{
    if (!pdinfo->vector[v]) {
	return -1;
    }

    if (ci != DIFF && ci != LDIFF && ci != SDIFF) {
	return -1;
    }

    if (ci == SDIFF && !dataset_is_seasonal(pdinfo)) {
	return -1;
    }

    return get_transform(ci, v, 0, pZ, pdinfo, VNAMELEN - 3);
}

/**
 * xpxgenr: 
 * @vi: ID number in dataset of first source variable.
 * @vj: ID number in dataset of second source variable.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 *
 * Creates the cross product of variables @vi and @vj if this 
 * variable does not already exist.
 *
 * Returns: the ID number of the cross-product variable, 
 * or -1 on error.
 */

int xpxgenr (int vi, int vj, double ***pZ, DATAINFO *pdinfo)
{
    if (vi == vj) {
	if (gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[vi])) {
	    return -1;
	}
    }

    return get_transform(SQUARE, vi, vj, pZ, pdinfo, VNAMELEN - 3);
}

static int 
get_starting_length (const int *list, DATAINFO *pdinfo, int trim)
{
    int width = VNAMELEN - 3 - trim;
    int len, maxlen = 0;
    int i, j;

    for (i=1; i<=list[0]; i++) {
	len = strlen(pdinfo->varname[list[i]]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    if (maxlen <= width) {
	/* no problem: generated names will fit in VNAMELEN - 3 chars */
	return VNAMELEN - 3;
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

    if (len < VNAMELEN - 3) {
	len = VNAMELEN - 3;
    } else if (len > VNAMELEN - 1) {
	len = VNAMELEN - 1;
    }

    return len;
}

/**
 * list_loggenr:
 * @list: on entry, list of variables to process; on exit,
 * holds the ID numbers of the generated variables.
 * @pZ: pointer to data array.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set the natural logs of the
 * variables given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_loggenr (int *list, double ***pZ, DATAINFO *pdinfo)
{
    int tnum, i, v;
    int startlen;
    int n_ok = 0;

    startlen = get_starting_length(list, pdinfo, 2);

    for (i=1; i<=list[0]; i++) {
	v = list[i];

	if (v == 0 || !pdinfo->vector[v]) {
	    continue; 
	}
	if (gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[v])) {
	    continue;
	}

	tnum = get_transform(LOGS, v, 0, pZ, pdinfo, startlen);

	if (tnum > 0) {
	    n_ok++;
	    list[i] = tnum;
	}
    }

    return (n_ok > 0)? 0 : E_LOGS;
}

static int *make_lags_list (int *list, int order, DATAINFO *pdinfo)
{
    int i, v, nl = 0;

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0 && pdinfo->vector[v]) {
	    nl += order;
	}
    }

    return gretl_list_new(nl);
}

/**
 * list_laggenr:
 * @plist: pointer to list of variables to process.  On exit
 * the list holds the ID numbers of the lag variables.
 * @order: number of lags to generate (or 0 for automatic).
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set @order lagged values of the 
 * variables given in @list.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int list_laggenr (int **plist, int order, double ***pZ, DATAINFO *pdinfo)
{
    int *list = *plist;
    int *laglist = NULL;
    int l, i, j;
    int startlen;

    if (order < 0) {
	return E_DATA;
    }

    if (order == 0) {
	order = default_lag_order(pdinfo);
    } 

    laglist = make_lags_list(list, order, pdinfo);
    if (laglist == NULL) {
	return E_ALLOC;
    }

    startlen = get_starting_length(list, pdinfo, (order > 9)? 3 : 2);

    j = 1;
    
    for (i=1; i<=list[0]; i++) {
	int lv, v = list[i];

	if (v == 0 || !pdinfo->vector[v]) {
	    continue;
	}

	for (l=1; l<=order; l++) {
	    lv = get_transform(LAGS, v, l, pZ, pdinfo, startlen);
#if TRDEBUG > 1
	    fprintf(stderr, "base var '%s', lag %d: lv = %d\n",
		    pdinfo->varname[v], l, lv);
#endif
	    if (lv < 0) {
		return 1;
	    }
#if TRDEBUG > 1
	    fprintf(stderr, "lag var name '%s', label '%s'\n",
		    pdinfo->varname[lv], VARLABEL(pdinfo, lv));
#endif
	    laglist[j++] = lv;
	}
    }

    free(*plist);
    *plist = laglist;

    return 0;
}

/**
 * default_lag_order:
 * @pdinfo: data information struct.
 *
 * Returns: the default lag order for generating lags, performing
 * autocorrelation test, and so on.
 *
 */

int default_lag_order (const DATAINFO *pdinfo)
{
    int order = 1;

    if (!dataset_is_panel(pdinfo)) {
	order = (pdinfo->pd < 52)? pdinfo->pd : 14;
    }

    return order;
}

/**
 * list_diffgenr:
 * @list: on entry, list of variables to process; on exit,
 * ID numbers of the generated variables.
 * @ci: must be %DIFF, %LDIFF or %SDIFF.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate differences of the variables in @list, and add them
 * to the data set.  If @ci is %DIFF these are ordinary first
 * differences; if @ci is %LDIFF they are log differences; and
 * if @ci is %SDIFF they are seasonal differences.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int list_diffgenr (int *list, int ci, double ***pZ, DATAINFO *pdinfo)
{
    int i, v, startlen;
    int tnum, err = 0;

    if (ci != DIFF && ci != LDIFF && ci != SDIFF) {
	return 1;
    }

    if (ci == SDIFF && !dataset_is_seasonal(pdinfo)) {
	return E_PDWRONG;
    }    

    startlen = get_starting_length(list, pdinfo, (ci == DIFF)? 2 : 3);
    
    for (i=1; i<=list[0] && !err; i++) {
	v = list[i];
	tnum = get_transform(ci, v, 0, pZ, pdinfo, startlen);
	if (tnum < 0) {
	    err = 1;
	} else {
	    list[i] = tnum;
	}
    }

    return err;
}

/**
 * list_xpxgenr:
 * @plist: pointer to list of variables to process.  On exit
 * the list holds the ID numbers of the squares (and possibly 
 * cross-products).
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: If OPT_O, both squares and cross-products are generated,
 * otherwise only squares.
 *
 * Generates and adds to the data set squares and (if @opt is %OPT_O) 
 * cross-products of the variables given in the list pointed to
 * by @plist.
 *
 * Returns: 0 on success, error code on error.
 */

int list_xpxgenr (int **plist, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt)
{
    int *list = *plist;
    int l0 = list[0];
    int *xpxlist = NULL;
    int tnum, i, j, k, vi, vj;
    int startlen;

    if (opt & OPT_O) {
	int maxterms = (l0 * l0 + l0) / 2;
    
	xpxlist = gretl_list_new(maxterms);
	if (xpxlist == NULL) {
	    return E_ALLOC;
	}
    } else {
	xpxlist = list;
    }

    startlen = get_starting_length(list, pdinfo, 3);
    xpxlist[0] = 0;

    k = 1;
    for (i=1; i<=l0; i++) {
	vi = list[i];

	if (vi == 0 || !pdinfo->vector[vi]) {
	    continue; 
	}
	if (gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[vi])) {
	    continue;
	}

	tnum = get_transform(SQUARE, vi, vi, pZ, pdinfo, startlen);
	if (tnum > 0) {
	    xpxlist[k++] = tnum;
	    xpxlist[0] += 1;
	}

	if (opt & OPT_O) {
	    for (j=i+1; j<=l0; j++) {
		vj = list[j];
		tnum = xpxgenr(vi, vj, pZ, pdinfo);
		if (tnum > 0) {
		    xpxlist[k++] = tnum;
		    xpxlist[0] += 1;
		}
	    }
	}
    }

    if (opt & OPT_O) {
	free(*plist);
	*plist = xpxlist;
    }

    return (xpxlist[0] > 0)? 0 : E_SQUARES;
}

int gettrend (double ***pZ, DATAINFO *pdinfo, int square)
{
    int idx, t, v = pdinfo->v;
    double x;

    idx = varindex(pdinfo, (square)? "timesq" : "time");

    if (idx < v) {
	return idx;
    }
    
    if (dataset_add_series(1, pZ, pdinfo)) {
	return 0; /* error: valid value cannot == 0 */
    }

    for (t=0; t<pdinfo->n; t++) {
	x = (double) t + 1;
	(*pZ)[v][t] = (square)? x * x : x;
    }

    if (square) {
	strcpy(pdinfo->varname[v], "timesq");
	strcpy(VARLABEL(pdinfo, v), _("squared time trend variable"));
    } else {
	strcpy(pdinfo->varname[v], "time");
	strcpy(VARLABEL(pdinfo, v), _("time trend variable"));
    }
	    
    return idx;
}
