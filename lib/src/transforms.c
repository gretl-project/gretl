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

#include "libgretl.h"

#define TRDEBUG 0

enum {
    VARS_IDENTICAL,
    X_HAS_MISSING,
    Y_HAS_MISSING,
    VARS_DIFFER
} varcomp_codes;

enum {
    INVERSE = NC + 1
};

static char *get_mangled_name_by_id (int v);

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
make_transform_varname (char *vname, const char *orig, int ci, 
			int aux, int len)
{
    *vname = '\0';


    if (ci == DIFF) {
	strcpy(vname, "d_");
	strncat(vname, orig, len - 2);
    } else if (ci == LDIFF) {
	strcpy(vname, "ld_");
	strncat(vname, orig, len - 3);
    } else if (ci == SDIFF) {
	strcpy(vname, "sd_");
	strncat(vname, orig, len - 3);
    } else if (ci == ORTHDEV) {
	strcpy(vname, "o_");
	strncat(vname, orig, len - 2);
    } else if (ci == LOGS) {
	strcpy(vname, "l_");
	strncat(vname, orig, len - 2);
    } else if (ci == SQUARE) {
	strcpy(vname, "sq_");
	strncat(vname, orig, len - 3);
    } else if (ci == LAGS) {
	char ext[6];

	if (aux >= 0) {
	    /* an actual lag */
	    sprintf(ext, "_%d", aux);
	} else {
	    /* in fact a lead */
	    sprintf(ext, "%d", -aux);
	}
	strncat(vname, orig, len - strlen(ext));
	strcat(vname, ext);
    } else if (ci == DUMMIFY) {
	char ext[6];

	sprintf(ext, "_%d", aux);
	strcpy(vname, "D");
	strncat(vname, orig, len - strlen(ext) - 1);
	strcat(vname, ext);
    } else if (ci == INVERSE) {
	strcpy(vname, "i_");
	strncat(vname, orig, len - 2);
    }

#if TRDEBUG
    fprintf(stderr, "make_transform_varname:\n"
	    "orig='%s', ci=%d, len=%d, vname='%s'\n", 
	    orig, ci, len, vname);
#endif

    return 0;
}

static int
make_transform_label (char *label, const char *parent,
		      int ci, int lag)
{
    int err = 0;

    if (ci == DIFF) {
	sprintf(label, _("= first difference of %s"), parent);
    } else if (ci == LDIFF) {
	sprintf(label, _("= log difference of %s"), parent);
    } else if (ci == SDIFF) {
	sprintf(label, _("= seasonal difference of %s"), parent);
    } else if (ci == LOGS) {
	sprintf(label, _("= log of %s"), parent);
    } else if (ci == SQUARE) {
	sprintf(label, _("= %s squared"), parent);
    } else if (ci == LAGS) {
	if (lag >= 0) {
	    sprintf(label, "= %s(t - %d)", parent, lag);
	} else {
	    sprintf(label, "= %s(t + %d)", parent, -lag);
	}
    } else if (ci == INVERSE) {
	sprintf(label, "= 1/%s", parent);
    } else {
	err = 1;
    }

    return err;
}

int is_standard_lag (int v, const DATAINFO *pdinfo, int *parent)
{
    int pv = 0, ret = 0;

    if (pdinfo == NULL || v <= 0 || v >= pdinfo->v) {
	return 0;
    }

    if (pdinfo->varinfo[v]->transform == LAGS) {
	pv = series_index(pdinfo, pdinfo->varinfo[v]->parent);
	pv = (pv < pdinfo->v)? pv : 0;
	ret = pdinfo->varinfo[v]->lag;
    }

    if (pv > 0) {
	int n = strcspn(pdinfo->varname[v], "_");

	if (!strncmp(pdinfo->varname[v], pdinfo->varname[pv], n)) {
	    if (parent != NULL) {
		*parent = pv;
	    }
	} else {
	    /* if the lag has a "custom" name, don't consider it
	       a child to parent variable pv */
	    ret = 0;
	}
    }

    return ret;
}

int is_standard_diff (int v, const DATAINFO *pdinfo, int *parent)
{
    int pv = 0, ret = 0;

    if (v <= 0 || v >= pdinfo->v) {
	return 0;
    }

    if (pdinfo->varinfo[v]->transform == DIFF) {
	pv = series_index(pdinfo, pdinfo->varinfo[v]->parent);
	pv = (pv < pdinfo->v)? pv : 0;
	if (pv > 0) {
	    if (parent != NULL) {
		*parent = pv;
	    }
	    ret = 1;
	}
    }

    return ret;
}

int is_dummy_child (int v, const DATAINFO *pdinfo, int *parent)
{
    int pv = pdinfo->v;
    int i = 0, ret = 0;

    if (pdinfo->varinfo[v]->transform == DUMMIFY) {
	pv = series_index(pdinfo, pdinfo->varinfo[v]->parent);
    } else if (!strncmp(pdinfo->varname[v], "dt_", 3)) {
	if (sscanf(pdinfo->varname[v] + 3, "%d", &i) && i > 1) {
	    pv = series_index(pdinfo, "dt_1");
	}
    } else if (!strncmp(pdinfo->varname[v], "du_", 3)) {
	if (sscanf(pdinfo->varname[v] + 3, "%d", &i) && i > 1) {
	    pv = series_index(pdinfo, "du_1");
	}
    }	

    if (pv < pdinfo->v) {
	*parent = pv;
	ret = 1;
    } else {
	*parent = 0;
    }

    return ret;
}

static void make_xp_varname (char *vname, int v1, int v2,
			     const DATAINFO *pdinfo,
			     int len)
{
    int n1 = strlen(pdinfo->varname[v1]);
    int n2 = strlen(pdinfo->varname[v2]);
    int cut = n1 + n2 + 1 - len;
    int cut1 = 0, cut2 = 0;

    if (cut > 0) {
	cut1 = cut2 = cut/2;
	if (cut % 2) {
	    cut2++;
	}
    }

    *vname = '\0';
    strncat(vname, pdinfo->varname[v1], n1 - cut1);
    strcat(vname, "_");
    strncat(vname, pdinfo->varname[v2], n2 - cut2);    
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
    int t1 = (lag > 0)? lag : 0;
    int t2 = pdinfo->n - 1;
    int t, lt, miss;

    for (t=0; t<pdinfo->n; t++) {
	lagvec[t] = NADBL;
    }

    if (dated_daily_data(pdinfo)) {
	for (t=t1; t<=t2; t++) {
	    lt = t - lag;
	    miss = 0;
	    while (na(Z[v][lt]) && lt > 0 && miss < 6) {
		lt--;
		miss++;
	    }
	    lagvec[t] = Z[v][lt];
	}
    } else { 
	/* the "standard" time-series case */
	for (t=t1; t<=t2; t++) {
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

	for (t=t1; t<=t2; t++) {
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
	xx = Z[v][t];
	if (na(xx) || xx <= 0.0) {
	    logvec[t] = NADBL;
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
	x1 = Z[v][t - t0];

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

/* orthogonal deviations */

static int get_orthdev (int v, double *xvec, const double **Z, 
			const DATAINFO *pdinfo)
{
    return orthdev_series(Z[v], xvec, pdinfo);
}

/* write square or cross-product into xvec */

static int get_xpx (int vi, int vj, double *xvec, const double **Z, 
		    const DATAINFO *pdinfo)
{
    double xit, xjt;
    int t;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xit = Z[vi][t];
	xjt = Z[vj][t];
	if (na(xit) || na(xjt)) {
	    xvec[t] = NADBL;
	} else {
	    xvec[t] = xit * xjt;
	}
    }

    return 0;
}

/* write reciprocal into xvec */

static int get_inverse (int v, double *xvec, const double **Z, 
			const DATAINFO *pdinfo)
{
    double xt;
    int t;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xt = Z[v][t];
	if (na(xt) || xt == 0.0) {
	    xvec[t] = NADBL;
	} else {
	    xvec[t] = 1.0 / xt;
	}
    }

    return 0;
}

/* write dummy for (v == value) into xvec */

static int get_discdum (int v, double val, double *xvec, const double **Z, 
			const DATAINFO *pdinfo)
{
    double xt;
    int t;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	xt = Z[v][t];
	if (na(xt)) {
	    xvec[t] = NADBL;
	} else {
	    xvec[t] = (xt == val)? 1.0 : 0.0;
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

/* Note: new behavior as of 2008-02-11 */

#define TR_OVERWRITE 1

static int transform_handle_duplicate (int ci, int lag, int v, 
				       const double *x, const char *label,
				       DATAINFO *pdinfo, double **Z,
				       int origv)
{
    int ret = VARNAME_DUPLICATE;
    int t, ok = 0;

    if (!strcmp(label, VARLABEL(pdinfo, v))) {
	/* labels identical, so OK? */
	ok = 1;
	
    }

#if TR_OVERWRITE
    if (!ok && v < origv) {
	ok = 1;
    }
#endif

    if (ok) {
#if TRDEBUG
	fprintf(stderr, "transform_handle_duplicate: updating var %d (%s)\n",
		v, pdinfo->varname[v]);
#endif
	for (t=0; t<pdinfo->n; t++) {
	    Z[v][t] = x[t];
	}
	if (*label != '\0') {
	    strcpy(VARLABEL(pdinfo, v), label);
	}
	pdinfo->varinfo[v]->transform = ci;
	pdinfo->varinfo[v]->lag = lag;
	pdinfo->varinfo[v]->flags = 0;
	ret = VAR_EXISTS_OK;
    }	

    return ret;
}

static int 
check_add_transform (int ci, int lag, int vnum, const double *x,
		     const char *vname, const char *label,
		     DATAINFO *pdinfo, double ***pZ,
		     int origv)
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
	} else if (chk == X_HAS_MISSING) { 
	    /* is this right? */
	    ret = VAR_EXISTS_OK;
	} else if (chk == Y_HAS_MISSING) {
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[vnum][t] = x[t];
	    }
	    ret = VAR_EXISTS_OK;
	} else {
	    ret = transform_handle_duplicate(ci, lag, vnum, x, label, pdinfo,
					     *pZ, origv);
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

static int get_lag_ID (int srcv, int lag, const DATAINFO *pdinfo)
{
    const char *vname = pdinfo->varname[srcv];
    VARINFO *vinfo;
    int i;

    for (i=1; i<pdinfo->v; i++) {
	vinfo = pdinfo->varinfo[i];
	if (vinfo->lag == lag && !strcmp(vinfo->parent, vname)) {
	    return i;
	}
    }

    return -1;
}

/* get_transform: create specified transformation of variable v if
   this variable does not already exist.

   The dummies resulting from DUMMIFY are automatically marked as
   discrete.

   origv is the number of variables in the dataset prior to the
   current round of adding.

   Return the ID number of the transformed var, or -1 on error.
*/

static int get_transform (int ci, int v, int aux, double x, 
			  double ***pZ, DATAINFO *pdinfo,
			  int startlen, int origv)
{
    char vname[VNAMELEN] = {0};
    char label[MAXLABEL] = {0};
    int vno = -1, err = 0;
    int len, lag = 0;
    const char *srcname;
    double *vx;

    vx = testvec(pdinfo->n);
    if (vx == NULL) {
	return -1;
    }

    if (ci == LAGS) {
	lag = aux;
	err = get_lag(v, lag, vx, (const double **) *pZ, pdinfo);
    } else if (ci == LOGS) {
	err = get_log(v, vx, (const double **) *pZ, pdinfo);
    } else if (ci == DIFF || ci == LDIFF || ci == SDIFF) {
	err = get_diff(v, vx, ci, (const double **) *pZ, pdinfo);
    } else if (ci == ORTHDEV) {
	err = get_orthdev(v, vx, (const double **) *pZ, pdinfo);
    } else if (ci == SQUARE) {
	/* "aux" = second variable number */
	err = get_xpx(v, aux, vx, (const double **) *pZ, pdinfo);
    } else if (ci == DUMMIFY) {
	/* "x" = value for dummy */
	err = get_discdum(v, x, vx, (const double **) *pZ, pdinfo);
    } else if (ci == INVERSE) {
	err = get_inverse(v, vx, (const double **) *pZ, pdinfo);
    }

    if (err) {
	return -1;
    }

    if (ci == LAGS && (vno = get_lag_ID(v, aux, pdinfo)) > 0) {
	/* special case: pre-existing lag */
	err = check_add_transform(ci, lag, vno, vx, vname, label, pdinfo, pZ, origv);
	if (err != VAR_EXISTS_OK) {
	    vno = -1;
	}
	return vno;
    }

    if (ci == SQUARE && v != aux) {
	sprintf(label, _("= %s times %s"), pdinfo->varname[v], 
		pdinfo->varname[aux]);
    } else if (ci == DUMMIFY) {
	sprintf(label, _("dummy for %s = %g"), pdinfo->varname[v], x);
    } else {
	make_transform_label(label, pdinfo->varname[v], ci, aux);
    }

    srcname = get_mangled_name_by_id(v);
    if (srcname == NULL) {
	srcname = pdinfo->varname[v];
    }

    for (len=startlen; len<=VNAMELEN; len++) {
	if (len == VNAMELEN) {
	    /* last resort: hack the name */
	    make_varname_unique(vname, 0, pdinfo);
	} else if (ci == SQUARE && v != aux) {
	    make_xp_varname(vname, v, aux, pdinfo, len);
	} else {
	    make_transform_varname(vname, srcname, ci, aux, len);
	}
	vno = series_index(pdinfo, vname);

	err = check_add_transform(ci, lag, vno, vx, vname, label, pdinfo, pZ, origv);
	if (err != VAR_EXISTS_OK && err != VAR_ADDED_OK) {
	    vno = -1;
	}
	if (err != VARNAME_DUPLICATE) {
	    /* either we're OK, or there was a fatal error */
	    break;
	}
    }

    if (!err && vno > 0) {
	if (ci == DUMMIFY) {
	    strcpy(pdinfo->varinfo[vno]->parent, pdinfo->varname[v]);
	    pdinfo->varinfo[vno]->transform = DUMMIFY;
	    set_var_discrete(pdinfo, vno, 1);
	} else if (ci == LAGS || ci == DIFF) {
	    strcpy(pdinfo->varinfo[vno]->parent, pdinfo->varname[v]);
	    pdinfo->varinfo[vno]->transform = ci;
	}

	if (ci == LAGS) {
	    pdinfo->varinfo[vno]->lag = aux;
	} else {
	    pdinfo->varinfo[vno]->lag = 0;
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

    if (lag > pdinfo->n || -lag > pdinfo->n) {
	gretl_errmsg_sprintf(_("Invalid lag order %d"), lag);
	lno = -1;
    } else if (lag == 0) {
	lno = v;
    } else {
	lno = get_transform(LAGS, v, lag, 0.0, pZ, pdinfo, 
			    VNAMELEN - 3, pdinfo->v);
    }

    return lno;
}

/**
 * laggenr_from_to: 
 * @v: ID number in dataset of source variable.
 * @minlag: minimum lag order.
 * @maxlag: maximum lag order.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 * @err: location to receive error code.
 *
 * Creates the specified lags of variable @v if they do not
 * already exist.
 *
 * Returns: list of lag variables, or %NULL or on error.
 */

int *laggenr_from_to (int v, int minlag, int maxlag, double ***pZ, 
		      DATAINFO *pdinfo, int *err)
{
    int *llist;
    int i, lv, nlags = -1;
    int p, reverse;

    if (maxlag < 0) {
	nlags = -maxlag + minlag + 1;
    } else if (maxlag > 0) {
	nlags = maxlag - minlag + 1;
    } 

    if (nlags <= 0) {
	*err = E_DATA;
	return NULL;
    }

    llist = gretl_list_new(nlags);
    if (llist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    reverse = (maxlag < minlag);
    p = minlag;

    for (i=0; i<nlags; i++) {
	lv = laggenr(v, p, pZ, pdinfo);
	if (lv < 0) {
	    *err = E_DATA;
	    free(llist);
	    llist = NULL;
	    break;
	}
	llist[i+1] = lv;
	if (reverse) {
	    p--;
	} else {
	    p++;
	}
    }

    return llist;
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
    return get_transform(LOGS, v, 0, 0.0, pZ, pdinfo, 
			 VNAMELEN - 3, pdinfo->v);
}

/**
 * invgenr: 
 * @v: ID number in dataset of source variable.
 * @pZ: pointer to data array.
 * @pdinfo: information on dataset.
 *
 * Creates the reciprocal of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the reciprocal, or -1 on error.
 */

int invgenr (int v, double ***pZ, DATAINFO *pdinfo)
{
    return get_transform(INVERSE, v, 0, 0.0, pZ, pdinfo, 
			 VNAMELEN - 3, pdinfo->v);
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
    if (ci != DIFF && ci != LDIFF && ci != SDIFF) {
	return -1;
    }

    if (ci == SDIFF && !dataset_is_seasonal(pdinfo)) {
	return -1;
    }

    return get_transform(ci, v, 0, 0.0, pZ, pdinfo, 
			 VNAMELEN - 3, pdinfo->v);
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

    return get_transform(SQUARE, vi, vj, 0.0, pZ, pdinfo, 
			 VNAMELEN - 3, pdinfo->v);
}

static int 
get_starting_length (const int *list, DATAINFO *pdinfo, int trim)
{
    int width = VNAMELEN - 3 - trim;
    const char *vni, *vnj;
    int len, maxlen = 0;
    int conflict = 0;
    int i, j;

    if (list[0] == 1) {
	return VNAMELEN - 1;
    }

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
	conflict = 0;
	for (i=1; i<=list[0] && !conflict; i++) {
	    vni = pdinfo->varname[list[i]];
	    for (j=i+1; j<=list[0] && !conflict; j++) {
		vnj = pdinfo->varname[list[j]];
		if (!strncmp(vni, vnj, len)) {
		    conflict = 1;
		}
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

struct mangled_name {
    int v;
    char *s;
};

static struct mangled_name *mnames;
static int n_mnames;

static void destroy_mangled_names (void)
{
    int i;

    for (i=0; i<n_mnames; i++) {
	free(mnames[i].s);
    }
    free(mnames);
    mnames = NULL;
    n_mnames = 0;
}

static int push_mangled_name (int v, const char *s)
{
    static struct mangled_name *tmp;

    tmp = realloc(mnames, (n_mnames + 1) * sizeof *tmp);

    if (tmp == NULL) {
	destroy_mangled_names();
	return E_ALLOC;
    }

    mnames = tmp;
    mnames[n_mnames].v = v;
    mnames[n_mnames].s = gretl_strdup(s);

    n_mnames++;

    return 0;
}

static char *get_mangled_name_by_id (int v)
{
    int i;

    for (i=0; i<n_mnames; i++) {
	if (mnames[i].v == v) {
	    return mnames[i].s;
	}
    }

    return NULL;
}

static int make_mangled_name (int v, const char *s, int nc)
{
    char tmp[16], sfx[16];
    int n, err;

    if (get_mangled_name_by_id(v)) {
	return 0;
    }

    *tmp = *sfx = '\0';
    /* we know the bit out here must be distinct? */
    strcat(sfx, s + nc);
    n = strlen(sfx);
    strncat(tmp, s, nc - n);
    strcat(tmp, sfx);

#if TRDEBUG
    fprintf(stderr, "mangled name: '%s' -> '%s'\n", s, tmp);
#endif

    err = push_mangled_name(v, tmp);

    return err;
}

static int 
transform_preprocess_list (int *list, double **Z, const DATAINFO *pdinfo,
			   int f)
{
    int maxc = VNAMELEN - 3;
    int longnames = 0;
    char **S = NULL;
    int i, v, ok;
    int err = 0;

    if (f == SQUARE || f == LDIFF || f == SDIFF) {
	/* 3-character prefixes */
	maxc--;
    } else if (f == DUMMIFY) {
	/* "D" plus suffix */
	maxc--;
    }

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	ok = 1;

	if (f == LOGS || f == SQUARE) {
	    if (v == 0) { /* FIXME?? */
		ok = 1; 
	    }
	    if (gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[v])) {
		ok = 0;
	    }
	} else if (f == LAGS) {
	    if (v == 0) {
		ok = 0;
	    }
	} else if (f == DUMMIFY) {
	    ok = 0; /* reverse burden of proof */
	    if (v > 0) {
		if (var_is_discrete(pdinfo, v)) {
		    /* pre-approved */
		    ok = 1;
		} else if (gretl_isdiscrete(0, pdinfo->n - 1, Z[v])) {
		    ok = 1;
		}
	    }
	}

	if (!ok) {
	    gretl_list_delete_at_pos(list, i--);
	} else if (strlen(pdinfo->varname[v]) > maxc) {
	    strings_array_add(&S, &longnames, pdinfo->varname[v]);
	}
    }

    if (list[0] == 0) {
	err = E_TYPES;
    }

    if (longnames > 0) {
	if (!err) {
	    int j, herr = 0;

	    for (i=0; i<longnames && !herr; i++) {
		for (j=i+1; j<longnames && !herr; j++) {
		    if (!strncmp(S[i], S[j], maxc)) {
			v = series_index(pdinfo, S[i]);
			herr = make_mangled_name(v, S[i], maxc);
			v = series_index(pdinfo, S[j]);
			herr += make_mangled_name(v, S[j], maxc);
		    }
		}
	    }
	}

	free_strings_array(S, longnames);
    }

    return err;
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
    int origv = pdinfo->v;
    int tnum, i, j, v;
    int startlen;
    int l0 = 0;
    int err;

    err = transform_preprocess_list(list, *pZ, pdinfo, LOGS);
    if (err) {
	return err;
    }

    startlen = get_starting_length(list, pdinfo, 2);

    j = 1;
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	tnum = get_transform(LOGS, v, 0, 0.0, pZ, pdinfo, startlen,
			     origv);
	if (tnum > 0) {
	    list[j++] = tnum;
	    l0++;
	}
    }

    list[0] = l0;

    destroy_mangled_names();

    return (l0 > 0)? 0 : E_LOGS;
}

static int *make_lags_list (int *list, int order, DATAINFO *pdinfo)
{
    int i, v, nl = 0;

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0) {
	    nl += order;
	}
    }

    return gretl_list_new(nl);
}

/**
 * list_laggenr:
 * @plist: on entry, pointer to list of variables to process.  On exit
 * the list holds the ID numbers of the lag variables.
 * @order: number of lags to generate (or 0 for automatic).
 * @pZ: pointer to data array.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set @order lagged values of the 
 * variables given in the list pointed to by @plist.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int list_laggenr (int **plist, int order, double ***pZ, DATAINFO *pdinfo)
{
    int origv = pdinfo->v;
    int *list = *plist;
    int *laglist = NULL;
    int l, i, j, v, lv;
    int startlen, l0 = 0;
    int err;

    if (order < 0) {
	gretl_errmsg_sprintf(_("Invalid lag order %d"), order);
	return E_DATA;
    }

    if (order == 0) {
	order = default_lag_order(pdinfo);
    } 

    err = transform_preprocess_list(list, *pZ, pdinfo, LAGS);
    if (err) {
	return err;
    }

    laglist = make_lags_list(list, order, pdinfo);
    if (laglist == NULL) {
	destroy_mangled_names();
	return E_ALLOC;
    }

    startlen = get_starting_length(list, pdinfo, (order > 9)? 3 : 2);

    j = 1;
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	for (l=1; l<=order; l++) {
	    lv = get_transform(LAGS, v, l, 0.0, pZ, pdinfo, startlen, origv);
#if TRDEBUG > 1
	    fprintf(stderr, "base var '%s', lag %d: lv = %d\n",
		    pdinfo->varname[v], l, lv);
#endif
	    if (lv > 0) {
#if TRDEBUG > 1
		fprintf(stderr, "lag var name '%s', label '%s'\n",
			pdinfo->varname[lv], VARLABEL(pdinfo, lv));
#endif
		laglist[j++] = lv;
		l0++;
	    }
	}
    }

    destroy_mangled_names();

    laglist[0] = l0;

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
    int origv = pdinfo->v;
    int i, v, startlen;
    int tnum, l0 = 0;
    int err;

    if (list[0] == 0) {
	return 0;
    }

    if (ci != DIFF && ci != LDIFF && ci != SDIFF) {
	return 1;
    }

    if (ci == SDIFF && !dataset_is_seasonal(pdinfo)) {
	return E_PDWRONG;
    } 

    err = transform_preprocess_list(list, *pZ, pdinfo, ci);
    if (err) {
	return err;
    }

    startlen = get_starting_length(list, pdinfo, (ci == DIFF)? 2 : 3);
    
    for (i=1; i<=list[0] && !err; i++) {
	v = list[i];
	tnum = get_transform(ci, v, 0, 0.0, pZ, pdinfo, startlen, origv);
	if (tnum < 0) {
	    err = 1;
	} else {
	    list[i] = tnum;
	    l0++;
	}
    }

    list[0] = l0;

    destroy_mangled_names();

    return err;
}

/**
 * list_orthdev:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate orthogonal deviations of the variables in @list, and add
 * them to the data set.
 *
 * Returns: 0 on success, error code on error.
 */

int list_orthdev (int *list, double ***pZ, DATAINFO *pdinfo)
{
    int origv = pdinfo->v;
    int i, v, startlen;
    int tnum, l0 = 0;
    int err;

    if (list[0] == 0) {
	return 0;
    }

    if (!dataset_is_panel(pdinfo)) {
	return E_PDWRONG;
    } 

    err = transform_preprocess_list(list, *pZ, pdinfo, ORTHDEV);
    if (err) {
	return err;
    }

    startlen = get_starting_length(list, pdinfo, 2);
    
    for (i=1; i<=list[0] && !err; i++) {
	v = list[i];
	tnum = get_transform(ORTHDEV, v, 0, 0.0, pZ, pdinfo, startlen, origv);
	if (tnum < 0) {
	    err = 1;
	} else {
	    list[i] = tnum;
	    l0++;
	}
    }

    list[0] = l0;

    destroy_mangled_names();

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
    int origv = pdinfo->v;
    int *list = *plist;
    int *xpxlist = NULL;
    int tnum, i, j, k, vi, vj;
    int startlen, l0;
    int err;

    err = transform_preprocess_list(list, *pZ, pdinfo, SQUARE);
    if (err) {
	return err;
    }

    l0 = list[0];

    if (opt & OPT_O) {
	int maxterms = (l0 * l0 + l0) / 2;
    
	xpxlist = gretl_list_new(maxterms);
	if (xpxlist == NULL) {
	    destroy_mangled_names();
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
	tnum = get_transform(SQUARE, vi, vi, 0.0, pZ, pdinfo, startlen,
			     origv);
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

    destroy_mangled_names();

    if (opt & OPT_O) {
	free(*plist);
	*plist = xpxlist;
    }

    return (xpxlist[0] > 0)? 0 : E_SQUARES;
}

#define DUMDEBUG 0

#define skip_j (x, xs) (!na(xs) && (xs == x))

static int real_list_dumgenr (int **plist, double ***pZ, DATAINFO *pdinfo,
			      double oddval, gretlopt opt)
{
    int origv = pdinfo->v;
    int *list = *plist;
    int *tmplist = NULL;
    double *x = NULL;
    int i, j, t, n;
    int startlen;
    int err;

    err = transform_preprocess_list(list, *pZ, pdinfo, DUMMIFY);
    if (err) {
	return err;
    }

    tmplist = gretl_null_list();
    if (tmplist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    
    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    startlen = get_starting_length(list, pdinfo, 3);

    for (i=1; i<=list[0] && !err; i++) {
	int vi = list[i];
	int nuniq, tnum;
	int jmin, jmax;
	double xt;

	n = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    xt = (*pZ)[vi][t];
	    if (!na(xt)) {
		x[n++] = xt;
	    }
	}

	if (n == 0) {
	    continue;
	}

	qsort(x, n, sizeof *x, gretl_compare_doubles);
	nuniq = count_distinct_values(x, n);

	if (nuniq == 1) {
	    continue;
	}

	rearrange_id_array(x, nuniq, n);

	jmin = (opt & OPT_F)? 1 : 0;
	jmax = (opt & OPT_L)? nuniq - 1 : nuniq;

#if DUMDEBUG 
	fprintf(stderr, "variable %d has %d distinct values\n", vi, nuniq);
	if (opt & OPT_F) fprintf(stderr, "skipping lowest value\n");
	if (opt & OPT_L) fprintf(stderr, "skipping highest value\n");
	fprintf(stderr, "jskip = %d\n", jskip);
#endif

	for (j=jmin; j<jmax && !err; j++) {
	    if (x[j] != oddval) {
		tnum = get_transform(DUMMIFY, vi, j+1, x[j], pZ, pdinfo, 
				     startlen, origv);
#if DUMDEBUG   
		fprintf(stderr, "VALUE = %g, tnum = %d\n", x[j], tnum);
#endif
		if (tnum > 0) {
		    tmplist = gretl_list_append_term(&tmplist, tnum);
		    if (tmplist == NULL) {
			err = E_ALLOC;
		    } 
		} else {
		    err = E_DATA;
		}
	    }
	}
    }

    if (!err && tmplist[0] == 0) {
	gretl_errmsg_set(_("dummify: no suitable variables were found"));
	err = E_DATA;
    }

    free(x);

 bailout:

    if (!err) {
	free(*plist);
	*plist = tmplist;
#if DUMDEBUG   
	printlist(tmplist, "output list");
#endif
    } else {
	free(tmplist);
    }

    destroy_mangled_names();

    return err;
}

/**
 * list_dumgenr:
 * @plist: pointer to list of variables to process; on exit
 * the list holds the ID numbers of the generated dummies.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: can be %OPT_F to drop the first value, %OPT_L to drop
 * the last, or %OPT_NONE.
 *
 * For each of the variables given in the list to which @plist
 * points, generates and adds to the data set k dummy variables 
 * coding for the k distinct values of the variable in question.
 * All these variables must have already been marked as discrete.
 * If the %OPT_F or %OPT_L option is given, either the first or
 * the last value of each variable is taken as the "base" and is
 * not given a dummy encoding (that is, only k - 1 dummies are
 * added for each variable).
 *
 * Returns: 0 on success, error code on error.
*/

int list_dumgenr (int **plist, double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt)
{
    return real_list_dumgenr(plist, pZ, pdinfo, NADBL, opt);
}

/**
 * dumgenr_with_oddval:
 * @plist: pointer to list of variables to process; on exit
 * the list holds the ID numbers of the generated dummies.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @oddval: value which should be skipped when encoding the
 * input values as dummies.
 *
 * For each of the variables given in the list to which @plist
 * points, generates and adds to the data set k dummy variables 
 * coding for the k distinct values of the variable in question.
 * All these variables must have already been marked as discrete.
 * if @oddval is not #NADBL, it is treated as the omitted
 * category and only k - 1 dummies are added for each variable).
 *
 * Returns: 0 on success, error code on error.
*/

int dumgenr_with_oddval (int **plist, double ***pZ, DATAINFO *pdinfo,
			 double oddval)
{
    return real_list_dumgenr(plist, pZ, pdinfo, oddval, OPT_NONE);
}

/**
 * list_makediscrete:
 * @list: list of variables to process.
 * @pdinfo: data information struct.
 * @opt: if %OPT_R, reverse the operation.
 *
 * Sets the variables given in @list as discrete, unless
 * opt is %OPT_R, in which case the variables are set as
 * continuous.
 *
 * Returns: 0 on success, error code on error.
 */

int list_makediscrete (const int *list, DATAINFO *pdinfo, gretlopt opt)
{
    int disc = !(opt & OPT_R);
    int i, v, err = 0;

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0) {
	    set_var_discrete(pdinfo, v, disc);
	}
    }

    return err;
}

int gettrend (double ***pZ, DATAINFO *pdinfo, int square)
{
    int idx, t, v = pdinfo->v;
    double x;

    idx = series_index(pdinfo, (square)? "timesq" : "time");

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


