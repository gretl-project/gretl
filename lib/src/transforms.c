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

/* newlag is a library global, used when auto-generating a lag in the
   context of reading a regression list (interact.c)
*/
int newlag;

int vars_identical (const double *x, const double *y, int n)
{
    int t;

    for (t=0; t<n; t++) {
	if (floatneq(x[t], y[t])) {
	    return 0;
	}
    }

    return 1;
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
    } else if (cmd == RHODIFF) {
	strncat(vname, orig, len - 1);
	strcat(vname, "#");
    }

    return 0;
}

static void
make_transform_label (char *label, const char *parent,
		      int cmd, int lag)
{
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
    } else if (cmd == RHODIFF) {
	sprintf(label, _("= rho-differenced %s"), parent);
    }
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

    if (pdinfo->time_series == STACKED_CROSS_SECTION) {
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
    if (pdinfo->time_series == STACKED_TIME_SERIES) {
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
	if (pdinfo->time_series == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}
	x0 = Z[v][t];
	if (pdinfo->time_series == STACKED_CROSS_SECTION) {
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
	if (vars_identical(x, (*pZ)[vnum], pdinfo->n)) {
	    /* and it is just what we want */
	    ret = VAR_EXISTS_OK;
	} else {
	    if (!strcmp(label, VARLABEL(pdinfo, vnum))) {
		/* labels match: update the values */
		for (t=0; t<pdinfo->n; t++) {
		    (*pZ)[vnum][t] = x[t];
		}
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
			  double ***pZ, DATAINFO *pdinfo)
{
    char vname[VNAMELEN];
    char label[MAXLABEL];
    int vno, len, err = 0;
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

    if (err) {
	return -1;
    }

    if (ci == SQUARE && v != aux) {
	sprintf(label, _("= %s times %s"), pdinfo->varname[v], 
		pdinfo->varname[aux]);
    } else {
	make_transform_label(label, pdinfo->varname[v], ci, aux);
    }
    
    for (len=8; len<=VNAMELEN; len++) {

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

    lno = get_transform(LOGS, v, lag, pZ, pdinfo);

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
    return get_transform(LOGS, v, 0, pZ, pdinfo);
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

    return get_transform((ldiff)? LDIFF : DIFF, v, 0, pZ, pdinfo);
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

    return get_transform(SQUARE, vi, vj, pZ, pdinfo);
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

int list_loggenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int lognum, i, v;
    int n_ok = 0;

    for (i=1; i<=list[0]; i++) {
	v = list[i];

	if (v == 0 || !pdinfo->vector[v]) {
	    continue; 
	}
	if (isdummy((*pZ)[v], pdinfo->t1, pdinfo->t2)) {
	    continue;
	}

	lognum = loggenr(v, pZ, pdinfo);
	if (lognum > 0) {
	    n_ok++;
	}
    }

    return (n_ok > 0)? 0 : E_LOGS;
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

int list_laggenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int lagnum, l, i, v;
    int maxlag = pdinfo->pd;

    /* play safe with panel data */
    if (dataset_is_panel(pdinfo)) {
	maxlag = 1;
    }
    
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v == 0 || !pdinfo->vector[v]) {
	    continue;
	}
	for (l=1; l<=maxlag; l++) {
	    lagnum = laggenr(v, l, pZ, pdinfo);
	    if (lagnum < 0) {
		return 1;
	    }
	}
    }

    return 0;
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

int list_diffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    
    for (i=1; i<=list[0]; i++) {
	if (diffgenr(list[i], pZ, pdinfo, 0) < 0) {
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

int list_ldiffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    
    for (i=1; i<=list[0]; i++) {
	if (diffgenr(list[i], pZ, pdinfo, 1) < 0) {
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

int list_xpxgenr (const LIST list, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt)
{
    int xnum, i, j, vi, vj;
    int n_ok = 0;

    for (i=1; i<=list[0]; i++) {
	vi = list[i];

	if (vi == 0 || !pdinfo->vector[vi]) {
	    continue; 
	}
	if (isdummy((*pZ)[vi], pdinfo->t1, pdinfo->t2)) {
	    continue;
	}

	xnum = xpxgenr(vi, vi, pZ, pdinfo);
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

#undef RHODEBUG

/* rhodiff: a "legacy" function */

/**
 * rhodiff:
 * @param: please see the gretl help on rhodiff() for syntax.
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set rho-differenced versions
 * of the variables given in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int rhodiff (char *param, const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i, j, maxlag, p, t, t1, nv;
    int v = pdinfo->v, n = pdinfo->n;
    char parmbit[VNAMELEN];
    double *rhot;

#ifdef RHODEBUG
    fprintf(stderr, "rhodiff: param = '%s'\n", param);
#endif

    maxlag = count_fields(param);
    rhot = malloc(maxlag * sizeof *rhot);
    if (rhot == NULL) {
	return E_ALLOC;
    }

    if (maxlag > pdinfo->t1) {
	t1 = maxlag;
    } else {
	t1 = pdinfo->t1;
    }

#ifdef RHODEBUG
    fprintf(stderr, "rhodiff: maxlag = %d, t1 = %d\n", maxlag, t1);
#endif

    /* parse "param" string */
    j = strlen(param);
    p = 0;
    for (i=0; i<j; i++) {
	if ((i == 0 || param[i] == ' ') && i < (j - 1)) {
	    sscanf(param + i + (i? 1 : 0), "%8s", parmbit); 
#ifdef RHODEBUG
	    fprintf(stderr, "rhodiff: parmbit = '%s'\n", parmbit);
#endif
	    if (isalpha((unsigned char) parmbit[0])) {
		nv = varindex(pdinfo, parmbit);
		if (nv == v) {
		    free(rhot);
		    return E_UNKVAR;
		}
		rhot[p] = get_xvalue(nv, (const double **) *pZ, pdinfo);
	    } else {
		rhot[p] = dot_atof(parmbit);
	    }
	    p++;
	}
    }

    if (dataset_add_vars(list[0], pZ, pdinfo)) {
	return E_ALLOC;
    }

    for (i=1; i<=list[0]; i++) {
	int vr = v + i - 1;
	double xx;

	j = list[i];

	make_transform_varname(pdinfo->varname[vr], 
			       pdinfo->varname[j], 
			       RHODIFF, 0, 8);
	make_transform_label(VARLABEL(pdinfo, vr), pdinfo->varname[j],
			     RHODIFF, 0);

	for (t=0; t<n; t++) {
	    if (t < t1 || t > pdinfo->t2) {
		(*pZ)[vr][t] = NADBL;
		continue;
	    }
	    xx = (*pZ)[j][t];
	    if (na(xx)) {
		(*pZ)[vr][t] = NADBL;
		continue;
	    }
	    for (p=0; p<maxlag; p++) {
		if (na((*pZ)[j][t-p-1])) {
		    xx = NADBL;
		    break;
		} else {
		    xx -= rhot[p] * (*pZ)[j][t-p-1];
		}
	    }
	    (*pZ)[vr][t] = xx;
	}
    }

    free(rhot);

    return 0;
}

