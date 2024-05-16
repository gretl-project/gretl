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
#include "uservar.h"
#include "gretl_func.h"
#include "matrix_extra.h"

#define TRDEBUG 0

/**
 * SECTION:transforms
 * @short_description: standard transformations of series in the dataset
 * @title: Transformations
 * @include: libgretl.h
 *
 * Functions to generate standard transformations (logs, lags, first
 * differences and so on) of series in the dataset.
 */

enum {
    VARS_IDENTICAL,
    NEW_HAS_MISSING,
    OLD_HAS_MISSING,
    VARS_DIFFER
} varcomp_codes;

enum {
    INVERSE = NC + 1,
    RESAMPLE,
    HFDIFF,
    HFLDIFF
};

static char *get_mangled_name_by_id (int v);

/* Note: in this comparison @test is the newly generated series
   and @cand is the previously existing candidate series.
*/

static int check_vals (const double *test, const double *cand, int n)
{
    int t, ret = VARS_IDENTICAL;

    for (t=0; t<n; t++) {
	if (test[t] == cand[t] || (na(test[t]) && na(cand[t]))) {
	    /* values the same, or both NA */
	    continue;
	}
	/* values differ at obs @t */
	if (na(test[t]) && !na(cand[t])) {
	    if (ret == VARS_IDENTICAL || ret == NEW_HAS_MISSING) {
		ret = NEW_HAS_MISSING;
	    } else {
		ret = VARS_DIFFER;
	    }
	} else if (!na(test[t]) && na(cand[t])) {
	    if (ret == VARS_IDENTICAL || ret == OLD_HAS_MISSING) {
		ret = OLD_HAS_MISSING;
	    } else {
		ret = VARS_DIFFER;
	    }
	} else {
	    ret = VARS_DIFFER;
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
	char ext[12];

	if (aux >= 0) {
	    /* an actual lag (or contemporaneous value) */
	    sprintf(ext, "_%d", aux);
	    strncat(vname, orig, len - strlen(ext));
	    strcat(vname, ext);
	} else {
	    /* a lead: in naming the output series we'll
	       try to avoid appending a plain digit to
	       a series name that already ends in a digit,
	       since this may get confusing
	    */
	    int n1, n2;

	    sprintf(ext, "%d", -aux);
	    n1 = strlen(orig);
	    n2 = strlen(ext);
	    if (n1 + n2 + 2 < VNAMELEN) {
		if (isdigit(orig[n1-1])) {
		    /* separate the digits */
		    sprintf(vname, "%s_f%d", orig, -aux);
		} else {
		    sprintf(vname, "%s%d", orig, -aux);
		}
	    } else {
		/* Oh, well, can't be helped */
		strncat(vname, orig, len - n2);
		strcat(vname, ext);
	    }
	}
    } else if (ci == DUMMIFY) {
	char ext[16];

	if (aux < 0) {
	    sprintf(ext, "_m%d", abs(aux));
	} else {
	    sprintf(ext, "_%d", aux);
	}
	strcpy(vname, "D");
	strncat(vname, orig, len - strlen(ext) - 1);
	strcat(vname, ext);
    } else if (ci == INVERSE) {
	strcpy(vname, "i_");
	strncat(vname, orig, len - 2);
    } else if (ci == STDIZE) {
	if (aux < 0) {
	    strcpy(vname, "c_");
	} else {
	    strcpy(vname, "s_");
	}
	strncat(vname, orig, len - 2);
    } else if (ci == RESAMPLE) {
	strcpy(vname, "rs_");
	strncat(vname, orig, len - 3);
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
		      int ci, int aux)
{
    gchar *s = NULL;
    int err = 0;

    if (ci == DIFF) {
	s = g_strdup_printf(_("= first difference of %s"), parent);
    } else if (ci == LDIFF) {
	s = g_strdup_printf(_("= log difference of %s"), parent);
    } else if (ci == SDIFF) {
	s = g_strdup_printf(_("= seasonal difference of %s"), parent);
    } else if (ci == LOGS) {
	s = g_strdup_printf(_("= log of %s"), parent);
    } else if (ci == SQUARE) {
	s = g_strdup_printf(_("= %s squared"), parent);
    } else if (ci == LAGS) {
	/* @aux = lag */
	if (aux >= 0) {
	    s = g_strdup_printf("= %s(t - %d)", parent, aux);
	} else {
	    s = g_strdup_printf("= %s(t + %d)", parent, -aux);
	}
    } else if (ci == INVERSE) {
	s = g_strdup_printf("= 1/%s", parent);
    } else if (ci == STDIZE) {
	/* @aux = dfcorr */
	if (aux < 0) {
	    s = g_strdup_printf(_("= centered %s"), parent);
	} else {
	    s = g_strdup_printf(_("= standardized %s"), parent);
	}
    } else if (ci == RESAMPLE) {
	s = g_strdup_printf(_("= resampled %s"), parent);
    } else if (ci == HFDIFF) {
	s = g_strdup_printf(_("= high-frequency difference of %s"), parent);
    } else if (ci == HFLDIFF) {
	s = g_strdup_printf(_("= high-frequency log difference of %s"), parent);
    } else {
	err = 1;
    }

    if (!err) {
	strcpy(label, gretl_utf8_truncate(s, MAXLABEL-1));
	g_free(s);
    }

    return err;
}

/**
 * standard_lag_of:
 * @v: ID number of series to test.
 * @parent: ID of potential parent series.
 * @dset: dataset information.
 *
 * Returns: the lag order of series @v, if it is marked as
 * a lag of @parent, otherwise 0.
 */

int standard_lag_of (int v, int parent, const DATASET *dset)
{
    int pv = 0, ret = 0;

    if (dset == NULL || v <= 0 || v >= dset->v) {
	return 0;
    }

    if (series_get_transform(dset, v) == LAGS) {
	pv = series_get_parent_id(dset, v);
	if (pv == parent) {
	    ret = series_get_lag(dset, v);
	}
    }

    return ret;
}

/**
 * is_standard_diff:
 * @v: ID number of variable to test.
 * @dset: dataset information.
 * @parent: location to receive ID number of parent variable,
 * or NULL.
 *
 * Returns: 1 if the variable @v is marked as being the first
 * difference of some "parent" variable in the dataset,
 * otherwise 0.
 */

int is_standard_diff (int v, const DATASET *dset, int *parent)
{
    int pv = 0, ret = 0;

    if (v <= 0 || v >= dset->v) {
	return 0;
    }

    if (series_get_transform(dset, v) == DIFF) {
	pv = series_get_parent_id(dset, v);
	if (pv > 0) {
	    if (parent != NULL) {
		*parent = pv;
	    }
	    ret = 1;
	}
    }

    return ret;
}

static void make_xp_varname (char *vname, int v1, int v2,
			     const DATASET *dset,
			     int len)
{
    int n1 = strlen(dset->varname[v1]);
    int n2 = strlen(dset->varname[v2]);
    int cut = n1 + n2 + 1 - len;
    int cut1 = 0, cut2 = 0;

    if (cut > 0) {
	cut1 = cut2 = cut/2;
	if (cut % 2) {
	    cut2++;
	}
    }

    *vname = '\0';
    strncat(vname, dset->varname[v1], n1 - cut1);
    strcat(vname, "_");
    strncat(vname, dset->varname[v2], n2 - cut2);
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

/**
 * gretl_transforms_cleanup:
 *
 * Called by libgretl_cleanup(). Frees any memory allocated
 * as workspace for the creation of transformed variables.
 */

void gretl_transforms_cleanup (void)
{
    testvec(0);
}

/* write lagged values of variable v into xlag */

static int get_lag (int v, int lag, double *xlag,
		    const DATASET *dset)
{
    const double *x = dset->Z[v];
    int t1 = (lag > 0)? lag : 0;
    int t2 = dset->n - 1;
    int t, s;

    for (t=0; t<dset->n; t++) {
	xlag[t] = NADBL;
    }

    for (t=t1; t<=t2; t++) {
	s = t - lag;
	if (dataset_is_panel(dset) &&
	    t / dset->pd != s / dset->pd) {
	    continue;
	}
	if (s >= 0 && s < dset->n) {
	    xlag[t] = x[s];
	}
    }

    return 0;
}

/* write log of variable v into logvec */

static int get_log (int v, double *logvec, const DATASET *dset)
{
    double xx;
    int t, err = 0;

    for (t=dset->t1; t<=dset->t2 && !err; t++) {
	xx = dset->Z[v][t];
	if (na(xx) || xx <= 0.0) {
	    logvec[t] = NADBL;
	    set_gretl_warning(W_GENMISS);
	} else {
	    logvec[t] = log(xx);
	}
    }

    return err;
}

/* write some sort of difference of variable v into diffvec */

static int get_diff (int v, double *diffvec, int ci,
		     const DATASET *dset)
{
    double x0, x1;
    int t, t0, t1;

    t0 = (ci == SDIFF)? dset->pd : 1;
    t1 = (dset->t1 > t0)? dset->t1 : t0;

    for (t=t1; t<=dset->t2; t++) {
	if (dset->structure == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, dset)) {
	    continue;
	}

	x0 = dset->Z[v][t];
	x1 = dset->Z[v][t - t0];

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

/* write a matrix containing high-frequency differences of
   the series in @list over the current sample range
*/

static gretl_matrix *get_hf_diffs (const int *list,
				   int ci, double mult,
				   const DATASET *dset,
				   int *err)
{
    gretl_matrix *dX = NULL;
    double **Z = dset->Z;
    double x0, x1, dti;
    int n = list[0];
    int v1 = list[1];
    int i, vi, s, t, T;

    T = dset->t2 - dset->t1 + 1;
    dX = gretl_zero_matrix_new(T, n);

    if (dX == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (t=T-1; t>=0; t--) {
	s = t + dset->t1;
	for (i=0; i<n; i++) {
	    vi = list[i+1];
	    x0 = Z[vi][s];
	    if (i < n-1) {
		x1 = Z[list[i+2]][s];
	    } else if (s > 0) {
		x1 = Z[v1][s-1];
	    } else {
		x1 = NADBL;
	    }
	    dti = NADBL;
	    if (!na(x0) && !na(x1)) {
		if (ci == DIFF) {
		    dti = mult * (x0 - x1);
		} else if (x0 > 0 && x1 > 0) {
		    dti = mult * log(x0/x1);
		}
	    }
	    gretl_matrix_set(dX, t, i, dti);
	}
    }

    return dX;
}

/* orthogonal deviations */

static int get_orthdev (int v, double *xvec, const DATASET *dset)
{
    return orthdev_series(dset->Z[v], xvec, dset);
}

/* write square or cross-product into xvec */

static int get_xpx (int vi, int vj, double *xvec, const DATASET *dset)
{
    double xit, xjt;
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
	xit = dset->Z[vi][t];
	xjt = dset->Z[vj][t];
	if (na(xit) || na(xjt)) {
	    xvec[t] = NADBL;
	} else {
	    xvec[t] = xit * xjt;
	}
    }

    return 0;
}

/* write reciprocal into xvec */

static int get_inverse (int v, double *xvec, const DATASET *dset)
{
    double xt;
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
	xt = dset->Z[v][t];
	if (na(xt) || xt == 0.0) {
	    xvec[t] = NADBL;
	} else {
	    xvec[t] = 1.0 / xt;
	}
    }

    return 0;
}

/* write standardized series v into targ */

static int get_std (int v, int dfc, double *targ, const DATASET *dset)
{
    double *x = dset->Z[v];

    return standardize_series(x, targ, dfc, dset);
}

/* write resampled series into rsvec */

static int get_resampled (int v, double *rsvec,
			  const DATASET *dset,
			  const int *z)
{
    const double *x = dset->Z[v];
    int t, j, s;

    if (dataset_is_panel(dset)) {
	int n = panel_sample_size(dset);
	int i, T = dset->pd;

	s = dset->t1;
	for (i=0; i<n; i++) {
	    j = z[i] * T;
	    for (t=0; t<T; t++) {
		rsvec[s++] = x[j + t];
	    }
	}
    } else {
	s = 0;
	for (t=dset->t1; t<=dset->t2; t++) {
	    j = z[s++];
	    rsvec[t] = x[j];
	}
    }

    return 0;
}

/* write dummy for (v == value) into xvec */

static int get_discdum (int v, double val, double *xvec,
			const DATASET *dset)
{
    double xt;
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
	xt = dset->Z[v][t];
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
				       const double *x,
				       const char *label,
				       DATASET *dset,
				       int parent,
				       int origv)
{
    const char *vlabel = series_get_label(dset, v);
    int ret = VARNAME_DUPLICATE;
    int no_overwrite = 0;
    int t, ok = 0;

    if (vlabel != NULL && *label && *vlabel) {
	if (!strcmp(label, vlabel)) {
	    /* labels identical, so OK? */
	    ok = 1;
	} else {
	    /* labels differ */
	    no_overwrite = 1;
	}
    }

    if (!no_overwrite) {
	int vp = series_get_parent_id(dset, v);

	no_overwrite = parent != vp;
    }

#if TR_OVERWRITE
    if (!ok && v < origv && !no_overwrite) {
	ok = 1;
    }
#endif

    if (ok) {
#if TRDEBUG
	fprintf(stderr, "transform_handle_duplicate: updating var %d (%s)\n",
		v, dset->varname[v]);
#endif
	for (t=0; t<dset->n; t++) {
	    dset->Z[v][t] = x[t];
	}
	if (*label != '\0') {
	    series_set_label(dset, v, label);
	}
	series_set_transform(dset, v, ci);
	series_set_lag(dset, v, lag);
	series_zero_flags(dset, v);
	ret = VAR_EXISTS_OK;
    }

#if TRDEBUG
    if (!ok) {
	fprintf(stderr, "transform_handle_duplicate: collision with var %d (%s)\n",
		v, dset->varname[v]);
    }
#endif

    return ret;
}

static int
check_add_transform (int ci, int lag, int vnum, const double *x,
		     const char *vname, const char *label,
		     DATASET *dset, int parent, int origv)
{
    int t, ret = VAR_ADDED_OK;

    if (vnum < dset->v) {
	/* a variable of this name already exists */
	int chk = check_vals(x, dset->Z[vnum], dset->n);

	/* heuristic: we'll assume that if the variables have the same
	   name, and differ only in respect of one being more complete
	   than the other (having valid values where the other has
	   missing values), then we should not create a new variable,
	   but rather use the more complete series.

	   Watch out for cases where this is not the desired behavior!
	*/

	if (chk == VARS_IDENTICAL) {
	    ret = VAR_EXISTS_OK;
	} else if (chk == NEW_HAS_MISSING) {
	    /* is this right? */
	    ret = VAR_EXISTS_OK;
	} else if (chk == OLD_HAS_MISSING) {
	    for (t=0; t<dset->n; t++) {
		dset->Z[vnum][t] = x[t];
	    }
	    ret = VAR_EXISTS_OK;
	} else {
	    ret = transform_handle_duplicate(ci, lag, vnum, x, label, dset,
					     parent, origv);
	}
    } else {
	/* no var of this name, working from scratch */
	if (dataset_add_series(dset, 1)) {
	    ret = VAR_ADD_FAILED;
	} else {
	    strcpy(dset->varname[vnum], vname);
	    series_set_label(dset, vnum, label);
	    for (t=0; t<dset->n; t++) {
		dset->Z[vnum][t] = x[t];
	    }
	}
    }

    return ret;
}

#define RESPECT_LAG_LEVEL 1

static int get_lag_ID (int srcv, int lag, const DATASET *dset)
{
#if RESPECT_LAG_LEVEL
    int fd = gretl_function_depth();
#endif
    const char *parent, *vname = dset->varname[srcv];
    int i;

    for (i=1; i<dset->v; i++) {
#if RESPECT_LAG_LEVEL
	if (series_get_stack_level(dset, i) != fd) {
	    continue;
	}
#endif
	if (lag == series_get_lag(dset, i)) {
	    parent = series_get_parent_name(dset, i);
	    if (parent != NULL && !strcmp(vname, parent)) {
		return i;
	    }
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
			  DATASET *dset, int startlen, int origv,
			  int *idxvec)
{
    char vname[VNAMELEN] = {0};
    char label[MAXLABEL] = {0};
    int vno = -1, err = 0;
    int len, lag = 0;
    const char *srcname;
    double *vx;

    vx = testvec(dset->n);
    if (vx == NULL) {
	return -1;
    }

    if (ci == LAGS) {
	lag = aux;
	err = get_lag(v, lag, vx, dset);
    } else if (ci == LOGS) {
	err = get_log(v, vx, dset);
    } else if (ci == DIFF || ci == LDIFF || ci == SDIFF) {
	err = get_diff(v, vx, ci, dset);
    } else if (ci == ORTHDEV) {
	err = get_orthdev(v, vx, dset);
    } else if (ci == SQUARE) {
	/* "aux" = second variable number */
	err = get_xpx(v, aux, vx, dset);
    } else if (ci == DUMMIFY) {
	/* "x" = value for dummy */
	err = get_discdum(v, x, vx, dset);
    } else if (ci == INVERSE) {
	err = get_inverse(v, vx, dset);
    } else if (ci == STDIZE) {
	err = get_std(v, aux, vx, dset);
    } else if (ci == RESAMPLE) {
	err = get_resampled(v, vx, dset, idxvec);
    }

    if (err) {
	return -1;
    }

    if (ci == LAGS && (vno = get_lag_ID(v, aux, dset)) > 0) {
	/* special case: pre-existing lag */
	err = check_add_transform(ci, lag, vno, vx, vname, label, dset, v, origv);
	if (err != VAR_EXISTS_OK) {
	    vno = -1;
	}
	return vno;
    }

    if (ci == DUMMIFY || (ci == SQUARE && v != aux)) {
	/* special cases in respect of labeling */
	gchar *s = NULL;

	if (ci == SQUARE && v != aux) {
	    s = g_strdup_printf(_("= %s times %s"), dset->varname[v],
				dset->varname[aux]);
	} else if (is_string_valued(dset, v)) {
	    const char *sval = series_get_string_for_value(dset, v, x);

	    if (sval != NULL && *sval != '\0') {
		s = g_strdup_printf(_("dummy for %s = '%s'"), dset->varname[v], sval);
	    } else {
		s = g_strdup_printf(_("dummy for %s = %g"), dset->varname[v], x);
	    }
	} else {
	    s = g_strdup_printf(_("dummy for %s = %g"), dset->varname[v], x);
	}
	strcpy(label, gretl_utf8_truncate(s, MAXLABEL-1));
	g_free(s);
    } else {
	make_transform_label(label, dset->varname[v], ci, aux);
    }

    srcname = get_mangled_name_by_id(v);
    if (srcname == NULL) {
	srcname = dset->varname[v];
    }

#if TRDEBUG
    fprintf(stderr, "srcname = '%s', startlen = %d, v = %d\n",
	    srcname, startlen, v);
#endif

    for (len=startlen; len<=VNAMELEN; len++) {
	if (*vname != '\0' && len == VNAMELEN) {
	    /* last resort: hack the name */
	    make_varname_unique(vname, 0, dset);
	} else if (ci == SQUARE && v != aux) {
	    make_xp_varname(vname, v, aux, dset, len);
	} else if (ci == DUMMIFY && series_is_integer_valued(dset, v)) {
	    make_transform_varname(vname, srcname, ci, (int) x, len);
	} else {
	    make_transform_varname(vname, srcname, ci, aux, len);
	}
	vno = series_index(dset, vname);
	err = check_add_transform(ci, lag, vno, vx, vname, label, dset, v, origv);
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
	    series_set_parent(dset, vno, dset->varname[v]);
	    series_set_transform(dset, vno, ci);
	    series_set_discrete(dset, vno, 1);
	} else if (ci == LAGS || ci == DIFF || ci == LOGS) {
	    series_set_parent(dset, vno, dset->varname[v]);
	    series_set_transform(dset, vno, ci);
	}

	if (ci == LAGS) {
	    series_set_lag(dset, vno, aux);
	} else {
	    series_set_lag(dset, vno, 0);
	}
    }

    return vno;
}

/**
 * laggenr:
 * @v: ID number in dataset of source variable.
 * @lag: the order of the lag to create.
 * @dset: dataset struct.
 *
 * Creates the specified lag of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the lagged variable, or -1 on error.
 */

int laggenr (int v, int lag, DATASET *dset)
{
    int lno;

    if (lag > dset->n || -lag > dset->n) {
	gretl_errmsg_sprintf(_("Invalid lag order %d"), lag);
	lno = -1;
    } else if (v == 0) {
	gretl_errmsg_set(_("The constant cannot be lagged"));
	lno = -1;
    } else if (lag == 0) {
	lno = v;
    } else {
	lno = get_transform(LAGS, v, lag, 0.0, dset,
			    VNAMELEN - 3, dset->v, NULL);
    }

    return lno;
}

/**
 * laggenr_from_to:
 * @v: ID number in dataset of source variable.
 * @fromlag: start of lag/lead range.
 * @tolag: end of lag/lead range.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Creates the specified lags of variable @v if they do not
 * already exist.
 *
 * Returns: list of lag variables, or NULL or on error.
 */

int *laggenr_from_to (int v, int fromlag, int tolag,
		      DATASET *dset, int *err)
{
    int *llist;
    int i, lv, nlags = -1;
    int p, count_down = 0;

    if (fromlag > tolag) {
	nlags = fromlag - tolag + 1;
	count_down = 1;
    } else {
	nlags = tolag - fromlag + 1;
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

    p = fromlag;

    for (i=0; i<nlags; i++) {
	lv = laggenr(v, -p, dset);
	if (lv < 0) {
	    *err = E_DATA;
	    free(llist);
	    llist = NULL;
	    break;
	}
	llist[i+1] = lv;
	if (count_down) {
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
 * @dset: dataset struct.
 *
 * Creates the natural log of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the log variable, or -1 on error.
 */

int loggenr (int v, DATASET *dset)
{
    return get_transform(LOGS, v, 0, 0.0, dset,
			 VNAMELEN - 3, dset->v, NULL);
}

/**
 * invgenr:
 * @v: ID number in dataset of source variable.
 * @dset: dataset struct.
 *
 * Creates the reciprocal of variable @v if this variable does
 * not already exist.
 *
 * Returns: the ID number of the reciprocal, or -1 on error.
 */

int invgenr (int v, DATASET *dset)
{
    return get_transform(INVERSE, v, 0, 0.0, dset,
			 VNAMELEN - 3, dset->v, NULL);
}

/**
 * diffgenr:
 * @v: ID number in dataset of source variable.
 * @ci: DIFF (first difference), LDIFF (log difference) or SDIFF
 * (seasonal difference).
 * @dset: dataset struct.
 *
 * Creates the first difference (or log- or seasonal difference,
 * depending on the value of @ci) of variable @v, if the
 * differenced variable does not already exist.
 *
 * Returns: the ID number of the differenced variable, or -1 on error.
 */

int diffgenr (int v, int ci, DATASET *dset)
{
    if (ci != DIFF && ci != LDIFF && ci != SDIFF) {
	return -1;
    }

    if (ci == SDIFF && !dataset_is_seasonal(dset)) {
	return -1;
    }

    return get_transform(ci, v, 0, 0.0, dset,
			 VNAMELEN - 3, dset->v, NULL);
}

/**
 * xpxgenr:
 * @vi: ID number in dataset of first source variable.
 * @vj: ID number in dataset of second source variable.
 * @dset: dataset struct.
 *
 * Creates the cross product of variables @vi and @vj if this
 * variable does not already exist.
 *
 * Returns: the ID number of the cross-product variable,
 * or -1 on error.
 */

int xpxgenr (int vi, int vj, DATASET *dset)
{
    if (vi == vj) {
	if (gretl_isdummy(dset->t1, dset->t2, dset->Z[vi])) {
	    return -1;
	}
    }

    return get_transform(SQUARE, vi, vj, 0.0, dset,
			 VNAMELEN - 3, dset->v, NULL);
}

static int
get_starting_length (const int *list, DATASET *dset, int trim)
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
	len = strlen(dset->varname[list[i]]);
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
	    vni = dset->varname[list[i]];
	    for (j=i+1; j<=list[0] && !conflict; j++) {
		vnj = dset->varname[list[j]];
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
    char tmp[32], sfx[32];
    int n, err;

    if (get_mangled_name_by_id(v)) {
	return 0;
    }

    *tmp = *sfx = '\0';
    /* we know the bit out here (+nc) must be distinct? */
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
transform_preprocess_list (int *list, const DATASET *dset, int f)
{
    int maxc = VNAMELEN - 3;
    int longnames = 0;
    char **S = NULL;
    int i, v, ok;
    int err = 0;

    if (f == SQUARE || f == LDIFF || f == SDIFF || f == RESAMPLE) {
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
	    if (gretl_isdummy(dset->t1, dset->t2, dset->Z[v])) {
		ok = 0;
	    }
	} else if (f == LAGS) {
	    if (v == 0) {
		gretl_errmsg_set(_("The constant cannot be lagged"));
		ok = 0;
	    }
	} else if (f == DUMMIFY) {
	    ok = v > 0 && accept_as_discrete(dset, v, 0);
	}

	if (!ok) {
	    gretl_list_delete_at_pos(list, i--);
	} else if (strlen(dset->varname[v]) > maxc) {
	    strings_array_add(&S, &longnames, dset->varname[v]);
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
			v = series_index(dset, S[i]);
			herr = make_mangled_name(v, S[i], maxc);
			v = series_index(dset, S[j]);
			herr += make_mangled_name(v, S[j], maxc);
		    }
		}
	    }
	}

	strings_array_free(S, longnames);
    }

    return err;
}

/**
 * list_loggenr:
 * @list: on entry, list of variables to process; on exit,
 * holds the ID numbers of the generated variables.
 * @dset: dataset struct.
 *
 * Generates and adds to the data set the natural logs of the
 * variables given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_loggenr (int *list, DATASET *dset)
{
    int origv = dset->v;
    int tnum, i, j, v;
    int startlen;
    int l0 = 0;
    int err;

    err = transform_preprocess_list(list, dset, LOGS);
    if (err) {
	return err;
    }

    startlen = get_starting_length(list, dset, 2);

    j = 1;
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	tnum = get_transform(LOGS, v, 0, 0.0, dset, startlen,
			     origv, NULL);
	if (tnum > 0) {
	    list[j++] = tnum;
	    l0++;
	}
    }

    list[0] = l0;

    destroy_mangled_names();

    return (l0 > 0)? 0 : E_LOGS;
}

/**
 * list_stdgenr:
 * @list: on entry, list of variables to process; on exit,
 * holds the ID numbers of the generated variables.
 * @dset: dataset struct.
 *
 * Generates and adds to the dataset standardized versions
 * of the series given in @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_stdgenr (int *list, DATASET *dset, gretlopt opt)
{
    int origv = dset->v;
    int tnum, i, j, v;
    int startlen;
    int dfc = 1;
    int l0 = 0;
    int err;

    err = transform_preprocess_list(list, dset, STDIZE);
    if (err) {
	return err;
    }

    if (opt & OPT_C) {
	dfc = -1;
    } else if (opt & OPT_N) {
	dfc = 0;
    }

    startlen = get_starting_length(list, dset, 2);

    j = 1;
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	tnum = get_transform(STDIZE, v, dfc, 0.0, dset, startlen,
			     origv, NULL);
	if (tnum > 0) {
	    list[j++] = tnum;
	    l0++;
	}
    }

    list[0] = l0;

    destroy_mangled_names();

    return (l0 > 0)? 0 : E_DATA;
}

static int *make_lags_list (int *list, int order)
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

static int lag_wanted (int p, const gretl_matrix *v, int n)
{
    int ret = 1;

    if (v != NULL) {
	int i;

	ret = 0; /* reverse the assumption */
	for (i=0; i<n; i++) {
	    if (v->val[i] == p) {
		ret = 1;
		break;
	    }
	}
    }

    return ret;
}

static int process_hf_lags_input (int hf_min, int hf_max,
				  DATASET *dset,
				  int compfac,
				  int *lf_min, int *lf_max,
				  int *skip_first,
				  int *n_terms)
{
    if (dset->pd != 1 && dset->pd != 12 && dset->pd != 4) {
	return E_PDWRONG;
    }

    *lf_min = (int) ceil(hf_min / (double) compfac);
    *lf_max = (int) ceil(hf_max / (double) compfac);
    *skip_first = (hf_min - 1) % compfac;
    *n_terms = hf_max - hf_min + 1;

    if (*skip_first < 0) {
	/* handle leads */
	*skip_first = compfac + *skip_first;
    }

#if 0
    fprintf(stderr, "hf_min = %d, lf_min = %d\n", hf_min, *lf_min);
    fprintf(stderr, "hf_max = %d, lf_max = %d\n", hf_max, *lf_max);
    fprintf(stderr, "skip_first = %d\n", *skip_first);
    fprintf(stderr, "n_terms = %d\n", *n_terms);
#endif

    return 0;
}

/**
 * list_laggenr:
 * @plist: on entry, pointer to list of variables to process.  On exit
 * the list holds the ID numbers of the lag variables.
 * @lmin: minimum lag to include (defaults to 1).
 * @lmax: maximum lag to include (or 0 for automatic).
 * @lvec: (alternative to @lmin, @lmax) vector holding lags to generate.
 * @dset: dataset struct.
 * @compfac: compaction factor (for MIDAS lists only, otherwise give 0).
 * @opt: may contain OPT_L to order the list by lag rather than
 * by variable.
 *
 * Generates and adds to the data set lagged values of the
 * variables given in the list pointed to by @plist.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int list_laggenr (int **plist, int lmin, int lmax,
		  const gretl_matrix *lvec,
		  DATASET *dset, int compfac,
		  gretlopt opt)
{
    int origv = dset->v;
    int *list = *plist;
    int *laglist = NULL;
    int l, i, j, v, lv;
    int startlen, l0 = 0;
    int skip_first = 0;
    int veclen = 0;
    int n_terms = 0;
    int err = 0;

    if (compfac < 0) {
	return E_INVARG;
    } else if (compfac > 0) {
	/* MIDAS: we want @lmin and @lmax, not @lvec */
	if (lvec != NULL) {
	    return E_INVARG;
	}
	if (!gretl_is_midas_list(list, dset)) {
	    gretl_warnmsg_set(_("The argument does not seem to be a MIDAS list"));
	}
    }

    if (lvec != NULL) {
	/* got a vector of lags in @lvec */
	veclen = gretl_vector_get_length(lvec);
	if (veclen == 0) {
	    err = E_INVARG;
	} else {
	    lmin = lvec->val[0];
	    lmax = lvec->val[veclen-1];
	    n_terms = veclen;
	}
    } else if (compfac > 0) {
	/* the MIDAS case */
	int hmin = lmin, hmax = lmax;

	if (hmax < hmin) {
	    err = E_INVARG;
	} else {
	    err = process_hf_lags_input(hmin, hmax, dset, compfac, &lmin,
					&lmax, &skip_first, &n_terms);
	}
    } else {
	/* non-MIDAS, using order = @lmax, ignoring @lmin */
	if (lmax < 0 || lmax > dset->n) {
	    gretl_errmsg_sprintf(_("Invalid lag order %d"), lmax);
	    return E_INVARG;
	} else if (lmax == 0) {
	    lmax = default_lag_order(dset);
	}
	lmin = 1;
	n_terms = lmax;
    }

    if (err) {
	return err;
    }

    /* Note: at this point @lmin and @lmax are defined in terms of
       the frequency of @dset, even when we're handling MIDAS data.
       In the regular case, @n_terms holds the number of lag terms
       wanted per series in the input @list, but in the MIDAS case
       @n_terms is the total number of hf lag terms wanted (since
       the input @list actually represents a single series).
    */

    err = transform_preprocess_list(list, dset, LAGS);
    if (err) {
	return err;
    }

    if (compfac) {
	laglist = gretl_list_new(n_terms);
    } else {
	laglist = make_lags_list(list, n_terms);
    }

    if (laglist == NULL) {
	destroy_mangled_names();
	return E_ALLOC;
    }

    startlen = get_starting_length(list, dset, (lmax > 9)? 3 : 2);

    j = 1;

    if (compfac > 0) {
	/* MIDAS high-frequency lags, by lags */
	int hfpmax = n_terms + skip_first;
	int mp, hfp = 0;

	for (l=lmin; l<=lmax; l++) {
	    for (i=1; i<=list[0]; i++) {
		hfp++;
		if (hfp <= skip_first) {
		    continue;
		} else if (hfp > hfpmax) {
		    break;
		}
		v = list[i];
		lv = get_transform(LAGS, v, l, 0.0, dset,
				   startlen, origv, NULL);
		if (lv > 0) {
		    mp = series_get_midas_period(dset, v);
		    series_set_midas_period(dset, lv, mp);
		    mp = series_get_midas_freq(dset, v);
		    series_set_midas_freq(dset, lv, mp);
		    laglist[j++] = lv;
		    l0++;
		}
	    }
	}
    } else if (opt & OPT_L) {
	/* order the list by lags */
	for (l=lmin; l<=lmax; l++) {
	    if (!lag_wanted(l, lvec, veclen)) {
		continue;
	    }
	    for (i=1; i<=list[0]; i++) {
		v = list[i];
		lv = get_transform(LAGS, v, l, 0.0, dset,
				   startlen, origv, NULL);
		if (lv > 0) {
		    laglist[j++] = lv;
		    l0++;
		}
	    }
	}
    } else {
	/* order by variable */
	for (i=1; i<=list[0]; i++) {
	    v = list[i];
	    for (l=lmin; l<=lmax; l++) {
		if (!lag_wanted(l, lvec, veclen)) {
		    continue;
		}
		lv = get_transform(LAGS, v, l, 0.0, dset,
				   startlen, origv, NULL);
		if (lv > 0) {
		    laglist[j++] = lv;
		    l0++;
		}
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
 * @dset: data information struct.
 *
 * Returns: the default lag order for generating lags, performing
 * autocorrelation test, and so on.
 */

int default_lag_order (const DATASET *dset)
{
    int order = 1;

    if (!dataset_is_panel(dset)) {
	order = (dset->pd < 52)? dset->pd : 14;
    }

    return order;
}

/**
 * list_diffgenr:
 * @list: on entry, list of variables to process; on exit,
 * ID numbers of the generated variables.
 * @ci: must be DIFF, LDIFF or SDIFF.
 * @dset: pointer to dataset struct.
 *
 * Generate differences of the variables in @list, and add them
 * to the data set.  If @ci is DIFF these are ordinary first
 * differences; if @ci is LDIFF they are log differences; and
 * if @ci is SDIFF they are seasonal differences.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int list_diffgenr (int *list, int ci, DATASET *dset)
{
    int origv = dset->v;
    int i, v, startlen;
    int tnum, l0 = 0;
    int err;

    if (list[0] == 0) {
	return 0;
    }

    if (ci != DIFF && ci != LDIFF && ci != SDIFF) {
	return 1;
    }

    if (ci == SDIFF && !dataset_is_seasonal(dset)) {
	return E_PDWRONG;
    }

    err = transform_preprocess_list(list, dset, ci);
    if (err) {
	return err;
    }

    startlen = get_starting_length(list, dset, (ci == DIFF)? 2 : 3);

    for (i=1; i<=list[0] && !err; i++) {
	v = list[i];
	tnum = get_transform(ci, v, 0, 0.0, dset,
			     startlen, origv, NULL);
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

static int check_hf_difflist (const int *list,
			      const DATASET *dset,
			      int ci, int *n_add)
{
    char vname[VNAMELEN];
    int i, v, li, err = 0;

    for (i=1; i<=list[0] && !err; i++) {
	li = list[i];
	make_transform_varname(vname, dset->varname[li],
			       ci, 0, VNAMELEN - 8);
	v = current_series_index(dset, vname);
	if (v > 0) {
	    *n_add -= 1;
	} else if (gretl_is_user_var(vname)) {
	    gretl_errmsg_sprintf(_("%s: collides with existing object name"),
				 vname);
	    err = E_TYPES;
	}
    }

    return err;
}

/**
 * hf_list_diffgenr:
 * @list: on entry, midas list of variables to process; on exit,
 * ID numbers of the generated variables.
 * @ci: must be DIFF or LDIFF.
 * @parm: optional scalar multiplier.
 * @dset: pointer to dataset struct.
 *
 * Generate high-frequency differences of the variables in @list
 * and add them to the data set.  If @ci is DIFF these are ordinary
 * first differences; if @ci is LDIFF they are log differences.
 * If @parm is not NA then the new series are all multiplied by
 * the specified value (as in multiplication of log-differences
 * by 100).
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int hf_list_diffgenr (int *list, int ci, double parm, DATASET *dset)
{
    gretl_matrix *dX = NULL;
    int n_add, hfci = 0;
    int set_midas = 1;
    int err = 0;

    if (list[0] == 0) {
	return 0;
    }

    if (ci == DIFF) {
	hfci = HFDIFF;
    } else if (ci == LDIFF) {
	hfci = HFLDIFF;
    } else {
	return E_INVARG;
    }

    if (!dataset_is_time_series(dset)) {
	return E_PDWRONG;
    }

    n_add = list[0];
    err = check_hf_difflist(list, dset, ci, &n_add);
    if (err) {
	return err;
    }

    if (!gretl_is_midas_list(list, dset)) {
	gretl_warnmsg_set(_("The argument does not seem to be a MIDAS list"));
	set_midas = 0;
    }

    if (na(parm)) {
	/* set to default */
	parm = 1.0;
    }

    dX = get_hf_diffs(list, ci, parm, dset, &err);

    if (!err) {
	char vname[VNAMELEN];
	char label[MAXLABEL];
	const char *parent;
	int i, s, t, n = dX->cols;
	int vi, v0 = dset->v;
	int offset = 0;

	err = dataset_add_series(dset, n_add);

	if (!err) {
	    for (i=0; i<n; i++) {
		parent = dset->varname[list[i+1]];
		make_transform_varname(vname, parent, ci, 0,
				       VNAMELEN - 8);
		vi = current_series_index(dset, vname);
		if (vi < 0) {
		    vi = v0 + offset;
		    offset++;
		    strcpy(dset->varname[vi], vname);
		    make_transform_label(label, parent, hfci, 0);
		    series_record_label(dset, vi, label);
		}
		s = 0;
		for (t=dset->t1; t<=dset->t2; t++) {
		    dset->Z[vi][t] = gretl_matrix_get(dX, s++, i);
		}
		list[i+1] = vi;
	    }
	    list[0] = n;
	}
	gretl_matrix_free(dX);
    }

    if (!err && set_midas) {
	gretl_list_set_midas(list, dset);
    }

    return err;
}

/**
 * list_orthdev:
 * @list: list of variables to process.
 * @dset: dataset struct.
 *
 * Generate orthogonal deviations of the variables in @list, and add
 * them to the data set.
 *
 * Returns: 0 on success, error code on error.
 */

int list_orthdev (int *list, DATASET *dset)
{
    int origv = dset->v;
    int i, v, startlen;
    int tnum, l0 = 0;
    int err;

    if (list[0] == 0) {
	return 0;
    }

    if (!dataset_is_panel(dset)) {
	return E_PDWRONG;
    }

    err = transform_preprocess_list(list, dset, ORTHDEV);
    if (err) {
	return err;
    }

    startlen = get_starting_length(list, dset, 2);

    for (i=1; i<=list[0] && !err; i++) {
	v = list[i];
	tnum = get_transform(ORTHDEV, v, 0, 0.0, dset,
			     startlen, origv, NULL);
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
 * @dset: dataset struct.
 * @opt: If OPT_O, both squares and cross-products are generated,
 * otherwise only squares.
 *
 * Generates and adds to the data set squares and (if @opt is OPT_O)
 * cross-products of the variables given in the list pointed to
 * by @plist.
 *
 * Returns: 0 on success, error code on error.
 */

int list_xpxgenr (int **plist, DATASET *dset, gretlopt opt)
{
    int origv = dset->v;
    int *list = *plist;
    int *xpxlist = NULL;
    int tnum, i, j, k, vi, vj;
    int startlen, l0;
    int err;

    err = transform_preprocess_list(list, dset, SQUARE);
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

    startlen = get_starting_length(list, dset, 3);
    xpxlist[0] = 0;

    k = 1;
    for (i=1; i<=l0; i++) {
	vi = list[i];
	tnum = get_transform(SQUARE, vi, vi, 0.0, dset,
			     startlen, origv, NULL);
	if (tnum > 0) {
	    xpxlist[k++] = tnum;
	    xpxlist[0] += 1;
	}
	if (opt & OPT_O) {
	    for (j=i+1; j<=l0; j++) {
		vj = list[j];
		tnum = xpxgenr(vi, vj, dset);
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

int list_resample (int *list, DATASET *dset)
{
    int origv = dset->v;
    int *z = NULL;
    int tnum, i, j, vi;
    int n, startlen, l0;
    int k1, k2;
    int err;

    err = transform_preprocess_list(list, dset, RESAMPLE);
    if (err) {
	return err;
    }

    if (dataset_is_panel(dset)) {
	k1 = dset->t1 / dset->pd;
	k2 = dset->t2 / dset->pd;
    } else {
	k1 = dset->t1;
	k2 = dset->t2;
    }

    n = k2 - k1 + 1;

    z = malloc(n * sizeof *z);
    if (z == NULL) {
	return E_ALLOC;
    }

    /* generate uniform drawings from sample range */
    gretl_rand_int_minmax(z, n, k1, k2);

    l0 = list[0];
    startlen = get_starting_length(list, dset, 3);
    list[0] = 0;

    j = 1;
    for (i=1; i<=l0; i++) {
	vi = list[i];
	tnum = get_transform(RESAMPLE, vi, 0, 0.0, dset,
			     startlen, origv, z);
	if (tnum > 0) {
	    list[j++] = tnum;
	    list[0] += 1;
	}
    }

    free(z);
    destroy_mangled_names();

    return (list[0] > 0)? 0 : E_DATA;
}

/**
 * list_dropcoll:
 * @list: on entry, list of variables to process; on exit,
 * the original list minus collinear terms.
 * @dset: dataset struct.
 *
 * Drop collinear variables from @list.
 *
 * Returns: 0 on success, error code on error.
 */

int list_dropcoll (int *list, double eps, DATASET *dset)
{
    gretl_matrix *R = NULL;
    gretl_matrix *X = NULL;
    const double *x;
    int i, T, k;
    int err = 0;

    if (list == NULL) {
	return E_DATA;
    } else if (list[0] < 2) {
	/* no-op */
	return 0;
    }

    if (eps < 0) {
	return E_DATA;
    } else if (na(eps)) {
	eps = R_DIAG_MIN;
    }

    for (i=list[0]; i>0; i--) {
	x = dset->Z[list[i]];
	if (gretl_iszero(dset->t1, dset->t2, x)) {
	    err = gretl_list_delete_at_pos(list, i);
	}
    }

    X = gretl_matrix_data_subset(list, dset, dset->t1, dset->t2,
				 M_MISSING_SKIP, &err);
    if (err) {
	return err;
    }

    T = X->rows;
    k = X->cols;

    if (T < k) {
	gretl_errmsg_sprintf(_("A minimum of %d observations is required"), k);
	err = E_TOOFEW;
    } else {
	R = gretl_matrix_alloc(k, k);
	if (R == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_QR_decomp(X, R);
	}
    }

    if (!err) {
	double rii;
	int j = 1;

	for (i=0; i<k && !err; i++, j++) {
	    rii = gretl_matrix_get(R, i, i);
	    if (fabs(rii) < eps) {
		err = gretl_list_delete_at_pos(list, j--);
	    }
	}
    }

    gretl_matrix_free(X);
    gretl_matrix_free(R);

    return err;
}

#define DUMDEBUG 0

static int real_list_dumgenr (int **plist, DATASET *dset,
			      double oddval, gretlopt opt)
{
    int origv = dset->v;
    int *list = *plist;
    int *tmplist = NULL;
    double *x = NULL;
    int i, j, t, n;
    int startlen;
    int err;

    err = transform_preprocess_list(list, dset, DUMMIFY);
    if (err) {
	return err;
    }

    tmplist = gretl_null_list();
    if (tmplist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    startlen = get_starting_length(list, dset, 3);

    for (i=1; i<=list[0] && !err; i++) {
	int vi = list[i];
	int nuniq, tnum;
	int jmin, jmax;
	double xt;

	n = 0;
	for (t=dset->t1; t<=dset->t2; t++) {
	    xt = dset->Z[vi][t];
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
        fprintf(stderr, "oddval = %g\n", oddval);
#endif

	for (j=jmin; j<jmax && !err; j++) {
	    if (na(oddval) || x[j] != oddval) {
                tnum = get_transform(DUMMIFY, vi, j+1, x[j], dset,
                                     startlen, origv, NULL);
#if DUMDEBUG
                fprintf(stderr, "%s: for value %g tnum = %d\n",
                        dset->varname[vi], x[j], tnum);
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

    if (!err && !(opt & OPT_A) && tmplist[0] == 0) {
	/* don't fail here if doing auto-dummify (OPT_A) */
	gretl_errmsg_set(_("dummify: no suitable variables were found"));
#if DUMDEBUG
	printlist(list, "list for which dummify failed");
#endif
	err = E_DATA;
    }

    free(x);

 bailout:

    if (!err && tmplist[0] > 0) {
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
 * @dset: dataset struct.
 * @opt: can include OPT_F to drop the first value, OPT_L to drop
 * the last value.
 *
 * For each of the variables given in the list to which @plist
 * points, generates and adds to the data set k dummy variables
 * coding for the k distinct values of the variable in question.
 * All these variables must have already been marked as discrete.
 * If the OPT_F or OPT_L option is given, either the first or
 * the last value of each variable is taken as the "base" and is
 * not given a dummy encoding (that is, only k - 1 dummies are
 * added for each variable).
 *
 * Returns: 0 on success, error code on error.
*/

int list_dumgenr (int **plist, DATASET *dset, gretlopt opt)
{
    return real_list_dumgenr(plist, dset, NADBL, opt);
}

/**
 * dumgenr_with_oddval:
 * @plist: pointer to list of variables to process; on exit
 * the list holds the ID numbers of the generated dummies.
 * @dset: dataset struct.
 * @oddval: value which should be skipped when encoding the
 * input values as dummies.
 *
 * For each of the variables given in the list to which @plist
 * points, generates and adds to the data set k dummy variables
 * coding for the k distinct values of the variable in question.
 * All these variables must have already been marked as discrete.
 * if @oddval is not %NADBL, it is treated as the omitted
 * category and only k - 1 dummies are added for each variable).
 *
 * Returns: 0 on success, error code on error.
*/

int dumgenr_with_oddval (int **plist, DATASET *dset, double oddval)
{
    return real_list_dumgenr(plist, dset, oddval, OPT_NONE);
}

/**
 * auto_dummify_list:
 * @plist: pointer to list of variables to process.
 * @dset: dataset struct.
 *
 * Produces a list in which any coded series in the original
 * list are replaced by a set of dummy variables. If the
 * original list contains no such series it is not altered.
 * In generating dummy series, the omitted category is the
 * minimum value of the coding variable.
 *
 * Returns: 0 on success, error code on error.
*/

int auto_dummify_list (int **plist, DATASET *dset)
{
    int *orig = *plist;
    int *list = NULL;
    int *dlist = NULL;
    int i, vi;
    int err = 0;

    /* empty starting point for full list */
    list = gretl_null_list();

    for (i=1; i<=orig[0] && !err; i++) {
	vi = orig[i];
	if (series_is_coded(dset, vi)) {
	    /* coded series: split into dummies */
	    dlist = gretl_list_new(1);
	    dlist[1] = vi;
	    err = real_list_dumgenr(&dlist, dset, NADBL, OPT_F | OPT_A);
	    if (!err) {
		err = gretl_list_append_list(&list, dlist);
	    }
	    free(dlist);
	} else if (!in_gretl_list(list, vi)) {
	    /* non-coded, not already present */
	    list = gretl_list_append_term(&list, vi);
	    if (list == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (err) {
	free(list);
    } else if (list != *plist) {
	/* replace original list */
	free(*plist);
	*plist = list;
    }

    return err;
}

/**
 * list_makediscrete:
 * @list: list of variables to process.
 * @dset: data information struct.
 * @opt: if OPT_R, reverse the operation.
 *
 * Sets the variables given in @list as discrete, unless
 * opt is OPT_R, in which case the variables are set as
 * continuous.
 *
 * Returns: 0 on success, error code on error.
 */

int list_makediscrete (const int *list, DATASET *dset, gretlopt opt)
{
    int disc = !(opt & OPT_R);
    int i, v, err = 0;

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0) {
	    series_set_discrete(dset, v, disc);
	}
    }

    return err;
}
