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
#include "gretl_func.h"
#include "gretl_scalar.h"
#include "gretl_string_table.h"
#include "libset.h"
#include "dbread.h"

#define DDEBUG 0
#define FULLDEBUG 0

static int dataset_changed;

int check_dataset_is_changed (void)
{
    int ret = dataset_changed;

    dataset_changed = 0;
    return ret;
}

void set_dataset_is_changed (void)
{
    dataset_changed = 1;
}

/**
 * free_Z:
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Does a deep free on the data matrix.
 */

void free_Z (double **Z, DATAINFO *pdinfo)
{
    if (Z != NULL && pdinfo != NULL) {
	int i;

#if DDEBUG
	fprintf(stderr, "Freeing Z (%p): %d vars\n", (void *) Z, 
		pdinfo->v);
#endif

	for (i=0; i<pdinfo->v; i++) {
	    free(Z[i]);
	}
	free(Z);
    }
}

/**
 * dataset_destroy_obs_markers:
 * @pdinfo: data information struct.
 *
 * Frees any allocated observation markers for @pdinfo.
 */

void dataset_destroy_obs_markers (DATAINFO *pdinfo)
{
    int i;

    if (pdinfo->S != NULL) {
	for (i=0; i<pdinfo->n; i++) { 
	   free(pdinfo->S[i]); 
	}
	free(pdinfo->S);
	pdinfo->S = NULL;
	pdinfo->markers = NO_MARKERS;
    } 
}

static void free_sorted_markers (DATAINFO *pdinfo, int v)
{
    VARINFO *vinfo = pdinfo->varinfo[v];
    int i;

    if (vinfo->sorted_markers != NULL) {
	for (i=0; i<pdinfo->n; i++) {
	    free(vinfo->sorted_markers[i]);
	}
	free(vinfo->sorted_markers);
	vinfo->sorted_markers = NULL;
    }    
}

void dataset_destroy_panel_info (DATAINFO *pdinfo)
{
    if (pdinfo->paninfo != NULL) {
	free(pdinfo->paninfo->unit);
	free(pdinfo->paninfo->period);
	free(pdinfo->paninfo->padmask);
	free(pdinfo->paninfo);
	pdinfo->paninfo = NULL;
    }
}

static void free_varinfo (DATAINFO *pdinfo, int v)
{
    free_sorted_markers(pdinfo, v);
    free(pdinfo->varinfo[v]);
}

void set_sorted_markers (DATAINFO *pdinfo, int v, char **S)
{
    free_sorted_markers(pdinfo, v);
    pdinfo->varinfo[v]->sorted_markers = S;
}

/**
 * clear_datainfo:
 * @pdinfo: data information struct.
 * @code: either %CLEAR_FULL or %CLEAR_SUBSAMPLE.
 *
 * Frees the allocated content of a data information struct.
 */

void clear_datainfo (DATAINFO *pdinfo, int code)
{
    int i;

    if (pdinfo == NULL) return;

    if (pdinfo->S != NULL) {
	dataset_destroy_obs_markers(pdinfo);
    } 

    if (pdinfo->submask != NULL) {
	free_subsample_mask(pdinfo->submask);
	pdinfo->submask = NULL;
    }

    if (pdinfo->restriction != NULL) {
	free(pdinfo->restriction);
	pdinfo->restriction = NULL;
    }	

    if (pdinfo->paninfo != NULL) {
	dataset_destroy_panel_info(pdinfo);
    }    

    /* if this is not a sub-sample datainfo, free varnames, labels, etc. */
    if (code == CLEAR_FULL) {
	if (pdinfo->varname != NULL) {
	    for (i=0; i<pdinfo->v; i++) {
		free(pdinfo->varname[i]); 
	    }
	    free(pdinfo->varname);
	    pdinfo->varname = NULL;
	}
	if (pdinfo->varinfo != NULL) {
	    for (i=0; i<pdinfo->v; i++) {
		free_varinfo(pdinfo, i);
	    }
	    free(pdinfo->varinfo);
	    pdinfo->varinfo = NULL;
	}
	if (pdinfo->descrip != NULL) {
	    free(pdinfo->descrip);
	    pdinfo->descrip = NULL;
	}

	maybe_free_full_dataset(pdinfo);

	pdinfo->v = pdinfo->n = 0;
    } 
}

/**
 * destroy_dataset:
 * @Z: data array.
 * @pdinfo: dataset information struct.
 *
 * Frees all resources associated with @Z and @pdinfo.
 */

void destroy_dataset (double **Z, DATAINFO *pdinfo)
{
    free_Z(Z, pdinfo);

    if (pdinfo != NULL) {
	clear_datainfo(pdinfo, CLEAR_FULL); 
	free(pdinfo);
    }
}

/**
 * copy_dataset_obs_info:
 * @targ: target dataset information struct.
 * @src: source dataset information struct.
 *
 * Sets the "date" or observations information in @targ to that
 * found in @src.
 */

void copy_dataset_obs_info (DATAINFO *targ, const DATAINFO *src)
{
    strcpy(targ->stobs, src->stobs);
    strcpy(targ->endobs, src->endobs);
    targ->sd0 = src->sd0;
    targ->pd = src->pd;
    targ->structure = src->structure;
    targ->decpoint = src->decpoint;
}

/**
 * dataset_obs_info_default:
 * @pdinfo: dataset information struct.
 *
 * Sets the "date" or observations information in @pdinfo to a
 * simple default of cross-sectional data, observations 1 to n,
 * where n is the %n element (number of observations) in @pdinfo.
 */

void dataset_obs_info_default (DATAINFO *pdinfo)
{
    strcpy(pdinfo->stobs, "1");
    sprintf(pdinfo->endobs, "%d", pdinfo->n);
    pdinfo->sd0 = 1.0;
    pdinfo->pd = 1;
    pdinfo->structure = CROSS_SECTION;
    pdinfo->decpoint = '.';
}

/**
 * dataset_allocate_obs_markers:
 * @pdinfo: dataset information struct
 *
 * Allocates space in @pdinfo for strings indentifying the
 * observations and initializes all of the markers to empty
 * strings.  Note that These strings have a fixed maximum 
 * length of #OBSLEN - 1.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_allocate_obs_markers (DATAINFO *pdinfo)
{
    char **S = NULL;
    int err = 0;

    if (pdinfo->S == NULL) {
	/* not already allocated */
	S = strings_array_new_with_length(pdinfo->n, OBSLEN);
	if (S == NULL) {
	    err = E_ALLOC;
	} else {
	    pdinfo->S = S;
	}
    }

    if (pdinfo->S != NULL) {
	pdinfo->markers = REGULAR_MARKERS;
    }

    return err;
}

/**
 * dataset_allocate_panel_info:
 * @pdinfo: dataset information struct
 *
 * Allocates space in @pdinfo for two indices representing
 * the unit or group and time-period, respectively, of each 
 * observation in a panel data set.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_allocate_panel_info (DATAINFO *pdinfo)
{
    PANINFO *pan;
    int i;

#if DDEBUG
    fprintf(stderr, "dataset_allocate_panel_info: pdinfo at %p\n",
	    (void *) pdinfo);
#endif

    /* just in case, clean out any previous stuff */
    dataset_destroy_panel_info(pdinfo);  

    pan = malloc(sizeof *pan);
    if (pan == NULL) {
	return E_ALLOC;
    }

    pan->unit = pan->period = NULL;
    pan->padmask = NULL;

    pan->unit = malloc(pdinfo->n * sizeof *pan->unit);
    pan->period = malloc(pdinfo->n * sizeof *pan->period);

    if (pan->unit == NULL || pan->period == NULL) {
	free(pan->unit);
	free(pan->period);
	free(pan);
	return E_ALLOC;
    }

    for (i=0; i<pdinfo->n; i++) {
	pan->unit[i] = pan->period[i] = -1; /* deliberately invalid */
    }

    pan->nunits = 0;
    pan->Tmin = pan->Tmax = 0;
    pan->olen = 0;

    pdinfo->paninfo = pan;

    return 0;
}

static int resize_panel_info (DATAINFO *pdinfo, int n)
{
    PANINFO *pan = pdinfo->paninfo;
    int *unit;
    int *period;

    unit = realloc(pan->unit, n * sizeof *unit);
    if (unit == NULL) {
	return E_ALLOC;
    }
    pan->unit = unit;

    period = realloc(pan->period, n * sizeof *period);
    if (period == NULL) {
	return E_ALLOC;
    }
    pan->period = period;

    if (n > pdinfo->n) {
	int uadd = (n - pdinfo->n) / pdinfo->pd;
	int i, t, j = pdinfo->paninfo->nunits;
	int s = pdinfo->n;

	for (i=0; i<uadd; i++) {
	    for (t=0; t<pdinfo->pd; t++) {
		pan->unit[s + i] = j;
		pan->period[s + i] = t;
		s++;
	    }
	    j++;
	}
    }

    pdinfo->paninfo->nunits = n / pdinfo->pd;

    return 0;
}

/**
 * dataset_add_default_panel_indices:
 * @pdinfo: dataset information struct.
 *
 * Adds a pair of indices for panel unit and panel period.
 * The default is that both are zero-based and increase
 * consecutively, per unit and per period respectively.
 * This function assumes a balanced panel.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_add_default_panel_indices (DATAINFO *pdinfo)
{
    int nunits, T;
    int i, s, t, err;

    if (pdinfo->n % pdinfo->pd != 0) {
	fprintf(stderr, "dataset_add_default_panel_indices: "
		"error: pdinfo->n %% pdinfo->pd = %d\n", 
		pdinfo->n % pdinfo->pd);	
	return E_DATA;
    }

    err = dataset_allocate_panel_info(pdinfo);

    if (!err) {
	char test[32];

	T = pdinfo->pd;
	nunits = pdinfo->n / pdinfo->pd;
	s = 0;
	for (i=0; i<nunits; i++) {
	    for (t=0; t<T; t++) {
		pdinfo->paninfo->unit[s] = i;
		pdinfo->paninfo->period[s] = t;
		s++;
	    }
	}
	pdinfo->paninfo->nunits = nunits;
	pdinfo->paninfo->Tmin = T;
	pdinfo->paninfo->Tmax = T;

	sprintf(test, "%d", T);
	pdinfo->paninfo->olen = strlen(test);
    }

    return err;
}

/**
 * dataset_finalize_panel_indices:
 * @pdinfo: dataset information struct.
 *
 * Having already added a pair of indices for panel unit 
 * and panel period, check these for consistency and
 * calculate the number of panel units and the minimum
 * and maximum observations per unit.  If it turns out
 * there's only one unit, or only one period, in the
 * dataset, then it's not really a panel: we destroy
 * the panel info and return %E_PDWRONG.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_finalize_panel_indices (DATAINFO *pdinfo)
{
    PANINFO *pan = pdinfo->paninfo;
    char test[32];
    int u0, ubak;
    int i, Ti, len;

    if (pan == NULL) {
	return E_DATA;
    }

    pan->nunits = 1;
    pan->Tmax = 0;
    pan->Tmin = 999999;
    pan->olen = 0;

    u0 = ubak = pan->unit[0];
    Ti = 0;

    /* basic validity check */
    for (i=0; i<pdinfo->n; i++) {
	if (pan->unit[i] < 0 || pan->period[i] < 0) {
	    gretl_errmsg_set("Panel index information is corrupted");
	    return E_DATA;
	}
	sprintf(test, "%d", pan->period[i] + 1);
	len = strlen(test);
	if (len > pan->olen) {
	    pan->olen = len;
	}
    }

    for (i=0; i<pdinfo->n; i++) {
	if (pan->unit[i] == ubak) {
	    Ti++;
	} else {
	    pan->nunits += 1;
	    if (Ti > pan->Tmax) {
		pan->Tmax = Ti;
	    }
	    if (Ti < pan->Tmin) {
		pan->Tmin = Ti;
	    }	    
	    Ti = 1;
	    ubak = pan->unit[i];
	}
    }

    if (pan->nunits == 1 || pan->Tmax < 2) {
	dataset_destroy_panel_info(pdinfo);
	pdinfo->structure = CROSS_SECTION;
	return E_PDWRONG;
    }

#if DDEBUG
    fprintf(stderr, "dataset_finalize_panel_indices:\n "
	    "nunits=%d, Tmin=%d, Tmax=%d, olen=%d\n", pan->nunits,
	    pan->Tmin, pan->Tmax, pan->olen);
#endif

    return 0;
}

static void gretl_varinfo_init (VARINFO *vinfo)
{
    vinfo->label[0] = '\0';
    vinfo->display_name[0] = '\0';
    vinfo->parent[0] = '\0';
    vinfo->flags = 0;
    vinfo->transform = 0;
    vinfo->lag = 0;
    vinfo->compact_method = COMPACT_NONE;
    vinfo->line_width = 1;
    vinfo->sorted_markers = NULL;
    vinfo->stack_level = gretl_function_depth();
}

/**
 * copy_varinfo:
 * @targ: target to which to copy.
 * @src: source to copy from.
 *
 * Copies all relevant information from @src to @targ.
 */

void copy_varinfo (VARINFO *targ, const VARINFO *src)
{
    if (src == NULL || targ == NULL) {
	return;
    }

    strcpy(targ->label, src->label);
    strcpy(targ->display_name, src->display_name);
    strcpy(targ->parent, src->parent);
    targ->flags = src->flags;
    targ->transform = src->transform;
    targ->lag = src->lag;
    targ->compact_method = src->compact_method;
    targ->stack_level = src->stack_level;
    targ->line_width = src->line_width;
    targ->sorted_markers = NULL;
}

/**
 * dataset_allocate_varnames:
 * @pdinfo: dataset information struct.
 *
 * Given a blank @pdinfo, which should have been obtained using
 * datainfo_new(), allocate space for the names of variables.
 * The @v member of @pdinfo (number of variables) must be
 * set before calling this function.
 * 
 * Returns: 0 on sucess, %E_ALLOC on failure.
 */

int dataset_allocate_varnames (DATAINFO *pdinfo)
{
    int i, j, v = pdinfo->v;
    int err = 0;

    pdinfo->varname = strings_array_new_with_length(v, VNAMELEN);
    if (pdinfo->varname == NULL) {
	return E_ALLOC;
    }

    pdinfo->varinfo = malloc(v * sizeof *pdinfo->varinfo);
    if (pdinfo->varinfo == NULL) {
	free(pdinfo->varname);
	return E_ALLOC;
    }

    for (i=0; i<v; i++) {
	pdinfo->varinfo[i] = malloc(sizeof **pdinfo->varinfo);
	if (pdinfo->varinfo[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(pdinfo->varinfo[j]);
	    }
	    free(pdinfo->varinfo);
	    pdinfo->varinfo = NULL;
	    err = E_ALLOC;
	    break;
	} else {
	    gretl_varinfo_init(pdinfo->varinfo[i]);
	}
    }

    if (!err) {
	strcpy(pdinfo->varname[0], "const");
	strcpy(VARLABEL(pdinfo, 0), _("auto-generated constant"));
    } 

    return err;
}

/**
 * datainfo_new:
 *
 * Creates a new data information struct pointer from scratch,
 * properly initialized as empty.
 * 
 * Returns: pointer to data information struct, or NULL on error.
 */

DATAINFO *datainfo_new (void)
{
    DATAINFO *dinfo;

    dinfo = malloc(sizeof *dinfo);
    if (dinfo == NULL) {
	return NULL;
    }

    dinfo->v = 0;
    dinfo->n = 0;
    dinfo->pd = 1;
    dinfo->sd0 = 1.0;
    dinfo->t1 = 0;
    dinfo->t2 = 0;
    dinfo->stobs[0] = '\0';
    dinfo->endobs[0] = '\0';

    dinfo->varname = NULL;
    dinfo->varinfo = NULL;    
    dinfo->paninfo = NULL;

    dinfo->markers = NO_MARKERS;  
    dinfo->delim = ',';
    dinfo->decpoint = '.';

    dinfo->S = NULL;
    dinfo->descrip = NULL;
    dinfo->submask = NULL;
    dinfo->restriction = NULL;

    dinfo->structure = CROSS_SECTION;

    return dinfo;
}

/**
 * create_new_dataset:
 * @pZ: pointer to data matrix.
 * @nvar: number of variables.
 * @nobs: number of observations per variable 
 * @markers: 1 if there are case markers for the observations, 0
 * otherwise.
 *
 * Creates a new data information struct corresponding to a given
 * data matrix.
 * 
 * Returns: pointer to data information struct, or %NULL on error.
 */

DATAINFO *
create_new_dataset (double ***pZ, int nvar, int nobs, int markers)
{
    DATAINFO *pdinfo;

    pdinfo = malloc(sizeof *pdinfo);
    if (pdinfo == NULL) return NULL;

    pdinfo->v = nvar;
    pdinfo->n = nobs;

    *pZ = NULL;

    if (start_new_Z(pZ, pdinfo, 0)) {
	free(pdinfo);
	return NULL;
    }

    if (markers) {
	if (dataset_allocate_obs_markers(pdinfo)) {
	    free_datainfo(pdinfo);
	    return NULL;
	}
    } 

    dataset_obs_info_default(pdinfo);
    pdinfo->descrip = NULL;

    return pdinfo;
}

DATAINFO *create_auxiliary_dataset (double ***pZ, int nvar, int nobs)
{
    DATAINFO *pdinfo;

    pdinfo = create_new_dataset(pZ, nvar, nobs, 0);

    if (pdinfo != NULL) {
	pdinfo->descrip = gretl_strdup("auxiliary");
    }

    return pdinfo;
}

/**
 * allocate_Z:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information struct.
 *
 * Allocates the two-dimensional array to which @pZ points,
 * based on the %v (number of variables) and %n (number of
 * observations) members of @pdinfo.  The variable at 
 * position 0 is initialized to all 1s; other variables
 * are initialized to #NADBL.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int allocate_Z (double ***pZ, const DATAINFO *pdinfo)
{
    double **Z;
    int i, t;
    int err = 0;

    if (*pZ != NULL) {
	free(*pZ);
    }

    Z = doubles_array_new(pdinfo->v, pdinfo->n);

    if (Z == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<pdinfo->v; i++) {
	    for (t=0; t<pdinfo->n; t++) {
		Z[i][t] = (i == 0)? 1.0 : NADBL;
	    }
	}
    }

    *pZ = Z;

    return err;
}

/**
 * start_new_Z:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @resample: 1 if we're sub-sampling from a full data set, 0 otherwise.
 *
 * Initializes the data matrix pointed to by @pZ (adding the constant in
 * position 0) and the data information struct @pdinfo.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 */

int start_new_Z (double ***pZ, DATAINFO *pdinfo, int resample)
{
    if (allocate_Z(pZ, pdinfo)) {
	return E_ALLOC;
    }

    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    if (resample) {
	pdinfo->varname = NULL;
	pdinfo->varinfo = NULL;
    } else if (dataset_allocate_varnames(pdinfo)) {
	free_Z(*pZ, pdinfo);
	*pZ = NULL;
	return E_ALLOC;
    }

    pdinfo->S = NULL;
    pdinfo->markers = NO_MARKERS;
    pdinfo->delim = ',';
    pdinfo->descrip = NULL;
    pdinfo->paninfo = NULL;
    pdinfo->submask = NULL;
    pdinfo->restriction = NULL;
    
    return 0;
}

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

/* Allow for the possibility of centered seasonal dummies: usually
   xon = 1 and xoff = 0, but in the centered case xon = 1 - 1/pd 
   and xoff = -1/pd.
*/

static int get_xon_xoff (const double *x, int n, int pd, double *xon, double *xoff)
{
    double cfac = 1.0 / pd;
    double xc = 1.0 - cfac, yc = -cfac;
    double x0 = 999, y0 = 999;
    int t, ret = 1;

    for (t=0; t<n && ret; t++) {
	if (x[t] == 1.0) {
	    if (x0 == 999) x0 = 1.0;
	    else if (x[t] != x0) ret = 0;
	} else if (x[t] == 0.0) {
	    if (y0 == 999) y0 = 0.0;
	    else if (x[t] != y0) ret = 0;
	} else if (x[t] == xc) {
	    if (x0 == 999) x0 = xc;
	    else if (x[t] != x0) ret = 0;
	} else if (x[t] == yc) {
	    if (y0 == 999) y0 = yc;
	    else if (x[t] != y0) ret = 0;
	} else {
	    ret = 0;
	}
    }

    if (ret) {
	*xon = x0;
	*xoff = y0;
    }

    return ret;
}

static int real_periodic_dummy (const double *x, int n,
				int *pd, int *offset, 
				double *pxon, double *pxoff)
{
    double xon = 1.0, xoff = 0.0;
    int onbak = 0;
    int gap = 0;
    int trail = 0;
    int t, m = n - 1, ret = 1;

    if (!get_xon_xoff(x, n, *pd, &xon, &xoff)) {
	return 0;
    }

    *pd = -1;
    *offset = -1;
    trail = 0;

    /* find number of trailing "off" values */
    for (t=n-1; t>0; t--) {
	if (x[t] == xoff) {
	    trail++;
	} else {
	    if (x[t] == xon) {
		m = t;
	    } else {
		ret = 0;
	    }
	    break;
	}
    }

    /* check for dummyhood and periodicity */
    for (t=0; t<=m && ret; t++) {
	if (x[t] == xoff) {
	    onbak = 0;
	    gap++;
	} else if (x[t] == xon) {
	    if (onbak) {
		ret = 0;
	    } else if (*offset < 0) {
		*offset = gap;
	    } else if (*pd < 0) {
		*pd = gap + 1;
		if (*pd < *offset + 1) {
		    ret = 0;
		}
	    } else if (gap != *pd - 1) {
		ret = 0;
	    } else if (gap < trail) {
		ret = 0;
	    }
	    gap = 0;
	    onbak = 1;
	} else {
	    ret = 0;
	    break;
	}
    }

    if (ret && pxon != NULL && pxoff != NULL) {
	*pxon = xon;
	*pxoff = xoff;
    }

    return ret;
}

/**
 * is_periodic_dummy:
 * @x: array to examine.
 * @pdinfo: pointer to dataset information struct.
 *
 * Returns: 1 if @x is a periodic dummy variable,
 * 0 otherwise.
 */

int is_periodic_dummy (const double *x, const DATAINFO *pdinfo)
{
    int offset, pd = pdinfo->pd;

    return real_periodic_dummy(x, pdinfo->n, &pd, &offset, NULL, NULL);
}

static int is_linear_trend (const double *x, int n)
{
    int t, ret = 1;
    
    for (t=1; t<n; t++) {
	if (x[t] != x[t-1] + 1.0) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

static int is_quadratic_trend (const double *x, int n)
{
    double t2;
    int t, ret = 1;
    
    for (t=0; t<n; t++) {
	t2 = (t + 1) * (t + 1);
	if (x[t] != t2) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * is_trend_variable:
 * @x: array to examine.
 * @n: number of elements in array.
 *
 * Returns: 1 if @x is a simple linear trend variable, with each
 * observation equal to the preceding observation plus 1, or
 * if @x is a quadratic trend starting at 1 for the first 
 * observation in the data set, and 0 otherwise.
 */

int is_trend_variable (const double *x, int n)
{
    int ret = 0;

    if (is_linear_trend(x, n)) {
	ret = 1;
    } else if (is_quadratic_trend(x, n)) {
	ret = 1;
    }

    return ret;
}

static void 
maybe_extend_trends (double **Z, const DATAINFO *pdinfo, int oldn)
{
    int i, t;

    for (i=1; i<pdinfo->v; i++) {
	if (is_linear_trend(Z[i], oldn)) {
	    for (t=oldn; t<pdinfo->n; t++) {
		Z[i][t] = Z[i][t-1] + 1.0;
	    }
	} else if (is_quadratic_trend(Z[i], oldn)) {
	    for (t=oldn; t<pdinfo->n; t++) {
		Z[i][t] = (t + 1) * (t + 1);
	    }
	}	    
    }
}

static void 
maybe_extend_dummies (double **Z, const DATAINFO *pdinfo, int oldn)
{
    int pd = pdinfo->pd;
    double xon = 1.0, xoff = 0.0;
    int offset;
    int i, t;

    for (i=1; i<pdinfo->v; i++) {
	if (real_periodic_dummy(Z[i], oldn, &pd, &offset, &xon, &xoff)) {
	    for (t=oldn; t<pdinfo->n; t++) {
		Z[i][t] = ((t - offset) % pd)? xoff : xon;
	    }
	}
    }
}

/**
 * dataset_add_observations:
 * @newobs: number of observations to add.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: use %OPT_A to attempt to recognize and
 * automatically extend simple deterministic variables such 
 * as a time trend and periodic dummy variables; 
 * use %OPT_D to drop any observation markers rather than
 * expanding the set of markers and padding it out with
 * dummy values.
 *
 * Extends all series in the dataset by the specified number of
 * extra observations.  The added values are initialized to
 * the missing value code, #NADBL, with the exception of
 * simple deterministic variables when %OPT_A is given.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_add_observations (int newobs, double ***pZ, DATAINFO *pdinfo,
			      gretlopt opt)
{
    double *x;
    int oldn = pdinfo->n;
    int i, t, bign;

    if (newobs <= 0) {
	return 0;
    }

    if (dataset_is_panel(pdinfo) && newobs % pdinfo->pd != 0) {
	return E_PDWRONG;
    }

    bign = pdinfo->n + newobs;

    for (i=0; i<pdinfo->v; i++) {
	x = realloc((*pZ)[i], bign * sizeof *x);
	if (x == NULL) {
	    return E_ALLOC;
	}
	(*pZ)[i] = x;
	for (t=pdinfo->n; t<bign; t++) {
	    (*pZ)[i][t] = (i == 0)? 1.0 : NADBL;
	}	    
    }
    
    if (pdinfo->markers && pdinfo->S != NULL) {
	if (opt & OPT_D) {
	    dataset_destroy_obs_markers(pdinfo);
	} else {
	    if (reallocate_markers(pdinfo, bign)) {
		return E_ALLOC;
	    }
	    for (t=oldn; t<bign; t++) {
		sprintf(pdinfo->S[t], "%d", t + 1);
	    }
	}
    }

    if (pdinfo->paninfo != NULL) {
	if (resize_panel_info(pdinfo, bign)) {
	    return E_ALLOC;
	}
    }
    
    if (pdinfo->t2 == pdinfo->n - 1) {
	pdinfo->t2 = bign - 1;
    }

    pdinfo->n = bign;

    if (opt & OPT_A) {
	maybe_extend_trends(*pZ, pdinfo, oldn);
	maybe_extend_dummies(*pZ, pdinfo, oldn);
    }

    /* does daily data need special handling? */
    ntodate(pdinfo->endobs, bign - 1, pdinfo);

    return 0;
}

/**
 * dataset_drop_observations:
 * @n: number of observations to drop.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Deletes @n observations from the end of each series in the 
 * dataset.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_drop_observations (int n, double ***pZ, DATAINFO *pdinfo)
{
    double *x;
    int i, newn;

    if (n <= 0) {
	return 0;
    }

    if (dataset_is_panel(pdinfo) && n % pdinfo->pd != 0) {
	return E_PDWRONG;
    }

    newn = pdinfo->n - n;

    for (i=0; i<pdinfo->v; i++) {
	x = realloc((*pZ)[i], newn * sizeof *x);
	if (x == NULL) {
	    return E_ALLOC;
	}
	(*pZ)[i] = x;
    }
    
    if (pdinfo->markers && pdinfo->S != NULL) {
	if (reallocate_markers(pdinfo, newn)) {
	    return E_ALLOC;
	}
    }

    if (pdinfo->paninfo != NULL) {
	if (resize_panel_info(pdinfo, newn)) {
	    return E_ALLOC;
	}
    }    

    if (pdinfo->t2 > newn - 1) {
	pdinfo->t2 = newn - 1;
    }

    pdinfo->n = newn;

    /* does daily data need special handling? */
    ntodate(pdinfo->endobs, newn - 1, pdinfo);

    return 0;
}

/**
 * dataset_shrink_obs_range:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Truncates the range of observations in the dataset, based on
 * the current values of the %t1 and %t2 members of @pdinfo.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_shrink_obs_range (double ***pZ, DATAINFO *pdinfo)
{
    int head = pdinfo->t1;
    int tail = pdinfo->n - 1 - pdinfo->t2;
    int err = 0;

    if (head > 0) {
	int mvsize, rem = pdinfo->n - head;
	double **Z = *pZ;
	int i;

	mvsize = rem * sizeof **Z;
	for (i=0; i<pdinfo->v; i++) {
	    memmove(Z[i], Z[i] + head, mvsize);
	}

	if (pdinfo->markers && pdinfo->S != NULL) {
	    for (i=0; i<head; i++) {
		free(pdinfo->S[i]);
	    }
	    mvsize = rem * sizeof *pdinfo->S;
	    memmove(pdinfo->S, pdinfo->S + head, mvsize);
	}

	if (pdinfo->paninfo != NULL) {
	    mvsize = rem * sizeof *pdinfo->paninfo->unit;
	    memmove(pdinfo->paninfo->unit, pdinfo->paninfo->unit + head, mvsize);
	    memmove(pdinfo->paninfo->period, pdinfo->paninfo->period + head, mvsize);
	}

	/* FIXME panel data */

	if (pdinfo->structure == CROSS_SECTION) {
	    ntodate(pdinfo->stobs, 0, pdinfo);
	} else {
	    ntodate(pdinfo->stobs, pdinfo->t1, pdinfo);
	    pdinfo->sd0 = get_date_x(pdinfo->pd, pdinfo->stobs);
	}
	pdinfo->t1 = 0;
	pdinfo->t2 -= head;
	pdinfo->n -= head;
	ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);
    }

    if (tail > 0) {
	err = dataset_drop_observations(tail, pZ, pdinfo);
    }

    return err;
}

static int 
dataset_expand_varinfo (int newvars, DATAINFO *pdinfo)
{
    char **varname = NULL;
    VARINFO **varinfo = NULL;
    int v = pdinfo->v;
    int bigv = v + newvars;
    int i;

    varname = realloc(pdinfo->varname, bigv * sizeof *varname);
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
	varinfo = realloc(pdinfo->varinfo, bigv * sizeof *varinfo);
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

    pdinfo->v += newvars;

    return 0;
}

static int real_add_series (int newvars, double *x,
			    double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    int v = pdinfo->v;
    int n = pdinfo->n;
    int i, t;
    int err = 0;

    if (newvars == 0) {
	/* no-op */
	return 0;
    }

    newZ = realloc(*pZ, (v + newvars) * sizeof *newZ); 

#if DDEBUG
    fprintf(stderr, "real_add_series: add %d vars, Z = %p\n",
	    newvars, (void *) newZ);
#endif
	    

    if (newZ == NULL) {
	err = E_ALLOC;
    } else {
	*pZ = newZ;
    }
    
    if (!err) {
	if (newvars == 1 && x != NULL) {
	    /* new var is pre-allocated */
	    newZ[v] = x;
	} else {
	    for (i=0; i<newvars && !err; i++) {
		newZ[v+i] = malloc(n * sizeof **newZ);
		if (newZ[v+i] == NULL) {
		    err = E_ALLOC;
		} else {
		    for (t=0; t<n; t++) {
			newZ[v+i][t] = 0.0;
		    }
		}
	    }
	}
    }

    if (!err) {
	err = dataset_expand_varinfo(newvars, pdinfo);
    }

    return err;
}

/**
 * dataset_add_series:
 * @newvars: number of series to add.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Adds space for the specified number of additional series
 * to the dataset.  It is the caller's responsibility to
 * initialize the numerical values of the new series.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int 
dataset_add_series (int newvars, double ***pZ, DATAINFO *pdinfo)
{
    return real_add_series(newvars, NULL, pZ, pdinfo);
}

/**
 * dataset_add_allocated_series:
 * @x: one-dimensional data array.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Adds @x as an additional series in the dataset.
 * The array @x is not copied; it should be treated as
 * belonging to @pZ after this operation.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int 
dataset_add_allocated_series (double *x, double ***pZ, DATAINFO *pdinfo)
{
    return real_add_series(1, x, pZ, pdinfo);
}

/**
 * dataset_add_series_as:
 * @x: array to be added.
 * @newname: name to give the new variable.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Adds to the dataset a new series with name @newname and
 * values given by @x.  The new variable is added at one
 * level "deeper" (in terms of function execution) than the
 * current level.  This is for use with user-defined functions.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_series_as (double *x, const char *newname,
			   double ***pZ, DATAINFO *pdinfo)
{
    int v, t, err = 0;

    if (pdinfo->varinfo == NULL) {
	gretl_errmsg_set(_("Please open a data file first"));
	return 1;
    }

#if DDEBUG
    fprintf(stderr, "dataset_add_series_as: incoming Z=%p, name='%s'\n",
	    (void *) *pZ, newname);
#endif

    err = real_add_series(1, NULL, pZ, pdinfo);

    if (!err) {
	v = pdinfo->v - 1;
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[v][t] = x[t];
	}
	strcpy(pdinfo->varname[v], newname);
	STACK_LEVEL(pdinfo, v) += 1;
    }

    return err;
}

/**
 * dataset_copy_variable_as:
 * @v: index number of variable to copy.
 * @newname: name to give the copy.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Makes a copy of variable @v under the name @newname.
 * The copy exists in a variable namespace one level "deeper"
 * (in terms of function execution) than the variable being copied. 
 * This is for use with user-defined functions: a variable
 * supplied to a function as an argument is copied into the
 * function's namespace under the name it was given as a
 * parameter.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_copy_variable_as (int v, const char *newname,
			      double ***pZ, DATAINFO *pdinfo)
{
    int t, err;

    err = real_add_series(1, NULL, pZ, pdinfo);

    if (!err) {
	int vnew = pdinfo->v - 1;

	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[vnew][t] = (*pZ)[v][t];
	}
	strcpy(pdinfo->varname[vnew], newname);
	STACK_LEVEL(pdinfo, vnew) += 1;
	/* FIXME other varinfo stuff? */
#if 0
	fprintf(stderr, "copied var %d ('%s', level %d) as var %d ('%s', level %d): ",
		v, pdinfo->varname[v], STACK_LEVEL(pdinfo, v),
		vnew, newname, STACK_LEVEL(pdinfo, vnew));
	fprintf(stderr, "Z[%d][0] = %g\n", vnew, (*pZ)[vnew][0]);
#endif
    }

    return err;
}

enum {
    DROP_NORMAL,
    DROP_SPECIAL
};

/* DROP_SPECIAL is used when deleting variables from the "full" shadow
   of a sub-sampled dataset, after deleting those same variables from
   the sub-sampled version: in that case we don't mess with the
   pointer elements of the DATAINFO struct, because these are shared
   between the full and sub-sampled versions.
*/

static int 
shrink_dataset_to_size (double ***pZ, DATAINFO *pdinfo, int nv, int drop)
{
    double **newZ;

#if DDEBUG
    fprintf(stderr, "shrink_dataset_to_size: nv=%d, got pZ at %p, pdinfo at %p\n"
	    " drop = %s\n", nv, (void *) pZ, (void *) pdinfo, 
	    (drop == DROP_NORMAL)? "DROP_NORMAL" : "DROP_SPECIAL");
#endif

    if (drop == DROP_NORMAL) {
	char **varname;
	VARINFO **varinfo;
    
	varname = realloc(pdinfo->varname, nv * sizeof *varname);
	if (varname == NULL) {
	    return E_ALLOC;
	}
	pdinfo->varname = varname;

	varinfo = realloc(pdinfo->varinfo, nv * sizeof *varinfo);
	if (varinfo == NULL) {
	    return E_ALLOC;
	}
	pdinfo->varinfo = varinfo;
    }

    newZ = realloc(*pZ, nv * sizeof *newZ); 
    if (newZ == NULL) {
	return E_ALLOC;
    }
    *pZ = newZ;

    pdinfo->v = nv;

    return 0;
}

static int vars_renumbered (const int *list, DATAINFO *pdinfo,
			    int dmin)
{
    int i, ndel = 0;

    for (i=dmin; i<pdinfo->v; i++) {
	if (in_gretl_list(list, i)) {
	    ndel++;
	} else if (ndel > 0 && !var_is_hidden(pdinfo, i)) {
	    return 1;
	}
    }

    return 0;
}

int overwrite_err (const char *name)
{
    gretl_errmsg_sprintf("The variable %s is read-only", name);

    return E_DATA;
}

/**
 * series_is_parent:
 * @pdinfo: dataset information.
 * @v: ID number of variable to test.
 * 
 * Returns: 1 if variable @v is "parent" to a transformed
 * variable (e.g. a log, lag or difference), othewise 0.
 */

int series_is_parent (const DATAINFO *pdinfo, int v)
{
    const char *s = pdinfo->varname[v];
    int i;

    for (i=1; i<pdinfo->v; i++) {
	if (i != v && !strcmp(s, pdinfo->varinfo[i]->parent)) {
	    return 1;
	}
    }

    return 0;
}

/**
 * dataset_rename_series:
 * @pdinfo: dataset information.
 * @v: ID number of the variable to be renamed.
 * @name: new name to give the variable.
 * 
 * Returns: 0 on sucess, non-zero on error.
 */

int dataset_rename_series (DATAINFO *pdinfo, int v, const char *name)
{
    if (v < 0 || v >= pdinfo->v) {
	return E_DATA;
    }

    if (object_is_const(pdinfo->varname[v]) ||
	series_is_parent(pdinfo, v)) {
	return overwrite_err(pdinfo->varname[v]);
    }

    if (strcmp(pdinfo->varname[v], name)) {
	strcpy(pdinfo->varname[v], name);
	set_dataset_is_changed();
    }

    return 0;
}

static int real_drop_listed_vars (int *list, double ***pZ, 
				  DATAINFO *pdinfo, int *renumber,
				  int drop, PRN *prn)
{
    double **Z = *pZ;
    int oldv = pdinfo->v, vmax = pdinfo->v;
    char vname[VNAMELEN] = {0};
    int d0, d1;
    int delmin = oldv;
    int i, v, ndel = 0; 
    int err = 0;

    if (renumber != NULL) {
	*renumber = 0;
    }

    if (list == NULL || list[0] == 0) {
	/* no-op */
	return 0;
    }

    d0 = list[0];

    check_variable_deletion_list(list, pdinfo);
    d1 = list[0];
    if (prn != NULL && d1 == 1) {
	strcpy(vname, pdinfo->varname[list[1]]);
    }

    if (d1 == 0) {
	goto finish;
    }

#if DDEBUG
    fprintf(stderr, "real_drop_listed_variables: dropping %d vars:\n",
	    list[0]);
#endif

    /* check that no vars to be deleted are marked "const", and do
       some preliminary accounting while we're at it */

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0 && v < oldv) {
	    if (object_is_const(pdinfo->varname[v])) {
		return overwrite_err(pdinfo->varname[v]);
	    }
	    if (v < delmin) {
		delmin = v;
	    }
	    ndel++;
	}
    }

    if (ndel == 0) {
	return 0;
    }

    if (renumber != NULL) {
	*renumber = vars_renumbered(list, pdinfo, delmin);
    }    

#if DDEBUG
    fprintf(stderr, "real_drop_listed_variables: lowest ID of deleted var"
	    " = %d\n", delmin);
#endif

    /* free and set to NULL all the vars to be deleted */
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v > 0 && v < oldv) {
	    free(Z[v]);
	    Z[v] = NULL;
	    if (drop == DROP_NORMAL) {
		free(pdinfo->varname[v]);
		free(pdinfo->varinfo[v]);
	    }
	}
    }

    /* repack pointers if necessary */

    for (v=1; v<vmax; v++) {
	if (Z[v] == NULL) {
	    int gap = 1;

	    for (i=v+1; i<vmax; i++) {
		if (Z[i] == NULL) {
		    gap++;
		} else {
		    break;
		}
	    }

	    if (i < vmax) {
		vmax -= gap;
		for (i=v; i<vmax; i++) {
		    if (drop == DROP_NORMAL) {
			pdinfo->varname[i] = pdinfo->varname[i + gap];
			pdinfo->varinfo[i] = pdinfo->varinfo[i + gap];
		    }
		    Z[i] = Z[i + gap];
		}		    
	    } else {
		/* deleting all subsequent vars: done */
		break;
	    }
	}
    }

    err = shrink_dataset_to_size(pZ, pdinfo, oldv - ndel, drop);

 finish:

    /* report results, if appropriate */

    if (!err && prn != NULL) {
	if (d0 == d1) {
	    if (gretl_messages_on()) {
		if (*vname != '\0') {
		    pprintf(prn, "Deleted %s", vname);
		} else {
		    pprintf(prn, "Deleted %d variables", d1);
		}
		pputc(prn, '\n');
	    }
	} else {
	    if (d1 == 0) {
		pputs(prn, "No variables deleted");
	    } else if (*vname != '\0') {
		pprintf(prn, "Deleted %s", vname);
	    } else {
		pprintf(prn, "Deleted %d variables", d1);
	    }
	    pputs(prn, " (some data were in use)\n");
	    
	}
    }

    return err;
}

static int *make_dollar_list (DATAINFO *pdinfo, int *err)
{
    int *list = NULL;
    int i;

    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->varname[i][0] == '$') {
	    list = gretl_list_append_term(&list, i);
	    if (list == NULL) {
		*err = E_ALLOC;
		break;
	    }
	}
    }

    return list;
}

/**
 * dataset_drop_listed_variables:
 * @list: list of variable to drop, by ID number.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @renumber: location for return of information on whether
 * remaining variables have been renumbered as a result, or
 * %NULL.
 * @prn: pointer to printing struct.
 *
 * Deletes the variables given in @list from the dataset.  Remaining
 * variables may have their ID numbers changed as a consequence. If
 * @renumber is not %NULL, this location receives 1 in case variables
 * have been renumbered, 0 otherwise.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_drop_listed_variables (int *list, double ***pZ, 
				   DATAINFO *pdinfo, 
				   int *renumber,
				   PRN *prn)
{
    int *dlist = NULL;
    int free_dlist = 0;
    int lastvar[2];
    int err = 0;

    if (pdinfo->n == 0 || pdinfo->v == 0) {
	return E_NODATA;
    }

    if (list == NULL) {
	/* signal to drop internal "$" variables */
	dlist = make_dollar_list(pdinfo, &err);
	if (err) {
	    return err;
	} else if (dlist == NULL) {
	    /* no-op */
	    return 0;
	}
	free_dlist = 1;
    } else if (list[0] == 0) {
	/* signal to drop last variable */
	lastvar[0] = 1;
	lastvar[1] = pdinfo->v - 1;
	dlist = lastvar;
    } else {
	dlist = list;
    }

    err = real_drop_listed_vars(dlist, pZ, pdinfo, renumber,
				DROP_NORMAL, prn);

    if (dlist[0] > 0) {
	if (!err) {
	    err = gretl_lists_revise(dlist, 0);
	}

	if (!err && complex_subsampled()) {
	    double ***fZ = fetch_full_Z();
	    DATAINFO *fdinfo = fetch_full_datainfo();

	    err = real_drop_listed_vars(dlist, fZ, fdinfo, NULL,
					DROP_SPECIAL, NULL);
	    reset_full_Z(fZ);
	}
    }

    if (free_dlist) {
	free(dlist);
    }

    return err;
}

/**
 * dataset_drop_variable:
 * @v: ID number of variable to drop.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Deletes variable @v from the dataset.  
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_drop_variable (int v, double ***pZ, DATAINFO *pdinfo) 
{
    int list[2] = {1, v};

    if (v <= 0 || v >= pdinfo->v) {
	return E_DATA;
    }

    return dataset_drop_listed_variables(list, pZ, pdinfo, NULL, NULL);
}

/**
 * dataset_destroy_hidden_variables:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @vmin: do not drop variables with ID numbers less than this.
 *
 * Deletes from the dataset any "hidden" variables that have
 * been added automatically (for example, auto-generated variables
 * used for the x-axis in graph plotting), and that have ID
 * numbers greater than or equal to @vmin.  Never deletes the
 * automatically generated constant (ID number 0).
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_destroy_hidden_variables (double ***pZ, DATAINFO *pdinfo,
				      int vmin)
{
    int i, nhid = 0;
    int err = 0;

    if (vmin <= 1) vmin = 1;

    for (i=vmin; i<pdinfo->v; i++) {
	if (var_is_hidden(pdinfo, i)) {
	    nhid++;
	}
    }

    if (nhid > 0) {
	int *list = gretl_list_new(nhid);

	if (list == NULL) {
	    err = 1;
	} else {
	    int j = 1;

	    for (i=vmin; i<pdinfo->v; i++) {
		if (var_is_hidden(pdinfo, i)) {
		    list[j++] = i;
		}
	    }	    
	    err = dataset_drop_listed_variables(list, pZ, pdinfo, 
						NULL, NULL);
	    free(list);
	}
    }

    return err;
}

/* intended for use with newly imported data: trash any 
   series that contain nothing but NAs
*/

int maybe_prune_dataset (double ***pZ, DATAINFO **ppdinfo, void *p)
{
    DATAINFO *pdinfo = *ppdinfo;
    int allmiss, prune = 0, err = 0;
    int i, t;

    for (i=1; i<pdinfo->v; i++) {
	allmiss = 1;
	for (t=0; t<pdinfo->n; t++) {
	    if (!na((*pZ)[i][t])) {
		allmiss = 0;
		break;
	    }
	}
	if (allmiss) {
	    prune = 1;
	    break;
	}
    }

    if (prune) {
	char *mask = calloc(pdinfo->v, 1);
	double **newZ = NULL;
	DATAINFO *newinfo = NULL;
	int ndrop = 0;

	if (mask == NULL) {
	    return E_ALLOC;
	}

	for (i=1; i<pdinfo->v; i++) {
	    allmiss = 1;
	    for (t=0; t<pdinfo->n; t++) {
		if (!na((*pZ)[i][t])) {
		    allmiss = 0;
		    break;
		}
	    }
	    if (allmiss) {
		mask[i] = 1;
		ndrop++;
	    }
	}

	newinfo = datainfo_new();
	if (newinfo == NULL) {
	    err = E_ALLOC;
	} else {
	    newinfo->v = pdinfo->v - ndrop;
	    newinfo->n = pdinfo->n;
	    err = start_new_Z(&newZ, newinfo, 0);
	}

	if (!err) {
	    gretl_string_table *st = (gretl_string_table *) p;
	    size_t ssize = pdinfo->n * sizeof **newZ;
	    int k = 1;

	    for (i=1; i<pdinfo->v; i++) {
		if (!mask[i]) {
		    memcpy(newZ[k], (*pZ)[i], ssize);
		    strcpy(newinfo->varname[k], pdinfo->varname[i]);
		    strcpy(VARLABEL(newinfo, k), VARLABEL(pdinfo, i));
		    if (st != NULL && k < i) {
			gretl_string_table_reset_column_id(st, i, k);
		    }
		    k++;
		}
	    }

	    destroy_dataset(*pZ, pdinfo);
	    *pZ = newZ;
	    *ppdinfo = newinfo;

	    fprintf(stderr, "pruned dataset to %d variables\n", newinfo->v);
	}

	free(mask);
    }

    return err;
}

/* apparatus for sorting entire dataset */

typedef struct spoint_t_ spoint_t;

struct spoint_t_ {
    int obsnum;
    double val;
};

static int compare_vals_up (const void *a, const void *b)
{
    const spoint_t *pa = (const spoint_t *) a;
    const spoint_t *pb = (const spoint_t *) b;
     
    return (pa->val > pb->val) - (pa->val < pb->val);
}

static int compare_vals_down (const void *a, const void *b)
{
    const spoint_t *pa = (const spoint_t *) a;
    const spoint_t *pb = (const spoint_t *) b;
     
    return (pa->val < pb->val) - (pa->val > pb->val);
}

int dataset_sort_by (int v, double **Z, DATAINFO *pdinfo, gretlopt opt)
{
    spoint_t *sv = NULL;
    double *x = NULL;
    char **S = NULL;
    int i, t;
    int err = 0;

    if (v < 1 || v >= pdinfo->v) {
	return E_DATA;
    }

    sv = malloc(pdinfo->n * sizeof *sv);
    if (sv == NULL) {
	return E_ALLOC;
    }

    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) {
	free(sv);
	return E_ALLOC;
    }    

    if (pdinfo->S != NULL) {
	S = strings_array_new_with_length(pdinfo->n, OBSLEN);
	if (S == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    for (t=0; t<pdinfo->n; t++) {
	sv[t].obsnum = t;
	sv[t].val = Z[v][t];
    }

    if (opt & OPT_D) {
	/* descending */
	qsort(sv, pdinfo->n, sizeof *sv, compare_vals_down);
    } else {
	qsort(sv, pdinfo->n, sizeof *sv, compare_vals_up);
    }

    for (i=1; i<pdinfo->v; i++) {
	for (t=0; t<pdinfo->n; t++) {
	    x[t] = Z[i][sv[t].obsnum];
	}
	for (t=0; t<pdinfo->n; t++) {
	    Z[i][t] = x[t];
	}
    }

    if (S != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    strcpy(S[t], pdinfo->S[sv[t].obsnum]);
	}
	free_strings_array(pdinfo->S, pdinfo->n);
	pdinfo->S = S;
    }

 bailout:

    free(sv);
    free(x);

    return err;
}

static int dataset_sort (const char *s, double **Z, DATAINFO *pdinfo,
			 gretlopt opt)
{
    char fmt[10];
    char vname[VNAMELEN];
    int v;

    if (dataset_is_time_series(pdinfo) ||
	dataset_is_panel(pdinfo)) {
	gretl_errmsg_set("You can only do this with undated data");
	return E_DATA;
    }

    s += strspn(s, " \t");

    sprintf(fmt, "%%%ds", VNAMELEN - 1);

    if (sscanf(s, fmt, vname) != 1) {
	return E_DATA;
    }

    v = series_index(pdinfo, vname);
    if (v == pdinfo->v) {
	return E_UNKVAR;
    }

    if (v < 1) {
	return E_DATA;
    }
    
    return dataset_sort_by(v, Z, pdinfo, opt);
}

/**
 * dataset_drop_last_variables:
 * @delvars: number of variables to be dropped.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Deletes from the dataset the number @delvars of variables 
 * that were added most recently (that have the highest ID numbers).
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_drop_last_variables (int delvars, double ***pZ, DATAINFO *pdinfo)
{
    int newv = pdinfo->v - delvars;
    int i, err = 0;

    if (delvars <= 0) {
	return 0;
    }

#if FULLDEBUG
    fprintf(stderr, "*** dataset_drop_last_variables: dropping %d, newv = %d\n",
	    delvars, newv);
#endif

    if (newv < 1) {
	fprintf(stderr, "dataset_drop_last_vars: pdinfo->v = %d, delvars = %d "
		" -> newv = %d\n (pdinfo = %p)\n", pdinfo->v, delvars,
		newv, (void *) pdinfo);
	return E_DATA;
    }

#if FULLDEBUG
    for (i=0; i<pdinfo->v; i++) {
	if ((*pZ)[i] == NULL) {
	    fprintf(stderr, "var %d (%s, level %d, val = NULL) %s\n", 
		    i, pdinfo->varname[i], STACK_LEVEL(pdinfo, i), 
		    (i >= newv)? "deleting" : "");
	} else {
	    fprintf(stderr, "var %d (%s, level %d, val[0] = %g) %s\n", 
		    i, pdinfo->varname[i], STACK_LEVEL(pdinfo, i), 
		    (*pZ)[i][0], (i >= newv)? "deleting" : "");
	}
    }
#endif

    for (i=newv; i<pdinfo->v; i++) {
	free(pdinfo->varname[i]);
	free_varinfo(pdinfo, i);
	free((*pZ)[i]);
	(*pZ)[i] = NULL;
    }

    err = shrink_dataset_to_size(pZ, pdinfo, newv, DROP_NORMAL);

    if (!err) {
	err = gretl_lists_revise(NULL, newv);
    }

    if (!err && complex_subsampled()) {
	double ***fZ = fetch_full_Z();
	DATAINFO *fdinfo = fetch_full_datainfo();

	/* 
	   Context: we're deleting @delvars variables at the end of Z,
	   leaving @newv variables.  The dataset is currently
	   subsampled.

	   Question: should we delete any variables at the end of
	   fullZ to keep the two arrays in sync?

	   If @newv < fdinfo->v, this must mean that at least some of
	   the extra vars we're deleting from the current sub-sampled
	   Z have already been synced to fullZ, so we should do the
	   deletion from fullZ.
	*/

	if (newv < fdinfo->v) { 
#if FULLDEBUG
	    fprintf(stderr, "prior fdinfo->v = %d: shrinking fullZ to %d vars\n", 
		    fdinfo->v, newv);
#endif
	    for (i=newv; i<fdinfo->v; i++) {
		free((*fZ)[i]);
		(*fZ)[i] = NULL;
	    }
	    err = shrink_dataset_to_size(fZ, fdinfo, newv, DROP_SPECIAL);
	    reset_full_Z(fZ);
	} 
    }

    return err;
}

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

static int get_stack_param_val (const char *s, const double **Z,
				const DATAINFO *pdinfo)
{
    int val = -1;

    if (isdigit(*s)) {
	val = atoi(s);
    } else {
	char vname[VNAMELEN];
	int i, len = strcspn(s, " -");

	if (len > VNAMELEN - 1) {
	    len = VNAMELEN - 1;
	}
	*vname = '\0';
	strncat(vname, s, len);
	if (gretl_is_scalar(vname)) {
	    val = gretl_scalar_get_value(vname);
	} else {
	    i = series_index(pdinfo, vname);
	    if (i < pdinfo->v) {
		val = (int) Z[i][0];
	    }
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
	    *err = E_PARSE;
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
	    *err = E_PARSE;
	} else {
	    len = get_stack_param_val(p + 9, Z, pdinfo);
	    if (len < 0 || len > pdinfo->n) {
		*err = E_DATA;
	    }
	}
    }

    return len;
}

/* Apparatus for stacking variables (e.g. in case of panel
   data that were read in "wrongly").
*/

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

/**
 * dataset_stack_variables:
 * @vname: name for new variable, to be produced by stacking.
 * @line: instructions for stacking existing variables.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @prn: printing apparatus.
 *
 * Really for internal use.  Don't worry about it.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_stack_variables (const char *vname, const char *line,
			     double ***pZ, DATAINFO *pdinfo, 
			     PRN *prn)
{
    char vn1[VNAMELEN], vn2[VNAMELEN];
    char format[16];
    char *p, *s = NULL, *scpy = NULL;
    int *vnum = NULL;
    double *bigx = NULL;
    int i, v1 = 0, v2 = 0, nv = 0;
    int maxok, offset;
    int oldn, bign, genv;
    int err = 0;

    /* copy full "stack(...)" spec for later reference */
    scpy = gretl_strdup(line);
    if (scpy == NULL) {
	return E_ALLOC;
    }

    line += 6; /* skip "stack(" */
    if (*line == ',') {
	free(scpy);
	return E_PARSE;
    }

    /* copy active portion of line (we use strtok below) */
    s = gretl_strdup(line);
    if (s == NULL) {
	free(scpy);
	return E_ALLOC;
    }

    p = strrchr(s, ')');
    if (p == NULL) {
	err = E_PARSE;
	goto bailout;
    }

    /* end of active portion of line */
    *p = '\0';

    genv = series_index(pdinfo, vname);

    /* do we have a range of vars? */
    sprintf(format, "%%%d[^.]..%%%ds", VNAMELEN-1, VNAMELEN-1);
    if (sscanf(s, format, vn1, vn2) == 2) {
	if (isdigit(*vn1) && isdigit(*vn2)) {
	    v1 = atoi(vn1);
	    v2 = atoi(vn2);
	} else {
	    v1 = series_index(pdinfo, vn1);
	    v2 = series_index(pdinfo, vn2);
	}
	if (v1 >= 0 && v2 > v1 && v2 < pdinfo->v) {
	    nv = v2 - v1 + 1;
	} else {
	    fputs("stack vars: range is invalid\n", stderr);
	    err = E_DATA;
	}
    } else {
	/* or a comma-separated list of vars? */
	char *p = s;

	while (*p) {
	    if (*p == ',') nv++;
	    p++;
	}
	nv++;

	if (nv < 2) {
	    return E_PARSE;
	}

	vnum = malloc(nv * sizeof *vnum);
	if (vnum == NULL) {
	    err = E_ALLOC;
	}

	for (i=0; i<nv && !err; i++) {
	    p = strtok((i == 0)? s : NULL, ",");
	    while (*p == ' ') p++;
	    if (isdigit(*p)) {
		v1 = atoi(p);
	    } else {
		v1 = series_index(pdinfo, p);
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

#if PDEBUG
    fprintf(stderr, "offset = %d, maxok = %d\n", offset, maxok);
#endif

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

	    ok = pdinfo->n - missing_tail((*pZ)[j], pdinfo->n);

	    if (ok > maxok) {
		maxok = ok;
	    }
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

#if PDEBUG
    fprintf(stderr, "bign = %d, allocating bigx (oldn = %d)\n", bign, pdinfo->n);
#endif

    /* allocate stacked series */
    bigx = malloc(bign * sizeof *bigx);
    if (bigx == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* extend length of all series? */
    oldn = pdinfo->n;
    if (bign > oldn) {
	err = dataset_add_observations(bign - oldn, pZ, pdinfo, OPT_NONE);
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
	    bigx[bigt] = (*pZ)[j][t];
	    if (pdinfo->S != NULL && bigt != t) {
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
	err = dataset_add_allocated_series(bigx, pZ, pdinfo);
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
	strcpy(pdinfo->varname[genv], vname);
	make_stack_label(VARLABEL(pdinfo, genv), scpy);
	pprintf(prn, "%s %s %s (ID %d)\n", 
		(genv == pdinfo->v - 1)? _("Generated") : _("Replaced"),
		_("series"), vname, genv);
    }

 bailout:

    free(vnum);
    free(s);
    free(scpy);

    return err;
}

static int found_log_parent (const char *s, char *targ)
{
    int len = gretl_namechar_spn(s);

    if (len < VNAMELEN && s[len] == ')') {
	sscanf(s, "%15[^)]", targ);
	return 1;
    }

    return 0;
}

/**
 * is_log_variable:
 * @i: ID number of variable.
 * @pdinfo: dataset information.
 * @parent: location to which to write the name of the
 * "parent" variable if any.
 *
 * Tries to determine if the variable with ID number @i is
 * the logarithm of some other variable.
 *
 * Returns: 1 if variable @i appears to be a log, else 0.
 */

int is_log_variable (int i, const DATAINFO *pdinfo, char *parent)
{
    const char *s = VARLABEL(pdinfo, i);

    *parent = '\0';

    if (s != NULL && *s != '\0') {
	if (sscanf(s, "= log of %15s", parent) == 1) {
	    return 1;
	} else if (!strncmp(s, "log(", 4)) {
	    return found_log_parent(s + 4, parent);
	} else {
	    s += strcspn(s, "=");
	    if (!strncmp(s, "=log(", 5)) {
		return found_log_parent(s + 5, parent);
	    }
	}
    }

    return 0;
}

/**
 * set_var_discrete:
 * @pdinfo: pointer to data information struct.
 * @i: index number of variable.
 * @s: non-zero to mark variable as discrete, zero to 
 * mark as not discrete.
 *
 * Mark a variable as being discrete or not.
 */

void set_var_discrete (DATAINFO *pdinfo, int i, int s) 
{
    if (i > 0 && i < pdinfo->v) {
	int flags = pdinfo->varinfo[i]->flags;

	if (s && !(flags & VAR_DISCRETE)) {
	    pdinfo->varinfo[i]->flags |= VAR_DISCRETE;
	    set_dataset_is_changed();
	} else if (!s && (flags & VAR_DISCRETE)) {
	    pdinfo->varinfo[i]->flags &= ~VAR_DISCRETE;
	    set_dataset_is_changed();
	}
    }
}

/**
 * set_var_hidden:
 * @pdinfo: pointer to data information struct.
 * @i: index number of variable.
 *
 * Mark a variable as being "hidden" (an automatically
 * generated variable that will not be shown in the main
 * GUI window).
 */

void set_var_hidden (DATAINFO *pdinfo, int i) 
{
    if (i > 0 && i < pdinfo->v) {
	pdinfo->varinfo[i]->flags |= VAR_HIDDEN;
    }
}

/**
 * var_set_linewidth:
 * @pdinfo: pointer to data information struct.
 * @i: index number of variable.
 * @w: with of plot line.
 *
 * Set the line width for use when this variable is displayed
 * in a line graph.
 */

void var_set_linewidth (DATAINFO *pdinfo, int i, int w) 
{
    if (w >= 1 && w <= 32 && i > 0 && i < pdinfo->v) {
	pdinfo->varinfo[i]->line_width = w;
    }
}

/**
 * var_get_linewidth:
 * @pdinfo: pointer to data information struct.
 * @i: index number of variable.
 *
 * Returns: the line width set for use when graphing
 * variable @i.
 */

int var_get_linewidth (const DATAINFO *pdinfo, int i) 
{
    if (i > 0 && i < pdinfo->v) {
	return pdinfo->varinfo[i]->line_width;
    } else {
	return 0;
    }
}

int var_set_description (DATAINFO *pdinfo, int i,
			 const char *s) 
{
    char *targ = pdinfo->varinfo[i]->label;

    if (strcmp(targ, s)) {
	*targ = 0;
	strncat(targ, s, MAXLABEL - 1);
	set_dataset_is_changed();
    }

    return 0;
}

int var_set_display_name (DATAINFO *pdinfo, int i,
			  const char *s) 
{
    char *targ = pdinfo->varinfo[i]->display_name;

    if (strcmp(targ, s)) {
	*targ = 0;
	strncat(targ, s, MAXDISP - 1);
	set_dataset_is_changed();
    }

    return 0;
}

int var_set_compact_method (DATAINFO *pdinfo, int i,
			    int method) 
{
    int orig = COMPACT_METHOD(pdinfo, i);

    if (method != orig) {
	COMPACT_METHOD(pdinfo, i) = method;
	set_dataset_is_changed();
    }

    return 0;
}

const char *var_get_graph_name (const DATAINFO *pdinfo, int i)
{
    const char *ret = pdinfo->varname[i];

    if (pdinfo->varinfo != NULL && pdinfo->varinfo[i] != NULL) {
	ret = pdinfo->varinfo[i]->display_name;
	if (ret[0] == '\0') {
	    ret = pdinfo->varname[i];
	}
    } 

    return ret;
}

static int add_obs (int n, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int err = 0;

    if (complex_subsampled()) {
	pprintf(prn, _("The data set is currently sub-sampled.\n"));
	err = E_DATA;
    } else if (n <= 0) {
	err = E_PARSE;
    } else {
	err = dataset_add_observations(n, pZ, pdinfo, OPT_A);
	if (!err) {
	    pprintf(prn, _("Dataset extended by %d observations"), n);
	    pputc(prn, '\n');
	}
    }

    return err;
}

int dataset_op_from_string (const char *s)
{
    int op = DS_NONE;

    if (!strcmp(s, "addobs")) {
	op = DS_ADDOBS;
    } else if (!strcmp(s, "compact")) {
	op = DS_COMPACT;
    } else if (!strcmp(s, "expand")) {
	op = DS_EXPAND;
    } else if (!strcmp(s, "transpose")) {
	op = DS_TRANSPOSE;
    } else if (!strcmp(s, "delete")) {
	op = DS_DELETE;
    } else if (!strcmp(s, "keep")) {
	op = DS_KEEP;
    } else if (!strcmp(s, "sortby")) {
	op = DS_SORTBY;
    } else if (!strcmp(s, "dsortby")) {
	op = DS_DSORTBY;
    } else if (!strcmp(s, "resample")) {
	op = DS_RESAMPLE;
    } else if (!strcmp(s, "restore")) {
	op = DS_RESTORE;
    } else if (!strcmp(s, "clear")) {
	op = DS_CLEAR;
    }

    return op;
}

static int dataset_int_param (const char **ps, int op, 
			      double **Z, DATAINFO *pdinfo, 
			      int *err)
{
    const char *s = *ps;
    char test[32];
    int k = 0;

    if ((op == DS_COMPACT || op == DS_EXPAND) &&
	!dataset_is_time_series(pdinfo)) {
	*err = E_PDWRONG;
	return 0;
    }

    *test = '\0';
    sscanf(s, "%31s", test);
    *ps += strlen(test);

    k = gretl_int_from_string(test, err);
    if (*err) {
	return 0;
    }

    if (k <= 0 || (op == DS_RESAMPLE && k < 1)) {
	*err = E_DATA;
    } else if (op == DS_COMPACT) {
	int ok = 0;

	if (pdinfo->pd == 12 && (k == 4 || k == 1)) {
	    ok = 1;
	} else if (pdinfo->pd == 4 && k == 1) {
	    ok = 1;
	} else if (pdinfo->pd == 52 && k == 12) {
	    ok = 1;
	} else if (dated_daily_data(pdinfo) && (k == 52 || k == 12)) {
	    ok = 1;
	}

	if (!ok) {
	    *err = E_PDWRONG;
	    gretl_errmsg_set("This conversion is not supported");
	}
    } else if (op == DS_EXPAND) {
	int ok = 0;

	if (pdinfo->pd == 1 && (k == 4 || k == 12)) {
	    ok = 1;
	} else if (pdinfo->pd == 4 && k == 12) {
	    ok = 1;
	} 

	if (!ok) {
	    *err = E_PDWRONG;
	    gretl_errmsg_set("This conversion is not supported");
	}
    }

    return k;
}

static int compact_data_set_wrapper (const char *s, double ***pZ, 
				     DATAINFO *pdinfo, int k)
{
    CompactMethod method = COMPACT_AVG;

    s += strspn(s, " ");

    if (!strcmp(s, "sum")) {
	method = COMPACT_SUM;
    } else if (!strcmp(s, "first") || !strcmp(s, "sop")) {
	method = COMPACT_SOP;
    } else if (!strcmp(s, "last") || !strcmp(s, "eop")) {
	method = COMPACT_EOP;
    } else if (!strcmp(s, "avg") || !strcmp(s, "average")) {
	method = COMPACT_AVG;
    } else if (*s != '\0') {
	return E_PARSE;
    }

    return compact_data_set(pZ, pdinfo, k, method, 0, 0);
}

static unsigned int resample_seed;

unsigned int get_resampling_seed (void)
{
    return resample_seed;
}

/* resample the dataset by observation, with replacement */

int dataset_resample (int n, unsigned int seed,
		      double ***pZ, DATAINFO *pdinfo)
{
    double **RZ = NULL;
    DATAINFO *rinfo = NULL;
    char **S = NULL;
    int T = sample_size(pdinfo);
    int v = pdinfo->v;
    int i, j, s, t;
    int err = 0;

    if (v < 2) {
	return E_DATA;
    }

    rinfo = datainfo_new();
    if (rinfo == NULL) {
	return E_ALLOC;
    }

    RZ = malloc(v * sizeof *RZ);
    if (RZ == NULL) {
	free(rinfo);
	return E_ALLOC;
    }

    for (i=0; i<v; i++) {
	RZ[i] = NULL;
    }

    rinfo->v = v;

    j = 0;
    for (i=0; i<pdinfo->v && !err; i++) {
	RZ[j] = malloc(n * sizeof **RZ);
	if (RZ[j] != NULL && i == 0) {
	    for (t=0; t<n; t++) {
		RZ[j][t] = 1.0;
	    }
	}
	if (RZ[j] == NULL) {
	    err = E_ALLOC;
	}
	j++;
    }

    if (err) {
	goto bailout;
    }

    if (pdinfo->markers == REGULAR_MARKERS) {
	S = strings_array_new_with_length(n, OBSLEN);
    }

    if (seed > 0) {
	resample_seed = seed;
	gretl_rand_set_seed(seed);
    } else {
	resample_seed = gretl_rand_get_seed();
    }

    for (t=0; t<n; t++) {
	s = gretl_rand_int_max(T) + pdinfo->t1;
	j = 1;
	for (i=1; i<pdinfo->v; i++) {
	    RZ[j][t] = (*pZ)[i][s];
	    j++;
	}
	if (S != NULL) {
	    strcpy(S[t], pdinfo->S[s]);
	}
    }

    if (S != NULL) {
	rinfo->S = S;
	rinfo->markers = REGULAR_MARKERS;
    } 

    rinfo->varname = pdinfo->varname;
    rinfo->varinfo = pdinfo->varinfo;
    rinfo->descrip = pdinfo->descrip;

    rinfo->n = n;
    rinfo->t1 = 0;
    rinfo->t2 = n - 1;
    dataset_obs_info_default(rinfo);

    set_dataset_resampled(rinfo);

 bailout:

    if (err) {
	free_Z(RZ, rinfo);
	clear_datainfo(rinfo, CLEAR_SUBSAMPLE);
	free(rinfo);
    } else {
	backup_full_dataset(*pZ, pdinfo);
	*pZ = RZ;
	*pdinfo = *rinfo;
	free(rinfo);
    }

    return err;
}

int modify_dataset (int op, const int *list, const char *s, 
		    double ***pZ, DATAINFO *pdinfo, 
		    PRN *prn)
{
    static int resampled;
    int k = 0, err = 0;

    if (pZ == NULL || *pZ == NULL || pdinfo == NULL) {
	return E_NODATA;
    }

    if (op == DS_CLEAR) {
	/* must be handled by the calling program */
	return E_NOTIMP;
    }

    if (gretl_function_depth() > 0) {
	gretl_errmsg_set(_("The 'dataset' command is not available within functions"));
	return 1;
    }

    if (op != DS_RESAMPLE && op != DS_RESTORE && gretl_looping()) {
	pputs(prn, _("Sorry, this command is not available in loop mode\n"));
	return 1;
    }

    if (op == DS_RESAMPLE && resampled) {
	/* repeated "resample": implicitly restore first */
	err = restore_full_sample(pZ, pdinfo, NULL);
	if (err) {
	    return err;
	} else {
	    resampled = 0;
	}
    }

    if (op != DS_RESTORE && complex_subsampled()) {
	gretl_errmsg_set(_("The data set is currently sub-sampled"));
	return 1;
    }

    if (op == DS_ADDOBS || op == DS_COMPACT || 
	op == DS_EXPAND || op == DS_RESAMPLE) {
	k = dataset_int_param(&s, op, *pZ, pdinfo, &err);
	if (err) {
	    return err;
	}
    }

    if (op == DS_ADDOBS) {
	err = add_obs(k, pZ, pdinfo, prn);
    } else if (op == DS_COMPACT) {
	err = compact_data_set_wrapper(s, pZ, pdinfo, k);
    } else if (op == DS_EXPAND) {
	err = expand_data_set(pZ, pdinfo, k);
    } else if (op == DS_TRANSPOSE) {
	err = transpose_data(pZ, pdinfo);
    } else if (op == DS_SORTBY) {
	err = dataset_sort(s, *pZ, pdinfo, OPT_NONE);
    } else if (op == DS_DSORTBY) {
	err = dataset_sort(s, *pZ, pdinfo, OPT_D);
    } else if (op == DS_RESAMPLE) {
	err = dataset_resample(k, 0, pZ, pdinfo);
	if (!err) {
	    resampled = 1;
	}
    } else if (op == DS_RESTORE) {
	if (resampled) {
	    err = restore_full_sample(pZ, pdinfo, NULL);
	    resampled = 0;
	} else {
	    pprintf(prn, _("dataset restore: dataset is not resampled\n"));
	    err = E_DATA;
	}
    } else if (op == DS_DELETE) {
	pprintf(prn, "dataset delete: not ready yet\n");
    } else if (op == DS_KEEP) {
	pprintf(prn, "dataset keep: not ready yet\n");
    } else {
	err = E_PARSE;
    }

    return err;
}

int dataset_get_structure (const DATAINFO *pdinfo)
{
    if (pdinfo == NULL || pdinfo->n == 0) {
	return DATA_NONE;
    } else if (dataset_is_panel(pdinfo)) {
	return DATA_PANEL;
    } else if (dataset_is_time_series(pdinfo)) {
	return DATA_TS;
    } else {
	return DATA_XSECT;
    }
}

/**
 * dataset_purge_missing_rows:
 * @Z: data array.
 * @pdinfo: pointer to data information struct.
 * 
 * Removes empty rows from the dataset -- that is, observations
 * at which there are no non-missing values.  This is intended
 * for daily data only.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_purge_missing_rows (double **Z, DATAINFO *pdinfo)
{
    int new_n, missrow, totmiss = 0;
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    char **S = NULL;
    double *Zi = NULL;
    size_t sz;
    int i, t, s;
    int err = 0;

    for (t=0; t<pdinfo->n; t++) {
	missrow = 1;
	for (i=1; i<pdinfo->v; i++) {
	    if (!na(Z[i][t])) {
		missrow = 0;
		break;
	    }
	}
	if (missrow) {
	    totmiss++;
	    if (t < pdinfo->t1) {
		t1--;
	    }
	    if (t < pdinfo->t2) {
		t2--;
	    }
	}
    }

    if (totmiss == 0) {
	/* no-op */
	return 0;
    }

    if (dated_daily_data(pdinfo) && pdinfo->S == NULL) {
	err = dataset_allocate_obs_markers(pdinfo);
	if (!err) {
	    for (t=0; t<pdinfo->n; t++) {
		calendar_date_string(pdinfo->S[t], t, pdinfo);
	    }
	}
    }

    for (t=0; t<pdinfo->n; t++) {
	missrow = 1;
	for (i=1; i<pdinfo->v; i++) {
	    if (!na(Z[i][t])) {
		missrow = 0;
		break;
	    }
	}
	if (missrow) {
	    sz = (pdinfo->n - t) * sizeof **Z;
	    for (i=1; i<pdinfo->v; i++) {
		memmove(Z[i] + t, Z[i] + t + 1, sz);
	    }
	    if (pdinfo->S != NULL) {
		free(pdinfo->S[t]);
		for (s=t; s<pdinfo->n - 1; s++) {
		    pdinfo->S[s] = pdinfo->S[s+1];
		}
	    }
	}
    } 

    new_n = pdinfo->n - totmiss;

    for (i=1; i<pdinfo->v; i++) {
	Zi = realloc(Z[i], new_n * sizeof *Zi);
	if (Zi == NULL) {
	    err = E_ALLOC;
	} else {
	    Z[i] = Zi;
	}
    }

    if (!err && pdinfo->S != NULL) {
	S = realloc(pdinfo->S, new_n * sizeof *S);
	if (S == NULL) {
	    err = E_ALLOC;
	} else {
	    pdinfo->S = S;
	    if (dated_daily_data(pdinfo)) {
		strcpy(pdinfo->stobs, pdinfo->S[0]);
		strcpy(pdinfo->endobs, pdinfo->S[new_n-1]);
		pdinfo->sd0 = get_epoch_day(pdinfo->stobs);
	    }
	}
    }

    pdinfo->n = new_n;
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    return err;
}
