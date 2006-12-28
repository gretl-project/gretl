/*
 *  Copyright (c) 2005 by Allin Cottrell
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

#define DDEBUG 0

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
	free(pdinfo->submask);
	pdinfo->submask = NULL;
    }

    pdinfo->submode = 0;

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

	/* added Sat Dec  4 12:19:26 EST 2004 */
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
 * the unit or group and time-period, respectively,  of each 
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
	    strcpy(gretl_errmsg, "Panel index information is corrupted");
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
    vinfo->flags = 0;
    vinfo->compact_method = COMPACT_NONE;
    vinfo->stack_level = 0;
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
    targ->flags = src->flags;
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
    dinfo->submode = 0;

    dinfo->data = NULL;

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

    pdinfo->markers = (unsigned char) markers;
    if (pdinfo->markers) {
	if (dataset_allocate_obs_markers(pdinfo)) {
	    free_datainfo(pdinfo);
	    return NULL;
	}
    } 

    dataset_obs_info_default(pdinfo);
    pdinfo->descrip = NULL;

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
    pdinfo->data = NULL;
    pdinfo->paninfo = NULL;
    pdinfo->submask = NULL;
    pdinfo->submode = 0;
    
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
	if (var_is_series(pdinfo, i)) {
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
	if (var_is_series(pdinfo, i)) {
	    x = realloc((*pZ)[i], newn * sizeof *x);
	    if (x == NULL) {
		return E_ALLOC;
	    }
	    (*pZ)[i] = x;
	}
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
	    if (var_is_series(pdinfo, i)) {
		memmove(Z[i], Z[i] + head, mvsize);
	    }
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
 * dataset_add_scalars:
 * @n: number of scalars to add.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Allocates space for @n new scalar members of the dataset.
 * The added variables are initialized to zero.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_scalars (int n, double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    int v = pdinfo->v;
    int i, err = 0;

    newZ = realloc(*pZ, (v + n) * sizeof *newZ);  

    if (newZ == NULL) {
	err = E_ALLOC;
    } else {
	*pZ = newZ;
    }

    if (!err) {
	for (i=0; i<n; i++) {
	    newZ[v+i] = NULL;
	}
	for (i=0; i<n; i++) {
	    newZ[v+i] = malloc(sizeof **newZ);
	    if (newZ[v+i] == NULL) {
		err = E_ALLOC;
		break;
	    }
	    newZ[v+i][0] = 0.0;
	}
    }

    if (!err) {
	err = dataset_expand_varinfo(n, pdinfo);
    }

    if (!err) {
	for (i=0; i<n; i++) {
	    set_var_scalar(pdinfo, v + i, 1);
	}
    }

    return err;
}

/**
 * dataset_add_scalar:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Allocates space for a new scalar member of the dataset.
 * The added variable is initialized to zero.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo)
{
    return dataset_add_scalars(1, pZ, pdinfo);
}

/**
 * dataset_add_scalar_as:
 * @x: scalar value.
 * @newname: name to give the new variable.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Adds to the dataset a new scalar with name @newname and
 * value given by @x.  The new variable is added at one
 * level "deeper" (in terms of function execution) than the
 * current level.  This is for use with user-defined functions.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_scalar_as (double x, const char *newname,
			   double ***pZ, DATAINFO *pdinfo)
{
    int v, err = 0;

    if (pdinfo->varinfo == NULL) {
	strcpy(gretl_errmsg, _("Please open a data file first"));
	return 1;
    }

    err = dataset_add_scalar(pZ, pdinfo);

    if (!err) {
	v = pdinfo->v - 1;
	(*pZ)[v][0] = x;
	strcpy(pdinfo->varname[v], newname);
	STACK_LEVEL(pdinfo, v) += 1;
    }

    return err;
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
	strcpy(gretl_errmsg, _("Please open a data file first"));
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

    if (var_is_series(pdinfo, v)) {
	err = real_add_series(1, NULL, pZ, pdinfo);
    } else {
	err = dataset_add_scalar(pZ, pdinfo);
    }

    if (!err) {
	int vnew = pdinfo->v - 1;

	if (var_is_series(pdinfo, v)) {
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[vnew][t] = (*pZ)[v][t];
	    }
	} else {
	    (*pZ)[vnew][0] = (*pZ)[v][0];
	}
	strcpy(pdinfo->varname[vnew], newname);
	STACK_LEVEL(pdinfo, vnew) += 1;
	 /* FIXME other varinfo stuff?? */
    }

    return err;
}

enum {
    DROP_NORMAL,
    DROP_SPECIAL
};

/* DROP_SPECIAL is used when deleting variables from the "full" shadow
   of a sub-sampled dataset: in that case we don't mess with the
   pointer elements of the DATAINFO struct, because these are shared
   between the full and sub-sampled datasets.
*/

static int 
shrink_dataset_to_size (double ***pZ, DATAINFO *pdinfo, int nv, int drop)
{
    double **newZ;

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

int overwrite_err (const DATAINFO *pdinfo, int v)
{
    sprintf(gretl_errmsg, "The variable %s is read-only", 
	    pdinfo->varname[v]);

    return E_DATA;
}

static int real_drop_listed_vars (const int *list, double ***pZ, 
				  DATAINFO *pdinfo, int *renumber,
				  int drop)
{
    double **Z = *pZ;
    int oldv = pdinfo->v, vmax = pdinfo->v;
    int delmin = oldv;
    int i, v, ndel = 0; 
    int err;

    if (list == NULL) {
	/* no-op */
	return 0;
    }

    if (renumber != NULL) {
	*renumber = 0;
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
	    if (var_is_const(pdinfo, v)) {
		return overwrite_err(pdinfo, v);
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

#if DDEBUG
    fprintf(stderr, "real_drop_listed_variables: lowest ID of deleted var"
	    " = %d\n", delmin);
#endif

#if 0 /* not yet: stuff for recreating saved lists in the face of the
	 renumbering of variables (see below) */
    if (1) {
	int ren = oldv - ndel - delmin;

	if (ren > 0) {
	    /* if the vars "to be renumbered" are actually hidden, we're OK? */
	    for (i=delmin; i<vmax; i++) {
		if (var_is_hidden(pdinfo, i) && !in_gretl_list(list, i)) {
		    ren--;
		}
	    }
	    if (ren > 0) {
		fprintf(stderr, "drop vars: renumbering needed (ren=%d)\n", ren);
	    }
	}
    }
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
			if (!var_is_hidden(pdinfo, i + gap)) {
			    if (renumber != NULL) {
				*renumber = 1;
			    }
			}
			pdinfo->varname[i] = pdinfo->varname[i + gap];
			pdinfo->varinfo[i] = pdinfo->varinfo[i + gap];
		    }
		    Z[i] = Z[i + gap];
		}		    
	    } else {
		/* deleting all subsequent vars */
		break;
	    }
	}
    }

    err = shrink_dataset_to_size(pZ, pdinfo, oldv - ndel, drop);

    if (!err) {
	/* for now we'll just remove from saved lists any vars
	   that have been deleted OR renumbered */
	gretl_lists_prune(delmin);
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
				   int *renumber)
{
    int *dlist = NULL;
    int free_dlist = 0;
    int lastvar[2];
    int err = 0;

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
				DROP_NORMAL);

    if (!err && complex_subsampled()) {
	double ***fZ = fetch_full_Z();
	DATAINFO *fdinfo = fetch_full_datainfo();

	err = real_drop_listed_vars(dlist, fZ, fdinfo, NULL,
				    DROP_SPECIAL);
	reset_full_Z(fZ);
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

    return dataset_drop_listed_variables(list, pZ, pdinfo, NULL);
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
	    err = dataset_drop_listed_variables(list, pZ, pdinfo, NULL);
	    free(list);
	}
    }

    return err;
}

static int 
real_drop_last_vars (int delvars, double ***pZ, DATAINFO *pdinfo,
		     int drop)
{
    int i, v = pdinfo->v; 
    int newv = v - delvars;

#if DDEBUG
    fprintf(stderr, "real_drop_last_vars: dropping %d\n", delvars);
#endif

    for (i=newv; i<v; i++) {
	if (drop == DROP_NORMAL) {
	    free(pdinfo->varname[i]);
	    free_varinfo(pdinfo, i);
	}
	free((*pZ)[i]);
	(*pZ)[i] = NULL;
    }

    return shrink_dataset_to_size(pZ, pdinfo, newv, drop);
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
    int err;

    if (delvars <= 0) {
	return 0;
    }

    if (pdinfo->v <= 1) {
	return E_DATA;
    }

    err = real_drop_last_vars(delvars, pZ, pdinfo, DROP_NORMAL);

    if (!err && complex_subsampled()) {
	double ***fZ = fetch_full_Z();
	DATAINFO *fdinfo = fetch_full_datainfo();

	err = real_drop_last_vars(delvars, fZ, fdinfo, DROP_SPECIAL);
	reset_full_Z(fZ);
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
	return E_SYNTAX;
    }

    /* copy active portion of line (we use strtok below) */
    s = gretl_strdup(line);
    if (s == NULL) {
	free(scpy);
	return E_ALLOC;
    }

    p = strrchr(s, ')');
    if (p == NULL) {
	err = E_SYNTAX;
	goto bailout;
    }

    /* end of active portion of line */
    *p = '\0';

    genv = varindex(pdinfo, vname);

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
	/* or a comma-separated list of vars? */
	char *p = s;

	while (*p) {
	    if (*p == ',') nv++;
	    p++;
	}
	nv++;

	if (nv < 2) {
	    return E_SYNTAX;
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

	    if (var_is_series(pdinfo, j)) {
		ok = pdinfo->n - missing_tail((*pZ)[j], pdinfo->n);
	    } else {
		ok = 1;
	    }
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
	    if (var_is_series(pdinfo, j)) {
		bigx[bigt] = (*pZ)[j][t];
	    } else {
		bigx[bigt] = (*pZ)[j][0];
	    }
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
    int len = gretl_varchar_spn(s);

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
	if (s) {
	    pdinfo->varinfo[i]->flags |= VAR_DISCRETE;
	} else {
	    pdinfo->varinfo[i]->flags &= ~VAR_DISCRETE;
	}
    }
}

/**
 * set_var_scalar:
 * @pdinfo: pointer to data information struct.
 * @i: index number of variable.
 * @s: non-zero to mark variable as a scalar, zero to 
 * mark as not scalar (i.e. a series).
 *
 * Mark a variable as being a scalar or not.
 */

void set_var_scalar (DATAINFO *pdinfo, int i, int s) 
{
    if (i > 0 && i < pdinfo->v) {
	if (s) {
	    pdinfo->varinfo[i]->flags |= VAR_SCALAR;
	} else {
	    pdinfo->varinfo[i]->flags &= ~VAR_SCALAR;
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
