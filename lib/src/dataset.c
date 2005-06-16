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

/**
 * free_Z:
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Does a deep free on the data matrix.
 */

void free_Z (double **Z, DATAINFO *pdinfo)
{
    int i;

    if (Z == NULL || pdinfo == NULL) return;

    for (i=0; i<pdinfo->v; i++) {
	free(Z[i]);
    }
    free(Z);
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
	if (pdinfo->descrip) {
	    free(pdinfo->descrip);
	    pdinfo->descrip = NULL;
	}
	if (pdinfo->vector) {
	    free(pdinfo->vector);
	    pdinfo->vector = NULL;
	}

	maybe_free_full_dataset(pdinfo);

	/* added Sat Dec  4 12:19:26 EST 2004 */
	pdinfo->v = pdinfo->n = 0;
    } 
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

static char **allocate_obslen_strings (int n)
{
    char **S;
    int j, t;

    S = malloc(n * sizeof *S);
    if (S == NULL) return NULL;

    for (t=0; t<n; t++) {
	S[t] = malloc(OBSLEN);
	if (S[t] == NULL) {
	    for (j=0; j<t; j++) {
		free(S[j]);
	    }
	    free(S);
	    return NULL;
	}
	S[t][0] = '\0';
    }

    return S;
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
	S = allocate_obslen_strings(pdinfo->n);
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

static void gretl_varinfo_init (VARINFO *vinfo)
{
    vinfo->label[0] = '\0';
    vinfo->display_name[0] = '\0';
    vinfo->compact_method = COMPACT_NONE;
    vinfo->stack_level = 0;
    vinfo->sorted_markers = NULL;

    if (gretl_executing_function()) {
	vinfo->stack_level = gretl_function_stack_depth();
    } 
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
    
    pdinfo->varname = malloc(v * sizeof *pdinfo->varname);
    if (pdinfo->varname == NULL) {
	return E_ALLOC;
    }

    pdinfo->varinfo = malloc(v * sizeof *pdinfo->varinfo);
    if (pdinfo->varinfo == NULL) {
	free(pdinfo->varname);
	return E_ALLOC;
    }

    pdinfo->vector = malloc(v * sizeof *pdinfo->vector);
    if (pdinfo->vector == NULL) {
	free(pdinfo->varname);
	free(pdinfo->varinfo);
	return E_ALLOC;
    }

    for (i=0; i<v; i++) {
	pdinfo->varname[i] = malloc(VNAMELEN);
	if (pdinfo->varname[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(pdinfo->varname[j]);
	    }
	    free(pdinfo->varname);
	    pdinfo->varname = NULL;
	    err = E_ALLOC;
	    break;
	} else {
	    pdinfo->varname[i][0] = '\0';
	}

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

	pdinfo->vector[i] = 1;
    }

    if (!err) {
	strcpy(pdinfo->varname[0], "const");
	strcpy(VARLABEL(pdinfo, 0), _("auto-generated constant"));
    } else {
	free(pdinfo->vector);
	pdinfo->vector = NULL;
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

    dinfo->markers = NO_MARKERS;  
    dinfo->delim = ',';
    dinfo->decpoint = '.';

    dinfo->S = NULL;
    dinfo->descrip = NULL;
    dinfo->vector = NULL;
    dinfo->submask = NULL;
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
 * position 0 is initialized to all 1s; no other variables
 * are initialized.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int allocate_Z (double ***pZ, const DATAINFO *pdinfo)
{
    double **Z;
    int i, j, t;
    int err = 0;

    if (*pZ != NULL) {
	free(*pZ);
    }

    Z = malloc(pdinfo->v * sizeof *Z);

    if (Z == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<pdinfo->v && !err; i++) {
	    Z[i] = malloc(pdinfo->n * sizeof **Z);
	    if (Z[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(Z[j]);
		}
		free(Z);
		Z = NULL;
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	for (t=0; t<pdinfo->n; t++) {
	    Z[0][t] = 1.0; 
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
 * Returns: 0 on successful completion, 1 on error.
 */

int start_new_Z (double ***pZ, DATAINFO *pdinfo, int resample)
{
    if (allocate_Z(pZ, pdinfo)) return 1;

    pdinfo->t1 = 0; 
    pdinfo->t2 = pdinfo->n - 1;

    if (resample) {
	pdinfo->varname = NULL;
	pdinfo->varinfo = NULL;
    } else if (dataset_allocate_varnames(pdinfo)) {
	return 1;
    }

    pdinfo->S = NULL;
    pdinfo->markers = NO_MARKERS;
    pdinfo->delim = ',';
    pdinfo->descrip = NULL;
    pdinfo->data = NULL;
    pdinfo->submask = NULL;
    
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

static int real_periodic_dummy (const double *x, int n,
				int *pd, int *offset,
				int *trail)
{
    int onebak = 0;
    int gap = 0;
    int t, m = n - 1, ret = 1;

    *pd = -1;
    *offset = -1;
    *trail = 0;

    /* find number of trailing zeros */
    for (t=n-1; t>0; t--) {
	if (x[t] == 0.0) {
	    *trail += 1;
	} else {
	    if (x[t] == 1.0) {
		m = t;
	    } else {
		ret = 0;
	    }
	    break;
	}
    }

    /* check for dummyhood and periodicity */
    for (t=0; t<=m && ret; t++) {
	if (x[t] == 0.0) {
	    onebak = 0;
	    gap++;
	} else if (x[t] == 1.0) {
	    if (onebak) {
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
	    } else if (gap < *trail) {
		ret = 0;
	    }
	    gap = 0;
	    onebak = 1;
	} else {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * is_periodic_dummy:
 * @x: array to examine.
 * @n: number of elements in array.
 *
 * Returns: 1 if @x is a periodic dummy variable,
 * 0 otherwise.
 */

int is_periodic_dummy (const double *x, int n)
{
    int pd, offset, trail;

    return real_periodic_dummy(x, n, &pd, &offset, &trail);
}

/**
 * is_trend_variable:
 * @x: array to examine.
 * @n: number of elements in array.
 *
 * Returns: 1 if @x is a simple linear trend variable, with each
 * observation equal to the preceding observation plus 1, 
 * 0 otherwise.
 */

int is_trend_variable (const double *x, int n)
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

static void 
maybe_extend_trends (double **Z, const DATAINFO *pdinfo, int oldn)
{
    int i, t;

    for (i=1; i<pdinfo->v; i++) {
	if (is_trend_variable(Z[i], oldn)) {
	    for (t=oldn; t<pdinfo->n; t++) {
		Z[i][t] = Z[i][t-1] + 1.0;
	    }
	}
    }
}

static void 
maybe_extend_dummies (double **Z, const DATAINFO *pdinfo, int oldn)
{
    int pd, offset, trail;
    int i, t;

    for (i=0; i<pdinfo->v; i++) {
	if (real_periodic_dummy(Z[i], oldn, &pd, &offset, &trail)) {
	    int s = pd - offset;

	    for (t=oldn; t<pdinfo->n; t++) {
		Z[i][t] = (s % pd)? 0.0 : 1.0;
		s++;
	    }
	}
    }
}

/**
 * dataset_add_observations:
 * @newobs: number of observations to add.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Extends all series in the dataset by the specified number of
 * extra observations, and initializes all the added values to
 * the missing value code, #NADBL.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_observations (int newobs, double ***pZ, DATAINFO *pdinfo)
{
    double *x;
    int oldn = pdinfo->n;
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
	for (t=oldn; t<bign; t++) {
	    sprintf(pdinfo->S[t], "%d", t + 1);
	}
    }
    
    if (pdinfo->t2 == pdinfo->n - 1) {
	pdinfo->t2 = bign - 1;
    }

    pdinfo->n = bign;

    maybe_extend_trends(*pZ, pdinfo, oldn);
    maybe_extend_dummies(*pZ, pdinfo, oldn);

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
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_drop_observations (int n, double ***pZ, DATAINFO *pdinfo)
{
    double *x;
    int i, newn;

    if (n <= 0) return 0;

    newn = pdinfo->n - n;

    for (i=0; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
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
    
    if (pdinfo->t2 > newn - 1) {
	pdinfo->t2 = newn - 1;
    }

    pdinfo->n = newn;

    /* does daily data need special handling? */
    ntodate(pdinfo->endobs, newn - 1, pdinfo);

    return 0;
}

static int 
dataset_expand_varinfo (int newvars, DATAINFO *pdinfo)
{
    char **varname = NULL;
    char *vector = NULL;
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

    vector = realloc(pdinfo->vector, bigv);
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

static int real_dataset_add_series (int newvars, double *x,
				    double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    int i, n = pdinfo->n, v = pdinfo->v;
    int err = 0;

    newZ = realloc(*pZ, (v + newvars) * sizeof *newZ);  

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
    return real_dataset_add_series(newvars, NULL, pZ, pdinfo);
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
    return real_dataset_add_series(1, x, pZ, pdinfo);
}

/**
 * dataset_add_scalar:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Allocates space for a new scalar member of the dataset.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo)
{
    double **newZ;
    int n = pdinfo->n, v = pdinfo->v; 
    int err = 0;

    newZ = realloc(*pZ, (v + 1) * sizeof *newZ);  

    if (newZ == NULL) {
	err = E_ALLOC;
    } else {
	*pZ = newZ;
    }

    if (!err) {
	newZ[v] = malloc(n * sizeof **newZ);
	if (newZ[v] == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = dataset_expand_varinfo(1, pdinfo);
    }

    if (!err) {
	pdinfo->vector[v] = 0;
    }

    return err;
}

/**
 * dataset_scalar_to_vector:
 * @v: index number of variable to process.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Expands an existing scalar member of a dataset to a
 * full-length vector.  All values are initialized to
 * the missing value code, #NABDL.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_scalar_to_vector (int v, double ***pZ, DATAINFO *pdinfo)
{
    double *tmp;
    int t, err = 0;

    tmp = realloc((*pZ)[v], pdinfo->n * sizeof *tmp);

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	/* initialize all vals to missing */
	for (t=0; t<pdinfo->n; t++) {
	    tmp[t] = NADBL;
	}
	(*pZ)[v] = tmp;
	pdinfo->vector[v] = 1;
    }

    return err;
}

/**
 * dataset_add_scalar_as:
 * @numstr: string representation of numeric value.
 * @newname: name to give the new variable.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Adds to the dataset a new scalar with name @newname and
 * value given by @numstr.  The new variable is added at one
 * level "deeper" (in terms of function execution) than the
 * current level.  This is for use with user-defined functions,
 * where a numeric string is given as a function argument.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_add_scalar_as (const char *numstr, const char *newname,
			   double ***pZ, DATAINFO *pdinfo)
{
    int err = 0;

    if (pdinfo->varinfo == NULL) {
	strcpy(gretl_errmsg, _("Please open a data file first"));
	return 1;
    }

    err = dataset_add_scalar(pZ, pdinfo);
    if (!err) {
	int vnew = pdinfo->v - 1;

	(*pZ)[vnew][0] = atof(numstr);
	strcpy(pdinfo->varname[vnew], newname);
	STACK_LEVEL(pdinfo, vnew) += 1;
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

    if (pdinfo->vector[v]) {
	err = real_dataset_add_series(1, NULL, pZ, pdinfo);
    } else {
	err = dataset_add_scalar(pZ, pdinfo);
    }

    if (!err) {
	int vnew = pdinfo->v - 1;

	if (pdinfo->vector[v]) {
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

int dataset_drop_listed_variables (const int *list, double ***pZ, 
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
		    if (renumber != NULL && !is_hidden_variable(i + gap, pdinfo)) {
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

/**
 * dataset_destroy_hidden_variables:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Deletes from the dataset all "hidden" variables that have
 * been added automatically (for example, auto-generated variables
 * used for the x-axis in graph plotting).  Does not delete the
 * automatically generated constant (ID number 0).
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int dataset_destroy_hidden_variables (double ***pZ, DATAINFO *pdinfo)
{
    int i, nhid = 0;
    int err = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (is_hidden_variable(i, pdinfo)) {
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
		if (is_hidden_variable(i, pdinfo)) {
		    hidlist[j++] = i;
		}
	    }	    
	    err = dataset_drop_listed_variables(hidlist, pZ, pdinfo, NULL);
	    free(hidlist);
	}
    }

    return err;
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
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @newvar: name for new variable, produced by stacking
 * @s: instructions for stacking existing variables.
 *
 * Really for internal use.  Don't worry about it.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int dataset_stack_variables (double ***pZ, DATAINFO *pdinfo, 
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
	err = dataset_add_observations(bign - oldn, pZ, pdinfo);
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

/**
 * is_hidden_variable:
 * @i: ID number of variable.
 * @pdinfo: dataset information.
 *
 * Used in various contexts to screen a list of variables being 
 * presented to the user.
 *
 * Returns: 1 if variable @i is a "hidden", automatically
 * generated variable, otherwise 0.  
 */

int is_hidden_variable (int i, const DATAINFO *pdinfo)
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
