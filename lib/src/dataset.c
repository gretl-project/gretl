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

static int real_dataset_add_series (int newvars, double *x,
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
