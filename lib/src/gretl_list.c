/*
 *  Copyright (c) Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_list.h"

/**
 * gretl_list_new:
 * @nterms: the maximum number of elements to be stored in the list.
 * 
 * Creates a newly allocated list with space for @nterms elements,
 * besides the leading element, which in a gretl list always
 * holds a count of the number of elements that follow.  This
 * leading element is initialized appropriately.  For example, if
 * @nterms = 4, space for 5 integers is allocated and the first
 * element of the array is set to 4.
 *
 * Returns: the newly allocated list, or %NULL on failure.
 */

int *gretl_list_new (int nterms)
{
    int *list = NULL;
    int i;

    list = malloc((nterms + 1) * sizeof *list);

    if (list != NULL) {
	list[0] = nterms;
	for (i=1; i<=nterms; i++) {
	    list[i] = 0;
	}
    }

    return list;
}

/**
 * gretl_null_list:
 * 
 * Creates a newly allocated "list" with only one member, 
 * which is set to zero.
 *
 * Returns: the newly allocated list, or %NULL on failure.
 */

int *gretl_null_list (void)
{
    int *list = malloc(sizeof *list);

    if (list != NULL) {
	list[0] = 0;
    }

    return list;
}

/**
 * gretl_list_copy:
 * @src: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Returns: an allocated copy @src.
 */

int *gretl_list_copy (const int *src)
{
    int *targ = NULL;
    int i;

    if (src != NULL && src[0] != 0) {
	targ = malloc((src[0] + 1) * sizeof *targ);
	if (targ != NULL) {
	    for (i=0; i<=src[0]; i++) {
		targ[i] = src[i];
	    }
	}
    }

    return targ;
}

/**
 * gretl_list_from_string:
 * @liststr: string representation of list of integers.
 *
 * Reads a string containing a list of integers and constructs
 * an array of these integers.  The first element is the number
 * of integers that follow.
 *
 * Returns: the allocated array.
 */

int *gretl_list_from_string (const char *liststr)
{
    const char *s = liststr;
    char numstr[8];
    int *list;
    int n = 0;

    while (*s) {
	while (*s == ' ') s++;
	if (sscanf(s, "%7s", numstr)) {
	    n++;
	    s += strlen(numstr);
	}
    }

    if (n == 0) {
	return NULL;
    }

    list = malloc((n + 1) * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    list[0] = n;

    s = liststr;
    n = 1;
    while (*s) {
	while (*s == ' ') s++;
	if (sscanf(s, "%7s", numstr)) {
	    list[n++] = atoi(numstr);
	    s += strlen(numstr);

	}
    }    

    return list;
}

/**
 * in_gretl_list:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @k: integer to test.
 *
 * Checks whether @k is present among the members of @list,
 * in position 1 or higher.
 *
 * Returns: the position of @k in @list, or 0 if @k is not
 * present. 
 */

int in_gretl_list (const int *list, int k)
{
    int i, ret = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == k) {
	    ret = i;
	    break;
	}
    }

    return ret;
}

/**
 * rearrange_list:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Checks @list for the presence of a constant term (ID 0), and 
 * if present, moves it to position 2 in @list.  This is
 * designed for producing a canonical version of a regression
 * list, with the constant (if any) preceding any other
 * regressors (position 1 holds the dependent variable in such
 * a list).
 */

void rearrange_list (int *list)
{
    int i, v;

    for (v=list[0]; v>2; v--) {
        if (list[v] == 0)  {
	    for (i=v; i>2; i--) {
		list[i] = list[i-1];
	    }
	    list[2] = 0;
	    return;
        }
    }
}

/**
 * gretl_list_delete_at_pos:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @pos: position at which to delete list element.
 *
 * Deletes the element at position @pos from @list and moves any
 * remaining elements forward.  Decrements the value of the first,
 * counter, element of @list.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_list_delete_at_pos (int *list, int pos)
{
    int i, err = 0;

    if (pos < 1 || pos > list[0]) {
	err = 1;
    } else {
	for (i=pos; i<list[0]; i++) {
	    list[i] = list[i + 1];
	}

	list[list[0]] = 0;
	list[0] -= 1;
    }

    return 0;
}

/**
 * gretl_list_add:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @add: list of variables to be added.
 * @err: location to receive error code.
 *
 * Creates a list containing the union of elements of @orig 
 * and the elements of @add.  If one or more elements of
 * @add were already present in @orig, the error code is
 * %E_ADDDUP.
 *
 * Returns: new list on success, %NULL on error.
 */

int *gretl_list_add (const int *orig, const int *add, int *err)
{
    int i, j, k;
    int *big;
    const int norig = orig[0];
    const int nadd = add[0];

    *err = 0;

    big = malloc((norig + nadd + 1) * sizeof *big);
    if (big == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<=norig; i++) {
	big[i] = orig[i];
    }

    k = orig[0];

    for (i=1; i<=nadd; i++) {
	int match = 0;

	for (j=1; j<=norig; j++) {
	    if (add[i] == orig[j]) {
		/* a "new" var was already present */
		free(big);
		*err = E_ADDDUP;
		return NULL;
	    }
	}
	if (!match) {
	    big[0] += 1;
	    big[++k] = add[i];
	}
    }

    if (big[0] == norig) {
	free(big);
	*err = E_NOADD;
	return NULL;
    }

    return big;
}

/**
 * gretl_list_omit_last:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @err: location to receive error code.
 *
 * Creates a list containing all but the last elements of @orig.
 *
 * Returns: new list on success, %NULL on error.
 */

int *gretl_list_omit_last (const int *orig, int *err)
{
    int *list = NULL;
    int i;

    *err = 0;

    if (orig[0] < 2) {
	*err = E_NOVARS;
    }

    /* can't handle compound lists */
    if (*err == 0) {
	for (i=1; i<=orig[0]; i++) {
	    if (orig[i] == LISTSEP) {
		*err = 1;
		break;
	    }
	}
    }

    if (*err == 0) {
	list = malloc(orig[0] * sizeof *list);
	if (list == NULL) {
	   *err = E_ALLOC;
	} 
    }
    
    if (list != NULL) {
	list[0] = orig[0] - 1;
	for (i=1; i<orig[0]; i++) {
	    list[i] = orig[i];
	}
    }

    return list;
}

static int list_count (const int *list)
{
    int i, k = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    break;
	} else {
	    k++;
	}
    }

    return k;
}

/**
 * gretl_list_omit:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @omit: list of variables to drop.
 * @err: pointer to receive error code.
 *
 * Creates a list containing the elements of @orig that are not
 * present in @omit. 
 *
 * Returns: new list on success, %NULL on error.
 */

int *gretl_list_omit (const int *orig, const int *omit, int *err)
{
    int i, j, k;
    int *smal;
    const int nomit = omit[0];
    const int norig = list_count(orig);

    *err = 0;

    /* check for spurious "omissions" */
    for (i=1; i<=nomit; i++) {
	int pos = in_gretl_list(orig, omit[i]);

	if (pos <= 1) {
	    sprintf(gretl_errmsg, _("Variable %d was not in the original list"),
		    omit[i]);
	    *err = 1;
	    return NULL;
	}
    }

    /* check for attempt to omit all vars */
    if (nomit == norig - 1) {
	*err = E_NOVARS;
	return NULL;
    }

    smal = malloc((norig - nomit + 1) * sizeof *smal);
    if (smal == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    smal[0] = norig - nomit;
    smal[1] = orig[1];

    k = 1;
    for (i=2; i<=norig; i++) {
        int match = 0;

        for (j=1; j<=nomit; j++) {
            if (orig[i] == omit[j]) {
                match = 1; /* matching var: omit it */
		break;
            }
        }
        if (!match) { /* var is not in omit list: keep it */
            smal[++k] = orig[i];
        }
    }

    return smal;
}

/**
 * gretl_list_diff:
 * @targ: target list (must be pre-allocated).
 * @biglist: inclusive list.
 * @sublist: subset of biglist.
 *
 * Fills out @targ with the elements of @biglist, from position 2 
 * onwards, that are not present in @sublist.  It is assumed that 
 * the variable ID number in position 1 (dependent variable) is the 
 * same in both lists.  It is an error if, from position 2 on, 
 * @sublist is not a proper subset of @biglist.  See also 
 * #gretl_list_diff_new.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_list_diff (int *targ, const int *biglist, const int *sublist)
{
    int i, j, k = 0;
    int match, err = 0;

    targ[0] = biglist[0] - sublist[0];

    if (targ[0] <= 0) {
	err = 1;
    } else {
	for (i=2; i<=biglist[0]; i++) {
	    match = 0;
	    for (j=2; j<=sublist[0]; j++) {
		if (sublist[j] == biglist[i]) {
		    match = 1;
		    break;
		}
	    }
	    if (!match) {
		targ[++k] = biglist[i];
	    }
	}
    }

    return err;
}

/**
 * gretl_list_diff_new:
 * @biglist: inclusive list.
 * @sublist: subset of biglist.
 *
 * Returns: a newly allocated list including the elements of @biglist,
 * from position 2 onwards, that are not present in @sublist, or %NULL on 
 * failure.  It is assumed that the variable ID number in position 1
 * (dependent variable) is the same in both lists.  It is an error
 * if, from position 2 on, @sublist is not a proper subset of @biglist.
 */

int *gretl_list_diff_new (const int *biglist, const int *sublist)
{
    int *targ = NULL;
    int i, j, n, k = 0;
    int match;

    n = biglist[0] - sublist[0];

    if (n <= 0) {
	return NULL;
    }

    targ = gretl_list_new(n);

    for (i=2; i<=biglist[0]; i++) {
	match = 0;
	for (j=2; j<=sublist[0]; j++) {
	    if (sublist[j] == biglist[i]) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    targ[++k] = biglist[i];
	}
    }

    return targ;
}

/**
 * list_members_replaced:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 * @pdinfo: dataset information.
 * @ref_id: ID number of reference #MODEL.
 *
 * Checks whether any variable in @list has been redefined via
 * gretl's %genr command since a previous model (identified by
 * @ref_id) was estimated.
 *
 * Returns: 1 if any variables have been replaced, 0 otherwise.
 */

int list_members_replaced (const int *list, const DATAINFO *pdinfo,
			   int ref_id)
{
    const char *label;
    char rword[16];
    int j, mc, repl, err = 0;

    if (ref_id == 0) {
	mc = get_model_count();
    } else {
	mc = ref_id;
    }

    for (j=1; j<=list[0]; j++) {
	if (list[j] == LISTSEP) {
	    continue;
	}
	label = VARLABEL(pdinfo, list[j]);
	*rword = '\0';
	sscanf(label, "%15s", rword);
	if (!strcmp(rword, _("Replaced"))) {
	    repl = 0;
	    sscanf(label, "%*s %*s %*s %d", &repl);
	    if (repl >= mc) {
		err = E_VARCHANGE;
		break;
	    }
	}
    }

    return err;
}

/**
 * gretl_list_has_const:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 *
 * Returns: 1 if the constant (variable ID 0) is found in @list,
 * in position 2 or higher.  This corresponds to determining
 * whether of not a set of regressors contains an intercept.
 */

int gretl_list_has_const (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
        if (list[i] == 0) {
	    return 1;
	}
    }

    return 0;
}

/**
 * gretl_list_duplicates:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 * @ci: index of gretl command (for context).
 *
 * Checks whether or not a gretl list contains duplicated elements.
 * Exactly what counts as duplication depends on the context of the
 * command in which @list will be used, which is given by @ci.
 *
 * Returns: 1 in case of erroneous duplication, 0 otherwise.
 */

int gretl_list_duplicates (const int *list, GretlCmdIndex ci)
{
    int i, j, start = 2;
    int ret = 0;

    if (ci == ARCH) {
	start = 3;
    } else if (ci == LAGS && list[0] > 1 && list[2] == LISTSEP) {
	start = 3;
    } else if (ci == TSLS || ci == AR || ci == ARMA || 
	       ci == SCATTERS || ci == MPOLS || ci == GARCH) {
	for (i=2; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
    }
    
    for (i=start; i<list[0] && !ret; i++) {
	for (j=start+1; j<=list[0] && !ret; j++) {
	    if (i != j && list[i] == list[j]) {
		ret = list[i];
	    }
	}
    }

    return ret;
}

/**
 * full_var_list:
 * @pdinfo: dataset information.
 * @nvars: location for return of number of elements in full list.
 *
 * Creates a newly allocated list including all variables in the
 * dataset that are not scalars and are not hidden variables.
 * The return value is %NULL in case either (a) allocation of
 * memory failed, or (b) the resulting list would be empty.
 * The caller can distinguish between these possibilities by
 * examining the value returned in @nvars, which will be zero if
 * and only if the resulting list would be empty.  If this is
 * not of interest to the caller, @nvars may be given as %NULL.
 *
 * Returns: the allocated list, or %NULL.
 */

int *full_var_list (const DATAINFO *pdinfo, int *nvars)
{
    int i, j, nv = 0;
    int *list = NULL;

    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i] && !is_hidden_variable(i, pdinfo)) {
	    nv++;
	}
    }

    if (nvars != NULL) {
	*nvars = nv;
    }
    
    if (nv > 0) {
	list = gretl_list_new(nv);
    }

    if (list != NULL) {
	j = 1;
	for (i=1; i<pdinfo->v; i++) {
	    if (pdinfo->vector[i] && !is_hidden_variable(i, pdinfo)) {
		list[j++] = i;
	    }
	}
    }	    

    return list;
}
