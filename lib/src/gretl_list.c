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

/* gretl_list_new:
 *
 */

int *gretl_list_new (int nterms)
{
    int *list;
    int i;

    list = malloc((nterms + 1) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = nterms;

    for (i=1; i<=nterms; i++) list[i] = 0;

    return list;
}

/* in_gretl_list:
 *
 */

int in_gretl_list (const int *list, int k)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == k) return i;
    }

    return 0;
}

/* checks a list for a constant term (ID 0), and if present, 
   moves it to the first dep var position (pos 2)
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

/* gretl_list_delete_at_pos:
 * @list: original list.
 * @pos: position at which to delete list element.
 *
 * deletes the element at position @pos from @list and moves any
 * remaining elements forward.
 *
 * Returns: 0 on success, 1 on error.
 * 
 */

int gretl_list_delete_at_pos (int *list, int pos)
{
    int i;

    if (pos < 1) return 1;

    for (i=pos; i<list[0]; i++) {
	list[i] = list[i + 1];
    }

    list[list[0]] = 0;
    list[0] -= 1;

    return 0;
}

/* gretl_list_add:
 * @orig: original list.
 * @add: list of variables to add.
 * @err: pointer to receive error code.
 *
 * creates a list containing the union of elements of @orig 
 * and the elements of @add. 
 *
 * Returns: new list on success, %NULL on error.
 * 
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
 * @orig: original list.
 * @err: pointer to receive error code.
 *
 * creates a list containing all but the last elements of @orig.
 *
 * Returns: new list on success, %NULL on error.
 * 
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
 * @orig: original list.
 * @omit: list of variables to drop.
 * @err: pointer to receive error code.
 *
 * creates a list containing the elements of @orig that are not
 * present in @omit. 
 *
 * Returns: new list on success, %NULL on error.
 * 
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
 * fills out list @targ with the elements of biglist that are not
 * present in sublist.
 * 
 */

void gretl_list_diff (int *targ, const int *biglist, const int *sublist)
{
    int i, j, k = 0;
    int match;

    targ[0] = biglist[0] - sublist[0];

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

/* Check if any var in list has been replaced via genr since a
   previous model (ref_id) was estimated.  Expects the "label"
   in datainfo to be of the form "Replaced after model <count>".
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
