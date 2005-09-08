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
#include "gretl_func.h"
#include "libset.h"

#define LDEBUG 0

typedef struct saved_list_ saved_list;

struct saved_list_ {
    int *list;
    char *name;
    int level;
};

static saved_list **list_stack;
static int n_lists;

static saved_list *saved_list_new (const int *list, const char *name)
{
    saved_list *sl = malloc(sizeof *sl);

    if (sl != NULL) {
	if (gretl_executing_function()) {
	    sl->level = gretl_function_stack_depth();
	} else {
	    sl->level = 0;
	}
	if (list != NULL && list[0] > 0) {
	    sl->list = gretl_list_copy(list);
	} else {
	    sl->list = gretl_null_list();
	}
	if (sl->list == NULL) {
	    free(sl);
	    sl = NULL;
	} else {
	    sl->name = gretl_strndup(name, 31);
	    if (sl->name == NULL) {
		free(sl->list);
		free(sl);
		sl = NULL;
	    }
	}
    }

    return sl;
}

static void free_saved_list (saved_list *sl)
{
    free(sl->list);
    free(sl->name);
    free(sl);
}

static saved_list *get_saved_list_by_name (const char *name)
{
    int fsd = 0;
    int i;

    if (gretl_executing_function()) {
	fsd = gretl_function_stack_depth();
    }

    for (i=0; i<n_lists; i++) {
	if (!strcmp(name, list_stack[i]->name) && 
	    fsd == list_stack[i]->level) {
	    return list_stack[i];
	}
    }

    return NULL;
}

/**
 * get_list_by_name:
 * @name: the name of the list to be found.
 *
 * Looks up @name in the stack of saved lists (if any) and
 * retrieves the associated list.
 *
 * Returns: the list, or %NULL if the lookup fails. 
 */

int *get_list_by_name (const char *name)
{
    int *ret = NULL;
    saved_list *sl = get_saved_list_by_name(name);

    if (sl != NULL) {
	ret = sl->list;
    }

    return ret;
}

static int real_remember_list (const int *list, const char *name, 
			       int force_new, PRN *prn)
{
    saved_list *orig = NULL;
    int err = 0;

#if LDEBUG
    fprintf(stderr, "remember_list (in): name='%s', force_new=%d,"
	    " n_lists=%d\n", name, force_new, n_lists);
#endif

    if (!force_new) {
	orig = get_saved_list_by_name(name);
    }

    if (orig != NULL) {
	/* replace existing list of same name */
	free(orig->list);
	orig->list = gretl_list_copy(list);
	if (orig->list == NULL) {
	    pprintf(prn, "Out of memory replacing list '%s'\n", name);
	    err = E_ALLOC;
	} else if (gretl_messages_on()) {
	    pprintf(prn, "Replaced list '%s'\n", name);
	}
    } else {
	saved_list **lstack;

	lstack = realloc(list_stack, (n_lists + 1) * sizeof *lstack);
	if (lstack == NULL) {
	    return E_ALLOC;
	}
	list_stack = lstack;

	list_stack[n_lists] = saved_list_new(list, name);
	if (list_stack[n_lists] == NULL) {
	    pprintf(prn, "Out of memory adding list '%s'\n", name);
	    err = E_ALLOC;
	} else {
	    if (gretl_messages_on()) {
		pprintf(prn, "Added list '%s'\n", name);
	    }
	    n_lists++;
	}
    }

#if LDEBUG
    fprintf(stderr, "remember_list (out): n_lists=%d, ", n_lists);
    fprintf(stderr, "list_stack[%d]=%p\n", n_lists - 1, 
	    (void *) list_stack[n_lists - 1]);
#endif

    return err;
}

/**
 * remember_list:
 * @list: array of integers, the first element being a count
 * of the following elements.
 * @name: name to be given to the list.
 * @prn: printing struct.
 *
 * Adds @list to the stack of saved lists and associates it
 * with @name, unless there is already a list with the given
 * name in which case the original list is replaced.  A status
 * message is printed to @prn.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int remember_list (const int *list, const char *name, PRN *prn)
{
    return real_remember_list(list, name, 0, prn);
}

/**
 * copy_named_list:
 * @orig: the name of the original list.
 * @new: the name to be given to the copy.
 *
 * If a saved list is found by the name @orig, a copy of
 * this list is added to the stack of saved lists under the
 * name @new.  This is intended for use when a list is given
 * as the argument to a user-defined function: it is copied
 * under the name assigned by the function's parameter list.
 *
 * Returns: 0 on success, non-zero on error.
 */

int copy_named_list_as (const char *orig, const char *new)
{
    saved_list *sl;
    int err = 0;

    sl = get_saved_list_by_name(orig);
    if (sl == NULL) {
	err = 1;
    } else {
	err = real_remember_list(sl->list, new, 1, NULL);
	if (!err) {
	    /* for use in functions */
	    sl = list_stack[n_lists - 1];
	    sl->level += 1;
	}
    }

    return err;
}

/**
 * destroy_saved_lists_at_level:
 * @level: stack level of function execution.
 *
 * Destroys and removes from the stack of saved lists all
 * lists that were created at the given @level.  This is 
 * part of the cleanup that is performed when a user-defined
 * function terminates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int destroy_saved_lists_at_level (int level)
{
    saved_list **lstack;
    int i, j, nl = 0;
    int err = 0;

    for (i=0; i<n_lists; i++) {
	if (list_stack[i] == NULL) {
	    break;
	}
	if (list_stack[i]->level == level) {
	    free_saved_list(list_stack[i]);
	    for (j=i; j<n_lists - 1; j++) {
		list_stack[j] = list_stack[j+1];
	    }
	    list_stack[n_lists - 1] = NULL;
	} else {
	    nl++;
	}
    }

    if (nl < n_lists) {
	n_lists = nl;
	lstack = realloc(list_stack, nl * sizeof *list_stack);
	if (lstack == NULL) {
	    err = E_ALLOC;
	} else {
	    list_stack = lstack;
	}
    }

    return err;
}

/**
 * gretl_lists_cleanup:
 *
 * Frees all resources associated with the internal
 * apparatus for saving and retrieving named lists.
 */

void gretl_lists_cleanup (void)
{
    int i;

    for (i=0; i<n_lists; i++) {
	free_saved_list(list_stack[i]);
    }

    free(list_stack);
    list_stack = NULL;
    n_lists = 0;
}

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
 * gretl_list_to_string:
 * @list: array of integers.
 * 
 * Prints the given @list of integers into a newly
 * allocated string, separated by single spaces and with
 * one leading space.  This function is designed to handle 
 * positive integers in a range that is sensible for ID 
 * numbers of variables, typically with three digits or less, 
 * and will fail if the list contains any numbers greater 
 * than 999.
 *
 * Returns: The string representation of the list on success,
 * or %NULL on failure.
 */

char *gretl_list_to_string (const int *list)
{
    char *buf;
    char numstr[8];
    int len, i, err = 0;

    len = 4 * (list[0] + 1);
    if (len > MAXLINE - 32) {
	return NULL;
    }

    buf = malloc(len);
    if (buf == NULL) {
	return NULL;
    }

    *buf = '\0';
    for (i=1; i<=list[0]; i++) {
	if (abs(list[i] > 999)) {
	    err = 1;
	    break;
	}
	sprintf(numstr, " %d", list[i]);
	strcat(buf, numstr);
    }

    if (err) {
	free(buf);
	buf = NULL;
    }

    return buf;
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
 * gretl_list_purge_const:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Checks @list from position 1 onward for the presence of a constant
 * (the variable with ID number 0).  If this is found, is is deleted
 * from list (that is, any following elements are moved forward by one
 * and list[0] is decremented by 1.
 *
 * Returns: 1 if the constant was found and deleted, else 0.
 */

int gretl_list_purge_const (int *list)
{
    int i, gotc = 0;
    int l0 = list[0];

    /* handle the case where the constant comes last; if it's
       the only element behind the list separator, remove both
       the constant and the separator */
    if (list[l0] == 0) {
	gotc = 1;
	list[0] -= 1;
	if (list[l0 - 1] == LISTSEP) {
	    list[l0 - 1] = 0;
	    list[0] -= 1;
	}
    } else {
	for (i=1; i<l0; i++) {
	    if (list[i] == 0) {
		for ( ; i<l0; i++) {
		    list[i] = list[i+1];
		}
		list[l0] = 0;
		list[0] -= 1;
		gotc = 1;
		break;
	    }
	}
    }

    return gotc;
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
 * gretl_list_add_list:
 * @targ: location of list to which @src should be added.
 * @src: list to be added to @targ.
 *
 * Adds @src onto the end of @targ.  The length of @targ becomes the
 * sum of the lengths of the two original lists.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int gretl_list_add_list (int **targ, const int *src)
{
    int *big;
    int n1 = (*targ)[0];
    int n2 = src[0];
    int i, err = 0;

    big = realloc(*targ, (n1 + n2 + 1) * sizeof *big);
    if (big == NULL) {
	err = E_ALLOC;
    } else {
	big[0] = n1 + n2;
	for (i=1; i<=src[0]; i++) {
	    big[n1 + i] = src[i];
	}
	*targ = big;
    }

    return err;
}

/**
 * gretl_list_insert_list:
 * @targ: location of list into which @src should be inserted.
 * @src: list to be inserted.
 * @pos: zero-based position at which @src should be inserted.
 *
 * Inserts @src into @targ at @pos.  The length of @targ becomes the
 * sum of the lengths of the two original lists.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_insert_list (int **targ, const int *src, int pos)
{
    int *big;
    int n1 = (*targ)[0];
    int n2 = src[0];
    int bign = n1 + n2;
    int i, err = 0;

    if (pos > n1 + 1) {
	return 1;
    }

    big = realloc(*targ, (bign + 1) * sizeof *big);
    if (big == NULL) {
	err = E_ALLOC;
    } else {
	big[0] = bign;
	for (i=bign; i>=pos+n2; i--) {
	    big[i] = big[i-n2];
	}
	for (i=1; i<=src[0]; i++) {
	    big[pos+i-1] = src[i];
	}
	*targ = big;
    }

    return err;
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
 * gretl_list_has_separator:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 *
 * Returns: 1 if @list contains the separator for compound
 * lists, #LISTSEP, else 0.  The search begins at position 2.
 */

int gretl_list_has_separator (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
	    return 1;
	}
    }

    return 0;
}

/**
 * gretl_list_split_on_separator:
 * @list: source list.
 * @plist1: pointer to accept first sub-list.
 * @plist2: pointer to accept second sub-list.
 *
 * If @list contains the list separator, #LISTSEP, creates two
 * sub-lists, one containing the elements of @list preceding
 * the separator and one containing the elements following
 * the separator.  The sub-lists are newly allocated, and assigned
 * as the content of @plist1 and @plist2 respectively.
 *
 * Returns: 0 on success, %E_ALLOC is memory allocation fails,
 * or %E_DATA if @list does not contain a separator.
 */

int gretl_list_split_on_separator (const int *list, int **plist1, int **plist2)
{
    int *list1 = NULL, *list2 = NULL;
    int i, n = -1;

    for (i=1; i<list[0] && n<0; i++) {
	if (list[i] == LISTSEP) {
	    n = i;
	}
    }

    if (n < 0) {
	return 1;
    }

    list1 = gretl_list_new(n - 1);
    if (list1 == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<n; i++) {
	list1[i] = list[i];
    }

    list2 = gretl_list_new(list[0] - n);
    if (list2 == NULL) {
	free(list1);
	return E_ALLOC;
    }

    for (i=1; i<=list2[0]; i++) {
	list2[i] = list[i + n];
    }

    *plist1 = list1;
    *plist2 = list2;
    
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
 * Returns: the ID number of the first duplicated variable found,
 * or -1 in case of no duplication.
 */

int gretl_list_duplicates (const int *list, GretlCmdIndex ci)
{
    int i, j, start = 2;
    int ret = -1;

    if (ci == ARCH) {
	start = 3;
    } else if (ci == ARMA) {
	for (i=list[0]-1; i>2; i--) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
    } else if (ci == LAGS && list[0] > 1 && list[2] == LISTSEP) {
	start = 3;
    } else if (ci == TSLS || ci == AR || ci == SCATTERS || 
	       ci == MPOLS || ci == GARCH) {
	for (i=2; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
    }
    
    for (i=start; i<list[0] && ret < 0; i++) {
	for (j=i+1; j<=list[0] && ret < 0; j++) {
	    if (list[i] == list[j]) {
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
