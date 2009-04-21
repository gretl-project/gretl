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
#include "gretl_list.h"
#include "gretl_func.h"
#include "libset.h"
#include "gretl_xml.h"

#include <errno.h>

#define LDEBUG 0

typedef struct saved_list_ saved_list;

struct saved_list_ {
    char name[VNAMELEN];
    int *list;
    int level;
};

enum {
    LIST_OMIT_TIGHT,
    LIST_OMIT_SLOPPY
};

static saved_list **list_stack;
static int n_lists;

static int *real_gretl_list_omit (const int *orig, const int *omit, 
				  int minpos, int mode, int *err);

static saved_list *saved_list_new (const int *list, const char *name)
{
    saved_list *sl = malloc(sizeof *sl);

    if (sl != NULL) {
	sl->level = gretl_function_depth();
	if (list != NULL && list[0] > 0) {
	    sl->list = gretl_list_copy(list);
	} else {
	    sl->list = gretl_null_list();
	}
	if (sl->list == NULL) {
	    free(sl);
	    sl = NULL;
	} else {
	    *sl->name = 0;
	    strncat(sl->name, name, VNAMELEN - 1);
	}
    }

    return sl;
}

static void free_saved_list (saved_list *sl)
{
    if (sl != NULL) {
	free(sl->list);
	free(sl);
    }
}

static saved_list *get_saved_list_by_name (const char *name)
{
    saved_list *sl = NULL;
    int fsd = gretl_function_depth();
    int i;

    for (i=0; i<n_lists; i++) {
	if (!strcmp(name, list_stack[i]->name) && 
	    fsd == list_stack[i]->level) {
	    sl = list_stack[i];
	    break;
	}
    }

    return sl;
}

/**
 * n_saved_lists:
 *
 * Returns: the number of saved (named) lists currently defined.
 */

int n_saved_lists (void)
{
    return n_lists;
}

/**
 * max_varno_in_saved_lists:
 *
 * Returns: the highest ID number of a variable referenced
 * in a saved (named) list.
 */

int max_varno_in_saved_lists (void)
{
    int *list;
    int i, j, vmax = 0;

    for (i=0; i<n_lists; i++) {
	list = list_stack[i]->list;
	for (j=1; j<=list[0]; j++) {
	    if (list[j] > vmax) {
		vmax = list[j];
	    }
	}
    }    

    return vmax;
}

/**
 * get_list_name_by_index:
 * @idx: 0-based index into array of saved lists.
 *
 * Returns: the name of the specified saved list, or %NULL if
 * @idx is out of bounds.
 */

const char *get_list_name_by_index (int idx)
{
    if (idx >= 0 && idx < n_lists) {
	return list_stack[idx]->name;
    } else {
	return NULL;
    }
}

/**
 * get_list_by_name:
 * @name: the name of the list to be found.
 *
 * Looks up @name in the stack of saved lists, at the current level
 * of function execution, and retrieves the associated list.
 *
 * Returns: the list, or %NULL if the lookup fails. 
 */

int *get_list_by_name (const char *name)
{
    int *ret = NULL;
    saved_list *sl;

    if (name != NULL) {
	sl = get_saved_list_by_name(name);
	if (sl != NULL) {
	    ret = sl->list;
	}
    }

    return ret;
}

/**
 * append_to_list_by_name:
 * @targ: the name of the target list.
 * @add: list to add.
 *
 * If @targ is the name of a saved list, append the list
 * @add to it.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int append_to_list_by_name (const char *targ, const int *add)
{
    saved_list *sl = get_saved_list_by_name(targ);
    int err = 0;

    if (sl == NULL) {
	err = E_UNKVAR;
    } else {
	int *tmp = gretl_list_copy(sl->list);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_list_add_list(&tmp, add);
	    if (!err) {
		free(sl->list);
		sl->list = tmp;
	    }
	}
    } 

    return err;
}

/**
 * subtract_from_list_by_name:
 * @targ: the name of the target list.
 * @sub: sub-list to remove.
 *
 * If @targ is the name of a saved list, remove from @targ
 * any elements of @sub that it contains.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int subtract_from_list_by_name (const char *targ, const int *sub)
{
    saved_list *sl = get_saved_list_by_name(targ);
    int err = 0;

    if (sl == NULL) {
	err = E_UNKVAR;
    } else {
	int *tmp = real_gretl_list_omit(sl->list, sub, 1, 
					LIST_OMIT_SLOPPY,
					&err);

	if (!err) {
	    free(sl->list);
	    sl->list = tmp;
	}
    } 

    return err;
}

/**
 * replace_list_by_name:
 * @targ: the name of the target list.
 * @new: replacement list
 *
 * If @targ is the name of a saved list, replace the
 * list of that name with @new.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int replace_list_by_name (const char *targ, const int *new)
{
    saved_list *sl = get_saved_list_by_name(targ);
    int err = 0;

    if (sl == NULL) {
	err = E_UNKVAR;
    } else {
	int *tmp = gretl_list_copy(new);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    free(sl->list);
	    sl->list = tmp;
	}	    
    }

    return err;
}

static int real_remember_list (const int *list, const char *name, 
			       int force_new, PRN *prn)
{
    saved_list *orig = NULL;
    int err = 0;

    if (list == NULL) {
	return E_DATA;
    }

#if LDEBUG
    fprintf(stderr, "remember_list (in): name='%s', force_new=%d,"
	    " n_lists=%d\n", name, force_new, n_lists);
#endif

    /* Note: 'force_new' means that we'll add a new list even if it
       has the same name as an existing one.  This makes sense only if
       we're copying a list in the context of running a user-defined
       function, since in that case the new list will exist at a
       different "stack level" from any prior list.
    */

    if (!force_new) {
	orig = get_saved_list_by_name(name);
    }

    if (orig != NULL) {
	/* replace existing list of same name */
	free(orig->list);
	orig->list = gretl_list_copy(list);
	if (orig->list == NULL) {
	    err = E_ALLOC;
	} else if (gretl_messages_on() && !gretl_looping_quietly()) {
	    pprintf(prn, _("Replaced list '%s'\n"), name);
	}
    } else {
	saved_list **lstack;

	err = check_varname(name);
	if (err) {
	    return err;
	}

	lstack = realloc(list_stack, (n_lists + 1) * sizeof *lstack);
	if (lstack == NULL) {
	    return E_ALLOC;
	}
	list_stack = lstack;

	list_stack[n_lists] = saved_list_new(list, name);
	if (list_stack[n_lists] == NULL) {
	    err = E_ALLOC;
	} else {
	    if (gretl_messages_on()) {
		const char *realname = list_stack[n_lists]->name;

		if (list[0] > 0) {
		    pprintf(prn, _("Added list '%s'\n"), realname);
		}
		if (strlen(realname) < strlen(name)) {
		    pprintf(prn, _("Warning: the name was truncated to %d characters\n"), 
			    VNAMELEN - 1);
		}
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
 * Adds a copy of @list to the stack of saved lists and associates
 * it with @name, unless there is already a list with the given
 * name in which case the original list is replaced.  A status
 * message is printed to @prn.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int remember_list (const int *list, const char *name, PRN *prn)
{
    return real_remember_list(list, name, 0, prn);
}

static int destroy_saved_list (saved_list *sl)
{
    saved_list **lstack;
    int found = 0;
    int i, j, err = 0;

    for (i=0; i<n_lists; i++) {
	if (list_stack[i] == NULL) {
	    break;
	}
	if (list_stack[i] == sl) {
	    free_saved_list(sl);
	    for (j=i; j<n_lists - 1; j++) {
		list_stack[j] = list_stack[j+1];
	    }
	    list_stack[n_lists - 1] = NULL;
	    found = 1;
	    break;
	} 
    }

    if (!found) {
	err = E_DATA;
    } else {
	n_lists--;
	if (n_lists == 0) {
	    free(list_stack);
	    list_stack = NULL;
	} else {
	    lstack = realloc(list_stack, n_lists * sizeof *list_stack);
	    if (lstack == NULL) {
		err = E_ALLOC;
	    } else {
		list_stack = lstack;
	    } 
	}
    }   

    return err;
}

/**
 * rename_saved_list:
 * @orig: the original name of the list.
 * @new: the new name to be given.
 *
 * Renames a saved list from @orig to @new.  If there is
 * already a list called @new, it is destroyed.
 *
 * Returns: 0 on success, non-zero on error.
 */

int rename_saved_list (const char *orig, const char *new)
{
    saved_list *sl0, *sl1;
    int err = 0;

    sl0 = get_saved_list_by_name(orig);
    if (sl0 == NULL) {
	err = 1;
    } else {
	/* is there already a list called @new? */
	sl1 = get_saved_list_by_name(new);
	if (sl1 != NULL) {
	    err = destroy_saved_list(sl1);
	}
	*sl0->name = '\0';
	strncat(sl0->name, new, VNAMELEN - 1);
    } 

    return err;
}

/**
 * copy_named_list_as:
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
	    sl = list_stack[n_lists - 1];
	    sl->level += 1;
	}
    }

    return err;
}

/**
 * named_list_lower_level:
 * @name: the name of the list.
 *
 * If a saved list is found by the name @name, at the
 * current level of function execution, lower its level
 * by 1.  This is intended for use when a list is 
 * returned by a user-defined function: it is shifted to
 * the level of the caller.
 *
 * Returns: 0 on success, non-zero on error.
 */ 

int named_list_lower_level (const char *name)
{
    saved_list *sl;
    int err = 0;

    sl = get_saved_list_by_name(name);
    if (sl == NULL) {
	err = E_DATA;
    } else {
	sl->level -= 1;
    }

    return err;
}

/**
 * create_named_null_list:
 * @name: the name to be given to the list.
 *
 * Creates an empty list under the given @name and adds it
 * to the stack of saved lists at the next level of
 * function execution.  This is intended for use when a null 
 * list is given as an argument to a user-defined function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int create_named_null_list (const char *name)
{
    saved_list *sl;
    int *list;
    int err;

    list = gretl_null_list();
    if (list == NULL) {
	return E_ALLOC;
    }

    err = real_remember_list(list, name, 1, NULL);
    if (!err) {
	sl = list_stack[n_lists - 1];
	sl->level += 1;
    } else {
	free(list);
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
	    i--;
	} else {
	    nl++;
	}
    }

    if (nl < n_lists) {
	n_lists = nl;
	if (nl == 0) {
	    free(list_stack);
	    list_stack = NULL;
	} else {
	    lstack = realloc(list_stack, nl * sizeof *list_stack);
	    if (lstack == NULL) {
		err = E_ALLOC;
	    } else {
		list_stack = lstack;
	    }
	}
    }

    return err;
}

int delete_list_by_name (const char *name)
{
    saved_list *sl;

    sl = get_saved_list_by_name(name);
    if (sl == NULL) {
	return E_UNKVAR;
    } else {
	return destroy_saved_list(sl);
    }
}

static int var_is_deleted (const int *dlist, int dmin, int i)
{
    int v = dmin + i - 1;

    if (dlist != NULL) {
	return in_gretl_list(dlist, v);
    } else {
	return (v >= dmin);
    }
}

/**
 * gretl_lists_revise:
 * @dlist: list of variables to be deleted (or %NULL).
 * @dmin: lowest ID number of deleted var (referenced only
 * if @dlist is %NULL).
 *
 * Goes through any saved lists, adjusting the ID numbers
 * they contain to reflect the deletion from the dataset of
 * certain variables: those referenced in @dlist, if given, 
 * or if @dlist is %NULL, those variables with IDs greater 
 * than or equal to @dmin.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int gretl_lists_revise (const int *dlist, int dmin)
{
    int *list, *maplist;
    int lmax = 0;
    int i, j, k;

#if 0
    fprintf(stderr, "gretl_lists_revise: dlist = %p, dmin = %d\n", 
	    (void *) dlist, dmin);
#endif

    if (dlist != NULL) {
	/* determine lowest deleted ID */
	dmin = dlist[1];
	for (i=2; i<=dlist[0]; i++) {
	    if (dlist[i] > 0 && dlist[i] < dmin) {
		dmin = dlist[i];
	    }
	}
    }

    /* find highest ID ref'd in any saved list */
    for (j=0; j<n_lists; j++) {
	list = list_stack[j]->list;
	for (i=1; i<=list[0]; i++) {
	    if (list[i] > lmax) {
		lmax = list[i];
	    }
	}
    }

    if (lmax < dmin) {
	/* nothing to be done */
	return 0;
    }

    /* make mapping from old to new IDs */

    maplist = gretl_list_new(lmax - dmin + 1);
    if (maplist == NULL) {
	return E_ALLOC;
    }

    j = dmin;

    for (i=1; i<=maplist[0]; i++) {
	if (var_is_deleted(dlist, dmin, i)) {
	    maplist[i] = -1;
	} else {
	    maplist[i] = j++;
	}
    }

    /* use mapping to revise saved lists */
    for (j=0; j<n_lists; j++) {
	list = list_stack[j]->list;
	for (i=list[0]; i>0; i--) {
	    k = list[i] - dmin + 1;
	    if (k >= 1) {
		if (maplist[k] == -1) {
		    gretl_list_delete_at_pos(list, i);
		} else {
		    list[i] = maplist[k];
		}
	    }
	}
    }

    free(maplist);

    return 0;
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
 * element of the array is set to 4.  The other elements of 
 * the list are initialized to 0.
 *
 * Returns: the newly allocated list, or %NULL on failure.
 */

int *gretl_list_new (int nterms)
{
    int *list = NULL;
    int i;
    
    if (nterms < 0) {
	return NULL;
    }

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
 * gretl_consecutive_list_new:
 * @lmin: starting value for consecutive list elements.
 * @lmax: ending value.
 * 
 * Creates a newly allocated list whose elements run from
 * @lmin to @lmax consecutively.
 *
 * Returns: the newly allocated list, or %NULL on failure.
 */

int *gretl_consecutive_list_new (int lmin, int lmax)
{
    int *list = NULL;
    int i, n;

    n = lmax - lmin + 1;
    if (n <= 0) {
	return NULL;
    }

    list = gretl_list_new(n);

    if (list != NULL) {
	for (i=0; i<n; i++) {
	    list[i+1] = lmin + i;
	}
    }

    return list;
}

/**
 * gretl_list_resize:
 * @oldlist: pointer to list to be resized.
 * @nterms: the new maximum number of elements for the list.
 * 
 * Resizes the content of @oldlist to hold @nterms, and adjusts
 * the first element to reflect the new size.  If the new
 * list is longer than the old, the extra elements are initialized 
 * to zero.
 *
 * Returns: the resized list, or %NULL on failure.
 */

int *gretl_list_resize (int **oldlist, int nterms)
{
    int *list = NULL;
    int i, oldn = 0;

    if (nterms < 0 || oldlist == NULL) {
	return NULL;
    }

    if (*oldlist != NULL) {
	oldn = (*oldlist)[0];
	if (nterms == oldn) {
	    /* no-op */
	    return *oldlist;
	}
    }

    list = realloc(*oldlist, (nterms + 1) * sizeof *list);

    if (list != NULL) {
	list[0] = nterms;
	*oldlist = list;
	for (i=oldn+1; i<=list[0]; i++) {
	    list[i] = 0;
	}
    } else {
	free(*oldlist);
	*oldlist = NULL;
    }

    return list;
}

/**
 * gretl_list_append_term:
 * @plist: pointer to list to be augmented.
 * @v: the term to be added.
 * 
 * Resizes (or allocates from scratch) the content of @plist,
 * so that it can hold one extra element, and sets the last
 * element to @v.
 *
 * Returns: the augmented list, or %NULL on failure.
 */

int *gretl_list_append_term (int **plist, int v)
{
    int *list = NULL;

    if (*plist == NULL) {
	list = gretl_list_new(1);
	if (list != NULL) {
	    list[1] = v;
	}
    } else {
	int oldn = (*plist)[0];

	list = realloc(*plist, (oldn + 2) * sizeof *list);
	if (list != NULL) {
	    list[0] += 1;
	    list[list[0]] = v;
	} else {
	    free(*plist);
	}
    }

    *plist = list;

    return list;
}

/**
 * gretl_list_sort:
 * @list: list to be sorted.
 * 
 * Sorts the elements from position 1 to the end of @list
 * in ascending order.
 *
 * Returns: the sorted list.
 */

int *gretl_list_sort (int *list)
{
    int i, sorted = 1;

    for (i=1; i<list[0]; i++) {
	if (list[i] > list[i+1]) {
	    sorted = 0;
	    break;
	}
    }

    if (!sorted) {
	qsort(list + 1, list[0], sizeof *list, gretl_compare_ints);
    }

    return list;
}

/**
 * gretl_list_cmp:
 * @list1: gretl list.
 * @list2: gretl list.
 * 
 * Returns: 0 if @list1 and @list2 have identical content,
 * otherwise 1.
 */

int gretl_list_cmp (const int *list1, const int *list2)
{
    int i;

    for (i=0; i<=list1[0]; i++) {
	if (list2[i] != list1[i]) {
	    return 1;
	}
    }

    return 0;
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
 * Returns: an allocated copy @src (or %NULL if @src is %NULL).
 */

int *gretl_list_copy (const int *src)
{
    int *targ = NULL;
    int i;

    if (src != NULL) {
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
 * gretl_list_copy_from_pos:
 * @src: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Returns: an allocated copy @src from position @pos onward
 * (or %NULL on failure).
 */

int *gretl_list_copy_from_pos (const int *src, int pos)
{
    int *targ = NULL;
    int i, n;

    if (src != NULL && ((n = src[0] - pos + 1) > 0)) {
	targ = gretl_list_new(n);
	if (targ != NULL) {
	    for (i=1; i<=n; i++) {
		targ[i] = src[i+pos-1];
	    }
	}
    }

    return targ;
}

/**
 * gretl_list_from_string:
 * @str: string representation of list of integers.
 * @err: location to receive error code.
 *
 * Reads a string containing a list of integers, separated by
 * spaces and/or commas and possibly wrapped in parentheses,
 * and constructs an array of these integers.  The first 
 * element is the number of integers that follow.
 * This function supports an abbreviation for consecutive
 * (increasing) integers in the list, using the notation, e.g., 
 * "1-4" as shorthand for "1 2 3 4".
 *
 * Returns: the allocated array, or %NULL on failure.
 */

int *gretl_list_from_string (const char *str, int *err)
{
    char *p, *q, *s, *next;
    int i, r1, r2, rg;
    int *list;
    int n = 0;

    if (str == NULL) {
	*err = E_DATA;
	return NULL;
    }

    /* 'p' marks the memory to be freed */
    p = s = gretl_strdup(str);
    if (s == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *err = 0;

    /* strip white space at both ends */
    while (isspace(*s)) s++;
    tailstrip(s);

    /* strip parentheses, if present */
    if (*s == '(') {
	n = strlen(s);
	if (s[n-1] != ')') {
	    /* got opening paren but no close */
	    *err = E_PARSE;
	    return NULL;
	}
	s[n-1] = '\0';
	s++;
	while (isspace(*s)) s++;
	tailstrip(s);
    }

    q = s; /* copy relevant starting point */

    charsub(s, ',', ' ');

    errno = 0;

    /* first pass: figure out the number of values
       in the list, checking for errors as we go
    */

    n = 0; /* value counter */

    while (*s && !*err) {
	r1 = strtol(s, &next, 10);
	if (errno || next == s) {
	    *err = E_PARSE;
	} else {
	    s = next;
	    if (*s == '-') {
		/* hyphen indicating range? */
		s++;
		r2 = strtol(s, &next, 10);
		if (errno || next == s) {
		    *err = E_PARSE;
		} else if (r2 < r1) {
		    *err = E_PARSE;
		} else {
		    n += r2 - r1 + 1;
		}
		s = next;
	    } else {
		/* single numerical value */
		n++;
	    }
	}
    }

    if (*err || n == 0) {
	free(p);
	return NULL;
    }

    list = gretl_list_new(n);
    if (list == NULL) {
	*err = E_ALLOC;
	free(p);
	return NULL;
    }

    /* second pass: fill out the list (no error
       checking should be needed at this stage) 
    */

    s = q; /* back to start of string */
    n = 1; /* list position indicator */

    while (*s) {
	r1 = strtol(s, &s, 10);
	if (*s == '-') {
	    s++;
	    r2 = strtol(s, &s, 10);
	    rg = r2 - r1 + 1;
	    for (i=0; i<rg; i++) {
		list[n++] = r1 + i;
	    }
	} else {
	    list[n++] = r1;
	}
    } 

    free(p);

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
 * than 9998.
 *
 * Returns: The string representation of the list on success,
 * or %NULL on failure.
 */

char *gretl_list_to_string (const int *list)
{
    char *buf;
    char numstr[8];
    int len, i, err = 0;

    len = 6 * (list[0] + 1);
    if (len > MAXLINE - 32) {
	return NULL;
    }

    buf = malloc(len);
    if (buf == NULL) {
	return NULL;
    }

    *buf = '\0';

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    strcat(buf, " ;");
	} else {
	    sprintf(numstr, " %d", list[i]);
	    strcat(buf, numstr);
	}
    }

    if (err) {
	free(buf);
	buf = NULL;
    }

    return buf;
}

/**
 * gretl_list_to_lags_string:
 * @list: array of integers.
 * @err: location to receive error code.
 * 
 * Prints the given @list of integers into a newly
 * allocated string, separated by commas.  Will fail
 * if the list contains any numbers greater than 998.
 *
 * Returns: The string representation of the list on success,
 * or %NULL on failure.
 */

char *gretl_list_to_lags_string (const int *list, int *err)
{
    char *buf;
    char numstr[8];
    int len, i;

    len = 4 * (list[0] + 1);

    if (len > MAXLINE - 32) {
	*err = E_DATA;
	return NULL;
    }

    buf = malloc(len);
    if (buf == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *buf = '\0';

    for (i=1; i<=list[0]; i++) {
	if (abs(list[i] >= 999)) {
	    *err = E_DATA;
	    break;
	} else {
	    if (i == 1) {
		sprintf(numstr, "%d", list[i]);
	    } else {
		sprintf(numstr, ",%d", list[i]);
	    }
	    strcat(buf, numstr);
	}
    }

    if (*err) {
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
    int i;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == k) {
		return i;
	    }
	}
    }

    return 0;
}

static void reglist_move_const (int *list, int k)
{
    int i, cnum = list[k];

    for (i=k; i>2; i--) {
	list[i] = list[i-1];
    }

    list[2] = cnum;
}

/**
 * reglist_check_for_const:
 * @list: regression list suitable for use with a gretl
 * model (should not contain #LISTSEP).
 * @Z: data array.
 * @pdinfo: dataset information.
 *
 * Checks @list for an intercept term (a variable all of
 * whose valid values in sample are 1).  If such a variable
 * is present, it is moved to position 2 in the list.
 *
 * Returns: 1 if the list contains an intercept, else 0.
 */

int reglist_check_for_const (int *list, const double **Z,
			     const DATAINFO *pdinfo)
{
    int cpos = gretl_list_const_pos(list, 2, Z, pdinfo);
    int ret = 0;

    if (cpos > 1) {
	ret = 1;
    }

    if (cpos > 2) {
	reglist_move_const(list, cpos);
    }

    return ret;
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

    return err;
}

/**
 * gretl_list_purge_const:
 * @list: list of variable ID numbers.
 * @Z: data array.
 * @pdinfo: dataset information.
 *
 * Checks @list from position 1 onward for the presence of a 
 * variable whose valid values in sample all equal 1.0.  If 
 * such a variable is found, it is deleted from @list (that is, 
 * any following elements are moved forward by one and list[0] 
 * is decremented by 1).
 *
 * Returns: 1 if a constant was found and deleted, else 0.
 */

int gretl_list_purge_const (int *list, const double **Z,
			    const DATAINFO *pdinfo)
{
    int i, gotc = 0;
    int l0 = list[0];

    /* handle the case where the constant comes last; if it's
       the only element behind the list separator, remove both
       the constant and the separator */

    if (list[l0] == 0 || true_const(list[l0], Z, pdinfo)) {
	gotc = 1;
	list[0] -= 1;
	if (list[l0 - 1] == LISTSEP) {
	    list[l0 - 1] = 0;
	    list[0] -= 1;
	}
    } else {
	for (i=1; i<l0; i++) {
	    if (list[i] == 0 || true_const(list[i], Z, pdinfo)) {
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

static int *real_gretl_list_union (const int *orig, const int *add, 
				   int *err, int duperr)
{
    int i, j, k;
    int *big;
    const int norig = orig[0];
    const int nadd = add[0];

    *err = 0;

    big = gretl_list_new(norig + nadd);
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
		if (duperr) {
		    free(big);
		    *err = E_ADDDUP;
		    return NULL;
		} else {
		    match = 1;
		}
	    }
	}
	if (!match) {
	    big[0] += 1;
	    big[++k] = add[i];
	}
    }

    if (duperr && big[0] == norig) {
	free(big);
	*err = E_NOADD;
	return NULL;
    }

    return big;
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
    return real_gretl_list_union(orig, add, err, 1);
}

/**
 * gretl_list_union:
 * @l1: list of integers.
 * @l2: list of integers.
 * @err: location to receive error code.
 *
 * Creates a list holding the union of @l1 and @l2.
 *
 * Returns: new list on success, %NULL on error.
 */

int *gretl_list_union (const int *l1, const int *l2, int *err)
{
    return real_gretl_list_union(l1, l2, err, 0);
}

/**
 * gretl_list_intersection:
 * @l1: list of integers.
 * @l2: list of integers.
 * @err: location to receive error code.
 *
 * Creates a list holding the intersection of @l1 and @l2.
 *
 * Returns: new list on success, %NULL on error.
 */

int *gretl_list_intersection (const int *l1, const int *l2, int *err)
{
    int *ret = NULL;
    int i, j;
    int n = 0;

    for (i=1; i<=l1[0]; i++) {
	for (j=1; j<=l2[0]; j++) {
	    if (l2[j] == l1[i]) {
		n++;
		break;
	    }
	}
    }

    if (n == 0) {
	ret = gretl_null_list();
    } else {
	ret = gretl_list_new(n);
	if (ret != NULL) {
	    n = 1;
	    for (i=1; i<=l1[0]; i++) {
		for (j=1; j<=l2[0]; j++) {
		    if (l2[j] == l1[i]) {
			ret[n++] = l1[i];
			break;
		    }
		}
	    }
	}
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    }

    return ret;
}

/**
 * gretl_list_omit_last:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @err: location to receive error code.
 *
 * Creates a list containing all but the last element of @orig.
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

static int *real_gretl_list_omit (const int *orig, const int *omit, 
				  int minpos, int mode, int *err)
{
    const int nomit = omit[0];
    const int norig = list_count(orig);
    int real_nomit = nomit;
    int *smal = NULL;
    int i, j, k;

    *err = 0;

    /* check for spurious "omissions" */
    for (i=1; i<=nomit; i++) {
	int pos = in_gretl_list(orig, omit[i]);

	if (pos < minpos) {
	    if (mode == LIST_OMIT_TIGHT) {
		sprintf(gretl_errmsg, _("Variable %d was not in the original list"),
			omit[i]);
		*err = 1;
		return NULL;
	    } else {
		real_nomit--;
	    }
	}
    }

    if (minpos > 1) {
	/* it's an error to attempt to omit all vars */
	if (nomit == norig - 1) {
	    *err = E_NOVARS;
	    return NULL;
	}
    }

    if (real_nomit == norig) {
	/* omitting all */
	smal = gretl_null_list();
	if (smal == NULL) {
	    *err = E_ALLOC;
	}
	return smal;
    }

    smal = gretl_list_new(norig - real_nomit);
    if (smal == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=1; i<minpos; i++) {
	/* copy unaltered portion */
	smal[i] = orig[i];
    }

    k = minpos;
    for (i=minpos; i<=norig; i++) {
        int match = 0;

        for (j=1; j<=nomit; j++) {
            if (orig[i] == omit[j]) {
                match = 1; /* matching var: omit it */
		break;
            }
        }
        if (!match) { /* var is not in omit list: keep it */
            smal[k++] = orig[i];
        }
    }

    return smal;
}

/**
 * gretl_list_omit:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @omit: list of variables to drop.
 * @minpos: minimum position to check.  This should be 2 for a regular
 * regression list, to skip the dependent var in position 1; but in
 * other contexts it may be 1 to start from the first element of @orig.
 * @err: pointer to receive error code.
 *
 * Creates a list containing the elements of @orig that are not
 * present in @omit.
 *
 * Returns: new list on success, %NULL on error.
 */

int *gretl_list_omit (const int *orig, const int *omit, int minpos, int *err)
{
    return real_gretl_list_omit(orig, omit, minpos, LIST_OMIT_TIGHT, err);
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
    int i, j, k, n;
    int match, err = 0;

    n = biglist[0] - sublist[0];
    targ[0] = n;

    if (n <= 0) {
	err = 1;
    } else {
	k = 1;
	for (i=2; i<=biglist[0]; i++) {
	    match = 0;
	    for (j=2; j<=sublist[0]; j++) {
		if (sublist[j] == biglist[i]) {
		    match = 1;
		    break;
		}
	    }
	    if (!match) {
		if (k <= n) {
		    targ[k++] = biglist[i];
		} else {
		    err = 1;
		}
	    }
	}
    }

    return err;
}

/**
 * gretl_list_diff_new:
 * @biglist: inclusive list.
 * @sublist: subset of biglist.
 * @minpos: position in lists at which to start.
 *
 * Returns: a newly allocated list including the elements of @biglist,
 * from position @minpos onwards, that are not present in @sublist, 
 * again from @minpos onwards, or %NULL on failure.  Note that
 * comparison stops whenever a list separator is found; i.e. only
 * the pre-separator portions of the lists are compared.
 */

int *gretl_list_diff_new (const int *biglist, const int *sublist,
			  int minpos)
{
    int *targ = NULL;
    int i, j, bi;
    int match;

    if (biglist == NULL || sublist == NULL) {
	return NULL;
    }

    targ = gretl_null_list();
    if (targ == NULL) {
	return NULL;
    }

    for (i=minpos; i<=biglist[0]; i++) {
	bi = biglist[i];
	if (bi == LISTSEP) {
	    break;
	}
	match = 0;
	for (j=minpos; j<=sublist[0]; j++) {
	    if (sublist[j] == LISTSEP) {
		break;
	    } else if (sublist[j] == bi) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    /* but is this var already accounted for? */
	    for (j=1; j<=targ[0]; j++) {
		if (targ[j] == bi) {
		    match = 1;
		    break;
		}
	    }
	}
	if (!match) {
	    targ = gretl_list_append_term(&targ, biglist[i]);
	    if (targ == NULL) {
		break;
	    }
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
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_add_list (int **targ, const int *src)
{
    int *big;
    int n1 = (*targ)[0];
    int n2 = src[0];
    int i, err = 0;

    if (targ == NULL || *targ == NULL) {
	return E_DATA;
    }

    if (src == NULL || *src == 0) {
	/* no-op */
	return 0;
    }

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
 * gretl_list_insert_list_minus:
 * @targ: location of list into which @src should be inserted.
 * @src: list to be inserted.
 * @pos: zero-based position at which @src should be inserted.
 *
 * Inserts @src into @targ at @pos.  The length of @targ becomes the
 * sum of the lengths of the two original lists minus one.  This
 * can be useful if we were expecting to insert a single variable
 * but found we had to insert a list instead.  Insertion of @src
 * overwrites any entries in @targ beyond @pos (the expectation is
 * that this function will be called in the process of assembling
 * @targ, in left-to-right mode).
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_insert_list_minus (int **targ, const int *src, int pos)
{
    int *big;
    int n1 = (*targ)[0];
    int n2 = src[0];
    int bign = n1 - 1 + n2;
    int i, err = 0;

    if (pos > n1 + 1) {
	return 1;
    }

    big = realloc(*targ, (bign + 1) * sizeof *big);
    if (big == NULL) {
	err = E_ALLOC;
    } else {
	big[0] = bign;
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
    const char *errmsg = N_("Can't do this: some vars in original "
			    "model have been redefined");
    const char *label;
    char rword[16];
    int j, mc, repl;

    if (ref_id == 0) {
	mc = get_model_count();
    } else {
	mc = ref_id;
    }

    for (j=1; j<=list[0]; j++) {
	if (list[j] == LISTSEP) {
	    continue;
	}
	if (list[j] >= pdinfo->v) {
	    strcpy(gretl_errmsg, _(errmsg));
	    return E_DATA;
	}
	label = VARLABEL(pdinfo, list[j]);
	*rword = '\0';
	sscanf(label, "%15s", rword);
	if (!strcmp(rword, _("Replaced"))) {
	    repl = 0;
	    sscanf(label, "%*s %*s %*s %d", &repl);
	    if (repl >= mc) {
		strcpy(gretl_errmsg, _(errmsg));
		return E_DATA;
	    }
	}
    }

    return 0;
}

/**
 * gretl_list_const_pos:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 * @minpos: position in @list at which to start the search (>= 1).
 *
 * Checks @list for the presence, in position @minpos or higher, of
 * a variable whose valid values in sample all equal 1.  This usually
 * amounts to checking whether a list of regressors includes
 * an intercept term.
 * 
 * Returns: The list position of the const, or 0 if none is
 * found.
 */

int gretl_list_const_pos (const int *list, int minpos, const double **Z, 
			  const DATAINFO *pdinfo)
{
    int i;

    if (minpos < 1) {
	return 0;
    }

    /* we give preference to the "official" const... */
    for (i=minpos; i<=list[0]; i++) {
        if (list[i] == 0) {
	    return i;
	}
    }

    /* ... but if it's not found */
    for (i=minpos; i<=list[0]; i++) {
        if (true_const(list[i], Z, pdinfo)) {
	    return i;
	}
    }

    return 0;
}

/**
 * gretl_list_separator_position:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 *
 * Returns: if @list contains the separator for compound
 * lists, #LISTSEP, the position in @list at which this is found,
 * else 0.  The search begins at position 1.
 */

int gretl_list_separator_position (const int *list)
{
    int i;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) {
		return i;
	    }
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
 * lists, #LISTSEP, else 0.  The search begins at position 1.
 */

int gretl_list_has_separator (const int *list)
{
    return gretl_list_separator_position(list) > 0;
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
    int i, n = 0;
    int err = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    n = i;
	    break;
	}
    }

    if (n == 0) {
	err = E_PARSE;
    }

    if (n > 1) {
	list1 = gretl_list_new(n - 1);
	if (list1 == NULL) {
	    return E_ALLOC;
	}
	for (i=1; i<n; i++) {
	    list1[i] = list[i];
	}
    }

    if (n < list[0]) {
	list2 = gretl_list_new(list[0] - n);
	if (list2 == NULL) {
	    free(list1);
	    return E_ALLOC;
	}

	for (i=1; i<=list2[0]; i++) {
	    list2[i] = list[i + n];
	}
    }

    *plist1 = list1;
    *plist2 = list2;
    
    return err;
}

/**
 * gretl_lists_join_with_separator:
 * @list1: first sub-list.
 * @list2: second sub-list.
 *
 * Concatenates the content of @list2 onto @list1, after first
 * appending the list separator.  It is acceptable that @list1
 * be %NULL, in which case the returned list is just @list2
 * with the separator prepended.  But it is not acceptable that
 * @list2 be null; in that this function returns %NULL.
 *
 * Returns: alllcated list  on success or %NULL on failure.
 */

int *gretl_lists_join_with_separator (const int *list1, const int *list2)
{
    int *biglist;
    int i, j, n;

    if (list2 == NULL) {
	return NULL;
    }

    n = (list1 != NULL)? list1[0] : 0;


    n += list2[0] + 1;
    biglist = gretl_list_new(n);

    if (biglist == NULL) {
	return NULL;
    }

    j = 1;

    if (list1 != NULL) {
	for (i=1; i<=list1[0]; i++) {
	    biglist[j++] = list1[i];
	}
    }

    biglist[j++] = LISTSEP;

    for (i=1; i<=list2[0]; i++) {
	biglist[j++] = list2[i];
    }    

    return biglist;
}

static int real_list_dup (const int *list, int start, int stop)
{
    int i, j, ret = -1;

    for (i=start; i<stop && ret<0; i++) {
	for (j=i+1; j<=stop && ret<0; j++) {
	    if (list[i] == list[j]) {
		ret = list[i];
	    }
	}
    }

    return ret;
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
    int multi = 0;
    int start = 2;
    int i, ret = -1;

    if (ci == COINT) {
	start = 1;
    } else if (ci == ARCH) {
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
    } else if (ci == AR || ci == SCATTERS || ci == MPOLS || ci == GARCH) {
	for (i=2; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
    } else if (ci == IVREG || ci == HECKIT || ci == EQUATION) {
	multi = 1;
	for (i=2; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
	ret = real_list_dup(list, start, list[0]);
	if (ret == -1) {
	    ret = real_list_dup(list, 2, start - 2);
	}
    } else if (ci == VAR || ci == VECM || ci == COINT2) {
	multi = 1;
	for (i=1; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i+1;
		break;
	    }
	}
	ret = real_list_dup(list, start, list[0]);
	if (ret == -1) {
	    ret = real_list_dup(list, 1, start - 2);
	}
    } else if (ci == ARBOND) {
	int stop = 0;

	multi = 1;
	for (i=2; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		start = i;
		break;
	    }
	}
	for (i=list[0]-1; i>=2; i--) {
	    if (list[i] == LISTSEP) {
		stop = i;
		break;
	    }
	}

	if (stop == start) {
	    ret = real_list_dup(list, start + 1, list[0]);
	} else {
	    ret = real_list_dup(list, start + 1, stop - 1);
	    if (ret == -1) {
		ret = real_list_dup(list, stop + 1, list[0]);
	    }
	}
	multi = 1;
    } 

    if (!multi) {
	ret = real_list_dup(list, start, list[0]);
    }

    return ret;
}

/**
 * gretl_lists_share_members:
 * @list1: 
 * @list2: 
 *
 * Returns: 1 if there are any elements in common between @list1
 * and @list2, otherwise 0.
 */

int gretl_lists_share_members (const int *list1, const int *list2)
{
    int i;

    if (list1 == NULL || list2 == NULL) {
	return 0;
    }

    for (i=1; i<=list1[0]; i++) {
	if (in_gretl_list(list2, list1[i])) {
	    return 1;
	}
    }

    return 0;
}

/**
 * full_var_list:
 * @pdinfo: dataset information.
 * @nvars: location for return of number of elements in full list.
 *
 * Creates a newly allocated list including all variables in the
 * dataset that are not scalars, are not hidden variables, and
 * are accessible at the current level of function execution.
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
    int fsd = gretl_function_depth();
    int i, j, nv = 0;
    int *list = NULL;

    for (i=1; i<pdinfo->v; i++) {
	if (!var_is_hidden(pdinfo, i) &&
	    STACK_LEVEL(pdinfo, i) == fsd) {
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
	    if (!var_is_hidden(pdinfo, i) &&
		STACK_LEVEL(pdinfo, i) == fsd) {
		list[j++] = i;
	    }
	}
    }	    

    return list;
}

/**
 * gretl_list_position:
 * @v: ID number of variable to find.
 * @list: list to search.
 *
 * Scans @list, from position 1 onward, for variable @v.
 * In case the list contains a list separator, #LISTSEP,
 * only the portion prior to the separator is considered.
 *
 * Returns: the 0-based position, or 0 if the variable is 
 * not found.
 */

int gretl_list_position (int v, const int *list)
{
    int i, lmax = list[0];

    /* handle lists with separator: search only first part */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    lmax = i - 1;
	    break;
	}
    }
	    
    for (i=lmax; i>=1; i--) {
	if (v == list[i]) return i;
    }

    return 0;
}

/**
 * gretl_list_is_consecutive:
 * @list: list to check.
 *
 * Returns: 1 if the elements of @list, from position 1 onward,
 * are consecutive integer values, else 0.
 */

int gretl_list_is_consecutive (const int *list)
{
    int i, ret = 1;

    for (i=2; i<=list[0]; i++) {
	if (list[i] != list[i-1] + 1) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * gretl_list_build:
 * @s: string list specification.
 * @pdinfo: dataset information.
 * @err: location to receive error code
 *
 * Builds a list based on the specification in @s, which may include
 * the ID numbers of variables, the names of variables, and/or the
 * names of previously defined lists (all separated by spaces).
 *
 * Returns: the constructed list, or %NULL on failure.
 */

int *gretl_list_build (const char *s, const DATAINFO *pdinfo, int *err)
{
    char test[32];
    int *list = NULL;
    int *nlist;
    int i, v, len, nf;

    list = gretl_null_list();
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    nf = count_fields(s);
    
    for (i=0; i<nf && !*err; i++) {
	s += strspn(s, " ");
	len = strcspn(s, " ");
	if (len > 31) {
	    *err = E_PARSE;
	} else {
	    *test = 0;
	    strncat(test, s, len);

	    /* valid elements: integers, varnames, named lists */

	    if (isdigit(*test)) {
		v = positive_int_from_string(test);
		if (v >= 0) {
		    list = gretl_list_append_term(&list, v);
		} else {
		    *err = E_PARSE;
		}
	    } else {
		v = series_index(pdinfo, test);
		if (v < pdinfo->v) {
		    list = gretl_list_append_term(&list, v);
		} else {
		    nlist = get_list_by_name(test);
		    if (nlist != NULL) {
			*err = gretl_list_add_list(&list, nlist);
		    } else {
			*err = E_UNKVAR;
		    }
		}
	    }

	    if (list == NULL) {
		*err = E_ALLOC;
	    }
	}
	s += len;
    }

    if (*err) {
	free(list);
	list = NULL;
    }

    return list;
}

int load_user_lists_file (const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    int i, nl, err = 0;

    xmlKeepBlanksDefault(0);

    err = gretl_xml_open_doc_root(fname, "gretl-lists", &doc, &node);
    if (err) {
	return err;
    }

    if (!gretl_xml_get_prop_as_int(node, "count", &nl)) {
	err = E_DATA;
    } else if (nl <= 0) {
	err = E_DATA;
    }

    if (!err) {
	list_stack = malloc(nl * sizeof *list_stack);
	if (list_stack == NULL) {
	    err = E_ALLOC;
	} else {
	    n_lists = nl;
	    for (i=0; i<nl; i++) {
		list_stack[i] = NULL;
	    }
	    for (i=0; i<nl && !err; i++) {
		list_stack[i] = malloc(sizeof **list_stack);
		if (list_stack[i] == NULL) {
		    err = E_ALLOC;
		} else {
		    list_stack[i]->name[0] = '\0';
		    list_stack[i]->list = NULL;
		    list_stack[i]->level = 0;
		}
	    }
	}
    }

    if (!err) {
	xmlNodePtr cur = node->xmlChildrenNode;
	char *lname;

	i = 0;
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (XUC) "list")) {
		if (!gretl_xml_get_prop_as_string(cur, "name", &lname)) {
		    err = E_DATA;
		} else {
		    strncat(list_stack[i]->name, lname, VNAMELEN - 1);
		    free(lname);
		    list_stack[i]->list = 
			gretl_xml_node_get_list(cur, doc, &err);
		    i++;
		}
	    }
	    cur = cur->next;
	}
    }

    if (err && list_stack != NULL) {
	gretl_lists_cleanup();
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

    return err;
}

/**
 * gretl_serialize_lists:
 * @fname: name of file to which output should be written.
 *
 * Prints an XML representation of the current saved lists,
 * if any.
 *
 * Returns: 0 on success, or if there are no saved lists, 
 * non-zero code on error.
 */

int gretl_serialize_lists (const char *fname)
{
    FILE *fp;
    int i;

    if (n_lists == 0) {
	return 0;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_xml_header(fp); 

    fprintf(fp, "<gretl-lists count=\"%d\">\n", n_lists);

    for (i=0; i<n_lists; i++) {
	gretl_xml_put_named_list(list_stack[i]->name, 
				 list_stack[i]->list, 
				 fp);
    }

    fputs("</gretl-lists>\n", fp);

    fclose(fp);

    return 0;
}

/**
 * gretl_list_print:
 * @lname: name of list.
 * @pdinfo:
 * @prn:
 * 
 * Prints to @prn the given @list of variables, by name.
 */

void gretl_list_print (const char *lname, const DATAINFO *pdinfo,
		       PRN *prn)
{
    const int *list = get_list_by_name(lname);
    int testlen = 62;
    int i, li, len = 0;

    if (list == NULL) {
	pprintf(prn, _("Unknown variable '%s'"), lname);
	pputc(prn, '\n');
    } else if (list[0] == 0) {
	pputs(prn, "null\n");
    } else {
	for (i=1; i<=list[0]; i++) {
	    li = list[i];
	    if (li == LISTSEP) {
		len += pputs(prn, "; ");
	    } else if (li < 0 || li >= pdinfo->v) {
		len += pputs(prn, "?? ");
	    } else {
		len += pprintf(prn, "%s ", pdinfo->varname[li]);
		if (len > testlen && i < list[0]) {
		    pputs(prn, "\\\n "); 
		    len = 1;
		}
	    }
	}
	pputc(prn, '\n');
    }
}

