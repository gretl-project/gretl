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
#include "gretl_midas.h"
#include "libset.h"
#include "uservar.h"

#include <errno.h>
#include <glib.h>

#define LDEBUG 0

/**
 * SECTION:gretl_list
 * @short_description: handling of lists of integers
 * @title: Lists
 * @include: libgretl.h
 *
 * Lists of integers are used in many contexts in libgretl, e.g.
 * for holding the ID numbers of variables in a regression
 * specification. A gretl "list" is simply an array of ints
 * following a definite convention: the value at position 0
 * gives the number of elements that follow. The total number
 * of ints in the list foo is therefore foo[0] + 1, and reading
 * the substantive members of foo involves looping from
 * position 1 to position foo[0].
 */

/**
 * LISTSEP:
 *
 * Symbolic name for the separator used in gretl lists; this
 * corresponds to a semicolon in the string representation of
 * a list.
 */

/**
 * gretl_is_midas_list:
 * @list: the list array.
 * @dset: pointer to dataset.
 *
 * Returns: 1 if @list has been set as a MIDAS list in an
 * approved manner, 0 otherwise.
 */

int gretl_is_midas_list (const int *list, const DATASET *dset)
{
    int ret = 0;

    if (list != NULL && list[0] > 2) {
	int i, m, p;

	if (!series_is_midas_anchor(dset, list[1])) {
	    return 0;
	}

	m = series_get_midas_period(dset, list[1]);
	if (!is_valid_midas_frequency_ratio(dset, m) || list[0] != m) {
	    return 0;
	}

	ret = 1;
	for (i=2; i<=list[0] && ret; i++) {
	    p = series_get_midas_period(dset, list[i]);
	    if (p != m - 1) {
		ret = 0;
	    }
	    m = p;
	}
    }

    return ret;
}

/**
 * gretl_list_set_midas:
 * @list: the list array.
 * @dset: pointer to dataset.
 *
 * Attempts to set the MIDAS flag on the members of @list.
 *
 * 0 on success, non-zero code on failure.
 */

int gretl_list_set_midas (const int *list, DATASET *dset)
{
    int err = 0;

    if (list != NULL) {
	int i, m = list[0];

	if (!is_valid_midas_frequency_ratio(dset, m)) {
	    err = E_INVARG;
	}

	for (i=1; i<=list[0] && !err; i++) {
	    if (list[i] < 1 || list[i] >= dset->v) {
		err = E_INVARG;
	    }
	}

	if (!err) {
	    int freq = get_midas_frequency(dset, m);

	    series_set_midas_anchor(dset, list[1]);
	    series_set_midas_period(dset, list[1], m);
	    series_set_midas_freq(dset, list[1], freq);

	    for (i=2; i<=list[0]; i++) {
		series_set_midas_period(dset, list[i], m - i + 1);
		series_set_midas_freq(dset, list[i], freq);
	    }
	}
    }

    return err;
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
 * Returns: the newly allocated list, or NULL on failure.
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
 * gretl_list_array_new:
 * @nlists: the number of lists to create.
 * @nterms: the maximum number of elements to be stored in each list.
 *
 * Creates an array of newly allocated lists, each of which as described
 * in connection with gretl_list_new().
 *
 * Returns: the newly allocated lists, or NULL on failure.
 */

int **gretl_list_array_new (int nlists, int nterms)
{
    int **lists = NULL;
    int i;

    if (nlists < 0) {
	return NULL;
    }

    lists = malloc(nlists * sizeof *lists);

    if (lists != NULL) {
	for (i=0; i<nlists; i++) {
	    lists[i] = gretl_list_new(nterms);
	}
    }

    return lists;
}

/**
 * gretl_list_array_free:
 * @lists: array of gretl lists.
 * @nlists: the number of lists in @lists.
 *
 * Frees all the lists in @lists and also the top-level pointer.
 */

void gretl_list_array_free (int **lists, int nlists)
{
    if (lists != NULL) {
	int i;

	for (i=0; i<nlists; i++) {
	    free(lists[i]);
	}
	free(lists);
    }
}

/**
 * gretl_consecutive_list_new:
 * @lmin: starting value for consecutive list elements.
 * @lmax: ending value.
 *
 * Creates a newly allocated list whose elements run from
 * @lmin to @lmax consecutively.
 *
 * Returns: the newly allocated list, or NULL on failure.
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
 * Returns: the resized list, or NULL on failure.
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
 * Returns: the augmented list, or NULL on failure.
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

    if (list1 == NULL && list2 != NULL) {
	return 1;
    } else if (list1 != NULL && list2 == NULL) {
	return 1;
    } else if (list1 == NULL && list2 == NULL) {
	return 0;
    }

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
 * Returns: the newly allocated list, or NULL on failure.
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
 * Returns: an allocated copy @src (or NULL if @src is NULL).
 */

int *gretl_list_copy (const int *src)
{
    int *targ = NULL;

    if (src != NULL) {
	int n = src[0] + 1;

	targ = malloc(n * sizeof *targ);
	if (targ != NULL) {
	    memcpy(targ, src, n * sizeof *targ);
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
 * (or NULL on failure).
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
 * Returns: the allocated array, or NULL on failure.
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

    /* strip parentheses or braces, if present */
    if (*s == '(' || *s == '{') {
	char close = (*s == '(')? ')' : '}';

	n = strlen(s);
	if (s[n-1] != close) {
	    /* got opening grouping character but no close */
	    *err = E_PARSE;
	    return NULL;
	}
	s[n-1] = '\0';
	s++;
	while (isspace(*s)) s++;
	tailstrip(s);
    }

    q = s; /* copy relevant starting point */

    gretl_charsub(s, ',', ' ');

    errno = 0;

    /* first pass: figure out the number of values
       in the list, checking for errors as we go
    */

    n = 0; /* value counter */

    while (*s && !*err) {
	s += strspn(s, " ");
	if (n > 0 && *s == ';') {
	    /* list separator */
	    n++;
	    s++;
	    continue;
	}
	r1 = strtol(s, &next, 10);
	if (errno || next == s) {
	    fprintf(stderr, "gretl_list_from_string: '%s'\n", s);
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
	s += strspn(s, " ");
	if (*s == ';') {
	    list[n++] = LISTSEP;
	    s++;
	    continue;
	}
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
 * gretl_list_from_varnames:
 * @str: string holding space-separated series names.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Returns: an allocated gretl list holding the ID numbers of
 * the series named in @str, or NULL on failure.
 */

int *gretl_list_from_varnames (const char *str,
			       const DATASET *dset,
			       int *err)
{
    int *list = NULL;
    char **S;
    int n = 0;

    if (str == NULL || count_fields(str, NULL) < 1) {
	*err = E_DATA;
	return NULL;
    }

    S = gretl_string_split(str, &n, NULL);
    if (S == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    list = gretl_list_new(n);

    if (list == NULL) {
	*err = E_ALLOC;
    } else {
	int i, vi;

	for (i=0; i<n; i++) {
	    if (!strcmp(S[i], "0")) {
		vi = 0;
	    } else {
		vi = current_series_index(dset, S[i]);
	    }
	    if (vi < 0) {
		*err = E_UNKVAR;
		break;
	    } else {
		list[i+1] = vi;
	    }
	}
    }

    strings_array_free(S, n);

    return list;
}

#define LN_10 2.30258509299404590

static int integer_length (int k)
{
    double x = k;
    int len = 0;

    if (x < 0) {
	x = fabs(x);
	len = 1;
    }

    if (x < 10) {
	len += 1;
    } else {
	len += (int) ceil(log(x)/LN_10);
	len += (k % 10 == 0);
    }

    return len;
}

/**
 * gretl_list_to_numeric_string:
 * @list: array of integers.
 *
 * Prints the given @list of integers into a newly
 * allocated string, separated by single spaces and with
 * one leading space.
 *
 * Returns: The string representation of the list on success,
 * or NULL on failure.
 */

char *gretl_list_to_numeric_string (const int *list)
{
    char *buf;
    int i, len = 1;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    len += 2;
	} else {
	    len += integer_length(list[i]) + 1;
	}
    }

    if (len > MAXLINE - 32) {
	/* string would be too long for command line */
	return NULL;
    }

    buf = malloc(len);

    if (buf != NULL) {
	char numstr[16];

	*buf = '\0';
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) {
		strcat(buf, " ;");
	    } else {
		sprintf(numstr, " %d", list[i]);
		strcat(buf, numstr);
	    }
	}
    }

    return buf;
}

/**
 * gretl_list_to_string:
 * @list: array of integers.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Returns: allocated string representation of @list, with ID
 * numbers cashed out as series names (and with one leading
 * space), or NULL on failure. The list separator #LISTSEP,
 * is accepted in the incoming @list, otherwise all terms
 * must be integers in the range 0 to the greatest current
 * series ID within @dset.
 */

char *gretl_list_to_string (const int *list,
			    const DATASET *dset,
			    int *err)
{
    char *buf = NULL;
    int len = 1;
    int i, vi;

    if (list == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (list[0] == 0) {
	return gretl_strdup("");
    }

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (vi == LISTSEP) {
	    len += 2;
	} else if (vi >= 0 && vi < dset->v) {
	    len += strlen(dset->varname[vi]) + 1;
	} else {
	    *err = E_DATA;
	    return NULL;
	}
    }

    buf = calloc(len, 1);

    if (buf == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == LISTSEP) {
		strcat(buf, " ;");
	    } else {
		strcat(buf, " ");
		strcat(buf, dset->varname[vi]);
	    }
	}
    }

    return buf;
}

static int get_excluded_id (const int *list, int n)
{
    int i;

    for (i=1; i<n; i++) {
	if (!in_gretl_list(list, i)) {
	    return i;
	}
    }

    return 0;
}

/**
 * gretl_list_to_compact_string:
 * @list: array of integers.
 * @dset: pointer to dataset.
 * @argstyle: boolean for function-argument variant.
 * @err: location to receive error code.
 *
 * Designed to produce a reasonably compact string representation
 * of the @list argument. If the list is short enough the task is
 * handed off to gretl_list_to_string(). Otherwise we try to use
 * either the "dataset" macro, or the ".." notation for a list of
 * contiguous series IDs. The fallback for a large list, if these
 * tactics do not seem applicable, is to compose a string holding
 * the series ID numbers.
 *
 * If @argstyle is non-zero the returned string will suitable for
 * giving as an argument to the userspace deflist function.
 *
 * Returns: allocated string representation of @list, or NULL
 * on error.
 */

char *gretl_list_to_compact_string (const int *list,
				    const DATASET *dset,
				    int argstyle,
				    int *err)
{
    char *buf = NULL;
    int *mylist = NULL;
    char cstr[8];
    char sep;
    int buflen;
    int nv = 0;      /* number of 'real' series */
    int v0 = -1;     /* first non-zero series ID */
    int vn = -1;     /* last non-zero series ID */
    int contig = 1;  /* contiguity flag */
    int vprev = -1;  /* ID of previous series */
    int idmax = 0;   /* greatest series ID */
    int cpos = 0;    /* position of constant */
    int i, vi;

    if (argstyle) {
	strcpy(cstr, "const,");
	sep = ',';
    } else {
	strcpy(cstr, "const ");
	sep = ' ';
    }

    if (list == NULL) {
	*err = E_DATA;
	return NULL;
    } else if (list[0] == 0) {
	/* empty list */
	return gretl_strdup("");
    } else if (list[0] <= 6) {
	/* short list: don't bother with the fancy stuff */
	buf = gretl_list_to_string(list, dset, err);
	if (*err) {
	    return buf;
	} else if (argstyle) {
	    return gretl_charsub(g_strchug(buf), ' ', ',');
	} else {
	    return g_strchug(buf);
	}
    }

    /* analyse the incoming list */
    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (vi == LISTSEP) {
	    gretl_errmsg_set("Invalid list element");
	    *err = E_INVARG;
	    return NULL;
	} else if (vi == 0) {
	    /* record the const position, if present */
	    cpos = i;
	} else {
	    if (v0 < 0) {
		/* record the first "real" series */
		v0 = vi;
	    }
	    if (vi > idmax) {
		/* highest-numbered series in list */
		idmax = vi;
	    }
	    if (contig && vprev > 0 && vi != vprev + 1) {
		/* list is not contiguous */
		contig = 0;
	    }
	    vprev = vi;
	}
    }

    if (cpos > 0) {
	/* delete the constant */
	mylist = gretl_list_copy(list);
	gretl_list_delete_at_pos(mylist, cpos);
    } else {
	mylist = (int *) list;
    }

    nv = mylist[0];
    vn = mylist[nv];

    if (nv >= dset->v - 2) {
	/* an almost-complete list: use the "dataset" macro,
	   maybe with const prepended and/or one exclusion
	*/
	if (nv == dset->v - 1) {
	    /* no "real" series excluded */
	    buflen = 6 * (cpos > 0) + 9;
	    buf = calloc(buflen, 1);
	    sprintf(buf, "%sdataset", cpos > 0 ? cstr : "");
	} else {
	    /* one "real" series excluded */
	    int vx = get_excluded_id(mylist, dset->v);
	    const char *sx = dset->varname[vx];

	    buflen = 12 + 6 * (cpos > 0) + strlen(sx);
	    buf = calloc(buflen, 1);
	    sprintf(buf, "%sdataset-%s", cpos > 0 ? cstr : "", sx);
	}
    } else if (contig) {
	/* use ".." notation for contiguous series IDs from
	   @v0 to @vn, possibly with const prepended
	*/
	const char *s0 = dset->varname[v0];
	const char *sn = dset->varname[vn];

	buflen = 4 + 6 * (cpos > 0) + strlen(s0) + strlen(sn);
	buf = calloc(buflen, 1);
	sprintf(buf, "%s%s..%s", cpos > 0 ? cstr : "", s0, sn);
    } else {
	/* use straight numeric IDs */
	int numlen = (int) floor(log10((double) idmax));
	char tmp[16];

	buflen = 2 + 2 * (cpos > 0) + nv * (2 + numlen);
	buf = calloc(buflen, 1);
	sprintf(buf, "%s%c", cpos > 0 ? "0" : "", cpos > 0 ? sep : '\0');
	for (i=1; i<=mylist[0]; i++) {
	    sprintf(tmp, "%d%c", mylist[i], sep);
	    strcat(buf, tmp);
	}
	g_strchomp(buf);
    }

    if (mylist != list) {
	free(mylist);
    }

    return buf;
}

/**
 * gretl_list_to_vector:
 * @list: array of integers.
 * @err: location to receive error code.
 *
 * Returns: allocated representation of @list as a row vector
 * or NULL on failure.
 */

gretl_matrix *gretl_list_to_vector (const int *list, int *err)
{
    gretl_vector *v = NULL;

    if (list == NULL) {
	*err = E_DATA;
    } else {
	int i, n = list[0];

	if (n == 0) {
	    v = gretl_null_matrix_new();
	    if (v == NULL) {
		*err = E_ALLOC;
	    }
	} else if (n > 0) {
	    v = gretl_vector_alloc(n);
	    if (v == NULL) {
		*err = E_ALLOC;
	    } else {
		for (i=0; i<n; i++) {
		    v->val[i] = list[i+1];
		}
	    }
	} else {
	    *err = E_DATA;
	}
    }

    return v;
}

/**
 * gretl_list_get_names:
 * @list: array of integers.
 * @dset: dataset information.
 * @err: location to receive error code.
 *
 * Prints the names of the members of @list of integers into
 * a newly allocated string, separated by commas.
 *
 * Returns: allocated string on success or NULL on failure.
 */

char *gretl_list_get_names (const int *list, const DATASET *dset,
			    int *err)
{
    char *buf = NULL;
    int len = 0;
    int i, vi;

    if (list == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (list[0] == 0) {
	return gretl_strdup("");
    }

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (vi < 0 || vi >= dset->v) {
	    len += strlen("unknown") + 1;
	} else {
	    len += strlen(dset->varname[vi]) + 1;
	}
    }

    buf = malloc(len);
    if (buf == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *buf = '\0';

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (vi < 0 || vi >= dset->v) {
	    strcat(buf, "unknown");
	} else {
	    strcat(buf, dset->varname[vi]);
	}
	if (i < list[0]) {
	    strcat(buf, ",");
	}
    }

    return buf;
}

/**
 * gretl_list_get_names_array:
 * @list: array of integers.
 * @dset: dataset information.
 * @err: location to receive error code.
 *
 * Returns: An array of strings holding the names of the
 * members of @list, or NULL on failure.
 */

char **gretl_list_get_names_array (const int *list,
				   const DATASET *dset,
				   int *err)
{
    char **S = NULL;
    int i, vi, n;

    if (list == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (list[0] == 0) {
	return NULL;
    }

    n = list[0];

    S = strings_array_new(n);
    if (S == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<n; i++) {
	vi = list[i+1];
	if (vi < 0 || vi >= dset->v) {
	    S[i] = gretl_strdup("unknown");
	} else {
	    S[i] = gretl_strdup(dset->varname[vi]);
	}
	if (S[i] == NULL) {
	    *err = E_ALLOC;
	    strings_array_free(S, n);
	    S = NULL;
	    break;
	}
    }

    return S;
}

/**
 * gretl_list_to_lags_string:
 * @list: array of integers.
 * @err: location to receive error code.
 *
 * Prints the given @list of integers into a newly
 * allocated string, enclosed by braces and separated by commas.
 * Will fail if @list contains any numbers greater than 998.
 *
 * Returns: The string representation of the list on success,
 * or NULL on failure.
 */

char *gretl_list_to_lags_string (const int *list, int *err)
{
    char *buf;
    char numstr[16];
    int len, i;

    len = 4 * (list[0] + 1) + 2;

    if (len > MAXLINE - 32) {
	*err = E_DATA;
	return NULL;
    }

    buf = calloc(len, 1);
    if (buf == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=1; i<=list[0]; i++) {
	if (abs(list[i] >= 999)) {
	    *err = E_DATA;
	    break;
	} else {
	    if (i == 1) {
		sprintf(numstr, "{%d", list[i]);
	    } else {
		sprintf(numstr, ",%d", list[i]);
	    }
	    strcat(buf, numstr);
	}
    }
    strcat(buf, "}");

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
 * @dset: dataset struct.
 *
 * Checks @list for an intercept term (a variable all of
 * whose valid values in sample are 1).  If such a variable
 * is present, it is moved to position 2 in the list.
 *
 * Returns: 1 if the list contains an intercept, else 0.
 */

int reglist_check_for_const (int *list, const DATASET *dset)
{
    int cpos = gretl_list_const_pos(list, 2, dset);
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
 * @dset: dataset struct.
 *
 * Checks @list from position 1 onward for the presence of a
 * variable whose valid values in sample all equal 1.0.  If
 * such a variable is found, it is deleted from @list (that is,
 * any following elements are moved forward by one and list[0]
 * is decremented by 1).
 *
 * Returns: 1 if a constant was found and deleted, else 0.
 */

int gretl_list_purge_const (int *list, const DATASET *dset)
{
    int i, gotc = 0;
    int l0 = list[0];

    /* handle the case where the constant comes last; if it's
       the only element behind the list separator, remove both
       the constant and the separator */

    if (list[l0] == 0 || true_const(list[l0], dset)) {
	gotc = 1;
	list[0] -= 1;
	if (list[l0 - 1] == LISTSEP) {
	    list[l0 - 1] = 0;
	    list[0] -= 1;
	}
    } else {
	for (i=1; i<l0; i++) {
	    if (list[i] == 0 || true_const(list[i], dset)) {
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
 * gretl_list_min_max:
 * @list: gretl list.
 * @lmin: location to receive minimum value.
 * @lmax: location to receive maximum value.
 *
 * Reads @list from position 1 onward and writes to @lmin
 * and @lmax the minimum and maximum values among the elements
 * of the list. Reading stops at the end of the list or when
 * a list separator is encountered.
 *
 * Returns: 0 on successful completion, error code if the
 * list is NULL or empty.
 */

int gretl_list_min_max (const int *list, int *lmin, int *lmax)
{
    if (list == NULL || list[0] == 0) {
	return E_DATA;
    } else {
	int i;

	*lmin = *lmax = list[1];

	for (i=2; i<=list[0]; i++) {
	    if (list[i] < *lmin) {
		*lmin = list[i];
	    }
	    if (list[i] > *lmax) {
		*lmax = list[i];
	    }
	}

	return 0;
    }
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
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_add (const int *orig, const int *add, int *err)
{
    int n_orig = orig[0];
    int n_add = add[0];
    int i, j, k;
    int *big;

    *err = 0;

    big = gretl_list_new(n_orig + n_add);
    if (big == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<=n_orig; i++) {
	big[i] = orig[i];
    }

    k = orig[0];

    for (i=1; i<=n_add; i++) {
	for (j=1; j<=n_orig; j++) {
	    if (add[i] == orig[j]) {
		/* a "new" var was already present */
		free(big);
		*err = E_ADDDUP;
		return NULL;
	    }
	}
	big[0] += 1;
	big[++k] = add[i];
    }

    if (big[0] == n_orig) {
	free(big);
	big = NULL;
	*err = E_NOADD;
    }

    return big;
}

/**
 * gretl_list_plus:
 * @l1: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @l2: list of variables to be added.
 * @err: location to receive error code.
 *
 * Creates a list containing all elements of @l1 followed
 * by all elements of @l2. This differs from gretl_list_union()
 * in that some elements may end up being repeated in the
 * returned list.
 *
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_plus (const int *l1, const int *l2, int *err)
{
    int n1 = l1[0];
    int n2 = l2[0];
    int i, j;
    int *ret;

    ret = gretl_list_new(n1 + n2);
    if (ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    j = 1;

    for (i=1; i<=n1; i++) {
	ret[j++] = l1[i];
    }

    for (i=1; i<=n2; i++) {
	ret[j++] = l2[i];
    }

    return ret;
}

static int list_mark_duplicates (const int *l1, int *l2)
{
    int n_dups = 0;
    int i, j;

    for (i=1; i<=l2[0]; i++) {
	if (l2[i] == -1) {
	    continue;
	}
	if (in_gretl_list(l1, l2[i])) {
	    /* element already present in l1 */
	    n_dups++;
	    l2[i] = -1;
	} else {
	    /* not present in @l1, but check for duplicates of
	       this element in @l2
	    */
	    for (j=1; j<=l2[0]; j++) {
		if (j != i && l2[j] == l2[i]) {
		    n_dups++;
		    l2[j] = -1;
		}
	    }
	}
    }

    return n_dups;
}

/**
 * gretl_list_union:
 * @l1: list of integers.
 * @l2: list of integers.
 * @err: location to receive error code.
 *
 * Creates a list holding the union of @l1 and @l2.
 *
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_union (const int *l1, const int *l2, int *err)
{
    int *ret, *lcopy;
    int n_orig = l1[0];
    int n_add = l2[0];

    *err = 0;

    lcopy = gretl_list_copy(l2);

    if (lcopy == NULL) {
	*err = E_ALLOC;
	return NULL;
    } else {
	n_add -= list_mark_duplicates(l1, lcopy);
    }

    if (n_add == 0) {
	ret = gretl_list_copy(l1);
    } else {
	ret = gretl_list_new(n_orig + n_add);
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    } else if (n_add > 0) {
	int i, j;

	for (i=1; i<=n_orig; i++) {
	    ret[i] = l1[i];
	}
	j = n_orig + 1;
	for (i=1; i<=lcopy[0]; i++) {
	    if (lcopy[i] != -1) {
		ret[j++] = lcopy[i];
	    }
	}
    }

    free(lcopy);

    return ret;
}

/* return a copy of @list from which any duplicated elements
   have been deleted
*/

static int *gretl_list_uniq (const int *list)
{
    int *ret = gretl_list_copy(list);
    int i, j;

    for (j=list[0]; j>0; j--) {
	for (i=1; i<j; i++) {
	    if (list[j] == list[i]) {
		gretl_list_delete_at_pos(ret, j);
	    }
	}
    }

    return ret;
}

/**
 * gretl_list_merge_list:
 * @targ: pointer to original list.
 * @src: list to merge.
 *
 * Appends to *targ any unique elements of @src that are not
 * already present in *targ. @targ must not be NULL, but it's
 * OK if the list to which it points is NULL.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_merge_list (int **targ, const int *src)
{
    int *tmp = NULL;
    int *scpy = NULL;
    int n1 = 0, n2 = 0;
    int i, j;
    int err = 0;

    if (targ == NULL) {
	return E_INVARG;
    } else if (src == NULL || src[0] == 0) {
	/* no-op */
	return 0;
    }

    if (src[0] > 1) {
	src = scpy = gretl_list_uniq(src);
    }

    if (*targ == NULL) {
	/* special case */
	*targ = gretl_list_copy(src);
	if (*targ == NULL) {
	    err = E_ALLOC;
	}
	free(scpy);
	return err;
    }

    n1 = (*targ)[0];
    for (i=1; i<=src[0]; i++) {
	if (!in_gretl_list(*targ, src[i])) {
	    n2++;
	}
    }
    if (n2 == 0) {
	/* nothing to be appended */
	free(scpy);
	return 0;
    }

    tmp = realloc(*targ, (n1 + n2 + 1) * sizeof *tmp);
    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	j = n1 + 1;
	for (i=1; i<=src[0]; i++) {
	    if (!in_gretl_list(tmp, src[i])) {
		tmp[j++] = src[i];
	    }
	}
	tmp[0] += n2;
	*targ = tmp;
    }

    free(scpy);

    return err;
}

/**
 * gretl_list_intersection:
 * @l1: list of integers.
 * @l2: list of integers.
 * @err: location to receive error code.
 *
 * Creates a list holding the intersection of @l1 and @l2.
 *
 * Returns: new list on success, NULL on error.
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

static void name_xprod_term (char *vname, int vi, int vj,
			     int di, const DATASET *dset)
{
    const char *si = dset->varname[vi];
    const char *sj = dset->varname[vj];
    int ilen = strlen(si);
    int jlen = strlen(sj);
    int totlen = ilen + jlen + 2;
    char numstr[16];

    sprintf(numstr, "%d", di);
    totlen += strlen(numstr);

    if (totlen >= VNAMELEN) {
	int decr = 1 + totlen - VNAMELEN;

	while (decr > 0) {
	    if (ilen > jlen) {
		ilen--;
	    } else {
		jlen--;
	    }
	    decr--;
	}
    }

    sprintf(vname, "%.*s_%.*s_%s", ilen, si, jlen, sj, numstr);
}

static void set_xprod_label (int v, int vi, int vj,
			     double val, DATASET *dset)
{
    const char *si = dset->varname[vi];
    const char *sj = dset->varname[vj];
    char label[MAXLABEL];

    sprintf(label, "interaction of %s and (%s == %g)", si, sj, val);
    series_record_label(dset, v, label);
}

static int nonneg_integer_series (const DATASET *dset, int v)
{
    double xt;
    int t;

    for (t=dset->t1; t<=dset->t2; t++) {
	xt = dset->Z[v][t];
	if (!na(xt) && (xt != floor(xt) || xt < 0)) {
	    return 0;
	}
    }

    return 1;
}

/**
 * gretl_list_product:
 * @X: list of integers (representing discrete variables).
 * @Y: list of integers.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Creates a list holding the Cartesian product of @X and @Y.
 *
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_product (const int *X, const int *Y,
			 DATASET *dset, int *err)
{
    int *ret = NULL;
    gretl_matrix *xvals;
    char vname[VNAMELEN];
    const double *x, *y;
    int *x_is_int = NULL;
    int newv, n_old = 0;
    int n, vi, vj;
    int i, j, k, t;

    if (X == NULL || Y == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (X[0] == 0 || Y[0] == 0) {
	ret = gretl_null_list();
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
	return ret;
    }

    x_is_int = gretl_list_new(X[0]);
    if (x_is_int == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* check the X list for discreteness */

    for (j=1; j<=X[0] && !*err; j++) {
	vj = X[j];
	if (nonneg_integer_series(dset, vj)) {
	    x_is_int[j] = 1;
	} else if (!series_is_discrete(dset, vj)) {
	    gretl_errmsg_sprintf(_("The variable '%s' is not discrete"),
				 dset->varname[vj]);
	    *err = E_DATA;
	}
    }

    if (*err) {
	free(x_is_int);
	return NULL;
    }

    n = sample_size(dset);
    newv = dset->v;

    for (j=1; j<=X[0] && !*err; j++) {
	vj = X[j];
	x = dset->Z[vj];
	xvals = gretl_matrix_values(x + dset->t1, n, OPT_S, err);
	if (!*err) {
	    *err = dataset_add_series(dset, Y[0] * xvals->rows);
	    if (!*err) {
		for (i=1; i<=Y[0] && !*err; i++) {
		    vi = Y[i];
		    y = dset->Z[vi];
		    for (k=0; k<xvals->rows && !*err; k++) {
			int v, oldv, iik;
			double xik;

			xik = gretl_vector_get(xvals, k);
			iik = x_is_int[j] ? (int) xik : (k + 1);
			name_xprod_term(vname, vi, vj, iik, dset);
			oldv = current_series_index(dset, vname);
			if (oldv > 0) {
			    /* reuse existing series of the same name */
			    v = oldv;
			    n_old++;
			} else {
			    /* make a new series */
			    v = newv++;
			}
			for (t=dset->t1; t<=dset->t2; t++) {
			    if (na(x[t]) || na(xik)) {
				dset->Z[v][t] = NADBL;
			    } else {
				dset->Z[v][t] = (x[t] == xik)? y[t] : 0;
			    }
			}
			gretl_list_append_term(&ret, v);
			if (ret == NULL) {
			    *err = E_ALLOC;
			} else {
			    if (v != oldv) {
				strcpy(dset->varname[v], vname);
			    }
			    set_xprod_label(v, vi, vj, xik, dset);
			}
		    }
		}
	    }
	    gretl_matrix_free(xvals);
	}
    }

    free(x_is_int);

    if (n_old > 0) {
	/* we added more series than were actually needed */
	dataset_drop_last_variables(dset, n_old);
    }

    return ret;
}

/**
 * gretl_list_omit_last:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @err: location to receive error code.
 *
 * Creates a list containing all but the last element of @orig,
 * which must not contain #LISTSEP and must contain at least
 * two members.
 *
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_omit_last (const int *orig, int *err)
{
    int *list = NULL;
    int i;

    *err = 0;

    if (orig[0] < 2) {
	*err = E_NOVARS;
    } else {
	for (i=1; i<=orig[0]; i++) {
	    if (orig[i] == LISTSEP) {
		/* can't handle compound lists */
		*err = 1;
		break;
	    }
	}
    }

    if (!*err) {
	list = gretl_list_new(orig[0] - 1);
	if (list == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=1; i<orig[0]; i++) {
		list[i] = orig[i];
	    }
	}
    }

    return list;
}

static int list_count (const int *list)
{
    int i, k = 0;

    if (list == NULL) {
	return 0;
    }

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
 * @minpos: minimum position to check. This should be 2 for a regular
 * regression list, to skip the dependent var in position 1; but in
 * other contexts it may be 1 to start from the first element of @orig.
 * @err: pointer to receive error code.
 *
 * Creates a list containing the elements of @orig that are not
 * present in @omit.  It is an error if the @omit list contains
 * members that are not present in @orig, and also if the @omit
 * list contains duplicated elements.
 *
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_omit (const int *orig, const int *omit,
		      int minpos, int *err)
{
    int n_omit = omit[0];
    int n_orig = list_count(orig);
    int *ret = NULL;
    int i, j, k;

    if (n_omit > n_orig) {
	*err = E_DATA;
	return NULL;
    }

    *err = 0;

    /* check for spurious "omissions" */
    for (i=1; i<=omit[0]; i++) {
	k = in_gretl_list(orig, omit[i]);
	if (k < minpos) {
	    gretl_errmsg_sprintf(_("Variable %d was not in the original list"),
				 omit[i]);
	    *err = 1;
	    return NULL;
	}
    }

    ret = gretl_list_new(n_orig - n_omit);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else if (n_omit < n_orig) {
	int match;

	k = 1;
	for (i=1; i<=n_orig; i++) {
	    if (i < minpos) {
		ret[k++] = orig[i];
	    } else {
		match = 0;
		for (j=1; j<=n_omit; j++) {
		    if (orig[i] == omit[j]) {
			/* matching var: omit it */
			match = 1;
			break;
		    }
		}
		if (!match) {
		    /* var is not in omit list: keep it */
		    ret[k++] = orig[i];
		}
	    }
	}
    }

    return ret;
}

/**
 * gretl_list_drop:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @drop: list of variables to drop.
 * @err: pointer to receive error code.
 *
 * Creates a list containing the elements of @orig that are not
 * present in @drop.  Unlike gretl_list_omit(), processing always
 * starts from position 1 in @orig, and it is not an error if
 * some members of @drop are not present in @orig, or if some
 * members of @drop are duplicated.
 *
 * Returns: new list on success, NULL on error.
 */

int *gretl_list_drop (const int *orig, const int *drop, int *err)
{
    int *lcopy = NULL;
    int *ret = NULL;
    int n_omit = 0;
    int i, k;

    *err = 0;

    lcopy = gretl_list_copy(orig);

    if (lcopy == NULL) {
	*err = E_ALLOC;
	return NULL;
    } else {
	/* see how many terms we're omitting */
	for (i=1; i<=drop[0]; i++) {
	    k = in_gretl_list(lcopy, drop[i]);
	    if (k > 0) {
		n_omit++;
		lcopy[k] = -1;
	    }
	}
    }

    if (n_omit == 0) {
	ret = lcopy;
    } else {
	ret = gretl_list_new(orig[0] - n_omit);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else if (n_omit < orig[0]) {
	    k = 1;
	    for (i=1; i<=orig[0]; i++) {
		if (lcopy[i] >= 0) {
		    ret[k++] = orig[i];
		}
	    }
	}
	free(lcopy);
    }

    return ret;
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
 * again from @minpos onwards, or NULL on failure.  Note that
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
 * gretl_list_append_list:
 * @targ: location of list to which @src should be added.
 * @src: list to be added to @targ.
 *
 * Adds @src onto the end of @targ.  The length of @targ becomes the
 * sum of the lengths of the two original lists. Note that the result
 * may contain duplicated elements.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_append_list (int **targ, const int *src)
{
    int *tmp;
    int i, n1, n2;
    int err = 0;

    if (targ == NULL || *targ == NULL) {
	return E_INVARG;
    } else if (src == NULL || src[0] == 0) {
	/* no-op */
	return 0;
    }

    n1 = (*targ)[0];
    n2 = src[0];

    tmp = realloc(*targ, (n1 + n2 + 1) * sizeof *tmp);

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	tmp[0] = n1 + n2;
	for (i=1; i<=src[0]; i++) {
	    tmp[n1 + i] = src[i];
	}
	*targ = tmp;
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
 * gretl_list_sublist:
 * @list: the source list.
 * @pos0: the starting position.
 * @pos1: the ending position.
 *
 * Returns: a newly allocated sublist containing elements @pos0
 * to @pos1 of the source.
 */

int *gretl_list_sublist (const int *list, int pos0, int pos1)
{
    int n = pos1 - pos0 + 1;
    int *ret = gretl_list_new(n);

    if (n > 0 && ret != NULL) {
	int i, j = 1;

	for (i=pos0; i<=pos1; i++) {
	    ret[j++] = list[i];
	}
    }

    return ret;
}

/**
 * gretl_list_select:
 * @list: the source list.
 * @sel: the selection list.
 *
 * Returns: a newly allocated list containing the elements
 * of @list that are selected by @sel.
 */

int *gretl_list_select (const int *list, const int *sel)
{
    int *ret = gretl_list_new(sel[0]);
    int i, j = 1;

    for (i=1; i<=sel[0]; i++) {
        ret[j++] = list[sel[i]];
    }

    return ret;
}

/**
 * list_members_replaced:
 * @pmod: the model whose list is to be tested.
 * @dset: dataset information.
 *
 * Checks whether any variable used in @pmod has been redefined
 * since the model in question was estimated.
 *
 * Returns: non-zero if any variables have been replaced, 0 otherwise.
 */

int list_members_replaced (const MODEL *pmod, const DATASET *dset)
{
    const char *errmsg = N_("Can't do this: some vars in original "
			    "model have been redefined");
    int i, vi;

    if (pmod->list == NULL) {
	return 0;
    }

    for (i=1; i<=pmod->list[0]; i++) {
	vi = pmod->list[i];
	if (vi == LISTSEP) {
	    continue;
	}
	if (vi >= dset->v) {
	    gretl_errmsg_set(_(errmsg));
	    return E_DATA;
	}
	if (series_get_mtime(dset, vi) > pmod->esttime) {
	    gretl_errmsg_set(_(errmsg));
	    return E_DATA;
	}
    }

    return 0;
}

/**
 * gretl_list_const_pos:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 * @minpos: position in @list at which to start the search (>= 1).
 * @dset: dataset struct.
 *
 * Checks @list for the presence, in position @minpos or higher, of
 * a variable whose valid values in sample all equal 1.  This usually
 * amounts to checking whether a list of regressors includes
 * an intercept term.
 *
 * Returns: The list position of the const, or 0 if none is
 * found.
 */

int gretl_list_const_pos (const int *list, int minpos,
			  const DATASET *dset)
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
        if (true_const(list[i], dset)) {
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
 * @plist1: pointer to accept first sub-list, or NULL.
 * @plist2: pointer to accept second sub-list, or NULL.
 *
 * If @list contains the list separator, #LISTSEP, creates two
 * sub-lists, one containing the elements of @list preceding
 * the separator and one containing the elements following
 * the separator.  The sub-lists are newly allocated, and assigned
 * as the content of @plist1 and @plist2 respectively. Note, however,
 * that one or other of the sublists can be discarded by passing
 * NULL as the second or third argument.
 *
 * Returns: 0 on success, %E_ALLOC is memory allocation fails,
 * or %E_DATA if @list does not contain a separator.
 */

int gretl_list_split_on_separator (const int *list,
				   int **plist1,
				   int **plist2)
{
    int *list1 = NULL, *list2 = NULL;
    int i, n = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    n = i;
	    break;
	}
    }

    if (n == 0) {
	return E_PARSE;
    }

    if (plist1 != NULL) {
	if (n > 1) {
	    list1 = gretl_list_new(n - 1);
	    if (list1 == NULL) {
		return E_ALLOC;
	    }
	    for (i=1; i<n; i++) {
		list1[i] = list[i];
	    }
	}
	*plist1 = list1;
    }

    if (plist2 != NULL) {
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
	*plist2 = list2;
    }

    return 0;
}

/**
 * gretl_lists_join_with_separator:
 * @list1: first sub-list.
 * @list2: second sub-list.
 *
 * Concatenates the content of @list2 onto @list1, after first
 * appending the list separator.  It is acceptable that @list1
 * be NULL, in which case the returned list is just @list2
 * with the separator prepended.  But it is not acceptable that
 * @list2 be null; in that this function returns NULL.
 *
 * Returns: alllcated list on success or NULL on failure.
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

    if (ci == COINT || ci == ANOVA || ci == DELEET) {
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
	int seppos = 0;

	multi = 1;
	for (i=1; i<list[0]; i++) {
	    if (list[i] == LISTSEP) {
		seppos = i;
		break;
	    }
	}
	if (seppos) {
	    /* check each sublist */
	    ret = real_list_dup(list, 1, seppos - 1);
	    if (ret == -1) {
		ret = real_list_dup(list, seppos + 1, list[0]);
	    }
	} else {
	    /* just one list to examine */
	    ret = real_list_dup(list, 1, list[0]);
	}
    } else if (ci == DPANEL) {
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
    } else if (ci == BIPROBIT) {
	multi = 1;
	if (list[1] == list[2]) {
	    ret = 1;
	}
	if (ret == -1) {
	    for (i=3; i<list[0]; i++) {
		if (list[i] == LISTSEP) {
		    start = i+1;
		    break;
		}
	    }
	    ret = real_list_dup(list, start, list[0]);
	    if (ret == -1) {
		ret = real_list_dup(list, 3, start - 2);
	    }
	}
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
 * Returns: the number of elements that are in common between
 * @list1 and @list2.
 */

int gretl_lists_share_members (const int *list1, const int *list2)
{
    int i, n = 0;

    if (list1 != NULL && list2 != NULL) {
	for (i=1; i<=list1[0]; i++) {
	    if (in_gretl_list(list2, list1[i])) {
		n++;
	    }
	}
    }

    return n;
}

/**
 * gretl_list_n_distinct_members:
 * @list: list to test.
 *
 * Returns: the count of distinct elements in list from position
 * 1 onward, not counting #LISTSEP if present.
 */

int gretl_list_n_distinct_members (const int *list)
{
    int i, j, n = list[0];

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    n--;
	} else {
	    for (j=2; j<i; j++) {
		if (list[i] == list[j]) {
		    n--;
		    break;
		}
	    }
	}
    }

    return n;
}

/**
 * full_var_list:
 * @dset: dataset information.
 * @nvars: location for return of number of elements in full list.
 *
 * Creates a newly allocated list including all series in the
 * dataset that are not hidden variables, and are accessible
 * at the current level of function execution.
 * The return value is NULL in case either (a) allocation of
 * memory failed, or (b) the resulting list would be empty.
 * The caller can distinguish between these possibilities by
 * examining the value returned in @nvars, which will be zero if
 * and only if the resulting list would be empty.  If this is
 * not of interest to the caller, @nvars may be given as NULL.
 *
 * Returns: the allocated list, or NULL.
 */

int *full_var_list (const DATASET *dset, int *nvars)
{
    int fsd = gretl_function_depth();
    int i, j, nv = 0;
    int *list = NULL;

    if (dset == NULL) {
	if (nvars != NULL) {
	    *nvars = 0;
	}
	return NULL;
    }

    for (i=1; i<dset->v; i++) {
	if (!series_is_hidden(dset, i) &&
	    series_get_stack_level(dset, i) == fsd) {
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
	for (i=1; i<dset->v; i++) {
	    if (!series_is_hidden(dset, i) &&
		series_get_stack_level(dset, i) == fsd) {
		list[j++] = i;
	    }
	}
    }

    return list;
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
 * @dset: dataset information.
 * @err: location to receive error code
 *
 * Builds a list based on the specification in @s, which may include
 * the ID numbers of variables, the names of variables, and/or the
 * names of previously defined lists (all separated by spaces).
 *
 * Returns: the constructed list, or NULL on failure.
 */

int *gretl_list_build (const char *s, const DATASET *dset, int *err)
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

    nf = count_fields(s, NULL);

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
		v = series_index(dset, test);
		if (v < dset->v) {
		    list = gretl_list_append_term(&list, v);
		} else {
		    nlist = get_list_by_name(test);
		    if (nlist != NULL) {
			*err = gretl_list_append_list(&list, nlist);
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

/**
 * gretl_list_print:
 * @list: list to print.
 * @dset: dataset information.
 * @prn: gretl printing struct.
 *
 * Prints to @prn the given @list of variables, by name
 * if @dset is non-NULL otherwise by ID number.
 */

void gretl_list_print (const int *list, const DATASET *dset,
		       PRN *prn)
{
    int testlen = 62;
    int i, li, len = 0;

    if (list == NULL) {
	pputs(prn, "null\n");
    } else if (list[0] == 0) {
	pputs(prn, "empty\n");
    } else {
	for (i=1; i<=list[0]; i++) {
	    li = list[i];
	    if (li == LISTSEP) {
		len += pputs(prn, "; ");
	    } else if (dset == NULL) {
		len += pprintf(prn, "%d ", li);
	    } else if (li < 0 || li >= dset->v) {
		len += pputs(prn, "?? ");
	    } else {
		len += pprintf(prn, "%s ", dset->varname[li]);
		if (len > testlen && i < list[0]) {
		    pputs(prn, "\\\n ");
		    len = 1;
		}
	    }
	}
	pputc(prn, '\n');
    }
}

/**
 * varname_match_list:
 * @dset: pointer to dataset information.
 * @pattern: pattern to be matched.
 * @err: location to receive error code.
 *
 * Returns: a list of ID numbers of variables whose names
 * match @pattern, or NULL if there are no matches.
 */

int *varname_match_list (const DATASET *dset, const char *pattern,
			 int *err)
{
    GPatternSpec *pspec;
    int *list = NULL;
    int i, fd, n = 0;

    if (dset == NULL || dset->v == 0) {
	return NULL;
    }

    fd = gretl_function_depth();

    pspec = g_pattern_spec_new(pattern);

    for (i=1; i<dset->v; i++) {
	if (fd == 0 || fd == series_get_stack_level(dset, i)) {
	    if (g_pattern_match_string(pspec, dset->varname[i])) {
		n++;
	    }
	}
    }

    if (n > 0) {
	list = gretl_list_new(n);
	if (list == NULL) {
	    *err = E_ALLOC;
	} else {
	    int j = 1;

	    for (i=1; i<dset->v; i++) {
		if (fd == 0 || fd == series_get_stack_level(dset, i)) {
		    if (g_pattern_match_string(pspec, dset->varname[i])) {
			list[j++] = i;
		    }
		}
	    }
	}
    }

    g_pattern_spec_free(pspec);

    return list;
}

/**
 * ellipsis_list:
 * @dset: pointer to dataset information.
 * @v1: index of first variable.
 * @v2: index of last variable.
 * @err: location to receive error code.
 *
 * Returns: a list of ID numbers of variables running
 * from @v1 to @v2.
 */

int *ellipsis_list (const DATASET *dset, int v1, int v2, int *err)
{
    int *list = NULL;
    int i, fd, n = 0;

    if (dset == NULL || dset->v == 0) {
	return NULL;
    }

    fd = gretl_function_depth();

    for (i=v1; i<=v2; i++) {
	if (fd == 0 || fd == series_get_stack_level(dset, i)) {
	    n++;
	}
    }

    if (n > 0) {
	list = gretl_list_new(n);
	if (list == NULL) {
	    *err = E_ALLOC;
	} else {
	    int j = 1;

	    for (i=v1; i<=v2; i++) {
		if (fd == 0 || fd == series_get_stack_level(dset, i)) {
		    list[j++] = i;
		}
	    }
	}
    }

    return list;
}

/**
 * gretl_list_from_vector:
 * @v: matrix (must be a vector).
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Tries to interpret the matrix @v as a list of ID
 * numbers of series. This can work only if @v is a
 * vector, and all its elements have integer values
 * k satisfying 0 <= k < v, where v is the number
 * of series in @dset. In the special case where @v
 * is a null matrix, an empty list is returned.
 *
 * Returns: a gretl list, or NULL on failure.
 */

int *gretl_list_from_vector (const gretl_matrix *v,
			     const DATASET *dset,
			     int *err)
{
    int *list = NULL;

    if (gretl_is_null_matrix(v)) {
	list = gretl_null_list();
	if (list == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	int i, vi, k = gretl_vector_get_length(v);

	if (k == 0) {
	    *err = E_TYPES;
	} else {
	    for (i=0; i<k && !*err; i++) {
		vi = gretl_int_from_double(v->val[i], err);
		if (!*err && (vi >= dset->v || (vi < 0 && vi != LISTSEP))) {
		    gretl_errmsg_sprintf("list from vector: series ID %d "
					 "is out of bounds", vi);
		    *err = E_UNKVAR;
		}
	    }
	    if (!*err) {
		list = gretl_list_new(k);
		if (list == NULL) {
		    *err = E_ALLOC;
		} else {
		    for (i=0; i<k; i++) {
			list[i+1] = (int) v->val[i];
		    }
		}
	    }
	}
    }

    return list;
}

static int colnames_ok_for_series (const char **S, int n)
{
    int i, j, err = 0;

    for (i=0; i<n; i++) {
	err = check_varname(S[i]);
	if (!err && gretl_is_user_var(S[i])) {
	    gretl_errmsg_sprintf("'%s': name conflicts with different type", S[i]);
	    err = E_TYPES;
	} else {
	    for (j=0; j<n; j++) {
		if (j != i && strcmp(S[j], S[i]) == 0) {
		    gretl_errmsg_sprintf("'%s': name is not unique", S[i]);
		    err = 1;
		}
	    }
	}
	if (err) {
	    return 0;
	}
    }

    return 1;
}

static int try_list_vname (char *chkname,
			   const char *pfx,
			   int j, int n)
{
    gchar *tmp = g_strdup_printf("%s%0*d", pfx, n, j);
    int err = check_varname(tmp);

    if (!err && gretl_is_user_var(tmp)) {
	gretl_errmsg_sprintf("'%s': name conflicts with different type", tmp);
	err = E_TYPES;
    }
    if (!err) {
	strcpy(chkname, tmp);
    }
    g_free(tmp);

    return err;
}

/**
 * gretl_list_from_matrix:
 * @X: matrix.
 * @prefix: prefix for series names, or NULL.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Tries to add the columns of @X to @dset as series, and
 * if successful constructs a list holding the added series.
 * numbers of series. This is possible only if the rows of
 * @X can be matched to observations in @dset.
 *
 * Returns: allocated list, or NULL on failure.
 */

/* There's actually no possibility of overflow below */
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-overflow"
#endif

int *gretl_list_from_matrix (const gretl_matrix *X,
			     const char *prefix,
			     DATASET *dset,
			     int *err)
{
    int *list = NULL;
    const char **S = NULL;
    int orig_v = dset->v;

    if (gretl_is_null_matrix(X)) {
	list = gretl_null_list();
	if (list == NULL) {
	    *err = E_ALLOC;
	}
	return list;
    }

    S = gretl_matrix_get_colnames(X);
    if (S == NULL && prefix == NULL) {
	*err = E_INVARG;
    } else if (S != NULL) {
	/* check the colnames */
	if (!colnames_ok_for_series(S, X->cols)) {
	    *err = E_INVARG;
	}
    } else {
	/* check the prefix */
	if (strlen(prefix) > VNAMELEN - 3) {
	    *err = E_INVARG;
	}
    }

    if (!*err) {
	char chkname[VNAMELEN];
	int mt1 = gretl_matrix_get_t1(X);
	int mt2 = gretl_matrix_get_t2(X);
	int s1 = 0, s2 = 0;
	int n_add = X->cols;
	int j, slen = 0;
	int ok = 0;

	if (mt2 > 0) {
	    if (X->rows != mt2 - mt1 + 1) {
		; /* time indices not usable! */
	    } else if (mt2 < dset->t1 || mt1 > dset->t2) {
		; /* no overlap of ranges */
	    } else {
		s1 = dset->t1 - mt1;
		s2 = X->rows - s1;
		ok = 1;
	    }
	} else if (X->rows == sample_size(dset)) {
	    /* length matches current sample */
	    s1 = 0;
	    s2 = X->rows - 1;
	    ok = 1;
	} else if (X->rows == dset->n) {
	    /* length matches full series length */
	    s1 = dset->t1;
	    s2 = dset->t2;
	    ok = 1;
	}
	if (!ok) {
	    gretl_errmsg_set("matrix to list: data ranges could not be matched");
	    *err = E_DATA;
	} else {
	    list = gretl_list_new(X->cols);
	    if (list == NULL) {
		*err = E_ALLOC;
	    }
	}

	if (S == NULL) {
	    slen = (int) floor(log10(X->cols)) + 1;
	    if (slen > 6) {
		*err = E_DATA;
	    }
	}

	/* first pass, check the putative series names */
	for (j=0; j<X->cols && !*err; j++) {
	    if (S != NULL) {
		strcpy(chkname, S[j]);
	    } else {
		*err = try_list_vname(chkname, prefix, j+1, slen);
	    }
	    if (!*err && gretl_is_series(chkname, dset)) {
		/* an existing series, decrement the count
		   of series to be added to @dset
		*/
		n_add--;
	    }
	}

	if (!*err && n_add > 0) {
	    *err = dataset_add_NA_series(dset, n_add);
	}
	if (!*err) {
	    int vnew = orig_v;
	    int vj, t, s;

	    for (j=0; j<X->cols && !*err; j++) {
		if (S != NULL) {
		    strcpy(chkname, S[j]);
		} else if (prefix != NULL) {
		    sprintf(chkname, "%s%0*d", prefix, slen, j+1);
		}
		vj = current_series_index(dset, chkname);
		if (vj < 0) {
		    vj = vnew++;
		    strcpy(dset->varname[vj], chkname);
		}
		for (t=dset->t1, s=s1; t<=dset->t2 && s<=s2; t++, s++) {
		    if (s < 0) {
			dset->Z[vj][t] = NADBL;
		    } else {
			dset->Z[vj][t] = gretl_matrix_get(X, s, j);
		    }
		}
		list[j+1] = vj;
	    }
	}
    }

    if (*err) {
	dataset_drop_last_variables(dset, dset->v > orig_v);
	free(list);
	list = NULL;
    }

    return list;
}

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#pragma GCC diagnostic pop
#endif

/**
 * gretl_auxlist_from_vector:
 * @v: source vector.
 * @err: location to receive error code.
 *
 * Unlike gretl_list_from_vector() this function does not
 * require that the elements of @v are valid series IDs. They
 * may represent a set of orders.
 *
 * Returns: a newly allocated gretl list containing the values
 * in @v, or NULL on failure. Note that it is an error if
 * @v is NULL, or is not a vector.
 */

int *gretl_auxlist_from_vector (const gretl_vector *v, int *err)
{
    int i, n = gretl_vector_get_length(v);
    int *list = NULL;

    if (n == 0) {
	*err = E_DATA;
    } else {
	list = gretl_list_new(n);
	if (list == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		list[i+1] = (int) v->val[i];
	    }
	}
    }

    return list;
}

/**
 * varname_match_any:
 * @dset: pointer to dataset information.
 * @pattern: pattern to be matched.
 *
 * Returns: 1 if at least one variable in the dataset has a
 * name that matches @pattern, otherwise 0.
 */

int varname_match_any (const DATASET *dset, const char *pattern)
{
    GPatternSpec *pspec;
    int i, fd, ret = 0;

    fd = gretl_function_depth();

    pspec = g_pattern_spec_new(pattern);

    for (i=1; i<dset->v; i++) {
	if (fd == 0 || fd == series_get_stack_level(dset, i)) {
	    if (g_pattern_match_string(pspec, dset->varname[i])) {
		ret = 1;
		break;
	    }
	}
    }

    g_pattern_spec_free(pspec);

    return ret;
}
