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
#include "libset.h"
#include "gretl_func.h"
#include "uservar.h"
#include "gretl_array.h"
#include "gretl_string_table.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif

#include <errno.h>
#include <glib.h>

enum {
    ST_QUOTED  = 1 << 0,
    ST_ALLINTS = 1 << 1,
    ST_ALLDBLS = 1 << 2
};

struct series_table_ {
    int n_strs;       /* number of strings in table */
    char **strs;      /* saved strings */
    GHashTable *ht;   /* hash table for quick lookup */
    int flags;        /* status flags (above) */
};

struct gretl_string_table_ {
    int *cols_list;       /* list of included columns */
    series_table **cols;  /* per-column tables (see above) */
    char *extra;          /* extra information, if any */
};

#define st_quoted(t) (t->flags & ST_QUOTED)
#define all_ints(t)  (t->flags & ST_ALLINTS)
#define all_dbls(t)  (t->flags & ST_ALLDBLS)
#define all_num(t)   (t->flags & (ST_ALLINTS | ST_ALLDBLS))

static series_table *series_table_alloc (void)
{
    series_table *st = malloc(sizeof *st);

    if (st != NULL) {
	st->strs = NULL;
	st->n_strs = 0;
	st->ht = g_hash_table_new(g_str_hash, g_str_equal);
	st->flags = 0;
    }

    return st;
}

static gretl_string_table *gretl_string_table_alloc (void)
{
    gretl_string_table *gst = malloc(sizeof *gst);

    if (gst != NULL) {
	gst->cols_list = NULL;
	gst->cols = NULL;
	gst->extra = NULL;
    }

    return gst;
}

/**
 * gretl_string_table_new:
 * @list: list of series IDs whose values are to be
 * given a string representation, or NULL.
 *
 * These values in @list should correspond to the 0-based indices
 * of the series in question within the dataset.  For example,
 * if strings are to be recorded for variables 2, 5 and 10 the
 * @list argument would be {3, 2, 5, 10}. If NULL is passed for
 * @list the return value is an initialized, empty string table.
 *
 * Returns: pointer to a newly allocated string table or NULL
 * on failure.
 */

gretl_string_table *gretl_string_table_new (const int *list)
{
    gretl_string_table *gst;
    int ncols = 0;
    int err = 0;

    gst = gretl_string_table_alloc();
    if (gst == NULL) {
	return NULL;
    }

    if (list != NULL && list[0] > 0) {
	gst->cols_list = gretl_list_copy(list);
	if (gst->cols_list == NULL) {
	    err = E_ALLOC;
	} else {
	    ncols = list[0];
	}
    }

    if (ncols > 0) {
	gst->cols = malloc(ncols * sizeof *gst->cols);
	if (gst->cols == NULL) {
	    err = E_ALLOC;
	} else {
	    int i, j;

	    for (i=0; i<ncols && !err; i++) {
		gst->cols[i] = series_table_alloc();
		if (gst->cols[i] == NULL) {
		    err = E_ALLOC;
		    for (j=0; j<i; j++) {
			free(gst->cols[j]);
		    }
		    free(gst->cols);
		}
	    }
	}
    }

    if (err) {
	free(gst->cols_list);
	free(gst);
	gst = NULL;
    }

    return gst;
}

static int series_table_get_index (const series_table *st,
				   const char *s)
{
    int ret = 0;

    if (st != NULL && st->ht != NULL && s != NULL) {
	gpointer p = g_hash_table_lookup(st->ht, s);

	if (p != NULL) {
	    ret = GPOINTER_TO_INT(p);
	}
    }

    return ret;
}

/**
 * series_table_get_value:
 * @st: a gretl series table.
 * @s: the string to look up.
 *
 * Returns: the numerical value associated with @s in the
 * given series table, or #NADBL in case there is no match.
 */

double series_table_get_value (series_table *st, const char *s)
{
    int k = series_table_get_index(st, s);

    return (k > 0)? (double) k : NADBL;
}

/**
 * series_table_get_string:
 * @st: a gretl series table.
 * @val: the numerical value to look up.
 *
 * Returns: the string associated with @val in the
 * given series table, or NULL in case there is no match.
 */

const char *series_table_get_string (series_table *st, double val)
{
    const char *ret = NULL;

    if (!na(val)) {
	int k = (int) lrint(val);

	if (k > 0 && k <= st->n_strs) {
	    ret = st->strs[k-1];
	}
    }

    return ret;
}

/**
 * series_table_map:
 * @st_from: gretl series table.
 * @st_to: gretl series table.
 *
 * Constructs a mapping from the integer codes in @st_from
 * to those in @st_to. For example, if the string "foo"
 * has code 3 in @st_from and code 12 in @st_to, then
 * element 3 in the mapping array will have value 12.
 * For any strings in @st_from that are not matched
 * in @st_to, the associated element of the map is set
 * to -1.
 *
 * Element 0 of the map holds the number of following
 * elements, which is the same as the number of strings in
 * @st_from.
 *
 * Returns: allocated array of int or NULL in case of failure.
 */

int *series_table_map (series_table *st_from, series_table *st_to)
{
    int *map = NULL;
    int n = st_from->n_strs;

    map = gretl_list_new(n);

    if (map != NULL) {
	const char *s1;
	int i, i2;

	for (i=0; i<n; i++) {
	    s1 = st_from->strs[i];
	    i2 = series_table_get_index(st_to, s1);
	    map[i+1] = i2 == 0 ? -1 : i2;
	}
    }

    return map;
}

/**
 * series_table_get_strings:
 * @st: a gretl series table.
 * @n_strs: location to receive the number of strings, or NULL.
 *
 * Returns: the array of strings associated with @st. These
 * should not be modified in any way.
 */

char **series_table_get_strings (series_table *st, int *n_strs)
{
    if (st != NULL) {
	if (n_strs != NULL) {
	    *n_strs = st->n_strs;
	}
	return st->strs;
    } else {
	return NULL;
    }
}

int series_table_get_n_strings (series_table *st)
{
    if (st != NULL) {
	return st->n_strs;
    } else {
	return 0;
    }
}

static char *get_unquoted (const char *s)
{
    char *tmp = NULL;
    int n = strlen(s);

    if (s[n-1] == s[0]) {
	tmp = gretl_strndup(s+1, n-2);
    }

    return tmp;
}

/**
 * series_table_add_string:
 * @st: a gretl series table.
 * @s: new string to add.
 *
 * Returns: the index of the new string within the table, or
 * -1 on failure.
 */

int series_table_add_string (series_table *st, const char *s)
{
    char *tmp = NULL;
    int n, err;

    if (s == NULL) {
	return -1;
    }

    if (*s == '"' || *s == '\'') {
	tmp = get_unquoted(s);
    }

    if (tmp != NULL) {
	st->flags |= ST_QUOTED;
	err = strings_array_add(&st->strs, &st->n_strs, tmp);
	free(tmp);
    } else {
	err = strings_array_add(&st->strs, &st->n_strs, s);
    }

    if (err) {
	n = -1;
    } else {
	n = st->n_strs;
	g_hash_table_insert(st->ht, (gpointer) st->strs[n-1],
			    GINT_TO_POINTER(n));
    }

    return n;
}

/**
 * series_table_add_strings:
 * @st: a gretl series table.
 * @S: array of new strings to add.
 * @ns: number of elements of @S.
 *
 * Appends the strings in @S to @st.

 * Returns: 0 on successful completion, or error code on error.
 */

int series_table_add_strings (series_table *st,
			      const char **S,
			      int ns)
{
    int oldn, newn;
    int i, j, err = 0;

    if (S == NULL || ns <= 0) {
	return 0;
    }

    oldn = st->n_strs;
    newn = oldn + ns;
    st->strs = strings_array_realloc_with_length(&st->strs,
						 oldn, newn,
						 0);
    if (st->strs == NULL) {
	err = E_ALLOC;
    } else {
	st->n_strs = newn;
	for (i=0, j=oldn; i<ns; i++, j++) {
	    st->strs[j] = gretl_strdup(S[i]);
	    g_hash_table_insert(st->ht, (gpointer) st->strs[j],
				GINT_TO_POINTER(j+1));
	}
    }

    return err;
}

series_table *series_table_new (char **strs, int n_strs, int *err)
{
    series_table *st = series_table_alloc();
    int i;

    if (st == NULL || strs == NULL) {
	*err = E_ALLOC;
    } else {
	st->n_strs = n_strs;
	st->strs = strs;
	for (i=0; i<n_strs; i++) {
	    if (st->strs[i] == NULL) {
		fprintf(stderr, "series_table_new: str %d is NULL\n", i);
		*err = E_DATA;
	    } else {
		g_hash_table_insert(st->ht, (gpointer) st->strs[i],
				    GINT_TO_POINTER(i+1));
	    }
	}
    }

    return st;
}

series_table *series_table_copy (series_table *st)
{
    series_table *ret = NULL;

    if (st != NULL) {
	ret = series_table_alloc();
    }

    if (ret != NULL) {
	char **S = strings_array_dup(st->strs, st->n_strs);
	int i;

	if (S == NULL) {
	    series_table_destroy(ret);
	    ret = NULL;
	} else {
	    ret->n_strs = st->n_strs;
	    ret->strs = S;
	    for (i=0; i<ret->n_strs; i++) {
		g_hash_table_insert(ret->ht, (gpointer) ret->strs[i],
				    GINT_TO_POINTER(i+1));
	    }
	}
    }

    return ret;
}

static series_table *
gretl_string_table_add_column (gretl_string_table *gst, int colnum)
{
    series_table **cols;
    int *newlist;
    int n, err = 0;

    newlist = gretl_list_append_term(&gst->cols_list, colnum);
    if (newlist == NULL) {
	return NULL;
    }

    n = gst->cols_list[0];

    cols = realloc(gst->cols, n * sizeof *cols);
    if (cols == NULL) {
	err = E_ALLOC;
    } else {
	gst->cols = cols;
	cols[n-1] = series_table_alloc();
	if (cols[n-1] == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	gst->cols_list[0] -= 1;
	return NULL;
    } else {
	return cols[n-1];
    }
}

/**
 * gretl_string_table_index:
 * @gst: a gretl string table.
 * @s: the string to look up or add.
 * @col: index of the column to be accessed or created.
 * @addcol: non-zero to indicate that a column should be
 * added if it's not already present.
 * @prn: gretl printer (or %NULL).
 *
 * This function has two main uses: for lookup in the context of
 * a completed string table, or for constructing such a
 * table (with @addcol non-zero).  The returned index reflects
 * any additions to the table that may be required (if column
 * @col does not already exist, or if string @s is not already
 * stored for column @col).
 *
 * Returns: the 1-based index of @s within the column of
 * @st that has index @col, if available, otherwise 0.
 */

int
gretl_string_table_index (gretl_string_table *gst, const char *s,
			  int col, int addcol, PRN *prn)
{
    series_table *st = NULL;
    char *tmp = NULL;
    int i, idx = 0;

    if (gst == NULL) {
	return idx;
    }

    if (gst->cols_list != NULL) {
	for (i=1; i<=gst->cols_list[0]; i++) {
	    if (gst->cols_list[i] == col) {
		st = gst->cols[i-1];
		break;
	    }
	}
    }

    if (*s == '"') {
	tmp = get_unquoted(s);
    }

    if (st != NULL) {
	/* there's a table for this column already */
	idx = series_table_get_index(st, tmp != NULL ? tmp : s);
    } else if (addcol) {
	/* no table for this column yet: start one now */
	st = gretl_string_table_add_column(gst, col);
	if (st != NULL) {
	    pprintf(prn, _("variable %d: translating from strings to "
			   "code numbers\n"), col);
	}
    }

    if (idx == 0 && st != NULL) {
	idx = series_table_add_string(st, s);
    }

    free(tmp);

    return idx;
}

/* Used in the context of deletion of "empty" variables from
   an imported dataset: the index of a given "column" in
   a string table is adjusted to match the new position of
   the variable in question.
*/

int gretl_string_table_reset_column_id (gretl_string_table *gst,
					int oldid, int newid)
{
    if (gst != NULL) {
	int i;

	for (i=1; i<=gst->cols_list[0]; i++) {
	    if (gst->cols_list[i] == oldid) {
		gst->cols_list[i] = newid;
		return 0;
	    }
	}
    }

    return E_DATA;
}

series_table *gretl_string_table_detach_col (gretl_string_table *gst,
					     int col)
{
    series_table *st = NULL;

    if (gst != NULL) {
	int pos = in_gretl_list(gst->cols_list, col);
	int i, n = gst->cols_list[0];

	if (pos > 0) {
	    st = gst->cols[pos-1];
	    for (i=pos-1; i<n-1; i++) {
		gst->cols[i] = gst->cols[i+1];
	    }
	    gst->cols[n-1] = NULL;
	    gretl_list_delete_at_pos(gst->cols_list, pos);
	}
    }

    return st;
}

int in_string_table (gretl_string_table *gst, int id)
{
    if (gst != NULL) {
	return in_gretl_list(gst->cols_list, id);
    } else {
	return 0;
    }
}

int *string_table_copy_list (gretl_string_table *gst)
{
    if (gst != NULL) {
	return gretl_list_copy(gst->cols_list);
    } else {
	return NULL;
    }
}

int string_table_replace_list (gretl_string_table *gst,
			       int *newlist)
{
    if (gst != NULL) {
	/* FIXME pruning? */
	free(gst->cols_list);
	gst->cols_list = newlist;
    }

    return E_DATA;
}

/**
 * series_table_destroy:
 * @st: series string table.
 *
 * Frees all resources associated with @st.
 */

void series_table_destroy (series_table *st)
{
    if (st != NULL) {
	strings_array_free(st->strs, st->n_strs);
	if (st->ht != NULL) {
	    g_hash_table_destroy(st->ht);
	}
	free(st);
    }
}

/**
 * series_table_free_shallow:
 * @st: series string table.
 *
 * Does a "shallow"free on @st, without freeing the
 * array of strings to which the table has a reference.
 */

void series_table_free_shallow (series_table *st)
{
    if (st != NULL) {
	if (st->ht != NULL) {
	    g_hash_table_destroy(st->ht);
	}
	free(st);
    }
}

/**
 * gretl_string_table_destroy:
 * @gst: gretl string table.
 *
 * Frees all resources associated with @gst.
 */

void gretl_string_table_destroy (gretl_string_table *gst)
{
    int i, ncols;

    if (gst == NULL) return;

    ncols = (gst->cols_list != NULL)? gst->cols_list[0] : 0;

    for (i=0; i<ncols; i++) {
	series_table_destroy(gst->cols[i]);
    }
    free(gst->cols);

    free(gst->cols_list);

    if (gst->extra != NULL) {
	free(gst->extra);
    }

    free(gst);
}

/* Given a series_table in which all the strings are just
   representations of integers, write the integer values
   into the series and destroy the table, while marking
   the series as "coded". Or, if all the strings are valid
   doubles, convert to numeric.
*/

static void series_commute_string_table (DATASET *dset, int i,
					 series_table *st)
{
    if (dset != NULL && i > 0 && i < dset->v) {
	const char *s;
	double val;
	int t;

	for (t=0; t<dset->n; t++) {
	    val = dset->Z[i][t];
	    if (!na(val)) {
		s = series_table_get_string(st, val);
		if (all_ints(st)) {
		    dset->Z[i][t] = (double) atoi(s);
		} else {
		    dset->Z[i][t] = atof(s);
		}
	    }
	}
	if (all_ints(st)) {
	    series_set_flag(dset, i, VAR_DISCRETE);
	    if (!gretl_isdummy(0, dset->n - 1, dset->Z[i])) {
		series_set_flag(dset, i, VAR_CODED);
	    }
	} else {
	    series_unset_flag(dset, i, VAR_DISCRETE);
	}
	series_table_destroy(st);
    }
}

void series_table_print (DATASET *dset, int i, PRN *prn)
{
    series_table *st = series_get_string_table(dset, i);
    int j;

    if (st != NULL) {
        pprintf(prn, _("String code table for variable %d (%s):\n"),
                i, dset->varname[i]);
	pputc(prn, '\n');
        for (j=0; j<st->n_strs; j++) {
            pprintf(prn, "%5d \"%s\"\n", j+1, st->strs[j]);
        }
        pputc(prn, '\n');
    }
}

/**
 * gretl_string_table_finalize:
 * @gst: gretl string table.
 * @dset: dataset.
 *
 * Attaches the string-value information in @gst to the dataset
 * @dset. However, if one or more of the series referenced by @gst
 * are deemed to be integer codes or misclassified numeric data,
 * their "series tables" are commuted into numeric form.
 *
 * Returns: the number of series tables actually attached.
 */

int gretl_string_table_finalize (gretl_string_table *gst, DATASET *dset)
{
    series_table *st;
    int i, vi, ret = 0;

    if (gst == NULL || gst->cols_list == NULL ||
        dset == NULL || dset->varinfo == NULL) {
	return 0;
    }

    for (i=0; i<gst->cols_list[0]; i++) {
        st = gst->cols[i];
        if (st != NULL) {
            vi = gst->cols_list[i+1];
            if (all_num(st)) {
                series_commute_string_table(dset, vi, st);
            } else {
                series_attach_string_table(dset, vi, st);
                ret++;
            }
            gst->cols[i] = NULL;
        }
    }

    return ret;
}

static int string_is_double (const char *s)
{
    char *test;

    errno = 0;
    strtod(s, &test);
    return errno == 0 && *test == '\0';
}

/**
 * gretl_string_table_validate:
 * @gst: gretl string table.
 * @opt: may include OPT_S to indicate that the string
 * table was constructed from a spreadsheet file, in which
 * case certain cells may be explicitly marked as holding
 * strings, which information may be not be reliable.
 *
 * Checks that the "string values" in @gst are not in fact
 * undigested quasi-numerical values. We run this on
 * imported CSV data to ensure we don't produce
 * misleading results. In the spreadsheet case we're on
 * the lookout for purely numerical data wrongly identified
 * as strings; we don't treat that as fatal but mark the
 * series table(s) for conversion to numeric.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_string_table_validate (gretl_string_table *gst,
				 gretlopt opt)
{
    const char *test = "0123456789.,";
    int ssheet = (opt & OPT_S);
    int i, ncols = 0;
    int err = 0;

    if (gst != NULL && gst->cols_list != NULL) {
	ncols = gst->cols_list[0];
    }

    for (i=0; i<ncols; i++) {
	series_table *st = gst->cols[i];
	const char *s;
	int nint = 0;
	int ndbl = 0;
	int j, myerr = E_DATA;

	if (st_quoted(st)) {
	    myerr = 0;
	}

	for (j=0; j<st->n_strs; j++) {
	    s = st->strs[j];
	    if (st_quoted(st)) {
		if (integer_string(s)) {
		    nint++;
		}
	    } else {
		/* not quoted */
		if (ssheet && (*s == '\0' || string_is_double(s))) {
		    /* could really be numeric? (2020-12-25) */
		    ndbl++;
		    continue;
		} else if (*s == '-' || *s == '+') {
		    s++;
		}
		if (strspn(s, test) < strlen(s)) {
		    /* not quasi-numeric */
		    myerr = 0;
		    break;
		}
	    }
	}

	if (nint == st->n_strs) {
	    /* treat as integer codes */
	    st->flags |= ST_ALLINTS;
	} else if (ndbl == st->n_strs) {
	    /* all really numeric */
	    st->flags |= ST_ALLDBLS;
	    myerr = 0;
	}

	if (myerr) {
	    err = myerr;
	    break;
	}
    }

    return err;
}

/**
 * gretl_string_table_save:
 * @gst: gretl string table.
 * @dset: dataset information (for names of variables).
 *
 * Attaches the content of @gst to @dset.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_string_table_save (gretl_string_table *gst, DATASET *dset)
{
    series_table *st;
    int i, vi, ncols = 0;

    if (gst == NULL || dset->varinfo == NULL) {
	return E_DATA;
    }

    ncols = (gst->cols_list != NULL)? gst->cols_list[0] : 0;

    for (i=0; i<ncols; i++) {
	st = gst->cols[i];
	if (st != NULL) {
	    vi = gst->cols_list[i+1];
	    st = gst->cols[i];
	    series_attach_string_table(dset, vi, st);
	    gst->cols[i] = NULL;
	}
    }

    return 0;
}

/**
 * gretl_string_table_add_extra:
 * @gst: gretl string table.
 * @prn: gretl printer.
 *
 * Steals the printing buffer from @prn and adds it to @gst.
 * The buffer will be appended when @gst is printed via
 * gretl_string_table_finalize().
 */

void gretl_string_table_add_extra (gretl_string_table *gst, PRN *prn)
{
    if (gst != NULL && prn != NULL) {
	gst->extra = gretl_print_steal_buffer(prn);
    }
}

/* apparatus for built-in strings */

struct built_in_string_ {
    char name[VNAMELEN];
    gchar *s;
};

typedef struct built_in_string_ built_in_string;

static built_in_string built_ins[] = {
    { "gretldir", NULL },
    { "dotdir",   NULL },
    { "workdir",  NULL },
    { "gnuplot",  NULL },
    { "x12a",     NULL },
    { "x12adir",  NULL },
    { "tramo",    NULL },
    { "tramodir", NULL },
    { "seats",    NULL },
    { "pkgdir",   NULL },
    { "lang",     NULL },
    { "logfile",  NULL }
};

void builtin_strings_cleanup (void)
{
    int i, n = sizeof built_ins / sizeof built_ins[0];

    for (i=0; i<n; i++) {
	g_free(built_ins[i].s);
    }
}

/**
 * gretl_insert_builtin_string:
 * @name: the name of the string to be added or replaced.
 * @s: the value for this string variable.
 *
 * Inserts value @s for string @name in gretl's table
 * of built-in string variables.
 */

void gretl_insert_builtin_string (const char *name, const char *s)
{
    int i, n = sizeof built_ins / sizeof built_ins[0];
    int m, gui = gretl_in_gui_mode();

    for (i=0; i<n; i++) {
	if (!strcmp(name, built_ins[i].name)) {
	    g_free(built_ins[i].s);
	    if (s == NULL) {
		built_ins[i].s = NULL;
	    } else if (gui && !g_utf8_validate(s, -1, NULL)) {
		/* handle non-ASCII Windows paths */
		gsize bytes;
		gchar *u;

		u = g_locale_to_utf8(s, -1, NULL, &bytes, NULL);
		if (u != NULL) {
		    m = strlen(u);
		    if (u[m-1] == SLASH) {
			u[m-1] = '\0';
		    }
		}
		built_ins[i].s = u;
	    } else {
		m = strlen(s);
		if (s[m-1] == SLASH) {
		    /* drop trailing dir separator for paths */
		    built_ins[i].s = g_strndup(s, m - 1);
		} else {
		    built_ins[i].s = g_strdup(s);
		}
	    }
	    return;
	}
    }
}

char *get_built_in_string_by_name (const char *name)
{
    int i, n = sizeof built_ins / sizeof built_ins[0];

    for (i=0; i<n; i++) {
	if (!strcmp(name, built_ins[i].name)) {
	    char *s = built_ins[i].s;

	    return s != NULL ? s : "";
	}
    }

    return NULL;
}

/* Try to recode the content of a local file or web resource
   to UTF-8. Can be tricky since we don't know the original
   encoding of the content.
*/

static gchar *recode_content (gchar *orig, const char *codeset,
			      int *err)
{
    const gchar *charset = NULL;
    GError *gerr = NULL;
    gsize wrote = 0;
    gchar *tr;

    if (codeset != NULL) {
	/* the user specified the source encoding */
	tr = g_convert(orig, -1, "UTF-8", codeset,
		       NULL, &wrote, &gerr);
    } else if (g_get_charset(&charset)) {
	/* we're in a UTF-8 locale, so we know that
	   g_locale_to_utf8 won't do the job; so guess
	   the content is iso-8859-something?
	*/
	tr = g_convert(orig, -1, "UTF-8", "ISO-8859-15",
		       NULL, &wrote, &gerr);
    } else {
	/* try assuming the material is in the locale
	   encoding */
	tr = g_locale_to_utf8(orig, -1, NULL, &wrote, &gerr);
	if (gerr != NULL) {
	    /* failed: try iso-8859-15? */
	    g_error_free(gerr);
	    gerr = NULL;
	    tr = g_convert(orig, -1, "UTF-8", "ISO-8859-15",
			   NULL, &wrote, &gerr);
	}
    }

    if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	*err = E_DATA;
	g_error_free(gerr);
    }

    g_free(orig);

    return tr;
}

static int shell_grab (const char *arg, char **sout)
{
    int err = 0;

    if (arg == NULL || *arg == '\0') {
	return E_PARSE;
    }

    if (!libset_get_bool(SHELL_OK)) {
	gretl_errmsg_set(_("The shell command is not activated."));
	return 1;
    }

    gretl_shell_grab(arg, sout);

    if (sout != NULL && *sout != NULL) {
	char *content = *sout;

	if (!g_utf8_validate(content, -1, NULL)) {
	    content = recode_content(content, NULL, &err);
	    *sout = content;
	}

	if (content != NULL) {
	    /* trim trailing newline */
	    int n = strlen(content);

	    if (content[n-1] == '\n') {
		content[n-1] = '\0';
	    }
	}
    }

    return err;
}

char *gretl_backtick (const char *arg, int *err)
{
    char *val = NULL;

    *err = shell_grab(arg, &val);

    if (!*err && val == NULL) {
	val = gretl_strdup("");
	if (val == NULL) {
	    *err = E_ALLOC;
	}
    }

    return val;
}

char *gretl_getenv (const char *key, int *defined, int *err)
{
    char *test = getenv(key);
    char *val = NULL;

    if (test == NULL) {
	*defined = 0;
	val = gretl_strdup("");
    } else {
	*defined = 1;
	val = gretl_strdup(test);
    }

    if (val == NULL) {
	*err = E_ALLOC;
    }

    return val;
}

char *retrieve_date_string (int t, const DATASET *dset, int *err)
{
    char *ret = NULL;

    if (t <= 0 || t > dset->n) {
	*err = E_DATA;
    } else if (dset->S != NULL) {
	ret = gretl_strdup(dset->S[t-1]);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	char datestr[OBSLEN] = {0};

	ntolabel(datestr, t - 1, dset);
	ret = gretl_strdup(datestr);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

/* returns a gretl_array of strings on success */

gretl_array *retrieve_date_strings (const gretl_vector *v,
				    const DATASET *dset,
				    int *err)
{
    gretl_array *ret = NULL;
    char *s = NULL;
    int i, t, n;

    n = gretl_vector_get_length(v);
    if (n == 0) {
	*err = E_INVARG;
    } else {
	ret = gretl_array_new(GRETL_TYPE_STRINGS, n, err);
    }

    for (i=0; i<n && !*err; i++) {
	t = gretl_int_from_double(v->val[i], err);
	if (!*err) {
	    s = retrieve_date_string(t, dset, err);
	}
	if (!*err) {
	    gretl_array_set_data(ret, i, s);
	}
    }

    if (*err && ret != NULL) {
	gretl_array_destroy(ret);
	ret = NULL;
    }

    return ret;
}

static int is_web_resource (const char *s)
{
    if (!strncmp(s, "http://", 7) ||
	!strncmp(s, "https://", 8) ||
	!strncmp(s, "ftp://", 6)) {
	return 1;
    } else {
	return 0;
    }
}

static gchar *gzipped_file_get_content (const char *fname,
					int *err)
{
    gzFile fz = gretl_gzopen(fname, "rb");
    gchar *ret = NULL;

    if (fz == NULL) {
	*err = E_FOPEN;
    } else {
	size_t len = 0;
	int chk;

	while (gzgetc(fz) > 0) {
	    len++;
	}
	if (len > 0) {
	    gzrewind(fz);
	    ret = g_try_malloc(len + 1);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    } else {
		chk = gzread(fz, ret, len);
		if (chk <= 0) {
		    *err = E_DATA;
		}
		ret[len] = '\0';
	    }
	} else {
	    ret = g_strdup("");
	}
	gzclose(fz);
    }

    return ret;
}

static gchar *regular_file_get_content (const char *fname,
					int *err)
{
    GError *gerr = NULL;
    gchar *ret = NULL;
    size_t len = 0;
    int done = 0;

#ifdef WIN32
    /* g_file_get_contents() requires a UTF-8 filename */
    if (!g_utf8_validate(fname, -1, NULL)) {
	gchar *fconv;
	gsize wrote = 0;

	fconv = g_locale_to_utf8(fname, -1, NULL, &wrote, &gerr);
	if (fconv != NULL) {
	    g_file_get_contents(fconv, &ret, &len, &gerr);
	    g_free(fconv);
	}
	done = 1;
    }
#endif
    if (!done) {
	g_file_get_contents(fname, &ret, &len, &gerr);
    }

    if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	*err = E_FOPEN;
	g_error_free(gerr);
    }

    return ret;
}

char *retrieve_file_content (const char *fname, const char *codeset,
			     int *err)
{
    char *ret = NULL;
    gchar *content = NULL;
    size_t len = 0;
    gssize sz;

    if (fname == NULL || *fname == '\0') {
	*err = E_INVARG;
    } else {
        len = strlen(fname);
        if (len >= MAXLEN) {
            gretl_errmsg_sprintf(_("filename too long (%d bytes)"), (int) len);
            *err = E_INVARG;
        }
        len = 0;
    }

    if (*err) {
	return NULL;
    }

    if (is_web_resource(fname)) {
#ifdef USE_CURL
	content = retrieve_public_file_as_buffer(fname, &len, err);
#else
	gretl_errmsg_set(_("Internet access not supported"));
	*err = E_DATA;
#endif
    } else {
	char fullname[MAXLEN] = {0};

	strcpy(fullname, fname);
	gretl_addpath(fullname, 0);
	if (is_gzipped(fullname)) {
	    content = gzipped_file_get_content(fullname, err);
	} else {
	    content = regular_file_get_content(fullname, err);
	}
    }

    sz = (len > 0)? len : -1;
    if (content != NULL && !g_utf8_validate(content, sz, NULL)) {
	content = recode_content(content, codeset, err);
    }

    if (content != NULL) {
	if (*err == 0) {
	    ret = gretl_strdup(content);
	}
	g_free(content);
    }

    return ret;
}

/* inserting string into format portion of (s)printf command:
   double any backslashes to avoid breakage of Windows paths
*/

static char *mod_strdup (const char *s)
{
    char *ret = NULL;
    int i, n = strlen(s);
    int bs = 0;

    for (i=0; i<n; i++) {
	if (s[i] == '\\' && (i == n - 1 || s[i+1] != '\\')) {
	    bs++;
	}
    }

    ret = malloc(n + 1 + bs);
    if (ret == NULL) {
	return NULL;
    }

    if (bs == 0) {
	strcpy(ret, s);
    } else {
	int j = 0;

	for (i=0; i<n; i++) {
	    if (s[i] == '\\' && (i == n - 1 || s[i+1] != '\\')) {
		ret[j++] = '\\';
	    }
	    ret[j++] = s[i];
	}
	ret[j] = '\0';
    }

    return ret;
}

static char *maybe_get_subst (char *name, int *n, int quoted,
			      int *freeit)
{
    char *s = NULL;
    int k = *n - 1;
    char *ret = NULL;

    while (k >= 0) {
	s = (char *) get_string_by_name(name);
	if (s != NULL) {
	    *n = k + 1;
	    ret = s;
	    break;
	}
	name[k--] = '\0';
    }

    if (ret != NULL) {
	if (quoted && strchr(ret, '\\')) {
	    ret = mod_strdup(ret);
	    *freeit = 1;
	}
    }

    return ret;
}

static void too_long (void)
{
    gretl_errmsg_sprintf(_("Maximum length of command line "
			   "(%d bytes) exceeded\n"), MAXLINE);
}

int substitute_named_strings (char *line, int *subst)
{
    char sname[VNAMELEN];
    int len = strlen(line);
    char *sub, *tmp, *s = line;
    int bs = 0, in_format = 0;
    int freeit;
    int n, m, err = 0;

    *subst = 0;

    if (*s == '#' || strchr(s, '@') == NULL) {
	return 0;
    }

    if (!strncmp(line, "printf", 6) || !strncmp(line, "sprintf", 7)) {
	s = strchr(s, '"');
	if (s == NULL) {
	    /* no format string */
	    return E_PARSE;
	}
	s++;
	in_format = 1;
    }

    while (*s && !err) {
	if (in_format) {
	    if (*s == '"' && (bs % 2 == 0)) {
		/* reached end of (s)printf format string */
		in_format = 0;
	    }
	    if (*s == '\\') {
		bs++;
	    } else {
		bs = 0;
	    }
	}
	if (*s == '@') {
	    n = gretl_namechar_spn(s + 1);
	    if (n > 0) {
		if (n >= VNAMELEN) {
		    n = VNAMELEN - 1;
		}
		*sname = '\0';
		strncat(sname, s + 1, n);
		freeit = 0;
		sub = maybe_get_subst(sname, &n, in_format, &freeit);
		if (sub != NULL) {
		    m = strlen(sub);
		    if (len + m + 2 >= MAXLINE) {
			too_long();
			err = 1;
			break;
		    }
		    tmp = gretl_strdup(s + n + 1);
		    if (tmp == NULL) {
			err = E_ALLOC;
		    } else {
			strcpy(s, sub);
			strcpy(s + m, tmp);
			free(tmp);
			len += m - (n + 1);
			s += m - 1;
			*subst = 1;
		    }
		    if (freeit) {
			free(sub);
		    }
		}
	    }
	}
	s++;
    }

    return err;
}
