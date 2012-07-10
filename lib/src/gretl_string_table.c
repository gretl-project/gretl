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
#include "gretl_string_table.h"

#include <glib.h>

struct _series_table {
    int n_strs;       /* number of strings in table */
    char **strs;      /* saved strings */
    GHashTable *ht;   /* hash table for quick lookup */
};

struct _gretl_string_table {
    int *cols_list;       /* list of included columns */
    series_table **cols;  /* per-column tables (see above) */
    char *extra;          /* extra information, if any */
};

static series_table *series_table_alloc (void)
{
    series_table *st = malloc(sizeof *st);

    if (st != NULL) {
	st->strs = NULL;
	st->n_strs = 0;
	st->ht = g_hash_table_new(g_str_hash, g_str_equal);
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
 * given a string representation. or NULL. 
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
    gpointer p = g_hash_table_lookup(st->ht, s);
    int ret = 0;

    if (p != NULL) {
	ret = GPOINTER_TO_INT(p);
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
	int k = val;

	if (k > 0 && k <= st->n_strs) {
	    ret = st->strs[k-1];
	}
    }

    return ret;
}

/**
 * series_table_get_strings:
 * @st: a gretl series table.
 * @n_strs: location to receive the number of strings.
 *
 * Returns: the array of strings associated with @st. These
 * should not be modified in any way.
 */

const char **series_table_get_strings (series_table *st, int *n_strs)
{
    *n_strs = st->n_strs;
    return (const char **) st->strs;
}

static int series_table_add_string (series_table *st, 
				    const char *s)
{
    int n, err;

    err = strings_array_add(&st->strs, &st->n_strs, s);

    if (err) {
	n = -1;
    } else {
	n = st->n_strs;
	g_hash_table_insert(st->ht, (gpointer) st->strs[n-1], 
			    GINT_TO_POINTER(n));
    }

    return n;
}

series_table *series_table_new (char **strs, int n_strs)
{
    series_table *st = series_table_alloc();
    int i;

    if (st != NULL) {
	st->n_strs = n_strs;
	st->strs = strs;
	for (i=0; i<n_strs; i++) {
	    g_hash_table_insert(st->ht, (gpointer) st->strs[i], 
				GINT_TO_POINTER(i+1));
	}
    }

    return st;
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
    int i, idx = 0;

    if (gst == NULL) return idx;

    for (i=1; i<=gst->cols_list[0]; i++) {
	if (gst->cols_list[i] == col) {
	    st = gst->cols[i-1];
	    break;
	}
    }

    if (st != NULL) {
	/* there's a table for this column already */
	idx = series_table_get_index(st, s);
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

#define SAVE_STRING_TABLES 1

/**
 * gretl_string_table_print:
 * @gst: gretl string table.
 * @dset: dataset information (for names of variables).
 * @fname: name of the datafile to which the table pertains.
 * @prn: gretl printer (or %NULL).
 *
 * Prints table @gst to a file named string_table.txt in the
 * user's working directory.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_string_table_print (gretl_string_table *gst, DATASET *dset,
			      const char *fname, PRN *prn)
{
    series_table *st;
    const char *fshort;
    char stname[MAXLEN];
    FILE *fp;
    int i, j, ncols = 0;
    int err = 0;

    if (gst == NULL) {
	return E_DATA;
    }

    strcpy(stname, "string_table.txt");
    gretl_path_prepend(stname, gretl_workdir());

    fp = gretl_fopen(stname, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fshort = strrchr(fname, SLASH);
    if (fshort != NULL) {
	fprintf(fp, "%s\n", fshort + 1);
    } else {
	fprintf(fp, "%s\n", fname);
    }

    ncols = (gst->cols_list != NULL)? gst->cols_list[0] : 0;

    if (ncols > 0) {
	fputc('\n', fp);
	fputs(_("One or more non-numeric variables were found.\n"
		"Gretl cannot handle such variables directly, so they\n"
		"have been given numeric codes as follows.\n\n"), fp);
	if (gst->extra != NULL) {
	    fputs(_("In addition, some mappings from numerical values to string\n"
		    "labels were found, and are printed below.\n\n"), fp);
	}
    }
    
    for (i=0; i<ncols; i++) {
	int vi = gst->cols_list[i+1];

	st = gst->cols[i];
	if (i > 0) {
	    fputc('\n', fp);
	}
	fprintf(fp, _("String code table for variable %d (%s):\n"), 
		vi, dset->varname[vi]);
	for (j=0; j<st->n_strs; j++) {
	    fprintf(fp, "%3d = '%s'\n", j+1, st->strs[j]);
	}
#if SAVE_STRING_TABLES
	if (dset->varinfo != NULL) {
	    series_attach_string_table(dset, vi, st);
	    gst->cols[i] = NULL;
	}
#endif
    }

    if (gst->extra != NULL) {
	fputs(gst->extra, fp);
    }

    pprintf(prn, _("String code table written to\n %s\n"), stname);

    fclose(fp);
    set_string_table_written();

    return err;
}

/**
 * gretl_string_table_add_extra:
 * @gst: gretl string table.
 * @prn: gretl printer.
 *
 * Steals the printing buffer from @prn and adds it to @gst.
 * The buffer will be appended when @gst is printed via
 * gretl_string_table_print().
 */

void gretl_string_table_add_extra (gretl_string_table *gst, PRN *prn)
{
    if (gst != NULL && prn != NULL) {
	gst->extra = gretl_print_steal_buffer(prn);
    }
}

/* below: saving of user-defined strings */

typedef struct saved_string_ saved_string;

struct saved_string_ {
    char name[VNAMELEN];
    int level;
    char *s;
};

static int n_saved_strings;
static saved_string *saved_strings;

static saved_string built_ins[] = {
    { "gretldir", 0, NULL },
    { "dotdir",   0, NULL },
    { "workdir",  0, NULL },
    { "gnuplot",  0, NULL },
    { "x12a",     0, NULL },
    { "x12adir",  0, NULL },
    { "tramo",    0, NULL },
    { "tramodir", 0, NULL },
    { "seats",    0, NULL },
    { "shelldir", 0, NULL },
    { "Rbin",     0, NULL },
    { "Rlib",     0, NULL },
    { "pkgdir",   0, NULL }
};

#ifdef WIN32
static saved_string dsep = { "dirsep", 0, "\\" };
#else
static saved_string dsep = { "dirsep", 0, "/" };
#endif

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
    int n = sizeof built_ins / sizeof built_ins[0];
    int i, m;

    for (i=0; i<n; i++) {
	if (!strcmp(name, built_ins[i].name)) {
	    free(built_ins[i].s);
	    if (s == NULL) {
		built_ins[i].s = NULL;
	    } else {
		m = strlen(s);
		if (s[m-1] == SLASH) {
		    /* drop trailing dir separator for paths */
		    built_ins[i].s = gretl_strndup(s, m - 1);
		} else {
		    built_ins[i].s = gretl_strdup(s);
		}
	    }
	    return;
	}
    }
}

static void gretl_free_builtin_strings (void)
{
    int i, n = sizeof built_ins / sizeof built_ins[0];

    for (i=0; i<n; i++) {
	free(built_ins[i].s);
    }    
}

static saved_string *get_saved_string_by_name (const char *name,
					       int *builtin)
{
    int i, d;

    if (builtin != NULL && !strcmp(name, "dirsep")) {
	*builtin = 1;
	return &dsep;
    }

    if (builtin != NULL) {
	int n = sizeof built_ins / sizeof built_ins[0];

	for (i=0; i<n; i++) {
	    if (!strcmp(name, built_ins[i].name)) {
		*builtin = 1;
		return &built_ins[i];
	    }
	}
    }

    d = gretl_function_depth();

    for (i=0; i<n_saved_strings; i++) {
	if (saved_strings[i].level == d &&
	    !strcmp(name, saved_strings[i].name)) {
	    return &saved_strings[i];
	}
    }

    return NULL;
}

/**
 * get_string_by_name:
 * @name: the name of the string variable to access.
 *
 * Returns: the value of string variable @name, or %NULL
 * if there is no such variable.
 */

const char *get_string_by_name (const char *name)
{
    int i, n, d;

    if (!strcmp(name, "dirsep")) {
	return dsep.s;
    } else if (!strcmp(name, "vcvtype")) {
	return last_model_get_vcv_type();
    }

    n = sizeof built_ins / sizeof built_ins[0];

    for (i=0; i<n; i++) {
	if (!strcmp(name, built_ins[i].name)) {
	    return built_ins[i].s;
	}
    }

    d = gretl_function_depth();

    for (i=0; i<n_saved_strings; i++) {
	if (saved_strings[i].level == d &&
	    !strcmp(name, saved_strings[i].name)) {
	    return saved_strings[i].s;
	}
    }

    return NULL;
}

static saved_string *add_named_string (const char *name)
{
    int n = n_saved_strings;
    saved_string *S;

    S = realloc(saved_strings, (n + 1) * sizeof *S);
    if (S == NULL) {
	return NULL;
    }

    strcpy(S[n].name, name);
    S[n].level = gretl_function_depth();
    S[n].s = NULL;
    saved_strings = S;
    n_saved_strings += 1;

    return &S[n];
}

/**
 * add_string_as:
 * @s: string value to be added.
 * @name: the name of the string variable to add.
 *
 * Adds @s to the saved array of string variables 
 * under the name @name.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int add_string_as (const char *s, const char *name)
{
    int n = n_saved_strings;
    saved_string *S;
    int err = 0;

    if (name == NULL || s == NULL) {
	return E_DATA;
    }

    S = realloc(saved_strings, (n + 1) * sizeof *S);
    if (S == NULL) {
	return E_ALLOC;
    }

    saved_strings = S;

    S[n].s = gretl_strdup(s);
    if (S[n].s == NULL) {
	err = E_ALLOC;
    } else {  
	strcpy(S[n].name, name);
	S[n].level = gretl_function_depth() + 1;
	n_saved_strings += 1;
    }

    return err;
}

static int destroy_saved_string (saved_string *S)
{
    int i, j, ns = n_saved_strings - 1;
    int err = 0;

    for (i=0; i<n_saved_strings; i++) {
	if (&saved_strings[i] == S) {
	    free(saved_strings[i].s);
	    for (j=i; j<ns; j++) {
		saved_strings[j] = saved_strings[j+1];
	    }
	    break;
	}
    } 

    if (ns == 0) {
	free(saved_strings);
	saved_strings = NULL;
    } else {
	saved_string *tmp = realloc(saved_strings, ns * sizeof *tmp);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    saved_strings = tmp;
	}
    }

    n_saved_strings = ns;

    return err;
}

/**
 * delete_saved_string:
 * @name: the name of the string variable to delete.
 * @prn: gretl printer (or %NULL).
 *
 * Deletes the string variable called @name.
 *
 * Returns: 0 on success, non-zero on failure (e.g.
 * attempting to delete a built-in string variable).
 */

int delete_saved_string (const char *name, PRN *prn)
{
    saved_string *S;
    int err, builtin = 0;

    S = get_saved_string_by_name(name, &builtin);
    if (S == NULL) {
	err = E_UNKVAR;
    } else if (builtin) {
	gretl_errmsg_sprintf(_("You cannot delete '%s'"), name);
	err = E_DATA;
    } else {
	err = destroy_saved_string(S);
	if (!err && prn != NULL && gretl_messages_on()) {
	    pprintf(prn, _("Deleted string %s"), name);
	    pputc(prn, '\n');
	}
    }

    return err;
}

void saved_strings_cleanup (void)
{
    int i;

    for (i=0; i<n_saved_strings; i++) {
	free(saved_strings[i].s);
    }

    free(saved_strings);
    saved_strings = NULL;
    n_saved_strings = 0;

    gretl_free_builtin_strings();
}

void destroy_user_strings (void)
{
    int i;

    for (i=0; i<n_saved_strings; i++) {
	free(saved_strings[i].s);
    }

    free(saved_strings);
    saved_strings = NULL;
    n_saved_strings = 0;
}

/* called on exiting a user-defined function to clean
   up any strings defined therein */

int destroy_saved_strings_at_level (int d)
{
    int i, ndel = 0;
    int err = 0;

    for (i=0; i<n_saved_strings; i++) {
	if (saved_strings[i].level == d) {
	    ndel++;
	}
    }

    if (ndel > 0) {
	int nnew = n_saved_strings - ndel;

	if (nnew == 0) {
	    for (i=0; i<n_saved_strings; i++) {
		free(saved_strings[i].s);
	    }
	    free(saved_strings);
	    saved_strings = NULL;
	    n_saved_strings = 0;
	} else {	    
	    saved_string *S = malloc(nnew * sizeof *S);
	    int j = 0;

	    if (S == NULL) {
		err = E_ALLOC;
	    } else {
		for (i=0; i<n_saved_strings; i++) {
		    if (saved_strings[i].level == d) {
			free(saved_strings[i].s);
		    } else {
			strcpy(S[j].name, saved_strings[i].name);
			S[j].level = saved_strings[i].level;
			S[j].s = saved_strings[i].s;
			j++;
		    }
		}
		free(saved_strings);
		saved_strings = S;
		n_saved_strings = nnew;
	    }
	}
    }

    return err;
}

#if 0
static void copy_unescape (char *targ, const char *src, int n)
{
    int c, i = 0;

    for (i=0; i<n; i++) {
	if (src[i] == '\\') {
	    c = src[i+1];
	    if (c == '\\' || c == '"') {
		*targ++ = c;
		i++;
	    } else {
		*targ++ = src[i];
	    }
	} else {
	    *targ++ = src[i];
	}
    }	

    *targ = '\0';
}
#endif

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
	/* trim trailing newline */
	int n = strlen(*sout);

	if ((*sout)[n-1] == '\n') {
	    (*sout)[n-1] = '\0';
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

#if 0
static char *strim (char *s)
{
    int n;

    tailstrip(s);
    n = strlen(s);
    if (n > 0 && s[n-1] == '"') {
	s[n-1] = '\0';
    }

    return s;
}
#endif

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

	ntodate(datestr, t - 1, dset);
	ret = gretl_strdup(datestr);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    } 

    return ret;
}

char *retrieve_file_content (const char *fname, int *err)
{
    char *ret = NULL;

    if (fname == NULL || *fname == '\0') {
	*err = E_DATA;
    } else {
	char fullname[FILENAME_MAX];
	GError *gerr = NULL;
	gsize len = 0;

	*fullname = '\0';
	strncat(fullname, fname, FILENAME_MAX - 1);
	gretl_addpath(fullname, 0);

	g_file_get_contents(fullname, &ret, &len, &gerr);

	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	    *err = E_FOPEN;
	}
    } 

    return ret;
}

/**
 * gretl_is_string:
 * @sname: string to test.
 *
 * Returns: 1 if @sname os the name of a currently defined
 * string variable, otherwise 0.
 */

int gretl_is_string (const char *sname)
{
    saved_string *str;
    int builtin = 0;

    if (*sname == '@' && *(sname + 1) != '@') {
	sname++;
    }
    
    str = get_saved_string_by_name(sname, &builtin);

    return (str != NULL && str->s != NULL);
}

int is_user_string (const char *sname)
{
    int i, d;

    if (*sname == '@' && *(sname + 1) != '@') {
	sname++;
    }

    d = gretl_function_depth();

    for (i=0; i<n_saved_strings; i++) {
	if (saved_strings[i].level == d &&
	    !strcmp(sname, saved_strings[i].name)) {
	    return 1;
	}
    }

    return 0;
}

int save_named_string (const char *name, const char *s, PRN *prn)
{
    saved_string *str;
    int builtin = 0;
    int add = 0;
    int err = 0;

    str = get_saved_string_by_name(name, &builtin);

    if (str != NULL && builtin) {
	if (prn == NULL) {
	    gretl_errmsg_sprintf(_("You cannot overwrite '%s'\n"), name);
	} else {
	    pprintf(prn, _("You cannot overwrite '%s'\n"), name);
	}
	return E_DATA;
    }	

    if (str != NULL) {
	free(str->s);
	str->s = NULL;
    } else {
	str = add_named_string(name);
	if (str == NULL) {
	    return E_ALLOC;
	}
	add = 1;
    } 

    str->s = gretl_strdup(s);
    if (str->s == NULL) {
	err = E_ALLOC;
    }

    if (prn != NULL && !err && gretl_messages_on() && 
	!gretl_looping_quietly() && *s != '\0') {
	if (add) {
	    pprintf(prn, _("Generated string %s\n"), name); 
	} else {
	    pprintf(prn, _("Replaced string %s\n"), name);
	}
    }

    return err;
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
    saved_string *str;
    int builtin = 0;
    int k = *n - 1;
    char *ret = NULL;

    while (k >= 0) {
	str = get_saved_string_by_name(name, &builtin);
	if (str != NULL && str->s != NULL) {
	    *n = k + 1;
	    ret = str->s;
	    break;
	}
	name[k--] = '\0';
    }

    if (ret != NULL) {
	if (quoted) {
	    if (strchr(ret, '\\')) {
		ret = mod_strdup(ret);
		*freeit = 1;
	    }
	} 
    }

    return ret;
}

static void too_long (void)
{
    gretl_errmsg_sprintf(_("Maximum length of command line "
			   "(%d bytes) exceeded\n"), MAXLINE);
}

#define var_context(s,i) (i > 8 && !strncmp(s - 9, "isstring(", 9))

int substitute_named_strings (char *line, int *subst)
{
    char sname[VNAMELEN];
    int len = strlen(line);
    char *sub, *tmp, *s = line;
    int bs = 0, pf = 0, quoted = 0;
    int freeit;
    int i, n, m, err = 0;

    if (*s == '#' || strchr(s, '@') == NULL) {
	return 0;
    }

    if (!strncmp(line, "sscanf", 6)) {
	/* when scanning, let @foo be handled as a variable (FIXME?) */
	return 0;
    }    

    if (!strncmp(line, "printf", 6) || !strncmp(line, "sprintf", 7)) {
	pf = 1;
	s = strchr(s, '"');
	if (s == NULL) {
	    /* no format string */
	    return E_PARSE;
	}
	s++;
	quoted = 1;
    }

    i = s - line;

    while (*s && !err) {
	if (pf) {
	    if (*s == '"' && (bs % 2 == 0)) {
		/* reached end of format string: stop substituting */
		break;
	    }
	    if (*s == '\\') {
		bs++;
	    } else {
		bs = 0;
	    }
	}
	if (*s == '@' && !var_context(s, i)) {
	    n = gretl_namechar_spn(s + 1);
	    if (n > 0) {
		if (n >= VNAMELEN) {
		    n = VNAMELEN - 1;
		}
		*sname = '\0';
		strncat(sname, s + 1, n);
		freeit = 0;
		sub = maybe_get_subst(sname, &n, quoted, &freeit);
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
			i += m - 1;
			*subst = 1;
		    }
		    if (freeit) {
			free(sub);
		    }
		}
	    }
	}
	s++;
	i++;
    }

    return err;
}
