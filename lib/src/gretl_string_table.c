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
#include "gretl_string_table.h"
#ifdef USE_CURL
# include "gretl_www.h"
#endif

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

    if (gst->cols_list != NULL) {
	for (i=1; i<=gst->cols_list[0]; i++) {
	    if (gst->cols_list[i] == col) {
		st = gst->cols[i-1];
		break;
	    }
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
		"These variables have been given numeric codes as follows.\n\n"), fp);
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
	if (dset->varinfo != NULL) {
	    series_attach_string_table(dset, vi, st);
	    gst->cols[i] = NULL;
	}
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
	vi = gst->cols_list[i+1];
	st = gst->cols[i];
	series_attach_string_table(dset, vi, st);
	gst->cols[i] = NULL;
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
 * gretl_string_table_print().
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
    char *s;
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
    { "shelldir", NULL },
    { "Rbin",     NULL },
    { "Rlib",     NULL },
    { "pkgdir",   NULL }
};

#ifdef WIN32
static built_in_string dsep = { "dirsep", "\\" };
#else
static built_in_string dsep = { "dirsep", "/" };
#endif

void builtin_strings_cleanup (void)
{
    int i, n = sizeof built_ins / sizeof built_ins[0];

    for (i=0; i<n; i++) {
	free(built_ins[i].s);
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

char *get_built_in_string_by_name (const char *name)
{
    if (!strcmp(name, "dirsep")) {
	return dsep.s;
    } else {
	int i, n = sizeof built_ins / sizeof built_ins[0];

	for (i=0; i<n; i++) {
	    if (!strcmp(name, built_ins[i].name)) {
		return built_ins[i].s;
	    }
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

	ntodate(datestr, t - 1, dset);
	ret = gretl_strdup(datestr);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
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

char *retrieve_file_content (const char *fname, const char *codeset,
			     int *err)
{
    char *content = NULL;
    size_t len = 0;

    if (fname == NULL || *fname == '\0') {
	*err = E_DATA;
    } else if (is_web_resource(fname)) {
#ifdef USE_CURL
	content = retrieve_public_file_as_buffer(fname, &len, err);
	if (!*err && !g_utf8_validate(content, len, NULL)) {
	    content = recode_content(content, codeset, err);
	}
#else
	gretl_errmsg_set(_("Internet access not supported"));
	*err = E_DATA;
#endif
    } else {
	char fullname[FILENAME_MAX];
	GError *gerr = NULL;

	*fullname = '\0';
	strncat(fullname, fname, FILENAME_MAX - 1);
	gretl_addpath(fullname, 0);

	g_file_get_contents(fullname, &content, &len, &gerr);

	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    *err = E_FOPEN;
	    g_error_free(gerr);
	} else if (!g_utf8_validate(content, len, NULL)) {
	    content = recode_content(content, codeset, err);
	}
    }

    return content;
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

#define PRINTF_SPECIAL 0

int substitute_named_strings (char *line, int *subst)
{
    char sname[VNAMELEN];
    int len = strlen(line);
    char *sub, *tmp, *s = line;
    int bs = 0, in_format = 0;
    int freeit;
    int i, n, m, err = 0;

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

    i = s - line;

    while (*s && !err) {
	if (in_format) {
	    if (*s == '"' && (bs % 2 == 0)) {
		/* reached end of (s)printf format string */
		in_format = 0;
#if PRINTF_SPECIAL
		break;
#endif
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
