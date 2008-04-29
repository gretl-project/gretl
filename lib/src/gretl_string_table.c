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

typedef struct _col_table col_table;

struct _col_table {
    int idx;
    int n_strs;
    char **strs;
    GHashTable *hash;
};

struct _gretl_string_table {
    int n_cols;
    col_table **cols;
};

static void clean_up_codevars (void);

static col_table *col_table_new (int colnum)
{
    col_table *ct = malloc(sizeof *ct);

    if (ct != NULL) {
	ct->strs = NULL;
	ct->n_strs = 0;
	ct->idx = colnum;
	ct->hash = g_hash_table_new(g_str_hash, g_str_equal);
    }

    return ct;
}

gretl_string_table *gretl_string_table_new (int *err)
{
    gretl_string_table *st = malloc(sizeof *st);

    if (st == NULL) {
	*err = E_ALLOC;
    } else {
	st->cols = NULL;
	st->n_cols = 0;
    }

    return st;
}

gretl_string_table *string_table_new_from_cols_list (int *list)
{
    gretl_string_table *st;
    int ncols = list[0];
    int i, j;

    st = malloc(sizeof *st);
    if (st == NULL) return NULL;

    st->cols = malloc(ncols * sizeof *st->cols);
    if (st->cols == NULL) {
	free(st);
	st = NULL;
    } else {
	st->n_cols = ncols;
	for (i=0; i<ncols; i++) {
	    st->cols[i] = col_table_new(list[i+1]);
	    if (st->cols[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(st->cols[j]);
		}
		free(st->cols);
		free(st);
		st = NULL;
	    } 
	}
    }

    return st;
}

static int col_table_get_index (const col_table *ct, const char *s)
{
    gpointer p = g_hash_table_lookup(ct->hash, s);
    int ret = 0;

    if (p != NULL) {
	ret = GPOINTER_TO_INT(p);
    }

    return ret;
}

static int 
col_table_add_string (col_table *ct, const char *s)
{
    char **strs;
    int n = ct->n_strs + 1;
    int ret = n;

    strs = realloc(ct->strs, n * sizeof *strs);
    if (strs == NULL) {
	ret = -1;
    } else {
	ct->strs = strs;
	strs[n-1] = gretl_strdup(s);

	if (strs[n-1] == NULL) {
	    ret = -1;
	} else {
	    ct->n_strs += 1;
	    g_hash_table_insert(ct->hash, (gpointer) strs[n-1], 
				GINT_TO_POINTER(n));
	}
    }

    return ret;
}

static col_table *
gretl_string_table_add_column (gretl_string_table *st, int colnum)
{
    col_table **cols;
    int n = st->n_cols + 1;

    cols = realloc(st->cols, n * sizeof *cols);
    if (cols == NULL) return NULL;

    st->cols = cols;
    cols[n-1] = col_table_new(colnum);
    if (cols[n-1] == NULL) return NULL;

    st->n_cols += 1;

    return cols[n-1];
}

int 
gretl_string_table_index (gretl_string_table *st, const char *s, int col,
			  int addcol, PRN *prn)
{
    col_table *ct = NULL;
    int i, idx = 0;

    if (st == NULL) return idx;

    for (i=0; i<st->n_cols; i++) {
	if ((st->cols[i])->idx == col) {
	    ct = st->cols[i];
	    break;
	}
    }

    if (ct != NULL) {
	/* there's a table for this column already */
	idx = col_table_get_index(ct, s);
    } else if (addcol) {
	/* no table for this column yet: start one now */
	ct = gretl_string_table_add_column(st, col);
	if (ct != NULL) {
	    pprintf(prn, M_("variable %d: translating from strings to "
			    "code numbers\n"), 
		    col);
	}
    }

    if (idx == 0 && ct != NULL) {
	idx = col_table_add_string(ct, s);
    }

    return idx;
}

static void col_table_destroy (col_table *ct)
{
    int i;

    if (ct == NULL) return;

    for (i=0; i<ct->n_strs; i++) {
	free(ct->strs[i]);
    }
    free(ct->strs);

    if (ct->hash != NULL) {
	g_hash_table_destroy(ct->hash);
    }

    free(ct);
}

void gretl_string_table_destroy (gretl_string_table *st)
{
    int i;

    if (st == NULL) return;

    for (i=0; i<st->n_cols; i++) {
	col_table_destroy(st->cols[i]);
    }
    free(st->cols);

    free(st);
}

int gretl_string_table_print (gretl_string_table *st, DATAINFO *pdinfo,
			      const char *fname, PRN *prn)
{
    int i, j;
    const col_table *ct;
    const char *fshort;
    char stname[MAXLEN];
    FILE *fp;
    int err = 0;

    if (st == NULL) return E_DATA;

    strcpy(stname, "string_table.txt");
    gretl_path_prepend(stname, gretl_work_dir());

    fp = gretl_fopen(stname, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fshort = strrchr(fname, SLASH);
    if (fshort != NULL) {
	fprintf(fp, "%s\n\n", fshort + 1);
    } else {
	fprintf(fp, "%s\n\n", fname);
    }

    fputs(M_("One or more non-numeric variables were found.\n"
	     "Gretl cannot handle such variables directly, so they\n"
	     "have been given numeric codes as follows.\n\n"), fp);

    for (i=0; i<st->n_cols; i++) {
	ct = st->cols[i];
	if (!err) {
	    fprintf(fp, M_("String code table for variable %d (%s):\n"), 
		    ct->idx, pdinfo->varname[ct->idx]);
	} else {
	    pprintf(prn, M_("String code table for variable %d (%s):\n"), 
		    ct->idx, pdinfo->varname[ct->idx]);
	}
	for (j=0; j<ct->n_strs; j++) {
	    if (!err) {
		fprintf(fp, "%3d = '%s'\n", j+1, ct->strs[j]);
	    } else {
		pprintf(prn, "%3d = '%s'\n", j+1, ct->strs[j]);
	    }
	}
    }

    if (fp != NULL) {
	pprintf(prn, M_("String code table written to\n %s\n"), stname);
	fclose(fp);
	set_string_table_written();
    }

    return err;
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
    { "shelldir", 0, NULL }
};

#ifdef WIN32
static saved_string dsep = { "dirsep", 0, "\\" };
#else
static saved_string dsep = { "dirsep", 0, "/" };
#endif

void gretl_insert_builtin_string (const char *name, const char *s)
{
    int n = sizeof built_ins / sizeof built_ins[0];
    int i, m;

    for (i=0; i<n; i++) {
	if (!strcmp(name, built_ins[i].name)) {
	    free(built_ins[i].s);
	    m = strlen(s);
	    if (s[m-1] == SLASH) {
		/* drop trailing dir separator for paths */
		built_ins[i].s = gretl_strndup(s, m - 1);
	    } else {
		built_ins[i].s = gretl_strdup(s);
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

/* access a saved string by name */

char *get_string_by_name (const char *name)
{
    int i, n, d;

    if (!strcmp(name, "dirsep")) {
	return dsep.s;
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

int delete_saved_string (const char *name, PRN *prn)
{
    saved_string *S;
    int err, builtin = 0;

    S = get_saved_string_by_name(name, &builtin);
    if (S == NULL) {
	err = E_UNKVAR;
    } else if (builtin) {
	sprintf(gretl_errmsg, _("You cannot delete '%s'"), name);
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

    clean_up_codevars();
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
	strcpy(gretl_errmsg, _("The shell command is not activated."));
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

char *gretl_getenv (const char *key, int *err)
{
    char *test = NULL;
    char *val = NULL;

    test = getenv(key);
    if (test == NULL) {
	val = gretl_strdup("");
    } else {
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

char *retrieve_date_string (int t, const DATAINFO *pdinfo, int *err)
{
    char datestr[OBSLEN] = {0};
    char *ret = NULL;

    if (t > 0 && t <= pdinfo->n) {
	ntodate(datestr, t - 1, pdinfo);
	ret = gretl_strdup(datestr);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    } else {
	*err = E_DATA;
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
	addpath(fullname, NULL, 0);

	g_file_get_contents(fullname, &ret, &len, &gerr);

	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	    *err = E_FOPEN;
	}
    } 

    return ret;
}

int string_is_defined (const char *sname)
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
	    sprintf(gretl_errmsg, _("You cannot overwrite '%s'\n"), name);
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

    if (prn != NULL && !err && gretl_messages_on() && *s != '\0') {
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
    sprintf(gretl_errmsg, _("Maximum length of command line "
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

    /* when using the genr apparatus for a string, let @foo be handled 
       as a variable (FIXME?) */
    if (!strncmp(line, "string ", 7)) {
	return 0;
    }

    if (!strncmp(line, "sscanf", 6)) {
	/* when scanning, let @foo be handled as a variable */
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
	    n = gretl_varchar_spn(s + 1);
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

/* setting of names of imported variables that should be treated
   as string codes */

static char **codevars;
static int n_codevars;

static void clean_up_codevars (void)
{
    if (codevars != NULL && n_codevars > 0) {
	free_strings_array(codevars, n_codevars);
	codevars = NULL;
	n_codevars = 0;
    }
}

int is_codevar (const char *s)
{
    int i;

    for (i=0; i<n_codevars; i++) {
	if (!strcmp(s, codevars[i])) {
	    return 1;
	}
    }

    return 0;
}

int set_codevars (const char *s)
{
    char chunk[32];
    const char *p;
    int err = 0;

    p = strstr(s, "codevars");
    if (p != NULL) {
	s = p + 9;
    }

    *chunk = '\0';
    sscanf(s, "%31s", chunk);

    if (*chunk == '\0') {
	err = E_DATA;
    } else {
	clean_up_codevars();
	if (strcmp(chunk, "null")) {
	    codevars = gretl_string_split(s, &n_codevars);
	    if (codevars == NULL) {
		err = E_ALLOC;
	    } 
	}
    }

    return err;
}
