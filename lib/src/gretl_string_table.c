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

/* for use in the context of the "printf" and "sprintf" commands:
   access a saved string by name */

char *get_named_string (const char *name)
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

static int append_to_newstr (char **newstr, char *s, int offset)
{
    char *tmp;
    int n;

    if (*newstr != NULL) {
	n = strlen(*newstr) + strlen(s) + 1;
    } else {
	n = strlen(s) + 1;
    }

    tmp = malloc(n - offset);

    if (tmp == NULL) {
	return E_ALLOC;
    }

    if (*newstr != NULL) {
	strcpy(tmp, *newstr);
	free(*newstr);
    } else {
	*tmp = '\0';
    }

    strcat(tmp, s + offset);
    *newstr = tmp;

    return 0;
}

static int append_to_saved_string (saved_string *str, char *s)
{
    char *tmp;
    int n;

    if (str->s != NULL) {
	n = strlen(str->s) + strlen(s) + 1;
    } else {
	n = strlen(s) + 1;
    }

    tmp = malloc(n);

    if (tmp == NULL) {
	return E_ALLOC;
    }

    if (str->s != NULL) {
	strcpy(tmp, str->s);
	free(str->s);
    } else {
	*tmp = '\0';
    }

    strcat(tmp, s);
    str->s = tmp;
    
    return 0;
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

static char *retrieve_string_var (const char **pline, 
				  const char *s,
				  int *err)
{
    char sname[VNAMELEN] = {0};
    char *p, *ret = NULL;
    int n, spn;

    spn = n = gretl_varchar_spn(s);
    if (n >= VNAMELEN) {
	n = VNAMELEN - 1;
    }

    strncat(sname, s, n);
    p = get_named_string(sname);

    if (p == NULL) {
	*err = E_UNKVAR;
    } else {
	ret = gretl_strdup(p);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    *pline = s + spn;
	}
    }

    return ret;
}

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

    arg += strspn(arg, " \t");

    if (strchr(arg, '@')) {
	char tmp[MAXLEN];

	strcpy(tmp, arg);
	substitute_named_strings(tmp);
	gretl_shell_grab(tmp, sout);
    } else {
	gretl_shell_grab(arg, sout);
    }

    if (sout != NULL && *sout != NULL) {
	/* trim trailing newline */
	int n = strlen(*sout);

	if ((*sout)[n-1] == '\n') {
	    (*sout)[n-1] = '\0';
	}
    }

    return err;
}

static char *gretl_backtick (const char **pline, int *err)
{
    const char *s = *pline + 2;
    char *arg = NULL;
    char *val = NULL;
    int n = 0, pcount = 1;

    while (*s && pcount) {
	if (*s == '(') pcount++;
	else if (*s == ')') pcount--;
	s++;
	n++;
    }

    if (pcount || n < 2) {
	*err = E_PARSE;
	return NULL;
    }

    arg = gretl_strndup(*pline + 2, n - 1);
    if (arg == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *pline = s;
    *err = shell_grab(arg, &val);

    free(arg);
    
    return val;
}

static char *gretl_getenv (const char **pline, int *err)
{
    const char *s = *pline;
    char *p = strchr(s, '"');
    char *key = NULL;
    char *test = NULL;
    char *val = NULL;
    int n;

    if (p == NULL || *(p+1) != ')') {
	*err = E_PARSE;
	return NULL;
    }

    n = strcspn(s, "\"");
    key = gretl_strndup(s, n);
    if (key == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *pline = p + 2;

    test = getenv(key);
    if (test == NULL) {
	val = gretl_strdup("");
    } else {
	val = gretl_strdup(test);
    }

    if (val == NULL) {
	*err = E_ALLOC;
    }

    free(key);
    
    return val;
}

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

static char *gretl_strstr (const char **pline, int *err)
{
    char *haystack = NULL, *needle = NULL;
    const char *s, *p = *pline;
    char *tmp, *ret = NULL;
    int n;

    p += strspn(p, " \t\r\n");
    s = p;

    if (*p == '"') {
	s++;
	p++;
	p = strchr(p, '"');
	if (p != NULL) {
	    p = strchr(p, ',');
	}
    } else {
	p = strchr(p, ',');
    }

    if (p == NULL) {
	*err = E_PARSE;
    } else {
	n = p - s;
	haystack = gretl_strndup(s, n);
	p++;
	p += strspn(p, " \t\r\n");
	s = p;
	if (*p == '"') {
	    p++;
	    s++;
	    p = strchr(p, '"');
	    if (p != NULL) {
		p = strrchr(p, ')');
	    }
	} else {
	    p = strrchr(p, ')');
	}
	if (p == NULL) {
	    *err = E_PARSE;
	} else {
	    n = p - s;
	    needle = gretl_strndup(s, n);
	}
    }

    if (haystack != NULL && needle != NULL) {

	strim(haystack);
	if (*haystack == '@') {
	    tmp = get_named_string(haystack + 1);
	    if (tmp != NULL) {
		free(haystack);
		haystack = gretl_strdup(tmp);
	    } 
	}

	strim(needle);
	if (*needle == '@') {
	    tmp = get_named_string(needle + 1);
	    if (tmp != NULL) {
		free(needle);
		needle = gretl_strdup(tmp);
	    }
	}
    } 

    if (haystack != NULL && needle != NULL) {
	tmp = strstr(haystack, needle);
	if (tmp != NULL) {
	    ret = gretl_strdup(tmp);
	} else {
	    ret = gretl_strdup("");
	}
    }	

    if (!*err && p != NULL) {
	*pline = ++p;
    }

    free(haystack);
    free(needle);
	
    return ret;
}

static char *retrieve_date_string (const char **pline, 
				   double ***pZ, DATAINFO *pdinfo,
				   int *err)
{
    const char *s = *pline;
    int p = 1, n = 0;
    char *ret = NULL;

    while (*s) {
	if (*s == '(') p++;
	else if (*s == ')') p--;
	if (p == 0) {
	    break;
	}
	n++;
	s++;
    }

    if (p == 0 && n > 0) {
	char datestr[OBSLEN] = {0};
	int t;

	ret = gretl_strndup(*pline, n);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    *pline = s + 1;
	    t = generate_scalar(ret, pZ, pdinfo, err);
	    if (!*err && t > 0 && t <= pdinfo->n) {
		ntodate(datestr, t - 1, pdinfo);
		free(ret);
		ret = gretl_strdup(datestr);
	    } else if (!*err) {
		*err = E_DATA;
	    }
	    if (*err) {
		free(ret);
		ret = NULL;
	    }
	}
    } else {
	*err = E_PARSE;
    }

    return ret;
}

static char *retrieve_arg_name (const char **pline, int *err)
{
    char argvar[VNAMELEN];
    char *ret = NULL;

    if (sscanf(*pline, "%15[^)])", argvar) == 1) {
	ret = gretl_func_get_arg_name(argvar);
	*pline += strlen(argvar);
	if (**pline == ')') *pline += 1;
    } else {
	*err = E_PARSE;
    }

    return ret;
}

static char *retrieve_file_content (const char **pline, int *err)
{
    char fmt[16], fname[FILENAME_MAX];
    char *ret = NULL;

    if (**pline == '"') {
	*pline += 1;
    }

    sprintf(fmt, "%%%d[^)])", FILENAME_MAX - 1);
    *fname = '\0';

    if (strchr(*pline, '@')) {
	char tmp[MAXLINE];

	strcpy(tmp, *pline);
	substitute_named_strings(tmp);
	sscanf(tmp, fmt, fname);
    } else {
	sscanf(*pline, fmt, fname);
    }

    if (*fname != '\0') {
	int n = strlen(fname);
	GError *gerr = NULL;
	gsize len = 0;

	if (fname[n-1] == '"') {
	    fname[n-1] = '\0';
	}

	addpath(fname, NULL, 0);

	*pline = strchr(*pline, ')') + 1;
	g_file_get_contents(fname, &ret, &len, &gerr);
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	    *err = E_FOPEN;
	}
    } else {
	*err = E_PARSE;
    }

    return ret;
}

static char *get_string_element (const char **pline, 
				 double ***pZ, DATAINFO *pdinfo,
				 int *err)
{
    const char *line = *pline;
    const char *s;
    char *cpy;
    int closed = 0;
    int bs = 0;
    int n = 0;

    line += strspn(line, " \t");
    *pline = line;

    if (!strncmp(line, "null", 4) && string_is_blank(line + 4)) {
	*pline += 4;
	return gretl_strdup("");
    }

    if (*line == '@') {
	/* should be saved string variable */
	return retrieve_string_var(pline, line + 1, err);
    }

    if (!strncmp(line, "getenv(\"", 8)) {
	*pline += 8;
	return gretl_getenv(pline, err);
    }

    if (!strncmp(line, "argname(", 8)) {
	*pline += 8;
	return retrieve_arg_name(pline, err);
    }

    if (!strncmp(line, "obslabel(", 9)) {
	*pline += 9;
	return retrieve_date_string(pline, pZ, pdinfo, err);
    }    

    if (!strncmp(line, "readfile(", 9)) {
	*pline += 9;
	return retrieve_file_content(pline, err);
    }  

    if (!strncmp(line, "strstr(", 7)) {
	*pline += 7;
	return gretl_strstr(pline, err);
    }    

    if (!strncmp(line, "$(", 2)) {
	return gretl_backtick(pline, err);
    }    

    if (*line != '"') {
	*err = E_PARSE;
	return NULL;
    }

    line++;
    s = line;
    while (*s) {
	/* allow for escaped quotes */
	if (*s == '\\') {
	    bs++;
	} else if (*s == '"' && (bs % 2 == 0)) {
	    closed = 1;
	    break;
	}
	if (*s != '\\') {
	    bs = 0;
	}
	s++;
	n++;
    }

    if (!closed) {
	fprintf(stderr, "quote not closed\n");
	*err = E_PARSE;
	return NULL;
    }

    cpy = malloc(n + 1);

    if (cpy == NULL) {
	*err = E_ALLOC;
    } else {
	copy_unescape(cpy, line, n);
    }

    s++; /* eat closing quote */
    *pline = s;

    return cpy;
}

int string_is_defined (const char *sname)
{
    saved_string *str;
    int builtin = 0;

    if (*sname == '@' && *(sname + 1) != '@') {
	sname++;
    }
    
    str = get_saved_string_by_name(sname, &builtin);

    return (str != NULL && str->s != NULL && str->s[0]);
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

/* for use in "sprintf" command */

int save_named_string (const char *name, const char *s, PRN *prn)
{
    saved_string *str;
    int builtin = 0;

    if (s == NULL) {
	return E_DATA;
    }

    str = get_saved_string_by_name(name, &builtin);
    
    if (str != NULL && builtin) {
	pprintf(prn, _("You cannot overwrite '%s'\n"), name);
	return E_DATA;
    }

    if (str == NULL) {
	str = add_named_string(name);
	if (str == NULL) {
	    return E_ALLOC;
	}
    }

    if (str->s != NULL) {
	free(str->s);
    }

    str->s = gretl_strdup(s);
    if (str->s == NULL) {
	return E_ALLOC;
    }

    if (gretl_messages_on()) {
	if (str->s[0] == '\0') {
	    pprintf(prn, _("Saved empty string as '%s'\n"), name);
	} else {
	    pprintf(prn, _("Saved string as '%s'\n"), name);
	}
    }

    return 0;
}

/* allow a string element to be of the form "str + n", to take
   an offset into str (provided the offset is not greater
   than the length of str)
*/

static int get_plus_mod (char *s1, const char **pline, int *plus)
{
    const char *s = *pline;
    int err = 0;

    s += strspn(s, " \t\r\n");
    if (*s == '+') {
	s++;
	s += strspn(s, " \t\r\n");
	if (isdigit(*s)) {
	    char *end = NULL;

	    *plus = (int) strtol(s, &end, 10);
	    if (*plus < 0 || *plus > strlen(s1)) {
		err = E_DATA;
	    } else {
		*pline = end;
	    }
	}
    }

    return err;
}

/* respond to commands of the forms:

     string <name> = "<s1>" "<s2>" ... "<sn>"
     string <name> += "<s2>" "<s2>"  ... "<sn>"
*/

int process_string_command (const char *line, double ***pZ,
			    DATAINFO *pdinfo, PRN *prn)
{
    saved_string *str;
    char *newstr = NULL;
    char targ[VNAMELEN];
    int builtin = 0;
    int n, add = 0;
    int err = 0;

    /* skip "string" plus any following space */
    line += 6;
    line += strspn(line, " \t");

    n = gretl_varchar_spn(line);
    if (n == 0 || n >= VNAMELEN) {
	return E_PARSE;
    }

    *targ = '\0';
    strncat(targ, line, n);
    line += n;

    /* eat any space before operator */
    line += strspn(line, " \t");

    if (*line == '\0') {
	/* just a call to echo an existing string? */
	str = get_saved_string_by_name(targ, &builtin);
	if (str == NULL) {
	    return E_UNKVAR;
	} else {
	    pprintf(prn, " %s\n", str->s);
	    return 0;
	}
    }

    /* operator must be '=' or '+=' */
    if (!strncmp(line, "+=", 2)) {
	add = 1;
    } else if (*line != '=') {
	return E_PARSE;
    }

    line += (add)? 2 : 1;

    /* set up the target */
    str = get_saved_string_by_name(targ, &builtin);

    if (str != NULL && builtin) {
	pprintf(prn, _("You cannot overwrite '%s'\n"), targ);
	return E_DATA;
    }	

    if (str == NULL) {
	if (add) {
	    return E_UNKVAR;
	} else {
	    str = add_named_string(targ);
	    if (str == NULL) {
		return E_ALLOC;
	    }
	}
    } 

    /* add strings(s) to target */
    while (!err && *line != '\0') {
	char *s1 = NULL;
	int plus = 0;

	s1 = get_string_element(&line, pZ, pdinfo, &err);
	if (!err) {
	    err = get_plus_mod(s1, &line, &plus);
	    if (!err) {
		err = append_to_newstr(&newstr, s1, plus);
	    }
	    free(s1);
	}
    }

    if (!err) {
	if (add) {
	    err = append_to_saved_string(str, newstr);
	    free(newstr);
	} else {
	    free(str->s);
	    str->s = newstr;
	}
    }

    if (!err && gretl_messages_on()) {
	if (str->s[0] == '\0') {
	    pprintf(prn, _("Saved empty string as '%s'\n"), targ);
	} else {
	    pprintf(prn, _("Saved string as '%s'\n"), targ);
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

int substitute_named_strings (char *line)
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

    if (!strncmp(line, "string", 6)) {
	/* when defining a string, let @foo be handled as a variable */
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
