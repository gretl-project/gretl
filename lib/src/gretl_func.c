/*
 *  Copyright (c) 2005 by Allin Cottrell
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
#include "monte_carlo.h"
#include "version.h"
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_xml.h"
#include "cmd_private.h"

#define FN_DEBUG 0
#define PKG_DEBUG 0

typedef struct fn_param_ fn_param;
typedef struct fncall_ fncall;

struct fn_param_ {
    char *name;
    char type;
    char flags;
    double deflt;
    double min;
    double max;
};

struct ufunc_ {
    char name[FN_NAMELEN];
    int pkgID;
    char *help;
    int private;
    int n_lines;
    char **lines;
    int n_params;
    fn_param *params;
    int rettype;
    char *retname;
    int in_use;
};

struct fnpkg_ {
    int ID;
    float minver;
    FuncDataReq dreq;
    char *fname;
    char *author;
    char *version;
    char *date;
    char *descrip;
};

enum {
    ARG_OPTIONAL = 1 << 0,
    ARG_CONST    = 1 << 1
};

#define ref_type(t) (t == ARG_REF_SCALAR || \
		     t == ARG_REF_SERIES || \
		     t == ARG_REF_MATRIX)

static int n_ufuns;
static ufunc **ufuns;
static ufunc *current_ufun;

static int n_pkgs;
static fnpkg **pkgs;

static int drop_function_vars = 1;

static void real_user_function_help (ufunc *fun, fnpkg *pkg, PRN *prn);

/* record of state, and communication of state with outside world */

static int compiling;
static int fn_executing;
static int fn_comment;

int gretl_compiling_function (void)
{
    return compiling;
}

static void set_compiling_on (void)
{
    fn_comment = 0;
    compiling = 1;
}

static void set_compiling_off (void)
{
    compiling = 0;
}

int gretl_function_depth (void)
{
    return fn_executing;
}

/* Switch to delete, or not, all the "private" variables internal to a
   function when execution terminates.  This is used in nls.c: if a
   user function will be called repeatedly in the context of
   NLS/MLE/GMM we turn off deletion temporarily, which saves repeated
   allocation/deallocation of storage for the variables.
*/

void set_drop_function_vars (int s)
{
    drop_function_vars = (s != 0)? 1 : 0;
}

int repeating_function_exec (void)
{
    return !drop_function_vars;
}

static void set_executing_on (ufunc *fun)
{
    fn_executing++;
    fun->in_use += 1;
#if FN_DEBUG
    fprintf(stderr, "set_executing_on: fun = %s, fn_executing=%d\n", 
	    fun->name, fn_executing);
#endif
}

static void set_executing_off (ufunc *fun)
{
    fn_executing--;
    fun->in_use -= 1;
#if FN_DEBUG
    fprintf(stderr, "set_executing_off: fun=%s, fn_executing=%d\n",
	    fun->name, fn_executing);
#endif
}

/* general info accessors */

int n_user_functions (void)
{
    return n_ufuns;
}

const ufunc *get_user_function_by_index (int idx)
{
    return (idx < 0 || idx >= n_ufuns)? NULL : ufuns[idx];
}

int fn_n_params (const ufunc *fun)
{
    return fun->n_params;
}

int fn_param_type (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? 0 :
	fun->params[i].type;
}

const char *fn_param_name (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NULL :
	fun->params[i].name;
}

double fn_param_default (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].deflt;
}    

double fn_param_minval (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].min;
}

double fn_param_maxval (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].max;
}

int fn_param_optional (const ufunc *fun, int i)
{
    int t;

    if (i < 0 || i >= fun->n_params) return 0;

    t = fun->params[i].type;

    return ((ref_type(t) || t == ARG_LIST) && 
	    (fun->params[i].flags & ARG_OPTIONAL));
}

int user_func_get_return_type (const ufunc *fun)
{
    if (fun == NULL) {
	return ARG_NONE;
    } else {
	return fun->rettype;
    }
}

const char *user_function_name_by_index (int i)
{
    if (i >= 0 && i < n_ufuns) {
	return ufuns[i]->name;
    } else {
	return NULL;
    }
}

int user_function_index_by_name (const char *name)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(name, ufuns[i]->name)) {
	    return i;
	}
    }

    return -1;
}

ufunc *get_user_function_by_name (const char *name)
{
    ufunc *fun = NULL;
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(name, (ufuns[i])->name)) {
	    fun = ufuns[i];
	    break;
	}
    }

#if FN_DEBUG
    if (fun != NULL) {
	fprintf(stderr, "get_user_function_by_name: name = '%s' (n_ufuns = %d);"
		" found match\n", name, n_ufuns);
    }
#endif

    return fun;
}

/* constructors and destructors */

static fn_param *allocate_params (int n)
{
    fn_param *params;
    int i;

    params = malloc(n * sizeof *params);
    if (params == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	params[i].name = NULL;
	params[i].type = 0;
	params[i].flags = 0;
	params[i].deflt = NADBL;
	params[i].min = NADBL;
	params[i].max = NADBL;
    }

    return params;
}

static ufunc *ufunc_new (void)
{
    ufunc *fun = malloc(sizeof *fun);

    if (fun == NULL) {
	return NULL;
    }

    fun->name[0] = '\0';
    fun->pkgID = 0;
    fun->private = 0;

    fun->help = NULL;

    fun->n_lines = 0;
    fun->lines = NULL;

    fun->n_params = 0;
    fun->params = NULL;

    fun->rettype = ARG_NONE;
    fun->retname = NULL;

    fun->in_use = 0;

    return fun;
}

static void free_params_array (fn_param *params, int n)
{
    int i;

    if (params == NULL) return;

    for (i=0; i<n; i++) {
	free(params[i].name);
    }
    free(params);
}

static void clear_ufunc_data (ufunc *fun)
{
    free(fun->help);
    free_strings_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    free(fun->retname);
    
    fun->help = NULL;
    fun->lines = NULL;
    fun->params = NULL;

    fun->n_lines = 0;
    fun->n_params = 0;

    fun->rettype = ARG_NONE;
    fun->retname = NULL;
    
    fun->in_use = 0;
}

static void free_ufunc (ufunc *fun)
{
    free(fun->help);
    free_strings_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    free(fun->retname);

    free(fun);
}

static int add_allocated_ufunc (ufunc *fun)
{
    int nf = n_ufuns;
    ufunc **myfuns;

    myfuns = realloc(ufuns, (nf + 1) * sizeof *myfuns);

    if (myfuns == NULL) {
	return E_ALLOC;
    }

    ufuns = myfuns;
    ufuns[nf] = fun;
    n_ufuns++;

#if PKG_DEBUG
    fprintf(stderr, "add_allocated_ufunc: name '%s', n_ufuns = %d\n",
	    fun->name, n_ufuns);
#endif

    return 0;
}

static ufunc *add_ufunc (void)
{
    ufunc *fun = ufunc_new();

    if (fun != NULL) {
	if (add_allocated_ufunc(fun)) {
	    free_ufunc(fun);
	    fun = NULL;
	}
    }

    return fun;
}

static void ufuncs_destroy (void)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	free_ufunc(ufuns[i]);
    }
    free(ufuns);

    ufuns = NULL;
    n_ufuns = 0;
}

/* handling of XML function packages */

enum {
    FUNCS_INFO,
    FUNCS_LOAD,
    FUNCS_CODE
};

static const char *arg_type_string (int type)
{
    if (type == ARG_BOOL) {
	return "bool";
    } else if (type == ARG_INT) {
	return "int";
    } else if (type == ARG_SCALAR) {
	return "scalar";
    } else if (type == ARG_SERIES) {
	return "series";
    } else if (type == ARG_LIST) {
	return "list";
    } else if (type == ARG_MATRIX) {
	return "matrix";
    } else if (type == ARG_REF_SCALAR) {
	return "scalar *";
    } else if (type == ARG_REF_SERIES) {
	return "series *";
    } else if (type == ARG_REF_MATRIX) {
	return "matrix *";
    } else {
	return "unknown";
    }
}

static int arg_type_from_string (const char *s)
{
    if (!strncmp(s, "bool", 4)) return ARG_BOOL;
    if (!strcmp(s, "int"))      return ARG_INT;
    if (!strcmp(s, "scalar"))   return ARG_SCALAR;
    if (!strcmp(s, "series"))   return ARG_SERIES;
    if (!strcmp(s, "list"))     return ARG_LIST;
    if (!strcmp(s, "matrix"))   return ARG_MATRIX;

    if (!strcmp(s, "scalar *"))  return ARG_REF_SCALAR;
    if (!strcmp(s, "series *"))  return ARG_REF_SERIES;
    if (!strcmp(s, "matrix *"))  return ARG_REF_MATRIX;

    if (!strcmp(s, "scalarref"))  return ARG_REF_SCALAR;
    if (!strcmp(s, "seriesref"))  return ARG_REF_SERIES;
    if (!strcmp(s, "matrixref"))  return ARG_REF_MATRIX;

    return 0;
}

static int field_to_type (const char *s)
{
#if FN_DEBUG
    fprintf(stderr, "field_to_type: looking at '%s'\n", s);
#endif

    if (isdigit(*s)) {
	return atoi(s);
    } else {
	return arg_type_from_string(s);
    }
}    

#define scalar_arg(t) (t == ARG_BOOL || t == ARG_INT || t == ARG_SCALAR)

static int func_read_params (xmlNodePtr node, ufunc *fun)
{
    xmlNodePtr cur;
    char *field;
    int n, err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n) || n < 0) {
	fprintf(stderr, "Couldn't read param count\n");
	return E_DATA;
    }  

    if (n == 0) {
	return 0;
    }

    fun->params = allocate_params(n);
    if (fun->params == NULL) {
	return E_ALLOC;
    } 

    fun->n_params = n;

    gretl_push_c_numeric_locale();

    cur = node->xmlChildrenNode;
    n = 0;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "param")) {
	    if (gretl_xml_get_prop_as_string(cur, "name", &field)) {
		fun->params[n].name = field;
	    } else {
		err = E_DATA;
	    }
	    if (gretl_xml_get_prop_as_string(cur, "type", &field)) {
		fun->params[n].type = field_to_type(field);
		free(field);
		if (scalar_arg(fun->params[n].type)) {
		    gretl_xml_get_prop_as_double(cur, "default", 
						 &fun->params[n].deflt);
		}
		if (fun->params[n].type == ARG_INT) {
		    gretl_xml_get_prop_as_double(cur, "min", 
						 &fun->params[n].min);
		    gretl_xml_get_prop_as_double(cur, "max", 
						 &fun->params[n].max);
		}
		if (gretl_xml_get_prop_as_bool(cur, "optional")) {
		    fun->params[n].flags |= ARG_OPTIONAL;
		}
		if (gretl_xml_get_prop_as_bool(cur, "const")) {
		    fun->params[n].flags |= ARG_CONST;
		}
	    } else {
		err = E_DATA;
	    }
	    n++;
	}	    
	cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    if (!err && n != fun->n_params) {
	err = E_DATA;
    }

    return err;
}

static int func_read_return (xmlNodePtr node, ufunc *fun)
{
    char *field;
    int err = 0;

    if (gretl_xml_get_prop_as_string(node, "name", &field)) {
	fun->retname = field;
    } else {
	fprintf(stderr, "return value: couldn't get name\n");
	err = E_DATA;
    }

    if (!err && gretl_xml_get_prop_as_string(node, "type", &field)) {
	fun->rettype = field_to_type(field);
	free(field);
    } else {
	fprintf(stderr, "return value: couldn't get type\n");
	err = E_DATA;
    }    

    return err;
}

static int func_read_code (xmlNodePtr node, xmlDocPtr doc, ufunc *fun,
			   PRN *prn)
{
    char line[MAXLINE];
    char *buf, *s;
    int err = 0;

    buf = (char *) xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
    if (buf == NULL) {
	return 1;
    }

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && !err) {
	if (string_is_blank(line)) {
	    continue;
	}
	if (prn != NULL) {
	    pprintf(prn, "  %s\n", line);
	} else {
	    s = line;
	    while (isspace(*s)) s++;
	    err = strings_array_add(&fun->lines, &fun->n_lines, s);
	}
    }

    bufgets_finalize(buf);

    free(buf);

    return err;
}

static void print_opt_flags (fn_param *param, PRN *prn)
{
    if (param->flags & ARG_OPTIONAL) {
	pputs(prn, "[null]");
    }
}

static void print_deflt_min_max (fn_param *param, PRN *prn)
{
    double x = param->min;
    double y = param->max;
    double z = param->deflt;

    if (na(x) && na(y) && na(z)) {
	return; /* no-op */
    } else if (na(x) && na(y) && !na(z)) {
	pprintf(prn, "[%g]", z);
	return;
    }

    pputc(prn, '[');
    if (na(x)) {
	pputc(prn, ':');
    } else {
	pprintf(prn, "%g:", x);
    }
    if (na(y)) {
	pputc(prn, ':');
    } else {
	pprintf(prn, "%g:", y);
    }  
    if (!na(z)) {
	pprintf(prn, "%g", z);
    }   
    pputc(prn, ']');
}

static void print_function_start (ufunc *fun, PRN *prn)
{
    const char *s;
    int i;

    gretl_push_c_numeric_locale();

    pprintf(prn, "function %s ", fun->name);
    for (i=0; i<fun->n_params; i++) {
	if (i == 0) {
	    pputc(prn, '(');
	}
	if (fun->params[i].flags & ARG_CONST) {
	    pputs(prn, "const ");
	}
	s = arg_type_string(fun->params[i].type);
	if (s[strlen(s) - 1] == '*') {
	    pprintf(prn, "%s%s", s, fun->params[i].name);
	} else {
	    pprintf(prn, "%s %s", s, fun->params[i].name);
	}
	if (fun->params[i].type == ARG_BOOL) {
	    if (!na(fun->params[i].deflt)) {
		pprintf(prn, "[%g]", fun->params[i].deflt);
	    }
	} else if (scalar_arg(fun->params[i].type)) {
	    print_deflt_min_max(&fun->params[i], prn);
	} else if (ref_type(fun->params[i].type)) {
	    print_opt_flags(&fun->params[i], prn);
	}
	if (i == fun->n_params - 1) {
	    pputc(prn, ')');
	} else {
	    pputs(prn, ", ");
	}
    }
    pputc(prn, '\n');

    gretl_pop_c_numeric_locale();
}

static void print_function_end (ufunc *fun, PRN *prn)
{
    const char *typestr;

    if (fun->rettype != ARG_NONE) {
	typestr = arg_type_string(fun->rettype);
	pprintf(prn, "  return %s %s\n", typestr, fun->retname);
    }

    pputs(prn, "end function\n");
}

static int read_ufunc_from_xml (xmlNodePtr node, xmlDocPtr doc, fnpkg *pkg, 
				int task, PRN *prn)
{
    ufunc *fun = ufunc_new();
    xmlNodePtr cur;
    char *fname;
    int err = 0;

    if (fun == NULL) {
	return E_ALLOC;
    }

    if (!gretl_xml_get_prop_as_string(node, "name", &fname)) {
	free_ufunc(fun);
	return E_DATA;
    }

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: got function name '%s'\n",
	    fname);
#endif

    *fun->name = 0;
    strncat(fun->name, fname, FN_NAMELEN - 1);
    free(fname);

    gretl_xml_get_prop_as_int(node, "private", &fun->private);
    if (pkg != NULL) {
	fun->pkgID = pkg->ID;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "help")) {
	    gretl_xml_node_get_string(cur, doc, &fun->help);
	} else if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    err = func_read_params(cur, fun);
	    if (err) {
		fprintf(stderr, "%s: error parsing function parameters\n",
			fun->name);
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "return")) {
	    err = func_read_return(cur, fun);
	    if (err) {
		fprintf(stderr, "%s: error parsing function return value\n",
			fun->name);
	    }
	} else if (task != FUNCS_INFO && 
		   !xmlStrcmp(cur->name, (XUC) "code")) {
	    if (task == FUNCS_CODE) {
		print_function_start(fun, prn);
	    }
	    err = func_read_code(cur, doc, fun, prn);
	}
	cur = cur->next;
    }

    if (task == FUNCS_LOAD) {
	if (!err) {
	   err = add_allocated_ufunc(fun);
	}
	if (err) {
	    free_ufunc(fun);
	}
    } else if (task == FUNCS_INFO) {
	if (!fun->private) {
	    real_user_function_help(fun, NULL, prn);
	}
	free_ufunc(fun);
    } else if (task == FUNCS_CODE) {
	print_function_end(fun, prn);
	free_ufunc(fun);
    }

#if FN_DEBUG
    if (err) {
	fprintf(stderr, "read_ufunc_from_xml: error reading spec\n");
    } 
#endif

    return err;
}

static void adjust_indent (const char *line, int *this_indent,
			   int *next_indent)
{
    int ti = *next_indent;
    int ni = *next_indent;

    if (!strncmp(line, "loop", 4)) {
	ni++;
    } else if (!strncmp(line, "if", 2)) {
	ni++;
    } else if (!strncmp(line, "nls", 3)) {
	ni++;
    } else if (!strncmp(line, "mle", 3)) {
	ni++;
    } else if (!strncmp(line, "end", 3)) {
	ti--;
	ni--;
    } else if (!strncmp(line, "else", 4)) {
	ni = ti;
	ti--;
    } 

    *this_indent = ti;
    *next_indent = ni;
}

static int write_function_xml (const ufunc *fun, FILE *fp)
{
    int this_indent = 0;
    int next_indent = 0;
    int i, j;

    fprintf(fp, "<gretl-function name=\"%s\" private=\"%d\">\n", 
	    fun->name, fun->private);

    if (fun->help != NULL) {
	gretl_xml_put_tagged_string("help", fun->help, fp);
    }

    if (fun->n_params > 0) {

	gretl_push_c_numeric_locale();

	fprintf(fp, " <params count=\"%d\">\n", fun->n_params);
	for (i=0; i<fun->n_params; i++) {
	    fprintf(fp, "  <param name=\"%s\" type=\"%s\"",
		    fun->params[i].name, 
		    arg_type_string(fun->params[i].type));
	    if (!na(fun->params[i].min)) {
		fprintf(fp, " min=\"%g\"", fun->params[i].min);
	    }
	    if (!na(fun->params[i].max)) {
		fprintf(fp, " max=\"%g\"", fun->params[i].max);
	    }
	    if (!na(fun->params[i].deflt)) {
		fprintf(fp, " default=\"%g\"", fun->params[i].deflt);
	    }
	    if (fun->params[i].flags & ARG_OPTIONAL) {
		fputs(" optional=\"true\"", fp);
	    }
	    if (fun->params[i].flags & ARG_CONST) {
		fputs(" const=\"true\"", fp);
	    }
	    fputs("/>\n", fp);
	}
	fputs(" </params>\n", fp);

	gretl_pop_c_numeric_locale();
    }

    if (fun->rettype != ARG_NONE) {
	fprintf(fp, " <return name=\"%s\" type=\"%s\"/>\n", fun->retname, 
		arg_type_string(fun->rettype));
    }

    fputs("<code>", fp);
    for (i=0; i<fun->n_lines; i++) {
	adjust_indent(fun->lines[i], &this_indent, &next_indent);
	for (j=0; j<this_indent; j++) {
	    fputs("  ", fp);
	}
	gretl_xml_put_raw_string(fun->lines[i], fp);
	fputc('\n', fp);
    }
    fputs("</code>\n", fp);

    fputs("</gretl-function>\n", fp);
    
    return 0;
}

int gretl_function_print_code (int i, PRN *prn)
{
    int this_indent = 0;
    int next_indent = 0;
    ufunc *fun;
    int j;
   
    if (i < 0 || i >= n_ufuns) {
	return 1;
    }

    fun = ufuns[i];

    print_function_start(fun, prn);

    for (i=0; i<fun->n_lines; i++) {
	adjust_indent(fun->lines[i], &this_indent, &next_indent);
	for (j=0; j<=this_indent; j++) {
	    pputs(prn, "  ");
	}
	pputs(prn, fun->lines[i]);
	pputc(prn, '\n');
    }  

    print_function_end(fun, prn);

    return 0;
}

int gretl_function_set_info (int i, const char *help)
{
    int err = 0;

    if (i >= 0 && i < n_ufuns) {
	free(ufuns[i]->help);
	if (help != NULL) {
	    ufuns[i]->help = gretl_strdup(help);
	} else {
	    ufuns[i]->help = NULL;
	}
    } else {
	err = 1;
    }

    return err;
}

static fnpkg *ufunc_get_parent_package (const ufunc *fun)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (fun->pkgID == pkgs[i]->ID) {
	    return pkgs[i];
	}
    }

    return NULL;
}

static char *make_pkgname (const char *fname)
{
    char *p = strrchr(fname, SLASH);
    char *ret;

    if (p != NULL) {
	ret = gretl_strdup(p + 1);
    } else {
	ret = gretl_strdup(fname);
    }

    if (ret == NULL) {
	return NULL;
    }

    p = strstr(ret, ".gfn");
    if (p != NULL) {
	*p = 0;
    }

    return ret;
}

int gretl_function_get_info (int i, const char *key, char const **value)
{
    if (i < 0 || i >= n_ufuns) {
	return E_DATA;
    }

    if (!strcmp(key, "help")) {
	*value = ufuns[i]->help;
    } else {
	fnpkg *pkg = ufunc_get_parent_package(ufuns[i]);

	if (pkg == NULL) {
	    *value = NULL;
	} else if (!strcmp(key, "author")) {
	    *value = pkg->author;
	} else if (!strcmp(key, "version")) {
	    *value = pkg->version;
	} else if (!strcmp(key, "date")) {
	    *value = pkg->date;
	} else if (!strcmp(key, "pkgdesc")) {
	    *value = pkg->descrip;
	}
    }

    return 0;
}

void gretl_function_set_private (int i, int priv)
{
    if (i >= 0 && i < n_ufuns) {
	ufuns[i]->private = (priv != 0);
    }
}

static void get_version_string (float ver, char *vstr)
{
    gretl_push_c_numeric_locale();
    sprintf(vstr, "%.2f", (double) ver);
    gretl_pop_c_numeric_locale();

    vstr[4] = vstr[3];
    vstr[3] = '.';
    vstr[5] = 0;
}

int write_function_package (fnpkg *pkg,
			    const char *fname,
			    int pub, 
			    const int *privlist,
			    const char *author,
			    const char *version,
			    const char *date,
			    const char *descrip,
			    FuncDataReq dreq,
			    float minver)
{
    char *pkgname;
    FILE *fp;
    int saveas = 0;
    int i, fi;
    int err = 0;

    if (n_ufuns == 0) {
	fprintf(stderr, "No functions are defined\n");
	return 0;
    }

    if (author == NULL || version == NULL || date == NULL || 
	descrip == NULL || pub < 0) {
	strcpy(gretl_errmsg, "Function information is incomplete");
	return E_DATA;
    }

    if (pkg != NULL && strcmp(fname, pkg->fname)) {
	saveas = 1;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	sprintf(gretl_errmsg, _("Couldn't open %s"), fname);
	return E_FOPEN;
    }

    gretl_xml_header(fp);    

    fputs("<gretl-functions>\n", fp);
    fputs("<gretl-function-package", fp);

    pkgname = make_pkgname(fname);
    if (pkgname != NULL) {
	fprintf(fp, " name=\"%s\"", pkgname);
	free(pkgname);
    } 

    if (dreq == FN_NEEDS_TS) {
	fprintf(fp, " %s=\"true\"", NEEDS_TS);
    } else if (dreq == FN_NEEDS_QM) {
	fprintf(fp, " %s=\"true\"", NEEDS_QM);
    } else if (dreq == FN_NEEDS_PANEL) {
	fprintf(fp, " %s=\"true\"", NEEDS_PANEL);
    }

    if (minver > 0) {
	char vstr[6];

	get_version_string(minver, vstr);
	fprintf(fp, " minver=\"%s\"", vstr);
    }

    fputs(">\n", fp);

    gretl_xml_put_tagged_string("author", author, fp);
    gretl_xml_put_tagged_string("version", version, fp);
    gretl_xml_put_tagged_string("date", date, fp);
    gretl_xml_put_tagged_string("description", descrip, fp);

    write_function_xml(ufuns[pub], fp);

    if (privlist != NULL) {
	for (i=1; i<=privlist[0]; i++) {
	    fi = privlist[i];
	    if (fi >= 0 && fi < n_ufuns) {
		write_function_xml(ufuns[fi], fp);
	    }
	}
    }

    fputs("</gretl-function-package>\n", fp);
    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    if (pkg != NULL && !saveas) {
	/* existing package, name has not changed: 
	   update package info */
	if (strcmp(author, pkg->author)) {
	    free(pkg->author);
	    pkg->author = gretl_strdup(author);
	}
	if (strcmp(version, pkg->version)) {
	    free(pkg->version);
	    pkg->version = gretl_strdup(version);
	}
	if (strcmp(date, pkg->date)) {
	    free(pkg->date);
	    pkg->date = gretl_strdup(date);
	}
	if (strcmp(descrip, pkg->descrip)) {
	    free(pkg->descrip);
	    pkg->descrip = gretl_strdup(descrip);
	}
	pkg->dreq = dreq;
	pkg->minver = minver;

	if (pkg->author == NULL || pkg->version == NULL ||
	    pkg->date == NULL || pkg->descrip == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

int function_package_get_info (const char *fname,
			       fnpkg **ppkg,
			       int *pub, 
			       int **privlist,
			       char **author,
			       char **version,
			       char **date,
			       char **descrip,
			       FuncDataReq *dreq,
			       float *minver)
{
    fnpkg *pkg = NULL;
    int pubnum = -1;
    int npriv = 0;
    int *list = NULL;
    int i, j, err = 0;

    if (n_pkgs == 0 || n_ufuns == 0) {
	fprintf(stderr, "function_package_get_info: no functions loaded\n");
	return 1;
    }

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    pkg = pkgs[i];
	    break;
	}
    }

    if (pkg == NULL) {
	fprintf(stderr, "No package associated with '%s'\n", fname);
	return 1;
    }

    if (ppkg != NULL) {
	*ppkg = pkg;
    }
    if (author != NULL) {
	*author = gretl_strdup(pkg->author);
    }
    if (date != NULL) {
	*date = gretl_strdup(pkg->date);
    }
    if (version != NULL) {
	*version = gretl_strdup(pkg->version);
    }
    if (descrip != NULL) {
	*descrip = gretl_strdup(pkg->descrip);
    }
    if (dreq != NULL) {
	*dreq = pkg->dreq;
    }
    if (minver != NULL) {
	*minver = pkg->minver;
    }

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkgID == pkg->ID) {
	    if (ufuns[i]->private) {
		npriv++;
	    } else {
		pubnum = i;
	    }
	}
    }

    if (!err && pub != NULL && pubnum >= 0) {
	*pub = pubnum;
    }

    if (!err && privlist != NULL && npriv > 0) {
	list = gretl_list_new(npriv);
	if (list != NULL) {
	    j = 1;
	    for (i=0; i<n_ufuns; i++) {
		if (ufuns[i]->pkgID == pkg->ID && ufuns[i]->private) {
		    list[j++] = i;
		}
	    }	    
	    *privlist = list;
	} else {
	    err = E_ALLOC;
	}
    }
	    
    return err;
}

static void print_function_package (fnpkg *pkg, FILE *fp)
{
    int i;

    fputs("<gretl-function-package>\n", fp);

    if (pkg->author != NULL) {
	gretl_xml_put_tagged_string("author", pkg->author, fp);
    }
    if (pkg->version != NULL) {
	gretl_xml_put_tagged_string("version", pkg->version, fp);
    }
    if (pkg->date != NULL) {
	gretl_xml_put_tagged_string("date", pkg->date, fp);
    }
    if (pkg->descrip != NULL) {
	gretl_xml_put_tagged_string("description", pkg->descrip, fp);
    }

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkgID == pkg->ID) {
	    write_function_xml(ufuns[i], fp);
	}
    }

    fputs("</gretl-function-package>\n", fp);
}

/* dump all currently defined functions */

int write_user_function_file (const char *fname)
{
    FILE *fp;
    int i;

    if (n_ufuns == 0) {
	return 0;
    }

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_xml_header(fp);  
    fputs("<gretl-functions>\n", fp);

    for (i=0; i<n_pkgs; i++) {
	print_function_package(pkgs[i], fp);
    }

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkgID == 0) {
	    write_function_xml(ufuns[i], fp);
	}
    }

    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    return 0;
}

static int function_package_add (fnpkg *pkg)
{
    fnpkg **tmp;
    int err = 0;

    tmp = realloc(pkgs, (n_pkgs + 1) * sizeof *tmp);
    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	pkgs = tmp;
	pkgs[n_pkgs] = pkg;
	n_pkgs++;
    }

    return err;
}

static fnpkg *function_package_new (const char *fname)
{
    fnpkg *pkg = malloc(sizeof *pkg);

    if (pkg == NULL) {
	return NULL;
    }

    pkg->ID = n_pkgs + 1; /* FIXME? */
    pkg->author = NULL;
    pkg->version = NULL;
    pkg->date = NULL;
    pkg->descrip = NULL;
    pkg->dreq = 0;
    pkg->minver = 0;

    pkg->fname = gretl_strdup(fname);
    if (pkg->fname == NULL) {
	free(pkg);
	pkg = NULL;
    }

    return pkg;
}

static void function_package_free (fnpkg *pkg)
{
    if (pkg != NULL) {
	int i;

	for (i=0; i<n_ufuns; i++) {
	    if (ufuns[i]->pkgID == pkg->ID) {
		ufuns[i]->pkgID = 0;
	    }
	}

	free(pkg->fname);
	free(pkg->author);
	free(pkg->version);
	free(pkg->date);
	free(pkg->descrip);
	free(pkg);
    }
}

static void packages_destroy (void)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	function_package_free(pkgs[i]);
    }
    free(pkgs);

    pkgs = NULL;
    n_pkgs = 0;
}

static float version_float_from_string (const char *s)
{
    int maj, min, pl;

    sscanf(s, "%d.%d.%d", &maj, &min, &pl);
    return maj + min / 10.0 + pl / 100.0;
}

static int 
read_user_function_package (xmlDocPtr doc, xmlNodePtr node, 
			    const char *fname, int task, 
			    PRN *prn, char **pname)
{
    xmlNodePtr cur;
    fnpkg *pkg;
    char *vstr = NULL;
    int err = 0;

    pkg = function_package_new(fname);
    if (pkg == NULL) {
	return E_ALLOC;
    }

    if (pname != NULL) {
	gretl_xml_get_prop_as_string(node, "name", pname);
    }

    if (gretl_xml_get_prop_as_bool(node, NEEDS_TS)) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_QM)) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_PANEL)) {
	pkg->dreq = FN_NEEDS_PANEL;
    }

    if (gretl_xml_get_prop_as_string(node, "minver", &vstr)) {
	pkg->minver = version_float_from_string(vstr);
	free(vstr);
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "author")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->author);
	} else if (!xmlStrcmp(cur->name, (XUC) "version")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->version);
	} else if (!xmlStrcmp(cur->name, (XUC) "date")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->date);
	} else if (!xmlStrcmp(cur->name, (XUC) "description")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->descrip);
	}
	cur = cur->next;
    }

    if (task == FUNCS_INFO) {
	if (pname != NULL && *pname != NULL) {
	    pprintf(prn, "Package: %s\n", *pname);
	} else {
	    pprintf(prn, "Package: %s\n", fname);
	}
	pprintf(prn, "Author: %s\n", (pkg->author)? pkg->author : "unknown");
	pprintf(prn, "Version: %s\n", (pkg->version)? pkg->version : "unknown");
	pprintf(prn, "Date: %s\n", (pkg->date)? pkg->date : "unknown");
	pprintf(prn, "Description: %s\n", (pkg->descrip)? pkg->descrip : "none");
    }

    /* now get specific functions from this package */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
	    if (prn != NULL) {
		pputc(prn, '\n');
	    }
	    err = read_ufunc_from_xml(cur, doc, pkg, task, prn);
	}
	cur = cur->next;
    }

    if (err || prn != NULL) {
	function_package_free(pkg);
    } else {
	err = function_package_add(pkg);
    }

#if PKG_DEBUG
    fprintf(stderr, "read_user_function_package:\n"
	    " err = %d reading '%s'\n", err, fname);
#endif

    return err;
}

/* if prn is non-NULL, we're just reading the contents
   of this file in order to display them */

static int real_read_user_function_file (const char *fname, int task, PRN *prn,
					 char **pname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    int err = 0;

    xmlKeepBlanksDefault(0);

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return err;
    }

    /* first get any function packages from this file */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
	    err = read_user_function_package(doc, cur, fname, task, prn, pname);
	} 
	cur = cur->next;
    }

    /* then get any unpackaged functions */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
	    if (prn != NULL) {
		pputc(prn, '\n');
	    }
	    err = read_ufunc_from_xml(cur, doc, NULL, task, prn);
	}
	cur = cur->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

#if PKG_DEBUG
    fprintf(stderr, "real_read_user_function_file:\n"
	    " returning %d for '%s'\n", err, fname);
#endif

    return err;
}

int user_function_file_is_loaded (const char *fname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    return 1;
	}
    }

    return 0;
}

/* read functions from file into gretl's workspace */

int load_user_function_file (const char *fname)
{
    return real_read_user_function_file(fname, FUNCS_LOAD, NULL, NULL);
}

/* read specific function info from file, but do not
   load into workspace */

int get_function_file_info (const char *fname, PRN *prn, char **pname)
{
    return real_read_user_function_file(fname, FUNCS_INFO, prn, pname);
}

int get_function_file_code (const char *fname, PRN *prn, char **pname)
{
    return real_read_user_function_file(fname, FUNCS_CODE, prn, pname);
}

char *get_function_file_header (const char *fname, int *err)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr sub;
    char *descrip = NULL;

    xmlKeepBlanksDefault(0);

    *err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (*err) {
	return NULL;
    }

    node = node->xmlChildrenNode;
    while (node != NULL) {
	if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
	    sub = node->xmlChildrenNode;
	    while (sub != NULL) {
		if (!xmlStrcmp(sub->name, (XUC) "description")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, &descrip);
		    break;
		}
		sub = sub->next;
	    }
	    if (descrip != NULL) break;
	}
	node = node->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

    if (descrip == NULL) {
	descrip = gretl_strdup(_("No description available"));
    }

    return descrip;
}

static char *function_name_from_line (const char *line, char *name)
{
    if (*line == '#' || !strncmp(line, "(*", 2)) {
	*name = '\0';
    } else {
	const char *p = strchr(line, '=');
	char *q;

	if (p == NULL) {
	    p = line;
	} else {
	    p++;
	}

	sscanf(p, "%31s", name);
	q = strchr(name, '(');
	if (q != NULL) {
	    *q = 0;
	}
#if FN_DEBUG
	fprintf(stderr, "function_name_from_line: line='%s', got '%s'\n", 
		line, name);
#endif
    }

    return name;
}

int gretl_is_user_function (const char *line)
{
    int ret = 0;

    if (n_ufuns > 0 && !string_is_blank(line)) {
	char name[FN_NAMELEN];

#if FN_DEBUG > 1
	fprintf(stderr, "gretl_is_user_function: testing '%s'\n", line);
#endif
	function_name_from_line(line, name);
	if (get_user_function_by_name(name) != NULL) {
	    ret = 1;
	}
    }

    return ret;
}

int gretl_is_public_user_function (const char *name)
{
    ufunc *fun = get_user_function_by_name(name);

    if (fun != NULL && !fun->private) {
	return 1;
    } else {
	return 0;
    }
}

int gretl_get_user_function (const char *line)
{
    int ret = 0;

    if (n_ufuns > 0 && !string_is_blank(line)) {
	char name[FN_NAMELEN];

#if FN_DEBUG > 1
	fprintf(stderr, "gretl_is_user_function: testing '%s'\n", line);
#endif
	function_name_from_line(line, name);
	if (get_user_function_by_name(name) != NULL) {
	    ret = 1;
	} else if (function_from_string(name)) {
	    ret = 1;
	}
    }

    return ret;
}

static void delete_ufunc_from_list (ufunc *fun)
{
    if (n_ufuns == 0 || fun == NULL) {
	/* "can't happen" */
	return;
    }

    if (n_ufuns == 1) {
	free_ufunc(fun);
	free(ufuns);
	ufuns = NULL;
	n_ufuns = 0;
    } else {
	int i, gotit = 0;

	for (i=0; i<n_ufuns; i++) {
	    if (gotit) {
		ufuns[i-1] = ufuns[i];
	    }
	    if (!gotit && ufuns[i] == fun) {
		gotit = 1;
		free_ufunc(fun);
	    }
	}
	
	ufuns = realloc(ufuns, (n_ufuns - 1) * sizeof *ufuns);
	n_ufuns--;
    }
}

static int check_func_name (const char *fname, ufunc **pfun, PRN *prn)
{
    int i, err = 0;

    if (!isalpha((unsigned char) *fname)) {
	strcpy(gretl_errmsg, "function names must start with a letter");
	err = 1;
    } else if (gretl_command_number(fname)) {
	sprintf(gretl_errmsg, "'%s' is the name of a gretl command",
		fname);
	err = 1;
    } else {
	for (i=0; i<n_ufuns; i++) {
	    if (!strcmp(fname, ufuns[i]->name)) {
		pprintf(prn, "Redefining function '%s'\n", fname);
		if (pfun != NULL) {
		    clear_ufunc_data(ufuns[i]);
		    *pfun = ufuns[i];
		} else {
		    delete_ufunc_from_list(ufuns[i]);
		}
		break;
	    }
	}
    }
    
    return err;
}

static int maybe_delete_function (const char *fname)
{
    ufunc *fun = get_user_function_by_name(fname);
    int err = 0;

    if (fun == NULL) {
	; /* no-op */
    } else if (fun->in_use) {
	sprintf(gretl_errmsg, "%s: function is in use", fname);
	err = 1;
    } else {
	delete_ufunc_from_list(fun);
    } 

    return err;
}

static int comma_count (const char *s)
{
    int nc = 0;

    while (*s) {
	if (*s == ',') nc++;
	s++;
    }

    return nc;
}

static int read_deflt_min_max (char *s, fn_param *param,
			       int *namelen)
{
    char *p = strchr(s, '[');
    int err = 0;

    if (p != NULL) {
	double x, y, z;

	if (param->type == ARG_BOOL) {
	    if (sscanf(p, "[%lf]", &x) == 1) {
		param->deflt = x;
	    } else {
		err = E_DATA;
	    }
	} else {
	    if (sscanf(p, "[%lf:%lf:%lf]", &x, &y, &z) == 3) {
		param->min = x;
		param->max = y;
		param->deflt = z;
	    } else if (sscanf(p, "[%lf::%lf]", &x, &y) == 2) {
		param->min = x;
		param->deflt = y;
	    } else if (sscanf(p, "[:%lf:%lf]", &x, &y) == 2) {
		param->max = x;
		param->deflt = y;
	    } else if (sscanf(p, "[%lf]", &x) == 1) {
		param->deflt = x;
	    } else {
		err = E_DATA;
	    }
	}
	if (!err) {
	    *namelen = p - s;
	}
    }

    return err;
}

static int read_param_option (char *s, fn_param *param,
			      int *namelen)
{
    char *p = strstr(s, "[null]");

    if (p != NULL) {
	param->flags |= ARG_OPTIONAL;
	*namelen -= 6;
    }

    return 0;
}

static int parse_function_param (char *s, fn_param *param, int i)
{
    char tstr[22] = {0};
    char *name;
    int type, len;
    int err = 0;

    while (isspace(*s)) s++;

    if (!strncmp(s, "const ", 6)) {
	param->flags |= ARG_CONST;
	s += 6;
	while (isspace(*s)) s++;
    }

    /* get arg or return type */
    len = strcspn(s, " ");
    if (len > 21) {
	err = E_PARSE;
    } else {
	strncat(tstr, s, len);
	s += len;
	while (isspace(*s)) s++;
	if (*s == '*') {
	    strcat(tstr, " *");
	    s++;
	}
	type = arg_type_from_string(tstr);
	if (type == 0) {
	    sprintf(gretl_errmsg, "Unrecognized data type '%s'", tstr);
	    err = E_PARSE;
	} 
    }

    if (err) {
	return err;
    }
    
    while (isspace(*s)) s++;
    len = strcspn(s, " ");
    if (len == 0) {
	sprintf(gretl_errmsg, "parameter %d: name is missing", i);
	err = E_PARSE;
    }

    if (err) {
	return err;
    }    

    if (scalar_arg(type)) {
	param->type = type;
	err = read_deflt_min_max(s, param, &len);
    }

    if (ref_type(type) || type == ARG_LIST) {
	param->type = type;
	err = read_param_option(s, param, &len);
    }    

    if (!err) {
	name = gretl_strndup(s, len);
	if (name == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	param->type = type;
	param->name = name;
    }

#if FN_DEBUG
    if (!err) {
	fprintf(stderr, " param[%d] = '%s', ptype = %d\n", 
		i, name, type);
	fprintf(stderr, "min=%g, max=%g, deflt=%g\n", 
		param->min, param->max, param->deflt);
    }
#endif

    return err;
}

static void arg_tail_strip (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>=0; i--) {
	if (isspace(s[i]) || s[i] == ')') {
	    s[i] = '\0';
	} else {
	    break;
	}
    }
}

static int parse_fn_definition (char *fname, 
				fn_param **pparams,
				int *pnp,
				const char *str, 
				ufunc **pfun, 
				PRN *prn)
{
    fn_param *params = NULL;
    char fmt[8] = "%31s";
    char *p, *s = NULL;
    int i, len, np = 0;
    int err = 0;

    while (isspace(*str)) str++;

    /* get the function name */
    len = strcspn(str, " (");
    if (len == 0) {
	err = E_PARSE;
    } else if (len < FN_NAMELEN - 1) {
	sprintf(fmt, "%%%ds", len);
    }

    if (!err) {
	if (sscanf(str, fmt, fname) != 1) {
	    err = E_PARSE;
	}
    }

    if (!err) {
	err = check_func_name(fname, pfun, prn);
    }

    if (!err) {
	str += len;
    }

    if (*str == '\0') {
	/* void function */
	return 0;
    }

    /* move to next bit and make a copy */
    if (!err) {
	str += strspn(str, " (");
	if (*str == 0) {
	    err = E_PARSE;
	} else {
	    s = gretl_strdup(str);
	    if (s == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!strcmp(s, ")")) {
	/* void function "foo()" */
	free(s);
	return 0;
    }

    if (!err) {
	/* strip trailing ')' and space */
	arg_tail_strip(s);
	np = comma_count(s) + 1;
	if (np == 1 && !strcmp(s, "void")) {
	    free(s);
	    return 0;
	}
	params = allocate_params(np);
	if (params == NULL) {
	    err = E_ALLOC;
	}
	if (!err) {
	    charsub(s, ',', 0);
	}
    }

    if (!err) {
	p = s;
	for (i=0; i<np && !err; i++) {
	    err = parse_function_param(p, &params[i], i);
	    p += strlen(p) + 1;
	}
    }

    free(s);

    if (err) {
	free_params_array(params, np);
    } else {
	*pparams = params;
	*pnp = np;
    }
    
    return err;
}

int gretl_start_compiling_function (const char *line, PRN *prn)
{
    ufunc *fun = NULL;
    fn_param *params = NULL;
    int nf, n_params = 0;
    char fname[FN_NAMELEN];
    char extra[8];
    int err = 0;

    nf = sscanf(line, "function %31s %7s", fname, extra);
    if (nf <= 0) {
	return E_PARSE;
    } 

    if (nf == 2) {
	if (!strcmp(extra, "clear") || !strcmp(extra, "delete")) {
	    maybe_delete_function(fname);
	    return 0;
	}
    } 

    /* the following takes care of replacing an existing function
       of the same name, if any */

    err = parse_fn_definition(fname, &params, &n_params,
			      line + 8, &fun, prn);

    if (!err && fun == NULL) {
	fun = add_ufunc();
	if (fun == NULL) {
	    free_params_array(params, n_params);
	    err = E_ALLOC;
	}
    }

    if (!err) {
	strcpy(fun->name, fname);
	fun->params = params;
	fun->n_params = n_params;
	current_ufun = fun;
	set_compiling_on();
    } else {
	current_ufun = NULL;
    }
    
    return err;
}

static int add_function_return (ufunc *fun, const char *line)
{
    char s1[16], s2[VNAMELEN];
    int type;
    int err = 0;

    if (fun->rettype != ARG_NONE) {
	sprintf(gretl_errmsg, "Function %s: return value is already defined",
		fun->name);
	return 1;
    }

    if (sscanf(line, "%15s %15s", s1, s2) != 2) {
	return E_PARSE;
    }

    type = field_to_type(s1);
    if (type == 0) {
	return E_PARSE;
    } 

    err = check_varname(s2);

    if (!err) {
	fun->retname = gretl_strdup(s2);
	if (fun->retname == NULL) {
	    err = E_ALLOC;
	} else {
	    fun->rettype = type;
	}
    }

    return err;
}

static int end_of_function (const char *s)
{
    int ret = 0;

    if (!strncmp(s, "end ", 4)) {
	char word[9];

	if (sscanf(s + 4, "%8s", word) && !strcmp(word, "function")) {
	    ret = 1;
	}
    }

    return ret;
}

int gretl_function_append_line (const char *line)
{
    ufunc *fun = current_ufun;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "gretl_function_append_line: '%s'\n", line);
#endif

    if (fun == NULL) {
#if FN_DEBUG
	fprintf(stderr, " fun == NULL!\n");
#endif
	return 1;
    }

    if (string_is_blank(line)) {
	return 0;
    }

    if (!strncmp(line, "/*", 2)) {
	fn_comment = 1;
    } else if (!strncmp(line, "*/", 2)) {
	fn_comment = 0;
    }

    if (fn_comment) {
	return strings_array_add(&fun->lines, &fun->n_lines, line);
    }

    if (end_of_function(line)) {
	if (fun->n_lines == 0) {
	    sprintf(gretl_errmsg, "%s: empty function", fun->name);
	    delete_ufunc_from_list(fun);
	    err = 1;
	}
	set_compiling_off();
	return err;
    }

    if (!strncmp(line, "quit", 4)) {
	/* abort compilation */
	delete_ufunc_from_list(fun);
	set_compiling_off();
	return err;
    }

    if (!strncmp(line, "function", 8)) {
	strcpy(gretl_errmsg, "You can't define a function within a function");
	return 1;
    }

    if (!strncmp(line, "return ", 7)) {
	err = add_function_return(fun, line + 7);
    } else {  
	err = strings_array_add(&fun->lines, &fun->n_lines, line);
    }

    return err;
}

int update_function_from_script (const char *fname, int idx)
{
    char line[MAXLINE];
    char *s;
    FILE *fp;
    int gotfn = 0;
    int err = 0;

    if (idx < 0 || idx >= n_ufuns) {
	return E_DATA;
    }

    fp = fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fprintf(stderr, "Going to update function id %d from %s\n",
	    idx, fname);

    while (fgets(line, sizeof line, fp) && !err) {
	s = line;
	while (*s == ' ') s++;
	tailstrip(s);
	if (!strncmp(s, "function ", 9)) {
	    if (gotfn) {
		err = 1;
	    } else {
		gotfn = 1;
		err = gretl_start_compiling_function(s, NULL);
		strcpy(gretl_errmsg, "Error compiling function");
	    }
	} else {
	    err = gretl_function_append_line(s);
	    strcpy(gretl_errmsg, "Error compiling function");
	}
    }

    fclose(fp);

    if (!err && current_ufun != NULL) {
	int ichk = user_function_index_by_name(current_ufun->name);

	if (ichk != idx) {
	    strcpy(gretl_errmsg, "Function name has been changed!");
	    fprintf(stderr, "idx = %d, but user_function_index_by_name() "
		    "gives %d for '%s'\n", idx, ichk, current_ufun->name);
	    err = 1;
	}
    }

    return err;
}

/* Given a named list of variables supplied as an argument to a
   function, copy the list under the name assigned by the function,
   and make the variables referenced in that list accessible to the
   function.  If the "const" flag is set for the parameter, mark the
   variables in question as const.
*/

static int localize_list (const char *oldname, fn_param *fp,
			  DATAINFO *pdinfo)
{
    const int *list;
    int i, err;

    list = get_list_by_name(oldname);

    if (list == NULL) {
	err = E_DATA;
    } else {
	err = copy_named_list_as(oldname, fp->name);
    }

    if (!err) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] != 0) {
		STACK_LEVEL(pdinfo, list[i]) += 1;
		if (fp->flags & ARG_CONST) {
		    set_var_const(pdinfo, list[i]);
		}
	    }
	}
    }

    return err;
}

/* Scalar function arguments only: if the arg is not supplied, use the
   default that is contained in the function specification, if any.
*/

static int add_scalar_arg_default (fn_param *param, double ***pZ,
				   DATAINFO *pdinfo)
{
    double x;

    if (na(param->deflt)) {
	/* should be impossible here, but... */
	return E_DATA;
    }

    if (param->type == ARG_BOOL || param->type == ARG_INT) {
	x = floor(param->deflt);
    } else {
	x = param->deflt;
    }
    
    return dataset_add_scalar_as(x, param->name, pZ, pdinfo);
}

static int allocate_function_args (ufunc *fun,
				   int argc, fnargs *args, 
				   double ***pZ,
				   DATAINFO *pdinfo)
{
    fn_param *fp;
    int xi = 0, Xi = 0, Mi = 0, li = 0;
    int i, err = 0;

    /* regular arguments, passed by value */
    
    for (i=0; i<fun->n_params && !err; i++) {
	fp = &fun->params[i];
	if (scalar_arg(fp->type)) {
	    if (xi >= args->nx) {
		err = add_scalar_arg_default(fp, pZ, pdinfo);
	    } else {
		err = dataset_add_scalar_as(args->x[xi++], fp->name, 
					    pZ, pdinfo);
	    } 
	} else if (fp->type == ARG_SERIES) {
	    err = dataset_add_series_as(args->X[Xi++], fp->name, 
					pZ, pdinfo);
	} else if (fp->type == ARG_MATRIX) {
	    err = copy_matrix_as(args->M[Mi++], fp->name);
	} else if (fp->type == ARG_LIST) {
	    if (li < args->nl) {
		err = localize_list(args->lists[li++], fp, pdinfo);
	    } else {
		err = create_named_null_list(fp->name);
	    } 
	} 
    }

    xi = Mi = 0;

    /* "pointer" parameters, passed by reference: make these vars
       visible at function level, under their parameter names
    */

    for (i=0; i<argc && !err; i++) {
	fp = &fun->params[i];
	if (ref_type(fp->type)) {
	    if (args->types[i] == ARG_REF_SCALAR ||
		args->types[i] == ARG_REF_SERIES) {
		err = strings_array_add(&args->upnames, &args->nnames, 
					pdinfo->varname[args->refv[xi]]);
		if (!err) {
		    STACK_LEVEL(pdinfo, args->refv[xi]) += 1;
		    strcpy(pdinfo->varname[args->refv[xi]], fp->name);
		    xi++;
		}
	    } else if (args->types[i] == ARG_REF_MATRIX) {
		err = strings_array_add(&args->upnames, &args->nnames, 
					user_matrix_get_name(args->refm[Mi]));
		if (!err) {
		    user_matrix_adjust_level(args->refm[Mi], 1);
		    user_matrix_set_name(args->refm[Mi], fp->name);
		    Mi++;
		}
	    } 
	} 
    }

    return err;
}

int check_function_needs (const DATAINFO *pdinfo, FuncDataReq dreq,
			  float minver)
{
    static float thisver = 0.0;

    if (thisver == 0.0) {
	thisver = version_float_from_string(GRETL_VERSION);
    }

    if (minver > thisver) {
	char vstr[6];

	get_version_string(minver, vstr);
	sprintf(gretl_errmsg, "This function needs gretl version %s", vstr);
	return 1;
    }

    if ((dreq == FN_NEEDS_TS) && 
	!dataset_is_time_series(pdinfo)) {
	strcpy(gretl_errmsg, "This function needs time-series data");
	return 1;
    }

    if ((dreq == FN_NEEDS_PANEL) && 
	!dataset_is_panel(pdinfo)) {
	strcpy(gretl_errmsg, "This function needs panel data");
	return 1;
    }

    if ((dreq == FN_NEEDS_QM) && 
	(!dataset_is_time_series(pdinfo) || 
	 (pdinfo->pd != 4 || pdinfo->pd != 12))) {
	strcpy(gretl_errmsg, "This function needs quarterly or monthly data");
	return 1;
    } 

    return 0;
}

static int 
maybe_check_function_needs (const DATAINFO *pdinfo,
				 const ufunc *fun)
{
    const fnpkg *pkg = ufunc_get_parent_package(fun);

    if (pkg == NULL) {
	return 0;
    } else {
	return check_function_needs(pdinfo, pkg->dreq, pkg->minver);
    }
}

enum {
    GET_PTR,
    GET_COPY
};

static double 
get_scalar_return (const char *vname, double **Z, DATAINFO *pdinfo,
		   int *err)
{
    int v = varindex(pdinfo, vname);

    if (v < pdinfo->v && var_is_scalar(pdinfo, v)) {
	return Z[v][0];
    } else if (v < pdinfo->v) {
	*err = E_TYPES;
    } else {
	*err = E_UNKVAR;
    }

    return NADBL;
}

static double *
get_series_return (const char *vname, double **Z, DATAINFO *pdinfo,
		   int action, int *err)
{
    int v = varindex(pdinfo, vname);
    double *x = NULL;

    if (v < pdinfo->v && var_is_series(pdinfo, v)) {
	if (action == GET_COPY) {
	    x = copyvec(Z[v], pdinfo->n);
	    if (x == NULL) {
		*err = E_ALLOC;
	    }
	} else {
	    x = Z[v];
	}
    } else if (v < pdinfo->v) {
	*err = E_TYPES;
    } else {
	*err = E_UNKVAR;
    }

    return x;
}

static gretl_matrix *
get_matrix_return (const char *mname, int action, int *err)
{
    gretl_matrix *m = get_matrix_by_name(mname);
    gretl_matrix *ret = NULL;

    if (m != NULL) {
	if (action == GET_COPY) {
	    ret = gretl_matrix_copy(m);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    }
	} else {
	    return ret = m;
	}
    } else {
	*err = E_UNKVAR;
    }

    return ret;
}

static int unlocalize_list (const char *listname, DATAINFO *pdinfo)
{
    const int *list = get_list_by_name(listname);
    int d = gretl_function_depth();
    int i, err = 0;

    if (list == NULL) {
	err = E_DATA;
    } else {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] != 0 && list[i] < pdinfo->v) {
		if (STACK_LEVEL(pdinfo, list[i]) == d) {
		    STACK_LEVEL(pdinfo, list[i]) -= 1;
		}
		unset_var_const(pdinfo, list[i]);
	    }
	}
    }

    return err;
}

static char *
get_list_return (const char *lname, DATAINFO *pdinfo, int *err)
{
    
    char *ret = gretl_strdup(lname);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	*err = unlocalize_list(lname, pdinfo);
	if (!*err) {
	    *err = named_list_lower_level(lname);
	}
    }

    return ret;
}

static int 
function_assign_returns (ufunc *u, fnargs *args, int argc, int rtype, 
			 double **Z, DATAINFO *pdinfo, 
			 void *ret, PRN *prn, int *perr)
{
    fn_param *fp;
    int vi = 0, mi = 0, li = 0;
    int i, j, err = 0;

    if (*perr == 0) {
	/* direct return value */
	if (rtype == ARG_SCALAR) {
	    *(double *) ret = get_scalar_return(u->retname, Z, pdinfo, &err);
	} else if (rtype == ARG_SERIES) {
	    *(double **) ret = get_series_return(u->retname, Z, pdinfo, GET_COPY, &err);
	} else if (rtype == ARG_MATRIX) {
	    *(gretl_matrix **) ret = get_matrix_return(u->retname, GET_COPY, &err);
	} else if (rtype == ARG_LIST) {
	    *(char **) ret = get_list_return(u->retname, pdinfo, &err);
	}

	if (err == E_UNKVAR) {
	    pprintf(prn, "Function %s did not provide the specified return value\n",
		    u->name);
	}
	*perr = err;
    }

    /* "indirect return" values: these should be restored even if the
       function bombed
    */

    j = 0;
    for (i=0; i<argc; i++) {
	fp = &u->params[i];
	if (ref_type(fp->type)) {
	    if (args->types[i] == ARG_REF_SCALAR ||
		args->types[i] == ARG_REF_SERIES) {
		STACK_LEVEL(pdinfo, args->refv[vi]) -= 1;
		strcpy(pdinfo->varname[args->refv[vi]], args->upnames[j++]);
		vi++;
	    } else if (args->types[i] == ARG_REF_MATRIX) {
		user_matrix_adjust_level(args->refm[mi], -1);
		user_matrix_set_name(args->refm[mi], args->upnames[j++]);
		mi++;
	    }
	} else if (fp->type == ARG_LIST) {
	    if (li < args->nl) {
		unlocalize_list(args->lists[li++], pdinfo);
	    } 
	}
    }

    return err;
}

static int stop_fncall (ufunc *u, double ***pZ, DATAINFO *pdinfo,
			int orig_v)
{
    int d = gretl_function_depth();
    int anyerr = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "stop_fncall: terminating call to "
	    "function '%s' at depth %d\n", u->name, d);
#endif

    anyerr = destroy_saved_lists_at_level(d);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_saved_lists_at_level(%d): err = %d\n", d, err);
#endif
    }

    anyerr = destroy_user_matrices_at_level(d);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_user_matrices_at_level(%d): err = %d\n", d, err);
#endif
    }  

    /* below: delete some or all local variables, taking care not to
       delete any 'local' vars that have been "promoted" to caller
       level via their inclusion in a returned list
    */

    if (drop_function_vars) {
	/* delete all local variables */
	int i, delv = 0;

	for (i=orig_v; i<pdinfo->v; i++) {
	    if (STACK_LEVEL(pdinfo, i) == d) {
		delv++;
	    }
	}

	if (delv == pdinfo->v - orig_v) {
	    anyerr = dataset_drop_last_variables(delv, pZ, pdinfo);
	    if (anyerr && !err) {
		err = anyerr;
	    }
	} else {
	    for (i=orig_v; i<pdinfo->v; i++) {
		if (STACK_LEVEL(pdinfo, i) == d) {
		    anyerr = dataset_drop_variable(i, pZ, pdinfo);
		    if (anyerr && !err) {
			err = anyerr;
		    }
		}
	    }
	}
    } else {
	/* delete only variables copied from arguments */
	fn_param *fp;
	int i, v;

	for (i=0; i<u->n_params && !err; i++) {
	    fp = &u->params[i];
	    if (scalar_arg(fp->type) || fp->type == ARG_SERIES) {
		v = varindex(pdinfo, fp->name);
		if (STACK_LEVEL(pdinfo, v) == d) {
		    anyerr = dataset_drop_variable(v, pZ, pdinfo);
		    if (anyerr && !err) {
			err = anyerr;
		    }
		}
	    } 
	}		
    }

    pop_program_state();
    set_executing_off(u);

    return err;
}

static void start_fncall (ufunc *u)
{
    set_executing_on(u);
    push_program_state();
    set_gretl_echo(0);
    set_gretl_messages(0);
}

static char funcerr_msg[256];

const char *get_funcerr_message (void)
{
    return funcerr_msg;
}

static void set_funcerr_message (ufunc *u, const char *s)
{
    int n;

    sprintf(funcerr_msg, _("Error message from %s:\n"), u->name);
    n = strlen(funcerr_msg);
    strncat(funcerr_msg, s, 255 - n);
}

int gretl_function_exec (ufunc *u, fnargs *args, int rtype,
			 double ***pZ, DATAINFO *pdinfo,
			 void *ret, PRN *prn)
{
    ExecState state;
    MODEL **models = NULL;
    char line[MAXLINE];
    CMD cmd;
    fn_param *fp;
    int started = 0;
    int funcerr = 0;
    int argc, i, j;
    int err = 0;

    int orig_v = pdinfo->v;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;

    *funcerr_msg = '\0';

    err = maybe_check_function_needs(pdinfo, u);
    if (err) {
	return err;
    }

    /* precautions */
    cmd.list = NULL;
    cmd.param = NULL;
    cmd.extra = NULL;
    cmd.linfo = NULL;
    state.cmd = NULL;
    state.models = NULL;

    argc = args->nx + args->nX + args->nM + args->nl +
	args->nrefv + args->nrefm + args->nnull;

#if FN_DEBUG
    fprintf(stderr, "gretl_function_exec: argc = %d\n", argc);
    fprintf(stderr, "u->n_params = %d\n", u->n_params);
#endif

    j = 0;
    for (i=0; i<argc && !err; i++) {
	fp = &u->params[i];
	if ((fp->flags & ARG_OPTIONAL) && args->types[i] == ARG_NONE) {
	    ; /* this is OK */
	} else if (scalar_arg(fp->type) && args->types[i] == ARG_SCALAR) {
	    ; /* this is OK too */
	} else if (fp->type != args->types[i]) {
	    pprintf(prn, "argv[%d] is of wrong type (got %s, should be %s)\n", 
		    i, arg_type_string(args->types[i]), 
		    arg_type_string(fp->type));
	    err = E_TYPES;
	} else if (fp->type == ARG_SCALAR) {
	    if ((!na(fp->min) && args->x[j] < fp->min) ||
		(!na(fp->max) && args->x[j] > fp->max)) {
		pprintf(prn, "argv[%d]: scalar value %g out of bounds\n", 
			i, args->x[j]);
		err = E_INVARG;
	    }
	    j++;
	}
    }

    for (i=argc; i<u->n_params && !err; i++) {
	/* do we have defaults for any empty args? */
	fp = &u->params[i];
	if (!(fp->flags & ARG_OPTIONAL) && na(fp->deflt)) {
	    pprintf(prn, "%s: not enough arguments\n", u->name);
	    err = E_ARGS;
	}
    }

    if (!err) {
	err = allocate_function_args(u, argc, args, pZ, pdinfo);
    }

    if (err) {
	/* get out before allocating further storage */
	return err;
    }    

    models = allocate_working_models(2);
    if (models == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_cmd_init(&cmd);
    }

    if (!err) {
	*line = '\0';
	gretl_exec_state_init(&state, FUNCTION_EXEC, line, &cmd, 
			      models, prn);
	if (pdinfo->submode) {
	    state.subinfo = pdinfo;
	}
    }

    if (!err) {
	start_fncall(u);
	started = 1;
    }

    /* get function lines in sequence and check, parse, execute */

    for (i=0; i<u->n_lines && !err; i++) {
	strcpy(line, u->lines[i]);
	err = maybe_exec_line(&state, pZ, &pdinfo, &funcerr);
	if (funcerr) {
	    pprintf(prn, "%s: %s\n", u->name, state.cmd->param);
	    set_funcerr_message(u, state.cmd->param);
	}
	if (gretl_execute_loop()) { 
	    err = gretl_loop_exec(&state, pZ, &pdinfo);
	    if (err) {
		fprintf(stderr, "function_exec: breaking on error in loop\n");
		break;
	    }
	}
    }

    /* restore the sample that was in place on entry */

    if (complex_subsampled()) {
	if (state.subinfo == NULL) {
	    /* we weren't sub-sampled on entry: easy */
	    restore_full_sample(pZ, &pdinfo, NULL);
	} else if (state.subinfo != pdinfo) {
	    /* we're differently sub-sampled: complex */
	    restore_full_sample(pZ, &pdinfo, &state);
	} 
    } 

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    function_assign_returns(u, args, argc, rtype, *pZ, pdinfo,
			    ret, prn, &err);

    gretl_exec_state_clear(&state);

    if (started) {
	int stoperr = stop_fncall(u, pZ, pdinfo, orig_v);

	if (stoperr && !err) {
	    err = stoperr;
	}
    }

#if FN_DEBUG
    fprintf(stderr, "gretl_function_exec: err = %d\n", err);
#endif

    return err;    
}

void gretl_functions_cleanup (void)
{
    ufuncs_destroy();
    packages_destroy();
}

static void real_user_function_help (ufunc *fun, fnpkg *pkg, PRN *prn)
{
    int i;

    if (pkg == NULL) {
	pkg = ufunc_get_parent_package(fun);
    }

    pprintf(prn, "function %s\n\n", fun->name);

    if (pkg != NULL) {
	pprintf(prn, "Author: %s\n", pkg->author? pkg->author : "unknown");
	pprintf(prn, "Version: %s\n", pkg->version? pkg->version : "unknown");
	pprintf(prn, "Date: %s\n\n", pkg->date? pkg->date : "unknown");
    }

    if (fun->n_params > 0) {
	pprintf(prn, "Parameters:\n");
	for (i=0; i<fun->n_params; i++) {
	    pprintf(prn, " %s (%s)\n", 
		    fun->params[i].name, arg_type_string(fun->params[i].type));
	}
	pputc(prn, '\n');
    } else {
	pputs(prn, "Parameters: none\n\n");
    }

    if (fun->rettype != ARG_NONE) {
	pprintf(prn, "Return value: %s (%s)\n\n", 
		fun->retname, arg_type_string(fun->rettype));
    } else {
	pputs(prn, "Return value: none\n\n");
    }
	
    if (fun->help != NULL) {
	pprintf(prn, "Help text:\n%s\n\n", fun->help);
    }
}

int user_function_help (const char *fnname, PRN *prn)
{
    ufunc *fun = get_user_function_by_name(fnname);
    int err = 0;

    if (fun == NULL) {
	pprintf(prn, _("\"%s\" is not defined.\n"), fnname);
	err = 1;
    } else {
	real_user_function_help(fun, NULL, prn);
    }

    return err;
}
