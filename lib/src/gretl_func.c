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
#include "monte_carlo.h"
#include "version.h"
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_xml.h"
#include "cmd_private.h"
#include "gretl_string_table.h"

#define FN_DEBUG 0
#define PKG_DEBUG 0
#define EXEC_DEBUG 0
#define UDEBUG 0

typedef struct fn_param_ fn_param;

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
    int *in_use;
    fnargs *args;
};

struct fnpkg_ {
    int ID;
    char name[FN_NAMELEN];
    char *fname;
    char *author;
    char *version;
    char *date;
    char *descrip;
    float minver;
    FuncDataReq dreq;
    ufunc *iface;
    ufunc **priv;
    int n_priv;
};

struct fnarg {
    int type;
    char *upname;  /* name of arg at caller level */
    union {
	int idnum;
	double x;
	double *px;
	gretl_matrix *m;
	user_matrix *um;
	char *str;
    } val;
};

struct fnargs_ {
    int argc;
    struct fnarg **arg;
};

enum {
    ARG_OPTIONAL = 1 << 0,
    ARG_CONST    = 1 << 1
};

static int n_ufuns;
static ufunc **ufuns;
static ufunc *current_ufun;

static int n_pkgs;
static fnpkg **pkgs;

static int drop_function_vars = 1;

static void real_user_function_help (ufunc *fun, int ci, PRN *prn);
static void delete_ufunc_from_list (ufunc *fun);
static fnpkg *function_package_new (const char *fname);
static int function_package_add (fnpkg *pkg);
static int function_package_remove_by_ID (int ID);

/* record of state, and communication of state with outside world */

static int compiling;
static int fn_executing;

int gretl_compiling_function (void)
{
    return compiling;
}

static void set_compiling_on (void)
{
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

/* handling of function arguments */

struct fnarg *fn_arg_new (int type, void *p, int *err)
{
    struct fnarg *arg;

    arg = malloc(sizeof *arg);
    if (arg == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    arg->type = type;
    arg->upname = NULL;
    
    if (type == GRETL_TYPE_NONE) {
	arg->val.x = 0;
    } else if (type == GRETL_TYPE_DOUBLE) {
	arg->val.x = *(double *) p;
    } else if (type == GRETL_TYPE_SERIES) {
	arg->val.px = (double *) p;
    } else if (type == GRETL_TYPE_MATRIX) {
	arg->val.m = (gretl_matrix *) p;
    } else if (type == GRETL_TYPE_LIST ||
	       type == GRETL_TYPE_STRING) {
	arg->val.str = (char *) p;
    } else if (type == GRETL_TYPE_SCALAR_REF ||
	       type == GRETL_TYPE_SERIES_REF ||
	       type == GRETL_TYPE_UVAR) {
	arg->val.idnum = * (int *) p;
    } else if (type == GRETL_TYPE_MATRIX_REF) {
	arg->val.um = (user_matrix *) p;
    } else {
	*err = E_TYPES;
	free(arg);
	arg = NULL;
    }

    return arg;
}

fnargs *fn_args_new (void)
{
    fnargs *args = malloc(sizeof *args);

    if (args == NULL) {
	return NULL;
    }

    args->argc = 0;
    args->arg = NULL;

    return args;
}

/* note: this is not supposed to be a "deep free" (in case the arg
   carries a pointer member); that is handled in geneval.c 
*/

static void free_fn_arg (struct fnarg *arg)
{
    free(arg->upname);
    free(arg);
}

void fn_args_free (fnargs *args)
{
    int i;

    for (i=0; i<args->argc; i++) {
	free_fn_arg(args->arg[i]);
    }
    
    free(args->arg);
    free(args);
}

int push_fn_arg (fnargs *args, int type, void *p)
{
    int err = 0;

#if 0
    fprintf(stderr, "push_fn_arg: starting on type %d\n", type);
#endif

    if (args == NULL) {
	err = E_DATA;
    } else {
	struct fnarg **arg;
	int n = args->argc + 1;

	arg = realloc(args->arg, n * sizeof *arg);
	if (arg == NULL) {
	    err = E_ALLOC;
	} else {
	    args->arg = arg;
	    args->arg[n-1] = fn_arg_new(type, p, &err);
	}

	if (!err) {
	    args->argc = n;
	}
    }

    return err;
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
    gretl_list_append_term(&fun->in_use, fn_executing);
#if EXEC_DEBUG
    fprintf(stderr, "set_executing_on: fun = %s, fn_executing=%d\n", 
	    fun->name, fn_executing);
#endif
}

static int use_list_delete_last (ufunc *fun)
{
    if (fun->in_use == NULL || fun->in_use[0] == 0) {
	/* shouldn't happen */
	return E_DATA;
    }

    if (fun->in_use[0] == 1) {
	free(fun->in_use);
	fun->in_use = NULL;
    } else {
	int *newlist = gretl_list_new(fun->in_use[0] - 1);
	int i;

	if (newlist == NULL) {
	    return E_ALLOC;
	}

	for (i=1; i<fun->in_use[0]; i++) {
	    newlist[i] = fun->in_use[i];
	}
	free(fun->in_use);
	fun->in_use = newlist;
    }

    return 0;
}

static void set_executing_off (ufunc *fun)
{
    use_list_delete_last(fun);
    fn_executing--;
#if EXEC_DEBUG
    fprintf(stderr, "set_executing_off: fun=%s, fn_executing=%d\n",
	    fun->name, fn_executing);
#endif
}

/* general info accessors */

int n_free_functions (void)
{
    int i, n = 0;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkgID == 0) {
	    n++;
	}
    }

    return n;
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

    return ((gretl_ref_type(t) || t == GRETL_TYPE_LIST) && 
	    (fun->params[i].flags & ARG_OPTIONAL));
}

int user_func_get_return_type (const ufunc *fun)
{
    if (fun == NULL) {
	return GRETL_TYPE_NONE;
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

static int fname_idx;

const char *next_free_function_name (void)
{
    const char *ret = NULL;
    ufunc *fun;

    if (n_ufuns == 0) {
	fname_idx = 0;
	return NULL;
    }

    while (fname_idx < n_ufuns) {
	fun = ufuns[fname_idx++];
	if (fun->pkgID == 0) {
	    ret = fun->name;
	    break;
	}
    }

    return ret;
}

void function_names_init (void)
{
    fname_idx = 0;
}

int current_func_pkgID (void)
{
    ufunc *fun;
    int i, f0;

    for (i=0; i<n_ufuns; i++) {
	fun = ufuns[i];
	if (fun->in_use != NULL) {
	    f0 = fun->in_use[0];
	    if (fun->in_use[f0] == fn_executing) {
		return fun->pkgID;
	    }
	}
    }

    return 0;
}

static ufunc *currently_called_function (void)
{
    ufunc *fun;
    int i, f0;

    for (i=0; i<n_ufuns; i++) {
	fun = ufuns[i];
	if (fun->in_use != NULL) {
	    f0 = fun->in_use[0];
	    if (fun->in_use[f0] == fn_executing) {
		return fun;
	    }
	}
    }

    return NULL;
}

ufunc *get_user_function_by_name (const char *name)
{
    ufunc *fun = NULL;
    int i, ID;

    if (n_ufuns == 0) {
	return NULL;
    }

    ID = current_func_pkgID();

    /* if we're currently exec'ing a function from a package,
       look first for functions from the same package */

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(name, (ufuns[i])->name)) {
	    fun = ufuns[i];
	    if (fun->pkgID == ID || ID == 0) {
		break;
	    } else {
		fun = NULL;
	    }
	}
    }

    if (ID > 0 && fun == NULL) {
	/* fall back on unpackaged functions */
	for (i=0; i<n_ufuns; i++) {
	    if (!strcmp(name, (ufuns[i])->name)) {
		fun = ufuns[i];
		if (fun->pkgID == 0) {
		    break;
		} else {
		    fun = NULL;
		}
	    }
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

    fun->rettype = GRETL_TYPE_NONE;
    fun->retname = NULL;

    fun->in_use = NULL;
    fun->args = NULL;

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

    fun->rettype = GRETL_TYPE_NONE;
    fun->retname = NULL;
    
    fun->in_use = NULL;
    fun->args = NULL;
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

static ufunc *add_ufunc (const char *fname)
{
    ufunc *fun = ufunc_new();

    if (fun != NULL) {
	strncat(fun->name, fname, FN_NAMELEN - 1);
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

static const char *arg_type_string (int t)
{
    switch (t) {
    case GRETL_TYPE_BOOL:       return "bool";
    case GRETL_TYPE_INT:        return "int";
    case GRETL_TYPE_DOUBLE:     return "scalar";
    case GRETL_TYPE_SERIES:     return "series";
    case GRETL_TYPE_LIST:       return "list";
    case GRETL_TYPE_MATRIX:     return "matrix";
    case GRETL_TYPE_SCALAR_REF: return "scalar *";
    case GRETL_TYPE_SERIES_REF: return "series *";
    case GRETL_TYPE_MATRIX_REF: return "matrix *";
    case GRETL_TYPE_STRING:     return "string";	
    }

    return "unknown";
}

static const char *arg_type_xml_string (int t)
{
    if (t == GRETL_TYPE_SCALAR_REF) {
	return "scalarref";
    } else if (t == GRETL_TYPE_SERIES_REF) {
	return "seriesref";
    } else if (t == GRETL_TYPE_MATRIX_REF) {
	return "matrixref";
    } else {
	return arg_type_string(t);
    }
}

/* FIXME updating and placement of these function */

static int arg_type_from_string (const char *s)
{
    if (!strncmp(s, "bool", 4)) return GRETL_TYPE_BOOL;
    if (!strcmp(s, "int"))      return GRETL_TYPE_INT;
    if (!strcmp(s, "scalar"))   return GRETL_TYPE_DOUBLE;
    if (!strcmp(s, "series"))   return GRETL_TYPE_SERIES;
    if (!strcmp(s, "list"))     return GRETL_TYPE_LIST;
    if (!strcmp(s, "matrix"))   return GRETL_TYPE_MATRIX;

    if (!strcmp(s, "scalar *"))  return GRETL_TYPE_SCALAR_REF;
    if (!strcmp(s, "series *"))  return GRETL_TYPE_SERIES_REF;
    if (!strcmp(s, "matrix *"))  return GRETL_TYPE_MATRIX_REF;

    if (!strcmp(s, "scalarref"))  return GRETL_TYPE_SCALAR_REF;
    if (!strcmp(s, "seriesref"))  return GRETL_TYPE_SERIES_REF;
    if (!strcmp(s, "matrixref"))  return GRETL_TYPE_MATRIX_REF;

    if (!strcmp(s, "string")) return GRETL_TYPE_STRING;

    return 0;
}

/* backward compatibility for nasty old numeric type references */

static int arg_type_from_int (const char *s)
{
    int i = atoi(s);

    switch (i) {
    case 1: return GRETL_TYPE_DOUBLE;
    case 2: return GRETL_TYPE_SERIES;
    case 3: return GRETL_TYPE_LIST;
    case 4: return GRETL_TYPE_MATRIX;
    case 5: return GRETL_TYPE_BOOL;
    case 6: return GRETL_TYPE_INT;
    case 7: return GRETL_TYPE_SCALAR_REF;
    case 8: return GRETL_TYPE_SERIES_REF;
    case 9: return GRETL_TYPE_MATRIX_REF;
    }

    return GRETL_TYPE_NONE;
}

static int field_to_type (const char *s)
{
#if FN_DEBUG
    fprintf(stderr, "field_to_type: looking at '%s'\n", s);
#endif

    if (isdigit(*s)) {
	return arg_type_from_int(s);
    } else {
	return arg_type_from_string(s);
    }
}    

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
		if (gretl_scalar_type(fun->params[n].type)) {
		    gretl_xml_get_prop_as_double(cur, "default", 
						 &fun->params[n].deflt);
		}
		if (fun->params[n].type == GRETL_TYPE_INT) {
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

static int func_read_code (xmlNodePtr node, xmlDocPtr doc, ufunc *fun)
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
	s = line;
	while (isspace(*s)) s++;
	tailstrip(s);
	err = strings_array_add(&fun->lines, &fun->n_lines, s);
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
	if (fun->params[i].type == GRETL_TYPE_BOOL) {
	    if (!na(fun->params[i].deflt)) {
		pprintf(prn, "[%g]", fun->params[i].deflt);
	    }
	} else if (gretl_scalar_type(fun->params[i].type)) {
	    print_deflt_min_max(&fun->params[i], prn);
	} else if (gretl_ref_type(fun->params[i].type) || 
		   fun->params[i].type == GRETL_TYPE_LIST) {
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

    if (fun->rettype != GRETL_TYPE_NONE) {
	typestr = arg_type_string(fun->rettype);
	pprintf(prn, "  return %s %s\n", typestr, fun->retname);
    }

    pputs(prn, "end function\n");
}

static int unload_package_by_ID (int ID)
{
    int i;

#if PKG_DEBUG
    fprintf(stderr, "unload_package_by_ID: unloading %d\n", ID);
#endif

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkgID == ID) {
	    delete_ufunc_from_list(ufuns[i]);
	}
    }

    function_package_remove_by_ID(ID);

    return 0;
}

static int attach_ufunc_to_package (ufunc *fun, fnpkg *pkg)
{
    int err = 0;

    if (fun->private) {
	ufunc **priv;

	priv = realloc(pkg->priv, (pkg->n_priv + 1) * sizeof *priv);
	if (priv == NULL) {
	    err = E_ALLOC;
	} else {
	    pkg->priv = priv;
	    pkg->priv[pkg->n_priv] = fun;
	    pkg->n_priv += 1;
	}
    } else if (pkg->iface != NULL) {
	err = E_DATA;
    } else {
	pkg->iface = fun;
    }

#if PKG_DEBUG
    fprintf(stderr, "attach_ufunc_to_package: id = %d, "
	    "private = %d, err = %d\n", fun->pkgID, 
	    fun->private, err);
#endif

    return err;
}

static int read_ufunc_from_xml (xmlNodePtr node, xmlDocPtr doc, fnpkg *pkg)
{
    ufunc *fun = ufunc_new();
    xmlNodePtr cur;
    char *tmp;
    int err = 0;

    if (fun == NULL) {
	return E_ALLOC;
    }

    if (!gretl_xml_get_prop_as_string(node, "name", &tmp)) {
	free_ufunc(fun);
	return E_DATA;
    }

    strncat(fun->name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (pkg != NULL) {
	fun->pkgID = pkg->ID;
    }

    gretl_xml_get_prop_as_int(node, "private", &fun->private);

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: name '%s', private = %d\n",
	    fun->name, fun->private);
#endif

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
	} else if (!xmlStrcmp(cur->name, (XUC) "code")) {
	    err = func_read_code(cur, doc, fun);
	}
	cur = cur->next;
    }

    if (!err) {
	if (pkg != NULL) {
	    err = attach_ufunc_to_package(fun, pkg);
	} else {
	    err = add_allocated_ufunc(fun);
	}
    }

    if (err) {
	free_ufunc(fun);
    } 

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: returning %d\n", err);
#endif

    return err;
}

static int wordmatch (const char *s, const char *test)
{
    int n = strlen(test);
    
    return (!strncmp(s, test, n) && (s[n] == '\0' || isspace(s[n])));
}

void adjust_indent (const char *s, int *this_indent, int *next_indent)
{
    int ti = *next_indent;
    int ni = *next_indent;

    if (wordmatch(s, "loop")) {
	ni++;
    } else if (wordmatch(s, "if")) {
	ni++;
    } else if (wordmatch(s, "nls")) {
	ni++;
    } else if (wordmatch(s, "mle")) {
	ni++;
    } else if (wordmatch(s, "gmm")) {
	ni++;
    } else if (wordmatch(s, "function")) {
	ni++;
    } else if (wordmatch(s, "restrict")) {
	ni++;
    } else if (wordmatch(s, "end") ||
	       wordmatch(s, "endif") ||
	       wordmatch(s, "endloop")) {
	ti--;
	ni--;
    } else if (wordmatch(s, "else") ||
	       wordmatch(s, "elif")) {
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
		    arg_type_xml_string(fun->params[i].type));
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

    if (fun->rettype != GRETL_TYPE_NONE) {
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

static int real_function_print_code (ufunc *fun, PRN *prn)
{
    int this_indent = 0;
    int next_indent = 0;
    int i, j;
   
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

int gretl_function_print_code (int i, PRN *prn)
{
    if (i < 0 || i >= n_ufuns) {
	return 1;
    } else {
	real_function_print_code(ufuns[i], prn);
	return 0;
    }
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

    p = strrchr(ret, '-');
    if (p == NULL) {
	p = strstr(ret, ".gfn");
    }

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

static fnpkg *new_package_with_funcs (const char *fname, int pub, 
				      const int *privlist)
{
    fnpkg *pkg = function_package_new(fname);
    int i;

    if (pkg == NULL) {
	return NULL;
    } 

    if (privlist != NULL && privlist[0] > 0) {
	pkg->priv = malloc(privlist[0] * sizeof *pkg->priv);
	if (pkg->priv == NULL) {
	    free(pkg);
	    return NULL;
	}
	for (i=1; i<=privlist[0]; i++) {
	    pkg->priv[i-1] = ufuns[privlist[i]];
	}
	pkg->n_priv = privlist[0];
    }

    pkg->iface = ufuns[pub];

    return pkg;
}

/* Three cases should be distinguished below: (a) saving a totally new
   package; (b) saving changes to an existing package with no change
   to the version; and (c) saving changes to an existing package with
   a version change.  

   In case (a) the "pkg" argument will be NULL, otherwise it will be
   non-NULL.  We can tell cases (b) and (c) apart by checking whether
   the given "fname" differs from the existing package's filename: if
   so, we're in case (c), which we mark using "saveas".
*/

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
    int newpkg = 0;
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

    if (pkg == NULL) {
	pkg = new_package_with_funcs(fname, pub, privlist);
	if (pkg == NULL) {
	    return E_ALLOC;
	} 
	newpkg = 1;
    } else if (strcmp(fname, pkg->fname)) {
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

    if (newpkg || saveas) {
	pkg->ID = (int) time(NULL);
    } 

    fprintf(fp, " ID=\"%d\"", pkg->ID);

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

    ufuns[pub]->pkgID = pkg->ID;
    write_function_xml(ufuns[pub], fp);

    if (privlist != NULL) {
	for (i=1; i<=privlist[0]; i++) {
	    fi = privlist[i];
	    if (fi >= 0 && fi < n_ufuns) {
		ufuns[fi]->pkgID = pkg->ID;
		write_function_xml(ufuns[fi], fp);
	    }
	}
    }

    fputs("</gretl-function-package>\n", fp);
    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    if (newpkg) {
	pkg->author = gretl_strdup(author);
	pkg->version = gretl_strdup(version);
	pkg->date = gretl_strdup(date);
	pkg->descrip = gretl_strdup(descrip);
    } else {
	if (strcmp(fname, pkg->fname)) {
	    free(pkg->fname);
	    pkg->fname = gretl_strdup(fname);
	}
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
    } 

    if (pkg->author == NULL || pkg->version == NULL ||
	pkg->date == NULL || pkg->descrip == NULL ||
	pkg->fname == NULL) {
	err = E_ALLOC;
    } else {
	pkg->dreq = dreq;
	pkg->minver = minver;
	if (newpkg) {
	    err = function_package_add(pkg);
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

#if PKG_DEBUG
    fprintf(stderr, "Found package, ID = %d\n", pkg->ID);
#endif

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
#if PKG_DEBUG
	fprintf(stderr, "checking func %i, pkgID = %d\n", i, 
		ufuns[i]->pkgID);
#endif
	if (ufuns[i]->pkgID == pkg->ID) {
	    if (ufuns[i]->private) {
		npriv++;
	    } else {
		pubnum = i;
	    }
	}
    }

#if PKG_DEBUG
    fprintf(stderr, "npriv = %d, pubnum = %d\n", npriv, pubnum);
#endif

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

enum {
    PKG_FREE_PKG,
    PKG_FREE_ALL
};

static void function_package_free (fnpkg *pkg, int mode)
{
    if (pkg != NULL) {
	int i;

	for (i=0; i<n_ufuns; i++) {
	    if (ufuns[i]->pkgID == pkg->ID) {
		ufuns[i]->pkgID = 0;
	    }
	}

	if (mode == PKG_FREE_ALL && pkg->iface != NULL) {
	    free_ufunc(pkg->iface);
	}

	if (pkg->priv != NULL) {
	    if (mode == PKG_FREE_ALL) {
		for (i=0; i<pkg->n_priv; i++) {
		    free_ufunc(pkg->priv[i]);
		}
	    }
	    free(pkg->priv);
	}

	free(pkg->fname);
	free(pkg->author);
	free(pkg->version);
	free(pkg->date);
	free(pkg->descrip);
	free(pkg);
    }
}

static int function_package_remove_by_ID (int ID)
{
    int i, j, err = E_DATA;

    for (i=0; i<n_pkgs; i++) {
	if (pkgs[i]->ID == ID) {
	    function_package_free(pkgs[i], PKG_FREE_PKG);
	    for (j=i; j<n_pkgs-1; j++) {
		pkgs[j] = pkgs[j+1];
	    }
	    err = 0;
	    break;
	}
    }

    if (err) {
	return err;
    }

    if (n_pkgs == 1) {
	free(pkgs);
	pkgs = NULL;
	n_pkgs = 0;
    } else {
	fnpkg **tmp;

	tmp = realloc(pkgs, (n_pkgs - 1) * sizeof *tmp);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    pkgs = tmp;
	    n_pkgs--;
	}
    }

    return err;
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

static int pkg_name_from_filename (fnpkg *pkg)
{
    const char *p = strrchr(pkg->fname, SLASH);
    int n, err = 0;

    if (p != NULL) {
	p++;
    } else {
	p = pkg->fname;
    }

    n = strcspn(p, "-");

    if (n == strlen(p) && has_suffix(p, ".gfn")) {
	n -= 4;
    }

    if (n >= FN_NAMELEN) {
	fprintf(stderr, "pkg_name_from_filename: name is too long\n");
	n = FN_NAMELEN - 1;
	err = 1;
    }

    strncat(pkg->name, p, n);

#if PKG_DEBUG
    fprintf(stderr, "pkg_name_from_filename: calculated '%s'\n", pkg->name);
#endif

    return err;
}

static int fname_is_tmpfile (const char *fname)
{
    const char *p = strrchr(fname, SLASH);

    if (p == NULL) {
	p = fname;
    } else {
	p++;
    }

    return strncmp(p, "dltmp.", 6) == 0;
}

static fnpkg *function_package_new (const char *fname)
{
    fnpkg *pkg = malloc(sizeof *pkg);

    if (pkg == NULL) {
	return NULL;
    }

    pkg->ID = n_pkgs + 1; /* FIXME? */
    *pkg->name = '\0';
    pkg->author = NULL;
    pkg->version = NULL;
    pkg->date = NULL;
    pkg->descrip = NULL;
    pkg->dreq = 0;
    pkg->minver = 0;

    pkg->iface = NULL;
    pkg->priv = NULL;
    pkg->n_priv = 0;

    pkg->fname = gretl_strdup(fname);
    if (pkg->fname == NULL) {
	free(pkg);
	pkg = NULL;
    } else if (!fname_is_tmpfile(fname)) {
	pkg_name_from_filename(pkg);
    }

    return pkg;
}

static void packages_destroy (void)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	function_package_free(pkgs[i], PKG_FREE_PKG);
    }
    free(pkgs);

    pkgs = NULL;
    n_pkgs = 0;
}

static void maybe_clear_out_duplicate (ufunc *fun)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fun->name, ufuns[i]->name)) {
	    if (fun->pkgID == ufuns[i]->pkgID) {
		fprintf(stderr, "Redefining function '%s'\n", fun->name);
		delete_ufunc_from_list(ufuns[i]);
	    } else if (!fun->private && !ufuns[i]->private) {
		fprintf(stderr, "Redefining function '%s'\n", fun->name);
		if (ufuns[i]->pkgID == 0) {
		    delete_ufunc_from_list(ufuns[i]);
		} else {
		    unload_package_by_ID(ufuns[i]->pkgID);
		}
	    } 
	}
    }
}

static int real_load_package (fnpkg *pkg)
{
    int i, err;

#if PKG_DEBUG
    fprintf(stderr, "real_load_package: name='%s', ID = %d\n",
	    pkg->name, pkg->ID);
#endif

    err = function_package_add(pkg);

    if (!err && pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv && !err; i++) {
	    maybe_clear_out_duplicate(pkg->priv[i]);
	    err = add_allocated_ufunc(pkg->priv[i]);
	}
    }

    if (!err && pkg->iface != NULL) {
	maybe_clear_out_duplicate(pkg->iface);
	err = add_allocated_ufunc(pkg->iface);
    }

    return err;
}

static float version_float_from_string (const char *s)
{
    int maj, min, pl;

    sscanf(s, "%d.%d.%d", &maj, &min, &pl);
    return maj + min / 10.0 + pl / 100.0;
}

static void print_package_info (const fnpkg *pkg, PRN *prn)
{
    pprintf(prn, "Package: %s\n", (*pkg->name)? pkg->name : "unknown");
    pprintf(prn, "Author: %s\n", (pkg->author)? pkg->author : "unknown");
    pprintf(prn, "Version: %s\n", (pkg->version)? pkg->version : "unknown");
    pprintf(prn, "Date: %s\n", (pkg->date)? pkg->date : "unknown");
    pprintf(prn, "Description: %s\n", (pkg->descrip)? pkg->descrip : "none");

    pputc(prn, '\n');
    real_user_function_help(pkg->iface, 0, prn);
}

static void print_package_code (const fnpkg *pkg, PRN *prn)
{
    int i;

    if (pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv; i++) {
	    real_function_print_code(pkg->priv[i], prn);
	    pputc(prn, '\n');
	}
    }

    if (pkg->iface != NULL) {
	real_function_print_code(pkg->iface, prn);
    }
}

static void check_package_name (const char *correct, char *pname)
{
    if (pname != NULL && strcmp(correct, pname)) {
	fprintf(stderr, "Got old package name '%s', should be '%s'\n",
		pname, correct);
    }
}

/* allocate a fnpkg structure and read from XML file into it */

static fnpkg * 
real_read_package (xmlDocPtr doc, xmlNodePtr node, const char *fname, int *err)
{
    xmlNodePtr cur;
    fnpkg *pkg;
    char *tmp = NULL;
    int id;

#if PKG_DEBUG
    fprintf(stderr, "real_read_package: fname='%s'\n", fname);
#endif

    pkg = function_package_new(fname);
    if (pkg == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    gretl_xml_get_prop_as_string(node, "name", &tmp);

    if (tmp == NULL) {
	*err = E_DATA;
	function_package_free(pkg, PKG_FREE_PKG);
	return NULL;
    }

    if (fname_is_tmpfile(fname)) {
	strncat(pkg->name, tmp, FN_NAMELEN - 1);
    } else {
	check_package_name(pkg->name, tmp);
    }

    free(tmp);

    if (gretl_xml_get_prop_as_bool(node, NEEDS_TS)) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_QM)) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_PANEL)) {
	pkg->dreq = FN_NEEDS_PANEL;
    }

    if (gretl_xml_get_prop_as_string(node, "minver", &tmp)) {
	pkg->minver = version_float_from_string(tmp);
	free(tmp);
    }

    if (gretl_xml_get_prop_as_int(node, "ID", &id)) {
	pkg->ID = id;
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

    cur = node->xmlChildrenNode;
    while (cur != NULL && !*err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
	    *err = read_ufunc_from_xml(cur, doc, pkg);
	}
	cur = cur->next;
    }

#if PKG_DEBUG
    fprintf(stderr, "real_read_package: err = %d\n", *err);
#endif

    return pkg;
}

/* read aggregated file that may contain one or more function packages
   and one or more "loose" functions */

int read_session_functions_file (const char *fname)
{
    fnpkg *pkg = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    int err = 0;

#if PKG_DEBUG
    fprintf(stderr, "read_session_functions_file: starting on '%s'\n", fname);
#endif

    xmlKeepBlanksDefault(0);

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return err;
    }

    /* first get any function packages from this file */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
	    pkg = real_read_package(doc, cur, fname, &err);
	    if (!err) {
		err = real_load_package(pkg);
	    }
	} 
	cur = cur->next;
    }

    /* then get any unpackaged functions */
    if (!err) {
	cur = node->xmlChildrenNode;
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
		err = read_ufunc_from_xml(cur, doc, NULL);
	    }
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

#if PKG_DEBUG
    fprintf(stderr, "read_session_functions_file: returning %d\n", err);
#endif

    return err;
}

/* parse an XML function package file and return an allocated
   package struct with the functions attached */

static fnpkg *read_package_file (const char *fname, int *err)
{
    fnpkg *pkg = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;

#if PKG_DEBUG
    fprintf(stderr, "read_function_package: starting on '%s'\n", fname);
#endif

    xmlKeepBlanksDefault(0);

    *err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (*err) {
	return NULL;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
	    pkg = real_read_package(doc, cur, fname, err);
	    break;
	} 
	cur = cur->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

#if PKG_DEBUG
    fprintf(stderr, "read_function_package: err = %d\n", *err);
#endif

    return pkg;
}

int function_package_is_loaded (const char *fname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    return 1;
	}
    }

    return 0;
}

const char *function_package_description (const char *fname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    return pkgs[i]->descrip;
	}
    }

    return NULL;
}

static fnpkg *get_loaded_package (const char *fname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    return pkgs[i];
	}
    }

    return NULL;
}

/* read functions from file into gretl's workspace */

int load_user_function_file (const char *fname)
{
    fnpkg *pkg = NULL;
    int err = 0;

    if (get_loaded_package(fname)) {
	/* no-op */
	return 0;
    }

    pkg = read_package_file(fname, &err);
    if (err) {
	return err;
    } 

    return real_load_package(pkg);
}

/* Retrieve summary info or code listing for a function package,
   identified by its filename.  If the package is loaded in memory, we
   read from the package file.
*/

static int 
real_get_function_file_info (const char *fname, PRN *prn, char **pname,
			     int task)
{
    fnpkg *pkg = NULL;
    int free_pkg = 0;
    int err = 0;

    pkg = get_loaded_package(fname);

#if PKG_DEBUG
    fprintf(stderr, "real_get_function_file_info: get_loaded_package gave %p\n",
	    (void *) pkg);
#endif

    if (pkg == NULL) {
	pkg = read_package_file(fname, &err);
	free_pkg = 1;
    }

    if (!err) {
	*pname = gretl_strdup(pkg->name);
	if (task == FUNCS_INFO) {
	    print_package_info(pkg, prn);
	} else {
	    print_package_code(pkg, prn);
	}
	if (free_pkg) {
	    function_package_free(pkg, PKG_FREE_ALL);
	}
    }

    return err;
}

int get_function_file_info (const char *fname, PRN *prn, char **pname)
{
    return real_get_function_file_info(fname, prn, pname, FUNCS_INFO);
}

int get_function_file_code (const char *fname, PRN *prn, char **pname)
{
    return real_get_function_file_info(fname, prn, pname, FUNCS_CODE);
}

/* Just read the header from a function package file -- this is
   used when displaying the available packages.  We write the
   version number (as a string) into *pver.
*/

char *get_function_file_header (const char *fname, char **pver,
				int *err)
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
		} else if (!xmlStrcmp(sub->name, (XUC) "version")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, pver);
		}
		if (descrip != NULL && *pver != NULL) {
		    break;
		}
		sub = sub->next;
	    }
	    if (descrip != NULL && *pver != NULL) {
		break;
	    }
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
    if (*pver == NULL) {
	*pver = gretl_strdup("unknown");
    }

    if (descrip == NULL || *pver == NULL) {
	*err = 1;
    }

    return descrip;
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

static int maybe_delete_function (const char *fname, PRN *prn)
{
    ufunc *fun = get_user_function_by_name(fname);
    int err = 0;

    if (fun == NULL) {
	; /* no-op */
    } else if (fun->in_use != NULL) {
	sprintf(gretl_errmsg, "%s: function is in use", fname);
	err = 1;
    } else if (fun->pkgID != 0) {
	sprintf(gretl_errmsg, "%s: function belongs to package", fname);
	err = 1;
    } else {
	delete_ufunc_from_list(fun);
	if (gretl_messages_on()) {
	    pprintf(prn, _("Deleted function '%s'\n"), fname);
	}
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

	if (param->type == GRETL_TYPE_BOOL) {
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

    if (gretl_scalar_type(type)) {
	param->type = type;
	err = read_deflt_min_max(s, param, &len);
    }

    if (gretl_ref_type(type) || type == GRETL_TYPE_LIST) {
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
    char *p, *s = NULL;
    int i, len, np = 0;
    int err = 0;

    while (isspace(*str)) str++;

    /* get the length of the function name */
    len = strcspn(str, " (");
    if (len == 0) {
	return E_PARSE;
    } 

    if (*fname == '\0') {
	char fmt[8] = "%31s";

	if (len < FN_NAMELEN - 1) {
	    sprintf(fmt, "%%%ds", len);
	}
	if (sscanf(str, fmt, fname) != 1) {
	    err = E_PARSE;
	}
	if (!err) {
	    err = check_func_name(fname, pfun, prn);
	}
    }

    if (err) {
	return err;
    }

    str += len;

    if (*str == '\0') {
	/* void function */
	return 0;
    }

    /* move to next bit and make a copy */
    str += strspn(str, " (");
    if (*str == 0) {
	err = E_PARSE;
    } else {
	s = gretl_strdup(str);
	if (s == NULL) {
	    return E_ALLOC;
	}
    }

    if (!strcmp(s, ")")) {
	/* void function "foo()" */
	free(s);
	return 0;
    }

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
	    return maybe_delete_function(fname, prn);
	}
    } 

    /* the following takes care of replacing an existing function
       of the same name, if any */

    *fname = '\0';
    err = parse_fn_definition(fname, &params, &n_params,
			      line + 8, &fun, prn);

    if (!err && fun == NULL) {
	fun = add_ufunc(fname);
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

    if (fun->rettype != GRETL_TYPE_NONE) {
	sprintf(gretl_errmsg, "Function %s: return value is already defined",
		fun->name);
	return 1;
    }

    if (sscanf(line, "%15s %15s", s1, s2) != 2) {
	return E_PARSE;
    }

    type = field_to_type(s1);

#if FN_DEBUG
    fprintf(stderr, "add_function_return: s1='%s', s2='%s'\n", s1, s2);
    fprintf(stderr, "field_to_type on '%s' gives %d\n", s1, type);
#endif

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

static int real_function_append_line (const char *line, ufunc *fun)
{
    int editing = 1;
    int err = 0;

    if (fun == NULL) {
	fun = current_ufun;
	editing = 0;
    } 

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

    if (end_of_function(line)) {
	if (fun->n_lines == 0) {
	    sprintf(gretl_errmsg, "%s: empty function", fun->name);
	    err = 1;
	}
	set_compiling_off();
    } else if (!strncmp(line, "quit", 4)) {
	/* abort compilation */
	if (!editing) {
	    delete_ufunc_from_list(fun);
	}
	set_compiling_off();
    } else if (!strncmp(line, "function", 8)) {
	strcpy(gretl_errmsg, "You can't define a function within a function");
	err = 1;
    } else if (!strncmp(line, "return ", 7)) {
	err = add_function_return(fun, line + 7);
    } else {  
	err = strings_array_add(&fun->lines, &fun->n_lines, line);
    }

    if (err && !editing) {
	delete_ufunc_from_list(fun);
	set_compiling_off();
    }	

    return err;
}

int gretl_function_append_line (const char *line)
{
    return real_function_append_line(line, NULL);
}

static int extract_funcname (char *name, const char *s)
{
    int n = strcspn(s, " (");

    if (n == 0 || n >= FN_NAMELEN) {
	return E_DATA;
    } else {
	*name = '\0';
	strncat(name, s, n);
	return 0;
    }
}

/* Called from GUI window within the package-editing apparatus: the
   content of the script-editing window is dumped to file and we're
   passed the filename, along with the index number of the function
   interface that is being edited.  We parse the new function
   definition, but hold off replacing the original definition until we
   know there are no compilation errors.
*/

int update_function_from_script (const char *fname, int idx)
{
    char line[MAXLINE];
    ufunc *fun, *orig;
    char *s;
    FILE *fp;
    int gotfn = 0;
    int err = 0;

    if (idx < 0 || idx >= n_ufuns) {
	return E_DATA;
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fun = ufunc_new();
    if (fun == NULL) {
	fclose(fp);
	return E_ALLOC;
    }

    orig = ufuns[idx];

    fprintf(stderr, "Going to update function id %d '%s' from %s\n",
	    idx, orig->name, fname);

    while (fgets(line, sizeof line, fp) && !err) {
	s = line;
	while (*s == ' ') s++;
	tailstrip(s);
	if (!strncmp(s, "function ", 9)) {
	    if (gotfn || extract_funcname(fun->name, s + 9)) {
		err = 1;
	    } else if (strcmp(fun->name, orig->name)) {
		err = 1;
		strcpy(gretl_errmsg, 
		       _("You can't change the name of a function here"));
	    } else {
		gotfn = 1;
		err = parse_fn_definition(fun->name, &fun->params,
					  &fun->n_params, s + 8, 
					  NULL, NULL);
		if (err) {
		    strcpy(gretl_errmsg, _("Error compiling function"));
		}
	    }
	} else {
	    err = real_function_append_line(s, fun);
	    if (err) {
		strcpy(gretl_errmsg, _("Error compiling function"));
	    }
	}
    }

    fclose(fp);

    if (!err) {
	/* actually replace function content */
	free_strings_array(orig->lines, orig->n_lines);
	orig->n_lines = fun->n_lines;
	orig->lines = fun->lines;
	fun->lines = NULL;

	free_params_array(orig->params, orig->n_params);
	orig->n_params = fun->n_params;
	orig->params = fun->params;
	fun->params = NULL;

	orig->rettype = fun->rettype;
	free(orig->retname);
	orig->retname = fun->retname;
	fun->retname = NULL;
    } else {
	/* just trash the attempt */
	free_strings_array(fun->lines, fun->n_lines);
	free_params_array(fun->params, fun->n_params);
	free(fun->retname);
    }

    free(fun);

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
    int level = fn_executing + 1;
    int i, err;

    list = get_list_by_name(oldname);

    if (list == NULL) {
	err = E_DATA;
    } else {
	err = copy_named_list_as(oldname, fp->name);
    }

#if UDEBUG
    fprintf(stderr, "localize_list: localizing '%s' as '%s'\n",
	    oldname, fp->name);
#endif

    if (!err) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] != 0) {
		STACK_LEVEL(pdinfo, list[i]) = level;
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

    if (param->type == GRETL_TYPE_BOOL || param->type == GRETL_TYPE_INT) {
	x = floor(param->deflt);
    } else {
	x = param->deflt;
    }
    
    return dataset_add_scalar_as(x, param->name, pZ, pdinfo);
}

static int allocate_function_args (ufunc *fun,
				   fnargs *args, 
				   double ***pZ,
				   DATAINFO *pdinfo)
{
    struct fnarg *arg;
    fn_param *fp;
    int i, err = 0;

    for (i=0; i<args->argc && !err; i++) {
	arg = args->arg[i];
	fp = &fun->params[i];

	if (gretl_scalar_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_UVAR) {
		err = dataset_copy_variable_as(arg->val.idnum, fp->name,
					       pZ, pdinfo);
	    } else if (arg->type == GRETL_TYPE_NONE) {
		err = add_scalar_arg_default(fp, pZ, pdinfo);
	    } else {
		err = dataset_add_scalar_as(arg->val.x, fp->name, 
					    pZ, pdinfo);
	    }
	} else if (fp->type == GRETL_TYPE_SERIES) {
	    if (arg->type == GRETL_TYPE_UVAR) {
		err = dataset_copy_variable_as(arg->val.idnum, fp->name,
					       pZ, pdinfo);
	    } else {
		err = dataset_add_series_as(arg->val.px, fp->name, 
					    pZ, pdinfo);
	    }	    
	} else if (fp->type == GRETL_TYPE_MATRIX) {
	    err = copy_matrix_as(arg->val.m, fp->name);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    if (arg->type == GRETL_TYPE_NONE) {
		err = create_named_null_list(fp->name);
	    } else {
		err = localize_list(arg->val.str, fp, pdinfo);
	    }
	} else if (fp->type == GRETL_TYPE_STRING) {
	    err = add_string_as(arg->val.str, fp->name);
	} else if (fp->type == GRETL_TYPE_SCALAR_REF ||
		   fp->type == GRETL_TYPE_SERIES_REF) {
	    int v = arg->val.idnum;

	    arg->upname = gretl_strdup(pdinfo->varname[v]);
	    if (arg->upname == NULL) {
		err = E_ALLOC;
	    } else {
		STACK_LEVEL(pdinfo, v) += 1;
		strcpy(pdinfo->varname[v], fp->name);
	    }
	} else if (fp->type == GRETL_TYPE_MATRIX_REF) {
	    user_matrix *u = arg->val.um;

	    arg->upname = gretl_strdup(user_matrix_get_name(u));
	    if (arg->upname == NULL) {
		err = E_ALLOC;
	    } else {	    
		user_matrix_adjust_level(u, 1);
		user_matrix_set_name(u, fp->name);
	    }
	}

	if (arg->type == GRETL_TYPE_UVAR && !err) {
	    arg->upname = gretl_strdup(pdinfo->varname[arg->val.idnum]);
	}
    }

    /* now for any parameters withut matching arguments */

    for (i=args->argc; i<fun->n_params && !err; i++) {
	fp = &fun->params[i];
	if (gretl_scalar_type(fp->type)) {
	    err = add_scalar_arg_default(fp, pZ, pdinfo);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    err = create_named_null_list(fp->name);
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

static int unlocalize_list (const char *lname, DATAINFO *pdinfo)
{
    const int *list = get_list_by_name(lname);
    int d = gretl_function_depth();
    int i, vi, err = 0;

#if UDEBUG
    fprintf(stderr, "unlocalize_list: '%s', function depth = %d\n", lname, d);
    printlist(list, lname);
#endif	    

    if (list == NULL) {
	err = E_DATA;
    } else {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi != 0 && vi < pdinfo->v) {
		if (STACK_LEVEL(pdinfo, vi) == d) {
		    STACK_LEVEL(pdinfo, vi) -= 1;
		}
		unset_var_const(pdinfo, vi);
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
#if UDEBUG
	fprintf(stderr, "get_list_return: calling unlocalize_list\n");
#endif
	*err = unlocalize_list(lname, pdinfo);
	if (!*err) {
	    *err = named_list_lower_level(lname);
	}
    }

    return ret;
}

static void 
maybe_set_return_description (ufunc *u, int rtype, DATAINFO *pdinfo, 
			      char **descrip)
{
    if (rtype == GRETL_TYPE_DOUBLE || rtype == GRETL_TYPE_SERIES) {
	int v = varindex(pdinfo, u->retname);

	if (v < pdinfo->v) {
	    *descrip = gretl_strdup(VARLABEL(pdinfo, v));
	}
    }
}

static int 
function_assign_returns (ufunc *u, fnargs *args, int rtype, 
			 double **Z, DATAINFO *pdinfo, 
			 void *ret, char **descrip, PRN *prn, 
			 int *perr)
{
    struct fnarg *arg;
    fn_param *fp;
    int i, err = 0;

#if UDEBUG
    fprintf(stderr, "function_assign_returns: retname = '%s', rtype = %d\n", 
	    u->retname, rtype);
#endif

    if (*perr == 0) {
	/* direct return value */
	if (rtype == GRETL_TYPE_DOUBLE) {
	    *(double *) ret = get_scalar_return(u->retname, Z, pdinfo, &err);
	} else if (rtype == GRETL_TYPE_SERIES) {
	    *(double **) ret = get_series_return(u->retname, Z, pdinfo, GET_COPY, &err);
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    *(gretl_matrix **) ret = get_matrix_return(u->retname, GET_COPY, &err);
	} else if (rtype == GRETL_TYPE_LIST) {
	    *(char **) ret = get_list_return(u->retname, pdinfo, &err);
	}

	if (err == E_UNKVAR) {
	    pprintf(prn, "Function %s did not provide the specified return value\n",
		    u->name);
	}

	*perr = err;

	if (!err && descrip != NULL) {
	    maybe_set_return_description(u, rtype, pdinfo, descrip);
	}
    }

    /* "indirect return" values: these should be restored even if the
       function bombed
    */

    for (i=0; i<args->argc; i++) {
	arg = args->arg[i];
	fp = &u->params[i];
	if (gretl_ref_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_SCALAR_REF ||
		arg->type == GRETL_TYPE_SERIES_REF) {
		int v = arg->val.idnum;

		STACK_LEVEL(pdinfo, v) -= 1;
		strcpy(pdinfo->varname[v], arg->upname);
	    } else if (arg->type == GRETL_TYPE_MATRIX_REF) {
		user_matrix *u = arg->val.um;

		user_matrix_adjust_level(u, -1);
		user_matrix_set_name(u, arg->upname);
	    }
	} else if (fp->type == GRETL_TYPE_LIST) {
#if UDEBUG
	    fprintf(stderr, "function_assign_returns: calling unlocalize_list on '%s'\n",
		    fp->name);
#endif	    
	    unlocalize_list(fp->name, pdinfo);
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

    u->args = NULL;

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

    anyerr = destroy_saved_strings_at_level(d);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_saved_strings_at_level(%d): err = %d\n", d, err);
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
		    anyerr = dataset_drop_variable(i--, pZ, pdinfo);
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
	    if (gretl_scalar_type(fp->type) || fp->type == GRETL_TYPE_SERIES) {
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

static void start_fncall (ufunc *u, fnargs *args)
{
    set_executing_on(u);
    u->args = args;
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

static void func_exec_callback (ExecState *s, double ***pZ,
				DATAINFO *pdinfo)
{
    int ci = s->cmd->ci;

    if (ci == VAR || ci == VECM) {
	maybe_stack_var(s->var, s->cmd);
    } else if (ci == END && !strcmp(s->cmd->param, "restrict")) {
	maybe_stack_var(s->var, s->cmd);
    }
}

static int uvar_scalar (struct fnarg *arg, const DATAINFO *pdinfo)
{
    if (arg->type == GRETL_TYPE_UVAR) {
	int v = arg->val.idnum;

	if (v >= 0 && v < pdinfo->v) {
	    return var_is_scalar(pdinfo, v);
	}
    }

    return 0;
}

static int uvar_series (struct fnarg *arg, const DATAINFO *pdinfo)
{
    if (arg->type == GRETL_TYPE_UVAR) {
	int v = arg->val.idnum;

	if (v >= 0 && v < pdinfo->v) {
	    return var_is_series(pdinfo, v);
	}
    }

    return 0;
}

static double arg_get_double_val (struct fnarg *arg,
				  const double **Z)
{
    if (gretl_scalar_type(arg->type)) {
	return arg->val.x;
    } else if (arg->type == GRETL_TYPE_UVAR) {
	return Z[arg->val.idnum][0];
    } else {
	return NADBL;
    }
}

static int check_function_args (ufunc *u, fnargs *args, 
				const double **Z,
				const DATAINFO *pdinfo,
				PRN *prn)
{
    struct fnarg *arg;
    fn_param *fp;
    double x;
    int i, err = 0;

    for (i=0; i<args->argc && !err; i++) {
	arg = args->arg[i];
	fp = &u->params[i];

	if ((fp->flags & ARG_OPTIONAL) && arg->type == GRETL_TYPE_NONE) {
	    ; /* this is OK */
	} else if (gretl_scalar_type(fp->type) && arg->type == GRETL_TYPE_DOUBLE) {
	    ; /* OK */
	} else if (gretl_scalar_type(fp->type) && uvar_scalar(arg, pdinfo)) {
	    ; /* OK */
	} else if (fp->type == GRETL_TYPE_SERIES && uvar_series(arg, pdinfo)) {
	    ; /* OK */
	} else if (fp->type != arg->type) {
	    pprintf(prn, "argv[%d] is of wrong type (got %s, should be %s)\n", 
		    i, arg_type_string(arg->type), arg_type_string(fp->type));
	    err = E_TYPES;
	}

	if (!err && fp->type == GRETL_TYPE_DOUBLE) {
	    x = arg_get_double_val(arg, Z);
	    if ((!na(fp->min) && x < fp->min) ||
		(!na(fp->max) && x > fp->max)) {
		pprintf(prn, "argv[%d]: scalar value %g out of bounds\n", i, x);
		err = E_INVARG;
	    }
	}
    }

    for (i=args->argc; i<u->n_params && !err; i++) {
	/* do we have defaults for any empty args? */
	fp = &u->params[i];
	if (!(fp->flags & ARG_OPTIONAL) && na(fp->deflt)) {
	    pprintf(prn, "%s: not enough arguments\n", u->name);
	    err = E_ARGS;
	}
    }

    return err;
}

static void fn_state_init (CMD *cmd, ExecState *state)
{
    cmd->list = NULL;
    cmd->param = NULL;
    cmd->extra = NULL;
    cmd->linfo = NULL;

    state->cmd = NULL;
    state->models = NULL;
}

int gretl_function_exec (ufunc *u, fnargs *args, int rtype,
			 double ***pZ, DATAINFO *pdinfo,
			 void *ret, char **descrip, 
			 PRN *prn)
{
    ExecState state;
    MODEL **models = NULL;
    char line[MAXLINE];
    CMD cmd;
    int orig_v = pdinfo->v;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int started = 0;
    int i, err = 0;
 
    *funcerr_msg = '\0';

#if FN_DEBUG
    fprintf(stderr, "gretl_function_exec: starting\n");
#endif

    err = maybe_check_function_needs(pdinfo, u);
    if (err) {
	return err;
    }

    /* precaution */
    fn_state_init(&cmd, &state);

#if FN_DEBUG
    fprintf(stderr, "gretl_function_exec: argc = %d\n", args->argc);
    fprintf(stderr, "u->n_params = %d\n", u->n_params);
#endif

    err = check_function_args(u, args, (const double **) *pZ, pdinfo, prn);

    if (!err) {
	err = allocate_function_args(u, args, pZ, pdinfo);
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
	state.callback = func_exec_callback;
    }

    if (!err) {
	start_fncall(u, args);
	started = 1;
    }

    /* get function lines in sequence and check, parse, execute */

    for (i=0; i<u->n_lines && !err; i++) {
	strcpy(line, u->lines[i]);
	err = maybe_exec_line(&state, pZ, &pdinfo);
	if (state.funcerr) {
#if UDEBUG
	    fprintf(stderr, "funcerr: gretl_function_exec: i=%d, line: '%s'\n", 
		    i, line);
#endif
	    pprintf(prn, "%s: %s\n", u->name, state.cmd->param);
	    set_funcerr_message(u, state.cmd->param);
	}
	if (gretl_execute_loop()) { 
#if UDEBUG
	    fprintf(stderr, "gretl_function_exec: calling gretl_loop_exec\n");
#endif
	    err = gretl_loop_exec(&state, pZ, &pdinfo);
	    if (state.funcerr) {
		pprintf(prn, "%s: %s\n", u->name, state.cmd->param);
		set_funcerr_message(u, state.cmd->param);
	    }
	    if (err) {
		fprintf(stderr, "function_exec: breaking on error in loop\n");
		break;
	    }
#if UDEBUG
	    fprintf(stderr, "gretl_function_exec: gretl_loop_exec done, err = %d\n", err);
#endif
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

    function_assign_returns(u, args, rtype, *pZ, pdinfo,
			    ret, descrip, prn, &err);

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

/* look up name of supplied argument based on name of variable
   inside function */

char *gretl_func_get_arg_name (const char *argvar)
{
    ufunc *u = currently_called_function();
    char *ret = NULL;

    if (u != NULL && u->args != NULL) {
	int i, n = u->args->argc;

	for (i=0; i<n; i++) {
	    if (!strcmp(argvar, u->params[i].name)) {
		if (u->args->arg[i]->upname != NULL) {
		    ret = gretl_strdup(u->args->arg[i]->upname); 
		}
		break;
	    }
	}
    }

    if (ret == NULL) {
	ret = gretl_strdup("");
    }

    return ret;
}

void gretl_functions_cleanup (void)
{
    ufuncs_destroy();
    packages_destroy();
}

static void real_user_function_help (ufunc *fun, int ci, PRN *prn)
{
    fnpkg *pkg = NULL;
    int i;

    if (ci == HELP) {
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

    if (fun->rettype != GRETL_TYPE_NONE) {
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
	real_user_function_help(fun, HELP, prn);
    }

    return err;
}
