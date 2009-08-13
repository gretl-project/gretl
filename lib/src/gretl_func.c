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
#include "gretl_scalar.h"
#include "flow_control.h"

#include <glib.h>

#define FNPARSE_DEBUG 0 /* debug parsing of function code */
#define EXEC_DEBUG 0    /* debugging of function execution */
#define UDEBUG 0        /* debug handling of args */
#define PKG_DEBUG 0     /* debug handling of function packages */
#define FN_DEBUG 0      /* miscellaneous debugging */
#define DDEBUG 0        /* debug the debugger */

typedef struct obsinfo_ obsinfo;

struct obsinfo_ {
    int structure;
    int pd;
    int t1, t2;
    char changed;
    char stobs[OBSLEN];
};

/* structure representing a parameter of a user-defined function */

typedef struct fn_param_ fn_param;

struct fn_param_ {
    char *name;     /* the name of the parameter */
    char *descrip;  /* its description */
    char type;      /* its type */
    char flags;     /* additional information (e.g. "const" flag) */
    double deflt;   /* default value */
    double min;     /* minimum value (scalar parameters only) */
    double max;     /* maximum value (scalar parameters only) */
};

/* structure representing a call to a user-defined function */

typedef struct fncall_ fncall;

struct fncall_ {
    ufunc *fun;    /* the function called */
    fnargs *args;  /* argument array */
    int *ptrvars;
    int *listvars;
    obsinfo obs;
};

/* structure representing a user-defined function */

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
    char *retname; /* used only with "old-style" function syntax */
    int debug;
};

/* structure representing a function "package" */

struct fnpkg_ {
    int ID;
    char name[FN_NAMELEN];
    char *fname;
    char *author;
    char *version;
    char *date;
    char *descrip;
    char *sample; /* sample caller script */
    int minver;
    FuncDataReq dreq;
    ufunc *iface;
    ufunc **priv;
    int n_priv;
};

enum {
    ARG_OPTIONAL = 1 << 0,
    ARG_CONST    = 1 << 1
};

/* structure representing an argument to a user-defined function */

struct fnarg {
    char type;
    char flags;       /* ARG_OPTIONAL, ARG_CONST as appropriate */
    const char *name; /* name as function param */
    char *upname;     /* name of supplied arg at caller level */
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

static int n_ufuns;
static ufunc **ufuns;
static ufunc *current_ufun;
static GList *callstack;

static int n_pkgs;
static fnpkg **pkgs;

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
    arg->flags = 0;
    arg->name = NULL;
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
	       type == GRETL_TYPE_USCALAR ||
	       type == GRETL_TYPE_USERIES) {
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

static fncall *fncall_new (ufunc *fun)
{
    fncall *call = malloc(sizeof *call);

    if (call != NULL) {
	call->fun = fun;
	call->ptrvars = NULL;
	call->listvars = NULL;
    }

    return call;
}

static void fncall_free (fncall *call)
{
    if (call != NULL) {
	free(call->ptrvars);
	free(call->listvars);
	free(call);
    }
}

static void set_listargs_from_call (fncall *call, DATAINFO *pdinfo)
{
    int i, v;

    if (pdinfo == NULL) {
	return;
    }

    for (i=1; i<pdinfo->v; i++) {
	unset_var_listarg(pdinfo, i);
    }

    if (call != NULL && call->listvars != NULL) {
	for (i=1; i<=call->listvars[0]; i++) {
	    v = call->listvars[i];
#if UDEBUG
	    fprintf(stderr, "setting listarg status on var %d (%s)\n",
		    v, pdinfo->varname[v]);
#endif
	    set_var_listarg(pdinfo, v);
	}
    }
}

static void set_executing_off (fncall *call, DATAINFO *pdinfo)
{
    fncall *popcall = NULL;

    fn_executing--;

    callstack = g_list_remove(callstack, call);
#if EXEC_DEBUG
    fprintf(stderr, "set_executing_off: removing call to %s, depth now %d\n",
	    call->fun->name, g_list_length(callstack));
#endif
    fncall_free(call);

    if (fn_executing > 0) {
	GList *tmp = g_list_last(callstack);
	
	popcall = tmp->data;
    } else {
	g_list_free(callstack);
	callstack = NULL;
    }

    if (pdinfo != NULL) {
	set_listargs_from_call(popcall, pdinfo);
    }
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

const char *fn_param_descrip (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NULL :
	fun->params[i].descrip;
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

    return ((gretl_ref_type(t) || 
	     t == GRETL_TYPE_LIST ||
	     t == GRETL_TYPE_STRING) && 
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

static fncall *current_function_call (void)
{
    if (callstack != NULL) {
	GList *tmp = g_list_last(callstack);

	return tmp->data;
    } else {
	return NULL;
    }
}

static ufunc *currently_called_function (void)
{
    fncall *call = current_function_call();

    if (call != NULL) {
	return call->fun;
    } else {
	return NULL;
    }
}

static int function_in_use (ufunc *fun)
{
    GList *tmp = callstack;
    fncall *call;

    while (tmp != NULL) {
	call = tmp->data;
	if (fun == call->fun) {
	    return 1;
	}
	tmp = tmp->next;
    }

    return 0;
}

int current_func_pkgID (void)
{
    ufunc *f = currently_called_function();

    return (f == NULL)? 0 : f->pkgID;
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
	if (!strcmp(name, ufuns[i]->name)) {
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
	    if (!strcmp(name, ufuns[i]->name)) {
		fun = ufuns[i];
		if (fun->pkgID == 0) {
		    break;
		} else {
		    fun = NULL;
		}
	    }
	}
    }	

#if FN_DEBUG > 1
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
	params[i].descrip = NULL;
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

    fun->debug = 0;

    return fun;
}

static void free_params_array (fn_param *params, int n)
{
    int i;

    if (params == NULL) return;

    for (i=0; i<n; i++) {
	free(params[i].name);
	free(params[i].descrip);
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
    case GRETL_TYPE_MATRIX:     return "matrix";	
    case GRETL_TYPE_LIST:       return "list";
    case GRETL_TYPE_SCALAR_REF: return "scalar *";
    case GRETL_TYPE_SERIES_REF: return "series *";
    case GRETL_TYPE_MATRIX_REF: return "matrix *";
    case GRETL_TYPE_STRING:     return "string";
    case GRETL_TYPE_VOID:       return "void";
    case GRETL_TYPE_NONE:       return "null";
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
    if (!strcmp(s, "matrix"))   return GRETL_TYPE_MATRIX;
    if (!strcmp(s, "list"))     return GRETL_TYPE_LIST;
    if (!strcmp(s, "string")) return GRETL_TYPE_STRING;

    if (!strcmp(s, "scalar *"))  return GRETL_TYPE_SCALAR_REF;
    if (!strcmp(s, "series *"))  return GRETL_TYPE_SERIES_REF;
    if (!strcmp(s, "matrix *"))  return GRETL_TYPE_MATRIX_REF;

    return 0;
}

#define ok_return_type(r) (r == GRETL_TYPE_DOUBLE || \
	                   r == GRETL_TYPE_SERIES || \
	                   r == GRETL_TYPE_MATRIX || \
	                   r == GRETL_TYPE_LIST || \
                           r == GRETL_TYPE_STRING || \
                           r == GRETL_TYPE_VOID)

static int return_type_from_string (const char *s)
{
    int t;

    if (!strcmp(s, "void")) {
	/* not OK as arg type, but OK as return */
	t = GRETL_TYPE_VOID;
    } else {
	t = arg_type_from_string(s);
    }

    return (ok_return_type(t))? t : 0;
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

static int param_field_to_type (const char *s)
{
#if FNPARSE_DEBUG
    fprintf(stderr, "param_field_to_type: looking at '%s'\n", s);
#endif

    if (isdigit(*s)) {
	return arg_type_from_int(s);
    } else {
	return arg_type_from_string(s);
    }
}    

static int func_read_params (xmlNodePtr node, xmlDocPtr doc,
			     ufunc *fun)
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
		break;
	    }
	    if (gretl_xml_get_prop_as_string(cur, "type", &field)) {
		fun->params[n].type = param_field_to_type(field);
		free(field);
		if (gretl_scalar_type(fun->params[n].type)) {
		    gretl_xml_get_prop_as_double(cur, "default", 
						 &fun->params[n].deflt);
		    if (fun->params[n].type != GRETL_TYPE_BOOL) {
			gretl_xml_get_prop_as_double(cur, "min", 
						     &fun->params[n].min);
			gretl_xml_get_prop_as_double(cur, "max", 
						     &fun->params[n].max);
		    }
		}
		if (gretl_xml_get_prop_as_bool(cur, "optional")) {
		    fun->params[n].flags |= ARG_OPTIONAL;
		}
		if (gretl_xml_get_prop_as_bool(cur, "const")) {
		    fun->params[n].flags |= ARG_CONST;
		}
	    } else {
		err = E_DATA;
		break;
	    }
	    gretl_xml_child_get_string(cur, doc, "description", 
				       &fun->params[n].descrip);
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

/* used only with old-style function definitions */

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
	fun->rettype = param_field_to_type(field);
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

static void print_param_description (fn_param *param, PRN *prn)
{
    if (param->descrip != NULL && *param->descrip != '\0') {
	pprintf(prn, " \"%s\"", param->descrip);
    }
}

static void print_min_max_deflt (fn_param *param, PRN *prn)
{
    if (na(param->min) && na(param->max) && na(param->deflt)) {
	return; /* no-op */
    } else if (na(param->min) && na(param->max)) {
	/* default value only */
	pprintf(prn, "[%g]", param->deflt);
	return;
    }

    pputc(prn, '[');
    if (!na(param->min)) pprintf(prn, "%g", param->min);
    pputc(prn, ':');
    if (!na(param->max)) pprintf(prn, "%g", param->max);
    pputc(prn, ':');
    if (!na(param->deflt)) pprintf(prn, "%g", param->deflt);
    pputc(prn, ']');
}

static void print_function_start (ufunc *fun, PRN *prn)
{
    const char *s;
    int i;

    if (fun->rettype == GRETL_TYPE_NONE) {
	pprintf(prn, "function void %s ", fun->name);
    } else {
	const char *typestr = arg_type_string(fun->rettype);

	pprintf(prn, "function %s %s ", typestr, fun->name);
    }

    gretl_push_c_numeric_locale();

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
	    print_min_max_deflt(&fun->params[i], prn);
	} else if (gretl_ref_type(fun->params[i].type) || 
		   fun->params[i].type == GRETL_TYPE_LIST ||
		   fun->params[i].type == GRETL_TYPE_STRING) {
	    print_opt_flags(&fun->params[i], prn);
	}
	print_param_description(&fun->params[i], prn);
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
    if (fun->retname != NULL) {
	/* support old-style */
	pprintf(prn, "  return %s\n", fun->retname);
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
	/* FIXME allow multiple interfaces */
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

/* for now we'll support old-style function definitions on
   reading from gfn XML */

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

    if (gretl_xml_get_prop_as_string(node, "type", &tmp)) {
	fun->rettype = return_type_from_string(tmp);
	free(tmp);
    } else {
	fun->rettype = GRETL_TYPE_VOID;
    }

    gretl_xml_get_prop_as_int(node, "private", &fun->private);

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: name '%s', type %d, private = %d\n",
	    fun->name, fun->rettype, fun->private);
#endif

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "help")) {
	    gretl_xml_node_get_string(cur, doc, &fun->help);
	} else if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    err = func_read_params(cur, doc, fun);
	    if (err) {
		fprintf(stderr, "%s: error parsing function parameters\n",
			fun->name);
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "return")) {
	    /* support old-style functions */
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
    } else if (wordmatch(s, "system")) {
	ni++;
    } else if (wordmatch(s, "foreign")) {
	ni++;
    } else if (wordmatch(s, "kalman")) {
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

/* ensure use of canonical forms "endif", "endloop" */

static void maybe_correct_line (char *line)
{
    char *p = strstr(line, "end if");

    if (p == NULL) {
	p = strstr(line, "end loop");
    }

    if (p != NULL && (p == line || *(p-1) == ' ')) {
	shift_string_left(p + 4, 1);
    }
}

/* on writing functions as XML, use new-style function syntax */

static int write_function_xml (const ufunc *fun, FILE *fp)
{
    int rtype = fun->rettype;
    int this_indent = 0;
    int next_indent = 0;
    int i, j;

    if (rtype == GRETL_TYPE_NONE) {
	rtype = GRETL_TYPE_VOID;
    }

    fprintf(fp, "<gretl-function name=\"%s\" type=\"%s\" private=\"%d\">\n", 
	    fun->name, arg_type_string(rtype), fun->private);

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
	    if (fun->params[i].descrip != NULL) {
		fputs(">\n", fp);
		gretl_xml_put_tagged_string("description", 
					    fun->params[i].descrip,
					    fp);
		fputs("  </param>\n", fp);
	    } else {
		fputs("/>\n", fp);
	    }
	}
	fputs(" </params>\n", fp);

	gretl_pop_c_numeric_locale();
    }

    fputs("<code>", fp);

    for (i=0; i<fun->n_lines; i++) {
	adjust_indent(fun->lines[i], &this_indent, &next_indent);
	for (j=0; j<this_indent; j++) {
	    fputs("  ", fp);
	}
	maybe_correct_line(fun->lines[i]);
	gretl_xml_put_raw_string(fun->lines[i], fp);
	fputc('\n', fp);
    }

    if (fun->rettype != GRETL_TYPE_NONE && fun->retname != NULL) {
	/* old-style */
	fprintf(fp, "  return %s\n", fun->retname);
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
	} else if (!strcmp(key, "sample")) {
	    *value = pkg->sample;
	} else if (!strcmp(key, "pkgname")) {
	    *value = pkg->name;
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

static char *get_version_string (char *s, int v)
{
    int x, y, z;

    x = v / 10000;
    y = (v - x * 10000) / 100;
    z = v % 100;

    sprintf(s, "%d.%d.%d", x, y, z);
    return s;
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

static char *trim_script (char *s)
{
    while (isspace(*s)) s++;
    return tailstrip(s);
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
			    char *sample,
			    FuncDataReq dreq,
			    int minver)
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
	char vstr[8];

	fprintf(fp, " minver=\"%s\"", get_version_string(vstr, minver));
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

    if (sample != NULL) {
 	/* escapes may be needed; also perhaps trimming */
	char *s = trim_script(sample);

	fputs("<sample-script>\n", fp);
	gretl_xml_put_raw_string(s, fp);
	fputs("\n</sample-script>\n", fp);	
    }

    fputs("</gretl-function-package>\n", fp);
    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    if (newpkg) {
	pkg->author = gretl_strdup(author);
	pkg->version = gretl_strdup(version);
	pkg->date = gretl_strdup(date);
	pkg->descrip = gretl_strdup(descrip);
	pkg->sample = gretl_strdup(sample);
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
	if (sample != NULL && 
	    (pkg->sample == NULL || strcmp(sample, pkg->sample))) {
	    free(pkg->sample);
	    pkg->sample = gretl_strdup(sample);
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
			       char **sample,
			       FuncDataReq *dreq,
			       int *minver)
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
    if (sample != NULL) {
	*sample = gretl_strdup(pkg->sample);
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

    if (pkg->sample != NULL) {
	/* escapes may be needed; also perhaps trimming */
	char *s = trim_script(pkg->sample);

	fputs("<sample-script>\n", fp);
	gretl_xml_put_raw_string(s, fp);
	fputs("\n</sample-script>\n", fp);	
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
	free(pkg->sample);
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
    pkg->sample = NULL;
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

static int version_number_from_string (const char *s)
{
    int x, y, z;

    sscanf(s, "%d.%d.%d", &x, &y, &z);
    return 10000 * x + 100 * y + z;
}

static void print_package_info (const fnpkg *pkg, PRN *prn)
{
    pprintf(prn, "Package: %s\n", (*pkg->name)? pkg->name : "unknown");
    pprintf(prn, "Author: %s\n", (pkg->author)? pkg->author : "unknown");
    pprintf(prn, "Version: %s\n", (pkg->version)? pkg->version : "unknown");
    pprintf(prn, "Date: %s\n", (pkg->date)? pkg->date : "unknown");
    pputs(prn, "Description: ");
    pputs(prn, (pkg->descrip)? pkg->descrip : "none");

    pputs(prn, "\n\n");
    real_user_function_help(pkg->iface, 0, prn);

    if (pkg->sample != NULL) {
	pputs(prn, "Sample script:\n");
	pputs(prn, pkg->sample);
	pputc(prn, '\n');
    }
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
	pkg->minver = version_number_from_string(tmp);
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
	} else if (!xmlStrcmp(cur->name, (XUC) "sample-script")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->sample);
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
	strcpy(gretl_errmsg, _("Function names must start with a letter"));
	err = 1;
    } else if (gretl_command_number(fname)) {
	sprintf(gretl_errmsg, _("'%s' is the name of a gretl command"),
		fname);
	err = 1;
    } else if (function_lookup(fname)) {
	sprintf(gretl_errmsg, _("'%s' is the name of a built-in function"),
		fname);
	err = 1;
    } else {
	for (i=0; i<n_ufuns; i++) {
	    if (!strcmp(fname, ufuns[i]->name)) {
		if (libset_get_bool(VERBOSE_INCLUDE)) {
		    pprintf(prn, _("Redefining function '%s'"), fname);
		    pputc(prn, '\n');
		}
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
    } else if (function_in_use(fun)) {
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
    int quoted = 0;
    int nc = 0;

    while (*s) {
	if (*s == '"') {
	    quoted = !quoted;
	} else if (!quoted && *s == ',') {
	    nc++;
	}
	s++;
    }

    return nc;
}

static int read_min_max_deflt (char **ps, fn_param *param)
{
    char *p = *ps;
    double x, y, z;
    int err = 0;

    if (param->type == GRETL_TYPE_BOOL) {
	if (sscanf(p, "[%lf]", &x) == 1) {
	    param->deflt = x;
	} else {
	    err = E_PARSE;
	}
    } else {
	if (sscanf(p, "[%lf:%lf:%lf]", &x, &y, &z) == 3) {
	    param->min = x;
	    param->max = y;
	    param->deflt = z;
	} else if (sscanf(p, "[%lf::%lf]", &x, &y) == 2) {
	    param->min = x;
	    param->deflt = y;
	} else if (sscanf(p, "[%lf:%lf:]", &x, &y) == 2) {
	    param->min = x;
	    param->max = y;
	} else if (sscanf(p, "[:%lf:%lf]", &x, &y) == 2) {
	    param->max = x;
	    param->deflt = y;
	} else if (sscanf(p, "[%lf]", &x) == 1) {
	    param->deflt = x;
	} else {
	    err = E_PARSE;
	}
    }

    if (!err) {
	p = strchr(p, ']');
	if (p == NULL) {
	    err = E_PARSE;
	} else {
	    *ps = p + 1;
	}
    }

    return err;
}

static int read_param_option (char **ps, fn_param *param)
{
    int err = E_PARSE;

#if FNPARSE_DEBUG
    fprintf(stderr, "read_param_option: got '%s'\n", *ps);
#endif

    if (!strncmp(*ps, "[null]", 6)) {
	param->flags |= ARG_OPTIONAL;
	err = 0;
	*ps += 6;
    }

    return err;
}

static int read_param_comment (char **ps, fn_param *param)
{
    char *p = *ps + 1;
    int len = 0;
    int err = E_PARSE;

    while (*p) {
	if (*p == '"') {
	    err = 0;
	    break;
	}
	len++;
	p++;
    }

    if (!err && len > 0) {
	p = *ps + 1;
	param->descrip = gretl_strndup(p, len);
	if (param->descrip == NULL) {
	    err = E_ALLOC;
	} else {
	    *ps = p + len + 1;
	}
    }

    return err;
}

static int parse_function_param (char *s, fn_param *param, int i)
{
    char tstr[22] = {0};
    char *name;
    int type, len;
    int err = 0;

#if FNPARSE_DEBUG
    fprintf(stderr, "parse_function_param: s = '%s'\n", s);
#endif

    while (isspace(*s)) s++;

    if (!strncmp(s, "const ", 6)) {
	param->flags |= ARG_CONST;
	s += 6;
	while (isspace(*s)) s++;
    }

    /* get parameter type */

    len = gretl_namechar_spn(s);
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
	if (*tstr == '\0') {
	    sprintf(gretl_errmsg, "Expected a type identifier");
	    err = E_PARSE;
	} else {
	    type = arg_type_from_string(tstr);
	    if (type == 0) {
		sprintf(gretl_errmsg, "Unrecognized data type '%s'", tstr);
		err = E_PARSE;
	    }
	} 
    }

    if (err) {
	return err;
    }
    
    while (isspace(*s)) s++;
    len = gretl_namechar_spn(s);
    if (len == 0) {
	sprintf(gretl_errmsg, "parameter %d: name is missing", i + 1);
	err = E_PARSE;
    } else {
	name = gretl_strndup(s, len);
	if (name == NULL) {
	    err = E_ALLOC;
	} else if (gretl_reserved_word(name)) {
	    free(name);
	    err = E_DATA;
	}
    }	

    if (err) {
	return err;
    }

    param->type = type;

    s += len;
    s += strspn(s, " ");

    if (gretl_scalar_type(type)) {
	if (*s == '[') { 
	    err = read_min_max_deflt(&s, param);
	}
    }

    if (gretl_ref_type(type) || 
	type == GRETL_TYPE_LIST ||
	type == GRETL_TYPE_STRING) {
	if (*s == '[') { 
	    err = read_param_option(&s, param);
	}
    } 

    s += strspn(s, " ");

    if (*s == '"') {
	err = read_param_comment(&s, param);
    } 

    if (!err && *s != '\0') {
	/* got trailing unparseable stuff */
	err = E_PARSE;
    }

    if (!err) {
	param->name = name;
    } else {
	free(name);
	free(param->descrip);
	param->descrip = NULL;
    }

#if FNPARSE_DEBUG
    if (!err) {
	fprintf(stderr, " param[%d] = '%s', ptype = %d\n", 
		i, name, type);
	fprintf(stderr, "  min=%g, max=%g, deflt=%g\n", 
		param->min, param->max, param->deflt);
	fprintf(stderr, "  comment = '%s'\n", param->descrip); 
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

/* Here we're parsing what follows "function ".  The
   traditional expectation is <funcname> (<args>),
   but in a new-style definition we can have
   <return-type> <funcname> (<args>)
*/

static int parse_fn_definition (char *fname, 
				fn_param **pparams,
				int *pnp,
				int *rettype,
				const char *str, 
				ufunc **pfun, 
				PRN *prn)
{
    fn_param *params = NULL;
    char *p, *s = NULL;
    int i, len, np = 0;
    int err = 0;

    /* skip to next word */
    while (isspace(*str)) str++;

    /* and find its length */
    len = gretl_namechar_spn(str);

    if (len == 0 || len >= FN_NAMELEN) {
	return E_PARSE;
    }

    if (len <= 8) {
	/* see if we have a type specifier: is so, record
	   it and move on to the next word
	*/
	char typeword[9] = {0};
	int t;

	strncat(typeword, str, len);
	t = return_type_from_string(typeword);

	if (t > 0) {
	    if (!ok_return_type(t)) {
		gretl_errmsg_set("Invalid return type for function");
		return E_TYPES;
	    } else {		
		*rettype = t;
		str += len;
		str += strspn(str, " ");
		len = gretl_namechar_spn(str);
		if (len == 0 || len >= FN_NAMELEN) {
		    return E_PARSE;
		}
	    }
	}
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
	int quoted = 0;

	p = s;
	while (*p) {
	    if (*p == '"') {
		quoted = !quoted;
	    } else if (!quoted && *p == ',') {
		*p = '\0';
	    }
	    p++;
	}
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
    int rettype = 0;
    char fname[FN_NAMELEN];
    char s1[FN_NAMELEN];
    char s2[FN_NAMELEN];
    int err = 0;

    nf = sscanf(line, "function %31s %31s", s1, s2);
    if (nf <= 0) {
	return E_PARSE;
    } 

    if (nf == 2) {
	if (!strcmp(s2, "clear") || !strcmp(s2, "delete")) {
	    return maybe_delete_function(s1, prn);
	}
    } 

    /* the following takes care of replacing an existing function
       of the same name, if any */

    *fname = '\0';
    err = parse_fn_definition(fname, &params, &n_params, &rettype,
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
	fun->rettype = rettype;
	current_ufun = fun;
	set_compiling_on();
    } else {
	current_ufun = NULL;
    }
    
    return err;
}

static int parse_function_return (ufunc *fun, const char *line)
{
    const char *s = line + 6; /* skip "return" */
    char s1[16], s2[VNAMELEN];
    int type, n;
    int err = 0;

#if FNPARSE_DEBUG
    fprintf(stderr, "parse_function_return: line = '%s'\n", line);
#endif

    n = sscanf(s, "%15s %15s", s1, s2);

    type = return_type_from_string(s1);

    if (type > 0 && fun->retname != NULL) {
	sprintf(gretl_errmsg, "%s: return value is already defined",
		fun->name);
	return E_PARSE;
    }

    if (fun->rettype == GRETL_TYPE_NONE) {
	/* old-style: return type should be specified inline */
	if (type == 0) {
	    gretl_errmsg_sprintf("%s: missing a valid return type\n", 
				 fun->name);
	    err = E_TYPES;
	} else {
	    err = check_varname(s2);
	    if (!err) {
		fun->retname = gretl_strdup(s2);
		if (fun->retname == NULL) {
		    err = E_ALLOC;
		} else {
		    fun->rettype = type;
		}
	    }
	}
    } else {
	/* return type was pre-defined: treat as normal line */
	err = strings_array_add(&fun->lines, &fun->n_lines, line);
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

#define function_return_line(s) (strncmp(s, "return", 6) == 0 && \
	                         (*(s + 6) == ' ' || *(s + 6) == '\0'))

#define bare_quote(p,s)   (*p == '"' && (p-s==0 || *(p-1) != '\\'))
#define starts_comment(p) (*p == '/' && *(p+1) == '*')
#define ends_comment(p)   (*p == '*' && *(p+1) == '/')

static int ignore_line (ufunc *fun)
{
    int i, quoted = 0, ignore = 0;
    char *s, *p;

    for (i=0; i<fun->n_lines; i++) {
	s = p = fun->lines[i];
	while (*p) {
	    if (!quoted && !ignore && *p == '#') {
		break;
	    }
	    if (!ignore && bare_quote(p, s)) {
		quoted = !quoted;
	    }
	    if (!quoted) {
		if (starts_comment(p)) {
		    ignore = 1;
		    p += 2;
		} else if (ends_comment(p)) {
		    ignore = 0;
		    p += 2;
		    p += strspn(p, " ");
		}
	    }
	    if (*p) {
		p++;
	    }
	}
    }

    return ignore;
}

#define NEEDS_IF(c) (c == ELSE || c == ELIF || c == ENDIF)

/* Rather minimal check for syntatic validity of "compiled" function.
   FIXME: it would be good to check here for messed up block structure
   too (e.g. "system" without "end system").
*/

static int check_function_structure (ufunc *fun)
{
    CMD cmd;
    int ifdepth = 0;
    int i, err = 0;

    gretl_cmd_init(&cmd);

    for (i=0; i<fun->n_lines && !err; i++) {
	get_command_index(fun->lines[i], &cmd);
	if (cmd.ci == FUNC) {
	    gretl_errmsg_set("You can't define a function within a function");
	    err = E_PARSE;
	} else if (cmd.ci == IF) {
	    ifdepth++;
	} else if (NEEDS_IF(cmd.ci) && ifdepth == 0) {
	    gretl_errmsg_sprintf("%s: unbalanced if/else/endif", fun->name);
	    err = E_PARSE;
	} else if (cmd.ci == ENDIF) {
	    ifdepth--;
	} 
    }

    if (!err && ifdepth != 0) {
	gretl_errmsg_sprintf("%s: unbalanced if/else/endif", fun->name);
	err = E_PARSE;
    }

    gretl_cmd_free(&cmd);

    return err;
}

static int real_function_append_line (const char *line, ufunc *fun)
{
    int editing = 1;
    int err = 0;

    if (fun == NULL) {
	fun = current_ufun;
	editing = 0;
    } 

#if FNPARSE_DEBUG
    fprintf(stderr, "gretl_function_append_line: '%s'\n", line);
#endif

    if (fun == NULL) {
#if FN_DEBUG
	fprintf(stderr, " fun == NULL!\n");
#endif
	return 1;
    }

    if (string_is_blank(line)) {
	err = strings_array_add(&fun->lines, &fun->n_lines, "");
    } else if (end_of_function(line) && !ignore_line(fun)) {
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
	return 0; /* handled */
    } else if (function_return_line(line) && !ignore_line(fun)) {
	err = parse_function_return(fun, line);
    } else {  
	err = strings_array_add(&fun->lines, &fun->n_lines, line);
    }

    if (err && !editing) {
	set_compiling_off();
    }	

    if (!err && !compiling) {
	/* finished composing function */
	err = check_function_structure(fun);
    }

    if (err && !editing) {
	delete_ufunc_from_list(fun);
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

static int fndef_maybe_append_next (char *s, FILE *fp)
{
    int n = strlen(s);

    if (s[n-1] == ',' || s[n-1] == '\\') {
	/* definition continues? */
	char line[MAXLINE];

	if (fgets(line, sizeof line, fp) == NULL) {  
	    return E_PARSE;
	}

	if (s[n-1] == '\\') {
	    s[n-1] = '\0';
	    n--;
	}	

	tailstrip(line);   
	n += strlen(line);

	if (n >= MAXLINE) {
	    return E_DATA;
	}

	strcat(s, line);
    }

    return 0;
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
		err = fndef_maybe_append_next(s, fp);
		if (!err) {
		    err = parse_fn_definition(fun->name, &fun->params,
					      &fun->n_params, &fun->rettype,
					      s + 8, NULL, NULL);
		}
	    }
	} else {
	    err = real_function_append_line(s, fun);
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
   function.
*/

static int localize_list (fncall *call, const char *oldname, 
			  fn_param *fp, DATAINFO *pdinfo)
{
    const int *list;
    int level = fn_executing + 1;
    int i, vi, err;

    list = get_list_by_name(oldname);

    if (list == NULL) {
	err = E_DATA;
    } else {
	/* note: it's OK if the two names are the same,
	   since their "levels" will be different */
	err = copy_named_list_as(oldname, fp->name);
    }

#if UDEBUG
    fprintf(stderr, "localize_list: localizing '%s' as '%s'\n",
	    oldname, fp->name);
#endif

    if (!err) {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0) {
		if (!in_gretl_list(call->listvars, vi)) {
		    gretl_list_append_term(&call->listvars, vi);
		}
		STACK_LEVEL(pdinfo, vi) = level;
	    }
	}
    }

    return err;
}

static void boolify_local_var (const char *vname)
{
    double x = gretl_scalar_get_value(vname);

    if (x != 0.0 && !na(x)) {
	gretl_scalar_set_value(vname, 1.0);
    }
}

static void maybe_set_arg_const (struct fnarg *arg, fn_param *fp)
{
    if ((fp->flags & ARG_CONST) || object_is_const(fp->name)) {
	arg->name = fp->name;
	arg->flags |= ARG_CONST;
    }
}

static int localize_const_matrix (struct fnarg *arg, fn_param *fp)
{
    user_matrix *u = get_user_matrix_by_data(arg->val.m);
    int err = 0;

    if (u == NULL) {
	/* the const argument is an anonymous matrix */
	err = matrix_add_as_shell(arg->val.m, fp->name);
    } else {
	/* the const argument is a named matrix */
	arg->upname = gretl_strdup(user_matrix_get_name(u));
	if (arg->upname == NULL) {
	    err = E_ALLOC;
	} else {
	    user_matrix_adjust_level(u, 1);
	    user_matrix_set_name(u, fp->name);
	}
    }

    if (!err) {
	arg->name = fp->name;
	arg->flags |= ARG_CONST;
    }

    return err;
}

static int localize_matrix_ref (struct fnarg *arg, fn_param *fp)
{
    user_matrix *u = arg->val.um;

    arg->upname = gretl_strdup(user_matrix_get_name(u));
    if (arg->upname == NULL) {
	return E_ALLOC;
    }

    user_matrix_adjust_level(u, 1);
    user_matrix_set_name(u, fp->name);

    maybe_set_arg_const(arg, fp);

    return 0;
}

static int localize_scalar_ref (fncall *call, struct fnarg *arg, 
				fn_param *fp, DATAINFO *pdinfo)
{
    int i = arg->val.idnum;
    const char *s = gretl_scalar_get_name(i);

    if (s == NULL) {
	return E_DATA;
    } 

    arg->upname = gretl_strdup(s);
    if (arg->upname == NULL) {
	return E_ALLOC;
    } 

    gretl_scalar_set_local_name(i, fp->name);
    maybe_set_arg_const(arg, fp);

    return 0;
}

static int localize_series_ref (fncall *call, struct fnarg *arg, 
				fn_param *fp, DATAINFO *pdinfo)
{
    int v = arg->val.idnum;

    arg->upname = gretl_strdup(pdinfo->varname[v]);
    if (arg->upname == NULL) {
	return E_ALLOC;
    } 

    STACK_LEVEL(pdinfo, v) += 1;
    strcpy(pdinfo->varname[v], fp->name);

    if (!in_gretl_list(call->ptrvars, v)) {
	gretl_list_append_term(&call->ptrvars, v);
    }

    maybe_set_arg_const(arg, fp);

    return 0;
}

/* Scalar function arguments only: if the arg is not supplied, use the
   default that is contained in the function specification, if any.
*/

static int add_scalar_arg_default (fn_param *param)
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
    
    return gretl_scalar_add_as_arg(param->name, x);
}

static void fncall_finalize_listvars (fncall *call)
{
    int i, v;

    for (i=call->listvars[0]; i>0; i--) {
	v = call->listvars[i];
	if (in_gretl_list(call->ptrvars, v)) {
	    gretl_list_delete_at_pos(call->listvars, i);
	}
    }

    if (call->listvars[0] == 0) {
	free(call->listvars);
	call->listvars = NULL;
    }
}

static int allocate_function_args (fncall *call,
				   double ***pZ,
				   DATAINFO *pdinfo)
{
    ufunc *fun = call->fun;
    fnargs *args = call->args;
    struct fnarg *arg;
    fn_param *fp;
    int i, err = 0;

    for (i=0; i<args->argc && !err; i++) {
	arg = args->arg[i];
	fp = &fun->params[i];

	if (gretl_scalar_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_USCALAR) {
		double x = gretl_scalar_get_value_by_index(arg->val.idnum);

		err = gretl_scalar_add_as_arg(fp->name, x);
	    } else if (arg->type == GRETL_TYPE_NONE) {
		err = add_scalar_arg_default(fp);
	    } else if (arg->type == GRETL_TYPE_MATRIX) {
		/* "cast" to scalar */
		err = gretl_scalar_add_as_arg(fp->name, arg->val.m->val[0]);
	    } else {
		err = gretl_scalar_add_as_arg(fp->name, arg->val.x);    
	    }
	    if (!err && fp->type == GRETL_TYPE_BOOL) {
		boolify_local_var(fp->name);
	    }
	} else if (fp->type == GRETL_TYPE_SERIES) {
	    if (arg->type == GRETL_TYPE_USERIES) {
		err = dataset_copy_variable_as(arg->val.idnum, fp->name,
					       pZ, pdinfo);
	    } else {
		err = dataset_add_series_as(arg->val.px, fp->name, 
					    pZ, pdinfo);
	    }	    
	} else if (fp->type == GRETL_TYPE_MATRIX) {
	    if (fp->flags & ARG_CONST) {
		err = localize_const_matrix(arg, fp);
	    } else {
		err = copy_matrix_as(arg->val.m, fp->name);
	    }
	} else if (fp->type == GRETL_TYPE_LIST) {
	    if (arg->type == GRETL_TYPE_NONE) {
		err = create_named_null_list(fp->name);
	    } else {
		err = localize_list(call, arg->val.str, fp, pdinfo);
	    }
	} else if (fp->type == GRETL_TYPE_STRING) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = add_string_as(arg->val.str, fp->name);
	    }
	} else if (fp->type == GRETL_TYPE_SCALAR_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_scalar_ref(call, arg, fp, pdinfo);
	    }
	} else if (fp->type == GRETL_TYPE_SERIES_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_series_ref(call, arg, fp, pdinfo);
	    }
	} else if (fp->type == GRETL_TYPE_MATRIX_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_matrix_ref(arg, fp);
	    }
	}

	if (!err) {
	    if (arg->type == GRETL_TYPE_USERIES) {
		arg->upname = gretl_strdup(pdinfo->varname[arg->val.idnum]);
	    } else if (arg->type == GRETL_TYPE_USCALAR) {
		arg->upname = gretl_strdup(gretl_scalar_get_name(arg->val.idnum));
	    }	
	}	
    }

    /* now for any parameters without matching arguments */

    for (i=args->argc; i<fun->n_params && !err; i++) {
	fp = &fun->params[i];
	if (gretl_scalar_type(fp->type)) {
	    err = add_scalar_arg_default(fp);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    err = create_named_null_list(fp->name);
	}
    }

    if (call->listvars != NULL) {
	if (err) {
	    free(call->listvars);
	    call->listvars = NULL;
	} else {
	    fncall_finalize_listvars(call);
	}
    }

    if (!err) {
	set_listargs_from_call(call, pdinfo);
    }

    return err;
}

int check_function_needs (const DATAINFO *pdinfo, FuncDataReq dreq,
			  int minver)
{
    static int thisver = 0;

    if (thisver == 0) {
	thisver = version_number_from_string(GRETL_VERSION);
    }

    if (minver > thisver) {
	char vstr[8];

	get_version_string(vstr, minver);
	sprintf(gretl_errmsg, "This function needs gretl version %s", vstr);
	return 1;
    }

    if (dreq == FN_NEEDS_TS && 
	(pdinfo == NULL || !dataset_is_time_series(pdinfo))) {
	strcpy(gretl_errmsg, "This function needs time-series data");
	return 1;
    }

    if (dreq == FN_NEEDS_PANEL && 
	(pdinfo == NULL || !dataset_is_panel(pdinfo))) {
	strcpy(gretl_errmsg, "This function needs panel data");
	return 1;
    }

    if (dreq == FN_NEEDS_QM && 
	(pdinfo == NULL || !dataset_is_time_series(pdinfo) || 
	 (pdinfo->pd != 4 || pdinfo->pd != 12))) {
	strcpy(gretl_errmsg, "This function needs quarterly or monthly data");
	return 1;
    } 

    return 0;
}

static int maybe_check_function_needs (const DATAINFO *pdinfo,
				       const ufunc *fun)
{
    const fnpkg *pkg = ufunc_get_parent_package(fun);

    if (pkg == NULL) {
	return 0;
    } else {
	return check_function_needs(pdinfo, pkg->dreq, pkg->minver);
    }
}

static double get_scalar_return (const char *vname, int *err)
{
    if (gretl_is_scalar(vname)) {
	return gretl_scalar_get_value(vname);
    } else {
	*err = E_UNKVAR; /* FIXME */
	return NADBL;
    }
}

static double *get_series_return (const char *vname, double **Z, 
				  DATAINFO *pdinfo, int copy, 
				  int *err)
{
    int v = series_index(pdinfo, vname);
    double *x = NULL;

    if (!copy && v == 0) {
	copy = 1;
    }

    if (v >= 0 && v < pdinfo->v) {
	if (copy) {
	    x = copyvec(Z[v], pdinfo->n);
	    if (x == NULL) {
		*err = E_ALLOC;
	    }
	} else {
	    x = Z[v];
	    Z[v] = NULL;
	}
    } else {
	*err = E_UNKVAR;
    }

    return x;
}

static gretl_matrix *get_matrix_return (const char *mname, int copy, int *err)
{
    gretl_matrix *ret = NULL;

    if (copy) {
	gretl_matrix *m = get_matrix_by_name(mname);

	if (m != NULL) {
	    ret = gretl_matrix_copy(m);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    }
	} 
    } else {
	ret = steal_matrix_by_name(mname);
    }

    if (ret == NULL && !*err) {
	*err = E_UNKVAR;
    }

    return ret;
}

enum {
    LIST_DIRECT_RETURN,
    LIST_WAS_ARG
};

/* Deal with a list that exists at the level of a user-defined
   function whose execution is now terminating.  Note that this list
   may be the direct return value of the function, or it may have been
   given to the function as an argument.  This is flagged in the
   @status argument.  
*/

static int unlocalize_list (const char *lname, int status,
			    double **Z, DATAINFO *pdinfo)
{
    const int *list = get_list_by_name(lname);
    int d = gretl_function_depth();
    int upd = d - 1;
    int i, vi;

#if UDEBUG
    fprintf(stderr, "unlocalize_list: '%s', function depth = %d\n", lname, d);
    printlist(list, lname);
    fprintf(stderr, " pdinfo = %p, pdinfo->v = %d\n", (void *) pdinfo, pdinfo->v);
#endif	    

    if (list == NULL) {
	return E_DATA;
    }

    /* Note, AC, 2009-04-04.  The following is tricky, but I think it
       _may_ now be right. If the list we're looking at was given as a
       function argument we simply shunt all its members to the prior
       STACK_LEVEL.  But if the list is the direct return value from
       the function, we need to overwrite any variables at caller
       level that have been redefined within the function.
    */

    if (status == LIST_DIRECT_RETURN) {
	int t, j, overwrite;
	const char *vname;
	
	for (i=1; i<=list[0]; i++) {
	    overwrite = 0;
	    vi = list[i];
	    vname = pdinfo->varname[vi];
	    unset_var_listarg(pdinfo, vi);
	    if (vi > 0 && vi < pdinfo->v && STACK_LEVEL(pdinfo, vi) == d) {
		for (j=1; j<pdinfo->v; j++) { 
		    if (STACK_LEVEL(pdinfo, j) == upd && 
			!strcmp(pdinfo->varname[j], vname)) {
			overwrite = 1;
			break;
		    }
		}
		if (overwrite) {
		    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
			Z[j][t] = Z[vi][t];
		    }
		} else {
		    STACK_LEVEL(pdinfo, vi) = upd;
		}
	    }
#if UDEBUG
	    fprintf(stderr, " list-member var %d, '%s': ", vi, vname); 
	    if (overwrite) {
		fprintf(stderr, "found match in caller, overwrote var %d\n", j);
	    } else {
		fprintf(stderr, "no match in caller\n");
	    }
#endif
	}
    } else {
	/* LIST_WAS_ARG */
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    unset_var_listarg(pdinfo, vi);
	    STACK_LEVEL(pdinfo, vi) = upd;
	}
    }

    return 0;
}

static char *
get_string_return (const char *sname, double **Z, DATAINFO *pdinfo, 
		   int *err)
{
    
    const char *s = get_string_by_name(sname);
    char *ret = NULL;

    if (s == NULL) {
	*err = E_DATA;
    } else {
	ret = gretl_strdup(s);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

static void 
maybe_set_return_description (ufunc *u, int rtype, DATAINFO *pdinfo, 
			      char **descrip)
{
    if (rtype == GRETL_TYPE_SERIES) {
	int v = series_index(pdinfo, u->retname);

	if (v < pdinfo->v) {
	    *descrip = gretl_strdup(VARLABEL(pdinfo, v));
	}
    }
}

#define types_match(pt,rt) ((pt==GRETL_TYPE_SERIES_REF && rt==GRETL_TYPE_SERIES) || \
                            (pt==GRETL_TYPE_MATRIX_REF && rt==GRETL_TYPE_MATRIX))

static int is_pointer_arg (ufunc *u, fnargs *args, int rtype)
{
    int i;

    for (i=0; i<args->argc; i++) {
	if (types_match(u->params[i].type, rtype)) {
	    if (!strcmp(u->params[i].name, u->retname)) {
		return 1;
	    }
	}
    }

    return 0;
}

#define needs_datainfo(t) (t == GRETL_TYPE_SERIES || \
                           t == GRETL_TYPE_LIST || \
                           t == GRETL_TYPE_SERIES_REF)

static int 
function_assign_returns (ufunc *u, fnargs *args, int rtype, 
			 double **Z, DATAINFO *pdinfo, 
			 void *ret, char **descrip, PRN *prn, 
			 int *perr)
{
    struct fnarg *arg;
    fn_param *fp;
    int copy, i, err = 0;

#if UDEBUG
    fprintf(stderr, "function_assign_returns: retname = '%s', rtype = %d\n", 
	    u->retname, rtype);
#endif

    if (*perr == 0) {
	/* first we work on the value directly returned by the
	   function (but only if there's no error) */
	if (needs_datainfo(rtype) && pdinfo == NULL) {
	    /* "can't happen" */
	    err = E_DATA;
	} else if (rtype == GRETL_TYPE_DOUBLE) {
	    *(double *) ret = get_scalar_return(u->retname, &err);
	} else if (rtype == GRETL_TYPE_SERIES) {
	    copy = is_pointer_arg(u, args, rtype);
	    *(double **) ret = get_series_return(u->retname, Z, pdinfo, copy, &err);
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    copy = is_pointer_arg(u, args, rtype);
	    *(gretl_matrix **) ret = get_matrix_return(u->retname, copy, &err);
	} else if (rtype == GRETL_TYPE_LIST) {
	    /* note: in this case the job is finished in
	       stop_fncall(); here we just adjust the info on the
	       listed variables so they don't get deleted
	    */
	    err = unlocalize_list(u->retname, LIST_DIRECT_RETURN, Z, pdinfo);
	} else if (rtype == GRETL_TYPE_STRING) {
	    *(char **) ret = get_string_return(u->retname, Z, pdinfo, &err);
	}

	if (err == E_UNKVAR) {
	    pprintf(prn, "Function %s did not provide the specified return value\n",
		    u->name);
	}

	*perr = err;

	if (!err && pdinfo != NULL && descrip != NULL) {
	    maybe_set_return_description(u, rtype, pdinfo, descrip);
	}
    }

    /* "indirect return" values and other pointerized args: 
       these should be handled even if the function bombed.
    */

    for (i=0; i<args->argc; i++) {
	arg = args->arg[i];
	fp = &u->params[i];
	if (needs_datainfo(fp->type) && pdinfo == NULL) {
	    err = E_DATA;
	} else if (gretl_ref_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_SERIES_REF) {
		int v = arg->val.idnum;

		STACK_LEVEL(pdinfo, v) -= 1;
		strcpy(pdinfo->varname[v], arg->upname);
	    } else if (arg->type == GRETL_TYPE_SCALAR_REF) {
		gretl_scalar_restore_name(arg->val.idnum, arg->upname);
	    } else if (arg->type == GRETL_TYPE_MATRIX_REF) {
		user_matrix *u = arg->val.um;

		user_matrix_adjust_level(u, -1);
		user_matrix_set_name(u, arg->upname);
	    }
	} else if (fp->type == GRETL_TYPE_MATRIX && arg->upname != NULL) {
	    user_matrix *u = get_user_matrix_by_data(arg->val.m);

	    user_matrix_adjust_level(u, -1);
	    user_matrix_set_name(u, arg->upname);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    unlocalize_list(fp->name, LIST_WAS_ARG, Z, pdinfo);
	}
    }

#if UDEBUG
    fprintf(stderr, "function_assign_returns: returning %d\n", err);
#endif

    return err;
}

static void record_obs_info (obsinfo *o, DATAINFO *pdinfo)
{
    o->changed = 0;

    if (pdinfo != NULL) {
	o->structure = pdinfo->structure;
	o->pd = pdinfo->pd;
	o->t1 = pdinfo->t1;
	o->t2 = pdinfo->t2;
	strcpy(o->stobs, pdinfo->stobs);
    }
}

static int restore_obs_info (obsinfo *o, double ***pZ, DATAINFO *pdinfo)
{
    char line[128];
    gretlopt opt = OPT_NONE;

    if (o->structure == CROSS_SECTION) {
	opt = OPT_X;
    } else if (o->structure == TIME_SERIES) {
	opt = OPT_T;
    } else if (o->structure == STACKED_TIME_SERIES) {
	opt = OPT_S;
    } else if (o->structure == SPECIAL_TIME_SERIES) {
	opt = OPT_N;
    } 

    sprintf(line, "setobs %d %s", o->pd, o->stobs);

    return set_obs(line, pZ, pdinfo, opt);
}

static int stop_fncall (fncall *call, int rtype, void *ret,
			double ***pZ, DATAINFO *pdinfo,
			int orig_v)
{
    int i, d = gretl_function_depth();
    int delv, anyerr = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "stop_fncall: terminating call to "
	    "function '%s' at depth %d, pdinfo->v = %d\n", 
	    call->fun->name, d, (pdinfo != NULL)? pdinfo->v : 0);
#endif

    call->args = NULL;

    anyerr = destroy_user_scalars_at_level(d);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_user_scalars_at_level(%d): err = %d\n", d, err);
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

    /* below: delete variables local to the function, taking care not to
       delete any 'local' vars that have been "promoted" to caller
       level via their inclusion in a returned list
    */

    if (pdinfo != NULL) {
	for (i=orig_v, delv=0; i<pdinfo->v; i++) {
	    if (STACK_LEVEL(pdinfo, i) == d) {
		delv++;
	    }
	}
	if (delv > 0) {
	    if (delv == pdinfo->v - orig_v) {
		/* deleting all added variables */
		anyerr = dataset_drop_last_variables(delv, pZ, pdinfo);
		if (anyerr && !err) {
		    err = anyerr;
		}
	    } else {
		for (i=pdinfo->v-1; i>=orig_v; i--) {
		    if (STACK_LEVEL(pdinfo, i) == d) {
			anyerr = dataset_drop_variable(i, pZ, pdinfo);
			if (anyerr && !err) {
			    err = anyerr;
			}
		    } 
		}
	    }    
	}
    }

    /* direct list return: write the possibly revised list to the
       return pointer.  Note that we can't do this earlier, because
       the ID numbers of the variables in the return list may be
       changed due to the deletion of function-local variables.
    */
    if (!err && rtype == GRETL_TYPE_LIST) {
	int *lret = gretl_list_copy(get_list_by_name(call->fun->retname));

	if (lret != NULL) {
	    *(int **) ret = lret;
	} else {
	    err = E_ALLOC;
	}
    }    

    anyerr = destroy_saved_lists_at_level(d);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_saved_lists_at_level(%d): err = %d\n", d, err);
#endif
    }    

    pop_program_state();

    if (pdinfo != NULL && call->obs.changed) {
	restore_obs_info(&call->obs, pZ, pdinfo);
    }

    set_executing_off(call, pdinfo);

    return err;
}

static int start_fncall (fncall *call, DATAINFO *pdinfo, PRN *prn)
{
    fn_executing++;
    push_program_state();

    callstack = g_list_append(callstack, call);
#if EXEC_DEBUG
    fprintf(stderr, "start_fncall: added call to %s, depth now %d\n", 
	    call->fun->name, g_list_length(callstack));
#endif

    record_obs_info(&call->obs, pdinfo);

    if (gretl_debugging_on() || call->fun->debug) {
	set_gretl_echo(1);
	set_gretl_messages(1);
	pprintf(prn, "*** executing function %s\n", call->fun->name);
    } else {
	set_gretl_echo(0);
	set_gretl_messages(0);
    }

    return 0;
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

static double arg_get_double_val (struct fnarg *arg)
{
    if (arg->type == GRETL_TYPE_USCALAR) {
	return gretl_scalar_get_value_by_index(arg->val.idnum);
    } else if (gretl_scalar_type(arg->type)) {
	return arg->val.x;
    } else if (arg->type == GRETL_TYPE_MATRIX) {
	return arg->val.m->val[0];
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
	} else if (gretl_scalar_type(fp->type) && arg->type == GRETL_TYPE_USCALAR) {
	    ; /* OK */
	} else if (fp->type == GRETL_TYPE_SERIES && arg->type == GRETL_TYPE_USERIES) {
	    ; /* OK */
	} else if (gretl_scalar_type(fp->type) && 
		   arg->type == GRETL_TYPE_MATRIX &&
		   gretl_matrix_is_scalar(arg->val.m)) {
	    ; /* OK */
	} else if (fp->type != arg->type) {
	    pprintf(prn, _("%s: argument %d is of the wrong type (is %s, should be %s)\n"), 
		    u->name, i + 1, arg_type_string(arg->type), arg_type_string(fp->type));
	    err = E_TYPES;
	}

	if (!err && fp->type == GRETL_TYPE_DOUBLE) {
	    x = arg_get_double_val(arg);
	    if ((!na(fp->min) && x < fp->min) ||
		(!na(fp->max) && x > fp->max)) {
		pprintf(prn, _("%s, argument %d: value %g is out of bounds\n"), 
			u->name, i + 1, x);
		err = E_INVARG;
	    }
	}
    }

    for (i=args->argc; i<u->n_params && !err; i++) {
	/* do we have defaults for any empty args? */
	fp = &u->params[i];
	if (!(fp->flags & ARG_OPTIONAL) && na(fp->deflt)) {
	    pprintf(prn, _("%s: not enough arguments\n"), u->name);
	    err = E_ARGS;
	}
    }

    return err;
}

/**
 * user_function_set_debug:
 * @name: the name of the function.
 * @debug: boolean, if 1 then start debugging function, if
 * 0 then stop debugging.
 *
 * Enables or disables debugging for a user-defined function.
 *
 * Returns: 0 on success, non-zero if no function is specified
 * or if the function is not found.
 */

int user_function_set_debug (const char *name, int debug)
{
    ufunc *fun;

    if (name == NULL || *name == '\0') {
	return E_ARGS;
    } 

    fun = get_user_function_by_name(name);

    if (fun == NULL) {
	return E_UNKVAR;
    } else {
	fun->debug = debug;
	return 0;
    }
}

#define debug_cont(c) (c->ci == FUNDEBUG && (c->opt & OPT_C))
#define debug_next(c) (c->ci == FUNDEBUG && (c->opt & OPT_N))

#define set_debug_cont(c) (c->ci = FUNDEBUG, c->opt = OPT_C)
#define set_debug_next(c) (c->ci = FUNDEBUG, c->opt = OPT_N)

/* loop for stepping through function commands; returns
   1 if we're debugging, otherwise 0 
*/

static int debug_command_loop (ExecState *state, 
			       double ***pZ,
			       DATAINFO *datainfo,
			       DEBUG_READLINE get_line,
			       DEBUG_OUTPUT put_func,
			       int *errp)
{
    int brk = 0, err = 0;

    state->flags |= DEBUG_EXEC;

    while (!brk) {
#if DDEBUG
	fprintf(stderr, "--- debug_command_loop calling get_line\n"); 
#endif
	*errp = (*get_line)(state);
	if (*errp) {
	    /* the debugger failed: get out right away */
	    state->flags &= ~DEBUG_EXEC;
	    return 0;
	}

	err = parse_command_line(state->line, state->cmd, 
				 pZ, datainfo);
	if (err) {
	    if (!strcmp(state->line, "c")) {
		/* short for 'continue' */
		set_debug_cont(state->cmd);
		err = 0;
	    } else if (!strcmp(state->line, "n")) {
		/* short for 'next' */
		set_debug_next(state->cmd);
		err = 0;
	    }
	}

	if (err || debug_cont(state->cmd) || debug_next(state->cmd)) {
	    brk = 1;
	} else {
	    /* execute interpolated command */
	    err = gretl_cmd_exec(state, pZ, datainfo);
	    if (put_func != NULL) {
#if DDEBUG		
		fprintf(stderr, "--- debug_command_loop calling put_func\n"); 
#endif
		(*put_func)(state);
	    }		
	}
    }

    if (!debug_next(state->cmd)) {
	state->flags &= ~DEBUG_EXEC;
    }

    return 1 + debug_next(state->cmd);
}

#define void_function(f) (f->rettype == 0 || f->rettype == GRETL_TYPE_VOID)

static int handle_return_statement (ufunc *fun,
				    ExecState *state,
				    double ***pZ,
				    DATAINFO *pdinfo)
{
    const char *s = state->line + 6; /* skip "return" */
    int err = 0;

#if EXEC_DEBUG
    fprintf(stderr, "%s: return: s = '%s'\n", fun->name, s);
#endif

    s += strspn(s, " ");

    if (*s == '\0' && void_function(fun)) {
	/* plain "return" from void function: OK */
	return 0;
    } else if (*s == '\0' && !void_function(fun)) {
	gretl_errmsg_sprintf("%s: return value is missing", fun->name);
	err = E_TYPES;
    } else if (*s != '\0' && void_function(fun)) {
	gretl_errmsg_sprintf("%s: non-null return value '%s' is not valid", 
			     fun->name, s);
	err = E_TYPES;
    } else {
	int len = gretl_namechar_spn(s);

	if (len == strlen(s)) {
	    /* returning a named variable */
	    fun->retname = gretl_strndup(s, len);
	} else {
	    const char *typestr = arg_type_string(fun->rettype);
	    char formula[MAXLINE];
	    
	    sprintf(formula, "%s $retval=%s", typestr, s);
	    err = generate(formula, pZ, pdinfo, OPT_P, NULL);
	    if (!err) {
		fun->retname = gretl_strdup("$retval");
	    }
	}
    }

    if (!err && fun->retname == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static void fn_state_init (CMD *cmd, ExecState *state, int *indent0)
{
    cmd->list = NULL;
    cmd->param = NULL;
    cmd->extra = NULL;
    cmd->linfo = NULL;

    state->cmd = NULL;
    state->models = NULL;
    state->submask = NULL;

    *indent0 = gretl_if_state_record();
}

int gretl_function_exec (ufunc *u, fnargs *args, int rtype,
			 double ***pZ, DATAINFO *pdinfo,
			 void *ret, char **descrip, 
			 PRN *prn)
{
    DEBUG_READLINE get_line = NULL;
    DEBUG_OUTPUT put_func = NULL;
    ExecState state;
    fncall *call = NULL;
    MODEL **models = NULL;
    char line[MAXLINE];
    CMD cmd;
    int orig_v = 0;
    int orig_t1 = 0;
    int orig_t2 = 0;
    int indent0, started = 0;
    int retline = -1;
    int debugging = u->debug;
    int i, err = 0;

    *funcerr_msg = '\0';

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: starting %s\n", u->name);
#endif

    err = maybe_check_function_needs(pdinfo, u);
    if (err) {
	return err;
    }

    if (pdinfo != NULL) {
	orig_v = pdinfo->v;
	orig_t1 = pdinfo->t1;
	orig_t2 = pdinfo->t2;
    }

    /* precaution */
    fn_state_init(&cmd, &state, &indent0);

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: argc = %d\n", args->argc);
    fprintf(stderr, "u->n_params = %d\n", u->n_params);
#endif

    call = fncall_new(u);
    if (call == NULL) {
	fprintf(stderr, "fncall_new() returned NULL\n");
	return E_ALLOC;
    }

    err = check_function_args(u, args, (pZ != NULL)? (const double **) *pZ : NULL,
			      pdinfo, prn);

    if (!err) {
	call->args = args;
	err = allocate_function_args(call, pZ, pdinfo);
    }

    if (err) {
	/* get out before allocating further storage */
	fncall_free(call);
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
	if (pdinfo != NULL && pdinfo->submask != NULL) {
	    state.submask = copy_datainfo_submask(pdinfo);
	}
	state.callback = func_exec_callback;
    }

    if (!err) {
	err = start_fncall(call, pdinfo, prn);
	if (!err) {
	    started = 1;
	}
    }

#if EXEC_DEBUG
    fprintf(stderr, "start_fncall: err = %d\n", err);
#endif

    if (debugging) {
	get_line = get_debug_read_func();
	put_func = get_debug_output_func();
	if (get_line != NULL) {
	    debugging = 2;
	}
    }

    /* get function lines in sequence and check, parse, execute */

    for (i=0; i<u->n_lines && !err; i++) {

	if (*u->lines[i] == '\0') {
	    continue;
	}

	strcpy(line, u->lines[i]);

	if (debugging) {
	    pprintf(prn, "%s> %s\n", u->name, line);
	} else if (gretl_echo_on()) {
	    pprintf(prn, "? %s\n", line);
	}

	err = maybe_exec_line(&state, pZ, pdinfo);

	if (!err && state.cmd->ci == FUNCRET) {
	    err = handle_return_statement(u, &state, pZ, pdinfo); 
	    retline = i;
	    break;
	}

	if (debugging > 1 && state.cmd->ci > 0 && 
	    !gretl_compiling_loop() && !state.cmd->context) {
	    if (put_func != NULL) {
		pprintf(prn, "-- debugging %s, line %d --\n", u->name, i + 1);
		(*put_func)(&state);
	    } else {
		pprintf(prn, "-- debugging %s, line %d --\n", u->name, i + 1);
	    }
	    debugging = debug_command_loop(&state, pZ, pdinfo,
					   get_line, put_func, 
					   &err);
	}

	if (err) {
	    if (*gretl_errmsg == '\0') {
		gretl_errmsg_sprintf("error in function %s\n"
				     "> %s", u->name, line);
	    }
	    fprintf(stderr, "error on line %d of function %s\n", i+1, u->name);
	    fprintf(stderr, "> %s\n", line);
	}

	if (state.funcerr) {
	    pprintf(prn, "%s: %s\n", u->name, state.cmd->param);
	    set_funcerr_message(u, state.cmd->param);
	}

	if (state.cmd->ci == SETOBS) {
	    /* set flag for reverting on exit */
	    call->obs.changed = 1;
	}

	if (gretl_execute_loop()) { 
#if EXEC_DEBUG
	    fprintf(stderr, "gretl_function_exec: calling gretl_loop_exec\n");
#endif
	    err = gretl_loop_exec(&state, pZ, pdinfo);
	    if (state.funcerr) {
		pprintf(prn, "%s: %s\n", u->name, state.cmd->param);
		set_funcerr_message(u, state.cmd->param);
	    }
	    if (err) {
		fprintf(stderr, "function_exec: breaking on error %d in loop\n", err);
		fprintf(stderr, "error on line %d of function %s\n", i+1, u->name);
		break;
	    }
#if EXEC_DEBUG
	    fprintf(stderr, "gretl_function_exec: gretl_loop_exec done, err = %d\n", err);
#endif
	}
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: %s: finished main exec, err = %d, pdinfo->v = %d\n", 
	    u->name, err, (pdinfo != NULL)? pdinfo->v : 0);
#endif

    if (pdinfo != NULL) {
	/* restore the sample that was in place on entry */
	if (complex_subsampled()) {
	    if (state.submask == NULL) {
		/* we were not sub-sampled on entry */
		restore_full_sample(pZ, pdinfo, NULL);
	    } else if (submask_cmp(state.submask, pdinfo->submask)) {
		/* we were sub-sampled differently on entry */
		restore_full_sample(pZ, pdinfo, NULL);
		restrict_sample_from_mask(state.submask, pZ, pdinfo, OPT_NONE);
	    } 
	}
	pdinfo->t1 = orig_t1;
	pdinfo->t2 = orig_t2;
    }

    if (err || retline >= 0) {
	gretl_if_state_clear();
    } else {
	err = gretl_if_state_check(indent0);
    }

    function_assign_returns(u, args, rtype, (pZ != NULL)? *pZ : NULL, 
			    pdinfo, ret, descrip, prn, &err);

    gretl_exec_state_clear(&state);

    if (started) {
	int stoperr = stop_fncall(call, rtype, ret, pZ, pdinfo, orig_v);

	if (stoperr && !err) {
	    err = stoperr;
	}
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: err = %d\n", err);
#endif

    return err;    
}

/* look up name of supplied argument based on name of variable
   inside function */

char *gretl_func_get_arg_name (const char *argvar, int *err)
{
    fncall *call = current_function_call();
    char *ret = NULL;

    *err = E_DATA;

    if (call != NULL && call->args != NULL) {
	ufunc *u = call->fun;
	fnargs *args = call->args;
	int i, n = args->argc;

	for (i=0; i<n; i++) {
	    if (!strcmp(argvar, u->params[i].name)) {
		*err = 0;
		if (args->arg[i]->upname != NULL) {
		    ret = gretl_strdup(args->arg[i]->upname); 
		    if (ret == NULL) {
			*err = E_ALLOC;
		    }
		}
		break;
	    }
	}
    }

    return ret;
}

int object_is_const (const char *name)
{
    fncall *call = current_function_call();

    if (call != NULL && call->args != NULL) {
	fnargs *args = call->args;
	int i, n = args->argc;

	for (i=0; i<n; i++) {
	    const char *aname = args->arg[i]->name;

	    if (aname != NULL && !strcmp(name, aname)) {
		return args->arg[i]->flags & ARG_CONST;
	    }
	}
    }

    return 0;
}

void sample_range_get_extrema (const DATAINFO *pdinfo, int *t1, int *t2)
{
    fncall *call = current_function_call();

    if (call != NULL) {
	*t1 = call->obs.t1;
	*t2 = call->obs.t2;
    } else {
	*t1 = 0;
	*t2 = pdinfo->n - 1;
    }
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

    if (fun->rettype != GRETL_TYPE_NONE && fun->rettype != GRETL_TYPE_VOID) {
	pprintf(prn, "Return value: %s\n\n", arg_type_string(fun->rettype));
    } else {
	pputs(prn, "Return value: none\n\n");
    }
	
    if (fun->help != NULL) {
	pputs(prn, "Help text:\n");
	pputs(prn, fun->help);
	pprintf(prn, "\n\n");
    }

    if (pkg != NULL && pkg->sample != NULL) {
	pputs(prn, "Sample script:\n");
	pputs(prn, pkg->sample);
	pprintf(prn, "\n\n");
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
