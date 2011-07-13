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

#define FULL_XML_HEADERS

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
#include "gretl_bundle.h"
#include "flow_control.h"
#include "kalman.h"

#include <glib.h>

#define FNPARSE_DEBUG 0 /* debug parsing of function code */
#define EXEC_DEBUG 0    /* debugging of function execution */
#define UDEBUG 0        /* debug handling of args */
#define PKG_DEBUG 0     /* debug handling of function packages */
#define FN_DEBUG 0      /* miscellaneous debugging */
#define DDEBUG 0        /* debug the debugger */

#define INT_USE_XLIST (-999)

typedef struct fn_param_ fn_param;
typedef struct obsinfo_ obsinfo;
typedef struct fncall_ fncall;

/* structure representing a parameter of a user-defined function */

struct fn_param_ {
    char *name;     /* the name of the parameter */
    char type;      /* its type */
    char *descrip;  /* its description */
    char **labels;  /* value labels, if applicable */
    int nlabels;    /* number of value labels */
    char flags;     /* additional information (e.g. "const" flag) */
    double deflt;   /* default value */
    double min;     /* minimum value (scalar parameters only) */
    double max;     /* maximum value (scalar parameters only) */
    double step;    /* step increment (scalars only) */
};

/* structure representing sample information at start of
   a function call */

struct obsinfo_ {
    int structure;      /* time-series, etc. */
    int pd;             /* data frequency */
    int t1, t2;         /* starting and ending observations */
    char changed;       /* sample has been changed within the function call? */
    char stobs[OBSLEN]; /* string representation of starting obs */
};

/* structure representing a call to a user-defined function */

struct fncall_ {
    ufunc *fun;    /* the function called */
    fnargs *args;  /* argument array */
    int *ptrvars;  /* list of pointer arguments */
    int *listvars; /* list of series included in a list argument */
    char *retname; /* name of return value (or dummy string) */
    obsinfo obs;   /* sample info */
};

enum {
    UFUN_PRIVATE  = 1 << 0,
    UFUN_PLUGIN   = 1 << 1,
    UFUN_NOPRINT  = 1 << 2  
};

enum {
    UFUN_BUNDLE_PRINT = 1,
    UFUN_BUNDLE_PLOT,
    UFUN_BUNDLE_TEST,
    UFUN_BUNDLE_FCAST,
    UFUN_BUNDLE_EXTRA,
    UFUN_GUI_MAIN,
    UFUN_GUI_PRECHECK
};    

/* structure representing a user-defined function */

struct ufunc_ {
    char name[FN_NAMELEN]; /* identifier */
    fnpkg *pkg;            /* pointer to parent package, or NULL */
    int pkg_role;          /* printer, plotter, etc. */
    int flags;             /* private, plugin, etc. */
    int n_lines;           /* number of lines of code */
    char **lines;          /* array of lines of code */
    int n_params;          /* number of parameters */
    fn_param *params;      /* parameter info array */
    int rettype;           /* return type (if any) */
    int debug;             /* are we debugging this function? */
};

/* structure representing a function package */

struct fnpkg_ {
    int ID;           /* unique ID (time stamp) */
    char name[FN_NAMELEN]; /* package name */
    char *fname;      /* filename */
    char *author;     /* author's name */
    char *version;    /* package version string */
    char *date;       /* last revision date */
    char *descrip;    /* package description */
    char *help;       /* package help text */
    char *gui_help;   /* GUI-specific help (optional) */
    char *sample;     /* sample caller script */
    char *label;      /* for use in GUI menus */
    int minver;       /* minimum required gretl version */
    FuncDataReq dreq; /* data requirement */
    ufunc **pub;      /* pointers to public interfaces */
    ufunc **priv;     /* pointers to private functions */
    int n_pub;        /* number of public functions */
    int n_priv;       /* number of private functions */
};

/* acceptable types for parameters of user-defined functions */

#define ok_function_arg_type(t) (t == GRETL_TYPE_BOOL ||	\
				 t == GRETL_TYPE_INT ||		\
				 t == GRETL_TYPE_OBS ||		\
				 t == GRETL_TYPE_DOUBLE ||	\
				 t == GRETL_TYPE_SERIES ||	\
				 t == GRETL_TYPE_LIST ||	\
				 t == GRETL_TYPE_MATRIX ||	\
				 t == GRETL_TYPE_STRING ||	\
				 t == GRETL_TYPE_SCALAR_REF ||	\
				 t == GRETL_TYPE_SERIES_REF ||	\
				 t == GRETL_TYPE_MATRIX_REF ||	\
				 t == GRETL_TYPE_BUNDLE_REF ||  \
	                         t == GRETL_TYPE_BUNDLE)

enum {
    ARG_OPTIONAL = 1 << 0,
    ARG_CONST    = 1 << 1
};

/* structure representing an argument to a user-defined function */

struct fnarg {
    char type;           /* argument type */
    char flags;          /* ARG_OPTIONAL, ARG_CONST as appropriate */
    const char *name;    /* name as function param */
    char *upname;        /* name of supplied arg at caller level */
    union {
	int idnum;        /* named series arg (series ID) */
	double x;         /* scalar arg */
	double *px;       /* anonymous series arg */
	gretl_matrix *m;  /* matrix arg */
	user_matrix *um;  /* named user-matrix arg */
	char *str;        /* string arg */
	gretl_bundle *b;  /* anonymous bundle pointer */
    } val;
};

struct fnargs_ {
    int argc;           /* count of arguments */
    struct fnarg **arg; /* array of arguments */
};

static int n_ufuns;         /* number of user-defined functions in memory */
static ufunc **ufuns;       /* array of pointers to user-defined functions */
static ufunc *current_fdef; /* pointer to function currently being defined */
static GList *callstack;    /* stack of function calls */
static int n_pkgs;          /* number of loaded function packages */
static fnpkg **pkgs;        /* array of pointers to loaded packages */

static int function_package_record (fnpkg *pkg);
static void function_package_free (fnpkg *pkg);
static int version_number_from_string (const char *s);

/* record of state, and communication of state with outside world */

static int compiling;    /* boolean: are we compiling a function currently? */
static int fn_executing; /* depth of function call stack */

#define function_is_private(f) (f->flags & UFUN_PRIVATE)
#define function_is_plugin(f)  (f->flags & UFUN_PLUGIN)
#define function_is_noprint(f) (f->flags & UFUN_NOPRINT)

#define null_return(t) (t == 0 || t == GRETL_TYPE_VOID)

struct flag_and_key {
    int flag;
    const char *key;
};

static struct flag_and_key pkg_lookups[] = {
    { UFUN_BUNDLE_PRINT, BUNDLE_PRINT },
    { UFUN_BUNDLE_PLOT,  BUNDLE_PLOT },
    { UFUN_BUNDLE_TEST,  BUNDLE_TEST },
    { UFUN_BUNDLE_FCAST, BUNDLE_FCAST },
    { UFUN_BUNDLE_EXTRA, BUNDLE_EXTRA },
    { UFUN_GUI_MAIN,     GUI_MAIN },
    { UFUN_GUI_PRECHECK, GUI_PRECHECK },
    { -1,                NULL }
};

static int pkg_key_get_role (const char *key)
{
    int i;

    for (i=0; pkg_lookups[i].flag > 0; i++) {
	if (!strcmp(key, pkg_lookups[i].key)) {
	    return pkg_lookups[i].flag;
	}
    }
    
    return 0;
}

static const char *pkg_role_get_key (int flag)
{
    int i;

    for (i=0; pkg_lookups[i].flag > 0; i++) {
	if (flag == pkg_lookups[i].flag) {
	    return pkg_lookups[i].key;
	}
    }
    
    return NULL;
}

static void set_function_private (ufunc *u, gboolean s)
{
    if (s) {
	u->flags |= UFUN_PRIVATE;
    } else {
	u->flags &= ~UFUN_PRIVATE;
    }
}

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

/**
 * fn_arg_new:
 * @type: type of argument
 * @p: pointer to value for argument.
 * @err: location to receive error code.
 *
 * Allocates a new function argument of the specified @type
 * and assigns it the value given in @p.
 *
 * Returns: the allocated argument, or %NULL of failure.
 */

struct fnarg *fn_arg_new (GretlType type, void *p, int *err)
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
    } else if (type == GRETL_TYPE_INT || type == GRETL_TYPE_OBS) {
	arg->val.x = *(int *) p;
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
    } else if (type == GRETL_TYPE_BUNDLE_REF) {
	arg->val.str = (char *) p;
    } else if (type == GRETL_TYPE_BUNDLE) {
	arg->val.b = (gretl_bundle *) p;
    } else {
	*err = E_TYPES;
	free(arg);
	arg = NULL;
    }

    return arg;
}

/**
 * fn_args_new:
 *
 * Returns: a newly allocated, empty array of function
 * arguments, or %NULL of failure.
 */

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
    if (args != NULL) {
	int i;

	for (i=0; i<args->argc; i++) {
	    free_fn_arg(args->arg[i]);
	}
    
	free(args->arg);
	free(args);
    }
}

/**
 * push_fn_arg:
 * @args: existing array of function arguments.
 * @type: type of argument to add.
 * @p: pointer to value to add.
 *
 * Appends a new argument of the specified type and value to the
 * array @args.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_fn_arg (fnargs *args, GretlType type, void *p)
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
	call->retname = NULL;
    }

    return call;
}

static void fncall_free (fncall *call)
{
    if (call != NULL) {
	free(call->ptrvars);
	free(call->listvars);
	free(call->retname);
	free(call);
    }
}

static fnpkg *function_package_alloc (const char *fname)
{
    fnpkg *pkg = malloc(sizeof *pkg);

    if (pkg == NULL) {
	return NULL;
    }

    pkg->fname = gretl_strdup(fname);
    if (pkg->fname == NULL) {
	free(pkg);
	return NULL;
    }     

    pkg->ID = (int) time(NULL);
    pkg->name[0] = '\0';
    pkg->author = NULL;
    pkg->version = NULL;
    pkg->date = NULL;
    pkg->descrip = NULL;
    pkg->help = NULL;
    pkg->gui_help = NULL;
    pkg->sample = NULL;
    pkg->label = NULL;
    pkg->dreq = 0;
    pkg->minver = 0;
    
    pkg->pub = pkg->priv = NULL;
    pkg->n_pub = pkg->n_priv = 0;

    return pkg;
}

/* For function call @call, set the 'listarg' status of any named
   series provided to the called function via list arguments.  This is
   required for correct handling of namespaces.  
*/

static void set_listargs_from_call (fncall *call, DATASET *dset)
{
    int i, v;

    if (dset == NULL) {
	return;
    }

    for (i=1; i<dset->v; i++) {
	unset_var_listarg(dset, i);
    }

    if (call != NULL && call->listvars != NULL) {
	for (i=1; i<=call->listvars[0]; i++) {
	    v = call->listvars[i];
#if UDEBUG
	    fprintf(stderr, "setting listarg status on var %d (%s)\n",
		    v, dset->varname[v]);
#endif
	    set_var_listarg(dset, v);
	}
    }
}

static void set_executing_off (fncall *call, DATASET *dset, PRN *prn)
{
    int dbg = gretl_debugging_on();
    fncall *popcall = NULL;

    fn_executing--;

    callstack = g_list_remove(callstack, call);

#if EXEC_DEBUG
    fprintf(stderr, "set_executing_off: removing call to %s, depth now %d\n",
	    call->fun->name, g_list_length(callstack));
#endif

    if (dbg) {
	pprintf(prn, "*** exiting function %s, ", call->fun->name);
    }

    fncall_free(call);

    if (fn_executing > 0) {
	GList *tmp = g_list_last(callstack);
	
	popcall = tmp->data;
    } else {
	g_list_free(callstack);
	callstack = NULL;
	gretl_insert_builtin_string("pkgdir", NULL);
    }

    if (dset != NULL) {
	set_listargs_from_call(popcall, dset);
    }

    if (dbg) {
	if (popcall != NULL) {
	    pprintf(prn, "returning to %s\n", popcall->fun->name);
	} else {
	    pputs(prn, "returning to main\n");
	}
    }
}

/**
 * n_free_functions:
 *
 * Returns: the number of functions loaded in memory
 * that are not currently attached to any function package,
 * and are therefore available for packaging.
 */

int n_free_functions (void)
{
    int i, n = 0;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == NULL) {
	    n++;
	}
    }

    return n;
}

/**
 * get_user_function_by_index:
 * @idx: index number.
 *
 * Returns: pointer to the user-function that currently
 * occupies (0-based) slot @idx in the array of loaded
 * functions, or %NULL if @idx is out of bounds.
 */

const ufunc *get_user_function_by_index (int idx)
{
    return (idx < 0 || idx >= n_ufuns)? NULL : ufuns[idx];
}

/**
 * fn_n_params:
 * @fun: pointer to user-function.
 *
 * Returns: the number of parameters associated with @fun.
 */

int fn_n_params (const ufunc *fun)
{
    return fun->n_params;
}

/**
 * fn_param_type:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the type of parameter @i of function
 * @fun.
 */

int fn_param_type (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? 0 :
	fun->params[i].type;
}

/**
 * fn_param_name:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the name of parameter @i of function
 * @fun.
 */

const char *fn_param_name (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NULL :
	fun->params[i].name;
}

/**
 * fn_param_descrip:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the description of parameter @i of function
 * @fun (if any), otherwise %NULL.
 */

const char *fn_param_descrip (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NULL :
	fun->params[i].descrip;
}

/**
 * fn_param_value_labels:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * @n: location to receive number of labels.
 * 
 * Returns: the value-labels associated with parameter @i 
 * of function @fun (if any), otherwise %NULL.
 */

const char **fn_param_value_labels (const ufunc *fun, int i, 
				    int *n)
{
    if (i >= 0 && i < fun->n_params) {
	*n = fun->params[i].nlabels;
	return (const char **) fun->params[i].labels;
    } else {
	*n = 0;
	return NULL;
    }
}

/**
 * fn_param_default:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the default value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_default (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].deflt;
}    

/**
 * fn_param_minval:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the minimum value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_minval (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].min;
}

/**
 * fn_param_maxval:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the maximum value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_maxval (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].max;
}

/**
 * fn_param_step:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: the step value for (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

double fn_param_step (const ufunc *fun, int i)
{
    return (i < 0 || i >= fun->n_params)? NADBL :
	fun->params[i].step;
}

/**
 * fn_param_optional:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: 1 if parameter @i of function @fun is optional,
 * otherwise 0.
 */

int fn_param_optional (const ufunc *fun, int i)
{
    int t;

    if (i < 0 || i >= fun->n_params) {
	return 0;
    }

    t = fun->params[i].type;

    return ((gretl_ref_type(t) || 
	     t == GRETL_TYPE_LIST ||
	     t == GRETL_TYPE_STRING) && 
	    (fun->params[i].flags & ARG_OPTIONAL));
}

/**
 * fn_param_uses_xlist:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * 
 * Returns: 1 if parameter @i of function @fun is
 * designed to select an integer based on a gretl
 * model's list of regressors, otherwise 0.
 */

int fn_param_uses_xlist (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
	return 0;
    }

    return (fun->params[i].type == GRETL_TYPE_INT &&
	    fun->params[i].deflt == INT_USE_XLIST);
}

/**
 * user_func_get_return_type:
 * @fun: pointer to user-function.
 * 
 * Returns: the return type of function @fun.
 */

int user_func_get_return_type (const ufunc *fun)
{
    if (fun == NULL) {
	return GRETL_TYPE_NONE;
    } else {
	return fun->rettype;
    }
}

/**
 * user_func_is_noprint:
 * @fun: pointer to user-function.
 * 
 * Returns: 1 if the function is not designed to print anything.
 */

int user_func_is_noprint (const ufunc *fun)
{
    if (fun == NULL) {
	return 0;
    } else {
	return function_is_noprint(fun);
    }
}

/**
 * user_function_name_by_index:
 * @i: the position of a user-function in the array of
 * loaded functions.
 * 
 * Returns: the name of the function, or %NULL if
 * @i is out of bounds.
 */

const char *user_function_name_by_index (int i)
{
    if (i >= 0 && i < n_ufuns) {
	return ufuns[i]->name;
    } else {
	return NULL;
    }
}

/**
 * user_function_index_by_name:
 * @name: function name.
 * @pkg: reference function package.
 * 
 * Returns: the 0-based position of a function of name
 * @name belonging to package @pkg, or -1 if there is
 * no such function.
 */

int user_function_index_by_name (const char *name,
				 fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg && 
	    !strcmp(name, ufuns[i]->name)) {
	    return i;
	}
    }

    return -1;
}	

/* Apparatus used in the GUI selector for composing a new function
   package.  We want a list of the names of currently unpackaged
   functions: we first call function_names_init(), then keep calling
   next_available_function_name() until it returns %NULL.  The pointer
   argument @idxp provides a means to grab the "index number"
   (position in the current functions array) corresponding to the
   returned function name.
*/

static int fname_idx;

void function_names_init (void)
{
    fname_idx = 0;
}

const char *next_available_function_name (int *idxp)
{
    const char *ret = NULL;
    ufunc *fun;

    if (n_ufuns == 0) {
	fname_idx = 0;
	return NULL;
    }

    while (fname_idx < n_ufuns) {
	fun = ufuns[fname_idx++];
	if (fun->pkg == NULL) {
	    ret = fun->name;
	    *idxp = fname_idx - 1;
	    break;
	}
    }

    return ret;
}

/* end GUI function-name apparatus */

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

    return (call != NULL)? call->fun : NULL;
}

/* see if a function is currently employed in the call stack */

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

/**
 * get_ufunc_by_name:
 * @name: name to test.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given name and it is accessible in
 * context (i.e. it's not private to a package other than
 * the one that's currently active, if any), otherwise
 * %NULL.
 */

ufunc *get_user_function_by_name (const char *name)
{
    fnpkg *pkg = NULL;
    ufunc *fun;
    int i;

    if (n_ufuns == 0) {
	return NULL;
    }

    fun = currently_called_function();

    if (fun != NULL) {
	pkg = fun->pkg;
	fun = NULL;
    }

    /* First pass: if there's no active function package, match any
       non-private function, but if there is an active package match
       only its member functions.  Try to optimize by putting the
       cheapest comparisons first.
    */

    for (i=0; i<n_ufuns; i++) {
	if ((pkg == NULL && !function_is_private(ufuns[i])) || ufuns[i]->pkg == pkg) {
	    if (!strcmp(name, ufuns[i]->name)) { 
		fun = ufuns[i];
		break;
	    }
	}
    }

    /* Second pass, if the first didn't work: match unpackaged
       functions or functions from other packages, so long as they're
       not private.
    */

    if (fun == NULL && pkg != NULL) {
	for (i=0; i<n_ufuns; i++) {
	    if (!function_is_private(ufuns[i]) && !strcmp(name, ufuns[i]->name)) {
		fun = ufuns[i];
		break;
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

/**
 * get_function_from_package:
 * @funname: name if fuction to retrieve.
 * @pkg: function package.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given @funname that is associated with
 * function package @pkg, otherwise NULL.  This is used
 * in the gretl function package editor.
 */

ufunc *get_function_from_package (const char *funname, fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg && 
	    !strcmp(funname, ufuns[i]->name)) {
	    return ufuns[i];
	}
    }

    return NULL;
}

/* allocate and initialize a new array of @n parameters */

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
	params[i].descrip = NULL;
	params[i].labels = NULL;
	params[i].nlabels = 0;
	params[i].flags = 0;
	params[i].deflt = NADBL;
	params[i].min = NADBL;
	params[i].max = NADBL;
	params[i].step = NADBL;
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
    fun->pkg = NULL;
    fun->flags = 0;
    fun->pkg_role = 0;

    fun->n_lines = 0;
    fun->lines = NULL;

    fun->n_params = 0;
    fun->params = NULL;

    fun->rettype = GRETL_TYPE_NONE;

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
	free_strings_array(params[i].labels, params[i].nlabels);
    }
    free(params);
}

static void clear_ufunc_data (ufunc *fun)
{
    free_strings_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    
    fun->lines = NULL;
    fun->params = NULL;

    fun->n_lines = 0;
    fun->n_params = 0;

    fun->rettype = GRETL_TYPE_NONE;
}

static void ufunc_free (ufunc *fun)
{
    free_strings_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);

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

#if FN_DEBUG
    fprintf(stderr, "add_ufunc: '%s'\n", fname);
#endif

    if (fun != NULL) {
	strncat(fun->name, fname, FN_NAMELEN - 1);
	if (add_allocated_ufunc(fun)) {
	    ufunc_free(fun);
	    fun = NULL;
	}
    }

    return fun;
}

/**
 * gretl_arg_type_name:
 * @type: a gretl type.
 *
 * Returns: the name of gretl type that is valid as a 
 * function argument or return value.
 */

const char *gretl_arg_type_name (GretlType type)
{
    switch (type) {
    case GRETL_TYPE_BOOL:       return "bool";
    case GRETL_TYPE_INT:        return "int";
    case GRETL_TYPE_OBS:        return "obs";
    case GRETL_TYPE_DOUBLE:     return "scalar";
    case GRETL_TYPE_SERIES:     return "series";
    case GRETL_TYPE_USERIES:    return "series";
    case GRETL_TYPE_MATRIX:     return "matrix";	
    case GRETL_TYPE_LIST:       return "list";
    case GRETL_TYPE_BUNDLE:     return "bundle";	
    case GRETL_TYPE_SCALAR_REF: return "scalar *";
    case GRETL_TYPE_SERIES_REF: return "series *";
    case GRETL_TYPE_MATRIX_REF: return "matrix *";
    case GRETL_TYPE_BUNDLE_REF: return "bundle *";
    case GRETL_TYPE_STRING:     return "string";
    case GRETL_TYPE_VOID:       return "void";
    case GRETL_TYPE_NONE:       return "null";
    default:
	return "invalid";
    }
}

/* handling of XML function packages */

enum {
    FUNCS_INFO,
    FUNCS_LOAD,
    FUNCS_CODE
};

static const char *arg_type_xml_string (int t)
{
    if (t == GRETL_TYPE_SCALAR_REF) {
	return "scalarref";
    } else if (t == GRETL_TYPE_SERIES_REF) {
	return "seriesref";
    } else if (t == GRETL_TYPE_MATRIX_REF) {
	return "matrixref";
    } else if (t == GRETL_TYPE_BUNDLE_REF) {
	return "bundleref";
    } else {
	return gretl_arg_type_name(t);
    }
}

GretlType gretl_type_from_string (const char *s)
{
    if (!strcmp(s, "bool"))     return GRETL_TYPE_BOOL;
    if (!strcmp(s, "boolean"))  return GRETL_TYPE_BOOL;
    if (!strcmp(s, "int"))      return GRETL_TYPE_INT;
    if (!strcmp(s, "obs"))      return GRETL_TYPE_OBS;
    if (!strcmp(s, "scalar"))   return GRETL_TYPE_DOUBLE;
    if (!strcmp(s, "series"))   return GRETL_TYPE_SERIES;
    if (!strcmp(s, "matrix"))   return GRETL_TYPE_MATRIX;
    if (!strcmp(s, "list"))     return GRETL_TYPE_LIST;
    if (!strcmp(s, "string"))   return GRETL_TYPE_STRING;
    if (!strcmp(s, "bundle"))   return GRETL_TYPE_BUNDLE;

    if (!strcmp(s, "scalar *"))  return GRETL_TYPE_SCALAR_REF;
    if (!strcmp(s, "series *"))  return GRETL_TYPE_SERIES_REF;
    if (!strcmp(s, "matrix *"))  return GRETL_TYPE_MATRIX_REF;
    if (!strcmp(s, "bundle *"))  return GRETL_TYPE_BUNDLE_REF;

    if (!strcmp(s, "scalarref"))  return GRETL_TYPE_SCALAR_REF;
    if (!strcmp(s, "seriesref"))  return GRETL_TYPE_SERIES_REF;
    if (!strcmp(s, "matrixref"))  return GRETL_TYPE_MATRIX_REF;
    if (!strcmp(s, "bundleref"))  return GRETL_TYPE_BUNDLE_REF;

    return 0;
}

static int return_type_from_string (const char *s)
{
    int t;

    if (!strcmp(s, "void")) {
	/* not OK as arg type, but OK as return */
	t = GRETL_TYPE_VOID;
    } else {
	t = gretl_type_from_string(s);
    }

    return (ok_function_return_type(t))? t : 0;
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
	/* backward compat */
	return arg_type_from_int(s);
    } else {
	return gretl_type_from_string(s);
    }
}    

/* read the parameter info for a function from XML file */

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
	    fn_param *param = &fun->params[n++];

	    if (gretl_xml_get_prop_as_string(cur, "name", &field)) {
		param->name = field;
	    } else {
		err = E_DATA;
		break;
	    }
	    if (gretl_xml_get_prop_as_string(cur, "type", &field)) {
		param->type = param_field_to_type(field);
		free(field);
		if (gretl_scalar_type(param->type)) {
		    gretl_xml_get_prop_as_double(cur, "default", 
						 &param->deflt);
		    if (param->type != GRETL_TYPE_BOOL) {
			gretl_xml_get_prop_as_double(cur, "min", &param->min);
			gretl_xml_get_prop_as_double(cur, "max", &param->max);
			gretl_xml_get_prop_as_double(cur, "step", &param->step);
		    }
		}
		if (gretl_xml_get_prop_as_bool(cur, "optional")) {
		    param->flags |= ARG_OPTIONAL;
		}
		if (gretl_xml_get_prop_as_bool(cur, "const")) {
		    param->flags |= ARG_CONST;
		}
	    } else {
		err = E_DATA;
		break;
	    }
	    gretl_xml_child_get_string(cur, doc, "description", 
				       &param->descrip);
	    gretl_xml_child_get_strings_array(cur, doc, "labels",
					      &param->labels,
					      &param->nlabels);
	}	    
	cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    if (!err && n != fun->n_params) {
	err = E_DATA;
    }

    return err;
}

/* we use this only with old-style function definitions */

static int func_read_return (xmlNodePtr node, ufunc *fun,
			     char **retname)
{
    char *field = NULL;
    int err = 0;

    if (gretl_xml_get_prop_as_string(node, "name", &field)) {
	*retname = field;
	field = NULL;
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

/* read the actual code lines from the XML representation of a 
   function */

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

static void print_param_labels (fn_param *param, PRN *prn)
{
    int i;

    pputs(prn, " {");

    for (i=0; i<param->nlabels; i++) {
	pprintf(prn, "\"%s\"", param->labels[i]);
	if (i < param->nlabels - 1) {
	    pputs(prn, ", ");
	}
    }

    pputc(prn, '}');
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

    /* minimum */
    if (!na(param->min)) pprintf(prn, "%g", param->min);
    pputc(prn, ':');

    /* maximum */
    if (!na(param->max)) pprintf(prn, "%g", param->max);
    pputc(prn, ':');

    /* default */
    if (!na(param->deflt)) pprintf(prn, "%g", param->deflt);

    if (!na(param->step)) {
	/* step */
	pputc(prn, ':');
	pprintf(prn, "%g", param->step);
    }

    pputc(prn, ']');
}

/* free @fun and also remove it from the list of loaded
   functions */

static void ufunc_unload (ufunc *fun)
{
    int i, j, found = 0;

    if (n_ufuns == 0 || fun == NULL) {
	/* "can't happen" */
	return;
    }

    /* remove this function from the array of loaded functions */

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i] == fun) {
	    for (j=i; j<n_ufuns-1; j++) {
		ufuns[j] = ufuns[j+1];
	    }
	    found = 1;
	    break;
	}
    }

    ufunc_free(fun);

    if (found) {
	n_ufuns--;
    }
}

/* remove a package from the current listing and free it; 
   if @full is non-zero, also unload any member functions
*/

static void real_function_package_unload (fnpkg *pkg, int full)
{
    int i, j, found = 0;

#if PKG_DEBUG
    fprintf(stderr, "function_package_unload: unloading '%s'\n", pkg->name);
#endif

    /* unload children? */

    if (full) {
	for (i=0; i<n_ufuns; i++) {
	    if (ufuns[i]->pkg == pkg) {
		ufunc_unload(ufuns[i]);
	    }
	}  
    }  

    /* remove this package from the array of loaded packages */

    for (i=0; i<n_pkgs; i++) {
	if (pkgs[i] == pkg) {
	    for (j=i; j<n_pkgs-1; j++) {
		pkgs[j] = pkgs[j+1];
	    }
	    found = 1;
	    break;
	}
    }

    /* free the package itself */
    function_package_free(pkg);

    if (found) {
	n_pkgs--;
    }
}

static void function_package_unload_full (fnpkg *pkg)
{
    real_function_package_unload(pkg, 1);
}

/* Append a pointer to @fun to the array of child-pointers in @pkg: we
   do this when reading the definition of a packaged function from an
   XML file.  Note that this action does not add the functions to the
   array of loaded functions -- that's done separately, if we're
   loading the package 'for real'.
*/

static int attach_ufunc_to_package (ufunc *fun, fnpkg *pkg)
{
    ufunc **uf;
    int n, err = 0;

    if (function_is_private(fun)) {
	n = pkg->n_priv;
	uf = realloc(pkg->priv, (n + 1) * sizeof *uf);
	if (uf == NULL) {
	    err = E_ALLOC;
	} else {
	    pkg->priv = uf;
	    pkg->priv[n] = fun;
	    pkg->n_priv += 1;
	}
    } else {
	n = pkg->n_pub;
	uf = realloc(pkg->pub, (n + 1) * sizeof *uf);
	if (uf == NULL) {
	    err = E_ALLOC;
	} else {
	    pkg->pub = uf;
	    pkg->pub[n] = fun;
	    pkg->n_pub += 1;
	}	
    } 

#if PKG_DEBUG
    fprintf(stderr, "attach_ufunc_to_package: package = %s, "
	    "private = %d, err = %d\n", pkg->name, 
	    function_is_private(fun), err);
#endif

    return err;
}

/* We got an old-style function definition from XML,
   in which the return type and name were recorded
   in the 'header'.  Having recorded the return type
   we now append a return statement.
*/

static int ufunc_add_return_statement (ufunc *fun,
				       const char *retname)
{
    char *s = malloc(8 + strlen(retname));
    int err = 0;

    if (s == NULL) {
	err = E_ALLOC;
    } else {
	sprintf(s, "return %s", retname);
	err = strings_array_add(&fun->lines, &fun->n_lines, s);
	free(s);
    }

    return err;
}

/* read a single user-function definition from XML file: if the
   function is a child of a package, the @pkg argument will
   be non-NULL
*/

static int read_ufunc_from_xml (xmlNodePtr node, xmlDocPtr doc, fnpkg *pkg)
{
    ufunc *fun = ufunc_new();
    xmlNodePtr cur;
    char *retname = NULL;
    char *tmp;
    int err = 0;

    if (fun == NULL) {
	return E_ALLOC;
    }

    if (!gretl_xml_get_prop_as_string(node, "name", &tmp)) {
	ufunc_free(fun);
	return E_DATA;
    }

    strncat(fun->name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (pkg != NULL) {
	fun->pkg = pkg;
    }

    if (gretl_xml_get_prop_as_string(node, "type", &tmp)) {
	fun->rettype = return_type_from_string(tmp);
	free(tmp);
    } else {
	fun->rettype = GRETL_TYPE_VOID;
    }

    if (gretl_xml_get_prop_as_bool(node, "private")) {
	fun->flags |= UFUN_PRIVATE;
    }

    if (gretl_xml_get_prop_as_bool(node, "plugin-wrapper")) {
	fun->flags |= UFUN_PLUGIN;
    }

    if (gretl_xml_get_prop_as_bool(node, "no-print")) {
	fun->flags |= UFUN_NOPRINT;
    }

    if (gretl_xml_get_prop_as_string(node, "pkg-role", &tmp)) {
	fun->pkg_role = pkg_key_get_role(tmp);
	free(tmp);
    }

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: name '%s', type %d\n"
	    " private = %d, plugin = %d\n", fun->name, fun->rettype, 
	    function_is_private(fun), function_is_plugin(fun));
#endif

    if (pkg == NULL && (function_is_private(fun) || function_is_plugin(fun))) {
	fprintf(stderr, "unpackaged function: can't be private or plugin\n");
	ufunc_free(fun);
	return E_DATA;
    }	

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "help")) {
	    /* backward compatibility: help used to be attached to
	       a package's public interface */
	    if (pkg->help != NULL) {
		free(pkg->help);
	    }
	    gretl_xml_node_get_string(cur, doc, &pkg->help);
	} else if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    err = func_read_params(cur, doc, fun);
	    if (err) {
		fprintf(stderr, "%s: error parsing function parameters\n",
			fun->name);
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "return")) {
	    /* support old-style functions */
	    err = func_read_return(cur, fun, &retname);
	    if (err) {
		fprintf(stderr, "%s: error parsing function return value\n",
			fun->name);
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "code")) {
	    err = func_read_code(cur, doc, fun);
	}
	cur = cur->next;
    }

    if (retname != NULL) {
	/* backward compat */
	if (!err) {
	    err = ufunc_add_return_statement(fun, retname);
	}
	free(retname);
    }

    if (!err) {
	if (pkg != NULL) {
	    /* function belongs to a package */
	    err = attach_ufunc_to_package(fun, pkg);
	} else {
	    /* reading a standalone function from session file */
	    err = add_allocated_ufunc(fun);
	}
    }

    if (err) {
	ufunc_free(fun);
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

    if (!strncmp(s, "catch ", 6)) {
	s += 6;
	s += strspn(s, " ");
    }

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
	shift_string_left(p + 3, 1);
    }
}

#define parm_has_children(p) (p->descrip != NULL || p->nlabels > 0)

/* write out a single user-defined function as XML, according to
   gretlfunc.dtd */

static int write_function_xml (ufunc *fun, FILE *fp)
{
    int rtype = fun->rettype;
    int this_indent = 0;
    int next_indent = 0;
    int i, j;

    if (rtype == GRETL_TYPE_NONE) {
	rtype = GRETL_TYPE_VOID;
    }

    fprintf(fp, "<gretl-function name=\"%s\" type=\"%s\"", 
	    fun->name, gretl_arg_type_name(rtype));

    if (function_is_private(fun)) {
	fputs(" private=\"1\"", fp);
    }    
    if (function_is_plugin(fun)) {
	fputs(" plugin-wrapper=\"1\"", fp);
    }
    if (function_is_noprint(fun)) {
	fputs(" no-print=\"1\"", fp);
    }

    if (fun->pkg_role) {
	fprintf(fp, " pkg-role=\"%s\"", pkg_role_get_key(fun->pkg_role));
    }

    fputs(">\n", fp);

    if (fun->n_params > 0) {

	gretl_push_c_numeric_locale();

	fprintf(fp, " <params count=\"%d\">\n", fun->n_params);
	for (i=0; i<fun->n_params; i++) {
	    fn_param *param = &fun->params[i];

	    fprintf(fp, "  <param name=\"%s\" type=\"%s\"",
		    param->name, arg_type_xml_string(param->type));
	    if (!na(param->min)) {
		fprintf(fp, " min=\"%g\"", param->min);
	    }
	    if (!na(param->max)) {
		fprintf(fp, " max=\"%g\"", param->max);
	    }
	    if (!na(param->deflt)) {
		fprintf(fp, " default=\"%g\"", param->deflt);
	    }
	    if (!na(param->step)) {
		fprintf(fp, " step=\"%g\"", param->step);
	    }
	    if (param->flags & ARG_OPTIONAL) {
		fputs(" optional=\"true\"", fp);
	    }
	    if (param->flags & ARG_CONST) {
		fputs(" const=\"true\"", fp);
	    }
	    if (parm_has_children(param)) {
		fputs(">\n", fp); /* terminate opening tag */
		if (param->descrip != NULL) {
		    gretl_xml_put_tagged_string("description", 
						param->descrip,
						fp);
		} 
		if (param->nlabels > 0) {
		    gretl_xml_put_strings_array_quoted("labels", 
						       (const char **) param->labels, 
						       param->nlabels, fp);
		}
		fputs("  </param>\n", fp);
	    } else {
		fputs("/>\n", fp); /* terminate opening tag */
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

    fputs("</code>\n", fp);

    fputs("</gretl-function>\n", fp);
    
    return 0;
}

/* script-style output */

static void print_function_start (ufunc *fun, PRN *prn)
{
    const char *s;
    int i;

    if (fun->rettype == GRETL_TYPE_NONE) {
	pprintf(prn, "function void %s ", fun->name);
    } else {
	const char *typestr = gretl_arg_type_name(fun->rettype);

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
	s = gretl_arg_type_name(fun->params[i].type);
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
	if (fun->params[i].nlabels > 0) {
	    print_param_labels(&fun->params[i], prn);
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

/**
 * gretl_function_print_code:
 * @u: pointer to user-function.
 * @prn: printing struct.
 *
 * Prints out function @fun to @prn, script-style.
 *
 * Returns: 0 on success, non-zero if @fun is %NULL.
 */

int gretl_function_print_code (ufunc *u, PRN *prn)
{
    int this_indent = 0;
    int next_indent = 0;
    int i, j;

    if (u == NULL) {
	return E_DATA;
    }
   
    print_function_start(u, prn);

    for (i=0; i<u->n_lines; i++) {
	adjust_indent(u->lines[i], &this_indent, &next_indent);
	for (j=0; j<=this_indent; j++) {
	    pputs(prn, "  ");
	}
	pputs(prn, u->lines[i]);
	pputc(prn, '\n');
    }

    pputs(prn, "end function\n");

    return 0;
}

/* construct a name for @pkg based on its filename member:
   take the basename and knock off ".gfn"
*/

static void name_package_from_filename (fnpkg *pkg)
{
    char *p = strrchr(pkg->fname, SLASH);
    int n;

    if (p != NULL) {
	p++;
    } else {
	p = pkg->fname;
    }

    n = strlen(p);
    if (has_suffix(p, ".gfn")) {
	n -= 4;
    }

    if (n > FN_NAMELEN - 1) {
	n = FN_NAMELEN - 1;
    }

    *pkg->name = '\0';
    strncat(pkg->name, p, n);

#if PKG_DEBUG
    fprintf(stderr, "filename '%s' -> pkgname '%s'\n",
	    pkg->fname, pkg->name);
#endif
}

/* name lookup for functions to be connected to a package,
   allowing for the possibility that they're already
   connected */

static ufunc *get_uf_array_member (const char *name, fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg || ufuns[i]->pkg == NULL) {
	    if (!strcmp(name, ufuns[i]->name)) {
		return ufuns[i];
	    }
	}
    }   

    return NULL;
}

static void check_special_comments (ufunc *fun)
{
    int i;

    for (i=0; i<fun->n_lines; i++) {
	if (strstr(fun->lines[i], "## plugin-wrapper ##")) {
	    fun->flags |= UFUN_PLUGIN;
	} else if (strstr(fun->lines[i], "## no-print ##")) {
	    fun->flags |= UFUN_NOPRINT;
	}
    }
}

/* Given an array of @n function names in @names, set up a
   corresponding array of pointers to the named functions in @pkg,
   Flag an error if any of the function names are bad, or if
   allocation fails.  If we're revising an existing function package
   we start by disconnecting the currently connected functions.  
*/

static int set_uf_array_from_names (fnpkg *pkg, char **names, 
				    int n, int priv)
{
    ufunc **uf;
    ufunc *fun;
    int i;

    /* error-checking first */
    
    for (i=0; i<n; i++) {
	if (get_uf_array_member(names[i], pkg) == NULL) {
	    return E_DATA;
	}
    }    

    /* then disconnect current functions, if any */

    if (priv) {
	for (i=0; i<pkg->n_priv; i++) {
	    pkg->priv[i]->pkg = NULL;
	    set_function_private(pkg->priv[i], FALSE);
	}
    } else {
	for (i=0; i<pkg->n_pub; i++) {
	    pkg->pub[i]->pkg = NULL;
	}
    }	

    /* then connect the specified functions */

    if (priv) {
	uf = realloc(pkg->priv, n * sizeof *uf);
    } else {
	uf = realloc(pkg->pub, n * sizeof *uf);
    } 

    if (uf == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	fun = get_uf_array_member(names[i], NULL);
	fun->pkg = pkg;
	set_function_private(fun, priv);
	check_special_comments(fun);
	uf[i] = fun;
    }

    if (priv) {
	pkg->priv = uf;
	pkg->n_priv = n;
    } else {
	pkg->pub = uf;
	pkg->n_pub = n;
    }

    return 0;
}

/**
 * function_package_connect_funcs:
 * @pkg: function package.
 * @pubnames: array of names of public functions.
 * @n_pub: number of strings in @pubnames.
 * @privnames: array of names of private functions (or %NULL).
 * @n_priv: number of strings in @privnames (may be 0).
 *
 * Looks up the functions named in @pubnames and @privnames
 * and adds pointers to these functions to @pkg, hence marking
 * the functions as belonging to @pkg.
 *
 * Returns: 0 on success, non-zero on error.
 */

int function_package_connect_funcs (fnpkg *pkg, 
				    char **pubnames, int n_pub,
				    char **privnames, int n_priv) 
{
    int err = 0;

    if (n_pub > 0) {
	err = set_uf_array_from_names(pkg, pubnames, n_pub, 0);
    }

    if (!err && n_priv > 0) {
	err = set_uf_array_from_names(pkg, privnames, n_priv, 1);
    }

    return err;
}

/**
 * function_package_new:
 * @fname: filename for package.
 * @pubnames: array of names of public functions.
 * @n_pub: number of strings in @pubnames.
 * @privnames: array of names of private functions (or %NULL).
 * @n_priv: number of strings in @privnames (may be 0).
 * @err: location to receive error code.
 *
 * Allocates a new package with filename-member @fname, including
 * the public and private functions named in @pubnames and
 * @privnames.  Note that this function does not cause the
 * package to be written to file; for that, see
 * write_function_package().
 *
 * Returns: pointer to package on success, %NULL on error.
 */

fnpkg *function_package_new (const char *fname, 
			     char **pubnames, int n_pub,
			     char **privnames, int n_priv, 
			     int *err)
{
    fnpkg *pkg = function_package_alloc(fname);

    if (pkg == NULL) {
	*err = E_ALLOC;
	return NULL;
    } 

    name_package_from_filename(pkg);

    *err = function_package_connect_funcs(pkg, pubnames, n_pub,
					  privnames, n_priv);

    if (!*err) {
	*err = function_package_record(pkg);
    }

    if (*err) {
	/* note: this does not free the packaged functions, if any */
	function_package_free(pkg);
	pkg = NULL;
    }

    return pkg;
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

static char *trim_text (char *s)
{
    while (isspace(*s)) s++;
    return tailstrip(s);
}

static int package_write_translatable_strings (fnpkg *pkg, PRN *prn)
{
    FILE *fp;
    gchar *trname;
    char **S = NULL;
    int i, n = 0;

    trname = g_strdup_printf("%s-i18n.c", pkg->name);

    fp = gretl_fopen(trname, "w");
    if (fp == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), trname);
	g_free(trname);
	return E_FOPEN;
    }

    if (pkg->pub != NULL) {
	int j, k;

	for (i=0; i<pkg->n_pub; i++) {
	    ufunc *fun = pkg->pub[i];

	    for (j=0; j<fun->n_params; j++) {
		fn_param *param = &fun->params[j];

		if (param->descrip != NULL) {
		    strings_array_add(&S, &n, param->descrip);
		}
		for (k=0; k<param->nlabels; k++) {
		    strings_array_add(&S, &n, param->labels[k]);
		}
	    }
	}
    }

    if (pkg->label != NULL || S != NULL) {
	fprintf(fp, "const char *%s_translations[] = {\n", pkg->name);
	if (pkg->label != NULL) {
	    fprintf(fp, "    N_(\"%s\"),\n", pkg->label);
	}	
	if (S != NULL) {
	    strings_array_sort(&S, &n, OPT_U);
	    for (i=0; i<n; i++) {
		fprintf(fp, "    N_(\"%s\")", S[i]);
		if (i < n-1) {
		    fputc(',', fp);
		}
		fputc('\n', fp);
	    }
	    free_strings_array(S, n);
	}
	fputs("};\n", fp);
    }

    fclose(fp);

    pprintf(prn, "Wrote translations file %s\n", trname);
    g_free(trname);

    return 0;
}

static int package_write_index (fnpkg *pkg, PRN *prn)
{
    gchar *idxname;
    FILE *fp;

    idxname = g_strdup_printf("%s.xml", pkg->name);

    fp = gretl_fopen(idxname, "w");
    if (fp == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), idxname);
	g_free(idxname);
	return E_FOPEN;
    }

    gretl_xml_header(fp);

    fputs("<gretl-addon", fp);

    fprintf(fp, " name=\"%s\"", pkg->name);

    if (pkg->dreq == FN_NEEDS_TS) {
	fprintf(fp, " %s=\"true\"", NEEDS_TS);
    } else if (pkg->dreq == FN_NEEDS_QM) {
	fprintf(fp, " %s=\"true\"", NEEDS_QM);
    } else if (pkg->dreq == FN_NEEDS_PANEL) {
	fprintf(fp, " %s=\"true\"", NEEDS_PANEL);
    } else if (pkg->dreq == FN_NODATA_OK) {
	fprintf(fp, " %s=\"true\"", NO_DATA_OK);
    }

    if (pkg->minver > 0) {
	char vstr[8];

	fprintf(fp, " minver=\"%s\"", get_version_string(vstr, pkg->minver));
    }

    fputs(">\n", fp);

    gretl_xml_put_tagged_string("author",  pkg->author, fp);
    gretl_xml_put_tagged_string("version", pkg->version, fp);
    gretl_xml_put_tagged_string("date",    pkg->date, fp);
    gretl_xml_put_tagged_string("description", pkg->descrip, fp);

    if (pkg->label != NULL) {
	gretl_xml_put_tagged_string("label", pkg->label, fp);
    }

    fputs("</gretl-addon>\n", fp);

    fclose(fp);

    pprintf(prn, "Wrote index file %s\n", idxname);
    g_free(idxname);

    return 0;
}

/* Write out a function package as XML. If @fp is NULL, then we write
   this as a package file in its own right, using the filename given
   in the package (the 'standalone' case), otherwise we're writing it
   as a component of the functions listing in a gretl session file,
   which already has an associated open stream.  
*/

static int real_write_function_package (fnpkg *pkg, FILE *fp)
{
    int standalone = (fp == NULL);
    int i, err = 0;

    if (standalone) {
	fp = gretl_fopen(pkg->fname, "w");
	if (fp == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't open %s"), pkg->fname);
	    return E_FOPEN;
	} else {
	    gretl_xml_header(fp);
	    fputs("<gretl-functions>\n", fp);
	}
    }

    fputs("<gretl-function-package", fp);

    name_package_from_filename(pkg); /* ? */

    fprintf(fp, " name=\"%s\"", pkg->name);
    fprintf(fp, " ID=\"%d\"", pkg->ID);

    if (pkg->dreq == FN_NEEDS_TS) {
	fprintf(fp, " %s=\"true\"", NEEDS_TS);
    } else if (pkg->dreq == FN_NEEDS_QM) {
	fprintf(fp, " %s=\"true\"", NEEDS_QM);
    } else if (pkg->dreq == FN_NEEDS_PANEL) {
	fprintf(fp, " %s=\"true\"", NEEDS_PANEL);
    } else if (pkg->dreq == FN_NODATA_OK) {
	fprintf(fp, " %s=\"true\"", NO_DATA_OK);
    }

    if (pkg->minver > 0) {
	char vstr[8];

	fprintf(fp, " minver=\"%s\"", get_version_string(vstr, pkg->minver));
    }

    fputs(">\n", fp);

    gretl_xml_put_tagged_string("author",  pkg->author, fp);
    gretl_xml_put_tagged_string("version", pkg->version, fp);
    gretl_xml_put_tagged_string("date",    pkg->date, fp);
    gretl_xml_put_tagged_string("description", pkg->descrip, fp);

    if (pkg->label != NULL) {
	gretl_xml_put_tagged_string("label", pkg->label, fp);
    }

    if (pkg->help != NULL) {
	fputs("<help>\n", fp);
	gretl_xml_put_raw_string(trim_text(pkg->help), fp);
	fputs("\n</help>\n", fp);
    } 

    if (pkg->gui_help != NULL) {
	fputs("<gui-help>\n", fp);
	gretl_xml_put_raw_string(trim_text(pkg->gui_help), fp);
	fputs("\n</gui-help>\n", fp);
    }     

    if (pkg->pub != NULL) {
	for (i=0; i<pkg->n_pub; i++) {
	    write_function_xml(pkg->pub[i], fp);
	}
    }    

    if (pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv; i++) {
	    write_function_xml(pkg->priv[i], fp);
	}
    }

    if (pkg->sample != NULL) {
	fputs("<sample-script>\n", fp);
	gretl_xml_put_raw_string(trim_text(pkg->sample), fp);
	fputs("\n</sample-script>\n", fp);	
    }

    fputs("</gretl-function-package>\n", fp);

    if (standalone) {
	fputs("</gretl-functions>\n", fp);
	fclose(fp);
    }

    return err;
}

/**
 * function_package_write_file:
 * @pkg: function package.
 *
 * Write out @pkg as an XML file, using the filename
 * recorded in the package.
 *
 * Returns: 0 on success, non-zero on error.
 */

int function_package_write_file (fnpkg *pkg)
{
    return real_write_function_package(pkg, NULL);
}

/* below: apparatus for constructing and saving a gfn function
   package from the command line
*/

static int is_unclaimed (const char *s, char **S, int n)
{
    int i;

    for (i=0; i<n; i++) {
	if (strcmp(s, S[i]) == 0) {
	    return 0;
	}
    }

    return 1;
}

static gchar *pkg_aux_content (const char *fname, int *err)
{
    gchar *ret = NULL;
    GError *gerr = NULL;
    gsize len = 0;

    g_file_get_contents(fname, &ret, &len, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
	*err = E_FOPEN;
    }

    return ret;
}

static int pkg_set_dreq (fnpkg *pkg, const char *s)
{
    int err = 0;

    if (!strcmp(s, NEEDS_TS)) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (!strcmp(s, NEEDS_QM)) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (!strcmp(s, NEEDS_PANEL)) {
	pkg->dreq = FN_NEEDS_PANEL;
    } else if (!strcmp(s, NO_DATA_OK)) {
	pkg->dreq = FN_NODATA_OK;
    } else {
	err = E_PARSE;
    }

    return err;
}

static int function_set_pkg_role (const char *name, fnpkg *pkg,
				  const char *attr, PRN *prn)
{
    ufunc *u;
    int role = pkg_key_get_role(attr);
    int i, j, err = 0;

    /* check that the function in question satisfies the
       requirements for its role, and if so, hook it up 
    */

    if (role == UFUN_GUI_PRECHECK) {
	for (i=0; i<pkg->n_priv; i++) {
	    if (!strcmp(name, pkg->priv[i]->name)) {
		u = pkg->priv[i];
		if (u->rettype != GRETL_TYPE_DOUBLE) {
		    pprintf(prn, "%s: must return a scalar\n", attr);
		    err = E_TYPES;
		} else if (u->n_params > 0) {
		    pprintf(prn, "%s: no parameters are allowed\n", attr);
		    err = E_TYPES;
		}
		if (!err) {
		    u->pkg_role = role;
		}		
		return err;
	    }
	}
	pprintf(prn, "%s: %s: no such private function\n", attr, name);
	return E_DATA;
    }

    for (i=0; i<pkg->n_pub; i++) {
	if (!strcmp(name, pkg->pub[i]->name)) {
	    u = pkg->pub[i];
	    if (0 && u->n_params == 0) {
		/* void functions can't be right? */
		pprintf(prn, "%s: is marked as public but is void\n", pkg->pub[i]->name);
		err = E_TYPES;
	    } else if (role == UFUN_GUI_MAIN) {
		if (u->rettype != GRETL_TYPE_BUNDLE && u->rettype != GRETL_TYPE_VOID) {
		    pprintf(prn, "%s: must return a bundle, or nothing\n", attr);
		    err = E_TYPES;
		}
	    } else {
		/* bundle-print, bundle-plot, etc. */
		for (j=0; j<u->n_params && !err; j++) {
		    if (j == 0 && u->params[j].type != GRETL_TYPE_BUNDLE_REF) {
			pprintf(prn, "%s: first param type must be %s\n",
				attr, gretl_arg_type_name(GRETL_TYPE_BUNDLE_REF));
			err = E_TYPES;
		    } else if (j == 1 && u->params[j].type != GRETL_TYPE_INT) {
			pprintf(prn, "%s: second param type must be %s\n",
				attr, gretl_arg_type_name(GRETL_TYPE_INT));
			err = E_TYPES;
		    } else if (j > 1 && !fn_param_optional(u, j) &&
			       na(fn_param_default(u, j))) {
			pprintf(prn, "%s: params beyond the second must be optional\n",
				attr);
			err = E_TYPES;
		    }
		}
	    }
	    if (!err) {
		u->pkg_role = role;
	    }
	    return err;
	}
    }

    pprintf(prn, "%s: %s: no such public function\n", attr, name);

    return E_DATA;
}

/* having assembled and checked the function-listing for a new
   package, now retrieve the additional information from the
   spec file
*/

static int new_package_info_from_spec (fnpkg *pkg, FILE *fp, PRN *prn)
{
    char *p, line[1024];
    gchar *tmp;
    int got = 0;
    int err = 0;

    while (fgets(line, sizeof line, fp) && !err) {
	if (*line == '#' || string_is_blank(line)) {
	    continue;
	}
	tailstrip(line);
	p = strchr(line, '=');
	if (p == NULL) {
	    continue;
	} else {
	    p++;
	    p += strspn(p, " ");
	    if (!strncmp(line, "author", 6)) {
		err = function_package_set_properties(pkg, "author", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "version", 7)) {
		err = function_package_set_properties(pkg, "version", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "date", 4)) {
		err = function_package_set_properties(pkg, "date", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "description", 11)) {
		err = function_package_set_properties(pkg, "description", p, NULL);
		if (!err) got++;
	    } else if (!strncmp(line, "label", 5)) {
		err = function_package_set_properties(pkg, "label", p, NULL);
	    } else if (!strncmp(line, "help", 4)) {
		if (has_suffix(p, ".pdf")) {
		    pprintf(prn, "Recording help reference %s\n", p);
		    tmp = g_strdup_printf("pdfdoc:%s", p);
		} else {
		    pprintf(prn, "Looking for help text in %s\n", p);
		    tmp = pkg_aux_content(p, &err);
		}
		if (!err) {
		    err = function_package_set_properties(pkg, "help", tmp, NULL);
		    if (!err) got++;
		    g_free(tmp);
		}
	    } else if (!strncmp(line, "gui-help", 8)) {
		pprintf(prn, "Looking for GUI help text in %s\n", p);
		tmp = pkg_aux_content(p, &err);
		if (!err) {
		    err = function_package_set_properties(pkg, "gui-help", tmp, NULL);
		    g_free(tmp);
		}
	    } else if (!strncmp(line, "sample-script", 13)) {
		pprintf(prn, "Looking for sample script in %s\n", p);
		tmp = pkg_aux_content(p, &err);
		if (!err) {
		    err = function_package_set_properties(pkg, "sample-script", tmp, NULL);
		    if (!err) got++;
		    g_free(tmp);
		}
	    } else if (!strncmp(line, "data-requirement", 16)) {
		err = pkg_set_dreq(pkg, p);
	    } else if (!strncmp(line, "min-version", 11)) {
		pkg->minver = version_number_from_string(p);
		got++;
	    } else {
		const char *key;
		int i;

		for (i=0; pkg_lookups[i].key != NULL; i++) {
		    key = pkg_lookups[i].key;
		    if (!strncmp(line, key, strlen(key))) {
			err = function_set_pkg_role(p, pkg, key, prn);
			if (!err) {
			    pprintf(prn, "%s function is %s, OK\n", key, p);
			}
			break;
		    }
		}
	    }
	}
    }

    if (!err && got < 7) {
	gretl_errmsg_set("Some required information was missing");
	err = E_DATA;
    }

    return err;
}

static fnpkg *new_pkg_from_spec_file (const char *gfnname, PRN *prn,
				      int *err)
{
    fnpkg *pkg = NULL;
    char *p, fname[FILENAME_MAX];
    char line[4096], cont[1024];
    FILE *fp;

    if (!has_suffix(gfnname, ".gfn")) {
	gretl_errmsg_set("Output must have extension \".gfn\"");
	*err = E_DATA;
	return NULL;
    }

    strcpy(fname, gfnname);
    p = strrchr(fname, '.');
    strcpy(p, ".spec");

    fp = gretl_fopen(fname, "r");

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	char **pubnames = NULL;
	char **privnames = NULL;
	int npub = 0, npriv = 0;
	ufunc *uf;
	int i;

	pprintf(prn, "Found spec file '%s'\n", fname);

	/* first pass: gather names of public functions */

	while (fgets(line, sizeof line, fp) && !*err) {
	    if (!strncmp(line, "public =", 8)) {
		while (ends_with_backslash(line)) {
		    charsub(line, '\\', '\0');
		    *cont = '\0';
		    fgets(cont, sizeof cont, fp);
		    strcat(line, cont);
		}
		tailstrip(line);
		pubnames = gretl_string_split(line + 8, &npub);
		if (npub == 0) {
		    *err = E_DATA;
		}
	    } 
	}

	if (!*err) {
	    pprintf(prn, "number of public interfaces = %d\n", npub);
	    for (i=0; i<npub && !*err; i++) {
		uf = get_user_function_by_name(pubnames[i]);
		pprintf(prn, " %s", pubnames[i]);
		if (uf == NULL) {
		    pputs(prn, ": *** not found");
		    gretl_errmsg_sprintf("'%s': no such function", pubnames[i]);
		    *err = E_DATA;
		} 
		pputc(prn, '\n');
	    }
	}

	/* note: any other functions currently loaded are assumed to be
	   private functions for this package */

	if (!*err) {
	    npriv = n_free_functions() - npub;
	    if (npriv < 0) {
		npriv = 0;
	    }
	}	

	if (!*err && npriv > 0) {
	    npriv = 0;
	    for (i=0; i<n_ufuns && !*err; i++) {
		if (ufuns[i]->pkg == NULL && is_unclaimed(ufuns[i]->name, pubnames, npub)) {
		    *err = strings_array_add(&privnames, &npriv, ufuns[i]->name);
		}
	    }
	}

	if (!*err && npriv > 0) {
	    pprintf(prn, "number of private functions = %d\n", npriv);
	    for (i=0; i<npriv; i++) {
		pprintf(prn, " %s\n", privnames[i]);
	    }
	}

	if (!*err) {
	    pkg = function_package_new(gfnname, pubnames, npub,
				       privnames, npriv, err);
	}

	free_strings_array(pubnames, npub);
	free_strings_array(privnames, npriv);

	if (!*err) {
	    rewind(fp);
	    *err = new_package_info_from_spec(pkg, fp, prn);
	}

	if (*err && pkg != NULL) {
	    real_function_package_unload(pkg, 0);
	    pkg = NULL;
	}

	fclose(fp);
    }

    return pkg;
}

static int cli_validate_package_file (const char *fname, PRN *prn)
{
    char dtdname[FILENAME_MAX];
    xmlDocPtr doc;
    xmlDtdPtr dtd;
    int err = 0;

    doc = xmlParseFile(fname);
    if (doc == NULL) {
	pprintf(prn, "Couldn't parse %s\n", fname);
	return 1;
    }

    sprintf(dtdname, "%sfunctions%cgretlfunc.dtd", gretl_home(), SLASH);
    dtd = xmlParseDTD(NULL, (const xmlChar *) dtdname); 

    if (dtd == NULL) {
	pputs(prn, "Couldn't open DTD to check package\n");
    } else {
	const char *pkgname = path_last_element(fname);
	xmlValidCtxtPtr cvp = xmlNewValidCtxt();

	if (cvp == NULL) {
	    pputs(prn, "Couldn't get an XML validation context\n");
	    xmlFreeDtd(dtd);
	    xmlFreeDoc(doc);
	    return 0;
	}

	cvp->userData = (void *) prn;
	cvp->error    = (xmlValidityErrorFunc) pprintf;
	cvp->warning  = (xmlValidityWarningFunc) pprintf;

	pprintf(prn, "Checking against %s\n", dtdname);

	if (!xmlValidateDtd(cvp, doc, dtd)) {
	    err = 1;
	} else {
	    pprintf(prn, _("%s: validated against DTD OK"), pkgname);
	    pputc(prn, '\n');
	}

	xmlFreeValidCtxt(cvp);
	xmlFreeDtd(dtd);
    } 

    xmlFreeDoc(doc);

    return err;
}

/**
 * create_and_write_function_package:
 * @fname: filename for function package.
 * @opt: may include OPT_I to write a package-index entry,
 * OPT_T to write translatable strings.
 * @prn: printer struct for feedback.
 *
 * Create a package based on the functions currently loaded, and
 * write it out as an XML file. Responds to the makepkg command.
 *
 * Returns: 0 on success, non-zero on error.
 */

int create_and_write_function_package (const char *fname, 
				       gretlopt opt,
				       PRN *prn)
{
    fnpkg *pkg = NULL;
    int err = 0;

    if (n_free_functions() == 0) {
	gretl_errmsg_set(_("No functions are available for packaging at present."));
	err = E_DATA;
    } else {
	pkg = new_pkg_from_spec_file(fname, prn, &err);
	if (pkg != NULL) {
	    err = function_package_write_file(pkg);
	    if (!err) {
		err = cli_validate_package_file(fname, prn);
		/* should we delete @fname ? */
	    }
	    if (!err) {
		if (opt & OPT_T) {
		    package_write_translatable_strings(pkg, prn);
		}
		if (opt & OPT_I) {
		    package_write_index(pkg, prn);
		}		
	    }
	}
    }

    return err;
}

/**
 * function_package_get_name:
 * @pkg: function package.
 *
 * Returns: the name of the package.
 */

const char *function_package_get_name (fnpkg *pkg)
{
    return (pkg != NULL)? pkg->name : NULL;
}

static int maybe_replace_string_var (char **svar, const char *src)
{
    if (*svar == NULL) {
	*svar = gretl_strdup(src);
    } else if (strcmp(*svar, src)) {
	free(*svar);
	*svar = gretl_strdup(src);
    }

    return (*svar == NULL)? E_ALLOC : 0;
}

static int is_string_property (const char *key)
{
    return !strcmp(key, "fname") ||
	!strcmp(key, "author")   ||
	!strcmp(key, "version")  ||
	!strcmp(key, "date")     ||
	!strcmp(key, "description") ||
	!strcmp(key, "label") ||
	!strcmp(key, "help") ||
	!strcmp(key, "gui-help") ||
	!strcmp(key, "sample-script");
}

/* varargs function for setting the properties of a function
   package: the settings take the form of a NULL-terminated
   set of (key, value) pairs.
*/

int function_package_set_properties (fnpkg *pkg, ...)
{
    va_list ap;
    const char *key;
    int i, err = 0;

    va_start(ap, pkg);

    for (i=1; !err; i++) {
	key = va_arg(ap, const char *);
	if (key == NULL) {
	    break;
	}
	if (is_string_property(key)) {
	    const char *cval = va_arg(ap, const char *);

	    if (!strcmp(key, "fname")) {
		err = maybe_replace_string_var(&pkg->fname, cval);
	    } else if (!strcmp(key, "author")) {
		err = maybe_replace_string_var(&pkg->author, cval);
	    } else if (!strcmp(key, "date")) {
		err = maybe_replace_string_var(&pkg->date, cval);
	    } else if (!strcmp(key, "version")) {
		err = maybe_replace_string_var(&pkg->version, cval);
	    } else if (!strcmp(key, "description")) {
		err = maybe_replace_string_var(&pkg->descrip, cval);
	    } else if (!strcmp(key, "label")) {
		err = maybe_replace_string_var(&pkg->label, cval);
	    } else if (!strcmp(key, "help")) {
		err = maybe_replace_string_var(&pkg->help, cval);
	    } else if (!strcmp(key, "gui-help")) {
		err = maybe_replace_string_var(&pkg->gui_help, cval);
	    } else if (!strcmp(key, "sample-script")) {
		err = maybe_replace_string_var(&pkg->sample, cval);
	    } 
	} else {
	    int ival = va_arg(ap, int);

	    if (!strcmp(key, "data-requirement")) {
		pkg->dreq = ival;
	    } else if (!strcmp(key, "min-version")) {
		pkg->minver = ival;
	    } 
	} 
    }

    va_end(ap);

    return err;
}

enum {
    PUBLIST,
    GUILIST,
    PRIVLIST
};

/* From a function package get a list of either its public or its
   private functions, in the form of a simple gretl list.  The
   elements of this list are just the positions of these functions in
   the current recorder array for user-functions.  This is fine if
   the information is to be used right away (as in the GUI function
   call dialog), but it should _not_ be assumed that the identifiers
   will remain valid indefinitely.
*/

static int *function_package_get_list (fnpkg *pkg, int code, int n)
{
    int *list = NULL;
    int subtract = 0;

    if (n > 0) {
 	list = gretl_list_new(n);
	if (list != NULL) {
	    int i, j = 1;

	    for (i=0; i<n_ufuns; i++) {
		if (ufuns[i]->pkg == pkg) {
		    int priv = function_is_private(ufuns[i]);

		    if (code == PRIVLIST && priv) {
			list[j++] = i;
		    } else if (code == PUBLIST && !priv) {
			list[j++] = i;
		    } else if (code == GUILIST && !priv) {
			if (ufuns[i]->pkg_role == UFUN_BUNDLE_PRINT) {
			    subtract = 1;
			} else {
			    list[j++] = i;
			}
		    }
		}
	    }	    
	} 
    }

    if (list != NULL && subtract) {
	list[0] -= 1;
	if (list[0] == 0) {
	    free(list);
	    list = NULL;
	}
    }

    return list;
}

static char *pkg_get_special_func (fnpkg *pkg, int role)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == role) {
	    return g_strdup(ufuns[i]->name);
	}
    }	    

    return NULL;
}

/* varargs function for retrieving the properties of a function
   package: the arguments after @pkg take the form of a
   NULL-terminated set of (key, pointer) pairs; values are written to
   the locations given by the pointers.
*/

int function_package_get_properties (fnpkg *pkg, ...)
{
    va_list ap;
    int npub = 0;
    int npriv = 0;
    int **plist;
    const char *key;
    void *ptr;
    char **ps;
    int *pi;
    int i, err = 0;

    g_return_val_if_fail(pkg != NULL, E_DATA);

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == pkg) {
	    if (function_is_private(ufuns[i])) {
		npriv++;
	    } else {
		npub++;
	    }
	}
    }

    va_start(ap, pkg);

    for (i=1; ; i++) {
	key = va_arg(ap, const char *);
	if (key == NULL) {
	    break;
	}
	ptr = va_arg(ap, void *);
	if (ptr == NULL) {
	    break;
	}
	if (!strcmp(key, "name")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->name);
	} else if (!strcmp(key, "author")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->author);
	} else if (!strcmp(key, "date")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->date);
	} else if (!strcmp(key, "version")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->version);
	} else if (!strcmp(key, "description")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->descrip);
	} else if (!strcmp(key, "help")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->help);
	} else if (!strcmp(key, "gui-help")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->gui_help);
	} else if (!strcmp(key, "sample-script")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->sample);
	} else if (!strcmp(key, "label")) {
	    ps = (char **) ptr;
	    *ps = g_strdup(pkg->label);
	} else if (!strcmp(key, "data-requirement")) {
	    pi = (int *) ptr;
	    *pi = pkg->dreq;
	} else if (!strcmp(key, "min-version")) {
	    pi = (int *) ptr;
	    *pi = pkg->minver;
	} else if (!strcmp(key, "publist")) {
	    plist = (int **) ptr;
	    *plist = function_package_get_list(pkg, PUBLIST, npub);
	} else if (!strcmp(key, "gui-publist")) {
	    plist = (int **) ptr;
	    *plist = function_package_get_list(pkg, GUILIST, npub);
	} else if (!strcmp(key, "privlist")) {
	    plist = (int **) ptr;
	    *plist = function_package_get_list(pkg, PRIVLIST, npriv);
	} else if (!strcmp(key, BUNDLE_PRINT)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_PRINT);
	} else if (!strcmp(key, BUNDLE_PLOT)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_PLOT);
	} else if (!strcmp(key, BUNDLE_TEST)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_TEST);
	} else if (!strcmp(key, BUNDLE_FCAST)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_FCAST);
	} else if (!strcmp(key, BUNDLE_EXTRA)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_BUNDLE_EXTRA);
	} else if (!strcmp(key, GUI_MAIN)) {	
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_GUI_MAIN);
	} else if (!strcmp(key, GUI_PRECHECK)) {
	    ps = (char **) ptr;
	    *ps = pkg_get_special_func(pkg, UFUN_GUI_PRECHECK);
	}	    
    }

    va_end(ap);

    return err;
}

/* quick check to see if there's a gross problem with a package */

static int validate_function_package (fnpkg *pkg)
{
    if (pkg->pub == NULL || pkg->author == NULL ||
	pkg->version == NULL || pkg->date == NULL ||
	pkg->descrip == NULL) {
	return 0;
    }

    return 1;
}

/* for session file use: dump all currently defined functions,
   packaged or not, into a single XML file */

int write_session_functions_file (const char *fname)
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

    /* first write any loaded function packages */

    for (i=0; i<n_pkgs; i++) {
	if (validate_function_package(pkgs[i])) {
	    real_write_function_package(pkgs[i], fp);
	}
    }

    /* then any unpackaged functions */

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkg == NULL) {
	    write_function_xml(ufuns[i], fp);
	}
    }

    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    return 0;
}

/* De-allocate a function package: this can be done in either of two
   modes.  If 'full' is non-zero then we destroy all functions that
   are children of the given package, otherwise we leave the
   function-children alone (just 'detaching' them from the parent
   package).
*/

static void real_function_package_free (fnpkg *pkg, int full)
{
    if (pkg != NULL) {
	int i;

	if (full) {
	    for (i=0; i<pkg->n_pub; i++) {
		ufunc_free(pkg->pub[i]);
	    }
	    for (i=0; i<pkg->n_priv; i++) {
		ufunc_free(pkg->priv[i]);
	    }
	} else {
	    for (i=0; i<n_ufuns; i++) {
		if (ufuns[i]->pkg == pkg) {
		    /* remove package info */
		    ufuns[i]->pkg = NULL;
		    set_function_private(ufuns[i], FALSE);
		}
	    }
	}

	free(pkg->pub);
	free(pkg->priv);
	free(pkg->fname);
	free(pkg->author);
	free(pkg->version);
	free(pkg->date);
	free(pkg->descrip);
	free(pkg->help);
	free(pkg->gui_help);
	free(pkg->sample);
	free(pkg->label);
	free(pkg);
    }
}

static void function_package_free (fnpkg *pkg)
{
    real_function_package_free(pkg, 0);
}

static void function_package_free_full (fnpkg *pkg)
{
    real_function_package_free(pkg, 1);
}

/* is the package with filename @fname already in memory? */

static fnpkg *get_loaded_pkg_by_filename (const char *fname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    return pkgs[i];
	}
    }

    return NULL;
}

/**
 * function_package_unload_by_filename:
 * @fname: package filename.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are not destroyed; they become available for inclusion in
 * other packages.
 */

void function_package_unload_by_filename (const char *fname)
{
    fnpkg *pkg = get_loaded_pkg_by_filename(fname);

    if (pkg != NULL) {
	real_function_package_unload(pkg, 0);
    }
}

/**
 * function_package_unload_full_by_filename:
 * @fname: package filename.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are also unloaded from memory.
 */

void function_package_unload_full_by_filename (const char *fname)
{
    fnpkg *pkg = get_loaded_pkg_by_filename(fname);

    if (pkg != NULL) {
	real_function_package_unload(pkg, 1);
    }
}

/* append a function package to the recorder array of loaded
   packages */

static int function_package_record (fnpkg *pkg)
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

static int broken_package_error (fnpkg *pkg)
{
    gretl_errmsg_sprintf("'%s': package contains "
			 "duplicated function names",
			 pkg->name);
    return E_DATA;
}

#define fn_redef_msg(s) fprintf(stderr, "Redefining function '%s'\n", s)

/* When loading a private function the only real conflict would be
   with a function of the same name owned by the same package.
   Obviously this shouldn't happen but we'll whack it if it does.
*/

static int load_private_function (ufunc *fun)
{
    const char *targ = fun->name;
    int i, err = 0;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(targ, ufuns[i]->name)) {
	    if (fun->pkg == ufuns[i]->pkg) {
		err = broken_package_error(fun->pkg);
		break;
	    }
	}
    }

    if (!err) {
	err = add_allocated_ufunc(fun);
    }

    return err;
}

/* When loading a public, packaged function we want to avoid conflicts
   with any non-private functions of the same name.  In case we get a
   conflict with a public member of another package we'd best delete
   the entire conflicting package.
*/

static int load_public_function (ufunc *fun)
{
    const char *targ = fun->name;
    int i, done = 0;
    int err = 0;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(targ, ufuns[i]->name)) {
	    if (ufuns[i]->pkg == fun->pkg) {
		/* name duplication in package */
		err = broken_package_error(fun->pkg);
	    } else if (ufuns[i]->pkg == NULL) {
		/* conflicting unpackaged function */
		ufunc_free(ufuns[i]);
		ufuns[i] = fun;
		done = 1;
	    } else if (!function_is_private(ufuns[i])) {
		/* got a conflicting package */
		fprintf(stderr, "unloading conflicting package '%s'\n",
			ufuns[i]->pkg->name);
		fprintf(stderr, "(filename '%s')\n", ufuns[i]->pkg->fname); 
		function_package_unload_full(ufuns[i]->pkg);
		break;
	    } 
	}
    }

    if (!err && !done) {
	err = add_allocated_ufunc(fun);
    }

    return err;
}

/* A 'real load' is in contrast to just reading some info from a
   package, as is done in various GUI contexts.  We do the real load
   in response to some GUI commands, the "include" command, and also
   when re-opening a gretl session file that contains function
   packages.
*/

static int real_load_package (fnpkg *pkg)
{
    int i, err = 0;

#if PKG_DEBUG
    fprintf(stderr, "real_load_package:\n loading '%s'\n", pkg->fname);
#endif

    if (pkg->priv != NULL) {
	for (i=0; i<pkg->n_priv && !err; i++) {
	    err = load_private_function(pkg->priv[i]);
	}
    }

    if (!err && pkg->pub != NULL) {
	for (i=0; i<pkg->n_pub && !err; i++) {
	    err = load_public_function(pkg->pub[i]);
	}
    } 

    if (!err) {
	/* add to array of loaded packages */
	err = function_package_record(pkg);
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
    char vstr[8];
    int i;

    if (pkg->minver > 0) {
	get_version_string(vstr, pkg->minver);
    } else {
	*vstr = '\0';
    }

    pprintf(prn, "Package: %s\n", (*pkg->name)? pkg->name : "unknown");
    pprintf(prn, "Author: %s\n", (pkg->author)? pkg->author : "unknown");
    pprintf(prn, "Version: %s\n", (pkg->version)? pkg->version : "unknown");
    pprintf(prn, "Date: %s\n", (pkg->date)? pkg->date : "unknown");
    pprintf(prn, "Required gretl version: %s\n", (*vstr)? vstr : "unknown");
    pputs(prn, "Description: ");
    pputs(prn, (pkg->descrip)? pkg->descrip : "none");

    pputs(prn, "\n\n");

    if (pkg->pub != NULL) {
	if (pkg->n_pub == 1) {
	    if (strcmp(pkg->pub[0]->name, pkg->name)) {
		pputs(prn, "Public interface: ");
		pprintf(prn, "%s()\n", pkg->pub[0]->name);
	    }
	} else {
	    pputs(prn, "Public interfaces:\n");
	    for (i=0; i<pkg->n_pub; i++) {
		pprintf(prn, "  %s()\n", pkg->pub[i]->name);
	    }
	}
	pputc(prn, '\n');
    }    

    if (pkg->help != NULL) {
	pputs(prn, "Help text:\n");
	if (has_suffix(pkg->help, ".pdf")) {
	    const char *s = strrchr(pkg->help, ':');

	    if (s != NULL) {
		pprintf(prn, "See %s", s + 1);
	    } else {
		pputs(prn, pkg->help);
	    }
	} else {	
	    pputs(prn, pkg->help);
	}
	pprintf(prn, "\n\n");
    }

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
	    gretl_function_print_code(pkg->priv[i], prn);
	    pputc(prn, '\n');
	}
    }

    if (pkg->pub != NULL) {
	for (i=0; i<pkg->n_pub; i++) {
	    gretl_function_print_code(pkg->pub[i], prn);
	    pputc(prn, '\n');
	}
    }
}

/* allocate a fnpkg structure and read from XML file into it */

static fnpkg * 
real_read_package (xmlDocPtr doc, xmlNodePtr node, const char *fname, 
		   int *err)
{
    xmlNodePtr cur;
    fnpkg *pkg;
    char *tmp = NULL;
    int id;

#if PKG_DEBUG
    fprintf(stderr, "real_read_package: fname='%s'\n", fname);
#endif

    pkg = function_package_alloc(fname);
    if (pkg == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    gretl_xml_get_prop_as_string(node, "name", &tmp);

    if (tmp == NULL) {
	fprintf(stderr, "real_read_package: package has no name\n");
	*err = E_DATA;
	function_package_free(pkg);
	return NULL;
    }

    strncat(pkg->name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (gretl_xml_get_prop_as_bool(node, NEEDS_TS)) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_QM)) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (gretl_xml_get_prop_as_bool(node, NEEDS_PANEL)) {
	pkg->dreq = FN_NEEDS_PANEL;
    } else if (gretl_xml_get_prop_as_bool(node, NO_DATA_OK)) {
	pkg->dreq = FN_NODATA_OK;
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
	} else if (!xmlStrcmp(cur->name, (XUC) "help")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->help);
	} else if (!xmlStrcmp(cur->name, (XUC) "gui-help")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->gui_help);
	} else if (!xmlStrcmp(cur->name, (XUC) "sample-script")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->sample);
	} else if (!xmlStrcmp(cur->name, (XUC) "label")) {
	    gretl_xml_node_get_trimmed_string(cur, doc, &pkg->label);
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
    fprintf(stderr, "read_function_package: got '%s'\n", fname);
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

/** 
 * function_package_is_loaded:
 * @fname: full path to gfn file.
 *
 * Returns: 1 if the function package with filename @fname is
 * loaded in memory, otherwise 0.
 */

int function_package_is_loaded (const char *fname)
{
    return (get_loaded_pkg_by_filename(fname) != NULL);
}

/** 
 * load_function_package_by_filename:
 * @fname: full path to gfn file.
 *
 * Loads the function package located by @fname into
 * memory, if possible. Supports gretl's "include" command
 * for gfn files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int load_function_package_by_filename (const char *fname)
{
    int err = 0;

    if (function_package_is_loaded(fname)) {
	/* already loaded: no-op */
	fprintf(stderr, "load_function_package_from_file:\n"
		" '%s' is already loaded\n", fname);
    } else {
	fnpkg *pkg = read_package_file(fname, &err);

	if (!err) {
	    err = real_load_package(pkg);
	} 
    }

    if (err) {
	fprintf(stderr, "load function package: failed on %s\n", fname);
    }

    return err;
}

/** 
 * get_function_package_by_filename:
 * @fname: gfn filename.
 * @err: location to receive error code.
 *
 * If the package whose filename is @fname is already loaded,
 * returns the package pointer, otherwise attempts to load the
 * package from file. Similar to load_function_package_from_file()
 * but returns the package pointer rather than just a status
 * code.
 *
 * Returns: package-pointer on success, NULL on error.
 */

fnpkg *get_function_package_by_filename (const char *fname, int *err)
{
    fnpkg *pkg = NULL;
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(fname, pkgs[i]->fname)) {
	    pkg = pkgs[i];
	    break;
	}
    }

    if (pkg == NULL) {
	pkg = read_package_file(fname, err);
	if (!*err) {
	    *err = real_load_package(pkg);
	    if (*err) {
		pkg = NULL;
	    }
	}
    } 	

    return pkg;
}

fnpkg *get_function_package_by_name (const char *pkgname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
	if (!strcmp(pkgname, pkgs[i]->name)) {
	    return pkgs[i];
	}
    }

    return NULL;
}

/* Retrieve summary info or code listing for a function package,
   identified by its filename.  This is called (indirectly) from the
   GUI (see below for the actual callbacks).
*/

static int 
real_print_function_package_data (const char *fname, PRN *prn, int task)
{
    fnpkg *pkg;
    int free_pkg = 0;
    int err = 0;

#if PKG_DEBUG
    fprintf(stderr, "real_get_function_file_info: fname = '%s'\n", fname);
#endif

    pkg = get_loaded_pkg_by_filename(fname);

    if (pkg == NULL) {
	/* package is not loaded, read it now */
	pkg = read_package_file(fname, &err);
	free_pkg = 1;
    }

    if (!err) {
	if (task == FUNCS_INFO) {
	    print_package_info(pkg, prn);
	} else {
	    print_package_code(pkg, prn);
	}
	if (free_pkg) {
	    function_package_free_full(pkg);
	}
    }

#if PKG_DEBUG
    fprintf(stderr, "real_get_function_file_info: err = %d (free_pkg = %d)\n", 
	    err, free_pkg);
#endif

    return err;
}

/* callback used in the  GUI function package browser */

int print_function_package_info (const char *fname, PRN *prn)
{
    return real_print_function_package_data(fname, prn, FUNCS_INFO);
}

/* callback used in the  GUI function package browser */

int print_function_package_code (const char *fname, PRN *prn)
{
    return real_print_function_package_data(fname, prn, FUNCS_CODE);
}

/* Read the header from a function package file -- this is used when
   displaying the available packages in the GUI.  We write the
   description into *pdesc and a string representation of version
   number into *pver.
*/

int get_function_file_header (const char *fname, char **pdesc, 
			      char **pver)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr sub;
    int err = 0;

    xmlKeepBlanksDefault(0);

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
	return err;
    }

    node = node->xmlChildrenNode;
    while (node != NULL) {
	if (!xmlStrcmp(node->name, (XUC) "gretl-function-package")) {
	    sub = node->xmlChildrenNode;
	    while (sub != NULL) {
		if (!xmlStrcmp(sub->name, (XUC) "description")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, pdesc);
		} else if (!xmlStrcmp(sub->name, (XUC) "version")) {
		    gretl_xml_node_get_trimmed_string(sub, doc, pver);
		}
		if (*pdesc != NULL && *pver != NULL) {
		    break;
		}
		sub = sub->next;
	    }
	    if (*pdesc != NULL && *pver != NULL) {
		/* already got what we want */
		break;
	    }
	}
	node = node->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

    if (*pdesc == NULL) {
	*pdesc = gretl_strdup(_("No description available"));
    }

    if (*pver == NULL) {
	*pver = gretl_strdup("unknown");
    }

    if (*pdesc == NULL || *pver == NULL) {
	err = E_ALLOC;
    }

    return err;
}

/**
 * gretl_is_public_user_function:
 * @name: name to test.
 *
 * Returns: 1 if @name is the name of a user-defined
 * function that is not a private member of a function
 * package, otherwise 0.
 */

int gretl_is_public_user_function (const char *name)
{
    ufunc *fun = get_user_function_by_name(name);

    return (fun != NULL && !function_is_private(fun));
}

static int check_func_name (const char *fname, ufunc **pfun, PRN *prn)
{
    int i, err = 0;

#if FN_DEBUG
    fprintf(stderr, "check_func_name: '%s'\n", fname);
#endif

    if (!isalpha((unsigned char) *fname)) {
	gretl_errmsg_set(_("Function names must start with a letter"));
	err = 1;
    } else if (gretl_command_number(fname)) {
	gretl_errmsg_sprintf(_("'%s' is the name of a gretl command"),
			     fname);
	err = 1;
    } else if (function_lookup(fname)) {
	gretl_errmsg_sprintf(_("'%s' is the name of a built-in function"),
			     fname);
	err = 1;
    } else {
	for (i=0; i<n_ufuns; i++) {
	    if (!function_is_private(ufuns[i]) && !strcmp(fname, ufuns[i]->name)) {
#if FN_DEBUG
		fprintf(stderr, "'%s': found an existing function of this name\n", fname);
#endif
		if (pfun != NULL) {
		    clear_ufunc_data(ufuns[i]);
		    *pfun = ufuns[i];
		} else {
		    ufunc_unload(ufuns[i]);
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
	gretl_errmsg_sprintf("%s: function is in use", fname);
	err = 1;
    } else if (fun->pkg != NULL) {
	gretl_errmsg_sprintf("%s: function belongs to package", fname);
	err = 1;
    } else {
	ufunc_unload(fun);
	if (gretl_messages_on()) {
	    pprintf(prn, _("Deleted function '%s'\n"), fname);
	}
    } 

    return err;
}

/* next: apparatus for parsing function definitions */

static int comma_count (const char *s)
{
    int quoted = 0;
    int braced = 0;
    int nc = 0;

    while (*s) {
	if (!quoted) {
	    if (*s == '{') {
		braced++;
	    } else if (*s == '}') {
		braced--;
	    }
	}
	if (*s == '"') {
	    quoted = !quoted;
	} else if (!quoted && !braced && *s == ',') {
	    nc++;
	}
	s++;
    }

    return nc;
}

static int check_parm_min_max (fn_param *p, const char *name, int *nvals)
{
    if (!na(p->min) && !na(p->max)) {
	if (p->min > p->max) {
	    gretl_errmsg_sprintf("%s: min > max", name);
	    return E_DATA;
	} else if (!na(p->step) && p->step > p->max - p->min) {
	    gretl_errmsg_sprintf("%s: step > max - min", name);
	    return E_DATA;
	}	    
	*nvals = (int) p->max - (int) p->min + 1;
    }

    if (!na(p->deflt)) {
	if (!na(p->min) && p->deflt < p->min) {
	    gretl_errmsg_sprintf("%s: default value out of bounds", name);
	    return E_DATA;
	} 
	if (!na(p->max) && p->deflt > p->max) {
	    gretl_errmsg_sprintf("%s: default value out of bounds", name);
	    return E_DATA;
	}
    } 

    return 0;
}

static int read_min_max_deflt (char **ps, fn_param *param, const char *name,
			       int *nvals)
{
    char *p = *ps;
    double x, y, z, s;
    int err = 0;

    gretl_push_c_numeric_locale();

    if (param->type == GRETL_TYPE_BOOL) {
	if (sscanf(p, "[%lf]", &x) == 1) {
	    param->deflt = x;
	} else {
	    err = E_PARSE;
	}
    } else if (!strncmp(p, "[$xlist]", 8)) {
	param->deflt = INT_USE_XLIST;
    } else {
	if (sscanf(p, "[%lf:%lf:%lf:%lf]", &x, &y, &z, &s) == 4) {
	    param->min = x;
	    param->max = y;
	    param->deflt = z;
	    param->step = s;
	} else if (sscanf(p, "[%lf:%lf:%lf]", &x, &y, &z) == 3) {
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

	if (!err) {
	    err = check_parm_min_max(param, name, nvals);
	}
    }

    gretl_pop_c_numeric_locale();

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

/* get the descriptive string for a function parameter */

static int read_param_descrip (char **ps, fn_param *param)
{
    char *p = *ps + 1; /* skip opening quote */
    int len = 0;
    int err = E_PARSE;

    while (*p) {
	if (*p == '"') {
	    /* OK, found closing quote */
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

/* get the value labels for a function parameter */

static int read_param_labels (char **ps, fn_param *param, 
			      const char *name, int nvals)
{
    char *p = *ps + 1; /* skip opening '{' */
    int len = 0;
    int err = E_PARSE;

    while (*p) {
	if (*p == '}') {
	    /* OK, found closing brace */
	    err = 0;
	    break;
	}
	len++;
	p++;
    }

    if (!err && len > 0) {
	char *tmp;

	p = *ps + 1;
	tmp = gretl_strndup(p, len);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    *ps = p + len + 1;
	    param->labels = gretl_string_split_quoted(tmp, &param->nlabels,
						      &err);
	    free(tmp);
	    if (!err && param->nlabels != nvals) {
		gretl_errmsg_sprintf("%s: found %d values but %d value-labels",
				     name, nvals, param->nlabels);
		err = E_DATA;
	    }
	}
    }

    return err;
}

static void trash_param_info (char *name, fn_param *param)
{
    free(name);
    free(param->descrip);
    param->descrip = NULL;
    if (param->nlabels > 0) {
	free_strings_array(param->labels, param->nlabels);
	param->labels = NULL;
	param->nlabels = 0;
    }
}

static int parse_function_param (char *s, fn_param *param, int i)
{
    char tstr[22] = {0};
    char *name;
    int type, len;
    int nvals = 0;
    int err = 0;

#if FNPARSE_DEBUG
    fprintf(stderr, "parse_function_param: s = '%s'\n", s);
#endif

    while (isspace(*s)) s++;

    /* pick up the "const" flag if present */

    if (!strncmp(s, "const ", 6)) {
	param->flags |= ARG_CONST;
	s += 6;
	while (isspace(*s)) s++;
    }

    /* get parameter type -- required */

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
	    gretl_errmsg_set("Expected a type identifier");
	    err = E_PARSE;
	} else {
	    type = gretl_type_from_string(tstr);
	    if (type == 0) {
		gretl_errmsg_sprintf("Unrecognized data type '%s'", tstr);
		err = E_PARSE;
	    } else if (!ok_function_arg_type(type)) {
		gretl_errmsg_sprintf("%s: invalid parameter type", tstr);
		err = E_INVARG;
	    }		
	} 
    }

    if (err) {
	return err;
    }

    /* now get the parameter name -- required */
    
    while (isspace(*s)) s++;

    len = gretl_namechar_spn(s);

    if (len == 0) {
	gretl_errmsg_sprintf("parameter %d: name is missing", i + 1);
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

    /* now scan for various optional extras: first we may have
       a [min:max:default] specification, or a [null] spec
       to allow no argument for "pointer"-type arguments
    */

    if (*s == '[') { 
	if (gretl_scalar_type(type)) {
	    err = read_min_max_deflt(&s, param, name, &nvals);
	} else if (gretl_ref_type(type) || 
		   type == GRETL_TYPE_LIST ||
		   type == GRETL_TYPE_STRING) {
	    err = read_param_option(&s, param);
	} else {
	    err = E_PARSE;
	}
	s += strspn(s, " ");
    } 

    /* then we may have a double-quoted descriptive string
       for the parameter */

    if (!err && *s == '"') {
	err = read_param_descrip(&s, param);
	s += strspn(s, " ");
    } 

    /* and finally we may have a set of value-labels enclosed
       in braces, but this is valid only if we have been able
       to determine a definite number of admissible values
       for the parameter
    */

    if (!err && *s == '{') {
	if (nvals == 0) {
	    err = E_PARSE;
	} else {
	    err = read_param_labels(&s, param, name, nvals);
	}
	s += strspn(s, " ");
    }    

    if (!err && *s != '\0') {
	/* got trailing unparseable stuff */
	err = E_PARSE;
    }

    if (!err) {
	param->name = name;
    } else {
	trash_param_info(name, param);
    }

#if FNPARSE_DEBUG
    if (!err) {
	fprintf(stderr, " param[%d] = '%s', ptype = %d\n", 
		i, name, type);
	fprintf(stderr, "  min=%g, max=%g, deflt=%g\n", 
		param->min, param->max, param->deflt);
	if (param->descrip != NULL) {
	    fprintf(stderr, "  comment = '%s'\n", param->descrip);
	}
	if (param->nlabels > 0) {
	    fprintf(stderr, "  value labels: %d\n", param->nlabels);
	}	
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

/* Here we're parsing what follows "function ".  The expectation
   is that we get <return-type> <funcname> (<args>), but we also
   support the old-style <funcname> (<args>).
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

#if FNPARSE_DEBUG > 1
    fprintf(stderr, "parse_fn_definition:\n '%s'\n", str);
#endif

    len = strlen(str);
    if (str[len-1] != ')') {
	/* somehow we didn't get a properly terminated 
	   param list */
	return E_PARSE;
    }

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
	    if (!ok_function_return_type(t)) {
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
	/* FIXME allow '(' then newline? */
	return E_PARSE;
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

#if FNPARSE_DEBUG
    fprintf(stderr, "function %s: looking for %d parameters\n",
	    fname, np);
#endif

    params = allocate_params(np);
    if (params == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	int quoted = 0;
	int braced = 0;

	p = s;
	while (*p) {
	    if (!quoted) {
		if (*p == '{') {
		    braced++;
		} else if (*p == '}') {
		    braced--;
		}
	    }
	    if (*p == '"') {
		quoted = !quoted;
	    } else if (!quoted && !braced && *p == ',') {
		*p = '\0';
	    }
	    p++;
	}
	p = s;
	if (braced != 0) {
	    err = E_PARSE;
	}
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

/**
 * gretl_start_compiling_function:
 * @line: command line.
 * @prn: printing struct for feedback.
 *
 * Responds to a command of the form "function ...".  In most
 * cases, embarks on compilation of a function, but this
 * also handles the construction "function foo delete".
 *
 * Returns: 0 on success, non-zero on error.
 */

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
    if (err) {
	pprintf(prn, "> %s\n", line);
    }

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
	current_fdef = fun;
	set_compiling_on();
    } else {
	current_fdef = NULL;
    }
    
    return err;
}

/* not reading from XML, but from script or other
   direct user input */

static int parse_function_return (ufunc *fun, const char *line)
{
    const char *s = line + 6; /* skip "return" */
    char s1[16], s2[VNAMELEN];
    int type, err = 0;

#if FNPARSE_DEBUG
    fprintf(stderr, "parse_function_return: line = '%s'\n", line);
#endif

    *s1 = *s2 = '\0';
    sscanf(s, "%15s %15s", s1, s2);
    type = return_type_from_string(s1);

    if (type > 0 && fun->rettype != GRETL_TYPE_NONE) {
	gretl_errmsg_sprintf("%s: return type is already defined",
			     fun->name);
	err = E_PARSE;
    } else if (fun->rettype != GRETL_TYPE_NONE) {
	/* new-style: the return type was pre-defined, and
	   we treat a "return " line as a normal line 
	*/
	err = strings_array_add(&fun->lines, &fun->n_lines, line);
    } else {
	/* old-style: we didn't get a return type at the start of the
	   definition, and it should be specified inline here
	*/
	if (type == 0) {
	    gretl_errmsg_sprintf("%s: missing a valid return type\n", 
				 fun->name);
	    err = E_TYPES;
	} else {
	    err = check_varname(s2);
	    if (!err) {
		err = ufunc_add_return_statement(fun, s2);
	    }
	    if (!err) {
		fun->rettype = type;
	    }
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

/* Rather minimal check for syntactic validity of "compiled" function.
   FIXME: it would be good to check here for messed up block structure
   too (e.g. "system" without "end system").
*/

static int check_function_structure (ufunc *fun)
{
    char line[MAXLINE];
    CMD cmd;
    int ifdepth = 0;
    int i, err = 0;

    gretl_cmd_init(&cmd);

#if 0
    fprintf(stderr, "checking function '%s'\n", fun->name);
#endif


    for (i=0; i<fun->n_lines && !err; i++) {
#if 0
	fprintf(stderr, "line[%d] = '%s'\n", i, fun->lines[i]);
#endif
	/* avoid losing comment lines */
	strcpy(line, fun->lines[i]);
	get_command_index(line, &cmd);
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
	fprintf(stderr, "After reading, ifdepth = %d\n", ifdepth);
	err = E_PARSE;
    }

    gretl_cmd_free(&cmd);

    return err;
}

/* Note: the input @fun may be NULL; in that case the assumption is
   that we're appending to a function-in-progress that is
   represented by the internal variable 'current_fdef'.  If
   the @fun argument is non-NULL that means we're editing an
   existing function.
*/

static int real_function_append_line (const char *line, ufunc *fun)
{
    int editing = 1;
    int err = 0;

    if (fun == NULL) {
	fun = current_fdef;
	editing = 0;
    } 

#if FNPARSE_DEBUG
    fprintf(stderr, "real_function_append_line: '%s'\n", line);
#endif

    if (fun == NULL) {
#if FN_DEBUG
	fprintf(stderr, "real_function_append_line: fun is NULL\n");
#endif
	return 1;
    }

    if (string_is_blank(line)) {
	err = strings_array_add(&fun->lines, &fun->n_lines, "");
    } else if (end_of_function(line) && !ignore_line(fun)) {
	if (fun->n_lines == 0) {
	    gretl_errmsg_sprintf("%s: empty function", fun->name);
	    err = 1;
	}
	set_compiling_off();
    } else if (!strncmp(line, "quit", 4)) {
	/* abort compilation */
	if (!editing) {
	    ufunc_unload(fun);
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

    if (!err && !compiling && !editing) {
	/* finished composing function */
	err = check_function_structure(fun);
    }

    if (err && !editing) {
	ufunc_unload(fun);
    }	

    return err;
}

/**
 * gretl_function_append_line:
 * @line: line of code to append.
 *
 * Continuation of definition of user-function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_function_append_line (const char *line)
{
    return real_function_append_line(line, NULL);
}

/* Given a line "function ..." get the function name, with
   some error checking.  The @s we are given here is at an 
   offset of 9 bytes into the line, skipping "function ".
*/

static int extract_funcname (char *name, const char *s)
{
    char word[FN_NAMELEN];
    int n, type, err = 0;

    s += strspn(s, " ");
    n = strcspn(s, " (");

    if (n == 0 || n > FN_NAMELEN - 1) {
	return E_DATA;
    }

    *word = '\0';
    strncat(word, s, n);

    if (!strcmp(word, "void")) {
	type = GRETL_TYPE_VOID;
    } else {
	type = gretl_type_from_string(word);
    }

    if (type == 0) {
	/* old-style */
	strcpy(name, word);
    } else if (!ok_function_return_type(type)) {
	err = E_DATA;
    } else {
	s += n;
	s += strspn(s, " ");
	n = strcspn(s, " (");
	if (n == 0 || n > FN_NAMELEN - 1) {
	    err = E_DATA;
	} else {
	    *name = '\0';
	    strncat(name, s, n);
	}
    }

    return err;
}

/* allow for line-continuation in a function definition */

static int fndef_append_next (char *s, FILE *fp)
{
    char line[MAXLINE];
    int n = strlen(s);

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
    } else {
	strcat(s, line);
	return 0;
    }
}

static int fndef_continues (const char *s)
{
    int n = strlen(s);

    return (s[n-1] == ',' || s[n-1] == '\\');
}

/* Called from the GUI package-editing apparatus.  The user is editing
   a function named @funname, and the content of the script-editing
   window has been dumped to file (given by @path).  We open the file
   and parse the new function definition, but hold off replacing the
   original definition until we know there are no (obvious)
   compilation errors.  We try to detect if the name of the function
   itself has been changed in the GUI editor window: that is not
   acceptable.

   The @pkg argument will be non-NULL if an existing package is being
   edited, but will be NULL when a new package is under construction
   and has not yet been saved,
*/

int update_function_from_script (const char *funname, const char *path,
				 fnpkg *pkg)
{
    char line[MAXLINE];
    ufunc *fun, *orig;
    char *s;
    FILE *fp;
    int gotfn = 0;
    int err = 0;

    if (pkg == NULL) {
	orig = get_user_function_by_name(funname);
    } else {
	orig = get_function_from_package(funname, pkg);
    }

    if (orig == NULL) {
	gretl_errmsg_sprintf("Couldn't find function '%s'", funname);
	return E_DATA;
    }

    fp = gretl_fopen(path, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    /* temporary function struct */
    fun = ufunc_new();
    if (fun == NULL) {
	fclose(fp);
	return E_ALLOC;
    }

    fprintf(stderr, "Updating function '%s' from %s\n", funname, path);

    while (fgets(line, sizeof line, fp) && !err) {
	s = line;
	while (*s == ' ') s++;
	tailstrip(s);
	if (!strncmp(s, "function ", 9)) {
	    if (gotfn || extract_funcname(fun->name, s + 9)) {
		err = E_DATA;
	    } else if (strcmp(fun->name, orig->name)) {
		err = E_DATA;
		gretl_errmsg_set(_("You can't change the name of a function here"));
		fprintf(stderr, "orig->name='%s', fun->name='%s'\n",
		       orig->name, fun->name);
	    } else {
		gotfn = 1;
		while (!err && fndef_continues(s)) {
		    err = fndef_append_next(s, fp);
		}
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
	err = check_function_structure(fun);
    }

    if (err) {
	/* didn't work: trash the attempt */
	free_strings_array(fun->lines, fun->n_lines);
	free_params_array(fun->params, fun->n_params);
    } else {	
	/* really replace the original function content */
	free_strings_array(orig->lines, orig->n_lines);
	orig->n_lines = fun->n_lines;
	orig->lines = fun->lines;
	fun->lines = NULL;

	free_params_array(orig->params, orig->n_params);
	orig->n_params = fun->n_params;
	orig->params = fun->params;
	fun->params = NULL;

	orig->rettype = fun->rettype;
    } 

    free(fun);

    return err;
}

/* next block: handling function arguments */

/* Given a named list of variables supplied as an argument to a
   function, copy the list under the name assigned by the function,
   and make the variables referenced in that list accessible to the
   function.
*/

static int localize_list (fncall *call, const char *oldname, 
			  fn_param *fp, DATASET *dset)
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
		STACK_LEVEL(dset, vi) = level;
	    }
	}
    }

#if UDEBUG
    fprintf(stderr, "localize_list: returning %d\n", err);
#endif

    return err;
}

/* coerce scalar value to boolean, as per parameter spec */

static void boolify_local_var (const char *vname)
{
    double x = gretl_scalar_get_value(vname);

    if (x != 0.0 && !na(x)) {
	gretl_scalar_set_value(vname, 1.0);
    }
}

static void maybe_set_arg_const (struct fnarg *arg, fn_param *fp)
{
    if (fp->flags & ARG_CONST) {
	/* param is marked CONST directly */
	arg->name = fp->name;
	arg->flags |= ARG_CONST;
    } else if (object_is_const(arg->upname)) {
	/* param is CONST by inheritance */
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
    int ulev = user_matrix_get_level(u);

    if (ulev == fn_executing + 1) {
	gretl_errmsg_set(_("Duplicated pointer argument: not allowed"));
	return E_DATA;
    }

    arg->upname = gretl_strdup(user_matrix_get_name(u));
    if (arg->upname == NULL) {
	return E_ALLOC;
    }

    user_matrix_adjust_level(u, 1);
    user_matrix_set_name(u, fp->name);

    maybe_set_arg_const(arg, fp);

    return 0;
}

static int localize_bundle_ref (fnargs *args, int i,
				fn_param *fp)
{
    struct fnarg *arg = args->arg[i];
    const char *alt, *bname = arg->val.str;
    int j, err = 0;

    for (j=0; j<i; j++) {
	alt = args->arg[j]->upname;
	if (alt != NULL && !strcmp(alt, bname)) {
	    gretl_errmsg_set(_("Duplicated pointer argument: not allowed"));
	    return E_DATA;
	}
    }

    arg->upname = gretl_strdup(bname);
    if (arg->upname == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_bundle_localize(bname, fp->name);
    }

    return err;
}

static int localize_scalar_ref (struct fnarg *arg, fn_param *fp)
{
    int i = arg->val.idnum;
    const char *s = gretl_scalar_get_name(i);

    if (s == NULL) {
	return E_DATA;
    } 

    if (gretl_scalar_get_level(i) == fn_executing + 1) {
	gretl_errmsg_set(_("Duplicated pointer argument: not allowed"));
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
				fn_param *fp, DATASET *dset)
{
    int v = arg->val.idnum;

    if (STACK_LEVEL(dset, v) == fn_executing + 1) {
	gretl_errmsg_set(_("Duplicated pointer argument: not allowed"));
	return E_DATA;
    }

    arg->upname = gretl_strdup(dset->varname[v]);
    if (arg->upname == NULL) {
	return E_ALLOC;
    } 

    STACK_LEVEL(dset, v) += 1;
    strcpy(dset->varname[v], fp->name);

    if (!in_gretl_list(call->ptrvars, v)) {
	gretl_list_append_term(&call->ptrvars, v);
    }

    maybe_set_arg_const(arg, fp);

    return 0;
}

/* Scalar function arguments only: if the arg is not supplied, use the
   default that is found in the function specification, if any.
*/

static int add_scalar_arg_default (fn_param *param)
{
    double x;

    if (na(param->deflt)) {
	/* should be impossible here, but... */
	return E_DATA;
    }

    if (param->type == GRETL_TYPE_BOOL || 
	param->type == GRETL_TYPE_INT ||
	param->type == GRETL_TYPE_OBS) {
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

static int allocate_function_args (fncall *call, DATASET *dset)
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
					       dset);
	    } else {
		err = dataset_add_series_as(arg->val.px, fp->name, 
					    dset);
	    }	    
	} else if (fp->type == GRETL_TYPE_MATRIX) {
	    if (fp->flags & ARG_CONST) {
		err = localize_const_matrix(arg, fp);
	    } else {
		err = copy_matrix_as(arg->val.m, fp->name);
	    }
	} else if (fp->type == GRETL_TYPE_BUNDLE) {
	    err = copy_bundle_arg_as(arg->val.b, fp->name);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    if (arg->type == GRETL_TYPE_NONE) {
		err = create_named_null_list(fp->name);
	    } else if (arg->type == GRETL_TYPE_USERIES) {
		err = create_named_singleton_list(arg->val.idnum, fp->name);
	    } else {
		err = localize_list(call, arg->val.str, fp, dset);
	    }
	} else if (fp->type == GRETL_TYPE_STRING) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = add_string_as(arg->val.str, fp->name);
	    }
	} else if (fp->type == GRETL_TYPE_SCALAR_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_scalar_ref(arg, fp);
	    }
	} else if (fp->type == GRETL_TYPE_SERIES_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_series_ref(call, arg, fp, dset);
	    }
	} else if (fp->type == GRETL_TYPE_MATRIX_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_matrix_ref(arg, fp);
	    }
	} else if (fp->type == GRETL_TYPE_BUNDLE_REF) {
	    if (arg->type != GRETL_TYPE_NONE) {
		err = localize_bundle_ref(args, i, fp);
	    }
	} 

	if (!err) {
	    if (arg->type == GRETL_TYPE_USERIES && fp->type != GRETL_TYPE_LIST) {
		arg->upname = gretl_strdup(dset->varname[arg->val.idnum]);
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
	set_listargs_from_call(call, dset);
    }

#if UDEBUG
    fprintf(stderr, "allocate_function args: returning %d\n", err);
#endif

    return err;
}

/**
 * check_function_needs:
 * @dset: pointer to dataset info.
 * @dreq: function data requirements flag.
 * @minver: function minimum program version requirement.
 *
 * Checks whether the requirements in @dreq and @minver
 * are jointly satisfied by the current dataset and gretl
 * version.
 *
 * Returns: 0 if all is OK, 1 otherwise.
 */

int check_function_needs (const DATASET *dset, FuncDataReq dreq,
			  int minver)
{
    static int thisver = 0;

    if (thisver == 0) {
	thisver = version_number_from_string(GRETL_VERSION);
    }

    if (minver > thisver) {
	char vstr[8];

	get_version_string(vstr, minver);
	gretl_errmsg_sprintf("This function needs gretl version %s", vstr);
	return 1;
    }

    if ((dset == NULL || dset->v == 0) && dreq != FN_NODATA_OK) {
	gretl_errmsg_set("This function needs a dataset in place");
	return 1;
    }

    if (dreq == FN_NEEDS_TS && !dataset_is_time_series(dset)) {
	gretl_errmsg_set("This function needs time-series data");
	return 1;
    }

    if (dreq == FN_NEEDS_PANEL && !dataset_is_panel(dset)) {
	gretl_errmsg_set("This function needs panel data");
	return 1;
    }

    if (dreq == FN_NEEDS_QM && 
	(!dataset_is_time_series(dset) || 
	 (dset->pd != 4 || dset->pd != 12))) {
	gretl_errmsg_set("This function needs quarterly or monthly data");
	return 1;
    } 

    return 0;
}

static int maybe_check_function_needs (const DATASET *dset,
				       const ufunc *fun)
{
    if (fun->pkg == NULL) {
	return 0;
    } else {
	return check_function_needs(dset, fun->pkg->dreq, 
				    fun->pkg->minver);
    }
}

/* next block: handling function return values */

static int handle_scalar_return (const char *vname, void *ptr)
{
    int err = 0;

    if (gretl_is_scalar(vname)) {
	*(double *) ptr = gretl_scalar_get_value(vname);
    } else {
	*(double *) ptr = NADBL;
	err = E_UNKVAR; /* FIXME */
    }

    return err;
}

static int handle_series_return (const char *vname, void *ptr,
				 DATASET *dset, int copy)
{
    int v = series_index(dset, vname);
    double *x = NULL;
    int err = 0;

    if (!copy && v == 0) {
	copy = 1;
    }
    if (v >= 0 && v < dset->v) {
	if (copy) {
	    x = copyvec(dset->Z[v], dset->n);
	    if (x == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    x = dset->Z[v];
	    dset->Z[v] = NULL;
	}
    } else {
	err = E_UNKVAR;
    }

    *(double **) ptr = x;

    return err;
}

static int handle_matrix_return (const char *name, void *ptr,
				 int copy)
{
    gretl_matrix *ret = NULL;
    int err = 0;

    if (copy) {
	gretl_matrix *m = get_matrix_by_name(name);

	if (m != NULL) {
	    ret = gretl_matrix_copy(m);
	    if (ret == NULL) {
		err = E_ALLOC;
	    }
	} 
    } else {
	ret = steal_matrix_by_name(name);
    }

    if (ret == NULL && !err) {
	err = E_UNKVAR;
    }

    *(gretl_matrix **) ptr = ret;

    return err;
}

static int handle_bundle_return (fncall *call, void *ptr, int copy)
{
    const char *name = call->retname;
    gretl_bundle *ret = NULL;
    int err = 0;

    if (copy) {
	gretl_bundle *b = get_gretl_bundle_by_name(name);

	if (b != NULL) {
	    ret = gretl_bundle_copy(b, &err);
	} 
    } else {
	ret = gretl_bundle_pull_from_stack(name, &err);
    }

    if (ret != NULL && call->fun->pkg != NULL &&
	gretl_function_depth() == 1) {
	gretl_bundle_set_creator(ret, call->fun->pkg->name);
    }

    *(gretl_bundle **) ptr = ret;

    return err;
}

enum {
    LIST_DIRECT_RETURN,
    LIST_WAS_ARG,
    LIST_WAS_TEMP
};

/* Deal with a list that exists at the level of a user-defined
   function whose execution is now terminating.  Note that this list
   may be the direct return value of the function, or it may have been
   given to the function as an argument, or it may have been constructed
   on the fly.  This is flagged by @status.  
*/

static int unlocalize_list (const char *lname, int status,
			    DATASET *dset)
{
    int *list = get_list_by_name(lname);
    int d = gretl_function_depth();
    int upd = d - 1;
    int i, vi;

#if UDEBUG
    fprintf(stderr, "unlocalize_list: '%s', function depth = %d\n", lname, d);
    printlist(list, lname);
    fprintf(stderr, " dset = %p, dset->v = %d\n", (void *) dset, dset->v);
#endif	    

    if (list == NULL) {
	return E_DATA;
    }

    /* 
       If the list we're looking at was given as a function argument
       we simply shunt all its members to the prior STACK_LEVEL.  But
       if the list is the direct return value from the function, we
       need to overwrite any variables at caller level that have been
       redefined within the function.  If any series have been
       redefined in that way we also need to adjust the list itself,
       replacing the ID numbers of local variables with the
       caller-level IDs.

       On the other hand, if the list was a temporary construction
       (a list parameter was required, but a single series was given as
       argument), we destroy the list.
    */

    if (status == LIST_DIRECT_RETURN) {
	int t, j, overwrite;
	const char *vname;
	
	for (i=1; i<=list[0]; i++) {
	    overwrite = 0;
	    vi = list[i];
	    vname = dset->varname[vi];
	    unset_var_listarg(dset, vi);
	    if (vi > 0 && vi < dset->v && STACK_LEVEL(dset, vi) == d) {
		for (j=1; j<dset->v; j++) { 
		    if (STACK_LEVEL(dset, j) == upd && 
			!strcmp(dset->varname[j], vname)) {
			overwrite = 1;
			break;
		    }
		}
		if (overwrite) {
		    /* replace data */
		    for (t=dset->t1; t<=dset->t2; t++) {
			dset->Z[j][t] = dset->Z[vi][t];
		    }
		    /* replace variable info */
		    strcpy(VARLABEL(dset, j), VARLABEL(dset, vi));
		    strcpy(DISPLAYNAME(dset, j), DISPLAYNAME(dset, vi));
		    /* replace ID number in list */
		    list[i] = j;
		} else {
		    STACK_LEVEL(dset, vi) = upd;
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
    } else if (status == LIST_WAS_ARG) {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    unset_var_listarg(dset, vi);
	    STACK_LEVEL(dset, vi) = upd;
	}
    } else if (status == LIST_WAS_TEMP) {
	delete_list_by_name(lname);
    }

    return 0;
}

static int handle_string_return (const char *sname, void *ptr)
{
    const char *s = get_string_by_name(sname);
    char *ret = NULL;
    int err = 0;

    if (s == NULL) {
	err = E_DATA;
    } else {
	ret = gretl_strdup(s);
	if (ret == NULL) {
	    err = E_ALLOC;
	}
    }

    *(char **) ptr = ret;

    return err;
}

static void 
maybe_set_return_description (fncall *call, int rtype, 
			      DATASET *dset, 
			      char **descrip)
{
    if (rtype == GRETL_TYPE_SERIES) {
	int v = series_index(dset, call->retname);

	if (v < dset->v) {
	    *descrip = gretl_strdup(VARLABEL(dset, v));
	}
    }
}

#define types_match(pt,rt) ((pt==GRETL_TYPE_SERIES_REF && rt==GRETL_TYPE_SERIES) || \
                            (pt==GRETL_TYPE_MATRIX_REF && rt==GRETL_TYPE_MATRIX))

static int is_pointer_arg (fncall *call, fnargs *args, int rtype)
{
    ufunc *u = call->fun;
    int i;

    if (call->retname == NULL) {
	return 0;
    }

    for (i=0; i<args->argc; i++) {
	if (types_match(u->params[i].type, rtype)) {
	    if (!strcmp(u->params[i].name, call->retname)) {
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
function_assign_returns (fncall *call, fnargs *args, int rtype, 
			 DATASET *dset, void *ret, char **descrip, 
			 PRN *prn, int *perr)
{
    ufunc *u = call->fun;
    struct fnarg *arg;
    fn_param *fp;
    int copy, i, err = 0;

#if UDEBUG
    fprintf(stderr, "function_assign_returns: rtype = %d, call->retname = %s\n", 
	    rtype, call->retname);
#endif

    /* check for missing return value */
    if (*perr == 0 && !null_return(rtype) && call->retname == NULL) {
	gretl_errmsg_sprintf("Function %s did not provide the specified return value",
			     u->name);
	*perr = err = E_UNKVAR;
    }

    if (*perr == 0) {
	/* first we work on the value directly returned by the
	   function (but only if there's no error) 
	*/
	if (needs_datainfo(rtype) && dset == NULL) {
	    /* "can't happen" */
	    err = E_DATA;
	} else if (rtype == GRETL_TYPE_DOUBLE) {
	    err = handle_scalar_return(call->retname, ret);
	} else if (rtype == GRETL_TYPE_SERIES) {
	    copy = is_pointer_arg(call, args, rtype);
	    err = handle_series_return(call->retname, ret, dset, copy);
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    copy = is_pointer_arg(call, args, rtype);
	    err = handle_matrix_return(call->retname, ret, copy);
	} else if (rtype == GRETL_TYPE_LIST) {
	    /* note: in this case the job is finished in
	       stop_fncall(); here we just adjust the info on the
	       listed variables so they don't get deleted
	    */
	    err = unlocalize_list(call->retname, LIST_DIRECT_RETURN, dset);
	} else if (rtype == GRETL_TYPE_BUNDLE) {
	    copy = is_pointer_arg(call, args, rtype);
	    err = handle_bundle_return(call, ret, copy);
	} else if (rtype == GRETL_TYPE_STRING) {
	    err = handle_string_return(call->retname, ret);
	} 

	if (err == E_UNKVAR) {
	    pprintf(prn, "Function %s did not provide the specified return value\n",
		    u->name);
	}

	*perr = err;

	if (!err && dset != NULL && descrip != NULL) {
	    maybe_set_return_description(call, rtype, dset, descrip);
	}
    }

    /* "indirect return" values and other pointerized args: 
       these should be handled even if the function bombed.
    */

    for (i=0; i<args->argc; i++) {
	arg = args->arg[i];
	fp = &u->params[i];
	if (needs_datainfo(fp->type) && dset == NULL) {
	    err = E_DATA;
	} else if (gretl_ref_type(fp->type)) {
	    if (arg->type == GRETL_TYPE_SERIES_REF) {
		int v = arg->val.idnum;

		STACK_LEVEL(dset, v) -= 1;
		strcpy(dset->varname[v], arg->upname);
	    } else if (arg->type == GRETL_TYPE_SCALAR_REF) {
		gretl_scalar_restore_name(arg->val.idnum, arg->upname);
	    } else if (arg->type == GRETL_TYPE_MATRIX_REF) {
		user_matrix *u = arg->val.um;

		user_matrix_adjust_level(u, -1);
		user_matrix_set_name(u, arg->upname);
	    } else if (arg->type == GRETL_TYPE_BUNDLE_REF) {
		gretl_bundle_unlocalize(fp->name, arg->upname);
	    }
	} else if (fp->type == GRETL_TYPE_MATRIX && arg->upname != NULL) {
	    user_matrix *u = get_user_matrix_by_data(arg->val.m);

	    user_matrix_adjust_level(u, -1);
	    user_matrix_set_name(u, arg->upname);
	} else if (fp->type == GRETL_TYPE_LIST) {
	    if (arg->type == GRETL_TYPE_LIST) {
		unlocalize_list(fp->name, LIST_WAS_ARG, dset);
	    } else {
		unlocalize_list(fp->name, LIST_WAS_TEMP, dset);
	    }
	} 
    }

#if UDEBUG
    fprintf(stderr, "function_assign_returns: returning %d\n", err);
#endif

    return err;
}

/* make a record of the sample information at the time a function is
   called */

static void record_obs_info (obsinfo *o, DATASET *dset)
{
    o->changed = 0;

    if (dset != NULL) {
	o->structure = dset->structure;
	o->pd = dset->pd;
	o->t1 = dset->t1;
	o->t2 = dset->t2;
	strcpy(o->stobs, dset->stobs);
    }
}

/* on function exit, restore the sample information that was in force
   on entry */

static int restore_obs_info (obsinfo *o, DATASET *dset)
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

    return set_obs(line, dset, opt);
}

/* do the basic housekeeping that is required when a function exits:
   destroy local variables, restore previous sample info, etc.
*/

static int stop_fncall (fncall *call, int rtype, void *ret,
			DATASET *dset, PRN *prn, int orig_v)
{
    int i, d = gretl_function_depth();
    int delv, anyerr = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "stop_fncall: terminating call to "
	    "function '%s' at depth %d, dset->v = %d\n", 
	    call->fun->name, d, (dset != NULL)? dset->v : 0);
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
       delete any local vars that have been "promoted" to caller
       level via their inclusion in a returned list
    */

    if (dset != NULL) {
	for (i=orig_v, delv=0; i<dset->v; i++) {
	    if (STACK_LEVEL(dset, i) == d) {
		delv++;
	    }
	}
	if (delv > 0) {
	    if (delv == dset->v - orig_v) {
		/* deleting all added variables */
		anyerr = dataset_drop_last_variables(delv, dset);
		if (anyerr && !err) {
		    err = anyerr;
		}
	    } else {
		for (i=dset->v-1; i>=orig_v; i--) {
		    if (STACK_LEVEL(dset, i) == d) {
			anyerr = dataset_drop_variable(i, dset);
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
	int *lret = gretl_list_copy(get_list_by_name(call->retname));

	if (lret != NULL) {
	    *(int **) ret = lret;
	} else {
	    err = E_ALLOC;
	}
    }    

    anyerr = destroy_saved_lists_at_level(d);

#if FN_DEBUG
    fprintf(stderr, "destroy_saved_lists_at_level(%d): err = %d\n", d, anyerr);
#endif


    if (anyerr && !err) {
	err = anyerr;
    }    

    /* if the function defined a Kalman filter, clean that up */
    delete_kalman(NULL);

    /* if any bundles were defined but not returned, clean up */
    destroy_saved_bundles_at_level(d);

    pop_program_state();

    if (dset != NULL && call->obs.changed) {
	restore_obs_info(&call->obs, dset);
    }

    set_executing_off(call, dset, prn);

    return err;
}

static void set_pkgdir (fnpkg *pkg)
{
    const char *p = strrchr(pkg->fname, SLASH);

    if (p != NULL) {
	char *pkgdir = gretl_strndup(pkg->fname, p - pkg->fname);

	gretl_insert_builtin_string("pkgdir", pkgdir);
	free(pkgdir);
    }
}

static int start_fncall (fncall *call, DATASET *dset, PRN *prn)
{
    fn_executing++;
    push_program_state();

    callstack = g_list_append(callstack, call);
#if EXEC_DEBUG
    fprintf(stderr, "start_fncall: added call to %s, depth now %d\n", 
	    call->fun->name, g_list_length(callstack));
#endif

    record_obs_info(&call->obs, dset);

    if (gretl_debugging_on() || call->fun->debug) {
	set_gretl_echo(1);
	set_gretl_messages(1);
	pprintf(prn, "*** executing function %s\n", call->fun->name);
    } else {
	set_gretl_echo(0);
	set_gretl_messages(0);
    }

    if (fn_executing == 1 && call->fun->pkg != NULL) {
	set_pkgdir(call->fun->pkg);
    } 

    return 0;
}

static void func_exec_callback (ExecState *s, void *ptr,
				GretlObjType type)
{
    int ci = s->cmd->ci;

    if (ci == GNUPLOT) {
	/* we permit "reach-back" into the GUI here */
	EXEC_CALLBACK gc = get_gui_callback();

	if (gc != NULL) {
	    gc(s, NULL, 0);
	}
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

static int check_function_args (ufunc *u, fnargs *args, PRN *prn)
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
	} else if (fp->type == GRETL_TYPE_LIST && arg->type == GRETL_TYPE_USERIES) {
	    ; /* OK */
	} else if (fp->type != arg->type) {
	    /* FIXME case where list is wanted and the arg is a series? */
	    pprintf(prn, _("%s: argument %d is of the wrong type (is %s, should be %s)\n"), 
		    u->name, i + 1, gretl_arg_type_name(arg->type), 
		    gretl_arg_type_name(fp->type));
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
			       DATASET *dset,
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
				 dset);
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
	    err = gretl_cmd_exec(state, dset);
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

static void set_function_error_message (int err, ufunc *u, 
					ExecState *state,
					const char *line)
{
    if (err == E_FUNCERR) {
	/* let the function writer set the message */
	gretl_errmsg_sprintf(_("Error message from %s():\n %s"), u->name,
			     state->cmd->param);
    } else {
	/* we'll handle this ourselves */
	const char *msg = gretl_errmsg_get();

	if (*msg == '\0') {
	    msg = errmsg_get_with_default(err);
	    gretl_errmsg_set(msg);
	} 

	if (*line != '\0') {
	    gretl_errmsg_sprintf("*** error in function %s\n> %s", u->name, line);
	} else {
	    gretl_errmsg_sprintf("*** error in function %s\n", u->name);
	}

	if (gretl_function_depth() > 1) {
	    GList *tmp = g_list_last(callstack);

	    if (tmp != NULL) {
		tmp = g_list_previous(tmp);
		if (tmp != NULL) {
		    fncall *call = tmp->data;

		    gretl_errmsg_sprintf(" called by function %s", call->fun->name);
		}
	    }
	}
    }
}

#define void_function(f) (f->rettype == 0 || f->rettype == GRETL_TYPE_VOID)

static int handle_return_statement (fncall *call,
				    ExecState *state,
				    DATASET *dset)
{
    const char *s = state->line + 6; /* skip "return" */
    ufunc *fun = call->fun;
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
	    call->retname = gretl_strndup(s, len);
	} else {
	    const char *typestr = gretl_arg_type_name(fun->rettype);
	    char formula[MAXLINE];
	    
	    sprintf(formula, "%s $retval=%s", typestr, s);
	    err = generate(formula, dset, OPT_P, NULL);
	    if (err) {
		set_function_error_message(err, fun, state, s);
	    } else {
		call->retname = gretl_strdup("$retval");
	    }
	}
    }

    if (!err && call->retname == NULL) {
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

static int do_debugging (ExecState *s)
{
    return s->cmd->ci > 0 && !s->cmd->context &&
	!gretl_compiling_loop();
}

#define CDEBUG 1

/* Construct bundle of arguments to send to C function;
   load plugin; send arguments; and try to retrieve
   return value, if any.
*/

static int handle_plugin_call (ufunc *u, fnargs *args,
			       DATASET *dset,
			       void *retval,
			       PRN *prn)
{
    int (*cfunc) (gretl_bundle *, PRN *);
    void *handle;
    gretl_bundle *b;
    int i, err = 0;

#if CDEBUG
    fprintf(stderr, "handle_plugin_call: function '%s'\n", u->name);
#endif

    if (u->pkg == NULL) {
	return E_DATA;
    }

    cfunc = get_packaged_C_function(u->pkg->name, u->name, &handle);

    if (cfunc == NULL) {
	fputs(I_("Couldn't load plugin function\n"), stderr);
	return E_FOPEN;
    }    

    b = gretl_bundle_new();

    if (b == NULL) {
	close_plugin(handle);
	return E_ALLOC;
    }

    for (i=0; i<args->argc && !err; i++) {
	fn_param *fp = &u->params[i];
	const char *key = fp->name;
	struct fnarg *arg = args->arg[i];

	if (!type_can_be_bundled(fp->type)) {
	    fprintf(stderr, "type %d: cannot be bundled\n", fp->type);
	    err = E_TYPES;
	    break;
	}

	if (gretl_scalar_type(fp->type)) {
	    double x;

	    if (arg->type == GRETL_TYPE_USCALAR) {
		x = gretl_scalar_get_value_by_index(arg->val.idnum);
	    } else if (arg->type == GRETL_TYPE_NONE) {
		x = fp->deflt;
	    } else if (arg->type == GRETL_TYPE_MATRIX) {
		x = arg->val.m->val[0];
	    } else {
		x = arg->val.x;
	    }
	    if (fp->type == GRETL_TYPE_INT || fp->type == GRETL_TYPE_OBS) {
		x = (int) x;
	    } else if (fp->type == GRETL_TYPE_BOOL) {
		x = (x != 0.0);
	    }
	    err = gretl_bundle_set_data(b, key, &x, GRETL_TYPE_DOUBLE, 0);
	} else if (fp->type == GRETL_TYPE_MATRIX) {
	    gretl_matrix *m = arg->val.m;
	    
	    err = gretl_bundle_set_data(b, key, m, fp->type, 0);
	} else if (fp->type == GRETL_TYPE_MATRIX_REF) {
	    gretl_matrix *m = user_matrix_get_matrix(arg->val.um);

	    err = gretl_bundle_set_data(b, key, m, fp->type, 0);
	} else if (fp->type == GRETL_TYPE_SERIES) {
	    int size = sample_size(dset);
	    double *px;

	    if (arg->type == GRETL_TYPE_USERIES) {
		px = dset->Z[arg->val.idnum];
	    } else {
		px = arg->val.px;
	    }
	    err = gretl_bundle_set_data(b, key, px + dset->t1, 
					GRETL_TYPE_SERIES, size);
	} else {
	    /* FIXME strings and maybe other types */
	    err = E_TYPES;
	}

#if CDEBUG
	fprintf(stderr, "arg[%d] (\"%s\") type = %d, err = %d\n", 
		i, key, fp->type, err);
#endif
    }

    if (!err) {
	err = (*cfunc) (b, prn);
    }

    close_plugin(handle);

    if (!err && u->rettype != GRETL_TYPE_VOID && retval != NULL) {
	GretlType type;
	int size;
	void *ptr;

	ptr = gretl_bundle_get_data(b, "retval", &type, &size, &err);

	if (!err) {
#if CDEBUG
	    fprintf(stderr, "%s: got return value of type %d\n", 
		    u->name, type);
#endif
	    if (type != u->rettype) {
		err = E_TYPES;
	    }
	}
	if (!err) {
	    if (type == GRETL_TYPE_DOUBLE) {
		double *px = ptr;
		
		*(double *) retval = *px;
	    } else if (type == GRETL_TYPE_MATRIX) {
		gretl_matrix *m = ptr;

		/* FIXME when to copy, when not to */
		*(gretl_matrix **) retval = gretl_matrix_copy(m);
		if (*(gretl_matrix **) retval == NULL) {
		    err = E_ALLOC;
		}
	    } else if (type == GRETL_TYPE_MATRIX_REF) {
		/* is this always right? ever right? */
		*(gretl_matrix **) retval = ptr;
	    }
	}
    }

    gretl_bundle_destroy(b);

    return err;
}

int gretl_function_exec (ufunc *u, fnargs *args, int rtype,
			 DATASET *dset, void *ret, 
			 char **descrip, PRN *prn)
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
    int indent0 = 0, started = 0;
    int retline = -1;
    int debugging = u->debug;
    int i, err = 0;

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: starting %s\n", u->name);
#endif

    err = maybe_check_function_needs(dset, u);
    if (err) {
	return err;
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: argc = %d\n", args->argc);
    fprintf(stderr, "u->n_params = %d\n", u->n_params);
#endif

    if (dset != NULL) {
	orig_v = dset->v;
	orig_t1 = dset->t1;
	orig_t2 = dset->t2;
    }

    if (!function_is_plugin(u)) {
	/* precaution */
	fn_state_init(&cmd, &state, &indent0);
	call = fncall_new(u);
	if (call == NULL) {
	    fprintf(stderr, "fncall_new() returned NULL\n");
	    return E_ALLOC;
	}
    }

    err = check_function_args(u, args, prn);

    if (function_is_plugin(u)) {
	if (!err) {
	    err = handle_plugin_call(u, args, dset, ret, prn);
	}
	return err;
    }

    if (!err) {
	call->args = args;
	err = allocate_function_args(call, dset);
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
	if (dset != NULL && dset->submask != NULL) {
	    state.submask = copy_datainfo_submask(dset, &err);
	}
	state.callback = func_exec_callback;
    }

    if (!err) {
	err = start_fncall(call, dset, prn);
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

	err = maybe_exec_line(&state, dset);

	if (!err && state.cmd->ci == FUNCRET) {
	    err = handle_return_statement(call, &state, dset);
	    if (i < u->n_lines - 1) {
		retline = i;
	    }
	    break;
	}

	if (debugging > 1 && do_debugging(&state)) {
	    if (put_func != NULL) {
		pprintf(prn, "-- debugging %s, line %d --\n", u->name, i + 1);
		(*put_func)(&state);
	    } else {
		pprintf(prn, "-- debugging %s, line %d --\n", u->name, i + 1);
	    }
	    debugging = debug_command_loop(&state, dset,
					   get_line, put_func, 
					   &err);
	}

	if (err) {
	    set_function_error_message(err, u, &state, line);
	    break;
	}

	if (state.cmd->ci == SETOBS) {
	    /* set flag for reverting on exit */
	    call->obs.changed = 1;
	}

	if (gretl_execute_loop()) { 
	    err = gretl_loop_exec(&state, dset);
	    if (err) {
		set_function_error_message(err, u, &state, state.line);
	    }
	}
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: %s: finished main exec, "
	    "err = %d, dset->v = %d\n", u->name, err, 
	    (dset != NULL)? dset->v : 0);
#endif

    if (dset != NULL) {
	/* restore the sample that was in place on entry */
	if (complex_subsampled()) {
	    if (state.submask == NULL) {
		/* we were not sub-sampled on entry */
		restore_full_sample(dset, NULL);
	    } else if (submask_cmp(state.submask, dset->submask)) {
		/* we were sub-sampled differently on entry */
		restore_full_sample(dset, NULL);
		restrict_sample_from_mask(state.submask, dset, OPT_NONE);
	    } 
	}
	dset->t1 = orig_t1;
	dset->t2 = orig_t2;
    }

    if (err) {
	if (gretl_function_depth() == 1) {
	    gretl_if_state_clear();
	} else {
	    gretl_if_state_reset(indent0);
	}
    } else if (retline >= 0) {
	/* we returned prior to the end of the function */
	gretl_if_state_reset(indent0);
    } else {
	err = gretl_if_state_check(indent0);
    }

    function_assign_returns(call, args, rtype, dset, ret, 
			    descrip, prn, &err);

    gretl_exec_state_clear(&state);

    if (started) {
	int stoperr = stop_fncall(call, rtype, ret, dset, prn,
				  orig_v);

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
		if (args->arg[i]->upname == NULL) {
		    ret = gretl_strdup("");
		} else {
		    ret = gretl_strdup(args->arg[i]->upname); 
		}
		if (ret == NULL) {
		    *err = E_ALLOC;
		}
		break;
	    }
	}
    }

    return ret;
}

/**
 * object_is_const:
 * @name: name of object (e.g. matrix).
 *
 * Checks whether the named object currently has 'const' status,
 * by virtue of its being made available as a const argument
 * to a user-defined function.
 *
 * Returns: non-zero if the object is const, 0 if it is not.
 */

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

/**
 * object_is_function_arg:
 * @name: name of object (e.g. matrix).
 *
 * Checks whether the named object has been made available
 * as a function argument.
 *
 * Returns: 1 if the object has arg status, 0 otherwise.
 */

int object_is_function_arg (const char *name)
{
    fncall *call = current_function_call();

    if (call != NULL && call->args != NULL) {
	fnargs *args = call->args;
	int i, n = args->argc;

	for (i=0; i<n; i++) {
	    const char *aname = args->arg[i]->name;

	    if (aname != NULL && !strcmp(name, aname)) {
		return 1;
	    }
	}
    }

    return 0;
}

/**
 * sample_range_get_extrema:
 * @dset: dataset info.
 * @t1: location to receive earliest possible starting observation.
 * @t2: location to receive latest possible ending observation.
 *
 * Fills out @t1 and @t2, making allowance for the possibility
 * that we're currently executing a function, on entry to 
 * which the sample range was restricted: within the function,
 * we are not allowed to overstep the bounds set on entry. 
 */

void sample_range_get_extrema (const DATASET *dset, int *t1, int *t2)
{
    fncall *call = current_function_call();

    if (call != NULL) {
	*t1 = call->obs.t1;
	*t2 = call->obs.t2;
    } else {
	*t1 = 0;
	*t2 = dset->n - 1;
    }
}

/**
 * gretl_functions_cleanup:
 * 
 * For internal use: frees all resources associated with
 * user-defined functions and function packages.
 */

void gretl_functions_cleanup (void)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	ufunc_free(ufuns[i]);
    }

    free(ufuns);
    ufuns = NULL;
    n_ufuns = 0;

    for (i=0; i<n_pkgs; i++) {
	function_package_free(pkgs[i]);
    }

    free(pkgs);
    pkgs = NULL;
    n_pkgs = 0;
}

/* generate help output for a packaged function, either for
   display on the console, or with markup for display in a
   GtkTextView window in the GUI (opt & OPT_M)
*/

static void real_user_function_help (ufunc *fun, gretlopt opt, PRN *prn)
{
    fnpkg *pkg = fun->pkg;
    int markup = (opt & OPT_M);
    int i;

    if (markup) {
	pprintf(prn, "<@hd1=\"%s\">\n\n", fun->name);
    } else {
	pprintf(prn, "%s\n\n", fun->name);
    }

    if (pkg != NULL) {
	if (markup) {
	    pprintf(prn, "<@hd1=\"Author\">: %s\n", pkg->author? pkg->author : "unknown");
	    pprintf(prn, "<@hd1=\"Version\">: %s (%s)\n\n", pkg->version? pkg->version : "unknown",
		    pkg->date? pkg->date : "unknown");
	} else {
	    pprintf(prn, "Author: %s\n", pkg->author? pkg->author : "unknown");
	    pprintf(prn, "Version: %s (%s)\n\n", pkg->version? pkg->version : "unknown",
		    pkg->date? pkg->date : "unknown");
	}
    }

    if ((opt & OPT_G) && pkg != NULL && pkg->gui_help != NULL) {
	/* GUI-specific help is preferred, and is available */
	if (markup) {
	    pputs(prn, "<@hd1=\"Help text\">:\n\n");
	    pputs(prn, "<mono>\n");
	} else {
	    pputs(prn, "Help text:\n");
	}   	
	pputs(prn, pkg->gui_help);
	if (markup) {
	    pputs(prn, "\n</mono>");
	}
	pputs(prn, "\n\n");
	return;
    }

    if (markup) {
	pputs(prn, "<@hd1=\"Parameters\">: ");
    } else {
	pputs(prn, "Parameters: ");
    }    

    if (fun->n_params > 0) {
	pputc(prn, '\n');
	for (i=0; i<fun->n_params; i++) {
	    pprintf(prn, " %s (%s)\n", 
		    fun->params[i].name, gretl_arg_type_name(fun->params[i].type));
	}
	pputc(prn, '\n');
    } else {
	pputs(prn, "none\n\n");
    }

    if (markup) {
	pputs(prn, "<@hd1=\"Return value\">: ");
    } else {
	pputs(prn, "Return value: ");
    }      

    if (fun->rettype != GRETL_TYPE_NONE && fun->rettype != GRETL_TYPE_VOID) {
	pprintf(prn, "%s\n\n", gretl_arg_type_name(fun->rettype));
    } else {
	pputs(prn, "none\n\n");
    }

    if (pkg != NULL && pkg->help != NULL) {
	if (markup) {
	    pputs(prn, "<@hd1=\"Help text\">:\n\n");
	    pputs(prn, "<mono>\n");
	} else {
	    pputs(prn, "Help text:\n");
	}   	
	if (has_suffix(pkg->help, ".pdf")) {
	    const char *s = strrchr(pkg->help, ':');

	    if (s != NULL) {
		if (markup) {
		    pprintf(prn, "See <@pdf=\"%s\">", s + 1);
		} else {
		    pprintf(prn, "See %s", s + 1);
		}
	    } else {
		pputs(prn, pkg->help);
	    }
	} else {	
	    pputs(prn, pkg->help);
	}
	if (markup) {
	    pputs(prn, "\n</mono>");
	}
	pputs(prn, "\n\n");
    }

    if (pkg != NULL && pkg->sample != NULL) {
	if (markup) {
	    pputs(prn, "<@hd1=\"Sample script\">:\n\n");
	    pputs(prn, "<code>\n");
	} else {
	    pputs(prn, "Sample script:\n\n");
	}
	pputs(prn, pkg->sample);
	if (markup) {
	    pputs(prn, "\n</code>\n");
	} else {
	    pprintf(prn, "\n\n");
	}
    }	
}

/**
 * user_function_help:
 * @fnname: name of function.
 * @opt: may include OPT_M for adding markup, OPT_G
 * for preferring GUI-specific help, if available.
 * @prn: printing struct.
 * 
 * Looks for a function named @fnname and prints
 * as much help information as possible.  
 *
 * Returns: 0 on success, non-zero if the function is not
 * found.
 */

int user_function_help (const char *fnname, gretlopt opt, PRN *prn)
{
    ufunc *fun = get_user_function_by_name(fnname);
    int err = 0;

    if (fun == NULL) {
	pprintf(prn, _("\"%s\" is not defined.\n"), fnname);
	err = 1;
    } else {
	real_user_function_help(fun, opt, prn);
    }

    return err;
}

/**
 * function_package_has_PDF_doc:
 * @pkg: function-package pointer.
 * @pdfname: location to receive basename of PDF file,
 * if applicable (or NULL if this is not wanted).
 * 
 * Checks whether @pkg is documented in the form of a PDF file.
 *
 * Returns: 1 if so (in which case the name of that file is
 * returned via @pdfname), 0 otherwise.
 */

int function_package_has_PDF_doc (fnpkg *pkg, char **pdfname)
{
    int ret = 0;

    if (pkg->help != NULL && !strncmp(pkg->help, "pdfdoc:", 7)) {
	ret = 1;
	if (pdfname != NULL) {
	    *pdfname = switch_ext_new(pkg->fname, "pdf");
	    if (*pdfname == NULL) {
		ret = 0;
	    }
	}
    }

    return ret;
}


