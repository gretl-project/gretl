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
#include "version.h"
#include "monte_carlo.h"
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_xml.h"
#include "cmd_private.h"
#include "gretl_string_table.h"
#include "gretl_typemap.h"
#include "gretl_zip.h"
#include "uservar.h"
#include "flow_control.h"
#include "system.h"
#include "genr_optim.h"
#include "gretl_foreign.h"
#include "gretl_mdconv.h"
#include "gen_public.h"
#include "addons_utils.h"
#include "gfn_translations.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#include <errno.h>
#include <glib.h>
#include <glib/gstdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define FNPARSE_DEBUG 0 /* debug parsing of function code */
#define EXEC_DEBUG 0    /* debugging of function execution */
#define ARGS_DEBUG 0    /* debug handling of args */
#define PKG_DEBUG 0     /* debug handling of function packages */
#define FN_DEBUG 0      /* miscellaneous debugging */
#define COMP_DEBUG 0    /* debug "compilation" */
#define CALL_DEBUG 0    /* debug handling of the callstack */
#define REC_DEBUG 0     /* debug recursion */

#define COMPILE_RECURSIVE 0 /* May 2024: too risky */

#define INT_USE_XLIST (-999)
#define INT_USE_MYLIST (-777)

typedef struct fn_param_ fn_param;
typedef struct fn_arg_ fn_arg;
typedef struct obsinfo_ obsinfo;

/* local symbol for structure representing a line of a
   user-defined function */
typedef struct stmt_ fn_line;

enum {
    FP_CONST    = 1 << 0, /* explicitly marked as "const" */
    FP_OPTIONAL = 1 << 1, /* marked as optional (null is OK) */
    FP_AUTO     = 1 << 2  /* marked as automatic (has a default value) */
};

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

#define param_is_const(p) (p->flags & FP_CONST)
#define param_is_optional(p) (p->flags & FP_OPTIONAL)
#define param_is_auto(p) (p->flags & FP_AUTO)

#define set_param_optional(p) (p->flags |= FP_OPTIONAL)
#define set_param_auto(p) (p->flags |= (FP_OPTIONAL | FP_AUTO))

typedef enum {
    LINE_IGNORE = 1 << 0,
    LINE_NOCOMP = 1 << 1  /* "don't compile genr" */
} ln_flags;

#define line_has_loop(l)  (l->ci == LOOP && l->ptr != NULL)
#define ignore_line(l)    (l->flags & LINE_IGNORE)
#define line_no_comp(l)   (l->flags & LINE_NOCOMP)

#define UNSET_VALUE (-1.0e200)
#define default_unset(p) (p->deflt == UNSET_VALUE)

/* structure representing sample information at start of
   a function call */

struct obsinfo_ {
    int structure;      /* time-series, etc. */
    int pd;             /* data frequency */
    int t1, t2;         /* starting and ending observations */
    int added;          /* number of observations added within function */
    char stobs[OBSLEN]; /* string representation of starting obs */
    int panel_pd;       /* panel time frequency, if applicable */
    double panel_sd0;   /* panel time starting point */
};

/* structure representing a call to a user-defined function */

typedef enum {
    FC_RECURSING = 1 << 0, /* call to f() has a call to f() in its ancestry */
    FC_PREV_MSGS = 1 << 1, /* record of "messages" setting before call */
    FC_PREV_ECHO = 1 << 2, /* record of "echo" setting before call */
    FC_PRESERVE  = 1 << 3  /* function call should be preserved */
} FCFlags;

/* Note: the FC_PRESERVE flag is set when we want to allow for repeated
   calls to a given function. In the non-recursive case we just reuse
   the stored fncall struct. In the recursive case we need to ensure
   there's a distinct call struct for each depth of function execution,
   and the function then acquires a GList of fncalls.
*/

struct linegen_ {
    int idx;
    GENERATOR *genr;
};

typedef struct linegen_ linegen;

struct fncall_ {
    ufunc *fun;      /* the function called */
    int argc;        /* argument count */
    int orig_v;      /* number of series defined on entry */
    fn_arg *args;    /* argument array */
    int *ptrvars;    /* list of pointer arguments */
    int *listvars;   /* list of series included in a list argument */
    GList *lists;    /* list of names of list arguments */
    char *retname;   /* name of return value (or dummy string) */
    GretlType rtype; /* return type (when not fixed in advance) */
    obsinfo obs;     /* sample info */
    FCFlags flags;   /* indicators for recursive call, etc */
    fn_line *line;   /* currently executing line */
    linegen *lgen;   /* array of per-line genrs */
    int n_lgen;      /* number of linegens */
};

/* structure representing a user-defined function */

struct ufunc_ {
    char name[FN_NAMELEN]; /* identifier */
    fnpkg *pkg;            /* pointer to parent package, or NULL */
    fncall *call;          /* pointer to current call, or NULL */
    GList *calls;          /* for use in recursion */
    int pkg_role;          /* printer, plotter, etc. */
    UfunAttrs flags;       /* private, plugin, etc. */
    int line_idx;          /* current line index (compiling) */
    int n_lines;           /* number of lines of code */
    fn_line *lines;        /* array of lines of code */
    int n_params;          /* number of parameters */
    fn_param *params;      /* parameter info array */
    int rettype;           /* return type (if any) */
};

/* structure representing a function package */

struct fnpkg_ {
    char name[FN_NAMELEN]; /* package name */
    char *fname;      /* filename */
    char *author;     /* author's name */
    char *email;      /* author's email address */
    char *version;    /* package version string */
    char *date;       /* last revision date */
    char *descrip;    /* package description */
    char *help;       /* package help text */
    char *gui_help;   /* GUI-specific help (optional) */
    char *Rdeps;      /* R dependencies (if any) */
    char *sample;     /* sample caller script */
    char *help_fname;     /* filename: package help text */
    char *gui_help_fname; /* filename: GUI-specific help text */
    char *sample_fname;   /* filename: sample caller script */
    char *tags;       /* tag string(s) */
    char *label;      /* for use in GUI menus */
    char *mpath;      /* menu path in GUI */
    int minver;       /* minimum required gretl version */
    guint8 uses_subdir; /* lives in subdirectory (0/1) */
    guint8 prechecked;  /* already checked for data requirement */
    guint8 data_access; /* wants access to full data range */
    DataReq dreq;     /* data requirement */
    int modelreq;     /* required model type, if applicable */
    ufunc **pub;      /* pointers to public interfaces */
    ufunc **priv;     /* pointers to private functions */
    int n_pub;        /* number of public functions */
    int n_priv;       /* number of private functions */
    guint8 overrides; /* number of overrides of built-in functions */
    char **datafiles; /* names of packaged data files */
    char **depends;   /* names of dependencies */
    char *provider;   /* name of "provider" package, if applicable */
    int n_files;      /* number of data files */
    int n_depends;    /* number of dependencies */
    void *editor;     /* for GUI use */
    void *trans;      /* translations */
};

/* acceptable types for parameters of user-defined functions */

#define ok_function_arg_type(t) (t == GRETL_TYPE_BOOL ||                \
                                 t == GRETL_TYPE_INT ||                 \
                                 t == GRETL_TYPE_OBS ||                 \
                                 t == GRETL_TYPE_DOUBLE ||              \
                                 t == GRETL_TYPE_SERIES ||              \
                                 t == GRETL_TYPE_LIST ||                \
                                 t == GRETL_TYPE_MATRIX ||              \
                                 t == GRETL_TYPE_STRING ||              \
                                 t == GRETL_TYPE_BUNDLE ||              \
                                 t == GRETL_TYPE_SCALAR_REF ||          \
                                 t == GRETL_TYPE_SERIES_REF ||          \
                                 t == GRETL_TYPE_MATRIX_REF ||          \
                                 t == GRETL_TYPE_BUNDLE_REF ||          \
                                 t == GRETL_TYPE_STRING_REF ||          \
                                 t == GRETL_TYPE_STRINGS ||             \
                                 t == GRETL_TYPE_MATRICES ||            \
                                 t == GRETL_TYPE_BUNDLES||              \
                                 t == GRETL_TYPE_LISTS ||               \
                                 t == GRETL_TYPE_ARRAYS ||              \
                                 t == GRETL_TYPE_STRINGS_REF ||         \
                                 t == GRETL_TYPE_MATRICES_REF ||        \
                                 t == GRETL_TYPE_BUNDLES_REF ||         \
                                 t == GRETL_TYPE_LISTS_REF ||           \
                                 t == GRETL_TYPE_ARRAYS_REF ||          \
                                 t == GRETL_TYPE_NUMERIC)

/* structure representing an argument to a user-defined function */

struct fn_arg_ {
    char type;            /* argument type */
    char shifted;         /* level was shifted for execution */
    char *upname;         /* name of supplied arg at caller level */
    user_var *uvar;       /* reference to "parent", if any */
    union {
        int idnum;        /* named series arg (series ID) */
        double x;         /* scalar arg */
        double *px;       /* anonymous series arg */
        gretl_matrix *m;  /* matrix arg */
        char *str;        /* string arg */
        int *list;        /* list arg */
        gretl_bundle *b;  /* anonymous bundle pointer */
        gretl_array *a;   /* array argument */
    } val;
};

static int n_ufuns;         /* number of user-defined functions in memory */
static ufunc **ufuns;       /* array of pointers to user-defined functions */
static ufunc *current_fdef; /* pointer to function currently being defined */
static GList *callstack;    /* stack of function calls */
static int n_pkgs;          /* number of loaded function packages */
static fnpkg **pkgs;        /* array of pointers to loaded packages */
static fnpkg *current_pkg;  /* pointer to package currently being edited */

static int function_package_record (fnpkg *pkg);
static void function_package_free (fnpkg *pkg);
static int load_function_package (const char *fname,
                                  gretlopt opt,
                                  GArray *pstack,
                                  PRN *prn,
                                  fnpkg **ppkg,
                                  int level);
static int ufunc_get_structure (ufunc *u);
#if CALL_DEBUG
static void print_callstack (ufunc *fun);
#endif

/* record of state, and communication of state with outside world */

static int compiling;    /* boolean: are we compiling a function currently? */
static int fn_executing; /* depth of function call stack */
static int compiling_python;
#ifdef HAVE_MPI
static char mpi_caller[FN_NAMELEN];
#endif

#define function_is_private(f)   (f->flags & UFUN_PRIVATE)
#define function_is_noprint(f)   (f->flags & UFUN_NOPRINT)
#define function_is_menu_only(f) (f->flags & UFUN_MENU_ONLY)
#define function_is_recursive(f) (f->flags & UFUN_RECURSES)

#define set_call_recursing(c)   (c->flags |= FC_RECURSING)
#define is_recursing(c)         (c->flags & FC_RECURSING)

struct flag_and_key {
    int flag;
    const char *key;
};

static struct flag_and_key pkg_lookups[] = {
    { UFUN_BUNDLE_PRINT,  BUNDLE_PRINT },
    { UFUN_BUNDLE_PLOT,   BUNDLE_PLOT },
    { UFUN_BUNDLE_TEST,   BUNDLE_TEST },
    { UFUN_BUNDLE_FCAST,  BUNDLE_FCAST },
    { UFUN_BUNDLE_EXTRA,  BUNDLE_EXTRA },
    { UFUN_GUI_MAIN,      GUI_MAIN },
    { UFUN_GUI_PRECHECK,  GUI_PRECHECK },
    { UFUN_PLOT_PRECHECK, PLOT_PRECHECK },
    { UFUN_LIST_MAKER,    LIST_MAKER },
    { UFUN_R_SETUP,       R_SETUP },
    { UFUN_UI_MAKER,      UI_MAKER },
    { -1,                 NULL }
};

#define pkg_aux_role(r) (r == UFUN_BUNDLE_PRINT || \
                         r == UFUN_BUNDLE_PLOT ||  \
                         r == UFUN_BUNDLE_TEST ||  \
                         r == UFUN_BUNDLE_FCAST || \
                         r == UFUN_BUNDLE_EXTRA || \
                         r == UFUN_R_SETUP || \
                         r == UFUN_UI_MAKER)

static int pkg_key_get_role (const char *key)
{
    int i;

    if (key != NULL && *key != '\0') {
        for (i=0; pkg_lookups[i].flag > 0; i++) {
            if (!strcmp(key, pkg_lookups[i].key)) {
                return pkg_lookups[i].flag;
            }
        }
    }

    return UFUN_ROLE_NONE;
}

const char *package_role_get_key (int flag)
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

int gretl_compiling_python (const char *line)
{
    if (compiling_python) {
        char s1[4], s2[8];

        if (sscanf(line, "%3s %7s", s1, s2) == 2 &&
            !strcmp(s1, "end") && !strcmp(s2, "foreign")) {
            compiling_python = 0;
        }
    }

    return compiling_python;
}

static void set_compiling_on (void)
{
    compiling = 1;
}

static void set_compiling_off (void)
{
    compiling = compiling_python = 0;
}

int gretl_function_depth (void)
{
    return fn_executing;
}

static void adjust_array_arg_type (fn_arg *arg)
{
    GretlType t = gretl_array_get_type(arg->val.a);

    if (arg->type == GRETL_TYPE_ARRAY_REF) {
        arg->type = gretl_type_get_ref_type(t);
    } else {
        arg->type = t;
    }
}

static void nullify_upname (fn_arg *arg)
{
    if (arg != NULL && arg->upname != NULL) {
        free(arg->upname);
        arg->upname = NULL;
    }
}

static int fn_arg_set_data (fn_arg *arg, const char *name,
                            user_var *uvar, GretlType type,
                            void *p)
{
    int err = 0;

    arg->type = type;
    arg->shifted = 0;
    arg->uvar = uvar;
    arg->upname = (name != NULL)? gretl_strdup(name) : NULL;

    if (type == GRETL_TYPE_NONE) {
        arg->val.x = 0;
    } else if (type == GRETL_TYPE_DOUBLE ||
               type == GRETL_TYPE_SCALAR_REF) {
        arg->val.x = *(double *) p;
    } else if (type == GRETL_TYPE_INT ||
               type == GRETL_TYPE_OBS) {
        arg->val.x = *(int *) p;
    } else if (type == GRETL_TYPE_SERIES) {
        arg->val.px = (double *) p;
    } else if (type == GRETL_TYPE_MATRIX ||
               type == GRETL_TYPE_MATRIX_REF) {
        arg->val.m = (gretl_matrix *) p;
    } else if (type == GRETL_TYPE_STRING ||
               type == GRETL_TYPE_STRING_REF) {
        arg->val.str = (char *) p;
    } else if (type == GRETL_TYPE_LIST) {
        arg->val.list = (int *) p;
    } else if (type == GRETL_TYPE_SERIES_REF ||
               type == GRETL_TYPE_USERIES) {
        arg->val.idnum = *(int *) p;
    } else if (type == GRETL_TYPE_BUNDLE ||
               type == GRETL_TYPE_BUNDLE_REF) {
        arg->val.b = (gretl_bundle *) p;
    } else if (type == GRETL_TYPE_ARRAY ||
               type == GRETL_TYPE_ARRAY_REF) {
        arg->val.a = (gretl_array *) p;
        adjust_array_arg_type(arg);
    } else {
        err = E_TYPES;
    }

    return err;
}

static int fncall_add_args_array (fncall *fc)
{
    int i, np = fc->fun->n_params;
    int err = 0;

    fc->args = malloc(np * sizeof *fc->args);

    if (fc->args == NULL) {
        err = E_ALLOC;
    } else {
        for (i=0; i<np; i++) {
            fc->args[i].type = 0;
            fc->args[i].shifted = 0;
            fc->args[i].upname = NULL;
            fc->args[i].uvar = NULL;
        }
    }

    return err;
}

static void fncall_clear_args_array (fncall *fc)
{
    int i, np = fc->fun->n_params;
    fn_arg *arg;

    for (i=0; i<np; i++) {
        arg = &fc->args[i];
        arg->type = 0;
        arg->shifted = 0;
        nullify_upname(arg);
        arg->uvar = NULL;
        arg->val.px = NULL;
    }

    fc->argc = 0;
}

static void fncall_destroy_args_array (fncall *fc)
{
    if (fc->args != NULL) {
	int i;

	for (i=0; i<fc->fun->n_params; i++) {
	    if (fc->args[i].upname != NULL) {
		free(fc->args[i].upname);
	    }
	}
	free(fc->args);
    }
}

static void maybe_set_param_const (fn_param *fp)
{
    fp->flags |= FP_CONST;
}

/**
 * push_function_arg:
 * @fc: pointer to function call.
 * @name: name of variable (or NULL for anonymous).
 * @uvar: reference to user_var or NULL.
 * @type: type of argument to add.
 * @value: pointer to value to add.
 *
 * Writes a new argument of the specified type and value into the
 * argument array of @fc.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_function_arg (fncall *fc, const char *name,
                       void *uvar, GretlType type,
                       void *value)
{
    int err = 0;

    if (fc == NULL || fc->fun == NULL) {
        err = E_DATA;
    } else if (fc->argc >= fc->fun->n_params) {
        fprintf(stderr, "function %s has %d parameters but argc = %d\n",
                fc->fun->name, fc->fun->n_params, fc->argc);
        err = E_DATA;
    } else if (fc->args == NULL) {
        err = fncall_add_args_array(fc);
    }

    if (!err) {
        err = fn_arg_set_data(&fc->args[fc->argc], name, uvar, type, value);
        fc->argc += 1;
    }

    return err;
}

/**
 * push_anon_function_arg:
 * @fc: pointer to function call.
 * @type: type of argument to add.
 * @value: pointer to value to add.
 *
 * Writes a new argument of the specified type and value into the
 * argument array of @fc.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_anon_function_arg (fncall *fc, GretlType type, void *value)
{
    return push_function_arg(fc, NULL, NULL, type, value);
}

/**
 * set_anon_function_arg:
 * @fc: pointer to function call.
 * @i: index of argument.
 * @type: type of argument.
 * @value: pointer to value to set.
 *
 * Sets type and value of (anonymous) argument @i in @fc. The index @i
 * must be >= 0 and less than the argc member of @fc.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int set_anon_function_arg (fncall *fc, int i, GretlType type, void *value)
{
    if (i >= 0 && i < fc->argc) {
        fn_arg_set_data(&fc->args[i], NULL, NULL, type, value);
        return 0;
    } else {
        return E_DATA;
    }
}

/**
 * push_function_args:
 * @fc: pointer to function call.
 *
 * Writes multiple entries into the argument array of @fc.
 * Each argument must be given in the form {type, value},
 * The list of entries must be terminated with -1.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int push_function_args (fncall *fc, ...)
{
    va_list ap;
    int argtype;
    void *value;
    int err = 0;

    va_start(ap, fc);
    while (err == 0) {
        argtype = va_arg(ap, int);
        if (argtype < 0) {
            /* reached the end of the args */
            break;
        }
        value = va_arg(ap, void *);
        err = push_function_arg(fc, NULL, NULL, argtype, value);
    }
    va_end(ap);

    return err;
}

/* fncall_new: if @preserve is non-zero, we don't destroy the fncall
   struct automatically on completion of function execution -- in
   which case the caller of this function is responsible for calling
   fncall_destroy() at a suitable point in the proceedings.
*/

fncall *fncall_new (ufunc *fun, int preserve)
{
    fncall *call = malloc(sizeof *call);

    if (call != NULL) {
        call->fun = fun;
        call->ptrvars = NULL;
        call->listvars = NULL;
        call->lists = NULL;
        call->retname = NULL;
        call->rtype = fun->rettype;
        call->argc = 0;
        call->args = NULL;
        call->flags = preserve ? FC_PRESERVE : 0;
        call->line = NULL;
        call->lgen = NULL;
        call->n_lgen = 0;
    }

    return call;
}

static int call_is_in_use (fncall *fc)
{
    if (callstack != NULL) {
        GList *tmp = g_list_last(callstack);

        while (tmp != NULL) {
            if (fc == tmp->data) {
                return 1;
            }
            tmp = tmp->prev;
        }
    }

    return 0;
}

static void clear_fncall_data (fncall *fc)
{
    if (fc->args != NULL) {
        fncall_clear_args_array(fc);
    }
    if (fc->ptrvars != NULL) {
        free(fc->ptrvars);
        fc->ptrvars = NULL;
    }
    if (fc->listvars != NULL) {
        free(fc->listvars);
        fc->listvars = NULL;
    }
    if (fc->lists != NULL) {
        g_list_free(fc->lists);
        fc->lists = NULL;
    }
    free(fc->retname);
    fc->retname = NULL;
    fc->line = NULL;
}

/* This function is called (only) from geneval.c, in the function
   eval_ufunc(). The primary role of the fncall struct is to serve as
   a vehicle for conveying the arguments to the user function in
   question, but it also carries state information about the call,
   including whether or not it is subject to recursion.
*/

fncall *user_func_get_fncall (ufunc *fun)
{
    GList *tmp = g_list_first(fun->calls);
    fncall *fc = NULL;

    while (tmp != NULL) {
        if (!call_is_in_use(tmp->data)) {
            fc = tmp->data;
            break;
        }
        tmp = tmp->next;
    }

    if (fc != NULL) {
        clear_fncall_data(fc);
        fun->call = fc;
    } else {
        fun->call = fncall_new(fun, 1);
        fun->calls = g_list_append(fun->calls, fun->call);
    }

#if REC_DEBUG
    fprintf(stderr, "call for %s, depth %d: %p (n_lgen %d)\n", fun->name,
            gretl_function_depth(), (void *) fun->call, fun->call->n_lgen);
#endif

    return fun->call;
}

void fncall_destroy (void *data)
{
    if (data != NULL) {
#if CALL_DEBUG
        fprintf(stderr, "destroying fncall at %p\n", data);
#endif
        fncall *call = data;

        fncall_destroy_args_array(call);
        free(call->ptrvars);
        free(call->listvars);
        free(call->retname);
        g_list_free(call->lists);
        if (call->n_lgen > 0) {
            int i;

            for (i=0; i<call->n_lgen; i++) {
                destroy_genr(call->lgen[i].genr);
            }
            free(call->lgen);
        }
        free(call);
    }
}

static void maybe_destroy_fncall (fncall **pcall)
{
    fncall *call = *pcall;

    if (call != NULL && !(call->flags & FC_PRESERVE)) {
        if (call->fun != NULL) {
            if (call->fun->calls != NULL) {
                call->fun->calls = g_list_remove(call->fun->calls, call);
            }
            if (call->fun->call == call) {
                call->fun->call = NULL;
            }
        }
        fncall_destroy(call);
        *pcall = NULL;
    }
}

/* For use with overloaded functions whose return type is
   'numeric'. The specific return type is not known
   until execution is completed; we can supply it here.
*/

GretlType fncall_get_return_type (fncall *call)
{
    if (call != NULL) {
        return call->rtype;
    }

    return GRETL_TYPE_NONE;
}

/* Portmanteau function to get a caller struct for a function named
   @funcname from a function package named @pkgname. We first check if
   the specified package is already in memory; if not we try to find
   it in the local filesystem, and load it into memory if successful.

   If/once the package is in fact loaded we look up the specified
   function; and if that's successful we allocate a caller struct and
   return it.

   The @pkgpath argument can be given as NULL, but if the path to the
   package is known to the caller and is provided via this argument
   this will speed up the look-up in case the package is not already
   loaded.
*/

fncall *get_pkg_function_call (const char *funcname,
                               const char *pkgname,
                               const char *pkgpath)
{
    fncall *fc = NULL;
    ufunc *uf = NULL;
    fnpkg *pkg;

    /* is the package already loaded? */
    pkg = get_function_package_by_name(pkgname);

    if (pkg == NULL) {
        /* no, so look it up */
        int err = 0;

        if (pkgpath != NULL) {
            /* path was supplied by caller */
            pkg = get_function_package_by_filename(pkgpath, &err);
        } else {
            /* we need to search */
            char *mypath;

            mypath = gretl_function_package_get_path(pkgname, PKG_ALL);
            if (mypath != NULL) {
                pkg = get_function_package_by_filename(mypath, &err);
                free(mypath);
            }
        }
    }

    if (pkg != NULL) {
        uf = get_function_from_package(funcname, pkg);
    }

    if (uf == NULL) {
        gretl_errmsg_sprintf(_("Couldn't find function %s"), funcname);
    } else {
        fc = fncall_new(uf, 0); /* FIXME second arg? */
    }

    return fc;
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

#if PKG_DEBUG
    fprintf(stderr, "function_package_alloc: fname='%s'\n", fname);
#endif

    pkg->name[0] = '\0';
    pkg->author = NULL;
    pkg->email = NULL;
    pkg->version = NULL;
    pkg->date = NULL;
    pkg->descrip = NULL;
    pkg->help = NULL;
    pkg->gui_help = NULL;
    pkg->Rdeps = NULL;
    pkg->sample = NULL;
    pkg->help_fname = NULL;
    pkg->gui_help_fname = NULL;
    pkg->sample_fname = NULL;
    pkg->tags = NULL;
    pkg->label = NULL;
    pkg->mpath = NULL;
    pkg->dreq = FN_NEEDS_DATA;
    pkg->modelreq = 0;
    pkg->minver = 0;
    pkg->uses_subdir = 0;
    pkg->prechecked = 0;
    pkg->data_access = 0;

    pkg->pub = pkg->priv = NULL;
    pkg->n_pub = pkg->n_priv = 0;
    pkg->overrides = 0;
    pkg->datafiles = NULL;
    pkg->n_files = 0;
    pkg->depends = NULL;
    pkg->n_depends = 0;
    pkg->provider = NULL;
    pkg->editor = NULL;
    pkg->trans = NULL;

    return pkg;
}

/* For function call @call, set the 'listarg' status of any named
   series provided to the called function via list arguments.  This is
   required for correct handling of namespaces.
*/

static void set_listargs_from_call (fncall *call, DATASET *dset)
{
    int i, vi;

    if (dset == NULL) {
        return;
    }

    for (i=1; i<dset->v; i++) {
        series_unset_flag(dset, i, VAR_LISTARG);
    }

    if (call != NULL && call->listvars != NULL) {
        for (i=1; i<=call->listvars[0]; i++) {
            vi = call->listvars[i];
#if ARGS_DEBUG
            fprintf(stderr, "setting listarg status on var %d (%s)\n",
                    vi, dset->varname[vi]);
#endif
            series_set_flag(dset, vi, VAR_LISTARG);
        }
    }
}

static void get_prior_call_for_function (ufunc *fun)
{
    GList *tmp = g_list_last(callstack);
    fncall *call;

    tmp = tmp->prev;
    while (tmp != NULL) {
        call = tmp->data;
        if (call->fun == fun) {
            /* reattach */
#if CALL_DEBUG
            fprintf(stderr, "reattach call at %p to function %s\n",
                    (void *) call, fun->name);
#endif
            fun->call = call;
            break;
        }
        tmp = tmp->prev;
    }
}

static void set_executing_off (fncall **pcall, DATASET *dset, PRN *prn)
{
    int dbg = gretl_debugging_on();
    fncall *thiscall = *pcall;
    fncall *popcall = NULL;
    ufunc *fun = thiscall->fun;

#if CALL_DEBUG
    fprintf(stderr, "set_executing_off, starting\n");
    print_callstack(fun);
#endif

    destroy_option_params_at_level(fn_executing);
    set_previous_depth(fn_executing);
    fn_executing--;

    if (is_recursing(thiscall)) {
        get_prior_call_for_function(fun);
    }

    callstack = g_list_remove(callstack, thiscall);

#if CALL_DEBUG || EXEC_DEBUG
    fprintf(stderr, "set_executing_off: removed call %p to %s, depth now %d\n",
            (void *) thiscall, fun->name, g_list_length(callstack));
#endif

    if (dbg) {
        pputs(prn, "*** ");
        bufspace(gretl_function_depth(), prn);
        pprintf(prn, "exiting function %s, ", fun->name);
    }

    maybe_destroy_fncall(pcall);

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

    if (popcall == NULL) {
        /* returning to main */
        switch_uservar_hash(0);
        if (dset != NULL) {
            series_ensure_level_zero(dset);
        }
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
 * n_user_functions:
 *
 * Returns: the number of hansl functions currently loaded in memory.
 */

int n_user_functions (void)
{
    return n_ufuns;
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
    const char *ret = NULL;

    if (i >= 0 && i < fun->n_params) {
        fn_param *fp = &fun->params[i];

        if (fp->descrip != NULL) {
            if (fun->pkg != NULL && fun->pkg->trans != NULL) {
                ret = get_gfn_translation(fun->pkg->trans, fp->descrip);
            } else {
                ret = fp->descrip;
            }
        }
    }

    return ret;
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

char **fn_param_value_labels (const ufunc *fun, int i, int *n)
{
    if (i >= 0 && i < fun->n_params) {
        fn_param *fp = &fun->params[i];

        if (fp->labels != NULL) {
            void *trans = fun->pkg == NULL ? NULL : fun->pkg->trans;
            char **S = calloc(fp->nlabels, sizeof *S);
            int j;

            if (S == NULL) {
                return NULL;
            }
            for (j=0; j<fp->nlabels; j++) {
                if (trans != NULL) {
                    S[j] = (char *) get_gfn_translation(fun->pkg->trans,
                                                        fp->labels[j]);
                } else {
                    S[j] = fp->labels[j];
                }
            }
            *n = fp->nlabels;
            return S;
        }
    }

    *n = 0;
    return NULL;
}

/**
 * fn_param_has_default:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if the (scalar) parameter @i of function @fun
 * (if any) has a default value set, otherwise 0.
 */

int fn_param_has_default (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
        return 0;
    } else {
        fn_param *fp = &fun->params[i];

        return !default_unset(fp);
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
    if (i < 0 || i >= fun->n_params) {
        return NADBL;
    } else {
        fn_param *fp = &fun->params[i];

        return default_unset(fp)? NADBL : fp->deflt;
    }
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

static int arg_may_be_optional (GretlType t)
{
    return gretl_ref_type(t) || gretl_is_array_type(t) ||
        t == GRETL_TYPE_SERIES ||
        t == GRETL_TYPE_MATRIX ||
        t == GRETL_TYPE_BUNDLE ||
        t == GRETL_TYPE_LIST ||
        t == GRETL_TYPE_STRING;
}

/**
 * fn_param_automatic:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if parameter @i of function @fun is marked "auto",
 * meaning that it is optional and has a default value,
 * otherwise 0.
 */

int fn_param_automatic (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
        return 0;
    }

    return arg_may_be_optional(fun->params[i].type) &&
        (fun->params[i].flags & FP_AUTO);
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
    if (i < 0 || i >= fun->n_params) {
        return 0;
    }

    return arg_may_be_optional(fun->params[i].type) &&
        (fun->params[i].flags & FP_OPTIONAL);
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
 * fn_param_uses_mylist:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if parameter @i of function @fun is
 * designed to select an integer based on a custom
 * list, constructed by the function, otherwise 0.
 */

int fn_param_uses_mylist (const ufunc *fun, int i)
{
    if (i < 0 || i >= fun->n_params) {
        return 0;
    }

    return (fun->params[i].type == GRETL_TYPE_INT &&
            fun->params[i].deflt == INT_USE_MYLIST);
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
 * user_func_is_menu_only:
 * @fun: pointer to user-function.
 *
 * Returns: 1 if the function is not designed to be called
 * other than via a GUI menu.
 */

int user_func_is_menu_only (const ufunc *fun)
{
    if (fun == NULL) {
        return 0;
    } else {
        return function_is_menu_only(fun);
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
 * @pkg: reference function package, or NULL.
 *
 * Looks up the 0-based index of the named function
 * in the current array of user-functions. If @pkg is
 * non-NULL the search for @name is confined to functions
 * that belong to @pkg; otherwise the index of the first
 * match is returned.
 *
 * Returns: 0-based index or -1 on failure.
 */

int user_function_index_by_name (const char *name,
                                 fnpkg *pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
        if ((pkg == NULL || ufuns[i]->pkg == pkg) &&
            !strcmp(name, ufuns[i]->name)) {
            return i;
        }
    }

    return -1;
}

static int fname_idx;

void function_names_init (void)
{
    fname_idx = 0;
}

/* Apparatus used in the GUI selector for composing a new function
   package, or editing an existing one.  In the first case we want a
   list of names of currently unpackaged functions (and the @pkg
   argument should be NULL); in the second the list should include
   both unpackaged functions and those belonging to the package
   in question (specified via the @pkg arg).

   The caller should first invoke function_names_init(), then keep
   calling next_available_function_name() until it returns %NULL. The
   pointer argument @idxp provides a means to grab the "index number"
   (position in the current functions array) corresponding to the
   returned function name.
*/

const char *next_available_function_name (fnpkg *pkg, int *idxp)
{
    const char *ret = NULL;
    ufunc *fun;

    if (n_ufuns == 0) {
        fname_idx = 0;
        return NULL;
    }

    while (fname_idx < n_ufuns) {
        fun = ufuns[fname_idx++];
        if (fun->pkg == NULL || fun->pkg == pkg) {
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

void current_function_info (char const **funcname,
                            char const **pkgname)
{
    ufunc *u = currently_called_function();

    if (u != NULL) {
        if (funcname != NULL) {
            *funcname = u->name;
        }
        if (pkgname != NULL && u->pkg != NULL) {
            *pkgname = u->pkg->name;
        }
    }
}

int function_is_executing (const char *funcname)
{
    int ret = 0;

    if (callstack != NULL) {
        GList *tmp = g_list_last(callstack);
        fncall *call;

        while (tmp != NULL) {
            call = tmp->data;
            if (!strcmp(call->fun->name, funcname)) {
                ret = 1;
                break;
            }
            tmp = tmp->prev;
        }
    }

    return ret;
}

fnpkg *get_active_function_package (gretlopt opt)
{
    ufunc *u = currently_called_function();

    if (u != NULL && u->pkg != NULL) {
        if (opt == OPT_NONE) {
            return u->pkg;
        } else if ((opt & OPT_O) && u->pkg->overrides) {
            /* in this case we're only interested if the
               package overrides any built-in functions
            */
            return u->pkg;
        }
    }

    return NULL;
}

const char *function_package_translate (const char *id)
{
    fnpkg *pkg = get_active_function_package(OPT_NONE);

    if (pkg != NULL && pkg->trans != NULL) {
        return get_gfn_translation(pkg->trans, id);
    } else {
        return id;
    }
}

void *function_package_translation (fnpkg *pkg)
{
    return pkg == NULL ? NULL : pkg->trans;
}

fnpkg *gretl_function_get_package (const ufunc *fun)
{
    return fun == NULL ? NULL : fun->pkg;
}

#if CALL_DEBUG

static void print_callstack (ufunc *fun)
{
    GList *tmp = callstack;
    fncall *call;
    char s[32];

    if (tmp == NULL) {
        fputs("callstack: empty\n", stderr);
        return;
    }

    fputs("callstack:\n", stderr);
    while (tmp != NULL) {
        call = tmp->data;
        *s = '\0';
        if (fun == call->fun) {
            strcpy(s, " *");
        }
        if (call->flags & FC_PRESERVE) {
            strcat(s, " preserve");
        }
        if (call->flags & FC_RECURSING) {
            strcat(s, " recursing");
        }
        fprintf(stderr, "  %s %p%s\n", call->fun->name, (void *) call, s);
        tmp = tmp->next;
    }
}

#endif

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

int gretl_function_recursing (void)
{
    if (callstack != NULL) {
        GList *tmp = g_list_last(callstack);
        fncall *call = tmp->data;

        return is_recursing(call);
    } else {
        return 0;
    }
}

#ifdef HAVE_MPI

static fnpkg *find_caller_package (const char *name)
{
    fnpkg *pkg = NULL;
    int i;

    for (i=0; i<n_ufuns; i++) {
        if (!strcmp(name, ufuns[i]->name)) {
            if (ufuns[i]->pkg != NULL) {
                ufuns[i]->pkg->prechecked = 1;
                pkg = ufuns[i]->pkg;
            }
            break;
        }
    }

    return pkg;
}

#endif

/**
 * get_user_function_by_name:
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
    fnpkg *pkg = current_pkg;
    ufunc *fun = NULL;
    int i;

    if (n_ufuns == 0) {
        return NULL;
    }

    if (pkg == NULL) {
        fun = currently_called_function();
        if (fun != NULL) {
            pkg = fun->pkg;
            fun = NULL;
        }
    }

#ifdef HAVE_MPI
    if (pkg == NULL && *mpi_caller != '\0') {
        pkg = find_caller_package(mpi_caller);
    }
#endif

    if (pkg != NULL) {
        /* There's an active function package: try first
           for functions that belong to the package.
        */
        for (i=0; i<pkg->n_pub; i++) {
            /* public members */
            if (!strcmp(name, pkg->pub[i]->name)) {
                fun = pkg->pub[i];
                break;
            }
        }
        if (fun == NULL) {
            /* private members */
            for (i=0; i<pkg->n_priv; i++) {
                if (!strcmp(name, pkg->priv[i]->name)) {
                    fun = pkg->priv[i];
                    break;
                }
            }
        }
        if (fun == NULL && pkg->provider != NULL) {
            /* functions shared by provider */
            fnpkg *prv = get_function_package_by_name(pkg->provider);

            if (prv != NULL) {
                for (i=0; i<prv->n_priv; i++) {
                    if (!strcmp(name, prv->priv[i]->name)) {
                        fun = prv->priv[i];
                        break;
                    }
                }
            }
        }
    }

    if (fun == NULL) {
        /* Match any function (WAS any non-private function) */
        for (i=0; i<n_ufuns; i++) {
            if (!function_is_private(ufuns[i]) &&
                !strcmp(name, ufuns[i]->name)) {
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
 * is_user_function:
 * @name: name to test.
 *
 * Returns: 1 if @name identifies a hansl function accessible
 * in the current context, otherwise 0.
 */

int is_user_function (const char *name)
{
    return get_user_function_by_name(name) != NULL;
}

/**
 * get_function_from_package:
 * @funname: name of function to retrieve.
 * @pkg: function package.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given @funname that is associated with
 * function package @pkg, otherwise NULL.  This is used
 * in the gretl function package editor.
 */

ufunc *get_function_from_package (const char *funname, fnpkg *pkg)
{
    if (pkg != NULL) {
        int i;

        for (i=0; i<pkg->n_pub; i++) {
            if (!strcmp(funname, pkg->pub[i]->name)) {
                return pkg->pub[i];
            }
        }
        for (i=0; i<pkg->n_priv; i++) {
            if (!strcmp(funname, pkg->priv[i]->name)) {
                return pkg->priv[i];
            }
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
        params[i].deflt = UNSET_VALUE;
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
    fun->call = NULL;
    fun->calls = NULL;
    fun->flags = 0;
    fun->pkg_role = 0;

    fun->n_lines = 0;
    fun->line_idx = 1;
    fun->lines = NULL;

    fun->n_params = 0;
    fun->params = NULL;

    fun->rettype = GRETL_TYPE_NONE;

    return fun;
}

static void free_lines_array (fn_line *lines, int n)
{
    fn_line *line;
    int i;

    if (lines == NULL) {
        return;
    }

#if COMP_DEBUG
    fprintf(stderr, "free_lines_array (n=%d)\n", n);
#endif

    for (i=0; i<n; i++) {
        line = &lines[i];
        free(line->s);
        if (line_has_loop(line)) {
#if COMP_DEBUG
            fprintf(stderr, "  L%d: has_loop at %p\n", i, (void *) line->ptr);
#endif
            gretl_loop_destroy(line->ptr);
        }
    }

    free(lines);
}

static void free_params_array (fn_param *params, int n)
{
    int i;

    if (params == NULL) return;

    for (i=0; i<n; i++) {
        free(params[i].name);
        free(params[i].descrip);
        strings_array_free(params[i].labels, params[i].nlabels);
    }
    free(params);
}

static void destroy_ufunc_calls (ufunc *fun)
{
    if (fun->calls != NULL) {
        g_list_free_full(fun->calls, fncall_destroy);
        fun->calls = NULL;
        fun->call = NULL;
    } else if (fun->call != NULL) {
        fncall_destroy(fun->call);
        fun->call = NULL;
    }
}

static void clear_ufunc_data (ufunc *fun)
{
    free_lines_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);

    fun->lines = NULL;
    fun->params = NULL;

    fun->n_lines = 0;
    fun->line_idx = 1;
    fun->n_params = 0;

    destroy_ufunc_calls(fun);

    fun->rettype = GRETL_TYPE_NONE;
}

static void ufunc_free (ufunc *fun)
{
    /* free_lines_array() handles any attached loops */
    free_lines_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    /* destroy_ufunc_calls() handles attached genrs */
    destroy_ufunc_calls(fun);
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

static int no_scalar_default (fn_param *fp)
{
    int ret = 0;

    if (default_unset(fp)) {
        ret = 1;
    } else if (fp->type != GRETL_TYPE_DOUBLE && na(fp->deflt)) {
        ret = 1;
    }

    return ret;
}

/* handling of XML function packages */

enum {
    FUNCS_INFO,
    FUNCS_LOAD,
    FUNCS_CODE,
    FUNCS_SAMPLE,
    FUNCS_HELP,
    FUNCS_QUERY
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
    } else if (t == GRETL_TYPE_STRING_REF) {
        return "stringref";
    } else if (t == GRETL_TYPE_STRINGS_REF) {
        return "stringsref";
    } else if (t == GRETL_TYPE_MATRICES_REF) {
        return "matricesref";
    } else if (t == GRETL_TYPE_BUNDLES_REF) {
        return "bundlesref";
    } else if (t == GRETL_TYPE_LISTS_REF) {
        return "listsref"; /* not actually allowed */
    } else {
        return gretl_type_get_name(t);
    }
}

static GretlType return_type_from_string (const char *s,
                                          int *err)
{
    GretlType t;

    if (!strcmp(s, "void")) {
        /* not OK as arg type, but OK as return */
        t = GRETL_TYPE_VOID;
    } else {
        t = gretl_type_from_string(s);
    }

    if (!ok_function_return_type(t)) {
        if (*s == '\0') {
            gretl_errmsg_sprintf(_("Missing function return type"));
        } else if (t == GRETL_TYPE_NONE) {
            gretl_errmsg_sprintf(_("Expected a function return type, found '%s'"),
                                 s);
        } else {
            gretl_errmsg_sprintf(_("%s: invalid return type for function"),
                                 s);
        }
        *err = E_TYPES;
    }

    return t;
}

static GretlType param_field_to_type (const char *s,
                                      const char *funname,
                                      int *err)
{
    GretlType t = gretl_type_from_string(s);

    if (!ok_function_arg_type(t)) {
        gretl_errmsg_sprintf(_("function %s: invalid parameter type '%s'"),
                             funname, s);
        *err = E_INVARG;
    }

    return t;
}

static gboolean special_int_default (fn_param *param,
                                     xmlNodePtr np)
{
    gboolean ret = FALSE;

    if (param->type == GRETL_TYPE_INT) {
        char *s = NULL;

        gretl_xml_get_prop_as_string(np, "default", &s);
        if (s != NULL) {
            if (strstr(s, "$mylist")) {
                param->deflt = INT_USE_MYLIST;
                ret = TRUE;
            } else if (strstr(s, "$xlist")) {
                param->deflt = INT_USE_XLIST;
                ret = TRUE;
            }
            free(s);
        }
    }

    return ret;
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
                param->type = param_field_to_type(field, fun->name, &err);
                free(field);
                if (special_int_default(param, cur)) {
                    ; /* handled */
                } else if (gretl_scalar_type(param->type)) {
                    double x;

                    if (gretl_xml_get_prop_as_double(cur, "default", &x)) {
                        param->deflt = x;
                    } else {
                        param->deflt = UNSET_VALUE;
                    }
                    if (param->type != GRETL_TYPE_BOOL) {
                        gretl_xml_get_prop_as_double(cur, "min", &param->min);
                        gretl_xml_get_prop_as_double(cur, "max", &param->max);
                        gretl_xml_get_prop_as_double(cur, "step", &param->step);
                    }
                }
                if (gretl_xml_get_prop_as_bool(cur, "auto")) {
                    set_param_auto(param);
                } else if (gretl_xml_get_prop_as_bool(cur, "optional")) {
                    set_param_optional(param);
                }
                if (gretl_xml_get_prop_as_bool(cur, "const")) {
                    maybe_set_param_const(param);
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

/* When we come here @fun may already have an array of lines allocated
   (if we know its size in advance). In that case @n_alloced will be
   non-NULL and will contain the number of lines allocated.  Otherwise
   the array of lines will be NULL initially, and @n_alloced will be
   NULL.

   In both cases fun->n_lines will be zero on the first call to this
   function, and it must be incremented each time a line is saved.
*/

static int push_function_line (ufunc *fun, char *s, int ci,
                               int donate, int *n_alloced)
{
    fn_line *ln = NULL;
    int n = fun->n_lines + 1;
    int err = 0;

    if (n_alloced == NULL || n > *n_alloced) {
        /* create or enlarge the array of lines */
        fn_line *lines = realloc(fun->lines, n * sizeof *lines);

        if (lines == NULL) {
            return E_ALLOC;
        }
        fun->lines = lines;
        fun->n_lines = n;
        if (n_alloced != NULL) {
            *n_alloced = n;
        }
    } else {
        /* we already have sufficient storage */
        fun->n_lines += 1;
    }

    ln = &fun->lines[n-1];
    ln->idx = fun->line_idx;

    if (donate) {
        /* take ownership of @s */
        ln->s = s;
    } else {
        /* make a copy of @s */
        ln->s = gretl_strdup(s);
        if (ln->s == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        ln->ptr = NULL;
        ln->ci = ci;
        ln->next = 0;
        /* be on the safe side wrt string substitution */
        ln->flags = strchr(ln->s, '$') ? LINE_NOCOMP : 0;
        if (ci == CMD_COMMENT) {
            ln->flags |= LINE_IGNORE;
        }
        fun->line_idx += 1;
    }

    return err;
}

/* Get a quick reading on the number of non-blank lines in @buf, so we
   can pre-allocate all the lines needed to store a function
   definition. This may be an overestimate if some "blank" lines
   actually contain spaces.
*/

static int buf_get_n_lines (const char *buf)
{
    const char *p = buf;
    const char *prev = NULL;
    int n = 1;

    while ((p = strchr(p, '\n')) != NULL) {
        if (prev != NULL && p - prev > 1) {
            n++;
        }
        prev = p++;
    }

    return n;
}

/* does a line that starts with "set" actually change
   the program state? heuristic
*/

static int changes_state (const char *line)
{
    if (string_is_blank(line)) {
	return 0;
    } else {
	char setarg[16];

	sscanf(line, "%15s", setarg);
	return strcmp(setarg, "stopwatch");
    }
}

/* Read the lines of hansl code from the <code> element of the XML
   representation of a function.
*/

static int func_read_code (xmlNodePtr node, xmlDocPtr doc, ufunc *fun)
{
    ExecState state = {0};
    CMD cmd = {0};
    char line[MAXLINE];
    char *buf, *s;
    gint8 uses_set = 0;
    gint8 has_flow = 0;
    int save_comments = 1; /* added 2023-03-16 */
    int n_lines = 0;
    int err = 0;

    buf = (char *) xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
    if (buf == NULL) {
        return 1;
    }

    n_lines = buf_get_n_lines(buf);
    bufgets_init(buf);
    gretl_cmd_init(&cmd);
    state.cmd = &cmd;

    /* pre-allocate storage */
    fun->lines = calloc(n_lines, sizeof *fun->lines);

    while (bufgets(line, sizeof line, buf)) {
        s = line;
        while (isspace(*s)) s++;
        state.line = s;
        err = get_command_index(&state, FUNC, save_comments);
        if (err) {
            break;
        }
        if (cmd.ci == SET && changes_state(s + 3)) {
            uses_set = 1;
        } else if (cmd.ci == IF || cmd.ci == LOOP) {
            has_flow = 1;
        }
        tailstrip(s);
        if (string_is_blank(s)) {
            fun->line_idx += 1;
        } else {
            err = push_function_line(fun, s, cmd.ci, 0, &n_lines);
        }
    }

    if (fun->n_lines < n_lines) {
        /* shrink to avoid a memory leak later */
        fun->lines = realloc(fun->lines, n_lines * sizeof *fun->lines);
    }

    if (uses_set) {
        fun->flags |= UFUN_USES_SET;
    }
    if (has_flow) {
        fun->flags |= UFUN_HAS_FLOW;
    }

    bufgets_finalize(buf);
    free(buf);

    if (!err && has_flow) {
        err = ufunc_get_structure(fun);
    }

    gretl_cmd_free(&cmd);

    return err;
}

static void print_opt_flags (fn_param *param, PRN *prn)
{
    if (param_is_auto(param)) {
        pputs(prn, "[auto]");
    } else if (param_is_optional(param)) {
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
    if (na(param->min) && na(param->max) && default_unset(param)) {
        return; /* no-op */
    } else if (na(param->min) && na(param->max)) {
        /* got a default value only? */
        if (!default_unset(param)) {
            if (na(param->deflt)) {
                pputs(prn, "[NA]");
            } else if (param->deflt == INT_USE_XLIST) {
                pputs(prn, "[$xlist]");
            } else if (param->deflt == INT_USE_MYLIST) {
                pputs(prn, "[$mylist]");
            } else {
                pprintf(prn, "[%g]", param->deflt);
            }
        }
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
    if (!default_unset(param)) {
        if (na(param->deflt)) {
            pputs(prn, "NA");
        } else {
            pprintf(prn, "%g", param->deflt);
        }
    }

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
        return;
    }

#if PKG_DEBUG
    fprintf(stderr, "ufunc_unload: name %s, pkg %p\n",
            fun->name, (void *) fun->pkg);
#endif

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

    if (fun->pkg != NULL && fun->pkg->overrides) {
        delete_function_override(fun->name, fun->pkg->name);
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
    fprintf(stderr, "function_package_unload: unloading '%s', full=%d\n",
            pkg->name, full);
#endif

    /* unload children? */
    if (full) {
        for (i=0; i<pkg->n_priv; i++) {
            ufunc_unload(pkg->priv[i]);
        }
        for (i=0; i<pkg->n_pub; i++) {
            ufunc_unload(pkg->pub[i]);
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

/* Append a pointer to @fun to the array of child-pointers in @pkg: we
   do this when reading the definition of a packaged function from an
   XML file.  Note that this action does not add the functions to the
   array of loaded functions -- that's done separately, if we're
   loading the package 'for real'.
*/

static int attach_ufunc_to_package (ufunc *fun, fnpkg *pkg)
{
    ufunc **ufs;
    int n, err = 0;

    if (function_is_private(fun)) {
        n = pkg->n_priv;
        ufs = realloc(pkg->priv, (n + 1) * sizeof *ufs);
        if (ufs == NULL) {
            err = E_ALLOC;
        } else {
            pkg->priv = ufs;
            pkg->priv[n] = fun;
            pkg->n_priv += 1;
        }
    } else {
        n = pkg->n_pub;
        ufs = realloc(pkg->pub, (n + 1) * sizeof *ufs);
        if (ufs == NULL) {
            err = E_ALLOC;
        } else {
            pkg->pub = ufs;
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

static void maybe_set_menu_only (ufunc *fun, fnpkg *pkg)
{
    if (pkg->mpath != NULL && strstr(pkg->mpath, "MODELWIN")) {
        if (fun->pkg_role == UFUN_GUI_MAIN) {
            fun->flags |= UFUN_MENU_ONLY;
        }
    }
}

/* read a single user-function definition from XML file: if the
   function is a child of a package, the @pkg argument will
   be non-NULL
*/

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
        ufunc_free(fun);
        return E_DATA;
    }

    strncat(fun->name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (pkg != NULL) {
        fun->pkg = pkg;
    }

    if (gretl_xml_get_prop_as_string(node, "type", &tmp)) {
        fun->rettype = return_type_from_string(tmp, &err);
        free(tmp);
    } else {
        fun->rettype = GRETL_TYPE_VOID;
    }

    if (err) {
        ufunc_free(fun);
        return err;
    }

    if (gretl_xml_get_prop_as_bool(node, "private")) {
        fun->flags |= UFUN_PRIVATE;
    }
    if (gretl_xml_get_prop_as_bool(node, "no-print")) {
        fun->flags |= UFUN_NOPRINT;
    }
    if (gretl_xml_get_prop_as_bool(node, "menu-only")) {
        fun->flags |= UFUN_MENU_ONLY;
    }
    if (gretl_xml_get_prop_as_string(node, "pkg-role", &tmp)) {
        fun->pkg_role = pkg_key_get_role(tmp);
        free(tmp);
    }
    if (!(fun->flags & UFUN_MENU_ONLY) && pkg != NULL) {
        maybe_set_menu_only(fun, pkg);
    }

#if PKG_DEBUG
    fprintf(stderr, "read_ufunc_from_xml: name '%s', type %d\n"
            " private = %d\n", fun->name, fun->rettype,
            function_is_private(fun));
#endif

    if (pkg == NULL && function_is_private(fun)) {
        fprintf(stderr, "unpackaged function: can't be private\n");
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
            pkg->help = gretl_xml_get_string(cur, doc);
        } else if (!xmlStrcmp(cur->name, (XUC) "params")) {
            err = func_read_params(cur, doc, fun);
            if (err) {
                fprintf(stderr, "%s: error parsing function parameters\n",
                        fun->name);
            }
        } else if (!xmlStrcmp(cur->name, (XUC) "return")) {
            gretl_errmsg_set(_("Old-style function definitions no longer supported"));
            err = E_DATA;
        } else if (!xmlStrcmp(cur->name, (XUC) "code")) {
            err = func_read_code(cur, doc, fun);
        }
        cur = cur->next;
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

static void print_param_limit (fn_param *param, int i, PRN *prn)
{
    if (param->type == GRETL_TYPE_INT) {
        if (i == 0) {
            pprintf(prn, " min=\"%.0f\"", param->min);
        } else if (i == 1) {
            pprintf(prn, " max=\"%.0f\"", param->max);
        } else {
            pprintf(prn, " default=\"%.0f\"", param->deflt);
        }
    } else {
        if (i == 0) {
            pprintf(prn, " min=\"%.15g\"", param->min);
        } else if (i == 1) {
            pprintf(prn, " max=\"%.15g\"", param->max);
        } else {
            pprintf(prn, " default=\"%.15g\"", param->deflt);
        }
    }
}

#define parm_has_children(p) (p->descrip != NULL || p->nlabels > 0)

/* write out a single user-defined function as XML, according to
   gretlfunc.dtd */

static int write_function_xml (ufunc *fun, PRN *prn, int mpi)
{

    int rtype = fun->rettype;
    int this_indent = 0;
    int next_indent = 0;
    int i;

    if (rtype == GRETL_TYPE_NONE) {
        rtype = GRETL_TYPE_VOID;
    }

    pprintf(prn, "<gretl-function name=\"%s\" type=\"%s\"",
            fun->name, gretl_type_get_name(rtype));

    if (function_is_private(fun)) {
        pputs(prn, " private=\"1\"");
    }
    if (function_is_noprint(fun)) {
        pputs(prn, " no-print=\"1\"");
    }
    if (function_is_menu_only(fun)) {
        pputs(prn, " menu-only=\"1\"");
    }

    if (fun->pkg_role) {
        pprintf(prn, " pkg-role=\"%s\"", package_role_get_key(fun->pkg_role));
    }

    pputs(prn, ">\n");

    if (fun->n_params > 0) {
        fn_param *param;

        gretl_push_c_numeric_locale();
        pprintf(prn, " <params count=\"%d\">\n", fun->n_params);
        for (i=0; i<fun->n_params; i++) {
            param = &fun->params[i];
            pprintf(prn, "  <param name=\"%s\" type=\"%s\"",
                    param->name, arg_type_xml_string(param->type));
            if (!na(param->min)) {
                print_param_limit(param, 0, prn);
            }
            if (!na(param->max)) {
                print_param_limit(param, 1, prn);
            }
            if (!default_unset(param)) {
                if (na(param->deflt)) {
                    pputs(prn, " default=\"NA\"");
                } else if (param->deflt == INT_USE_MYLIST) {
                    pputs(prn, " default=\"$mylist\"");
                } else if (param->deflt == INT_USE_XLIST) {
                    pputs(prn, " default=\"$xlist\"");
                } else {
                    print_param_limit(param, 2, prn);
                }
            }
            if (!na(param->step)) {
                pprintf(prn, " step=\"%g\"", param->step);
            }
            if (param_is_auto(param)) {
                pputs(prn, " auto=\"true\"");
            } else if (param_is_optional(param)) {
                pputs(prn, " optional=\"true\"");
            }
            if (param_is_const(param)) {
                pputs(prn, " const=\"true\"");
            }
            if (parm_has_children(param)) {
                pputs(prn, ">\n"); /* terminate opening tag */
                if (param->descrip != NULL) {
                    gretl_xml_put_tagged_string("description",
                                                param->descrip,
                                                prn);
                }
                if (param->nlabels > 0) {
                    gretl_xml_put_strings_array_quoted("labels",
                                                       (const char **) param->labels,
                                                       param->nlabels, prn);
                }
                pputs(prn, "  </param>\n");
            } else {
                pputs(prn, "/>\n"); /* terminate opening tag */
            }
        }
        pputs(prn, " </params>\n");
        gretl_pop_c_numeric_locale();
    }

    pputs(prn, "<code>");

    if (mpi) {
        /* minimal output of function code */
        for (i=0; i<fun->n_lines; i++) {
            gretl_xml_put_string(fun->lines[i].s, prn);
            pputc(prn, '\n');
        }
    } else {
        /* pay some attention to readability */
        for (i=0; i<fun->n_lines; i++) {
            if (i > 0 && fun->lines[i].idx - fun->lines[i-1].idx > 1) {
                /* reinstate single blank lines */
                pputc(prn, '\n');
            }
            /* basic indentation adjustment */
            adjust_indent(fun->lines[i].s, &this_indent, &next_indent);
            bufspace(2 * this_indent, prn);
            gretl_xml_put_string(fun->lines[i].s, prn);
            pputc(prn, '\n');
        }
    }

    pputs(prn, "</code>\n");
    pputs(prn, "</gretl-function>\n");

    return 0;
}

/* script-style output */

static void print_function_start (ufunc *fun, PRN *prn)
{
    const char *s;
    int i, pos = 0;

    if (fun->rettype == GRETL_TYPE_NONE) {
        pos += pprintf(prn, "function void %s ", fun->name);
    } else {
        const char *typestr = gretl_type_get_name(fun->rettype);

        pos += pprintf(prn, "function %s %s ", typestr, fun->name);
    }

    gretl_push_c_numeric_locale();

    if (fun->n_params == 0) {
        pputs(prn, "(void)");
    } else {
        pos += pputc(prn, '(');
    }

    for (i=0; i<fun->n_params; i++) {
        fn_param *fp = &fun->params[i];

        if (fp->flags & FP_CONST) {
            pputs(prn, "const ");
        }
        s = gretl_type_get_name(fp->type);
        if (s[strlen(s) - 1] == '*') {
            pprintf(prn, "%s%s", s, fp->name);
        } else {
            pprintf(prn, "%s %s", s, fp->name);
        }
        if (fp->type == GRETL_TYPE_BOOL) {
            if (!default_unset(fp) && !na(fp->deflt)) {
                pprintf(prn, "[%g]", fp->deflt); /* FIXME? */
            }
        } else if (gretl_scalar_type(fp->type)) {
            print_min_max_deflt(fp, prn);
        } else if (arg_may_be_optional(fp->type)) {
            print_opt_flags(fp, prn);
        }
        print_param_description(fp, prn);
        if (fp->nlabels > 0) {
            print_param_labels(fp, prn);
        }
        if (i == fun->n_params - 1) {
            pputc(prn, ')');
        } else {
            pputs(prn, ",\n");
            bufspace(pos, prn);
        }
    }

    pputc(prn, '\n');

    gretl_pop_c_numeric_locale();
}

/* The following is an attempt at generic indentation based
   on occurrences of '{' and '}', which should probably
   work for R code.
*/

static void foreign_adjust_indent (const char *s,
                                   int *this_indent,
                                   int *next_indent)
{
    int n = strlen(s);

    *this_indent = *next_indent;

    if (strncmp(s, "} else", 6) == 0) {
        *next_indent = *this_indent;
        *this_indent = *this_indent - 1;
    } else if (s[n-1] == '{') {
        *next_indent = *this_indent + 1;
    } else if (strcmp(s, "}") == 0 ||
               strncmp(s, "} else", 6) == 0) {
        *this_indent = *next_indent = *this_indent - 1;
    }
}

static void maybe_toggle_foreign (const char *s, int *in_foreign)
{
    if (!*in_foreign) {
        if (strncmp(s, "foreign ", 8) == 0) {
            *in_foreign = 1;
        }
    } else if (strncmp(s, "end foreign", 11) == 0) {
        *in_foreign = 0;
    }
}

/**
 * gretl_function_print_code:
 * @u: pointer to user-function.
 * @tabwidth: number of spaces per "tab" (logical indent).
 * @prn: printing struct.
 *
 * Prints out function @fun to @prn, script-style. Note
 * that this function may be called for several functions
 * in turn, with a single @prn as target.
 *
 * Returns: 0 on success, non-zero if @fun is %NULL.
 */

int gretl_function_print_code (ufunc *u, int tabwidth, PRN *prn)
{
    PRN *ptmp;
    const char *buf;
    int this_indent = 0;
    int next_indent = 0;
    int in_foreign = 0;
    int i, j;

    if (u == NULL) {
        return E_DATA;
    }

    if (tabwidth == 0) {
        tabwidth = 2;
    }

    ptmp = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    print_function_start(u, ptmp);

    for (i=0; i<u->n_lines; i++) {
        maybe_toggle_foreign(u->lines[i].s, &in_foreign);
        if (in_foreign) {
            foreign_adjust_indent(u->lines[i].s, &this_indent, &next_indent);
        } else {
            adjust_indent(u->lines[i].s, &this_indent, &next_indent);
        }
        for (j=0; j<=this_indent; j++) {
            bufspace(tabwidth, ptmp);
        }
        pputs(ptmp, u->lines[i].s);
        if (i < u->n_lines - 1) {
            if (u->lines[i+1].idx > u->lines[i].idx + 1) {
                pputc(ptmp, '\n');
            }
        }
        pputc(ptmp, '\n');
    }

    pputs(ptmp, "end function\n");

    /* now refine the indentation */
    buf = gretl_print_get_buffer(ptmp);
    normalize_hansl(buf, tabwidth, prn);
    gretl_print_destroy(ptmp);

    return 0;
}

/* called from kalman.c for time-varying matrices */

char **gretl_function_retrieve_code (ufunc *u, int *nlines)
{
    char **S = NULL;
    fn_line *line;
    int i, j = 0;

    for (i=0; i<u->n_lines; i++) {
        line = &(u->lines[i]);
        if (!ignore_line(line)) {
            j++;
        }
    }

    if (j > 0) {
        S = strings_array_new(j);
    }

    if (S != NULL) {
        *nlines = j;
        j = 0;
        for (i=0; i<u->n_lines; i++) {
            line = &(u->lines[i]);
            if (!ignore_line(line)) {
                S[j++] = line->s;
            }
        }
    }

    return S;
}

/* construct a name for @pkg based on its filename member:
   take the basename and knock off ".gfn"
*/

static void name_package_from_filename (fnpkg *pkg)
{
    char *p = strrslash(pkg->fname);
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
        if (strstr(fun->lines[i].s, "## no-print ##")) {
            fun->flags |= UFUN_NOPRINT;
        } else if (strstr(fun->lines[i].s, "## menu-only ##")) {
            fun->flags |= UFUN_MENU_ONLY;
        }
    }
}

/* before revising the function members of an existing
   package, detach all its current functions
*/

static void package_disconnect_funcs (fnpkg *pkg)
{
    int i;

    if (pkg->pub != NULL) {
        for (i=0; i<pkg->n_pub; i++) {
            pkg->pub[i]->pkg = NULL;
        }
        free(pkg->pub);
        pkg->pub = NULL;
        pkg->n_pub = 0;
    }

    if (pkg->priv != NULL) {
        for (i=0; i<pkg->n_priv; i++) {
            pkg->priv[i]->pkg = NULL;
            set_function_private(pkg->priv[i], FALSE);
        }
        free(pkg->priv);
        pkg->priv = NULL;
        pkg->n_priv = 0;
    }
}

/* Given an array of @n function names in @names, set up a
   corresponding array of pointers to the named functions in @pkg,
   Flag an error if any of the function names are bad, or if
   allocation fails.
*/

static int set_uf_array_from_names (fnpkg *pkg, char **names,
                                    int n, int priv)
{
    ufunc **uf = NULL;
    ufunc *fun;
    int i, j;

    /* check the supplied names */
    for (i=0; i<n; i++) {
        for (j=0; j<i; j++) {
            if (!strcmp(names[j], names[i])) {
                gretl_errmsg_sprintf(_("Duplicated function name '%s'"), names[i]);
                return E_DATA;
            }
        }
        if (get_uf_array_member(names[i], NULL) == NULL) {
            fprintf(stderr, "%s: function not found!\n", names[i]);
            return E_DATA;
        }
    }

    /* allocate storage */
    if (n > 0) {
        uf = malloc(n * sizeof *uf);
        if (uf == NULL) {
            return E_ALLOC;
        }
    }

    /* connect the specified functions */
    for (i=0; i<n; i++) {
        fun = get_uf_array_member(names[i], NULL);
        fun->pkg = pkg;
        set_function_private(fun, priv);
        check_special_comments(fun);
        uf[i] = fun;
    }

    /* set up the package info */
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
    int err;

    /* clear the decks first */
    package_disconnect_funcs(pkg);

    err = set_uf_array_from_names(pkg, pubnames, n_pub, 0);

    if (!err) {
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
    fnpkg *pkg = NULL;

    if (n_pub <= 0) {
        /* we need at least one public function */
        *err = E_DATA;
    } else {
        pkg = function_package_alloc(fname);
        if (pkg == NULL) {
            *err = E_ALLOC;
        }
    }

    if (*err) {
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

static char *trim_text (char *s)
{
    while (isspace(*s)) s++;
    return tailstrip(s);
}

static int addon_write_translatable_strings (fnpkg *pkg, PRN *prn)
{
    FILE *fp;
    gchar *trname;
    char **S = NULL;
    int trans = 0;
    int i, n = 0;

    /* open a file to contain the results */
    trname = g_strdup_printf("%s-i18n.c", pkg->name);
    fp = gretl_fopen(trname, "wb");
    if (fp == NULL) {
        gretl_errmsg_sprintf(_("Couldn't open %s"), trname);
        g_free(trname);
        return E_FOPEN;
    }

    if (pkg->pub != NULL) {
        /* pick up any parameter descriptions and/or parameter-value
           enumeration strings
        */
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
        /* we got some relevant content */
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
            strings_array_free(S, n);
        }
        fputs("};\n", fp);
        trans = 1;
    }

    fclose(fp);

    if (!trans) {
        gretl_remove(trname);
    } else {
        pprintf(prn, "Wrote translations file %s\n", trname);
    }

    g_free(trname);

    return 0;
}

/* Note: return 0 on failure; 1 on writing _something_ but no
   T_()-marked strings; 2 on writing some T_()-marked strings.
*/

int package_write_translatables (fnpkg *pkg, PRN *prn)
{
    gchar *enstr;
    char **S = NULL;
    const char *s;
    ufunc *fun;
    fn_param *param;
    int T_strs = 0;
    int i, j, k;
    int ns = 0;
    int ret = 0;

    if (pkg->pub != NULL) {
        /* pick up any parameter descriptions and/or parameter-value
           enumeration strings
        */
        for (i=0; i<pkg->n_pub; i++) {
            fun = pkg->pub[i];
            for (j=0; j<fun->n_params; j++) {
                param = &fun->params[j];
                if (param->descrip != NULL) {
                    strings_array_add_uniq(&S, &ns, param->descrip, NULL);
                }
                for (k=0; k<param->nlabels; k++) {
                    strings_array_add_uniq(&S, &ns, param->labels[k], NULL);
                }
            }
            /* also check lines of code for @fun, looking for T_() */
            for (j=0; j<fun->n_lines; j++) {
                s = fun->lines[j].s;
                while ((s = strstr(s, "T_(")) != NULL) {
                    enstr = get_translatable_content(&s);
                    if (enstr != NULL) {
                        strings_array_add_uniq(&S, &ns, enstr, NULL);
                        g_free(enstr);
                        T_strs++;
                    }
                }
            }
        }
    }

    if (pkg->priv != NULL) {
        for (i=0; i<pkg->n_priv; i++) {
            fun = pkg->priv[i];
            for (j=0; j<fun->n_lines; j++) {
                s = fun->lines[j].s;
                 while ((s = strstr(s, "T_(")) != NULL) {
                    enstr = get_translatable_content(&s);
                    if (enstr != NULL) {
                        strings_array_add_uniq(&S, &ns, enstr, NULL);
                        g_free(enstr);
                        T_strs++;
                    }
                }
            }
        }
    }

    if (pkg->label != NULL || S != NULL) {
        /* we got some relevant content */
        pputs(prn, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        pputs(prn, "<translation lang=\"xx\">\n");
        if (pkg->label != NULL) {
            pprintf(prn, " <msg en=\"%s\"></msg>\n", pkg->label);
        }
        if (S != NULL) {
            for (i=0; i<ns; i++) {
                pprintf(prn, " <msg en=\"%s\"></msg>\n", S[i]);
            }
            strings_array_free(S, ns);
        }
        pputs(prn, "</translation>\n");
        ret = 1 + (T_strs > 0);
    }

    return ret;
}

static int package_write_index (fnpkg *pkg, PRN *inprn)
{
    PRN *prn;
    gchar *idxname;
    int err = 0;

    idxname = g_strdup_printf("%s.xml", pkg->name);

    prn = gretl_print_new_with_filename(idxname, &err);
    if (prn == NULL) {
        g_free(idxname);
        return err;
    }

    gretl_xml_header(prn);

    pprintf(prn, "<gretl-addon name=\"%s\"", pkg->name);

    if (pkg->dreq == FN_NEEDS_TS) {
        pprintf(prn, " %s=\"true\"", NEEDS_TS);
    } else if (pkg->dreq == FN_NEEDS_QM) {
        pprintf(prn, " %s=\"true\"", NEEDS_QM);
    } else if (pkg->dreq == FN_NEEDS_PANEL) {
        pprintf(prn, " %s=\"true\"", NEEDS_PANEL);
    } else if (pkg->dreq == FN_NODATA_OK) {
        pprintf(prn, " %s=\"true\"", NO_DATA_OK);
    }

    if (pkg->modelreq > 0) {
        pprintf(prn, " model-requirement=\"%s\"",
                gretl_command_word(pkg->modelreq));
    }

    if (pkg->minver > 0) {
        char vstr[10];

        pprintf(prn, " minver=\"%s\"",
                gretl_version_string(vstr, pkg->minver));
    }

    if (pkg->uses_subdir) {
        pputs(prn, " lives-in-subdir=\"true\"");
    }

    if (pkg->data_access) {
        pputs(prn, " wants-data-access=\"true\"");
    }

    pputs(prn, ">\n");

    gretl_xml_put_tagged_string("author",  pkg->author, prn);
    gretl_xml_put_tagged_string("version", pkg->version, prn);
    gretl_xml_put_tagged_string("date",    pkg->date, prn);
    gretl_xml_put_tagged_string("description", pkg->descrip, prn);

    if (pkg->help_fname != NULL) {
	gretl_xml_put_tagged_string("helpfile", pkg->help_fname, prn);
    } else if (pkg->help != NULL && !strncmp(pkg->help, "pdfdoc:", 7)) {
	gretl_xml_put_tagged_string("helpfile", pkg->help + 7, prn);
    }

    if (pkg->tags != NULL) {
        gretl_xml_put_tagged_string("tags", pkg->tags, prn);
    }
    if (pkg->label != NULL) {
        gretl_xml_put_tagged_string("label", pkg->label, prn);
    }
    if (pkg->mpath != NULL) {
        gretl_xml_put_tagged_string("menu-attachment", pkg->mpath, prn);
    }

    pputs(prn, "</gretl-addon>\n");

    gretl_print_destroy(prn);

    pprintf(inprn, "Wrote index file %s\n", idxname);
    g_free(idxname);

    return 0;
}

/* Write out a function package as XML. If @fp is NULL, then we write
   this as a package file in its own right, using the filename given
   in the package (the 'standalone' case), otherwise we're writing it
   as a component of the functions listing in a gretl session file,
   which already has an associated open stream.
*/

static int real_write_function_package (fnpkg *pkg, PRN *prn, int mpi)
{
    int standalone = (prn == NULL);
    int i, err = 0;

    if (standalone) {
        prn = gretl_print_new_with_filename(pkg->fname, &err);
        if (prn == NULL) {
            return err;
        } else {
            gretl_xml_header(prn);
            pputs(prn, "<gretl-functions>\n");
        }
    }

    pputs(prn, "<gretl-function-package");

    if (pkg->name[0] == '\0') {
        name_package_from_filename(pkg);
    }

    pprintf(prn, " name=\"%s\"", pkg->name);

    if (pkg->dreq == FN_NEEDS_TS) {
        pprintf(prn, " %s=\"true\"", NEEDS_TS);
    } else if (pkg->dreq == FN_NEEDS_QM) {
        pprintf(prn, " %s=\"true\"", NEEDS_QM);
    } else if (pkg->dreq == FN_NEEDS_PANEL) {
        pprintf(prn, " %s=\"true\"", NEEDS_PANEL);
    } else if (pkg->dreq == FN_NODATA_OK) {
        pprintf(prn, " %s=\"true\"", NO_DATA_OK);
    }

    if (pkg->modelreq > 0) {
        pprintf(prn, " model-requirement=\"%s\"",
                gretl_command_word(pkg->modelreq));
    }

    if (pkg->minver > 0) {
        char vstr[10];

        pprintf(prn, " minver=\"%s\"", gretl_version_string(vstr, pkg->minver));
    }
    if (pkg->uses_subdir) {
        pprintf(prn, " lives-in-subdir=\"true\"");
    }
    if (pkg->data_access) {
        pprintf(prn, " wants-data-access=\"true\"");
    }

    pputs(prn, ">\n");

    if (pkg->email != NULL && *pkg->email != '\0') {
        gretl_xml_put_tagged_string_plus("author", pkg->author,
                                         "email", pkg->email,
                                         prn);
    } else {
        gretl_xml_put_tagged_string("author", pkg->author, prn);
    }
    gretl_xml_put_tagged_string("version", pkg->version, prn);
    gretl_xml_put_tagged_string("date",    pkg->date, prn);
    gretl_xml_put_tagged_string("description", pkg->descrip, prn);

    if (pkg->tags != NULL) {
        gretl_xml_put_tagged_string("tags", pkg->tags, prn);
    }
    if (pkg->label != NULL) {
        gretl_xml_put_tagged_string("label", pkg->label, prn);
    }

    if (pkg->mpath != NULL) {
        gretl_xml_put_tagged_string("menu-attachment", pkg->mpath, prn);
    }
    if (pkg->help != NULL && !mpi) {
        if (pkg->help_fname != NULL) {
            pprintf(prn, "<help filename=\"%s\">\n", pkg->help_fname);
        } else {
            pputs(prn, "<help>\n");
        }
        gretl_xml_put_string(trim_text(pkg->help), prn);
        pputs(prn, "\n</help>\n");
    }

    if (pkg->gui_help != NULL && !mpi) {
        if (pkg->gui_help_fname != NULL) {
            pprintf(prn, "<gui-help filename=\"%s\">\n",
                              pkg->gui_help_fname);
        } else {
            pputs(prn, "<gui-help>\n");
        }
        gretl_xml_put_string(trim_text(pkg->gui_help), prn);
        pputs(prn, "\n</gui-help>\n");
    }

    if (pkg->datafiles != NULL) {
        gretl_xml_put_strings_array("data-files",
                                    (const char **) pkg->datafiles,
                                    pkg->n_files, prn);
    }

    if (pkg->depends != NULL) {
        gretl_xml_put_strings_array("depends",
                                    (const char **) pkg->depends,
                                    pkg->n_depends, prn);
    }
    if (pkg->provider != NULL) {
        gretl_xml_put_tagged_string("provider", pkg->provider, prn);
    }
    if (pkg->Rdeps != NULL) {
        gretl_xml_put_tagged_string("R-depends", pkg->Rdeps, prn);
    }
    if (pkg->pub != NULL) {
        for (i=0; i<pkg->n_pub; i++) {
            write_function_xml(pkg->pub[i], prn, mpi);
        }
    }
    if (pkg->priv != NULL) {
        for (i=0; i<pkg->n_priv; i++) {
            write_function_xml(pkg->priv[i], prn, mpi);
        }
    }
    if (pkg->sample != NULL && !mpi) {
        if (pkg->sample_fname != NULL) {
            pprintf(prn, "<sample-script filename=\"%s\">\n",
                              pkg->sample_fname);
        } else {
            pputs(prn, "<sample-script>\n");
        }
        gretl_xml_put_string(trim_text(pkg->sample), prn);
        pputs(prn, "\n</sample-script>\n");
    }

    if (pkg->trans != NULL) {
        write_translation(pkg->trans, prn);
    }

    pputs(prn, "</gretl-function-package>\n");

    if (standalone) {
        pputs(prn, "</gretl-functions>\n");
        gretl_print_destroy(prn);
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
    return real_write_function_package(pkg, NULL, 0);
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

static void strip_cr (gchar *s)
{
    while (*s) {
        if (*s == 0x0d && *(s+1) == 0x0a) {
            /* CR + LF -> LF */
            memmove(s, s+1, strlen(s));
            s++;
        }
        s++;
    }
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
    } else if (strchr(ret, '\r')) {
        strip_cr(ret);
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

static int pkg_set_modelreq (fnpkg *pkg, const char *s)
{
    int ci = gretl_command_number(s);

    if (ci > 0) {
        pkg->modelreq = ci;
        return 0;
    } else {
        return E_PARSE;
    }
}

static int pkg_set_datafiles (fnpkg *pkg, const char *s)
{
    int err = 0;

    if (string_is_blank(s)) {
        err = E_DATA;
    } else {
        int n = 0;

        pkg->datafiles = gretl_string_split(s, &n, NULL);
        if (pkg->datafiles == NULL) {
            pkg->n_files = 0;
            err = E_ALLOC;
        } else {
            pkg->n_files = n;
        }
    }

    if (!err) {
        pkg->uses_subdir = 1;
    }

    return err;
}

static int pkg_set_depends (fnpkg *pkg, const char *s)
{
    int err = 0;

    if (string_is_blank(s)) {
        err = E_DATA;
    } else {
        int n = 0;

        pkg->depends = gretl_string_split(s, &n, NULL);
        if (pkg->depends == NULL) {
            pkg->n_depends = 0;
            err = E_ALLOC;
        } else {
            pkg->n_depends = n;
        }
    }

    return err;
}

static int pkg_boolean_from_string (const char *s)
{
    if (!strcmp(s, "true")) {
        return 1;
    } else {
        return 0;
    }
}

static int pkg_remove_role (fnpkg *pkg, int role)
{
    ufunc *u;
    int i;

    for (i=0; i<pkg->n_priv; i++) {
        u = pkg->priv[i];
        if (u->pkg_role == role) {
            u->pkg_role = UFUN_ROLE_NONE;
            return 0;
        }
    }
    for (i=0; i<pkg->n_pub; i++) {
        u = pkg->pub[i];
        if (u->pkg_role == role) {
            u->pkg_role = UFUN_ROLE_NONE;
            return 0;
        }
    }

    return E_DATA;
}

static int valid_list_maker (ufunc *u)
{
    if (u->n_params == 0) {
        return 1; /* OK */
    } else if (u->n_params == 1 &&
               u->params[0].type == GRETL_TYPE_LIST) {
        return 1; /* OK */
    } else {
        return 0;
    }
}

static int valid_plot_checker (ufunc *u)
{
    return (u->n_params == 1 &&
	    u->params[0].type == GRETL_TYPE_BUNDLE_REF);
}

#define PRIVATE_ONLY(r) (r == UFUN_GUI_PRECHECK ||	\
			 r == UFUN_PLOT_PRECHECK ||	\
                         r == UFUN_LIST_MAKER ||	\
                         r == UFUN_R_SETUP ||		\
                         r == UFUN_UI_MAKER)

static ufunc *get_package_member_by_name (const char *name,
                                          fnpkg *pkg,
                                          int public_ok,
                                          int private_ok)
{
    int i;

    if (public_ok) {
        for (i=0; i<pkg->n_pub; i++) {
            if (!strcmp(name, pkg->pub[i]->name)) {
                return pkg->pub[i];
            }
        }
    }
    if (private_ok) {
        for (i=0; i<pkg->n_priv; i++) {
            if (!strcmp(name, pkg->priv[i]->name)) {
                return pkg->priv[i];
            }
        }
    }

    return NULL;
}

int function_set_package_role (const char *name, fnpkg *pkg,
                               const char *attr, PRN *prn)
{
    ufunc *u = NULL;
    int role = pkg_key_get_role(attr);
    int err = 0;

    if (name == NULL) {
        /* removing a role altogether */
        pkg_remove_role(pkg, role);
        return 0;
    } else if (role == UFUN_ROLE_NONE) {
        /* canceling a function's role */
        u = get_package_member_by_name(name, pkg, 1, 1);
        if (u != NULL) {
            u->pkg_role = role;
            return 0;
        } else {
            return E_DATA;
        }
    }

    /* Check that the function identified by @name satisfies the
       requirements for the role identified by @attr: if so, hook
       the function up in this role.
    */

    if (PRIVATE_ONLY(role)) {
        /* function roles that must be private */
        u = get_package_member_by_name(name, pkg, 0, 1);
        if (u == NULL) {
            pprintf(prn, "%s: %s: no such private function\n", attr, name);
            return E_DATA;
        }
        if (role == UFUN_GUI_PRECHECK) {
            if (u->rettype != GRETL_TYPE_DOUBLE) {
                pprintf(prn, "%s: must return a scalar\n", attr);
                err = E_TYPES;
            } else if (u->n_params > 0) {
                pprintf(prn, "%s: no parameters are allowed\n", attr);
                err = E_TYPES;
            }
	} else if (role == UFUN_PLOT_PRECHECK) {
            if (u->rettype != GRETL_TYPE_MATRIX) {
                pprintf(prn, "%s: must return a matrix (vector)\n", attr);
                err = E_TYPES;
            } else if (!valid_plot_checker(u)) {
                pprintf(prn, "%s: must have a single bundle parameter\n", attr);
                err = E_TYPES;
            }
        } else if (role == UFUN_LIST_MAKER) {
            if (u->rettype != GRETL_TYPE_LIST) {
                pprintf(prn, "%s: must return a list\n", attr);
                err = E_TYPES;
            } else if (!valid_list_maker(u)) {
                pprintf(prn, "%s: must have 0 parameters or a single "
                        "list parameter\n", attr);
                err = E_TYPES;
            }
        } else if (role == UFUN_R_SETUP) {
            if (u->rettype != GRETL_TYPE_VOID) {
                pprintf(prn, "%s: should not return anything\n", attr);
                err = E_TYPES;
            } else if (u->n_params > 0) {
                pprintf(prn, "%s: no parameters are allowed\n", attr);
                err = E_TYPES;
            }
        } else if (role == UFUN_UI_MAKER) {
            if (u->rettype != GRETL_TYPE_BUNDLE) {
                pprintf(prn, "%s: should return a bundle\n", attr);
                err = E_TYPES;
            } else if (u->n_params > 0) {
                pprintf(prn, "%s: no parameters are allowed\n", attr);
                err = E_TYPES;
            }
        }
        if (!err) {
            u->pkg_role = role;
        }
        return err; /* private cases handled */
    }

    /* all other special-role functions may be public */
    u = get_package_member_by_name(name, pkg, 1, 1);
    if (u == NULL) {
        pprintf(prn, "%s: %s: no such function\n", attr, name);
        return E_DATA;
    }

    if (role == UFUN_GUI_MAIN) {
        ; /* OK, type doesn't matter */
    } else {
        /* bundle-print, bundle-plot, etc: these cases require
           some further checks
        */
        int fcast = (role == UFUN_BUNDLE_FCAST);
        int j, pmin = fcast ? 2 : 1;
        GretlType pjt;

        if (u->n_params == 0) {
            pprintf(prn, "%s: must take a %s argument\n", attr,
                    gretl_type_get_name(GRETL_TYPE_BUNDLE_REF));
            err = E_TYPES;
        }
        for (j=0; j<u->n_params && !err; j++) {
            pjt = u->params[j].type;
            if (j == 0 && pjt != GRETL_TYPE_BUNDLE_REF) {
                pprintf(prn, "%s: first param type must be %s\n",
                        attr, gretl_type_get_name(GRETL_TYPE_BUNDLE_REF));
                err = E_TYPES;
            } else if (j == 1 && pjt != GRETL_TYPE_INT) {
                pprintf(prn, "%s: second param type must be %s\n",
                        attr, gretl_type_get_name(GRETL_TYPE_INT));
                err = E_TYPES;
            } else if (j == 2 && fcast && pjt != GRETL_TYPE_INT) {
                pprintf(prn, "%s: third param type must be %s\n",
                        attr, gretl_type_get_name(GRETL_TYPE_INT));
                err = E_TYPES;
            } else if (j > pmin && !fn_param_optional(u, j) &&
                       na(fn_param_default(u, j))) {
                pprintf(prn, "%s: extra params must be optional\n",
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

/* Called from the GUI package editor to check whether
   a given function of name @name can be shown as a
   candidate for a specified GUI-special @role.
*/

int function_ok_for_package_role (const char *name,
                                  int role)
{
    ufunc *u = NULL;
    int i, err = 0;

    if (name == NULL || role == UFUN_ROLE_NONE) {
        return 0;
    }

    for (i=0; i<n_ufuns; i++) {
        if (!strcmp(name, ufuns[i]->name)) {
            u = ufuns[i];
            break;
        }
    }

    if (u == NULL) {
        return 0;
    }

    if (role == UFUN_GUI_PRECHECK) {
        if (u->rettype != GRETL_TYPE_DOUBLE) {
            err = E_TYPES;
        } else if (u->n_params > 0) {
            err = E_TYPES;
        }
        return !err; /* found */
    } else if (role == UFUN_PLOT_PRECHECK) {
	if (u->rettype != GRETL_TYPE_MATRIX) {
	    err = E_TYPES;
	} else if (!valid_plot_checker(u)) {
	    err = E_TYPES;
	}
	return !err; /* found */
    } else if (role == UFUN_LIST_MAKER) {
        if (u->rettype != GRETL_TYPE_LIST) {
            err = E_TYPES;
        } else if (!valid_list_maker(u)) {
            err = E_TYPES;
        }
        return !err; /* found */
    } else if (role == UFUN_R_SETUP) {
        if (u->rettype != GRETL_TYPE_VOID) {
            err = E_TYPES;
        } else if (u->n_params > 0) {
            err = E_TYPES;
        }
        return !err; /* found */
    } else if (role == UFUN_UI_MAKER) {
        if (u->rettype != GRETL_TYPE_BUNDLE) {
            err = E_TYPES;
        } else if (u->n_params > 0) {
            err = E_TYPES;
        }
        return !err; /* found */
    }

    if (role == UFUN_GUI_MAIN) {
        ; /* OK, we don't mind what type it is */
    } else {
        /* bundle-print, bundle-plot, etc. */
        if (u->n_params == 0) {
            err = E_TYPES;
        }
        for (i=0; i<u->n_params && !err; i++) {
            if (i == 0 && u->params[i].type != GRETL_TYPE_BUNDLE_REF) {
                err = E_TYPES;
            } else if (i == 1 && u->params[i].type != GRETL_TYPE_INT) {
                err = E_TYPES;
            } else if (i > 1 && !fn_param_optional(u, i) &&
                       na(fn_param_default(u, i))) {
                err = E_TYPES;
            }
        }
    }

    return !err;
}

static int pkg_set_funcs_attribute (fnpkg *pkg, const char *s,
                                    int flag)
{
    char **S;
    int ns = 0, err = 0;

    S = gretl_string_split(s, &ns, NULL);
    if (ns == 0) {
        err = E_DATA;
    } else {
        int i, j, match;

        for (i=0; i<ns && !err; i++) {
            match = 0;
            for (j=0; j<pkg->n_pub; j++) {
                if (!strcmp(S[i], pkg->pub[j]->name)) {
                    pkg->pub[j]->flags |= flag;
                    match = 1;
                    break;
                }
            }
            if (!match) {
                err = E_DATA;
            }
        }

        strings_array_free(S, ns);
    }

    return err;
}

static void function_package_set_auxfile (fnpkg *pkg,
                                          const char *id,
                                          const char *fname)
{
    gchar *test = NULL;

    /* Maybe set the source filename for an element
       of a function package read from file, but only
       if it's not the standard, default filename for
       the element in question.
    */

    if (!strcmp(id, "help-fname")) {
        test = g_strdup_printf("%s_help.txt", pkg->name);
    } else if (!strcmp(id, "gui-help-fname")) {
        test = g_strdup_printf("%s_gui_help.txt", pkg->name);
    } else if (!strcmp(id, "sample-fname")) {
        test = g_strdup_printf("%s_sample.inp", pkg->name);
    }

    if (test != NULL) {
        if (strcmp(fname, test)) {
           function_package_set_properties(pkg, id, fname, NULL);
        }
        g_free(test);
    }
}

static int check_for_pdf_ref (const char *s, fnpkg *pkg,
                              int *err)
{
    int len, ret = 0;

    if (!strncmp(s, "pdfdoc:", 7)) {
        s += 7;
    }

    len = strlen(s);
    if (len < 64 && strchr(s, ' ') == NULL && has_suffix(s, ".pdf")) {
        char *tmp = gretl_strndup(s, len - 4);

        if (strcmp(tmp, pkg->name)) {
            gretl_errmsg_sprintf(_("PDF doc should be called %s.pdf"), pkg->name);
            *err = E_DATA;
        } else {
            ret = 1;
        }
        free(tmp);
    }

    return ret;
}

static int is_pdf_ref (const char *s)
{
    if (!strncmp(s, "pdfdoc:", 7)) {
        s += 7;
    }
    return strlen(s) < 64 && strchr(s, ' ') == NULL &&
        has_suffix(s, ".pdf");
}

static int validate_R_depends (const char *s)
{
    int err = 0;

    if (strncmp(s, "R ", 2)) {
        err = E_DATA;
    } else {
        /* skip "R " */
        s += 2;
        s += strcspn(s, " ");
        if (*s == '\0') {
            err = E_DATA;
        } else {
            /* skip R version */
            s += strcspn(s, " ");
            while (*s && !err) {
                s += strspn(s, " ");
                if (*s == '\0') break;
                /* R package? */
                s += strcspn(s, " ");
                s += strspn(s, " ");
                if (*s == '\0') {
                    /* no version given? */
                    err = E_DATA;
                } else {
                    /* skip package version */
                    s += strcspn(s, " ");
                }
            }
        }
    }

    if (err) {
        gretl_errmsg_set(_("Invalid R-depends line"));
    }

    return err;
}

static int function_package_set_version (fnpkg *pkg, const char *vstr)
{
    int err = 0;

    if (pkg->version != NULL) {
	free(pkg->version);
	pkg->version = NULL;
    }
    if (!strcmp(vstr, "@VERSION@")) {
        /* auto-versioning for gretl addons */
	if (is_gretl_addon(pkg->name)) {
	    pkg->version = gretl_strdup(GRETL_VERSION);
	} else {
	    err = E_DATA;
	}
    } else {
	pkg->version = gretl_strdup(vstr);
    }

    return err;
}

/* Having assembled and checked the function-listing for a new
   package, now retrieve the additional information from the
   spec file (named by @fname, opened as @fp).
*/

static int new_package_info_from_spec (fnpkg *pkg, const char *fname,
                                       FILE *fp, gretlopt opt,
                                       PRN *prn)
{
    const char *okstr, *failed;
    gchar *currdir = NULL;
    int quiet = (opt & OPT_Q);
    char *p, line[1024];
    int got = 0;
    int err = 0;

    if (opt & OPT_G) {
        /* GUI use */
        okstr = "<@ok>\n";
        failed = "<@fail>\n";
    } else {
        okstr = "OK\n";
        failed = "failed\n";
    }

    if (strrslash(fname) != NULL) {
        /* directory change needed */
        char dirname[FILENAME_MAX];

        strcpy(dirname, fname);
        p = strrslash(dirname);
        *p = '\0';
        currdir = g_get_current_dir();
        err = gretl_chdir(dirname);
    }

#if PKG_DEBUG
    fprintf(stderr, "new_package_info_from_spec\n");
#endif

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
	    if (!strncmp(line, "public", 6)) {
		continue;
	    } else if (!strncmp(line, "author", 6)) {
                err = function_package_set_properties(pkg, "author", p, NULL);
                if (!err) got++;
            } else if (!strncmp(line, "email", 5)) {
                err = function_package_set_properties(pkg, "email", p, NULL);
            } else if (!strncmp(line, "version", 7)) {
                err = function_package_set_version(pkg, p);
                if (!err) got++;
            } else if (!strncmp(line, "date", 4)) {
                err = function_package_set_properties(pkg, "date", p, NULL);
                if (!err) got++;
            } else if (!strncmp(line, "description", 11)) {
                err = function_package_set_properties(pkg, "description", p, NULL);
                if (!err) got++;
            } else if (!strncmp(line, "tags", 4)) {
                err = function_package_set_properties(pkg, "tags", p, NULL);
            } else if (!strncmp(line, "label", 5)) {
                err = function_package_set_properties(pkg, "label", p, NULL);
            } else if (!strncmp(line, "menu-attachment", 15)) {
                err = function_package_set_properties(pkg, "menu-attachment", p, NULL);
            } else if (!strncmp(line, "provider", 8)) {
                if (!quiet) {
                    pprintf(prn, "Recording provider %s\n", p);
                }
                err = function_package_set_properties(pkg, "provider", p, NULL);
            } else if (!strncmp(line, "help", 4)) {
                gchar *hstr = NULL;
                int pdfdoc;

                pdfdoc = check_for_pdf_ref(p, pkg, &err);
                if (pdfdoc) {
                    if (!quiet) {
                        pprintf(prn, "Recording help reference %s\n", p);
                    }
                    hstr = g_strdup_printf("pdfdoc:%s", p);
                } else if (!err) {
                    if (!quiet) {
                        pprintf(prn, "Looking for help text in %s... ", p);
                    }
                    hstr = pkg_aux_content(p, &err);
                    if (err) {
                        pputs(prn, failed);
                    } else {
                        pputs(prn, okstr);
                    }
                }
                if (!err) {
                    err = function_package_set_properties(pkg, "help", hstr, NULL);
                    if (!err) {
                        got++;
                        if (!pdfdoc) {
                            function_package_set_auxfile(pkg, "help-fname", p);
                        }
                    }
                }
                g_free(hstr);
            } else if (!strncmp(line, "gui-help", 8)) {
                gchar *ghstr = NULL;

                if (!quiet) {
                    pprintf(prn, "Looking for GUI help text in %s... ", p);
                }
                ghstr = pkg_aux_content(p, &err);
                if (err) {
                    pputs(prn, failed);
                } else {
                    pputs(prn, okstr);
                    err = function_package_set_properties(pkg, "gui-help", ghstr, NULL);
                    if (!err) {
                        function_package_set_auxfile(pkg, "gui-help-fname", p);
                    }
                }
                g_free(ghstr);
            } else if (!strncmp(line, "sample-script", 13)) {
                gchar *script = NULL;

                if (!quiet) {
                    pprintf(prn, "Looking for sample script in %s... ", p);
                }
                script = pkg_aux_content(p, &err);
                if (err) {
                    pputs(prn, failed);
                } else {
                    pputs(prn, okstr);
                    err = function_package_set_properties(pkg, "sample-script", script, NULL);
                    if (!err) {
                        got++;
                        function_package_set_auxfile(pkg, "sample-fname", p);
                    }
                }
                g_free(script);
            } else if (!strncmp(line, "translations", 12)) {
                if (!quiet) {
                    pprintf(prn, "Looking for translations in %s... ", p);
                }
                pkg->trans = read_translations_file(p, &err);
                if (err) {
                    pputs(prn, failed);
                } else {
                    pputs(prn, okstr);
                }
            } else if (!strncmp(line, "data-files", 10)) {
                if (!quiet) {
                    pprintf(prn, "Recording data-file list: %s\n", p);
                }
                err = pkg_set_datafiles(pkg, p);
            } else if (!strncmp(line, "depends", 7)) {
                if (!quiet) {
                    pprintf(prn, "Recording dependency list: %s\n", p);
                }
                err = pkg_set_depends(pkg, p);
            } else if (!strncmp(line, "R-depends", 9)) {
                err = validate_R_depends(p);
                if (!err) {
                    if (!quiet) {
                        pprintf(prn, "Recording R dependency list: %s\n", p);
                    }
                    err = function_package_set_properties(pkg, "R-depends", p, NULL);
                }
            } else if (!strncmp(line, "data-requirement", 16)) {
                err = pkg_set_dreq(pkg, p);
            } else if (!strncmp(line, "model-requirement", 17)) {
                err = pkg_set_modelreq(pkg, p);
            } else if (!strncmp(line, "min-version", 11)) {
                pkg->minver = gretl_version_number(p);
                got++;
            } else if (!strncmp(line, "lives-in-subdir", 15)) {
                pkg->uses_subdir = pkg_boolean_from_string(p);
            } else if (!strncmp(line, "wants-data-access", 17)) {
                pkg->data_access = pkg_boolean_from_string(p);
            } else if (!strncmp(line, "no-print", 8)) {
                err = pkg_set_funcs_attribute(pkg, p, UFUN_NOPRINT);
            } else if (!strncmp(line, "menu-only", 9)) {
                err = pkg_set_funcs_attribute(pkg, p, UFUN_MENU_ONLY);
            } else {
                const char *key;
		int found = 0;
                int i;

                for (i=0; pkg_lookups[i].key != NULL; i++) {
                    key = pkg_lookups[i].key;
                    if (!strncmp(line, key, strlen(key))) {
                        err = function_set_package_role(p, pkg, key, prn);
                        if (!err && !quiet) {
                            pprintf(prn, "%s function is %s, %s", key, p, okstr);
                        }
			found = 1;
                        break;
                    }
                }
		if (!found) {
		    fprintf(stderr, "*** unrecognized '%s' ***\n", line);
		}
            }
        }
    }

    if (!err && pkg->provider != NULL) {
        /* in case provider isn't registered as a dependency */
        err = strings_array_prepend_uniq(&pkg->depends, &pkg->n_depends,
                                         pkg->provider);
    }

    if (currdir != NULL) {
        /* go back where we came from */
        gretl_chdir(currdir);
        g_free(currdir);
    }

    if (!err && got < 7) {
        gretl_errmsg_set(_("Some required information was missing"));
        err = E_DATA;
    }

    return err;
}

static int is_public_funcs_line (const char *line, int *offset)
{
    const char *s = line;

    if (!strncmp(s, "public", 6)) {
        s += 6;
        s += strspn(s, " ");
        if (*s == '=') {
            *offset = 1 + s - line;
            return 1;
        }
    }

    return 0;
}

static fnpkg *new_pkg_from_spec_file (const char *gfnname, gretlopt opt,
                                      PRN *prn, int *err)
{
    fnpkg *pkg = NULL;
    char fname[FILENAME_MAX];
    char line[4096], cont[1024];
    int quiet = (opt & OPT_Q);
    FILE *fp;

    if (!has_suffix(gfnname, ".gfn")) {
        gretl_errmsg_set(_("Output must have extension \".gfn\""));
        *err = E_DATA;
        return NULL;
    }

    switch_ext(fname, gfnname, "spec");
    fp = gretl_fopen(fname, "r"); /* "rb" ? */

    if (fp == NULL) {
        *err = E_FOPEN;
    } else {
        char **pubnames = NULL;
        char **privnames = NULL;
        int npub = 0, npriv = 0;
        ufunc *uf;
        int i;

        if (!quiet) {
            pprintf(prn, _("Found spec file '%s'\n"), fname);
        }

        /* first pass: gather names of public functions */

        while (fgets(line, sizeof line, fp) && !*err) {
            int offset = 0;

            if (is_public_funcs_line(line, &offset)) {
                while (ends_with_backslash(line)) {
                    gretl_charsub(line, '\\', '\0');
                    *cont = '\0';
                    if (fgets(cont, sizeof cont, fp)) {
                        strcat(line, cont);
                    } else {
                        *err = E_PARSE;
                    }
                }
                if (!*err) {
                    tailstrip(line);
                    pubnames = gretl_string_split(line + offset, &npub, NULL);
                }
                break;
            }
        }

        if (npub == 0) {
            pprintf(prn, _("no specification of public functions was found\n"));
            *err = E_DATA;
        }

        if (!*err) {
            if (!quiet) {
                pprintf(prn, _("number of public interfaces = %d\n"), npub);
            }
            for (i=0; i<npub && !*err; i++) {
                uf = get_user_function_by_name(pubnames[i]);
                if (!quiet) {
                    pprintf(prn, " %s", pubnames[i]);
                }
                if (uf == NULL) {
                    if (!quiet) {
                        pputs(prn, _(": *** not found"));
                    }
                    gretl_errmsg_sprintf(_("'%s': no such function"), pubnames[i]);
                    *err = E_DATA;
                }
                if (!quiet) {
                    pputc(prn, '\n');
                }
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

        if (!quiet && !*err && npriv > 0) {
            pprintf(prn, "number of private functions = %d\n", npriv);
            for (i=0; i<npriv; i++) {
                pprintf(prn, " %s\n", privnames[i]);
            }
        }

        if (!*err) {
            pkg = function_package_new(gfnname, pubnames, npub,
                                       privnames, npriv, err);
        }

        strings_array_free(pubnames, npub);
        strings_array_free(privnames, npriv);

        if (!*err) {
            rewind(fp);
            *err = new_package_info_from_spec(pkg, fname, fp, opt, prn);
        }

        if (*err && pkg != NULL) {
            real_function_package_unload(pkg, 0);
            pkg = NULL;
        }

        fclose(fp);
    }

    return pkg;
}

static int cli_validate_package_file (const char *fname,
                                      gretlopt opt, PRN *prn)
{
    char dtdname[FILENAME_MAX];
    xmlDocPtr doc = NULL;
    xmlDtdPtr dtd = NULL;
    int err;

    err = gretl_xml_open_doc_root(fname, NULL, &doc, NULL);
    if (err) {
        pprintf(prn, "Couldn't parse %s\n", fname);
        return 1;
    }

    *dtdname = '\0';

    if (opt & OPT_D) {
        const char *dpath = get_optval_string(MAKEPKG, OPT_D);

        if (dpath != NULL && *dpath != '\0') {
            strcat(dtdname, dpath);
        }
    } else {
        sprintf(dtdname, "%sfunctions%cgretlfunc.dtd", gretl_home(), SLASH);
    }

    if (*dtdname != '\0') {
        dtd = xmlParseDTD(NULL, (const xmlChar *) dtdname);
    }

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
        cvp->error    = (xmlValidityErrorFunc) pprintf2;
        cvp->warning  = (xmlValidityWarningFunc) pprintf2;

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

static int get_gfn_info_for_zip (const char *fname,
                                 int *pdfdoc,
                                 char ***datafiles,
                                 int *n_datafiles)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    int found = 0;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return err;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL && found < 2) {
        if (!xmlStrcmp(cur->name, (XUC) "help")) {
            char *s = NULL;

            gretl_xml_node_get_trimmed_string(cur, doc, &s);
            *pdfdoc = is_pdf_ref(s);
            free(s);
            found++;
        } else if (!xmlStrcmp(cur->name, (XUC) "data-files")) {
            *datafiles =
                gretl_xml_get_strings_array(cur, doc, n_datafiles,
                                            0, &err);
            found++;
        } else if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
            /* we've overshot */
            found = 2;
        }
        cur = cur->next;
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    return err;
}

static int cli_build_zip_package (const char *fname,
                                  const char *gfnname,
                                  fnpkg *pkg,
                                  PRN *prn)
{
    char **datafiles = NULL;
    int n_datafiles = 0;
    int pdfdoc = 0;
    int freeit = 0;
    int err = 0;

    if (pkg != NULL) {
        datafiles = pkg->datafiles;
        n_datafiles = pkg->n_files;
        pdfdoc = is_pdf_ref(pkg->help);
    } else {
        /* grope in the gfn file */
        err = get_gfn_info_for_zip(gfnname,
                                   &pdfdoc,
                                   &datafiles,
                                   &n_datafiles);
        freeit = 1;
    }

    if (!err) {
        err = package_make_zipfile(gfnname,
                                   pdfdoc,
                                   datafiles,
                                   n_datafiles,
                                   NULL,
                                   fname,
                                   OPT_NONE,
                                   prn);
    }

    if (freeit) {
        strings_array_free(datafiles, n_datafiles);
    }

    return err;
}

static int should_rebuild_gfn (const char *gfnname)
{
    char testname[FILENAME_MAX];
    struct stat b1, b2, b3, b4;
    int err;

    err = gretl_stat(gfnname, &b1);
    if (err) {
        /* gfn not found: so have to rebuild */
        return 1;
    }

    switch_ext(testname, gfnname, "inp");
    err = gretl_stat(testname, &b2);
    if (err) {
        /* no corresponding inp: can't rebuild */
        return 0;
    }

    switch_ext(testname, gfnname, "spec");
    err = gretl_stat(testname, &b3);
    if (err) {
        /* no corresponding spec: can't rebuild */
        return 0;
    }

    switch_ext(testname, gfnname, "xml");
    err = gretl_stat(testname, &b4);
    if (!err && b4.st_mtime > b1.st_mtime) {
        /* xml is newer than gfn */
        return 1;
    }

    if (b2.st_mtime > b1.st_mtime ||
        b3.st_mtime > b1.st_mtime) {
        /* inp or spec is newer than gfn */
        return 1;
    }

    return 0;
}

static int gfn_write_translatables_to_file (fnpkg *pkg, PRN *prn)
{
    gchar *trname = g_strdup_printf("%s_xx.xml", pkg->name);
    FILE *fp = gretl_fopen(trname, "wb");
    int ret = 0;

    if (fp == NULL) {
        gretl_errmsg_sprintf(_("Couldn't open %s"), trname);
    } else {
        PRN *fprn = gretl_print_new_with_stream(fp);

        ret = package_write_translatables(pkg, fprn);
        if (ret == 0) {
            gretl_remove(trname);
        } else {
            pprintf(prn, "Wrote translations file %s\n", trname);
        }
        gretl_print_destroy(fprn);
    }

    g_free(trname);
    return ret;
}

/**
 * create_and_write_function_package:
 * @fname: filename for function package.
 * @opt: may include OPT_I to write a package-index entry,
 * OPT_T to write translatable strings (for addons). OPT_G
 * is used internally when this function is called from the
 * GUI, to inflect the informative printout.
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
    char gfnname[FILENAME_MAX];
    fnpkg *pkg = NULL;
    int build_gfn = 1;
    int build_zip = 0;
    int err = 0;

    if (has_suffix(fname, ".zip")) {
        /* building a zip package */
        switch_ext(gfnname, fname, "gfn");
        build_gfn = should_rebuild_gfn(gfnname);
        build_zip = 1;
    } else {
        /* just building a gfn file */
        strcpy(gfnname, fname);
    }

    if (build_gfn && n_free_functions() == 0) {
        gretl_errmsg_set(_("No functions are available for packaging at present."));
        err = E_DATA;
    } else if (build_gfn) {
        pkg = new_pkg_from_spec_file(gfnname, opt, prn, &err);
        if (pkg != NULL) {
            err = function_package_write_file(pkg);
            if (!err) {
                err = cli_validate_package_file(gfnname, opt, prn);
                /* should we delete @gfnname ? */
            }
            if (!err) {
                if (opt & OPT_T) {
                    if (is_gretl_addon(pkg->name)) {
                        /* addons are a special case */
                        addon_write_translatable_strings(pkg, prn);
                    } else {
                        gfn_write_translatables_to_file(pkg, prn);
                    }
                }
                if (opt & OPT_I) {
                    package_write_index(pkg, prn);
                }
            }
        }
    }

    if (!err && build_zip) {
        err = cli_build_zip_package(fname, gfnname, pkg, prn);
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

/**
 * function_package_get_version:
 * @pkg: function package.
 *
 * Returns: the name of the package.
 */

double function_package_get_version (fnpkg *pkg)
{
    if (pkg == NULL) {
        return NADBL;
    } else {
        return dot_atof(pkg->version);
    }
}

static int maybe_replace_string_var (char **svar, const char *src)
{
    if (src == NULL) {
        gretl_errmsg_set(_("string value is missing"));
        return E_DATA;
    } else {
        free(*svar);
        *svar = gretl_strdup(src);
        return (*svar == NULL)? E_ALLOC : 0;
    }
}

/* unlike the case above, here we'll accept NULL for @src, to wipe
   out an existing string var */

static int maybe_replace_optional_string_var (char **svar,
                                              const char *src)
{
    if (src == NULL) {
        free(*svar);
        *svar = NULL;
        return 0;
    } else {
        free(*svar);
        *svar = gretl_strdup(src);
        return (*svar == NULL)? E_ALLOC : 0;
    }
}

/* Called (indirectly) from GUI function packager. Note that we're
   careful not to touch UFUN_PRIVATE on the uf->flags side, since
   this flag is not represented in @attrs.
*/

static void pkg_set_gui_attrs (fnpkg *pkg, const unsigned char *attrs)
{
    ufunc *uf;
    int i, r;

    for (r=1; r<UFUN_GUI_PRECHECK; r++) {
        uf = NULL;
        for (i=0; i<n_ufuns; i++) {
            if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == r) {
                uf = ufuns[i];
                break;
            }
        }
        if (uf != NULL) {
            if (attrs[r-1] & UFUN_NOPRINT) {
                uf->flags |= UFUN_NOPRINT;
            } else {
                uf->flags &= ~UFUN_NOPRINT;
            }
            if (attrs[r-1] & UFUN_MENU_ONLY) {
                uf->flags |= UFUN_MENU_ONLY;
            } else {
                uf->flags &= ~UFUN_MENU_ONLY;
            }
        }
    }
}

/* Called (indirectly) from GUI function packager. Note that we scrub
   UFUN_PRIVATE from the user function flags passed to the packager:
   this flag is handled by a different mechanism.
*/

static void pkg_get_gui_attrs (fnpkg *pkg, unsigned char *attrs)
{
    ufunc *uf;
    UfunAttrs a;
    int i, r;

    for (r=1; r<UFUN_GUI_PRECHECK; r++) {
        uf = NULL;
        a = 0;
        for (i=0; i<n_ufuns; i++) {
            if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == r) {
                uf = ufuns[i];
                break;
            }
        }
        if (uf != NULL) {
            if (uf->flags & UFUN_NOPRINT) {
                a |= UFUN_NOPRINT;
            }
            if (uf->flags & UFUN_MENU_ONLY) {
                a |= UFUN_MENU_ONLY;
            }
        }
        attrs[r-1] = a;
    }
}

static char **pkg_strvar_pointer (fnpkg *pkg, const char *key,
                                  int *optional)
{
    *optional = 0;

    if (!strcmp(key, "fname")) {
        return &pkg->fname;
    } else if (!strcmp(key, "author")) {
        return &pkg->author;
    } else if (!strcmp(key, "email")) {
        return &pkg->email;
    } else if (!strcmp(key, "version")) {
        return &pkg->version;
    } else if (!strcmp(key, "date")) {
        return &pkg->date;
    } else if (!strcmp(key, "description")) {
        return &pkg->descrip;
    } else if (!strcmp(key, "help")) {
        return &pkg->help;
    } else if (!strcmp(key, "sample-script")) {
        return &pkg->sample;
    }

    *optional = 1;

    if (!strcmp(key, "tags")) {
        return &pkg->tags; /* FIXME should be non-optional */
    } else if (!strcmp(key, "label")) {
        return &pkg->label;
    } else if (!strcmp(key, "menu-attachment")) {
        return &pkg->mpath;
    } else if (!strcmp(key, "gui-help")) {
        return &pkg->gui_help;
    } else if (!strcmp(key, "R-depends")) {
        return &pkg->Rdeps;
    } else if (!strcmp(key, "help-fname")) {
        return &pkg->help_fname;
    } else if (!strcmp(key, "gui-help-fname")) {
        return &pkg->gui_help_fname;
    } else if (!strcmp(key, "sample-fname")) {
        return &pkg->sample_fname;
    } else if (!strcmp(key, "provider")) {
        return &pkg->provider;
    }

    return NULL;
}

/* varargs function for setting the properties of a function
   package: the settings take the form of a NULL-terminated
   set of (key, value) pairs.
*/

int function_package_set_properties (fnpkg *pkg, ...)
{
    va_list ap;
    const char *key;
    char **sptr;
    int optional;
    int err = 0;

    va_start(ap, pkg);

    while (err == 0) {
        key = va_arg(ap, const char *);
        if (key == NULL) {
            break;
        }

        sptr = pkg_strvar_pointer(pkg, key, &optional);

        if (sptr != NULL) {
            const char *sval = va_arg(ap, const char *);

            if (optional) {
                err = maybe_replace_optional_string_var(sptr, sval);
            } else {
                err = maybe_replace_string_var(sptr, sval);
            }
            if (!err && !strcmp(key, "help")) {
                if (!strncmp(sval, "pdfdoc", 6) ||
                    is_pdf_ref(sval)) {
                    pkg->uses_subdir = 1;
                }
            }
        } else if (!strcmp(key, "gui-attrs")) {
            const unsigned char *np = va_arg(ap, const unsigned char *);

            pkg_set_gui_attrs(pkg, np);
        } else {
            int ival = va_arg(ap, int);

            if (!strcmp(key, "data-requirement")) {
                pkg->dreq = ival;
            } else if (!strcmp(key, "model-requirement")) {
                pkg->modelreq = ival;
            } else if (!strcmp(key, "min-version")) {
                pkg->minver = ival;
            } else if (!strcmp(key, "lives-in-subdir")) {
                pkg->uses_subdir = (ival != 0);
            } else if (!strcmp(key, "wants-data-access")) {
                pkg->data_access = (ival != 0);
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
    int j = 0;

    if (n > 0) {
        list = gretl_list_new(n);
        if (list != NULL) {
            int i, priv, menu_only;

            for (i=0; i<n_ufuns; i++) {
                if (ufuns[i]->pkg == pkg) {
                    priv = function_is_private(ufuns[i]);
                    menu_only = function_is_menu_only(ufuns[i]);
                    if (code == PRIVLIST && priv) {
                        list[++j] = i;
                    } else if (code == PUBLIST && !priv) {
                        list[++j] = i;
                    } else if (code == GUILIST && !priv && !menu_only &&
                               !pkg_aux_role(ufuns[i]->pkg_role)) {
                        /* in the GUI list of public functions, don't
                           display post-processing functions
                        */
                        list[++j] = i;
                    }
                }
            }
        }
    }

    if (list != NULL) {
        if (j == 0) {
            free(list);
            list = NULL;
        } else {
            list[0] = j;
        }
    }

    return list;
}

static gchar *pkg_get_special_func_name (fnpkg *pkg, UfunRole role)
{
    int i;

    for (i=0; i<pkg->n_pub; i++) {
        if (pkg->pub[i]->pkg_role == role) {
            return g_strdup(pkg->pub[i]->name);
        }
    }
    for (i=0; i<pkg->n_priv; i++) {
        if (pkg->priv[i]->pkg_role == role) {
            return g_strdup(pkg->priv[i]->name);
        }
    }

    return NULL;
}

static int pkg_get_special_func_id (fnpkg *pkg, UfunRole role)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i]->pkg == pkg && ufuns[i]->pkg_role == role) {
            return i;
        }
    }

    return -1;
}

static void handle_optional_string (char **ps, const char *src)
{
    if (src == NULL) {
        *ps = NULL;
    } else {
        *ps = g_strdup(src);
    }
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
    UfunRole role;
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
        } else if (!strcmp(key, "email")) {
            ps = (char **) ptr;
            if (pkg->email != NULL) {
                *ps = g_strdup(pkg->email);
            } else {
                *ps = g_strdup(""); /* ? */
            }
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
            handle_optional_string(ps, pkg->gui_help);
        } else if (!strcmp(key, "R-depends")) {
            ps = (char **) ptr;
            handle_optional_string(ps, pkg->Rdeps);
        } else if (!strcmp(key, "sample-script")) {
            ps = (char **) ptr;
            *ps = g_strdup(pkg->sample);
        } else if (!strcmp(key, "help-fname")) {
            ps = (char **) ptr;
            handle_optional_string(ps, pkg->help_fname);
        } else if (!strcmp(key, "gui-help-fname")) {
            ps = (char **) ptr;
            handle_optional_string(ps, pkg->gui_help_fname);
        } else if (!strcmp(key, "sample-fname")) {
            ps = (char **) ptr;
            handle_optional_string(ps, pkg->sample_fname);
        } else if (!strcmp(key, "tags")) {
            ps = (char **) ptr;
            *ps = g_strdup(pkg->tags);
        } else if (!strcmp(key, "label")) {
            ps = (char **) ptr;
            *ps = g_strdup(pkg->label);
        } else if (!strcmp(key, "menu-attachment")) {
            ps = (char **) ptr;
            *ps = g_strdup(pkg->mpath);
        } else if (!strcmp(key, "provider")) {
            ps = (char **) ptr;
            *ps = g_strdup(pkg->provider);
        } else if (!strcmp(key, "data-requirement")) {
            pi = (int *) ptr;
            *pi = pkg->dreq;
        } else if (!strcmp(key, "model-requirement")) {
            pi = (int *) ptr;
            *pi = pkg->modelreq;
        } else if (!strcmp(key, "min-version")) {
            pi = (int *) ptr;
            *pi = pkg->minver;
        } else if (!strcmp(key, "lives-in-subdir")) {
            pi = (int *) ptr;
            *pi = pkg->uses_subdir;
        } else if (!strcmp(key, "wants-data-access")) {
            pi = (int *) ptr;
            *pi = pkg->data_access;
        } else if (!strcmp(key, "publist")) {
            plist = (int **) ptr;
            *plist = function_package_get_list(pkg, PUBLIST, npub);
        } else if (!strcmp(key, "gui-publist")) {
            plist = (int **) ptr;
            *plist = function_package_get_list(pkg, GUILIST, npub);
        } else if (!strcmp(key, "privlist")) {
            plist = (int **) ptr;
            *plist = function_package_get_list(pkg, PRIVLIST, npriv);
        } else if (!strcmp(key, "gui-main-id")) {
            pi = (int *) ptr;
            *pi = pkg_get_special_func_id(pkg, UFUN_GUI_MAIN);
        } else if (!strcmp(key, "gui-attrs")) {
            unsigned char *s = (unsigned char *) ptr;

            pkg_get_gui_attrs(pkg, s);
        } else if ((role = pkg_key_get_role(key))) {
            ps = (char **) ptr;
            *ps = pkg_get_special_func_name(pkg, role);
        }
    }

    va_end(ap);

    return err;
}

/* don't tamper with return value! */

const char *function_package_get_string (fnpkg *pkg,
                                         const char *id)
{
    if (pkg == NULL || id == NULL) {
        return NULL;
    } else if (!strcmp(id, "fname")) {
        return pkg->fname;
    } else if (!strcmp(id, "help-fname")) {
        return pkg->help_fname;
    } else if (!strcmp(id, "gui-help-fname")) {
        return pkg->gui_help_fname;
    } else if (!strcmp(id, "sample-fname")) {
        return pkg->sample_fname;
    } else if (!strcmp(id, "sample-script")) {
        return pkg->sample;
    } else if (!strcmp(id, "help")) {
        return pkg->help;
    } else if (!strcmp(id, "gui-help")) {
        return pkg->gui_help;
    } else if (!strcmp(id, "R-depends")) {
        return pkg->Rdeps;
    } else {
        return NULL;
    }
}

char **function_package_get_data_files (fnpkg *pkg, int *n)
{
    char **S = NULL;

    *n = 0;
    if (pkg->datafiles != NULL) {
        S = strings_array_dup(pkg->datafiles, pkg->n_files);
        if (S != NULL) {
            *n = pkg->n_files;
        }
    }

    return S;
}

char **function_package_get_depends (fnpkg *pkg, int *n)
{
    char **S = NULL;

    *n = 0;
    if (pkg->depends != NULL) {
        S = strings_array_dup(pkg->depends, pkg->n_depends);
        if (S != NULL) {
            *n = pkg->n_depends;
        }
    }

    return S;
}

int function_package_set_data_files (fnpkg *pkg, char **S, int n)
{
    int err = 0;

    if (pkg->datafiles != NULL) {
        strings_array_free(pkg->datafiles, pkg->n_files);
        pkg->datafiles = NULL;
        pkg->n_files = 0;
    }

    if (n > 0) {
        if (S == NULL) {
            err = E_DATA;
        } else {
            pkg->datafiles = strings_array_dup(S, n);
            if (pkg->datafiles == NULL) {
                err = E_ALLOC;
            } else {
                pkg->n_files = n;
                pkg->uses_subdir = 1;
            }
        }
    }

    return err;
}

int function_package_set_depends (fnpkg *pkg, char **S, int n)
{
    int err = 0;

    if (pkg->depends != NULL) {
        strings_array_free(pkg->depends, pkg->n_depends);
        pkg->depends = NULL;
        pkg->n_depends = 0;
    }

    if (n > 0) {
        if (S == NULL) {
            err = E_DATA;
        } else {
            pkg->depends = strings_array_dup(S, n);
            if (pkg->depends == NULL) {
                err = E_ALLOC;
            } else {
                pkg->n_depends = n;
            }
        }
    }

    if (!err && pkg->provider != NULL) {
        err = strings_array_prepend_uniq(&pkg->depends,
                                         &pkg->n_depends,
                                         pkg->provider);
    }

    return err;
}

/* quick check to see if there's a gross problem with a package,
   in the context of considering packing it into a gretl
   session file
*/

static int validate_function_package (fnpkg *pkg)
{
    if (pkg->pub == NULL || pkg->author == NULL ||
        pkg->version == NULL || pkg->date == NULL ||
        pkg->descrip == NULL) {
        return 0;
    } else if (pkg->name[0] == '\0') {
        return 0;
    }

    return 1;
}

/* for session file use, and also MPI: dump all currently defined
   functions, packaged or not, into a single XML file */

int write_loaded_functions_file (const char *fname, int mpicall)
{
    PRN *prn;
    int i, err = 0;

    if (n_ufuns == 0) {
        return 0;
    }

    prn = gretl_print_new_with_filename(fname, &err);
    if (prn == NULL) {
        return err;
    }

    gretl_xml_header(prn);
    pputs(prn, "<gretl-functions>\n");

#ifdef HAVE_MPI
    if (mpicall) {
        /* if we're launching MPI, record the name of the
           currently executing function, if any
        */
        ufunc *u = currently_called_function();

        if (u != NULL) {
            pprintf(prn, "<caller>%s</caller>\n", u->name);
        }
    }
#endif

    /* write any loaded function packages */
    for (i=0; i<n_pkgs; i++) {
        if (validate_function_package(pkgs[i])) {
            real_write_function_package(pkgs[i], prn, mpicall);
        }
    }

    /* then any unpackaged functions */
    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i]->pkg == NULL) {
            write_function_xml(ufuns[i], prn, mpicall);
        }
    }

    pputs(prn, "</gretl-functions>\n");

    gretl_print_destroy(prn);

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

        if (pkg->datafiles != NULL && pkg->n_files > 0) {
            strings_array_free(pkg->datafiles, pkg->n_files);
        }
        if (pkg->depends != NULL && pkg->n_depends > 0) {
            strings_array_free(pkg->depends, pkg->n_depends);
        }
        if (pkg->trans != NULL) {
            destroy_translation(pkg->trans);
        }

        free(pkg->pub);
        free(pkg->priv);
        free(pkg->fname);
        free(pkg->author);
        free(pkg->email);
        free(pkg->version);
        free(pkg->date);
        free(pkg->descrip);
        free(pkg->help);
        free(pkg->gui_help);
        free(pkg->Rdeps);
        free(pkg->sample);
        free(pkg->help_fname);
        free(pkg->gui_help_fname);
        free(pkg->sample_fname);
        free(pkg->tags);
        free(pkg->label);
        free(pkg->mpath);
        free(pkg->provider);
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

static fnpkg *get_loaded_pkg_by_filename (const char *fname,
                                          const char **version)
{
    int i;

    if (fname == NULL) {
        return NULL;
    }

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(fname, pkgs[i]->fname)) {
            if (version != NULL) {
                *version = pkgs[i]->version;
            }
            return pkgs[i];
        }
    }

    return NULL;
}

static fnpkg *get_loaded_pkg_by_name (const char *pkgname)
{
    int i;

    if (pkgname == NULL) {
        return NULL;
    }

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(pkgname, pkgs[i]->name)) {
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
    fnpkg *pkg = get_loaded_pkg_by_filename(fname, NULL);

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
    fnpkg *pkg = get_loaded_pkg_by_filename(fname, NULL);

    if (pkg != NULL) {
        real_function_package_unload(pkg, 1);
    }
}

/**
 * function_package_unload_full:
 * @pkgname: package name.
 *
 * Unloads the specified function package from memory, if it
 * is currently loaded.  The functions 'owned' by the package
 * are also unloaded from memory.
 *
 * Returns: 1 if the specified package was loaded and is
 * new unloaded, otherwise 0.
 */

int function_package_unload_full (const char *pkgname)
{
    fnpkg *pkg = get_loaded_pkg_by_name(pkgname);

    if (pkg != NULL) {
        real_function_package_unload(pkg, 1);
        return 1;
    } else {
        return 0;
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
    gretl_errmsg_sprintf(_("'%s': package contains "
                         "duplicated function names"),
                         pkg->name);
    return E_DATA;
}

#if USE_RLIB

/* Try defining and calling a trivial R function to check that
   this functionality is actually working. FIXME: if it fails,
   can we diagnose why exactly that happens?
*/

static int verify_libR_functionality (PRN *prn)
{
    const char *buf =
        "<gretl-function name=\"Rlib_check_\" type=\"scalar\">\n"
        "<code>foreign language=R\n"
        "plusone &lt;- function(x) {\n"
        "x + 1\n"
        "}\n"
        "end foreign\n"
        "catch scalar a = R.plusone(1)\n"
        "return $error\n"
        "</code>\n"
        "</gretl-function>\n";
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    int err;

    err = gretl_xml_read_buffer(buf, "gretl-function", &doc, &node);
    if (!err) {
        err = read_ufunc_from_xml(node, doc, NULL);
    }
    if (!err) {
        /* @ival should be 0 on success */
        int ival = generate_int("Rlib_check_()", NULL, &err);

        if (!err && ival) {
            fprintf(stderr, "Rlib_check_ call: return value %d\n", ival);
            err = ival;
        }
    }
    if (doc != NULL) {
        xmlFreeDoc(doc);
    }
    if (err) {
        pputs(prn, _("Warning: R_functions could not be enabled\n"));
    }

    return err;
}

/* For an R-dependent package that wants to define one or more R
   functions: extract the relevant lines of code from its dedicated
   R-setup function (which must be private) and send them to the R
   library for "compilation". If this fails then we can't use R_lib.
*/

static int package_run_R_setup (ufunc *fun, PRN *prn)
{
    char *buf = NULL;
    PRN *fprn = NULL;
    int i, err;

    fprintf(stderr, "start package_run_R_setup\n");

    /* first screen for nominal availability of libR */
    err = libset_set_bool(R_FUNCTIONS, 1);

    if (!err) {
        /* then check that a simple example works */
        err = verify_libR_functionality(prn);
        if (err) {
            libset_set_bool(R_FUNCTIONS, 0);
            return 1;
        }
    }

    if (!err) {
        fprn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    }

    if (fprn != NULL) {
        const char *s;

        for (i=0; i<fun->n_lines; i++) {
            s = fun->lines[i].s;
            /* skip "foreign" and "end foreign" */
            if (!strncmp(s, "foreign lan", 11) ||
                !strncmp(s, "end foreign", 11)) {
                continue;
            } else if (fun->lines[i].flags & LINE_IGNORE) {
                continue;
            } else {
                pputs(fprn, s);
                pputc(fprn, '\n');
            }
        }
        buf = gretl_print_steal_buffer(fprn);
        gretl_print_destroy(fprn);
    }

    if (buf != NULL) {
        err = execute_R_buffer(buf, NULL, OPT_NONE, NULL);
        free(buf);
    }

    fprintf(stderr, "package_run_R_setup: err = %d\n", err);

    return err;
}

#endif /* USE_RLIB */

static void do_override (ufunc *fun)
{
    gretl_warnmsg_sprintf(_("'%s' is the name of a built-in function"),
                          fun->name);
    install_function_override(fun->name, fun->pkg->name, fun);
    fun->pkg->overrides += 1;
}

/* When loading a private function the only real conflict would be
   with a function of the same name owned by the same package.
   Obviously this shouldn't happen but we'll whack it if it does.
*/

static int load_private_function (fnpkg *pkg, int i, PRN *prn)
{
    ufunc *fun = pkg->priv[i];
    int j, err = 0;

    for (j=0; j<n_ufuns; j++) {
        if (!strcmp(fun->name, ufuns[j]->name)) {
            if (pkg == ufuns[j]->pkg) {
                return broken_package_error(pkg);
            }
        }
    }

#if USE_RLIB
    if (fun->pkg_role == UFUN_R_SETUP) {
        err = package_run_R_setup(fun, prn);
    }
#endif

    if (!err) {
        err = add_allocated_ufunc(fun);
        if (!err && function_lookup(fun->name)) {
            do_override(fun);
        }
    }

    return err;
}

/* When loading a public, packaged function we want to avoid conflicts
   with any non-private functions of the same name.  In case we get a
   conflict with a public member of an already loaded distinct package
   we'll flag an error (otherwise things will get too confusing).
*/

static int load_public_function (fnpkg *pkg, int i)
{
    ufunc *fun = pkg->pub[i];
    int j, done = 0;
    int err = 0;

    for (j=0; j<n_ufuns; j++) {
        if (!strcmp(fun->name, ufuns[j]->name)) {
            if (pkg == ufuns[j]->pkg) {
                /* name duplication in single package */
                err = broken_package_error(pkg);
            } else if (ufuns[j]->pkg == NULL) {
                /* conflicting unpackaged function */
                ufunc_free(ufuns[j]);
                ufuns[j] = fun;
                done = 1;
            } else if (!function_is_private(ufuns[j])) {
                /* got a conflicting package */
                gretl_errmsg_sprintf(_("The function %s is already defined "
                                     "by package '%s'"), fun->name,
                                     ufuns[j]->pkg->name);
                err = E_DATA;
                break;
            }
        }
    }

    if (!err && !done && function_lookup(fun->name)) {
        if (!strcmp(pkg->name, "PanelTools") &&
            !strcmp(fun->name, "pxmean")) {
            /* give this a pass: just use the built-in pxmean() */
            return 0;
        } else if (strcmp(fun->name, "bkw")) {
            /* for now, don't throw an error on loading Lee Adkins'
               bkw package */
            gretl_errmsg_sprintf(_("Package loading error:\n"
                                   "%s %s defines a function %s(), "
                                   "which conflicts with a built-in\n"
                                   "function in gretl %s. Please check for an "
                                   "updated version of %s."),
                                 pkg->name, pkg->version, fun->name,
                                 GRETL_VERSION, pkg->name);
            err = E_DATA;
        }
    }

    if (!err && !done) {
        err = add_allocated_ufunc(fun);
    }

    return err;
}

static int pkg_in_stack (const char *name, GArray *pstack)
{
    char *s;
    int i;

    for (i=0; i<pstack->len; i++) {
        s = g_array_index(pstack, char *, i);
        if (!strcmp(name, s)) {
            return 1;
        }
    }

    return 0;
}

#define TS_REQUIRED(d) (d == FN_NEEDS_TS || d == FN_NEEDS_QM)

static int data_requirements_conflict (fnpkg *a, fnpkg *b)
{
    if (a->dreq == FN_NEEDS_PANEL && TS_REQUIRED(b->dreq)) {
        return 1;
    } else if (b->dreq == FN_NEEDS_PANEL && TS_REQUIRED(a->dreq)) {
        return 1;
    }
    return 0;
}

static int load_gfn_dependencies (fnpkg *pkg, GArray *pstack, int level)
{
    int err = 0;

    if (pkg->depends != NULL) {
        char *pkgpath;
        int i;

        fprintf(stderr, "*** load_gfn_dependencies (%d) for %s (level %d) ***\n",
                pkg->n_depends, pkg->name, level);

        for (i=0; i<pkg->n_depends && !err; i++) {
            const char *depname = pkg->depends[i];
            fnpkg *dep = NULL;

            if (get_function_package_by_name(depname) != NULL) {
                fprintf(stderr, " %s: already loaded\n", depname);
                ; /* OK, already loaded */
            } else if (pkg_in_stack(depname, pstack)) {
                fprintf(stderr, " %s: found in stack\n", depname);
                ; /* don't go into infinite loop! */
            } else {
                fprintf(stderr, " trying for %s\n", depname);
                pkgpath = gretl_function_package_get_path(depname, PKG_ALL);
                if (pkgpath == NULL) {
                    err = E_DATA;
                    gretl_errmsg_sprintf(_("%s: dependency %s was not found"),
                                         pkg->name, depname);
                } else {
                    err = load_function_package(pkgpath, OPT_NONE,
                                                pstack, NULL, &dep,
                                                level + 1);
                    free(pkgpath);
                    if (!err && level == 0) {
                        err = data_requirements_conflict(pkg, dep);
                    }
                    if (!err) {
                        g_array_append_val(pstack, depname);
                    }
                }
            }
        }
    }

    return err;
}

/* A 'real load' is in contrast to just reading some info from a
   package, as is done in various GUI contexts.  We do the real load
   in response to some GUI commands, the "include" command, and also
   when re-opening a gretl session file that contains function
   packages.
*/

static int real_load_package (fnpkg *pkg, GArray *pstack, PRN *prn,
                              int level)
{
    int i, err = 0;

#if PKG_DEBUG
    fprintf(stderr, "real_load_package:\n loading '%s'\n", pkg->fname);
#endif

    gretl_error_clear();

    if (pstack != NULL) {
        err = load_gfn_dependencies(pkg, pstack, level);
    }

    if (!err && pkg->pub != NULL) {
        for (i=0; i<pkg->n_pub && !err; i++) {
            err = load_public_function(pkg, i);
        }
    }

    if (!err && pkg->priv != NULL) {
        for (i=0; i<pkg->n_priv && !err; i++) {
            err = load_private_function(pkg, i, prn);
        }
    }

    if (!err && pkg->provider != NULL) {
        /* check that provider really got loaded */
        if (get_function_package_by_name(pkg->provider) == NULL) {
            gretl_errmsg_sprintf(_("Provider package %s is not loaded\n"),
                                 pkg->provider);
            err = E_DATA;
        }
    }

    if (!err) {
        /* add to array of loaded packages */
        err = function_package_record(pkg);
    }

    return err;
}

static const char *data_needs_string (DataReq dr)
{
    if (dr == FN_NEEDS_TS) {
        return N_("Time-series data");
    } else if (dr == FN_NEEDS_QM) {
        return N_("Quarterly or monthly data");
    } else if (dr == FN_NEEDS_PANEL) {
        return N_("Panel data");
    } else if (dr == FN_NODATA_OK) {
        return N_("none");
    } else {
        return N_("some sort of dataset");
    }
}

static void pkg_print_depends (const fnpkg *pkg,  PRN *prn)
{
    int rdep = pkg->Rdeps != NULL;
    int i, n = pkg->n_depends + rdep;

    pputs(prn, "<@itl=\"Dependencies\">: ");
    for (i=0; i<pkg->n_depends; i++) {
        pprintf(prn, "%s%s", pkg->depends[i],
                (i < n-1)? ", " : "\n");
    }
    if (rdep) {
        pprintf(prn, "%s%s\n", pkg->n_depends ? "," : "",
                pkg->Rdeps);
    }
}

int help_text_is_markdown (const fnpkg *pkg, int gui_help)
{
    const char *fname;
    const char *text;

    text = gui_help ? pkg->gui_help : pkg->help;
    if (text == NULL) {
        return 0;
    }

    fname = gui_help ? pkg->gui_help_fname : pkg->help_fname;
    if (fname != NULL && has_suffix(fname, ".md")) {
        /* A filename matching "*.md" should clinch it. */
        return 1;
    }

    if (pkg->date != NULL) {
        /* Use of markdown in gfn help text was first supported in
           gretl 2023b (2023-07-21); we assume that packages with
           an earlier last-modified date can't be using it.
        */
        int y, m, d;

        if (sscanf(pkg->date, "%4d-%02d-%02d", &y, &m, &d) == 3) {
            if (10000*y + 100*m + d < 20230721) {
                return 0;
            }
        }
    }

    /* We might have anonymous help text (no filename specified),
       in which case we fall back on a simple heuristic.
    */
    return simple_markdown_detect(text);
}

static void print_package_info (const fnpkg *pkg, const char *fname, PRN *prn)
{
    char vstr[8];
    int remote, pdfdoc;
    int md = 0;
    int i;

    remote = (strstr(fname, "dltmp.") != NULL);

    if (!remote && g_path_is_absolute(fname)) {
        pprintf(prn, "<@itl=\"File\">: %s\n\n", fname);
    }

    if (pkg->name[0] == '\0' || pkg->author == NULL ||
        pkg->minver <= 0 || pkg->descrip == NULL ||
        pkg->version == NULL || pkg->date == NULL ||
        pkg->help == NULL) {
        pprintf(prn, _("\nBroken package! Basic information is missing\n"));
        return;
    }

    gretl_version_string(vstr, pkg->minver);
    pdfdoc = is_pdf_ref(pkg->help);

    pprintf(prn, "<@itl=\"Package\">: %s %s (%s)\n", pkg->name, pkg->version,
            pkg->date);
    pprintf(prn, "<@itl=\"Author\">: %s\n", pkg->author);
    if (pkg->email != NULL && *pkg->email != '\0') {
        pprintf(prn, "<@itl=\"Email\">: %s\n", pkg->email);
    }
    pprintf(prn, "<@itl=\"Required gretl version\">: %s\n", vstr);
    pprintf(prn, "<@itl=\"Data requirement\">: %s\n", _(data_needs_string(pkg->dreq)));
    pprintf(prn, "<@itl=\"Description\">: %s\n", gretl_strstrip(pkg->descrip));
    if (pkg->n_depends > 0 || pkg->Rdeps != NULL) {
        pkg_print_depends(pkg, prn);
    }
    if (pkg->provider != NULL) {
        pprintf(prn, "<@itl=\"Provider\">: %s\n", pkg->provider);
    }

    if (pdfdoc) {
        const char *s = strrchr(pkg->help, ':');

        if (remote) {
            pprintf(prn, "<@itl=\"Documentation\">: %s\n\n", s + 1);
        } else {
            gchar *localpdf = g_strdup(fname);
            gchar *p = strrchr(localpdf, '.');

            *p = '\0';
            strcat(p, ".pdf");
            pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n\n", localpdf);
            g_free(localpdf);
        }
    } else {
        pputc(prn, '\n');
    }

    if (pkg->pub != NULL) {
	/* note: pkg->pub is non-NULL only if the package is loaded */
        if (pkg->n_pub == 1) {
            if (strcmp(pkg->pub[0]->name, pkg->name)) {
                pputs(prn, "<@itl=\"Public interface\">: ");
                pputs(prn, "<mono>\n");
                pprintf(prn, "%s()\n", pkg->pub[0]->name);
                pputs(prn, "</mono>\n\n");
            }
        } else {
            pputs(prn, "<@itl=\"Public interfaces\">:\n\n");
            pputs(prn, "<mono>\n");
            for (i=0; i<pkg->n_pub; i++) {
                pprintf(prn, "  %s()\n", pkg->pub[i]->name);
            }
            pputs(prn, "</mono>\n\n");
        }
    }

    md = help_text_is_markdown(pkg, 0);

    if (!pdfdoc) {
        if (md) {
            md_to_gretl(pkg->help, prn);
        } else {
            pputs(prn, "<@itl=\"Help text\">:\n\n");
            pputs(prn, "<mono>\n");
            pputs(prn, pkg->help);
            pputs(prn, "\n\n");
            pputs(prn, "</mono>\n");
        }
    }

    if (remote && pkg->sample != NULL) {
        pputs(prn, "<@itl=\"Sample script\">:\n\n");
        pputs(prn, "<code>\n");
        pputs(prn, pkg->sample);
        pputs(prn, "</code>\n");
        pputc(prn, '\n');
    }
}

/* start apparatus for augmenting "pkg query" to cover contents
   of the "examples" directory of a function package, as of
   February/March 2024.
*/

static int resource_wanted (const char *s)
{
    return has_suffix(s, ".inp") ||
	has_suffix(s, ".gdt") ||
	has_suffix(s, ".gdtb");
}

static int pkg_has_examples (const fnpkg *pkg)
{
    int i;

    for (i=0; i<pkg->n_files; i++) {
        if (!strcmp(pkg->datafiles[i], "examples")) {
            return 1;
        }
    }

    return 0;
}

static int pkg_has_datafiles (const fnpkg *pkg)
{
    int i;

    for (i=0; i<pkg->n_files; i++) {
        if (strlen(pkg->datafiles[i]) > 4 &&
	    has_suffix(pkg->datafiles[i], ".gdt")) {
            return 1;
        }
    }

    return 0;
}

/* helper for compare_resources() below */

static int same_subdir (const char *sa, const char *sb)
{
    int n = strcspn(sa, SLASHSTR);

    if (strcspn(sb, SLASHSTR) != n) {
	return 0;
    } else {
	return strncmp(sa, sb, n) == 0;
    }
}

/* qsort callback for resource filenames */

static int compare_resources (const void *a, const void *b)
{
    const char *sa = *(const char **) a;
    const char *sb = *(const char **) b;
    int ac = strchr(sa, SLASH) != NULL;
    int bc = strchr(sb, SLASH) != NULL;

    if (ac + bc == 1) {
	/* sort top-level files before subdir contents */
	return ac ? 1 : -1;
    } else {
	/* sort scripts before data files, per directory */
	int same_dir = (ac + bc == 0) ? 1 : same_subdir(sa, sb);

	if (same_dir) {
	    ac = has_suffix(sa, ".inp");
	    bc = has_suffix(sb, ".inp");
	    if (ac + bc == 1) {
		return ac ? -1 : 1;
	    }
	}
	/* just alphabetize */
	return strcmp(sa, sb);
    }
}

static char **package_resources_info (const fnpkg *pkg,
                                      const char *fname,
                                      char resdir[],
                                      int *pn_res)
{
    char **resources = NULL;
    gchar *pkg_base = NULL;
    const gchar *dname;
    const gchar *subdir;
    char tmp[128];
    gchar *savedir = NULL;
    GDir *dir = NULL;
    GDir *child = NULL;
    int path_len = 0;
    int n_res = 0;

    path_len = strlen(fname) - strlen(pkg->name) - 4;
    pkg_base = g_strndup(fname, path_len);

    gretl_build_path(resdir, pkg_base, "examples", NULL);
    dir = gretl_opendir(resdir);
    if (dir == NULL) {
	g_free(pkg_base);
	return NULL;
    }

    /* cd into the examples directory to make life easier */
    savedir = g_get_current_dir();
    gretl_chdir(resdir);

    /* first get a count of relevant files */
    while ((dname = g_dir_read_name(dir)) != NULL) {
	if (g_file_test(dname, G_FILE_TEST_IS_DIR)) {
	    child = gretl_opendir(dname);
	    if (child != NULL) {
		while ((dname = g_dir_read_name(child)) != NULL) {
		    if (resource_wanted(dname)) {
			n_res++;
		    }
		}
		g_dir_close(child);
	    }
	} else if (resource_wanted(dname)) {
	    n_res++;
	}
    }

    if (n_res > 0) {
	/* if there are any relevant files, create an array */
	int i = 0;

	resources = strings_array_new(n_res);
	g_dir_rewind(dir);

	while ((dname = g_dir_read_name(dir)) != NULL) {
	    if (g_file_test(dname, G_FILE_TEST_IS_DIR)) {
		child = gretl_opendir(dname);
		if (child != NULL) {
		    subdir = dname;
		    while ((dname = g_dir_read_name(child)) != NULL) {
			if (resource_wanted(dname)) {
			    sprintf(tmp, "%s%c%s", subdir, SLASH, dname);
			    resources[i++] = gretl_strdup(tmp);
			}
		    }
		    g_dir_close(child);
		}
	    } else if (resource_wanted(dname)) {
		resources[i++] = gretl_strdup(dname);
	    }
	}
    }

    if (dir != NULL) {
        g_dir_close(dir);
    }

    if (savedir != NULL) {
	/* return to the original working directory */
	gretl_chdir(savedir);
	g_free(savedir);
    }

    *pn_res = n_res;
    g_free(pkg_base);

    if (resources != NULL) {
	qsort(resources, n_res, sizeof *resources, compare_resources);
    }

    return resources;
}

static char **package_datafiles_info (const fnpkg *pkg,
                                      const char *fname,
                                      char resdir[],
                                      int *pn_res)
{
    char **resources = NULL;
    const gchar *dname;
    GDir *dir = NULL;
    int path_len = 0;
    int n_res = 0;

    path_len = strlen(fname) - strlen(pkg->name) - 5;
    resdir[0] = '\0';
    strncat(resdir, fname, path_len);

    dir = gretl_opendir(resdir);
    if (dir == NULL) {
	return NULL;
    }

    /* first get a count of relevant files */
    while ((dname = g_dir_read_name(dir)) != NULL) {
	if (resource_wanted(dname)) {
	    n_res++;
	}
    }

    if (n_res > 0) {
	/* if there are any relevant files, create an array */
	int i = 0;

	resources = strings_array_new(n_res);
	g_dir_rewind(dir);
	while ((dname = g_dir_read_name(dir)) != NULL) {
	    if (resource_wanted(dname)) {
		resources[i++] = gretl_strdup(dname);
	    }
	}
    }

    if (dir != NULL) {
        g_dir_close(dir);
    }

    *pn_res = n_res;

    if (resources != NULL) {
	qsort(resources, n_res, sizeof *resources, compare_resources);
    }

    return resources;
}

/* end apparatus for augmenting "pkg query" to "examples" */

static void plain_print_package_info (const fnpkg *pkg,
                                      const char *fname,
                                      PRN *prn)
{
    char vstr[8];
    int has_examples = 0;
    int i;

    if (g_path_is_absolute(fname) && !strstr(fname, "dltmp.")) {
        pprintf(prn, "File: %s\n", fname);
    }
    if (pkg->name[0] == '\0' || pkg->author == NULL ||
        pkg->minver <= 0 || pkg->descrip == NULL ||
        pkg->version == NULL || pkg->date == NULL ||
        pkg->help == NULL) {
        pprintf(prn, _("Broken package! Basic information is missing\n"));
        return;
    }

    gretl_version_string(vstr, pkg->minver);

    pprintf(prn, "Package: %s %s (%s)\n", pkg->name, pkg->version, pkg->date);
    pprintf(prn, "Author: %s\n", pkg->author);
    if (pkg->email != NULL && *pkg->email != '\0') {
        pprintf(prn, "Email: %s\n", pkg->email);
    }
    pprintf(prn, "Required gretl version: %s\n", vstr);
    pprintf(prn, "Data requirement: %s\n", _(data_needs_string(pkg->dreq)));
    pprintf(prn, "Description: %s\n", gretl_strstrip(pkg->descrip));
    if (pkg->n_depends > 0) {
        pputs(prn, "Dependencies: ");
        for (i=0; i<pkg->n_depends; i++) {
            pprintf(prn, "%s%s", pkg->depends[i],
                    (i < pkg->n_depends-1)? ", " : "\n");
        }
    }
    if (pkg->provider != NULL) {
        pprintf(prn, "Provider: %s\n", pkg->provider);
    }

    has_examples = pkg_has_examples(pkg);

    if (has_examples || pkg_has_datafiles(pkg)) {
        char resdir[MAXLEN];
        char **resources;
        int n_res = 0;

	if (has_examples) {
	    resources = package_resources_info(pkg, fname, resdir, &n_res);
	} else {
	    resources = package_datafiles_info(pkg, fname, resdir, &n_res);
	}
        if (resources != NULL) {
            pprintf(prn, "Resources in %s:\n", resdir);
            for (i=0; i<n_res; i++) {
                pprintf(prn, "%s\n", resources[i]);
            }
            strings_array_free(resources, n_res);
        }
    }

    pputc(prn, '\n');
}

static void real_bundle_package_info (const fnpkg *pkg,
                                      const char *fname,
                                      gretl_bundle *b)
{
    int has_examples = 0;
    int err = 0;

    if (g_path_is_absolute(fname) && !strstr(fname, "dltmp.")) {
        gretl_bundle_set_string(b, "file", fname);
    }
    if (pkg->name[0] == '\0' || pkg->author == NULL ||
        pkg->minver <= 0 || pkg->descrip == NULL ||
        pkg->version == NULL || pkg->date == NULL ||
        pkg->help == NULL) {
        gretl_bundle_set_int(b, "broken", 1);
        return;
    }

    gretl_bundle_set_string(b, "name", pkg->name);
    gretl_bundle_set_string(b, "author", pkg->author);
    /* for backward compatibility, we need a numerical package
       version in the bundle returned by "pkg query"
    */
    if (is_gretl_addon(pkg->name)) {
	gretl_bundle_set_scalar(b, "version",
				gretl_version_number(pkg->version));
    } else {
	gretl_bundle_set_scalar(b, "version", dot_atof(pkg->version));
    }
    gretl_bundle_set_string(b, "date", pkg->date);
    if (pkg->email != NULL && *pkg->email != '\0') {
        gretl_bundle_set_string(b, "email", pkg->email);
    }
    gretl_bundle_set_int(b, "gretl_version", pkg->minver);
    gretl_bundle_set_string(b, "description", gretl_strstrip(pkg->descrip));
    gretl_bundle_set_int(b, "gretl_version", pkg->minver);
    gretl_bundle_set_string(b, "data_requirement", data_needs_string(pkg->dreq));

    if (pkg->n_depends > 0) {
        gretl_array *D = gretl_array_from_strings(pkg->depends,
                                                  pkg->n_depends,
                                                  1, &err);
        gretl_bundle_donate_data(b, "depends", D, GRETL_TYPE_ARRAY, 0);
    }
    if (pkg->provider != NULL) {
        gretl_bundle_set_string(b, "provider", pkg->provider);
    }

    has_examples = pkg_has_examples(pkg);

    if (has_examples || pkg_has_datafiles(pkg)) {
        char resdir[MAXLEN];
        char **resources;
        int n_res = 0;

	if (has_examples) {
	    resources = package_resources_info(pkg, fname, resdir, &n_res);
	} else {
	    resources = package_datafiles_info(pkg, fname, resdir, &n_res);
	}
        if (resources != NULL) {
            gretl_array *R = gretl_array_from_strings(resources, n_res, 0, &err);

            if (!err) {
		gretl_bundle_set_string(b, "resource_dir", resdir);
                gretl_bundle_donate_data(b, "resources", R, GRETL_TYPE_ARRAY, 0);
            }
        }
    }
}

static void reflow_package_help (const char *help, PRN *prn)
{
    char *rem, *prev, line[2048];
    int i, lmax = 80;

    bufgets_init(help);

    while (bufgets(line, sizeof line, help)) {
        if (strlen(line) <= lmax) {
            pputs(prn, line);
        } else {
            rem = line;
            prev = NULL;
            while (strlen(rem) > lmax) {
                for (i=lmax-1; i>0; i--) {
                    if (rem[i] == ' ') {
                        rem[i] = '\0';
                        pprintf(prn, "%s\n", rem);
                        rem = rem + i + 1;
                        break;
                    }
                }
                if (rem == line || rem == prev) {
                    /* let's not get into an infinite loop */
                    break;
                }
                prev = rem;
            }
            if (*rem != '\0') {
                pputs(prn, rem);
            }
        }
    }

    bufgets_finalize(help);
    pputs(prn, "\n\n");
}

/* Simple plain-text help output, or converted markdown if applicable
   -- that is, we're in GUI mode and @pkg has help text in
   markdown. Being in GUI mode is signalled by @pbuf being non-NULL.
*/

static void print_package_help (const fnpkg *pkg,
                                const char *fname,
                                char **pbuf,
                                PRN *prn)
{
    pprintf(prn, "%s %s (%s), %s\n", pkg->name, pkg->version,
            pkg->date, pkg->author);
    pputs(prn, gretl_strstrip(pkg->descrip));
    pputs(prn, "\n\n");

    if (pbuf != NULL && help_text_is_markdown(pkg, 0)) {
        PRN *myprn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);

         /* FIXME should the next line really be needed? */
        pprintf(myprn, "<@bld=\"%s\">\n\n", pkg->name);
        md_to_gretl(pkg->help, myprn);
        *pbuf = gretl_print_steal_buffer(myprn);
        gretl_print_destroy(myprn);
    } else {
         reflow_package_help(pkg->help, prn);
    }
}

static void print_package_code (const fnpkg *pkg,
                                int tabwidth,
                                PRN *prn)
{
    int i;

    pprintf(prn, "# hansl code from package %s %s (%s)\n\n",
            pkg->name, pkg->version, pkg->date);

    if (pkg->priv != NULL) {
        pputs(prn, "# private functions\n\n");
        for (i=0; i<pkg->n_priv; i++) {
            gretl_function_print_code(pkg->priv[i], tabwidth, prn);
            pputc(prn, '\n');
        }
    }

    if (pkg->pub != NULL) {
        pputs(prn, "# public functions\n\n");
        for (i=0; i<pkg->n_pub; i++) {
            gretl_function_print_code(pkg->pub[i], tabwidth, prn);
            pputc(prn, '\n');
        }
    }
}

/* allocate a fnpkg structure and read from XML file into it */

static fnpkg *
real_read_package (xmlDocPtr doc, xmlNodePtr node,
                   const char *fname, int get_funcs,
                   int *err)
{
    xmlNodePtr cur;
    fnpkg *pkg;
    char *tmp = NULL;

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

    if (gretl_xml_get_prop_as_string(node, "model-requirement", &tmp)) {
        pkg->modelreq = gretl_command_number(tmp);
        free(tmp);
    }

    if (gretl_xml_get_prop_as_string(node, "minver", &tmp)) {
        pkg->minver = gretl_version_number(tmp);
        free(tmp);
    }

    pkg->uses_subdir = gretl_xml_get_prop_as_bool(node, "lives-in-subdir");
    pkg->data_access = gretl_xml_get_prop_as_bool(node, "wants-data-access");

    cur = node->xmlChildrenNode;

    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "author")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->author);
            gretl_xml_get_prop_as_string(cur, "email", &pkg->email);
        } else if (!xmlStrcmp(cur->name, (XUC) "version")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->version);
        } else if (!xmlStrcmp(cur->name, (XUC) "date")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->date);
        } else if (!xmlStrcmp(cur->name, (XUC) "description")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->descrip);
        } else if (!xmlStrcmp(cur->name, (XUC) "help")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->help);
            gretl_xml_get_prop_as_string(cur, "filename", &pkg->help_fname);
        } else if (!xmlStrcmp(cur->name, (XUC) "gui-help")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->gui_help);
            gretl_xml_get_prop_as_string(cur, "filename", &pkg->gui_help_fname);
        } else if (!xmlStrcmp(cur->name, (XUC) "R-depends")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->Rdeps);
        } else if (!xmlStrcmp(cur->name, (XUC) "sample-script")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->sample);
            gretl_xml_get_prop_as_string(cur, "filename", &pkg->sample_fname);
        } else if (!xmlStrcmp(cur->name, (XUC) "tags")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->tags);
        } else if (!xmlStrcmp(cur->name, (XUC) "label")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->label);
        } else if (!xmlStrcmp(cur->name, (XUC) "menu-attachment")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->mpath);
        } else if (!xmlStrcmp(cur->name, (XUC) "provider")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &pkg->provider);
        } else if (!xmlStrcmp(cur->name, (XUC) "data-files")) {
            pkg->datafiles =
                gretl_xml_get_strings_array(cur, doc, &pkg->n_files,
                                            0, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "depends")) {
            pkg->depends =
                gretl_xml_get_strings_array(cur, doc, &pkg->n_depends,
                                            0, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "translation")) {
            pkg->trans = read_translation_element(cur, doc);
        }

        cur = cur->next;
    }

    if (get_funcs) {
        cur = node->xmlChildrenNode;
        while (cur != NULL && !*err) {
            if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
                *err = read_ufunc_from_xml(cur, doc, pkg);
            }
            cur = cur->next;
        }
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
#ifdef HAVE_MPI
    int get_caller = 0;
#endif
    int err = 0;

#if PKG_DEBUG
    fprintf(stderr, "read_session_functions_file: starting on '%s'\n", fname);
#endif

    err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (err) {
        return err;
    }

#ifdef HAVE_MPI
    if (gretl_mpi_rank() >= 0) {
        get_caller = 1;
    }
#endif

    /* get any function packages from this file */
    cur = node->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
            pkg = real_read_package(doc, cur, fname, 1, &err);
            if (!err) {
                err = real_load_package(pkg, NULL, NULL, 0);
            }
        }
#ifdef HAVE_MPI
        else if (get_caller && !xmlStrcmp(cur->name, (XUC) "caller")) {
            /* are these functions being loaded from within a
               function that's calling mpi?
            */
            char *caller = gretl_xml_get_string(cur, doc);

            if (caller != NULL) {
                strcpy(mpi_caller, caller);
                free(caller);
            }
            get_caller = 0;
        }
#endif
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
    }

#if PKG_DEBUG
    fprintf(stderr, "read_session_functions_file: returning %d\n", err);
#endif

    return err;
}

/* Parse an XML function package file and return an allocated
   package struct with the functions attached. Note that this
   function does not actually "load" the package (making it
   available to users) and does not check dependencies.
*/

static fnpkg *read_package_file (const char *fname,
                                 int get_funcs,
                                 int *err)
{
    fnpkg *pkg = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;

#if PKG_DEBUG
    fprintf(stderr, "read_package_file: got '%s'\n", fname);
#endif

    node = gretl_xml_get_gfn(fname, &doc, err);
    if (!*err) {
        pkg = real_read_package(doc, node, fname, get_funcs, err);
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

#if PKG_DEBUG
    fprintf(stderr, "read_function_package: err = %d\n", *err);
#endif

    return pkg;
}

static char **read_package_strings_array (const char *fname,
                                          const char *tag,
                                          int *ns)
{
    char **S = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    int err;

    node = gretl_xml_get_gfn(fname, &doc, &err);

    if (node != NULL) {
        cur = node->xmlChildrenNode;
        while (cur != NULL) {
            if (!xmlStrcmp(cur->name, (XUC) tag)) {
                S = gretl_xml_get_strings_array(cur, doc, ns, 0, &err);
                break;
            }
            cur = cur->next;
        }
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    return S;
}

char **package_peek_dependencies (const char *fname, int *ndeps)
{
    return read_package_strings_array(fname, "depends", ndeps);
}

/**
 * function_package_is_loaded:
 * @fname: full path to gfn file.
 * @version: location to receive version info, or NULL.
 *
 * Returns: 1 if the function package with filename @fname is
 * loaded in memory, otherwise 0.
 */

int function_package_is_loaded (const char *fname,
                                const char **version)
{
    return (get_loaded_pkg_by_filename(fname, version) != NULL);
}

/**
 * gfn_is_loaded:
 * @gfnname: basename of gfn file, including suffix.
 *
 * Returns: 1 if the function package with basename @gfnname is
 * loaded in memory, otherwise 0.
 */

int gfn_is_loaded (const char *gfnname)
{
    int ret = 0;

    if (strchr(gfnname, SLASH) == NULL && has_suffix(gfnname, ".gfn")) {
        /* plain basename of gfn file */
        gchar *test = g_strdup(gfnname);
        gchar *p = strrchr(test, '.');

        *p = '\0';
        if (get_function_package_by_name(test)) {
            ret = 1;
        }
        g_free(test);
    }

    return ret;
}

static int not_mpi_duplicate (void)
{
#ifdef HAVE_MPI
    return gretl_mpi_rank() < 1;
#else
    return 0;
#endif
}

static fnpkg *check_for_loaded (const char *fname, gretlopt opt)
{
    fnpkg *pkg = get_loaded_pkg_by_filename(fname, NULL);

    if (pkg != NULL) {
        if (opt & OPT_F) {
            /* force re-reading from file */
            real_function_package_unload(pkg, 1);
            pkg = NULL;
        }
    }

    return pkg;
}

/**
 * load_function_package:
 * @fname: full path to gfn file.
 * @opt: may include OPT_F to force loading even when
 * the package is already loaded.
 * @prn: gretl printer.
 * @ppkg: optional pointer to receive fnpkg pointer.
 *
 * Loads the function package located by @fname into
 * memory, if possible. Supports gretl's "include" command
 * for gfn files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

static int load_function_package (const char *fname,
                                  gretlopt opt,
                                  GArray *pstack,
                                  PRN *prn,
                                  fnpkg **ppkg,
                                  int level)
{
    fnpkg *pkg;
    int err = 0;

    pkg = check_for_loaded(fname, opt);
    if (pkg != NULL) {
        return 0;
    }

    pkg = read_package_file(fname, 1, &err);

    if (!err) {
        if (ppkg != NULL) {
            *ppkg = pkg;
        }
        if (pkg->Rdeps != NULL) {
            err = check_R_depends(pkg->name, pkg->Rdeps, prn);
        }
    }

    if (!err) {
        /* Let's double-check that we don't have a
           colliding package (it would have to be
           with a different filename).
        */
        fnpkg *oldpkg;

        oldpkg = get_function_package_by_name(pkg->name);
        if (oldpkg != NULL) {
            real_function_package_unload(oldpkg, 1);
        }
        err = real_load_package(pkg, pstack, prn, level);
    }

    if (err) {
	gchar *tmp = g_strdup_printf(_("Error loading %s\n"), fname);

	gretl_errmsg_prepend(tmp, err);
	g_free(tmp);
    } else if (pkg != NULL && prn != NULL && not_mpi_duplicate()) {
        pprintf(prn, "%s %s, %s (%s)\n", pkg->name, pkg->version,
                pkg->date, pkg->author);
    }

    return err;
}

/**
 * include_gfn:
 * @fname: full path to gfn file.
 * @opt: may include OPT_F to force loading even when
 * the package is already loaded.
 * @prn: gretl printer.
 *
 * Loads the function package located by @fname into
 * memory, if possible. Supports gretl's "include" command
 * for gfn files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int include_gfn (const char *fname, gretlopt opt, PRN *prn)
{
    GArray *pstack;
    gchar *p, *pkgname;
    int err;

    pkgname = g_path_get_basename(fname);
    p = strrchr(pkgname, '.');
    if (p != NULL) {
        *p = '\0';
    }

    pstack = g_array_new(FALSE, FALSE, sizeof(char *));
    g_array_append_val(pstack, pkgname);
    err = load_function_package(fname, opt, pstack, prn, NULL, 0);
    g_array_free(pstack, TRUE);
    g_free(pkgname);

    return err;
}

/**
 * get_function_package_by_filename:
 * @fname: gfn filename.
 * @err: location to receive error code.
 *
 * If the package whose filename is @fname is already loaded,
 * returns the package pointer, otherwise attempts to load the
 * package from file. Similar to include_gfn() but returns
 * the package pointer rather than just a status code.
 *
 * Returns: package-pointer on success, NULL on error.
 */

fnpkg *get_function_package_by_filename (const char *fname, int *err)
{
    fnpkg *pkg = NULL;
    int i, myerr = 0;

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(fname, pkgs[i]->fname)) {
            pkg = pkgs[i];
            break;
        }
    }

    if (pkg == NULL) {
        PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);

        myerr = include_gfn(fname, OPT_NONE, prn);
        if (!myerr) {
            for (i=0; i<n_pkgs; i++) {
                if (!strcmp(fname, pkgs[i]->fname)) {
                    pkg = pkgs[i];
                    break;
                }
            }
        } else if (prn != NULL) {
            const char *buf = gretl_print_get_buffer(prn);

            gretl_errmsg_set(buf);
        }
        gretl_print_destroy(prn);
    }

    if (err != NULL) {
        *err = myerr;
    }

    return pkg;
}

/**
 * get_function_package_by_name:
 * @pkgname: name of function package.
 *
 * Returns: pointer to function package if a package named
 * @pkgname is already in memory, otherwise NULL.
 */

fnpkg *get_function_package_by_name (const char *pkgname)
{
    int i;

    if (has_suffix(pkgname, ".gfn") || has_suffix(pkgname, ".zip")) {
        /* just in case: strip off the extension */
        gchar *tmp = g_strdup(pkgname);
        gchar *p = strrchr(tmp, '.');
        fnpkg *ret = NULL;

        *p = '\0';
        for (i=0; i<n_pkgs; i++) {
            if (!strcmp(tmp, pkgs[i]->name)) {
                ret = pkgs[i];
                break;
            }
        }
        g_free(tmp);
        return ret;
    } else {
        for (i=0; i<n_pkgs; i++) {
            if (!strcmp(pkgname, pkgs[i]->name)) {
                return pkgs[i];
            }
        }
        return NULL;
    }
}

/**
 * get_function_package_path_by_name:
 * @pkgname: name of function package.
 *
 * Returns: the full path to the named function package if it
 * is already in memory, otherwise NULL.
 */

const char *get_function_package_path_by_name (const char *pkgname)
{
    int i;

    for (i=0; i<n_pkgs; i++) {
        if (!strcmp(pkgname, pkgs[i]->name)) {
            return pkgs[i]->fname;
        }
    }

    return NULL;
}

static int gfn_version_fail (const char *fname, PRN *prn)
{
    FILE *fp = gretl_fopen(fname, "r");
    int err = 0;

    if (fp != NULL) {
        char line[128];

        if (fgets(line, sizeof line, fp) &&
            strchr(line, '<') == NULL &&
            strstr(line, "requires gretl") != NULL) {
            gretl_errmsg_set(gretl_strstrip(line));
            err = 1;
        }
        fclose(fp);
    }

    return err;
}

static char *sample_script_from_xml (const char *pkgname, int *err)
{
    char *ret = NULL;
    gchar *gfnname = NULL;
    char fullname[MAXLEN];
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;

    gfnname = g_strdup_printf("%s.gfn", pkgname);
    *err = get_full_filename(gfnname, fullname, OPT_I);

    if (!*err) {
        node = gretl_xml_get_gfn(fullname, &doc, err);
    }

    if (!*err) {
	xmlNodePtr cur = node->xmlChildrenNode;

	while (cur != NULL) {
            if (!xmlStrcmp(cur->name, (XUC) "sample-script")) {
                gretl_xml_node_get_trimmed_string(cur, doc, &ret);
                break;
            }
            cur = cur->next;
        }
    }

    g_free(gfnname);
    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return ret;
}

int grab_package_sample (const char *pkgname, char **pscript)
{
    int err = 0;

    if (pkgname == NULL) {
	err = E_DATA;
    } else {
	fnpkg *pkg = get_loaded_pkg_by_name(pkgname);
	char *ss = NULL;

	if (pkg != NULL) {
	    ss = gretl_strdup(pkg->sample);
	} else {
	    ss = sample_script_from_xml(pkgname, &err);
	}
	*pscript = ss;
	if (ss == NULL) {
	    gretl_errmsg_sprintf(_("Couldn't find sample script for '%s'"), pkgname);
	    err = E_DATA;
	}
    }

    return err;
}

/* Retrieve summary info, sample script, or code listing for a
   function package, identified by its filename.  This is called
   (indirectly) from the GUI -- see below for the actual callbacks.
*/

static int real_print_gfn_data (const char *fname,
                                char **pbuf, PRN *prn,
                                int tabwidth, int task,
                                gretl_bundle *b)
{
    fnpkg *pkg = NULL;
    int free_pkg = 0;
    int err = 0;

#if PKG_DEBUG
    fprintf(stderr, "real_print_gfn_data: fname='%s', task %d\n", fname, task);
#endif

    if (task == FUNCS_INFO && prn != NULL && strstr(fname, "dltmp.")) {
        if (gfn_version_fail(fname, prn)) {
            return 1;
        }
    }

    pkg = get_loaded_pkg_by_filename(fname, NULL);

#if PKG_DEBUG
    fprintf(stderr, "real_print_gfn_data: pkg=%p\n", (void *) pkg);
#endif

    if (pkg == NULL) {
        /* package is not loaded, read it now */
        int get_funcs = 1;

        if (task == FUNCS_HELP || task == FUNCS_INFO ||
            task == FUNCS_SAMPLE || task == FUNCS_QUERY) {
            get_funcs = 0;
        }
        pkg = read_package_file(fname, get_funcs, &err);
        free_pkg = 1;
    }

    if (!err) {
        if (task == FUNCS_HELP) {
            print_package_help(pkg, fname, pbuf, prn);
        } else if (task == FUNCS_INFO) {
            print_package_info(pkg, fname, prn);
        } else if (task == FUNCS_QUERY) {
            if (b != NULL) {
                real_bundle_package_info(pkg, fname, b);
            } else {
                plain_print_package_info(pkg, fname, prn);
            }
        } else if (task == FUNCS_SAMPLE) {
            pputs(prn, pkg->sample);
        } else if (task == FUNCS_CODE) {
            /* called by print_function_package_code() */
            print_package_code(pkg, tabwidth, prn);
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

/* callback used in the GUI function package browser */

int print_function_package_info (const char *fname, int gui_mode,
                                 PRN *prn)
{
    int mode = gui_mode ? FUNCS_INFO : FUNCS_QUERY;

    return real_print_gfn_data(fname, NULL, prn, 0, mode, NULL);
}

/* callback used by "pkg" command with query action + --quiet */

int bundle_function_package_info (const char *fname, gretl_bundle *b)
{
    return real_print_gfn_data(fname, NULL, NULL, 0, FUNCS_QUERY, b);
}

/* callback used in the GUI function package browser (gui/datafiles.c) */

int print_function_package_code (const char *fname, int tabwidth,
                                 PRN *prn)
{
    return real_print_gfn_data(fname, NULL, prn, tabwidth, FUNCS_CODE, NULL);
}

/* callback used in the GUI function package browser */

int print_function_package_sample (const char *fname, int tabwidth,
                                   PRN *prn)
{
    return real_print_gfn_data(fname, NULL, prn, tabwidth, FUNCS_SAMPLE, NULL);
}

/* callback used via command line */

int print_function_package_help (const char *fname, char **pbuf, PRN *prn)
{
    return real_print_gfn_data(fname, pbuf, prn, 0, FUNCS_HELP, NULL);
}

static void maybe_fix_broken_date (char **pdate)
{
    char *date = *pdate;

    if (strlen(date) < 10 && date[4] == '-') {
        int y, m, d;

        if (sscanf(date, "%d-%d-%d", &y, &m, &d) == 3) {
            char tmp[12];

            sprintf(tmp, "%d-%02d-%02d", y, m, d);
            free(*pdate);
            *pdate = gretl_strdup(tmp);
        }
    }
}

static int get_pdfdoc_status (xmlNodePtr sub,
                              xmlDocPtr doc)
{
    char *tmp = NULL;
    int ret = 0;

    gretl_xml_node_get_trimmed_string(sub, doc, &tmp);
    if (tmp != NULL) {
        if (!strncmp(tmp, "pdfdoc", 6) || is_pdf_ref(tmp)) {
            ret = 1;
        }
        free(tmp);
    }

    return ret;
}

/* Read the header from a function package file -- this is used when
   displaying the available packages in the GUI (see gui/datafiles.c).
   We write the description into *pdesc, a string representation of
   version number into *pver, date into *pdate, and author's name into
   *pauthor.

   The *pdfdoc pointer can be used to collect an int value indicating
   whether the package contains PDF documentation (1/0), or it can be
   set to NULL. None of the other arguments can be NULL.
*/

int get_function_file_header (const char *fname, char **pdesc,
                              char **pver, char **pdate,
                              char **pauthor, int *pdfdoc)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    int docdone = 0;
    int ndone = 0;
    int err = 0;

    if (pdesc == NULL || pver == NULL || pdate == NULL || pauthor == NULL) {
        fprintf(stderr, "get_function_file_header: missing parameter\n");
        return E_DATA;
    }

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return err;
    }

    if (pdfdoc == NULL) {
        /* not wanted, so count it as "done" already */
        docdone = 1;
    } else {
        *pdfdoc = 0;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "description")) {
            gretl_xml_node_get_trimmed_string(cur, doc, pdesc);
            ndone++;
        } else if (!xmlStrcmp(cur->name, (XUC) "version")) {
            gretl_xml_node_get_trimmed_string(cur, doc, pver);
            ndone++;
        } else if (!xmlStrcmp(cur->name, (XUC) "date")) {
            gretl_xml_node_get_trimmed_string(cur, doc, pdate);
            if (*pdate != NULL) {
                maybe_fix_broken_date(pdate);
            }
            ndone++;
        } else if (!xmlStrcmp(cur->name, (XUC) "author")) {
            gretl_xml_node_get_trimmed_string(cur, doc, pauthor);
            ndone++;
        } else if (!docdone && !xmlStrcmp(cur->name, (XUC) "help")) {
            *pdfdoc = get_pdfdoc_status(cur, doc);
            docdone = 1;
        }
        if (docdone && ndone == 4) {
            break;
        }
        cur = cur->next;
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    if (*pdesc == NULL) {
        *pdesc = gretl_strdup(_("No description available"));
    }
    if (*pver == NULL) {
        *pver = gretl_strdup("unknown");
    }
    if (*pdate == NULL) {
        *pdate = gretl_strdup("unknown");
    }
    if (*pauthor == NULL) {
        *pauthor = gretl_strdup("unknown");
    }

    if (*pdesc == NULL || *pver == NULL ||
        *pdate == NULL || *pauthor == NULL) {
        err = E_ALLOC;
    }

    return err;
}

int package_has_menu_attachment (const char *fname,
                                 char **pkgname,
                                 char **attach,
                                 char **label)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    char *tmp = NULL;
    int skipit = 0;
    int got_label = 0;
    int got_attach = 0;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return 0;
    }

    if (pkgname != NULL) {
        gretl_xml_get_prop_as_string(node, "name", &tmp);
        if (tmp != NULL && !strcmp(tmp, "ridge")) {
            /* don't give "ridge" a menu attachment */
            free(tmp);
            skipit = 1;
        } else {
            *pkgname = tmp;
        }
    }

    if (skipit) {
        xmlFreeDoc(doc);
        return 0;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "menu-attachment")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            if (tmp != NULL) {
                got_attach = 1;
                if (attach != NULL) {
                    *attach = tmp;
                } else {
                    free(tmp);
                }
            }
        } else if (!xmlStrcmp(cur->name, (XUC) "label")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            if (tmp != NULL) {
                got_label = 1;
                if (label != NULL) {
                    *label = tmp;
                } else {
                    free(tmp);
                }
            }
        } else if (!xmlStrcmp(cur->name, (XUC) "help")) {
            /* we've overshot */
            break;
        }
        if (got_attach && got_label) {
            break;
        }
        cur = cur->next;
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    return got_attach && got_label;
}

int package_needs_zipping (const char *fname,
                           int *pdfdoc,
                           char ***datafiles,
                           int *n_files)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;
    char *tmp = NULL;
    int retmax = 1;
    int ret = 0;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);
    if (err) {
        return 0;
    }

    if (datafiles != NULL) {
        retmax = 2;
    }

    cur = node->xmlChildrenNode;
    while (cur != NULL) {
        if (!xmlStrcmp(cur->name, (XUC) "help")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            if (tmp != NULL && !strncmp(tmp, "pdfdoc:", 7)) {
                if (pdfdoc != NULL) {
                    *pdfdoc = 1;
                }
                ret++;
            }
            free(tmp);
        } else if (!xmlStrcmp(cur->name, (XUC) "data-files")) {
            if (datafiles != NULL) {
                *datafiles = gretl_xml_get_strings_array(cur, doc, n_files,
                                                         0, &err);
            }
            ret++;
        } else if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
            /* we've overshot */
            break;
        }
        if (ret == retmax) {
            break;
        }
        cur = cur->next;
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    return ret;
}

double gfn_file_get_version (const char *fname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr node;
    double version = NADBL;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);

    if (!err) {
        xmlNodePtr cur = node->xmlChildrenNode;

        while (cur != NULL) {
            if (!xmlStrcmp(cur->name, (XUC) "version")) {
                gretl_xml_node_get_double(cur, doc, &version);
                break;
            }
            cur = cur->next;
        }
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    return version;
}

const char *try_for_label_translation (const char *label,
                                       const char *fname,
                                       const char *lang)
{
    const char *ret = label;
    xmlDocPtr doc = NULL;
    xmlNodePtr node;
    int err = 0;

    node = gretl_xml_get_gfn(fname, &doc, &err);

    if (!err) {
        xmlNodePtr cur = node->xmlChildrenNode;
        xmlNodePtr msg;
        char *trlang = NULL;
        char *tr = NULL;

        while (cur != NULL) {
            if (!xmlStrcmp(cur->name, (XUC) "translation")) {
                if (gretl_xml_get_prop_as_string(cur, "lang", &trlang)) {
                    if (!strncmp(lang, trlang, 2)) {
                        msg = cur->xmlChildrenNode;
                        while (msg != NULL) {
                            if (gretl_xml_node_get_string(msg, doc, &tr)) {
                                if (!strcmp(label, tr)) {
                                    ret = tr;
                                    break;
                                }
                            }
                            msg = msg->next;
                        }
                    } else {
                        break;
                    }
                }
                break;
            }
            cur = cur->next;
        }
    }

    if (doc != NULL) {
        xmlFreeDoc(doc);
    }

    return ret;
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

void set_current_function_package (fnpkg *pkg)
{
    current_pkg = pkg;
}

static int skip_private_func (ufunc *ufun)
{
    int skip = 0;

    if (function_is_private(ufun)) {
        /* skip it, unless we're "authorized" to edit it,
           i.e. coming from the gui package editor */
        skip = 1;
        if (ufun->pkg != NULL && ufun->pkg == current_pkg) {
            skip = 0;
        }
    }

    return skip;
}

static int check_function_name (const char *name, ufunc **pfun,
                                const DATASET *dset, PRN *prn)
{
    int i, err = 0;

#if FN_DEBUG
    fprintf(stderr, "check_function_name: '%s'\n", name);
#endif

    if (!isalpha((unsigned char) *name)) {
        gretl_errmsg_set(_("Function names must start with a letter"));
        err = 1;
    } else if (gretl_reserved_word(name)) {
        err = 1;
    } else if (function_lookup(name)) {
        gretl_errmsg_sprintf(_("'%s' is the name of a built-in function"), name);
        err = 1;
    } else if (gretl_is_user_var(name)) {
        gretl_errmsg_sprintf(_("'%s' is the name of a variable"), name);
        err = 1;
    } else if (gretl_is_series(name, dset)) {
        gretl_errmsg_sprintf(_("'%s' is the name of a variable"), name);
        err = 1;
    } else {
        /* @name is OK, now check for existing function of the same name */
        for (i=0; i<n_ufuns; i++) {
            if (!skip_private_func(ufuns[i]) && !strcmp(name, ufuns[i]->name)) {
#if FN_DEBUG
                fprintf(stderr, "'%s': found an existing function of this name\n", name);
#endif
                if (ufuns[i]->pkg != NULL && ufuns[i]->pkg != current_pkg) {
                    /* don't allow overwriting */
                    gretl_errmsg_sprintf(_("The function %s is already defined "
                                         "by package '%s'"), name,
                                         ufuns[i]->pkg->name);
                    err = E_DATA;
                } else if (pfun != NULL) {
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
        gretl_errmsg_sprintf(_("%s: function is in use"), fname);
        err = 1;
    } else if (fun->pkg != NULL) {
        gretl_errmsg_sprintf(_("%s: function belongs to package"), fname);
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

static int check_parm_min_max (fn_param *p, const char *name,
                               int *nvals)
{
    int err = 0;

    if (p->type != GRETL_TYPE_DOUBLE && na(p->deflt)) {
        gretl_errmsg_sprintf(_("%s: parameters of type %s cannot have a "
                             "default value of NA"), name,
                             gretl_type_get_name(p->type));
        err = E_DATA;
    }

    if (!err && !na(p->min) && !na(p->max)) {
        if (p->min > p->max) {
            gretl_errmsg_sprintf("%s: min > max", name);
            err = E_DATA;
        } else if (!na(p->step) && p->step > p->max - p->min) {
            gretl_errmsg_sprintf("%s: step > max - min", name);
            err = E_DATA;
        }
        if (!err) {
            *nvals = (int) p->max - (int) p->min + 1;
        }
    }

    if (!err && !na(p->deflt) && !default_unset(p)) {
        if (!na(p->min) && p->deflt < p->min) {
            gretl_errmsg_sprintf(_("%s: default value out of bounds"), name);
            err = E_DATA;
        } else if (!na(p->max) && p->deflt > p->max) {
            gretl_errmsg_sprintf(_("%s: default value out of bounds"), name);
            err = E_DATA;
        }
    }

    return err;
}

static int colon_count (const char *p)
{
    int n = 0;

    while (*p && *p != ']') {
        if (*p == ':') {
            n++;
        }
        p++;
    }

    return n;
}

#define VALDEBUG 0

static int read_min_max_deflt (char **ps, fn_param *param,
                               const char *name, int *nvals)
{
    char *p = *ps;
    int err = 0;

    gretl_push_c_numeric_locale();

    errno = 0;

    if (!strcmp(p, "[null]")) {
        gretl_errmsg_set(_("'null' is not a valid default for scalars"));
        err = E_TYPES;
    } else if (!strcmp(p, "[auto]")) {
        /* FIXME */
        err = E_TYPES;
    } else if (param->type == GRETL_TYPE_BOOL) {
        if (!strncmp(p, "[TRUE]", 6)) {
            param->deflt = 1;
        } else if (!strncmp(p, "[FALSE]", 7)) {
            param->deflt = 0;
        } else if (sscanf(p, "[%lf]", &param->deflt) != 1) {
            err = E_PARSE;
        }
    } else if (!strncmp(p, "[$xlist]", 8)) {
        param->deflt = INT_USE_XLIST;
    } else if (!strncmp(p, "[$mylist]", 9)) {
        param->deflt = INT_USE_MYLIST;
    } else {
        int i, len, nf = colon_count(p) + 1;
        char *test, valstr[32];
        char *q = p + 1;
        double x[4] = {NADBL, NADBL, UNSET_VALUE, NADBL};

        for (i=0; i<nf && !err; i++) {
            len = strcspn(q, ":]");
            if (len > 31) {
                err = E_PARSE;
            } else if (len > 0) {
                *valstr = '\0';
                strncat(valstr, q, len);
#if VALDEBUG > 1
                fprintf(stderr, "valstr(%d) = '%s'\n", i, valstr);
#endif
                if (!strcmp(valstr, "NA")) {
                    x[i] = NADBL;
                } else {
                    x[i] = strtod(valstr, &test);
                    if (*test != '\0' || errno) {
                        err = E_PARSE;
                    }
                }
            }
            q += len + 1;
        }

        if (!err) {
            if (nf == 1) {
                /* we take a single value with no colons as
                   indicating a default */
                param->deflt = x[0];
            } else {
                param->min = x[0];
                param->max = x[1];
                param->deflt = x[2];
                param->step = x[3];
            }
#if VALDEBUG
            fprintf(stderr, "min %g, max %g, deflt %g, step %g\n",
                    param->min, param->max, param->deflt, param->step);
#endif
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
    int err = 0;

#if FNPARSE_DEBUG
    fprintf(stderr, "read_param_option: got '%s'\n", *ps);
#endif

    if (!strncmp(*ps, "[auto]", 6)) {
        set_param_auto(param);
        *ps += 6;
    } else if (!strncmp(*ps, "[null]", 6)) {
        set_param_optional(param);
        *ps += 6;
    } else {
        gretl_errmsg_sprintf(_("got invalid field '%s'"), *ps);
        err = E_PARSE;
    }

    return err;
}

/* get the descriptive string for a function parameter */

static int read_param_descrip (char **ps, fn_param *param)
{
    int offset = 1;
    char *p = *ps + offset;
    int len = 0;
    int err = E_PARSE;

    while (*p) {
        if (*p == '"') {
            /* found closing quote */
            err = 0;
            break;
        }
        len++;
        p++;
    }

    if (!err && len > 0) {
        p = *ps + offset;
        param->descrip = gretl_strndup(p, len);
        if (param->descrip == NULL) {
            err = E_ALLOC;
        } else {
            *ps = p + len + 1;
        }
    }

    return err;
}

/* Get the value labels for a function parameter. The syntactic
   element we're looking at starts with "{".
*/

static int read_param_labels (char **ps, fn_param *param,
                              const char *name, int nvals)
{
    char *p = *ps;
    int offset = 1;
    int len = 0;
    int err = E_PARSE;

    p += offset;

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

        p = *ps + offset;
        tmp = gretl_strndup(p, len);

        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            *ps = p + len + 1;
            param->labels = gretl_string_split_quoted(tmp, &param->nlabels,
                                                      " ,", &err);
            free(tmp);
            if (!err && param->nlabels != nvals) {
                gretl_errmsg_sprintf(_("%s: found %d values but %d value-labels"),
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
        strings_array_free(param->labels, param->nlabels);
        param->labels = NULL;
        param->nlabels = 0;
    }
}

static int parse_function_param (char *s, fn_param *param, int i)
{
    GretlType type = GRETL_TYPE_NONE;
    char tstr[22] = {0};
    char *name;
    int len, nvals = 0;
    int const_flag = 0;
    int err = 0;

    while (isspace(*s)) s++;

#if FNPARSE_DEBUG
    fprintf(stderr, "parse_function_param: s = '%s'\n", s);
#endif

    /* pick up the "const" flag if present */
    if (!strncmp(s, "const ", 6)) {
        const_flag = 1;
        s += 6;
        while (isspace(*s)) s++;
    }

    /* get parameter type, required */
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
            gretl_errmsg_set(_("Expected a type identifier"));
            err = E_PARSE;
        } else {
            type = gretl_type_from_string(tstr);
            if (type == 0) {
                gretl_errmsg_sprintf(_("Unrecognized data type '%s'"), tstr);
                err = E_PARSE;
            } else if (!ok_function_arg_type(type)) {
                gretl_errmsg_sprintf(_("%s: invalid parameter type"), tstr);
                err = E_INVARG;
            }
        }
    }

    if (err) {
        return err;
    }

    /* get the required parameter name */
    while (isspace(*s)) s++;
    len = gretl_namechar_spn(s);
    if (len == 0) {
        gretl_errmsg_sprintf(_("parameter %d: name is missing"), i + 1);
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
    if (const_flag) {
        maybe_set_param_const(param);
    }

    s += len;
    s += strspn(s, " ");

    /* now scan for various optional extras: first we may have
       a [min:max:default] specification (scalars only), or a
       [null|auto] spec to allow no argument for "pointer"-type
       parameters (broadly interpreted)
    */

    if (*s == '[') {
        if (gretl_scalar_type(type)) {
            err = read_min_max_deflt(&s, param, name, &nvals);
        } else if (arg_may_be_optional(type)) {
            err = read_param_option(&s, param);
        } else {
            gretl_errmsg_sprintf(_("'%s': error scanning default value"),
                                 name);
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
        fprintf(stderr, " param[%d] = '%s', ptype = %d (%s)\n",
                i, name, type, gretl_type_get_name(type));
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
        if (isspace(s[i])) {
            s[i] = '\0';
        } else if (s[i] == ')') {
            s[i] = '\0';
            break;
        } else {
            break;
        }
    }
}

static int parse_function_parameters (const char *str,
                                      fn_param **pparams,
                                      int *pnp,
                                      const char *fname,
                                      PRN *prn)
{
    fn_param *params = NULL;
    char *p, *s;
    int i, j, np = 0;
    int err = 0;

    s = gretl_strdup(str);
    if (s == NULL) {
        return E_ALLOC;
    }

    if (!strcmp(s, ")")) {
        /* we got a void function "foo()" */
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

    for (i=0; i<np && !err; i++) {
        for (j=i+1; j<np && !err; j++) {
            if (!strcmp(params[i].name, params[j].name)) {
                gretl_errmsg_sprintf(_("%s: duplicated parameter name '%s'"),
                                     fname, params[i].name);
                err = E_DATA;
            }
        }
    }

    if (err) {
        free_params_array(params, np);
    } else {
        *pparams = params;
        *pnp = np;
    }

    return err;
}

static int get_two_words (const char *s, char *w1, char *w2,
                          int *err)
{
    int nf = 0;

    if (strncmp(s, "function ", 9)) {
        *err = E_PARSE;
    } else {
        int n;

        *w1 = *w2 = '\0';

        s += 9;
        s += strspn(s, " ");

        /* @w1 should be a return type, except in the case
           of "function foo delete"
        */
        n = strcspn(s, " ");

        if (n == 0 || n >= FN_NAMELEN) {
            *err = E_PARSE;
        } else {
            strncat(w1, s, n);
            nf++;
            s += n;
            s += strspn(s, " ");
            n = strcspn(s, " (");
        }

        /* @w2 should generally be the function name */
        if (n > 0 && n < FN_NAMELEN) {
            strncat(w2, s, n);
            nf++;
        } else if (n >= FN_NAMELEN) {
            gretl_errmsg_set(_("Identifier exceeds the maximum of 31 characters"));
            *err = E_PARSE;
        }
    }

    return nf;
}

/**
 * gretl_start_compiling_function:
 * @line: command line.
 * @prn: printing struct for feedback.
 *
 * Responds to a statement of the form "function ...".  In most
 * cases, embarks on compilation of a function, but this
 * also handles the construction "function foo delete".
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_start_compiling_function (const char *line,
                                    const DATASET *dset,
                                    PRN *prn)
{
    ufunc *fun = NULL;
    fn_param *params = NULL;
    int nf, n_params = 0;
    int rettype = 0;
    char s1[FN_NAMELEN];
    char s2[FN_NAMELEN];
    char *p = NULL, *name = NULL;
    int err = 0;

    if (gretl_function_depth() > 0) {
        return E_FNEST;
    }

    nf = get_two_words(line, s1, s2, &err);

    if (err) {
        return err;
    } else if (nf < 2) {
        gretl_errmsg_set(_("A function definition must have a return type and name"));
        return E_PARSE;
    }

    if (!strcmp(s2, "clear") || !strcmp(s2, "delete")) {
        return maybe_delete_function(s1, prn);
    }

    /* If we didn't get a special such as "function foo delete",
       then @s1 should contain the return type and @s2 the
       name.
    */

    name = s2;
    rettype = return_type_from_string(s1, &err);

    if (!err) {
        /* note: this handles a name collision */
        err = check_function_name(name, &fun, dset, prn);
    }

    if (!err) {
        /* now for the args bit */
        p = strchr(line, '(');
        if (p == NULL || strchr(p, ')') == NULL) {
            err = E_PARSE;
        } else {
	    /* skip the left paren and any following spaces */
            p++;
	    p += strspn(p, " \t");
        }
    }

    if (!err) {
        err = parse_function_parameters(p, &params, &n_params,
                                        name, prn);
    }

    if (err) {
        pprintf(prn, "> %s\n", line);
    }

    if (!err && fun == NULL) {
        fun = add_ufunc(name);
        if (fun == NULL) {
            free_params_array(params, n_params);
            err = E_ALLOC;
        }
    }

#if COMP_DEBUG
    fprintf(stderr, "started compiling function %s (err = %d)\n",
            name, err);
#endif

    if (!err) {
        strcpy(fun->name, name);
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

static void python_check (const char *line)
{
    char s1[8], s2[16];

    if (sscanf(line, "%7s %15s", s1, s2) == 2) {
        if (!strcmp(s1, "foreign") && strstr(s2, "ytho")) {
            compiling_python = 1;
        } else if (!strcmp(s1, "end") && !strcmp(s2, "foreign")) {
            compiling_python = 0;
        }
    }
}

static int ufunc_get_structure (ufunc *u)
{
    int ret, recurses = 0;

    ret = statements_get_structure(u->lines, u->n_lines,
				   u->name, &recurses);
    if (recurses) {
#if 0
	fprintf(stderr, "%s calls itself\n", u->name);
#endif
	u->flags |= UFUN_RECURSES;
    }

    return ret;
}

#define NEEDS_IF(c) (c == ELSE || c == ELIF || c == ENDIF)

#define FLOW_CI(c) (c == IF || c == ELSE || c == ELIF || c == ENDIF)

/**
 * gretl_function_append_line:
 * @s: pointer to execution state.
 *
 * Continuation of definition of user-function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_function_append_line (ExecState *state)
{
    static int ifdepth;
    static int last_flow;
    CMD *cmd = state->cmd;
    char *line = state->line;
    char *origline = NULL;
    char *s;
    ufunc *fun = current_fdef;
    int blank = 0;
    int ignore = 0;
    int i, err = 0;

    if (fun == NULL) {
        fprintf(stderr, "gretl_function_append_line: fun is NULL\n");
        return 1;
    }

    blank = string_is_blank(line);

#if FNPARSE_DEBUG
    if (fun->line_idx == 1) {
        fputc('\n', stderr);
    }
    if (blank) {
        fprintf(stderr, "%s: append blank line (idx = %d)\n",
                fun->name, fun->line_idx);
    } else {
        fprintf(stderr, "%s: append line '%s' (idx = %d)\n",
                fun->name, line, fun->line_idx);
    }
#endif

    /* note: avoid losing comment lines */
    origline = gretl_strdup(line);

    if (!blank) {
	s = line + strspn(line, " \t");
	state->line = s;    /* temporary! */
        err = get_command_index(state, FUNC, 0);
	state->line = line; /* reset to original */
    }
    if (blank || err) {
        goto next_step;
    }

    /* carry out some basic structural checks */
    if (cmd->ci == QUIT) {
        gretl_errmsg_sprintf(_("%s: \"quit\" cannot be used in a function"),
                             fun->name);
        err = E_PARSE;
        cmd->ci = 0;
    } else if (cmd->flags & CMD_ENDFUN) {
        if (fun->n_lines == 0) {
            gretl_errmsg_sprintf(_("%s: empty function"), fun->name);
            err = 1;
        }
        set_compiling_off();
    } else if (cmd->ci == SET && changes_state(s + 3)) {
        fun->flags |= UFUN_USES_SET;
    } else if (cmd->ci == FUNC) {
        err = E_FNEST;
    } else if (FLOW_CI(cmd->ci)) {
        if (cmd->ci == IF) {
            fun->flags |= UFUN_HAS_FLOW;
            ifdepth++;
        } else if (ifdepth == 0) {
            gretl_errmsg_sprintf(_("%s: unbalanced if/else/endif"), fun->name);
            err = E_PARSE;
        } else if (cmd->ci == ELSE && last_flow == ELSE) {
            gretl_errmsg_sprintf(_("%s: unbalanced if/else/endif"), fun->name);
            err = E_PARSE;
        } else if (cmd->ci == ENDIF) {
            ifdepth--;
        }
        last_flow = cmd->ci;
    } else if (cmd->ci == LOOP) {
        fun->flags |= UFUN_HAS_FLOW;
    } else if (cmd->ci < 0) {
        ignore = 1;
    }

  next_step:

    if (err) {
        set_compiling_off();
    }

    if (compiling) {
        if (blank) {
            /* just increment the index (placeholder) */
            fun->line_idx += 1;
        } else {
            /* actually add the line */
            i = fun->n_lines;
            err = push_function_line(fun, origline, cmd->ci, 1, NULL);
            if (!err) {
                origline = NULL; /* successfully donated */
                if (ignore) {
                    fun->lines[i].flags |= LINE_IGNORE;
                } else {
                    python_check(line);
                }
            }
        }
    } else {
        /* finished compilation */
        if (!err && ifdepth != 0) {
            gretl_errmsg_sprintf(_("%s: unbalanced if/else/endif"), fun->name);
            err = E_PARSE;
        }
#if COMP_DEBUG
        fprintf(stderr, "finished compiling function %s\n", fun->name);
#endif
        /* reset static var */
        ifdepth = 0;
        last_flow = 0;
    }

    if (!err && !compiling && (fun->flags & UFUN_HAS_FLOW)) {
        ufunc_get_structure(fun);
    }

    free(origline);
    cmd->flags &= ~CMD_ENDFUN;

    if (err) {
        ufunc_unload(fun);
    }

    return err;
}

/* next block: handling function arguments */

/* Given a list of variables supplied as an argument to a function,
   copy the list under the name assigned by the function and
   make the variables referenced in that list accessible within
   the function. We also handle the case where the given arg
   was not a real list, but we can make a list out of it for the
   purposes of the function.
*/

/* Missing or "null" arg -> gives an empty list, rather
   than no list. This is traditional behavior though it's
   (now) inconsistent with what happens with other types
   of argument.
*/

static int add_empty_list (fncall *call, fn_param *fp,
                           DATASET *dset)
{
    int err = 0;

    if (dset == NULL || dset->n == 0) {
        err = E_NODATA;
    } else {
        int tmp[] = {0};

        copy_list_as_arg(fp->name, tmp, &err);
    }

    return err;
}

static void localize_list_members (fncall *call, int *list,
                                   const char *lname,
                                   DATASET *dset)
{
    int i, vi, level = fn_executing + 1;

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (vi > 0 && vi < dset->v) {
            if (!in_gretl_list(call->listvars, vi)) {
                gretl_list_append_term(&call->listvars, vi);
            }
            series_set_stack_level(dset, vi, level);
        }
    }
}

static int localize_list (fncall *call, fn_arg *arg,
                          fn_param *fp, DATASET *dset)
{
    int *list = NULL;
    int err = 0;

    if (arg->type == GRETL_TYPE_LIST) {
        /* actual list arg -> copy to function level */
        list = arg->val.list;
        err = copy_as_arg(fp->name, GRETL_TYPE_LIST, list);
        call->lists = g_list_prepend(call->lists, fp->name);
    } else if (arg->type == GRETL_TYPE_USERIES) {
        /* series arg -> becomes a singleton list */
        int tmp[] = {1, arg->val.idnum};

        list = copy_list_as_arg(fp->name, tmp, &err);
        nullify_upname(arg);
    } else {
        /* "can't happen" */
        err = E_DATA;
    }

    if (!err && list == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
        localize_list_members(call, list, fp->name, dset);
    }

#if ARGS_DEBUG
    fprintf(stderr, "localize_list (%s): returning %d\n",
            fp->name, err);
#endif

    return err;
}

/* For the series with index number @ID, find the name of
   the list parameter of the current function call under
   which this series was made accessible to the function.
   We call this only if we know already that series @ID
   was in fact passed to the function via a list.
*/

const char *series_get_list_parent (int ID)
{
    fncall *call = current_function_call();
    const char *ret = NULL;

    if (call != NULL && call->lists != NULL) {
        GList *L = call->lists;
        int *list;

        while (L != NULL) {
            list = get_list_by_name(L->data);
            if (list != NULL && in_gretl_list(list, ID)) {
                ret = L->data;
                break;
            }
            L = L->next;
        }
    }

    return ret;
}

static int localize_bundled_lists (fncall *call, fn_arg *arg,
                                   fn_param *fp, DATASET *dset)
{
    gretl_bundle *b = arg->val.b;
    GList *ll = gretl_bundle_get_lists(b);
    int err = 0;

    if (ll != NULL) {
        GList *lli = g_list_first(ll);

        while (lli != NULL) {
            localize_list_members(call, lli->data, NULL, dset);
            lli = g_list_next(lli);
        }
        g_list_free(ll);
    }

    return err;
}

static void *arg_get_data (fn_arg *arg, GretlType type)
{
    void *data = NULL;

    if (type == 0) {
        type = arg->type;
    }

    if (type == GRETL_TYPE_MATRIX) {
        data = arg->val.m;
    } else if (type == GRETL_TYPE_BUNDLE) {
        data = arg->val.b;
    } else if (type == GRETL_TYPE_STRING) {
        data = arg->val.str;
    } else if (gretl_is_array_type(type)) {
        data = arg->val.a;
    }

    return data;
}

/* handle the case where the GUI passed an anonymous
   object in "pointerized" form */

static int localize_object_as_shell (fn_arg *arg,
                                     fn_param *fp)
{
    GretlType type = gretl_type_get_plain_type(arg->type);
    void *data = arg_get_data(arg, type);

    return arg_add_as_shell(fp->name, type, data);
}

static int localize_const_object (fncall *call, int i, fn_param *fp)
{
    fn_arg *arg = &call->args[i];
    void *data = arg_get_data(arg, 0);
    int err = 0;

    if (data == NULL) {
        err = E_DATA;
    } else if (arg->uvar == NULL) {
        /* the const argument is an anonymous object */
        err = arg_add_as_shell(fp->name, arg->type, data);
    } else {
        /* a named object: in the simplest case we just
           adjust the level and name of the object for
           the duration of the function call, but note
           that we can do this only once for any given
           object.
        */
        user_var *thisvar = arg->uvar;
        int j, done = 0;

        for (j=0; j<i; j++) {
            if (call->args[j].uvar == thisvar) {
                /* we've already localized this one! */
                err = arg_add_as_shell(fp->name, arg->type, data);
                done = 1;
                break;
            }
        }
        if (!done) {
            user_var_adjust_level(thisvar, 1);
            user_var_set_name(thisvar, fp->name);
            arg->shifted = 1;
        }
    }

    return err;
}

static int localize_series_ref (fncall *call, fn_arg *arg,
                                fn_param *fp, DATASET *dset)
{
    int level = fn_executing + 1;
    int v = arg->val.idnum;

    series_set_stack_level(dset, v, level);
    strcpy(dset->varname[v], fp->name);

    if (!in_gretl_list(call->ptrvars, v)) {
        gretl_list_append_term(&call->ptrvars, v);
    }

    return 0;
}

static int argval_get_int (fn_param *param, double x, int *err)
{
    int ret = gretl_int_from_double(x, err);

    if (*err) {
        gretl_errmsg_sprintf(_("%s: expected an integer but found %g"),
                             param->name, x);
    }

    return ret;
}

static int real_add_scalar_arg (fn_param *param, double x)
{
    int err = 0;

    if (na(x)) {
        /* always allow NA for scalar args (?) */
        err = copy_as_arg(param->name, GRETL_TYPE_DOUBLE, &x);
    } else {
        if (param->type == GRETL_TYPE_BOOL) {
            if (x != 0.0) {
                x = 1.0;
            }
        } else if (param->type == GRETL_TYPE_INT ||
                   param->type == GRETL_TYPE_OBS) {
            x = argval_get_int(param, x, &err);
        }

        if (!err) {
            if ((!na(param->min) && x < param->min) ||
                (!na(param->max) && x > param->max)) {
                gretl_errmsg_sprintf(_("%s: argument value %g is out of bounds"),
                                     param->name, x);
                err = E_DATA;
            } else {
                err = copy_as_arg(param->name, GRETL_TYPE_DOUBLE, &x);
            }
        }
    }

    return err;
}

/* Scalar function arguments only: if the arg is not supplied, use the
   default that is found in the function specification, if any.
*/

static int add_scalar_arg_default (fn_param *param)
{
    if (default_unset(param)) {
        /* should be impossible here, but... */
        return E_DATA;
    } else {
        return real_add_scalar_arg(param, param->deflt);
    }
}

/* Handle the case of a scalar parameter for which a 1 x 1 matrix
   was given as argument.
*/

static int do_matrix_scalar_cast (fn_arg *arg, fn_param *param)
{
    gretl_matrix *m = arg->val.m;

    return real_add_scalar_arg(param, m->val[0]);
}

/* Handle the case of a matrix parameter for which a scalar
   was given as argument.
*/

static int do_scalar_matrix_cast (fn_arg *arg, fn_param *param)
{
    gretl_matrix *m = gretl_matrix_alloc(1, 1);
    int err = 0;

    if (m == NULL) {
        err = E_ALLOC;
    } else {
        m->val[0] = arg->val.x;
        err = copy_as_arg(param->name, GRETL_TYPE_MATRIX, m);
        gretl_matrix_free(m);
    }

    return err;
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

static int upnames_match (fn_arg *ai, fn_arg *aj)
{
    return ai->upname != NULL && aj->upname != NULL &&
        strcmp(ai->upname, aj->upname) == 0;
}

/* Note: a caller cannot be allowed to supply a given variable in
   pointer form for more than one argument slot in a function call
   (although ordinary arguments may be repeated).
*/

static int duplicated_pointer_arg_check (fncall *call)
{
    fn_arg *ai, *aj;
    int i, j, err = 0;

    for (i=0; i<call->argc && !err; i++) {
        ai = &call->args[i];
        if (gretl_ref_type(ai->type)) {
            for (j=i+1; j<call->argc && !err; j++) {
                aj = &call->args[j];
                if (gretl_ref_type(aj->type) && upnames_match(ai, aj)) {
                    gretl_errmsg_set(_("Duplicated pointer argument: not allowed"));
                    err = E_DATA;
                }
            }
        }
    }

    return err;
}

static int process_series_arg (fncall *call, fn_param *fp,
                               fn_arg *arg, DATASET *dset)
{
    if (arg->type == GRETL_TYPE_USERIES) {
        /* an existing named series */
        if (param_is_const(fp) && arg->val.idnum != 0) {
            /* we can pass it by reference, if the series
	       in question is not the one named "const",
	       with ID number 0
	    */
            return localize_series_ref(call, arg, fp, dset);
        } else {
            /* pass it by value */
            return dataset_copy_series_as(dset, arg->val.idnum, fp->name);
        }
    } else if (arg->val.px != NULL) {
        /* an on-the-fly series */
        return dataset_add_series_as(dset, arg->val.px, fp->name);
    } else {
        return E_DATA;
    }
}

/* for matrix, bundle, array, string */

static int process_object_arg (fncall *call, int i, fn_param *fp)
{
    fn_arg *arg = &call->args[i];

    if (param_is_const(fp)) {
        /* we can pass it by reference */
        return localize_const_object(call, i, fp);
    } else {
        /* pass it by value */
        return copy_as_arg(fp->name, arg->type,
                           arg_get_data(arg, 0));
    }
}

/* for *matrix, *bundle, *array, *string */

static int process_object_ref_arg (fn_arg *arg, fn_param *fp)
{
    if (arg->upname != NULL) {
        return user_var_localize(arg->upname, fp->name, fp->type);
    } else {
        return localize_object_as_shell(arg, fp);
    }
}

/* Note that if we reach here we've already successfully negotiated
   check_function_args(): now we're actually making the argument
   objects (if any) available within the function.
*/

static int allocate_function_args (fncall *call, DATASET *dset)
{
    ufunc *fun = call->fun;
    fn_arg *arg;
    fn_param *fp;
    int i, err;

    err = duplicated_pointer_arg_check(call);

    for (i=0; i<call->argc && !err; i++) {
        arg = &call->args[i];
        fp = &fun->params[i];

#if ARGS_DEBUG
        fprintf(stderr, "arg[%d], param type %s (%s), arg type %s (%s)\n",
                i, gretl_type_get_name(fp->type), fp->name,
                gretl_type_get_name(arg->type), arg->upname);
        fprintf(stderr, "  const %s, optional %s\n",
                param_is_const(fp) ? "yes" : "no",
                param_is_optional(fp) ? "yes" : "no");
#endif
        if (arg->upname != NULL && object_is_const(arg->upname, -1) &&
            gretl_ref_type(fp->type) && !param_is_const(fp)) {
            const char *caller = NULL;

            current_function_info(&caller, NULL);
            gretl_errmsg_sprintf(_("%s() tries to pass const argument "
                                 "%s to %s() in pointer form with no\n"
                                 "const guarantee\n"),
                                 caller, fp->name, call->fun->name);
            err = E_INVARG;
            break;
        }

        if (arg->type == GRETL_TYPE_NONE) {
            if (gretl_scalar_type(fp->type)) {
                err = add_scalar_arg_default(fp);
            } else if (fp->type == GRETL_TYPE_LIST) {
                err = add_empty_list(call, fp, dset);
            }
            continue;
        }

        if (fp->type == GRETL_TYPE_NUMERIC) {
            /* param supports overloading */
            if (gretl_scalar_type(arg->type)) {
                err = real_add_scalar_arg(fp, arg->val.x);
            } else if (gretl_is_series_type(arg->type)) {
                process_series_arg(call, fp, arg, dset);
            } else if (arg->type == GRETL_TYPE_MATRIX) {
                err = process_object_arg(call, i, fp);
            } else if (arg->type == GRETL_TYPE_LIST) {
                err = localize_list(call, arg, fp, dset);
            }
            continue;
        }

        if (gretl_scalar_type(fp->type)) {
            if (arg->type == GRETL_TYPE_MATRIX) {
                err = do_matrix_scalar_cast(arg, fp);
            } else {
                err = real_add_scalar_arg(fp, arg->val.x);
            }
        } else if (fp->type == GRETL_TYPE_SERIES) {
            err = process_series_arg(call, fp, arg, dset);
        } else if (fp->type == GRETL_TYPE_MATRIX &&
                   arg->type == GRETL_TYPE_DOUBLE) {
            err = do_scalar_matrix_cast(arg, fp);
        } else if (fp->type == GRETL_TYPE_MATRIX ||
                   fp->type == GRETL_TYPE_BUNDLE ||
                   fp->type == GRETL_TYPE_STRING ||
                   gretl_array_type(fp->type)) {
            err = process_object_arg(call, i, fp);
        } else if (fp->type == GRETL_TYPE_LIST) {
            err = localize_list(call, arg, fp, dset);
        } else if (fp->type == GRETL_TYPE_SERIES_REF) {
            err = localize_series_ref(call, arg, fp, dset);
        } else if (gretl_ref_type(fp->type)) {
            err = process_object_ref_arg(arg, fp);
        }
        if (!err && (fp->type == GRETL_TYPE_BUNDLE ||
                     fp->type == GRETL_TYPE_BUNDLE_REF)) {
            err = localize_bundled_lists(call, arg, fp, dset);
        }
    }

    /* now for any trailing parameters without matching arguments */

    for (i=call->argc; i<fun->n_params && !err; i++) {
        fp = &fun->params[i];
        if (gretl_scalar_type(fp->type)) {
            err = add_scalar_arg_default(fp);
        } else if (fp->type == GRETL_TYPE_LIST) {
            err = add_empty_list(call, fp, dset);
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

#if ARGS_DEBUG
    fprintf(stderr, "allocate_function args: returning %d\n", err);
#endif

    return err;
}

/**
 * check_function_needs:
 * @dset: pointer to dataset info.
 * @dreq: function data requirements flag.
 * @minver: function minimum program version requirement.
 * @pkg: the package to which the function belongs.
 *
 * Checks whether the requirements in @dreq and @minver
 * are jointly satisfied by the current dataset and gretl
 * version.
 *
 * Returns: 0 if all is OK; E_DATA if the current dataset
 * (or lack thereof) does not satisfy @dreq; 1 otherwise.
 */

int check_function_needs (const DATASET *dset, DataReq dreq,
                          int minver, fnpkg *pkg)
{
    static int thisver = 0;

    if (thisver == 0) {
        thisver = gretl_version_number(GRETL_VERSION);
    }

    if (minver > thisver) {
        char vstr[8];

        gretl_version_string(vstr, minver);
        if (pkg != NULL) {
            gretl_errmsg_sprintf(_("The package %s needs gretl version %s"),
                                 pkg->name, vstr);
        }
        return 1;
    }

    if ((dset == NULL || dset->v == 0) && dreq != FN_NODATA_OK) {
        if (pkg != NULL) {
            gretl_errmsg_sprintf(_("The package %s needs a dataset in place"),
                                 pkg->name);
        }
        return E_DATA;
    }

    if (dreq == FN_NEEDS_TS && !dataset_is_time_series(dset)) {
        if (pkg != NULL) {
            gretl_errmsg_sprintf(_("The package %s needs time-series data"),
                                 pkg->name);
        }
        return E_DATA;
    }

    if (dreq == FN_NEEDS_PANEL && !dataset_is_panel(dset)) {
        if (pkg != NULL) {
            gretl_errmsg_sprintf(_("The package %s needs panel data"),
                                 pkg->name);
        }
        return E_DATA;
    }

    if (dreq == FN_NEEDS_QM && !quarterly_or_monthly(dset)) {
         if (pkg != NULL) {
            gretl_errmsg_sprintf(_("The package %s needs quarterly or monthly data"),
                                 pkg->name);
        }
        return E_DATA;
    }

    if (pkg != NULL) {
        pkg->prechecked = 1;
    }

    return 0;
}

/**
 * package_version_ok:
 * @minver: function package minimum program version requirement.
 * @reqstr: location to write required version string (should be
 * at least 8 bytes), or NULL.
 *
 * Returns: 1 if the running version of gretl satisfies the
 * minimum version requirement in @minver; otherwise returns 0,
 * in which case the string representation of @minver is written
 * @reqstr if this is non-NULL.
 */

int package_version_ok (int minver, char *reqstr)
{
    static int thisver = 0;
    int ret = 0;

    if (thisver == 0) {
        thisver = gretl_version_number(GRETL_VERSION);
    }

    ret = thisver >= minver;

    if (!ret && reqstr != NULL) {
        gretl_version_string(reqstr, minver);
    }

    return ret;
}

static int maybe_check_function_needs (const DATASET *dset,
                                       const ufunc *fun)
{
    if (fun->pkg == NULL || fun->pkg->prechecked) {
        return 0;
    }

#ifdef HAVE_MPI
    if (*mpi_caller != '\0') {
        /* let's assume we're OK */
        return 0;
    }
#endif

    return check_function_needs(dset, fun->pkg->dreq,
                                fun->pkg->minver, fun->pkg);
}

/* next block: handling function return values */

static int handle_scalar_return (fncall *call, void *ptr, int rtype)
{
    static double xret;
    const char *vname = call->retname;
    int err = 0;

    if (gretl_is_scalar(vname)) {
        xret = gretl_scalar_get_value(vname, NULL);
    } else {
        gretl_matrix *m = get_matrix_by_name(vname);

        if (gretl_matrix_is_scalar(m)) {
            xret = m->val[0];
        } else {
            user_var *uv = get_user_var_by_name(vname);
            GretlType t = user_var_get_type(uv);

            gretl_errmsg_sprintf(_("Function %s did not provide the specified return value\n"
                                 "(expected %s, got %s)"), call->fun->name,
                                 gretl_type_get_name(rtype), gretl_type_get_name(t));
            err = E_TYPES;
            xret = NADBL;
        }
    }

    *(double **) ptr = &xret;

    return err;
}

static int handle_series_return (const char *vname, void *ptr,
                                 DATASET *dset, int copy,
                                 char **label, series_table **stab)
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

    if (!err) {
        if (label != NULL) {
            const char *vstr = series_get_label(dset, v);

            if (vstr != NULL) {
                *label = gretl_strdup(vstr);
            } else {
                *label = gretl_strdup("");
            }
        }
        if (stab != NULL && is_string_valued(dset, v)) {
            series_table *st = series_get_string_table(dset, v);

            *stab = series_table_copy(st);
        }
    }

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

static GretlType fix_numeric_return_type (fncall *call,
                                          DATASET *dset)
{
    GretlType t = 0;
    user_var *uv;

    uv = user_var_get_value_and_type(call->retname, &t);
    if (uv != NULL) {
        call->rtype = t;
    } else {
        int v = current_series_index(dset, call->retname);

        if (v >= 0) {
            t = call->rtype = GRETL_TYPE_SERIES;
        }
    }

    return t;
}

static int handle_bundle_return (fncall *call, void *ptr, int copy)
{
    const char *name = call->retname;
    gretl_bundle *ret = NULL;
    int err = 0;

    if (copy) {
        gretl_bundle *b = get_bundle_by_name(name);

        if (b != NULL) {
            ret = gretl_bundle_copy(b, &err);
        }
    } else {
        ret = gretl_bundle_pull_from_stack(name, &err);
    }

    if (ret != NULL && call->fun->pkg != NULL &&
        gretl_function_depth() == 1) {
        if (call->fun->pkg_role != UFUN_BUNDLE_FCAST) {
            gretl_bundle_set_creator(ret, call->fun->pkg->name);
        }
    }

    *(gretl_bundle **) ptr = ret;

    return err;
}

static int handle_array_return (fncall *call, void *ptr, int copy)
{
    const char *name = call->retname;
    gretl_array *ret = NULL;
    int err = 0;

    if (copy) {
        gretl_array *a = get_array_by_name(name);

        if (a != NULL) {
            ret = gretl_array_copy(a, &err);
        }
    } else {
        ret = gretl_array_pull_from_stack(name, &err);
    }

    *(gretl_array **) ptr = ret;

    return err;
}

static void replace_caller_series (int targ, int src, DATASET *dset)
{
    const char *vlabel = series_get_label(dset, src);
    int t;

    /* replace data values */
    for (t=dset->t1; t<=dset->t2; t++) {
        dset->Z[targ][t] = dset->Z[src][t];
    }

    /* replace variable info? */
    series_set_label(dset, targ, vlabel == NULL ? "" : vlabel);
    series_set_display_name(dset, targ, series_get_display_name(dset, src));
}

/* Deal with a list, named @lname, that exists at the level of a
   function whose execution is now terminating.  Note that this list
   may be the direct return value of the function (in which case @ret
   will be non-NULL), or it may have been passed to the function as an
   argument (in which case @arg will be non-NULL), or it may have been
   constructed on the fly.

   If the list was given as a function argument we simply shunt all
   its members to the caller's stack level.  But if it's the direct
   return value we need to overwrite any variables at caller level
   that have been redefined within the function.  If any series have
   been redefined in that way we also need to adjust the list itself,
   replacing the ID numbers of local variables with the caller-level
   IDs -- all supposing the direct return value is being assigned by
   the caller.

   If the list was a temporary construction (a list parameter was
   required, but a single series was given as argument), we just
   destroy it.
*/

static int unlocalize_list (fncall *call, const char *lname,
                            fn_arg *arg, void *ret,
                            DATASET *dset)
{
    int *list = get_list_by_name(lname);
    int d = gretl_function_depth();
    int upd = d - 1;
    int i, vi;

#if ARGS_DEBUG
    fprintf(stderr, "\n*** unlocalize_list '%s', function %s depth = %d\n",
            lname, call->fun->name, d);
    printlist(list, lname);
    fprintf(stderr, " dset = %p, dset->v = %d\n", (void *) dset, dset->v);
    fprintf(stderr, " list is direct return value? %s\n", arg == NULL ? "yes" : "no");
#endif

    if (list == NULL) {
        return E_DATA;
    }

    if (arg == NULL && ret == NULL) {
        ; /* no-op: we'll trash the list later */
    } else if (arg == NULL) {
        /* handle list as direct return value */
        int j, lev, overwrite;
        const char *vname;

        for (i=1; i<=list[0]; i++) {
            overwrite = 0;
            vi = list[i];
            vname = dset->varname[vi];
            if (in_gretl_list(call->ptrvars, vi)) {
                /* 2021-06-28: detect a somewhat anomalous case */
                fprintf(stderr, "*** return list '%s' contains series %s, passed in "
                        "pointer form ***\n", lname, vname);
                continue;
            }
            series_unset_flag(dset, vi, VAR_LISTARG);
            if (vi > 0 && vi < dset->v && series_get_stack_level(dset, vi) == d) {
                for (j=1; j<dset->v; j++) {
                    lev = series_get_stack_level(dset, j);
                    if (lev == upd && !strcmp(dset->varname[j], vname)) {
                        overwrite = 1;
                        break;
                    }
                    if (lev == d && j < vi && series_is_listarg(dset, j, NULL) &&
                        !strcmp(dset->varname[j], vname)) {
                        overwrite = 1;
                        break;
                    }
                }
                if (overwrite) {
                    replace_caller_series(j, vi, dset);
                    /* replace ID number in list */
                    list[i] = j;
                } else {
                    series_set_stack_level(dset, vi, upd);
                }
            }
#if ARGS_DEBUG
            fprintf(stderr, " list-member var %d, '%s': ", vi, vname);
            if (overwrite) {
                fprintf(stderr, "found match in caller, overwrote var %d\n", j);
            } else {
                fprintf(stderr, "no match in caller\n");
            }
#endif
        }
    } else {
        /* list was given as argument to function */
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == LISTSEP) {
                continue;
            }
            if (series_is_listarg(dset, vi, NULL)) {
                series_unset_flag(dset, vi, VAR_LISTARG);
                series_set_stack_level(dset, vi, upd);
            }
        }
        if (arg->type != GRETL_TYPE_LIST) {
            /* the list was constructed on the fly */
            user_var_delete_by_name(lname, NULL);
        }
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

static int check_sub_object_return (fn_arg *arg, fn_param *fp)
{
    GretlType type = gretl_type_get_plain_type(arg->type);
    void *orig = arg_get_data(arg, type);
    void *curr = user_var_get_value_by_name(fp->name);
    int err = 0;

    if (curr == NULL) {
        fprintf(stderr, "sub_object return: value has disappeared!\n");
        err = E_DATA;
    } else if (curr != orig) {
        /* reattachment is needed */
        GretlType uptype = user_var_get_type(arg->uvar);
        void *updata = user_var_get_value(arg->uvar);

        if (uptype == GRETL_TYPE_ARRAY) {
            gretl_array *a = updata;
            int i, n = gretl_array_get_length(a);
            int found = 0;

            for (i=0; i<n; i++) {
                if (gretl_array_get_data(a, i) == orig) {
                    gretl_array_set_data(a, i, curr);
                    found = 1;
                }
            }
            if (!found) {
                err = E_DATA;
            }
        } else if (uptype == GRETL_TYPE_BUNDLE) {
            fprintf(stderr, "sub_object return: uptype = bundle, not handled\n");
            err = E_DATA;
        } else {
            /* no other types handled */
            err = E_TYPES;
        }
    }

    return err;
}

static int is_pointer_arg (fncall *call, int rtype)
{
    ufunc *u = call->fun;

    if (call->retname != NULL) {
        fn_param *fp;
        int i;

        for (i=0; i<call->argc; i++) {
            fp = &u->params[i];
            if ((fp->type == gretl_type_get_ref_type(rtype) ||
                 (fp->flags & FP_CONST)) &&
                strcmp(fp->name, call->retname) == 0) {
                return 1;
            }
        }
    }

    return 0;
}

/* On completing execution of a function, move a series at the current
   level up to the level of the caller. This must be done for series
   that were passed to the function by reference rather than being
   copied -- that is, series passed in explicit pointer form and also
   series "loaned" by the caller via the "const" qualifier.
*/

static void push_series_to_caller (ufunc *u, fn_arg *arg, DATASET *dset)
{
    int v = arg->val.idnum;

    if (v == 0) {
	/* this is a no-op */
	return;
    } else if (arg->upname == NULL) {
	/* an internal error! */
        fprintf(stderr, "ERROR in push_series_to_caller: arg->upname is NULL\n");
        return;
    }

    series_decrement_stack_level(dset, v);
    if (series_get_stack_level(dset, v) < 0) {
        fprintf(stderr, "@@@ After decrement in %s, stack level=%d for %s @@@\n",
                u->name, series_get_stack_level(dset, v), arg->upname);
        /* set_stack_level(dset, v, 0); */
    }
    strcpy(dset->varname[v], arg->upname);
}

static void push_object_to_caller (fn_arg *arg)
{
    user_var_adjust_level(arg->uvar, -1);
    user_var_set_name(arg->uvar, arg->upname);
}

#define null_return(t) (t == GRETL_TYPE_VOID || t == GRETL_TYPE_NONE)

#define needs_dataset(t) (t == GRETL_TYPE_SERIES || \
                          t == GRETL_TYPE_LIST ||   \
                          t == GRETL_TYPE_SERIES_REF || \
                          t == GRETL_TYPE_LISTS || \
                          t == GRETL_TYPE_LISTS_REF)

static int
function_assign_returns (fncall *call, int rtype,
                         DATASET *dset, void *ret,
                         char **label, series_table **stab,
                         PRN *prn, int *perr)
{
    ufunc *u = call->fun;
    int i, err = 0;

#if ARGS_DEBUG
    fprintf(stderr, "function_assign_returns: rtype = %s, call->retname = %s\n",
            gretl_type_get_name(rtype), call->retname);
#endif

    if (*perr == 0 && !null_return(rtype) && call->retname == NULL) {
        /* missing return value */
        gretl_errmsg_sprintf(_("Function %s did not provide the specified return value\n"
                             "(expected %s)"), u->name, gretl_type_get_name(rtype));
        *perr = err = E_UNKVAR;
    } else if (*perr == 0 && needs_dataset(rtype) && dset == NULL) {
        /* "can't happen" */
        *perr = err = E_DATA;
    }

    if (*perr == 0) {
        /* first we work on the value directly returned by the
           function (but only if there's no error)
        */
        int copy = is_pointer_arg(call, rtype);

        if (rtype == GRETL_TYPE_NUMERIC) {
            rtype = fix_numeric_return_type(call, dset);
        }

        if (rtype == GRETL_TYPE_DOUBLE) {
            err = handle_scalar_return(call, ret, rtype);
        } else if (rtype == GRETL_TYPE_SERIES) {
            err = handle_series_return(call->retname, ret, dset, copy, label, stab);
        } else if (rtype == GRETL_TYPE_MATRIX) {
            err = handle_matrix_return(call->retname, ret, copy);
        } else if (rtype == GRETL_TYPE_LIST) {
            /* note: in this case the job is finished in
               stop_fncall(); here we just adjust the info on the
               listed variables so they don't get deleted
            */
            err = unlocalize_list(call, call->retname, NULL, ret, dset);
        } else if (rtype == GRETL_TYPE_BUNDLE) {
            err = handle_bundle_return(call, ret, copy);
        } else if (gretl_array_type(rtype)) {
            err = handle_array_return(call, ret, copy);
        } else if (rtype == GRETL_TYPE_STRING) {
            err = handle_string_return(call->retname, ret);
        }

        if (err == E_UNKVAR) {
            pprintf(prn, _("Function %s did not provide the specified return value\n"),
                    u->name);
        }

        *perr = err;
    }

    /* "Indirect return" values and other pointerized args:
       these should be handled even if the function bombed.
    */

    for (i=0; i<call->argc; i++) {
        fn_arg *arg = &call->args[i];
        fn_param *fp = &u->params[i];
        int ierr = 0;

        if (needs_dataset(fp->type) && dset == NULL) {
            ierr = E_DATA;
        } else if (gretl_ref_type(fp->type)) {
            if (arg->type == GRETL_TYPE_SERIES_REF) {
                push_series_to_caller(u, arg, dset);
            } else if (arg->upname != NULL) {
                push_object_to_caller(arg);
            } else if (arg->uvar != NULL) {
                ierr = check_sub_object_return(arg, fp);
            } else {
                ; /* pure "shell" object: no-op */
            }
        } else if (arg->type == GRETL_TYPE_USERIES &&
                   (fp->type == GRETL_TYPE_SERIES ||
                    fp->type == GRETL_TYPE_NUMERIC) &&
                   param_is_const(fp)) {
            push_series_to_caller(u, arg, dset);
        } else if ((fp->type == GRETL_TYPE_MATRIX ||
                    fp->type == GRETL_TYPE_BUNDLE ||
                    fp->type == GRETL_TYPE_STRING ||
                    fp->type == GRETL_TYPE_NUMERIC ||
                    gretl_array_type(fp->type))) {
            if (arg->shifted) {
                /* non-pointerized const object argument,
                   which we renamed and shifted
                */
                push_object_to_caller(arg);
            }
        } else if (fp->type == GRETL_TYPE_LIST) {
            unlocalize_list(call, fp->name, arg, NULL, dset);
        }
        if (ierr) {
            *perr = err = ierr;
        }
    }

#if ARGS_DEBUG
    fprintf(stderr, "function_assign_returns: returning %d\n", err);
#endif

    return err;
}

/* make a record of the sample information at the time a function is
   called */

static void record_obs_info (obsinfo *oi, DATASET *dset)
{
    oi->added = 0;
    oi->structure = dset->structure;
    oi->pd = dset->pd;
    oi->t1 = dset->t1;
    oi->t2 = dset->t2;
    strcpy(oi->stobs, dset->stobs);
    oi->panel_pd = dset->panel_pd;
    oi->panel_sd0 = dset->panel_sd0;
}

/* On function exit, restore the observations information
   ("setobs" stuff) that was in force on entry.
*/

static void restore_obs_info (obsinfo *oi, DATASET *dset)
{
    dset->structure = oi->structure;
    dset->pd = oi->pd;
    dset->t1 = oi->t1;
    dset->t2 = oi->t2;
    strcpy(dset->stobs, oi->stobs);
    dset->panel_pd = oi->panel_pd;
    dset->panel_sd0 = oi->panel_sd0;
}

static void push_verbosity (fncall *call)
{
    if (!is_recursing(call)) {
        if (gretl_messages_on()) {
            call->flags |= FC_PREV_MSGS;
        } else {
            call->flags &= ~FC_PREV_MSGS;
        }
        if (gretl_echo_on()) {
            call->flags |= FC_PREV_ECHO;
        } else {
            call->flags &= ~FC_PREV_ECHO;
        }
    }
}

static void pop_verbosity (fncall *call)
{
    if (!is_recursing(call)) {
        set_gretl_messages(call->flags & FC_PREV_MSGS);
        set_gretl_echo(call->flags & FC_PREV_ECHO);
    }
}

/* do the basic housekeeping that is required when a function exits:
   destroy local variables, restore previous sample info, etc.
*/

static int stop_fncall (fncall *call, int rtype, void *ret,
                        DATASET *dset, PRN *prn, int redir_level)
{
    int i, d = gretl_function_depth();
    int delv, anyerr = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "stop_fncall: terminating call to "
            "function '%s' at depth %d, dset->v = %d\n",
            call->fun->name, d, (dset != NULL)? dset->v : 0);
#endif

    /* below: delete series local to the function, taking care not to
       delete any that have been "promoted" to caller level via their
       inclusion in a returned list; also catch any "listarg" series
       that are no longer attached to a list
    */

    if (dset != NULL) {
        for (i=call->orig_v, delv=0; i<dset->v; i++) {
            if (series_get_stack_level(dset, i) == d) {
                delv++;
            }
        }
        if (delv > 0) {
            if (delv == dset->v - call->orig_v) {
                /* deleting all added series */
                anyerr = dataset_drop_last_variables(dset, delv);
                if (anyerr && !err) {
                    err = anyerr;
                }
            } else {
                for (i=dset->v-1; i>=call->orig_v; i--) {
                    if (series_get_stack_level(dset, i) == d) {
                        anyerr = dataset_drop_variable(i, dset);
                        if (anyerr && !err) {
                            err = anyerr;
                        }
                    }
                }
            }
        }
        if (call->listvars != NULL) {
            int vi;

            for (i=1; i<=call->listvars[0]; i++) {
                vi = call->listvars[i];
                series_unset_flag(dset, vi, VAR_LISTARG);
                series_set_stack_level(dset, vi, d - 1);
            }
        }
    }

    /* direct list return: write the possibly revised list to the
       return pointer.  Note that we can't do this earlier, because
       the ID numbers of the variables in the return list may be
       changed due to the deletion of function-local variables.
    */

    if (!err && rtype == GRETL_TYPE_LIST && ret != NULL) {
        int *lret = gretl_list_copy(get_list_by_name(call->retname));

        if (lret != NULL) {
            *(int **) ret = lret;
        } else {
            err = E_ALLOC;
        }
    }

    /* now we're ready to trash function-local vars */
    anyerr = destroy_user_vars_at_level(d);

    if (anyerr && !err) {
        err = anyerr;
    }

    /* if any anonymous equations system was defined: clean up */
    delete_anonymous_equation_system(d);

    if (call->fun->flags & UFUN_USES_SET) {
        pop_program_state();
    } else {
        pop_verbosity(call);
    }

    if (dset != NULL && dset->v > 0) {
        restore_obs_info(&call->obs, dset);
    }

    set_executing_off(&call, dset, prn);

    if (print_redirection_level(prn) > redir_level) {
        gretl_errmsg_set(_("Incorrect use of 'outfile' in function"));
        err = 1;
    }

    if (call != NULL && call->retname != NULL) {
        /* we're surely done with this now: avoid potential leak */
        free(call->retname);
        call->retname = NULL;
    }

    return err;
}

/* reset any saved uservar addresses in context of loops or
   genrs attached to @u */

static void reset_saved_uservars (ufunc *u, int on_recurse)
{
    fn_line *line;
    int i;

#if COMP_DEBUG
    if (on_recurse) {
        fprintf(stderr, "on recursion of %s, reset_saved_uservars\n", u->name);
    } else {
        fprintf(stderr, "at exit from %s, reset_saved_uservars\n", u->name);
    }
#endif

    for (i=0; i<u->n_lines; i++) {
        line = &(u->lines[i]);
        if (line_has_loop(line)) {
            loop_reset_uvars(line->ptr);
        }
    }

#if 1
    /* 2022-09-11: we shouldn't have to do genr_reset_uvars (below)
       any more, with genrs now attached to fncall rather than line of
       function?
    */
    if (on_recurse) return;
#endif

    if (u->call != NULL && u->call->lgen != NULL) {
        for (i=0; i<u->call->n_lgen; i++) {
            genr_reset_uvars(u->call->lgen[i].genr);
        }
    }
}

static void set_pkgdir (fnpkg *pkg)
{
    const char *p = strrslash(pkg->fname);

    if (p != NULL) {
        char *pkgdir = gretl_strndup(pkg->fname, p - pkg->fname);

        gretl_insert_builtin_string("pkgdir", pkgdir);
        free(pkgdir);
    }
}

static int start_fncall (fncall *call, DATASET *dset, PRN *prn)
{
    GList *tmp = g_list_last(callstack);
    fncall *prevcall;

    while (tmp != NULL) {
        prevcall = tmp->data;
        if (prevcall->fun == call->fun) {
            set_call_recursing(call);
            reset_saved_uservars(call->fun, 1); /* 2022-08-27 */
            break;
        }
        tmp = tmp->prev;
    }

    set_previous_depth(fn_executing);
    fn_executing++;

    if (call->fun->flags & UFUN_USES_SET) {
        push_program_state();
    } else {
        push_verbosity(call);
    }

    callstack = g_list_append(callstack, call);

#if EXEC_DEBUG
    fprintf(stderr, "start_fncall: %s, depth now %d, recursing %d\n",
            call->fun->name, g_list_length(callstack), is_recursing(call));
#endif

    switch_uservar_hash(fn_executing);

    if (dset != NULL && dset->v > 0) {
        record_obs_info(&call->obs, dset);
    }

    if (gretl_debugging_on()) {
        set_gretl_echo(1);
        set_gretl_messages(1);
        pputs(prn, "*** ");
        bufspace(gretl_function_depth(), prn);
        pprintf(prn, "executing function %s\n", call->fun->name);
    } else {
        set_gretl_echo(0);
        set_gretl_messages(0);
    }

    if (fn_executing == 1 && call->fun->pkg != NULL) {
        set_pkgdir(call->fun->pkg);
    }

    return 0;
}

static int func_exec_callback (ExecState *s, void *ptr,
                               GretlObjType type)
{
    if (s->cmd->ci == FLUSH || is_plotting_command(s->cmd)) {
        /* we permit "reach-back" into the GUI for these */
        EXEC_CALLBACK gc = get_gui_callback();

        if (gc != NULL) {
            gc(s, NULL, 0);
        }
    }

    return 0;
}

static double arg_get_double_val (fn_arg *arg)
{
    if (gretl_scalar_type(arg->type)) {
        return arg->val.x;
    } else if (arg->type == GRETL_TYPE_MATRIX) {
        return arg->val.m->val[0];
    } else {
        return NADBL;
    }
}

static int check_function_args (fncall *call, PRN *prn)
{
    ufunc *u = call->fun;
    fn_arg *arg;
    fn_param *fp;
    double x;
    int i, err = 0;

    for (i=0; i<call->argc && !err; i++) {
        int scalar_check = 1;

        arg = &call->args[i];
        fp = &u->params[i];

        if (param_is_optional(fp) && arg->type == GRETL_TYPE_NONE) {
            ; /* this is OK */
        } else if (gretl_scalar_type(fp->type) && arg->type == GRETL_TYPE_DOUBLE) {
            ; /* OK: types match */
        } else if (fp->type == GRETL_TYPE_SERIES && arg->type == GRETL_TYPE_USERIES) {
            ; /* OK: types match */
        } else if (gretl_scalar_type(fp->type) &&
                   arg->type == GRETL_TYPE_MATRIX &&
                   gretl_matrix_is_scalar(arg->val.m)) {
            ; /* OK: either types match or we can convert */
        } else if (fp->type == GRETL_TYPE_MATRIX && arg->type == GRETL_TYPE_DOUBLE) {
            ; /* OK, we can convert */
        } else if (fp->type == GRETL_TYPE_LIST && arg->type == GRETL_TYPE_USERIES) {
            ; /* OK, we'll handle it (create singleton list) */
        } else if (fp->type == GRETL_TYPE_LIST && arg->type == GRETL_TYPE_NONE) {
            ; /* OK (special: "null" was passed as argument) */
        } else if (gretl_scalar_type(fp->type) && arg->type == GRETL_TYPE_NONE &&
                   !no_scalar_default(fp)) {
            scalar_check = 0; /* OK: arg missing but we have a default value */
        } else if (fp->type == GRETL_TYPE_NUMERIC && NUMERIC_TYPE(arg->type)) {
            ; /* OK, for overloaded param */
        } else if (fp->type != arg->type) {
            pprintf(prn, _("%s: argument %d is of the wrong type (is %s, should be %s)\n"),
                    u->name, i + 1, gretl_type_get_name(arg->type),
                    gretl_type_get_name(fp->type));
            err = E_TYPES;
        }

        if (!err && scalar_check && fp->type == GRETL_TYPE_DOUBLE) {
            /* check for out-of-bounds scalar value */
            x = arg_get_double_val(arg);
            if ((!na(fp->min) && x < fp->min) ||
                (!na(fp->max) && x > fp->max)) {
                pprintf(prn, _("%s, argument %d: value %g is out of bounds\n"),
                        u->name, i + 1, x);
                err = E_INVARG;
            }
        }
    }

    for (i=call->argc; i<u->n_params && !err; i++) {
        /* do we have defaults for any empty args? */
        fp = &u->params[i];
        if (!param_is_optional(fp) && no_scalar_default(fp)) {
            pprintf(prn, _("%s: not enough arguments\n"), u->name);
            err = E_ARGS;
        }
    }

#if ARGS_DEBUG
    fprintf(stderr, "check_function_args: err = %d\n", err);
#endif

    return err;
}

static const char *get_funcerr_message (ExecState *state)
{
    const char *p = state->cmd->param;

    if (p != NULL && strlen(p) == gretl_namechar_spn(p)) {
        const char *s = get_string_by_name(p);

        if (s != NULL) {
            return s;
        }
    } else if (p == NULL) {
        return "none";
    }

    return p;
}

static void set_func_error_message (int err, ufunc *u,
                                    ExecState *state,
                                    const char *cmdline,
                                    fn_line *fline)
{
    if (err == E_FUNCERR) {
        /* let the function writer set the message? */
        const char *msg = get_funcerr_message(state);

        if (msg != NULL && strcmp(msg, "none")) {
            gretl_errmsg_sprintf(_("Error message from %s():\n %s"),
                                 u->name, msg);
            return;
        }
    }

    if (err == E_STOP) {
        ; /* no-op */
    } else {
        /* we'll handle this here */
        const char *msg = gretl_errmsg_get();
        int showline = 1;

        if (*msg == '\0') {
            msg = errmsg_get_with_default(err);
            gretl_errmsg_set(msg);
        }

        if (*cmdline == '\0' || strncmp(cmdline, "funcerr(", 8) == 0 ||
            strncmp(cmdline, "errorif(", 8) == 0) {
            showline = 0;
        }

        if (fline == NULL) {
            /* indicates that we don't have a real line number
               available (the error occurred inside a loop)
            */
            if (showline && err != E_FNEST) {
                gretl_errmsg_sprintf(_("*** error within loop in function %s\n> %s"),
                                     u->name, cmdline);
            } else {
                gretl_errmsg_sprintf(_("*** error within loop in function %s\n"),
                                     u->name);
            }
        } else {
            if (showline && err != E_FNEST) {
		/* 2024-05-12: last arg here was cmdline */
                gretl_errmsg_sprintf(_("*** error in function %s, line %d\n> %s"),
                                     u->name, fline->idx, fline->s);
            } else {
                gretl_errmsg_sprintf(_("*** error in function %s, line %d\n> %s\n"),
                                     u->name, fline->idx, fline->s);
            }
        }

        if (gretl_function_depth() > 1) {
            GList *tmp = g_list_last(callstack);

            if (tmp != NULL) {
                tmp = g_list_previous(tmp);
                if (tmp != NULL) {
                    fncall *call = tmp->data;

                    gretl_errmsg_sprintf(_(" called by function %s"), call->fun->name);
                }
            }
        }
    }
}

static GENERATOR *fnline_get_genr (fncall *call, fn_line *line)
{
    GENERATOR *genr = NULL;
    int i;

    for (i=0; i<call->n_lgen; i++) {
        if (call->lgen[i].idx == line->idx) {
            genr = call->lgen[i].genr;
            break;
        }
    }

#if REC_DEBUG
    fprintf(stderr, "?? seek genr for call %p line idx %d: got %p (n_lgen %d)\n",
            (void *) call, line->idx, (void *) genr, call->n_lgen);
#endif

    return genr;
}

static int fnline_set_genr (fncall *call, fn_line *line, GENERATOR *genr)
{
    linegen *lgen;
    int n = call->n_lgen;
    int err = 0;

    lgen = realloc(call->lgen, (n + 1) * sizeof *lgen);
    if (lgen == NULL) {
        err = E_ALLOC;
    } else {
        call->lgen = lgen;
        call->lgen[n].idx = line->idx;
        call->lgen[n].genr = genr;
        call->n_lgen = n + 1;
#if REC_DEBUG
        fprintf(stderr, "++ adding genr %p for call %p line idx %d (n_lgen %d)\n",
                (void *) *pgen, (void *) call, line->idx, call->n_lgen);
#endif
    }

    return err;
}

static int generate_return_value2 (fncall *call,
                                   ExecState *state,
                                   DATASET *dset,
                                   const char *s,
                                   fn_line *line)
{
    gchar *formula = g_strdup_printf("$retval=%s", s);
    int done = 0;
    int err = 0;

    if (line != NULL) {
        /* try compilation */
        GENERATOR *genr;

        genr = genr_compile(formula, dset, call->fun->rettype,
                            OPT_P, state->prn, &err);
        if (!err && genr != NULL) {
            /* succeeded */
            fnline_set_genr(call, line, genr);
            done = 1;
        } else if (err == E_EQN) {
            /* failed, but not fatally */
            fprintf(stderr, "did NOT attach genr at %p\n", line->ptr);
            line->flags |= LINE_NOCOMP;
            err = 0;
        }
    }
    if (!done && !err) {
        err = generate(formula, dset, call->fun->rettype,
                       OPT_P, state->prn);
    }
    g_free(formula);

    return err;
}

static int generate_return_value1 (fncall *call,
                                   ExecState *state,
                                   DATASET *dset,
                                   const char *s)
{
    gchar *formula = g_strdup_printf("$retval=%s", s);
    int err;

    err = generate(formula, dset, call->fun->rettype, OPT_P, state->prn);
    g_free(formula);

    return err;
}

static int handle_return_statement (fncall *call,
                                    ExecState *state,
                                    DATASET *dset,
                                    int gencomp,
                                    fn_line *line)
{
    ufunc *fun = call->fun;
    const char *s = NULL;
    int err = 0;

    if (line != NULL) {
        GENERATOR *genr = fnline_get_genr(call, line);

        if (genr != NULL) {
#if REC_DEBUG
            fprintf(stderr, "%s: handle_return: exec compiled genr %p\n",
                    fun->name, genr);
#endif
            /* FIXME: when exactly is this reset required? */
            genr_reset_uvars(genr);
            err = execute_genr(genr, dset, state->prn);
            goto last_step;
        } else {
            err = parse_command_line(state, dset, NULL);
            if (err) {
                goto last_step;
            }
        }
    }

    s = state->cmd->vstart;

#if EXEC_DEBUG
    fprintf(stderr, "%s: handle_return_statement '%s'\n", fun->name, s);
#endif

    if (fun->rettype == GRETL_TYPE_VOID) {
        if (s == NULL || *s == '\0') {
            ; /* plain "return" from void function: OK */
        } else if (line == NULL) {
            gretl_errmsg_sprintf(_("%s: non-null return value '%s' is not valid"),
                                 fun->name, s);
            err = E_TYPES;
        } else {
            gretl_errmsg_sprintf(_("%s, line %d: non-null return value '%s' is not valid"),
                                 fun->name, line->idx, s);
            err = E_TYPES;
        }
        return err; /* handled */
    }

    if (s == NULL || *s == '\0') {
        if (line == NULL) {
            gretl_errmsg_sprintf(_("%s: return value is missing"), fun->name);
        } else {
            gretl_errmsg_sprintf(_("%s, line %d: return value is missing"),
                                 fun->name, line->idx);
        }
        err = E_TYPES;
    } else {
        if (gretl_namechar_spn(s) == strlen(s)) {
            /* returning a named variable, simple */
            call->retname = gretl_strdup(s);
        } else {
            if (line != NULL && gencomp && !cmd_subst(state->cmd)) {
                /* allow trying compilation */
                err = generate_return_value2(call, state, dset, s, line);
            } else {
                /* play it safe */
                err = generate_return_value1(call, state, dset, s);
            }
        }
    }

 last_step:

    if (err) {
        /* FIXME case of line == NULL? */
        if (s == NULL && line != NULL) {
            s = line->s;
        }
        set_func_error_message(err, fun, state, s, line);
    } else if (call->retname == NULL) {
        call->retname = gretl_strdup("$retval");
    }

    return err;
}

int current_function_size (void)
{
    ufunc *u = currently_called_function();

    return (u != NULL)? u->n_lines : 0;
}

static gchar *return_line;

/* to be called when a "return" statement occurs within a
   loop that's being called by a function
*/

int set_function_should_return (const char *line)
{
    if (gretl_function_depth() > 0) {
        g_free(return_line);
        return_line = g_strdup(line);
        return 0;
    } else {
        gretl_errmsg_set(_("return: can only be used in a function"));
        return E_PARSE;
    }
}

static int get_return_line (ExecState *state)
{
    if (return_line == NULL) {
        return 0;
    } else {
        strcpy(state->line, return_line);
        g_free(return_line);
        return_line = NULL;
        return 1;
    }
}

#define ONE_WORKSPACE 1

/* Under "ONE_WORKSPACE" we use a shared command line (workspace) for
   function calls.
*/

static ExecState *make_func_exec_state (CMD *cmd,
					DATASET *dset,
                                        PRN *prn,
                                        int *err)
{
#if ONE_WORKSPACE
    static char *shared_line;
#endif
    char *line = NULL;
    int free_line = 0;
    ExecState *state;
    MODEL *model;

#if ONE_WORKSPACE
    if (shared_line == NULL) {
	line = shared_line = calloc(MAXLINE, 1);
    } else {
	line = shared_line;
    }
#else
    /* always allocate a distinct @line */
    line = calloc(MAXLINE, 1);
    free_line = 1;
#endif

    if (line == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    state = calloc(1, sizeof *state);
    model = allocate_working_model();

    if (state == NULL || model == NULL) {
        *err = E_ALLOC;
    } else {
        *err = gretl_cmd_init(cmd);
    }

    if (!*err) {
        gretl_exec_state_init(state, FUNCTION_EXEC, line, cmd, model, prn);
        if (dset != NULL) {
            if (dset->submask != NULL) {
                state->submask = copy_dataset_submask(dset, err);
            }
            state->padded = dset->padmask != NULL;
        }
        state->callback = func_exec_callback;
	state->free_line = free_line;
    }

    return state;
}

static void restore_dataset_for_caller(ExecState *state,
                                       DATASET *dset,
                                       int orig_t1,
                                       int orig_t2,
                                       int orig_n)
{
    if (complex_subsampled()) {
        if (state->submask == NULL) {
            /* we were not sub-sampled on entry */
            restore_full_sample(dset, NULL);
        } else if (submask_cmp(state->submask, dset->submask)) {
            /* we were sub-sampled differently on entry */
            gretlopt opt = state->padded ? OPT_B : OPT_NONE;

            restore_full_sample(dset, NULL);
            restrict_sample_from_mask(state->submask, dset, opt);
        }
    } else if (dset->n > orig_n) {
        /* some observations were added inside the function; note
           that this is not allowed if the dataset is subsampled
           on entry
        */
        dataset_drop_observations(dset, dset->n - orig_n);
    }

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;
}

#define do_if_check(c) (c == IF || c == ELIF || c == ELSE || c == ENDIF)
#define may_have_genr(c) (c == IF || c == ELIF || c == GENR)

int gretl_function_exec_full (fncall *call, int rtype, DATASET *dset,
                              void *ret, char **descrip,
                              series_table **stab, PRN *prn)
{
    ufunc *u = call->fun;
    ExecState *state = NULL;
    GENERATOR *genr = NULL;
    CMD cmd = {0};
    void *ptr = NULL;
    int orig_n = 0;
    int orig_t1 = 0;
    int orig_t2 = 0;
    int indent0 = 0, started = 0;
    int redir_level = 0;
    int retline = -1;
    int gencomp = 0;
    int n_saved = 0;
    int i, j, err = 0;

#if COMP_DEBUG || EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: starting %s (depth %d, recursing %d)\n",
            u->name, gretl_function_depth(), is_recursing(call));
#endif

    err = maybe_check_function_needs(dset, u);
    if (err) {
        maybe_destroy_fncall(&call);
        return err;
    }

    /* record incoming dataset dimensions */
    if (dset != NULL) {
        call->orig_v = dset->v;
        orig_n = dset->n;
        orig_t1 = dset->t1;
        orig_t2 = dset->t2;
    } else {
        call->orig_v = 0;
    }

    /* precaution aganst early error */
    indent0 = gretl_if_state_record();

    err = check_function_args(call, prn);
    if (!err) {
        err = allocate_function_args(call, dset);
    }
    if (err) {
        /* get out before allocating further storage */
        maybe_destroy_fncall(&call);
        return err;
    }

    state = make_func_exec_state(&cmd, dset, prn, &err);

#if EXEC_DEBUG
    fprintf(stderr, "after prepare_func_exec_state: err = %d\n", err);
#endif

    if (!err) {
        start_fncall(call, dset, prn);
        started = 1;
        redir_level = print_redirection_level(prn);
    }

    /* when should we try to compile genrs, loops? */
#if COMPILE_RECURSIVE
    /* as of 2024-05-08 this seems to be too risky */
    gencomp = gretl_iterating() && !get_loop_renaming();
#else
    /* add clause to cut out functions that recurse */
    gencomp = gretl_iterating() && !get_loop_renaming() &&
	!function_is_recursive(u);
#endif

    /* get function lines in sequence and check, parse, execute */

    for (i=0; i<u->n_lines && !err; i++) {
        fn_line *fline = &u->lines[i];
        int this_gencomp = gencomp && !line_no_comp(fline);

        if (gretl_echo_on()) {
            pprintf(prn, "? %s\n", fline->s);
        }
        if (ignore_line(fline)) {
            continue;
        }
        strcpy(state->line, fline->s); /* needed? */
        call->line = fline;

        if (this_gencomp && may_have_genr(fline->ci)) {
            genr = fnline_get_genr(call, fline);
        } else {
            genr = NULL;
        }

#if 0
        fprintf(stderr, " line %d: '%s' (no_comp %d, genr %p)\n",
                i, fline->s, line_no_comp(fline), (void *) genr);
#endif

        /* check and adjust the if-state */
        if (do_if_check(fline->ci)) {
            if (fline->ci == ELSE || fline->ci == ENDIF) {
                state->cmd->ci = fline->ci;
                flow_control(state, NULL, NULL);
            } else if (genr != NULL) {
                state->cmd->ci = fline->ci;
                flow_control(state, dset, &genr);
            } else {
                ptr = this_gencomp ? &genr : NULL;
                err = maybe_exec_line(state, dset, ptr);
                if (genr != NULL) {
                    fnline_set_genr(call, fline, genr);
                } else if (ptr != NULL) {
                    fline->flags |= LINE_NOCOMP;
                }
            }
            if (err) {
                goto err_next;
            }
            if (gretl_if_state_false() && fline->next > 0) {
                /* skip to next relevant statement */
                i = fline->next - 1;
            }
            continue;
        }

        if (fline->ci == LOOP) {
            if (fline->ptr != NULL) {
                state->loop = fline->ptr;
                err = gretl_loop_exec(state, dset);
                n_saved++;
            } else {
                /* assemble then execute the loop */
                if (this_gencomp && !gretl_function_recursing()) {
                    ptr = &fline->ptr;
                } else {
                    ptr = NULL;
                }
                for (j=i; j<=fline->next && !err; j++) {
                    strcpy(state->line, u->lines[j].s);
                    err = maybe_exec_line(state, dset, ptr);
                }
                if (!err) {
                    state->loop = fline->ptr;
		    err = gretl_loop_exec(state, dset);
                }
            }
            if (err) {
                set_func_error_message(err, u, state, state->line, NULL);
                break;
            } else if (get_return_line(state)) {
                /* a "return" statement was encountered in the loop */
                err = handle_return_statement(call, state, dset,
                                              this_gencomp, NULL);
                if (i < u->n_lines) {
                    retline = i;
                }
                break;
            }
            /* skip past the loop lines */
            i = fline->next;
            continue;
        } else if (fline->ci == FUNCRET) {
            err = handle_return_statement(call, state, dset, this_gencomp, fline);
            if (i < u->n_lines) {
                retline = i;
            }
            break;
        } else if (genr != NULL) {
            err = execute_genr(genr, dset, prn);
            n_saved++;
        } else {
            ptr = (this_gencomp && may_have_genr(fline->ci))? &genr : NULL;
            err = maybe_exec_line(state, dset, ptr);
            if (genr != NULL) {
                fnline_set_genr(call, fline, genr);
            } else if (ptr != NULL) {
                fline->flags |= LINE_NOCOMP;
            }
        }

    err_next:

        if (err) {
            set_func_error_message(err, u, state, state->line, fline);
            break;
        }
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec: %s: finished main exec, "
            "err = %d, dset->v = %d\n", u->name, err,
            (dset != NULL)? dset->v : 0);
#endif

    if (dset != NULL) {
        /* restore the sample that was in place on entry */
        restore_dataset_for_caller(state, dset, orig_t1, orig_t2, orig_n);
    }

    if (err || (retline >= 0)) {
        /* hit an error, or returned prior to the end of the function */
        gretl_if_state_reset(indent0);
    } else {
        err = gretl_if_state_check(indent0);
    }

    function_assign_returns(call, rtype, dset, ret, descrip, stab, prn, &err);

    if (gencomp || n_saved > 0) {
        reset_saved_uservars(call->fun, 0);
    }

    gretl_exec_state_clear(state);
    free(state);

    if (started) {
        int stoperr = stop_fncall(call, rtype, ret, dset, prn, redir_level);

        if (stoperr && !err) {
            err = stoperr;
        }
    } else {
        maybe_destroy_fncall(&call);
    }

#if EXEC_DEBUG
    fprintf(stderr, "gretl_function_exec (%s) finished: err = %d\n",
            u->name, err);
#endif
#if COMP_DEBUG
    if (err) {
        fprintf(stderr, "*** %s: exiting with err = %d ***\n", u->name, err);
    }
#endif

    return err;
}

int gretl_function_exec (fncall *call, int rtype, DATASET *dset,
                         void *ret, PRN *prn)
{
    return gretl_function_exec_full(call, rtype, dset, ret,
                                    NULL, NULL, prn);
}

/* look up name of supplied argument based on name of variable
   inside function */

char *gretl_func_get_arg_name (const char *argvar, int *err)
{
    fncall *call = current_function_call();
    char *ret = NULL;

    *err = E_DATA;

    if (call != NULL) {
        ufunc *u = call->fun;
        int i, n = call->argc;

        for (i=0; i<n; i++) {
            if (!strcmp(argvar, u->params[i].name)) {
                *err = 0;
                if (call->args[i].upname != NULL) {
                    ret = gretl_strdup(call->args[i].upname);
                } else {
                    ret = gretl_strdup("");
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
 * @name: name of object.
 * @vnum: ID number (specific to objects of type series).
 *
 * Checks whether the named object currently has 'const' status,
 * by virtue of its being made available as a const argument
 * to a user-defined function. Note that @name is the name by
 * which the object is known within the function, which will
 * likely differ from its name in the caller.
 *
 * Returns: non-zero if the object is const, 0 if it is not.
 */

int object_is_const (const char *name, int vnum)
{
    fncall *call = current_function_call();
    int ret = 0;

    if (call != NULL) {
        if (name != NULL) {
            fn_param *param;
            int i;

            for (i=0; i<call->fun->n_params; i++) {
                param = &call->fun->params[i];
                if (!strcmp(name, param->name)) {
                    ret = param_is_const(param);
                    break;
                }
            }
        }
        if (!ret && vnum > 0 && vnum < call->orig_v) {
            /* We're looking at a series that is not local to
               the called function, but not yet identified as
               read-only. It probably _should_ be read-only
               unless it was given in pointer form.
               Note: this check added 2018-10-18.
            */
            if (!in_gretl_list(call->ptrvars, vnum)) {
                ret = 1;
            }
        }
    }

    return ret;
}

/**
 * object_is_function_arg:
 * @name: name of object (e.g. matrix).
 *
 * Checks whether the named object has been made available
 * as a function argument.
 *
 * Returns: non-zero (the gretl type of the object) if the
 * object has arg status, 0 otherwise.
 */

int object_is_function_arg (const char *name)
{
    fncall *call = current_function_call();

    if (call != NULL) {
        fn_param *params = call->fun->params;
        int i;

        for (i=0; i<call->fun->n_params; i++) {
            if (!strcmp(name, params[i].name)) {
                return params[i].type;
            }
        }
    }

    return 0;
}

static int allow_full_data;

void allow_full_data_access (int s)
{
    if (s > 0) {
        allow_full_data = 1;
    } else if (s == 0) {
        allow_full_data = 0;
    }
}

/**
 * sample_range_get_extrema:
 * @dset: dataset info.
 * @t1: location to receive earliest possible starting
 * observation, or NULL.
 * @t2: location to receive latest possible ending observation,
 * or NULL.
 *
 * Fills out @t1 and/or @t2, making allowance for the possibility
 * that we're currently executing a function, on entry to
 * which the sample range was restricted: within the function,
 * we are not allowed to overstep the bounds set on entry.
 */

void sample_range_get_extrema (const DATASET *dset, int *t1, int *t2)
{
    int tvals[2] = {0, dset->n - 1};

    if (!allow_full_data) {
        fncall *call = current_function_call();

        if (call != NULL) {
            tvals[0] = call->obs.t1;
            tvals[1] = call->obs.t2 + call->obs.added;
            /* FIXME this needs to be smarter */
            if (tvals[1] > dset->n - 1) {
                tvals[1] = dset->n - 1;
            }
        }
    }

    if (t1 != NULL) {
        *t1 = tvals[0];
    }
    if (t2 != NULL) {
        *t2 = tvals[1];
    }
}

void extend_function_sample_range (int addobs)
{
    fncall *call = current_function_call();

    if (call != NULL) {
        call->obs.added += addobs;
    }
}

/**
 * series_is_accessible_in_function:
 * @ID: series ID number (0 = constant).
 *
 * Returns: 1 if the series with ID number @ID is accessible
 * by number at the current level of function execution,
 * otherwise 0.
 */

int series_is_accessible_in_function (int ID, const DATASET *dset)
{
    /* This is harder to get right than one might think, and
       for now (as of 2018-10-18) we will not attempt to
       screen access to series in this way. However, we do
       strive to ensure that series are treated as read-only
       when they ought to be. See geneval.c for (potential)
       uses of this check, currently def'd out.
    */
#if 1
    return 1; /* allow all */
#else
    fncall *fc = current_function_call();
    int ret = 1;

    if (fc != NULL) {
        /* assume not accessible without contrary evidence */
        ret = 0;
        if (ID == 0 || ID == LISTSEP) {
            /* the constant, always OK; also LISTSEP */
            ret = 1;
        } else if (ID >= fc->orig_v) {
            /* the series post-dates the start of execution */
            ret = 1;
        } else if (in_gretl_list(fc->listvars, ID)) {
            /* series was given in a list argument */
            ret = 1;
        } else if (in_gretl_list(fc->ptrvars, ID)) {
            /* series was given in "pointer" form */
            ret = 1;
        }
    }

    if (ret == 0) {
        gretl_errmsg_sprintf("Series %d (%s) is not accessible",
                             ID, dset->varname[ID]);
    }

    return ret;
#endif
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

static void show_gfn_help (fnpkg *pkg, int gui_help,
                           int markup, PRN *prn)
{
    const char *text;
    int mdconv;

    text = gui_help ? pkg->gui_help : pkg->help;
    mdconv = help_text_is_markdown(pkg, gui_help);

    if (mdconv) {
        /* help is in markdown */
        md_to_gretl(text, prn);
    } else {
        if (markup) {
            pputs(prn, "<@itl=\"Help text\">:\n\n");
            pputs(prn, "<mono>\n");
        } else {
            pputs(prn, "Help text:\n");
        }
        pputs(prn, text);
        if (markup) {
            pputs(prn, "\n</mono>");
        }
        pputs(prn, "\n\n");
    }
}

/* Generate help output for a packaged function, either for
   display on the console, or with markup for display in a
   GtkTextView window in the GUI (opt & OPT_M).
*/

static void real_user_function_help (ufunc *fun, gretlopt opt, PRN *prn)
{
    fnpkg *pkg = fun->pkg;
    int markup = (opt & OPT_M);
    int i;

    if (markup) {
        pprintf(prn, "<@itl=\"Function\">: %s\n", fun->name);
    } else {
        pprintf(prn, "Function: %s\n", fun->name);
    }

    if (pkg != NULL) {
        if (markup) {
            pprintf(prn, "<@itl=\"Package\">: %s %s (%s)\n", pkg->name, pkg->version,
                    pkg->date);
            pprintf(prn, "<@itl=\"Author\">: %s\n", pkg->author? pkg->author : "unknown");
            if (pkg->email != NULL && *pkg->email != '\0') {
                pprintf(prn, "<@itl=\"Email\">: %s\n", pkg->email);
            }
        } else {
            pprintf(prn, "Package: %s %s (%s)\n", pkg->name, pkg->version,
                    pkg->date);
            pprintf(prn, "Author: %s\n", pkg->author? pkg->author : "unknown");
            if (pkg->email != NULL && *pkg->email != '\0') {
                pprintf(prn, "Email:  %s\n", pkg->email);
            }
        }
        pputc(prn, '\n');
    }

    if ((opt & OPT_G) && pkg != NULL && pkg->gui_help != NULL) {
        /* GUI-specific help is preferred, and is available */
        show_gfn_help(pkg, 1, markup, prn);
        return;
    }

    if (markup) {
        pputs(prn, "<@itl=\"Parameters\">: ");
    } else {
        pputs(prn, "Parameters: ");
    }

    if (fun->n_params > 0) {
        pputc(prn, '\n');
        for (i=0; i<fun->n_params; i++) {
            pprintf(prn, " %s (%s",
                    fun->params[i].name, gretl_type_get_name(fun->params[i].type));
            if (fun->params[i].descrip != NULL) {
                pprintf(prn, ": %s)\n", fun->params[i].descrip);
            } else {
                pputs(prn, ")\n");
            }
        }
        pputc(prn, '\n');
    } else {
        pputs(prn, "none\n\n");
    }

    if (markup) {
        pputs(prn, "<@itl=\"Return value\">: ");
    } else {
        pputs(prn, "Return value: ");
    }

    if (fun->rettype != GRETL_TYPE_NONE && fun->rettype != GRETL_TYPE_VOID) {
        pprintf(prn, "%s\n\n", gretl_type_get_name(fun->rettype));
    } else {
        pputs(prn, "none\n\n");
    }

    if (pkg != NULL && pkg->help != NULL) {
        if (is_pdf_ref(pkg->help)) {
            gchar *pdfname = g_strdup(pkg->fname);
            gchar *p = strrchr(pdfname, '.');

            *p = '\0';
            strcat(p, ".pdf");
            if (markup) {
                pprintf(prn, "<@itl=\"Documentation\">: <@adb=\"%s\">\n\n", pdfname);
            } else {
                pprintf(prn, "See %s\n\n", pdfname);
            }
            g_free(pdfname);
        } else {
            show_gfn_help(pkg, 0, markup, prn);
        }
    }

    if (pkg != NULL && pkg->sample != NULL) {
        if (markup) {
            pputs(prn, "<@itl=\"Sample script\">:\n\n");
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

    if (pkg->help != NULL && !strncmp(pkg->help, "pdfdoc", 6)) {
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

/**
 * function_package_has_gui_help:
 * @pkg: function-package pointer.
 *
 * Checks whether @pkg has GUI-specific help text.
 *
 * Returns: 1 if so, 0 otherwise.
 */

int function_package_has_gui_help (fnpkg *pkg)
{
    return pkg->gui_help != NULL && !string_is_blank(pkg->gui_help);
}

void function_package_set_editor (fnpkg *pkg, void *editor)
{
    if (pkg != NULL) {
        pkg->editor = editor;
    }
}

void *function_package_get_editor (fnpkg *pkg)
{
    return pkg == NULL ? NULL : pkg->editor;
}

/**
 * delete_function_package:
 * @gfnname: full path to gfn file.
 *
 * Deletes @gfnname, taking care of deleting its enclosing
 * specific subdirectory and all of its contents if applicable.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int delete_function_package (const char *gfnpath)
{
    char *p = strrslash(gfnpath);
    gchar *pkgname = NULL;
    gchar *pkgdir = NULL;
    gchar *pkgsub = NULL;
    int err = 0;

    errno = 0;

    if (p != NULL) {
        /* @pkgname: strip extension from @gfnpath */
        pkgname = g_strdup(p + 1);
        p = strrchr(pkgname, '.');
        if (p != NULL) {
            *p = '\0';
        }
        /* @pkgdir: directory in which the gfn lives */
        pkgdir = g_strdup(gfnpath);
        p = strrslash(pkgdir);
        *p = '\0';
        p = strrslash(pkgdir);
        /* @pkgsub: subdir in which the gfn lives */
        if (p != NULL) {
            pkgsub = g_strdup(p + 1);
        }
    }

    if (pkgname != NULL && pkgdir != NULL && pkgsub != NULL &&
        !strcmp(pkgname, pkgsub)) {
        /* We should delete the tree @pkgdir only if the last
           directory in that path, @pkgsub, compares equal to
           the basename of the package, @pkgname. For example:
           foo(.gfn) lives in <path-to-functions>/foo/.
        */
        err = gretl_deltree(pkgdir);
        if (err) {
            gretl_errmsg_sprintf(_("Couldn't delete %s"), pkgdir);
        }
    } else {
        /* just delete the .gfn file itself */
        err = gretl_remove(gfnpath);

        if (err) {
            gretl_errmsg_sprintf(_("Couldn't delete %s"), gfnpath);
        }
    }

    if (err) {
        fprintf(stderr, "failure in delete_function_package: gfnpath '%s'\n",
                gfnpath);
        fprintf(stderr, " pkgname '%s', pkgdir '%s', pkgsub '%s'\n",
                pkgname, pkgdir, pkgsub);
    }

    g_free(pkgname);
    g_free(pkgdir);
    g_free(pkgsub);

    return err;
}

/**
 * uninstall_function_package:
 * @package: name of package (with or without .gfn extension).
 * @opt: may include OPT_P to "purge" the package.
 * @prn: gretl printer.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int uninstall_function_package (const char *package, gretlopt opt,
                                PRN *prn)
{
    gchar *gfnname = NULL;
    gchar *pkgname = NULL;
    char *p, fname[MAXLEN];
    int err;

    if (has_suffix(package, ".gfn")) {
        gfnname = g_strdup(package);
        pkgname = g_strdup(package);
        p = strrchr(pkgname, '.');
        *p = '\0';
    } else {
        gfnname = g_strdup_printf("%s.gfn", package);
        pkgname = g_strdup(package);
    }

    *fname = '\0';
    err = get_full_filename(gfnname, fname, OPT_I);

    if (!err && !gretl_file_exists(fname)) {
        gretl_errmsg_sprintf(_("Couldn't find %s"), gfnname);
        err = E_FOPEN;
    }

    if (!err) {
        function_package_unload_full_by_filename(fname);
        if (opt & OPT_P) {
            err = delete_function_package(fname);
        }
    }

    if (!err && gretl_messages_on()) {
        if (opt & OPT_P) {
            pprintf(prn, _("Removed %s\n"), pkgname);
        } else {
            pprintf(prn, _("Unloaded %s\n"), pkgname);
        }
    }

    g_free(gfnname);
    g_free(pkgname);

    return err;
}

/* full apparatus for hansl indentation follows */

static void get_cmdword (const char *s, char *word)
{
    if (!strncmp(s, "catch ", 6)) {
        s += 6;
    }

    if (sscanf(s, "%*s <- %8s", word) != 1) {
        sscanf(s, "%8s", word);
    }
}

#define bare_quote(p,s) (*p == '"' && (p-s==0 || *(p-1) != '\\'))
#define starts_ccmt(p)  (*p == '/' && *(p+1) == '*')
#define ends_ccmt(p)    (*p == '*' && *(p+1) == '/')

static void check_for_comment (const char *s, int *incomm)
{
    const char *p = s;
    int commbak = *incomm;
    int quoted = 0;

    while (*p) {
        if (!quoted && !*incomm && *p == '#') {
            break;
        }
        if (!*incomm && bare_quote(p, s)) {
            quoted = !quoted;
        }
        if (!quoted) {
            if (starts_ccmt(p)) {
                *incomm = 1;
                p += 2;
            } else if (ends_ccmt(p)) {
                *incomm = 0;
                p += 2;
                p += strspn(p, " ");
            }
        }
        if (*p) {
            p++;
        }
    }

    if (*incomm && commbak) {
        /* on the second or subsequent line of a C-style comment */
        *incomm = 2;
    }
}

/* determine whether a given line is subject to continuation,
   i.e. ends with backslash or comma, but not in a comment
*/

static int line_broken (const char *s)
{
    int ret = 0;

    if (*s != '\0') {
        int i, n = strlen(s);

        for (i=n-1; i>=0; i--) {
            if (s[i] == '\\' || s[i] == ',') {
                ret = 1;
            } else if (!ret && !isspace(s[i])) {
                break;
            } else if (ret && s[i] == '#') {
                ret = 0;
                break;
            }
        }
    }

    return ret;
}

static void strip_trailing_whitespace (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>=0; i--) {
        if (s[i] == '\n' || s[i] == ' ' || s[i] == '\t') {
            s[i] = '\0';
        } else {
            break;
        }
    }

    strcat(s, "\n");
}

/* determine position of unmatched left parenthesis,
   when applicable, if @s starts a function definition
*/

static int left_paren_offset (const char *s)
{
    const char *p = strchr(s, '(');

    if (p != NULL && strchr(p, ')') == NULL) {
        return p - s;
    } else {
        return 0;
    }
}

static int wordmatch (const char *s, const char *test)
{
    int n = strlen(test);

    return (!strncmp(s, test, n) && (s[n] == '\0' || isspace(s[n])));
}

void adjust_indent (const char *s, int *this_indent, int *next_indent)
{
    const char *block_starts[] = {
        "loop", "if", "nls", "mle", "gmm", "mpi", "plot",
        "function", "restrict", "system", "foreign",
        "outfile", "gpbuild", NULL
    };
    int ti = *next_indent;
    int ni = *next_indent;
    int i, matched = 0;

    if (*s == '\0') {
        *this_indent = *next_indent;
        return;
    }

    /* skip "catch" if present */
    if (!strncmp(s, "catch ", 6)) {
        s += 6;
        s += strspn(s, " ");
    }

    for (i=0; block_starts[i] != NULL && !matched; i++) {
        if (wordmatch(s, block_starts[i])) {
            matched = 1;
            ni++;
        }
    }

    if (!matched) {
        if (wordmatch(s, "end") ||
            wordmatch(s, "endif") ||
            wordmatch(s, "endloop")) {
            ti--;
            ni--;
        } else if (wordmatch(s, "else") ||
                   wordmatch(s, "elif")) {
            ni = ti;
            ti--;
        }
    }

    *this_indent = ti;
    *next_indent = ni;
}

void normalize_hansl (const char *buf, int tabwidth, PRN *prn)
{
    char *line;
    char word[9];
    const char *ins;
    size_t llen = 1024;
    int this_indent = 0;
    int next_indent = 0;
    int continuation = 0;
    int incomment = 0;
    int inforeign = 0;
    int lp_pos = 0;
    int lp_zero = 0;
    int nsp;

    line = malloc(llen);
    line[0] = '\0';
    bufgets_init(buf);

    while (safe_bufgets(&line, &llen, buf)) {
        int handled = 0;
	int lbreak = 0;

        strip_trailing_whitespace(line);

        if (string_is_blank(line)) {
            pputc(prn, '\n');
            continue;
        }

        check_for_comment(line, &incomment);
        ins = line + strspn(line, " \t");

        if (!incomment) {
	    lbreak = line_broken(line);
            *word = '\0';
            get_cmdword(ins, word);
            if (!strcmp(word, "foreign")) {
                inforeign = 1;
            } else if (inforeign) {
                if (!strncmp(ins, "end foreign", 11)) {
                    inforeign = 0;
                } else {
                    /* insert foreign line as is */
                    pputs(prn, line);
                    handled = 1;
                }
            } else {
                if (!strcmp(word, "function")) {
		    /* record position of left parenthesis */
                    lp_pos = left_paren_offset(ins);
                } else if (lp_pos > 0 && strchr(ins, ')') != NULL) {
		    /* flag end of function signature */
                    lp_zero = 1;
                }
		adjust_indent(word, &this_indent, &next_indent);
            }
        }

        if (!handled) {
            nsp = this_indent * tabwidth;
            if (incomment == 2) {
		/* C-style comment continuation */
		nsp += 3;
            } else if (continuation) {
                if (lp_pos > 0) {
		    /* in a function signature */
                    nsp = lp_pos + 1;
                } else {
		    /* small offset for continuation */
                    nsp += 2;
                }
            }
            /* insert required number of spaces */
            bufspace(nsp, prn);
            /* insert line, with leading space pruned */
            pputs(prn, ins);
        }

	/* set parameters for next line */
	this_indent = next_indent;
	continuation = lbreak;
        if (lp_zero) {
	    /* scrub the left-paren record */
            lp_zero = lp_pos = 0;
        }
    }

    bufgets_finalize(buf);
    free(line);
}
