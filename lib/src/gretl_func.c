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
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"
#include "gretl_xml.h"

#define CALLSTACK_DEPTH 8
#define FN_NAMELEN 32
#define FN_DEBUG 0

typedef struct fn_param_ fn_param;
typedef struct fn_return_ fn_return;
typedef struct fncall_ fncall;
typedef struct fnpkg_ fnpkg;

enum {
    FN_PARAMS,
    FN_RETURNS
};

struct fn_param_ {
    char *name;
    char type;
    double deflt;
    double min;
    double max;
};

struct fn_return_ {
    char *name;
    char type;
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
    int n_returns;
    fn_return *returns;
};

struct fncall_ {
    ufunc *fun;
    int lnum;
    int argc;
    char **argv;
    char **assv;
    int *asslist;
};

struct fnpkg_ {
    int ID;
    FuncDataReq dreq;
    char *fname;
    char *author;
    char *version;
    char *date;
    char *descrip;
};

static int n_ufuns;
static ufunc **ufuns;
static ufunc *current_ufun;

static int n_pkgs;
static fnpkg **pkgs;

static fncall **callstack;

static void free_fncall (fncall *call);
static int real_add_fn_line (ufunc *fun, const char *s);
static void real_user_function_help (ufunc *fun, fnpkg *pkg, PRN *prn);

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

int gretl_executing_function (void)
{
    return fn_executing;
}

static void set_executing_on (fncall *call)
{
    fn_executing++;
#if FN_DEBUG
    fprintf(stderr, "set_executing_on: fun=%s, fn_executing=%d\n", 
	    call->fun->name, fn_executing);
#endif
}

static void set_executing_off (fncall *call)
{
    fn_executing--;
#if FN_DEBUG
    fprintf(stderr, "set_executing_off: fun=%s, fn_executing=%d\n",
	    call->fun->name, fn_executing);
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

int user_func_get_return_types (const ufunc *fun,
				int *n_returns,
				char **return_types)
{
    char *rtypes = NULL;
    int i, nr = 0;
    int err = 0;

    if (fun == NULL) {
	return E_DATA;
    }

    nr = fun->n_returns;
    *n_returns = nr;

    if (nr > 0) {
	rtypes = malloc(nr);
	if (rtypes == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<nr; i++) {
		rtypes[i] = fun->returns[i].type;
	    }
	    *return_types = rtypes;
	}
    } else {
	*return_types = NULL;
    }
	
    return err;
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

/* function call stack mechanism */

static int callstack_init (void)
{
    int i, err = 0;

    if (callstack != NULL) {
	return 0;
    }

    callstack = malloc(CALLSTACK_DEPTH * sizeof *callstack);
    if (callstack == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<CALLSTACK_DEPTH; i++) {
	    callstack[i] = NULL;
	}
    }

    return err;
}

static void callstack_destroy (void)
{
    int i;

    if (callstack != NULL) {
	for (i=0; i<CALLSTACK_DEPTH; i++) {
	    if (callstack[i] != NULL) {
		free_fncall(callstack[i]);
	    }
	}
    }

    free(callstack);
    callstack = NULL;
}

int gretl_function_stack_depth (void)
{
    int i, n = 0;

    if (callstack == NULL) {
	callstack_init();
    }

    if (callstack != NULL) {
	for (i=0; i<CALLSTACK_DEPTH; i++) {
	    if (callstack[i] != NULL) n++;
	    else break;
	}
    }

    return n;
}

static int push_fncall (fncall *call, const DATAINFO *pdinfo)
{
    int i, nc;

    if (callstack == NULL && callstack_init()) {
	return E_ALLOC;
    }

    nc = gretl_function_stack_depth();
    if (nc == CALLSTACK_DEPTH) {
	strcpy(gretl_errmsg, "Function call stack depth exceeded");
	return 1;
    }

    for (i=nc; i>0; i--) {
	callstack[i] = callstack[i-1];
    }	

    callstack[0] = call;

    set_executing_on(call);
    push_program_state(pdinfo);

    return 0;
}

static void copy_values_to_assignee (int targ, int src, double **Z, 
				     const DATAINFO *pdinfo)
{
    int t, n = (var_is_series(pdinfo, targ))? pdinfo->n : 1;

#if FN_DEBUG
    fprintf(stderr, "copy_values_to_assignee: copying from Z[%d] (%s) to Z[%d] (%s)\n",
	    src, pdinfo->varname[src], targ, pdinfo->varname[targ]);
#endif

    for (t=0; t<n; t++) {
	Z[targ][t] = Z[src][t];
    }
}

static int new_var_assignment (fncall *call, double ***pZ, DATAINFO *pdinfo, 
			       int nc, int *locals)
{
    ufunc *cfun = call->fun;
    int i, j, src, targ;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "\n*** start assignment of return values ***\n");
#endif

    if (call->asslist == NULL) {
#if FN_DEBUG
	fprintf(stderr, "function call: no assigments\n");
#endif
	return 0;
    }

#if FN_DEBUG
    printlist(call->asslist, "call->asslist");
    fprintf(stderr, "number of assignments = %d, number of returns = %d\n",
	    call->asslist[0], cfun->n_returns);
#endif

    if (cfun->n_returns < call->asslist[0]) {
	fprintf(stderr, "bad function call: assignments > returns\n");
	return 1;
    }

    for (i=0; i<call->asslist[0]; i++) {
	targ = call->asslist[i+1];
	src = -1;

#if FN_DEBUG
	fprintf(stderr, " assv[%d] = '%s'\n", i, call->assv[i]);
	fprintf(stderr, " return[%d] = '%s'\n", i, cfun->returns[i].name);
#endif
	if (call->assv[i] == NULL || *call->assv[i] == '\0') {
	    fprintf(stderr, " got a null assignment\n");
	    continue;
	}

	/* find the source variable */
	for (j=1; j<pdinfo->v; j++) {
	    if (STACK_LEVEL(pdinfo, j) != nc) {
		continue;
	    }
	    if (!strcmp(pdinfo->varname[j], cfun->returns[i].name)) {
#if FN_DEBUG
		fprintf(stderr, " identified source as var %d (%s)\n",
			j, pdinfo->varname[j]);
#endif
		src = j;
		break;
	    }
	}

	if (src > 0) {
	    if (targ > 0) {
		/* copy values to pre-existing var at caller level */
		copy_values_to_assignee(targ, src, *pZ, pdinfo);
	    } else if (locals != NULL) {
		/* rename variable as caller desired and mark as global */
		int pos;

		strcpy(pdinfo->varname[src], call->assv[i]);
		STACK_LEVEL(pdinfo, src) -= 1; 
		if ((pos = in_gretl_list(locals, src))) {
		    gretl_list_delete_at_pos(locals, pos);
		}
	    }
	} else {
	    /* assigning a matrix, not a variable? */
	    gretl_matrix *S = get_matrix_by_name(cfun->returns[i].name,
						 pdinfo);


	    if (S != NULL) {
		gretl_matrix *T;
#if FN_DEBUG
		fprintf(stderr, " identified source as matrix '%s' (%p)\n",
			cfun->returns[i].name, (void *) S);
#endif
		T = get_matrix_by_name_at_level(call->assv[i], nc - 1, pdinfo); 
		if (T != NULL) {
#if FN_DEBUG
		    fprintf(stderr, " identified target as matrix '%s' (%p)\n",
			    call->assv[i], (void *) T);
#endif
		    /* copy values to pre-existing matrix at caller level */
		    err = gretl_matrix_copy_values(T, S);
		} else {
#if FN_DEBUG
		    fprintf(stderr, " naming matrix '%s' (%p) as '%s' (level = %d)\n",
			    cfun->returns[i].name, (void *) S, call->assv[i], nc-1);
#endif
		    /* rename matrix as caller desired and mark as global */
		    err = user_matrix_set_name_and_level(S, call->assv[i], nc - 1);
		}
	    }
	}
    }

#if FN_DEBUG
    fprintf(stderr, "*** done assignment of returns ***\n\n");
#endif

    return err;
}

static int *make_locals_list (const DATAINFO *pdinfo, int nc, int *err)
{
    int *locals = NULL;
    int i, nlocal = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (STACK_LEVEL(pdinfo, i) == nc) {
	    nlocal++;
	}
    }

    if (nlocal > 0) {
	locals = gretl_list_new(nlocal);
	if (locals == NULL) {
	    *err = E_ALLOC;
	} else {
	    int j = 1;

	    for (i=1; i<pdinfo->v; i++) {
		if (STACK_LEVEL(pdinfo, i) == nc) {
		    locals[j++] = i;
		}
	    }
	}
    }

    return locals;
}

static int fn_offers_matrix (fncall *call)
{
    int i;

    for (i=0; i< call->fun->n_returns; i++) {
	if (call->fun->returns[i].type == ARG_MATRIX) {
	    return 1;
	}
    }

    return 0;
}

static int 
destroy_or_assign_local_vars (fncall *call, double ***pZ, DATAINFO *pdinfo, 
			      int nc)
{
    int *locals = NULL;
    int anyerr = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "destroy_or_assign_local_vars: depth = %d\n", nc);
#endif

    locals = make_locals_list(pdinfo, nc, &err);

    if (err || (locals == NULL && !fn_offers_matrix(call))) {
#if FN_DEBUG
	fprintf(stderr, "locals = %p, err = %d, returning\n", 
		(void *) locals, err);
#endif
	destroy_saved_lists_at_level(nc);
	return err;
    }    

    err = new_var_assignment(call, pZ, pdinfo, nc, locals);

    if (locals != NULL) {
	anyerr = dataset_drop_listed_variables(locals, pZ, pdinfo, NULL);
	if (anyerr && !err) {
	    err = anyerr;
#if FN_DEBUG
	    fprintf(stderr, "dataset_drop_listed_variables: err = %d\n", err);
#endif
	}
	free(locals);
    }

    anyerr = destroy_saved_lists_at_level(nc);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_saved_lists_at_level(%d): err = %d\n", nc, err);
#endif
    }

    anyerr = destroy_user_matrices_at_level(nc);
    if (anyerr && !err) {
	err = anyerr;
#if FN_DEBUG
	fprintf(stderr, "destroy_user_matrices_at_level(%d): err = %d\n", nc, err);
#endif
    }    

#if FN_DEBUG
    fprintf(stderr, "destroy_or_assign_local_vars: returning err = %d\n", err);
#endif

    return err;
}

static int unstack_fncall (double ***pZ, DATAINFO **ppdinfo)
{
    fncall *call;   
    int i, nc;
    int err = 0;

    if (callstack == NULL) {
	return 1;
    }
    
    call = callstack[0];

    nc = gretl_function_stack_depth();

#if FN_DEBUG
    fprintf(stderr, "unstack_fncall: terminating call to "
	    "function '%s' at depth %d\n", 
	    call->fun->name, nc);
#endif

    pop_program_state(pZ, ppdinfo);

    err = destroy_or_assign_local_vars(call, pZ, *ppdinfo, nc);

    set_executing_off(call);
    free_fncall(call);

    for (i=0; i<nc; i++) {
	if (i == nc - 1) {
	    callstack[i] = NULL;
	} else {
	    callstack[i] = callstack[i+1];
	}
    }

    return err;
}

static int function_is_on_stack (ufunc *func)
{
    int i;

    if (callstack != NULL) {
	for (i=0; i<CALLSTACK_DEPTH; i++) {
	    if (callstack[i] == NULL) {
		break;
	    }
	    if ((callstack[i])->fun == func) {
		return 1;
	    }
	}
    }

    return 0;
}

static fncall *current_call (void)
{
    if (callstack == NULL) {
	return NULL;
    } else {
	return callstack[0];
    }
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
	params[i].deflt = NADBL;
	params[i].min = NADBL;
	params[i].max = NADBL;
    }

    return params;
}

static fn_return *allocate_returns (int n)
{
    fn_return *returns;
    int i;

    returns = malloc(n * sizeof *returns);
    if (returns == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	returns[i].name = NULL;
	returns[i].type = 0;
    }

    return returns;
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

    fun->n_returns = 0;
    fun->returns = NULL;

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

static void free_returns_array (fn_return *returns, int n)
{
    int i;

    if (returns == NULL) return;

    for (i=0; i<n; i++) {
	free(returns[i].name);
    }
    free(returns);
}

static void clear_ufunc_data (ufunc *fun)
{
    free(fun->help);
    free_strings_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    free_returns_array(fun->returns, fun->n_returns);
    
    fun->help = NULL;
    fun->lines = NULL;
    fun->params = NULL;
    fun->returns = NULL;

    fun->n_lines = 0;
    fun->n_returns = 0;
    fun->n_params = 0;
}

static void free_ufunc (ufunc *fun)
{
    free(fun->help);
    free_strings_array(fun->lines, fun->n_lines);
    free_params_array(fun->params, fun->n_params);
    free_returns_array(fun->returns, fun->n_returns);

    free(fun);
}

static fncall *fncall_new (ufunc *fun, int argc, char **argv,
			   int *asslist, char **assv)
{
    fncall *call = malloc(sizeof *call);

    if (call == NULL) {
	free_strings_array(argv, argc);
	free_strings_array(assv, asslist[0]);
	free(asslist);
	return NULL;
    }

    call->fun = fun;
    call->lnum = 0;

    call->argc = argc;
    call->argv = argv;

    call->asslist = asslist;
    call->assv = assv;

    return call;
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
    } else {
	return "unknown";
    }
}

static int field_to_type (const char *s)
{
#if FN_DEBUG
    fprintf(stderr, "field_to_type: looking at '%s'\n", s);
#endif

    if (isdigit(*s)) {
	return atoi(s);
    }

    if (!strncmp(s, "bool", 4)) return ARG_BOOL;
    if (!strcmp(s, "int"))      return ARG_INT;
    if (!strcmp(s, "scalar"))   return ARG_SCALAR;
    if (!strcmp(s, "series"))   return ARG_SERIES;
    if (!strcmp(s, "list"))     return ARG_LIST;
    if (!strcmp(s, "matrix"))   return ARG_MATRIX;

    return 0;
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

static int func_read_returns (xmlNodePtr node, ufunc *fun)
{
    xmlNodePtr cur;
    char *field;
    int n, err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n) || n < 0) {
	fprintf(stderr, "Couldn't read return values count\n");
	return E_DATA;
    }  

    if (n == 0) {
	return 0;
    }

    fun->returns = allocate_returns(n);
    if (fun->returns == NULL) {
	return E_ALLOC;
    } 

    fun->n_returns = n;

    cur = node->xmlChildrenNode;
    n = 0;
    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "return")) {
	    if (gretl_xml_get_prop_as_string(cur, "name", &field)) {
		fun->returns[n].name = field;
	    } else {
		fprintf(stderr, "return %d: couldn't get name\n", n);
		err = E_DATA;
	    }
	    if (gretl_xml_get_prop_as_string(cur, "type", &field)) {
		fun->returns[n].type = field_to_type(field);
		free(field);
	    } else {
		fprintf(stderr, "return %d: couldn't get type\n", n);
		err = E_DATA;
	    }
	    n++;
	}
	cur = cur->next;
    }

    if (!err && n != fun->n_returns) {
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
	    err = real_add_fn_line(fun, s);
	}
    }

    free(buf);

    return err;
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
    const char *typestr;
    int i;

    gretl_push_c_numeric_locale();

    pprintf(prn, "function %s ", fun->name);
    for (i=0; i<fun->n_params; i++) {
	if (i == 0) {
	    pputc(prn, '(');
	}
	typestr = arg_type_string(fun->params[i].type);
	pprintf(prn, "%s %s", typestr, fun->params[i].name);
	if (fun->params[i].type == ARG_BOOL) {
	    if (!na(fun->params[i].deflt)) {
		pprintf(prn, "[%g]", fun->params[i].deflt);
	    }
	} else if (scalar_arg(fun->params[i].type)) {
	    print_deflt_min_max(&fun->params[i], prn);
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
    int i;

    for (i=0; i<fun->n_returns; i++) {
	if (i == 0) {
	    pputs(prn, "  return ");
	}
	typestr = arg_type_string(fun->returns[i].type);
	pprintf(prn, "%s %s", typestr, fun->returns[i].name);
	if (i < fun->n_returns - 1) {
	    pputs(prn, ", ");
	} else {
	    pputc(prn, '\n');
	}
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

#if FN_DEBUG
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
	} else if (!xmlStrcmp(cur->name, (XUC) "returns")) {
	    err = func_read_returns(cur, fun);
	    if (err) {
		fprintf(stderr, "%s: error parsing function return values\n",
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
	    fputs("/>\n", fp);
	}
	fputs(" </params>\n", fp);

	gretl_pop_c_numeric_locale();
    }

    if (fun->n_returns > 0) {
	fprintf(fp, " <returns count=\"%d\">\n", fun->n_returns);
	for (i=0; i<fun->n_returns; i++) {
	    fprintf(fp, "  <return name=\"%s\" type=\"%s\"/>\n",
		    fun->returns[i].name, 
		    arg_type_string(fun->returns[i].type));
	}
	fputs(" </returns>\n", fp);
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

int write_function_package (const char *fname,
			    const int *privlist, 
			    const int *publist,
			    const char *author,
			    const char *version,
			    const char *date,
			    const char *descrip,
			    FuncDataReq dreq)
{
    char *pkgname;
    FILE *fp;
    int i, fi;

    if (n_ufuns == 0) {
	fprintf(stderr, "No functions are defined\n");
	return 0;
    }

    if (author == NULL || version == NULL || date == NULL || 
	descrip == NULL || publist == NULL || publist[0] == 0) {
	strcpy(gretl_errmsg, "Function information is incomplete");
	return E_DATA;
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
	fputs(" needs-time-series=\"true\"", fp);
    } else if (dreq == FN_NEEDS_QM) {
	fputs(" needs-qm=\"true\"", fp);
    } else if (dreq == FN_NEEDS_PANEL) {
	fputs(" needs-panel=\"true\"", fp);
    }

    fputs(">\n", fp);

    gretl_xml_put_tagged_string("author", author, fp);
    gretl_xml_put_tagged_string("version", version, fp);
    gretl_xml_put_tagged_string("date", date, fp);
    gretl_xml_put_tagged_string("description", descrip, fp);

    if (privlist != NULL) {
	for (i=1; i<=privlist[0]; i++) {
	    fi = privlist[i];
	    if (fi >= 0 && fi < n_ufuns) {
		write_function_xml(ufuns[fi], fp);
	    }
	}
    }

    for (i=1; i<=publist[0]; i++) {
	fi = publist[i];
	if (fi >= 0 && fi < n_ufuns) {
	    write_function_xml(ufuns[fi], fp);
	}
    }
	    
    fputs("</gretl-function-package>\n", fp);
    fputs("</gretl-functions>\n", fp);

    fclose(fp);

    return 0;
}

int function_package_get_info (const char *fname,
			       int **privlist, 
			       int **publist,
			       char **author,
			       char **version,
			       char **date,
			       char **descrip,
			       FuncDataReq *dreq)
{
    fnpkg *pkg = NULL;
    int npriv = 0, npub = 0;
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

    for (i=0; i<n_ufuns; i++) {
	if (ufuns[i]->pkgID == pkg->ID) {
	    if (ufuns[i]->private) {
		npriv++;
	    } else {
		npub++;
	    }
	}
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

    if (!err && publist != NULL && npub > 0) {
	list = gretl_list_new(npub);
	if (list != NULL) {
	    j = 1;
	    for (i=0; i<n_ufuns; i++) {
		if (ufuns[i]->pkgID == pkg->ID && !ufuns[i]->private) {
		    list[j++] = i;
		}
	    }
	    *publist = list;
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

static int 
read_user_function_package (xmlDocPtr doc, xmlNodePtr node, 
			    const char *fname, int task, 
			    PRN *prn, char **pname)
{
    xmlNodePtr cur;
    fnpkg *pkg;
    int err = 0;

    pkg = function_package_new(fname);
    if (pkg == NULL) {
	return E_ALLOC;
    }

    if (pname != NULL) {
	gretl_xml_get_prop_as_string(node, "name", pname);
    }

    if (gretl_xml_get_prop_as_bool(node, "needs-time-series")) {
	pkg->dreq = FN_NEEDS_TS;
    } else if (gretl_xml_get_prop_as_bool(node, "needs-qm")) {
	pkg->dreq = FN_NEEDS_QM;
    } else if (gretl_xml_get_prop_as_bool(node, "needs-panel")) {
	pkg->dreq = FN_NEEDS_PANEL;
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

#if FN_DEBUG
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

#if FN_DEBUG
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

static void free_fncall (fncall *call)
{
    free_strings_array(call->argv, call->argc);
    free_strings_array(call->assv, call->asslist[0]);
    free(call->asslist);

    free(call);
}

static ufunc *get_ufunc_by_name (const char *name)
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
	fprintf(stderr, "get_ufunc_by_name: name = '%s' (n_ufuns = %d);"
		" found match\n", name, n_ufuns);
    }
#endif

    return fun;
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
	if (get_ufunc_by_name(name) != NULL) {
	    ret = 1;
	}
    }

    return ret;
}

int gretl_is_public_user_function (const char *name)
{
    ufunc *fun = get_ufunc_by_name(name);

    if (fun != NULL && !fun->private) {
	return 1;
    } else {
	return 0;
    }
}

int is_user_matrix_function (const char *word)
{
    ufunc *fun = get_ufunc_by_name(word);

    if (fun != NULL) {
	if (fun->n_returns == 1 && fun->returns[0].type == ARG_MATRIX) {
	    return 1;
	}
    }

    return 0;
}

int gretl_get_user_function (const char *line, char **fnname)
{
    int ret = 0;

    if (n_ufuns > 0 && !string_is_blank(line)) {
	char name[FN_NAMELEN];

#if FN_DEBUG > 1
	fprintf(stderr, "gretl_is_user_function: testing '%s'\n", line);
#endif
	function_name_from_line(line, name);
	if (get_ufunc_by_name(name) != NULL) {
	    free(*fnname);
	    *fnname = gretl_strdup(name);
	    if (*fnname != NULL) {
		ret = 1;
	    }
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

/* count fields that may be separated by spaces or commas,
   including blank comma-delimited fields */

static int special_count_fields (char *s)
{
    int n, nf = 0;

    while (*s) {
	while (*s == ' ') s++;
	n = strcspn(s, " ,");
	if (n > 0) {
	    nf++;
	    s += n;
	    while (*s == ' ') s++;
	    if (*s == ',') s++;
	} else if (*s == ',') {
	    nf++;
	    s++;
	}
	while (*s == ' ') s++;
    }

    return nf;
}

static char *get_next_arg (char **parg, char *s, int maxlen)
{
    int m, n;

    while (*s == ' ') s++;

    n = strcspn(s, " ,");
    
    if (n > 0) {
	m = (n > maxlen)? maxlen : n;
	*parg = gretl_strndup(s, m);
	s += n;
	while (*s == ' ') s++;
	if (*s == ',') s++;
    } else if (*s == ',') {
	*parg = gretl_strdup("");
	s++;
    }

    while (*s == ' ') s++;

    return s;
}

static int arg_type_from_string (const char *s)
{
    int ret = 0;

    if (!strncmp(s, "bool", 4)) {
	ret = ARG_BOOL;
    } else if (!strcmp(s, "int")) {
	ret = ARG_INT;
    } else if (!strcmp(s, "scalar")) {
	ret = ARG_SCALAR;
    } else if (!strcmp(s, "series")) {
	ret = ARG_SERIES;
    } else if (!strcmp(s, "list")) {
	ret = ARG_LIST;
    } else if (!strcmp(s, "matrix")) {
	ret = ARG_MATRIX;
    }

    return ret;
}

/* Parse line and return an allocated array of strings consisting of
   the space- or comma-separated fields in line, each one truncated if
   necessary to a maximum of maxlen characters.
*/ 

static char **
get_separated_fields (const char *line, int *nfields, int maxlen, 
		      char **passtypes, int *err)
{
    char **fields = NULL;
    char *asstypes = NULL;
    char *cpy, *s;
    int realnf;
    int i, n, nf;

    *nfields = 0;

    if (string_is_blank(line)) {
	*err = E_PARSE;
	return NULL;
    }

    cpy = gretl_strdup(line);
    if (cpy == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* clean up the copied string */
    s = cpy;
    s += strspn(s, " ");
    if (*s == '(') {
	s++;
    }
    tailstrip(s);
    n = strlen(s);
    if (s[n-1] == ')') {
	s[n-1] = '\0';
    } 

    realnf = nf = special_count_fields(s);

    if (nf > 0) {
	if (passtypes != NULL) {
	    asstypes = calloc(nf, 1);
	    if (asstypes == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (!*err) {
	    fields = strings_array_new(nf);
	    if (fields == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (!*err) {
	    int atype = 0;
	    int j = 0;

	    for (i=0; i<nf && !*err; i++) {
		if (asstypes != NULL && atype) {
		    asstypes[j] = atype;
		}
		s = get_next_arg(&fields[j], s, maxlen);
		if (fields[j] == NULL) {
		    *err = E_ALLOC;
		} else if ((atype = arg_type_from_string(fields[j]))) {
		    if (asstypes == NULL || i == nf - 1) {
			*err = E_PARSE;
		    } else {
			free(fields[j]);
			fields[j] = NULL;
		    } 
		} else {
		    j++;
		}
	    }
	    realnf = j;
	}
    }

    if (!*err && realnf < nf) {
	if (realnf == 0) {
	    *err = E_PARSE;
	} else {
	    char **realf = realloc(fields, realnf * sizeof *fields);

	    if (realf == NULL) {
		*err = E_ALLOC;
	    } else {
		fields = realf;
		nf = realnf;
	    }
	}
    }

    if (*err) {
	free_strings_array(fields, nf);
	fields = NULL;
	free(asstypes);
    } else {
	*nfields = nf;
	if (passtypes != NULL) {
	    *passtypes = asstypes;
	}
    }

    free(cpy);

    return fields;
}

static int 
parse_function_args_etc (const char *s, int *argc, char ***pargv,
			 int **asslist, char ***passv, char **passtypes)
{
    char **argv = NULL;
    char **assign = NULL;
    const char *p;
    int na = 0, err = 0;

    *argc = 0;
    *pargv = NULL;
    *passv = NULL;

    if ((p = strchr(s, '=')) != NULL && *s != '=') {
	char *astr = gretl_strndup(s, p - s);

	/* seems we have a left-hand side assignment */
	assign = get_separated_fields(astr, &na, VNAMELEN - 1, 
				      passtypes, &err);
	free(astr);
	s = p + 1;
    } else {
	s++;
    }

    if (err) {
	return err;
    }

    if (na > 0) {
	*asslist = gretl_list_new(na);
    } else {
	*asslist = gretl_null_list();
    }

    if (asslist == NULL) {
	err = E_ALLOC;
	free_strings_array(assign, na);
    } else {
	*passv = assign;
    }

    if (!err) {
	/* skip over function name and spaces before args,
	   or deal with the case where the args are parenthesized
	*/
	const char *p = strchr(s, '(');

	if (p != NULL) {
	    s = p;
	} else {
	    s += strspn(s, " ");
	    s += strcspn(s, " ");
	    s += strspn(s, " ");
	}

	if (*s != '\0') {
#if FN_DEBUG
	    fprintf(stderr, "function_args: looking at '%s'\n", s);
#endif
	    argv = get_separated_fields(s, &na, 32, NULL, &err);
	    if (!err) {
		*pargv = argv;
		*argc = na;
	    }
	}
    }

    return err;
}

static int maybe_delete_function (const char *fname)
{
    ufunc *fun = get_ufunc_by_name(fname);
    int err = 0;

    if (fun == NULL) {
	; /* no-op */
    } else if (function_is_on_stack(fun)) {
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

static int 
parse_function_element (char *s, fn_param *param, fn_return *ret, int i)
{
    char tstr[8] = {0};
    char *name;
    int type, len;
    int err = 0;

    while (isspace(*s)) s++;

    /* get arg or return type */
    len = strcspn(s, " ");
    if (len > 7) {
	err = E_PARSE;
    } else {
	strncat(tstr, s, len);
	type = arg_type_from_string(tstr);
	if (type == 0) {
	    sprintf(gretl_errmsg, "Unrecognized data type '%s'", tstr);
	    err = E_PARSE;
	} else if (type == ARG_LIST && ret != NULL) {
	    strcpy(gretl_errmsg, "A function cannot return a list");
	    err = E_DATA;
	}
    }

    if (err) {
	return err;
    }

    s += len;
    while (isspace(*s)) s++;
    len = strcspn(s, " ");
    if (len == 0) {
	if (param != NULL) {
	    sprintf(gretl_errmsg, "parameter %d: name is missing", i);
	} else {
	    sprintf(gretl_errmsg, "return value %d: name is missing", i);
	}
	err = E_PARSE;
    }

    if (err) {
	return err;
    }    

    if (param != NULL && scalar_arg(type)) {
	param->type = type;
	err = read_deflt_min_max(s, param, &len);
    }

    if (!err) {
	name = gretl_strndup(s, len);
	if (name == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	if (param != NULL) {
	    param->type = type;
	    param->name = name;
	} else {
	    ret->type = type;
	    ret->name = name;
	}
    }

#if FN_DEBUG
    if (!err) {
	fprintf(stderr, " %s[%d] = '%s', ptype = %d\n", 
		(param != NULL)? "parm" : "ret",
		i, name, type);
	if (param != NULL) {
	    fprintf(stderr, "min=%g, max=%g, deflt=%g\n", 
		    param->min, param->max, param->deflt);
	}
    }
#endif

    return err;
}

static int 
parse_fn_definition_or_returns (char *fname, 
				fn_param **pparams,
				fn_return **preturns,
				int *pnp,
				const char *str, 
				ufunc **pfun, 
				PRN *prn)
{
    fn_param *params = NULL;
    fn_return *returns = NULL;
    char *p, *s = NULL;
    int i, len, np = 0;
    int which, err = 0;

    if (fname != NULL) {
	which = FN_PARAMS;
    } else {
	which = FN_RETURNS;
    }

    while (isspace(*str)) str++;

    if (which == FN_PARAMS) {
	char fmt[8] = "%31s";

	/* get the function name */
	len = strcspn(str, " (");
	if (len == 0) {
	    err = E_PARSE;
	} else if (len < FN_NAMELEN -1) {
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

    if (!err) {
	/* strip trailing ')' and space */
	tailstrip(s);
	len = strlen(s);
	if (s[len-1] == ')') {
	    s[len-1] = 0;
	}
	np = comma_count(s) + 1;
	if (which == FN_PARAMS) {
	    params = allocate_params(np);
	    if (params == NULL) {
		err = E_ALLOC;
	    }
	} else {
	    returns = allocate_returns(np);
	    if (returns == NULL) {
		err = E_ALLOC;
	    }
	}
	if (!err) {
	    charsub(s, ',', 0);
	}
    }

    if (!err) {
	p = s;
	for (i=0; i<np && !err; i++) {
	    if (params != NULL) {
		err = parse_function_element(p, &params[i], NULL, i);
	    } else {
		err = parse_function_element(p, NULL, &returns[i], i);
	    }
	    p += strlen(p) + 1;
	}
    }

    free(s);

    if (err) {
	free_params_array(params, np);
	free_returns_array(returns, np);
    } else {
	if (pparams != NULL) {
	    *pparams = params;
	} else if (preturns != NULL) {
	    *preturns = returns;
	}
	*pnp = np;
    }
    
    return err;
}

int gretl_start_compiling_function (const char *line, PRN *prn)
{
    ufunc *fun = NULL;
    fn_param *params;
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

    err = parse_fn_definition_or_returns(fname, &params, NULL, &n_params,
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

static int create_function_return_list (ufunc *fun, const char *line)
{
    int err = 0;

    if (fun->returns != NULL) {
	sprintf(gretl_errmsg, "Function %s: return value is already defined",
		fun->name);
	return 1;
    }

    err = parse_fn_definition_or_returns(NULL, NULL, &fun->returns, 
					 &fun->n_returns, line,
					 NULL, NULL);

    return err;
}

static int real_add_fn_line (ufunc *fun, const char *s)
{
    char **lines;
    int nl = fun->n_lines;
    int err = 0;

#if FN_DEBUG > 1
    fprintf(stderr, "real_add_fn_line: '%s'\n", s);
    fprintf(stderr, "currently fun->n_lines = %d\n", nl);
#endif

    lines = realloc(fun->lines, (nl + 1) * sizeof *lines);

    if (lines == NULL) {
	err = E_ALLOC;
    } else {
	fun->lines = lines;
	fun->lines[nl] = gretl_strdup(s);
	if (fun->lines[nl] == NULL) {
	    err = E_ALLOC;
	} else {
	    fun->n_lines += 1;
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
	err = create_function_return_list(fun, line + 7);
    } else {    
	err = real_add_fn_line(fun, line);
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

#define LOCAL_COPY_LISTS 1

#ifdef LOCAL_COPY_LISTS

/* expand by one element the mapping (in the form of two parallel
   lists) between the ID numbers of "outer scope" variables and
   the ID numbers of the corresponding local copies of those
   variables */

static int 
grow_vars_mapping (int outv, int inv, int **p_outlist, int **p_inlist)
{
    int *outlist = NULL, *inlist = NULL;
    int n;

    if (*p_outlist != NULL) {
	n = (*p_outlist)[0] + 2;
    } else {
	n = 2;
    }

    outlist = realloc(*p_outlist, n * sizeof *outlist);
    if (outlist == NULL) {
	return 1;
    }

    inlist = realloc(*p_inlist, n * sizeof *inlist);
    if (inlist == NULL) {
	*p_outlist = outlist;
	return 1;
    }

    outlist[0] = inlist[0] = n - 1;
    
    outlist[n-1] = outv;
    inlist[n-1] = inv;

    *p_outlist = outlist;
    *p_inlist = inlist;

    return 0;
}

/* look up a variable by its "outer" ID number: if we've already made
   a local copy of the variable, return the ID number of the copy
   (otherwise return 0) */

static int find_mapped_var (int v, const int *outlist, const int *inlist)
{
    int i, match = 0;

    if (outlist != NULL && inlist != NULL) {
	for (i=1; i<=outlist[0]; i++) {
	    if (v == outlist[i]) {
		match = inlist[i];
		break;
	    }
	}
    }

    return match;
}

/* given a named list of variables supplied as an argument to a function,
   make a revised version of the list in which the original variable
   ID numbers are replaced by the ID numbers of local copies of those
   variables */

static int localize_list (const char *oldname, const char *newname,
			  double ***pZ, DATAINFO *pdinfo,
			  int **p_outlist, int **p_inlist)
{
    int *orig = NULL;
    int *new = NULL;
    int origv = pdinfo->v;
    int orig0 = 0;
    int i, v, match, err = 0;

    if (!strcmp(oldname, "null")) {
	new = gretl_null_list();
    } else {
	orig = get_list_by_name(oldname);
	if (orig == NULL) {
	    return 1;
	}
	new = gretl_list_copy(orig);
	orig0 = orig[0];
    }

    if (new == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<=orig0 && !err; i++) {
	v = orig[i];
	if (v == pdinfo->v) {
	    /* unknown variable */
	    err = E_DATA;
	} else if ((match = find_mapped_var(v, *p_outlist, *p_inlist))) {
	    /* a local copy of this variable already exists */
	    new[i] = match;
	} else {
	    /* add a local copy of the variable */
	    err = dataset_copy_variable_as(v, pdinfo->varname[v], 
					   pZ, pdinfo);
	    if (!err) {
		new[i] = pdinfo->v - 1;
		err = grow_vars_mapping(v, pdinfo->v - 1, p_outlist, p_inlist);
	    }
	}
    }
	
    if (!err) {
	/* save the new list at appropriate stacking level */
	err = stack_localized_list_as(new, newname);
    } else {
	dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);
	free(new);
    }

    return err;
}

#endif

/* Check number of arguments to function.  If there are trailing
   parameters that have specified default values, it's OK if the
   argument count is short (the defaults will be used).  An excess
   argument count is always an error.
*/

static int check_argc (int argc, ufunc *fun) 
{
    int err = 1;

    if (argc == fun->n_params) {
	err = 0;
    } else if (argc < fun->n_params) {
	int i, nreq = fun->n_params;

	for (i=fun->n_params-1; i>=0; i--) {
	    if (!na(fun->params[i].deflt)) {
		nreq--;
	    } else {
		break;
	    }
	}

	if (argc >= nreq) {
	    err = 0;
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
    char defstr[32];

    if (na(param->deflt)) {
	/* should be impossible here, but... */
	return E_DATA;
    }

    if (param->type == ARG_BOOL || param->type == ARG_INT) {
	sprintf(defstr, "%g", floor(param->deflt));
    } else {
	sprintf(defstr, "%g", param->deflt);
    }
    
    return dataset_add_scalar_as(defstr, param->name, pZ, pdinfo);
}

static int check_and_allocate_function_args (ufunc *fun,
					     int argc, char **argv, 
					     double ***pZ,
					     DATAINFO *pdinfo)
{
    int *outlist = NULL;
    int *inlist = NULL;
    int i, v, err;

    err = check_argc(argc, fun);
    if (err) {
	sprintf(gretl_errmsg, _("Number of arguments (%d) does not "
				"match the number of\nparameters for "
				"function %s (%d)"),
		argc, fun->name, fun->n_params);
	return err;
    }

    for (i=0; i<fun->n_params && !err; i++) {
#if FN_DEBUG
	fprintf(stderr, "fn argv[%d]: arg='%s', param.name='%s' param.type=%d\n", 
		i, argv[i], fun->params[i].name, fun->params[i].type);
#endif
	if (scalar_arg(fun->params[i].type)) {
	    if (i >= argc) {
		err = add_scalar_arg_default(&fun->params[i], pZ, pdinfo);
	    } else if (numeric_string(argv[i])) {
		err = dataset_add_scalar_as(argv[i], fun->params[i].name, 
					    pZ, pdinfo);
	    } else {
		v = varindex(pdinfo, argv[i]);
		if (var_is_scalar(pdinfo, v)) {
		    err = dataset_copy_variable_as(v, fun->params[i].name, 
						   pZ, pdinfo);
		    if (!err) {
			err = grow_vars_mapping(v, pdinfo->v - 1, &outlist, &inlist);
		    }
		} else {
		    sprintf(gretl_errmsg, "argument %d (%s): not a scalar", 
			    i+1, argv[i]);
		    err = 1;
		} 
	    }
	} else if (fun->params[i].type == ARG_SERIES) {
	    if (numeric_string(argv[i])) {
		v = atoi(argv[i]);
	    } else {
		v = varindex(pdinfo, argv[i]);
	    }
	    if (v < pdinfo->v && var_is_series(pdinfo, v)) {
		err = dataset_copy_variable_as(v, fun->params[i].name, 
					       pZ, pdinfo);
		if (!err) {
		    err = grow_vars_mapping(v, pdinfo->v - 1, &outlist, &inlist);
		}
	    } else {
		sprintf(gretl_errmsg, "argument %d (%s): not a series", 
			i+1, argv[i]);
		err = 1;
	    } 
	} else if (fun->params[i].type == ARG_LIST) {
	    if (get_list_by_name(argv[i]) != NULL || !strcmp(argv[i], "null")) {
#ifdef LOCAL_COPY_LISTS
		err = localize_list(argv[i], fun->params[i].name, pZ, pdinfo,
				    &outlist, &inlist);
#else
		err = copy_named_list_as(argv[i], fun->params[i].name);
#endif
	    } else {
		sprintf(gretl_errmsg, "argument %d (%s): not a list", i+1, argv[i]);
		err = 1;
	    }
	} else if (fun->params[i].type == ARG_MATRIX) {
	    if (get_matrix_by_name(argv[i], pdinfo) != NULL) {
		err = copy_named_matrix_as(argv[i], fun->params[i].name);
#if FN_DEBUG
		fprintf(stderr, "done copy_named_matrix_as, '%s' -> '%s', err = %d\n", 
			argv[i], fun->params[i].name, err);
#endif
	    }
	} else {
	    /* impossible */
	    err = 1;
	}
    }

    if (inlist != NULL) {
	free(inlist);
    }
    if (outlist != NULL) {
	free(outlist);
    }

    return err;
}

static int check_function_assignments (ufunc *fun, int *asslist, 
				       char **assv, char *asstypes,
				       DATAINFO *pdinfo)
{
    int i, v, err = 0;

    if (asslist[0] > fun->n_returns) {
	sprintf(gretl_errmsg, _("Number of assignments (%d) exceeds the "
				"number of values returned by\n%s (%d)"), 
		asslist[0], fun->name, fun->n_returns);
	err = 1;
    }

    for (i=0; i<asslist[0] && !err; i++) {
#if FN_DEBUG
	fprintf(stderr, "checking assignment %d to '%s'\n", i, assv[i]);
#endif
	if (assv[i] == NULL || *assv[i] == '\0') {
	    /* placeholder, non-assignment */
	    continue;
	}
	v = varindex(pdinfo, assv[i]);
	if (v < pdinfo->v) {
#if FN_DEBUG
	    fprintf(stderr, " variable '%s' has ID %d, series = %d\n", assv[i], 
		    v, var_is_series(pdinfo, v));
#endif
	    if ((fun->returns[i].type == ARG_SCALAR && var_is_series(pdinfo, v)) ||
		(fun->returns[i].type == ARG_SERIES && var_is_scalar(pdinfo, v))) {
		sprintf(gretl_errmsg, "%s: wrong type for assignment", assv[i]);
		err = 1;
	    } else if (fun->returns[i].type == ARG_MATRIX) {
		sprintf(gretl_errmsg, "%s: wrong type for assignment", assv[i]);
		err = 1;
	    } else {
		asslist[i+1] = v;
	    }
	} else if (fun->returns[i].type == ARG_LIST) {
	    fprintf(stderr, "requested return of list as '%s'\n", assv[i]);
	    /* FIXME? */
	} else {
	    /* new variable */
	    if (asstypes != NULL && asstypes[i] && asstypes[i] != fun->returns[i].type) {
		sprintf(gretl_errmsg, "%s: wrong type for assignment", assv[i]);
		err = 1;
	    }
	    if (!err) {
		err = check_varname(assv[i]);
	    }
	}
    }

    return err;
}

static int check_function_data_needs (const DATAINFO *pdinfo,
				      const ufunc *fun)
{
    const fnpkg *pkg;

    pkg = ufunc_get_parent_package(fun);
    if (pkg == NULL) {
	return 0;
    }

    if ((pkg->dreq == FN_NEEDS_TS) && 
	!dataset_is_time_series(pdinfo)) {
	return 1;
    }

    if ((pkg->dreq == FN_NEEDS_PANEL) && 
	!dataset_is_panel(pdinfo)) {
	return 1;
    }

    if ((pkg->dreq == FN_NEEDS_QM) && 
	(!dataset_is_time_series(pdinfo) || 
	 (pdinfo->pd != 4 || pdinfo->pd != 12))) {
	return 1;
    } 

    return 0;
}

int gretl_function_start_exec (const char *line, const char *fname, 
			       double ***pZ, DATAINFO *pdinfo)
{
    char **argv = NULL;
    char **assv = NULL;
    char *asstypes = NULL;
    int *asslist = NULL;
    ufunc *fun;
    fncall *call = NULL;
    int argc = 0;
    int err = 0;

#if FN_DEBUG
    fprintf(stderr, "gretl_function_start_exec: line='%s'\n"
	    " fname='%s'\n", line, fname);
#endif

    fun = get_ufunc_by_name(fname);

    if (fun == NULL) {
	return 1;
    }

    err = check_function_data_needs(pdinfo, fun);
    if (err) {
	return err;
    }

    err = parse_function_args_etc(line, &argc, &argv, &asslist, &assv, &asstypes);
    if (err) {
	return err;
    }

    err = check_and_allocate_function_args(fun, argc, argv, pZ, pdinfo);
    
    if (!err && asslist[0] > 0) {
	err = check_function_assignments(fun, asslist, assv, asstypes, pdinfo);
    }

    free(asstypes);

    if (err) {
	free_strings_array(argv, argc);
	free_strings_array(assv, asslist[0]);
	free(asslist);
	return 1;
    }

    call = fncall_new(fun, argc, argv, asslist, assv);
    if (call == NULL) {
	return E_ALLOC;
    } 

    err = push_fncall(call, pdinfo);
    if (err) {
	free_fncall(call);
    }

#if FN_DEBUG
    fprintf(stderr, "gretl_function_start_exec: returning err = %d\n", err);
#endif

    return err;
}

char *gretl_function_get_line (char *line, int len,
			       double ***pZ, DATAINFO **ppdinfo,
			       int *err)
{
    fncall *call = current_call();
    const char *src;
    int unstack = 0;

    *err = 0;

    if (call == NULL || call->fun == NULL) {
#if FN_DEBUG
	fprintf(stderr, "gretl_function_get_line: returning NULL\n"); 
#endif
	return NULL;
    }

#if FN_DEBUG
    fprintf(stderr, "gretl_function_get_line: current fun='%s'\n", 
	    call->fun->name);
#endif

    if (call->lnum > call->fun->n_lines - 1) {
	/* finished executing */
	unstack = 1;
    } else {
	src = call->fun->lines[call->lnum];
	if (!strncmp(src, "exit", 4)) {
	    /* terminate execution */
	    unstack = 1;
	}
    } 

    if (unstack) {
	*err = unstack_fncall(pZ, ppdinfo);
	return "";
    }	

    call->lnum += 1;
    strcpy(line, src);

    return line;
}

void gretl_functions_cleanup (void)
{
    callstack_destroy();
    ufuncs_destroy();
    packages_destroy();
}

void gretl_function_stop_on_error (double ***pZ, DATAINFO **ppdinfo, PRN *prn)
{
    gretl_function_flagged_error(NULL, prn);

    callstack_destroy();

    if (fn_executing > 0) {
	libset_restore_state_zero(pZ, ppdinfo);
    }

    fn_executing = 0;
}

int gretl_function_flagged_error (const char *s, PRN *prn)
{
    fncall *call;

    if (!fn_executing) {
	return 0;
    }

    call = current_call();

    if (s != NULL && *s != '\0') {
	pprintf(prn, "%s: %s\n", call->fun->name, s);
    } else {
	pprintf(prn, _("Error condition in execution of function %s"),
		(callstack[0])->fun->name);
	pputc(prn, '\n');
    }

    return 1;
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

    if (fun->n_returns > 0) {
	pprintf(prn, "Return values:\n");
	for (i=0; i<fun->n_returns; i++) {
	    pprintf(prn, " %s (%s)\n", 
		    fun->returns[i].name, arg_type_string(fun->returns[i].type));
	}
	pputc(prn, '\n');
    } else {
	pputs(prn, "Return values: none\n\n");
    }
	
    if (fun->help != NULL) {
	pprintf(prn, "Help text:\n%s\n\n", fun->help);
    }
}

int user_function_help (const char *fnname, PRN *prn)
{
    ufunc *fun = get_ufunc_by_name(fnname);
    int err = 0;

    if (fun == NULL) {
	pprintf(prn, _("\"%s\" is not defined.\n"), fnname);
	err = 1;
    } else {
	real_user_function_help(fun, NULL, prn);
    }

    return err;
}

