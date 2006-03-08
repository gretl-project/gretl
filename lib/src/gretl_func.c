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

#define CALLSTACK_DEPTH 8

#define FN_DEBUG 0

typedef struct ufunc_ ufunc;
typedef struct fncall_ fncall;

#define FN_NAMELEN 32

enum {
    ARG_NONE = 0,
    ARG_SCALAR,
    ARG_SERIES,
    ARG_LIST,
    ARG_MATRIX
};

struct ufunc_ {
    char name[FN_NAMELEN];
    int n_lines;
    char **lines;
    int n_returns;
    char **returns;
    char *rtype;
    int n_params;
    char **params;
    char *ptype;
};

struct fncall_ {
    ufunc *fun;
    int lnum;
    int argc;
    char **argv;
    char **assv;
    int *asslist;
};

static int n_ufuns;
static ufunc **ufuns;

static fncall **callstack;

static void free_fncall (fncall *call);

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
    int t, n = (pdinfo->vector[targ])? pdinfo->n : 1;

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
	    call->asslist[0], call->fun->n_returns);
#endif

    if (call->fun->n_returns < call->asslist[0]) {
	fprintf(stderr, "bad function call: assignments > returns\n");
	return 1;
    }

    for (i=0; i<call->asslist[0]; i++) {
	targ = call->asslist[i+1];
	src = -1;

#if FN_DEBUG
	fprintf(stderr, " assv[%d] = '%s'\n", i, call->assv[i]);
	fprintf(stderr, " return[%d] = '%s'\n", i, call->fun->returns[i]);
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
	    if (!strcmp(pdinfo->varname[j], call->fun->returns[i])) {
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
	    gretl_matrix *S = get_matrix_by_name(call->fun->returns[i],
						 pdinfo);


	    if (S != NULL) {
		gretl_matrix *T;
#if FN_DEBUG
		fprintf(stderr, " identified source as matrix '%s' (%p)\n",
			call->fun->returns[i], (void *) S);
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
		    fprintf(stderr, " renaming matrix '%s' (%p) as '%s'\n",
			    call->fun->returns[i], (void *) S, call->assv[i]);
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
	if (call->fun->rtype[i] == ARG_MATRIX) {
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

/* constructors */

static ufunc *ufunc_new (void)
{
    ufunc *func = malloc(sizeof *func);

    if (func == NULL) {
	return NULL;
    }

    func->name[0] = '\0';

    func->n_lines = 0;
    func->lines = NULL;

    func->n_returns = 0;
    func->returns = NULL;
    func->rtype = NULL;

    func->n_params = 0;
    func->params = NULL;
    func->ptype = NULL;

    return func;
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

static ufunc *add_ufunc (void)
{
    int nf = n_ufuns;
    ufunc **myfuns;

    myfuns = realloc(ufuns, (nf + 1) * sizeof *myfuns);
    if (myfuns == NULL) {
	return NULL;
    }
    ufuns = myfuns;

    ufuns[nf] = ufunc_new();
    if (ufuns[nf] == NULL) {
	return NULL;
    }

    n_ufuns++;

    return ufuns[nf];
}

/* destructors */

static void free_ufunc (ufunc *fun)
{
    free_strings_array(fun->lines, fun->n_lines);
    free_strings_array(fun->returns, fun->n_returns);
    free_strings_array(fun->params, fun->n_params);
    free(fun->ptype);
    free(fun->rtype);

    free(fun);
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

static int check_func_name (const char *fname, PRN *prn)
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
		delete_ufunc_from_list(ufuns[i]);
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

static int is_type_name (const char *s)
{
    int ret = ARG_NONE;

    if (!strcmp(s, "scalar")) {
	ret = ARG_SCALAR;
    } else if (!strcmp(s, "series")) {
	ret = ARG_SERIES;
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
	    fields = create_strings_array(nf);
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
		} else if ((atype = is_type_name(fields[j]))) {
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

static int get_param_type (const char *s)
{
    int ret = 0;

    if (!strcmp(s, "scalar")) {
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

static int comma_count (const char *s)
{
    int nc = 0;

    while (*s) {
	if (*s == ',') nc++;
	s++;
    }

    return nc;
}

static int allocate_parmv_ptype (char ***pparmv, char **pptype, int n)
{
    char **parmv;
    char *ptype;
    int err = 0;

    parmv = create_strings_array(n);
    if (parmv == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	ptype = malloc(n * sizeof *ptype);
	if (ptype == NULL) {
	    free(parmv);
	    err = E_ALLOC;
	}
    }

    if (!err) {
	*pparmv = parmv;
	*pptype = ptype;
    }

    return err;
}

enum {
    FN_PARAMS,
    FN_RETURNS
};

static int parse_fn_element (char *s, char **parmv, char *ptype, int i,
			     int which)
{
    char tstr[8] = {0};
    int n, len;
    int err = 0;

    while (isspace((unsigned char) *s)) s++;
    len = strcspn(s, " ");
    if (len < 7) {
	n = len;
    } else {
	n = 7;
    }

    strncat(tstr, s, n);
    ptype[i] = get_param_type(tstr);
    if (ptype[i] == 0) {
	sprintf(gretl_errmsg, "Unrecognized data type '%s'", tstr);
	err = E_PARSE;
    } else if (ptype[i] == ARG_LIST && which == FN_RETURNS) {
	/* FIXME? */
	strcpy(gretl_errmsg, "A function cannot return a list");
	err = 1;
    }

    if (!err) {
	s += len;
	while (isspace((unsigned char) *s)) s++;
	len = strcspn(s, " ");
	if (len == 0) {
	    if (which == FN_PARAMS) {
		sprintf(gretl_errmsg, "parameter %d: name is missing", i);
	    } else {
		sprintf(gretl_errmsg, "return value %d: name is missing", i);
	    }
	    err = E_PARSE;
	}
    }

    if (!err) {
	if (len < 31) {
	    n = len;
	} else {
	    n = 31;
	}
	parmv[i] = gretl_strndup(s, n);
	if (parmv[i] == NULL) {
	    err = E_ALLOC;
	}
    }

#if FN_DEBUG
    if (!err) {
	fprintf(stderr, " %s[%d] = '%s', ptype = %d\n", 
		(which == FN_PARAMS)? "parm" : "ret",
		i, parmv[i], ptype[i]);
    }
#endif

    return err;
}

static int 
parse_fn_definition_or_returns (char *fname, char ***pparmv, int *pnp,
				char **pptype, const char *str,
				PRN *prn)
{
    char **parmv = NULL;
    char *ptype = NULL;
    char *p, *s = NULL;
    int i, len, np = 0;
    int which, err = 0;

    if (fname != NULL) {
	which = FN_PARAMS;
    } else {
	which = FN_RETURNS;
    }

    while (isspace((unsigned char) *str)) str++;

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
	    err = check_func_name(fname, prn);
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
	err = allocate_parmv_ptype(&parmv, &ptype, np);
	charsub(s, ',', 0);
    }

    if (!err) {
	p = s;
	for (i=0; i<np && !err; i++) {
	    err = parse_fn_element(p, parmv, ptype, i, which);
	    p += strlen(p) + 1;
	}
    }

    free(s);

    if (err) {
	free_strings_array(parmv, np);
	free(ptype);
    } else {
	*pparmv = parmv;
	*pptype = ptype;
	*pnp = np;
    }
    
    return err;
}

int gretl_start_compiling_function (const char *line, PRN *prn)
{
    char **params = NULL;
    char *ptype = NULL;
    int n_params = 0;
    char fname[FN_NAMELEN];
    char extra[8];
    int nf;
    ufunc *fun = NULL;
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

    err = parse_fn_definition_or_returns(fname, &params, &n_params, 
					 &ptype, line + 8, prn);

    if (!err) {
	fun = add_ufunc();
	if (fun == NULL) {
	    free_strings_array(params, n_params);
	    free(ptype);
	    err = E_ALLOC;
	}
    }

    if (!err) {
	strcpy(fun->name, fname);
	fun->params = params;
	fun->n_params = n_params;
	fun->ptype = ptype;
	set_compiling_on();
    }
    
    return err;
}

static ufunc *get_latest_ufunc (void)
{
    if (n_ufuns > 0) {
	return ufuns[n_ufuns - 1];
    } else {
	return NULL;
    }
}

static int create_function_return_list (ufunc *fun, const char *line)
{
    int err = 0;

    if (fun->returns != NULL) {
	sprintf(gretl_errmsg, "Function %s: return value is already defined",
		fun->name);
	return 1;
    }

    err = parse_fn_definition_or_returns(NULL, &fun->returns, 
					 &fun->n_returns, 
					 &fun->rtype, 
					 line,
					 NULL);

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
    ufunc *fun = get_latest_ufunc();
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

static int check_and_allocate_function_args (ufunc *fun,
					     int argc, char **argv, 
					     double ***pZ,
					     DATAINFO *pdinfo)
{
    int *outlist = NULL, *inlist = NULL;
    int i, v, err = 0;

    if (argc != fun->n_params) {
	sprintf(gretl_errmsg, _("Number of arguments (%d) does not "
				"match the number of\nparameters for "
				"function %s (%d)"),
		argc, fun->name, fun->n_params);
	err = 1;
    }

    for (i=0; i<argc && !err; i++) {
#if FN_DEBUG
	fprintf(stderr, "fn argv[%d]: arg='%s', fun->param='%s' fun->ptype=%d\n", 
		i, argv[i], fun->params[i], fun->ptype[i]);
#endif
	if (fun->ptype[i] == ARG_SCALAR) {
	    if (numeric_string(argv[i])) {
		err = dataset_add_scalar_as(argv[i], fun->params[i], pZ, pdinfo);
	    } else {
		v = varindex(pdinfo, argv[i]);
		if (v < pdinfo->v && !pdinfo->vector[v]) {
		    err = dataset_copy_variable_as(v, fun->params[i], pZ, pdinfo);
		    if (!err) {
			err = grow_vars_mapping(v, pdinfo->v - 1, &outlist, &inlist);
		    }
		} else {
		    sprintf(gretl_errmsg, "argument %d (%s): not a scalar", i+1, argv[i]);
		    err = 1;
		} 
	    }
	} else if (fun->ptype[i] == ARG_SERIES) {
	    if (numeric_string(argv[i])) {
		v = atoi(argv[i]);
	    } else {
		v = varindex(pdinfo, argv[i]);
	    }
	    if (v < pdinfo->v && pdinfo->vector[v]) {
		err = dataset_copy_variable_as(v, fun->params[i], pZ, pdinfo);
		if (!err) {
		    err = grow_vars_mapping(v, pdinfo->v - 1, &outlist, &inlist);
		}
	    } else {
		sprintf(gretl_errmsg, "argument %d (%s): not a series", i+1, argv[i]);
		err = 1;
	    } 
	} else if (fun->ptype[i] == ARG_LIST) {
	    if (get_list_by_name(argv[i]) != NULL || !strcmp(argv[i], "null")) {
#ifdef LOCAL_COPY_LISTS
		err = localize_list(argv[i], fun->params[i], pZ, pdinfo,
				    &outlist, &inlist);
#else
		err = copy_named_list_as(argv[i], fun->params[i]);
#endif
	    } else {
		sprintf(gretl_errmsg, "argument %d (%s): not a list", i+1, argv[i]);
		err = 1;
	    }
	} else if (fun->ptype[i] == ARG_MATRIX) {
	    if (get_matrix_by_name(argv[i], pdinfo) != NULL) {
		err = copy_named_matrix_as(argv[i], fun->params[i]);
#if FN_DEBUG
		fprintf(stderr, "done copy_named_matrix_as, '%s' -> '%s', err = %d\n", 
			argv[i], fun->params[i], err);
#endif
	    }
	} else {
	    /* impossible */
	    err = 1;
	}
    }

    free(outlist);
    free(inlist);

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
	    fprintf(stderr, " variable '%s' has ID %d, vector = %d\n", assv[i], 
		    v, pdinfo->vector[v]);
#endif
	    if ((fun->rtype[i] == ARG_SCALAR && pdinfo->vector[v]) ||
		(fun->rtype[i] == ARG_SERIES && !pdinfo->vector[v])) {
		sprintf(gretl_errmsg, "%s: wrong type for assignment", assv[i]);
		err = 1;
	    } else if (fun->rtype[i] == ARG_MATRIX) {
		sprintf(gretl_errmsg, "%s: wrong type for assignment", assv[i]);
		err = 1;
	    } else {
		asslist[i+1] = v;
	    }
	} else if (fun->rtype[i] == ARG_LIST) {
	    fprintf(stderr, "requested return of list as '%s'\n", assv[i]);
	    /* FIXME? */
	} else {
	    /* new variable */
	    if (asstypes != NULL && asstypes[i] && asstypes[i] != fun->rtype[i]) {
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

int gretl_function_start_exec (const char *line, double ***pZ,
			       DATAINFO *pdinfo)
{
    char **argv = NULL;
    char **assv = NULL;
    char *asstypes = NULL;
    int *asslist = NULL;
    char fname[FN_NAMELEN];
    ufunc *fun;
    fncall *call = NULL;
    int argc = 0;
    int err = 0;

    function_name_from_line(line, fname);
    fun = get_ufunc_by_name(fname);

    if (fun == NULL) {
	return 1;
    }

    err = parse_function_args_etc(line, &argc, &argv, &asslist, &assv, &asstypes);
    if (err) {
	return err;
    }

    if (argc > 0) {
	err = check_and_allocate_function_args(fun, argc, argv, pZ, pdinfo);
    }
    
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
