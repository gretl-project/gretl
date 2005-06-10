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

#define CALLSTACK_DEPTH 8

#define FN_DEBUG 0

typedef struct ufunc_ ufunc;
typedef struct fncall_ fncall;

struct ufunc_ {
    char name[32];
    int is_macro;
    int n_lines;
    char **lines;
    int n_returns;
    char **returns;
    int n_params;
    char **params;
};

struct fncall_ {
    ufunc *fun;
    int lnum;
    int argc;
    char **argv;
    int assc;
    char **assv;
};

static int n_ufuns;
static ufunc **ufuns;

static fncall **callstack;

static void free_fncall (fncall *call);

/* record of state, and communication of state with outside world */

static int compiling;
static int fn_executing;
static int macro_executing;

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

int gretl_executing_macro (void)
{
    return macro_executing;
}

int gretl_executing_function_or_macro (void)
{
    return (fn_executing || macro_executing);
}

static void set_executing_on (fncall *call)
{
    if (call->fun->is_macro) {
	macro_executing++;
    } else {
	fn_executing++;
    }
#if FN_DEBUG
    fprintf(stderr, "set_executing_on: fun=%s, macro_executing=%d, "
	    "fn_executing=%d\n", call->fun->name, macro_executing,
	    fn_executing);
#endif
}

static void set_executing_off (fncall *call)
{
    if (call->fun->is_macro) {
	macro_executing--;
    } else {
	fn_executing--;
    }
#if FN_DEBUG
    fprintf(stderr, "set_executing_off: fun=%s, macro_executing=%d, "
	    "fn_executing=%d\n", call->fun->name, macro_executing,
	    fn_executing);
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

static int real_gretl_function_stack_depth (void)
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

int gretl_function_stack_depth (void)
{
    return real_gretl_function_stack_depth();
}

int gretl_macro_stack_depth (void)
{
    return real_gretl_function_stack_depth();
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

    if (!call->fun->is_macro) {
	push_program_state(pdinfo);
    }

    return 0;
}

static int 
destroy_local_vars (fncall *call, double ***pZ, DATAINFO *pdinfo, int nc)
{
    int n_returns = call->fun->n_returns;
    char **returns = call->fun->returns;
    int saves = call->assc;
    int *locals = NULL;
    int i, nlocal = 0;
    int err = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (STACK_LEVEL(pdinfo, i) == nc) {
	    nlocal++;
	}
    }

    if (nlocal > 0) {
	locals = gretl_list_new(nlocal);
	if (locals == NULL) {
	    err = 1;
	} else {
	    locals[0] = 0;
	}
    }

    if (locals != NULL) {
	int k, v = 0, j = 1;
	int wanted;

	for (i=1; i<pdinfo->v; i++) {
	    if (STACK_LEVEL(pdinfo, i) != nc) {
		continue;
	    }

	    wanted = 0;

	    if (saves > 0) {
		for (k=0; k<n_returns; k++) {
		    if (!strcmp(pdinfo->varname[i], returns[k])) {
#if FN_DEBUG
			fprintf(stderr, "%s: this var is wanted by caller\n",
				pdinfo->varname[i]);
#endif
			wanted = 1;
			break;
		    }
		}
	    }

	    if (wanted) {
		/* rename variable as caller desired */
		/* FIXME: what if var already exists at caller level? */
		strcpy(pdinfo->varname[i], call->assv[v++]);
		STACK_LEVEL(pdinfo, i) -= 1; 
		saves--;
	    } else {
#if FN_DEBUG
		fprintf(stderr, "local variable %d (%s) "
			"marked for deletion\n", i,
			pdinfo->varname[i]);
#endif
		locals[j++] = i;
		locals[0] += 1;
	    }
	}

	err = dataset_drop_listed_variables(locals, pZ, pdinfo, NULL);
	free(locals);
    }

    if (!err) {
	err = destroy_saved_lists_at_level(nc);
    }

    return err;
}

static int unstack_fncall (double ***pZ, DATAINFO *pdinfo)
{
    fncall *call;   
    int i, nc;
    int err = 0;

    if (callstack == NULL) {
	return 1;
    }
    
    call = callstack[0];

    /* FIXME */
    if (call->fun->is_macro) {
	nc = gretl_macro_stack_depth();
    } else {
	nc = gretl_function_stack_depth();
    }

#if FN_DEBUG
    fprintf(stderr, "unstack_fncall: terminating call to "
	    "function '%s' at depth %d\n", 
	    call->fun->name, nc);
#endif

    err = destroy_local_vars(call, pZ, pdinfo, nc);

    set_executing_off(call);

    if (!call->fun->is_macro) {
	pop_program_state(pdinfo);
    }

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

    for (i=0; i<CALLSTACK_DEPTH; i++) {
	if (callstack[i] == NULL) {
	    break;
	}
	if ((callstack[i])->fun == func) {
	    return 1;
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
    func->is_macro = 0;

    func->n_lines = 0;
    func->lines = NULL;

    func->n_returns = 0;
    func->returns = NULL;

    func->n_params = 0;
    func->params = NULL;

    return func;
}

static fncall *fncall_new (ufunc *fun, int argc, char **argv,
			   int assc, char **assv)
{
    fncall *call = malloc(sizeof *call);

    if (call == NULL) {
	if (argc > 0) {
	    int i;

	    for (i=0; i<argc; i++) {
		free(argv[i]);
	    }
	    free(argv);
	}
	if (assc > 0) {
	    int i;

	    for (i=0; i<assc; i++) {
		free(assv[i]);
	    }
	    free(assv);
	}
	return NULL;
    }

    call->fun = fun;
    call->lnum = 0;

    call->argc = argc;
    call->argv = argv;

    call->assc = assc;
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
    free_strings_array(call->assv, call->assc);

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

	if (p == NULL) {
	    p = line;
	} else {
	    p++;
	}

	sscanf(p, "%31s", name);
#if FN_DEBUG
	fprintf(stderr, "function_name_from_line: got '%s'\n", name);
#endif
    }

    return name;
}

int gretl_is_user_function (const char *line)
{
    int ret = 0;

    if (n_ufuns > 0 && !string_is_blank(line)) {
	char name[32];

#if FN_DEBUG
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

enum function_name_returns {
    FN_NAME_OK,
    FN_NAME_BAD,
    FN_NAME_TAKEN
};

static int check_func_name (const char *fname)
{
    int i;

    if (!isalpha((unsigned char) *fname)) {
	strcpy(gretl_errmsg, "function names must start with a letter");
	return FN_NAME_BAD;
    }

    if (gretl_command_number(fname)) {
	sprintf(gretl_errmsg, "'%s' is the name of a gretl command",
		fname);
	return FN_NAME_BAD;
    }

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fname, (ufuns[i])->name)) {
	    sprintf(gretl_errmsg, "'%s': function is already defined",
		    fname);
	    return FN_NAME_TAKEN;
	}
    }

    return FN_NAME_OK;
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

static char **get_comma_separated_strings (const char *s, int *argc)
{
    char **argv = NULL;
    int i, na, err = 0;

    *argc = 0;

    /* count comma-separated arguments */
    na = comma_count(s) + 1;

    argv = malloc(na * sizeof *argv);
    if (argv == NULL) {
	return NULL;
    }

    for (i=0; i<na; i++) {
	char *arg;
	int len;

	if (i < na - 1) {
	    len = strcspn(s, ",");
	} else {
	    len = strlen(s);
	}

	arg = gretl_strndup(s, len);
	if (arg == NULL) {
	    na = i;
	    err = 1;
	    break;
	}
	argv[i] = arg;

	s += len;
	while (*s) {
	    if (*s != ',' && *s != ' ') break;
	    s++;
	}
    }

    if (err) {
	for (i=0; i<na; i++) {
	    free(argv[i]);
	}
	free(argv);
	argv = NULL;
    } else {
	*argc = na;
    }

    return argv;
}

/* Parse line and return an allocated array of strings consisting of
   the space- or comma-separated fields in line, each one truncated if
   necessary to a maximum of 8 characters.
*/ 

static char **get_separated_fields_8 (const char *line, int *nfields)
{
    char **fields = NULL;
    char *cpy, *s;
    char str[9];
    int i, j, nf, err = 0;

    *nfields = 0;

    if (string_is_blank(line)) {
	return NULL;
    }

    cpy = gretl_strdup(line);
    if (cpy == NULL) {
	return NULL;
    }

    s = cpy;
    charsub(s, ',', ' ');

    nf = count_fields(s);
    fields = malloc(nf * sizeof *fields);
    if (fields == NULL) {
	free(s);
	return NULL;
    }

    for (i=0; i<nf && !err; i++) {
	sscanf(s, "%8s", str);
	fields[i] = gretl_strdup(str);
	if (fields[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(fields[j]);
	    }
	    free(fields);
	    fields = NULL;
	    err = 1;
	} else {
	    s += strcspn(s, " ");
	    s += strspn(s, " ");
	}
    }

    if (!err) {
	*nfields = nf;
    }

    free(cpy);

    return fields;
}

static char **parse_assignment (char *s, int *na)
{
    char **vnames;
    int n, err = 0;

    /* clean the string up */
    s += strspn(s, " ");
    if (*s == '(') {
	s++;
    }
    tailstrip(s);
    n = strlen(s);
    if (s[n-1] == ')') {
	s[n-1] = '\0';
    }

    vnames = get_separated_fields_8(s, na);

    if (vnames != NULL) {
	int i;

	for (i=0; i<*na && !err; i++) {
#if FN_DEBUG
	    fprintf(stderr, "assignee %d: '%s'\n", i, vnames[i]);
#endif
	    err = check_varname(vnames[i]);
	}
    }

    if (err) {
	free_strings_array(vnames, *na);
	vnames = NULL;
    }

    return vnames;
}

static int 
parse_function_args_etc (const char *s, int *argc, char ***pargv,
			 int *assc, char ***passv, int is_macro)
{
    char **argv = NULL;
    char **assign = NULL;
    const char *p;
    int na, err = 0;

    *argc = 0;
    *pargv = NULL;
    *assc = 0;
    *passv = NULL;

    if ((p = strchr(s, '=')) != NULL && *s != '=') {
	char *astr = gretl_strndup(s, p - s);

	/* seems we have a left-hand side assignment */
	assign = parse_assignment(astr, &na);
	if (assign == NULL) {
	    err = 1;
	} else {
	    *passv = assign;
	    *assc = na;
	}
	free(astr);
	s = p + 1;
    } else {
	s++;
    }

    if (!err) {
	/* skip over function name and spaces before args */
	s += strspn(s, " ");
	s += strcspn(s, " ");
	s += strspn(s, " ");

	if (*s != '\0') {
#if FN_DEBUG
	    fprintf(stderr, "function_args: looking at '%s'\n", s);
#endif
	    if (is_macro) {
		argv = get_comma_separated_strings(s, &na);
	    } else {
		argv = get_separated_fields_8(s, &na);
	    }
	    if (argv == NULL) {
		err = 1;
	    } else {
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

static int parse_func_name_and_params (char *name, char ***params,
				       int *n_params, const char *line)
{
    char **pm1 = NULL;
    char **pm2 = NULL;
    char *fname, *param;
    const char *p, *q;
    int np = 0, err = 0;

    p = strchr(line, '(');
    q = strchr(line, ')');

    if ((p != NULL && q == NULL) ||
	(p == NULL && q != NULL) ||
	q - p < 0) {
	err = E_PARSE;
    }

    if (!err) {
	fname = gretl_word_strdup(line, &p);
	if (fname == NULL) {
	    err = E_ALLOC;
	} else {
	    *name = '\0';
	    strncat(name, fname, 31);
	    free(fname);
	}
    }

    while (*p && !err) {
	param = gretl_word_strdup(p, &p);
	if (param != NULL) {
	    pm2 = realloc(pm1, (np + 1) * sizeof *pm1);
	    if (pm2 == NULL) {
		err = E_ALLOC;
	    } else {
		pm2[np] = param;
		pm1 = pm2;
		np++;
	    }	    
	}
    }

    *params = pm1;
    *n_params = np;

    return err;
}

int gretl_start_compiling_function (const char *line)
{
    char **params = NULL;
    int n_params = 0;
    char fname[32];
    char extra[8];
    int nm, nf;
    int name_status;
    int is_macro = 0;
    ufunc *fun = NULL;

    /* for now (May 2005), "function" means macro and
       "newfunc" means function, but this will change */

    nm = sscanf(line, "function %31s %7s", fname, extra);
    nf = sscanf(line, "newfunc %31s %7s", fname, extra);

    if (nm <= 0 && nf <= 0) {
	return E_PARSE;
    } 

    if (nm == 2 || nf == 2) {
	if (!strcmp(extra, "clear") || !strcmp(extra, "delete")) {
	    maybe_delete_function(fname);
	    return 0;
	}
    }

    if (nm > 0) {
	is_macro = 1;
    }

    if (!is_macro) {
	parse_func_name_and_params(fname, &params, &n_params, line + 8);
    }

    name_status = check_func_name(fname);
    if (name_status == FN_NAME_BAD || name_status == FN_NAME_TAKEN) {
	return 1;
    }

    fun = add_ufunc();
    if (fun == NULL) {
	free_strings_array(params, n_params);
	return E_ALLOC;
    }

    strcpy(fun->name, fname);
    fun->is_macro = is_macro;

    fun->params = params;
    fun->n_params = n_params;

    set_compiling_on();
    
    return 0;
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
    int nr, err = 0;
#if FN_DEBUG
    int i;
#endif

    line += strspn(line, " ");

    fun->returns = get_separated_fields_8(line, &nr);

    if (fun->returns != NULL) {
	fun->n_returns = nr;
#if FN_DEBUG
	fprintf(stderr, "done create_function_return_list:\n");
	for (i=0; i<nr; i++) {
	    fprintf(stderr, " return %d = '%s'\n", i, fun->returns[i]);
	}
#endif
    } else {
	err = E_ALLOC;
    }

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

/* FIXME headers */
extern int dataset_copy_variable_as (int v, const char *newname,
				     double ***pZ, DATAINFO *pdinfo);
int dataset_add_scalar_as (const char *numstr, const char *newname,
			   double ***pZ, DATAINFO *pdinfo);

static int check_and_allocate_function_args (ufunc *fun,
					     int argc, char **argv, 
					     double ***pZ,
					     DATAINFO *pdinfo)
{
    int i, v, err = 0;

    if (argc != fun->n_params) {
	sprintf(gretl_errmsg, _("Number of arguments (%d) does not "
				"match the number of\nparameters for "
				"function %s (%d)"),
		argc, fun->name, fun->n_params);
	err = 1;
    }

    for (i=0; i<argc && !err; i++) {
	if (numeric_string(argv[i])) {
	    err = dataset_add_scalar_as(argv[i], fun->params[i], pZ, pdinfo);
#if FN_DEBUG
	    fprintf(stderr, "fn args: new scalar '%s' = %g: err = %d\n",
		    fun->params[i], atof(argv[i]), err);
#endif
	} else {
	    v = varindex(pdinfo, argv[i]);
	    if (v < pdinfo->v) {
		err = dataset_copy_variable_as(v, fun->params[i], pZ, pdinfo);
#if FN_DEBUG
                fprintf(stderr, "fn args: copy var %d as '%s': err = %d\n", 
			v, fun->params[i], err);
#endif
	    } else if (get_list_by_name(argv[i]) != NULL) {
		err = copy_named_list_as(argv[i], fun->params[i]);
#if FN_DEBUG
		fprintf(stderr, "fn args: copy list '%s' as '%s': err = %d\n", 
			argv[i], fun->params[i], err);
#endif
	    } else {
		fprintf(stderr, "argument %d (%s): not a known variable\n", 
			i + 1, argv[i]);
		err = 1;
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
    char fname[32];
    ufunc *fun;
    fncall *call = NULL;
    int argc = 0;
    int assc = 0;
    int err = 0;

    function_name_from_line(line, fname);
    fun = get_ufunc_by_name(fname);

    if (fun == NULL) {
	return 1;
    }

    err = parse_function_args_etc(line, &argc, &argv, &assc, &assv,
				  fun->is_macro);
    if (err) {
	return E_ALLOC;
    }

    if (assc > fun->n_returns) {
	sprintf(gretl_errmsg, _("Number of assignments (%d) exceeds the "
				"number of values returned by\n%s (%d)"), 
		assc, fun->name, fun->n_returns);
	free_strings_array(argv, argc);
	free_strings_array(assv, assc);
	return 1;
    }

    if (!fun->is_macro) {
	if (argc > 0) {
	    err = check_and_allocate_function_args(fun, argc, argv, pZ, pdinfo);
	}
	if (err) {
	    free_strings_array(argv, argc);
	    free_strings_array(assv, assc);
	    return 1;
	}
    }

    call = fncall_new(fun, argc, argv, assc, assv);
    if (call == NULL) {
	return E_ALLOC;
    } 

    err = push_fncall(call, pdinfo);

    if (err) {
	free_fncall(call);
    }

    return err;
}

static int 
safe_strncat (char *targ, const char *src, int n, int maxlen)
{
    if (n == 0) {
	n = strlen(src);
    }

    if (strlen(targ) + n >= maxlen) {
	return 1;
    }

    strncat(targ, src, n);

    return 0;
}

static int 
substitute_dollar_terms (char *targ, const char *src, 
			 int maxlen, int argc, 
			 const char **argv)
{
    int pos, err = 0;

    *targ = '\0';

    while (*src && !err) {
	int len;

	if (strchr(src, '$') == NULL) {
	    err = safe_strncat(targ, src, 0, maxlen);
	    break;
	}

	len = strcspn(src, "$");

	if (len > 0) {
	    err = safe_strncat(targ, src, len, maxlen);
	    if (err) {
		break;
	    }
	    src += len;
	}

	/* got a positional parameter? */
	if (*(src+1) && sscanf(src, "$%d", &pos)) {
	    if (pos <= argc) {
		/* make the substitution */
		err = safe_strncat(targ, argv[pos - 1], 0, maxlen);
	    }
	    src++;
	    while (isdigit((unsigned char) *src)) src++;
	} else {
	    err = safe_strncat(targ, "$", 1, maxlen);
	    src++;
	} 
    }

    return err;
}

char *gretl_function_get_line (char *line, int len,
			       double ***pZ, DATAINFO *pdinfo)
{
    fncall *call = current_call();
    const char *src;

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
	unstack_fncall(pZ, pdinfo);
	return "";
    } 

    src = call->fun->lines[call->lnum];
    if (!strncmp(src, "exit", 4)) {
	/* terminate execution */
	unstack_fncall(pZ, pdinfo);
	return "";
    } 

    call->lnum += 1;

    if (call->fun->is_macro) {
	if (!string_is_blank(src)) {
	    int err = substitute_dollar_terms(line, src, len, call->argc, 
					      (const char **) call->argv); 
	    if (err) {
		sprintf(gretl_errmsg,
			_("Maximum length of command line "
			  "(%d bytes) exceeded\n"), MAXLEN);
		unstack_fncall(pZ, pdinfo);
		return NULL;
	    }
	}
# if FN_DEBUG
	fprintf(stderr, "function_get_line, $substitution: \n"
		" before: '%s'\n  after: '%s'\n", src, line);
# endif
    } else {
	strcpy(line, src);
    }

    return line;
}

void gretl_functions_cleanup (void)
{
    callstack_destroy();
    ufuncs_destroy();
}

void gretl_function_stop_on_error (DATAINFO *pdinfo, PRN *prn)
{
    gretl_function_flagged_error(NULL, prn);

    callstack_destroy();

    if (fn_executing > 0) {
	libset_restore_state_zero(pdinfo);
    }

    fn_executing = 0;
    macro_executing = 0;
}

int gretl_function_flagged_error (const char *s, PRN *prn)
{
    fncall *call;

    if (!fn_executing && !macro_executing) {
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
