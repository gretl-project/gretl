/*
 *  Copyright (c) by Allin Cottrell
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

#define CALLSTACK_DEPTH 8

#define FN_DEBUG 0

typedef struct ufunc_ ufunc;
typedef struct fncall_ fncall;

struct ufunc_ {
    char name[32];
    int n_lines;
    char **lines;
};

struct fncall_ {
    ufunc *fun;
    int lnum;
    int argc;
    char **argv;
};

static int n_ufuns;
static ufunc **ufuns;

static fncall **callstack;

static void free_fncall (fncall *call);

/* record of state, and communication of state with outside world */

static int compiling;
static int executing;

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
    return executing;
}

static void set_executing_on (void)
{
    executing++;
}

static void set_executing_off (void)
{
    executing--;
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

    for (i=0; i<CALLSTACK_DEPTH; i++) {
	if (callstack[i] != NULL) n++;
	else break;
    }

    return n;
}

static int push_fncall (fncall *call)
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

    set_executing_on();

    return 0;
}

static int 
destroy_local_vars (double ***pZ, DATAINFO *pdinfo, int nc)
{
    int i, nlocal = 0;
    int err = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (STACK_LEVEL(pdinfo, i) == nc) {
	    nlocal++;
	}
    }

    if (nlocal > 0) {
	int *locals = malloc((nlocal + 1) * sizeof *locals);

	if (locals == NULL) {
	    err = 1;
	} else {
	    int j = 1;

	    locals[0] = nlocal;
	    for (i=1; i<pdinfo->v; i++) {
		if (STACK_LEVEL(pdinfo, i) == nc) {
#if FN_DEBUG
                    fprintf(stderr, "local variable %d (%s) "
			    "marked for deletion\n", i,
			    pdinfo->varname[i]);
#endif
		    locals[j++] = i;
		}
	    }
	    err = dataset_drop_listed_vars(locals, pZ, pdinfo, NULL);
	}
	free(locals);
    }

    return err;
}

static int unstack_fncall (double ***pZ, DATAINFO *pdinfo)
{
    int i, nc;
    int err = 0;

    if (callstack == NULL) return 1;

    nc = gretl_function_stack_depth();

#if FN_DEBUG
    fprintf(stderr, "unstack_fncall: terminating call to "
	    "function '%s' at depth %d\n", 
	    (callstack[0])->fun->name, nc);
#endif

    free_fncall(callstack[0]);
    err = destroy_local_vars(pZ, pdinfo, nc);

    for (i=0; i<nc; i++) {
	if (i == nc - 1) {
	    callstack[i] = NULL;
	} else {
	    callstack[i] = callstack[i+1];
	}
    }

    set_executing_off();

    return err;
}

static int function_is_on_stack (ufunc *func)
{
    int i;

    for (i=0; i<CALLSTACK_DEPTH; i++) {
	if (callstack[i] == NULL) break;
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

    if (func == NULL) return NULL;

    func->name[0] = '\0';
    func->n_lines = 0;
    func->lines = NULL;

    return func;
}

static fncall *fncall_new (ufunc *fun, int argc, char **argv)
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
	return NULL;
    }

    call->fun = fun;
    call->lnum = 0;
    call->argc = argc;
    call->argv = argv;

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
    int i;

    for (i=0; i<fun->n_lines; i++) {
	free(fun->lines[i]);
    }
    free(fun->lines);

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
    int i;

    for (i=0; i<call->argc; i++) {
	free(call->argv[i]);
    }
    free(call->argv);

    free(call);
}

static ufunc *get_ufunc_by_name (const char *name)
{
    ufunc *fun = NULL;
    int i;

#if FN_DEBUG
    fprintf(stderr, "get_ufunc_by_name: name = '%s' (n_ufuns = %d)\n",
	    name, n_ufuns);
#endif

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(name, (ufuns[i])->name)) {
	    fun = ufuns[i];
	    break;
	}
    }

    return fun;
}

int gretl_is_user_function (const char *s)
{
    int ret = 0;

    if (n_ufuns > 0 && !string_is_blank(s)) {
	char name[32];

	if (sscanf(s, "%31s", name) &&
	    get_ufunc_by_name(name) != NULL) {
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

    /* or should we overwrite? */

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

static char **parse_args (const char *s, int *argc, int *err)
{
    char **argv = NULL;
    int i, na;

    *argc = 0;
    *err = 0;

    /* skip over function name */
    while (*s) {
	if (*s == ' ') break;
	s++;
    }

    /* skip over spaces to args */
    while (*s) {
	if (*s != ' ') break;
	s++;
    } 

    if (*s == '\0') {
	/* got no args */
	return NULL;
    }

    /* count comma-separated arguments */
    na = comma_count(s) + 1;

    argv = malloc(na * sizeof *argv);
    if (argv == NULL) {
	*err = 1;
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
	    *err = 1;
	    break;
	}
	argv[i] = arg;

	s += len;
	while (*s) {
	    if (*s != ',' && *s != ' ') break;
	    s++;
	}
    }

    if (*err) {
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

static int maybe_delete_function (const char *fname)
{
    ufunc *fun = get_ufunc_by_name(fname);
    int err = 0;

    if (fun == NULL) {
	err = 1;
    } else if (function_is_on_stack(fun)) {
	sprintf(gretl_errmsg, "%s: function is in use", fname);
	err = 1;
    } else {
	delete_ufunc_from_list(fun);
    } 

    return err;
}

int gretl_start_compiling_function (const char *line)
{
    char fname[32];
    int name_status;
    ufunc *fun = NULL;

    if (!sscanf(line, "function %31s", fname)) {
	return E_PARSE;
    }

    name_status = check_func_name(fname);

    if (name_status == FN_NAME_BAD) {
	return 1;
    } 

    if (name_status == FN_NAME_TAKEN) {
	if (strstr(line, "delete")) {
	    return maybe_delete_function(fname);
	} else {
	    return 1;
	}
    }

    fun = add_ufunc();
    if (fun == NULL) {
	return E_ALLOC;
    }

    strcpy(fun->name, fname);

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

int gretl_function_append_line (const char *line)
{
    ufunc *fun = get_latest_ufunc();
    char **lines;
    int nl, err = 0;

    if (fun == NULL) return 1;

    if (string_is_blank(line)) {
	return 0;
    }

    if (!strncmp(line, "end ", 4)) {
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

    nl = fun->n_lines;
    lines = realloc(fun->lines, (nl + 1) * sizeof *lines);
    if (lines == NULL) {
	return E_ALLOC;
    }

    fun->lines = lines;

    fun->lines[nl] = gretl_strdup(line);
    if (fun->lines[nl] == NULL) {
	return E_ALLOC;
    }

    fun->n_lines += 1;

    return err;
}

int gretl_function_start_exec (const char *line)
{
    char **argv;
    char fname[32];
    ufunc *fun;
    fncall *call;
    int argc;
    int err = 0;

    sscanf(line, "%31s", fname);
    fun = get_ufunc_by_name(fname);

    if (fun == NULL) {
	return 1;
    }

    argv = parse_args(line + 1, &argc, &err);

    if (err) {
	return E_ALLOC;
    }

    call = fncall_new(fun, argc, argv);
    if (call == NULL) {
	return E_ALLOC;
    } 

    err = push_fncall(call);
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
	if ((err = safe_strncat(targ, src, len, maxlen))) {
	    break;
	}
	src += len;

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
    int err = 0;

    if (call == NULL || call->fun == NULL) {
	return NULL;
    }

    if (call->lnum > call->fun->n_lines - 1) {
	/* finished executing */
	unstack_fncall(pZ, pdinfo);
	return NULL;
    } 

    src = call->fun->lines[call->lnum];
    if (!strncmp(src, "exit", 4)) {
	/* terminate execution */
	unstack_fncall(pZ, pdinfo);
	return NULL;
    } 

    call->lnum += 1;

    if (!string_is_blank(src)) {
	err = substitute_dollar_terms(line, src, len, call->argc, 
				      (const char **) call->argv); 
    }

    if (err) {
	sprintf(gretl_errmsg,
		_("Maximum length of command line "
		  "(%d bytes) exceeded\n"), MAXLEN);
	unstack_fncall(pZ, pdinfo);
	return NULL;
    }

#if FN_DEBUG
	fprintf(stderr, "function_get_line, $substitution: \n"
		" before: '%s'\n  after: '%s'\n", src, line);
#endif

    return line;
}

void gretl_functions_cleanup (void)
{
    callstack_destroy();
    ufuncs_destroy();
}

void gretl_function_error (void)
{
    callstack_destroy();
    executing = 0;
}
