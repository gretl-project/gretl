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

/* monte_carlo.c - loop procedures */

#include "libgretl.h" 
#include "monte_carlo.h"
#include "libset.h"
#include "compat.h"
#include "cmd_private.h"
#include "var.h"
#include "objstack.h"
#include "gretl_func.h"
#include "uservar.h"
#include "flow_control.h"
#include "system.h"
#include "genparse.h"
#include "gretl_string_table.h"

#include <time.h>
#include <unistd.h>

#define LOOP_DEBUG 0
#define SUBST_DEBUG 0

#include <gmp.h>

typedef mpf_t bigval;

enum loop_types {
    COUNT_LOOP,
    WHILE_LOOP,
    INDEX_LOOP,
    DATED_LOOP,
    FOR_LOOP,
    EACH_LOOP
};

#define DEFAULT_NOBS 512

#define indexed_loop(l) (l->type == INDEX_LOOP || \
                         l->type == DATED_LOOP || \
			 l->type == EACH_LOOP)

typedef struct {
    int lineno;    /* location: line number in loop */
    int n;         /* number of repetitions */
    int nvars;     /* number of variables */
    char **names;  /* names of vars to print */
    bigval *sum;   /* running sum of values */
    bigval *ssq;   /* running sum of squares */
    double *xbak;  /* previous values */
    int *diff;     /* indicator for difference */
    char *na;      /* indicator for NAs in calculation */
} LOOP_PRINT;  

/* below: used only in "progressive" loops */ 

typedef struct {
    int lineno;             /* location: line number in loop */
    int n;                  /* number of repetitions */
    int nc;                 /* number of coefficients */
    MODEL *model0;          /* copy of initial model */
    bigval *bigarray;       /* global pointer array */
    bigval *sum_coeff;      /* sums of coefficient estimates */
    bigval *ssq_coeff;      /* sums of squares of coeff estimates */
    bigval *sum_sderr;      /* sums of estimated std. errors */
    bigval *ssq_sderr;      /* sums of squares of estd std. errs */
    double *cbak;           /* previous values of coeffs */
    double *sbak;           /* previous values of std. errs */
    int *cdiff;             /* indicator for difference in coeff */
    int *sdiff;             /* indicator for difference in s.e. */
} LOOP_MODEL;

typedef struct {
    int lineno;     /* location: line number in loop */ 
    int n;          /* number of observations */
    int nvars;      /* number of variables to store */
    char **names;   /* names of vars to print */
    char *fname;    /* filename for output */
    gretlopt opt;   /* formatting option */
    DATASET *dset;  /* temporary data storage */
} LOOP_STORE;

enum loop_flags {
    LOOP_PROGRESSIVE = 1 << 0,
    LOOP_VERBOSE     = 1 << 1,
    LOOP_QUIET       = 1 << 2,
    LOOP_DELVAR      = 1 << 3
};

struct controller_ {
    double val;            /* evaluated value */
    char vname[VNAMELEN];  /* name of (scalar) variable, if used */
    int vsign;             /* 1 or -1, if vname is used */
    char *expr;            /* expression to pass to genr, if used */
};

typedef struct controller_ controller;

enum loop_command_codes {
    LOOP_CMD_GENR    = 1 << 0, /* compiled "genr" */
    LOOP_CMD_LIT     = 1 << 1, /* literal printing */
    LOOP_CMD_NODOL   = 1 << 2, /* no $-substitution this line */
    LOOP_CMD_NOSUB   = 1 << 3, /* no @-substitution this line */
    LOOP_CMD_NOOPT   = 1 << 4, /* no option flags in this line */
    LOOP_CMD_CATCH   = 1 << 5, /* "catch" flag present */
    LOOP_CMD_COND    = 1 << 6  /* compiled conditional */
};

struct loop_command_ {
    char *line;
    int ci;
    char flags;
    GENERATOR *genr;
};

typedef struct loop_command_ loop_command;

typedef struct LOOPSET_ LOOPSET;

struct LOOPSET_ {
    /* basic characteristics */
    char type;
    char flags;
    int level;
    int err;

    /* iterations */
    int itermax;
    int iter;
    int index;

    /* index/foreach control variables */
    char idxname[VNAMELEN];
    int idxval;
    char listname[VNAMELEN];

    /* break signal */
    char brk;

    /* control structures */
    controller init;
    controller test;
    controller delta;
    controller final;

    /* numbers of various subsidiary objects */
    int n_cmds;
    int n_models;
    int n_loop_models;
    int n_prints;
    int n_children;

    /* subsidiary objects */
    loop_command *cmds;   /* saved command info */
    char **eachstrs;      /* for use with "foreach" loop */
    MODEL **models;       /* regular model pointers */
    int *model_lines;     
    LOOP_MODEL *lmodels;
    LOOP_PRINT *prns;
    LOOP_STORE store;
    LOOPSET *parent;
    LOOPSET **children;
    int parent_line;
};

#define loop_is_progressive(l)  (l->flags & LOOP_PROGRESSIVE)
#define loop_set_progressive(l) (l->flags |= LOOP_PROGRESSIVE)
#define loop_is_verbose(l)      (l->flags & LOOP_VERBOSE)
#define loop_set_verbose(l)     (l->flags |= LOOP_VERBOSE)
#define loop_is_quiet(l)        (l->flags & LOOP_QUIET)
#define loop_set_quiet(l)       (l->flags |= LOOP_QUIET)

#define model_print_deferred(o) (o & OPT_F)

static void controller_init (controller *clr);
static int gretl_loop_prepare (LOOPSET *loop);
static void loop_model_free (LOOP_MODEL *lmod);
static void loop_print_free (LOOP_PRINT *lprn);
static void loop_store_free (LOOP_STORE *lstore);
static int extend_loop_dataset (LOOP_STORE *lstore);
static void controller_free (controller *clr);

static int 
make_dollar_substitutions (char *str, int maxlen,
			   const LOOPSET *loop,
			   const DATASET *dset,
			   int *subst,
			   gretlopt opt);

#define LOOP_BLOCK 32

/* record of state, and communication of state with outside world */

static LOOPSET *currloop;

static int compile_level;
static int loop_execute;

int gretl_compiling_loop (void)
{
    return compile_level;
}

int gretl_execute_loop (void)
{
    return loop_execute;
}

/* For indexed loops: get a value from a loop "limit" element (lower
   or upper).  If we got the name of a scalar variable at setup time,
   look up its current value (and modify the sign if wanted).  Or if
   we got a "genr" expression, evaluate it.  Otherwise we should have
   got a numerical constant at setup, in which case we just return
   that value.  
*/

static double controller_get_val (controller *clr,
				  LOOPSET *loop,
				  DATASET *dset, 
				  int *err)
{
    if (clr->vname[0] != '\0') {
	if (gretl_is_scalar(clr->vname)) {
	    clr->val = gretl_scalar_get_value(clr->vname, NULL) * clr->vsign;
	} else {
	    gretl_errmsg_sprintf(_("'%s': not a scalar"), clr->vname);
	} 
    } else if (clr->expr != NULL) {
	int done = 0;

	if (strchr(clr->expr, '@')) {
	    /* the expression needs string substitution? */
	    int subst = 0;
	    char expr[32];

	    *expr = '\0';
	    strncat(expr, clr->expr, 31);
	    *err = substitute_named_strings(expr, &subst);
	    if (!*err && subst) {
		clr->val = generate_scalar(expr, dset, err);
		done = 1;
	    }
	}
	if (!done && !*err && strchr(clr->expr, '$')) {
	    /* the expression needs dollar substitution? */
	    int subst = 0;
	    char expr[32];

	    *expr = '\0';
	    strncat(expr, clr->expr, 31);
	    *err = make_dollar_substitutions(expr, 31, loop, dset, 
					     &subst, OPT_T);
	    if (!*err && subst) {
		clr->val = generate_scalar(expr, dset, err);
		done = 1;
	    }	    
	}
	if (!*err && !done) {
	    clr->val = generate_scalar(clr->expr, dset, err);
	}
    }

#if LOOP_DEBUG
    fprintf(stderr, "controller_get_val: vname='%s', expr='%s', val=%g, err=%d\n", 
	    clr->vname, clr->expr, clr->val, *err);
#endif

    return clr->val;
}

/* apply initialization in case of for-loop */

static void
forloop_init (LOOPSET *loop, DATASET *dset, int *err)
{
    const char *expr = loop->init.expr;

    if (expr != NULL) {
	*err = generate(expr, dset, OPT_Q, NULL);
	if (*err) {
	    gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
				 expr);
	}
    }
}

/* evaluate boolean condition in for-loop or while-loop */

static int 
loop_testval (LOOPSET *loop, DATASET *dset, int *err)
{
    const char *expr = loop->test.expr;
    int ret = 1;

    if (expr != NULL) {
	double x = generate_scalar(expr, dset, err);

	if (!*err && na(x)) {
	    *err = E_DATA;
	    ret = 0;
	} else {
	    ret = x;
	}
	if (*err) {
	    gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
				 expr);
	}
    }

    return ret;
}

/* evaluate third expression in for-loop, if any */

static void 
loop_delta (LOOPSET *loop, DATASET *dset, int *err)
{
    const char *expr = loop->delta.expr;

    if (expr != NULL) {
	*err = generate(expr, dset, OPT_Q, NULL);
	if (*err) {
	    gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
				 expr);
	}
    }
}

static gretlopt get_loop_opts (ExecState *s, int *err)
{
    gretlopt opt = OPT_NONE;

    if (compile_level > 0) {
	/* options will not have been parsed out already */
	opt = get_gretl_options(s->line, err);
    } else {
	opt = s->cmd->opt;
    }

    return opt;
}

static void set_loop_opts (LOOPSET *loop, gretlopt opt)
{
    if (opt & OPT_P) {
	loop_set_progressive(loop);
    }
    if (opt & OPT_V) {
	loop_set_verbose(loop);
    }
    if (opt & OPT_Q) {
	loop_set_quiet(loop);
    }
}

#define plain_model_ci(c) (MODEL_COMMAND(c) && \
                           c != NLS && \
                           c != MLE && \
                           c != GMM)

/**
 * ok_in_loop:
 * @ci: command index.
 *
 * Returns: 1 if the given command is acceptable inside the loop construct,
 * 0 otherwise.
 */

int ok_in_loop (int c)
{
    /* here are the commands we _don't_ currently allow */

    if (c == CORRGM ||
	c == CUSUM ||
	c == FUNC ||
	c == HURST ||
	c == INCLUDE ||
	c == LEVERAGE ||
	c == NULLDATA ||
	c == RMPLOT ||
	c == RUN ||
	c == SETMISS ||
	c == VIF)  {
	return 0;
    }

    return 1;
}

static int loop_attach_child (LOOPSET *loop, LOOPSET *child)
{
    LOOPSET **children;
    int nc = loop->n_children;

    children = realloc(loop->children, (nc + 1) * sizeof *children);
    if (children == NULL) {
	return E_ALLOC;
    } 

    loop->children = children;
    loop->children[nc] = child;
    child->parent = loop;
    child->parent_line = loop->n_cmds;
    child->level = loop->level + 1;

#if LOOP_DEBUG
    fprintf(stderr, "child loop %p has parent %p\n", 
	    (void *) child, (void *) child->parent);
#endif

    loop->n_children += 1;

    return 0;
}

static void loop_store_init (LOOP_STORE *lstore)
{
    lstore->lineno = -1;
    lstore->n = 0;
    lstore->nvars = 0;
    lstore->names = NULL;
    lstore->fname = NULL;
    lstore->opt = OPT_NONE;
    lstore->dset = NULL;
}

static void gretl_loop_init (LOOPSET *loop)
{
#if LOOP_DEBUG
    fprintf(stderr, "gretl_loop_init: initing loop at %p\n", (void *) loop);
#endif

    loop->flags = 0;
    loop->level = 0;

    loop->itermax = 0;
    loop->iter = 0;
    loop->err = 0;
    *loop->idxname = 0;
    loop->idxval = 0;
    loop->brk = 0;
    *loop->listname = '\0';

    controller_init(&loop->init);
    controller_init(&loop->test);
    controller_init(&loop->delta);
    controller_init(&loop->final);

    loop->n_cmds = 0;
    loop->n_models = 0;
    loop->n_loop_models = 0;
    loop->n_prints = 0;

    loop->cmds = NULL;
    loop->model_lines = NULL;

    loop->eachstrs = NULL;

    loop->models = NULL;
    loop->lmodels = NULL;
    loop->prns = NULL;

    loop_store_init(&loop->store);

    loop->parent = NULL;
    loop->children = NULL;
    loop->n_children = 0;
    loop->parent_line = 0;
}

static LOOPSET *gretl_loop_new (LOOPSET *parent)
{
    LOOPSET *loop = malloc(sizeof *loop);

    if (loop == NULL) {
	return NULL;
    }

    gretl_loop_init(loop);

    if (parent != NULL) {
	int err = loop_attach_child(parent, loop);

	if (err) {
	    free(loop);
	    loop = NULL;
	} 
    }
	
    return loop;
}

static void gretl_loop_destroy (LOOPSET *loop)
{
    int i;

    if (loop == NULL) {
	return;
    }

    for (i=0; i<loop->n_children; i++) {
	gretl_loop_destroy(loop->children[i]);
	loop->children[i] = NULL;
    }

    controller_free(&loop->init);
    controller_free(&loop->test);
    controller_free(&loop->delta);
    controller_free(&loop->final);

    if (loop->cmds != NULL) {
	for (i=0; i<loop->n_cmds; i++) {
	    free(loop->cmds[i].line);
	    if (loop->cmds[i].genr != NULL) {
		destroy_genr(loop->cmds[i].genr);
	    }	    
	}
	free(loop->cmds);
    }

    free(loop->model_lines);
    free(loop->models);

    if (loop->eachstrs != NULL) {
	strings_array_free(loop->eachstrs, loop->itermax);
    }

    if (loop->lmodels != NULL) {
	for (i=0; i<loop->n_loop_models; i++) {
#if LOOP_DEBUG
	    fprintf(stderr, "freeing loop->lmodels[%d]\n", i);
#endif
	    loop_model_free(&loop->lmodels[i]);
	}
	free(loop->lmodels);
    }

    if (loop->prns != NULL) {
	for (i=0; i<loop->n_prints; i++) { 
	    loop_print_free(&loop->prns[i]);
	}
	free(loop->prns);
    }

    loop_store_free(&loop->store);

    if (loop->children != NULL) {
	free(loop->children);
    }

    if (loop->flags & LOOP_DELVAR) {
	user_var_delete_by_name(loop->idxname, NULL);
    }

    free(loop);
}

static int parse_as_while_loop (LOOPSET *loop, const char *s)
{
    int err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_while_loop: cond = '%s'\n", s);
#endif

    if (s == NULL || *s == '\0') {
	err = E_PARSE;
    } else {
	loop->type = WHILE_LOOP;
	loop->test.expr = gretl_strdup(s);
	if (loop->test.expr == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static int loop_attach_index_var (LOOPSET *loop, const char *vname,
				  DATASET *dset)
{
    int err = 0;

    if (gretl_is_scalar(vname)) {
	strcpy(loop->idxname, vname);
	gretl_scalar_set_value_authorized(vname, loop->init.val);
    } else {
	char genline[64];
	
	if (na(loop->init.val)) {
	    sprintf(genline, "scalar %s = NA", vname);
	} else {
	    gretl_push_c_numeric_locale();
	    sprintf(genline, "scalar %s = %g", vname, loop->init.val);
	    gretl_pop_c_numeric_locale();
	}
	err = generate(genline, dset, OPT_Q, NULL);
	if (!err) {
	    strcpy(loop->idxname, vname);
	    loop->flags |= LOOP_DELVAR;
	}
    } 

    return err;
}

/* for a loop control expression such as "j=start..end", get the
   initial or final value from the string @s (we also use this to get
   the count for a simple count loop).
*/

static int index_get_limit (LOOPSET *loop, controller *clr, const char *s, 
			    DATASET *dset)
{
    int v, err = 0;

    if (integer_string(s)) {
	/* plain numerical value */
	clr->val = atoi(s);
    } else {
	if (*s == '-') {
	    /* negative of variable? */
	    clr->vsign = -1;
	    s++;
	}
	if (gretl_is_scalar(s)) {
	    *clr->vname = '\0';
	    strncat(clr->vname, s, VNAMELEN - 1);
	    clr->val = (int) gretl_scalar_get_value(s, NULL);
	} else if ((v = current_series_index(dset, s)) >= 0) {
	    /* found a series by the name of s */
	    gretl_errmsg_sprintf(_("'%s': not a scalar"), s);
	} else if (loop->parent != NULL && strlen(s) == gretl_namechar_spn(s)) {
	    /* potentially valid varname, but unknown at present */
	    *clr->vname = '\0';
	    strncat(clr->vname, s, VNAMELEN - 1);
	} else {
	    /* expression to be evaluated to scalar? */
	    clr->expr = gretl_strdup(s);
	    if (clr->expr == NULL) {
		err = E_ALLOC;
	    }
	} 
    }

    return err;
}

#define maybe_date(s) (strchr(s, ':') || strchr(s, '/'))

static int parse_as_indexed_loop (LOOPSET *loop,
				  DATASET *dset,
				  const char *lvar, 
				  const char *start,
				  const char *end)
{
    int err = 0;

    /* starting and ending values: the order in which we try 
       for valid values is: dates, numeric constants,
       named scalars, scalar expressions.
    */

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_indexed_loop: start='%s', end='%s'\n", start, end);
#endif

    if (maybe_date(start)) {
	loop->init.val = dateton(start, dset);
	if (loop->init.val < 0) {
	    err = E_DATA;
	} else {
	    loop->init.val += 1;
	    loop->final.val = dateton(end, dset);
	    if (loop->final.val < 0) {
		err = E_DATA;
	    } else {
		loop->final.val += 1;
		loop->type = DATED_LOOP;
	    }
	}
    } else {
	err = index_get_limit(loop, &loop->init, start, dset);
	if (!err) {
	    err = index_get_limit(loop, &loop->final, end, dset);
	}
	if (!err) {
	    loop->type = INDEX_LOOP;
	}
    }

    if (!err) {
	err = loop_attach_index_var(loop, lvar, dset);
    }

#if LOOP_DEBUG
    fprintf(stderr, "indexed_loop: init.val=%g, final.val=%g, err=%d\n",
	    loop->init.val, loop->final.val, err);
#endif

    return err;
}

static int parse_as_count_loop (LOOPSET *loop, 
				DATASET *dset,
				const char *s)
{
    int err;

    err = index_get_limit(loop, &loop->final, s, dset);

    if (!err) {
	loop->init.val = 1;
	loop->type = COUNT_LOOP;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_count_loop: init.val=%g, final.val=%g\n",
	    loop->init.val, loop->final.val);
#endif

    return err;
}

static int set_forloop_element (char *s, LOOPSET *loop, int i)
{
    controller *clr = (i == 0)? &loop->init :
	(i == 1)? &loop->test : &loop->delta;
    int len, err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "set_forloop_element: i=%d: '%s'\n", i, s);
#endif

    if (s == NULL || *s == '\0') {
	/* an empty "for" field */
	if (i == 1) {
	    /* test is implicitly always true */
	    clr->val = 1;
	} else {
	    /* no-op */
	    clr->val = 0;
	}
	return 0;
    }

    clr->expr = gretl_strdup(s);
    if (clr->expr == NULL) {
	err = E_ALLOC;
    }

    if (!err && i == 0) {
	/* initialization: look for varname for possible substitution */
	err = extract_varname(clr->vname, s, &len);
    } 

#if LOOP_DEBUG
    fprintf(stderr, " expr='%s', vname='%s'\n", clr->expr, clr->vname);
#endif

    return err;
}

static int allocate_each_strings (LOOPSET *loop, int n)
{
    loop->eachstrs = strings_array_new(n);

    return (loop->eachstrs == NULL)? E_ALLOC : 0;
}

static int list_vars_to_strings (LOOPSET *loop, const int *list,
				 const DATASET *dset)
{
    int i, vi;
    int err;

#if LOOP_DEBUG
    fprintf(stderr, "list_vars_to_strings: adding %d strings\n", list[0]);
#endif

    err = allocate_each_strings(loop, list[0]);

    for (i=0; i<list[0] && !err; i++) {
	vi = list[i+1];
	if (vi < 0 || vi >= dset->v) {
	    err = E_DATA;
	} else {
	    loop->eachstrs[i] = gretl_strdup(dset->varname[vi]);
	    if (loop->eachstrs[i] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

/* At loop runtime, check the named list and insert the names (or
   numbers) of the variables as "eachstrs"; flag an error if the list
   has disappeared. We also have to handle the case where the name
   of the loop-controlling list is subject to $-substitution.
*/

static int loop_list_refresh (LOOPSET *loop, const DATASET *dset)
{
    int *list = NULL;
    int err = 0;

    if (strchr(loop->listname, '$') != NULL) {
	/* $-string substitution required */
	char lname[VNAMELEN];

	strcpy(lname, loop->listname);
	err = make_dollar_substitutions(lname, VNAMELEN, loop, 
					dset, NULL, OPT_T);
	if (!err) {
	    list = get_list_by_name(lname);
	}
    } else if (*loop->listname == '@') {
	/* @-string substitution required */
	const char *s = get_string_by_name(loop->listname + 1);

	if (s != NULL && strlen(s) < VNAMELEN) {
	    list = get_list_by_name(s);
	}
    } else {
	/* no string substitution needed */
	list = get_list_by_name(loop->listname);
    }

    if (loop->eachstrs != NULL) {
	strings_array_free(loop->eachstrs, loop->itermax);
	loop->eachstrs = NULL;
    }

    loop->itermax = loop->final.val = 0;

    if (list == NULL) {
	if (!err) {
	    err = E_UNKVAR;
	}
    } else if (list[0] > 0) {
	err = list_vars_to_strings(loop, list, dset);
	if (!err) {
	    loop->final.val = list[0];
	}
    }

    return err;
}

static int find_list_in_parentage (LOOPSET *loop, const char *s)
{
    char fmt[16], lname[VNAMELEN];
    int i;

    sprintf(fmt, "list %%%d[^ =]", VNAMELEN-1);

    while ((loop = loop->parent) != NULL) {
	for (i=0; i<loop->n_cmds; i++) {
	    if (sscanf(loop->cmds[i].line, fmt, lname)) {
		if (!strcmp(lname, s)) {
		    return 1;
		}
	    }
	}
    }

    return 0;
}

/* We're looking at a "foreach" loop with just one field after the
   index variable, so it's most likely a loop over a list.

   We begin by looking for a currently existing named list, but if
   this fails we don't give up immediately.  If we're working on an
   embedded loop, the list may be created within a parent loop whose
   commands have not yet been executed, so we search upward among the
   ancestors of this loop (if any) for a relevant list-creation
   command.

   Even if we find an already-existing list, we do not yet fill out
   the variable-name (or variable-number) strings: these will be set
   when the loop is actually run, since the list may have changed in
   the meantime.
*/

static int list_loop_setup (LOOPSET *loop, char *s, int *nf)
{
    int *list;
    int err = 0;

    while (isspace(*s)) s++;
    tailstrip(s);

    if (*s == '@') {
	/* tricksy: got a list-name that needs string subst? */
	*loop->listname = '\0';
	strncat(loop->listname, s, VNAMELEN - 1);
	*nf = 0;
	return 0;
    }

    list = get_list_by_name(s);

#if LOOP_DEBUG
    fprintf(stderr, "list_loop_setup: s = '%s'\n", s);
    printlist(list, "get_list_by_name");
#endif

    if (list == NULL && !find_list_in_parentage(loop, s)) {
	err = E_UNKVAR;
    } else {
	*loop->listname = '\0';
	strncat(loop->listname, s, VNAMELEN - 1);
	*nf = (list != NULL)? list[0] : 0;
    } 

    return err;
}

enum {
    DOTTED_LIST,
    WILDCARD_LIST
};

static int
each_strings_from_list_of_vars (LOOPSET *loop, const DATASET *dset, 
				char *s, int *pnf, int type)
{
    int *list = NULL;
    int err = 0;

    if (type == WILDCARD_LIST) {
	s += strspn(s, " \t");
	list = varname_match_list(dset, s, &err);
    } else {
	char vn1[VNAMELEN], vn2[VNAMELEN];
	char fmt[16];
	
	gretl_delchar(' ', s);
	sprintf(fmt, "%%%d[^.]..%%%ds", VNAMELEN-1, VNAMELEN-1);
	

	if (sscanf(s, fmt, vn1, vn2) != 2) {
	    err = E_PARSE;
	} else {
	    int v1 = current_series_index(dset, vn1);
	    int v2 = current_series_index(dset, vn2);

	    if (v1 < 0 || v2 < 0) {
		err = E_UNKVAR;
	    } else if (v2 - v1 + 1 <= 0) {
		err = E_DATA;
	    } else {
		list = gretl_consecutive_list_new(v1, v2);
		if (list == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
	if (err) {
	    *pnf = 0;
	}
    }

    if (list != NULL) {
	int i, vi;

	err = allocate_each_strings(loop, list[0]);
	if (!err) {
	    for (i=1; i<=list[0] && !err; i++) {
		vi = list[i];
		loop->eachstrs[i-1] = gretl_strdup(dset->varname[vi]);
		if (loop->eachstrs[i-1] == NULL) {
		    strings_array_free(loop->eachstrs, list[0]);
		    loop->eachstrs = NULL;
		    err = E_ALLOC;
		}
	    }
	}
	if (!err) {
	    *pnf = list[0];
	}
	free(list);
    }

    return err;
}

/* in context of "foreach" loop, split a string variable by
   both spaces and newlines */

static int count_each_fields (const char *s)
{
    int nf = 0;

    if (s != NULL && *s != '\0') {
	const char *p;

	s += strspn(s, " ");

	if (*s != '\0' && *s != '\n') {
	    s++;
	    nf++;
	}

	while (*s) {
	    p = strpbrk(s, " \n");
	    if (p != NULL) {
		s = p + strspn(p, " \n");
		if (*s) {
		    nf++;
		}
	    } else {
		break;
	    }
	}
    }
	    
    return nf;
}

static int
parse_as_each_loop (LOOPSET *loop, DATASET *dset, char *s)
{
    char ivar[VNAMELEN] = {0};
    int done = 0;
    int nf, err = 0;

    /* we're looking at the string that follows "loop foreach" */
    if (*s == '\0') {
	return E_PARSE;
    }
    
    s += strspn(s, " "); /* skip any spaces */

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: s = '%s'\n", s);
#endif 

    /* get the index variable name (as in "foreach i") */
    if (gretl_scan_varname(s, ivar) != 1) {
	return E_PARSE;
    } 

    s += strlen(ivar);
    nf = count_each_fields(s);

#if LOOP_DEBUG
    fprintf(stderr, " number of fields = %d\n", nf);
#endif 

    if (nf == 0) {
	return E_PARSE;
    }

    if (nf <= 3 && strstr(s, "..") != NULL) {
	/* range of values, foo..quux */
	err = each_strings_from_list_of_vars(loop, dset, s, &nf,
					     DOTTED_LIST);
	done = 1;
    } else if (nf == 1 && strchr(s, '*')) {
	err = each_strings_from_list_of_vars(loop, dset, s, &nf,
					     WILDCARD_LIST);
	done = (err == 0);
    }	

    if (!done && nf == 1) {
	/* try for a named list? */
	err = list_loop_setup(loop, s, &nf);
	done = (err == 0);
    }

    if (!done) {
	/* simple array of strings: allow for quoted substrings */
	loop->eachstrs = gretl_string_split_quoted(s, &nf, NULL, &err);
    }

    if (!err) {
	loop->type = EACH_LOOP;
	loop->init.val = 1;
	loop->final.val = nf;
	err = loop_attach_index_var(loop, ivar, dset);
    }   

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: final.val=%g\n", loop->final.val);
#endif 

    return err;
}

/* try to parse out (expr1; expr2; expr3) */

static int parse_as_for_loop (LOOPSET *loop, char *s)
{
    char *tmp, *q;
    int i, j, len;
    int sc = 0;
    int err = 0;

    s += strcspn(s, "(");
    if (*s != '(') {
	return E_PARSE;
    }

    s++;
    q = strrchr(s, ')');
    if (q == NULL) {
	return E_PARSE;
    }

    len = q - s;
    if (len < 2) { /* minimal OK string is ";;" */
	return E_PARSE;
    }

    tmp = malloc(len + 1);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    for (j=0; j<3 && s!=q && !err; j++) {
	/* make a compressed copy of field j */
	i = 0;
	while (s != q) {
	    if (*s == ';') {
		sc++;
		s++;
		break; /* onto next field */
	    }
	    if (*s != ' ') {
		tmp[i++] = *s;
	    }
	    s++;
	}
	tmp[i] = '\0';
	err = set_forloop_element(tmp, loop, j);
    }  

    if (!err && (sc != 2 || s != q)) {
	/* we've reached the reached rightmost ')' but have not
	   found two semi-colons */
	err = E_PARSE;
    }

    free(tmp);

    if (!err) {
	loop->type = FOR_LOOP;
    }

    return err;
}

static int parse_first_loopline (char *s, LOOPSET *loop, 
				 DATASET *dset)
{
    char lvar[VNAMELEN], rvar[VNAMELEN], op[VNAMELEN];
    char fmt[32];
    int err = 0;

    /* skip preliminary string */
    while (isspace(*s)) s++;
    if (!strncmp(s, "loop", 4)) {
	s += 4;
	while (isspace(*s)) s++;
    }

    /* syntactic slop: accept "for i=lo..hi" -> "i=lo..hi" */
    if (!strncmp(s, "for ", 4) && !strchr(s, ';')) {
	s += 4;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_first_loopline: '%s'\n", s);
#endif

    sprintf(fmt, "%%%d[^= ] = %%%d[^.]..%%%ds", VNAMELEN-1, VNAMELEN-1, 
	    VNAMELEN-1);

    if (!strncmp(s, "foreach ", 8)) {
	err = parse_as_each_loop(loop, dset, s + 8);
    } else if (sscanf(s, fmt, lvar, op, rvar) == 3) {
	err = parse_as_indexed_loop(loop, dset, lvar, op, rvar);
    } else if (!strncmp(s, "for", 3)) {
	err = parse_as_for_loop(loop, s + 4);
    } else if (!strncmp(s, "while", 5)) {
	err = parse_as_while_loop(loop, s + 6);
    } else if (gretl_scan_varname(s, lvar) == 1) {
	err = parse_as_count_loop(loop, dset, lvar);
    } else {
	printf("parse_first_loopline: failed on '%s'\n", s);
	gretl_errmsg_set(_("No valid loop condition was given."));
	err = 1;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_first_loopline: returning %d\n", err);
#endif

    return err;
}

/**
 * start_new_loop:
 * @s: loop specification line.
 * @inloop: current loop struct pointer, or %NULL.
 * @dset: dataset struct.
 * @opt: options associated with new loop.
 * @nested: location to receive info on whether a new
 * loop was created, nested within the input loop.
 * @err: location to receive error code.
 *
 * Create a new LOOPSET based on the input line; this may or
 * may not be a child of @inloop.
 *
 * Returns: loop pointer on successful completion, %NULL on error.
 */

static LOOPSET *start_new_loop (char *s, LOOPSET *inloop, 
				DATASET *dset,
				gretlopt opt,
				int *nested,
				int *err)
{
    LOOPSET *loop = NULL;

    gretl_error_clear();

#if LOOP_DEBUG
    fprintf(stderr, "start_new_loop: inloop=%p, line='%s'\n", 
	    (void *) inloop, s);
#endif

#if 0
    if (inloop != NULL && (opt & OPT_P)) {
	/* don't allow nesting of "progressive" loops? */
	fprintf(stderr, "inner loop with OPT_P\n");
    }
#endif

    if (inloop == NULL || compile_level <= inloop->level) {
	loop = gretl_loop_new(NULL);
    } else {
	loop = gretl_loop_new(inloop);
	*nested = 1;
    } 

    if (loop == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	*err = E_ALLOC;
	return NULL;
    }

#if LOOP_DEBUG
    fprintf(stderr, " added loop at %p (%s)\n", (void *) loop,
	    (*nested)? "nested" : "independent");
#endif

    *err = parse_first_loopline(s, loop, dset);

    if (!*err) {
	*err = gretl_loop_prepare(loop);
    }

    if (*err) {
	free(loop->cmds);
	free(loop);
	loop = NULL;
    } 

    return loop;
}

#if LOOP_DEBUG
# define MAX_FOR_TIMES  10
#else
# define MAX_FOR_TIMES  5000000
#endif

static int loop_count_too_high (LOOPSET *loop)
{
    int nt = loop->iter + 1;

    if (loop->type == FOR_LOOP) {
	if (nt > MAX_FOR_TIMES) {
	    gretl_errmsg_sprintf(_("Reached maximum iterations, %d"),
				 MAX_FOR_TIMES);
	    loop->err = 1;
	}
    } else {
	int maxit = libset_get_int(LOOP_MAXITER);

	if (maxit > 0 && nt > maxit) {
	    gretl_errmsg_sprintf(_("Reached maximum iterations, %d"),
				 maxit);
	    loop->err = 1;
	}
    }

    return loop->err;
}

/**
 * loop_condition:
 * @loop: pointer to loop commands struct.
 * @dset: data information struct.
 * @err: location to receive error code.
 *
 * Check whether a loop continuation condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

static int 
loop_condition (LOOPSET *loop, DATASET *dset, int *err)
{
    int ok = 0;

    if (loop->brk) {
	/* got "break" comand */
	loop->brk = 0;
	ok = 0;
    } else if (loop->type == COUNT_LOOP || indexed_loop(loop)) {
	if (loop->iter < loop->itermax) {
	    ok = 1;
	    if (indexed_loop(loop) && loop->iter > 0) {
		loop->idxval += 1;
		gretl_scalar_set_value_authorized(loop->idxname, loop->idxval);
	    }
	}
    } else if (!loop_count_too_high(loop)) {
	/* more complex forms of control (for, while) */
	if (loop->type == FOR_LOOP) {
	    if (loop->iter > 0) {
		loop_delta(loop, dset, err);
	    }
	    ok = loop_testval(loop, dset, err);
	} else if (loop->type == WHILE_LOOP) {
	    ok = loop_testval(loop, dset, err);
	}
    }

    return ok;
}

static void controller_init (controller *clr)
{
    clr->val = NADBL;
    clr->vname[0] = 0;
    clr->vsign = 1;
    clr->expr = NULL;
}

static void controller_free (controller *clr)
{
    if (clr->expr != NULL) {
	free(clr->expr);
	clr->expr = NULL;
    }
}

static void loop_cmds_init (LOOPSET *loop, int i1, int i2)
{
    int i;

    for (i=i1; i<i2; i++) {
	loop->cmds[i].line = NULL;
	loop->cmds[i].ci = 0;
	loop->cmds[i].flags = 0;
	loop->cmds[i].genr = NULL;
    }
}

static int gretl_loop_prepare (LOOPSET *loop)
{
    mpf_set_default_prec(256);

    /* allocate some initial lines/commands for loop */
    loop->cmds = malloc(LOOP_BLOCK * sizeof *loop->cmds); 
    
    if (loop->cmds == NULL) {
	return E_ALLOC;
    } else {
	loop_cmds_init(loop, 0, LOOP_BLOCK);
    } 

    return 0;
}

static void loop_model_free (LOOP_MODEL *lmod)
{
    int i, n;

#if LOOP_DEBUG
    fprintf(stderr, "loop_model_free: lmod at %p, model0 at %p\n",
	    (void *) lmod, (void *) lmod->model0);
#endif

    n = 4 * lmod->model0->ncoeff;

    for (i=0; i<n; i++) {
	mpf_clear(lmod->bigarray[i]);
    }

    free(lmod->bigarray);
    free(lmod->cbak);
    free(lmod->cdiff);

    gretl_model_free(lmod->model0);
}

/* Reset the loop model */

static void loop_model_zero (LOOP_MODEL *lmod, int started)
{
    int i, bnc = 4 * lmod->nc;

#if LOOP_DEBUG
    fprintf(stderr, "loop_model_zero: %p\n", (void *) lmod);
#endif

    for (i=0; i<bnc; i++) {
	if (started) {
	    mpf_set_d(lmod->bigarray[i], 0.0);
	} else {
	    mpf_init(lmod->bigarray[i]);
	}
    }

    for (i=0; i<lmod->nc; i++) {
	lmod->cbak[i] = lmod->sbak[i] = NADBL;
	lmod->cdiff[i] = lmod->sdiff[i] = 0;
    }

    lmod->n = 0;
}

/* Set everything in lmod to 0/null in case of failure */

static void loop_model_init (LOOP_MODEL *lmod, int lno)
{
    lmod->lineno = lno;
    lmod->nc = 0;
    lmod->model0 = NULL;
    lmod->bigarray = NULL;
    lmod->cbak = NULL;
    lmod->cdiff = NULL;
}

/* Start up a LOOP_MODEL struct: copy @pmod into place and
   allocate storage */

static int loop_model_start (LOOP_MODEL *lmod, const MODEL *pmod)
{
    int nc = pmod->ncoeff;
    int err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "init: copying model at %p\n", (void *) pmod);
#endif

    lmod->model0 = gretl_model_copy(pmod);
    if (lmod->model0 == NULL) {
	return E_ALLOC;
    }

    lmod->nc = nc;

    lmod->bigarray = malloc(nc * 4 * sizeof *lmod->bigarray);
    if (lmod->bigarray == NULL) {
	return E_ALLOC;
    }

    lmod->sum_coeff = lmod->bigarray;
    lmod->ssq_coeff = lmod->sum_coeff + nc;
    lmod->sum_sderr = lmod->ssq_coeff + nc;
    lmod->ssq_sderr = lmod->sum_sderr + nc;

    lmod->cbak = malloc(nc * 2 * sizeof *lmod->cbak);
    if (lmod->cbak == NULL) {
	err = E_ALLOC;
    } else {
	lmod->sbak = lmod->cbak + nc;
    }

    if (!err) {
	lmod->cdiff = malloc(nc * 2 * sizeof *lmod->cdiff);
	if (lmod->cdiff == NULL) {
	    err = E_ALLOC;
	} else {
	    lmod->sdiff = lmod->cdiff + nc;
	}
    }

    if (!err) {
	loop_model_zero(lmod, 0);
#if LOOP_DEBUG
	fprintf(stderr, " model copied to %p, returning 0\n", 
		(void *) lmod->model0);
#endif
    }

    if (err) {
	free(lmod->bigarray);
	free(lmod->cbak);
	free(lmod->cdiff);
    }

    return err;
}

static void loop_print_free (LOOP_PRINT *lprn)
{
    int i;

    for (i=0; i<lprn->nvars; i++) {
	mpf_clear(lprn->sum[i]);
	mpf_clear(lprn->ssq[i]);
    }

    strings_array_free(lprn->names, lprn->nvars);

    free(lprn->sum);
    free(lprn->ssq);
    free(lprn->xbak);
    free(lprn->diff);
    free(lprn->na);
}

static void loop_print_zero (LOOP_PRINT *lprn, int started)
{
    int i;

    lprn->n = 0;

    for (i=0; i<lprn->nvars; i++) { 
	if (started) {
	    mpf_set_d(lprn->sum[i], 0.0);
	    mpf_set_d(lprn->ssq[i], 0.0);
	} else {
	    mpf_init(lprn->sum[i]);
	    mpf_init(lprn->ssq[i]);
	}
	lprn->xbak[i] = NADBL;
	lprn->diff[i] = 0;
	lprn->na[i] = 0;
    }
}

/* allocate and initialize @lprn, based on the number of
   elements in @namestr */

static int loop_print_start (LOOP_PRINT *lprn, const char *namestr)
{
    int i, nv;

    lprn->names = gretl_string_split(namestr, &lprn->nvars, NULL);
    if (lprn->names == NULL) {
	return E_ALLOC;
    }

    nv = lprn->nvars;

    for (i=0; i<nv; i++) {
	if (!gretl_is_scalar(lprn->names[i])) {
	    gretl_errmsg_sprintf(_("'%s': not a scalar"), lprn->names[i]);
	    strings_array_free(lprn->names, lprn->nvars);
	    lprn->names = NULL;
	    lprn->nvars = 0;
	    return E_DATA;
	}
    }

    lprn->sum = malloc(nv * sizeof *lprn->sum);
    if (lprn->sum == NULL) goto cleanup;

    lprn->ssq = malloc(nv * sizeof *lprn->ssq);
    if (lprn->ssq == NULL) goto cleanup;

    lprn->xbak = malloc(nv * sizeof *lprn->xbak);
    if (lprn->xbak == NULL) goto cleanup;

    lprn->diff = malloc(nv * sizeof *lprn->diff);
    if (lprn->diff == NULL) goto cleanup;

    lprn->na = malloc(nv);
    if (lprn->na == NULL) goto cleanup;

    loop_print_zero(lprn, 0);

    return 0;

 cleanup:

    strings_array_free(lprn->names, lprn->nvars);
    lprn->names = NULL;
    lprn->nvars = 0;
    
    free(lprn->sum);
    free(lprn->ssq);
    free(lprn->xbak);
    free(lprn->diff);
    free(lprn->na);

    lprn->sum = NULL;
    lprn->ssq = NULL;
    lprn->xbak = NULL;
    lprn->diff = NULL;
    lprn->na = NULL;

    return E_ALLOC;
}

static void loop_print_init (LOOP_PRINT *lprn, int lno)
{
    lprn->lineno = lno;
    lprn->nvars = 0;
    lprn->names = NULL;
    lprn->sum = NULL;
    lprn->ssq = NULL;
    lprn->xbak = NULL;
    lprn->diff = NULL;
    lprn->na = NULL;
}

static LOOP_PRINT *get_loop_print_by_line (LOOPSET *loop, int lno, int *err)
{
    LOOP_PRINT *prns;
    int i, np = loop->n_prints;

    for (i=0; i<np; i++) {
	if (loop->prns[i].lineno == lno) {
	    return &loop->prns[i];
	}
    }

    prns = realloc(loop->prns, (np + 1) * sizeof *prns);
    if (prns == NULL) {
	*err = E_ALLOC;
	return NULL;
    } else {
	loop->prns = prns;
    }

    loop_print_init(&loop->prns[np], lno);
    loop->n_prints += 1;

    return &loop->prns[np];
}

static void loop_store_free (LOOP_STORE *lstore)
{
    destroy_dataset(lstore->dset);
    lstore->dset = NULL;

    strings_array_free(lstore->names, lstore->nvars);
    lstore->nvars = 0;
    lstore->names = NULL;

    free(lstore->fname);
    lstore->fname = NULL;

    lstore->lineno = -1;
    lstore->n = 0;
    lstore->opt = OPT_NONE;
}

static int loop_store_set_filename (LOOP_STORE *lstore, 
				    const char *fname,
				    gretlopt opt)
{
    if (fname == NULL || *fname == '\0') {
	return E_ARGS;
    }

    lstore->fname = gretl_strdup(fname);
    if (lstore->fname == NULL) {
	return E_ALLOC;
    }

    if (opt == OPT_NONE) {
	opt = data_save_opt_from_suffix(lstore->fname);
    }

    lstore->opt = opt;    

    return 0;
}

/* check, allocate and initialize loop data storage */

static int loop_store_start (LOOPSET *loop, const char *names, 
			     const char *fname, gretlopt opt)
{
    LOOP_STORE *lstore = &loop->store;
    int i, n, err = 0;

    if (names == NULL || *names == '\0') {
	gretl_errmsg_set("'store' list is empty");
	return E_DATA;
    }

    lstore->names = gretl_string_split(names, &lstore->nvars, NULL);
    if (lstore->names == NULL) {
	return E_ALLOC;
    }

    err = loop_store_set_filename(lstore, fname, opt);
    if (err) {
	return err;
    }

    n = (loop->itermax > 0)? loop->itermax : DEFAULT_NOBS;

    lstore->dset = create_auxiliary_dataset(lstore->nvars + 1, n, 0);
    if (lstore->dset == NULL) {
	return E_ALLOC;
    }
    
#if LOOP_DEBUG
    fprintf(stderr, "loop_store_init: created sZ, v = %d, n = %d\n",
	    lstore->dset->v, lstore->dset->n);
#endif

    for (i=0; i<lstore->nvars && !err; i++) {
	const char *s = lstore->names[i];

	if (!gretl_is_scalar(s)) {
	    gretl_errmsg_sprintf(_("'%s': not a scalar"), s);
	    err = E_DATA;
	} else {
	    strcpy(lstore->dset->varname[i+1], s);
	}
    }

    return err;
}

static int loop_store_update (LOOPSET *loop, int lno,
			      const char *names, const char *fname,
			      gretlopt opt)
{
    LOOP_STORE *lstore = &loop->store;
    int i, t, err = 0;

    if (lstore->lineno >= 0 && lstore->lineno != lno) {
	gretl_errmsg_set("Only one 'store' command is allowed in a "
			 "progressive loop");
	return E_DATA;
    }

    if (lstore->dset == NULL) {
	/* not started yet */
	err = loop_store_start(loop, names, fname, opt);
	if (err) {
	    return err;
	}
    }

    lstore->lineno = lno;
    t = lstore->n;

    if (t >= lstore->dset->n) {
	if (extend_loop_dataset(lstore)) {
	    err = E_ALLOC;
	}
    }

    for (i=0; i<lstore->nvars && !err; i++) {
	lstore->dset->Z[i+1][t] = 
	    gretl_scalar_get_value(lstore->names[i], &err);
    }

    if (!err) {
	lstore->n += 1;
    }

    return err;
}

/* See if we already have a model recorder in place for the command on
   line @lno of the loop.  If so, return it, else create a new one and
   return it.
*/

static MODEL *get_model_record_by_line (LOOPSET *loop, int lno, int *err)
{
    MODEL **models, *pmod;
    int *modlines;
    int n = loop->n_models;
    int i;

    for (i=0; i<n; i++) {
	if (lno == loop->model_lines[i]) {
	    return loop->models[i];
	}
    }

    modlines = realloc(loop->model_lines, (n + 1) * sizeof *modlines);
    if (modlines == NULL) {
	*err = E_ALLOC;
	return NULL;
    } else {
	loop->model_lines = modlines;
    }

    models = realloc(loop->models, (n + 1) * sizeof *models);
    if (models == NULL) {
	*err = E_ALLOC;
	return NULL;
    } else {
	loop->models = models;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    loop->model_lines[n] = lno;
    pmod->ID = n + 1;
    loop->models[n] = pmod;
    loop->n_models += 1;

    return pmod;
}

int model_is_in_loop (const MODEL *pmod)
{
    LOOPSET *loop = currloop;
    int i;

    while (loop != NULL) {
	for (i=0; i<loop->n_models; i++) {
	    if (pmod == loop->models[i]) {
		return 1;
	    }
	}
	loop = loop->parent;
    }

    return 0;
}

/* See if we already have a LOOP_MODEL in place for the command
   on line @lno of the loop.  If so, return it, else create
   a new LOOP_MODEL and return it.
*/

static LOOP_MODEL *
get_loop_model_by_line (LOOPSET *loop, int lno, int *err)
{
    LOOP_MODEL *lmods;
    int n = loop->n_loop_models;
    int i;

#if LOOP_DEBUG
    fprintf(stderr, "get_loop_model_by_line: loop->n_loop_models = %d\n",
	    loop->n_loop_models);
#endif

    for (i=0; i<n; i++) {
	if (loop->lmodels[i].lineno == lno) {
	    return &loop->lmodels[i];
	}
    }

    lmods = realloc(loop->lmodels, (n + 1) * sizeof *loop->lmodels);
    if (lmods == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    loop->lmodels = lmods;
    loop_model_init(&loop->lmodels[n], lno);
    loop->n_loop_models += 1;

    return &loop->lmodels[n];
}

#define realdiff(x,y) (fabs((x)-(y)) > 2.0e-13)

/* Update the info stored in LOOP_MODEL based on the results in pmod.
   If this is the first use we have to do some allocation first.
*/

static int loop_model_update (LOOP_MODEL *lmod, MODEL *pmod)
{
    mpf_t m;
    int j, err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "loop_model_update: lmod = %p, pmod = %p\n", 
	    (void *) lmod, (void *) pmod);
#endif

    if (lmod == NULL) {
	fprintf(stderr, "loop_model_update: got NULL loop model\n");
	return E_DATA;
    }

    if (lmod->nc == 0) {
	/* not started yet */
	err = loop_model_start(lmod, pmod);
	if (err) {
	    return err;
	}
    } else if (pmod->ncoeff != lmod->nc) {
	gretl_errmsg_set(_("progressive loop: model must be of constant size"));
	return E_DATA;
    }

    mpf_init(m);

    for (j=0; j<pmod->ncoeff; j++) {
	mpf_set_d(m, pmod->coeff[j]);
	mpf_add(lmod->sum_coeff[j], lmod->sum_coeff[j], m); 
	mpf_mul(m, m, m);
	mpf_add(lmod->ssq_coeff[j], lmod->ssq_coeff[j], m);

	mpf_set_d(m, pmod->sderr[j]);
	mpf_add(lmod->sum_sderr[j], lmod->sum_sderr[j], m);
	mpf_mul(m, m, m);
	mpf_add(lmod->ssq_sderr[j], lmod->ssq_sderr[j], m);
	if (!na(lmod->cbak[j]) && realdiff(pmod->coeff[j], lmod->cbak[j])) {
	    lmod->cdiff[j] = 1;
	}
	if (!na(lmod->sbak[j]) && realdiff(pmod->sderr[j], lmod->sbak[j])) {
	    lmod->sdiff[j] = 1;
	}
	lmod->cbak[j] = pmod->coeff[j];
	lmod->sbak[j] = pmod->sderr[j];
    }

    mpf_clear(m);

    lmod->n += 1;

#if LOOP_DEBUG
    fprintf(stderr, "loop_model_update: returning %d\n", err);
#endif

    return err;
}

/* Update the LOOP_PRINT struct @lprn using the current values of the
   specified variables. If this is the first use we need to do some
   allocation first.
*/

static int loop_print_update (LOOP_PRINT *lprn, const char *names) 
{
    mpf_t m;
    double x;
    int i, err = 0;

    if (lprn->names == NULL) {
	/* not started yet */
	err = loop_print_start(lprn, names);
	if (err) {
	    return err;
	}
    }

    mpf_init(m);
    
    for (i=0; i<lprn->nvars; i++) {
	if (lprn->na[i]) {
	    continue;
	}
	x = gretl_scalar_get_value(lprn->names[i], &err);
	if (err) {
	    break;
	}
	if (na(x)) {
	    lprn->na[i] = 1;
	    continue;
	}
	mpf_set_d(m, x); 
	mpf_add(lprn->sum[i], lprn->sum[i], m);
	mpf_mul(m, m, m);
	mpf_add(lprn->ssq[i], lprn->ssq[i], m);
	if (!na(lprn->xbak[i]) && realdiff(x, lprn->xbak[i])) {
	    lprn->diff[i] = 1;
	}
	lprn->xbak[i] = x;
    }

    mpf_clear(m);

    lprn->n += 1;

    return err;
}

static int add_more_loop_commands (LOOPSET *loop)
{
    int nb = 1 + (loop->n_cmds + 1) / LOOP_BLOCK;
    int totcmds = nb * LOOP_BLOCK;
    loop_command *cmds;

    /* in case we ran out of space */
    cmds = realloc(loop->cmds, totcmds * sizeof *cmds); 
    
    if (cmds == NULL) {
	return E_ALLOC;
    }

    loop->cmds = cmds;
    loop_cmds_init(loop, loop->n_cmds, totcmds);

    return 0;
}  

static int real_append_line (ExecState *s, LOOPSET *loop)
{
    const char *flagstr = NULL;
    int n = loop->n_cmds;
    int len, err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "real_append_line: s->line = '%s'\n", s->line);
#endif

    if ((n + 1) % LOOP_BLOCK == 0) {
	if (add_more_loop_commands(loop)) {
	    return E_ALLOC;
	}
    }

    len = strlen(s->line);

    if (s->cmd->opt) {
	flagstr = print_flags(s->cmd->opt, s->cmd->ci);
	len += strlen(flagstr);
	if (len >= MAXLINE) {
	    gretl_errmsg_set("loop: line is too long");
	    err = 1;
	}
    }

    if (!err) {
	if (flagstr != NULL) {
	    loop->cmds[n].line = malloc(len + 1);
	} else {
	    loop->cmds[n].line = gretl_strdup(s->line);
	}
	if (loop->cmds[n].line == NULL) {
	    err = E_ALLOC;
	} else {
	    if (flagstr != NULL) {
		sprintf(loop->cmds[n].line, "%s%s", s->line, flagstr);
	    }
	    compress_spaces(loop->cmds[n].line);
	}
    }

    if (!err) {
	if (s->cmd->ci == PRINT && (!loop_is_progressive(loop) || strchr(s->line, '"'))) {
	    /* printing a literal string, not a variable's value */
	    loop->cmds[n].flags |= LOOP_CMD_LIT;
	} 
	loop->cmds[n].ci = s->cmd->ci;
	loop->n_cmds += 1;
    }

#if LOOP_DEBUG
    fprintf(stderr, "loop %p: n_cmds=%d, line[%d]='%s', ci=%d\n",
	    (void *) loop, loop->n_cmds, n, loop->cmds[n].line,
	    loop->cmds[n].ci);
#endif

    return err;
}  

static void destroy_loop_stack (void)
{
    LOOPSET *loop = currloop;

    /* find the origin of the stack */
    while (loop->parent != NULL) {
	loop = loop->parent;
    }

    /* and destroy recursively */
    gretl_loop_destroy(loop);

    compile_level = 0;
}

/**
 * gretl_loop_append_line:
 * @s: program execution state.
 * @dset: dataset struct.
 *
 * Add command line to accumulated loop buffer.  This may be
 * called "starting from cold", in which case the "line"
 * member of @s will have been parsed and any options
 * extracted to s->cmd->opt.  But it may also be called in
 * the process of ongoing loop compilation, and in that
 * case option flags will not have been processed already.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_loop_append_line (ExecState *s, DATASET *dset)
{
    LOOPSET *loop = currloop;
    LOOPSET *newloop = currloop;
    int err = 0;

    warnmsg(s->prn); /* catch "end loop" if present */
    gretl_error_clear();

#if LOOP_DEBUG > 1
    fprintf(stderr, "gretl_loop_append_line: currloop = %p, line = '%s'\n", 
	    (void *) loop, s->line);
#endif

    if (!ok_in_loop(s->cmd->ci)) {
	gretl_errmsg_set(_("Sorry, this command is not available in loop mode"));
	fprintf(stderr, "ci = %d (%s)\n", s->cmd->ci, s->line);
	destroy_loop_stack();
	return E_NOTIMP;
    }

    if (s->cmd->ci == LOOP) {
	/* starting from scratch */
	gretlopt opt = get_loop_opts(s, &err);
	int nested = 0;

	if (!err) {
	    newloop = start_new_loop(s->line, loop, dset, 
				     opt, &nested, &err);
#if LOOP_DEBUG
	    fprintf(stderr, "got LOOP: newloop at %p (err = %d)\n", 
		    (void *) newloop, err);
#endif
	    if (newloop == NULL) {
		return err;
	    } else {
		set_loop_opts(newloop, opt);
		compile_level++;
		if (!nested) {
		    currloop = newloop;
		    return 0; /* done */
		}
	    }
	}
    } else if (s->cmd->ci == ENDLOOP) {
	/* got to the end */
	compile_level--;
#if LOOP_DEBUG
	fprintf(stderr, "got ENDLOOP, compile_level now %d\n",
		compile_level);
#endif
	if (compile_level == 0) {
	    /* set flag to run the loop */
	    loop_execute = 1;
	} else {
	    /* back up a level */
	    newloop = loop->parent;
	}
    } 

    if (!err && loop != NULL && s->cmd->ci != ENDLOOP) {
	err = real_append_line(s, loop);
    }

    if (err) {
	if (loop != NULL) {
	    gretl_loop_destroy(loop);
	    compile_level = 0;
	}
    } else {
	currloop = newloop;
    }

    return err;
}

static void print_loop_coeff (const DATASET *dset, 
			      const LOOP_MODEL *lmod, 
			      int i, PRN *prn)
{
    char pname[VNAMELEN];
    char tmp[NAMETRUNC];
    mpf_t c1, c2, m, sd1, sd2;
    unsigned long ln = lmod->n;

    mpf_init(c1);
    mpf_init(c2);
    mpf_init(m);
    mpf_init(sd1);
    mpf_init(sd2);

    mpf_div_ui(c1, lmod->sum_coeff[i], ln);
    if (lmod->cdiff[i] == 0) {
	mpf_set_d(sd1, 0.0);
    } else {
	mpf_mul(m, c1, c1);
	mpf_mul_ui(m, m, ln);
	mpf_sub(m, lmod->ssq_coeff[i], m);
	mpf_div_ui(sd1, m, ln);
	if (mpf_cmp_d(sd1, 0.0) > 0) {
	    mpf_sqrt(sd1, sd1);
	} else {
	    mpf_set_d(sd1, 0.0);
	}
    }

    mpf_div_ui(c2, lmod->sum_sderr[i], ln);
    if (lmod->sdiff[i] == 0) {
	mpf_set_d(sd2, 0.0);
    } else {
	mpf_mul(m, c2, c2);
	mpf_mul_ui(m, m, ln);
	mpf_sub(m, lmod->ssq_sderr[i], m);
	mpf_div_ui(sd2, m, ln);
	if (mpf_cmp_d(sd2, 0.0) > 0) {
	    mpf_sqrt(sd2, sd2);
	} else {
	    mpf_set_d(sd2, 0.0);
	}
    }

    gretl_model_get_param_name(lmod->model0, dset, i, pname);
    maybe_trim_varname(tmp, pname);
    pprintf(prn, "%*s", 15, tmp); /* FIXME length */
    pprintf(prn, "%#14g %#14g %#14g %#14g\n", mpf_get_d(c1), mpf_get_d(sd1), 
	    mpf_get_d(c2), mpf_get_d(sd2));

    mpf_clear(c1);
    mpf_clear(c2);
    mpf_clear(m);
    mpf_clear(sd1);
    mpf_clear(sd2);
}

static void loop_model_print (LOOP_MODEL *lmod, const DATASET *dset, 
			      PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    int i;

    ntodate(startdate, lmod->model0->t1, dset);
    ntodate(enddate, lmod->model0->t2, dset);

    pputc(prn, '\n');
    pprintf(prn, _("%s estimates using the %d observations %s-%s\n"),
	    _(estimator_string(lmod->model0, prn)), lmod->model0->nobs, 
	    startdate, enddate);
    print_model_vcv_info(lmod->model0, dset, prn);
    pprintf(prn, _("Statistics for %d repetitions\n"), lmod->n); 
    pprintf(prn, _("Dependent variable: %s\n\n"), 
	    gretl_model_get_depvar_name(lmod->model0, dset));

    pputs(prn, _("                     mean of      std. dev. of     mean of"
		 "     std. dev. of\n"
		 "                    estimated      estimated"
		 "      estimated      estimated\n"
		 "      Variable     coefficients   coefficients   std. errors"
		 "    std. errors\n\n"));

    for (i=0; i<lmod->model0->ncoeff; i++) {
	print_loop_coeff(dset, lmod, i, prn);
    }

    pputc(prn, '\n');
}

static void loop_print_print (LOOP_PRINT *lprn, PRN *prn)
{
    bigval mean, m, sd;
    int len, maxlen = 7;
    int i, n;
    const char *s;

    if (lprn == NULL) {
	return;
    }

    n = lprn->n;

    mpf_init(mean);
    mpf_init(m);
    mpf_init(sd);

    for (i=0; i<lprn->nvars; i++) {
	len = strlen(lprn->names[i]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    pprintf(prn, _("Statistics for %d repetitions\n"), n); 
    pputc(prn, '\n');
    bufspace(maxlen + 1, prn);

    len = get_utf_width(_("mean"), 14);
    pprintf(prn, "%*s ", len, _("mean"));

    len = get_utf_width(_("std. dev"), 14);
    pprintf(prn, "%*s\n", len, _("std. dev"));
    
    for (i=0; i<lprn->nvars; i++) {
	s = lprn->names[i];
	if (lprn->na[i]) {
	    pprintf(prn, "%*s", maxlen + 1, s);
	    pprintf(prn, "%14s %14s\n", "NA   ", "NA   ");
	    continue;
	}
	mpf_div_ui(mean, lprn->sum[i], (unsigned long) n);
	if (lprn->diff[i] == 0) {
	    mpf_set_d(sd, 0.0);
	} else {
	    mpf_mul(m, mean, mean);
	    mpf_mul_ui(m, m, (unsigned long) n);
	    mpf_sub(sd, lprn->ssq[i], m);
	    mpf_div_ui(sd, sd, (unsigned long) n);
	    if (mpf_cmp_d(sd, 0.0) > 0) {
		mpf_sqrt(sd, sd);
	    } else {
		mpf_set_d(sd, 0.0);
	    }
	}
	pprintf(prn, "%*s", maxlen + 1, s);
	pprintf(prn, "%#14g %#14g\n", mpf_get_d(mean), mpf_get_d(sd));
    }

    mpf_clear(mean);
    mpf_clear(m);
    mpf_clear(sd);
 
    pputc(prn, '\n');
}

static int loop_store_save (LOOP_STORE *lstore, PRN *prn)
{
    int *list;
    int err = 0;

    list = gretl_consecutive_list_new(1, lstore->dset->v - 1);
    if (list == NULL) {
	return E_ALLOC;
    }

    lstore->dset->t2 = lstore->n - 1;
    pprintf(prn, _("store: using filename %s\n"), lstore->fname);
    err = write_data(lstore->fname, list, lstore->dset, lstore->opt, prn);

    if (err) {
	pprintf(prn, _("write of data file failed\n"));
    }

    free(list);

    return err;
}

#define loop_literal(l,i) (l->cmds[i].flags & LOOP_CMD_LIT)

/**
 * print_loop_results:
 * @loop: pointer to loop struct.
 * @dset: data information struct.
 * @prn: gretl printing struct.
 *
 * Print out the results after completion of the loop @loop.
 */

static void print_loop_results (LOOPSET *loop, const DATASET *dset, 
				PRN *prn)
{
    int iters = loop->iter;
    int i, j = 0, k = 0;

    if (loop->type != COUNT_LOOP && !(loop_is_quiet(loop))) {
	pprintf(prn, _("\nNumber of iterations: %d\n\n"), iters);
    }

    for (i=0; i<loop->n_cmds; i++) {
	gretlopt opt = OPT_NONE;
	int ci = loop->cmds[i].ci;

#if LOOP_DEBUG
	fprintf(stderr, "print_loop_results: loop command %d: %s\n", 
		i, loop->cmds[i].line);
#endif

	if (plain_model_ci(ci)) {
	    char linecpy[MAXLINE];

	    strcpy(linecpy, loop->cmds[i].line);
	    opt = get_gretl_options(linecpy, NULL);
	}	    

	if (!loop_is_progressive(loop) && ci == OLS) {
	    if (model_print_deferred(opt)) {
		MODEL *pmod = loop->models[j++];

		set_model_id(pmod);
		printmodel(pmod, dset, opt, prn);
	    }	    
	}

	if (loop_is_progressive(loop)) {
	    if (plain_model_ci(ci) && !(opt & OPT_Q)) {
		loop_model_print(&loop->lmodels[j], dset, prn);
		loop_model_zero(&loop->lmodels[j], 1);
		j++;
	    } else if (ci == PRINT && !loop_literal(loop, i)) {
		loop_print_print(&loop->prns[k], prn);
		loop_print_zero(&loop->prns[k], 1);
		k++;
	    } else if (ci == STORE) {
		loop_store_save(&loop->store, prn);
	    }
	}
    }
}

static int substitute_dollar_targ (char *str, int maxlen, 
				   const LOOPSET *loop,
				   const DATASET *dset,
				   int *subst)
{
    char insert[32], targ[VNAMELEN + 3] = {0};
    char *p, *ins, *q, *s;
    int targlen, inslen, idx = 0;
    int incr, cumlen = 0;
    int err = 0;

#if SUBST_DEBUG
    fprintf(stderr, "subst_dollar_targ:\n original: '%s'\n", str);
#endif

    /* construct the target for substitution */

    if (loop->type == FOR_LOOP) {
	if (!gretl_is_scalar(loop->init.vname)) {
	    /* nothing to substitute */
	    return 0;
	} 
	sprintf(targ, "$%s", loop->init.vname);
	targlen = strlen(targ);
    } else if (indexed_loop(loop)) {
	sprintf(targ, "$%s", loop->idxname);
	targlen = strlen(targ);
	idx = loop->init.val + loop->iter;
    } else {
	/* shouldn't be here! */
	return 1;
    }

#if SUBST_DEBUG
    fprintf(stderr, " target = '%s', idx = %d\n", targ, idx);
#endif

    if (strstr(str, targ) == NULL) {
	/* nothing to be done */
	return 0;
    }

    ins = insert;

    /* prepare the substitute string */

    if (loop->type == FOR_LOOP) {
	double x = gretl_scalar_get_value(loop->init.vname, NULL);

	if (na(x)) {
	    strcpy(insert, "NA");
	} else {
	    sprintf(insert, "%g", x);
	}
    } else if (loop->type == INDEX_LOOP) {
	sprintf(insert, "%d", idx);
    } else if (loop->type == DATED_LOOP) {
	/* note: ntodate is 0-based */
	ntodate(insert, idx - 1, dset);
    } else if (loop->type == EACH_LOOP) {
	ins = loop->eachstrs[idx - 1];
    } 

    inslen = strlen(ins);
    incr = inslen - targlen;
    if (incr > 0) {
	/* substitution will lengthen the string */
	cumlen = strlen(str);
    }

    q = malloc(strlen(strstr(str, targ)));
    if (q == NULL) {
	err = E_ALLOC;
    }

    /* crawl along str, replacing targ with ins */

    s = str;
    while ((p = strstr(s, targ)) != NULL && !err) {
	if (is_gretl_accessor(p)) {
	    s++;
	    continue;
	}
	if (incr > 0) {
	    cumlen += incr;
	    if (cumlen >= maxlen) {
		/* substitution would cause overflow */
		err = (maxlen == VNAMELEN)? E_UNKVAR : E_TOOLONG;
		break;
	    }
	}
	strcpy(q, p + targlen);
	strcpy(p, ins);
	strcpy(p + inslen, q);
	if (subst != NULL) {
	    *subst = 1;
	}
	s++; /* += strlen(ins)? */
    }

    free(q);

#if SUBST_DEBUG
    fprintf(stderr, " after: '%s'\n", str);
#endif

    return err;
}

static int extend_loop_dataset (LOOP_STORE *lstore)
{
    double *x;
    int oldn = lstore->dset->n;
    int n = oldn + DEFAULT_NOBS;
    int i, t;

    for (i=0; i<lstore->dset->v; i++) {
	x = realloc(lstore->dset->Z[i], n * sizeof *x);
	if (x == NULL) {
	    return E_ALLOC;
	}
	lstore->dset->Z[i] = x;
	for (t=oldn; t<n; t++) {
	    lstore->dset->Z[i][t] = (i == 0)? 1.0 : NADBL;
	}	    
    }
    
    lstore->dset->n = n;
    lstore->dset->t2 = n - 1;

    ntodate(lstore->dset->endobs, n - 1, lstore->dset);

    return 0;
}

static void progressive_loop_zero (LOOPSET *loop)
{
    int i;

    for (i=0; i<loop->n_loop_models; i++) {
	loop_model_free(&loop->lmodels[i]);
    }

    loop->lmodels = NULL;
    loop->n_loop_models = 0;

    for (i=0; i<loop->n_prints; i++) { 
	loop_print_free(&loop->prns[i]);
    }

    loop->prns = NULL;
    loop->n_prints = 0;

    loop_store_free(&loop->store);
}

static int top_of_loop (LOOPSET *loop, DATASET *dset)
{
    int err = 0;

    loop->iter = 0;

    if (loop->listname[0] != '\0') {
	err = loop_list_refresh(loop, dset);
    } else if (loop->type == INDEX_LOOP) {
	loop->init.val = controller_get_val(&loop->init, loop, dset, &err);
    } else if (loop->type == FOR_LOOP) {
	forloop_init(loop, dset, &err);
    }

    if (!err && (loop->type == COUNT_LOOP || indexed_loop(loop))) {
	loop->final.val = controller_get_val(&loop->final, loop, dset, &err);
	if (na(loop->init.val) || na(loop->final.val)) {
	    gretl_errmsg_set(_("error evaluating loop condition"));
	    fprintf(stderr, "loop: got NA for init and/or final value\n");
	    err = E_DATA;
	} else {
	    loop->itermax = loop->final.val - loop->init.val + 1;
#if LOOP_DEBUG
	    fprintf(stderr, "*** itermax = %g - %g + 1 = %d\n",
		    loop->final.val, loop->init.val, loop->itermax);
#endif
	}
    }

    if (!err) {
	if (indexed_loop(loop)) {
	    loop->idxval = loop->init.val;
	    gretl_scalar_set_value_authorized(loop->idxname, loop->idxval);
	}

	/* initialization, in case this loop is being run more than
	   once (i.e. it's embedded in an outer loop) */

	if (loop_is_progressive(loop)) {
	    progressive_loop_zero(loop);
	} else {
	    free(loop->models);
	    loop->models = NULL;
	    loop->n_models = 0;
	}
    }

    return err;
}

static void 
print_loop_progress (const LOOPSET *loop, const DATASET *dset,
		     PRN *prn)
{
    int i = loop->init.val + loop->iter;

    if (loop->type == INDEX_LOOP) {
	pprintf(prn, "loop: %s = %d\n\n", loop->idxname, i);
    } else if (loop->type == DATED_LOOP) {
	char obs[OBSLEN];

	ntodate(obs, i - 1, dset);
	pprintf(prn, "loop: %s = %s\n\n", loop->idxname, obs);
    }
}

static const LOOPSET *
subst_loop_in_parentage (const LOOPSET *loop)
{
    while ((loop = loop->parent) != NULL) {
	if (indexed_loop(loop) || loop->type == FOR_LOOP) break;
    }

    return loop;
}

static int 
make_dollar_substitutions (char *str, int maxlen,
			   const LOOPSET *loop,
			   const DATASET *dset,
			   int *subst,
			   gretlopt opt)
{
    int err = 0;

    /* if (opt & OPT_T) we're just processing a variable name, at the top
       of a loop, so we can skip to the "parentage" bit 
    */

    if (!(opt & OPT_T) && (indexed_loop(loop) || loop->type == FOR_LOOP)) {
	err = substitute_dollar_targ(str, maxlen, loop, dset, subst);
    }

    while (!err && (loop = subst_loop_in_parentage(loop)) != NULL) {
	err = substitute_dollar_targ(str, maxlen, loop, dset, subst);
    }

    return err;
}

int scalar_is_read_only_index (const char *name)
{
    const LOOPSET *loop = currloop;

    while (loop != NULL) {
	if (indexed_loop(loop) && !strcmp(name, loop->idxname)) {
	    return 1;
	}
	loop = loop->parent;
    }

    return 0;
}

static LOOPSET *get_child_loop_by_line (LOOPSET *loop, int lno)
{
    int i;

    for (i=0; i<loop->n_children; i++) {
	if (loop->children[i]->parent_line == lno) {
	    return loop->children[i];
	}
    }

    return NULL;
}

static int add_loop_genr (LOOPSET *loop, int lno, 
			  const char *line, 
			  DATASET *dset)
{
    int err = 0;

    loop->cmds[lno].genr = genr_compile(line, dset, OPT_NONE, &err);
    if (!err) {
	loop->cmds[lno].flags |= LOOP_CMD_GENR;
    }

    return err;
}

static int loop_print_save_model (MODEL *pmod, DATASET *dset,
				  PRN *prn, ExecState *s)
{
    int err = pmod->errcode;
    
    if (!err) {
	set_gretl_errno(0);
	if (!(s->cmd->opt & OPT_Q)) {
	    printmodel(pmod, dset, s->cmd->opt, prn);
	}
	attach_subsample_to_model(pmod, dset);
	s->pmod = maybe_stack_model(pmod, s->cmd, prn, &err);
	if (!err && s->callback != NULL && *s->cmd->savename != '\0' &&
	    gretl_in_gui_mode()) {
	    s->callback(s, s->pmod, GRETL_OBJ_EQN);
	}
    } 

    return err;
}

/* get the next command for a loop by pulling a line off the
   stack of loop commands.
*/

static int loop_next_command (char *targ, LOOPSET *loop, int *pj)
{
    int ret = 1, j = *pj + 1;

    if (j < loop->n_cmds) {
	strcpy(targ, loop->cmds[j].line);
	*pj = j;
    } else {
	ret = 0;
    }

    return ret;
}

#define genr_compiled(l,j) (l->cmds[j].flags & LOOP_CMD_GENR)
#define loop_cmd_nodol(l,j) (l->cmds[j].flags & LOOP_CMD_NODOL)
#define loop_cmd_nosub(l,j) (l->cmds[j].flags & LOOP_CMD_NOSUB)
#define loop_cmd_noopt(l,j) (l->cmds[j].flags & LOOP_CMD_NOOPT)
#define loop_cmd_catch(l,j) (l->cmds[j].flags & LOOP_CMD_CATCH)
#define conditional_compiled(l,j) (l->cmds[j].flags & LOOP_CMD_COND)

static int loop_process_error (LOOPSET *loop, int j, int err, PRN *prn)
{
#if LOOP_DEBUG
    fprintf(stderr, "loop_process_error: j=%d, catch=%d\n",
	    j, loop_cmd_catch(loop, j));
#endif
    if (loop_cmd_catch(loop, j)) {
	err = 0;
    } else if (!libset_get_bool(HALT_ON_ERR)) {
	errmsg(err, prn);
	err = 0;
    }

    return err;
}

/* Based on the stored flags in the loop-line record, set
   or unset some flags for the command parser: this can 
   reduce the amount of work the parser has to do on each
   iteration of a loop.
*/

static inline void loop_info_to_cmd (LOOPSET *loop, int j,
				     CMD *cmd)
{
    if (loop_is_progressive(loop)) {
	cmd->flags |= CMD_PROG;
    } else {
	cmd->flags &= ~CMD_PROG;
    }

    if (loop_cmd_nosub(loop, j)) {
	/* tell parser not to try for @-substitution */
	cmd->flags |= CMD_NOSUB;
    } else {
	cmd->flags &= ~CMD_NOSUB;
    }
    
    if (loop_cmd_noopt(loop, j)) {
	/* tell parser not to try for option flags */
	cmd->flags |= CMD_NOOPT;
	/* since the cmd->opt is (or can be) persistent, we
	   should zero it here */
	cmd->opt = OPT_NONE;
    } else {
	cmd->flags &= ~CMD_NOOPT;
    }

    if (loop_cmd_catch(loop, j)) {
	cmd->flags |= CMD_CATCH;
    } else if (!cmd->context) {
	cmd->flags &= ~CMD_CATCH;
    }
}

/* Based on the parsed info in @cmd, maybe modify some flags in
   the current loop-line record.
*/

static inline void cmd_info_to_loop (LOOPSET *loop, int j,
				     CMD *cmd, int *subst)
{
    if (!loop_cmd_nosub(loop, j)) {
	/* this loop line has not already been marked as
	   free of @-substitution
	*/
	if (cmd_subst(cmd)) {
	    *subst = 1;
	} else {
	    /* record: no @-substitution in this line */
	    loop->cmds[j].flags |= LOOP_CMD_NOSUB;
	}
    }

    if (loop_cmd_nosub(loop, j)) {
	if (cmd->opt == OPT_NONE) {
	    /* record: no options are present on this line */
	    loop->cmds[j].flags |= LOOP_CMD_NOOPT;
	} 
    }

    if (cmd->flags & CMD_CATCH) {
	loop->cmds[j].flags |= LOOP_CMD_CATCH;
    }
}

static int loop_delete_object (CMD *cmd, PRN *prn)
{
    int err = 0;

    if (cmd->list != NULL && cmd->list[0] > 0) {
	pputs(prn, _("You cannot delete series in this context\n"));
	err = 1;
    } else if (gretl_is_scalar(cmd->param)) {
	pputs(prn, _("You cannot delete scalars in this context\n"));
	err = 1;
    } else {
	err = gretl_delete_var_by_name(cmd->param, prn);
    }

    return err;
}

static int loop_report_error (LOOPSET *loop, int err, 
			      const char *errline,
			      ExecState *state,
			      PRN *prn)
{
    int fd = gretl_function_depth();

    if (err) {
	if (fd == 0) {
	    errmsg(err, prn);
	    if (*errline != '\0') {
		pprintf(prn, ">> %s\n", errline);
	    }
	}
    } else if (loop->err) {
	if (fd == 0) {
	    errmsg(loop->err, prn);
	}
	err = loop->err;
    }

    if (fd > 0 && err && *errline != '\0') {
	strcpy(state->line, errline);
    }

    return err;
}

static int conditional_line (LOOPSET *loop, int j)
{
    return loop->cmds[j].ci == IF || loop->cmds[j].ci == ELIF;
}

#define COMPILE_IF 1 /* should be OK now? */

#if COMPILE_IF 

static int do_compile_conditional (LOOPSET *loop, int j)
{
    if ((loop->cmds[j].ci == IF || loop->cmds[j].ci == ELIF) &&
	loop_cmd_nodol(loop, j) && loop_cmd_nosub(loop, j)) {
	return 1;
    } else {
	return 0;
    }
}

#endif	

static int block_model (CMD *cmd)
{
    return cmd->ci == END && 
	(!strcmp(cmd->param, "mle") || 
	 !strcmp(cmd->param, "nls") ||
	 !strcmp(cmd->param, "gmm"));
}

#define not_ok_in_progloop(c) (NEEDS_MODEL_CHECK(c) || \
			       c == NLS ||  \
			       c == MLE ||  \
			       c == GMM)

int gretl_loop_exec (ExecState *s, DATASET *dset) 
{
    LOOPSET *loop = currloop;
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    MODEL *pmod;
    LOOP_MODEL *lmod;
    LOOP_PRINT *lprn;
    char errline[MAXLINE];
    int indent0;
    int show_activity = 0;
    int j, err = 0;

    /* for the benefit of the caller: register the fact that execution
       of this loop is now under way */
    loop_execute = 0;

    if (loop == NULL) {
	pputs(prn, "Got a NULL loop\n");
	set_loop_off();
	return 1;
    }

    indent0 = gretl_if_state_record();

    set_loop_on(loop_is_quiet(loop), loop_is_progressive(loop));

#if LOOP_DEBUG
    fprintf(stderr, "loop_exec: loop = %p\n", (void *) loop);
#endif

    err = top_of_loop(loop, dset);
    if (err) {
	*errline = '\0';
    }

    show_activity = show_activity_func_installed();
    
    while (!err && loop_condition(loop, dset, &err)) {
#if LOOP_DEBUG
	fprintf(stderr, "*** top of loop: iter = %d\n", loop->iter);
#endif
	j = -1;

	pmod = NULL;
	lmod = NULL;
	lprn = NULL;

	if (gretl_echo_on() && indexed_loop(loop) && !loop_is_quiet(loop)) {
	    print_loop_progress(loop, dset, prn);
	}

	while (!err && loop_next_command(line, loop, &j)) {
	    int subst = 0;

#if LOOP_DEBUG
	    fprintf(stderr, " j=%d, line='%s', compiled = %d\n", j, line,
	            genr_compiled(loop, j));
#endif
	    strcpy(errline, line);

	    if (!gretl_if_state_false() && genr_compiled(loop, j)) {
		/* If the current line already has "compiled genr" status,
		   we should be able to skip several steps that are
		   potentially quite time-consuming.
		*/
		if (gretl_echo_on() && !loop_is_quiet(loop)) {
		    pprintf(prn, "? %s\n", line);
		}
		err = execute_genr(loop->cmds[j].genr, dset, prn);
		if (err) {
		    err = loop_process_error(loop, j, err, prn);
		}
		if (err) {
		    break;
		} else {
		    continue;
		}
	    }

	    if (!loop_cmd_nodol(loop, j)) {
		if (strchr(line, '$')) {
		    /* handle loop-specific $-string substitution */
		    err = make_dollar_substitutions(line, MAXLINE, loop, 
						    dset, &subst, OPT_NONE);
		    if (err) {
			break;
		    } else if (!subst) {
			loop->cmds[j].flags |= LOOP_CMD_NODOL;
		    }
		} else {
		    loop->cmds[j].flags |= LOOP_CMD_NODOL;
		}
	    }

	    /* transcribe loop -> cmd */
	    loop_info_to_cmd(loop, j, cmd);

	    /* call the full command parser, with special treatment
	       for "if" or "elif" conditions that may be already
	       compiled, or that should now be compiled
	    */
#if COMPILE_IF
	    if (conditional_compiled(loop, j)) {
		err = parse_command_line(line, cmd, dset, &loop->cmds[j].genr);
	    } else if (do_compile_conditional(loop, j)) {
		GENERATOR *ifgen = NULL;

		err = parse_command_line(line, cmd, dset, &ifgen);
		if (ifgen != NULL) {
		    loop->cmds[j].genr = ifgen;
		    loop->cmds[j].flags |= LOOP_CMD_COND;
		}
	    } else {
		err = parse_command_line(line, cmd, dset, NULL);
	    }
#else
	    err = parse_command_line(line, cmd, dset, NULL);
#endif

#if LOOP_DEBUG
	    fprintf(stderr, "    after: '%s'\n", line);
	    fprintf(stderr, "    cmd->savename = '%s'\n", cmd->savename);
	    fprintf(stderr, "    err from parse_command_line: %d\n", err);
#endif

	    if (err) {
		err = loop_process_error(loop, j, err, prn);
		if (err) {
		    break;
		} else {
		    continue;
		}
	    } else if (cmd->ci < 0) {
		if (conditional_line(loop, j)) {
		    cmd_info_to_loop(loop, j, cmd, &subst);
		}
		continue;
	    } else {
		gretl_exec_state_transcribe_flags(s, cmd);
		cmd_info_to_loop(loop, j, cmd, &subst);
	    }

	    if (gretl_echo_on()) {
		if (s->cmd->ci == ENDLOOP) {
		    if (indexed_loop(loop)) {
			pputc(prn, '\n');
		    }
		} else if (!loop_is_quiet(loop)) {
		    echo_command(cmd, dset, line, prn);
		}
	    }

	    /* now branch based on the command index: some commands
	       require special treatment in loop context
	    */

	    if (cmd->ci == LOOP) {
		currloop = get_child_loop_by_line(loop, j);
		if (currloop == NULL) {
		    currloop = loop;
		    fprintf(stderr, "Got a LOOP command, don't know what to do!\n");
		    err = 1;
		} else {
		    err = gretl_loop_exec(s, dset);
		}
	    } else if (cmd->ci == BREAK) {
		loop->brk = 1;
		break;
	    } else if (cmd->ci == ENDLOOP) {
		; /* implicit break */
	    } else if (plain_model_ci(cmd->ci)) {
		/* model may need special handling */
		if (loop_is_progressive(loop) && !(cmd->opt & OPT_Q)) {
		    lmod = get_loop_model_by_line(loop, j, &err);
		} else if (model_print_deferred(cmd->opt)) {
		    pmod = get_model_record_by_line(loop, j, &err);
		}
		/* estimate the model called for */
		if (!err) {
		    err = gretl_cmd_exec(s, dset);
		}
		if (!err) {
		    int moderr = check_gretl_errno();

		    if (moderr) {
			if (loop_is_progressive(loop) || model_print_deferred(cmd->opt)) {
			    err = moderr;
			} else {
			    errmsg(moderr, prn);
			}
		    } else if (loop_is_progressive(loop) && !(cmd->opt & OPT_Q)) {
			err = loop_model_update(lmod, s->model);
			set_as_last_model(s->model, GRETL_OBJ_EQN);
		    } else if (model_print_deferred(cmd->opt)) {
			swap_models(s->model, pmod);
			pmod->ID = j + 1;
			set_as_last_model(pmod, GRETL_OBJ_EQN);
			model_count_minus();
		    } else {
			loop_print_save_model(s->model, dset, prn, s);
		    }
		}
	    } else if (cmd->ci == PRINT && *cmd->param == '\0' &&
		       loop_is_progressive(loop)) {
		lprn = get_loop_print_by_line(loop, j, &err);
		if (!err) {
		    err = loop_print_update(lprn, cmd->parm2);
		}
	    } else if (cmd->ci == STORE && loop_is_progressive(loop)) {
		err = loop_store_update(loop, j, cmd->parm2, cmd->param,
					cmd->opt);
	    } else if (loop_is_progressive(loop) && not_ok_in_progloop(cmd->ci)) {
		gretl_errmsg_sprintf(_("%s: not implemented in 'progressive' loops"),
				     cmd->word);
		err = 1;
	    } else if (cmd->ci == GENR) {
		if (subst || (cmd->opt & OPT_O)) {
		    /* can't use a "compiled" genr if string substitution
		       has been done, since the genr expression will not
		       be constant 
		    */
		    if (!loop_is_verbose(loop)) {
			cmd->opt |= OPT_Q;
		    }
		    err = generate(line, dset, cmd->opt, prn);
		} else {
		    err = add_loop_genr(loop, j, line, dset);
		    if (!err) {
			err = execute_genr(loop->cmds[j].genr, dset, prn);
		    }
		}
	    } else if (cmd->ci == DELEET && !(cmd->opt & (OPT_F | OPT_T))) {
		err = loop_delete_object(cmd, prn);
	    } else {
		err = gretl_cmd_exec(s, dset);
		if (!err && !check_gretl_errno() && block_model(cmd)) {
		    /* NLS, etc. */
		    loop_print_save_model(s->model, dset, prn, s);
		}
	    }

	    if (err && (cmd->flags & CMD_CATCH)) {
		set_gretl_errno(err);
		cmd->flags ^= CMD_CATCH;
		err = 0;
	    }
	} /* end execution of commands within loop */

	if (err) {
	    if (!libset_get_bool(HALT_ON_ERR)) {
		/* print error message but keep going */
		errmsg(err, prn);
		err = 0;
		/* should we clear the "if state" here? */
	    } else {
		gretl_if_state_clear();
	    }
	} 

	if (!err) {
	    if (loop->brk) {
		gretl_if_state_reset(indent0);
	    } else {
		err = gretl_if_state_check(indent0);
	    }
	}

	if (!err && !loop->brk) {
	    loop->iter += 1;
	    if (show_activity && (loop->iter % 10 == 0)) {
		show_activity_callback();
	    }
	}

    } /* end iterations of loop */

    cmd->flags &= ~CMD_NOSUB;
    cmd->flags &= ~CMD_NOOPT;

    if (loop->brk) {
	/* turn off break flag */
	loop->brk = 0;
    }

    if (err || loop->err) {
	err = loop_report_error(loop, err, errline, s, prn);
    }

    if (!err && loop->iter > 0) {
	print_loop_results(loop, dset, prn); 
    }

    if (loop->n_models > 0) {
	/* we need to update models[0] */
	GretlObjType type;
	void *ptr = get_last_model(&type);
	int i;

	if (type == GRETL_OBJ_EQN && s->model != ptr) {
	    swap_models(s->model, loop->models[loop->n_models - 1]);
	    set_as_last_model(s->model, GRETL_OBJ_EQN);
	}
	for (i=0; i<loop->n_models; i++) {
	    gretl_model_free(loop->models[i]);
	}
    }

    if (err && gretl_function_depth() > 0) {
	; /* leave 'line' alone */
    } else if (line != NULL) {
	*line = '\0';
    } 

    /* be sure to clear some loop-special parser flags */
    cmd->flags &= ~CMD_PROG;

    if (loop->parent == NULL) {
	/* reached top of stack: clean up */
	gretl_loop_destroy(loop);
	currloop = NULL;
	set_loop_off();
    }

    return process_command_error(cmd, err);
}
