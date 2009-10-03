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

/*  monte_carlo.c - loop procedures */

#include "libgretl.h" 
#include "monte_carlo.h"
#include "libset.h"
#include "compat.h"
#include "cmd_private.h"
#include "var.h"
#include "objstack.h"
#include "gretl_func.h"
#include "gretl_scalar.h"
#include "flow_control.h"

#include <time.h>
#include <unistd.h>

#define LOOP_DEBUG 0
#define SUBST_DEBUG 0

#if LOOP_DEBUG
# undef ENABLE_GMP
#endif

#if defined(ENABLE_GMP)
# include <gmp.h>
  typedef mpf_t bigval;
#elif defined(HAVE_LONG_DOUBLE)
  typedef long double bigval;
#else
  typedef double bigval;
#endif

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
    int linenum;   /* location: line number in loop */
    int n;         /* number of repetitions */
    int *list;     /* list of vars to print */
    bigval *sum;   /* running sum of values */
    bigval *ssq;   /* running sum of squares */
    double *xbak;  /* previous values */
    int *diff;     /* indicator for difference */
} LOOP_PRINT;  

/* below: used only in "progressive" loops */ 

typedef struct {
    int linenum;            /* location: line number in loop */
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
    int linenum;      /* location: line number in loop */ 
    int n;            /* number of observations */
    char *fname;      /* filename for output */
    gretlopt opt;     /* formatting option */
    double **Z;       /* data storage */
    DATAINFO *dinfo;  /* data info */
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
    int n_genrs;
    int n_children;

    /* subsidiary objects */
    char **lines;
    int *ci;
    char **eachstrs;
    GENERATOR **genrs;
    MODEL **models;
    int *model_lines;
    LOOP_MODEL *lmodels;
    LOOP_PRINT *prns;
    LOOP_STORE store;
    LOOPSET *parent;
    LOOPSET **children;
    int parent_line;
};

#define loop_is_progressive(l) (l->flags & LOOP_PROGRESSIVE)
#define loop_set_progressive(l) (l->flags |= LOOP_PROGRESSIVE)
#define loop_is_verbose(l) (l->flags & LOOP_VERBOSE)
#define loop_set_verbose(l) (l->flags |= LOOP_VERBOSE)
#define loop_is_quiet(l) (l->flags & LOOP_QUIET)
#define loop_set_quiet(l) (l->flags |= LOOP_QUIET)

#define is_list_loop(l) (l->listname[0] != '\0')

#define model_print_deferred(o) (o & OPT_F)

static void controller_init (controller *clr);
static int gretl_loop_prepare (LOOPSET *loop);
static void loop_model_free (LOOP_MODEL *lmod);
static void loop_print_free (LOOP_PRINT *lprn);
static void loop_store_free (LOOP_STORE *lstore);
static int extend_loop_dataset (LOOP_STORE *lstore);
static void controller_free (controller *clr);

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
   look up its current value (and modify the sign if wanted).
   Otherwise we should have got a numerical constant at setup, in
   which case we just return that value.
*/

static double controller_get_val (controller *clr)
{
    if (clr->vname[0] != '\0') {
	if (gretl_is_scalar(clr->vname)) {
	    clr->val = gretl_scalar_get_value(clr->vname) * clr->vsign;
	} else {
	    gretl_errmsg_sprintf(_("'%s': not a scalar"), clr->vname);
	} 
    } 

#if LOOP_DEBUG
    fprintf(stderr, "controller_get_val: vname='%s', returning %g\n", 
	    clr->vname, clr->val);
#endif

    return clr->val;
}

/* apply initialization in case of for-loop */

static void
forloop_init (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo, int *err)
{
    const char *expr = loop->init.expr;

    if (expr != NULL) {
	*err = generate(expr, pZ, pdinfo, OPT_Q, NULL);
	if (*err) {
	    gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
				 expr);
	}
    }
}

/* evaluate boolean condition in for-loop or while-loop */

static int 
loop_testval (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo, int *err)
{
    const char *expr = loop->test.expr;
    int ret = 1;

    if (expr != NULL) {
	double x = generate_scalar(expr, pZ, pdinfo, err);

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
loop_delta (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo, int *err)
{
    const char *expr = loop->delta.expr;

    if (expr != NULL) {
	*err = generate(expr, pZ, pdinfo, OPT_Q, NULL);
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
	c == DATA ||
	c == DELEET ||
	c == EQNPRINT ||
	c == FOREIGN ||
	c == FUNC ||
	c == HURST ||
	c == INCLUDE ||
	c == LEVERAGE ||
	c == MODELTAB ||
	c == NULLDATA ||
	c == OPEN ||
	c == RMPLOT ||
	c == RUN ||
	c == SCATTERS ||
	c == SETMISS ||
	c == SETOBS ||
	c == TABPRINT ||
	c == VIF ||
	c == XCORRGM) {
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
    child->parent_line = loop->n_cmds + 1;
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
    lstore->linenum = -1;
    lstore->n = 0;
    lstore->fname = NULL;
    lstore->opt = OPT_NONE;
    lstore->Z = NULL;
    lstore->dinfo = NULL;
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

    loop->n_genrs = 0;
    loop->genrs = NULL;

    loop->lines = NULL;
    loop->ci = NULL;
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

    for (i=0; i<loop->n_children; i++) {
	gretl_loop_destroy(loop->children[i]);
    }

    controller_free(&loop->init);
    controller_free(&loop->test);
    controller_free(&loop->delta);
    controller_free(&loop->final);

    if (loop->lines != NULL) {
	for (i=0; i<loop->n_cmds; i++) {
	    free(loop->lines[i]);
	}
	free(loop->lines);
    }

    if (loop->ci != NULL) { 
	free(loop->ci);
    }

    if (loop->model_lines != NULL) { 
	free(loop->model_lines);
    }    

    if (loop->eachstrs != NULL) {
	free_strings_array(loop->eachstrs, loop->itermax);
    }

    if (loop->models != NULL) {
	free(loop->models);
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

    for (i=0; i<loop->n_genrs; i++) {
	destroy_genr(loop->genrs[i]);
    }
    free(loop->genrs);    

    if (loop->children != NULL) {
	free(loop->children);
    }

    if (loop->flags & LOOP_DELVAR) {
	gretl_scalar_delete(loop->idxname, NULL);
    }

    free(loop);
}

static int parse_as_while_loop (LOOPSET *loop,
				double ***pZ, DATAINFO *pdinfo,
				const char *s)
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
				  const DATAINFO *pdinfo)
{
    int err = 0;

    if (gretl_is_scalar(vname)) {
	strcpy(loop->idxname, vname);
	gretl_scalar_set_value(vname, loop->init.val);
    } else {
	err = gretl_scalar_add_with_check(vname, loop->init.val, pdinfo);
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
			    const double **Z, const DATAINFO *pdinfo)
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
	    clr->val = (int) gretl_scalar_get_value(s);
	} else if ((v = current_series_index(pdinfo, s)) >= 0) {
	    /* found a series by the name of s */
	    gretl_errmsg_sprintf(_("'%s': not a scalar"), s);
	} else if (loop->parent != NULL && strlen(s) == gretl_namechar_spn(s)) {
	    /* potentially valid varname, but unknown at present */
	    *clr->vname = '\0';
	    strncat(clr->vname, s, VNAMELEN - 1);
	} else {
	    gretl_errmsg_sprintf(_("Undefined variable '%s' in loop condition."), s);
	    err = E_UNKVAR;
	}
    }

    return err;
}

#define maybe_date(s) (strchr(s, ':') || strchr(s, '/'))

static int parse_as_indexed_loop (LOOPSET *loop,
				  const double **Z,
				  const DATAINFO *pdinfo,
				  const char *lvar, 
				  const char *start,
				  const char *end)
{
    int err = 0;

    /* starting and ending values: try for dates first, then
       numeric constants, then variables */

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_indexed_loop: start='%s', end='%s'\n", start, end);
#endif

    if (maybe_date(start)) {
	loop->init.val = dateton(start, pdinfo);
	if (loop->init.val < 0) {
	    err = E_DATA;
	} else {
	    loop->init.val += 1;
	    loop->final.val = dateton(end, pdinfo);
	    if (loop->final.val < 0) {
		err = E_DATA;
	    } else {
		loop->final.val += 1;
		loop->type = DATED_LOOP;
	    }
	}
    } else {
	err = index_get_limit(loop, &loop->init, start, Z, pdinfo);
	if (!err) {
	    err = index_get_limit(loop, &loop->final, end, Z, pdinfo);
	}
	if (!err) {
	    loop->type = INDEX_LOOP;
	}
    }

    if (!err) {
	err = loop_attach_index_var(loop, lvar, pdinfo);
    }

#if LOOP_DEBUG
    fprintf(stderr, "indexed_loop: init.val=%g, final.val=%g, err=%d\n",
	    loop->init.val, loop->final.val, err);
#endif

    return err;
}

static int parse_as_count_loop (LOOPSET *loop, 
				const double **Z,
				const DATAINFO *pdinfo,
				const char *s)
{
    int err;

    err = index_get_limit(loop, &loop->final, s, Z, pdinfo);

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

static int 
set_forloop_element (char *s, LOOPSET *loop,
		     double ***pZ, DATAINFO *pdinfo,
		     int i)
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
				 const DATAINFO *pdinfo)
{
    int i, vi;
    int err;

#if LOOP_DEBUG
    fprintf(stderr, "list_vars_to_strings: adding %d strings\n", list[0]);
#endif

    err = allocate_each_strings(loop, list[0]);

    for (i=0; i<list[0] && !err; i++) {
	vi = list[i+1];
	if (vi < 0 || vi >= pdinfo->v) {
	    err = E_DATA;
	} else {
	    loop->eachstrs[i] = gretl_strdup(pdinfo->varname[vi]);
	    if (loop->eachstrs[i] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

/* At loop runtime, check the named list and insert the names (or
   numbers) of the variables as "eachstrs"; flag an error if the list
   has disappeared.
*/

static int loop_list_refresh (LOOPSET *loop, const DATAINFO *pdinfo)
{
    int *list = get_list_by_name(loop->listname);
    int err = 0;

    if (loop->eachstrs != NULL) {
	free_strings_array(loop->eachstrs, loop->itermax);
	loop->eachstrs = NULL;
    }

    loop->itermax = loop->final.val = 0;

    if (list == NULL) {
	err = E_UNKVAR;
    } else if (list[0] > 0) {
	err = list_vars_to_strings(loop, list, pdinfo);
	if (!err) {
	    loop->final.val = list[0];
	}
    }

    return err;
}

static int find_list_in_parentage (LOOPSET *loop, const char *s)
{
    char lname[VNAMELEN];
    int i;

    while ((loop = loop->parent) != NULL) {
	for (i=0; i<loop->n_cmds; i++) {
	    if (sscanf(loop->lines[i], "list %15[^ =]", lname)) {
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
each_strings_from_list_of_vars (LOOPSET *loop, const DATAINFO *pdinfo, 
				char *s, int *pnf, int type)
{
    int *list = NULL;
    int err = 0;

    if (type == WILDCARD_LIST) {
	s += strspn(s, " \t");
	list = varname_match_list(pdinfo, s);
	if (list == NULL) {
	    err = 1;
	}
    } else {
	char vn1[VNAMELEN], vn2[VNAMELEN];
	
	delchar(' ', s);

	if (sscanf(s, "%15[^.]..%15s", vn1, vn2) != 2) {
	    err = E_PARSE;
	} else {
	    int v1 = current_series_index(pdinfo, vn1);
	    int v2 = current_series_index(pdinfo, vn2);

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
		loop->eachstrs[i-1] = gretl_strdup(pdinfo->varname[vi]);
		if (loop->eachstrs[i-1] == NULL) {
		    free_strings_array(loop->eachstrs, list[0]);
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

static int
parse_as_each_loop (LOOPSET *loop, const DATAINFO *pdinfo, char *s)
{
    char ivar[VNAMELEN] = {0};
    int done = 0;
    int i, nf, err = 0;

    /* we're looking at the string that follows "loop foreach" */
    if (*s == '\0') {
	return E_PARSE;
    }
    
    s += strspn(s, " "); /* skip any spaces */

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: s = '%s'\n", s);
#endif 

    if (sscanf(s, "%15s", ivar) != 1) {
	return E_PARSE;
    } 

    s += strlen(ivar);
    nf = count_fields(s);

    if (nf == 0) {
	return E_PARSE;
    }

    if (nf <= 3 && strstr(s, "..") != NULL) {
	/* range of values, foo..quux */
	err = each_strings_from_list_of_vars(loop, pdinfo, s, &nf,
					     DOTTED_LIST);
	done = 1;
    } else if (nf == 1 && strchr(s, '*')) {
	err = each_strings_from_list_of_vars(loop, pdinfo, s, &nf,
					     WILDCARD_LIST);
	done = (err == 0);
    }	

    if (!done && nf == 1) {
	/* try for a named list? */
	err = list_loop_setup(loop, s, &nf);
	done = (err == 0);
    }

    if (!done) {
	/* simple list of values */
	err = allocate_each_strings(loop, nf);

	for (i=0; i<nf && !err; i++) {
	    int len;

	    while (isspace((unsigned char) *s)) s++;
	    len = strcspn(s, " ");

	    loop->eachstrs[i] = gretl_strndup(s, len);
	    if (loop->eachstrs[i] == NULL) {
		free_strings_array(loop->eachstrs, nf);
		loop->eachstrs = NULL;
		err = E_ALLOC;
	    } else {
		s += len;
	    }
	}
    }

    if (!err) {
	loop->type = EACH_LOOP;
	loop->init.val = 1;
	loop->final.val = nf;
	err = loop_attach_index_var(loop, ivar, pdinfo);
    }   

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: final.val=%g\n", loop->final.val);
#endif 

    return err;
}

/* try to parse out (expr1; expr2; expr3) */

static int 
parse_as_for_loop (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo, char *s)
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
	err = set_forloop_element(tmp, loop, pZ, pdinfo, j);
    }  

    if (!err && (sc != 2 || s != q)) {
	/* must have two semi-colons; must have reached rightmost ')' */
	err = E_PARSE;
    }

    free(tmp);

    if (!err) {
	loop->type = FOR_LOOP;
    }

    return err;
}

static int parse_first_loopline (char *s, LOOPSET *loop, 
				 double ***pZ, DATAINFO *pdinfo)
{
    char lvar[VNAMELEN], rvar[VNAMELEN], op[VNAMELEN];
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

    if (!strncmp(s, "foreach ", 8)) {
	err = parse_as_each_loop(loop, pdinfo, s + 8);
    } else if (sscanf(s, "%15[^= ] = %15[^.]..%15s", lvar, op, rvar) == 3) {
	err = parse_as_indexed_loop(loop, (const double **) *pZ, pdinfo, 
				    lvar, op, rvar);
    } else if (!strncmp(s, "for", 3)) {
	err = parse_as_for_loop(loop, pZ, pdinfo, s + 4);
    } else if (!strncmp(s, "while", 5)) {
	err = parse_as_while_loop(loop, pZ, pdinfo, s + 6);
    } else if (sscanf(s, "%15s", lvar) == 1) {
	err = parse_as_count_loop(loop, (const double **) *pZ, pdinfo, lvar);
    } else {
	printf("parse_first_loopline: failed on '%s'\n", s);
	gretl_errmsg_set(_("No valid loop condition was given."));
	err = 1;
    }

    return err;
}

/**
 * start_new_loop:
 * @s: loop specification line.
 * @inloop: current loop struct pointer, or %NULL.
 * @pZ: pointer to data array.
 * @pdinfo: data information struct.
 * @pZ: pointer to data array.
 * @nested: location to receive info on whether a new
 * loop was created, nested within the input loop.
 * @err: location to receive error code.
 *
 * Create a new LOOPSET based on the input line; this may or
 * may not be a child of @inloop.
 *
 * Returns: loop pointer on successful completion, %NULL on error.
 */

static LOOPSET *
start_new_loop (char *s, LOOPSET *inloop, 
		double ***pZ, DATAINFO *pdinfo,
		int *nested, int *err)
{
    LOOPSET *loop = NULL;

    gretl_error_clear();

#if LOOP_DEBUG
    fprintf(stderr, "start_new_loop: inloop=%p, line='%s'\n", 
	    (void *) inloop, s);
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

    *err = parse_first_loopline(s, loop, pZ, pdinfo);

    if (!*err) {
	*err = gretl_loop_prepare(loop);
    }

    if (*err) {
	free(loop->lines);
	free(loop->ci);
	free(loop);
	loop = NULL;
    } 

    return loop;
}

#if LOOP_DEBUG
# define MAX_FOR_TIMES  10
#else
# define MAX_FOR_TIMES  100000
#endif

static int loop_count_too_high (LOOPSET *loop)
{
    static int max_iters = 0;
    int nt = loop->iter + 1;

    if (loop->type == FOR_LOOP) {
	if (nt >= MAX_FOR_TIMES) {
	    gretl_errmsg_sprintf(_("Reached maximum iterations, %d"),
				 MAX_FOR_TIMES);
	    loop->err = 1;
	}
    } else {
	if (max_iters == 0) {
	    max_iters = libset_get_int(LOOP_MAXITER);
	}
	if (nt >= max_iters) {
	    gretl_errmsg_sprintf(_("Warning: no convergence after %d iterations"),
				 max_iters);
	    loop->err = 1;
	}
    }

    return loop->err;
}

/**
 * loop_condition:
 * @loop: pointer to loop commands struct.
 * @pZ: pointer to data array.
 * @pdinfo: data information struct.
 * @err: location to receive error code.
 *
 * Check whether a looping condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

static int 
loop_condition (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo, int *err)
{
    int ok = 0;

    if (loop->brk) {
	/* got "break" comand */
	loop->brk = 0;
	ok = 0;
    } else if (loop->type == COUNT_LOOP || indexed_loop(loop)) {
	if (loop->iter < loop->itermax) {
	    ok = 1;
	}
    } else {
	/* more complex forms of control (for, while) */
	if (loop_count_too_high(loop)) {
	    /* safeguard against infinite loops */
	    ok = 0;
	} else if (loop->type == FOR_LOOP) {
	    if (loop->iter > 0) {
		loop_delta(loop, pZ, pdinfo, err);
	    }
	    ok = loop_testval(loop, pZ, pdinfo, err);
	} else if (loop->type == WHILE_LOOP) {
	    ok = loop_testval(loop, pZ, pdinfo, err);
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

static int gretl_loop_prepare (LOOPSET *loop)
{
#ifdef ENABLE_GMP
    mpf_set_default_prec(256);
#endif

    /* allocate some initial lines/commands for loop */
    loop->lines = malloc(LOOP_BLOCK * sizeof *loop->lines); 
    loop->ci = malloc(LOOP_BLOCK * sizeof *loop->ci);
    
    if (loop->lines == NULL || loop->ci == NULL) {
	return 1;
    }

    return 0;
}

static void loop_model_free (LOOP_MODEL *lmod)
{
#if LOOP_DEBUG
    fprintf(stderr, "loop_model_free: lmod at %p, model0 at %p\n",
	    (void *) lmod, (void *) lmod->model0);
#endif

#ifdef ENABLE_GMP
    int i, n = 4 * lmod->model0->ncoeff;

    for (i=0; i<n; i++) {
	mpf_clear(lmod->bigarray[i]);
    }
#endif

    free(lmod->bigarray);
    free(lmod->cbak);
    free(lmod->sbak);
    free(lmod->cdiff);
    free(lmod->sdiff);

    gretl_model_free(lmod->model0);
}

/* Reset the loop model */

static void loop_model_zero (LOOP_MODEL *lmod)
{
    int i;

    lmod->n = 0;

    for (i=0; i<lmod->nc; i++) {
#ifdef ENABLE_GMP
	mpf_init(lmod->sum_coeff[i]);
	mpf_init(lmod->ssq_coeff[i]);
	mpf_init(lmod->sum_sderr[i]);
	mpf_init(lmod->ssq_sderr[i]);
#else
	lmod->sum_coeff[i] = lmod->ssq_coeff[i] = 0.0;
	lmod->sum_sderr[i] = lmod->ssq_sderr[i] = 0.0;
#endif
	lmod->cbak[i] = lmod->sbak[i] = NADBL;
	lmod->cdiff[i] = lmod->sdiff[i] = 0;
    }
}

/* Set everything in lmod to 0/null in case of failure */

static void loop_model_init (LOOP_MODEL *lmod, int lno)
{
    lmod->linenum = lno;
    lmod->nc = 0;
    lmod->model0 = NULL;
    lmod->bigarray = NULL;
    lmod->cbak = lmod->sbak = NULL;
    lmod->cdiff = lmod->sdiff = NULL;
}

/* Start up a LOOP_MODEL struct: copy pmod into place and
   allocate storage */

static int loop_model_start (LOOP_MODEL *lmod, const MODEL *pmod)
{
    int nc;

#if LOOP_DEBUG
    fprintf(stderr, "init: copying model at %p\n", (void *) pmod);
#endif

    lmod->nc = nc = pmod->ncoeff;

    lmod->model0 = gretl_model_copy(pmod);
    if (lmod->model0 == NULL) {
	return E_ALLOC;
    }

    lmod->bigarray = malloc(nc * 4 * sizeof *lmod->bigarray);
    if (lmod->bigarray == NULL) return E_ALLOC;

    lmod->sum_coeff = lmod->bigarray;
    lmod->ssq_coeff = lmod->sum_coeff + nc;
    lmod->sum_sderr = lmod->ssq_coeff + nc;
    lmod->ssq_sderr = lmod->sum_sderr + nc;

    lmod->cbak = malloc(nc * sizeof *lmod->cbak);
    if (lmod->cbak == NULL) goto cleanup;

    lmod->sbak = malloc(nc * sizeof *lmod->sbak);
    if (lmod->sbak == NULL) goto cleanup;

    lmod->cdiff = malloc(nc * sizeof *lmod->cdiff);
    if (lmod->cdiff == NULL) goto cleanup;

    lmod->sdiff = malloc(nc * sizeof *lmod->sdiff);
    if (lmod->sdiff == NULL) goto cleanup;

    loop_model_zero(lmod);

#if LOOP_DEBUG
    fprintf(stderr, " model copied to %p, returning 0\n", 
	    (void *) lmod->model0);
#endif

    return 0;

 cleanup:
    free(lmod->bigarray);
    free(lmod->cbak);
    free(lmod->sbak);
    free(lmod->cdiff);
    free(lmod->sdiff);

    return E_ALLOC;
}

static void loop_print_free (LOOP_PRINT *lprn)
{
#ifdef ENABLE_GMP
    int i;

    for (i=0; i<lprn->list[0]; i++) {
	mpf_clear(lprn->sum[i]);
	mpf_clear(lprn->ssq[i]);
    }
#endif

    free(lprn->sum);
    free(lprn->ssq);
    free(lprn->list);
    free(lprn->xbak);
    free(lprn->diff);
}

static void loop_print_zero (LOOP_PRINT *lprn)
{
    int i;

    lprn->n = 0;

    for (i=0; i<lprn->list[0]; i++) { 
#ifdef ENABLE_GMP
	mpf_init(lprn->sum[i]);
	mpf_init(lprn->ssq[i]);
#else
	lprn->sum[i] = lprn->ssq[i] = 0.0;
#endif
	lprn->xbak[i] = NADBL;
	lprn->diff[i] = 0;
    }
}

static void loop_print_init (LOOP_PRINT *lprn, int lno)
{
    lprn->linenum = lno;
    lprn->list = NULL;
    lprn->sum = NULL;
    lprn->ssq = NULL;
    lprn->xbak = NULL;
    lprn->diff = NULL;
}

/* allocate and initialize @lprn, based on the number of
   elements in @list */

static int loop_print_start (LOOP_PRINT *lprn, const int *list)
{
    int k = list[0];

    lprn->list = gretl_list_copy(list);
    if (lprn->list == NULL) return 1;

    lprn->sum = malloc(k * sizeof *lprn->sum);
    if (lprn->sum == NULL) goto cleanup;

    lprn->ssq = malloc(k * sizeof *lprn->ssq);
    if (lprn->ssq == NULL) goto cleanup;

    lprn->xbak = malloc(k * sizeof *lprn->xbak);
    if (lprn->xbak == NULL) goto cleanup;

    lprn->diff = malloc(k * sizeof *lprn->diff);
    if (lprn->diff == NULL) goto cleanup;

    loop_print_zero(lprn);

    return 0;

 cleanup:

    free(lprn->list);
    free(lprn->sum);
    free(lprn->ssq);
    free(lprn->xbak);
    free(lprn->diff);

    return 1;
}

static LOOP_PRINT *get_loop_print_by_line (LOOPSET *loop, int lno, int *err)
{
    LOOP_PRINT *prns;
    int i, np = loop->n_prints;

    for (i=0; i<np; i++) {
	if (loop->prns[i].linenum == lno) {
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
    destroy_dataset(lstore->Z, lstore->dinfo);
    lstore->Z = NULL;
    lstore->dinfo = NULL;
    free(lstore->fname);
    lstore->fname = NULL;

    lstore->linenum = -1;
    lstore->n = 0;
    lstore->opt = OPT_NONE;
}

static int loop_store_set_filename (LOOPSET *loop, const char *fname,
				    gretlopt opt)
{
    if (fname == NULL || *fname == '\0') {
	return E_ARGS;
    }

    loop->store.fname = gretl_strdup(fname);
    if (loop->store.fname == NULL) {
	return E_ALLOC;
    }

    if (opt == OPT_NONE) {
	opt = data_save_opt_from_suffix(loop->store.fname);
    }

    loop->store.opt = opt;    

    return 0;
}

/* check, allocate and initialize loop data storage */

static int loop_store_start (LOOPSET *loop, const int *list, 
			     const char *fname, DATAINFO *pdinfo,
			     gretlopt opt)
{
    int i, n, err = 0;

    if (list == NULL || list[0] == 0) {
	gretl_errmsg_set("'store' list is empty");
	return E_DATA;
    }

    err = loop_store_set_filename(loop, fname, opt);
    if (err) {
	return err;
    }

    n = (loop->itermax > 0)? loop->itermax : DEFAULT_NOBS;

    loop->store.dinfo = create_auxiliary_dataset(&loop->store.Z, list[0] + 1, n);
    if (loop->store.dinfo == NULL) {
	return E_ALLOC;
    }
    
#if LOOP_DEBUG
    fprintf(stderr, "loop_store_init: created sZ, v = %d, n = %d\n",
	    loop->store.dinfo->v, loop->store.dinfo->n);
#endif

    for (i=1; i<=list[0] && !err; i++) {
	const char *s = gretl_scalar_get_name(list[i]);

	if (s == NULL) {
	    err = E_DATA;
	} else {
	    strcpy(loop->store.dinfo->varname[i], s);
	}
    }

    return err;
}

static int loop_store_update (LOOPSET *loop, int lno,
			      const int *list, const char *fname,
			      const double **Z, DATAINFO *pdinfo,
			      gretlopt opt)
{
    int i, t, err = 0;

    if (loop->store.linenum >= 0 && loop->store.linenum != lno) {
	gretl_errmsg_set("Only one 'store' command is allowed in a "
			 "progressive loop");
	return E_DATA;
    }

    if (loop->store.Z == NULL) {
	/* not started yet */
	err = loop_store_start(loop, list, fname, pdinfo, opt);
	if (err) {
	    return err;
	}
    }

    loop->store.linenum = lno;
    t = loop->store.n;

    if (t >= loop->store.dinfo->n) {
	if (extend_loop_dataset(&loop->store)) {
	    return E_ALLOC;
	}
    }

    for (i=1; i<=list[0]; i++) {
	loop->store.Z[i][t] = gretl_scalar_get_value_by_index(list[i]);
    }

    loop->store.n += 1;

    return 0;
}

/* See if we already have a model recorder in place for the command on
   line @lno of the loop.  If so, return it, else create a new one and
   return it.
*/

static MODEL *get_model_record_by_line (LOOPSET *loop, int lno, int *err)
{
    MODEL **models, *pmod;
    int *modlines;
    int nm = loop->n_models;
    int i;

    for (i=0; i<nm; i++) {
	if (lno == loop->model_lines[i]) {
	    return loop->models[i];
	}
    }

    modlines = realloc(loop->model_lines, (nm + 1) * sizeof *modlines);
    if (modlines == NULL) {
	*err = E_ALLOC;
	return NULL;
    } else {
	loop->model_lines = modlines;
    }

    models = realloc(loop->models, (nm + 1) * sizeof *models);
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

    loop->model_lines[nm] = lno;
    pmod->ID = nm + 1;
    loop->models[nm] = pmod;
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
    int nlm = loop->n_loop_models;
    int i;

#if LOOP_DEBUG
    fprintf(stderr, "get_loop_model_by_line: loop->n_loop_models = %d\n",
	    loop->n_loop_models);
#endif

    for (i=0; i<nlm; i++) {
	if (loop->lmodels[i].linenum == lno) {
	    return &loop->lmodels[i];
	}
    }

    lmods = realloc(loop->lmodels, (nlm + 1) * sizeof *loop->lmodels);
    if (lmods == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    loop->lmodels = lmods;
    loop_model_init(&loop->lmodels[nlm], lno);
    loop->n_loop_models += 1;

    return &loop->lmodels[nlm];
}

#define realdiff(x,y) (fabs((x)-(y)) > 2.0e-13)

/* Update the info stored in LOOP_MODEL based on the results in pmod.
   If this is the first use we have to do some allocation first.
*/

static int loop_model_update (LOOP_MODEL *lmod, MODEL *pmod)
{
#ifdef ENABLE_GMP
    mpf_t m;
#endif
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
    }

#ifdef ENABLE_GMP
    mpf_init(m);
#endif

    for (j=0; j<pmod->ncoeff; j++) {
#ifdef ENABLE_GMP
	mpf_set_d(m, pmod->coeff[j]);
	mpf_add(lmod->sum_coeff[j], lmod->sum_coeff[j], m); 
	mpf_mul(m, m, m);
	mpf_add(lmod->ssq_coeff[j], lmod->ssq_coeff[j], m);

	mpf_set_d(m, pmod->sderr[j]);
	mpf_add(lmod->sum_sderr[j], lmod->sum_sderr[j], m);
	mpf_mul(m, m, m);
	mpf_add(lmod->ssq_sderr[j], lmod->ssq_sderr[j], m);
#else
	lmod->sum_coeff[j] += pmod->coeff[j];
	lmod->ssq_coeff[j] += pmod->coeff[j] * pmod->coeff[j];
	lmod->sum_sderr[j] += pmod->sderr[j];
	lmod->ssq_sderr[j] += pmod->sderr[j] * pmod->sderr[j];
#endif
	if (!na(lmod->cbak[j]) && realdiff(pmod->coeff[j], lmod->cbak[j])) {
	    lmod->cdiff[j] = 1;
	}
	if (!na(lmod->sbak[j]) && realdiff(pmod->sderr[j], lmod->sbak[j])) {
	    lmod->sdiff[j] = 1;
	}
	lmod->cbak[j] = pmod->coeff[j];
	lmod->sbak[j] = pmod->sderr[j];
    }

#ifdef ENABLE_GMP
    mpf_clear(m);
#endif

    lmod->n += 1;

#if LOOP_DEBUG
    fprintf(stderr, "loop_model_update: returning %d\n", err);
#endif

    return err;
}

/* Update the LOOP_PRINT struct lprn using the current values from the
   variables in list.  If this is the fist use we need to do some
   allocation first.
*/

static int loop_print_update (LOOP_PRINT *lprn,
			      const int *list, const double **Z, 
			      const DATAINFO *pdinfo)
{
#ifdef ENABLE_GMP
    mpf_t m;
#endif
    int j, vj;
    double x;
    int err = 0;

    if (lprn->list == NULL) {
	/* not started yet */
	err = loop_print_start(lprn, list);
	if (err) {
	    return err;
	}
    }

#ifdef ENABLE_GMP
    mpf_init(m);
#endif
    
    for (j=0; j<list[0]; j++) {
	vj = list[j+1];
	x = gretl_scalar_get_value_by_index(vj);
#ifdef ENABLE_GMP
	mpf_set_d(m, x); 
	mpf_add(lprn->sum[j], lprn->sum[j], m);
	mpf_mul(m, m, m);
	mpf_add(lprn->ssq[j], lprn->ssq[j], m);
#else
	lprn->sum[j] += x;
	lprn->ssq[j] += x * x;
#endif
	if (!na(lprn->xbak[j]) && realdiff(x, lprn->xbak[j])) {
	    lprn->diff[j] = 1;
	}
	lprn->xbak[j] = x;
    }

#ifdef ENABLE_GMP
    mpf_clear(m);
#endif

    lprn->n += 1;

    return err;
}

static int add_more_loop_lines (LOOPSET *loop)
{
    int nb = 1 + (loop->n_cmds + 1) / LOOP_BLOCK;
    char **lines;
    int *ci;
    
    lines = realloc(loop->lines, (nb * LOOP_BLOCK) * sizeof *lines); 
    ci = realloc(loop->ci, (nb * LOOP_BLOCK) * sizeof *ci);
    
    if (lines == NULL || ci == NULL) {
	return 1;
    }

    loop->lines = lines;
    loop->ci = ci;

    return 0;
}  

static int real_append_line (ExecState *s, LOOPSET *loop)
{
    const char *flagstr = NULL;
    int nc = loop->n_cmds;
    int len, err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "real_append_line: s->line = '%s'\n", s->line);
#endif

    if ((nc + 1) % LOOP_BLOCK == 0) {
	if (add_more_loop_lines(loop)) {
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
	    loop->lines[nc] = malloc(len + 1);
	} else {
	    loop->lines[nc] = gretl_strdup(s->line);
	}
	if (loop->lines[nc] == NULL) {
	    err = E_ALLOC;
	} else if (flagstr != NULL) {
	    sprintf(loop->lines[nc], "%s%s", s->line, flagstr);
	}
    }

    if (!err) {
	if (s->cmd->ci == PRINT && (!loop_is_progressive(loop) || strchr(s->line, '"'))) {
	    loop->ci[nc] = 0;
	} else {
	    loop->ci[nc] = s->cmd->ci;
	}
	loop->n_cmds += 1;
    }

#if LOOP_DEBUG
    fprintf(stderr, "loop %p: n_cmds=%d, line[%d]='%s', ci=%d\n",
	    (void *) loop, loop->n_cmds, nc, loop->lines[nc],
	    loop->ci[nc]);
#endif

    return err;
}  

/**
 * gretl_loop_append_line:
 * @s: program execution state.
 * @pZ: pointer to data matrix.
 * @pdinfo: dataset information.
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

int gretl_loop_append_line (ExecState *s, double ***pZ, 
			    DATAINFO *pdinfo)
{
    LOOPSET *loop = currloop;
    LOOPSET *newloop = currloop;
    int err = 0;

    gretl_error_clear();

#if LOOP_DEBUG > 1
    fprintf(stderr, "gretl_loop_append_line: currloop = %p, line = '%s'\n", 
	    (void *) loop, s->line);
#endif

    if (s->cmd->ci == LOOP) {
	gretlopt opt = get_loop_opts(s, &err);
	int nested = 0;

	if (!err) {
	    newloop = start_new_loop(s->line, loop, pZ, pdinfo, 
				     &nested, &err);
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

static void print_loop_coeff (const DATAINFO *pdinfo, 
			      const LOOP_MODEL *lmod, 
			      int i, PRN *prn)
{
    char pname[VNAMELEN];
#ifdef ENABLE_GMP
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

    gretl_model_get_param_name(lmod->model0, pdinfo, i, pname);
    pprintf(prn, "%*s", VNAMELEN - 1, pname);
    pprintf(prn, "%#14g %#14g %#14g %#14g\n", mpf_get_d(c1), mpf_get_d(sd1), 
	    mpf_get_d(c2), mpf_get_d(sd2));

    mpf_clear(c1);
    mpf_clear(c2);
    mpf_clear(m);
    mpf_clear(sd1);
    mpf_clear(sd2);
#else /* non-GMP */
    bigval m1, m2, sd1, sd2;
    int n = lmod->n;

    m1 = lmod->sum_coeff[i] / n;
    if (lmod->cdiff[i] == 0) {
	sd1 = 0.0;
    } else {
	sd1 = (lmod->ssq_coeff[i] - n * m1 * m1) / n;
	sd1 = (sd1 <= 0.0)? 0.0 : sqrt((double) sd1);
    }

    m2 = lmod->sum_sderr[i] / n;
    if (lmod->sdiff[i] == 0) {
	sd2 = 0.0;
    } else {
	sd2 = (lmod->ssq_sderr[i] - n * m2 * m2) / n;
	sd2 = (sd2 <= 0.0)? 0 : sqrt((double) sd2);
    }

    gretl_model_get_param_name(lmod->model0, pdinfo, i, pname);
    pprintf(prn, "%*s", VNAMELEN - 1, pname);
    pprintf(prn, "%#14g %#14g %#14g %#14g\n", (double) m1, (double) sd1, 
	    (double) m2, (double) sd2);
#endif
}

static void print_loop_model (LOOP_MODEL *lmod, const DATAINFO *pdinfo, 
			      PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    int i;

    ntodate(startdate, lmod->model0->t1, pdinfo);
    ntodate(enddate, lmod->model0->t2, pdinfo);

    pputc(prn, '\n');
    pprintf(prn, _("%s estimates using the %d observations %s-%s\n"),
	    _(estimator_string(lmod->model0, prn)), lmod->model0->nobs, 
	    startdate, enddate);
    print_model_vcv_info(lmod->model0, prn);
    pprintf(prn, _("Statistics for %d repetitions\n"), lmod->n); 
    pprintf(prn, _("Dependent variable: %s\n\n"), 
	    gretl_model_get_depvar_name(lmod->model0, pdinfo));

    pputs(prn, _("                     mean of      std. dev. of     mean of"
		 "     std. dev. of\n"
		 "                    estimated      estimated"
		 "      estimated      estimated\n"
		 "      Variable     coefficients   coefficients   std. errors"
		 "    std. errors\n\n"));

    for (i=0; i<lmod->model0->ncoeff; i++) {
	print_loop_coeff(pdinfo, lmod, i, prn);
    }

    pputc(prn, '\n');
}

static void print_loop_prn (LOOP_PRINT *lprn, const DATAINFO *pdinfo, 
			    PRN *prn)
{
    bigval mean, m, sd;
    int i, vi, n;
    const char *s;

    if (lprn == NULL) {
	return;
    }

    n = lprn->n;

    pprintf(prn, _("Statistics for %d repetitions\n"), n); 
    pputs(prn, "    ");
    pputs(prn, _("   Variable     mean         std. dev.\n"));

#ifdef ENABLE_GMP
    mpf_init(mean);
    mpf_init(m);
    mpf_init(sd);
    
    for (i=0; i<lprn->list[0]; i++) {
	vi = lprn->list[i+1];
	s = gretl_scalar_get_name(vi);
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
	pprintf(prn, "%*s", VNAMELEN - 1, s);
	pprintf(prn, "%#14g %#14g\n", mpf_get_d(mean), mpf_get_d(sd));
    }

    mpf_clear(mean);
    mpf_clear(m);
    mpf_clear(sd);
#else
    for (i=0; i<lprn->list[0]; i++) {
	vi = lprn->list[i+1];
	s = gretl_scalar_get_name(vi);
	mean = lprn->sum[i] / n;
	if (lprn->diff[i] == 0) {
	    sd = 0.0;
	} else {
	    m = (lprn->ssq[i] - n * mean * mean) / n;
	    sd = (m < 0)? 0 : sqrt((double) m);
	}
	pprintf(prn, "%*s", VNAMELEN - 1, s);
	pprintf(prn, "%#14g %#14g\n", (double) mean, (double) sd);
    }
#endif
    pputc(prn, '\n');
}

static int loop_store_save (LOOP_STORE *lstore, PRN *prn)
{
    int *list;
    int err = 0;

    list = gretl_consecutive_list_new(1, lstore->dinfo->v - 1);
    if (list == NULL) {
	return E_ALLOC;
    }

    lstore->dinfo->t2 = lstore->n - 1;

    pprintf(prn, _("store: using filename %s\n"), lstore->fname);

    err = write_data(lstore->fname, list, (const double **) lstore->Z, 
		     lstore->dinfo, lstore->opt, 0);

    if (!err) {
	pprintf(prn, _("Data written OK.\n"));
    } else {
	pprintf(prn, _("write of data file failed\n"));
    }

    free(list);

    return err;
}

/**
 * print_loop_results:
 * @loop: pointer to loop struct.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Print out the results after completion of the loop @loop.
 */

static void print_loop_results (LOOPSET *loop, const DATAINFO *pdinfo, 
				PRN *prn)
{
    char linecpy[MAXLINE];
    int iters = loop->iter;
    int i, j, k;

    if (loop->type != COUNT_LOOP && !(loop_is_quiet(loop))) {
	pprintf(prn, _("\nNumber of iterations: %d\n\n"), iters);
    }

    j = 0;
    k = 0;

    for (i=0; i<loop->n_cmds; i++) {
	gretlopt opt = OPT_NONE;

#if LOOP_DEBUG
	fprintf(stderr, "print_loop_results: loop command %d (i=%d): %s\n", 
		i+1, i, loop->lines[i]);
#endif

	if (plain_model_ci(loop->ci[i])) {
	    strcpy(linecpy, loop->lines[i]);
	    opt = get_gretl_options(linecpy, NULL);
	}	    

	if (!loop_is_progressive(loop) && loop->ci[i] == OLS) {
	    if (model_print_deferred(opt)) {
		MODEL *pmod = loop->models[j++];

		set_model_id(pmod);
		printmodel(pmod, pdinfo, opt, prn);
	    }	    
	}

	if (loop_is_progressive(loop)) {
	    if (plain_model_ci(loop->ci[i]) && !(opt & OPT_Q)) {
		print_loop_model(&loop->lmodels[j], pdinfo, prn);
		loop_model_zero(&loop->lmodels[j]);
		j++;
	    } else if (loop->ci[i] == PRINT) {
		print_loop_prn(&loop->prns[k], pdinfo, prn);
		loop_print_zero(&loop->prns[k]);
		k++;
	    } else if (loop->ci[i] == STORE) {
		loop_store_save(&loop->store, prn);
	    }
	}
    }
}

static int 
substitute_dollar_targ (char *str, const LOOPSET *loop,
			const double **Z, const DATAINFO *pdinfo,
			int *subst)
{
    char targ[VNAMELEN + 3] = {0};
    char ins[32];
    char *p, *pins;
    int targlen;
    double forval = 0;
    int idx = 0;
    int err = 0;

#if SUBST_DEBUG
    fprintf(stderr, "subst_dollar_targ:\n original: '%s'\n", str);
#endif

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
	return 1;
    }

#if SUBST_DEBUG
    fprintf(stderr, " target = '%s', idx = %d\n", targ, idx);
#endif

    if (strstr(str, targ) == NULL) {
	/* nothing to be done */
	return 0;
    }

    pins = ins;

    /* prepare substitute */
    if (loop->type == FOR_LOOP) {
	forval = gretl_scalar_get_value(loop->init.vname);
	/* the rest is handled below */
    } else if (loop->type == INDEX_LOOP) {
	sprintf(ins, "%d", idx);
    } else if (loop->type == DATED_LOOP) {
	ntodate(ins, idx, pdinfo);
    } else if (loop->type == EACH_LOOP) {
	pins = loop->eachstrs[idx - 1];
    }  

    while ((p = strstr(str, targ)) != NULL) {
	char *q = malloc(strlen(p));

	if (q == NULL) {
	    err = 1;
	    break;
	}

	strcpy(q, p + targlen);

	if (loop->type == FOR_LOOP) {
	    sprintf(ins, "%g", forval);
	    if (p - str > 0 && *(p - 1) == '[' && *(p + targlen) == ']') {
		/* got an obs-type string, on the pattern [$lvar] */
		int t = dateton(ins, pdinfo);

		if (t < 0) {
		    t = atoi(ins) - 1;
		}
		sprintf(ins, "%d", t);
	    } 
	} 

	strcpy(p, pins);
	strcpy(p + strlen(pins), q);
	free(q);

	*subst = 1;
    }

#if SUBST_DEBUG
    fprintf(stderr, " after: '%s'\n", str);
#endif

    return err;
}

static int extend_loop_dataset (LOOP_STORE *lstore)
{
    double *x;
    int oldn = lstore->dinfo->n;
    int n = oldn + DEFAULT_NOBS;
    int i, t;

    for (i=0; i<lstore->dinfo->v; i++) {
	x = realloc(lstore->Z[i], n * sizeof *x);
	if (x == NULL) {
	    return E_ALLOC;
	}
	lstore->Z[i] = x;
	for (t=oldn; t<n; t++) {
	    lstore->Z[i][t] = (i == 0)? 1.0 : NADBL;
	}	    
    }
    
    lstore->dinfo->n = n;
    lstore->dinfo->t2 = n - 1;

    ntodate(lstore->dinfo->endobs, n - 1, lstore->dinfo);

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

static int 
top_of_loop (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
    int err = 0;

    loop->iter = 0;

    if (is_list_loop(loop)) {
	err = loop_list_refresh(loop, pdinfo);
    } else if (loop->type == INDEX_LOOP) {
	loop->init.val = controller_get_val(&loop->init);
    } else if (loop->type == FOR_LOOP) {
	forloop_init(loop, pZ, pdinfo, &err);
    }

    if (!err && (loop->type == COUNT_LOOP || indexed_loop(loop))) {
	loop->final.val = controller_get_val(&loop->final);
	if (na(loop->init.val) || na(loop->final.val)) {
	    err = 1;
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
	    gretl_scalar_set_value(loop->idxname, loop->idxval);
	}

	/* initialization in case this loop is being run more than
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
print_loop_progress (const LOOPSET *loop, const DATAINFO *pdinfo,
		     PRN *prn)
{
    int i = loop->init.val + loop->iter;

    if (loop->type == INDEX_LOOP) {
	pprintf(prn, "loop: %s = %d\n\n", loop->idxname, i);
    } else if (loop->type == DATED_LOOP) {
	char obs[OBSLEN];

	ntodate(obs, i, pdinfo);
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
make_dollar_substitutions (char *str, const LOOPSET *loop,
			   const double **Z, const DATAINFO *pdinfo,
			   int *subst)
{
    int err = 0;

    if (indexed_loop(loop) || loop->type == FOR_LOOP) {
	err = substitute_dollar_targ(str, loop, Z, pdinfo, subst);
    }

    while (!err && (loop = subst_loop_in_parentage(loop)) != NULL) {
	err = substitute_dollar_targ(str, loop, Z, pdinfo, subst);
    }

    return err;
}

/* Try to determine if a "list" command within a loop modifies the
   list that is controlling the loop: if so, we'll have to arrange to
   have the list of variable names refreshed at the top of the loop.
*/

static int modifies_loop_list (const LOOPSET *loop, const char *s)
{
    int ret = 0;

    if (!strncmp(s, "list", 4)) {
	char lname[VNAMELEN];

	s += 4;
	s += strspn(s, " ");
	*lname = '\0';
	sscanf(s, "%15[^ +=-]", lname);
	ret = !strcmp(lname, loop->listname);
    }

    return ret;
}

static int maybe_refresh_list (CMD *cmd, const LOOPSET *loop,
			       const char *line)
{
    int ret = 0;

    if (cmd->ci == GENR) {
	ret = modifies_loop_list(loop, line);
    }

    return ret;
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

static int add_loop_genr (LOOPSET *loop, GENERATOR *genr, int lno)
{
    GENERATOR **genrs;
    int n = loop->n_genrs;

    genrs = realloc(loop->genrs, (n+1) * sizeof *genrs);
    if (genrs == NULL) {
	return E_ALLOC;
    } 

    loop->genrs = genrs;
    loop->genrs[n] = genr;
    loop->n_genrs = n + 1;
    genr_set_loopline(genr, lno);

    return 0;
}

static GENERATOR *get_loop_genr_by_line (LOOPSET *loop, int lno, 
					 const char *line, 
					 double ***pZ, DATAINFO *pdinfo,
					 int *err)
{
    GENERATOR *genr;
    int i, ll;

    for (i=0; i<loop->n_genrs; i++) {
	ll = genr_get_loopline(loop->genrs[i]);
	if (ll == lno) {
	    return loop->genrs[i];
	}
    }

    genr = genr_compile(line, pZ, pdinfo, OPT_NONE, err);

    if (!*err) {
	*err = add_loop_genr(loop, genr, lno);
    }

    return genr;
}

/* get the next command for a loop by pulling a line off the
   stack of loop commands.
*/

static int next_command (char *targ, LOOPSET *loop, int *pj)
{
    int ret = 1, j = *pj;

    if (j < loop->n_cmds) {
	strcpy(targ, loop->lines[j++]);
	*pj = j;
    } else {
	ret = 0;
    }

    return ret;
}

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

int gretl_loop_exec (ExecState *s, double ***pZ, DATAINFO *pdinfo) 
{
    GENERATOR *genr;
    LOOPSET *loop = currloop;
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    MODEL *pmod;
    LOOP_MODEL *lmod;
    LOOP_PRINT *lprn;
    char errline[MAXLINE];
    int indent0, subst, lrefresh;
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

    set_loop_on(loop_is_quiet(loop)); /* libset.c */

#if LOOP_DEBUG
    fprintf(stderr, "loop_exec: loop = %p\n", (void *) loop);
#endif

    err = top_of_loop(loop, pZ, pdinfo);
    
    while (!err && loop_condition(loop, pZ, pdinfo, &err)) {
#if LOOP_DEBUG
	fprintf(stderr, "top of loop: iter = %d\n", loop->iter);
#endif
	j = lrefresh = subst = 0;

	pmod = NULL;
	lmod = NULL;
	lprn = NULL;

	if (gretl_echo_on() && indexed_loop(loop) && !loop_is_quiet(loop)) {
	    print_loop_progress(loop, pdinfo, prn);
	}

	while (!err && next_command(line, loop, &j)) {
#if LOOP_DEBUG
	    fprintf(stderr, " j=%d, line='%s'\n", j, line);
#endif
	    strcpy(errline, line);
	    err = make_dollar_substitutions(line, loop, 
					    (const double **) *pZ,
					    pdinfo, &subst);
	    if (err) {
		break;
	    }

	    if (loop_is_progressive(loop)) {
		cmd->flags |= CMD_PROG;
	    }

	    /* We already have the "ci" index recorded, but here
	       we do some further parsing. 
	    */
	    err = parse_command_line(line, cmd, pZ, pdinfo);

	    if (cmd->ci < 0) {
		continue;
	    } else if (err) {
		if (libset_get_bool(HALT_ON_ERR)) {
		    break;
		} else {
		    /* try soldiering on */
		    errmsg(err, prn);
		    err = 0;
		    continue;
		}
	    }

	    if (!subst && cmd_subst(cmd)) {
		subst = 1;
	    }

	    if (is_list_loop(loop) && maybe_refresh_list(cmd, loop, line)) {
		lrefresh = 1;
	    }

	    if (gretl_echo_on()) {
		if (s->cmd->ci == ENDLOOP) {
		    if (indexed_loop(loop)) {
			pputc(prn, '\n');
		    }
		} else if (!loop_is_quiet(loop)) {
		    echo_cmd(cmd, pdinfo, line, CMD_BATCH_MODE, prn);
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
		    err = gretl_loop_exec(s, pZ, pdinfo);
		}
	    } else if (cmd->ci == BREAK) {
		loop->brk = 1;
	    } else if (cmd->ci == ENDLOOP) {
		; /* implicit break */
	    } else if (cmd->ci == FREQ) {
		err = freqdist(cmd->list[1], (const double **) *pZ, pdinfo, 0, 
			       cmd->opt, prn);
	    } else if (plain_model_ci(cmd->ci)) {
		/* model may need special handling */
		if (loop_is_progressive(loop) && !(cmd->opt & OPT_Q)) {
		    lmod = get_loop_model_by_line(loop, j, &err);
		} else if (model_print_deferred(cmd->opt)) {
		    pmod = get_model_record_by_line(loop, j, &err);
		}
		/* estimate the model called for */
		if (!err) {
		    err = gretl_cmd_exec(s, pZ, pdinfo);
		}
		if (!err) {
		    int moderr = check_gretl_errno();

		    if (moderr) {
			if (loop_is_progressive(loop) || 
			    model_print_deferred(cmd->opt)) {
			    err = moderr;
			} else {
			    errmsg(moderr, prn);
			}
		    } else if (loop_is_progressive(loop) && !(cmd->opt & OPT_Q)) {
			err = loop_model_update(lmod, s->models[0]);
			set_as_last_model(s->models[0], GRETL_OBJ_EQN);
		    } else if (model_print_deferred(cmd->opt)) {
			swap_models(s->models[0], pmod);
			pmod->ID = j;
			set_as_last_model(pmod, GRETL_OBJ_EQN);
			model_count_minus();
		    } else {
			if (!(cmd->opt & OPT_Q)) {
			    printmodel(s->models[0], pdinfo, cmd->opt, prn);
			}
			set_as_last_model(s->models[0], GRETL_OBJ_EQN);
		    }
		}
	    } else if (cmd->ci == PRINT && *cmd->param == '\0' &&
		       loop_is_progressive(loop)) {
		lprn = get_loop_print_by_line(loop, j, &err);
		if (!err) {
		    loop_print_update(lprn, cmd->list, (const double **) *pZ, 
				      pdinfo);
		}
	    } else if (cmd->ci == STORE && loop_is_progressive(loop)) {
		err = loop_store_update(loop, j, cmd->list, cmd->param,
					(const double **) *pZ, 
					pdinfo, cmd->opt);
	    } else if (loop_is_progressive(loop) && not_ok_in_progloop(cmd->ci)) {
		gretl_errmsg_sprintf(_("%s: not implemented in 'progressive' loops"),
				     cmd->word);
		err = 1;
	    } else if (cmd->ci == GENR) {
		if (!loop_is_verbose(loop)) {
		    cmd->opt |= OPT_Q;
		}
		if (subst || (cmd->opt & OPT_U)) {
		    /* can't use a "compiled" genr if string substitution
		       has been done, since the genr expression will not
		       be constant */
		    err = generate(line, pZ, pdinfo, cmd->opt, prn);
		} else {
		    genr = get_loop_genr_by_line(loop, j, line, pZ, pdinfo, &err);
		    if (!err) {
			err = execute_genr(genr, pZ, pdinfo, OPT_L, prn);
		    }
		}
	    } else {
		err = gretl_cmd_exec(s, pZ, pdinfo);
		if (!err && block_model(cmd)) {
		    /* NLS, etc. */
		    printmodel(s->models[0], pdinfo, cmd->opt, prn);
		    set_as_last_model(s->models[0], GRETL_OBJ_EQN);
		}
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
	    err = gretl_if_state_check(indent0);
	}

	if (!err) {
	    loop->iter += 1;
	    if (indexed_loop(loop)) {
		loop->idxval += 1;
		gretl_scalar_set_value(loop->idxname, loop->idxval);
	    } else if (lrefresh) {
		/* added 2008-01-11, AC */
		loop_list_refresh(loop, pdinfo);
	    }
	}

    } /* end iterations of loop */

    if (err) {
	if (!s->funcerr) {
	    errmsg(err, prn);
	    pprintf(prn, ">> %s\n", errline);
	}
    } else if (loop->err) {
	errmsg(loop->err, prn);
	err = loop->err;
    }

    if (!err && loop->iter > 0) {
	print_loop_results(loop, pdinfo, prn); 
    }

    if (loop->n_models > 0) {
	/* we need to update models[0] */
	GretlObjType type;
	void *ptr = get_last_model(&type);
	int i;

	if (type == GRETL_OBJ_EQN && s->models[0] != ptr) {
	    swap_models(s->models[0], loop->models[loop->n_models - 1]);
	    set_as_last_model(s->models[0], GRETL_OBJ_EQN);
	}
	for (i=0; i<loop->n_models; i++) {
	    gretl_model_free(loop->models[i]);
	}
    }

    if (line != NULL) {
	*line = '\0';
    } 

    if (loop->parent == NULL) {
	/* reached top of stack: clean up */
	gretl_loop_destroy(loop);
	currloop = NULL;
	set_loop_off();
    }

    return (libset_get_bool(HALT_ON_ERR))? err : 0;
}



