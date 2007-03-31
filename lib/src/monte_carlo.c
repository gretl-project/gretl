/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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

/*  monte_carlo.c - loop simulation procedures */

#include "libgretl.h" 
#include "monte_carlo.h"
#include "loop_private.h"
#include "libset.h"
#include "compat.h"
#include "cmd_private.h"
#include "var.h"
#include "objstack.h"
#include "gretl_func.h"

#include <time.h>
#include <unistd.h>

#define LOOP_DEBUG 0
#define SUBST_DEBUG 0
#define IDX_DEBUG 0

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

enum loop_val_codes {
    LOOP_VAL_BAD = 0,
    LOOP_VAL_UNDEF = -9999
};

#define indexed_loop(l) (l->type == INDEX_LOOP || \
                         l->type == DATED_LOOP || \
			 l->type == EACH_LOOP)

typedef struct {
    int *list;
    bigval *sum;
    bigval *ssq;
} LOOP_PRINT;  

/* below: used for special "progressive" loop */ 

typedef struct {
    int nc;                 /* number of coefficients */
    MODEL *model0;          /* copy of initial model */
    bigval *sum_coeff;      /* sums of coefficient estimates */
    bigval *ssq_coeff;      /* sums of squares of coeff estimates */
    bigval *sum_sderr;      /* sums of estimated std. errors */
    bigval *ssq_sderr;      /* sums of squares of estd std. errs */
} LOOP_MODEL;

enum loop_flags {
    LOOP_PROGRESSIVE = 1 << 0,
    LOOP_VERBOSE     = 1 << 1,
    LOOP_QUIET       = 1 << 2
};

struct controller_ {
    double val;
    int vnum;
    char vname[VNAMELEN];
    int vsign;
    char *expr;
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

    /* control variables */
    char ichar;
    int ival;
    char brk;

    controller init;
    controller test;
    controller delta;
    controller final;

    /* numbers of various subsidiary objects */
    int n_cmds;
    int n_models;
    int n_loop_models;
    int n_prints;

    char **lines;
    int *ci;
    char **eachstrs;
    MODEL **models;
    LOOP_MODEL *lmodels;
    LOOP_PRINT *prns;
    char storefile[MAXLEN];
    gretlopt storeopt;
    double **sZ;
    DATAINFO *sdinfo;
    LOOPSET *parent;
    LOOPSET **children;
    int n_children;
};

#define loop_is_progressive(l) (l->flags & LOOP_PROGRESSIVE)
#define loop_set_progressive(l) (l->flags |= LOOP_PROGRESSIVE)
#define loop_is_verbose(l) (l->flags & LOOP_VERBOSE)
#define loop_set_verbose(l) (l->flags |= LOOP_VERBOSE)
#define loop_is_quiet(l) (l->flags & LOOP_QUIET)
#define loop_set_quiet(l) (l->flags |= LOOP_QUIET)

static void gretl_loop_init (LOOPSET *loop);
static int prepare_loop_for_action (LOOPSET *loop);
static void free_loop_model (LOOP_MODEL *lmod);
static void free_loop_print (LOOP_PRINT *lprn);
static void loop_store_free (LOOPSET *loop);
static void controller_free (controller *clr);
static void print_loop_model (LOOP_MODEL *lmod, int loopnum,
			      const DATAINFO *pdinfo, PRN *prn);
static void print_loop_coeff (const DATAINFO *pdinfo, const LOOP_MODEL *lmod, 
			      int c, int n, PRN *prn);
static void print_loop_prn (LOOP_PRINT *lprn, int n,
			    const DATAINFO *pdinfo, PRN *prn);
static int save_loop_store (LOOPSET *loop, PRN *prn);
static void set_active_loop (LOOPSET *loop);

#define LOOP_BLOCK 32
#define N_LOOP_INDICES 5

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

/* --------------------------------------------------*/

static double controller_evaluate_expr (const char *expr, 
					const char *vname,
					int vnum,
					double ***pZ,
					DATAINFO *pdinfo)
{
    double x = NADBL;
    int iftest = 0;
    int v, err;

    gretl_error_clear();

    if (!strcmp(vname, "iftest")) {
	iftest = 1;
	x = generate_scalar(expr, pZ, pdinfo, &err);
    } else {
	err = generate(expr, pZ, pdinfo, OPT_P, NULL);
    }

#if LOOP_DEBUG
    fprintf(stderr, "controller_evaluate_expr: expr = '%s', genr err = %d\n", 
	    expr, err);
#endif

    if (!err && !iftest) {
	if (vnum > 0) {
	    v = vnum;
	} else {
	    v = varindex(pdinfo, vname);
	}
	if (v < pdinfo->v) {
	    x = (*pZ)[v][0];
	    if (na(x)) {
		err = 1;
	    }
	} else {
	    err = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("error evaluating loop condition"));
    }

    return x;
}

/* get a value from a loop controller element, either via a
   "genr"-type expression, or from a specified scalar variable, or
   from a numeric constant */

static double 
loop_controller_get_value (controller *clr, double ***pZ, DATAINFO *pdinfo)
{
    double ret = NADBL;

    if (clr->expr != NULL) {
	ret = controller_evaluate_expr(clr->expr, clr->vname, clr->vnum,
				       pZ, pdinfo);
	clr->val = ret;
    } else if (clr->vnum > 0) {
	ret = (*pZ)[clr->vnum][0] * clr->vsign;
	clr->val = ret;
    } else {
	ret = clr->val;
    } 

#if LOOP_DEBUG
    fprintf(stderr, "loop_controller_get_value: vnum = %d, returning %g\n", 
	    clr->vnum, ret);
#endif

    return ret;
}

static double loop_initval (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
#if LOOP_DEBUG
    fprintf(stderr, "getting loop_initval...\n");
#endif
    return loop_controller_get_value(&loop->init, pZ, pdinfo);
}

static double loop_finalval (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
#if LOOP_DEBUG
    fprintf(stderr, "getting loop_finalval...\n");
#endif
    return loop_controller_get_value(&loop->final, pZ, pdinfo);
}

static double loop_testval (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
#if LOOP_DEBUG
    fprintf(stderr, "getting loop_testval...\n");
#endif
    return loop_controller_get_value(&loop->test, pZ, pdinfo);
}

static double loop_delta (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
    double x = loop_controller_get_value(&loop->delta, pZ, pdinfo);

#if LOOP_DEBUG
    fprintf(stderr, "loop_incrval, returning %g\n", x);
#endif
    return x;
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

#define OK_LOOP_MODEL(c) (c == ARMA || c == CORC || c == GARCH || \
                          c == HCCM || c == HILU || c == HSK || \
                          c == LAD || c == OLS || c == TSLS || \
                          c == PWE || c == WLS || c == MPOLS || \
                          c == TOBIT)

/**
 * ok_in_loop:
 * @ci: command index.
 *
 * Returns: 1 if the given command is acceptable inside the loop construct,
 * 0 otherwise.
 */

int ok_in_loop (int c)
{
    if (c == APPEND ||
        c == GENR ||
	c == LOOP ||
	c == STORE ||
	c == PRINT ||
	c == PRINTF ||
	c == PVALUE ||
	c == SMPL ||
	c == IF ||
	c == ELSE ||
	c == ENDIF ||
	c == BREAK ||
	c == ENDLOOP) { 
	return 1;
    }

    /* "simple_commands" */
    if (c == ADF || 
	c == COINT || 
	c == COINT2 || 
	c == CORR ||
	c == CRITERIA || 
	c == DIFF || 
	c == HURST ||	
	c == KPSS ||
        c == LABELS ||
	c == LAGS || 
	c == LDIFF || 
	c == LOGS ||
	c == MEANTEST || 
	c == MULTIPLY || 
	c == OUTFILE ||
	c == PCA ||
	c == REMEMBER ||
        c == RENAME || 
	c == RHODIFF ||
	c == RUNS || 
        c == SET ||
	c == SETINFO ||
	c == SHELL || 
	c == SPEARMAN || 
        c == STRING ||
	c == SQUARE || 
	c == SUMMARY ||
	c == VARLIST ||
	c == VARTEST ||
	c == XTAB ) {
	return 1;
    }

    /* frequencies --- no graph allowed */
    if (c == FREQ) {
	return 1;
    }

    /* modeling commands */
    if (OK_LOOP_MODEL(c)) {
	return 1;
    }

    /* vector models */
    if (c == VAR || c == VECM) {
	return 1;
    }
    
    /* nonlinear models */
    if (c == NLS || c == MLE || c == GMM || c == END) {
	return 1;
    }

    /* basic model tests */
    if (c == ADD || c == OMIT || c == TESTUHAT) {
	return 1;
    }

    return 0;
}

static int loop_attach_child (LOOPSET *loop, LOOPSET *child)
{
    LOOPSET **children;
    int nc = loop->n_children;

    children = realloc(loop->children, (nc + 1) * sizeof *children);
    if (children == NULL) {
	return 1;
    } 

    loop->children = children;
    loop->children[nc] = child;
    child->parent = loop;
    child->level = loop->level + 1;

#if LOOP_DEBUG
    fprintf(stderr, "child loop %p has parent %p\n", 
	    (void *) child, (void *) child->parent);
#endif

    loop->n_children += 1;

    return 0;
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

    if (loop->eachstrs != NULL) {
	for (i=0; i<loop->itermax; i++) {
	    free(loop->eachstrs[i]);
	}
	free(loop->eachstrs);
    }    

    if (loop->models != NULL) {
	free(loop->models);
    } 

    if (loop->lmodels != NULL) {
	for (i=0; i<loop->n_loop_models; i++) {
#if LOOP_DEBUG
	    fprintf(stderr, "freeing loop->lmodels[%d]\n", i);
#endif
	    free_loop_model(&loop->lmodels[i]);
	}
	free(loop->lmodels);
    }

    if (loop->prns != NULL) {
	for (i=0; i<loop->n_prints; i++) { 
	    free_loop_print(&loop->prns[i]);
	}
	free(loop->prns);
    }

    loop_store_free(loop);

    if (loop->children != NULL) {
	free(loop->children);
    }

    free(loop);
}

static int 
ok_loop_var (const LOOPSET *loop, const DATAINFO *pdinfo, const char *vname)
{
    int v = varindex(pdinfo, vname);

    if (v >= pdinfo->v) {
	if (loop->parent == NULL) {
	    sprintf(gretl_errmsg, 
		    _("Undefined variable '%s' in loop condition."), vname);
	    v = LOOP_VAL_BAD;
	} else {
	    /* the variable may be defined in an outer loop,
	       and will become available at runtime? */
	    v = LOOP_VAL_UNDEF;
	}
    } else if (var_is_series(pdinfo, v)) {
	strcpy(gretl_errmsg, _("The loop control variable "
	       "must be a scalar"));
	v = LOOP_VAL_BAD;
    }
	
    return v;
}

static int 
controller_set_var (controller *clr, LOOPSET *loop, const DATAINFO *pdinfo,
		    const char *s)
{
    int v, err = 0;
    int vsign = 1;

    if (*s == '-' && isalpha((unsigned char) *(s + 1))) {
	/* negative sign in front of varname? */
	vsign = -1;
	v = ok_loop_var(loop, pdinfo, s + 1);
    } else {
	v = ok_loop_var(loop, pdinfo, s);
    }

    if (v == LOOP_VAL_BAD) {
	err = 1;
    } else if (v == LOOP_VAL_UNDEF) {
	clr->vnum = LOOP_VAL_UNDEF;
	*clr->vname = '\0';
	strncat(clr->vname, s, VNAMELEN - 1);
    } else {
	clr->vnum = v;
	clr->vsign = vsign;
    }

#if LOOP_DEBUG
    fprintf(stderr, " controller set var: vnum=%d, vsign=%d\n", 
	    clr->vnum, clr->vsign);
#endif

    return err;
}

static int parse_as_while_loop (LOOPSET *loop,
				double ***pZ, DATAINFO *pdinfo,
				const char *s)
{
    char *expr = NULL;
    double x = NADBL;
    int err = 0;

    expr = gretl_strdup(s);
    x = controller_evaluate_expr(expr, "iftest", 0, pZ, pdinfo);

    if (na(x)) {
	err = E_DATA;
	free(expr);
    }

    if (!err) {
	loop->test.expr = expr;
	strcpy(loop->test.vname, "iftest");
	loop->type = WHILE_LOOP;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_while_loop: cond = '%s', x = %g, err = %d\n", 
	    s, x, err);
#endif

    return err;
}

static const char ichars[N_LOOP_INDICES] = "ijklm"; 

static int loop_index_char_pos (int c)
{
    int i;

    for (i=0; i<N_LOOP_INDICES; i++) {
	if (c == ichars[i]) return i;
    }

    return -1;
}

static int bad_ichar (char c)
{
    int err = 0;

    if (loop_index_char_pos(c) < 0) {
	sprintf(gretl_errmsg, _("The index in a 'for' loop must be a "
				"single character in the range '%c' to '%c'"),
		ichars[0], ichars[N_LOOP_INDICES - 1]); 
	err = E_DATA;
    }

    return err;
}

#define maybe_date(s) (strchr(s, ':') || strchr(s, '/'))

static int parse_as_indexed_loop (LOOPSET *loop,
				  const double **Z,
				  const DATAINFO *pdinfo,
				  char ichar,
				  const char *lvar, 
				  const char *start,
				  const char *end)
{
    int nstart = -1, nend = -1;
    int dated = 0;
    int err = 0;

    if (lvar != NULL) {
	if (strlen(lvar) > 1) {
	    ichar = 'x';
	} else {
	    ichar = *lvar;
	}
    }

    err = bad_ichar(ichar);

    /* starting and ending values: try for dates first, then
       numeric constants, then variables */

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_indexed_loop: start='%s'\n", start);
#endif

    if (!err) {
	if (maybe_date(start)) {
	    nstart = dateton(start, pdinfo);
	    if (nstart >= 0) {
		dated = 1;
	    } else {
		err = E_DATA;
	    }
	} else if (numeric_string(start)) {
	    nstart = atoi(start);
#if LOOP_DEBUG
	    fprintf(stderr, "numeric string: nstart = %d\n", nstart);
#endif
	} else {
	    err = controller_set_var(&loop->init, loop, pdinfo, start);
	}
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_indexed_loop: end='%s'\n", end);
#endif

    if (!err) {
	if (dated) {
	    nend = dateton(end, pdinfo);
	    if (nend < 0) {
		err = E_DATA;
	    }
	} else if (numeric_string(end)) {
	    nend = atoi(end);
#if LOOP_DEBUG
	    fprintf(stderr, "numeric string: nend = %d\n", nend);
#endif
	} else {
	    err = controller_set_var(&loop->final, loop, pdinfo, end);
	}
    }

    /* if the starting and ending values are constants, check the
       range right now */

    if (!err && loop->init.vnum == 0 && 
	loop->final.vnum == 0 && nend <= nstart) {
	strcpy(gretl_errmsg, _("Ending value for loop index must be greater "
			       "than starting value."));
	err = E_DATA;
    }

    if (!err) {
	if (loop->init.vnum == 0) {
	    loop->init.val = nstart;
	} 
	if (loop->final.vnum == 0) {
	    loop->final.val = nend;
	} 
	if (dated) {
	    loop->type = DATED_LOOP;
	} else {
	    loop->type = INDEX_LOOP;
	}
	loop->ichar = ichar;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_indexed_loop: init.val=%g, final.val=%g, "
	    "final.vnum=%d\n", loop->init.val, loop->final.val,
	    loop->final.vnum);
#endif

    return err;
}

static int parse_as_count_loop (LOOPSET *loop, 
				const double **Z,
				const DATAINFO *pdinfo,
				const char *count)
{
    int nt = -1;
    int err = 0;

    /* try for a numeric value, or failing that, a variable */

    if (numeric_string(count)) {
	nt = atoi(count); 
    } else { 
	err = controller_set_var(&loop->final, loop, pdinfo, count);
    }

    if (!err && loop->final.vnum == 0 && nt <= 0) {
	strcpy(gretl_errmsg, _("Loop count must be positive."));
	err = E_DATA;
    }

    if (!err) {
	loop->init.val = 1;
	if (loop->final.vnum == 0) {
	    loop->final.val = nt;
	}
	loop->type = COUNT_LOOP;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_count_loop: init.val=%g, final.val=%g, "
	    "final.vnum=%d\n", loop->init.val, loop->final.val,
	    loop->final.vnum);
#endif

    return err;
}

/* before testing the "delta" element in a putative "for" loop, save
   the state of the variable in question so we can restore it after
   running the test */

static double *save_delta_state (int v, double **Z, DATAINFO *pdinfo)
{
    double *x0;
    int i, n;

    n = (var_is_series(pdinfo, v))? pdinfo->n : 1;

    x0 = malloc(n * sizeof *x0);
    if (x0 != NULL) {
	for (i=0; i<n; i++) {
	    x0[i] = Z[v][i];
	}
    }

    return x0;
}

static void 
restore_delta_state (int v, double *x0, double **Z, DATAINFO *pdinfo)
{
    int i, n;

    n = (var_is_series(pdinfo, v))? pdinfo->n : 1;

    for (i=0; i<n; i++) {
	Z[v][i] = x0[i];
    }

    free(x0);
}

static int 
test_forloop_element (char *s, LOOPSET *loop,
		      double ***pZ, DATAINFO *pdinfo,
		      int i)
{
    char vname[VNAMELEN] = {0};
    double x = NADBL;
    int len, err = 0;

    if (s == NULL) {
	/* FIXME: allow empty "for" fields? */
	return E_PARSE;
    }

    if (i == 1) {
	/* middle term: Boolean test */
	strcpy(vname, "iftest");
	x = controller_evaluate_expr(s, vname, 0, pZ, pdinfo);
    } else {
	/* treat as regular "genr" expression */
	len = gretl_varchar_spn(s);
	if (len < VNAMELEN) {
	    strncat(vname, s, len);
	    if (i == 2) {
		/* "increment" expression: we'll have to undo whatever
		   the test evaluation does */
		int v = varindex(pdinfo, vname);

		if (v < pdinfo->v) {
		    double *x0 = save_delta_state(v, *pZ, pdinfo);

		    if (x0 == NULL) {
			err = E_ALLOC;
		    } else {
			x = controller_evaluate_expr(s, vname, v, pZ, pdinfo);
			restore_delta_state(v, x0, *pZ, pdinfo);
		    }
		} else {
		    err = E_UNKVAR;
		}
	    } else {
		x = controller_evaluate_expr(s, vname, 0, pZ, pdinfo);
	    }	
	} else {
	    err = E_UNKVAR;
	}
    }

#if LOOP_DEBUG
    fprintf(stderr, "test_forloop_element: i=%d: '%s', x = %g\n", i, s, x);
#endif

    if (na(x)) {
	err = 1;
    }

    if (!err) {
	controller *clr = (i == 0)? &loop->init :
	    (i == 1)? &loop->test : &loop->delta;

	clr->expr = gretl_strdup(s);
	if (clr->expr == NULL) {
	    err = 1;
	} else {
	    strcpy(clr->vname, vname);
	    if (i != 1) {
		clr->vnum = varindex(pdinfo, vname);
	    }
	}	
    }	

    return err;
}

static int allocate_each_strings (LOOPSET *loop, int n)
{
    loop->eachstrs = strings_array_new(n);

    return (loop->eachstrs == NULL)? E_ALLOC : 0;
}

static int each_strings_from_named_list (LOOPSET *loop, const DATAINFO *pdinfo,
					 char *s, int *nf)
{
    int *list;
    int err = 0;

    while (isspace(*s)) s++;

    tailstrip(s);
    list = get_list_by_name(s);

    if (list == NULL) {
	err = E_UNKVAR;
    } else {
	err = allocate_each_strings(loop, list[0]);
    }

    /* when cashing out list-members, use varnames rather than numbers
    */
    if (!err) {
	int i, li;

	for (i=0; i<list[0] && !err; i++) {
	    li = list[i+1];
	    if (li < 0 || li >= pdinfo->v) {
		err = 1;
	    } else {
		loop->eachstrs[i] = gretl_strdup(pdinfo->varname[li]);
	    }
	    if (!err && loop->eachstrs[i] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (err && loop->eachstrs != NULL) {
	free_strings_array(loop->eachstrs, list[0]);
	loop->eachstrs = NULL;
    }

    if (!err) {
	*nf = list[0];
    }

    return err;
}

static int
each_strings_from_list_of_vars (LOOPSET *loop, const DATAINFO *pdinfo, 
				char *s, int *nf)
{
    char vn1[VNAMELEN], vn2[VNAMELEN];
    int v1, v2;
    int err = 0;

    *nf = 0;

    delchar(' ', s);

    if (sscanf(s, "%15[^.]..%15s", vn1, vn2) != 2) {
	err = E_PARSE;
    } else {
	v1 = varindex(pdinfo, vn1);
	v2 = varindex(pdinfo, vn2);

	if (v1 < 0 || v2 < 0 || v1 >= pdinfo->v || v2 >= pdinfo->v) {
	    err = 1;
	} else {
	    *nf = v2 - v1 + 1;
	    if (*nf <= 0) {
		err = 1;
	    }
	}

	if (!err) {
	    err = allocate_each_strings(loop, *nf);
	}

	if (!err) {
	    int i;

	    for (i=v1; i<=v2 && !err; i++) {
		loop->eachstrs[i-v1] = gretl_strdup(pdinfo->varname[i]);
		if (loop->eachstrs[i-v1] == NULL) {
		    free_strings_array(loop->eachstrs, *nf);
		    loop->eachstrs = NULL;
		    err = E_ALLOC;
		}
	    }
	}
    }
    
    return err;
}

static int
parse_as_each_loop (LOOPSET *loop, const DATAINFO *pdinfo, char *s)
{
    char ivar[8] = {0};
    char ichar = 0;
    int i, nf, err = 0;

    if (*s == '\0') {
	return E_PARSE;
    }

    /* we're looking at the string that follows "loop foreach" */

    s++; /* skip space */

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: s = '%s'\n", s);
#endif 

    /* should be something like "i " */
    if (sscanf(s, "%7s", ivar) != 1) {
	err = E_PARSE;
    } else if (strlen(ivar) > 1) {
	err = E_PARSE;
    } else {
	ichar = *ivar;
	err = bad_ichar(ichar);
    }

    if (err) {
	return err;
    }

    s += strlen(ivar);
    nf = count_fields(s);

    if (nf == 0) {
	return E_PARSE;
    }

    loop->ichar = ichar;
    
    if (nf <= 3 && strstr(s, "..") != NULL) {
	/* range of values, foo..quux */
	err = each_strings_from_list_of_vars(loop, pdinfo, s, &nf);
    } else if (nf == 1) {
	/* named list? */
	err = each_strings_from_named_list(loop, pdinfo, s, &nf);
    } else {
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
    }   

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: final.val=%g\n", loop->final.val);
#endif 

    return err;
}

static int 
parse_as_for_loop (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo, char *s)
{
    char *forstr = NULL;
    char *forbits[3];
    char *p = strchr(s, '(');
    int i, len, err = 0;

    if (p == NULL) {
	err = 1;
    } else {
	len = strcspn(p, ")");
	if (len < 4) {
	    err = E_PARSE;
	} else if ((forstr = malloc(len + 1)) == NULL) {
	    err = E_ALLOC;
	} else {
	    i = 0;
	    /* make compressed copy of string */
	    p++;
	    while (*p) {
		if (*p == ')') break;
		if (*p != ' ') {
		    forstr[i++] = *p;
		}
		p++;
	    }
	    forstr[i] = '\0';
	    /* split terms separated by ';' */
	    for (i=0; i<3 && !err; i++) {
		forbits[i] = strtok((i == 0)? forstr : NULL, ";");
		err = test_forloop_element(forbits[i], loop, 
					   pZ, pdinfo, i);
	    }
	    free(forstr);
	}
    }

    if (!err) {
	loop->type = FOR_LOOP;
    }

    return err;
}

static int parse_first_loopline (char *s, LOOPSET *loop, 
				 double ***pZ, DATAINFO *pdinfo)
{
    char lvar[VNAMELEN], rvar[VNAMELEN], op[VNAMELEN];
    char ichar;
    int err = 0;

    /* skip preliminary string */
    while (isspace(*s)) s++;
    if (!strncmp(s, "loop", 4)) {
	s += 4;
	while (isspace(*s)) s++;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_first_loopline: '%s'\n", s);
#endif

    if (sscanf(s, "%c = %15[^.]..%15s", &ichar, op, rvar) == 3) {
	err = parse_as_indexed_loop(loop, (const double **) *pZ, pdinfo,
				    ichar, NULL, op, rvar);
    } else if (sscanf(s, "for %15[^= ] = %15[^.]..%15s", lvar, op, rvar) == 3) {
	err = parse_as_indexed_loop(loop, (const double **) *pZ, pdinfo, 
				    0, lvar, op, rvar);
    } else if (!strncmp(s, "foreach", 7)) {
	err = parse_as_each_loop(loop, pdinfo, s + 7);
    } else if (!strncmp(s, "for", 3)) {
	err = parse_as_for_loop(loop, pZ, pdinfo, s + 3);
    } else if (!strncmp(s, "while", 5)) {
	err = parse_as_while_loop(loop, pZ, pdinfo, s + 6);
    } else if (sscanf(s, "%15s", lvar) == 1) {
	err = parse_as_count_loop(loop, (const double **) *pZ, pdinfo, lvar);
    } else {
	printf("parse_first_loopline: failed on '%s'\n", s);
	strcpy(gretl_errmsg, _("No valid loop condition was given."));
	err = 1;
    }

    if (!err && !na(loop->init.val) && !na(loop->final.val)) {
	int nt = loop->final.val - loop->init.val + 1;

	if (loop->type != FOR_LOOP && nt <= 0) {
	    strcpy(gretl_errmsg, _("Loop count missing or invalid"));
	    err = 1;
	}
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
	/* allocate loop->lines, loop->ci */
	*err = prepare_loop_for_action(loop);
    }

    if (*err) {
	free(loop->lines);
	free(loop->ci);
	free(loop);
	loop = NULL;
    } 

    return loop;
}

#define DEFAULT_MAX_ITER 250
#define MAX_FOR_TIMES  100000

static int get_max_iters (void)
{
    static int ml = 0;

    if (ml == 0) {
	char *mlstr = getenv("GRETL_MAX_ITER");

	if (mlstr != NULL && sscanf(mlstr, "%d", &ml)) ;
	else ml = DEFAULT_MAX_ITER;
    }

    return ml;
}

static int loop_count_too_high (LOOPSET *loop)
{
    static int max_iters = 0;
    int nt = loop->iter + 1;

    if (loop->type == FOR_LOOP) {
	/* FIXME? */
	if (nt >= MAX_FOR_TIMES) {
	    sprintf(gretl_errmsg, _("Reached maximum interations, %d"),
		    MAX_FOR_TIMES);
	    loop->err = 1;
	}
    } else {  
	if (max_iters == 0) {
	    max_iters = get_max_iters();
	}

	if (nt >= max_iters) {
	    sprintf(gretl_errmsg, _("Warning: no convergence after %d interations"),
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
 *
 * Check whether a looping condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

static int 
loop_condition (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
    int ok = 0;

    if (loop->brk) {
	/* got "break" comand */
	loop->brk = 0;
	ok = 0;
    } else if (loop->type == COUNT_LOOP || indexed_loop(loop)) {
#if LOOP_DEBUG
	fprintf(stderr, "** COUNT_LOOP: iter = %d, itermax = %d\n",
		loop->iter, loop->itermax);
#endif
	/* simple count loop or indexed loop */
	if (loop->iter < loop->itermax) {
	    ok = 1;
	}
    } else {
	/* more complex forms of control... */
	if (loop_count_too_high(loop)) {
	    /* safeguard against infinite loops */
	    ok = 0;
	} else if (loop->type == FOR_LOOP) {
	    if (loop->iter > 0) {
		loop_delta(loop, pZ, pdinfo);
	    }
	    ok = loop_testval(loop, pZ, pdinfo);
	} else if (loop->type == WHILE_LOOP) {
	    ok = loop_testval(loop, pZ, pdinfo);
	}
    }

    if (!ok) {
	if (loop->iter == 0) {
	    strcpy(gretl_errmsg, _("Loop condition not satisfied at first round"));
	    loop->err = 1;
	    loop->brk = 1;
	}
    }

    return ok;
}

static void controller_init (controller *clr)
{
    clr->val = NADBL;
    clr->vnum = 0;
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
    loop->ichar = 0;
    loop->ival = 0;
    loop->brk = 0;

    controller_init(&loop->init);
    controller_init(&loop->test);
    controller_init(&loop->delta);
    controller_init(&loop->final);

    loop->n_cmds = 0;
    loop->n_models = 0;
    loop->n_loop_models = 0;
    loop->n_prints = 0;

    loop->lines = NULL;
    loop->ci = NULL;

    loop->eachstrs = NULL;

    loop->models = NULL;
    loop->lmodels = NULL;
    loop->prns = NULL;

    loop->storefile[0] = '\0';
    loop->storeopt = OPT_NONE;
    loop->sZ = NULL;
    loop->sdinfo = NULL;

    loop->parent = NULL;
    loop->children = NULL;
    loop->n_children = 0;
}

static int prepare_loop_for_action (LOOPSET *loop)
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

static void loop_model_zero (LOOP_MODEL *lmod)
{
    int i;

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
    }
}

/**
 * loop_model_init:
 * @lmod: pointer to struct to initialize.
 * @pmod: model to take as basis.
 *
 * Initialize a #LOOP_MODEL struct, based on @pmod.
 *
 * Returns: 0 on successful completion, %E_ALLOC on error.
 */

static int loop_model_init (LOOP_MODEL *lmod, const MODEL *pmod)
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

    lmod->sum_coeff = malloc(nc * sizeof *lmod->sum_coeff);
    if (lmod->sum_coeff == NULL) return E_ALLOC;

    lmod->ssq_coeff = malloc(nc * sizeof *lmod->ssq_coeff);
    if (lmod->ssq_coeff == NULL) goto cleanup;

    lmod->sum_sderr = malloc(nc * sizeof *lmod->sum_sderr);
    if (lmod->sum_sderr == NULL) goto cleanup;

    lmod->ssq_sderr = malloc(nc * sizeof *lmod->ssq_sderr);
    if (lmod->ssq_sderr == NULL) goto cleanup;

    loop_model_zero(lmod);

#if LOOP_DEBUG
    fprintf(stderr, " model copied to %p, returning 0\n", 
	    (void *) lmod->model0);
#endif

    return 0;

 cleanup:
    free(lmod->ssq_coeff);
    free(lmod->sum_sderr);
    free(lmod->ssq_sderr);

    return E_ALLOC;
}

static void loop_print_zero (LOOP_PRINT *lprn)
{
    int i;

    for (i=0; i<lprn->list[0]; i++) { 
#ifdef ENABLE_GMP
	mpf_init(lprn->sum[i]);
	mpf_init(lprn->ssq[i]);
#else
	lprn->sum[i] = lprn->ssq[i] = 0.0;
#endif
    }
}

/**
 * loop_print_init:
 * @lprn: pointer to struct to initialize.
 * @list: list of variables to be printed.
 *
 * Initialize a #LOOP_PRINT struct.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

static int loop_print_init (LOOP_PRINT *lprn, const int *list)
{
    lprn->list = gretl_list_copy(list);
    if (lprn->list == NULL) return 1;

    lprn->sum = malloc(list[0] * sizeof *lprn->sum);
    if (lprn->sum == NULL) goto cleanup;

    lprn->ssq = malloc(list[0] * sizeof *lprn->ssq);
    if (lprn->ssq == NULL) goto cleanup;

    loop_print_zero(lprn);

    return 0;

 cleanup:
    free(lprn->list);
    free(lprn->sum);
    free(lprn->ssq);

    return 1;
}

static void loop_store_free (LOOPSET *loop)
{
    destroy_dataset(loop->sZ, loop->sdinfo);
    loop->sZ = NULL;
    loop->sdinfo = NULL;
}

/**
 * loop_store_init:
 * @loop: pointer to loop struct.
 * @fname: name of file in which to store data.
 * @list: list of variables to be stored (written to file).
 * @pdinfo: data information struct.
 * @opt: option for data storage format.
 *
 * Set up @loop for saving a set of variables.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

static int loop_store_init (LOOPSET *loop, const char *fname, 
			    const int *list, DATAINFO *pdinfo,
			    gretlopt opt)
{
    int i;

    if (loop->sZ != NULL) {
	strcpy(gretl_errmsg, "Only one 'store' command is allowed in a "
	       "progressive loop");
	return E_DATA;
    }

    loop->sdinfo = create_new_dataset(&loop->sZ, list[0] + 1, 
				      loop->itermax, 0);
    if (loop->sdinfo == NULL) {
	return E_ALLOC;
    }

    loop->storefile[0] = '\0';
    strncat(loop->storefile, fname, MAXLEN - 1);
    loop->storeopt = opt;

    for (i=1; i<loop->sdinfo->v; i++) {
	char *p;

	strcpy(loop->sdinfo->varname[i], pdinfo->varname[list[i]]);
	strcpy(VARLABEL(loop->sdinfo, i), VARLABEL(pdinfo, list[i]));
	if ((p = strstr(VARLABEL(loop->sdinfo, i), "(scalar)"))) {
	    *p = 0;
	}
    }

    return 0;
}

static int add_model_record (LOOPSET *loop)
{
    MODEL **models;
    int nm = loop->n_models;
    int err = 0;

    models = realloc(loop->models, (nm + 1) * sizeof *models);
    if (models == NULL) {
	err = E_ALLOC;
    } else {
	loop->models = models;
	loop->models[nm] = gretl_model_new();
	if (loop->models[nm] == NULL) {
	    err = E_ALLOC;
	} else {
	    loop->models[nm]->ID = nm + 1;
	}
    }

    if (!err) {
	loop->n_models += 1;
    }

    return err;
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

static int add_loop_model (LOOPSET *loop)
{
    int nlm = loop->n_loop_models;
    int err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "add_loop_model: loop->n_loop_models = %d\n",
	    loop->n_loop_models);
#endif

    loop->lmodels = realloc(loop->lmodels, (nlm + 1) * sizeof *loop->lmodels);
    if (loop->lmodels == NULL) {
	err = 1;
    } else {
	loop->n_loop_models += 1;
    }

    return err;
}

/**
 * update_loop_model:
 * @loop: pointer to loop struct.
 * @i: index number of the #LOOP_MODEL within @loop.
 * @pmod: contains estimates from the current iteration.
 *
 * Update a #LOOP_MODEL belonging to @loop, based on the results
 * in @pmod.
 */

static void update_loop_model (LOOPSET *loop, int i, MODEL *pmod)
{
    LOOP_MODEL *lmod;
    int j;
#ifdef ENABLE_GMP
    mpf_t m;

    mpf_init(m);
#endif

    lmod = &loop->lmodels[i];

#if LOOP_DEBUG
    fprintf(stderr, "update_loop_model: i=%d\n", i);
    fprintf(stderr, "pmod = %p, lmod = %p\n", (void *) pmod, (void *) lmod);
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
    }

#ifdef ENABLE_GMP
    mpf_clear(m);
#endif

#if LOOP_DEBUG
    fprintf(stderr, "update_loop_model: returning 0\n");
#endif
}

static int add_loop_print (LOOPSET *loop, const int *list)
{
    LOOP_PRINT *prns;
    int np = loop->n_prints;
    int err = 0;

    prns = realloc(loop->prns, (np + 1) * sizeof *prns);
    if (prns == NULL) {
	return 1;
    }

    loop->prns = prns;

    if (loop_print_init(&loop->prns[np], list)) { 
	fprintf(stderr, "Failed to initialize print struct for loop\n");
	err = E_ALLOC;
    }

    if (!err) {
	loop->n_prints += 1;
    }

    return err;
}

/**
 * update_loop_print:
 * @loop: pointer to loop struct.
 * @cmdnum: sequential index number of the command within @loop.
 * @list: list of variables to be printed.
 * @pZ: pointer to data matrix.
 * @pdinfo: pointer to data information struct.
 *
 * Update a #LOOP_PRINT belonging to @loop, based on the current
 * data values.
 */

static void update_loop_print (LOOPSET *loop, int i, 
			       const int *list, double ***pZ, 
			       const DATAINFO *pdinfo)
{
    LOOP_PRINT *lprn = &loop->prns[i];
    int j, t;
#ifdef ENABLE_GMP
    mpf_t m;

    mpf_init(m);
#endif
    
    for (j=1; j<=list[0]; j++) {
	if (var_is_series(pdinfo, list[j])) {
	    t = pdinfo->t1;
	} else {
	    t = 0;
	}
#ifdef ENABLE_GMP
	mpf_set_d(m, (*pZ)[list[j]][t]); 
	mpf_add(lprn->sum[j-1], lprn->sum[j-1], m);
	mpf_mul(m, m, m);
	mpf_add(lprn->ssq[j-1], lprn->ssq[j-1], m);
#else
	lprn->sum[j-1] += (*pZ)[list[j]][t];
	lprn->ssq[j-1] += (*pZ)[list[j]][t] * (*pZ)[list[j]][t];
#endif
    }

#ifdef ENABLE_GMP
    mpf_clear(m);
#endif
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
    int i, j, k;

    if (loop->type != COUNT_LOOP && !(loop_is_quiet(loop))) {
	pprintf(prn, _("\nNumber of iterations: %d\n\n"), loop->iter);
    }

    j = 0;
    k = 0;

    for (i=0; i<loop->n_cmds; i++) {
#if LOOP_DEBUG
	fprintf(stderr, "print_loop_results: loop command %d (i=%d): %s\n", 
		i+1, i, loop->lines[i]);
#endif
	if (!loop_is_progressive(loop) && loop->ci[i] == OLS) {
	    gretlopt opt;

	    strcpy(linecpy, loop->lines[i]);
	    opt = get_gretl_options(linecpy, NULL);
	    
	    if (opt & OPT_P) {
		/* deferred printing of model was requested */
		MODEL *pmod = loop->models[j++];

		set_model_id(pmod);
		printmodel(pmod, pdinfo, opt, prn);
	    }	    
	}

	if (loop_is_progressive(loop)) {
	    if (OK_LOOP_MODEL(loop->ci[i])) {
		print_loop_model(&loop->lmodels[j], 
				 loop->itermax, pdinfo, prn);
		loop_model_zero(&loop->lmodels[j]);
		j++;
	    } else if (loop->ci[i] == PRINT) {
		print_loop_prn(&loop->prns[k], 
			       loop->itermax, pdinfo, prn);
		loop_print_zero(&loop->prns[k]);
		k++;
	    } else if (loop->ci[i] == STORE) {
		save_loop_store(loop, prn);
	    }
	}
    }
}

static void free_loop_model (LOOP_MODEL *lmod)
{
#if LOOP_DEBUG
    fprintf(stderr, "free_loop_model: lmod at %p, model0 at %p\n",
	    (void *) lmod, (void *) lmod->model0);
#endif

#ifdef ENABLE_GMP
    int i;

    for (i=0; i<lmod->model0->ncoeff; i++) {
	mpf_clear(lmod->sum_coeff[i]);
	mpf_clear(lmod->sum_sderr[i]);
	mpf_clear(lmod->ssq_coeff[i]);
	mpf_clear(lmod->ssq_sderr[i]);
    }
#endif

    free(lmod->sum_coeff);
    free(lmod->sum_sderr);
    free(lmod->ssq_coeff);
    free(lmod->ssq_sderr);

    gretl_model_free(lmod->model0);
}

static void free_loop_print (LOOP_PRINT *lprn)
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
    int nc = loop->n_cmds;
    int err = 0;

#if LOOP_DEBUG
    fprintf(stderr, "real_append_line: s->line = '%s'\n", s->line);
#endif

    if ((nc + 1) % LOOP_BLOCK == 0) {
	if (add_more_loop_lines(loop)) {
	    gretl_errmsg_set(_("Out of memory!"));
	    return E_ALLOC;
	}
    }

    loop->lines[nc] = malloc(MAXLEN);
    if (loop->lines[nc] == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return E_ALLOC;
    }

    if (s->cmd->ci == PRINT && loop->type != COUNT_LOOP) {
	/* fixme: what's going on here? */
	loop->ci[nc] = 0;
    } else {
	loop->ci[nc] = s->cmd->ci;
    }

    loop->lines[nc][0] = '\0';

    if (s->cmd->opt) {
	const char *flagstr = print_flags(s->cmd->opt, s->cmd->ci);

	if (strlen(s->line) + strlen(flagstr) >= MAXLEN) {
	    strcpy(gretl_errmsg, "loop: line is too long");
	    err = 1;
	} else {
	    sprintf(loop->lines[nc], "%s%s", s->line, flagstr);
	}
    } else {
	strcpy(loop->lines[nc], s->line);
    }

    loop->n_cmds += 1;

#if LOOP_DEBUG
    fprintf(stderr, "loop %p: n_cmds = %d, line[%d] = '%s'\n",
	    (void *) loop, loop->n_cmds, nc, loop->lines[nc]);
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

static void print_loop_model (LOOP_MODEL *lmod, int loopnum,
			      const DATAINFO *pdinfo, PRN *prn)
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
    pprintf(prn, _("Statistics for %d repetitions\n"), loopnum); 
    pprintf(prn, _("Dependent variable: %s\n\n"), 
	    pdinfo->varname[lmod->model0->list[1]]);

    pputs(prn, _("                     mean of      std. dev. of     mean of"
		 "     std. dev. of\n"
		 "                    estimated      estimated"
		 "      estimated      estimated\n"
		 "      Variable     coefficients   coefficients   std. errors"
		 "    std. errors\n\n"));

    for (i=0; i<lmod->model0->ncoeff; i++) {
	print_loop_coeff(pdinfo, lmod, i, loopnum, prn);
    }

    pputc(prn, '\n');
}

static void print_loop_coeff (const DATAINFO *pdinfo, 
			      const LOOP_MODEL *lmod, 
			      int c, int n, PRN *prn)
{
#ifdef ENABLE_GMP
    mpf_t c1, c2, m, sd1, sd2;
    unsigned long ln = n;

    mpf_init(c1);
    mpf_init(c2);
    mpf_init(m);
    mpf_init(sd1);
    mpf_init(sd2);

    mpf_div_ui(c1, lmod->sum_coeff[c], ln);
    mpf_mul(m, c1, c1);
    mpf_mul_ui(m, m, ln);
    mpf_sub(m, lmod->ssq_coeff[c], m);
    mpf_div_ui(sd1, m, ln);
    if (mpf_cmp_d(sd1, 0.0) > 0) {
	mpf_sqrt(sd1, sd1);
    } else {
	mpf_set_d(sd1, 0.0);
    }

    mpf_div_ui(c2, lmod->sum_sderr[c], ln);
    mpf_mul(m, c2, c2);
    mpf_mul_ui(m, m, ln);
    mpf_sub(m, lmod->ssq_sderr[c], m);
    mpf_div_ui(sd2, m, ln);
    if (mpf_cmp_d(sd2, 0.0) > 0) {
	mpf_sqrt(sd2, sd2);
    } else {
	mpf_set_d(sd2, 0.0);
    }

    pprintf(prn, " %*s ", VNAMELEN, 
	    pdinfo->varname[lmod->model0->list[c+2]]);
    pprintf(prn, "%#14g %#14g %#14g %#14g\n", mpf_get_d(c1), mpf_get_d(sd1), 
	    mpf_get_d(c2), mpf_get_d(sd2));

    mpf_clear(c1);
    mpf_clear(c2);
    mpf_clear(m);
    mpf_clear(sd1);
    mpf_clear(sd2);
#else /* non-GMP */
    bigval m1, m2, var1, var2, sd1, sd2;
    
    m1 = lmod->sum_coeff[c] / n;
    var1 = (lmod->ssq_coeff[c] - n * m1 * m1) / n;
    sd1 = (var1 <= 0.0)? 0.0 : sqrt((double) var1);

    m2 = lmod->sum_sderr[c] / n;
    var2 = (lmod->ssq_sderr[c] - n * m2 * m2) / n;
    sd2 = (var2 <= 0.0)? 0 : sqrt((double) var2);

    pprintf(prn, " %*s ", VNAMELEN, 
	   pdinfo->varname[lmod->model0->list[c+2]]);
    pprintf(prn, "%#14g %#14g %#14g %#14g\n", (double) m1, (double) sd1, 
	    (double) m2, (double) sd2);
#endif
}

static void print_loop_prn (LOOP_PRINT *lprn, int n,
			    const DATAINFO *pdinfo, PRN *prn)
{
    int i;
    bigval mean, m, sd;

    if (lprn == NULL) return;

    pputs(prn, _("   Variable     mean         std. dev.\n"));

#ifdef ENABLE_GMP
    mpf_init(mean);
    mpf_init(m);
    mpf_init(sd);
    
    for (i=1; i<=lprn->list[0]; i++) {
	mpf_div_ui(mean, lprn->sum[i-1], (unsigned long) n);
	mpf_mul(m, mean, mean);
	mpf_mul_ui(m, m, (unsigned long) n);
	mpf_sub(sd, lprn->ssq[i-1], m);
	mpf_div_ui(sd, sd, (unsigned long) n);
	if (mpf_cmp_d(sd, 0.0) > 0) {
	    mpf_sqrt(sd, sd);
	} else {
	    mpf_set_d(sd, 0.0);
	}
	pprintf(prn, " %*s ", VNAMELEN, pdinfo->varname[lprn->list[i]]);
	pprintf(prn, "%#14g %#14g\n", mpf_get_d(mean), mpf_get_d(sd));
    }

    mpf_clear(mean);
    mpf_clear(m);
    mpf_clear(sd);
#else
    for (i=1; i<=lprn->list[0]; i++) {
	mean = lprn->sum[i-1] / n;
	m = (lprn->ssq[i-1] - n * mean * mean) / n;
	sd = (m < 0)? 0 : sqrt((double) m);
	pprintf(prn, " %*s ", VNAMELEN, pdinfo->varname[lprn->list[i]]);
	pprintf(prn, "%#14g %#14g\n", (double) mean, (double) sd);
    }
#endif
    pputc(prn, '\n');
}

static int save_loop_store (LOOPSET *loop, PRN *prn)
{
    char fname[MAXLEN];
    int *list;
    int err = 0;

    list = gretl_consecutive_list_new(1, loop->sdinfo->v - 1);
    if (list == NULL) {
	return E_ALLOC;
    }

    /* organize filename */
    if (loop->storefile[0] == '\0') {
	sprintf(fname, "%sloopdata.gdt", gretl_user_dir());	
    } else {
	strcpy(fname, loop->storefile);
    }

    if (loop->storeopt == OPT_NONE &&
	strchr(fname, '.') == NULL) {
	strcat(fname, ".gdt");
    }

    err = write_data(fname, list, (const double **) loop->sZ, 
		     loop->sdinfo, loop->storeopt, NULL);

    free(list);

    return err;
}

static int 
substitute_dollar_targ (char *str, const LOOPSET *loop,
			const double **Z, const DATAINFO *pdinfo)
{
    char targ[VNAMELEN + 3] = {0};
    char *p;
    int targlen;
    int idx = 0;
    int err = 0;

#if SUBST_DEBUG
    fprintf(stderr, "subst_dollar_targ:\n original: '%s'\n", str);
#endif

    if (loop->type == FOR_LOOP) {
	sprintf(targ, "$%s", pdinfo->varname[loop->init.vnum]);
	targlen = strlen(targ);
    } else if (indexed_loop(loop)) {
	targ[0] = '$';
	targ[1] = loop->ichar;
	targlen = 2;
	idx = loop->init.val + loop->iter;
    } else {
	return 1;
    }

#if SUBST_DEBUG
    fprintf(stderr, " target = '%s'\n", targ);
#endif

    while ((p = strstr(str, targ)) != NULL) {
	char ins[32];
	char *pins = ins;
	char *q;

	q = malloc(strlen(p));
	if (q == NULL) {
	    err = 1;
	    break;
	}

	strcpy(q, p + targlen);

	if (loop->type == FOR_LOOP) {
	    sprintf(ins, "%g", Z[loop->init.vnum][0]); /* scalar */

	    if (p - str > 0 && *(p - 1) == '[' && *(p + targlen) == ']') {
		/* got an obs-type string, on the pattern [$lvar] */
		int t = dateton(ins, pdinfo);

		if (t < 0) {
		    t = atoi(ins) - 1;
		}
		sprintf(ins, "%d", t);
	    } 
	} else if (loop->type == INDEX_LOOP) {
	    sprintf(ins, "%d", idx);
	} else if (loop->type == DATED_LOOP) {
	    ntodate(ins, idx, pdinfo);
	} else if (loop->type == EACH_LOOP) {
	    pins = loop->eachstrs[idx - 1];
	}

	strcpy(p, pins);
	strcpy(p + strlen(pins), q);
	free(q);	
    }

#if SUBST_DEBUG
    fprintf(stderr, " after: '%s'\n", str);
#endif

    return err;
}

static void update_loop_store (const int *list, LOOPSET *loop,
			       const double **Z, DATAINFO *pdinfo)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (var_is_series(pdinfo, list[i])) { 
	    loop->sZ[i][loop->iter] = Z[list[i]][pdinfo->t1 + 1]; /* ?? */
	} else {
	    loop->sZ[i][loop->iter] = Z[list[i]][0];
	}
    }
}

static int clr_attach_var (controller *clr, const DATAINFO *pdinfo)
{
    int v = varindex(pdinfo, clr->vname);
    int err = 0;

    if (v >= pdinfo->v) {
	sprintf(gretl_errmsg, 
		_("Undefined variable '%s' in loop condition."), clr->vname);
	err = 1;
    } else {
	clr->vnum = v;
    } 

    return err;
}

/* Connect up any loop control variables that were undefined at
   compile time, and flag an error if any are still undefined
   at run time.
*/

static int 
connect_loop_control_vars (LOOPSET *loop, const DATAINFO *pdinfo)
{
    int err = 0;

    if (loop->init.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->init, pdinfo);
    }

    if (!err && loop->test.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->test, pdinfo);
    }

    if (!err && loop->delta.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->delta, pdinfo);
    }

    if (!err && loop->final.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->final, pdinfo);
    }

    return err;
}

static void progressive_loop_zero (LOOPSET *loop)
{
    int i;

    for (i=0; i<loop->n_loop_models; i++) {
	free_loop_model(&loop->lmodels[i]);
    }

    loop->n_loop_models = 0;

    for (i=0; i<loop->n_prints; i++) { 
	free_loop_print(&loop->prns[i]);
    }

    loop->n_prints = 0;

    loop_store_free(loop);
}

static int 
top_of_loop (LOOPSET *loop, double ***pZ, DATAINFO *pdinfo)
{
    int err;

    loop->iter = 0;

    err = connect_loop_control_vars(loop, pdinfo);

    if (!err) {
	if (loop->init.vnum > 0) {
	    loop_initval(loop, pZ, pdinfo);
	} 

	if (loop->type == COUNT_LOOP || indexed_loop(loop)) {
	    double maxval = loop_finalval(loop, pZ, pdinfo);

	    if (na(loop->init.val) || na(maxval)) {
		err = 1;
	    } else {
		loop->itermax = maxval - loop->init.val + 1;
#if LOOP_DEBUG
		fprintf(stderr, "*** itermax = %g - %g + 1 = %d\n",
			maxval, loop->init.val, loop->itermax);
#endif
	    }
	}
    } 

    if (!err) {
	if (indexed_loop(loop)) {
	    loop->ival = loop->init.val;
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

	set_active_loop(loop);
    }

    return err;
}

static void 
print_loop_progress (const LOOPSET *loop, const DATAINFO *pdinfo,
		     PRN *prn)
{
    int i = loop->init.val + loop->iter;

    if (loop->type == INDEX_LOOP) {
	pprintf(prn, "loop: %c = %d\n\n", loop->ichar, i);
    } else if (loop->type == DATED_LOOP) {
	char obs[OBSLEN];

	ntodate(obs, i, pdinfo);
	pprintf(prn, "loop: %c = %s\n\n", loop->ichar, obs);
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
			   const double **Z, const DATAINFO *pdinfo)
{
    int err = 0;

    if (indexed_loop(loop) || loop->type == FOR_LOOP) {
	err = substitute_dollar_targ(str, loop, Z, pdinfo);
    }

    while (!err && (loop = subst_loop_in_parentage(loop)) != NULL) {
	err = substitute_dollar_targ(str, loop, Z, pdinfo);
    }

    return err;
}

/* get the next command for a loop by pulling a line off the
   stack of loop commands.
*/

static int next_command (char *targ, LOOPSET *loop, int *j)
{
    int ret = 1;

    if (*j < loop->n_cmds) {
	if (!isprint(loop->lines[*j][0])) {
	    fprintf(stderr, "bad line %d from loop at %p\n", *j, (void *) loop);
	    exit(EXIT_FAILURE);
	}
	strcpy(targ, loop->lines[*j]);
	*j += 1;
    } else {
	ret = 0;
    }

    return ret;
}

int gretl_loop_exec (ExecState *s, double ***pZ, DATAINFO **ppdinfo) 
{
    LOOPSET *loop = currloop;
    DATAINFO *pdinfo = *ppdinfo;
    MODEL **models = s->models;
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    char errline[MAXLINE];
    int mod_id = 0;
    int err = 0;

    /* for the benefit of the caller: register the fact that execution
       of this loop is now under way */
    loop_execute = 0;

    if (loop == NULL) {
	pputs(prn, "Got a NULL loop\n");
	return 1;
    }

    if (loop->n_cmds == 0) {
	pputs(prn, _("No commands in loop\n"));
	return 0;
    }

    set_loop_on(); /* libset.c */

#if LOOP_DEBUG
    fprintf(stderr, "loop_exec: loop = %p\n", (void *) loop);
#endif

    err = top_of_loop(loop, pZ, pdinfo);
    
    while (!err && loop_condition(loop, pZ, pdinfo)) {
	int modnum = 0;
	int lmodnum = 0;
	int printnum = 0;
	int childnum = 0;
	int j;

#if LOOP_DEBUG
	fprintf(stderr, "top of loop: iter = %d\n", loop->iter);
#endif

	if (gretl_echo_on() && indexed_loop(loop) &&
	    !loop_is_quiet(loop)) {
	    print_loop_progress(loop, pdinfo, prn);
	}

	j = 0;
	while (!err && next_command(line, loop, &j)) {

	    pdinfo = *ppdinfo;
	    set_active_loop(loop);

#if LOOP_DEBUG
	    fprintf(stderr, " j=%d, line='%s'\n", j, line);
#endif
	    strcpy(errline, line);

	    err = make_dollar_substitutions(line, loop, 
					    (const double **) *pZ,
					    pdinfo);
	    if (err) {
		break;
	    }

	    /* We already have the "ci" index recorded, but here
	       we do some further parsing. 
	    */
	    err = parse_command_line(line, cmd, pZ, pdinfo);

	    if (cmd->ci < 0) {
		continue;
	    } else if (err) {
		break;
	    }

	    if (gretl_echo_on() && indexed_loop(loop)) {
		if (s->cmd->ci == ENDLOOP) {
		    pputc(prn, '\n');
		} else if (!loop_is_quiet(loop)) {
		    echo_cmd(cmd, pdinfo, line, 0, prn);
		}
	    }

	    switch (cmd->ci) {

	    case LOOP:
		if (childnum < loop->n_children) {
		    currloop = loop->children[childnum++];
		    err = gretl_loop_exec(s, pZ, ppdinfo);
		    pdinfo = *ppdinfo;
		} else {
		    fprintf(stderr, "Got a LOOP command, don't know what to do!\n");
		    err = 1;
		}
		break;

	    case BREAK:
		loop->brk = 1;
		break;

	    case ENDLOOP:
		break;

	    case FREQ:
		/* note: no graphs in loops */
		err = freqdist(cmd->list[1], (const double **) *pZ, pdinfo, 0, 
			       cmd->opt, prn);
		break;

	    case ARMA:
	    case CORC:
	    case GARCH:
	    case HCCM:
	    case HILU:
	    case HSK:
	    case LAD:
	    case MPOLS:
	    case OLS:
	    case TOBIT:
	    case TSLS:
	    case PWE:
	    case WLS:
		/* if this is the first time round, allocate space
		   for each loop model */
		if (loop->iter == 0) {
		    if (loop_is_progressive(loop)) {
			err = add_loop_model(loop);
		    } else if (cmd->opt & OPT_P) {
			err = add_model_record(loop);
		    }
		    if (err) {
			break;
		    }
		} 

		/* estimate the model called for, using models[0] */
		clear_model(models[0]);

		if (cmd->ci == OLS || cmd->ci == WLS || cmd->ci == HCCM) {
		    *models[0] = lsq(cmd->list, pZ, pdinfo, cmd->ci, cmd->opt);
		} else if (cmd->ci == LAD) {
		    *models[0] = lad(cmd->list, pZ, pdinfo);
		} else if (cmd->ci == MPOLS) {
		    *models[0] = mp_ols(cmd->list, (const double **) *pZ, pdinfo);
		} else if (cmd->ci == HSK) {
		    *models[0] = hsk_func(cmd->list, pZ, pdinfo);
		} else if (cmd->ci == ARMA) {
		    *models[0] = arma(cmd->list, (const double **) *pZ, pdinfo, 
				      cmd->opt, prn);
		} else if (cmd->ci == TOBIT) {
		    *models[0] = tobit_model(cmd->list, pZ, pdinfo,
					     (cmd->opt & OPT_V)? prn : NULL);
		} else if (cmd->ci == TSLS) {
		    *models[0] = tsls_func(cmd->list, TSLS, pZ, pdinfo, cmd->opt);
		} else if (cmd->ci == GARCH) {
		    *models[0] = garch(cmd->list, pZ, pdinfo, cmd->opt, prn);
		} else if (cmd->ci == CORC || cmd->ci == HILU || cmd->ci == PWE) {
		    double rho = estimate_rho(cmd->list, pZ, pdinfo, cmd->ci,
					      &err, cmd->opt, prn);
		    if (err) {
			break;
		    }
		    *models[0] = ar1_lsq(cmd->list, pZ, pdinfo, cmd->ci, cmd->opt, rho);
		}

		if ((err = models[0]->errcode)) {
		    break;
		} 

		if (loop_is_progressive(loop)) {
		    int m = lmodnum++;

		    if (loop->iter == 0) {
			err = loop_model_init(&loop->lmodels[m], models[0]);
			if (err) {
			    gretl_errmsg_set(_("Failed to initialize model for loop\n"));
			    break;
			}
		    }
		    update_loop_model(loop, m, models[0]);
		    set_as_last_model(models[0], GRETL_OBJ_EQN);
		} else if (cmd->opt & OPT_P) {
		    /* deferred printing of model results */
		    int m = modnum++;

		    swap_models(models[0], loop->models[m]);
		    loop->models[m]->ID = j;
		    set_as_last_model(loop->models[m], GRETL_OBJ_EQN);
		    model_count_minus();
		} else {
		    models[0]->ID = ++mod_id;
		    printmodel(models[0], pdinfo, cmd->opt, prn);
		    set_as_last_model(models[0], GRETL_OBJ_EQN);
		}
		break;

	    case ADD:
	    case OMIT:
	    case GMM:
	    case MLE:
	    case NLS:
		if (loop_is_progressive(loop)) {
		    err = 1;
		} else {
		    goto cmd_exec;
		}
		break;

	    case PRINT:
		if (cmd->param[0] != '\0') {
		    goto cmd_exec;
		} else if (loop_is_progressive(loop)) {
		    int p = printnum++;

		    if (loop->iter == 0) {
			err = add_loop_print(loop, cmd->list);
			if (err) {
			    break;
			}
		    }
		    update_loop_print(loop, p, cmd->list, pZ, pdinfo);
		} else {
		    err = printdata(cmd->list, cmd->extra,
				    (const double **) *pZ, pdinfo, 
				    cmd->opt, prn);
		}
		break;

	    case STORE:
		if (loop_is_progressive(loop)) {
		    if (loop->iter == 0) {
			err = loop_store_init(loop, cmd->param, cmd->list, 
					      pdinfo, cmd->opt);
			if (err) {
			    break;
			}
		    }
		    update_loop_store(cmd->list, loop, (const double **) *pZ, 
				      pdinfo);
		} else {
		    goto cmd_exec;
		}
		break;

	    default: 
	    cmd_exec:
		if (cmd->ci == GENR && !loop_is_verbose(loop)) {
		    cmd->opt |= OPT_Q;
		}
		err = gretl_cmd_exec(s, pZ, ppdinfo);
		pdinfo = *ppdinfo;
		break;

	    } /* end switch on specific command number */

	} /* end execution of commands within loop */

	if (err && get_halt_on_error() == 0) {
	    errmsg(err, prn);
	    err = 0;
	}

	loop->iter += 1;

	if (indexed_loop(loop)) {
	    loop->ival += 1;
	}

    } /* end iterations of loop */

    if (err) {
	errmsg(err, prn);
	pprintf(prn, ">> %s\n", errline);
    } else if (loop->err) {
	errmsg(loop->err, prn);
	err = loop->err;
    }

    if (!err && loop->iter > 0) {
	print_loop_results(loop, pdinfo, prn); 
    }

    if (loop->n_models > 0) {
	/* need to update models[0] */
	GretlObjType type;
	void *ptr = get_last_model(&type);
	int i;

	if (type == GRETL_OBJ_EQN && models[0] != ptr) {
	    swap_models(models[0], loop->models[loop->n_models - 1]);
	    set_as_last_model(models[0], GRETL_OBJ_EQN);
	}
	for (i=0; i<loop->n_models; i++) {
	    gretl_model_free(loop->models[i]);
	}
    }

    if (line != NULL) {
	*line = '\0';
    } 

    set_active_loop(loop->parent);
    set_loop_off();

    if (loop->parent == NULL) {
	/* reached top of stack: clean up */
	gretl_loop_destroy(loop);
	currloop = NULL;
    }

    if (get_halt_on_error()) {
	return err;
    } else {
	return 0;
    }
}

static const LOOPSET *
parent_with_ichar (const LOOPSET *loop, int c)
{
    while ((loop = loop->parent) != NULL) {
	if (indexed_loop(loop) && loop->ichar == c) {
	    return loop;
	}
    }

    return NULL;
}

static const LOOPSET *
indexed_loop_set_or_get (LOOPSET *loop, int set, int c)
{
    static LOOPSET *active_loop;
    const LOOPSET *ret = NULL;

    if (set) {
	active_loop = loop;
    } else if (active_loop != NULL) {
	if (indexed_loop(active_loop) && active_loop->ichar == c) {
	    ret = active_loop;
	} else {
	    ret = parent_with_ichar(active_loop, c);
	}
    }

#if IDX_DEBUG > 1
    if (set) {
	fprintf(stderr, "indexed_loop_record: set active_loop = %p\n",
		(void *) loop);
    } else {
	fprintf(stderr, "indexed_loop_record: returning %p for c='%c'\n",
		(void *) ret, c);
    }	
#endif

    return ret;
}

static void set_active_loop (LOOPSET *loop)
{
    indexed_loop_set_or_get(loop, 1, 0);
}

int is_active_index_loop_char (int c)
{
    return indexed_loop_set_or_get(NULL, 0, c) != NULL;
}

int loop_scalar_read (int c)
{
    const LOOPSET *loop = indexed_loop_set_or_get(NULL, 0, c);

    return (loop == NULL)? -1 : loop->ival;
}

#define IFDEBUG 0

/* if-then stuff - conditional execution */

int if_eval (const char *line, double ***pZ, DATAINFO *pdinfo)
{
    double val = NADBL;
    int err, ret = -1;

#if IFDEBUG
    fprintf(stderr, "if_eval: line = '%s'\n", line);
#endif

    if (!strncmp(line, "if", 2)) {
	line += 2;
    } else if (!strncmp(line, "elif", 4)) {
	/* dead code, for now */
	line += 4;
    }

    while (*line == ' ') line++;

    val = generate_scalar(line, pZ, pdinfo, &err);

#if IFDEBUG
    if (err) {
	fprintf(stderr, "if_eval: generate returned %d\n", err);
    }
#endif

    if (err) {
	strcpy(gretl_errmsg, _("error evaluating 'if'"));
    } else if (na(val)) {
	strcpy(gretl_errmsg, _("indeterminate condition for 'if'"));
    } else {
	ret = (int) val;
    }

#if IFDEBUG
    fprintf(stderr, "if_eval: returning %d\n", ret);
#endif

    return ret;
}

#if IFDEBUG
static const char *ifstr (int c)
{
    if (c == SET_FALSE) return "SET_FALSE";
    if (c == SET_TRUE)  return "SET_TRUE";
    if (c == SET_ELSE)  return "SET_ELSE";
    if (c == SET_ELIF)  return "SET_ELIF";
    if (c == SET_ENDIF) return "SET_ENDIF";
    if (c == IS_FALSE)  return "IS_FALSE";
    if (c == DOINDENT)  return "DOINDENT";
    if (c == UNINDENT)  return "UNINDENT";
    if (c == GETINDENT) return "GETINDENT";
    if (c == RELAX)     return "RELAX";
    return "UNKNOWN";
}
#endif

#define IF_DEPTH 32

int ifstate (int code)
{
    static unsigned char T[IF_DEPTH];
    static unsigned char got_if[IF_DEPTH];
    static unsigned char got_else[IF_DEPTH];
    static unsigned char indent;
    int i, ret = 0;

#if IFDEBUG
    fprintf(stderr, "ifstate: code = %s\n", ifstr(code));
#endif

    if (code == RELAX) {
	indent = 0;
    } else if (code == DOINDENT) {
	return ++indent;
    } else if (code == UNINDENT) {
	return --indent;
    } else if (code == GETINDENT) {
	return indent;
    } else if (code == SET_FALSE || code == SET_TRUE) {
	indent++;
	if (indent >= IF_DEPTH) {
	    fprintf(stderr, "IF depth (%d) exceeded\n", IF_DEPTH);
	    ret = 1;
	} else {
	    T[indent] = (code == SET_TRUE);
	    got_if[indent] = 1;
	    got_else[indent] = 0;
	}
    } else if (code == SET_ELSE) {
	if (got_else[indent] || !got_if[indent]) {
	    strcpy(gretl_errmsg, "Unmatched \"else\"");
	    ret = 1; 
	} else {
	    T[indent] = !T[indent];
	    got_else[indent] = 1;
	}
    } else if (code == SET_ELIF) {
	if (got_else[indent] || !got_if[indent]) {
	    strcpy(gretl_errmsg, "Unmatched \"elif\"");
	    ret = 1; 
	} else {
	    T[indent] = !T[indent];
	}
    } else if (code == SET_ENDIF) {
	if (!got_if[indent] || indent == 0) {
	    strcpy(gretl_errmsg, "Unmatched \"endif\"");
	    ret = 1; 
	} else {
	    got_if[indent] = 0;
	    got_else[indent] = 0;
	    indent--;
	}
    } else if (code == IS_FALSE) {
	for (i=1; i<=indent; i++) {
	    if (T[i] == 0) {
		ret = 1;
		break;
	    }
	}
    }

#if IFDEBUG
    fprintf(stderr, "ifstate: returning %d (indent = %d)\n", ret, indent);
#endif

    return ret;
}
