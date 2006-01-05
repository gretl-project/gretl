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

/*  monte_carlo.c - loop simulation procedures
*/  

#include "libgretl.h" 
#include "loop_private.h"
#include "libset.h"
#include "compat.h"
#include "cmd_private.h"
#include "var.h"
#include "objstack.h"

#include <time.h>
#include <unistd.h>

#define LOOP_DEBUG 0

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

enum loop_scalar_codes {
    LOOP_SCALAR_READ,
    LOOP_SCALAR_PUT,
    LOOP_SCALAR_INCR
};

#define indexed_loop(l) (l->type == INDEX_LOOP || \
                         l->type == DATED_LOOP || \
			 l->type == EACH_LOOP)

typedef struct {
    int ID;
    int *list;
    bigval *sum;
    bigval *ssq;
} LOOP_PRINT;  

/* below: used for special "progressive" loop */ 

typedef struct {
    int ID;                 /* ID number for model */
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
};

typedef struct controller_ controller;

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
    char ineq;
    char brk;

    controller left;
    controller right;
    controller init;
    controller incr;

    int incrsgn;

    /* numbers of various subsidiary objects */
    int ncmds;
    int nmod;
    int nprn;
    int nstore;

    int next_model;
    int next_print;
    char **lines;
    int *ci;
    char **eachstrs;
    MODEL **models;
    LOOP_MODEL *lmodels;
    LOOP_PRINT *prns;
    char storefile[MAXLEN];
    char **storename;
    char **storelbl;
    double *storeval;
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
static void print_loop_model (LOOP_MODEL *lmod, int loopnum,
			      const DATAINFO *pdinfo, PRN *prn);
static void print_loop_coeff (const DATAINFO *pdinfo, const LOOP_MODEL *lmod, 
			      int c, int n, PRN *prn);
static void print_loop_prn (LOOP_PRINT *lprn, int n,
			    const DATAINFO *pdinfo, PRN *prn);
static int print_loop_store (LOOPSET *loop, PRN *prn);
static int get_prnnum_by_id (LOOPSET *loop, int id);
static int get_modnum_by_cmdnum (LOOPSET *loop, int cmdnum);
static void set_active_loop (LOOPSET *loop);
static int loop_scalar_index (int c, int opt, int put);

#define LOOP_BLOCK 32
#define N_LOOP_INDICES 5

static double 
loop_controller_get_value (controller *clr, const double **Z)
{
    double ret = NADBL;

    if (clr->vnum > 0) {
	ret = Z[clr->vnum][0] * clr->vsign;
    } else {
	ret = clr->val;
    } 

#if LOOP_DEBUG
    fprintf(stderr, "loop_controller_get_value: vnum = %d, returning %g\n", 
	    clr->vnum, ret);
#endif

    return ret;
}

static double loop_lval (LOOPSET *loop, const double **Z)
{
#if LOOP_DEBUG
    fprintf(stderr, "getting loop_lval...\n");
#endif
    return loop_controller_get_value(&loop->left, Z);
}

static double loop_rval (LOOPSET *loop, const double **Z)
{
#if LOOP_DEBUG
    fprintf(stderr, "getting loop_rval...\n");
#endif
    return loop_controller_get_value(&loop->right, Z);
}

static double loop_incrval (LOOPSET *loop, const double **Z)
{
    double x = loop_controller_get_value(&loop->incr, Z);

    if (!na(x)) {
	x *= loop->incrsgn;
    }

#if LOOP_DEBUG
    fprintf(stderr, "loop_incrval, returning %g\n", x);
#endif

    return x;
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

/**
 * ok_in_loop:
 * @ci: command index.
 * @loop: pointer to loop structure
 *
 * Returns: 1 if the given command is acceptable inside the loop construct,
 * 0 otherwise.
 */

int ok_in_loop (int ci, const LOOPSET *loop)
{
    int ok = 0;

    if (ci == GENR ||
	ci == LOOP ||
	ci == STORE ||
	ci == PRINT ||
	ci == PRINTF ||
	ci == PVALUE ||
	ci == SMPL ||
	ci == IF ||
	ci == ELSE ||
	ci == ENDIF ||
	ci == BREAK ||
	ci == ENDLOOP) { 
	ok = 1;
    }

    /* "simple_commands" */
    else if (ci == ADF || 
	     ci == COINT || 
	     ci == COINT2 || 
	     ci == CORR ||
	     ci == CRITERIA || 
	     ci == CRITICAL || 
	     ci == DIFF || 
	     ci == HURST ||	
	     ci == KPSS ||
	     ci == LABEL ||
	     ci == LAGS || 
	     ci == LDIFF || 
	     ci == LOGS ||
	     ci == MEANTEST || 
	     ci == MULTIPLY || 
	     ci == OUTFILE ||
	     ci == PCA ||
	     ci == RHODIFF ||
	     ci == RUNS || 
	     ci == SIM ||
	     ci == SPEARMAN || 
	     ci == SQUARE || 
	     ci == SUMMARY ||
	     ci == VARTEST) {
	ok = 1;
    }

    /* modeling commands */
    else if (ci == ARMA ||
	     ci == CORC ||
	     ci == GARCH ||
	     ci == HCCM ||
	     ci == HILU ||
	     ci == HSK ||
	     ci == LAD ||
	     ci == OLS ||
	     ci == PWE ||
	     ci == WLS) {
	return 1;
    }

    /* vector models */
    else if (ci == VAR || ci == VECM) {
	ok = 1;
    }
    
    /* nonlinear models */
    else if (ci == NLS || ci == MLE || ci == END) {
	ok = 1;
    }

    /* basic model tests */
    else if (ci == ADD || ci == OMIT) {
	ok = 1;
    }

    return ok;
}

static int loop_attach_child (LOOPSET *loop, LOOPSET *child)
{
    LOOPSET **children;
    int nc = loop->n_children + 1;

    children = realloc(loop->children, nc * sizeof *loop->children);
    if (children == NULL) {
	return 1;
    } 

    loop->children = children;
    loop->children[nc - 1] = child;
    child->parent = loop;
    child->level = loop->level + 1;

#if LOOP_DEBUG
    fprintf(stderr, "child loop %p has parent %p\n", 
	    (void *) child, (void *) child->parent);
#endif

    loop->n_children += 1;

    return 0;
}

static LOOPSET *gretl_loop_new (LOOPSET *parent, int loopstack)
{
    LOOPSET *loop;

    loop = malloc(sizeof *loop);
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

int opstr_to_op (const char *s)
{
    int op = LOOP_VAL_BAD;

    if (strstr(s, ">=")) {
	op = OP_GTE;
    } else if (strstr(s, ">")) {
	op = OP_GT;
    } else if (strstr(s, "<=")) {
	op = OP_LTE;
    } else if (strstr(s, "<")) {
	op = OP_LT;
    } 
    
    return op;
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
    } else if (v == 0 || pdinfo->vector[v]) {
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
				const DATAINFO *pdinfo,
				char *lvar, char *rvar, 
				char *opstr)
{
    int err = 0;

    loop->ineq = opstr_to_op(opstr);

    if (loop->ineq == LOOP_VAL_BAD) {
	err = 1;
    }

    if (!err) {
	err = controller_set_var(&loop->left, loop, pdinfo, lvar);
    }

    if (!err) {
	if (numeric_string(rvar)) {
	    loop->right.val = dot_atof(rvar);
	} else {
	    err = controller_set_var(&loop->right, loop, pdinfo, rvar);
	}
    }

    if (!err) {
	loop->type = WHILE_LOOP;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_while_loop: left.var=%d, right.var=%d, right.val=%g\n", 
	    loop->left.vnum, loop->right.vnum, loop->right.val);
    fprintf(stderr, " (lvar='%s', rvar='%s')\n", lvar, rvar);
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
	err = 1;
    }

    return err;
}

#define maybe_date(s) (strchr(s, ':') || strchr(s, '/'))

static int parse_as_indexed_loop (LOOPSET *loop,
				  const DATAINFO *pdinfo,
				  const double **Z,
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
		err = 1;
	    }
	} else if (numeric_string(start)) {
	    nstart = atoi(start);
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
		err = 1;
	    }
	} else if (numeric_string(end)) {
	    nend = atoi(end);
	} else {
	    err = controller_set_var(&loop->right, loop, pdinfo, end);
	}
    }

    /* if the starting and ending values are constants, check the
       range right now */

    if (!err && loop->init.vnum == 0 && 
	loop->right.vnum == 0 && nend <= nstart) {
	strcpy(gretl_errmsg, _("Ending value for loop index must be greater "
			       "than starting value."));
	err = 1;
    }

    if (!err) {
	if (loop->init.vnum == 0) {
	    loop->init.val = nstart;
	} 

	if (loop->right.vnum == 0) {
	    loop->right.val = nend;
	} 

	if (dated) {
	    loop->type = DATED_LOOP;
	} else {
	    loop->type = INDEX_LOOP;
	}

	loop->ichar = ichar;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_indexed_loop: init.val=%g, right.val=%g, "
	    "right.vnum=%d\n", loop->init.val, loop->right.val,
	    loop->right.vnum);
#endif

    return err;
}

static int parse_as_count_loop (LOOPSET *loop, 
				const DATAINFO *pdinfo,
				const double **Z,
				const char *count)
{
    int nt = -1;
    int err = 0;

    /* try for a numeric value, or failing that, a variable */

    if (numeric_string(count)) {
	nt = atoi(count); 
    } else { 
	err = controller_set_var(&loop->right, loop, pdinfo, count);
    }

    if (!err && loop->right.vnum == 0 && nt <= 0) {
	strcpy(gretl_errmsg, _("Loop count must be positive."));
	err = 1;
    }

    if (!err) {
	loop->init.val = 1;
	if (loop->right.vnum == 0) {
	    loop->right.val = nt;
	}
	loop->type = COUNT_LOOP;
    }

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_count_loop: init.val=%g, right.val=%g, "
	    "right.vnum=%d\n", loop->init.val, loop->right.val,
	    loop->right.vnum);
#endif

    return err;
}

static int 
test_forloop_element (const char *s, LOOPSET *loop,
		      DATAINFO *pdinfo, double ***pZ,
		      int i)
{
    char lhs[9], opstr[3], rhs[9];
    int ngot, err = 0;

    if (s == NULL) return 1;

    if (i == 0) {
	ngot = sscanf(s, "%8[^=]=%8s", lhs, rhs) + 1;
	strcpy(opstr, "=");
    } else {
	ngot = sscanf(s, "%8[^-+*/=<>]%2[-+*/=<>]%8[^-+*/=<>]", 
		      lhs, opstr, rhs);
    }

#if LOOP_DEBUG
    fprintf(stderr, "read forloop element, i=%d: '%s'\n", i, s);
    fprintf(stderr, " got lhs='%s', opstr='%s', rhs='%s'\n",
	    lhs, opstr, rhs);
#endif

    if (ngot != 3) {
	err = E_PARSE;
    } else {
	int v = varindex(pdinfo, lhs);

#if LOOP_DEBUG
	fprintf(stderr, " lhs: varindex = %d (pdinfo->v = %d)\n", v, pdinfo->v);
#endif
	/* examine the LHS */
	if (i == 0) {
	    if (v == pdinfo->v) {
		err = dataset_add_scalar(pZ, pdinfo);
		if (err) {
		    strcpy(gretl_errmsg, _("Out of memory!"));
		} else {
		    strcpy(pdinfo->varname[v], lhs);
		    (*pZ)[v][0] = 0.0;
		    loop->left.vnum = v;
		}
	    } else {
		loop->left.vnum = ok_loop_var(loop, pdinfo, lhs);
		if (loop->left.vnum == LOOP_VAL_BAD) {
		    err = 1;
		} 
	    }
	} else if (v != loop->left.vnum) {
	    /* the LHS var must be the same in all three "for" fields */
	    printf("error in test_forloop_element: i=%d, lhs='%s', v=%d\n", 
                   i, lhs, v);
	    strcpy(gretl_errmsg, _("No valid loop condition was given."));
	    err = 1;
	}
	    
	if (!err) {
	    /* examine the RHS */
	    controller *clr;

	    clr = (i == 0)? &loop->init :
		(i == 1)? &loop->right : &loop->incr;

	    if (numeric_string(rhs)) {
		clr->val = dot_atof(rhs);
	    } else {
		err = controller_set_var(clr, loop, pdinfo, rhs);
	    } 
	}
	
	/* examine operator(s) */
	if (!err) {
	    if (i == 1) {
		loop->ineq = opstr_to_op(opstr);
		if (loop->ineq == LOOP_VAL_BAD) {
		    err = 1;
		}
	    } else if (i == 2) {
		if (!strcmp(opstr, "-=")) {
		    loop->incrsgn = -1;
		} else if (strcmp(opstr, "+=")) {
		    sprintf(gretl_errmsg, "Invalid operator '%s'", opstr);
		    err = 1;
		}
	    }
	}
    }

    return err;
}

static void destroy_each_strings (LOOPSET *loop, int n)
{
    int i;

    for (i=0; i<n; i++) {
	free(loop->eachstrs[i]);
    }
    
    free(loop->eachstrs);
    loop->eachstrs = NULL;
}

static int allocate_each_strings (LOOPSET *loop, int n)
{
    int i, err = 0;

    loop->eachstrs = malloc(n * sizeof *loop->eachstrs);

    if (loop->eachstrs == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<n; i++) {
	    loop->eachstrs[i] = NULL;
	}
    }

    return err;
}

static int each_strings_from_named_list (LOOPSET *loop, const DATAINFO *pdinfo,
					 char *s, int *nf)
{
    int *list;
    int err = 0;

    chopstr(s);
    list = get_list_by_name(s);
    
    if (list == NULL) {
	err = 1;
    } else {
	err = allocate_each_strings(loop, list[0]);
    }

#if 1 
    /* when cashing out list-members, use varnames rather than numbers */
    if (!err) {
	int i, li;

	for (i=1; i<=list[0] && !err; i++) {
	    li = list[i];
	    if (li < 0 || li >= pdinfo->v) {
		err = 1;
	    } else {
		loop->eachstrs[i-1] = gretl_strdup(pdinfo->varname[li]);
		if (loop->eachstrs[i-1] == NULL) {
		    err = 1;
		}
	    }
	}
    }
#else
    if (!err) {
	char numstr[16];
	int i, li;

	for (i=1; i<=list[0] && !err; i++) {
	    li = list[i];
	    if (abs(li) > 9999999) {
		err = 1;
	    } else {
		sprintf(numstr, "%d", li);
		loop->eachstrs[i-1] = gretl_strdup(numstr);
		if (loop->eachstrs[i-1] == NULL) {
		    err = 1;
		}
	    }
	}
    }
#endif

    if (err && loop->eachstrs != NULL) {
	destroy_each_strings(loop, list[0]);
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

    if (sscanf(s, "%8[^.]..%8s", vn1, vn2) != 2) {
	err = 1;
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
		    destroy_each_strings(loop, *nf);
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
	return 1;
    }

    /* we're looking at the string that follows "loop foreach" */

    s++;

    /* should be something like "i " */
    if (sscanf(s, "%7s", ivar) != 1) {
	err = 1;
    } else if (strlen(ivar) > 1) {
	err = 1;
    } else {
	ichar = *ivar;
	err = bad_ichar(ichar);
    }

    if (err) {
	return 1;
    }

    s += strlen(ivar);
    nf = count_fields(s);

    if (nf == 0) {
	return 1;
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
		destroy_each_strings(loop, nf);
		err = E_ALLOC;
	    } else {
		s += len;
	    }
	}
    }

    if (!err) {
	loop->type = EACH_LOOP;
	loop->init.val = 1;
	loop->right.val = nf;
    }   

#if LOOP_DEBUG
    fprintf(stderr, "parse_as_each_loop: right.val=%g\n", loop->right.val);
#endif 

    return err;
}

static int 
parse_as_for_loop (LOOPSET *loop,
		   DATAINFO *pdinfo, double ***pZ,
		   char *s)
{
    char *p = strchr(s, '(');
    int err = 0;

    if (p == NULL) {
	err = 1;
    } else {
	char *forstr = NULL;
	int len = strcspn(p, ")");

	if (len < 4 || (forstr = malloc(len)) == NULL) {
	    err = 1;
	} else {
	    char *forbits[3];
	    int i = 0;

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
					   pdinfo, pZ, i);
	    }

	    free(forstr);
	}
    }

    if (!err) {
	loop->type = FOR_LOOP;
    }

    return err;
}

/**
 * parse_loopline:
 * @line: command line.
 * @ploop: current loop struct pointer, or %NULL.
 * @loopstack: stacking level for the loop.
 * @pdinfo: data information struct.
 * @pZ: pointer to data array.
 *
 * Parse a line specifying a loop condition.
 *
 * Returns: loop pointer on successful completion, %NULL on error.
 */

static LOOPSET *
parse_loopline (char *line, LOOPSET *ploop, int loopstack,
		DATAINFO *pdinfo, double ***pZ)
{
    LOOPSET *loop;
    char lvar[VNAMELEN], rvar[VNAMELEN], op[VNAMELEN];
    gretlopt opt = OPT_NONE;
    char ichar;
    int err = 0;

    *gretl_errmsg = '\0';

#if LOOP_DEBUG
    fprintf(stderr, "parse_loopline: ploop = %p, loopstack = %d\n"
	    " line: '%s'\n", (void *) ploop, loopstack, line);
#endif

    while (isspace((unsigned char) *line)) {
	line++;
    }

    if (strncmp(line, "loop", 4)) {
	printf("parse_loopline: line didn't begin with 'loop': '%s'\n", line);
	strcpy(gretl_errmsg, _("No valid loop condition was given."));
	return NULL;
    }

    if (ploop == NULL) {
	/* starting from scratch */
#if LOOP_DEBUG
	fprintf(stderr, "parse_loopline: starting from scratch\n");
#endif
	loop = gretl_loop_new(NULL, 0);
	if (loop == NULL) {
	    gretl_errmsg_set(_("Out of memory!"));
	    return NULL;
	}
    } else if (loopstack > ploop->level) {
	/* have to nest this loop */
#if LOOP_DEBUG
	fprintf(stderr, "parse_loopline: adding child\n");
#endif
	loop = gretl_loop_new(ploop, loopstack);
	if (loop == NULL) {
	    gretl_errmsg_set(_("Out of memory!"));
	    return NULL;
	} else {
	    opt = get_gretl_options(line, &err);
	    if (err) {
		gretl_loop_destroy(loop);
		return NULL;
	    } 
	}
    } else {
	/* shouldn't happen: need error message? */
	loop = ploop;
    }

    line += 4; /* "loop" */
    while (isspace((unsigned char) *line)) {
	line++;
    }

#if LOOP_DEBUG
    fprintf(stderr, "ready for parse: '%s'\n", line);
#endif

    /* try parsing the loop line in various ways */

    if (sscanf(line, "while %8[^ <>=]%8[ <>=] %8s", lvar, op, rvar) == 3) {
	err = parse_as_while_loop(loop, pdinfo, lvar, rvar, op);
    }

    else if (sscanf(line, "%c = %8[^.]..%8s", &ichar, op, rvar) == 3) {
	err = parse_as_indexed_loop(loop, pdinfo, (const double **) *pZ, 
				    ichar, NULL, op, rvar);
    }	

    else if (sscanf(line, "for %8[^= ] = %8[^.]..%8s", lvar, op, rvar) == 3) {
	err = parse_as_indexed_loop(loop, pdinfo, (const double **) *pZ, 
				    0, lvar, op, rvar);
    }

    else if (!strncmp(line, "foreach", 7)) {
	err = parse_as_each_loop(loop, pdinfo, line + 7);
    }    

    else if (!strncmp(line, "for", 3)) {
	err = parse_as_for_loop(loop, pdinfo, pZ, line + 3);
    }

    else if (sscanf(line, "%8s", lvar) == 1) {
	err = parse_as_count_loop(loop, pdinfo, (const double **) *pZ, 
				  lvar);
    }

    /* out of options, complain */
    else {
	printf("parse_loopline: failed on '%s'\n", line);
	strcpy(gretl_errmsg, _("No valid loop condition was given."));
	err = 1;
    }

    if (!err && !na(loop->init.val) && !na(loop->right.val)) {
	int nt = loop->right.val - loop->init.val;

	if (loop->type != FOR_LOOP && nt <= 0) {
	    strcpy(gretl_errmsg, _("Loop count missing or invalid\n"));
	    err = 1;
	}
    }

    if (!err) {
	/* allocates loop->lines, loop->ci */
	err = prepare_loop_for_action(loop);
    }

    if (err) {
	if (loop != ploop) {
	    free(loop->lines);
	    free(loop->ci);
	    free(loop);
	}
	loop = NULL;
    }

    if (!err && opt != OPT_NONE) {
	set_loop_opts(loop, opt);
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
	/* FIXME */
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

static int 
eval_numeric_condition (int op, double xl, double xr)
{
    int cont = 0;

#if LOOP_DEBUG
    char opstr[4];

    if (op == OP_GT) strcpy(opstr, "GT");
    else if (op == OP_GTE) strcpy(opstr, "GTE");
    else if (op == OP_LTE) strcpy(opstr, "LTE");
    else if (op == OP_LT) strcpy(opstr, "LT");
    else strcpy(opstr, "??");
    fprintf(stderr, "** eval_numeric_condition: xl=%g, xr=%g, "
	    "op=%s\n", xl, xr, opstr);
#endif

    if (op == OP_GT) {
	cont = (xl > xr);
    } else if (op == OP_GTE) {
	cont = (xl >= xr);
    } else if (op == OP_LTE) {
	cont = (xl <= xr);
    } else if (op == OP_LT) {
	cont = (xl < xr);
    }

    return cont;
}

/**
 * loop_condition:
 * @loop: pointer to loop commands struct.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Check whether a looping condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

static int 
loop_condition (LOOPSET *loop, double **Z, DATAINFO *pdinfo)
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
	} else {
	    double lval, rval = loop_rval(loop, (const double **) Z);

	    if (loop->type == FOR_LOOP && loop->iter > 0) {
		Z[loop->left.vnum][0] += loop_incrval(loop, (const double **) Z);
		lval = Z[loop->left.vnum][0];
	    } else {
		lval = loop_lval(loop, (const double **) Z);
	    }

	    ok = eval_numeric_condition(loop->ineq, lval, rval);
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
    loop->brk = 0;

    controller_init(&loop->left);
    controller_init(&loop->right);
    controller_init(&loop->init);
    controller_init(&loop->incr);

    loop->incrsgn = 1;
    loop->ineq = 0;

    loop->ncmds = 0;
    loop->nmod = 0;
    loop->nprn = 0;
    loop->nstore = 0;

    loop->next_model = 0;
    loop->next_print = 0;

    loop->lines = NULL;
    loop->ci = NULL;

    loop->eachstrs = NULL;

    loop->models = NULL;
    loop->lmodels = NULL;
    loop->prns = NULL;

    loop->storefile[0] = '\0';
    loop->storename = NULL;
    loop->storelbl = NULL;
    loop->storeval = NULL;

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

void gretl_loop_destroy (LOOPSET *loop)
{
    int i;

    for (i=0; i<loop->n_children; i++) {
	gretl_loop_destroy(loop->children[i]);
    }

    if (loop->lines != NULL) {
	for (i=0; i<loop->ncmds; i++) {
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
	for (i=0; i<loop->nmod; i++) {
	    gretl_model_free(loop->models[i]);
	}
	free(loop->models);
    } 

    if (loop->lmodels != NULL) {
	for (i=0; i<loop->nmod; i++) {
	    free_loop_model(&loop->lmodels[i]);
	}
	free(loop->lmodels);
    }

    if (loop->prns != NULL) {
	for (i=0; i<loop->nprn; i++) { 
	    free_loop_print(&loop->prns[i]);
	}
	free(loop->prns);
    }

    if (loop->storename != NULL) {
	for (i=0; i<loop->nstore; i++) {
	    free(loop->storename[i]);
	}
	free(loop->storename);
    }

    if (loop->storelbl != NULL) {
	for (i=0; i<loop->nstore; i++) {
	    free(loop->storelbl[i]);
	}
	free(loop->storelbl);
    }

    if (loop->storeval != NULL) { 
	free(loop->storeval);
    }

    if (loop->children != NULL) {
	free(loop->children);
    }

    free(loop);
}

/**
 * loop_model_init:
 * @lmod: pointer to struct to initialize.
 * @pmod: model to take as basis.
 * @id: ID number to assign to @lmod.
 *
 * Initialize a #LOOP_MODEL struct, based on @pmod.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

static int loop_model_init (LOOP_MODEL *lmod, const MODEL *pmod,
			    int id)
{
    int i, nc = pmod->ncoeff;

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

    for (i=0; i<nc; i++) {
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

    lmod->ID = id;

    return 0;

 cleanup:
    free(lmod->ssq_coeff);
    free(lmod->sum_sderr);
    free(lmod->ssq_sderr);

    return E_ALLOC;
}

/**
 * loop_print_init:
 * @lprn: pointer to struct to initialize.
 * @list: list of variables to be printed.
 * @id: ID number to assign to @lprn.
 *
 * Initialize a #LOOP_PRINT struct.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

static int loop_print_init (LOOP_PRINT *lprn, const int *list, int id)
{
    int i;

    lprn->list = gretl_list_copy(list);
    if (lprn->list == NULL) return 1;

    lprn->sum = malloc(list[0] * sizeof *lprn->sum);
    if (lprn->sum == NULL) goto cleanup;

    lprn->ssq = malloc(list[0] * sizeof *lprn->ssq);
    if (lprn->ssq == NULL) goto cleanup;

    for (i=0; i<list[0]; i++) { 
#ifdef ENABLE_GMP
	mpf_init(lprn->sum[i]);
	mpf_init(lprn->ssq[i]);
#else
	lprn->sum[i] = lprn->ssq[i] = 0.0;
#endif
    }

    lprn->ID = id;

    return 0;

 cleanup:
    free(lprn->list);
    free(lprn->sum);
    free(lprn->ssq);

    return 1;
}

/**
 * loop_store_init:
 * @loop: pointer to loop struct.
 * @fname: name of file in which to store data.
 * @list: list of variables to be stored (written to file).
 * @pdinfo: data information struct.
 *
 * Set up @loop for saving a set of variables.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

static int loop_store_init (LOOPSET *loop, const char *fname, 
			    const int *list, DATAINFO *pdinfo)
{
    int i, tot = list[0] * loop->itermax;

    loop->storefile[0] = '\0';
    strncat(loop->storefile, fname, MAXLEN - 1);

    loop->storename = malloc(list[0] * sizeof *loop->storename);
    if (loop->storename == NULL) return 1;

    loop->storelbl = malloc(list[0] * sizeof *loop->storelbl);
    if (loop->storelbl == NULL) goto cleanup;

    loop->storeval = malloc(tot * sizeof *loop->storeval);
    if (loop->storeval == NULL) goto cleanup;

    for (i=0; i<list[0]; i++) {
	char *p;

	loop->storename[i] = malloc(VNAMELEN);
	if (loop->storename[i] == NULL) goto cleanup;

	strcpy(loop->storename[i], pdinfo->varname[list[i+1]]);

	loop->storelbl[i] = malloc(MAXLABEL);
	if (loop->storelbl[i] == NULL) goto cleanup;

	strcpy(loop->storelbl[i], VARLABEL(pdinfo, list[i+1]));
	if ((p = strstr(loop->storelbl[i], "(scalar)"))) {
	    *p = 0;
	}
    }

    return 0;

 cleanup:
    free(loop->storename);
    free(loop->storelbl);
    free(loop->storeval);

    return 1;
}

static int add_loop_model_record (LOOPSET *loop, int cmdnum)
{
    MODEL **lmods;
    int err = 0;
    int nm = loop->nmod + 1;

    lmods = realloc(loop->models, nm * sizeof *lmods);
    if (lmods == NULL) {
	err = 1;
    } else {
	loop->models = lmods;
	loop->models[loop->nmod] = gretl_model_new();
	if (loop->models[loop->nmod] == NULL) {
	    err = 1;
	} else {
	    (loop->models[loop->nmod])->ID = cmdnum;
	}
    }

    if (!err) {
	loop->nmod += 1;
    }

    return err;
}

static int add_loop_model (LOOPSET *loop)
{
    int nm = loop->nmod + 1;
    int err = 0;

    loop->lmodels = realloc(loop->lmodels, nm * sizeof *loop->lmodels);
    if (loop->lmodels == NULL) {
	err = 1;
    } else {
	loop->nmod += 1;
    }

    return err;
}

/**
 * update_loop_model:
 * @loop: pointer to loop struct.
 * @cmdnum: sequential index number of the command within @loop.
 * @pmod: contains estimates from the current iteration.
 *
 * Update a #LOOP_MODEL belonging to @loop, based on the results
 * in @pmod.
 *
 * Returns: 0 on successful completion.
 */

static int update_loop_model (LOOPSET *loop, int cmdnum, MODEL *pmod)
{
    int j, i = get_modnum_by_cmdnum(loop, cmdnum);
    LOOP_MODEL *lmod;
#ifdef ENABLE_GMP
    mpf_t m;

    mpf_init(m);
#endif

    lmod = &loop->lmodels[i];

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

    return 0;
}

static int add_loop_print (LOOPSET *loop, const int *list, int cmdnum)
{
    LOOP_PRINT *prns;
    int np = loop->nprn + 1;
    int err = 0;

    prns = realloc(loop->prns, np * sizeof *prns);
    if (prns == NULL) {
	return 1;
    }

    loop->prns = prns;

    if (loop_print_init(&loop->prns[loop->nprn], list, cmdnum)) { 
	strcpy(gretl_errmsg, _("Failed to initalize print struct for loop\n"));
	err = 1;
    }

    if (!err) {
	loop->nprn += 1;
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
 *
 * Returns: 0 on successful completion.
 */

static int update_loop_print (LOOPSET *loop, int cmdnum, 
			      const int *list, double ***pZ, 
			      const DATAINFO *pdinfo)
{
    int j, t, i = get_prnnum_by_id(loop, cmdnum);
    LOOP_PRINT *lprn = &loop->prns[i];
#ifdef ENABLE_GMP
    mpf_t m;

    mpf_init(m);
#endif
    
    for (j=1; j<=list[0]; j++) {
	if (pdinfo->vector[list[j]]) t = pdinfo->t1;
	else t = 0;
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

    return 0;
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
    int i;

    if (loop->type != COUNT_LOOP && !(loop_is_quiet(loop))) {
	pprintf(prn, _("\nNumber of iterations: %d\n\n"), loop->iter);
    }

    for (i=0; i<loop->ncmds; i++) {
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
		MODEL *pmod = loop->models[loop->next_model];

		set_model_id(pmod);
		printmodel(pmod, pdinfo, opt, prn);
		loop->next_model += 1;
	    }	    
	}

	if (loop_is_progressive(loop)) {
	    if (loop->ci[i] == OLS || loop->ci[i] == LAD ||
		loop->ci[i] == HSK || loop->ci[i] == HCCM || 
		loop->ci[i] == WLS) {
		print_loop_model(&loop->lmodels[loop->next_model], 
				 loop->itermax, pdinfo, prn);
		loop->next_model += 1;
	    } else if (loop->ci[i] == PRINT) {
		print_loop_prn(&loop->prns[loop->next_print], 
			       loop->itermax, pdinfo, prn);
		loop->next_print += 1;
	    } else if (loop->ci[i] == STORE) {
		print_loop_store(loop, prn);
	    }
	}
    }
}

static void free_loop_model (LOOP_MODEL *lmod)
{
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
    int nb = 1 + (loop->ncmds + 1) / LOOP_BLOCK;
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

/**
 * add_to_loop:
 * @line: command line.
 * @ci: command index number.
 * @opt: option flag(s) associated with the command.
 * @pdinfo: dataset information.
 * @pZ: pointer to data matrix.
 * @loop: #LOOPSET to which command should be added
 * @loopstack: pointer to integer stacking level.
 * @looprun: pointer to integer switch for running loop.
 *
 * Add line and command index to accumulated loop buffer.
 *
 * Returns: pointer to loop struct on success, %NULL on failure.
 */

LOOPSET *add_to_loop (char *line, int ci, gretlopt opt,
		      DATAINFO *pdinfo, double ***pZ,
		      LOOPSET *loop, int *loopstack, int *looprun)
{
    LOOPSET *lret = loop;

    *gretl_errmsg = '\0';

#if LOOP_DEBUG
    fprintf(stderr, "add_to_loop: loop = %p, loopstack = %d, line = '%s'\n", 
	    (void *) loop, *loopstack, line);
#endif

    if (ci == LOOP) {
	/* starting a new loop, possibly inside another */
	lret = parse_loopline(line, loop, *loopstack, pdinfo, pZ);
	if (lret == NULL) {
	    if (*gretl_errmsg == '\0') {
		gretl_errmsg_set(_("No valid loop condition was given."));
		gretl_errno = E_PARSE;
	    }
	    goto bailout;
	} else {
	    set_loop_opts(lret, opt);
	    *loopstack += 1;
	}
    } else if (ci == ENDLOOP) {
	*loopstack -= 1;
	if (*loopstack == 0) {
	    *looprun = 1;
	} else {
	    lret = loop->parent;
	}
    } 

    if (loop != NULL) {
	int nc = loop->ncmds;

	if ((nc + 1) % LOOP_BLOCK == 0) {
	    if (add_more_loop_lines(loop)) {
		gretl_errmsg_set(_("Out of memory!"));
		goto bailout;
	    }
	}

	loop->lines[nc] = malloc(MAXLEN);
	if (loop->lines[nc] == NULL) {
	    gretl_errmsg_set(_("Out of memory!"));
	    goto bailout;
	}

	top_n_tail(line);

	if (ci == PRINT && loop->type != COUNT_LOOP) {
	    /* fixme: what's going on here? */
	    loop->ci[nc] = 0;
	} else {
	    loop->ci[nc] = ci;
	}

	loop->lines[nc][0] = '\0';

	if (opt) {
	    const char *flagstr = print_flags(opt, ci);

	    if (strlen(line) + strlen(flagstr) >= MAXLEN) {
		goto bailout;
	    } else {
		sprintf(loop->lines[nc], "%s%s", line, flagstr);
	    }
	} else {
	    strcpy(loop->lines[nc], line);
	}

	loop->ncmds += 1;

#if LOOP_DEBUG
	fprintf(stderr, "loop: ncmds=%d, line[%d] = '%s'\n",
		loop->ncmds, nc, loop->lines[nc]);
#endif
    }

    return lret;

 bailout:
    
    if (loop != NULL) {
	gretl_loop_destroy(loop);
    }

    return NULL;
}

/* FIXME: additional model info such as robust std errs method */

static void print_loop_model (LOOP_MODEL *lmod, int loopnum,
			      const DATAINFO *pdinfo, PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    int i;

    ntodate(startdate, lmod->model0->t1, pdinfo);
    ntodate(enddate, lmod->model0->t2, pdinfo);

    pputc(prn, '\n');
    pprintf(prn, _("%s estimates using the %d observations %s-%s\n"),
	    _(estimator_string(lmod->model0->ci, prn)), lmod->model0->nobs, 
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

    pprintf(prn, " %3d) %8s ", lmod->model0->list[c+2], 
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

    pprintf(prn, " %3d) %8s ", lmod->model0->list[c+2], 
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
	pprintf(prn, " %8s ", pdinfo->varname[lprn->list[i]]);
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
	pprintf(prn, " %8s ", pdinfo->varname[lprn->list[i]]);
	pprintf(prn, "%#14g %#14g\n", (double) mean, (double) sd);
    }
#endif
    pputc(prn, '\n');
}

static int print_loop_store (LOOPSET *loop, PRN *prn)
{
    int i, t;
    FILE *fp;
    char gdtfile[MAXLEN], infobuf[1024];
    char *xmlbuf = NULL;
    time_t writetime;

    /* organize filename */
    if (loop->storefile[0] == '\0') {
	sprintf(gdtfile, "%sloopdata.gdt", gretl_user_dir());	
    } else {
	strcpy(gdtfile, loop->storefile);
    }

    if (strchr(gdtfile, '.') == NULL) {
	strcat(gdtfile, ".gdt");
    }

    fp = gretl_fopen(gdtfile, "w");
    if (fp == NULL) return 1;

    writetime = time(NULL);

    pprintf(prn, _("printing %d values of variables to %s\n"), 
	    loop->itermax, gdtfile);

    fprintf(fp, "<?xml version=\"1.0\"?>\n"
	    "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
	    "<gretldata name=\"%s\" frequency=\"1\" "
	    "startobs=\"1\" endobs=\"%d\" ", 
	    gdtfile, loop->itermax);

    fprintf(fp, "type=\"cross-section\">\n");

    sprintf(infobuf, "%s %s", _("simulation data written"),
	    print_time(&writetime)); 
    xmlbuf = gretl_xml_encode(infobuf);
    fprintf(fp, "<description>\n%s\n</description>\n", xmlbuf);
    free(xmlbuf);

    gretl_push_c_numeric_locale();

    /* print info on variables */
    fprintf(fp, "<variables count=\"%d\">\n", loop->nstore);

    for (i=0; i<loop->nstore; i++) {
	xmlbuf = gretl_xml_encode(loop->storename[i]);
	fprintf(fp, "<variable name=\"%s\"", xmlbuf);
	free(xmlbuf);
	xmlbuf = gretl_xml_encode(loop->storelbl[i]);
	fprintf(fp, "\n label=\"%s\"/>\n", xmlbuf);
	free(xmlbuf);
    }

    fputs("</variables>\n", fp);

    /* print actual data */
    fprintf(fp, "<observations count=\"%d\" labels=\"false\">\n",
	    loop->itermax);

    for (t=0; t<loop->itermax; t++) {
	double x;

	fputs("<obs>", fp);
	for (i=0; i<loop->nstore; i++) {
	    x = loop->storeval[loop->itermax * i + t];
	    if (na(x)) {
		fputs("NA ", fp);
	    } else {
		fprintf(fp, "%g ", x);
	    }
	}
	fputs("</obs>\n", fp);
    }

    fprintf(fp, "</observations>\n</gretldata>\n");

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static int get_prnnum_by_id (LOOPSET *loop, int id)
{
    int i;

    for (i=0; i<loop->nprn; i++) {
	if (loop->prns[i].ID == id) return i;
    }
    return -1;
}

/**
 * get_modnum_by_cmdnum:
 * @loop: pointer to loop struct.
 * @cmdnum: sequential index of command within @loop.
 *
 * Determine the ID number of a model within a "while" loop construct.
 *
 * Returns: model ID number, or -1 in case of no match.
 */

static int get_modnum_by_cmdnum (LOOPSET *loop, int cmdnum)
{
    int i;

    if (loop_is_progressive(loop)) {
	for (i=0; i<loop->nmod; i++) {
	    if (loop->lmodels[i].ID == cmdnum) {
		return i;
	    }
	}
    } else {
	for (i=0; i<loop->nmod; i++) {
	    if ((loop->models[i])->ID == cmdnum) {
		return i;
	    }
	}
    }

    return -1;
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

#if LOOP_DEBUG
    fprintf(stderr, "subst_dollar_targ:\n original: '%s'\n", str);
#endif

    if (loop->type == FOR_LOOP) {
	sprintf(targ, "$%s", pdinfo->varname[loop->left.vnum]);
	targlen = strlen(targ);
    } else if (indexed_loop(loop)) {
	targ[0] = '$';
	targ[1] = loop->ichar;
	targlen = 2;
	idx = loop->init.val + loop->iter;
    } else {
	return 1;
    }

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
	    sprintf(ins, "%g", Z[loop->left.vnum][0]); /* scalar */

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

#if LOOP_DEBUG
    fprintf(stderr, " after: '%s'\n", str);
#endif

    return err;
}

static void 
loop_add_storevals (const int *list, LOOPSET *loop,
		    const double **Z, DATAINFO *pdinfo)
{
    int i, sv;

    for (i=0; i<list[0]; i++) {
	sv = i * loop->itermax + loop->iter;
	if (pdinfo->vector[list[i+1]]) { 
	    loop->storeval[sv] = Z[list[i+1]][pdinfo->t1 + 1];
	} else {
	    loop->storeval[sv] = Z[list[i+1]][0];
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

    if (loop->left.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->left, pdinfo);
    }

    if (!err && loop->right.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->right, pdinfo);
    }

    if (!err && loop->init.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->init, pdinfo);
    }

    if (!err && loop->incr.vnum == LOOP_VAL_UNDEF) {
	err = clr_attach_var(&loop->incr, pdinfo);
    }

    return err;
}

static int 
top_of_loop (LOOPSET *loop, double **Z, const DATAINFO *pdinfo)
{
    int err;

    loop->iter = 0;

    err = connect_loop_control_vars(loop, pdinfo);

    if (!err) {
	if (loop->init.vnum > 0) {
	    loop->init.val = 
		loop_controller_get_value(&loop->init, (const double **) Z);
	} 

	if (loop->type == FOR_LOOP) {
	    Z[loop->left.vnum][0] = loop->init.val;
	}

	if (loop->type == COUNT_LOOP || indexed_loop(loop)) {
	    double maxval = loop_rval(loop, (const double **) Z);

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
	    loop_scalar_index(loop->ichar, LOOP_SCALAR_PUT, 
			      (int) loop->init.val);
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

/* Below, in loop_exec, as of October 29, 2005, the handling of the
   multiple models is probably not right yet -- needs more work in
   light of the new object naming/saving mechanism in gretl.
*/

int loop_exec (LOOPSET *loop, char *line,
	       double ***pZ, DATAINFO **ppdinfo, 
	       MODEL **models, PRN *prn)
{
    CMD cmd;
    MODEL *lastmod = models[0];
    GRETL_VAR *var;
    char errline[MAXLINE];
    char linecpy[MAXLINE];
    int m = 0;
    int modnum = 0;
    int err = 0;

    if (loop == NULL) {
	pputs(prn, "Got a NULL loop\n");
	return 1;
    }

    if (loop->ncmds == 0) {
	pputs(prn, _("No commands in loop\n"));
	return 0;
    }

    err = gretl_cmd_init(&cmd);
    if (err) {
	return err;
    }

    set_loop_on();

#if LOOP_DEBUG
    fprintf(stderr, "loop_exec: loop = %p\n", (void *) loop);
#endif

    err = top_of_loop(loop, *pZ, *ppdinfo);

    while (!err && loop_condition(loop, *pZ, *ppdinfo)) {
	int childnum = 0;
	int j;

#if LOOP_DEBUG
	fprintf(stderr, "top of loop: iter = %d\n", loop->iter);
#endif

	if (gretl_echo_on() && indexed_loop(loop)) {
	    print_loop_progress(loop, *ppdinfo, prn);
	}

	for (j=0; !err && j<loop->ncmds; j++) {
#if LOOP_DEBUG
	    fprintf(stderr, "loop->lines[%d] = '%s'\n", j, loop->lines[j]);
#endif
	    strcpy(linecpy, loop->lines[j]);
	    strcpy(errline, loop->lines[j]);

	    err = make_dollar_substitutions(linecpy, loop, 
					    (const double **) *pZ,
					    *ppdinfo);
	    if (err) break;

	    /* We already have the "ci" index recorded, but this line
	       will do some checking that hasn't been done earlier.
	    */

	    err = parse_command_line(linecpy, &cmd, pZ, *ppdinfo);

	    if (cmd.ci < 0) {
		continue;
	    } else if (err) {
		break;
	    }

	    if (gretl_echo_on() && indexed_loop(loop)) {
		if (cmd.ci == ENDLOOP) {
		    pputc(prn, '\n');
		} else {
		    echo_cmd(&cmd, *ppdinfo, linecpy, 0, prn);
		}
	    }

	    switch (cmd.ci) {

	    case ADF: 
	    case COINT2: 
	    case COINT: 
	    case CORR:
	    case CRITERIA: 
	    case CRITICAL: 
	    case DIFF: 
	    case HURST:	
	    case KPSS:
	    case LABEL:
	    case LAGS: 
	    case LDIFF: 
	    case LOGS:
	    case MEANTEST: 
	    case MULTIPLY: 
	    case OUTFILE:
	    case PCA:
	    case RHODIFF:
	    case RUNS: 
	    case SIM:
	    case SPEARMAN: 
	    case SQUARE: 
	    case SUMMARY:
	    case VARTEST: 
		err = simple_commands(&cmd, linecpy, pZ, *ppdinfo, prn);
		break;

	    case LOOP:
		err = loop_exec(loop->children[childnum++], NULL,
				pZ, ppdinfo, models, prn);
		break;

	    case BREAK:
		loop->brk = 1;
		break;

	    case ENDLOOP:
		/* no-op */
		break;

	    case GENR:
		err = generate(linecpy, pZ, *ppdinfo, cmd.opt);
		if (loop_is_verbose(loop) && !err) { 
		    print_gretl_msg(prn);
		}
		break;

	    case ARMA:
	    case CORC:
	    case GARCH:
	    case HCCM:
	    case HILU:
	    case HSK:
	    case LAD:
	    case OLS:
	    case PWE:
	    case WLS:
		/* if this is the first time round, allocate space
		   for each loop model */
		if (loop->iter == 0) {
		    if (loop_is_progressive(loop)) {
			err = add_loop_model(loop);
		    } else if (cmd.opt & OPT_P) {
			err = add_loop_model_record(loop, j);
		    }
		    if (err) {
			break;
		    }
		} 

		/* estimate the model called for, using models[0] */
		clear_model(models[0]);

		if (cmd.ci == OLS || cmd.ci == WLS || cmd.ci == HCCM) {
		    *models[0] = lsq(cmd.list, pZ, *ppdinfo, cmd.ci, cmd.opt, 0.0);
		} else if (cmd.ci == LAD) {
		    *models[0] = lad(cmd.list, pZ, *ppdinfo);
		} else if (cmd.ci == HSK) {
		    *models[0] = hsk_func(cmd.list, pZ, *ppdinfo);
		} else if (cmd.ci == ARMA) {
		    *models[0] = arma(cmd.list, (const double **) *pZ, *ppdinfo, 
				      cmd.opt, prn);
		} else if (cmd.ci == GARCH) {
		    *models[0] = garch(cmd.list, pZ, *ppdinfo, cmd.opt, prn);
		} else if (cmd.ci == CORC || cmd.ci == HILU || cmd.ci == PWE) {
		    double rho = estimate_rho(cmd.list, pZ, *ppdinfo, cmd.ci,
					      &err, cmd.opt, prn);
		    if (err) {
			break;
		    }
		    *models[0] = lsq(cmd.list, pZ, *ppdinfo, cmd.ci, cmd.opt, rho);
		}

		if ((err = models[0]->errcode)) {
		    break;
		} 

		if (loop_is_progressive(loop)) {
		    if (loop->iter == 0 && loop_model_init(&loop->lmodels[loop->nmod - 1], 
							   models[0], j)) { 
			gretl_errmsg_set(_("Failed to initialize model for loop\n"));
			err = 1;
			break;
		    } else if (update_loop_model(loop, j, models[0])) {
			gretl_errmsg_set(_("Failed to add results to loop model\n"));
			err = 1;
			break;
		    }
		    set_as_last_model(models[0], EQUATION);
		} else if (cmd.opt & OPT_P) {
		    /* deferred printing of model results */
		    m = get_modnum_by_cmdnum(loop, j);
		    swap_models(&models[0], &loop->models[m]); /* direction? */
		    loop->models[m]->ID = j;
		    lastmod = loop->models[m];
		    set_as_last_model(loop->models[m], EQUATION); /* FIXME? */
		    model_count_minus();
		} else {
		    models[0]->ID = ++modnum; /* FIXME? */
		    printmodel(models[0], *ppdinfo, cmd.opt, prn);
		    lastmod = models[0];
		    set_as_last_model(models[0], EQUATION);
		}
		break;

	    case ADD:
	    case OMIT:
		if (loop_is_progressive(loop)) {
		    err = 1;
		    break;
		}
		/* FIXME: this needs work, and should only be allowed
		   under certain conditions */
		clear_model(models[1]);
		if (cmd.ci == ADD || cmd.ci == ADDTO) {
		    err = add_test(cmd.list, models[0], models[1], 
				   pZ, *ppdinfo, cmd.opt, prn);
		} else {
		    err = omit_test(cmd.list, models[0], models[1],
				    pZ, *ppdinfo, cmd.opt, prn);
		}
		if (err) {
		    errmsg(err, prn);
		    clear_model(models[1]);
		} else {
		    swap_models(&models[0], &models[1]);
		    lastmod = models[0];
		    set_as_last_model(models[0], EQUATION);
		    clear_model(models[1]);
		}
		break;	

	    case MLE:
	    case NLS:
		if (loop_is_progressive(loop) || (cmd.opt & OPT_P)) {
		    err = 1;
		} else {
		    err = nls_parse_line(cmd.ci, linecpy, (const double **) *pZ, 
					 *ppdinfo, prn);
		    if (err) {
			errmsg(err, prn);
		    } else {
			gretl_cmd_set_context(&cmd, cmd.ci);
		    }
		}
		break;

	    case END:
		if (!strcmp(cmd.param, "nls") || !strcmp(cmd.param, "mle")) {
		    clear_model(models[0]);
		    *models[0] = nls(pZ, *ppdinfo, cmd.opt, prn);
		    if ((err = (models[0])->errcode)) {
			errmsg(err, prn);
		    } else {
			printmodel(models[0], *ppdinfo, cmd.opt, prn);
			lastmod = models[0];
			set_as_last_model(models[0], EQUATION);
		    }
		} else {
		    err = 1;
		}
		break;

	    case PRINT:
		if (cmd.param[0] != '\0') {
		    err = simple_commands(&cmd, linecpy, pZ, *ppdinfo, prn);
		} else if (loop_is_progressive(loop)) {
		    if (loop->iter == 0) {
			if ((err = add_loop_print(loop, cmd.list, j))) {
			    break;
			}
		    }
		    if (update_loop_print(loop, j, cmd.list, pZ, *ppdinfo)) {
			gretl_errmsg_set(_("Failed to add values to print loop\n"));
			err = 1;
		    }
		} else {
		    err = printdata(cmd.list, (const double **) *pZ, *ppdinfo, 
				    cmd.opt, prn);
		}
		break;

	    case PRINTF:
		err = do_printf(linecpy, pZ, *ppdinfo, prn);
		break;

	    case SMPL:
		if (cmd.opt) {
		    err = restore_full_sample(pZ, ppdinfo, cmd.opt);
		    if (err) {
			errmsg(err, prn);
			break;
		    } else {
			err = restrict_sample(linecpy, pZ, ppdinfo, 
					      cmd.list, cmd.opt);
		    }
		} else if (!strcmp(linecpy, "smpl full") ||
			   !strcmp(linecpy, "smpl --full")) {
		    err = restore_full_sample(pZ, ppdinfo, OPT_C);
		} else { 
		    err = set_sample(linecpy, (const double **) *pZ,
				     *ppdinfo);
		}

		if (err) {
		    errmsg(err, prn);
		} else if (1 || gretl_echo_on()) {
		    print_smpl(*ppdinfo, get_full_length_n(), prn);
		}
		break;

	    case STORE:
		if (loop_is_progressive(loop)) {
		    if (loop->iter == 0) {
			loop->nstore = cmd.list[0];
			if (loop_store_init(loop, cmd.param, cmd.list, *ppdinfo)) {
			    err = 1;
			}
		    }
		    if (!err) {
			loop_add_storevals(cmd.list, loop, 
					   (const double **) *pZ, 
					   *ppdinfo);
		    }
		} else {
		    simple_commands(&cmd, linecpy, pZ, *ppdinfo, prn);
		}
		break;

	    case PVALUE:
		batch_pvalue(linecpy, (const double **) *pZ, *ppdinfo, prn, OPT_NONE);
		break;

	    case VAR:
	    case VECM:
		if (cmd.ci == VAR) {
		    var = full_VAR(atoi(cmd.param), cmd.list, pZ, *ppdinfo, 
				   cmd.opt, prn);
		} else {
		    var = vecm(atoi(cmd.param), atoi(cmd.extra), cmd.list, 
			   pZ, *ppdinfo, cmd.opt, prn);
		}
		if (var == NULL) {
		    err = 1;
		} else {
		    set_as_last_model(var, VAR);
		}
		break;

	    default: 
		/* not reachable (since commands were screened in advance) */
		pprintf(prn, _("command: '%s'\nThis is not available in a loop.\n"),
			linecpy);
		err = 1;
		break;

	    } /* end switch on specific command number */

	} /* end execution of commands within loop */

	if (err && get_halt_on_error() == 0) {
	    errmsg(err, prn);
	    err = 0;
	}

	loop->iter += 1;

	if (indexed_loop(loop)) {
	    loop_scalar_index(loop->ichar, LOOP_SCALAR_INCR, 1);
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
	print_loop_results(loop, *ppdinfo, prn); 
    }

    if (lastmod != models[0]) {
	/* to get commands that reference models[0], after the loop
	   has finished, to come out right
	*/
	swap_models(&models[0], &loop->models[m]);
    }

    gretl_cmd_free(&cmd);

    if (line != NULL) {
	*line = '\0';
    } 

    set_active_loop(loop->parent);
    set_loop_off();

    if (get_halt_on_error()) {
	return err;
    } else {
	return 0;
    }
}

static int ichar_in_parentage (const LOOPSET *loop, int c)
{
    while ((loop = loop->parent) != NULL) {
	if (indexed_loop(loop) && loop->ichar == c) {
	    return 1;
	}
    }

    return 0;
}

/* apparatus relating to retrieval of loop index value in genr */

static int loop_scalar_index (int c, int opt, int put)
{
    static int idx[N_LOOP_INDICES];
    int i = loop_index_char_pos(c);
    int ret = -1;

#if LOOP_DEBUG
    fprintf(stderr, "loop_scalar_index: c='%c', opt=%d, put=%d\n", 
	    c, opt, put);
#endif

    if (i >= 0) {
	if (opt == LOOP_SCALAR_PUT) {
	    idx[i] = put;
	} else if (opt == LOOP_SCALAR_INCR) {
	    idx[i] += put;
	}
	ret = idx[i];
    }

#if LOOP_DEBUG
    fprintf(stderr, "loop_scalar_index: returning %d\n", ret);
#endif

    return ret;
}

int loop_scalar_read (int c)
{
    return loop_scalar_index(c, LOOP_SCALAR_READ, 0);
}

static int indexed_loop_record (LOOPSET *loop, int set, int test)
{
    static LOOPSET *active_loop;
    int ret = 0;

    if (set) {
	active_loop = loop;
    } else if (active_loop != NULL) {
	if (indexed_loop(active_loop) && active_loop->ichar == test) {
	    ret = 1;
	} else if (ichar_in_parentage(active_loop, test)) {
	    ret = 1;
	}
    }

#if LOOP_DEBUG
    if (set) {
	fprintf(stderr, "indexed_loop_record: set active_loop = %p\n",
		(void *) loop);
    } else {
	fprintf(stderr, "indexed_loop_record: returning %d for test='%c'\n",
		ret, test);
    }	
#endif

    return ret;
}

static void set_active_loop (LOOPSET *loop)
{
    indexed_loop_record(loop, 1, 0);
}

int is_active_index_loop_char (int c)
{
    return indexed_loop_record(NULL, 0, c);
}

#define IFDEBUG 0

/* if-then stuff - conditional execution */

int if_eval (const char *line, double ***pZ, DATAINFO *pdinfo)
{
    char formula[MAXLEN];
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

    sprintf(formula, "__iftest=%s", line);

    err = generate(formula, pZ, pdinfo, OPT_P);
#if IFDEBUG
    fprintf(stderr, "if_eval: generate returned %d\n", err);
#endif

    if (err) {
	strcpy(gretl_errmsg, _("error evaluating 'if'"));
    } else {
	int v = varindex(pdinfo, "iftest");

	if (v < pdinfo->v) {
	    double val = (*pZ)[v][0];

	    if (na(val)) {
		strcpy(gretl_errmsg, _("indeterminate condition for 'if'"));
	    } else {
		ret = (int) val;
	    }
	    dataset_drop_last_variables(1, pZ, pdinfo);
	}
    }

#if IFDEBUG
    fprintf(stderr, "if_eval: returning %d\n", ret);
#endif

    return ret;
}

#define IF_DEPTH 9

int ifstate (int code)
{
    static unsigned char T[IF_DEPTH];
    static unsigned char got_if[IF_DEPTH];
    static unsigned char got_else[IF_DEPTH];
    static unsigned char indent;

    if (code == RELAX) {
	indent = 0;
    } else if (code == SET_FALSE || code == SET_TRUE) {
	indent++;
	if (indent >= IF_DEPTH) {
	    fprintf(stderr, "if depth (%d) exceeded\n", IF_DEPTH);
	    return 1; /* too deeply nested */
	}
	T[indent] = (code == SET_TRUE);
	got_if[indent] = 1;
	got_else[indent] = 0;
    } else if (code == SET_ELSE) {
	if (got_else[indent] || !got_if[indent]) {
	    strcpy(gretl_errmsg, "Unmatched \"else\"");
	    return 1; 
	}
	T[indent] = !T[indent];
	got_else[indent] = 1;
    } else if (code == SET_ENDIF) {
	if (!got_if[indent] || indent == 0) {
	    strcpy(gretl_errmsg, "Unmatched \"endif\"");
	    return 1; 
	}
	got_if[indent] = 0;
	got_else[indent] = 0;
	indent--;
    } else if (code == IS_FALSE) {
	int i;

	for (i=1; i<=indent; i++) {
	    if (T[i] == 0) return 1;
	}
    }

    return 0;
}
