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
#include "gretl_private.h"
#include "libset.h"

#include <time.h>
#include <unistd.h>

#undef LOOP_DEBUG

#if defined(ENABLE_GMP)
# include <gmp.h>
typedef mpf_t bigval;
#elif defined(HAVE_LONG_DOUBLE)
typedef long double bigval;
#else
typedef double bigval;
#endif

enum inequalities {
    GT = 1,
    LT
};

enum loop_types {
    COUNT_LOOP,
    WHILE_LOOP,
    INDEX_LOOP,
    FOR_LOOP
};

typedef struct {
    int ID;
    int *list;
    bigval *sum;
    bigval *ssq;
} LOOP_PRINT;   

typedef struct {
    int ID;                      /* ID number for model */
    int ci;                      /* command index for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by ID number */
    int ifc;                     /* = 1 if the equation includes a constant,
                                    else = 0 */
    bigval *sum_coeff;      /* sums of coefficient estimates */
    bigval *ssq_coeff;      /* sums of squares of coeff estimates */
    bigval *sum_sderr;      /* sums of estimated std. errors */
    bigval *ssq_sderr;      /* sums of squares of estd std. errs */
} LOOP_MODEL;

struct LOOPSET_ {
    char type;
    int level;
    int err;
    int ntimes;
    int index;
    int index_start;
    int lvar;
    int rvar;
    double rval;
    int ineq;
    int ncmds;
    int nmod;
    int nprn;
    int nstore;
    int next_model;
    int next_print;
    char **lines;
    int *ci;
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
static int print_loop_store (LOOPSET *loop, PRN *prn, PATHS *ppaths);
static int get_prnnum_by_id (LOOPSET *loop, int id);
static int get_modnum_by_cmdnum (LOOPSET *loop, int cmdnum);

#define LOOP_BLOCK 32

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
    if (ci == OLS || 
	ci == GENR ||
	ci == LOOP ||
	ci == STORE ||
	ci == PRINT ||
	ci == PRINTF ||
	ci == PVALUE ||
	ci == SIM ||
	ci == SMPL ||
	ci == SUMMARY ||
	ci == IF ||
	ci == ELSE ||
	ci == ENDIF ||
	ci == ENDLOOP) { 
	return 1;
    }

    if (loop->type == COUNT_LOOP && 
	(ci == LAD || ci == HSK || ci == HCCM || ci == WLS)) {
	return 1;
    }

    return 0;
}

/* ......................................................  */

LOOPSET *loop_get_parent (LOOPSET *loop)
{
    if (loop != NULL) {
	return loop->parent;
    } else {
	return NULL;
    }
}

/* ......................................................  */

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

#ifdef LOOP_DEBUG
    printf("child loop %p has parent %p\n", 
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

#ifdef LOOP_DEBUG
	printf("gretl_loop_new: loop = %p, parent = %p\n", 
	       (void *) loop, (void *) parent);
#endif	

	if (err) {
	    free(loop);
	    loop = NULL;
	} 
    }
	
    return loop;
}

static int parse_as_while_loop (LOOPSET *loop,
				const DATAINFO *pdinfo,
				char *lvar, char *rvar, 
				char *op)
{
    int v;
    int err = 0;

    if (strstr(op, ">")) {
	loop->ineq = GT;
    } else {
	loop->ineq = LT;
    }

    v = varindex(pdinfo, lvar);
    if (v > 0 && v < pdinfo->v) {
	loop->lvar = v;
    } else {
	sprintf(gretl_errmsg, 
		_("Undefined variable '%s' in loop condition."), lvar);
	err = 1;
    }

    if (!err && (isdigit((unsigned char) *rvar) || *rvar == '.')) { 
	/* numeric rvalue?  FIXME: too restrictive?  */
	if (check_atof(rvar)) {
	    err = 1;
	} else {
	    loop->rval = atof(rvar);
	}
    } else if (!err) { /* otherwise try a varname */
	v = varindex(pdinfo, rvar);
	if (v > 0 && v < pdinfo->v) {
	    loop->rvar = v;
	} else {
	    sprintf(gretl_errmsg, 
		    _("Undefined variable '%s' in loop condition."), rvar);
	    loop->lvar = 0;
	    err = 1;
	}
    }

    if (!err) {
	loop->type = WHILE_LOOP;
    }

    return err;
}

static int parse_as_indexed_loop (LOOPSET *loop,
				  char *lvar, int start,
				  int end)
{
    int err = 0;

    if (strcmp(lvar, "i")) {
	sprintf(gretl_errmsg, 
		_("The index variable in a 'for' loop must be the "
		  "special variable 'i'"));
	err = 1;
    }
    if (!err && end <= start) {
	sprintf(gretl_errmsg, _("Ending value for loop index must be greater "
				"than starting value."));
	err = 1;
    }
    if (!err) {
	/* initialize loop index to starting value */
	loop->index_start = start;
	loop->lvar = 0;
	loop->rvar = 0;
	loop->ntimes = end;
	loop->type = INDEX_LOOP;
    }

    return err;
}

static int parse_as_simple_count_loop (LOOPSET *loop, int n)
{
    int err = 0;

    if (n <= 0) {
	strcpy(gretl_errmsg, _("Loop count must be positive."));
	err = 1;
    } else {
	loop->ntimes = n;
	loop->type = COUNT_LOOP;
    }

    return err;
}

static int parse_as_named_count_loop (LOOPSET *loop, 
				      const DATAINFO *pdinfo,
				      const double **Z,
				      const char *lvar)
{
    int v = varindex(pdinfo, lvar);
    int err = 0;

    if (v > 0 && v < pdinfo->v && pdinfo->vector[v] == 0) {
	int n = Z[v][0];

	if (n <= 0) {
	    strcpy(gretl_errmsg, _("Loop count must be positive."));
	    err = 1;
	} else {
	    loop->ntimes = n;
	    loop->type = COUNT_LOOP;
	}
    }

    return err;
}

/**
 * parse_loopline:
 * @line: command line.
 * @ploop: current loop struct pointer, or %NULL.
 * @loopstack: stacking level for the loop.
 * @pdinfo: data information struct.
 * @Z: data array.
 *
 * Parse a line specifying a loop condition.
 *
 * Returns: loop pointer on successful completion, %NULL on error.
 */

static LOOPSET *
parse_loopline (char *line, LOOPSET *ploop, int loopstack,
		DATAINFO *pdinfo, const double **Z)
{
    LOOPSET *loop;
    char lvar[VNAMELEN], rvar[VNAMELEN], op[8];
    int start, end;
    int n, err = 0;

#ifdef LOOP_DEBUG
    printf("parse_loopline: ploop = %p, loopstack = %d\n",
	   (void *) ploop, loopstack);
#endif

    if (ploop == NULL) {
	/* starting from scratch */
#ifdef LOOP_DEBUG
	printf("parse_loopline: starting from scratch\n");
#endif
	loop = gretl_loop_new(NULL, 0);
	if (loop == NULL) {
	    return NULL;
	}
    } else if (loopstack > ploop->level) {
	/* have to nest this loop */
#ifdef LOOP_DEBUG
	printf("parse_loopline: adding child\n");
#endif
	loop = gretl_loop_new(ploop, loopstack);
	if (loop == NULL) {
	    gretl_errmsg_set(_("Out of memory!"));
	    return NULL;
	}
    } else {
	/* shouldn't happen: need error message? */
	loop = ploop;
    }

    *gretl_errmsg = '\0';
    
    err = prepare_loop_for_action(loop);
    if (err) {
	goto bailout;
    }

    if (sscanf(line, "loop while %[^ <>]%[ <>] %s", lvar, op, rvar) == 3) {
	err = parse_as_while_loop(loop, pdinfo, lvar, rvar, op);
    }

    else if (sscanf(line, "loop for %[^= ] = %d..%d", lvar, &start, &end) == 3) {
	err = parse_as_indexed_loop(loop, lvar, start, end);
    }

    else if (sscanf(line, "loop %d", &n) == 1) {
	err = parse_as_simple_count_loop(loop, n);
    }

    else if (sscanf(line, "loop %8s", lvar) == 1) {
	err = parse_as_named_count_loop(loop, pdinfo, Z, lvar);
    }

    /* out of options, complain */
    else {
	strcpy(gretl_errmsg, _("No valid loop condition was given."));
	err = 1;
    }

    if (!err && loop->lvar == 0 && loop->ntimes < 2) {
	strcpy(gretl_errmsg, _("Loop count missing or invalid\n"));
	err = 1;
    }

 bailout:

    if (err) {
	if (loop != ploop) {
	    free(loop);
	}
	loop = NULL;
    }

    return loop;
}

#define DEFAULT_MAX_ITER 250

static int get_maxloop (void)
{
    static int ml = 0;

    if (ml == 0) {
	char *mlstr = getenv("GRETL_MAX_ITER");

	if (mlstr != NULL && sscanf(mlstr, "%d", &ml)) ;
	else ml = DEFAULT_MAX_ITER;
    }

    return ml;
}

/**
 * loop_condition:
 * @k: in case of a simple count loop, the number of iterations so far.
 * @loop: pointer to loop commands struct.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Check whether a looping condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

static int 
loop_condition (int k, LOOPSET *loop, double **Z, DATAINFO *pdinfo)
{
    static int maxloop = 0;
    int cont = 0;
    int oldtimes = loop->ntimes;

    if (maxloop == 0) {
	maxloop = get_maxloop();
    }

    /* an inequality between variables */
    if (loop->rvar > 0) {
	loop->ntimes += 1;
	if (loop->ntimes >= maxloop) {
	    sprintf(gretl_errmsg, _("Warning: no convergence after %d interations"),
		    maxloop);
	    loop->err = 1;
	    return 0; /* safety measure */
	}
	if (loop->ineq == GT) {
	    cont = (Z[loop->lvar][0] > Z[loop->rvar][0]);
	} else {
	    cont = (Z[loop->lvar][0] < Z[loop->rvar][0]);
	}
    } 

    /* a 'for' indexed loop */
    else if (loop->type == INDEX_LOOP) {  
	if (loop->index <= loop->ntimes) {
	    loop->index += 1;
	    cont = 1;
	}
    }

    /* inequality between a var and a number */
    else if (loop->lvar) {
	loop->ntimes += 1;
	if (loop->ntimes >= maxloop) {
	    sprintf(gretl_errmsg, _("Warning: no convergence after %d interations"),
		    maxloop);
	    loop->err = 1;
	    return 0; /* safety measure */
	}
	if (loop->ineq == GT) {
	    cont = (Z[loop->lvar][0] > loop->rval);
	} else {
	    cont = (Z[loop->lvar][0] < loop->rval);
	}
    }

    /* a simple count loop */
    else {
	if (k < loop->ntimes) cont = 1;
    }

    if (!cont && oldtimes == 0) {
	strcpy(gretl_errmsg, _("Loop condition not satisfied at first round"));
	loop->err = 1;
	loop->ntimes = 0;
    }

    return cont;
}

/* ......................................................  */

static void gretl_loop_init (LOOPSET *loop)
{
#ifdef LOOP_DEBUG
    printf("gretl_loop_init: initing loop at %p\n", (void *) loop);
#endif

    loop->level = 0;

    loop->ntimes = 0;
    loop->index = 0;
    loop->index_start = 0;
    loop->err = 0;
    loop->lvar = 0;
    loop->rvar = 0;
    loop->rval = 0;
    loop->ineq = 0;

    loop->ncmds = 0;
    loop->nmod = 0;
    loop->nprn = 0;
    loop->nstore = 0;

    loop->next_model = 0;
    loop->next_print = 0;

    loop->lines = NULL;
    loop->ci = NULL;
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

    if (loop->lines) {
	for (i=0; i<loop->ncmds; i++)
	    if (loop->lines[i]) free (loop->lines[i]);
	free(loop->lines);
	loop->lines = NULL;
    }

    if (loop->ci) { 
	free(loop->ci);
    }

    if (loop->lmodels) {
	for (i=0; i<loop->nmod; i++) {
	    free_loop_model(&loop->lmodels[i]);
	}
	free(loop->lmodels);
	loop->lmodels = NULL;
    }

    if (loop->models) {
	for (i=0; i<loop->nmod; i++) {
	    free_model(loop->models[i]);
	    loop->models[i] = NULL;
	}
	free(loop->models);
	loop->models = NULL;
    } 
   
    if (loop->prns) {
	for (i=0; i<loop->nprn; i++) { 
	    free_loop_print(&loop->prns[i]);
	}
	free(loop->prns);
	loop->prns = NULL;
    }

    if (loop->storename) {
	for (i=0; i<loop->nstore; i++) {
	    if (loop->storename[i]) {
		free(loop->storename[i]);
	    }
	}
	free(loop->storename);
	loop->storename = NULL;
    }

    if (loop->storelbl) {
	for (i=0; i<loop->nstore; i++) {
	    if (loop->storelbl[i]) {
		free(loop->storelbl[i]);
	    }
	}
	free(loop->storelbl);
	loop->storelbl = NULL;
    }

    if (loop->storeval) { 
	free(loop->storeval);
	loop->storeval = NULL;
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
    int i, ncoeff = pmod->ncoeff;

    lmod->sum_coeff = malloc(ncoeff * sizeof *lmod->sum_coeff);
    if (lmod->sum_coeff == NULL) return 1;

    lmod->ssq_coeff = malloc(ncoeff * sizeof *lmod->ssq_coeff);
    if (lmod->ssq_coeff == NULL) goto cleanup;

    lmod->sum_sderr = malloc(ncoeff * sizeof *lmod->sum_sderr);
    if (lmod->sum_sderr == NULL) goto cleanup;

    lmod->ssq_sderr = malloc(ncoeff * sizeof *lmod->ssq_sderr);
    if (lmod->ssq_sderr == NULL) goto cleanup;

    lmod->list = copylist(pmod->list);
    if (lmod->list == NULL) goto cleanup;

    for (i=0; i<ncoeff; i++) {
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

    lmod->ncoeff = ncoeff;
    lmod->t1 = pmod->t1;
    lmod->t2 = pmod->t2;
    lmod->nobs = pmod->nobs;
    lmod->ID = id;
    lmod->ci = pmod->ci;

    return 0;

 cleanup:
    free(lmod->ssq_coeff);
    free(lmod->sum_sderr);
    free(lmod->ssq_sderr);

    return 1;
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

static int loop_print_init (LOOP_PRINT *lprn, const LIST list, int id)
{
    int i;

    lprn->list = copylist(list);
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
			    const LIST list, DATAINFO *pdinfo)
{
    int i, tot = list[0] * loop->ntimes;

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

static int add_loop_model (LOOPSET *loop, int cmdnum)
{
    int err = 0;
    int nm = loop->nmod + 1;

    if (loop->type != COUNT_LOOP) { 
	/* a conditional loop */
	loop->models = realloc(loop->models, nm * sizeof(MODEL *));
	if (loop->models == NULL) {
	    err = 1;
	} else {
	    loop->models[loop->nmod] = gretl_model_new();
	    if (loop->models[loop->nmod] == NULL) {
		err = 1;
	    } else {
		(loop->models[loop->nmod])->ID = cmdnum;
	    }
	}
    } else { 
	/* looping a fixed number of times */
	loop->lmodels = realloc(loop->lmodels, nm * sizeof *loop->lmodels);
	if (loop->lmodels == NULL) {
	    err = 1;
	}
    }

    if (!err) {
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

static int add_loop_print (LOOPSET *loop, const LIST list, int cmdnum)
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
			      const LIST list, double ***pZ, 
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
 * @ppaths: path information struct.
 *
 * Print out the results after completion of the loop @loop.
 *
 */

static void print_loop_results (LOOPSET *loop, const DATAINFO *pdinfo, 
				PRN *prn, PATHS *ppaths)
{
    int i, j;
    gretlopt opt;
    MODEL *pmod = NULL;

    if (loop->lvar && loop->type == INDEX_LOOP) {
	pprintf(prn, _("\nNumber of iterations: %d\n\n"), loop->ntimes);
    }

    for (i=0; i<loop->ncmds; i++) {
#ifdef LOOP_DEBUG
	fprintf(stderr, "loop command %d (i=%d): %s\n\n", i+1, i, loop->lines[i]);
#endif
	catchflags(loop->lines[i], &opt);

	if (loop->lvar && loop->ci[i] == OLS) {
	    double dfadj, sqrta;

	    pmod = loop->models[loop->next_model];

	    set_model_id(pmod);

	    /* std. errors are asymptotic; degrees of freedom
	       correction is not wanted */
	    dfadj = (double) pmod->dfd / pmod->nobs;
	    sqrta = sqrt(dfadj);
	    pmod->sigma = sqrt((1.0 / pmod->nobs) * pmod->ess);
	    for (j=0; j<pmod->ncoeff; j++) {
		pmod->sderr[j] *= sqrta;
	    }

	    printmodel(pmod, pdinfo, opt, prn);

#if 0
	    if (opt & OPT_O) { /* FIXME */
		if (pmod->vcv) {
		    int nc = pmod->ncoeff;
		    int nt = nc * (nc + 1) / 2;

		    for (j=0; j<nt; j++) {
			pmod->vcv[j] *= dfadj;
		    }
		} else {
		    makevcv(pmod);
		}
		outcovmx(pmod, pdinfo, prn);
	    }
#endif

	    loop->next_model += 1;	    
	}
	else if (loop->ci[i] == OLS || loop->ci[i] == LAD ||
		 loop->ci[i] == HSK || loop->ci[i] == HCCM || 
		 loop->ci[i] == WLS) {
	    print_loop_model(&loop->lmodels[loop->next_model], 
			     loop->ntimes, pdinfo, prn);
	    loop->next_model += 1;
	}
	else if (loop->ci[i] == PRINT) {
	    print_loop_prn(&loop->prns[loop->next_print], 
			   loop->ntimes, pdinfo, prn);
	    loop->next_print += 1;
	}
	else if (loop->ci[i] == STORE) {
	    print_loop_store(loop, prn, ppaths);
	}
    }
}

/* ......................................................  */

static void free_loop_model (LOOP_MODEL *lmod)
{
#ifdef ENABLE_GMP
    int i;

    for (i=0; i<lmod->ncoeff; i++) {
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
    free(lmod->list);
}

/* ......................................................  */

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
 * @oflags: option flag(s) associated with the command.
 * @pdinfo: dataset information.
 * @Z: data matrix.
 * @loopstack: pointer to integer stacking level.
 * @pointer to integer switch for running loop.
 *
 * Add line and command index to accumulated loop buffer.
 *
 * Returns: pointer to loop struct on success, %NULL on failure.
 */

LOOPSET *add_to_loop (char *line, int ci, gretlopt opt,
		      DATAINFO *pdinfo, const double **Z,
		      LOOPSET *loop, int *loopstack, int *looprun)
{
    LOOPSET *lret = loop;

#ifdef LOOP_DEBUG
    printf("add_to_loop: loop = %p, loopstack = %d, line = '%s'\n", 
	   (void *) loop, *loopstack, line);
#endif

    if (ci == LOOP) {
	lret = parse_loopline(line, loop, *loopstack, pdinfo, Z);
	if (lret == NULL) {
	    gretl_errmsg_set("error: loop is NULL");
	} else {
	    *loopstack += 1;
#ifdef LOOP_DEBUG
	    printf("after LOOP: loop = %p, parent = %p, loopstack = %d\n", 
		   (void *) lret, (void *) lret->parent, *loopstack);
#endif
	}
    } else if (ci == ENDLOOP) {
	*loopstack -= 1;
	if (*loopstack == 0) {
	    *looprun = 1;
	} else {
	    lret = loop->parent;
#ifdef LOOP_DEBUG
	    printf("after ENDLOOP: loop = %p, parent = %p, loopstack = %d\n", 
		   (void *) lret, (void *) loop_get_parent(lret), *loopstack);
#endif
	}
    } 

    if (loop != NULL) {
	int nc = loop->ncmds;

	if ((nc + 1) % LOOP_BLOCK == 0) {
	    if (add_more_loop_lines(loop)) {
		gretl_loop_destroy(loop);
		loop = NULL;
		goto bailout;
	    }
	}

	loop->lines[nc] = malloc(MAXLEN);
	if (loop->lines[nc] == NULL) {
	    gretl_loop_destroy(loop);
	    loop = NULL;
	    goto bailout;
	}

	top_n_tail(line);

	if (ci == PRINT && loop->type != COUNT_LOOP) {
	    loop->ci[nc] = 0;
	} else {
	    loop->ci[nc] = ci;
	}

	loop->lines[nc][0] = '\0';
	/* FIXME magic number 24 */
	strncat(loop->lines[nc], line, MAXLEN - 24);

	loop->ncmds += 1;

	if (opt) {
	    const char *flagstr = print_flags(opt, ci);

	    strcat(loop->lines[nc], flagstr);
	}
    }

 bailout:

    return lret;
}

/* ......................................................... */ 

static void print_loop_model (LOOP_MODEL *lmod, int loopnum,
			      const DATAINFO *pdinfo, PRN *prn)
{
    int i;
    char startdate[OBSLEN], enddate[OBSLEN];

    ntodate(startdate, lmod->t1, pdinfo);
    ntodate(enddate, lmod->t2, pdinfo);

    pprintf(prn, _("%s estimates using the %d observations %s-%s\n"),
	    _(estimator_string(lmod->ci, prn->format)), lmod->t2 - lmod->t1 + 1, 
	    startdate, enddate);
    pprintf(prn, _("Statistics for %d repetitions\n"), loopnum); 
    pprintf(prn, _("Dependent variable: %s\n\n"), 
	    pdinfo->varname[lmod->list[1]]);

    pputs(prn, _("                     mean of      std. dev. of     mean of"
		 "     std. dev. of\n"
		 "                    estimated      estimated"
		 "      estimated      estimated\n"
		 "      Variable     coefficients   coefficients   std. errors"
		 "    std. errors\n\n"));

    for (i=0; i<lmod->ncoeff; i++) {
	print_loop_coeff(pdinfo, lmod, i, loopnum, prn);
    }
    pputs(prn, "\n");
}

/* ......................................................... */ 

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

    pprintf(prn, " %3d) %8s ", lmod->list[c+2], 
	   pdinfo->varname[lmod->list[c+2]]);

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

    pprintf(prn, " %3d) %8s ", lmod->list[c+2], 
	   pdinfo->varname[lmod->list[c+2]]);

    pprintf(prn, "%#14g %#14g %#14g %#14g\n", (double) m1, (double) sd1, 
	    (double) m2, (double) sd2);
#endif
}

/* ......................................................... */ 

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
    pputs(prn, "\n");
}

/* ......................................................... */ 

static int print_loop_store (LOOPSET *loop, PRN *prn, PATHS *ppaths)
{
    int i, t;
    FILE *fp;
    char gdtfile[MAXLEN], infobuf[1024];
    char *xmlbuf = NULL;
    time_t writetime;

    /* organize filename */
    if (loop->storefile[0] == '\0') {
	sprintf(gdtfile, "%sloopdata.gdt", ppaths->userdir);	
    } else if (slashpos(loop->storefile) == 0) { 
	/* no path given (FIXME) */
	sprintf(gdtfile, "%s%s", ppaths->userdir, loop->storefile);
    } else {
	strcpy(gdtfile, loop->storefile);
    }

    if (strchr(gdtfile, '.') == NULL) {
	strcat(gdtfile, ".gdt");
    }

    fp = fopen(gdtfile, "w");
    if (fp == NULL) return 1;

    writetime = time(NULL);

    pprintf(prn, _("printing %d values of variables to %s\n"), 
	    loop->ntimes, gdtfile);

    fprintf(fp, "<?xml version=\"1.0\"?>\n"
	    "<!DOCTYPE gretldata SYSTEM \"gretldata.dtd\">\n\n"
	    "<gretldata name=\"%s\" frequency=\"1\" "
	    "startobs=\"1\" endobs=\"%d\" ", 
	    gdtfile, loop->ntimes);

    fprintf(fp, "type=\"cross-section\">\n");

    sprintf(infobuf, "%s %s", _("simulation data written"),
	    print_time(&writetime)); 
    xmlbuf = gretl_xml_encode(infobuf);
    fprintf(fp, "<description>\n%s\n</description>\n", xmlbuf);
    free(xmlbuf);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

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
	    loop->ntimes);

    for (t=0; t<loop->ntimes; t++) {
	double x;

	fputs("<obs>", fp);
	for (i=0; i<loop->nstore; i++) {
	    x = loop->storeval[loop->ntimes*i + t];
	    if (na(x)) {
		fputs("NA ", fp);
	    } else {
		fprintf(fp, "%g ", x);
	    }
	}
	fputs("</obs>\n", fp);
    }

    fprintf(fp, "</observations>\n</gretldata>\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);
    return 0;
}

/* ......................................................... */ 

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

    for (i=0; i<loop->nmod; i++) {
	if (loop->lvar != 0 && (loop->models[i])->ID == cmdnum) {
	    return i;
	}
	if (loop->lvar == 0 && loop->lmodels[i].ID == cmdnum) {
	    return i;
	}
    }

    return -1;
}

/**
 * get_cmd_ci:
 * @line: command line.
 * @command: pointer to gretl command struct.
 *
 * Parse @line and assign to @command->ci the index number of
 * the command embedded in @line.
 */

void get_cmd_ci (const char *line, CMD *command)
{
    /* allow for leading spaces */
    while (isspace(*line)) line++;

    if (sscanf(line, "%s", command->cmd) != 1 || 
	*line == '(' || *line == '#') {
	command->nolist = 1;
	command->ci = -1;
	return;
    }
    if ((command->ci = gretl_command_number(command->cmd)) == 0) {
	command->errcode = 1;
	sprintf(gretl_errmsg, _("command \"%s\" not recognized"), 
		command->cmd);
	return;
    }    
} 

static void substitute_dollar_i (char *str, int index)
{
    char *p;

    while ((p = strstr(str, "$i")) != NULL) {
	char ins[8];
	char *q;

	q = malloc(strlen(p));
	strcpy(q, p + 2);
	sprintf(ins, "%d", index - 1);
	strcpy(p, ins);
	strcpy(p + strlen(ins), q);
	free(q);	
    }
}

int loop_exec (LOOPSET *loop, char *line,
	       double ***pZ, DATAINFO **ppdinfo, 
	       MODEL **models, PATHS *paths, 
	       int *echo_off, PRN *prn)
{
    DATAINFO *pdinfo = *ppdinfo;
    CMD cmd;
    int lround = 0;
    int ignore = 0;
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

    gretl_set_text_pause(0);

#ifdef LOOP_DEBUG
    printf("loop_exec: loop = %p\n", (void *) loop);
#endif

    loop->index = loop->index_start;

    while (!err && loop_condition(lround, loop, *pZ, pdinfo)) {
	int childnum = 0;
	int j;

#ifdef LOOP_DEBUG
	printf("top of loop: lround = %d\n", lround);
#endif

	if (loop->type == INDEX_LOOP && !(*echo_off)) {
	    pprintf(prn, "loop: i = %d\n\n", loop->index - 1);
	}

	for (j=0; !err && j<loop->ncmds; j++) {
	    char linecpy[MAXLINE];
	    static MODEL *tmpmodel;

#ifdef LOOP_DEBUG
	    printf("loop->lines[%d] = '%s'\n", j, loop->lines[j]);
#endif
	    strcpy(linecpy, loop->lines[j]);

	    err = catchflags(linecpy, &cmd.opt);
	    if (err) {
		break;
	    }

	    substitute_dollar_i(linecpy, loop->index);

	    getcmd(linecpy, pdinfo, &cmd, &ignore, pZ, NULL);

	    if (cmd.ci < 0) continue;

	    if (cmd.errcode) {
		err = cmd.errcode;
		break;
	    }

	    if (!(*echo_off) && loop->type == INDEX_LOOP) {
		echo_cmd(&cmd, pdinfo, linecpy, 0, 1, prn);
	    }

	    switch (cmd.ci) {

	    case LOOP:
		err = loop_exec(loop->children[childnum++], NULL,
				pZ, ppdinfo, models, paths,
				echo_off, prn);
		break;

	    case ENDLOOP:
		/* no-op */
		break;

	    case GENR:
		err = generate(pZ, pdinfo, linecpy, tmpmodel);
		break;

	    case SIM:
		err = simulate(linecpy, pZ, pdinfo);
		break;	

	    case OLS:
	    case WLS:
	    case LAD:
	    case HSK:
	    case HCCM:
		/* if this is the first time round the loop, allocate space
		   for each loop model */
		if (lround == 0 && loop->type != INDEX_LOOP) {
		    err = add_loop_model(loop, j);
		    if (err) {
			break;
		    }
		} /* end of basic round 0 setup */

		/* estimate the model called for */
		clear_model(models[0]);

		if (cmd.ci == OLS || cmd.ci == WLS) {
		    *models[0] = lsq(cmd.list, pZ, pdinfo, cmd.ci, cmd.opt, 0.0);
		} else if (cmd.ci == LAD) {
		    *models[0] = lad(cmd.list, pZ, pdinfo);
		} else if (cmd.ci == HSK) {
		    *models[0] = hsk_func(cmd.list, pZ, pdinfo);
		} else if (cmd.ci == HCCM) {
		    *models[0] = hccm_func(cmd.list, pZ, pdinfo);
		}

		if ((err = (models[0])->errcode)) {
		    break;
		}

		if (loop->type == INDEX_LOOP) {
		    (models[0])->ID = lround + 1;
		    printmodel(models[0], pdinfo, cmd.opt, prn); 
		    tmpmodel = models[0];
		} else if (loop->type != COUNT_LOOP) { /* conditional loop */
		    /* deal with model estimate for "while" loop */
		    int m = get_modnum_by_cmdnum(loop, j);

		    swap_models(&models[0], &loop->models[m]);
		    (loop->models[m])->ID = j;
		    tmpmodel = loop->models[m];
		    model_count_minus();
		} else { 
		    /* looping a fixed number of times */
		    if (lround == 0 && loop_model_init(&loop->lmodels[loop->nmod - 1], 
						       models[0], j)) { 
			gretl_errmsg_set(_("Failed to initialize model for loop\n"));
			err = 1;
			break;
		    } else if (update_loop_model(loop, j, models[0])) {
			gretl_errmsg_set(_("Failed to add results to loop model\n"));
			err = 1;
			break;
		    }
		    tmpmodel = models[0];
		}
		break;

	    case PRINT:
		if (cmd.param[0] != '\0') {
		    simple_commands(&cmd, linecpy, pZ, pdinfo, paths, prn);
		    break;
		}
		if (loop->type != COUNT_LOOP) {
		    printdata(cmd.list, pZ, pdinfo, cmd.opt, prn);
		    break;
		}
		if (lround == 0) {
		    if ((err = add_loop_print(loop, cmd.list, j))) {
			break;
		    }
		}
		if (update_loop_print(loop, j, cmd.list, pZ, pdinfo)) {
		    gretl_errmsg_set(_("Failed to add values to print loop\n"));
		    err = 1;
		}
		break;

	    case PRINTF:
		err = do_printf(linecpy, pZ, pdinfo, models[0], prn);
		break;

	    case SMPL:
		if (cmd.opt) {
		    if (restore_full_sample(pZ, ppdinfo, cmd.opt)) {
			err = 1;
			break;
		    }
		    if (restrict_sample(linecpy, pZ, ppdinfo, 
					cmd.list, cmd.opt)) {
			err = 1;
			break;
		    }
		    pdinfo = *ppdinfo;
		} else {
		    /* FIXME */
		    gretl_errmsg_set(_("loop: only the '-o' and '-r' forms of the smpl "
				       " command may be used.\n"));  
		    err = 1;
		}
		break;

	    case STORE:
		if (lround == 0) {
		    loop->nstore = cmd.list[0];
		    if (loop_store_init(loop, cmd.param, cmd.list, pdinfo)) {
			err = 1;
		    }
		} else {
		    int i;

		    for (i=0; i<cmd.list[0]; i++) {
			if (pdinfo->vector[cmd.list[i+1]]) { 
			    loop->storeval[i * loop->ntimes + lround] = 
				(*pZ)[cmd.list[i+1]][pdinfo->t1 + 1];
			} else {
			    loop->storeval[i * loop->ntimes + lround] = 
				(*pZ)[cmd.list[i+1]][0];
			}
		    }
		}	
		break;

	    case PVALUE:
		batch_pvalue(loop->lines[j], *pZ, pdinfo, prn);
		break;

	    case SUMMARY:
		if (loop->type == COUNT_LOOP) {
		    gretl_errmsg_set( _("The summary command is not available in "
					"this sort of loop.\n"));
		    err = 1;
		} else {
		    GRETLSUMMARY *summ;

		    summ = summary(cmd.list, pZ, pdinfo, prn);
		    if (summ == NULL) {
			gretl_errmsg_set(_("generation of summary stats failed\n"));
			err = 1;
		    } else {
			print_summary(summ, pdinfo, prn);
			free_summary(summ);
		    }
		}	    
		break; 

	    default: 
		/* not reachable */
		pprintf(prn, _("command: '%s'\nThis is not available in a loop.\n"),
			linecpy);
		err = 1;
		break;

	    } /* end switch on command number */
	} /* end list of commands within loop */

	lround++;

    } /* end iterations of loop */

    if (err) {
	print_gretl_errmsg(prn);
    } else if (loop->err) {
	print_gretl_errmsg(prn);
	err = loop->err;
    }

    if (!err && lround > 0) {
	if (loop->type != INDEX_LOOP) {
	    print_loop_results(loop, pdinfo, prn, paths); 
	}
    }

    if (line != NULL) {
	*line = '\0';
    }

    gretl_cmd_free(&cmd);

    return err;
}

/* ifthen stuff - conditional execution */

int if_eval (const char *line, double ***pZ, DATAINFO *pdinfo)
{
    char formula[MAXLEN];
    int err, ret = -1;

    /* + 2 below to omit "if" */
    sprintf(formula, "__iftest=%s", line + 2);
    err = generate(pZ, pdinfo, formula, NULL);
    if (!err) {
	int v = varindex(pdinfo, "iftest");
	
	if (v < pdinfo->v) {
	    ret = (*pZ)[v][0];
	    dataset_drop_vars(1, pZ, pdinfo);
	}
    }

    return ret;
}

int ifstate (int code)
{
    static unsigned char T[9];
    static unsigned char got_if[9];
    static unsigned char got_else[9];
    static unsigned char indent;

    if (code == RELAX) {
	indent = 0;
    }
    else if (code == SET_FALSE || code == SET_TRUE) {
	indent++;
	if (indent > 8) {
	    return 1; /* too deeply nested */
	}
	T[indent] = (code == SET_TRUE);
	got_if[indent] = 1;
	got_else[indent] = 0;
    }
    else if (code == SET_ELSE) {
	if (got_else[indent] || !got_if[indent]) {
	    sprintf(gretl_errmsg, "Unmatched \"else\"");
	    return 1; 
	}
	T[indent] = !T[indent];
	got_else[indent] = 1;
    }
    else if (code == SET_ENDIF) {
	if (!got_if[indent] || indent == 0) {
	    sprintf(gretl_errmsg, "Unmatched \"endif\"");
	    return 1; 
	}
	got_if[indent] = 0;
	got_else[indent] = 0;
	indent--;
    }
    else if (code == IS_FALSE) {
	int i;

	for (i=1; i<=indent; i++) {
	    if (T[i] == 0) return 1;
	}
    }

    return 0;
}


