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

#include <time.h>
#include <unistd.h>

enum inequalities {
    GT = 1,
    LT
};

static void gretl_loop_init (LOOPSET *loop);
static int monte_carlo_init (LOOPSET *loop);
static void monte_carlo_free (LOOPSET *loop);
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
	ci == ENDLOOP) 
	return 1;

    if (loop->type == COUNT_LOOP && 
	(ci == LAD || ci == HSK || ci == HCCM || ci == WLS)) 
	return 1;

    return 0;
}

/* ......................................................  */

static int loop_attach_child (LOOPSET *loop, LOOPSET *child)
{
    LOOPSET **children;
    int n;

    n = loop->n_children + 1;
    children = realloc(loop->children, n * sizeof *loop->children);
    if (children == NULL) {
	return 1;
    } 

    loop->children = children;
    loop->children[n - 1] = child;
    child->parent = loop;
    child->level = loop->level + 1;

    loop->n_children += 1;

    return 0;
}

static LOOPSET *gretl_loop_new (LOOPSET *ploop, int loopstack)
{
    LOOPSET *loop;

    loop = malloc(sizeof *loop);
    if (loop == NULL) {
	return NULL;
    }

    gretl_loop_init(loop);

    if (ploop != NULL) {
	int err = loop_attach_child(ploop, loop);

	if (err) {
	    free(loop);
	    loop = NULL;
	} 
    }
	
    return loop;
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

LOOPSET *parse_loopline (char *line, LOOPSET *ploop, int loopstack,
			 DATAINFO *pdinfo, const double **Z)
{
    LOOPSET *loop;
    char lvar[VNAMELEN], rvar[VNAMELEN], op[8];
    int start, end;
    int n, v;
    int err = 0;

    if (ploop == NULL) {
	/* starting from scratch */
	loop = gretl_loop_new(NULL, 0);
	if (loop == NULL) {
	    return NULL;
	}
    } else if (loopstack > ploop->level) {
	/* have to nest this loop */
	loop = gretl_loop_new(ploop, loopstack);
	if (loop == NULL) {
	    return NULL;
	}
    } else {
	/* shouldn't happen: need error message? */
	loop = ploop;
    }

    *gretl_errmsg = '\0';
    monte_carlo_init(loop);

    /* try parsing as a while loop */
    if (sscanf(line, "loop while %[^ <>]%[ <>] %s", lvar, op, rvar) == 3) {
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
    }

    /* or try parsing as a for loop */
    else if (sscanf(line, "loop for %[^= ] = %d..%d", lvar, &start, &end) == 3) {
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
	    /* initialize special genr index to starting value */
	    genr_scalar_index(1, start - 1);
	    loop->lvar = INDEXNUM;
	    loop->rvar = 0;
	    loop->ntimes = end;
	    loop->type = FOR_LOOP;
	}
    }

    /* or as a simple count loop */
    else if (sscanf(line, "loop %d", &n) == 1) {
	if (n <= 0) {
	    strcpy(gretl_errmsg, _("Loop count must be positive."));
	    err = 1;
	} else {
	    loop->ntimes = n;
	    loop->type = COUNT_LOOP;
	}
    }

    /* or as a count loop with named scalar max */
    else if (sscanf(line, "loop %8s", lvar) == 1) {
	v = varindex(pdinfo, lvar);
	if (v > 0 && v < pdinfo->v && pdinfo->vector[v] == 0) {
	    n = Z[v][0];
	    if (n <= 0) {
		strcpy(gretl_errmsg, _("Loop count must be positive."));
		err = 1;
	    } else {
		loop->ntimes = n;
		loop->type = COUNT_LOOP;
	    }
	}
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

int loop_condition (int k, LOOPSET *loop, double **Z, DATAINFO *pdinfo)
{
    static int maxloop = 0;
    int t, cont = 0;
    int oldtimes = loop->ntimes;

    if (maxloop == 0) maxloop = get_maxloop();

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

    /* a 'for' loop */
    else if (loop->lvar == INDEXNUM) {  
	t = genr_scalar_index(0, 0); /* fetch index */
	if (t < loop->ntimes) {
	    genr_scalar_index(2, 1); /* increment index */
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

LOOPSET *gretl_loop_terminate (LOOPSET *loop)
{
    LOOPSET *nextloop = loop->parent;

    /* free allocated storage in loop */
    monte_carlo_free(loop);

    /* free loop struct pointer itself */
    free(loop);

    /* FIXME: do more */

    return nextloop;
}

/* ......................................................  */

static void gretl_loop_init (LOOPSET *loop)
{
    loop->level = 0;

    loop->ntimes = 0;
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

static int monte_carlo_init (LOOPSET *loop)
{
    gretl_loop_init(loop);

#ifdef ENABLE_GMP
    mpf_set_default_prec(256);
#endif

    loop->lines = malloc(32 * sizeof *loop->lines); 
    loop->ci = malloc(32 * sizeof *loop->ci);
    
    if (loop->lines == NULL || loop->ci == NULL) {
	return 1;
    }

    return 0;
}

static void monte_carlo_free (LOOPSET *loop)
{
    int i;

    if (loop->lines) {
	for (i=0; i<loop->ncmds; i++)
	    if (loop->lines[i]) free (loop->lines[i]);
	free(loop->lines);
	loop->lines = NULL;
    }
    if (loop->ci) 
	free(loop->ci);
    if (loop->lmodels) {
	for (i=0; i<loop->nmod; i++) 
	    free_loop_model(&loop->lmodels[i]);
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
	for (i=0; i<loop->nprn; i++) 
	    free_loop_print(&loop->prns[i]);
	free(loop->prns);
	loop->prns = NULL;
    }
    if (loop->storename) {
	for (i=0; i<loop->nstore; i++)
	    if (loop->storename[i])
		free(loop->storename[i]);
	free(loop->storename);
	loop->storename = NULL;
    }
    if (loop->storelbl) {
	for (i=0; i<loop->nstore; i++)
	    if (loop->storelbl[i])
		free(loop->storelbl[i]);
	free(loop->storelbl);
	loop->storelbl = NULL;
    }
    if (loop->storeval) { 
	free(loop->storeval);
	loop->storeval = NULL;
    }
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

int loop_model_init (LOOP_MODEL *lmod, const MODEL *pmod,
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

int loop_print_init (LOOP_PRINT *lprn, const LIST list, int id)
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

int loop_store_init (LOOPSET *loop, const char *fname, 
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

int add_loop_model (LOOPSET *loop, int cmdnum)
{
    int err = 0;
    int nm = loop->nmod + 1;

    loop->nmod += 1;

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

int update_loop_model (LOOPSET *loop, int cmdnum, MODEL *pmod)
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

int add_loop_print (LOOPSET *loop, const LIST list, int cmdnum)
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


int update_loop_print (LOOPSET *loop, int cmdnum, 
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

void print_loop_results (LOOPSET *loop, const DATAINFO *pdinfo, 
			 PRN *prn, PATHS *ppaths)
{
    int i, j;
    gretlopt opt;
    MODEL *pmod = NULL;

    if (loop->lvar && loop->lvar != INDEXNUM) {
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

/**
 * add_to_loop:
 * @loop: pointer to loop struct.
 * @line: command line.
 * @ci: command index number.
 * @oflags: option flag(s) associated with the command.
 *
 * Add line and command index to accumulated loop buffer.
 *
 * Returns: 0 on successful completion.
 */

int add_to_loop (LOOPSET *loop, char *line, int ci,
		 gretlopt oflags)
{
    int i = loop->ncmds;

    loop->ncmds += 1;

    loop->lines[i] = malloc(MAXLEN);
    if (loop->lines[i] == NULL) return E_ALLOC;

    top_n_tail(line);

    if (ci == PRINT && loop->type != COUNT_LOOP) {
	loop->ci[i] = 0;
    } else {
	loop->ci[i] = ci;
    }

    strncpy(loop->lines[i], line, MAXLEN - 4);

    if (oflags) {
	const char *flagstr = print_flags(oflags, ci);

	strcat(loop->lines[i], flagstr);
    }

    return 0;
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

int get_modnum_by_cmdnum (LOOPSET *loop, int cmdnum)
{
    int i;

    for (i=0; i<loop->nmod; i++) {
	if (loop->lvar && (loop->models[i])->ID == cmdnum) return i;
	if (loop->lvar == 0 && loop->lmodels[i].ID == cmdnum) return i;
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
	if (indent > 8) return 1; /* too deeply nested */
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

	for (i=1; i<=indent; i++) if (T[i] == 0) return 1;
    }
    return 0;
}


