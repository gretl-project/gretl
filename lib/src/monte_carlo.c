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
#include "internal.h"
#include <unistd.h>

enum inequalities {
    GT = 1,
    LT
};

extern int genr_scalar_index (int opt, int put);

static int monte_carlo_init (LOOPSET *ploop);
static void free_loop_model (LOOP_MODEL *plmod);
static void free_loop_print (LOOP_PRINT *pprn);
static void print_loop_model (LOOP_MODEL *plmod, int loopnum,
			       const DATAINFO *pdinfo, PRN *prn);
static void print_loop_coeff (const DATAINFO *pdinfo, const LOOP_MODEL *plmod, 
			       int c, int n, PRN *prn);
static void print_loop_prn (LOOP_PRINT *pprn, int n,
			     const DATAINFO *pdinfo, PRN *prn);
static int print_loop_store (LOOPSET *ploop, PRN *prn, PATHS *ppaths,
			      char *loopstorefile);
static int get_prnnum_by_id (LOOPSET *ploop, int id);

/**
 * ok_in_loop:
 * @ci: command index.
 *
 * Returns: 1 if the given command is acceptable inside a loop construct,
 * 0 otherwise.
 */

int ok_in_loop (int ci)
{
    if (ci == OLS || 
	ci == LAD ||
	ci == HSK ||
	ci == HCCM ||
	ci == GENR ||
	ci == STORE ||
	ci == PRINT ||
	ci == PVALUE ||
	ci == SIM ||
	ci == SMPL ||
	ci == SUMMARY ||
	ci == IF ||
	ci == ELSE ||
	ci == ENDIF ||
	ci == ENDLOOP) 
	return 1;
    return 0;
}

/**
 * parse_loopline:
 * @line: command line.
 * @ploop: loop commands struct.
 * @pdinfo: data information struct.
 *
 * Parse a line specifying a loop condition.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int parse_loopline (char *line, LOOPSET *ploop, DATAINFO *pdinfo)
{
    char lvar[9], rvar[9], op[8];
    int start, end;
    int n, v;
    int err = 0;

    gretl_errmsg[0] = '\0';
    monte_carlo_init(ploop);

    /* try parsing as a while loop */
    if (sscanf(line, "loop while %[^ <>]%[ <>] %s", lvar, op, rvar) == 3) {
	if (strstr(op, ">")) ploop->ineq = GT;
	else ploop->ineq = LT;
	v = varindex(pdinfo, lvar);
	if (v > 0 && v < pdinfo->v) ploop->lvar = v;
	else {
	    sprintf(gretl_errmsg, 
		    _("Undefined variable '%s' in loop condition."), lvar);
	    err = 1;
	}
	if (!err && (isdigit((unsigned char) rvar[0]) 
	    || rvar[0] == '.')) { /* numeric rvalue? */
	    ploop->rval = atof(rvar);
	} else if (!err) { /* otherwise try a varname */
	    v = varindex(pdinfo, rvar);
	    if (v > 0 && v < pdinfo->v) ploop->rvar = v;
	    else {
		sprintf(gretl_errmsg, 
			_("Undefined variable '%s' in loop condition."), rvar);
		ploop->lvar = 0;
		err = 1;
	    }
	}
	ploop->type = WHILE_LOOP;
	return err;
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
	    ploop->lvar = INDEXNUM;
	    ploop->rvar = 0;
	    ploop->ntimes = end;
	    ploop->type = FOR_LOOP;
	}
	return err;
    }

    /* or just as a simple count loop */
    else if (sscanf(line, "loop %d", &n) == 1) {
	if (n <= 0) {
	    strcpy(gretl_errmsg, _("Loop count must be positive."));
	    err = 1;
	} else {
	    ploop->ntimes = n;
	    ploop->type = COUNT_LOOP;
	}
	return err;
    }

    /* out of options, complain */
    strcpy(gretl_errmsg, _("No valid loop condition was given."));
    return 1;
}

/**
 * loop_condition:
 * @k: in case of a simple count loop, the number of iterations so far.
 * @ploop: loop commands struct.
 * @Z: data matrix.
 * @pdinfo: data information struct.
 *
 * Check whether a looping condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

int loop_condition (int k, LOOPSET *ploop, double **Z, DATAINFO *pdinfo)
{
    int t;
    const int LOOPMAX = 80000; /* bodge: safety measure */

    if (ploop->lvar && ploop->ntimes > LOOPMAX) 
	return 0;

    /* case of an inequality between variables */
    if (ploop->rvar > 0) {
	ploop->ntimes += 1;
	if (ploop->ineq == GT) {
	    if (Z[ploop->lvar][0] > Z[ploop->rvar][0])
		return 1;
	    else return 0;
	} else {
	    if (Z[ploop->lvar][0] < Z[ploop->rvar][0])
		return 1;
	    else return 0;
	}
    } 
    else if (ploop->lvar == INDEXNUM) {  /* a 'for' loop */
	t = genr_scalar_index(0, 0); /* fetch index */
	if (t == ploop->ntimes) return 0;
	else {
	    genr_scalar_index(2, 1); /* increment index */
	    return 1;
	}
    }
    else if (ploop->lvar) {
    /* case of inequality between a var and a number */
	ploop->ntimes += 1;
	if (ploop->ineq == GT) {
	    if (Z[ploop->lvar][0] > ploop->rval) return 1;
	    else return 0;
	} else {
	    if (Z[ploop->lvar][0] < ploop->rval) return 1;
	    else return 0;
	}
    }
    else {
    /* case of a simple count */
	if (k < ploop->ntimes) return 1;
	else return 0;
    }
}

/* ......................................................  */

static int monte_carlo_init (LOOPSET *ploop)
{
    ploop->ntimes = 0;
    ploop->lvar = 0;
    ploop->rvar = 0;
    ploop->rval = 0;
    ploop->ineq = 0;
    ploop->ncmds = 0;
    ploop->nmod = 0;
    ploop->nprn = 0;
    ploop->nstore = 0;
    ploop->next_model = 0;
    ploop->next_print = 0;
    ploop->lines = malloc(32 * sizeof *ploop->lines); 
    ploop->ci = malloc(32 * sizeof *ploop->ci);
    ploop->models = NULL;
    ploop->lmodels = NULL;
    ploop->prns = NULL;
    return 0;
}

/**
 * monte_carlo_free:
 * @ploop: loop commands struct.
 *
 * Free allocated elements of @ploop.
 *
 */

void monte_carlo_free (LOOPSET *ploop)
{
    int i;

    if (ploop->lines) {
	for (i=0; i<ploop->ncmds; i++)
	    if (ploop->lines[i]) free (ploop->lines[i]);
	free(ploop->lines);
	ploop->lines = NULL;
    }
    if (ploop->ci) 
	free(ploop->ci);
    if (ploop->lmodels) {
	for (i=0; i<ploop->nmod; i++) 
	    free_loop_model(&ploop->lmodels[i]);
	free(ploop->lmodels);
	ploop->lmodels = NULL;
    }
    if (ploop->models) {
	for (i=0; i<ploop->nmod; i++) {
	    free_model(ploop->models[i]);
	    ploop->models[i] = NULL;
	}
	free(ploop->models);
	ploop->models = NULL;
    }    
    if (ploop->prns) {
	for (i=0; i<ploop->nprn; i++) 
	    free_loop_print(&ploop->prns[i]);
	free(ploop->prns);
	ploop->prns = NULL;
    }
    if (ploop->storename) {
	for (i=0; i<ploop->nstore; i++)
	    if (ploop->storename[i])
		free(ploop->storename[i]);
	free(ploop->storename);
	ploop->storename = NULL;
    }
    if (ploop->storelbl) {
	for (i=0; i<ploop->nstore; i++)
	    if (ploop->storelbl[i])
		free(ploop->storelbl[i]);
	free(ploop->storelbl);
	ploop->storelbl = NULL;
    }
    if (ploop->storeval) { 
	free(ploop->storeval);
	ploop->storeval = NULL;
    }
}

/**
 * loop_model_init:
 * @plmod: pointer to struct to initialize.
 * @pmod: model to take as basis.
 * @id: ID number to assign to @plmod.
 *
 * Initialize a #LOOP_MODEL struct, based on @pmod.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int loop_model_init (LOOP_MODEL *plmod, const MODEL *pmod,
		     int id)
{
    int i, ncoeff = pmod->ncoeff;

    plmod->sum_coeff = malloc((ncoeff + 1) * sizeof *plmod->sum_coeff);
    if (plmod->sum_coeff == NULL) return 1;

    plmod->ssq_coeff = malloc((ncoeff + 1) * sizeof *plmod->ssq_coeff);
    if (plmod->ssq_coeff == NULL) goto cleanup;

    plmod->sum_sderr = malloc((ncoeff + 1) * sizeof *plmod->sum_sderr);
    if (plmod->sum_sderr == NULL) goto cleanup;

    plmod->ssq_sderr = malloc((ncoeff + 1) * sizeof *plmod->ssq_sderr);
    if (plmod->ssq_sderr == NULL) goto cleanup;

    plmod->list = NULL;
    copylist(&plmod->list, pmod->list);
    if (plmod->list == NULL) goto cleanup;

    for (i=0; i<=ncoeff; i++) {
#ifdef ENABLE_GMP
	mpf_init(plmod->sum_coeff[i]);
	mpf_init(plmod->ssq_coeff[i]);
	mpf_init(plmod->sum_sderr[i]);
	mpf_init(plmod->ssq_sderr[i]);
#else
	plmod->sum_coeff[i] = plmod->ssq_coeff[i] = 0.0;
	plmod->sum_sderr[i] = plmod->ssq_sderr[i] = 0.0;
#endif
    }

    plmod->ncoeff = ncoeff;
    plmod->t1 = pmod->t1;
    plmod->t2 = pmod->t2;
    plmod->nobs = pmod->nobs;
    plmod->ID = id;
    plmod->ci = pmod->ci;
    return 0;

 cleanup:
    free(plmod->ssq_coeff);
    free(plmod->sum_sderr);
    free(plmod->ssq_sderr);
    return 1;
}

/**
 * loop_print_init:
 * @pprn: pointer to struct to initialize.
 * @list: list of variables to be printed.
 * @id: ID number to assign to @pprn.
 *
 * Initialize a #LOOP_PRINT struct.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int loop_print_init (LOOP_PRINT *pprn, const LIST list, int id)
{
    int i;

    pprn->list = NULL;

    copylist(&pprn->list, list);
    if (pprn->list == NULL) return 1;

    pprn->sum = malloc(list[0] * sizeof *pprn->sum);
    if (pprn->sum == NULL) goto cleanup;

    pprn->ssq = malloc(list[0] * sizeof *pprn->ssq);
    if (pprn->ssq == NULL) goto cleanup;

    for (i=0; i<list[0]; i++) { 
#ifdef ENABLE_GMP
	mpf_init(pprn->sum[i]);
	mpf_init(pprn->ssq[i]);
#else
	pprn->sum[i] = pprn->ssq[i] = 0.0;
#endif
    }

    pprn->ID = id;

    return 0;

 cleanup:
    free(pprn->list);
    free(pprn->sum);
    free(pprn->ssq);
    return 1;
}

/**
 * loop_store_init:
 * @ploop: pointer to loop struct.
 * @list: list of variables to be stored (written to file).
 * @pdinfo: data information struct.
 *
 * Set up @ploop for saving a set of variables.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int loop_store_init (LOOPSET *ploop, const LIST list, DATAINFO *pdinfo)
{
    int i, tot = list[0] * ploop->ntimes;

    ploop->storename = malloc(list[0] * sizeof *ploop->storename);
    if (ploop->storename == NULL) return 1;

    ploop->storelbl = malloc(list[0] * sizeof *ploop->storelbl);
    if (ploop->storelbl == NULL) goto cleanup;

    ploop->storeval = malloc(tot * sizeof *ploop->storeval);
    if (ploop->storeval == NULL) goto cleanup;

    for (i=0; i<list[0]; i++) {
	char *p;

	ploop->storename[i] = malloc(9);
	if (ploop->storename[i] == NULL) goto cleanup;

	strcpy(ploop->storename[i], pdinfo->varname[list[i+1]]);

	ploop->storelbl[i] = malloc(MAXLABEL);
	if (ploop->storelbl[i] == NULL) goto cleanup;

	strcpy(ploop->storelbl[i], pdinfo->label[list[i+1]]);
	if ((p = strstr(ploop->storelbl[i], "(scalar)"))) {
	    *p = 0;
	}
    }

    return 0;
 cleanup:
    free(ploop->storename);
    free(ploop->storelbl);
    free(ploop->storeval);
    return 1;
}

/**
 * update_loop_model:
 * @ploop: pointer to loop struct.
 * @cmdnum: sequential index number of the command within @ploop.
 * @pmod: contains estimates from the current iteration.
 *
 * Update a #LOOP_MODEL belonging to @ploop, based on the results
 * in @pmod.
 *
 * Returns: 0 on successful completion.
 */

int update_loop_model (LOOPSET *ploop, int cmdnum, MODEL *pmod)
{
    int j, i = get_modnum_by_cmdnum(ploop, cmdnum);
    LOOP_MODEL *plmod;
#ifdef ENABLE_GMP
    mpf_t m;

    mpf_init(m);
#endif

    plmod = &ploop->lmodels[i];

    for (j=1; j<=pmod->ncoeff; j++) {
#ifdef ENABLE_GMP
	mpf_set_d(m, pmod->coeff[j]);
	mpf_add(plmod->sum_coeff[j], plmod->sum_coeff[j], m);
	mpf_mul(m, m, m);
	mpf_add(plmod->ssq_coeff[j], plmod->ssq_coeff[j], m);

	mpf_set_d(m, pmod->sderr[j]);
	mpf_add(plmod->sum_sderr[j], plmod->sum_sderr[j], m);
	mpf_mul(m, m, m);
	mpf_add(plmod->ssq_sderr[j], plmod->ssq_sderr[j], m);
#else
	plmod->sum_coeff[j] += pmod->coeff[j];
	plmod->ssq_coeff[j] += pmod->coeff[j] * pmod->coeff[j];
	plmod->sum_sderr[j] += pmod->sderr[j];
	plmod->ssq_sderr[j] += pmod->sderr[j] * pmod->sderr[j];
#endif
    }

#ifdef ENABLE_GMP
    mpf_clear(m);
#endif

    return 0;
}

/**
 * update_loop_print:
 * @ploop: pointer to loop struct.
 * @cmdnum: sequential index number of the command within @ploop.
 * @list: list of variables to be printed.
 * @pZ: pointer to data matrix.
 * @pdinfo: pointer to data information struct.
 *
 * Update a #LOOP_PRINT belonging to @ploop, based on the current
 * data values.
 *
 * Returns: 0 on successful completion.
 */


int update_loop_print (LOOPSET *ploop, int cmdnum, 
		       const LIST list, double ***pZ, 
		       const DATAINFO *pdinfo)
{
    int j, t, i = get_prnnum_by_id(ploop, cmdnum);
    LOOP_PRINT *pprn = &ploop->prns[i];
#ifdef ENABLE_GMP
    mpf_t m;

    mpf_init(m);
#endif
    
    for (j=1; j<=list[0]; j++) {
	if (pdinfo->vector[list[j]]) t = pdinfo->t1;
	else t = 0;
#ifdef ENABLE_GMP
	mpf_set_d(m, (*pZ)[list[j]][t]); 
	mpf_add(pprn->sum[j-1], pprn->sum[j-1], m);
	mpf_mul(m, m, m);
	mpf_add(pprn->ssq[j-1], pprn->ssq[j-1], m);
#else
	pprn->sum[j-1] += (*pZ)[list[j]][t];
	pprn->ssq[j-1] += (*pZ)[list[j]][t] * (*pZ)[list[j]][t];
#endif
    }

#ifdef ENABLE_GMP
    mpf_clear(m);
#endif

    return 0;
}

/**
 * print_loop_results:
 * @ploop: pointer to loop struct.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 * @ppaths: path information struct.
 * @model_count: pointer to count of models estimated so far.
 * @loopstorefile: name of file into which to save data (or NULL).
 *
 * Print out the results after completion of the loop @ploop.
 *
 */

void print_loop_results (LOOPSET *ploop, const DATAINFO *pdinfo, 
			 PRN *prn, PATHS *ppaths, int *model_count,
			 char *loopstorefile)
{
    int i, j;
    MODEL *pmod = NULL;

    if (ploop->lvar && ploop->lvar != INDEXNUM) 
	pprintf(prn, _("\nNumber of iterations: %d\n\n"), ploop->ntimes);

    for (i=0; i<ploop->ncmds; i++) {
	/*  pprintf(prn, "loop command %d: %s\n\n", i+1, ploop->lines[i]); */
	if (ploop->lvar) {
	    if (ploop->ci[i] == OLS) {
		pmod = ploop->models[ploop->next_model];

		*model_count += 1;
		pmod->ID = *model_count;

		/* std. errors are asymptotic; degrees of freedom
		 correction is not wanted */
		if (pmod->vcv) free(pmod->vcv);
		pmod->vcv = NULL;

		pmod->sigma = sqrt((1.0 / pmod->nobs) * pmod->ess);
		makevcv(pmod);
		for (j=1; j<=pmod->ncoeff; j++) {
		    pmod->sderr[j] *= 
			sqrt((double) pmod->dfd /(double) pmod->nobs);
		}
		printmodel(pmod, pdinfo, prn);
		if (pmod->correct) {
		    /* -o flag was given */
		    outcovmx(pmod, pdinfo, 0, prn);
		}
		ploop->next_model += 1;	    
		continue;
	    }
	}
	if (ploop->ci[i] == OLS) {
	    print_loop_model(&ploop->lmodels[ploop->next_model], 
			     ploop->ntimes, pdinfo, prn);
	    ploop->next_model += 1;
	}
	else if (ploop->ci[i] == PRINT) {
	    print_loop_prn(&ploop->prns[ploop->next_print], 
			   ploop->ntimes, pdinfo, prn);
	    ploop->next_print += 1;
	}
	else if (ploop->ci[i] == STORE) {
	    print_loop_store(ploop, prn, ppaths, loopstorefile);
	}
    }
}

/* ......................................................  */

static void free_loop_model (LOOP_MODEL *plmod)
{
#ifdef ENABLE_GMP
    int i;

    for (i=0; i<=plmod->ncoeff; i++) {
	mpf_clear(plmod->sum_coeff[i]);
	mpf_clear(plmod->sum_sderr[i]);
	mpf_clear(plmod->ssq_coeff[i]);
	mpf_clear(plmod->ssq_sderr[i]);
    }
#endif
    free(plmod->sum_coeff);
    free(plmod->sum_sderr);
    free(plmod->ssq_coeff);
    free(plmod->ssq_sderr);
    free(plmod->list);
}

/* ......................................................  */

static void free_loop_print (LOOP_PRINT *pprn)
{
#ifdef ENABLE_GMP
    int i;

    for (i=0; i<pprn->list[0]; i++) {
	mpf_clear(pprn->sum[i]);
	mpf_clear(pprn->ssq[i]);
    }
#endif
    free(pprn->sum);
    free(pprn->ssq);
    free(pprn->list);    
}

/**
 * add_to_loop:
 * @ploop: pointer to loop struct.
 * @line: command line.
 * @ci: command index number.
 * @opt: option flag associated with the command.
 *
 * Add line and command index to accumulated loop buffer.
 *
 * Returns: 0 on successful completion.
 */

int add_to_loop (LOOPSET *ploop, char *line, int ci,
		 int opt)
{
    int i = ploop->ncmds;

    ploop->ncmds += 1;

    ploop->lines[i] = malloc(MAXLEN);
    if (ploop->lines[i] == NULL) return E_ALLOC;

    top_n_tail(line);

    if (ci == PRINT && ploop->type != COUNT_LOOP) {
	ploop->ci[i] = 0;
    } else {
	ploop->ci[i] = ci;
    }

    strncpy(ploop->lines[i], line, MAXLEN - 4);

    if (opt) {
	char flagstr[4];

	sprintf(flagstr, " -%c", getflag(opt));
	strcat(ploop->lines[i], flagstr);
    }

    return 0;
}

/* ......................................................... */ 

static const char *estimator_string (int ci)
{
    if (ci == OLS || ci == VAR) return N_("OLS");
    else if (ci == WLS) return N_("WLS"); 
    else if (ci == ARCH) return N_("WLS (ARCH)");
    else if (ci == TSLS) return N_("TSLS");
    else if (ci == HSK) return N_("Heteroskedasticity-corrected");
    else if (ci == AR) return N_("AR");
    else if (ci == LAD) return N_("LAD");
    else if (ci == HCCM) return N_("HCCM");
    else if (ci == PROBIT) return N_("Probit");
    else if (ci == LOGIT) return N_("Logit");
    else if (ci == POOLED) return N_("Pooled OLS");
    else if (ci == CORC) return N_("Cochrane-Orcutt");
    else if (ci == HILU) return N_("Hildreth-Lu");
    else return "";
}

/* ......................................................... */ 

static void print_loop_model (LOOP_MODEL *plmod, int loopnum,
			      const DATAINFO *pdinfo, PRN *prn)
{
    int i, nc = plmod->ncoeff;
    char startdate[9], enddate[9];
    int t1 = plmod->t1, t2 = plmod->t2;

    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(prn, _("%s estimates using the %d observations %s-%s\n"),
	    _(estimator_string(plmod->ci)), t2 - t1 + 1, 
	    startdate, enddate);
    pprintf(prn, _("Statistics for %d repetitions\n"), loopnum); 
    pprintf(prn, _("Dependent variable: %s\n\n"), 
	   pdinfo->varname[plmod->list[1]]);

    pputs(prn, _("                     mean of      std. dev. of     mean of"
	    "     std. dev. of\n"
	    "                    estimated      estimated"
	    "      estimated      estimated\n"
	    "      Variable     coefficients   coefficients   std. errors"
	    "    std. errors\n\n"));

    for (i=1; i<=nc; i++) {
	print_loop_coeff(pdinfo, plmod, i, loopnum, prn);
    }
    pputs(prn, "\n");
}

/* ......................................................... */ 

static void print_loop_coeff (const DATAINFO *pdinfo, 
			      const LOOP_MODEL *plmod, 
			      int c, int n, PRN *prn)
{
#ifdef ENABLE_GMP
    bigval c1, c2, m, sd1, sd2;

    mpf_init(c1);
    mpf_init(c2);
    mpf_init(m);
    mpf_init(sd1);
    mpf_init(sd2);

    mpf_div_ui(c1, plmod->sum_coeff[c], (unsigned long) n);
    mpf_mul(m, c1, c1);
    mpf_mul_ui(m, m, (unsigned long) n);
    mpf_sub(m, plmod->ssq_coeff[c], m);
    mpf_div_ui(sd1, m, (unsigned long) n);
    if (mpf_cmp_d(sd1, 0.0) > 0) {
	mpf_sqrt(sd1, sd1);
    } else {
	mpf_set_d(sd1, 0.0);
    }

    mpf_div_ui(c2, plmod->sum_sderr[c], (unsigned long) n);
    mpf_mul(m, c2, c2);
    mpf_mul_ui(m, m, (unsigned long) n);
    mpf_sub(m, plmod->ssq_sderr[c], m);
    mpf_div_ui(sd2, m, (unsigned long) n);
    if (mpf_cmp_d(sd2, 0.0) > 0) {
	mpf_sqrt(sd2, sd2);
    } else {
	mpf_set_d(sd2, 0.0);
    }

    pprintf(prn, " %3d) %8s ", plmod->list[c+1], 
	   pdinfo->varname[plmod->list[c+1]]);

    pprintf(prn, "%#14g %#14g %#14g %#14g\n", mpf_get_d(c1), mpf_get_d(sd1), 
	    mpf_get_d(c2), mpf_get_d(sd2));

    mpf_clear(c1);
    mpf_clear(c2);
    mpf_clear(m);
    mpf_clear(sd1);
    mpf_clear(sd2);
#else /* non-GMP */
    bigval m1, m2, var1, var2, sd1, sd2;
    
    m1 = plmod->sum_coeff[c] / n;
    var1 = (plmod->ssq_coeff[c] - n * m1 * m1) / n;
    sd1 = (var1 <= 0.0)? 0.0 : sqrt((double) var1);

    m2 = plmod->sum_sderr[c] / n;
    var2 = (plmod->ssq_sderr[c] - n * m2 * m2) / n;
    sd2 = (var2 <= 0.0)? 0 : sqrt((double) var2);

    pprintf(prn, " %3d) %8s ", plmod->list[c+1], 
	   pdinfo->varname[plmod->list[c+1]]);

    pprintf(prn, "%#14g %#14g %#14g %#14g\n", (double) m1, (double) sd1, 
	    (double) m2, (double) sd2);
#endif
}

/* ......................................................... */ 

static void print_loop_prn (LOOP_PRINT *pprn, int n,
			    const DATAINFO *pdinfo, PRN *prn)
{
    int i;
    bigval mean, m, sd;

    if (pprn == NULL) return;

    pputs(prn, _("   Variable     mean         std. dev.\n"));

#ifdef ENABLE_GMP
    mpf_init(mean);
    mpf_init(m);
    mpf_init(sd);
    
    for (i=1; i<=pprn->list[0]; i++) {
	mpf_div_ui(mean, pprn->sum[i-1], (unsigned long) n);
	mpf_mul(m, mean, mean);
	mpf_mul_ui(m, m, (unsigned long) n);
	mpf_sub(sd, pprn->ssq[i-1], m);
	mpf_div_ui(sd, sd, (unsigned long) n);
	if (mpf_cmp_d(sd, 0.0) > 0) {
	    mpf_sqrt(sd, sd);
	} else {
	    mpf_set_d(sd, 0.0);
	}
	pprintf(prn, " %8s ", pdinfo->varname[pprn->list[i]]);
	pprintf(prn, "%#14g %#14g\n", mpf_get_d(mean), mpf_get_d(sd));
    }

    mpf_clear(mean);
    mpf_clear(m);
    mpf_clear(sd);
#else
    for (i=1; i<=pprn->list[0]; i++) {
	mean = pprn->sum[i-1] / n;
	m = (pprn->ssq[i-1] - n * mean * mean) / n;
	sd = (m < 0)? 0 : sqrt((double) m);
	pprintf(prn, " %8s ", pdinfo->varname[pprn->list[i]]);
	pprintf(prn, "%#14g %#14g\n", (double) mean, (double) sd);
    }
#endif
    pputs(prn, "\n");
}

/* ......................................................... */ 

static int print_loop_store (LOOPSET *ploop, PRN *prn, PATHS *ppaths,
			     char *loopstorefile)
{
    int i, t, dot;
    FILE *fd;
    char prefix[MAXLEN], datfile[MAXLEN], hdrfile[MAXLEN], lblfile[MAXLEN];
    time_t writetime;

    /* organize filenames */
    sprintf(datfile, "%s%s", ppaths->userdir, loopstorefile);
    dot = dotpos(datfile);

    *prefix = 0;
    strncat(prefix, datfile, dot);
    sprintf(hdrfile, "%s.hdr", prefix);
    sprintf(lblfile, "%s.lbl", prefix);

    /* print to header file */
    writetime = time(NULL);
    fd = fopen(hdrfile, "w");
    if (fd == NULL) return 1;

    pprintf(prn, _("printing data header info to %s\n"), hdrfile);
    fprintf(fd, "(*\n %s %s\n*)\n", _("simulation data written"), 
	    print_time(&writetime));

    for (i=0; i<ploop->nstore; i++) {
	fprintf(fd, "%s ", ploop->storename[i]);
    }

    fprintf(fd, ";\n1 1 %d\nBYOBS\n", ploop->ntimes);
    fclose(fd);

    /* print to label file */
    fd = fopen(lblfile, "w");
    if (fd == NULL) return 1;

    for (i=0; i<ploop->nstore; i++) {
	fprintf(fd, "%s %s\n", ploop->storename[i], ploop->storelbl[i]);
    }
    fclose(fd);

    /* print to data file */
    fd = fopen(datfile, "w");    
    if (fd == NULL) return 1;

    pprintf(prn, _("printing %d values of variables to %s\n"), 
	    ploop->ntimes, datfile);

    for (t=0; t<ploop->ntimes; t++) {
	for (i=0; i<ploop->nstore; i++) {
	    fprintf(fd, "%g ", ploop->storeval[ploop->ntimes*i + t]);
	}
	fputc('\n', fd);
    }

    fclose(fd);
    return 0;
}

/* ......................................................... */ 

static int get_prnnum_by_id (LOOPSET *ploop, int id)
{
    int i;

    for (i=0; i<ploop->nprn; i++) {
	if (ploop->prns[i].ID == id) return i;
    }
    return -1;
}

/**
 * get_modnum_by_cmdnum:
 * @ploop: pointer to loop struct.
 * @cmdnum: sequential index of command within @ploop.
 *
 * Determine the ID number of a model within a "while" loop construct.
 *
 * Returns: model ID number, or -1 in case of no match.
 */

int get_modnum_by_cmdnum (LOOPSET *ploop, int cmdnum)
{
    int i;

    for (i=0; i<ploop->nmod; i++) {
	if (ploop->lvar && (ploop->models[i])->ID == cmdnum) return i;
	if (ploop->lvar == 0 && ploop->lmodels[i].ID == cmdnum) return i;
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
    if (sscanf(line, "%s", command->cmd) != 1 || line[0] == '(') {
	command->nolist = 1;
	command->ci = -1;
	return;
    }
    if ((command->ci = command_number(command->cmd)) == 0) {
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
    sprintf(formula, "iftest=%s", line + 2);
    err = generate(pZ, pdinfo, formula, 0, NULL, 1);
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


