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

extern int _command_number (const char *cmd);
extern int genr_scalar_index (int opt, int put);

static int monte_carlo_init (LOOPSET *ploop);
static void _free_loop_model (LOOP_MODEL *plmod);
static void _free_loop_print (LOOP_PRINT *pprn);
static void _print_loop_model (LOOP_MODEL *plmod, const int loopnum,
			       const DATAINFO *pdinfo, print_t *prn);
static void _print_loop_coeff (const DATAINFO *pdinfo, const LOOP_MODEL *plmod, 
			       const int c, const int n, print_t *prn);
static void _print_loop_prn (LOOP_PRINT *pprn, const int n,
			     const DATAINFO *pdinfo, print_t *prn);
static int _print_loop_store (LOOPSET *ploop, print_t *prn, PATHS *ppaths,
			      char *loopstorefile);
static int get_prnnum_by_id (LOOPSET *ploop, const int id);

/* ......................................................  */

int ok_in_loop (int ci)
{
    if (ci == OLS || 
	ci == GENR ||
	ci == STORE ||
	ci == PRINT ||
	ci == SMPL ||
	ci == SUMMARY ||
	ci == ENDLOOP) 
	return 1;
    return 0;
}

/* ......................................................  */

int parse_loopline (char *line, LOOPSET *ploop, DATAINFO *pdinfo)
{
    char lvar[9], rvar[9], op[8];
    int start, end;
    int n, v;

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
		    "Undefined variable '%s' in loop condition.", lvar);
	    return 1;
	}
	if (isdigit((unsigned char) rvar[0]) 
	    || rvar[0] == '.') { /* numeric rvalue? */
	    ploop->rval = atof(rvar);
	    return 0; /* OK, done */
	} /* otherwise try a varname */
	v = varindex(pdinfo, rvar);
	if (v > 0 && v < pdinfo->v) {
	    ploop->rvar = v;
	    return 0;
	} else {
	    sprintf(gretl_errmsg, 
		    "Undefined variable '%s' in loop condition.", rvar);
	    ploop->lvar = 0;
	    return 1;
	}
    }

    /* or try parsing as a for loop */
    else if (sscanf(line, "loop for %[^= ] = %d..%d", lvar, &start, &end) == 3) {
	if (strcmp(lvar, "i")) {
	    sprintf(gretl_errmsg, 
		    "The index variable in a 'for' loop must be the "
		    "special variable 'i'");
	    return 1;
	}
	if (end <= start) {
	    sprintf(gretl_errmsg, "Ending value for loop index must be greater "
		    "than starting value.");
	    return 1;
	}
	/* initialize special genr index to starting value */
	genr_scalar_index(1, start - 1);
	ploop->lvar = INDEXNUM;
	ploop->rvar = 0;
	ploop->ntimes = end;
	return 0;
    }

    /* or just as a simple count loop */
    else if (sscanf(line, "loop %d", &n) == 1) {
	ploop->ntimes = n;
	return 0;
    }

    /* out of options, complain */
    strcpy(gretl_errmsg, "No valid loop condition was given.");
    return 1;
}

/* ......................................................  */

int loop_condition (int k, LOOPSET *ploop, double *Z, DATAINFO *pdinfo)
{
    int n = pdinfo->n, t = pdinfo->t2;
    const int LOOPMAX = 2500; /* safety measure */

    if (ploop->lvar && ploop->ntimes > LOOPMAX) 
	return 0;

    /* case of an inequality between variables */
    if (ploop->rvar > 0) {
	ploop->ntimes += 1;
	if (ploop->ineq == GT) {
	    if (Z(ploop->lvar, t) > Z(ploop->rvar, t))
		return 1;
	    else return 0;
	} else {
	    if (Z(ploop->lvar, t) < Z(ploop->rvar, t))
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
	    if (Z(ploop->lvar, t) > ploop->rval) return 1;
	    else return 0;
	} else {
	    if (Z(ploop->lvar, t) < ploop->rval) return 1;
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
    ploop->lines = malloc(32 * sizeof(char *)); 
    ploop->ci = malloc(32 * sizeof(int));
    ploop->models = NULL;
    ploop->lmodels = NULL;
    ploop->prns = NULL;
    return 0;
}

/* ......................................................  */

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
	    _free_loop_model(&ploop->lmodels[i]);
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
	    _free_loop_print(&ploop->prns[i]);
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

/* ......................................................  */

int loop_model_init (LOOP_MODEL *plmod, const MODEL *pmod,
		     const int id)
{
    int i, ncoeff = pmod->ncoeff;

    plmod->sum_coeff = malloc((ncoeff + 1) * sizeof(long double));
    if (plmod->sum_coeff == NULL) return 1;
    plmod->ssq_coeff = malloc((ncoeff + 1) * sizeof(long double));
    if (plmod->ssq_coeff == NULL) goto cleanup;
    plmod->sum_sderr = malloc((ncoeff + 1) * sizeof(long double));
    if (plmod->sum_sderr == NULL) goto cleanup;
    plmod->ssq_sderr = malloc((ncoeff + 1) * sizeof(long double));
    if (plmod->ssq_sderr == NULL) goto cleanup;
    for (i=0; i<=ncoeff; i++) {
	plmod->sum_coeff[i] = plmod->ssq_coeff[i] = 0;
	plmod->sum_sderr[i] = plmod->ssq_sderr[i] = 0;
    }
    plmod->ncoeff = ncoeff;
    plmod->t1 = pmod->t1;
    plmod->t2 = pmod->t2;
    plmod->nobs = pmod->nobs;
    plmod->list = NULL;
    copylist(&plmod->list, pmod->list);
    plmod->ID = id;
    return 0;

 cleanup:
    free(plmod->ssq_coeff);
    free(plmod->sum_sderr);
    free(plmod->ssq_sderr);
    return 1;
}

/* ......................................................  */

int loop_print_init (LOOP_PRINT *pprn, const int *list, const int id)
{
    int i;

    pprn->list = NULL;
    copylist(&pprn->list, list);
    if (pprn->list == NULL) return 1;
    pprn->sum = malloc(list[0] * sizeof(long double));
    if (pprn->sum == NULL) goto cleanup;
    pprn->ssq = malloc(list[0] * sizeof(long double));
    if (pprn->ssq == NULL) goto cleanup;
    for (i=0; i<list[0]; i++) 
	pprn->sum[i] = pprn->ssq[i] = 0;
    pprn->ID = id;
    return 0;
 cleanup:
    free(pprn->list);
    free(pprn->sum);
    free(pprn->ssq);
    return 1;
}

/* ......................................................  */

int loop_store_init (LOOPSET *ploop, const int *list, DATAINFO *pdinfo)
{
    int i, tot;

    ploop->storename = malloc(list[0] * sizeof(char *));
    if (ploop->storename == NULL) return 1;
    ploop->storelbl = malloc(list[0] * sizeof(char *));
    if (ploop->storelbl == NULL) goto cleanup;
    tot = list[0] * ploop->ntimes;
    ploop->storeval = malloc(tot * sizeof(double));
    if (ploop->storeval == NULL) goto cleanup;
    for (i=0; i<list[0]; i++) {
	ploop->storename[i] = malloc(9);
	if (ploop->storename[i] == NULL) goto cleanup;
	strcpy(ploop->storename[i], pdinfo->varname[list[i+1]]);
	ploop->storelbl[i] = malloc(MAXLABEL);
	if (ploop->storelbl[i] == NULL) goto cleanup;
	strcpy(ploop->storelbl[i], pdinfo->label[list[i+1]]);
    }
    return 0;
 cleanup:
    free(ploop->storename);
    free(ploop->storelbl);
    free(ploop->storeval);
    return 1;
}

/* ......................................................  */

int update_loop_model (LOOPSET *ploop, const int cmdnum,
		       MODEL *pmod, const DATAINFO *pdinfo, 
		       const int oflag)
{
    int j, i = get_modnum_by_id(ploop, cmdnum);
    LOOP_MODEL *plmod;

    plmod = &ploop->lmodels[i];
    for (j=1; j<=pmod->ncoeff; j++) {
	plmod->sum_coeff[j] += pmod->coeff[j];
	plmod->ssq_coeff[j] += pmod->coeff[j] * pmod->coeff[j];
	plmod->sum_sderr[j] += pmod->sderr[j];
	plmod->ssq_sderr[j] += pmod->sderr[j] * pmod->sderr[j];
    }
    return 0;
}

/* ......................................................  */

int update_loop_print (LOOPSET *ploop, const int cmdnum, 
		       const int *list, double **pZ, 
		       const int n, const int t)
{
    int j, i = get_prnnum_by_id(ploop, cmdnum);
    LOOP_PRINT *pprn = &ploop->prns[i];
    
    for (j=1; j<=list[0]; j++) {
	pprn->sum[j-1] += (*pZ)[list[j]*n + t];
	pprn->ssq[j-1] += (*pZ)[list[j]*n + t] * (*pZ)[list[j]*n + t];
    }
    return 0;
}

/* ......................................................  */

void print_loop_results (LOOPSET *ploop, const DATAINFO *pdinfo, 
			 print_t *prn, PATHS *ppaths, int *model_count,
			 char *loopstorefile)
{
    int i, j;
    MODEL *pmod = NULL;

    if (ploop->lvar && ploop->lvar != INDEXNUM) 
	pprintf(prn, "\nNumber of iterations: %d\n\n", ploop->ntimes);

    for (i=0; i<ploop->ncmds; i++) {
	/*  pprintf(prn, "loop command %d: %s\n\n", i+1, ploop->lines[i]); */
	if (ploop->lvar) {
	    if (ploop->ci[i] == OLS) {
		pmod = ploop->models[ploop->next_model];
		*model_count += 1;
		pmod->ID = *model_count;
		/* std. errors are asymptotic; degrees of freedom
		 correction is not wanted */
		free(pmod->vcv);
		pmod->vcv = NULL;
		pmod->sigma = sqrt((1.0/pmod->nobs) * pmod->ess);
		makevcv(pmod);
		for (j=1; j<=pmod->ncoeff; j++)
		    pmod->sderr[j] *= 
			sqrt((double) pmod->dfd /(double) pmod->nobs);
		printmodel(pmod, pdinfo, prn);
		if (pmod->correct) /* -o flag was given */
		    outcovmx(pmod, pdinfo, 0, prn);
		ploop->next_model += 1;	    
		continue;
	    }
	}
	if (ploop->ci[i] == OLS) {
	    _print_loop_model(&ploop->lmodels[ploop->next_model], 
			     ploop->ntimes, pdinfo, prn);
	    ploop->next_model += 1;
	}
	else if (ploop->ci[i] == PRINT) {
	    _print_loop_prn(&ploop->prns[ploop->next_print], 
			   ploop->ntimes, pdinfo, prn);
	    ploop->next_print += 1;
	}
	else if (ploop->ci[i] == STORE) {
	    _print_loop_store(ploop, prn, ppaths, loopstorefile);
	}
    }
}

/* ......................................................  */

static void _free_loop_model (LOOP_MODEL *plmod)
{
    free(plmod->sum_coeff);
    free(plmod->sum_sderr);
    free(plmod->ssq_coeff);
    free(plmod->ssq_sderr);
    free(plmod->list);
}

/* ......................................................  */

static void _free_loop_print (LOOP_PRINT *pprn)
{
    free(pprn->list);
    free(pprn->sum);
    free(pprn->ssq);
}

/* ......................................................  */

int add_to_loop (LOOPSET *ploop, char *line, const int ci,
		 double **pZ, DATAINFO *pdinfo, const int opt)
/* for loop command, add line and command index to accumulated buffer */
{
    int i = ploop->ncmds;

    ploop->ncmds += 1;
    ploop->lines[i] = malloc(MAXLEN);
    if (ploop->lines[i] == NULL) return E_ALLOC;
    top_n_tail(line);
    ploop->ci[i] = ci;
    strncpy(ploop->lines[i], line, MAXLEN - 4);
    if (opt) {
	char flagstr[4];

	sprintf(flagstr, " -%c", getflag(opt));
	strcat(ploop->lines[i], flagstr);
    }
    return 0;
}

/* ......................................................... */ 

static void _print_loop_model (LOOP_MODEL *plmod, const int loopnum,
			       const DATAINFO *pdinfo, print_t *prn)
{
    int i, nc = plmod->ncoeff;
    char startdate[8];
    char enddate[8];
    int t1 = plmod->t1, t2 = plmod->t2;

    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(prn, "OLS estimates using the %d observations %s-%s\n",
	   t2-t1+1, startdate, enddate);
    pprintf(prn, "Statistics for %d repetitions\n", loopnum); 
    pprintf(prn, "Dependent variable: %s\n\n", 
	   pdinfo->varname[plmod->list[1]]);

    pprintf(prn, "                     mean of      std. dev. of     mean of"
	    "     std. dev. of\n"
	    "                    estimated      estimated"
	    "      estimated      estimated\n"
	    "      Variable     coefficients   coefficients   std. errors"
	    "    std. errors\n\n");
    for (i=1;i<=nc; i++) 
	_print_loop_coeff(pdinfo, plmod, i, loopnum, prn);
    pprintf(prn, "\n");
}

/* ......................................................... */ 

static void _print_loop_coeff (const DATAINFO *pdinfo, 
			       const LOOP_MODEL *plmod, 
			       const int c, const int n, print_t *prn)
{
    long double m1, m2, var1, var2, sd1, sd2;

    m1 = plmod->sum_coeff[c] / n;
    m1 = plmod->sum_coeff[c] / n;
    var1 = (plmod->ssq_coeff[c] - n * m1 * m1) / n;
    sd1 = (var1 < 0)? 0 : sqrt((double) var1);
    m2 = plmod->sum_sderr[c] / n;
    var2 = (plmod->ssq_sderr[c] - n * m2 * m2) / n;
    sd2 = (var2 < 0)? 0 : sqrt((double) var2);

    pprintf(prn, " %3d) %8s ", plmod->list[c+1], 
	   pdinfo->varname[plmod->list[c+1]]);
    pprintf(prn, "%14f %14f %14f %14f\n", (double) m1, (double) sd1, 
	    (double) m2, (double) sd2);
}

/* ......................................................... */ 

static void _print_loop_prn (LOOP_PRINT *pprn, const int n,
			     const DATAINFO *pdinfo, print_t *prn)
{
    int i;
    long double mean, var, sd;

    if (pprn == NULL) return;

    pprintf(prn, "   Variable     mean         std. dev.\n");
    for (i=1; i<=pprn->list[0]; i++) {
	mean = pprn->sum[i-1] / n;
	var = (pprn->ssq[i-1] - n * mean * mean) / n;
	sd = (var < 0)? 0 : sqrt((double) var);
	pprintf(prn, " %8s ", pdinfo->varname[pprn->list[i]]);
	pprintf(prn, "%14f %14f\n", (double) mean, (double) sd);
    }
    pprintf(prn, "\n");
}

/* ......................................................... */ 

static int _print_loop_store (LOOPSET *ploop, print_t *prn, PATHS *ppaths,
			      char *loopstorefile)
{
    int i, t, dot;
    FILE *fd;
    char prefix[MAXLEN], datfile[MAXLEN], hdrfile[MAXLEN], lblfile[MAXLEN];
    time_t writetime;

    /* organize filenames */
    sprintf(datfile, "%s%s", ppaths->userdir, loopstorefile);
    dot = dotpos(datfile);
    clear(prefix, MAXLEN);
    strncpy(prefix, datfile, dot);
    clear(hdrfile, MAXLEN); 
    sprintf(hdrfile, "%s.hdr", prefix);
    clear(lblfile, MAXLEN);
    sprintf(lblfile, "%s.lbl", prefix);

    /* print to header file */
    writetime = time(NULL);
    fd = fopen(hdrfile, "w");
    if (fd == NULL) return 1;
    pprintf(prn, "printing data header info to %s\n", hdrfile);
    fprintf(fd, "(*\n simulation data written %s\n*)\n", 
	    ctime(&writetime));
    for (i=0; i<ploop->nstore; i++)
	fprintf(fd, "%s ", ploop->storename[i]);
    fprintf(fd, ";\n1 1 %d\nBYOBS\n", ploop->ntimes);
    fclose(fd);

    /* print to label file */
    fd = fopen(lblfile, "w");
    if (fd == NULL) return 1;
    for (i=0; i<ploop->nstore; i++)
	fprintf(fd, "%s %s\n", ploop->storename[i], ploop->storelbl[i]);
    fclose(fd);

    /* print to data file */
    fd = fopen(datfile, "w");    
    if (fd == NULL) return 1;
    pprintf(prn, "printing %d values of variables to %s\n", 
	    ploop->ntimes, datfile);
    for (t=0; t<ploop->ntimes; t++) {
	for (i=0; i<ploop->nstore; i++) 
	    fprintf(fd, "%f ", ploop->storeval[ploop->ntimes*i + t]);
	fprintf(fd, "\n");
    }
    fclose(fd);
    return 0;
}

/* ......................................................... */ 

static int get_prnnum_by_id (LOOPSET *ploop, const int id)
{
    int i;

    for (i=0; i<ploop->nprn; i++) 
	if (ploop->prns[i].ID == id) return i;
    return -1;
}

/* ......................................................... */ 

int get_modnum_by_id (LOOPSET *ploop, const int id)
{
    int i;

    for (i=0; i<ploop->nmod; i++) {
	if (ploop->lvar && (ploop->models[i])->ID == id) return i;
	if (ploop->lvar == 0 && ploop->lmodels[i].ID == id) return i;
    }
    return -1;
}

/* ......................................................... */

void get_cmd_ci (const char *line, CMD *command)
{
    if (sscanf(line, "%s", command->cmd) != 1 || line[0] == '(') {
	command->nolist = 1;
	command->ci = -1;
	return;
    }
    if ((command->ci = _command_number(command->cmd)) == 0) {
	command->errcode = 1;
	sprintf(gretl_errmsg, "command \"%s\" not recognized", 
		command->cmd);
	return;
    }    
} 
