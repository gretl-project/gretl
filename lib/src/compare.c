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

/*
    compare.c - gretl model comparison procedures
*/

#include "libgretl.h"
#include "gretl_private.h"
#include "gretl_matrix.h"

#ifdef WIN32
# include <windows.h>
#endif

/* ........................................................... */

void difflist (int *biglist, int *smalist, int *targ)
{
    int i, j, k, match;

    targ[0] = biglist[0] - smalist[0];
    k = 1;

    for (i=2; i<=biglist[0]; i++) {
	match = 0;
	for (j=2; j<=smalist[0]; j++) {
	    if (smalist[j] == biglist[i]) {
		match = 1;
		break;
	    }
	}
	if (!match) {
	    targ[k++] = biglist[i];
	}
    }
}

/* ........................................................... */

static COMPARE 
add_or_omit_compare (const MODEL *pmodA, const MODEL *pmodB, int add)
{
    COMPARE cmp;
    const MODEL *umod, *rmod;
    int i;	

    if (add) {
	umod = pmodB;
	rmod = pmodA;
    } else {
	umod = pmodA;
	rmod = pmodB;
    }

    cmp.m1 = pmodA->ID;
    cmp.m2 = pmodB->ID;
    cmp.ci = pmodA->ci;

    cmp.F = cmp.chisq = cmp.trsq = 0.0;
    cmp.score = 0;

    if (gretl_model_get_int(pmodA, "robust") || pmodA->ci == HCCM) {
	cmp.robust = 1;
    } else {
	cmp.robust = 0;
    }

    cmp.dfn = umod->ncoeff - rmod->ncoeff;
    cmp.dfd = umod->dfd;

    if (LIMDEP(cmp.ci)) {
	cmp.chisq = 2.0 * (umod->lnL - rmod->lnL);
	return cmp;
    }

    if (cmp.ci == OLS && !cmp.robust) {
	cmp.F = ((rmod->ess - umod->ess) / umod->ess) * cmp.dfd / cmp.dfn;
    }

    for (i=0; i<8; i++) { 
	if (pmodB->criterion[i] < pmodA->criterion[i]) {
	    cmp.score++;
	}
    }

    return cmp;
}

/* ........................................................... */

static int *
full_model_list (const MODEL *pmod, const int *inlist, int *ppos)
/* reconstitute full varlist for WLS and AR models */
{
    int i, len, pos = 0;
    int *flist = NULL;

    if (pmod->ci != WLS && pmod->ci != AR) 
	return NULL;

    if (pmod->ci == WLS) { 
	len = inlist[0] + 2;
    } else {
	pos = pmod->arinfo->arlist[0] + 1;
	len = pos + inlist[0] + 2;
    }

    flist = malloc(len * sizeof *flist);
    if (flist == NULL) return NULL;

    if (pmod->ci == WLS) { 
	flist[0] = len - 1;
	flist[1] = pmod->nwt;
	for (i=1; i<=inlist[0]; i++) {
	    flist[i+1] = inlist[i];
	}
    }
    else if (pmod->ci == AR) {
	flist[0] = len - 2;
	for (i=1; i<pos; i++) {
	    flist[i] = pmod->arinfo->arlist[i];
	}
	flist[pos] = LISTSEP;
	for (i=1; i<=inlist[0]; i++) {
	    flist[pos+i] = inlist[i];
	}
    }

    *ppos = pos;

    return flist;
}

/* ........................................................... */

static int be_quiet (gretlopt opt)
{
    if ((opt & OPT_A) || (opt & OPT_Q)) return 1;
    else return 0;
}

/* ........................................................... */

static MODEL replicate_estimator (const MODEL *orig, int **plist,
				  double ***pZ, DATAINFO *pdinfo,
				  gretlopt lsqopt, PRN *prn)
{
    MODEL rep;
    double rho = 0.0;
    int *list = *plist;
    int pos, mc = get_model_count();

    gretl_model_init(&rep);

    if (orig->ci == CORC || orig->ci == HILU || orig->ci == PWE) {
	rep.errcode = hilu_corc(&rho, list, pZ, pdinfo, NULL, 1, orig->ci, prn);
    }
    else if (orig->ci == WLS || orig->ci == AR) {
	int *full_list = full_model_list(orig, list, &pos);

	free(list);
	if (full_list == NULL) {
	    rep.errcode = E_ALLOC;
	} else {
	    list = *plist = full_list;
	}
    }

    if (rep.errcode) return rep;

    switch (orig->ci) {

    case AR:
	rep = ar_func(list, pos, pZ, pdinfo, lsqopt, prn);
	break;
    case ARCH:
	rep = arch(orig->order, list, pZ, pdinfo, NULL, lsqopt, prn);
	break;
    case LOGIT:
    case PROBIT:
	rep = logit_probit(list, pZ, pdinfo, orig->ci);
	break;
    case TOBIT:
	rep = tobit_model(list, pZ, pdinfo, NULL);
	break;
    case LAD:
	rep = lad(list, pZ, pdinfo);
	break;
    case LOGISTIC: 
	{
	    char lmaxstr[32];
	    double lmax;

	    lmax = gretl_model_get_double(orig, "lmax");
	    sprintf(lmaxstr, "lmax=%g", lmax);
	    rep = logistic_model(list, pZ, pdinfo, lmaxstr);
	}
	break;
    default:
	/* handles OLS, WLS, HSK, HCCM, etc. */
	if (gretl_model_get_int(orig, "robust")) lsqopt |= OPT_R;
	rep = lsq(list, pZ, pdinfo, orig->ci, lsqopt, rho);
	break;
    }

    /* check that we got the same sample as the original */
    if (!rep.errcode && rep.nobs < orig->nobs) {
	rep.errcode = E_MISS;
    }

    /* if the model count went up for an aux regression,
       bring it back down */
    if (be_quiet(lsqopt) && get_model_count() > mc) {
	model_count_minus();
    }

    return rep;
}

/**
 * auxreg:
 * @addvars: list of variables to add to original model (or NULL)
 * @orig: pointer to original model.
 * @new: pointer to new (modified) model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @aux_code: code indicating what sort of aux regression to run.
 * @test: hypothesis test results struct.
 * @opt: can contain options flags (--quiet, --vcv).
 * @prn: gretl printing struct.
 *
 * Run an auxiliary regression, in order to test a given set of added
 * variables, or to test for non-linearity (squares, logs).
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int auxreg (LIST addvars, MODEL *orig, MODEL *new, 
	    double ***pZ, DATAINFO *pdinfo, int aux_code, 
	    GRETLTEST *test, gretlopt opt, PRN *prn)
{
    COMPARE add;  
    MODEL aux;
    int *newlist = NULL, *tmplist = NULL;
    int i, j, t, listlen;
    const int n = pdinfo->n, orig_nvar = pdinfo->v; 
    double trsq = 0.0;
    int newvars = 0, err = 0;

    if (!command_ok_for_model(ADD, orig->ci)) {
	return E_NOTIMP;
    }

    if (aux_code != AUX_ADD && (orig->ci == LOGISTIC || orig->ci == LAD)) {
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(orig->list, pdinfo, orig->ID);
    if (err) return err;

    /* if adding specified vars, build the list */
    if (addvars != NULL) {
	newlist = gretl_list_add(orig->list, addvars, &err);
	if (err) return err;
    }

    /* temporarily re-impose the sample that was in force when the
       original model was estimated */
    exchange_smpl(orig, pdinfo);

    gretl_model_init(&aux);

    if (addvars == NULL) {
	listlen = orig->list[0] - orig->ifc;
	tmplist = gretl_list_new(listlen - 1);
	if (tmplist == NULL) {
	    err = E_ALLOC;
	} else {
	    j = 2;
	    for (i=1; i<=tmplist[0]; i++) {
		if (orig->list[j] == 0) j++;
		tmplist[i] = orig->list[j++];
	    }
	    /* no cross-products yet */
	    if (aux_code == AUX_SQ) { 
		/* add squares of original variables */
		newvars = xpxgenr(tmplist, pZ, pdinfo, 0, 0);
		if (newvars < 0) {
		    fprintf(stderr, "gretl: generation of squares failed\n");
		    free(tmplist);
		    err = E_SQUARES;
		}
	    }
	    else if (aux_code == AUX_LOG) { /* add logs of orig vars */
		newvars = logs(tmplist, pZ, pdinfo);
		if (newvars < 0) {
		    fprintf(stderr, "gretl: generation of logs failed\n");
		    free(tmplist);
		    err = E_LOGS;
		}
	    }	    
	    /* now construct an "addvars" list including all the
	       vars that were just generated -- re-use tmplist */
	    if (!err) {
		tmplist = realloc(tmplist, (newvars + 2) * sizeof *tmplist);
		if (tmplist == NULL) {
		    err = E_ALLOC;
		} else {
		    tmplist[0] = pdinfo->v - orig_nvar;
		    for (i=1; i<=tmplist[0]; i++) { 
			tmplist[i] = i + orig_nvar - 1;
		    }
		    newlist = gretl_list_add(orig->list, tmplist, &err);
		}
	    }
	} /* tmplist != NULL */
    }

    /* ADD: run an augmented regression, matching the original
       estimation method */
    if (!err && aux_code == AUX_ADD) {
	*new = replicate_estimator(orig, &newlist, pZ, pdinfo, opt, prn);
	if (new->errcode) {
	    err = new->errcode;
	    free(newlist);
	    if (addvars == NULL) {
		free(tmplist); 
	    }
	    clear_model(new);
	}
    } /* end if AUX_ADD */

    /* non-linearity test? Run auxiliary regression here -- 
       Replace depvar with uhat from orig */
    else if (!err && (aux_code == AUX_SQ || aux_code == AUX_LOG)) {
	int df = 0;

	/* grow data set to accommodate new dependent var */
	if (dataset_add_vars(1, pZ, pdinfo)) {
	    err = E_ALLOC;
	} else {
	    for (t=0; t<n; t++) {
		(*pZ)[pdinfo->v - 1][t] = orig->uhat[t];
	    }
	    newlist[1] = pdinfo->v - 1;

	    aux = lsq(newlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	    if (aux.errcode) {
		err = aux.errcode;
		fprintf(stderr, "auxiliary regression failed\n");
		free(newlist);
		if (addvars == NULL) free(tmplist); 
	    } else {
		aux.aux = aux_code;
		printmodel(&aux, pdinfo, opt, prn);
		trsq = aux.rsq * aux.nobs;

		if (test != NULL) {
		    gretl_test_init(test);
		    df = newlist[0] - orig->list[0];
		    strcpy(test->type, (aux_code == AUX_SQ)?
			    N_("Non-linearity test (squares)") :
			    N_("Non-linearity test (logs)"));
		    strcpy(test->h_0, N_("relationship is linear"));
		    test->teststat = GRETL_TEST_TR2;
		    test->dfn = df;
		    test->value = trsq;
		    test->pvalue = chisq(trsq, df);
		}
	    } /* ! aux.errcode */
	    clear_model(&aux);
	    /* shrink for uhat */
	    dataset_drop_vars(1, pZ, pdinfo);
	}
    }

    if (!err) {
	if (aux_code == AUX_ADD) {
	    new->aux = aux_code;
	}

	add = add_or_omit_compare(orig, new, 1);
	add.trsq = trsq;

	if (aux_code == AUX_ADD && !(opt & OPT_Q) && 
	    new->ci != AR && new->ci != ARCH) {
	    printmodel(new, pdinfo, opt, prn);
	}

	if (addvars != NULL) {
	    difflist(new->list, orig->list, addvars);
	    if (gretl_model_get_int(orig, "robust") || orig->ci == HCCM) {
		add.F = robust_omit_F(addvars, new);
	    }
	    gretl_print_add(&add, addvars, pdinfo, aux_code, prn, opt);
	} else {
	    add.dfn = newlist[0] - orig->list[0];
	    gretl_print_add(&add, tmplist, pdinfo, aux_code, prn, opt);
	}

	free(newlist);
	if (addvars == NULL) free(tmplist); 

	/* trash any extra variables generated (squares, logs) */
	if (pdinfo->v > orig_nvar) {
	    dataset_drop_vars(pdinfo->v - orig_nvar, pZ, pdinfo);
	}
    }

    /* put back into pdinfo what was there on input */
    exchange_smpl(orig, pdinfo);

    return err;
}

static int omit_index (int i, const int *list, const MODEL *pmod)
{
    /* return pos. in model coeff array of ith var to omit */
    int k = 0;

    if (list != NULL) {
	/* omitting a specific list of vars */
	int j, match = 0;

	for (j=2; j<=pmod->list[0]; j++) {
	    if (in_gretl_list(list, pmod->list[j])) {
		if (match == i) {
		    k = j - 2;
		    break;
		}
		match++;
	    }
	}
    } else {
	/* omitting all but the constant */
	k = i + pmod->ifc;
    }

    return k;
}

#undef FDEBUG

/**
 * robust_omit_F:
 * @list: list of variables to omit (or NULL).
 * @pmod: model to be tested.
 *
 * Simple form of Wald F-test for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated F-value, or #NADBL on failure.
 * 
 */

double robust_omit_F (const int *list, MODEL *pmod)
{
    int q, err = 0;
    gretl_matrix *sigma = NULL;
    gretl_vector *br = NULL;
    double F;
    int i, j, ii, jj, idx;
#ifdef FDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
#endif

    if (list != NULL) {
	q = list[0];
    } else {
	q = pmod->list[0] - 1 - pmod->ifc;
    } 
    
    sigma = gretl_matrix_alloc(q, q);
    br = gretl_column_vector_alloc(q);

    if (sigma == NULL || br == NULL) {
	gretl_matrix_free(sigma);
	gretl_matrix_free(br);
	return NADBL;
    }

    for (i=0; i<q; i++) {
	ii = omit_index(i, list, pmod);
	gretl_vector_set(br, i, pmod->coeff[ii]);
	for (j=0; j<=i; j++) {
	    jj = omit_index(j, list, pmod);
	    idx = ijton(ii, jj, pmod->ncoeff);
	    gretl_matrix_set(sigma, i, j, pmod->vcv[idx]);
	    if (i != j) {
		gretl_matrix_set(sigma, j, i, pmod->vcv[idx]);
	    }
	}
    }

#ifdef FDEBUG
    gretl_matrix_print(sigma, "sigma", prn);
#endif

    err = gretl_invert_symmetric_matrix(sigma);

#ifdef FDEBUG
    pprintf(prn, "invert returned %d\n", err);
#endif

    if (!err) {
	F = gretl_scalar_b_prime_X_b(br, sigma, &err);
    }

    if (err) {
	F = NADBL;
    } else {
	F /= q;
    }

    gretl_matrix_free(sigma);
    gretl_matrix_free(br);

#ifdef FDEBUG
    gretl_print_destroy(prn);
#endif

    return F;
}

/**
 * omit_test:
 * @omitvars: list of variables to omit from original model.
 * @orig: pointer to original model.
 * @new: pointer to new (modified) model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: can contain option flags (--quiet, --vcv).
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after removing a list of 
 * specified variables.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test (LIST omitvars, MODEL *orig, MODEL *new, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    COMPARE omit;
    int *tmplist;
    int maxlag = 0, t1 = pdinfo->t1;
    int err = 0;

    if (!command_ok_for_model(OMIT, orig->ci))
	return E_NOTIMP;

    /* check that vars to omit have not been redefined */
    if ((err = list_members_replaced(orig->list, pdinfo, orig->ID)))
	return err;

    /* create list for test model */
    tmplist = gretl_list_omit(orig->list, omitvars, &err);
    if (tmplist == NULL) {
	return err;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    exchange_smpl(orig, pdinfo);

    if (orig->ci == AR) { 
	maxlag = orig->arinfo->arlist[orig->arinfo->arlist[0]];
    } else if (orig->ci == ARCH) {
	maxlag = orig->order;
    }
    pdinfo->t1 = orig->t1 - maxlag; /* FIXME: problem */

    if (orig->ci == CORC || orig->ci == HILU) {
	pdinfo->t1 -= 1;
    }

    *new = replicate_estimator(orig, &tmplist, pZ, pdinfo, opt, prn);

    if (new->errcode) {
	pprintf(prn, "%s\n", gretl_errmsg);
	free(tmplist);
	err = new->errcode; 
    }

    if (!err) {
	if (orig->ci == LOGIT || orig->ci == PROBIT) {
	    new->aux = AUX_OMIT;
	}

	omit = add_or_omit_compare(orig, new, 0);

	if (!(opt & OPT_Q) && orig->ci != AR && orig->ci != ARCH) {
	    printmodel(new, pdinfo, opt, prn); 
	}

	difflist(orig->list, new->list, omitvars);

	if (gretl_model_get_int(orig, "robust") || orig->ci == HCCM) {
	    omit.F = robust_omit_F(omitvars, orig);
	}

	gretl_print_omit(&omit, omitvars, pdinfo, prn, opt); 

	free(tmplist);

	if (orig->ci == LOGIT || orig->ci == PROBIT) {
	    new->aux = AUX_NONE;
	}
    }

    pdinfo->t1 = t1;

    /* put back into pdinfo what was there on input */
    exchange_smpl(orig, pdinfo);

    return err;
}

static int ljung_box (int varno, int order, const double **Z, 
		      DATAINFO *pdinfo, double *lb)
{
    double *x, *y, *acf;
    int k, l, nobs, n = pdinfo->n; 
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];

    list[0] = 1;
    list[1] = varno;
    adjust_t1t2(NULL, list, &t1, &t2, Z, NULL);
    nobs = t2 - t1 + 1;

    x = malloc(n * sizeof *x);
    y = malloc(n * sizeof *y);
    acf = malloc((order + 1) * sizeof *acf);

    if (x == NULL || y == NULL || acf == NULL)
	return E_ALLOC;    

    for (l=1; l<=order; l++) {
	for (t=t1+l; t<=t2; t++) {
	    k = t - (t1+l);
	    x[k] = Z[varno][t];
	    y[k] = Z[varno][t-l];
	}
	acf[l] = gretl_corr(nobs-l, x, y);
    }

    /* compute Ljung-Box statistic */
    *lb = 0;
    for (t=1; t<=order; t++) { 
	*lb += acf[t] * acf[t] / (nobs - t);
    }
    *lb *= nobs * (nobs + 2.0);

    free(x);
    free(y);
    free(acf);

    return 0;
}

void gretl_test_init (GRETLTEST *test)
{
    test->type[0] = 0;
    test->h_0[0] = 0;
    test->param[0] = 0;
}

/**
 * reset_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @test: hypothesis test results struct.
 *
 * Ramsey's RESET test for model specification.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int reset_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		PRN *prn, GRETLTEST *test)
{
    int *newlist;
    MODEL aux;
    int i, t, v = pdinfo->v; 
    double RF;
    int err = 0;

    if (pmod->ci != OLS) return E_OLSONLY;

    gretl_model_init(&aux);

    if (pmod->ncoeff + 2 >= pdinfo->t2 - pdinfo->t1)
	return E_DF;

    newlist = malloc((pmod->list[0] + 3) * sizeof *newlist);
    if (newlist == NULL) {
	err = E_ALLOC;
    } else {
	newlist[0] = pmod->list[0] + 2;
	for (i=1; i<=pmod->list[0]; i++) {
	    newlist[i] = pmod->list[i];
	}
	if (dataset_add_vars(2, pZ, pdinfo)) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add yhat^2, yhat^3 to data set */
	for (t = pmod->t1; t<=pmod->t2; t++) {
	    double xx = pmod->yhat[t];

	    (*pZ)[v][t] = xx * xx;
	    (*pZ)[v+1][t] = xx * xx * xx;
	}
	strcpy(pdinfo->varname[v], "yhat^2");
	strcpy(pdinfo->varname[v+1], "yhat^3");
	newlist[pmod->list[0] + 1] = v;
	newlist[pmod->list[0] + 2] = v + 1;
    }

    if (!err) {
	aux = lsq(newlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    } 

    if (!err) {
	aux.aux = AUX_RESET;
	printmodel(&aux, pdinfo, OPT_NONE, prn);
	RF = ((pmod->ess - aux.ess) / 2) / (aux.ess / aux.dfd);

	pprintf(prn, "\n%s: F = %f,\n", _("Test statistic"), RF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		2, aux.dfd, RF, fdist(RF, 2, aux.dfd));

	if (test != NULL) {
	    gretl_test_init(test);
	    strcpy(test->type, N_("RESET test for specification"));
	    strcpy(test->h_0, N_("specification is adequate"));
	    test->teststat = GRETL_TEST_RESET;
	    test->dfn = 2;
	    test->dfd = aux.dfd;
	    test->value = RF;
	    test->pvalue = fdist(RF, test->dfn, test->dfd);
	}
    }

    free(newlist);
    dataset_drop_vars(2, pZ, pdinfo); 
    clear_model(&aux); 

    return err;
}

/* Below: apparatus for generating standard errors that are robust in 
   face of general serial correlation (see Wooldridge, Introductory
   Econometrics, chapter 12)
*/

static double get_vhat (double *ahat, int g, int t1, int t2)
{
    int t, h;
    double weight, a_cross_sum;
    double vhat;

    vhat = 0.0;

    for (t=t1; t<=t2; t++) {
	vhat += ahat[t] * ahat[t];
    }

    for (h=1; h<=g; h++) {
	weight = 1.0 - (double) h / (g + 1);
	a_cross_sum = 0.0;
	for (t=h+t1; t<=t2; t++) {
	    a_cross_sum += ahat[t] * ahat[t-h];
	}
	vhat += 2.0 * weight * a_cross_sum;
    }

    return vhat;
}

static int autocorr_standard_errors (MODEL *pmod, double ***pZ, 
				     DATAINFO *pdinfo, PRN *prn)
{
    int *auxlist = NULL;
    double *ahat = NULL;
    double *robust = NULL;
    double *tmp;
    int i, j, g;
    int aux = AUX_NONE, order = 0;
    MODEL auxmod;

    auxlist = malloc(pmod->list[0] * sizeof *auxlist);
    ahat = malloc(pdinfo->n * sizeof *ahat);
    robust = malloc(pmod->ncoeff * sizeof *robust);

    if (auxlist == NULL || ahat == NULL || robust == NULL) {
	free(auxlist);
	free(ahat);
	free(robust);
	return E_ALLOC;
    }

    g = get_hac_lag(pmod->nobs);

    auxlist[0] = pmod->list[0] - 1;

    gretl_model_init(&auxmod);

    /* loop across the indep vars in the original model */
    for (i=2; i<=pmod->list[0]; i++) {
	double vhat = 0;
	double sderr;
	int k, t;

	/* set the given indep var as the dependent */
	auxlist[1] = pmod->list[i];

	k = 2;
	for (j=2; j<=pmod->list[0]; j++) {
	    /* add other indep vars as regressors */
	    if (pmod->list[j] == auxlist[1]) continue;
	    auxlist[k++] = pmod->list[j];
	}

	auxmod = lsq(auxlist, pZ, pdinfo, OLS, OPT_A, 0.0);

	if (auxmod.errcode) {
	    fprintf(stderr, "Error estimating auxiliary model, code=%d\n", 
		    auxmod.errcode);
	    pmod->sderr[i-2] = NADBL;
	} else {
	    /* compute robust standard error */
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		ahat[t] = pmod->uhat[t] * auxmod.uhat[t];
	    }
	    vhat = get_vhat(ahat, g, pmod->t1, pmod->t2);
	    sderr = pmod->sderr[i-2] / pmod->sigma;
	    sderr = sderr * sderr;
	    sderr *= sqrt(vhat);
	    robust[i-2] = sderr;
	}

	clear_model(&auxmod);
    }

    /* save original model data */
    tmp = pmod->sderr;
    aux = pmod->aux;
    order = pmod->order;

    /* adjust data for SC-robust version */
    pmod->sderr = robust;
    pmod->aux = AUX_SCR;
    pmod->order = g;

    /* print original model, showing robust std errors */
    printmodel(pmod, pdinfo, OPT_NONE, prn);  

    /* reset the original model data */
    pmod->sderr = tmp;
    pmod->aux = aux;
    pmod->order = order;

    free(auxlist);
    free(ahat);
    free(robust);
	
    return 0;
}

/**
 * autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @test: hypothesis test results struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero. Gives TR^2 and LMF test statistics.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int autocorr_test (MODEL *pmod, int order, 
		   double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn, GRETLTEST *test)
{
    int *newlist;
    MODEL aux;
    int i, k, t, n = pdinfo->n, v = pdinfo->v; 
    double trsq, LMF, lb, pval = 1.0;
    int err = 0;

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == LOGISTIC) 
	return E_NOTIMP;

    if (dataset_is_panel(pdinfo)) {
	void *handle;
	int (*panel_autocorr_test)(MODEL *, int, 
				   double **, DATAINFO *, 
				   PRN *, GRETLTEST *);

	panel_autocorr_test = get_plugin_function("panel_autocorr_test", 
						  &handle);
	if (panel_autocorr_test == NULL) {
	    return 1;
	}

	err = panel_autocorr_test(pmod, order, *pZ, pdinfo,
				  prn, NULL);
	close_plugin(handle);
	return err;
    }

    exchange_smpl(pmod, pdinfo);
    gretl_model_init(&aux);

    if (order <= 0) order = pdinfo->pd;

    if (pmod->ncoeff + order >= pdinfo->t2 - pdinfo->t1)
	return E_DF;

    k = order + 1;
    newlist = malloc((pmod->list[0] + k) * sizeof *newlist);

    if (newlist == NULL) {
	err = E_ALLOC;
    } else {
	newlist[0] = pmod->list[0] + order;
	for (i=2; i<=pmod->list[0]; i++) newlist[i] = pmod->list[i];
	if (dataset_add_vars(1, pZ, pdinfo)) {
	    k = 0;
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add uhat to data set */
	for (t=0; t<n; t++)
	    (*pZ)[v][t] = NADBL;
	for (t = pmod->t1; t<= pmod->t2; t++)
	    (*pZ)[v][t] = pmod->uhat[t];
	strcpy(pdinfo->varname[v], "uhat");
	strcpy(VARLABEL(pdinfo, v), _("residual"));
	/* then lags of same */
	for (i=1; i<=order; i++) {
	    int lnum;

	    lnum = laggenr(v, i, 1, pZ, pdinfo);
	    if (lnum < 0) {
		sprintf(gretl_errmsg, _("lagging uhat failed"));
		err = E_LAGS;
	    } else {
		newlist[pmod->list[0] + i] = lnum;
	    }
	}
    }

    if (!err) {
	newlist[1] = v;
	/*  printlist(newlist); */
	aux = lsq(newlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    } 

    if (!err) {
	aux.aux = AUX_AR;
	aux.order = order;
	printmodel(&aux, pdinfo, OPT_NONE, prn);
	trsq = aux.rsq * aux.nobs;
	LMF = (aux.rsq/(1.0 - aux.rsq)) * 
	    (aux.nobs - pmod->ncoeff - order)/order; 

	pprintf(prn, "\n%s: LMF = %f,\n", _("Test statistic"), LMF);
	pval = fdist(LMF, order, aux.nobs - pmod->ncoeff - order);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, aux.nobs - pmod->ncoeff - order, LMF, pval);

	pprintf(prn, "\n%s: TR^2 = %f,\n", 
		_("Alternative statistic"), trsq);
	pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", 	_("with p-value"), 
		_("Chi-square"), order, trsq, chisq(trsq, order));

	/* add Ljung-Box Q' */
	if (ljung_box(v, order, (const double **) *pZ, pdinfo, &lb) == 0) {
	    pprintf(prn, "Ljung-Box Q' = %g %s = P(%s(%d) > %g) = %.3g\n", 
		    lb, _("with p-value"), _("Chi-square"), order,
		    lb, chisq(lb, order));
	}

	if (test != NULL) {
	    gretl_test_init(test);
	    strcpy(test->type, N_("LM test for autocorrelation up to order %s"));
	    strcpy(test->h_0, N_("no autocorrelation"));
	    sprintf(test->param, "%d", order);
	    test->teststat = GRETL_TEST_LMF;
	    test->dfn = order;
	    test->dfd = aux.nobs - pmod->ncoeff - order;
	    test->value = LMF;
	    test->pvalue = fdist(LMF, test->dfn, test->dfd);
	}
    }

    free(newlist);
    dataset_drop_vars(k, pZ, pdinfo); 
    clear_model(&aux); 

    if (pval < 0.05 && !gretl_model_get_int(pmod, "robust")) {
	autocorr_standard_errors(pmod, pZ, pdinfo, prn);
    }

    exchange_smpl(pmod, pdinfo);

    return err;
}

/**
 * chow_test:
 * @line: command line for parsing.
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @test: hypothesis test results struct.
 *
 * Tests the given model for structural stability (Chow test).
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test (const char *line, MODEL *pmod, double ***pZ,
	       DATAINFO *pdinfo, PRN *prn, GRETLTEST *test)
{
    int *chowlist = NULL;
    int newvars = pmod->list[0] - 1;
    int i, t, v = pdinfo->v, n = pdinfo->n;
    char chowdate[OBSLEN], s[OBSLEN];
    MODEL chow_mod;
    double F;
    int split = 0, err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (has_missing_obs(pmod)) {
	return E_DATA;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    exchange_smpl(pmod, pdinfo);

    gretl_model_init(&chow_mod);

    if (sscanf(line, "%*s %8s", chowdate) != 1) 
	err = E_PARSE;
    else {
	split = dateton(chowdate, pdinfo) - 1;
	if (split <= 0 || split >= pdinfo->n) 
	    err = E_SPLIT;
    }

    if (!err) {
	/* take the original regression list, add a split dummy
	   and interaction terms. */
	if (pmod->ifc == 0) newvars++;

	if (dataset_add_vars(newvars, pZ, pdinfo)) {
	    newvars = 0;
	    err = E_ALLOC;
	} else {
	    chowlist = malloc((pmod->list[0] + newvars + 1) * sizeof *chowlist);
	    if (chowlist == NULL) 
		err = E_ALLOC;
	}
    }

    if (!err) {
	chowlist[0] = pmod->list[0] + newvars;
	for (i=1; i<=pmod->list[0]; i++) { 
	    chowlist[i] = pmod->list[i];
	}

	/* generate the split variable */
	for (t=0; t<n; t++) 
	    (*pZ)[v][t] = (double) (t > split); 
	strcpy(pdinfo->varname[v], "splitdum");
	strcpy(VARLABEL(pdinfo, v), _("dummy variable for Chow test"));
	chowlist[pmod->list[0] + 1] = v;

	/* and the interaction terms */
	for (i=1; i<newvars; i++) {
	    int orig = i + 1 + pmod->ifc;

	    for (t=0; t<n; t++) {
		(*pZ)[v+i][t] = (*pZ)[v][t] * 
		    (*pZ)[pmod->list[orig]][t];
	    }
	    strcpy(s, pdinfo->varname[pmod->list[orig]]); 
	    gretl_trunc(s, 5);
	    strcpy(pdinfo->varname[v+i], "sd_");
	    strcat(pdinfo->varname[v+i], s);
	    sprintf(VARLABEL(pdinfo, v+i), "splitdum * %s", 
		    pdinfo->varname[pmod->list[orig]]);
	    chowlist[pmod->list[0]+1+i] = v+i;
	}

	chow_mod = lsq(chowlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (chow_mod.errcode) {
	    err = chow_mod.errcode;
	    errmsg(err, prn);
	} else {
	    chow_mod.aux = AUX_CHOW;
	    printmodel(&chow_mod, pdinfo, OPT_NONE, prn);
	    F = (pmod->ess - chow_mod.ess) * chow_mod.dfd / 
		(chow_mod.ess * newvars);
	    pprintf(prn, _("\nChow test for structural break at observation %s:\n"
		    "  F(%d, %d) = %f with p-value %f\n\n"), chowdate,
		    newvars, chow_mod.dfd, F, 
		    fdist(F, newvars, chow_mod.dfd)); 

	    if (test != NULL) {
		gretl_test_init(test);
		strcpy(test->type, N_("Chow test for structural break at "
			"observation %s"));
		strcpy(test->param, chowdate);
		strcpy(test->h_0, N_("no structural break"));
		test->teststat = GRETL_TEST_F;
		test->dfn = newvars;
		test->dfd = chow_mod.dfd;
		test->value = F;
		test->pvalue = fdist(F, newvars, chow_mod.dfd);
	    }
	}
	clear_model(&chow_mod);
    }

    /* clean up extra variables */
    dataset_drop_vars(newvars, pZ, pdinfo);
    free(chowlist);

    exchange_smpl(pmod, pdinfo);    

    return err;
}

/* ........................................................... */

static double vprime_M_v (double *v, double *M, int n)
     /* compute v'Mv, for symmetric M */
{
    int i, j, jmin, jmax, k;
    double xx, val = 0.0;

    k = jmin = 0;
    for (i=0; i<n; i++) {
	xx = 0.0;
	for (j=jmin; j<n; j++) {
	    xx += v[j] * M[k++];
	}
	val += xx * v[i];
	jmin++;
    }

    jmax = 1;
    for (i=1; i<n; i++) {
	k = i;
	xx = 0.0;
	for (j=0; j<jmax; j++) {
	    xx += v[j] * M[k];
	    k += n - j - 1;
	}
	val += xx * v[i];
	jmax++;
    }

    return val;
}

/**
 * cusum_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @ppaths: path information struct.
 * @test: hypothesis test results struct.
 *
 * Tests the given model for parameter stability (CUSUM test).
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int cusum_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, PRN *prn, 
		PATHS *ppaths, GRETLTEST *test)
{
    int n_est, i, j, t;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int xno, yno = pmod->list[1];
    const int T = pmod->t2 - pmod->t1 + 1;
    const int K = pmod->ncoeff;
    MODEL cum_mod;
    char cumdate[OBSLEN];
    double xx, yy, sigma, hct, wbar = 0.0;
    double *cresid = NULL, *W = NULL, *xvec = NULL;
    FILE *fq = NULL;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (has_missing_obs(pmod)) {
	return E_DATA;
    }

    n_est = T - K;
    /* set sample based on model to be tested */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t1 + K - 1;    

    cresid = malloc(n_est * sizeof *cresid);
    W = malloc(n_est * sizeof *W);
    xvec = malloc(K * sizeof *xvec);

    if (cresid == NULL || W == NULL || xvec == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	for (j=0; j<n_est && !err; j++) {
	    cum_mod = lsq(pmod->list, pZ, pdinfo, OLS, OPT_C, 0.0);
	    err = cum_mod.errcode;
	    if (err) {
		errmsg(err, prn);
	    } else {
		t = pdinfo->t2 + 1;
		yy = 0.0;
		for (i=0; i<K; i++) {
		    xno = cum_mod.list[i+2];
		    xvec[i] = (*pZ)[xno][t];
		    yy += cum_mod.coeff[i] * (*pZ)[xno][t];
		}
		cresid[j] = (*pZ)[yno][t] - yy;
		cum_mod.ci = CUSUM;
		makevcv(&cum_mod);
		xx = vprime_M_v(xvec, cum_mod.vcv, K);
		cresid[j] /= sqrt(1.0 + xx);
		/*  printf("w[%d] = %g\n", t, cresid[j]); */
		wbar += cresid[j];
		clear_model(&cum_mod);
		pdinfo->t2 += 1;
	    }
	    clear_model(&cum_mod); 
	}
    }

    if (!err) {
	wbar /= T - K;
	pprintf(prn, "\n%s\n\n",
		_("CUSUM test for stability of parameters"));
	pprintf(prn, _("mean of scaled residuals = %g\n"), wbar);
	sigma = 0;
	for (j=0; j<n_est; j++) {
	    xx = (cresid[j] - wbar);
	    sigma += xx * xx;
	}
	sigma /= T - K - 1;
	sigma = sqrt(sigma);
	pprintf(prn, _("sigmahat                 = %g\n\n"), sigma);

	xx = 0.948 * sqrt((double) (T-K));
	yy = 2.0 * xx / (T-K);

	pputs(prn, _("Cumulated sum of scaled residuals\n"
		"('*' indicates a value outside of 95%% confidence band):\n\n"));
    
	for (j=0; j<n_est; j++) {
	    W[j] = 0.0;
	    for (i=0; i<=j; i++) W[j] += cresid[i];
	    W[j] /= sigma;
	    t = pmod->t1 + K + j;
	    ntodate(cumdate, t, pdinfo);
	    /* FIXME printing of number below? */
	    pprintf(prn, " %s %9.3f %s\n", cumdate, W[j],
		    (fabs(W[j]) > xx + (j+1)*yy)? "*" : "");
	}
	hct = (sqrt((double) (T-K)) * wbar) / sigma;
	pprintf(prn, _("\nHarvey-Collier t(%d) = %g with p-value %.4g\n\n"), 
		T-K-1, hct, tprob(hct, T-K-1));

	if (test != NULL) {
	    gretl_test_init(test);
	    strcpy(test->type, N_("CUSUM test for parameter stability"));
	    strcpy(test->h_0, N_("no change in parameters"));
	    test->teststat = GRETL_TEST_HARVEY_COLLIER;
	    test->dfn = T-K-1;
	    test->value = hct;
	    test->pvalue = tprob(hct, T-K-1);
	}

#ifdef ENABLE_NLS
        setlocale(LC_NUMERIC, "C");
#endif
	/* plot with 95% confidence bands, if not batch mode */
	if (prn->fp == NULL && gnuplot_init(ppaths, PLOT_CUSUM, &fq) == 0) {
	    fputs("# CUSUM test\n", fq);
	    fprintf(fq, "set xlabel \"%s\"\n", I_("Observation"));
	    fputs("set xzeroaxis\n", fq);
	    fprintf(fq, "set title \"%s\"\n",
		    /* xgettext:no-c-format */
		    I_("CUSUM plot with 95% confidence band"));
	    fputs("set nokey\n", fq);
	    fprintf(fq, "plot %f+%f*x w l 1, \\\n", xx - K*yy, yy);
	    fprintf(fq, "%f-%f*x w l 1, \\\n", -xx + K*yy, yy);
	    fputs("'-' using 1:2 w lp\n", fq);
	    for (j=0; j<n_est; j++) { 
		t = pmod->t1 + K + j;
		fprintf(fq, "%d %f\n", t, W[j]);
	    }
	    fputs("e\n", fq);

	    fclose(fq);

	    err = gnuplot_display(ppaths);
	}
#ifdef ENABLE_NLS
        setlocale(LC_NUMERIC, "");
#endif
    }

    /* restore sample */
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;
    
    free(cresid);
    free(W);
    free(xvec);

    return err;
}

/**
 * hausman_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Tests the given pooled model for fixed and random effects.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int hausman_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		  PRN *prn) 
{
    if (pmod->ci != POOLED) {
	pputs(prn, _("This test is only relevant for pooled models\n"));
	return 1;
    }

    if (pmod->ifc == 0) {
	pputs(prn, _("This test requires that the model contains a constant\n"));
	return 1;
    }

    if (!balanced_panel(pdinfo)) {
	pputs(prn, _("Sorry, can't do this test on an unbalanced panel.\n"
		"You need to have the same number of observations\n"
		"for each cross-sectional unit"));
	return 1;
    } else {
	void *handle;
	void (*panel_diagnostics)(MODEL *, double ***, DATAINFO *, PRN *);

	panel_diagnostics = get_plugin_function("panel_diagnostics", &handle);
	if (panel_diagnostics == NULL) {
	    return 1;
	}
	(*panel_diagnostics) (pmod, pZ, pdinfo, prn);
	close_plugin(handle);
    }
    return 0;
}

/**
 * add_leverage_values_to_dataset:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @m: matrix containing leverage values.
 * @opt: option flag: combination of SAVE_LEVERAGE, SAVE_INFLUENCE,
 * and SAVE_DFFITS.
 *
 * Adds to the working dataset one or more series calculated by
 * the gretl test for leverage/influence of data points.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int add_leverage_values_to_dataset (double ***pZ, DATAINFO *pdinfo,
				    gretl_matrix *m, unsigned char opt)
{
    int t1, t2;
    int addvars = 0;

    if (opt & SAVE_LEVERAGE) addvars++;
    if (opt & SAVE_INFLUENCE) addvars++;
    if (opt & SAVE_DFFITS) addvars++;

    if (dataset_add_vars(addvars, pZ, pdinfo)) {
	strcpy(gretl_errmsg, _("Out of memory adding series"));
	return 1;
    }

    t1 = gretl_matrix_get_int(m);
    t2 = t1 + gretl_matrix_rows(m);

    /* add leverage? */
    if (opt & SAVE_LEVERAGE) {
	int t, v = pdinfo->v - addvars;
	int j = 0;

	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t >= t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = gretl_matrix_get(m, j++, 0);
	    }
	}
	strcpy(pdinfo->varname[v], "lever");
	make_varname_unique(pdinfo->varname[v], v, pdinfo);
	strcpy(VARLABEL(pdinfo, v), "leverage values");
    }

    /* add influence? */
    if (opt & SAVE_INFLUENCE) {
	int t, v = pdinfo->v - (addvars - 1);
	int j = 0;

	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t >= t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = gretl_matrix_get(m, j++, 1);
	    }
	}	
	strcpy(pdinfo->varname[v], "influ");
	make_varname_unique(pdinfo->varname[v], v, pdinfo);
	strcpy(VARLABEL(pdinfo, v), "influence values");
    }

    /* add DFFITS? */
    if (opt & SAVE_DFFITS) {
	int t, v = pdinfo->v - (addvars - 2);
	int j = 0;

	for (t=0; t<pdinfo->n; t++) {
	    double s, h;

	    if (t < t1 || t >= t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		/* s = studentized residuals */
		h = gretl_matrix_get(m, j, 0);
		s = gretl_matrix_get(m, j, 2);
		if (na(h) || na(s)) {
		    (*pZ)[v][t] = NADBL;
		} else {
		    (*pZ)[v][t] = s * sqrt(h / (1.0 - h));
		}
		j++;
	    }
	}	
	strcpy(pdinfo->varname[v], "dffits");
	make_varname_unique(pdinfo->varname[v], v, pdinfo);
	strcpy(VARLABEL(pdinfo, v), "DFFITS values");
    }

    return 0;
}

/**
 * leverage_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @ppaths: path information struct (should be NULL if a graph
 * is not wanted).
 * @oflag: if non-zero, add calculated series to data set.
 *
 * Tests the data used in the given model for points with
 * high leverage and influence on the estimates
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int leverage_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn, PATHS *ppaths, gretlopt oflag)
{
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, double ***, 
				     const DATAINFO *, PRN *, PATHS *);
    gretl_matrix *m;
    int err = 0;

    if (pmod->ci != OLS) return E_OLSONLY;

    model_leverage = get_plugin_function("model_leverage", &handle);
    if (model_leverage == NULL) {
	return 1;
    }

    m = (*model_leverage)(pmod, pZ, pdinfo, prn, ppaths);
    if (m == NULL) {
	err = 1;
    } else {
	if (oflag) {
	    err = add_leverage_values_to_dataset(pZ, pdinfo, m, 
						 SAVE_LEVERAGE |
						 SAVE_INFLUENCE| 
						 SAVE_DFFITS);
	}
	gretl_matrix_free(m);
    }

    close_plugin(handle);

    return err;
}

/**
 * vif_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Calculates and displays the Variance Inflation Factors for
 * the independent variables in the given model.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int vif_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    void *handle;
    int (*print_vifs) (MODEL *, double ***, DATAINFO *, PRN *);
    int err;

    *gretl_errmsg = '\0';

    print_vifs = get_plugin_function("print_vifs", &handle);
    if (print_vifs == NULL) {
	return 1;
    }

    err = (*print_vifs)(pmod, pZ, pdinfo, prn);

    close_plugin(handle);

    if (err && *gretl_errmsg == '\0') {
	gretl_errmsg_set(_("Command failed"));
    }

    return err;
}

int make_mp_lists (const LIST list, const char *str,
		   int **reglist, int **polylist)
{
    int i, pos;

    pos = atoi(str);

    *reglist = malloc(pos * sizeof **polylist);
    *polylist = malloc((list[0] - pos + 2) * sizeof **reglist);

    if (*reglist == NULL || *polylist == NULL) {
	free(*reglist);
	free(*polylist);
	return 1;
    }
    
    (*reglist)[0] = pos - 1;
    for (i=1; i<pos; i++) {
	(*reglist)[i] = list[i];
    }

    (*polylist)[0] = list[0] - pos;
    for (i=1; i<=(*polylist)[0]; i++) {
	(*polylist)[i] = list[i+pos];
    }

    return 0;
}

/**
 * mp_ols:
 * @list: specification of variables to use
 * @pos: string rep. of integer position in list at which
 * the regular list of variables ends and a list of polynomial
 * terms begins (or empty string in case of no polynomial terms)
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * 
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int mp_ols (const LIST list, const char *pos,
	    double ***pZ, DATAINFO *pdinfo, 
	    PRN *prn) 
{
    void *handle = NULL;
    int (*mplsq)(const int *, const int *, double ***, 
		 DATAINFO *, PRN *, char *, mp_results *);
    const int *reglist = NULL;
    int *polylist = NULL, *tmplist = NULL;
    mp_results *mpvals = NULL;
    int nc, err = 0;

    mplsq = get_plugin_function("mplsq", &handle);
    if (mplsq == NULL) return 1;

    if (!err && *pos) { /* got a list of polynomial terms? */
	err = make_mp_lists(list, pos, &tmplist, &polylist);
	if (err) {
	    pputs(prn, _("Failed to parse mp_ols command\n"));
	}
	reglist = tmplist;
    } 

    if (!err && !*pos) {
	reglist = list;
    }

    nc = list[0] - 1;
    if (polylist != NULL) nc--;

    mpvals = gretl_mp_results_new(nc);
    if (mpvals == NULL || allocate_mp_varnames(mpvals)) {
	pprintf(prn, "%s\n", _("Out of memory!"));
	err = 1;
    }

    if (!err) {
	err = (*mplsq)(reglist, polylist, pZ, pdinfo, prn, 
		       gretl_errmsg, mpvals); 
    }

    if (!err) {
	print_mpols_results(mpvals, pdinfo, prn);
    }

    close_plugin(handle);

    free(polylist);
    free(tmplist);
    free_gretl_mp_results(mpvals);

    return err;
}

static int varmatch (int *sumvars, int test)
{
    int j;

    for (j=1; j<=sumvars[0]; j++) {
	if (sumvars[j] == test) return 1;
    }

    return 0;
}

static 
void fill_sum_var (double **Z, int n, int v, int vrepl, int vfirst)
{
    int t;

    for (t=0; t<n; t++) {
	Z[v][t] = Z[vrepl][t] - Z[vfirst][t];
    }
}

static 
int make_sum_test_list (MODEL *pmod, double **Z, DATAINFO *pdinfo,
			int *tmplist, int *sumvars, int vstart)
{
    int repl = 0;
    int newv = vstart;
    int testcoeff = 0;
    int nnew = sumvars[0] - 1;
    int i;

    tmplist[0] = pmod->list[0];
    tmplist[1] = pmod->list[1];

    for (i=2; i<=pmod->list[0]; i++) {
	if (nnew > 0 && varmatch(sumvars, pmod->list[i])) {
	    if (repl) {
		fill_sum_var(Z, pdinfo->n, newv, pmod->list[i], sumvars[1]);
		tmplist[i] = newv++;
		nnew--;
	    } else {
		tmplist[i] = pmod->list[i];
		testcoeff = i;
		repl = 1;
	    }
	} else {
	    tmplist[i] = pmod->list[i];
	}
    }

    if (nnew == 0) return testcoeff;
    else return -1;
}

/**
 * sum_test:
 * @sumvars: specification of variables to use.
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * 
 * Calculates the sum of the coefficients, relative to the given model, 
 * for the variables given in @sumvars.  Prints this estimate along 
 * with its standard error.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int sum_test (LIST sumvars, MODEL *pmod, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn)
{
    int *tmplist;
    int add = sumvars[0] - 1;
    int v = pdinfo->v;
    int testcoeff;
    MODEL summod;
    PRN *nullprn;
    int err = 0;

    if (sumvars[0] < 2) {
	pprintf(prn, _("Invalid input\n"));
	return E_DATA;
    }

    if (!command_ok_for_model(COEFFSUM, pmod->ci)) 
	return E_NOTIMP;

    tmplist = malloc((pmod->list[0] + 1) * sizeof *tmplist);
    if (tmplist == NULL) return E_ALLOC;

    if (dataset_add_vars(add, pZ, pdinfo)) {
	free(tmplist);
	return E_ALLOC;
    }

    nullprn = gretl_print_new(GRETL_PRINT_NULL, NULL);

    testcoeff = make_sum_test_list(pmod, *pZ, pdinfo, tmplist, sumvars, v);

    if (testcoeff < 0) {
	pprintf(prn, _("Invalid input\n"));
	free(tmplist);
	dataset_drop_vars(add, pZ, pdinfo);
	return E_DATA;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    exchange_smpl(pmod, pdinfo);

    gretl_model_init(&summod);

    if (!err) {
	summod = replicate_estimator(pmod, &tmplist, pZ, pdinfo, OPT_A, 
				     nullprn);
	if (summod.errcode) {
	    pprintf(prn, "%s\n", gretl_errmsg);
	    err = summod.errcode; 
	} else {
	    int i;

	    pprintf(prn, "\n%s: ", _("Variables"));
	    for (i=1; i<=sumvars[0]; i++) {
		pprintf(prn, "%s ", pdinfo->varname[sumvars[i]]);
	    }
	    /* FIXME: check indexing of summod.coeff[] below */
	    pprintf(prn, "\n   %s = %g\n", _("Sum of coefficients"), 
		    summod.coeff[testcoeff - 2]);
	    if (!na(summod.sderr[testcoeff - 2])) {
		double tval;

		pprintf(prn, "   %s = %g\n", _("Standard error"),
			summod.sderr[testcoeff - 2]);
		tval = summod.coeff[testcoeff - 2] / 
		    summod.sderr[testcoeff - 2];
		pprintf(prn, "   t(%d) = %g ", summod.dfd, tval);
		pprintf(prn, _("with p-value = %g\n"), 
			tprob(tval, summod.dfd));
	    }
	}
    }

    free(tmplist);
    clear_model(&summod);
    dataset_drop_vars(add, pZ, pdinfo);
    gretl_print_destroy(nullprn);

    /* put back into pdinfo what was there on input */
    exchange_smpl(pmod, pdinfo);

    return err;
}
