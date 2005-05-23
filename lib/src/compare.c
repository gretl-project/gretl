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
#include "libset.h"
#include "gretl_matrix.h"
#include "missing_private.h"

#undef WDEBUG

struct COMPARE {
    int cmd;       /* ADD or OMIT */
    int m1;        /* ID for first model */
    int m2;        /* ID for second model */
    int ci;        /* estimator code for the first model */
    int dfn;       /* numerator degrees of freedom */
    int dfd;       /* denominator degrees of freedom */ 
    double F;      /* F test statistic */
    double chisq;  /* Chi-square test statistic */
    double trsq;   /* T*R^2 test statistic */
    int score;     /* "cases correct" for discrete models */
    int robust;    /* = 1 when robust vcv is in use, else 0 */
};

/* Given a list of variables, check them against the independent
   variables included in a model, and construct a mask with 1s in
   positions where there is a match, 0s otherwise.  If the test list
   is NULL, match all variables except the constant.
*/

static char *
mask_from_test_list (const int *test, const MODEL *pmod)
{
    char *mask;
    int i, j;

    mask = calloc(pmod->ncoeff, 1);
    if (mask == NULL) {
	return NULL;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	if (test != NULL) {
	    for (j=1; j<=test[0]; j++) {
		if (pmod->list[i+2] == test[j]) {
		    mask[i] = 1;
		}
	    }
	} else {
	    if (pmod->list[i+2] != 0) {
		mask[i] = 1;
	    }
	}
    }

    return mask;
}

enum {
    CHI_SQUARE_FORM,
    F_FORM
};

/* Wald (chi-square or F) test for a set of zero restrictions on the
   parameters of a given model, based on the covariance matrix of the
   unrestricted model. Suitable for use where the original model is
   estimated by FGLS or IV.  Note that if list is NULL, we do an
   automatic test, for the significance of all vars but the constant.
*/

static double wald_test (const int *list, MODEL *pmod, int form)
{
    char *mask = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    double w = NADBL;
    int err = 0;

    mask = mask_from_test_list(list, pmod);
    if (mask == NULL) {
	err = 1;
    } 

    if (!err) {
	C = gretl_vcv_matrix_from_model(pmod, mask);
	if (C == NULL) {
	    err = 1;
	} 
    }

    if (!err) {
	b = gretl_coeff_vector_from_model(pmod, mask);
	if (b == NULL) {
	    err = 1;
	} 
    }  

    if (!err) {
#if WDEBUG
	gretl_matrix_print(C, "Wald VCV matrix", NULL);
	gretl_matrix_print(b, "Wald coeff vector", NULL);
#endif
	err = gretl_invert_symmetric_matrix(C);
    }

    if (!err) {
	w = gretl_scalar_b_prime_X_b(b, C, &err);
    }

    if (form == F_FORM && !err) {
	w /= gretl_vector_get_length(b);
    }

#if WDEBUG
    fprintf(stderr, "Wald test: %s = %g\n", 
	    (form == F_FORM)? "F" : "Chi^2", w);
#endif

    free(mask);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return w;
}

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
    return wald_test(list, pmod, F_FORM);
}

/* ----------------------------------------------------- */

static int 
add_diffvars_to_test (ModelTest *test, const int *list, 
		      const DATAINFO *pdinfo)
{
    char *vnames;
    int i, len = 0;
    int err = 0;

    for (i=1; i<=list[0]; i++) {
	len += strlen(pdinfo->varname[list[i]]) + 1;
    }

    vnames = malloc(len);

    if (vnames == NULL) {
	err = 1;
    } else {
	*vnames = '\0';
	for (i=1; i<=list[0]; i++) {
	    strcat(vnames, pdinfo->varname[list[i]]);
	    if (i < list[0]) {
		strcat(vnames, " ");
	    }
	}
	model_test_set_allocated_param(test, vnames);
    }

    return err;
}

static void
gretl_print_compare (const struct COMPARE *cmp, const int *diffvars, 
		     MODEL *new, const DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn)
{
    ModelTest *test = NULL;
    double pval = NADBL;
    int stat_ok, i;

    if (cmp->ci == LAD) {
	return;
    }

    stat_ok = !na(cmp->F) || !na(cmp->chisq);

    if (stat_ok && (opt & OPT_S)) {
	test = new_test_on_model(new, (cmp->cmd == OMIT)? 
				 GRETL_TEST_OMIT : GRETL_TEST_ADD);
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, _("Comparison of Model %d and Model %d:\n"), 
		cmp->m1, cmp->m2);
    } 

    if (stat_ok) {
	pputs(prn, _("\n  Null hypothesis: the regression parameters are "
		     "zero for the variables\n\n"));
	for (i=1; i<=diffvars[0]; i++) {
	    pprintf(prn, "    %s\n", pdinfo->varname[diffvars[i]]);	
	}
    }

    if (!na(cmp->F)) {
	pprintf(prn, "\n  %s: %s(%d, %d) = %g, ", _("Test statistic"), 
		(cmp->robust)? _("Robust F") : "F",
		cmp->dfn, cmp->dfd, cmp->F);
	pval = fdist(cmp->F, cmp->dfn, cmp->dfd);
	pprintf(prn, _("with p-value = %g\n"), pval);
	record_test_result(cmp->F, pval, (cmp->cmd == OMIT)? "omit" : "add");

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_F);
	    model_test_set_dfn(test, cmp->dfn);
	    model_test_set_dfd(test, cmp->dfd);
	    model_test_set_value(test, cmp->F);
	    model_test_set_pvalue(test, pval);
	}
    } else if (!na(cmp->chisq)) {
	pprintf(prn, "\n  %s:%s%s(%d) = %g, ",  
		(LIMDEP(cmp->ci))? _("Test statistic") : 
		_("Asymptotic test statistic"),
		(LIMDEP(cmp->ci))? " " : "\n    ",
		_("Chi-square"), cmp->dfn, cmp->chisq);
	pval = chisq(cmp->chisq, cmp->dfn);
	pprintf(prn, _("with p-value = %g\n\n"), pval);
	record_test_result(cmp->chisq, pval, (cmp->cmd == OMIT)? "omit" : "add");

	if (test != NULL) {
	    model_test_set_teststat(test, (LIMDEP(cmp->ci))? 
				    GRETL_STAT_LR : GRETL_STAT_WALD_CHISQ);
	    model_test_set_dfn(test, cmp->dfn);
	    model_test_set_value(test, cmp->chisq);
	    model_test_set_pvalue(test, pval);
	}
    } 	

    if (test != NULL) {
	add_diffvars_to_test(test, diffvars, pdinfo);
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, _("  Of the %d model selection statistics, %d "), 
		C_MAX, cmp->score);
	if (cmp->score == 1) {
	    pputs(prn, _("has improved.\n"));
	} else {
	    pputs(prn, _("have improved.\n\n"));
	}
    }
}

static struct COMPARE 
add_or_omit_compare (MODEL *pmodA, MODEL *pmodB, int add,
		     const int *testvars)
{
    struct COMPARE cmp;
    MODEL *umod, *rmod;
    int i;	

    if (add) {
	umod = pmodB;
	rmod = pmodA;
	cmp.cmd = ADD;
    } else {
	umod = pmodA;
	rmod = pmodB;
	cmp.cmd = OMIT;
    }

    cmp.m1 = pmodA->ID;
    cmp.m2 = pmodB->ID;
    cmp.ci = pmodA->ci;

    cmp.F = cmp.chisq = cmp.trsq = NADBL;
    cmp.score = 0;
    cmp.robust = 0;
    
    cmp.dfn = umod->ncoeff - rmod->ncoeff;
    cmp.dfd = umod->dfd;

    /* FIXME TSLS (F or chi-square?) */

    if (gretl_model_get_int(pmodA, "robust") || pmodA->ci == HCCM) {
	cmp.F = robust_omit_F(testvars, umod);
	cmp.robust = 1;
    } else if (LIMDEP(cmp.ci)) {
	cmp.chisq = 2.0 * (umod->lnL - rmod->lnL);
    } else if (cmp.ci == OLS) {
	cmp.F = ((rmod->ess - umod->ess) / umod->ess) * cmp.dfd / cmp.dfn;
    } else if (cmp.dfn > 1) {
	cmp.chisq = wald_test(testvars, umod, CHI_SQUARE_FORM);
    }

    for (i=0; i<C_MAX; i++) { 
	if (na(pmodB->criterion[i]) || na(pmodA->criterion[i])) {
	    continue;
	}
	if (pmodB->criterion[i] < pmodA->criterion[i]) {
	    cmp.score++;
	}
    }

    return cmp;
}

static int get_tsls_pos (const int *list)
{
    int i, pos = 0;
    
    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    pos = i;
	    break;
	}
    }

    return pos;
}

/* reconstitute full varlist for WLS, POISSON, AR and TSLS models */

static int *
full_model_list (const MODEL *pmod, const int *inlist, int *ppos)
{
    int i, len, pos = 0;
    int *flist = NULL;

    if (pmod->ci != WLS && pmod->ci != POISSON && 
	pmod->ci != AR && pmod->ci != TSLS) {
	return NULL;
    }

    if (pmod->ci == WLS) { 
	len = inlist[0] + 2;
    } else if (pmod->ci == POISSON) {
	len = inlist[0] + 3;
    } else if (pmod->ci == TSLS) {
	pos = get_tsls_pos(pmod->list);
	len = inlist[0] + 2 + pmod->list[0] - pos; 
    } else {
	pos = pmod->arinfo->arlist[0] + 1;
	len = pos + inlist[0] + 2;
    }

    flist = malloc(len * sizeof *flist);
    if (flist == NULL) {
	return NULL;
    }

    if (pmod->ci == WLS) { 
	flist[0] = len - 1;
	flist[1] = pmod->nwt;
	for (i=1; i<=inlist[0]; i++) {
	    flist[i+1] = inlist[i];
	}
    } else if (pmod->ci == POISSON) {
	int offvar = gretl_model_get_int(pmod, "offset_var");

	flist[0] = len - 1;
	for (i=1; i<=inlist[0]; i++) {
	    flist[i] = inlist[i];
	}
	flist[flist[0] - 1] = LISTSEP;
	flist[flist[0]] = offvar;
    } else if (pmod->ci == TSLS) {
	int j = 1;

	flist[0] = len - 1;
	for (i=1; i<=inlist[0]; i++) {
	    flist[j++] = inlist[i];
	}
	for (i=pos; i<=pmod->list[0]; i++) {
	    flist[j++] = pmod->list[i];
	}
	pos = get_tsls_pos(flist);
    } else if (pmod->ci == AR) {
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

#define be_quiet(o) ((o & OPT_A) || (o & OPT_Q))

#undef SMPL_DEBUG

static MODEL replicate_estimator (MODEL *orig, int **plist,
				  double ***pZ, DATAINFO *pdinfo,
				  gretlopt lsqopt, PRN *prn)
{
    MODEL rep;
    double rho = 0.0;
    int *list = *plist;
    int pos, mc = get_model_count();
    int repci = orig->ci;

    gretl_model_init(&rep);

    if (orig->ci == CORC || orig->ci == HILU || orig->ci == PWE) {
	gretlopt hlopt = OPT_NONE;

	if (orig->ci == HILU && gretl_model_get_int(orig, "no-corc")) {
	    hlopt = OPT_B;
	}
	rho = estimate_rho(list, pZ, pdinfo, orig->ci, 
			   &rep.errcode, hlopt, prn);
    } else if (gretl_model_get_int(orig, "unit_weights")) {
	/* panel model with per-unit weights */
	lsqopt |= OPT_W;
	repci = POOLED;
    } else if (orig->ci == WLS || orig->ci == AR || orig->ci == TSLS ||
	       (orig->ci == POISSON && gretl_model_get_int(orig, "offset_var"))) {
	int *full_list = full_model_list(orig, list, &pos);

	free(list);
	if (full_list == NULL) {
	    rep.errcode = E_ALLOC;
	} else {
	    list = *plist = full_list;
	}
    }

    if (rep.errcode) {
	return rep;
    }

    switch (orig->ci) {

    case AR:
	rep = ar_func(list, pos, pZ, pdinfo, lsqopt, prn);
	break;
    case ARCH:
	rep = arch_model(list, orig->order, pZ, pdinfo, lsqopt, prn);
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
    case POISSON:
	rep = poisson_model(list, pZ, pdinfo, NULL);
	break;
    case TSLS:
	rep = tsls_func(list, pos, pZ, pdinfo, lsqopt);
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
	if (gretl_model_get_int(orig, "robust")) {
	    lsqopt |= OPT_R;
	}
	if (gretl_model_get_int(orig, "hc_version") == 4) {
	    repci = HCCM;
	}
	rep = lsq(list, pZ, pdinfo, repci, lsqopt, rho);
	break;
    }

#if SMPL_DEBUG
    fprintf(stderr, "replicate_estimator:\n"
	    " orig: t1=%d, t2=%d, nobs = %d\n"
	    " rep:  t1=%d, t2=%d, nobs = %d\n",
	    orig->t1, orig->t2, orig->nobs,
	    rep.t1, rep.t2, rep.nobs);
#endif

    /* check that we got the same sample as the original */
    if (!rep.errcode) {
	if (rep.nobs < orig->nobs) {
	    rep.errcode = E_MISS;
	} else if (rep.nobs > orig->nobs) {
	    rep.errcode = E_DATA;
	}
    } 

    /* if the model count went up for an aux regression,
       bring it back down */
    if (be_quiet(lsqopt) && get_model_count() > mc) {
	model_count_minus();
    }

    return rep;
}

static int
real_nonlinearity_test (MODEL *pmod, int *list,
			double ***pZ, DATAINFO *pdinfo,
			int aux_code, gretlopt opt, 
			PRN *prn)
{
    MODEL aux;
    int t, err = 0;

    /* grow data set to accommodate new dependent var */
    if (dataset_add_vars(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[pdinfo->v - 1][t] = pmod->uhat[t];
    }

    /* replace the dependent var */
    list[1] = pdinfo->v - 1;

    aux = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (aux.errcode) {
	err = aux.errcode;
	fprintf(stderr, "auxiliary regression failed\n");
    } else {
	double trsq = aux.rsq * aux.nobs;
	int df = list[0] - pmod->list[0];
	double pval = chisq(trsq, df);

	aux.aux = aux_code;
	printmodel(&aux, pdinfo, opt, prn);

	pprintf(prn, "\n%s: TR^2 = %g,\n", _("Test statistic"), trsq);
	pprintf(prn, _("with p-value = prob(Chi-square(%d) > %g) = %g\n\n"), 
		df, trsq, pval);

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = new_test_on_model(pmod, (aux_code == AUX_SQ)?
				     GRETL_TEST_SQUARES : GRETL_TEST_LOGS);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_TR2);
		model_test_set_dfn(test, df);
		model_test_set_value(test, trsq);
		model_test_set_pvalue(test, chisq(trsq, df));
	    }
	}

	record_test_result(trsq, pval, "non-linearity");
    } 

    clear_model(&aux);

    return err;
}

/**
 * nonlinearity_test:
 * @pmod: pointer to original model.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @aux_code: %AUX_SQ for squares or %AUX_LOG for logs
 * @opt: if contains %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Run an auxiliary regression to test @pmod for non-linearity,
 * via the addition of either squares or logs of the original
 * indepdendent variables.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int nonlinearity_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		       int aux_code, gretlopt opt, PRN *prn) 
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *tmplist = NULL;
    const int orig_nvar = pdinfo->v; 
    int err = 0;

    if (!command_ok_for_model(ADD, pmod->ci)) {
	return E_NOTIMP;
    }

    if (pmod->ci == LOGISTIC || pmod->ci == LAD) {
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod->list, pdinfo, pmod->ID);
    if (err) {
	return err;
    }

    /* re-impose the sample that was in force when the original model
       was estimated */
    impose_model_smpl(pmod, pdinfo);

    /* add squares or logs */
    tmplist = augment_regression_list(pmod->list, aux_code, 
				      pZ, pdinfo);
    if (tmplist == NULL) {
	return E_ALLOC;
    } else if (tmplist[0] == pmod->list[0]) {
	/* no vars were added */
	if (aux_code == AUX_SQ) {
	    fprintf(stderr, "gretl: generation of squares failed\n");
	    err = E_SQUARES;
	} else if (aux_code == AUX_LOG) {
	    fprintf(stderr, "gretl: generation of logs failed\n");
	    err = E_LOGS;
	}
    }

    if (!err) {
	err = real_nonlinearity_test(pmod, tmplist, pZ, pdinfo, aux_code, 
				     opt, prn);
    }
	
    /* trash any extra variables generated (squares, logs) */
    dataset_drop_vars(pdinfo->v - orig_nvar, pZ, pdinfo);

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

    free(tmplist);

    return err;
}

/**
 * add_test:
 * @addvars: list of variables to add to original model.
 * @orig: pointer to original model.
 * @new: pointer to receive new model, with vars added.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: can contain %OPT_Q (quiet) to suppress printing
 * of the new model, %OPT_O to print covariance matrx.
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after adding the specified
 * variables.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int add_test (const int *addvars, MODEL *orig, MODEL *new, 
	      double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn)
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *tmplist = NULL;
    const int orig_nvar = pdinfo->v; 
    int save_test = 0;
    int err = 0;

    if (orig == NULL || orig->list == NULL) {
	return 1;
    }

    if (!command_ok_for_model(ADD, orig->ci)) {
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(orig->list, pdinfo, orig->ID);
    if (err) {
	return err;
    }

    /* create augmented regression list */
    tmplist = gretl_list_add(orig->list, addvars, &err);
    if (err) {
	return err;
    }

    /* impose as sample range the sample range in force
       when the original model was estimated */
    impose_model_smpl(orig, pdinfo);

    if (opt & OPT_S) {
	/* don't pass OPT_S to replicate_estimator() */
	save_test = 1;
	opt &= ~OPT_S;
    }

    /* Run augmented regression, matching the original estimation
       method; use OPT_Z to suppress the elimination of perfectly
       collinear variables.
    */
    *new = replicate_estimator(orig, &tmplist, pZ, pdinfo, 
			       (opt | OPT_Z), prn);

    if (new->errcode) {
	pprintf(prn, "%s\n", gretl_errmsg);
	err = new->errcode; 
    }

    if (!err) {

	new->aux = AUX_ADD;

	if (!(opt & OPT_Q) && new->ci != AR && new->ci != ARCH) {
	    printmodel(new, pdinfo, opt, prn);
	}

	if (new->nobs == orig->nobs) {
	    struct COMPARE cmp;
	    int *addlist;

	    addlist = gretl_list_diff_new(new->list, orig->list);
	    cmp = add_or_omit_compare(orig, new, 1, addlist);

	    if (save_test) {
		opt |= OPT_S;
	    }

	    gretl_print_compare(&cmp, addlist, orig, pdinfo, opt, prn);
	    free(addlist);
	}
    }

    /* trash any extra variables generated (squares, logs) */
    dataset_drop_vars(pdinfo->v - orig_nvar, pZ, pdinfo);

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

    free(tmplist);

    return err;
}

/**
 * omit_test:
 * @omitvars: list of variables to omit from original model.
 * @orig: pointer to original model.
 * @new: pointer to receive new model, with vars omitted.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt: can contain %OPT_Q (quiet) to suppress printing
 * of the new model, %OPT_O to print covariance matrx.
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after removing a list of 
 * specified variables.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test (const int *omitvars, MODEL *orig, MODEL *new, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *tmplist = NULL;
    int maxlag = 0;
    int save_test = 0;
    int err = 0;

    if (orig == NULL || orig->list == NULL) {
	return 1;
    }

    if (!command_ok_for_model(OMIT, orig->ci)) {
	return E_NOTIMP;
    }

    /* check that vars to omit have not been redefined */
    if ((err = list_members_replaced(orig->list, pdinfo, orig->ID))) {
	return err;
    }

    /* create list for test model */
    if (omitvars == NULL) {
	/* special: just drop the last variable */
	tmplist = gretl_list_omit_last(orig->list, &err);
    } else {
	tmplist = gretl_list_omit(orig->list, omitvars, &err);
    }

    if (tmplist == NULL) {
	return err;
    }

    /* impose as sample range the estimation range of the 
       original model */
    pdinfo->t1 = orig->t1;
    pdinfo->t2 = orig->t2;

    /* set the mask for missing obs within the sample range, based
       on the original model */
    set_reference_missmask(orig);

    if (orig->ci == AR) { 
	maxlag = orig->arinfo->arlist[orig->arinfo->arlist[0]];
    } else if (orig->ci == ARCH) {
	maxlag = orig->order;
    }

    pdinfo->t1 = orig->t1 - maxlag; /* FIXME: problem? */

    if (orig->ci == CORC || orig->ci == HILU) {
	pdinfo->t1 -= 1;
    }

    if (opt & OPT_S) {
	/* don't pass OPT_S to replicate estimator */
	save_test = 1;
	opt &= ~OPT_S;
    }

    *new = replicate_estimator(orig, &tmplist, pZ, pdinfo, opt, prn);

    if (new->errcode) {
	pprintf(prn, "%s\n", gretl_errmsg);
	err = new->errcode; 
    }

    if (!err) {
	if (orig->ci == LOGIT || orig->ci == PROBIT) {
	    new->aux = AUX_OMIT;
	}

	if (!(opt & OPT_Q) && orig->ci != AR && orig->ci != ARCH) {
	    printmodel(new, pdinfo, opt, prn); 
	}	

	if (new->nobs == orig->nobs && omitvars != NULL) {
	    struct COMPARE cmp;
	    int *omitlist;

	    omitlist = gretl_list_diff_new(orig->list, new->list);
	    cmp = add_or_omit_compare(orig, new, 0, omitlist);

	    if (save_test) {
		opt |= OPT_S;
	    }

	    gretl_print_compare(&cmp, omitlist, orig, pdinfo, opt, prn); 
	    free(omitlist);
	}

	if (orig->ci == LOGIT || orig->ci == PROBIT) {
	    new->aux = AUX_NONE;
	}
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

    /* reset the missing mask */
    set_reference_missmask(NULL);

    free(tmplist);

    return err;
}

static int ljung_box (int varno, int order, const double **Z, 
		      DATAINFO *pdinfo, double *lb)
{
    double *x, *y, *acf;
    int k, m, nobs, n = pdinfo->n; 
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int list[2];

    list[0] = 1;
    list[1] = varno;

    varlist_adjust_sample(list, &t1, &t2, Z);
    nobs = t2 - t1 + 1;

    x = malloc(n * sizeof *x);
    y = malloc(n * sizeof *y);
    acf = malloc((order + 1) * sizeof *acf);

    if (x == NULL || y == NULL || acf == NULL) {
	return E_ALLOC;
    }

    for (m=1; m<=order; m++) {
	for (t=t1+m; t<=t2; t++) {
	    k = t - (t1 + m);
	    x[k] = Z[varno][t];
	    y[k] = Z[varno][t-m];
	}
	acf[m] = gretl_corr(0, nobs - m - 1, x, y);
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

/**
 * reset_test:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if contains %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Carries out Ramsey's RESET test for model specification.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int reset_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn)
{
    int *newlist;
    MODEL aux;
    int i, t, v = pdinfo->v; 
    double RF;
    int err = 0;

    if (pmod->ci != OLS) return E_OLSONLY;

    gretl_model_init(&aux);

    if (pmod->ncoeff + 2 >= pdinfo->t2 - pdinfo->t1) {
	return E_DF;
    }

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
	double pval;

	aux.aux = AUX_RESET;
	printmodel(&aux, pdinfo, OPT_NONE, prn);
	RF = ((pmod->ess - aux.ess) / 2) / (aux.ess / aux.dfd);
	pval = fdist(RF, 2, aux.dfd);

	pprintf(prn, "\n%s: F = %f,\n", _("Test statistic"), RF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		2, aux.dfd, RF, pval);

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = new_test_on_model(pmod, GRETL_TEST_RESET);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_RESET);
		model_test_set_dfn(test, 2);
		model_test_set_dfd(test, aux.dfd);
		model_test_set_value(test, RF);
		model_test_set_pvalue(test, pval);
	    }	    
	}

	record_test_result(RF, pval, "RESET");
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
 * @opt: if flags include %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero. Gives TR^2 and LMF test statistics.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int autocorr_test (MODEL *pmod, int order, 
		   double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn)
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *newlist = NULL;
    MODEL aux;
    int i, t, n = pdinfo->n, v = pdinfo->v; 
    double trsq, LMF, lb, pval = 1.0;
    int err = 0;

    if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == LOGISTIC) 
	return E_NOTIMP;

    if (pmod->missmask != NULL) {
	return E_DATA;
    }

    if (dataset_is_panel(pdinfo)) {
	void *handle;
	int (*panel_autocorr_test)(MODEL *, int, 
				   double **, DATAINFO *, 
				   PRN *, ModelTest *);

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

    /* impose original sample range */
    impose_model_smpl(pmod, pdinfo);

    gretl_model_init(&aux);

    if (order <= 0) order = pdinfo->pd;

    if (pmod->ncoeff + order >= pdinfo->t2 - pdinfo->t1) {
	return E_DF;
    }

    newlist = malloc((pmod->list[0] + order + 1) * sizeof *newlist);

    if (newlist == NULL) {
	err = E_ALLOC;
    } else {
	newlist[0] = pmod->list[0] + order;
	for (i=2; i<=pmod->list[0]; i++) {
	    newlist[i] = pmod->list[i];
	}
	if (dataset_add_vars(1, pZ, pdinfo)) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add uhat to data set */
	for (t=0; t<n; t++) {
	    (*pZ)[v][t] = pmod->uhat[t];
	}
	strcpy(pdinfo->varname[v], "uhat");
	strcpy(VARLABEL(pdinfo, v), _("residual"));
	/* then lags of same */
	for (i=1; i<=order; i++) {
	    int lnum;

	    lnum = laggenr(v, i, pZ, pdinfo);
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

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = new_test_on_model(pmod, GRETL_TEST_AUTOCORR);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LMF);
		model_test_set_dfn(test, order);
		model_test_set_dfd(test, aux.nobs - pmod->ncoeff - order);
		model_test_set_order(test, order);
		model_test_set_value(test, LMF);
		model_test_set_pvalue(test, pval);
	    }	    
	}

	record_test_result(LMF, pval, "autocorrelation");
    }

    free(newlist);
    dataset_drop_vars(pdinfo->v - v, pZ, pdinfo); 
    clear_model(&aux); 

    if (pval < 0.05 && !gretl_model_get_int(pmod, "robust")) {
	autocorr_standard_errors(pmod, pZ, pdinfo, prn);
    }

    /* reset sample as it was */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

    return err;
}

/**
 * chow_test:
 * @line: command line for parsing.
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability (Chow test).
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test (const char *line, MODEL *pmod, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *chowlist = NULL;
    int newvars = pmod->list[0] - 1;
    int i, t, v = pdinfo->v;
    char chowdate[OBSLEN], s[OBSLEN];
    MODEL chow_mod;
    double F;
    int split = 0, err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    impose_model_smpl(pmod, pdinfo);

    gretl_model_init(&chow_mod);

    if (sscanf(line, "%*s %8s", chowdate) != 1) {
	err = E_PARSE;
    } else {
	split = dateton(chowdate, pdinfo) - 1;
	if (split <= 0 || split >= pdinfo->n) { 
	    err = E_SPLIT;
	}
    }

    if (!err) {
	/* take the original regression list, add a split dummy
	   and interaction terms. */
	if (pmod->ifc == 0) {
	    newvars++;
	}
	if (dataset_add_vars(newvars, pZ, pdinfo)) {
	    newvars = 0;
	    err = E_ALLOC;
	} else {
	    chowlist = malloc((pmod->list[0] + newvars + 1) * sizeof *chowlist);
	    if (chowlist == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	chowlist[0] = pmod->list[0] + newvars;
	for (i=1; i<=pmod->list[0]; i++) { 
	    chowlist[i] = pmod->list[i];
	}

	/* generate the split variable */
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[v][t] = (double) (t > split); 
	}
	strcpy(pdinfo->varname[v], "splitdum");
	strcpy(VARLABEL(pdinfo, v), _("dummy variable for Chow test"));
	chowlist[pmod->list[0] + 1] = v;

	/* and the interaction terms */
	for (i=1; i<newvars; i++) {
	    int orig = i + 1 + pmod->ifc;
	    int plo = pmod->list[orig];

	    for (t=0; t<pdinfo->n; t++) {
		if (model_missing(pmod, t)) {
		    (*pZ)[v+i][t] = NADBL;
		} else {
		    (*pZ)[v+i][t] = (*pZ)[v][t] * (*pZ)[plo][t];
		}
	    }
	    strcpy(s, pdinfo->varname[plo]); 
	    gretl_trunc(s, 5);
	    strcpy(pdinfo->varname[v+i], "sd_");
	    strcat(pdinfo->varname[v+i], s);
	    sprintf(VARLABEL(pdinfo, v+i), "splitdum * %s", pdinfo->varname[plo]);
	    chowlist[pmod->list[0] + 1 + i] = v + i;
	}

	chow_mod = lsq(chowlist, pZ, pdinfo, OLS, OPT_A, 0.0);

	if (chow_mod.errcode) {
	    err = chow_mod.errcode;
	    errmsg(err, prn);
	} else {
	    double pval;

	    chow_mod.aux = AUX_CHOW;
	    printmodel(&chow_mod, pdinfo, OPT_NONE, prn);
	    F = (pmod->ess - chow_mod.ess) * chow_mod.dfd / 
		(chow_mod.ess * newvars);
	    pval = fdist(F, newvars, chow_mod.dfd);
	    pprintf(prn, _("\nChow test for structural break at observation %s:\n"
		    "  F(%d, %d) = %f with p-value %f\n\n"), chowdate,
		    newvars, chow_mod.dfd, F, pval);

	    if (opt & OPT_S) {
		ModelTest *test;

		test = new_test_on_model(pmod, GRETL_TEST_CHOW);
		if (test != NULL) {
		    model_test_set_teststat(test, GRETL_STAT_F);
		    model_test_set_param(test, chowdate);
		    model_test_set_dfn(test, newvars);
		    model_test_set_dfd(test, chow_mod.dfd);
		    model_test_set_value(test, F);
		    model_test_set_pvalue(test, pval);
		}	  
	    }

	    record_test_result(F, pval, "Chow");
	}
	clear_model(&chow_mod);
    }

    /* clean up extra variables */
    dataset_drop_vars(newvars, pZ, pdinfo);
    free(chowlist);

    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

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
 * @opt: if flags include %OPT_S, save results of test to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for parameter stability (CUSUM test).
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int cusum_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn) 
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int n_est, i, j, t;
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
	    cum_mod = lsq(pmod->list, pZ, pdinfo, OLS, (OPT_C | OPT_A), 0.0);
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
	double pval = NADBL;

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
	    for (i=0; i<=j; i++) {
		W[j] += cresid[i];
	    }
	    W[j] /= sigma;
	    t = pmod->t1 + K + j;
	    ntodate(cumdate, t, pdinfo);
	    /* FIXME printing of number below? */
	    pprintf(prn, " %s %9.3f %s\n", cumdate, W[j],
		    (fabs(W[j]) > xx + (j+1)*yy)? "*" : "");
	}

	hct = (sqrt((double) (T-K)) * wbar) / sigma;
	pval = t_pvalue_2(hct, T-K-1);
	pprintf(prn, _("\nHarvey-Collier t(%d) = %g with p-value %.4g\n\n"), 
		T-K-1, hct, pval);

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = new_test_on_model(pmod, GRETL_TEST_CUSUM);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_HARVEY_COLLIER);
		model_test_set_dfn(test, T-K-1);
		model_test_set_value(test, hct);
		model_test_set_pvalue(test, pval);
	    }
	}

	record_test_result(hct, pval, "Harvey-Collier");

#ifdef ENABLE_NLS
        setlocale(LC_NUMERIC, "C");
#endif
	/* plot with 95% confidence bands, if not batch mode */
	if (!gretl_in_batch_mode() &&
	    gnuplot_init(PLOT_CUSUM, &fq) == 0) {
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

	    err = gnuplot_make_graph();
	}

#ifdef ENABLE_NLS
        setlocale(LC_NUMERIC, "");
#endif
    }

    /* restore sample */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;
    
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
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Tests the given pooled model for fixed and random effects.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int hausman_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn) 
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
	void (*panel_diagnostics)(MODEL *, double ***, DATAINFO *, 
				  gretlopt, PRN *);

	panel_diagnostics = get_plugin_function("panel_diagnostics", &handle);
	if (panel_diagnostics == NULL) {
	    return 1;
	}
	(*panel_diagnostics) (pmod, pZ, pdinfo, opt, prn);
	close_plugin(handle);
    }
    return 0;
}

/**
 * add_leverage_values_to_dataset:
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @m: matrix containing leverage values.
 * @flags: option flags: combination of %SAVE_LEVERAGE, %SAVE_INFLUENCE,
 * and %SAVE_DFFITS.
 *
 * Adds to the working dataset one or more series calculated by
 * the gretl test for leverage/influence of data points.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int add_leverage_values_to_dataset (double ***pZ, DATAINFO *pdinfo,
				    gretl_matrix *m, unsigned char flags)
{
    int t1, t2;
    int addvars = 0;

    if (flags & SAVE_LEVERAGE) addvars++;
    if (flags & SAVE_INFLUENCE) addvars++;
    if (flags & SAVE_DFFITS) addvars++;

    if (dataset_add_vars(addvars, pZ, pdinfo)) {
	strcpy(gretl_errmsg, _("Out of memory adding series"));
	return 1;
    }

    t1 = gretl_matrix_get_int(m);
    t2 = t1 + gretl_matrix_rows(m);

    /* add leverage? */
    if (flags & SAVE_LEVERAGE) {
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
    if (flags & SAVE_INFLUENCE) {
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
    if (flags & SAVE_DFFITS) {
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
 * @opt: if %OPT_S, add calculated series to data set.
 * @prn: gretl printing struct.
 *
 * Tests the data used in the given model for points with
 * high leverage and influence on the estimates.
 * 
 * Returns: 0 on successful completion, error code on error.
 *
 */

int leverage_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn)
{
    void *handle;
    gretl_matrix *(*model_leverage) (const MODEL *, double ***, 
				     const DATAINFO *, PRN *, int);
    gretl_matrix *m;
    int err = 0;

    if (pmod->ci != OLS) return E_OLSONLY;

    model_leverage = get_plugin_function("model_leverage", &handle);
    if (model_leverage == NULL) {
	return 1;
    }

    m = (*model_leverage)(pmod, pZ, pdinfo, prn, 0);
    if (m == NULL) {
	err = 1;
    } else {
	if (opt & OPT_S) {
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

int make_mp_lists (const int *list, const char *str,
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

int mp_ols (const int *list, const char *pos,
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

static int varmatch (const int *sumvars, int test)
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
			int *tmplist, const int *sumvars, int newv)
{
    int repl = 0;
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

    return (nnew == 0)? testcoeff : -1;
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

int sum_test (const int *sumvars, MODEL *pmod, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn)
{
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *tmplist = NULL;
    const int oldv = pdinfo->v;
    int testcoeff;
    MODEL summod;
    PRN *nullprn;
    int err = 0;

    if (sumvars[0] < 2) {
	pprintf(prn, _("Invalid input\n"));
	return E_DATA;
    }

    if (!command_ok_for_model(COEFFSUM, pmod->ci)) {
	return E_NOTIMP;
    }

    tmplist = malloc((pmod->list[0] + 1) * sizeof *tmplist);
    if (tmplist == NULL) {
	return E_ALLOC;
    }

    if (dataset_add_vars(sumvars[0] - 1, pZ, pdinfo)) {
	free(tmplist);
	return E_ALLOC;
    }

    nullprn = gretl_print_new(GRETL_PRINT_NULL);

    testcoeff = make_sum_test_list(pmod, *pZ, pdinfo, tmplist, sumvars, oldv);

    if (testcoeff < 0) {
	pprintf(prn, _("Invalid input\n"));
	free(tmplist);
	dataset_drop_vars(pdinfo->v - oldv, pZ, pdinfo);
	return E_DATA;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    impose_model_smpl(pmod, pdinfo);

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
		double tval, pval;

		pprintf(prn, "   %s = %g\n", _("Standard error"),
			summod.sderr[testcoeff - 2]);
		tval = summod.coeff[testcoeff - 2] / 
		    summod.sderr[testcoeff - 2];
		pval = coeff_pval(&summod, tval, summod.dfd);
		pprintf(prn, "   t(%d) = %g ", summod.dfd, tval);
		pprintf(prn, _("with p-value = %g\n"), pval);
		record_test_result(tval, pval, "sum");
	    }
	}
    }

    free(tmplist);
    clear_model(&summod);

    dataset_drop_vars(pdinfo->v - oldv, pZ, pdinfo);
    gretl_print_destroy(nullprn);

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

    return err;
}
