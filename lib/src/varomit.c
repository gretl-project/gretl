

static int *gretl_VAR_get_exo_list (const GRETL_VAR *var, int *err)
{
    int *vlist, *elist;
    int nendo, nexo;
    int i, j;

    if (var->models == NULL) {
	*err = E_DATA;
	return NULL;
    }

    vlist = var->models[0]->list;

    /* the endogenous vars start in position 2, and there are
       order * neqns such terms */

    nendo = var->order * var->neqns;
    nexo = vlist[0] - 1 - nendo;
    if (nexo == 0) {
	*err = E_DATA;
	return NULL;
    }

    elist = gretl_list_new(nexo);
    if (elist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    j = 1;
    for (i=nendo+2; i<=vlist[0]; i++) {
	elist[j++] = vlist[i];
    }

    return elist;
}

static int *rebuild_VAR_list (const GRETL_VAR *orig, int *exolist, int *err)
{
    int *list = NULL;
    int i, j = 1;

    list = gretl_list_new(orig->neqns + exolist[0]);
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<orig->neqns; i++) {
	list[j++] = orig->models[i]->list[1];
    }

    for (i=1; i<=exolist[0]; i++) {
	list[j++] = exolist[i];
    }    

    return list;
}

static int gretl_VAR_real_omit_test (const GRETL_VAR *orig,
				     const int *exolist0,
				     const GRETL_VAR *new,
				     const int *exolist1,
				     const DATAINFO *pdinfo,
				     PRN *prn)
{
    int *omitlist;
    double LR, pval;
    int df, err = 0;

    omitlist = gretl_list_diff_new(exolist0, exolist1);
    if (omitlist == NULL) {
	return E_ALLOC;
    }

    LR = orig->T * (new->ldet - orig->ldet);
    df = orig->neqns * omitlist[0];
    pval = chisq(LR, df);
    
    pputs(prn, _("\n  Null hypothesis: the regression parameters are "
		 "zero for the variables\n\n"));
    for (i=1; i<=omitlist[0]; i++) {
	pprintf(prn, "    %s\n", pdinfo->varname[omitlist[i]]);	
    }

    pprintf(prn, "\n  %s:%s%s(%d) = %g, ", _("Test statistic"), 
	    _("Chi-square"), df, LR);
    pprintf(prn, _("with p-value = %g\n\n"), pval);

    free(omitlist);

    return err;
}

/**
 * gretl_VAR_omit_test:
 * @omitvars: list of variables to omit from original model.
 * @var: pointer to original VAR.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 *
 * Re-estimates a given VAR after removing the variables
 * specified in @omitvars, and reports per-equation F-tests
 * and system-wide LR tests for the null hypothesis that
 * the omitted variables have zero parameters.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int gretl_VAR_omit_test (const int *omitvars, const GRETL_VAR *orig, 
			 double ***pZ, DATAINFO *pdinfo, 
			 PRN *prn)
{
    GRETL_VAR *var = NULL;
    gretlopt opt = OPT_NONE; /* FIXME: need to pickup opts from orig */
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *exolist = NULL;
    int *tmplist = NULL;
    int *varlist = NULL;
    int err = 0;

    if (orig == NULL) {
	return 1;
    }

    if (omitvars == NULL || omitvars[0] == 0) {
	return E_PARSE;
    }

    /* recreate the exog vars list for original VAR */
    exolist = gretl_VAR_get_exo_list(orig, &err);
    if (exolist == NULL) {
	return err;
    }

    /* create exogenous vars list for test VAR */
    tmplist = gretl_list_omit(exolist, omitvars, &err);
    if (tmplist == NULL) {
	goto bailout;
    }

    /* recreate full input VAR list for test VAR */
    varlist = rebuild_VAR_list(orig, tmplist, &err);
    if (varlist == NULL) {
	goto bailout;
    }

    /* impose as sample range the estimation range of the 
       original model */
    pdinfo->t1 = orig->t1;
    pdinfo->t2 = orig->t2;

    var = gretl_VAR(orig->order, list, pZ, pdinfo, opt, NULL, &err);

    /* now, if var is non-NULL, do the actual test(s) */
    if (var != NULL) {
	err = gretl_VAR_real_omit_test(orig, exolist, var, tmplist,
				       pdinfo, prn);
	gretl_VAR_free(var);
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

 bailout:

    free(exolist);
    free(tmplist);
    free(varlist);

    return err;
}
