#define VO_DEBUG 0

const int *gretl_VAR_get_exo_list (const GRETL_VAR *var)
{
    return var->xlist;
}

/* Based on the specification stored in the VAR struct, reconstitute
   the list that was intially passed to the gretl_VAR() function to
   set up the system.
*/

static int *rebuild_VAR_list (const GRETL_VAR *orig, int *exolist, int *err)
{
    int *list = NULL;
    int lsep = (exolist[0] > 0);
    int i, j = 1;

    list = gretl_list_new(orig->neqns + exolist[0] + lsep);
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<orig->neqns; i++) {
	list[j++] = orig->models[i]->list[1];
    }

    if (lsep) {
	list[j++] = LISTSEP;
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
    int i, df, err = 0;

#if VO_DEBUG
    fprintf(stderr, "gretl_VAR_real_omit_test: about to diff lists\n");
    printlist(exolist0, "exolist0");
    printlist(exolist1, "exolist1");
#endif

    omitlist = gretl_list_diff_new(exolist0, exolist1, 1);
    if (omitlist == NULL) {
	return E_ALLOC;
    }

    LR = orig->T * (new->ldet - orig->ldet);
    df = orig->neqns * omitlist[0];
    pval = chisq_cdf_comp(LR, df);
    
    pputs(prn, _("\n  Null hypothesis: the regression parameters are "
		 "zero for the variables\n\n"));
    for (i=1; i<=omitlist[0]; i++) {
	pprintf(prn, "    %s\n", pdinfo->varname[omitlist[i]]);	
    }

    pprintf(prn, "\n  %s: %s(%d) = %g, ", _("Test statistic"), 
	    _("Chi-square"), df, LR);
    pprintf(prn, _("with p-value = %g\n\n"), pval);

    free(omitlist);

    return err;
}

/**
 * gretl_VAR_omit_test:
 * @omitvars: list of variables to omit from original model.
 * @orig: pointer to original VAR.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Re-estimates a given VAR after removing the variables
 * specified in @omitvars, and reports per-equation F-tests
 * and system-wide LR tests for the null hypothesis that
 * the omitted variables have zero parameters.
 * 
 * Returns: restricted VAR on sucess, %NULL on error.
 */

GRETL_VAR *gretl_VAR_omit_test (const int *omitvars, const GRETL_VAR *orig, 
				double ***pZ, DATAINFO *pdinfo, 
				PRN *prn, int *err)
{
    GRETL_VAR *var = NULL;
    gretlopt opt = OPT_NONE;
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    const int *exolist = NULL;
    int *tmplist = NULL;
    int *varlist = NULL;
    int c0, c1;

    *err = 0;

    if (orig == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (omitvars == NULL || omitvars[0] == 0) {
	*err = E_PARSE;
	return NULL;
    }

    /* recreate the exog vars list for original VAR */
    exolist = gretl_VAR_get_exo_list(orig);
    if (exolist == NULL) {
	*err = E_DATA;
	return NULL;
    }

#if VO_DEBUG
    printlist(exolist, "original exolist");
#endif

    c0 = gretl_list_const_pos(exolist, 1, (const double **) *pZ, pdinfo);
    if (c0 > 0) {
	c1 = !gretl_list_const_pos(omitvars, 1, (const double **) *pZ, pdinfo);
    } else {
	c1 = 0;
    }

    /* create exogenous vars list for test VAR */
    tmplist = gretl_list_omit(exolist, omitvars, 1, err);
    if (tmplist == NULL) {
	goto bailout;
    }

#if VO_DEBUG
    fprintf(stderr, "co = %d, c1 = %d\n", c0, c1);
    printlist(tmplist, "exog vars list for test VAR");
#endif

    /* recreate full input VAR list for test VAR */
    varlist = rebuild_VAR_list(orig, tmplist, err);
    if (varlist == NULL) {
	goto bailout;
    }

#if VO_DEBUG
    printlist(varlist, "full list for test VAR");
#endif

    /* If the original VAR did not include a constant, we need to
       pass OPT_N to the test VAR to prevent the addition of a
       constant.  We also need to pass OPT_N in case the constant was
       present originally but is now to be omitted.
    */
    if (c0 == 0 || c1 == 0) {
	opt = OPT_N;
    }

    /* impose as sample range the estimation range of the 
       original model */
    pdinfo->t1 = orig->t1;
    pdinfo->t2 = orig->t2;

    var = gretl_VAR(orig->order, varlist, pZ, pdinfo, opt, prn, err);

    /* now, if var is non-NULL, do the actual test(s) */
    if (var != NULL) {
	*err = gretl_VAR_real_omit_test(orig, exolist, var, tmplist,
					pdinfo, prn);
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

 bailout:

    free(tmplist);
    free(varlist);

    return var;
}
