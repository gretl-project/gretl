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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* panel data plugin for gretl */

#include "libgretl.h"
#include "f2c.h"
#include "clapack_double.h"

typedef struct {
    int ns;
    double sigma_e;
    double H;
    double *bdiff;
    double *sigma;
} hausman_t;

/* #define PDEBUG 1 */

/* .................................................................. */

static void print_panel_coeff (MODEL *pmod, MODEL *panelmod,
			       DATAINFO *pdinfo, int i, 
			       PRN *prn)
{
    char numstr[18];

    sprintf(numstr, "(%.5g)", panelmod->sderr[i]);
    pprintf(prn, "%*s: %14.5g %15s\n", VNAMELEN,
	    pdinfo->varname[pmod->list[i+2]],
	    panelmod->coeff[i], numstr);
}

/* .................................................................. */

static double group_means_variance (MODEL *pmod, 
				    double **Z, DATAINFO *pdinfo,
				    double ***groupZ, DATAINFO **ginfo,
				    int nunits, int T) 
{
    int i, j, k, t, start, *list;
    double xx;
    MODEL meanmod;

#ifdef PDEBUG
    fprintf(stderr, "group_means_variance: creating dataset\n"
	    " nvar=%d, nobs=%d, *groupZ=%p\n", pmod->list[0],
	    nunits, (void *)*groupZ);
#endif

    *ginfo = create_new_dataset(groupZ, pmod->list[0], nunits, 0);
    if (*ginfo == NULL) return NADBL;

    list = malloc((pmod->list[0] + 1) * sizeof *list);
    if (list == NULL) {
	clear_datainfo(*ginfo, CLEAR_FULL);
	free(*ginfo);
	return NADBL;
    }

#ifdef PDEBUG
    fprintf(stderr, "gmv: *groupZ=%p\n", (void *) *groupZ);
#endif

    list[0] = pmod->list[0];
    k = 1;
    for (j=1; j<=list[0]; j++) { /* the variables */
	if (pmod->list[j] == 0) {
	    list[j] = 0;
	    continue;
	}
	list[j] = k;
	start = 0;
	for (i=0; i<nunits; i++) { /* the observations */
	    xx = 0.0;
	    if (pdinfo->time_series == 2) {
		for (t=start; t<start+T; t++) {
		    xx += Z[pmod->list[j]][t];
		}
		start += T;
	    } else {
		for (t=start; t<pdinfo->n; t += nunits) {
		    xx += Z[pmod->list[j]][t];
		}
		start++;
	    }
	    xx /= (double) T;
	    (*groupZ)[k][i] = xx;
	}
	k++;
    }

#ifdef PDEBUG
    fprintf(stderr, "group_means_variance: about to run OLS\n");
    printlist(list, NULL);
    fprintf(stderr, "*groupZ=%p, ginfo=%p\n", (void *)*groupZ, (void *)*ginfo);
#endif

    meanmod = lsq(list, groupZ, *ginfo, OLS, OPT_A, 0.0);
#ifdef PDEBUG
    fprintf(stderr, "gmv: lsq errcode was %d\n", meanmod.errcode);
#endif
    if (meanmod.errcode) xx = NADBL;
    else xx = meanmod.sigma * meanmod.sigma;

    clear_model(&meanmod, NULL);
    free(list);
#ifdef PDEBUG
    fprintf(stderr, "gmv: done freeing stuff\n");
#endif

    return xx;
}

/* .................................................................. */

static void vcv_slopes (hausman_t *haus, MODEL *pmod, int nunits, 
			int subt)
{
    int i, j, k = 0, idx;

    for (i=0; i<haus->ns; i++) {
	for (j=i; j<haus->ns; j++) {
	    idx = ijton(i+1, j+1, pmod->ncoeff);
#ifdef PDEBUG
	    fprintf(stderr, "setting sigma[%d] using vcv[%d] (%d, %d) = %g\n",
		    k, idx, i+2, j+2, pmod->vcv[idx]);
#endif
	    if (subt) haus->sigma[k++] -= pmod->vcv[idx];
	    else haus->sigma[k++] = pmod->vcv[idx];
	}
    }
}

/* .................................................................. */

static int bXb (hausman_t *haus)
{
    char uplo = 'L'; 
    integer nrhs = 1;
    integer n = haus->ns;
    integer ldb = n;
    integer info = 0;
    integer *ipiv;
    double *x;
    int i;

    /* make a copy of haus->bdiff first */
    x = malloc(n * sizeof *x);
    if (x == NULL) return E_ALLOC;

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) {
	free(x);
	return E_ALLOC;
    }

    for (i=0; i<haus->ns; i++) {
	x[i] = haus->bdiff[i];
    }

    /* solve for X-inverse * b */
    dspsv_(&uplo, &n, &nrhs, haus->sigma, ipiv, x, &ldb, &info);
    if (info > 0) {
	fprintf(stderr, "Hausman sigma matrix is singular\n");
    } else if (info < 0) {
	fprintf(stderr, "Illegal entry in Hausman sigma matrix\n");
    } else {
	haus->H = 0.0;
	for (i=0; i<haus->ns; i++) {
	    haus->H += x[i] * haus->bdiff[i];
	}
    }

    free(x);
    free(ipiv);

    return info;
}

/* .................................................................. */

static double LSDV (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		    int nunits, int T, hausman_t *haus, PRN *prn) 
{
    int i, t, oldv = pdinfo->v, start;
    int dvlen = pmod->list[0] + nunits;
    int *dvlist;
    double var, F;
    MODEL lsdv;

    /* We can be assured there's an intercept in the original
       regression */

    dvlist = malloc(dvlen * sizeof *dvlist);
    if (dvlist == NULL) return NADBL;
    if (dataset_add_vars(nunits - 1, pZ, pdinfo)) {
	free(dvlist);
	return NADBL;
    }

    start = 0;
    for (i=0; i<nunits-1; i++) {
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[oldv+i][t] = 0.0;
	}
	if (pdinfo->time_series == STACKED_TIME_SERIES) {
	    for (t=start; t<start+T; t++) {
		(*pZ)[oldv+i][t] = 1.0;
	    }
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) {
		(*pZ)[oldv+i][t] = 1.0;
	    }
	    start++;
	}
    }

    dvlist[0] = dvlen - 1;

    for (i=1; i<=pmod->list[0]; i++) {
	dvlist[i] = pmod->list[i];
    }
    for (i=1; i<nunits; i++) {
	dvlist[pmod->list[0] + i] = oldv + i - 1;
    }

#ifdef PDEBUG
    fprintf(stderr, "LSDV: about to run OLS\n");
#endif

    lsdv = lsq(dvlist, pZ, pdinfo, OLS, OPT_A, 0.0);

    if (lsdv.errcode) {
	var = NADBL;
	pputs(prn, _("Error estimating fixed effects model\n"));
	errmsg(lsdv.errcode, prn);
    } else {
	haus->sigma_e = lsdv.sigma;
	var = lsdv.sigma * lsdv.sigma;
	pputs(prn, 
	      _("                          Fixed effects estimator\n"
		"          allows for differing intercepts by cross-sectional "
		"unit\n"
		"         (slope standard errors in parentheses, a_i = "
		"intercepts)\n\n"));

	/* skip the constant in this printing */
	for (i=1; i<pmod->list[0] - 1; i++) {
	    print_panel_coeff(&lsdv, &lsdv, pdinfo, i, prn);
	    haus->bdiff[i-1] = lsdv.coeff[i];
	}

	for (i=0; i<nunits; i++) {
	    /* print per-unit intercept estimates */
	    char dumstr[VNAMELEN];
	    double x;

	    if (i == nunits - 1) x = lsdv.coeff[0];
	    else x = lsdv.coeff[i + pmod->list[0] - 1] + lsdv.coeff[0];
	    sprintf(dumstr, "a_%d", i + 1);
	    pprintf(prn, "%*s: %14.4g\n", VNAMELEN, dumstr, x);
	}

	pprintf(prn, _("\nResidual variance: %g/(%d - %d) = %g\n"), 
		lsdv.ess, pdinfo->n, lsdv.ncoeff, var);
	F = (pmod->ess - lsdv.ess) * lsdv.dfd /
	    (lsdv.ess * (nunits - 1.0));
	pprintf(prn, _("Joint significance of unit dummy variables:\n"
		       " F(%d, %d) = %g with p-value %g\n"), nunits - 1,
		lsdv.dfd, F, fdist(F, nunits - 1, lsdv.dfd));
	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the pooled OLS model\nis adequate, in favor of the fixed "
		     "effects alternative.)\n\n"));
	makevcv(&lsdv);
	vcv_slopes(haus, &lsdv, nunits, 0);
    }

    clear_model(&lsdv, NULL);
    dataset_drop_vars(nunits - 1, pZ, pdinfo);
    free(dvlist);

    return var;
}

/* .................................................................. */

static int random_effects (MODEL *pmod, double **Z, DATAINFO *pdinfo, 
			   double **groupZ, double theta, int nunits, int T, 
			   hausman_t *haus, PRN *prn)
{
    double **reZ;
    DATAINFO *reinfo;
    MODEL remod;
    int *relist;
    int i, j, k, t, err = 0;

    reinfo = create_new_dataset(&reZ, pmod->list[0], pdinfo->n, 0);
    if (reinfo == NULL) return E_ALLOC;

    relist = malloc((pmod->list[0] + 1) * sizeof *relist);
    if (relist == NULL) {
	free_Z(reZ, reinfo);
	clear_datainfo(reinfo, CLEAR_FULL);
	free(reinfo);
	return E_ALLOC;
    }

    relist[0] = pmod->list[0];

    /* create transformed variables */
    k = 1;
    for (i=1; i<=relist[0]; i++) {
	if (pmod->list[i] == 0) {
	    relist[i] = 0;
	    continue;
	}
	relist[i] = k;
	j = 0;
	if (pdinfo->time_series == STACKED_TIME_SERIES) { 
	    for (t=0; t<pdinfo->n; t++) {
		if (t && (t % T == 0)) j++; 
		reZ[k][t] = Z[pmod->list[i]][t] 
		    - theta * groupZ[k][j];
	    }
	} else { /* stacked cross sections */
	    for (t=0; t<pdinfo->n; t++) {
		if (t && t % nunits == 0) j = 0; /* FIXME ?? */
		reZ[k][t] = Z[pmod->list[i]][t] 
		    - theta * groupZ[k][j];
		j++;
	    }
	}
	k++;
    }

    for (t=0; t<pdinfo->n; t++) reZ[0][t] = 1.0 - theta;

#ifdef PDEBUG
    fprintf(stderr, "random_effects: about to run OLS\n");
#endif

    remod = lsq(relist, &reZ, reinfo, OLS, OPT_A, 0.0);
    if ((err = remod.errcode)) {
	pputs(prn, _("Error estimating random effects model\n"));
	errmsg(err, prn);
    } else {
	pputs(prn,
	      _("                         Random effects estimator\n"
		"           allows for a unit-specific component to the "
		"error term\n"
		"                     (standard errors in parentheses)\n\n"));
	for (i=0; i<relist[0] - 1; i++) {
	    print_panel_coeff(pmod, &remod, pdinfo, i, prn);
	    if (i > 0) haus->bdiff[i-1] -= remod.coeff[i];
	}
	makevcv(&remod);
	vcv_slopes(haus, &remod, nunits, 1);
    }

    clear_model(&remod, NULL);
    free_Z(reZ, reinfo);
    clear_datainfo(reinfo, CLEAR_FULL);
    free(reinfo);
    free(relist);    

    return err;
}

/* .................................................................. */

static int breusch_pagan_LM (MODEL *pmod, DATAINFO *pdinfo, 
			     int nunits, int T, PRN *prn)
{
    double *ubar, LM, eprime = 0.0;
    int i, t, start = 0;

    ubar = malloc(nunits * sizeof *ubar);
    if (ubar == NULL) return E_ALLOC;

    for (i=0; i<nunits; i++) {
	ubar[i] = 0.0;
	if (pdinfo->time_series == STACKED_TIME_SERIES) {
	    for (t=start; t<start+T; t++) 
		ubar[i] += pmod->uhat[t];
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) 
		ubar[i] += pmod->uhat[t];
	    start++;
	}
	ubar[i] /= (double) T;
	eprime += ubar[i] * ubar[i];
    }
#ifdef PDEBUG
    fprintf(stderr,  "breusch_pagan: found ubars\n");
#endif

    pputs(prn, _("\nMeans of pooled OLS residuals for cross-sectional "
		 "units:\n\n"));
    for (i=0; i<nunits; i++) {
	pprintf(prn, _(" unit %2d: %13.5g\n"), 
		i + 1, ubar[i]);
    }
    free(ubar);

    LM = (double) pdinfo->n/(2.0*(T - 1.0)) * 
	pow((T * T * eprime/pmod->ess) - 1.0, 2);
    pprintf(prn, _("\nBreusch-Pagan test statistic:\n"
		   " LM = %g with p-value = prob(chi-square(1) > %g) = %g\n"), 
	    LM, LM, chisq(LM, 1));
    pputs(prn, _("(A low p-value counts against the null hypothesis that "
		 "the pooled OLS model\nis adequate, in favor of the random "
		 "effects alternative.)\n\n"));
    return 0;
}

/* .................................................................. */

static int do_hausman_test (hausman_t *haus, PRN *prn)
{
#ifdef PDEBUG
    int i, ns = haus->ns;
    int nterms = (ns * ns + ns) / 2;

    for (i=0; i<ns; i++)
 	pprintf(prn, "b%d_FE - beta%d_RE = %g\n", i, i, haus->bdiff[i]);
    pprintf(prn, "\n");

    for (i=0; i<nterms; i++)
 	pprintf(prn, "vcv_diff[%d] = %g\n", i, haus->sigma[i]);
#endif

    if (bXb(haus)) { 
	pputs(prn, _("Error attempting to invert vcv difference matrix\n"));
	return 1;
    }
    if (haus->H < 0) 
	pputs(prn, _("\nHausman test matrix is not positive definite (this "
		     "result may be treated as\n\"fail to reject\" the random effects "
		     "specification).\n"));
    else {
	pprintf(prn, _("\nHausman test statistic:\n"
		       " H = %g with p-value = prob(chi-square(%d) > %g) = %g\n"),
		haus->H, haus->ns, haus->H, chisq(haus->H, haus->ns));
	pputs(prn, _("(A low p-value counts against the null hypothesis that "
		     "the random effects\nmodel is consistent, in favor of the fixed "
		     "effects model.)\n"));
    }

    return 0;
}

/* .................................................................. */

int panel_diagnostics (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		       PRN *prn)
{
    int nunits, ns, T;
    double var1, var2, theta;
    double **groupZ = NULL;
    DATAINFO *ginfo = NULL;
    hausman_t haus;

    if (pmod->ifc == 0) return 1;

    if (get_panel_structure(pdinfo, &nunits, &T))
	return 1;

    if (nunits > pmod->ncoeff) {
	ns = haus.ns = pmod->ncoeff - 1;
	haus.bdiff = malloc(haus.ns * sizeof *haus.bdiff);
	if (haus.bdiff == NULL) return E_ALLOC;
	haus.sigma = malloc(((ns * ns + ns) / 2) * sizeof *haus.sigma);
	if (haus.sigma == NULL) return E_ALLOC; 
    }   
    
    pprintf(prn, _("      Diagnostics: assuming a balanced panel with %d "
		   "cross-sectional units\n "
		   "                        observed over %d periods\n\n"), 
	    nunits, T);

    var2 = LSDV(pmod, pZ, pdinfo, nunits, T, &haus, prn);

#ifdef PDEBUG
    fprintf(stderr, "panel_diagnostics: LSDV gave %g\n", var2);
#endif

    breusch_pagan_LM(pmod, pdinfo, nunits, T, prn);

#ifdef PDEBUG
    fprintf(stderr, "panel_diagnostics: done breusch_pagan_LM()\n");
#endif
    
    if (nunits > pmod->ncoeff && var2 > 0) {
	var1 = group_means_variance(pmod, *pZ, pdinfo, 
				    &groupZ, &ginfo, nunits, T);
	if (var1 < 0) 
	    pputs(prn, _("Couldn't estimate group means regression\n"));
	else {
	    pprintf(prn, _("Residual variance for group means "
			   "regression: %g\n\n"), var1);    
	    theta = 1.0 - sqrt(var2 / (T * var1));
	    random_effects(pmod, *pZ, pdinfo, groupZ, theta, nunits, T, 
			   &haus, prn);
	    do_hausman_test(&haus, prn);
	}
	free_Z(groupZ, ginfo);
	clear_datainfo(ginfo, CLEAR_FULL);
	free(ginfo);
	free(haus.bdiff);
	free(haus.sigma);
    }

    return 0;
}

/* .................................................................. */

static void panel_copy_var (double **targZ, DATAINFO *targinfo, int targv,
			    double *src, DATAINFO *srcinfo, int srcv,
			    int order)
{
    int t, j = 0;

    for (t=srcinfo->t1; t<=srcinfo->t2; t++) {
	if (t % srcinfo->pd >= order) {
	    targZ[targv][j++] = src[t];
	} 
    }

    if (srcv == -1) {
	strcpy(targinfo->varname[targv], "uhat");
	strcpy(VARLABEL(targinfo, targv), _("residual"));
    } else {
	strcpy(targinfo->varname[targv], srcinfo->varname[srcv]);
	strcpy(VARLABEL(targinfo, targv), VARLABEL(srcinfo, srcv));
    }
}

static void make_reduced_data_info (DATAINFO *targ, DATAINFO *src, int order)
{
    targ->pd = src->pd - order;
    ntodate(targ->stobs, src->t1 + order, src);
    targ->sd0 = obs_str_to_double(targ->stobs); 
    targ->time_series = src->time_series;
}

static void panel_lag (double **tmpZ, DATAINFO *tmpinfo, 
		       double *src, DATAINFO *srcinfo, 
		       int v, int order, int lag)
{
    int t, j = 0;

    for (t=srcinfo->t1; t<=srcinfo->t2; t++) {
	if (t % srcinfo->pd >= order) {
	    tmpZ[v][j++] = src[t - lag];
	}
    }

    sprintf(tmpinfo->varname[v], "uhat_%d", lag);
    *VARLABEL(tmpinfo, v) = 0;
}

/* - do some sanity checks
   - create a local copy of the required portion of the data set,
     skipping the obs that will be missing
   - copy in the lags of uhat
   - estimate the aux model
   - destroy the temporary data set
*/

int panel_autocorr_test (MODEL *pmod, int order, 
			 double **Z, DATAINFO *pdinfo, 
			 PRN *prn, GRETLTEST *test)
{
    int *aclist;
    double **tmpZ;
    DATAINFO *tmpinfo;
    MODEL aux;
    double trsq, LMF;
    int i, nv, nunits, nobs, err = 0;
    int sn = pdinfo->t2 - pdinfo->t1 + 1;

    /* basic checks */
    if (order <= 0) order = 1;
    if (order > pdinfo->pd - 1) return E_DF;
    if (pmod->ncoeff + order >= sn) return E_DF;

    if (pdinfo->time_series != STACKED_TIME_SERIES ||
	!balanced_panel(pdinfo)) { 
        return E_DATA;
    }

    /* get number of cross-sectional units */
    nunits = sn / pdinfo->pd;

    /* we lose "order" observations for each unit */
    nobs = sn - nunits * order;

    /* the required number of variables */
    nv = pmod->list[0] + order;

    /* create temporary reduced dataset */
    tmpinfo = create_new_dataset(&tmpZ, nv, nobs, 0);
    if (tmpinfo == NULL) return E_ALLOC;

    make_reduced_data_info(tmpinfo, pdinfo, order);

#ifdef PDEBUG
    fprintf(stderr, "Created data set, n=%d, pd=%d, vars=%d, stobs='%s'\n", 
	    tmpinfo->n, tmpinfo->pd, tmpinfo->v, tmpinfo->stobs);
#endif

    /* allocate the auxiliary regression list */
    aclist = malloc((nv + 1) * sizeof *aclist);
    if (aclist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	int k;

	aclist[0] = pmod->list[0] + order;
	/* copy model uhat to position 1 in temp data set */
	aclist[1] = 1;
	panel_copy_var(tmpZ, tmpinfo, 1,
		       &pmod->uhat[0], pdinfo, -1,
		       order);
	/* copy across the original indep vars, making
	   the new regression list while we're at it */
	k = 2;
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == 0) { /* the constant */
		aclist[i] = 0;
	    } else {
		aclist[i] = k;
		panel_copy_var(tmpZ, tmpinfo, k, 
			       &Z[pmod->list[i]][0], pdinfo, pmod->list[i], 
			       order);
		k++;
	    }
	}
    }

    if (!err) {
	int v = pmod->list[0] - 1;

	/* add lags of uhat to temp data set */
	for (i=1; i<=order; i++) {
	    panel_lag(tmpZ, tmpinfo, &pmod->uhat[0], 
		      pdinfo, v + i, order, i);
	    aclist[v + i + 1] = v + i;
	}
    }

    if (!err) {
	aux = lsq(aclist, &tmpZ, tmpinfo, OLS, OPT_D | OPT_A, 0.0);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    } 

    if (!err) {
	aux.aux = AUX_AR;
	aux.order = order;
	printmodel(&aux, tmpinfo, prn);
	trsq = aux.rsq * aux.nobs;
	LMF = (aux.rsq/(1.0 - aux.rsq)) * 
	    (aux.nobs - pmod->ncoeff - order)/order; 

	pprintf(prn, "\n%s: LMF = %f,\n", _("Test statistic"), LMF);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"), 
		order, aux.nobs - pmod->ncoeff - order, LMF,
		fdist(LMF, order, aux.nobs - pmod->ncoeff - order));

	pprintf(prn, "\n%s: TR^2 = %f,\n", 
		_("Alternative statistic"), trsq);
	pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", 	_("with p-value"), 
		_("Chi-square"), order, trsq, chisq(trsq, order));

	if (test != NULL) {
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

    free(aclist);
    clear_model(&aux, tmpinfo); 

    free_Z(tmpZ, tmpinfo);
    clear_datainfo(tmpinfo, CLEAR_FULL);
    free(tmpinfo);

    return err;
}

int switch_panel_orientation (double **Z, DATAINFO *pdinfo)
{
    double **tmpZ;
    int i, j, k, t, nvec;
    int nunits = pdinfo->pd;
    int nperiods = pdinfo->n / nunits;
    char **markers = NULL;

    tmpZ = malloc((pdinfo->v - 1) * sizeof *tmpZ);
    if (tmpZ == NULL) return E_ALLOC;

    /* allocate temporary data matrix */
    j = 0;
    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
	    tmpZ[j] = malloc(pdinfo->n * sizeof **tmpZ);
	    if (tmpZ[j] == NULL) {
		for (i=0; i<j; i++) free(tmpZ[i]);
		free(tmpZ);
		return E_ALLOC;
	    }
	    j++;
	} 
    }
    nvec = j;

    /* allocate marker space if relevant */
    if (pdinfo->S != NULL) {
	markers = malloc(pdinfo->n * sizeof *markers);
	if (markers != NULL) {
	    for (t=0; t<pdinfo->n; t++) {
		markers[t] = malloc(OBSLEN);
		if (markers[t] == NULL) {
		    free(markers);
		    markers = NULL;
		    break;
		} else {
		    strcpy(markers[t], pdinfo->S[t]);
		}
	    }
	}
    }

    /* copy the data (vectors) across */
    j = 0;
    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
	    for (t=0; t<pdinfo->n; t++) {
		tmpZ[j][t] = Z[i][t];
	    }
	    j++;
	} 
    }

    /* copy the data back in transformed order: construct a set of
       time series for each unit in turn -- and do markers if present */
    for (k=0; k<nunits; k++) {
	j = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (!pdinfo->vector[i]) continue;
	    for (t=0; t<nperiods; t++) {
		Z[i][k * nperiods + t] = tmpZ[j][k + nunits * t];
	    }
	    j++;
	}
	if (markers != NULL) {
	    for (t=0; t<nperiods; t++) {
		strcpy(pdinfo->S[k * nperiods + t], markers[k + nunits * t]);
	    }
	}
    }

    /* change the datainfo setup */
    pdinfo->time_series = STACKED_TIME_SERIES;
    pdinfo->pd = nperiods;
    if (nperiods < 9) {
	strcpy(pdinfo->stobs, "1:1");
    } else {
	strcpy(pdinfo->stobs, "1:01");
    }
    pdinfo->sd0 = obs_str_to_double(pdinfo->stobs);
    ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);

    /* clean up */
    for (i=0; i<nvec; i++) free(tmpZ[i]);
    free(tmpZ);
    if (markers != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    free(markers[t]);
	}
	free(markers);
    }

    return 0;
}



