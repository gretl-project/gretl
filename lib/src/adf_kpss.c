/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#define ADF_DEBUG 0
#define KPSS_DEBUG 0

#include "libgretl.h" 
#include "transforms.h"

enum {
    ADF_EG_TEST   = 1 << 0,
    ADF_PRINT_ACK = 1 << 1
} adf_flags;

static int *
adf_prepare_vars (int order, int varno, int nseas, int *d0,
		  double ***pZ, DATAINFO *pdinfo)
{
    int nl = 5 + nseas;
    int i, j, orig_t1 = pdinfo->t1;
    int *list;
    int err = 0;

    if (varno == 0) {
	return NULL;
    }

    list = gretl_list_new(nl + order);
    if (list == NULL) {
	return NULL;
    }

    /* temporararily reset sample */
    pdinfo->t1 = 0;

    /* generate first difference of the given variable */
    list[1] = diffgenr(varno, DIFF, pZ, pdinfo);
    if (list[1] < 0) {
	pdinfo->t1 = orig_t1;
	free(list);
	return NULL;
    }	

    /* generate lag of given var */
    list[2] = laggenr(varno, 1, pZ, pdinfo); 
    if (list[2] < 0) {
	pdinfo->t1 = orig_t1;
	free(list);
	return NULL;
    }

    /* undo reset sample */
    pdinfo->t1 = orig_t1;

    /* generate lags of difference for augmented test */
    j = 3;
    for (i=1; i<=order && !err; i++) {
	int lnum = laggenr(list[1], i, pZ, pdinfo);

	if (lnum < 0) {
	    fprintf(stderr, "Error generating lag variable\n");
	    err = 1;
	} else {
	    list[j++] = lnum;
	} 
    } 

    if (nseas > 0 && !err) {
	*d0 = dummy(pZ, pdinfo, 0); /* should we center these? */
	if (*d0 < 0) {
	    fprintf(stderr, "Error generating seasonal dummies\n");
	    err = 1;
	} 
    }

    if (err) {
	free(list);
	list = NULL;
    } 

#if ADF_DEBUG
    printlist(list, "adf initial list");
#endif

    return list;
}

static void 
print_adf_results (int order, int auto_order, double DFt, double pv, const MODEL *dfmod,
		   int dfnum, const char *vname, int *blurb_done,
		   unsigned char flags, int i, int nseas, PRN *prn)
{
    const char *models[] = {
	"(1 - L)y = (a-1)*y(-1) + e",
	"(1 - L)y = b0 + (a-1)*y(-1) + e",
	"(1 - L)y = b0 + b1*t + (a-1)*y(-1) + e",
	"(1 - L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + e"
    };
    const char *aug_models[] = {
	"(1 - L)y = (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + b1*t + (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + ... + e"
    };
    const char *teststrs[] = {
	N_("test without constant"),
	N_("test with constant"),
	N_("with constant and trend"),
	N_("with constant and quadratic trend")
    };

    char pvstr[48];

    if (prn == NULL) return;

    if (na(pv)) {
	sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
    } else {
	sprintf(pvstr, "%s %.4g", 
		(order > 0)? _("asymptotic p-value") : _("p-value"), 
		pv);
    } 

    if (*blurb_done == 0) {
	if (order > 0 && !auto_order) {
	    pprintf(prn, _("\nAugmented Dickey-Fuller tests, order %d, for %s\n"),
		    order, vname);
	} else {
	    pprintf(prn, _("\nDickey-Fuller tests for %s\n"), vname);
	}
	pprintf(prn, _("sample size %d\n"), dfmod->nobs);
	pputs(prn, _("unit-root null hypothesis: a = 1"));
	pputs(prn, "\n\n");
	*blurb_done = 1;
    }

    pprintf(prn, "   %s ", _(teststrs[i]));
    if (nseas > 0 && i > 0) {
	pputs(prn, _("plus seasonal dummies"));
    }
    pputc(prn, '\n');

    if (!(flags & ADF_EG_TEST)) {
	pprintf(prn, "   %s: %s\n", _("model"), 
		(order > 0)? aug_models[i] : models[i]);
	if (auto_order) {
	    pprintf(prn, "   %s %d\n", _("lag order:"), order);
	}
    }

    pprintf(prn, "   %s: %g\n"
	    "   %s: t = %g\n"
	    "   %s\n",
	    _("estimated value of (a - 1)"), dfmod->coeff[dfnum],
	    _("test statistic"), DFt,
	    pvstr);	
}

static int auto_adjust_order (int *list, int order_max,
			      double ***pZ, DATAINFO *pdinfo,
			      PRN *prn)
{
    MODEL kmod;
    double tstat, pval = 1.0;
    int l0 = list[0];
    int k;

    for (k=order_max; k>0; k--) {
	int pos = k;

	if (list[l0] == 0) {
	    /* constructed list includes constant */
	    pos++;
	}

	kmod = lsq(list, pZ, pdinfo, OLS, OPT_A);

	if (kmod.errcode) {
	    clear_model(&kmod);
	    fprintf(stderr, "adf: model failed in auto_adjust_order()\n");
	    k = -1;
	    break;
	}

#if ADF_DEBUG
	printmodel(&kmod, pdinfo, OPT_NONE, prn);
#endif

	tstat = kmod.coeff[pos] / kmod.sderr[pos];
	clear_model(&kmod);
	pval = normal_pvalue_2(tstat);

	if (pval > 0.10) {
#if ADF_DEBUG
	    pprintf(prn, "\nauto_adjust_order: lagged difference not "
		    "significant at order %d (t = %g)\n\n", k, tstat);
#endif
	    if (k == 1) {
		k = 0;
		break;
	    } else {
		gretl_list_delete_at_pos(list, k + 2);
	    }
	} else {
#if ADF_DEBUG
	    pprintf(prn, "\nauto_adjust_order: lagged difference is "
		    "significant at order %d (t = %g)\n\n", k, tstat);
#endif
	    break;
	}
    }

    return k;
}

/* targ must be big enough to accept all of src! */

static void copy_list_values (int *targ, const int *src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }
}

static double df_pvalue_from_plugin (double tau, int n, int niv, int itv)
{
    char datapath[FILENAME_MAX];
    void *handle;
    double (*mackinnon_pvalue)(double, int, int, int, char *);
    double pval = NADBL;
    static int nodata;
    
    if (nodata) {
	return pval;
    }

    mackinnon_pvalue = get_plugin_function("mackinnon_pvalue", &handle);
    if (mackinnon_pvalue == NULL) {
	nodata = 1;
        return pval;
    }

    strcpy(datapath, gretl_lib_path());
#ifdef WIN32
    append_dir(datapath, "plugins");
#endif

    pval = (*mackinnon_pvalue)(tau, n, niv, itv, datapath);

#if ADF_DEBUG
    fprintf(stderr, "getting pval: tau=%g, n=%d, niv=%d, itv=%d: pval=%g\n",
	    tau, n, niv, itv, pval);
#endif

    close_plugin(handle);

    if (*datapath == '\0') {
	nodata = 1;
    } 

    return pval;
}

#define test_opt_not_set(o) (!(o & OPT_N) && !(o & OPT_C) && \
                             !(o & OPT_T) && !(o & OPT_R))

static int test_wanted (int test, gretlopt opt)
{
    int ret = 0;

    switch (test) {
    case UR_NO_CONST:
	ret = (opt & OPT_N);
	break;
    case UR_CONST:
	ret = (opt & OPT_C);
	break;
    case UR_TREND:
	ret = (opt & OPT_T);
	break;
    case UR_QUAD_TREND:
	ret = (opt & OPT_R);
	break;
    default:
	break;
    }

    return ret;
}

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags,
			  PRN *prn)
{
    MODEL dfmod;

    int orig_nvars = pdinfo->v;
    int blurb_done = 0;
    int auto_order = 0;
    int order_max = 0;
    int *list;
    int *biglist = NULL;
    double DFt = NADBL;
    double pv = NADBL;
    int i, d0 = 0;
    int nseas = 0;
    int err = 0;

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: got order = %d\n", order);
#endif

    if (opt & OPT_F) {
	/* difference the variable before testing */
	varno = diffgenr(varno, DIFF, pZ, pdinfo);
	if (varno < 0) {
	    return E_DATA;
	}
    }

    if (order < 0) {
	auto_order = 1;
	order = -order;
    }

    order_max = order;

    if ((opt & OPT_D) && pdinfo->pd > 1) {
	nseas = pdinfo->pd - 1;
    }

    list = adf_prepare_vars(order, varno, nseas, &d0, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    if (auto_order) {
	list[0] = order + 5;
	biglist = gretl_list_copy(list);
	if (biglist == NULL) {
	    free(list);
	    return E_ALLOC;
	}
    }

    gretl_model_init(&dfmod);

    if (test_opt_not_set(opt)) {
	/* the default set of models */
	opt |= (OPT_C | OPT_T | OPT_R);
    }

    for (i=UR_NO_CONST; i<UR_MAX; i++) {
	int j, k, dfnum = (i > UR_NO_CONST);

	if (!test_wanted(i, opt)) {
	    continue;
	}

	if (auto_order) {
	    /* re-establish max order before testing down */
	    order = order_max;
	    copy_list_values(list, biglist);
	}

	list[0] = 1 + order + i;

	/* list[1] and list[2], plus the "order" lags, are in common
	   for all models */

	if (i >= UR_TREND) {
	    k = 3 + order;
	    list[k] = gettrend(pZ, pdinfo, 0);
	    if (list[k] == 0) {
		err = E_ALLOC;
		goto bailout;
	    }
	}

	if (i == UR_QUAD_TREND) {
	    k = 4 + order;
	    list[k] = gettrend(pZ, pdinfo, 1);
	    if (list[k] == 0) {
		err = E_ALLOC;
		goto bailout;
	    }
	}

	if (i != UR_NO_CONST) {
	    k = list[0];
	    list[0] += nseas;
	    /* stick constant on end of list */
	    list[list[0]] = 0;
	    /* preceded by seasonal dummies if wanted */
	    for (j=0; j<nseas; j++) {
		list[k++] = d0 + j;
	    }	    
	} 

	if (auto_order) {
	    order = auto_adjust_order(list, order_max, pZ, pdinfo, prn);
	    if (order < 0) {
		err = 1;
		clear_model(&dfmod);
		goto bailout;
	    }
	}

	dfmod = lsq(list, pZ, pdinfo, OLS, OPT_A);
	if (dfmod.errcode) {
	    fprintf(stderr, "adf_test: dfmod.errcode = %d\n", 
		    dfmod.errcode);
	    err = dfmod.errcode;
	    clear_model(&dfmod);
	    goto bailout;
	}

	DFt = dfmod.coeff[dfnum] / dfmod.sderr[dfnum];

	pv = df_pvalue_from_plugin(DFt, 
				   /* use asymptotic p-value for augmented case */
				   (order > 0)? 0 : dfmod.nobs, 
				   niv, i);

	if (!(opt & OPT_Q)) {
	    print_adf_results(order, auto_order, DFt, pv, &dfmod, dfnum, 
			      pdinfo->varname[varno], &blurb_done, flags, 
			      i - 1, nseas, prn);
	}

	if (opt & OPT_V) {
	    /* verbose */
	    dfmod.aux = (order > 0)? AUX_ADF : AUX_DF;
	    if (!na(pv)) {
		gretl_model_set_int(&dfmod, "dfnum", dfnum + 2);
		gretl_model_set_double(&dfmod, "dfpval", pv);
	    }
	    printmodel(&dfmod, pdinfo, OPT_NONE, prn);
	} else if (!(opt & OPT_Q)) {
	    pputc(prn, '\n');
	}

	clear_model(&dfmod);
    }

    if (!err) {
	if (!(flags & ADF_EG_TEST)) {
	    record_test_result(DFt, pv, "Dickey-Fuller");
	}
	if ((flags & ADF_PRINT_ACK) && !(opt & OPT_Q)) {
	    pputs(prn, _("P-values based on MacKinnon (JAE, 1996)\n"));
	}	
    }

 bailout:

    free(list);

    if (biglist != NULL) {
	free(biglist);
    }

    dataset_drop_last_variables(pdinfo->v - orig_nvars, pZ, pdinfo);

    return err;
}

/**
 * adf_test:
 * @order: lag order for the test.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller test for 
 * a unit root.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int adf_test (int order, int varno, double ***pZ,
	      DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    return real_adf_test(varno, order, 1, pZ, pdinfo, opt, 
			 ADF_PRINT_ACK, prn);
}

/**
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the KPSS test for 
 * stationarity.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int kpss_test (int order, int varno, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    MODEL KPSSmod;
    int list[4];
    int hastrend = 0;
    double s2 = 0.0;
    double cumsum = 0.0, cumsum2 = 0.0;
    double teststat;
    double *autocov;
    double et;

    int i, t;
    int t1, t2, T;

    /* sanity check */
    if (order < 0 || varno <= 0 || varno >= pdinfo->v) {
	return 1;
    }

    if (opt & OPT_F) {
	/* difference the variable before testing */
	varno = diffgenr(varno, DIFF, pZ, pdinfo);
	if (varno < 0) {
	    return E_DATA;
	}
    }

    if (opt & OPT_T) {
	hastrend = 1;
    }

    list[0] = (2 + hastrend);
    list[1] = varno;
    list[2] = 0;

    if (hastrend) {
	list[3] = gettrend(pZ, pdinfo, 0);
	if (list[3] == 0) {
	    return E_ALLOC;
	}
    }

    /* OPT_M: reject missing values within sample range */
    KPSSmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_M);
    if (KPSSmod.errcode) {
	clear_model(&KPSSmod);
	return KPSSmod.errcode;
    }

    t1 = KPSSmod.t1;
    t2 = KPSSmod.t2;
    T = KPSSmod.nobs;

    if (opt & OPT_V) {
	KPSSmod.aux = AUX_KPSS;
	printmodel(&KPSSmod, pdinfo, OPT_NONE, prn);
    }
  
    autocov = malloc(order * sizeof *autocov);
    if (autocov == NULL) {
	return E_ALLOC;
    }
  
    for (i=0; i<order; i++) {
	autocov[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	et = KPSSmod.uhat[t];
	if (na(et)) {
	    continue;
	}
	cumsum += et;
	cumsum2 += cumsum * cumsum;
	s2 += et * et;
	for (i=0; i<order; i++) {
	    int s = i + 1;

	    if (t - s >= t1) {
		autocov[i] += et * KPSSmod.uhat[t - s];
	    }
	}
#if KPSS_DEBUG
	fprintf(stderr, "%d: %#12.4g %#12.4g %#12.4g %#12.4g \n", 
		t, et, KPSSmod.uhat[t-1], s2, cumsum2);
#endif
    }

    for (i=0; i<order; i++) {
	double wt = 1.0 - ((double) (i + 1)) / (order + 1);

	s2 += 2.0 * wt * autocov[i];
    }

    s2 /= T;
    teststat = cumsum2 / (s2 * T * T);

    if (opt & OPT_V) {
	pprintf(prn, "  %s: %g\n", _("Robust estimate of variance"), s2);
	pprintf(prn, "  %s: %g\n", _("Sum of squares of cumulated residuals"), 
		cumsum2);
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, _("\nKPSS test for %s %s\n\n"), pdinfo->varname[varno],
		(hastrend)? _("(including trend)") : _("(without trend)"));
	pprintf(prn, _("Lag truncation parameter = %d\n"), order);
	pprintf(prn, "%s = %g\n\n", _("Test statistic"), teststat);
	pprintf(prn, "%*s    ", TRANSLATED_WIDTH(_("Critical values")), " ");
	pprintf(prn, "%g%%      %g%%    %g%%      %g%%\n", 10.0, 5.0, 2.5, 1.0);
	if (hastrend) {
	    pprintf(prn, "%s: %.3f   %.3f   %.3f   %.3f\n\n", 
		    _("Critical values"), 0.119, 0.146, 0.176, 0.216);
	} else {
	    pprintf(prn, "%s: %.3f   %.3f   %.3f   %.3f\n\n", 
		    _("Critical values"), 0.347, 0.463, 0.574, 0.739);
	}
    }

    record_test_result(teststat, NADBL, "KPSS");
    clear_model(&KPSSmod);

    free(autocov);

    return 0;
}

/**
 * coint:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: if OPT_N, do not an include a constant in the
 *       cointegrating regression.
 * @prn: gretl printing struct.
 *
 * Carries out the Engle-Granger test for cointegration.  
 *
 * Returns: 0 on successful completion.
 */

int coint (int order, const int *list, double ***pZ, 
	   DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int i, t, n, nv, l0 = list[0];
    int cpos, ifc, cnum = 0;
    MODEL cmod;
    int *cointlist = NULL;

    cpos = gretl_list_const_pos(list, 2, (const double **) *pZ, pdinfo);
    if (cpos > 0) {
	cnum = list[cpos];
	ifc = 1;
    } else {
	ifc = 0;
    }

    if (order <= 0 || list[0] - ifc < 2) {
	strcpy(gretl_errmsg, "coint: needs a positive lag order "
	       "and at least two variables");
	return 1;
    }

    gretl_model_init(&cmod);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	if (list[i] == cnum) {
	    continue;
	}
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		i, pdinfo->varname[list[i]]);
	real_adf_test(list[i], order, 1, pZ, pdinfo, OPT_NONE, 
		      ADF_EG_TEST, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (!ifc && !(opt & OPT_N)) {
	/* add const to coint regression list */
	cointlist = malloc((l0 + 2) * sizeof *cointlist);
	if (cointlist == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<=l0; i++) {
	    cointlist[i] = list[i];
	}
	cointlist[l0 + 1] = 0;
	cointlist[0] += 1;
    } else {
	cointlist = gretl_list_copy(list);
	if (cointlist == NULL) {
	    return E_ALLOC;
	}
    }

    pprintf(prn, _("Step %d: cointegrating regression\n"), l0 + 1);
    
    cmod = lsq(cointlist, pZ, pdinfo, OLS, OPT_NONE); 
    cmod.aux = AUX_COINT;
    printmodel(&cmod, pdinfo, OPT_NONE, prn);

    /* add residuals from cointegrating regression to data set */
    n = pdinfo->n;
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }
    nv = pdinfo->v - 1;

    for (t=0; t<cmod.t1; t++) {
	(*pZ)[nv][t] = NADBL;
    }
    for (t=cmod.t1; t<=cmod.t2; t++) {
	(*pZ)[nv][t] = cmod.uhat[t];
    }
    for (t=cmod.t2+1; t<n; t++) {
	(*pZ)[nv][t] = NADBL;
    }

    strcpy(pdinfo->varname[nv], "uhat");

    pputc(prn, '\n');
    pprintf(prn, _("Step %d: Dickey-Fuller test on residuals\n"), l0 + 2);

    /* Run (A)DF test on the residuals */
    real_adf_test(pdinfo->v - 1, order, 1 + cmod.ncoeff - cmod.ifc, 
		  pZ, pdinfo, OPT_N, ADF_EG_TEST | ADF_PRINT_ACK, prn);

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		 "(a) The unit-root hypothesis is not rejected for the individual"
		 " variables.\n(b) The unit-root hypothesis is rejected for the "
		 "residuals (uhat) from the \n    cointegrating regression.\n"));

    /* clean up and get out */
    clear_model(&cmod);
    free(cointlist);
    dataset_drop_last_variables(1, pZ, pdinfo);

    return 0;
}
