/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#define ADF_DEBUG 1
#define KPSS_DEBUG 0

#include "libgretl.h" 
#include "transforms.h"

enum {
    ADF_EG_TEST   = 1 << 0,
    ADF_EG_RESIDS = 1 << 1
} adf_flags;

/* replace y with demeaned or detrended y */

static int GLS_demean_detrend (double *y, int T, int test)
{
    gretl_matrix *ya = NULL;
    gretl_matrix *Za = NULL;
    gretl_matrix *b = NULL;
    double c;
    int t, zcols;
    int err = 0;

    zcols = (test == UR_TREND)? 2 : 1;

    ya = gretl_column_vector_alloc(T);
    Za = gretl_matrix_alloc(T, zcols);
    b = gretl_column_vector_alloc(zcols);

    if (ya == NULL || Za == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    c = (test == UR_CONST)? (1.0 - 7.0/T) : (1.0 - 13.5/T);

    ya->val[0] = y[0];
    for (t=1; t<T; t++) {
	ya->val[t] = y[t] - y[t-1] * c;
    }

    gretl_matrix_set(Za, 0, 0, 1);
    if (zcols == 2) {
	gretl_matrix_set(Za, t, 1, 1);
    }

    for (t=1; t<T; t++) {
	gretl_matrix_set(Za, t, 0, 1 - c);
	if (zcols == 2) {
	    gretl_matrix_set(Za, t, 1, t+1 - t*c);
	}
    }

    err = gretl_matrix_ols(ya, Za, b, NULL, NULL, NULL);

    if (!err) {
	for (t=0; t<T; t++) {
	    y[t] -= b->val[0];
	    if (zcols == 2) {
		y[t] -= b->val[1] * (t+1);
	    }
	}
    }    

 bailout:

    gretl_matrix_free(ya);
    gretl_matrix_free(Za);
    gretl_matrix_free(b);
    
    return err;
}

static int *
adf_prepare_vars (int order, int varno, int nseas, int *d0,
		  double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt)
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

    /* GLS adjustment wanted? */
    if (opt & OPT_G) {
	int t, v = pdinfo->v;

	err = dataset_add_series(1, pZ, pdinfo);
	if (!err) {
	    for (t=0; t<=pdinfo->t2; t++) {
		(*pZ)[v][t] = (*pZ)[varno][t];
	    }
	    err = GLS_demean_detrend((*pZ)[v], pdinfo->t2 + 1, UR_CONST); /* FIXME */
	}
	if (err) {
	    free(list);
	    return NULL;
	}
	strcpy(pdinfo->varname[v], "yd");
	varno = v;
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
		   unsigned char flags, int i, int niv, int nseas, 
		   gretlopt opt, PRN *prn)
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
    const char *urcstrs[] = {
	"nc", "c", "ct", "ctt"
    };

    char pvstr[48];

    if (prn == NULL) return;

    if (!(opt & OPT_G)) {
	i--;
    }

    if (na(pv)) {
	sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
    } else {
	sprintf(pvstr, "%s %.4g", 
		(order > 0)? _("asymptotic p-value") : _("p-value"), 
		pv);
    } 

    if (*blurb_done == 0) {
	if (flags & ADF_EG_RESIDS) {
	    pputc(prn, '\n');
	    pprintf(prn, _("lag order %d\n"), order);
	} else if (flags & ADF_EG_TEST) {
	    if (order > 0) {
		pprintf(prn, _("\nAugmented Dickey-Fuller test, order %d, for %s\n"),
			order, vname);
	    } else {
		pprintf(prn, _("\nDickey-Fuller test for %s\n"), vname);
	    }
	} else {
	    if (order > 0 && !auto_order) {
		pprintf(prn, _("\nAugmented Dickey-Fuller tests, order %d, for %s\n"),
			order, vname);
	    } else {
		pprintf(prn, _("\nDickey-Fuller tests for %s\n"), vname);
	    }
	} 
	pprintf(prn, _("sample size %d\n"), dfmod->nobs);
	pputs(prn, _("unit-root null hypothesis: a = 1"));
	pputs(prn, "\n\n");
	*blurb_done = 1;
    }

    if (!(flags & ADF_EG_RESIDS)) {
	pprintf(prn, "   %s ", _(teststrs[i]));
	if (opt & OPT_G) {
	    pputs(prn, "(GLS) ");
	}
	if (nseas > 0 && i > 0) {
	    pputs(prn, _("plus seasonal dummies"));
	}
	pputc(prn, '\n');
    }

    if (!(flags & ADF_EG_TEST)) {
	pprintf(prn, "   %s: %s\n", _("model"), 
		(order > 0)? aug_models[i] : models[i]);
	if (auto_order) {
	    pprintf(prn, "   %s %d\n", _("lag order:"), order);
	}
	if (!na(dfmod->rho)) {
	    pprintf(prn, "   %s: %.3f\n", _("1st-order autocorrelation coeff. for e"), 
		    dfmod->rho);
	}
    }

    pprintf(prn, "   %s: %g\n"
	    "   %s: tau_%s(%d) = %g\n"
	    "   %s\n",
	    _("estimated value of (a - 1)"), dfmod->coeff[dfnum],
	    _("test statistic"), urcstrs[i], niv, DFt,
	    pvstr);	
}

static int auto_adjust_order (int *list, int order_max,
			      double ***pZ, DATAINFO *pdinfo,
			      PRN *prn)
{
    MODEL kmod;
    double tstat, pval;
    int k, pos;

    for (k=order_max; k>0; k--) {

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

	pos = k + kmod.ifc;
	tstat = kmod.coeff[pos] / kmod.sderr[pos];
	clear_model(&kmod);
	pval = normal_pvalue_2(tstat);

	if (pval > 0.10) {
#if ADF_DEBUG
	    pprintf(prn, "\nauto_adjust_order: lagged difference not "
		    "significant at order %d (t = %g)\n\n", k, tstat);
#endif
	    gretl_list_delete_at_pos(list, k + 2);
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

static int engle_granger_itv (gretlopt opt)
{
    int itv = UR_CONST;

    if (opt & OPT_N) {
	itv = UR_NO_CONST;
    } else if (opt & OPT_T) {
	itv = UR_TREND;
    } else if (opt & OPT_R) {
	itv = UR_QUAD_TREND;
    }

    return itv;
}

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags,
			  PRN *prn)
{
    MODEL dfmod;
    gretlopt eg_opt = OPT_NONE;
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

    if (flags & ADF_EG_RESIDS) {
	/* final step of Engle-Granger test: the (A)DF test regression
	   will contain no deterministic terms, but the selection of the
	   p-value is based on the deterministic terms in the cointegrating
	   regression, represented by "eg_opt".
	 */
	eg_opt = opt;
	opt = OPT_N;
    }

    if (opt & OPT_F) {
	/* difference the variable before testing */
	varno = diffgenr(varno, DIFF, pZ, pdinfo);
	if (varno < 0) {
	    return E_DATA;
	}
    }

    if (opt & OPT_E) {
	auto_order = 1;
    }

    if (order < 0) {
	auto_order = 1;
	order = -order;
    }

    order_max = order;

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: order = %d, auto_order = %d\n", order, auto_order);
#endif

    if ((opt & OPT_D) && pdinfo->pd > 1) {
	nseas = pdinfo->pd - 1;
    }

    list = adf_prepare_vars(order, varno, nseas, &d0, pZ, pdinfo, opt);
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
	/* do the full default set of models */
	opt |= (OPT_C | OPT_T | OPT_R);
    }

    for (i=UR_NO_CONST; i<UR_MAX; i++) {
	int j, k, dfnum = (i > UR_NO_CONST);
	int itv = i;

	if (!test_wanted(i, opt)) {
	    continue;
	}

	if (auto_order) {
	    /* re-establish max order before testing down */
	    order = order_max;
	    copy_list_values(list, biglist);
	}

	if (opt & OPT_G) {
	    fprintf(stderr, "DF-GLS, skipping const, trend\n");
	    list[0] = 2 + order;
	    dfnum--;
	    itv--;
	    goto skipdet;
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

    skipdet:

	if (auto_order) {
	    order = auto_adjust_order(list, order_max, pZ, pdinfo, prn);
	    if (order < 0) {
		err = 1;
		clear_model(&dfmod);
		goto bailout;
	    }
	}

	printlist(list, "final ADF regression list");

	dfmod = lsq(list, pZ, pdinfo, OLS, OPT_A);
	if (dfmod.errcode) {
	    fprintf(stderr, "adf_test: dfmod.errcode = %d\n", 
		    dfmod.errcode);
	    err = dfmod.errcode;
	    clear_model(&dfmod);
	    goto bailout;
	}

	DFt = dfmod.coeff[dfnum] / dfmod.sderr[dfnum];

	if (flags & ADF_EG_RESIDS) {
	    itv = engle_granger_itv(eg_opt);
	} 

	pv = df_pvalue_from_plugin(DFt, 
				   /* use asymptotic p-value for augmented case */
				   (order > 0)? 0 : dfmod.nobs, 
				   niv, itv);

	if (!(opt & OPT_Q)) {
	    print_adf_results(order, auto_order, DFt, pv, &dfmod, dfnum, 
			      pdinfo->varname[varno], &blurb_done, flags,
			      itv, niv, nseas, opt, prn);
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
	if (!(flags & ADF_EG_TEST) || (flags & ADF_EG_RESIDS)) {
	    record_test_result(DFt, pv, "Dickey-Fuller");
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
 * @list: list of variables to test.
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

int adf_test (int order, const int *list, double ***pZ,
	      DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int i, err = 0;

    for (i=1; i<=list[0] && !err; i++) {
	err = real_adf_test(list[i], order, 1, pZ, pdinfo, opt, 
			    0, prn);
    }

    return err;
}

static int 
real_kpss_test (int order, int varno, double ***pZ,
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
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @list: list of variables to test.
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

int kpss_test (int order, const int *list, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int i, err = 0;

    for (i=1; i<=list[0] && !err; i++) {
	err = real_kpss_test(order, list[i], pZ, pdinfo,
			     opt, prn);
    }

    return err;
}

static int *make_coint_list (const int *list, int detcode, int *nv, 
			     double ***pZ, DATAINFO *pdinfo,
			     int *err)
{
    int *clist = NULL;
    int ifc = 0;
    int i, j = 1;

    /* does the incoming list contain a constant? */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    ifc = 1;
	    break;
	}
    }

    /* check for sufficient arguments */
    *nv = list[0] - ifc;
    if (*nv < 2) {
	*err = E_ARGS;
	return NULL;
    }

    /* allocate list for cointegrating regression */
    clist = gretl_list_new(*nv + detcode - 1);
    if (clist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* transcribe original vars */
    for (i=1; i<=list[0]; i++) {
	if (list[i] != 0) {
	    clist[j++] = list[i];
	}
    }    

    /* add trend, if wanted */
    if (detcode >= UR_TREND) {
	clist[j] = gettrend(pZ, pdinfo, 0);
	if (clist[j++] == 0) {
	    *err = E_ALLOC;
	} 
    }

    /* add trend-squared, if wanted */
    if (!*err && detcode == UR_QUAD_TREND) {
	clist[j] = gettrend(pZ, pdinfo, 1);
	if (clist[j++] == 0) {
	    *err = E_ALLOC;
	}
    }

    /* add const, if wanted */
    if (!*err && detcode != UR_NO_CONST) {
	clist[j] = 0;
    } 

    return clist;
}

static int 
coint_check_opts (gretlopt opt, int *detcode, gretlopt *adf_opt)
{
    if (opt & OPT_N) {
	if ((opt & OPT_T) || (opt & OPT_R)) {
	    return E_BADOPT;
	} else {
	    *detcode = UR_NO_CONST;
	    *adf_opt = OPT_N;
	}
    } else if (opt & OPT_T) {
	if (opt & OPT_R) {
	    return E_BADOPT;
	} else {
	    *detcode = UR_TREND;
	    *adf_opt = OPT_T;
	}
    } else if (opt & OPT_R) {
	*detcode = UR_QUAD_TREND;
	*adf_opt = OPT_R;
    }

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
    gretlopt adf_opt = OPT_C;
    MODEL cmod;
    int detcode = UR_CONST;
    int i, nv, k = 0;
    int step = 1;
    int *clist = NULL;
    int err = 0;

    err = coint_check_opts(opt, &detcode, &adf_opt);
    if (err) {
	return err;
    }

    clist = make_coint_list(list, detcode, &nv, pZ, pdinfo, &err);
    if (err) {
	return err;
    }

    gretl_model_init(&cmod);

    if (!(opt & OPT_S)) {
	/* test all candidate vars for unit root */
	for (i=1; i<=nv; i++) {
	    pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		    step++, pdinfo->varname[clist[i]]);
	    real_adf_test(clist[i], order, 1, pZ, pdinfo, adf_opt, 
			  ADF_EG_TEST, prn);
	}
    }

    pprintf(prn, _("Step %d: cointegrating regression\n"), step++);

    cmod = lsq(clist, pZ, pdinfo, OLS, OPT_NONE);
    err = cmod.errcode;
    if (err) {
	goto bailout;
    }

    cmod.aux = AUX_COINT;
    printmodel(&cmod, pdinfo, OPT_NONE, prn);

    /* add residuals from cointegrating regression to data set */
    err = dataset_add_allocated_series(cmod.uhat, pZ, pdinfo);
    if (err) {
	goto bailout;
    }

    k = pdinfo->v - 1;
    strcpy(pdinfo->varname[k], "uhat");
    cmod.uhat = NULL;

    pputc(prn, '\n');
    pprintf(prn, _("Step %d: Dickey-Fuller test on residuals\n"), step);

    /* Run (A)DF test on the residuals */
    real_adf_test(k, order, nv, pZ, pdinfo, adf_opt, 
		  ADF_EG_TEST | ADF_EG_RESIDS, prn); 

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		 "(a) The unit-root hypothesis is not rejected for the individual"
		 " variables.\n(b) The unit-root hypothesis is rejected for the "
		 "residuals (uhat) from the \n    cointegrating regression.\n"));

 bailout:
    
    clear_model(&cmod);
    free(clist);
    if (k > 0) {
	dataset_drop_variable(k, pZ, pdinfo);
    }

    return err;
}
