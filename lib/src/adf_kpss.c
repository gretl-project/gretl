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

#define ADF_DEBUG 0
#define KPSS_DEBUG 0

#include "libgretl.h" 
#include "transforms.h"

/**
 * SECTION:adf_kpss
 * @short_description: unit root and cointegration tests
 * @title: ADF, KPSS, Engle-Granger 
 * @include: libgretl.h
 *
 * Implementations of the (Augmented) Dickey-Fuller test
 * and the Kwiatkowski, Phillips, Schmidt and Shin test for
 * the presence of a unit root in a time series, along with
 * the Engle-Granger test for cointegration of two or more
 * time series.
 *
 * The Johansen cointegration test is also provided in
 * libgretl; see johansen_test() and johansen_test_simple().
 */

typedef struct adf_info_ adf_info;
typedef struct kpss_info_ kpss_info;

struct adf_info_ {
    int T;
    int df;
    int order;
    double b;
    double tau;
    double pval;
};

struct kpss_info_ {
    int T;
    double test;
    double pval;
};

enum {
    UR_NO_CONST = 1,
    UR_CONST,
    UR_TREND,
    UR_QUAD_TREND,
    UR_MAX
} AdfCode;

enum {
    ADF_EG_TEST   = 1 << 0,
    ADF_EG_RESIDS = 1 << 1,
    ADF_PANEL     = 1 << 2
} adf_flags;

/* replace y with demeaned or detrended y */

static int GLS_demean_detrend (double *y, int T, int test)
{
    gretl_matrix *ya = NULL;
    gretl_matrix *Za = NULL;
    gretl_matrix *b = NULL;
    double c, b0, b1 = 0;
    int t, zcols;
    int err = 0;

    zcols = (test == UR_TREND)? 2 : 1;

    if (T - zcols <= 0) {
	return E_DF;
    }

    ya = gretl_column_vector_alloc(T);
    Za = gretl_matrix_alloc(T, zcols);
    b = gretl_column_vector_alloc(zcols);

    if (ya == NULL || Za == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    c = (test == UR_CONST)? (1.0 - 7.0/T) : (1.0 - 13.5/T);

    gretl_vector_set(ya, 0, y[0] /* (1 - c) * y[0] ?? */);
    for (t=1; t<T; t++) {
	gretl_vector_set(ya, t, y[t] - c * y[t-1]);
    }

    gretl_matrix_set(Za, 0, 0, 1);
    if (zcols == 2) {
	gretl_matrix_set(Za, 0, 1, 1);
    }

    for (t=1; t<T; t++) {
	gretl_matrix_set(Za, t, 0, 1 - c);
	if (zcols == 2) {
	    gretl_matrix_set(Za, t, 1, t+1 - t*c);
	}
    }

    err = gretl_matrix_ols(ya, Za, b, NULL, NULL, NULL);

    if (!err) {
	b0 = gretl_vector_get(b, 0);
	if (zcols == 2) {
	    b1 = gretl_vector_get(b, 1);
	}
	for (t=0; t<T; t++) {
	    y[t] -= b0;
	    if (zcols == 2) {
		y[t] -= b1 * (t+1);
	    }
	}
    }    

 bailout:

    gretl_matrix_free(ya);
    gretl_matrix_free(Za);
    gretl_matrix_free(b);
    
    return err;
}

/* generate the various differences and lags required for
   the ADF test */

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

    if (opt & OPT_G) {
	/* GLS adjustment wanted */
	int test = (opt & OPT_T)? UR_TREND : UR_CONST;
	int t, v = pdinfo->v;

	err = dataset_add_series(1, pZ, pdinfo);
	if (!err) {
	    int T, offset = 0;

	    for (t=0; t<=pdinfo->t2; t++) {
		(*pZ)[v][t] = (*pZ)[varno][t];
		if (na((*pZ)[v][t])) {
		    offset = t+1;
		}
	    }
	    if (offset < pdinfo->t1 - order) {
		offset = pdinfo->t1 - order;
	    }
	    T = pdinfo->t2 - offset + 1;
	    err = GLS_demean_detrend((*pZ)[v] + offset, T, test);
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

static double *df_gls_ct_cval (int T)
{
    /* Elliott, Rothenberg and Stock (1996), Table 1 */
    static double df_gls_ct_cvals[4][4] = {
	/* 10%     5%    2.5%    1% */
	{ -2.89, -3.19, -3.46, -3.77 }, /* T = 50  */
	{ -2.74, -3.03, -3.29, -3.58 }, /* T = 100 */
	{ -2.64, -2.93, -3.18, -3.46 }, /* T = 200 */
	{ -2.57, -2.89, -3.15, -3.48 }  /* \infty  */
    };
    int i = 3;

    if (T <= 50){
	i = 0;
    } else if (T <= 100){
	i = 1;
    } else if (T <= 200){
	i = 2;
    } 

    return df_gls_ct_cvals[i];
}

/* display an F-test for the joint significance of the lagged
   \delta y terms in ADF test */

static void show_lags_test (MODEL *pmod, int order, PRN *prn)
{
    int *llist = gretl_list_new(order);
    double F;
    int i;

    if (llist != NULL) {
	for (i=0; i<order; i++) {
	    /* lagged differences */
	    llist[i+1] = pmod->list[pmod->ifc + 3 + i];
	}

	F = wald_omit_F(llist, pmod);

	if (!na(F)) {
	    pprintf(prn, "   %s: F(%d, %d) = %.3f [%.4f]\n",
		    _("lagged differences"), order, pmod->dfd, F, 
		    snedecor_cdf_comp(order, pmod->dfd, F));
	}

	free(llist);
    }
}

static void DF_header (const char *s, int p, int pmax,
		       gretlopt opt, PRN *prn)
{
    pputc(prn, '\n');

    if (p <= 0) {
	if (opt & OPT_G) {
	    pprintf(prn, _("Dickey-Fuller (GLS) test for %s\n"), s);
	} else {
	    pprintf(prn, _("Dickey-Fuller test for %s\n"), s);
	}	
    } else {
	if (opt & OPT_G) {
	    pprintf(prn, _("Augmented Dickey-Fuller (GLS) test for %s\n"), s);
	} else {
	    pprintf(prn, _("Augmented Dickey-Fuller test for %s\n"), s);
	}
	if (p == 1) {
	    pprintf(prn, _("including one lag of (1-L)%s"), s);
	} else {
	    pprintf(prn, _("including %d lags of (1-L)%s"), p, s);
	}
	if (pmax >= p) {
	    pputc(prn, ' ');
	    pprintf(prn, _("(max was %d)"), pmax);
	}
	pputc(prn, '\n');
    }
}

static const char *DF_model_string (int i)
{
    const char *models[] = {
	"(1-L)y = (a-1)*y(-1) + e",
	"(1-L)y = b0 + (a-1)*y(-1) + e",
	"(1-L)y = b0 + b1*t + (a-1)*y(-1) + e",
	"(1-L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + e"
    };

    if (i >= 0 && i < 4) {
	return models[i];
    } else {
	return "";
    }
}

static const char *ADF_model_string (int i)
{
    const char *models[] = {
	"(1-L)y = (a-1)*y(-1) + ... + e",
	"(1-L)y = b0 + (a-1)*y(-1) + ... + e",
	"(1-L)y = b0 + b1*t + (a-1)*y(-1) + ... + e",
	"(1-L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + ... + e"
    };

    if (i >= 0 && i < 4) {
	return models[i];
    } else {
	return "";
    }
}

static const char *DF_test_string (int i)
{
    const char *tests[] = {
	N_("test without constant"),
	N_("test with constant"),
	N_("with constant and trend"),
	N_("with constant and quadratic trend")
    };

    if (i >= 0 && i < 4) {
	return tests[i];
    } else {
	return "";
    }
}

static void 
print_adf_results (int order, int pmax, double DFt, double pv, 
		   MODEL *dfmod, int dfnum, const char *vname, 
		   int *blurb_done, unsigned char flags, int i, 
		   int niv, int nseas, gretlopt opt, PRN *prn)
{
    const char *urcstrs[] = {
	"nc", "c", "ct", "ctt"
    };
    char pvstr[48];
    char taustr[16];

    if (prn == NULL) return;

    /* convert test-type to 0-base */
    i--;

    if (na(pv)) {
	sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
    } else {
	int asy = (order > 0 || (opt & OPT_G));

	sprintf(pvstr, "%s %.4g", 
		(asy)? _("asymptotic p-value") : _("p-value"), 
		pv);
    } 

    if (*blurb_done == 0) {
	DF_header(vname, order, pmax, opt, prn);
	pprintf(prn, _("sample size %d\n"), dfmod->nobs);
	if (flags & ADF_PANEL) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, _("unit-root null hypothesis: a = 1"));
	    pputs(prn, "\n\n");
	}
	*blurb_done = 1;
    }

    if (!(flags & ADF_EG_RESIDS)) {
	pprintf(prn, "   %s ", _(DF_test_string(i)));
	if (nseas > 0 && i > 0) {
	    pputs(prn, _("plus seasonal dummies"));
	}
	pputc(prn, '\n');
    }

    pprintf(prn, "   %s: %s\n", _("model"), 
	    (order > 0)? ADF_model_string(i) : DF_model_string(i));

    if (!na(dfmod->rho)) {
	pprintf(prn, "   %s: %.3f\n", _("1st-order autocorrelation coeff. for e"), 
		dfmod->rho);
    }

    if (order > 1) {
	show_lags_test(dfmod, order, prn);
    }	

    if (opt & OPT_G) {
	strcpy(taustr, "tau");
    } else {
	sprintf(taustr, "tau_%s(%d)", urcstrs[i], niv);
    }

    pprintf(prn, "   %s: %g\n"
	    "   %s: %s = %g\n",
	    _("estimated value of (a - 1)"), dfmod->coeff[dfnum],
	    _("test statistic"), taustr, DFt);

    if ((opt & OPT_G) && i+1 == UR_TREND) {
	const double *c = df_gls_ct_cval(dfmod->nobs);

	pprintf(prn, "\n   %*s    ", TRANSLATED_WIDTH(_("Critical values")), " ");
	pprintf(prn, "%g%%     %g%%     %g%%     %g%%\n", 10.0, 5.0, 2.5, 1.0);
	pprintf(prn, "   %s: %.2f   %.2f   %.2f   %.2f\n", 
		_("Critical values"), c[0], c[1], c[2], c[3]);
    } else {
	pprintf(prn, "   %s\n", pvstr);
    } 
}

static int auto_adjust_order (int *list, int order_max,
			      double **Z, DATAINFO *pdinfo,
			      PRN *prn)
{
    MODEL kmod;
    double tstat, pval;
    int k, pos;

    for (k=order_max; k>0; k--) {
	kmod = lsq(list, Z, pdinfo, OLS, OPT_A);
	if (kmod.errcode) {
	    fprintf(stderr, "auto_adjust_order: k = %d, err = %d\n", k,
		    kmod.errcode);
	    clear_model(&kmod);
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

/**
 * get_urc_pvalue:
 * @tau: test statistic.
 * @n: sample size (or 0 for asymptotic result).
 * @niv: number of potentially cointegrated variables
 * (1 for simple unit-root test).
 * @itv: code: 1, 2, 3, 4 for nc, c, ct, ctt models
 * respectively.
 * @opt: give OPT_G if GLS adjustment was applied in
 * the test from which @tau was obtained.
 *
 * Retrieves the p-value for @tau from the Dickey–Fuller 
 * unit-root test or the Engle–Granger cointegration 
 * test, as per James MacKinnon (1996).
 *
 * Returns: p-value, or %NADBL on failure.
 */

double get_urc_pvalue (double tau, int n, int niv, int itv,
		       gretlopt opt)
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

    if ((opt & OPT_G) && itv == UR_CONST) {
	itv = UR_NO_CONST;
    }

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

static int gettrend (double ***pZ, DATAINFO *pdinfo, int square)
{
    int idx, t, v = pdinfo->v;
    double x;

    idx = series_index(pdinfo, (square)? "timesq" : "time");

    if (idx < v) {
	return idx;
    }

    if (dataset_add_series(1, pZ, pdinfo)) {
	return 0; /* error: valid value cannot == 0 */
    }

    for (t=0; t<pdinfo->n; t++) {
	x = (double) t + 1; 
	(*pZ)[v][t] = (square)? x * x : x;
    }

    if (square) {
	strcpy(pdinfo->varname[v], "timesq");
	strcpy(VARLABEL(pdinfo, v), _("squared time trend variable"));
    } else {
	strcpy(pdinfo->varname[v], "time");
	strcpy(VARLABEL(pdinfo, v), _("time trend variable"));
    }
	    
    return idx;
}

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags,
			  adf_info *ainfo, PRN *prn)
{
    MODEL dfmod;
    gretlopt eg_opt = OPT_NONE;
    gretlopt df_mod_opt = OPT_A | OPT_Z;
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

    if (gretl_isconst(pdinfo->t1, pdinfo->t2, (*pZ)[varno])) {
	gretl_errmsg_sprintf(_("%s is a constant"), pdinfo->varname[varno]);
	return E_DATA;
    }    

    if (opt & OPT_E) {
	/* testing down */
	auto_order = 1;
	order_max = order;
    }

    if (order < 0) {
	/* testing down: backward compatibility */
	auto_order = 1;
	order = order_max = -order;
    }

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: order = %d, auto_order = %d\n", order, auto_order);
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

#if 0
    if (flags & ADF_PANEL) {
	/* testing for unit root with panel data: compute standard
	   errors (and so the ADF t-statistic) without a degrees-of-
	   freedom adjustment in the ADF regression
	*/
	df_mod_opt |= OPT_N;
    }
#endif

    if (opt & OPT_F) {
	/* difference the variable before testing */
	int t1 = pdinfo->t1;

	pdinfo->t1 = 0;
	varno = diffgenr(varno, DIFF, pZ, pdinfo);
	pdinfo->t1 = t1;
	if (varno < 0) {
	    return E_DATA;
	}
    }

    if ((opt & OPT_D) && pdinfo->pd > 1) {
	/* add seasonal dummies */
	nseas = pdinfo->pd - 1;
    }

    if (test_opt_not_set(opt)) {
	/* default model(s) */
	if (opt & OPT_G) {
	    opt |= OPT_C;
	} else {
	    opt |= (OPT_C | OPT_T);
	}
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
	    /* DF-GLS: skip const, trend */
	    list[0] = 2 + order;
	    dfnum--;
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
	    order = auto_adjust_order(list, order_max, *pZ, pdinfo, prn);
	    if (order < 0) {
		err = 1;
		clear_model(&dfmod);
		goto bailout;
	    }
	}

#if ADF_DEBUG
	printlist(list, "final ADF regression list");
#endif

	dfmod = lsq(list, *pZ, pdinfo, OLS, df_mod_opt);
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

	if (getenv("DFGLS_NO_PVALUE")) {
	    /* to speed up monte carlo stuff */
	    pv = NADBL;
	} else if ((opt & OPT_G) && itv == UR_TREND) {
	    /* DF-GLS with trend: MacKinnon p-values won't work */
	    pv = NADBL;
	} else {
	    /* Use asymp. p-value in augmented case; also the
	       finite-sample MacKinnon p-values are not correct
	       in the GLS case.
	    */
	    int asymp = (order > 0 || (opt & OPT_G));

	    pv = get_urc_pvalue(DFt, (asymp)? 0 : dfmod.nobs, 
				niv, itv, opt);
	}

	if (ainfo != NULL) {
	    ainfo->T = dfmod.nobs;
	    ainfo->df = dfmod.dfd;
	    ainfo->order = order;
	    ainfo->b = dfmod.coeff[dfnum];
	    ainfo->tau = DFt;
	    ainfo->pval = pv;
	}

	if (!(opt & OPT_Q) && !(flags & ADF_PANEL)) {
	    print_adf_results(order, order_max, DFt, pv, &dfmod, dfnum, 
			      pdinfo->varname[varno], &blurb_done, flags,
			      itv, niv, nseas, opt, prn);
	}

	if ((opt & OPT_V) && !(flags & ADF_PANEL)) {
	    /* verbose */
	    dfmod.aux = (order > 0)? AUX_ADF : AUX_DF;
	    if (!na(pv)) {
		gretl_model_set_int(&dfmod, "dfnum", dfnum);
		gretl_model_set_double(&dfmod, "dfpval", pv);
	    }
	    printmodel(&dfmod, pdinfo, OPT_NONE, prn);
	} else if (!(opt & OPT_Q) && !(flags & (ADF_EG_RESIDS | ADF_PANEL))) {
	    pputc(prn, '\n');
	}

	clear_model(&dfmod);
    }

    if (!err) {
	if (!(flags & (ADF_EG_TEST | ADF_PANEL)) || (flags & ADF_EG_RESIDS)) {
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

/* print critical value for ADF-IPS or KPSS */

static void print_critical_values (double *a, double *cv, 
				   int ci, PRN *prn)
{
    const char *label = N_("Critical values");
    int figs = (ci == ADF)? 2 : 3;

    pprintf(prn, "%*s    ", TRANSLATED_WIDTH(_(label)), " ");
    pprintf(prn, "%g%%      %g%%      %g%%\n", 
	    100*a[0], 100*a[1], 100*a[2]);
    pprintf(prn, "%s: %.*f   %.*f   %.*f\n", 
	    _(label), figs, cv[0], figs, cv[1], figs, cv[2]);
}

static int panel_adjust_ADF_opt (gretlopt *opt)
{
    int err = 0;

    /* Has the user selected an option governing the 
       deterministic terms to be included? If so, don't 
       mess with it.
    */
    if (*opt & (OPT_N | OPT_C | OPT_R | OPT_T)) {
	; /* no-op */
    } else {
	/* panel default: test with constant */
	*opt |= OPT_C;
    }

    return err;
}

static int DF_index (gretlopt opt)
{
    if (opt & OPT_N) {
	return 0;
    } else if (opt & OPT_C) {
	return 1;
    } else if (opt & OPT_T) {
	return 2;
    } else {
	return 3;
    }
}

/* See Im, Pesaran and Shin, "Testing for unit roots in
   heterogeneous panels", Journal of Econometrics 115 (2003),
   53-74.
*/

static int do_IPS_test (double tbar, int n, const int *Ti, 
			int order, const int *Oi,
			const DATAINFO *pdinfo,
			gretlopt opt, PRN *prn)
{
    int (*get_IPS_critvals) (int, int, int, double *);
    int (*IPS_tbar_moments) (int, double *, double *);
    int (*IPS_tbar_rho_moments) (int, int, int, double *, double *);
    void *handle = NULL;
    int Tmin = Ti[1], Tmax = Ti[1];
    int i, T, err = 0;

    for (i=2; i<=n; i++) {
	if (Ti[i] > Tmax) {
	    Tmax = Ti[i];
	}
	if (Ti[i] < Tmin) {
	    Tmin = Ti[i];
	}
    }

    if (Oi != NULL || order > 0) {
	/* non-zero lag order: use IPS's W_{tbar} statistic */
	double E, V, Wtbar, Etbar = 0, Vtbar = 0;
	int order_i;

	IPS_tbar_rho_moments = get_plugin_function("IPS_tbar_rho_moments", &handle);

	if (IPS_tbar_rho_moments != NULL) {
	    for (i=0; i<n && !err; i++) {
		T = Ti[i+1];
		order_i = (Oi != NULL)? Oi[i+1] : order;
		err = IPS_tbar_rho_moments(order_i, T, (opt & OPT_T), &E, &V);
		Etbar += E;
		Vtbar += V;
	    }

	    if (!err) {
		Etbar /= n;
		Vtbar /= n;
		Wtbar = sqrt(n) * (tbar - Etbar) / sqrt(Vtbar);
		pprintf(prn, "N = %d, Tmin = %d, Tmax = %d\n", n, Tmin, Tmax);
		pprintf(prn, "Im-Pesaran-Shin W_tbar = %g [%.4f]\n", Wtbar, 
			normal_pvalue_1(-Wtbar));
	    }
	}
    } else if (Tmax > Tmin) {
	/* sample sizes differ: use IPS's Z_{tbar} */
	double E, V, Ztbar, Etbar = 0, Vtbar = 0;

	IPS_tbar_moments = get_plugin_function("IPS_tbar_moments", &handle);

	if (IPS_tbar_moments != NULL) {
	    for (i=0; i<n && !err; i++) {
		T = Ti[i+1];
		err = IPS_tbar_moments(T, &E, &V);
		Etbar += E;
		Vtbar += V;
	    }

	    if (!err) {
		Etbar /= n;
		Vtbar /= n;
		Ztbar = sqrt(n) * (tbar - Etbar) / sqrt(Vtbar);
		pprintf(prn, "N = %d, Tmin = %d, Tmax = %d\n", n, Tmin, Tmax);
		pprintf(prn, "Im-Pesaran-Shin Z_tbar = %g [%.4f]\n", Ztbar, 
			normal_pvalue_1(-Ztbar));
	    }
	}
    } else {
	/* simple case: use tbar with exact critical values */
	pprintf(prn, "N,T = (%d,%d)\n", n, Tmax);
	pprintf(prn, "Im-Pesaran-Shin t-bar = %g\n", tbar);

	get_IPS_critvals = get_plugin_function("get_IPS_critvals", &handle);

	if (get_IPS_critvals != NULL) {
	    double a[] = { 0.1, 0.05, 0.01 };
	    double cv[3];
		
	    err = (*get_IPS_critvals) (n, Tmax, (opt & OPT_T), cv);
	    if (!err) {
		print_critical_values(a, cv, ADF, prn);
	    }
	}
    }

    if (handle != NULL) {
	close_plugin(handle);
    }

    return err;
}

/* See In Choi, "Unit root tests for panel data", Journal of
   International Money and Finance 20 (2001), 249-272.
*/

static void do_choi_test (double ppv, double zpv, double lpv, 
			  int n, PRN *prn)
{
    double P = -2 * ppv;
    double Z = zpv / sqrt((double) n);
    int tdf = 5 * n + 4;
    double k = (3.0*tdf)/(M_PI*M_PI*n*(5*n+2));
    double L = sqrt(k) * lpv;

    pprintf(prn, "%s\n", _("Choi meta-tests:"));
    pprintf(prn, "   %s(%d) = %g [%.4f]\n", _("Inverse chi-square"),
	    2*n, P, chisq_cdf_comp(2*n, P));
    pprintf(prn, "   %s = %g [%.4f]\n", _("Inverse normal test"),
	    Z, normal_pvalue_1(-Z));
    pprintf(prn, "   %s: t(%d) = %g [%.4f]\n", _("Logit test"),
	    tdf, L, student_pvalue_1(tdf, -L));
}

static int panel_DF_test (int v, int order, 
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, PRN *prn)
{
    int u0 = pdinfo->t1 / pdinfo->pd;
    int uN = pdinfo->t2 / pdinfo->pd;
    int quiet = (opt & OPT_Q);
    int verbose = (opt & OPT_V);
    double ppv = 0.0, zpv = 0.0, lpv = 0.0;
    double pval, tbar = 0.0;
    int *Ti = NULL, *Oi = NULL;
    int i, n, err;

    err = panel_adjust_ADF_opt(&opt);
    if (err) {
	return err;
    }

    if (opt & OPT_G) {
	/* GLS option: can't do IPS t-bar test */
	tbar = NADBL;
    } else {
	Ti = gretl_list_new(uN - u0 + 1);
	if (Ti == NULL) {
	    return E_ALLOC;
	}
	if ((opt & OPT_E) || order < 0) {
	    /* testing down: lag order may vary by unit */
	    Oi = gretl_list_new(uN - u0 + 1);
	    if (Oi == NULL) {
		free(Ti);
		return E_ALLOC;
	    }
	}	    
    }

    if (!quiet) {
	int j = DF_index(opt);

	DF_header(pdinfo->varname[v], (opt & OPT_E)? 0 : order, 0, opt, prn);
	pprintf(prn, "   %s ", _(DF_test_string(j)));
	pputc(prn, '\n');
	pprintf(prn, "   %s: %s\n\n", _("model"), 
		(order > 0)? ADF_model_string(j) : DF_model_string(j));
    }

    /* number of units in sample range */
    n = uN - u0 + 1;

    /* run a Dickey-Fuller test for each unit and record the
       results */

    for (i=u0; i<=uN && !err; i++) {
	adf_info ainfo;

	pdinfo->t1 = i * pdinfo->pd;
	pdinfo->t2 = pdinfo->t1 + pdinfo->pd - 1;
	err = series_adjust_sample((*pZ)[v], &pdinfo->t1, &pdinfo->t2);

	if (!err) {
	    err = real_adf_test(v, order, 1, pZ, pdinfo, opt, 
				ADF_PANEL, &ainfo, prn);
	    if (!err && verbose) {
		pprintf(prn, "%s %d, T = %d, %s = %d\n", _("Unit"), i + 1, 
			ainfo.T, _("lag order"), ainfo.order);
		pprintf(prn, "   %s: %g\n"
			"   %s = %g", 
			_("estimated value of (a - 1)"), ainfo.b,
			_("test statistic"), ainfo.tau);
		if (na(ainfo.pval)) {
		    pputs(prn, "\n\n");
		} else {
		    pprintf(prn, " [%.4f]\n\n", ainfo.pval);
		}
	    }	    
	}

	if (!err) {
	    if (Ti != NULL) {
		Ti[i-u0+1] = ainfo.T;
	    }
	    if (Oi != NULL) {
		Oi[i-u0+1] = ainfo.order;
	    }	    
	    if (!na(tbar)) {
		tbar += ainfo.tau;
	    }
	    pval = ainfo.pval;
	    if (na(pval)) {
		ppv = zpv = lpv = NADBL;
	    } else if (!na(ppv)) {
		ppv += log(pval);
		zpv += normal_cdf_inverse(pval);
		lpv += log(pval / (1-pval));
	    }
	}
    }

    /* process the results as per Im-Pesaran-Shin and/or Choi */

    if (!err) {
	pprintf(prn, "%s\n\n", _("H0: all groups have unit root"));
	if (!na(tbar)) {
	    tbar /= n;
	    do_IPS_test(tbar, n, Ti, order, Oi, pdinfo, opt, prn);
	}
	if (!na(ppv)) {
	    pputc(prn, '\n');
	    do_choi_test(ppv, zpv, lpv, n, prn);
	}
	pputc(prn, '\n');
    }

    free(Ti);
    free(Oi);

    return err;
}

/**
 * levin_lin_test:
 * @vnum: ID number of variable to test.
 * @plist: list of ADF lag orders.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Levin-Lin-Chu test
 * for a unit root in panel data. 
 *
 * The list @plist should contain either a single lag order
 * to be applied to all units, or a set of unit-specific
 * orders; in the latter case the length of the list must
 * equal the number of panel units in the current sample
 * range. (This is a gretl list: the first element holds
 * a count of the number of elements following.)
 *
 * By default a test with constant is performed, but the
 * (mutually exclusive) options OPT_N and OPT_T in @opt switch to
 * the case of no constant or constant plus trend respectively.
 * The OPT_Q flag may be used to suppress printed output.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int levin_lin_test (int vnum, const int *plist,
		    double **Z, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn)
{
    int (*real_levin_lin) (int, const int *, double **,
			   DATAINFO *, gretlopt, PRN *);
    void *handle;
    int panelmode;
    int err = 0;

    panelmode = multi_unit_panel_sample(pdinfo);

    if (!panelmode || incompatible_options(opt, OPT_N | OPT_T)) {
	return E_BADOPT;
    }

    real_levin_lin = get_plugin_function("real_levin_lin", &handle);

    if (real_levin_lin == NULL) {
	fputs(I_("Couldn't load plugin function\n"), stderr);
	err = E_FOPEN;
    } else {
	int save_t1 = pdinfo->t1;
	int save_t2 = pdinfo->t2;
 
	err = (*real_levin_lin) (vnum, plist, Z, pdinfo, opt, prn);
	close_plugin(handle);

	pdinfo->t1 = save_t1;
	pdinfo->t2 = save_t2;
    }

    return err;
}

/**
 * adf_test:
 * @order: lag order for the (augmented) test.
 * @list: list of variables to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller 
 * test for a unit root. 
 *
 * By default two tests are performed, one for a model
 * including a constant and one including a linear trend. The 
 * deterministic components of the model can be controlled via
 * flags in @opt as follows: OPT_N, omit the constant; OPT_C,
 * run just one test using the constant; OPT_T, one test including
 * linear trend; OPT_R, one test including a quadratic trend;
 * OPT_D, include seasonal dummy variables.
 *
 * Additional flags that may be given in @opt include: 
 * OPT_V for verbose operation; OPT_F to apply first-differencing
 * before testing; OPT_G for GLS preprocessing as in Elliott, Rothenberg
 * and Stock (incompatible with OPT_N, OPT_R, OPT_D); OPT_E to
 * "test down" from a given maximum lag order (see the entry for
 * "adf" in the Gretl Command Reference for details).
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int adf_test (int order, const int *list, double ***pZ,
	      DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int panelmode;
    int err;

    /* GLS incompatible with no const, quadratic trend or seasonals */
    err = incompatible_options(opt, OPT_N | OPT_R | OPT_G);
    if (!err) {
	err = incompatible_options(opt, OPT_D | OPT_G);
    }

    if (!err && (opt & OPT_G)) {
	/* under GLS, have to choose between cases */
	err = incompatible_options(opt, OPT_C | OPT_T);
    }

    panelmode = multi_unit_panel_sample(pdinfo);

    if (panelmode) {
	err = panel_DF_test(list[1], order, pZ, pdinfo, opt, prn);
    } else {
	/* regular time series case */
	int i, v, vlist[2] = {1, 0};

	for (i=1; i<=list[0] && !err; i++) {
	    vlist[1] = v = list[i];
	    err = list_adjust_sample(vlist, &pdinfo->t1, &pdinfo->t2, 
				     (const double **) *pZ);
	    if (!err && order == -1) {
		/* default to L_{12}: see G. W. Schwert, "Tests for Unit Roots:
		   A Monte Carlo Investigation", Journal of Business and
		   Economic Statistics, 7(2), 1989, pp. 5-17.
		*/
		int T = pdinfo->t2 - pdinfo->t1 + 1;

		order = 12.0 * pow(T/100.0, 0.25);
	    }
	    if (!err) {
		err = real_adf_test(v, order, 1, pZ, pdinfo, opt, 
				    0, NULL, prn);
	    }
	    pdinfo->t1 = save_t1;
	    pdinfo->t2 = save_t2;
	}
    }

    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

    return err;
}

/* See Peter S. Sephton, "Response surface estimates of the KPSS 
   stationarity test", Economics Letters 47 (1995) 255-261.

   The estimates below of \beta_\infty and \beta_1 allow the
   construction of better critical values for finite samples
   than the values given in the original KPSS article.
*/

static void kpss_parms (double a, int trend, double *b)
{
    const double b0_level[] = { 0.74375, 0.46119, 0.34732 };
    const double b1_level[] = { -0.99187, 0.45911, 0.20695 };
    const double b0_trend[] = { 0.21778, 0.14795, 0.119298 };
    const double b1_trend[] = { -0.235089, 0.035327, 0.100804 };
    int i = (a == .01)? 0 : (a == .05)? 1 : 2;

    if (trend) {
	b[0] = b0_trend[i];
	b[1] = b1_trend[i];
    } else {
	b[0] = b0_level[i];
	b[1] = b1_level[i];
    }	
}

static double kpss_critval (double alpha, int T, int trend)
{
    double b[2];

    kpss_parms(alpha, trend, b);
    return b[0] + b[1]/T;
}

#define PV_GT10 1.1
#define PV_LT01 -1.0

static double kpss_interp (double s, int T, int trend)
{
    double c10, c05, c01;
    double pv;

    c10 = kpss_critval(.10, T, trend);
    if (s < c10) {
	return PV_GT10;
    }

    c01 = kpss_critval(.01, T, trend);
    if (s > c01) {
	return PV_LT01;
    }  

    /* OK, p-value must lie between .01 and .10 */

    c05 = kpss_critval(.05, T, trend);
    if (s > c05) {
	pv = .01 + .04 * (c01 - s) / (c01 - c05);
    } else {
	pv = .05 + .05 * (c05 - s) / (c05 - c10);
    }

    return pv;
}

static int 
real_kpss_test (int order, int varno, double ***pZ,
		DATAINFO *pdinfo, gretlopt opt, 
		kpss_info *kinfo, PRN *prn)
{
    MODEL KPSSmod;
    int list[4];
    int hastrend = 0;
    double et, s2 = 0.0;
    double cumsum = 0.0, cumsum2 = 0.0;
    double teststat, pval = NADBL;
    double *autocov;
    int t1, t2, T;
    int i, t;

    /* sanity check */
    if (varno <= 0 || varno >= pdinfo->v) {
	return E_DATA;
    }

    if (gretl_isconst(pdinfo->t1, pdinfo->t2, (*pZ)[varno])) {
	gretl_errmsg_sprintf(_("%s is a constant"), pdinfo->varname[varno]);
	return E_DATA;
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
    KPSSmod = lsq(list, *pZ, pdinfo, OLS, OPT_A | OPT_M);
    if (KPSSmod.errcode) {
	clear_model(&KPSSmod);
	return KPSSmod.errcode;
    }

    t1 = KPSSmod.t1;
    t2 = KPSSmod.t2;
    T = KPSSmod.nobs;

    if (order <= 0) {
	order = 4.0 * pow(T / 100.0, 0.25);
    }

    if (kinfo == NULL && (opt & OPT_V)) {
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

    if (s2 <= 0.0) {
	teststat = pval = NADBL;
    } else {
	teststat = cumsum2 / (s2 * T * T);
	pval = kpss_interp(teststat, T, hastrend);
    }

    if (kinfo != NULL) {
	/* storing info for panel test */
	kinfo->T = T;
	kinfo->test = teststat;
	kinfo->pval = pval;
    } else {
	/* testing individual time series */
	if (pval == PV_GT10 || pval == PV_LT01) {
	    pval = NADBL;
	}
	if (opt & OPT_V) {
	    pprintf(prn, "  %s: %g\n", _("Robust estimate of variance"), s2);
	    pprintf(prn, "  %s: %g\n", _("Sum of squares of cumulated residuals"), 
		    cumsum2);
	}
    }

    if (!(opt & OPT_Q)) {
	double a[] = { 0.1, 0.05, 0.01 };
	double cv[3];

	cv[0] = kpss_critval(a[0], T, hastrend);
	cv[1] = kpss_critval(a[1], T, hastrend);
	cv[2] = kpss_critval(a[2], T, hastrend);

	pprintf(prn, _("\nKPSS test for %s %s\n\n"), pdinfo->varname[varno],
		(hastrend)? _("(including trend)") : _("(without trend)"));
	pprintf(prn, "T = %d\n", T);
	pprintf(prn, _("Lag truncation parameter = %d\n"), order);
	pprintf(prn, "%s = %g\n\n", _("Test statistic"), teststat);
	print_critical_values(a, cv, KPSS, prn);
	if (pval == PV_GT10) {
	    pprintf(prn, "%s > .10\n", _("P-value"));
	    pval = NADBL; /* invalidate for record_test_result */
	} else if (pval == PV_LT01) {
	    pprintf(prn, "%s < .01\n", _("P-value"));
	    pval = NADBL;
	} else if (!xna(pval)) {
	    pprintf(prn, "%s %.3f\n", _("Interpolated p-value"), pval);
	}
	pputc(prn, '\n');
    }

    if (kinfo == NULL) {
	record_test_result(teststat, pval, "KPSS");
    }

    clear_model(&KPSSmod);
    free(autocov);

    return 0;
}

static int panel_kpss_test (int order, int v, 
			    double ***pZ, DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn)
{
    kpss_info kinfo;
    int u0 = pdinfo->t1 / pdinfo->pd;
    int uN = pdinfo->t2 / pdinfo->pd;
    int n = uN - u0 + 1;
    int verbose = (opt & OPT_V);
    double ppv = 0.0, zpv = 0.0, lpv = 0.0;
    int gt_10 = 0, lt_01 = 0;
    double pval;
    int i, err = 0;

    /* run a KPSS test for each unit and record the
       results */

    pprintf(prn, _("\nKPSS test for %s %s\n"), pdinfo->varname[v],
	    (opt & OPT_T)? _("(including trend)") : _("(without trend)"));
    pprintf(prn, _("Lag truncation parameter = %d\n"), order);
    pputc(prn, '\n');

    for (i=u0; i<=uN && !err; i++) {
	pdinfo->t1 = i * pdinfo->pd;
	pdinfo->t2 = pdinfo->t1 + pdinfo->pd - 1;
	err = series_adjust_sample((*pZ)[v], &pdinfo->t1, &pdinfo->t2);
	if (!err) {
	    err = real_kpss_test(order, v, pZ, pdinfo, opt | OPT_Q, &kinfo, prn);
	    if (!err && verbose) {
		pprintf(prn, "Unit %d, T = %d\n", i + 1, kinfo.T);
		if (na(kinfo.pval)) {
		    pputs(prn, "\n\n");
		} else {
		    pprintf(prn, "test = %g, ", kinfo.test);
		    if (kinfo.pval == PV_GT10) {
			pprintf(prn, "%s > .10\n", _("p-value"));
		    } else if (kinfo.pval == PV_LT01) {
			pprintf(prn, "%s < .01\n", _("p-value"));
		    } else {
			pprintf(prn, "%s %.3f\n", _("interpolated p-value"), kinfo.pval);
		    }
		    pputc(prn, '\n');
		}
	    }
	}

	if (!err) {
	    pval = kinfo.pval;

	    if (pval == PV_GT10) {
		gt_10++;
		if (lt_01 == 0) {
		    /* record lower bound */
		    pval = .10;
		} else {
		    pval = NADBL;
		}
	    } else if (pval == PV_LT01) {
		lt_01++;
		if (gt_10 == 0) {
		    /* record upper bound */
		    pval = .01;
		} else {
		    pval = NADBL;
		}
	    }

	    if (xna(pval)) {
		ppv = zpv = lpv = NADBL;
	    } else if (!na(ppv)) {
		ppv += log(pval);
		zpv += normal_cdf_inverse(pval);
		lpv += log(pval / (1-pval));
	    }
	}
    }

    if (!err && !na(ppv)) {
	/* process the results as per Choi, as best we can */
	pprintf(prn, "%s\n\n", _("H0: all groups are stationary"));
	do_choi_test(ppv, zpv, lpv, n, prn);
	if (gt_10 > 0) {
	    pputs(prn, "   Note: these are LOWER BOUNDS "
		  "on the true p-values\n");
	    pprintf(prn, "   (Individual p-values > .10, and recorded as .10: %d)\n",
		    gt_10);
	} else if (lt_01 > 0) {
	    pputs(prn, "   Note: these are UPPER BOUNDS "
		  "on the true p-values\n");
	    pprintf(prn, "   (Individual p-values < .01, and recorded as .01: %d)\n",
		    lt_01);
	} 
	pputc(prn, '\n');
    } else {
	pprintf(prn, "Choi test: cannot be calculated\n");
    }

    return err;
}

/**
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @list: list of variables to test.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the KPSS test for 
 * stationarity. Flags that may be given in @opt include:
 * OPT_T to include a linear trend; OPT_F to apply 
 * first-differencing before testing; OPT_V for verbose
 * operation.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int kpss_test (int order, const int *list, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int save_t1 = pdinfo->t1;
    int save_t2 = pdinfo->t2;
    int err = 0;

    if (multi_unit_panel_sample(pdinfo)) {
	err = panel_kpss_test(order, list[1], pZ, pdinfo, opt, prn);
    } else {
	/* regular time series case */
	int i, v, vlist[2] = {1, 0};

	for (i=1; i<=list[0] && !err; i++) {
	    v = list[i];
	    vlist[1] = v;
	    err = list_adjust_sample(vlist, &pdinfo->t1, &pdinfo->t2, 
				     (const double **) *pZ);
	    if (!err) {
		err = real_kpss_test(order, v, pZ, pdinfo, opt, NULL, prn);
	    }
	    pdinfo->t1 = save_t1;
	    pdinfo->t2 = save_t2;
	}
    }

    pdinfo->t1 = save_t1;
    pdinfo->t2 = save_t2;

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

    if (opt & OPT_E) {
	*adf_opt |= OPT_E;
    }

    return 0;
}

/* Engle-Granger: try to ensure a uniform sample for the individual
   (A)DF tests and the test on the cointegrating regression
*/

static int coint_set_sample (const int *list, int nv, int order,
			     double **Z, DATAINFO *pdinfo)
{
    int anymiss;
    int i, t;

    for (t=pdinfo->t1; t<pdinfo->t2; t++) {
	anymiss = 0;
	for (i=1; i<=nv; i++) {
	    if (na(Z[list[i]][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (!anymiss) {
	    break;
	}
    }

    pdinfo->t1 = t + order + 1;
    
    for (t=pdinfo->t2; t>pdinfo->t1; t--) {
	anymiss = 0;
	for (i=1; i<=nv; i++) {
	    if (na(Z[list[i]][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (!anymiss) {
	    break;
	}
    }

    pdinfo->t2 = t;

    return 0;
}

/**
 * engle_granger_test:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Carries out the Engle-Granger test for cointegration. 
 * Flags that may be given in @opt include: OPT_N, do
 * not an include a constant in the cointegrating regression;
 * OPT_T include constant and linear trend; OPT_R, include
 * quadratic trend; OPT_S, skip DF tests for individual variables; 
 * OPT_E, test down from maximum lag order (see the entry for
 * "adf" in the Gretl Command Reference for details).
 *
 * Returns: 0 on successful completion, non-zero code
 * on error.
 */

int engle_granger_test (int order, const int *list, double ***pZ, 
			DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int test_t1, test_t2;
    gretlopt adf_opt = OPT_C;
    MODEL cmod;
    int detcode = UR_CONST;
    int i, nv, k = 0;
    int step = 1;
    int *clist = NULL;
    int err = 0;

    if (multi_unit_panel_sample(pdinfo)) {
	gretl_errmsg_set("Sorry, this command is not yet available "
			 "for panel data");
	return E_DATA;
    }

    err = coint_check_opts(opt, &detcode, &adf_opt);
    if (err) {
	return err;
    }

    clist = make_coint_list(list, detcode, &nv, pZ, pdinfo, &err);
    if (err) {
	return err;
    }

    /* backward compatibility: let a negative lag order
       indicate that we should test down */
    if (order < 0) {
	order = -order;
	adf_opt |= OPT_E;
    }

    gretl_model_init(&cmod);

    if (!(opt & OPT_S)) {
	/* test all candidate vars for unit root */
	coint_set_sample(clist, nv, order, *pZ, pdinfo);
	for (i=1; i<=nv; i++) {
	    if (step == 1) {
		pputc(prn, '\n');
	    }
	    pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		    step++, pdinfo->varname[clist[i]]);
	    real_adf_test(clist[i], order, 1, pZ, pdinfo, adf_opt, 
			  ADF_EG_TEST, NULL, prn);
	}
    }

    if (step == 1) {
	pputc(prn, '\n');
    }
    pprintf(prn, _("Step %d: cointegrating regression\n"), step++);

    test_t1 = pdinfo->t1;
    test_t2 = pdinfo->t2;

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    cmod = lsq(clist, *pZ, pdinfo, OLS, OPT_NONE);
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

    pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
	    step, pdinfo->varname[k]);

    pdinfo->t1 = test_t1;
    pdinfo->t2 = test_t2;

    /* Run (A)DF test on the residuals */
    real_adf_test(k, order, nv, pZ, pdinfo, adf_opt, 
		  ADF_EG_TEST | ADF_EG_RESIDS, NULL, prn); 

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		 "(a) The unit-root hypothesis is not rejected for the individual"
		 " variables.\n(b) The unit-root hypothesis is rejected for the "
		 "residuals (uhat) from the \n    cointegrating regression.\n"));
    pputc(prn, '\n');

 bailout:
    
    clear_model(&cmod);
    free(clist);
    if (k > 0) {
	dataset_drop_variable(k, pZ, pdinfo);
    }

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return err;
}
