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

#include "libgretl.h" 
#include "var.h"  
#include "vartest.h"
#include "matrix_extra.h"
#include "libset.h"

int gretl_VAR_normality_test (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->E == NULL || var->S == NULL) {
	err = 1;
    } else {
	err = multivariate_normality_test(var->E, var->S, prn);
    }

    return err;
}

int gretl_VAR_autocorrelation_test (GRETL_VAR *var, int order, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    int i, err = 0;

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
	err = autocorr_test(var->models[i], order, pZ, pdinfo,
			    OPT_Q | OPT_S, prn);
	gretl_model_test_print(var->models[i], 0, prn);
	gretl_model_destroy_tests(var->models[i]);
    }

    return err;
}

int gretl_VAR_arch_test (GRETL_VAR *var, int order, 
			 DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
	err = arch_test(var->models[i], order, pdinfo, OPT_NONE, prn);
    }

    return err;
}

static void 
form_C0j (const GRETL_VAR *var, gretl_matrix *C0j, 
	  gretl_matrix *et, gretl_matrix *ej,
	  int j)
{
    int i, t;

    gretl_matrix_zero(C0j);

    for (t=j; t<var->T; t++) {
	/* load e_t and e_{t-j} */
	for (i=0; i<var->neqns; i++) {
	    et->val[i] = gretl_matrix_get(var->E, t, i);
	    ej->val[i] = gretl_matrix_get(var->E, t-j, i);
	}
	/* add e_t' * e_{t-j} */
	gretl_matrix_multiply_mod(et, GRETL_MOD_TRANSPOSE,
				  ej, GRETL_MOD_NONE,
				  C0j, GRETL_MOD_CUMULATE);
    }

    gretl_matrix_divide_by_scalar(C0j, var->T);
}

/* See Johansen (1995 book) pp. 21-22 */

int VAR_portmanteau_test (GRETL_VAR *var)
{
    gretl_matrix *C00 = NULL;
    gretl_matrix *C0j = NULL;
    gretl_matrix *et = NULL;
    gretl_matrix *ej = NULL;
    gretl_matrix *L = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *Tmp = NULL;
    int n = var->neqns;
    double trj, LB = 0.0;
    int s, j, err = 0;

    /* we'll do this just for unrestricted VARs */
    if (var->ci == VECM && jrank(var) < var->neqns) {
	return 0;
    }

    /* any guidance on the order for this test? */
    s = var->T / 4;
    if (s > 48) s = 48;

    clear_gretl_matrix_err();

    C00 = gretl_matrix_alloc(n, n);
    C0j = gretl_matrix_alloc(n, n);
    et = gretl_matrix_alloc(1, n);
    ej = gretl_matrix_alloc(1, n);
    L = gretl_matrix_alloc(n, n);
    R = gretl_matrix_alloc(n, n);
    Tmp = gretl_matrix_alloc(n, n);

    err = get_gretl_matrix_err();

    form_C0j(var, C00, et, ej, 0);

    if (!err) {
	err = gretl_invert_symmetric_matrix(C00);
    }

    for (j=1; j<=s && !err; j++) {
	form_C0j(var, C0j, et, ej, j);
	gretl_matrix_multiply(C0j, C00, L);
	gretl_matrix_multiply_mod(C0j, GRETL_MOD_TRANSPOSE,
				  C00, GRETL_MOD_NONE,
				  R, GRETL_MOD_NONE);
	gretl_matrix_multiply(L, R, Tmp);
	trj = gretl_matrix_trace(Tmp, &err);
	LB += (1.0 / (var->T - j)) * trj;
    }

    if (!err) {
	LB *= var->T * (var->T + 2);
	var->LB = LB;
	var->LBs = s;
    } else {
	var->LB = NADBL;
	var->LBs = 0;
    }

    gretl_matrix_free(C00);
    gretl_matrix_free(C0j);
    gretl_matrix_free(et);
    gretl_matrix_free(ej);
    gretl_matrix_free(L);
    gretl_matrix_free(R);
    gretl_matrix_free(Tmp);

    return err;
}

/* make and record residuals for LR test on last lag */

int last_lag_LR_prep (GRETL_VAR *var, int ifc)
{
    gretl_matrix *X = NULL;
    gretl_matrix *B = NULL;
    double x;
    int g = var->ncoeff - var->neqns;
    int lag, lmax;
    int i, t, k, err = 0;

    if (var->F == NULL) {
	var->F = gretl_matrix_alloc(var->T, var->neqns);
	if (var->F == NULL) {
	    return E_ALLOC;
	}
    }   

    X = gretl_matrix_alloc(var->T, g);
    B = gretl_matrix_alloc(g, var->neqns);
    if (X == NULL || B == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    lag = (ifc)? 0 : 1;
    lmax = ifc + var->neqns * var->order;
    k = 0;

    for (i=0; i<var->ncoeff; i++) {
	if (lag != var->order) {
	    for (t=0; t<var->T; t++) {
		x = gretl_matrix_get(var->X, t, i);
		gretl_matrix_set(X, t, k, x);
	    }
	    k++;
	}
	lag = (i >= lmax)? 0 : (lag < var->order)? lag + 1 : 1;
    }

    err = gretl_matrix_multi_ols(var->Y, X, B, var->F, NULL);

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(B);

    return err;
}

int VAR_LR_lag_test (GRETL_VAR *var)
{
    double ldet;
    int err = 0;

    ldet = gretl_VAR_ldet(var, &err);

    if (!err) {
	double ll, AIC, BIC, HQC;
	int T = var->T;
	int g = var->neqns;
	int m = var->ncoeff - g;
	int k = g * m;

	var->LR = T * (ldet - var->ldet);

	ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	AIC = (-2.0 * ll + 2.0 * k) / T;
	BIC = (-2.0 * ll + k * log(T)) / T;
	HQC = (-2.0 * ll + 2.0 * k * log(log(T))) / T;
	var->Ivals[0] = AIC;
	var->Ivals[1] = BIC;
	var->Ivals[2] = HQC;
    }

    /* we're done with this set of residuals */
    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

static void gretl_VAR_print_lagsel (gretl_matrix *lltab,
				    gretl_matrix *crittab,
				    int *best_row,
				    PRN *prn)
{
    int maxlag = gretl_matrix_rows(crittab);
    double x;
    int i, j;

    pprintf(prn, _("VAR system, maximum lag order %d"), maxlag);
    pputs(prn, "\n\n");

    pputs(prn, _("The asterisks below indicate the best (that is, minimized) values\n"
	  "of the respective information criteria, AIC = Akaike criterion,\n"
	  "BIC = Schwartz Bayesian criterion and HQC = Hannan-Quinn criterion."));
    pputs(prn, "\n\n");

    pputs(prn, _("lags        loglik    p(LR)       AIC          BIC          HQC"));
    pputs(prn, "\n\n");

    for (i=0; i<maxlag; i++) {
	pprintf(prn, "%4d", i + 1);
	x = gretl_matrix_get(lltab, i, 0);
	pprintf(prn, "%14.5f", x);
	if (i > 0) {
	    x = gretl_matrix_get(lltab, i, 1);
	    pprintf(prn, "%9.5f", x);
	} else {
	    pputs(prn, "         ");
	}
	for (j=0; j<N_IVALS; j++) {
	    x = gretl_matrix_get(crittab, i, j);
	    pprintf(prn, "%12.6f", x);
	    if (i == best_row[j]) {
		pputc(prn, '*');
	    } else {
		pputc(prn, ' ');
	    }
	}
	pputc(prn, '\n');
    }
}

/* apparatus for selecting the optimal lag length for a VAR */

int VAR_do_lagsel (GRETL_VAR *var, const double **Z, 
		   const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *crittab = NULL;
    gretl_matrix *lltab = NULL;

    int p = var->order;
    int r = p - 1;
    int T = var->T;
    int n = var->neqns;

    /* initialize the "best" at the longest lag */
    double best[N_IVALS] = { var->AIC, var->BIC, var->HQC };
    int best_row[N_IVALS] = { r, r, r };
    double crit[N_IVALS];
    double LRtest;
    double ldet = NADBL;
    int cols0;
    int j, m = 0;
    int err = 0;

    if (p < 2) {
	return 0;
    }

    if (var->F != NULL) {
	gretl_matrix_free(var->F);
    }

    var->F = gretl_matrix_alloc(T, n);
    if (var->F == NULL) {
	return E_ALLOC;
    }

    crittab = gretl_matrix_alloc(p, N_IVALS);
    lltab = gretl_matrix_alloc(p, 2);
    if (crittab == NULL || lltab == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* # of cols in X that are not Y lags */
    cols0 = var->ncoeff - p * n; 

    for (j=1; j<p && !err; j++) {
	int jxcols = cols0 + j * n;

	VAR_fill_X(var, j, Z, pdinfo);

	gretl_matrix_reuse(var->X, T, jxcols);
	gretl_matrix_reuse(var->B, jxcols, n);

	err = gretl_matrix_multi_ols(var->Y, var->X, var->B, 
				     var->F, NULL);

	if (!err) {
	    ldet = gretl_VAR_ldet(var, &err);
	}

	if (!err) {
	    double ll;
	    int q = var->ncoeff - (n * (p - j));
	    int c, k = n * q;

	    ll = -(n * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	    crit[0] = (-2.0 * ll + 2.0 * k) / T;               /* AIC */
	    crit[1] = (-2.0 * ll + k * log(T)) / T;            /* BIC */
	    crit[2] = (-2.0 * ll + 2.0 * k * log(log(T))) / T; /* HQC */

	    gretl_matrix_set(lltab, m, 0, ll);
	    if (j == 1) {
		gretl_matrix_set(lltab, m, 1, 0);
	    } else {
		LRtest = 2.0 * (ll - gretl_matrix_get(lltab, m-1, 0));
		gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, n * n));
	    }	
	    
	    for (c=0; c<N_IVALS; c++) {
		gretl_matrix_set(crittab, m, c, crit[c]);
		if (crit[c] < best[c]) {
		    best[c] = crit[c];
		    best_row[c] = m;
		}
	    }
	
	    m++;
	}
    }

    if (!err) {
	gretl_matrix_set(lltab, m, 0, var->ll);
	LRtest = 2.0 * (var->ll - gretl_matrix_get(lltab, m - 1, 0));
	gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, n * n));
	gretl_matrix_set(crittab, m, 0, var->AIC);
	gretl_matrix_set(crittab, m, 1, var->BIC);
	gretl_matrix_set(crittab, m, 2, var->HQC);
	gretl_VAR_print_lagsel(lltab, crittab, best_row, prn);
    }

    bailout:

    gretl_matrix_free(crittab);
    gretl_matrix_free(lltab);

    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

/* (X'X)^{-1} * X'\Omega X * (X'X)^{-1} : right now
   we're supporting only HC0 and HC1 */

static int VAR_robust_vcv (GRETL_VAR *var, gretl_matrix *V,
			   MODEL *pmod, int hcv, int k)
{
    gretl_matrix *XOX = NULL;
    double xij, xti, xtj, utk;
    int T = var->T;
    int g = var->ncoeff;
    int i, j, t;

    XOX = gretl_matrix_alloc(g, g);
    if (XOX == NULL) {
	return E_ALLOC;
    }

    /* form X' \Omega X */
    for (i=0; i<g; i++) {
	for (j=i; j<g; j++) {
	    xij = 0.0;
	    for (t=0; t<T; t++) {
		xti = gretl_matrix_get(var->X, t, i);
		xtj = gretl_matrix_get(var->X, t, j);
		utk = gretl_matrix_get(var->E, t, k);
		xij += utk * utk * xti * xtj;
	    }
	    if (hcv > 0) {
		/* cheating here, for now */
		xij *= (double) T / (T - g);
	    }
	    gretl_matrix_set(XOX, i, j, xij);
	    if (i != j) {
		gretl_matrix_set(XOX, j, i, xij);
	    }
	}
    }

    gretl_matrix_qform(var->XTX, GRETL_MOD_TRANSPOSE, XOX,
		       V, GRETL_MOD_NONE);

    gretl_model_set_int(pmod, "hc", 1);
    if (hcv > 0) {
	gretl_model_set_int(pmod, "hc_version", 1);
    }

    gretl_matrix_free(XOX);

    return 0;
}

/* Run the various per-equation omit tests (all lags of each var in
   turn, last lag of all vars) using the Wald method.  We also
   add the standard errors to the models here, since we have the
   covariance matrix to hand.
*/

int VAR_wald_omit_tests (GRETL_VAR *var, int ifc)
{
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    int hcv = libset_get_int(HC_VERSION);
    int p = var->order;
    int n = var->neqns;
    int g = var->ncoeff;
    int dim = (p > n)? p : n;
    int i, j, k, m = 0;
    int err = 0;

    if (ifc && var->robust && g - 1 > dim) {
	/* need bigger arrays for robust overall F-test */
	dim = g - 1;
    }

    V = gretl_matrix_alloc(g, g);
    C = gretl_matrix_alloc(dim, dim);
    b = gretl_column_vector_alloc(dim);

    if (V == NULL || C == NULL || b == NULL) {
	return E_ALLOC;
    }     

    for (i=0; i<n && !err; i++) {
	MODEL *pmod = var->models[i];
	int ii, jj, jpos, ipos = ifc;
	double w, vij;

	gretl_matrix_reuse(V, g, g);

	if (var->robust) {
	    err = VAR_robust_vcv(var, V, pmod, hcv, i);
	} else {
	    gretl_matrix_copy_values(V, var->XTX);
	    gretl_matrix_multiply_by_scalar(V, pmod->sigma * pmod->sigma);
	}
	
	if (!err) {
	    /* set (possibly robust) standard errors */
	    for (j=0; j<g; j++) {
		vij = gretl_matrix_get(V, j, j);
		pmod->sderr[j] = sqrt(vij);
	    }
	}

	/* exclusion of each var, all lags */

	gretl_matrix_reuse(C, p, p);
	gretl_matrix_reuse(b, p, 1);

	for (j=0; j<n && !err; j++) {
	    double w = NADBL;

	    gretl_matrix_extract_matrix(C, V, ipos, ipos, GRETL_MOD_NONE);
	    for (k=0; k<p; k++) {
		b->val[k] = pmod->coeff[k + ipos];
	    }
	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		var->Fvals[m++] = w / p;
	    }

	    ipos += p;
	}

	/* exclusion of last lag, all vars? */

	if (p > 1) {
	    gretl_matrix_reuse(C, n, n);
	    gretl_matrix_reuse(b, n, 1);

	    ipos = ifc + p - 1;
	    for (ii=0; ii<n; ii++) {
		jpos = ifc + p - 1;
		for (jj=0; jj<n; jj++) {
		    vij = gretl_matrix_get(V, ipos, jpos);
		    gretl_matrix_set(C, ii, jj, vij);
		    jpos += p;
		}
		b->val[ii] = pmod->coeff[ipos];
		ipos += p;
	    }

	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		var->Fvals[m++] = w / n;
	    }
	}

	/* exclusion of all but const? */

	if (ifc && var->robust) {
	    gretl_matrix_reuse(C, g-1, g-1);
	    gretl_matrix_reuse(b, g-1, 1);

	    gretl_matrix_extract_matrix(C, V, 1, 1, GRETL_MOD_NONE);
	    for (k=0; k<g-1; k++) {
		b->val[k] = pmod->coeff[k+1];
	    }
	    err = gretl_invert_symmetric_matrix(C);
	    if (!err) {
		w = gretl_scalar_qform(b, C, &err);
	    }
	    if (!err) {
		pmod->fstt = w / (g-1);
	    }
	}
    }

    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return err;
}

#define VO_DEBUG 0

const int *gretl_VAR_get_exo_list (const GRETL_VAR *var)
{
    return var->xlist;
}

/* Based on the specification stored in the original VAR struct, 
   plus a new exogenous list, constitute a full VAR list.
*/

static int *build_VAR_list (const GRETL_VAR *var, int *exolist, int *err)
{
    return list_composite(var->ylist, exolist, var->rlist);
}

static int gretl_VAR_real_omit_test (const GRETL_VAR *orig,
				     const GRETL_VAR *new,
				     const DATAINFO *pdinfo,
				     PRN *prn)
{
    int *omitlist;
    double LR, pval;
    int i, df, err = 0;

#if VO_DEBUG
    fprintf(stderr, "gretl_VAR_real_omit_test: about to diff lists\n");
    printlist(orig->xlist, "orig xlist");
    printlist(new->xlist, "new xlist");
#endif

    if (new->xlist == NULL) {
	omitlist = gretl_list_copy(orig->xlist);
    } else {
	omitlist = gretl_list_diff_new(orig->xlist, new->xlist, 1);
    }

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
 * @Z: data array.
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
				const double **Z, DATAINFO *pdinfo, 
				PRN *prn, int *err)
{
    GRETL_VAR *var = NULL;
    gretlopt opt = OPT_NONE;
    int smpl_t1 = pdinfo->t1;
    int smpl_t2 = pdinfo->t2;
    int *tmplist = NULL;
    int *varlist = NULL;
    int c1 = 0;

    *err = 0;

    if (orig == NULL || orig->xlist == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (omitvars == NULL || omitvars[0] == 0) {
	*err = E_PARSE;
	return NULL;
    }

#if VO_DEBUG
    printlist(orig->xlist, "original xlist");
#endif

    if (orig->ifc) {
	c1 = !gretl_list_const_pos(omitvars, 1, Z, pdinfo);
    } 

    /* create reduced exogenous vars list for test VAR */
    tmplist = gretl_list_omit(orig->xlist, omitvars, 1, err);
    if (tmplist == NULL) {
	goto bailout;
    }

#if VO_DEBUG
    fprintf(stderr, "c1 = %d\n", c1);
    printlist(tmplist, "exog vars list for test VAR");
#endif

    /* create full input VAR list for test VAR */
    varlist = build_VAR_list(orig, tmplist, err);
    if (varlist == NULL) {
	goto bailout;
    }

#if VO_DEBUG
    printlist(varlist, "full list for test VAR");
#endif

    if (orig->detflags & DET_SEAS) {
	opt |= OPT_D;
    }

    if (orig->detflags & DET_TREND) {
	opt |= OPT_T;
    }

    /* If the original VAR did not include a constant, we need to
       pass OPT_N to the test VAR to suppress the constant.
       We also need to pass OPT_N in case the constant was
       present originally but is now to be omitted.
    */
    if (orig->ifc == 0 || c1 == 0) {
	opt |= OPT_N;
    }

    /* impose as sample range the estimation range of the 
       original VAR */
    pdinfo->t1 = orig->t1;
    pdinfo->t2 = orig->t2;

    var = gretl_VAR(orig->order, varlist, Z, pdinfo, opt, prn, err);

    /* now, if var is non-NULL, do the actual test(s) */
    if (var != NULL) {
	*err = gretl_VAR_real_omit_test(orig, var, pdinfo, prn);
    }

    /* put back into pdinfo what was there on input */
    pdinfo->t1 = smpl_t1;
    pdinfo->t2 = smpl_t2;

 bailout:

    free(tmplist);
    free(varlist);

    return var;
}

