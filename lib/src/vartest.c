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
#include "johansen.h"
#include "vartest.h"
#include "matrix_extra.h"
#include "libset.h"
#include "qr_estimate.h"
#include "varprint.h"

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
				    DATASET *dset, gretlopt opt,
				    PRN *prn)
{
    MODEL *pmod;
    gretl_matrix *tests, *pvals;
    double lb;
    int i, quiet = (opt & OPT_Q);
    int err = 0;

    if (order == 0) {
	order = dset->pd;
    }

    tests = gretl_column_vector_alloc(var->neqns);
    pvals = gretl_column_vector_alloc(var->neqns);

    if (tests == NULL || pvals == NULL) {
	err = E_ALLOC;
    }

    for (i=0; i<var->neqns && !err; i++) {
	pmod = var->models[i];
	if (!quiet) {
	    pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
	}
	if (pmod->list[0] == 1) {
	    /* only the dependent variable is recorded */
	    lb = ljung_box(order, pmod->t1, pmod->t2, pmod->uhat, &err);
	    if (!err) {
		tests->val[i] = lb;
		pvals->val[i] = chisq_cdf_comp(order, lb);
		if (!(opt & OPT_Q)) {
		    pprintf(prn, "Ljung-Box Q' = %g %s = P(%s(%d) > %g) = %.3g\n", 
			    lb, _("with p-value"), _("Chi-square"), order,
			    lb, pvals->val[i]);
		    pputc(prn, '\n');
		}
	    }
	} else {
	    if (quiet) {
		err = autocorr_test(pmod, order, dset, OPT_Q, prn);
	    } else {
		err = autocorr_test(pmod, order, dset, OPT_Q | OPT_S, prn);
	    }
	    if (!err) {
		tests->val[i] = get_last_test_statistic(NULL);
		pvals->val[i] = get_last_pvalue(NULL);
		if (!quiet) {
		    gretl_model_test_print(var->models[i], 0, prn);
		    gretl_model_destroy_tests(var->models[i]);
		}
	    }
	}
    }

    if (!err) {
	record_matrix_test_result(tests, pvals);
    } else {
	gretl_matrix_free(tests);
	gretl_matrix_free(pvals);
    }

    return err;
}

int gretl_VAR_arch_test (GRETL_VAR *var, int order, 
			 DATASET *dset, gretlopt opt,
			 PRN *prn)
{
    gretl_matrix *tests, *pvals;
    int i, err = 0;

    if (order == 0) {
	order = dset->pd;
    }

    tests = gretl_column_vector_alloc(var->neqns);
    pvals = gretl_column_vector_alloc(var->neqns);

    if (tests == NULL || pvals == NULL) {
	err = E_ALLOC;
    } else {
	pprintf(prn, "%s %d\n\n", _("Test for ARCH of order"), 
		order);
    }

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
	/* add OPT_M for multi-equation output */
	err = arch_test(var->models[i], order, dset, 
			opt | OPT_M, prn);
	if (!err) {
	    tests->val[i] = get_last_test_statistic(NULL);
	    pvals->val[i] = get_last_pvalue(NULL);
	}
    }

    if (!err) {
	record_matrix_test_result(tests, pvals);
    } else {
	gretl_matrix_free(tests);
	gretl_matrix_free(pvals);
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
    gretl_matrix_block *B;
    gretl_matrix *C00, *C0j;
    gretl_matrix *et, *ej;
    gretl_matrix *L, *R, *Tmp;
    int k, n = var->neqns;
    double trj, LB = 0.0;
    int s, j, err = 0;

    var->LB = NADBL;
    var->LBs = 0;

    /* we'll do this only for unrestricted VARs */
    if (var->ci == VECM && jrank(var) < var->neqns) {
	return 0;
    }

    /* any guidance on the order for this test? */
    s = var->T / 4;
    if (s > 48) s = 48;

    k = var->order + (var->ci == VECM);
    if (s - k <= 0) {
	/* no degrees of freedom */
	return 0;
    }

    B = gretl_matrix_block_new(&C00, n, n,
			       &C0j, n, n,
			       &et,  1, n,
			       &ej,  1, n,
			       &L,   n, n,
			       &R,   n, n,
			       &Tmp, n, n,
			       NULL);

    if (B == NULL) {
	return E_ALLOC;
    }

    form_C0j(var, C00, et, ej, 0);
    err = gretl_invert_symmetric_matrix(C00);

    for (j=1; j<=s && !err; j++) {
	form_C0j(var, C0j, et, ej, j);
	gretl_matrix_multiply(C0j, C00, L);
	gretl_matrix_multiply_mod(C0j, GRETL_MOD_TRANSPOSE,
				  C00, GRETL_MOD_NONE,
				  R, GRETL_MOD_NONE);
	gretl_matrix_multiply(L, R, Tmp);
	trj = gretl_matrix_trace(Tmp);
	LB += (1.0 / (var->T - j)) * trj;
    }

    if (!err) {
	LB *= var->T * (var->T + 2);
	var->LB = LB;
	var->LBs = s;
    }

    gretl_matrix_block_destroy(B);

    return err;
}

int VAR_LR_lag_test (GRETL_VAR *var, const gretl_matrix *E)
{
    double test_ldet;
    int err = 0;

    test_ldet = gretl_VAR_ldet(var, E, &err);

    if (!err) {
	double ll, AIC, BIC, HQC;
	int T = var->T;
	int g = var->neqns;
	int m = var->ncoeff - g;
	int k = g * m;

	var->LR = T * (test_ldet - var->ldet);

	ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * test_ldet;
	AIC = (-2.0 * ll + 2.0 * k) / T;
	BIC = (-2.0 * ll + k * log(T)) / T;
	HQC = (-2.0 * ll + 2.0 * k * log(log(T))) / T;
	var->Ivals[0] = AIC;
	var->Ivals[1] = BIC;
	var->Ivals[2] = HQC;
    }

    return err;
}

static void gretl_VAR_print_lagsel (int minlag,
				    gretl_matrix *lltab,
				    gretl_matrix *crittab,
				    int *best_row,
				    PRN *prn)
{
    int maxlag, nrows = gretl_matrix_rows(crittab);
    double x;
    int i, j;

    maxlag = nrows + minlag - 1;

    pprintf(prn, _("VAR system, maximum lag order %d"), maxlag);
    pputs(prn, "\n\n");

    pputs(prn, _("The asterisks below indicate the best (that is, minimized) values\n"
	  "of the respective information criteria, AIC = Akaike criterion,\n"
	  "BIC = Schwarz Bayesian criterion and HQC = Hannan-Quinn criterion."));
    pputs(prn, "\n\n");

    pputs(prn, _("lags        loglik    p(LR)       AIC          BIC          HQC"));
    pputs(prn, "\n\n");

    for (i=0; i<nrows; i++) {
	pprintf(prn, "%4d", minlag + i);
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
    pputc(prn, '\n');
}

static int lagsel_get_min_lag (int p, int *err)
{
    int m = get_optval_int(VAR, OPT_M, err);

    if (!*err && (m < 0 || m > p)) {
	*err = E_INVARG;
    }

    return m;
}

/* apparatus for selecting the optimal lag length for a VAR */

int VAR_do_lagsel (GRETL_VAR *var, const DATASET *dset, 
		   gretlopt opt, PRN *prn)
{
    gretl_matrix *crittab = NULL;
    gretl_matrix *lltab = NULL;
    gretl_matrix *E = NULL;
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
    int cols0, minlag = 1;
    int nrows, use_QR = 0;
    int j, m = 0;
    int err = 0;

    /* number of cols in X that are not Y lags */
    cols0 = var->ncoeff - p * n;

    if (opt & OPT_M) {
	minlag = lagsel_get_min_lag(p, &err);
    }

    if (p < minlag + 1) {
	return 0;
    }

    E = gretl_matrix_alloc(T, n);
    if (E == NULL) {
	return E_ALLOC;
    }

    nrows = p - minlag + 1;
    crittab = gretl_matrix_alloc(nrows, N_IVALS);
    lltab = gretl_matrix_alloc(nrows, 2);

    if (crittab == NULL || lltab == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (getenv("VAR_USE_QR") != NULL) {
	use_QR = 1;
    }

    for (j=minlag; j<p && !err; j++) {
	int jxcols = cols0 + j * n;

	if (jxcols == 0) {
	    gretl_matrix_copy_values(E, var->Y);
	} else {
	    VAR_fill_X(var, j, dset);

	    gretl_matrix_reuse(var->X, T, jxcols);
	    gretl_matrix_reuse(var->B, jxcols, n);

	    if (use_QR) {
		err = gretl_matrix_QR_ols(var->Y, var->X, var->B, 
					  E, NULL, NULL);
	    } else {
		err = gretl_matrix_multi_ols(var->Y, var->X, var->B, 
					     E, NULL);
	    }
	}

	if (!err) {
	    ldet = gretl_VAR_ldet(var, E, &err);
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
	    if (j == minlag) {
		gretl_matrix_set(lltab, m, 1, 0);
	    } else {
		LRtest = 2.0 * (ll - gretl_matrix_get(lltab, m-1, 0));
		gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(n * n, LRtest));
	    }	
	    
	    for (c=0; c<N_IVALS; c++) {
		gretl_matrix_set(crittab, m, c, crit[c]);
		if (crit[c] < best[c]) {
		    best[c] = crit[c];
		    best_row[c] = m;
		}
	    }
	
	    m++; /* increment table row */
	}
    }

    if (!err) {
	gretl_matrix_set(lltab, m, 0, var->ll);
	LRtest = 2.0 * (var->ll - gretl_matrix_get(lltab, m - 1, 0));
	gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(n * n, LRtest));
	gretl_matrix_set(crittab, m, 0, var->AIC);
	gretl_matrix_set(crittab, m, 1, var->BIC);
	gretl_matrix_set(crittab, m, 2, var->HQC);
	if (!(opt & OPT_S)) {
	    gretl_VAR_print_lagsel(minlag, lltab, crittab, best_row, prn);
	}
	record_matrix_test_result(crittab, NULL);
	crittab = NULL;
    }

 bailout:

    gretl_matrix_free(crittab);
    gretl_matrix_free(lltab);
    gretl_matrix_free(E);

    return err;
}

static gretl_matrix *VAR_get_hvec (const gretl_matrix *X,
				   const gretl_matrix *XTX,
				   int *err)
{
    gretl_matrix *hvec;
    gretl_matrix *Xt;
    int t, T = X->rows;
    int k = X->cols;
    double ht;

    Xt = gretl_matrix_alloc(1, k);
    hvec = gretl_column_vector_alloc(T);

    if (Xt == NULL || hvec == NULL) {
	gretl_matrix_free(Xt);
	gretl_matrix_free(hvec);
	*err = E_ALLOC;
	return NULL;
    }

    for (t=0; t<T && !*err; t++) {
	gretl_matrix_get_row(X, t, Xt);
	ht = gretl_scalar_qform(Xt, XTX, err);
	gretl_vector_set(hvec, t, ht);
    }

    gretl_matrix_free(Xt);

    return hvec;
}

static gretl_matrix *var_hac_xox (GRETL_VAR *var, int k,
				  VCVInfo *vi, int *err)
{
    gretl_matrix *XOX = NULL;
    gretl_matrix *uk;
    int t;

    uk = gretl_column_vector_alloc(var->T);

    if (uk == NULL) {
	*err = E_ALLOC;
    } else {
	for (t=0; t<var->T; t++) {
	    uk->val[t] = gretl_matrix_get(var->E, t, k);
	}
	XOX = HAC_XOX(uk, var->X, vi, err);
	gretl_matrix_free(uk);
    }

    return XOX;
}

static gretl_matrix *var_hc_xox (GRETL_VAR *var, int k,
				 int hcv, int *err)
{
    gretl_matrix *XOX;
    gretl_matrix *hvec = NULL;
    int T = var->T;
    int g = var->ncoeff;

    XOX = gretl_matrix_alloc(g, g);

    if (XOX == NULL) {
	*err = E_ALLOC;
    } else if (hcv > 1) {
	hvec = VAR_get_hvec(var->X, var->XTX, err);
    }

    if (!*err) {
	/* form X' \Omega X */
	double xti, xtj, xij;
	double utk, u2, ht;
	int i, j, t;

	for (i=0; i<g; i++) {
	    for (j=i; j<g; j++) {
		xij = 0.0;
		for (t=0; t<T; t++) {
		    utk = gretl_matrix_get(var->E, t, k);
		    u2 = utk * utk;
		    if (hcv > 1) {
			ht = gretl_vector_get(hvec, t);
			u2 /= 1.0 - ht;
			if (hcv > 2) {
			    u2 /= 1.0 - ht;
			}
		    }
		    xti = gretl_matrix_get(var->X, t, i);
		    xtj = gretl_matrix_get(var->X, t, j);
		    xij += u2 * xti * xtj;
		}
		if (hcv == 1) {
		    xij *= (double) T / (T - g);
		} 
		gretl_matrix_set(XOX, i, j, xij);
		if (i != j) {
		    gretl_matrix_set(XOX, j, i, xij);
		}
	    }
	}
    }

    gretl_matrix_free(hvec);

    return XOX;
}

/* (X'X)^{-1} * X'\Omega X * (X'X)^{-1} */

static int VAR_robust_vcv (GRETL_VAR *var, gretl_matrix *V,
			   MODEL *pmod, int hcv, int k)
{
    gretl_matrix *XOX = NULL;
    VCVInfo vi = {0};
    int err = 0;

    if (var->robust == VAR_HAC) {
	XOX = var_hac_xox(var, k, &vi, &err);
    } else {
	XOX = var_hc_xox(var, k, hcv, &err);
    }

    if (!err) {
	gretl_matrix_qform(var->XTX, GRETL_MOD_TRANSPOSE, XOX,
			   V, GRETL_MOD_NONE);

	if (var->robust == VAR_HAC) {
	    gretl_model_set_full_vcv_info(pmod, VCV_HAC, vi.vmin,
					  vi.order, vi.flags,
					  vi.bw);
	} else {
	    gretl_model_set_vcv_info(pmod, VCV_HC, hcv);
	}
    }

    gretl_matrix_free(XOX);

    return err;
}

/* Run the various per-equation omit tests (all lags of each var in
   turn, last lag of all vars) using the Wald method.  We also
   add the standard errors to the models here, since we have the
   per-equation covariance matrices to hand.
*/

int VAR_wald_omit_tests (GRETL_VAR *var)
{
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    int hcv = libset_get_int(HC_VERSION);
    int p = (var->lags != NULL)? var->lags[0] : var->order;
    int n = var->neqns;
    int g = var->ncoeff;
    int dim = (p > n)? p : n;
    int i, j, k, m = 0;
    int any_F_err = 0;
    int err = 0;

    if (var->ifc && var->robust && g - 1 > dim) {
	/* need bigger arrays for robust overall F-test */
	dim = g - 1;
    }

    V = gretl_matrix_alloc(g, g);
    C = gretl_matrix_alloc(dim, dim);
    b = gretl_column_vector_alloc(dim);

    if (V == NULL || C == NULL || b == NULL) {
	return E_ALLOC;
    }     

    for (i=0; i<var->neqns && !err; i++) {
	MODEL *pmod = var->models[i];
	int ii, jj, jpos, ipos = var->ifc;
	double w, vij;
	int F_err = 0;

	gretl_matrix_reuse(V, g, g);

	if (var->robust) {
	    err = VAR_robust_vcv(var, V, pmod, hcv, i);
	} else {
	    gretl_matrix_copy_values(V, var->XTX);
	    gretl_matrix_multiply_by_scalar(V, pmod->sigma * pmod->sigma);
	}

	if (err) {
	    break;
	}
	
	/* set (possibly robust) standard errors */
	for (j=0; j<g; j++) {
	    vij = gretl_matrix_get(V, j, j);
	    pmod->sderr[j] = sqrt(vij);
	}

	/* exclusion of each var, all lags */

	gretl_matrix_reuse(C, p, p);
	gretl_matrix_reuse(b, p, 1);

	for (j=0; j<var->neqns; j++) {
	    double w = NADBL;

	    gretl_matrix_extract_matrix(C, V, ipos, ipos, GRETL_MOD_NONE);
	    for (k=0; k<p; k++) {
		b->val[k] = pmod->coeff[k + ipos];
	    }
	    F_err = gretl_invert_symmetric_matrix(C);
	    if (!F_err) {
		w = gretl_scalar_qform(b, C, &F_err);
	    }
	    if (F_err) {
		any_F_err = 1;
		var->Fvals[m++] = NADBL;
	    } else {
		var->Fvals[m++] = w / p;
	    }
	    ipos += p;
	}

	/* exclusion of last lag, all vars? */

	if (p > 1) {
	    gretl_matrix_reuse(C, n, n);
	    gretl_matrix_reuse(b, n, 1);

	    ipos = var->ifc + p - 1;
	    for (ii=0; ii<n; ii++) {
		jpos = var->ifc + p - 1;
		for (jj=0; jj<n; jj++) {
		    vij = gretl_matrix_get(V, ipos, jpos);
		    gretl_matrix_set(C, ii, jj, vij);
		    jpos += p;
		}
		b->val[ii] = pmod->coeff[ipos];
		ipos += p;
	    }

	    F_err = gretl_invert_symmetric_matrix(C);
	    if (!F_err) {
		w = gretl_scalar_qform(b, C, &F_err);
	    }
	    if (F_err) {
		any_F_err = 1;
		var->Fvals[m++] = NADBL;
	    } else {
		var->Fvals[m++] = w / n;
	    }
	}

	/* exclusion of all but const? */

	if (var->ifc && var->robust) {
	    gretl_matrix_reuse(C, g-1, g-1);
	    gretl_matrix_reuse(b, g-1, 1);

	    gretl_matrix_extract_matrix(C, V, 1, 1, GRETL_MOD_NONE);
	    for (k=0; k<g-1; k++) {
		b->val[k] = pmod->coeff[k+1];
	    }
	    F_err = gretl_invert_symmetric_matrix(C);
	    if (!F_err) {
		w = gretl_scalar_qform(b, C, &F_err);
	    }
	    if (F_err) {
		any_F_err = 1;
		pmod->fstt = NADBL;
	    } else {
		pmod->fstt = w / (g-1);
	    }
	}
    }

    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    if (!err && any_F_err) {
	fprintf(stderr, "*** Warning: some F-tests could not be computed\n");
    }

    return err;
}

#define VO_DEBUG 0

const int *gretl_VAR_get_exo_list (const GRETL_VAR *var)
{
    return (var == NULL)? NULL : var->xlist;
}

const int *gretl_VAR_get_endo_list (const GRETL_VAR *var)
{
    return (var == NULL)? NULL : var->ylist;
}

/* Based on the specification stored in the original VAR struct, 
   plus a new exogenous list, constitute a full VAR list.
*/

static int *build_VAR_list (const GRETL_VAR *var, int *exolist, int *err)
{
    return VAR_list_composite(var->ylist, exolist, var->rlist);
}

static int gretl_VAR_real_omit_test (const GRETL_VAR *orig,
				     const GRETL_VAR *new,
				     const DATASET *dset,
				     gretlopt opt,
				     PRN *prn)
{
    int *omitlist = NULL;
    double LR, pval;
    int df, nr = 0;
    int err = 0;

#if VO_DEBUG
    fprintf(stderr, "gretl_VAR_real_omit_test: about to diff lists\n");
    printlist(orig->xlist, "orig xlist");
    printlist(new->xlist, "new xlist");
#endif

    if (orig->xlist != NULL) {
	if (new->xlist == NULL) {
	    omitlist = gretl_list_copy(orig->xlist);
	} else {
	    omitlist = gretl_list_diff_new(orig->xlist, new->xlist, 1);
	}
	if (omitlist == NULL) {
	    return E_ALLOC;
	}
	nr = omitlist[0];
    }
    if (opt & OPT_E) {
	/* omitting seasonals */
	nr += dset->pd + 1;
    }
    if (opt & OPT_T) {
	/* omitting trend */
	nr++;
    }    

    LR = orig->T * (new->ldet - orig->ldet);
    df = orig->neqns * nr;
    pval = chisq_cdf_comp(df, LR);

    record_test_result(LR, pval, _("omit"));

    pprintf(prn, "%s:\n", _("Test on the original VAR"));

    print_add_omit_null(omitlist, dset, opt | OPT_S, prn);

    pprintf(prn, "  %s: %s(%d) = %g, ", _("LR test"), 
	    _("Chi-square"), df, LR);
    pprintf(prn, _("with p-value = %g\n"), pval);

    free(omitlist);

    return err;
}

static int VAR_omit_check (GRETL_VAR *var, const int *omitlist,
			   gretlopt opt)
{
    int nx = omitlist == NULL ? 0 : omitlist[0];
    int err = 0;

    if (nx == 0 && !(opt & (OPT_T | OPT_E))) {
	/* nothing to be omitted */
	err = E_NOOMIT;
    } else if ((opt & OPT_T) && !(var->detflags & DET_TREND)) {
	/* can't have the --trend option with no auto-trend */
	err = E_BADOPT;
    } else if ((opt & OPT_E) && !(var->detflags & DET_SEAS)) {
	/* can't have the --seasonals option with no auto-seasonals */
	err = E_BADOPT;
    }    

    return err;
}

/**
 * gretl_VAR_omit_test:
 * @var: pointer to original VAR. 
 * @omitlist: list of variables to omit from original model.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Re-estimates a given VAR after removing the variables
 * specified in @omitlist, and reports per-equation F-tests
 * and system-wide LR tests for the null hypothesis that
 * the omitted variables have zero parameters.
 * 
 * Returns: restricted VAR on success, %NULL on error.
 */

GRETL_VAR *gretl_VAR_omit_test (GRETL_VAR *var, const int *omitlist, 
				DATASET *dset, gretlopt opt,
				PRN *prn, int *err)
{
    GRETL_VAR *vnew = NULL;
    gretlopt varopt = OPT_NONE;
    int smpl_t1 = dset->t1;
    int smpl_t2 = dset->t2;
    int *tmplist = NULL;
    int *varlist = NULL;
    int c1 = 0;

    if (var == NULL) {
	*err = E_DATA;
	return NULL;
    }

    *err = VAR_omit_check(var, omitlist, opt);
    if (*err) {
	return NULL;
    }

    if (var->ifc) {
	c1 = !gretl_list_const_pos(omitlist, 1, dset);
    }

    if (omitlist != NULL && omitlist[0] > 0) {
	/* create reduced exogenous vars list for test VAR */
	tmplist = gretl_list_omit(var->xlist, omitlist, 1, err);
	if (tmplist == NULL) {
	    goto bailout;
	}
    } else if (var->xlist != NULL) {
	tmplist = gretl_list_copy(var->xlist);
	if (tmplist == NULL) {
	    goto bailout;
	}	
    }

    /* create full input VAR list for test VAR */
    varlist = build_VAR_list(var, tmplist, err);
    if (varlist == NULL) {
	goto bailout;
    }

    if ((var->detflags & DET_SEAS) && !(opt & OPT_E)) {
	varopt |= OPT_D;
    }

    if ((var->detflags & DET_TREND) && !(opt & OPT_T)) {
	varopt |= OPT_T;
    }

    /* If the original VAR did not include a constant, we need to
       pass OPT_N to the test VAR to suppress the constant.
       We also need to pass OPT_N in case the constant was
       present originally but is now to be omitted.
    */
    if (var->ifc == 0 || c1 == 0) {
	varopt |= OPT_N;
    }

    /* impose as sample range the estimation range of the 
       original VAR */
    dset->t1 = var->t1;
    dset->t2 = var->t2;

    vnew = gretl_VAR(var->order, var->lags, varlist, dset, 
		     varopt, NULL, err);

    if (vnew != NULL) {
	/* do the actual test(s) */
	*err = gretl_VAR_real_omit_test(var, vnew, dset, opt, prn);
	if (!*err && prn != NULL) {
	    gretl_VAR_print(vnew, dset, OPT_NONE, prn);
	}	
    }

    /* put back into dset what was there on input */
    dset->t1 = smpl_t1;
    dset->t2 = smpl_t2;

 bailout:

    free(tmplist);
    free(varlist);

    return vnew;
}

static int set_wald_R_ones (gretl_matrix *R, int k, int r0, int neq,
			    int stride)
{
    if (k < 0) {
	return E_DATA;
    } else {
	int i, rmax = r0 + neq;

	for (i=r0; i<rmax; i++) {
	    gretl_matrix_set(R, i, k, 1.0);
	    k += stride;
	}
	return 0;
    }
}

/**
 * gretl_VAR_wald_omit_test:
 * @var: pointer to original. 
 * @omitlist: list of variables to omit from original model.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Runs a Wald test for the joint omission of the exogenous
 * variables specified in @omitlist.
 * 
 * Returns: 0 on successful completion, non-zero error code
 * otherwise.
 */

int gretl_VAR_wald_omit_test (GRETL_VAR *var, const int *omitlist, 
			      DATASET *dset, gretlopt opt,
			      PRN *prn)
{
    gretl_matrix *B = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *RB = NULL;
    gretl_matrix *RVR = NULL;
    double x, test = NADBL;
    int i, j, neq, eqk, row0;
    int pos, nr = 0, nx_omit = 0;
    int err = 0;

    if (var == NULL || var->B == NULL || 
	var->S == NULL || var->XTX == NULL) {
	return E_DATA;
    }

    err = VAR_omit_check(var, omitlist, opt);
    if (err) {
	return err;
    }

    if (omitlist != NULL) {
	nx_omit = omitlist[0];
    }

    B = gretl_matrix_vectorize_new(var->B);
    if (B == NULL) {
	return E_ALLOC;
    }

    neq = var->neqns;

    if (nx_omit > 0) {
	/* omitting user-specified exogenous vars */
	nr += nx_omit * neq;
    }
    if (opt & OPT_T) {
	/* omit trend */
	nr += neq;
    }
    if (opt & OPT_E) {
	/* omit seasonals */
	nr += (dset->pd - 1) * neq;
    }

    S = gretl_zero_matrix_new(neq, neq);
    R = gretl_zero_matrix_new(nr, B->rows);
    RB = gretl_matrix_alloc(nr, 1);
    RVR = gretl_matrix_alloc(nr, nr);

    if (S == NULL || R == NULL || RB == NULL || RVR == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<neq; i++) {
	x = gretl_matrix_get(var->S, i, i);
	gretl_matrix_set(S, i, i, x);
    }

    V = gretl_matrix_kronecker_product_new(S, var->XTX, &err);
    if (err) {
	goto bailout;
    }

    eqk = var->B->rows;
    pos = var->ifc + neq * var->order;
    row0 = 0;

    if (nx_omit > 0) {
	/* @omitlist holds ID numbers of exog variables */ 
	int vi;

	for (i=1; i<=omitlist[0]; i++) {
	    vi = omitlist[i];
	    if (vi == 0) {
		j = 0;
	    } else {
		j = in_gretl_list(var->xlist, vi);
		j += pos - 1;
	    }
	    set_wald_R_ones(R, j, row0, neq, eqk);
	    row0 += neq;
	}	
    }
    if (opt & OPT_E) {
	/* omit auto-seasonals */
	for (i=1; i<dset->pd; i++) {
	    j = pos + i - 1;
	    set_wald_R_ones(R, j, row0, neq, eqk);
	    row0 += neq;
	}
    }   
    if (opt & OPT_T) {
	/* omit auto-trend (comes as last param) */
	j = eqk - 1;
	set_wald_R_ones(R, j, row0, neq, eqk);
	row0 += neq;
    }

    gretl_matrix_multiply(R, B, RB);
    gretl_matrix_qform(R, GRETL_MOD_NONE, V,
		       RVR, GRETL_MOD_NONE);
    err = gretl_invert_symmetric_matrix(RVR);

    if (!err) {
	test = gretl_scalar_qform(RB, RVR, &err);
	if (!err) {
	    double pval = chisq_cdf_comp(nr, test);

	    record_test_result(test, pval, _("omit"));

	    if (!(opt & OPT_I)) {
		/* not silent */
		pprintf(prn, "%s:\n", _("Test on VAR"));
		print_add_omit_null(omitlist, dset, opt | OPT_S, prn);
		pprintf(prn, "  %s: %s(%d) = %g, %s %g\n\n",  _("Wald test"),
			_("Chi-square"), nr, test,
			_("p-value"), pval);
	    }
	}
    }

 bailout:

    gretl_matrix_free(B);
    gretl_matrix_free(S);
    gretl_matrix_free(V);
    gretl_matrix_free(R);
    gretl_matrix_free(RB);
    gretl_matrix_free(RVR);

    return err;
}
