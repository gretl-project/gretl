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

/* discrete.c for gretl: logit and probit models */

#include "libgretl.h"
#include "libset.h"

#define TINY 1.0e-13

static int neginv (const double *xpx, double *diag, int nv);
static int cholesky_decomp (double *xpx, int nv);

#define LPDEBUG 0

static double logit (double x)
{
    double l = 1.0 / (1.0 + exp(-x));

#if LPDEBUG
    if (x > 40 || x < -40) {
	fprintf(stderr, "x = %g, logit(x) = %g\n", x, l);
    }
#endif

    return l;
}

static double logit_pdf (double x)
{
    double l, z = exp(-x);

    l = z / ((1.0 + z) * (1.0 + z));

#if LPDEBUG
    if (x > 40 || x < -40) {
	fprintf(stderr, "x = %g, logit_pdf(x) = %g\n", x, l);
    }
#endif

    if (x < 0 && isnan(l)) {
#if LPDEBUG
	fprintf(stderr, "logit_pdf(): x = %g, forcing l to zero\n", x);
#endif
	l = 0;
    }

    return l;
}

static void Lr_chisq (MODEL *pmod, double **Z)
{
    int t, zeros, ones = 0, m = pmod->nobs;
    double Lr, chisq;
    
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}
	if (Z[pmod->list[1]][t] == 1.0) {
	    ones++;
	} 
    }

    zeros = m - ones;

    Lr = (double) ones * log((double) ones / (double) m);
    Lr += (double) zeros * log((double) zeros /(double) m);

    chisq = 2.0 * (pmod->lnL - Lr);
    gretl_model_set_double(pmod, "chisq", chisq);
    
    /* McFadden pseudo-R^2 */
    pmod->rsq = 1.0 - pmod->lnL / Lr;
    pmod->adjrsq = NADBL;
}

static double 
logit_probit_llhood (const double *y, const MODEL *pmod, int opt)
{
    double q, lnL = 0.0;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->yhat[t])) {
	    continue;
	}
	q = 2.0 * y[t] - 1.0;
	if (opt == LOGIT) {
	    lnL += log(logit(q * pmod->yhat[t]));
	} else {
	    lnL += log(normal_cdf(q * pmod->yhat[t]));
	}
    }

    return lnL;
}

static int add_slopes_to_model (MODEL *pmod, double fbx)
{
    double *slopes;
    size_t ssize = pmod->ncoeff * sizeof *slopes;
    int i;

    slopes = malloc(ssize);

    if (slopes == NULL) {
	return 1;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	if (pmod->list[i+2] == 0) {
	    continue;
	}
	slopes[i] = pmod->coeff[i] * fbx;
    }

    if (gretl_model_set_data(pmod, "slopes", slopes, 
			     MODEL_DATA_DOUBLE_ARRAY,
			     ssize)) {
	free(slopes);
	return 1;
    }

    return 0;
}

static double *hess_wts (MODEL *pmod, const double **Z, int ci) 
{
    int t, tw, n = pmod->t2 - pmod->t1 + 1;
    double q, bx, xx;
    double *w;

    w = malloc(n * sizeof *w);
    if (w == NULL) {
	return NULL;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	tw = t - pmod->t1;
	if (model_missing(pmod, t)) {
	    w[tw] = NADBL;
	    continue;
	}

	q = 2.0 * Z[pmod->list[1]][t] - 1.0;
	bx = pmod->yhat[t];

	if (ci == LOGIT) {
	    w[tw] = -1.0 * logit(bx) * (1.0 - logit(bx));
	} else {
	    xx = (q * normal_pdf(q * bx)) / normal_cdf(q * bx);
	    w[tw] = -xx * (xx + bx);
	}
    }

    return w;
}

static double *hessian (MODEL *pmod, const double **Z, int ci) 
{
    int i, j, li, lj, m, t;
    const int l0 = pmod->list[0];
    double xx, *wt, *xpx;

    i = l0 - 1;
    m = i * (i + 1) / 2;

    xpx = malloc(m * sizeof *xpx);
    if (xpx == NULL) {
	return NULL;
    }

    wt = hess_wts(pmod, Z, ci);
    if (wt == NULL) {
	free(xpx);
	return NULL;
    }

    m = 0;
    for (i=2; i<=l0; i++) {
	li = pmod->list[i];
	for (j=i; j<=l0; j++) {
	    lj = pmod->list[j];
	    xx = 0.0;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!model_missing(pmod, t)) {
		    xx += wt[t - pmod->t1] * Z[li][t] * Z[lj][t];
		}
	    }
	    if (floateq(xx, 0.0) && li == lj) {
		free(xpx);
		free(wt);
		return NULL;
	    }
	    xpx[m++] = -xx;
	}
    }

    free(wt);

    return xpx; 
}

static int 
compute_QML_vcv (MODEL *pmod, const double **Z)
{
    gretl_matrix *G = NULL;
    gretl_matrix *H = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;

    const double *y = Z[pmod->list[1]];
    const double *xi;

    double x;
    int k = pmod->ncoeff;
    int T = pmod->nobs;
    int i, j, t, gt, err = 0;

    G = gretl_matrix_alloc(k, T);
    H = gretl_matrix_alloc(k, k);
    S = gretl_matrix_alloc(k, k);
    V = gretl_matrix_alloc(k, k);

    if (G == NULL || H == NULL || S == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* compute gradient or score matrix */
    for (i=0; i<k; i++) {
	xi = Z[pmod->list[i+2]];
	gt = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->yhat[t])) {
		/* score for obs t */
		if (pmod->ci == LOGIT) {
		    x = (y[t] - logit(pmod->yhat[t])) * xi[t];
		} else {
		    double c = normal_cdf(pmod->yhat[t]);

		    x = (y[t] - c) * normal_pdf(pmod->yhat[t]) * xi[t] /
			(c * (1.0 - c));
		}
		gretl_matrix_set(G, i, gt++, x);
	    }
	}
    }

    /* transcribe Hessian from model */
    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    x = pmod->xpx[ijton(i, j, k)];
	    gretl_matrix_set(H, i, j, x);
	    if (i != j) {
		gretl_matrix_set(H, j, i, x);
	    }
	}
    }   

    gretl_invert_symmetric_matrix(H);
    gretl_matrix_multiply_by_scalar(H, -1.0);

    /* form S = GG' */
    gretl_matrix_multiply_mod(G, GRETL_MOD_NONE,
			      G, GRETL_MOD_TRANSPOSE,
			      S, GRETL_MOD_NONE);

    /* form sandwich: V = H^{-1} S H^{-1} */
    gretl_matrix_qform(H, GRETL_MOD_NONE, S,
		       V, GRETL_MOD_NONE);

    pmod->vcv = malloc((k * (k + 1) / 2) * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* set VCV and std err values */
    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    x = gretl_matrix_get(V, i, j);
	    pmod->vcv[ijton(i, j, k)] = x;
	    if (i == j) {
		pmod->sderr[i] = sqrt(x);
	    }
	}
    }

    gretl_model_set_int(pmod, "ml_vcv", VCV_QML);

 bailout:

    gretl_matrix_free(G);
    gretl_matrix_free(H);
    gretl_matrix_free(S);
    gretl_matrix_free(V);

    return err;
}

static MODEL ordered_model (const int *list, int ci, double ***pZ, 
			    DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    MODEL omod;
    void *handle;
    MODEL (* ordered_estimate) (const int *, int, double ***, DATAINFO *, 
				gretlopt, PRN *);

    gretl_error_clear();

    ordered_estimate = get_plugin_function("ordered_estimate", &handle);
    if (ordered_estimate == NULL) {
	gretl_model_init(&omod);
	omod.errcode = E_FOPEN;
	return omod;
    }

    omod = (*ordered_estimate) (list, ci, pZ, pdinfo, opt, prn);

    close_plugin(handle);

    set_model_id(&omod);

    return omod;
}

/* calculate standard errors etc using the Hessian */

static int logit_probit_vcv (MODEL *dmod, gretlopt opt, const double **Z)
{
    int i, err = 0;

    if (dmod->vcv != NULL) {
	free(dmod->vcv);
	dmod->vcv = NULL;
    }

    if (dmod->xpx != NULL) {
	free(dmod->xpx);
    }

    dmod->xpx = hessian(dmod, Z, dmod->ci);
    if (dmod->xpx == NULL) {
	strcpy(gretl_errmsg, _("Failed to construct Hessian matrix"));
	return E_ALLOC;
    } 

    if (opt & OPT_R) {
	err = compute_QML_vcv(dmod, Z);
	if (err) {
	    return err;
	}
    } else {    
	/* obtain negative inverse of Hessian */
	double *xpx = NULL, *diag = NULL;

	cholesky_decomp(dmod->xpx, dmod->ncoeff); 

	diag = malloc(dmod->ncoeff * sizeof *diag); 
	if (diag == NULL) {
	    return E_ALLOC;
	}

	xpx = copyvec(dmod->xpx, dmod->ncoeff * (dmod->ncoeff + 1) / 2);
	if (xpx == NULL) {
	    free(diag);
	    return E_ALLOC;
	}

	neginv(xpx, diag, dmod->ncoeff);

	for (i=0; i<dmod->ncoeff; i++) {
	    dmod->sderr[i] = sqrt(diag[i]);
	}

	free(diag);
	free(xpx);
    }

    return err;
}

#if 0 /* not yet */

static double 
matrix_probit_llhood (const gretl_matrix *y,
		      const gretl_matrix *yhat)
{
    double q, ll = 0.0;
    int t;

    for (t=0; t<y->rows; t++) {
	q = 2.0 * y->val[t] - 1.0;
	ll += log(normal_cdf(q * yhat->val[t]));
    }

    return ll;
}

int gretl_matrix_probit (const gretl_matrix *X,
			 const gretl_matrix *y,
			 gretl_matrix **pb)
{
    gretl_matrix *bX = NULL;
    gretl_matrix *by = NULL;
    gretl_matrix *bi = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *yhat = NULL;
    double ll, lldiff, llbak = -1.0e9;
    double wt, xti, yt, fx, Fx;
    double tol = 1.0e-9;
    int nobs = X->rows;
    int nc = X->cols;
    int itermax = 250;
    int i, t, iter;
    int err = 0;

    bX = gretl_matrix_alloc(nobs, nc);
    by = gretl_column_vector_alloc(nobs);
    bi = gretl_column_vector_alloc(nc);
    yhat = gretl_column_vector_alloc(nobs);
    b = gretl_column_vector_alloc(nc);

    if (bX == NULL || by == NULL || bi == NULL || b == NULL) {
	gretl_matrix_free(bX);
	gretl_matrix_free(by);
	gretl_matrix_free(bi);
	gretl_matrix_free(b);
	return E_ALLOC;
    }

    err = gretl_matrix_ols(y, X, b, NULL, NULL, NULL);

    gretl_matrix_multiply(X, b, yhat);

    for (iter=0; iter<itermax; iter++) {
	/* construct BRMR dataset */
	for (t=0; t<nobs; t++) {
	    yt = yhat->val[t];
	    fx = normal_pdf(yt);
	    Fx = normal_cdf(yt);

	    if (Fx < 1.0) {
		wt = 1.0 / sqrt(Fx * (1.0 - Fx));
	    } else {
		wt = 0.0;
	    }

	    by->val[t] = wt * (y->val[t] - Fx);
	    wt *= fx;
	    for (i=0; i<nc; i++) {
		xti = gretl_matrix_get(X, t, i);
		gretl_matrix_set(bX, t, i, wt * xti);
	    }
	}

	ll = matrix_probit_llhood(y, yhat);

	lldiff = fabs(ll - llbak);
	if (lldiff < tol) {
	    break; 
	}

	llbak = ll;

	err = gretl_matrix_ols(by, bX, bi, NULL, NULL, NULL);
	if (err) {
	    fprintf(stderr, "logit_probit: err = %d\n", err);
	    if (iter > 0) {
		err = E_NOCONV;
	    }
	    break;
	}

	/* update coefficient estimates */
	gretl_matrix_add_to(b, bi);

	/* calculate yhat */
	gretl_matrix_multiply(X, b, yhat);
    }

    if (!err && lldiff > tol) {
	err = E_NOCONV;
    }

    gretl_matrix_free(bX);
    gretl_matrix_free(by);
    gretl_matrix_free(bi);
    gretl_matrix_free(yhat);

    if (!err) {
	*pb = b;
    } else {
	gretl_matrix_free(b);
    }

    return err;
}

#endif

/* BRMR, Davidson and MacKinnon, ETM, p. 461 */

static int do_BRMR (const int *list, MODEL *dmod, int ci,
		    double *beta, const double **Z, 
		    const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *b = NULL;
    const double *yvar = Z[list[1]];
    double lldiff, llbak = -1.0e9;
    double wt, yt, fx, Fx;
    double tol = 1.0e-9;
    int itermax = 250;
    int nc = dmod->ncoeff;
    int i, s, t, iter;
    int err = 0;

    X = gretl_matrix_alloc(dmod->nobs, nc);
    y = gretl_column_vector_alloc(dmod->nobs);
    b = gretl_column_vector_alloc(nc);

    if (X == NULL || y == NULL || b == NULL) {
	gretl_matrix_free(X);
	gretl_matrix_free(y);
	gretl_matrix_free(b);
	return E_ALLOC;
    }

    for (iter=0; iter<itermax; iter++) {
	/* construct BRMR dataset */
	s = 0;
	for (t=dmod->t1; t<=dmod->t2; t++) {
	    yt = dmod->yhat[t];
	    if (!na(yt)) {
		if (ci == LOGIT) {
		    fx = logit_pdf(yt);
		    Fx = logit(yt);
		} else {
		    fx = normal_pdf(yt);
		    Fx = normal_cdf(yt);
		}

		if (Fx < 1.0) {
		    wt = 1.0 / sqrt(Fx * (1.0 - Fx));
		} else {
		    wt = 0.0;
		}

		gretl_vector_set(y, s, wt * (yvar[t] - Fx));
		wt *= fx;
		for (i=0; i<nc; i++) {
		    gretl_matrix_set(X, s, i, wt * Z[list[i+2]][t]);
		}
		s++;
	    }
	}

	if (s != dmod->nobs) {
	    fprintf(stderr, "logit_probit: data error\n");
	    err = E_DATA;
	    break;
	}

	dmod->lnL = logit_probit_llhood(yvar, dmod, ci);

	lldiff = fabs(dmod->lnL - llbak);
	if (lldiff < tol) {
	    break; 
	}

	if (na(dmod->lnL)) {
	    pprintf(prn, _("Iteration %d: log likelihood = NA"), iter);	
	} else {
	    pprintf(prn, _("Iteration %d: log likelihood = %#.12g"), 
		    iter, dmod->lnL);
	}
	pputc(prn, '\n');

	llbak = dmod->lnL;

	err = gretl_matrix_ols(y, X, b, NULL, NULL, NULL);
	if (err) {
	    fprintf(stderr, "logit_probit: err = %d\n", err);
	    if (iter > 0) {
		err = E_NOCONV;
	    }
	    break;
	}

	/* update coefficient estimates */
	for (i=0; i<nc; i++) {
	    beta[i] += b->val[i];
	}

	/* calculate yhat */
	for (t=dmod->t1; t<=dmod->t2; t++) {
	    if (!na(dmod->yhat[t])) {
		dmod->yhat[t] = 0.0;
		for (i=0; i<nc; i++) {
		    dmod->yhat[t] += beta[i] * Z[list[i+2]][t];
		}
	    }
	}
    }

    if (!err && lldiff > tol) {
	err = E_NOCONV;
    }

    if (!err) {
	gretl_model_set_int(dmod, "iters", iter);
	pputc(prn, '\n');
    }

    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(b);

    return err;
}

static MODEL 
binary_logit_probit (const int *list, double ***pZ, DATAINFO *pdinfo, 
		     int ci, gretlopt opt, PRN *prn)
{
    int i, t, depvar = list[1];
    int nx = list[0] - 1;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    int *act_pred = NULL;
    MODEL dmod;
    double xx, fbx, f, F;
    double *xbar = NULL;
    double *beta = NULL;

    gretl_model_init(&dmod);
    
    /* allocate space for means of indep vars, coeffs */
    xbar = malloc(nx * sizeof *xbar);
    beta = malloc(nx * sizeof *beta);
    if (xbar == NULL || beta == NULL) {
	dmod.errcode = E_ALLOC;
	goto bailout;
    }

    /* space for actual/predicted matrix */
    act_pred = malloc(4 * sizeof *act_pred);
    if (act_pred != NULL) {
	for (i=0; i<4; i++) {
	    act_pred[i] = 0;
	}
    }

    varlist_adjust_sample(list, &pdinfo->t1, &pdinfo->t2, 
			  (const double **) *pZ);

    dmod = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (dmod.errcode == 0 && dmod.ncoeff != nx) {
	dmod.errcode = E_DATA;
    }

    if (dmod.errcode) {
	goto bailout;
    }

    for (i=0; i<dmod.ncoeff; i++) {
	beta[i] = dmod.coeff[i];
	xbar[i] = 0.0;
	for (t=dmod.t1; t<=dmod.t2; t++) {
	    if (!model_missing(&dmod, t)) {
		xbar[i] += (*pZ)[list[i+2]][t];
	    }
	}
	xbar[i] /= dmod.nobs;
    }

    dmod.errcode = do_BRMR(list, &dmod, ci, beta, (const double **) *pZ, 
			   pdinfo, prn);
    if (dmod.errcode) {
	goto bailout;
    }

    /* transcribe coefficients */
    for (i=0; i<dmod.ncoeff; i++) {
	dmod.coeff[i] = beta[i];
    }

    dmod.lnL = logit_probit_llhood((*pZ)[depvar], &dmod, ci);
    Lr_chisq(&dmod, *pZ);
    dmod.ci = ci;

    /* calculate standard errors etc */
    dmod.errcode = logit_probit_vcv(&dmod, opt, (const double **) *pZ);
    if (dmod.errcode) {
	goto bailout;
    }

    /* apparatus for calculating slopes at means */
    xx = 0.0;
    for (i=0; i<dmod.ncoeff; i++) {
	xx += dmod.coeff[i] * xbar[i];
    }

    if (ci == LOGIT) {
	fbx = logit_pdf(xx);
#if LPDEBUG
	fprintf(stderr, "xx = %.8g, fbx = %.8g\n", xx, fbx);
#endif
    } else {
	fbx = normal_pdf(xx);
    }

    if (opt & OPT_P) {
	gretl_model_set_int(&dmod, "show-pvals", 1);
    } else if (add_slopes_to_model(&dmod, fbx)) {
	dmod.errcode = E_ALLOC;
	goto bailout;
    }

    /* calculate additional statistics */
    xx = 0.0;
    for (t=dmod.t1; t<=dmod.t2; t++) {
	double xb = dmod.yhat[t];
	int yt;

	if (model_missing(&dmod, t)) {
	    continue;
	}

	yt = (int) (*pZ)[depvar][t];
	xx += yt;

	if (act_pred != NULL) {
	    i = 2 * yt + (xb > 0.0);
	    act_pred[i] += 1;
	}

	if (dmod.ci == LOGIT) {
	    dmod.yhat[t] = exp(xb) / (1.0 + exp(xb)); 
	    dmod.uhat[t] = yt - dmod.yhat[t];
	} else {
	    f = normal_pdf(xb);
	    F = normal_cdf(xb);
	    dmod.yhat[t] = F; 
	    dmod.uhat[t] = yt ? f/F : -f/(1-F);
	}
    }

    xx /= dmod.nobs;
    dmod.ybar = xx;
    dmod.sdy = fbx;

    if (act_pred != NULL) {
	gretl_model_set_data(&dmod, "discrete_act_pred", act_pred, 
			     MODEL_DATA_INT_ARRAY, 
			     4 * sizeof *act_pred);
    }

    mle_criteria(&dmod, 0);

    if (opt & OPT_A) {
	dmod.aux = AUX_AUX;
    } else {
	dmod.ID = model_count_plus();
    }

 bailout:

    free(xbar);
    free(beta);

    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    return dmod;
}

static int 
ordered_model_ok (double **Z, const DATAINFO *pdinfo, int v)
{
    if (!var_is_discrete(pdinfo, v)) {
	sprintf(gretl_errmsg, "The variable '%s' is not discrete",
		pdinfo->varname[v]);
	return 0;
    } else {
	double m = gretl_min(pdinfo->t1, pdinfo->t2, Z[v]);

	if (m != 0.0) {
	    sprintf(gretl_errmsg, "The minimum value of '%s' is not 0",
		    pdinfo->varname[v]);
	    return 0;
	}
    }

    return 1;
}

/**
 * logit_probit:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: command index: if = %LOGIT, perform logit regression, otherwise
 * perform probit regression.
 * @opt: if includes %OPT_R form robust (QML) estimates of standard
 * errors and covariance matrix, in binary case; if %OPT_P arrange for
 * printing of p-values, not slopes at mean; if %OPT_A treat as an
 * auxiliary regression.
 * @prn: printing struct in case additional information is
 * wanted (%OPT_V).
 *
 * Computes estimates of the discrete model specified by @list,
 * using an estimator determined by the value of @ci.  In the
 * binary case, uses the BRMR auxiliary regression; see Davidson 
 * and MacKinnon.  If the dependent variable is not binary but 
 * is discrete and has a minimum value of 0, we do ordered 
 * logit/probit.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL logit_probit (const int *list, double ***pZ, DATAINFO *pdinfo, 
		    int ci, gretlopt opt, PRN *prn)
{
    int yv = list[1];

    if (gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[yv])) {
	return binary_logit_probit(list, pZ, pdinfo, ci, opt, prn);
    } else if (ordered_model_ok(*pZ, pdinfo, yv)) {
	return ordered_model(list, ci, pZ, pdinfo, opt, prn);
    } else {
	MODEL dmod;

	gretl_model_init(&dmod);
	dmod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, _("The dependent variable '%s' is not a 0/1 "
				"variable.\n"), pdinfo->varname[yv]);
	return dmod;
    }
}

/* Solves for diagonal elements of X'X inverse matrix.
   X'X must be Cholesky-decomposed already.
*/

static int neginv (const double *xpx, double *diag, int n)
{
    int kk, l, m, k, i, j;
    const int nxpx = n * (n + 1) / 2;
    double *tmp;
    double d, e;

    tmp = malloc((n + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<=n; i++) {
	tmp[i] = 0.0;
    }

    kk = 0;

    for (l=1; l<=n-1; l++) {
        d = xpx[kk];
        tmp[l] = d;
        e = d * d;
        m = 0;
        if (l > 1) {
	    for (j=1; j<=l-1; j++) {
		m += n - j;
	    }
	}
        for (i=l+1; i<=n; i++) {
            d = 0.0;
            k = i + m - 1;
            for (j=l; j<=i-1; j++) {
                d += tmp[j] * xpx[k];
                k += n - j;
            }
            d = -d * xpx[k];
            tmp[i] = d;
            e += d * d;
        }
        kk += n + 1 - l;
        diag[l-1] = e;
    }

    diag[n - 1] = xpx[nxpx - 1] * xpx[nxpx - 1];

    free(tmp);

    return 0;
}

/* Cholesky decomposition of X'X */

static int cholesky_decomp (double *xpx, int nv)
{
    int i, j, k, kk, l, jm1;
    double e, d, d1, test, xx;

    e = 1.0 / sqrt(xpx[0]);
    xpx[0] = e;

    for (i=1; i<nv; i++) {
	xpx[i] *= e;
    }

    kk = nv;

    for (j=2; j<=nv; j++) {
	/* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d += xx * xx;
            k += nv-l;
        }
        test = xpx[kk] - d;
        if (test / xpx[kk] < TINY) {
	    return 1;
	}
        e = 1 / sqrt(test);
        xpx[kk] = e;
        /* off-diagonal elements */
        for (i=j+1; i<=nv; i++) {
            kk++;
            d = 0.0;
            k = j - 1;
            for (l=1; l<=jm1; l++) {
                d += xpx[k] * xpx[k-j+i];
                k += nv - l;
            }
            xpx[kk] = (xpx[kk] - d) * e;
        }
        kk++;
    }

    return 0; 
}

/* logistic model: doesn't exactly belong here but it seems like
   the most suitable place for it 
*/

int logistic_ymax_lmax (const double *y, const DATAINFO *pdinfo,
			double *ymax, double *lmax)
{
    int t;

    *ymax = 0.0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(y[t])) {
	    continue;
	}
	if (y[t] <= 0.0) {
	    strcpy(gretl_errmsg, _("Illegal non-positive value of the "
				   "dependent variable"));
	    return 1;
	}
	if (y[t] > *ymax) {
	    *ymax = y[t];
	}
    }

    if (*ymax < 1.0) {
	*lmax = 1.0;
    } else if (*ymax < 100.0) {
	*lmax = 100.0;
    } else {
	*lmax = 1.1 * *ymax;
    }
	    
    return 0;
}

static double real_get_lmax (const double *y, const DATAINFO *pdinfo,
			     const char *lmstr)
{
    double lmax, ymax = 0.0;
    int err;

    err = logistic_ymax_lmax(y, pdinfo, &ymax, &lmax);

    if (err) {
	return NADBL;
    }

    if (lmstr != NULL && *lmstr != '\0') {
	lmax = atof(lmstr + 5);
	if (lmax <= ymax) {
	    gretl_errmsg_set(_("Invalid value for the maximum of the "
			       "dependent variable"));
	    lmax = NADBL;
	}
    }	
	    
    return lmax;
}

static int make_logistic_depvar (double ***pZ, DATAINFO *pdinfo, 
				 int dv, double lmax)
{
    int t, v = pdinfo->v;

    if (dataset_add_series(1, pZ, pdinfo)) {
	return 1;
    }

    for (t=0; t<pdinfo->n; t++) {
	double p = (*pZ)[dv][t];

	if (na(p)) {
	    (*pZ)[v][t] = NADBL;
	} else {
	    (*pZ)[v][t] = log(p / (lmax - p));
	}
    }

    return 0;
}

static int rewrite_logistic_stats (const double **Z, const DATAINFO *pdinfo,
				   MODEL *pmod, int dv, double lmax)
{
    int t;
    double x;

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, Z[dv]);

    /* make the VCV matrix before messing with the model stats */
    makevcv(pmod, pmod->sigma);

    pmod->ess = 0.0;
    for (t=0; t<pdinfo->n; t++) {
	x = pmod->yhat[t];
	if (na(x)) {
	    continue;
	}
	pmod->yhat[t] = lmax / (1.0 + exp(-x));
	pmod->uhat[t] = Z[dv][t] - pmod->yhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	x = Z[dv][t];
	if (!na(x)) {
	    pmod->tss += (x - pmod->ybar) * (x - pmod->ybar);
	}
    }

    pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    pmod->list[1] = dv;
    gretl_model_set_double(pmod, "lmax", lmax);
    pmod->ci = LOGISTIC;
    ls_criteria(pmod);

    return 0;
}

/**
 * logistic_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @param: may contain "ymax=value" for user setting of the
 * asymptotic maximum of the dependent variable.
 *
 * Estimate the model given in @list using the logistic transformation
 * of the dependent variable.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL logistic_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		      const char *param) 
{
    double lmax;
    int *llist = NULL;
    int dv = list[1];
    MODEL lmod;

    gretl_model_init(&lmod); 

    llist = gretl_list_copy(list);
    if (llist == NULL) {
	lmod.errcode = E_ALLOC;
	return lmod;
    }

    lmax = real_get_lmax((*pZ)[dv], pdinfo, param);

    if (na(lmax)) {
	lmod.errcode = E_DATA;
    } else if (lmax == 0.0) {
	lmod.errcode = E_CANCEL;
    }

    if (!lmod.errcode) {
	if (make_logistic_depvar(pZ, pdinfo, dv, lmax)) {
	    lmod.errcode = E_ALLOC;	
	}
    }

    if (lmod.errcode) {
	free(llist);
	return lmod;
    }

    llist[1] = pdinfo->v - 1;

    lmod = lsq(llist, pZ, pdinfo, OLS, OPT_A);
    if (!lmod.errcode) {
	rewrite_logistic_stats((const double **) *pZ, pdinfo, &lmod,
			       dv, lmax);
	set_model_id(&lmod);
    }

    dataset_drop_last_variables(1, pZ, pdinfo);
    free(llist);
    
    return lmod;
}

static double table_prob (double a, double b, double c, double d,
			  double num, double nf)
{
    double den, P;

    den = nf * x_factorial(a) * x_factorial(b) * 
	x_factorial(c) * x_factorial(d);
    P = num / den;

#if 0
    fprintf(stderr, "\n%6g %6g\n", a, b);
    fprintf(stderr, "%6g %6g\n", c, d);
    fprintf(stderr, " Pi = %g\n", P);
#endif

    return P;
}

/* for use with cross-tabulation (xtab command) in the
   2 x 2 case
*/

int fishers_exact_test (const Xtab *tab, PRN *prn)
{
    double a, b, c, d, E0;
    double P0, Pi, PL, PR, P2;
    double num, nf;

    a = tab->f[0][0];
    b = tab->f[0][1];
    c = tab->f[1][0];
    d = tab->f[1][1];

    E0 = (tab->rtotal[0] * tab->ctotal[0]) / (double) tab->n;

    num = x_factorial(a + b) * x_factorial(c + d) * x_factorial(a + c) * 
	x_factorial(b + d);
    nf = x_factorial(tab->n);
    
    /* Probability of the observed table */
    PL = PR = P2 = P0 = table_prob(a, b, c, d, num, nf);

    while (a > 0 && d > 0) {
	a -= 1; d -= 1;
	c += 1; b += 1;
	Pi = table_prob(a, b, c, d, num, nf);
	if (Pi <= P0 || tab->f[0][0] > E0) {
	    PL += Pi;
	}
	if (Pi <= P0) {
	    P2 += Pi;
	}
    }

    a = tab->f[0][0];
    b = tab->f[0][1];
    c = tab->f[1][0];
    d = tab->f[1][1];

    while (c > 0 && b > 0) {
	c -= 1; b -= 1;
	a += 1; d += 1;
	Pi = table_prob(a, b, c, d, num, nf);
	if (Pi <= P0 || tab->f[0][0] < E0) {
	    PR += Pi;
	}
	if (Pi <= P0) {
	    P2 += Pi;
	}
    } 

    pprintf(prn, "\n%s:\n", _("Fisher's Exact Test"));
    pprintf(prn, "  Left:   P-value = %g\n", PL);
    pprintf(prn, "  Right:  P-value = %g\n", PR);
    pprintf(prn, "  2-Tail: P-value = %g\n", P2);
    pputc(prn, '\n');

    return 0;
}
