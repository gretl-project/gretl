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
#include "missing_private.h"
#include "gretl_bfgs.h"

#include <errno.h>

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

    if (chisq < 0) {
	pmod->rsq = pmod->adjrsq = pmod->chisq = NADBL;
    } else {
	pmod->chisq = chisq;
	/* McFadden pseudo-R^2 */
	pmod->rsq = 1.0 - pmod->lnL / Lr;
	pmod->adjrsq = 1.0 - (pmod->lnL - pmod->ncoeff) / Lr;
    }
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
			     GRETL_TYPE_DOUBLE_ARRAY,
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

    err = gretl_model_write_vcv(pmod, V);

    if (!err) {
	gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
	pmod->opt |= OPT_R;
    }

 bailout:

    gretl_matrix_free(G);
    gretl_matrix_free(H);
    gretl_matrix_free(S);
    gretl_matrix_free(V);

    return err;
}

static MODEL ordered_model (int *list, int ci, double ***pZ, 
			    DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    MODEL omod;
    void *handle;
    MODEL (* ordered_estimate) (int *, int, double ***, DATAINFO *, 
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

/*
  If min1 > max0, then there exists a separating hyperplane between
  all the zeros and all the ones; in this case, the likelihood has
  no maximum (despite having a supremum at 1) and no MLE exists.
*/

static int perfect_pred_check (const double *y, MODEL *dmod)
{
    double max0 = -1.0e200;
    double min1 = 1.0e200;
    double yht;
    int t;

    for (t=dmod->t1; t<=dmod->t2; t++) {
	yht = dmod->yhat[t];
	if (!na(yht)) {
	    if (y[t] == 0 && yht > max0) {
		max0 = yht;
	    } 
	    if (y[t] == 1 && yht < min1) {
		min1 = yht;
	    }
	}
    }

    return (min1 > max0);
}

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

	if (perfect_pred_check(yvar, dmod)) {
	    gretl_errmsg_sprintf("Perfect prediction detected at iteration %d;\nno MLE exists",
				 iter);
	    err = E_NOCONV;
	    break;
	}

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

#if LPDEBUG > 1
	gretl_matrix_print(y, "y, in BRMR");
	gretl_matrix_print(X, "X, in BRMR");
	gretl_matrix_print(b, "b, from BRMR");
#endif

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

static char *classifier_check (int *list, double **Z, DATAINFO *pdinfo,
			       gretlopt opt, PRN *prn, int *err)
{
    char *mask = NULL;
    const double *y = Z[list[1]];
    int i, v, ni, t;

    for (i=2; i<=list[0]; i++) {
	v = list[i];
	if (v == 0) {
	    continue;
	}
	ni = gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[v]);
	if (ni > 0) {
	    int same[2] = {0};
	    int diff[2] = {0};
	    int maskval, pc = 1;

	    for (t=pdinfo->t1; t<=pdinfo->t2 && pc; t++) {
		if (Z[v][t] > 0) {
		    if (y[t] > 0) {
			same[0] += 1;
		    } else {
			diff[0] += 1;
		    }
		} else {
		    if (y[t] == 0) {
			same[1] += 1;
		    } else {
			diff[1] += 1;
		    }
		}
		if (same[0] && diff[0] && same[1] && diff[1]) {
		    pc = 0;
		}
	    }

	    if (pc) {
		if (!(same[0] && diff[0])) {
		    maskval = 1;
		    pprintf(prn, "Note: %s != 0 predicts %s perfectly\n",
			    pdinfo->varname[v], (same[0])? "success" : "failure");
		} else {
		    maskval = 0;
		    pprintf(prn, "Note: %s == 0 predicts %s perfectly\n",
			    pdinfo->varname[v], (diff[1])? "success" : "failure");
		}

		if (mask == NULL) {
		    mask = malloc(pdinfo->n + 1);
		    if (mask == NULL) {
			*err = E_ALLOC;
		    } else {
			memset(mask, '0', pdinfo->n);
			mask[pdinfo->n] = 0;
		    }
		}

		if (mask != NULL) {
		    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
			if (Z[v][t] == maskval) {
			    mask[t-pdinfo->t1] = '1';
			}
		    }
		    pprintf(prn, "%d observations not used\n", ni);
		}
		gretl_list_delete_at_pos(list, i--);
	    }
	}
    }

    return mask;
}

static MODEL 
binary_logit_probit (const int *inlist, double ***pZ, DATAINFO *pdinfo, 
		     int ci, gretlopt opt, PRN *prn)
{
    int *list = NULL;
    char *mask = NULL;
    int i, t, depvar, nx;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    int *act_pred = NULL;
    MODEL dmod;
    double xx, fbx, f, F;
    double *xbar = NULL;
    double *beta = NULL;

    /* FIXME do we need to insist on a constant in this sort
       of model? */

    gretl_model_init(&dmod);

    list = gretl_list_copy(inlist);
    if (list == NULL) {
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

    mask = classifier_check(list, *pZ, pdinfo, opt, prn,
			    &dmod.errcode);
    if (dmod.errcode) {
	goto bailout;
    }

    depvar = list[1];
    nx = list[0] - 1;

    /* allocate space for means of indep vars, coeffs */
    xbar = malloc(nx * sizeof *xbar);
    beta = malloc(nx * sizeof *beta);
    if (xbar == NULL || beta == NULL) {
	dmod.errcode = E_ALLOC;
	goto bailout;
    }

    if (mask != NULL) {
	dmod.errcode = copy_to_reference_missmask(mask);
    }

    if (!dmod.errcode) {
	dmod = lsq(list, pZ, pdinfo, OLS, OPT_A);
	if (dmod.errcode == 0 && dmod.ncoeff != nx) {
	    dmod.errcode = E_DATA;
	}
    }

    if (dmod.errcode) {
	goto bailout;
    }

#if LPDEBUG
    printmodel(&dmod, pdinfo, OPT_NONE, prn);
#endif

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

    dmod.errcode = do_BRMR(list, &dmod, ci, beta, 
			   (const double **) *pZ, pdinfo, 
			   (opt & OPT_V)? prn : NULL);
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
	dmod.opt |= OPT_P;
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
			     GRETL_TYPE_INT_ARRAY, 
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
    free(list);
    free(mask);

    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    return dmod;
}

/* struct for holding multinomial logit info */

typedef struct mnl_info_ mnl_info;

struct mnl_info_ {
    gretl_matrix_block *B;
    gretl_matrix *y;  /* dependent variable */
    gretl_matrix *X;  /* regressors */
    gretl_matrix *b;  /* coefficients, matrix form */
    gretl_matrix *Xb; /* coeffs times regressors */
    gretl_matrix *P;  /* probabilities */
    double *theta;    /* coeffs for BFGS */
    int n;            /* number of categories (excluding base) */
    int k;            /* number of coeffs per category */
    int npar;         /* total number of parameters */
    int T;            /* number of observations */
};

static void mnl_info_destroy (mnl_info *mnl)
{
    if (mnl != NULL) {
	gretl_matrix_block_destroy(mnl->B);
	free(mnl->theta);
	free(mnl);
    }
}

static mnl_info *mnl_info_new (int n, int k, int T)
{
    mnl_info *mnl = malloc(sizeof *mnl);
    int i;

    if (mnl != NULL) {
	mnl->npar = k * n;
	mnl->theta = malloc(mnl->npar * sizeof *mnl->theta);
	if (mnl->theta == NULL) {
	    free(mnl);
	    return NULL;
	}
	mnl->B = gretl_matrix_block_new(&mnl->y, T, 1,
					&mnl->X, T, k,
					&mnl->b, k, n,
					&mnl->Xb, T, n,
					&mnl->P, T, n,
					NULL);
	if (mnl->B == NULL) {
	    free(mnl->theta);
	    free(mnl);
	    return NULL;
	} else {
	    for (i=0; i<mnl->npar; i++) {
		mnl->theta[i] = 0.0;
	    }
	    mnl->n = n;
	    mnl->k = k;
	    mnl->T = T;
	}
    }

    return mnl;
}

/* compute loglikelihood for multinomial logit */

static double mn_logit_loglik (const double *theta, void *ptr)
{
    mnl_info *mnl = (mnl_info *) ptr;
    double x, xti, ll = 0.0;
    int npar = mnl->npar;
    int yt, i, t;

    errno = 0;

    for (i=0; i<npar; i++) {
	mnl->b->val[i] = theta[i];
    }

    gretl_matrix_multiply(mnl->X, mnl->b, mnl->Xb);

    for (t=0; t<mnl->T && !errno; t++) {
	x = 1.0;
	for (i=0; i<mnl->n; i++) {
	    /* sum row i of exp(Xb) */
	    xti = exp(gretl_matrix_get(mnl->Xb, t, i));
	    x += xti;
	}
	ll -= log(x);
	yt = gretl_vector_get(mnl->y, t);
	if (yt > 0) {
	    ll += gretl_matrix_get(mnl->Xb, t, yt-1);
	}
    }

    if (errno != 0) {
	ll = NADBL;
    }

    return ll;
}

static int mn_logit_score (double *theta, double *s, int npar, BFGS_CRIT_FUNC ll, 
			   void *ptr)
{
    mnl_info *mnl = (mnl_info *) ptr;
    double x, xti, pti, p, g;
    int i, j, k, t, yt;

    errno = 0;

    for (i=0; i<npar; i++) {
	s[i] = 0.0;
    }

    for (t=0; t<mnl->T && !errno; t++) {
	x = 1.0;
	for (i=0; i<mnl->n; i++) {
	    /* sum row i of exp(Xb) */
	    xti = exp(gretl_matrix_get(mnl->Xb, t, i));
	    x += xti;
	    gretl_matrix_set(mnl->P, t, i, xti);
	}

	if (!errno) {
	    yt = gretl_vector_get(mnl->y, t);
	    k = 0;
	    for (i=0; i<mnl->n; i++) {
		pti = gretl_matrix_get(mnl->P, t, i) / x;
		gretl_matrix_set(mnl->P, t, i, pti);
		p = (i == (yt-1)) - pti;
		for (j=0; j<mnl->k; j++) {
		    g = p * gretl_matrix_get(mnl->X, t, j);
		    s[k++] += g;
		}
	    }
	}
    }

    return 1;
}

/* multinomial logit; form the negative inverse of the analytical
   Hessian */

static gretl_matrix *mnl_hess_inverse (mnl_info *mnl, int *err)
{
    gretl_matrix_block *B;
    gretl_matrix *H;
    gretl_matrix *x;
    gretl_matrix *xx;
    gretl_matrix *hjk;
    double xti, ptj, ptk;
    int r, c;
    int i, j, k, t;

    H = gretl_zero_matrix_new(mnl->npar, mnl->npar);
    if (H == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    B = gretl_matrix_block_new(&x, 1, mnl->k,
			       &xx, mnl->k, mnl->k,
			       &hjk, mnl->k, mnl->k,
			       NULL);
    if (B == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(H);
	return NULL;
    }

    r = c = 0;

    for (j=0; j<mnl->n; j++) {
	for (k=0; k<=j; k++) {
	    gretl_matrix_zero(hjk);
	    for (t=0; t<mnl->T; t++) {
		for (i=0; i<mnl->k; i++) {
		    xti = gretl_matrix_get(mnl->X, t, i);
		    gretl_vector_set(x, i, xti);
		}
		gretl_matrix_multiply_mod(x, GRETL_MOD_TRANSPOSE,
					  x, GRETL_MOD_NONE,
					  xx, GRETL_MOD_NONE);
		ptj = gretl_matrix_get(mnl->P, t, j);
		ptk = gretl_matrix_get(mnl->P, t, k);
		gretl_matrix_multiply_by_scalar(xx, ptj * ((j == k) - ptk));
		gretl_matrix_add_to(hjk, xx);
	    }
	    gretl_matrix_inscribe_matrix(H, hjk, r, c, GRETL_MOD_NONE);
	    if (j != k) {
		gretl_matrix_inscribe_matrix(H, hjk, c, r, GRETL_MOD_NONE);
	    }
	    c += mnl->k;
	}
	r += mnl->k;
	c = 0;
    }

    gretl_matrix_block_destroy(B);

    *err = gretl_invert_symmetric_matrix(H);

    if (*err) {
	gretl_matrix_free(H);
	H = NULL;
    }

    return H;
}

static gretl_matrix *mnl_opg_matrix (mnl_info *mnl, int *err)
{
    gretl_matrix *G;
    gretl_matrix *GG;
    double p, g;
    int yt, i, j, k, t;

    G = gretl_matrix_alloc(mnl->T, mnl->npar);
    if (G == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    GG = gretl_matrix_alloc(mnl->npar, mnl->npar);
    if (GG == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(G);
	return NULL;
    }

    for (t=0; t<mnl->T; t++) {
	yt = gretl_vector_get(mnl->y, t);
	k = 0;
	for (i=0; i<mnl->n; i++) {
	    p = (i == (yt-1)) - gretl_matrix_get(mnl->P, t, i);
	    for (j=0; j<mnl->k; j++) {
		g = p * gretl_matrix_get(mnl->X, t, j);
		gretl_matrix_set(G, t, k++, g);
	    }
	}
    }

    gretl_matrix_multiply_mod(G, GRETL_MOD_TRANSPOSE,
			      G, GRETL_MOD_NONE,
			      GG, GRETL_MOD_NONE);

    gretl_matrix_free(G);

    return GG;
}

static int mnl_add_variance_matrix (MODEL *pmod, mnl_info *mnl,
				    gretlopt opt)
{
    gretl_matrix *H = NULL;
    gretl_matrix *G = NULL;
    gretl_matrix *V = NULL;
    int err = 0;

    H = mnl_hess_inverse(mnl, &err);
    if (err) {
	return err;
    }

    if (opt & OPT_R) {
	G = mnl_opg_matrix(mnl, &err);
	if (!err) {
	    V = gretl_matrix_alloc(mnl->npar, mnl->npar);
	    if (V == NULL) {
		err = E_ALLOC;
	    } else {
		err = gretl_matrix_qform(H, GRETL_MOD_NONE, G,
					 V, GRETL_MOD_NONE);
	    }
	}
    }

    if (!err) {
	if (opt & OPT_R) {
	    err = gretl_model_write_vcv(pmod, V);
	    if (!err) {
		gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
		pmod->opt |= OPT_R;
	    }
	} else {
	    err = gretl_model_write_vcv(pmod, H);
	    if (!err) {
		gretl_model_set_vcv_info(pmod, VCV_ML, VCV_HESSIAN);
	    }
	}
    }

    gretl_matrix_free(H);
    gretl_matrix_free(G);
    gretl_matrix_free(V);

    return err;
}

/* Construct 'yhat' and 'uhat'.  Maybe this is too simple-minded; we
   just find, for each observation, the y-value for which the
   probability is maximized and set that as yhat[t].  We then set
   uhat[t] as a binary "hit" (residual = 0) or "miss" (residual = 1).
   Constructing a "quantitative" residual as y[t] - yhat[t] seems
   spurious, since there's no meaningful metric for the "distance"
   between y and yhat when y is an unordered response.
*/

static void mn_logit_yhat (MODEL *pmod, mnl_info *mnl,
			   const gretl_matrix *yvals)
{
    double p, pmax;
    int i, s, t, iymax;
    int ncorrect = 0;

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->yhat[t])) {
	    continue;
	}
	pmax = 1.0;
	iymax = 0;
	for (i=0; i<mnl->n; i++) {
	    p = gretl_matrix_get(mnl->Xb, s, i);
	    if (p > pmax) {
		pmax = p;
		iymax = i + 1;
	    }
	}
	if (iymax == (int) mnl->y->val[s]) {
	    ncorrect++;
	    pmod->uhat[t] = 0;
	} else {
	    pmod->uhat[t] = 1;
	}
	if (yvals != NULL) {
	    pmod->yhat[t] = yvals->val[iymax];
	} else {
	    pmod->yhat[t] = iymax;
	}
	s++;
    }

    gretl_model_set_int(pmod, "correct", ncorrect);
}

struct sorter {
    double x;
    int t;
};

/* In case the dependent variable is not in canonical form for
   multinomial logit, construct a transformed version */

static int make_canonical_depvar (MODEL *pmod, const double *y, 
				  gretl_matrix *yvec)
{
    struct sorter *s;
    double nexty, bady;
    int i, t, n = pmod->nobs;

    s = malloc(n * sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    s[i].x = y[t];
	    s[i].t = i;
	    i++;
	}
    }

    qsort(s, n, sizeof *s, gretl_compare_doubles);

    /* normalize to a minimum of zero */
    if (s[0].x != 0) {
	double ymin = s[0].x;

	for (i=0; i<n && s[i].x == ymin; i++) {
	    s[i].x = 0.0;
	}
    }

    /* ensure that the sorted values increase by steps of one */
    for (i=1; i<n; i++) {
	if (s[i].x != s[i-1].x) {
	    nexty = s[i-1].x + 1;
	    if (s[i].x != nexty) {
		bady = s[i].x;
		while (i < n && s[i].x == bady) {
		    s[i++].x = nexty;
		}
		i--; /* compensate for outer i++ */
	    } 
	} 
    }

    /* write canonical version of y into yvec */
    for (i=0; i<n; i++) {
	yvec->val[s[i].t] = s[i].x;
    }

    free(s);

    return 0;
}

/* multinomial logit: count the distinct values taken on by the
   dependent variable.  If the variable does not take the form of a
   sequence of consecutive integers with a base of zero, flag this
   by returning via @yvals the vector of distinct values.
*/

static int mn_value_count (const double *y, MODEL *pmod, gretl_vector **yvals)
{
    double *sy = malloc(pmod->nobs * sizeof *sy);
    int T = pmod->nobs;
    int n = 0;

    if (sy == NULL) {
	pmod->errcode = E_ALLOC;
    } else {
	gretl_vector *v;
	int t, s = 0;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t])) {
		sy[s++] = y[t];
	    }
	}

	v = gretl_matrix_values(sy, T, &pmod->errcode);
	if (pmod->errcode) {
	    free(sy);
	    return 0;
	}

	if (v->val[0] != 0.0) {
	    /* we'll need v later */
	    *yvals = v;
	}

	n = gretl_vector_get_length(v);

	for (t=0; t<n && !pmod->errcode; t++) {
	    if (v->val[t] != floor(v->val[t])) {
		pmod->errcode = E_DATA;
		gretl_errmsg_set("logit: the dependent variable must form a sequence of "
				 "integers");
	    } else if (t > 0 && v->val[t] != v->val[t-1] + 1) {
		/* we'll need v later */
		*yvals = v;
	    }
	}
	if (pmod->errcode) {
	    *yvals = NULL;
	    n = 0;
	} else if (*yvals == v) {
	    /* avoid freeing v */
	    v = NULL;
	}
	gretl_matrix_free(v);
	free(sy);
    }

    return n;
}

const char *mn_logit_coeffsep (const MODEL *pmod, const DATAINFO *pdinfo, int i)
{
    static char ret[32];
    const char *vname = gretl_model_get_depvar_name(pmod, pdinfo);
    const int *list = (const int *) gretl_model_get_data(pmod, "yvals");
    int val = (list != NULL)? list[i+1] : i;

    sprintf(ret, "%s = %d", vname, val);
    return ret;
}

/* transcribe multinomial logit results into @pmod, which was
   initialized via OLS, and add covariance matrix
*/

static void mnl_finish (mnl_info *mnl, MODEL *pmod,
			const gretl_matrix *yvals,
			const DATAINFO *pdinfo,
			gretlopt opt)
{
    int i;

    pmod->errcode = gretl_model_write_coeffs(pmod, mnl->theta, mnl->npar);

    if (!pmod->errcode) {
	pmod->errcode = mnl_add_variance_matrix(pmod, mnl, opt);
    }

    if (!pmod->errcode) {
	pmod->ci = LOGIT;
	pmod->lnL = mn_logit_loglik(mnl->theta, mnl);
	mle_criteria(pmod, 0);

	gretl_model_set_int(pmod, "multinom", mnl->n);
	gretl_model_set_int(pmod, "cblock", mnl->k);

	pmod->ess = pmod->sigma = NADBL;
	pmod->fstt = pmod->rsq = pmod->adjrsq = NADBL;

	gretl_model_allocate_params(pmod, mnl->npar);
    }

    if (!pmod->errcode) {
	int j, vj, k = 0;

	for (i=0; i<mnl->n; i++) {
	    for (j=0; j<mnl->k; j++) {
		vj = pmod->list[j+2];
		strcpy(pmod->params[k++], pdinfo->varname[vj]);
	    }
	}
	mn_logit_yhat(pmod, mnl, yvals);
	set_model_id(pmod);
    }

    if (!pmod->errcode && yvals != NULL) {
	int *list = gretl_list_new(yvals->rows);

	if (list != NULL) {
	    for (i=0; i<yvals->rows; i++) {
		list[i+1] = yvals->val[i];
	    }
	    gretl_model_set_list_as_data(pmod, "yvals", list);
	}
    }

    if (!pmod->errcode) {
	/* add overall likelihood ratio test */
	int t, ni, df = pmod->ncoeff;
	double LR, L0 = 0.0;

	if (pmod->ifc) {
	    df -= mnl->n;
	}

	for (i=0; i<=mnl->n; i++) {
	    ni = 0;
	    for (t=0; t<mnl->T; t++) {
		if (mnl->y->val[t] == i) {
		    ni++;
		}
	    }
	    if (ni > 0) {
		if (pmod->ifc) {
		    L0 += ni * log((double) ni / mnl->T);
		} else {
		    L0 += ni * log(1.0 / (mnl->n + 1));
		}
	    }
	}

	LR = 2.0 * (pmod->lnL - L0);
	gretl_model_set_double(pmod, "wald", LR);
	gretl_model_set_int(pmod, "wald_df", df);
    }

    if (pmod->xpx != NULL) {
	free(pmod->xpx);
	pmod->xpx = NULL;
    }

    pmod->opt |= OPT_M;
}

/* multinomial logit */

static MODEL mnl_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn)
{
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;
    MODEL mod;
    mnl_info *mnl;
    gretl_vector *yvals = NULL;
    int n, k = list[0] - 1;
    int i, vi, t, s;

    /* we'll start with OLS to flush out data issues */
    mod = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (mod.errcode) {
	return mod;
    }

    n = mn_value_count((*pZ)[list[1]], &mod, &yvals);
    if (mod.errcode) {
	return mod;
    }

    n--; /* exclude the first value */

    if (n * k > mod.nobs) {
	mod.errcode = E_DF;
	return mod;
    }

    mnl = mnl_info_new(n, k, mod.nobs);
    if (mnl == NULL) {
	mod.errcode = E_ALLOC;
	return mod;
    }

    if (yvals != NULL) {
	/* the dependent variable needs transforming */
	mod.errcode = make_canonical_depvar(&mod, (*pZ)[list[1]], mnl->y);
	if (mod.errcode) {
	    goto bailout;
	}
    } else {
	s = 0;
	for (t=mod.t1; t<=mod.t2; t++) {
	    if (!na(mod.yhat[t])) {
		mnl->y->val[s++] = (*pZ)[list[1]][t];
	    }
	}
    }

    for (i=0; i<k; i++) {
	vi = list[i+2];
	s = 0;
	for (t=mod.t1; t<=mod.t2; t++) {
	    if (!na(mod.yhat[t])) {
		gretl_matrix_set(mnl->X, s++, i, (*pZ)[vi][t]);
	    }
	}
    } 

    mod.errcode = BFGS_max(mnl->theta, mnl->npar, maxit, 0.0, 
			   &fncount, &grcount, mn_logit_loglik, C_LOGLIK,
			   mn_logit_score, mnl, 
			   (prn != NULL)? OPT_V : OPT_NONE, prn);

    if (!mod.errcode) {
	mnl_finish(mnl, &mod, yvals, pdinfo, opt);
    }  

 bailout:  

    mnl_info_destroy(mnl);
    gretl_matrix_free(yvals);

    return mod;
}

static int 
discrete_depvar_ok (double **Z, const DATAINFO *pdinfo, int v)
{
    if (!var_is_discrete(pdinfo, v) && 
	!gretl_is_oprobit_ok(pdinfo->t1, pdinfo->t2, Z[v])) {
	sprintf(gretl_errmsg, _("The variable '%s' is not discrete"),
		pdinfo->varname[v]);
	return 0;
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

MODEL logit_probit (int *list, double ***pZ, DATAINFO *pdinfo, 
		    int ci, gretlopt opt, PRN *prn)
{
    PRN *vprn = NULL;
    int yv = list[1];

    if (ci == LOGIT && (opt & OPT_M)) {
	vprn = (opt & OPT_V)? prn : NULL;
	return mnl_model(list, pZ, pdinfo, opt, vprn);
    } else if (gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[yv])) {
	return binary_logit_probit(list, pZ, pdinfo, ci, opt, prn);
    } else if (discrete_depvar_ok(*pZ, pdinfo, yv)) {
	vprn = (opt & OPT_V)? prn : NULL;
	return ordered_model(list, ci, pZ, pdinfo, opt, vprn);
    } else {
	MODEL dmod;

	gretl_model_init(&dmod);
	dmod.errcode = E_UNSPEC;
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
    double x, ess, sigma;
    int t;

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, Z[dv]);

    /* make the VCV matrix before messing with the model stats */
    makevcv(pmod, pmod->sigma);

    ess = 0.0;

    for (t=0; t<pdinfo->n; t++) {
	x = pmod->yhat[t];
	if (na(x)) {
	    continue;
	}
	pmod->yhat[t] = lmax / (1.0 + exp(-x));
	pmod->uhat[t] = Z[dv][t] - pmod->yhat[t];
	ess += pmod->uhat[t] * pmod->uhat[t];
    }

    sigma = sqrt(ess / pmod->dfd);

    pmod->list[1] = dv;
    gretl_model_set_double(pmod, "lmax", lmax);
    gretl_model_set_double(pmod, "ess_orig", ess);
    gretl_model_set_double(pmod, "sigma_orig", sigma);
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

static double choose (double n, double k)
{
    double c = 1.0;
    int i;

    for (i=0; i<k; i++) {
	c *= (n - i) / (k - i);
    }

    return c;
}

static double table_prob (double a, double b, double c, double d,
			  double n)
{
    double p1 = choose(a+b, a);
    double p2 = choose(c+d, c);
    double p3 = choose(n, a+c);
    double P = p1 * p2 / p3;

#if 0
    fprintf(stderr, "\ntable_prob: a=%g, b=%g, c=%g, d=%g, n=%g\n", 
	    a, b, c, d, n);
    fprintf(stderr, " p1=%g, p2=%g, p3=%g; P = %g\n", 
	    p1, p2, p3, P);
#endif

    return P;
}

/* for use with cross-tabulation (xtab command) in the
   2 x 2 case
*/

int fishers_exact_test (const Xtab *tab, PRN *prn)
{
    double a, b, c, d, n, E0;
    double P0, Pi, PL, PR, P2;

    a = tab->f[0][0];
    b = tab->f[0][1];
    c = tab->f[1][0];
    d = tab->f[1][1];
    n = tab->n;

    E0 = (tab->rtotal[0] * tab->ctotal[0]) / n;

    /* Probability of the observed table */
    PL = PR = P2 = P0 = table_prob(a, b, c, d, n);

    while (a > 0 && d > 0) {
	a -= 1; d -= 1;
	c += 1; b += 1;
	Pi = table_prob(a, b, c, d, n);
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
	Pi = table_prob(a, b, c, d, n);
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

/**
 * intreg:
 * @list: high/low (2 variables) plus list of regressors.
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
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL intreg (int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt, 
	      PRN *prn)
{
    MODEL intmod;
    void *handle;
    MODEL (* interval_estimate) (int *, double ***, DATAINFO *, gretlopt, 
				 PRN *);

    gretl_error_clear();

    interval_estimate = get_plugin_function("interval_estimate", &handle);
    if (interval_estimate == NULL) {
	gretl_model_init(&intmod);
	intmod.errcode = E_FOPEN;
	return intmod;
    }

    intmod = (*interval_estimate) (list, pZ, pdinfo, opt, prn);

    close_plugin(handle);

    set_model_id(&intmod);

    return intmod;
}

