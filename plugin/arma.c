/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* Native gretl code for ARMA estimation.  Much of the code here
   was contributed by Riccardo "Jack" Lucchetti, the rest is due
   to Allin Cottrell; thanks also to Stephen Moshier for cephes.
*/

#include "libgretl.h"
#include "internal.h"

#include "bhhh_max.h"

#include "../cephes/polrt.c"

#undef ARMA_DEBUG

static void add_arma_varnames (MODEL *pmod, const DATAINFO *pdinfo)
{
    int p = pmod->list[1];
    int q = pmod->list[2];
    int r = pmod->list[0] - 4;
    int i, j, np = 2 + p + q + r;

    pmod->params = malloc(np * sizeof pmod->params);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    pmod->nparams = np;

    for (i=0; i<np; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    pmod->errcode = E_ALLOC;
	    return;
	}
    }

    strcpy(pmod->params[0], pdinfo->varname[pmod->list[4]]);
    strcpy(pmod->params[1], pdinfo->varname[0]);

    j = 2;
    for (i=0; i<p; i++) {
	const char *depvar = pmod->params[0];
	size_t n = strlen(depvar);
	
	if (n < VNAMELEN - 4) {
	    sprintf(pmod->params[j++], "%s(-%d)", depvar, i + 1);
	} else {
	    sprintf(pmod->params[j++], "y(-%d)", i + 1);
	}
    }
    for (i=0; i<q; i++) {
	sprintf(pmod->params[j++], "e(-%d)", i + 1);
    }
    for (i=0; i<r; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[pmod->list[5+i]]);
    }    
}

static int ma_out_of_bounds (int q, const double *ma_coeff)
{
    double *temp = NULL, *tmp2 = NULL;
    double re, im, rt;
    cmplx *roots = NULL;
    int i, err = 0;

    temp  = malloc((q+1) * sizeof *temp);
    tmp2  = malloc((q+1) * sizeof *tmp2);
    roots = malloc(q * sizeof *roots);

    if (temp == NULL || tmp2 == NULL || roots == NULL) {
	free(temp);
	free(tmp2);
	free(roots);
	return 1;
    }

    temp[0] = 1.0;
    /* ma_coeff is 1-based */
    for (i=1; i<=q; i++){
	temp[i] = ma_coeff[i];
    }
    polrt(temp, tmp2, q, roots);

    for (i=0; i<q; i++) {
	re = roots[i].r;
	im = roots[i].i;
	rt = re * re + im * im;
	if (rt > DBL_EPSILON && rt <= 1.0) {
	    fprintf(stderr, "MA root %d = %g\n", i, rt);
	    err = 1;
	    break;
	}
    }

    free(temp);
    free(tmp2);
    free(roots);

    return err;
}

static int arma_ll (double *coeff, 
		    const double **X, double **Z, 
		    model_info *arma,
		    int do_score)
{
    int i, j, t;
    int t1 = model_info_get_t1(arma);
    int t2 = model_info_get_t2(arma);
    int n = t2 - t1 + 1;
    int p, q, r;
    const double K = 1.41893853320467274178; /* ln(sqrt(2*pi)) + 0.5 */
    const double *y = X[0];
    double **series = model_info_get_series(arma);
    double *e = series[0];
    double **de = series + 1;
    const double *ar_coeff, *ma_coeff, *reg_coeff;
    double ll, s2 = 0.0;
    int err = 0;

    model_info_get_pqr(arma, &p, &q, &r);
    ar_coeff = coeff;
    ma_coeff = coeff + p;
    reg_coeff = coeff + p + q;

    if (ma_out_of_bounds(q, ma_coeff)) {
	fputs("MA estimate(s) out of bounds\n", stderr);
	return 1;
    }

    /* update forecast errors */

    for (t=t1; t<=t2; t++) {

	e[t] = y[t] - coeff[0];

	for (i=1; i<=p; i++) {
	    e[t] -= ar_coeff[i] * y[t-i];
	}

	for (i=1; i<=q; i++) {
	    if (t - i >= t1) {
		e[t] -= ma_coeff[i] * e[t-i];
	    }
	}

	for (i=1; i<=r; i++) {
	    e[t] -= reg_coeff[i] * X[i][t];
	}

	s2 += e[t] * e[t];
    }

    /* get error variance and log-likelihood */

    s2 /= (double) n;

    ll = -n * (0.5 * log(s2) + K);
    model_info_set_ll(arma, ll, do_score);

    if (do_score) {
	int col, nc = p + q + r + 1;
	double x;

	for (t=t1; t<=t2; t++) {

	    col = 0;

	    /* the constant term (de_c) */
	    de[col][t] = -1.0;
	    for (i=1; i<=q; i++) {
		de[col][t] -= ma_coeff[i] * de[col][t-i];
	    }

	    /* AR terms (de_a) */
	    for (j=0; j<p; j++) {
		col++;
		if (t >= col) {
		    de[col][t] = -y[t-col];
		    for (i=1; i<=q; i++) {
			de[col][t] -= ma_coeff[i] * de[col][t-i];
		    }
		}
	    }

	    /* MA terms (de_m) */
	    for (j=0; j<q; j++) {
		int lag = ++col - p;

		if (t >= lag) {
		    de[col][t] = -e[t-lag];
		    for (i=1; i<=q; i++) {
			de[col][t] -= ma_coeff[i] * de[col][t-i];
		    }
		}
	    }

	    /* ordinary regressors */
	    for (j=0; j<r; j++) {
		col++;
		de[col][t] = -X[j+1][t]; 
		for (i=1; i<=q; i++) {
		    de[col][t] -= ma_coeff[i] * de[col][t-i];
		}
	    }

	    /* update OPG data set */
	    x = e[t] / s2; /* sqrt(s2)? */
	    for (i=0; i<nc; i++) {
		Z[i+1][t] = -de[i][t] * x;
	    }
	}
    }

    if (isnan(ll)) err = 1;

    return err;
}

static cmplx *arma_roots (int p, int q, const double *coeff) 
     /*
       Given an ARMA process $A(L) y_t = C(L) \epsilon_t$, returns the 
       roots of the two polynomials;

       Syntax:
       p: order of A(L)
       q: order of C(L)
       coeff: p+q+1 vector of coefficients (element 0 is the constant
       and is ignored)
       returns: the p + q roots (AR part first)
     */
{
    int maxlag = (p > q)? p : q;
    const double *ar_coeff = coeff;
    const double *ma_coeff = coeff + p;
    double *temp = NULL, *temp2 = NULL;
    cmplx *roots = NULL;
    int j;

    temp  = malloc((maxlag+1) * sizeof *temp);
    temp2 = malloc((maxlag+1) * sizeof *temp2);
    roots = malloc((p + q) * sizeof *roots);

    if (temp == NULL || temp2 == NULL || roots == NULL) {
	free(temp);
	free(temp2);
	free(roots);
	return NULL;
    }

    temp[0] = 1.0;

    /* A(L) */
    for (j=1; j<=p; j++){
	temp[j] = -ar_coeff[j];
    }
    polrt(temp, temp2, p, roots);

    /* C(L) */
    for (j=1; j<=q; j++){
	temp[j] = ma_coeff[j];
    }
    polrt(temp, temp2, q, roots + p);

    free(temp);
    free(temp2);

    return roots;
}

static void rewrite_arma_model_stats (MODEL *pmod, model_info *arma,
				      const int *list, const double *y, 
				      const double *theta, int nc,
				      const DATAINFO *pdinfo)
{
    int i, t;
    int p = list[1], q = list[2], r = list[0] - 4;
    double **series = model_info_get_series(arma);
    const double *e = series[0];
    double mean_error;

    pmod->ci = ARMA;
    pmod->ifc = 1;

    pmod->lnL = model_info_get_ll(arma);

    pmod->dfn = p + q + r;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = nc;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = theta[i];
    }

    copylist(&pmod->list, list);

    pmod->ybar = _esl_mean(pmod->t1, pmod->t2, y);
    pmod->sdy = _esl_stddev(pmod->t1, pmod->t2, y);

    mean_error = pmod->ess = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = e[t];
	pmod->yhat[t] = y[t] - pmod->uhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	mean_error += pmod->uhat[t];
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->tss += (y[t] - pmod->ybar) * (y[t] - pmod->ybar);
    }

    if (pmod->tss > pmod->ess) {
	pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / 
	    (pmod->dfn * pmod->ess);
    } else {
	pmod->fstt = NADBL;
    }

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    /* AIC, as per X-12-ARIMA */
    pmod->criterion[C_AIC] = -2.0 * pmod->lnL + 2.0 * (pmod->ncoeff + 1);

    /* BIC, as per X-12-ARIMA */
    pmod->criterion[C_BIC] = -2.0 * pmod->lnL + (pmod->ncoeff + 1) * log(pmod->nobs);
}

static void remove_const (int *list)
{
    int i, j;

    for (i=5; i<=list[0]; i++) {
	if (list[i] == 0) {
	    for (j=i; j<list[0]; j++) {
		list[j] = list[j+1];
	    }
	    list[0] -= 1;
	    break;
	}
    }
}

static int check_arma_list (int *list)
{
    int err = 0;

    if (list[0] < 4) err = 1;

    /* for now we'll accept ARMA (4,4) at max */
    else if (list[1] < 0 || list[1] > 4) err = 1;
    else if (list[2] < 0 || list[2] > 4) err = 1;
    else if (list[1] + list[2] == 0) err = 1;

    /* remove const from list of regressors (ARMAX), since it
       is added automatically */
    if (list[0] > 4) {
	remove_const(list);
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    }
    
    return err;
}

static const double **make_armax_X (int *list, const double **Z)
{
    const double **X;
    int nv = list[0] - 4;
    int v, i;

    X = malloc((nv + 1) * sizeof *X);
    if (X == NULL) return NULL;

    X[0] = Z[list[4]];

    for (i=1; i<=nv; i++) {
	v = list[i + 4];
	X[i] = Z[v];
    }

    return X;
}

/* below: run an initial OLS to get starting values for the
   AR coefficients */

static int ar_init_by_ols (const int *list, double *coeff,
			   const double **Z, const DATAINFO *pdinfo,
			   int arma_t1)
{
    int an = pdinfo->t2 - arma_t1 + 1;
    int p = list[1];
    int q = list[2];
    int r = list[0] - 4;
    int ynum = list[4];
    int av = p + r + 2;
    double **aZ = NULL;
    DATAINFO *ainfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int i, j, t, err = 0;

    gretl_model_init(&armod, NULL);  

    alist = malloc((av + 1) * sizeof *alist);
    if (alist == NULL) return 1;

    alist[0] = av;
    alist[1] = 1;
    alist[2] = 0;
    for (i=0; i<p; i++) {
	alist[i + 3] = i + 2;
    }
    for (i=0; i<r; i++) {
	alist[i + 3 + p] = i + p + 2;
    }

    ainfo = create_new_dataset(&aZ, av, an, 0);
    if (ainfo == NULL) {
	free(alist);
	return 1;
    }
    
    for (t=0; t<an; t++) {
	int j, s;

	for (i=0; i<=p; i++) {
	    s = t + arma_t1 - i;
	    aZ[i+1][t] = Z[ynum][s];
	}
	for (i=0; i<r; i++) {
	    j = list[i+5];
	    s = t + arma_t1;
	    aZ[i+p+2][t] = Z[j][s];
	}
    }

    armod = lsq(alist, &aZ, ainfo, OLS, OPT_A, 0.0);
    err = armod.errcode;
    if (!err) {
	j = 0;
	for (i=0; i<armod.ncoeff; i++) {
	    if (i == p + 1) j += q; /* leave space for MA coeffs */
	    coeff[j++] = armod.coeff[i];
	}
	for (i=0; i<q; i++) {
	    /* squeeze in some zeros */
	    coeff[i+p+1] = 0.0;
	} 
    }

    free(alist);
    free_Z(aZ, ainfo);
    clear_datainfo(ainfo, CLEAR_FULL);
    free(ainfo);

    clear_model(&armod, NULL);

    return err;
}

static int adjust_sample (DATAINFO *pdinfo, const double **Z, const int *list,
			  int *arma_t1, int *arma_t2)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int p = list[1];
    int q = list[2];
    int r = list[0] - 4;
    int maxlag = (p > q)? p : q;
    int an, i, v, t, t1min = 0;
    int anymiss;

    for (t=0; t<=pdinfo->t2; t++) {
	anymiss = 0;
	for (i=4; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) t1min++;
        else break;
    }

    t1min += maxlag;
    if (t1 < t1min) t1 = t1min;

    for (t=pdinfo->t2; t>=t1; t--) {
	anymiss = 0;
	for (i=4; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) t2--;
        else break;
    }

    for (t=t1-p; t<t2; t++) {
	for (i=4; i<=list[0]; i++) {
	    if (t < t1 && i > 4) continue;
	    v = list[i];
	    if (na(Z[v][t])) {
		char msg[64];

		sprintf(msg, _("Missing value encountered for "
                           "variable %d, obs %d"), v, t + 1);
		gretl_errmsg_set(msg);
		return 1;
	    }
	}
    }

    an = t2 - t1 + 1;
    if (an <= p + q + r + 1) return 1; 

    *arma_t1 = t1;
    *arma_t2 = t2;

    return 0;
}

static model_info *
set_up_arma_info (int t1, int t2, int p, int q, int r)
{
    model_info *arma = model_info_new();

    if (arma == NULL) return NULL;

    model_info_set_opts(arma, PRESERVE_OPG_MODEL);
    model_info_set_tol(arma, 1.0e-6);
    model_info_set_pqr(arma, p, q, r);
    model_info_set_n_series(arma, 2 + p + q + r);
    model_info_set_t1_t2(arma, t1, t2);

    return arma;
}

MODEL arma_model (int *list, const double **Z, DATAINFO *pdinfo, 
		  PRN *prn)
{
    int nc, v, p, q, r, maxlag;
    int arma_t1, arma_t2;
    int err = 0;
    double *coeff;
    const double **X = NULL;
    MODEL armod;
    model_info *arma;

    gretl_model_init(&armod, NULL);  

    if (check_arma_list(list)) {
	armod.errcode = E_UNSPEC;
	return armod;
    }

    v = list[4]; /* dependent var */
    p = list[1]; /* AR order */
    q = list[2]; /* MA order */
    maxlag = (p > q)? p : q; /* maximum lag in the model */

    r = list[0] - 4; /* number of ordinary regressors */

    /* adjust sample? */
    if (adjust_sample(pdinfo, Z, list, &arma_t1, &arma_t2)) {
        armod.errcode = E_DATA;
        return armod;
    }

    /* number of coefficients */
    nc = 1 + p + q + r;

    /* coefficient vector (plus error variance) */
    coeff = malloc(nc * sizeof *coeff);
    if (coeff == NULL) {
	armod.errcode = E_ALLOC;
	return armod;
    }

    arma = set_up_arma_info(arma_t1, arma_t2, p, q, r);
    if (arma == NULL) {
	free(coeff);
	armod.errcode = E_ALLOC;
	return armod;
    }

    /* initialize the coefficients: AR and regression part by OLS, 
       MA at 0 */
    err = ar_init_by_ols(list, coeff, Z, pdinfo, arma_t1);
    if (err) {
	free(coeff);
	model_info_free(arma);
	armod.errcode = E_ALLOC;
	return armod;
    }	

    /* construct virtual dataset for dep var, real regressors */
    X = make_armax_X(list, Z);
    if (X == NULL) {
	armod.errcode = E_ALLOC;
	free(coeff);
	return armod;
    }

    err = bhhh_max(arma_ll, X, coeff, arma, prn);

    if (err) {
	armod.errcode = E_NOCONV;
    } else {
	MODEL *pmod = model_info_capture_OPG_model(arma);
	double *theta = model_info_get_theta(arma);
	cmplx *roots;

	rewrite_arma_model_stats(pmod, arma, list, Z[v], theta, nc,
				 pdinfo);

	/* compute and save polynomial roots */
	roots = arma_roots(p, q, theta);
	if (roots != NULL) {
	    gretl_model_set_data(pmod, "roots", roots,
				 (p + q) * sizeof *roots);
	}

	add_arma_varnames(pmod, pdinfo);

	armod = *pmod;
    }

    free(coeff);
    free(X);
    model_info_free(arma);

    return armod;
}
