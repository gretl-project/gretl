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

#include "libgretl.h"
#include "internal.h"
#include "../cephes/polrt.c"

#define DEFAULT_MAX_ITER 1000

static void arma_coeff_name (char *s, const DATAINFO *pdinfo,
			     const MODEL *pmod, int i)
{
    int j, p = pmod->list[1];

    if (i == 0) {
	strcpy(s, pdinfo->varname[pmod->list[4]]);
	return;
    }

    if (i == 1 && pmod->ifc) {
	strcpy(s, pdinfo->varname[0]);
	return;
    }

    if (pmod->ifc) j = i - 1;
    else j = i;

    if (j - p < 1) {
	const char *depvar = pmod->params[0];
	size_t n = strlen(depvar);
	
	if (n < VNAMELEN - 4) {
	    sprintf(s, "%s(-%d)", depvar, j);
	} else {
	    sprintf(s, "y(-%d)", j);
	}
    } else {
	sprintf(s, "e(-%d)", j - p);
    }
}

static void add_arma_varnames (MODEL *pmod, const DATAINFO *pdinfo)
{
    int i, np = 2 + pmod->list[1] + pmod->list[2];

    pmod->params = malloc(np * sizeof pmod->params);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    for (i=0; i<np; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->errcode = E_ALLOC;
	    return;
	}
    }

    for (i=0; i<np; i++) { 
	arma_coeff_name(pmod->params[i], pdinfo, pmod, i);
    }
}

static int get_maxiter (void)
{
    static int ml = 0;

    if (ml == 0) {
        char *mlstr = getenv("GRETL_MAX_ITER");

        if (mlstr != NULL && sscanf(mlstr, "%d", &ml)) ;
        else ml = DEFAULT_MAX_ITER;
    }

    return ml;
}

static double update_fcast_errs (double *e, const double *y, 
				 const double *coeff, 
				 int n, int p, int q)
{
    int i, k, t;
    double s2 = 0.0;
    int SampleSize = 0;
    int maxlag = (p>q) ? p : q;

    for (t=0; t<n; t++) {

	if (na(y[t])){
	  e[t] = NADBL;
	  continue;
	}

	e[t] = y[t] - coeff[0];

	for (i=0; i<p; i++) {
	    k = t - (i + 1);
	    if (k < 0 || na(y[k])) continue;
	    e[t] -= coeff[i+1] * y[k];
	}

	for (i=0; i<q; i++) {
	    k = t - (i + 1);
	    if (k < 0) continue;
	    e[t] -= coeff[i+p+1] * e[k];
	}

	if(t>=maxlag){
	  s2 += e[t] * e[t];
	  SampleSize++;
	}
    }

    s2 /= (double) SampleSize;

    return s2;
}

static void update_error_partials (double *e, const double *y, 
				   double **Z, const double *coeff,
				   int n, int p, int q)
{
    int t, col, i, t2;

    for (t=0; t<n; t++) {
	if (na(y[t])) continue;

	/* the constant term */
	col = 2;
	Z[col][t] = -1.0;
	for (i=0; i<q; i++) {
	    if (t > i) Z[col][t] -= coeff[i+p+1] * Z[col][t-i-1];
	}

	/* AR terms */
	if (p > 0) {
	    for (col=3; col <= p+3; col++) {
		if (t>=col-2) {
		    Z[col][t] = -y[t-col+2];
		    for (i=0; i<q; i++) {
			if (t > i) Z[col][t] -= coeff[i+p+1] * Z[col][t-i-1];
		    }
		}
	    }
	}

	/* MA terms */
	if (q > 0){
	    for (col=p+3; col <= p+q+3; col++) {
		t2 = col - p - 2;
		if (t>=t2) {
		    Z[col][t] = -e[t-t2];
		    for (i=0; i<q; i++) {
			if (t > i) Z[col][t] -= coeff[i+p+1] * Z[col][t-i-1];
		    }
		}
	    }
	}
    }
}

static void update_ll_partials (double *e, double **Z, double s2,
				int n, int p, int q)
{
    int i, j, t;
    int maxlag = (p>q) ? p : q;
    int nc = p + q + 1;

    for (i=0; i<nc; i++) {
	j = 3 + p + q + i;
	for (t=0; t<maxlag; t++) {
	    Z[j][t] = NADBL;
	}
	for (t=maxlag; t<n; t++) {
	    if (na(e[t])) continue;
	    Z[j][t] = -Z[j-nc][t] * e[t] / s2;
	}
    }
}

static double get_ll (const double *e, int n, int maxlag)
{
    int t;
    double ll = 0.0, s2;
    double SampleSize = n - maxlag;
    const double K = 1.41893853320467274178;
    /* that's ln(sqrt(2*pi)) + 0.5 */

    for (t=maxlag; t<n; t++) {
	if (na(e[t])) continue;
	ll += e[t] * e[t];
    }

    s2 = ll / SampleSize;
    ll = -SampleSize * (0.5 * log(s2) + K);

    return ll;
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
    int maxlag = (p>q) ? p : q;

    double *temp, *temp2;
    cmplx  *roots;
    int j;

    temp =  malloc((maxlag+1) * sizeof *temp);
    temp2 = malloc((maxlag+1) * sizeof *temp2);
    roots = malloc((p + q) * sizeof *roots);

    if (temp == NULL || temp2 == NULL || roots == NULL) {
	return NULL;
    }

    temp[0] = 1.0;

    /* A(L) */
    for (j=1; j<p+1; j++){
	temp[j] = -coeff[j];
    }

    polrt(temp, temp2, p, roots);

    /* C(L) */
    for (j=1; j<q+1; j++){
	temp[j] = coeff[j+p];
    }

    polrt(temp, temp2, q, roots + p);

    free(temp);
    free(temp2);

    return roots;
}

static int reallocate_hat_series (MODEL *pmod, int n)
{
    double *tmp;

    tmp = realloc(pmod->uhat, n * sizeof *pmod->uhat);
    if (tmp != NULL) {
	pmod->uhat = tmp;
    } else {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    tmp = realloc(pmod->yhat, n * sizeof *pmod->yhat);
    if (tmp != NULL) {
	pmod->yhat = tmp;
    } else {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    return 0;
}

static void rewrite_arma_model_stats (MODEL *pmod, const double *coeff,
				      const int *list, const double *y, 
				      const double *e, 
				      const DATAINFO *pdinfo)
{
    int i, t;
    int p = list[1], q = list[2];
    int realt1 = pmod->t1 + pdinfo->t1;
    int realt2 = pmod->t2 + pdinfo->t1;
    double mean_error;

    pmod->ci = ARMA;
    pmod->ifc = 1;

    pmod->dfn = p + q;
    pmod->dfd = pmod->nobs - pmod->dfn;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = coeff[i];
    }

    copylist(&pmod->list, list);

    pmod->ybar = _esl_mean(realt1, realt2, y);
    pmod->sdy = _esl_stddev(realt1, realt2, y);

    /* if model was estimated on a sub-sample, we need to expand the
       residual and fitted series (which are always full-length) */
    if (pdinfo->t1 + (pdinfo->n - 1 - pdinfo->t2) > 0) {
	if (reallocate_hat_series(pmod, pdinfo->n)) {
	    return;
	}
    }

    mean_error = pmod->ess = 0.0;
    for (t=0; t<pdinfo->n; t++) {
	if (t < realt1 || t > realt2) {
	    pmod->uhat[t] = pmod->yhat[t] = NADBL;
	} else {
	    pmod->uhat[t] = e[t - realt1];
	    pmod->yhat[t] = y[t] - pmod->uhat[t];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    mean_error += pmod->uhat[t];
	}
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=realt1; t<=realt2; t++) {
	if (!na(y[t])) {
	    pmod->tss += (y[t] - pmod->ybar) * (y[t] - pmod->ybar);
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

    /* get the model start/end dates in sync with the main dataset */
    pmod->t1 = realt1;
    pmod->t2 += pdinfo->t1;

    /* AIC, as per X-12-ARIMA */
    pmod->criterion[0] = -2.0 * pmod->lnL + 2.0 * (pmod->ncoeff + 1);
    /* BIC, as per X-12-ARIMA */
    pmod->criterion[1] = -2.0 * pmod->lnL + (pmod->ncoeff + 1) * log(pmod->nobs);
}

static int check_arma_list (const int *list)
{
    int err = 0;

    if (list[0] != 4) err = 1;

    /* for now we'll accept ARMA (4,4) at max */
    else if (list[1] < 0 || list[1] > 4) err = 1;
    else if (list[2] < 0 || list[2] > 4) err = 1;
    else if (list[1] + list[2] == 0) err = 1;

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    }
    
    return err;
}

/* below: run an initial OLS to get starting values for the
   coefficients, in the pure AR case.  But the pure AR case
   doesn't seem to work right anyway... FIXME */

static int ar_init_by_ols (int v, int p, double *coeff,
			   const double **Z, const DATAINFO *pdinfo)
{
    int an = pdinfo->t2 - pdinfo->t1 + 1;
    double **aZ = NULL;
    DATAINFO *ainfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int i, t, err = 0;

    gretl_model_init(&armod, pdinfo);  

    alist = malloc((p + 3) * sizeof *alist);
    if (alist == NULL) return 1;

    alist[0] = p + 2;
    alist[1] = 1;
    alist[2] = 0;
    for (i=0; i<p; i++) {
	alist[i + 3] = i + 2;
    }

    ainfo = create_new_dataset(&aZ, p + 2, an, 0);
    if (ainfo == NULL) {
	free(alist);
	return 1;
    }
    ainfo->extra = 1; 
    
    for (t=0; t<an; t++) {
	for (i=0; i<=p; i++) {
	    int k = t + pdinfo->t1 - i;

	    if (k < 0) aZ[i+1][t] = NADBL;
	    else aZ[i+1][t] = Z[v][k];
	}
    }

    armod = lsq(alist, &aZ, ainfo, OLS, 1, 0.0);
    err = armod.errcode;
    if (!err) {
	for (i=0; i<armod.ncoeff; i++) {
	    coeff[i] = armod.coeff[i];
	}
    }
    
    free(alist);
    free_Z(aZ, ainfo);
    clear_datainfo(ainfo, CLEAR_FULL);
    free(ainfo);

    clear_model(&armod, NULL);

    return err;
}

MODEL arma_model (int *list, const double **Z, DATAINFO *pdinfo, 
		  PRN *prn)
{
    int an = pdinfo->t2 - pdinfo->t1 + 1;
    int nc, v, p, q, maxlag;
    int iters, itermax;
    int subiters, subitermax;
    int i, t;
    double s2, tol, crit, ll = 0.0;
    double steplength, ll_prev = -1.0e+8;
    double *e, *coeff, *d_coef;
    const double *y;
    double **aZ;
    DATAINFO *ainfo;
    int *alist;
    MODEL armod;

    gretl_model_init(&armod, pdinfo);  

    if (check_arma_list(list)) {
	armod.errcode = E_UNSPEC;
	return armod;
    }

    v = list[4]; /* the dependent variable */
    p = list[1]; /* AR order */
    q = list[2]; /* MA order */
    maxlag = (p>q) ? p : q; /* maximum lag in the model */

    /* number of coefficients */
    nc = 1 + p + q;

    alist = malloc((nc + 2) * sizeof *alist);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	return armod;
    }

    alist[0] = nc + 1;
    alist[1] = 0; /* dep var is constant, in OPG */
    for (i=0; i<=p+q; i++) {
	alist[i+2] = 3 + p + q + i;
    }

    /* coefficient vector */
    coeff = malloc(nc * sizeof *coeff);
    if (coeff == NULL) {
	free(alist);
	armod.errcode = E_ALLOC;
	return armod;
    }

    /* and its change on each iteration */
    d_coef = malloc(nc * sizeof *d_coef);
    if (d_coef == NULL) {
	free(alist);
	free(coeff);
	armod.errcode = E_ALLOC;
	return armod;
    }

    /* temporary dataset: we need two series for each coefficient,
       plus one for the forecast errors and one for the const */
    ainfo = create_new_dataset(&aZ, nc * 2 + 2, an, 0);
    if (ainfo == NULL) {
	free(alist);
	free(coeff);
	free(d_coef);
	armod.errcode = E_ALLOC;
	return armod;
    }
    ainfo->extra = 1; /* ? can't remember what this does */

    /* initialize the coefficients */
    /* AR part by OLS, MA at 0's */
    ar_init_by_ols(v, p, coeff, Z, pdinfo);
    for (i=0; i<q; i++) coeff[i+p+1] = 0.0; 

    /* initialize forecast errors and derivatives */
    for (t=0; t<an; t++) {
	for (i=1; i<=nc+1; i++) {
	    aZ[i][t] = 0.0;
	}
    }

    /* forecast errors */
    e = aZ[1];
    /* dependent variable */
    y = &Z[v][pdinfo->t1];

    crit = 1.0;
    tol = 1.0e-6;
    iters = 0;
    itermax = get_maxiter();
    subitermax = 20;
    steplength = 0.125;

    /* generate one-step forecast errors */
    s2 = update_fcast_errs(e, y, coeff, an, p, q);

    /* calculate log-likelihood */
    ll = get_ll(e, an, maxlag);

    while (crit > tol && iters++ < itermax && !isnan(ll)) {

	pprintf(prn, "Iteration %d\n", iters);
        pprintf(prn, "  log likelihood = %g\n", ll);

        /* partials of e wrt coeffs */
	update_error_partials(e, y, aZ, coeff, an, p, q);

	/* partials of l wrt coeffs */
	update_ll_partials(e, aZ, s2, an, p, q);

	/* OPG regression */
	clear_model(&armod, NULL);
	armod = lsq(alist, &aZ, ainfo, OLS, 1, 0.0);
	if (armod.errcode) {
	    goto arma_bailout;
	}

	/* compute direction */
	for (i=0; i<armod.ncoeff; i++) {
	    d_coef[i] = armod.coeff[i];
	}

	/* update parameters */
	for (i=0; i<armod.ncoeff; i++) {
	  coeff[i] += d_coef[i] * steplength;
	}

	/* compute log-likelihood at new point */
	s2 = update_fcast_errs(e, y, coeff, an, p, q);
	ll = get_ll(e, an, maxlag);

	/* subiterations for best steplength */
	subiters = 0;
	while (steplength > 1.0e-08 && subiters++ < subitermax && ll < ll_prev) {

	    /* if we've gone down, halve steplength and go back */
	    steplength *= 0.5;
	    for (i=0; i<armod.ncoeff; i++) {
		coeff[i] -= d_coef[i] * steplength;
	    }

	    /* compute loglikelihood again */
	    s2 = update_fcast_errs(e, y, coeff, an, p, q);
	    ll = get_ll(e, an, maxlag);
	}

	if (subiters > subitermax) {
	    gretl_errmsg_set(_("Maximum number of subiterations reached: quitting"));
	    armod.errcode = E_UNSPEC;
	    goto arma_bailout;
	}

	pprintf(prn, _("  constant        = %.8g (gradient = %#.6g)\n"), 
		coeff[0], armod.coeff[0]);
	for (i=0; i<p; i++) {
	    pprintf(prn, _("  ar%d coefficient = %.8g (gradient = %#.6g)\n"), 
		    i + 1, coeff[1+i], armod.coeff[1+i]);
	}
	for (i=0; i<q; i++) {
	    pprintf(prn, _("  ma%d coefficient = %.8g (gradient = %#.6g)\n"), 
		    i + 1, coeff[1+p+i], armod.coeff[1+p+i]);
	}
	pprintf(prn, _("  steplength = %.6g (subiterations = %d)\n"), 
		steplength, subiters);

	/* update criterion */
	crit = ll - ll_prev;
	pprintf(prn, _("  criterion = %g\n\n"), crit);

	/* try and re-double the steplength for next iteration */
	if (subiters == 1 && steplength < 4.0){
	    steplength *= 2.0;
	} else if (subiters == 0) {
	    /* time to quit */
	    break;
	}

	ll_prev = ll;

    }

    if (crit > tol) {
	armod.errcode = E_NOCONV;
    } else {
	cmplx *roots;

	y = Z[v];
	armod.lnL = ll; 
	rewrite_arma_model_stats(&armod, coeff, list, y, e, pdinfo);

	/* compute and save polynomial roots */
	roots = arma_roots(p, q, coeff);
	if (roots != NULL) {
	    gretl_model_set_data(&armod, "roots", roots,
				 (p + q) * sizeof *roots);
	}
    }

    if (!armod.errcode) {
	add_arma_varnames(&armod, pdinfo);
    }

 arma_bailout:

    free(alist);
    free(coeff);
    free(d_coef);
    free_Z(aZ, ainfo);
    clear_datainfo(ainfo, CLEAR_FULL);
    free(ainfo);

    return armod;
}



    

    

    
    
