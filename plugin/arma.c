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

#include "../cephes/polrt.c"

#undef ARMA_DEBUG

#define DEFAULT_MAX_ITER 1000
#define SUBITERMAX         40
#define MINSTEPLEN    1.0e-08     

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
    char *mistr = getenv("GRETL_MAX_ITER");
    int mi = DEFAULT_MAX_ITER;

    if (mistr != NULL) {
	if (!sscanf(mistr, "%d", &mi)) {
	    mi = DEFAULT_MAX_ITER;
	}
    }

    return mi;
}

static double update_fcast_errs (double *e, const double *y, 
				 const double *coeff, 
				 const DATAINFO *ainfo,
				 int p, int q)
{
    int i, t;
    double s2 = 0.0;
    const double *ar_coeff = coeff;
    const double *ma_coeff = coeff + p;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {

	e[t] = y[t] - coeff[0];

#ifdef ARMA_DEBUG
	if (t < ainfo->t1 + 3) {
	    fprintf(stderr, "update_fcast_errs: initializing e[%d] = %g - %g = %g\n", 
		    t, y[t], coeff[0], e[t]);
	}
#endif

	for (i=1; i<=p; i++) {
	    e[t] -= ar_coeff[i] * y[t-i];
	}

	for (i=1; i<=q; i++) {
	    if (t - i >= ainfo->t1) {
		e[t] -= ma_coeff[i] * e[t-i];
	    }
	}

#ifdef ARMA_DEBUG
	if (t < ainfo->t1 + 3) {
	    fprintf(stderr, "update_fcast_errs: after processing: e[%d] = %g\n", 
		    t, e[t]);
	}
#endif

	s2 += e[t] * e[t];
    }

    s2 /= (double) (ainfo->t2 - ainfo->t1 + 1);

#ifdef ARMA_DEBUG
    fprintf(stderr, "update_fcast_errs: SampleSize = %d, s2 = %g\n", 
	    ainfo->t2 - ainfo->t1 + 1, s2);
#endif

    return s2;
}

static void update_error_partials (double *e, const double *y, 
				   double **Z, const double *ma_coeff,
				   const DATAINFO *ainfo,
				   int p, int q)
{
    int t, col, i, t2;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {

	/* the constant term (de_c) */
	col = 2;
	Z[col][t] = -1.0;
	for (i=1; i<=q; i++) {
	    Z[col][t] -= ma_coeff[i] * Z[col][t-i];
	}

	/* AR terms (de_a) */
	if (p > 0) {
	    for (col=3; col <= p+3; col++) {
		if (t >= col-2) {
		    Z[col][t] = -y[t-col+2];
		    for (i=1; i<=q; i++) {
			Z[col][t] -= ma_coeff[i] * Z[col][t-i];
		    }
		}
	    }
	}

	/* MA terms (de_m) */
	if (q > 0) {
	    for (col=p+3; col <= p+q+3; col++) {
		t2 = col - p - 2;
		if (t >= t2) {
		    Z[col][t] = -e[t-t2];
		    for (i=1; i<=q; i++) {
			Z[col][t] -= ma_coeff[i] * Z[col][t-i];
		    }
		}
	    }
	}
    }
}

static void update_ll_partials (double *e, double **Z, double s2,
				const DATAINFO *ainfo, int p, int q)
{
    int i, t, col;
    int nc = p + q + 1;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	double x = e[t] / s2;

	for (i=0; i<nc; i++) {
	    col = 3 + p + q + i;
	    Z[col][t] = -Z[col-nc][t] * x;
	}
    }
}

static double get_ll (const double *e, const DATAINFO *ainfo)
{
    int t;
    double ll = 0.0, s2;
    int SampleSize = ainfo->t2 - ainfo->t1 + 1;
    const double K = 1.41893853320467274178; /* ln(sqrt(2*pi)) + 0.5 */

    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	ll += e[t] * e[t];
    }

    s2 = ll / SampleSize;
    ll = -SampleSize * (0.5 * log(s2) + K);

#ifdef ARMA_DEBUG
    fprintf(stderr, "get_ll: s2 = %g, ll = %g\n", s2, ll);
#endif

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

static void rewrite_arma_model_stats (MODEL *pmod, const double *coeff,
				      const int *list, const double *y, 
				      const double *e, 
				      const DATAINFO *pdinfo)
{
    int i, t;
    int p = list[1], q = list[2];
    double mean_error;

    pmod->ci = ARMA;
    pmod->ifc = 1;

    pmod->dfn = p + q;
    pmod->dfd = pmod->nobs - pmod->dfn;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = coeff[i];
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
    pmod->criterion[1] = -2.0 * pmod->lnL + 2.0 * (pmod->ncoeff + 1);
    /* BIC, as per X-12-ARIMA */
    pmod->criterion[4] = -2.0 * pmod->lnL + (pmod->ncoeff + 1) * log(pmod->nobs);
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
   AR coefficients */

static int ar_init_by_ols (int v, int p, double *coeff,
			   const double **Z, const DATAINFO *pdinfo,
			   int arma_t1)
{
    int an = pdinfo->t2 - arma_t1 + 1;
    double **aZ = NULL;
    DATAINFO *ainfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int i, t, err = 0;

    gretl_model_init(&armod, NULL);  

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
    
    for (t=0; t<an; t++) {
	for (i=0; i<=p; i++) {
	    int s = t + arma_t1 - i;

	    aZ[i+1][t] = Z[v][s];
	}
    }

    armod = lsq(alist, &aZ, ainfo, OLS, OPT_A, 0.0);
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

#ifdef ARMA_DEBUG
    fprintf(stderr, "\n");
    for (i=0; i<armod.ncoeff; i++) {
	fprintf(stderr, "ar_init_by_ols: coeff[%d] = %g\n", 
		i, coeff[i]);
    }
#endif

    return err;
}

#ifdef ARMA_DEBUG
static void make_tmp_varnames (DATAINFO *ainfo, int p, int q)
{
    int i;

    strcpy(ainfo->varname[1], "e");
    strcpy(ainfo->varname[2], "de_c");
    strcpy(ainfo->varname[2+p+q+1], "dl_c");

    for (i=0; i<p; i++) {
	sprintf(ainfo->varname[3+i], "de_a%d", i + 1);
	sprintf(ainfo->varname[3+p+q+1+i], "dl_a%d", i + 1);
    }

    for (i=0; i<q; i++) {
	sprintf(ainfo->varname[3+p+i], "de_m%d", i + 1);
	sprintf(ainfo->varname[3+2*p+q+1+i], "dl_m%d", i + 1);
    }   
}
#endif

static int adjust_sample (DATAINFO *pdinfo, const double *y,
                          int p, int q, int v,
			  int *arma_t1, int *arma_t2)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int maxlag = (p>q)? p : q;
    int an, t, t1min = 0;

    for (t=0; t<=pdinfo->t2; t++) {
        if (na(y[t])) t1min++;
        else break;
    }

    t1min += maxlag;
    if (t1 < t1min) t1 = t1min;

    for (t=pdinfo->t2; t>=t1; t--) {
        if (na(y[t])) t2--;
        else break;
    }

    for (t=t1-p; t<t2; t++) {
        if (na(y[t])) {
            char msg[64];

            sprintf(msg, _("Missing value encountered for "
                           "variable %d, obs %d"), v, t + 1);
	    gretl_errmsg_set(msg);
            return 1;
        }
    }

    an = t2 - t1 + 1;
    if (an <= p + q + 1) return 1; 

    *arma_t1 = t1;
    *arma_t2 = t2;

    return 0;
}

MODEL arma_model (int *list, const double **Z, DATAINFO *pdinfo, 
		  PRN *prn)
{
    int nc, v, p, q, maxlag;
    int subiters, iters, itermax;
    int arma_t1, arma_t2;
    int i, t;
    double s2, tol, crit, ll = 0.0;
    double steplength, ll_prev = -1.0e+8;
    double *e, *coeff, *d_coef;
    const double *y;
    double **aZ;
    DATAINFO *ainfo;
    int *alist;
    MODEL armod;
#ifdef ARMA_DEBUG
    PRN *errprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
#endif

    gretl_model_init(&armod, NULL);  

    if (check_arma_list(list)) {
	armod.errcode = E_UNSPEC;
	return armod;
    }

    v = list[4]; /* dependent var */
    y = Z[v];

    p = list[1]; /* AR order */
    q = list[2]; /* MA order */
    maxlag = (p>q) ? p : q; /* maximum lag in the model */

    /* adjust sample? */
    if (adjust_sample(pdinfo, y, p, q, v, &arma_t1, &arma_t2)) {
        armod.errcode = E_DATA;
        return armod;
    }

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
    ainfo = create_new_dataset(&aZ, nc * 2 + 2, pdinfo->n, 0);
    if (ainfo == NULL) {
	free(alist);
	free(coeff);
	free(d_coef);
	armod.errcode = E_ALLOC;
	return armod;
    }
    ainfo->t1 = arma_t1;
    ainfo->t2 = arma_t2;

    /* initialize the coefficients: AR part by OLS, MA at 0 */
    ar_init_by_ols(v, p, coeff, Z, pdinfo, ainfo->t1);
    for (i=0; i<q; i++) coeff[i+p+1] = 0.0;

    /* initialize forecast errors and derivatives */
    for (t=0; t<ainfo->n; t++) {
	for (i=1; i<ainfo->v; i++) {
	    aZ[i][t] = 0.0;
	}
    }

    /* forecast errors */
    e = aZ[1];

#ifdef ARMA_DEBUG
    make_tmp_varnames(ainfo, p, q);
#endif 

    crit = 1.0;
    tol = 1.0e-6;
    iters = 0;
    itermax = get_maxiter();
    steplength = 0.125;

    /* generate one-step forecast errors */
    s2 = update_fcast_errs(e, y, coeff, ainfo, p, q);

    /* calculate log-likelihood */
    ll = get_ll(e, ainfo);

    while (crit > tol && iters++ < itermax && !isnan(ll)) {

	pprintf(prn, "Iteration %d\n", iters);
        pprintf(prn, "  log likelihood = %g\n", ll);

#ifdef ARMA_DEBUG
	fprintf(stderr, "arma main loop: iters = %d, ll = %g, steplength = %g\n", 
		iters, ll, steplength);
#endif

        /* partials of e wrt coeffs */
	update_error_partials(e, y, aZ, coeff + p, ainfo, p, q);

	/* partials of l wrt coeffs */
	update_ll_partials(e, aZ, s2, ainfo, p, q);

	/* OPG regression */
	clear_model(&armod, NULL);
	armod = lsq(alist, &aZ, ainfo, OLS, OPT_A, 0.0);
	if (armod.errcode) {
	    goto arma_bailout;
	}

#ifdef ARMA_DEBUG
	armod.ID = 0;
	printmodel(&armod, ainfo, errprn);
#endif

	/* compute direction */
	for (i=0; i<armod.ncoeff; i++) {
	    d_coef[i] = armod.coeff[i];
	}

	/* update parameters */
	for (i=0; i<armod.ncoeff; i++) {
	    coeff[i] += d_coef[i] * steplength;
	}

	/* compute log-likelihood at new point */
	s2 = update_fcast_errs(e, y, coeff, ainfo, p, q);
	ll = get_ll(e, ainfo);

	/* subiterations for best steplength */
	subiters = 0;
	while (steplength > MINSTEPLEN && subiters++ < SUBITERMAX && ll < ll_prev) {

#ifdef ARMA_DEBUG
	    fprintf(stderr, "arma sub-loop 1: "
		    "subiters = %d, steplength = %g, ll = %g\n",
		    subiters, steplength, ll);
#endif

	    /* if we've gone down, halve steplength and go back */
	    steplength *= 0.5;
	    for (i=0; i<armod.ncoeff; i++) {
		coeff[i] -= d_coef[i] * steplength;
	    }

	    /* compute loglikelihood again */
	    s2 = update_fcast_errs(e, y, coeff, ainfo, p, q);
	    ll = get_ll(e, ainfo);

#ifdef ARMA_DEBUG
	    fprintf(stderr, "arma sub-loop 2: "
		    "subiters = %d, steplength = %g, ll = %g\n",
		    subiters, steplength, ll);
#endif
	}

	if (subiters > SUBITERMAX) {
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
	if (subiters == 1 && steplength < 4.0) {
	    steplength *= 2.0;
	} else if (subiters == 0) {
	    /* time to quit */
	    break;
	}

	ll_prev = ll;

    }

#ifdef ARMA_DEBUG
    gretl_print_destroy(errprn);
#endif

    if (crit > tol) {
	armod.errcode = E_NOCONV;
    } else {
	cmplx *roots;

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
