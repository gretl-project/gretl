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
				 const double **X,
				 const double *coeff, 
				 const DATAINFO *ainfo,
				 int p, int q, int r)
{
    int i, t;
    double s2 = 0.0;
    const double *ar_coeff = coeff;
    const double *ma_coeff = coeff + p;
    const double *reg_coeff = coeff + p + q;

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

	for (i=1; i<=r; i++) {
	    e[t] -= reg_coeff[i] * X[i-1][t];
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

static void update_error_partials (double *e, 
				   const double *y, const double **X,
				   double **Z, const double *ma_coeff,
				   const DATAINFO *ainfo,
				   int p, int q, int r)
{
    int t, col, i, j, t2;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {

	col = 2;

	/* the constant term (de_c) */
	Z[col][t] = -1.0;
	for (i=1; i<=q; i++) {
	    Z[col][t] -= ma_coeff[i] * Z[col][t-i];
	}

	/* AR terms (de_a) */
	for (j=0; j<p; j++) {
	    col++;
	    if (t >= col-2) {
		Z[col][t] = -y[t-col+2];
		for (i=1; i<=q; i++) {
		    Z[col][t] -= ma_coeff[i] * Z[col][t-i];
		}
	    }
	}

	/* MA terms (de_m) */
	for (j=0; j<q; j++) {
	    col++;
	    t2 = col - p - 2;
	    if (t >= t2) {
		Z[col][t] = -e[t-t2];
		for (i=1; i<=q; i++) {
		    Z[col][t] -= ma_coeff[i] * Z[col][t-i];
		}
	    }
	}

	/* ordinary regressors */
	for (j=0; j<r; j++) {
	    col++;
	    Z[col][t] = -X[j][t]; 
	    for (i=1; i<=q; i++) {
		Z[col][t] -= ma_coeff[i] * Z[col][t-i];
	    }
	}	    
    }
}

static void update_ll_partials (double *e, double **Z, double s2,
				const DATAINFO *ainfo, 
				int p, int q, int r)
{
    int i, t, col;
    int nc = p + q + r + 1;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	double x = e[t] / s2;

	for (i=0; i<nc; i++) {
	    col = 3 + p + q + r + i;
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
    int p = list[1], q = list[2], r = list[0] - 4;
    double mean_error;

    pmod->ci = ARMA;
    pmod->ifc = 1;

    pmod->dfn = p + q + r;
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

    X = malloc(nv * sizeof *X);
    if (X == NULL) return NULL;

    for (i=0; i<nv; i++) {
	v = list[i + 5];
#ifdef ARMA_DEBUG
	fprintf(stderr, "make_armax_X: setting X[%d] -> Z[%d]\n", i, v);
#endif
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

#ifdef ARMA_DEBUG
    fprintf(stderr, "\n");
    for (i=0; i<=p + q + r; i++) {
	fprintf(stderr, "ar_init_by_ols: coeff[%d] = %g\n", 
		i, coeff[i]);
    }
#endif

    clear_model(&armod, NULL);

    return err;
}

#ifdef ARMA_DEBUG
static void make_tmp_varnames (DATAINFO *ainfo, int p, int q, int r)
{
    int i;
    int nt = p + q + r;

    strcpy(ainfo->varname[1], "e");
    strcpy(ainfo->varname[2], "de_c");
    strcpy(ainfo->varname[2+nt+1], "dl_c");

    for (i=0; i<p; i++) {
	sprintf(ainfo->varname[3+i], "de_a%d", i + 1);
	sprintf(ainfo->varname[3+nt+1+i], "dl_a%d", i + 1);
    }

    for (i=0; i<q; i++) {
	sprintf(ainfo->varname[3+p+i], "de_m%d", i + 1);
	sprintf(ainfo->varname[3+nt+1+p+i], "dl_m%d", i + 1);
    }

    for (i=0; i<r; i++) {
	sprintf(ainfo->varname[3+p+q+i], "de_x%d", i + 1);
	sprintf(ainfo->varname[3+nt+1+p+q+i], "dl_x%d", i + 1);
    }   
}
#endif

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

MODEL arma_model (int *list, const double **Z, DATAINFO *pdinfo, 
		  PRN *prn)
{
    int nc, v, p, q, r, maxlag;
    int subiters, iters, itermax;
    int arma_t1, arma_t2;
    int i, t;
    double s2, tol, crit, ll = 0.0;
    double steplength, ll_prev = -1.0e+8;
    double *e, *coeff, *d_coef;
    const double *y;
    double **aZ = NULL;
    const double **X = NULL;
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
    maxlag = (p > q)? p : q; /* maximum lag in the model */

    /* number of ordinary regressors */
    r = list[0] - 4;

    /* adjust sample? */
    if (adjust_sample(pdinfo, Z, list, &arma_t1, &arma_t2)) {
        armod.errcode = E_DATA;
        return armod;
    }

    /* number of coefficients */
    nc = 1 + p + q + r;

    alist = malloc((nc + 2) * sizeof *alist);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	return armod;
    }

    alist[0] = nc + 1;
    alist[1] = 0; /* dep var is constant, in OPG */
    for (i=0; i<=p+q+r; i++) {
	alist[i+2] = 3 + p + q + r + i;
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

    /* initialize the coefficients: AR and regression part by OLS, 
       MA at 0 */
    ar_init_by_ols(list, coeff, Z, pdinfo, ainfo->t1);

    /* initialize forecast errors and derivatives */
    for (t=0; t<ainfo->n; t++) {
	for (i=1; i<ainfo->v; i++) {
	    aZ[i][t] = 0.0;
	}
    }

    /* construct virtual dataset for real regressors */
    if (r > 0) {
	X = make_armax_X(list, Z);
	if (X == NULL) {
	    armod.errcode = E_ALLOC;
	    free(alist);
	    free(coeff);
	    free(d_coef);
	    return armod;
	}
    }

    /* forecast errors */
    e = aZ[1];

#ifdef ARMA_DEBUG
    make_tmp_varnames(ainfo, p, q, r);
#endif 

    crit = 1.0;
    tol = 1.0e-6;
    iters = 0;
    itermax = get_maxiter();
    steplength = 0.125;

    /* generate one-step forecast errors */
    s2 = update_fcast_errs(e, y, X, coeff, ainfo, p, q, r);

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
	update_error_partials(e, y, X, aZ, coeff + p, ainfo, p, q, r);

	/* partials of l wrt coeffs */
	update_ll_partials(e, aZ, s2, ainfo, p, q, r);

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
	s2 = update_fcast_errs(e, y, X, coeff, ainfo, p, q, r);
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
	    s2 = update_fcast_errs(e, y, X, coeff, ainfo, p, q, r);
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
	int qr_bak = get_use_qr();

	/* run OPG once more using QR, to get VCV matrix */
	clear_model(&armod, NULL);
	set_use_qr(1);
	armod = lsq(alist, &aZ, ainfo, OLS, OPT_A, 0.0);
	set_use_qr(qr_bak);

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
    free(X);

    return armod;
}
