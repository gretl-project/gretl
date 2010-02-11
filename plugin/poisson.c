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

/* models for count data */

#include "libgretl.h"
#include "matrix_extra.h"
#include "libset.h"
#include "bhhh_max.h"

#define PDEBUG 0
#define PR2DEBUG 0

#define POISSON_TOL 1.0e-10 
#define POISSON_MAX_ITER 100 

typedef struct negbin_info_ negbin_info;

struct negbin_info_ {
    int k;               /* length of parameter vector */
    gretl_vector *y;     /* dependent variable */
    gretl_matrix *X;     /* regressors */
    gretl_matrix *theta; /* parameter vector */
    gretl_matrix *mu;    /* exp(X\beta) */
    gretl_matrix *llt;   /* per-observation likelihood */
    gretl_matrix *G;     /* score matrix */
    gretl_matrix *V;     /* covariance matrix */
    PRN *prn;            /* verbose printer */
};

static void negbin_free (negbin_info *nbinfo)
{
    gretl_matrix_replace(&nbinfo->y, NULL);
    gretl_matrix_replace(&nbinfo->X, NULL);
    gretl_matrix_replace(&nbinfo->theta, NULL);
    gretl_matrix_replace(&nbinfo->mu, NULL);
    gretl_matrix_replace(&nbinfo->llt, NULL);
    gretl_matrix_replace(&nbinfo->G, NULL);
    gretl_matrix_replace(&nbinfo->V, NULL);
}

static int negbin_init (negbin_info *nbinfo, MODEL *pmod,
			const double **Z, gretlopt opt, 
			PRN *prn)
{
    int n = pmod->nobs;
    int k = pmod->ncoeff + 1;
    int i, s, t, v;

    nbinfo->y = gretl_column_vector_alloc(n);
    nbinfo->X = gretl_matrix_alloc(n, k);
    nbinfo->theta = gretl_column_vector_alloc(k + 1);
    nbinfo->mu = gretl_column_vector_alloc(n);
    nbinfo->llt = gretl_column_vector_alloc(n);
    nbinfo->G = gretl_matrix_alloc(n, k);
    nbinfo->V= gretl_matrix_alloc(k, k);

    if (nbinfo->y == NULL || nbinfo->X == NULL ||  
	nbinfo->theta == NULL || nbinfo->mu == NULL ||
	nbinfo->llt == NULL || nbinfo->G == NULL || 
	nbinfo->V == NULL) {
	negbin_free(nbinfo);
	return E_ALLOC;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	v = pmod->list[1];
	nbinfo->y->val[s] = Z[v][t];
	for (i=0; i<k; i++) {
	    v = pmod->list[i+2];
	    gretl_matrix_set(nbinfo->X, s, i, Z[v][t]);
	}
	s++;
    }

    nbinfo->k = k;

    for (i=0; i<k; i++) {
	nbinfo->theta->val[i] = 0.01;
    }

    nbinfo->theta->val[k] = 1.0;

    nbinfo->prn = (opt & OPT_V)? prn : NULL;

    return 0;
}

static double negbin_ll (negbin_info *nbinfo, 
			 gretl_matrix *beta, 
			 int k, double alpha)
{
    double eps = 1.0e-9;
    double *ll = nbinfo->llt->val;
    double *mu = nbinfo->mu->val;
    double *y = nbinfo->y->val;
    int t, T = nbinfo->y->rows;
    double psi, mpp;
    double lltot = 0.0;

    gretl_matrix_reuse(beta, k, 1);
    gretl_matrix_multiply(nbinfo->X, beta, nbinfo->mu);
    gretl_matrix_reuse(beta, k+1, 1);

    for (t=0; t<T; t++) {
	mu[t] = eps + exp(mu[t]);
	psi = eps + mu[t] / alpha;
	mpp = mu[t] + psi;
	ll[t] = log_gamma_function(y[t] + psi) 
	    - log_gamma_function(psi);
	ll[t] -= log_gamma_function(y[t] + 1.0);
	ll[t] += psi * log(psi / mpp);
	ll[t] += y[t] * log(mu[t] / mpp);
	lltot += ll[t];
    }

    return lltot;
}

static int negbin_score (void)
{
    /* to be written */
    return 1;
}

static double negbin_callback (double *theta, 
			       gretl_matrix *G, 
			       void *data,
			       int do_score,
			       int *err)
{
    negbin_info *nbinfo = (negbin_info *) data;
    int k = nbinfo->k;
    gretl_matrix *beta = nbinfo->theta;
    double alpha = beta->val[k];
    double ll = NADBL;

    *err = E_NOTIMP;

    if (!*err) {
	ll = negbin_ll(nbinfo, beta, k, alpha);
    }

    if (!*err && do_score) {
	*err = negbin_score();
    }

    return ll;
}

static int do_negbin (MODEL *pmod, const double **Z, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    negbin_info nbinfo;
    double tol = libset_get_double(BHHH_TOLER);
    int iters, err = 0;

    err = negbin_init(&nbinfo, pmod, Z, opt, prn);

    err = E_NOTIMP; /* till ready */

    if (!err) {
	gretlopt bhhh_opt = OPT_NONE;

	if (opt & OPT_V) {
	    bhhh_opt |= OPT_V;
	}

	err = bhhh_max(nbinfo.theta->val, nbinfo.k, nbinfo.G,
		       negbin_callback, tol, &iters,
		       &nbinfo, nbinfo.V, bhhh_opt, nbinfo.prn);
    }

    negbin_free(&nbinfo);

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }
    
    return pmod->errcode;
}

/* sandwich of hessian and OPG */

static int poisson_robust_vcv (MODEL *pmod, gretl_matrix *G)
{
    gretl_matrix *H = NULL;
    gretl_matrix *GG = NULL;
    gretl_matrix *V = NULL;
    int k = G->cols;
    int err = 0;

    H = gretl_vcv_matrix_from_model(pmod, NULL, &err);

    if (!err) {
	GG = gretl_matrix_alloc(k, k);
	V = gretl_matrix_alloc(k, k);
	if (GG == NULL || V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	gretl_matrix_multiply_mod(G, GRETL_MOD_TRANSPOSE,
				  G, GRETL_MOD_NONE,
				  GG, GRETL_MOD_NONE);    
	/* form sandwich: V = H^{-1} GG' H^{-1} */
	err = gretl_matrix_qform(H, GRETL_MOD_NONE,
				 GG, V, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, V);
    }

    if (!err) {
	gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
	pmod->opt |= OPT_R;
    }

    gretl_matrix_free(V);
    gretl_matrix_free(GG);
    gretl_matrix_free(H);

    return err;
} 

/* Overdispersion test via augmented OPG regression: see
   Davidson and McKinnon, ETM, section 11.5 */

static int overdispersion_test (MODEL *pmod, const double **Z,
				gretlopt opt)
{
    const double *mu = pmod->yhat;
    const double *y = Z[pmod->list[1]];
    gretl_matrix *u, *G, *b, *e;
    double mt, mxt, zt;
    int n = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t, v;
    int err = 0;

    u = gretl_unit_matrix_new(n, 1);
    G = gretl_matrix_alloc(n, k+1);
    b = gretl_matrix_alloc(k+1, 1);
    e = gretl_matrix_alloc(n, 1);

    if (u == NULL || G == NULL || b == NULL || e == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* Construct gradient matrix, G, with an extra column
       holding z_t = (y_t - exp(X_t\beta))^2 - y_t.  Under 
       the null of no overdispersion, (n - SSR) from the
       artificial regression follows \chi^2(1) asymptotically.
    */

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(y[t]) || na(mu[t])) {
	    continue;
	}
	mt = y[t] - mu[t];
	for (i=0; i<k; i++) {
	    v = pmod->list[i+2];
	    mxt = mt * Z[v][t];
	    gretl_matrix_set(G, s, i, mxt);
	}
	zt = mt * mt - y[t];
	gretl_matrix_set(G, s, k, zt);
	s++;
    } 

    if (opt & OPT_R) {
	gretl_matrix_reuse(G, n, k);
	err = poisson_robust_vcv(pmod, G);
	gretl_matrix_reuse(G, n, k+1);
    }

    if (!err) {
	err = gretl_matrix_ols(u, G, b, NULL, e, NULL);

	if (!err) {
	    double X2 = e->rows;

	    for (i=0; i<e->rows; i++) {
		X2 -= e->val[i] * e->val[i];
	    }

	    if (X2 > 0) {
		gretl_model_set_double(pmod, "overdisp", X2);
	    }
	}
    }  

 bailout:

    gretl_matrix_free(u);
    gretl_matrix_free(G);
    gretl_matrix_free(b);
    gretl_matrix_free(e);

    return err;
}

static double poisson_ll (const double *y, const double *mu, 
			  int t1, int t2)
{
    double loglik = 0.0;
    double lytfact, llt;
    int t;

    for (t=t1; t<=t2; t++) {
	if (na(y[t]) || na(mu[t])) {
	    continue;
	}
	lytfact = log_x_factorial(y[t]);
	if (na(lytfact)) {
	    loglik = NADBL;
	    break;
	}
	llt = -mu[t] + y[t] * log(mu[t]) - lytfact;
	loglik += llt;
    }  

    return loglik;
}

static void add_pseudoR2 (MODEL *pmod, const double *y, const double *offset, 
			  double offmean)
{
    double llt, ll0 = 0.0;
    double K, lytfact;
    double ybar = gretl_mean(pmod->t1, pmod->t2, y);
    int use_offset = (offset != NULL);
    int t;

    if (use_offset) {
	K = ybar * (log(ybar/offmean) - 1.0);
    } else {
	K = ybar * (log(ybar) - 1.0);
    }

#if PR2DEBUG
    fprintf(stderr, "pseudoR2: K = %g\n", K);
#endif

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(y[t]) || (use_offset && na(offset[t]))) {
	    continue;
	}

	lytfact = log_x_factorial(y[t]);
	if (na(lytfact)) {
	    break;
	}

	llt = K - lytfact;

	if (use_offset) {
	    llt += y[t] * log(offset[t]); 
	}

#if PR2DEBUG
	fprintf(stderr, "ll[%d] = %g\n", t, llt);
#endif
	ll0 += llt;
    }  

#if PR2DEBUG
    fprintf(stderr, "ll0 = %g\n", ll0);
#endif

    if (na(ll0)) {
	pmod->rsq = pmod->adjrsq = NADBL;
    } else {
	int k = pmod->ncoeff; /* FIXME? - pmod->ifc */

	pmod->rsq = 1.0 - (pmod->lnL / ll0);
	pmod->adjrsq = 1.0 - ((pmod->lnL - k) / ll0);
    }
}

static int 
transcribe_poisson_results (MODEL *targ, MODEL *src, const double *y, 
			    int iter, int offvar, const double *offset, 
			    double offmean)
{
    int i, t;
    int err = 0;

    targ->ci = POISSON;
    
    gretl_model_set_int(targ, "iters", iter);

    if (offvar > 0) {
	gretl_model_set_int(targ, "offset_var", offvar);
    }

    targ->ess = 0.0;

    for (t=targ->t1; t<=targ->t2; t++) {
	if (na(targ->yhat[t])) {
	    targ->uhat[t] = NADBL;
	} else {
	    targ->uhat[t] = y[t] - targ->yhat[t];
	    targ->ess += targ->uhat[t] * targ->uhat[t];
	}
    }

    targ->sigma = sqrt(targ->ess / targ->dfd);

    for (i=0; i<targ->ncoeff; i++) {
	targ->sderr[i] = src->sderr[i] / src->sigma;
    }

    targ->lnL = poisson_ll(y, targ->yhat, targ->t1, targ->t2);

    add_pseudoR2(targ, y, offset, offmean);

#if PDEBUG
    fprintf(stderr, "log-likelihood = %g\n", targ->lnL);
#endif

    mle_criteria(targ, 0); 

    /* mask invalid statistics */
    targ->fstt = targ->chisq = NADBL;

    /* make the covariance matrix */
    if (makevcv(src, 1.0)) {
	err = 1;
    } else {
	if (targ->vcv != NULL) {
	    free(targ->vcv);
	}
	targ->vcv = src->vcv;
	src->vcv = NULL;
    }   

    return err;
}

static double *get_offset (MODEL *pmod, int offvar, double **Z,
			   double *offmean)
{
    double *offset = NULL;
    int t, err = 0;

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	} else if (na(Z[offvar][t])) {
	    err = 1;
	} else if (Z[offvar][t] < 0.0) {
	    err = 1;
	} 
    }

    if (err == 0) {
	offset = Z[offvar];
	*offmean = gretl_mean(pmod->t1, pmod->t2, offset);
    }

    return offset;
}

static int 
do_poisson (MODEL *pmod, int offvar, double ***pZ, DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn)
{
    int origv = pdinfo->v;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int iter = 0;
    double crit = 1.0;
    double *offset = NULL;
    double offmean = NADBL;
    double *y;
    double *wgt;
    double *depvar;
    MODEL tmpmod;
    int *local_list = NULL;
    int i, t;

    gretl_model_init(&tmpmod);

    /* set the sample to that of the initial OLS model */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    local_list = gretl_list_new(pmod->list[0] + 1);
    if (local_list == NULL) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    if (dataset_add_series(2, pZ, pdinfo)) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    if (offvar > 0) {
	offset = get_offset(pmod, offvar, *pZ, &offmean);
	if (offset == NULL) {
	    pmod->errcode = E_DATA;
	    goto bailout;
	}
    }

    /* the original dependent variable */
    y = (*pZ)[pmod->list[1]];

    /* weighting variable (first newly added var) */
    local_list[1] = origv;
    wgt = (*pZ)[origv];

    /* dependent variable for GNR (second newly added var) */
    local_list[2] = origv + 1;
    depvar = (*pZ)[origv + 1];
    
    for (i=3; i<=local_list[0]; i++) { 
	/* original independent vars */
	local_list[i] = pmod->list[i-1];
    }    

    pmod->coeff[0] = log(pmod->ybar);
    if (offvar > 0) {
	pmod->coeff[0] -= log(offmean);
    }

    for (i=1; i<pmod->ncoeff; i++) { 
	pmod->coeff[i] = 0.0;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    depvar[t] = NADBL;
	    wgt[t] = NADBL;
	} else {
	    pmod->yhat[t] = pmod->ybar;
	    if (offvar > 0) {
		pmod->yhat[t] *= offset[t] / offmean;
	    }
	    depvar[t] = y[t] / pmod->yhat[t] - 1.0;
	    wgt[t] = pmod->yhat[t];
	}
    }

    pputc(prn, '\n');

    while (iter < POISSON_MAX_ITER && crit > POISSON_TOL) {

	iter++;

	tmpmod = lsq(local_list, pZ, pdinfo, WLS, OPT_A);

	if (tmpmod.errcode) {
	    fprintf(stderr, "poisson_estimate: lsq returned %d\n", 
		    tmpmod.errcode);
	    pmod->errcode = tmpmod.errcode;
	    break;
	}

	crit = tmpmod.nobs * tmpmod.rsq;

	pprintf(prn, "%s %3d\tcrit = %g\n", _("iteration"), iter, crit);

	for (i=0; i<tmpmod.ncoeff; i++) { 
	    pmod->coeff[i] += tmpmod.coeff[i];
#if PDEBUG
	    fprintf(stderr, "coeff[%d] = %g,\tgrad[%d] = %g\n", 
		    i, pmod->coeff[i], i, tmpmod.coeff[i]);
#endif
	}

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (na(pmod->uhat[t])) {
		continue;
	    }
	    pmod->yhat[t] *= exp(tmpmod.yhat[t]);
	    depvar[t] = y[t] / pmod->yhat[t] - 1;
	    wgt[t] = pmod->yhat[t];
	}

	if (crit > POISSON_TOL) {
	    clear_model(&tmpmod);
	}
    }

    pputc(prn, '\n');

    if (crit > POISSON_TOL) {
	pmod->errcode = E_NOCONV;
    } 

    if (pmod->errcode == 0) {
	transcribe_poisson_results(pmod, &tmpmod, y, iter, offvar, 
				   offset, offmean);
	overdispersion_test(pmod, (const double **) *pZ, opt);
    }

 bailout:

    clear_model(&tmpmod);
    free(local_list);
    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return pmod->errcode;
}

int count_data_estimate (MODEL *pmod, int ci, int offvar, 
		     double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn) 
{
    int err = 0;

    if (!gretl_iscount(pmod->t1, pmod->t2, (*pZ)[pmod->list[1]])) {
	gretl_errmsg_set(_("poisson: the dependent variable must be count data"));
	err = pmod->errcode = E_DATA;
    } else if (ci == POISSON) {
	err = do_poisson(pmod, offvar, pZ, pdinfo, opt, prn);
    } else {
	err = do_negbin(pmod, (const double **) *pZ, pdinfo, opt, prn);
    }

    return err;
}

