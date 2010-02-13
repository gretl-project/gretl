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

/* models for count data, Poisson and Negative Binomial */

#include "libgretl.h"
#include "matrix_extra.h"
#include "libset.h"
#include "bhhh_max.h"
#include "gretl_bfgs.h"
#include "../../cephes/libprob.h"

#include <errno.h>

#define PDEBUG 0
#define PR2DEBUG 0

#define POISSON_TOL 1.0e-10 
#define POISSON_MAX_ITER 100 

#define NEGBIN_USE_BFGS 1

typedef struct negbin_info_ negbin_info;
typedef struct offset_info_ offset_info;

struct negbin_info_ {
    int type;              /* variance type: 1 or 2 */
    double ll;             /* loglikelihood */
    int k;                 /* number of regressors */
    double *theta;         /* params array, length k + 1 */
    gretl_matrix_block *B; /* workspace */
    gretl_vector *y;       /* dependent variable */
    gretl_matrix *X;       /* regressors */
    gretl_matrix *beta;    /* coeffs on regressors */
    gretl_matrix *mu;      /* exp(X\beta) */
    gretl_matrix *llt;     /* per-observation likelihood */
    gretl_matrix *G;       /* score matrix */
    gretl_matrix *V;       /* covariance matrix */
    gretl_vector *offset;  /* offset/exposure vector */
    PRN *prn;              /* verbose printer */
};

struct offset_info_ {
    int vnum;        /* ID number of offset variable */
    const double *x; /* series in dataset */
    double mean;     /* sample mean */
};

static void add_pseudoR2 (MODEL *pmod, const double *y, 
			  offset_info *oinfo);

static void negbin_free (negbin_info *nbinfo)
{
    gretl_matrix_block_destroy(nbinfo->B);
    free(nbinfo->theta);
    gretl_matrix_free(nbinfo->offset);
}

static int negbin_init (negbin_info *nbinfo, MODEL *pmod,
			const double **Z, offset_info *oinfo,
			gretlopt opt, PRN *prn)
{
    int n = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t, v;

    /* NEGBIN2 is the default */
    nbinfo->type = (opt & OPT_M)? 1 : 2;

    nbinfo->B = NULL;
    nbinfo->offset = NULL;

    nbinfo->theta = malloc((k + 1) * sizeof *nbinfo->theta);
    if (nbinfo->theta == NULL) {
	return E_ALLOC;
    }

    if (oinfo != NULL) {
	nbinfo->offset = gretl_column_vector_alloc(n);
	if (nbinfo->offset == NULL) {
	    return E_ALLOC;
	}
    }

    nbinfo->B = gretl_matrix_block_new(&nbinfo->y, n, 1,
				       &nbinfo->X, n, k,
				       &nbinfo->beta, k, 1,
				       &nbinfo->mu, n, 1,
				       &nbinfo->llt, n, 1,
				       &nbinfo->G, n, k+1,
				       &nbinfo->V, k+1, k+1,
				       NULL);
    if (nbinfo->B == NULL) {
	return E_ALLOC;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	v = pmod->list[1];
	nbinfo->y->val[s] = Z[v][t];
	if (oinfo != NULL) {
	    nbinfo->offset->val[s] = oinfo->x[t];
	}
	for (i=0; i<k; i++) {
	    v = pmod->list[i+2];
	    gretl_matrix_set(nbinfo->X, s, i, Z[v][t]);
	}
	s++;
    }

    nbinfo->k = k;

    for (i=0; i<k; i++) {
	/* initialize using Poisson estimates */
	nbinfo->theta[i] = pmod->coeff[i];
    }
    nbinfo->theta[k] = 1.0; /* FIXME smart alpha initialization? */

    nbinfo->ll = NADBL;
    nbinfo->prn = (opt & OPT_V)? prn : NULL;

    return 0;
}

/* calculate negative binomial loglikelihood, and score if
   wanted: supports NEGBIN1 and NEGBIN2
*/

static double negbin_callback (double *theta, 
			       gretl_matrix *G, 
			       void *data,
			       int do_score,
			       int *errp)
{
    double eps = 1.111e-16;
    negbin_info *nbinfo = (negbin_info *) data;
    int k = nbinfo->k;
    gretl_matrix *X = nbinfo->X;
    double alpha = theta[k];
    double *ll = nbinfo->llt->val;
    double *mu = nbinfo->mu->val;
    double *y = nbinfo->y->val;
    double psi, mpp;
    int i, t, T = nbinfo->y->rows;
    int err = 0;

    for (i=0; i<k; i++) {
	nbinfo->beta->val[i] = theta[i];
    }

    gretl_matrix_multiply(nbinfo->X, nbinfo->beta, nbinfo->mu);

    nbinfo->ll = 0.0;

    errno = 0;

    for (t=0; t<T; t++) {
	mu[t] = eps + exp(mu[t]);
	if (nbinfo->offset != NULL) {
	    mu[t] *= nbinfo->offset->val[t];
	}
	if (nbinfo->type == 1) {
	    psi = eps + mu[t]/alpha;
	} else {
	    psi = eps + 1/alpha;
	}
	mpp = mu[t] + psi;
	ll[t] = ln_gamma(y[t] + psi) - ln_gamma(psi);
	ll[t] -= ln_gamma(y[t] + 1.0);
	ll[t] += psi * log(psi/mpp) + y[t] * log(mu[t]/mpp);
	nbinfo->ll += ll[t];
    }

    if (errno || get_cephes_errno()) {
	err = E_NAN;
	nbinfo->ll = NADBL;
    }

    if (!err && do_score) {
	double dpsi_da, dpsi_dmu, dmu_dbi;
	double dl_dpsi, dl_dmu, dl_dbi, dl_da;
	double a2 = alpha * alpha;

	for (t=0; t<T; t++) {
	    if (nbinfo->type == 1) {
		psi = eps + mu[t]/alpha;
		dpsi_da = -eps - mu[t]/a2;
		dpsi_dmu = eps + 1/alpha;
	    } else {
		psi = eps + 1/alpha;
		dpsi_da = -eps - 1/a2;
		dpsi_dmu = 0;
	    }	    

	    mpp = mu[t] + psi;

	    dl_dpsi = digamma(psi + y[t]) - digamma(psi) 
		- log(1 + mu[t]/psi) + (mu[t]/mpp) - (y[t]/mpp);

	    dl_da = dl_dpsi * dpsi_da;
	    if (nbinfo->type == 2) {
		dl_da *= alpha;
	    }
	    
	    dl_dmu = y[t] * psi/mu[t] / mpp - (psi/mpp);

	    for (i=0; i<k; i++) {
		dmu_dbi = mu[t] * gretl_matrix_get(X, t, i);
		dl_dbi = dl_dpsi * dpsi_dmu * dmu_dbi + dl_dmu * dmu_dbi;
		gretl_matrix_set(G, t, i, dl_dbi);
	    }
	    gretl_matrix_set(G, t, k, dl_da);
	}
    }

    *errp = err;

    return nbinfo->ll;
}

/* glue needed when using BFGS or forming Hessian */

static double negbin_loglik (const double *theta, void *data)
{
    int err = 0;

    return negbin_callback((double *) theta, NULL, data, 0, &err);
}

/* glue needed when using BFGS */

static int negbin_score (double *theta, double *g, int k, BFGS_CRIT_FUNC ll, 
			 void *data)
{
    negbin_info *nbinfo = (negbin_info *) data;
    int i, t, T = nbinfo->y->rows;
    double gti;
    int err = 0;

    negbin_callback(theta, nbinfo->G, data, 1, &err);

    for (i=0; i<k; i++) {
	g[i] = 0.0;
	for (t=0; t<T; t++) {
	    gti = gretl_matrix_get(nbinfo->G, t, i);
	    g[i] += gti;
	}
    }

    return err;
}

/* QML sandwich VCV */

static int negbin_robust_vcv (MODEL *pmod, negbin_info *nbinfo)
{
    gretl_matrix *H = NULL;
    gretl_matrix *GG = NULL;
    int nc = nbinfo->k + 1;
    int err = 0;

    H = numerical_hessian(nbinfo->theta, nc, negbin_loglik, 
			  nbinfo, &err);

    if (!err) {
	GG = gretl_matrix_alloc(nc, nc);
	if (GG == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	gretl_matrix_multiply_mod(nbinfo->G, GRETL_MOD_TRANSPOSE,
				  nbinfo->G, GRETL_MOD_NONE,
				  GG, GRETL_MOD_NONE);
	/* form sandwich */
	err = gretl_matrix_qform(H, GRETL_MOD_NONE,
				 GG, nbinfo->V, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, nbinfo->V);
	if (!err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
	    pmod->opt |= OPT_R;
	}
    }
    
    gretl_matrix_free(H);
    gretl_matrix_free(GG);

    return err;
}

static int negbin_hessian_vcv (MODEL *pmod, negbin_info *nbinfo)
{
    gretl_matrix *H;
    int nc = nbinfo->k + 1;
    int err = 0;

    H = numerical_hessian(nbinfo->theta, nc, negbin_loglik, 
			  nbinfo, &err);

    if (!err) {
	err = gretl_model_write_vcv(pmod, H);
	if (!err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_HESSIAN);
	}
    }
    
    gretl_matrix_free(H);

    return err;
}

static int 
transcribe_negbin_results (MODEL *pmod, negbin_info *nbinfo, 
			   const double **Z, const DATAINFO *pdinfo, 
			   offset_info *oinfo, int fncount, int grcount,
			   gretlopt opt)
{
    int nc = nbinfo->k + 1;
    int i, s, t, err = 0;

    pmod->ci = NEGBIN;

    gretl_model_set_int(pmod, "fncount", fncount);
    gretl_model_set_int(pmod, "grcount", grcount);

    if (oinfo != NULL) {
	gretl_model_set_int(pmod, "offset_var", oinfo->vnum);
    }

    pmod->ess = 0.0;

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    pmod->yhat[t] = nbinfo->mu->val[s];
	    pmod->uhat[t] = nbinfo->y->val[s] - pmod->yhat[t];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    s++;
	}
    }

    if (!err) {
	err = gretl_model_allocate_params(pmod, nc);
	if (!err) {
	    int v;

	    for (i=0; i<nbinfo->k; i++) {
		v = pmod->list[i+2];
		strcpy(pmod->params[i], pdinfo->varname[v]);
	    }
	    strcpy(pmod->params[nc-1], "alpha");
	}
    }

    pmod->dfd -= 1;
    pmod->dfn += 1;

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    err = gretl_model_write_coeffs(pmod, nbinfo->theta, nc);
    
    if (!err) {
	if (opt & OPT_R) {
	    err = negbin_robust_vcv(pmod, nbinfo);
	} else {
	    err = negbin_hessian_vcv(pmod, nbinfo);
	} 
    }

    if (!err) {
	pmod->lnL = nbinfo->ll;
	mle_criteria(pmod, 0); 
	/* mask invalid statistics */
	pmod->fstt = pmod->chisq = NADBL;
	pmod->rsq = pmod->adjrsq = NADBL;
	pmod->ess = NADBL;
	pmod->sigma = NADBL;
	if (opt & OPT_M) {
	    /* NEGBIN1 */
	    pmod->opt |= OPT_M;
	}
    }

    return err;
}

static int do_negbin (MODEL *pmod, offset_info *oinfo,
		      const double **Z, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn)
{
    negbin_info nbinfo;
    gretlopt max_opt = OPT_NONE;
    double toler;
    int fncount = 0;
    int grcount = 0;
    int err = 0;

    err = negbin_init(&nbinfo, pmod, Z, oinfo, opt, prn);
    if (err) {
	goto bailout;
    }

    if (opt & OPT_V) {
	max_opt |= OPT_V;
    }

    if (getenv("NEGBIN_USE_BHHH")) {
	toler = libset_get_double(BHHH_TOLER);
	err = bhhh_max(nbinfo.theta, nbinfo.k + 1, nbinfo.G,
		       negbin_callback, toler, &fncount, &grcount,
		       &nbinfo, nbinfo.V, max_opt, nbinfo.prn);
    } else {	
	int maxit;

	BFGS_defaults(&maxit, &toler, NEGBIN);

	err = BFGS_max(nbinfo.theta, nbinfo.k + 1, maxit, toler, 
		       &fncount, &grcount, negbin_loglik, C_LOGLIK,
		       negbin_score, &nbinfo, NULL, max_opt, nbinfo.prn);
    }

    if (!err) {
	err = transcribe_negbin_results(pmod, &nbinfo, Z, pdinfo, 
					oinfo, fncount, grcount, 
					opt);
    }

    bailout:

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

static void add_pseudoR2 (MODEL *pmod, const double *y, 
			  offset_info *oinfo)
{
    double llt, ll0 = 0.0;
    double K, lytfact;
    double ybar = gretl_mean(pmod->t1, pmod->t2, y);
    int t;

    if (oinfo != NULL) {
	K = ybar * (log(ybar/oinfo->mean) - 1.0);
    } else {
	K = ybar * (log(ybar) - 1.0);
    }

#if PR2DEBUG
    fprintf(stderr, "pseudoR2: K = %g\n", K);
#endif

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->yhat[t])) {
	    continue;
	}

	lytfact = log_x_factorial(y[t]);
	if (na(lytfact)) {
	    ll0 = NADBL;
	    break;
	}

	llt = K - lytfact;

	if (oinfo != NULL && y[t] > 0) {
	    llt += y[t] * log(oinfo->x[t]); 
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
			    int iter, offset_info *oinfo)
{
    int i, t, err = 0;

    targ->ci = POISSON;
    
    gretl_model_set_int(targ, "iters", iter);

    if (oinfo != NULL) {
	gretl_model_set_int(targ, "offset_var", oinfo->vnum);
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

    add_pseudoR2(targ, y, oinfo);

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

static int do_poisson (MODEL *pmod, offset_info *oinfo, 
		       double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    int origv = pdinfo->v;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int iter = 0;
    double crit = 1.0;
    double *y, *wgt;
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

    /* the original dependent variable */
    y = (*pZ)[pmod->list[1]];

    /* weighting variable (first newly added var) */
    local_list[1] = origv;
    wgt = (*pZ)[origv];

    /* dependent variable for GNR (second newly added var) */
    local_list[2] = origv + 1;
    depvar = (*pZ)[origv + 1];
    
    for (i=3; i<=local_list[0]; i++) { 
	/* the original independent vars */
	local_list[i] = pmod->list[i-1];
    }    

    pmod->coeff[0] = log(pmod->ybar);
    if (oinfo != NULL) {
	pmod->coeff[0] -= log(oinfo->mean);
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
	    if (oinfo != NULL) {
		pmod->yhat[t] *= oinfo->x[t] / oinfo->mean;
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

    if (pmod->errcode == 0 && !(opt & OPT_A)) {
	transcribe_poisson_results(pmod, &tmpmod, y, iter, oinfo);
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

/* Check the offset series (if wanted) and retrieve its mean
   if the series is OK.   
*/

static int get_offset_info (MODEL *pmod, offset_info *oinfo)
{
    int t, err = 0;

    oinfo->mean = 0.0;

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	} else if (na(oinfo->x[t])) {
	    err = E_MISSDATA;
	} else if (oinfo->x[t] < 0.0) {
	    err = E_DATA;
	} else {
	    oinfo->mean += oinfo->x[t];
	}
    }

    if (!err) {
	oinfo->mean /= pmod->nobs;
	if (oinfo->mean == 0.0) {
	    err = E_DATA;
	}
    }

    return err;
}

/* the incoming @pmod has been estimated via OLS */

int count_data_estimate (MODEL *pmod, int ci, int offvar, 
			 double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn) 
{
    offset_info oinfo_t, *oinfo = NULL;
    int err = 0;

    if (offvar > 0) {
	/* handle the offset variable */
	oinfo_t.vnum = offvar;
	oinfo_t.x = (*pZ)[offvar];
	err = get_offset_info(pmod, &oinfo_t);
	if (err) {
	    pmod->errcode = err;
	    return err;
	} else {
	    oinfo = &oinfo_t;
	}
    } 

    if (ci == NEGBIN) {
	/* use auxiliary poisson to initialize the estimates */
	err = do_poisson(pmod, oinfo, pZ, pdinfo, OPT_A, NULL);
	if (!err) {
	    err = do_negbin(pmod, oinfo, (const double **) *pZ, 
			    pdinfo, opt, prn);
	}
    } else {
	err = do_poisson(pmod, oinfo, pZ, pdinfo, opt, prn);
    }

    return err;
}

