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
#include "version.h"
#include "matrix_extra.h"
#include "libset.h"
#include "gretl_bfgs.h"
#include "../../cephes/libprob.h"

#include <errno.h>

#define PDEBUG 0
#define PR2DEBUG 0

#define POISSON_TOL 1.0e-10 
#define POISSON_MAX_ITER 100 

typedef struct negbin_info_ negbin_info;
typedef struct offset_info_ offset_info;

enum {
    SCORE_UPDATE_MU = 1
};

struct negbin_info_ {
    int type;              /* variance type: 1 or 2 */
    int flags;             /* control info */
    double ll;             /* loglikelihood */
    int k;                 /* number of regressors */
    int T;                 /* number of observations */
    double *theta;         /* params array, length k + 1 */
    gretl_matrix_block *B; /* workspace */
    gretl_vector *y;       /* dependent variable */
    gretl_matrix *X;       /* regressors */
    gretl_matrix *beta;    /* coeffs on regressors */
    gretl_matrix *mu;      /* exp(X\beta) */
    gretl_matrix *llt;     /* per-observation likelihood */
    gretl_matrix *G;       /* score matrix */
    gretl_vector *offset;  /* offset/exposure vector */
    PRN *prn;              /* verbose printer */
};

struct offset_info_ {
    int vnum;        /* ID number of offset variable */
    const double *x; /* series in dataset */
    double mean;     /* sample mean */
};

static void negbin_free (negbin_info *nbinfo)
{
    gretl_matrix_block_destroy(nbinfo->B);
    free(nbinfo->theta);
    gretl_matrix_free(nbinfo->offset);
}

static int negbin_init (negbin_info *nbinfo, MODEL *pmod,
			const DATASET *dset, offset_info *oinfo,
			gretlopt opt, PRN *prn)
{
    int n = pmod->nobs;
    int k = pmod->ncoeff;
    int np = k + 1;
    int i, s, t, v;

    /* NegBin2 is the default */
    nbinfo->type = (opt & OPT_M)? 1 : 2;
    nbinfo->flags = 0;

    nbinfo->B = NULL;
    nbinfo->offset = NULL;

    nbinfo->theta = malloc(np * sizeof *nbinfo->theta);
    
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
				       &nbinfo->G, n, np,
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
	nbinfo->y->val[s] = dset->Z[v][t];
	if (oinfo != NULL) {
	    nbinfo->offset->val[s] = oinfo->x[t];
	}
	for (i=0; i<k; i++) {
	    v = pmod->list[i+2];
	    gretl_matrix_set(nbinfo->X, s, i, dset->Z[v][t]);
	}
	s++;
    }

    nbinfo->k = k;
    nbinfo->T = n;

    for (i=0; i<k; i++) {
	/* initialize using Poisson estimates */
	nbinfo->theta[i] = pmod->coeff[i];
    }
    nbinfo->theta[k] = 1.0; /* FIXME smart alpha initialization? */

    nbinfo->ll = NADBL;
    nbinfo->prn = (opt & OPT_V)? prn : NULL;

    return 0;
}

static int negbin_update_mu (negbin_info *nbinfo, const double *theta)
{
    double *mu = nbinfo->mu->val;
    int i, t, err = 0;

    for (i=0; i<nbinfo->k; i++) {
	nbinfo->beta->val[i] = theta[i];
    }

    gretl_matrix_multiply(nbinfo->X, nbinfo->beta, nbinfo->mu);

    for (t=0; t<nbinfo->T && !err; t++) {
	mu[t] = exp(mu[t]);
	if (mu[t] == 0) {
	    err = E_NAN;
	} else if (nbinfo->offset != NULL) {
	    mu[t] *= nbinfo->offset->val[t];
	}
    }

    return err;
}

static double negbin_loglik (const double *theta, void *data)
{
    negbin_info *nbinfo = (negbin_info *) data;
    double alpha = theta[nbinfo->k];
    double *ll = nbinfo->llt->val;
    double *mu = nbinfo->mu->val;
    double *y = nbinfo->y->val;
    double psi = 0, lgpsi = 0;
    double mpp, rat;
    int t, err = 0;

    if (alpha <= 0) {
	return NADBL;
    }

    err = negbin_update_mu(nbinfo, theta);
    if (err) {
	return NADBL;
    }

    nbinfo->ll = 0.0;
    errno = 0;

    if (nbinfo->type == 2) {
	/* in NegBin2 psi is the same for all obs, so it can 
	   be dealt with outside the loop */
	psi = 1/alpha;
	lgpsi = ln_gamma(psi);
    }

    for (t=0; t<nbinfo->T; t++) {
	if (nbinfo->type == 1) {
	    psi = mu[t]/alpha;
	    lgpsi = ln_gamma(psi);
	}

	mpp = mu[t] + psi;
	rat = psi/mpp;

	ll[t] = ln_gamma(y[t] + psi) - lgpsi - ln_gamma(y[t] + 1.0);
	ll[t] += psi * log(rat) + y[t] * log(1-rat);
	nbinfo->ll += ll[t];
    }

    if (errno || get_cephes_errno()) {
	nbinfo->ll = NADBL;
    }

    return nbinfo->ll;
}

static int negbin_score (double *theta, double *g, int np, BFGS_CRIT_FUNC ll, 
			 void *data)
{
    negbin_info *nbinfo = (negbin_info *) data;
    double dpsi_dmu, dmu_dbi, dpsi_da = 0;
    double dl_dpsi, dl_dmu, dl_da;
    const double *y = nbinfo->y->val;
    const double *mu = nbinfo->mu->val;
    double alpha = theta[nbinfo->k];
    double a2 = alpha * alpha;
    double psi = 0, dgpsi = 0;
    double mpp, gti;
    int i, t, err = 0;

    if (nbinfo->flags == SCORE_UPDATE_MU) {
	negbin_update_mu(nbinfo, theta);
    }

    if (g != NULL) {
	for (i=0; i<np; i++) {
	    g[i] = 0.0;
	}
    }    

    if (nbinfo->type == 1) {
	dpsi_dmu = 1/alpha;
    } else {
	psi = 1/alpha;
	dpsi_dmu = 0;
	dgpsi = digamma(psi);
	dpsi_da = -1/a2;
    }	 

    for (t=0; t<nbinfo->T; t++) {
	if (nbinfo->type == 1) {
	    psi = mu[t]/alpha;
	    dgpsi = digamma(psi);
	    dpsi_da = -mu[t]/a2;
	}	    

	mpp = mu[t] + psi;
	dl_dpsi = digamma(psi + y[t]) - dgpsi
	    - log(1 + mu[t]/psi) - (y[t] - mu[t])/mpp;
	dl_da = dl_dpsi * dpsi_da;
	dl_dmu = y[t]/mu[t] - (psi + y[t])/mpp;

	for (i=0; i<np; i++) {
	    if (i < nbinfo->k) {
		dmu_dbi = mu[t] * gretl_matrix_get(nbinfo->X, t, i);
		gti = (dl_dpsi * dpsi_dmu + dl_dmu) * dmu_dbi;
	    } else {
		gti = dl_da;
	    }
	    gretl_matrix_set(nbinfo->G, t, i, gti);
	    if (g != NULL) {
		g[i] += gti;
	    }
	}
    }

    return err;
}

static gretl_matrix *negbin_init_H (negbin_info *nbinfo)
{
    gretl_matrix *H = NULL;
    int err;

    nbinfo->flags = SCORE_UPDATE_MU;
    err = negbin_score(nbinfo->theta, NULL, nbinfo->k + 1, 
		       NULL, nbinfo);
    nbinfo->flags = 0;

    if (!err) {
	H = gretl_matrix_GG_inverse(nbinfo->G, &err);
    }

    return H;
}

static gretl_matrix *negbin_hessian_inverse (negbin_info *nbinfo, 
					     int *err)
{
    gretl_matrix *H;
    int np = nbinfo->k + 1;

    nbinfo->flags = SCORE_UPDATE_MU;
    H = hessian_inverse_from_score(nbinfo->theta, np,
				   negbin_score, NULL,
				   nbinfo, err);
    nbinfo->flags = 0;

    return H;
}

static int negbin_model_add_vcv (MODEL *pmod, negbin_info *nbinfo,
				 const DATASET *dset, gretlopt opt)
{
    gretl_matrix *H = NULL;
    int err = 0;

    if (opt & OPT_G) {
	err = gretl_model_add_OPG_vcv(pmod, nbinfo->G);
    } else {
	H = negbin_hessian_inverse(nbinfo, &err);
	if (!err) {
	    if (opt & OPT_R) {
		err = gretl_model_add_QML_vcv(pmod, pmod->ci, 
					      H, nbinfo->G,
					      dset, opt);
	    } else {
		err = gretl_model_add_hessian_vcv(pmod, H);
	    }
	}
    }

    gretl_matrix_free(H);

    return err;
}

static int 
transcribe_negbin_results (MODEL *pmod, negbin_info *nbinfo, 
			   const DATASET *dset, offset_info *oinfo, 
			   int fncount, int grcount,
			   gretlopt opt)
{
    int nc = nbinfo->k + 1;
    int i, s, t, err = 0;

    pmod->ci = NEGBIN;

    if (grcount > 0) {
	gretl_model_set_int(pmod, "fncount", fncount);
	gretl_model_set_int(pmod, "grcount", grcount);
    } else {
	gretl_model_set_int(pmod, "iters", fncount);
    }

    if (oinfo != NULL) {
	gretl_model_set_int(pmod, "offset_var", oinfo->vnum);
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    pmod->yhat[t] = nbinfo->mu->val[s];
	    pmod->uhat[t] = nbinfo->y->val[s] - pmod->yhat[t];
	    s++;
	}
    }

    if (!err) {
	err = gretl_model_allocate_param_names(pmod, nc);
	if (!err) {
	    int v;

	    for (i=0; i<nbinfo->k; i++) {
		v = pmod->list[i+2];
		gretl_model_set_param_name(pmod, i, dset->varname[v]);
	    }
	    gretl_model_set_param_name(pmod, nc-1, "alpha");
	}
    }

    pmod->dfd -= 1;
    pmod->dfn += 1;

    err = gretl_model_write_coeffs(pmod, nbinfo->theta, nc);
    
    if (!err) {
	err = negbin_model_add_vcv(pmod, nbinfo, dset, opt);
    }

    if (!err) {
	pmod->lnL = nbinfo->ll;
	mle_criteria(pmod, 0); 
	/* mask invalid statistics (FIXME chisq?) */
	pmod->fstt = pmod->chisq = NADBL;
	pmod->rsq = pmod->adjrsq = NADBL;
	pmod->ess = NADBL;
	pmod->sigma = NADBL;
	if (opt & OPT_M) {
	    /* NegBin1 */
	    pmod->opt |= OPT_M;
	}
    }

    return err;
}

static int do_negbin (MODEL *pmod, offset_info *oinfo,
		      DATASET *dset, gretlopt opt, 
		      PRN *prn)
{
    gretlopt max_opt = (opt & OPT_V)? OPT_V : OPT_NONE;
    gretl_matrix *H = NULL;
    negbin_info nbinfo;
    double toler;
    int maxit = 100;
    int fncount = 0;
    int grcount = 0;
    int use_newton = 0;
    int err = 0;

    if (libset_get_int(GRETL_OPTIM) == OPTIM_NEWTON) {
	use_newton = 1;
    }

    err = negbin_init(&nbinfo, pmod, dset, oinfo, opt, prn);

    if (!err && !use_newton) {
	/* initialize BFGS curvature */
	H = negbin_init_H(&nbinfo);
    }

    if (!err && use_newton) {
	double crittol = 1.0e-7;
	double gradtol = 1.0e-7;

	nbinfo.flags = SCORE_UPDATE_MU;
	err = newton_raphson_max(nbinfo.theta, nbinfo.k + 1, maxit, 
				 crittol, gradtol, &fncount, 
				 C_LOGLIK, negbin_loglik, 
				 negbin_score, NULL, 
				 &nbinfo, max_opt, nbinfo.prn);
	nbinfo.flags = 0;
    } else if (!err) {
	BFGS_defaults(&maxit, &toler, NEGBIN);
	err = BFGS_max(nbinfo.theta, nbinfo.k + 1, maxit, toler, 
		       &fncount, &grcount, negbin_loglik, C_LOGLIK,
		       negbin_score, &nbinfo, H, max_opt, 
		       nbinfo.prn);
	gretl_matrix_free(H);
    }

    if (!err) {
	err = transcribe_negbin_results(pmod, &nbinfo, dset, 
					oinfo, fncount, grcount, 
					opt);
    }

    negbin_free(&nbinfo);

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }
    
    return pmod->errcode;
}

/* sandwich of hessian and OPG, or clustered */

static int poisson_robust_vcv (MODEL *pmod, 
			       const gretl_matrix *G,
			       const DATASET *dset, 
			       gretlopt opt)
{
    gretl_matrix *H;
    int err = 0;

    H = gretl_vcv_matrix_from_model(pmod, NULL, &err);

    if (!err) {
	err = gretl_model_add_QML_vcv(pmod, POISSON, 
				      H, G, dset, opt);
    }

    gretl_matrix_free(H);

    return err;
} 

/* Overdispersion test via augmented OPG regression: see
   Davidson and McKinnon, ETM, section 11.5 */

static int overdispersion_test (MODEL *pmod, const DATASET *dset,
				gretlopt opt)
{
    const double *mu = pmod->yhat;
    const double *y = dset->Z[pmod->list[1]];
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
	    mxt = mt * dset->Z[v][t];
	    gretl_matrix_set(G, s, i, mxt);
	}
	zt = mt * mt - y[t];
	gretl_matrix_set(G, s, k, zt);
	s++;
    } 

    if (opt & OPT_R) {
	gretl_matrix_reuse(G, n, k);
	err = poisson_robust_vcv(pmod, G, dset, opt);
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

/* If the Poisson specification does not include a constant term we
   can't use the nice short-cut that's available otherwise, so here
   we run an auxiliary regression to get a proper starting point for
   the iteratively weighted least squares routine.
*/

static int get_noconst_estimates (MODEL *pmod, DATASET *dset, int tmpno)
{
    int yno = pmod->list[1];
    double *tmp = dset->Z[tmpno];
    double *y = dset->Z[yno];
    MODEL logmod;
    int *tmplist;
    int i, t, err = 0;

    tmplist = gretl_list_copy(pmod->list);
    if (tmplist == NULL) {
	return E_ALLOC;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    tmp[t] = NADBL;
	} else {
	    tmp[t] = log(y[t] + 1);
	}
    }

    tmplist[1] = tmpno;
    logmod = lsq(tmplist, dset, OLS, OPT_A);
    err = logmod.errcode;

    if (!err) {
	for (i=0; i<pmod->ncoeff; i++) {
	    pmod->coeff[i] = logmod.coeff[i];
	}
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t])) {
		pmod->yhat[t] = exp(logmod.yhat[t]);
	    }
	}
    }

    free(tmplist);
    clear_model(&logmod);    

    return err;
}

static int do_poisson (MODEL *pmod, offset_info *oinfo, 
		       DATASET *dset, gretlopt opt, 
		       PRN *prn)
{
    int origv = dset->v;
    int orig_t1 = dset->t1;
    int orig_t2 = dset->t2;
    int iter = 0;
    double crit = 1.0;
    double *y, *wgt;
    double *depvar;
    MODEL tmpmod;
    int *local_list = NULL;
    int i, t;

    gretl_model_init(&tmpmod, NULL);

    /* set the sample to that of the initial OLS model */
    dset->t1 = pmod->t1;
    dset->t2 = pmod->t2;

    /* list for iterated WLS */
    local_list = gretl_list_new(pmod->list[0] + 1);
    if (local_list == NULL) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    if (dataset_add_series(dset, 2)) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    /* the original dependent variable */
    y = dset->Z[pmod->list[1]];

    /* weighting variable (the first newly added var) */
    local_list[1] = origv;
    wgt = dset->Z[origv];

    /* dependent variable for GNR (the second newly added var) */
    local_list[2] = origv + 1;
    depvar = dset->Z[origv + 1];

#if PDEBUG
    strcpy(dset->varname[origv], "GNRweight");
    strcpy(dset->varname[origv+1], "GNRy");
#endif

    for (i=3; i<=local_list[0]; i++) { 
	/* insert the original independent vars */
	local_list[i] = pmod->list[i-1];
    }

    if (pmod->ifc) {
	pmod->coeff[0] = log(pmod->ybar);
	if (oinfo != NULL) {
	    pmod->coeff[0] -= log(oinfo->mean);
	}
	for (i=1; i<pmod->ncoeff; i++) { 
	    pmod->coeff[i] = 0.0;
	}
    } else {
	pmod->errcode = get_noconst_estimates(pmod, dset, origv);
	if (pmod->errcode) {
	    goto bailout;
	}
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    depvar[t] = NADBL;
	    wgt[t] = NADBL;
	} else {
	    if (pmod->ifc) {
		pmod->yhat[t] = pmod->ybar;
		if (oinfo != NULL) {
		    pmod->yhat[t] *= oinfo->x[t] / oinfo->mean;
		}
	    } else if (oinfo != NULL) {
		pmod->yhat[t] *= oinfo->x[t];
	    }
	    depvar[t] = y[t] / pmod->yhat[t] - 1.0;
	    wgt[t] = pmod->yhat[t];
	}
    }

    if (prn != NULL) 
	pputc(prn, '\n');

    while (iter < POISSON_MAX_ITER && crit > POISSON_TOL) {

	iter++;

	/* weighted least squares */
	tmpmod = lsq(local_list, dset, WLS, OPT_A);

	if (tmpmod.errcode) {
	    fprintf(stderr, "poisson_estimate: lsq returned %d\n", 
		    tmpmod.errcode);
	    pmod->errcode = tmpmod.errcode;
	    break;
	}

#if PDEBUG
	printmodel(&tmpmod, dset, OPT_S, prn);
#endif

	if (pmod->ifc) {
	    crit = tmpmod.nobs * tmpmod.rsq;
	} else {
	    /* use sum of squared gradients? */
	    crit = 0.0;
	    for (i=0; i<tmpmod.ncoeff; i++) { 
		crit += tmpmod.coeff[i] * tmpmod.coeff[i];
	    }
	}

	if (prn != NULL) {
	    double ll = poisson_ll(y, pmod->yhat, pmod->t1, pmod->t2);

	    print_iter_info(iter, ll, C_LOGLIK, 0, 
			    NULL, NULL, NADBL, prn);
	}

	for (i=0; i<tmpmod.ncoeff; i++) { 
	    pmod->coeff[i] += tmpmod.coeff[i];
	}

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t])) {
		pmod->yhat[t] *= exp(tmpmod.yhat[t]);
		depvar[t] = y[t] / pmod->yhat[t] - 1;
		wgt[t] = pmod->yhat[t];
	    }
	}

	if (crit > POISSON_TOL) {
	    clear_model(&tmpmod);
	}
    }

    if (crit > POISSON_TOL) {
	pmod->errcode = E_NOCONV;
    } 

    if (pmod->errcode == 0 && !(opt & OPT_A)) {
	transcribe_poisson_results(pmod, &tmpmod, y, iter, oinfo);
	overdispersion_test(pmod, dset, opt);
    }

 bailout:

    clear_model(&tmpmod);
    free(local_list);
    dataset_drop_last_variables(dset, dset->v - origv);

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

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
			 DATASET *dset, gretlopt opt, 
			 PRN *prn) 
{
    offset_info oinfo_t, *oinfo = NULL;
    PRN *vprn = NULL;
    int err = 0;

    if (offvar > 0) {
	/* handle the offset variable */
	oinfo_t.vnum = offvar;
	oinfo_t.x = dset->Z[offvar];
	err = get_offset_info(pmod, &oinfo_t);
	if (err) {
	    pmod->errcode = err;
	    return err;
	} else {
	    oinfo = &oinfo_t;
	}
    } 

    if (opt & OPT_C) {
	/* cluster implies robust */
	opt |= OPT_R;
    }

    if (opt & OPT_V) {
	vprn = prn;
    }

    if (ci == NEGBIN) {
	/* use auxiliary poisson to initialize the estimates */
	err = do_poisson(pmod, oinfo, dset, OPT_A, NULL);
	if (!err) {
	    err = do_negbin(pmod, oinfo, dset, opt, vprn);
	}
    } else {
	err = do_poisson(pmod, oinfo, dset, opt, vprn);
    }

    return err;
}

