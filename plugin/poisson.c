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

/* initialize NegBin via Poisson estimates? */
#define POISSON_INIT 1

typedef struct count_info_ count_info;

enum {
    SCORE_UPDATE_MU = 1
};

struct count_info_ {
    int ci;                /* POISSON or NEGBIN */
    int nbtype;            /* NEGBIN only: 1 or 2 */
    int flags;             /* NEGBIN only: control info */
    int k;                 /* number of covariates (including constant) */
    int n;                 /* number of observations */
    double ll;             /* loglikelihood */
    int offvar;            /* ID of offset variable, or 0 */
    double omean;          /* mean of offset variable */
    gretl_matrix *y;       /* dependent variable */
    gretl_matrix *X;       /* covariates */
    gretl_vector *offset;  /* offset variable, or NULL */
    gretl_matrix *logoff;  /* log of offset variable, or NULL */
    gretl_matrix *b;       /* coeffs on covariates */
    gretl_matrix *Xb;      /* X \times \beta (+ offset) */
    gretl_matrix *mu;      /* exp(X\beta) */
    gretl_matrix *HX;      /* workspace */
    gretl_matrix *V;       /* covariance matrix */
    gretl_matrix *G;       /* score matrix */
    double *theta;         /* for NEGBIN */
    PRN *prn;              /* verbose printer */
};

static int poisson_hessian (double *theta, gretl_matrix *H,
			    void *data);

static void negbin_free (count_info *cinfo)
{
    free(cinfo->theta);
    gretl_matrix_free(cinfo->offset);
}

static int cinfo_allocate (count_info *cinfo,
			   gretl_matrix_block **pB)
{
    gretl_matrix_block *B;
    int n = cinfo->n;
    int k = cinfo->k;

    if (cinfo->ci == POISSON) {
	B = gretl_matrix_block_new(&cinfo->y, n, 1,
				   &cinfo->X, n, k,
				   &cinfo->b, k, 1,
				   &cinfo->Xb, n, 1,
				   &cinfo->mu, n, 1,
				   &cinfo->HX, n, k,
				   &cinfo->V, k, k,
				   NULL);
    } else {
	int np = k + 1;

	B = gretl_matrix_block_new(&cinfo->y, n, 1,
				   &cinfo->X, n, k,
				   &cinfo->b, k, 1,
				   &cinfo->mu, n, 1,
				   &cinfo->G, n, np,
				   NULL);
    }

    if (B == NULL) {
	return E_ALLOC;
    } else {
	*pB = B;
	return 0;
    }
}

/* copy data (y, X and offset, if applicable) in matrix form */

static int cinfo_add_data (count_info *cinfo, MODEL *pmod,
			   const DATASET *dset)
{
    const double *y = dset->Z[pmod->list[1]];
    const double *off = NULL;
    int i, s, t, vi;

    if (cinfo->offvar > 0) {
	off = dset->Z[cinfo->offvar];
    }

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	cinfo->y->val[s] = y[t];
	for (i=2; i<=pmod->list[0]; i++) {
	    vi = pmod->list[i];
	    gretl_matrix_set(cinfo->X, s, i-2, dset->Z[vi][t]);
	}
	if (off != NULL) {
	    if (cinfo->ci == NEGBIN) {
		cinfo->offset->val[s] = off[t];
	    } else {
		cinfo->logoff->val[s] = log(off[t]);
	    }
	}
	s++;
    }

    return 0;
}

/* If the count specification does not include a constant term we
   can't use the nice short-cut that's available otherwise, so here
   we run an auxiliary regression to get a proper starting point for
   the iteratively weighted least squares routine.
*/

static int get_noconst_estimates (count_info *cinfo)
{
    gretl_matrix *lny1;
    int i, t, err = 0;

    /* borrow a vector of the right length */
    lny1 = (cinfo->mu != NULL)? cinfo->mu : cinfo->Xb;

    for (t=0; t<cinfo->n; t++) {
	lny1->val[t] = log(cinfo->y->val[t] + 1);
    }

    err = gretl_matrix_ols(lny1, cinfo->X, cinfo->b,
			   NULL, NULL, NULL);

    if (cinfo->ci == NEGBIN) {
	for (i=0; i<cinfo->k; i++) {
	    cinfo->theta[i] = cinfo->b->val[i];
	}
    }

    return err;
}

static int count_coeffs_init (count_info *cinfo, MODEL *pmod)
{
    double *targ;
    int i, err = 0;

    targ = (cinfo->ci == NEGBIN)? cinfo->theta : cinfo->b->val;

    if (pmod->ifc) {
	targ[0] = log(pmod->ybar);
	if (cinfo->offvar > 0) {
	    targ[0] -= log(cinfo->omean);
	}
	for (i=1; i<pmod->ncoeff; i++) {
	    targ[i] = 0.0;
	}
    } else {
	err = get_noconst_estimates(cinfo);
    }

    return err;
}

static int negbin_init (count_info *cinfo,
			gretl_matrix_block **pB,
			MODEL *pmod, const DATASET *dset,
			int offvar, double omean,
			gretlopt opt, PRN *prn)
{
    int i, k = pmod->ncoeff;
    int err = 0;

    cinfo->ci = NEGBIN;
    cinfo->nbtype = (opt & OPT_M)? 1 : 2;
    cinfo->flags = 0;
    cinfo->k = k;
    cinfo->n = pmod->nobs;

    cinfo->offset = NULL;
    cinfo->offvar = offvar;
    cinfo->omean = omean;

    cinfo->theta = malloc((k+1) * sizeof *cinfo->theta);
    if (cinfo->theta == NULL) {
	return E_ALLOC;
    }

    if (offvar > 0) {
	cinfo->offset = gretl_column_vector_alloc(cinfo->n);
	if (cinfo->offset == NULL) {
	    return E_ALLOC;
	}
    }

    err = cinfo_allocate(cinfo, pB);
    if (err) {
	return err;
    }

    cinfo_add_data(cinfo, pmod, dset);

#if POISSON_INIT /* initialize via Poisson */
    for (i=0; i<k; i++) {
	cinfo->theta[i] = pmod->coeff[i];
    }
#else
    /* or simple initialization */
    count_coeffs_init(cinfo, pmod);
#endif
    cinfo->theta[k] = 1.0; /* smarter alpha initialization? */

    cinfo->ll = NADBL;
    cinfo->prn = (opt & OPT_V)? prn : NULL;

    return 0;
}

static int negbin_update_mu (count_info *cinfo, const double *theta)
{
    double *mu = cinfo->mu->val;
    int i, t, err = 0;

    for (i=0; i<cinfo->k; i++) {
	cinfo->b->val[i] = theta[i];
    }

    gretl_matrix_multiply(cinfo->X, cinfo->b, cinfo->mu);

    for (t=0; t<cinfo->n && !err; t++) {
	mu[t] = exp(mu[t]);
	if (mu[t] == 0) {
	    err = E_NAN;
	} else if (cinfo->offset != NULL) {
	    mu[t] *= cinfo->offset->val[t];
	}
    }

    return err;
}

static double negbin_loglik (const double *theta, void *data)
{
    count_info *cinfo = (count_info *) data;
    double alpha = theta[cinfo->k];
    double *mu = cinfo->mu->val;
    double *y = cinfo->y->val;
    double psi = 0, lgpsi = 0;
    double mpp, rat, llt;
    int t, err = 0;

    if (alpha <= 0) {
	return NADBL;
    }

    err = negbin_update_mu(cinfo, theta);
    if (err) {
	return NADBL;
    }

    cinfo->ll = 0.0;
    errno = 0;

    if (cinfo->nbtype == 2) {
	/* in NegBin2 psi is the same for all obs, so it can
	   be dealt with outside the loop */
	psi = 1/alpha;
	lgpsi = lngamma(psi);
    }

    for (t=0; t<cinfo->n; t++) {
	if (cinfo->nbtype == 1) {
	    psi = mu[t]/alpha;
	    lgpsi = lngamma(psi);
	}
	mpp = mu[t] + psi;
	rat = psi/mpp;
	llt = lngamma(y[t] + psi) - lgpsi - lngamma(y[t] + 1.0);
	llt += psi * log(rat) + y[t] * log(1-rat);
	cinfo->ll += llt;
    }

    if (errno || get_cephes_errno()) {
	cinfo->ll = NADBL;
    }

    return cinfo->ll;
}

static int negbin_score (double *theta, double *g, int np,
			 BFGS_CRIT_FUNC ll, void *data)
{
    count_info *cinfo = (count_info *) data;
    double dpsi_dmu, dmu_dbi, dpsi_da = 0;
    double dl_dpsi, dl_dmu, dl_da;
    const double *y = cinfo->y->val;
    const double *mu = cinfo->mu->val;
    double alpha = theta[cinfo->k];
    double a2 = alpha * alpha;
    double psi = 0, dgpsi = 0;
    double mpp, gti;
    int i, t, err = 0;

    if (cinfo->flags == SCORE_UPDATE_MU) {
	negbin_update_mu(cinfo, theta);
    }

    if (g != NULL) {
	for (i=0; i<np; i++) {
	    g[i] = 0.0;
	}
    }

    if (cinfo->nbtype == 1) {
	dpsi_dmu = 1/alpha;
    } else {
	psi = 1/alpha;
	dpsi_dmu = 0;
	dgpsi = digamma(psi);
	dpsi_da = -1/a2;
    }

    for (t=0; t<cinfo->n; t++) {
	if (cinfo->nbtype == 1) {
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
	    if (i < cinfo->k) {
		dmu_dbi = mu[t] * gretl_matrix_get(cinfo->X, t, i);
		gti = (dl_dpsi * dpsi_dmu + dl_dmu) * dmu_dbi;
	    } else {
		gti = dl_da;
	    }
	    gretl_matrix_set(cinfo->G, t, i, gti);
	    if (g != NULL) {
		g[i] += gti;
	    }
	}
    }

    return err;
}

static int negbin2_hessian (double *theta, gretl_matrix *H,
			    void *data)
{
    count_info *cinfo = (count_info *) data;
    double *mu = cinfo->mu->val;
    double *y = cinfo->y->val;
    double a = theta[cinfo->k];
    double am1, am2;
    double amin2, amin3;
    double hrs, xtr, xts;
    int k = cinfo->k;
    int t, r, s, j;
    int err = 0;

    gretl_matrix_zero(H);

    amin2 = pow(a, -2);
    amin3 = pow(a, -3);

    for (t=0; t<cinfo->n; t++) {
	am1 = 1 + a * mu[t];
	am2 = am1 * am1;
	/* beta/beta terms */
	for (r=0; r<k; r++) {
	    xtr = gretl_matrix_get(cinfo->X, t, r);
	    for (s=0; s<=r; s++) {
		xts = gretl_matrix_get(cinfo->X, t, s);
		hrs = gretl_matrix_get(H, r, s);
		hrs += (mu[t] * (1+a*y[t]) * xtr * xts) / am2;
		gretl_matrix_set(H, r, s, hrs);
		if (r != s) {
		    gretl_matrix_set(H, s, r, hrs);
		}
	    }
	}
	/* beta/alpha terms */
	for (r=0; r<k; r++) {
	    xtr = gretl_matrix_get(cinfo->X, t, r);
	    hrs = gretl_matrix_get(H, r, k);
	    hrs += (mu[t] * (y[t] - mu[t]) * xtr) / am2;
	    gretl_matrix_set(H, r, cinfo->k, hrs);
	    gretl_matrix_set(H, cinfo->k, r, hrs);
	}
	/* alpha/alpha term */
	hrs = gretl_matrix_get(H, k, k);
	for (j=0; j<y[t]; j++) {
	    hrs += pow(j/(1+a*j), 2);
	}
	hrs += 2 * amin3 * log(am1);
	hrs -= (2 * amin2 * mu[t]) / am1;
	hrs -= ((y[t] + 1/a) * mu[t] * mu[t]) / am2;
	gretl_matrix_set(H, k, k, hrs);
    }

    return err;
}

static gretl_matrix *negbin_init_H (count_info *cinfo)
{
    gretl_matrix *H = NULL;
    int err;

    cinfo->flags = SCORE_UPDATE_MU;
    err = negbin_score(cinfo->theta, NULL, cinfo->k + 1,
		       NULL, cinfo);
    cinfo->flags = 0;

    if (!err) {
	H = gretl_matrix_GG_inverse(cinfo->G, &err);
    }

    return H;
}

static gretl_matrix *negbin_hessian_inverse (count_info *cinfo,
					     int *err)
{
    gretl_matrix *H = NULL;
    int np = cinfo->k + 1;

    if (cinfo->nbtype == 2) {
	H = gretl_matrix_alloc(np, np);
	if (H == NULL) {
	    *err = E_ALLOC;
	} else {
	    negbin2_hessian(cinfo->theta, H, cinfo);
	    *err = gretl_invert_symmetric_matrix(H);
	}
    } else {
	cinfo->flags = SCORE_UPDATE_MU;
	H = hessian_inverse_from_score(cinfo->theta, np,
				       negbin_score, NULL,
				       cinfo, err);
	cinfo->flags = 0;
    }

    return H;
}

static int negbin_model_add_vcv (MODEL *pmod, count_info *cinfo,
				 const DATASET *dset, gretlopt opt)
{
    gretl_matrix *H = NULL;
    int err = 0;

    if (opt & OPT_G) {
	err = gretl_model_add_OPG_vcv(pmod, cinfo->G, NULL);
    } else {
	H = negbin_hessian_inverse(cinfo, &err);
	if (!err) {
	    if (opt & OPT_R) {
		err = gretl_model_add_QML_vcv(pmod, pmod->ci,
					      H, cinfo->G,
					      dset, opt, NULL);
	    } else {
		err = gretl_model_add_hessian_vcv(pmod, H);
	    }
	}
    }

    gretl_matrix_free(H);

    return err;
}

static int negbin_fill_hatvars (MODEL *pmod, count_info *cinfo)
{
    double *y = cinfo->y->val;
    double u, nb_par, x, l = 0;
    int s, t, k = cinfo->k;
    int err = 0;

    if (cinfo->nbtype == 1) {
	nb_par = cinfo->theta[k];
	l = log1p(nb_par);
    } else {
	nb_par = 1.0 / cinfo->theta[k];
    }

    pmod->ess = 0.0;

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
	if (na(pmod->yhat[t])) {
	    continue;
	}
	pmod->yhat[t] = cinfo->mu->val[s];
	u = y[s] - pmod->yhat[t];
	pmod->ess += u*u;

	/* NOTE: uhat holds the generalized residuals */
	if (cinfo->nbtype == 1) {
	    x = pmod->yhat[t] / nb_par;
	    pmod->uhat[t] = x * (digamma(y[s] + x) - digamma(x) - l);
	} else {
	    x = nb_par / (pmod->yhat[t] + nb_par);
	    pmod->uhat[t] = u * x;
	}
	s++;
    }

    return err;
}

static int
transcribe_negbin_results (MODEL *pmod, count_info *cinfo,
			   const DATASET *dset, int fncount,
			   int grcount, gretlopt opt)
{
    int nc = cinfo->k + 1;
    int i, err = 0;

    pmod->ci = cinfo->ci;

    if (grcount > 0) {
	gretl_model_set_int(pmod, "fncount", fncount);
	gretl_model_set_int(pmod, "grcount", grcount);
    } else {
	gretl_model_set_int(pmod, "iters", fncount);
    }

    if (cinfo->offvar) {
	gretl_model_set_int(pmod, "offset_var", cinfo->offvar);
    }

#if 0
    /* messes up the overall chi-square test */
    pmod->dfd -= 1;
    pmod->dfn += 1;
#endif

    err = negbin_fill_hatvars(pmod, cinfo);
    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    if (!err) {
	err = gretl_model_allocate_param_names(pmod, nc);
	if (!err) {
	    int vi;

	    for (i=0; i<cinfo->k; i++) {
		vi = pmod->list[i+2];
		gretl_model_set_param_name(pmod, i, dset->varname[vi]);
	    }
	    gretl_model_set_param_name(pmod, nc-1, "alpha");
	}
    }

    err = gretl_model_write_coeffs(pmod, cinfo->theta, nc);

    if (!err) {
	err = negbin_model_add_vcv(pmod, cinfo, dset, opt);
    }

    if (!err) {
	pmod->lnL = cinfo->ll;
	mle_criteria(pmod, 0);
	pmod->chisq = wald_omit_chisq(NULL, pmod);
	pmod->fstt = NADBL;
	pmod->rsq = pmod->adjrsq = NADBL;
	if (opt & OPT_M) {
	    /* NegBin1 */
	    pmod->opt |= OPT_M;
	}
    }

    return err;
}

static int do_negbin (MODEL *pmod, int offvar, double omean,
		      DATASET *dset, gretlopt opt, PRN *prn)
{
    gretlopt maxopt = (opt & OPT_V) | OPT_U;
    gretl_matrix_block *B = NULL;
    count_info cinfo = {0};
    double toler;
    int maxit = 100;
    int fncount = 0;
    int grcount = 0;
    int err = 0;

    err = negbin_init(&cinfo, &B, pmod, dset, offvar, omean, opt, prn);

    if (!err && cinfo.nbtype == 2) {
	void *hfunc = negbin2_hessian;
	double crittol = 1.0e-7;
	double gradtol = 1.0e-7;

	err = newton_raphson_max(cinfo.theta, cinfo.k + 1, maxit,
				 crittol, gradtol, &fncount,
				 C_LOGLIK, negbin_loglik,
				 negbin_score, hfunc,
				 &cinfo, maxopt, cinfo.prn);
    } else if (!err) {
	gretl_matrix *H = negbin_init_H(&cinfo);

	BFGS_defaults(&maxit, &toler, NEGBIN);
	err = BFGS_max(cinfo.theta, cinfo.k + 1, maxit, toler,
		       &fncount, &grcount, negbin_loglik, C_LOGLIK,
		       negbin_score, &cinfo, H, maxopt,
		       cinfo.prn);
	gretl_matrix_free(H);
    }

    if (!err) {
	err = transcribe_negbin_results(pmod, &cinfo, dset,
					fncount, grcount, opt);
    }

    gretl_matrix_block_destroy(B);
    negbin_free(&cinfo);

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
	err = gretl_model_add_QML_vcv(pmod, POISSON, H, G,
				      dset, opt, NULL);
    }

    gretl_matrix_free(H);

    return err;
}

/* Overdispersion test via augmented OPG regression: see
   Davidson and McKinnon, ETM, section 11.5. In the process
   we add a robust covariance matrix if wanted.
*/

static int overdispersion_test (MODEL *pmod, count_info *cinfo,
				DATASET *dset, gretlopt opt)
{
    const double *mu = cinfo->mu->val;
    const double *y = cinfo->y->val;
    gretl_matrix *u, *G, *b, *e;
    double mt, xti, zt;
    int n = cinfo->n;
    int k = cinfo->k;
    int i, t;
    int err = 0;

    u = gretl_unit_matrix_new(n, 1);
    G = gretl_matrix_alloc(n, k+1);
    b = gretl_matrix_alloc(k+1, 1);
    e = cinfo->Xb; /* borrowed */

    if (u == NULL || G == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* Construct gradient matrix, G, with an extra column
       holding z_t = (y_t - exp(X_t\beta))^2 - y_t.  Under
       the null of no overdispersion, (n - SSR) from the
       artificial regression follows \chi^2(1) asymptotically.
    */

    for (t=0; t<cinfo->n; t++) {
	mt = y[t] - mu[t];
	for (i=0; i<k; i++) {
	    xti = gretl_matrix_get(cinfo->X, t, i);
	    gretl_matrix_set(G, t, i, mt * xti);
	}
	zt = mt * mt - y[t];
	gretl_matrix_set(G, t, k, zt);
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

    return err;
}

static void add_pseudoR2 (MODEL *pmod, count_info *cinfo)
{
    const double *y = cinfo->y->val;
    const double *logoff = NULL;
    double llt, ll0 = 0.0;
    double K, lytfact;
    int t;

    if (cinfo->offvar) {
	logoff = cinfo->logoff->val;
	K = pmod->ybar * (log(pmod->ybar/cinfo->omean) - 1.0);
    } else {
	K = pmod->ybar * (log(pmod->ybar) - 1.0);
    }

    for (t=0; t<cinfo->n; t++) {
	lytfact = log_x_factorial(y[t]);
	if (na(lytfact)) {
	    ll0 = NADBL;
	    break;
	}
	llt = K - lytfact;
	if (logoff != NULL && y[t] > 0) {
	    llt += y[t] * logoff[t];
	}
	ll0 += llt;
    }

    if (na(ll0)) {
	pmod->rsq = pmod->adjrsq = NADBL;
    } else {
	int k = cinfo->k; /* FIXME? - pmod->ifc */

	pmod->rsq = 1.0 - (pmod->lnL / ll0);
	pmod->adjrsq = 1.0 - ((pmod->lnL - k) / ll0);
    }
}

static int transcribe_poisson_results (MODEL *pmod,
				       count_info *cinfo,
				       int iter)
{
    double *Xb = cinfo->Xb->val;
    double *y = cinfo->y->val;
    int i, s, t, err = 0;

    pmod->ci = cinfo->ci;
    gretl_model_set_int(pmod, "iters", iter);

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = cinfo->b->val[i];
    }

    if (cinfo->offvar) {
	gretl_model_set_int(pmod, "offset_var", cinfo->offvar);
    }

    pmod->ess = 0.0;

    for (t=pmod->t1, s=0; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	pmod->yhat[t] = exp(Xb[s]);
	pmod->uhat[t] = y[s] - pmod->yhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	s++;
    }

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    poisson_hessian(cinfo->b->val, cinfo->V, cinfo);
    err = gretl_invert_symmetric_matrix(cinfo->V);
    if (!err) {
	err = gretl_model_add_hessian_vcv(pmod, cinfo->V);
    }

    pmod->lnL = cinfo->ll;
    add_pseudoR2(pmod, cinfo);
    mle_criteria(pmod, 0);

    /* mask missing statistics */
    pmod->fstt = pmod->chisq = NADBL;

    return err;
}

/* Check the offset series (if present) and retrieve its mean
   if the series is OK. Here we allow for the possibility
   that the offset variable has fewer valid observations than
   the dependent variable and regular regressors; if that's
   the case we need to recalculate the dep. var. statistics.
*/

static int get_offset_info (MODEL *pmod, DATASET *dset,
			    int offvar, double *pmean)
{
    const double *x = dset->Z[offvar];
    double omean = 0.0;
    int ndrop = 0;
    int t, err = 0;

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	} else if (x[t] < 0.0) {
	    err = E_DATA;
	} else if (na(x[t])) {
	    pmod->uhat[t] = pmod->yhat[t] = NADBL;
	    if (pmod->missmask == NULL) {
		model_add_missmask(pmod, dset->n);
	    }
	    if (pmod->missmask == NULL) {
		err = E_ALLOC;
	    } else {
		pmod->missmask[t] = '1';
	    }
	    pmod->nobs -= 1;
	    pmod->dfd -= 1;
	    ndrop++;
	} else {
	    omean += x[t];
	}
    }

    if (!err && ndrop > 0 && pmod->nobs < pmod->ncoeff) {
	/* we don't have enough observations left */
	err = E_DF;
    }

    if (!err && ndrop > 0) {
	/* recalculate dep. var. statistics */
	const double *y = dset->Z[pmod->list[1]];
	double v;

	pmod->ybar = gretl_restricted_mean(pmod->t1, pmod->t2, y,
					   pmod->uhat, OP_NEQ, NADBL);
	v = gretl_restricted_variance(pmod->t1, pmod->t2, y,
				      pmod->uhat, OP_NEQ, NADBL);
	pmod->sdy = sqrt(v);
    }

    if (!err) {
	omean /= pmod->nobs;
	if (omean == 0.0) {
	    err = E_DATA;
	} else {
	    *pmean = omean;
	}
    }

    return err;
}

static void poisson_update_mu (count_info *cinfo, const double *theta)
{
    int i;

    if (theta != cinfo->b->val) {
	for (i=0; i<cinfo->k; i++) {
	    cinfo->b->val[i] = theta[i];
	}
    }
    if (cinfo->logoff != NULL) {
	gretl_matrix_copy_values(cinfo->Xb, cinfo->logoff);
	gretl_matrix_multiply_mod(cinfo->X, GRETL_MOD_NONE,
				  cinfo->b, GRETL_MOD_NONE,
				  cinfo->Xb, GRETL_MOD_CUMULATE);
    } else {
	gretl_matrix_multiply(cinfo->X, cinfo->b, cinfo->Xb);
    }

    for (i=0; i<cinfo->n; i++) {
	cinfo->mu->val[i] = exp(cinfo->Xb->val[i]);
    }
}

static double poisson_loglik (const double *theta, void *data)
{
    count_info *cinfo = (count_info *) data;
    double *y  = cinfo->y->val;
    double *Xb = cinfo->Xb->val;
    double *mu = cinfo->mu->val;
    double llt;
    int t;

    poisson_update_mu(cinfo, theta);
    cinfo->ll = 0.0;
    errno = 0;

    for (t=0; t<cinfo->n; t++) {
	llt = -mu[t] + y[t] * Xb[t] - lngamma(y[t] + 1);
	cinfo->ll += llt;
    }

    if (errno) {
	cinfo->ll = NADBL;
    }

    return cinfo->ll;
}

static int poisson_score (double *theta, double *g, int np,
			  BFGS_CRIT_FUNC ll, void *data)
{
    count_info *cinfo = (count_info *) data;
    double *y  = cinfo->y->val;
    double *mu = cinfo->mu->val;
    double dev;
    int i, t, err = 0;

    for (i=0; i<np; i++) {
	g[i] = 0.0;
    }

    for (t=0; t<cinfo->n; t++) {
	dev = y[t] - mu[t];
	for (i=0; i<np; i++) {
	    g[i] += dev * gretl_matrix_get(cinfo->X, t, i);
	}
    }

    return err;
}

static int poisson_hessian (double *theta, gretl_matrix *H,
			    void *data)
{
    count_info *cinfo = (count_info *) data;
    double *mu = cinfo->mu->val;
    double hxti;
    int i, t, err = 0;

    for (t=0; t<cinfo->n; t++) {
	for (i=0; i<cinfo->k; i++) {
	    hxti = gretl_matrix_get(cinfo->X, t, i);
	    gretl_matrix_set(cinfo->HX, t, i, mu[t] * hxti);
	}
    }

    gretl_matrix_multiply_mod(cinfo->X, GRETL_MOD_TRANSPOSE,
			      cinfo->HX, GRETL_MOD_NONE,
			      H, GRETL_MOD_NONE);

    return err;
}

static int do_poisson (MODEL *pmod, int offvar, double omean,
		       DATASET *dset, gretlopt opt, PRN *prn)
{
    count_info cinfo = {0};
    gretl_matrix_block *B = NULL;
    int iters = 0;
    int err = 0;

    cinfo.ci = POISSON;
    cinfo.n = pmod->nobs;
    cinfo.k = pmod->ncoeff;
    cinfo.offvar = offvar;
    cinfo.omean = omean;
    cinfo.logoff = NULL;
    err = cinfo_allocate(&cinfo, &B);

    if (!err && offvar > 0) {
	cinfo.logoff = gretl_matrix_alloc(cinfo.n, 1);
	if (cinfo.logoff == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	cinfo_add_data(&cinfo, pmod, dset);
    }

    if (!err) {
	err = count_coeffs_init(&cinfo, pmod);
    }

    if (!err) {
	gretlopt maxopt = (opt & OPT_V) | OPT_U;
	double crittol = 1.0e-7;
	double gradtol = 1.0e-7;
	int maxit = 200;

	if (opt & OPT_A) {
	    /* doing prep for negbin */
	    maxopt = OPT_NONE;
	}
	err = newton_raphson_max(cinfo.b->val, cinfo.k, maxit,
				 crittol, gradtol, &iters,
				 C_LOGLIK, poisson_loglik,
				 poisson_score, poisson_hessian,
				 &cinfo, maxopt, prn);
    }

    if (!err) {
	if (opt & OPT_A) {
	    /* we just want the coefficient estimates */
	    int i;

	    for (i=0; i<pmod->ncoeff; i++) {
		pmod->coeff[i] = cinfo.b->val[i];
	    }
	} else {
	    transcribe_poisson_results(pmod, &cinfo, iters);
	    overdispersion_test(pmod, &cinfo, dset, opt);
	    /* do this after settling the variance matrix */
	    pmod->chisq = wald_omit_chisq(NULL, pmod);
	}
    }

    if (err) {
	pmod->errcode = err;
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_free(cinfo.logoff);

    return err;
}

/* the incoming @pmod has been estimated via OLS */

int count_data_estimate (MODEL *pmod, int ci, int offvar,
			 DATASET *dset, gretlopt opt,
			 PRN *prn)
{
    PRN *vprn = (opt & OPT_V)? prn : NULL;
    double omean = NADBL;
    int err = 0;

    if (offvar > 0) {
	/* handle the offset variable */
	err = get_offset_info(pmod, dset, offvar, &omean);
	if (err) {
	    pmod->errcode = err;
	    return err;
	}
    }

    if (opt & OPT_C) {
	/* cluster implies robust */
	opt |= OPT_R;
    }

    if (ci == NEGBIN) {
#if POISSON_INIT
	/* use auxiliary poisson to initialize the estimates */
	err = do_poisson(pmod, offvar, omean, dset, OPT_A, NULL);
	if (!err) {
	    err = do_negbin(pmod, offvar, omean, dset, opt, vprn);
	}
#else
	err = do_negbin(pmod, offvar, omean, dset, opt, vprn);
#endif
    } else {
	err = do_poisson(pmod, offvar, omean, dset, opt, vprn);
    }

    return err;
}
