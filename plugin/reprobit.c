/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2010 Allin Cottrell and Riccardo "Jack" Lucchetti
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

#include "libgretl.h"
#include "version.h"
#include "gretl_model.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"
#include "gretl_normal.h"
#include "libset.h"
#include "gretl_mt.h"

#ifdef _OPENMP
# include <omp.h>
#endif

typedef struct reprob_container_ reprob_container;

struct reprob_container_ {
    const int *list;         /* model specification */
    int depvar;		     /* location of y in array Z */
    int npar;		     /* no. of parameters */
    double ll;		     /* log-likelihood */
    double scale;            /* ln(sigma^2_u) */
    MODEL *pmod;             /* pointer to model struct */

    int N;	             /* number of included units */
    int T;                   /* time dimension of panel */
    int Tmax;                /* max usable time-series observations */
    int *unit_obs;           /* array of effective unit T's */
    int *unit_start;         /* array of starting indices for units */
    int nobs;                /* total number of observations */
    int qp;                  /* number of quadrature points */
    int parallel;

    int *y;	             /* dependent var (0/1) */
    gretl_matrix *X;	     /* main eq. regressors */
    gretl_matrix *R;         /* holds inverse Mills ratios */

    gretl_matrix_block *B;   /* holder for the following */
    gretl_matrix *ndx;       /* index function */
    gretl_matrix *nodes;     /* Gauss-Hermite quadrature nodes */
    gretl_matrix *wts;       /* Gauss-Hermite quadrature weights */
    gretl_matrix *P;         /* probabilities (by individual and qpoints) */
    gretl_matrix *lik;       /* probabilities (by individual) */
    gretl_vector *beta;      /* parameters (excluding log of variance 
				of individual effect) */
    gretl_matrix *qi;        /* storage for score caculation */
};

reprob_container *rep_container_new (const int *list)
{
    reprob_container *C = malloc(sizeof *C);

    if (C != NULL) {
	C->list = list;
	C->depvar = list[1];
	C->npar = list[0];
	C->ll = NADBL;
	C->N = 0;
	C->nobs = 0;
	C->parallel = 0;
	C->qi = NULL;

	C->unit_obs = NULL;
	C->unit_start = NULL;
	
	C->y = NULL;
	C->X = NULL;
	C->R = NULL;
	C->B = NULL;
    }

    return C;
}

static void rep_container_destroy (reprob_container *C)
{
    if (C != NULL) {
	free(C->y);
	free(C->unit_obs);
	gretl_matrix_free(C->X);
	gretl_matrix_free(C->R);
	gretl_matrix_block_destroy(C->B);
	gretl_matrix_free(C->qi);
	free(C);
    }
}

/* based on pooled probit, figure how many cross-sectional
   units are included, and construct an array holding the
   number of observations used for each included unit
*/

static int reprobit_obs_accounts (reprob_container *C,
				  MODEL *pmod,
				  const DATASET *dset)
{
    int unit, ubak = -1;
    int i, t, s;

    C->nobs = pmod->nobs;
    C->T = dset->pd;

    /* count the included units */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    unit = t / C->T;
	    if (unit != ubak) {
		C->N += 1;
	    }
	    ubak = unit;
	}
    }

    C->unit_obs = malloc(2 * C->N * sizeof *C->unit_obs);
    if (C->unit_obs == NULL) {
	return E_ALLOC;
    }

    C->unit_start = C->unit_obs + C->N;

    for (i=0; i<C->N; i++) {
	C->unit_obs[i] = C->unit_start[i] = 0;
    }

    C->Tmax = 0;
    ubak = -1;
    i = s = 0;

    /* build the obs-count array */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    unit = t / C->T;
	    if (t > pmod->t1 && unit != ubak) {
		if (C->unit_obs[i] > C->Tmax) {
		    C->Tmax = C->unit_obs[i];
		}
		C->unit_obs[++i] = 1;
		C->unit_start[i] = s;
	    } else {
		C->unit_obs[i] += 1;
	    }
	    ubak = unit;
	    s++;
	}
    }

    return 0;
}

/*  Use pooled probit as a starting point; an initial 
    guess of rho is obtained by a rough estimate of the
    sample autocorrelation of (generalized) residuals; 
    sigma_a and beta are then scaled accordingly.
*/

static int params_init_from_pooled (MODEL *pmod, 
				    reprob_container *C,
				    double *theta)
{
    double rho, sigma2_a, e, elag = 0.0;
    double num = 0.0, den = 0.0;
    int k = C->npar - 1;
    int i, t;

    if (pmod->ncoeff != k) {
	fprintf(stderr, "reprobit init: n. of coeffs doesn't match (%d vs %d)\n", 
		pmod->ncoeff, k);
	return E_NONCONF;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (t / C->T > (t-1) / C->T) {
	    /* first obs of a new unit */
	    elag = pmod->uhat[t];
	    elag = na(elag) ? 0.0 : elag;
	    continue;
	}
	e = pmod->uhat[t];
	if (!na(e)) {
	    num += e * elag;
	    den += e * e;
	    elag = e;
	}
    }

    rho = (num > 0.0) ? num/den : 0.01;
    sigma2_a = rho / (1-rho);

#if 1
    fprintf(stderr, "num = %g, den = %g\n", num, den);
    fprintf(stderr, "Initial rho = %g (sigma2_a = %g)\n", rho, sigma2_a);
#endif

    for (i=0; i<k; i++) {
	C->beta->val[i] = theta[i] = pmod->coeff[i] / (1-rho);
    }

    C->scale = theta[k] = log(sigma2_a);

    return 0;
}

static int rep_container_fill (reprob_container *C,
			       MODEL *pmod,
			       DATASET *dset, 
			       int quadpoints,
			       double **ptheta)
{
    gretl_matrix *tmp = NULL;
    int vj, k = C->npar - 1;
    int i, j, t, s;
    double x;
    int err;

    err = reprobit_obs_accounts(C, pmod, dset);
    if (err) {
	return err;
    }

    C->pmod = pmod;
    C->qp = quadpoints;

    C->y = malloc(C->nobs * sizeof *C->y);
    C->X = gretl_matrix_alloc(C->nobs, k);
    C->R = gretl_matrix_alloc(C->nobs, C->qp);

    if (C->y == NULL || C->X == NULL || C->R == NULL) {
	return E_ALLOC;
    }

    *ptheta = malloc(C->npar * sizeof **ptheta);
    if (*ptheta == NULL) {
	return E_ALLOC;
    }

    C->B = gretl_matrix_block_new(&C->ndx, C->nobs, 1,
				  &C->P, C->N, C->qp,
				  &C->lik, C->N, 1,
				  &C->beta, k, 1,
				  &C->nodes, 1, C->qp,
				  &C->wts, C->qp, 1,
				  NULL);
    if (C->B == NULL) {
	return E_ALLOC;
    }

    /* initialize the parameter vector */
    err = params_init_from_pooled(pmod, C, *ptheta);

    /* write the data into C->y and C->X, skipping
       any observations with missing values */
    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    C->y[s] = (dset->Z[C->depvar][t] != 0);
	    for (j=0; j<k; j++) {
		vj = C->list[j+2];
		x = dset->Z[vj][t];
		gretl_matrix_set(C->X, s, j, x);
	    }
	    s++;
	}
    }

    /* form and transcribe the quadrature matrix */
    tmp = gretl_gauss_hermite_matrix_new(C->qp, &err);
    if (!err) {
	for (i=0; i<C->qp; i++) {
	    gretl_vector_set(C->nodes, i, gretl_matrix_get(tmp, i, 0));
	    gretl_vector_set(C->wts, i, gretl_matrix_get(tmp, i, 1));
	}
	gretl_matrix_free(tmp);
    }

#if defined(_OPENMP)
    if (C->nobs * C->qp > 900 && gretl_get_omp_threads() > 1) {
	/* FIXME threshold? */
	C->parallel = 1;
    }
#endif

    if (!C->parallel) {
	C->qi = gretl_matrix_alloc(1, C->qp);
	if (C->qi == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static void update_ndx (reprob_container *C, const double *theta)
{
    int i;

    for (i=0; i<C->npar-1; i++) {
	C->beta->val[i] = theta[i];
    }

    gretl_matrix_multiply(C->X, C->beta, C->ndx);

    C->scale = exp(theta[C->npar-1]/2.0);
}

static double quick_dot_product (const double *x,
				 const double *y,
				 int n)
{
    double ret = 0.0;
    int i;

    for (i=0; i<n; i++) {
	ret += x[i] * y[i];
    }

    return ret;
}

static int reprobit_score (double *theta, double *g, int npar, 
			   BFGS_CRIT_FUNC ll, void *p)
{
    reprob_container *C = (reprob_container *) p;
    gretl_matrix *Q = C->P; /* re-use existing storage */
    const double *nodes = C->nodes->val;
    int i, j, k, t;
    int done = 0;
    int err = 0;

    k = C->npar - 1;
    update_ndx(C, theta);

    /* form the Q and R matrices */

#if defined(_OPENMP)
#pragma omp parallel for private(i, j, t) if (C->parallel)
#endif    
    for (i=0; i<C->N; i++) {
	int Ti = C->unit_obs[i];
	int t0 = C->unit_start[i];
	double x, qij, node, ndxi;
	int sign;

	for (j=0; j<C->qp; j++) {
	    node = C->scale * nodes[j];
	    qij = 1.0;
	    for (t=0; t<Ti; t++) {
		ndxi = C->ndx->val[t0+t];
		sign = C->y[t0+t] ? 1 : -1;
		x = sign * (ndxi + node);
		qij *= normal_cdf(x);
		x = sign * invmills(-x);
		gretl_matrix_set(C->R, t0+t, j, x);
	    }
	    gretl_matrix_set(Q, i, j, qij);
	}
    }

    gretl_matrix_multiply(Q, C->wts, C->lik);

    for (i=0; i<C->npar; i++) {
	g[i] = 0.0;
    }

#if defined(_OPENMP)
    if (C->parallel) {
	double *qi, *all_qi = NULL;
	
#pragma omp parallel private(i, j, t, qi)
	{
	    int nt = omp_get_num_threads();
	    int tid = omp_get_thread_num();

	    if (tid == 0) {
		/* master: allocate workspace for all */
		all_qi = malloc(nt * C->qp * sizeof *all_qi);
		if (all_qi == NULL) {
		    err = E_ALLOC;
		}
	    }

#pragma omp barrier	    
	    qi = err ? NULL : all_qi + tid * C->qp;

#pragma omp for
	    for (i=0; i<C->N; i++) {
		int ii, Ti = C->unit_obs[i];
		int t0 = C->unit_start[i];
		double x, qij, rtj;

		if (qi == NULL) {
		    continue;
		}

		for (ii=0; ii<=k; ii++) {
		    for (j=0; j<C->qp; j++) {
			x = qi[j] = 0.0;
			qij = gretl_matrix_get(Q, i, j);
			if (ii == k) {
			    x = C->scale * nodes[j];
			}
			for (t=0; t<Ti; t++) {
			    if (ii < k) {
				x = gretl_matrix_get(C->X, t0+t, ii);
			    }
			    rtj = gretl_matrix_get(C->R, t0+t, j);
			    qi[j] += x * rtj * qij;
			}
			qi[j] /= C->lik->val[i];
		    }
		    x = quick_dot_product(qi, C->wts->val, C->qp);
#pragma omp atomic		
		    g[ii] += x;
		}
	    }
	}
	free(all_qi);
	done = 1;
    }
#endif /* _OPENMP */

    if (!done) {
	/* single-threaded variant */
	double *qi = C->qi->val;
	
	for (i=0; i<C->N; i++) {
	    int ii, Ti = C->unit_obs[i];
	    int t0 = C->unit_start[i];
	    double x, qij, rtj;

	    for (ii=0; ii<=k; ii++) {
		for (j=0; j<C->qp; j++) {
		    x = qi[j] = 0.0;
		    qij = gretl_matrix_get(Q, i, j);
		    if (ii == k) {
			x = C->scale * nodes[j];
		    }
		    for (t=0; t<Ti; t++) {
			if (ii < k) {
			    x = gretl_matrix_get(C->X, t0+t, ii);
			}
			rtj = gretl_matrix_get(C->R, t0+t, j);
			qi[j] += x * rtj * qij;
		    }
		    qi[j] /= C->lik->val[i];
		}
		g[ii] += quick_dot_product(qi, C->wts->val, C->qp);
	    }
	}
    }

    g[k] /= 2;
    
    return err;
}

static double reprobit_ll (const double *theta, void *p)
{
    reprob_container *C = (reprob_container *) p;
    double x, pij, node;
    int i, j, t;
    int err = 0;

    if (theta[C->npar-1] < -9.0) {
	fprintf(stderr, "reprobit_ll: scale too small\n");
	return NADBL;
    }

    update_ndx(C, theta);
    gretl_matrix_zero(C->P);

#if defined(_OPENMP)
#pragma omp parallel for private(i, j, t, node, pij, x) if (C->parallel)
#endif
    for (i=0; i<C->N; i++) {
	int Ti = C->unit_obs[i];
	int t0 = C->unit_start[i];

	for (j=0; j<C->qp; j++) {
	    node = gretl_vector_get(C->nodes, j);
	    pij = 1.0;
	    for (t=0; t<Ti; t++) {
		x = C->ndx->val[t0+t] + C->scale * node;
		/* the probability */
		pij *= normal_cdf(C->y[t0+t] ? x : -x);
		if (pij < 1.0e-200) {
		    break;
		}
	    }
	    gretl_matrix_set(C->P, i, j, pij);
	}
    }

    err = gretl_matrix_multiply(C->P, C->wts, C->lik);

#if 0
    if (1) {
	/* for analysing the behavior of f(x) */
	static int iter;
	char pname[8];

	sprintf(pname, "P%d.mat", iter);
	gretl_matrix_write_as_text(C->P, pname, 0);
	iter++;
    }
#endif

    if (!err) {
	C->ll = 0.0;
	for (i=0; i<C->N; i++) {
	    C->ll += log(C->lik->val[i]);
	}
    } else {
	C->ll = NADBL;
    }

    return C->ll;
}

static int add_rho_LR_test (MODEL *pmod, double LR)
{
    ModelTest *test = model_test_new(GRETL_TEST_RE);
    int err = 0;

    if (test == NULL) {
	err = E_ALLOC;
    } else {
        model_test_set_teststat(test, GRETL_STAT_LR);
        model_test_set_dfn(test, 1);
        model_test_set_value(test, LR);
        model_test_set_pvalue(test, chisq_cdf_comp(1, LR));
        maybe_add_test_to_model(pmod, test);
    } 

    return err;
}

static int transcribe_reprobit (MODEL *pmod, reprob_container *C,
				double *theta, const DATASET *dset)
{
    gretl_matrix *Hinv;
    double sigma, LR;
    int Tmin = C->nobs;
    int npar = C->npar;
    int k = pmod->ncoeff;
    int i, vi, err = 0;

    Hinv = gretl_zero_matrix_new(npar, npar);

    if (Hinv == NULL) {
	err = E_ALLOC;
    } else {
	err = hessian_from_score(theta, Hinv, reprobit_score, 
				 reprobit_ll, C);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(Hinv);
	if (err) {
	    fprintf(stderr, "hessian_inverse_from_score: failed (err = %d)\n", 
		    err);
	} 
    }

    if (err) {
	/* try generic inverse */
	err = gretl_invert_symmetric_indef_matrix(Hinv);
	if (err) {
	    fprintf(stderr, "hessian_inverse_from_score: failed again (err = %d)\n", 
		    err);
	    gretl_matrix_free(Hinv);
	    Hinv = NULL;
	} else {
	    gretl_model_set_int(pmod, "non-pd-hess", 1);
	}
    }

    if (!err) {
	err = gretl_model_allocate_param_names(pmod, C->npar);
	if (!err) {
	    for (i=0; i<k; i++) {
		vi = pmod->list[i+2];
		gretl_model_set_param_name(pmod, i, dset->varname[vi]);
	    }
	    gretl_model_set_param_name(pmod, k, "lnsigma2");
	}
    }

    if (!err) {
	err = gretl_model_write_coeffs(pmod, theta, C->npar);
    }
    
    if (!err) {
	err = gretl_model_add_hessian_vcv(pmod, Hinv);
    }

    gretl_matrix_free(Hinv);

    if (err) {
	return err;
    }

    /* mask R-squared  */
    pmod->rsq = pmod->adjrsq = NADBL;

    LR = 2.0 * (C->ll - pmod->lnL);
    if (LR < 0) {
	LR = 0;
    }
    add_rho_LR_test(pmod, LR);

    /* set the likelihood-based measures */
    pmod->lnL = C->ll;
    mle_criteria(pmod, 0);

    pmod->sigma = sigma = exp(theta[k]/2);
    pmod->rho = 1 - 1/(1 + sigma * sigma);

    /* note: this may need more thought */
    binary_model_hatvars(pmod, C->ndx, C->y, OPT_E);

    /* panel observation stats */
    for (i=0; i<C->N; i++) {
	if (C->unit_obs[i] < Tmin) {
	    Tmin = C->unit_obs[i];
	}
    }

    gretl_model_set_int(pmod, "quadpoints", C->qp);
    gretl_model_set_int(pmod, "n_included_units", C->N);
    gretl_model_set_int(pmod, "Tmin", Tmin);
    gretl_model_set_int(pmod, "Tmax", C->Tmax);

    pmod->opt |= OPT_E;

    return err;
}

MODEL reprobit_estimate (const int *list, DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    MODEL mod;
    int err = 0;

    /* estimate pooled model as starting point */
    mod = binary_probit(list, dset, OPT_P | OPT_X, prn);

    if ((err = mod.errcode) != 0) {
	pputs(prn, "reprobit: error estimating initial probit\n");
    } else if (mod.list[0] < list[0]) {
	/* something was dropped */
	gretl_errmsg_set("Couldn't estimate initial probit");
    }

    if (!err) {
	/* do the actual reprobit stuff */
	reprob_container *C;
	double *theta = NULL;
	int quadpoints = 32;
	int maxit = libset_get_int(BFGS_MAXITER);
	int fcount = 0;
	gretlopt maxopt;
	
	if (opt & OPT_G) {
	    int qp = get_optval_int(mod.ci, OPT_G, &err);

	    if (qp >= 3 && qp <= 100) {
		quadpoints = qp;
	    }
	}

	if (!err) {
	    C = rep_container_new(mod.list);
	    if (C == NULL) {
		err = E_ALLOC;
	    } else {
		err = rep_container_fill(C, &mod, dset, quadpoints, &theta);
	    }
	}

	if (err) {
	    free(theta);
	    goto bailout;
	}	

	if (libset_get_int(GRETL_OPTIM) == OPTIM_NEWTON) {
	    double crittol = 1.0e-06;
	    double gradtol = 1.0e-05;
	    int quiet = opt & OPT_Q;

	    maxopt = (opt & OPT_V) | OPT_U;
	    err = newton_raphson_max(theta, C->npar, maxit, 
				     crittol, gradtol, &fcount, C_LOGLIK, 
				     reprobit_ll, reprobit_score, NULL, 
				     C, maxopt, quiet ? NULL : prn);
	} else {
	    int gcount = 0;

	    maxopt = opt | OPT_U;
	    err = BFGS_max(theta, C->npar, maxit, 1.0e-9, 
			   &fcount, &gcount, reprobit_ll, C_LOGLIK, 
			   reprobit_score, C, NULL, maxopt, prn);
	}

	if (!err) {
	    transcribe_reprobit(&mod, C, theta, dset);
	}

	rep_container_destroy(C);
	free(theta);
    }

 bailout:

    if (err && mod.errcode == 0) {
	mod.errcode = err;
    }

    if (!err && (opt & OPT_R)) {
	pputs(prn, "Robust variance estimation is not supported for "
	      "random-effects probit.\n");
    }

    return mod;    
}
