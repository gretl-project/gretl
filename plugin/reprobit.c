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
#include "gretl_model.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"

typedef struct reprob_container_ reprob_container;

struct reprob_container_ {
    int *list;               /* model specification */
    int depvar;		     /* location of y in array Z */
    int npar;		     /* no. of parameters */
    double ll;		     /* log-likelihood */
    MODEL *pmod;             /* pointer to model struct */

    int N;	             /* number of included units */
    int Tmax;                /* maximum tim-series observations */
    int *unit_obs;           /* array of effective unit T's */
    int nobs;                /* total number of observations */
    int qp;                  /* number of quadrature points */

    int *y;	             /* dependent var (0/1) */
    gretl_matrix *X;	     /* main eq. regressors */
    gretl_matrix *R;         /* holds inverse Mills ratios */

    gretl_matrix_block *B;   /* holder for the following */
    gretl_matrix *ndx;       /* index function */
    gretl_matrix *gh_nodes;  /* Gauss-Hermite quadrature nodes */
    gretl_matrix *gh_wts;    /* Gauss-Hermite quadrature weights */
    gretl_matrix *P;         /* probabilities (by individual and qpoints) */
    gretl_matrix *lik;       /* probabilities (by individual) */
    gretl_vector *theta;     /* parameters (including log of variance 
				of individual effect) */
    gretl_matrix *qi;        /* storage for score caculation */
};

reprob_container *rep_container_new (const int *list)
{
    reprob_container *C = malloc(sizeof *C);

    if (C != NULL) {
	C->list = gretl_list_copy(list);
	if (C->list == NULL) {
	    free(C);
	    return NULL;
	}

	C->depvar = C->list[1];
	C->npar = list[0];
	C->ll = NADBL;
	C->N = 0;
	C->nobs = 0;

	C->unit_obs = NULL;
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
	free(C->list);
	free(C->y);
	free(C->unit_obs);
	gretl_matrix_free(C->X);
	gretl_matrix_free(C->R);
	gretl_matrix_block_destroy(C->B);
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
    int i, t;

    C->nobs = pmod->nobs;

    /* count the included units */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    unit = t / dset->pd;
	    if (unit != ubak) {
		C->N += 1;
	    }
	    ubak = unit;
	}
    }

    C->unit_obs = malloc(C->N * sizeof *C->unit_obs);
    if (C->unit_obs == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<C->N; i++) {
	C->unit_obs[i] = 0;
    }

    C->Tmax = 0;
    ubak = -1;
    i = 0;

    /* build the obs-count array */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    unit = t / dset->pd;
	    if (t > pmod->t1 && unit != ubak) {
		if (C->unit_obs[i] > C->Tmax) {
		    C->Tmax = C->unit_obs[i];
		}
		C->unit_obs[++i] = 1;
	    } else {
		C->unit_obs[i] += 1;
	    }
	    ubak = unit;
	}
    }

    return 0;
}

static int params_init_from_pooled (MODEL *pmod, 
				    gretl_vector *par)
{
    /* 
       use pooled probit as a starting point; an initial 
       guess of rho is obtained by a rough estimate of the
       sample autocorrelation of (generalized) residuals. 
       sigma_a and beta are scaled accordingly
     */

    int i, s, t;
    int k = par->rows - 1;

    if (pmod->ncoeff != k) {
	fprintf(stderr, "What??? %d vs %d\n", pmod->ncoeff, k);
	return E_NONCONF;
    }

    s = 0;
    double rho, sigma2_a, e, elag = 0.0;
    double num = 0.0;
    double den = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	e = pmod->uhat[t];
	if (!xna(e)) {
	    num += e * elag;
	    den += e * e;
	    elag = e;
	}
    }

    rho = (num>0) ? num/den : 0.001;
    sigma2_a = rho / (1-rho);

#if 0
    fprintf(stderr, "num = %g, den = %g\n", num, den);
    fprintf(stderr, "Initial rho = %g (sigma2_a = %g)\n", rho, sigma2_a);
#endif

    for (i=0; i<k; i++) {
	par->val[i] = pmod->coeff[i] / (1-rho);
    }

    par->val[k] = log(sigma2_a);
    return 0;
}

static int rep_container_fill (reprob_container *C,
			       MODEL *pmod,
			       DATASET *dset, 
			       int quadpoints)
{
    gretl_matrix *tmp = NULL;
    int vj, k = C->npar;
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
    C->X = gretl_matrix_alloc(C->nobs, k-1);
    C->R = gretl_matrix_alloc(C->nobs, C->qp);

    if (C->y == NULL || C->X == NULL || C->R == NULL) {
	return E_ALLOC;
    }

    C->B = gretl_matrix_block_new(&C->ndx, C->nobs, 1,
				  &C->P, C->N, C->qp,
				  &C->lik, C->N, 1,
				  &C->theta, k, 1,
				  &C->gh_nodes, 1, C->qp,
				  &C->gh_wts, C->qp, 1,
				  &C->qi, 1, C->qp,
				  NULL);
    if (C->B == NULL) {
	return E_ALLOC;
    }

    /* initialize the parameters vector */

    err = params_init_from_pooled(pmod, C->theta);

    /* write the data into C->y and C->X, skipping
       any observations with missing values */

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    C->y[s] = (dset->Z[C->depvar][t] != 0);
	    for (j=0; j<k-1; j++) {
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
	    gretl_vector_set(C->gh_nodes, i, gretl_matrix_get(tmp, i, 0));
	    gretl_vector_set(C->gh_wts, i, gretl_matrix_get(tmp, i, 1));
	}
	gretl_matrix_free(tmp);
    }

    return err;
}

static void update_ndx (reprob_container *C, double *scale)
{
    gretl_matrix_reuse(C->theta, C->npar-1, 1);
    gretl_matrix_multiply(C->X, C->theta, C->ndx);
    gretl_matrix_reuse(C->theta, C->npar, 1);

    *scale = exp(C->theta->val[C->npar-1]/2.0);
}

static int reprobit_score (double *theta, double *g, int npar, 
			   BFGS_CRIT_FUNC ll, void *p)
{
    reprob_container *C = (reprob_container *) p;
    gretl_matrix *Q, *qi;
    const double *nodes = C->gh_nodes->val;
    double scale, x, qij, ndxi, node;
    int i, j, k, s, t, h = C->qp;
    int sign = 1;
    int err = 0;

    k = C->npar - 1;
    Q = C->P;   /* re-use existing storage of right size */
    qi = C->qi; /* reduce verbosity below */ 

    update_ndx(C, &scale);

    /* form the Q and R matrices */
    s = 0;
    for (i=0; i<C->N; i++) {
	int Ti = C->unit_obs[i];

	for (j=0; j<h; j++) {
	    node = scale * nodes[j];
	    qij = 1.0;
	    for (t=0; t<Ti; t++) {
		ndxi = C->ndx->val[s+t];
		sign = C->y[s+t] ? 1 : -1;
		x = sign * (ndxi + node);
		qij *= normal_cdf(x);
		x = sign * invmills(-x);
		gretl_matrix_set(C->R, s+t, j, x);
	    }
	    gretl_matrix_set(Q, i, j, qij);
	}
	s += Ti;
    }

    gretl_matrix_multiply(Q, C->gh_wts, C->lik);

    for (i=0; i<C->npar; i++) {
	g[i] = 0.0;
    }

    s = 0;
    for (i=0; i<C->N; i++) {
	int ii, Ti = C->unit_obs[i];
	double rtj, tmp;

	for (ii=0; ii<=k; ii++) {
	    for (j=0; j<h; j++) {
		x = qi->val[j] = 0.0;
		qij = gretl_matrix_get(Q, i, j);
		if (ii == k) {
		    x = scale * nodes[j];
		}
		for (t=0; t<Ti; t++) {
                    if (ii < k) {
		        x = gretl_matrix_get(C->X, s+t, ii);
                    }
 		    rtj = gretl_matrix_get(C->R, s+t, j);
		    qi->val[j] += x * rtj * qij;
		}
		qi->val[j] /= C->lik->val[i];
	    }
            tmp = gretl_vector_dot_product(qi, C->gh_wts, &err);
	    g[ii] += (ii < k)? tmp : tmp * scale;
 	}
	s += Ti;
    }
    
    return err;
}

static double reprobit_ll (const double *theta, void *p)
{
    reprob_container *C = (reprob_container *) p;
    double scale, x, pit, node;
    int i, j, t, s, h = C->qp;
    int err = 0;

    update_ndx(C, &scale);
    gretl_matrix_zero(C->P);

    s = 0;
    for (i=0; i<C->N; i++) {
	int Ti = C->unit_obs[i];

	for (j=0; j<h; j++) {
	    node = gretl_vector_get(C->gh_nodes, j);
	    pit = 1.0;
	    for (t=0; t<Ti; t++) {
		x = C->ndx->val[s+t] + scale * node;
		/* the probability */
		pit *= normal_cdf(C->y[s+t] ? x : -x);
		if (pit < 1.0e-30) {
		    break;
		}
	    }
	    gretl_matrix_set(C->P, i, j, pit);
	}
	s += Ti;
    }	    

    err = gretl_matrix_multiply(C->P, C->gh_wts, C->lik);

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

static int transcribe_reprobit (MODEL *pmod, reprob_container *C)
{
    gretl_matrix *Hinv;
    int Tmin = C->nobs, Tmax = 0;
    int i, k = pmod->ncoeff;
    double sigma, LR;
    int err = 0;

    Hinv = hessian_inverse_from_score(C->theta->val, C->npar,
				      reprobit_score,
				      reprobit_ll,
				      C, &err);

    for (i=0; i<k; i++) {
	pmod->coeff[i] = C->theta->val[i];
	if (Hinv != NULL) {
	    pmod->sderr[i] = sqrt(gretl_matrix_get(Hinv, i, i));
	} else {
	    pmod->sderr[i] = NADBL;
	}
    }

    gretl_matrix_free(Hinv);

    pmod->rsq = pmod->adjrsq = NADBL;
    LR = 2.0 * (C->ll - pmod->lnL);
    /* LR test for var(u) = 0 */
    fprintf(stderr, "LR = %g\n", LR);
    pmod->lnL = C->ll;
    mle_criteria(pmod, 1);

    for (i=0; i<C->N; i++) {
	if (C->unit_obs[i] < Tmin) {
	    Tmin = C->unit_obs[i];
	}
	if (C->unit_obs[i] > Tmax) {
	    Tmax = C->unit_obs[i];
	}
    }

    pmod->sigma = sigma = exp(C->theta->val[k]/2);
    pmod->rho = 1 - 1/(1 + sigma * sigma);

    /* check that this is doing the right thing */
    binary_model_hatvars(pmod, C->ndx, C->y, OPT_E);

    gretl_model_set_int(pmod, "n_included_units", C->N);
    gretl_model_set_int(pmod, "Tmin", Tmin);
    gretl_model_set_int(pmod, "Tmax", Tmax);

    pmod->opt |= OPT_E;

    return err;
}

static int reprobit_init (MODEL *pmod, const int *list,
			  DATASET *dset, PRN *prn)
{
    int *tmplist = gretl_list_copy(list);

    if (tmplist == NULL) {
	gretl_model_init(pmod);
	pmod->errcode = E_ALLOC;
    } else {
	*pmod = binary_probit(tmplist, dset, OPT_A | OPT_P | OPT_X, prn);
	if (pmod->errcode) {
	    fprintf(stderr, "reprobit_estimate: error %d in "
		    "initial probit\n", pmod->errcode);
	}
	free(tmplist);
    }

    return pmod->errcode;
}

MODEL reprobit_estimate (const int *list, DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    MODEL mod;
    int err;

    err = reprobit_init(&mod, list, dset, prn);

    if (!err) {
	/* do the actual reprobit stuff */
	reprob_container *C;
	int fc, gc;
	int quadpoints;
	
	if (opt & OPT_G) {
	    quadpoints = get_optval_int(mod.ci, OPT_G, &err);
	} else {
	    quadpoints = 32;
	}

	C = rep_container_new(list);
	if (C == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	rep_container_fill(C, &mod, dset, quadpoints);
#if 1
	err = BFGS_max(C->theta->val, C->npar, 100, 1.0e-9, 
		       &fc, &gc, reprobit_ll, C_LOGLIK, 
		       reprobit_score, C, NULL, opt, prn);
#else
	double crittol = 1.0e-06;
	double gradtol = 1.0e-05;
	gretlopt maxopt = opt & OPT_V;
	int quiet = opt & OPT_Q;
	int maxit = 1000;

	err = newton_raphson_max(C->theta->val, C->npar, maxit, 
				 crittol, gradtol, &fc, C_LOGLIK, 
				 reprobit_ll, reprobit_score, NULL, 
				 C, maxopt, quiet ? NULL : prn);
#endif
	
	if (!err) {
	    transcribe_reprobit(&mod, C);
	}

	rep_container_destroy(C);
    }

 bailout:

    if (err && mod.errcode == 0) {
	mod.errcode = err;
    }

    return mod;    
}
