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
    int *unit_obs;           /* array of effective unit T's */
    int nobs;                /* total number of observations */
    int qp;                  /* number of quadrature points */

    int *y;	             /* dependent var (0/1) */
    gretl_matrix *X;	     /* main eq. regressors */

    gretl_matrix_block *B;   /* holder for the following */
    gretl_matrix *ndx;       /* index function */
    gretl_matrix *gh_nodes;  /* Gauss-Hermite quadrature nodes */
    gretl_matrix *gh_wts;    /* Gauss-Hermite quadrature weights */
    gretl_matrix *P;         /* probabilities (by individual and qpoints) */
    gretl_matrix *lik;       /* probabilities (by individual) */
    gretl_vector *theta;     /* parameters (including log of variance 
				of individual effect) */
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
	    unit = (int) floor(t / (double) dset->pd);
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

    ubak = -1;
    i = 0;

    /* build the obs-count array */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    unit = (int) floor(t / (double) dset->pd);
	    if (t > pmod->t1 && unit != ubak) {
		C->unit_obs[++i] = 1;
	    } else {
		C->unit_obs[i] += 1;
	    }
	    ubak = unit;
	}
    }

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
    if (C->y == NULL || C->X == NULL) {
	return E_ALLOC;
    }

    C->B = gretl_matrix_block_new(&C->ndx, C->nobs, 1,
				  &C->P, C->N, C->qp,
				  &C->lik, C->N, 1,
				  &C->theta, k, 1,
				  &C->gh_nodes, 1, C->qp,
				  &C->gh_wts, C->qp, 1,
				  NULL);
    if (C->B == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k-1; i++) {
	/* use pooled probit as a starting point */
	C->theta->val[i] = pmod->coeff[i];
    }
    C->theta->val[k-1] = log(0.5); /* rho = 0.5 */

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

static double reprobit_ll (const double *theta, void *p)
{
    reprob_container *C = (reprob_container *) p;
    double lambda, x, pit, a;
    int i, j, t, s, h = C->qp;
    int err = 0;

    /* note: theta and C->theta->val are the same object */
    gretl_matrix_reuse(C->theta, C->npar - 1, 1);
    err = gretl_matrix_multiply(C->X, C->theta, C->ndx);
    gretl_matrix_reuse(C->theta, C->npar, 1);

    lambda = exp(theta[C->npar-1]/2.0);
    gretl_matrix_zero(C->P);

    s = 0;
    for (i=0; i<C->N; i++) {
	int Ti = C->unit_obs[i];

	for (j=0; j<h; j++) {
	    a = gretl_vector_get(C->gh_nodes, j);
	    pit = 1.0;
	    for (t=0; t<Ti; t++) {
		x = gretl_vector_get(C->ndx, s+t) + lambda * a;
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

static void transcribe_reprobit (MODEL *pmod, reprob_container *C)
{
    int Tmin = C->nobs, Tmax = 0;
    int i, t, k = C->npar - 1;

    for (i=0; i<k; i++) {
	pmod->coeff[i] = C->theta->val[i];
	pmod->sderr[i] = 0.0;
    }

    pmod->rsq = pmod->adjrsq = NADBL;
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

    /* check that this is doing the right thing */
    binary_model_hatvars(pmod, C->ndx, C->y, OPT_E);

    gretl_model_set_int(pmod, "n_included_units", C->N);
    gretl_model_set_int(pmod, "Tmin", Tmin);
    gretl_model_set_int(pmod, "Tmax", Tmax);

    pmod->opt |= OPT_E;
}

MODEL reprobit_estimate (const int *list, DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    MODEL mod;
    int *tmplist = NULL;
    int err = 0;

    gretl_model_init(&mod);

    tmplist = gretl_list_copy(list);
    if (tmplist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }    

    /* baseline: estimate via pooled probit */
    mod = binary_probit(tmplist, dset, OPT_A | OPT_P | OPT_X, prn);
    if (mod.errcode) {
	err = mod.errcode;
	fprintf(stderr, "real_panel_model: error %d in intial OLS\n", 
		mod.errcode);
	goto bailout;
    } 

    free(tmplist);

    if (!err) {
	/* do the actual reprobit stuff */
	reprob_container *C;
	int quadpoints = 32;
	int fc, gc;

	C = rep_container_new(list);
	if (C == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	rep_container_fill(C, &mod, dset, quadpoints);

	err = BFGS_max(C->theta->val, C->npar, 100, 1.0e-9, 
		       &fc, &gc, reprobit_ll, C_LOGLIK, 
		       NULL, C, NULL, opt, prn);
	
	if (!err) {
	    pprintf(prn, "estimate of ln(var(u)) = %g\n", 
		    C->theta->val[C->npar-1]);
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
