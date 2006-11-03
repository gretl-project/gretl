/* 
 * Copyright (C) 2004 Riccardo Lucchetti and Allin Cottrell
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
#include "libset.h"

#define ODEBUG 1

#define OPROBIT_TOL 1.0e-12

typedef struct op_container_ op_container;

struct op_container_ {
    int type;         /* model type */
    int *y;           /* dependent variable */
    double **Z;       /* data */
    int *list;
    int M;            /* max of y */
    int t1;           /* beginning of sample */
    int t2;           /* end of sample */
    int nobs;         /* number of observations */
    int nx;           /* number of explanatory variables */
    int k;            /* total number of parameters */
    double *ndx;      /* index variable */
    double *dP;       /* probabilities */
    MODEL *pmod;      /* model struct, initially containing OLS */
    int ascore;       /* 0 for numerical score, 1 analytical */
    double **G;       /* score matrix by observation */
    double *g;        /* total score vector */
};

static double distfunc (double x, int type)
{
    switch (type) {
    case PROBIT:
	return normal_cdf(x);
    case LOGIT:
	return 1.0 / (1.0 + exp(-x));
    default:
	return NADBL;
    }
}

static double densfunc (double x, int type)
{
    double tmp;

    switch (type) {
    case PROBIT:
	return normal_pdf(x);
    case LOGIT:
	tmp = 1.0 + exp(-x);
	return (tmp - 1.0) / (tmp * tmp);
    default:
	return NADBL;
    }
}

static void op_container_destroy (op_container *OC)
{
    free(OC->y);
    free(OC->ndx);
    free(OC->dP);
    free(OC->list);
    doubles_array_free(OC->G, OC->k);
    free(OC->g);

    free(OC);
}

static op_container *op_container_new (int type, double **Z, MODEL *pmod, 
				       int ascore)
{
    int i, t, nx, vy;
    int nobs = pmod->nobs;

    op_container *OC = malloc(sizeof *OC);
    if (OC == NULL) {
	return NULL;
    }

    OC->type = type;

    OC->Z = Z;
    OC->pmod = pmod;
    OC->t1 = pmod->t1;
    OC->t2 = pmod->t2;
    OC->nobs = nobs;
    OC->k = pmod->ncoeff;

    vy = pmod->list[1];
    OC->M = (int) gretl_max(OC->t1, OC->t2, Z[vy]);
    nx = OC->k - (OC->M - 1);
    OC->nx = nx;
    OC->ascore = ascore;

    OC->y = NULL;
    OC->ndx = NULL;
    OC->dP = NULL;
    OC->list = NULL;
    OC->G = NULL;
    OC->g = NULL;

    OC->y = malloc(nobs * sizeof *OC->y);
    OC->ndx = malloc(nobs * sizeof *OC->ndx);
    OC->dP = malloc(nobs * sizeof *OC->dP);
    OC->list = gretl_list_new(1 + nx);
    OC->G = doubles_array_new(OC->k, nobs);
    OC->g = malloc(OC->k * sizeof *OC->g);

    if (OC->y == NULL || OC->ndx == NULL || 
	OC->dP == NULL || OC->list == NULL ||
	OC->G == NULL || OC->g == NULL) {
	op_container_destroy(OC);
	return NULL;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    OC->y[i++] = (int) Z[vy][t];
	}
    }

#if ODEBUG
    fprintf(stderr, "nobs = %d\n", OC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", OC->t1, OC->t2);
    fprintf(stderr, "k = %d\n", OC->k);
    fprintf(stderr, "nx = %d\n", OC->nx);
    fprintf(stderr, "Max(y) = %d\n", OC->M);
#endif

    OC->list[1] = vy;
    for (i=1; i<=nx; i++) {
	OC->list[i+1] = pmod->list[i+1];
    }

#if ODEBUG > 1
    printlist(OC->list, "list, in op_container_new");
#endif

    return OC;
}

static int compute_probs (const double *theta, op_container *OC)
{
    double m0, m1, ystar0 = 0.0, ystar1 = 0.0;
    int M = OC->M;
    int nx = OC->nx;
    double P0, P1, dP;
    int i, t, yt, v;
    int type = OC->type;

    if (OC->ascore) {
	for (i=0; i<OC->k; i++) {
	    OC->g[i] = 0.0;
	}
    }

    for (t=0; t<OC->nobs; t++) {
	yt = OC->y[t];
	if (yt == 0) {
	    ystar1 = OC->ndx[t];
	} else if (yt == 1) {
	    m1 = theta[nx];
	    ystar0 = OC->ndx[t];
	    ystar1 = OC->ndx[t] + m1;
	} else {
	    m0 = theta[nx + yt - 2];
	    ystar0 = OC->ndx[t] + m0;
	    if (yt < M) {
		m1 = theta[nx + yt - 1];
		ystar1 = OC->ndx[t] + m1;
	    }
	} 

	P0 = (yt == 0)? 0.0 : distfunc(ystar0, type);
	P1 = (yt == M)? 1.0 : distfunc(ystar1, type);
	dP = P1 - P0;

	if (dP > 1.0e-09) {
	    OC->dP[t] = dP;
	} else {
	    fprintf(stderr, "t:%4d y=%d, ndx = %10.6f, dP = %9.7f\n", 
		    t, yt, OC->ndx[t], dP);
	    return 1;
	} 

	if (OC->ascore) {
	    double mills0 = (yt == 0)? 0.0 : densfunc(ystar0, type) / dP;
	    double mills1 = (yt == M)? 0.0 : densfunc(ystar1, type) / dP;
	    double dm = mills1 - mills0;

	    for (i=0; i<nx; i++) {
		v = OC->list[i+2];
		OC->G[i][t] = -dm * OC->Z[v][t];
		OC->g[i] += OC->G[i][t];
	    }

	    for (i=nx; i<OC->k; i++) {
		OC->G[i][t] = 0.0;
		if (i == (nx + yt - 2)) {
		    OC->G[i][t] = -mills0;
		    OC->g[i] += OC->G[i][t];
		}
		if (i == (nx + yt - 1)) {
		    OC->G[i][t] = mills1;
		    OC->g[i] += OC->G[i][t];
		}
	    }
	}
    }

    return 0;
}

static int bad_cutpoints (const double *theta, const op_container *OC)
{
    int i;
    
    if (theta[OC->nx] <= 0.0) {
	return 1;
    }
    
    for (i=OC->nx+1; i<OC->k; i++) {
	if (theta[i] <= theta[i-1]) {
	    return 1;
	}
    }
	
    return 0;
}

static double op_loglik (const double *theta, void *ptr)
{
    double x, ll = 0.0;
    int err;
    int i, s, t, v;
    
    op_container *OC = (op_container *) ptr;
    
    if (bad_cutpoints(theta, OC)) {
#if ODEBUG > 1
	fputs("Non-increasing cutpoint!\n", stderr);
#endif
	return NADBL;
    }

    s = 0;
    for (t=OC->t1; t<=OC->t2; t++) {
	if (na(OC->pmod->uhat[t])) {
	    continue;
	}
	x = 0.0;
	for (i=0; i<OC->nx; i++) {
	    v = OC->list[i+2];
	    x -= theta[i] * OC->Z[v][t];
	}
	OC->ndx[s++] = x;
    }
    
    err = compute_probs(theta, OC);
    if (err) {
	ll = NADBL;
    } else {
	for (t=0; t<OC->nobs; t++) {
	    ll += log(OC->dP[t]);
	}
    }

#if ODEBUG > 1
    fprintf(stderr, "ll = %16.10f\n", ll);
#endif

    return ll;
}

static int op_score (double *theta, double *s, int npar, BFGS_LL_FUNC ll, 
		     void *ptr)
{
    op_container *OC = (op_container *) ptr;
    int i;

    for (i=0; i<npar; i++) {
	s[i] = OC->g[i];
    }

    return 1;
}

static int ihess_from_ascore (double *theta, op_container *OC, gretl_matrix *inH) 
{
    int i, j, err, k = OC->k;
    double smal = 1.0e-07;  /* "small" is some sort of macro on win32 */
    double smal2 = 2.0 * smal;
    double x, ll;

    double *g0 = malloc(k * sizeof *g0);

    if (g0 == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	theta[i] -= smal;
	ll = op_loglik(theta, OC);
	for (j=0; j<k; j++) {
	    g0[j] = OC->g[j];
	}
	theta[i] += smal2;
	ll = op_loglik(theta, OC);
	for (j=0; j<k; j++) {
	    x = (OC->g[j] - g0[j]) / smal2;
	    gretl_matrix_set(inH, i, j, -x);
	}
	theta[i] -= smal;
    }

    /* make symmetric */

    for (i=0; i<k; i++) {
	for (j=i+1; j<k; j++) {
	    x = 0.5 * (gretl_matrix_get(inH, i, j) + gretl_matrix_get(inH, j, i));
	    gretl_matrix_set(inH, i, j, x);
	    gretl_matrix_set(inH, j, i, x);
	}
    }

    free(g0);

    err = gretl_invert_symmetric_matrix(inH);

    return err;
}

static int 
numerical_ihess (op_container *OC, double *theta, gretl_matrix *invH)
{
    double vij, *V;
    int i, j, k, npar = OC->k;

    if (invH == NULL) {
	return E_ALLOC;
    }

    V = numerical_hessian(theta, npar, op_loglik, OC);
    if (V == NULL) {
	return E_ALLOC;
    }

    k = 0;
    for (i=0; i<npar; i++) {
	for (j=i; j<npar; j++) {
	    vij = V[k++];
	    gretl_matrix_set(invH, i, j, vij);
	    if (i != j) {
		gretl_matrix_set(invH, j, i, vij);
	    }
	}
    }

    free(V);

    return 0;
}

static void fill_model (MODEL *pmod, const DATAINFO *pdinfo, 
			op_container *OC, double *theta, 
			gretl_matrix *invH)
{
    int npar = OC->k;
    int nx = OC->nx;
    double x;
    int i, j, k, t, v;

    pmod->t1 = OC->t1;
    pmod->t2 = OC->t2;
    pmod->nobs = OC->nobs;
    pmod->ci = OC->type;
    gretl_model_set_int(pmod, "ordered", 1);

    pmod->ncoeff = npar;

    if (pmod->vcv == NULL) {
	pmod->vcv = malloc(npar * (npar+1) / 2 * sizeof *pmod->vcv);
	if (pmod->vcv == NULL) {
	    pmod->errcode = E_ALLOC;
	    return;
	}
    }

    k = 0;
    for (i=0; i<npar; i++) {
	pmod->coeff[i] = theta[i];
#if ODEBUG > 1
	fprintf(stderr,"theta[%d] = %12.8f\n", i, theta[i]);
#endif
	for (j=0; j<=i; j++) {
	    x = gretl_matrix_get(invH, i, j);
	    pmod->vcv[ijton(i,j,npar)] = x;
	    if (i == j) {
		pmod->sderr[i] = sqrt(x);
	    }
	}	    
    }

    for (t=OC->t1; t<=OC->t2; t++) {
	if (na(OC->pmod->uhat[t])) {
	    continue;
	}
	x = 0.0;
	for (i=0; i<OC->nx; i++) {
	    v = OC->list[i+2];
	    x -= theta[i] * OC->Z[v][t];
	}
	/* - \hat{beta}'x */
	pmod->yhat[t] = x;
	pmod->uhat[t] = NADBL; /* can we do something sensible here? */
    }

    pmod->lnL = op_loglik(theta, OC);

    gretl_model_allocate_params(pmod, npar);
    if (pmod->errcode == 0) {
	for (i=0; i<nx; i++) {
	    strcpy(pmod->params[i], pdinfo->varname[OC->list[i+2]]);
	}
	k = 1;
	for (i=nx; i<npar; i++) {
	    sprintf(pmod->params[i], "cut(%d)", k++);
	}
    }
}

/* Main ordered estimation function */

static int do_ordered (int ci, double **Z, DATAINFO *pdinfo, MODEL *pmod,
		       PRN *prn)
{
    int i, npar;
    gretl_matrix *iH = NULL;
    double *theta = NULL;
    int err;

    /* BFGS apparatus */
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;
    int ascore = 1;

    pprintf(prn, "Analytical score = %d\n", ascore);

    op_container *OC = op_container_new(ci, Z, pmod, ascore);
    if (OC == NULL) {
	return E_ALLOC;
    }

    npar = OC->k;

    theta = malloc(npar * sizeof *theta);
    if (theta == NULL) {
	op_container_destroy(OC);
	return E_ALLOC;
    }

    for (i=0; i<OC->nx; i++) {
	theta[i] = pmod->coeff[i];
    }

    /* 
       the following is very heuristic, has no sound theoretical 
       justification but seems to work in practice. So much
       for the Scientific Method (TM) !!!
    */
 
    for (i=OC->nx; i<npar; i++) {
	theta[i] = 0.5 * pmod->coeff[i];
    }


#if ODEBUG
    fputs("Starting values:\n", stderr);
    for (i=0; i<npar; i++) {
	fprintf(stderr, "theta[%d] = %12.8f\n", i, theta[i]);
    }
    fprintf(stderr, "\ninitial loglikelihood = %12.8f\n", 
	    op_loglik(theta, OC));
#endif

    err = BFGS_max(theta, npar, maxit, OPROBIT_TOL, 
		   &fncount, &grcount, op_loglik, 
		   (OC->ascore)? op_score : NULL, 
		   OC, (prn != NULL)? OPT_V : OPT_NONE,
		   prn);

    if (err) {
	return err;
    } else {
	fprintf(stderr, "Number of iterations = %d (%d)\n", fncount, grcount);
    }

    for (i=0; i<npar; i++) {
	pmod->coeff[i] = theta[i];
    }

    iH = gretl_matrix_alloc(npar, npar);
    if(ascore) {
	err = ihess_from_ascore(theta, OC, iH);
    } else {
	err = numerical_ihess(OC, theta, iH);
    }

#if ODEBUG > 1
    gretl_matrix_print(iH, "Inverse Hessian");
#endif

    fill_model(pmod, pdinfo, OC, theta, iH);

    free(theta);
    gretl_matrix_free(iH);
    op_container_destroy(OC);

    return 0;
}

/* the driver function for the plugin: note, if prn is non-NULL,
   the verbose option has been selected
*/

MODEL ordered_estimate (const int *list, int ci, double ***pZ, DATAINFO *pdinfo,
			PRN *prn) 
{
    MODEL model;
    int *newlist = NULL;
    int *dumlist = NULL;
    int i, k, nv;

    gretl_model_init(&model);

    dumlist = malloc(2 * sizeof *dumlist);
    if (dumlist == NULL) {
	model.errcode = E_ALLOC;
	return model;
    }

    dumlist[0] = 1;
    dumlist[1] = list[1];

    model.errcode = list_dumgenr(&dumlist, pZ, pdinfo, OPT_F);
    if (model.errcode) {
	free(dumlist);
	return model;
    }

    /* we drop 2 dummies, not just the first one */
    nv = list[0] + dumlist[0] - 1;

    newlist = gretl_list_new(nv);
    if (newlist == NULL) {
	model.errcode = E_ALLOC;
	goto bailout;
    }

    k = 1;
    for (i=1; i<=list[0]; i++) {
	newlist[k++] = list[i];
    }
    for (i=2; i<=dumlist[0]; i++) {
	newlist[k++] = dumlist[i];
    }

    /* run initial OLS */
    model = lsq(newlist, pZ, pdinfo, OLS, OPT_A);
    if (model.errcode) {
	goto bailout;
    }

#if ODEBUG > 1
    pprintf(prn, "oprobit_estimate: initial OLS\n");
    printmodel(&model, pdinfo, OPT_NONE, prn);
#endif

    /* do the actual ordered probit analysis */
    if (model.errcode == 0) {
	model.errcode = do_ordered(ci, *pZ, pdinfo, &model, prn);
    }

 bailout:

    free(dumlist);
    free(newlist);

    return model;
}
