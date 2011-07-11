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

#include "libgretl.h"
#include "libset.h"
#include "gretl_bfgs.h"

#define INTDEBUG 0

enum {
    INT_LOW,   /* no lower bound */
    INT_MID,   /* both bounds */
    INT_HIGH,  /* no upper bound */
    INT_POINT  /* point value */
};

typedef struct int_container_ int_container;

struct int_container_ {
    MODEL *pmod;      /* model struct, initially containing OLS */
    int hiv, lov;     /* variable numbers for limit series */
    double ll;        /* loglikelihood */
    double *dspace;   /* workspace */
    double *hi, *lo;  /* limit series for dependent variable */
    int *obstype;     /* dependent variable classifier */
    int typecount[4]; /* no. of obs by category */
    double **X;       /* regressors */
    int *list;        /* dependent vars plus regular regressors */
    int t1;           /* beginning of sample */
    int t2;           /* end of sample */
    int nobs;         /* number of observations */
    int nx;           /* number of explanatory variables */
    int k;            /* total number of parameters */
    double *theta;    /* real parameter estimates */
    double *ndx;      /* index variable */
    double *uhat;     /* generalized residuals */
    double *dP;       /* probabilities */
    double *f0;       /* normalized density at min */
    double *f1;       /* normalized density at max */
    double **G;       /* score matrix by observation */
    double *g;        /* total score vector */
};

static void int_container_destroy (int_container *IC)
{
    free(IC->dspace);

    free(IC->theta);
    free(IC->g);
    
    free(IC->obstype);
    free(IC->list);

    doubles_array_free(IC->X, IC->nx);
    doubles_array_free(IC->G, IC->k);

    free(IC);
}

static int_container *int_container_new (int *list, DATASET *dset, 
					 MODEL *mod)
{
    int_container *IC;
    int nobs, i, s, t, vlo = list[1], vhi = list[2];
    double x0, x1;

    IC = malloc(sizeof *IC);
    if (IC == NULL) {
	return NULL;
    }

    IC->lov = vlo;
    IC->hiv = vhi;
    IC->pmod = mod;

    IC->t1 = mod->t1;
    IC->t2 = mod->t2;
    nobs = IC->nobs = mod->nobs;
    IC->nx = mod->ncoeff;
    IC->k = mod->ncoeff + 1; /* beta + sigma */
    for (i=0; i<4; i++) {
	IC->typecount[i] = 0;
    }

    IC->dspace = NULL;

    IC->theta = IC->g = NULL;
    IC->X = IC->G = NULL;
    
    IC->obstype = NULL;
    IC->list = NULL;

    /* doubles arrays, length nobs */
    IC->dspace = malloc(7 * nobs * sizeof *IC->dspace);

    /* doubles arrays, length k */
    IC->theta = malloc(IC->k * sizeof *IC->theta);
    IC->g = malloc(IC->k * sizeof *IC->g);

    /* two-dimensional doubles arrays */
    IC->X = doubles_array_new(IC->nx, nobs);
    IC->G = doubles_array_new(IC->k, nobs);

    /* int arrays */
    IC->obstype = malloc(nobs * sizeof *IC->obstype);
    IC->list = gretl_list_new(2 + IC->nx);

    if (IC->dspace == NULL ||
	IC->theta == NULL || IC->g == NULL || 
	IC->X == NULL || IC->G == NULL || 
	IC->obstype == NULL || IC->list == NULL) {
	int_container_destroy(IC);
	return NULL;
    }

    IC->hi = IC->dspace;
    IC->lo = IC->hi + nobs;
    IC->ndx = IC->lo + nobs;
    IC->uhat = IC->ndx + nobs;
    IC->dP = IC->uhat + nobs;
    IC->f0 = IC->dP + nobs;
    IC->f1 = IC->f0 + nobs;

    s = 0;
    for (t=mod->t1; t<=mod->t2; t++) {
	if (!na(mod->uhat[t])) {

	    x0 = dset->Z[vlo][t];
	    x1 = dset->Z[vhi][t];
	    IC->lo[s] = x0;
	    IC->hi[s] = x1;

	    if (na(x0)) {
		IC->obstype[s] = INT_LOW;
	    } else if (na(x1)) {
		IC->obstype[s] = INT_HIGH;
	    } else if (x0 == x1) {
		IC->obstype[s] = INT_POINT;
	    } else {
		IC->obstype[s] = INT_MID;
	    }

	    for (i=0; i<IC->nx; i++) {
		IC->X[i][s] = dset->Z[list[i+3]][t];
	    }

	    s++;
	}
    }

    IC->list[1] = vlo;
    IC->list[2] = vhi;
    for (i=0; i<IC->nx; i++) {
	IC->list[i+3] = mod->list[i+2];
    }

    for (i=0; i<IC->nx; i++) {
	IC->theta[i] = mod->coeff[i];
    }
    IC->theta[i] = log(mod->sigma);

#if INTDEBUG > 2
    fprintf(stderr, "nobs = %d\n", IC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", IC->t1, IC->t2);
    fprintf(stderr, "k = %d\n", IC->k);
    fprintf(stderr, "nx = %d\n", IC->nx);
    for (i=0; i<=IC->nx; i++) {
	fprintf(stderr, "theta[%2d] = %10.6f\n", i, IC->theta[i]);
    }
#endif

    return IC;
}

static int create_midpoint_y (int *list, DATASET *dset, 
			      int **initlist)
{
    int mpy = dset->v;
    double *lo, *hi, *mid;
    double x0, x1;
    int i, t, err;

    err = dataset_add_series(1, dset);
    if (err) {
	return err;
    }

    lo = dset->Z[list[1]];
    hi = dset->Z[list[2]];
    mid = dset->Z[mpy];

    for (t=dset->t1; t<=dset->t2 && !err; t++) {
	x0 = lo[t];
	x1 = hi[t];

	if (na(x0)) {
	    mid[t] = x1;
	} else if (na(x1)) {
	    mid[t] = x0;
	} else if (x0 > x1) {
	    gretl_errmsg_sprintf(_("Obs %d: lower bound (%g) "
				   "exceeds upper (%g)"), t + 1,
				 x0, x1);
	    err = E_DATA;
	} else {
	    mid[t] = 0.5 * (x0 + x1);
	}
    }

    if (err) {
	return err;
    }

    *initlist = gretl_list_new(list[0] - 1);

    if (*initlist == NULL) {
	err = E_ALLOC;
    } else {
	(*initlist)[1] = mpy;
	for (i=3; i<=list[0]; i++) {
	    (*initlist)[i-1] = list[i];
	}

#if INTDEBUG > 1
	for (t=dset->t1; t<=dset->t2; t++) {
	    fprintf(stderr, "%2d: %g\n", t, mid[t]);
	}
#endif
    }

    return err;
}

static void loglik_prelim (const double *theta, int_container *IC)
{
    int i, t, k = IC->k;
    double ndxt, z0, z1, x0, x1;
    double sigma = exp(theta[k-1]);

    for (t=0; t<IC->nobs; t++) {
	ndxt = 0.0;
	for (i=0; i<IC->nx; i++) {
	    ndxt += theta[i] * IC->X[i][t];
	}

	IC->ndx[t] = ndxt;

	x0 = IC->lo[t];
	x1 = IC->hi[t];

	switch (IC->obstype[t]) {
	case INT_LOW:
	    z1 = (x1 - ndxt)/sigma;
	    IC->dP[t] = normal_cdf(z1);
	    IC->f0[t] = 0;
	    IC->f1[t] = normal_pdf(z1) / IC->dP[t];
	    break;
	case INT_HIGH:
	    z0 = (x0 - ndxt)/sigma;
	    IC->dP[t] = normal_cdf_comp(z0);
	    IC->f0[t] = normal_pdf(z0) / IC->dP[t];
	    IC->f1[t] = 0;
	    break;
	case INT_MID:
	    z0 = (x0 - ndxt)/sigma;
	    z1 = (x1 - ndxt)/sigma;
	    IC->dP[t] = normal_cdf(z1) - normal_cdf(z0);
	    IC->f0[t] = normal_pdf(z0) / IC->dP[t];
	    IC->f1[t] = normal_pdf(z1) / IC->dP[t];
	    break;
	case INT_POINT:
	    z0 = (x0 - ndxt)/sigma;
	    IC->dP[t] = normal_pdf(z0)/sigma;
	    IC->f0[t] = IC->f1[t] = 0;
	}
    }
}

static double int_loglik (const double *theta, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    double x0, x1, z0, z1;
    double derivs = 0.0, derivb = 0.0;
    double sigma, ndxt, ll = 0.0;
    int i, t, k = IC->k;

    sigma = exp(theta[k-1]);

    for (i=0; i<k; i++) {
	IC->g[i] = 0;
    }

    loglik_prelim(theta, IC);

    for (t=0; t<IC->nobs; t++) {
	x0 = IC->lo[t];
	x1 = IC->hi[t];
	ndxt = IC->ndx[t];

	switch (IC->obstype[t]) {
	case INT_LOW:
	    z1 = (x1 - ndxt)/sigma;
	    derivb = -IC->f1[t]/sigma;
	    derivs = -IC->f1[t]*z1;
	    break;
	case INT_HIGH:
	    z0 = (x0 - ndxt)/sigma;
	    derivb = IC->f0[t]/sigma;
	    derivs = IC->f0[t]*z0;
	    break;
	case INT_MID:
	    z0 = (x0 - ndxt)/sigma;
	    z1 = (x1 - ndxt)/sigma;
	    derivb = (IC->f0[t] - IC->f1[t])/sigma;
	    derivs = (IC->f0[t]*z0 - IC->f1[t]*z1);
	    break;
	case INT_POINT:
	    z0 = (x0 - ndxt)/sigma;
	    derivb = z0/sigma;
	    derivs = z0*z0 - 1;
	}

	ll += log(IC->dP[t]);

	for (i=0; i<IC->nx; i++) {
	    IC->G[i][t] = derivb * IC->X[i][t];
	    IC->g[i] += IC->G[i][t]; 
	}

	IC->G[k-1][t] = derivs;
	IC->g[k-1] += derivs;
    }
    
#if INTDEBUG > 1
    fprintf(stderr, "ll = %16.10f\n", ll);

    for (i=0; i<IC->k; i++) {
	fprintf(stderr, "g[%d] = %16.10f\t", i, IC->g[i]);
    }
    fputc('\n', stderr);
#endif

    return ll;
}

static int int_score (double *theta, double *s, int npar, BFGS_CRIT_FUNC ll, 
		      void *ptr)
{
    int_container *IC = (int_container *) ptr;
    int i;

    for (i=0; i<npar; i++) {
	s[i] = IC->g[i];
    }

    return 0;
}

int int_ahess (double *theta, gretl_matrix *V, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    double z0 = 0, z1 = 0, mu = 0;
    double q0 = 0, q1 = 0, nu = 0;
    double x, vij, x0, x1;    
    double f0, f1, Hss = 0;
    double sigma, ndxt, lambda = 0;
    int i, j, t, k = IC->k;
    int err = 0;

    sigma = exp(theta[k-1]);
    loglik_prelim(theta, IC);

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    gretl_matrix_set(V, i, j, 0);
	}
    }

#if INTDEBUG > 1
    fprintf(stderr, "int_ahess:\n");
    for (i=0; i<k; i++) {
	fprintf(stderr, "theta[%2d] = %20.16f\n", i, theta[i]);
    }
#endif

    for (t=0; t<IC->nobs; t++) {
	/* ndx, dP, f0, f1 etc should already be there */
	x0 = IC->lo[t];
	x1 = IC->hi[t];
	ndxt = IC->ndx[t];
	f0 = IC->f0[t];
	f1 = IC->f1[t];

	switch (IC->obstype[t]) {
	case INT_LOW:
	    z1 = (x1 - ndxt)/sigma;
	    mu     = -f1/sigma;
	    lambda = mu * z1;
	    q1     = (z1*z1 - 1.0);
	    nu     = mu * q1;
	    break;
	case INT_HIGH:
	    z0 = (x0 - ndxt)/sigma;
	    mu     = f0/sigma;
	    lambda = mu * z0;
	    q0     = (z0*z0 - 1.0);
	    nu     = mu * q0;
	    break;
	case INT_MID:
	    z0 = (x0 - ndxt)/sigma;
	    z1 = (x1 - ndxt)/sigma;
	    mu     = (f0 - f1)/sigma;
	    lambda = (f0*z0 - f1*z1)/sigma;
	    q0     = (z0*z0 - 1.0);
	    q1     = (z1*z1 - 1.0);
	    nu     = (f0 * q0 - f1 * q1)/sigma;
	    break;
	case INT_POINT:
	    z0 = (x0 - ndxt)/sigma;
	}

	if (IC->obstype[t] == INT_POINT) {
	    x = 1.0 / (sigma*sigma);
	} else {
	    x = (mu*mu - lambda/sigma);
	}

	for (i=0; i<IC->nx; i++) {
	    x0 = IC->X[i][t];
	    for (j=i; j<IC->nx; j++) {
		x1 = IC->X[j][t];
		vij = gretl_matrix_get(V, i, j);
		vij += x * x0 * x1;
		gretl_matrix_set(V, i, j, vij);
	    }
	}

	if (IC->obstype[t] == INT_POINT) {
	    x = 2 * z0 / sigma;
	} else {
	    x = (mu*lambda*sigma - nu);
	}

	for (i=0; i<IC->nx; i++) {
	    x0 = IC->X[i][t];
	    vij = gretl_matrix_get(V, i, k-1);
	    vij += x * x0;
	    gretl_matrix_set(V, i, k-1, vij);
	}

	if (IC->obstype[t] == INT_POINT) {
	    x = 2 * (z0*z0);
	    Hss += x;
	} else {
	    x = lambda*sigma;
	    Hss += x*(x+1) - (f0 * q0 * z0 - f1 * q1 *z1);
	}
    }

    gretl_matrix_set(V, k-1, k-1, Hss);

    for (i=0; i<k; i++) {
	for (j=i; j<k; j++) {
	    vij = gretl_matrix_get(V, i, j);
	    gretl_matrix_set(V, j, i, vij);
	}
    }

#if INTDEBUG > 1
    gretl_matrix_print(V, "Hessian (analytical)");
#endif

    return err;
}

static int 
int_ahess_inverse (double *theta, gretl_matrix *V, void *ptr)
{
    int err = int_ahess(theta, V, ptr);

    if (!err) {
	err = gretl_invert_symmetric_matrix(V);
    }

    return err;
}

static gretl_matrix *intreg_sandwich (int_container *IC,
				      gretl_matrix *H,
				      int *err)
{
    gretl_matrix *G, *S, *V;
    double x;
    int k = IC->k;
    int i, j, jj = 0;

    G = gretl_matrix_alloc(IC->nobs, k);
    S = gretl_matrix_alloc(k, k);
    V = gretl_matrix_alloc(k, k);

    if (G == NULL || S == NULL || V == NULL) {
	gretl_matrix_free(V);
	V = NULL;
	*err = E_ALLOC;
	goto bailout;
    }

    /* form the G matrix */
    for (j=0; j<k; j++) {
	for (i=0; i<IC->nobs; i++) {
	    x = IC->G[j][i];
	    G->val[jj++] = x;
	}
    }

    /* form S = G'G */
    *err = gretl_matrix_multiply_mod(G, GRETL_MOD_TRANSPOSE,
				     G, GRETL_MOD_NONE,
				     S, GRETL_MOD_NONE);

    if (!*err) {
	/* form sandwich: V = H^{-1} S H^{-1} */
	*err = gretl_matrix_qform(H, GRETL_MOD_NONE, S,
				  V, GRETL_MOD_NONE);
    }

 bailout:
    
    gretl_matrix_free(G);
    gretl_matrix_free(S);
    gretl_matrix_free(H);

    return V;
}

static gretl_matrix *intreg_VCV (int_container *IC, gretlopt opt,
				 int *err)
{
    gretl_matrix *H = gretl_zero_matrix_new(IC->k, IC->k);

    if (H == NULL) {
	*err = E_ALLOC;
    } else {
	*err = int_ahess_inverse(IC->theta, H, IC);
    }

    if (!*err && (opt & OPT_R)) {
	return intreg_sandwich(IC, H, err);
    } else {
	return H;
    }
}

static void int_compute_gresids (int_container *IC)
{
    int t, k = IC->k;
    double sigma = exp(IC->theta[k-1]);

    for (t=0; t<IC->nobs; t++) {
	switch (IC->obstype[t]) {
	case INT_LOW:
	    IC->uhat[t] = -sigma * IC->f1[t];
	    break;
	case INT_HIGH:
	    IC->uhat[t] = sigma * IC->f0[t];
	    break;
	case INT_MID:
	    IC->uhat[t] = sigma * (IC->f0[t] - IC->f1[t]);
	    break;
	case INT_POINT:
	    IC->uhat[t] = IC->lo[t] - IC->ndx[t];
	}
    }
} 

/*
  Wald test for zeroing all coefficient apart from the constant,
  which is assumed to be the first explanatory variable
*/

static double chisq_overall_test (int_container *IC)
{
    gretl_vector *b;
    gretl_matrix *V;
    int k = IC->nx - 1;
    int i, j, n, err;
    double ret = NADBL;

    b = gretl_vector_alloc(k);
    V = gretl_matrix_alloc(k, k);

    if (b == NULL || V == NULL) {
	gretl_vector_free(b);
	gretl_matrix_free(V);
	return NADBL;
    }

    n = IC->nx;
    for (i=0; i<k; i++) {
	gretl_vector_set(b, i, IC->pmod->coeff[i+1]);
	for (j=i; j<k; j++) {
	    gretl_matrix_set(V, i, j, IC->pmod->vcv[n]);
	    gretl_matrix_set(V, j, i, IC->pmod->vcv[n]);
	    n++;
	}
    }

    err = gretl_invert_symmetric_matrix(V);

    if (!err) {
	ret = gretl_scalar_qform(b, V, &err);
    }

    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return ret;
} 

/* 
   Here we exploit the following properties of sub-support normal 
   moments: if x ~ N(0,1) and A = [a,b], then

   E(x^n | x \in A) = (n-1) E(x^{n-2} | x \in A) + 
   + \frac{a^{n-1}\phi(a) - b^{n-1}\phi(b)}{\Phi(b)-\Phi(a)}

   We have:

   E(x^0 | x \in A) = 1 (trivially) 
   E(x^1 | x \in A) = f0 - f1
   E(x^2 | x \in A) = 1 + lo*f0 - hi*f1
   ...

   Returns the matrix of the orthogonality conditions (conditional
   moments); the non-zero ones are stored in the two pointers sm3
   and sm4.
*/

static gretl_matrix *cond_moments (int_container *IC, double *sm3, 
				   double *sm4)
{
    gretl_matrix *ret = NULL;
    int n = IC->nobs;  
    int k = IC->k;  
    int noc = k + 2;  
    double a, b, phi0, phi1, m1, m2, u, x, y;
    double m3 = 0.0, m4 = 0.0;
    double sigma = exp(IC->theta[k - 1]);
    int i, t;

    *sm3 = 0.0;
    *sm4 = 0.0;

    ret = gretl_matrix_alloc(n, noc);
    if (ret == NULL) {
	return NULL;
    }

    for (t=0; t<n; t++) {
	switch (IC->obstype[t]) {
	case INT_LOW:
	    b = (IC->hi[t] - IC->ndx[t])/sigma;
	    phi1 = y = IC->f1[t];
	    m1 = -phi1;
	    m2 = 1 - (y *= b);
	    m3 = 2*m1 - (y *= b);
	    m4 = 3*m2 - (y *= b);
	    break;
	case INT_HIGH:
	    a = (IC->lo[t] - IC->ndx[t])/sigma;
	    m1 = phi0 = x = IC->f0[t];
	    m2 = 1 + (x *= a);
	    m3 = 2*m1 + (x *= a);
	    m4 = 3*m2 + (x *= a);
	    break;
	case INT_MID:
	    a = (IC->lo[t] - IC->ndx[t])/sigma;
	    b = (IC->hi[t] - IC->ndx[t])/sigma;
	    phi0 = x = IC->f0[t];
	    phi1 = y = IC->f1[t];
	    m1 = phi0 - phi1;
	    m2 = 1 + (x *= a) - (y *= b);
	    m3 = 2*m1 + (x *= a) - (y *= b);
	    m4 = 3*m2 + (x *= a) - (y *= b);
	    break;
	case INT_POINT:
	    u = IC->uhat[t] / sigma;
	    m3 = u*u*u;
	    m4 = m3*u;
	}

	m4 -= 3;
	*sm3 += m3;
	*sm4 += m4;

	for (i=0; i<k; i++) {
	    gretl_matrix_set(ret, t, i, IC->G[i][t]);
	}

	gretl_matrix_set(ret, t, k, m3);
	gretl_matrix_set(ret, t, k+1, m4);
    }

    return ret;
}

/*
  Conditional moment test for normality: in practice, a regression
  where the dep. var is 1 and the matrix of explanatory variables
  is [ G | sk | ku ], where:

  * G is the by-obs score matrix (should already be in IC by now)
  * sk[t] is an unbiased estimator of \epsilon_t^3
  * ku[t] is an unbiased estimator of (\epsilon_t^4 - 3)

  The test statistic is the explained sum of squares.
*/

static int intreg_normtest (int_container *IC, double *teststat)
{
    gretl_matrix *condmom;
    gretl_matrix *GG = NULL;
    gretl_matrix *g = NULL;
    double skew, kurt;
    int noc = IC->k + 2;
    int err = 0;

    condmom = cond_moments(IC, &skew, &kurt);
    if (condmom == NULL) {
	return E_ALLOC;
    }

    GG = gretl_matrix_XTX_new(condmom);
    g = gretl_zero_matrix_new(1, noc);

    if (GG == NULL || g == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_invert_symmetric_matrix(GG);
	if (!err) {
	    gretl_vector_set(g, noc-2, skew);
	    gretl_vector_set(g, noc-1, kurt);
	    *teststat = gretl_scalar_qform(g, GG, &err);
	} else {
	    gretl_matrix_print(condmom, "conditional moment matrix");
	    *teststat = NADBL;
	}
    }

    gretl_matrix_free(condmom);
    gretl_matrix_free(GG);
    gretl_matrix_free(g);

    return err;
}

static int fill_intreg_model (int_container *IC, gretl_matrix *V,
			      int fncount, int grcount,
			      const DATASET *dset,
			      gretlopt opt)
{
    MODEL *pmod = IC->pmod;
    double x, ndx, u;
    int i, j, k = IC->k;
    int obstype, vtype;
    int err = 0;

    pmod->ci = (opt & OPT_T)? TOBIT : INTREG;
    pmod->lnL = IC->ll;
    mle_criteria(pmod, 1);
    pmod->sigma = exp(IC->theta[k-1]);

    for (i=0; i<IC->nx; i++) {
	pmod->coeff[i] = IC->theta[i];
    }

    err = gretl_model_write_vcv(pmod, V, IC->nx);
    if (err) {
	return err;
    }

    vtype = (opt & OPT_R)? VCV_QML : VCV_HESSIAN;
    gretl_model_set_vcv_info(pmod, VCV_ML, vtype);

    /* get the s.e. of sigma via the delta method */
    x = gretl_matrix_get(V, k-1, k-1);
    x *= pmod->sigma * pmod->sigma;
    gretl_model_set_double(pmod, "se_sigma", sqrt(x));

    int_compute_gresids(IC);

    j = 0;

    for (i=IC->t1; i<=IC->t2; i++) {
	if (!na(pmod->uhat[i])) {
	    obstype = IC->obstype[j];
	    IC->typecount[obstype] += 1;
	    ndx = IC->ndx[j];
	    u = IC->uhat[j];
	    pmod->uhat[i] = u;
	    pmod->yhat[i] = ndx;
	    j++;
	}
    }

    /* mask invalid statistics */
    pmod->ybar = pmod->sdy = NADBL;
    pmod->ess = pmod->tss = NADBL;
    pmod->fstt = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;

    if (grcount > 0) {
	/* BFGS */
	gretl_model_set_int(pmod, "fncount", fncount);
	gretl_model_set_int(pmod, "grcount", grcount);
    } else {
	/* Newton */
	gretl_model_set_int(pmod, "iters", fncount);
    }

    gretl_model_set_int(pmod, "n_left", IC->typecount[0]);
    gretl_model_set_int(pmod, "n_right", IC->typecount[2]);

    if (!(opt & OPT_T)) {
	gretl_model_set_int(pmod, "lovar", IC->lov);
	gretl_model_set_int(pmod, "hivar", IC->hiv);
	gretl_model_set_int(pmod, "n_both", IC->typecount[1]);
	gretl_model_set_int(pmod, "n_point", IC->typecount[3]);
    }

    if (opt & OPT_R) {
	pmod->opt |= OPT_R;
    }

    if (IC->pmod->ifc && IC->nx > 1) {
	pmod->chisq = chisq_overall_test(IC);
    } else {
	pmod->chisq = NADBL;
    }

    if (!(opt & OPT_T)) {
	/* the variable at list[1] will be deleted */
	pmod->list[1] = 0;
    }

    return 0;
}

static int do_interval (int *list, DATASET *dset, MODEL *mod, 
			gretlopt opt, PRN *prn) 
{
    int_container *IC;
    gretl_matrix *V = NULL;
    int maxit, fncount, grcount = 0;
    double toler, normtest = NADBL;
    gretlopt maxopt = opt & OPT_V;
    int use_bfgs = 0;
    int err = 0;

    IC = int_container_new(list, dset, mod);
    if (IC == NULL) {
	return E_ALLOC;
    }

    BFGS_defaults(&maxit, &toler, INTREG);

    if (libset_get_int(GRETL_OPTIM) == OPTIM_BFGS) {
	use_bfgs = 1;
    }    

    if (use_bfgs) {
	err = BFGS_max(IC->theta, IC->k, maxit, toler, 
		       &fncount, &grcount, int_loglik, C_LOGLIK,
		       int_score, IC, NULL, maxopt, prn);
    } else {
	double crittol = 1.0e-07;
	double gradtol = 1.0e-07;
	
	err = newton_raphson_max(IC->theta, IC->k, maxit, crittol, gradtol,
				 &fncount, C_LOGLIK, int_loglik, 
				 int_score, int_ahess, IC, maxopt, prn);
    }

    if (!err) {
	IC->ll = int_loglik(IC->theta, IC);
    }

    if (!err) {
	V = intreg_VCV(IC, opt, &err);
    }

    if (!err) {
	err = fill_intreg_model(IC, V, fncount, grcount, dset, opt);
    }

    if (!err) {
	err = intreg_normtest(IC, &normtest);
	if (!err) {
	    err = gretl_model_add_normality_test(IC->pmod, normtest);
	}
	/* don't let this be a show-stopper */
	err = 0;
    }

    gretl_matrix_free(V);
    int_container_destroy(IC);

    return err;
}

/* if the list contains a constant, ensure that it appears
   as the first regressor, in position 3 */

static void maybe_reposition_const (int *list, const DATASET *dset)
{
    int cpos = gretl_list_const_pos(list, 4, dset);
    int i;

    if (cpos > 0) {
	for (i=cpos; i>3; i--) {
	    list[i] = list[i-1];
	}
	list[3] = 0;
    }	
}

MODEL interval_estimate (int *list, DATASET *dset,
			 gretlopt opt, PRN *prn) 
{
    MODEL model;
    int *initlist = NULL;

    gretl_model_init(&model);
    
    if (list[0] > 3) {
	maybe_reposition_const(list, dset);
    }

    /* create extra series for model initialization and
       corresponding regression list
    */
    model.errcode = create_midpoint_y(list, dset, &initlist);
    if (model.errcode) {
	return model;
    }

    /* run initial OLS */
    model = lsq(initlist, dset, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "interval_estimate: initial OLS failed\n");
	free(initlist);
	return model;
    }

#if INTDEBUG
    pprintf(prn, "interval_estimate: initial OLS\n");
    printmodel(&model, dset, OPT_S, prn);
#endif

    /* clean up midpoint-y */
    dataset_drop_last_variables(1, dset);
    free(initlist);

    /* do the actual analysis */
    model.errcode = do_interval(list, dset, &model, opt, prn);

    clear_model_xpx(&model);

    return model;
}

static int tobit_add_lo_hi (MODEL *pmod, double llim, double rlim,
			    DATASET *dset, int **plist)
{
    double *y, *lo, *hi;
    int *list;
    int lv, hv;
    int i, t, err;

    err = dataset_add_series(2, dset);
    if (err) {
	return err;
    }

    lv = dset->v - 2;
    hv = dset->v - 1;

    lo = dset->Z[lv];
    hi = dset->Z[hv];
    y = dset->Z[pmod->list[1]];

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	if (na(y[t])) {
	    lo[t] = hi[t] = NADBL;
	} else if (y[t] > llim && y[t] < rlim) {
	    lo[t] = hi[t] = y[t];
	} else if (y[t] >= rlim) {
	    lo[t] = rlim;
	    hi[t] = NADBL;
	} else if (y[t] <= llim) {
	    lo[t] = NADBL;
	    hi[t] = llim;
	} 
    }

    list = gretl_list_new(pmod->list[0] + 1);

    if (list == NULL) {
	err = E_ALLOC;
    } else {
	list[1] = lv;
	list[2] = hv;
	for (i=3; i<=list[0]; i++) {
	    list[i] = pmod->list[i-1];
	}
	*plist = list;
    }

    return err;
}

MODEL tobit_via_intreg (int *list, double llim, double rlim,
			DATASET *dset, gretlopt opt, 
			PRN *prn) 
{
    MODEL model;
    int *ilist = NULL;
    int origv = dset->v;

    /* run initial OLS */
    model = lsq(list, dset, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "intreg: initial OLS failed\n");
	return model;
    }

#if INTDEBUG
    pprintf(prn, "tobit_via_intreg: initial OLS\n");
    printmodel(&model, dset, OPT_S, prn);
#endif

    model.errcode = tobit_add_lo_hi(&model, llim, rlim,
				    dset, &ilist);

    if (!model.errcode) {
	/* do the actual analysis */
	model.errcode = do_interval(ilist, dset, &model, 
				    opt | OPT_T, prn);
    }

    clear_model_xpx(&model);

    if (!model.errcode) {
	if (opt & OPT_L) {
	    model.opt |= OPT_L;
	    gretl_model_set_double(&model, "llimit", llim);
	}
	if ((opt & OPT_M) && !na(rlim)) {
	    model.opt |= OPT_M;
	    gretl_model_set_double(&model, "rlimit", rlim);
	}	
    }

    /* clean up extra data */
    dataset_drop_last_variables(dset->v - origv, dset);
    free(ilist);

    return model;
}
