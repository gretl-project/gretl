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

static int_container *int_container_new (int *list, double **Z, 
					 DATAINFO *pdinfo, MODEL *mod)
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

	    x0 = Z[vlo][t];
	    x1 = Z[vhi][t];
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
		IC->X[i][s] = Z[list[i+3]][t];
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

#if INTDEBUG > 1
    fprintf(stderr, "nobs = %d\n", IC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", IC->t1, IC->t2);
    fprintf(stderr, "k = %d\n", IC->k);
    fprintf(stderr, "nx = %d\n", IC->nx);
#endif

    return IC;
}

static int create_midpoint_y (int *list, double ***pZ, DATAINFO *pdinfo, 
			      int **initlist)
{
    int n = list[0];
    int mpy = pdinfo->v;
    int lv, hv;
    double x0, x1;
    int i, t, err;

    err = dataset_add_series(1, pZ, pdinfo);
    if (err) {
	return err;
    }

    lv = list[1];
    hv = list[2];

    for (t=pdinfo->t1; t<=pdinfo->t2 && !err; t++) {
	x0 = (*pZ)[lv][t];
	x1 = (*pZ)[hv][t];

	if (na(x0)) {
	    (*pZ)[mpy][t] = x1;
	} else if (na(x1)) {
	    (*pZ)[mpy][t] = x0;
	} else if (x0 > x1) {
	    gretl_errmsg_sprintf(_("Obs %d: lower bound (%g) "
				   "exceeds upper (%g)"), t + 1,
				 x0, x1);
	    err = E_DATA;
	} else {
	    (*pZ)[mpy][t] = 0.5 * (x0 + x1);
	}
    }

    if (err) {
	return err;
    }

    *initlist = gretl_list_new(n - 1);

    if (*initlist == NULL) {
	err = E_ALLOC;
    } else {
	(*initlist)[1] = mpy;
	for (i=3; i<=n; i++) {
	    (*initlist)[i-1] = list[i];
	}

#if INTDEBUG > 1
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    fprintf(stderr, "%2d: %g\n", t, (*pZ)[mpy][t]);
	}
#endif
    }

    return err;
}

static double int_loglik (const double *theta, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    double x0, x1, z0, z1;
    double ndx, ll = 0.0;
    double derivs = 0.0, derivb = 0.0;
    int i, t;
    int k = IC->k;
    double sigma = exp(theta[k-1]);

    for (i=0; i<k; i++) {
	IC->g[i] = 0;
    }

    for (t=0; t<IC->nobs; t++) {
	x0 = IC->lo[t];
	x1 = IC->hi[t];

	ndx = 0;
	for (i=0; i<IC->nx; i++) {
	    ndx += theta[i] * IC->X[i][t];
	}

	IC->ndx[t] = ndx;

	switch (IC->obstype[t]) {
	case INT_LOW:
	    z1 = (x1 - ndx)/sigma;
	    IC->dP[t] = normal_cdf(z1);
	    IC->f1[t] = normal_pdf(z1) / IC->dP[t];
	    derivb = -IC->f1[t]/sigma;
	    derivs = -IC->f1[t]*z1;
	    break;
	case INT_HIGH:
	    z0 = (x0 - ndx)/sigma;
	    IC->dP[t] = normal_cdf_comp(z0);
	    IC->f0[t] = normal_pdf(z0) / IC->dP[t];
	    derivb = IC->f0[t]/sigma;
	    derivs = IC->f0[t]*z0;
	    break;
	case INT_MID:
	    z0 = (x0 - ndx)/sigma;
	    z1 = (x1 - ndx)/sigma;
	    IC->dP[t] = normal_cdf(z1) - normal_cdf(z0);
	    IC->f0[t] = normal_pdf(z0) / IC->dP[t];
	    IC->f1[t] = normal_pdf(z1) / IC->dP[t];
	    derivb = (IC->f0[t] - IC->f1[t])/sigma;
	    derivs = (IC->f0[t]*z0 - IC->f1[t]*z1);
	    break;
	case INT_POINT:
	    z0 = (x0 - ndx)/sigma;
	    IC->dP[t] = normal_pdf(z0)/sigma;
	    derivb = z0/sigma;
	    derivs = z0*z0 - 1;
	}

	ll += log(IC->dP[t]);

	for (i=0; i<IC->nx; i++) {
	    IC->G[i][t] = derivb * IC->X[i][t];
	    IC->g[i] += IC->G[i][t];
	}

	IC->G[k-1][t] = derivs;
	IC->g[k-1] += IC->G[k-1][t];

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

    return 1;
}

static gretl_matrix *int_hess (int_container *IC, int *err)
{
    double x, smal = 1.0e-07;
    gretl_matrix *V = NULL;
    double *q, *g, *gplus, *gminus, *hss;
    int k = IC->k;
    int i, j, n;

    q = malloc(k * sizeof *q); 
    g = malloc(k * sizeof *g); 
    gplus = malloc(k * sizeof *gplus); 
    gminus = malloc(k * sizeof *gminus);
    hss = malloc(k * k * sizeof *hss);

    if (q == NULL || g == NULL || gplus == NULL ||
	gminus == NULL || hss == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }	

    for (i=0; i<k; i++) {
	g[i] = IC->g[i];
	q[i] = IC->theta[i];
    }

    n = 0;
    for (i=0; i<k; i++) {
	q[i] += smal;
	int_loglik(q, IC);
	for (j=0; j<k; j++) {
	    gplus[j] = IC->g[j];
	}

	q[i] -= 2 * smal;
	int_loglik(q, IC);
	for (j=0; j<k; j++) {
	    gminus[j] = IC->g[j];
	}

	for (j=0; j<k; j++) {
	    hss[n++] = (gplus[j] - gminus[j]) / (2 * smal);
	}
    }

    for (i=0; i<k; i++) {
	for (j=i; j<k; j++) {
	    x = -0.5 * (hss[i+k*j] + hss[i*k+j]);
	    gretl_matrix_set(V, i, j, x);
	    gretl_matrix_set(V, j, i, x);
	}
    }

    *err = gretl_invert_symmetric_matrix(V);
#if INTDEBUG
    gretl_matrix_print(V, "V = -H^{-1}");
#endif

 bailout:

    free(q);
    free(g);
    free(gplus);
    free(gminus);
    free(hss);

    return V;
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
    gretl_matrix *H = int_hess(IC, err);

    if (!*err && (opt & OPT_R)) {
	return intreg_sandwich(IC, H, err);
    } else {
	return H;
    }
}

static void int_compute_gresids (int_container *IC)
{
    int t, k = IC->k;
    double ndx, sigma = exp(IC->theta[k-1]);

    for (t=0; t<IC->nobs; t++) {
	ndx = IC->ndx[t];

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
	    IC->uhat[t] = IC->lo[t] - ndx;
	}
    }
} 

static double chisq_overall_test (int_container *IC)
{
    /*
      Wald test for zeroing all coefficient apart from the constant,
      which is assumed to be the first explanatory variable
    */

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

static gretl_matrix *cond_moments (int_container *IC, double *sm3, 
				   double *sm4)
{
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

    int n = IC->nobs;  
    int k = IC->k;  
    int noc = k + 2;  
    double *ndx = IC->ndx;
    double *lo = IC->lo; 
    double *hi = IC->hi;
    double *f0 = IC->f0;
    double *f1 = IC->f1;
    gretl_matrix *ret = NULL;
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
	    b = (hi[t] - ndx[t])/sigma;
	    phi1 = y = f1[t];
	    m1 = -phi1;
	    m2 = 1 - (y *= b);
	    m3 = 2*m1 - (y *= b);
	    m4 = 3*m2 - (y *= b);
	    break;
	case INT_HIGH:
	    a = (lo[t] - ndx[t])/sigma;
	    m1 = phi0 = x = f0[t];
	    m2 = 1 + (x *= a);
	    m3 = 2*m1 + (x *= a);
	    m4 = 3*m2 + (x *= a);
	    break;
	case INT_MID:
	    a = (lo[t] - ndx[t])/sigma;
	    b = (hi[t] - ndx[t])/sigma;

	    phi0 = x = f0[t];
	    phi1 = y = f1[t];
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

	for(i=0; i<k; i++) {
	    gretl_matrix_set(ret, t, i, IC->G[i][t]);
	}

	gretl_matrix_set(ret, t, k, m3);
	gretl_matrix_set(ret, t, k+1, m4);
    }

    return ret;
}

static int intreg_normtest (int_container *IC, double *teststat)
{
    /*
      Conditonal moment test for normality: in practice, a regression
      where the dep. var is 1 and the matrix of explanatory variables
      is [ G | sk | ku ], where:

      * G is the by-obs score matrix (should already be in IC by now)
      * sk[t] is an unbiased estimator of \epsilon_t^3
      * ku[t] is an unbiased estimator of (\epsilon_t^4 - 3)

      The test statistic is the explained sum of squares.
    */

    gretl_matrix *condmom;
    double skew, kurt;
    gretl_matrix *GG = NULL;
    gretl_matrix *g = NULL;
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
	goto bailout;
    }

    err = gretl_invert_symmetric_matrix(GG);
    if (!err) {
	gretl_vector_set(g, noc-2, skew);
	gretl_vector_set(g, noc-1, kurt);
	*teststat = gretl_scalar_qform(g, GG, &err);
    } else {
	gretl_matrix_print(condmom, "conditional moment matrix");
	*teststat = NADBL;
    }

 bailout:

    gretl_matrix_free(condmom);
    gretl_matrix_free(GG);
    gretl_matrix_free(g);

    return err;
}

static int add_norm_test_to_model (MODEL *pmod, double X2)
{
    ModelTest *test = model_test_new(GRETL_TEST_NORMAL);
    int err = 0;

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_NORMAL_CHISQ);
        model_test_set_dfn(test, 2);
        model_test_set_value(test, X2);
        model_test_set_pvalue(test, chisq_cdf_comp(2, X2));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = 1;
    }

    return err;
}

static int fill_intreg_model (int_container *IC, gretl_matrix *V,
			      gretlopt opt, const DATAINFO *pdinfo)
{
    MODEL *pmod = IC->pmod;
    double x, ndx, u;
    int i, j, n, m, k = IC->k, nx = IC->nx;
    int obstype;

    pmod->ci = INTREG;

    pmod->lnL = IC->ll;
    mle_criteria(pmod, 1);

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
    }

    m =  k * (k - 1) / 2;

    pmod->vcv = malloc(m * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return E_ALLOC;
    }

    n = 0;
    for (i=0; i<nx; i++) {
	pmod->coeff[i] = IC->theta[i];
	for (j=i; j<nx; j++) {
	    x = gretl_matrix_get(V, i, j);
	    if (i == j) {
		pmod->sderr[i] = sqrt(x);
	    }
	    pmod->vcv[n++] = x;
	}
    }

    pmod->sigma = exp(IC->theta[k-1]);

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

    pmod->ybar = pmod->sdy = NADBL;
    pmod->ess = pmod->tss = NADBL;
    pmod->fstt = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;

    gretl_model_set_int(pmod, "lovar", IC->lov);
    gretl_model_set_int(pmod, "hivar", IC->hiv);
    gretl_model_set_int(pmod, "n_left", IC->typecount[0]);
    gretl_model_set_int(pmod, "n_both", IC->typecount[1]);
    gretl_model_set_int(pmod, "n_right", IC->typecount[2]);
    gretl_model_set_int(pmod, "n_point", IC->typecount[3]);

    if (opt & OPT_R) {
	gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
	pmod->opt |= OPT_R;
    }

    if (IC->pmod->ifc && IC->nx > 1) {
	pmod->chisq = chisq_overall_test(IC);
    } else {
	pmod->chisq = NADBL;
    }

    /* the variable at list[1] will be deleted */
    pmod->list[1] = 0;

    return 0;
}

static int do_interval (int *list, double **Z, DATAINFO *pdinfo, 
			MODEL *mod, gretlopt opt, PRN *prn) 
{
    int_container *IC;
    gretl_matrix *V = NULL;
    int maxit, fncount, grcount;
    double toler, normtest = NADBL;
    int err = 0;

    IC = int_container_new(list, Z, pdinfo, mod);
    if (IC == NULL) {
	return E_ALLOC;
    }

    BFGS_defaults(&maxit, &toler, INTREG);

    err = BFGS_max(IC->theta, IC->k, maxit, toler, 
		   &fncount, &grcount, int_loglik, C_LOGLIK,
		   int_score, IC, opt & OPT_V, prn);

    if (!err) {
	IC->ll = int_loglik(IC->theta, IC);
    }

    if (!err) {
	V = intreg_VCV(IC, opt, &err);
    }

    if (!err) {
	err = fill_intreg_model(IC, V, opt, pdinfo);
    }

    if (!err) {
	err = intreg_normtest(IC, &normtest);
	if (!err) {
	    err = add_norm_test_to_model(IC->pmod, normtest);
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

static void maybe_reposition_const (int *list, const double **Z,
				    const DATAINFO *pdinfo)
{
    int cpos = gretl_list_const_pos(list, 4, Z, pdinfo);
    int i;

    if (cpos > 0) {
	for (i=cpos; i>3; i--) {
	    list[i] = list[i-1];
	}
	list[3] = 0;
    }	
}

MODEL interval_estimate (int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn) 
{
    MODEL model;
    int *initlist = NULL;

    gretl_model_init(&model);
    
    if (list[0] > 3) {
	maybe_reposition_const(list, (const double **) *pZ, pdinfo);
    }

    /* create extra variable for model initialization and
       corresponding list
    */
    model.errcode = create_midpoint_y(list, pZ, pdinfo, &initlist);
    if (model.errcode) {
	return model;
    }

    /* run initial OLS */
    model = lsq(initlist, pZ, pdinfo, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "interval_estimate: initial OLS failed\n");
	free(initlist);
	return model;
    }

#if INTDEBUG
    pprintf(prn, "interval_estimate: initial OLS\n");
    printmodel(&model, pdinfo, OPT_S, prn);
#endif

    /* clean up midpoint-y */
    dataset_drop_last_variables(1, pZ, pdinfo);
    free(initlist);

    /* do the actual analysis */
    model.errcode = do_interval(list, *pZ, pdinfo, &model, opt, prn);

    clear_model_xpx(&model);

    return model;
}
