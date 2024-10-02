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
#include "version.h"
#include "libset.h"
#include "gretl_bfgs.h"
#include "gretl_normal.h"
#include "matrix_extra.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#define INTDEBUG 0
#define BIGDEBUG 0

enum {
    INT_LOW,   /* no lower bound */
    INT_MID,   /* both bounds */
    INT_HIGH,  /* no upper bound */
    INT_POINT, /* point value */
    INT_FPOINT /* forced to point value */
};

typedef struct int_container_ int_container;

struct int_container_ {
    MODEL *pmod;      /* model struct, initially containing OLS */
    int hiv, lov;     /* variable numbers for limit series */
    double ll;        /* loglikelihood */
    double *dspace;   /* workspace */
    double *hi, *lo;  /* limit series for dependent variable */
    int *obstype;     /* dependent variable classifier */
    int typecount[5]; /* no. of obs by category */
    gretl_matrix *X;  /* regressors */
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
    gretl_matrix *G;  /* score matrix by observation */
    double *g;        /* total score vector */
};

static void int_container_destroy (int_container *IC)
{
    free(IC->dspace);

    free(IC->theta);
    free(IC->g);

    free(IC->obstype);
    free(IC->list);

    gretl_matrix_free(IC->X);
    gretl_matrix_free(IC->G);

    free(IC);
}

static int_container *int_container_new (int *list, DATASET *dset,
					 MODEL *mod)
{
    int_container *IC;
    int nobs, i, s, t, vlo = list[1], vhi = list[2];
    double x0, x1, xti;

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
    for (i=0; i<5; i++) {
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

    /* big matrices */
    IC->X = gretl_matrix_alloc(nobs, IC->nx);
    IC->G = gretl_matrix_alloc(nobs, IC->k);

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
		xti = dset->Z[list[i+3]][t];
		gretl_matrix_set(IC->X, s, i, xti);
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

    err = dataset_add_series(dset, 1);
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
    double ndxt, z0, z1, x0, x1, xti;
    double sigma = exp(theta[k-1]);

#if defined(_OPENMP)
    int par_t = IC->nobs >= 2000;
#pragma omp parallel for private(t, i, ndxt, xti, x0, x1, z0, z1) if (par_t)
#endif
    for (t=0; t<IC->nobs; t++) {
	ndxt = 0.0;
	for (i=0; i<IC->nx; i++) {
	    xti = gretl_matrix_get(IC->X, t, i);
	    ndxt += theta[i] * xti;
	}

	IC->ndx[t] = ndxt;
	x0 = IC->lo[t];
	x1 = IC->hi[t];

	/* reset forced observations to mid just in case */
	if (IC->obstype[t] == INT_FPOINT) {
	    IC->obstype[t] = INT_MID;
	}

	switch (IC->obstype[t]) {
	case INT_LOW:
	    z1 = (x1 - ndxt)/sigma;
	    IC->dP[t] = normal_cdf(z1);
	    IC->f0[t] = 0;
	    IC->f1[t] = invmills(-z1);
	    break;
	case INT_HIGH:
	    z0 = (x0 - ndxt)/sigma;
	    IC->dP[t] = normal_cdf_comp(z0);
	    IC->f0[t] = invmills(z0);
	    IC->f1[t] = 0;
	    break;
	case INT_MID:
	    z0 = (x0 - ndxt)/sigma;
	    z1 = (x1 - ndxt)/sigma;
	    IC->dP[t] = normal_cdf(z1) - normal_cdf(z0);
	    if (IC->dP[t] < 1.0e-12) {
		fprintf(stderr, "obs %d forced to point\n", t);
		IC->obstype[t] = INT_FPOINT;
		IC->dP[t] = normal_pdf(z0)/sigma;
		IC->f0[t] = IC->f1[t] = 0;
	    } else {
		IC->f0[t] = normal_pdf(z0) / IC->dP[t];
		IC->f1[t] = normal_pdf(z1) / IC->dP[t];
	    }
	    break;
	case INT_POINT:
	    z0 = (x0 - ndxt)/sigma;
	    IC->dP[t] = normal_pdf(z0)/sigma;
	    IC->f0[t] = IC->f1[t] = 0;
	}
    }
}

static double interval_loglik (const double *theta, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    double x0, x1, z0, z1, gti;
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
	case INT_FPOINT:
	    z0 = (x0 - ndxt)/sigma;
	    derivb = z0/sigma;
	    derivs = z0*z0 - 1;
	}

	ll += log(IC->dP[t]);

	for (i=0; i<IC->nx; i++) {
	    gti = derivb * gretl_matrix_get(IC->X, t, i);
	    gretl_matrix_set(IC->G, t, i, gti);
	    IC->g[i] += gti;
	}

	gretl_matrix_set(IC->G, t, k-1, derivs);
	IC->g[k-1] += derivs;
    }

#if BIGDEBUG
    fprintf(stderr, "*** interval loglik: ll = %g ***\n", ll);
#endif

#if INTDEBUG > 1
    fprintf(stderr, "ll = %16.10f\n", ll);

    for (i=0; i<IC->k; i++) {
	fprintf(stderr, "g[%d] = %16.10f\t", i, IC->g[i]);
    }
    fputc('\n', stderr);
#endif

    return ll;
}

static int interval_score (double *theta, double *s, int npar,
			   BFGS_CRIT_FUNC ll, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    int i;

    for (i=0; i<npar; i++) {
	s[i] = IC->g[i];
    }

    return 0;
}

int interval_hessian (double *theta, gretl_matrix *V, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    double z0 = 0, z1 = 0, mu = 0;
    double q0 = 0, q1 = 0, nu = 0;
    double x, vij, x0, x1;
    double f0, f1, Hss = 0;
    double sigma, ndxt, lambda = 0;
    int i, j, t, k = IC->k;
#if defined(_OPENMP)
    int parhess = (IC->nx > 400);
#endif
    int err = 0;

    sigma = exp(theta[k-1]);
    loglik_prelim(theta, IC);

    gretl_matrix_zero(V);

#if BIGDEBUG
    fprintf(stderr, "*** interval hessian ***\n");
#endif

#if INTDEBUG > 1
    fprintf(stderr, "interval_hessian:\n");
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
	case INT_FPOINT:
	    z0 = (x0 - ndxt)/sigma;
	}

	if ((IC->obstype[t] == INT_POINT) ||
	    (IC->obstype[t] == INT_FPOINT)) {
	    x = 1.0 / (sigma*sigma);
	} else {
	    x = (mu*mu - lambda/sigma);
	}

#if defined(_OPENMP)
#pragma omp parallel for private(j, i, x0, x1, vij) if (parhess)
#endif
	for (j=0; j<IC->nx; j++) {
	    x0 = gretl_matrix_get(IC->X, t, j);
	    for (i=j; i<IC->nx; i++) {
		x1 = gretl_matrix_get(IC->X, t, i);
		vij = gretl_matrix_get(V, j, i);
		vij += x * x0 * x1;
		gretl_matrix_set(V, j, i, vij);
	    }
	}

	if ((IC->obstype[t] == INT_POINT) ||
	    (IC->obstype[t] == INT_FPOINT)) {
	    x = 2 * z0 / sigma;
	} else {
	    x = (mu*lambda*sigma - nu);
	}

	for (i=0; i<IC->nx; i++) {
	    x0 = gretl_matrix_get(IC->X, t, i);
	    vij = gretl_matrix_get(V, i, k-1);
	    vij += x * x0;
	    gretl_matrix_set(V, i, k-1, vij);
	}

	if ((IC->obstype[t] == INT_POINT) ||
	    (IC->obstype[t] == INT_FPOINT)) {
	    x = 2 * (z0*z0);
	    Hss += x;
	} else {
	    x = lambda*sigma;
	    Hss += x*x - (f0 * q0 * z0 - f1 * q1 * z1);
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

static gretl_matrix *interval_hessian_inverse (double *theta,
					       void *ptr,
					       int *err)
{
    int_container *IC = ptr;
    gretl_matrix *H;

    H = gretl_zero_matrix_new(IC->k, IC->k);

    if (H == NULL) {
	*err = E_ALLOC;
    } else {
	*err = interval_hessian(theta, H, ptr);
    }

    if (!*err) {
	*err = gretl_invert_symmetric_matrix(H);
    }

    return H;
}

static int interval_trim_vcv (MODEL *pmod, int k,
			      const gretl_matrix *V)
{
    gretl_matrix *Vt;
    double vij;
    int n = V->rows;
    int i, j;

    Vt = gretl_matrix_copy(V);
    if (Vt == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_reuse(Vt, k, k);

    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    vij = pmod->vcv[ijton(i, j, n)];
	    gretl_matrix_set(Vt, i, j, vij);
	}
    }

    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    vij = gretl_matrix_get(Vt, i, j);
	    pmod->vcv[ijton(i, j, k)] = vij;
	}
    }

    gretl_matrix_free(Vt);

    return 0;
}

static int intreg_model_add_vcv (MODEL *pmod,
				 int_container *IC,
				 const DATASET *dset,
				 gretlopt opt,
				 double *x)
{
    gretl_matrix *H = NULL;
    gretl_matrix *V = NULL;
    int err = 0;

#if BIGDEBUG
    fprintf(stderr, "*** intreg: add vcv ***\n");
#endif

    if (opt & OPT_G) {
	err = gretl_model_add_OPG_vcv(pmod, IC->G, &V);
    } else {
	H = interval_hessian_inverse(IC->theta, IC, &err);
	if (!err) {
	    if (opt & OPT_R) {
		err = gretl_model_add_QML_vcv(pmod, pmod->ci, H, IC->G,
					      dset, opt, &V);
	    } else {
		V = H;
		err = gretl_model_add_hessian_vcv(pmod, H);
	    }
	}
    }

    if (!err) {
	int k = IC->k - 1;

	*x = gretl_model_get_vcv_element(pmod, k, k, V->rows);
	err = interval_trim_vcv(pmod, IC->nx, V);
	if (!err) {
	    /* note: donates @V to @pmod */
	    gretl_model_set_matrix_as_data(pmod, "full_vcv", V);
	} else if (V != H) {
	    gretl_matrix_free(V);
	}
    }

    if (opt & OPT_R) {
	/* In the robust case we're finished with H, which was
	   just one ingredient in the covariance matrix
	*/
	gretl_matrix_free(H);
    }

    return err;
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
	case INT_FPOINT:
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
    double gti, m3 = 0.0, m4 = 0.0;
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
	case INT_FPOINT:
	    u = IC->uhat[t] / sigma;
	    m3 = u*u*u;
	    m4 = m3*u;
	}

	m4 -= 3;
	*sm3 += m3;
	*sm4 += m4;

	for (i=0; i<k; i++) {
	    gti = gretl_matrix_get(IC->G, t, i);
	    gretl_matrix_set(ret, t, i, gti);
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

static int fill_intreg_model (int_container *IC,
			      int fncount, int grcount,
			      const DATASET *dset,
			      gretlopt opt)
{
    MODEL *pmod = IC->pmod;
    double ndx, u, x = 0;
    int i, j, k = IC->k;
    int obstype;
    int err = 0;

    pmod->ci = (opt & OPT_T)? TOBIT : INTREG;
    pmod->lnL = IC->ll;
    mle_criteria(pmod, 1);
    pmod->sigma = exp(IC->theta[k-1]);

    for (i=0; i<IC->nx; i++) {
	pmod->coeff[i] = IC->theta[i];
    }

    err = intreg_model_add_vcv(pmod, IC, dset, opt, &x);
    if (err) {
	goto bailout;
    }

    if (!na(x)) {
	/* get the s.e. of sigma via the delta method */
	 x *= pmod->sigma * pmod->sigma;
	 gretl_model_set_double(pmod, "se_sigma", sqrt(x));
    }

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
	gretl_model_set_int(pmod, "n_fpoint", IC->typecount[4]);
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

 bailout:

    if (err && pmod->errcode == 0) {
	pmod->errcode = err;
    }

    return err;
}

static int do_interval (int *list, DATASET *dset, MODEL *mod,
			gretlopt opt, PRN *prn)
{
    int_container *IC;
    int maxit, fncount = 0, grcount = 0;
    double toler, normtest = NADBL;
    gretlopt maxopt = (opt & OPT_V) | OPT_U;
    int optim, use_bfgs = 0;
    int err = 0;

    IC = int_container_new(list, dset, mod);
    if (IC == NULL) {
	return E_ALLOC;
    }

    BFGS_defaults(&maxit, &toler, INTREG);

    optim = libset_get_int(GRETL_OPTIM);
    if (optim == OPTIM_BFGS) {
	use_bfgs = 1;
    } else if (optim == OPTIM_AUTO && mod->ncoeff > 600) {
	/* avoid bogging down */
	use_bfgs = 1;
    }

#if INTDEBUG || BIGDEBUG
    fprintf(stderr, "do_interval: use_bfgs = %d\n", use_bfgs);
#endif

    if (use_bfgs) {
	err = BFGS_max(IC->theta, IC->k, maxit, toler,
		       &fncount, &grcount, interval_loglik, C_LOGLIK,
		       interval_score, IC, NULL, maxopt, prn);
    } else {
	double crittol = 1.0e-07;
	double gradtol = 1.0e-07;

	err = newton_raphson_max(IC->theta, IC->k, maxit, crittol, gradtol,
				 &fncount, C_LOGLIK, interval_loglik,
				 interval_score, interval_hessian,
				 IC, maxopt, prn);
    }

#if INTDEBUG || BIGDEBUG
    fprintf(stderr, "after maximization, err = %d\n", err);
#endif

    if (!err) {
	IC->ll = interval_loglik(IC->theta, IC);
    }

     if (!err) {
	err = fill_intreg_model(IC, fncount, grcount, dset, opt);
    }

    if (!err) {
	err = intreg_normtest(IC, &normtest);
	if (!err) {
	    err = gretl_model_add_normality_test(IC->pmod, normtest);
	}
	/* don't let this be a show-stopper */
	err = 0;
    }

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

MODEL interval_estimate (const int *list, DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    MODEL model;
    int *mylist = gretl_list_copy(list);
    int *initlist = NULL;

    gretl_model_init(&model, NULL);

    if (mylist[0] > 3) {
	maybe_reposition_const(mylist, dset);
    }

    /* create extra series for model initialization and
       corresponding regression list
    */
    model.errcode = create_midpoint_y(mylist, dset, &initlist);
    if (model.errcode) {
        free(mylist);
	return model;
    }

    /* run initial OLS */
    model = lsq(initlist, dset, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "interval_estimate: initial OLS failed\n");
	free(initlist);
        free(mylist);
	return model;
    }

#if INTDEBUG
    pprintf(prn, "interval_estimate: initial OLS\n");
    printmodel(&model, dset, OPT_S, prn);
#endif

    /* clean up midpoint-y */
    dataset_drop_last_variables(dset, 1);
    free(initlist);

    if (opt & OPT_C) {
	opt |= OPT_R;
    }

    /* do the actual analysis */
    model.errcode = do_interval(mylist, dset, &model, opt, prn);

    clear_model_xpx(&model);
    free(mylist);

    return model;
}

static int tobit_add_lo_hi (MODEL *pmod, double llim, double rlim,
			    DATASET *dset, int **plist)
{
    double *y, *lo, *hi;
    int *list;
    int lv, hv;
    int i, t, err;

    err = dataset_add_series(dset, 2);
    if (err) {
	return err;
    }

    if (na(llim)) {
	llim = -1.0e300;
    }
    if (na(rlim)) {
	rlim = 1.0e300;
    }

    lv = dset->v - 2;
    hv = dset->v - 1;
    lo = dset->Z[lv];
    hi = dset->Z[hv];
    y = dset->Z[pmod->list[1]];

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	if (na(y[t])) {
	    /* y missing */
	    lo[t] = hi[t] = NADBL;
	} else if (y[t] > llim && y[t] < rlim) {
	    /* y in uncensored region */
	    lo[t] = hi[t] = y[t];
	} else if (y[t] >= rlim) {
	    /* y at or above right bound */
	    lo[t] = rlim;
	    hi[t] = NADBL;
	} else if (y[t] <= llim) {
	    /* y at or below left bound */
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

static int tobit_depvar_check (const int *list)
{
    int i, yvar = list[1];

    for (i=2; i<=list[0]; i++) {
	if (list[i] == yvar) {
	    gretl_errmsg_set(_("tobit: the dependent variable cannot be "
			     "included as a regressor"));
	    return E_DATA;
	}
    }

    return 0;
}

MODEL tobit_via_intreg (int *list, double llim, double rlim,
			DATASET *dset, gretlopt opt,
			PRN *prn)
{
    MODEL model;
    int *ilist = NULL;
    int origv = dset->v;
    int err;

    err = tobit_depvar_check(list);
    if (err) {
	gretl_model_init(&model, NULL);
	model.errcode = err;
	return model;
    }

    /* run initial OLS */
    model = lsq(list, dset, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "intreg: initial OLS failed\n");
	return model;
    }

#if INTDEBUG
    pprintf(prn, "tobit: llim=%g, rlim=%g\n", llim, rlim);
    pprintf(prn, "tobit_via_intreg: initial OLS\n");
    printmodel(&model, dset, OPT_S, prn);
#endif

    model.errcode = tobit_add_lo_hi(&model, llim, rlim,
				    dset, &ilist);

    if (opt & OPT_C) {
	opt |= OPT_R;
    }

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
    dataset_drop_last_variables(dset, dset->v - origv);
    free(ilist);

    return model;
}
