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

/* The algorithm here was contributed by Riccardo "Jack" Lucchetti; 
   C coding by Allin Cottrell. 
*/

#include "libgretl.h"
#include "bhhh_max.h"
#include "libset.h"

#undef TDEBUG 

#define USE_BFGS

#define TOBIT_TOL 1.0e-10 /* calibrated against William Greene */

typedef struct tob_container_ tob_container;

struct tob_container_ {
    int k;
    int n;
    int do_score;
    double ll;

    const double **X;       /* data */
    double **G;
    double *g;
    double *theta;

    /* workspace */
    double *e;
    double *f;
    double *P;
    double *ystar;
};

static void tob_container_destroy (tob_container *TC)
{
    if (TC != NULL) {
	doubles_array_free(TC->G, TC->k);

	free(TC->g);
	free(TC->theta);
	free(TC->e);
	free(TC->f);
	free(TC->P);
	free(TC->ystar);

	free(TC);
    }
}

static tob_container *tob_container_new (int k, MODEL *pmod, const double **X, 
					 int do_score)
{
    tob_container *TC = malloc(sizeof *TC);
    int i, n = pmod->nobs;

    if (TC == NULL) {
	return NULL;
    }
    
    TC->G = NULL;
    TC->g = NULL;
    TC->theta = NULL;

    TC->e = NULL;
    TC->f = NULL;
    TC->P = NULL;
    TC->ystar = NULL;

    TC->k = k;
    TC->n = n;
    TC->X = X;

    TC->do_score = do_score;

    TC->G = doubles_array_new(k, n);
    TC->g = malloc(k * sizeof *TC->g);
    TC->theta = malloc(k * sizeof *TC->theta);

    TC->e = malloc(n * sizeof *TC->e);
    TC->f = malloc(n * sizeof *TC->f);
    TC->P = malloc(n * sizeof *TC->P);
    TC->ystar = malloc(n * sizeof *TC->ystar);

    if (TC->G == NULL || TC->g == NULL || TC->theta == NULL ||
	TC->e == NULL || TC->f == NULL || TC->P == NULL ||
	TC->ystar == NULL) {
	tob_container_destroy(TC);
	TC = NULL;
    } else {
	for (i=0; i<k; i++) {
	    TC->theta[i] = pmod->coeff[i];
	}
    }

    return TC;
}

#ifndef USE_BFGS

/* BHHH: Compute log likelihood, and the score matrix if do_score is
   non-zero */

static int tobit_ll (double *theta, const double **X, double **Z, 
		     model_info *tobit, int do_score)
{

    const double *y = X[1];
    double **series = model_info_get_series(tobit);
    double *e = series[0];
    double *f = series[1];
    double *P = series[2];
    double *ystar = series[3];
    int k = model_info_get_k(tobit);
    int n = model_info_get_n(tobit);
    
    double siginv = theta[k-1]; /* inverse of variance */
    double ll, llt;
    int i, t;

    if (siginv < 0.0) {
	fprintf(stderr, "tobit_ll: got a negative variance\n");
	return 1;
    } 

    /* calculate ystar, e, f, and P vectors */
    for (t=0; t<n; t++) {
	ystar[t] = theta[0];
	for (i=1; i<k-1; i++) {
	    ystar[t] += theta[i] * X[i+1][t]; /* coeff * orig data */
	}
	e[t] = y[t] * siginv - ystar[t];
	f[t] = siginv * normal_pdf(e[t]);
	P[t] = normal_cdf(ystar[t]);
#ifdef TDEBUG
	fprintf(stderr, "tobit_ll: e[%d]=%g, f[%d]=%g, P[%d]=%g\n",
		t, e[t], t, f[t], t, P[t]);
#endif
    }

    /* compute loglikelihood for each obs, cumulate into ll */
    ll = 0.0;
    for (t=0; t<n; t++) {
	if (y[t] == 0.0) {
	    llt = 1.0 - P[t];
	} else {
	    llt = f[t];
	}
	if (llt == 0.0) {
	    fprintf(stderr, "tobit_ll: L[%d] is zero\n", t);
	    return 1;
	}
	ll += log(llt);
    }

    model_info_set_ll(tobit, ll, do_score);

    if (do_score) {
	int i, gi, xi;
	double den, tail;

	for (t=0; t<n; t++) {
	    for (i=0; i<k; i++) {

		/* set the indices into the data arrays */
		gi = i + 1;
		xi = (i == 0)? 0 : i + 1;

		if (y[t] == 0.0) {
		    /* score if y is censored */
		    if (xi < k) {
			den = normal_pdf(ystar[t]);
			tail = 1.0 - P[t];
			Z[gi][t] = -den / tail * X[xi][t];
		    } else {
			Z[gi][t] = 0.0;
		    } 
		} else {
		    /* score if y is not censored */
		    if (xi < k) {
			Z[gi][t] = e[t] * X[xi][t];
		    } else {
			Z[gi][t] = e[t] * -y[t];
		    }
		    if (xi == k) {
			Z[gi][t] += 1.0 / siginv;
		    }
		}
	    }
	}
    }

    return 0;
}

#endif

static double t_loglik (const double *theta, void *ptr)
{
    double ll = NADBL;
    
    tob_container *TC = (tob_container *) ptr;
    const double **X = TC->X;
    const double *y = X[1];
    int k = TC->k;
    int n = TC->n;
    int do_score = TC->do_score;
    double siginv = theta[k-1];  /* inverse of variance */

    double *e = TC->e;
    double *f = TC->f;
    double *P = TC->P;
    double *ystar = TC->ystar;
    double llt;
    int i, t;

    if (siginv < 0.0) {
	fprintf(stderr, "tobit_ll: got a negative variance\n");
	return NADBL;
    } 

    /* calculate ystar, e, f, and P vectors */
    for (t=0; t<n; t++) {
	ystar[t] = theta[0];
	for (i=1; i<k-1; i++) {
	    ystar[t] += theta[i] * X[i+1][t]; /* coeff * orig data */
	}
	e[t] = y[t] * siginv - ystar[t];
	f[t] = siginv * normal_pdf(e[t]);
	P[t] = normal_cdf(ystar[t]);
#if 0
	fprintf(stderr, "t_loglik: e[%d]=%g, f[%d]=%g, P[%d]=%g\n",
		t, e[t], t, f[t], t, P[t]);
#endif
    }

    /* compute loglikelihood for each obs, cumulate into ll */
    ll = 0.0;
    for (t=0; t<n; t++) {
	if (y[t] == 0.0) {
	    llt = 1.0 - P[t];
	} else {
	    llt = f[t];
	}
	if (llt == 0.0) {
#if 0
	    fprintf(stderr, "tobit_ll: L[%d] is zero\n", t);
#endif
	    return NADBL;
	}
	ll += log(llt);
    }

    if (do_score) {
	double *score = TC->g;
	double **Z = TC->G;

	int i, gi, xi;
	double den, tail;

	for (i=0; i<k; i++) {
	    score[i] = 0.0;
	}

	for (t=0; t<n; t++) {
	    for (i=0; i<k; i++) {
		/* set the indices into the data arrays */
		gi = i + 1;
		xi = (i == 0)? 0 : i + 1;

		if (y[t] == 0.0) {
		    /* score if y is censored */
		    if (i < (k-1)) {
			den = normal_pdf(ystar[t]);
			tail = 1.0 - P[t];
			Z[i][t] = -den / tail * X[xi][t];
			score[i] += Z[i][t];
		    } else {
			Z[i][t] = 0.0;
		    } 
		} else {
		    /* score if y is not censored */
		    if (i < (k-1)) {
			Z[i][t] = e[t] * X[xi][t];
		    } else {
			Z[i][t] = e[t] * -y[t];
		    }
		    if (i == (k-1)) {
			Z[i][t] += 1.0 / siginv;
		    }
		    score[i] += Z[i][t];
		}
		
	    }
	}
    }

    TC->ll = ll;

    return ll;
}

static int t_score (double *theta, double *s, int npar, BFGS_LL_FUNC ll, 
		    void *ptr)
{
    tob_container *TC = (tob_container *) ptr;
    int i;

    for (i=0; i<npar; i++) {
	s[i] = TC->g[i];
    }

    return 1;
}

static gretl_matrix *tobit_opg_vcv (const tob_container *TC) 
{

    int n = TC->n;
    int k = TC->k;
    double x;
    int i, j, t;
    double **G = TC->G;
    gretl_matrix *V;

    V = gretl_matrix_alloc(k,k);
    if (V == NULL) {
	return NULL;
    }

    for (i=0; i<k; i++) {
	for (j=i; j<k; j++) {
	    x = 0.0;
	    for (t=0; t<n; t++) {
		x += G[i][t] * G[j][t];
	    }
	    gretl_matrix_set(V, i, j, x);
	    gretl_matrix_set(V, j, i, x);
	}
    }

    gretl_invert_symmetric_matrix(V);

    return V;
}

/* Transcribe the VCV matrix into packed triangular form */

static int pack_vcv (MODEL *pmod, gretl_matrix *v, double scale)
{
    const int nv = pmod->ncoeff;
    const int nterms = nv * (nv + 1) / 2;
    double x;
    int i, j, k;

    if (pmod->vcv == NULL) {
	pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
	if (pmod->vcv == NULL) {
	    return E_ALLOC;
	}
    } 

    for (i=0; i<nv; i++) {
	for (j=0; j<=i; j++) {
	    k = ijton(i, j, nv);
	    x = gretl_matrix_get(v, i, j);
	    pmod->vcv[k] = x;
	    if (scale != 1.0) {
		pmod->vcv[k] /= scale * scale;
	    }
	}
    }

    return 0;
}

static double recompute_tobit_ll (const MODEL *pmod, const double *y)
{
    double lt, ll = 0.0;
    int t, s = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	if (y[s++] == 0.0) {
	    lt = normal_cdf(-pmod->yhat[t] / pmod->sigma);
	} else {
	    lt = (1.0 / pmod->sigma) * normal_pdf(pmod->uhat[t] / pmod->sigma);
	}
	ll += log(lt);
    }

    return ll;
}

static int add_norm_test_to_model (MODEL *pmod, double X2)
{
    ModelTest *test = model_test_new(GRETL_TEST_NORMAL);
    int err = 0;

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_NORMAL_CHISQ);
        model_test_set_dfn(test, 2);
        model_test_set_value(test, X2);
        model_test_set_pvalue(test, chisq_cdf_comp(X2, 2));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = 1;
    }

    return err;
}

/* Chesher-Irish test for normality based on the generalized
   residuals.  See Chesher and Irish, "Residual analysis in the
   grouped and censored normal linear model," Journal of Econometrics,
   vol. 34(1-2), 1987, pp. 33-61.  Note that the presentation of
   this test in William Greene's "Econometric Analysis", 4e, is
   badly broken.
*/

static int chesher_irish_test (MODEL *pmod, const double **X)
{
    double **cZ = NULL;
    DATAINFO *cinfo;
    MODEL mod;
    const double *y, *e, *Xb;
    double s, s2, e2, li;
    double es, ndxs, ndxs2;
    int *list;
    int i, t, k, nv, uncens;
    int err = 0;

    k = pmod->ncoeff;
    nv = 4 + k;

    list = gretl_list_new(nv);
    if (list == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nv; i++) {
	list[i+1] = i;
    }

    cinfo = create_new_dataset(&cZ, nv, pmod->nobs, 0);
    if (cinfo == NULL) {
	free(list);
	return E_ALLOC;
    }

    e = pmod->uhat + pmod->t1;
    Xb = pmod->yhat + pmod->t1;
    y = X[1];
    s = pmod->sigma;
    s2 = s * s;

    for (t=0; t<cinfo->n; t++) {
	es = e[t] / s;
	ndxs = Xb[t] / s;
	e2 = es * es;
	li = normal_pdf(ndxs) / (1 - normal_cdf(ndxs));
	ndxs2 = ndxs * ndxs;
	uncens = (y[t] > 0);

	for (i=1; i<nv; i++) {
	    if (i == 1) {
		cZ[i][t] = e[t];
	    } else if (i <= k) {
		cZ[i][t] = e[t] * X[i][t];
	    } else if (i == k + 1) {
		if (uncens) {
		    /* use the generalized residual here */
		    cZ[i][t] = (es * es - 1) / (2 * s2);
		} else {
		    cZ[i][t] = ndxs * li / (2 * s2);
		}
	    } else {
		if (uncens) {
		    if (i == k + 2) {
			/* skewness */
			cZ[i][t] = es * e2;
		    } else {
			/* excess kurtosis */
			cZ[i][t] = e2 * e2 - 3;
		    }
		} else {
		    if (i == k + 2) {
			/* skewness */
			cZ[i][t] = es * (2 + ndxs2);

		    } else {
			/* excess kurtosis */
			cZ[i][t] = li * ndxs * (3 + ndxs2);
		    }
		}
	    }
	}
    }

    mod = lsq(list, &cZ, cinfo, OLS, OPT_A);
    if (!mod.errcode) {
	add_norm_test_to_model(pmod, mod.nobs - mod.ess);
    }

    clear_model(&mod);
    free(list);
    destroy_dataset(cZ, cinfo);
    
    return err;
}

static double gen_res (double yt, double ndxt, double sig)
{
    double ret;

    if (yt > 0.0) {
	ret = yt - ndxt;
    } else {
	double std = ndxt / sig;
	double F = 1.0 - normal_cdf(std);

	ret = -sig * normal_pdf(std) / F;
    }

#if 0 /* do we need/want this? */
    ret /= sig * sig;
#endif

    return ret;
}

/* Taking the original OLS model as a basis, re-write the statistics
   to reflect the Tobit results.
*/

static int write_tobit_stats (MODEL *pmod, double *theta, int ncoeff,
			      double sigma, double ll, const double **X,
			      gretl_matrix *VCV, double scale, int iters)
{
    int i, t, s, cenc = 0;
    const double *y = X[1];

    for (i=0; i<ncoeff; i++) {
	pmod->coeff[i] = theta[i];
	pmod->sderr[i] = sqrt(gretl_matrix_get(VCV, i, i));
	if (scale != 1.0) {
	    pmod->coeff[i] /= scale;
	    pmod->sderr[i] /= scale;
	}
    }

    pmod->sigma = sigma;

    if (scale != 1.0) {
	pmod->sigma /= scale;	
	pmod->ybar /= scale;
	pmod->sdy /= scale;
    }

    pmod->ess = 0.0;
    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	double yt;

	if (na(pmod->uhat[t])) {
	    continue;
	}

	pmod->yhat[t] = pmod->coeff[0];
	for (i=1; i<ncoeff; i++) {
	    pmod->yhat[t] += pmod->coeff[i] * X[i + 1][s];
	}

	yt = y[s];
	if (yt > 0 && scale != 1.0) {
	    yt /= scale;
	}

	pmod->uhat[t] = gen_res(yt, pmod->yhat[t], pmod->sigma);
	pmod->ess += pmod->uhat[t] * pmod->uhat[t]; /* Is this meaningful? */

	if (yt == 0.0) cenc++;
	s++;
    }

    if (scale != 1.0) {
	pmod->lnL = recompute_tobit_ll(pmod, y);
    } else {
	pmod->lnL = ll;
    }

    chesher_irish_test(pmod, X);
    pmod->fstt = pmod->rsq = pmod->adjrsq = NADBL;
    mle_criteria(pmod, 1);
    pack_vcv(pmod, VCV, scale);
    pmod->ci = TOBIT;

    gretl_model_set_int(pmod, "censobs", cenc);
    gretl_model_set_int(pmod, "iters", iters);

    return 0;
}

#ifndef USE_BFGS

static model_info *
tobit_model_info_init (int model_obs, int bign, int k, int n_series)
{
    double tol = get_bhhh_toler();
    model_info *tobit;

    if (na(tol)) {
	tol = TOBIT_TOL;
    }

    tobit = model_info_new(k, 0, model_obs - 1, bign, tol);

    if (tobit != NULL) {
	model_info_set_opts(tobit, FULL_VCV_MATRIX);
	model_info_set_n_series(tobit, n_series);
    }

    return tobit;
}

#endif

static int transform_tobit_VCV (gretl_matrix *VCV, double sigma,
				const double *theta, int k)
{
    gretl_matrix *J = NULL;
    gretl_matrix *kk = NULL;
    int i, j;

    /* Jacobian mat. for transforming VCV from Olsen to slopes + variance */
    J = gretl_zero_matrix_new(k, k);
    if (J == NULL) {
	return E_ALLOC;
    }

    kk = gretl_matrix_copy(VCV);
    if (kk == NULL) {
	gretl_matrix_free(J);
	return E_ALLOC;
    }    

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    if (i == j && i < k-1) {
		/* upper left diagonal component */
		gretl_matrix_set(J, i, j, sigma);
	    } else if (j == k-1 && i < j) {
		/* right-hand column */
		gretl_matrix_set(J, i, j, -sigma * theta[i]);
	    } else if (j == k-1 && i == j) {
		/* bottom right-hand element */
		gretl_matrix_set(J, i, j, -sigma * sigma);
	    }
	}
    }

    gretl_matrix_qform(J, GRETL_MOD_NONE, kk,
		       VCV, GRETL_MOD_NONE);

    gretl_matrix_free(J);
    gretl_matrix_free(kk);

    return 0;
}

static int tobit_add_variance (MODEL *pmod, int k)
{
    double *b;

    b = realloc(pmod->coeff, k * sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    /* initialize variance */
    b[k-1] = 1.0;
    pmod->coeff = b;

    return 0;
}

/* Main Tobit function: two variants */

#ifdef USE_BFGS

static int do_tobit (double **Z, DATAINFO *pdinfo, MODEL *pmod,
		     double scale, int missvals, PRN *prn)
{
    tob_container *TC = NULL;
    int fncount, grcount;
    double **X;
    gretl_matrix *VCV = NULL;
    double sigma;
    int i, k, n;
    int nv = pmod->list[0];
    int err = 0;

    k = pmod->ncoeff + 1; /* add the variance */
    n = pmod->nobs;

    /* set of pointers into original data, or a reduced copy 
       if need be */
    X = data_array_from_model(pmod, Z, missvals);
    if (X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = tobit_add_variance(pmod, k);
    if (err) {
	goto bailout;
    }

    TC = tob_container_new(k, pmod, (const double **) X, 1);
    if (TC == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = BFGS_max(TC->theta, TC->k, 1000, TOBIT_TOL, 
		   &fncount, &grcount, t_loglik, 
		   (1 ? t_score : NULL), 
		   TC, (prn != NULL)? OPT_V : OPT_NONE,
		   prn);
    if (err) {
	goto bailout;
    }

    /* recover estimate of variance */
    sigma = 1.0 / TC->theta[k-1]; 

    /* recover slope estimates */
    for (i=0; i<k-1; i++) {
	TC->theta[i] *= sigma;
    }

    /* get estimate of variance matrix for Olsen parameters */
    VCV = tobit_opg_vcv(TC);

    if (VCV == NULL) {
	err = E_ALLOC;
    } else {
	/* Do Jacobian thing */
	err = transform_tobit_VCV(VCV, sigma, TC->theta, k);
    }

    if (!err) {
	write_tobit_stats(pmod, TC->theta, k-1, sigma, TC->ll, (const double **) X, 
			  VCV, scale, fncount);
    }

 bailout:

    if (missvals) {
	doubles_array_free(X, nv);
    } else {
	free(X);
    }

    tob_container_destroy(TC);

    return err;
}

#else /* BHHH */

static int do_tobit (double **Z, DATAINFO *pdinfo, MODEL *pmod,
		     double scale, int missv, PRN *prn)
{
    double **X;
    double *theta = NULL;
    gretl_matrix *VCV = NULL;
    double sigma;
    int i, k, n;
    int nv = pmod->list[0];
    int err = 0;

    model_info *tobit = NULL;
    int n_series = 4;
    double ll;

    k = pmod->ncoeff + 1; /* add the variance */
    n = pmod->nobs;

    /* set of pointers into original data, or a reduced copy 
       if need be */
    X = data_array_from_model(pmod, Z, missv);
    if (X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = tobit_add_variance(pmod, k);
    if (err) {
	goto bailout;
    }

    tobit = tobit_model_info_init(pmod->nobs, pdinfo->n, k, n_series);
    if (tobit == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = bhhh_max(tobit_ll, (const double **) X, pmod->coeff, tobit, prn);
    if (err) {
	goto bailout;
    }

    theta = model_info_get_theta(tobit);

    /* recover estimate of variance */
    sigma = 1.0 / theta[k-1]; 

    /* recover slope estimates */
    for (i=0; i<k-1; i++) {
	theta[i] *= sigma;
    }

    /* get estimate of variance matrix for Olsen parameters */
    VCV = model_info_get_VCV(tobit);
    gretl_invert_symmetric_matrix(VCV);
    gretl_matrix_divide_by_scalar(VCV, n);

    /* Do Jacobian thing */
    err = transform_tobit_VCV(VCV, sigma, theta, k);

    ll = model_info_get_ll(tobit);
    write_tobit_stats(pmod, theta, k-1, sigma, ll, (const double **) X, 
		      VCV, scale, model_info_get_iters(tobit));

 bailout:

    if (missv) {
	doubles_array_free(X, nv);
    } else {
	free(X);
    }

    if (tobit != NULL) {
	model_info_free(tobit);
    }

    return err;
}

#endif /* BFGS vs BHHH */

static double tobit_depvar_scale (const MODEL *pmod, int *miss)
{
    double ut, umax = 0.0;
    double scale = 1.0;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    *miss = 1;
	} else {
	    ut = fabs(pmod->uhat[t]);
	    if (ut > umax) {
		umax = ut;
	    }
	}
    }

    if (umax > 5.0) {
	/* normal pdf will lose accuracy beyond, say, z = 5 */
	scale = 5.0 / umax;
    }

    return scale;
}

/* the driver function for the plugin */

MODEL tobit_estimate (const int *list, double ***pZ, DATAINFO *pdinfo,
		      PRN *prn) 
{
    MODEL model;
    double *y;
    double scale = 1.0;
    int missvals = 0;
    int t;

    /* run initial OLS: OPT_M would ban missing obs */
    model = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (model.errcode) {
	return model;
    }

    /* handle scale issues */
    scale = tobit_depvar_scale(&model, &missvals);
    if (scale != 1.0) {
	y = (*pZ)[model.list[1]];
	for (t=0; t<pdinfo->n; t++) {
	    if (!na(y[t])) {
		y[t] *= scale;
	    }
	}
	clear_model(&model);
	model = lsq(list, pZ, pdinfo, OLS, OPT_A);
    }

#ifdef TDEBUG
    pprintf(prn, "tobit_estimate: initial OLS\n");
    printmodel(&model, pdinfo, OPT_NONE, prn);
#endif

    /* do the actual Tobit analysis */
    if (model.errcode == 0) {
	model.errcode = do_tobit(*pZ, pdinfo, &model, scale, 
				 missvals, prn);
    }

    if (scale != 1.0) {
	y = (*pZ)[model.list[1]];
	for (t=0; t<pdinfo->n; t++) {
	    if (!na(y[t])) {
		y[t] /= scale;
	    }
	}
    }

    return model;
}
