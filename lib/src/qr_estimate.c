/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "qr_estimate.h"
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"
#include "internal.h"

#define QR_RCOND_MIN 1e-15 /* experiment with this? */
#define QR_SMALL 1e-15     /* SSR < this counts as zero */

/* In fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are
   contiguous.
*/

static double get_tss (const double *y, int n, int ifc)
{
    double ymean = 0.0;
    double x, tss = 0.0;
    int i;

    if (ifc) ymean = _esl_mean(0, n-1, y);

    for (i=0; i<n; i++) {
	x = y[i] - ymean;
	tss += x * x;
    }

    return tss;
}

static void qr_compute_r_squared (MODEL *pmod, const double *y, int n)
{
    int t1 = pmod->t1;

    if (pmod->rho) t1++;

    if (pmod->dfd > 0) {
	if (pmod->ifc) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	    pmod->adjrsq = 1 - (pmod->ess * (n - 1) / den);
	} else {
	    double alt = corrrsq(n, y + t1, pmod->yhat + t1);

	    if (na(alt)) {
		pmod->rsq = pmod->adjrsq = NADBL;
	    } else {
		pmod->rsq = alt;
		pmod->adjrsq = 1 - ((1 - alt) * (n - 1) / pmod->dfd);
	    }
	}
	pmod->fstt = (pmod->tss - pmod->ess) * pmod->dfd / 
	    (pmod->ess * pmod->dfn);
    } else {
	pmod->rsq = 1.0;
	pmod->fstt = NADBL;
    }
}

static int get_vcv_index (MODEL *pmod, int i, int j, int n)
{
    int k, vi, vj;

    if (pmod->ifc) {
	vi = 1 + (i + 1) % n;
	vj = 1 + (j + 1) % n;
    } else {
	vi = i + 1;
	vj = j + 1;
    }

    k = ijton(vi, vj, n);

    return k;
}

static int qr_make_vcv (MODEL *pmod, gretl_matrix *v, int robust)
{
    const int nv = pmod->ncoeff;
    const int nterms = nv * (nv + 1) / 2;
    double x;
    int i, j, k;

    pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) return 1;  

    for (i=0; i<nv; i++) {
	for (j=0; j<=i; j++) {
	    k = get_vcv_index(pmod, i, j, nv);
	    x = gretl_matrix_get(v, i, j);
	    if (!robust) {
		x *= pmod->sigma * pmod->sigma;
	    }
	    pmod->vcv[k] = x;
	}
    }

    return 0;
}

static void get_resids_and_SSR (MODEL *pmod, const double **Z,
				gretl_matrix *y, double ypy, 
				int fulln)
{
    int t, i = 0;
    int dwt = gretl_model_get_int(pmod, "wt_dummy");

    if (dwt) dwt = pmod->nwt;

    pmod->ess = 0.0;

    if (pmod->rho) {
	for (t=0; t<fulln; t++) {
	    if (t <= pmod->t1 || t > pmod->t2) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    } else {
		double x = Z[pmod->list[1]][t];

		x -= pmod->rho * Z[pmod->list[1]][t-1];
		x -= y->val[i];
		pmod->uhat[t] = x;
		pmod->ess += x * x;
		i++;
	    }
	}
    } else if (pmod->nwt) {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    } else {
		double x = Z[pmod->list[1]][t];

		if (dwt && Z[dwt][t] == 0.0) {
		    pmod->yhat[t] = NADBL;
		} else {
		    if (!dwt) x *= Z[pmod->nwt][t];
		    pmod->yhat[t] = y->val[i];
		    pmod->uhat[t] = x - y->val[i];
		    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
		    i++;
		}
	    }
	}
    } else {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    } else {
		pmod->yhat[t] = y->val[i];
		pmod->uhat[t] = Z[pmod->list[1]][t] - y->val[i];
		pmod->ess += pmod->uhat[t] * pmod->uhat[t];
		i++;
	    }
	}
    }

    /* if SSR is small enough, treat it as zero */
    if (pmod->ess < QR_SMALL && pmod->ess > (-QR_SMALL)) {
	pmod->ess = 0.0;
    } 
}

static gretl_matrix *make_data_X (const MODEL *pmod, const double **Z)
{
    gretl_matrix *X;
    int i, t, j = 0;
    int start = (pmod->ifc)? 3 : 2;

    X = gretl_matrix_alloc(pmod->nobs, pmod->ncoeff);
    if (X == NULL) return NULL;

    /* copy independent vars into matrix X */
    for (i=start; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    X->val[j++] = Z[pmod->list[i]][t];
	}
    }

    /* insert constant */
    if (pmod->ifc) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    X->val[j++] = 1.0;
	}
    }

    return X;
}

/* Calculate W(t)-transpose * W(t-lag) */

static void wtw (gretl_matrix *wt, gretl_matrix *X, 
		 int n, int t, int lag)
{
    int i, j;
    double xi, xj;

    for (i=0; i<n; i++) {
	xi = X->val[mdx(X, t, i)];
	for (j=0; j<n; j++) {
	    xj = X->val[mdx(X, t - lag, j)];
	    wt->val[mdx(wt, i, j)] = xi * xj;
	}
    }
}

/* Calculate the Newey-West HAC covariance matrix.  Algorithm and
   (basically) notation taken from Davidson and MacKinnon (DM), 
   Econometric Theory and Methods, chapter 9.
*/

static int qr_make_hac (MODEL *pmod, const double **Z, 
			gretl_matrix *xpxinv)
{
    gretl_matrix *vcv = NULL, *wtj = NULL, *gammaj = NULL;
    gretl_matrix *X;
    int m = pmod->nobs;
    int n = pmod->ncoeff;
    int p, i, j, t;
    double weight, uu;
    double *uhat = pmod->uhat + pmod->t1;
    int err = 0;

    X = make_data_X(pmod, Z);
    if (X == NULL) return 1;

    /* get the user's preferred maximum lag setting */
    p = get_hac_lag(m);
    gretl_model_set_int(pmod, "hac_lag", p);

    vcv = gretl_matrix_alloc(n, n);
    wtj = gretl_matrix_alloc(n, n);
    gammaj = gretl_matrix_alloc(n, n);
    if (vcv == NULL || wtj == NULL || gammaj == NULL) {
	err = 1;
	goto bailout;
    }

    gretl_matrix_zero(vcv);

    for (j=0; j<=p; j++) {
	/* cumulate running sum of Gamma-hat terms */
	gretl_matrix_zero(gammaj);
	for (t=j; t<m; t++) {
	    /* W(t)-transpose * W(t-j) */
	    wtw(wtj, X, n, t, j);
	    uu = uhat[t] * uhat[t-j];
	    gretl_matrix_multiply_by_scalar(wtj, uu);
	    /* DM equation (9.36), p. 363 */
	    gretl_matrix_add_to(gammaj, wtj);
	}

	if (j > 0) {
	    /* Gamma(j) = Gamma(j) + Gamma(j)-transpose */
	    gretl_matrix_add_self_transpose(gammaj);
	    weight = 1.0 - (double) j / (p + 1.0);
	    /* multiply by Newey-West weight */
	    gretl_matrix_multiply_by_scalar(gammaj, weight);
	}

	/* DM equation (9.38), p. 364 */
	gretl_matrix_add_to(vcv, gammaj);
    }

    gretl_matrix_multiply_mod(xpxinv, GRETL_MOD_TRANSPOSE,
			      vcv, GRETL_MOD_NONE,
			      wtj);
    gretl_matrix_multiply(wtj, xpxinv, vcv);

    /* vcv now holds HAC */
    for (i=0; i<n; i++) {
	double x = gretl_matrix_get(vcv, i, i);

	j = (pmod->ifc)? (i + 1) % n : i;
	pmod->sderr[j] = sqrt(x);
    }

    /* Transcribe vcv into triangular representation */
    err = qr_make_vcv(pmod, vcv, 1);

 bailout:

    gretl_matrix_free(wtj);
    gretl_matrix_free(gammaj);
    gretl_matrix_free(vcv);
    gretl_matrix_free(X);

    return err;
}

static int qr_make_hccme (MODEL *pmod, const double **Z, 
			  gretl_matrix *xpxinv)
{
    gretl_matrix *X;
    gretl_matrix *diag = NULL;
    gretl_matrix *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL;
    int m = pmod->nobs; 
    int n = pmod->list[0] - 1;
    int i, j = 0;
    int err = 0;

    X = make_data_X(pmod, Z);
    if (X == NULL) return 1;

    diag = gretl_diagonal_matrix(pmod->uhat, m, GRETL_MOD_SQUARE);
    if (diag == NULL) {
	err = 1;
	goto bailout;
    }   

    tmp1 = gretl_matrix_alloc(n, m);
    tmp2 = gretl_matrix_alloc(n, n);
    tmp3 = gretl_matrix_alloc(n, n);
    if (tmp1 == NULL || tmp2 == NULL || tmp3 == NULL) {
	err = 1;
	goto bailout;
    }  

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
			      diag, GRETL_MOD_NONE,
			      tmp1);
    gretl_matrix_multiply(tmp1, X, tmp2);
    gretl_matrix_multiply(xpxinv, tmp2, tmp3); 
    gretl_matrix_multiply(tmp3, xpxinv, tmp2);

    /* tmp2 now holds HCCM */
    for (i=0; i<n; i++) {
	double x = tmp2->val[mdx(tmp2, i, i)];

	j = (pmod->ifc)? (i + 1) % n : i;
	pmod->sderr[j] = sqrt(x);
    }

    err = qr_make_vcv(pmod, tmp2, 1);

    bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(diag);
    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);
    gretl_matrix_free(tmp3);

    return err;
}

static int qr_make_regular_vcv (MODEL *pmod, gretl_matrix *v)
{
    int i, k, n = pmod->ncoeff;
    int err = 0;

    for (i=0; i<n; i++) {
	double x = v->val[mdx(v, i, i)];

	k = (pmod->ifc)? (i + 1) % n : i;
	pmod->sderr[k] = pmod->sigma * sqrt(x);
    }

    err = qr_make_vcv(pmod, v, 0);

    return err;
}

static double get_model_data (MODEL *pmod, const double **Z, 
			      gretl_matrix *Q, gretl_matrix *y)
{
    int i, j, t, start;
    double x, ypy = 0.0;
    int t1 = pmod->t1;
    int dwt = gretl_model_get_int(pmod, "wt_dummy");

    if (pmod->rho) t1++;

    start = (pmod->ifc)? 3 : 2;

    if (dwt) dwt = pmod->nwt;

    /* copy independent vars into matrix Q */
    j = 0;
    for (i=start; i<=pmod->list[0]; i++) {
	for (t=t1; t<=pmod->t2; t++) {
	    x = Z[pmod->list[i]][t];
	    if (dwt) {
		if (Z[dwt][t] == 0.0) continue;
	    } else if (pmod->nwt) {
		x *= Z[pmod->nwt][t];
	    } else if (pmod->rho && pmod->list[i] != 0) {
		x -= pmod->rho * Z[pmod->list[i]][t-1];
	    }
	    Q->val[j++] = x;
	}
    }

    /* insert constant last (numerical issues) */
    if (pmod->ifc) {
	for (t=t1; t<=pmod->t2; t++) {
	    if (dwt) {
		if (Z[dwt][t] == 0.0) continue;
	    } if (pmod->nwt) {
		Q->val[j++] = Z[pmod->nwt][t];
	    } else {
		Q->val[j++] = 1.0;
	    }
	}
    }

    /* copy dependent variable into y vector */
    j = 0;
    for (t=t1; t<=pmod->t2; t++) {
	x = Z[pmod->list[1]][t];
	if (dwt) {
	    if (Z[dwt][t] == 0.0) continue;
	} else if (pmod->nwt) {
	    x *= Z[pmod->nwt][t];
	} else if (pmod->rho) {
	    x -= pmod->rho * Z[pmod->list[1]][t-1];
	}
	y->val[j++] = x;
	ypy += x * x;
    }

    /* fetch tss based on the (possibly transformed) y */
    pmod->tss = get_tss(y->val, pmod->nobs, pmod->ifc);

    return ypy;
}

static void save_coefficients (MODEL *pmod, gretl_matrix *b,
			       int n)
{
    pmod->coeff = gretl_matrix_steal_data(b);

    if (pmod->ifc) {
	int i;
	double tmp = pmod->coeff[n-1];

	for (i=n-1; i>0; i--) pmod->coeff[i] = pmod->coeff[i-1];
	pmod->coeff[0] = tmp;
    }
}

int gretl_qr_regress (MODEL *pmod, const double **Z, int fulln,
		      unsigned long opts)
{
    integer info, lwork;
    integer m, n, lda;
    integer *iwork;
    gretl_matrix *Q, *y;
    gretl_matrix *R = NULL, *g = NULL, *b = NULL;
    gretl_matrix *xpxinv = NULL;
    doublereal *tau, *work;
    doublereal rcond, ypy;
    char uplo = 'U';
    char diag = 'N';
    char norm = '1';
    int i, j;
    int err = 0;

    m = pmod->nobs;               /* # of rows = # of observations */
    lda = m;                      /* leading dimension of Q */
    n = pmod->list[0] - 1;        /* # of cols = # of variables */

    Q = gretl_matrix_alloc(m, n);
    y = gretl_matrix_alloc(m, 1);

    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);

    work = malloc(sizeof *work);
    iwork = malloc(n * sizeof *iwork);

    if (Q == NULL || y == NULL || tau == NULL || work == NULL ||
	iwork == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    ypy = get_model_data(pmod, Z, Q, y);

    /* do a workspace size query */
    lwork = -1;
    info = 0;
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* set up optimally sized work array */
    lwork = (integer) work[0];
    work = realloc(work, (size_t) lwork * sizeof *work);
    if (work == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* run actual QR factorization */
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* check reciprocal condition number of R */
    dtrcon_(&norm, &uplo, &diag, &n, Q->val, &lda, &rcond, work, 
	    iwork, &info);
    if (info != 0) {
	fprintf(stderr, "dtrcon: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    if (rcond < QR_RCOND_MIN) {
	fprintf(stderr, "dtrcon: rcond = %g, but min is %g\n", rcond,
		QR_RCOND_MIN);
	err = E_SINGULAR;
	goto qr_cleanup;
    }

    /* allocate temporary auxiliary matrices */
    R = gretl_matrix_alloc(n, n);
    g = gretl_matrix_alloc(n, 1);
    b = gretl_matrix_alloc(n, 1);
    if (R == NULL || y == NULL || g == NULL || b == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* allocate storage in model struct */
    pmod->sderr = malloc(n * sizeof *pmod->sderr);
    pmod->yhat = malloc(fulln * sizeof *pmod->yhat);
    pmod->uhat = malloc(fulln * sizeof *pmod->uhat);
    if (pmod->sderr == NULL || pmod->yhat == NULL || 
	pmod->uhat == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* invert R */
    dtrtri_(&uplo, &diag, &n, Q->val, &lda, &info);
    if (info != 0) {
	fprintf(stderr, "dtrtri: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* copy the upper triangular R-inverse out of Q */
    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    if (i <= j) {
		gretl_matrix_set(R, i, j, 
				 gretl_matrix_get(Q, i, j));
	    } else {
		gretl_matrix_set(R, i, j, 0.0);
	    }
	}
    }

    /* obtain the real "Q" matrix */
    dorgqr_(&m, &n, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dorgqr: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* make "g" into gamma-hat */    
    gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
			      y, GRETL_MOD_NONE, g);

    /* OLS coefficients */
    gretl_matrix_multiply(R, g, b);
    save_coefficients(pmod, b, n);

    /* write vector of fitted values into y */
    gretl_matrix_multiply(Q, g, y);    

    /* get vector of residuals and SSR */
    get_resids_and_SSR(pmod, Z, y, ypy, fulln);

    /* standard error of regression */
    if (m - n > 0) {
	pmod->sigma = sqrt(pmod->ess / (m - n));
    } else {
	pmod->sigma = 0.0;
    }

    /* create (X'X)' */
    xpxinv = gretl_matrix_alloc(n, n);
    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      xpxinv);

    /* VCV and standard errors */
    if (opts & OPT_R) { 
	gretl_model_set_int(pmod, "robust", 1);
	if (opts & OPT_T) {
	    qr_make_hac(pmod, Z, xpxinv);
	} else {
	    qr_make_hccme(pmod, Z, xpxinv);
	}
    } else {
	qr_make_regular_vcv(pmod, xpxinv);
    }

    /* get R^2 and F-stat */
    qr_compute_r_squared(pmod, Z[pmod->list[1]], m);

 qr_cleanup:
    gretl_matrix_free(Q);
    gretl_matrix_free(y);

    free(tau); free(work);
    free(iwork);

    gretl_matrix_free(R);
    gretl_matrix_free(g);
    gretl_matrix_free(b);
    gretl_matrix_free(xpxinv);

    pmod->errcode = err;

    return err;    
}
