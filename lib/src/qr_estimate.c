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
	    pmod->rsq = corrrsq(n, y + t1, pmod->yhat + t1);
	    pmod->adjrsq = 
		1 - ((1 - pmod->rsq) * (n - 1) / pmod->dfd);
	}
	pmod->fstt = (pmod->tss - pmod->ess) * pmod->dfd / 
	    (pmod->ess * pmod->dfn);
    } else {
	pmod->rsq = 1.0;
	pmod->fstt = NADBL;
    }
}

static int qr_make_vcv (MODEL *pmod, gretl_matrix *v)
{
    const int nv = pmod->ncoeff;
    const int nterms = nv * (nv + 1) / 2;
    double x;
    int i, j, k;

    pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) return 1;  

    for (i=0; i<nv; i++) {
	for (j=0; j<=i; j++) {
	    k = ijton(i+1, j+1, nv);
	    x = gretl_matrix_get(v, i, j);
	    x *= pmod->sigma * pmod->sigma;
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

		x *= Z[pmod->nwt][t];
		pmod->yhat[t] = y->val[i];
		pmod->uhat[t] = x - y->val[i];
		pmod->ess += pmod->uhat[t] * pmod->uhat[t];
		i++;
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

static double get_model_data (MODEL *pmod, const double **Z, 
			      gretl_matrix *Q, gretl_matrix *y)
{
    int i, j, t;
    double x, ypy = 0.0;
    int t1 = pmod->t1;

    if (pmod->rho) t1++;

    /* copy independent vars into matrix Q */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	for (t=t1; t<=pmod->t2; t++) {
	    x = Z[pmod->list[i]][t];
	    if (pmod->nwt) {
		x *= Z[pmod->nwt][t];
	    } else if (pmod->rho && pmod->list[i] != 0) {
		x -= pmod->rho * Z[pmod->list[i]][t-1];
	    }
	    Q->val[j++] = x;
	}
    }

    /* copy dependent variable into y vector */
    j = 0;
    for (t=t1; t<=pmod->t2; t++) {
	x = Z[pmod->list[1]][t];
	if (pmod->nwt) {
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

int gretl_qr_regress (MODEL *pmod, const double **Z, int fulln)
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
    pmod->coeff = malloc(n * sizeof *pmod->coeff);
    pmod->sderr = malloc(n * sizeof *pmod->sderr);
    pmod->yhat = malloc(fulln * sizeof *pmod->yhat);
    pmod->uhat = malloc(fulln * sizeof *pmod->uhat);
    if (pmod->coeff == NULL || pmod->sderr == NULL || 
	pmod->yhat == NULL || pmod->uhat == NULL) {
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
    gretl_matmult_mod(Q, GRETL_MOD_TRANSPOSE,
		      y, GRETL_MOD_NONE, g);

    /* OLS coefficients */
    gretl_matmult(R, g, b);
    for (i=0; i<n; i++) {
	pmod->coeff[i] = b->val[i];
    }

    /* write vector of fitted values into y */
    gretl_matmult(Q, g, y);    

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
    gretl_matmult_mod(R, GRETL_MOD_NONE,
		      R, GRETL_MOD_TRANSPOSE,
		      xpxinv);

    /* get standard errors */
    for (i=0; i<n; i++) {
	double x = gretl_matrix_get(xpxinv, i, i);

	pmod->sderr[i] = pmod->sigma * sqrt(x);
    }

    /* set up covar matrix (triangular) */
    qr_make_vcv(pmod, xpxinv);

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
