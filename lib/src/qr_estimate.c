#include "libgretl.h"
#include "gretl_matrix.h"

/* In fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are
   contiguous.
   
*/

#define TINY 1e-9 /* was 1.0e-13, produced poor results on NIST Filip
		       test */

int gretl_qr_regress (MODEL *pmod, double **Z, int fulln)
{
    integer info, lwork;
    integer m, n, lda;
    integer *iwork;
    gretl_matrix *Q, *y;
    gretl_matrix *R = NULL, *g = NULL, *b = NULL;
    gretl_matrix *xpxinv = NULL;
    doublereal *tau, *work;
    doublereal rcond;
    char uplo = 'U';
    char diag = 'N';
    char norm = '1';
    int i, j, t;
    int err = 0;

    m = pmod->t2 - pmod->t1 + 1;  /* # of rows = # of observations */
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

    /* copy independent var values into Q */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    Q->val[j++] = Z[pmod->list[i]][t];
	}
    }

    /* fill out dependent variable vector */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	gretl_matrix_set(y, t, 0, Z[pmod->list[1]][t]);
    }

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

    /* check condition number of R */
    dtrcon_(&norm, &uplo, &diag, &n, Q->val, &lda, &rcond, work, 
	    iwork, &info);
    if (info != 0) {
	fprintf(stderr, "dtrcon: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    fprintf(stderr, "rcond = %g\n", rcond);
    if (rcond < TINY) {
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
    pmod->coeff = malloc((n + 1) * sizeof *pmod->coeff);
    pmod->sderr = malloc((n + 1) * sizeof *pmod->sderr);
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

    /* copy the upper triangular R out of Q */
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

    /* write vector of fitted values into y */
    gretl_matmult(Q, g, y);

    /* get vector of residuals and SSR */
    pmod->ess = 0.0;
    i = 0;
    for (t=0; t<fulln; t++) {
	if (t < pmod->t1 || t > pmod->t2) {
	    pmod->yhat[t] = pmod->uhat[t] == NADBL;
	} else {
	    pmod->yhat[t] = y->val[i];
	    pmod->uhat[t] = Z[pmod->list[1]][t] - y->val[i];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    i++;
	}
    }
    pmod->sigma = sqrt(pmod->ess / (m - n));

    printf("SSR = %g, SE = %g\n", pmod->ess, pmod->sigma);

    /* OLS coefficients */
    gretl_matmult(R, g, b);
    for (i=0; i<n; i++) {
	pmod->coeff[i+1] = b->val[i];
    }

    /* create (X'X)' */
    xpxinv = gretl_matrix_alloc(n, n);
    gretl_matmult_mod(R, GRETL_MOD_NONE,
		      R, GRETL_MOD_TRANSPOSE,
		      xpxinv);

    /* get standard errors */
    for (i=0; i<n; i++) {
	double x = gretl_matrix_get(xpxinv, i, i);

	pmod->sderr[i+1] = pmod->sigma * sqrt(x);
    }

 qr_cleanup:
    gretl_matrix_free(Q);
    gretl_matrix_free(y);

    free(tau); free(work);
    free(iwork);

    gretl_matrix_free(R);
    gretl_matrix_free(g);
    gretl_matrix_free(b);
    gretl_matrix_free(xpxinv);

    return err;    
}
