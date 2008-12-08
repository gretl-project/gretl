/* lad.c -- Least Absolute Deviation regression for gretl */

#include "libgretl.h"

static double toler = 1.0e-9;

static void col_ (double *v1, double *v2, double amlt, 
		  int m1, int iout)
{
    int i;

    /* Parameter adjustments */
    --v2;
    --v1;

    for (i=1; i<=m1; i++) {
	if (i != iout) {
	    v1[i] -= v2[i] * amlt;
	}
    }
}

/*    Based on SUBROUTINE L1(M,N,TOLER,X,A,B)
C
C  ***************************************************************
C  *  COMMUNICATIONS OF THE ASSOCIATION FOR COMPUTING MACHINERY  *
C  *  ALGORITHM 478                                              *
C  *  SOLUTION OF AN OVERDETERMINED SYSTEM OF EQUATIONS IN THE   *
C  *  L1 NORM -- I. BARRODALE AND F.D.K. ROBERTS                 *
C  ***************************************************************
C
C  ****************************************************************
C  * THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX METHOD OF *
C  * LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION TO THE OVER-  *
C  * DETERMINED SYSTEM OF LINEAR EQUATIONS, A*X = B.
C  *                                                              *
C  *                         PARAMETERS                           *
C  *                                                              *
C  * M     - NUMBER OF EQUATIONS.  NOTE THAT THIS IS THE N IN     *
C  *         SUBROUTINE LAD.                                      *
C  * N     - NUMBER OF UNKNOWNS (M.GE.N).                         *
C  * TOLER - A SMALL POSITIVE TOLERANCE. THE ROUTINE REGARDS ANY  *
C  *         QUANTITY AS ZERO UNLES ITS MAGNITUDE EXCEEDS TOLER.  *
C  *         EMPIRICAL EVIDENCE SUGGESTS TOLER = 10**(-D*2/3)     *
C  *         WHERE D REPRESENTS THE NUMBER OF DECIMAL DIGITS OF   *
C  *         ACCURACY AVAILABLE.                                  *
C  * X     - ONE DIMENSIONAL REAL ARRAY.  ON EXIT, THIS ARRAY     *
C  *         CONTAINS A SOLUTION TO THE L1 PROBLEM.               *
C  *                                                              *
C  * THE ORIGINAL CONTENTS OF THE ARRAYS A AND B ARE DESTROYED BY *
C  * THIS ROUTINE.                                                *
C  *                                                              *
C  * ON EXIT FROM THE SUBROUTINE, THE ARRAY A CONTAINS THE        *
C  * FOLLOWING INFORMATION.                                       *
C  *                                                              *
C  * A(M+1,N+1)  THE MINIMUM SUM OF THE ABSOLUTE VALUES OF THE    *
C  *             RESIDUALS.                                       *
C  * A(M+1,N+2)  THE RANK OF THE MATRIX OF COEFFICIENTS.          *
C  * A(M+2,N+1)  EXIT CODE WITH VALUES                            *
C  *             0 - OPTIMAL SOLUTION - PROBABLY NON-UNIQUE       *
C  *             1 - UNIQUE OPTIMAL SOLUTION                      *
C  *             2 - PREMATURE TERMINATION DUE TO ROUNDING ERRORS *
C  * A(M+2,N+2)  NUMBER OF SIMPLEX ITERATIONS PERFORMED.          *
C  *                                                              *
C  ****************************************************************
C */

/* Modifications to the above: "toler" is not entered as a parameter,
   it is "global" to this translation unit.  Regression residuals
   are found in the array e on exit.  The following function was
   translated from the Fortran by f2c, then rendered into slightly more
   idiomatic C by me.  Allin Cottrell, September 2002.
*/

static int l1_ (int m, int n, 
		double *a, double *b, 
		double *x, double *e)
{
    const double big = 1e200;
    double amin, amax;
    int test = 0, stage;
    double sum;
    int iout = 0;
    int i, j, k, l;
    int m1, n1, m2, n2, kount;
    double d, pivot;
    int kl, in = 0, kr, *ls;
    int nrows = m + 2;

    ls = malloc(m * sizeof *ls);

    /* Parameter adjustments */
    --e;
    --x;
    --b;
    a -= (nrows + 1);

    m1 = m + 1;
    n1 = n + 1;
    m2 = m + 2;
    n2 = n + 2;

    for (j = 1; j <= n; ++j) {
	a[m2 + j * nrows] = (double) j;
	x[j] = 0.0;
    }

    for (i = 1; i <= m; ++i) {
	a[i + n2 * nrows] = (double) (n + i);
	a[i + n1 * nrows] = b[i];
	if (b[i] < 0.0) {
	    for (j = 1; j <= n2; ++j) {
		a[i + j * nrows] = -a[i + j * nrows];
	    }
	}
	e[i] = 0.0;
    }

    for (j = 1; j <= n1; ++j) {
	sum = 0.0;
	for (i = 1; i <= m; ++i) {
	    sum += a[i + j * nrows];
	}
	a[m1 + j * nrows] = sum;
    }

    stage = 1;
    kount = 0;
    kr = 1;
    kl = 1;

L70:
    amax = -1.0;
    for (j = kr; j <= n; ++j) {
	if (fabs(a[m2 + j * nrows]) > (double) n) {
	    continue;
	}
	d = fabs(a[m1 + j * nrows]);
	if (d <= amax) {
	    continue;
	}
	amax = d;
	in = j;
    }

    if (a[m1 + in * nrows] < 0.0) {
	for (i = 1; i <= m2; ++i) {
	    a[i + in * nrows] = -a[i + in * nrows];
	}
    }

L100:
    k = 0;
    for (i = kl; i <= m; ++i) {
	d = a[i + in * nrows];
	if (d > toler) {
	    ++k;
	    b[k] = a[i + n1 * nrows] / d;
	    ls[k - 1] = i;
	    test = 1;
	}
    }

L120:
    if (k <= 0) {
	test = 0;
	goto L150;
    }

    amin = big;
    for (i = 1; i <= k; ++i) {
	if (b[i] < amin) {
	    j = i;
	    amin = b[i];
	    iout = ls[i - 1];
	}
    }
    b[j] = b[k];
    ls[j - 1] = ls[k - 1];
    --k;

L150:
    if (test || !stage) {
	goto L170;
    }
    for (i = 1; i <= m2; ++i) {
	d = a[i + kr * nrows];
	a[i + kr * nrows] = a[i + in * nrows];
	a[i + in * nrows] = d;
    }
    ++kr;
    goto L260;

L170:
    if (test) {
	goto L180;
    }
    a[m2 + n1 * nrows] = 2.0;
    goto L350;

L180:
    pivot = a[iout + in * nrows];
    if (a[m1 + in * nrows] - pivot - pivot <= toler) {
	goto L200;
    }
    for (j = kr; j <= n1; ++j) {
	d = a[iout + j * nrows];
	a[m1 + j * nrows] = a[m1 + j * nrows] - d - d;
	a[iout + j * nrows] = -d;
    }
    a[iout + n2 * nrows] = -a[iout + n2 * nrows];
    goto L120;

L200:
    for (j = kr; j <= n1; ++j) {
	if (j != in) {
	    a[iout + j * nrows] /= pivot;
	}
    }
    for (j = kr; j <= n1; ++j) {
	if (j != in) {
	    col_(&a[j * nrows + 1], &a[in * nrows + 1], a[iout + j * nrows], 
		 m1, iout);
	}
    }
    for (i = 1; i <= m1; ++i) {
	if (i != iout) {
	    a[i + in * nrows] = -a[i + in * nrows] / pivot;
	}
    }

    a[iout + in * nrows] = 1.0 / pivot;
    d = a[iout + n2 * nrows];
    a[iout + n2 * nrows] = a[m2 + in * nrows];
    a[m2 + in * nrows] = d;
    ++kount;

    if (!stage) {
	goto L270;
    }

    ++kl;
    for (j = kr; j <= n2; ++j) {
	d = a[iout + j * nrows];
	a[iout + j * nrows] = a[kount + j * nrows];
	a[kount + j * nrows] = d;
    }

L260:
    if (kount + kr != n1) {
	goto L70;
    }

    stage = 0;

L270:
    amax = -big;
    for (j = kr; j <= n; ++j) {
	d = a[m1 + j * nrows];
	if (d < 0.0) {
	    if (d > -2.0) {
		continue;
	    }
	    d = -d - 2.0;
	}
	if (d <= amax) {
	    continue;
	}
	amax = d;
	in = j;
    }

    if (amax > toler) {
	if (a[m1 + in * nrows] > 0.0) {
	    goto L100;
	}
	for (i = 1; i <= m2; ++i) {
	    a[i + in * nrows] = -a[i + in * nrows];
	}
	a[m1 + in * nrows] -= 2.0;
	goto L100;
    }

    l = kl - 1;
    for (i = 1; i <= l; ++i) {
	if (a[i + n1 * nrows] < 0.0) {
	    for (j = kr; j <= n2; ++j) {
		a[i + j * nrows] = -a[i + j * nrows];
	    }
	}
    }

    a[m2 + n1 * nrows] = 0.0;
    if (kr != 1) {
	goto L350;
    }

    for (j = 1; j <= n; ++j) {
	d = fabs(a[m1 + j * nrows]);
	if (d <= toler || 2.0 - d <= toler) {
	    goto L350;
	}
    }
    a[m2 + n1 * nrows] = 1.0;

L350:
    for (i = 1; i <= m; ++i) {
	k = (int) a[i + n2 * nrows];
	d = a[i + n1 * nrows];
	if (k <= 0) {
	    k = -k;
	    d = -d;
	}
	if (i >= kl) {
	    k -= n;
	    e[k] = d;
	} else {
	    x[k] = d;
	}
    }

    a[m2 + n2 * nrows] = (double) kount;
    a[m1 + n2 * nrows] = (double) (n1 - kr);

    sum = 0.0;
    for (i = kl; i <= m; ++i) {
	sum += a[i + n1 * nrows];
    }

    a[m1 + n1 * nrows] = sum;

    free(ls);

    return 0;
} 

static int missobs_before (const MODEL *pmod, int t)
{
    int i, c = 0;

    for (i=0; i<t; i++) {
	if (pmod->missmask[i + pmod->t1] == '1') {
	    c++;
	}
    }

    return c;
}

static void 
adjust_sample_for_missing (int *sample, int n, const MODEL *pmod)
{
    int i;

    for (i=0; i<n; i++) {
	sample[i] += missobs_before(pmod, sample[i]);
    }
}

#define ITERS 500
#define RESAMPLE_RESIDUALS 0

#if RESAMPLE_RESIDUALS

/* populate dependent var using resampled residuals */

static void
make_data_arrays (MODEL *pmod, double **Z,
		  double *a, double *b,
		  const int *sample,
		  int nrows, int k, int m)
{
    int i, j, v, t;

    /* we need to do this on each iteration because the "a" array is
       overwritten by the LAD calculations
    */
    for (j=0; j<k; j++) {
	v = pmod->list[j+2];
	t = pmod->t1;
	for (i=0; i<m; i++) {
	    while (model_missing(pmod, t)) {
		t++;
	    }
	    a[i + j * nrows] = Z[v][t++];
	}
    }

    t = pmod->t1;

    for (i=0; i<m; i++) {
	while (model_missing(pmod, t)) {
	    t++;
	}
	b[i] = a[i + k * nrows] = 
	    pmod->yhat[t++] + pmod->uhat[sample[i]];
    }
}

#else

/* populate both y and X using resampled data rows */

static void
make_data_arrays (MODEL *pmod, double **Z,
		  double *a, double *b,
		  const int *sample,
		  int nrows, int k, int m)
{
    int i, j, v, t;

    for (i=0; i<m; i++) {
	t = sample[i];
	for (j=0; j<k; j++) {
	    v = pmod->list[j+2];
	    a[i + j * nrows] = Z[v][t];
	}
	v = pmod->list[1];
	b[i] = a[i + k * nrows] = Z[v][t];
    }
}

#endif

/* obtain bootstrap estimates of LAD covariance matrix */

static int bootstrap_vcv (MODEL *pmod, double **Z,
			  double *a, double *b, double *e, double *x,
			  int m, int n, int dim)
{
    double **coeffs = NULL;
    double *meanb = NULL;
    int *sample = NULL;
    double xi, xj;
    int i, j, k;
    int nvcv, nrows = m + 2;
    int err = 0;

    /* note: new_vcv sets all entries to zero */
    err = gretl_model_new_vcv(pmod, &nvcv);
    if (err) {
	return err;
    }

    /* an array for each coefficient */
    coeffs = malloc(pmod->ncoeff * sizeof *coeffs);
    if (coeffs == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* a scalar for each coefficient mean */
    meanb = malloc(pmod->ncoeff * sizeof *meanb);
    if (meanb == NULL) {
	err = E_ALLOC;
	goto bailout;
    }    

    /* each array has length ITERS */
    for (i=0; i<pmod->ncoeff; i++) {
	coeffs[i] = malloc(ITERS * sizeof **coeffs);
	if (coeffs[i] == NULL) {
	    for (k=0; k<i; k++) {
		free(coeffs[k]);
	    }
	    free(coeffs);
	    coeffs = NULL;
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    /* sample array has length pmod->nobs */
    sample = malloc(m * sizeof *sample);
    if (sample == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (k=0; k<ITERS; k++) {

	/* create random sample index array */
	for (i=0; i<m; i++) {
	    sample[i] = pmod->t1 + gretl_rand_int_max(m);
	}

	if (pmod->missmask != NULL) {
	    adjust_sample_for_missing(sample, m, pmod);
	}

	/* initialize arrays */
	for (i=0; i<dim; i++) {
	    a[i] = 0.0;
	}
	for (i=0; i<m; i++) {
	    e[i] = b[i] = 0.0;
	}
	for (i=0; i<n; i++) {
	    x[i] = 0.0;
	}

	make_data_arrays(pmod, Z, a, b, sample, nrows, n, m);

	/* estimate LAD model and store coeffs */
	l1_(m, n, a, b, x, e);

	for (i=0; i<n; i++) {
	    coeffs[i][k] = x[i];
	}
    }

    /* find means of coeff estimates */
    for (i=0; i<pmod->ncoeff; i++) {
	double bbar = 0.0;

	for (k=0; k<ITERS; k++) {
	   bbar += coeffs[i][k];
	} 
	meanb[i] = bbar / ITERS;
    }    

    /* find variances and covariances */
    for (i=0; i<pmod->ncoeff; i++) {
	double vi = 0.0;

	for (k=0; k<ITERS; k++) {
	    xi = coeffs[i][k] - meanb[i];
	    vi += xi * xi;
	    for (j=0; j<=i; j++) {
		xj = coeffs[j][k] - meanb[j];
		pmod->vcv[ijton(i, j, pmod->ncoeff)] += xi * xj;
	    }
	}
	pmod->sderr[i] = sqrt(vi / ITERS);
    }

    for (i=0; i<nvcv; i++) {
	pmod->vcv[i] /= ITERS;
    }

 bailout:

    free(sample);

    if (coeffs != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    free(coeffs[i]);
	}
	free(coeffs);
    }

    free(meanb);

    return 0;
}

static void lad_loglik (MODEL *pmod)
{
    double tau = 0.5, R = 0.0;
    int n = pmod->nobs;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	R += pmod->uhat[t] * (tau - (pmod->uhat[t] < 0));
    }

    pmod->lnL = n * (log(tau * (1-tau)) - 1 - log(R/n));
    mle_criteria(pmod, 0);
}

int lad_driver (MODEL *pmod, double **Z, DATAINFO *pdinfo)
{
    double *a = NULL, *b = NULL, *e = NULL, *x = NULL;
    int i, j, t, m, n, nrows, dim;
    int yno = pmod->list[1];
    int ladcode;

    m = pmod->nobs;
    n = pmod->list[0] - 1;

    nrows = m + 2;
    dim = nrows * (n + 2);

    a = malloc(dim * sizeof *a);
    x = malloc(n * sizeof *x);
    e = malloc(m * sizeof *e);
    b = malloc(m * sizeof *b);

    if (a == NULL || x == NULL || e == NULL || b == NULL) {
	free(a);
	free(x);
	free(e);
	free(b);
	return 1;
    }

    /* initialize arrays */
    for (i=0; i<dim; i++) {
	a[i] = 0.0;
    }
    for (i=0; i<m; i++) {
	e[i] = b[i] = 0.0;
    }
    for (i=0; i<n; i++) {
	x[i] = 0.0;
    }

    /* populate data array */
    for (j=0; j<n; j++) {
	int v = pmod->list[j+2];

	t = pmod->t1;
	for (i=0; i<m; i++) {
	    while (model_missing(pmod, t)) {
		t++;
	    }
	    a[i + j * nrows] = Z[v][t++];
	}
    }

    t = pmod->t1;
    for (i=0; i<m; i++) {
	while (model_missing(pmod, t)) {
	    t++;
	}
	b[i] = a[i + n * nrows] = Z[yno][t++];
    }

    l1_(m, n, a, b, x, e);

    /* handle case where exit code indicates numeric error */
    ladcode = (int) a[m + 1 + n * nrows];
    if (ladcode == 2) {
	pmod->errcode = E_SINGULAR;
    } else if (ladcode == 0) {
	gretl_model_set_int(pmod, "nonunique", 1);
    }

    if (pmod->errcode == 0) {
	double SAR = a[m + n * nrows];

	for (i=0; i<n; i++) {
	    pmod->coeff[i] = x[i];
	}

	pmod->ess = 0.0;
	for (i=0; i<m; i++) {
	    t = i + pmod->t1;
	    pmod->yhat[t] = Z[yno][t] - e[i];
	    pmod->uhat[t] = e[i];
	    pmod->ess += e[i] * e[i];
	}

	/* sum of absolute residuals */
	gretl_model_set_double(pmod, "ladsum", SAR);

	/* median of dependent variable */
	gretl_model_add_y_median(pmod, Z[yno]);

	/* set various stats to missing value */
	pmod->rsq = pmod->adjrsq = NADBL;
	pmod->fstt = pmod->chisq = NADBL;

	/* LaPlace errors: equivalent of standard error is sum of
	   absolute residuals over nobs */
	pmod->sigma = SAR / pmod->nobs; 

	lad_loglik(pmod);

	if (bootstrap_vcv(pmod, Z, a, b, e, x, m, n, dim)) {
	    pmod->errcode = E_ALLOC;
	}
    }

    pmod->ci = LAD;

    free(a);
    free(x);
    free(e);
    free(b);

    return 0;
}
