/* lad.c -- Least Absolute Deviation regression for gretl */

#include "libgretl.h"

static double toler = 1.0e-9;

static int l1_ (int m, int n, 
		double *a, double *b, 
		double *x, double *e);

static int col_(double *v1, double *v2, double amlt, 
		int m1, int iout);

static int 
bootstrap_stderrs (MODEL *pmod, double **Z,
		   double *a, double *b, double *e, double *x,
		   int m, int n, int dim);


int lad_driver (MODEL *pmod, double **Z, DATAINFO *pdinfo)
{
    double *a = NULL, *b = NULL, *e = NULL, *x = NULL;
    int i, j, k, t, m, n, nrows, dim;
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
    k = 0;
    for (j = 1; j <= n; j++) {
	int lj1 = pmod->list[j+1];

	t = pmod->t1;
	for (i = 0; i < m; i++) {
	    while (model_missing(pmod, t)) {
		t++;
	    }
	    a[i + k * nrows] = Z[lj1][t++];
	}
	k++;
    }

    t = pmod->t1;
    for (i = 0; i < m; i++) {
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
    } else {
	gretl_model_set_int(pmod, "ladcode", ladcode);
    }

    if (pmod->errcode == 0) {

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

	/* sum of absolute residuals (abuse of "rho") */
	pmod->rho = a[m + n * nrows];

	/* set ess-based stats to missing value */
	pmod->rsq = NADBL;
	pmod->adjrsq = NADBL;
	pmod->fstt = NADBL;

	/* LaPlace errors: equivalent of standard error is sum of
	   absolute residuals over nobs */
	pmod->sigma = pmod->rho / pmod->nobs; 

	if (bootstrap_stderrs (pmod, Z, a, b, e, x, m, n, dim)) {
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
    static double big = 1e200;

    /* System generated locals */
    int i1, i2;

    /* Local variables */
    double amin, amax;
    int test = 0, stage;
    double zero, one, two, sum;
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

    zero = 0.;
    one = 1.;
    two = 2.;
    m1 = m + 1;
    n1 = n + 1;
    m2 = m + 2;
    n2 = n + 2;
    i1 = n;
    for (j = 1; j <= i1; ++j) {
	a[m2 + j * nrows] = (double) j;
	x[j] = zero;
/* L10: */
    }
    i1 = m;
    for (i = 1; i <= i1; ++i) {
	a[i + n2 * nrows] = (double) (n + i);
	a[i + n1 * nrows] = b[i];
	if (b[i] >= zero) {
	    goto L30;
	}
	i2 = n2;
	for (j = 1; j <= i2; ++j) {
	    a[i + j * nrows] = -a[i + j * nrows];
/* L20: */
	}
L30:
	e[i] = zero;
/* L40: */
    }

    i1 = n1;
    for (j = 1; j <= i1; ++j) {
	sum = zero;
	i2 = m;
	for (i = 1; i <= i2; ++i) {
	    sum += a[i + j * nrows];
/* L50: */
	}
	a[m1 + j * nrows] = sum;
/* L60: */
    }

    stage = 1;
    kount = 0;
    kr = 1;
    kl = 1;
L70:
    amax = -one;
    i1 = n;
    for (j = kr; j <= i1; ++j) {
	if (fabs(a[m2 + j * nrows]) > (double) n) {
	    goto L80;
	}
	d = fabs(a[m1 + j * nrows]);
	if (d <= amax) {
	    goto L80;
	}
	amax = d;
	in = j;
L80:
	;
    }
    if (a[m1 + in * nrows] >= zero) {
	goto L100;
    }
    i1 = m2;
    for (i = 1; i <= i1; ++i) {
	a[i + in * nrows] = -a[i + in * nrows];
/* L90: */
    }

L100:
    k = 0;
    i1 = m;
    for (i = kl; i <= i1; ++i) {
	d = a[i + in * nrows];
	if (d <= toler) {
	    goto L110;
	}
	++k;
	b[k] = a[i + n1 * nrows] / d;
	ls[k - 1] = i;
	test = 1;
L110:
	;
    }
L120:
    if (k > 0) {
	goto L130;
    }
    test = 0;
    goto L150;
L130:
    amin = big;
    i1 = k;
    for (i = 1; i <= i1; ++i) {
	if (b[i] >= amin) {
	    goto L140;
	}
	j = i;
	amin = b[i];
	iout = ls[i - 1];
L140:
	;
    }
    b[j] = b[k];
    ls[j - 1] = ls[k - 1];
    --k;

L150:
    if (test || ! stage) {
	goto L170;
    }
    i1 = m2;
    for (i = 1; i <= i1; ++i) {
	d = a[i + kr * nrows];
	a[i + kr * nrows] = a[i + in * nrows];
	a[i + in * nrows] = d;
/* L160: */
    }
    ++kr;
    goto L260;
L170:
    if (test) {
	goto L180;
    }
    a[m2 + n1 * nrows] = two;
    goto L350;
L180:
    pivot = a[iout + in * nrows];
    if (a[m1 + in * nrows] - pivot - pivot <= toler) {
	goto L200;
    }
    i1 = n1;
    for (j = kr; j <= i1; ++j) {
	d = a[iout + j * nrows];
	a[m1 + j * nrows] = a[m1 + j * nrows] - d - d;
	a[iout + j * nrows] = -d;
/* L190: */
    }
    a[iout + n2 * nrows] = -a[iout + n2 * nrows];
    goto L120;

L200:
    i1 = n1;
    for (j = kr; j <= i1; ++j) {
	if (j == in) {
	    goto L210;
	}
	a[iout + j * nrows] /= pivot;
L210:
	;
    }
    i1 = n1;
    for (j = kr; j <= i1; ++j) {
	if (j == in) {
	    goto L220;
	}
	col_(&a[j * nrows + 1], &a[in * nrows + 1], a[iout + j * nrows], m1, iout);
L220:
	;
    }
    i1 = m1;
    for (i = 1; i <= i1; ++i) {
	if (i == iout) {
	    goto L240;
	}
	a[i + in * nrows] = -a[i + in * nrows] / pivot;
L240:
	;
    }
    a[iout + in * nrows] = one / pivot;
    d = a[iout + n2 * nrows];
    a[iout + n2 * nrows] = a[m2 + in * nrows];
    a[m2 + in * nrows] = d;
    ++kount;
    if (! stage) {
	goto L270;
    }

    ++kl;
    i1 = n2;
    for (j = kr; j <= i1; ++j) {
	d = a[iout + j * nrows];
	a[iout + j * nrows] = a[kount + j * nrows];
	a[kount + j * nrows] = d;
/* L250: */
    }
L260:
    if (kount + kr != n1) {
	goto L70;
    }

    stage = 0;

L270:
    amax = -big;
    i1 = n;
    for (j = kr; j <= i1; ++j) {
	d = a[m1 + j * nrows];
	if (d >= zero) {
	    goto L280;
	}
	if (d > -two) {
	    goto L290;
	}
	d = -d - two;
L280:
	if (d <= amax) {
	    goto L290;
	}
	amax = d;
	in = j;
L290:
	;
    }
    if (amax <= toler) {
	goto L310;
    }
    if (a[m1 + in * nrows] > zero) {
	goto L100;
    }
    i1 = m2;
    for (i = 1; i <= i1; ++i) {
	a[i + in * nrows] = -a[i + in * nrows];
/* L300: */
    }
    a[m1 + in * nrows] -= two;
    goto L100;

L310:
    l = kl - 1;
    i1 = l;
    for (i = 1; i <= i1; ++i) {
	if (a[i + n1 * nrows] >= zero) {
	    goto L330;
	}
	i2 = n2;
	for (j = kr; j <= i2; ++j) {
	    a[i + j * nrows] = -a[i + j * nrows];
/* L320: */
	}
L330:
	;
    }
    a[m2 + n1 * nrows] = zero;
    if (kr != 1) {
	goto L350;
    }
    i1 = n;
    for (j = 1; j <= i1; ++j) {
	d = fabs(a[m1 + j * nrows]);
	if (d <= toler || two - d <= toler) {
	    goto L350;
	}
/* L340: */
    }
    a[m2 + n1 * nrows] = one;
L350:
    i1 = m;
    for (i = 1; i <= i1; ++i) {
	k = (int) a[i + n2 * nrows];
	d = a[i + n1 * nrows];
	if (k > 0) {
	    goto L360;
	}
	k = -k;
	d = -d;
L360:
	if (i >= kl) {
	    goto L370;
	}
	x[k] = d;
	goto L380;
L370:
	k -= n;
	e[k] = d;
L380:
	;
    }
    a[m2 + n2 * nrows] = (double) kount;
    a[m1 + n2 * nrows] = (double) (n1 - kr);
    sum = zero;
    i1 = m;
    for (i = kl; i <= i1; ++i) {
	sum += a[i + n1 * nrows];
/* L390: */
    }
    a[m1 + n1 * nrows] = sum;

    free(ls);

    return 0;
} 

static int col_(double *v1, double *v2, double amlt, 
		int m1, int iout)
{
    int i;

    /* Parameter adjustments */
    --v2;
    --v1;

    for (i = 1; i <= m1; ++i) {
	if (i != iout) {
	    v1[i] -= v2[i] * amlt;
	}
    }
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

/* obtain bootstrap estimates of LAD standard errors */

static int 
bootstrap_stderrs (MODEL *pmod, double **Z,
		   double *a, double *b, double *e, double *x,
		   int m, int n, int dim)
{
    double **coeffs = NULL;
    int i, j, k, r, *sample = NULL;
    int yno = pmod->list[1];
    int nrows = m + 2;

    /* an array for each coefficient */
    coeffs = malloc(pmod->ncoeff * sizeof *coeffs);
    if (coeffs == NULL) {
	return 1;
    }

    /* each array has length ITERS + 1 */
    for (i=0; i<pmod->ncoeff; i++) {
	coeffs[i] = malloc((ITERS + 1) * sizeof **coeffs);
	if (coeffs[i] == NULL) {
	    for (k=0; k<i; k++) {
		free(coeffs[k]);
	    }
	    free(coeffs);
	    return 1;
	}
    }

    /* sample array has length pmod->nobs */
    sample = malloc(m * sizeof *sample);
    if (sample == NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    free(coeffs[i]);
	}
	free(coeffs);
	return 1;
    }

    for (k=0; k<ITERS; k++) {

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
	
	/* create random sample index array */
	for (i=0; i<m; i++) {
	    sample[i] = pmod->t1 + gretl_rand_int_max(m);
	}

	if (pmod->missmask != NULL) {
	    adjust_sample_for_missing(sample, m, pmod);
	}

	/* populate data array using the sample */
	r = 0;
	for (j = 1; j <= n; j++) {
	    for (i = 0; i < m; i++) {
		a[i + r * nrows] = Z[pmod->list[j+1]][sample[i]];
	    }
	    r++;
	}

	for (i = 0; i < m; i++) {
	    b[i] = a[i + n * nrows] = Z[yno][sample[i]];
	}

	/* estimate LAD model and store coeffs */
	l1_(m, n, a, b, x, e);

	for (i=0; i<n; i++) {
	    coeffs[i][k] = x[i];
	}
    }

    /* initialize means and standard deviations */
    for (i=0; i<pmod->ncoeff; i++) {
	coeffs[i][ITERS] = 0.0;
	pmod->sderr[i] = 0.0;
    }

    /* find means of coeff estimates */
    for (i=0; i<pmod->ncoeff; i++) {
	for (k=0; k<ITERS; k++) {
	   coeffs[i][ITERS] += coeffs[i][k];
	} 
	coeffs[i][ITERS] /= ITERS;
    }    

    /* find standard deviations */
    for (i=0; i<pmod->ncoeff; i++) {
	for (k=0; k<ITERS; k++) {
	   pmod->sderr[i] += (coeffs[i][k] - coeffs[i][ITERS]) *
	       (coeffs[i][k] - coeffs[i][ITERS]);
	}
	pmod->sderr[i] = sqrt(pmod->sderr[i] / ITERS);
    }

    free(sample);

    for (i=0; i<pmod->ncoeff; i++) {
	free(coeffs[i]);
    }

    free(coeffs);

    return 0;
}

