/* ladrho.c */

#include <stdlib.h>
#include "libgretl.h"

#define TRUE_ 1
#define FALSE_ 0

static double toler = 1.0e-9;

static int lad_calculate (int m, int n, 
			  double *a, double *b, 
			  double *x, double *e);

static int col_(double *v1, double *v2, double amlt, 
		int m1, int iout);


int lad_driver (MODEL *pmod, double **Z, DATAINFO *pdinfo, PRN *prn)
{
    int iter;
    double *a = NULL, *b = NULL, *e = NULL, *x = NULL;
    int i, j, k, m, n, nrows, dim;
    int irank;
    double sumre;
    int iexcode;

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
    for (j = 1; j <= n; ++j) {
	for (i = 1; i <= m; ++i) {
	    a[i + k * nrows - 1] = Z[pmod->list[j+1]][i - 1 + pmod->t1];
	}
	k++;
    }

    for (i = 1; i <= m; ++i) {
	b[i - 1] = a[i + (n + 1) * nrows - (nrows + 1)] =
	    Z[pmod->list[1]][i - 1 + pmod->t1];
    }

    lad_calculate(m, n, a, b, x, e);

    for (i=0; i<n; i++) {
	pmod->coeff[i+1] = x[i];
    }

    irank = (int) a[m + 1 + (n + 2) * nrows - (nrows + 1)];
    iter = (int) a[m + 2 + (n + 2) * nrows - (nrows + 1)];
    iexcode = (int) a[m + 2 + (n + 1) * nrows - (nrows + 1)];
    sumre = a[m + 1 + (n + 1) * nrows - (nrows + 1)];

    pprintf(prn, "LAD: irank=%d, iter=%d, iexcode=%d\n",
	    irank, iter, iexcode);

    pmod->ess = sumre;
    pmod->sigma = sumre / pmod->nobs;

    free(a);
    free(x);
    free(e);
    free(b);

    return 0;
}

static int lad_calculate (int m, int n, 
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

    stage = TRUE_;
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
	test = TRUE_;
L110:
	;
    }
L120:
    if (k > 0) {
	goto L130;
    }
    test = FALSE_;
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

    stage = FALSE_;

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
