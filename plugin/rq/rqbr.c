/* rqbr.f -- translated by f2c (version 20061008) */

#include "gretl_f2c.h"
#include <stdio.h>

/* Table of constant values */

#define RMAX 1000

static doublereal c_b15 = 1.;

static double d_sign (doublereal *a, doublereal *b)
{
    double x;
    x = (*a >= 0 ? *a : - *a);
    return ( *b >= 0 ? x : -x);
}

/* Output from Public domain Ratfor, version 1.0 */
int rqbr_(integer *m, integer *nn, integer *m5, integer *n3, 
	integer *n4, doublereal *a, doublereal *b, doublereal *t, doublereal *
	toler, integer *ift, doublereal *x, doublereal *e, integer *s, 
	doublereal *wa, doublereal *wb, integer *nsol, integer *ndsol, 
	doublereal *sol, doublereal *dsol, integer *lsol, integer *h__, 
	doublereal *qn, doublereal *cutoff, doublereal *ci, doublereal *tnmat,
	 doublereal *big, logical *lci1)
{
    /* System generated locals */
    integer h_dim1, h_offset, sol_dim1, sol_offset, a_dim1, a_offset, wa_dim1,
	     wa_offset, dsol_dim1, dsol_offset, i1, i2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i, j, k, l, n;
    static doublereal a1, b1;
    static integer m1, m2, m3, m4, n1, n2;
    static doublereal t0, t1;
    static integer kd, kl, in, kr;
    static doublereal tn, dif, min__, max__, aux;
    static logical lup;
    static doublereal sum, tnt;
    static integer out;
    static logical lci2;
    static doublereal half;
    static logical iend;
    static doublereal told;
    static logical init, skip;
    static doublereal smax, tnew;
    static logical test;
    static integer idxcf;
    static logical stage;
    static integer kount;
    static doublereal pivot;

    integer rcount = 0;

    /* Parameter adjustments */
    --wb;
    --s;
    --e;
    --b;
    tnmat -= 5;
    ci -= 5;
    --qn;
    --x;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    wa_dim1 = *m5;
    wa_offset = 1 + wa_dim1;
    wa -= wa_offset;
    h_dim1 = *nn;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    sol_dim1 = *n3;
    sol_offset = 1 + sol_dim1;
    sol -= sol_offset;
    dsol_dim1 = *m;
    dsol_offset = 1 + dsol_dim1;
    dsol -= dsol_offset;

    /* Function Body */
    n = *nn;
    *ift = 0;
    wa[*m + 2 + (*nn + 1) * wa_dim1] = 1.;
    if (*m5 != *m + 5) {
	*ift = 3;
    }
    if (*n3 != n + 3) {
	*ift = 4;
    }
    if (*n4 != n + 4) {
	*ift = 5;
    }
    if ((doublereal) (*m) <= 0. || (doublereal) n <= 0.) {
	*ift = 6;
    }
    if ((doublereal) (*ift) <= 2.) {
	half = .5;
	iend = TRUE_;
	lci2 = FALSE_;
	lup = TRUE_;
	skip = FALSE_;
	idxcf = 0;
	tnew = 0.;
	tn = 0.;
	m1 = *m + 1;
	n1 = n + 1;
	n2 = n + 2;
	m2 = *m + 2;
	m3 = *m + 3;
	m4 = *m + 4;
	i1 = n;
	for (j = 1; j <= i1; ++j) {
	    x[j] = 0.;
	}
	i1 = *m;
	for (i = 1; i <= i1; ++i) {
	    e[i] = 0.;
	}
	if (*t < 0. || *t > 1.) {
	    t0 = 1. / ((real) (*m) * 2.);
	    t1 = 1. - t0;
	    *t = t0;
	    iend = FALSE_;
	    *lci1 = FALSE_;
	}
L23016:
	i1 = *m;
	for (i = 1; i <= i1; ++i) {
	    k = 1;
	    i2 = *nn;
	    for (j = 1; j <= i2; ++j) {
		if (k <= *nn) {
		    if (j == idxcf) {
			skip = TRUE_;
		    } else {
			skip = FALSE_;
		    }
		    if (! skip) {
			wa[i + k * wa_dim1] = a[i + j * a_dim1];
			++k;
		    }
		}
	    }
	    wa[i + *n4 * wa_dim1] = (doublereal) (n + i);
	    wa[i + n2 * wa_dim1] = b[i];
	    if (idxcf != 0) {
		wa[i + *n3 * wa_dim1] = tnew * a[i + idxcf * a_dim1];
	    } else {
		wa[i + *n3 * wa_dim1] = 0.;
	    }
	    wa[i + n1 * wa_dim1] = wa[i + n2 * wa_dim1] - wa[i + *n3 * 
		    wa_dim1];
	    if (wa[i + n1 * wa_dim1] < 0.) {
		i2 = *n4;
		for (j = 1; j <= i2; ++j) {
		    wa[i + j * wa_dim1] = -wa[i + j * wa_dim1];
		}
	    }
	}
	i1 = n;
	for (j = 1; j <= i1; ++j) {
	    wa[m4 + j * wa_dim1] = (doublereal) j;
	    wa[m2 + j * wa_dim1] = 0.;
	    wa[m3 + j * wa_dim1] = 0.;
	    i2 = *m;
	    for (i = 1; i <= i2; ++i) {
		aux = d_sign(&c_b15, &wa[m4 + j * wa_dim1]) * wa[i + j * 
			wa_dim1];
		wa[m2 + j * wa_dim1] += aux * (1. - d_sign(&c_b15, &wa[i + *
			n4 * wa_dim1]));
		wa[m3 + j * wa_dim1] += aux * d_sign(&c_b15, &wa[i + *n4 * 
			wa_dim1]);
	    }
	    wa[m3 + j * wa_dim1] *= 2.;
	}
	dif = 0.;
	init = FALSE_;
	if (! lci2) {
	    i1 = n;
	    for (k = 1; k <= i1; ++k) {
		wa[*m5 + k * wa_dim1] = 0.;
		i2 = *m;
		for (i = 1; i <= i2; ++i) {
		    wa[*m5 + k * wa_dim1] += a[i + k * a_dim1];
		}
		wa[*m5 + k * wa_dim1] /= (real) (*m);
	    }
	}
	*lsol = 1;
	kount = 0;
L23045:
	i1 = n;
	for (j = 1; j <= i1; ++j) {
	    wa[m1 + j * wa_dim1] = wa[m2 + j * wa_dim1] + wa[m3 + j * wa_dim1]
		     * *t;
	}
	if (! init) {
	    stage = TRUE_;
	    kr = 1;
	    kl = 1;
	    goto L30;
	}
L23052:
	stage = FALSE_;
L23055:
	max__ = -(*big);
	i1 = n;
	for (j = kr; j <= i1; ++j) {
	    d__ = wa[m1 + j * wa_dim1];
	    if (d__ < 0.) {
		if (d__ > -2.) {
		    goto L23058;
		}
		d__ = -d__ - 2.;
	    }
	    if (d__ > max__) {
		max__ = d__;
		in = j;
	    }
L23058:
	    ;
	}
	if (max__ <= *toler) {
	    if (++rcount > RMAX) {
		fprintf(stderr, "max=%g, toler=%g, rcount=%d: stopping\n", 
			max__, *toler, rcount);
		*ift = 2;
		return 0;
	    }
	    goto L23054;
	}
	if (wa[m1 + in * wa_dim1] <= 0.) {
	    i1 = m4;
	    for (i = 1; i <= i1; ++i) {
		wa[i + in * wa_dim1] = -wa[i + in * wa_dim1];
	    }
	    wa[m1 + in * wa_dim1] += -2.;
	    wa[m2 + in * wa_dim1] += -2.;
	}
L23072:
	k = 0;
	i1 = *m;
	for (i = kl; i <= i1; ++i) {
	    d__ = wa[i + in * wa_dim1];
	    if (d__ > *toler) {
		++k;
		wb[k] = wa[i + n1 * wa_dim1] / d__;
		s[k] = i;
		test = TRUE_;
	    }
	}
L23079:
	if (k <= 0) {
	    test = FALSE_;
	} else {
	    min__ = *big;
	    i1 = k;
	    for (i = 1; i <= i1; ++i) {
		if (wb[i] < min__) {
		    j = i;
		    min__ = wb[i];
		    out = s[i];
		}
	    }
	    wb[j] = wb[k];
	    s[j] = s[k];
	    --k;
	}
	if (! test && stage) {
	    goto L23081;
	}
	if (! test) {
	    goto L23047;
	}
	pivot = wa[out + in * wa_dim1];
	if (wa[m1 + in * wa_dim1] - pivot - pivot <= *toler) {
	    goto L10;
	}
	i1 = *n3;
	for (j = kr; j <= i1; ++j) {
	    d__ = wa[out + j * wa_dim1];
	    wa[m1 + j * wa_dim1] = wa[m1 + j * wa_dim1] - d__ - d__;
	    wa[m2 + j * wa_dim1] = wa[m2 + j * wa_dim1] - d__ - d__;
	    wa[out + j * wa_dim1] = -d__;
	}
	wa[out + *n4 * wa_dim1] = -wa[out + *n4 * wa_dim1];
	goto L23079;
L23081:
	i1 = m4;
	for (i = 1; i <= i1; ++i) {
	    d__ = wa[i + kr * wa_dim1];
	    wa[i + kr * wa_dim1] = wa[i + in * wa_dim1];
	    wa[i + in * wa_dim1] = d__;
	}
	++kr;
	goto L20;
L10:
	i1 = *n3;
	for (j = kr; j <= i1; ++j) {
	    if (j != in) {
		wa[out + j * wa_dim1] /= pivot;
	    }
	}
	i1 = m3;
	for (i = 1; i <= i1; ++i) {
	    if (i != out) {
		d__ = wa[i + in * wa_dim1];
		i2 = *n3;
		for (j = kr; j <= i2; ++j) {
		    if (j != in) {
			wa[i + j * wa_dim1] -= d__ * wa[out + j * wa_dim1];
		    }
		}
	    }
	}
	i1 = m3;
	for (i = 1; i <= i1; ++i) {
	    if (i != out) {
		wa[i + in * wa_dim1] = -wa[i + in * wa_dim1] / pivot;
	    }
	}
	wa[out + in * wa_dim1] = 1. / pivot;
	d__ = wa[out + *n4 * wa_dim1];
	wa[out + *n4 * wa_dim1] = wa[m4 + in * wa_dim1];
	wa[m4 + in * wa_dim1] = d__;
	++kount;
	if (! stage) {
	    goto L23074;
	}
	++kl;
	i1 = *n4;
	for (j = kr; j <= i1; ++j) {
	    d__ = wa[out + j * wa_dim1];
	    wa[out + j * wa_dim1] = wa[kount + j * wa_dim1];
	    wa[kount + j * wa_dim1] = d__;
	}
L20:
	if (kount + kr == n1) {
	    goto L23057;
	}
L30:
	max__ = -1.;
	i1 = n;
	for (j = kr; j <= i1; ++j) {
	    if ((d__1 = wa[m4 + j * wa_dim1], abs(d__1)) <= (doublereal) n) {
		d__ = (d__1 = wa[m1 + j * wa_dim1], abs(d__1));
		if (d__ > max__) {
		    max__ = d__;
		    in = j;
		}
	    }
	}
	if (wa[m1 + in * wa_dim1] < 0.) {
	    i1 = m4;
	    for (i = 1; i <= i1; ++i) {
		wa[i + in * wa_dim1] = -wa[i + in * wa_dim1];
	    }
	}
	goto L23072;
L23074:
	goto L23055;
L23057:
	goto L23052;
L23054:
	if (kr == 1) {
	    i1 = n;
	    for (j = 1; j <= i1; ++j) {
		d__ = (d__1 = wa[m1 + j * wa_dim1], abs(d__1));
		if (d__ <= *toler || 2. - d__ <= *toler) {
		    *ift = 1;
		    wa[m2 + (*nn + 1) * wa_dim1] = 0.;
		    goto L80;
		}
	    }
	}
L80:
	kount = 0;
	sum = 0.;
	if (! lci2) {
	    i1 = kl - 1;
	    for (i = 1; i <= i1; ++i) {
		k = (integer) (wa[i + *n4 * wa_dim1] * d_sign(&c_b15, &wa[
			i + *n4 * wa_dim1]));
		x[k] = wa[i + n1 * wa_dim1] * d_sign(&c_b15, &wa[i + *n4 *
			 wa_dim1]);
	    }
	}
	i1 = n;
	for (i = 1; i <= i1; ++i) {
	    kd = (integer) ((d__1 = wa[m4 + i * wa_dim1], abs(d__1)) - n);
	    dsol[kd + *lsol * dsol_dim1] = wa[m1 + i * wa_dim1] / 2. + 1.;
	    if (wa[m4 + i * wa_dim1] < 0.) {
		dsol[kd + *lsol * dsol_dim1] = 1. - dsol[kd + *lsol * 
			dsol_dim1];
	    }
	    if (! lci2) {
		sum += x[i] * wa[*m5 + i * wa_dim1];
		sol[i + 3 + *lsol * sol_dim1] = x[i];
		h__[i + *lsol * h_dim1] = kd;
	    }
	}
	i1 = *m;
	for (i = kl; i <= i1; ++i) {
	    kd = (integer) ((d__1 = wa[i + *n4 * wa_dim1], abs(d__1)) - n);
	    if (wa[i + *n4 * wa_dim1] < 0.) {
		dsol[kd + *lsol * dsol_dim1] = 0.;
	    }
	    if (wa[i + *n4 * wa_dim1] > 0.) {
		dsol[kd + *lsol * dsol_dim1] = 1.;
	    }
	}
	if (! lci2) {
	    sol[*lsol * sol_dim1 + 1] = smax;
	    sol[*lsol * sol_dim1 + 2] = sum;
	    sum = 0.;
	    i1 = *m;
	    for (j = kl; j <= i1; ++j) {
		d__ = wa[j + n1 * wa_dim1] * d_sign(&c_b15, &wa[j + *n4 * 
			wa_dim1]);
		sum += d__ * (smax + half * (d_sign(&c_b15, &d__) - 1.));
	    }
	    sol[*lsol * sol_dim1 + 3] = sum;
	    i1 = *m;
	    for (i = 1; i <= i1; ++i) {
		dsol[i + (*lsol + 1) * dsol_dim1] = dsol[i + *lsol * 
			dsol_dim1];
	    }
	}
	if (lci2) {
	    a1 = 0.;
	    i1 = *m;
	    for (i = 1; i <= i1; ++i) {
		a1 += a[i + idxcf * a_dim1] * (dsol[i + *lsol * dsol_dim1]
			 + *t - 1.);
	    }
	    tn = a1 / sqrt(qn[idxcf] * *t * (1. - *t));
	    if (abs(tn) < *cutoff) {
		if (lup) {
		    smax = *big;
		} else {
		    smax = -(*big);
		}
		i1 = kl - 1;
		for (i = 1; i <= i1; ++i) {
		    k = (integer) (wa[i + *n4 * wa_dim1] * d_sign(&c_b15, &
			    wa[i + *n4 * wa_dim1]));
		    sol[k + sol_dim1] = wa[i + n2 * wa_dim1] * d_sign(&
			    c_b15, &wa[i + *n4 * wa_dim1]);
		    sol[k + (sol_dim1 << 1)] = wa[i + *n3 * wa_dim1] * 
			    d_sign(&c_b15, &wa[i + *n4 * wa_dim1]) / tnew;
		}
		i1 = *m;
		for (i = kl; i <= i1; ++i) {
		    a1 = 0.;
		    b1 = 0.;
		    k = (integer) (wa[i + *n4 * wa_dim1] * d_sign(&c_b15, &
			    wa[i + *n4 * wa_dim1]) - n);
		    l = 1;
		    i2 = n;
		    for (j = 1; j <= i2; ++j) {
			if (j == idxcf) {
			    ++l;
			}
			a1 += a[k + l * a_dim1] * sol[j + sol_dim1];
			b1 += a[k + l * a_dim1] * sol[j + (sol_dim1 << 1)];
			++l;
		    }
		    tnt = (b[k] - a1) / (a[k + idxcf * a_dim1] - b1);
		    if (lup) {
			if (tnt > tnew) {
			    if (tnt < smax) {
				smax = tnt;
				out = i;
			    }
			}
		    } else {
			if (tnt < tnew) {
			    if (tnt > smax) {
				smax = tnt;
				out = i;
			    }
			}
		    }
		}
		if (lup) {
		    told = tnew;
		    tnew = smax + *toler;
		    ci[(idxcf << 2) + 3] = told - *toler;
		    tnmat[(idxcf << 2) + 3] = tn;
		    if (! (tnew < *big - *toler)) {
			ci[(idxcf << 2) + 3] = *big;
			ci[(idxcf << 2) + 4] = *big;
			tnmat[(idxcf << 2) + 3] = tn;
			tnmat[(idxcf << 2) + 4] = tn;
			lup = FALSE_;
			goto L70;
		    }
		} else {
		    told = tnew;
		    tnew = smax - *toler;
		    ci[(idxcf << 2) + 2] = told + *toler;
		    tnmat[(idxcf << 2) + 2] = tn;
		    if (! (tnew > -(*big) + *toler)) {
			ci[(idxcf << 2) + 2] = -(*big);
			ci[(idxcf << 2) + 1] = -(*big);
			tnmat[(idxcf << 2) + 2] = tn;
			tnmat[(idxcf << 2) + 1] = tn;
			lup = TRUE_;
			goto L60;
		    }
		}
		i1 = *m;
		for (i = 1; i <= i1; ++i) {
		    wa[i + *n3 * wa_dim1] = wa[i + *n3 * wa_dim1] / told *
			     tnew;
		    wa[i + n1 * wa_dim1] = wa[i + n2 * wa_dim1] - wa[i 
			    + *n3 * wa_dim1];
		}
		i1 = *n3;
		for (j = kr; j <= i1; ++j) {
		    d__ = wa[out + j * wa_dim1];
		    wa[m1 + j * wa_dim1] = wa[m1 + j * wa_dim1] - d__ - d__;
		    wa[m2 + j * wa_dim1] = wa[m2 + j * wa_dim1] - d__ - d__;
		    wa[out + j * wa_dim1] = -d__;
		}
		wa[out + *n4 * wa_dim1] = -wa[out + *n4 * wa_dim1];
		init = TRUE_;
	    } else {
		if (lup) {
		    ci[(idxcf << 2) + 4] = tnew - *toler;
		    tnmat[(idxcf << 2) + 4] = tn;
		    lup = FALSE_;
		    goto L70;
		} else {
		    ci[(idxcf << 2) + 1] = tnew + *toler;
		    tnmat[(idxcf << 2) + 1] = tn;
		    lup = TRUE_;
		    goto L60;
		}
	    }
	}
	if (iend && ! lci2) {
	    goto L40;
	}
	if (! lci2) {
	    init = TRUE_;
	    ++(*lsol);
	    i1 = *m;
	    for (i = 1; i <= i1; ++i) {
		s[i] = 0;
	    }
	    i1 = n;
	    for (j = 1; j <= i1; ++j) {
		x[j] = 0.;
	    }
	    smax = 2.;
	    i1 = n;
	    for (j = 1; j <= i1; ++j) {
		b1 = wa[m3 + j * wa_dim1];
		a1 = (-2. - wa[m2 + j * wa_dim1]) / b1;
		b1 = -wa[m2 + j * wa_dim1] / b1;
		if (a1 >= *t) {
		    if (a1 < smax) {
			smax = a1;
			dif = (b1 - a1) / 2.;
		    }
		}
		if (b1 > *t) {
		    if (b1 < smax) {
			smax = b1;
			dif = (b1 - a1) / 2.;
		    }
		}
	    }
	    tnt = smax + *toler * (abs(dif) + 1.);
	    if (tnt >= t1 + *toler) {
		iend = TRUE_;
	    }
	    *t = tnt;
	    if (iend) {
		*t = t1;
	    }
	}
	goto L23045;
L23047:
	wa[m2 + (*nn + 1) * wa_dim1] = 2.;
	*ift = 2;
	goto L50;
L40:
	if (*lsol > 2) {
	    sol[sol_dim1 + 1] = 0.;
	    sol[sol_dim1 + 2] = 0.;
	    sol[sol_dim1 + 3] = 0.;
	    sol[*lsol * sol_dim1 + 1] = 1.;
	    sol[*lsol * sol_dim1 + 2] = 0.;
	    sol[*lsol * sol_dim1 + 3] = 0.;
	    i1 = *m;
	    for (i = 1; i <= i1; ++i) {
		dsol[i + dsol_dim1] = 1.;
		dsol[i + *lsol * dsol_dim1] = 0.;
		dsol[i + (*lsol + 1) * dsol_dim1] = 0.;
	    }
	}
	l = kl - 1;
	i1 = l;
	for (i = 1; i <= i1; ++i) {
	    if (wa[i + n1 * wa_dim1] < 0.) {
		i2 = *n4;
		for (j = kr; j <= i2; ++j) {
		    wa[i + j * wa_dim1] = -wa[i + j * wa_dim1];
		}
	    }
	}
L50:
	sum = 0.;
	if (! lci2) {
	    i1 = *m;
	    for (i = kl; i <= i1; ++i) {
		k = (integer) (wa[i + *n4 * wa_dim1] * d_sign(&c_b15, &wa[
			i + *n4 * wa_dim1]));
		d__ = wa[i + n1 * wa_dim1] * d_sign(&c_b15, &wa[i + *n4 * 
			wa_dim1]);
		sum += d__ * d_sign(&c_b15, &d__) * (half + d_sign(&c_b15, &
			d__) * (*t - half));
		k -= n;
		e[k] = d__;
	    }
	    wa[m2 + n2 * wa_dim1] = (doublereal) kount;
	    wa[m1 + n2 * wa_dim1] = (doublereal) (n1 - kr);
	    wa[m1 + n1 * wa_dim1] = sum;
	}
	if (wa[m2 + (*nn + 1) * wa_dim1] == 2.) {
	    goto L23018;
	}
	if (! (*lci1)) {
	    goto L23018;
	}
	if (lci2) {
	    goto L23234;
	}
	lci2 = TRUE_;
	n = *nn - 1;
	n1 = n + 1;
	n2 = n + 2;
	*n3 = n + 3;
	*n4 = n + 4;
L60:
	++idxcf;
	if (! (idxcf > *nn)) {
	    goto L23236;
	}
	goto L23018;
L23236:
L70:
	if (! lup) {
	    goto L23238;
	}
	tnew = x[idxcf] + *toler;
	told = tnew;
	ci[(idxcf << 2) + 3] = x[idxcf];
	tnmat[(idxcf << 2) + 3] = 0.;
	goto L23239;
L23238:
	tnew = x[idxcf] - *toler;
	told = tnew;
	ci[(idxcf << 2) + 2] = x[idxcf];
	tnmat[(idxcf << 2) + 2] = 0.;
L23239:
L23234:
	goto L23016;
L23018:
	i1 = *m;
	for (i = 1; i <= i1; ++i) {
	    dsol[i + *lsol * dsol_dim1] = dsol[i + (*lsol + 1) * 
		    dsol_dim1];
	}
    }

    return 0;
}

