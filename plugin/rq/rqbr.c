/* rqbr.f -- translated by f2c (version 20061008) */

#include "gretl_f2c.h"

/* Table of constant values */

static doublereal c_b15 = 1.;

static double d_sign (doublereal *a, doublereal *b)
{
    double x;
    x = (*a >= 0 ? *a : - *a);
    return ( *b >= 0 ? x : -x);
}

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int rqbr_(integer *m, integer *nn, integer *m5, integer *n3, 
	integer *n4, doublereal *a, doublereal *b, doublereal *t, doublereal *
	toler, integer *ift, doublereal *x, doublereal *e, integer *s, 
	doublereal *wa, doublereal *wb, integer *nsol, integer *ndsol, 
	doublereal *sol, doublereal *dsol, integer *lsol, integer *h__, 
	doublereal *qn, doublereal *cutoff, doublereal *ci, doublereal *tnmat,
	 doublereal *big, logical *lci1)
{
    /* System generated locals */
    integer h_dim1, h_offset, sol_dim1, sol_offset, a_dim1, a_offset, wa_dim1,
	     wa_offset, dsol_dim1, dsol_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, l, n;
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
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    x[j] = 0.;
/* L23010: */
	}
/* L23011: */
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    e[i__] = 0.;
/* L23012: */
	}
/* L23013: */
	if (*t < 0. || *t > 1.) {
	    t0 = 1. / ((real) (*m) * 2.);
	    t1 = 1. - t0;
	    *t = t0;
	    iend = FALSE_;
	    *lci1 = FALSE_;
	}
L23016:
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = 1;
	    i__2 = *nn;
	    for (j = 1; j <= i__2; ++j) {
		if (k <= *nn) {
		    if (j == idxcf) {
			skip = TRUE_;
		    } else {
			skip = FALSE_;
		    }
		    if (! skip) {
			wa[i__ + k * wa_dim1] = a[i__ + j * a_dim1];
			++k;
		    }
		}
/* L23021: */
	    }
/* L23022: */
	    wa[i__ + *n4 * wa_dim1] = (doublereal) (n + i__);
	    wa[i__ + n2 * wa_dim1] = b[i__];
	    if (idxcf != 0) {
		wa[i__ + *n3 * wa_dim1] = tnew * a[i__ + idxcf * a_dim1];
	    } else {
		wa[i__ + *n3 * wa_dim1] = 0.;
	    }
	    wa[i__ + n1 * wa_dim1] = wa[i__ + n2 * wa_dim1] - wa[i__ + *n3 * 
		    wa_dim1];
	    if (wa[i__ + n1 * wa_dim1] < 0.) {
		i__2 = *n4;
		for (j = 1; j <= i__2; ++j) {
		    wa[i__ + j * wa_dim1] = -wa[i__ + j * wa_dim1];
/* L23033: */
		}
/* L23034: */
	    }
/* L23019: */
	}
/* L23020: */
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    wa[m4 + j * wa_dim1] = (doublereal) j;
	    wa[m2 + j * wa_dim1] = 0.;
	    wa[m3 + j * wa_dim1] = 0.;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		aux = d_sign(&c_b15, &wa[m4 + j * wa_dim1]) * wa[i__ + j * 
			wa_dim1];
		wa[m2 + j * wa_dim1] += aux * (1. - d_sign(&c_b15, &wa[i__ + *
			n4 * wa_dim1]));
		wa[m3 + j * wa_dim1] += aux * d_sign(&c_b15, &wa[i__ + *n4 * 
			wa_dim1]);
/* L23037: */
	    }
/* L23038: */
	    wa[m3 + j * wa_dim1] *= 2.;
/* L23035: */
	}
/* L23036: */
	dif = 0.;
	init = FALSE_;
	if (! lci2) {
	    i__1 = n;
	    for (k = 1; k <= i__1; ++k) {
		wa[*m5 + k * wa_dim1] = 0.;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    wa[*m5 + k * wa_dim1] += a[i__ + k * a_dim1];
/* L23043: */
		}
/* L23044: */
		wa[*m5 + k * wa_dim1] /= (real) (*m);
/* L23041: */
	    }
/* L23042: */
	}
	*lsol = 1;
	kount = 0;
L23045:
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    wa[m1 + j * wa_dim1] = wa[m2 + j * wa_dim1] + wa[m3 + j * wa_dim1]
		     * *t;
/* L23048: */
	}
/* L23049: */
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
	i__1 = n;
	for (j = kr; j <= i__1; ++j) {
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
/* L23059: */
	if (max__ <= *toler) {
	    goto L23054;
	}
	if (wa[m1 + in * wa_dim1] <= 0.) {
	    i__1 = m4;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		wa[i__ + in * wa_dim1] = -wa[i__ + in * wa_dim1];
/* L23070: */
	    }
/* L23071: */
	    wa[m1 + in * wa_dim1] += -2.;
	    wa[m2 + in * wa_dim1] += -2.;
	}
L23072:
	k = 0;
	i__1 = *m;
	for (i__ = kl; i__ <= i__1; ++i__) {
	    d__ = wa[i__ + in * wa_dim1];
	    if (d__ > *toler) {
		++k;
		wb[k] = wa[i__ + n1 * wa_dim1] / d__;
		s[k] = i__;
		test = TRUE_;
	    }
/* L23075: */
	}
/* L23076: */
L23079:
	if (k <= 0) {
	    test = FALSE_;
	} else {
	    min__ = *big;
	    i__1 = k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (wb[i__] < min__) {
		    j = i__;
		    min__ = wb[i__];
		    out = s[i__];
		}
/* L23084: */
	    }
/* L23085: */
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
	i__1 = *n3;
	for (j = kr; j <= i__1; ++j) {
	    d__ = wa[out + j * wa_dim1];
	    wa[m1 + j * wa_dim1] = wa[m1 + j * wa_dim1] - d__ - d__;
	    wa[m2 + j * wa_dim1] = wa[m2 + j * wa_dim1] - d__ - d__;
	    wa[out + j * wa_dim1] = -d__;
/* L23094: */
	}
/* L23095: */
	wa[out + *n4 * wa_dim1] = -wa[out + *n4 * wa_dim1];
/* L23080: */
	goto L23079;
L23081:
	i__1 = m4;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__ = wa[i__ + kr * wa_dim1];
	    wa[i__ + kr * wa_dim1] = wa[i__ + in * wa_dim1];
	    wa[i__ + in * wa_dim1] = d__;
/* L23096: */
	}
/* L23097: */
	++kr;
	goto L20;
L10:
	i__1 = *n3;
	for (j = kr; j <= i__1; ++j) {
	    if (j != in) {
		wa[out + j * wa_dim1] /= pivot;
	    }
/* L23098: */
	}
/* L23099: */
	i__1 = m3;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ != out) {
		d__ = wa[i__ + in * wa_dim1];
		i__2 = *n3;
		for (j = kr; j <= i__2; ++j) {
		    if (j != in) {
			wa[i__ + j * wa_dim1] -= d__ * wa[out + j * wa_dim1];
		    }
/* L23106: */
		}
/* L23107: */
	    }
/* L23102: */
	}
/* L23103: */
	i__1 = m3;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ != out) {
		wa[i__ + in * wa_dim1] = -wa[i__ + in * wa_dim1] / pivot;
	    }
/* L23110: */
	}
/* L23111: */
	wa[out + in * wa_dim1] = 1. / pivot;
	d__ = wa[out + *n4 * wa_dim1];
	wa[out + *n4 * wa_dim1] = wa[m4 + in * wa_dim1];
	wa[m4 + in * wa_dim1] = d__;
	++kount;
	if (! stage) {
	    goto L23074;
	}
	++kl;
	i__1 = *n4;
	for (j = kr; j <= i__1; ++j) {
	    d__ = wa[out + j * wa_dim1];
	    wa[out + j * wa_dim1] = wa[kount + j * wa_dim1];
	    wa[kount + j * wa_dim1] = d__;
/* L23116: */
	}
/* L23117: */
L20:
	if (kount + kr == n1) {
	    goto L23057;
	}
L30:
	max__ = -1.;
	i__1 = n;
	for (j = kr; j <= i__1; ++j) {
	    if ((d__1 = wa[m4 + j * wa_dim1], abs(d__1)) <= (doublereal) n) {
		d__ = (d__1 = wa[m1 + j * wa_dim1], abs(d__1));
		if (d__ > max__) {
		    max__ = d__;
		    in = j;
		}
	    }
/* L23120: */
	}
/* L23121: */
	if (wa[m1 + in * wa_dim1] < 0.) {
	    i__1 = m4;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		wa[i__ + in * wa_dim1] = -wa[i__ + in * wa_dim1];
/* L23128: */
	    }
/* L23129: */
	}
/* L23073: */
	goto L23072;
L23074:
/* L23056: */
	goto L23055;
L23057:
/* L23053: */
	goto L23052;
L23054:
	if (kr == 1) {
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		d__ = (d__1 = wa[m1 + j * wa_dim1], abs(d__1));
		if (d__ <= *toler || 2. - d__ <= *toler) {
		    *ift = 1;
		    wa[m2 + (*nn + 1) * wa_dim1] = 0.;
		    goto L80;
		}
/* L23132: */
	    }
/* L23133: */
	}
L80:
	kount = 0;
	sum = 0.;
	if (! lci2) {
	    i__1 = kl - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		k = (integer) (wa[i__ + *n4 * wa_dim1] * d_sign(&c_b15, &wa[
			i__ + *n4 * wa_dim1]));
		x[k] = wa[i__ + n1 * wa_dim1] * d_sign(&c_b15, &wa[i__ + *n4 *
			 wa_dim1]);
/* L23138: */
	    }
/* L23139: */
	}
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    kd = (integer) ((d__1 = wa[m4 + i__ * wa_dim1], abs(d__1)) - n);
	    dsol[kd + *lsol * dsol_dim1] = wa[m1 + i__ * wa_dim1] / 2. + 1.;
	    if (wa[m4 + i__ * wa_dim1] < 0.) {
		dsol[kd + *lsol * dsol_dim1] = 1. - dsol[kd + *lsol * 
			dsol_dim1];
	    }
	    if (! lci2) {
		sum += x[i__] * wa[*m5 + i__ * wa_dim1];
		sol[i__ + 3 + *lsol * sol_dim1] = x[i__];
		h__[i__ + *lsol * h_dim1] = kd;
	    }
/* L23140: */
	}
/* L23141: */
	i__1 = *m;
	for (i__ = kl; i__ <= i__1; ++i__) {
	    kd = (integer) ((d__1 = wa[i__ + *n4 * wa_dim1], abs(d__1)) - n);
	    if (wa[i__ + *n4 * wa_dim1] < 0.) {
		dsol[kd + *lsol * dsol_dim1] = 0.;
	    }
	    if (wa[i__ + *n4 * wa_dim1] > 0.) {
		dsol[kd + *lsol * dsol_dim1] = 1.;
	    }
/* L23146: */
	}
/* L23147: */
	if (! lci2) {
	    sol[*lsol * sol_dim1 + 1] = smax;
	    sol[*lsol * sol_dim1 + 2] = sum;
	    sum = 0.;
	    i__1 = *m;
	    for (j = kl; j <= i__1; ++j) {
		d__ = wa[j + n1 * wa_dim1] * d_sign(&c_b15, &wa[j + *n4 * 
			wa_dim1]);
		sum += d__ * (smax + half * (d_sign(&c_b15, &d__) - 1.));
/* L23154: */
	    }
/* L23155: */
	    sol[*lsol * sol_dim1 + 3] = sum;
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dsol[i__ + (*lsol + 1) * dsol_dim1] = dsol[i__ + *lsol * 
			dsol_dim1];
/* L23156: */
	    }
/* L23157: */
	}
	if (lci2) {
	    a1 = 0.;
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a1 += a[i__ + idxcf * a_dim1] * (dsol[i__ + *lsol * dsol_dim1]
			 + *t - 1.);
/* L23160: */
	    }
/* L23161: */
	    tn = a1 / sqrt(qn[idxcf] * *t * (1. - *t));
	    if (abs(tn) < *cutoff) {
		if (lup) {
		    smax = *big;
		} else {
		    smax = -(*big);
		}
		i__1 = kl - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    k = (integer) (wa[i__ + *n4 * wa_dim1] * d_sign(&c_b15, &
			    wa[i__ + *n4 * wa_dim1]));
		    sol[k + sol_dim1] = wa[i__ + n2 * wa_dim1] * d_sign(&
			    c_b15, &wa[i__ + *n4 * wa_dim1]);
		    sol[k + (sol_dim1 << 1)] = wa[i__ + *n3 * wa_dim1] * 
			    d_sign(&c_b15, &wa[i__ + *n4 * wa_dim1]) / tnew;
/* L23166: */
		}
/* L23167: */
		i__1 = *m;
		for (i__ = kl; i__ <= i__1; ++i__) {
		    a1 = 0.;
		    b1 = 0.;
		    k = (integer) (wa[i__ + *n4 * wa_dim1] * d_sign(&c_b15, &
			    wa[i__ + *n4 * wa_dim1]) - n);
		    l = 1;
		    i__2 = n;
		    for (j = 1; j <= i__2; ++j) {
			if (j == idxcf) {
			    ++l;
			}
			a1 += a[k + l * a_dim1] * sol[j + sol_dim1];
			b1 += a[k + l * a_dim1] * sol[j + (sol_dim1 << 1)];
			++l;
/* L23170: */
		    }
/* L23171: */
		    tnt = (b[k] - a1) / (a[k + idxcf * a_dim1] - b1);
		    if (lup) {
			if (tnt > tnew) {
			    if (tnt < smax) {
				smax = tnt;
				out = i__;
			    }
			}
		    } else {
			if (tnt < tnew) {
			    if (tnt > smax) {
				smax = tnt;
				out = i__;
			    }
			}
		    }
/* L23168: */
		}
/* L23169: */
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
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    wa[i__ + *n3 * wa_dim1] = wa[i__ + *n3 * wa_dim1] / told *
			     tnew;
		    wa[i__ + n1 * wa_dim1] = wa[i__ + n2 * wa_dim1] - wa[i__ 
			    + *n3 * wa_dim1];
/* L23190: */
		}
/* L23191: */
		i__1 = *n3;
		for (j = kr; j <= i__1; ++j) {
		    d__ = wa[out + j * wa_dim1];
		    wa[m1 + j * wa_dim1] = wa[m1 + j * wa_dim1] - d__ - d__;
		    wa[m2 + j * wa_dim1] = wa[m2 + j * wa_dim1] - d__ - d__;
		    wa[out + j * wa_dim1] = -d__;
/* L23192: */
		}
/* L23193: */
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
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s[i__] = 0;
/* L23200: */
	    }
/* L23201: */
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		x[j] = 0.;
/* L23202: */
	    }
/* L23203: */
	    smax = 2.;
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
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
/* L23204: */
	    }
/* L23205: */
	    tnt = smax + *toler * (abs(dif) + 1.);
	    if (tnt >= t1 + *toler) {
		iend = TRUE_;
	    }
	    *t = tnt;
	    if (iend) {
		*t = t1;
	    }
	}
/* L23046: */
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
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dsol[i__ + dsol_dim1] = 1.;
		dsol[i__ + *lsol * dsol_dim1] = 0.;
		dsol[i__ + (*lsol + 1) * dsol_dim1] = 0.;
/* L23220: */
	    }
/* L23221: */
	}
	l = kl - 1;
	i__1 = l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (wa[i__ + n1 * wa_dim1] < 0.) {
		i__2 = *n4;
		for (j = kr; j <= i__2; ++j) {
		    wa[i__ + j * wa_dim1] = -wa[i__ + j * wa_dim1];
/* L23226: */
		}
/* L23227: */
	    }
/* L23222: */
	}
/* L23223: */
L50:
	sum = 0.;
	if (! lci2) {
	    i__1 = *m;
	    for (i__ = kl; i__ <= i__1; ++i__) {
		k = (integer) (wa[i__ + *n4 * wa_dim1] * d_sign(&c_b15, &wa[
			i__ + *n4 * wa_dim1]));
		d__ = wa[i__ + n1 * wa_dim1] * d_sign(&c_b15, &wa[i__ + *n4 * 
			wa_dim1]);
		sum += d__ * d_sign(&c_b15, &d__) * (half + d_sign(&c_b15, &
			d__) * (*t - half));
		k -= n;
		e[k] = d__;
/* L23230: */
	    }
/* L23231: */
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
/* L23017: */
	goto L23016;
L23018:
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dsol[i__ + *lsol * dsol_dim1] = dsol[i__ + (*lsol + 1) * 
		    dsol_dim1];
/* L23242: */
	}
/* L23243: */
    }
    return 0;
} /* rqbr_ */

