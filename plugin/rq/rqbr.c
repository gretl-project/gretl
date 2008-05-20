/* rqbr.f -- translated by f2c (version 20061008) then
   somewhat cleaned up by Allin Cottrell.
   Original RATFOR code by Roger Koenker.
*/

#include <stdio.h>
#include <math.h>

#define sgn(x) (((x) >= 0)? 1.0 : -1.0)

int rqbr_(int m, int nn, double *a, double *b, double tau, 
	  double tol, double *x, double *e, 
	  int *s, double *wa, double *wb, 
	  double *sol, double *dsol, int *h, 
	  double *qn, double cutoff, double *ci, double *tnmat,
	  double big, int rmax, int ci1)
{
    static double d, a1, b1;
    static int i, j, k, l;
    static int m1, m2, m3, m4, n1, n2;
    static int kd, kl, in, kr;
    static double t1, tn, dif, dmin, dmax, aux;
    static double sum, tnt, told;
    static char init, skip, iend;
    static double smax, tnew;
    static char lup, test, stage;
    static int idxcf, kount, out;
    static double pivot;
    int n = nn;
    int n3 = n + 3, n4 = n + 4, m5 = m + 5;
    int sdim, ci2 = 0;
    int lsol, ift = 0, rcount = 0;

    /* Parameter adjustments (ugh!) */
    --wb;
    --s;
    --e;
    --b;
    tnmat -= 5;
    ci -= 5;
    --qn;
    --x;
 
    a -= (1 + m);
    wa -= (1 + m5);
    h -= (1 + nn);
    sdim = n3;
    sol -= (1 + sdim);
    dsol -= (1 + m);

    /* Function Body */
    n = nn;
    ift = 0;
    wa[m + 2 + (nn + 1) * m5] = 1.0;

    iend = 1;
    ci2 = 0;
    lup = 1;
    skip = 0;
    idxcf = 0;
    tnew = 0;
    tn = 0;

    m1 = m + 1;
    n1 = n + 1;
    n2 = n + 2;
    m2 = m + 2;
    m3 = m + 3;
    m4 = m + 4;

    for (j = 1; j <= n; ++j) {
	x[j] = 0;
    }
    for (i = 1; i <= m; ++i) {
	e[i] = 0;
    }

 L23016:
    for (i = 1; i <= m; ++i) {
	k = 1;
	for (j = 1; j <= nn; ++j) {
	    if (k <= nn) {
		if (j == idxcf) {
		    skip = 1;
		} else {
		    skip = 0;
		}
		if (! skip) {
		    wa[i + k * m5] = a[i + j * m];
		    ++k;
		}
	    }
	}
	wa[i + n4 * m5] = n + i;
	wa[i + n2 * m5] = b[i];
	if (idxcf != 0) {
	    wa[i + n3 * m5] = tnew * a[i + idxcf * m];
	} else {
	    wa[i + n3 * m5] = 0;
	}
	wa[i + n1 * m5] = wa[i + n2 * m5] - wa[i + n3 * m5];
	if (wa[i + n1 * m5] < 0) {
	    for (j = 1; j <= n4; ++j) {
		wa[i + j * m5] = -wa[i + j * m5];
	    }
	}
    }
    for (j = 1; j <= n; ++j) {
	wa[m4 + j * m5] = j;
	wa[m2 + j * m5] = 0;
	wa[m3 + j * m5] = 0;
	for (i = 1; i <= m; ++i) {
	    aux = sgn(wa[m4 + j * m5]) * wa[i + j * m5];
	    wa[m2 + j * m5] += aux * (1. - sgn(wa[i + n4 * m5]));
	    wa[m3 + j * m5] += aux * sgn(wa[i + n4 * m5]);
	}
	wa[m3 + j * m5] *= 2.;
    }
    dif = 0;
    init = 0;
    if (!ci2) {
	for (k = 1; k <= n; ++k) {
	    wa[m5 + k * m5] = 0;
	    for (i = 1; i <= m; ++i) {
		wa[m5 + k * m5] += a[i + k * m];
	    }
	    wa[m5 + k * m5] /= m;
	}
    }
    lsol = 1;
    kount = 0;
 L23045:
    for (j = 1; j <= n; ++j) {
	wa[m1 + j * m5] = wa[m2 + j * m5] + wa[m3 + j * m5] * tau;
    }
    if (! init) {
	stage = 1;
	kr = 1;
	kl = 1;
	goto L30;
    }
 L23052:
    stage = 0;
 L23055:
    dmax = -(big);
    for (j = kr; j <= n; ++j) {
	d = wa[m1 + j * m5];
	if (d < 0) {
	    if (d > -2.) {
		goto L23058;
	    }
	    d = -d - 2.;
	}
	if (d > dmax) {
	    dmax = d;
	    in = j;
	}
    L23058:
	;
    }
    if (dmax <= tol) {
	if (++rcount > rmax) {
	    fprintf(stderr, "max=%g, tol=%g, rcount=%d: stopping\n", 
		    dmax, tol, rcount);
	    ift = 3;
	    return ift;
	}
	goto L23054;
    }
    if (wa[m1 + in * m5] <= 0) {
	for (i = 1; i <= m4; ++i) {
	    wa[i + in * m5] = -wa[i + in * m5];
	}
	wa[m1 + in * m5] += -2.;
	wa[m2 + in * m5] += -2.;
    }
 L23072:
    k = 0;
    for (i = kl; i <= m; ++i) {
	d = wa[i + in * m5];
	if (d > tol) {
	    ++k;
	    wb[k] = wa[i + n1 * m5] / d;
	    s[k] = i;
	    test = 1;
	}
    }
 L23079:
    if (k <= 0) {
	test = 0;
    } else {
	dmin = big;
	for (i = 1; i <= k; ++i) {
	    if (wb[i] < dmin) {
		j = i;
		dmin = wb[i];
		out = s[i];
	    }
	}
	wb[j] = wb[k];
	s[j] = s[k];
	--k;
    }
    if (!test && stage) {
	goto L23081;
    }
    if (!test) {
	goto L23047;
    }
    pivot = wa[out + in * m5];
    if (wa[m1 + in * m5] - pivot - pivot <= tol) {
	goto L10;
    }
    for (j = kr; j <= n3; ++j) {
	d = wa[out + j * m5];
	wa[m1 + j * m5] = wa[m1 + j * m5] - d - d;
	wa[m2 + j * m5] = wa[m2 + j * m5] - d - d;
	wa[out + j * m5] = -d;
    }
    wa[out + n4 * m5] = -wa[out + n4 * m5];
    goto L23079;
 L23081:
    for (i = 1; i <= m4; ++i) {
	d = wa[i + kr * m5];
	wa[i + kr * m5] = wa[i + in * m5];
	wa[i + in * m5] = d;
    }
    ++kr;
    goto L20;
 L10:
    for (j = kr; j <= n3; ++j) {
	if (j != in) {
	    wa[out + j * m5] /= pivot;
	}
    }
    for (i = 1; i <= m3; ++i) {
	if (i != out) {
	    d = wa[i + in * m5];
	    for (j = kr; j <= n3; ++j) {
		if (j != in) {
		    wa[i + j * m5] -= d * wa[out + j * m5];
		}
	    }
	}
    }
    for (i = 1; i <= m3; ++i) {
	if (i != out) {
	    wa[i + in * m5] = -wa[i + in * m5] / pivot;
	}
    }
    wa[out + in * m5] = 1.0 / pivot;
    d = wa[out + n4 * m5];
    wa[out + n4 * m5] = wa[m4 + in * m5];
    wa[m4 + in * m5] = d;
    ++kount;
    if (!stage) {
	goto L23074;
    }
    ++kl;
    for (j = kr; j <= n4; ++j) {
	d = wa[out + j * m5];
	wa[out + j * m5] = wa[kount + j * m5];
	wa[kount + j * m5] = d;
    }
 L20:
    if (kount + kr == n1) {
	goto L23057;
    }
 L30:
    dmax = -1.;
    for (j = kr; j <= n; ++j) {
	if (fabs(wa[m4 + j * m5]) <= n) {
	    d = fabs(wa[m1 + j * m5]);
	    if (d > dmax) {
		dmax = d;
		in = j;
	    }
	}
    }
    if (wa[m1 + in * m5] < 0) {
	for (i = 1; i <= m4; ++i) {
	    wa[i + in * m5] = -wa[i + in * m5];
	}
    }
    goto L23072;
 L23074:
    goto L23055;
 L23057:
    goto L23052;
 L23054:
    if (kr == 1) {
	for (j = 1; j <= n; ++j) {
	    d = fabs(wa[m1 + j * m5]);
	    if (d <= tol || 2.0 - d <= tol) {
		ift = 1;
		wa[m2 + (nn + 1) * m5] = 0;
		goto L80;
	    }
	}
    }
 L80:
    kount = 0;
    sum = 0;
    if (!ci2) {
	for (i = 1; i <= kl - 1; ++i) {
	    k = (int) (wa[i + n4 * m5] * sgn(wa[i + n4 * m5]));
	    x[k] = wa[i + n1 * m5] * sgn(wa[i + n4 * m5]);
	}
    }
    for (i = 1; i <= n; ++i) {
	kd = fabs(wa[m4 + i * m5]) - n;
	dsol[kd + lsol * m] = wa[m1 + i * m5] / 2.0 + 1.0;
	if (wa[m4 + i * m5] < 0) {
	    dsol[kd + lsol * m] = 1. - dsol[kd + lsol * m];
	}
	if (!ci2) {
	    sum += x[i] * wa[m5 + i * m5];
	    sol[i + 3 + lsol * sdim] = x[i];
	    h[i + lsol * nn] = kd;
	}
    }
    for (i = kl; i <= m; ++i) {
	kd = fabs(wa[i + n4 * m5]) - n;
	if (wa[i + n4 * m5] < 0) {
	    dsol[kd + lsol * m] = 0;
	}
	if (wa[i + n4 * m5] > 0) {
	    dsol[kd + lsol * m] = 1.;
	}
    }
    if (!ci2) {
	sol[lsol * sdim + 1] = smax;
	sol[lsol * sdim + 2] = sum;
	sum = 0;
	for (j = kl; j <= m; ++j) {
	    d = wa[j + n1 * m5] * sgn(wa[j + n4 * m5]);
	    sum += d * (smax + .5 * (sgn(d) - 1.0));
	}
	sol[lsol * sdim + 3] = sum;
	for (i = 1; i <= m; ++i) {
	    dsol[i + (lsol + 1) * m] = dsol[i + lsol * m];
	}
    }
    if (ci2) {
	a1 = 0;
	for (i = 1; i <= m; ++i) {
	    a1 += a[i + idxcf * m] * (dsol[i + lsol * m] + tau - 1);
	}
	tn = a1 / sqrt(qn[idxcf] * tau * (1. - tau));
	if (fabs(tn) < cutoff) {
	    smax = (lup)? big : -big;
	    for (i = 1; i <= kl - 1; ++i) {
		k = (int) (wa[i + n4 * m5] * sgn(wa[i + n4 * m5]));
		sol[k + sdim] = wa[i + n2 * m5] * sgn(wa[i + n4 * m5]);
		sol[k + (sdim << 1)] = wa[i + n3 * m5] * 
		    sgn(wa[i + n4 * m5]) / tnew;
	    }
	    for (i = kl; i <= m; ++i) {
		a1 = 0;
		b1 = 0;
		k = (int) (wa[i + n4 * m5] * sgn(wa[i + n4 * m5]) - n);
		l = 1;
		for (j = 1; j <= n; ++j) {
		    if (j == idxcf) {
			++l;
		    }
		    a1 += a[k + l * m] * sol[j + sdim];
		    b1 += a[k + l * m] * sol[j + (sdim << 1)];
		    ++l;
		}
		tnt = (b[k] - a1) / (a[k + idxcf * m] - b1);
		if (lup) {
		    if (tnt > tnew && tnt < smax) {
			smax = tnt;
			out = i;
		    }
		} else if (tnt < tnew && tnt > smax) {
		    smax = tnt;
		    out = i;
		}
	    }
	    if (lup) {
		told = tnew;
		tnew = smax + tol;
		ci[(idxcf << 2) + 3] = told - tol;
		tnmat[(idxcf << 2) + 3] = tn;
		if (! (tnew < big - tol)) {
		    ci[(idxcf << 2) + 3] = big;
		    ci[(idxcf << 2) + 4] = big;
		    tnmat[(idxcf << 2) + 3] = tn;
		    tnmat[(idxcf << 2) + 4] = tn;
		    lup = 0;
		    goto L70;
		}
	    } else {
		told = tnew;
		tnew = smax - tol;
		ci[(idxcf << 2) + 2] = told + tol;
		tnmat[(idxcf << 2) + 2] = tn;
		if (! (tnew > -(big) + tol)) {
		    ci[(idxcf << 2) + 2] = -(big);
		    ci[(idxcf << 2) + 1] = -(big);
		    tnmat[(idxcf << 2) + 2] = tn;
		    tnmat[(idxcf << 2) + 1] = tn;
		    lup = 1;
		    goto L60;
		}
	    }
	    for (i = 1; i <= m; ++i) {
		wa[i + n3 * m5] = wa[i + n3 * m5] / told * tnew;
		wa[i + n1 * m5] = wa[i + n2 * m5] - wa[i + n3 * m5];
	    }
	    for (j = kr; j <= n3; ++j) {
		d = wa[out + j * m5];
		wa[m1 + j * m5] = wa[m1 + j * m5] - d - d;
		wa[m2 + j * m5] = wa[m2 + j * m5] - d - d;
		wa[out + j * m5] = -d;
	    }
	    wa[out + n4 * m5] = -wa[out + n4 * m5];
	    init = 1;
	} else {
	    if (lup) {
		ci[(idxcf << 2) + 4] = tnew - tol;
		tnmat[(idxcf << 2) + 4] = tn;
		lup = 0;
		goto L70;
	    } else {
		ci[(idxcf << 2) + 1] = tnew + tol;
		tnmat[(idxcf << 2) + 1] = tn;
		lup = 1;
		goto L60;
	    }
	}
    }
    if (iend && !ci2) {
	goto L40;
    }
    if (!ci2) {
	init = 1;
	lsol++;
	for (i = 1; i <= m; ++i) {
	    s[i] = 0;
	}
	for (j = 1; j <= n; ++j) {
	    x[j] = 0;
	}
	smax = 2.0;
	for (j = 1; j <= n; ++j) {
	    b1 = wa[m3 + j * m5];
	    a1 = (-2.0 - wa[m2 + j * m5]) / b1;
	    b1 = -wa[m2 + j * m5] / b1;
	    if (a1 >= tau) {
		if (a1 < smax) {
		    smax = a1;
		    dif = (b1 - a1) / 2.0;
		}
	    }
	    if (b1 > tau) {
		if (b1 < smax) {
		    smax = b1;
		    dif = (b1 - a1) / 2.0;
		}
	    }
	}
	tnt = smax + tol * (fabs(dif) + 1.0);
	if (tnt >= t1 + tol) {
	    iend = 1;
	}
	tau = tnt;
	fprintf(stderr, "set tau = tnt = %g\n", tau);
	if (iend) {
	    tau = t1;
	}
    }
    goto L23045;
 L23047:
    wa[m2 + (nn + 1) * m5] = 2.0;
    ift = 2;
    goto L50;
 L40:
    if (lsol > 2) {
	sol[sdim + 1] = 0;
	sol[sdim + 2] = 0;
	sol[sdim + 3] = 0;
	sol[lsol * sdim + 1] = 1.;
	sol[lsol * sdim + 2] = 0;
	sol[lsol * sdim + 3] = 0;
	for (i = 1; i <= m; ++i) {
	    dsol[i + m] = 1.;
	    dsol[i + lsol * m] = 0;
	    dsol[i + (lsol + 1) * m] = 0;
	}
    }
    l = kl - 1;
    for (i = 1; i <= l; ++i) {
	if (wa[i + n1 * m5] < 0) {
	    for (j = kr; j <= n4; ++j) {
		wa[i + j * m5] = -wa[i + j * m5];
	    }
	}
    }
 L50:
    sum = 0;
    if (!ci2) {
	for (i = kl; i <= m; ++i) {
	    k = (int) (wa[i + n4 * m5] * sgn(wa[i + n4 * m5]));
	    d = wa[i + n1 * m5] * sgn(wa[i + n4 * m5]);
	    sum += d * sgn(d) * (.5 + sgn(d) * (tau - .5));
	    k -= n;
	    e[k] = d;
	}
	wa[m2 + n2 * m5] = kount;
	wa[m1 + n2 * m5] = n1 - kr;
	wa[m1 + n1 * m5] = sum;
    }
    if (wa[m2 + (nn + 1) * m5] == 2.0) {
	goto L23018;
    }
    if (!ci1) {
	goto L23018;
    }
    if (ci2) {
	goto L23234;
    }
    ci2 = 1;
    n = nn - 1;
    n1 = n + 1;
    n2 = n + 2;
    n3 = n + 3;
    n4 = n + 4;
 L60:
    ++idxcf;
    if (!(idxcf > nn)) {
	goto L23236;
    }
    goto L23018;
 L23236:
 L70:
    if (!lup) {
	goto L23238;
    }
    tnew = x[idxcf] + tol;
    told = tnew;
    ci[(idxcf << 2) + 3] = x[idxcf];
    tnmat[(idxcf << 2) + 3] = 0;
    goto L23239;
 L23238:
    tnew = x[idxcf] - tol;
    told = tnew;
    ci[(idxcf << 2) + 2] = x[idxcf];
    tnmat[(idxcf << 2) + 2] = 0;
 L23239:
 L23234:
    goto L23016;
 L23018:
    for (i = 1; i <= m; ++i) {
	dsol[i + lsol * m] = dsol[i + (lsol + 1) * m];
    }

    return ift;
}

