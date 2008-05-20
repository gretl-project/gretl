/* rqbr.f -- translated by f2c (version 20061008) then
   somewhat cleaned up by Allin Cottrell.
   Original RATFOR code by Roger Koenker.
*/

#include <stdio.h>
#include <math.h>

#define sgn(x) (((x) >= 0)? 1.0 : -1.0)

int rqbr_(int n, int pp, double *x, double *y, double tau, 
	  double tol, double *coeff, double *resid, 
	  int *s, double *wa, double *wb, 
	  double *sol, double *dsol, int *h, 
	  double *qn, double cutoff, double *ci, double *tnmat,
	  double big, int rmax, int ci1)
{
    static double d, a1, b1;
    static int i, j, k, l;
    static int n1, n2, n3, n4, p1, p2;
    static int kd, kl, in, kr;
    static double t1, tn, dif, dmin, dmax, aux;
    static double sum, tnt, told;
    static char init, iend;
    static double smax, tnew;
    static char lup, test, stage;
    static int idxcf, kount, out;
    static double pivot;
    int p = pp;
    int p3 = p + 3, p4 = p + 4, n5 = n + 5;
    int sdim, ci2 = 0;
    int lsol, ift = 0, rcount = 0;

    /* Parameter adjustments (ugh!) */
    --wb;
    --s;
    --resid;
    --y;
    tnmat -= 5;
    ci -= 5;
    --qn;
    --coeff;
 
    x -= (1 + n);
    wa -= (1 + n5);
    h -= (1 + pp);
    sdim = p3;
    sol -= (1 + sdim);
    dsol -= (1 + n);

    /* Function Body */
    p = pp;
    ift = 0;
    wa[n + 2 + (pp + 1) * n5] = 1.0;

    iend = 1;
    ci2 = 0;
    lup = 1;
    idxcf = 0;
    tnew = 0;
    tn = 0;

    p1 = p + 1;
    p2 = p + 2;
    n1 = n + 1;
    n2 = n + 2;
    n3 = n + 3;
    n4 = n + 4;

    for (j = 1; j <= p; ++j) {
	coeff[j] = 0;
    }
    for (i = 1; i <= n; ++i) {
	resid[i] = 0;
    }

 L23016:
    for (i = 1; i <= n; ++i) {
	k = 1;
	for (j = 1; j <= pp; ++j) {
	    if (k <= pp && j != idxcf) {
		wa[i + k * n5] = x[i + j * n];
		++k;
	    }
	}
	wa[i + p4 * n5] = p + i;
	wa[i + p2 * n5] = y[i];
	wa[i + p3 * n5] = (idxcf == 0)? 0 : tnew * x[i + idxcf * n];
	wa[i + p1 * n5] = wa[i + p2 * n5] - wa[i + p3 * n5];
	if (wa[i + p1 * n5] < 0) {
	    for (j = 1; j <= p4; ++j) {
		wa[i + j * n5] = -wa[i + j * n5];
	    }
	}
    }
    for (j = 1; j <= p; ++j) {
	wa[n4 + j * n5] = j;
	wa[n2 + j * n5] = 0;
	wa[n3 + j * n5] = 0;
	for (i = 1; i <= n; ++i) {
	    aux = sgn(wa[n4 + j * n5]) * wa[i + j * n5];
	    wa[n2 + j * n5] += aux * (1. - sgn(wa[i + p4 * n5]));
	    wa[n3 + j * n5] += aux * sgn(wa[i + p4 * n5]);
	}
	wa[n3 + j * n5] *= 2.0;
    }
    dif = 0;
    init = 0;
    if (!ci2) {
	for (k = 1; k <= p; ++k) {
	    wa[n5 + k * n5] = 0;
	    for (i = 1; i <= n; ++i) {
		wa[n5 + k * n5] += x[i + k * n];
	    }
	    wa[n5 + k * n5] /= n;
	}
    }
    lsol = 1;
    kount = 0;
 L23045:
    for (j = 1; j <= p; ++j) {
	wa[n1 + j * n5] = wa[n2 + j * n5] + wa[n3 + j * n5] * tau;
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
    for (j = kr; j <= p; ++j) {
	d = wa[n1 + j * n5];
	if (d < 0) {
	    if (d > -2.0) {
		continue;
	    }
	    d = -d - 2.0;
	}
	if (d > dmax) {
	    dmax = d;
	    in = j;
	}
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
    if (wa[n1 + in * n5] <= 0) {
	for (i = 1; i <= n4; ++i) {
	    wa[i + in * n5] = -wa[i + in * n5];
	}
	wa[n1 + in * n5] += -2.;
	wa[n2 + in * n5] += -2.;
    }
 L23072:
    k = 0;
    for (i = kl; i <= n; ++i) {
	d = wa[i + in * n5];
	if (d > tol) {
	    ++k;
	    wb[k] = wa[i + p1 * n5] / d;
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
    pivot = wa[out + in * n5];
    if (wa[n1 + in * n5] - pivot - pivot <= tol) {
	goto L10;
    }
    for (j = kr; j <= p3; ++j) {
	d = wa[out + j * n5];
	wa[n1 + j * n5] = wa[n1 + j * n5] - d - d;
	wa[n2 + j * n5] = wa[n2 + j * n5] - d - d;
	wa[out + j * n5] = -d;
    }
    wa[out + p4 * n5] = -wa[out + p4 * n5];
    goto L23079;
 L23081:
    for (i = 1; i <= n4; ++i) {
	d = wa[i + kr * n5];
	wa[i + kr * n5] = wa[i + in * n5];
	wa[i + in * n5] = d;
    }
    ++kr;
    goto L20;
 L10:
    for (j = kr; j <= p3; ++j) {
	if (j != in) {
	    wa[out + j * n5] /= pivot;
	}
    }
    for (i = 1; i <= n3; ++i) {
	if (i != out) {
	    d = wa[i + in * n5];
	    for (j = kr; j <= p3; ++j) {
		if (j != in) {
		    wa[i + j * n5] -= d * wa[out + j * n5];
		}
	    }
	}
    }
    for (i = 1; i <= n3; ++i) {
	if (i != out) {
	    wa[i + in * n5] = -wa[i + in * n5] / pivot;
	}
    }
    wa[out + in * n5] = 1.0 / pivot;
    d = wa[out + p4 * n5];
    wa[out + p4 * n5] = wa[n4 + in * n5];
    wa[n4 + in * n5] = d;
    ++kount;
    if (!stage) {
	goto L23074;
    }
    ++kl;
    for (j = kr; j <= p4; ++j) {
	d = wa[out + j * n5];
	wa[out + j * n5] = wa[kount + j * n5];
	wa[kount + j * n5] = d;
    }
 L20:
    if (kount + kr == p1) {
	goto L23057;
    }
 L30:
    dmax = -1.;
    for (j = kr; j <= p; ++j) {
	if (fabs(wa[n4 + j * n5]) <= p) {
	    d = fabs(wa[n1 + j * n5]);
	    if (d > dmax) {
		dmax = d;
		in = j;
	    }
	}
    }
    if (wa[n1 + in * n5] < 0) {
	for (i = 1; i <= n4; ++i) {
	    wa[i + in * n5] = -wa[i + in * n5];
	}
    }
    goto L23072;
 L23074:
    goto L23055;
 L23057:
    goto L23052;
 L23054:
    if (kr == 1) {
	for (j = 1; j <= p; ++j) {
	    d = fabs(wa[n1 + j * n5]);
	    if (d <= tol || 2.0 - d <= tol) {
		ift = 1;
		wa[n2 + (pp + 1) * n5] = 0;
		goto L80;
	    }
	}
    }
 L80:
    kount = 0;
    sum = 0;
    if (!ci2) {
	for (i = 1; i <= kl - 1; ++i) {
	    k = (int) (wa[i + p4 * n5] * sgn(wa[i + p4 * n5]));
	    coeff[k] = wa[i + p1 * n5] * sgn(wa[i + p4 * n5]);
	}
    }
    for (i = 1; i <= p; ++i) {
	kd = fabs(wa[n4 + i * n5]) - p;
	dsol[kd + lsol * n] = wa[n1 + i * n5] / 2.0 + 1.0;
	if (wa[n4 + i * n5] < 0) {
	    dsol[kd + lsol * n] = 1. - dsol[kd + lsol * n];
	}
	if (!ci2) {
	    sum += coeff[i] * wa[n5 + i * n5];
	    sol[i + 3 + lsol * sdim] = coeff[i];
	    h[i + lsol * pp] = kd;
	}
    }
    for (i = kl; i <= n; ++i) {
	kd = fabs(wa[i + p4 * n5]) - p;
	if (wa[i + p4 * n5] < 0) {
	    dsol[kd + lsol * n] = 0;
	}
	if (wa[i + p4 * n5] > 0) {
	    dsol[kd + lsol * n] = 1.;
	}
    }
    if (!ci2) {
	sol[lsol * sdim + 1] = smax;
	sol[lsol * sdim + 2] = sum;
	sum = 0;
	for (j = kl; j <= n; ++j) {
	    d = wa[j + p1 * n5] * sgn(wa[j + p4 * n5]);
	    sum += d * (smax + .5 * (sgn(d) - 1.0));
	}
	sol[lsol * sdim + 3] = sum;
	for (i = 1; i <= n; ++i) {
	    dsol[i + (lsol + 1) * n] = dsol[i + lsol * n];
	}
    }
    if (ci2) {
	a1 = 0;
	for (i = 1; i <= n; ++i) {
	    a1 += x[i + idxcf * n] * (dsol[i + lsol * n] + tau - 1);
	}
	tn = a1 / sqrt(qn[idxcf] * tau * (1 - tau));
	if (fabs(tn) < cutoff) {
	    smax = (lup)? big : -big;
	    for (i = 1; i <= kl - 1; ++i) {
		k = (int) (wa[i + p4 * n5] * sgn(wa[i + p4 * n5]));
		sol[k + sdim] = wa[i + p2 * n5] * sgn(wa[i + p4 * n5]);
		sol[k + (sdim << 1)] = wa[i + p3 * n5] * 
		    sgn(wa[i + p4 * n5]) / tnew;
	    }
	    for (i = kl; i <= n; ++i) {
		a1 = b1 = 0;
		k = (int) (wa[i + p4 * n5] * sgn(wa[i + p4 * n5]) - p);
		l = 1;
		for (j = 1; j <= p; ++j) {
		    if (j == idxcf) {
			++l;
		    }
		    a1 += x[k + l * n] * sol[j + sdim];
		    b1 += x[k + l * n] * sol[j + (sdim << 1)];
		    ++l;
		}
		tnt = (y[k] - a1) / (x[k + idxcf * n] - b1);
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
	    told = tnew;
	    if (lup) {
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
	    for (i = 1; i <= n; ++i) {
		wa[i + p3 * n5] = wa[i + p3 * n5] / told * tnew;
		wa[i + p1 * n5] = wa[i + p2 * n5] - wa[i + p3 * n5];
	    }
	    for (j = kr; j <= p3; ++j) {
		d = wa[out + j * n5];
		wa[n1 + j * n5] = wa[n1 + j * n5] - d - d;
		wa[n2 + j * n5] = wa[n2 + j * n5] - d - d;
		wa[out + j * n5] = -d;
	    }
	    wa[out + p4 * n5] = -wa[out + p4 * n5];
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
	for (i = 1; i <= n; ++i) {
	    s[i] = 0;
	}
	for (j = 1; j <= p; ++j) {
	    coeff[j] = 0;
	}
	smax = 2.0;
	for (j = 1; j <= p; ++j) {
	    b1 = wa[n3 + j * n5];
	    a1 = (-2.0 - wa[n2 + j * n5]) / b1;
	    b1 = -wa[n2 + j * n5] / b1;
	    if (a1 >= tau) {
		if (a1 < smax) {
		    smax = a1;
		    dif = (b1 - a1) / 2.0;
		}
	    }
	    if (b1 > tau && b1 < smax) {
		smax = b1;
		dif = (b1 - a1) / 2.0;
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
    wa[n2 + (pp + 1) * n5] = 2.0;
    ift = 2;
    goto L50;
 L40:
    if (lsol > 2) {
	sol[sdim + 1] = 0;
	sol[sdim + 2] = 0;
	sol[sdim + 3] = 0;
	sol[lsol * sdim + 1] = 1;
	sol[lsol * sdim + 2] = 0;
	sol[lsol * sdim + 3] = 0;
	for (i = 1; i <= n; ++i) {
	    dsol[i + n] = 1.;
	    dsol[i + lsol * n] = 0;
	    dsol[i + (lsol + 1) * n] = 0;
	}
    }
    l = kl - 1;
    for (i = 1; i <= l; ++i) {
	if (wa[i + p1 * n5] < 0) {
	    for (j = kr; j <= p4; ++j) {
		wa[i + j * n5] = -wa[i + j * n5];
	    }
	}
    }
 L50:
    sum = 0;
    if (!ci2) {
	for (i = kl; i <= n; ++i) {
	    k = (int) (wa[i + p4 * n5] * sgn(wa[i + p4 * n5]));
	    d = wa[i + p1 * n5] * sgn(wa[i + p4 * n5]);
	    sum += d * sgn(d) * (.5 + sgn(d) * (tau - .5));
	    k -= p;
	    resid[k] = d;
	}
	wa[n2 + p2 * n5] = kount;
	wa[n1 + p2 * n5] = p1 - kr;
	wa[n1 + p1 * n5] = sum;
    }
    if (wa[n2 + (pp + 1) * n5] == 2.0) {
	goto L23018;
    }
    if (!ci1) {
	goto L23018;
    }
    if (ci2) {
	goto L23234;
    }
    ci2 = 1;
    p = pp - 1;
    p1 = p + 1;
    p2 = p + 2;
    p3 = p + 3;
    p4 = p + 4;
 L60:
    ++idxcf;
    if (idxcf <= pp) {
	goto L23236;
    }
    goto L23018;
 L23236:
 L70:
    if (!lup) {
	goto L23238;
    }
    tnew = coeff[idxcf] + tol;
    told = tnew;
    ci[(idxcf << 2) + 3] = coeff[idxcf];
    tnmat[(idxcf << 2) + 3] = 0;
    goto L23239;
 L23238:
    tnew = coeff[idxcf] - tol;
    told = tnew;
    ci[(idxcf << 2) + 2] = coeff[idxcf];
    tnmat[(idxcf << 2) + 2] = 0;
 L23239:
 L23234:
    goto L23016;
 L23018:
    for (i = 1; i <= n; ++i) {
	dsol[i + lsol * n] = dsol[i + (lsol + 1) * n];
    }

    return ift;
}

