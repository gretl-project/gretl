/* rqbr.{r,f,c}
   Original RATFOR code by Roger Koenker.
   Translated to C by f2c (version 20061008).
   Somewhat cleaned up by Allin Cottrell.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sgn(x) (((x) >= 0)? 1.0 : -1.0)

int rqbr_ (int n, int pp, double *x, double *y, double tau, 
	   double tol, double *coeff, double *resid, 
	   int *s, double *wa, double *wb, 
	   double *sol, double *dsol, int *h, 
	   double *qn, double cutoff, double *ci, double *tnmat,
	   double big, int rmax, int ci1,
	   void (*callback)(void))
{
    static double d, a1, b1;
    int i, j, k, l, jj;
    int n1, n2, n3, n4, p1, p2;
    int kd, kl = 0, in = 0, kr = 0;
    double tn, dmin, dmax, aux;
    double sum, tnt, told;
    char init, iend;
    double tnew, smax = 0;
    char lup, stage1, test = 0;
    int idxcf, kount, out = 0;
    int main_iters = 0;
    double pivot;
    int p = pp;
    int p3 = p + 3, p4 = p + 4, n5 = n + 5;
    int sdim, ci2 = 0;
    int ift = 0, rcount = 0;

    /* Parameter adjustments (ugh!) */
    --wb;
    --s;
    --resid;
    --y;
    --qn;
    --coeff;

    if (ci != NULL) {
	ci -= 5;
    }
    if (tnmat != NULL) {
	tnmat -= 5;
    }
 
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

 looptop:

    if (callback != NULL && (main_iters++ % 10 == 0)) {
	callback();
    }

    for (i = 1; i <= n; ++i) {
	k = 1;
	for (j = 1; j <= pp; ++j) {
	    if (k <= pp && j != idxcf) {
		wa[i + k * n5] = x[i + j * n];
		k++;
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
	    wa[n2 + j * n5] += aux * (1 - sgn(wa[i + p4 * n5]));
	    wa[n3 + j * n5] += aux * sgn(wa[i + p4 * n5]);
	}
	wa[n3 + j * n5] *= 2.0;
    }
    init = 0;
    if (!ci2) {
	/* compute the p column means */
	for (k = 1; k <= p; ++k) {
	    wa[n5 + k * n5] = 0;
	    for (i = 1; i <= n; ++i) {
		wa[n5 + k * n5] += x[i + k * n];
	    }
	    wa[n5 + k * n5] /= n;
	}
    }
    kount = 0;
 L23045:
    /* compute new marginal costs */
    for (j = 1; j <= p; ++j) {
	wa[n1 + j * n5] = wa[n2 + j * n5] + wa[n3 + j * n5] * tau;
    }
    if (!init) {
	/* stage 1: determine the vector to enter the basis */
	stage1 = 1;
	kr = kl = 1;
	goto L30;
    }
 L23052:
    /* stage 2: determine the vector to enter the basis */
    stage1 = 0;
 L23055:
    dmax = -big;
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
	jj = in * n5;
	for (i = 1; i <= n4; ++i) {
	    wa[i + jj] = -wa[i + jj];
	}
	wa[n1 + jj] += -2.0;
	wa[n2 + jj] += -2.0;
    }
 L23072:
    /* determine the vector to leave the basis */
    k = 0;
    for (i = kl; i <= n; ++i) {
	d = wa[i + in * n5];
	if (d > tol) {
	    k++;
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
	k--;
    }
    /* check for linear dependence in stage 1 */
    if (!test && stage1) {
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
	wa[n1 + j * n5] -= d + d;
	wa[n2 + j * n5] -= d + d;
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
    kr++;
    goto L20;
 L10:
    /* pivot on wa(out, in) */
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
    kount++;
    if (!stage1) {
	goto L23055;
    }
    /* interchange rows in stage 1 */
    kl++;
    for (j = kr; j <= p4; ++j) {
	d = wa[out + j * n5];
	wa[out + j * n5] = wa[kount + j * n5];
	wa[kount + j * n5] = d;
    }
 L20:
    if (kount + kr == p1) {
	goto L23052;
    }
 L30:
    dmax = -1;
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
 L23054:
    if (kr == 1) {
	for (j = 1; j <= p; ++j) {
	    d = fabs(wa[n1 + j * n5]);
	    if (d <= tol || 2.0 - d <= tol) {
		ift = 1;
		wa[n2 + (pp + 1) * n5] = 0;
		break;
	    }
	}
    }
    kount = 0;
    sum = 0;
    if (!ci2) {
	jj = p4 * n5;
	for (i = 1; i <= kl - 1; ++i) {
	    k = (int) (wa[i + jj] * sgn(wa[i + jj]));
	    coeff[k] = wa[i + p1 * n5] * sgn(wa[i + jj]);
	}
    }
    for (i = 1; i <= p; ++i) {
	kd = fabs(wa[n4 + i * n5]) - p;
	dsol[kd + n] = wa[n1 + i * n5] / 2.0 + 1.0;
	if (wa[n4 + i * n5] < 0) {
	    dsol[kd + n] = 1 - dsol[kd + n];
	}
	if (!ci2) {
	    sum += coeff[i] * wa[n5 + i * n5];
	    sol[i + 3 + sdim] = coeff[i];
	    h[i + pp] = kd;
	}
    }
    for (i = kl; i <= n; ++i) {
	kd = fabs(wa[i + p4 * n5]) - p;
	if (wa[i + p4 * n5] < 0) {
	    dsol[kd + n] = 0;
	}
	if (wa[i + p4 * n5] > 0) {
	    dsol[kd + n] = 1.0;
	}
    }
    if (!ci2) {
	sol[sdim + 1] = smax;
	sol[sdim + 2] = sum;
	sum = 0;
	for (j = kl; j <= n; ++j) {
	    d = wa[j + p1 * n5] * sgn(wa[j + p4 * n5]);
	    sum += d * (smax + .5 * (sgn(d) - 1.0));
	}
	sol[sdim + 3] = sum;
	for (i = 1; i <= n; ++i) {
	    dsol[i + 2 * n] = dsol[i + n];
	}
    } else {
	/* ci2: compute next theta */
	a1 = 0;
	for (i = 1; i <= n; ++i) {
	    a1 += x[i + idxcf * n] * (dsol[i + n] + tau - 1);
	}
	tn = a1 / sqrt(qn[idxcf] * tau * (1 - tau));
	if (fabs(tn) < cutoff) {
	    jj = p4 * n5;
	    smax = (lup)? big : -big;
	    for (i = 1; i <= kl - 1; ++i) {
		k = (int) (wa[i + jj] * sgn(wa[i + jj]));
		sol[k + sdim] = wa[i + p2 * n5] * sgn(wa[i + jj]);
		sol[k + (sdim << 1)] = wa[i + p3 * n5] * sgn(wa[i + jj]) / tnew;
	    }
	    for (i = kl; i <= n; ++i) {
		a1 = b1 = 0;
		k = (int) (wa[i + jj] * sgn(wa[i + jj]) - p);
		l = 1;
		for (j = 1; j <= p; ++j) {
		    if (j == idxcf) {
			l++;
		    }
		    a1 += x[k + l * n] * sol[j + sdim];
		    b1 += x[k + l * n] * sol[j + (sdim << 1)];
		    l++;
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
	    jj = idxcf << 2;
	    if (lup) {
		tnew = smax + tol;
		ci[jj + 3] = told - tol;
		tnmat[jj + 3] = tn;
		if (! (tnew < big - tol)) {
		    ci[jj + 3] = big;
		    ci[jj + 4] = big;
		    tnmat[jj + 3] = tn;
		    tnmat[jj + 4] = tn;
		    lup = 0;
		    goto L70;
		}
	    } else {
		tnew = smax - tol;
		ci[jj + 2] = told + tol;
		tnmat[jj + 2] = tn;
		if (! (tnew > -big + tol)) {
		    ci[jj + 2] = -big;
		    ci[jj + 1] = -big;
		    tnmat[jj + 2] = tn;
		    tnmat[jj + 1] = tn;
		    lup = 1;
		    goto L60;
		}
	    }
	    /* update the new marginal cost */
	    for (i = 1; i <= n; ++i) {
		wa[i + p3 * n5] = wa[i + p3 * n5] / told * tnew;
		wa[i + p1 * n5] = wa[i + p2 * n5] - wa[i + p3 * n5];
	    }
	    for (j = kr; j <= p3; ++j) {
		d = wa[out + j * n5];
		wa[n1 + j * n5] -= d + d;
		wa[n2 + j * n5] -= d + d;
		wa[out + j * n5] = -d;
	    }
	    wa[out + p4 * n5] = -wa[out + p4 * n5];
	    init = 1;
	} else if (lup) {
	    jj = idxcf << 2;
	    ci[jj + 4] = tnew - tol;
	    tnmat[jj + 4] = tn;
	    lup = 0;
	    goto L70;
	} else {
	    jj = idxcf << 2;
	    ci[jj + 1] = tnew + tol;
	    tnmat[jj + 1] = tn;
	    lup = 1;
	    goto L60;
	}
    }
    if (iend && !ci2) {
	goto L40;
    }
    goto L23045;
 L23047:
    wa[n2 + (pp + 1) * n5] = 2.0;
    ift = 2;
    goto L50;
 L40:
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
	jj = p4 * n5;
	for (i = kl; i <= n; ++i) {
	    k = (int) (wa[i + jj] * sgn(wa[i + jj]));
	    d = wa[i + p1 * n5] * sgn(wa[i + jj]);
	    sum += d * sgn(d) * (.5 + sgn(d) * (tau - .5));
	    k -= p;
	    resid[k] = d;
	}
	wa[n2 + p2 * n5] = kount;
	wa[n1 + p2 * n5] = p1 - kr;
	wa[n1 + p1 * n5] = sum;
    }
    if (!ci1 || wa[n2 + (pp + 1) * n5] == 2.0) {
	goto getout;
    }
    if (ci2) {
	goto looptop;
    }
    ci2 = 1;
    p = pp - 1;
    p1 = p + 1;
    p2 = p + 2;
    p3 = p + 3;
    p4 = p + 4;
 L60:
    if (++idxcf > pp) {
	goto getout;
    }
 L70:
    if (lup) {
	jj = idxcf << 2;
	tnew = coeff[idxcf] + tol;
	told = tnew;
	ci[jj + 3] = coeff[idxcf];
	tnmat[jj + 3] = 0;
    } else {
	jj = idxcf << 2;
	tnew = coeff[idxcf] - tol;
	told = tnew;
	ci[jj + 2] = coeff[idxcf];
	tnmat[jj + 2] = 0;
    }
    goto looptop;

 getout:
    /* restore the original value of dual when ci is true */
    for (i = 1; i <= n; ++i) {
	dsol[i + n] = dsol[i + 2 * n];
    }

    return ift;
}

