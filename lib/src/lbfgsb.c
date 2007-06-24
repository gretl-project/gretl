/* This code is ased on L-BFGS-B (version 2.1); the core
   is an f2c translation of routines.f

   See http://www.ece.northwestern.edu/~nocedal/lbfgsb.html

*/

#include "libgretl.h"
#include "nls_private.h"

#include <time.h>

#define abs(x) ((x) >= 0 ? (x) : -(x))
#ifndef min
# define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
# define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#define STPMIN 0.0
#define FTOL .001
#define GTOL .9
#define XTOL .1

static double ddot_(int n, double *dx, double *dy)
{
    double dtemp = 0;
    int i;
    
    for (i=1; i<=n; i++) {
	dtemp += dx[i] * dy[i];
    }

    return dtemp;
}

static void dscal_(int n, double da, double *dx)
{
    int i;

    for (i=1; i<=n; i++) {
	dx[i] *= da;
    }    
} 

static void dpofa_(double *a, int lda, int n, int *info)
{
    int j, k;
    double s, t;
    int jm1;

    for (j = 1; j <= n; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	for (k = 1; k <= jm1; ++k) {
	    t = a[k + j * lda] - ddot_(k - 1, &a[k * lda], &a[j * lda]);
	    t /= a[k + k * lda];
	    a[k + j * lda] = t;
	    s += t * t;
	}
L20:
	s = a[j + j * lda] - s;
	if (s <= 0.) {
	    return;
	}
	a[j + j * lda] = sqrt(s);
    }

    *info = 0;
}

static void dcopy_(int n, double *dx, double *dy)
{
    int i;

    for (i=1; i<=n; i++) {
	dy[i] = dx[i];
    }
}

static void daxpy_(int n, double da, double *dx, double *dy)
{
    int i;

    if (n <= 0 || da == 0) {
	return;
    }

    for (i=1; i<=n; i++) {
	dy[i] += da * dx[i];
    }
} 

static void dtrsl_(double *t, int ldt, int n, 
		   double *b, int job, int *info)
{
    int j, jj, K = 1;
    double temp;

    for (*info = 1; *info <= n; ++(*info)) {
	if (t[*info + *info * ldt] == 0.) {
	    return;
	}
    }

    *info = 0;

    if (job % 10 != 0) {
	K = 2;
    }
    if (job % 100 / 10 != 0) {
	K += 2;
    }

    if (K == 1) {
	b[1] /= t[ldt + 1];
	for (j = 2; j <= n; ++j) {
	    temp = -b[j - 1];
	    daxpy_(n - j + 1, temp, &t[j - 1 + (j - 1) * ldt], &b[j-1]);
	    b[j] /= t[j + j * ldt];
	}
    } else if (K == 2) {
	b[n] /= t[n + n * ldt];
	for (jj = 2; jj <= n; ++jj) {
	    j = n - jj + 1;
	    temp = -b[j + 1];
	    daxpy_(j, temp, &t[(j + 1) * ldt], b);
	    b[j] /= t[j + j * ldt];
	}
    } else if (K == 3) {
	b[n] /= t[n + n * ldt];
	for (jj = 2; jj <= n; ++jj) {
	    j = n - jj + 1;
	    b[j] -= ddot_(jj - 1, &t[j + j * ldt], &b[j]);
	    b[j] /= t[j + j * ldt];
	}
    } else if (K == 4) {
	b[1] /= t[ldt + 1];
	for (j = 2; j <= n; ++j) {
	    b[j] -= ddot_(j - 1, &t[j * ldt], b);
	    b[j] /= t[j + j * ldt];
	}
    }
} 

static void active_(int n, double *l, double *u, 
		    int *nbd, double *x, int *iwhere,  
		    int *prjctd, int *cnstnd, int *boxed)
{
    int i, nbdd;

    nbdd = 0;
    *prjctd = 0;
    *cnstnd = 0;
    *boxed = 1;

    for (i = 1; i <= n; ++i) {
	if (nbd[i] > 0) {
	    if (nbd[i] <= 2 && x[i] <= l[i]) {
		if (x[i] < l[i]) {
		    *prjctd = 1;
		    x[i] = l[i];
		}
		++nbdd;
	    } else if (nbd[i] >= 2 && x[i] >= u[i]) {
		if (x[i] > u[i]) {
		    *prjctd = 1;
		    x[i] = u[i];
		}
		++nbdd;
	    }
	}
    }

    for (i = 1; i <= n; ++i) {
	if (nbd[i] != 2) {
	    *boxed = 0;
	}
	if (nbd[i] == 0) {
	    iwhere[i] = -1;
	} else {
	    *cnstnd = 1;
	    if (nbd[i] == 2 && u[i] - l[i] <= 0.) {
		iwhere[i] = 3;
	    } else {
		iwhere[i] = 0;
	    }
	}
    }
}

static void bmv_(int m, double *sy, double *wt, int col, 
		 double *v, double *p, int *info)
{
    int i, k, icol;
    double sum;

    if (col == 0) {
	return;
    }

    p[col + 1] = v[col + 1];

    for (i = 2; i <= col; ++i) {
	icol = col + i;
	sum = 0.;
	for (k = 1; k <= i - 1; ++k) {
	    sum += sy[i + k * m] * v[k] / sy[k + k * m];
	}
	p[icol] = v[icol] + sum;
    }

    dtrsl_(wt, m, col, &p[col], 11, info);
    if (*info != 0) {
	return;
    }

    for (i = 1; i <= col; ++i) {
	p[i] = v[i] / sqrt(sy[i + i * m]);
    }

    dtrsl_(wt, m, col, &p[col], 1, info);
    if (*info != 0) {
	return;
    }
    
    for (i = 1; i <= col; ++i) {
	p[i] = -p[i] / sqrt(sy[i + i * m]);
    }

    for (i = 1; i <= col; ++i) {
	sum = 0.;
	for (k = i + 1; k <= col; ++k) {
	    sum += sy[k + i * m] * p[col + k] / sy[i + i * m];
	}
	p[i] += sum;
    }
}

static void hpsolb_(int n, double *t, int *iorder, int iheap)
{
    int i, j, k;
    double out, ddum;
    int indxin, indxou;

    if (iheap == 0) {
	for (k = 2; k <= n; ++k) {
	    ddum = t[k];
	    indxin = iorder[k];
	    i = k;
L10:
	    if (i > 1) {
		j = i / 2;
		if (ddum < t[j]) {
		    t[i] = t[j];
		    iorder[i] = iorder[j];
		    i = j;
		    goto L10;
		}
	    }
	    t[i] = ddum;
	    iorder[i] = indxin;
	}
    }

    if (n > 1) {
	i = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[n];
	indxin = iorder[n];
L30:
	j = i + i;
	if (j <= n - 1) {
	    if (t[j + 1] < t[j]) {
		++j;
	    }
	    if (t[j] < ddum) {
		t[i] = t[j];
		iorder[i] = iorder[j];
		i = j;
		goto L30;
	    }
	}
	t[i] = ddum;
	iorder[i] = indxin;
	t[n] = out;
	iorder[n] = indxou;
    }
}

static int 
cauchy_(int n, double *x, double *l, double *u, int *nbd, 
	double *g, int *iorder, int *iwhere, double *t, 
	double *d, double *xcp, int m, 
	double *wy, double *ws, double *sy, double *wt, 
	double theta, int col, int head, double *p, 
	double *c, double *wbp, double *v, int *nint, 
	double *sg, double *yg, double sbgnrm, int *info, 
	double epsmch)
{
    double d1;

    int i, j;
    double f1, f2, dt, tj, tj0;
    double tu = 0.0, tl = 0.0;
    int ibp;
    double dtm;
    double wmc, wmp, wmw;
    int col2;
    double dibp;
    int iter;
    double zibp, tsum, dibp2;
    int bnded;
    double neggi;
    int nfree;
    double bkmin;
    int nleft;
    double f2_org__;
    int nbreak, ibkmin;
    int pointr;
    int xlower, xupper;

    if (sbgnrm <= 0.) {
	dcopy_(n, x, xcp);
	return 0;
    }

    bnded = 1;
    nfree = n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = col << 1;
    f1 = 0.;

    for (i = 1; i <= col2; ++i) {
	p[i] = 0.;
    }

    for (i = 1; i <= n; ++i) {
	neggi = -g[i];
	if (iwhere[i] != 3 && iwhere[i] != -1) {
	    if (nbd[i] <= 2) {
		tl = x[i] - l[i];
	    }
	    if (nbd[i] >= 2) {
		tu = u[i] - x[i];
	    }
	    xlower = nbd[i] <= 2 && tl <= 0.;
	    xupper = nbd[i] >= 2 && tu <= 0.;
	    iwhere[i] = 0;
	    if (xlower) {
		if (neggi <= 0.) {
		    iwhere[i] = 1;
		}
	    } else if (xupper) {
		if (neggi >= 0.) {
		    iwhere[i] = 2;
		}
	    } else {
		if (abs(neggi) <= 0.) {
		    iwhere[i] = -3;
		}
	    }
	}
	pointr = head;
	if (iwhere[i] != 0 && iwhere[i] != -1) {
	    d[i] = 0.;
	} else {
	    d[i] = neggi;
	    f1 -= neggi * neggi;
	    for (j = 1; j <= col; ++j) {
		p[j] += wy[i + pointr * n] * neggi;
		p[col + j] += ws[i + pointr * n] * neggi;
		pointr = pointr % m + 1;
	    }
	    if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.) {
		++nbreak;
		iorder[nbreak] = i;
		t[nbreak] = tl / (-neggi);
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else if (nbd[i] >= 2 && neggi > 0.) {
		++nbreak;
		iorder[nbreak] = i;
		t[nbreak] = tu / neggi;
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else {
		--nfree;
		iorder[nfree] = i;
		if (abs(neggi) > 0.) {
		    bnded = 0;
		}
	    }
	}
    }

    if (theta != 1.) {
	dscal_(col, theta, &p[col]);
    }

    dcopy_(n, x, xcp);
    if (nbreak == 0 && nfree == n + 1) {
	return 0;
    }

    for (j = 1; j <= col2; ++j) {
	c[j] = 0.;
    }

    f2 = -theta * f1;
    f2_org__ = f2;
    if (col > 0) {
	bmv_(m, sy, wt, col, p, v, info);
	if (*info != 0) {
	    return 0;
	}
	f2 -= ddot_(col2, v, p);
    }
    dtm = -f1 / f2;
    tsum = 0.;
    *nint = 1;

    if (nbreak == 0) {
	goto L888;
    }
    nleft = nbreak;
    iter = 1;
    tj = 0.;

 L777:

    tj0 = tj;

    if (iter == 1) {
	tj = bkmin;
	ibp = iorder[ibkmin];
    } else {
	if (iter == 2) {
	    if (ibkmin != nbreak) {
		t[ibkmin] = t[nbreak];
		iorder[ibkmin] = iorder[nbreak];
	    }
	}
	hpsolb_(nleft, t, iorder, iter - 2);
	tj = t[nleft];
	ibp = iorder[nleft];
    }

    dt = tj - tj0;

    if (dtm < dt) {
	goto L888;
    }

    tsum += dt;
    --nleft;
    ++iter;
    dibp = d[ibp];
    d[ibp] = 0.;
    if (dibp > 0.) {
	zibp = u[ibp] - x[ibp];
	xcp[ibp] = u[ibp];
	iwhere[ibp] = 2;
    } else {
	zibp = l[ibp] - x[ibp];
	xcp[ibp] = l[ibp];
	iwhere[ibp] = 1;
    }
    if (nleft == 0 && nbreak == n) {
	dtm = dt;
	goto L999;
    }

    ++(*nint);
    d1 = dibp;
    dibp2 = d1 * d1;
    f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
    f2 -= theta * dibp2;

    if (col > 0) {
	daxpy_(col2, dt, p, c);
	pointr = head;
	for (j = 1; j <= col; ++j) {
	    wbp[j] = wy[ibp + pointr * n];
	    wbp[col + j] = theta * ws[ibp + pointr * n];
	    pointr = pointr % m + 1;
	}
	bmv_(m, sy, wt, col, wbp, v, info);
	if (*info != 0) {
	    return 0;
	}
	wmc = ddot_(col2, c, v);
	wmp = ddot_(col2, p, v);
	wmw = ddot_(col2, wbp, v);
	d1 = -dibp;
	daxpy_(col2, d1, wbp, p);
	f1 += dibp * wmc;
	f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }

    d1 = epsmch * f2_org__;
    f2 = max(d1,f2);
    if (nleft > 0) {
	dtm = -f1 / f2;
	goto L777;
    } else if (bnded) {
	f1 = 0.;
	f2 = 0.;
	dtm = 0.;
    } else {
	dtm = -f1 / f2;
    }

 L888:
    if (dtm <= 0.) {
	dtm = 0.;
    }
    tsum += dtm;

    daxpy_(n, tsum, d, xcp);

 L999:
    if (col > 0) {
	daxpy_(col2, dtm, p, c);
    }
    return 0;
}

static void 
cmprlb_(int n, int m, double *x, 
	double *g, double *ws, double *wy, double *sy, 
	double *wt, double *z, double *r, double *wa, 
	int *index, double theta, int col, int head, 
	int nfree, int cnstnd, int *info)
{
    int i, j, k;
    double a1, a2;
    int pointr;

    if (!cnstnd && col > 0) {
	for (i = 1; i <= n; ++i) {
	    r[i] = -g[i];
	}
    } else {
	for (i = 1; i <= nfree; ++i) {
	    k = index[i];
	    r[i] = -(theta) * (z[k] - x[k]) - g[k];
	}
	bmv_(m, sy, wt, col, &wa[(m << 1)], wa, info);
	if (*info != 0) {
	    *info = -8;
	    return;
	}
	pointr = head;
	for (j = 1; j <= col; ++j) {
	    a1 = wa[j];
	    a2 = theta * wa[col + j];
	    for (i = 1; i <= nfree; ++i) {
		k = index[i];
		r[i] = r[i] + wy[k + pointr * n] * a1 + 
		    ws[k + pointr * n] * a2;
	    }
	    pointr = pointr % m + 1;
	}
    }
} 

static void errclb_(int n, int m, double factr, 
		    double *l, double *u, int *nbd, 
		    char *task, int *info, int *k)
{
    int i;

    if (n <= 0) {
	strcpy(task, "ERROR: N .LE. 0");
    }
    if (m <= 0) {
	strcpy(task, "ERROR: M .LE. 0");
    }
    if (factr < 0.) {
	strcpy(task, "ERROR: FACTR .LT. 0");
    }

    for (i = 1; i <= n; ++i) {
	if (nbd[i] < 0 || nbd[i] > 3) {
	    strcpy(task, "ERROR: INVALID NBD");
	    *info = -6;
	    *k = i;
	} else if (nbd[i] == 2) {
	    if (l[i] > u[i]) {
		strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
		*info = -7;
		*k = i;
	    }
	}
    }
} 

static void 
formk_(int n, int nsub, int *ind, int nenter, int ileave, 
       int *indx2, int iupdat, int updatd, double *wn, double *wn1, 
       int m, double *ws, double *wy, double *sy, double theta, 
       int col, int head, int *info)
{
    int i, k, k1, m2, is, js, iy, jy, is1, js1, col2, dend, pend;
    double temp1, temp2, temp3, temp4;
    int upcl, ipntr, jpntr, dbegin, pbegin;

    m2 = m << 1;

    if (updatd) {
	if (iupdat > m) {
	    for (jy = 1; jy <= m - 1; ++jy) {
		js = m + jy;
		dcopy_(m - jy, &wn1[jy + (jy + 1) * m2], 
		       &wn1[jy + jy * m2 - 1]);
		dcopy_(m - jy, &wn1[js + (js + 1) * m2], 
		       &wn1[js + js * m2 - 1]);
		dcopy_(m - 1, &wn1[m + 1 + (jy + 1) * m2], 
		       &wn1[m + jy * m2]);
	    }
	}
	pbegin = 1;
	pend = nsub;
	dbegin = nsub + 1;
	dend = n;
	iy = col;
	is = m + col;
	ipntr = head + col - 1;
	if (ipntr > m) {
	    ipntr -= m;
	}
	jpntr = head;
	for (jy = 1; jy <= col; ++jy) {
	    js = m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    for (k = pbegin; k <= pend; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    for (k = dbegin; k <= dend; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * n] * ws[k1 + jpntr * n];
		temp3 += ws[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    wn1[iy + jy * m2] = temp1;
	    wn1[is + js * m2] = temp2;
	    wn1[is + jy * m2] = temp3;
	    jpntr = jpntr % m + 1;
	}

	jy = col;
	jpntr = head + col - 1;
	if (jpntr > m) {
	    jpntr -= m;
	}
	ipntr = head;
	for (i = 1; i <= col; ++i) {
	    is = m + i;
	    temp3 = 0.;
	    for (k = pbegin; k <= pend; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    ipntr = ipntr % m + 1;
	    wn1[is + jy * m2] = temp3;
	}
	upcl = col - 1;
    } else {
	upcl = col;
    }

    ipntr = head;
    for (iy = 1; iy <= upcl; ++iy) {
	is = m + iy;
	jpntr = head;
	for (jy = 1; jy <= iy; ++jy) {
	    js = m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    temp4 = 0.;
	    for (k = 1; k <= nenter; ++k) {
		k1 = indx2[k];
		temp1 += wy[k1 + ipntr * n] * wy[k1 + jpntr * n];
		temp2 += ws[k1 + ipntr * n] * ws[k1 + jpntr * n];
	    }
	    for (k = ileave; k <= n; ++k) {
		k1 = indx2[k];
		temp3 += wy[k1 + ipntr * n] * wy[k1 + jpntr * n];
		temp4 += ws[k1 + ipntr * n] * ws[k1 + jpntr * n];
	    }
	    wn1[iy + jy * m2] = wn1[iy + jy * m2] + temp1 - temp3;
	    wn1[is + js * m2] = wn1[is + js * m2] - temp2 + temp4;
	    jpntr = jpntr % m + 1;
	}
	ipntr = ipntr % m + 1;
    }

    ipntr = head;
    for (is = m + 1; is <= m + upcl; ++is) {
	jpntr = head;
	for (jy = 1; jy <= upcl; ++jy) {
	    temp1 = 0.;
	    temp3 = 0.;
	    for (k = 1; k <= nenter; ++k) {
		k1 = indx2[k];
		temp1 += ws[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    for (k = ileave; k <= n; ++k) {
		k1 = indx2[k];
		temp3 += ws[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    if (is <= jy + m) {
		wn1[is + jy * m2] = wn1[is + jy * m2] + temp1 - temp3;
	    } else {
		wn1[is + jy * m2] = wn1[is + jy * m2] - temp1 + temp3;
	    }
	    jpntr = jpntr % m + 1;
	}
	ipntr = ipntr % m + 1;
    }

    for (iy = 1; iy <= col; ++iy) {
	is = col + iy;
	is1 = m + iy;
	for (jy = 1; jy <= iy; ++jy) {
	    js = col + jy;
	    js1 = m + jy;
	    wn[jy + iy * m2] = wn1[iy + jy * m2] / theta;
	    wn[js + is * m2] = wn1[is1 + js1 * m2] * theta;
	}
	for (jy = 1; jy <= iy - 1; ++jy) {
	    wn[jy + is * m2] = -wn1[is1 + jy * m2];
	}
	for (jy = iy; jy <= col; ++jy) {
	    wn[jy + is * m2] = wn1[is1 + jy * m2];
	}
	wn[iy + iy * m2] += sy[iy + iy * m];
    }

    dpofa_(wn, m2, col, info);
    if (*info != 0) {
	*info = -1;
	return;
    }

    col2 = col << 1;

    for (js = col + 1; js <= col2; ++js) {
	dtrsl_(wn, m2, col, &wn[js * m2], 11, info);
    }

    for (is = col + 1; is <= col2; ++is) {
	for (js = is; js <= col2; ++js) {
	    wn[is + js * m2] += ddot_(col, &wn[is * m2], &wn[js * m2]);
	}
    }

    dpofa_(&wn[col + col * m2], m2, col, info);
    if (*info != 0) {
	*info = -2;
    }
} 

static void formt_(int m, double *wt, double *sy, 
		   double *ss, int col, double theta, 
		   int *info)
{
    int i, j, k, k1;
    double ddum;

    for (j = 1; j <= col; ++j) {
	wt[j * m + 1] = theta * ss[j * m + 1];
    }

    for (i = 2; i <= col; ++i) {
	for (j = i; j <= col; ++j) {
	    k1 = min(i,j) - 1;
	    ddum = 0.;
	    for (k = 1; k <= k1; ++k) {
		ddum += sy[i + k * m] * sy[j + k * m] / sy[k + k * m];
	    }
	    wt[i + j * m] = ddum + theta * ss[i + j * m];
	}
    }

    dpofa_(wt, m, col, info);
    if (*info != 0) {
	*info = -3;
    }
} 

static void
freev_(int n, int *nfree, int *index, 
       int *nenter, int *ileave, int *indx2, int *iwhere, 
       int *wrk, int updatd, int cnstnd, int iter)
{
    int i, k, iact;

    *nenter = 0;
    *ileave = n + 1;

    if (iter > 0 && cnstnd) {
	for (i = 1; i <= *nfree; ++i) {
	    k = index[i];
	    if (iwhere[k] > 0) {
		--(*ileave);
		indx2[*ileave] = k;
	    }
	}
	for (i = *nfree + 1; i <= n; ++i) {
	    k = index[i];
	    if (iwhere[k] <= 0) {
		++(*nenter);
		indx2[*nenter] = k;
	    }
	}
    }

    *wrk = *ileave < n + 1 || *nenter > 0 || updatd;

    *nfree = 0;
    iact = n + 1;
    for (i = 1; i <= n; ++i) {
	if (iwhere[i] <= 0) {
	    ++(*nfree);
	    index[*nfree] = i;
	} else {
	    --iact;
	    index[iact] = i;
	}
    }
} 

static int dcstep_(double *stx, double *fx, double *dx, 
		   double *sty, double *fy, double *dy, 
		   double *stp, double fp, double dp, 
		   int *brackt, double stpmin, double stpmax)
{
    double d1, d2, d3;
    double p, q, r, s, sgnd, stpc, stpf, stpq, gamma, theta;

    sgnd = dp * (*dx / abs(*dx));

    if (fp > *fx) {
	theta = (*fx - fp) * 3. / (*stp - *stx) + *dx + dp;
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(dp);
	s = max(d1, d2);
	d1 = theta / s;
	gamma = s * sqrt(d1 * d1 - *dx / s * (dp / s));
	if (*stp < *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dx + theta;
	q = gamma - *dx + gamma + dp;
	r = p / q;
	stpc = *stx + r * (*stp - *stx);
	stpq = *stx + *dx / ((*fx - fp) / (*stp - *stx) + *dx) / 2. * (*stp - *stx);
	if ((d1 = stpc - *stx, abs(d1)) < (d2 = stpq - *stx, abs(d2))) {
	    stpf = stpc;
	} else {
	    stpf = stpc + (stpq - stpc) / 2.;
	}
	*brackt = 1;
    } else if (sgnd < 0.) {
	theta = (*fx - fp) * 3. / (*stp - *stx) + *dx + dp;
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(dp);
	s = max(d1,d2);
	d1 = theta / s;
	gamma = s * sqrt(d1 * d1 - *dx / s * (dp / s));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - dp + theta;
	q = gamma - dp + gamma + *dx;
	r = p / q;
	stpc = *stp + r * (*stx - *stp);
	stpq = *stp + dp / (dp - *dx) * (*stx - *stp);
	if ((d1 = stpc - *stp, abs(d1)) > (d2 = stpq - *stp, abs(d2))) {
	    stpf = stpc;
	} else {
	    stpf = stpq;
	}
	*brackt = 1;
    } else if (abs(dp) < abs(*dx)) {
	theta = (*fx - fp) * 3. / (*stp - *stx) + *dx + dp;

	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(dp);
	s = max(d1,d2);

	d3 = theta / s;
	d1 = 0., d2 = d3 * d3 - *dx / s * (dp / s);
	gamma = s * sqrt((max(d1,d2)));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - dp + theta;
	q = gamma + (*dx - dp) + gamma;
	r = p / q;
	if (r < 0. && gamma != 0.) {
	    stpc = *stp + r * (*stx - *stp);
	} else if (*stp > *stx) {
	    stpc = stpmax;
	} else {
	    stpc = stpmin;
	}

	stpq = *stp + dp / (dp - *dx) * (*stx - *stp);

	if (*brackt) {
	    if ((d1 = stpc - *stp, abs(d1)) < (d2 = stpq - *stp, abs(d2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    if (*stp > *stx) {
		d1 = *stp + (*sty - *stp) * .66;
		stpf = min(d1,stpf);
	    } else {
		d1 = *stp + (*sty - *stp) * .66;
		stpf = max(d1,stpf);
	    }
	} else {
	    if ((d1 = stpc - *stp, abs(d1)) > (d2 = stpq - *stp, abs(d2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    stpf = min(stpmax, stpf);
	    stpf = max(stpmin, stpf);
	}
    } else {
	if (*brackt) {
	    theta = (fp - *fy) * 3. / (*sty - *stp) + *dy + dp;
	    d1 = abs(theta), d2 = abs(*dy), d1 = max(d1,d2), d2 = abs(dp);
	    s = max(d1,d2);
	    d1 = theta / s;
	    gamma = s * sqrt(d1 * d1 - *dy / s * (dp / s));
	    if (*stp > *sty) {
		gamma = -gamma;
	    }
	    p = gamma - dp + theta;
	    q = gamma - dp + gamma + *dy;
	    r = p / q;
	    stpc = *stp + r * (*sty - *stp);
	    stpf = stpc;
	} else if (*stp > *stx) {
	    stpf = stpmax;
	} else {
	    stpf = stpmin;
	}
    }

    if (fp > *fx) {
	*sty = *stp;
	*fy = fp;
	*dy = dp;
    } else {
	if (sgnd < 0.) {
	    *sty = *stx;
	    *fy = *fx;
	    *dy = *dx;
	}
	*stx = *stp;
	*fx = fp;
	*dx = dp;
    }

    *stp = stpf;

    return 0;
} 

static int dcsrch_(double f, double g, double *stp, 
		   double ftol, double gtol, double xtol, 
		   double stpmin, double stpmax, 
		   char *task, int *isave, double *dsave)
{
    double d1;

    double fm, gm, fx, fy, gx, gy, fxm, fym, gxm, gym, stx, sty;
    double finit, ginit, width, ftest, gtest, stmin, stmax, width1;
    int stage, brackt;

    if (!strcmp(task, "START")) {
	if (*stp < stpmin) {
	    strcpy(task, "ERROR: STP .LT. STPMIN");
	}
	if (*stp > stpmax) {
	    strcpy(task, "ERROR: STP .GT. STPMAX");
	}
	if (g >= 0.) {
	    strcpy(task, "ERROR: INITIAL G .GE. ZERO");
	}
	if (ftol < 0.) {
	    strcpy(task, "ERROR: FTOL .LT. ZERO");
	}
	if (gtol < 0.) {
	    strcpy(task, "ERROR: GTOL .LT. ZERO");
	}
	if (xtol < 0.) {
	    strcpy(task, "ERROR: XTOL .LT. ZERO");
	}
	if (stpmin < 0.) {
	    strcpy(task, "ERROR: STPMIN .LT. ZERO");
	}
	if (stpmax < stpmin) {
	    strcpy(task, "ERROR: STPMAX .LT. STPMIN");
	}
	if (!strncmp(task, "ERROR", 5)) {
	    return 1;
	}

	brackt = 0;
	stage = 1;
	finit = f;
	ginit = g;
	gtest = ftol * ginit;
	width = stpmax - stpmin;
	width1 = width / .5;

	stx = 0.;
	fx = finit;
	gx = ginit;
	sty = 0.;
	fy = finit;
	gy = ginit;
	stmin = 0.;
	stmax = *stp + *stp * 4.;
	strcpy(task, "FG");
	goto L1000;
    } else {
	if (isave[1] == 1) {
	    brackt = 1;
	} else {
	    brackt = 0;
	}
	stage = isave[2];
	ginit = dsave[1];
	gtest = dsave[2];
	gx = dsave[3];
	gy = dsave[4];
	finit = dsave[5];
	fx = dsave[6];
	fy = dsave[7];
	stx = dsave[8];
	sty = dsave[9];
	stmin = dsave[10];
	stmax = dsave[11];
	width = dsave[12];
	width1 = dsave[13];
    }

    ftest = finit + *stp * gtest;
    if (stage == 1 && f <= ftest && g >= 0.) {
	stage = 2;
    }

    if (brackt && (*stp <= stmin || *stp >= stmax)) {
	strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
    }
    if (brackt && stmax - stmin <= xtol * stmax) {
	strcpy(task, "WARNING: XTOL TEST SATISFIED");
    }
    if (*stp == stpmax && f <= ftest && g <= gtest) {
	strcpy(task, "WARNING: STP = STPMAX");
    }
    if (*stp == stpmin && (f > ftest || g >= gtest)) {
	strcpy(task, "WARNING: STP = STPMIN");
    }

    if (f <= ftest && abs(g) <= gtol * (-ginit)) {
	strcpy(task, "CONVERGENCE");
    }

    if (!strncmp(task, "WARN", 4) || !strncmp(task, "CONV", 4)) {
	goto L1000;
    }

     if (stage == 1 && f <= fx && f > ftest) {
	fm = f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, fm, gm, &brackt,
		stmin, stmax);
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, 
		stmin, stmax);
    }

    if (brackt) {
	if ((d1 = sty - stx, abs(d1)) >= width1 * .66) {
	    *stp = stx + (sty - stx) * .5;
	}
	width1 = width;
	width = (d1 = sty - stx, abs(d1));
    }

    if (brackt) {
	stmin = min(stx,sty);
	stmax = max(stx,sty);
    } else {
	stmin = *stp + (*stp - stx) * 1.1;
	stmax = *stp + (*stp - stx) * 4.;
    }

    *stp = max(*stp, stpmin);
    *stp = min(*stp, stpmax);

    if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin 
	    <= xtol * stmax) {
	*stp = stx;
    }

    strcpy(task, "FG");

L1000:

    if (brackt) {
	isave[1] = 1;
    } else {
	isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;

    return 0;
} 

static void lnsrlb_(int n, double *l, double *u, 
		    int *nbd, double *x, double f, double *fold, 
		    double *gd, double *gdold, double *g, double *d, 
		    double *r, double *t, double *z, double *stp, 
		    double *dnorm, double *dtd, double *xstep, double *stpmx, 
		    int iter, int *ifun, int *iback, int *nfgv, 
		    int *info, char *task, int boxed, int cnstnd, 
		    char *csave, int *isave, double *dsave)
{
    double a1, a2, d1; 
    int i;

    if (!strncmp(task, "FG_LN", 5)) {
	goto L556;
    }

    *dtd = ddot_(n, d, d);
    *dnorm = sqrt(*dtd);

    *stpmx = 1.0e10;

    if (cnstnd) {
	if (iter == 0) {
	    *stpmx = 1.;
	} else {
	    for (i = 1; i <= n; ++i) {
		a1 = d[i];
		if (nbd[i] != 0) {
		    if (a1 < 0. && nbd[i] <= 2) {
			a2 = l[i] - x[i];
			if (a2 >= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx < a2) {
			    *stpmx = a2 / a1;
			}
		    } else if (a1 > 0. && nbd[i] >= 2) {
			a2 = u[i] - x[i];
			if (a2 <= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx > a2) {
			    *stpmx = a2 / a1;
			}
		    }
		}
	    }
	}
    }

    if (iter == 0 && !boxed) {
	d1 = 1. / *dnorm;
	*stp = min(d1, *stpmx);
    } else {
	*stp = 1.;
    }

    dcopy_(n, x, t);
    dcopy_(n, g, r);
    *fold = f;
    *ifun = 0;
    *iback = 0;
    strcpy(csave, "START");

L556:
    *gd = ddot_(n, g, d);
    if (*ifun == 0) {
	*gdold = *gd;
	if (*gd >= 0.) {
	    *info = -4;
	    return;
	}
    }

    dcsrch_(f, *gd, stp, FTOL, GTOL, XTOL, STPMIN, *stpmx, csave, 
	    isave, dsave);

    *xstep = *stp * *dnorm;

    if (strncmp(csave, "CONV", 4) && strncmp(csave, "WARN", 4)) {
	strcpy(task, "FG_LNSRCH");
	++(*ifun);
	++(*nfgv);
	*iback = *ifun - 1;
	if (*stp == 1.) {
	    dcopy_(n, z, x);
	} else {
	    for (i = 1; i <= n; ++i) {
		x[i] = *stp * d[i] + t[i];
	    }
	}
    } else {
	strcpy(task, "NEW_X");
    }
} 

static int matupd_(int n, int m, double *ws, 
		   double *wy, double *sy, double *ss, double *d, 
		   double *r, int *itail, int iupdat, int *col, 
		   int *head, double *theta, double rr, double dr, 
		   double *stp, double *dtd)
{
    int j, jmax;
    int pointr;

    if (iupdat <= m) {
	*col = iupdat;
	*itail = (*head + iupdat - 2) % m + 1;
    } else {
	*itail = *itail % m + 1;
	*head = *head % m + 1;
    }

    dcopy_(n, d, &ws[*itail * n]);
    dcopy_(n, r, &wy[*itail * n]);

    *theta = rr / dr;

    if (iupdat > m) {
	jmax = *col - 1;
	for (j = 1; j <= jmax; ++j) {
	    dcopy_(j, &ss[(j + 1) * m + 1], &ss[j * m]);
	    dcopy_(*col - j, &sy[j + (j + 1) * m], &sy[j + j * m - 1]);
	}
    }

    pointr = *head;
    jmax = *col - 1;
    for (j = 1; j <= jmax; ++j) {
	sy[*col + j * m] = ddot_(n, d, &wy[pointr * n]);
	ss[j + *col * m] = ddot_(n, &ws[pointr * n], d);
	pointr = pointr % m + 1;
    }
    if (*stp == 1.) {
	ss[*col + *col * m] = *dtd;
    } else {
	ss[*col + *col * m] = *stp * *stp * *dtd;
    }

    sy[*col + *col * m] = dr;

    return 0;
} 

static void projgr_(int n, double *l, double *u, 
		    int *nbd, double *x, double *g, 
		    double *sbgnrm)
{
    double gi, d1, d2;
    int i;

    *sbgnrm = 0.;

    for (i = 1; i <= n; ++i) {
	gi = g[i];

	if (nbd[i] != 0) {
	    if (gi < 0.) {
		if (nbd[i] >= 2) {
		    d1 = x[i] - u[i];
		    gi = max(d1, gi);
		}
	    } else {
		if (nbd[i] <= 2) {
		    d1 = x[i] - l[i];
		    gi = min(d1, gi);
		}
	    }
	}

	d1 = *sbgnrm, d2 = abs(gi);
	*sbgnrm = max(d1, d2);
    }
} 

static int 
subsm_(int n, int m, int nsub, int *ind, 
       double *l, double *u, int *nbd, double *x, 
       double *d, double *ws, double *wy, double theta, 
       int col, int head, int *iword, double *wv, 
       double *wn, int *info)
{
    int i, j, k, m2;
    int js, jy, col2, ibd = 0;
    double dk, temp1, temp2, alpha;
    int pointr;

    if (nsub <= 0) {
	return 0;
    }

    pointr = head;
    for (i = 1; i <= col; ++i) {
	temp1 = 0.;
	temp2 = 0.;
	for (j = 1; j <= nsub; ++j) {
	    k = ind[j];
	    temp1 += wy[k + pointr * n] * d[j];
	    temp2 += ws[k + pointr * n] * d[j];
	}
	wv[i] = temp1;
	wv[col + i] = theta * temp2;
	pointr = pointr % m + 1;
    }

    m2 = m << 1;
    col2 = col << 1;

    dtrsl_(wn, m2, col2, wv, 11, info);
    if (*info != 0) {
	return 0;
    }
    for (i = 1; i <= col; ++i) {
	wv[i] = -wv[i];
    }
    dtrsl_(wn, m2, col2, wv, 1, info);
    if (*info != 0) {
	return 0;
    }

    pointr = head;
    for (jy = 1; jy <= col; ++jy) {
	js = col + jy;
	for (i = 1; i <= nsub; ++i) {
	    k = ind[i];
	    d[i] = d[i] + wy[k + pointr * n] * wv[jy] / theta 
		    + ws[k + pointr * n] * wv[js];
	}
	pointr = pointr % m + 1;
    }
    for (i = 1; i <= nsub; ++i) {
	d[i] /= theta;
    }

    alpha = 1.;
    temp1 = alpha;
    for (i = 1; i <= nsub; ++i) {
	k = ind[i];
	dk = d[i];
	if (nbd[k] != 0) {
	    if (dk < 0. && nbd[k] <= 2) {
		temp2 = l[k] - x[k];
		if (temp2 >= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha < temp2) {
		    temp1 = temp2 / dk;
		}
	    } else if (dk > 0. && nbd[k] >= 2) {
		temp2 = u[k] - x[k];
		if (temp2 <= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha > temp2) {
		    temp1 = temp2 / dk;
		}
	    }
	    if (temp1 < alpha) {
		alpha = temp1;
		ibd = i;
	    }
	}
    }

    if (alpha < 1.) {
	dk = d[ibd];
	k = ind[ibd];
	if (dk > 0.) {
	    x[k] = u[k];
	    d[ibd] = 0.;
	} else if (dk < 0.) {
	    x[k] = l[k];
	    d[ibd] = 0.;
	}
    }

    for (i = 1; i <= nsub; ++i) {
	k = ind[i];
	x[k] += alpha * d[i];
    }

    if (alpha < 1.) {
	*iword = 1;
    } else {
	*iword = 0;
    }
    return 0;
} 

static void timer_ (double *ttime)
{
    clock_t tim = clock();
    *ttime = (double) tim / CLOCKS_PER_SEC;
}

static double dpmeps_(void)
{
    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;

    double ret_val;
    double a, b;
    int i, it;
    double beta;
    int irnd;
    double temp, temp1, betah;
    int ibeta, negep;
    double tempa;
    int itemp;
    double betain;

    a = one;
    b = one;
L10:
    a += a;
    temp = a + one;
    temp1 = temp - a;
    if (temp1 - one == zero) {
	goto L10;
    }
L20:
    b += b;
    temp = a + b;
    itemp = (int) (temp - a);
    if (itemp == 0) {
	goto L20;
    }
    ibeta = itemp;
    beta = (double) ibeta;
    it = 0;
    b = one;
L30:
    ++it;
    b *= beta;
    temp = b + one;
    temp1 = temp - b;
    if (temp1 - one == zero) {
	goto L30;
    }
    irnd = 0;
    betah = beta / two;
    temp = a + betah;
    if (temp - a != zero) {
	irnd = 1;
    }
    tempa = a + beta;
    temp = tempa + betah;
    if (irnd == 0 && temp - tempa != zero) {
	irnd = 2;
    }
    negep = it + 3;
    betain = one / beta;
    a = one;
    for (i = 1; i <= negep; ++i) {
	a *= betain;
    }
L50:
    temp = one + a;
    if (temp - one != zero) {
	goto L60;
    }
    a *= beta;
    goto L50;
L60:
    ret_val = a;
    if (ibeta == 2 || irnd == 0) {
	goto L70;
    }
    a = a * (one + a) / two;
    temp = one + a;
    if (temp - one != zero) {
	ret_val = a;
    }
L70:
    return ret_val;
} 

static int mainlb_(int n, int m, double *x, 
		   double *l, double *u, int *nbd, double *f, double *g,
		   double reltol, double factr, double pgtol, 
		   double *ws, double *wy, double *sy, double *ss, 
		   double *yy, double *wt, 
		   double *wn, double *snd, double *z, double *r, 
		   double *d, double *t, double *wa, double *sg, 
		   double *sgo, double *yg, double *ygo, int *index, 
		   int *iwhere, int *indx2, char *task, char *csave, 
		   int *lsave, int *isave, double *dsave)
{
    double d1, d2;
    int i, k;
    double gd, dr, rr, dtd;
    double stp, cpu1, cpu2;
    double ddum, dnorm, gdold;
    double time, time1, time2;
    double xstep, stpmx;
    double epsmch, cachyt;
    double sbtime, sbgnrm, lnscht;
    double fold, theta;
    int col, wrk, head;
    int nact = 0;
    int info, iback, nfree;
    int nfgv, ifun, iter, nint;
    int itail = 0, iword = 0, ileave = 0;
    int boxed, nskip;
    int itfile = 0;
    int updatd, iupdat;
    int prjctd, cnstnd;
    int nintol, nenter = 0;

    /* fortran indexing fudges */

    --indx2;
    --iwhere;
    --index;
    --t;
    --d;
    --r;
    --z;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --ygo;
    --yg;
    --sgo;
    --sg;
    --wa;

    snd -= 1 + 2 * m;
    wn -= 1 + 2 * m;
    wt -= 1 + m;
    yy -= 1 + m;
    ss -= 1 + m;
    sy -= 1 + m;
    wy -= 1 + n;
    ws -= 1 + n;

    if (!strcmp(task, "START")) {
	timer_(&time1);

	epsmch = dpmeps_();

	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;

	iter = 0;
	nfgv = 0;
	nint = 0;
	nintol = 0;
	nskip = 0;
	nfree = n;

	if (reltol == 0 && na(reltol)) {
	    reltol = factr * epsmch;
	}

	cachyt = 0.;
	sbtime = 0.;
	lnscht = 0.;

	info = 0;

	errclb_(n, m, factr, l, u, nbd, task, &info, &k);
	if (!strncmp(task, "ERROR", 5)) {
	    return E_DATA;
	}

	active_(n, l, u, nbd, x, iwhere, &prjctd, &cnstnd, &boxed);
    } else {
	prjctd = lsave[0];
	cnstnd = lsave[1];
	boxed = lsave[2];
	updatd = lsave[3];
	nintol = isave[0];
	itfile = isave[2];
	iback = isave[3];
	nskip = isave[4];
	head = isave[5];
	col = isave[6];
	itail = isave[7];
	iter = isave[8];
	iupdat = isave[9];
	nint = isave[11];
	nfgv = isave[12];
	info = isave[13];
	ifun = isave[14];
	iword = isave[15];
	nfree = isave[16];
	nact = isave[17];
	ileave = isave[18];
	nenter = isave[19];
	theta = dsave[0];
	fold = dsave[1];
	reltol = dsave[2];
	dnorm = dsave[3];
	epsmch = dsave[4];
	cpu1 = dsave[5];
	cachyt = dsave[6];
	sbtime = dsave[7];
	lnscht = dsave[8];
	time1 = dsave[9];
	gd = dsave[10];
	stpmx = dsave[11];
	sbgnrm = dsave[12];
	stp = dsave[13];
	gdold = dsave[14];
	dtd = dsave[15];

	if (!strncmp(task, "FG_LN", 5)) {
	    goto L666;
	}
	if (!strncmp(task, "NEW_X", 5)) {
	    goto L777;
	}
	if (!strncmp(task, "FG_ST", 5)) {
	    goto L111;
	}
	if (!strncmp(task, "STOP", 4)) {
	    goto L999;
	}
    }

    strcpy(task, "FG_START");

    goto L1000;

L111:
    nfgv = 1;
    projgr_(n, l, u, nbd, x, g, &sbgnrm);

    if (sbgnrm <= pgtol) {
	strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }

    /* start of main loop  */

L222:
    iword = -1;

    if (!cnstnd && col > 0) {
	dcopy_(n, x, z);
	wrk = updatd;
	nint = 0;
	goto L333;
    }

    timer_(&cpu1);

    cauchy_(n, x, l, u, nbd, g, indx2, iwhere, 
	    t, d, z, m, wy, ws, sy, wt,
	    theta, col, head, wa, 
	    &wa[(m << 1)], &wa[(m << 2)], &wa[m * 6], &nint, 
	    sg, yg, sbgnrm, &info, epsmch);

    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	timer_(&cpu2);
	cachyt = cachyt + cpu2 - cpu1;
	goto L222;
    }

    timer_(&cpu2);
    cachyt = cachyt + cpu2 - cpu1;
    nintol += nint;

    freev_(n, &nfree, index, &nenter, &ileave, indx2, iwhere, 
	   &wrk, updatd, cnstnd, iter);
    nact = n - nfree;

L333:
    if (nfree == 0 || col == 0) {
	goto L555;
    }

    timer_(&cpu1);

    if (wrk) {
	formk_(n, nfree, index, nenter, ileave, indx2, iupdat, 
	       updatd, wn, snd, m, ws, wy, sy, theta, col, 
	       head, &info);
    }

    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }

    cmprlb_(n, m, x, g, ws, wy, sy, wt, z, r, wa, 
	    index, theta, col, head, nfree, cnstnd, &info);
    if (info != 0) {
	goto L444;
    }

    subsm_(n, m, nfree, index, l, u, nbd, z, r, 
	   ws, wy, theta, col, head, &iword, 
	   wa, wn, &info);

 L444:
    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    sbtime = sbtime + cpu2 - cpu1;
L555:

    for (i = 1; i <= n; ++i) {
	d[i] = z[i] - x[i];
    }

    timer_(&cpu1);

L666:
    lnsrlb_(n, l, u, nbd, x, *f, &fold, &gd, &gdold, g,
	    d, r, t, z, &stp, &dnorm, &dtd, &xstep, 
	    &stpmx, iter, &ifun, &iback, &nfgv, &info, 
	    task, boxed, cnstnd, csave, &isave[20], 
	    &dsave[15]);

    if (info != 0 || iback >= 20) {
	dcopy_(n, t, x);
	dcopy_(n, r, g);
	*f = fold;
	if (col == 0) {
	    if (info == 0) {
		info = -9;
		--nfgv;
		--ifun;
		--iback;
	    }
	    strcpy(task, "ABNORMAL_TERMINATION_IN_LNSRCH");
	    ++iter;
	    goto L999;
	} else {
	    if (info == 0) {
		--nfgv;
	    }
	    info = 0;
	    col = 0;
	    head = 1;
	    theta = 1.;
	    iupdat = 0;
	    updatd = 0;
	    strcpy(task, "RESTART_FROM_LNSRCH");
	    timer_(&cpu2);
	    lnscht = lnscht + cpu2 - cpu1;
	    goto L222;
	}
    } else if (!strncmp(task, "FG_LN", 5)) {
	goto L1000;
    } else {
	timer_(&cpu2);
	lnscht = lnscht + cpu2 - cpu1;
	++iter;
	projgr_(n, l, u, nbd, x, g, &sbgnrm);
	goto L1000;
    }

L777:

    if (sbgnrm <= pgtol) {
	strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }

    d1 = abs(fold), d2 = abs(*f), d1 = max(d1,d2);
    ddum = max(d1, 1.);
    if (fold - *f <= reltol * ddum) {
	sprintf(task, "CONVERGENCE: REL_REDUCTION_OF_F <= %g", reltol);
	if (iback >= 10) {
	    info = -5;
	}
	goto L999;
    }

    for (i = 1; i <= n; ++i) {
	r[i] = g[i] - r[i];
    }
    rr = ddot_(n, r, r);

    if (stp == 1.) {
	dr = gd - gdold;
	ddum = -gdold;
    } else {
	dr = (gd - gdold) * stp;
	dscal_(n, stp, d);
	ddum = -gdold * stp;
    }

    if (dr <= epsmch * ddum) {
	++nskip;
	updatd = 0;
	goto L888;
    }

    updatd = 1;
    ++iupdat;

    matupd_(n, m, ws, wy, sy, ss, d, r, &itail, iupdat, &col, 
	    &head, &theta, rr, dr, &stp, &dtd);

    formt_(m, wt, sy, ss, col, theta, &info);
    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	goto L222;
    }

L888:

    /* end of main loop  */
    goto L222;

L999:
    timer_(&time2);
    time = time2 - time1;

L1000:
    lsave[0] = prjctd;
    lsave[1] = cnstnd;
    lsave[2] = boxed;
    lsave[3] = updatd;

    isave[0] = nintol;
    isave[2] = itfile;
    isave[3] = iback;
    isave[4] = nskip;
    isave[5] = head;
    isave[6] = col;
    isave[7] = itail;
    isave[8] = iter;
    isave[9] = iupdat;
    isave[11] = nint;
    isave[12] = nfgv;
    isave[13] = info;
    isave[14] = ifun;
    isave[15] = iword;
    isave[16] = nfree;
    isave[17] = nact;
    isave[18] = ileave;
    isave[19] = nenter;

    dsave[0] = theta;
    dsave[1] = fold;
    dsave[2] = reltol;
    dsave[3] = dnorm;
    dsave[4] = epsmch;
    dsave[5] = cpu1;
    dsave[6] = cachyt;
    dsave[7] = sbtime;
    dsave[8] = lnscht;
    dsave[9] = time1;
    dsave[10] = gd;
    dsave[11] = stpmx;
    dsave[12] = sbgnrm;
    dsave[13] = stp;
    dsave[14] = gdold;
    dsave[15] = dtd;

    return 0;
} 

static int setulb_(int n, int m, double *x, 
		   double *l, double *u, int *nbd, double *f, double *g,
		   double reltol, double factr, double pgtol, 
		   double *wa, int *iwa,
		   char *task, char *csave, int *lsave, 
		   int *isave, double *dsave)
{
    int ld, lr, lt, lz, lwa, lsg, lyg, lwn, lss;
    int lws, lwt, lsy, lwy, lyy, lsnd, lsgo, lygo;

    if (!strcmp(task, "START")) {
	isave[0] = m * n;
	isave[1] = m * m;
	isave[2] = m * m << 2;
	isave[3] = 0;
	isave[4] = isave[3] + isave[0];
	isave[5] = isave[4] + isave[0];
	isave[6] = isave[5] + isave[1];
	isave[7] = isave[6] + isave[1];
	isave[8] = isave[7] + isave[1];
	isave[9] = isave[8] + isave[1];
	isave[10] = isave[9] + isave[2];
	isave[11] = isave[10] + isave[2];
	isave[12] = isave[11] + n;
	isave[13] = isave[12] + n;
	isave[14] = isave[13] + n;
	isave[15] = isave[14] + n;
	isave[16] = isave[15] + (m << 3);
	isave[17] = isave[16] + m;
	isave[18] = isave[17] + m;
	isave[19] = isave[18] + m;
    }

    lws = isave[3];
    lwy = isave[4];
    lsy = isave[5];
    lss = isave[6];
    lyy = isave[7];
    lwt = isave[8];
    lwn = isave[9];
    lsnd = isave[10];
    lz = isave[11];
    lr = isave[12];
    ld = isave[13];
    lt = isave[14];
    lwa = isave[15];
    lsg = isave[16];
    lsgo = isave[17];
    lyg = isave[18];
    lygo = isave[19];

    int err;

    err = mainlb_(n, m, x, l, u, nbd, f, g, reltol, factr, pgtol, 
		  &wa[lws], &wa[lwy], &wa[lsy], &wa[lss], &wa[lyy], &wa[lwt], 
		  &wa[lwn], &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], 
		  &wa[lwa], &wa[lsg], &wa[lsgo], &wa[lyg], &wa[lygo], 
		  iwa, &iwa[n], &iwa[(n << 1)], task, csave, 
		  lsave, &isave[21], dsave);

    return err;
}

/* description of same of the key variables below, adapted
   from the FORTRAN file routines.f:

   factr (double):

   factr >= 0 is specified by the user.  The iteration
   will stop when

   (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch

   where epsmch is the machine precision, which is automatically
   generated by the code. Typical values for factr: 1.d+12 for
   low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
   high accuracy.

   pgtol (double):

   pgtol >= 0 is specified by the user.  The iteration
   will stop when

     max{|proj g_i | i = 1, ..., n} <= pgtol

   where pg_i is the ith component of the projected gradient.   

   dsave (double array):
   On exit with task = NEW_X, the following information is available:

   dsave[0] = current 'theta' in the BFGS matrix
   dsave[1] = f(x) in the previous iteration
   dsave[2] = factr*epsmch
   dsave[3] = 2-norm of the line search direction vector
   dsave[4] = the machine precision epsmch generated by the code
   dsave[6] = the accumulated time spent on searching for Cauchy points
   dsave[7] = the accumulated time spent on subspace minimization
   dsave[8] = the accumulated time spent on line search
   dsave[10] = the slope of the line search function at the current 
               point of line search
   dsave[11] = the maximum relative step length imposed in line search
   dsave[12] = the infinity norm of the projected gradient
   dsave[13] = the relative step length in the line search
   dsave[14] = the slope of the line search function at the starting 
               point of the line search
   dsave[15] = the square of the 2-norm of the line search direction 
               vector
*/

int BFGS_test (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn)
{
    double *g = NULL;
    double *l = NULL;
    double *u = NULL;
    double *wa = NULL;
    int *nbd = NULL;
    int *iwa = NULL;

    int i, m, wadim;
    char task[60];
    char csave[60];
    double f, factr, pgtol;
    double dsave[29];
    int isave[44];
    int lsave[4];
    int err = 0;

    BFGS_get_user_values(b, n, &maxit, &reltol, opt, prn);

    /*
      m: the number of corrections used in the limited memory matrix.
      It is not altered by the routine.  Values of m < 3 are not
      recommended, and large values of m can result in excessive
      computing time. The range 3 <= m <= 20 is recommended.
    */
    m = 10; /* was initially set to 5 */

    wadim = (2*m+4)*n + 12*m*m + 12*m;

    g = malloc(n * sizeof *g);
    l = malloc(n * sizeof *l);
    u = malloc(n * sizeof *u);
    wa = malloc(wadim * sizeof *wa);
    nbd = malloc(n * sizeof *nbd);
    iwa = malloc(3*n * sizeof *iwa);

    if (g == NULL || l == NULL || u == NULL ||
	wa == NULL || nbd == NULL || iwa == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (gradfunc == NULL) {
	gradfunc = BFGS_numeric_gradient;
    }

    /* Convergence criteria: we don't use "reltol" here because I
       haven't yet figured out exactly how it relates to the stopping
       criteria used in this implementation. See the long comment
       above this function.
    */
    factr = 1e5;  /* ?? */
    pgtol = 0;

    /* Bounds on the parameters: for now we just set them all to be
       less than some ridiculously large number */
    for (i=0; i<n; i++) {
	nbd[i] = 3; /* case 3: upper bound only */
	u[i] = NADBL / 100;
    }	

    /* Start the iteration by initializing 'task' */
    strcpy(task, "START");

    while (1) {
	/* Call the L-BFGS-B code */
	setulb_(n, m, b, l, u, nbd, &f, g, reltol, factr, pgtol, wa, iwa, 
		task, csave, lsave, isave, dsave);

	if (!strncmp(task, "FG", 2)) {
	    double minusf;

	    /* Compute function value, f */
	    minusf = cfunc(b, data);
	    if (!na(minusf)) {
		f = -minusf;
	    }
	    *fncount += 1;

	    /* Compute gradient, g */
	    gradfunc(b, g, n, cfunc, data);
	    reverse_gradient(g, n);
	    *grcount += 1;

	} else if (!strncmp(task, "NEW_X", 5)) {
	    /* The optimizer has produced a new set of parameter values */
	    if (isave[33] >= maxit) {
		strcpy(task, "STOP: TOTAL NO. of f AND g "
		       "EVALUATIONS EXCEEDS LIMIT");
		err = E_NOCONV;
	    } else if (opt & OPT_V) {
		reverse_gradient(g, n);
		print_iter_info(isave[29], -f, crittype, n, b, g, dsave[13], prn);
		reverse_gradient(g, n);
	    }
	} else {
	    /* study up on other possible "task" contents! */
	    fprintf(stderr, "got task = '%s'\n", task);
	    break;
	}
    }

    if (opt & OPT_V) {
	pputs(prn, _("\n--- FINAL VALUES: \n"));
	reverse_gradient(g, n); /* ?? */
	print_iter_info(isave[29], -f, crittype, n, b, g, dsave[13], prn);
	pputc(prn, '\n');
    }

 bailout:

    free(g);
    free(l);
    free(u);
    free(wa);
    free(nbd);
    free(iwa);

    return err;
}
