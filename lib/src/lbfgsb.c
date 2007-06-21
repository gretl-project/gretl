/* This code is ased on L-BFGS-B (version 2.1); the core
   is an f2c translation of routines.f

   See http://www.ece.northwestern.edu/~nocedal/lbfgsb.html

*/

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

/* dot product of two vectors -- uses unrolled loops for increments
   equal to one.  jack dongarra, linpack, 3/11/78. */

static double ddot_(int n, double *dx, double *dy)
{
    int i, m, mp1;
    double dtemp = 0;

    /* Parameter adjustments */
    --dy;
    --dx;

    if (n <= 0) {
	return dtemp;
    }

    m = n % 5;
    if (m == 0) {
	goto L40;
    }
    for (i = 1; i <= m; ++i) {
	dtemp += dx[i] * dy[i];
    }
    if (n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 5) {
	dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * 
	    dy[i + 1] + dx[i + 2] * dy[i + 2] + dx[i + 3] * 
	    dy[i + 3] + dx[i + 4] * dy[i + 4];
    }
L60:

    return dtemp;
} /* ddot_ */

/* scales a vector by a constant -- uses unrolled loops for increment
   equal to one.  jack dongarra, linpack, 3/11/78. 
*/

static void dscal_(int n, double da, double *dx)
{
    int i, m, mp1;

    /* Parameter adjustments */
    --dx;

    if (n <= 0) {
	return;
    }

    m = n % 5;
    if (m == 0) {
	goto L40;
    }
    for (i = 1; i <= m; ++i) {
	dx[i] = da * dx[i];
    }
    if (n < 5) {
	return;
    }
 L40:
    mp1 = m + 1;
    for (i = mp1; i <= n; i += 5) {
	dx[i] = da * dx[i];
	dx[i + 1] = da * dx[i + 1];
	dx[i + 2] = da * dx[i + 2];
	dx[i + 3] = da * dx[i + 3];
	dx[i + 4] = da * dx[i + 4];
    }
} /* dscal_ */

static void dpofa_(double *a, int lda, int n, int *info)
{
    int j, k;
    double s, t;
    int jm1;

    /* Parameter adjustments */
    a -= 1 + lda;

    for (j = 1; j <= n; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	for (k = 1; k <= jm1; ++k) {
	    t = a[k + j * lda] - ddot_(k - 1, &a[k * lda + 1], &a[j * lda + 1]);
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
} /* dpofa_ */

static void dcopy_(int n, double *dx, double *dy)
{
    int i, m, mp1;

    /* Parameter adjustments */
    --dy;
    --dx;

    if (n <= 0) {
	return;
    }

    m = n % 7;

    if (m != 0) {
	for (i = 1; i <= m; ++i) {
	    dy[i] = dx[i];
	}
	if (n < 7) {
	    return;
	}
    }

    mp1 = m + 1;
    for (i = mp1; i <= n; i += 7) {
	dy[i] = dx[i];
	dy[i + 1] = dx[i + 1];
	dy[i + 2] = dx[i + 2];
	dy[i + 3] = dx[i + 3];
	dy[i + 4] = dx[i + 4];
	dy[i + 5] = dx[i + 5];
	dy[i + 6] = dx[i + 6];
    }
} /* dcopy_ */

/* constant times a vector plus a vector -- uses unrolled loops for
   increments equal to one.  jack dongarra, linpack, 3/11/78.
*/

static void daxpy_(int n, double da, double *dx, double *dy)
{
    int i, m, mp1;

    /* Parameter adjustments */
    --dy;
    --dx;

    if (n <= 0 || da == 0) {
	return;
    }

    m = n % 4;

    if (m != 0) {
	for (i = 1; i <= m; ++i) {
	    dy[i] += da * dx[i];
	}
	if (n < 4) {
	    return;
	}
    }

    mp1 = m + 1;
    for (i = mp1; i <= n; i += 4) {
	dy[i] += da * dx[i];
	dy[i + 1] += da * dx[i + 1];
	dy[i + 2] += da * dx[i + 2];
	dy[i + 3] += da * dx[i + 3];
    }
} /* daxpy_ */

/* dtrsl solves systems of the form t * x = b or trans(t) * x = b,
   where t is a triangular matrix of order n
*/

static void dtrsl_(double *t, int ldt, int n, 
		   double *b, int job, int *info)
{
    int j, jj, case__ = 1;
    double temp;

    /* Parameter adjustments */
    t -= 1 + ldt;
    --b;

    for (*info = 1; *info <= n; ++(*info)) {
	if (t[*info + *info * ldt] == 0.) {
	    return;
	}
    }

    *info = 0;

    /* determine the task */
    if (job % 10 != 0) {
	case__ = 2;
    }
    if (job % 100 / 10 != 0) {
	case__ += 2;
    }

    if (case__ == 1) {
	/* solve t*x=b for t lower triangular */
	b[1] /= t[ldt + 1];
	for (j = 2; j <= n; ++j) {
	    temp = -b[j - 1];
	    daxpy_(n - j + 1, temp, &t[j + (j - 1) * ldt], &b[j]);
	    b[j] /= t[j + j * ldt];
	}
    } else if (case__ == 2) {
	/* solve t*x=b for t upper triangular. */
	b[n] /= t[n + n * ldt];
	for (jj = 2; jj <= n; ++jj) {
	    j = n - jj + 1;
	    temp = -b[j + 1];
	    daxpy_(j, temp, &t[(j + 1) * ldt + 1], &b[1]);
	    b[j] /= t[j + j * ldt];
	}
    } else if (case__ == 3) {
	/* solve trans(t)*x=b for t lower triangular. */
	b[n] /= t[n + n * ldt];
	for (jj = 2; jj <= n; ++jj) {
	    j = n - jj + 1;
	    b[j] -= ddot_(jj - 1, &t[j + 1 + j * ldt], &b[j + 1]);
	    b[j] /= t[j + j * ldt];
	}
    } else if (case__ == 4) {
	/* solve trans(t)*x=b for t upper triangular. */
	b[1] /= t[ldt + 1];
	for (j = 2; j <= n; ++j) {
	    b[j] -= ddot_(j - 1, &t[j * ldt + 1], &b[1]);
	    b[j] /= t[j + j * ldt];
	}
    }
} /* dtrsl_ */

static int active_(int n, double *l, double *u, 
		   int *nbd, double *x, int *iwhere,  
		   int *prjctd, int *cnstnd, int *boxed)
{
    int i, nbdd;

    /* Parameter adjustments */
    --iwhere;
    --x;
    --nbd;
    --u;
    --l;

    nbdd = 0;
    *prjctd = 0;
    *cnstnd = 0;
    *boxed = 1;

    /* Project the initial x to the easible set if necessary. */
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

    /* Initialize iwhere and assign values to cnstnd and boxed. */
    for (i = 1; i <= n; ++i) {
	if (nbd[i] != 2) {
	    *boxed = 0;
	}
	if (nbd[i] == 0) {
	    /* this variable is always free */
	    iwhere[i] = -1;
	    /* otherwise set x(i)=mid(x(i), u(i), l(i)). */
	} else {
	    *cnstnd = 1;
	    if (nbd[i] == 2 && u[i] - l[i] <= 0.) {
		/* this variable is always fixed */
		iwhere[i] = 3;
	    } else {
		iwhere[i] = 0;
	    }
	}
    }
    return 0;
} /* active_ */

static void bmv_(int m, double *sy, double *wt, int col, 
		 double *v, double *p, int *info)
{
    /* System generated locals */
    int wt_offset;

    int i, k, icol;
    double sum;

    /* Parameter adjustments */
    wt_offset = 1 + m;
    wt -= wt_offset;
    sy -= 1 + m;
    --p;
    --v;

    if (col == 0) {
	return;
    }

    /*     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ] */
    /*                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ]. */
    /* 	solve Jp2=v2+LD^(-1)v1. */
    p[col + 1] = v[col + 1];
    for (i = 2; i <= col; ++i) {
	icol = col + i;
	sum = 0.;
	for (k = 1; k <= i - 1; ++k) {
	    sum += sy[i + k * m] * v[k] / sy[k + k * m];
	}
	p[icol] = v[icol] + sum;
    }

    /* Solve the triangular system */
    dtrsl_(&wt[wt_offset], m, col, &p[col + 1], 11, info);
    if (*info != 0) {
	return;
    }
    /* solve D^(1/2)p1=v1. */
    for (i = 1; i <= col; ++i) {
	p[i] = v[i] / sqrt(sy[i + i * m]);
    }

    /*     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ] */
    /*                    [  0         J'           ] [ p2 ]   [ p2 ]. */
    /*       solve J^Tp2=p2. */
    dtrsl_(&wt[wt_offset], m, col, &p[col + 1], 1, info);
    if (*info != 0) {
	return;
    }
    
    /*       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2) */
    /*                 =-D^(-1/2)p1+D^(-1)L'p2. */
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
} /* bmv_ */

static int hpsolb_(int n, double *t, int *iorder, int iheap)
{
    int i, j, k;
    double out, ddum;
    int indxin, indxou;

    /* Parameter adjustments */
    --iorder;
    --t;

    if (iheap == 0) {
	/* Rearrange the elements t(1) to t(n) to form a heap. */
	for (k = 2; k <= n; ++k) {
	    ddum = t[k];
	    indxin = iorder[k];
	    /* Add ddum to the heap. */
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

    /* Assign to 'out' the value of t(1), the least member of the heap, 
       and rearrange the remaining members to form a heap as 
       elements 1 to n-1 of t. */
    if (n > 1) {
	i = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[n];
	indxin = iorder[n];
	/* Restore the heap */
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
	/* Put the least member in t(n). */
	t[n] = out;
	iorder[n] = indxou;
    }
    return 0;
} /* hpsolb_ */

static int 
cauchy_(int n, double *x, double *l, double *u, int *nbd, 
	double *g, int *iorder, int *iwhere, double *t, 
	double *d__, double *xcp, int m, 
	double *wy, double *ws, double *sy, double *wt, 
	double theta, int col, int head, double *p, 
	double *c__, double *wbp, double *v, int *nint, 
	double *sg, double *yg, double sbgnrm, int *info, 
	double epsmch)
{
    /* System generated locals */
    int wy_offset, ws_offset, sy_offset, wt_offset;
    double d__1;

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

    /* Parameter adjustments */
    --xcp;
    --d__;
    --t;
    --iwhere;
    --iorder;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --yg;
    --sg;
    --v;
    --wbp;
    --c__;
    --p;

    wt_offset = 1 + m;
    wt -= wt_offset;
    sy_offset = 1 + m;
    sy -= sy_offset;
    ws_offset = 1 + n;
    ws -= ws_offset;
    wy_offset = 1 + n;
    wy -= wy_offset;

    if (sbgnrm <= 0.) {
	dcopy_(n, &x[1], &xcp[1]);
	return 0;
    }

    bnded = 1;
    nfree = n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = col << 1;
    f1 = 0.;

    /* We set p to zero and build it up as we determine d. */
    for (i = 1; i <= col2; ++i) {
	p[i] = 0.;
    }

    /* In the following loop we determine for each variable its bound */
    /* status and its breakpoint, and update p accordingly. */
    /* Smallest breakpoint is identified. */

    for (i = 1; i <= n; ++i) {
	neggi = -g[i];
	if (iwhere[i] != 3 && iwhere[i] != -1) {
	    /* if x(i) is not a constant and has bounds, */
	    /* compute the difference between x(i) and its bounds. */
	    if (nbd[i] <= 2) {
		tl = x[i] - l[i];
	    }
	    if (nbd[i] >= 2) {
		tu = u[i] - x[i];
	    }
	    /* If a variable is close enough to a bound */
	    /* we treat it as at bound. */
	    xlower = nbd[i] <= 2 && tl <= 0.;
	    xupper = nbd[i] >= 2 && tu <= 0.;
	    /* reset iwhere(i). */
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
	    d__[i] = 0.;
	} else {
	    d__[i] = neggi;
	    f1 -= neggi * neggi;
	    /* calculate p := p - W'e_i* (g_i). */
	    for (j = 1; j <= col; ++j) {
		p[j] += wy[i + pointr * n] * neggi;
		p[col + j] += ws[i + pointr * n] * neggi;
		pointr = pointr % m + 1;
	    }
	    if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.) {
		/* x(i) + d(i) is bounded; compute t(i). */
		++nbreak;
		iorder[nbreak] = i;
		t[nbreak] = tl / (-neggi);
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else if (nbd[i] >= 2 && neggi > 0.) {
		/* x(i) + d(i) is bounded; compute t(i). */
		++nbreak;
		iorder[nbreak] = i;
		t[nbreak] = tu / neggi;
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else {
		/* x(i) + d(i) is not bounded. */
		--nfree;
		iorder[nfree] = i;
		if (abs(neggi) > 0.) {
		    bnded = 0;
		}
	    }
	}
    }

    /* The indices of the nonzero components of d are now stored */
    /* in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
    /* The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
    if (theta != 1.) {
	/* complete the initialization of p for theta not= one. */
	dscal_(col, theta, &p[col + 1]);
    }

    /* Initialize GCP xcp = x. */
    dcopy_(n, &x[1], &xcp[1]);
    if (nbreak == 0 && nfree == n + 1) {
	/* is a zero vector, return with the initial xcp as GCP. */
	return 0;
    }

    /* Initialize c = W'(xcp - x) = 0. */
    for (j = 1; j <= col2; ++j) {
	c__[j] = 0.;
    }

    /* Initialize derivative f2. */
    f2 = -theta * f1;
    f2_org__ = f2;
    if (col > 0) {
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	f2 -= ddot_(col2, &v[1], &p[1]);
    }
    dtm = -f1 / f2;
    tsum = 0.;
    *nint = 1;

    /* If there are no breakpoints, locate the GCP and return. */
    if (nbreak == 0) {
	goto L888;
    }
    nleft = nbreak;
    iter = 1;
    tj = 0.;

/* ------------------- the beginning of the loop ------------------------- */
L777:
    /* Find the next smallest breakpoint; */

    /* compute dt = t(nleft) - t(nleft + 1). */
    tj0 = tj;

    if (iter == 1) {
	/* Since we already have the smallest breakpoint we need not do */
	/* heapsort yet. Often only one breakpoint is used and the */
	/* cost of heapsort is avoided. */
	tj = bkmin;
	ibp = iorder[ibkmin];
    } else {
	if (iter == 2) {
	    /* Replace the already used smallest breakpoint with the */
	    /* breakpoint numbered nbreak > nlast, before heapsort call. */
	    if (ibkmin != nbreak) {
		t[ibkmin] = t[nbreak];
		iorder[ibkmin] = iorder[nbreak];
	    }
	    /* Update heap structure of breakpoints */
	    /* (if iter=2, initialize heap). */
	}
	hpsolb_(nleft, &t[1], &iorder[1], iter - 2);
	tj = t[nleft];
	ibp = iorder[nleft];
    }

    dt = tj - tj0;

    /* If a minimizer is within this interval, locate the GCP and return. */
    if (dtm < dt) {
	goto L888;
    }

    /* Otherwise fix one variable and reset the corresponding
       component of d to zero. */
    tsum += dt;
    --nleft;
    ++iter;
    dibp = d__[ibp];
    d__[ibp] = 0.;
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
	/* all n variables are fixed, return with xcp as GCP. */
	dtm = dt;
	goto L999;
    }

    /* Update the derivative information. */
    ++(*nint);
    /* Computing 2nd power */
    d__1 = dibp;
    dibp2 = d__1 * d__1;
    /* Update f1 and f2. */
    /* temporarily set f1 and f2 for col=0. */
    f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
    f2 -= theta * dibp2;
    if (col > 0) {
	/* update c = c + dt*p. */
	daxpy_(col2, dt, &p[1], &c__[1]);
	/* choose wbp, */
	/* the row of W corresponding to the breakpoint encountered. */
	pointr = head;
	for (j = 1; j <= col; ++j) {
	    wbp[j] = wy[ibp + pointr * n];
	    wbp[col + j] = theta * ws[ibp + pointr * n];
	    pointr = pointr % m + 1;
	}
	/* compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	wmc = ddot_(col2, &c__[1], &v[1]);
	wmp = ddot_(col2, &p[1], &v[1]);
	wmw = ddot_(col2, &wbp[1], &v[1]);
	/* update p = p - dibp*wbp. */
	d__1 = -dibp;
	daxpy_(col2, d__1, &wbp[1], &p[1]);
	/* complete updating f1 and f2 while col > 0. */
	f1 += dibp * wmc;
	f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }

    /* Computing MAX */
    d__1 = epsmch * f2_org__;
    f2 = max(d__1,f2);
    if (nleft > 0) {
	dtm = -f1 / f2;
	goto L777;
	/* to repeat the loop for unsearched intervals. */
    } else if (bnded) {
	f1 = 0.;
	f2 = 0.;
	dtm = 0.;
    } else {
	dtm = -f1 / f2;
    }

/* ------------------- the end of the loop ------------------------------- */

L888:
    if (dtm <= 0.) {
	dtm = 0.;
    }
    tsum += dtm;

    /* Move free variables (i.e., the ones w/o breakpoints) and */
    /* the variables whose breakpoints haven't been reached. */
    daxpy_(n, tsum, &d__[1], &xcp[1]);

L999:
    /* Update c = c + dtm*p = W'(x^c - x) */
    /* which will be used in computing r = Z'(B(x^c - x) + g). */
    if (col > 0) {
	daxpy_(col2, dtm, &p[1], &c__[1]);
    }
    return 0;
} /* cauchy_ */

static int 
cmprlb_(int n, int m, double *x, 
	double *g, double *ws, double *wy, double *sy, 
	double *wt, double *z, double *r, double *wa, 
	int *index, double theta, int col, int head, 
	int nfree, int cnstnd, int *info)
{
    /* System generated locals */
    int ws_offset, wy_offset, sy_offset, wt_offset;

    int i, j, k;
    double a1, a2;
    int pointr;

    /* Parameter adjustments */
    --index;
    --r;
    --z;
    --g;
    --x;
    --wa;
    wt_offset = 1 + m;
    wt -= wt_offset;
    sy_offset = 1 + m;
    sy -= sy_offset;
    wy_offset = 1 + n;
    wy -= wy_offset;
    ws_offset = 1 + n;
    ws -= ws_offset;

    if (!cnstnd && col > 0) {
	for (i = 1; i <= n; ++i) {
	    r[i] = -g[i];
	}
    } else {
	for (i = 1; i <= nfree; ++i) {
	    k = index[i];
	    r[i] = -(theta) * (z[k] - x[k]) - g[k];
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wa[(m << 1) + 1], 
	     &wa[1], info);
	if (*info != 0) {
	    *info = -8;
	    return 0;
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
    return 0;
} /* cmprlb_ */

static int errclb_(int n, int m, double factr, 
		   double *l, double *u, int *nbd, 
		   char *task, int *info, int *k)
{
    int i;

    /* Parameter adjustments */
    --nbd;
    --u;
    --l;

    if (n <= 0) {
	strcpy(task, "ERROR: N .LE. 0");
    }
    if (m <= 0) {
	strcpy(task, "ERROR: M .LE. 0");
    }
    if (factr < 0.) {
	strcpy(task, "ERROR: FACTR .LT. 0");
    }

    /* Check the validity of the arrays nbd(i), u(i), and l(i). */
    for (i = 1; i <= n; ++i) {
	if (nbd[i] < 0 || nbd[i] > 3) {
	    /* return */
	    strcpy(task, "ERROR: INVALID NBD");
	    *info = -6;
	    *k = i;
	}
	if (nbd[i] == 2) {
	    if (l[i] > u[i]) {
		/* return */
		strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
		*info = -7;
		*k = i;
	    }
	}
    }

    return 0;
} /* errclb_ */

static int 
formk_(int n, int nsub, int *ind, int nenter, int ileave, 
       int *indx2, int iupdat, int updatd, double *wn, double *wn1, 
       int m, double *ws, double *wy, double *sy, double theta, 
       int col, int head, int *info)
{
    /* System generated locals */
    int wn_offset, wn1_offset, ws_offset, wy_offset, sy_offset;

    int i, k, k1, m2, is, js, iy, jy, is1, js1, col2, dend, pend;
    int upcl, two_m;
    double temp1, temp2, temp3, temp4;
    int ipntr, jpntr, dbegin, pbegin;

    /* Parameter adjustments */
    --indx2;
    --ind;
    sy_offset = 1 + m;
    sy -= sy_offset;
    wy_offset = 1 + n;
    wy -= wy_offset;
    ws_offset = 1 + n;
    ws -= ws_offset;
    two_m = 2 * m;
    wn1_offset = 1 + two_m;
    wn1 -= wn1_offset;
    wn_offset = 1 + two_m;
    wn -= wn_offset;

    if (updatd) {
	if (iupdat > m) {
	    /* shift old part of WN1. */
	    for (jy = 1; jy <= m - 1; ++jy) {
		js = m + jy;
		dcopy_(m - jy, &wn1[jy + 1 + (jy + 1) * two_m], 
		       &wn1[jy + jy * two_m]);
		dcopy_(m - jy, &wn1[js + 1 + (js + 1) * two_m], 
		       &wn1[js + js * two_m]);
		dcopy_(m - 1, &wn1[m + 2 + (jy + 1) * two_m], 
		       &wn1[m + 1 + jy * two_m]);
	    }
	}
	/* put new rows in blocks (1,1), (2,1) and (2,2). */
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
	    /* compute element jy of row 'col' of Y'ZZ'Y */
	    for (k = pbegin; k <= pend; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    /* compute elements jy of row 'col' of L_a and S'AA'S */
	    for (k = dbegin; k <= dend; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * n] * ws[k1 + jpntr * n];
		temp3 += ws[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    wn1[iy + jy * two_m] = temp1;
	    wn1[is + js * two_m] = temp2;
	    wn1[is + jy * two_m] = temp3;
	    jpntr = jpntr % m + 1;
	}

	/* put new column in block (2,1). */
	jy = col;
	jpntr = head + col - 1;
	if (jpntr > m) {
	    jpntr -= m;
	}
	ipntr = head;
	for (i = 1; i <= col; ++i) {
	    is = m + i;
	    temp3 = 0.;
	    /* compute element i of column 'col' of R_z */
	    for (k = pbegin; k <= pend; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * n] * wy[k1 + jpntr * n];
	    }
	    ipntr = ipntr % m + 1;
	    wn1[is + jy * two_m] = temp3;
	}
	upcl = col - 1;
    } else {
	upcl = col;
    }

    /* modify the old parts in blocks (1,1) and (2,2) due to changes */
    /* in the set of free variables. */
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
	    wn1[iy + jy * two_m] = wn1[iy + jy * two_m] + temp1 - temp3;
	    wn1[is + js * two_m] = wn1[is + js * two_m] - temp2 + temp4;
	    jpntr = jpntr % m + 1;
	}
	ipntr = ipntr % m + 1;
    }

    /* modify the old parts in block (2,1). */
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
		wn1[is + jy * two_m] = wn1[is + jy * two_m] + temp1 - temp3;
	    } else {
		wn1[is + jy * two_m] = wn1[is + jy * two_m] - temp1 + temp3;
	    }
	    jpntr = jpntr % m + 1;
	}
	ipntr = ipntr % m + 1;
    }

    /* Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] */
    /*                                 [-L_a +R_z        S'AA'S*theta] */
    m2 = m << 1;
    for (iy = 1; iy <= col; ++iy) {
	is = col + iy;
	is1 = m + iy;
	for (jy = 1; jy <= iy; ++jy) {
	    js = col + jy;
	    js1 = m + jy;
	    wn[jy + iy * two_m] = wn1[iy + jy * two_m] / theta;
	    wn[js + is * two_m] = wn1[is1 + js1 * two_m] * theta;
	}
	for (jy = 1; jy <= iy - 1; ++jy) {
	    wn[jy + is * two_m] = -wn1[is1 + jy * two_m];
	}
	for (jy = iy; jy <= col; ++jy) {
	    wn[jy + is * two_m] = wn1[is1 + jy * two_m];
	}
	wn[iy + iy * two_m] += sy[iy + iy * m];
    }

    /* Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] */
    /*                                [(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
    /* first Cholesky factor (1,1) block of wn to get LL' */
    /* with L' stored in the upper triangle of wn. */
    dpofa_(&wn[wn_offset], m2, col, info);
    if (*info != 0) {
	*info = -1;
	return 0;
    }

    /* then form L^-1(-L_a'+R_z') in the (1,2) block. */
    col2 = col << 1;
    for (js = col + 1; js <= col2; ++js) {
	dtrsl_(&wn[wn_offset], m2, col, &wn[js * two_m + 1], 11, info);
    }

    /* Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
    /* upper triangle of (2,2) block of wn. */
    for (is = col + 1; is <= col2; ++is) {
	for (js = is; js <= col2; ++js) {
	    wn[is + js * two_m] += 
		ddot_(col, &wn[is * two_m + 1], &wn[js * two_m + 1]);
	}
    }

    /* Cholesky factorization of (2,2) block of wn. */
    dpofa_(&wn[col + 1 + (col + 1) * two_m], m2, col, info);
    if (*info != 0) {
	*info = -2;
	return 0;
    }

    return 0;
} /* formk_ */

static int formt_(int m, double *wt, double *sy, 
		  double *ss, int col, double theta, 
		  int *info)
{
    int i, j, k, k1, wt_offset;
    double ddum;

    /* Parameter adjustments */
    ss -= 1 + m;
    sy -= 1 + m;
    wt_offset = 1 + m;
    wt -= wt_offset;

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

    /* Cholesky factorize T to J*J' with */
    /* J' stored in the upper triangle of wt. */
    dpofa_(&wt[wt_offset], m, col, info);
    if (*info != 0) {
	*info = -3;
    }
    return 0;
} /* formt_ */

static int 
freev_(int n, int *nfree, int *index, 
       int *nenter, int *ileave, int *indx2, int *iwhere, 
       int *wrk, int updatd, int cnstnd, int iter)
{
    int i, k, iact;

    /* Parameter adjustments */
    --iwhere;
    --indx2;
    --index;

    *nenter = 0;
    *ileave = n + 1;

    if (iter > 0 && cnstnd) {
	/* count the entering and leaving variables. */
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

    /* Find the index set of free and active variables at the GCP. */
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
    return 0;
} /* freev_ */

static int dcstep_(double *stx, double *fx, double *dx, 
		   double *sty, double *fy, double *dy, 
		   double *stp, double fp, double dp, 
		   int *brackt, double stpmin, double stpmax)
{
    double d1, d2, d3;
    double p, q, r, s, sgnd, stpc, stpf, stpq, gamma, theta;

    sgnd = dp * (*dx / abs(*dx));

    /* First case: A higher function value. The minimum is bracketed. */
    /* If the cubic step is closer to stx than the quadratic step, the */
    /* cubic step is taken, otherwise the average of the cubic and */
    /* quadratic steps is taken. */

    if (fp > *fx) {
	theta = (*fx - fp) * 3. / (*stp - *stx) + *dx + dp;
	/* Computing MAX */
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(dp);
	s = max(d1, d2);
	/* Computing 2nd power */
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

	/* Second case: A lower function value and derivatives of opposite */
	/* sign. The minimum is bracketed. If the cubic step is farther from */
	/* stp than the secant step, the cubic step is taken, otherwise the */
	/* secant step is taken. */
    } else if (sgnd < 0.) {
	theta = (*fx - fp) * 3. / (*stp - *stx) + *dx + dp;
	/* Computing MAX */
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(dp);
	s = max(d1,d2);
	/* Computing 2nd power */
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

	/* Third case: A lower function value, derivatives of the same sign, */
	/* and the magnitude of the derivative decreases. */
    } else if (abs(dp) < abs(*dx)) {
	/* The cubic step is computed only if the cubic tends to infinity */
	/* in the direction of the step or if the minimum of the cubic */
	/* is beyond stp. Otherwise the cubic step is defined to be the */
	/* secant step. */
	theta = (*fx - fp) * 3. / (*stp - *stx) + *dx + dp;

	/* Computing MAX */
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(dp);
	s = max(d1,d2);

	/* The case gamma = 0 only arises if the cubic does not tend
	   to infinity in the direction of the step. */
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
	    /* A minimizer has been bracketed. If the cubic step is */
	    /* closer to stp than the secant step, the cubic step is */
	    /* taken, otherwise the secant step is taken. */
	    if ((d1 = stpc - *stp, abs(d1)) < (d2 = stpq - *stp, abs(d2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    if (*stp > *stx) {
		/* Computing MIN */
		d1 = *stp + (*sty - *stp) * .66;
		stpf = min(d1,stpf);
	    } else {
		/* Computing MAX */
		d1 = *stp + (*sty - *stp) * .66;
		stpf = max(d1,stpf);
	    }
	} else {
	    /* A minimizer has not been bracketed. If the cubic step is */
	    /* farther from stp than the secant step, the cubic step is */
	    /* taken, otherwise the secant step is taken. */
	    if ((d1 = stpc - *stp, abs(d1)) > (d2 = stpq - *stp, abs(d2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    stpf = min(stpmax, stpf);
	    stpf = max(stpmin, stpf);
	}

	/* Fourth case: A lower function value, derivatives of the same sign, */
	/* and the magnitude of the derivative does not decrease. If the */
	/* minimum is not bracketed, the step is either stpmin or stpmax, */
	/* otherwise the cubic step is taken. */
    } else {
	if (*brackt) {
	    theta = (fp - *fy) * 3. / (*sty - *stp) + *dy + dp;
	    /* Computing MAX */
	    d1 = abs(theta), d2 = abs(*dy), d1 = max(d1,d2), d2 = abs(dp);
	    s = max(d1,d2);
	    /* Computing 2nd power */
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

    /* Update the interval which contains a minimizer. */
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

    /* Compute the new step. */
    *stp = stpf;

    return 0;
} /* dcstep_ */

static int dcsrch_(double f, double g, double *stp, 
		   double ftol, double gtol, double xtol, 
		   double stpmin, double stpmax, 
		   char *task, int *isave, double *dsave)
{
    /* System generated locals */
    double d__1;

    double fm, gm, fx, fy, gx, gy, fxm, fym, gxm, gym, stx, sty;
    double finit, ginit, width, ftest, gtest, stmin, stmax, width1;
    int stage, brackt;

    /* Parameter adjustments */
    --dsave;
    --isave;

    if (!strcmp(task, "START")) {
	/* Check the input arguments for errors. */
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
	/* Exit if there are errors on input. */
	if (!strncmp(task, "ERROR", 5)) {
	    return 1;
	}
	/* Initialize local variables. */
	brackt = 0;
	stage = 1;
	finit = f;
	ginit = g;
	gtest = ftol * ginit;
	width = stpmax - stpmin;
	width1 = width / .5;

	/* The variables stx, fx, gx contain the values of the step,
	   function, and derivative at the best step.
	   The variables sty, fy, gy contain the value of the step,
	   function, and derivative at sty.
	   The variables stp, f, g contain the values of the step,
	   function, and derivative at stp. 
	*/
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
	/* Restore local variables. */
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

    /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
       algorithm enters the second stage. */
    ftest = finit + *stp * gtest;
    if (stage == 1 && f <= ftest && g >= 0.) {
	stage = 2;
    }

    /* Test for warnings */
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

    /* Test for convergence. */
    if (f <= ftest && abs(g) <= gtol * (-ginit)) {
	strcpy(task, "CONVERGENCE");
    }

    /* Test for termination. */
    if (!strncmp(task, "WARN", 4) || !strncmp(task, "CONV", 4)) {
	goto L1000;
    }

    /* A modified function is used to predict the step during the 
       first stage if a lower function value has been obtained but 
       the decrease is not sufficient. */
    if (stage == 1 && f <= fx && f > ftest) {
	/* Define the modified function and derivative values. */
	fm = f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
	/* Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, fm, gm, &brackt,
		stmin, stmax);
	/* Reset the function and derivative values for f. */
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
	/* Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, 
		stmin, stmax);
    }

    /* Decide if a bisection step is needed. */
    if (brackt) {
	if ((d__1 = sty - stx, abs(d__1)) >= width1 * .66) {
	    *stp = stx + (sty - stx) * .5;
	}
	width1 = width;
	width = (d__1 = sty - stx, abs(d__1));
    }

    /* Set the minimum and maximum steps allowed for stp. */
    if (brackt) {
	stmin = min(stx,sty);
	stmax = max(stx,sty);
    } else {
	stmin = *stp + (*stp - stx) * 1.1;
	stmax = *stp + (*stp - stx) * 4.;
    }

    /* Force the step to be within the bounds stpmax and stpmin. */
    *stp = max(*stp, stpmin);
    *stp = min(*stp, stpmax);

    /* If further progress is not possible, let stp be the best */
    /* point obtained during the search. */
    if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin 
	    <= xtol * stmax) {
	*stp = stx;
    }

    /* Obtain another function and derivative. */
    strcpy(task, "FG");

L1000:

    /* Save local variables. */
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
} /* dcsrch_ */

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

    /* Parameter adjustments */
    --z;
    --t;
    --r;
    --d;
    --g;
    --x;
    --nbd;
    --u;
    --l;
    --isave;
    --dsave;

    if (!strncmp(task, "FG_LN", 5)) {
	goto L556;
    }

    *dtd = ddot_(n, &d[1], &d[1]);
    *dnorm = sqrt(*dtd);

    /* Determine the maximum step length */
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
	/* Computing MIN */
	d1 = 1. / *dnorm;
	*stp = min(d1, *stpmx);
    } else {
	*stp = 1.;
    }

    dcopy_(n, &x[1], &t[1]);
    dcopy_(n, &g[1], &r[1]);
    *fold = f;
    *ifun = 0;
    *iback = 0;
    strcpy(csave, "START");

L556:
    *gd = ddot_(n, &g[1], &d[1]);
    if (*ifun == 0) {
	*gdold = *gd;
	if (*gd >= 0.) {
	    /* the directional derivative >=0. 
	       Line search is impossible. */
	    *info = -4;
	    return;
	}
    }

    dcsrch_(f, *gd, stp, FTOL, GTOL, XTOL, STPMIN, *stpmx, csave, 
	    &isave[1], &dsave[1]);

    *xstep = *stp * *dnorm;

    if (strncmp(csave, "CONV", 4) && strncmp(csave, "WARN", 4)) {
	strcpy(task, "FG_LNSRCH");
	++(*ifun);
	++(*nfgv);
	*iback = *ifun - 1;
	if (*stp == 1.) {
	    dcopy_(n, &z[1], &x[1]);
	} else {
	    for (i = 1; i <= n; ++i) {
		x[i] = *stp * d[i] + t[i];
	    }
	}
    } else {
	strcpy(task, "NEW_X");
    }
} /* lnsrlb_ */

static int matupd_(int n, int m, double *ws, 
		   double *wy, double *sy, double *ss, double *d__, 
		   double *r__, int *itail, int iupdat, int *col, 
		   int *head, double *theta, double rr, double dr, 
		   double *stp, double *dtd)
{
    /* System generated locals */
    int ws_offset, wy_offset, sy_offset, ss_offset;

    int j, jmax;
    int pointr;

    /* Parameter adjustments */
    --r__;
    --d__;
    ss_offset = 1 + m;
    ss -= ss_offset;
    sy_offset = 1 + m;
    sy -= sy_offset;
    wy_offset = 1 + n;
    wy -= wy_offset;
    ws_offset = 1 + n;
    ws -= ws_offset;

    if (iupdat <= m) {
	*col = iupdat;
	*itail = (*head + iupdat - 2) % m + 1;
    } else {
	*itail = *itail % m + 1;
	*head = *head % m + 1;
    }

    /* Update matrices WS and WY. */
    dcopy_(n, &d__[1], &ws[*itail * n + 1]);
    dcopy_(n, &r__[1], &wy[*itail * n + 1]);

    /* Set theta = yy/ys. */
    *theta = rr / dr;

    /* Form the middle matrix in B. */
    /* update the upper triangle of SS, */
    /* and the lower triangle of SY: */
    if (iupdat > m) {
	/* move old information */
	jmax = *col - 1;
	for (j = 1; j <= jmax; ++j) {
	    dcopy_(j, &ss[(j + 1) * m + 2], &ss[j * m + 1]);
	    dcopy_(*col - j, &sy[j + 1 + (j + 1) * m], &sy[j + j * m]);
	}
    }

    /* add new information: the last row of SY */
    /* and the last column of SS: */
    pointr = *head;
    jmax = *col - 1;
    for (j = 1; j <= jmax; ++j) {
	sy[*col + j * m] = ddot_(n, &d__[1], &wy[pointr * n + 1]);
	ss[j + *col * m] = ddot_(n, &ws[pointr * n + 1], &d__[1]);
	pointr = pointr % m + 1;
    }
    if (*stp == 1.) {
	ss[*col + *col * m] = *dtd;
    } else {
	ss[*col + *col * m] = *stp * *stp * *dtd;
    }

    sy[*col + *col * m] = dr;

    return 0;
} /* matupd_ */

static int projgr_(int n, double *l, double *u, 
		   int *nbd, double *x, double *g, 
		   double *sbgnrm)
{
    double gi, d1, d2;
    int i;

    /* Parameter adjustments */
    --g;
    --x;
    --nbd;
    --u;
    --l;

    *sbgnrm = 0.;

    for (i = 1; i <= n; ++i) {
	gi = g[i];

	if (nbd[i] != 0) {
	    if (gi < 0.) {
		if (nbd[i] >= 2) {
		    /* Computing MAX */
		    d1 = x[i] - u[i];
		    gi = max(d1, gi);
		}
	    } else {
		if (nbd[i] <= 2) {
		    /* Computing MIN */
		    d1 = x[i] - l[i];
		    gi = min(d1, gi);
		}
	    }
	}

	/* Computing MAX */
	d1 = *sbgnrm, d2 = abs(gi);
	*sbgnrm = max(d1, d2);
    }
    return 0;
} /* projgr_ */

static int 
subsm_(int n, int m, int nsub, int *ind, 
       double *l, double *u, int *nbd, double *x, 
       double *d__, double *ws, double *wy, double theta, 
       int col, int head, int *iword, double *wv, 
       double *wn, int *info)
{
    /* System generated locals */
    int ws_offset, wy_offset, wn_offset;

    int i, j, k, m2;
    double dk;
    int js, jy, col2, ibd = 0;
    double temp1, temp2, alpha;
    int pointr;

    /* Parameter adjustments */
    --d__;
    --x;
    --nbd;
    --u;
    --l;
    wn_offset = 1 + 2 * m;
    wn -= wn_offset;
    --wv;
    wy_offset = 1 + n;
    wy -= wy_offset;
    ws_offset = 1 + n;
    ws -= ws_offset;
    --ind;

    if (nsub <= 0) {
	return 0;
    }

    /* Compute wv = W'Zd. */
    pointr = head;
    for (i = 1; i <= col; ++i) {
	temp1 = 0.;
	temp2 = 0.;
	for (j = 1; j <= nsub; ++j) {
	    k = ind[j];
	    temp1 += wy[k + pointr * n] * d__[j];
	    temp2 += ws[k + pointr * n] * d__[j];
	}
	wv[i] = temp1;
	wv[col + i] = theta * temp2;
	pointr = pointr % m + 1;
    }

    /* Compute wv:=K^(-1)wv. */
    m2 = m << 1;
    col2 = col << 1;
    dtrsl_(&wn[wn_offset], m2, col2, &wv[1], 11, info);
    if (*info != 0) {
	return 0;
    }
    for (i = 1; i <= col; ++i) {
	wv[i] = -wv[i];
    }
    dtrsl_(&wn[wn_offset], m2, col2, &wv[1], 1, info);
    if (*info != 0) {
	return 0;
    }

    /* Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
    pointr = head;
    for (jy = 1; jy <= col; ++jy) {
	js = col + jy;
	for (i = 1; i <= nsub; ++i) {
	    k = ind[i];
	    d__[i] = d__[i] + wy[k + pointr * n] * wv[jy] / theta 
		    + ws[k + pointr * n] * wv[js];
	}
	pointr = pointr % m + 1;
    }
    for (i = 1; i <= nsub; ++i) {
	d__[i] /= theta;
    }

    /* Backtrack to the feasible region. */
    alpha = 1.;
    temp1 = alpha;
    for (i = 1; i <= nsub; ++i) {
	k = ind[i];
	dk = d__[i];
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
	dk = d__[ibd];
	k = ind[ibd];
	if (dk > 0.) {
	    x[k] = u[k];
	    d__[ibd] = 0.;
	} else if (dk < 0.) {
	    x[k] = l[k];
	    d__[ibd] = 0.;
	}
    }

    for (i = 1; i <= nsub; ++i) {
	k = ind[i];
	x[k] += alpha * d__[i];
    }

    if (alpha < 1.) {
	*iword = 1;
    } else {
	*iword = 0;
    }
    return 0;
} /* subsm_ */

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

    /* determine ibeta, beta ala malcolm. */
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
    /* determine it, irnd. */
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
    /* determine dpmeps. */
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
} /* dpmeps_ */

static int mainlb_(int n, int m, double *x, 
		   double *l, double *u, int *nbd, double *f, double *g,
		   double factr, double pgtol, double *ws, double *wy,
		   double *sy, double *ss, double *yy, double *wt, 
		   double *wn, double *snd, double *z, double *r, 
		   double *d, double *t, double *wa, double *sg, 
		   double *sgo, double *yg, double *ygo, int *index, 
		   int *iwhere, int *indx2, char *task, char *csave, 
		   int *lsave, int *isave, double *dsave)
{
    /* System generated locals */
    int ws_offset, wy_offset, sy_offset;
    int ss_offset, yy_offset, wt_offset;
    int wn_offset, snd_offset;
    double d1, d2;

    int i, k;
    double gd, dr, rr, dtd, tol;
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
    int boxed;
    int itail = 0, iword = 0, ileave = 0;
    int nskip;
    int itfile = 0;
    int updatd, iupdat;
    int prjctd, cnstnd;
    int nenter = 0;
    int nintol;

    /* Parameter adjustments */
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

    snd_offset = 1 + 2 * m;
    snd -= snd_offset;
    wn_offset = 1 + 2 * m;
    wn -= wn_offset;
    wt_offset = 1 + m;
    wt -= wt_offset;
    yy_offset = 1 + m;
    yy -= yy_offset;
    ss_offset = 1 + m;
    ss -= ss_offset;
    sy_offset = 1 + m;
    sy -= sy_offset;
    wy_offset = 1 + n;
    wy -= wy_offset;
    ws_offset = 1 + n;
    ws -= ws_offset;

    --lsave;
    --isave;
    --dsave;

    if (!strcmp(task, "START")) {
	timer_(&time1);

	/* Generate the current machine precision */
	epsmch = dpmeps_();

	/* Initialize counters and scalars when task='START' */

	/* limited memory BFGS matrices */
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;

	/* operation counts */
	iter = 0;
	nfgv = 0;
	nint = 0;
	nintol = 0;
	nskip = 0;
	nfree = n;

	/* stopping tolerance */
	tol = factr * epsmch;

	/* measuring running time */
	cachyt = 0.;
	sbtime = 0.;
	lnscht = 0.;

	/* 'info' records the termination information */
	info = 0;

	/* Check the input arguments for errors */
	errclb_(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k);
	if (!strncmp(task, "ERROR", 5)) {
	    return E_DATA;
	}

	/* Initialize iwhere and project x onto the feasible set */
	active_(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], &prjctd, 
		&cnstnd, &boxed);
    } else {
	/* restore local variables. */
	prjctd = lsave[1];
	cnstnd = lsave[2];
	boxed = lsave[3];
	updatd = lsave[4];
	nintol = isave[1];
	itfile = isave[3];
	iback = isave[4];
	nskip = isave[5];
	head = isave[6];
	col = isave[7];
	itail = isave[8];
	iter = isave[9];
	iupdat = isave[10];
	nint = isave[12];
	nfgv = isave[13];
	info = isave[14];
	ifun = isave[15];
	iword = isave[16];
	nfree = isave[17];
	nact = isave[18];
	ileave = isave[19];
	nenter = isave[20];
	theta = dsave[1];
	fold = dsave[2];
	tol = dsave[3];
	dnorm = dsave[4];
	epsmch = dsave[5];
	cpu1 = dsave[6];
	cachyt = dsave[7];
	sbtime = dsave[8];
	lnscht = dsave[9];
	time1 = dsave[10];
	gd = dsave[11];
	stpmx = dsave[12];
	sbgnrm = dsave[13];
	stp = dsave[14];
	gdold = dsave[15];
	dtd = dsave[16];

	/* After returning from the driver go to the point where execution
	   is to resume. */
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

    /* Compute f0 and g0 */
    strcpy(task, "FG_START");

    /* return to the driver to calculate f and g; reenter at 111 */
    goto L1000;

L111:
    nfgv = 1;
    /* Compute the infinity norm of the (-) projected gradient. */
    projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);

    if (sbgnrm <= pgtol) {
	/* terminate the algorithm. */
	strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }

/* ----------------- the beginning of the loop -------------------------- */

L222:
    iword = -1;

    if (!cnstnd && col > 0) {
	/* skip the search for GCP */
	dcopy_(n, &x[1], &z[1]);
	wrk = updatd;
	nint = 0;
	goto L333;
    }

    /* Compute the Generalized Cauchy Point (GCP). */

    timer_(&cpu1);

    cauchy_(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], 
	    &t[1], &d[1], &z[1], m, &wy[wy_offset], &ws[ws_offset], 
	    &sy[sy_offset], &wt[wt_offset], theta, col, head, &wa[1], 
	    &wa[(m << 1) + 1], &wa[(m << 2) + 1], &wa[m * 6 + 1], &nint, 
	    &sg[1], &yg[1], sbgnrm, &info, epsmch);

    if (info != 0) {
	/* singular triangular system detected; refresh the lbfgs memory. */
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

    /* Count the entering and leaving variables for iter > 0; */
    /* find the index set of free and active variables at the GCP. */
    freev_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iwhere[1], 
	   &wrk, updatd, cnstnd, iter);
    nact = n - nfree;

L333:
    /* If there are no free variables or B=theta*I, then */
    /* skip the subspace minimization. */
    if (nfree == 0 || col == 0) {
	goto L555;
    }

    /* Subspace minimization. */

    timer_(&cpu1);

    /*     Form  the LEL^T factorization of the indefinite */
    /*       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
    /*                     [L_a -R_z           theta*S'AA'S ] */
    /*       where     E = [-I  0] */
    /*                     [ 0  I] */
    if (wrk) {
	formk_(n, nfree, &index[1], nenter, ileave, &indx2[1], iupdat, 
	       updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset],
	       &wy[wy_offset], &sy[sy_offset], theta, col, head, &info);
    }

    if (info != 0) {
	/* nonpositive definiteness in Cholesky factorization; 
	   refresh the lbfgs memory and restart the iteration. */
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

    /* compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) 
       from 'cauchy'). */
    cmprlb_(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], 
	    &sy[sy_offset], &wt[wt_offset], &z[1], &r[1], &wa[1], 
	    &index[1], theta, col, head, nfree, cnstnd, &info);
    if (info != 0) {
	goto L444;
    }

    /* call the direct method */
    subsm_(n, m, nfree, &index[1], &l[1], &u[1], &nbd[1], &z[1], &r[1], 
	   &ws[ws_offset], &wy[wy_offset], theta, col, head, &iword, 
	   &wa[1], &wn[wn_offset], &info);

 L444:
    if (info != 0) {
	/* singular triangular system detected; refresh the lbfgs
	   memory and restart the iteration.
	*/
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

    /* Line search and optimality tests */

    /* Generate the search direction d:=z-x */
    for (i = 1; i <= n; ++i) {
	d[i] = z[i] - x[i];
    }

    timer_(&cpu1);

L666:
    lnsrlb_(n, &l[1], &u[1], &nbd[1], &x[1], *f, &fold, &gd, &gdold, &g[1],
	    &d[1], &r[1], &t[1], &z[1], &stp, &dnorm, &dtd, &xstep, 
	    &stpmx, iter, &ifun, &iback, &nfgv, &info, task, boxed, cnstnd, 
	    csave, &isave[22], &dsave[17]);

    if (info != 0 || iback >= 20) {
	/* restore the previous iterate */
	dcopy_(n, &t[1], &x[1]);
	dcopy_(n, &r[1], &g[1]);
	*f = fold;
	if (col == 0) {
	    /* abnormal termination */
	    if (info == 0) {
		info = -9;
		/* restore the actual number of f and g evaluations etc. */
		--nfgv;
		--ifun;
		--iback;
	    }
	    strcpy(task, "ABNORMAL_TERMINATION_IN_LNSRCH");
	    ++iter;
	    goto L999;
	} else {
	    /* refresh the lbfgs memory and restart the iteration */
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
	/* return to the driver for calculating f and g; reenter at 666 */
	goto L1000;
    } else {
	/* calculate the quantities related to the new X */
	timer_(&cpu2);
	lnscht = lnscht + cpu2 - cpu1;
	++iter;
	/* Compute the infinity norm of the projected (-)gradient */
	projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
	goto L1000;
    }

L777:

    /* Test for termination. */
    if (sbgnrm <= pgtol) {
	/* terminate the algorithm. */
	strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }

    /* Computing MAX */
    d1 = abs(fold), d2 = abs(*f), d1 = max(d1,d2);
    ddum = max(d1,1.);
    if (fold - *f <= tol * ddum) {
	/* terminate the algorithm. */
	strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
	if (iback >= 10) {
	    info = -5;
	}
	/* i.e., to issue a warning if iback>10 in the line search. */
	goto L999;
    }

    /* Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's. */
    for (i = 1; i <= n; ++i) {
	r[i] = g[i] - r[i];
    }
    rr = ddot_(n, &r[1], &r[1]);

    if (stp == 1.) {
	dr = gd - gdold;
	ddum = -gdold;
    } else {
	dr = (gd - gdold) * stp;
	dscal_(n, stp, &d[1]);
	ddum = -gdold * stp;
    }

    if (dr <= epsmch * ddum) {
	/* skip the L-BFGS update. */
	++nskip;
	updatd = 0;
	goto L888;
    }

    /* Update the L-BFGS matrix. */

    updatd = 1;
    ++iupdat;

    /* Update matrices WS and WY and form the middle matrix in B. */
    matupd_(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], 
	    &ss[ss_offset], &d[1], &r[1], &itail, iupdat, &col, 
	    &head, &theta, rr, dr, &stp, &dtd);

    /* Form the upper half of the pds T = theta*SS + L*D^(-1)*L'; 
       Store T in the upper triangular of the array wt;
       Cholesky factorize T to J*J' with
       J' stored in the upper triangular of wt. 
    */

    formt_(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], 
	   col, theta, &info);
    if (info != 0) {
	/* nonpositive definiteness in Cholesky factorization;
	   refresh the lbfgs memory and restart the iteration. */
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	goto L222;
    }
    /*     Now the inverse of the middle matrix in B is */
    /*       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ] */
    /*       [ -L*D^(-1/2)   J ] [  0        J'          ] */

L888:

/* -------------------- the end of the loop ----------------------------- */
    goto L222;

L999:
    timer_(&time2);
    time = time2 - time1;

L1000:
    /* Save local variables. */
    lsave[1] = prjctd;
    lsave[2] = cnstnd;
    lsave[3] = boxed;
    lsave[4] = updatd;
    isave[1] = nintol;
    isave[3] = itfile;
    isave[4] = iback;
    isave[5] = nskip;
    isave[6] = head;
    isave[7] = col;
    isave[8] = itail;
    isave[9] = iter;
    isave[10] = iupdat;
    isave[12] = nint;
    isave[13] = nfgv;
    isave[14] = info;
    isave[15] = ifun;
    isave[16] = iword;
    isave[17] = nfree;
    isave[18] = nact;
    isave[19] = ileave;
    isave[20] = nenter;
    dsave[1] = theta;
    dsave[2] = fold;
    dsave[3] = tol;
    dsave[4] = dnorm;
    dsave[5] = epsmch;
    dsave[6] = cpu1;
    dsave[7] = cachyt;
    dsave[8] = sbtime;
    dsave[9] = lnscht;
    dsave[10] = time1;
    dsave[11] = gd;
    dsave[12] = stpmx;
    dsave[13] = sbgnrm;
    dsave[14] = stp;
    dsave[15] = gdold;
    dsave[16] = dtd;

    return 0;
} /* mainlb_ */

static int setulb_(int n, int m, double *x, 
		   double *l, double *u, int *nbd, double *f, double *g,
		   double factr, double pgtol, double *wa, int *iwa,
		   char *task, char *csave, int *lsave, 
		   int *isave, double *dsave)
{
    int ld, lr, lt, lz, lwa, lsg, lyg, lwn, lss, lws, 
	lwt, lsy, lwy, lyy, lsnd, lsgo, lygo;

    /* Parameter adjustments */
    --iwa;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --wa;
    --lsave;
    --isave;
    --dsave;

    if (!strcmp(task, "START")) {
	isave[1] = m * n;
	isave[2] = m * m;
	isave[3] = m * m << 2;
	isave[4] = 1;
	isave[5] = isave[4] + isave[1];
	isave[6] = isave[5] + isave[1];
	isave[7] = isave[6] + isave[2];
	isave[8] = isave[7] + isave[2];
	isave[9] = isave[8] + isave[2];
	isave[10] = isave[9] + isave[2];
	isave[11] = isave[10] + isave[3];
	isave[12] = isave[11] + isave[3];
	isave[13] = isave[12] + n;
	isave[14] = isave[13] + n;
	isave[15] = isave[14] + n;
	isave[16] = isave[15] + n;
	isave[17] = isave[16] + (m << 3);
	isave[18] = isave[17] + m;
	isave[19] = isave[18] + m;
	isave[20] = isave[19] + m;
    }

    lws = isave[4];
    lwy = isave[5];
    lsy = isave[6];
    lss = isave[7];
    lyy = isave[8];
    lwt = isave[9];
    lwn = isave[10];
    lsnd = isave[11];
    lz = isave[12];
    lr = isave[13];
    ld = isave[14];
    lt = isave[15];
    lwa = isave[16];
    lsg = isave[17];
    lsgo = isave[18];
    lyg = isave[19];
    lygo = isave[20];

    int err;

    err = mainlb_(n, m, &x[1], &l[1], &u[1], &nbd[1], f, &g[1], factr, pgtol, 
		  &wa[lws], &wa[lwy], &wa[lsy], &wa[lss], &wa[lyy], &wa[lwt], 
		  &wa[lwn], &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], 
		  &wa[lwa], &wa[lsg], &wa[lsgo], &wa[lyg], &wa[lygo], 
		  &iwa[1], &iwa[n + 1], &iwa[(n << 1) + 1], task, csave, 
		  &lsave[1], &isave[22], &dsave[1]);

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

    /* the number of limited memory corrections stored:
       3 < m < 15 is recommended */
    m = 5;

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

    /* Convergence criteria: we don't use "reltol" because I haven't
       yet figured out exactly how it relates to the stopping
       criteria used here. */
    factr = 100.;
    pgtol = 0.; /* should set this here? */

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
	setulb_(n, m, b, l, u, nbd, &f, g, factr, pgtol, wa, iwa, task,
		csave, lsave, isave, dsave);

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
	    for (i=0; i<n; i++) {
		g[i] = -g[i];
	    }	
	    *grcount += 1;
	} else if (!strncmp(task, "NEW_X", 5)) {
	    /* The optimizer has produced a new set of parameter values */

	    if (isave[33] >= maxit) {
		strcpy(task, "STOP: TOTAL NO. of f AND g "
		       "EVALUATIONS EXCEEDS LIMIT");
		err = E_NOCONV;
	    } else if (dsave[12] <= (abs(f) + 1.) * 1e-10) {
		/* Terminate if  |proj g|/(1+|f|) < 1.0e-10, where
		   "proj g" denoted the projected gradient -- or could
		   rely on built-in mechanism? */
		strcpy(task, "STOP: THE PROJECTED GRADIENT IS "
		       "SUFFICIENTLY SMALL");
	    } else if (opt & OPT_V) {
		reverse_gradient(g, n);
		print_iter_info(isave[29], -f, crittype, n, b, g, dsave[13], prn);
		reverse_gradient(g, n);
	    }
	} else {
	    /* ?? study up on other possible "task" contents! */
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
