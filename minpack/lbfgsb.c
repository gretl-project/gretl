/* 
   routines.f -- translated by f2c (version 20061008).

   See http://www.ece.northwestern.edu/~nocedal/lbfgsb.html
*/

#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef min
# define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
# define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

static int c1 = 1;
static int c11 = 11;
static double c_b173 = .001;
static double c_b174 = .9;
static double c_b175 = .1;
static double c_b176 = 1.0e-15; /* min. step value: was 0. */

/* Builtin functions */
double sqrt(double);

static int dcopy_(int *n, double *dx, int *incx, 
		  double *dy, int *incy)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    int i, m, ix, iy, mp1;

    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
    }
    return 0;
 L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i1 = m;
    for (i = 1; i <= i1; ++i) {
	dy[i] = dx[i];
    }
    if (*n < 7) {
	return 0;
    }
 L40:
    mp1 = m + 1;
    i1 = *n;
    for (i = mp1; i <= i1; i += 7) {
	dy[i] = dx[i];
	dy[i + 1] = dx[i + 1];
	dy[i + 2] = dx[i + 2];
	dy[i + 3] = dx[i + 3];
	dy[i + 4] = dx[i + 4];
	dy[i + 5] = dx[i + 5];
	dy[i + 6] = dx[i + 6];
    }
    return 0;
} /* dcopy_ */

static double ddot_(int *n, double *dx, int *incx, double *dy, 
		    int *incy)
{
    /* System generated locals */
    int i1;
    double ret_val;

    /* Local variables */
    int i, m, ix, iy, mp1;
    double dtemp;

    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
    }
    ret_val = dtemp;
    return ret_val;
 L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i1 = m;
    for (i = 1; i <= i1; ++i) {
	dtemp += dx[i] * dy[i];
    }
    if (*n < 5) {
	goto L60;
    }
 L40:
    mp1 = m + 1;
    i1 = *n;
    for (i = mp1; i <= i1; i += 5) {
	dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[
								   i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 
														   4] * dy[i + 4];
    }
 L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

static int daxpy_(int *n, double *da, double *dx, 
		  int *incx, double *dy, int *incy)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    int i, m, ix, iy, mp1;

    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
    }
    return 0;
 L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i1 = m;
    for (i = 1; i <= i1; ++i) {
	dy[i] += *da * dx[i];
    }
    if (*n < 4) {
	return 0;
    }
 L40:
    mp1 = m + 1;
    i1 = *n;
    for (i = mp1; i <= i1; i += 4) {
	dy[i] += *da * dx[i];
	dy[i + 1] += *da * dx[i + 1];
	dy[i + 2] += *da * dx[i + 2];
	dy[i + 3] += *da * dx[i + 3];
    }
    return 0;
} /* daxpy_ */

static int dscal_(int *n, double *da, double *dx, int *incx)
{
    /* System generated locals */
    int i1, i2;

    /* Local variables */
    int i, m, mp1, nincx;

    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }
    nincx = *n * *incx;
    i1 = nincx;
    i2 = *incx;
    for (i = 1; i2 < 0 ? i >= i1 : i <= i1; i += i2) {
	dx[i] = *da * dx[i];
    }
    return 0;
 L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i2 = m;
    for (i = 1; i <= i2; ++i) {
	dx[i] = *da * dx[i];
    }
    if (*n < 5) {
	return 0;
    }
 L40:
    mp1 = m + 1;
    i2 = *n;
    for (i = mp1; i <= i2; i += 5) {
	dx[i] = *da * dx[i];
	dx[i + 1] = *da * dx[i + 1];
	dx[i + 2] = *da * dx[i + 2];
	dx[i + 3] = *da * dx[i + 3];
	dx[i + 4] = *da * dx[i + 4];
    }
    return 0;
} /* dscal_ */

static int dtrsl_(double *t, int *ldt, int *n, 
		  double *b, int *job, int *info)
{
    /* System generated locals */
    int t_dim1, t_offset, i1, i2;

    /* Local variables */
    int j, jj, case__;
    double temp;

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i1 = *n;
    for (*info = 1; *info <= i1; ++(*info)) {
	if (t[*info + *info * t_dim1] == 0.) {
	    goto L150;
	}
    }
    *info = 0;
    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch (case__) {
    case 1:  goto L20;
    case 2:  goto L50;
    case 3:  goto L80;
    case 4:  goto L110;
    }
 L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L40;
    }
    i1 = *n;
    for (j = 2; j <= i1; ++j) {
	temp = -b[j - 1];
	i2 = *n - j + 1;
	daxpy_(&i2, &temp, &t[j + (j - 1) * t_dim1], &c1, &b[j], &c1);
	b[j] /= t[j + j * t_dim1];
    }
 L40:
    goto L140;
 L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L70;
    }
    i1 = *n;
    for (jj = 2; jj <= i1; ++jj) {
	j = *n - jj + 1;
	temp = -b[j + 1];
	daxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c1, &b[1], &c1);
	b[j] /= t[j + j * t_dim1];
    }
 L70:
    goto L140;
 L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L100;
    }
    i1 = *n;
    for (jj = 2; jj <= i1; ++jj) {
	j = *n - jj + 1;
	i2 = jj - 1;
	b[j] -= ddot_(&i2, &t[j + 1 + j * t_dim1], &c1, &b[j + 1], &c1);
	b[j] /= t[j + j * t_dim1];
    }
 L100:
    goto L140;
 L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L130;
    }
    i1 = *n;
    for (j = 2; j <= i1; ++j) {
	i2 = j - 1;
	b[j] -= ddot_(&i2, &t[j * t_dim1 + 1], &c1, &b[1], &c1);
	b[j] /= t[j + j * t_dim1];
    }
 L130:
 L140:
 L150:
    return 0;
} /* dtrsl_ */

static int dpofa_(double *a, int *lda, int *n, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i1, i2, i3;

    /* Local variables */
    int j, k;
    double s, t;
    int jm1;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i1 = *n;
    for (j = 1; j <= i1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i2 = jm1;
	for (k = 1; k <= i2; ++k) {
	    i3 = k - 1;
	    t = a[k + j * a_dim1] - ddot_(&i3, &a[k * a_dim1 + 1], &c1, &
					  a[j * a_dim1 + 1], &c1);
	    t /= a[k + k * a_dim1];
	    a[k + j * a_dim1] = t;
	    s += t * t;
	}
    L20:
	s = a[j + j * a_dim1] - s;
	if (s <= 0.) {
	    goto L40;
	}
	a[j + j * a_dim1] = sqrt(s);
    }
    *info = 0;
 L40:
    return 0;
} /* dpofa_ */

static int active_(int *n, double *l, double *u, 
		   int *nbd, double *x, int *iwhere, 
		   int *prjctd, int *cnstnd, int *boxed)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    int i, nbdd;

    /* Parameter adjustments */
    --iwhere;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    nbdd = 0;
    *prjctd = 0;
    *cnstnd = 0;
    *boxed = 1;
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
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
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
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
    return 0;
} /* active_ */

static int bmv_(int *m, double *sy, double *wt, int *col, 
		double *v, double *p, int *info)
{
    /* System generated locals */
    int sy_dim1, sy_offset, wt_dim1, wt_offset, i1, i2;

    /* Local variables */
    int i, k, icol;
    double sum;

    /* Parameter adjustments */
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    --p;
    --v;

    /* Function Body */
    if (*col == 0) {
	return 0;
    }
    p[*col + 1] = v[*col + 1];
    i1 = *col;
    for (i = 2; i <= i1; ++i) {
	icol = *col + i;
	sum = 0.;
	i2 = i - 1;
	for (k = 1; k <= i2; ++k) {
	    sum += sy[i + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
	}
	p[icol] = v[icol] + sum;
    }
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c11, info);
    if (*info != 0) {
	return 0;
    }
    i1 = *col;
    for (i = 1; i <= i1; ++i) {
	p[i] = v[i] / sqrt(sy[i + i * sy_dim1]);
    }
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c1, info);
    if (*info != 0) {
	return 0;
    }
    i1 = *col;
    for (i = 1; i <= i1; ++i) {
	p[i] = -p[i] / sqrt(sy[i + i * sy_dim1]);
    }
    i1 = *col;
    for (i = 1; i <= i1; ++i) {
	sum = 0.;
	i2 = *col;
	for (k = i + 1; k <= i2; ++k) {
	    sum += sy[k + i * sy_dim1] * p[*col + k] / sy[i + i * 
							  sy_dim1];
	}
	p[i] += sum;
    }
    return 0;
} /* bmv_ */

static int cmprlb_(int *n, int *m, double *x, 
		   double *g, double *ws, double *wy, double *sy, 
		   double *wt, double *z, double *r, double *wa, 
		   int *index, double *theta, int *col, int *head, 
		   int *nfree, int *cnstnd, int *info)
{
    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	wt_dim1, wt_offset, i1, i2;

    /* Local variables */
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
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;

    /* Function Body */
    if (! (*cnstnd) && *col > 0) {
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	    r[i] = -g[i];
	}
    } else {
	i1 = *nfree;
	for (i = 1; i <= i1; ++i) {
	    k = index[i];
	    r[i] = -(*theta) * (z[k] - x[k]) - g[k];
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wa[(*m << 1) + 1], &wa[
									     1], info);
	if (*info != 0) {
	    *info = -8;
	    return 0;
	}
	pointr = *head;
	i1 = *col;
	for (j = 1; j <= i1; ++j) {
	    a1 = wa[j];
	    a2 = *theta * wa[*col + j];
	    i2 = *nfree;
	    for (i = 1; i <= i2; ++i) {
		k = index[i];
		r[i] = r[i] + wy[k + pointr * wy_dim1] * a1 + ws[k + 
								 pointr * ws_dim1] * a2;
	    }
	    pointr = pointr % *m + 1;
	}
    }
    return 0;
} /* cmprlb_ */

static int errclb_(int *n, int *m, double *factr, 
		   double *l, double *u, int *nbd, char *task, int *info,
		   int *k)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    int i;

    /* Parameter adjustments */
    --nbd;
    --u;
    --l;

    /* Function Body */
    if (*n <= 0) {
	strcpy(task, "ERROR: N .LE. 0");
    }
    if (*m <= 0) {
	strcpy(task, "ERROR: M .LE. 0");
    }
    if (*factr < 0.) {
	strcpy(task, "ERROR: FACTR .LT. 0");
    }
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	if (nbd[i] < 0 || nbd[i] > 3) {
	    strcpy(task, "ERROR: INVALID NBD");
	    *info = -6;
	    *k = i;
	}
	if (nbd[i] == 2) {
	    if (l[i] > u[i]) {
		strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
		*info = -7;
		*k = i;
	    }
	}
    }
    return 0;
} /* errclb_ */

static int formk_(int *n, int *nsub, int *ind, int *
		  nenter, int *ileave, int *indx2, int *iupdat, int *
		  updatd, double *wn, double *wn1, int *m, double *ws, 
		  double *wy, double *sy, double *theta, int *col, 
		  int *head, int *info)
{
    /* System generated locals */
    int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset, 
	wy_dim1, wy_offset, sy_dim1, sy_offset, i1, i2, i3;

    /* Local variables */
    int i, k, k1, m2, is, js, iy, jy, is1, js1, col2, dend, pend;
    int upcl;
    double temp1, temp2, temp3, temp4;
    int ipntr, jpntr, dbegin, pbegin;

    /* Parameter adjustments */
    --indx2;
    --ind;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    wn1_dim1 = 2 * *m;
    wn1_offset = 1 + wn1_dim1;
    wn1 -= wn1_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;

    /* Function Body */
    if (*updatd) {
	if (*iupdat > *m) {
	    i1 = *m - 1;
	    for (jy = 1; jy <= i1; ++jy) {
		js = *m + jy;
		i2 = *m - jy;
		dcopy_(&i2, &wn1[jy + 1 + (jy + 1) * wn1_dim1], &c1, &wn1[
									  jy + jy * wn1_dim1], &c1);
		i2 = *m - jy;
		dcopy_(&i2, &wn1[js + 1 + (js + 1) * wn1_dim1], &c1, &wn1[
									  js + js * wn1_dim1], &c1);
		i2 = *m - 1;
		dcopy_(&i2, &wn1[*m + 2 + (jy + 1) * wn1_dim1], &c1, &wn1[
									  *m + 1 + jy * wn1_dim1], &c1);
	    }
	}
	pbegin = 1;
	pend = *nsub;
	dbegin = *nsub + 1;
	dend = *n;
	iy = *col;
	is = *m + *col;
	ipntr = *head + *col - 1;
	if (ipntr > *m) {
	    ipntr -= *m;
	}
	jpntr = *head;
	i1 = *col;
	for (jy = 1; jy <= i1; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    i2 = pend;
	    for (k = pbegin; k <= i2; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    i2 = dend;
	    for (k = dbegin; k <= i2; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    wn1[iy + jy * wn1_dim1] = temp1;
	    wn1[is + js * wn1_dim1] = temp2;
	    wn1[is + jy * wn1_dim1] = temp3;
	    jpntr = jpntr % *m + 1;
	}
	jy = *col;
	jpntr = *head + *col - 1;
	if (jpntr > *m) {
	    jpntr -= *m;
	}
	ipntr = *head;
	i1 = *col;
	for (i = 1; i <= i1; ++i) {
	    is = *m + i;
	    temp3 = 0.;
	    i2 = pend;
	    for (k = pbegin; k <= i2; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    ipntr = ipntr % *m + 1;
	    wn1[is + jy * wn1_dim1] = temp3;
	}
	upcl = *col - 1;
    } else {
	upcl = *col;
    }
    ipntr = *head;
    i1 = upcl;
    for (iy = 1; iy <= i1; ++iy) {
	is = *m + iy;
	jpntr = *head;
	i2 = iy;
	for (jy = 1; jy <= i2; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    temp4 = 0.;
	    i3 = *nenter;
	    for (k = 1; k <= i3; ++k) {
		k1 = indx2[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
	    }
	    i3 = *n;
	    for (k = *ileave; k <= i3; ++k) {
		k1 = indx2[k];
		temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
	    }
	    wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
	    wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
	    jpntr = jpntr % *m + 1;
	}
	ipntr = ipntr % *m + 1;
    }
    ipntr = *head;
    i1 = *m + upcl;
    for (is = *m + 1; is <= i1; ++is) {
	jpntr = *head;
	i2 = upcl;
	for (jy = 1; jy <= i2; ++jy) {
	    temp1 = 0.;
	    temp3 = 0.;
	    i3 = *nenter;
	    for (k = 1; k <= i3; ++k) {
		k1 = indx2[k];
		temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    i3 = *n;
	    for (k = *ileave; k <= i3; ++k) {
		k1 = indx2[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    if (is <= jy + *m) {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] + temp1 - 
		    temp3;
	    } else {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] - temp1 + 
		    temp3;
	    }
	    jpntr = jpntr % *m + 1;
	}
	ipntr = ipntr % *m + 1;
    }
    m2 = *m << 1;
    i1 = *col;
    for (iy = 1; iy <= i1; ++iy) {
	is = *col + iy;
	is1 = *m + iy;
	i2 = iy;
	for (jy = 1; jy <= i2; ++jy) {
	    js = *col + jy;
	    js1 = *m + jy;
	    wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
	    wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
	}
	i2 = iy - 1;
	for (jy = 1; jy <= i2; ++jy) {
	    wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
	}
	i2 = *col;
	for (jy = iy; jy <= i2; ++jy) {
	    wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
	}
	wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
    }
    dpofa_(&wn[wn_offset], &m2, col, info);
    if (*info != 0) {
	*info = -1;
	return 0;
    }
    col2 = *col << 1;
    i1 = col2;
    for (js = *col + 1; js <= i1; ++js) {
	dtrsl_(&wn[wn_offset], &m2, col, &wn[js * wn_dim1 + 1], &c11, info);
    }
    i1 = col2;
    for (is = *col + 1; is <= i1; ++is) {
	i2 = col2;
	for (js = is; js <= i2; ++js) {
	    wn[is + js * wn_dim1] += ddot_(col, &wn[is * wn_dim1 + 1], &c1, 
					   &wn[js * wn_dim1 + 1], &c1);
	}
    }
    dpofa_(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
    if (*info != 0) {
	*info = -2;
	return 0;
    }
    return 0;
} /* formk_ */

static int formt_(int *m, double *wt, double *sy, 
		  double *ss, int *col, double *theta, int *info)
{
    /* System generated locals */
    int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i1, 
	i2, i3;

    /* Local variables */
    int i, j, k, k1;
    double ddum;

    /* Parameter adjustments */
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;

    /* Function Body */
    i1 = *col;
    for (j = 1; j <= i1; ++j) {
	wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
    }
    i1 = *col;
    for (i = 2; i <= i1; ++i) {
	i2 = *col;
	for (j = i; j <= i2; ++j) {
	    k1 = min(i,j) - 1;
	    ddum = 0.;
	    i3 = k1;
	    for (k = 1; k <= i3; ++k) {
		ddum += sy[i + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k + 
								       k * sy_dim1];
	    }
	    wt[i + j * wt_dim1] = ddum + *theta * ss[i + j * ss_dim1];
	}
    }
    dpofa_(&wt[wt_offset], m, col, info);
    if (*info != 0) {
	*info = -3;
    }
    return 0;
} /* formt_ */

static int freev_(int *n, int *nfree, int *index, 
		  int *nenter, int *ileave, int *indx2, int *iwhere, 
		  int *wrk, int *updatd, int *cnstnd, 
		  int *iter)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    int i, k, iact;

    /* Parameter adjustments */
    --iwhere;
    --indx2;
    --index;

    /* Function Body */
    *nenter = 0;
    *ileave = *n + 1;
    if (*iter > 0 && *cnstnd) {
	i1 = *nfree;
	for (i = 1; i <= i1; ++i) {
	    k = index[i];
	    if (iwhere[k] > 0) {
		--(*ileave);
		indx2[*ileave] = k;
	    }
	}
	i1 = *n;
	for (i = *nfree + 1; i <= i1; ++i) {
	    k = index[i];
	    if (iwhere[k] <= 0) {
		++(*nenter);
		indx2[*nenter] = k;
	    }
	}
    }
    *wrk = *ileave < *n + 1 || *nenter > 0 || *updatd;
    *nfree = 0;
    iact = *n + 1;
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
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

static int hpsolb_(int *n, double *t, int *iorder, 
		   int *iheap)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    int i, j, k;
    double out, ddum;
    int indxin, indxou;

    /* Parameter adjustments */
    --iorder;
    --t;

    /* Function Body */
    if (*iheap == 0) {
	i1 = *n;
	for (k = 2; k <= i1; ++k) {
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
    if (*n > 1) {
	i = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[*n];
	indxin = iorder[*n];
    L30:
	j = i + i;
	if (j <= *n - 1) {
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
	t[*n] = out;
	iorder[*n] = indxou;
    }
    return 0;
} /* hpsolb_ */


static int matupd_(int *n, int *m, double *ws, 
		   double *wy, double *sy, double *ss, double *d, 
		   double *r, int *itail, int *iupdat, int *col, 
		   int *head, double *theta, double *rr, double *dr, 
		   double *stp, double *dtd)
{
    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	ss_dim1, ss_offset, i1, i2;

    /* Local variables */
    int j;
    int pointr;

    /* Parameter adjustments */
    --r;
    --d;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;

    /* Function Body */
    if (*iupdat <= *m) {
	*col = *iupdat;
	*itail = (*head + *iupdat - 2) % *m + 1;
    } else {
	*itail = *itail % *m + 1;
	*head = *head % *m + 1;
    }
    dcopy_(n, &d[1], &c1, &ws[*itail * ws_dim1 + 1], &c1);
    dcopy_(n, &r[1], &c1, &wy[*itail * wy_dim1 + 1], &c1);
    *theta = *rr / *dr;
    if (*iupdat > *m) {
	i1 = *col - 1;
	for (j = 1; j <= i1; ++j) {
	    dcopy_(&j, &ss[(j + 1) * ss_dim1 + 2], &c1, &ss[j * ss_dim1 + 1]
		   , &c1);
	    i2 = *col - j;
	    dcopy_(&i2, &sy[j + 1 + (j + 1) * sy_dim1], &c1, &sy[j + j * 
								 sy_dim1], &c1);
	}
    }
    pointr = *head;
    i1 = *col - 1;
    for (j = 1; j <= i1; ++j) {
	sy[*col + j * sy_dim1] = ddot_(n, &d[1], &c1, &wy[pointr * 
							  wy_dim1 + 1], &c1);
	ss[j + *col * ss_dim1] = ddot_(n, &ws[pointr * ws_dim1 + 1], &c1, &
				       d[1], &c1);
	pointr = pointr % *m + 1;
    }
    if (*stp == 1.) {
	ss[*col + *col * ss_dim1] = *dtd;
    } else {
	ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
    }
    sy[*col + *col * sy_dim1] = *dr;
    return 0;
} /* matupd_ */

static int projgr_(int *n, double *l, double *u, 
		   int *nbd, double *x, double *g, double *sbgnrm)
{
    /* System generated locals */
    int i1;
    double d1, d2;

    /* Local variables */
    int i;
    double gi;

    /* Parameter adjustments */
    --g;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    *sbgnrm = 0.;
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	gi = g[i];
	if (nbd[i] != 0) {
	    if (gi < 0.) {
		if (nbd[i] >= 2) {
		    d1 = x[i] - u[i];
		    gi = max(d1,gi);
		}
	    } else {
		if (nbd[i] <= 2) {
		    d1 = x[i] - l[i];
		    gi = min(d1,gi);
		}
	    }
	}
	d1 = *sbgnrm, d2 = abs(gi);
	*sbgnrm = max(d1,d2);
    }
    return 0;
} /* projgr_ */

static int subsm_(int *n, int *m, int *nsub, int *
		  ind, double *l, double *u, int *nbd, double *x, 
		  double *d, double *ws, double *wy, double *theta, 
		  int *col, int *head, int *iword, double *wv, 
		  double *wn, int *info)
{
    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, wn_dim1, wn_offset, i1, 
	i2;

    /* Local variables */
    static int ibd;
    int i, j, k, m2;
    double dk;
    int js, jy, col2;
    double temp1, temp2, alpha;
    int pointr;

    /* Parameter adjustments */
    --d;
    --x;
    --nbd;
    --u;
    --l;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    --wv;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    --ind;

    /* Function Body */
    if (*nsub <= 0) {
	return 0;
    }
    pointr = *head;
    i1 = *col;
    for (i = 1; i <= i1; ++i) {
	temp1 = 0.;
	temp2 = 0.;
	i2 = *nsub;
	for (j = 1; j <= i2; ++j) {
	    k = ind[j];
	    temp1 += wy[k + pointr * wy_dim1] * d[j];
	    temp2 += ws[k + pointr * ws_dim1] * d[j];
	}
	wv[i] = temp1;
	wv[*col + i] = *theta * temp2;
	pointr = pointr % *m + 1;
    }
    m2 = *m << 1;
    col2 = *col << 1;
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c11, info);
    if (*info != 0) {
	return 0;
    }
    i1 = *col;
    for (i = 1; i <= i1; ++i) {
	wv[i] = -wv[i];
    }
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c1, info);
    if (*info != 0) {
	return 0;
    }
    pointr = *head;
    i1 = *col;
    for (jy = 1; jy <= i1; ++jy) {
	js = *col + jy;
	i2 = *nsub;
	for (i = 1; i <= i2; ++i) {
	    k = ind[i];
	    d[i] = d[i] + wy[k + pointr * wy_dim1] * wv[jy] / *theta 
		+ ws[k + pointr * ws_dim1] * wv[js];
	}
	pointr = pointr % *m + 1;
    }
    i1 = *nsub;
    for (i = 1; i <= i1; ++i) {
	d[i] /= *theta;
    }
    alpha = 1.;
    temp1 = alpha;
    i1 = *nsub;
    for (i = 1; i <= i1; ++i) {
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
    i1 = *nsub;
    for (i = 1; i <= i1; ++i) {
	k = ind[i];
	x[k] += alpha * d[i];
    }
    if (alpha < 1.) {
	*iword = 1;
    } else {
	*iword = 0;
    }
    return 0;
} /* subsm_ */

static int dcstep_(double *stx, double *fx, double *dx, 
		   double *sty, double *fy, double *dy, double *stp, 
		   double *fp, double *dp, int *brackt, double *stpmin, 
		   double *stpmax)
{
    /* System generated locals */
    double d1, d2, d3;

    /* Local variables */
    double p, q, r, s, sgnd, stpc, stpf, stpq, gamma, theta;

    sgnd = *dp * (*dx / abs(*dx));
    if (*fp > *fx) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(
								  *dp);
	s = max(d1,d2);
	d1 = theta / s;
	gamma = s * sqrt(d1 * d1 - *dx / s * (*dp / s));
	if (*stp < *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dx + theta;
	q = gamma - *dx + gamma + *dp;
	r = p / q;
	stpc = *stx + r * (*stp - *stx);
	stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp 
									- *stx);
	if ((d1 = stpc - *stx, abs(d1)) < (d2 = stpq - *stx, abs(d2)))
	    {
		stpf = stpc;
	    } else {
	    stpf = stpc + (stpq - stpc) / 2.;
	}
	*brackt = 1;
    } else if (sgnd < 0.) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(
								  *dp);
	s = max(d1,d2);
	d1 = theta / s;
	gamma = s * sqrt(d1 * d1 - *dx / s * (*dp / s));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma - *dp + gamma + *dx;
	r = p / q;
	stpc = *stp + r * (*stx - *stp);
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if ((d1 = stpc - *stp, abs(d1)) > (d2 = stpq - *stp, abs(d2)))
	    {
		stpf = stpc;
	    } else {
	    stpf = stpq;
	}
	*brackt = 1;
    } else if (abs(*dp) < abs(*dx)) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
	d1 = abs(theta), d2 = abs(*dx), d1 = max(d1,d2), d2 = abs(
								  *dp);
	s = max(d1,d2);
	d3 = theta / s;
	d1 = 0., d2 = d3 * d3 - *dx / s * (*dp / s);
	gamma = s * sqrt((max(d1,d2)));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma + (*dx - *dp) + gamma;
	r = p / q;
	if (r < 0. && gamma != 0.) {
	    stpc = *stp + r * (*stx - *stp);
	} else if (*stp > *stx) {
	    stpc = *stpmax;
	} else {
	    stpc = *stpmin;
	}
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if (*brackt) {
	    if ((d1 = stpc - *stp, abs(d1)) < (d2 = stpq - *stp, abs(
								     d2))) {
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
	    if ((d1 = stpc - *stp, abs(d1)) > (d2 = stpq - *stp, abs(
								     d2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    stpf = min(*stpmax,stpf);
	    stpf = max(*stpmin,stpf);
	}
    } else {
	if (*brackt) {
	    theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
	    d1 = abs(theta), d2 = abs(*dy), d1 = max(d1,d2), d2 = 
		abs(*dp);
	    s = max(d1,d2);
	    d1 = theta / s;
	    gamma = s * sqrt(d1 * d1 - *dy / s * (*dp / s));
	    if (*stp > *sty) {
		gamma = -gamma;
	    }
	    p = gamma - *dp + theta;
	    q = gamma - *dp + gamma + *dy;
	    r = p / q;
	    stpc = *stp + r * (*sty - *stp);
	    stpf = stpc;
	} else if (*stp > *stx) {
	    stpf = *stpmax;
	} else {
	    stpf = *stpmin;
	}
    }
    if (*fp > *fx) {
	*sty = *stp;
	*fy = *fp;
	*dy = *dp;
    } else {
	if (sgnd < 0.) {
	    *sty = *stx;
	    *fy = *fx;
	    *dy = *dx;
	}
	*stx = *stp;
	*fx = *fp;
	*dx = *dp;
    }
    *stp = stpf;
    return 0;
} /* dcstep_ */

static int dcsrch_(double *f, double *g, double *stp, 
		   double *ftol, double *gtol, double *xtol, double *
		   stpmin, double *stpmax, char *task, int *isave, double *
		   dsave)
{
    /* System generated locals */
    double d1;

    /* Local variables */
    double fm, gm, fx, fy, gx, gy, fxm, fym, gxm, gym, stx, sty;
    int stage;
    double finit, ginit, width, ftest, gtest, stmin, stmax, width1;
    int brackt;

    /* Parameter adjustments */
    --dsave;
    --isave;

    /* Function Body */
    if (strncmp(task, "START", 5) == 0) {
	if (*stp < *stpmin) {
	    strcpy(task, "ERROR: STP .LT. STPMIN");
	}
	if (*stp > *stpmax) {
	    strcpy(task, "ERROR: STP .GT. STPMAX");
	}
	if (*g >= 0.) {
	    strcpy(task, "ERROR: INITIAL G .GE. ZERO");
	}
	if (*ftol < 0.) {
	    strcpy(task, "ERROR: FTOL .LT. ZERO");
	}
	if (*gtol < 0.) {
	    strcpy(task, "ERROR: GTOL .LT. ZERO");
	}
	if (*xtol < 0.) {
	    strcpy(task, "ERROR: XTOL .LT. ZERO");
	}
	if (*stpmin < 0.) {
	    strcpy(task, "ERROR: STPMIN .LT. ZERO");
	}
	if (*stpmax < *stpmin) {
	    strcpy(task, "ERROR: STPMAX .LT. STPMIN");
	}
	if (strncmp(task, "ERROR", 5) == 0) {
	    return 0;
	}
	brackt = 0;
	stage = 1;
	finit = *f;
	ginit = *g;
	gtest = *ftol * ginit;
	width = *stpmax - *stpmin;
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
    if (stage == 1 && *f <= ftest && *g >= 0.) {
	stage = 2;
    }
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
	strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
	strcpy(task, "WARNING: XTOL TEST SATISFIED");
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
	strcpy(task, "WARNING: STP = STPMAX");
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
	strcpy(task, "WARNING: STP = STPMIN");
    }
    if (*f <= ftest && abs(*g) <= *gtol * (-ginit)) {
	strcpy(task, "CONVERGENCE");
    }
    if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4) == 0) {
	goto L1000;
    }
    if (stage == 1 && *f <= fx && *f > ftest) {
	fm = *f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = *g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
		stmin, &stmax);
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
		stmax);
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
    *stp = max(*stp,*stpmin);
    *stp = min(*stp,*stpmax);
    if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin 
	<= *xtol * stmax) {
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
} /* dcsrch_ */

static int lnsrlb_(int *n, double *l, double *u, 
		   int *nbd, double *x, double *f, double *fold, 
		   double *gd, double *gdold, double *g, double *d, 
		   double *r, double *t, double *z, double *stp, 
		   double *dnorm, double *dtd, double *xstep, double *
		   stpmx, int *iter, int *ifun, int *iback, int *nfgv, 
		   int *info, char *task, int *boxed, int *cnstnd, char *
		   csave, int *isave, double *dsave)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    int i;
    double a1, a2;

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

    /* Function Body */
    if (strncmp(task, "FG_LN", 5) == 0) {
	goto L556;
    }
    *dtd = ddot_(n, &d[1], &c1, &d[1], &c1);
    *dnorm = sqrt(*dtd);
    *stpmx = 1e10;
    if (*cnstnd) {
	if (*iter == 0) {
	    *stpmx = 1.;
	} else {
	    i1 = *n;
	    for (i = 1; i <= i1; ++i) {
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
    if (*iter == 0 && ! (*boxed)) {
	d1 = 1. / *dnorm;
	*stp = min(d1,*stpmx);
    } else {
	*stp = 1.;
    }
    dcopy_(n, &x[1], &c1, &t[1], &c1);
    dcopy_(n, &g[1], &c1, &r[1], &c1);
    *fold = *f;
    *ifun = 0;
    *iback = 0;
    strcpy(csave, "START");
 L556:
    *gd = ddot_(n, &g[1], &c1, &d[1], &c1);
    if (*ifun == 0) {
	*gdold = *gd;
	if (*gd >= 0.) {
	    *info = -4;
	    return 0;
	}
    }
    dcsrch_(f, gd, stp, &c_b173, &c_b174, &c_b175, &c_b176, stpmx, csave, &
	    isave[1], &dsave[1]);
    *xstep = *stp * *dnorm;
    if (strncmp(csave, "CONV", 4) != 0 && strncmp(csave, "WARN", 4) != 0) {
	strcpy(task, "FG_LNSRCH");
	++(*ifun);
	++(*nfgv);
	*iback = *ifun - 1;
	if (*stp == 1.) {
	    dcopy_(n, &z[1], &c1, &x[1], &c1);
	} else {
	    i1 = *n;
	    for (i = 1; i <= i1; ++i) {
		x[i] = *stp * d[i] + t[i];
	    }
	}
    } else {
	strcpy(task, "NEW_X");
    }
    return 0;
} /* lnsrlb_ */

static void timer_(double *ttime)
{
    clock_t tim = clock();
    *ttime = (double) tim / CLOCKS_PER_SEC;
} /* timer_ */

static double dpmeps_(void)
{
    /* Initialized data */

    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;

    /* System generated locals */
    int i1;
    double ret_val;

    /* Local variables */
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
    i1 = negep;
    for (i = 1; i <= i1; ++i) {
	a *= betain;
	/* L40: */
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

static int cauchy_(int *n, double *x, double *l, 
		   double *u, int *nbd, double *g, int *iorder, int *
		   iwhere, double *t, double *d, double *xcp, int *m, 
		   double *wy, double *ws, double *sy, double *wt, 
		   double *theta, int *col, int *head, double *p, 
		   double *c, double *wbp, double *v, int *nint, 
		   double *sg, double *yg, double *sbgnrm, 
		   int *info, double *epsmch)
{
    /* System generated locals */
    int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset, 
	wt_dim1, wt_offset, i1, i2;
    double d1;

    /* Local variables */
    static double tu, tl;
    int i, j;
    double f1, f2, dt, tj, tj0;
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
    --d;
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
    --c;
    --p;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;

    /* Function Body */
    if (*sbgnrm <= 0.) {
	dcopy_(n, &x[1], &c1, &xcp[1], &c1);
	return 0;
    }
    bnded = 1;
    nfree = *n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = *col << 1;
    f1 = 0.;
    i1 = col2;
    for (i = 1; i <= i1; ++i) {
	p[i] = 0.;
    }
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
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
	pointr = *head;
	if (iwhere[i] != 0 && iwhere[i] != -1) {
	    d[i] = 0.;
	} else {
	    d[i] = neggi;
	    f1 -= neggi * neggi;
	    i2 = *col;
	    for (j = 1; j <= i2; ++j) {
		p[j] += wy[i + pointr * wy_dim1] * neggi;
		p[*col + j] += ws[i + pointr * ws_dim1] * neggi;
		pointr = pointr % *m + 1;
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
    if (*theta != 1.) {
	dscal_(col, theta, &p[*col + 1], &c1);
    }
    dcopy_(n, &x[1], &c1, &xcp[1], &c1);
    if (nbreak == 0 && nfree == *n + 1) {
	return 0;
    }
    i1 = col2;
    for (j = 1; j <= i1; ++j) {
	c[j] = 0.;
    }
    f2 = -(*theta) * f1;
    f2_org__ = f2;
    if (*col > 0) {
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	f2 -= ddot_(&col2, &v[1], &c1, &p[1], &c1);
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
	i1 = iter - 2;
	hpsolb_(&nleft, &t[1], &iorder[1], &i1);
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
    if (nleft == 0 && nbreak == *n) {
	dtm = dt;
	goto L999;
    }
    ++(*nint);
    d1 = dibp;
    dibp2 = d1 * d1;
    f1 = f1 + dt * f2 + dibp2 - *theta * dibp * zibp;
    f2 -= *theta * dibp2;
    if (*col > 0) {
	daxpy_(&col2, &dt, &p[1], &c1, &c[1], &c1);
	pointr = *head;
	i1 = *col;
	for (j = 1; j <= i1; ++j) {
	    wbp[j] = wy[ibp + pointr * wy_dim1];
	    wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
	    pointr = pointr % *m + 1;
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	wmc = ddot_(&col2, &c[1], &c1, &v[1], &c1);
	wmp = ddot_(&col2, &p[1], &c1, &v[1], &c1);
	wmw = ddot_(&col2, &wbp[1], &c1, &v[1], &c1);
	d1 = -dibp;
	daxpy_(&col2, &d1, &wbp[1], &c1, &p[1], &c1);
	f1 += dibp * wmc;
	f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }
    d1 = *epsmch * f2_org__;
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
    daxpy_(n, &tsum, &d[1], &c1, &xcp[1], &c1);
 L999:
    if (*col > 0) {
	daxpy_(&col2, &dtm, &p[1], &c1, &c[1], &c1);
    }
    return 0;
} /* cauchy_ */

static int mainlb_(int *n, int *m, double *x, 
		   double *l, double *u, int *nbd, double *f, double 
		   *g, double *factr, double *pgtol, double *ws, double *
		   wy, double *sy, double *ss, double *yy, double *wt, 
		   double *wn, double *snd, double *z, double *r, 
		   double *d, double *t, double *wa, double *sg, 
		   double *sgo, double *yg, double *ygo, int *index, 
		   int *iwhere, int *indx2, char *task, char *
		   csave, int *lsave, int *isave, double *dsave)
{
    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	ss_dim1, ss_offset, yy_dim1, yy_offset, wt_dim1, wt_offset, 
	wn_dim1, wn_offset, snd_dim1, snd_offset, i1;
    double d1, d2;

    /* Local variables */
    int i, k;
    double gd, dr, rr, dtd;
    int col;
    double tol;
    int wrk;
    double stp, cpu1, cpu2;
    int head;
    double fold;
    double ddum;
    int info;
    double time;
    int nfgv, ifun, iter, nint;
    double time1, time2;
    int iback;
    double gdold;
    int nfree;
    int boxed;
    double theta;
    double dnorm;
    int nskip;
    double xstep, stpmx;
    int ileave;
    double cachyt;
    double epsmch;
    int updatd;
    double sbtime;
    int prjctd;
    int iupdat;
    int cnstnd;
    double sbgnrm;
    int nenter;
    double lnscht;
    int nintol;

    static int iword, itail, nact, itfile;

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
    snd_dim1 = 2 * *m;
    snd_offset = 1 + snd_dim1;
    snd -= snd_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    yy_dim1 = *m;
    yy_offset = 1 + yy_dim1;
    yy -= yy_offset;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    --lsave;
    --isave;
    --dsave;

    /* Function Body */
    if (strncmp(task, "START", 5) == 0) {
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
	nfree = *n;
#if 0
	tol = *factr * epsmch;
#else
	tol = *factr / 100.0;
#endif
	cachyt = 0.;
	sbtime = 0.;
	lnscht = 0.;
	info = 0;
	errclb_(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k);
	if (strncmp(task, "ERROR", 5) == 0) {
	    return 0;
	}
	active_(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], &prjctd, 
		&cnstnd, &boxed);
    } else {
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
	if (strncmp(task, "FG_LN", 5) == 0) {
	    goto L666;
	}
	if (strncmp(task, "NEW_X", 5) == 0) {
	    goto L777;
	}
	if (strncmp(task, "FG_ST", 5) == 0) {
	    goto L111;
	}
	if (strncmp(task, "STOP", 4) == 0) {
	    goto L999;
	}
    }
    strcpy(task, "FG_START");
    goto L1000;
 L111:
    nfgv = 1;
    projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
    if (sbgnrm <= *pgtol) {
	strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }
 L222:
    iword = -1;
    if (! cnstnd && col > 0) {
	dcopy_(n, &x[1], &c1, &z[1], &c1);
	wrk = updatd;
	nint = 0;
	goto L333;
    }
    timer_(&cpu1);
    cauchy_(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], &t[
									      1], &d[1], &z[1], m, &wy[wy_offset], &ws[ws_offset], &sy[
																       sy_offset], &wt[wt_offset], &theta, &col, &head, &wa[1], &wa[(*m 
																								     << 1) + 1], &wa[(*m << 2) + 1], &wa[*m * 6 + 1], &nint, &sg[1], &
	    yg[1], &sbgnrm, &info, &epsmch);
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
    freev_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iwhere[1], &
	   wrk, &updatd, &cnstnd, &iter);
    nact = *n - nfree;
 L333:
    if (nfree == 0 || col == 0) {
	goto L555;
    }
    timer_(&cpu1);
    if (wrk) {
	formk_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iupdat, &
	       updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset], &
	       wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
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
    cmprlb_(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], &sy[sy_offset]
	    , &wt[wt_offset], &z[1], &r[1], &wa[1], &index[1], &theta, &
	    col, &head, &nfree, &cnstnd, &info);
    if (info != 0) {
	goto L444;
    }
    subsm_(n, m, &nfree, &index[1], &l[1], &u[1], &nbd[1], &z[1], &r[1], &
	   ws[ws_offset], &wy[wy_offset], &theta, &col, &head, &iword, &wa[1]
	   , &wn[wn_offset], &info);
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
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	d[i] = z[i] - x[i];
    }
    timer_(&cpu1);
 L666:
    lnsrlb_(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1], &
	    d[1], &r[1], &t[1], &z[1], &stp, &dnorm, &dtd, &xstep, &
	    stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed, &cnstnd, 
	    csave, &isave[22], &dsave[17]);
    if (info != 0 || iback >= 20) {
	dcopy_(n, &t[1], &c1, &x[1], &c1);
	dcopy_(n, &r[1], &c1, &g[1], &c1);
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
    } else if (strncmp(task, "FG_LN", 5) == 0) {
	goto L1000;
    } else {
	timer_(&cpu2);
	lnscht = lnscht + cpu2 - cpu1;
	++iter;
	projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
	goto L1000;
    }
 L777:
    if (sbgnrm <= *pgtol) {
	strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }
    d1 = abs(fold), d2 = abs(*f), d1 = max(d1,d2);
    ddum = max(d1,1.);
    if (fold - *f <= tol * ddum) {
	strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
	if (iback >= 10) {
	    info = -5;
	}
	goto L999;
    }
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	r[i] = g[i] - r[i];
    }
    rr = ddot_(n, &r[1], &c1, &r[1], &c1);
    if (stp == 1.) {
	dr = gd - gdold;
	ddum = -gdold;
    } else {
	dr = (gd - gdold) * stp;
	dscal_(n, &stp, &d[1], &c1);
	ddum = -gdold * stp;
    }
    if (dr <= epsmch * ddum) {
	++nskip;
	updatd = 0;
	goto L888;
    }
    updatd = 1;
    ++iupdat;
    matupd_(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], &ss[
								      ss_offset], &d[1], &r[1], &itail, &iupdat, &col, &head, &
	    theta, &rr, &dr, &stp, &dtd);
    formt_(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col, &theta, &
	   info);
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
    goto L222;
 L999:
    timer_(&time2);
    time = time2 - time1;
 L1000:
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

int setulb_(int *n, int *m, double *x, 
	    double *l, double *u, int *nbd, double *f, double *g, 
	    double *factr, double *pgtol, double *wa, int *iwa, 
	    char *task, char *csave, int *lsave, 
	    int *isave, double *dsave)
{
    /* System generated locals */
    int i1;

    /* Local variables */
    static int l1, l2, l3, ld, lr, lt, lz, lwa, lsg, lyg, lwn, lss, lws, 
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

    /* Function Body */
    if (strncmp(task, "START", 5) == 0) {
	isave[1] = *m * *n;
	i1 = *m;
	isave[2] = i1 * i1;
	i1 = *m;
	isave[3] = i1 * i1 << 2;
	isave[4] = 1;
	isave[5] = isave[4] + isave[1];
	isave[6] = isave[5] + isave[1];
	isave[7] = isave[6] + isave[2];
	isave[8] = isave[7] + isave[2];
	isave[9] = isave[8] + isave[2];
	isave[10] = isave[9] + isave[2];
	isave[11] = isave[10] + isave[3];
	isave[12] = isave[11] + isave[3];
	isave[13] = isave[12] + *n;
	isave[14] = isave[13] + *n;
	isave[15] = isave[14] + *n;
	isave[16] = isave[15] + *n;
	isave[17] = isave[16] + (*m << 3);
	isave[18] = isave[17] + *m;
	isave[19] = isave[18] + *m;
	isave[20] = isave[19] + *m;
    }
    l1 = isave[1];
    l2 = isave[2];
    l3 = isave[3];
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
    mainlb_(n, m, &x[1], &l[1], &u[1], &nbd[1], f, &g[1], factr, pgtol, &wa[
									    lws], &wa[lwy], &wa[lsy], &wa[lss], &wa[lyy], &wa[lwt], &wa[lwn], 
	    &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa], &wa[lsg],
	    &wa[lsgo], &wa[lyg], &wa[lygo], &iwa[1], &iwa[*n + 1], &iwa[(*n 
									 << 1) + 1], task, csave, &lsave[1], &isave[22], &dsave[1]
	    );
    return 0;
} /* setulb_ */
