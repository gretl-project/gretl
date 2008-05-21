/* rqfnb.f -- translated by f2c (version 20050501) */

#include "gretl_f2c.h"

/* Table of constant values */

static doublereal c_b4 = 1.;
static integer one = 1;
static doublereal zero = 0.;
static doublereal c_b13 = -1.;

/* lapack/blas functions called below */

extern int dsyr_(char *, integer *, doublereal *, doublereal *, integer *, 
		 doublereal *, integer *, ftnlen);

extern int dposv_(char *, integer *, integer *, doublereal *, integer *, 
		  doublereal *, integer *, integer *, ftnlen);

extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
		  integer *, doublereal *, integer *, doublereal *, doublereal *, 
		  integer *, ftnlen);

extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);

extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);

extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, 
		  integer *);

extern int dpotrs_(char *, integer *, integer *, doublereal *, integer *, doublereal *, 
		   integer *, integer *, ftnlen);

/* principal private function */

static int lpfnb_(integer *, integer *, doublereal *, 
		  doublereal *, doublereal *, doublereal *, doublereal *, 
		  doublereal *, doublereal *, doublereal *, doublereal *, 
		  doublereal *, doublereal *, doublereal *, doublereal *, 
		  doublereal *, doublereal *, doublereal *, doublereal *, 
		  doublereal *, doublereal *, doublereal *, integer *, integer *);

/* Output from Public domain Ratfor, version 1.0 */
int rqfnb_(integer *n, integer *p, doublereal *a, doublereal *y, 
	   doublereal *rhs, doublereal *d__, doublereal *u, doublereal *beta,
	   doublereal *eps, doublereal *wn, doublereal *wp, integer *nit, 
	   integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, wn_dim1, wn_offset, wp_dim1, wp_offset;

    /* Parameter adjustments */
    wn_dim1 = *n;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    --u;
    --d__;
    --y;
    wp_dim1 = *p;
    wp_offset = 1 + wp_dim1;
    wp -= wp_offset;
    --rhs;
    a_dim1 = *p;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --nit;

    /* Function Body */
    lpfnb_(n, p, &a[a_offset], &y[1], &rhs[1], &d__[1], &u[1], beta, eps, &wn[
	    wn_dim1 + 1], &wn[(wn_dim1 << 1) + 1], &wp[wp_dim1 + 1], &wn[
	    wn_dim1 * 3 + 1], &wn[(wn_dim1 << 2) + 1], &wn[wn_dim1 * 5 + 1], &
	    wn[wn_dim1 * 6 + 1], &wp[(wp_dim1 << 1) + 1], &wn[wn_dim1 * 7 + 1]
	    , &wn[(wn_dim1 << 3) + 1], &wn[wn_dim1 * 9 + 1], &wp[wp_dim1 * 3 
	    + 1], &wp[(wp_dim1 << 2) + 1], &nit[1], info);
    return 0;
} /* rqfnb_ */

static int stepy_(integer *n, integer *p, doublereal *a, doublereal 
		  *d, doublereal *b, doublereal *ada, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ada_dim1, ada_offset;

    /* Local variables */
    static integer i, j, k, pp;

    /* Parameter adjustments */
    --d;
    ada_dim1 = *p;
    ada_offset = 1 + ada_dim1;
    ada -= ada_offset;
    --b;
    a_dim1 = *p;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    pp = *p * *p;
    for (j = 1; j <= *p; ++j) {
	for (k = 1; k <= *p; ++k) {
	    ada[j + k * ada_dim1] = 0.;
	}
    }
    for (i = 1; i <= *n; ++i) {
	dsyr_("U", p, &d[i], &a[i * a_dim1 + 1], &one, &ada[ada_offset],
	       p, (ftnlen)1);
    }

    dposv_("U", p, &one, &ada[ada_offset], p, &b[1], p, info, (ftnlen) 1);

    return 0;
} /* stepy_ */


static int lpfnb_(integer *n, integer *p, doublereal *a, doublereal *c__, 
		  doublereal *b, doublereal *d__, doublereal *u, doublereal *beta,
		  doublereal *eps, doublereal *x, doublereal *s, doublereal *y, 
		  doublereal *z__, doublereal *w, doublereal *dx, doublereal *ds, 
		  doublereal *dy, doublereal *dz, doublereal *dw, doublereal *dr, 
		  doublereal *rhs, doublereal *ada, integer *nit, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ada_dim1, ada_offset, i1;
    doublereal d1, d2;

    /* Local variables */
    static doublereal g;
    static integer i, pp;
    static doublereal mu, gap;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal dsdw, dxdz;
    static doublereal deltad, deltap;

    /* Parameter adjustments */
    --dr;
    --dw;
    --dz;
    --ds;
    --dx;
    --w;
    --z__;
    --s;
    --x;
    --u;
    --d__;
    --c__;
    ada_dim1 = *p;
    ada_offset = 1 + ada_dim1;
    ada -= ada_offset;
    --rhs;
    --dy;
    --y;
    --b;
    a_dim1 = *p;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --nit;

    /* Function Body */
    nit[1] = 0;
    nit[2] = 0;
    nit[3] = *n;
    pp = *p * *p;
    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &c__[1], &one, &zero, &y[1], &
	    one, (ftnlen)1);
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	d__[i] = 1.;
    }
    stepy_(n, p, &a[a_offset], &d__[1], &y[1], &ada[ada_offset], info);
    if (*info != 0) {
	return 0;
    }
    dcopy_(n, &c__[1], &one, &s[1], &one);
    dgemv_("T", p, n, &c_b13, &a[a_offset], p, &y[1], &one, &c_b4, &s[1], &
	    one, (ftnlen)1);
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	if ((d1 = s[i], abs(d1)) < *eps) {
	    /* Computing MAX */
	    d1 = s[i];
	    z__[i] = max(d1,0.) + *eps;
	    /* Computing MAX */
	    d1 = -s[i];
	    w[i] = max(d1,0.) + *eps;
	} else {
	    /* Computing MAX */
	    d1 = s[i];
	    z__[i] = max(d1,0.);
	    /* Computing MAX */
	    d1 = -s[i];
	    w[i] = max(d1,0.);
	}
	s[i] = u[i] - x[i];
    }
    gap = ddot_(n, &z__[1], &one, &x[1], &one) + ddot_(n, &w[1], &one, &s[
	    1], &one);
L23008:
    if (gap > *eps && nit[1] < 50) {
	++nit[1];
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	    d__[i] = 1. / (z__[i] / x[i] + w[i] / s[i]);
	    ds[i] = z__[i] - w[i];
	    dz[i] = d__[i] * ds[i];
	}
	dcopy_(p, &b[1], &one, &dy[1], &one);
	dgemv_("N", p, n, &c_b13, &a[a_offset], p, &x[1], &one, &c_b4, &dy[1]
		, &one, (ftnlen)1);
	dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dz[1], &one, &c_b4, &dy[1]
		, &one, (ftnlen)1);
	dcopy_(p, &dy[1], &one, &rhs[1], &one);
	stepy_(n, p, &a[a_offset], &d__[1], &dy[1], &ada[ada_offset], info);
	if (*info != 0) {
	    return 0;
	}
	dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &one, &c_b13, &ds[
		1], &one, (ftnlen)1);
	deltap = 1e20;
	deltad = 1e20;
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	    dx[i] = d__[i] * ds[i];
	    ds[i] = -dx[i];
	    dz[i] = -z__[i] * (dx[i] / x[i] + 1.);
	    dw[i] = -w[i] * (ds[i] / s[i] + 1.);
	    if (dx[i] < 0.) {
		d1 = deltap, d2 = -x[i] / dx[i];
		deltap = min(d1,d2);
	    }
	    if (ds[i] < 0.) {
		d1 = deltap, d2 = -s[i] / ds[i];
		deltap = min(d1,d2);
	    }
	    if (dz[i] < 0.) {
		d1 = deltad, d2 = -z__[i] / dz[i];
		deltad = min(d1,d2);
	    }
	    if (dw[i] < 0.) {
		d1 = deltad, d2 = -w[i] / dw[i];
		deltad = min(d1,d2);
	    }
	}
	/* Computing MIN */
	d1 = *beta * deltap;
	deltap = min(d1,1.);
	/* Computing MIN */
	d1 = *beta * deltad;
	deltad = min(d1,1.);
	if (min(deltap,deltad) < 1.) {
	    ++nit[2];
	    mu = ddot_(n, &x[1], &one, &z__[1], &one) + ddot_(n, &s[1], &
		    one, &w[1], &one);
	    g = mu + deltap * ddot_(n, &dx[1], &one, &z__[1], &one) + 
		    deltad * ddot_(n, &dz[1], &one, &x[1], &one) + deltap * 
		    deltad * ddot_(n, &dz[1], &one, &dx[1], &one) + deltap *
		     ddot_(n, &ds[1], &one, &w[1], &one) + deltad * ddot_(n,
		     &dw[1], &one, &s[1], &one) + deltap * deltad * ddot_(n,
		     &ds[1], &one, &dw[1], &one);
	    /* Computing 3rd power */
	    d1 = g / mu;
	    mu = mu * (d1 * (d1 * d1)) / (doublereal) (*n << 1);
	    i1 = *n;
	    for (i = 1; i <= i1; ++i) {
		dr[i] = d__[i] * (mu * (1 / s[i] - 1 / x[i]) + dx[i]
			 * dz[i] / x[i] - ds[i] * dw[i] / s[i]);
	    }
	    dswap_(p, &rhs[1], &one, &dy[1], &one);
	    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dr[1], &one, &c_b4,
		    &dy[1], &one, (ftnlen)1);
	    dpotrs_("U", p, &one, &ada[ada_offset], p, &dy[1], p, info, 
		    (ftnlen)1);
	    dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &one, &zero, 
		   &u[1], &one, (ftnlen)1);
	    deltap = 1e20;
	    deltad = 1e20;
	    i1 = *n;
	    for (i = 1; i <= i1; ++i) {
		dxdz = dx[i] * dz[i];
		dsdw = ds[i] * dw[i];
		dx[i] = d__[i] * (u[i] - z__[i] + w[i]) - dr[i];
		ds[i] = -dx[i];
		dz[i] = -z__[i] + (mu - z__[i] * dx[i] - dxdz) / x[
			i];
		dw[i] = -w[i] + (mu - w[i] * ds[i] - dsdw) / s[i];
		if (dx[i] < 0.) {
		    /* Computing MIN */
		    d1 = deltap, d2 = -x[i] / dx[i];
		    deltap = min(d1, d2);
		}
		if (ds[i] < 0.) {
		    /* Computing MIN */
		    d1 = deltap, d2 = -s[i] / ds[i];
		    deltap = min(d1, d2);
		}
		if (dz[i] < 0.) {
		    /* Computing MIN */
		    d1 = deltad, d2 = -z__[i] / dz[i];
		    deltad = min(d1, d2);
		}
		if (dw[i] < 0.) {
		    /* Computing MIN */
		    d1 = deltad, d2 = -w[i] / dw[i];
		    deltad = min(d1, d2);
		}
	    }
	    /* Computing MIN */
	    d1 = *beta * deltap;
	    deltap = min(d1,1.);
	    /* Computing MIN */
	    d1 = *beta * deltad;
	    deltad = min(d1,1.);
	}
	daxpy_(n, &deltap, &dx[1], &one, &x[1], &one);
	daxpy_(n, &deltap, &ds[1], &one, &s[1], &one);
	daxpy_(p, &deltad, &dy[1], &one, &y[1], &one);
	daxpy_(n, &deltad, &dz[1], &one, &z__[1], &one);
	daxpy_(n, &deltad, &dw[1], &one, &w[1], &one);
	gap = ddot_(n, &z__[1], &one, &x[1], &one) + ddot_(n, &w[1], &one, 
		&s[1], &one);
	goto L23008;
    }
    daxpy_(n, &c_b13, &w[1], &one, &z__[1], &one);
    dswap_(n, &z__[1], &one, &x[1], &one);
    return 0;
} /* lpfnb_ */


