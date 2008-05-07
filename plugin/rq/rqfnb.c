/* rqfnb.f -- translated by f2c (version 20050501).
*/

#include "gretl_f2c.h"

/* Table of constant values */

static doublereal c_b4 = 1.;
static integer c__1 = 1;
static doublereal c_b6 = 0.;
static doublereal c_b13 = -1.;

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int rqfnb_(integer *n, integer *p, doublereal *a, doublereal 
	*y, doublereal *rhs, doublereal *d__, doublereal *u, doublereal *beta,
	 doublereal *eps, doublereal *wn, doublereal *wp, integer *nit, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, wn_dim1, wn_offset, wp_dim1, wp_offset;

    /* Local variables */
    extern /* Subroutine */ int lpfnb_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);

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

/* Subroutine */ int lpfnb_(integer *n, integer *p, doublereal *a, doublereal 
	*c__, doublereal *b, doublereal *d__, doublereal *u, doublereal *beta,
	 doublereal *eps, doublereal *x, doublereal *s, doublereal *y, 
	doublereal *z__, doublereal *w, doublereal *dx, doublereal *ds, 
	doublereal *dy, doublereal *dz, doublereal *dw, doublereal *dr, 
	doublereal *rhs, doublereal *ada, integer *nit, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ada_dim1, ada_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal g;
    static integer i__, pp;
    static doublereal mu, gap;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal dsdw, dxdz;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), stepy_(integer *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *);
    static doublereal deltad, deltap;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);

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
    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &c__[1], &c__1, &c_b6, &y[1], &
	    c__1, (ftnlen)1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = 1.;
/* L23000: */
    }
/* L23001: */
    stepy_(n, p, &a[a_offset], &d__[1], &y[1], &ada[ada_offset], info);
    if (*info != 0) {
	return 0;
    }
    dcopy_(n, &c__[1], &c__1, &s[1], &c__1);
    dgemv_("T", p, n, &c_b13, &a[a_offset], p, &y[1], &c__1, &c_b4, &s[1], &
	    c__1, (ftnlen)1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = s[i__], abs(d__1)) < *eps) {
/* Computing MAX */
	    d__1 = s[i__];
	    z__[i__] = max(d__1,0.) + *eps;
/* Computing MAX */
	    d__1 = -s[i__];
	    w[i__] = max(d__1,0.) + *eps;
	} else {
/* Computing MAX */
	    d__1 = s[i__];
	    z__[i__] = max(d__1,0.);
/* Computing MAX */
	    d__1 = -s[i__];
	    w[i__] = max(d__1,0.);
	}
	s[i__] = u[i__] - x[i__];
/* L23004: */
    }
/* L23005: */
    gap = ddot_(n, &z__[1], &c__1, &x[1], &c__1) + ddot_(n, &w[1], &c__1, &s[
	    1], &c__1);
L23008:
    if (gap > *eps && nit[1] < 50) {
	++nit[1];
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__[i__] = 1. / (z__[i__] / x[i__] + w[i__] / s[i__]);
	    ds[i__] = z__[i__] - w[i__];
	    dz[i__] = d__[i__] * ds[i__];
/* L23010: */
	}
/* L23011: */
	dcopy_(p, &b[1], &c__1, &dy[1], &c__1);
	dgemv_("N", p, n, &c_b13, &a[a_offset], p, &x[1], &c__1, &c_b4, &dy[1]
		, &c__1, (ftnlen)1);
	dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dz[1], &c__1, &c_b4, &dy[1]
		, &c__1, (ftnlen)1);
	dcopy_(p, &dy[1], &c__1, &rhs[1], &c__1);
	stepy_(n, p, &a[a_offset], &d__[1], &dy[1], &ada[ada_offset], info);
	if (*info != 0) {
	    return 0;
	}
	dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &c__1, &c_b13, &ds[
		1], &c__1, (ftnlen)1);
	deltap = 1e20;
	deltad = 1e20;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dx[i__] = d__[i__] * ds[i__];
	    ds[i__] = -dx[i__];
	    dz[i__] = -z__[i__] * (dx[i__] / x[i__] + 1.);
	    dw[i__] = -w[i__] * (ds[i__] / s[i__] + 1.);
	    if (dx[i__] < 0.) {
/* Computing MIN */
		d__1 = deltap, d__2 = -x[i__] / dx[i__];
		deltap = min(d__1,d__2);
	    }
	    if (ds[i__] < 0.) {
/* Computing MIN */
		d__1 = deltap, d__2 = -s[i__] / ds[i__];
		deltap = min(d__1,d__2);
	    }
	    if (dz[i__] < 0.) {
/* Computing MIN */
		d__1 = deltad, d__2 = -z__[i__] / dz[i__];
		deltad = min(d__1,d__2);
	    }
	    if (dw[i__] < 0.) {
/* Computing MIN */
		d__1 = deltad, d__2 = -w[i__] / dw[i__];
		deltad = min(d__1,d__2);
	    }
/* L23014: */
	}
/* L23015: */
/* Computing MIN */
	d__1 = *beta * deltap;
	deltap = min(d__1,1.);
/* Computing MIN */
	d__1 = *beta * deltad;
	deltad = min(d__1,1.);
	if (min(deltap,deltad) < 1.) {
	    ++nit[2];
	    mu = ddot_(n, &x[1], &c__1, &z__[1], &c__1) + ddot_(n, &s[1], &
		    c__1, &w[1], &c__1);
	    g = mu + deltap * ddot_(n, &dx[1], &c__1, &z__[1], &c__1) + 
		    deltad * ddot_(n, &dz[1], &c__1, &x[1], &c__1) + deltap * 
		    deltad * ddot_(n, &dz[1], &c__1, &dx[1], &c__1) + deltap *
		     ddot_(n, &ds[1], &c__1, &w[1], &c__1) + deltad * ddot_(n,
		     &dw[1], &c__1, &s[1], &c__1) + deltap * deltad * ddot_(n,
		     &ds[1], &c__1, &dw[1], &c__1);
/* Computing 3rd power */
	    d__1 = g / mu;
	    mu = mu * (d__1 * (d__1 * d__1)) / (doublereal) (*n << 1);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dr[i__] = d__[i__] * (mu * (1 / s[i__] - 1 / x[i__]) + dx[i__]
			 * dz[i__] / x[i__] - ds[i__] * dw[i__] / s[i__]);
/* L23026: */
	    }
/* L23027: */
	    dswap_(p, &rhs[1], &c__1, &dy[1], &c__1);
	    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dr[1], &c__1, &c_b4, &
		    dy[1], &c__1, (ftnlen)1);
	    dpotrs_("U", p, &c__1, &ada[ada_offset], p, &dy[1], p, info, (
		    ftnlen)1);
	    dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &c__1, &c_b6, &
		    u[1], &c__1, (ftnlen)1);
	    deltap = 1e20;
	    deltad = 1e20;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dxdz = dx[i__] * dz[i__];
		dsdw = ds[i__] * dw[i__];
		dx[i__] = d__[i__] * (u[i__] - z__[i__] + w[i__]) - dr[i__];
		ds[i__] = -dx[i__];
		dz[i__] = -z__[i__] + (mu - z__[i__] * dx[i__] - dxdz) / x[
			i__];
		dw[i__] = -w[i__] + (mu - w[i__] * ds[i__] - dsdw) / s[i__];
		if (dx[i__] < 0.) {
/* Computing MIN */
		    d__1 = deltap, d__2 = -x[i__] / dx[i__];
		    deltap = min(d__1,d__2);
		}
		if (ds[i__] < 0.) {
/* Computing MIN */
		    d__1 = deltap, d__2 = -s[i__] / ds[i__];
		    deltap = min(d__1,d__2);
		}
		if (dz[i__] < 0.) {
/* Computing MIN */
		    d__1 = deltad, d__2 = -z__[i__] / dz[i__];
		    deltad = min(d__1,d__2);
		}
		if (dw[i__] < 0.) {
/* Computing MIN */
		    d__1 = deltad, d__2 = -w[i__] / dw[i__];
		    deltad = min(d__1,d__2);
		}
/* L23028: */
	    }
/* L23029: */
/* Computing MIN */
	    d__1 = *beta * deltap;
	    deltap = min(d__1,1.);
/* Computing MIN */
	    d__1 = *beta * deltad;
	    deltad = min(d__1,1.);
	}
	daxpy_(n, &deltap, &dx[1], &c__1, &x[1], &c__1);
	daxpy_(n, &deltap, &ds[1], &c__1, &s[1], &c__1);
	daxpy_(p, &deltad, &dy[1], &c__1, &y[1], &c__1);
	daxpy_(n, &deltad, &dz[1], &c__1, &z__[1], &c__1);
	daxpy_(n, &deltad, &dw[1], &c__1, &w[1], &c__1);
	gap = ddot_(n, &z__[1], &c__1, &x[1], &c__1) + ddot_(n, &w[1], &c__1, 
		&s[1], &c__1);
	goto L23008;
    }
/* L23009: */
    daxpy_(n, &c_b13, &w[1], &c__1, &z__[1], &c__1);
    dswap_(n, &z__[1], &c__1, &x[1], &c__1);
    return 0;
} /* lpfnb_ */

/* Subroutine */ int stepy_(integer *n, integer *p, doublereal *a, doublereal 
	*d__, doublereal *b, doublereal *ada, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ada_dim1, ada_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, pp;
    extern /* Subroutine */ int dsyr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dposv_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *, ftnlen);

    /* Parameter adjustments */
    --d__;
    ada_dim1 = *p;
    ada_offset = 1 + ada_dim1;
    ada -= ada_offset;
    --b;
    a_dim1 = *p;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    pp = *p * *p;
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *p;
	for (k = 1; k <= i__2; ++k) {
	    ada[j + k * ada_dim1] = 0.;
/* L23040: */
	}
/* L23041: */
/* L23038: */
    }
/* L23039: */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dsyr_("U", p, &d__[i__], &a[i__ * a_dim1 + 1], &c__1, &ada[ada_offset]
		, p, (ftnlen)1);
/* L23042: */
    }
/* L23043: */
    dposv_("U", p, &c__1, &ada[ada_offset], p, &b[1], p, info, (ftnlen)1);
    return 0;
} /* stepy_ */

