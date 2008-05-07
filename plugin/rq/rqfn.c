/* rqfn.f -- translated by Ratfor, version 1.0, then f2c (version 20050501).
*/

#include "gretl_f2c.h"

static doublereal c_b4 = 1.;
static integer c__1 = 1;
static doublereal c_b6 = 0.;
static doublereal c_b20 = -1.;

int rqfn_(integer *n, integer *p, doublereal *a, doublereal *
	y, doublereal *rhs, doublereal *d__, doublereal *u, doublereal *beta, 
	doublereal *eps, doublereal *wn, doublereal *wp, doublereal *aa, 
	integer *nit, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, wn_dim1, wn_offset, wp_dim1, wp_offset, aa_dim1,
	     aa_offset;

    /* Local variables */
    extern /* Subroutine */ int fna_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

    /* Parameter adjustments */
    wn_dim1 = *n;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    --u;
    --d__;
    --y;
    aa_dim1 = *p;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    wp_dim1 = *p;
    wp_offset = 1 + wp_dim1;
    wp -= wp_offset;
    --rhs;
    a_dim1 = *p;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --nit;

    /* Function Body */
    fna_(n, p, &a[a_offset], &y[1], &rhs[1], &d__[1], &u[1], beta, eps, &wn[
	    wn_dim1 + 1], &wn[(wn_dim1 << 1) + 1], &wp[wp_dim1 + 1], &wn[
	    wn_dim1 * 3 + 1], &wn[(wn_dim1 << 2) + 1], &wn[wn_dim1 * 5 + 1], &
	    wn[wn_dim1 * 6 + 1], &wp[(wp_dim1 << 1) + 1], &wn[wn_dim1 * 7 + 1]
	    , &wn[(wn_dim1 << 3) + 1], &wn[wn_dim1 * 9 + 1], &wn[wn_dim1 * 10 
	    + 1], &wp[wp_dim1 * 3 + 1], &wp[(wp_dim1 << 2) + 1], &aa[
	    aa_offset], &nit[1], info);
    return 0;
} /* rqfn_ */

/* Subroutine */ int fna_(integer *n, integer *p, doublereal *a, doublereal *
	c__, doublereal *b, doublereal *d__, doublereal *u, doublereal *beta, 
	doublereal *eps, doublereal *x, doublereal *s, doublereal *y, 
	doublereal *z__, doublereal *w, doublereal *dx, doublereal *ds, 
	doublereal *dy, doublereal *dz, doublereal *dw, doublereal *dsdw, 
	doublereal *dxdz, doublereal *rhs, doublereal *ada, doublereal *aa, 
	integer *nit, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ada_dim1, ada_offset, aa_dim1, aa_offset, i__1, 
	    i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal g;
    static integer i__, j;
    static doublereal cx;
    static integer pp;
    static doublereal by, mu, uw, uz, rdg, mua;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal acomp;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), stepy_(integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal deltad, deltap;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dtrtrs_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);

    /* Parameter adjustments */
    --dxdz;
    --dsdw;
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
    aa_dim1 = *p;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
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
    stepy_(n, p, &a[a_offset], &d__[1], &y[1], &aa[aa_offset], info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *p;
	for (j = 1; j <= i__2; ++j) {
	    ada[i__ + j * ada_dim1] = 0.;
/* L23004: */
	}
/* L23005: */
	ada[i__ + i__ * ada_dim1] = 1.;
/* L23002: */
    }
/* L23003: */
    dtrtrs_("U", "T", "N", p, p, &aa[aa_offset], p, &ada[ada_offset], p, info,
	     (ftnlen)1, (ftnlen)1, (ftnlen)1);
    dcopy_(&pp, &ada[ada_offset], &c__1, &aa[aa_offset], &c__1);
    dcopy_(n, &c__[1], &c__1, &s[1], &c__1);
    dgemv_("T", p, n, &c_b20, &a[a_offset], p, &y[1], &c__1, &c_b4, &s[1], &
	    c__1, (ftnlen)1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = 1.;
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
/* L23006: */
    }
/* L23007: */
    cx = ddot_(n, &c__[1], &c__1, &x[1], &c__1);
    by = ddot_(p, &b[1], &c__1, &y[1], &c__1);
    uw = dasum_(n, &w[1], &c__1);
    uz = dasum_(n, &z__[1], &c__1);
    rdg = cx - by + uw;
L23010:
    if (rdg > *eps) {
	++nit[1];
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__[i__] = 1. / (z__[i__] / x[i__] + w[i__] / s[i__]);
	    ds[i__] = z__[i__] - w[i__];
	    dx[i__] = d__[i__] * ds[i__];
/* L23012: */
	}
/* L23013: */
	dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dx[1], &c__1, &c_b6, &dy[1]
		, &c__1, (ftnlen)1);
	dcopy_(p, &dy[1], &c__1, &rhs[1], &c__1);
	stepy_(n, p, &a[a_offset], &d__[1], &dy[1], &ada[ada_offset], info);
	if (*info != 0) {
	    return 0;
	}
	dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &c__1, &c_b20, &ds[
		1], &c__1, (ftnlen)1);
	deltap = 1e20;
	deltad = 1e20;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dx[i__] = d__[i__] * ds[i__];
	    ds[i__] = -dx[i__];
	    dz[i__] = -z__[i__] * (dx[i__] / x[i__] + 1.);
	    dw[i__] = w[i__] * (dx[i__] / s[i__] - 1.);
	    dxdz[i__] = dx[i__] * dz[i__];
	    dsdw[i__] = ds[i__] * dw[i__];
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
/* L23016: */
	}
/* L23017: */
/* Computing MIN */
	d__1 = *beta * deltap;
	deltap = min(d__1,1.);
/* Computing MIN */
	d__1 = *beta * deltad;
	deltad = min(d__1,1.);
	if (deltap * deltad < 1.) {
	    ++nit[2];
	    acomp = ddot_(n, &x[1], &c__1, &z__[1], &c__1) + ddot_(n, &s[1], &
		    c__1, &w[1], &c__1);
	    g = acomp + deltap * ddot_(n, &dx[1], &c__1, &z__[1], &c__1) + 
		    deltad * ddot_(n, &dz[1], &c__1, &x[1], &c__1) + deltap * 
		    deltad * ddot_(n, &dz[1], &c__1, &dx[1], &c__1) + deltap *
		     ddot_(n, &ds[1], &c__1, &w[1], &c__1) + deltad * ddot_(n,
		     &dw[1], &c__1, &s[1], &c__1) + deltap * deltad * ddot_(n,
		     &ds[1], &c__1, &dw[1], &c__1);
	    mu = acomp / (doublereal) (*n << 1);
	    mua = g / (doublereal) (*n << 1);
/* Computing 3rd power */
	    d__1 = mua / mu;
	    mu *= d__1 * (d__1 * d__1);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dz[i__] = d__[i__] * (mu * (1 / s[i__] - 1 / x[i__]) + dx[i__]
			 * dz[i__] / x[i__] - ds[i__] * dw[i__] / s[i__]);
/* L23028: */
	    }
/* L23029: */
	    dswap_(p, &rhs[1], &c__1, &dy[1], &c__1);
	    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dz[1], &c__1, &c_b4, &
		    dy[1], &c__1, (ftnlen)1);
	    dpotrs_("U", p, &c__1, &ada[ada_offset], p, &dy[1], p, info, (
		    ftnlen)1);
	    daxpy_(p, &c_b20, &dy[1], &c__1, &rhs[1], &c__1);
	    dgemv_("T", p, n, &c_b4, &a[a_offset], p, &rhs[1], &c__1, &c_b6, &
		    dw[1], &c__1, (ftnlen)1);
	    deltap = 1e20;
	    deltad = 1e20;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dx[i__] = dx[i__] - dz[i__] - d__[i__] * dw[i__];
		ds[i__] = -dx[i__];
		dz[i__] = mu / x[i__] - z__[i__] * dx[i__] / x[i__] - z__[i__]
			 - dxdz[i__] / x[i__];
		dw[i__] = mu / s[i__] - w[i__] * ds[i__] / s[i__] - w[i__] - 
			dsdw[i__] / s[i__];
		if (dx[i__] < 0.) {
/* Computing MIN */
		    d__1 = deltap, d__2 = -x[i__] / dx[i__];
		    deltap = min(d__1,d__2);
		} else {
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
/* L23030: */
	    }
/* L23031: */
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
	cx = ddot_(n, &c__[1], &c__1, &x[1], &c__1);
	by = ddot_(p, &b[1], &c__1, &y[1], &c__1);
	uw = dasum_(n, &w[1], &c__1);
	uz = dasum_(n, &z__[1], &c__1);
	rdg = cx - by + uw;
	goto L23010;
    }
/* L23011: */
    daxpy_(n, &c_b20, &w[1], &c__1, &z__[1], &c__1);
    dswap_(n, &z__[1], &c__1, &x[1], &c__1);
    return 0;
} /* fna_ */

