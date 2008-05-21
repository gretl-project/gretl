/* rqfn.f -- translated by f2c (version 20061008) */

#include "gretl_f2c.h"

/* Table of constant values */

static doublereal c_b4 = 1.;
static integer c__1 = 1;
static doublereal c_b6 = 0.;
static doublereal mone = -1.;

/* blas calls */

extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
			integer *);
extern int dgemv_(char *, integer *, integer *, 
		  doublereal *, doublereal *, integer *, doublereal *, integer *, 
		  doublereal *, doublereal *, integer *, ftnlen);
extern doublereal dasum_(integer *, doublereal *, integer *);
extern int dcopy_(integer *, doublereal *, integer *, 
		  doublereal *, integer *);
extern int  dswap_(integer *, doublereal *, integer *, doublereal *, 
		   integer *);
extern int daxpy_(integer *, doublereal *, doublereal *, integer *, 
		  doublereal *, integer *);
extern int stepy_(integer *, integer *, doublereal *, doublereal *, 
		  doublereal *, doublereal *, integer *);
extern int dpotrs_(char *, integer *, integer *, doublereal *, integer *, 
		   doublereal *, integer *, integer *, ftnlen);
extern int dtrtrs_(char *, char *, char *, integer *, integer *, 
		   doublereal *, integer *, doublereal *, integer *, 
		   integer *, ftnlen, ftnlen, ftnlen);

/* Output from Public domain Ratfor, version 1.0 */

static int fna_(integer *n, integer *p, doublereal *a, doublereal *c__, 
		doublereal *b, doublereal *d, doublereal *u, doublereal *beta, 
		doublereal *eps, doublereal *x, doublereal *s, doublereal *y, 
		doublereal *z__, doublereal *w, doublereal *dx, doublereal *ds, 
		doublereal *dy, doublereal *dz, doublereal *dw, doublereal *dsdw, 
		doublereal *dxdz, doublereal *rhs, doublereal *ada, doublereal *aa, 
		integer *nit, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ada_dim1, ada_offset, aa_dim1, aa_offset, i1, 
	    i2;
    doublereal d1, d2;

    static doublereal g;
    static integer i, j;
    static doublereal cx;
    static integer pp;
    static doublereal by, mu, uw, uz, rdg, mua;
    static doublereal acomp;
    static doublereal deltad, deltap;

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
    --d;
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
    stepy_(n, p, &a[a_offset], &d[1], &y[1], &aa[aa_offset], info);
    if (*info != 0) {
	return 0;
    }
    i1 = *p;
    for (i = 1; i <= i1; ++i) {
	i2 = *p;
	for (j = 1; j <= i2; ++j) {
	    ada[i + j * ada_dim1] = 0.;
	}
	ada[i + i * ada_dim1] = 1.;
    }
    dtrtrs_("U", "T", "N", p, p, &aa[aa_offset], p, &ada[ada_offset], p, info,
	     (ftnlen)1, (ftnlen)1, (ftnlen)1);
    dcopy_(&pp, &ada[ada_offset], &c__1, &aa[aa_offset], &c__1);
    dcopy_(n, &c__[1], &c__1, &s[1], &c__1);
    dgemv_("T", p, n, &mone, &a[a_offset], p, &y[1], &c__1, &c_b4, &s[1], &
	    c__1, (ftnlen)1);
    i1 = *n;
    for (i = 1; i <= i1; ++i) {
	d[i] = 1.;
	if ((d1 = s[i], abs(d1)) < *eps) {
	    d1 = s[i];
	    z__[i] = max(d1,0.) + *eps;
	    d1 = -s[i];
	    w[i] = max(d1,0.) + *eps;
	} else {
	    d1 = s[i];
	    z__[i] = max(d1,0.);
	    d1 = -s[i];
	    w[i] = max(d1,0.);
	}
	s[i] = u[i] - x[i];
    }
    cx = ddot_(n, &c__[1], &c__1, &x[1], &c__1);
    by = ddot_(p, &b[1], &c__1, &y[1], &c__1);
    uw = dasum_(n, &w[1], &c__1);
    uz = dasum_(n, &z__[1], &c__1);
    rdg = cx - by + uw;
L23010:
    if (rdg > *eps) {
	++nit[1];
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	    d[i] = 1. / (z__[i] / x[i] + w[i] / s[i]);
	    ds[i] = z__[i] - w[i];
	    dx[i] = d[i] * ds[i];
	}
	dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dx[1], &c__1, &c_b6, &dy[1]
		, &c__1, (ftnlen)1);
	dcopy_(p, &dy[1], &c__1, &rhs[1], &c__1);
	stepy_(n, p, &a[a_offset], &d[1], &dy[1], &ada[ada_offset], info);
	if (*info != 0) {
	    return 0;
	}
	dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &c__1, &mone, &ds[
		1], &c__1, (ftnlen)1);
	deltap = 1e20;
	deltad = 1e20;
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	    dx[i] = d[i] * ds[i];
	    ds[i] = -dx[i];
	    dz[i] = -z__[i] * (dx[i] / x[i] + 1.);
	    dw[i] = w[i] * (dx[i] / s[i] - 1.);
	    dxdz[i] = dx[i] * dz[i];
	    dsdw[i] = ds[i] * dw[i];
	    if (dx[i] < 0.) {
		d1 = deltap, d2 = -x[i] / dx[i];
		deltap = min(d1, d2);
	    }
	    if (ds[i] < 0.) {
		d1 = deltap, d2 = -s[i] / ds[i];
		deltap = min(d1, d2);
	    }
	    if (dz[i] < 0.) {
		d1 = deltad, d2 = -z__[i] / dz[i];
		deltad = min(d1, d2);
	    }
	    if (dw[i] < 0.) {
		d1 = deltad, d2 = -w[i] / dw[i];
		deltad = min(d1, d2);
	    }
	}
	d1 = *beta * deltap;
	deltap = min(d1,1.);
	d1 = *beta * deltad;
	deltad = min(d1,1.);
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
	    d1 = mua / mu;
	    mu *= d1 * (d1 * d1);
	    i1 = *n;
	    for (i = 1; i <= i1; ++i) {
		dz[i] = d[i] * (mu * (1 / s[i] - 1 / x[i]) + dx[i]
			 * dz[i] / x[i] - ds[i] * dw[i] / s[i]);
	    }
	    dswap_(p, &rhs[1], &c__1, &dy[1], &c__1);
	    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dz[1], &c__1, &c_b4, &
		    dy[1], &c__1, (ftnlen)1);
	    dpotrs_("U", p, &c__1, &ada[ada_offset], p, &dy[1], p, info, (
		    ftnlen)1);
	    daxpy_(p, &mone, &dy[1], &c__1, &rhs[1], &c__1);
	    dgemv_("T", p, n, &c_b4, &a[a_offset], p, &rhs[1], &c__1, &c_b6, &
		    dw[1], &c__1, (ftnlen)1);
	    deltap = 1e20;
	    deltad = 1e20;
	    i1 = *n;
	    for (i = 1; i <= i1; ++i) {
		dx[i] = dx[i] - dz[i] - d[i] * dw[i];
		ds[i] = -dx[i];
		dz[i] = mu / x[i] - z__[i] * dx[i] / x[i] - z__[i]
			 - dxdz[i] / x[i];
		dw[i] = mu / s[i] - w[i] * ds[i] / s[i] - w[i] - 
			dsdw[i] / s[i];
		if (dx[i] < 0.) {
		    d1 = deltap, d2 = -x[i] / dx[i];
		    deltap = min(d1,d2);
		} else {
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
	    d1 = *beta * deltap;
	    deltap = min(d1,1.);
	    d1 = *beta * deltad;
	    deltad = min(d1,1.);
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
    daxpy_(n, &mone, &w[1], &c__1, &z__[1], &c__1);
    dswap_(n, &z__[1], &c__1, &x[1], &c__1);
    return 0;
} /* fna_ */

int rqfn_(integer *n, integer *p, doublereal *a, doublereal *
	y, doublereal *rhs, doublereal *d, doublereal *u, doublereal *beta, 
	doublereal *eps, doublereal *wn, doublereal *wp, doublereal *aa, 
	integer *nit, integer *info)
{
    integer a_dim1, a_offset, wn_dim1, wn_offset, wp_dim1, wp_offset, aa_dim1,
	     aa_offset;

    wn_dim1 = *n;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    aa_dim1 = *p;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    wp_dim1 = *p;
    wp_offset = 1 + wp_dim1;
    wp -= wp_offset;
    a_dim1 = *p;
    a_offset = 1 + a_dim1;
    a -= a_offset;
 
    fna_(n, p, &a[a_offset], y, rhs, d, u, beta, eps, &wn[wn_dim1 + 1], 
	 &wn[(wn_dim1 << 1) + 1], &wp[wp_dim1 + 1], &wn[wn_dim1 * 3 + 1], 
	 &wn[(wn_dim1 << 2) + 1], &wn[wn_dim1 * 5 + 1], &wn[wn_dim1 * 6 + 1], 
	 &wp[(wp_dim1 << 1) + 1], &wn[wn_dim1 * 7 + 1], &wn[(wn_dim1 << 3) + 1], 
	 &wn[wn_dim1 * 9 + 1], &wn[wn_dim1 * 10 + 1], &wp[wp_dim1 * 3 + 1], 
	 &wp[(wp_dim1 << 2) + 1], &aa[aa_offset], nit, info);

    return 0;
} 
