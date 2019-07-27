/* rqfnb.f -- translated by f2c (version 20050501)
   and slightly cleaned up by Allin Cottrell
*/

#include "libgretl.h"
#include "gretl_f2c.h"

/* Table of constant values */

static double c_b4 = 1.;
static integer one = 1;
static double zero = 0.;
static double c_b13 = -1.;

/* lapack/blas functions called below */

extern int dsyr_ (char *, integer *, double *, double *, 
		  integer *, double *, integer *);

extern int dposv_ (char *, integer *, integer *, double *, integer *, 
		   double *, integer *, integer *);

extern int dgemv_ (char *, integer *, integer *, double *, double *, 
		   integer *, double *, integer *, double *, 
		   double *, integer *);

extern int dcopy_ (integer *, double *, integer *, double *, 
		   integer *);

extern int dswap_ (integer *, double *, integer *, double *, 
		   integer *);

extern int daxpy_ (integer *, double *, double *, integer *, 
		   double *, integer *);

extern int dpotrs_ (char *, integer *, integer *, double *, integer *, 
		    double *, integer *, integer *);

extern double ddot_ (integer *, double *, integer *, double *, 
		     integer *);

static int stepy_ (integer *n, integer *p, double *a, 
		   double *d, double *b, double *ada, 
		   integer *info)
{
    integer i, m = *p * *p;
    int attempt = 0;
    int err = 0;

 try_again:

    for (i=0; i<m; i++) {
	ada[i] = 0.0;
    }

    for (i=0; i<*n; i++) {
	dsyr_("U", p, &d[i], &a[i * *p], &one, ada, p);
    }

    if (attempt == 0) {
	dposv_("U", p, &one, ada, p, b, p, info);
	if (*info != 0) {
	    fprintf(stderr, "stepy: dposv gave info = %d\n", *info);
	    attempt = 1;
	    goto try_again;
	}
    } else {
	gretl_matrix A, B;

	gretl_matrix_init(&A);
	gretl_matrix_init(&B);

	A.rows = A.cols = *p;
	A.val = ada;
	B.rows = *p;
	B.cols = 1;
	B.val = b;

	err = gretl_LU_solve(&A, &B);
	if (err) {
	    fprintf(stderr, "stepy: gretl_LU_solve: err = %d\n", err);
	}
    }

    return err;
}

#define ITERSTEP 5

static int lpfnb_ (integer *n, integer *p, double *a, double *c__, 
		   double *b, double *d__, double *u, double *beta,
		   double *eps, double *x, double *s, double *y, 
		   double *z__, double *w, double *dx, double *ds, 
		   double *dy, double *dz, double *dw, double *dr, 
		   double *rhs, double *ada, integer *nit, integer *info,
		   void (*callback)(void))
{
    integer a_dim1 = *p, ada_dim1 = *p;
    integer a_offset = 1 + a_dim1, ada_offset = 1 + ada_dim1;
    double d1, d2;
    static double g;
    static integer i;
    static double mu, gap;
    static double dsdw, dxdz;
    static double deltad, deltap;
    int main_iters = 0;
    int err = 0;

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
    ada -= ada_offset;
    --rhs;
    --dy;
    --y;
    --b;
    a -= a_offset;
    --nit;

    /* Function Body */
    nit[1] = 0;
    nit[2] = 0;
    nit[3] = *n;
    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &c__[1], &one, &zero, &y[1],
	   &one);
    for (i = 1; i <= *n; ++i) {
	d__[i] = 1.;
    }
    err = stepy_(n, p, &a[a_offset], &d__[1], &y[1], &ada[ada_offset], info);
    if (err) {
	return err;
    }
    dcopy_(n, &c__[1], &one, &s[1], &one);
    dgemv_("T", p, n, &c_b13, &a[a_offset], p, &y[1], &one, &c_b4, &s[1],
	   &one);
    for (i = 1; i <= *n; ++i) {
	if ((d1 = s[i], fabs(d1)) < *eps) {
	    d1 = s[i];
	    z__[i] = max(d1,0.) + *eps;
	    d1 = -s[i];
	    w[i] = max(d1,0.) + *eps;
	} else {
	    d1 = s[i];
	    z__[i] = max(d1, 0.);
	    d1 = -s[i];
	    w[i] = max(d1, 0.);
	}
	s[i] = u[i] - x[i];
    }

    gap = ddot_(n, &z__[1], &one, &x[1], &one) + 
	ddot_(n, &w[1], &one, &s[1], &one);

looptop:

    if (callback != NULL && (main_iters++ % ITERSTEP == 0)) {
	callback();
    }

    if (gap > *eps && nit[1] < 50) {
	++nit[1];
	for (i = 1; i <= *n; ++i) {
	    d__[i] = 1. / (z__[i] / x[i] + w[i] / s[i]);
	    ds[i] = z__[i] - w[i];
	    dz[i] = d__[i] * ds[i];
	}

	dcopy_(p, &b[1], &one, &dy[1], &one);
	dgemv_("N", p, n, &c_b13, &a[a_offset], p, &x[1], &one, &c_b4, &dy[1],
		&one);
	dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dz[1], &one, &c_b4, &dy[1],
		&one);
	dcopy_(p, &dy[1], &one, &rhs[1], &one);
	err = stepy_(n, p, &a[a_offset], &d__[1], &dy[1], &ada[ada_offset], info);
	if (err) {
	    return err;
	}

	dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &one, &c_b13, 
	       &ds[1], &one);
	deltap = 1e20;
	deltad = 1e20;

	for (i = 1; i <= *n; ++i) {
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
	d1 = *beta * deltap;
	deltap = min(d1,1.);
	d1 = *beta * deltad;
	deltad = min(d1,1.);
	if (min(deltap,deltad) < 1.) {
	    ++nit[2];
	    mu = ddot_(n, &x[1], &one, &z__[1], &one) + 
		ddot_(n, &s[1], &one, &w[1], &one);
	    g = mu + deltap * ddot_(n, &dx[1], &one, &z__[1], &one) + 
		deltad * ddot_(n, &dz[1], &one, &x[1], &one) + 
		deltap * deltad * ddot_(n, &dz[1], &one, &dx[1], &one) + 
		deltap * ddot_(n, &ds[1], &one, &w[1], &one) + 
		deltad * ddot_(n, &dw[1], &one, &s[1], &one) + 
		deltap * deltad * ddot_(n, &ds[1], &one, &dw[1], &one);
	    d1 = g / mu;
	    mu = mu * (d1 * (d1 * d1)) / (double) (*n << 1);
	    for (i = 1; i <= *n; ++i) {
		dr[i] = d__[i] * (mu * (1 / s[i] - 1 / x[i]) + dx[i]
				  * dz[i] / x[i] - ds[i] * dw[i] / s[i]);
	    }
	    dswap_(p, &rhs[1], &one, &dy[1], &one);
	    dgemv_("N", p, n, &c_b4, &a[a_offset], p, &dr[1], &one, &c_b4,
		   &dy[1], &one);
	    dpotrs_("U", p, &one, &ada[ada_offset], p, &dy[1], p, info);
	    if (*info != 0) {
		fprintf(stderr, "lpfnb: dpotrs_ gave info = %d\n", *info);
	    }
	    dgemv_("T", p, n, &c_b4, &a[a_offset], p, &dy[1], &one, &zero, 
		   &u[1], &one);
	    deltap = 1e20;
	    deltad = 1e20;
	    for (i = 1; i <= *n; ++i) {
		dxdz = dx[i] * dz[i];
		dsdw = ds[i] * dw[i];
		dx[i] = d__[i] * (u[i] - z__[i] + w[i]) - dr[i];
		ds[i] = -dx[i];
		dz[i] = -z__[i] + (mu - z__[i] * dx[i] - dxdz) / x[i];
		dw[i] = -w[i] + (mu - w[i] * ds[i] - dsdw) / s[i];
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
	}

	daxpy_(n, &deltap, &dx[1], &one, &x[1], &one);
	daxpy_(n, &deltap, &ds[1], &one, &s[1], &one);
	daxpy_(p, &deltad, &dy[1], &one, &y[1], &one);
	daxpy_(n, &deltad, &dz[1], &one, &z__[1], &one);
	daxpy_(n, &deltad, &dw[1], &one, &w[1], &one);
	gap = ddot_(n, &z__[1], &one, &x[1], &one) + 
	    ddot_(n, &w[1], &one, &s[1], &one);

	goto looptop;
    }

    daxpy_(n, &c_b13, &w[1], &one, &z__[1], &one);
    dswap_(n, &z__[1], &one, &x[1], &one);

    return err;
} /* end of lpfnb_ */

int rqfnb_ (integer *n, integer *p, double *a, double *y, 
	    double *rhs, double *d, double *u, double *beta,
	    double *eps, double *wn, double *wp, integer *nit, 
	    integer *info, void (*callback)(void))
{
    integer a_dim = *p, wn_dim = *n, wp_dim = *p;
    int err;

    /* Parameter adjustments */
    wn -= 1 + wn_dim;
    wp -= 1 + wp_dim;
    a -= 1 + a_dim;

    err = lpfnb_(n, p, &a[a_dim + 1], y, rhs, d, u, beta, eps, 
		 &wn[wn_dim + 1], &wn[(wn_dim << 1) + 1], &wp[wp_dim + 1], 
		 &wn[wn_dim * 3 + 1], &wn[(wn_dim << 2) + 1], &wn[wn_dim * 5 + 1], 
		 &wn[wn_dim * 6 + 1], &wp[(wp_dim << 1) + 1], &wn[wn_dim * 7 + 1],
		 &wn[(wn_dim << 3) + 1], &wn[wn_dim * 9 + 1], &wp[wp_dim * 3 + 1], 
		 &wp[(wp_dim << 2) + 1], nit, info, callback);

    return err;
} 
