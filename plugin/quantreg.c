/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/* Please see the additional license, rq/uiuc-ncsa.txt: this license
   governs the ratfor ("rational fortran") code, authored by Roger
   Koenker, which is called by this module.
*/

#include "libgretl.h"
#include "gretl_f2c.h"

extern void rqfn_ (integer *n, integer *p, double *a, double *y,
		   double *rhs, double *d, double *u, double *beta,
		   double *eps, double *wn, double *wp, double *aa,
		   integer *nit, integer *info);

extern void rqfnb_ (integer *n, integer *p, double *a, double *y,
		    double *rhs, double *d, double *u, double *beta,
		    double *eps, double *wn, double *wp, 
		    integer *nit, integer *info);

#if 0 /* OMG! */
extern void rqbr_ (integer *m, integer *nn, integer *m5, integer *n3,
		   integer *n4, double *a, double *b, double *t, double *toler,
		   integer *ift, double *x, double *e, integer *s,
		   double *wa, double *wb, integer *nsol, integer *ndsol,
		   double *sol, double *dsol, integer *lsol, integer *h,
		   double *qn, double *cutoff, double *ci, double *tnmat,
		   double *big, logical lci1);
#endif

/* 
    z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n +
        5), as.integer(p + 3), as.integer(p + 4), as.double(x),
        as.double(y), as.double(tau), as.double(tol), flag = as.integer(1),
        coef = double(p), resid = double(n), integer(n), double((n +
            5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol),
        sol = double((p + 3) * nsol), dsol = double(n * ndsol),
        lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn),
        cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 *
            p), as.double(big), as.logical(lci1), PACKAGE = "quantreg")
*/

static double rq_loglik (MODEL *pmod, double tau)
{
    double R = 0.0;
    int n = pmod->nobs;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	R += pmod->uhat[t] * (tau - (pmod->uhat[t] < 0));
    }

    return n * (log(tau * (1-tau)) - 1 - log(R/n));
}

#if 0
static double rq_AIC (int n, double ll)
{
    int k = log(n);

    return AIC(ll, k);
}
#endif

/* Bandwidth selection for sparsity estimation, as per
   Hall and Sheather (1988, JRSS(B)); rate = O(n^{-1/3})
*/

static double rq_bandwidth (double tau, int n, double alpha)
{
    double x0 = normal_cdf_inverse(tau);
    double f0 = normal_pdf(x0);
    double b1 = pow(n, -1 / 3.0);
    double b2 = pow(normal_cdf_inverse(1 - alpha/2), 2 / 3.0);
    double b3 = (1.5 * f0 * f0) / (2 * x0 * x0 + 1);

    return b1 * b2 * pow(b3, 1 / 3.0);
}

static void rq_transcribe_results (MODEL *pmod, 
				   const gretl_matrix *y, 
				   double tau,
				   const double *b, int p,
				   const double *u, int n)
{
    double SAR = 0.0;
    int i, t;

    for (i=0; i<p; i++) {
	pmod->coeff[i] = b[i];
    }

    pmod->ess = 0.0;

    for (i=0, t=pmod->t1; i<n; i++, t++) {
	pmod->uhat[t] = u[i];
	pmod->yhat[t] = y->val[i] - u[i];
	SAR += fabs(u[i]);
	pmod->ess += u[i] * u[i];
    }

    /* sum of absolute residuals */
    gretl_model_set_double(pmod, "ladsum", SAR);

    /* set ess-based stats to missing value */
    pmod->rsq = NADBL;
    pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;

    /* LaPlace errors: equivalent of standard error is sum of
       absolute residuals over nobs */
    pmod->sigma = SAR / pmod->nobs; 
    pmod->lnL = rq_loglik(pmod, tau);
    pmod->ci = LAD;
}

static void workspace_init (const gretl_matrix *XT,
			    double tau, double *rhs,
			    double *d, double *u,
			    double *wn)
{
    double xsum;
    int p = XT->rows;
    int n = XT->cols;
    int n10 = n * 10;
    int i, t;

    for (i=0; i<p; i++) {
	xsum = 0.0;
	for (t=0; t<n; t++) {
	    xsum += gretl_matrix_get(XT, i, t);
	}
	rhs[i] = tau * xsum;
    }  

    for (i=0; i<n; i++) {
	d[i] = u[i] = 1.0;
	wn[i] = tau;
    }

    for (i=n; i<n10; i++) {
	wn[i] = 0.0;
    }
}

static int rq_run (gretl_matrix *y, gretl_matrix *XT, 
		   double tau, MODEL *pmod)
{
    integer n = y->rows;
    integer p = XT->rows;
    integer info, nit[3] = {0};
    double *rhs = NULL, *d = NULL, *u = NULL; 
    double *wn = NULL, *wp = NULL, *aa = NULL;
    double beta = .99995;
    double h, eps = 1.0e-7;
    int err = 0;

    /* workspace */
    rhs = malloc(p * sizeof *rhs);
    d =   malloc(n * sizeof *d);
    u =   malloc(n * sizeof *u);
    wn =  malloc(n * 10 * sizeof *wn);
    wp =  malloc(p * (p + 4) * sizeof *wp);
    aa =  malloc(p * p * sizeof *aa);

    if (rhs == NULL || d == NULL || u == NULL || wn == NULL ||
	wp == NULL || aa == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    workspace_init(XT, tau, rhs, d, u, wn);

    /* get coefficients (in wp) and residuals (in wn) */
    rqfn_(&n, &p, XT->val, y->val, rhs, d, u, &beta, &eps, wn, wp, 
	  aa, nit, &info);

    if (info != 0) {
	fprintf(stderr, "rqfn gave info = %d\n", info);
	err = E_DATA;
	goto bailout;
    }

    /* save coeffs, residuals, etc. */
    rq_transcribe_results(pmod, y, tau, wp, p, wn, n);

    /* bandwidth for variance calculation */
    h = rq_bandwidth(tau, n, 0.05);

    if (tau + h > 1) {
	fprintf(stderr, "rq variance: tau + h > 1\n");
	err = 1;
    } else if (tau - h < 0) {
	fprintf(stderr, "rq variance: tau - h < 0\n");
	err = 1;
    }

    if (!err) {
	/* compute covariance matrix via "rank inversion" */
	gretl_matrix *p1 = NULL;
	gretl_matrix *dyhat = NULL;
	gretl_matrix *fX = NULL;
	gretl_matrix *fXX = NULL;
	gretl_matrix *XTX = NULL;
	gretl_matrix *V = NULL;
	double macheps = 2.0e-16;
	double tau_h, x, eps23;
	int i, t;

	eps23 = pow(macheps, 2/3.0);

	p1 =    gretl_matrix_alloc(p, 1);
	dyhat = gretl_matrix_alloc(n, 1);
	fX =    gretl_matrix_alloc(p, n);
	fXX =   gretl_matrix_alloc(p, p);
	XTX =   gretl_matrix_alloc(p, p);
	V =     gretl_matrix_alloc(p, p);

	if (p1 == NULL || dyhat == NULL || fX == NULL ||
	    fXX == NULL || XTX == NULL || V == NULL) {
	    err = E_ALLOC;
	    goto vcv_bailout;
	}
	
	tau_h = tau + h;
	workspace_init(XT, tau_h, rhs, d, u, wn);
	rqfnb_(&n, &p, XT->val, y->val, rhs, d, u, &beta, &eps, wn, wp, 
	       nit, &info);
	if (info != 0) {
	    fprintf(stderr, "tau + h: info = %d\n", info);
	    goto vcv_bailout;
	}

	/* grab coeffs */
	for (i=0; i<p; i++) {
	    p1->val[i] = wp[i];
	}

	tau_h = tau - h;
	workspace_init(XT, tau_h, rhs, d, u, wn);
	rqfnb_(&n, &p, XT->val, y->val, rhs, d, u, &beta, &eps, wn, wp, 
	       nit, &info);
	if (info != 0) {
	    fprintf(stderr, "tau - h: info = %d\n", info);
	    goto vcv_bailout;
	}

	/* diff coeffs */
	for (i=0; i<p; i++) {
	    p1->val[i] -= wp[i];
	}

	gretl_matrix_multiply_mod(XT, GRETL_MOD_TRANSPOSE,
				  p1, GRETL_MOD_NONE,
				  dyhat, GRETL_MOD_NONE);

	for (t=0; t<n; t++) {
	    if (dyhat->val[t] <= 0) {
		fprintf(stderr, "dyhat[%d] is non-positive\n", t);
	    }
	    x = (2.0 * h) / (dyhat->val[t] - eps23);
	    dyhat->val[t] = (x > 0)? sqrt(x) : 0;
	    dyhat->val[t] = (x > 0)? x : 0; 
	    for (i=0; i<p; i++) {
		x = gretl_matrix_get(XT, i, t);
		x *= dyhat->val[t];
		gretl_matrix_set(fX, i, t, x);
	    }
	}

	gretl_matrix_multiply_mod(fX, GRETL_MOD_NONE,
				  XT, GRETL_MOD_TRANSPOSE,
				  fXX, GRETL_MOD_NONE);

	gretl_invert_symmetric_matrix(fXX);

	gretl_matrix_multiply_mod(XT, GRETL_MOD_NONE,
				  XT, GRETL_MOD_TRANSPOSE,
				  XTX, GRETL_MOD_NONE);

	gretl_matrix_qform(fXX, GRETL_MOD_NONE,
			   XTX, V, GRETL_MOD_NONE);

	gretl_matrix_multiply_by_scalar(V, tau * (1 - tau));

#if 0
	gretl_matrix_print(V, "VCV");
#endif

	gretl_model_write_vcv(pmod, V);

    vcv_bailout:

	gretl_matrix_free(p1);
	gretl_matrix_free(dyhat);
	gretl_matrix_free(fX);
	gretl_matrix_free(fXX);
	gretl_matrix_free(XTX);
	gretl_matrix_free(V);
    }

 bailout:

    free(rhs);
    free(d);
    free(u);
    free(wn);
    free(wp);
    free(aa);

    return err;
}

static int rq_make_matrices (MODEL *pmod,
			     double **Z, DATAINFO *pdinfo,
			     gretl_matrix **py,
			     gretl_matrix **pXT)
{
    int n = pdinfo->n;
    int p = pmod->list[0] - 1;
    int yno = pmod->list[1];
    gretl_matrix *XT = NULL;
    gretl_matrix *y = NULL;
    int i, t, v;
    int err = 0;

    XT = gretl_matrix_alloc(p, n);
    y = gretl_matrix_alloc(n, 1);

    if (XT == NULL || y == NULL) {
	gretl_matrix_free(y);
	gretl_matrix_free(XT);
	return E_ALLOC;
    }

    for (t=0; t<n; t++) {
	gretl_vector_set(y, t, Z[yno][t]);
    }

    for (i=0; i<p; i++) {
	v = pmod->list[i+2];
	for (t=0; t<n; t++) {
	    gretl_matrix_set(XT, i, t, Z[v][t]);
	}
    }

    if (err) {
	gretl_matrix_free(y);
	gretl_matrix_free(XT);
    } else {
	*py = y;
	*pXT = XT;
    }

    return err;
}

int rq_driver (const char *parm, MODEL *pmod,
	       double **Z, DATAINFO *pdinfo)
{
    gretl_matrix *y = NULL;
    gretl_matrix *XT = NULL;
    double tau = atof(parm); /* FIXME generalize this */
    int err = 0;

    if (tau < .01 || tau > .99) {
	gretl_errmsg_sprintf("quantreg: tau must be >= .01 and <= .99");
	err = E_DATA;
    }

    fprintf(stderr, "tau = %g\n", tau);
    printlist(pmod->list, "list in quantreg_driver");

    if (!err) {
	err = rq_make_matrices(pmod, Z, pdinfo, &y, &XT);
    }

    if (!err) {
	err = rq_run(y, XT, tau, pmod);
    }

    gretl_matrix_free(y);
    gretl_matrix_free(XT);

    if (err && pmod->errcode == 0) {
	pmod->errcode = err;
    }

    return err;
}
