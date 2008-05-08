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

/* Frisch-Newton algorithm */

extern void rqfn_ (integer *n, integer *p, double *a, double *y,
		   double *rhs, double *d, double *u, double *beta,
		   double *eps, double *wn, double *wp, double *aa,
		   integer *nit, integer *info);

extern void rqfnb_ (integer *n, integer *p, double *a, double *y,
		    double *rhs, double *d, double *u, double *beta,
		    double *eps, double *wn, double *wp, 
		    integer *nit, integer *info);

#if 0 /* OMG! */

/* Modified simplex, a la Barrodale-Roberts: this variant lets us get
   1-alpha confidence intervals using the rank-inversion method.
*/

extern void rqbr_ (integer *n,    /* number of observations */
		   integer *p,    /* number of parameters */
		   integer *n5,   /* n + 5 */
		   integer *p3,   /* p + 3 */
		   integer *p4,   /* p + 4 */
		   double *x,     /* matrix of independent vars */
		   double *y,     /* dependent var vector */
		   double *tau,   /* the desired quantile */
		   double *tol,   /* machine precision to the 2/3 */
		   integer *ift,  /* exit code: 0 = OK */
		   double *coeff, /* p-vector of parameter estimates */
		   double *resid, /* n-vector of residuals */
		   integer *s,    /* integer work array, length n */
		   double *wa,    /* real work array, size n5 * p4 */
		   double *wb,    /* and another, size n */
		   integer *nsol, /* estimated (row) dim. of primal solution array,
				     "say 3*n; minimum value = 2" */
		   integer *ndsol, /* estimated (row) dim. of dual solution array,
				      "say 3*n; minimum value = 2" */
		   double *sol,   /* primal solution array, size (p + 3) * nsol (?) */
		   double *dsol,  /* dual solution array, size n * ndsol */
		   integer *lsol, /* actual dimension of the solution arrays */
		   integer *h,    /* matrix of basic observations indices, size p * nsol */
		   double *qn,    /* vector of residual variances from projection of 
				     each column of x on the remaining columns */
		   double *cutoff,  /* critical point for N(0,1) */
		   double *ci,      /* matrix of confidence intervals, size 4 * p */
		   double *tnmat,   /* matrix of JGPK rank test statistics */
		   double *big,     /* "large positive finite floating-point number" */
		   logical *lci1);  /* do confidence intervals? */

/* 
   Koenker he say:

     Utilization:  If you just want a solution at a single quantile you
     needn't bother with sol, nsol, etc.  If you want all the solutions
     then set theta (what?) to something <0 and sol and dsol will return 
     all the estimated quantile solutions.

     The algorithm is a slightly modified version of algorithm as 229
     described in Koenker and Dorey, Computing Regression Quantiles,
     Applied Statistics, pp. 383-393.
*/

static int rq_fit_br (gretl_matrix *x, gretl_matrix *y, double tau, 
		      double alpha, /* for confidence intervals, default 0.1 */
		      int do_ci,    /* doing c.i.s? */
		      int iid,      /* assume i.i.d. errors? (R default yes) */
		      int tcrit)    /* use t rather than normal (R default yes) */
{
    integer p = gretl_matrix_cols(x);
    integer n = gretl_matrix_rows(x);
    integer n5 = n + 5;
    integer p3 = p + 3;
    integer p4 = p + 4;
    integer ift, *s = NULL, *h = NULL;
    double *coeff = NULL;
    double *resid = NULL;
    double *wa = NULL;
    double *wb = NULL;
    double *qn = NULL;
    double *sol = NULL;
    double *dsol = NULL;
    gretl_matrix *tnmat = NULL;
    gretl_matrix *ci = NULL;
    double tol = 1.0e-12; /* Machine$double.eps^(2/3); */
    double eps = tol;
    double big = NADBL;
    double cutoff = 0.0;
    integer nsol = 2;
    integer ndsol = 2;
    integer lsol = 0;
    logical lci1 = 0;
    int i;

    if (tau < 0 || tau > 1) {
	/* We'll not bother with the "whole tau process" */
	return 1;
    }

    if (p == 1) {
	/* just one indep var: don't do confidence intervals */
	do_ci = 0;
    }

    coeff = malloc(p * sizeof *coeff);
    resid = malloc(n * sizeof *resid);
    s     = malloc(n * sizeof *s);
    wa    = malloc(n5 * p4 * sizeof *wa);
    wb    = malloc(n * sizeof *wb);
    sol   = malloc(p3 * nsol * sizeof *sol); /* ?? */
    dsol  = malloc(n * ndsol * sizeof *dsol); /* ?? */
    qn    = malloc(p * sizeof *qn);
    h     = malloc(p * nsol * sizeof *h); /* ?? */

    ci = gretl_matrix_alloc(4, p);
    tnmat = gretl_matrix_alloc(4, p);

    for (i=0; i<p; i++) {
	qn[i] = 0.0;
    }

    if (do_ci) {
	/* doing confidence intervals */
	lci1 = 1;

#if 0 /* translate! */
	if (tcrit) {
	    /* Student t critical values */
	    cutoff = qt(1 - alpha/2, n - p);
	} else {
	    /* Normal critical values */
	    cutoff = qnorm(1 - alpha/2);
	}
#endif

	if (iid) {
	    /* assuming iid errors */
	    gretl_matrix *xtx = gretl_matrix_alloc(p, p);

	    gretl_matrix_multiply_mod(x, GRETL_MOD_TRANSPOSE,
				      x, GRETL_MOD_NONE,
				      xtx, GRETL_MOD_NONE);
	    for (i=0; i<p; i++) {
		qn[i] = 1 / gretl_matrix_get(xtx, i, i);
	    }
	    gretl_matrix_free(xtx);
	} else {
	    /* allowing for non-iid errors (more translation needed!) */
	    gretl_matrix *bdiff;
	    gretl_matrix *dyhat;
#if 0
	    double h = rq_bandwidth(tau, n, hs = TRUE);
#endif
	    double h, dyi;

	    bdiff = gretl_matrix_alloc(p, 1);
	    dyhat = gretl_matrix_alloc(n, 1);

	    rq_fit_br(x, y, tau + h, 0, 0, 0, 0);
	    for (i=0; i<p; i++) {
		bdiff->val[i] = coeff[i];
	    }

	    rq_fit_br(x, y, tau - h, 0, 0, 0, 0);
	    for (i=0; i<p; i++) {
		bdiff->val[i] -= coeff[i];
	    }

	    gretl_matrix_multiply(x, bdiff, dyhat);

#if 0
	    if (any(dyhat <= 0)) {
		pfis = (100 * sum(dyhat <= 0)) / n;
		warning(paste(pfis, "percent fis <= 0"));
	    }
#endif

	    for (i=0; i<n; i++) {
		dyi = (2 * h) / (dyhat->val[i] - eps);
		dyhat->val[i] = (dyi < eps)? eps : dyi;
	    }

#if 0
	    for (j=1; j<=p; j++) {
		qnj = lm(x[, j] ~ x[, -j] - 1, weights = f)$resid;
		qn[j] = sum(qnj * qnj);
	    }
#endif
	} 
    } 

    rqbr_(&n, &p, &n5, &p3, &p4, x->val, y->val, &tau, &tol, &ift,
	  coeff, resid, s, wa, wb, &nsol, &ndsol, sol, dsol,
	  &lsol, h, qn, &cutoff, ci->val, tnmat->val, &big, &lci1);

    if (ift != 0) {
	if (1) { /* ?? */
	    fprintf(stderr, "Solution may be nonunique\n");
	} else { 
	    fprintf(stderr, "Premature end - possible conditioning problem in x\n");
	}
    }

    if (do_ci) {
	/* interpolated confidence intervals */
	double c1j, c2j, c3j, c4j;
	double tn1, tn2, tn3, tn4;
	int j;

	for (j=0; j<p; j++) {
	    c1j = gretl_matrix_get(ci, 0, j);
	    c2j = gretl_matrix_get(ci, 1, j);
	    c3j = gretl_matrix_get(ci, 2, j);
	    c4j = gretl_matrix_get(ci, 3, j);
	    tn1 = gretl_matrix_get(tnmat, 0, j);
	    tn2 = gretl_matrix_get(tnmat, 1, j);
	    tn3 = gretl_matrix_get(tnmat, 2, j);
	    tn4 = gretl_matrix_get(tnmat, 3, j);

	    c3j += fabs(c4j - c3j) * (cutoff - fabs(tn3)) / fabs(tn4 - tn3);
	    c2j -= fabs(c1j - c2j) * (cutoff - fabs(tn2)) / fabs(tn1 - tn2);

	    gretl_matrix_set(ci, 2, j, c3j);
	    gretl_matrix_set(ci, 1, j, c2j);
	}

	/* look out for NAs? */

	/* intervals now in rows 1 and 2? */
#if 0
        coefficients = cbind(coef, t(Tci[2:3, ]));
	cnames = c("coefficients", "lower bd", "upper bd");
        residuals = y - x %*% coef;
#endif
    } 

    free(coeff);
    free(resid);
    free(s);
    free(wa);
    free(wb);
    free(sol);
    free(dsol);
    free(qn);
    free(h);

    gretl_matrix_free(ci);
    gretl_matrix_free(tnmat);

    return 0;
}

#endif

struct rq_info {
    integer n, p;
    double tau;
    double beta;
    double eps;
    double *rhs;
    double *d;
    double *u;
    double *wn;
    double *wp;
    double *aa;
    integer nit[3];
    integer info;
};

static void rq_info_free (struct rq_info *rq)
{
    free(rq->rhs);
    free(rq->d);
    free(rq->u);
    free(rq->wn);
    free(rq->wp);
    free(rq->aa);
}

static int rq_info_alloc (struct rq_info *rq, int n, int p,
			  double tau)
{
    int err = 0;

    rq->rhs = NULL;
    rq->d =   NULL;
    rq->u =   NULL;
    rq->wn =  NULL;
    rq->wp =  NULL;
    rq->aa =  NULL;

    rq->rhs = malloc(p * sizeof *rq->rhs);
    rq->d =   malloc(n * sizeof *rq->d);
    rq->u =   malloc(n * sizeof *rq->u);
    rq->wn =  malloc(n * 10 * sizeof *rq->wn);
    rq->wp =  malloc(p * (p + 4) * sizeof *rq->wp);
    rq->aa =  malloc(p * p * sizeof *rq->aa);

    if (rq->rhs == NULL || rq->d == NULL ||
	rq->u == NULL || rq->wn == NULL ||
	rq->wp == NULL || rq->aa == NULL) {
	rq_info_free(rq);
	err = E_ALLOC;
    } else {
	rq->n = n;
	rq->p = p;
	rq->tau = tau;
	rq->beta = .99995;
	rq->eps = 1.0e-7;
    }

    return err;
}

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
				   struct rq_info *rq)
{
    double *b = rq->wp;
    double *u = rq->wn;
    double SAR = 0.0;
    int i, t;

    for (i=0; i<rq->p; i++) {
	pmod->coeff[i] = b[i];
    }

    pmod->ess = 0.0;

    for (i=0, t=pmod->t1; i<rq->n; i++, t++) {
	pmod->uhat[t] = u[i];
	pmod->yhat[t] = y->val[i] - u[i];
	SAR += fabs(u[i]);
	pmod->ess += u[i] * u[i];
    }

    gretl_model_set_double(pmod, "tau", rq->tau);

    /* sum of absolute residuals */
    gretl_model_set_double(pmod, "ladsum", SAR);

    /* set ess-based stats to missing value */
    pmod->rsq = NADBL;
    pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;

    /* LaPlace errors: equivalent of standard error is sum of
       absolute residuals over nobs */
    pmod->sigma = SAR / pmod->nobs; 
    pmod->lnL = rq_loglik(pmod, rq->tau);
    mle_criteria(pmod, 0);
    pmod->ci = LAD;
}

static void workspace_init (const gretl_matrix *XT,
			    double tau, 
			    struct rq_info *ws)
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
	ws->rhs[i] = tau * xsum;
    }  

    for (i=0; i<n; i++) {
	ws->d[i] = ws->u[i] = 1.0;
	ws->wn[i] = tau;
    }

    for (i=n; i<n10; i++) {
	ws->wn[i] = 0.0;
    }
}

static int rq_VCV (MODEL *pmod, gretl_matrix *y,
		   gretl_matrix *XT, double tau,
		   struct rq_info *rq)
{
    gretl_matrix *p1 = NULL;
    gretl_matrix *dyhat = NULL;
    gretl_matrix *fX = NULL;
    gretl_matrix *fXX = NULL;
    gretl_matrix *XTX = NULL;
    gretl_matrix *V = NULL;
    double macheps = 2.0e-16;
    double h, tau_h, x, eps23;
    int n = rq->n;
    int p = rq->p;
    int i, t, err = 0;

    /* bandwidth for variance calculation */
    h = rq_bandwidth(tau, rq->n, 0.05);

    if (tau + h > 1) {
	fprintf(stderr, "rq variance: tau + h > 1\n");
	return E_DATA;
    } else if (tau - h < 0) {
	fprintf(stderr, "rq variance: tau - h < 0\n");
	return E_DATA;
    }    

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
	goto bailout;
    }
	
    tau_h = tau + h;
    workspace_init(XT, tau_h, rq);
    rqfnb_(&n, &p, XT->val, y->val, rq->rhs, rq->d, rq->u, &rq->beta, &rq->eps, 
	   rq->wn, rq->wp, rq->nit, &rq->info);
    if (rq->info != 0) {
	fprintf(stderr, "tau + h: info = %d\n", rq->info);
	err = E_DATA;
	goto bailout;
    }

    /* grab coeffs */
    for (i=0; i<p; i++) {
	p1->val[i] = rq->wp[i];
    }

    tau_h = tau - h;
    workspace_init(XT, tau_h, rq);
    rqfnb_(&n, &p, XT->val, y->val, rq->rhs, rq->d, rq->u, &rq->beta, &rq->eps, 
	   rq->wn, rq->wp, rq->nit, &rq->info);
    if (rq->info != 0) {
	fprintf(stderr, "tau - h: info = %d\n", rq->info);
	err = E_DATA;
	goto bailout;
    }

    /* diff coeffs */
    for (i=0; i<p; i++) {
	p1->val[i] -= rq->wp[i];
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

    gretl_model_write_vcv(pmod, V);

 bailout:

    gretl_matrix_free(p1);
    gretl_matrix_free(dyhat);
    gretl_matrix_free(fX);
    gretl_matrix_free(fXX);
    gretl_matrix_free(XTX);
    gretl_matrix_free(V);

    return err;
}

static int rq_run (gretl_matrix *y, gretl_matrix *XT, 
		   double tau, MODEL *pmod)
{
    struct rq_info rq;
    integer n = y->rows;
    integer p = XT->rows;
    int err = 0;

    err = rq_info_alloc(&rq, n, p, tau);
    if (err) {
	return err;
    }

    workspace_init(XT, tau, &rq);

    /* get coefficients (in wp) and residuals (in wn) */
    rqfn_(&n, &p, XT->val, y->val, rq.rhs, rq.d, rq.u, &rq.beta, &rq.eps, 
	  rq.wn, rq.wp, rq.aa, rq.nit, &rq.info);

    if (rq.info != 0) {
	fprintf(stderr, "rqfn gave info = %d\n", rq.info);
	err = E_DATA;
    }

    if (!err) {
	/* save coeffs, residuals, etc., before computing VCV */
	rq_transcribe_results(pmod, y, &rq);
    }

    if (!err) {
	/* compute and add covariance matrix via "rank inversion" */
	err = rq_VCV(pmod, y, XT, tau, &rq);
    }

    rq_info_free(&rq);

    return err;
}

/* Write y and X from pmod into gretl matrices, respecting the sample
   range.  Note that for the purposes of the underlying rq functions
   we want the X matrix in transpose form (XT).
*/

static int rq_make_matrices (MODEL *pmod,
			     double **Z, DATAINFO *pdinfo,
			     gretl_matrix **py,
			     gretl_matrix **pXT)
{
    int n = pmod->nobs;
    int p = pmod->ncoeff;
    int yno = pmod->list[1];
    gretl_matrix *XT = NULL;
    gretl_matrix *y = NULL;
    int i, s, t, v;
    int err = 0;

    XT = gretl_matrix_alloc(p, n);
    y = gretl_matrix_alloc(n, 1);

    if (XT == NULL || y == NULL) {
	gretl_matrix_free(y);
	gretl_matrix_free(XT);
	return E_ALLOC;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	gretl_vector_set(y, s++, Z[yno][t]);
    }

    for (i=0; i<p; i++) {
	v = pmod->list[i+2];
	s = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    gretl_matrix_set(XT, i, s++, Z[v][t]);
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
    gretl_matrix *X = NULL;
    double tau = atof(parm); /* FIXME generalize this */
    int err = 0;

    if (tau < .01 || tau > .99) {
	gretl_errmsg_sprintf("quantreg: tau must be >= .01 and <= .99");
	err = E_DATA;
    }

    if (!err) {
	err = rq_make_matrices(pmod, Z, pdinfo, &y, &X);
    }

    if (!err) {
	err = rq_run(y, X, tau, pmod);
    }

    if (!err) {
	gretl_model_add_y_median(pmod, Z[pmod->list[1]]);
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);

    if (err && pmod->errcode == 0) {
	pmod->errcode = err;
    }

    return err;
}
