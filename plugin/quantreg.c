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
   governs the Fortran code, authored by Roger Koenker, which is
   called by this module.
*/

#include "libgretl.h"
#include "gretl_f2c.h"

#include <errno.h>

#define QDEBUG 1

/* Frisch-Newton algorithm: we use this if we're not computing
   rank-inversion confidence intervals.
*/

extern int rqfnb_ (integer *n, integer *p, double *a, double *y,
		   double *rhs, double *d, double *u, double *beta,
		   double *eps, double *wn, double *wp, 
		   integer *nit, integer *info);

/* Modified simplex, a la Barrodale-Roberts: this variant lets us get
   1-alpha confidence intervals using the rank-inversion method.
*/

extern int rqbr_ (integer *n,    /* number of observations */
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
		  double *tnmat,   /* matrix of JGPK rank-test statistics */
		  double *big,     /* "large positive finite floating-point number" */
		  logical *lci1);  /* do confidence intervals? */

static int rq_fit_br (gretl_matrix *y, gretl_matrix *X,
		      double tau, double alpha, gretlopt opt,
		      MODEL *pmod, gretl_matrix *theta);

#define calc_eps23 (pow(2.0e-16, 2/3.0))

/* Bandwidth selection for sparsity estimation, as per
   Hall and Sheather (1988, JRSS(B)); rate = O(n^{-1/3})
*/

static double rq_bandwidth (double tau, int n)
{
    double alpha = 0.05;
    double x0 = normal_cdf_inverse(tau);
    double f0 = normal_pdf(x0);
    double b1 = pow(n, -1 / 3.0);
    double b2 = pow(normal_cdf_inverse(1 - alpha/2), 2 / 3.0);
    double b3 = (1.5 * f0 * f0) / (2 * x0 * x0 + 1);

    return b1 * b2 * pow(b3, 1 / 3.0);
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

static void rq_transcribe_results (MODEL *pmod, 
				   const gretl_matrix *y,
				   double tau,
				   const double *b,
				   const double *u)
{
    double SAR = 0.0;
    int n = y->rows;
    int i, t;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = b[i];
	pmod->sderr[i] = NADBL; /* will have to be replaced */
    }

    pmod->ess = 0.0;

    for (i=0, t=pmod->t1; i<n; i++, t++) {
	pmod->uhat[t] = u[i];
	pmod->yhat[t] = y->val[i] - u[i];
	SAR += fabs(u[i]);
	pmod->ess += u[i] * u[i];
    }

    gretl_model_set_int(pmod, "rq", 1);
    gretl_model_set_double(pmod, "tau", tau);

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
    mle_criteria(pmod, 0);
    pmod->ci = LAD;
}

/* extract the interpolated lower and upper bounds from ci into
   a new matrix and attach this to the model for printing 
*/

static int rq_attach_intervals (MODEL *pmod, gretl_matrix *ci,
				double alpha, gretlopt opt)
{
    gretl_matrix *rqci;
    double blo, bhi;
    int i, p = ci->cols;

    rqci = gretl_matrix_alloc(p, 2);
    if (rqci == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<p; i++) {
	blo = gretl_matrix_get(ci, 1, i);
	bhi = gretl_matrix_get(ci, 2, i);
	gretl_matrix_set(rqci, i, 0, blo);
	gretl_matrix_set(rqci, i, 1, bhi);
    }

    gretl_model_set_matrix_as_data(pmod, "coeff_intervals", rqci);
    gretl_model_set_double(pmod, "rq_alpha", alpha);

    if (opt & OPT_R) {
	gretl_model_set_int(pmod, "rq_nid", 1);
    }

    return 0;
}

static void rq_interpolate_intervals (gretl_matrix *ci,
				      gretl_matrix *tn,
				      double cut)
{
    double c1j, c2j, c3j, c4j;
    double tn1, tn2, tn3, tn4;
    int p = gretl_matrix_cols(ci);
    int j;

#if QDEBUG > 1
    gretl_matrix_print(ci, "ci (original, from rqbr)");
    gretl_matrix_print(tn, "tnmat");
#endif

    for (j=0; j<p; j++) {
	c1j = gretl_matrix_get(ci, 0, j);
	c2j = gretl_matrix_get(ci, 1, j);
	c3j = gretl_matrix_get(ci, 2, j);
	c4j = gretl_matrix_get(ci, 3, j);
	tn1 = gretl_matrix_get(tn, 0, j);
	tn2 = gretl_matrix_get(tn, 1, j);
	tn3 = gretl_matrix_get(tn, 2, j);
	tn4 = gretl_matrix_get(tn, 3, j);

	c3j += fabs(c4j - c3j) * (cut - fabs(tn3)) / fabs(tn4 - tn3);
	c2j -= fabs(c1j - c2j) * (cut - fabs(tn2)) / fabs(tn1 - tn2);

	/* Write the 1-alpha intervals into rows 1 and 2 
	   of the matrix ci */

	gretl_matrix_set(ci, 2, j, c3j);
	gretl_matrix_set(ci, 1, j, c2j);
    }
}

/* Robust version of rank-inversion confidence interval 
   calculation
*/

static int make_nid_qn (gretl_matrix *y, gretl_matrix *X, 
			double *qn, double tau, double eps)
{
    int n = gretl_matrix_rows(X);
    int p = gretl_matrix_cols(X);
    gretl_matrix *coeff = NULL;
    gretl_matrix *bdiff = NULL;
    gretl_matrix *f = NULL;
    gretl_matrix *xj = NULL;
    gretl_matrix *uj = NULL;
    gretl_matrix *Xcmp = NULL;
    double h, fi, ui, xik;
    int i, j, k, jj;
    int badf = 0;
    int err = 0;

    h = rq_bandwidth(tau, n);

    coeff = gretl_matrix_alloc(p, 1);
    bdiff = gretl_matrix_alloc(p, 1);
    f     = gretl_matrix_alloc(n, 1);
    xj    = gretl_matrix_alloc(n, 1);
    uj    = gretl_matrix_alloc(n, 1);
    Xcmp  = gretl_matrix_alloc(n, p - 1);

    if (coeff == NULL || bdiff == NULL || f == NULL ||
	xj == NULL || uj == NULL || Xcmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    rq_fit_br(y, X, tau + h, 0, OPT_NONE, NULL, coeff);
    for (j=0; j<p; j++) {
	bdiff->val[j] = coeff->val[j];
    }

    rq_fit_br(y, X, tau - h, 0, OPT_NONE, NULL, coeff);
    for (j=0; j<p; j++) {
	bdiff->val[j] -= coeff->val[j];
    }

#if QDEBUG
    fprintf(stderr, "make_nid_qn: bandwidth = %g\n", h);
    gretl_matrix_print(bdiff, "coeff diffs");
#endif

    gretl_matrix_multiply(X, bdiff, f);

    for (i=0; i<n; i++) {
	if (f->val[i] <= 0) {
	    badf++;
	}
    }

    if (badf > 0) {
	fprintf(stderr, "Warning: %g percent of fi's <= 0\n",
		100 * (double) badf / n);
    }

    for (i=0; i<n; i++) {
	fi = (2 * h) / (f->val[i] - eps);
	fi = (fi < eps)? eps : fi;
	f->val[i] = sqrt(fi); /* ?? */
    }

    /* Now set each qn[j] to the SSR from an f-weighted regression of 
       X_j on the other regressors.
    */

    gretl_matrix_reuse(coeff, p - 1, 1);

    for (j=0; j<p; j++) {
	for (i=0; i<n; i++) {
	    /* weighted dependent var */
	    xj->val[i] = f->val[i] * gretl_matrix_get(X, i, j);
	}

	jj = 0;
	for (k=0; k<p; k++) {
	    if (k != j) {
		/* weighted independent var */
		for (i=0; i<n; i++) {
		    xik = gretl_matrix_get(X, i, k);
		    gretl_matrix_set(Xcmp, i, jj, f->val[i] * xik);
		}
		jj++;
	    }
	}

	err = gretl_matrix_ols(xj, Xcmp, coeff, NULL, uj, NULL);

	if (!err) {
	    qn[j] = 0.0;
	    for (i=0; i<n; i++) {
		ui = uj->val[i] / f->val[i];
		qn[j] += ui * ui;
	    }
	}
    }

 bailout:

    gretl_matrix_free(coeff);
    gretl_matrix_free(bdiff);
    gretl_matrix_free(f);
    gretl_matrix_free(xj);
    gretl_matrix_free(uj);
    gretl_matrix_free(Xcmp);

    return err;
}

static gretl_matrix *get_XTX_inverse (const gretl_matrix *X, int *err)
{
    int k = min(X->rows, X->cols);
    GretlMatrixMod flag1, flag2;
    gretl_matrix *XTX;

    XTX = gretl_matrix_alloc(k, k);
    if (XTX == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    flag1 = (X->cols == k)? GRETL_MOD_TRANSPOSE : GRETL_MOD_NONE;
    flag2 = (X->cols == k)? GRETL_MOD_NONE : GRETL_MOD_TRANSPOSE;

    *err = gretl_matrix_multiply_mod(X, flag1,
				     X, flag2,
				     XTX, GRETL_MOD_NONE);

    if (!*err) {
	*err = gretl_invert_symmetric_matrix(XTX);
    }

    if (*err) {
	gretl_matrix_free(XTX);
	XTX = NULL;
    }    

    return XTX;
}

/* IID version of rank-inversion confidence interval 
   calculation
*/

static int make_iid_qn (const gretl_matrix *X, double *qn)
{
    int i, p = gretl_matrix_cols(X);
    gretl_matrix *XTX;
    int err = 0;

    XTX = get_XTX_inverse(X, &err);

    if (!err) {
	for (i=0; i<p; i++) {
	    qn[i] = 1 / gretl_matrix_get(XTX, i, i);
	}
	gretl_matrix_free(XTX);
    }

    return err;
}

/* OPT_I for intervals; OPT_N for no df; OPT_R for robust (not iid) */

static int rq_fit_br (gretl_matrix *y, gretl_matrix *X, 
		      double tau, double alpha, gretlopt opt,
		      MODEL *pmod, gretl_matrix *theta)
{
    integer p = gretl_matrix_cols(X);
    integer n = gretl_matrix_rows(X);
    integer n5 = n + 5;
    integer p3 = p + 3;
    integer p4 = p + 4;
    integer *s, *h, *ispace = NULL;
    double *rspace = NULL;
    double *coeff, *resid;
    double *wa, *wb;
    double *qn = NULL;
    double *sol, *dsol;
    gretl_matrix *tnmat = NULL;
    gretl_matrix *ci = NULL;
    size_t rsize, isize;
    double tol = calc_eps23;
    double big = NADBL;
    double cut = 0.0;
    integer nsol = 2;  /* set to min. allowed */
    integer ndsol = 2; /* set to min. allowed */
    integer ift, lsol = 0;
    int do_ci = (opt & OPT_I)? 1 : 0;
    int iid = !(opt & OPT_R);
    int i, err = 0;

#if QDEBUG
    fprintf(stderr, "rq_fit_br: alpha = %g, do_ci = %d, iid = %d, tcrit = %d\n", 
	    alpha, do_ci, iid, (opt & OPT_N)? 0 : 1);
#endif

    rsize = p + n + n5 * p4 + n + p + nsol * p3 + ndsol * n;
    isize = n + nsol * p;

    rspace = malloc(rsize * sizeof *rspace);
    if (rspace == NULL) {
	return E_ALLOC;
    }

    ispace = malloc(isize * sizeof *ispace);
    if (ispace == NULL) {
	free(rspace);
	return E_ALLOC;
    }

    coeff = rspace;
    resid = coeff + p;
    wa = resid + n;
    wb = wa + n5 * p4;
    qn = wb + n;
    sol = qn + p;
    dsol = sol + nsol * p3;

    s = ispace;
    h = s + n;

    ci = gretl_matrix_alloc(4, p);
    tnmat = gretl_matrix_alloc(4, p);

    if (ci == NULL || tnmat == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (do_ci) {
	if (opt & OPT_N) {
	    /* no df correction */
	    cut = normal_cdf_inverse(1 - alpha/2);
	    fprintf(stderr, "Normal cutoff = %g\n", cut);
	} else {
	    cut = student_cdf_inverse(n - p, 1 - alpha/2);
	    fprintf(stderr, "t cutoff = %g\n", cut);
	}

	if (iid) {
	    /* assuming iid errors */
	    err = make_iid_qn(X, qn);
	} else {
	    /* robust variant */
	    err = make_nid_qn(y, X, qn, tau, tol);
	}
    } 
    
    if (!err) {
	logical lci1 = do_ci;

	rqbr_(&n, &p, &n5, &p3, &p4, X->val, y->val, &tau, &tol, &ift,
	      coeff, resid, s, wa, wb, &nsol, &ndsol, sol, dsol,
	      &lsol, h, qn, &cut, ci->val, tnmat->val, &big, &lci1);

	if (ift == 1) {
	    fprintf(stderr, "Warning: solution may be non-unique\n");
	} else if (ift == 2){ 
	    fprintf(stderr, "Premature end: conditioning problem in X?\n");
	    err = E_NOCONV;
	}
    }

    if (!err && do_ci) {
	rq_interpolate_intervals(ci, tnmat, cut);
    } 

    if (!err && pmod != NULL) {
	rq_transcribe_results(pmod, y, tau, coeff, resid);
	if (do_ci) {
	    rq_attach_intervals(pmod, ci, alpha, opt);
	}
    }

    if (!err && theta != NULL) {
	for (i=0; i<p; i++) {
	    theta->val[i] = coeff[i];
	}
    }    

 bailout:

    free(rspace);
    free(ispace);
    gretl_matrix_free(ci);
    gretl_matrix_free(tnmat);

    return err;
}

/* wrapper struct for use with Frisch-Newton method */

struct rq_info {
    integer n, p;
    double tau;
    double beta;
    double eps;
    double *rspace;
    double *rhs;
    double *d;
    double *u;
    double *wn;
    double *wp;
    integer nit[3];
    integer info;
};

static void rq_info_free (struct rq_info *rq)
{
    free(rq->rspace);
}

/* allocate workspace to be fed to the Fortran functions */

static int rq_info_alloc (struct rq_info *rq, int n, int p,
			  double tau)
{
    int n10 = n * 10;
    int pp4 = p * (p + 4);
    size_t rsize = p + n + n + n10 + pp4;

    rq->rspace = malloc(rsize * sizeof *rq->rspace);

    if (rq->rspace == NULL) {
	return E_ALLOC;
    }

    rq->rhs = rq->rspace;
    rq->d   = rq->rhs + p;
    rq->u   = rq->d + n;
    rq->wn  = rq->u + n;
    rq->wp  = rq->wn + n10;

    rq->n = n;
    rq->p = p;
    rq->tau = tau;
    rq->beta = .99995;
    rq->eps = 1.0e-7;

    return 0;
}

/* Initialize the tau-dependent arrays rhs and wn (rhs is
   also X-depdendent) and set other entries to 0 or 1.
*/

static void rq_workspace_init (const gretl_matrix *XT,
			       double tau, 
			       struct rq_info *ws)
{
    double xsum;
    int p = gretl_matrix_rows(XT);
    int n = gretl_matrix_cols(XT);
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

/* IID version of asymptotic F-N covariance matrix */

static int rq_fn_iid_VCV (MODEL *pmod, gretl_matrix *y,
			  gretl_matrix *XT, double tau,
			  struct rq_info *rq)
{
    gretl_matrix *V = NULL;
    gretl_matrix *vx = NULL;
    gretl_matrix *vy = NULL;
    gretl_matrix *S0 = NULL;
    gretl_matrix *S1 = NULL;
    double h, sparsity, eps23;
    int i, df, pz = 0;
    integer vn, vp = 2;
    int err = 0;

    V = get_XTX_inverse(XT, &err);
    if (err) {
	return err;
    }

    h = ceil(rq->n * rq_bandwidth(tau, rq->n));
    h = (rq->p + 1 > h)? rq->p + 1 : h;
    vn = h + 1;
    eps23 = calc_eps23;

    for (i=0; i<rq->n; i++) {
	if (fabs(rq->wn[i]) < eps23) 
	    pz++;
    }

#if QDEBUG
    fprintf(stderr, "rq_fn_iid_VCV: eps = %g, h = %g, pz = %d\n", eps23, h, pz);
#endif

    vy = gretl_column_vector_alloc(vn);
    vx = gretl_matrix_alloc(2, vn);
    S0 = gretl_matrix_alloc(rq->n, 2);

    if (vy == NULL || vx == NULL || S0 == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* construct sequence of length h + 1 */

    for (i=0; i<vn; i++) {
	vy->val[i] = i + pz + 1;
    }

    /* construct artificial X (transpose) matrix, constant in
       first column
    */

    df = rq->n - rq->p;
    for (i=0; i<vn; i++) {
	gretl_matrix_set(vx, 0, i, 1.0);
	gretl_matrix_set(vx, 1, i, vy->val[i] / df);
    }

    /* sort the residuals by absolute magnitude */

    for (i=0; i<rq->n; i++) {
	gretl_matrix_set(S0, i, 0, rq->wn[i]);
	gretl_matrix_set(S0, i, 1, fabs(rq->wn[i]));
    }

    S1 = gretl_matrix_sort_by_column(S0, 1, &err);
    if (err) {
	goto bailout;
    }

    /* extract the desired range of residuals */

    for (i=0; i<vn; i++) {
	vy->val[i] = gretl_matrix_get(S1, pz + i, 0);
    }

    /* and re-sort, respecting sign */

    qsort(vy->val, vn, sizeof *vy->val, gretl_compare_doubles);

    /* run artificial L1 regression to get sparsity measure */

    rq_workspace_init(vx, 0.5, rq);
    rqfnb_(&vn, &vp, vx->val, vy->val, rq->rhs, 
	   rq->d, rq->u, &rq->beta, &rq->eps, 
	   rq->wn, rq->wp, rq->nit, &rq->info);

    if (rq->info != 0) {
	fprintf(stderr, "rq_fn_iid_VCV: rqfnb: info = %d\n", rq->info);
	err = E_DATA;
    } else {
	/* scale X'X-inverse appropriately */
	sparsity = rq->wp[1];
	h = sparsity * sparsity * tau * (1 - tau);
	gretl_matrix_multiply_by_scalar(V, h); 
	err = gretl_model_write_vcv(pmod, V);
    }

 bailout:

    gretl_matrix_free(V);
    gretl_matrix_free(vy);
    gretl_matrix_free(vx);
    gretl_matrix_free(S0);
    gretl_matrix_free(S1);

    return err;
}

/* Robust sandwich version of F-N covariance matrix */

static int rq_fn_nid_VCV (MODEL *pmod, gretl_matrix *y,
			  gretl_matrix *XT, double tau,
			  struct rq_info *rq)
{
    gretl_matrix *p1 = NULL;
    gretl_matrix *f = NULL;
    gretl_matrix *fX = NULL;
    gretl_matrix *fXX = NULL;
    gretl_matrix *XTX = NULL;
    gretl_matrix *V = NULL;
    double h, tau_h, x, eps23;
    int n = rq->n;
    int p = rq->p;
    int i, t, err = 0;

    h = rq_bandwidth(tau, rq->n);

    if (tau + h > 1) {
	fprintf(stderr, "rq_fn_nid_VCV: tau + h > 1\n");
	return E_DATA;
    } else if (tau - h < 0) {
	fprintf(stderr, "rq_fn_nid_VCV: tau - h < 0\n");
	return E_DATA;
    }    

    eps23 = calc_eps23;

    p1  = gretl_matrix_alloc(p, 1);
    f   = gretl_matrix_alloc(n, 1);
    fX  = gretl_matrix_alloc(p, n);
    fXX = gretl_matrix_alloc(p, p);
    XTX = gretl_matrix_alloc(p, p);
    V   = gretl_matrix_alloc(p, p);

    if (p1 == NULL || f == NULL || fX == NULL ||
	fXX == NULL || XTX == NULL || V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
	
    tau_h = tau + h;
    rq_workspace_init(XT, tau_h, rq);
    rqfnb_(&n, &p, XT->val, y->val, rq->rhs, 
	   rq->d, rq->u, &rq->beta, &rq->eps, 
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
    rq_workspace_init(XT, tau_h, rq);
    rqfnb_(&n, &p, XT->val, y->val, rq->rhs, 
	   rq->d, rq->u, &rq->beta, &rq->eps, 
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

#if QDEBUG
    fprintf(stderr, "rq_fn_nid_VCV: bandwidth = %g\n", h);
    gretl_matrix_print(p1, "coeff diffs");
#endif

    gretl_matrix_multiply_mod(XT, GRETL_MOD_TRANSPOSE,
			      p1, GRETL_MOD_NONE,
			      f, GRETL_MOD_NONE);

    for (t=0; t<n; t++) {
	if (f->val[t] <= 0) {
	    fprintf(stderr, "f[%d] is non-positive\n", t);
	}
	x = (2 * h) / (f->val[t] - eps23);
	f->val[t] = (x > 0)? sqrt(x) : 0;
	f->val[t] = (x > 0)? x : 0; 
	for (i=0; i<p; i++) {
	    x = gretl_matrix_get(XT, i, t);
	    x *= f->val[t];
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
    gretl_model_set_int(pmod, "rq_nid", 1);

 bailout:

    gretl_matrix_free(p1);
    gretl_matrix_free(f);
    gretl_matrix_free(fX);
    gretl_matrix_free(fXX);
    gretl_matrix_free(XTX);
    gretl_matrix_free(V);

    return err;
}

/* sub-driver for Frisch-Newton interior point variant */

static int rq_fit_fn (gretl_matrix *y, gretl_matrix *XT, 
		      double tau, gretlopt opt, MODEL *pmod)
{
    struct rq_info rq;
    integer n = y->rows;
    integer p = XT->rows;
    int err = 0;

    err = rq_info_alloc(&rq, n, p, tau);
    if (err) {
	return err;
    }

    rq_workspace_init(XT, tau, &rq);

    /* get coefficients (in wp) and residuals (in wn) */
    rqfnb_(&n, &p, XT->val, y->val, rq.rhs, rq.d, rq.u, &rq.beta, &rq.eps, 
	   rq.wn, rq.wp, rq.nit, &rq.info);

    if (rq.info != 0) {
	fprintf(stderr, "rqfn gave info = %d\n", rq.info);
	err = E_DATA;
    }

    if (!err) {
	/* save coeffs, residuals, etc., before computing VCV */
	rq_transcribe_results(pmod, y, tau, rq.wp, rq.wn);
    }

    if (!err) {
	/* compute and add covariance matrix */
	if (opt & OPT_R) {
	    err = rq_fn_nid_VCV(pmod, y, XT, tau, &rq);
	} else {
	    err = rq_fn_iid_VCV(pmod, y, XT, tau, &rq);
	}
    }

    rq_info_free(&rq);

    return err;
}

/* Write y and X from pmod into gretl matrices, respecting the sample
   range.  Note that for use with the Frisch-Newton version of the
   underlying rq function we want the X matrix in transpose form.
*/

static int rq_make_matrices (MODEL *pmod,
			     double **Z, DATAINFO *pdinfo,
			     gretl_matrix **py,
			     gretl_matrix **pX,
			     gretlopt opt)
{
    int n = pmod->nobs;
    int p = pmod->ncoeff;
    int yno = pmod->list[1];
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    int tr = !(opt & OPT_I);
    int i, s, t, v;
    int err = 0;

    y = gretl_matrix_alloc(n, 1);

    if (tr) {
	X = gretl_matrix_alloc(p, n);
    } else {
	X = gretl_matrix_alloc(n, p);
    }

    if (X == NULL || y == NULL) {
	gretl_matrix_free(y);
	gretl_matrix_free(X);
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
	    if (tr) {
		gretl_matrix_set(X, i, s++, Z[v][t]);
	    } else {
		gretl_matrix_set(X, s++, i, Z[v][t]);
	    }
	}
    }

    if (err) {
	gretl_matrix_free(y);
	gretl_matrix_free(X);
    } else {
	*py = y;
	*pX = X;
    }

    return err;
}

static double get_user_tau (const char *s, double **Z, DATAINFO *pdinfo,
			    int *err)
{
    char *test;
    double tau = NADBL;

    *err = E_DATA;
    errno = 0;

    /* try for a numerical value */
    tau = strtod(s, &test);
    if (!errno && *test == '\0') {
	/* fine, got one */
	*err = 0;
    } else {
	/* try for named scalar */
	int v = varindex(pdinfo, s);

	if (v < pdinfo->v && var_is_scalar(pdinfo, v)) {
	    tau = Z[v][0];
	    *err = 0;
	}
	errno = 0;
    }

    if (!*err && (tau < .01 || tau > .99)) {
	gretl_errmsg_sprintf("quantreg: tau must be >= .01 and <= .99");
	*err = E_DATA;
    }	

    return tau;
}

static int get_ci_alpha (double *a)
{
    double c = get_optval_double(QUANTREG, OPT_I);

    if (na(c)) {
	*a = 0.1;
    } else if (c > 1.0 && c < 100.0) {
	*a = 1 - c / 100;
    } else if (c < 0.1 || c > .999) {
	gretl_errmsg_sprintf("Confidence level out of bounds");
	return E_DATA;
    } else {
	*a = 1 - c;
    }

    return 0;
}

int rq_driver (const char *parm, MODEL *pmod,
	       double **Z, DATAINFO *pdinfo,
	       gretlopt opt, PRN *prn)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    double tau;
    int err = 0;

    /* FIXME: we should probably accept a matrix for tau.  Also, when
       we're doing confidence intervals we should accept an alpha
       argument (right now we're hard-wired to 90%).
    */

    tau = get_user_tau(parm, Z, pdinfo, &err);
    if (err) {
	return E_DATA;
    }

    if ((opt & OPT_I) && pmod->list[0] < 3) {
	gretl_errmsg_set("quantreg: can't do confidence intervals with "
			 "only one regressor");
	err = E_DATA;
    }

    if (!err) {
	err = rq_make_matrices(pmod, Z, pdinfo, &y, &X, opt);
    }

    if (!err) {
	if (opt & OPT_I) {
	    /* doing confidence intervals -> use Borrodale-Roberts */
	    double alpha;

	    err = get_ci_alpha(&alpha);
	    if (!err) {
		err = rq_fit_br(y, X, tau, alpha, opt, pmod, NULL);
	    }
	} else {
	    /* use Frisch-Newton */
	    err = rq_fit_fn(y, X, tau, opt, pmod);
	}
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
