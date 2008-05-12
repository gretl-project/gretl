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
#include "usermat.h"

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

#define calc_eps23 (pow(2.22045e-16, 2/3.0))

/* wrapper struct for use with Barrodale-Roberts */

struct br_info {
    integer n, p;
    integer n5, p3, p4;
    integer nsol, ndsol;
    double tau;
    double tol;
    double big;
    double cut;
    double *rspace;
    double *coeff, *resid;
    double *wa, *wb;
    double *qn;
    double *sol, *dsol;
    integer *ispace;
    integer *s, *h;
    gretl_matrix *ci;
    gretl_matrix *tnmat;
};

/* wrapper struct for use with Frisch-Newton method */

struct fn_info {
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

/* destructors for the above wrappers */

static void fn_info_free (struct fn_info *rq)
{
    free(rq->rspace);
}

static void br_info_free (struct br_info *rq)
{
    free(rq->rspace);
    free(rq->ispace);
    gretl_matrix_free(rq->ci);
    gretl_matrix_free(rq->tnmat);
}

static int real_br_calc (gretl_matrix *y, gretl_matrix *X,
			 double tau, struct br_info *rq,
			 int calc_ci);

/* Bandwidth selection for sparsity estimation, as per Hall, P. and
   Sheather, S. J. (1988), "On the distribution of a studentized
   quantile", Journal of the Royal Statistical Society, Series B, 50,
   381–391: rate = O(n^{-1/3})
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

/* attach the special matrix generated when we estimate
   the model for several values of tau */

static int rq_attach_multi_intervals (MODEL *pmod, 
				      gretl_vector *tauvec,
				      gretl_matrix *tbeta,
				      double alpha, gretlopt opt)
{
#if QDEBUG
    gretl_matrix_print(tauvec, "tauvec");
    gretl_matrix_print(tbeta, "tbeta");
#endif

    gretl_model_set_matrix_as_data(pmod, "rq_tauvec", tauvec);
    gretl_model_set_matrix_as_data(pmod, "rq_sequence", tbeta);
    gretl_model_set_double(pmod, "rq_alpha", alpha);

    if (opt & OPT_R) {
	gretl_model_set_int(pmod, "rq_nid", 1);
    }

    pmod->ci = LAD;
    gretl_model_set_int(pmod, "rq", 1);

    pmod->lnL = NADBL;
    pmod->ess = NADBL;
    pmod->rsq = NADBL;
    pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;

    pmod->criterion[C_AIC] = NADBL;
    pmod->criterion[C_BIC] = NADBL;
    pmod->criterion[C_HQC] = NADBL;

    return 0;
}

static void rq_interpolate_intervals (struct br_info *rq)
{
    double c1j, c2j, c3j, c4j;
    double tn1, tn2, tn3, tn4;
    int p = gretl_matrix_cols(rq->ci);
    int j;

    for (j=0; j<p; j++) {
	c1j = gretl_matrix_get(rq->ci, 0, j);
	c2j = gretl_matrix_get(rq->ci, 1, j);
	c3j = gretl_matrix_get(rq->ci, 2, j);
	c4j = gretl_matrix_get(rq->ci, 3, j);
	tn1 = gretl_matrix_get(rq->tnmat, 0, j);
	tn2 = gretl_matrix_get(rq->tnmat, 1, j);
	tn3 = gretl_matrix_get(rq->tnmat, 2, j);
	tn4 = gretl_matrix_get(rq->tnmat, 3, j);

	c3j += fabs(c4j - c3j) * (rq->cut - fabs(tn3)) / fabs(tn4 - tn3);
	c2j -= fabs(c1j - c2j) * (rq->cut - fabs(tn2)) / fabs(tn1 - tn2);

	/* Write the (1 - alpha) intervals into rows 1 and 2 
	   of the matrix ci */

	gretl_matrix_set(rq->ci, 2, j, c3j);
	gretl_matrix_set(rq->ci, 1, j, c2j);
    }
}

/* prep for robust version of rank-inversion confidence interval 
   calculation
*/

static int make_nid_qn (gretl_matrix *y, gretl_matrix *X, 
			struct br_info *rq)
{
    int n = rq->n;
    int p = rq->p;
    double eps = rq->tol;
    gretl_matrix *b = NULL;
    gretl_matrix *f = NULL;
    gretl_matrix *xj = NULL;
    gretl_matrix *uj = NULL;
    gretl_matrix *Xcmp = NULL;
    double h, fi, ui, xik;
    int i, j, k, jj;
    int badf = 0;
    int err = 0;

    h = rq_bandwidth(rq->tau, n);

    b    = gretl_matrix_alloc(p, 1);
    f    = gretl_matrix_alloc(n, 1);
    xj   = gretl_matrix_alloc(n, 1);
    uj   = gretl_matrix_alloc(n, 1);
    Xcmp = gretl_matrix_alloc(n, p - 1);

    if (b == NULL || f == NULL || xj == NULL || 
	uj == NULL || Xcmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    real_br_calc(y, X, rq->tau + h, rq, 0);
    for (j=0; j<p; j++) {
	b->val[j] = rq->coeff[j];
    }

    real_br_calc(y, X, rq->tau - h, rq, 0);
    for (j=0; j<p; j++) {
	b->val[j] -= rq->coeff[j];
    }

#if QDEBUG
    fprintf(stderr, "make_nid_qn: bandwidth = %g\n", h);
    gretl_matrix_print(b, "coeff diffs");
#endif

    gretl_matrix_multiply(X, b, f);

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
	f->val[i] = sqrt(fi);
    }

    /* Now set each qn[j] to the SSR from an f-weighted regression of 
       X_j on the other regressors.
    */

    gretl_matrix_reuse(b, p - 1, 1);

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

	err = gretl_matrix_ols(xj, Xcmp, b, NULL, uj, NULL);

	if (!err) {
	    rq->qn[j] = 0.0;
	    for (i=0; i<n; i++) {
		ui = uj->val[i] / f->val[i];
		rq->qn[j] += ui * ui;
	    }
	}
    }

 bailout:

    gretl_matrix_free(b);
    gretl_matrix_free(f);
    gretl_matrix_free(xj);
    gretl_matrix_free(uj);
    gretl_matrix_free(Xcmp);

    return err;
}

/* maybe this should be in gretl_matrix.c ? */

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

/* prep for IID version of rank-inversion confidence interval 
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

/* allocate workspace to be fed to the function rqbr */

static int br_info_alloc (struct br_info *rq, int n, int p,
			  double tau, double alpha,
			  gretlopt opt)
{
    size_t rsize, isize;

    rq->rspace = NULL;
    rq->ispace = NULL;
    rq->ci = NULL;
    rq->tnmat = NULL;

    rq->n5 = n + 5;
    rq->p3 = p + 3;
    rq->p4 = p + 4;
    rq->nsol = rq->ndsol = 2;

    rsize = p + n + rq->n5 * rq->p4 + n + p + rq->nsol * rq->p3 + rq->ndsol * n;
    isize = n + rq->nsol * p;

    rq->rspace = malloc(rsize * sizeof *rq->rspace);
    if (rq->rspace == NULL) {
	return E_ALLOC;
    }

    rq->ispace = malloc(isize * sizeof *rq->ispace);
    if (rq->ispace == NULL) {
	return E_ALLOC;
    }

    rq->ci = gretl_matrix_alloc(4, p);
    rq->tnmat = gretl_matrix_alloc(4, p);

    if (rq->ci == NULL || rq->tnmat == NULL) {
	return E_ALLOC;
    }

    /* real arrays */
    rq->coeff = rq->rspace;
    rq->resid = rq->coeff + p;
    rq->wa    = rq->resid + n;
    rq->wb    = rq->wa + rq->n5 * rq->p4;
    rq->qn    = rq->wb + n;
    rq->sol   = rq->qn + p;
    rq->dsol  = rq->sol + rq->nsol * rq->p3;

    /* integer arrays */
    rq->s = rq->ispace;
    rq->h = rq->s + n;

    rq->n = n;
    rq->p = p;
    rq->tau = tau;
    rq->tol = calc_eps23;
    rq->big = NADBL;

    if (opt & OPT_N) {
	/* no df correction */
	rq->cut = normal_cdf_inverse(1 - alpha/2);
    } else {
	rq->cut = student_cdf_inverse(n - p, 1 - alpha/2);
    }

    return 0;
}

/* Fortran rqbr changes some of its arguments: we need to
   reset these for the next round */

static void br_info_reinit (gretl_matrix *X, struct br_info *rq)
{
    rq->n = X->rows;
    rq->p = X->cols;

    rq->n5 = rq->n + 5;
    rq->p3 = rq->p + 3;
    rq->p4 = rq->p + 4;
    rq->nsol = rq->ndsol = 2;
}

/* Call Barrodale-Roberts Fortran code */

static int real_br_calc (gretl_matrix *y, gretl_matrix *X, 
			 double tau, struct br_info *rq,
			 int calc_ci)
{
    integer ift, lsol = 0;
    logical lci1 = calc_ci;
    int err = 0;

    rqbr_(&rq->n, &rq->p, &rq->n5, &rq->p3, &rq->p4, 
	  X->val, y->val, &tau, &rq->tol, &ift,
	  rq->coeff, rq->resid, rq->s, rq->wa, rq->wb, 
	  &rq->nsol, &rq->ndsol, rq->sol, rq->dsol,
	  &lsol, rq->h, rq->qn, &rq->cut, 
	  rq->ci->val, rq->tnmat->val, 
	  &rq->big, &lci1);

    if (ift == 1) {
	fprintf(stderr, "Warning: solution may be non-unique\n");
    } else if (ift == 2){ 
	fprintf(stderr, "Premature end: conditioning problem in X?\n");
	err = E_NOCONV;
    }

    return err;
}

/* allocate workspace to be fed to the function rqfnb */

static int fn_info_alloc (struct fn_info *rq, int n, int p,
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
   also X-dependent) and set other entries to 0 or 1.
*/

static void rq_workspace_init (const gretl_matrix *XT,
			       double tau, 
			       struct fn_info *ws)
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
			  struct fn_info *rq)
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
			  struct fn_info *rq)
{
    gretl_matrix *p1 = NULL;
    gretl_matrix *f = NULL;
    gretl_matrix *fX = NULL;
    gretl_matrix *fXX = NULL;
    gretl_matrix *XTX = NULL;
    gretl_matrix *V = NULL;
    double h, x, eps23;
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
	
    rq_workspace_init(XT, tau + h, rq);
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

    rq_workspace_init(XT, tau - h, rq);
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

static int 
write_tbeta_block (gretl_matrix *tbeta, int nt, double *coeff,
		   gretl_matrix *ci, int k)

{
    double blo, bhi;
    int p = ci->cols;
    int i;

    for (i=0; i<p; i++) {
	blo = gretl_matrix_get(ci, 1, i);
	bhi = gretl_matrix_get(ci, 2, i);
	gretl_matrix_set(tbeta, k, 0, coeff[i]);
	gretl_matrix_set(tbeta, k, 1, blo);
	gretl_matrix_set(tbeta, k, 2, bhi);
	k += nt;
    }

    return 0;
}

/* sub-driver for Barrodale-Roberts estimation, with confidence
   intervals
*/

static int rq_fit_br (gretl_matrix *y, gretl_matrix *X, 
		      gretl_vector *tauvec, gretlopt opt, 
		      MODEL *pmod)
{
    struct br_info rq;
    gretl_matrix *tbeta = NULL;
    integer n = y->rows;
    integer p = X->cols;
    double tau, alpha;
    int i, ntau;
    int err = 0;

    err = get_ci_alpha(&alpha);
    if (err) {
	/* return right away: don't call br_info_free! */
	return err;
    }

    ntau = gretl_vector_get_length(tauvec);
    tau = gretl_vector_get(tauvec, 0);

    err = br_info_alloc(&rq, n, p, tau, alpha, opt);

    if (ntau > 1) {
	tbeta = gretl_zero_matrix_new(p * ntau, 3);
	if (tbeta == NULL) {
	    err = E_ALLOC;
	} 
    } 

    for (i=0; i<ntau && !err; i++) {
	tau = rq.tau = gretl_vector_get(tauvec, i);

#if QDEBUG
	fprintf(stderr, "rq_fit_br: i = %d, tau = %g\n", i, tau);
#endif

	/* preliminary calculations relating to confidence intervals */
	if (opt & OPT_R) {
	    /* robust variant */
	    err = make_nid_qn(y, X, &rq);
	} else {
	    /* assuming iid errors */
	    err = make_iid_qn(X, rq.qn);
	}

	if (!err) {
	    /* get the actual estimates */
	    err = real_br_calc(y, X, tau, &rq, 1);
	}

	if (!err) {
	    /* post-process confidence intervals */
	    rq_interpolate_intervals(&rq);
	    if (ntau == 1) {
		/* done: just transcribe everything */
		rq_transcribe_results(pmod, y, tau, rq.coeff, rq.resid);
		rq_attach_intervals(pmod, rq.ci, alpha, opt);
	    } else {
		/* using multiple tau values */
		write_tbeta_block(tbeta, ntau, rq.coeff, rq.ci, i);
		br_info_reinit(X, &rq);
	    }
	}
    }

    if (tbeta != NULL) {
	/* multiple tau values */
	if (err) {
	    gretl_matrix_free(tbeta);
	} else {
	    err = rq_attach_multi_intervals(pmod, tauvec, tbeta, alpha, opt);
	    tauvec = NULL;
	}
    }

    gretl_vector_free(tauvec);
    br_info_free(&rq);

    return err;
}

/* sub-driver for Frisch-Newton interior point variant */

static int rq_fit_fn (gretl_matrix *y, gretl_matrix *XT, 
		      double tau, gretlopt opt, MODEL *pmod)
{
    struct fn_info rq;
    integer n = y->rows;
    integer p = XT->rows;
    int err = 0;

    err = fn_info_alloc(&rq, n, p, tau);
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

    fn_info_free(&rq);

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
	/* try for a named scalar */
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

static gretl_vector *get_user_tauvec (const char *s, double ***pZ, 
				      DATAINFO *pdinfo, int *err)
{
    gretl_vector *tau = NULL;
    int i, n;

    tau = generate_matrix(s, pZ, pdinfo, err);

    if (*err) {
	return NULL;
    }

    n = gretl_vector_get_length(tau);

    if (n == 0) {
	*err = E_DATA;
    } else {
	double p;

	for (i=0; i<n; i++) {
	    p = gretl_vector_get(tau, i);
	    if (p < .01 || p > .99) {
		gretl_errmsg_sprintf("quantreg: tau must be >= .01 and <= .99");
		*err = E_DATA;
	    }
	}
    }

    if (*err) {
	gretl_matrix_free(tau);
	tau = NULL;
    }

    return tau;
}

int rq_driver (const char *parm, MODEL *pmod,
	       double ***pZ, DATAINFO *pdinfo,
	       gretlopt opt, PRN *prn)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_vector *tauvec = NULL;
    double tau = NADBL;
    int err = 0;

    if ((opt & OPT_I) && pmod->list[0] < 3) {
	gretl_errmsg_set("quantreg: can't do confidence intervals with "
			 "only one regressor");
	err = E_DATA;
    }

    if (opt & OPT_I) {
	tauvec = get_user_tauvec(parm, pZ, pdinfo, &err);
    } else {
	tau = get_user_tau(parm, *pZ, pdinfo, &err);
    }

    if (!err) {
	err = rq_make_matrices(pmod, *pZ, pdinfo, &y, &X, opt);
    }

    if (!err) {
	if (opt & OPT_I) {
	    /* doing confidence intervals -> use Borrodale-Roberts */
	    err = rq_fit_br(y, X, tauvec, opt, pmod);
	} else {
	    /* use Frisch-Newton */
	    err = rq_fit_fn(y, X, tau, opt, pmod);
	}
    }

    if (!err) {
	gretl_model_add_y_median(pmod, (*pZ)[pmod->list[1]]);
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);

    if (err && pmod->errcode == 0) {
	pmod->errcode = err;
    }

    return err;
}
