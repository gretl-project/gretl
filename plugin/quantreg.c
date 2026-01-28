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
#include "version.h"
#include "gretl_f2c.h"
#include "usermat.h"
#include "matrix_extra.h"
#include "libset.h"

#include <errno.h>

#define QDEBUG 0

/* Frisch-Newton algorithm: we use this if we're not computing
   rank-inversion confidence intervals.
*/

extern int rqfnb_ (integer *n,    /* number of observations */
		   integer *p,    /* number of parameters */
		   double *a,     /* transposed matrix of regressors */
		   double *y,     /* dependent variable vector */
		   double *rhs,
		   double *d,
		   double *u,
		   double *beta,
		   double *eps,    /* tolerance */
		   double *wn,     /* work array, length n */
		   double *wp,     /* work array, length p */
		   integer *nit,   /* iteration counts */
		   integer *info, /* exit status */
		   void (*callback)(void));

/* Modified simplex, a la Barrodale-Roberts: this variant lets us get
   confidence intervals using the rank-inversion method.
*/

extern int rqbr_ (int n,         /* number of observations */
		  int p,         /* number of parameters */
		  double *x,     /* matrix of regressors */
		  double *y,     /* dependent var vector */
		  double tau,    /* the desired quantile */
		  double tol,    /* machine precision to the 2/3 */
		  double *coeff, /* p-vector of parameter estimates */
		  double *resid, /* n-vector of residuals */
		  int *s,        /* int work array, length n */
		  double *wa,    /* real work array, size (n+5) * (p+4) */
		  double *wb,    /* and another, size n */
		  double *sol,   /* primal solution array, size (p + 3) * 2 */
		  double *dsol,  /* dual solution array, size n * 2 */
		  int *h,        /* matrix of observation indices, size p * nsol */
		  double *qn,    /* vector of residual variances from projection of
				    each column of x on the remaining columns */
		  double cutoff, /* critical point for interval */
		  double *ci,    /* matrix of confidence intervals, size 4 * p */
		  double *tnmat, /* matrix of JGPK rank-test statistics */
		  double big,    /* "large positive floating-point number" */
		  int rmax,      /* max iterations (added for gretl) */
		  int ci1,       /* doing confidence intervals? */
		  void (*callback)(void));

/* machine precision to the 2/3 */
#define calc_eps23 (pow(2.22045e-16, 2/3.0))

/* wrapper struct for use with Barrodale-Roberts */

struct br_info {
    int warning;
    int rmax;
    int n, p;
    int n5, p3, p4;
    int nsol, ndsol;
    double tau;
    double tol;
    double big;
    double cut;
    double *rspace;
    double *coeff, *resid;
    double *wa, *wb;
    double *qn;
    double *sol, *dsol;
    int *ispace;
    int *s, *h;
    gretl_matrix *ci;
    gretl_matrix *tnmat;
    void (*callback)();
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
    double *resid;
    double *coeff;
    integer nit[3];
    integer info;
    void (*callback)();
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

static double hs_bandwidth (double tau, int n, int *err)
{
    double alpha = 0.05;
    double x0 = normal_cdf_inverse(tau);
    double f0 = normal_pdf(x0);
    double b1 = pow(n, -1 / 3.0);
    double b2 = pow(normal_cdf_inverse(1 - alpha/2), 2 / 3.0);
    double b3 = (1.5 * f0 * f0) / (2 * x0 * x0 + 1);
    double h = b1 * b2 * pow(b3, 1 / 3.0);

#if QDEBUG
    fprintf(stderr, "tau = %g, hs bandwidth = %g\n", tau, h);
#endif

    if (err != NULL) {
	if (tau + h > 1) {
	    gretl_errmsg_set("Hall-Sheather bandwidth is out of bounds");
	    fprintf(stderr, "hs_bandwidth: tau + h > 1\n");
	    *err = E_DATA;
	} else if (tau - h < 0) {
	    gretl_errmsg_set("Hall-Sheather bandwidth is out of bounds");
	    fprintf(stderr, "hs_bandwidth: tau - h < 0\n");
	    *err = E_DATA;
	}
    }

    return h;
}

/* Compute the loglikelihood for a quantile model */

static double rq_loglik (MODEL *pmod, double tau)
{
    double R = 0.0;
    int n = pmod->nobs;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    R += pmod->uhat[t] * (tau - (pmod->uhat[t] < 0));
	}
    }

    return n * (log(tau * (1-tau)) - 1 - log(R/n));
}

enum {
    RQ_STAGE_1,
    RQ_STAGE_2,
    RQ_LAD
};

/* Transcribe quantreg results into model struct.  Note: we have to be
   careful to replace or invalidate any values generated by the
   initial OLS run that are not appropriate in the quantreg context,
   so we don't serve up misleading values when the user calls a
   model-data accessor.
*/

static void rq_transcribe_results (MODEL *pmod,
				   const gretl_matrix *y,
				   double tau,
				   const double *b,
				   const double *u,
				   int stage)
{
    double SAR = 0.0;
    int i, s, t;

    if (stage == RQ_STAGE_1) {
	gretl_model_set_double(pmod, "tau", tau);
    }

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = b[i];
	if (stage == RQ_STAGE_1 || stage == RQ_LAD) {
	    pmod->sderr[i] = NADBL; /* these will be replaced later */
	}
    }

    pmod->ess = 0.0;

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    pmod->uhat[t] = u[s];
	    pmod->yhat[t] = y->val[s] - u[s];
	    SAR += fabs(u[s]);
	    pmod->ess += u[s] * u[s];
	    s++;
	}
    }

    /* sum of absolute residuals */
    gretl_model_set_double(pmod, "ladsum", SAR);

    /* invalidate unused stats */
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = pmod->chisq = NADBL;

    /* LaPlace errors: equivalent of standard error is sum of
       absolute residuals over nobs */
    pmod->sigma = SAR / pmod->nobs;
    pmod->lnL = rq_loglik(pmod, tau);
    mle_criteria(pmod, 0);
}

/* Extract the interpolated lower and upper bounds from the matrix ci
   into a new matrix and attach this to the model for printing.
*/

static int rq_attach_intervals (MODEL *pmod, struct br_info *rq,
				double alpha, gretlopt opt)
{
    gretl_matrix *ci = rq->ci;
    gretl_matrix *rqci;
    double blo, bhi;
    int i, p = ci->cols;
    int err = 0;

    rqci = gretl_matrix_alloc(p, 2);
    if (rqci == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<p; i++) {
	blo = gretl_matrix_get(ci, 1, i);
	bhi = gretl_matrix_get(ci, 2, i);
	if (na(blo) || na(bhi) || (blo == 0.0 && bhi == 0.0)) {
	    err = E_NOCONV;
	    break;
	}
	gretl_matrix_set(rqci, i, 0, blo);
	gretl_matrix_set(rqci, i, 1, bhi);
    }

    if (err) {
	gretl_matrix_free(rqci);
    } else {
	gretl_model_set_matrix_as_data(pmod, "coeff_intervals", rqci);
	gretl_model_set_double(pmod, "rq_alpha", alpha);
    }

    return err;
}

static void correct_multi_tau_model (MODEL *pmod)
{
    pmod->lnL = NADBL;
    pmod->ess = NADBL;
    pmod->sigma = NADBL;
    pmod->rsq = NADBL;
    pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;

    pmod->criterion[C_AIC] = NADBL;
    pmod->criterion[C_BIC] = NADBL;
    pmod->criterion[C_HQC] = NADBL;

    free(pmod->uhat);
    free(pmod->yhat);
    free(pmod->xpx);
    free(pmod->vcv);

    pmod->uhat = NULL;
    pmod->yhat = NULL;
    pmod->xpx = NULL;
    pmod->vcv = NULL;
}

/* Attach the special results matrix generated when we estimate the
   model for several values of tau. In this case we retain the
   OLS coefficients and standard errors for comparison, but
   otherwise we invalidate most statistics.
*/

static int rq_attach_multi_results (MODEL *pmod,
				    const gretl_vector *tauvec,
				    gretl_matrix *tbeta,
				    double alpha, gretlopt opt)
{
    gretl_vector *tcpy;

#if QDEBUG
    gretl_matrix_print(tauvec, "tauvec");
    gretl_matrix_print(tbeta, "tbeta");
#endif

    tcpy = gretl_matrix_copy(tauvec);
    gretl_model_set_matrix_as_data(pmod, "rq_tauvec", tcpy);
    gretl_model_set_matrix_as_data(pmod, "rq_sequence", tbeta);

    if (alpha > 0) {
	gretl_model_set_double(pmod, "rq_alpha", alpha);
    }

    correct_multi_tau_model(pmod);

    return 0;
}

static int rq_interpolate_intervals (struct br_info *rq)
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

	/* check for garbage */
	if (na(c3j) || na(c2j)) {
	    fprintf(stderr, "rq_interpolate_intervals: coeff %d: low = %g, "
		    "high = %g\n", j, c2j, c3j);
	    gretl_matrix_print(rq->ci, "rq->ci");
	    gretl_matrix_print(rq->tnmat, "rq->tnmat");
	    gretl_errmsg_set(_("Couldn't calculate confidence intervals "
			       "for this model"));
	    return E_NAN;
	}

	/* Write the (1 - alpha) intervals into rows 1 and 2
	   of the matrix ci.
	*/
	gretl_matrix_set(rq->ci, 2, j, c3j);
	gretl_matrix_set(rq->ci, 1, j, c2j);
    }

    return 0;
}

static void bad_f_count (const gretl_matrix *f)
{
    int n = gretl_vector_get_length(f);
    int i, badf = 0;

    for (i=0; i<n; i++) {
	if (f->val[i] <= 0) {
	    badf++;
	}
    }

    if (badf > 0) {
	fprintf(stderr, "Warning: %g percent of fi's <= 0\n",
		100 * (double) badf / n);
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
    int err = 0;

    h = hs_bandwidth(rq->tau, n, &err);
    if (err) {
	return err;
    }

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
    gretl_matrix_print(b, "coeff diffs");
#endif

    gretl_matrix_multiply(X, b, f);

    bad_f_count(f);

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

/* Get {X'X}^{-1}, handling both the B-R case, where the X matrix
   has the shape you'd expect, and the F-N case, where "X" is
   actually X-transpose.
*/

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

#if QDEBUG
    fprintf(stderr, "make_iid_qn: returning %d\n", err);
#endif

    return err;
}

/* Allocate workspace to be fed to the function rqbr, and
   initialize various things.
*/

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

    if (!(opt & OPT_L)) {
	/* doing confidence intervals */
	rq->ci = gretl_matrix_alloc(4, p);
	rq->tnmat = gretl_matrix_alloc(4, p);

	if (rq->ci == NULL || rq->tnmat == NULL) {
	    return E_ALLOC;
	}
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

    rq->warning = 0;
    rq->n = n;
    rq->p = p;
    rq->tau = tau;
    rq->tol = calc_eps23;
    rq->big = DBL_MAX / 100;
    rq->rmax = libset_get_int(RQ_MAXITER);

    if (opt & OPT_L) {
	/* doing simple LAD */
	rq->cut = 0;
    } else if (opt & OPT_N) {
	/* asymptotic: no df correction */
	rq->cut = normal_cdf_inverse(1 - alpha/2);
    } else {
	rq->cut = student_cdf_inverse(n - p, 1 - alpha/2);
    }

    if (show_activity_func_installed()) {
	rq->callback = show_activity_callback;
    } else {
	rq->callback = NULL;
    }

    return 0;
}

/* Call Barrodale-Roberts code */

static int real_br_calc (gretl_matrix *y, gretl_matrix *X,
			 double tau, struct br_info *rq,
			 int calc_ci)
{
    double *ci_val, *tnmat_val;
    int ret, err = 0;

#if QDEBUG
    fprintf(stderr, "real_br_calc: calling rqbr, calc_ci = %d\n", calc_ci);
#endif

    /* these inputs are not needed if we're not computing
       confidence intervals */
    ci_val = (rq->ci == NULL)? NULL : rq->ci->val;
    tnmat_val = (rq->tnmat == NULL)? NULL : rq->tnmat->val;

    ret = rqbr_(rq->n, rq->p, X->val, y->val, tau, rq->tol,
		rq->coeff, rq->resid, rq->s, rq->wa, rq->wb,
		rq->sol, rq->dsol, rq->h, rq->qn, rq->cut,
		ci_val, tnmat_val,
		rq->big, rq->rmax, calc_ci,
		rq->callback);

#if QDEBUG
    fprintf(stderr, "rqbr: ift = %d\n", ret);
#endif

    if (ret == 1) {
	rq->warning = 1;
	fprintf(stderr, "Warning: solution may be non-unique\n");
    } else if (ret == 2){
	fprintf(stderr, "Premature end: conditioning problem in X?\n");
	err = E_NOCONV;
    } else if (ret == 3) {
	gretl_errmsg_sprintf("Maximum number of iterations (%d) exceeded",
			     (int) rq->rmax);
	err = E_NOCONV;
    }

    return err;
}

/* allocate workspace for F-N algorithm */

static int fn_info_alloc (struct fn_info *rq, int n, int p,
			  double tau, gretlopt opt)
{
    int n10 = n * 10;
    int pp4, rp = p;
    size_t rsize;

    if (p == 1 && !(opt & OPT_R)) {
	/* just one regressor, doing "iid" standard errors:
	   we'll need workspace for p = 2 */
	rp = 2;
    }

    pp4 = rp * (rp + 4);
    rsize = rp + n + n + n10 + pp4;

    rq->rspace = malloc(rsize * sizeof *rq->rspace);

    if (rq->rspace == NULL) {
	return E_ALLOC;
    }

    rq->rhs = rq->rspace;
    rq->d   = rq->rhs + rp;
    rq->u   = rq->d + n;
    rq->resid = rq->u + n;
    rq->coeff = rq->resid + n10;

    rq->n = n;
    rq->p = p;
    rq->tau = tau;
    rq->beta = .99995;
    rq->eps = 1.0e-7;

    if (show_activity_func_installed()) {
	rq->callback = show_activity_callback;
    } else {
	rq->callback = NULL;
    }

    return 0;
}

/* Initialize the tau-dependent arrays rhs and resid (rhs is
   also X-dependent) and set other entries to 0 or 1.
*/

static void rq_workspace_init (struct fn_info *rq,
			       const gretl_matrix *XT,
			       double tau)
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
	rq->rhs[i] = tau * xsum;
    }

    for (i=0; i<n; i++) {
	rq->d[i] = rq->u[i] = 1.0;
	rq->resid[i] = tau;
    }

    for (i=n; i<n10; i++) {
	rq->resid[i] = 0.0;
    }
}

static int rq_call_FN (integer *n, integer *p, gretl_matrix *XT,
		       gretl_matrix *y, struct fn_info *rq,
		       double tau)
{
    rq_workspace_init(rq, XT, tau);

    return rqfnb_(n, p, XT->val, y->val, rq->rhs,
		  rq->d, rq->u, &rq->beta, &rq->eps,
		  rq->resid, rq->coeff, rq->nit, &rq->info,
		  rq->callback);
}

static int rq_write_variance (const gretl_matrix *V,
			      MODEL *pmod, double *se)
{
    double x;
    int i, err = 0;

    if (se != NULL) {
	for (i=0; i<V->cols; i++) {
	    x = gretl_matrix_get(V, i, i);
	    if (na(x) || x < 0) {
		se[i] = NADBL;
	    } else {
		se[i] = sqrt(x);
	    }
	}
    } else {
	err = gretl_model_write_vcv(pmod, V);
    }

    return err;
}

/* IID version of asymptotic F-N covariance matrix */

static int rq_fn_iid_VCV (MODEL *pmod, gretl_matrix *y,
			  gretl_matrix *XT, double tau,
			  struct fn_info *rq,
			  double *se)
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

    h = ceil(rq->n * hs_bandwidth(tau, rq->n, NULL));
    h = (rq->p + 1 > h)? rq->p + 1 : h;
    vn = h + 1;
    eps23 = calc_eps23;

    for (i=0; i<rq->n; i++) {
	if (fabs(rq->resid[i]) < eps23) {
	    pz++;
	}
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
	gretl_matrix_set(S0, i, 0, rq->resid[i]);
	gretl_matrix_set(S0, i, 1, fabs(rq->resid[i]));
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

    err = rq_call_FN(&vn, &vp, vx, vy, rq, 0.5);

    if (err) {
	fprintf(stderr, "rq_fn_iid_VCV: rqfn: info = %d\n", rq->info);
    } else {
	/* scale X'X-inverse appropriately */
	sparsity = rq->coeff[1];
	h = sparsity * sparsity * tau * (1 - tau);
	gretl_matrix_multiply_by_scalar(V, h);
	err = rq_write_variance(V, pmod, se);
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
			  struct fn_info *rq,
			  double *se)
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

    h = hs_bandwidth(tau, rq->n, &err);
    if (err) {
	return err;
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

    err = rq_call_FN(&n, &p, XT, y, rq, tau + h);
    if (err) {
	fprintf(stderr, "tau + h: info = %d\n", rq->info);
	goto bailout;
    }

    /* grab coeffs */
    for (i=0; i<p; i++) {
	p1->val[i] = rq->coeff[i];
    }

    err = rq_call_FN(&n, &p, XT, y, rq, tau - h);
    if (err) {
	fprintf(stderr, "tau - h: info = %d\n", rq->info);
	goto bailout;
    }

    /* diff coeffs */
    for (i=0; i<p; i++) {
	p1->val[i] -= rq->coeff[i];
    }

#if QDEBUG
    fprintf(stderr, "rq_fn_nid_VCV: bandwidth = %g\n", h);
    gretl_matrix_print(p1, "coeff diffs");
#endif

    gretl_matrix_multiply_mod(XT, GRETL_MOD_TRANSPOSE,
			      p1, GRETL_MOD_NONE,
			      f, GRETL_MOD_NONE);

    bad_f_count(f);

    for (t=0; t<n; t++) {
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

    err = rq_write_variance(V, pmod, se);

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
    int err = 0;
    double c;

    c = get_optval_double(QUANTREG, OPT_I, &err);
    if (err) {
	return err;
    }

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

/* write coefficients and lower/upper c.i. values for
   all parameters at a given value of tau */

static int
write_tbeta_block_br (gretl_matrix *tbeta, int nt, double *coeff,
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

/* write coefficients or standard errors for
   all parameters at a given value of tau */

static int
write_tbeta_block_fn (gretl_matrix *tbeta, int nt, double *x,
		      integer p, int k, int j)

{
    int i;

    for (i=0; i<p; i++) {
	if (na(x[i])) {
	    fprintf(stderr, "write_tbeta_block_fn: x[%d] = %g\n", i, x[i]);
	    return E_NAN;
	}
	gretl_matrix_set(tbeta, k, j, x[i]);
	k += nt;
    }

    return 0;
}

/* Sub-driver for Barrodale-Roberts estimation, with confidence
   intervals.
*/

static int rq_fit_br (gretl_matrix *y, gretl_matrix *X,
		      const gretl_vector *tauvec, gretlopt opt,
		      MODEL *pmod)
{
    struct br_info rq;
    gretl_matrix *tbeta = NULL;
    integer n = y->rows;
    integer p = X->cols;
    double tau, alpha = 0;
    int i, ntau;
    int err = 0;

    err = get_ci_alpha(&alpha);
    if (err) {
	return err;
    }

    ntau = gretl_vector_get_length(tauvec);
    tau = gretl_vector_get(tauvec, 0);

    err = br_info_alloc(&rq, n, p, tau, alpha, opt);

    if (!err && ntau > 1) {
	tbeta = gretl_zero_matrix_new(p * ntau, 3);
	if (tbeta == NULL) {
	    err = E_ALLOC;
	}
#if QDEBUG
	fprintf(stderr, "p = %d, ntau = %d, alpha = %g\n", p, ntau, alpha);
	fprintf(stderr, "tbeta = %d x %d\n", tbeta->rows, tbeta->cols);
#endif
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
	    err = rq_interpolate_intervals(&rq);
	}

	if (!err) {
	    if (ntau == 1) {
		/* done: put intervals onto the model */
		err = rq_attach_intervals(pmod, &rq, alpha, opt);
		if (!err) {
		    rq_transcribe_results(pmod, y, tau, rq.coeff,
					  rq.resid, RQ_STAGE_2);
		}
	    } else {
		/* using multiple tau values */
		err = write_tbeta_block_br(tbeta, ntau, rq.coeff, rq.ci, i);
	    }
	}
    }

    if (!err && rq.warning) {
	gretl_model_set_int(pmod, "nonunique", 1);
    }

    if (tbeta != NULL) {
	/* multiple tau values */
	if (err) {
	    gretl_matrix_free(tbeta);
	} else {
	    err = rq_attach_multi_results(pmod, tauvec, tbeta, alpha, opt);
	}
    }

    br_info_free(&rq);

    return err;
}

/* sub-driver for Frisch-Newton interior point variant */

static int rq_fit_fn (gretl_matrix *y, gretl_matrix *XT,
		      const gretl_vector *tauvec, gretlopt opt,
		      MODEL *pmod)
{
    struct fn_info rq;
    gretl_matrix *tbeta = NULL;
    double *se = NULL;
    integer n = y->rows;
    integer p = XT->rows;
    double tau;
    int i, ntau;
    int err = 0;

    ntau = gretl_vector_get_length(tauvec);
    tau = gretl_vector_get(tauvec, 0);

    err = fn_info_alloc(&rq, n, p, tau, opt);
    if (err) {
	return err;
    }

    if (ntau > 1) {
	tbeta = gretl_zero_matrix_new(p * ntau, 2);
	se = malloc(p * sizeof *se);
	if (tbeta == NULL || se == NULL) {
	    err = E_ALLOC;
	}
    }

    for (i=0; i<ntau && !err; i++) {
	tau = rq.tau = gretl_vector_get(tauvec, i);

#if QDEBUG
	fprintf(stderr, "rq_fit_fn: i = %d, tau = %g\n", i, tau);
#endif

	/* get coefficients and residuals */
	err = rq_call_FN(&n, &p, XT, y, &rq, tau);
	if (err) {
	    fprintf(stderr, "rqfn gave info = %d\n", rq.info);
	}

	if (!err) {
	    if (ntau == 1) {
		/* save coeffs, residuals, etc., before computing VCV */
		rq_transcribe_results(pmod, y, tau, rq.coeff, rq.resid,
				      RQ_STAGE_1);
	    } else {
		/* write coeffs for this tau value */
		write_tbeta_block_fn(tbeta, ntau, rq.coeff, p, i, 0);
	    }
	}

	if (!err) {
	    /* compute covariance matrix */
	    if (opt & OPT_R) {
		err = rq_fn_nid_VCV(pmod, y, XT, tau, &rq, se);
	    } else {
		err = rq_fn_iid_VCV(pmod, y, XT, tau, &rq, se);
	    }
	}

	if (!err && ntau > 1) {
	    /* write std errs for this tau */
	    write_tbeta_block_fn(tbeta, ntau, se, p, i, 1);
	}
    }

    if (tbeta != NULL) {
	/* multiple tau values */
	if (err) {
	    gretl_matrix_free(tbeta);
	} else {
	    err = rq_attach_multi_results(pmod, tauvec, tbeta, 0, opt);
	}
    }

    fn_info_free(&rq);
    free(se);

    return err;
}

/* Write y and X from pmod into gretl matrices, respecting the sample
   range.  Note that for use with the Frisch-Newton version of the
   underlying rq function we want the X matrix in transpose form,
   which is signaled by @tr = 1. Also note that @pmod was first
   estimated via OLS and its sample range may contain interior
   NAs, in which case we have to avoid them when filling the
   data matrices.
*/

static int rq_make_matrices (MODEL *pmod,
			     DATASET *dset,
			     gretl_matrix **py,
			     gretl_matrix **pX,
			     int tr)
{
    int n = pmod->nobs;
    int p = pmod->ncoeff;
    int yno = pmod->list[1];
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
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
	if (!na(pmod->uhat[t])) {
	    gretl_vector_set(y, s++, dset->Z[yno][t]);
	}
    }

    for (i=0; i<p; i++) {
	v = pmod->list[i+2];
	s = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t])) {
		if (tr) {
		    gretl_matrix_set(X, i, s++, dset->Z[v][t]);
		} else {
		    gretl_matrix_set(X, s++, i, dset->Z[v][t]);
		}
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

static int rq_transpose_X (gretl_matrix **pX)
{
    gretl_matrix *XT = gretl_matrix_copy_transpose(*pX);
    int err = 0;

    if (XT == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_free(*pX);
	*pX = XT;
    }

    return err;
}

static int check_user_tau (const gretl_matrix *tau, int *ntau)
{
    int err = 0;

    *ntau = gretl_vector_get_length(tau);

    if (*ntau == 0) {
	err = E_DATA;
    } else {
	double p;
	int i;

	for (i=0; i<*ntau; i++) {
	    p = gretl_vector_get(tau, i);
	    if (p < .01 || p > .99) {
		gretl_errmsg_sprintf("quantreg: tau must be >= .01 and <= .99");
		err = E_DATA;
	    }
	}
    }

    return err;
}

int rq_driver (const gretl_matrix *tau, MODEL *pmod,
	       DATASET *dset, gretlopt opt, PRN *prn)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    int ntau = 0;
    int err;

    err = check_user_tau(tau, &ntau);

    if (!err && (opt & OPT_I) && pmod->list[0] < 3) {
	gretl_errmsg_set("quantreg: can't do confidence intervals with "
			 "only one regressor");
	err = E_DATA;
    }

    if (!err) {
	if ((opt & OPT_I) && ntau == 1) {
	    /* doing intervals, only one tau: we'll do a first run
	       using F-N so the model will be equipped with standard
	       errors
	    */
	    err = rq_make_matrices(pmod, dset, &y, &X, 1);
	    if (!err) {
		err = rq_fit_fn(y, X, tau, opt, pmod);
	    }
	    if (!err) {
		/* flip the X matrix for use with B-R */
		err = rq_transpose_X(&X);
	    }
	} else {
	    int tr = !(opt & OPT_I);

	    err = rq_make_matrices(pmod, dset, &y, &X, tr);
	}
    }

    if (!err) {
	if (opt & OPT_I) {
	    /* doing confidence intervals -> use Barrodale-Roberts */
	    err = rq_fit_br(y, X, tau, opt, pmod);
	} else {
	    /* otherwise use Frisch-Newton */
	    err = rq_fit_fn(y, X, tau, opt, pmod);
	}
    }

    if (!err) {
	/* some common finishing touches */
	gretl_model_add_y_median(pmod, dset->Z[pmod->list[1]]);
	pmod->ci = QUANTREG;
	gretl_model_set_int(pmod, "rq", 1);
	if (opt & OPT_R) {
	    gretl_model_set_int(pmod, "rq_nid", 1);
	    pmod->opt |= OPT_R;
	}
	gretl_model_set_vcv_info(pmod, VCV_RQ, (opt & OPT_R)?
				 RQ_NID : RQ_ASY);
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);

    if (err && pmod->errcode == 0) {
	pmod->errcode = err;
    }

    return err;
}

/* restock the y and X matrices with a bootstrap sample */

static void rq_refill_matrices (MODEL *pmod,
				DATASET *dset,
				gretl_matrix *y,
				gretl_matrix *X,
				int *sample)
{
    int n = pmod->nobs;
    int p = pmod->ncoeff;
    int yno = pmod->list[1];
    int i, j, t, v;

    for (i=0; i<n; i++) {
	t = sample[i];
	gretl_vector_set(y, i, dset->Z[yno][t]);
    }

    for (j=0; j<p; j++) {
	v = pmod->list[j+2];
	for (i=0; i<n; i++) {
	    t = sample[i];
	    gretl_matrix_set(X, i, j, dset->Z[v][t]);
	}
    }
}

static int *good_observations_array (MODEL *pmod)
{
    int *g = malloc(pmod->nobs * sizeof *g);

    if (g != NULL) {
	int t, s = 0;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t])) {
		g[s++] = t;
	    }
	}
    }

    return g;
}

#define ITERS 500

/* obtain bootstrap estimates of LAD covariance matrix */

static int lad_bootstrap_vcv (MODEL *pmod, DATASET *dset,
			      gretl_matrix *y, gretl_matrix *X,
			      struct br_info *rq)
{
    double **coeffs = NULL;
    double *meanb = NULL;
    int *sample = NULL;
    int *goodobs = NULL;
    double xi, xj;
    int i, j, k;
    int nc = pmod->ncoeff;
    int nvcv, n = pmod->nobs;
    int err = 0;

    /* note: new_vcv sets all entries to zero */
    err = gretl_model_new_vcv(pmod, &nvcv);
    if (err) {
	return err;
    }

    /* an array for each coefficient */
    coeffs = doubles_array_new(nc, ITERS);

    /* a scalar for each coefficient mean */
    meanb = malloc(nc * sizeof *meanb);

    /* resampling array of length pmod->nobs */
    sample = malloc(n * sizeof *sample);

    if (coeffs == NULL || meanb == NULL || sample == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* apparatus for handling interior NAs */
    if (model_has_missing_obs(pmod)) {
	goodobs = good_observations_array(pmod);
	if (goodobs == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    for (k=0; k<ITERS && !err; k++) {
	/* create random sample index array */
	for (i=0; i<n; i++) {
	    j = gretl_rand_int_max(n);
	    if (goodobs != NULL) {
		sample[i] = goodobs[j];
	    } else {
		sample[i] = pmod->t1 + j;
	    }
	}

	rq_refill_matrices(pmod, dset, y, X, sample);

	/* re-estimate LAD model */
	err = real_br_calc(y, X, 0.5, rq, 0);

	if (!err) {
	    for (i=0; i<nc; i++) {
		coeffs[i][k] = rq->coeff[i];
	    }
	}
    }

    /* find means of coeff estimates */
    for (i=0; i<nc && !err; i++) {
	double bbar = 0.0;

	for (k=0; k<ITERS; k++) {
	   bbar += coeffs[i][k];
	}
	meanb[i] = bbar / ITERS;
    }

    /* find variances and covariances */
    for (i=0; i<nc && !err; i++) {
	double vi = 0.0;

	for (k=0; k<ITERS; k++) {
	    xi = coeffs[i][k] - meanb[i];
	    vi += xi * xi;
	    for (j=0; j<=i; j++) {
		xj = coeffs[j][k] - meanb[j];
		pmod->vcv[ijton(i, j, nc)] += xi * xj;
	    }
	}
	pmod->sderr[i] = sqrt(vi / ITERS);
    }

    if (!err) {
	for (i=0; i<nvcv; i++) {
	    pmod->vcv[i] /= ITERS;
	}
    }

 bailout:

    free(sample);
    free(meanb);
    doubles_array_free(coeffs, nc);

    if (goodobs != NULL) {
	free(goodobs);
    }

    return err;
}

static void lad_scrub_vcv (MODEL *pmod)
{
    int i;

    free(pmod->vcv);
    pmod->vcv = NULL;
    free(pmod->xpx);
    pmod->xpx = NULL;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->sderr[i] = NADBL;
    }
}

static int lad_fit_br (MODEL *pmod, DATASET *dset,
		       gretl_matrix *y, gretl_matrix *X,
		       gretlopt opt)
{
    struct br_info rq;
    integer n = y->rows;
    integer p = X->cols;
    int err = 0;

    err = br_info_alloc(&rq, n, p, 0.5, 0.0, OPT_L);

    if (!err) {
	/* get the actual estimates */
	err = real_br_calc(y, X, 0.5, &rq, 0);
    }

    if (!err) {
	rq_transcribe_results(pmod, y, 0.5, rq.coeff,
			      rq.resid, RQ_LAD);
	if (rq.warning) {
	    gretl_model_set_int(pmod, "nonunique", 1);
	}
    }

    if (!err) {
	if (opt & OPT_N) {
	    /* --no-vcv */
	    lad_scrub_vcv(pmod);
	} else {
	    err = lad_bootstrap_vcv(pmod, dset, y, X, &rq);
	}
    }

    br_info_free(&rq);

    return err;
}

/* Support the "lad" command: estimate with tau = 0.5, using
   the Barrodale-Roberts method, and add a bootstrapped
   covariance matrix.
*/

int lad_driver (MODEL *pmod, DATASET *dset, gretlopt opt)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    int err;

    err = rq_make_matrices(pmod, dset, &y, &X, 0);

    if (!err) {
	err = lad_fit_br(pmod, dset, y, X, opt);
    }

    if (!err) {
	gretl_model_add_y_median(pmod, dset->Z[pmod->list[1]]);
	pmod->ci = LAD;
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);

    if (err && pmod->errcode == 0) {
	pmod->errcode = err;
    }

    return err;
}
