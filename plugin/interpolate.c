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
 */

#include "libgretl.h"
#include "version.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"
#include "libset.h"

#define CL_DEBUG 0

enum {
    AGG_SUM, /* sum */
    AGG_AVG, /* average */
    AGG_EOP, /* end of period */
    AGG_SOP  /* start of period */
};

enum {
    R_ACF1,  /* "traditional" Chow-Lin */
    R_MLE,   /* GLS + MLE (Bournay-Laroque) */
    R_SSR,   /* using GLS criterion */
    R_FIXED, /* pre-specified by user */
};

struct chowlin {
    int n;
    double targ;
};

/* Callback for fzero(), as we adjust the coefficient @a so the
   theoretically derived ratio of polynomials in @a matches the
   empirical first-order autocorrelation of the OLS residuals
   (cl->targ). Return the residual.
*/

static double chow_lin_callback (double a, void *p)
{
    struct chowlin *cl = (struct chowlin *) p;
    double r, num, den, resid;

    if (a == 0) {
	r = 0;
    } else {
	/* Calculate the ratio of immediate off-diagonal
	   element of CVC' to diagonal element. Avoid use
	   of pow() since all we require are successive
	   integer powers of @a.
	*/
	double apow = a;
	int n = cl->n;
	int np = 2 * n;
	int i, coef = 1;

	num = 0.0;
	for (i=0; i<np-1; i++) {
	    num += coef * apow;
	    apow *= a;
	    coef += (i < n-1)? 1 : -1;
	}
	den = n;
	apow = a;
	for (i=1; i<n; i++) {
	    den += 2*(n-i) * apow;
	    apow *= a;
	}
	r = num/den;
    }

    resid = r - cl->targ;

#if CL_DEBUG
    fprintf(stderr, "chow_lin_callback: target %g, a %g residual %g\n",
	    cl->targ, a, resid);
#endif

    return resid;
}

typedef struct ar1data_ {
    const gretl_matrix *y;
    const gretl_matrix *X;
} ar1data;

/* BFGS callback for ar1_mle(); see Davidson and MacKinnon,
   ETM, pp. 435-6.
*/

static double ar1_loglik (const double *theta, void *data)
{
    ar1data *a = data;
    int n = a->y->rows, k = a->X->cols;
    double r = theta[0];
    double s = theta[1];
    const double *b = theta + 2;
    double onemr2 = 1.0 - r*r;
    double inv2s2 = 1.0 / (2*s*s);
    double ll1 = -0.5*n*LN_2_PI - n*log(s) + 0.5*log(onemr2);
    double u, yf2, Xb1, Xb = 0;
    int i, t;

    /* the first observation */
    for (i=0; i<k; i++) {
	Xb += gretl_matrix_get(a->X, 0, i) * b[i];
    }
    u = a->y->val[0] - Xb;
    yf2 = onemr2 * u * u;

    /* subsequent observations */
    for (t=1; t<n; t++) {
	Xb1 = Xb;
	Xb = 0;
	for (i=0; i<k; i++) {
	    Xb += gretl_matrix_get(a->X, t, i) * b[i];
	}
	u = a->y->val[t] - r*a->y->val[t-1] - Xb + r*Xb1;
	yf2 += u * u;
    }

    return ll1 - inv2s2 * yf2;
}

static int ar1_mle (const gretl_matrix *y,
		    const gretl_matrix *X,
		    const gretl_matrix *b,
		    double s, double *rho)
{
    struct ar1data_ a = {y, X};
    double *theta;
    int fc = 0, gc = 0;
    int i, nt, err;

    nt = X->cols + 2;
    theta = malloc(nt * sizeof *theta);
    if (theta == NULL) {
	return E_ALLOC;
    }

    theta[0] = *rho;
    theta[1] = s;
    for (i=0; i<X->cols; i++) {
	theta[i+2] = b->val[i];
    }

    err = BFGS_max(theta, nt, 300, 1.0e-10,
		   &fc, &gc, ar1_loglik, C_LOGLIK,
		   NULL, &a, NULL, OPT_NONE, NULL);

    if (err) {
	if (err == E_NOCONV) {
	    /* try taking the final value regardless */
	    err = 0;
	} else {
	    fprintf(stderr, "ar1_mle: BFGS_max gave err=%d "
		    "(incoming rho %g, final %g)\n",
		    err, *rho, theta[0]);
	}
    }

    if (!err) {
#if CL_DEBUG
	fprintf(stderr, "ar1_mle, rho %g -> %g\n", *rho, theta[0]);
#endif
	*rho = theta[0];
    }

    free(theta);

    return err;
}

/* Produce W = CVC' without storing C */

static void make_CVC (gretl_matrix *W, const gretl_matrix *VC,
		      int s, int agg)
{
    double wij;
    int N = W->rows;
    int i, j, k, ii;

    if (agg >= AGG_EOP) {
	/* select rows of VC */
	for (j=0; j<N; j++) {
	    ii = agg == AGG_SOP ? 0 : s-1;
	    for (i=0; i<N; i++) {
		wij = gretl_matrix_get(VC, ii, j);
		gretl_matrix_set(W, i, j, wij);
		gretl_matrix_set(W, j, i, wij);
		ii += s;
	    }
	}
    } else {
	/* cumulate rows of VC */
	for (j=0; j<N; j++) {
	    ii = 0;
	    for (i=0; i<N; i++) {
		wij = 0.0;
		for (k=0; k<s; k++) {
		    wij += gretl_matrix_get(VC, ii++, j);
		}
		gretl_matrix_set(W, i, j, wij);
		gretl_matrix_set(W, j, i, wij);
	    }
	}
    }
}

/* Build VC' without storing V or C as such */

static void make_VC (gretl_matrix *VC, int N,
		     int s, double a, int agg)
{
    int sN = s*N;
    int i, j, k;

    if (agg >= AGG_EOP) {
	k = agg == AGG_SOP ? 0 : s-1;
	for (j=0; j<N; j++) {
	    for (i=0; i<sN; i++) {
		gretl_matrix_set(VC, i, j, pow(a, abs(k-i)));
	    }
	    k += s;
	}
    } else {
	double vij;
	int kk = 0;

	for (j=0; j<N; j++) {
	    for (i=0; i<sN; i++) {
		vij = 0;
		for (k=0; k<s; k++) {
		    vij += pow(a, abs(kk+k-i));
		}
		gretl_matrix_set(VC, i, j, vij);
	    }
	    kk += s;
	}
    }
}

/* Make the counterpart to VC' that's required for extrapolation;
   see p. 375 in Chow and Lin 1971.
*/

static void make_EVC (gretl_matrix *EVC, int s,
		      double a, int agg)
{
    double vij;
    int m = EVC->rows;
    int N = EVC->cols;
    int sN = s * N;
    int i, j, k, p;

    if (agg >= AGG_EOP) {
	p = agg == AGG_SOP ? sN : s*(N-1)+1;
	for (j=0; j<N; j++) {
	    vij = pow(a, p);
	    for (i=0; i<m; i++) {
		gretl_matrix_set(EVC, i, j, vij);
		vij *= a;
	    }
	    p -= s;
	}
    } else {
	p = sN;
	for (j=0; j<N; j++) {
	    vij = 0;
	    for (k=0; k<s; k++) {
		vij += pow(a, p-k);
	    }
	    for (i=0; i<m; i++) {
		gretl_matrix_set(EVC, i, j, vij);
		vij *= a;
	    }
	    p -= s;
	}
    }
}

/* Multiply VC' into W*u and increment y by the result */

static void multiply_by_VC (gretl_matrix *y,
			    const gretl_matrix *VC,
			    const gretl_matrix *wu,
			    int s, int m, double a,
			    int agg)
{
    gretl_matrix *EVC = NULL;
    gretl_matrix ext = {0};
    int sN = y->rows;

    if (m > 0) {
	/* we have some extrapolation to do */
	sN -= m;
	gretl_matrix_reuse(y, sN, 1);
	gretl_matrix_init(&ext);
	ext.rows = m;
	ext.cols = 1;
	ext.val = y->val + sN;
	EVC = gretl_matrix_alloc(m, sN/s);
    }

    gretl_matrix_multiply_mod(VC, GRETL_MOD_NONE,
			      wu, GRETL_MOD_NONE,
			      y, GRETL_MOD_CUMULATE);

    if (m > 0) {
	make_EVC(EVC, s, a, agg);
	gretl_matrix_multiply_mod(EVC, GRETL_MOD_NONE,
				  wu, GRETL_MOD_NONE,
				  &ext, GRETL_MOD_CUMULATE);
	gretl_matrix_reuse(y, sN+m, 1);
	gretl_matrix_free(EVC);
    }
}

struct gls_info {
    const gretl_matrix *y0;
    gretl_matrix *CX;
    gretl_matrix *VC;
    gretl_matrix *W;
    gretl_matrix *Z;
    gretl_matrix *Tmp1;
    gretl_matrix *Tmp2;
    gretl_matrix *b;
    gretl_matrix *u;
    int s, det, agg;
    int method;
    double lnl;
};

/* Carry out enough of the Chow-Lin GLS calculations to permit
   computation of a figure of merit (loglikelihood or sum of squared
   GLS residuals), conditional on rho.  This is used if we already
   have a "final" rho value greater than zero, or if we're trying to
   find optimal rho via GLS.
*/

static double cl_gls_calc (const double *rho, void *data)
{
    struct gls_info *G = data;
    gretl_matrix *Wcpy;
    double a = *rho;
    double ldet, crit;
    double SSR = NADBL;
    int N = G->y0->rows;
    int err = 0;

    make_VC(G->VC, N, G->s, a, G->agg);
    make_CVC(G->W, G->VC, G->s, G->agg);
    if (G->method == R_SSR) {
	gretl_matrix_multiply_by_scalar(G->W, 1/(1.0-a*a));
    }
    Wcpy = gretl_matrix_copy(G->W);
    err = gretl_invert_symmetric_matrix(G->W);

    if (!err) {
	gretl_matrix_qform(G->CX, GRETL_MOD_TRANSPOSE,
			   G->W, G->Z, GRETL_MOD_NONE);
	err = gretl_invert_symmetric_matrix(G->Z);
    }

    if (!err) {
	/* GLS coefficients */
	gretl_matrix_multiply_mod(G->Z, GRETL_MOD_NONE,
				  G->CX, GRETL_MOD_TRANSPOSE,
				  G->Tmp1, GRETL_MOD_NONE);
	gretl_matrix_multiply(G->Tmp1, G->W, G->Tmp2);
	gretl_matrix_multiply(G->Tmp2, G->y0, G->b);

	/* GLS residuals */
	gretl_matrix_copy_values(G->u, G->y0);
	gretl_matrix_multiply_mod(G->CX, GRETL_MOD_NONE,
				  G->b, GRETL_MOD_NONE,
				  G->u, GRETL_MOD_DECREMENT);

	ldet = gretl_matrix_log_determinant(Wcpy, &err);
    }

    if (!err) {
	SSR = gretl_scalar_qform(G->u, G->W, &err);
	if (!err) {
	    G->lnl = -0.5*N - 0.5*N*LN_2_PI -N*log(SSR/N)/2 - ldet/2;
	}
    }

    gretl_matrix_free(Wcpy);

    if (err) {
	crit = G->lnl = NADBL;
    } else {
	crit = G->method == R_SSR ? -SSR : G->lnl;
    }

    return crit;
}

/* Driver for optimization of rho via GLS, either on the criterion
   of likelihood or the straight GLS criterion.  In line with other
   implementations we place an upper bound on rho short of 1.0.
*/

static int cl_gls_max (double *a, struct gls_info *G,
		       PRN *prn)
{
    double theta[] = {0.5, 0, 0.999};
    int err, iters = 0;

    err = gretl_gss(theta, 1.0e-12, &iters, cl_gls_calc,
		    G, OPT_NONE, prn);

#if CL_DEBUG
    fprintf(stderr, "cl_gls_max: err=%d, iters %d\n", err, iters);
#endif

    if (!err) {
	double r = theta[0];

	*a = r;
	if (G->method == R_SSR) {
	    gretl_matrix_multiply_by_scalar(G->W, 1/(1-r*r));
	}
    }

    return err;
}

/* Regressor matrix: we put in constant (if det > 0) plus linear trend
   (if det > 1) and squared trend (if det = 3), summed appropriately
   based on @s.

   If the user has supplied high-frequency covariates in @X, we
   compress them from column @det onward.

   Note: this version of the implicit C matrix assumes what Chow and
   Lin call "distribution", which is appropriate for flow variables.
*/

static void fill_CX (gretl_matrix *CX, int s, int det,
		     const gretl_matrix *X)
{
    double xt1, xt2;
    int i, j, k = 1;
    int t, r = 0;

    for (t=0; t<CX->rows; t++) {
	if (det > 0) {
	    gretl_matrix_set(CX, t, 0, s);
	    if (det > 1) {
		xt1 = xt2 = 0.0;
		for (i=0; i<s; i++) {
		    xt1 += k;
		    if (det > 2) {
			xt2 += k * k;
		    }
		    k++;
		}
		gretl_matrix_set(CX, t, 1, xt1);
		if (det > 2) {
		    gretl_matrix_set(CX, t, 2, xt2);
		}
	    }
	}
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		xt1 = 0.0;
		for (i=0; i<s; i++) {
		    xt1 += gretl_matrix_get(X, r + i, j);
		}
		gretl_matrix_set(CX, t, det+j, xt1);
	    }
	    r += s;
	}
    }
}

/* Variant of fill_CX() in which C is a selection matrix,
   for interpolation in the strict sense (stock variables).
*/

static void fill_CX2 (gretl_matrix *CX, int s, int det,
		      const gretl_matrix *X, int agg)
{
    double xkj;
    int i, j, r, t;

    gretl_matrix_zero(CX);
    r = (agg == AGG_SOP)? 0 : s-1;

    for (i=0; i<CX->rows; i++) {
	if (det > 0) {
	    gretl_matrix_set(CX, i, 0, 1);
	    if (det > 1) {
		t = r + 1;
		gretl_matrix_set(CX, i, 1, t);
		if (det > 2) {
		    gretl_matrix_set(CX, i, 2, t*t);
		}
	    }
	}
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		xkj = gretl_matrix_get(X, r, j);
		gretl_matrix_set(CX, i, det+j, xkj);
	    }
	}
	r += s;
    }
}

static void make_X_beta (gretl_vector *y, const double *b,
			 const gretl_matrix *X, int det)
{
    int i, j, t;

    for (i=0; i<y->rows; i++) {
	if (det > 0) {
	    y->val[i] = b[0];
	    if (det > 1) {
		t = i + 1;
		y->val[i] += b[1]*t;
		if (det > 2) {
		    y->val[i] += b[2]*t*t;
		}
	    }
	}
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		y->val[i] += b[det+j] * gretl_matrix_get(X, i, j);
	    }
	}
    }
}

/* first-order autocorrelation of residuals: see also
   rhohat() in estimate.c.
*/

static double acf_1 (const gretl_matrix *y,
		     const gretl_matrix *X,
		     const gretl_matrix *b,
		     const gretl_matrix *u)
{
    double rho, num = 0, den = 0;
    int t;

    for (t=0; t<u->rows; t++) {
	den += u->val[t] * u->val[t];
	if (t > 0) {
	    num += u->val[t] * u->val[t-1];
	}
    }

    if (num < 1.0e-9) {
	return 0;
    }

    rho = num / den;

    /* improve the initial estimate of @rho via ML */
    ar1_mle(y, X, b, sqrt(den / u->rows), &rho);

    return rho;
}

/* We come here if (a) the caller has specified a fixed rho
   of zero, or (b) we're doing the "traditional" Chow-Lin
   thing of estimating low-frequency rho from an initial
   OLS regression and converting to high-frequency.
*/

static int cl_ols (struct gls_info *G,
		   const gretl_matrix *y0,
		   const gretl_matrix *X,
		   gretl_matrix *y,
		   double *rho,
		   PRN *prn)
{
    gretl_matrix *u;
    double a = 0;
    int err;

    u = G->method == R_FIXED ? NULL : G->u;
    err = gretl_matrix_ols(y0, G->CX, G->b, NULL, u, NULL);
    if (err) {
	return err;
    }

    if (G->method == R_ACF1) {
	a = acf_1(y0, G->CX, G->b, G->u);
	if (a <= 0.0) {
	    a = 0; /* don't pursue negative @a */
	} else if (G->agg >= AGG_EOP) {
	    a = pow(a, 1.0/G->s);
	} else {
	    double bracket[] = {0, 0.9999};
	    struct chowlin cl = {G->s, a};

	    err = gretl_fzero(bracket, 1.0e-12,
			      chow_lin_callback, &cl,
			      &a, OPT_NONE, prn);
#if CL_DEBUG
	    fprintf(stderr, "gretl_fzero: err=%d, a=%g\n", err, a);
#endif
	}
    }

    if (a == 0) {
	/* GLS isn't needed */
	make_X_beta(y, G->b->val, X, G->det);
	if (G->agg == AGG_AVG) {
	    gretl_matrix_multiply_by_scalar(y, G->s);
	}
    }

    *rho = a;

    return err;
}

static void show_GLS_results (const struct gls_info *G,
			      double a, int det,
			      PRN *prn)
{
    const char *dnames[] = {"const", "trend", "trend^2"};
    int i;

    pputs(prn, "\nGLS coefficients:\n");
    for (i=0; i<G->b->rows; i++) {
	if (i < det) {
	    pprintf(prn, " %-8s", dnames[i]);
	} else {
	    pprintf(prn, " %c%-7d", 'X', i-det+1);
	}
	pprintf(prn, "%#g\n", G->b->val[i]);
    }
    pprintf(prn, " %-8s%#.8g\n", "rho", a);
    pprintf(prn, " loglikelihood = %.8g\n", G->lnl);
}

/**
 * chow_lin_disagg:
 * @Y0: N x k: holds the original data to be expanded.
 * @X: (optionally) holds covariates of Y at the higher frequency;
 * if these are supplied they supplement the deterministic
 * terms (if any) as signalled by @det.
 * @s: the expansion factor: 3 for quarterly to monthly,
 * 4 for annual to quarterly or 12 for annual to monthly.
 * @det: 0 for none, 1 for constant, 2 for linear trend, 3 for
 * quadratic trend.
 * @agg: aggregation type.
 * @err: location to receive error code.
 *
 * Distribute or interpolate via the method of Chow and Lin. See
 * See Gregory C. Chow and An-loh Lin, "Best Linear Unbiased
 * Interpolation, Distribution, and Extrapolation of Time Series
 * by Related Series", Review of Economics and Statistics, Vol. 53,
 * No. 4 (November 1971) pp. 372-375.
 *
 * If @X is given it must have @s * N rows.
 *
 * Returns: matrix containing the expanded series, or
 * NULL on failure.
 */

static gretl_matrix *chow_lin_disagg (const gretl_matrix *Y0,
				      const gretl_matrix *X,
				      int s, int det, int agg,
				      int method, double rho,
				      PRN *prn, int *perr)
{
    struct gls_info G = {0};
    gretl_matrix_block *B;
    gretl_matrix *Y;
    gretl_matrix *y0, *y;
    gretl_matrix my0, my;
    int ny = Y0->cols;
    int nx = det;
    int N = Y0->rows;
    int sN = s * N;
    int m = 0;
    int i, err = 0;

    if (X != NULL) {
	nx += X->cols;
	/* m: "extra" observations */
	m = X->rows - sN;
    }
    if (nx == 0) {
	/* nothing to work with! */
	*perr = E_ARGS;
	return NULL;
    }

    /* the return value */
    Y = gretl_zero_matrix_new(sN + m, ny);
    if (Y == NULL) {
	*perr = E_ALLOC;
	return NULL;
    }

    /* workspace */
    B = gretl_matrix_block_new(&G.CX, N, nx,
			       &G.VC, sN, N,
			       &G.W, N, N,
			       &G.b, nx, 1,
			       &G.u, N, 1,
			       &G.Z, nx, nx,
			       &G.Tmp1, nx, N,
			       &G.Tmp2, nx, N,
			       NULL);
    if (B == NULL) {
	*perr = E_ALLOC;
	gretl_matrix_free(Y);
	return NULL;
    }

    if (!na(rho)) {
	method = R_FIXED;
    }

    /* regressors: deterministic terms (as wanted), plus
       anything else the user has added
    */
    if (agg >= AGG_EOP) {
	fill_CX2(G.CX, s, det, X, agg);
    } else {
	fill_CX(G.CX, s, det, X);
    }
#if CL_DEBUG > 1
    gretl_matrix_print(G.CX, "CX");
#endif

    /* original y0 vector, length N */
    gretl_matrix_init(&my0);
    y0 = &my0;
    y0->rows = N;
    y0->cols = 1;
    y0->val = Y0->val;

    /* y vector for return, length sN + m */
    gretl_matrix_init(&my);
    y = &my;
    y->rows = sN + m;
    y->cols = 1;
    y->val = Y->val;

    /* hook up various non-malloc'd things */
    G.y0 = y0;
    G.s = s;
    G.det = det;
    G.agg = agg;
    G.method = method;

    for (i=0; i<ny && !err; i++) {
	double a = method == R_FIXED ? rho : 0.0;

	if (i > 0) {
	    /* pick up the current columns for reading and writing */
	    y0->val = Y0->val + i * N;
	    y->val = Y->val + i * (sN + m);
	}

	if (method == R_ACF1 || (method == R_FIXED && a <= 0)) {
	    err = cl_ols(&G, y0, X, y, &a, prn);
	    if (a == 0) {
		/* this iteration is handled */
		continue;
	    }
	}

	if (!err) {
	    if (method == R_MLE || method == R_SSR) {
		err = cl_gls_max(&a, &G, prn);
	    } else {
		cl_gls_calc(&a, &G);
	    }
	}

	if (!err) {
	    /* y = X*beta + V*C'*W*u */
	    make_X_beta(y, G.b->val, X, det);
	    gretl_matrix_reuse(G.Tmp1, N, 1);
	    gretl_matrix_multiply(G.W, G.u, G.Tmp1);
	    multiply_by_VC(y, G.VC, G.Tmp1, s, m, a, agg);
	    gretl_matrix_reuse(G.Tmp1, nx, N);

	    if (prn != NULL && gretl_messages_on()) {
		show_GLS_results(&G, a, det, prn);
	    }

	    if (agg == AGG_AVG) {
		gretl_matrix_multiply_by_scalar(y, s);
	    }
	}
    }

    *perr = err;
    gretl_matrix_block_destroy(B);

    return Y;
}

/* The method of F. T. Denton, "Adjustment of Monthly or Quarterly
   Series to Annual Totals: An Approach Based on Quadratic
   Minimization", Journal of the American Statistical Association
   Vol. 66, No. 333 (March 1971), pp. 99-102, proportional first
   difference variant, as modified by P. A. Cholette, "Adjusting
   Sub-annual Series to Yearly Benchmarks," Survey Methodology,
   Vol. 10, 1984, pp. 35-49.

   The solution method is based on Tommaso Di Fonzo and Marco Marini,
   "On the Extrapolation with the Denton Proportional Benchmarking
   Method", IMF Working Paper WP/12/169, 2012.
*/

static gretl_matrix *denton_pfd (const gretl_vector *y0,
				 const gretl_vector *p,
				 int s, int agg,
				 int *err)
{
    gretl_matrix *M;
    gretl_matrix *y;
    gretl_matrix *tmp;
    int N = y0->rows;
    int sN = p->rows;
    int sNN = sN + N;
    int i, j, k = 0;
    int offset;

    /* we need one big matrix, @M */
    M = gretl_zero_matrix_new(sNN, sNN);
    tmp = gretl_matrix_alloc(sN, N);
    y = gretl_matrix_alloc(sN, 1);

    if (M == NULL || tmp == NULL || y == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* In @M, create (D'D ~ diag(p)*J') | (J*diag(p) ~ 0);
       see di Fonzo and Marini, equation (4)
    */

    /* the upper left portion, D'D */
    for (i=0; i<sN; i++) {
	gretl_matrix_set(M, i, i, (i == 0 || i == sN-1)? 1 : 2);
	if (i > 0) {
	    gretl_matrix_set(M, i, i-1, -1);
	}
	if (i < sN-1) {
	    gretl_matrix_set(M, i, i+1, -1);
	}
    }

    /* the bottom and right portions, using @p */
    k = offset = (agg == AGG_EOP)? s-1 : 0;
    for (i=sN; i<sNN; i++) {
	if (agg >= AGG_EOP) {
	    gretl_matrix_set(M, i, offset, p->val[k]);
	    gretl_matrix_set(M, offset, i, p->val[k]);
	    k += s;
	} else {
	    for (j=offset; j<offset+s; j++) {
		gretl_matrix_set(M, i, j, p->val[k]);
		gretl_matrix_set(M, j, i, p->val[k]);
		k++;
	    }
	}
	offset += s;
    }

    *err = gretl_invert_symmetric_indef_matrix(M);

    if (*err) {
	gretl_matrix_free(y);
	y = NULL;
    } else {
	/* extract the relevant portion of M-inverse and
	   premultiply by (diag(p) ~ 0) | (0 ~ I)
	*/
	double mij;

	for (j=0; j<N; j++) {
	    for (i=0; i<sN; i++) {
		mij = gretl_matrix_get(M, i, j+sN);
		gretl_matrix_set(tmp, i, j, mij * p->val[i]);
	    }
	}
	gretl_matrix_multiply(tmp, y0, y);
    }

    if (agg == AGG_AVG) {
	gretl_matrix_multiply_by_scalar(y, s);
    }

    gretl_matrix_free(M);
    gretl_matrix_free(tmp);

    return y;
}

#if 0 /* not yet */

static int get_aggregation_type (const char *s, int *err)
{
    if (!strcmp(s, "sum")) {
	return 0;
    } else if (!strcmp(s, "avg")) {
	return 1;
    } else if (!strcmp(s, "last")) {
	return 2;
    } else if (!strcmp(s, "first")) {
	return 3;
    } else {
	*err = E_INVARG;
	return -1;
    }
}

static int get_tdisagg_method (const char *s, int *err)
{
    if (!strcmp(s, "chow-lin")) {
	return 0;
    } else if (!strcmp(s, "chow-lin-mle")) {
	return 1;
    } else if (!strcmp(s, "chow-lin-ssr")) {
	return 2;
    } else if (!strcmp(s, "denton")) {
	return 3;
    } else {
	*err = E_INVARG;
	return -1;
    }
}

static int tdisagg_get_options (const gretl_bundle *b)
{
    double rho;
    const char *str;
    int s, agg, method, det;
    int err = 0;

    if (gretl_bundle_has_key(b, "s")) {
	s = gretl_bundle_get_int(b, "s", &err);
	if (!err && s != 3 && s != 4 && s != 12) {
	    err = E_INVARG;
	}
    }
    if (!err && gretl_bundle_has_key(b, "agg")) {
	str = gretl_bundle_get_string(b, "agg", &err);
	if (!err) {
	    agg = get_aggregation_type(str, &err);
	}
    }
    if (!err && gretl_bundle_has_key(b, "method")) {
	str = gretl_bundle_get_string(b, "method", &err);
	if (!err) {
	    method = get_tdisagg_method(str, &err);
	}
    }
    if (!err && gretl_bundle_has_key(b, "det")) {
	det = gretl_bundle_get_int(b, "det", &err);
	if (!err && (det < 0 || det > 3)) {
	    err = E_INVARG;
	}
    }
    if (!err && gretl_bundle_has_key(b, "rho")) {
	rho = gretl_bundle_get_scalar(b, "rho", &err);
	if (!err && (rho <= -1.0 || rho >= 1.0)) {
	    err = E_INVARG;
	}
    }

    return err;
}

#endif /* not yet */

gretl_matrix *time_disaggregate (const gretl_matrix *Y0,
				 const gretl_matrix *X,
				 int s, int agg, int method,
				 int det, double rho,
				 PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;

    if (method < 3) {
	/* Chow-Lin variants */
	if (det < 0) {
	    det = X == NULL ? 2 : 1;
	}
	ret = chow_lin_disagg(Y0, X, s, det, agg, method, rho, prn, err);
    } else if (method == 3) {
	/* Modified Denton, proportional first differences */
	int ylen = gretl_vector_get_length(Y0);
	gretl_matrix *X0 = NULL;

	if (X == NULL) {
	    /* as per R, tempdisagg, use a constant if
	       X is not provided */
	    X0 = gretl_unit_matrix_new(s * ylen, 1);
	    X = X0;
	} else {
	    int xlen = gretl_vector_get_length(X);

	    if (ylen == 0 || xlen == 0 || xlen != s * ylen) {
		*err = E_INVARG;
	    }
	}
	if (!*err) {
	    ret = denton_pfd(Y0, X, s, agg, err);
	    gretl_matrix_free(X0);
	}
    } else {
	/* no other choices at present */
	*err = E_INVARG;
    }

    return ret;
}
