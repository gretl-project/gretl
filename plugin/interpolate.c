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

/* Note: at present this plugin only offers Chow-Lin interpolation.
   At some point we may want to add interpolation using the Kalman
   filter or some other more sophisticated variant(s).
*/

#include "libgretl.h"
#include "version.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"

#define CL_DEBUG 0
#define AR1_MLE 1

enum {
    AGG_AVG, /* average */
    AGG_SUM, /* sum */
    AGG_SOP, /* start of period */
    AGG_EOP  /* end of period */
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

#if AR1_MLE

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
#if 0 /* debugging */
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    gretlopt opt = OPT_V;
#else
    PRN *prn = NULL;
    gretlopt opt = OPT_NONE;
#endif
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
		   NULL, &a, NULL, opt, prn);

    if (err) {
	if (err == E_NOCONV) {
	    /* try taking the final value regardless?? */
	    err = 0;
	} else {
	    fprintf(stderr, "ar1_mle: BFGS_max gave err=%d (incoming rho %g, final %g)\n",
		    err, *rho, theta[0]);
	}
    }

    if (!err) {
#if 1 || CL_DEBUG
	/* let's show this for now */
	fprintf(stderr, "ar1_mle, rho %g -> %g\n", *rho, theta[0]);
#endif
	*rho = theta[0];
    }

    free(theta);
    gretl_print_destroy(prn);

    return err;
}

#endif /* AR1_MLE */

static double csum (int n, double a, int k)
{
    double s = 0.0;
    int i;

    for (i=0; i<n; i++) {
	s += pow(a, abs(k));
	k++;
    }

    return s;
}

/* Generate W = CVC' without storing the full C or V matrices. C is
   the matrix that transforms from higher to lower frequency by
   summation; V is the autocovariance matrix for AR(1) disturbances
   with autoregressive coefficient @a; @n is the expansion factor.
*/

static void make_CVC (gretl_matrix *W, int n, double a)
{
    double wij;
    int i, j, k, m;

    for (i=0; i<W->rows; i++) {
	m = 0;
	for (j=i; j<W->cols; j++) {
	    wij = 0.0;
	    for (k=0; k<n; k++) {
		wij += csum(n, a, m--);
	    }
	    gretl_matrix_set(W, i, j, wij);
	    gretl_matrix_set(W, j, i, wij);
	}
    }
}

/* Variant of make_CVC() in which C is the selection matrix for
   interpolation, selecting either the first or the last sub-period.
*/

static void make_CVC2 (gretl_matrix *W, int n, double a, int agg)
{
    double wij;
    int i, j, s;

    s = (agg == AGG_SOP)? 0 : 1;

    for (i=0; i<W->rows; i++) {
	gretl_matrix_set(W, i, i, 1.0);
	for (j=0; j<i; j++) {
	    wij = pow(a, abs(n*(j+s) - n*(i+s)));
	    gretl_matrix_set(W, i, j, wij);
	    gretl_matrix_set(W, j, i, wij);
	}
    }
}

/* Multiply VC' into u and increment yx by the result;
   again, without storing V or C'.
*/

static void mult_VC (gretl_matrix *yx, gretl_matrix *u,
		     int n, double a)
{
    int Tx = yx->rows;
    int T = u->rows;
    int t, j;

    for (t=0; t<Tx; t++) {
	for (j=0; j<T; j++) {
	    yx->val[t] += u->val[j] * csum(n, a, j * n - t);
	}
    }
}

/* Regressor matrix: we put in constant (if det > 0) plus linear trend
   (if det > 1) and squared trend (if det = 3), summed appropriately
   based on @n.

   If the user has supplied high-frequency covariates in @X, we
   compress them from column @det onward.

   Note: this version of the implicit C matrix assumes what Chow and
   Lin call "distribution", which is appropriate for flow variables.
*/

static void fill_CX (gretl_matrix *CX, int n, int det,
		     const gretl_matrix *X)
{
    double xt1, xt2;
    int i, j, k = 1;
    int t, s = 0;

    for (t=0; t<CX->rows; t++) {
	if (det > 0) {
	    gretl_matrix_set(CX, t, 0, n);
	    if (det > 1) {
		xt1 = xt2 = 0.0;
		for (i=0; i<n; i++) {
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
		for (i=0; i<n; i++) {
		    xt1 += gretl_matrix_get(X, s + i, j);
		}
		gretl_matrix_set(CX, t, det+j, xt1);
	    }
	    s += n;
	}
    }
}

/* Variant of fill_CX() in which C is a selection matrix,
   for interpolation in the strict sense (stock variables).
*/

static void fill_CX2 (gretl_matrix *CX, int n, int det,
		      const gretl_matrix *X, int agg)
{
    double xt;
    int j, s, t;

    gretl_matrix_zero(CX);
    s = (agg == AGG_SOP)? 0 : n-1;

    for (t=0; t<CX->rows; t++) {
	if (det > 0) {
	    gretl_matrix_set(CX, t, 0, 1);
	    if (det > 1) {
		gretl_matrix_set(CX, t, 1, s);
		if (det > 2) {
		    gretl_matrix_set(CX, t, 2, s*s);
		}
	    }
	}
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		xt = gretl_matrix_get(X, s, j);
		gretl_matrix_set(CX, t, det+j, xt);
	    }
	}
	s += n;
    }
}

static void make_Xx_beta (gretl_vector *y, const double *b,
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

#if AR1_MLE
    /* improve the estimate of @rho via ML? */
    ar1_mle(y, X, b, sqrt(den / u->rows), &rho);
#endif

    return rho;
}

/**
 * chow_lin_disagg:
 * @Y: T x k: holds the original data to be expanded.
 * @X: (optionally) holds covariates of @Y at the higher frequency;
 * if these are supplied they supplement the deterministic
 * terms (if any) as signalled by @det.
 * @xfac: the expansion factor: 3 for quarterly to monthly,
 * 4 for annual to quarterly or 12 for annual to monthly.
 * @det: 0 for none, 1 for constant, 2 for linear trend, 3 for
 * quadratic trend.
 * @err: location to receive error code.
 *
 * Distribute or interpolate via the method of Chow and Lin. See
 * See Gregory C. Chow and An-loh Lin, "Best Linear Unbiased
 * Interpolation, Distribution, and Extrapolation of Time Series
 * by Related Series", Review of Economics and Statistics, Vol. 53,
 * No. 4 (November 1971) pp. 372-375.
 *
 * If @X is given it must have T * @xfac rows.
 *
 * Returns: matrix containing the expanded series, or
 * NULL on failure.
 */

static gretl_matrix *chow_lin_disagg (const gretl_matrix *Y,
				      const gretl_matrix *X,
				      int xfac, int det,
				      PRN *prn, int *err)
{
    gretl_matrix_block *B;
    gretl_matrix *CX, *b, *u, *W, *Z;
    gretl_matrix *Tmp1, *Tmp2;
    gretl_matrix *Yx = NULL;
    gretl_matrix *y, *yx;
    gretl_matrix my, myx;
    int agg = AGG_AVG;
    int nx, ny = Y->cols;
    int T = Y->rows;
    int Tx = T * xfac;
    int i;

    /* Note: checks to the effect that xfac = 3, 4 or 12 and,
       if X is non-NULL, that X->rows = xfac * Y->rows,
       should have already been performed.
    */

    nx = det;
    if (X != NULL) {
	nx += X->cols;
    }

    if (nx == 0) {
	/* nothing to work with! */
	*err = E_ARGS;
	return NULL;
    }

    gretl_matrix_init(&my);
    gretl_matrix_init(&myx);

    /* the return value */
    Yx = gretl_zero_matrix_new(Tx, ny);
    if (Yx == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* block of low-frequency matrices */
    B = gretl_matrix_block_new(&CX, T, nx,
			       &W, T, T,
			       &b, nx, 1,
			       &u, T, 1,
			       &Z, nx, nx,
			       &Tmp1, nx, T,
			       &Tmp2, nx, T,
			       NULL);
    if (B == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(Yx);
	return NULL;
    }

    /* regressors: deterministic terms (as wanted), plus
       anything the user has added
    */
    if (agg >= AGG_SOP) {
	fill_CX2(CX, xfac, det, X, agg);
    } else {
	fill_CX(CX, xfac, det, X);
    }
#if CL_DEBUG > 1
    gretl_matrix_print(CX, "CX");
#endif

    y = &my;
    yx = &myx;
    y->rows = T;
    yx->rows = Tx;
    y->cols = yx->cols = 1;
    y->val = Y->val;
    yx->val = Yx->val;

    for (i=0; i<ny; i++) {
	double a = 0.0;

	if (i > 0) {
	    /* pick up the current column */
	    y->val = Y->val + i * T;
	    yx->val = Yx->val + i * Tx;
	}

	/* initial low-frequency OLS */
	*err = gretl_matrix_ols(y, CX, b, NULL, u, NULL);

	if (!*err) {
	    a = acf_1(y, CX, b, u);
#if CL_DEBUG
	    fprintf(stderr, "initial acf_1 = %g\n", a);
#endif
	    if (a <= 0.0) {
		/* don't pursue negative @a */
		make_Xx_beta(yx, b->val, X, det);
		gretl_matrix_multiply_by_scalar(yx, xfac);
		/* nothing more to do, this iteration */
		continue;
	    } else {
		double bracket[] = {0, 0.9999};
		struct chowlin cl = {xfac, a};

		*err = gretl_fzero(bracket, 1.0e-12,
				   chow_lin_callback, &cl,
				   &a, OPT_NONE, prn);
#if CL_DEBUG
		fprintf(stderr, "gretl_fzero: err=%d, a=%g\n", *err, a);
#endif
	    }
	}

	if (!*err) {
	    if (agg >= AGG_SOP) {
		make_CVC2(W, xfac, a, agg);
	    } else {
		make_CVC(W, xfac, a);
	    }
	    *err = gretl_invert_symmetric_matrix(W);
	}

	if (!*err) {
	    gretl_matrix_qform(CX, GRETL_MOD_TRANSPOSE,
			       W, Z, GRETL_MOD_NONE);
	    *err = gretl_invert_symmetric_matrix(Z);
	}

	if (!*err) {
	    /* GLS \hat{\beta} */
	    gretl_matrix_multiply_mod(Z, GRETL_MOD_NONE,
				      CX, GRETL_MOD_TRANSPOSE,
				      Tmp1, GRETL_MOD_NONE);
	    gretl_matrix_multiply(Tmp1, W, Tmp2);
	    gretl_matrix_multiply(Tmp2, y, b);

	    /* X(expanded) * \hat{\beta} */
	    make_Xx_beta(yx, b->val, X, det);

	    /* GLS residuals */
	    gretl_matrix_copy_values(u, y);
	    gretl_matrix_multiply_mod(CX, GRETL_MOD_NONE,
				      b, GRETL_MOD_NONE,
				      u, GRETL_MOD_DECREMENT);

	    /* yx = Xx*beta + V*C'*W*u */
	    gretl_matrix_reuse(Tmp1, T, 1);
	    gretl_matrix_multiply(W, u, Tmp1);
	    mult_VC(yx, Tmp1, xfac, a);
	    gretl_matrix_reuse(Tmp1, nx, T);

	    /* FIXME make this conditional */
	    gretl_matrix_multiply_by_scalar(yx, xfac);
	}
    }

    gretl_matrix_block_destroy(B);

    return Yx;
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
				 int s, int *err)
{
    gretl_matrix *M;
    gretl_matrix *ret;
    gretl_matrix *tmp;
    int N = y0->rows;
    int sN = p->rows;
    int sNN = sN + N;
    int i, j, k = 0;
    int offset;

    /* we need one big matrix, @M */
    M = gretl_zero_matrix_new(sNN, sNN);
    tmp = gretl_matrix_alloc(sN, N);
    ret = gretl_matrix_alloc(sN, 1);

    if (M == NULL || tmp == NULL || ret == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* In @M, create (D'D ~ diag(p)*J') | (J*diag(p) ~ 0);
       see di Fonzo and Marini, equation (4)
    */
    for (i=0; i<sN; i++) {
	/* upper left portion, D'D */
	gretl_matrix_set(M, i, i, (i == 0 || i == sN-1)? 1 : 2);
	if (i > 0) {
	    gretl_matrix_set(M, i, i-1, -1);
	}
	if (i < sN-1) {
	    gretl_matrix_set(M, i, i+1, -1);
	}
    }
    k = offset = 0;
    for (i=sN; i<sNN; i++) {
	/* bottom and right portions, using @p */
	for (j=offset; j<offset+s; j++) {
	    gretl_matrix_set(M, i, j, p->val[k]);
	    gretl_matrix_set(M, j, i, p->val[k]);
	    k++;
	}
	offset += s;
    }

    *err = gretl_invert_symmetric_indef_matrix(M);

    if (*err) {
	gretl_matrix_free(ret);
	ret = NULL;
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
	gretl_matrix_multiply(tmp, y0, ret);
    }

    gretl_matrix_free(M);
    gretl_matrix_free(tmp);

    return ret;
}

gretl_matrix *time_disaggregate (const gretl_matrix *Y,
				 const gretl_matrix *X,
				 int xfac, int det,
				 int method, PRN *prn,
				 int *err)
{
    if (method == 0) {
	/* Chow-Lin */
	return chow_lin_disagg(Y, X, xfac, det, prn, err);
    } else if (method == 1) {
	/* Modified Denton, proportional first differences */
	if (gretl_vector_get_length(Y) == 0 ||
	    gretl_vector_get_length(X) == 0) {
	    *err = E_INVARG;
	    return NULL;
	}
	return denton_pfd(Y, X, xfac, err);
    } else {
	/* no other options at present */
	*err = E_INVARG;
	return NULL;
    }
}
