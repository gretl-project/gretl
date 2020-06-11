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

struct chowlin {
    int n;
    double targ;
};

/* Callback for fzero(), as we adjust the coefficient @a so the
   theoretically derived ratio of polynomials in @a matches the
   empirical first-order autocorrelation of the OLS residuals
   (cl->targ). Return the negative of the squared residual.
*/

static double chow_lin_callback (double a, void *p)
{
    struct chowlin *cl = (struct chowlin *) p;
    double r, num, den, resid;

    if (a == 0) {
	r = 0;
    } else {
	if (cl->n == 3) {
	    num = a + 2*a*a + 3*pow(a, 3) + 2*pow(a, 4) + pow(a, 5);
	    den = 3 + 4*a + 2*a*a;
	} else {
	    /* n = 4: requires cl->targ > 0 */
	    num = a + 2*a*a + 3*pow(a, 3) + 4*pow(a, 4) + 3*pow(a, 5)
		+ 2*pow(a, 6) + pow(a, 7);
	    den = 4 + 6*a + 4*a*a + 2*pow(a, 3);
	}
	r = num/den;
    }

    resid = r - cl->targ;

#if CL_DEBUG > 1
    fprintf(stderr, "chow_lin_callback: target %g, a %g residual %g\n",
	    cl->targ, a, resid);
#endif

    return resid;
}

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

/* Generate W = CVC' without storing the full C or V matrices.  C is
   the selection matrix that transforms from higher frequency to lower
   frequency by summation; V is the autocovariance matrix for AR(1)
   disturbances with autoregressive coefficient @a; @n is the expansion
   factor.
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

/* Regressor matrix: by default we put in constant plus linear trend
   (det = 1) or quadratic trend (det = 2), summed appropriately based
   on @n.

   If the user has supplied high-frequency covariates in @X, we
   compress them from column 2 or 3 onward.

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
	gretl_matrix_set(CX, t, 0, n);
	xt1 = xt2 = 0.0;
	for (i=0; i<n; i++) {
	    xt1 += k;
	    if (det == 2) {
		xt2 += k * k;
	    }
	    k++;
	}
	gretl_matrix_set(CX, t, 1, xt1);
	if (det == 2) {
	    gretl_matrix_set(CX, t, 2, xt2);
	}
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		xt1 = 0.0;
		for (i=0; i<n; i++) {
		    xt1 += gretl_matrix_get(X, s + i, j);
		}
		gretl_matrix_set(CX, t, det+1+j, xt1);
	    }
	    s += n;
	}
    }
}

static void make_Xx_beta (gretl_vector *y, const double *b,
			  const gretl_matrix *X, int det)
{
    int i, j, t;

    for (i=0; i<y->rows; i++) {
	t = i + 1;
	y->val[i] = b[0] + b[1]*t;
	if (det == 2) {
	    y->val[i] += b[2]*t*t;
	}
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		y->val[i] += b[det+1+j] * gretl_matrix_get(X, i, j);
	    }
	}
    }
}

/* first-order autocorrelation of residuals */

static double acf_1 (const double *u, int T)
{
    double num = 0, den = 0;
    int t;

    for (t=0; t<T; t++) {
	den += u[t] * u[t];
	if (t > 0) {
	    num += u[t] * u[t-1];
	}
    }

    if (num < 1.0e-9) {
	return 0;
    }

    return num / den;
}

/**
 * chow_lin_interpolate:
 * @Y: T x k: holds the original data to be expanded.
 * @X: (optionally) holds covariates of Y at the higher frequency:
 * if these are supplied they supplement the default set of
 * regressors, namely, constant plus linear or quadratic trend.
 * @xfac: the expansion factor: 3 for quarterly to monthly
 * or 4 for annual to quarterly. Only these factors are
 * supported.
 * @det: 1 for linear, 2 for quadratic trend.
 * @err: location to receive error code.
 *
 * Interpolate, from annual to quarterly or quarterly to monthly,
 * via the Chow-Lin method. See Gregory C. Chow and An-loh Lin,
 * "Best Linear Unbiased Interpolation, Distribution, and
 * Extrapolation of Time Series by Related Series", The
 * Review of Economics and Statistics, Vol. 53, No. 4
 * (November 1971) pp. 372-375.
 *
 * If @X is given it must have T * @xfac rows.
 *
 * Returns: matrix containing the expanded series, or
 * NULL on failure.
 */

gretl_matrix *chow_lin_interpolate (const gretl_matrix *Y,
				    const gretl_matrix *X,
				    int xfac, int det,
				    int *err)
{
    gretl_matrix_block *B;
    gretl_matrix *CX, *b, *u, *W, *Z;
    gretl_matrix *Tmp1, *Tmp2;
    gretl_matrix *Yx = NULL;
    gretl_matrix *y, *yx;
    gretl_matrix my, myx;
    int nx, ny = Y->cols;
    int T = Y->rows;
    int Tx = T * xfac;
    int i;

    /* Note: checks to the effect that xfac = 3 or 4, and,
       if X is non-NULL, that X->rows = xfac * Y->rows,
       should have already been performed.
    */

    gretl_matrix_init(&my);
    gretl_matrix_init(&myx);

    nx = (det == 2)? 3 : 2;
    if (X != NULL) {
	nx += X->cols;
    }

    /* the return value */
    Yx = gretl_zero_matrix_new(Tx, ny);
    if (Yx == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

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

    /* regressors: constant and linear or quadratic trend,
       plus anything the user has added */
    fill_CX(CX, xfac, det, X);
#if CL_DEBUG
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

	/* initial OLS */
	*err = gretl_matrix_ols(y, CX, b, NULL, u, NULL);

	if (!*err) {
	    a = acf_1(u->val, T);
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
		double bracket[] = {0, 0.99};
		struct chowlin cl = {xfac, a};

		*err = gretl_fzero(bracket, 1.0e-12,
				   chow_lin_callback, &cl,
				   &a, OPT_NONE, NULL);
	    }
	}

	if (!*err) {
	    make_CVC(W, xfac, a);
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

	    gretl_matrix_multiply_by_scalar(yx, xfac);
	}
    }

    gretl_matrix_block_destroy(B);

    return Yx;
}
