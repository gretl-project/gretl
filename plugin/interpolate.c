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

struct chowlin {
    int n;
    double targ;
};

/* Callback for BFGS, as we adjust the coefficient a so that the
   theoretically derived ratio of polynomials in a matches the
   empirical first-order autocorrelation of the OLS residuals
   (cl->targ).  
*/

static double chow_lin_callback (const double *pa, void *p)
{
    struct chowlin *cl = (struct chowlin *) p;
    double a = *pa;
    double num, den, val;

    if (cl->n == 3) {
	num = a + 2*a*a + 3*pow(a, 3) + 2*pow(a, 4) + pow(a, 5);
	den = 3 + 4*a + 2*a*a;
    } else {
	/* n = 4 */
	num = a + 2*a*a + 3*pow(a, 3) + 4*pow(a, 4) + 3*pow(a, 5)
	    + 2*pow(a, 6) + pow(a, 7);
	den = 4 + 6*a + 4*a*a + 2*pow(a, 3);
    }

    val = num/den - cl->targ;

    return -val * val;
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

/* Generate CVC' without storing the full C or V matrices.  C is the
   selection matrix that transforms from higher frequency to lower
   frequency by summation; V is the autocovariance matrix for AR(1)
   disturbances with autoregressive coefficient a; n is the expansion
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

/* Regressor matrix: by default we put in constant plus 
   quadratic trend, summed appropriately based on @n.
   If the user has supplied high-frequency covariates,
   in @X, we compress them from column 3 onward.
*/

static void make_CX (gretl_matrix *CX, int n,
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
	    xt2 += k * k;
	    k++;
	}
	gretl_matrix_set(CX, t, 1, xt1);
	gretl_matrix_set(CX, t, 2, xt2);

	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		xt1 = 0.0;
		for (i=0; i<n; i++) {
		    xt1 += gretl_matrix_get(X, s + i, j);
		}
		gretl_matrix_set(CX, t, 3+j, xt1);
	    }
	    s += n;
	}
    }
}

static void make_Xx_beta (gretl_vector *y, const double *b,
			  const gretl_matrix *X)
{
    int i, j, t;

    for (i=0; i<y->rows; i++) {
	t = i + 1;
	y->val[i] = b[0] + b[1]*t + b[2]*t*t;
	if (X != NULL) {
	    for (j=0; j<X->cols; j++) {
		y->val[i] += b[3+j] * gretl_matrix_get(X, i, j);
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

    return num / den;
}

static int make_y_vectors (gretl_matrix **py,
			   gretl_matrix **pyx,
			   int T, int Tx)
{
    gretl_matrix *y = gretl_null_matrix_new();
    gretl_matrix *yx = gretl_null_matrix_new();

    if (y == NULL || yx == NULL) {
	free(y);
	free(yx);
	return E_ALLOC;
    } else {
	y->rows = T;
	yx->rows = Tx;
	y->cols = yx->cols = 1;
	*py = y;
	*pyx = yx;
	return 0;
    }
}

/**
 * chow_lin_interpolate:
 * @Y: T x k: holds the original data to be expanded. 
 * @X: (optionally) holds covariates of Y at the higher frequency:
 * if these are supplied they supplement the default set of
 * regressors, namely, constant plus quadratic trend.
 * @xfac: the expansion factor: 3 for quarterly to monthly
 * or 4 for annual to quarterly. Only these factors are
 * supported.
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
				    int xfac, int *err)
{
    gretl_matrix_block *B;
    gretl_matrix *CX, *b, *u, *W, *Z;
    gretl_matrix *Tmp1, *Tmp2;
    gretl_matrix *Yx = NULL;
    gretl_matrix *y, *yx;
    int nx = 3;
    int ny = Y->cols;
    int T = Y->rows;
    int Tx = T * xfac;
    int i;

    /* Note: checks to the effect that xfac = 3 or 4, and,
       if X is non-NULL, that X->rows = xfac * Y->rows,
       should have already been performed.
    */

    if (X != NULL) {
	nx += X->cols;
    }

    Yx = gretl_matrix_alloc(Tx, ny);
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

    /* regressors: constant and quadratic trend, plus
       anything the user has added */
    make_CX(CX, xfac, X);

    if (ny > 1) {
	/* Y has more than 1 column */
	*err = make_y_vectors(&y, &yx, T, Tx);
	if (*err) {
	    gretl_matrix_free(Yx);
	    gretl_matrix_block_destroy(B);
	    return NULL;
	}
    } else {
	y = (gretl_matrix *) Y; /* don't worry, it's really const */
	yx = Yx;
    }

    for (i=0; i<ny; i++) {
	double a = 0.0;

	if (ny > 1) {
	    /* pick up the current column */
	    y->val = Y->val + i * T;
	    yx->val = Yx->val + i * Tx;
	}

	/* initial OLS */
	*err = gretl_matrix_ols(y, CX, b, NULL, u, NULL);

	if (!*err) {
	    struct chowlin cl;
	    int c1, c2;

	    cl.n = xfac;
	    cl.targ = acf_1(u->val, T);
	    *err = BFGS_max(&a, 1, 50, 1.0e-12, &c1, &c2, 
			    chow_lin_callback, C_OTHER, NULL,
			    &cl, NULL, OPT_NONE, NULL);
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
	    make_Xx_beta(yx, b->val, X);

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

    if (ny > 1) {
	free(y);
	free(yx);
    }

    gretl_matrix_block_destroy(B);
    
    return Yx;
}
