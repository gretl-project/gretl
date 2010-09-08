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

#include "libgretl.h"
#include "matrix_extra.h"
#include "gretl_bfgs.h"

struct chowlin {
    int xfac;
    double targ;
};

/* Callback for BFGS, as we adjust a so that the
   theoretically derived ratio of polynomials in a
   matches the empirical first-order autocorrelation
   of the OLS residuals (cl->targ).
*/

static double chow_lin_callback (const double *pa, void *p)
{
    struct chowlin *cl = (struct chowlin *) p;
    double a = *pa;
    double num, den, val;

    if (cl->xfac == 3) {
	num = a + 2*a*a + 3*pow(a, 3) + 2*pow(a, 4) + pow(a, 5);
	den = 3 + 4*a + 2*a*a;
    } else {
	/* xfac = 4 */
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

/* generate CVC' without storing the full C or V matrices */

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

/* multiply VC' into u and increment yx by the result */

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

/* regressor matrix: constant plus quadratic trend,
   summed appropriately based on @n */

static void make_CX (gretl_matrix *X, int n)
{
    double xt1, xt2;
    int t, i, k = 1;

    for (t=0; t<X->rows; t++) {
	gretl_matrix_set(X, t, 0, n);
	xt1 = xt2 = 0.0;
	for (i=0; i<n; i++) {
	    xt1 += k;
	    xt2 += k * k;
	    k++;
	}
	gretl_matrix_set(X, t, 1, xt1);
	gretl_matrix_set(X, t, 2, xt2);
    }
}

static void make_Xx_beta (gretl_vector *y, const double *b)
{
    int i;

    for (i=1; i<=y->rows; i++) {
	y->val[i-1] = b[0] + b[1] * i + b[2] * i * i;
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

/*
   Interpolate, from annual to quarterly or quaterly to monthly,
   via Chow-Lin method. See Gregory C. Chow and An-loh Lin,
   "Best Linear Unbiased Interpolation, Distribution, and 
   Extrapolation of Time Series by Related Series",
   The Review of Economics and Statistics, Vol. 53, No. 4 
   (Nov., 1971) pp. 372-375.

   In this implementation the only regressors used are a
   constant and quadratic trend.

   @y holds the original data to be expanded.
   @xfac is the expansion factor: 3 for quarterly to monthly
   or 4 for annual to quarterly. Only these factors are
   supported.
*/

gretl_matrix *chow_lin_interpolate (const gretl_matrix *y, 
				    int xfac, int *err)
{
    gretl_matrix_block *B;
    gretl_matrix *X, *b, *u, *W, *Z;
    gretl_matrix *Tmp1, *Tmp2;
    gretl_matrix *yx = NULL;
    double a = 0.5;
    int T = y->rows;

    if (xfac != 3 && xfac != 4) {
	*err = E_DATA;
	return NULL;
    }	

    yx = gretl_column_vector_alloc(T * xfac);
    if (yx == NULL) {
	*err = E_ALLOC;
	return NULL;
    }	

    B = gretl_matrix_block_new(&X, T, 3,
			       &W, T, T,
			       &b, 3, 1,
			       &u, T, 1,
			       &Z, 3, 3,
			       &Tmp1, 3, T,
			       &Tmp2, 3, T,
			       NULL);
    if (B == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(yx);
	return NULL;
    }

    /* regressors: constant and quadratic trend */
    make_CX(X, xfac);

    /* initial OLS */
    *err = gretl_matrix_ols(y, X, b, NULL, u, NULL);

    if (!*err) {
	struct chowlin cl;
	int c1, c2;

	cl.xfac = xfac;
	/* first-order autocorrelation of residuals */
	cl.targ = acf_1(u->val, T);
	/* solve for a */
	*err = BFGS_max(&a, 1, 50, 1.0e-12, &c1, &c2, 
			chow_lin_callback, C_OTHER, NULL,
			&cl, NULL, OPT_NONE, NULL);
    }

    if (!*err) {
	make_CVC(W, xfac, a);
	*err = gretl_invert_symmetric_matrix(W);
    }

    if (!*err) {
	gretl_matrix_qform(X, GRETL_MOD_TRANSPOSE,
			   W, Z, GRETL_MOD_NONE);
	*err = gretl_invert_symmetric_matrix(Z);
    }  

    if (!*err) {
	/* GLS \hat{\beta} */
	gretl_matrix_multiply_mod(Z, GRETL_MOD_NONE,
				  X, GRETL_MOD_TRANSPOSE,
				  Tmp1, GRETL_MOD_NONE);
	gretl_matrix_multiply(Tmp1, W, Tmp2);
	gretl_matrix_multiply(Tmp2, y, b);

	/* Xx \hat{\beta} */
	make_Xx_beta(yx, b->val);

	/* GLS residuals */
	gretl_matrix_copy_values(u, y);
	gretl_matrix_multiply_mod(X, GRETL_MOD_NONE,
				  b, GRETL_MOD_NONE,
				  u, GRETL_MOD_DECREMENT);

	/* yx = Xx*beta + V*C'*W*u */
	gretl_matrix_reuse(Tmp1, T, 1);
	gretl_matrix_multiply(W, u, Tmp1);
	mult_VC(yx, Tmp1, xfac, a);
	gretl_matrix_multiply_by_scalar(yx, xfac);
    }

    gretl_matrix_block_destroy(B);
    
    return yx;
}



