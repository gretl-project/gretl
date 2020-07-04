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
#define LIMIT_R_SSR 1
#define SSR_LOGISTIC 1
#define RHOMAX 0.999

/* aggregation types */
enum {
    AGG_SUM, /* sum */
    AGG_AVG, /* average */
    AGG_EOP, /* end of period */
    AGG_SOP  /* start of period */
};

/* rho estimation methods */
enum {
    R_ACF1,  /* "traditional" Chow-Lin */
    R_MLE,   /* GLS + MLE (Bournay-Laroque) */
    R_SSR,   /* use GLS criterion (Barbone et al) */
    R_FIXED, /* pre-specified by user */
};

struct chowlin {
    int n;
    double targ;
};

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
    gretl_matrix *Wcpy;
    gretl_matrix *se;
    int s, det, agg;
    int method;
    int netvcv;
    int verbose;
    double lnl;
    double SSR;
    double s2;
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

/* BFGS callback for ar1_mle(); see Davidson and MacKinnon,
   ETM, pp. 435-6.
*/

static double ar1_loglik (const double *theta, void *data)
{
    struct gls_info *G = data;
    int n = G->y0->rows, k = G->CX->cols;
    double r = theta[0];
    double s = theta[1];
    const double *b = theta + 2;
    const double *y = G->y0->val;
    double onemr2 = 1.0 - r*r;
    double inv2s2 = 1.0 / (2*s*s);
    double ll1 = -0.5*n*LN_2_PI - n*log(s) + 0.5*log(onemr2);
    double u, yf2, Xb1, Xb = 0;
    int i, t;

    /* the first observation */
    for (i=0; i<k; i++) {
	Xb += gretl_matrix_get(G->CX, 0, i) * b[i];
    }
    u = y[0] - Xb;
    yf2 = onemr2 * u * u;

    /* subsequent observations */
    for (t=1; t<n; t++) {
	Xb1 = Xb;
	Xb = 0;
	for (i=0; i<k; i++) {
	    Xb += gretl_matrix_get(G->CX, t, i) * b[i];
	}
	u = y[t] - r*y[t-1] - Xb + r*Xb1;
	yf2 += u * u;
    }

    return ll1 - inv2s2 * yf2;
}

static void show_regression_results (const struct gls_info *G,
				     double a, int gls,
				     PRN *prn)
{
    const char *enames[] = {"OLS", "GLS"};
    const char *dnames[] = {"const", "trend", "trend^2"};
    const char *snames[] = {"rho", "SSR", "lnl"};
    int i, k = G->b->rows;

    if (G->se != NULL) {
	/* not yet the case for OLS */
	int err = 0, p = 0, n = k + 3;
	char tmp[16];
	gretl_matrix *cs = gretl_matrix_alloc(k, 2);
	gretl_matrix *adds = gretl_matrix_alloc(1, 3);
	gretl_array *names = gretl_array_new(GRETL_TYPE_STRINGS, n, &err);
	int df = G->W->rows - G->CX->cols;

	if (cs == NULL || adds == NULL || names == NULL) {
	    return;
	}

	for (i=0; i<k; i++) {
	    gretl_matrix_set(cs, i, 0, G->b->val[i]);
	    gretl_matrix_set(cs, i, 1, G->se->val[i]);
	}
	adds->val[0] = a;
	adds->val[1] = G->SSR;
	adds->val[2] = G->lnl;

	for (i=0; i<n; i++) {
	    if (i < G->det) {
		gretl_array_set_data(names, i, gretl_strdup(dnames[i]));
	    } else if (i < k) {
		sprintf(tmp, "X%d", i - G->det + 1);
		gretl_array_set_data(names, i, gretl_strdup(tmp));
	    } else {
		gretl_array_set_data(names, i, gretl_strdup(snames[p++]));
	    }
	}

	pprintf(prn, "  %s estimates:\n", enames[gls]);
	print_model_from_matrices(cs, adds, names, df, OPT_I, prn);

	gretl_matrix_free(cs);
	gretl_matrix_free(adds);
	gretl_array_destroy(names);
    } else {
	pprintf(prn, "\n%s coefficients:\n", enames[gls]);
	for (i=0; i<k; i++) {
	    if (i < G->det) {
		pprintf(prn, " %-8s", dnames[i]);
	    } else {
		pprintf(prn, " %c%-7d", 'X', i - G->det + 1);
	    }
	    pprintf(prn, "%#g", G->b->val[i]);
	    if (G->se != NULL) {
		pprintf(prn, " (%#g)", G->se->val[i]);
	    }
	    pputc(prn, '\n');
	}
	pprintf(prn, " %-8s%#.8g\n", "rho", a);
	if (!na(G->SSR)) {
	    pprintf(prn, " %-8s%#.8g\n", "SSR", G->SSR);
	}
	if (!na(G->lnl)) {
	    pprintf(prn, " %-8s%#.8g\n", "lnl", G->lnl);
	}
    }
}

static int ar1_mle (struct gls_info *G, double s, double *rho)
{
    double *theta;
    int fc = 0, gc = 0;
    int k = G->CX->cols;
    int i, nt, err;

    nt = k + 2;
    theta = malloc(nt * sizeof *theta);
    if (theta == NULL) {
	return E_ALLOC;
    }

    theta[0] = *rho;
    theta[1] = s;
    for (i=0; i<k; i++) {
	theta[i+2] = G->b->val[i];
    }

    err = BFGS_max(theta, nt, 300, 1.0e-10,
		   &fc, &gc, ar1_loglik, C_LOGLIK,
		   NULL, G, NULL, OPT_NONE, NULL);

    if (err) {
	if (err == E_NOCONV) {
	    /* try taking the final value regardless */
#if CL_DEBUG
	    fprintf(stderr, "ar1_mle: BFGS_max gave E_NOCONV, "
		    "continuing anyway\n");
#endif
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
	*rho = theta[0] > 0.999 ? 0.999 : theta[0];
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

/* Carry out enough of the Chow-Lin GLS calculations to permit
   computation of a figure of merit (loglikelihood or sum of squared
   GLS residuals), conditional on rho.  This is used if we already
   have a "final" rho value greater than zero, or if we're trying to
   find optimal rho via GLS.
*/

static double cl_gls_calc (const double *rho, void *data)
{
    struct gls_info *G = data;
    double a = *rho;
    double ldet, crit;
    double SSR = NADBL;
    int N = G->y0->rows;
    int err = 0;

#if LIMIT_R_SSR && !SSR_LOGISTIC
    if (G->method == R_MLE) {
	a = logistic_cdf(*rho);
    }
#else
    if (G->method == R_MLE || G->method == R_SSR) {
	a = logistic_cdf(*rho);
    }
#endif

    make_VC(G->VC, N, G->s, a, G->agg);
    make_CVC(G->W, G->VC, G->s, G->agg);
    if (G->netvcv) {
	gretl_matrix_multiply_by_scalar(G->W, 1/(1.0-a*a));
    }
    if (G->Wcpy == NULL) {
	G->Wcpy = gretl_matrix_copy(G->W);
    } else {
	gretl_matrix_copy_values(G->Wcpy, G->W);
    }
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

	ldet = gretl_matrix_log_determinant(G->Wcpy, &err);
    }

    if (!err) {
	G->SSR = SSR = gretl_scalar_qform(G->u, G->W, &err);
	if (!err) {
	    G->s2 = SSR / N;
	    G->lnl = -0.5*N - 0.5*N*LN_2_PI - N*log(G->s2)/2 - ldet/2;
	    /* df adjustment, for R-compatibility */
	    G->s2 = SSR / (N - G->CX->cols);
	}
    }

    if (err) {
	crit = G->lnl = NADBL;
    } else {
	crit = G->method == R_SSR ? SSR : G->lnl;
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
    gretlopt opt = OPT_NONE;
    double r = *a > 0 ? *a : 0.5;
    int fc = 0, gc = 0;
    int err;

    if (G->verbose > 1) {
	opt = OPT_V;
    }

    /* this can be made conditional on G->method == R_GSS */
    G->netvcv = 1;

#if LIMIT_R_SSR
    /* prevent R_SSR from pushing @r above RHOMAX */
    if (G->method == R_SSR) {
	gretl_matrix bounds;
# if SSR_LOGISTIC
	double bvals[] = {1, -10, -log(1.0/RHOMAX - 1)};
	double lrho = -log(1/r - 1);
	double *rptr = &lrho;
# else
	double bvals[] = {1, 0, RHOMAX};
	double *rptr = &r;
# endif

	gretl_matrix_init_full(&bounds, 1, 3, bvals);
	err = LBFGS_max(rptr, 1, 200, 1.0e-12,
			&fc, &gc, cl_gls_calc, C_SSR,
			NULL, NULL, G, &bounds, opt, prn);
# if SSR_LOGISTIC
	if (!err) {
	    r = logistic_cdf(lrho);
	}
# endif
    } else {
	/* G->method = R_MLE */
	double lrho = -log(1/r - 1);

	err = BFGS_max(&lrho, 1, 200, 1.0e-14,
		       &fc, &gc, cl_gls_calc, C_LOGLIK,
		       NULL, G, NULL, opt, prn);
	if (!err) {
	    r = logistic_cdf(lrho);
	}
    }
#else
    /* treat R_SSR and R_MLE symmetrically */
    int crit = G->method == R_SSR ? C_SSR : C_LOGLIK;
    gretlopt opt = G->method == R_SSR ? (opt | OPT_I) : opt;
    double lrho = -log(1/r - 1);

    err = BFGS_max(&lrho, 1, 200, 1.0e-12,
		   &fc, &gc, cl_gls_calc, crit,
		   NULL, G, NULL, opt, prn);
    if (!err) {
	r = logistic_cdf(lrho);
    }
#endif

#if CL_DEBUG
    fprintf(stderr, "cl_gls_max: method = %d, err = %d\n",
	    G->method, err);
#endif

    if (!err) {
	/* record standard errors of GLS coeffs */
	int i, k = G->Z->cols;

	if (G->se == NULL) {
	    G->se = gretl_matrix_alloc(k, 1);
	}
	for (i=0; i<k; i++) {
	    G->se->val[i] = sqrt(G->s2 * gretl_matrix_get(G->Z, i, i));
	}
    }

    if (!err) {
	*a = r;
	if (G->netvcv) {
	    /* restore full W for subsequent calculations */
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

static double acf_1 (struct gls_info *G)
{
    double rho, num = 0, den = 0;
    double *u = G->u->val;
    int t, T = G->u->rows;

    for (t=0; t<T; t++) {
	den += u[t] * u[t];
	if (t > 0) {
	    num += u[t] * u[t-1];
	}
    }

    if (num < 1.0e-9) {
	return 0;
    }

    rho = num / den;

#if CL_DEBUG
    fprintf(stderr, "initial rho from OLS residuals: %g\n", rho);
#endif

    /* improve the initial estimate of @rho via ML */
    ar1_mle(G, sqrt(den / T), &rho);

    return rho;
}

/* We come here if (a) the caller has specified a fixed rho
   of zero, or (b) we're estimating low-frequency rho from
   an initial OLS regression and converting to high-frequency.
*/

static int cl_ols (struct gls_info *G,
		   const gretl_matrix *X,
		   gretl_matrix *y,
		   double *rho,
		   PRN *prn)
{
    gretl_matrix *u;
    double a = 0;
    int err;

    u = G->method == R_FIXED ? NULL : G->u;
    err = gretl_matrix_ols(G->y0, G->CX, G->b, NULL, u, NULL);
    if (err) {
	return err;
    }

    if (G->method != R_FIXED) {
	a = acf_1(G);
	if (a <= 0.0) {
	    a = 0; /* don't pursue negative @a */
	} else if (G->agg >= AGG_EOP) {
	    a = pow(a, 1.0/G->s);
	} else if (a >= RHOMAX) {
	    ; /* leave it be! */
	} else {
	    double bracket[] = {0, RHOMAX};
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
	/* no autocorrelation: GLS isn't needed */
	make_X_beta(y, G->b->val, X, G->det);
	if (G->agg == AGG_AVG) {
	    gretl_matrix_multiply_by_scalar(y, G->s);
	}
	if (prn != NULL && G->verbose) {
	    show_regression_results(G, a, 0, prn);
	}
    }

    *rho = a;

    return err;
}

static int fill_chowlin_bundle (struct gls_info *G,
				double a,
				gretl_bundle *b)
{
    gretl_matrix *m;

    gretl_bundle_set_scalar(b, "rho", a);
    gretl_bundle_set_scalar(b, "lnl", G->lnl);
    if (!na(G->SSR)) {
	gretl_bundle_set_scalar(b, "SSR", G->SSR);
    }
    m = gretl_matrix_copy(G->b);
    gretl_bundle_donate_data(b, "coeff", m, GRETL_TYPE_MATRIX, 0);
    if (G->se != NULL) {
	gretl_bundle_donate_data(b, "stderr", G->se, GRETL_TYPE_MATRIX, 0);
	G->se = NULL;
    }

    return 0;
}

/**
 * chow_lin_disagg:
 * @Y0: N x k: holds the original data to be expanded.
 * @X: (optionally) holds covariates of Y at the higher frequency;
 * if these are supplied they supplement the deterministic
 * terms (if any) as signalled by @det.
 * @s: the expansion factor: 3 for quarterly to monthly,
 * 4 for annual to quarterly or 12 for annual to monthly.
 * @agg: aggregation type.
 * @method: method for estimating rho.
 * @det: 0 for none, 1 for constant, 2 for linear trend, 3 for
 * quadratic trend.
 * @rho: fixed rho value, if wanted.
 * @prn: printing struct pointer.
 * @err: location to receive error code.
 *
 * Distribute or interpolate via the method of Chow and Lin. See
 * See Gregory C. Chow and An-loh Lin, "Best Linear Unbiased
 * Interpolation, Distribution, and Extrapolation of Time Series
 * by Related Series", Review of Economics and Statistics, Vol. 53,
 * No. 4 (November 1971) pp. 372-375.
 *
 * If @X is given it must have at least @s * N rows.
 *
 * Returns: matrix containing the expanded series, or
 * NULL on failure.
 */

static gretl_matrix *chow_lin_disagg (const gretl_matrix *Y0,
				      const gretl_matrix *X,
				      int s, int agg, int method,
				      int det, double rho,
				      gretl_bundle *res,
				      int verbose, PRN *prn,
				      int *perr)
{
    struct gls_info G = {0};
    gretl_matrix_block *B;
    gretl_matrix *Y;
    gretl_matrix *y0, *y;
    gretl_matrix my0, my;
    double a = 0;
    int rho_given = 0;
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
	if (method == R_ACF1) {
	    method = R_FIXED;
	}
	/* don't mess with negative rho */
	rho = rho < 0 ? 0 : rho;
	rho_given = 1;
    }

    if (verbose == 0 && gretl_messages_on()) {
	verbose = 1;
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

    /* convenience pointers to columns of Y0, Y */
    y0 = gretl_matrix_init_full(&my0, N, 1, Y0->val);
    y = gretl_matrix_init_full(&my, sN + m, 1, Y->val);

    /* hook up various non-malloc'd things */
    G.y0 = y0;
    G.s = s;
    G.det = det;
    G.agg = agg;
    G.method = method;
    G.netvcv = 0;
    G.verbose = verbose;
    G.lnl = G.SSR = G.s2 = NADBL;
    G.Wcpy = G.se = NULL;

    for (i=0; i<ny && !err; i++) {
	a = rho_given ? rho : 0.0;

	if (i > 0) {
	    /* pick up the current columns for reading and writing */
	    y0->val = Y0->val + i * N;
	    y->val = Y->val + i * (sN + m);
	}

	if (!rho_given || a == 0) {
	    err = cl_ols(&G, X, y, &a, prn);
	    if (a == 0) {
		/* this iteration is handled */
		continue;
	    }
	}

	if (!err) {
	    if (method == R_MLE || method == R_SSR) {
		/* pursue further optimization of @a */
		err = cl_gls_max(&a, &G, prn);
	    } else {
		/* just calculate with the current @a */
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

	    if (prn != NULL && G.verbose) {
		show_regression_results(&G, a, 1, prn);
	    }

	    if (agg == AGG_AVG) {
		gretl_matrix_multiply_by_scalar(y, s);
	    }
	}
    }

    if (!err && ny == 1 && res != NULL) {
	fill_chowlin_bundle(&G, a, res);
    }

    *perr = err;
    gretl_matrix_block_destroy(B);
    gretl_matrix_free(G.Wcpy);
    gretl_matrix_free(G.se);

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

static gretl_matrix *real_tdisagg (const gretl_matrix *Y0,
				   const gretl_matrix *X,
				   int s, int agg, int method,
				   int det, double rho,
				   gretl_bundle *r,
				   int verbose, PRN *prn,
				   int *err)
{
    gretl_matrix *ret = NULL;

    if (method < 3) {
	/* Chow-Lin variants */
	if (det < 0) {
	    det = 1;
	}
	ret = chow_lin_disagg(Y0, X, s, agg, method, det, rho,
			      r, verbose, prn, err);
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

static int tdisagg_get_options (gretl_bundle *b,
				int *pagg, int *pmeth,
				int *pdet, double *pr,
				int *pverb)
{
    double rho = NADBL;
    const char *str;
    int agg = 0;
    int method = 0;
    int det = 1;
    int verbose = 0;
    int err = 0;

    if (gretl_bundle_has_key(b, "agg")) {
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
    if (!err && gretl_bundle_has_key(b, "verbose")) {
	verbose = gretl_bundle_get_int(b, "verbose", NULL);
    }

    if (!err) {
	*pagg = agg;
	*pmeth = method;
	*pdet = det;
	*pr = rho;
	*pverb = verbose;
    }

    return err;
}

gretl_matrix *time_disaggregate (const gretl_matrix *Y0,
				 const gretl_matrix *X,
				 int s, gretl_bundle *b,
				 gretl_bundle *r,
				 PRN *prn, int *err)
{
    int agg = 0, method = 0, det = 1;
    int verbose = 0;
    double rho = NADBL;

    if (b != NULL) {
	*err = tdisagg_get_options(b, &agg, &method, &det,
				   &rho, &verbose);
	if (*err) {
	    return NULL;
	}
    }

    return real_tdisagg(Y0, X, s, agg, method, det, rho,
			r, verbose, prn, err);
}

gretl_matrix *tdisagg_basic (const gretl_matrix *Y0,
			     const gretl_matrix *X,
			     int s, int agg, int *err)
{
    int method = 0, det = 1;
    double rho = NADBL;

    return real_tdisagg(Y0, X, s, agg, method, det, rho,
			NULL, 0, NULL, err);
}
