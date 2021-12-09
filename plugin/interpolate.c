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
#define SSR_LOGISTIC 1
#define RHOMAX 0.999

/* aggregation types */
enum {
    AGG_AVG, /* average */
    AGG_SUM, /* sum */
    AGG_EOP, /* end of period */
    AGG_SOP  /* start of period */
};

/* rho estimation methods */
enum {
    R_ACF1,  /* "traditional" Chow-Lin */
    R_MLE,   /* GLS + MLE (Bournay-Laroque) */
    R_SSR,   /* use GLS criterion (Barbone et al) */
    R_UROOT, /* Fernandez */
    R_FIXED  /* pre-specified by user */
};

typedef enum {
    CL_SAVE_SE = 1 << 0,
    CL_NETVCV  = 1 << 1,
    CL_VERBOSE = 1 << 2,
    CL_SHOWMAX = 1 << 3,
    CL_TRUNC   = 1 << 4
} CLflags;

struct chowlin {
    int n;
    double targ;
};

struct gls_info {
    const gretl_matrix *y0;
    const gretl_matrix *X;
    const char *yname;
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
    CLflags flags;
    double lnl;
    double SSR;
    double s2;
};

static const char *aggtype_names[] = {
    "average", "sum", "last", "first"
};

static const char *method_names[] = {
    "chow-lin", "chow-lin-mle", "chow-lin-ssr",
    "fernandez", "denton-pfd", "denton-afd", NULL
};

static int td_plot (const gretl_matrix *y0,
		    const gretl_matrix *y,
		    int s, int agg, int method,
		    DATASET *dset);

static int force_iteration (int m)
{
    return m == R_MLE || m == R_SSR;
}

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

static void maybe_print_depvar (const struct gls_info *G,
				PRN *prn)
{
    if (G->yname != NULL && prn != NULL) {
	pprintf(prn, "%s: %s\n", _("Dependent variable"), G->yname);
    }
}

static void show_regression_results (const struct gls_info *G,
				     double a, int gls,
				     PRN *prn)
{
    const char *dnames[] = {"const", "trend", "trend^2"};
    const char *snames[] = {"rho", "SSR", "lnl"};
    const char **Xnames = NULL;
    char tmp[16];
    int i, k = G->b->rows;
    int df = G->CX->rows - G->CX->cols;
    int T = G->CX->rows;
    int err = 0, p = 0, n = k + 3;

    gretl_matrix *cs = gretl_matrix_alloc(k, 2);
    gretl_matrix *adds = gretl_matrix_alloc(1, 3);
    gretl_array *names = gretl_array_new(GRETL_TYPE_STRINGS, n, &err);

    if (cs == NULL || adds == NULL || names == NULL) {
	return;
    }

    if (G->X != NULL) {
	Xnames = gretl_matrix_get_colnames(G->X);
    }

    for (i=0; i<k; i++) {
	gretl_matrix_set(cs, i, 0, G->b->val[i]);
	if (G->se != NULL) {
	    gretl_matrix_set(cs, i, 1, G->se->val[i]);
	} else {
	    gretl_matrix_set(cs, i, 1, NADBL);
	}
    }
    adds->val[0] = a;
    adds->val[1] = G->SSR;
    adds->val[2] = G->lnl;

    for (i=0; i<n; i++) {
	if (i < G->det) {
	    gretl_array_set_data(names, i, gretl_strdup(dnames[i]));
	} else if (i < k) {
	    if (Xnames != NULL) {
		gretl_array_set_data(names, i, gretl_strdup(Xnames[i-G->det]));
	    } else {
		sprintf(tmp, "X%d", i - G->det + 1);
		gretl_array_set_data(names, i, gretl_strdup(tmp));
	    }
	} else {
	    gretl_array_set_data(names, i, gretl_strdup(snames[p++]));
	}
    }

    if (G->method == R_UROOT) {
	pprintf(prn, "%s", _("GLS estimates"));
	pprintf(prn, " (fernandez) T = %d\n", T);
	maybe_print_depvar(G, prn);
    } else if (G->method == R_MLE || G->method == R_SSR) {
	pprintf(prn, "%s", _("Iterated GLS estimates"));
	pprintf(prn, " (%s) T = %d\n", method_names[G->method], T);
	maybe_print_depvar(G, prn);
	if (G->flags & CL_TRUNC) {
	    pprintf(prn, "%s\n", _("rho truncated to zero"));
	}
    } else if (a == 0) {
	pprintf(prn, "%s T = %d\n", _("OLS estimates"), T);
	maybe_print_depvar(G, prn);
    } else {
	pprintf(prn, "%s", _("GLS estimates"));
	pprintf(prn, " (%s) T = %d\n", (G->method == R_FIXED)?
		"fixed rho" : "chow-lin", T);
	maybe_print_depvar(G, prn);
    }

    print_model_from_matrices(cs, adds, names, df, OPT_I, prn);

    gretl_matrix_free(cs);
    gretl_matrix_free(adds);
    gretl_array_destroy(names);
}

static int ar1_mle (struct gls_info *G, double s, double *rho,
		    PRN *prn)
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
	*rho = theta[0] > 0.999 ? 0.999 : theta[0];
	if (G->flags & CL_SHOWMAX) {
	    pprintf(prn, "rho as revised via ML: %g\n", *rho);
	}
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

/* The counterpart of make_VC() above, when using Fernandez
   or when rho = 0.
*/

static void make_alt_VC (gretl_matrix *VC, int s, int agg,
			 int method)
{
    double vij, *vj, *cval;
    int N = VC->cols;
    int rows = VC->rows;
    size_t colsz = rows * sizeof *vj;
    int i, j, k;

    vj = malloc(colsz);
    gretl_matrix_zero(VC);

    /* (1) construct C' */
    k = (agg == AGG_EOP)? s-1 : 0;
    for (j=0; j<N; j++) {
	if (agg >= AGG_EOP) {
	    gretl_matrix_set(VC, k, j, 1);
	} else {
	    for (i=0; i<s; i++) {
		gretl_matrix_set(VC, i+k, j, 1);
	    }
	}
	k += s;
    }

    if (method == R_UROOT) {
	/* (2) premultiply by inv(D'D) */
	for (k=0; k<2; k++) {
	    cval = VC->val;
	    for (j=0; j<N; j++) {
		memcpy(vj, cval, colsz);
		vij = vj[rows-1];
		for (i=0; i<rows; i++) {
		    cval[i] = vij;
		    if (i < rows-1) {
			vij += vj[rows-i-2];
		    }
		}
		cval += rows;
	    }
	}
    }

    free(vj);
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

/* version of make_EVC() for use with fernandez */

static void make_EDDC (gretl_matrix *EVC,
		       gretl_matrix *VC)
{
    double vij;
    int i, j;

    for (j=0; j<EVC->cols; j++) {
	vij = gretl_matrix_get(VC, VC->rows-1, j);
	for (i=0; i<EVC->rows; i++) {
	    gretl_matrix_set(EVC, i, j, vij);
	}
    }
}

/* Multiply VC' into W*u and increment y by the result */

static void multiply_by_VC (gretl_matrix *y,
			    struct gls_info *G,
			    int m, double a)
{
    gretl_matrix *Wu = G->Tmp1;
    gretl_matrix *EVC = NULL;
    gretl_matrix ext = {0};
    int i, sN = y->rows;

    if (m > 0) {
	/* we need a different covariance matrix for
	   extrapolation
	*/
	sN -= m;
	gretl_matrix_reuse(y, sN, 1);
	gretl_matrix_init(&ext);
	ext.rows = m;
	ext.cols = 1;
	ext.val = y->val + sN;
	EVC = gretl_matrix_alloc(m, sN / G->s);
    }

    gretl_matrix_multiply_mod(G->VC, GRETL_MOD_NONE,
			      Wu, GRETL_MOD_NONE,
			      y, GRETL_MOD_CUMULATE);

    if (EVC != NULL) {
	if (G->method == R_UROOT) {
	    make_EDDC(EVC, G->VC);
	} else {
	    make_EVC(EVC, G->s, a, G->agg);
	}
	gretl_matrix_multiply_mod(EVC, GRETL_MOD_NONE,
				  Wu, GRETL_MOD_NONE,
				  &ext, GRETL_MOD_CUMULATE);
	gretl_matrix_reuse(y, sN+m, 1);
	gretl_matrix_free(EVC);
    } else if (m > 0) {
	/* fallback */
	for (i=0; i<m; i++) {
	    y->val[sN+i] = NADBL;
	}
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

    if (G->method == R_UROOT) {
	make_alt_VC(G->VC, G->s, G->agg, G->method);
    } else {
	if (!(G->flags & CL_TRUNC)) {
	    /* if we're on a second pass after truncation
	       to zero, @rho will not be given in
	       logistic form
	    */
#if SSR_LOGISTIC
	    if (G->method == R_MLE || G->method == R_SSR) {
		a = logistic_cdf(*rho);
	    }
#else
	    if (G->method == R_MLE) {
		a = logistic_cdf(*rho);
	    }
#endif
	}
	make_VC(G->VC, N, G->s, a, G->agg);
    }

    make_CVC(G->W, G->VC, G->s, G->agg);
    if ((G->flags & CL_NETVCV) && a > 0) {
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

/* record standard errors of GLS coeffs */

static void add_gls_se (struct gls_info *G)
{
    int i, k = G->Z->cols;

    if (G->se == NULL) {
	G->se = gretl_matrix_alloc(k, 1);
    }
    for (i=0; i<k; i++) {
	G->se->val[i] = sqrt(G->s2 * gretl_matrix_get(G->Z, i, i));
    }
}

static int cl_gls_one (double *a, struct gls_info *G)
{
    double crit = cl_gls_calc(a, G);
    int err = 0;

    if (na(crit)) {
	err = E_NOCONV; /* FIXME? */
    } else if (G->flags & CL_SAVE_SE) {
	add_gls_se(G);
    }

    if (G->flags & CL_NETVCV) {
	/* restore full W for subsequent calculations */
	double r = *a;
	gretl_matrix_multiply_by_scalar(G->W, 1/(1-r*r));
    }

    return err;
}

/* Driver for optimization of rho via GLS, either on the criterion
   of likelihood or the straight GLS criterion.  In line with other
   implementations we place an upper bound on rho short of 1.0.
*/

static int cl_gls_max (double *a, struct gls_info *G, PRN *prn)
{
    gretlopt opt = OPT_NONE;
    double r = *a > 0 ? *a : 0.5;
    int fc = 0, gc = 0;
    int err;

    if (G->flags & CL_SHOWMAX) {
	opt = OPT_V;
    }

    if (G->method == R_SSR) {
	/* prevent R_SSR from pushing @r above RHOMAX */
	gretl_matrix bounds;
#if SSR_LOGISTIC
	double bvals[] = {1, -20, -log(1.0/RHOMAX - 1)};
	double lrho = -log(1/r - 1);
	double *rptr = &lrho;
#else
	double bvals[] = {1, 0, RHOMAX};
	double *rptr = &r;
#endif

	gretl_matrix_init_full(&bounds, 1, 3, bvals);
	err = LBFGS_max(rptr, 1, 200, 1.0e-12,
			&fc, &gc, cl_gls_calc, C_SSR,
			NULL, NULL, G, &bounds, opt, prn);
#if SSR_LOGISTIC
	if (!err) {
	    r = logistic_cdf(lrho);
	}
#endif
    } else {
	/* G->method = R_MLE */
	double lrho = -log(1/r - 1);

	err = BFGS_max(&lrho, 1, 200, 1.0e-14,
		       &fc, &gc, cl_gls_calc, C_LOGLIK,
		       NULL, G, NULL, opt, prn);
	if (!err) {
	    r = logistic_cdf(lrho);
	} else if (err == E_NOCONV && lrho < -6) {
	    /* looks like we're butting up against 0 */
	    fprintf(stderr, "R_MLE: stopped at lrho = %g\n", lrho);
	    err = 0;
	    r = 0;
	    /* recalculate with rho = 0 */
	    G->flags |= CL_TRUNC;
	    cl_gls_calc(&r, G);
	}
    }

#if CL_DEBUG
    fprintf(stderr, "cl_gls_max: method = %s, err = %d\n",
	    method_names[G->method], err);
#endif

    if (!err && (G->flags & CL_SAVE_SE)) {
	add_gls_se(G);
    }

    if (!err) {
	*a = r;
	if ((G->flags & CL_NETVCV) && r > 0) {
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

static double acf_1 (struct gls_info *G, PRN *prn)
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

    rho = num / den;

    if (G->flags & CL_SHOWMAX) {
	pprintf(prn, "Initial rho from OLS residuals: %g\n", rho);
    }

    if (rho < 1.0e-6) {
	/* truncate */
	rho = 0;
    } else {
	/* improve the initial estimate of @rho via ML */
	ar1_mle(G, sqrt(den / T), &rho, prn);
    }

    return rho;
}

static int prepare_OLS_stats (struct gls_info *G)
{
    int k = G->CX->cols;
    int n = G->CX->rows;
    int err = 0;

    G->se = gretl_matrix_alloc(k, 1);
    if (G->se == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(G->CX, GRETL_MOD_TRANSPOSE,
			      G->CX, GRETL_MOD_NONE,
			      G->Z, GRETL_MOD_NONE);
    err = gretl_invert_symmetric_matrix(G->Z);

    if (!err) {
	const double *u = G->u->val;
	double s2, zii, n2 = 0.5*n;
	int i;

	G->SSR = 0.0;
	for (i=0; i<n; i++) {
	    G->SSR += u[i] * u[i];
	}
	s2 = G->SSR / (n - k);
	for (i=0; i<k; i++) {
	    zii = gretl_matrix_get(G->Z, i, i);
	    G->se->val[i] = sqrt(s2 * zii);
	}
	G->lnl = -n2 * (1 + LN_2_PI - log(n)) - n2 * log(G->SSR);
	if (G->agg == AGG_SUM) {
	    G->SSR /= G->s;
	}
    }

    return err;
}

/* We come here if (a) the caller has specified a fixed rho
   of zero, or (b) we're estimating low-frequency rho from
   an initial OLS regression and converting to high-frequency.
*/

static int cl_ols (struct gls_info *G,
		   const gretl_matrix *X,
		   gretl_matrix *y,
		   int m, double *rho,
		   PRN *prn)
{
    double a = 0;
    int err;

    err = gretl_matrix_ols(G->y0, G->CX, G->b, NULL, G->u, NULL);
    if (err) {
#if CL_DEBUG
	fprintf(stderr, "cl_ols: error %d from initial regression\n", err);
#endif
	return err;
    }

    if (G->method != R_FIXED) {
	a = acf_1(G, prn);
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

    if (a == 0 && !force_iteration(G->method)) {
	/* no autocorrelation: GLS isn't needed */
	int N = G->u->rows;
	int nx = G->b->rows;

	/* y = X\hat{\beta}_{OLS} */
	make_X_beta(y, G->b->val, X, G->det);

	/* add C'(CC')^{-1} * \hat{u}_{OLS} */
	make_alt_VC(G->VC, G->s, G->agg, R_ACF1);
	make_CVC(G->W, G->VC, G->s, G->agg);
	gretl_invert_symmetric_matrix(G->W);
	gretl_matrix_reuse(G->Tmp1, N, 1);
	gretl_matrix_multiply(G->W, G->u, G->Tmp1);
	multiply_by_VC(y, G, m, 0);
	gretl_matrix_reuse(G->Tmp1, nx, N);

	if (G->agg == AGG_AVG) {
	    gretl_matrix_multiply_by_scalar(y, G->s);
	}
	if (G->flags & CL_SAVE_SE) {
	    prepare_OLS_stats(G);
	}
	if (G->method == R_MLE || G->method == R_SSR) {
	    /* don't lie about the method we actually used */
	    G->method = R_ACF1;
	}
	if (prn != NULL && (G->flags & CL_VERBOSE)) {
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

    if (G->method == R_FIXED) {
	gretl_bundle_set_string(b, "method", "fixed rho");
    } else {
	gretl_bundle_set_string(b, "method",
				method_names[G->method]);
    }
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
 * @verbose: 0, 1 or 2.
 * @plot: currently just a boolean.
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
				      int verbose, int plot,
				      DATASET *dset,
				      PRN *prn, int *perr)
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

#if CL_DEBUG
    fprintf(stderr, "chow_lin_disagg; N=%d, s=%d, m=%d\n", N, s, m);
    fprintf(stderr, "Y0 is %d x %d\n",Y0->rows, Y0->cols);
    if (X != NULL) {
	fprintf(stderr, "X is %d x %d\n", X->rows, X->cols);
    }
#endif

    /* pointer to get varnames, if present */
    G.X = X;

    /* get (single) yname if present */
    if (Y0->cols == 1) {
	const char **S = gretl_matrix_get_colnames(Y0);

	if (S != NULL) {
	    G.yname = S[0];
	}
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
    G.lnl = G.SSR = G.s2 = NADBL;
    G.Wcpy = G.se = NULL;

    /* set some state flags */
    if (verbose) {
	G.flags |= (CL_VERBOSE | CL_SAVE_SE);
	if (verbose > 1) {
	    G.flags |= CL_SHOWMAX;
	}
    } else if (res != NULL) {
	G.flags |= CL_SAVE_SE;
    }
    if (method != R_UROOT) {
	G.flags |= CL_NETVCV;
    }

    for (i=0; i<ny && !err; i++) {
	a = rho_given ? rho : 0.0;

	if (i > 0) {
	    /* pick up the current columns for reading and writing */
	    y0->val = Y0->val + i * N;
	    y->val = Y->val + i * (sN + m);
	}

	if (method != R_UROOT && (a == 0 || !rho_given)) {
	    err = cl_ols(&G, X, y, m, &a, prn);
	    if (a == 0 && !force_iteration(method)) {
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
		err = cl_gls_one(&a, &G);
		/* adjust @a for reporting when doing Fernández */
		if (method == R_UROOT) {
		    a = 1;
		}
	    }
	}

	if (!err) {
	    /* final series: y = X*beta + V*C'*W*u */
	    make_X_beta(y, G.b->val, X, det);
	    gretl_matrix_reuse(G.Tmp1, N, 1);
	    gretl_matrix_multiply(G.W, G.u, G.Tmp1);
	    multiply_by_VC(y, &G, m, a);
	    gretl_matrix_reuse(G.Tmp1, nx, N);

	    if (prn != NULL && verbose) {
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
    if (!err && ny == 1 && plot) {
	int meth = method == R_FIXED ? R_ACF1 : method;

	td_plot(Y0, Y, s, agg, meth, dset);
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
   Vol. 66, No. 333 (March 1971), pp. 99-102, first difference
   variant, as modified by P. A. Cholette, "Adjusting Sub-annual
   Series to Yearly Benchmarks," Survey Methodology, Vol. 10, 1984,
   pp. 35-49.

   In the proportional case the solution method is based on Tommaso Di
   Fonzo and Marco Marini, "On the Extrapolation with the Denton
   Proportional Benchmarking Method", IMF Working Paper WP/12/169,
   2012. In the additive case, we employ the "soluzione generale" on
   page 3 of Di Fonzi's 2003 working paper, "Benchmarking de serie
   storiche economiche. Nota tecnica ed estenioni".
*/

static gretl_matrix *denton_fd (const gretl_matrix *Y0,
				const gretl_vector *p,
				int s, int agg,
				int afd, int plot,
				DATASET *dset,
				int *err)
{
    gretl_matrix *M;
    gretl_matrix *Y;
    gretl_matrix *tmp;
    gretl_matrix *DDp = NULL;
    gretl_matrix *y0, *y;
    gretl_matrix my0, my;
    double pk, mij;
    int N = Y0->rows;
    int ny = Y0->cols;
    int sN = s * N;
    int m = p->rows - sN;
    int sNm = sN + m;
    int dim = sNm + N;
    int i, j, ii, k = 0;
    int offset;

    /* we need one big matrix, @M */
    M = gretl_zero_matrix_new(dim, dim);
    /* plus some workspace */
    if (afd) {
	tmp = gretl_matrix_alloc(dim, 1);
    } else {
	tmp = gretl_matrix_alloc(sNm, N);
    }
    /* and the return matrix */
    Y = gretl_matrix_alloc(sNm, ny);

    if (M == NULL || tmp == NULL || Y == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (afd) {
	DDp = gretl_zero_matrix_new(dim, 1);
	if (DDp == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
    }

    /* convenience pointers to columns of Y0, Y */
    y0 = gretl_matrix_init_full(&my0, N, 1, Y0->val);
    y = gretl_matrix_init_full(&my, sNm, 1, Y->val);

    /* In @M, create (D'D ~ diag(p)*J') | (J*diag(p) ~ 0);
       see di Fonzo and Marini, equation (4), or in the
       afd case, (D'D ~ J') | (J ~ 0).
    */

    /* the upper left portion, D'D */
    for (i=0; i<sNm; i++) {
	gretl_matrix_set(M, i, i, (i == 0 || i == sNm-1)? 1 : 2);
	if (i > 0) {
	    gretl_matrix_set(M, i, i-1, -1);
	}
	if (i < sNm-1) {
	    gretl_matrix_set(M, i, i+1, -1);
	}
    }

    /* the bottom and right portions, using @p (or not) */
    k = offset = (agg == AGG_EOP)? s-1 : 0;
    for (i=sNm; i<dim; i++) {
	if (agg >= AGG_EOP) {
	    pk = afd ? 1 : p->val[k];
	    gretl_matrix_set(M, i, offset, pk);
	    gretl_matrix_set(M, offset, i, pk);
	    k += s;
	} else {
	    for (j=offset; j<offset+s; j++) {
		pk = afd ? 1 : p->val[k];
		gretl_matrix_set(M, i, j, pk);
		gretl_matrix_set(M, j, i, pk);
		k++;
	    }
	}
	offset += s;
    }

    *err = gretl_invert_symmetric_indef_matrix(M);

    for (ii=0; ii<ny && !*err; ii++) {
	if (ii > 0) {
	    /* pick up the current columns for reading and writing */
	    y0->val = Y0->val + ii * N;
	    y->val = Y->val + ii * sNm;
	}
	if (afd) {
	    /* form (D'Dp | y0) */
	    DDp->val[0] = p->val[0] - p->val[1];
	    for (i=1; i<sNm-1; i++) {
		DDp->val[i] = 2 * p->val[i] - p->val[i-1] - p->val[i+1];
	    }
	    DDp->val[sNm-1] = p->val[sNm-1] - p->val[sNm-2];
	    memcpy(DDp->val + sNm, y0->val, N * sizeof(double));
	    /* multiply M into (D'Dp | y0) */
	    gretl_matrix_multiply(M, DDp, tmp);
	    /* and extract the portion we want */
	    memcpy(y->val, tmp->val, sNm * sizeof(double));
	} else {
	    /* extract the relevant portion of M-inverse and
	       premultiply by (diag(p) ~ 0) | (0 ~ I)
	    */
	    for (j=0; j<N; j++) {
		for (i=0; i<sNm; i++) {
		    mij = gretl_matrix_get(M, i, j+sNm);
		    gretl_matrix_set(tmp, i, j, mij * p->val[i]);
		}
	    }
	    gretl_matrix_multiply(tmp, y0, y);
	}

	if (agg == AGG_AVG) {
	    gretl_matrix_multiply_by_scalar(y, s);
	}
    }

 bailout:

    gretl_matrix_free(M);
    gretl_matrix_free(tmp);
    gretl_matrix_free(DDp);

    if (*err) {
	gretl_matrix_free(Y);
	Y = NULL;
    } else if (ny == 1 && plot) {
	td_plot(Y0, Y, s, agg, afd ? 5 : 4, dset);
    }

    return Y;
}

static gretl_matrix *real_tdisagg (const gretl_matrix *Y0,
				   const gretl_matrix *X,
				   int s, int agg, int method,
				   int det, double rho,
				   gretl_bundle *r,
				   int verbose, int plot,
				   DATASET *dset,
				   PRN *prn, int *err)
{
    gretl_matrix *ret = NULL;

    if (method < 4) {
	/* Chow-Lin variants */
	if (det < 0) {
	    det = 1;
	}
	ret = chow_lin_disagg(Y0, X, s, agg, method, det, rho,
			      r, verbose, plot, dset, prn, err);
    } else if (method < 6) {
	/* Modified Denton, first differences */
	int ylen = Y0->rows;
	gretl_matrix *X0 = NULL;

	if (X == NULL) {
	    /* as per R, tempdisagg, use a constant if
	       X is not provided */
	    X0 = gretl_unit_matrix_new(s * ylen, 1);
	    X = X0;
	} else {
	    int xlen = gretl_vector_get_length(X);

	    if (ylen == 0 || xlen == 0 || xlen < s * ylen) {
		*err = E_INVARG;
	    }
	}
	if (!*err) {
	    int afd = (method == 5);

	    ret = denton_fd(Y0, X, s, agg, afd, plot, dset, err);
	    gretl_matrix_free(X0);
	}
    } else {
	/* no other choices at present */
	*err = E_INVARG;
    }

    if (ret != NULL && X != NULL) {
	if (ret->rows == X->rows && gretl_matrix_is_dated(X)) {
	    gretl_matrix_transcribe_obs_info(ret, X);
	}
    } else if (ret != NULL && gretl_matrix_get_t2(Y0) > 0) {
	/* experimental */
	int t1 = gretl_matrix_get_t1(Y0);
	
	gretl_matrix_set_t1(ret, t1);
	gretl_matrix_set_t2(ret, t1 + ret->rows - 1);
    }

    return ret;
}

static int get_aggregation_type (const char *s, int *err)
{
    if (!strcmp(s, "avg")) {
	return AGG_AVG;
    } else if (!strcmp(s, "sum")) {
	return AGG_SUM;
    } else if (!strcmp(s, "last")) {
	return AGG_EOP;
    } else if (!strcmp(s, "first")) {
	return AGG_SOP;
    } else {
	*err = E_INVARG;
	return -1;
    }
}

static int get_tdisagg_method (const char *s, int *err)
{
    int i;

    for (i=0; method_names[i] != NULL; i++) {
	if (!strcmp(s, method_names[i])) {
	    return i;
	}
    }

    /* backward compat for the moment */
    if (!strcmp(s, "denton")) {
	return R_UROOT + 1;
    } else if (!strcmp(s, "denton-additive")) {
	return R_UROOT + 2;
    }

    *err = E_INVARG;
    return -1;
}

static int tdisagg_get_options (gretl_bundle *b,
				int *pagg, int *pmeth,
				int *pdet, double *pr,
				int *pverb, int *pplot,
				PRN *prn)
{
    double rho = NADBL;
    const char *str;
    int agg = AGG_SUM; /* debatable */
    int method = 0;
    int det = 1;
    int verbose = 0;
    int plot = 0;
    int err = 0;

    if (gretl_bundle_has_key(b, "aggtype")) {
	str = gretl_bundle_get_string(b, "aggtype", &err);
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
    if (!err && gretl_bundle_has_key(b, "plot")) {
	plot = gretl_bundle_get_int(b, "plot", NULL);
    }

    if (!err && verbose < 0) {
	verbose = gretl_messages_on();
    }

    if (!err) {
	*pagg = agg;
	*pmeth = method;
	*pdet = det;
	*pr = rho;
	*pverb = verbose;
	*pplot = plot;
	if (verbose && method <= R_UROOT) {
	    pprintf(prn, "Aggregation type %s\n", aggtype_names[agg]);
	    if (!na(rho)) {
		pprintf(prn, "Input rho value %g\n", rho);
	    }
	}
    }

    return err;
}

static int td_matrix_check (const gretl_matrix *Y,
			    const gretl_matrix *X)
{
    int i, n;

    if (Y->is_complex) {
	return E_CMPLX;
    } else if (gretl_is_null_matrix(Y)) {
	return E_INVARG;
    } else if (X != NULL && X->is_complex) {
	return E_CMPLX;
    } else if (X != NULL && X->rows * X->cols == 0) {
	return E_INVARG;
    }

    n = Y->rows * Y->cols;
    for (i=0; i<n; i++) {
	if (na(Y->val[i])) {
	    return E_MISSDATA;
	}
    }

    if (X != NULL) {
	n = X->rows * X->cols;
	for (i=0; i<n; i++) {
	    if (na(X->val[i])) {
		return E_MISSDATA;
	    }
	}
    }

    return 0;
}

gretl_matrix *time_disaggregate (const gretl_matrix *Y0,
				 const gretl_matrix *X,
				 int s, gretl_bundle *b,
				 gretl_bundle *res,
				 DATASET *dset,
				 PRN *prn, int *err)
{
    int agg = 0, method = 0, det = 1;
    int verbose = 0, plot = 0;
    double rho = NADBL;

    *err = td_matrix_check(Y0, X);

    if (!*err && b != NULL) {
	*err = tdisagg_get_options(b, &agg, &method, &det,
				   &rho, &verbose, &plot,
				   prn);
    }

    if (*err) {
	return NULL;
    }

#if CL_DEBUG
    fprintf(stderr, "time_disaggregate: method = %s\n",
	    method_names[method]);
#endif

    return real_tdisagg(Y0, X, s, agg, method, det, rho,
			res, verbose, plot, dset, prn, err);
}

/* Add basic info to dummy dataset @hf, sufficient to get
   a time-series x-axis in td plot, if possible.
*/

static int set_hf_data_info (DATASET *hf, DATASET *lf, int s)
{
    int ok = 0;

    if (lf->pd == 1) {
	if (s == 4 || s == 12) {
	    /* annual to quarterly or monthly */
	    hf->pd = s;
	    hf->sd0 = lf->sd0 + (s == 4 ? 0.1 : 0.01);
	    ok = 1;
	}
    } else if (lf->pd == 4) {
	if (s == 3) {
	    /* quarterly to monthly */
	    hf->pd = 12;
	    hf->sd0 = lf->sd0 - 0.1 + 0.01;
	    ok = 1;
	}
    }

    if (ok) {
	hf->structure = TIME_SERIES;
	hf->t1 = lf->t1 * s;
	hf->t2 = lf->t2 * s;
	hf->n = lf->n * s;
    }

    return ok;
}

/* This plot -- which is invoked by giving a non-zero value
   under the "plot" key in the tdisagg options bundle -- is
   intended as a simple sanity check. It's a time series plot
   of the original low-frequency series (with repetition, and
   shown in "step" form) and the "final series" (scaled for
   comparability with the original when aggregation is by
   summation). The plot should offer the right-click "Zoom"
   option so one can inspect it in more detail.
*/

static int td_plot (const gretl_matrix *y0,
		    const gretl_matrix *y,
		    int s, int agg, int method,
		    DATASET *dset)
{
    DATASET hfd = {0};
    gretl_matrix *YY;
    gchar *title;
    int i, k, T = y0->rows;
    int sT = s * T;
    int t, sTm = y->rows;
    int mult = 1;
    int err = 0;

    YY = gretl_matrix_alloc(sTm, 2);
    if (YY == NULL) {
	return E_ALLOC;
    }

    /* original data in first column */
    for (i=0, t=0; i<T; i++) {
	for (k=0; k<s; k++) {
	    gretl_matrix_set(YY, t++, 0, y0->val[i]);
	}
    }

    if (agg == AGG_SUM) {
	mult = s;
    }

    /* final series in second column */
    for (t=0; t<sTm; t++) {
	gretl_matrix_set(YY, t, 1, (mult > 1)? mult*y->val[t] : y->val[t]);
	if (t >= sT) {
	    gretl_matrix_set(YY, t, 0, NADBL);
	}
    }

    if (dset != NULL) {
	/* Either Y0 or X was a dataset object: we should make
	   use of dataset info in the plot if we can. But is the
	   dataset of the higher or lower frequency?
	*/
	int ok, hf = (dset->pd > 1 && dset->pd != 52);

	if (hf && dset->pd == 4 && s == 3) {
	    hf = 0;
	}
	if (!hf) {
	    ok = set_hf_data_info(&hfd, dset, s);
	    dset = ok ? &hfd : NULL;
	}
    }

    title = g_strdup_printf("%s (%s)", _("Temporal disaggregation"),
			    method_names[method]);
    err = write_tdisagg_plot(YY, mult, title, dset);

    gretl_matrix_free(YY);
    g_free(title);

    return err;
}
