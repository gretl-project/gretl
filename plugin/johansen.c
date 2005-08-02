/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "pvalues.h"
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"
#include "var.h"

/* 
   Critical values for Johansen's likelihood ratio tests
   are computed using J. Doornik's gamma approximation
*/

/* Matrices for the trace test */

const double s_mTrace_m_coef[5][6] = {
/*  n^2     n        1    n==1     n==2  n^1/2 */
  {2,  -1.00,    0.07,   0.07,       0,     0},
  {2,   2.01,       0,   0.06,    0.05,     0},
  {2,   1.05,   -1.55,  -0.50,   -0.23,     0},
  {2,   4.05,    0.50,  -0.23,   -0.07,     0},
  {2,   2.85,   -5.10,  -0.10,   -0.06,  1.35}
};

const double s_mTrace_v_coef[5][6] = {
  {3,  -0.33,  -0.55,      0,        0,     0},
  {3,    3.6,   0.75,   -0.4,     -0.3,     0},
  {3,    1.8,      0,   -2.8,     -1.1,     0},
  {3,    5.7,    3.2,   -1.3,     -0.5,     0},
  {3,    4.0,    0.8,   -5.8,    -2.66,     0}
};

const double s_mTrace_m_time[5][7] = {
/* sqrt(n)/T   n/T  n^2/T^2   n==1/T     n==1     n==2     n==3 */
  {-0.101,   0.499,   0.896,  -0.562, 0.00229, 0.00662,       0}, 
  {     0,   0.465,   0.984,  -0.273,-0.00244,       0,       0}, 
  { 0.134,   0.422,    1.02,    2.17,-0.00182,       0,-0.00321}, 
  {0.0252,   0.448,    1.09,  -0.353,       0,       0,       0}, 
  {-0.819,   0.615,   0.896,    2.43, 0.00149,       0,       0}
};

const double s_mTrace_v_time[5][7] = {
  {-0.204,   0.980,    3.11,   -2.14,  0.0499, -0.0103,-0.00902}, 
  { 0.224,   0.863,    3.38,  -0.807,       0,       0, -0.0091}, 
  { 0.422,   0.734,    3.76,    4.32,-0.00606,       0,-0.00718}, 
  {     0,   0.836,    3.99,   -1.33,-0.00298,-0.00139,-0.00268}, 
  { -1.29,    1.01,    3.92,    4.67, 0.00484,-0.00127, -0.0199}
};

/* Matrices for the lambdamax test */

const double s_mMaxev_m_coef[5][5] = {
/*   n            1        n==1         n==2        n^1/2 */
  {6.0019,     -2.7558,     0.67185,     0.11490,     -2.7764},  
  {5.9498,     0.43402,    0.048360,    0.018198,     -2.3669},  
  {5.8271,     -1.6487,     -1.6118,    -0.25949,     -1.5666},  
  {5.8658,      2.5595,    -0.34443,   -0.077991,     -1.7552},  
  {5.6364,    -0.90531,     -3.5166,    -0.47966,    -0.21447}
}; 

const double s_mMaxev_v_coef[5][5] = {
  {1.8806,     -15.499,      1.1136,    0.070508,      14.714},  
  {2.2231,     -7.9064,     0.58592,   -0.034324,      12.058},  
  {2.0785,     -9.7846,     -3.3680,    -0.24528,      13.074},  
  {1.9955,     -5.5428,      1.2425,     0.41949,      12.841},  
  {2.0899,     -5.3303,     -7.1523,    -0.25260,      12.393}
}; 


static int
gamma_par_asymp (double tracetest, double lmaxtest, JohansenCode det, 
		 int N, double *pval)
{
    /*
      Asymptotic critical values for Johansen's LR tests via gamma approximation

      params:
      tracetest, lmaxtest: trace and lambdamax est. statistics
      det: index of setup of deterministic regressors 
        J_NO_CONST     = no constant
        J_REST_CONST   = restricted constant
        J_UNREST_CONST = unrestricted constant
        J_REST_TREND   = restricted trend
        J_UNREST_TREND = unrestricted trend
      N: cointegration rank under H0;
      pval: on output, array of pvalues, for the two tests;
    */
    
    double mt, vt, ml, vl;
    const double *tracem, *tracev, *lmaxm, *lmaxv;
    double g, x[7];
    int i;

    tracem = s_mTrace_m_coef[det];
    tracev = s_mTrace_v_coef[det];
    lmaxm = s_mMaxev_m_coef[det];
    lmaxv = s_mMaxev_v_coef[det];

    mt = vt = 0.0;
    ml = vl = 0.0;

    x[0] = N * N;
    x[1] = N;
    x[2] = 1.0;
    x[3] = (N == 1)? 1.0 : 0.0;
    x[4] = (N == 2)? 1.0 : 0.0;
    x[5] = sqrt((double) N);

    for (i=0; i<6; i++) {
	mt += x[i] * tracem[i];
	vt += x[i] * tracev[i];
	if (i) {
	    ml += x[i] * lmaxm[i-1];
	    vl += x[i] * lmaxv[i-1];
	}
    }

    g = gamma_dist(mt, vt, tracetest, 2);
    if (na(g)) {
	pval[0] = NADBL;
    } else {
	pval[0] = 1.0 - g;
	if (pval[0] < 0.0) {
	    pval[0] = 0.0;
	}
    }

    g = gamma_dist(ml, vl, lmaxtest, 2);
    if (na(g)) {
	pval[1] = NADBL;
    } else {
	pval[1] = 1.0 - g;
	if (pval[1] < 0.0) {
	    pval[1] = 0.0;
	}
    }

    return 0;
}

/* for comparison with Eviews */
double ll_adj;

struct eigval {
    double v;
    int idx;
};

static int inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da < *db) - (*da > *db);
}

static int 
johansen_normalize (JVAR *jv, gretl_matrix *evecs, const double *eigvals)
{
    gretl_matrix *a = NULL, *b = NULL;
    double x, den;
    int i, j, k = gretl_matrix_rows(jv->Svv);
    int err = 0;

    a = gretl_column_vector_alloc(k);
    b = gretl_column_vector_alloc(k);

    if (a == NULL || b == NULL) {
	gretl_matrix_free(a);
	gretl_matrix_free(b);
	return E_ALLOC;
    }

    for (j=0; j<k; j++) {
	/* select column from evecs */
	for (i=0; i<k; i++) {
	    x = gretl_matrix_get(evecs, i, j);
	    gretl_vector_set(a, i, x);
	}

	/* find value of appropriate denominator */
	gretl_matrix_multiply(jv->Svv, a, b);
	den = gretl_vector_dot_product(a, b, &err);

	if (!err) {
	    den = sqrt(den);
	    for (i=0; i<k; i++) {
		x = gretl_matrix_get(evecs, i, j);
		gretl_matrix_set(evecs, i, j, x / den);
	    }
	} 
    }

    gretl_matrix_free(a);
    gretl_matrix_free(b);

    return err;
}

enum {
    PRINT_ALPHA,
    PRINT_BETA
};

static void print_coint_eqns (JVAR *jv, const DATAINFO *pdinfo, PRN *prn)
{
    char s[16];
    int rows = gretl_matrix_rows(jv->Beta);
    int i, j;
    double r;

    pprintf(prn, "%s\n\n", _("Cointegrating vectors"));

    for (j=0; j<jv->rank; j++) {
	sprintf(s, "CV%d", j + 1);
	if (j == 0) {
	    pprintf(prn, "%22s", s);
	} else {
	    pprintf(prn, "%13s", s);
	}
    }
    pputc(prn, '\n');

    /* ? check the column ordering for rank > 1 */

    for (i=0; i<rows; i++) {
	if (i < jv->list[0]) {
	    sprintf(s, "%s(-1)", pdinfo->varname[jv->list[i+1]]);
	    pprintf(prn, "%-10s", s);
	} else if (jv->code == J_REST_CONST) {
	    pprintf(prn, "%-10s", "const");
	} else if (jv->code == J_REST_TREND) {
	    pprintf(prn, "%-10s", "trend");
	}
	for (j=0; j<jv->rank; j++) {
	    r = gretl_matrix_get(jv->Beta, j, j);
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(jv->Beta, i, j) / r);
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* print cointegrating vectors or adjustments, either "raw" or
   normalized */

static void print_beta_or_alpha (JVAR *jv, 
				 struct eigval *evals, int k, 
				 const DATAINFO *pdinfo, PRN *prn,
				 int which, int normalize)
{
    gretl_matrix *c = (which == PRINT_BETA)? jv->Beta : jv->Alpha;
    int rows = gretl_matrix_rows(c);
    int i, j, col;
    double x;

    if (normalize) {
	pprintf(prn, "\n%s\n", (which == PRINT_BETA)? 
		_("renormalized beta") :
		_("renormalized alpha"));
    } else {
	pprintf(prn, "\n%s\n", (which == PRINT_BETA)? 
		_("beta (cointegrating vectors)") : 
		_("alpha (adjustment vectors)"));
    }

    for (j=0; j<rows; j++) {
	if (j < jv->list[0]) {
	    pprintf(prn, "%-10s", pdinfo->varname[jv->list[j+1]]);
	} else if (jv->code == J_REST_CONST) {
	    pprintf(prn, "%-10s", "const");
	} else if (jv->code == J_REST_TREND) {
	    pprintf(prn, "%-10s", "trend");
	}
	for (i=0; i<k; i++) {
	    col = evals[i].idx;
	    if (normalize) {
		x = gretl_matrix_get(jv->Beta, i, col);
		if (which == PRINT_BETA) {
		    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, j, col) / x);
		} else {
		    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, j, col) * x);
		}
	    } else {
		pprintf(prn, "%#12.5g ", gretl_matrix_get(c, j, col));
	    }
	}
	pputc(prn, '\n');
    }
}

#if 0
static void normalize_alpha (JVAR *jv)
{
    double x, r;
    int i, j;

    for (i=0; i<jv->neqns; i++) {
	for (j=0; j<jv->rank; j++) {
	    r = gretl_matrix_get(jv->Beta, j, j);
	    x = gretl_matrix_get(jv->Alpha, i, j);
	    gretl_matrix_set(jv->Alpha, i, j, x * r);
	}
    }
}

static void normalize_beta (JVAR *jv)
{
    double x, r;
    int i, j;

    for (i=0; i<jv->neqns; i++) {
	for (j=0; j<jv->rank; j++) {
	    r = gretl_matrix_get(jv->Beta, j, j);
	    x = gretl_matrix_get(jv->Beta, i, j);
	    gretl_matrix_set(jv->Beta, i, j, x / r);
	}
    }
}
#endif

/* calculate alpha (adjustments) matrix as per Johansen, 1991, eqn
   2.8, p. 1554 */

static int compute_alpha (JVAR *jv, int h)
{
    gretl_matrix *alpha = NULL;
    gretl_matrix *tmp1 = NULL;
    gretl_matrix *tmp2 = NULL;
    int n = gretl_matrix_rows(jv->Beta);
    int err = 0;

    tmp1 = gretl_matrix_alloc(n, h);
    tmp2 = gretl_matrix_alloc(h, h);
    alpha = gretl_matrix_alloc(n, h);

    if (tmp1 == NULL || tmp2 == NULL || alpha == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	gretl_matrix_multiply(jv->Svv, jv->Beta, tmp1);
	gretl_matrix_multiply_mod(jv->Beta, GRETL_MOD_TRANSPOSE,
				  tmp1, GRETL_MOD_NONE,
				  tmp2);
	err = gretl_invert_general_matrix(tmp2);
    }

    if (!err) {
	gretl_matrix_multiply(jv->Beta, tmp2, tmp1);
	gretl_matrix_multiply(jv->Suv, tmp1, alpha);
    }

    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);

    if (!err) {
	jv->Alpha = alpha;
    } else {
	gretl_matrix_free(alpha);
    }

    return err;
}

static void print_lr_matrix (JVAR *jv, gretl_matrix *Zeta_0, int n, int h,
			     const DATAINFO *pdinfo, PRN *prn)
{
    int i, j;

    pprintf(prn, "%s\n", _("long-run matrix (alpha * beta')"));
    pprintf(prn, "%22s", pdinfo->varname[jv->list[1]]);

    for (j=1; j<n; j++) {
	if (j < jv->list[0]) {
	    pprintf(prn, "%13s", pdinfo->varname[jv->list[j+1]]);
	} else if (jv->code == J_REST_CONST) {
	    pprintf(prn, "%13s", "const");
	} else if (jv->code == J_REST_TREND) {
	    pprintf(prn, "%13s", "trend");
	}
    }

    pputc(prn, '\n');

    for (i=0; i<h; i++) {
	pprintf(prn, "%-10s", pdinfo->varname[jv->list[i+1]]);
	for (j=0; j<n; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(Zeta_0, i, j));
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* make a matrix containing the eigenvectors corresponding to
   the hmax largest eigenvalues (Hamilton's "A") */

static gretl_matrix *select_eigvecs (JVAR *jv, const gretl_matrix *evecs,
				     struct eigval *evals, int hmax)
{
    gretl_matrix *A;
    double x;
    int i, j, col;

    A = gretl_matrix_alloc(jv->neqns, hmax);

    if (A != NULL) {
	for (j=0; j<hmax; j++) {
	    col = evals[j].idx;
	    for (i=0; i<jv->neqns; i++) {
		x = gretl_matrix_get(evecs, i, col);
		gretl_matrix_set(A, i, j, x);
	    }
	}
    }

    return A;
}

static int compute_estimator_variance (JVAR *jv)
{
    gretl_matrix *V = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *D = NULL;
    gretl_matrix *tmp = NULL;
    double x;
    int k = gretl_matrix_cols(jv->Data);
    int T = jv->t2 - jv->t1 + 1;
    int p = jv->order - 1;
    int h = jv->rank;
    int n = jv->neqns;
    int pn = p * n;
    int dcols = k - (n - h);
    int i, j;
    int err = 0;

    X = gretl_matrix_alloc(T, n);
    tmp = gretl_matrix_alloc(T, h);
    D = gretl_matrix_alloc(T, dcols);

    if (X == NULL || tmp == NULL || D == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* copy the X_{t-1} data into X */
    for (j=0; j<n; j++) {
	for (i=0; i<T; i++) {
	    x = gretl_matrix_get(jv->Data, i, j + pn);
	    gretl_matrix_set(X, i, j, x);
	}
    }

    /* multiply X_{t-1} into \beta */
    gretl_matrix_multiply(X, jv->Beta, tmp);

    /* copy the \Delta X_{t-i} data into D */
    for (j=0; j<pn; j++) {
	for (i=0; i<T; i++) {
	    x = gretl_matrix_get(jv->Data, i, j);	
	    gretl_matrix_set(D, i, j, x);
	}
    }

    /* copy the beta-modified X_{t-1} into D */
    for (j=pn; j<pn+h; j++) {
	for (i=0; i<T; i++) {
	    x = gretl_matrix_get(tmp, i, j-pn);
	    gretl_matrix_set(D, i, j, x);
	}
    }

    /* add deterministic vars to D */
    for (j=pn+h; j<dcols; j++) {
	for (i=0; i<T; i++) {
	    x = gretl_matrix_get(jv->Data, i, j + (n - h));	
	    gretl_matrix_set(D, i, j, x);
	}
    }
	
    V = gretl_matrix_vcv(D);

    if (V == NULL) {
	err = 1;
    } else {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	int nc = pn + h;
	double a, vii;

	nc += (jv->code >= J_UNREST_CONST) + (jv->code == J_UNREST_TREND);
	jv->Ase = gretl_matrix_alloc(nc, n);
	if (jv->Ase == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	gretl_matrix_zero(jv->Ase);

	for (j=0; j<n; j++) {
	    a = gretl_matrix_get(jv->Omega, j, j) / T;
	    for (i=0; i<nc; i++) {
		vii = gretl_matrix_get(V, i, i);
		gretl_matrix_set(jv->Ase, i, j, sqrt(vii * a));
	    }
	}

#if 0
	for (j=0; j<n; j++) {
	    a = gretl_matrix_get(jv->Omega, j, j) / T;
	    for (i=0; i<n; i++) {
		double vii = 2.0 * gretl_matrix_get(jv->Omega, i, i);

		fprintf(stderr, "se[O](%d,%d) = %g\n", j+1, i+1, sqrt(vii * a));
	    }
	    fputc('\n', stderr);
	}
#endif	
    }

 bailout:

    gretl_matrix_free(D);
    gretl_matrix_free(X);
    gretl_matrix_free(tmp);
    gretl_matrix_free(V);

    return err;
}

/* compute Hamilton's Omega (Johansen 1991 calls it Lambda): the
   cross-equation variance matrix 
*/

static int compute_omega (JVAR *jv, const gretl_matrix *Z0)
{
    gretl_matrix *uhat = NULL;
    gretl_matrix *omega = NULL;
    gretl_matrix *tmp = NULL;
    int n = jv->neqns;
    int T = jv->t2 - jv->t1 + 1;
    int df_adj = 1; /* not really wanted? */
    int err = 0;

    uhat = gretl_matrix_alloc(n, T);
    omega = gretl_matrix_alloc(n, n);
    tmp = gretl_matrix_alloc(n, T);

    if (uhat == NULL || omega == NULL || tmp == NULL) {
	err = E_ALLOC;
    } else {
	int den = T;

	gretl_matrix_copy_values(uhat, jv->u);
	gretl_matrix_multiply(Z0, jv->v, tmp);
	gretl_matrix_subtract_from(uhat, tmp);

	gretl_matrix_multiply_mod(uhat, GRETL_MOD_NONE,
				  uhat, GRETL_MOD_TRANSPOSE,
				  omega);

	if (df_adj) {
	    /* Eviews 4.0 applies df adjustment here */
	    den -= 1 + jv->neqns * (jv->order - 1)
		+ (jv->code >= J_UNREST_CONST) 
		+ (jv->code == J_UNREST_TREND);
	}

	gretl_matrix_divide_by_scalar(omega, den);
    }

    gretl_matrix_free(tmp);

    if (!err) {
	jv->uhat = uhat;
	jv->Omega = omega;
    } else {
	gretl_matrix_free(omega);
	gretl_matrix_free(uhat);
    }

    return err;
}

#if 0 /* this seems to be broken (at least does not agree w Eviews) */

/* find \alpha(\alpha' \alpha)^{-1} \alpha' \mu ?? */

static int coint_consts_from_mu (const gretl_matrix *mu, const gretl_matrix *alpha)
{
    gretl_matrix *tmp1 = NULL;
    gretl_matrix *tmp2 = NULL;
    gretl_matrix *tmp3 = NULL;
    gretl_matrix *c = NULL;
    int h = gretl_matrix_cols(alpha);
    int n = gretl_matrix_rows(mu);
    int err = 0;

    tmp1 = gretl_matrix_alloc(h, h);
    tmp2 = gretl_matrix_alloc(h, 1);
    tmp3 = gretl_matrix_alloc(h, 1);
    c = gretl_matrix_alloc(n, 1);

    /* FIXME checking */

    gretl_matrix_multiply_mod(alpha, GRETL_MOD_TRANSPOSE,
			      alpha, GRETL_MOD_NONE,
			      tmp1);

    gretl_invert_symmetric_matrix(tmp1);

    gretl_matrix_multiply_mod(alpha, GRETL_MOD_TRANSPOSE,
			      mu, GRETL_MOD_NONE,
			      tmp2);

    gretl_matrix_multiply(tmp1, tmp2, tmp3);
    gretl_matrix_reuse(tmp1, 4, 1);
    gretl_matrix_multiply(alpha, tmp3, c);

    gretl_matrix_print(c, "consts in coint relations (?)", NULL);

    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);
    gretl_matrix_free(tmp3);
    gretl_matrix_free(c);
    
    return err;
}

#endif

static int print_omega (JVAR *jv, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *M;
    char s[32];
    int i, j;

    pprintf(prn, "%s\n\n", _("Cross-equation covariance matrix"));

    for (i=0; i<jv->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[jv->list[i+1]]);
	if (i == 0) {
	    pprintf(prn, "%22s", s);
	} else {
	    pprintf(prn, "%13s", s);
	}
    }
    pputc(prn, '\n');

    for (i=0; i<jv->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[jv->list[i+1]]);
	pprintf(prn, "%-10s", s);
	for (j=0; j<jv->neqns; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(jv->Omega, i, j));
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    M = gretl_matrix_copy(jv->Omega);
    if (M != NULL) {
	pprintf(prn, "determinant = %g\n", gretl_matrix_determinant(M));
	gretl_matrix_free(M);
    }
}

static double get_tss (JVAR *jv, double *xd, const double **Z, int i)
{
    double x, xbar, tss = 0.0;
    int T = jv->t2 - jv->t1 + 1;
    int t, v = jv->list[i+1];
    
    for (t=jv->t1; t<=jv->t2; t++) {
	xd[t - jv->t1] = Z[v][t] - Z[v][t - 1];
    }

    xbar = gretl_mean(0, T - 1, xd);

    if (na(xbar)) {
	tss = NADBL;
    } else {
	for (t=0; t<T; t++) {
	    x = xd[t];
	    x -= xbar;
	    tss += x * x;
	}
    }

    return tss;
}

static int 
print_vecm (JVAR *jv, const double **Z, const DATAINFO *pdinfo, PRN *prn)
{
    char s[32];
    double *ess = NULL;
    double *xd = NULL;
    double x;
    int doconst = 0, dotrend = 0;
    int T = jv->t2 - jv->t1 + 1;
    int rows = gretl_matrix_rows(jv->A);
    int p = jv->order - 1;
    int i, j, m, k;

    ess = malloc(jv->neqns * sizeof *ess);
    if (ess == NULL) {
	return E_ALLOC;
    }

    xd = malloc(T * sizeof *xd);
    if (xd == NULL) {
	free(ess);
	return E_ALLOC;
    }    

    doconst = jv->code >= J_REST_CONST;
    dotrend = jv->code == J_REST_TREND;

    pprintf(prn, "%s\n", _("Error correction equations"));
    pprintf(prn, "%s\n\n", _("standard errors in (); t-ratios in []"));

    for (i=0; i<jv->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[jv->list[i+1]]);
	if (i == 0) {
	    pprintf(prn, "%22s", s);
	} else {
	    pprintf(prn, "%13s", s);
	}
    }
    pputs(prn, "\n\n");

    m = k = 1;

    for (i=0; i<rows; i++) {
	if (i < jv->rank) {
	    sprintf(s, "CV%d", i + 1);
	    pprintf(prn, "%-10s", s);
	} else if (k <= jv->list[0]) {
	    sprintf(s, "d_%s(-%d)", pdinfo->varname[jv->list[k]], m++);
	    pprintf(prn, "%-10s", s);
	} else if (doconst) {
	    pprintf(prn, "%-10s", "const");
	    doconst = 0;
	} else if (dotrend) {
	    pprintf(prn, "%-10s", "trend");
	    dotrend = 0;
	}
	for (j=0; j<jv->neqns; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(jv->A, i, j));
	}
	if (i > 0) {
	    pputc(prn, '\n');
	    pprintf(prn, "%-10s", " ");
	    for (j=0; j<jv->neqns; j++) {
		sprintf(s, "(%#.5g)", gretl_matrix_get(jv->Ase, i-1, j));
		pprintf(prn, " %12s", s);
	    }
	    pputc(prn, '\n');
	    pprintf(prn, "%-10s", " ");
	    for (j=0; j<jv->neqns; j++) {
		double tr =  gretl_matrix_get(jv->A, i, j)/
		    gretl_matrix_get(jv->Ase, i-1, j);

		sprintf(s, "[%#.4g]", tr);
		pprintf(prn, " %12s", s);
	    }	    
	}	
	pputs(prn, "\n\n");

	if (i > 0 && i % p == 0) {
	    k++;
	    m = 1;
	}
    }

    /* per-equation SSR */
    pprintf(prn, "%-10s", "SSR");
    for (i=0; i<jv->neqns; i++) {
	double u;

	ess[i] = 0.0;
	for (j=0; j<T; j++) {
	    u = gretl_matrix_get(jv->uhat, i, j);
	    ess[i] += u * u;
	}
	pprintf(prn, "%#12.5g ", ess[i]);
    }    

    /* per-equation standard errors */
    pputc(prn, '\n');
    pprintf(prn, "%-10s", "SE");
    for (i=0; i<jv->neqns; i++) {
	x = gretl_matrix_get(jv->Omega, i, i);
	pprintf(prn, "%#12.5g ", sqrt(x));
    }

    /* per-equation R-squared */
    pputc(prn, '\n');
    pprintf(prn, "%-10s", "R-squared");
    for (i=0; i<jv->neqns; i++) {
	double tss = get_tss(jv, xd, Z, i);

	pprintf(prn, "%#12.5g ", 1.0 - ess[i] / tss);
    }
    pputs(prn, "\n\n");    

    print_omega(jv, pdinfo, prn);

    if (!na(jv->ll)) {
	int T = jv->t2 - jv->t1 + 1;
	int n = jv->neqns;
	int k = gretl_matrix_rows(jv->A) + jv->rank;
	double aic, bic;

	aic = (-2.0 * jv->ll + 2.0 * k * n) / T;
	bic = (-2.0 * jv->ll + log(T) * k * n) / T;

	pprintf(prn, "\nlog-likelihood = %g\n", jv->ll);
	pprintf(prn, "AIC = %g\n", aic);
	pprintf(prn, "BIC = %g\n", bic);

	aic = (-2.0 * ll_adj + 2.0 * k * n) / T;
	bic = (-2.0 * ll_adj + log(T) * k * n) / T;

	pputs(prn, "\nd.f. adjusted\n");
	pprintf(prn, "log-likelihood = %g\n", ll_adj);
	pprintf(prn, "AIC = %g\n", aic);
	pprintf(prn, "BIC = %g\n", bic);
    }

    free(ess);
    free(xd);

    return 0;
}

static int 
compute_vecm (JVAR *jv, const double **Z, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *Z0 = NULL;
    gretl_matrix *Zeta = NULL;
    gretl_matrix *tmp = NULL;
    int n = jv->neqns;
    int p = jv->order - 1;
    int nc = n * p + jv->rank;
    int i, err = 0;

    nc += (jv->code >= J_UNREST_CONST) + (jv->code == J_UNREST_TREND);

    /* full-size coefficient matrix */
    jv->A = gretl_matrix_alloc(nc, n);
    if (jv->A == NULL) {
	return E_ALLOC;
    }

    Z0 = gretl_matrix_alloc(n, n);
    Zeta = gretl_matrix_alloc(n, n);
    tmp = gretl_matrix_alloc(n, n);

    if (jv->code >= J_UNREST_CONST) {
	jv->mu = gretl_column_vector_alloc(n);
    }

    if (Z0 == NULL || Zeta == NULL || tmp == NULL) {
	err = E_ALLOC;
    } else {
	int j, k, arow = 0;
	double r, x;

	gretl_matrix_multiply_mod(jv->Beta, GRETL_MOD_NONE,
				  jv->Beta, GRETL_MOD_TRANSPOSE,
				  tmp);

	/* \Zeta_0 = \Sigma_{uv} A A' (Hamilton eq. 20.2.12) */
	gretl_matrix_multiply(jv->Suv, tmp, Z0);

	gretl_matrix_zero(jv->A);

	/* put adjustment (alpha) terms in first */
	for (i=0; i<jv->rank; i++) {
	    for (j=0; j<n; j++) {
		r = gretl_matrix_get(jv->Beta, i, i);
		x = gretl_matrix_get(jv->Alpha, j, i);
		gretl_matrix_set(jv->A, i, j, x * r);
	    }
	}

	for (i=0; i<p; i++) {
	    gretl_matrix_copy_values(Zeta, jv->Pi[i]);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(Zeta, tmp);
#if 0
	    gretl_matrix_print(Zeta, "Zeta", prn);
#endif
	    /* place coeffs in big matrix */
	    arow = jv->rank + i;
	    for (j=0; j<n; j++) {
		for (k=0; k<n; k++) {
		    x = gretl_matrix_get(Zeta, k, j);
		    gretl_matrix_set(jv->A, arow, k, x);
		}
		arow += p;
	    }
	}

	/* calculate constants (Johansen calls then \mu), as
	   per Hamilton, eq. 20.2.14, p. 637 (Hamilton calls
	   them \alpha) */

	if (jv->mu != NULL) {
	    gretl_matrix_copy_values(jv->mu, jv->Pi[i]);
	    gretl_matrix_reuse(tmp, n, 1);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(jv->mu, tmp);
	    gretl_matrix_reuse(tmp, n, n);
#if 0 /* broken */
	    coint_consts_from_mu(jv->mu, jv->Alpha);
#endif
	    arow = (jv->code == J_UNREST_TREND)? nc - 2 : nc - 1;
	    for (j=0; j<n; j++) {
		x = gretl_vector_get(jv->mu, j);
		gretl_matrix_set(jv->A, arow, j, x);
	    }
	}

	err = compute_omega(jv, Z0);

	if (!err) {
	    err = compute_estimator_variance(jv);
	}
    }

    if (!err) {
	print_vecm(jv, Z, pdinfo, prn);
    }

    gretl_matrix_free(Z0);
    gretl_matrix_free(Zeta);
    gretl_matrix_free(tmp);

    return err;
}

/* calculate and print the "long-run matrix", \alpha \beta', (what
   Hamilton calls \Zeta_0) */

static int 
compute_long_run_matrix (JVAR *jv, struct eigval *evals, int h, 
			 const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *Z0 = NULL;
    gretl_matrix *tmp = NULL;
    int n = jv->neqns;
    int err = 0;

    Z0 = gretl_matrix_alloc(h, n); 
    tmp = gretl_matrix_alloc(n, n);

    if (Z0 == NULL || tmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	gretl_matrix_multiply_mod(jv->Beta, GRETL_MOD_NONE,
				  jv->Beta, GRETL_MOD_TRANSPOSE,
				  tmp);

	/* \Zeta_0 = \Sigma_{uv} A A' */
	gretl_matrix_multiply(jv->Suv, tmp, Z0);
	print_lr_matrix(jv, Z0, n, h, pdinfo, prn);
    }

    gretl_matrix_free(Z0);
    gretl_matrix_free(tmp);
    
    return err;
}

static int
print_beta_and_alpha (JVAR *jv, struct eigval *evals, int h, 
		      const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    pprintf(prn, "\n%s", _("eigenvalue"));
    for (i=0; i<h; i++) {
	pprintf(prn, "%#12.5g ", evals[i].v);
    }
    pputc(prn, '\n');

    /* "raw" vectors */
    print_beta_or_alpha(jv, evals, h, pdinfo, prn, PRINT_BETA, 0);
    print_beta_or_alpha(jv, evals, h, pdinfo, prn, PRINT_ALPHA, 0);

    /* normalized versions */
    print_beta_or_alpha(jv, evals, h, pdinfo, prn, PRINT_BETA, 1);
    print_beta_or_alpha(jv, evals, h, pdinfo, prn, PRINT_ALPHA, 1);

    pputc(prn, '\n');
    
    return err;
}

static void print_test_case (JohansenCode jcode, PRN *prn)
{
    pputc(prn, '\n');

    switch (jcode) {
    case J_NO_CONST:
	pputs(prn, "Case 1: No constant");
	break;
    case J_REST_CONST:
	pputs(prn, "Case 2: Restricted constant");
	break;
    case J_UNREST_CONST:
	pputs(prn, "Case 3: Unrestricted constant");
	break;
    case J_REST_TREND:
	pputs(prn, "Case 4: Restricted trend, unrestricted constant");
	break;
    case J_UNREST_TREND:
	pputs(prn, "Case 5: Unrestricted trend and constant");
	break;
    default:
	break;
    }

    pputc(prn, '\n');
}

static char *safe_print_pval (double p, int i)
{
    static char pv[2][32];

    if (!na(p)) {
	sprintf(pv[i], "%6.4f", p);
    } else {
	strcpy(pv[i], "  NA  ");
    }

    return pv[i];
}

static int johansen_ll_init (JVAR *jv, int T)
{
    gretl_matrix *Suu;
    double T_2 = (double) T / 2.0;
    double ldet, n = jv->neqns;
    int err = 0;

    Suu = gretl_matrix_copy(jv->Suu);

    if (Suu == NULL) {
	err = E_ALLOC;
    } else {
	int T_k = T - (n * (jv->order - 1) + 1 + jv->rank);

	ldet = gretl_matrix_log_determinant(Suu);
	jv->ll = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;

	/* "df-adjusted ll" for comparison with Eviews */
	gretl_matrix_copy_values(Suu, jv->Suu);
	gretl_matrix_multiply_by_scalar(Suu, (double) T / T_k);
	ldet = gretl_matrix_log_determinant(Suu);
	ll_adj = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;

	gretl_matrix_free(Suu);
    }

    return err;
}

int johansen_eigenvals (JVAR *jv, const double **Z, const DATAINFO *pdinfo, 
			PRN *prn)
{
    int k = gretl_matrix_cols(jv->Suu);
    int kv = gretl_matrix_cols(jv->Svv);
    gretl_matrix *alpha = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *TmpL = NULL, *TmpR = NULL;
    gretl_matrix *Suu = NULL, *Svv = NULL;
    double *eigvals = NULL;
    int err = 0;

    TmpL = gretl_matrix_alloc(kv, k);
    TmpR = gretl_matrix_alloc(kv, kv);
    M = gretl_matrix_alloc(kv, kv);

    /* let's preserve Svv and Suu, so copy them before
       inverting */
    Svv = gretl_matrix_copy(jv->Svv);
    Suu = gretl_matrix_copy(jv->Suu);

    if (Suu == NULL || Svv == NULL || TmpL == NULL || 
	TmpR == NULL || M == NULL) {
	err = 1;
	goto eigenvals_bailout;
    }

    if (kv > k) {
	gretl_matrix_reuse(TmpR, k, kv);
    }

    /* calculate Suu^{-1} Suv */
    err = gretl_invert_general_matrix(Suu);
    if (!err) {
	err = gretl_matrix_multiply(Suu, jv->Suv, TmpR);
    }

    /* calculate Svv^{-1} Suv' */
    if (!err) {
	err = gretl_invert_general_matrix(Svv);
    }
    if (!err) {
	err = gretl_matrix_multiply_mod(Svv, GRETL_MOD_NONE,
					jv->Suv, GRETL_MOD_TRANSPOSE, 
					TmpL);
    }

    if (!err) {
	err = gretl_matrix_multiply(TmpL, TmpR, M);
    }

    if (err) {
	goto eigenvals_bailout;
    }

    if (kv > k) {
	gretl_matrix_reuse(TmpR, kv, kv);
    }

    eigvals = gretl_general_matrix_eigenvals(M, TmpR);

    if (eigvals != NULL) {
	int T = jv->t2 - jv->t1 + 1;
	double cumeig = 0.0;
	double *lambdamax = NULL, *trace = NULL;
	struct eigval *evals;
	double pval[2];
	int hmax, i;

	/* in case kv > k */
	hmax = k;
	k = kv;

	/* needs more work */
	if (jv->rank > 0) {
	    hmax = jv->rank;
	}

	trace = malloc(k * sizeof *trace);
	lambdamax = malloc(k * sizeof *lambdamax);
	evals = malloc(k * sizeof *evals);
	if (trace == NULL || lambdamax == NULL || evals == NULL) {
	    free(trace);
	    free(lambdamax);
	    free(evals);
	    err = 1;
	    goto eigenvals_bailout;
	}

	for (i=0; i<k; i++) {
#if JDEBUG
	    fprintf(stderr, "unsorted eigvals[%d] = %g\n", i, eigvals[i]);
#endif
	    evals[i].v = eigvals[i];
	    evals[i].idx = i;
	}

	qsort(evals, k, sizeof *evals, inverse_compare_doubles);

	for (i=hmax-1; i>=0; i--){
      	    lambdamax[i] = -T * log(1.0 - evals[i].v); 
	    cumeig += lambdamax[i];
 	    trace[i] = cumeig; 
	}

	print_test_case(jv->code, prn);
	johansen_ll_init(jv, T);

	/* first column shows cointegration rank under H0, 
	   second shows associated eigenvalue */
	pprintf(prn, "\n%s %s %s %s   %s  %s\n", _("Rank"), _("Eigenvalue"), 
		_("Trace test"), _("p-value"),
		_("Lmax test"), _("p-value"));

	for (i=0; i<hmax; i++) {
	    if (!na(jv->ll)) {
		jv->ll += lambdamax[i] / 2.0;
		ll_adj += lambdamax[i] / 2.0;
	    }
	    gamma_par_asymp(trace[i], lambdamax[i], jv->code, hmax - i, pval);
	    pprintf(prn, "%4d%#11.5g%#11.5g [%s]%#11.5g [%s]\n", \
		    i, evals[i].v, trace[i], safe_print_pval(pval[0], 0), 
		    lambdamax[i], safe_print_pval(pval[1], 1));
	}
	pputc(prn, '\n');

	/* TmpR holds the full set of eigenvectors: normalize before
	   shrinking to rank 
	*/
	johansen_normalize(jv, TmpR, eigvals);
	
	if (jv->rank == 0) {
	    jv->Beta = TmpR;
	    TmpR = NULL;
	    err = compute_alpha(jv, hmax);
	    print_beta_and_alpha(jv, evals, hmax, pdinfo, prn);
	    compute_long_run_matrix(jv, evals, hmax, pdinfo, prn);
	} else {
	    jv->Beta = select_eigvecs(jv, TmpR, evals, hmax);
	    err = compute_alpha(jv, hmax);
	    print_coint_eqns(jv, pdinfo, prn);
	    compute_vecm(jv, Z, pdinfo, prn);
	}

	free(eigvals);
	free(evals);
	free(lambdamax);
	free(trace);

    } else {
	pputs(prn, _("Failed to find eigenvalues\n"));
    }

 eigenvals_bailout:    

    gretl_matrix_free(TmpL);
    gretl_matrix_free(TmpR);
    gretl_matrix_free(M);
    gretl_matrix_free(alpha);
    gretl_matrix_free(Svv);
    gretl_matrix_free(Suu);

    return err;
}
