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

/* print cointegrating vectors or adjustments, either "raw" or
   normalized */

static void print_beta_or_alpha (JVAR *jv, 
				 const gretl_matrix *beta,
				 const gretl_matrix *alpha,
				 struct eigval *evals, int k, 
				 const DATAINFO *pdinfo, PRN *prn,
				 int normalized)
{
    const gretl_matrix *c = (alpha == NULL)? beta : alpha;
    int rows = gretl_matrix_rows(c);
    int i, j, col;
    double x;

    if (normalized) {
	pprintf(prn, "\n%s\n", (alpha == NULL)? 
		_("renormalized beta") :
		_("renormalized alpha"));
    } else {
	pprintf(prn, "\n%s\n", (alpha == NULL)? 
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
	    if (normalized) {
		x = gretl_matrix_get(beta, i, col);
		if (alpha == NULL) {
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

/* calculate alpha matrix as per Johansen, 1991, eqn 2.8 p. 1554 */

static gretl_matrix *
compute_alpha (JVAR *jv, const gretl_matrix *evecs, int h)
{
    gretl_matrix *alpha = NULL;
    gretl_matrix *tmp1 = NULL;
    gretl_matrix *tmp2 = NULL;
    int n = gretl_matrix_rows(evecs);
    int err = 0;

    tmp1 = gretl_matrix_alloc(n, h);
    tmp2 = gretl_matrix_alloc(h, h);
    alpha = gretl_matrix_alloc(n, h);

    if (tmp1 == NULL || tmp2 == NULL || alpha == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	gretl_matrix_multiply(jv->Svv, evecs, tmp1);
	gretl_matrix_multiply_mod(evecs, GRETL_MOD_TRANSPOSE,
				  tmp1, GRETL_MOD_NONE,
				  tmp2);
	err = gretl_invert_general_matrix(tmp2);
    }

    if (!err) {
	gretl_matrix_multiply(evecs, tmp2, tmp1);
	gretl_matrix_multiply(jv->Suv, tmp1, alpha);
    }

    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);

    if (err && alpha != NULL) {
	gretl_matrix_free(alpha);
	alpha = NULL;
    }

    return alpha;
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

static int compute_omega (JVAR *jv, const gretl_matrix *Z0, PRN *prn)
{
    gretl_matrix *resids = NULL;
    gretl_matrix *omega = NULL;
    gretl_matrix *tmp = NULL;
    int n = jv->neqns;
    int T = jv->t2 - jv->t1 + 1;
    int df_adj = 1;
    int i, err = 0;
    
    resids = gretl_matrix_alloc(n, T);
    omega = gretl_matrix_alloc(n, n);
    tmp = gretl_matrix_alloc(n, T);

    if (resids == NULL || omega == NULL || tmp == NULL) {
	err = E_ALLOC;
    } else {
	int den = T;

	gretl_matrix_copy_values(resids, jv->u);
	gretl_matrix_multiply(Z0, jv->v, tmp);
	gretl_matrix_subtract_from(resids, tmp);

	gretl_matrix_multiply_mod(resids, GRETL_MOD_NONE,
				  resids, GRETL_MOD_TRANSPOSE,
				  omega);

	if (df_adj) {
	    /* Eviews applies df adjustment here */
	    den -= 1 + jv->neqns * (jv->order - 1)
		+ (jv->code >= J_UNREST_CONST) 
		+ (jv->code == J_UNREST_TREND);
	}

	gretl_matrix_divide_by_scalar(omega, den);
	gretl_matrix_print(omega, "omega", prn);

	for (i=0; i<n; i++) {
	    double x = gretl_matrix_get(omega, i, i);

	    pprintf(prn, "SE[%d] = %g\n", i+1, sqrt(x));
	}

	pprintf(prn, "\ndeterminant = %g\n", gretl_matrix_determinant(omega));
    }

    gretl_matrix_free(resids);
    gretl_matrix_free(omega);
    gretl_matrix_free(tmp);

    return err;
}

#if 0 /* this seems to be broken (at least does not agree w Eviews) */

/* find \alpha(\alpha' \alpha)^{-1} \alpha' \mu */

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

static int compute_vecm (JVAR *jv, const gretl_matrix *alpha, 
			 const DATAINFO *pdinfo, PRN *prn)
{
    char msg[32];
    gretl_matrix *Z0 = NULL;
    gretl_matrix *Zeta = NULL;
    gretl_matrix *mu = NULL;
    gretl_matrix *tmp = NULL;
    int n = jv->neqns;
    int p = jv->order - 1;
    int i, err = 0;

    Z0 = gretl_matrix_alloc(n, n);
    Zeta = gretl_matrix_alloc(n, n);
    tmp = gretl_matrix_alloc(n, n);

    if (jv->code >= J_UNREST_CONST) {
	mu = gretl_column_vector_alloc(n);
    }

    if (Z0 == NULL || Zeta == NULL || tmp == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_multiply_mod(jv->A, GRETL_MOD_NONE,
				  jv->A, GRETL_MOD_TRANSPOSE,
				  tmp);

	/* \zeta_0 = \Sigma_{uv} A A' */
	gretl_matrix_multiply(jv->Suv, tmp, Z0);

	for (i=0; i<p; i++) {
	    gretl_matrix_copy_values(Zeta, jv->Pi[i]);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(Zeta, tmp);
	    sprintf(msg, "Zeta[%d]", i+1);
	    gretl_matrix_print(Zeta, msg, prn);
	}

	if (mu != NULL) {
	    gretl_matrix_copy_values(mu, jv->Pi[i]);
	    gretl_matrix_reuse(tmp, n, 1);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(mu, tmp);
	    gretl_matrix_print(mu, "mu", prn);
	    gretl_matrix_reuse(tmp, n, n);
#if 0
	    coint_consts_from_mu(mu, alpha);
#endif
	}

	compute_omega(jv, Z0, prn);
    }

    gretl_matrix_free(Z0);
    gretl_matrix_free(Zeta);
    gretl_matrix_free(mu);
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
	gretl_matrix_multiply_mod(jv->A, GRETL_MOD_NONE,
				  jv->A, GRETL_MOD_TRANSPOSE,
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
print_beta_and_alpha (JVAR *jv, const gretl_matrix *beta, const gretl_matrix *alpha,
		      struct eigval *evals, int h, const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    pprintf(prn, "\n%s", _("eigenvalue"));
    for (i=0; i<h; i++) {
	pprintf(prn, "%#12.5g ", evals[i].v);
    }
    pputc(prn, '\n');

    if (jv->rank == 0) {
	print_beta_or_alpha(jv, beta, NULL, evals, h, pdinfo, prn, 0);
	if (alpha != NULL) {
	    print_beta_or_alpha(jv, beta, alpha, evals, h, pdinfo, prn, 0);
	}
    }

    print_beta_or_alpha(jv, beta, NULL, evals, h, pdinfo, prn, 1);
    if (alpha != NULL) {
	print_beta_or_alpha(jv, beta, alpha, evals, h, pdinfo, prn, 1);
    }

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

char *safe_print_pval (double p, int i)
{
    static char pv[2][32];

    if (!na(p)) {
	sprintf(pv[i], "%6.4f", p);
    } else {
	strcpy(pv[i], "  NA  ");
    }

    return pv[i];
}

int johansen_eigenvals (JVAR *jv, const DATAINFO *pdinfo, PRN *prn)
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

	/* first column shows cointegration rank under H0, 
	   second shows associated eigenvalue */
	pprintf(prn, "\n%s %s %s %s   %s  %s\n", _("Rank"), _("Eigenvalue"), 
		_("Trace test"), _("p-value"),
		_("Lmax test"), _("p-value"));

	for (i=0; i<hmax; i++) {
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
	    jv->A = TmpR;
	    TmpR = NULL;
	    alpha = compute_alpha(jv, jv->A, hmax);
	    print_beta_and_alpha(jv, jv->A, alpha, evals, hmax, pdinfo, prn);
	    compute_long_run_matrix(jv, evals, hmax, pdinfo, prn);
	} else {
	    jv->A = select_eigvecs(jv, TmpR, evals, hmax);
	    alpha = compute_alpha(jv, jv->A, hmax);
	    print_beta_and_alpha(jv, jv->A, alpha, evals, hmax, pdinfo, prn);
	    compute_vecm(jv, alpha, pdinfo, prn);
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
