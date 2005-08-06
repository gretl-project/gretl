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

#if DF_ADJ
/* for comparison with Eviews */
double ll_adj;
#endif

static int 
johansen_normalize (JVAR *jv, gretl_matrix *evecs)
{
    gretl_matrix *a = NULL, *b = NULL;
    double x, den;
    int nv = gretl_matrix_rows(jv->Svv);
    int i, j, h;
    int err = 0;

    if (jv->rank > 0) {
	h = jv->rank;
    } else {
	h = nv;
    }

    a = gretl_column_vector_alloc(nv);
    b = gretl_column_vector_alloc(nv);

    if (a == NULL || b == NULL) {
	gretl_matrix_free(a);
	gretl_matrix_free(b);
	return E_ALLOC;
    }

    for (j=0; j<h; j++) {
	/* select column from evecs */
	for (i=0; i<nv; i++) {
	    x = gretl_matrix_get(evecs, i, j);
	    gretl_vector_set(a, i, x);
	}

	/* find value of appropriate denominator */
	gretl_matrix_multiply(jv->Svv, a, b);
	den = gretl_vector_dot_product(a, b, &err);

	if (!err) {
	    den = sqrt(den);
	    for (i=0; i<nv; i++) {
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

/* VECM: print the cointegrating vectors, up to the specified
   rank */

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
	    /* re-scale */
	    r = gretl_matrix_get(jv->Beta, j, j);
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(jv->Beta, i, j) / r);
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* print cointegrating vectors or adjustments, either "raw" or
   re-scaled */

static void print_beta_or_alpha (JVAR *jv, int k,
				 const DATAINFO *pdinfo, PRN *prn,
				 int which, int rescale)
{
    gretl_matrix *c = (which == PRINT_BETA)? jv->Beta : jv->Alpha;
    int rows = gretl_matrix_rows(c);
    int i, j;
    double x;

    if (rescale) {
	pprintf(prn, "\n%s\n", (which == PRINT_BETA)? 
		_("renormalized beta") :
		_("renormalized alpha"));
    } else {
	pprintf(prn, "\n%s\n", (which == PRINT_BETA)? 
		_("beta (cointegrating vectors)") : 
		_("alpha (adjustment vectors)"));
    }

    for (i=0; i<rows; i++) {
	if (i < jv->list[0]) {
	    pprintf(prn, "%-10s", pdinfo->varname[jv->list[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    pprintf(prn, "%-10s", "const");
	} else if (jv->code == J_REST_TREND) {
	    pprintf(prn, "%-10s", "trend");
	}
	for (j=0; j<k; j++) {
	    if (rescale) {
		x = gretl_matrix_get(jv->Beta, j, j);
		if (which == PRINT_BETA) {
		    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j) / x);
		} else {
		    pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j) * x);
		}
	    } else {
		pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
	    }
	}
	pputc(prn, '\n');
    }
}

#if 0
static void rescale_alpha (JVAR *jv)
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

static void rescale_beta (JVAR *jv)
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

/* Calculate alpha (adjustments) matrix as per Johansen, 1991, eqn
   2.8, p. 1554.
*/

static int compute_alpha (JVAR *jv)
{
    gretl_matrix *alpha = NULL;
    gretl_matrix *tmp1 = NULL;
    gretl_matrix *tmp2 = NULL;
    int nv = gretl_matrix_rows(jv->Svv);
    int h, n = jv->neqns;
    int err = 0;

    if (jv->rank > 0) {
	h = jv->rank;
    } else {
	h = nv;
    }

    tmp1 = gretl_matrix_alloc(nv, h);
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

/* For cointegration test, print the "long-run matrix," \alpha \beta'
   in Johansen's notation or \Zeta_0 in Hamilton's.
*/

static void print_lr_matrix (JVAR *jv, gretl_matrix *Zeta_0,
			     const DATAINFO *pdinfo, PRN *prn)
{
    int cols = gretl_matrix_rows(jv->Svv);
    int i, j;

    pprintf(prn, "%s\n", _("long-run matrix (alpha * beta')"));
    pprintf(prn, "%22s", pdinfo->varname[jv->list[1]]);

    for (j=1; j<=cols; j++) {
	if (j < jv->list[0]) {
	    pprintf(prn, "%13s", pdinfo->varname[jv->list[j+1]]);
	} else if (jv->code == J_REST_CONST) {
	    pprintf(prn, "%13s", "const");
	} else if (jv->code == J_REST_TREND) {
	    pprintf(prn, "%13s", "trend");
	}
    }

    pputc(prn, '\n');

    for (i=0; i<jv->neqns; i++) {
	pprintf(prn, "%-10s", pdinfo->varname[jv->list[i+1]]);
	for (j=0; j<cols; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(Zeta_0, i, j));
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* Calculate the variance of the AR coefficients and the "adjustments"
   (alpha) in the VECM, as per Johansen, 1991, Appendix C, p. 1574.
   We don't do the full Kronecker thing since we just need the
   diagonal elements to compute standard errors.
*/

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
    int nc = k - (n - h); /* reduced rank */
    int i, j;
    int err = 0;

    X = gretl_matrix_alloc(T, n);
    tmp = gretl_matrix_alloc(T, h);
    D = gretl_matrix_alloc(T, nc);

    if (X == NULL || tmp == NULL || D == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_print(jv->Data, "jv->Data", NULL);

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
	    x = gretl_matrix_get(tmp, i, j - pn);
	    gretl_matrix_set(D, i, j, x);
	}
    }

    /* add deterministic vars to D */
    for (j=pn+h; j<nc; j++) {
	for (i=0; i<T; i++) {
	    x = gretl_matrix_get(jv->Data, i, j + (n - h));	
	    gretl_matrix_set(D, i, j, x);
	}
    }

    V = gretl_matrix_vcv(D);

    gretl_matrix_print(D, "D", NULL);
    gretl_matrix_print(V, "V", NULL);
    

    if (V == NULL) {
	err = 1;
    } else {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	int se_rows = nc + (jv->code >= J_UNREST_CONST);
	int ndet = nc - (pn + h);
	double a, vii, se;

	/* matrix to hold standard errors */
	jv->Ase = gretl_matrix_alloc(se_rows, n);
	if (jv->Ase == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	gretl_matrix_zero(jv->Ase);

	/* get the standard errors off the diagonal; shift the std
	   errs of the adjustments into the leading places in Ase, and
	   scale them using beta (FIXME?) */
	for (j=0; j<n; j++) {
	    a = gretl_matrix_get(jv->Omega, j, j) / T;
	    k = jv->rank;
	    for (i=0; i<pn+h; i++) { 
		vii = gretl_matrix_get(V, i, i);
		se = sqrt(vii * a);
		if (k < jv->rank) {
		    se *= fabs(gretl_matrix_get(jv->Beta, k, k));
		} 
		gretl_matrix_set(jv->Ase, k, j, se);
		k = (i < pn - 1)? k + 1 : i - pn + 1;
	    }
	}

	/* deterministic vars: FIXME this is wrong */
	if (ndet > 0) {
	    for (j=0; j<n; j++) {
		a = gretl_matrix_get(jv->Omega, j, j) / T;
		k = pn + h + (jv->code >= J_UNREST_CONST);
		for (i=0; i<ndet; i++) {
		    vii = 2.0 * gretl_matrix_get(jv->Omega, i, i);
		    se = sqrt(vii * a);
#if 0 /* temporary */
		    gretl_matrix_set(jv->Ase, k++, j, se);
#else
		    gretl_matrix_set(jv->Ase, k++, j, 0.0); 
#endif
		}
	    }
	}
    }

 bailout:

    gretl_matrix_free(D);
    gretl_matrix_free(X);
    gretl_matrix_free(tmp);
    gretl_matrix_free(V);

    return err;
}

/* Compute Hamilton's Omega (Johansen 1991 calls it Lambda): the
   cross-equation variance matrix.  This is not yet set up to
   handle the restricted cases.
*/

static int compute_omega (JVAR *jv, const gretl_matrix *Z0)
{
    gretl_matrix *uhat = NULL;
    gretl_matrix *omega = NULL;
    gretl_matrix *tmp = NULL;
    int n = jv->neqns;
    int T = jv->t2 - jv->t1 + 1;
#if DF_ADJ
    int df_adj = 0; /* not really wanted? */
#endif
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
#if DF_ADJ
	if (df_adj) {
	    /* Eviews 4.0 applies df adjustment here */
	    den -= jv->nparam;
	}
#endif

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

/* VECM: print the cross-equation variance matrix */

static int print_omega (JVAR *jv, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *M;
    char s[32];
    int i, j;

    pprintf(prn, "%s\n\n", _("Cross-equation covariance matrix"));

    for (i=0; i<jv->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[jv->list[i+1]]);
	if (i == 0) {
	    pprintf(prn, "%25s", s);
	} else {
	    pprintf(prn, "%13s", s);
	}
    }
    pputc(prn, '\n');

    for (i=0; i<jv->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[jv->list[i+1]]);
	pprintf(prn, "%-13s", s);
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

/* Fetch the Total Sum of Squares for the first difference of a given
   variable */

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

/* VECM: print the results of estimation */

static int 
print_vecm (JVAR *jv, const double **Z, const DATAINFO *pdinfo, PRN *prn)
{
    char s[32];
    double *ess = NULL;
    double *xd = NULL;
    double x;
    int constrow = 0;
    int trendrow = 0;
    int T = jv->t2 - jv->t1 + 1;
    int rows = gretl_matrix_rows(jv->A);
    int p = jv->order - 1;
    int n = jv->neqns;
    int h = jv->rank;
    int d = 0;
    int i, j, k, m;

    ess = malloc(jv->neqns * sizeof *ess);
    if (ess == NULL) {
	return E_ALLOC;
    }

    xd = malloc(T * sizeof *xd);
    if (xd == NULL) {
	free(ess);
	return E_ALLOC;
    }   

    if (jv->code >= J_UNREST_CONST) {
	constrow = h + p * n;
    }
    
    if (jv->code == J_UNREST_TREND) {
	trendrow = rows - 1;
    }


    pprintf(prn, "%s\n", _("Error correction equations"));
    pprintf(prn, "%s\n\n", _("standard errors in (); t-ratios in []"));

    for (i=0; i<jv->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[jv->list[i+1]]);
	if (i == 0) {
	    pprintf(prn, "%25s", s);
	} else {
	    pprintf(prn, "%13s", s);
	}
    }
    pputs(prn, "\n\n");

    m = k = 1;

    for (i=0; i<rows; i++) {
	if (i < h) {
	    sprintf(s, "CV%d", i + 1);
	    pprintf(prn, "%-13s", s);
	} else if (k <= jv->list[0]) {
	    sprintf(s, "d_%s(-%d)", pdinfo->varname[jv->list[k]], m++);
	    pprintf(prn, "%-13s", s);
	} else if (i == constrow) {
	    pprintf(prn, "%-13s", "const");
	} else if (i == trendrow) {
	    pprintf(prn, "%-13s", "trend");
	} else {
	    sprintf(s, "S%d", ++d);
	    pprintf(prn, "%-13s", s);
	}

	/* print coefficients */
	for (j=0; j<n; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(jv->A, i, j));
	}
	pputc(prn, '\n');

	/* print standard errors and t-ratios if known (unknown std
	   errs are currently set to zero)
	*/
	x = gretl_matrix_get(jv->Ase, i, 0);
	if (x == 0.0) {
	    pputc(prn, '\n');
	    continue;
	}
	pprintf(prn, "%-13s", " ");
	for (j=0; j<n; j++) {
	    sprintf(s, "(%#.5g)", gretl_matrix_get(jv->Ase, i, j));
	    pprintf(prn, " %12s", s);
	}
	pputc(prn, '\n');
	pprintf(prn, "%-13s", " ");
	for (j=0; j<n; j++) {
	    double tr =  gretl_matrix_get(jv->A, i, j)/
		gretl_matrix_get(jv->Ase, i, j);

	    sprintf(s, "[%#.4g]", tr);
	    pprintf(prn, " %12s", s);
	}	    
	pputs(prn, "\n\n");

	if (i > 0 && i % p == 0) {
	    k++;
	    m = 1;
	}
    }

    /* per-equation SSR */
    pprintf(prn, "%-13s", "SSR");
    for (i=0; i<n; i++) {
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
    pprintf(prn, "%-13s", "SE");
    for (i=0; i<n; i++) {
	x = gretl_matrix_get(jv->Omega, i, i);
	pprintf(prn, "%#12.5g ", sqrt(x));
    }

    /* per-equation R-squared */
    pputc(prn, '\n');
    pprintf(prn, "%-13s", "R-squared");
    for (i=0; i<n; i++) {
	double tss = get_tss(jv, xd, Z, i);

	pprintf(prn, "%#12.5g ", 1.0 - ess[i] / tss);
    }
    pputs(prn, "\n\n");    

    print_omega(jv, pdinfo, prn);

    if (!na(jv->ll)) {
	int k = gretl_matrix_rows(jv->A) + jv->rank; /* is this right?? */
	double aic, bic;

	aic = (-2.0 * jv->ll + 2.0 * k * n) / T;
	bic = (-2.0 * jv->ll + log(T) * k * n) / T;

	pprintf(prn, "\nlog-likelihood = %g\n", jv->ll);
	pprintf(prn, "AIC = %g\n", aic);
	pprintf(prn, "BIC = %g\n", bic);
#if DF_ADJ
	aic = (-2.0 * ll_adj + 2.0 * k * n) / T;
	bic = (-2.0 * ll_adj + log(T) * k * n) / T;

	pputs(prn, "\nd.f. adjusted\n");
	pprintf(prn, "log-likelihood = %g\n", ll_adj);
	pprintf(prn, "AIC = %g\n", aic);
	pprintf(prn, "BIC = %g\n", bic);
#endif
    }

    free(ess);
    free(xd);

    return 0;
}

/* restricted case: copy out first column into \mu, and
   copy the square remainder of the matriz into \Zeta_0 
*/

static gretl_matrix *get_real_zeta_0 (JVAR *jv, gretl_matrix *C)
{
    gretl_matrix *Z0;
    int n = jv->neqns;
    int i, j;

    Z0 = gretl_matrix_alloc(n, n);
    if (Z0 == NULL) {
	Z0 = C;
    } else {
	double x;

	for (i=0; i<n; i++) {
	    x = gretl_matrix_get(C, i, 0);
	    gretl_vector_set(jv->mu, i, x);
	}
	for (j=0; j<n; j++) {
	    for (i=0; i<n; i++) {
		x = gretl_matrix_get(C, i, j + 1);
		gretl_matrix_set(Z0, i, j, x);
	    }
	}
	gretl_matrix_free(C);
    }

    return Z0;
}

/* main driver for computing VECM results */

static int 
compute_vecm (JVAR *jv, const double **Z, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *Z0 = NULL;
    gretl_matrix *Zeta = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *det = NULL;
    int n = jv->neqns;
    int p = jv->order - 1;
    int nv = gretl_matrix_cols(jv->Suv);
    int nc = n * p + jv->rank;
    int i, err = 0;

    nc += (jv->code >= J_UNREST_CONST) + jv->nseas + 
	(jv->code == J_UNREST_TREND);

    /* full-size coefficient matrix */
    jv->A = gretl_matrix_alloc(nc, n);
    if (jv->A == NULL) {
	return E_ALLOC;
    }

    Z0 = gretl_matrix_alloc(n, nv);
    Zeta = gretl_matrix_alloc(n, n);
    tmp = gretl_matrix_alloc(nv, nv);

    /* handle intercept */
    jv->mu = gretl_column_vector_alloc(n);
    if (jv->mu == NULL) {
	err = E_ALLOC;
    }

    /* handle seasonals and/or trend */
    if (!err && (jv->code == J_UNREST_TREND || jv->nseas > 0)) {
	det = gretl_column_vector_alloc(n);
	if (det == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && (Z0 == NULL || Zeta == NULL || tmp == NULL)) {
	err = E_ALLOC;
    } else {
	int j, k, arow = 0;
	double r, x;

	gretl_matrix_multiply_mod(jv->Beta, GRETL_MOD_NONE,
				  jv->Beta, GRETL_MOD_TRANSPOSE,
				  tmp);

	/* \Zeta_0 = \Sigma_{uv} A A' (Hamilton eq. 20.2.12) */
	gretl_matrix_multiply(jv->Suv, tmp, Z0);

	if (nv > n) {
	    /* restricted case */
	    Z0 = get_real_zeta_0(jv, Z0);
	}

	gretl_matrix_zero(jv->A);

	/* put adjustment (alpha) terms in first */
	for (i=0; i<jv->rank; i++) {
	    for (j=0; j<n; j++) {
		r = gretl_matrix_get(jv->Beta, i, i);
		x = gretl_matrix_get(jv->Alpha, j, i);
		gretl_matrix_set(jv->A, i, j, x * r);
	    }
	}

	/* then the coeffs on lagged \Delta X */
	for (i=0; i<p; i++) {
	    gretl_matrix_copy_values(Zeta, jv->Pi[i]);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(Zeta, tmp);
	    if (jv->code == J_REST_CONST) {
		gretl_matrix_multiply_mod(jv->mu, GRETL_MOD_NONE,
					  jv->Aux[i], GRETL_MOD_TRANSPOSE,
					  tmp);
		gretl_matrix_subtract_from(Zeta, tmp);
	    }
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

	arow = jv->rank + p * n;

	/* it's all vectors below */
	gretl_matrix_reuse(tmp, n, 1);

	/* constants (Johansen's \mu, Hamilton's \alpha),
	   as per Hamilton, eq. 20.2.14, p. 637 */
	if (jv->code >= J_UNREST_CONST) {
	    gretl_matrix_copy_values(jv->mu, jv->Pi[i]);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(jv->mu, tmp);
#if 0 /* broken */
	    coint_consts_from_mu(jv->mu, jv->Alpha);
#endif
	    for (j=0; j<n; j++) {
		x = gretl_vector_get(jv->mu, j);
		gretl_matrix_set(jv->A, arow, j, x);
	    }
	    arow++;
	    i++;
	}

	/* FIXME when seasonals are included, both the consts
	   and the seasonals are way off */

	/* seasonals, if any */
	for (k=0; k<jv->nseas; k++) {
	    gretl_matrix_copy_values(det, jv->Pi[i]);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(det, tmp);
	    for (j=0; j<n; j++) {
		x = gretl_vector_get(det, j);
		gretl_matrix_set(jv->A, arow, j, x);
	    }
	    arow++;
	    i++;
	}

	/* trend, if any */
	if (jv->code == J_UNREST_TREND) {
	    gretl_matrix_copy_values(det, jv->Pi[i]);
	    gretl_matrix_multiply(Z0, jv->Theta[i], tmp);
	    gretl_matrix_subtract_from(det, tmp);
	    for (j=0; j<n; j++) {
		x = gretl_vector_get(det, j);
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
    gretl_matrix_free(det);

    return err;
}

/* Calculate and print the "long-run matrix", \alpha \beta', (what
   Hamilton calls \Zeta_0) */

static int 
compute_long_run_matrix (JVAR *jv, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *Z0 = NULL;
    gretl_matrix *tmp = NULL;
    int nv = gretl_matrix_rows(jv->Svv);
    int n = jv->neqns;
    int err = 0;

    Z0 = gretl_matrix_alloc(n, nv); 
    tmp = gretl_matrix_alloc(nv, nv);

    if (Z0 == NULL || tmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	gretl_matrix_multiply_mod(jv->Beta, GRETL_MOD_NONE,
				  jv->Beta, GRETL_MOD_TRANSPOSE,
				  tmp);

	/* \Zeta_0 = \Sigma_{uv} A A' */
	gretl_matrix_multiply(jv->Suv, tmp, Z0);
	print_lr_matrix(jv, Z0, pdinfo, prn);
    }

    gretl_matrix_free(Z0);
    gretl_matrix_free(tmp);
    
    return err;
}

/* Print both "raw" and re-scaled versions of the beta and alpha
   matrices (cointegrating vectors and vectors of adjustments
   respectively).
*/

static int
print_beta_and_alpha (JVAR *jv, double *evals, int h,
		      const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    pprintf(prn, "\n%s", _("eigenvalue"));
    for (i=0; i<h; i++) {
	pprintf(prn, "%#12.5g ", evals[i]);
    }
    pputc(prn, '\n');

    /* "raw" vectors */
    print_beta_or_alpha(jv, h, pdinfo, prn, PRINT_BETA, 0);
    print_beta_or_alpha(jv, h, pdinfo, prn, PRINT_ALPHA, 0);

    /* re-scaled versions */
    print_beta_or_alpha(jv, h, pdinfo, prn, PRINT_BETA, 1);
    print_beta_or_alpha(jv, h, pdinfo, prn, PRINT_ALPHA, 1);

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

/* Calculate the "base" component of the log-likelihood for a VECM,
   prior to adding the component associated with the eigenvalues.
   Fiddle around with df adjustment for comparison with Eviews.
*/

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
	ldet = gretl_matrix_log_determinant(Suu);
	jv->ll = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;

#if DF_ADJ
	/* "df-adjusted ll" for comparison with Eviews */
	gretl_matrix_copy_values(Suu, jv->Suu);
	gretl_matrix_multiply_by_scalar(Suu, (double) T / (T - jv->nparam));
	ldet = gretl_matrix_log_determinant(Suu);
	ll_adj = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
#endif

	gretl_matrix_free(Suu);
    }

    return err;
}

/* Public entry point for both cointegration test and VECM
   estimation */

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

    /* We want to preserve Svv and Suu, so copy them before
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
	/* "re-expand" this matrix */
	gretl_matrix_reuse(TmpR, kv, kv);
    }

    eigvals = gretl_general_matrix_eigenvals(M, TmpR);
    if (eigvals == NULL) {
	pputs(prn, _("Failed to find eigenvalues\n"));
	err = E_ALLOC;
    } else {
	err = gretl_eigen_sort(eigvals, TmpR, jv->rank);
    }

    if (!err) {
	int T = jv->t2 - jv->t1 + 1;
	double cumeig = 0.0;
	double *lambdamax = NULL, *trace = NULL;
	double pval[2];
	int hmax, i;

	/* in case kv > k */
	hmax = k;
	k = kv;

	if (jv->rank > 0) {
	    hmax = jv->rank;
	}

	trace = malloc(hmax * sizeof *trace);
	lambdamax = malloc(hmax * sizeof *lambdamax);
	if (trace == NULL || lambdamax == NULL) {
	    free(trace);
	    free(lambdamax);
	    err = 1;
	    goto eigenvals_bailout;
	}

	for (i=hmax-1; i>=0; i--){
      	    lambdamax[i] = -T * log(1.0 - eigvals[i]); 
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
#if DF_ADJ
		ll_adj += lambdamax[i] / 2.0;
#endif
	    }
	    gamma_par_asymp(trace[i], lambdamax[i], jv->code, hmax - i, pval);
	    pprintf(prn, "%4d%#11.5g%#11.5g [%s]%#11.5g [%s]\n", \
		    i, eigvals[i], trace[i], safe_print_pval(pval[0], 0), 
		    lambdamax[i], safe_print_pval(pval[1], 1));
	}
	pputc(prn, '\n');

	/* normalize the eigenvectors and compute adjustments */
	johansen_normalize(jv, TmpR);
	jv->Beta = TmpR;
	TmpR = NULL;
	err = compute_alpha(jv);

	if (!err) {
	    if (jv->rank == 0) {
		print_beta_and_alpha(jv, eigvals, hmax, pdinfo, prn);
		compute_long_run_matrix(jv, pdinfo, prn);
	    } else {
		print_coint_eqns(jv, pdinfo, prn);
		compute_vecm(jv, Z, pdinfo, prn);
	    }
	}

	free(eigvals);
	free(lambdamax);
	free(trace);
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
