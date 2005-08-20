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

#define JDEBUG 0

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

/* normalize the eigenvectors based on the Svv matrix */

static int johansen_normalize (JohansenInfo *jv, gretl_matrix *evecs)
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

/* for cointegration test: print cointegrating vectors or adjustments,
   either "raw" or re-scaled */

static void print_beta_or_alpha (JohansenInfo *jv, int k,
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

/* Calculate alpha (adjustments) matrix as per Johansen, 1991, eqn
   2.8, p. 1554.
*/

static int compute_alpha (JohansenInfo *jv, int n)
{
    gretl_matrix *alpha = NULL;
    gretl_matrix *tmp1 = NULL;
    gretl_matrix *tmp2 = NULL;
    int h, nv = gretl_matrix_rows(jv->Svv);
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

static void print_lr_matrix (JohansenInfo *jv, gretl_matrix *Zeta_0,
			     int neqns, const DATAINFO *pdinfo, PRN *prn)
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

    for (i=0; i<neqns; i++) {
	pprintf(prn, "%-10s", pdinfo->varname[jv->list[i+1]]);
	for (j=0; j<cols; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(Zeta_0, i, j));
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* Restricted const or trend: copy out last column into \rho, and copy
   the square remainder of the matrix into \Zeta_0.
*/

static gretl_matrix *
get_real_zeta_0 (gretl_matrix *C, gretl_matrix *rho, int n)
{
    gretl_matrix *Z0 = gretl_matrix_alloc(n, n);
    int i, j;

    if (Z0 == NULL) {
	Z0 = C;
    } else {
	double x;

	for (i=0; i<n; i++) {
	    x = gretl_matrix_get(C, i, n);
	    gretl_vector_set(rho, i, x);
	}

	for (j=0; j<n; j++) {
	    for (i=0; i<n; i++) {
		x = gretl_matrix_get(C, i, j);
		gretl_matrix_set(Z0, i, j, x);
	    }
	}

	gretl_matrix_free(C);
    }

    return Z0;
}

/* Compute Hamilton's Omega (Johansen 1991 calls it Lambda): the
   cross-equation variance matrix.  Uses Hamilton's eq. 20.2.15,
   p. 638, or equation at top of p. 645 in restricted const case.
*/

static int compute_omega (GRETL_VAR *vecm)
{
    gretl_matrix *Z0 = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *uhat = NULL;
    gretl_matrix *omega = NULL;
    gretl_matrix *rho = NULL;
    int n = vecm->neqns;
    int nv = gretl_matrix_cols(vecm->jinfo->Suv);
    int T = vecm->T;
    int err = 0;

    Z0 = gretl_matrix_alloc(n, nv);
    tmp = gretl_matrix_alloc(nv, nv);

    if (Z0 == NULL || tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (restricted(vecm)) {
	rho = gretl_column_vector_alloc(n);
	if (rho == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    gretl_matrix_multiply_mod(vecm->jinfo->Beta, GRETL_MOD_NONE,
			      vecm->jinfo->Beta, GRETL_MOD_TRANSPOSE,
			      tmp);

    /* \Zeta_0 = \Sigma_{uv} A A' (Hamilton eq. 20.2.12) */
    gretl_matrix_multiply(vecm->jinfo->Suv, tmp, Z0);

    if (restricted(vecm)) {
	/* in the restricted cases "Z0" contains an extra column,
	   which has to be separated out */
	Z0 = get_real_zeta_0(Z0, rho, n);
    }

    uhat = gretl_matrix_alloc(n, T);
    omega = gretl_matrix_alloc(n, n);

    gretl_matrix_free(tmp);
    tmp = gretl_matrix_alloc(n, T);

    if (uhat == NULL || omega == NULL || tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_copy_values(uhat, vecm->jinfo->u);
    gretl_matrix_multiply(Z0, vecm->jinfo->v, tmp);
    gretl_matrix_subtract_from(uhat, tmp);

    if (restricted(vecm)) {
	gretl_matrix_multiply(rho, vecm->jinfo->w, tmp);
	gretl_matrix_subtract_from(uhat, tmp);
    }

    gretl_matrix_multiply_mod(uhat, GRETL_MOD_NONE,
			      uhat, GRETL_MOD_TRANSPOSE,
			      omega);
    gretl_matrix_divide_by_scalar(omega, T);

 bailout:

    if (!err) {
	vecm->S = omega;
    } else {
	gretl_matrix_free(omega);
    }

    gretl_matrix_free(Z0);
    gretl_matrix_free(tmp);
    gretl_matrix_free(uhat);
    gretl_matrix_free(rho);

    return err;
}

/* compute the EC terms and add them to the dataset, Z, so we
   can run OLS */

static int 
add_EC_terms_to_dataset (GRETL_VAR *vecm, double ***pZ, DATAINFO *pdinfo)
{
    const gretl_matrix *Beta = vecm->jinfo->Beta;
    int rank = jrank(vecm);
    int *list = vecm->jinfo->list;
    double xt, bxt, sb;
    int i, j, t, v = pdinfo->v;
    int id = gretl_VECM_id(vecm);
    int err;

    err = dataset_add_series(rank, pZ, pdinfo);

    if (!err) {
	for (j=0; j<rank; j++) {
	    sprintf(pdinfo->varname[v+j], "CV%d", j + 1);
	    sprintf(VARLABEL(pdinfo, v+j), "cointegrating vector %d from VECM %d", 
		    j + 1, id);
	    for (t=0; t<pdinfo->n; t++) {
		if (t < vecm->t1 || t > vecm->t2) {
		    (*pZ)[v+j][t] = NADBL;
		} else { 
		    bxt = 0.0;
		    /* beta * X(t-1) */
		    for (i=0; i<vecm->neqns; i++) {
			xt = (*pZ)[list[i+1]][t-1];
			sb = gretl_matrix_get(Beta, i, j);
			sb /= gretl_matrix_get(Beta, j, j);
			bxt += sb * xt;
		    }
		    /* restricted const or trend */
		    if (restricted(vecm)) {
			sb = gretl_matrix_get(Beta, i, j);
			sb /= gretl_matrix_get(Beta, j, j);
			if (jcode(vecm) == J_REST_TREND) {
			    sb *= t;
			}
			bxt += sb;
		    }
		    (*pZ)[v+j][t] = bxt;
		}
	    }
	}
    }
	
    return err;
}

/* Run OLS taking the betas, as calculated via the eigen-analysis,
   as given.  So obtain estimates and standard errors for the
   coefficients on the lagged differences and the unrestricted
   deterministic vars.  Construct full residuals matrix while
   we're at it.

   FIXME: when seasonals are included, we're not getting the
   same results as JMulTi for the constant (though the results
   are the same for the seasonal dummies themselves).
*/

static int build_VECM_models (GRETL_VAR *vecm, double ***pZ, DATAINFO *pdinfo,
			      PRN *prn)
{
    int i, j, k, v = pdinfo->v;
    int mt, t, r = vecm->jinfo->rank;
    int *biglist = vecm->jinfo->biglist;
    int err = 0;

    vecm->E = gretl_matrix_alloc(vecm->T, vecm->neqns);

    err = add_EC_terms_to_dataset(vecm, pZ, pdinfo);

    for (i=0; i<vecm->neqns && !err; i++) {
	biglist[1] = vecm->jinfo->difflist[i+1];
	k = biglist[0] - r + 1;
	for (j=0; j<r; j++) {
	    biglist[k++] = v + j;
	}
	*vecm->models[i] = lsq(biglist, pZ, pdinfo, OLS, OPT_N | OPT_Z, 0.0);
	err = vecm->models[i]->errcode;
	if (!err) {
	    vecm->models[i]->ID = i + 1;
	    vecm->models[i]->aux = AUX_VECM;
	    vecm->models[i]->adjrsq = NADBL;
	    if (vecm->E != NULL) {
		for (t=0; t<vecm->T; t++) {
		    mt = t + vecm->t1;
		    gretl_matrix_set(vecm->E, t, i, vecm->models[i]->uhat[mt]);
		}
	    }
	}
    }

    return err;
}

/* Calculate and print the "long-run matrix", \alpha \beta', (what
   Hamilton calls \Zeta_0) */

static int 
compute_long_run_matrix (JohansenInfo *jv, int n, const DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *Z0 = NULL;
    gretl_matrix *tmp = NULL;
    int nv = gretl_matrix_rows(jv->Svv);
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
	print_lr_matrix(jv, Z0, n, pdinfo, prn);
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
print_beta_and_alpha (JohansenInfo *jv, double *evals, int h,
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

/* Calculate the "base" component of the log-likelihood for a VECM,
   prior to adding the component associated with the eigenvalues.
*/

static int johansen_ll_init (GRETL_VAR *vecm)
{
    gretl_matrix *Suu;
    double T_2 = (double) vecm->T / 2.0;
    double ldet, n = vecm->neqns;
    int err = 0;

    Suu = gretl_matrix_copy(vecm->jinfo->Suu);

    if (Suu == NULL) {
	err = E_ALLOC;
    } else {
	ldet = gretl_matrix_log_determinant(Suu);
	vecm->ll = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
	gretl_matrix_free(Suu);
    }

    return err;
}

#if NOTYET 
/* Johansen (1995, p. 184):

   T(I - \beta c')S_{11}^{-1}(I - c \beta') Kr 
       (\alpha' \Omega^{-1} \a)^{-1}

   The following is not right yet: I haven't got the normalization
   correct, I think.

*/

static int beta_variance (GRETL_VAR *vecm)
{
    gretl_matrix *O;
    gretl_matrix *aOa;
    gretl_matrix *Svv;
    gretl_matrix *tmp;
    gretl_matrix *LH;
    gretl_matrix *alpha;
    gretl_matrix *beta;
    gretl_matrix *c, *IB1, *IB2;
    int r = jrank(vecm);
    int n = vecm->neqns; /* ?? */
    int i, j, err = 0;

    tmp = gretl_matrix_alloc(n, n);
    aOa = gretl_matrix_alloc(r, r);

    alpha = gretl_matrix_copy(vecm->jinfo->Alpha);
    for (i=0; i<n; i++) {
	for (j=0; j<r; j++) {
	    double x = gretl_matrix_get(vecm->jinfo->Beta, j, j);
	    
	    gretl_matrix_set(alpha, i, j, gretl_matrix_get(vecm->jinfo->Alpha, i, j) * x);
	}
    }

    gretl_matrix_print(alpha, "alpha", NULL);

    /* compute right hand inverse matrix */
    O = gretl_matrix_copy(vecm->S);
    gretl_invert_symmetric_matrix(O);

    gretl_matrix_print(vecm->S, "jv->S", NULL);
    gretl_matrix_print(O, "O", NULL);

    gretl_matrix_reuse(tmp, r, n);
    gretl_matrix_multiply_mod(alpha, GRETL_MOD_TRANSPOSE,
			      O, GRETL_MOD_NONE,
			      tmp);
    gretl_matrix_multiply(tmp, alpha, aOa);
    gretl_invert_symmetric_matrix(aOa);

    gretl_matrix_print(aOa, "aOa", NULL);

    /* compute left-hand matrix */
    Svv = gretl_matrix_copy(vecm->jinfo->Svv);
    gretl_invert_symmetric_matrix(Svv);

    gretl_matrix_print(vecm->jinfo->Svv, "Svv", NULL);
    gretl_matrix_print(Svv, "Svv(-1)", NULL);

    beta = gretl_matrix_copy(vecm->jinfo->Beta);
    for (i=0; i<n; i++) {
	for (j=0; j<r; j++) {
	    double x = gretl_matrix_get(vecm->jinfo->Beta, j, j);

	    gretl_matrix_set(beta, i, j, gretl_matrix_get(vecm->jinfo->Beta, i, j) / x);
	}
    }

    c = gretl_matrix_alloc(n, r);
    gretl_matrix_zero(c);
    for (i=0; i<r; i++) {
	gretl_matrix_set(c, i, i, 1.0);
    }

    gretl_matrix_print(c, "c", NULL);

    IB1 = gretl_identity_matrix_new(n);
    IB2 = gretl_identity_matrix_new(n);
    LH = gretl_matrix_alloc(n, n); /* nv? */

    gretl_matrix_print(beta, "beta", NULL);

    gretl_matrix_reuse(tmp, n, n);
    gretl_matrix_multiply_mod(beta, GRETL_MOD_NONE,
			      c, GRETL_MOD_TRANSPOSE,
			      tmp);
    gretl_matrix_subtract_from(IB1, tmp);

    gretl_matrix_print(IB1, "IB1", NULL);

    gretl_matrix_multiply_mod(c, GRETL_MOD_NONE,
			      beta, GRETL_MOD_TRANSPOSE,
			      tmp);
    gretl_matrix_subtract_from(IB2, tmp);

    gretl_matrix_print(IB2, "IB2", NULL);

    gretl_matrix_multiply(Svv, IB2, tmp);
    gretl_matrix_multiply(IB1, tmp, LH);

    gretl_matrix_print(LH, "LH", NULL);

    for (i=0; i<n; i++) {
	double lh = gretl_matrix_get(LH, i, i);
	double rh = gretl_matrix_get(aOa, 0, 0);

	fprintf(stderr, "diag[%d] = %g\n", i+1, sqrt(lh * rh / vecm->T));
    }
    
    gretl_matrix_free(O);
    gretl_matrix_free(aOa);
    gretl_matrix_free(tmp);
    gretl_matrix_free(c);
    gretl_matrix_free(IB1);
    gretl_matrix_free(IB2);
    gretl_matrix_free(Svv);
    gretl_matrix_free(LH);
    gretl_matrix_free(beta);
    gretl_matrix_free(alpha);

    return err;
}
#endif /* NOTYET */

static int 
compute_coint_test (GRETL_VAR *jvar, const double *eigvals, PRN *prn)
{
    int T = jvar->T;
    int rank = jrank(jvar);
    double cumeig = 0.0;
    double *lambdamax = NULL;
    double *trace = NULL;
    double pvals[2];
    int h, i;

    if (rank > 0) {
	h = rank;
    } else {
	h = jvar->neqns;
    }

    trace = malloc(h * sizeof *trace);
    lambdamax = malloc(h * sizeof *lambdamax);

    if (trace == NULL || lambdamax == NULL) {
	free(trace);
	free(lambdamax);
	return E_ALLOC;
    }

    for (i=h-1; i>=0; i--){
	lambdamax[i] = -T * log(1.0 - eigvals[i]); 
	cumeig += lambdamax[i];
	trace[i] = cumeig; 
    }

    if (rank == 0) {
	/* just doing cointegration test */
	pputc(prn, '\n');
	print_Johansen_test_case(jcode(jvar), prn);
    }

    pprintf(prn, "\n%s %s %s %s   %s  %s\n", _("Rank"), _("Eigenvalue"), 
	    _("Trace test"), _("p-value"),
	    _("Lmax test"), _("p-value"));	

    for (i=0; i<h; i++) {
	if (!na(jvar->ll)) {
	    jvar->ll += lambdamax[i] / 2.0;
	}
	gamma_par_asymp(trace[i], lambdamax[i], jcode(jvar), h - i, pvals);
	pprintf(prn, "%4d%#11.5g%#11.5g [%6.4f]%#11.5g [%6.4f]\n", \
		i, eigvals[i], trace[i], pvals[0], lambdamax[i], pvals[1]);
    }
    pputc(prn, '\n');

    free(lambdamax);
    free(trace);

    return 0;
}

/* Public entry point for both cointegration test and VECM
   estimation */

int johansen_analysis (GRETL_VAR *jvar, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *M = NULL;
    gretl_matrix *TmpL = NULL, *TmpR = NULL;
    gretl_matrix *Suu = NULL, *Svv = NULL;

    double *eigvals = NULL;

    int n = jvar->neqns;
    int nv = gretl_matrix_cols(jvar->jinfo->Svv);
    int rank = jrank(jvar);
    int err = 0;

    TmpL = gretl_matrix_alloc(nv, n);
    TmpR = gretl_matrix_alloc(nv, nv);
    M = gretl_matrix_alloc(nv, nv);

    /* We want to preserve Svv and Suu, so copy them before
       inverting */
    Svv = gretl_matrix_copy(jvar->jinfo->Svv);
    Suu = gretl_matrix_copy(jvar->jinfo->Suu);

    if (Suu == NULL || Svv == NULL || TmpL == NULL || 
	TmpR == NULL || M == NULL) {
	err = 1;
	goto eigenvals_bailout;
    }

    if (nv > n) {
	gretl_matrix_reuse(TmpR, n, nv);
    }

    /* calculate Suu^{-1} Suv */
    err = gretl_invert_general_matrix(Suu);
    if (!err) {
	err = gretl_matrix_multiply(Suu, jvar->jinfo->Suv, TmpR);
    }

    /* calculate Svv^{-1} Suv' */
    if (!err) {
	err = gretl_invert_general_matrix(Svv);
    }
    if (!err) {
	err = gretl_matrix_multiply_mod(Svv, GRETL_MOD_NONE,
					jvar->jinfo->Suv, GRETL_MOD_TRANSPOSE, 
					TmpL);
    }

    if (!err) {
	err = gretl_matrix_multiply(TmpL, TmpR, M);
    }

    if (err) {
	goto eigenvals_bailout;
    }

    if (nv > n) {
	/* "re-expand" this matrix */
	gretl_matrix_reuse(TmpR, nv, nv);
    }

    eigvals = gretl_general_matrix_eigenvals(M, TmpR);
    if (eigvals == NULL) {
	pputs(prn, _("Failed to find eigenvalues\n"));
	err = E_ALLOC;
    } else {
	err = gretl_eigen_sort(eigvals, TmpR, rank);
    }

    if (!err) {
	johansen_ll_init(jvar);

	compute_coint_test(jvar, eigvals, prn);

	/* normalize the eigenvectors */
	johansen_normalize(jvar->jinfo, TmpR);
	jvar->jinfo->Beta = TmpR;
	TmpR = NULL;

	if (!err) {
	    if (rank == 0) {
		/* just running cointegration test */
		err = compute_alpha(jvar->jinfo, n);
		if (!err) {
		    print_beta_and_alpha(jvar->jinfo, eigvals, n, pdinfo, prn);
		    compute_long_run_matrix(jvar->jinfo, n, pdinfo, prn);
		}
	    } else {
		/* estimating VECM */
		err = build_VECM_models(jvar, pZ, pdinfo, prn);
		if (!err) {
		    err = compute_omega(jvar);
#if NOTYET
		    /* testing */
		    compute_alpha(jvar->jinfo, n);
		    beta_variance(jvar);
#endif
		}
	    }
	}
    } 

 eigenvals_bailout:    

    gretl_matrix_free(TmpL);
    gretl_matrix_free(TmpR);
    gretl_matrix_free(M);
    gretl_matrix_free(Svv);
    gretl_matrix_free(Suu);

    free(eigvals);

    return err;
}
