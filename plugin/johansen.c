/*
 *  Copyright (c) by Allin Cottrell and Riccardo "Jack" Lucchetti
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
#include "var.h"

#define JDEBUG 0
#define OLDEIG 0

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

    g = gamma_cdf_comp(mt, vt, tracetest, 2);
    if (na(g)) {
	pval[0] = NADBL;
    } else {
	pval[0] = 1.0 - g;
	if (pval[0] < 0.0) {
	    pval[0] = 0.0;
	}
    }

    g = gamma_cdf_comp(ml, vl, lmaxtest, 2);
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

/* Remove a possible excess zero from the end of a floating point
   number printed to the given precision p: are we working around
   a bug in the C library?
*/

static void fix_xstr (char *s, int p)
{
    int n = strlen(s);

    if (n > p && strspn(s + n - p, "0") == p) {
	s[n-1] = 0;
    }
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
    char xstr[32];
    int i, j;
    double x, y;

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
	    x = gretl_matrix_get(c, i, j);
	    if (rescale) {
		y = gretl_matrix_get(jv->Beta, j, j);
		if (which == PRINT_BETA) {
		    x /= y;
		} else {
		    x *= y;
		}
	    }
	    if (x == -0.0) {
		x = 0.0;
	    }
	    sprintf(xstr, "%#.5g", x);
	    fix_xstr(xstr, 5);
	    pprintf(prn, "%12s ", xstr);
	}
	pputc(prn, '\n');
    }
}

/* Calculate \alpha (adjustments) matrix as per Johansen, 1991, eqn
   2.8, p. 1554.  Required for the cointegration test, but not
   needed when doing a VECM (in which case we get \alpha via
   build_VECM_models() below).
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
	err = gretl_matrix_qform(jv->Beta, GRETL_MOD_TRANSPOSE, jv->Svv,
				 tmp2, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(tmp2);
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
   in Johansen's notation, or \Zeta_0 in Hamilton's.
*/

static void print_lr_matrix (JohansenInfo *jv, gretl_matrix *Zeta_0,
			     int neqns, const DATAINFO *pdinfo, PRN *prn)
{
    int cols = gretl_matrix_rows(jv->Svv);
    int i, j;

    pprintf(prn, "%s\n", _("long-run matrix (alpha * beta')"));

    pprintf(prn, "%22s", pdinfo->varname[jv->list[1]]); /* N.B. */
    for (j=2; j<=jv->list[0]; j++) {
	pprintf(prn, "%13s", pdinfo->varname[jv->list[j]]);
    }

    if (jv->code == J_REST_CONST) {
	pprintf(prn, "%13s", "const");
    } else if (jv->code == J_REST_TREND) {
	pprintf(prn, "%13s", "trend");
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

/* Compute Hamilton's Omega (Johansen 1991 calls it Lambda): the
   cross-equation variance matrix.
*/

static int compute_omega (GRETL_VAR *vecm)
{
    if (vecm->S == NULL) {
	vecm->S = gretl_matrix_alloc(vecm->neqns, vecm->neqns);
    }

    if (vecm->S == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(vecm->E, GRETL_MOD_TRANSPOSE,
			      vecm->E, GRETL_MOD_NONE,
			      vecm->S, GRETL_MOD_NONE);

    gretl_matrix_divide_by_scalar(vecm->S, vecm->T);
    
    return 0;
}

static void gretl_matrix_I (gretl_matrix *A, int n)
{
    int i;

    gretl_matrix_zero(A);
    for (i=0; i<n; i++) {
	gretl_matrix_set(A, i, i, 1.0);
    }
}

/* compute the EC terms and add them to the dataset, Z, so we
   can estimate the VECM by OLS (conditional on \beta) */

static int 
add_EC_terms_to_dataset (GRETL_VAR *vecm, double ***pZ, DATAINFO *pdinfo,
			 int iter)
{
    const gretl_matrix *B = vecm->jinfo->Beta;
    int rank = jrank(vecm);
    int *list = vecm->jinfo->list;
    double xt, bxt, sb;
    int i, j, t, v = pdinfo->v;
    int id = gretl_VECM_id(vecm);
    int vj, err = 0;

    if (iter == 0) {
	err = dataset_add_series(rank, pZ, pdinfo);
    }

    if (!err) {
	char vname[VNAMELEN];

	for (j=0; j<rank; j++) {
	    sprintf(vname, "EC%d", j + 1);

	    if (iter > 0) {
		/* series already allocated */
		vj = varindex(pdinfo, vname);
	    } else {
		vj = v + j;
		strcpy(pdinfo->varname[vj], vname);
		make_varname_unique(pdinfo->varname[vj], vj, pdinfo);
		sprintf(VARLABEL(pdinfo, vj), "error correction term %d from VECM %d", 
			j + 1, id);
	    }

	    for (t=0; t<pdinfo->n; t++) {
		if (t < vecm->t1 || t > vecm->t2) {
		    (*pZ)[vj][t] = NADBL;
		} else { 
		    bxt = 0.0;
		    /* beta * X(t-1) */
		    for (i=0; i<vecm->neqns; i++) {
			xt = (*pZ)[list[i+1]][t-1];
			sb = gretl_matrix_get(B, i, j);
			sb /= gretl_matrix_get(B, j, j);
			bxt += sb * xt;
		    }
		    /* restricted const or trend */
		    if (restricted(vecm)) {
			sb = gretl_matrix_get(B, i, j);
			sb /= gretl_matrix_get(B, j, j);
			if (jcode(vecm) == J_REST_TREND) {
			    sb *= t;
			}
			bxt += sb;
		    }
		    (*pZ)[vj][t] = bxt;
		}
	    }
	}
    }
	
    return err;
}

/* After doing OLS estimation of the VECM conditional on \beta: copy
   the coefficients on the lagged differences (i.e. form the \Gamma
   matrices) so we can compute the VAR representation */

static void copy_coeffs_to_Gamma (MODEL *pmod, int i, gretl_matrix **G,
				  int maxlag, int nv)
{
    int j, k, h;
    double x;
    
    for (k=0; k<maxlag; k++) {
	h = k + pmod->ifc;
	/* successive lags (distinct \Gamma_i matrices) */
	for (j=0; j<nv; j++) {
	    /* successive \Delta x_j */
	    x = pmod->coeff[h];
	    gretl_matrix_set(G[k], i, j, x);
	    h += maxlag;
	}
    }
}

/* Again, after doing OLS estimation of the VECM conditional on \beta:
   copy the coefficients on the EC terms (\beta' X) into the \alpha
   matrix. */

static void copy_coeffs_to_Alpha (GRETL_VAR *vecm, int i, gretl_matrix *Alpha,
				  int maxlag)
{
    double x;
    const MODEL *pmod = vecm->models[i];
    /* position in coeff array of first \alpha term */
    int base = vecm->jinfo->nexo + gretl_matrix_rows(Alpha) * maxlag;
    int j;

    for (j=0; j<vecm->jinfo->rank; j++) {
	x = pmod->coeff[base + j];
	gretl_matrix_set(Alpha, i, j, x);
    }
}

/* Form the matrix \Pi = \alpha \beta': since \beta is augmented
   in the case of restricted constant or restricted trend, we
   may have to make a reduced copy */

static int form_Pi (GRETL_VAR *vecm, const gretl_matrix *Alpha,
		    gretl_matrix *Pi)
{
    gretl_matrix *Beta = vecm->jinfo->Beta;
    int err = 0, freeit = 0;

    if (gretl_matrix_rows(Beta) > vecm->neqns) {
	Beta = gretl_matrix_alloc(vecm->neqns, vecm->jinfo->rank);
	if (Beta == NULL) {
	    err = E_ALLOC;
	} else {
	    double x;
	    int i, j;

	    for (i=0; i<vecm->neqns; i++) {
		for (j=0; j<vecm->jinfo->rank; j++) {
		    x = gretl_matrix_get(vecm->jinfo->Beta, i, j);
		    gretl_matrix_set(Beta, i, j, x);
		}
	    }
	    freeit = 1;
	}
    }

    if (!err) {
	gretl_matrix_multiply_mod(Alpha, GRETL_MOD_NONE,
				  Beta, GRETL_MOD_TRANSPOSE,
				  Pi, GRETL_MOD_NONE);
    }

    if (freeit) {
	gretl_matrix_free(Beta);
    }

    return err;
}

/* VAR representation: transcribe the coefficient matrix A_i (for lag
   i) into its place in the full VAR coefficient matrix, A */

static void add_Ai_to_VAR_A (gretl_matrix *Ai, GRETL_VAR *vecm, int k)
{
    int i, j, offset = k * vecm->neqns;
    double x;

    for (i=0; i<vecm->neqns; i++) {
	for (j=0; j<vecm->neqns; j++) {
	    x = gretl_matrix_get(Ai, i, j);
	    gretl_matrix_set(vecm->A, i, j + offset, x);
	}
    }
}

/* Run OLS taking the betas, as calculated via the eigen-analysis,
   as given.  So obtain estimates and standard errors for the
   coefficients on the lagged differences and the unrestricted
   deterministic vars.  Construct full residuals matrix while
   we're at it.

   FIXME: when seasonals are included, we're not getting the
   same results as JMulTi for the constant (though the results
   are the same for the seasonal dummies themselves)?
*/

static int 
build_VECM_models (GRETL_VAR *vecm, double ***pZ, DATAINFO *pdinfo, int iter)
{
    gretl_matrix *Pi = NULL;
    gretl_matrix *A = NULL;
    gretl_matrix **G = NULL;

    int rv0 = pdinfo->v;
    int mt, t, r = vecm->jinfo->rank;
    int p = vecm->order;
    int nv = vecm->neqns;
    int *biglist = vecm->jinfo->biglist;
    gretlopt lsqopt = OPT_N | OPT_Z;
    int i, j, k;
    int err = 0;

    /* Note: "vecm->order" is actually the order of the VAR system,
       which corresponds to the number of lagged differences on the
       RHS of the VAR system.  We need that number of G matrices to
       hold the coefficients on those lagged differences.
    */

#if JDEBUG
    fprintf(stderr, "build_VECM_models: vecm->order = %d\n", vecm->order);
#endif

    if (p < 0) {
	return E_DATA;
    }

    /* for computing VAR representation */
    Pi = gretl_matrix_alloc(nv, nv);
    A = gretl_matrix_alloc(nv, nv);
    if (Pi == NULL || A == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (p > 0) {
	G = gretl_matrix_array_alloc_with_size(p, nv, nv);
	if (G == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}	
    }    

    if (vecm->jinfo->Alpha == NULL) {
	vecm->jinfo->Alpha = gretl_matrix_alloc(nv, r);
	if (vecm->jinfo->Alpha == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (iter > 0) {
	/* bootstrapping: EC terms already in dataset */
	rv0 -= r;
	lsqopt |= OPT_A;
    }

    err = add_EC_terms_to_dataset(vecm, pZ, pdinfo, iter);

    for (i=0; i<nv && !err; i++) {
	biglist[1] = vecm->jinfo->difflist[i+1];
	k = biglist[0] - r + 1;
	for (j=0; j<r; j++) {
	    biglist[k++] = rv0 + j;
	}
#if JDEBUG
	printlist(biglist, "build_VECM_models: biglist");
#endif
	*vecm->models[i] = lsq(biglist, pZ, pdinfo, OLS, lsqopt);
	err = vecm->models[i]->errcode;
	if (!err) {
	    vecm->models[i]->ID = i + 1;
	    vecm->models[i]->aux = AUX_VECM;
	    vecm->models[i]->adjrsq = NADBL;
	    if (p > 0) {
		copy_coeffs_to_Gamma(vecm->models[i], i, G, p, nv);
	    }
	    copy_coeffs_to_Alpha(vecm, i, vecm->jinfo->Alpha, p);
	    for (t=0; t<vecm->T; t++) {
		mt = t + vecm->t1;
		gretl_matrix_set(vecm->E, t, i, vecm->models[i]->uhat[mt]);
	    }
	    if (i == 0) {
		vecm->ncoeff = vecm->models[i]->ncoeff;
	    }
	} else {
	    fprintf(stderr, "build_VECM_models: error %d from lsq, eqn %d, iter %d\n",
		    err, i + 1, iter);
	}
    }

    if (!err) {
	/* \Pi = \alpha \beta' */
	err = form_Pi(vecm, vecm->jinfo->Alpha, Pi);
    }

    if (err) {
	goto bailout;
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->Alpha, "Alpha from models");
    gretl_matrix_print(Pi, "Pi");
    for (i=0; i<p; i++) {
	fprintf(stderr, "Gamma matrix, lag %d\n\n", i+1);
	gretl_matrix_print(G[i], NULL);
    } 
#endif

    if (p == 0) {
	gretl_matrix_I(A, nv);
	gretl_matrix_add_to(A, Pi);
	add_Ai_to_VAR_A(A, vecm, 0);
    } else {
	for (i=0; i<=p; i++) {
	    if (i == 0) {
		gretl_matrix_I(A, nv);
		gretl_matrix_add_to(A, Pi);
		gretl_matrix_add_to(A, G[0]);
	    } else if (i == p) {
		gretl_matrix_zero(A);
		gretl_matrix_subtract_from(A, G[i-1]);
	    } else {
		gretl_matrix_copy_values(A, G[i]);
		gretl_matrix_subtract_from(A, G[i-1]);
	    }
#if JDEBUG
	    fprintf(stderr, "A matrix, lag %d\n\n", i+1);
	    gretl_matrix_print(A, NULL);
#endif
	    add_Ai_to_VAR_A(A, vecm, i);
	}
    }

#if JDEBUG
    gretl_matrix_print(vecm->A, "vecm->A");
#endif

 bailout:

    gretl_matrix_free(Pi);
    gretl_matrix_array_free(G, p);
    gretl_matrix_free(A);

    return err;
}

/* Cointegration test: calculate and print the "long-run matrix",
   \alpha \beta', (what Hamilton calls \Zeta_0) */

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
				  tmp, GRETL_MOD_NONE);

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
print_beta_and_alpha (JohansenInfo *jv, gretl_matrix *evals, int h,
		      const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    pprintf(prn, "\n%s", _("eigenvalue"));
    for (i=0; i<h; i++) {
	pprintf(prn, "%#12.5g ", evals->val[i]);
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

/*
   renormalize \beta such that its uppermost submatrix of
   size rank * rank is the identity matrix:

   \beta' = [ I | *free elements* ]
*/

static int phillips_normalize_beta (GRETL_VAR *vecm)
{
    gretl_matrix *c = NULL;
    gretl_matrix *beta_c = NULL;

    int r = jrank(vecm);
    int n = gretl_matrix_rows(vecm->jinfo->Beta);
    int i, j, err = 0;

    double x;

    c = gretl_matrix_alloc(r, r);
    beta_c = gretl_matrix_alloc(n, r);
    if (c == NULL || beta_c == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<r; i++) {
        for (j=0; j<r; j++) {
	    x = gretl_matrix_get(vecm->jinfo->Beta, i, j);
	    gretl_matrix_set(c, i, j, x);
	}
    }

    /* form \beta_c = \beta c^{-1} */
    gretl_invert_general_matrix(c);
    gretl_matrix_multiply(vecm->jinfo->Beta, c, beta_c);

    /* correct rounding error: set true zeros in \beta_c */
    for (i=0; i<n; i++) {
	for (j=0; j<r; j++) {
	    if (i >= r) {
		if (gretl_matrix_get(beta_c, i, j) == -0) {
		    gretl_matrix_set(beta_c, i, j, 0);
		}
	    } else if (i == j) {
		gretl_matrix_set(beta_c, i, j, 1.0);
	    } else {
		gretl_matrix_set(beta_c, i, j, 0.0);
	    }
	}
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->Beta, "original beta");
    gretl_matrix_print(beta_c, "beta_c = beta * c^{-1}");
#endif

    gretl_matrix_copy_values(vecm->jinfo->Beta, beta_c);

 bailout:
    
    gretl_matrix_free(c);
    gretl_matrix_free(beta_c);

    return err;
}

/* VECM: compute the variance of the estimator of \beta, after doing
   Phillips normalization */

static int beta_variance (GRETL_VAR *vecm)
{
    gretl_matrix *O = NULL;
    gretl_matrix *aOa = NULL;
    gretl_matrix *HSH = NULL;

    int r = jrank(vecm);
    int m = gretl_matrix_cols(vecm->jinfo->Alpha);
    int n = gretl_matrix_rows(vecm->jinfo->Beta);
    int i, j, k, err = 0;

    double x;

    O = gretl_matrix_copy(vecm->S);
    aOa = gretl_matrix_alloc(m, m);
    HSH = gretl_matrix_alloc(n - r, n - r);

    if (O == NULL || aOa == NULL || HSH == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* compute \alpha' \Omega^{-1} \alpha */

    err = gretl_invert_symmetric_matrix(O);
    if (err) {
	goto bailout;
    }

    gretl_matrix_qform(vecm->jinfo->Alpha, GRETL_MOD_TRANSPOSE, O,
		       aOa, GRETL_MOD_NONE);

#if JDEBUG
    gretl_matrix_print(vecm->S, "vecm->S");
    gretl_matrix_print(O, "O = inverse(vecm->S)");
    gretl_matrix_print(vecm->jinfo->Alpha, "alpha_c");
    gretl_matrix_print(aOa, "aOa = alpha_c' * O * alpha_c");
#endif

    /* compute H'SH (just keep the south-east corner) */
    for (i=r; i<n; i++) {
	for (j=r; j<n; j++) {
	    x = gretl_matrix_get(vecm->jinfo->Svv, i, j);
	    gretl_matrix_set(HSH, i - r, j - r, x);
	}
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->Svv, "full Svv");
    gretl_matrix_print(HSH, "HSH = subset(Svv)");
#endif

    vecm->jinfo->Bvar = gretl_matrix_kronecker_product_new(aOa, HSH);
    if (vecm->jinfo->Bvar == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = gretl_invert_symmetric_matrix(vecm->jinfo->Bvar);
    if (err) {
	goto bailout;
    }

    gretl_matrix_divide_by_scalar(vecm->jinfo->Bvar, vecm->T);

    vecm->jinfo->Bse = gretl_matrix_alloc(n - r, r);
    if (vecm->jinfo->Bse == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    k = 0;
    for (j=0; j<r; j++) {
	/* cointegrating vector j */
        for (i=0; i<n-r; i++) {
	    x = gretl_matrix_get(vecm->jinfo->Bvar, k, k);
	    gretl_matrix_set(vecm->jinfo->Bse, i, j, sqrt(x));
	    k++;
	}
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->Bvar, "var(beta)");
    gretl_matrix_print(vecm->jinfo->Bse, "se(beta)");
#endif

 bailout:

    gretl_matrix_free(O);
    gretl_matrix_free(aOa);
    gretl_matrix_free(HSH);

    return err;
}

static int johansen_ll_calc (GRETL_VAR *jvar, const gretl_matrix *eigvals)
{
    gretl_matrix *Suu;
    double ldet, T_2 = (double) jvar->T / 2.0;
    int n = jvar->neqns;
    int h, i, err = 0;

    h = (jrank(jvar) > 0)? jrank(jvar) : n;

    Suu = gretl_matrix_copy(jvar->jinfo->Suu);

    if (Suu == NULL) {
	err = E_ALLOC;
    } else {
	ldet = gretl_matrix_log_determinant(Suu, &err);
	jvar->ll = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
	for (i=0; i<h; i++) {
	    jvar->ll -= T_2 * log(1.0 - eigvals->val[i]); 
	}
	gretl_matrix_free(Suu);
    }

    return err;
}

static int vecm_ll_stats (GRETL_VAR *vecm)
{
    gretl_matrix *S;
    int T = vecm->T;
    int g = vecm->neqns;
    int k = g * (vecm->order + 1);

    S = gretl_matrix_copy(vecm->S);
    if (S == NULL) {
	return E_ALLOC;
    } 

    vecm->ldet = gretl_vcv_log_determinant(S);
    gretl_matrix_free(S);

    /* FIXME: is k right (in all cases)? */
    k += vecm->jinfo->nexo;

    vecm->AIC = (-2.0 * vecm->ll + 2.0 * k * g) / T;
    vecm->BIC = (-2.0 * vecm->ll + log(T) * k * g) / T;
    vecm->HQC = (-2.0 * vecm->ll + 2.0 * log(log(T)) * k * g) / T;

    return 0;
}

static int 
compute_coint_test (GRETL_VAR *jvar, const gretl_matrix *evals, PRN *prn)
{
    int T = jvar->T;
    int n = jvar->neqns;
    double cumeig = 0.0;
    double *lmax = NULL;
    double *trace = NULL;
    double pvals[2];
    int i;

    trace = malloc(n * sizeof *trace);
    lmax = malloc(n * sizeof *lmax);

    if (trace == NULL || lmax == NULL) {
	free(trace);
	free(lmax);
	return E_ALLOC;
    }

    for (i=n-1; i>=0; i--){
	lmax[i] = -T * log(1.0 - evals->val[i]); 
	cumeig += lmax[i];
	trace[i] = cumeig; 
    }

    pputc(prn, '\n');
    print_Johansen_test_case(jcode(jvar), prn);
    pprintf(prn, "\n%s %s %s %s   %s  %s\n", _("Rank"), _("Eigenvalue"), 
	    _("Trace test"), _("p-value"),
	    _("Lmax test"), _("p-value"));	

    for (i=0; i<n; i++) {
	gamma_par_asymp(trace[i], lmax[i], jcode(jvar), n - i, pvals);
	pprintf(prn, "%4d%#11.5g%#11.5g [%6.4f]%#11.5g [%6.4f]\n", 
		i, evals->val[i], trace[i], pvals[0], lmax[i], pvals[1]);
    }
    pputc(prn, '\n');

    free(lmax);
    free(trace);

    return 0;
}

static int johansen_get_eigenvalues (gretl_matrix *Suu,
				     const gretl_matrix *Suv,
				     const gretl_matrix *Svv,
				     gretl_matrix *M,
				     gretl_matrix **evals,
				     int rank)
{
    gretl_matrix *Tmp = NULL;
    int n = Svv->cols;
    int err;

    err = gretl_invert_symmetric_matrix(Suu);
    if (err) {
	return err;
    }

    Tmp = gretl_matrix_alloc(n, n);
    if (Tmp == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_qform(Suv, GRETL_MOD_TRANSPOSE, 
		       Suu, Tmp, GRETL_MOD_NONE);

    *evals = gretl_gensymm_eigenvals(Tmp, Svv, M, &err);

    if (!err) {
	err = gretl_symmetric_eigen_sort(*evals, M, rank);
    }

    gretl_matrix_free(Tmp);

    return err;
}

/* Public entry point for cointegration test */

int johansen_coint_test (GRETL_VAR *jvar, const DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn)
{
    gretl_matrix *Suu = NULL;
    gretl_matrix *evals = NULL;

    int n = jvar->neqns;
    int m = gretl_matrix_cols(jvar->jinfo->Svv);
    int err = 0;

    jvar->jinfo->Beta = gretl_matrix_alloc(m, m);
    Suu = gretl_matrix_copy(jvar->jinfo->Suu);

    if (jvar->jinfo->Beta == NULL || Suu == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = johansen_get_eigenvalues(Suu, jvar->jinfo->Suv, jvar->jinfo->Svv,
				       jvar->jinfo->Beta, &evals, 0);
    }

    if (err) {
	pputs(prn, _("Failed to find eigenvalues\n"));
    } else {
	johansen_ll_calc(jvar, evals);
	compute_coint_test(jvar, evals, prn);

	if (!(opt & OPT_Q)) {
	    err = compute_alpha(jvar->jinfo, n);
	    if (!err) {
		print_beta_and_alpha(jvar->jinfo, evals, n, pdinfo, prn);
		compute_long_run_matrix(jvar->jinfo, n, pdinfo, prn);
	    }
	}
    }

    gretl_matrix_free(Suu);
    gretl_matrix_free(evals);

    return err;
}

static int johansen_prep_restriction (GRETL_VAR *jvar, 
				      gretl_matrix *Suu,
				      const gretl_matrix *D)
{
    gretl_matrix *Svv = NULL;
    gretl_matrix *Suv = NULL;
    int n = jvar->neqns;
    int m = D->cols;
    int err = 0;

    Svv = gretl_matrix_alloc(m, m);
    Suv = gretl_matrix_alloc(n, m);
    jvar->jinfo->Beta = gretl_matrix_alloc(D->rows, jrank(jvar));
    
    if (Svv == NULL || Suv == NULL || jvar->jinfo->Beta == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	/* calculate Svv <- D' Svv D */
	err = gretl_matrix_qform(D, GRETL_MOD_TRANSPOSE,
				 jvar->jinfo->Svv, 
				 Svv, GRETL_MOD_NONE);
    }

    if (!err) {
	/* Suv <- SuvD */
	err = gretl_matrix_multiply(jvar->jinfo->Suv, D, Suv);
    } 

    if (!err) {
	gretl_matrix_free(jvar->jinfo->Svv);
	jvar->jinfo->Svv = Svv;
	gretl_matrix_free(jvar->jinfo->Suv);
	jvar->jinfo->Suv = Suv;
    } else {
	gretl_matrix_free(Svv);
	gretl_matrix_free(Suv);
    }

    return err;
}

/* Public entry point for VECM estimation (with "Case 1" restriction
   on beta, if D != NULL) 
*/

int johansen_estimate (GRETL_VAR *jvar, const gretl_matrix *D,
		       double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    gretl_matrix *M = NULL;
    gretl_matrix *Suu = NULL;
    gretl_matrix *evals = NULL;

    int rank = jrank(jvar);
    int m, err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting johansen_estimate()\n\n");
#endif

    if (D != NULL) {
	m = gretl_matrix_cols(D);
    } else {
	m = gretl_matrix_cols(jvar->jinfo->Svv);
    }

    M = gretl_matrix_alloc(m, m);
    Suu = gretl_matrix_copy(jvar->jinfo->Suu);

    if (M == NULL || Suu == NULL) {
	err = E_ALLOC;
    }

    if (!err && D != NULL) {
	err = johansen_prep_restriction(jvar, Suu, D);
    }

    if (!err) {
	err = johansen_get_eigenvalues(Suu, jvar->jinfo->Suv, jvar->jinfo->Svv,
				       M, &evals, rank);
    }

    if (err) {
	pputs(prn, _("Failed to find eigenvalues\n"));
	goto bailout;
    } 

#if JDEBUG
    gretl_matrix_print(M, "raw eigenvector(s)");
#endif

    if (D != NULL) {
	err = gretl_matrix_multiply(D, M, jvar->jinfo->Beta);
    } else {
	jvar->jinfo->Beta = M;
	M = NULL;
    }

    if (!err) {
	int do_stderrs = rank < jvar->neqns;

	err = johansen_ll_calc(jvar, evals);

	if (!err) {
	    err = phillips_normalize_beta(jvar); 
	}
	if (!err) {
	    err = build_VECM_models(jvar, pZ, pdinfo, 0);
	}
	if (!err) {
	    err = compute_omega(jvar);
	}
	if (!err && do_stderrs && D == NULL) {
	    err = beta_variance(jvar);
	}
	if (!err) {
	    err = gretl_VAR_do_error_decomp(jvar->S, jvar->C);
	}
	if (!err) {
	    err = vecm_ll_stats(jvar);
	}
    } 

 bailout:    

    gretl_matrix_free(M);
    gretl_matrix_free(Suu);
    gretl_matrix_free(evals);

    return err;
}

/* Simplified version of the Johansen procedure, to be called in
   the process of computing bootstrap confidence intervals for
   impulse response functions.  We just have to do enough to
   generate the VAR representation.
*/

int 
johansen_boots_round (GRETL_VAR *jvar, double ***pZ, DATAINFO *pdinfo,
		      int iter)
{
    gretl_matrix *M = NULL;
    gretl_matrix *evals = NULL;
    int m = gretl_matrix_cols(jvar->jinfo->Svv);
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting johansen_bootstrap_round()\n\n");
#endif

    M = gretl_matrix_alloc(m, m);
    if (M == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = johansen_get_eigenvalues(jvar->jinfo->Suu, jvar->jinfo->Suv, 
				   jvar->jinfo->Svv, M, &evals, 
				   jrank(jvar));
    if (err) {
	goto bailout;
    }

#if JDEBUG
    gretl_matrix_print(M, "raw eigenvector(s)");
#endif

    if (!err) {
	if (jvar->jinfo->Beta == NULL) {
	    jvar->jinfo->Beta = gretl_matrix_copy(M);
	} else {
	    gretl_matrix_copy_values(jvar->jinfo->Beta, M);
	}
	if (jvar->jinfo->Beta == NULL) {
	    err = E_ALLOC;
	}
	if (!err) {
	    err = phillips_normalize_beta(jvar); 
	}
	if (!err) {
	    err = build_VECM_models(jvar, pZ, pdinfo, iter);
	}
	if (!err) {
	    err = compute_omega(jvar);
	}
    } 

 bailout:    

    gretl_matrix_free(M);
    gretl_matrix_free(evals);

    return err;
}

static int
johansen_LR_calc (GRETL_VAR *jvar, const gretl_matrix *evals, 
		  const gretl_matrix *D, int save_ll, PRN *prn)
{
    gretl_matrix *Suu;
    double llr = 0.0;
    double ldet = 0.0;
    double T_2 = (double) jvar->T / 2.0;
    int n = jvar->neqns;
    int h, i, err = 0;

    h = (jrank(jvar) > 0)? jrank(jvar) : n;

    Suu = gretl_matrix_copy(jvar->jinfo->Suu);

    if (Suu == NULL) {
	err = E_ALLOC;
    } else {
	ldet = gretl_matrix_log_determinant(Suu, &err);
    }

    if (!err) {
	llr = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
	for (i=0; i<h; i++) {
	    pprintf(prn, _("eigenvalue %d = %g\n"), i+1, evals->val[i]);
	    llr -= T_2 * log(1.0 - evals->val[i]); 
	}
	pputc(prn, '\n');
    }

    if (Suu != NULL) {
	gretl_matrix_free(Suu);
    }

    if (!err) {
	double x = 2.0 * (jvar->ll - llr);
	int nb = gretl_matrix_rows(jvar->jinfo->Beta);
	int df = h * (nb - gretl_matrix_cols(D));

	pprintf(prn, _("Unrestricted loglikelihood (lu) = %g\n"), jvar->ll);
	pprintf(prn, _("Restricted loglikelihood (lr) = %g\n"), llr);
	pprintf(prn, "2 * (lu - lr) = %g\n", x);
	pprintf(prn, _("P(Chi-Square(%d) > %g = %g\n"), df, x, 
		chisq_cdf_comp(x, df));

	if (save_ll) {
	    jvar->ll = llr;
	}
    }

    return err;
}

/* 
   Test of (homogeneous) linear restrictions on the cointegrating
   relations in a VECM.  This all needs verification and possibly
   fixing in the case where the "unrestricted" VECM includes a
   restricted constant or trend
*/

int vecm_beta_test (GRETL_VAR *jvar, const gretl_matrix *D, const DATAINFO *pdinfo, 
		    PRN *prn)
{
    gretl_matrix *M = NULL;
    gretl_matrix *Svv = NULL;
    gretl_matrix *Suv = NULL;
    gretl_matrix *Suu = NULL;
    gretl_matrix *evals = NULL;

    int n = jvar->neqns;
    int m = gretl_matrix_cols(D);
    int rank = jrank(jvar);
    int err = 0;

    M = gretl_matrix_alloc(m, m);
    Svv = gretl_matrix_alloc(m, m);
    Suv = gretl_matrix_alloc(n, m);
    Suu = gretl_matrix_copy(jvar->jinfo->Suu);

    if (M == NULL || Svv == NULL || Suv == NULL || Suu == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    pputs(prn, "\nTest of restrictions on cointegrating relations\n\n");

    gretl_matrix_print_to_prn(D, "Restriction matrix, D", prn);

    /* calculate Svv <- D' Svv D */
    gretl_matrix_qform(D, GRETL_MOD_TRANSPOSE,
		       jvar->jinfo->Svv, Svv, GRETL_MOD_NONE);

    gretl_matrix_print_to_prn(Svv, "D'SvvD", prn);

    if (!err) {
	/* Suv <- SuvD */
	err = gretl_matrix_multiply(jvar->jinfo->Suv, D, Suv);
    }

    gretl_matrix_print_to_prn(Suv, "SuvD", prn);

    err = johansen_get_eigenvalues(Suu, Suv, Svv, M, &evals, rank);

    if (!err) {
	gretl_matrix_print_to_prn(M, "M", prn);
	johansen_LR_calc(jvar, evals, D, 0, prn);
    } 

 bailout:    

    gretl_matrix_free(M);
    gretl_matrix_free(Suu);
    gretl_matrix_free(evals);
    gretl_matrix_free(Svv);
    gretl_matrix_free(Suv);

    return err;
}
