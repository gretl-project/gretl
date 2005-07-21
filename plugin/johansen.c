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
		 int N, int T, double *pval)
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
      T: sample size;
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

    if (T > 0) {
	double m2 = 0.0, v2 = 0.0;

	tracem = s_mTrace_m_time[det];
	tracev = s_mTrace_v_time[det];

	x[0] = sqrt((double) N) / T;
	x[1] = (double) N / T;
	x[2] = x[1] * x[1];
	x[3] = (N == 1)? 1.0 / T : 0.0;
	x[4] = (N == 1)? 1.0 : 0.0;
	x[5] = (N == 2)? 1.0 : 0.0;
	x[6] = (N == 3)? 1.0 : 0.0;

	for (i=0; i<7; i++) {
	    m2 += x[i] * tracem[i];
	    v2 += x[i] * tracev[i];
	}

	mt *= exp(m2);
	vt *= exp(v2);
    }

    g = gamma_dist(mt, vt, tracetest, 2);
    if (na(g)) {
	pval[0] = NADBL;
    } else {
	pval[0] = 1.0 - g;
    }

    g = gamma_dist(ml, vl, lmaxtest, 2);
    if (na(g)) {
	pval[1] = NADBL;
    } else {
	pval[1] = 1.0 - g;
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
johansen_normalize (gretl_matrix *A, const gretl_matrix *Svv, int j)
{
    gretl_matrix *a = NULL, *b = NULL;
    double x, den;
    int i, err = 0;
    int k = Svv->rows;

    a = gretl_matrix_alloc(k, 1);
    b = gretl_matrix_alloc(k, 1);

    if (a == NULL || b == NULL) {
	gretl_matrix_free(a);
	gretl_matrix_free(b);
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	x = gretl_matrix_get(A, i, j);
	gretl_matrix_set(a, i, 0, x);
    }

    gretl_matrix_multiply(Svv, a, b);

    den = gretl_matrix_dot_product(a, GRETL_MOD_TRANSPOSE,
				   b, GRETL_MOD_NONE, &err);

    if (!err) {
	den = sqrt(den);
	for (i=0; i<k; i++) {
	    x = gretl_matrix_get(A, i, j);
	    gretl_matrix_set(A, i, j, x / den);
	}
    } 
    
    gretl_matrix_free(a);
    gretl_matrix_free(b);

    return err;
}

static void 
print_coint_vecs (struct eigval *evals, const gretl_matrix *vr, 
		  int k, PRN *prn)
{
    int i, j, rows = vr->rows;

    for (i=0; i<k; i++) {
	double ev, x = 0.0;
	int col = evals[i].idx;

	if (k > 1) {
	    pprintf(prn, "(%d) %s = %g\n", i + 1, _("Eigenvalue"), evals[i].v);
	} else {
	    pprintf(prn, "%s = %g\n", _("Eigenvalue"), evals[i].v);
	}

	pprintf(prn, " %s: [ ", _("cointegrating vector"));
	for (j=0; j<rows; j++) {
	    pprintf(prn, "%10.5g ", gretl_matrix_get(vr, j, col));
	}
	pputs(prn, "]\n");

	pprintf(prn, " %s:         [ ", _("renormalized"));
	for (j=0; j<rows; j++) {
	    if (j == 0) {
		x = gretl_matrix_get(vr, j, col);
		pprintf(prn, "%10.5g ", 1.0);
	    } else {
		ev = gretl_matrix_get(vr, j, col) / x;
		pprintf(prn, "%10.5g ", ev);
	    }
	}
	pputs(prn, "]\n");

    }
    pputc(prn, '\n');
}

static void print_test_case (JohansenCode jcode, PRN *prn)
{
    pputc(prn, '\n');

    switch (jcode) {
    case J_NO_CONST:
	pputs(prn, "Case 1: no constant");
	break;
    case J_REST_CONST:
	pputs(prn, "Case 2: restricted constant");
	break;
    case J_UNREST_CONST:
	pputs(prn, "Case 3: unrestricted constant");
	break;
    case J_REST_TREND:
	pputs(prn, "Case 4: restricted trend");
	break;
    case J_UNREST_TREND:
	pputs(prn, "Case 5: unrestricted trend");
	break;
    default:
	break;
    }

    pputc(prn, '\n');
}

char *safe_print_pval (double p, int i)
{
    static char pv0[32], pv1[32];

    if (!na(p)) {
	sprintf((i == 0)? pv0 : pv1, "%6.4f", p);
    } else {
	strcpy((i == 0)? pv0 : pv1, "  NA  ");
    }

    return (i== 0)? pv0 : pv1;
}

int johansen_eigenvals (gretl_matrix *Suu, gretl_matrix *Svv, gretl_matrix *Suv,
			int T, JohansenCode jcode, PRN *prn)
{
    int k = gretl_matrix_cols(Suu);
    int kv = gretl_matrix_cols(Svv);
    gretl_matrix *M = NULL, *Svvcpy = NULL;
    gretl_matrix *TmpL = NULL, *TmpR = NULL;
    double *eigvals = NULL;
    int r = 0, err = 0;

    TmpL = gretl_matrix_alloc(kv, k);
    TmpR = gretl_matrix_alloc(kv, kv);
    M = gretl_matrix_alloc(kv, kv);
    Svvcpy = gretl_matrix_copy(Svv);

    if (Suu == NULL || Svv == NULL || Suv == NULL ||
	TmpL == NULL || TmpR == NULL ||
	M == NULL || Svvcpy == NULL) {
	err = 1;
	goto eigenvals_bailout;
    }

    if (kv > k) {
	gretl_matrix_reuse(TmpR, k, kv);
    }

    /* calculate Suu^{-1} Suv */
    err = gretl_invert_general_matrix(Suu);
    if (!err) {
	err = gretl_matrix_multiply(Suu, Suv, TmpR);
    }

    /* calculate Svv^{-1} Suv' */
    if (!err) {
	err = gretl_invert_general_matrix(Svv);
    }
    if (!err) {
	err = gretl_matrix_multiply_mod(Svv, GRETL_MOD_NONE,
					Suv, GRETL_MOD_TRANSPOSE, 
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
	gretl_matrix_zero(TmpR);
    }

    eigvals = gretl_general_matrix_eigenvals(M, TmpR);

    if (eigvals != NULL) {
	double cumeig = 0.0;
	double *lambdamax = NULL, *trace = NULL;
	struct eigval *evals;
	double pval[2];
	int i;

	/* in case kv > k */
	k = kv;

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
	    evals[i].v = eigvals[i];
	    evals[i].idx = i;
	}

	qsort(evals, k, sizeof *evals, inverse_compare_doubles);

	for (i=k-1; i>=0; i--){
      	    lambdamax[i] = -T * log(1.0 - evals[i].v); 
	    cumeig += lambdamax[i];
 	    trace[i] = cumeig; 
	}

	print_test_case(jcode, prn);

	/* first col shows cointegration rank under H0, 
	   second shows associated eigenvalue */
	pprintf(prn, "\n%s %s %s %s   %s %s\n", _("Rank"), _("Eigenvalue"), 
		_("Trace test"), _("p-value"),
		_("Lmax test"), _("p-value"));

	for (i=0; i<k; i++) {
	    /* second-last arg below was T, but Doornik does not use sample size */
	    gamma_par_asymp(trace[i], lambdamax[i], jcode, k-i, 0, pval);
	    pprintf(prn, "%4d%11.4f%11.4f [%s]%11.4f [%s]\n", \
		    i, evals[i].v, trace[i], safe_print_pval(pval[0], 0), 
		    lambdamax[i], safe_print_pval(pval[1], 1));
	    if (pval[0] < 0.05) {
		r++;
	    }
	}
	pputc(prn, '\n');

	/* tmpR holds the eigenvectors */
	johansen_normalize(TmpR, Svvcpy, 0);

	if (r > 0) {
	    pputs(prn, 
		    /* xgettext:no-c-format */
		    _("Cointegrating vectors (trace test, 5% significance level):")); 
	    pputc(prn, '\n');
	    print_coint_vecs(evals, TmpR, r, prn);
	} else {
	    pputs(prn, 
		  /* xgettext:no-c-format */
		  _("No cointegrating vectors (trace test, 5% significance level)"));
	    pputc(prn, '\n');
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
    gretl_matrix_free(Svvcpy);

    return err;
}
