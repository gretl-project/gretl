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
gamma_par_asymp (double tracetest, double lmaxtest, int det, 
		 int N, int T, double *pval)
{
    /*
      Asymptotic critical values for Johansen's LR tests via gamma approximation

      params:
      tracetest, lmaxtest: trace and lambdamax est. statistics
      det: deterministic trends; 
      0 = no constant
      1 = restricted constant
      2 = unrestricted constant
      4 = restricted trend
      5 = unrestricted trend
      N: cointegration rank under H0;
      T: sample size;
      pval: on output, array of pvalues, for the two tests;
    */
    
    double mt, vt, ml, vl, *x;
    const double *tracem, *tracev, *lmaxm, *lmaxv;
    int i;

    tracem = s_mTrace_m_coef[det];
    tracev = s_mTrace_v_coef[det];
    lmaxm = s_mMaxev_m_coef[det];
    lmaxv = s_mMaxev_v_coef[det];

    mt = vt = 0.0;
    ml = vl = 0.0;

    x = malloc(7 * sizeof *x);
    if (x == NULL) return 1;

    x[0] = N * N;
    x[1] = N;
    x[2] = 1.0;
    x[3] = (N == 1) ? 1.0 : 0.0 ;
    x[4] = (N == 2) ? 1.0 : 0.0 ;
    x[5] = sqrt((double) N);

    for (i=0; i<6; i++) {
	mt += x[i] * tracem[i];
	vt += x[i] * tracev[i];
	if(i){
	    ml += x[i] * lmaxm[i-1];
	    vl += x[i] * lmaxv[i-1];
	}
    }

    if (T > 0) {
	double m2 = 0.0, v2 = 0.0;

	tracem = s_mTrace_m_time[det];
	tracev = s_mTrace_v_time[det];

	x[0] = sqrt((double) N) / T;
	x[1] = N / T;
	x[2] = x[1] * x[1];
	x[3] = (N == 1) ? 1.0/ T : 0.0;
	x[4] = (N == 1) ? 1.0 : 0.0;
	x[5] = (N == 2) ? 1.0 : 0.0;
	x[6] = (N == 3) ? 1.0 : 0.0;

	for (i=0; i<7; i++) {
	    m2 += x[i] * tracem[i];
	    v2 += x[i] * tracev[i];
	}

	mt *= exp(m2);
	vt *= exp(v2);
    }

    free(x);

    pval[0] = 1.0 - gamma_dist(mt, vt, tracetest, 2);
    pval[1] = 1.0 - gamma_dist(ml, vl, lmaxtest, 2);

    return 0;
}

static int inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da < *db) - (*da > *db);
}

int johansen_eigenvals (const double **X, const double **Y, const double **Z, 
			int k, int T, int trends, PRN *prn)
{
    gretl_matrix *Suu, *Svv, *Suv;
    gretl_matrix *Inv, *TmpL, *TmpR, *M;
    double *eigvals;
    int err = 0;

    Suu = gretl_matrix_from_2d_array(X, k, k);
    Svv = gretl_matrix_from_2d_array(Y, k, k);
    Suv = gretl_matrix_from_2d_array(Z, k, k);

    Inv = gretl_matrix_alloc(k, k);
    TmpL = gretl_matrix_alloc(k, k);
    TmpR = gretl_matrix_alloc(k, k);
    M = gretl_matrix_alloc(k, k);

    /* calculate Suu^{-1} Suv */
    gretl_invert_general_matrix(Suu);
    gretl_matrix_multiply(Suu, Suv, TmpR);

    /* calculate Svv^{-1} Suv' */
    gretl_invert_general_matrix(Svv);
    gretl_matrix_multiply_mod(Svv, GRETL_MOD_NONE,
			      Suv, GRETL_MOD_TRANSPOSE, 
			      TmpL);

    gretl_matrix_multiply(TmpL, TmpR, M);

    eigvals = gretl_general_matrix_eigenvals(M);

    if (eigvals != NULL) {
	int i;
	double cumeig = 0.0;
	double *lambdamax = NULL, *trace = NULL;
	double pval[2];

	trace = malloc(k * sizeof *trace);
	lambdamax = malloc(k * sizeof *lambdamax);
	if (trace == NULL || lambdamax == NULL) {
	    free(trace);
	    free(lambdamax);
	    err = 1;
	    goto eigenvals_bailout;
	}

	qsort(eigvals, k, sizeof *eigvals, inverse_compare_doubles);

	for (i=k-1; i>=0; i--){
      	    lambdamax[i] = -T * log(1.0 - eigvals[i]); 
	    cumeig += lambdamax[i];
 	    trace[i] = cumeig; 
	}

	pputs(prn, _("\nRank Eigenvalue Trace test [p.val.]  Lmax test [p.val.]\n"));

	for (i=0; i<k; i++) {
	    gamma_par_asymp(trace[i], lambdamax[i], 2 , k-i, T, pval);
	    pprintf(prn, "%4d%11.4f%11.4f [%6.4f]%11.4f [%6.4f]\n", \
		    i, eigvals[i], trace[i], pval[0], lambdamax[i], pval[1]);
	}
	pputc(prn, '\n');

	free(eigvals);
	free(lambdamax);
	free(trace);

    } else {
	pputs(prn, _("Failed to find eigenvalues\n"));
    }

 eigenvals_bailout:    

    gretl_matrix_free(Svv);
    gretl_matrix_free(Suu);
    gretl_matrix_free(Suv);

    gretl_matrix_free(Inv);
    gretl_matrix_free(TmpL);
    gretl_matrix_free(TmpR);
    gretl_matrix_free(M);

    return err;
}
