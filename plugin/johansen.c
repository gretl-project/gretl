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
#include "libset.h"
#include "jprivate.h"

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
    double x[7];
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
	if (i > 0) {
	    ml += x[i] * lmaxm[i-1];
	    vl += x[i] * lmaxv[i-1];
	}
    }

    pval[0] = gamma_cdf_comp(mt, vt, tracetest, 2);
    pval[1] = gamma_cdf_comp(ml, vl, lmaxtest, 2);

    return 0;
}

/* Remove a possible excess zero from the end of a floating point
   number printed to the given precision p (working around a bug 
   in the C library).
*/

static void fix_xstr (char *s, int p)
{
    int n = strlen(s);

    if (n > p && strspn(s + n - p, "0") == p) {
	s[n-1] = 0;
    }
}

#define ABMIN 1.0e-15

/* for cointegration test: print cointegrating vectors or adjustments,
   either "raw" or re-scaled */

static void print_beta_or_alpha (GRETL_VAR *jvar, int k,
				 const DATAINFO *pdinfo, PRN *prn,
				 int job, int rescale)
{
    JohansenInfo *jv = jvar->jinfo;
    gretl_matrix *c = (job == V_BETA)? jv->Beta : jv->Alpha;
    int rows = gretl_matrix_rows(c);
    char xstr[32];
    int i, j;
    double x, y;

    if (rescale) {
	pprintf(prn, "\n%s\n", (job == V_BETA)? 
		_("renormalized beta") :
		_("renormalized alpha"));
    } else {
	pprintf(prn, "\n%s\n", (job == V_BETA)? 
		_("beta (cointegrating vectors)") : 
		_("alpha (adjustment vectors)"));
    }

    for (i=0; i<rows; i++) {
	if (i < jvar->ylist[0]) {
	    pprintf(prn, "%-10s", pdinfo->varname[jvar->ylist[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    pprintf(prn, "%-10s", "const");
	} else if (jv->code == J_REST_TREND) {
	    pprintf(prn, "%-10s", "trend");
	}
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(c, i, j);
	    if (rescale) {
		y = gretl_matrix_get(jv->Beta, j, j);
		if (job == V_BETA) {
		    x /= y;
		} else {
		    x *= y;
		}
	    }
	    if (x == -0.0 || fabs(x) < ABMIN) {
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

static int compute_alpha (JohansenInfo *jv)
{
    const gretl_matrix *B = jv->Beta;
    gretl_matrix *alpha = NULL;
    gretl_matrix *BSB = NULL;
    gretl_matrix *Tmp = NULL;
    int err = 0;

    BSB = gretl_matrix_alloc(B->cols, B->cols);
    Tmp = gretl_matrix_alloc(B->rows, B->cols);
    alpha = gretl_matrix_alloc(jv->S01->rows, B->cols);

    if (BSB == NULL || Tmp == NULL || alpha == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	err = gretl_matrix_qform(B, GRETL_MOD_TRANSPOSE, jv->S11,
				 BSB, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(BSB);
    }

    if (!err) {
	gretl_matrix_multiply(B, BSB, Tmp);
	gretl_matrix_multiply(jv->S01, Tmp, alpha);
    }

    gretl_matrix_free(BSB);
    gretl_matrix_free(Tmp);

    if (!err) {
	jv->Alpha = alpha;
    } else {
	gretl_matrix_free(alpha);
    }

    return err;
}

/* print the long-run matrix, \alpha \beta' */

static int print_long_run_matrix (GRETL_VAR *jvar, 
				  const DATAINFO *pdinfo, 
				  PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    gretl_matrix *Pi;
    double x;
    int i, j;

    Pi = gretl_matrix_alloc(jv->Alpha->rows, jv->Beta->rows);
    if (Pi == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(jv->Alpha, GRETL_MOD_NONE,
			      jv->Beta, GRETL_MOD_TRANSPOSE,
			      Pi, GRETL_MOD_NONE);

    pprintf(prn, "%s\n", _("long-run matrix (alpha * beta')"));

    pprintf(prn, "%22s", pdinfo->varname[jvar->ylist[1]]); /* N.B. */
    for (j=2; j<=jvar->ylist[0]; j++) {
	pprintf(prn, "%13s", pdinfo->varname[jvar->ylist[j]]);
    }

    if (jv->code == J_REST_CONST) {
	pprintf(prn, "%13s", "const");
    } else if (jv->code == J_REST_TREND) {
	pprintf(prn, "%13s", "trend");
    }    

    pputc(prn, '\n');

    for (i=0; i<Pi->rows; i++) {
	pprintf(prn, "%-10s", pdinfo->varname[jvar->ylist[i+1]]);
	for (j=0; j<Pi->cols; j++) {
	    x = gretl_matrix_get(Pi, i, j);
	    if (fabs(x) < 0.5e-14) {
		x = 0.0;
	    }
	    pprintf(prn, "%#12.5g ", x);
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    gretl_matrix_free(Pi);

    return 0;
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

/* After doing OLS estimation of the VECM conditional on \beta: copy
   the coefficients on the lagged differences (i.e. form the \Gamma
   matrices) so we can compute the VAR representation */

static void copy_coeffs_to_Gamma (GRETL_VAR *vecm, gretl_matrix **G)
{
    int i, j, k, h;
    double x;

    for (i=0; i<vecm->neqns; i++) {
	for (k=0; k<vecm->order; k++) {
	    h = k + vecm->ifc;
	    /* successive lags (distinct \Gamma_i matrices) */
	    for (j=0; j<vecm->neqns; j++) {
		/* successive \Delta x_j */
		x = gretl_matrix_get(vecm->B, h, i);
		gretl_matrix_set(G[k], i, j, x);
		h += vecm->order;
	    }
	}
    }
}

/* \Pi, as will be used in forming the VAR representation */

static void form_Pi (GRETL_VAR *v, gretl_matrix *Pi)
{
    gretl_matrix_multiply_mod(v->jinfo->Alpha, GRETL_MOD_NONE,
			      v->jinfo->Beta, GRETL_MOD_TRANSPOSE,
			      Pi, GRETL_MOD_NONE);
}

/* After doing OLS estimation of the VECM conditional on \beta:
   copy the coefficients on the EC terms (\beta' X) into the \alpha
   matrix.
*/

static int copy_to_alpha (GRETL_VAR *v)
{
    int rank = v->jinfo->rank;
    int pos = v->ncoeff - rank;
    double x;
    int i, j;

    if (v->jinfo->Alpha == NULL) {
	v->jinfo->Alpha = gretl_matrix_alloc(v->neqns, rank);
	if (v->jinfo->Alpha == NULL) {
	    return E_ALLOC;
	}
    }

    for (i=0; i<v->neqns; i++) {
	for (j=0; j<rank; j++) {
	    x = gretl_matrix_get(v->B, pos + j, i);
	    gretl_matrix_set(v->jinfo->Alpha, i, j, x);
	}
    }

    return 0;
}

/* VAR representation: transcribe the coefficient matrix A_i (for lag
   i) into its place in the full VAR coefficient matrix, A 
*/

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

enum {
    ESTIMATE_ALPHA = 0,
    NET_OUT_ALPHA
};

/* For the case where we want to show complete output: transcribe
   the matrix multi-OLS results into MODEL structs for printing
   and saving.
*/

static int 
transcribe_VECM_models (GRETL_VAR *vecm, 
			const double **Z, const DATAINFO *pdinfo,
			const gretl_matrix *XTX, int flag)
{
    MODEL *pmod;
    char **params = NULL;
    const char *dvname;
    const double *y;
    double x;
    int r = jrank(vecm);
    int yno, N = pdinfo->n;
    int i, j, jmax;
    int err = 0;

    jmax = (flag == ESTIMATE_ALPHA)? vecm->ncoeff :
	vecm->ncoeff - r;

    params = strings_array_new_with_length(vecm->ncoeff, VNAMELEN);
    if (params == NULL) {
	return E_ALLOC;
    }

    gretl_VAR_param_names(vecm, params, pdinfo);

    for (i=0; i<vecm->neqns; i++) {
	yno = vecm->ylist[i+1];
	y = Z[yno];
	dvname = pdinfo->varname[yno];

	pmod = vecm->models[i];
	pmod->ID = i + 1;
	pmod->ci = OLS;
	pmod->aux = AUX_VECM;

	pmod->full_n = N;
	pmod->nobs = vecm->T;
	pmod->t1 = vecm->t1;
	pmod->t2 = vecm->t2;
	pmod->ncoeff = vecm->ncoeff;
	pmod->dfd = pmod->nobs - pmod->ncoeff;
	pmod->ifc = vecm->ifc;
	pmod->dfn = pmod->ncoeff - pmod->ifc;

	err = gretl_model_allocate_storage(pmod);

	pmod->depvar = malloc(VNAMELEN);
	if (pmod->depvar != NULL) {
	    strcpy(pmod->depvar, "d_");
	    strncat(pmod->depvar, dvname, VNAMELEN - 3);
	}

	if (i == 0) {
	    pmod->params = params;
	} else {
	    pmod->params = strings_array_dup(params, vecm->ncoeff);
	}
	pmod->nparams = vecm->ncoeff;

	pmod->list = gretl_list_new(1);
	pmod->list[1] = yno;

	pmod->dfd = vecm->T;
	set_VAR_model_stats(vecm, i);
	pmod->dfd = vecm->T - pmod->ncoeff;

	for (j=0; j<jmax; j++) {
	    pmod->coeff[j] = gretl_matrix_get(vecm->B, j, i);
	    x = gretl_matrix_get(XTX, j, j);
	    pmod->sderr[j] = pmod->sigma * sqrt(x);
	}

	for (j=jmax; j<pmod->ncoeff; j++) {
	    pmod->coeff[j] = gretl_matrix_get(vecm->jinfo->Alpha, i, j-jmax);
	    pmod->sderr[j] = NADBL;
	}
    }

    return err;
}

/* The X (data) and B (coefficient) matrices may need expanding
   to take account of the EC terms */

static int vecm_check_size (GRETL_VAR *v, int flag)
{
    int err = 0;

    if (flag == ESTIMATE_ALPHA) {
	fprintf(stderr, "vecm_check_size: ncoeff: %d -> %d\n",
		v->ncoeff, v->ncoeff + jrank(v));
	v->ncoeff += jrank(v);
    } else if (v->ncoeff == 0) {
	return 0;
    }

    if (v->X == NULL) {
	v->X = gretl_matrix_alloc(v->T, v->ncoeff);
	if (v->X == NULL) {
	    err = E_ALLOC;
	}
    } else if (v->X->cols < v->ncoeff) {
	err = gretl_matrix_realloc(v->X, v->T, v->ncoeff);
    }

    if (err) {
	return err;
    }

    if (v->B == NULL) {
	v->B = gretl_matrix_alloc(v->ncoeff, v->neqns);
	if (v->B == NULL) {
	    err = E_ALLOC;
	}	
    } else if (v->B->rows < v->ncoeff) {
	err = gretl_matrix_realloc(v->B, v->ncoeff, v->neqns);
    }

    if (err) {
	return err;
    }

    if (v->Y->cols > v->neqns) {
	gretl_matrix_reuse(v->Y, -1, v->neqns);
    }

    return 0;
}

/* For estimating both alpha and Gamma: add the EC terms into
   the X data matrix */

static int add_EC_terms_to_X (GRETL_VAR *v, const double **Z)
{
    const gretl_matrix *B = v->jinfo->Beta;
    int rank = jrank(v);
    int col = v->ncoeff - rank;
    double xt, bxt, bij;
    int i, j, s, t;
    int err = 0;

    for (j=0; j<rank; j++) {
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    bxt = 0.0;

	    /* beta * X(t-1) */
	    for (i=0; i<v->neqns; i++) {
		xt = Z[v->ylist[i+1]][t-1];
		bij = gretl_matrix_get(B, i, j);
		bxt += bij * xt;
	    }

	    /* restricted const or trend */
	    if (restricted(v)) {
		bij = gretl_matrix_get(B, i, j);
		if (jcode(v) == J_REST_TREND) {
		    bij *= t;
		}
		bxt += bij;
	    }

	    gretl_matrix_set(v->X, s++, col, bxt);
	}
	col++;
    }
	
    return err;
}

/* preparing for OLS conditional on beta: construct the
   appropriate dependent variable matrix */

static int make_vecm_models_Y (GRETL_VAR *v, const double **Z, 
			       gretl_matrix *Pi, int flag)
{
    int i, s, t, vi, vj;
    double yti, xti;
    int err = 0;

    if (flag == ESTIMATE_ALPHA) {
	/* "Y" is composed of plain DYt */
	for (i=0; i<v->neqns; i++) {
	    vi = v->ylist[i+1];
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		yti = Z[vi][t] - Z[vi][t-1];
		gretl_matrix_set(v->Y, s++, i, yti);
	    }
	}
    } else {
	int j, p1 = v->jinfo->Beta->rows;

	err = compute_alpha(v->jinfo);
	if (err) {
	    return err;
	}

	form_Pi(v, Pi);

	/* form "Y" = DY_t - \Pi Y*_t */
	for (i=0; i<v->neqns; i++) {
	    vi = v->ylist[i+1];
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		/* first difference */
		yti = Z[vi][t] - Z[vi][t-1];
		for (j=0; j<p1; j++) {
		    if (j < v->neqns) {
			/* lagged level */
			vj = v->ylist[j+1];
			xti = Z[vj][t-1];
		    } else {
			xti = (jcode(v) == J_REST_TREND)? t : 1;
		    } 
		    yti -= gretl_matrix_get(Pi, i, j) * xti;
		}
		gretl_matrix_set(v->Y, s++, i, yti);
	    }
	}
    } 

    return err;
}

/* The following is designed to accommodate the case where alpha is
   restricted, so we can't just run OLS conditional on beta.

   DY_t = \Pi Y*_{t-1} + \sum_{i=1}^{k-1}\Gamma_i DY_{t-1} + ...

   Subtract \Pi Y*_t from both sides, call DY_t - \Pi Y*_t "Yt",
   and call the lagged DYs "X".  Regress "Y" on "X" to find
   estimates of the \Gammas, etc.  

   But the function also handles the case where \alpha will
   be estimated along with \Gamma.
*/

static int 
build_VECM_models (GRETL_VAR *v, const double **Z, const DATAINFO *pdinfo,
		   int flag, int iter)
{
    gretl_matrix *Pi = NULL;
    gretl_matrix *XTX = NULL;
    gretl_matrix *Ai = NULL;
    gretl_matrix **G = NULL;

    gretl_matrix *beta = v->jinfo->Beta;

    int order = v->order;
    int nc, n = v->neqns;
    int i, err;

    err = vecm_check_size(v, flag);
    if (err) {
	return err;
    }

    nc = v->ncoeff;

    Pi = gretl_matrix_alloc(n, beta->rows);
    if (Pi == NULL) {
	err = E_ALLOC;
    }

    if (!err && nc > 0) {
	Ai = gretl_matrix_alloc(n, n);
	if (Ai == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && order > 0) {
	G = gretl_matrix_array_alloc_with_size(order, n, n);
	if (G == NULL) {
	    err = E_ALLOC;
	}	
    }  

    /* FIXME bootstrapping !! */

    if (!err) {
	err = make_vecm_models_Y(v, Z, Pi, flag);
    }

    if (!err && flag == ESTIMATE_ALPHA) {
	err = add_EC_terms_to_X(v, Z);
    }

    if (!err) {
	if (nc > 0) {
	    /* run the regressions */
	    err = gretl_matrix_multi_ols(v->Y, v->X, v->B, v->E, &XTX);
	} else {
	    /* nothing to estimate, with alpha in hand */
	    gretl_matrix_copy_values(v->E, v->Y);
	    v->ncoeff = jrank(v);
	}
    }

    if (!err && order > 0) {
	copy_coeffs_to_Gamma(v, G);
    }

    if (!err && flag == ESTIMATE_ALPHA) {
	err = copy_to_alpha(v);
	if (!err) {
	    form_Pi(v, Pi);
	}
    }

    if (err) {
	goto bailout;
    }

    if (Pi->cols > n) {
	gretl_matrix_reuse(Pi, -1, n);
    }

    if (order == 0) {
	gretl_matrix_I(Ai, n);
	gretl_matrix_add_to(Ai, Pi);
	add_Ai_to_VAR_A(Ai, v, 0);
    } else {
	for (i=0; i<=order; i++) {
	    if (i == 0) {
		gretl_matrix_I(Ai, n);
		gretl_matrix_add_to(Ai, Pi);
		gretl_matrix_add_to(Ai, G[0]);
	    } else if (i == order) {
		gretl_matrix_zero(Ai);
		gretl_matrix_subtract_from(Ai, G[i-1]);
	    } else {
		gretl_matrix_copy_values(Ai, G[i]);
		gretl_matrix_subtract_from(Ai, G[i-1]);
	    }
#if JDEBUG
	    fprintf(stderr, "Ai matrix, lag %d\n\n", i+1);
	    gretl_matrix_print(Ai, NULL);
#endif
	    add_Ai_to_VAR_A(Ai, v, i);
	}
    }

#if JDEBUG
    gretl_matrix_print(v->A, "vecm->A");
#endif

    if (!err) {
	/* FIXME conditionality here? */
	transcribe_VECM_models(v, Z, pdinfo, XTX, flag);
    }

 bailout:
    
    gretl_matrix_free(Pi);
    gretl_matrix_free(XTX);
    gretl_matrix_free(Ai);
    gretl_matrix_array_free(G, order);

    return err;
}

/* Print both "raw" and re-scaled versions of the beta and alpha
   matrices (cointegrating vectors and vectors of adjustments
   respectively).
*/

static int
print_beta_and_alpha (GRETL_VAR *jvar, gretl_matrix *evals, int h,
		      const DATAINFO *pdinfo, PRN *prn)
{
    int i, err = 0;

    pprintf(prn, "\n%s", _("eigenvalue"));
    for (i=0; i<h; i++) {
	pprintf(prn, "%#12.5g ", evals->val[i]);
    }
    pputc(prn, '\n');

    /* "raw" vectors */
    print_beta_or_alpha(jvar, h, pdinfo, prn, V_BETA, 0);
    print_beta_or_alpha(jvar, h, pdinfo, prn, V_ALPHA, 0);

    /* re-scaled versions */
    print_beta_or_alpha(jvar, h, pdinfo, prn, V_BETA, 1);
    print_beta_or_alpha(jvar, h, pdinfo, prn, V_ALPHA, 1);

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
    err = gretl_invert_general_matrix(c);
    if (err) {
	fprintf(stderr, "phillips_normalize_beta: c is singular\n");
	goto bailout;
    }

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

static int col_normalize_beta (GRETL_VAR *vecm)
{
    gretl_matrix *B = vecm->jinfo->Beta;
    double x, den;
    int i, j;

    for (j=0; j<B->cols; j++) {
	den = gretl_matrix_get(B, j, j);
	if (den != 0.0) {
	    for (i=0; i<B->rows; i++) {
		x = gretl_matrix_get(B, i, j);
		gretl_matrix_set(B, i, j, x / den);
	    }
	}
    }

    return 0;
}

static int normalize_beta (GRETL_VAR *vecm, const gretl_matrix *R,
			   int *do_stderrs)
{
    if (R == NULL) {
	if (get_vecm_norm() == NORM_DIAG) {
	    if (do_stderrs != NULL) {
		*do_stderrs = 0;
	    }
	    return col_normalize_beta(vecm);
	} else {
	    return phillips_normalize_beta(vecm);
	}
    } else {
	gretl_matrix *B = vecm->jinfo->Beta;

	if (B->cols == 1) {
	    double den = B->val[0];
	    int i;

	    if (den != 0.0) {
		for (i=0; i<B->rows; i++) {
		    if (B->val[i] != 0) {
			B->val[i] /= den;
		    }
		}
	    }
	}
    } 

    return 0;
}

static gretl_matrix *make_H (const gretl_matrix *R, int *err)
{
    gretl_matrix *R1;
    gretl_matrix *H;
    double x;
    int n = R->cols;
    int i, j;

    /* FIXME how general is this? */

    R1 = gretl_zero_matrix_new(R->rows + 1, n);
    if (R1 == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    R1->val[0] = 1;

    for (i=1; i<R1->rows; i++) {
	for (j=0; j<n; j++) {
	    x = gretl_matrix_get(R, i-1, j);
	    gretl_matrix_set(R1, i, j, x);
	}
    }

    H = gretl_matrix_right_nullspace(R1, err);

    gretl_matrix_free(R1);
    
    return H;
}

/* VECM: compute the variance of the estimator of \beta, either
   under restriction R or after doing Phillips normalization.
   FIXME this is not right yet when R is given. 
*/

static int beta_variance (GRETL_VAR *vecm, const gretl_matrix *R)
{
    gretl_matrix *O = NULL;
    gretl_matrix *aOa = NULL;
    gretl_matrix *H = NULL;
    gretl_matrix *HSH = NULL;
    double x;
    int r = jrank(vecm);
    int m = gretl_matrix_cols(vecm->jinfo->Alpha);
    int n = gretl_matrix_rows(vecm->jinfo->Beta);
    int nh = n - r;
    int i, j, err = 0;

    if (R != NULL) {
	if (r > 1) {
	    /* not handled */
	    return 0;
	}
	H = make_H(R, &err);
	if (err) {
	    return err;
	}
	nh = H->cols;
    }    

    O = gretl_matrix_copy(vecm->S);
    aOa = gretl_matrix_alloc(m, m);
    HSH = gretl_matrix_alloc(nh, nh);

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

    /* compute H'*S11*H */

    if (H != NULL) {
	gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE,
			   vecm->jinfo->S11, 
			   HSH, GRETL_MOD_NONE);
    } else {
	/* phillips: just keep the south-east corner */
	for (i=r; i<n; i++) {
	    for (j=r; j<n; j++) {
		x = gretl_matrix_get(vecm->jinfo->S11, i, j);
		gretl_matrix_set(HSH, i - r, j - r, x);
	    }
	}
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->S11, "full S11");
    gretl_matrix_print(H, "H");
    gretl_matrix_print(HSH, "H'*S11*H");
#endif

    vecm->jinfo->Bvar = gretl_matrix_kronecker_product_new(aOa, HSH, &err);
    if (err) {
	goto bailout;
    }

    err = gretl_invert_symmetric_matrix(vecm->jinfo->Bvar);
    if (err) {
	goto bailout;
    }

    vecm->jinfo->Bse = gretl_zero_matrix_new(n, r);
    if (vecm->jinfo->Bse == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_divide_by_scalar(vecm->jinfo->Bvar, vecm->T);

    if (H != NULL) {
	gretl_matrix *V = gretl_matrix_alloc(H->rows, H->rows);

	if (V == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	gretl_matrix_qform(H, GRETL_MOD_NONE,
			   vecm->jinfo->Bvar,
			   V, GRETL_MOD_NONE);
	for (i=0; i<V->rows; i++) {
	    x = gretl_matrix_get(V, i, i);
	    vecm->jinfo->Bse->val[i] = sqrt(x);
	}
	gretl_matrix_free(vecm->jinfo->Bvar);
	vecm->jinfo->Bvar = V;
    } else {
	int k = 0;

	for (j=0; j<r; j++) {
	    /* cointegrating vector j */
	    for (i=r; i<n; i++) {
		x = gretl_matrix_get(vecm->jinfo->Bvar, k, k);
		gretl_matrix_set(vecm->jinfo->Bse, i, j, sqrt(x));
		k++;
	    }
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
    gretl_matrix_free(H);

    return err;
}

static int johansen_ll_calc (GRETL_VAR *jvar, const gretl_matrix *eigvals)
{
    gretl_matrix *S00;
    int n = jvar->neqns;
    int r = jrank(jvar);
    int h = (r > 0)? r : n;
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (S00 == NULL) {
	err = E_ALLOC;
    } else {
	double ldet = gretl_matrix_log_determinant(S00, &err);

	jvar->ll = n * (1.0 + LN_2_PI) + ldet;
	for (i=0; i<h; i++) {
	    jvar->ll += log(1.0 - eigvals->val[i]); 
	}
	jvar->ll *= -(jvar->T / 2.0);
	gretl_matrix_free(S00);
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

    k += vecm->jinfo->seasonals;
    if (jcode(vecm) >= J_UNREST_CONST) {
	k++;
    }
    if (jcode(vecm) == J_UNREST_TREND) {
	k++;
    }
    if (vecm->xlist != NULL) {
	k += vecm->xlist[0];
    } 

    /* FIXME */

#if 0
    fprintf(stderr, "k=%d; vecm->ncoeff=%d\n", k, vecm->ncoeff);
#endif

#if 1
    /* This doesn't agree with what we had before: 
       but what's wrong with it? */
    k = vecm->ncoeff;
#endif

    /* FIXME: is k right (in all cases)? */
    
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

static
int johansen_get_eigenvalues (gretl_matrix *S00,
			      const gretl_matrix *S01,
			      const gretl_matrix *S11,
			      gretl_matrix *M,
			      gretl_matrix **evals,
			      int rank)
{
    gretl_matrix *Tmp = NULL;
    int n = S11->cols;
    int err;

    err = gretl_invert_symmetric_matrix(S00);
    if (err) {
	return err;
    }

    Tmp = gretl_matrix_alloc(n, n);
    if (Tmp == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_qform(S01, GRETL_MOD_TRANSPOSE, 
		       S00, Tmp, GRETL_MOD_NONE);

    *evals = gretl_gensymm_eigenvals(Tmp, S11, M, &err);

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
    gretl_matrix *S00 = NULL;
    gretl_matrix *evals = NULL;

    int n = jvar->neqns;
    int m = gretl_matrix_cols(jvar->jinfo->S11);
    int err = 0;

    jvar->jinfo->Beta = gretl_matrix_alloc(m, m);
    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (jvar->jinfo->Beta == NULL || S00 == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = johansen_get_eigenvalues(S00, jvar->jinfo->S01, jvar->jinfo->S11,
				       jvar->jinfo->Beta, &evals, 0);
    }

    if (err) {
	pputs(prn, _("Failed to find eigenvalues\n"));
    } else {
	johansen_ll_calc(jvar, evals);
	compute_coint_test(jvar, evals, prn);

	if (!(opt & OPT_Q)) {
	    err = compute_alpha(jvar->jinfo);
	    if (!err) {
		print_beta_and_alpha(jvar, evals, n, pdinfo, prn);
		print_long_run_matrix(jvar, pdinfo, prn);
	    }
	}
    }

    gretl_matrix_free(S00);
    gretl_matrix_free(evals);

    return err;
}

static void set_beta_test_df (GRETL_VAR *jvar, const gretl_matrix *H)
{
    int r = jrank(jvar);
    int nb = jvar->jinfo->Beta->rows;

    jvar->jinfo->bdf = r * (nb - H->cols);
}

int
johansen_LR_calc (GRETL_VAR *jvar, const gretl_matrix *evals, 
		  const gretl_matrix *H, int job, PRN *prn)
{
    gretl_matrix *S00;
    double llr = 0.0;
    double ldet = 0.0;
    double T_2 = (double) jvar->T / 2.0;
    int n = jvar->neqns;
    int r = jrank(jvar);
    int h = (r > 0)? r : n;
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (S00 == NULL) {
	err = E_ALLOC;
    } else {
	ldet = gretl_matrix_log_determinant(S00, &err);
    }

    if (!err) {
	llr = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
	for (i=0; i<h; i++) {
	    pprintf(prn, _("eigenvalue %d = %g\n"), i+1, evals->val[i]);
	    llr -= T_2 * log(1.0 - evals->val[i]); 
	}
	pputc(prn, '\n');
    }

    if (S00 != NULL) {
	gretl_matrix_free(S00);
    }

    if (!err) {
	double x = 2.0 * (jvar->ll - llr);
	int nb = gretl_matrix_rows(jvar->jinfo->Beta);
	int df;

	if (job == V_BETA) {
	    df = h * (nb - H->cols);
	} else {
	    df = h * (n - H->cols);
	}

	pprintf(prn, _("Unrestricted loglikelihood (lu) = %g\n"), jvar->ll);
	pprintf(prn, _("Restricted loglikelihood (lr) = %g\n"), llr);
	pprintf(prn, "2 * (lu - lr) = %g\n", x);
	if (df > 0) {
	    pprintf(prn, _("P(Chi-Square(%d) > %g = %g\n"), df, x, 
		    chisq_cdf_comp(x, df));
	}
    }

    return err;
}

static int johansen_prep_restriction (GRETL_VAR *jvar, 
				      const gretl_matrix *R,
				      gretl_matrix **S01,
				      gretl_matrix **S11,
				      gretl_matrix **pH)
{
    gretl_matrix *H;
    int m, n = jvar->neqns;
    int err = 0;

    H = gretl_matrix_right_nullspace(R, &err);
    if (err) {
	return err;
    }

    *pH = H;
    m = gretl_matrix_cols(H);

    *S11 = gretl_matrix_alloc(m, m);
    *S01 = gretl_matrix_alloc(n, m);
    if (*S11 == NULL || *S01 == NULL) {
	return E_ALLOC;
    }

    jvar->jinfo->Beta = gretl_matrix_alloc(H->rows, jrank(jvar));
    if (jvar->jinfo->Beta == NULL) {
	return E_ALLOC;
    }

    /* calculate S11 <- H' S11 H */
    err = gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE,
			     jvar->jinfo->S11, 
			     *S11, GRETL_MOD_NONE);

    if (!err) {
	/* S01 <- S01H */
	err = gretl_matrix_multiply(jvar->jinfo->S01, H, *S01);
    } 

    return err;
}

static int use_switcher (gretlopt opt)
{
    if (opt & OPT_W) return 1;
    return (getenv("GRETL_SWITCHER") != NULL);
}

/* test for homogeneous restriction, either for a rank-1 system 
   or in common across the columns of beta (or alpha)
*/

int 
simple_restriction (GRETL_VAR *jvar,
		    const gretl_restriction_set *rset)
{
    const gretl_matrix *R, *q;
    int acols = rset_VECM_acols(rset);
    int bcols = rset_VECM_bcols(rset);
    int rcols = jvar->neqns;
    int ret = 1;

    if (acols > 0 && bcols > 0) {
	return 0;
    }

    if (bcols > 0) {
	rcols += restricted(jvar);
	R = rset_get_R_matrix(rset);
	q = rset_get_q_matrix(rset);
    } else {
	R = rset_get_Ra_matrix(rset);
	q = rset_get_qa_matrix(rset);
    }

    if (!gretl_is_zero_matrix(q)) {
	/* non-homogeneous */
	ret = 0;
    } else if (R->cols > rcols) {
	/* not common to all columns */
	ret = 0;
    }

    return ret;
}

/* driver for VECM estimation subject to "general" restrictions
   on beta and/or alpha */

static int johansen_estimate_general (GRETL_VAR *jvar, 
				      const gretl_restriction_set *rset,
				      const double **Z, const DATAINFO *pdinfo, 
				      gretlopt opt, PRN *prn)
{
    const gretl_matrix *R, *q;
    gretlopt vopt = OPT_F;
    int err;

    if (use_switcher(opt)) {
	vopt |= OPT_W;
    }

    err = general_vecm_analysis(jvar, rset, pdinfo, vopt, prn);

    if (!err && !(vopt & OPT_W)) {
	int vflag;

	vflag = (rset_VECM_acols(rset) > 0)? 
	    NET_OUT_ALPHA : ESTIMATE_ALPHA;
	err = build_VECM_models(jvar, Z, pdinfo, vflag, 0);
    }

    if (!err) {
	err = vecm_ll_stats(jvar);
    }

    if (!err && rset_VECM_bcols(rset) > 0) {
	R = rset_get_R_matrix(rset);
	q = rset_get_q_matrix(rset);

	gretl_matrix_free(jvar->jinfo->R);
	gretl_matrix_free(jvar->jinfo->q);

	jvar->jinfo->R = gretl_matrix_copy(R);
	jvar->jinfo->q = gretl_matrix_copy(q);

	if (jvar->jinfo->R == NULL || jvar->jinfo->q == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && rset_VECM_acols(rset) > 0) {
	R = rset_get_Ra_matrix(rset);
        q = rset_get_qa_matrix(rset);

	gretl_matrix_free(jvar->jinfo->Ra);
	gretl_matrix_free(jvar->jinfo->qa);

	jvar->jinfo->Ra = gretl_matrix_copy(R);
	jvar->jinfo->qa = gretl_matrix_copy(q);

	if (jvar->jinfo->Ra == NULL || jvar->jinfo->qa == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

/* Public entry point for VECM estimation.  If rset != NULL we're
   imposing a restriction on the cointegrating vectors (or possibly
   alpha); and in that case how we proceed depends on whether the
   restrictions can be handled by the modified eigen-system approach.
*/

int johansen_estimate (GRETL_VAR *jvar, 
		       const gretl_restriction_set *rset,
		       const double **Z, const DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    const gretl_matrix *R = NULL;

    gretl_matrix *H = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *S00 = NULL;
    gretl_matrix *S01 = NULL;
    gretl_matrix *S11 = NULL;
    gretl_matrix *evals = NULL;

    int genrest = 0; /* doing general restriction? */
    int rank = jrank(jvar);
    int m, err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting johansen_estimate(), rset = %p\n\n",
	    rset);
#endif

    if (rset != NULL) {
	genrest = !simple_restriction(jvar, rset);
    }

    if (rset_VECM_acols(rset) > 0) {
	pprintf(prn, "\"full\" restriction on VECM via alpha: "
		"not handled yet\n");
	return E_NOTIMP;
    }

    if (rset != NULL && !genrest) {
	R = rset_get_R_matrix(rset);
	err = johansen_prep_restriction(jvar, R, &S01, &S11, &H);
	if (err) {
	    goto bailout;
	} 
	m = gretl_matrix_cols(H);
    } else {
	S11 = jvar->jinfo->S11;
	S01 = jvar->jinfo->S01;
	m = gretl_matrix_cols(S11);
    }

    M = gretl_matrix_alloc(m, m);
    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (M == NULL || S00 == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = johansen_get_eigenvalues(S00, S01, S11, M, &evals, rank);
    }

    if (err) {
	pputs(prn, _("Failed to find eigenvalues\n"));
	goto bailout;
    } 

#if JDEBUG
    gretl_matrix_print(M, "raw eigenvector(s)");
#endif

    if (H != NULL) {
	err = gretl_matrix_multiply(H, M, jvar->jinfo->Beta);
	set_beta_test_df(jvar, H);
    } else {
	jvar->jinfo->Beta = M;
	M = NULL;
    }

    if (!err) {
	err = johansen_ll_calc(jvar, evals);
    }

    if (!genrest) {
	int do_stderrs = rank < jvar->neqns;

	if (!err) {
	    err = normalize_beta(jvar, R, &do_stderrs); 
	}
	if (!err) {
	    err = build_VECM_models(jvar, Z, pdinfo, ESTIMATE_ALPHA, 0);
	}
	if (!err) {
	    err = compute_omega(jvar);
	}
	if (!err && do_stderrs) {
	    /* FIXME case where R != NULL */
	    err = beta_variance(jvar, R);
	}
	if (!err) {
	    err = gretl_VAR_do_error_decomp(jvar->S, jvar->C);
	}
	if (!err) {
	    err = vecm_ll_stats(jvar);
	}
	if (!err && R != NULL) {
	    jvar->jinfo->R = gretl_matrix_copy(R);
	}
    }

 bailout:    

    gretl_matrix_free(H);
    gretl_matrix_free(M);
    gretl_matrix_free(S00);
    gretl_matrix_free(evals);

    if (H != NULL) {
	gretl_matrix_free(S01);
	gretl_matrix_free(S11);
    }

    if (!err && genrest) {
	err = johansen_estimate_general(jvar, rset, Z, pdinfo, opt, prn);
    }

    return err;
}

/* Simplified version of the Johansen procedure, to be called in
   the process of computing bootstrap confidence intervals for
   impulse response functions.  We just have to do enough to
   generate the VAR representation.
*/

/* FIXME case of restricted beta */

int 
johansen_boots_round (GRETL_VAR *jvar, const double **Z, 
		      const DATAINFO *pdinfo,
		      int iter)
{
    gretl_matrix *M = NULL;
    gretl_matrix *evals = NULL;
    int m = gretl_matrix_cols(jvar->jinfo->S11);
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting johansen_bootstrap_round()\n\n");
#endif

    M = gretl_matrix_alloc(m, m);
    if (M == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = johansen_get_eigenvalues(jvar->jinfo->S00, jvar->jinfo->S01, 
				   jvar->jinfo->S11, M, &evals, 
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
	    err = normalize_beta(jvar, NULL, NULL); 
	}
	if (!err) {
	    err = build_VECM_models(jvar, Z, pdinfo, ESTIMATE_ALPHA, 
				    iter);
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

void print_beta_alpha_Pi (GRETL_VAR *jvar,
			  const DATAINFO *pdinfo,
			  PRN *prn)
{
    int r = jrank(jvar);

    print_beta_or_alpha(jvar, r, pdinfo, prn, V_BETA, 0);
    print_beta_or_alpha(jvar, r, pdinfo, prn, V_ALPHA, 0);
    pputc(prn, '\n');
    print_long_run_matrix(jvar, pdinfo, prn);
}

/* compute and print beta, alpha and alpha*beta', in the context where
   we've tested a (common, homogeneous) restriction on beta,
   represented by H, and verbose output has been requested.
*/

static int show_beta_alpha_etc (GRETL_VAR *jvar,
				const gretl_matrix *H,
				const gretl_matrix *M,
				const DATAINFO *pdinfo,
				PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    int err = 0;

    gretl_matrix_multiply_mod(H, GRETL_MOD_NONE,
			      M, GRETL_MOD_NONE,
			      jv->Beta, GRETL_MOD_NONE);

    if (jv->rank == 1) { 
	/* and if r > 1? */
	double den = jv->Beta->val[0];

	if (!floateq(den, 0.0)) {
	    gretl_matrix_divide_by_scalar(jv->Beta, den);
	}
    }

    if (!err) {
	err = compute_alpha(jv);
    }

    if (!err) {
	print_beta_alpha_Pi(jvar, pdinfo, prn);
    }

    return err;
}

/* test for a common, homogeneous restriction on beta (only) */

static int vecm_beta_test (GRETL_VAR *jvar, 
			   const gretl_restriction_set *rset,
			   const DATAINFO *pdinfo, 
			   gretlopt opt,
			   PRN *prn)
{
    const gretl_matrix *R;
    gretl_matrix *H = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *S11 = NULL;
    gretl_matrix *S01 = NULL;
    gretl_matrix *S00 = NULL;
    gretl_matrix *evals = NULL;
    int m, n, rank;
    int err = 0;

    R = rset_get_R_matrix(rset);
    H = gretl_matrix_right_nullspace(R, &err);

    if (err) {
	return err;
    }

    n = jvar->neqns;
    rank = jrank(jvar);
    m = gretl_matrix_cols(H);

    M = gretl_matrix_alloc(m, m);
    S11 = gretl_matrix_alloc(m, m);
    S01 = gretl_matrix_alloc(n, m);
    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (M == NULL || S11 == NULL || S01 == NULL || S00 == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    pputs(prn, "\nTest of restrictions on cointegrating relations\n\n");

    if (opt & OPT_V) {
	gretl_matrix_print_to_prn(H, "Restriction matrix, H", prn);
    }

    /* calculate S11 <- H' S11 H */
    err = gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE,
			     jvar->jinfo->S11, S11, 
			     GRETL_MOD_NONE);

    if (opt & OPT_V) {
	gretl_matrix_print_to_prn(S11, "H'*S11*H", prn);
    }

    if (!err) {
	/* S01 <- S01*H */
	err = gretl_matrix_multiply(jvar->jinfo->S01, H, S01);
    }

    if (opt & OPT_V) {
	gretl_matrix_print_to_prn(S01, "S01*H", prn);
    }

    if (!err) {
	err = johansen_get_eigenvalues(S00, S01, S11, M, &evals, rank);
    }

    if (!err) {
	if (opt & OPT_V) {
	    gretl_matrix_print_to_prn(M, "M", prn);
	}
	johansen_LR_calc(jvar, evals, H, V_BETA, prn);
    } 

    if (!err && (opt & OPT_V)) {
	show_beta_alpha_etc(jvar, H, M, pdinfo, prn);
    }

 bailout:    

    gretl_matrix_free(H);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);
    gretl_matrix_free(S00);
    gretl_matrix_free(S11);
    gretl_matrix_free(S01);

    return err;
}

/* Test of linear restrictions on the cointegrating relations in a
   VECM.  If the restrictions are "simple" (homogeneous and in common)
   we do the test using the eigen-system approach.  If they are
   "general" restrictions we hand off to the specialized machinery in
   jrestrict.c.
*/

int vecm_test_restriction (GRETL_VAR *jvar, 
			   const gretl_restriction_set *rset,
			   const DATAINFO *pdinfo, 
			   gretlopt opt,
			   PRN *prn)
{
    int acols = rset_VECM_acols(rset);
    int err = 0;

#if 0
    if (acols > 0 && !use_switcher(opt)) {
	pprintf(prn, "alpha restriction: not handled yet\n");
	return E_NOTIMP;
    } 
#endif

#if 0
    if (alpha_restricted_VECM(jvar) && bcols > 0) {
	pprintf(prn, "Beta restriction for an alpha-restricted VECM: "
		"not handled yet\n");
	return E_NOTIMP;
    }
#endif	

    if (simple_restriction(jvar, rset)) {
	if (acols > 0) {
	    err = vecm_alpha_test(jvar, rset, pdinfo, opt, prn);
	} else {
	    err = vecm_beta_test(jvar, rset, pdinfo, opt, prn);
	}
    } else {
	/* "general" restriction */
	if (use_switcher(opt)) {
	    opt |= OPT_W;
	}
	err = general_vecm_analysis(jvar, rset, pdinfo, opt, prn);
    }

    return err;
}
