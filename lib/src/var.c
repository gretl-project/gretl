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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/* var.c - vector autoregressions */  

#include "libgretl.h" 
#include "var.h"  
#include "libset.h"

/* in transforms.c */
extern int 
real_list_laggenr (const int *list, double ***pZ, DATAINFO *pdinfo,
		   int maxlag, int **lagnums);

static gretl_matrix *
gretl_VAR_get_fcast_decomp (GRETL_VAR *var, int targ, int periods);

static gretl_matrix *irf_bootstrap (const GRETL_VAR *var, 
				    int targ, int shock, int periods,
				    const double **Z, 
				    const DATAINFO *pdinfo);

#define VAR_DEBUG 0

struct var_resids {
    int *levels_list;
    gretl_matrix *u;
    gretl_matrix *v;
    int t1, t2;
};

struct var_lists {
    int *detvars;
    int *stochvars;
    int *reglist;
    int *testlist;
    int **lagvlist;
};

enum {
    ADF_EG_TEST   = 1 << 0,
    ADF_PRINT_ACK = 1 << 1
} adf_flags;

#define TREND_FAILED 9999

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags, 
			  PRN *prn);

static void var_resids_init (struct var_resids *vr)
{
    vr->levels_list = NULL;
    vr->u = vr->v = NULL;
    vr->t1 = vr->t2 = 0;
}

static void var_resids_free (struct var_resids *vr)
{
    free(vr->levels_list);
    gretl_matrix_free(vr->u);
    gretl_matrix_free(vr->v);
}

static void pad_var_coeff_matrix (gretl_matrix *A, int neqns, int order)
{
    int i, j;
    int rowmax = neqns * order;

    for (i=neqns; i<rowmax; i++) {
	for (j=0; j<rowmax; j++) {
	    gretl_matrix_set(A, i, j, (j == i - neqns)? 1.0 : 0.0);
	}
    }
}

static MODEL **allocate_VAR_models (int neqns)
{
    MODEL **models = malloc(neqns * sizeof *models);
    int i, j;

    if (models != NULL) {
	for (i=0; i<neqns; i++) {
	    models[i] = gretl_model_new();
	    if (models[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(models[i]);
		}
		free(models);
		models = NULL;
	    }
	}
    }

    return models;
}

int gretl_VAR_add_coeff_matrix (GRETL_VAR *var)
{
    int n = var->neqns;
    int err = 0;

    if (var->order > 1) {
	n *= var->order;
    } 

    var->A = gretl_matrix_alloc(n, n);
    if (var->A == NULL) {
	err = E_ALLOC;
    } else {
	pad_var_coeff_matrix(var->A, var->neqns, var->order);
    }

    return err;
}

int gretl_VAR_add_C_matrix (GRETL_VAR *var)
{
    int n = var->neqns;
    int err = 0;

    if (var->order > 1) {
	n *= var->order;
    } 

    var->C = gretl_matrix_alloc(n, var->neqns);
    if (var->C == NULL) {
	err = 1;
    } else {
	gretl_matrix_zero(var->C);
    }

    return err;
}

static GRETL_VAR *gretl_VAR_new (int neqns, int order, const DATAINFO *pdinfo)
{
    GRETL_VAR *var;
    int err = 0;

    if (neqns == 0 || order == 0) {
	return NULL;
    }

    var = malloc(sizeof *var);
    if (var == NULL) {
	return NULL;
    }

    var->ci = VAR;
    var->err = 0;
    var->t1 = var->t2 = var->T = 0;

    var->neqns = neqns;
    var->order = order;
    var->ncoeff = 0;

    var->A = NULL;
    var->lambda = NULL;
    var->E = NULL;
    var->C = NULL;
    var->S = NULL;
    var->F = NULL;

    var->models = NULL;
    var->Fvals = NULL;
    var->name = NULL;

    var->ll = var->ldet = NADBL;
    var->AIC = var->BIC = NADBL;
    var->LR = NADBL;

    var->jinfo = NULL;

    err = gretl_VAR_add_coeff_matrix(var);

    if (!err) {
	err = gretl_VAR_add_C_matrix(var);
    }

    if (!err) {
	var->models = allocate_VAR_models(neqns);
	if (var->models == NULL) {
	    err = 1;
	} 
    } 

    if (!err) {
	int m = neqns * neqns + neqns;
	
	var->Fvals = malloc(m  * sizeof *var->Fvals);
	if (var->Fvals == NULL) {
	    err = 1;
	}
    }

    if (err) {
	gretl_VAR_free(var);
	var = NULL;
    }

    return var;
}

static void destroy_VAR_models (MODEL **models, int neqns)
{
    int i;

    if (models != NULL) {
	for (i=0; i<neqns; i++) {
	    clear_model(models[i]);
	    free(models[i]);
	}
	free(models);
    }
}

static void johansen_info_free (JohansenInfo *jv)
{
    free(jv->list);
    free(jv->difflist);
    free(jv->biglist);

    gretl_matrix_free(jv->u);
    gretl_matrix_free(jv->v);
    gretl_matrix_free(jv->w);

    gretl_matrix_free(jv->Suu);
    gretl_matrix_free(jv->Svv);
    gretl_matrix_free(jv->Suv);

    gretl_matrix_free(jv->Beta);
    gretl_matrix_free(jv->Alpha);
    gretl_matrix_free(jv->Bse);

    free(jv);
}

void gretl_VAR_free (GRETL_VAR *var)
{
    if (var == NULL) return;

    gretl_matrix_free(var->A);
    gretl_matrix_free(var->lambda);
    gretl_matrix_free(var->E);
    gretl_matrix_free(var->C);
    gretl_matrix_free(var->S);
    gretl_matrix_free(var->F);

    free(var->Fvals);
    free(var->name);

    if (var->models != NULL) {
	destroy_VAR_models(var->models, var->neqns);
    }

    if (var->jinfo != NULL) {
	johansen_info_free(var->jinfo);
    }

    free(var);
}

void gretl_VAR_free_unnamed (GRETL_VAR *var)
{
    if (var == NULL) return;

    if (var->name == NULL || *var->name == '\0') {
	gretl_VAR_free(var);
    }
}

/* FIXME: the following function has to be re-done for VECMs, since in
   the case of VECMs, var->models contains the VECM models where the
   dependent variable is the first difference.  This should be
   converted into a forecast for the level, I think.
*/

static int
gretl_VAR_add_forecast (GRETL_VAR *var, int t1, int t2, const double **Z, 
			const DATAINFO *pdinfo, gretlopt opt)
{
    const MODEL *pmod;
    gretl_matrix *F;
    double fti, xti;
    int i, j, k, s, t;
    int nf, ns, lag, vj, m;
    int staticfc, fcols;

    pmod = var->models[0];

    nf = t2 - t1 + 1;

    staticfc = (opt & OPT_S);
    if (staticfc) {
	fcols = var->neqns;
    } else {
	fcols = 2 * var->neqns;
    }

    /* rows = number of forecast periods; cols = 1 to hold forecast
       for each variable, plus 1 to hold variance for each variable
       if forecast is dynamic.
    */
    F = gretl_matrix_alloc(nf, fcols);
    if (F == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_zero(F);

    ns = var->order * var->neqns;

    for (t=t1; t<=t2; t++) {
	int miss = 0;

	s = t - t1;
	for (i=0; i<var->neqns; i++) {
	    pmod = var->models[i];
	    fti = 0.0;
	    lag = 1;
	    k = 0;
	    for (j=0; j<pmod->ncoeff; j++) {
		vj = pmod->list[j + 2];
		if (j < ns + pmod->ifc && vj > 0) {
		    /* stochastic var */
		    if (staticfc || s - lag < 0) {
			/* pre-forecast value */
			m = (j - pmod->ifc) / var->order;
			vj = var->models[m]->list[1];
			if (t - lag < 0) {
			    xti = NADBL;
			} else {
			    xti = Z[vj][t-lag];
			}
			if (na(xti)) {
			    miss = 1;
			}
		    } else {
			/* prior forecast value */
			xti = gretl_matrix_get(F, s - lag, k);
		    }
		    lag++;
		    if (lag > var->order) {
			lag = 1;
			k++;
		    }
		} else {
		    /* deterministic var: value from dataset */
		    xti = Z[vj][t];
		    if (na(xti)) {
			miss = 1;
		    }
		}
		if (miss) {
		    fti = NADBL;
		} else {
		    fti += pmod->coeff[j] * xti;
		}
	    }
	    gretl_matrix_set(F, s, i, fti);
	}
    }

    /* now get variances, if not static */
    if (!staticfc) {
	double vti;
	int totcol;

	for (i=0; i<var->neqns; i++) {
	    gretl_matrix *vd;

	    vd = gretl_VAR_get_fcast_decomp(var, i, nf);
	    if (vd != NULL) {
		totcol = gretl_matrix_cols(vd) - 1;
		for (s=0; s<nf; s++) {
		    vti = gretl_matrix_get(vd, s, totcol);
		    gretl_matrix_set(F, s, var->neqns + i, vti);
		}
		gretl_matrix_free(vd);
	    } else {
		for (s=0; s<nf; s++) {
		    gretl_matrix_set(F, s, var->neqns + i, NADBL);
		}
	    }
	}
    }

#if 0
    gretl_matrix_print(F, "var->F", NULL);
#endif

    gretl_matrix_set_int(F, t1);

    var->F = F;

    return 0;
}

const gretl_matrix *
gretl_VAR_get_forecast_matrix (GRETL_VAR *var, int t1, int t2, const double **Z, 
			       const DATAINFO *pdinfo, gretlopt opt)
{
    if (var->F != NULL) {
	int ncols, nf = t2 - t1 + 1;
	int ft1 = gretl_matrix_get_int(var->F);

	ncols = (opt & OPT_S)? var->neqns: 2 * var->neqns;

	if (nf == gretl_matrix_rows(var->F) && t1 == ft1 && 
	    ncols == gretl_matrix_cols(var->F)) {
	    ; /* already done, fine */
	} else {
	    gretl_matrix_free(var->F);
	    var->F = NULL;
	}
    }

    if (var->F == NULL) {
	gretl_VAR_add_forecast(var, t1, t2, Z, pdinfo, opt);
    }

    return var->F;
}

const gretl_matrix *
gretl_VAR_get_residual_matrix (const GRETL_VAR *var)
{
    return var->E;
}

static int
get_moments (const gretl_matrix *M, int row, double *skew, double *kurt)
{
    int j, n = gretl_matrix_cols(M);
    double xi, xbar, dev, var;
    double s = 0.0;
    double s2 = 0.0;
    double s3 = 0.0;
    double s4 = 0.0;
    int err = 0;
    
    for (j=0; j<n; j++) {
	s += gretl_matrix_get(M, row, j);
    }

    xbar = s / n;

    for (j=0; j<n; j++) {
	xi = gretl_matrix_get(M, row, j);
	dev = xi - xbar;
	s2 += dev * dev;
	s3 += pow(dev, 3);
	s4 += pow(dev, 4);
    }

    var = s2 / n;

    if (var > 0.0) {
	/* if variance is effectively zero, these should be undef'd */
	*skew = (s3 / n) / pow(s2 / n, 1.5);
	*kurt = ((s4 / n) / pow(s2 / n, 2));
    } else {
	*skew = *kurt = NADBL;
	err = 1;
    }

    return err;
}

static int 
real_VAR_normality_test (const gretl_matrix *E, const gretl_matrix *Sigma, 
			 PRN *prn)
{
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *tmp = NULL;

    /* convenience pointers: do not free! */
    gretl_matrix *H;
    gretl_vector *Z1;
    gretl_vector *Z2;

    double *evals = NULL;
    double x, skew, kurt;
    double X2 = NADBL;
    int n, p;
    int i, j;
    int err = 0;

    if (E == NULL) {
	err = 1;
	goto bailout;
    }

    p = gretl_matrix_cols(E);
    n = gretl_matrix_rows(E);

    S = gretl_matrix_copy(Sigma);
    V = gretl_vector_alloc(p);
    C = gretl_matrix_alloc(p, p);
    X = gretl_matrix_copy_transpose(E);
    R = gretl_matrix_alloc(p, n);
    tmp = gretl_matrix_alloc(p, p);

    if (S == NULL || V == NULL || C == NULL || X == NULL || 
	R == NULL || tmp == NULL) {
	err = 1;
	goto bailout;
    }

    for (i=0; i<p; i++) {
	x = gretl_matrix_get(S, i, i);
	gretl_vector_set(V, i, 1.0 / sqrt(x));
    }

    err = gretl_matrix_diagonal_sandwich(V, S, C);
    if (err) {
	goto bailout;
    }

    gretl_matrix_print(C, "\nResidual correlation matrix, C", prn);

    evals = gretl_symmetric_matrix_eigenvals(C, 1);
    if (evals == NULL) {
	goto bailout;
    }

    pputs(prn, "Eigenvalues of the correlation matrix:\n\n");
    for (i=0; i<p; i++) {
	pprintf(prn, " %10g\n", evals[i]);
    }
    pputc(prn, '\n');

    /* C should now contain eigenvectors of the original C:
       relabel as 'H' for perspicuity */
    H = C;
#if 0
    gretl_matrix_print(H, "Eigenvectors, H", prn);
#endif
    gretl_matrix_copy_values(tmp, H);

    /* make "tmp" into $H \Lambda^{-1/2}$ */
    for (i=0; i<p; i++) {
	for (j=0; j<p; j++) {
	    x = gretl_matrix_get(tmp, i, j);
	    x *= 1.0 / sqrt(evals[j]);
	    gretl_matrix_set(tmp, i, j, x);
	}
    }

    /* make S into $H \Lambda^{-1/2} H'$ */
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
			      H, GRETL_MOD_TRANSPOSE,
			      S);

    /* compute VX', in X (don't need to subtract means, because these
       are OLS residuals)
    */
    for (i=0; i<p; i++) {
	for (j=0; j<n; j++) {
	    x = gretl_matrix_get(X, i, j);
	    x *= gretl_vector_get(V, i);
	    gretl_matrix_set(X, i, j, x);
	}
    }

    /* finally, compute $R' = H  \Lambda^{-1/2} H' V X'$ 
       Doornik and Hansen, 1994, section 3 */
    gretl_matrix_multiply(S, X, R);

    /* Z_1 and Z_2 are p-vectors: use existing storage */
    Z1 = V;
    Z2 = gretl_matrix_reuse(tmp, p, 1);

    for (i=0; i<p && !err; i++) {
	get_moments(R, i, &skew, &kurt);
	if (na(skew) || na(kurt)) {
	    err = 1;
	} else {
	    double z1i = dh_root_b1_to_z1(skew, n);
	    double z2i = dh_b2_to_z2(skew * skew, kurt, n);

	    gretl_vector_set(Z1, i, z1i);
	    gretl_vector_set(Z2, i, z2i);
	}
    }
	
    if (!err) {
	X2 = gretl_vector_dot_product(Z1, Z1, &err);
	X2 += gretl_vector_dot_product(Z2, Z2, &err);
    }

    if (na(X2)) {
	pputs(prn, "Calculation of test statistic failed\n");
    } else {
	pputs(prn, "Test for multivariate normality of residuals\n");
	pprintf(prn, "Doornik-Hansen Chi-square(%d) = %g, ", 2 * p, X2);
	pprintf(prn, "with p-value = %g\n", chisq(X2, 2 * p));
    }

 bailout:

    gretl_matrix_free(S);
    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(X);
    gretl_matrix_free(R);
    gretl_matrix_free(tmp);

    free(evals);

    return err;
}

int gretl_VAR_normality_test (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->E == NULL || var->S == NULL) {
	err = 1;
    } else {
	err = real_VAR_normality_test(var->E, var->S, prn);
    }

    return err;
}

int gretl_VAR_autocorrelation_test (GRETL_VAR *var, int order, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    int i, err = 0;

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "Equation %d:\n", i + 1);
	err = autocorr_test(var->models[i], order, pZ, pdinfo,
			    OPT_Q | OPT_S, prn);
	gretl_model_test_print(var->models[i], 0, prn);
	gretl_model_destroy_tests(var->models[i]);
    }

    return err;
}

int gretl_VAR_arch_test (GRETL_VAR *var, int order, 
			 double ***pZ, DATAINFO *pdinfo,
			 PRN *prn)
{
    int i, err = 0;

    for (i=0; i<var->neqns && !err; i++) {
	pprintf(prn, "Equation %d:\n", i + 1);
	err = arch_test_simple(var->models[i], order, pZ, pdinfo,
			       prn);
	gretl_model_test_print(var->models[i], 0, prn);
	gretl_model_destroy_tests(var->models[i]);
    }

    return err;
}

int gretl_VAR_do_error_decomp (int g, const gretl_matrix *S,
			       gretl_matrix *C)
{
    gretl_matrix *tmp = NULL;
    double x;
    int i, j, err = 0;

    /* copy cross-equation covariance matrix */
    tmp = gretl_matrix_copy(S);
    if (tmp == NULL) {
	err = E_ALLOC;
    }

    /* lower-triangularize and decompose */
    if (!err) {
	for (i=0; i<g-1; i++) {
	    for (j=i+1; j<g; j++) {
		gretl_matrix_set(tmp, i, j, 0.0);
	    }
	}
	err = gretl_matrix_cholesky_decomp(tmp);
    }

    /* write the decomposition into the C matrix */
    if (!err) {
	for (i=0; i<g; i++) {
	    for (j=0; j<g; j++) {
		x = gretl_matrix_get(tmp, i, j);
		gretl_matrix_set(C, i, j, x);
	    }
	}
    }

    if (tmp != NULL) {
	gretl_matrix_free(tmp);
    }

    return err;
}

int gretl_VAR_get_variable_number (const GRETL_VAR *var, int k)
{
    if (var->ci == VECM) {
	return var->jinfo->list[k + 1];
    } else {
	return (var->models[k])->list[1];
    }
}

int gretl_VAR_get_n_equations (const GRETL_VAR *var)
{
    return var->neqns;
}

int gretl_VAR_get_t1 (const GRETL_VAR *var)
{
    return var->t1;
}

int gretl_VAR_get_t2 (const GRETL_VAR *var)
{
    return var->t2;
}

const MODEL *gretl_VAR_get_model (const GRETL_VAR *var, int i)
{
    if (i < var->neqns) {
	return var->models[i];
    } else {
	return NULL;
    }
}

static int periods_from_pd (int pd)
{
    int periods = 10;

    if (pd == 4) {
	/* quarterly: try 5 years */
	periods = 20;
    } else if (pd == 12) {
	/* monthly: two years */
	periods = 24;
    } else if (pd == 7 || pd == 6 || pd == 5) {
	/* daily: three weeks */
	periods = 3 * pd;
    } 

    return periods;
}

int default_VAR_horizon (const DATAINFO *pdinfo)
{
    int h = get_VAR_horizon();

    if (h <= 0) {
	h = periods_from_pd(pdinfo->pd);
    }

    return h;
}

static gretl_matrix *
gretl_VAR_get_point_responses (GRETL_VAR *var, int targ, int shock,
			       int periods) 
{
    int rows = var->neqns * var->order;
    gretl_matrix *rtmp = NULL;
    gretl_matrix *ctmp = NULL;
    gretl_matrix *resp = NULL;
    double rt;
    int t, err = 0;

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return NULL;
    }  

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return NULL;
    } 

    if (periods <= 0) {
	fprintf(stderr, "Invalid number of periods\n");
	return NULL;
    }

    resp = gretl_matrix_alloc(periods, 1);
    if (resp == NULL) {
	return NULL;
    }

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) {
	gretl_matrix_free(resp);
	return NULL;
    }

    ctmp = gretl_matrix_alloc(rows, var->neqns);
    if (ctmp == NULL) {
	free(resp);
	gretl_matrix_free(rtmp);
	return NULL;
    }

    for (t=0; t<periods && !err; t++) {
	if (t == 0) {
	    /* calculate initial estimated responses */
	    err = gretl_matrix_copy_values(rtmp, var->C);
	} else {
	    /* calculate further estimated responses */
	    err = gretl_matrix_multiply(var->A, rtmp, ctmp);
	    gretl_matrix_copy_values(rtmp, ctmp);
	}

	if (!err) {
	    rt = gretl_matrix_get(rtmp, targ, shock);
	    gretl_matrix_set(resp, t, 0, rt);
	}
    }

    gretl_matrix_free(rtmp);
    gretl_matrix_free(ctmp);

    return resp;    
}

/**
 * gretl_VAR_get_impulse_response:
 * @var: pointer to VAR struct.
 * @targ: index of the target or response variable.
 * @shock: index of the source or shock variable.
 * @periods: number of periods over which to compute the response.
 * @Z: data array (or %NULL).
 * @pdinfo: dataset information.
 *
 * Computes the response of @targ to a perturbation of @shock
 * in the context of @var: @targ and @shock are zero-based indices 
 * relative to the structure of @var.  For example if @targ = 0 and 
 * @shock = 1, we compute the response of the dependent variable in 
 * the first VAR equation to a perturbation of the variable that
 * appears as dependent in the second VAR equation.
 *
 * If @Z is %NULL, the response matrix returned is a column vector 
 * of length @periods, giving the point estimate of the response 
 * function.  If @Z is not %NULL, the response matrix returned
 * has three columns, containing the point estimate, the 0.025
 * and the 0.975 quantile, where the quantiles are based on 999
 * bootstrap replications, with resampling of the original
 * residuals with replacement.
 *
 * Returns: matrix containing the estimated impulse responses,
 * with or without a confidence interval depending on whether
 * or not @Z is provided.
 */

gretl_matrix *
gretl_VAR_get_impulse_response (GRETL_VAR *var, 
				int targ, int shock, int periods,
				const double **Z,
				const DATAINFO *pdinfo)
{
    gretl_matrix *point = NULL;
    gretl_matrix *full = NULL;
    gretl_matrix *ret = NULL;
    int i;

    point = gretl_VAR_get_point_responses(var, targ, shock, periods);

    if (Z == NULL) {
	/* no data matrix given: just return point estimate */
	ret = point;
    } else if (point != NULL) {
	full = irf_bootstrap(var, targ, shock, periods, Z, pdinfo);
	if (full != NULL) {
	    double p;

	    for (i=0; i<periods; i++) {
		p = gretl_matrix_get(point, i, 0);
		gretl_matrix_set(full, i, 0, p);
	    }
	}
	gretl_matrix_free(point);
	ret = full;
    }

    return ret;
}

gretl_matrix *
gretl_VAR_get_fcast_decomp (GRETL_VAR *var, int targ, int periods) 
{
    int i, t;
    int rows = var->neqns;
    gretl_matrix *ctmp = NULL, *idx = NULL, *vtmp = NULL;
    gretl_matrix *cic = NULL, *vt = NULL;
    gretl_matrix *vd = NULL;
    int err = 0;

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return NULL;
    } 

    if (periods <= 0) {
	fprintf(stderr, "Invalid number of periods\n");
	return NULL;
    }

    if (var->order > 1) {
	rows *= var->order;
    }

    vd = gretl_matrix_alloc(periods, var->neqns + 1);
    ctmp = gretl_matrix_alloc(var->neqns, rows);
    idx = gretl_matrix_alloc(var->neqns, var->neqns); 
    cic = gretl_matrix_alloc(rows, rows);
    vt = gretl_matrix_alloc(rows, rows);
    vtmp = gretl_matrix_alloc(rows, rows);

    if (vd == NULL || ctmp == NULL || idx == NULL ||
	cic == NULL || vt == NULL || vtmp == NULL) {
	gretl_matrix_free(vd);
	gretl_matrix_free(ctmp);
	gretl_matrix_free(idx);
	gretl_matrix_free(cic);
	gretl_matrix_free(vt);
	gretl_matrix_free(vtmp);
	return NULL;
    }

    for (i=0; i<var->neqns; i++) {
	double vti;

	/* make appropriate index matrix */
	gretl_matrix_zero(idx);
	gretl_matrix_set(idx, i, i, 1.0);

	for (t=0; t<periods && !err; t++) {

	    if (t == 0) {
		/* calculate initial variances */
		err = gretl_matrix_multiply_mod(idx, GRETL_MOD_NONE,
						var->C, GRETL_MOD_TRANSPOSE,
						ctmp);
		err = gretl_matrix_multiply(var->C, ctmp, cic);
		gretl_matrix_copy_values(vt, cic);
	    } else {
		/* calculate further variances */
		err = gretl_matrix_multiply_mod(vt, GRETL_MOD_NONE,
						var->A, GRETL_MOD_TRANSPOSE,
						vtmp);
		err = gretl_matrix_multiply(var->A, vtmp, vt);
		gretl_matrix_add_to(vt, cic);
	    }

	    if (err) break;

	    vti = gretl_matrix_get(vt, targ, targ);
	    gretl_matrix_set(vd, t, i, vti);
	}
    }

    /* normalize variance contributions as percentage shares */
    for (t=0; t<periods && !err; t++) {
	double vtot = 0.0;
	double vi;

	for (i=0; i<var->neqns; i++) {
	    vtot += gretl_matrix_get(vd, t, i);
	}

	for (i=0; i<var->neqns; i++) {
	    vi = gretl_matrix_get(vd, t, i);
	    gretl_matrix_set(vd, t, i, 100.0 * vi / vtot);
	}

	gretl_matrix_set(vd, t, var->neqns, sqrt(vtot));
    }

    gretl_matrix_free(ctmp);
    gretl_matrix_free(idx);
    gretl_matrix_free(cic);
    gretl_matrix_free(vt);
    gretl_matrix_free(vtmp);

    return vd;
}

static int gettrend (double ***pZ, DATAINFO *pdinfo, int square)
{
    int index;
    int t, n = pdinfo->n, v = pdinfo->v;
    double x;

    if (square) {
	index = varindex(pdinfo, "timesq");
    } else {
	index = varindex(pdinfo, "time");
    }

    if (index < v) {
	return index;
    }
    
    if (dataset_add_series(1, pZ, pdinfo)) {
	return TREND_FAILED;
    }

    for (t=0; t<n; t++) {
	x = (double) t + 1;
	(*pZ)[v][t] = (square)? x * x : x;
    }

    if (square) {
	strcpy(pdinfo->varname[v], "timesq");
	strcpy(VARLABEL(pdinfo, v), _("squared time trend variable"));
    } else {
	strcpy(pdinfo->varname[v], "time");
	strcpy(VARLABEL(pdinfo, v), _("time trend variable"));
    }
	    
    return index;
}

static void var_lists_free (struct var_lists *vl)
{
    if (vl->lagvlist != NULL && vl->stochvars != NULL) {
	int i, ns = vl->stochvars[0];

	for (i=0; i<ns; i++) {
	    free(vl->lagvlist[i]);
	}
	free(vl->lagvlist);
    }

    free(vl->detvars);
    free(vl->stochvars);
    free(vl->reglist);
    free(vl->testlist);
}

static int **lagvlist_construct (int nstoch, int order)
{
    int **lvlist;
    int i, j;

    lvlist = malloc(nstoch * sizeof *lvlist);
    if (lvlist == NULL) {
	return NULL;
    }

    for (i=0; i<nstoch; i++) {
	lvlist[i] = malloc((order + 1) * sizeof **lvlist);
	if (lvlist[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(lvlist[j]);
	    }
	    free(lvlist);
	    lvlist = NULL;
	}
	lvlist[i][0] = order;
    }

    return lvlist;
}

static int var_lists_init (struct var_lists *vl,
			   int ndet, int nstoch, 
			   int order)
{
    int nreg = 1 + ndet + nstoch * order;
    int ntest = nstoch + ndet;

    if (order > 0) {
	/* test max lag for missing values */
	ntest += nstoch;
    }

#if VAR_DEBUG
    fprintf(stderr, "var_lists_init: order = %d, nreg = %d, ntest = %d\n", 
	    order, nreg, ntest);
#endif

    vl->detvars = NULL;
    vl->stochvars = NULL;
    vl->reglist = NULL;
    vl->testlist = NULL;
    vl->lagvlist = NULL;

    vl->detvars = gretl_list_new(ndet);
    vl->stochvars = gretl_list_new(nstoch);
    vl->reglist = gretl_list_new(nreg);
    vl->testlist = gretl_list_new(ntest);

    if (vl->detvars == NULL || vl->stochvars == NULL ||
	vl->reglist == NULL || vl->testlist == NULL) {
	goto bailout;
    }

    vl->lagvlist = lagvlist_construct(nstoch, order);
    if (vl->lagvlist == NULL) {
	goto bailout;
    }

    return 0;
    
 bailout:

    var_lists_free(vl);

    return E_ALLOC;
}

int var_max_order (const int *list, const DATAINFO *pdinfo)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int nstoch = 0, ndet = 0;
    int gotsep = 0;
    int order = 1;
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (!gotsep) {
	    nstoch++;
	} else {
	    ndet++;
	}
    }

    order = (T - ndet) / nstoch;

    while (order > 0) {
	int t1 = (order > pdinfo->t1)? order : pdinfo->t1;

	T = pdinfo->t2 - t1 + 1;
	if (nstoch * order + ndet > T) {
	    order--;
	} else {
	    break;
	}
    }

    return order - 1;
}

static int centered_dummy (const DATAINFO *pdinfo, int v)
{
    if (!strcmp(VARLABEL(pdinfo, v), "centered periodic dummy")) {
	return 1;
    } else {
	return 0;
    }
}

/* Given an incoming regression list, separate it into deterministic
   components (constant, trend, dummy variables) and stochastic
   components, and construct a list for each sort of variable.  Also
   allocate a VAR regression list that is long enough to hold the
   deterministic vars plus order lags of each stochastic var.
*/

static int organize_var_lists (const int *list, const double **Z,
			       const DATAINFO *pdinfo, int order,
			       struct var_lists *vlists)
{
    int ndet = 0, nstoch = 0;
    int gotsep = 0;
    char *d;
    int i, j, k, li;
    
    d = calloc(list[0] + 1, 1);
    if (d == NULL) {
	return E_ALLOC;
    }

    /* figure out the lengths of the lists */
    for (i=1; i<=list[0]; i++) {
	li = list[i];
	if (li == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (gotsep || 
	    !strcmp(pdinfo->varname[li], "const") ||	   
	    !strcmp(pdinfo->varname[li], "time") ||
	    centered_dummy(pdinfo, li) ||
	    gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[li])) {
	    d[i] = 1;
	    ndet++;
	} else {
	    nstoch++;
	}
    }

    /* check for degrees of freedom */
    if (nstoch * order + ndet > pdinfo->t2 - pdinfo->t1 + 1) {
	free(d);
	return E_DF;
    }

    /* allocate the lists */
    if (var_lists_init(vlists, ndet, nstoch, order)) {
	free(d);
	return E_ALLOC;
    }

    /* fill out the detvars and stochvars lists */
    j = k = 1;
    for (i=1; i<=list[0]; i++) {
	if (list[i] != LISTSEP) {
	    if (d[i]) {
		vlists->detvars[j++] = list[i];
	    } else {
		vlists->stochvars[k++] = list[i];
	    }
	}
    }

    free(d);

#if VAR_DEBUG
    printlist(vlists->detvars, "deterministic vars");
    printlist(vlists->stochvars, "stochastic vars");
#endif

    return 0;
}

/* Compose a VAR regression list: it may be complete, or one variable
   may be omitted (to run an F-test), or the order may be one less
   than the full VAR order (again, for an F-test).  If test = 1, also
   compose a list containing the maximum lag of each stochastic
   variable, for use in testing for missing values.
*/

static int
compose_varlist (struct var_lists *vl, int depvar, int order, int omit, 
		 int test, const DATAINFO *pdinfo)
{
    int l0 = 1 + vl->detvars[0] + order * vl->stochvars[0];
    int i, j, pos;
    int err = 0;

#if VAR_DEBUG
    fprintf(stderr, "compose_varlist: order = %d\n", order);
#endif

    if (omit) {
	l0 -= order;
    } 

    vl->reglist[0] = l0;
    vl->reglist[1] = depvar;

    pos = 2;
    for (i=1; i<=vl->stochvars[0]; i++) {
	if (i != omit) {
	    /* insert order lags of the given var */
	    for (j=1; j<=order; j++) {
		vl->reglist[pos++] = vl->lagvlist[i-1][j];
	    }
	}
    }

    /* append the deterministic vars */
    for (i=1; i<=vl->detvars[0]; i++) {
	vl->reglist[pos++] = vl->detvars[i];
    }

#if VAR_DEBUG
    printlist(vl->reglist, "composed VAR list");
#endif

    if (test) {
	/* now build the test list (to screen missing values) */
	pos = 1;
	for (i=1; i<=vl->stochvars[0]; i++) {
	    vl->testlist[pos++] = vl->stochvars[i];
	    if (order > 0) {
		/* include max lag */
		vl->testlist[pos++] = vl->lagvlist[i-1][order];
	    }
	}

	for (i=1; i<=vl->detvars[0]; i++) {
	    vl->testlist[pos++] = vl->detvars[i];
	}

#if VAR_DEBUG
	printlist(vl->testlist, "composed test list");
#endif
    }

    return err;
}

static int add_model_data_to_var (GRETL_VAR *var, const MODEL *pmod, int k)
{
    int i, j;
    int v = 0, lag = 0;
    int start = pmod->ifc;
    int rowmax = var->neqns * var->order + start;
    int err = 0;

    if (k == 0) {
	/* first equation: set up storage for residuals */
	var->T = pmod->t2 - pmod->t1 + 1;
	var->E = gretl_matrix_alloc(var->T, var->neqns);
	if (var->E == NULL) {
	    err = 1;
	} else {
	    var->ifc = pmod->ifc;
	}
    }

    /* save residuals */
    if (!err) {
	for (i=0; i<var->T; i++) {
	    gretl_matrix_set(var->E, i, k, pmod->uhat[pmod->t1 + i]);
	}
    }	

    /* save coefficients */
    if (!err) {
	for (i=start; i<rowmax; i++) {
	    if ((i - start) % var->order == 0) {
		v++;
		lag = 1;
	    } else {
		lag++;
	    }
	    j = (lag - 1) * var->neqns + v - 1;
	    gretl_matrix_set(var->A, k, j, pmod->coeff[i]);
	}
    }

    return err;
}

static int gretl_VAR_add_roots (GRETL_VAR *var)
{
    gretl_matrix *CompForm = NULL;
    double *eigA = NULL;
    double x, y;
    int i, np, err = 0;

    if (var->A == NULL) {
	return 1;
    }

    np = gretl_matrix_rows(var->A);

    var->lambda = gretl_matrix_alloc(np, 2);
    if (var->lambda == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
	CompForm = gretl_matrix_copy(var->A);
	if (CompForm == NULL) {
	    err = E_ALLOC;
	}
    }

    /* save eigenvalues of companion form matrix in polar form */
    if (!err) {
        eigA = gretl_general_matrix_eigenvals(CompForm, NULL);
	if (eigA == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<np; i++) {
		x = eigA[i];
		y = eigA[np+i];
		gretl_matrix_set(var->lambda, i, 0, atan2(y, x));
		gretl_matrix_set(var->lambda, i, 1, sqrt(x * x + y * y));
	    }
#if 0
	    gretl_matrix_print(var->A, "Companion form matrix", NULL);
	    gretl_matrix_print(var->lambda, "Eigenvalues in polar form", NULL);
#endif
	}
    }

    free(eigA);
    gretl_matrix_free(CompForm);

    if (err) {
	gretl_matrix_free(var->lambda);
	var->lambda = NULL;
    }

    return err;
}

const gretl_matrix *gretl_VAR_get_roots (GRETL_VAR *var)
{
    if (var->lambda == NULL) {
	/* roots not computed yet */
	gretl_VAR_add_roots(var);
    }

    return var->lambda;
}

static int VAR_LR_lag_test (GRETL_VAR *var)
{
    gretl_matrix *S = NULL;
    double ldet = NADBL;
    int err = 0;

    S = gretl_matrix_alloc(var->neqns, var->neqns);
    if (S == NULL) {
	err = 1;
    }

    if (!err) {
	gretl_matrix_multiply_mod(var->F, GRETL_MOD_TRANSPOSE,
				  var->F, GRETL_MOD_NONE,
				  S);
	gretl_matrix_divide_by_scalar(S, var->T);
    }

    if (!err) {
	ldet = gretl_vcv_log_determinant(S);
	if (na(var->ldet)) {
	    err = 1;
	}
    }    

    if (!err) {
	var->LR = var->T * (ldet - var->ldet);
    }

    gretl_matrix_free(S);

    /* we're done with this set of residuals */
    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

/* per-equation F-tests for excluding variables and maximum
   lag */

static int VAR_compute_tests (MODEL *varmod, GRETL_VAR *var,
			      struct var_lists *vl,
			      double ***pZ, DATAINFO *pdinfo,
			      int i, int *k) 
{
    MODEL testmod;
    double F = NADBL;
    int robust = gretl_model_get_int(varmod, "robust");
    int depvar = vl->stochvars[i + 1];
    int *outlist = NULL;
    int j, err = 0;

    if (i == 0 && var->order > 1) {
	/* first equation: allocate residual matrix for likelihood
	   ratio test on maximum lag */
	var->F = gretl_matrix_alloc(var->T, var->neqns);
    } 

    if (robust) {
	outlist = malloc(varmod->list[0] * sizeof *outlist);
	if (outlist == NULL) {
	    return E_ALLOC;
	}
    }

    /* restrictions for all lags of specific variables */
    for (j=0; j<var->neqns && !err; j++) {

	compose_varlist(vl, depvar, var->order, j + 1, 0, pdinfo);	

	if (robust) {
	    gretl_list_diff(outlist, varmod->list, vl->reglist);
	    F = robust_omit_F(outlist, varmod);
	    if (na(F)) {
		err = 1;
	    }
	} else {
	    testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A, 0.0);
	    err = testmod.errcode;
	    if (!err) {
		F = ((testmod.ess - varmod->ess) / var->order) / 
		    (varmod->ess / varmod->dfd);
	    }
	    clear_model(&testmod);
	}

	if (!err) {
	    var->Fvals[*k] = F;
	    *k += 1;
	}
    }
    
    /* restrictions for last lag, all variables */
    if (!err && var->order > 1) {

	compose_varlist(vl, depvar, var->order - 1, 0, 0, pdinfo);	

	testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A, 0.0);
	err = testmod.errcode;

	if (!err) {
	    if (robust) {
		gretl_list_diff(outlist, varmod->list, vl->reglist);
		F = robust_omit_F(outlist, varmod);
		if (na(F)) {
		    err = 1;
		}
	    } else {
		F = ((testmod.ess - varmod->ess) / var->neqns) / 
		    (varmod->ess / varmod->dfd);
	    }
	}

	if (!err && var->F != NULL) {
	    /* record residuals for LR test */
	    int j, t = testmod.t1;

	    for (j=0; j<var->T; j++) {
		gretl_matrix_set(var->F, j, i, testmod.uhat[t++]);
	    }
	}

	clear_model(&testmod);

	if (!err) {
	    var->Fvals[*k] = F;
	    *k += 1;
	}
    }

    if (outlist != NULL) {
	free(outlist);
    }

    return err;
}

static int VAR_add_stats (GRETL_VAR *var)
{
    int err = 0;

    var->S = gretl_matrix_alloc(var->neqns, var->neqns);
    if (var->S == NULL) {
	err = 1;
    }

    if (!err) {
	gretl_matrix_multiply_mod(var->E, GRETL_MOD_TRANSPOSE,
				  var->E, GRETL_MOD_NONE,
				  var->S);
	gretl_matrix_divide_by_scalar(var->S, var->T);
    }

    if (!err) {
	var->ldet = gretl_vcv_log_determinant(var->S);
	if (na(var->ldet)) {
	    err = 1;
	}
    }    

    if (!err) {
	int T = var->T;
	int g = var->neqns;
	int k = var->ncoeff;

	var->ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * var->ldet;
	var->AIC = (-2.0 * var->ll + 2.0 * k * g) / T;
	var->BIC = (-2.0 * var->ll + log(T) * k * g) / T;
    }

    if (!err && var->F != NULL) {
	VAR_LR_lag_test(var);
    }

    return err;
}


/* construct the respective VAR lists by adding the appropriate
   number of lags ("order") to the variables in list 

   Say the list is "x_1 const time x_2 x_3", and the order is 2.
   Then the first list should be

   x_1 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

   the second:

   x_2 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

   and so on.

   Run the regressions and print the results.
*/

static GRETL_VAR *real_var (int order, const int *inlist, 
			    double ***pZ, DATAINFO *pdinfo,
			    gretlopt opt, int *err)
{
    GRETL_VAR *var = NULL;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    struct var_lists vlists;
    gretlopt lsqopt = OPT_A | OPT_Z;
    int i, k, neqns;

    if (order < 1) {
	fprintf(stderr, I_("Not much point in a zero-order \"VAR\" surely?\n"));
	*err = 1;
	return NULL;
    }

    if (opt & OPT_R) {
	lsqopt |= OPT_R;
    }

    *err = organize_var_lists(inlist, (const double **) *pZ, pdinfo, 
			      order, &vlists);
    if (*err) {
	return NULL;
    }

    /* generate the required lags */
    if (real_list_laggenr(vlists.stochvars, pZ, pdinfo, 
			  order, vlists.lagvlist)) {
	*err = E_ALLOC;
	goto var_bailout;
    }

    neqns = vlists.stochvars[0];    

    /* compose base VAR list (entry 1 will vary across equations);
       assemble test list for t1 and t2 while we're at it */
    *err = compose_varlist(&vlists, vlists.stochvars[1], 
			   order, 0, 1, pdinfo);
    if (*err) {
	*err = E_DATA;
	goto var_bailout;
    }

    /* sort out sample range */
    if (check_for_missing_obs(vlists.testlist, &pdinfo->t1, &pdinfo->t2,
			      (const double **) *pZ, NULL)) {
	*err = E_MISSDATA;
	goto var_bailout;
    }
    
    var = gretl_VAR_new(neqns, order, pdinfo);
    if (var == NULL) {
	*err = E_ALLOC;
	goto var_bailout;
    }

    k = 0;

    for (i=0; i<neqns && !*err; i++) {
	MODEL *pmod = var->models[i];

	compose_varlist(&vlists, vlists.stochvars[i + 1], 
			order, 0, 0, pdinfo);

	*pmod = lsq(vlists.reglist, pZ, pdinfo, VAR, lsqopt, 0.0);

	if (pmod->errcode) {
	    *err = pmod->errcode;
	} else {
	    pmod->aux = AUX_VAR;
	    pmod->ID = i + 1;
	}

	if (!*err) {
	    *err = add_model_data_to_var(var, pmod, i);
	}

	if (!*err) {
	    *err = VAR_compute_tests(pmod, var, &vlists, pZ, pdinfo, i, &k);
	}
    }

 var_bailout:

    var_lists_free(&vlists);

    /* reset sample range */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    if (!*err) {
	var->ncoeff = var->models[0]->ncoeff;
	var->t1 = var->models[0]->t1;
	var->t2 = var->models[0]->t2;
	var->T = var->t2 - var->t1 + 1;
	*err = VAR_add_stats(var);
    }

    if (!*err) {
	*err = gretl_VAR_do_error_decomp(var->neqns, var->S, var->C);
    }

    if (*err) {
	gretl_VAR_free(var);
	var = NULL;
    }

#if 0
    if (!*err) {
	gretl_matrix_print(var->A, "var->A", NULL);
    }
#endif

    return var;
}

static int *
maybe_expand_VAR_list (const int *list, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, int *err)
{
    int needsep, addconst = 0, addseas = 0;
    int *vlist = NULL;
    int i, l0, di0 = 0;

    if (!(opt & OPT_N)) {
	addconst = 1;
    }

    if (pdinfo->pd > 1 && (opt & OPT_D)) {
	addseas = 1;
    }

    if (!addconst && !addseas) {
	return NULL;
    }

    if (addseas) {
	di0 = dummy(pZ, pdinfo, 0);
	if (di0 == 0) {
	    *err = E_ALLOC;
	    return NULL;
	}
    } 

    needsep = !gretl_list_has_separator(list);
    l0 = list[0] + needsep;

    if (addconst) {
	l0 += 1;
    }
    if (addseas) {
	l0 += pdinfo->pd - 1;
    }

    vlist = gretl_list_new(l0);

    if (vlist == NULL) {
	*err = E_ALLOC;
    } else {
	int j = 1;

	for (i=1; i<=list[0]; i++) {
	    vlist[j++] = list[i];
	}
	if (needsep) {
	    vlist[j++] = LISTSEP;
	}    
	if (addseas) {
	    for (i=0; i<pdinfo->pd - 1; i++) {
		vlist[j++] = di0 + i;
	    }
	}
	if (addconst) {
	    vlist[j++] = 0;
	}
    }

    return vlist;
}

/**
 * simple_VAR:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: if includes %OPT_R, use robust VCV;
 *       if includes %OPT_I, print impulse responses;
 *       if includes %OPT_F, print forecast variance decompositions;
 *       if includes %OPT_D, add seasonal dummies;
 *       if includes %OPT_N, do not include a constant.
 * @prn: gretl printing struct.
 *
 * Estimate a vector autoregression (VAR) and print the results.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int simple_VAR (int order, int *list, double ***pZ, DATAINFO *pdinfo,
		gretlopt opt, PRN *prn)
{
    GRETL_VAR *var = NULL;
    int *vlist = NULL;
    int err = 0;
    
    gretl_list_purge_const(list);
    vlist = maybe_expand_VAR_list(list, pZ, pdinfo, opt, &err);

    if (!err) {
	var = real_var(order, (vlist != NULL)? vlist : list, 
		       pZ, pdinfo, opt, &err);
    }

    if (var != NULL) {
	gretl_VAR_print(var, pdinfo, opt, prn);
	gretl_VAR_free(var);
    }

    free(vlist);

    return err;
}

/**
 * full_VAR:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: if includes %OPT_R, use robust VCV;
 *       if includes %OPT_I, print impulse responses;
 *       if includes %OPT_F, print forecast variance decompositions;
 *       if includes %OPT_D, add seasonal dummies;
 *       if includes %OPT_N, do not include a constant.
 * @prn: gretl printing struct.
 *
 * Estimate a vector auto-regression (VAR), print and save
 * the results.
 *
 * Returns: pointer to VAR struct, which may be %NULL on error.
 */

GRETL_VAR *full_VAR (int order, int *list, double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn)
{
    GRETL_VAR *var = NULL;
    int *vlist = NULL;
    int err = 0;

    gretl_list_purge_const(list);
    vlist = maybe_expand_VAR_list(list, pZ, pdinfo, opt, &err);

    if (!err) {
	var = real_var(order, (vlist != NULL)? vlist : list, 
		       pZ, pdinfo, opt, &err);
    }

    if (var != NULL) {
	gretl_VAR_print(var, pdinfo, opt, prn);
    }

    free(vlist);

    return var;
}

static double df_pvalue_from_plugin (double tau, int n, int niv, int itv)
{
    char datapath[FILENAME_MAX];
    void *handle;
    double (*mackinnon_pvalue)(double, int, int, int, char *);
    double pval = NADBL;
    static int nodata;
    
    if (nodata) {
	return pval;
    }

    mackinnon_pvalue = get_plugin_function("mackinnon_pvalue", &handle);
    if (mackinnon_pvalue == NULL) {
	nodata = 1;
        return pval;
    }

    strcpy(datapath, gretl_lib_path());
#ifdef WIN32
    append_dir(datapath, "plugins");
#endif

    pval = (*mackinnon_pvalue)(tau, n, niv, itv, datapath);

#if 0
    fprintf(stderr, "getting pval: tau=%g, n=%d, niv=%d, itv=%d: pval=%g\n",
	    tau, n, niv, itv, pval);
#endif

    close_plugin(handle);

    if (*datapath == '\0') {
	nodata = 1;
    } 

    return pval;
}

/**
 * coint:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: if OPT_N, do not an include a constant in the
 *       cointegrating regression.
 * @prn: gretl printing struct.
 *
 * Test for cointegration.  
 *
 * Returns: 0 on successful completion.
 */

int coint (int order, const int *list, double ***pZ, 
	   DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int i, t, n, nv, l0 = list[0];
    int hasconst = gretl_list_has_const(list);
    MODEL cmod;
    int *cointlist = NULL;

    if (order <= 0 || list[0] - hasconst < 2) {
	strcpy(gretl_errmsg, "coint: needs a positive lag order "
	       "and at least two variables");
	return 1;
    }

    gretl_model_init(&cmod);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	if (list[i] == 0) {
	    continue;
	}
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		i, pdinfo->varname[list[i]]);
	real_adf_test(list[i], order, 1, pZ, pdinfo, OPT_NONE, 
		      ADF_EG_TEST, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (!hasconst && !(opt & OPT_N)) {
	/* add const to coint regression list */
	cointlist = malloc((l0 + 2) * sizeof *cointlist);
	if (cointlist == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<=l0; i++) {
	    cointlist[i] = list[i];
	}
	cointlist[l0 + 1] = 0;
	cointlist[0] += 1;
    } else {
	cointlist = gretl_list_copy(list);
	if (cointlist == NULL) {
	    return E_ALLOC;
	}
    }

    pprintf(prn, _("Step %d: cointegrating regression\n"), l0 + 1);
    
    cmod = lsq(cointlist, pZ, pdinfo, OLS, OPT_NONE, 0.0); 
    cmod.aux = AUX_COINT;
    printmodel(&cmod, pdinfo, OPT_NONE, prn);

    /* add residuals from cointegrating regression to data set */
    n = pdinfo->n;
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }
    nv = pdinfo->v - 1;

    for (t=0; t<cmod.t1; t++) {
	(*pZ)[nv][t] = NADBL;
    }
    for (t=cmod.t1; t<=cmod.t2; t++) {
	(*pZ)[nv][t] = cmod.uhat[t];
    }
    for (t=cmod.t2+1; t<n; t++) {
	(*pZ)[nv][t] = NADBL;
    }

    strcpy(pdinfo->varname[nv], "uhat");

    pputc(prn, '\n');
    pprintf(prn, _("Step %d: Dickey-Fuller test on residuals\n"), l0 + 2);

    /* Run (A)DF test on the residuals */
    real_adf_test(pdinfo->v - 1, order, 1 + cmod.ncoeff - cmod.ifc, 
		  pZ, pdinfo, OPT_N, ADF_EG_TEST | ADF_PRINT_ACK, prn);

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		 "(a) The unit-root hypothesis is not rejected for the individual"
		 " variables.\n(b) The unit-root hypothesis is rejected for the "
		 "residuals (uhat) from the \n    cointegrating regression.\n"));

    /* clean up and get out */
    clear_model(&cmod);
    free(cointlist);
    dataset_drop_last_variables(1, pZ, pdinfo);

    return 0;
}

static int *adf_prepare_vars (int order, int varno,
			      double ***pZ, DATAINFO *pdinfo)
{
    int i, orig_t1 = pdinfo->t1;
    int *list;
    int err = 0;

    if (varno == 0) {
	return NULL;
    }

    list = malloc((6 + order) * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    /* temporararily reset sample */
    pdinfo->t1 = 0;

    /* generate first difference of the given variable */
    list[1] = diffgenr(varno, pZ, pdinfo, 0);
    if (list[1] < 0) {
	pdinfo->t1 = orig_t1;
	free(list);
	return NULL;
    }	

    /* generate lag of given var */
    list[2] = laggenr(varno, 1, pZ, pdinfo); 
    if (list[2] < 0) {
	pdinfo->t1 = orig_t1;
	free(list);
	return NULL;
    }

    /* undo reset sample */
    pdinfo->t1 = orig_t1;

    /* generate lags of difference for augmented test */
    for (i=1; i<=order && !err; i++) {
	int lnum = laggenr(list[1], i, pZ, pdinfo);

	if (lnum < 0) {
	    fprintf(stderr, "Error generating lag variable\n");
	    err = 1;
	} else {
	    list[2 + i] = lnum;
	} 
    } 

    return list;
}

#define ADF_DEBUG 0

static int auto_adjust_order (int *list, int order_max,
			      double ***pZ, DATAINFO *pdinfo,
			      PRN *prn)
{
    MODEL kmod;
    double tstat, pval = 1.0;
    int i, k = order_max;

    for (k=order_max; k>0; k--) {
	int j = k;

	if (list[list[0]] == 0) j++;

	kmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);

	if (kmod.errcode) {
	    clear_model(&kmod);
	    fprintf(stderr, "adf: model failed in auto_adjust_order()\n");
	    k = -1;
	    break;
	}

#if ADF_DEBUG
	printmodel(&kmod, pdinfo, OPT_NONE, prn);
#endif

	tstat = kmod.coeff[j] / kmod.sderr[j];
	clear_model(&kmod);
	pval = normal_pvalue_2(tstat);

	if (pval > 0.10) {
#if ADF_DEBUG
	    fprintf(stderr, "auto_adjust_order: lagged difference not "
		    "significant at order %d (t = %g)\n", k, tstat);
#endif
	    if (k == 1) {
		k = 0;
		break;
	    } else {
		for (i=k+2; i<list[0]; i++) {
		    list[i] = list[i+1];
		}
		list[0] -= 1;
	    }
	} else {
#if ADF_DEBUG
	    fprintf(stderr, "auto_adjust_order: lagged difference is "
		    "significant at order %d (t = %g)\n", k, tstat);
#endif
	    break;
	}
    }

    return k;
}

static void copy_list_values (int *targ, const int *src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }
}

static void 
print_adf_results (int order, double DFt, double pv, const MODEL *dfmod,
		   int dfnum, const char *vname, int *blurb_done,
		   unsigned char flags, int i, PRN *prn)
{
    const char *models[] = {
	"(1 - L)y = (a-1)*y(-1) + e",
	"(1 - L)y = b0 + (a-1)*y(-1) + e",
	"(1 - L)y = b0 + b1*t + (a-1)*y(-1) + e",
	"(1 - L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + e"
    };
    const char *aug_models[] = {
	"(1 - L)y = (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + b1*t + (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + ... + e"
    };
    const char *teststrs[] = {
	N_("test without constant"),
	N_("test with constant"),
	N_("with constant and trend"),
	N_("with constant and quadratic trend")
    };

    char pvstr[48];

    if (prn == NULL) return;

    if (na(pv)) {
	sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
    } else {
	sprintf(pvstr, "%s %.4g", 
		(order > 0)? _("asymptotic p-value") : _("p-value"), 
		pv);
    } 

    if (*blurb_done == 0) {
	if (order > 0) {
	    pprintf(prn, _("\nAugmented Dickey-Fuller tests, order %d, for %s\n"),
		    order, vname);
	} else {
	    pprintf(prn, _("\nDickey-Fuller tests for %s\n"), vname);
	}
	pprintf(prn, _("sample size %d\n"), dfmod->nobs);
	pputs(prn, _("unit-root null hypothesis: a = 1"));
	pputs(prn, "\n\n");
	*blurb_done = 1;
    }

    pprintf(prn, "   %s\n", _(teststrs[i]));

    if (!(flags & ADF_EG_TEST)) {
	pprintf(prn, "   %s: %s\n", _("model"), 
		(order > 0)? aug_models[i] : models[i]);
    }

    pprintf(prn, "   %s: %g\n"
	    "   %s: t = %g\n"
	    "   %s\n",
	    _("estimated value of (a - 1)"), dfmod->coeff[dfnum],
	    _("test statistic"), DFt,
	    pvstr);	
}

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags,
			  PRN *prn)
{
    MODEL dfmod;

    int orig_nvars = pdinfo->v;
    int blurb_done = 0;
    int auto_order = 0;
    int order_max = 0;
    int *list;
    int *biglist = NULL;
    double DFt = NADBL;
    double pv = NADBL;
    char mask[4] = {0};
    int i, itv;
    int err = 0;

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: got order = %d\n", order);
#endif

    if (order < 0) {
	auto_order = 1;
	order = -order;
    }

    order_max = order;

    list = adf_prepare_vars(order, varno, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    if (auto_order) {
	int tmp = list[0];

	list[0] = order + 5;
	biglist = gretl_list_copy(list);
	if (biglist == NULL) {
	    free(list);
	    return E_ALLOC;
	}
	list[0] = tmp;
    }

    gretl_model_init(&dfmod);

    if (opt == OPT_NONE || opt == OPT_V) {
	/* default display */
	mask[1] = mask[2] = mask[3] = 1;
    } else {
	if (opt & OPT_N) {
	    /* nc model */
	    mask[0] = 1;
	}
	if (opt & OPT_C) {
	    /* c */
	    mask[1] = 1;
	}
	if (opt & OPT_T) {
	    /* ct */
	    mask[2] = 1;
	}
	if (opt & OPT_R) {
	    /* ctt */
	    mask[3] = 1;
	}
    }

    for (i=0; i<4; i++) {
	int dfnum = (i > 0);

	if (mask[i] == 0) {
	    continue;
	}

	if (auto_order) {
	    order = order_max;
	    copy_list_values(list, biglist);
	}

	list[0] = 2 + order + i;

	if (i > 0) {
	    list[list[0]] = 0;
	} 

	if (i >= 2) {
	    list[3 + order] = gettrend(pZ, pdinfo, 0);
	    if (list[3 + order] == TREND_FAILED) {
		err = E_ALLOC;
		goto bailout;
	    }
	}

	if (i > 2) {
	    list[4 + order] = gettrend(pZ, pdinfo, 1);
	    if (list[4 + order] == TREND_FAILED) {
		err = E_ALLOC;
		goto bailout;
	    }
	}

	if (auto_order) {
	    order = auto_adjust_order(list, order_max, pZ, pdinfo, prn);
	    if (order < 0) {
		err = 1;
		clear_model(&dfmod);
		goto bailout;
	    }
	}

	dfmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (dfmod.errcode) {
	    fprintf(stderr, "adf_test: dfmod.errcode = %d\n", 
		    dfmod.errcode);
	    err = dfmod.errcode;
	    clear_model(&dfmod);
	    goto bailout;
	}

	DFt = dfmod.coeff[dfnum] / dfmod.sderr[dfnum];

	itv = (i == 0)? UR_NO_CONST :
	    (i == 1)? UR_CONST : 
	    (i == 2)? UR_TREND :
	    UR_TREND_SQUARED;

	pv = df_pvalue_from_plugin(DFt, 
				   /* use asymptotic p-value for augmented case */
				   (order > 0)? 0 : dfmod.nobs, 
				   niv, itv);

	if (!(opt & OPT_Q)) {
	    print_adf_results(order, DFt, pv, &dfmod, dfnum, pdinfo->varname[varno],
			      &blurb_done, flags, i, prn);
	}

	if (opt & OPT_V) {
	    /* verbose */
	    dfmod.aux = (order > 0)? AUX_ADF : AUX_DF;
	    if (!na(pv)) {
		gretl_model_set_int(&dfmod, "dfnum", dfnum + 2);
		gretl_model_set_double(&dfmod, "dfpval", pv);
	    }
	    printmodel(&dfmod, pdinfo, OPT_NONE, prn);
	} else if (!(opt & OPT_Q)) {
	    pputc(prn, '\n');
	}

	clear_model(&dfmod);
    }

    if (!err) {
	if (!(flags & ADF_EG_TEST)) {
	    record_test_result(DFt, pv, "Dickey-Fuller");
	}
	if ((flags & ADF_PRINT_ACK) && !(opt & OPT_Q)) {
	    pputs(prn, _("P-values based on MacKinnon (JAE, 1996)\n"));
	}	
    }

 bailout:

    free(list);

    if (biglist != NULL) {
	free(biglist);
    }

    dataset_drop_last_variables(pdinfo->v - orig_nvars, pZ, pdinfo);

    return err;
}

/**
 * adf_test:
 * @order: lag order for the test.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller test for 
 * a unit root.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int adf_test (int order, int varno, double ***pZ,
	      DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    return real_adf_test(varno, order, 1, pZ, pdinfo, opt, 
			 ADF_PRINT_ACK, prn);
}

/**
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the KPSS test for 
 * stationarity.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int kpss_test (int order, int varno, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    MODEL KPSSmod;
    int list[4];
    int hastrend = 0;
    double s2 = 0.0;
    double cumsum = 0.0, cumsum2 = 0.0;
    double teststat;
    double *autocov;
    double et;

    int i, t;
    int t1, t2, T;

    /* sanity check */
    if (order < 0 || varno <= 0 || varno >= pdinfo->v) {
	return 1;
    }

    if (opt & OPT_T) {
	hastrend = 1;
    }

    list[0] = (2 + hastrend);
    list[1] = varno;
    list[2] = 0;
    if (hastrend) {
	list[3] = gettrend(pZ, pdinfo, 0);
    }

    /* OPT_M: reject missing values within sample range */
    KPSSmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_M, 0.0);
    if (KPSSmod.errcode) {
	clear_model(&KPSSmod);
	return KPSSmod.errcode;
    }

    t1 = KPSSmod.t1;
    t2 = KPSSmod.t2;
    T = KPSSmod.nobs;

    if (opt & OPT_V) {
	KPSSmod.aux = AUX_KPSS;
	printmodel(&KPSSmod, pdinfo, OPT_NONE, prn);
    }
  
    autocov = malloc(order * sizeof *autocov);
    if (autocov == NULL) {
	return E_ALLOC;
    }
  
    for (i=0; i<order; i++) {
	autocov[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	et = KPSSmod.uhat[t];
	if (na(et)) {
	    continue;
	}
	cumsum += et;
	cumsum2 += cumsum * cumsum;
	s2 += et * et;
	for (i=0; i<order; i++) {
	    int s = i + 1;

	    if (t - s >= t1) {
		autocov[i] += et * KPSSmod.uhat[t - s];
	    }
	}
#ifdef KPSS_DEBUG
	fprintf(stderr, "%d: %#12.4g %#12.4g %#12.4g %#12.4g \n", 
		t, et, KPSSmod.uhat[t-1], s2, cumsum2);
#endif
    }

    for (i=0; i<order; i++) {
	double wt = 1.0 - ((double) (i + 1)) / (order + 1);

	s2 += 2.0 * wt * autocov[i];
    }

    s2 /= T;
    teststat = cumsum2 / (s2 * T * T);

    if (opt & OPT_V) {
	pprintf(prn, "  %s: %g\n", _("Robust estimate of variance"), s2);
	pprintf(prn, "  %s: %g\n", _("Sum of squares of cumulated residuals"), 
		cumsum2);
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, _("\nKPSS test for %s %s\n\n"), pdinfo->varname[varno],
		(hastrend)? _("(including trend)") : _("(without trend)"));
	pprintf(prn, _("Lag truncation parameter = %d\n"), order);
	pprintf(prn, "%s = %g\n\n", _("Test statistic"), teststat);
	pprintf(prn, "		    10%%\t   5%%\t 2.5%%\t   1%%\n");
	if (hastrend) {
	    pprintf(prn, "%s: 0.119\t0.146\t0.176\t0.216\n\n", _("Critical values"));
	} else {
	    pprintf(prn, "%s: 0.347\t0.463\t0.574\t0.739\n\n", _("Critical values"));
	}
    }

    record_test_result(teststat, NADBL, "KPSS");
    clear_model(&KPSSmod);

    free(autocov);

    return 0;
}

static int allocate_johansen_sigmas (JohansenInfo *jv, int k)
{
    int vk = k;
    int err = 0;

    if (jv->code == J_REST_CONST || jv->code == J_REST_TREND) {
	vk++;
    }    

    jv->Suu = gretl_matrix_alloc(k, k);
    jv->Svv = gretl_matrix_alloc(vk, vk);
    jv->Suv = gretl_matrix_alloc(k, vk);

    if (jv->Suu == NULL || jv->Svv == NULL || jv->Suv == NULL) {
	gretl_matrix_free(jv->Suu);
	gretl_matrix_free(jv->Svv);
	gretl_matrix_free(jv->Suv);
	
	jv->Suu = NULL;
	jv->Svv = NULL;
	jv->Suv = NULL;

	err = E_ALLOC;
    } 

    return err;
}

static void 
print_johansen_sigmas (const JohansenInfo *jv, PRN *prn)
{
    int nr, nc;
    int i, j;

    pprintf(prn, "\n%s\n\n", _("Sample variance-covariance matrices for residuals"));

    nr = gretl_matrix_rows(jv->Suu);
    pprintf(prn, " %s\n\n", _("VAR system in first differences"));
    for (i=0; i<nr; i++) {
	for (j=0; j<nr; j++) {
	    pprintf(prn, "%#12.5g", gretl_matrix_get(jv->Suu, i, j));
	}
	pputc(prn, '\n');
    }

    nr = gretl_matrix_rows(jv->Svv);
    pprintf(prn, "\n %s\n\n", _("System with levels as dependent variable"));
    for (i=0; i<nr; i++) {
	for (j=0; j<nr; j++) {
	    pprintf(prn, "%#12.5g", gretl_matrix_get(jv->Svv, i, j));
	}
	pputc(prn, '\n');
    } 
    
    nr = gretl_matrix_rows(jv->Suv);
    nc = gretl_matrix_cols(jv->Suv);
    pprintf(prn, "\n %s\n\n", _("Cross-products"));
    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    pprintf(prn, "%#12.5g", gretl_matrix_get(jv->Suv, i, j));
	}
	pputc(prn, '\n');
    }     
}

static void
transcribe_uhat_to_matrix (const MODEL *pmod, gretl_matrix *u, int row)
{
    int j, cols = gretl_matrix_cols(u);
    int t = pmod->t1;

    for (j=0; j<cols; j++) {
	gretl_matrix_set(u, row, j, pmod->uhat[t++]);
    }
}

static void
transcribe_data_as_uhat (int v, const double **Z, gretl_matrix *u, 
			 int row, int t)
{
    int j, cols = gretl_matrix_cols(u);

    for (j=0; j<cols; j++) {
	gretl_matrix_set(u, row, j, Z[v][t++]);
    }
}

static int 
johansen_complete (GRETL_VAR *jvar, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    void *handle = NULL;
    int (*johansen) (GRETL_VAR *, double ***, DATAINFO *, PRN *);
    int err = 0;

    *gretl_errmsg = 0;
    
    johansen = get_plugin_function("johansen_analysis", &handle);

    if (johansen == NULL) {
	err = 1;
    } else {
	err = (* johansen) (jvar, pZ, pdinfo, prn);
	close_plugin(handle);
    }
    
    return err;
}

static int 
allocate_residual_matrices (struct var_resids *resids, int k,
			    JohansenCode jcode)
{
    int T = resids->t2 - resids->t1 + 1;
    int vk = k;
    int err = 0;

    if (jcode == J_REST_CONST || jcode == J_REST_TREND) {
	vk++;
    }

    resids->u = gretl_matrix_alloc(k, T);
    if (resids->u == NULL) {
	err = E_ALLOC;
    } else {
	resids->v = gretl_matrix_alloc(vk, T);
	if (resids->v == NULL) {
	    gretl_matrix_free(resids->u);
	    resids->u = NULL;
	    err = E_ALLOC;
	}
    }

    return err;
}

/* Create a "master list" for the models in a VECM: it contains all
   the required lags of the first differences of the endogenous vars,
   plus any deterministic vars, plus blank spaces for the dependent
   variable (which will vary across equations) and the Error
   Correction terms(s).  At the same time allocate a list that will
   contain the ID numbers of the first differences of the endogenous
   variables.
*/

static int make_johansen_VECM_lists (JohansenInfo *jv, int *varlist,
				     int neqns)
{
    int i, k = varlist[0] + jv->rank;
    int err = 0;

    jv->difflist = gretl_list_new(neqns);
    if (jv->difflist == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	jv->biglist = gretl_list_new(k);
	if (jv->biglist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=2; i<=varlist[0]; i++) {
		jv->biglist[i] = varlist[i];
	    }
	}
    }

#if VAR_DEBUG
    printlist(varlist, "make_johansen_VECM_lists: varlist");
    printlist(jv->biglist, "make_johansen_VECM_lists: jv->biglist");
#endif

    return err;
}

/* For Johansen analysis: estimate VAR in differences along with the
   other auxiliary regressions required to compute the relevant
   matrices of residuals, for concentration of the log-likelihood
*/

static int johansen_VAR (int order, const int *inlist, 
			 double ***pZ, DATAINFO *pdinfo,
			 GRETL_VAR *jvar, struct var_resids *resids, 
			 gretlopt opt, PRN *prn)
{
    struct var_lists vlists;
    MODEL jmod;
    int i, err = 0;

    err = organize_var_lists(inlist, (const double **) *pZ, pdinfo, 
			     order, &vlists);
    if (err) {
	return err;
    }

    /* generate the required lags, if any */
    if (order > 0) {
	if (real_list_laggenr(vlists.stochvars, pZ, pdinfo, 
			      order, vlists.lagvlist)) {
	    err = E_ALLOC;
	    goto var_bailout;
	}
    }

    jvar->neqns = vlists.stochvars[0];    

    /* compose base VAR list (entry 1 will vary across equations);
       assemble test list for setting t1 and t2 while we're at it 
    */
    err = compose_varlist(&vlists, vlists.stochvars[1], 
			  order, 0, 1, pdinfo);
    if (err) {
	err = E_DATA;
	goto var_bailout;
    }

    if (check_for_missing_obs(vlists.testlist, &pdinfo->t1, &pdinfo->t2,
			      (const double **) *pZ, NULL)) {
	err = E_MISSDATA;
	goto var_bailout;
    }

    resids->t1 = pdinfo->t1;
    resids->t2 = pdinfo->t2;
    err = allocate_residual_matrices(resids, jvar->neqns, jcode(jvar));
    if (err) {
	goto var_bailout;
    }

    if (jrank(jvar) > 0) {
	err = make_johansen_VECM_lists(jvar->jinfo, vlists.reglist, 
				       jvar->neqns);
	if (err) {
	    goto var_bailout;
	}	    
    }

    gretl_model_init(&jmod);

    if (opt & OPT_V) {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), order);
    }

    for (i=0; i<jvar->neqns; i++) {

	/* insert the appropriate dependent variable */
	vlists.reglist[1] = vlists.stochvars[i + 1];

	if (jrank(jvar) > 0) {
	    /* VECM: record ID number of first difference.  Note:
	       we'll need this information in johansen.c, for building
	       the final VECM models, even if the order of the present
	       VAR is zero.
	    */
	    jvar->jinfo->difflist[i+1] = vlists.reglist[1];
	}

	/* VAR in differences */
	if (vlists.reglist[0] == 1) {
	    /* degenerate model (nothing to concentrate out) */
	    transcribe_data_as_uhat(vlists.reglist[1], (const double **) *pZ,
				    resids->u, i, resids->t1);
	} else {
	    jmod = lsq(vlists.reglist, pZ, pdinfo, VAR, OPT_A | OPT_Z, 0.0);
	    if ((err = jmod.errcode)) {
		goto var_bailout;
	    }
	    if (opt & OPT_V) {
		jmod.aux = AUX_VAR;
		jmod.ID = i + 1;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    transcribe_uhat_to_matrix(&jmod, resids->u, i);
	    clear_model(&jmod);
	}

	/* y_{t-1} regressions */
	vlists.reglist[1] = resids->levels_list[i + 1];
	if (vlists.reglist[0] == 1) {
	    /* degenerate */
	    transcribe_data_as_uhat(vlists.reglist[1], (const double **) *pZ,
				    resids->v, i, resids->t1);
	} else {
	    jmod = lsq(vlists.reglist, pZ, pdinfo, VAR, OPT_A | OPT_Z, 0.0);
	    if ((err = jmod.errcode)) {
		goto var_bailout;
	    }	
	    if (opt & OPT_V) {
		jmod.aux = AUX_JOHANSEN;
		jmod.ID = -1;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    transcribe_uhat_to_matrix(&jmod, resids->v, i);
	    clear_model(&jmod);
	}
    }

    /* supplementary regressions for restricted cases */
    if (restricted(jvar)) {
	if (jcode(jvar) == J_REST_CONST) {
	    vlists.reglist[1] = 0;
	} else {
	    vlists.reglist[1] = gettrend(pZ, pdinfo, 0);
	}
	if (vlists.reglist[0] == 1) {
	    /* degenerate case */
	    transcribe_data_as_uhat(vlists.reglist[1], (const double **) *pZ,
				    resids->v, i, resids->t1);
	} else {	    
	    jmod = lsq(vlists.reglist, pZ, pdinfo, VAR, OPT_A, 0.0);
	    if ((err = jmod.errcode)) {
		goto var_bailout;
	    }
	    if (opt & OPT_V) {
		jmod.aux = AUX_JOHANSEN;
		jmod.ID = -1;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    transcribe_uhat_to_matrix(&jmod, resids->v, i);
	    clear_model(&jmod);
	}
    }     

    pputc(prn, '\n');

 var_bailout:

    var_lists_free(&vlists);

    return err;
}

static JohansenCode jcode_from_opt (gretlopt opt)
{
    JohansenCode jc = J_UNREST_CONST;

    if (opt & OPT_N) {
	jc = J_NO_CONST;
    } else if (opt & OPT_T) {
	jc = J_UNREST_TREND;
    } else if (opt & OPT_R) {
	jc = J_REST_CONST;
    } else if (opt & OPT_A) {
	jc = J_REST_TREND;
    }

    return jc;
}

static JohansenInfo *
johansen_info_new (const int *list, int rank, gretlopt opt)
{
    JohansenInfo *jv = malloc(sizeof *jv);

    if (jv != NULL) {
	jv->code = jcode_from_opt(opt);
	jv->list = gretl_list_copy(list);
	if (jv->list == NULL) {
	    free(jv);
	    jv = NULL;
	} else {
	    jv->ID = 0;

	    jv->u = NULL;
	    jv->v = NULL;
	    jv->w = NULL;

	    jv->Suu = NULL;
	    jv->Svv = NULL;
	    jv->Suv = NULL;

	    jv->Beta = NULL;
	    jv->Alpha = NULL;
	    jv->Bse = NULL;

	    jv->difflist = NULL;
	    jv->biglist = NULL;

	    jv->rank = rank;

	    jv->seasonals = 0;
	    jv->nexo = 0;
	}
    }

    return jv;
}

/* VECM with restricted const or trend: split the "v" residuals
   matrix into two components */

static int split_v_resids (GRETL_VAR *jvar, gretl_matrix *v)
{
    double x;
    int T = jvar->T;
    int n = jvar->neqns;
    int i, j;

    jvar->jinfo->v = gretl_matrix_alloc(n, T);
    if (jvar->jinfo->v == NULL) {
	return E_ALLOC;
    } 

    jvar->jinfo->w = gretl_vector_alloc(T);
    if (jvar->jinfo->w == NULL) {
	gretl_matrix_free(jvar->jinfo->v);
	jvar->jinfo->v = NULL;
	return E_ALLOC;
    } 

    for (i=0; i<n; i++) {
	for (j=0; j<T; j++) {
	    x = gretl_matrix_get(v, i, j);
	    gretl_matrix_set(jvar->jinfo->v, i, j, x);
	}
    }

    for (j=0; j<T; j++) {
	x = gretl_matrix_get(v, n, j);
	gretl_vector_set(jvar->jinfo->w, j, x);
    }  

    gretl_matrix_free(v);
	
    return 0;
}

static GRETL_VAR *
johansen_VAR_new (const int *list, int rank, int order, gretlopt opt)
{
    GRETL_VAR *var = malloc(sizeof *var);

    if (var != NULL) {
	var->jinfo = johansen_info_new(list, rank, opt);
	if (var->jinfo == NULL) {
	    free(var);
	    var = NULL;
	}
    }

    if (var != NULL) {
	var->ci = VECM;
	var->err = 0;

	var->neqns = 0;
	var->order = order;
	var->ncoeff = 0;

	var->A = NULL;
	var->lambda = NULL;
	var->E = NULL;
	var->C = NULL;
	var->S = NULL;
	var->F = NULL;

	var->models = NULL;
	var->Fvals = NULL;
	var->name = NULL;

	var->ll = var->ldet = NADBL;
	var->AIC = var->BIC = NADBL;
	var->LR = NADBL;
    }

    return var;
}

static GRETL_VAR *johansen_driver (int order, int rank, const int *list, 
				   const int *exolist, double ***pZ, DATAINFO *pdinfo,
				   gretlopt opt, PRN *prn)
{
    PRN *varprn = NULL;
    GRETL_VAR *jvar;

    int seasonals = (opt & OPT_D);
    struct var_resids resids;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    int orig_v = pdinfo->v;
    int *varlist = NULL;
    int nexo, di0 = 0, l0 = list[0];
    int i, k;

    jvar = johansen_VAR_new(list, rank, order - 1, opt);
    if (jvar == NULL) {
	return NULL;
    }    

    var_resids_init(&resids);

    if (order <= 0 || list[0] < 2) {
	strcpy(gretl_errmsg, "coint2: needs a positive lag order "
	       "and at least two variables");
	jvar->err = 1;
	return jvar;
    }

    /* nexo will hold total number of deterministic/exogenous vars */
    nexo = 0;

    if (seasonals) {
	if (pdinfo->pd > 1) {
#if 0
	    /* Center the dummies only if there's a constant in the VAR? 
	       This seems to be a bad idea.
	    */
	    int flag = (jcode(jvar) >= J_UNREST_CONST)? 1 : -1;
#else
	    int flag = 1;
#endif
	    jvar->jinfo->seasonals = pdinfo->pd - 1;
	    nexo += pdinfo->pd - 1;
	    di0 = dummy(pZ, pdinfo, flag);
	    if (di0 == 0) {
		jvar->err = E_ALLOC;
	    }
	} else {
	    fprintf(stderr, "seasonals option ignored\n");
	}
    }

    if (jvar->err) {
	return jvar;
    }

    if (jcode(jvar) >= J_UNREST_CONST) {
	nexo++;
    }
    if (jcode(jvar) == J_UNREST_TREND) {
	nexo++;
    }
    if (exolist != NULL) {
	nexo += exolist[0];
    }
    if (nexo > 0) {
	l0 += nexo + 1; /* include list separator */
    } 

    /* "levels_list" will hold the first lags of the endogenous variables,
       which figure as the dependent variables in the second set of
       Johansen preliminary regressions */
    resids.levels_list = gretl_list_new(list[0]);
    if (resids.levels_list == NULL) {
	jvar->err = E_ALLOC;
	goto johansen_exit;
    }

    /* full VAR list, including both endog and exog vars */
    varlist = gretl_list_new(l0);
    if (varlist == NULL) {
	jvar->err = E_ALLOC;
	goto johansen_exit;
    }

    /* try to respect the chosen sample period: don't limit the
       generation of lags unnecessarily */
    pdinfo->t1 -= (order - 1);
    if (pdinfo->t1 < 0) {
	pdinfo->t1 = 0;
    }

    /* put first lags of endog vars into "levels_list" */
    for (i=1; i<=list[0]; i++) {
	resids.levels_list[i] = laggenr(list[i], 1, pZ, pdinfo);
	if (resids.levels_list[i] < 0) {
	    jvar->err = E_ALLOC;
	    goto johansen_exit;
	}
    }

#if VAR_DEBUG
    printlist(resids.levels_list, "resids.levels_list (first lags)");
#endif

    /* put first differences of endog vars into VAR list */
    k = 1;
    for (i=1; i<=list[0]; i++) {
	varlist[k] = diffgenr(list[i], pZ, pdinfo, 0);
	if (varlist[k] < 0) {
	    jvar->err = E_ALLOC;
	    goto johansen_exit;
	} 
	k++;
    }

#if VAR_DEBUG
    printlist(varlist, "varlist (first differences)");
#endif

    if (nexo > 0) {
	/* add separator before exog vars */
	varlist[k++] = LISTSEP;
	jvar->jinfo->nexo = nexo;
    }

    if (exolist != NULL) {
	/* add specified exogenous variables to list */
	for (i=1; i<=exolist[0]; i++) {
	    varlist[k++] = exolist[i];
	}
    }	    

    if (seasonals) {
	/* add seasonal dummies to list */
	for (i=0; i<pdinfo->pd-1; i++) {
	    varlist[k++] = di0 + i;
	}
    }

    if (jcode(jvar) == J_UNREST_TREND) {
	/* add trend to VAR list */
	varlist[k++] = gettrend(pZ, pdinfo, 0);
    }	

    if (jcode(jvar) >= J_UNREST_CONST) {
	/* add the constant to the VAR list */
	varlist[k++] = 0;
    }

    if (opt & OPT_V) {
	varprn = prn;
    } 

    jvar->err = johansen_VAR(order - 1, varlist, pZ, pdinfo, 
			     jvar, &resids, opt, varprn); 

    if (!jvar->err) {
	k = gretl_matrix_rows(resids.u);
	jvar->err = allocate_johansen_sigmas(jvar->jinfo, k);
    }

    if (!jvar->err) {
	char stobs[OBSLEN], endobs[OBSLEN];

	jvar->t1 = resids.t1;
	jvar->t2 = resids.t2;
	jvar->T = jvar->t2 - jvar->t1 + 1;

	gretl_matrix_multiply_mod(resids.u, GRETL_MOD_NONE,
				  resids.u, GRETL_MOD_TRANSPOSE,
				  jvar->jinfo->Suu);
	gretl_matrix_multiply_mod(resids.v, GRETL_MOD_NONE,
				  resids.v, GRETL_MOD_TRANSPOSE,
				  jvar->jinfo->Svv);
	gretl_matrix_multiply_mod(resids.u, GRETL_MOD_NONE,
				  resids.v, GRETL_MOD_TRANSPOSE,
				  jvar->jinfo->Suv);

	gretl_matrix_divide_by_scalar(jvar->jinfo->Suu, jvar->T);
	gretl_matrix_divide_by_scalar(jvar->jinfo->Svv, jvar->T);
	gretl_matrix_divide_by_scalar(jvar->jinfo->Suv, jvar->T);

	if (jrank(jvar) == 0) {
	    pprintf(prn, "%s:\n", _("Johansen test"));
	    pprintf(prn, "%s = %d\n", _("Number of equations"), k);
	    pprintf(prn, "%s = %d\n", _("Lag order"), order);
	    pprintf(prn, "%s: %s - %s (T = %d)\n", _("Estimation period"),
		    ntodate(stobs, resids.t1, pdinfo), 
		    ntodate(endobs, resids.t2, pdinfo), jvar->T);
	}

	if (opt & OPT_V) {
	    print_johansen_sigmas(jvar->jinfo, prn);
	}

	if (jrank(jvar) > 0) {
	    /* steal the resids for VECM analysis */
	    jvar->jinfo->u = resids.u;
	    if (restricted(jvar)) {
		jvar->err = split_v_resids(jvar, resids.v);
	    } else {
		jvar->jinfo->v = resids.v;
	    }
	    resids.u = NULL;
	    resids.v = NULL;
	    /* and allocate models for filling out
	       by Johansen plugin */
	    if (!jvar->err) {
		jvar->models = allocate_VAR_models(jvar->neqns);
		if (jvar->models == NULL) {
		    jvar->err = E_ALLOC;
		}
	    }
	}

	/* now get johansen plugin to finish the job */
	if (!jvar->err) {
	    jvar->err = johansen_complete(jvar, pZ, pdinfo, prn);
	}
    } 

johansen_exit:

    var_resids_free(&resids);
    free(varlist);
    
    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    if (jvar->err || !(opt & OPT_S)) {
	dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
    }

    return jvar;
}

/**
 * johansen_test:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt:
 * @prn: gretl printing struct.
 *
 *
 * Returns: pointer to struct containing information on 
 * the test.
 */

GRETL_VAR *johansen_test (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn)
{
    return johansen_driver(order, 0, list, NULL, pZ, pdinfo, opt, prn);
}

/**
 * johansen_test_simple:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt:
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int johansen_test_simple (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar;
    int err;

    jvar = johansen_driver(order, 0, list, NULL, pZ, pdinfo, opt, prn);
    if (jvar == NULL) {
	err = E_ALLOC;
    } else {
	err = jvar->err;
    }

    if (jvar != NULL) {
	gretl_VAR_free(jvar);
    }

    return err;
}

/**
 * vecm:
 * @order: lag order for test.
 * @rank: pre-specified cointegration rank.
 * @list: list of variables to test for cointegration.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt:
 * @prn: gretl printing struct.
 *
 *
 * Returns: pointer to struct containing information on 
 * the VECM system.
 */

GRETL_VAR *vecm (int order, int rank, int *list, 
		 double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar = NULL;
    int *endo_list = NULL, *exo_list = NULL;
    const int *vecm_list = list;
    int err = 0;

    gretl_list_purge_const(list);

    if (gretl_list_has_separator(list)) {
	err = gretl_list_split_on_separator(list, &endo_list, &exo_list);
	if (err) {
	    return jvar;
	}
	vecm_list = endo_list;
    } 

    if (rank <= 0 || rank > list[0]) {
	sprintf(gretl_errmsg, _("vecm: rank %d is out of bounds"), rank);
	return jvar;
    }

    jvar = johansen_driver(order, rank, vecm_list, exo_list,
			   pZ, pdinfo, opt | OPT_S, prn);
    
    if (jvar != NULL && !jvar->err) {
	gretl_VAR_print(jvar, pdinfo, OPT_NONE, prn);
    }

    free(endo_list);
    free(exo_list);

    return jvar;
}

/**
 * vecm_simple:
 * @order: lag order for test.
 * @rank: pre-specified cointegration rank.
 * @list: list of variables to test for cointegration.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt:
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int vecm_simple (int order, int rank, int *list, 
		 double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar;
    int err = 0;

    jvar = vecm(order, rank, list, pZ, pdinfo, opt, prn);

    if (jvar == NULL) {
	err = E_ALLOC;
    } else {
	err = jvar->err;
    }

    if (jvar != NULL) {
	gretl_VAR_free(jvar);
    }

    return err;
}

/**
 * gretl_VAR_print_VCV:
 * @var:
 * @prn:
 *
 *
 * Returns:
 */

int gretl_VAR_print_VCV (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->S == NULL) {
	err = 1;
    } else {
	print_contemp_covariance_matrix(var->S, var->ldet, prn);
    }

    return err;
}

static void tex_print_double (double x, PRN *prn)
{
    char number[16];

    x = screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);

    if (x < 0.) {
	pprintf(prn, "$-$%s", number + 1);
    } else {
	pputs(prn, number);
    }
}

/* printing of impulse responses and variance decompositions */

#define IRF_ROW_MAX 4
#define VDC_ROW_MAX 5

enum {
    IRF,
    VDC
};

static void VAR_RTF_row_spec (int ncols, PRN *prn)
{
    int lcol = 800, colwid = 1600;
    int i, cellx = lcol;

    pputs(prn, "{\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");
    for (i=0; i<ncols; i++) {
	pprintf(prn, "\\cellx%d", cellx);
	cellx += colwid;
    }
    pputc(prn, '\n');
}

static void VAR_info_header_block (int code, int v, int block, 
				   const DATAINFO *pdinfo, 
				   PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[16];

    if (tex) {
	pputs(prn, "\\vspace{1em}\n\n");
	if (code == IRF) {
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    tex_escape(vname, pdinfo->varname[v]));
	} else {
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    tex_escape(vname, pdinfo->varname[v]));
	}
	if (block == 0) {
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, " (%s)\n\n", I_("continued"));
	}
	pprintf(prn, "\\vspace{1em}\n\n\\begin{longtable}{%s}\n",
		(code == IRF)? "rrrrr" : "rrrrrr");
    } else if (rtf) {
	pputs(prn, "\\par\n\n");
	if (code == IRF) {
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[v]);
	} else {
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    pdinfo->varname[v]);
	}
	if (block == 0) {
	    pputs(prn, "\\par\n\n");
	} else {
	    pprintf(prn, " (%s)\\par\n\n", I_("continued"));
	}
	/* FIXME */
	VAR_RTF_row_spec((code == IRF)? IRF_ROW_MAX : VDC_ROW_MAX, prn);
    } else {
	if (code == IRF) {	
	    pprintf(prn, _("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[v]);
	} else {
	    pprintf(prn, _("Decomposition of variance for %s"), 
		    pdinfo->varname[v]);
	}
	if (block == 0) {
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, " (%s)\n\n", _("continued"));
	}
    }

    /* first column: period number header */
    if (tex) {
	pprintf(prn, "%s & ", I_("period"));
    } else if (rtf) {
	pprintf(prn, "\\intbl \\qc %s\\cell ", I_("period"));
    } else {
	pputs(prn, _("period"));
    }
}

static void VAR_info_print_vname (int code, int i, int v, int endrow, 
				  const DATAINFO *pdinfo, PRN *prn)
{
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    char vname[16];

    if (tex) {
	pprintf(prn, " %s ", tex_escape(vname, pdinfo->varname[v]));
	if (endrow) {
	   pputs(prn, "\\\\");
	} else { 
	    pputs(prn, "& ");
	}
    } else if (rtf) {
	pprintf(prn, "\\qc %s\\cell", pdinfo->varname[v]);
	if (endrow) {
	    pputs(prn, " \\intbl \\row");
	} 
    } else {
	int w = (code == IRF)? 13 : 11;

	if (code == IRF && i == 0) w--; /* FIXME width for i18n */
	pprintf(prn, "%*s", w, pdinfo->varname[v]);
    }
}

static void VAR_info_print_period (int t, PRN *prn)
{
    if (tex_format(prn)) {
	pprintf(prn, "%d & ", t);
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\intbl \\qc %d\\cell ", t);
    } else {
	pprintf(prn, " %3d  ", t);
    }
}

static void VAR_info_end_row (PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\intbl \\row\n");
    } else {
	pputc(prn, '\n');
    }
}

static void VAR_info_end_table (PRN *prn)	
{
    if (tex_format(prn)) {
	pputs(prn, "\\end{longtable}\n\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "}\n");
    } else {
	pputc(prn, '\n');
    }
}

#define IRF_ROW_MAX 4
#define VDC_ROW_MAX 5

/**
 * gretl_VAR_print_impulse_response:
 * @var: pointer to VAR struct.
 * @shock:
 * @periods: number of periods over which to print response.
 * @pdinfo: dataset information.
 * @pause: if non-zero, pause between sections of output.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  int pause, PRN *prn)
{
    int i, t;
    int vsrc;
    int rows = var->neqns;
    gretl_matrix *rtmp, *ctmp;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (var->order > 1) {
	rows *= var->order;
    }

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return 1;
    }  

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) {
	return E_ALLOC;
    }

    ctmp = gretl_matrix_alloc(rows, var->neqns);
    if (ctmp == NULL) {
	gretl_matrix_free(rtmp);
	return E_ALLOC;
    }

    if (var->ci == VECM) {
	vsrc = var->jinfo->list[shock + 1];
    } else {
	vsrc = (var->models[shock])->list[1];
    }

    blockmax = var->neqns / IRF_ROW_MAX;
    if (var->neqns % IRF_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax && !err; block++) {
	int k, vtarg, endrow;
	double r;

	VAR_info_header_block(IRF, vsrc, block, pdinfo, prn);

	for (i=0; i<IRF_ROW_MAX; i++) {
	    k = IRF_ROW_MAX * block + i;
	    if (k >= var->neqns) {
		break;
	    }
	    if (var->ci == VECM) {
		vtarg = var->jinfo->list[k + 1];
	    } else {
		vtarg = (var->models[k])->list[1];
	    }
	    endrow = !(i < IRF_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(IRF, i, vtarg, endrow, pdinfo, prn);
	}

	if (tex || rtf) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}

	for (t=0; t<periods && !err; t++) {
	    VAR_info_print_period(t + 1, prn);
	    if (t == 0) {
		/* calculate initial estimated responses */
		err = gretl_matrix_copy_values(rtmp, var->C);
	    } else {
		/* calculate further estimated responses */
		err = gretl_matrix_multiply(var->A, rtmp, ctmp);
		gretl_matrix_copy_values(rtmp, ctmp);
	    }

	    if (err) break;

	    /* matrix rtmp holds the responses */

	    for (i=0; i<IRF_ROW_MAX; i++) {
		k = IRF_ROW_MAX * block + i;
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(rtmp, k, shock);
		if (tex) {
		    tex_print_double(r, prn);
		    if (i < IRF_ROW_MAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.5g\\cell ", r);
		} else {
		    pprintf(prn, "%#12.5g ", r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (rtmp != NULL) gretl_matrix_free(rtmp);
    if (ctmp != NULL) gretl_matrix_free(ctmp);

    return err;
}

int gretl_VAR_print_all_impulse_responses (GRETL_VAR *var, const DATAINFO *pdinfo, 
					   int horizon, PRN *prn)
{
    int i, pause = 0, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (plain_format(prn)) {
	pause = gretl_get_text_pause();
    } else if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_impulse_response(var, i, horizon, pdinfo, 
					       pause, prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }   

    return err;
}

#define VD_ROW_MAX 5

/**
 * gretl_VAR_print_fcast_decomp:
 * @var: pointer to VAR struct.
 * @targ:
 * @periods: number of periods over which to print decomposition.
 * @pdinfo: dataset information.
 * @pause: if non-zero, pause between sections of output.
 * @prn: gretl printing struct.
 *
 *
 * Returns: 0 on success, non-zero code on error.
 */

int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      int pause, PRN *prn)
{
    int i, t;
    int vtarg;
    gretl_matrix *vd = NULL;
    int block, blockmax;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return 1;
    } 

    vd = gretl_VAR_get_fcast_decomp(var, targ, periods);
    if (vd == NULL) {
	return E_ALLOC;
    }

    if (var->ci == VECM) {
	vtarg = var->jinfo->list[targ + 1];
    } else {
	vtarg = (var->models[targ])->list[1];
    }

    blockmax = (var->neqns + 1) / VDC_ROW_MAX;
    if ((var->neqns + 1) % VDC_ROW_MAX) {
	blockmax++;
    }

    for (block=0; block<blockmax; block++) {
	int k, vsrc, endrow;
	double r;

	VAR_info_header_block(VDC, vtarg, block, pdinfo, prn);

	for (i=0; i<VDC_ROW_MAX; i++) {
	    k = VDC_ROW_MAX * block + i - 1;
	    if (k < 0) {
		if (tex) {
		    pprintf(prn, " %s & ", I_("std. error"));
		} else if (rtf) {
		    pprintf(prn, " \\qc %s\\cell ", I_("std. error"));
		} else {
		    pprintf(prn, " %12s ", _("std. error"));
		}
		continue;
	    }
	    if (k >= var->neqns) {
		break;
	    }
	    if (var->ci == VECM) {
		vsrc = var->jinfo->list[k + 1];
	    } else {
		vsrc = (var->models[k])->list[1];
	    }
	    endrow = !(i < VDC_ROW_MAX - 1 && k < var->neqns - 1);
	    VAR_info_print_vname(VDC, i, vsrc, endrow, pdinfo, prn);
	}

	if (tex || rtf) {
	    pputc(prn, '\n');
	} else {
	    pputs(prn, "\n\n");
	}

	for (t=0; t<periods && !err; t++) {
	    VAR_info_print_period(t + 1, prn);
	    for (i=0; i<VDC_ROW_MAX; i++) {
		k = VDC_ROW_MAX * block + i - 1;
		if (k < 0) {
		    r = gretl_matrix_get(vd, t, var->neqns);
		    if (tex) {
			pprintf(prn, "%g & ", r);
		    } else if (rtf) {
			pprintf(prn, "\\qc %g\\cell", r);
		    } else {
			pprintf(prn, " %14g ", r);
		    }
		    continue;
		}
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(vd, t, k);
		if (tex) {
		    pprintf(prn, "$%.4f$", r);
		    if (i < VDC_ROW_MAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.4f\\cell", r);
		} else {
		    pprintf(prn, "%10.4f ", r);
		}
	    }

	    VAR_info_end_row(prn);
	}

	VAR_info_end_table(prn);

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (vd != NULL) {
	gretl_matrix_free(vd);
    }

    return err;
}

int gretl_VAR_print_all_fcast_decomps (GRETL_VAR *var, const DATAINFO *pdinfo, 
				       int horizon, PRN *prn)
{
    int i, pause = 0, err = 0;

    if (horizon <= 0) {
	horizon = default_VAR_horizon(pdinfo);
    }

    if (plain_format(prn)) {
	pause = gretl_get_text_pause();
    } else if (rtf_format(prn)) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    for (i=0; i<var->neqns && !err; i++) {
	err = gretl_VAR_print_fcast_decomp(var, i, horizon, pdinfo, 
					   pause, prn);
    }

    if (rtf_format(prn)) {
	pputs(prn, "}\n");
    }

    return err;
}

void print_Johansen_test_case (JohansenCode jcode, PRN *prn)
{
    const char *jcase[] = {
	N_("Case 1: No constant"),
	N_("Case 2: Restricted constant"),
	N_("Case 3: Unrestricted constant"),
	N_("Case 4: Restricted trend, unrestricted constant"),
	N_("Case 5: Unrestricted trend and constant")
    };

    if (jcode <= J_UNREST_TREND) {
	if (plain_format(prn)) {
	    pputs(prn, _(jcase[jcode]));
	} else {
	    pputs(prn, I_(jcase[jcode]));
	}
    }
}

static void 
print_VECM_coint_eqns (JohansenInfo *jv, const DATAINFO *pdinfo, PRN *prn)
{
    int rtf = rtf_format(prn);
    char s[16];
    int rows = gretl_matrix_rows(jv->Beta);
    int i, j;
    double x;

    pputs(prn, _("Cointegrating vectors"));
    if (jv->Bse != NULL) {
	pprintf(prn, " (%s)", _("standard errors in parentheses"));
    } 
    gretl_prn_newline(prn);
    gretl_prn_newline(prn);

    for (i=0; i<rows; i++) {
	char vname[16];

	if (i < jv->list[0]) {
	    sprintf(vname, "%s(-1)", pdinfo->varname[jv->list[i+1]]);
	} else if (jv->code == J_REST_CONST) {
	    strcpy(vname, "const");
	} else if (jv->code == J_REST_TREND) {
	    strcpy(vname, "trend");
	}
	if (rtf) {
	    pputs(prn, vname);
	} else {
	    pprintf(prn, "%-10s", vname);
	}

	/* coefficients */
	for (j=0; j<jv->rank; j++) {
	    x = gretl_matrix_get(jv->Beta, i, j);
	    if (jv->Bse == NULL) {
		x /= gretl_matrix_get(jv->Beta, j, j);
	    }
	    if (rtf) {
		pprintf(prn, "\t%#.5g ", x);
	    } else {
		pprintf(prn, "%#12.5g ", x);
	    }
	}
	gretl_prn_newline(prn);

	if (jv->Bse != NULL) {
	    /* standard errors */
	    if (rtf) {
		pputs(prn, "\t");
	    } else {
		pprintf(prn, "%11s", " ");
	    }
	    for (j=0; j<jv->rank; j++) {
		if (i < jv->rank) {
		    x = 0.0;
		} else {
		    x = gretl_matrix_get(jv->Bse, i - jv->rank, j);
		}
		sprintf(s, "(%#.5g)", x);
		if (rtf) {
		    pprintf(prn, "\t%s", s);
		} else {
		    pprintf(prn, "%12s ", s);
		}
	    }
	    gretl_prn_newline(prn);
	}
    }

    gretl_prn_newline(prn);
}

static void print_VECM_omega (GRETL_VAR *jvar, const DATAINFO *pdinfo, PRN *prn)
{
    int rtf = rtf_format(prn);
    int *list = jvar->jinfo->list;
    char s[32];
    int i, j;

    pprintf(prn, "%s\n", _("Cross-equation covariance matrix"));
    gretl_prn_newline(prn);

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[list[i+1]]);
	if (i == 0) {
	    if (rtf) {
		pprintf(prn, "\t\t%s", s);
	    } else {
		pprintf(prn, "%25s", s);
	    }
	} else {
	    if (rtf) {
		pprintf(prn, "\t%s", s);
	    } else {
		pprintf(prn, "%13s", s);
	    }
	}
    }
    gretl_prn_newline(prn);

    for (i=0; i<jvar->neqns; i++) {
	sprintf(s, "d_%s", pdinfo->varname[list[i+1]]);
	if (rtf) {
	    pputs(prn, s);
	    if (strlen(s) < 8) {
		pputc(prn, '\t');
	    }	    
	} else {
	    pprintf(prn, "%-13s", s);
	}
	for (j=0; j<jvar->neqns; j++) {
	    if (rtf) {
		pprintf(prn, "\t%#.5g", gretl_matrix_get(jvar->S, i, j));
	    } else {
		pprintf(prn, "%#12.5g ", gretl_matrix_get(jvar->S, i, j));
	    }
	}
	gretl_prn_newline(prn);
    }

    gretl_prn_newline(prn);

    pprintf(prn, "%s = %g", _("determinant"), exp(jvar->ldet));
    gretl_prn_newline(prn);
}

static void 
print_vecm_header_info (GRETL_VAR *vecm, PRN *prn)
{
    gretl_prn_newline(prn);
    pprintf(prn, "%s = %d", 
	    (plain_format(prn))? _("Cointegration rank") : I_("Cointegration rank"),
	    jrank(vecm));
    gretl_prn_newline(prn);
    print_Johansen_test_case(jcode(vecm), prn); 
}

/**
 * gretl_VAR_print:
 * @var: pointer to VAR struct.
 * @pdinfo: dataset information.
 * @opt: if includes %OPT_I, include impulse responses; if
 * includes %OPT_F, include forecast variance decompositions.
 * @prn: pointer to printing struct.
 *
 * Prints the models in @var, along with relevant F-tests and
 * possibly impulse responses and variance decompositions.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print (GRETL_VAR *var, const DATAINFO *pdinfo, gretlopt opt, 
		     PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    char Vstr[72];
    int vecm = var->ci == VECM;
    int dfd = var->models[0]->dfd;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int i, j, k, v;
    int pause = 0;

    if (prn == NULL) {
	return 0;
    }

    ntodate(startdate, var->t1, pdinfo);
    ntodate(enddate, var->t2, pdinfo);

    if (rtf) {
	pputs(prn, "{\\rtf1\\par\n\\qc ");
    }

    if (vecm) {
	if (tex || rtf) {
	    sprintf(Vstr, I_("VECM system, lag order %d"), var->order + 1);
	} else {
	    sprintf(Vstr, _("VECM system, lag order %d"), var->order + 1);
	}
    } else {
	if (tex || rtf) {
	    sprintf(Vstr, I_("VAR system, lag order %d"), var->order);
	} else {
	    sprintf(Vstr, _("VAR system, lag order %d"), var->order);
	}
    }

    if (tex) {
	pputs(prn, "\\begin{center}");
	pprintf(prn, "\n%s\\\\\n", Vstr);
	pprintf(prn, I_("%s estimates, observations %s--%s ($T=%d$)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, prn);
	}
	pputs(prn, "\n\\end{center}\n");
    } else if (rtf) {
	gretl_print_toggle_doc_flag(prn);
	pprintf(prn, "\n%s\\par\n", Vstr);
	pprintf(prn, I_("%s estimates, observations %s%-%s (T = %d)"),
		(vecm)? I_("Maximum likelihood") : I_("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, prn);
	}	
	pputs(prn, "\\par\n\n");
    } else {
	pause = gretl_get_text_pause();
	pprintf(prn, "\n%s\n", Vstr);
	pprintf(prn, _("%s estimates, observations %s-%s (T = %d)"),
		(vecm)? ("Maximum likelihood") : _("OLS"), startdate, enddate, var->T);
	if (vecm) {
	    print_vecm_header_info(var, prn);
	}
	pputs(prn, "\n\n");
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_coint_eqns(var, pdinfo, prn);
	} else {
	    print_VECM_coint_eqns(var->jinfo, pdinfo, prn);
	}
    }

    if (tex) {
	tex_print_VAR_ll_stats(var, prn);
    } else if (rtf) {
	pprintf(prn, "%s = %#g\\par\n", I_("Log-likelihood"), var->ll);
	pprintf(prn, "%s = %#g\\par\n", I_("Determinant of covariance matrix"), 
		exp(var->ldet));
	pprintf(prn, "%s = %.4f\\par\n", I_("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\\par\n", I_("BIC"), var->BIC);
    } else {
	pprintf(prn, "%s = %#g\n", _("Log-likelihood"), var->ll);
	pprintf(prn, "%s = %#g\n", _("Determinant of covariance matrix"), exp(var->ldet));
	pprintf(prn, "%s = %.4f\n", _("AIC"), var->AIC);
	pprintf(prn, "%s = %.4f\n", _("BIC"), var->BIC);
    }

    if (vecm) {
	pputc(prn, '\n');
    }

    k = 0;

    for (i=0; i<var->neqns; i++) {
	char Fstr[24];

	printmodel(var->models[i], pdinfo, OPT_NONE, prn);

	if (pause) {
	    scroll_pause();
	}

	if (vecm) {
	    continue;
	}

	if (tex) {
	    pputs(prn, "\n\\begin{center}\n");
	    pprintf(prn, "%s\\\\[1em]\n", I_("F-tests of zero restrictions"));
	    pputs(prn, "\\begin{tabular}{lll}\n");
	} else if (rtf) {
	    pprintf(prn, "%s:\\par\n\n", I_("F-tests of zero restrictions"));
	} else {
	    pprintf(prn, "  %s:\n", _("F-tests of zero restrictions"));
	}

	for (j=0; j<var->neqns; j++) {
	    v = (var->models[j])->list[1];
	    if (tex) {
		pprintf(prn, I_("All lags of %-8s "), pdinfo->varname[v]);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\\\\n", I_("p-value"), 
			fdist(var->Fvals[k], var->order, dfd));
	    } else if (rtf) {
		pprintf(prn, I_("All lags of %-8s "), pdinfo->varname[v]);
		pprintf(prn, "F(%d, %d) = %10g, ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\par\n", I_("p-value"), 
			fdist(var->Fvals[k], var->order, dfd));
	    } else {
		pputs(prn, "    ");
		pprintf(prn, _("All lags of %-8s "), pdinfo->varname[v]);
		sprintf(Fstr, "F(%d, %d)", var->order, dfd);
		pprintf(prn, "%12s = %#10.5g, ", Fstr, var->Fvals[k]);
		pprintf(prn, "%s %.4f\n", _("p-value"), 
			fdist(var->Fvals[k], var->order, dfd));
	    }
	    k++;
	}

	if (var->order > 1) {
	    if (tex) {
		pprintf(prn, I_("All vars, lag %-6d "), var->order);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\\\\n", I_("p-value"), 
			fdist(var->Fvals[k], var->neqns, dfd));
	    } else if (rtf) {
		pprintf(prn, I_("All vars, lag %-6d "), var->order);
		pprintf(prn, "F(%d, %d) = %10g, ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, "%s %.4f\\par\n", I_("p-value"), 
			fdist(var->Fvals[k], var->neqns, dfd));
	    } else {
		pputs(prn, "    ");
		pprintf(prn, _("All vars, lag %-6d "), var->order);
		sprintf(Fstr, "F(%d, %d)", var->neqns, dfd);
		pprintf(prn, "%12s = %#10.5g, ", Fstr, var->Fvals[k]);
		pprintf(prn, "%s %.4f\n", _("p-value"), 
			fdist(var->Fvals[k], var->neqns, dfd));
	    } 
	    k++;
	}

	if (tex) {
	    pputs(prn, "\\end{tabular}\n"
		  "\\end{center}\n\n"
		  "\\clearpage\n\n");
	} else if (rtf) {
	    pputs(prn, "\\par\\n\n");
	} else if (pause) {
	    scroll_pause();
	}
    }

    pputc(prn, '\n');

    /* global LR test on max lag */
    if (!na(var->LR)) {
	char h0str[64];
	char h1str[64];
	int df = var->neqns * var->neqns;

	if (tex || rtf) {
	    sprintf(h0str, I_("the longest lag is %d"), var->order - 1);
	    sprintf(h1str, I_("the longest lag is %d"), var->order);
	} else {
	    sprintf(h0str, _("the longest lag is %d"), var->order - 1);
	    sprintf(h1str, _("the longest lag is %d"), var->order);
	}	    

	if (tex) {
	    pprintf(prn, "\\noindent %s ---\\par\n", I_("For the system as a whole"));
	    pprintf(prn, "%s: %s\\par\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "%s: %s\\par\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "%s: $\\chi^2_{%d}$ = %.3f (%s %f)\\par\n",
		    I_("Likelihood ratio test"), 
		    df, var->LR, I_("p-value"), chisq(var->LR, df));
	} else if (rtf) {
	    pprintf(prn, "\\par %s\n", I_("For the system as a whole"));
	    pprintf(prn, "\\par %s: %s\n", I_("Null hypothesis"), h0str);
	    pprintf(prn, "\\par %s: %s\n", I_("Alternative hypothesis"), h1str);
	    pprintf(prn, "\\par %s: %s(%d) = %g (%s %f)\n",
		    I_("Likelihood ratio test"), I_("Chi-square"), 
		    df, var->LR, I_("p-value"), chisq(var->LR, df));
	} else {
	    pprintf(prn, "%s\n", _("For the system as a whole"));
	    pprintf(prn, "  %s: %s\n", _("Null hypothesis"), h0str);
	    pprintf(prn, "  %s: %s\n", _("Alternative hypothesis"), h1str);
	    pprintf(prn, "  %s: %s(%d) = %g (%s %f)\n",
		    _("Likelihood ratio test"), _("Chi-square"), 
		    df, var->LR, _("p-value"), chisq(var->LR, df));
	}
    }

    if (vecm) {
	if (tex_format(prn)) {
	    tex_print_VECM_omega(var, pdinfo, prn);
	} else {
	    print_VECM_omega(var, pdinfo, prn);
	}
    } else {
	pputc(prn, '\n');
    }

    if (opt & OPT_I) {
	gretl_VAR_print_all_impulse_responses(var, pdinfo, 0, prn);
    }

    if (opt & OPT_F) {
	gretl_VAR_print_all_fcast_decomps(var, pdinfo, 0, prn);
    }

    if (rtf) {
	pputs(prn, "}\n");
    }

    return 0;
}

void gretl_VAR_assign_name (GRETL_VAR *var)
{
    static int nvar = 0;

    if (var->name != NULL) {
	free(var->name);
    }

    var->name = malloc(10);

    if (var->name != NULL) {
	if (var->ci == VAR) {
	    sprintf(var->name, "%s %d", _("VAR"), ++nvar);
	} else {
	    sprintf(var->name, "%s %d", _("VECM"), gretl_VECM_id(var));
	}
    }
}

void gretl_VAR_assign_specific_name (GRETL_VAR *var, const char *name)
{
    if (var->name != NULL) {
	free(var->name);
    }

    var->name = gretl_strdup(name);
}

const char *gretl_VAR_get_name (const GRETL_VAR *var)
{
    return var->name;
}

int gretl_VAR_add_resids_to_dataset (GRETL_VAR *var, int eqnum,
				     double ***pZ, DATAINFO *pdinfo)
{
    MODEL *pmod = var->models[eqnum];
    int i, t;

    if (dataset_add_series(1, pZ, pdinfo)) return E_ALLOC;

    i = pdinfo->v - 1;

    for (t=0; t<pdinfo->n; t++) {
	if (t < pmod->t1 || t > pmod->t2) {
	    (*pZ)[i][t] = NADBL;
	} else {
	    (*pZ)[i][t] = pmod->uhat[t];
	}
    }

    sprintf(pdinfo->varname[i], "uhat%d", eqnum + 1);

    if (var->ci == VAR) {
	sprintf(VARLABEL(pdinfo, i), _("residual from VAR system, equation %d"), 
		eqnum + 1);
    } else {
	sprintf(VARLABEL(pdinfo, i), _("residual from VECM system, equation %d"), 
		eqnum + 1);
    }

    return 0;
}

int gretl_VAR_get_highest_variable (const GRETL_VAR *var,
				    const DATAINFO *pdinfo)
{
    int vmax = 0;

    if (var->models != NULL && var->neqns >= 1) {
	vmax = highest_numbered_var_in_model(var->models[0], pdinfo);
    }

    return vmax;
}

int gretl_VECM_id (GRETL_VAR *vecm)
{
    static int nvecm;

    if (vecm->jinfo->ID == 0) {
	vecm->jinfo->ID = ++nvecm;
    }

    return vecm->jinfo->ID;
}

#include "irfboot.c"
