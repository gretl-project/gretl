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
#include "varprint.h"
#include "libset.h"
#include "transforms.h"
#include "gretl_xml.h"

#define VAR_DEBUG 0
#define BDEBUG    0  /* for debugging bootstrap IRFs */

#define N_IVALS 3

static gretl_matrix *irf_bootstrap (const GRETL_VAR *var, 
				    int targ, int shock, int periods,
				    const double **Z, 
				    const DATAINFO *pdinfo);
static gretl_matrix *VAR_coeff_matrix_from_VECM (const GRETL_VAR *var);

struct var_lists {
    int *detvars;
    int *stochvars;
    int *reglist;
    int *testlist;
    int **lagvlist;
};

static int VAR_list_laggenr (const int *list, double ***pZ, DATAINFO *pdinfo,
			     int order, int **lagnums)
{
    int *tmplist = gretl_list_copy(list);
    int i, j, k;
    int err = 0;
    
    if (tmplist == NULL) {
	return E_ALLOC;
    }

    err = list_laggenr(&tmplist, order, pZ, pdinfo);

    if (!err && tmplist[0] < list[0] * order) {
	/* we didn't get all the lags we wanted */
	err = E_DATA;
    }

    if (!err) {
	k = 1;
	for (i=0; i<list[0]; i++) {
	    for (j=1; j<=order; j++) {
		lagnums[i][j] = tmplist[k++];
	    }
	}
    }

    free(tmplist);
    
    return err;
}

static int gretl_VAR_add_models (GRETL_VAR *var)
{
    int n = var->neqns;
    int i, err = 0;

    if (var->models != NULL) {
	for (i=0; i<n; i++) {
	    clear_model(var->models[i]);
	}
    } else {
	var->models = gretl_model_array_new(n);
	if (var->models == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static int VAR_add_companion_matrix (GRETL_VAR *var)
{
    int n = var->neqns * (var->order + var->ecm);
    int i, j, err = 0;

    if (var->A != NULL) {
	return 0;
    }

    var->A = gretl_matrix_alloc(n, n);

    if (var->A == NULL) {
	err = E_ALLOC;
    } else {
	for (i=var->neqns; i<n; i++) {
	    for (j=0; j<n; j++) {
		gretl_matrix_set(var->A, i, j, (j == i - var->neqns)? 1.0 : 0.0);
	    }
	}
    }

    return err;
}

static int VAR_add_cholesky_matrix (GRETL_VAR *var)
{
    int n = var->neqns * (var->order + var->ecm);
    int err = 0;

    if (var->C != NULL) {
	return 0;
    }    

    var->C = gretl_matrix_alloc(n, var->neqns);
    if (var->C == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_zero(var->C);
    }

    return err;
}

static int VAR_add_residuals_matrix (GRETL_VAR *var)
{
    int err = 0;

    if (var->E != NULL) {
	return 0;
    }      

    var->E = gretl_matrix_alloc(var->T, var->neqns);
    if (var->E == NULL) {
	err = E_ALLOC;
    } 

    return err;
}

static void gretl_VAR_zero (GRETL_VAR *var)
{
    var->err = 0;
    var->ifc = 0;
    var->ncoeff = 0;
    var->t1 = var->t2 = var->T = 0;
    
    var->A = NULL;
    var->lambda = NULL;
    var->E = NULL;
    var->C = NULL;
    var->S = NULL;
    var->F = NULL;

    var->models = NULL;
    var->Fvals = NULL;
    var->Ivals = NULL;
    var->name = NULL;

    var->ll = var->ldet = NADBL;
    var->AIC = var->BIC = var->HQC = NADBL;
    var->LR = NADBL;
}

static GRETL_VAR *gretl_VAR_new (int ci, int neqns, int order)
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

    var->ci = ci;
    var->refcount = 0;

    var->neqns = neqns;
    var->order = order;
    var->ecm = 0;

    gretl_VAR_zero(var);
    var->jinfo = NULL;

    err = VAR_add_companion_matrix(var);

    if (!err) {
	err = VAR_add_cholesky_matrix(var);
    }

    if (!err) {
	err = gretl_VAR_add_models(var);
    }

    if (!err) {
	int m = neqns * neqns + neqns;
	
	var->Fvals = malloc(m  * sizeof *var->Fvals);
	if (var->Fvals == NULL) {
	    err = 1;
	}
    }

    if (!err) {
	var->Ivals = malloc(N_IVALS  * sizeof *var->Ivals);
	if (var->Ivals == NULL) {
	    err = 1;
	}
    }

    if (err) {
	gretl_VAR_free(var);
	var = NULL;
    }

    return var;
}

static void johansen_info_free (JohansenInfo *jv)
{
    free(jv->list);
    free(jv->exolist);
    free(jv->difflist);
    free(jv->biglist);
    free(jv->levels_list);
    free(jv->varlist);

    gretl_matrix_free(jv->u);
    gretl_matrix_free(jv->v);

    gretl_matrix_free(jv->Suu);
    gretl_matrix_free(jv->Svv);
    gretl_matrix_free(jv->Suv);

    gretl_matrix_free(jv->Beta);
    gretl_matrix_free(jv->Alpha);
    gretl_matrix_free(jv->Bse);
    gretl_matrix_free(jv->Bvar);
    gretl_matrix_free(jv->D);

    free(jv);
}

void gretl_VAR_free (GRETL_VAR *var)
{
    if (var == NULL) return;

#if 0
    fprintf(stderr, "gretl_VAR_free: var = %p, refcount = %d\n",
	    (void *) var, var->refcount);
#endif

    var->refcount -= 1;
    if (var->refcount > 0) {
	return;
    }

    gretl_matrix_free(var->A);
    gretl_matrix_free(var->lambda);
    gretl_matrix_free(var->E);
    gretl_matrix_free(var->C);
    gretl_matrix_free(var->S);
    gretl_matrix_free(var->F);

    free(var->Fvals);
    free(var->Ivals);
    free(var->name);

    if (var->models != NULL) {
	gretl_model_array_destroy(var->models, var->neqns);
    }

    if (var->jinfo != NULL) {
	johansen_info_free(var->jinfo);
    }

    free(var);

#if 0
    fprintf(stderr, "gretl_VAR_free: done\n");
    fflush(stderr);
#endif
}

static int add_VAR_fcast_variance (GRETL_VAR *var, gretl_matrix *F,
				   int nf, int pre_obs)
{
    double vti;
    int i, s;
    int err = 0;

    for (i=0; i<var->neqns; i++) {
	gretl_matrix *vd;
	int totcol;

	vd = gretl_VAR_get_fcast_decomp(var, i, nf - pre_obs);
	if (vd != NULL) {
	    totcol = gretl_matrix_cols(vd) - 1;
	    for (s=0; s<nf; s++) {
		if (s < pre_obs) {
		    gretl_matrix_set(F, s, var->neqns + i, NADBL);
		} else {
		    vti = gretl_matrix_get(vd, s - pre_obs, totcol);
		    gretl_matrix_set(F, s, var->neqns + i, vti);
		}
	    }
	    gretl_matrix_free(vd);
	} else {
	    err = E_ALLOC;
	    for (s=0; s<nf; s++) {
		gretl_matrix_set(F, s, var->neqns + i, NADBL);
	    }
	}
    }

    return err;
}

static int
gretl_VAR_add_forecast (GRETL_VAR *var, int t0, int t1, int t2,
			const double **Z, const DATAINFO *pdinfo, 
			gretlopt opt)
{
    const MODEL *pmod;
    gretl_matrix *F;
    double fti, xti;
    int nf = t2 - t0 + 1;
    int staticfc = (opt & OPT_S);
    int i, j, k, s, t;
    int ns, lag, vj, m;
    int fcols;

    fcols = (staticfc)? var->neqns : 2 * var->neqns;

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

    for (t=t0, s=0; t<=t2; t++, s++) {
	int miss = 0;

	for (i=0; i<var->neqns; i++) {
	    pmod = var->models[i];
	    fti = 0.0;
	    lag = 1;
	    k = 0;
	    for (j=0; j<pmod->ncoeff; j++) {
		vj = pmod->list[j + 2];
		if (j < ns + pmod->ifc && vj > 0) {
		    /* stochastic var */
		    if (t < t1 || staticfc || s - lag < 0) {
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
	add_VAR_fcast_variance(var, F, nf, t1 - t0);
    }

    gretl_matrix_set_t1(F, t0);
    var->F = F;

#if 0
    gretl_matrix_print(F, "var->F");
#endif

    return 0;
}

static int
gretl_VECM_add_forecast (GRETL_VAR *var, int t0, int t1, int t2,
			 const double **Z, DATAINFO *pdinfo, 
			 gretlopt opt)
{
    gretl_matrix *F;
    gretl_matrix *B;

    int order = var->order + 1;
    int nexo = (var->jinfo->exolist != NULL)? var->jinfo->exolist[0] : 0;
    int nseas = var->jinfo->seasonals;
    int nf = t2 - t0 + 1;
    int staticfc = (opt & OPT_S);
    int i, j, k, vj, s, t;
    int fcols, d0 = 0;

    if (nseas > 0) {
	/* find out where the seasonal dummies are */
	d0 = dummy(NULL, pdinfo, 1);
	if (d0 < 0) {
	    return E_DATA;
	}
    }

    fcols = (staticfc)? var->neqns : 2 * var->neqns;
    F = gretl_matrix_alloc(nf, fcols);
    if (F == NULL) {
	return E_ALLOC;
    }

    B = VAR_coeff_matrix_from_VECM(var);
    if (B == NULL) {
	gretl_matrix_free(F);
	return E_ALLOC;
    }

    gretl_matrix_zero(F);

    for (t=t0, s=0; t<=t2; t++, s++) {

	for (i=0; i<var->neqns; i++) {
	    double y, bij, fti = 0.0;
	    int ft, col = 0;

	    /* unrestricted constant, if present */
	    if (var->ifc) {
		bij = gretl_matrix_get(B, i, col++);
		fti += bij;
	    }

	    /* lags of endogenous vars */
	    for (j=0; j<var->neqns; j++) {
		vj = var->jinfo->list[j+1];
		for (k=0; k<order; k++) {
		    if (t - k - 1 < 0) {
			fti = NADBL;
			break;
		    }			
		    bij = gretl_matrix_get(B, i, col++);
		    ft = s - k - 1;
		    if (t >= t1 && ft >= 0 && !staticfc) {
			/* use prior forecast if available */
			y = gretl_matrix_get(F, ft, j);
		    } else {
			y = Z[vj][t-k-1];
		    }
		    if (na(y)) {
			fti = NADBL;
			break;
		    } else {
			fti += bij * y;
		    }
		}
		if (na(fti)) {
		    break;
		}
	    }

	    if (na(fti)) goto set_fcast;

	    /* exogenous vars, if present */
	    for (j=0; j<nexo; j++) {
		vj = var->jinfo->exolist[j+1];
		bij = gretl_matrix_get(B, i, col++);
		if (na(Z[vj][t])) {
		    fti = NADBL;
		} else {
		    fti += bij * Z[vj][t];
		}
	    }

	    if (na(fti)) goto set_fcast;

	    /* seasonals, if present */
	    for (j=0; j<nseas; j++) {
		bij = gretl_matrix_get(B, i, col++);
		fti += bij * Z[d0+j][t];
	    }

	    if (jcode(var) == J_UNREST_TREND) {
		/* unrestricted trend */
		bij = gretl_matrix_get(B, i, col);
		fti += bij * (t + 1);
	    } else if (jcode(var) == J_REST_CONST) {
		/* restricted constant */
		fti += gretl_matrix_get(B, i, col);
	    } else if (jcode(var) == J_REST_TREND) {
		/* restricted trend */
		bij = gretl_matrix_get(B, i, col);
		fti += bij * t;
	    }
	    
	set_fcast:

	    gretl_matrix_set(F, s, i, fti);
	}
    }

    gretl_matrix_free(B);

    /* now get variances, if not static */
    if (!staticfc) {
	add_VAR_fcast_variance(var, F, nf, t1 - t0);
    }

    gretl_matrix_set_t1(F, t0);
    var->F = F;

#if 0
    gretl_matrix_print(F, "var->F");
#endif

    return 0;
}

const gretl_matrix *
gretl_VAR_get_forecast_matrix (GRETL_VAR *var, int t0, int t1, int t2,
			       const double **Z, DATAINFO *pdinfo, 
			       gretlopt opt)
{
    if (var->F != NULL) {
	int ncols, nf = t2 - t0 + 1;
	int ft1 = gretl_matrix_get_t1(var->F);

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
	if (var->ecm) {
	    gretl_VECM_add_forecast(var, t0, t1, t2, Z, pdinfo, opt);
	} else {
	    gretl_VAR_add_forecast(var, t0, t1, t2, Z, pdinfo, opt);
	}
    }

    return var->F;
}

const gretl_matrix *
gretl_VAR_get_residual_matrix (const GRETL_VAR *var)
{
    return var->E;
}

int gretl_VAR_normality_test (const GRETL_VAR *var, PRN *prn)
{
    int err = 0;

    if (var->E == NULL || var->S == NULL) {
	err = 1;
    } else {
	err = gretl_system_normality_test(var->E, var->S, prn);
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
	err = arch_test_simple(var->models[i], order, pZ, pdinfo, OPT_NONE, prn);
    }

    return err;
}

int gretl_VAR_do_error_decomp (const gretl_matrix *S,
			       gretl_matrix *C)
{
    int g = gretl_matrix_rows(S);
    gretl_matrix *tmp = NULL;
    double x;
    int i, j, err = 0;

    /* copy cross-equation covariance matrix (note: the C matrix has
       more rows than S)
    */
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

    /* write the decomposition (lower triangle) into the C matrix */
    if (!err) {
	for (i=0; i<g; i++) {
	    for (j=0; j<=i; j++) {
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
    int rows = var->neqns * (var->order + var->ecm);
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
	    /* initial estimated responses */
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
	/* no data matrix provided: just return point estimate */
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
    int rows = var->neqns * (var->order + var->ecm);
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
						ctmp, GRETL_MOD_NONE);
		err = gretl_matrix_multiply(var->C, ctmp, cic);
		gretl_matrix_copy_values(vt, cic);
	    } else {
		/* calculate further variances */
		err = gretl_matrix_multiply_mod(vt, GRETL_MOD_NONE,
						var->A, GRETL_MOD_TRANSPOSE,
						vtmp, GRETL_MOD_NONE);
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
	lvlist[i] = gretl_list_new(order);
	if (lvlist[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(lvlist[j]);
	    }
	    free(lvlist);
	    lvlist = NULL;
	}
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
   than the full VAR order (again, for an F-test).  If misstest = 1,
   also compose a list containing the maximum lag of each stochastic
   variable, for use in testing for missing values.
*/

static int
compose_varlist (struct var_lists *vl, int depvar, int order, int omit, 
		 int misstest, const DATAINFO *pdinfo)
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

    /* append the exogenous vars */
    for (i=1; i<=vl->detvars[0]; i++) {
	vl->reglist[pos++] = vl->detvars[i];
    }

#if VAR_DEBUG
    printlist(vl->reglist, "composed VAR list");
#endif

    if (misstest) {
	/* build a test list to screen missing values */
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
	var->ifc = pmod->ifc;
	var->T = pmod->t2 - pmod->t1 + 1;
	err = VAR_add_residuals_matrix(var);
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

static int VAR_add_roots (GRETL_VAR *var)
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

    /* save eigenvalues of companion form matrix */
    if (!err) {
        eigA = gretl_general_matrix_eigenvals(CompForm, 0, &err);
	if (!err) {
	    for (i=0; i<np; i++) {
		x = eigA[i];
		y = eigA[np + i];
		gretl_matrix_set(var->lambda, i, 0, x);
		gretl_matrix_set(var->lambda, i, 1, y);
	    }
#if 0
	    gretl_matrix_print(var->A, "Companion form matrix");
	    gretl_matrix_print(var->lambda, "Eigenvalues");
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
	VAR_add_roots(var);
    }

    return var->lambda;
}

static double gretl_VAR_ldet (GRETL_VAR *var, int *err)
{
    gretl_matrix *S = NULL;
    double ldet = NADBL;

    S = gretl_matrix_alloc(var->neqns, var->neqns);

    if (S == NULL) {
	*err = E_ALLOC;
    } else {
	gretl_matrix_multiply_mod(var->F, GRETL_MOD_TRANSPOSE,
				  var->F, GRETL_MOD_NONE,
				  S, GRETL_MOD_NONE);
	gretl_matrix_divide_by_scalar(S, var->T);
	ldet = gretl_vcv_log_determinant(S);
	if (na(ldet)) {
	    *err = 1;
	}
	gretl_matrix_free(S);
    }

    return ldet;
}

static int VAR_LR_lag_test (GRETL_VAR *var)
{
    double ldet;
    int err = 0;

    ldet = gretl_VAR_ldet(var, &err);

    if (!err) {
	double ll, AIC, BIC, HQC;
	int T = var->T;
	int g = var->neqns;
	int m = var->ncoeff - g;
	int k = g * m;

	var->LR = T * (ldet - var->ldet);

	ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	AIC = (-2.0 * ll + 2.0 * k) / T;
	BIC = (-2.0 * ll + k * log(T)) / T;
	HQC = (-2.0 * ll + 2.0 * k * log(log(T))) / T;
	var->Ivals[0] = AIC;
	var->Ivals[1] = BIC;
	var->Ivals[2] = HQC;
    }

    /* we're done with this set of residuals */
    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

static void gretl_VAR_print_lagsel (gretl_matrix *lltab,
				    gretl_matrix *crittab,
				    int *best_row,
				    PRN *prn)
{
    int maxlag = gretl_matrix_rows(crittab);
    double x;
    int i, j;

    pprintf(prn, _("VAR system, maximum lag order %d"), maxlag);
    pputs(prn, "\n\n");

    pputs(prn, _("The asterisks below indicate the best (that is, minimized) values\n"
	  "of the respective information criteria, AIC = Akaike criterion,\n"
	  "BIC = Schwartz Bayesian criterion and HQC = Hannan-Quinn criterion."));
    pputs(prn, "\n\n");

    pputs(prn, _("lags        loglik    p(LR)       AIC          BIC          HQC"));
    pputs(prn, "\n\n");

    for (i=0; i<maxlag; i++) {
	pprintf(prn, "%4d", i + 1);
	x = gretl_matrix_get(lltab, i, 0);
	pprintf(prn, "%14.5f", x);
	if (i > 0) {
	    x = gretl_matrix_get(lltab, i, 1);
	    pprintf(prn, "%9.5f", x);
	} else {
	    pputs(prn, "         ");
	}
	for (j=0; j<N_IVALS; j++) {
	    x = gretl_matrix_get(crittab, i, j);
	    pprintf(prn, "%12.6f", x);
	    if (i == best_row[j]) {
		pputc(prn, '*');
	    } else {
		pputc(prn, ' ');
	    }
	}
	pputc(prn, '\n');
    }
}

static int 
gretl_VAR_do_lagsel (GRETL_VAR *var, struct var_lists *vl,
		     double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    gretl_matrix *crittab = NULL;
    gretl_matrix *lltab = NULL;
    MODEL testmod;

    /* initialize the "best" at the longest lag */
    double best[N_IVALS] = { var->AIC, var->BIC, var->HQC };
    int r = var->order - 1;
    int best_row[N_IVALS] = { r, r, r };

    double crit[N_IVALS];
    double LRtest;
    int T = var->T;
    int g = var->neqns;

    double ldet = NADBL;
    int depvar;
    int i, j, t;
    int tm, m = 0;
    int err = 0;

    if (var->order < 2) {
	return 0;
    }

    if (var->F != NULL) {
	gretl_matrix_free(var->F);
    }

    var->F = gretl_matrix_alloc(var->T, var->neqns);
    if (var->F == NULL) {
	return E_ALLOC;
    }

    crittab = gretl_matrix_alloc(var->order, N_IVALS);
    lltab = gretl_matrix_alloc(var->order, 2);
    if (crittab == NULL || lltab == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (j=1; j<var->order && !err; j++) {
	for (i=0; i<var->neqns && !err; i++) {
	    depvar = vl->stochvars[i + 1];
	    compose_varlist(vl, depvar, j, 0, 0, pdinfo);
	    testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A);
	    err = testmod.errcode;
	    if (!err) {
		/* record residuals for equation i, order j */
		tm = testmod.t1;
		for (t=0; t<var->T; t++) {
		    gretl_matrix_set(var->F, t, i, testmod.uhat[tm++]);
		}
	    }
	    clear_model(&testmod);
	}
	if (!err) {
	    /* resids matrix should now be complete for order j */
	    ldet = gretl_VAR_ldet(var, &err);
	}
	if (!err) {
	    double ll;
	    int p = var->ncoeff - (g * (var->order - j));
	    int c, k = g * p;

	    ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * ldet;
	    crit[0] = (-2.0 * ll + 2.0 * k) / T;               /* AIC */
	    crit[1] = (-2.0 * ll + k * log(T)) / T;            /* BIC */
	    crit[2] = (-2.0 * ll + 2.0 * k * log(log(T))) / T; /* HQC */

	    gretl_matrix_set(lltab, m, 0, ll);
	    if (j == 1) {
		gretl_matrix_set(lltab, m, 1, 0);
	    } else {
		LRtest = 2.0 * (ll - gretl_matrix_get(lltab, m-1, 0));
		gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, g * g));
	    }	
	    
	    for (c=0; c<N_IVALS; c++) {
		gretl_matrix_set(crittab, m, c, crit[c]);
		if (crit[c] < best[c]) {
		    best[c] = crit[c];
		    best_row[c] = m;
		}
	    }
	
	    m++;
	}
    }

    if (!err) {
	gretl_matrix_set(lltab, m, 0, var->ll);
	LRtest = 2.0 * (var->ll - gretl_matrix_get(lltab, m - 1, 0));
	gretl_matrix_set(lltab, m, 1, chisq_cdf_comp(LRtest, g * g));
	gretl_matrix_set(crittab, m, 0, var->AIC);
	gretl_matrix_set(crittab, m, 1, var->BIC);
	gretl_matrix_set(crittab, m, 2, var->HQC);
	gretl_VAR_print_lagsel(lltab, crittab, best_row, prn);
    }

    bailout:

    gretl_matrix_free(crittab);
    gretl_matrix_free(lltab);

    gretl_matrix_free(var->F);
    var->F = NULL;

    return err;
}

/* Per-equation F-tests for excluding variables and for excluding the
   last lag, plus AIC and BIC comparison for last lag.
*/

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
	    testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A);
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

	testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A);
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
	    /* record residuals for LR test, AIC, BIC, HQC */
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

static int VAR_add_variance_matrix (GRETL_VAR *var)
{
    int err = 0;

    var->S = gretl_matrix_alloc(var->neqns, var->neqns);
    if (var->S == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	gretl_matrix_multiply_mod(var->E, GRETL_MOD_TRANSPOSE,
				  var->E, GRETL_MOD_NONE,
				  var->S, GRETL_MOD_NONE);
	gretl_matrix_divide_by_scalar(var->S, var->T);
    }

    return err;
}

static int VAR_add_stats (GRETL_VAR *var)
{
    int err;

    err = VAR_add_variance_matrix(var);

    if (!err) {
	var->ldet = gretl_vcv_log_determinant(var->S);
	if (na(var->ldet)) {
	    err = 1;
	}
    }    

    if (!err) {
	int T = var->T;
	int g = var->neqns;
	int p = var->ncoeff;
	int k = g * p;

	var->ll = -(g * T / 2.0) * (LN_2_PI + 1) - (T / 2.0) * var->ldet;
	var->AIC = (-2.0 * var->ll + 2.0 * k) / T;
	var->BIC = (-2.0 * var->ll + k * log(T)) / T;
	var->HQC = (-2.0 * var->ll + 2.0 * k * log(log(T))) / T;
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
			    struct var_lists *vl,
			    double ***pZ, DATAINFO *pdinfo,
			    gretlopt opt, int *err)
{
    GRETL_VAR *var = NULL;
    gretlopt lsqopt = OPT_A | OPT_Z;
    int lagtest = (opt & OPT_L);
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
			      order, vl);
    if (*err) {
	return NULL;
    }

    /* generate the required lags */
    if ((*err = VAR_list_laggenr(vl->stochvars, pZ, pdinfo, 
				 order, vl->lagvlist))) {
	goto var_bailout;
    }

    neqns = vl->stochvars[0];    

    /* compose base VAR list (entry 1 will vary across equations);
       assemble test list for t1 and t2 while we're at it */
    *err = compose_varlist(vl, vl->stochvars[1], 
			   order, 0, 1, pdinfo);
    if (*err) {
	*err = E_DATA;
	goto var_bailout;
    }

    /* sort out sample range */
    if (check_for_missing_obs(vl->testlist, &pdinfo->t1, &pdinfo->t2,
			      (const double **) *pZ, NULL)) {
	*err = E_MISSDATA;
	goto var_bailout;
    }

    /* allocate storage */
    var = gretl_VAR_new(VAR, neqns, order);
    if (var == NULL) {
	*err = E_ALLOC;
	goto var_bailout;
    }

    k = 0;

    for (i=0; i<neqns && !*err; i++) {
	MODEL *pmod = var->models[i];

	compose_varlist(vl, vl->stochvars[i + 1], order, 0, 0, pdinfo);

	*pmod = lsq(vl->reglist, pZ, pdinfo, VAR, lsqopt);

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
	    *err = VAR_compute_tests(pmod, var, vl, pZ, pdinfo, i, &k);
	}
    }

 var_bailout:

    if (!*err) {
	var->ncoeff = var->models[0]->ncoeff;
	var->t1 = var->models[0]->t1;
	var->t2 = var->models[0]->t2;
	var->T = var->t2 - var->t1 + 1;
	*err = VAR_add_stats(var);
    }

    if (!*err && !lagtest) {
	*err = gretl_VAR_do_error_decomp(var->S, var->C);
    }

    if (*err) {
	gretl_VAR_free(var);
	var = NULL;
    }

#if VAR_DEBUG
    if (!*err) {
	gretl_matrix_print(var->A, "var->A");
    }
#endif

    return var;
}

static int *
maybe_expand_VAR_list (const int *list, double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, int *err)
{
    int needsep;
    int addconst = 0;
    int addtrend = 0;
    int addseas = 0;
    int *vlist = NULL;
    int vt = 0, di0 = 0;
    int i, l0;

    if (!(opt & OPT_N)) {
	addconst = 1;
    }

    if (opt & OPT_T) {
	addtrend = 1;
    }    

    if (pdinfo->pd > 1 && (opt & OPT_D)) {
	addseas = 1;
    }

    if (!addconst && !addtrend && !addseas) {
	return NULL;
    }

    if (addtrend) {
	vt = gettrend(pZ, pdinfo, 0);
	if (vt == 0) {
	    *err = E_ALLOC;
	    return NULL;
	}
    }	

    if (addseas) {
	di0 = dummy(pZ, pdinfo, 0);
	if (di0 == 0) {
	    *err = E_ALLOC;
	    return NULL;
	}
    } 

    needsep = !gretl_list_has_separator(list);

    l0 = list[0] + needsep + addconst + addtrend;

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
	if (addtrend) {
	    vlist[j++] = vt;
	}
	if (addconst) {
	    vlist[j++] = 0;
	}
    }

    return vlist;
}

static void var_lists_null (struct var_lists *vl)
{
    vl->detvars = NULL;
    vl->stochvars = NULL;
    vl->reglist = NULL;
    vl->testlist = NULL;
    vl->lagvlist = NULL;
}

/**
 * gretl_VAR:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: if includes %OPT_R, use robust VCV;
 *       if includes %OPT_I, print impulse responses;
 *       if includes %OPT_F, print forecast variance decompositions;
 *       if includes %OPT_D, add seasonal dummies;
 *       if includes %OPT_N, do not include a constant.
 *       if includes %OPT_Q, do not show individual regressions.
 *       if includes %OPT_T, include a linear trend.
 *       if includes %OPT_L, test for optimal lag length (only).
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Estimate a vector auto-regression (VAR), print and save
 * the results.
 *
 * Returns: pointer to VAR struct, which may be %NULL on error.
 */

GRETL_VAR *gretl_VAR (int order, int *list, double ***pZ, DATAINFO *pdinfo,
		      gretlopt opt, PRN *prn, int *err)
{
    GRETL_VAR *var = NULL;
    struct var_lists vl;
    int *vlist = NULL;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;

    var_lists_null(&vl);

    gretl_list_purge_const(list, (const double **) *pZ, pdinfo);

    vlist = maybe_expand_VAR_list(list, pZ, pdinfo, opt, err);

    if (!*err) {
	var = real_var(order, (vlist != NULL)? vlist : list, 
		       &vl, pZ, pdinfo, opt, err);
    }
    
    if (var != NULL) {
	if (opt & OPT_L) {
	    gretl_VAR_do_lagsel(var, &vl, pZ, pdinfo, prn);
	    gretl_VAR_free(var);
	    var = NULL;
	} else {
	    gretl_VAR_print(var, pdinfo, opt, prn);
	}
    }

    var_lists_free(&vl);

    if (vlist != NULL) {
	free(vlist);
    }

    /* reset original sample range */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2; 

    return var;
}

static int allocate_johansen_sigmas (JohansenInfo *jv)
{
    int k = gretl_matrix_rows(jv->u);
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

int gretl_VECM_test_beta (GRETL_VAR *vecm, PRN *prn)
{
    void *handle = NULL;
    int (*vecm_beta_test) (GRETL_VAR *, PRN *);
    int err = 0;

    if (vecm->jinfo == NULL || vecm->jinfo->D == NULL) {
	return E_DATA;
    }    

    gretl_error_clear();
    
    vecm_beta_test = get_plugin_function("vecm_beta_test", &handle);
    
    if (vecm_beta_test == NULL) {
	err = 1;
    } else {
	err = (* vecm_beta_test) (vecm, prn);
	close_plugin(handle);
    }
    
    return err;    
}

static int 
johansen_complete (GRETL_VAR *jvar, double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn)
{
    void *handle = NULL;
    int (*johansen) (GRETL_VAR *, double ***, DATAINFO *, gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();
    
    johansen = get_plugin_function("johansen_analysis", &handle);

    if (johansen == NULL) {
	err = 1;
    } else {
	err = (* johansen) (jvar, pZ, pdinfo, opt, prn);
	close_plugin(handle);
    }
    
    return err;
}

static int 
allocate_johansen_residual_matrices (GRETL_VAR *jvar)
{
    int T = jvar->t2 - jvar->t1 + 1;
    int vk = jvar->neqns;
    int err = 0;

    if (jvar->jinfo->u == NULL) {
	jvar->jinfo->u = gretl_matrix_alloc(jvar->neqns, T);
	if (jvar->jinfo->u == NULL) {
	    return E_ALLOC;
	}
    }

    if (restricted(jvar)) {
	vk++;
    }

    if (gretl_matrix_rows(jvar->jinfo->v) < vk) {
	gretl_matrix_free(jvar->jinfo->v);
	jvar->jinfo->v = gretl_matrix_alloc(vk, T);
	if (jvar->jinfo->v == NULL) {
	    gretl_matrix_free(jvar->jinfo->u);
	    jvar->jinfo->u = NULL;
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
    int i, k, err = 0;

    if (jv->difflist != NULL) {
	/* already done */
	return 0;
    }

    jv->difflist = gretl_list_new(neqns);
    if (jv->difflist == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	k = varlist[0] + jv->rank;
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

static int johansen_VAR (GRETL_VAR *jvar, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn)
{
    struct var_lists vlists;
    struct var_lists *vl = &vlists;
    MODEL jmod;
    int i, err = 0;

    err = organize_var_lists(jvar->jinfo->varlist, (const double **) *pZ, pdinfo, 
			     jvar->order, &vlists);
    if (err) {
	return err;
    }

    /* generate the required lags, if any */
    if (jvar->order > 0) {
	if ((err = VAR_list_laggenr(vl->stochvars, pZ, pdinfo, 
				    jvar->order, vl->lagvlist))) {
	    goto var_bailout;
	}
    }

    jvar->neqns = vl->stochvars[0];    

    /* compose base VAR list (entry 1 will vary across equations);
       assemble test list for setting t1 and t2 while we're at it 
    */
    err = compose_varlist(&vlists, vl->stochvars[1], 
			  jvar->order, 0, 1, pdinfo);
    if (err) {
	err = E_DATA;
	goto var_bailout;
    }

    if (jvar->t2 == 0) {
	/* sample hasn't been determined yet */
	if (check_for_missing_obs(vl->testlist, &pdinfo->t1, &pdinfo->t2,
				  (const double **) *pZ, NULL)) {
	    err = E_MISSDATA;
	    goto var_bailout;
	}

	jvar->t1 = pdinfo->t1;
	jvar->t2 = pdinfo->t2;
	jvar->T = jvar->t2 - jvar->t1 + 1;
    }

    err = allocate_johansen_residual_matrices(jvar); 
    if (err) {
	goto var_bailout;
    }

    if (jrank(jvar) > 0) {
	err = make_johansen_VECM_lists(jvar->jinfo, vl->reglist, 
				       jvar->neqns);
	if (err) {
	    goto var_bailout;
	}	    
    }

    gretl_model_init(&jmod);

    if (opt & OPT_V) {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), jvar->order);
    }

    for (i=0; i<jvar->neqns; i++) {

	/* insert the appropriate dependent variable */
	vl->reglist[1] = vl->stochvars[i + 1];

	if (jrank(jvar) > 0) {
	    /* VECM: record ID number of first difference.  Note:
	       we'll need this information in johansen.c, for building
	       the final VECM models, even if the order of the present
	       VAR is zero.
	    */
	    jvar->jinfo->difflist[i+1] = vl->reglist[1];
	}

	/* VAR in differences */
	if (vl->reglist[0] == 1) {
	    /* degenerate model (nothing to concentrate out) */
	    transcribe_data_as_uhat(vl->reglist[1], (const double **) *pZ,
				    jvar->jinfo->u, i, jvar->t1);
	} else {
	    jmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A | OPT_Z);
	    if ((err = jmod.errcode)) {
		fprintf(stderr, "*** johansen_VAR: VAR in differences, eqn %d, "
			"lsq err %d\n", i+1, err);
		printlist(vl->reglist, "list for this model");
		goto var_bailout;
	    }
	    if (opt & OPT_V) {
		jmod.aux = AUX_VAR;
		jmod.ID = i + 1;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    transcribe_uhat_to_matrix(&jmod, jvar->jinfo->u, i);
	    if (i == 0) {
		jvar->ifc = jmod.ifc;
	    }
	    clear_model(&jmod);
	}

	/* y_{t-1} regressions */
	vl->reglist[1] = jvar->jinfo->levels_list[i + 1];
	if (vl->reglist[0] == 1) {
	    /* degenerate */
	    transcribe_data_as_uhat(vl->reglist[1], (const double **) *pZ,
				    jvar->jinfo->v, i, jvar->t1);
	} else {
	    jmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A | OPT_Z);
	    if ((err = jmod.errcode)) {
		fprintf(stderr, "johansen_VAR: y_{t-1} regression, eqn %d, lsq err %d\n",
			i+1, err);
		goto var_bailout;
	    }	
	    if (opt & OPT_V) {
		jmod.aux = AUX_JOHANSEN;
		jmod.ID = -1;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    transcribe_uhat_to_matrix(&jmod, jvar->jinfo->v, i);
	    clear_model(&jmod);
	}
    }

    /* supplementary regressions for restricted cases */
    if (restricted(jvar)) {
	if (jcode(jvar) == J_REST_CONST) {
	    vl->reglist[1] = 0;
	} else {
	    vl->reglist[1] = gettrend(pZ, pdinfo, 0);
	    if (vl->reglist[1] == 0) {
		err = E_ALLOC;
		goto var_bailout;
	    }
	}

	if (vl->reglist[0] == 1) {
	    /* degenerate case */
	    transcribe_data_as_uhat(vl->reglist[1], (const double **) *pZ,
				    jvar->jinfo->v, i, jvar->t1);
	} else {	    
	    jmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A);
	    if ((err = jmod.errcode)) {
		fprintf(stderr, "johansen_VAR: restriction regression, "
			"eqn %d, lsq err %d\n", i+1, err);
		goto var_bailout;
	    }
	    if (opt & OPT_V) {
		jmod.aux = AUX_JOHANSEN;
		jmod.ID = -1;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    transcribe_uhat_to_matrix(&jmod, jvar->jinfo->v, i);
	    clear_model(&jmod);
	}
    }     

    pputc(prn, '\n');

 var_bailout:

    var_lists_free(&vlists);

#if VAR_DEBUG
    fprintf(stderr, "johansen_VAR: returning err = %d\n", err);
#endif    

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
johansen_info_new (const int *list, const int *exolist, int rank, gretlopt opt)
{
    JohansenInfo *jv = malloc(sizeof *jv);

    if (jv == NULL) {
	return NULL;
    }

    jv->levels_list = NULL;
    jv->varlist = NULL;

    jv->list = gretl_list_copy(list);
    if (jv->list == NULL) {
	free(jv);
	return NULL;
    }

    if (exolist != NULL) {
	jv->exolist = gretl_list_copy(exolist);
	if (jv->exolist == NULL) {
	    free(jv->list);
	    free(jv);
	    return NULL;
	}
    } else {
	jv->exolist = NULL;
    }

    jv->code = jcode_from_opt(opt);

    jv->ID = 0;

    jv->u = NULL;
    jv->v = NULL;

    jv->Suu = NULL;
    jv->Svv = NULL;
    jv->Suv = NULL;

    jv->Beta = NULL;
    jv->Alpha = NULL;
    jv->Bse = NULL;
    jv->Bvar = NULL;
    jv->D = NULL;

    jv->difflist = NULL;
    jv->biglist = NULL;

    jv->rank = rank;

    jv->seasonals = 0;
    jv->nexo = 0;

    return jv;
}

/* We don't simply use gretl_VAR_new(), because depending on the
   intended use of this VAR we may not need all the features
   added by that function.
*/

static GRETL_VAR *
johansen_VAR_new (const int *list, const int *exolist, int rank, int order, gretlopt opt)
{
    GRETL_VAR *var = malloc(sizeof *var);

    if (var == NULL) {
	return NULL;
    }

    var->jinfo = johansen_info_new(list, exolist, rank, opt);
    if (var->jinfo == NULL) {
	free(var);
	return NULL;
    }

    var->ci = VECM;
    var->refcount = 0;

    var->neqns = 0;
    var->order = order;
    var->ecm = 1;

    gretl_VAR_zero(var);

    return var;
}

/* Allocate storage, assemble basic lists and add basic required vars
   to the dataset, for Johansen analysis.  Note that the order in
   which variables are added (both to the dataset and to the lists) is
   significant, for the purposes of the bootstrap IRF analysis in
   irfboot.c: changing the order here will break things there.
*/

static GRETL_VAR *
johansen_VAR_prepare (int order, int rank, const int *list, const int *exolist, 
		      double ***pZ, DATAINFO *pdinfo, gretlopt opt)
{
    GRETL_VAR *jvar;
    int seasonals = (opt & OPT_D);
    int orig_t1 = pdinfo->t1;
    int orig_v = pdinfo->v;
    int nexo, di0 = 0, l0 = list[0];
    int i, k;

    jvar = johansen_VAR_new(list, exolist, rank, order - 1, opt);
    if (jvar == NULL) {
	return NULL;
    }    

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
	    int center = 1;

	    jvar->jinfo->seasonals = pdinfo->pd - 1;
	    nexo += pdinfo->pd - 1;
	    di0 = dummy(pZ, pdinfo, center);
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
    jvar->jinfo->levels_list = gretl_list_new(list[0]);
    if (jvar->jinfo->levels_list == NULL) {
	jvar->err = E_ALLOC;
	goto bailout;
    }

    /* full VAR list, including both endog and exog vars */
    jvar->jinfo->varlist = gretl_list_new(l0);
    if (jvar->jinfo->varlist == NULL) {
	jvar->err = E_ALLOC;
	goto bailout;
    }

    /* try to respect the chosen sample period: don't limit the
       generation of lags unnecessarily */
    pdinfo->t1 -= (order - 1);
    if (pdinfo->t1 < 0) {
	pdinfo->t1 = 0;
    }

    /* put first lags of endog vars into "levels_list" */
    for (i=1; i<=list[0]; i++) {
	jvar->jinfo->levels_list[i] = laggenr(list[i], 1, pZ, pdinfo);
	if (jvar->jinfo->levels_list[i] < 0) {
	    jvar->err = E_ALLOC;
	    goto bailout;
	}
    }

#if VAR_DEBUG
    printlist(jvar->jinfo->levels_list, "jvar->jinfo->levels_list (first lags)");
#endif

    /* put first differences of endog vars into VAR list */
    k = 1;
    for (i=1; i<=list[0]; i++) {
	jvar->jinfo->varlist[k] = diffgenr(list[i], DIFF, pZ, pdinfo);
	if (jvar->jinfo->varlist[k] < 0) {
	    jvar->err = E_ALLOC;
	    goto bailout;
	} 
	k++;
    }

#if VAR_DEBUG
    printlist(jvar->jinfo->varlist, "jvar->jinfo->varlist (first differences)");
#endif

    if (nexo > 0) {
	/* add separator before exog vars */
	jvar->jinfo->varlist[k++] = LISTSEP;
	jvar->jinfo->nexo = nexo;
    }

    if (exolist != NULL) {
	/* add specified exogenous variables to list */
	for (i=1; i<=exolist[0]; i++) {
	    jvar->jinfo->varlist[k++] = exolist[i];
	}
    }	    

    if (seasonals) {
	/* add seasonal dummies to list */
	for (i=0; i<pdinfo->pd-1; i++) {
	    jvar->jinfo->varlist[k++] = di0 + i;
	}
    }

    if (jcode(jvar) == J_UNREST_TREND) {
	/* add trend to VAR list */
	int vt = gettrend(pZ, pdinfo, 0);

	if (vt == 0) {
	    jvar->err = E_ALLOC;
	    goto bailout;
	} else {
	    jvar->jinfo->varlist[k++] = vt;
	}
    }	

    if (jcode(jvar) >= J_UNREST_CONST) {
	/* add the constant to the VAR list */
	jvar->jinfo->varlist[k++] = 0;
    }

 bailout:

    if (jvar->err) {
	dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
    }

    /* reset in case we moved this back */
    pdinfo->t1 = orig_t1;

    return jvar;
}

/* Driver function for Johansen analysis.  An appropriately
   initialized "jvar" must have been set up already, as in
   prepare_johansen_VAR().
*/

static int
johansen_driver (GRETL_VAR *jvar, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    PRN *varprn = (opt & OPT_V)? prn : NULL;
    
    /* estimate the equations in first differences, plus the
       equations in first lag of levels */
    jvar->err = johansen_VAR(jvar, pZ, pdinfo, opt, varprn); 

    if (jvar->jinfo->Suu == NULL && !jvar->err) {
	jvar->err = allocate_johansen_sigmas(jvar->jinfo);
    }

    if (!jvar->err) {
	gretl_matrix_multiply_mod(jvar->jinfo->u, GRETL_MOD_NONE,
				  jvar->jinfo->u, GRETL_MOD_TRANSPOSE,
				  jvar->jinfo->Suu, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(jvar->jinfo->v, GRETL_MOD_NONE,
				  jvar->jinfo->v, GRETL_MOD_TRANSPOSE,
				  jvar->jinfo->Svv, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(jvar->jinfo->u, GRETL_MOD_NONE,
				  jvar->jinfo->v, GRETL_MOD_TRANSPOSE,
				  jvar->jinfo->Suv, GRETL_MOD_NONE);

	gretl_matrix_divide_by_scalar(jvar->jinfo->Suu, jvar->T);
	gretl_matrix_divide_by_scalar(jvar->jinfo->Svv, jvar->T);
	gretl_matrix_divide_by_scalar(jvar->jinfo->Suv, jvar->T);

	if (jrank(jvar) == 0) {
	    char stobs[OBSLEN], endobs[OBSLEN];

	    pprintf(prn, "%s:\n", _("Johansen test"));
	    pprintf(prn, "%s = %d\n", _("Number of equations"), jvar->neqns);
	    pprintf(prn, "%s = %d\n", _("Lag order"), jvar->order + 1);
	    pprintf(prn, "%s: %s - %s (T = %d)\n", _("Estimation period"),
		    ntodate(stobs, jvar->t1, pdinfo), 
		    ntodate(endobs, jvar->t2, pdinfo), jvar->T);
	}

	if (opt & OPT_V) {
	    print_johansen_sigmas(jvar->jinfo, prn);
	}

	if (jrank(jvar) > 0) {
	    /* estimating VECM, not just doing cointegration test */
	    jvar->err = gretl_VAR_add_models(jvar);
	    if (!jvar->err) {
		jvar->err = VAR_add_companion_matrix(jvar);
	    }
	    if (!jvar->err) {
		jvar->err = VAR_add_cholesky_matrix(jvar);
	    }
	    if (!jvar->err) {
		jvar->err = VAR_add_residuals_matrix(jvar);
	    }	    
	}

	/* Now get the johansen plugin to finish the job, unless we're
	   bootstrapping IRFs (OPT_B), in which case that is handled
	   by the special code in irfboot.c
	*/
	if (!jvar->err && !(opt & OPT_B)) {
	    jvar->err = johansen_complete(jvar, pZ, pdinfo, opt, prn);
	}
    } 

    return jvar->err;
}

static GRETL_VAR *
johansen_wrapper (int order, int rank, const int *list, const int *exolist, 
		  double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    int oldv = pdinfo->v;

    jvar = johansen_VAR_prepare(order, rank, list, exolist, pZ, pdinfo, opt);

    if (jvar != NULL && !jvar->err) {
	jvar->err = johansen_driver(jvar, pZ, pdinfo, opt, prn);
    }

    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    if (jvar->err || !(opt & OPT_S)) {
	dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo);
    }

    return jvar;
}

/**
 * johansen_test:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: %OPT_A: include constant plus restricted trend; %OPT_D:
 * include centered seasonals; %OPT_N: no constant; %OPT_R:
 * restricted constant; %OPT_T: constant and unrestricted trend
 * (note: default "case" is unrestricted constant).
 * %OPT_V: produce verbose results.
 * @prn: gretl printing struct.
 *
 * Carries out the Johansen test for cointegration.
 *
 * Returns: pointer to struct containing information on 
 * the test.
 */

GRETL_VAR *johansen_test (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn)
{
    return johansen_wrapper(order, 0, list, NULL, pZ, pdinfo, opt, prn);
}

/**
 * johansen_test_simple:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: %OPT_A: include constant plus restricted trend; %OPT_D:
 * include centered seasonals; %OPT_N: no constant; %OPT_R:
 * restricted constant; %OPT_T: constant and unrestricted trend
 * (note: default "case" is unrestricted constant);
 * %OPT_V: produce verbose results; %OPT_Q: just print the tests.
 * @prn: gretl printing struct.
 *
 * Carries out the Johansen test for cointegration and prints the
 * results (but unlike johansen_test(), does not return the
 * allocated results in VAR form).
 *
 * Returns: 0 on success, non-zero code on error.
 */

int johansen_test_simple (int order, const int *list, double ***pZ, 
			  DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar;
    int err;

    jvar = johansen_wrapper(order, 0, list, NULL, pZ, pdinfo, opt, prn);
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
 * @err: location to receive error code.
 *
 * Returns: pointer to struct containing information on 
 * the VECM system or %NULL on failure.
 */

GRETL_VAR *gretl_VECM (int order, int rank, int *list, 
		       double ***pZ, DATAINFO *pdinfo,
		       gretlopt opt, PRN *prn, int *err)
{
    GRETL_VAR *jvar = NULL;
    int *endo_list = NULL, *exo_list = NULL;
    const int *vecm_list = list;

    gretl_list_purge_const(list, (const double **) *pZ, pdinfo);

    if (gretl_list_has_separator(list)) {
	*err = gretl_list_split_on_separator(list, &endo_list, &exo_list);
	if (*err) {
	    return jvar;
	}
	vecm_list = endo_list;
    } 

    if (rank <= 0 || rank > list[0]) {
	sprintf(gretl_errmsg, _("vecm: rank %d is out of bounds"), rank);
	*err = E_DATA;
	return jvar;
    }

    jvar = johansen_wrapper(order, rank, vecm_list, exo_list,
			    pZ, pdinfo, opt | OPT_S, prn);
    
    if (jvar != NULL) {
	if (!jvar->err) {
	    gretl_VAR_print(jvar, pdinfo, opt, prn);
	} else {
	    *err = jvar->err;
	}
    } else {
	*err = 1;
    }

    free(endo_list);
    free(exo_list);

    return jvar;
}

int gretl_VAR_attach_restrictions (GRETL_VAR *var, gretl_matrix *D)
{
    if (var->jinfo == NULL) {
	return 1;
    }

    if (var->jinfo->D != NULL) {
	gretl_matrix_free(var->jinfo->D);
    }

    var->jinfo->D = D;

    return 0;
}

void gretl_VAR_set_name (GRETL_VAR *var, const char *name)
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

int gretl_VAR_do_irf (GRETL_VAR *var, const char *line,
		      const double **Z, const DATAINFO *pdinfo)
{
    int targ = -1, shock = 1;
    int h = 0, boot = 0;
    int err = 0;
    char *s;

    s = strstr(line, "--targ=");
    if (s != NULL) {
	targ = atoi(s + 7) - 1;
    }

    s = strstr(line, "--shock=");
    if (s != NULL) {
	shock = atoi(s + 8) - 1;
    }

    s = strstr(line, "--horizon=");
    if (s != NULL) {
	h = atoi(s + 10);
    } else {
	h = 20;
    }

    if (strstr(line, "--bootstrap") != NULL) {
	boot = 1;
    }

#if 0
    fprintf(stderr, "targ=%d, shock=%d, h=%d, boot=%d\n", 
	    targ, shock, h, boot);
#endif

    if (targ >= 0 && shock >= 0 && h > 0) {
        if (boot) {
	    err = gretl_VAR_plot_impulse_response(var, targ, shock, h,
						  Z, pdinfo);
	} else {
	    err = gretl_VAR_plot_impulse_response(var, targ, shock, h, 
						  NULL, pdinfo);
	}
    }

    return err;
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

int gretl_VECM_n_beta (const GRETL_VAR *vecm)
{
    int nb = 0;

    if (vecm->jinfo != NULL && vecm->jinfo->Beta != NULL) {
	nb = gretl_matrix_rows(vecm->jinfo->Beta);
    }

    return nb;
}

int gretl_VECM_rank (const GRETL_VAR *vecm)
{
    int r = 0;

    if (vecm->jinfo != NULL) {
	r = vecm->jinfo->rank;
    }

    return r;
}

const int *gretl_VECM_list (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL) {
	return vecm->jinfo->list;
    }

    return NULL;
}

double *gretl_VAR_get_series (const GRETL_VAR *var, const DATAINFO *pdinfo, 
			      int idx, const char *key, int *err)
{
    double *x = NULL;
    const char *msel;
    int t, col = 0;

    if (var == NULL || idx != M_UHAT) { /* FIXME generalize this */
	*err = E_BADSTAT;
	return NULL;
    }

    msel = strchr(key, '[');
    if (msel == NULL || sscanf(msel, "[,%d]", &col) != 1) {
	*err = E_PARSE;
    } else if (col <= 0 || col > var->neqns) {
	*err = E_DATA;
    }

    if (!*err) {
	x = malloc(pdinfo->n * sizeof *x);
	if (x == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	const MODEL *pmod = var->models[col-1];

	if (pmod == NULL || pmod->full_n != pdinfo->n) {
	    *err = E_DATA;
	    free(x);
	    x = NULL;
	} else {
	    for (t=0; t<pdinfo->n; t++) {
		x[t] = pmod->uhat[t];
	    }
	}
    }

    return x;    
}

#define vecm_matrix(i) (i == M_JALPHA || i == M_JBETA || \
                        i == M_JVBETA || i == M_JS00 || \
                        i == M_JS11 || i == M_JS01)

gretl_matrix *gretl_VAR_get_matrix (const GRETL_VAR *var, int idx, 
				    int *err)
{
    const gretl_matrix *src = NULL;
    gretl_matrix *M = NULL;

    if (var == NULL) {
	*err = E_BADSTAT;
	return NULL;
    }

    if (idx == M_UHAT) {
	src = var->E;
    } else if (idx == M_COEFF) {
	src = var->A;
    } else if (idx == M_VCV) {
	src = var->S;
    } else if (vecm_matrix(idx)) {
	if (var->jinfo != NULL) {
	    switch (idx) {
	    case M_JALPHA: 
		src = var->jinfo->Alpha;
		break;
	    case M_JBETA: 
		src = var->jinfo->Beta;
		break;
	    case M_JVBETA: 
		src = var->jinfo->Bvar;
		break;
	    case M_JS00:
		src = var->jinfo->Suu;
		break;
	    case M_JS11:
		src = var->jinfo->Svv;
		break;
	    case M_JS01:
		src = var->jinfo->Suv;
		break;
	    }
	}
    }

    if (src == NULL) {
	*err = E_BADSTAT;
    } else {
	M = gretl_matrix_copy(src);
	if (M == NULL) {
	    *err = E_ALLOC;
	}
    }

    return M;
}

static GRETL_VAR *gretl_VAR_rebuilder_new (void)
{
    GRETL_VAR *var = malloc(sizeof *var);

    if (var == NULL) {
	return NULL;
    }

    var->ci = 0;
    var->refcount = 0;
    var->neqns = 0;
    var->order = 0;
    var->ecm = 0;

    gretl_VAR_zero(var);

    var->jinfo = NULL;

    return var;
}

static int VAR_retrieve_jinfo (xmlNodePtr node, xmlDocPtr doc,
			       GRETL_VAR *var)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    JohansenInfo *jinfo;
    int *list = NULL;
    int ID, code, rank;
    int seas, nexo;
    int got = 0;
    int err = 0;

    got += gretl_xml_get_prop_as_int(node, "ID", &ID);
    got += gretl_xml_get_prop_as_int(node, "code", &code);
    got += gretl_xml_get_prop_as_int(node, "rank", &rank);
    got += gretl_xml_get_prop_as_int(node, "seasonals", &seas);
    got += gretl_xml_get_prop_as_int(node, "nexo", &nexo);

    if (got != 5) {
	return E_DATA;
    }

    got = 0;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "list")) {
	    list = gretl_xml_node_get_list(cur, doc, &err);
	    got = 1;
	    break;
	} 
	cur = cur->next;
    }

    if (!err && !got) {
	err = E_DATA;
    }

    if (err) {
	return err;
    }

    jinfo = johansen_info_new(list, NULL, rank, OPT_NONE);
    if (jinfo == NULL) {
	return E_ALLOC;
    }

    jinfo->ID = ID;
    jinfo->code = code;
    jinfo->rank = rank;
    jinfo->seasonals = seas;
    jinfo->nexo = nexo;

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "difflist")) {
	    jinfo->difflist = gretl_xml_node_get_list(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "exolist")) {
	    jinfo->exolist = gretl_xml_node_get_list(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "biglist")) {
	    jinfo->biglist = gretl_xml_node_get_list(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "levels_list")) {
	    jinfo->levels_list = gretl_xml_node_get_list(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "varlist")) {
	    jinfo->varlist = gretl_xml_node_get_list(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "u")) {
	    jinfo->u = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "v")) {
	    jinfo->v = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Suu")) {
	    jinfo->Suu = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Suv")) {
	    jinfo->Suv = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Beta")) {
	    jinfo->Beta = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Alpha")) {
	    jinfo->Alpha = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Bse")) {
	    jinfo->Bse = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "D")) {
	    jinfo->D = gretl_xml_get_matrix(cur, doc, &err);
	}
	cur = cur->next;
    }

    if (err) {
	johansen_info_free(jinfo);
    } else {
	var->jinfo = jinfo;
    }

    fprintf(stderr, "VAR_retrieve_jinfo: err = %d\n", err);

    return err;
}

static int VAR_retrieve_equations (xmlNodePtr node, xmlDocPtr doc,
				   MODEL **models, int neqns)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    int i = 0, err = 0;

    while (cur != NULL && !err) {
	MODEL *pmod;

	if (!xmlStrcmp(cur->name, (XUC) "gretl-model")) {
	    pmod = gretl_model_from_XML(cur, doc, &err);
	    if (pmod != NULL) {
		models[i++] = pmod;
	    }
	}
	cur = cur->next;
    }

    fprintf(stderr, "VAR_retrieve_equations: got %d models (neqns = %d), err = %d\n", 
	    i, neqns, err);

    return err;
}

GRETL_VAR *gretl_VAR_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    GRETL_VAR *var;
    MODEL *pmod;
    xmlNodePtr cur;
    int start, rowmax;
    int i, j, k;
    int n, got = 0;

    var = gretl_VAR_rebuilder_new();
    if (var == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    got += gretl_xml_get_prop_as_int(node, "ci", &var->ci);
    got += gretl_xml_get_prop_as_int(node, "neqns", &var->neqns);
    got += gretl_xml_get_prop_as_int(node, "order", &var->order);
    got += gretl_xml_get_prop_as_int(node, "ecm", &var->ecm);

    if (got < 4) {
	*err = E_DATA;
	goto bailout;
    } 

    var->models = malloc(var->neqns * sizeof *var->models);
    if (var->models == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }
    
    for (i=0; i<var->neqns; i++) {
	var->models[i] = NULL;
    }    

    gretl_push_c_numeric_locale();

    cur = node->xmlChildrenNode;

    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "Fvals")) {
	    var->Fvals = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Ivals")) {
	    var->Ivals = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "equations")) {
	    *err = VAR_retrieve_equations(cur, doc, var->models, var->neqns);
	} else if (!xmlStrcmp(cur->name, (XUC) "gretl-johansen")) {
	    *err = VAR_retrieve_jinfo(cur, doc, var);
	}
	cur = cur->next;
    } 

    gretl_pop_c_numeric_locale();

    if (*err) {
	goto bailout;
    }

    pmod = var->models[0];

    var->ncoeff = pmod->ncoeff;
    var->ifc = pmod->ifc;
    var->t1 = pmod->t1;
    var->t2 = pmod->t2;
    var->T = var->t2 - var->t1 + 1;

    start = pmod->ifc;
    rowmax = var->neqns * var->order + start;

    /* set up storage for residuals */
    *err = VAR_add_residuals_matrix(var);

    /* set up storage for coefficients */
    if (!*err) {
	*err = VAR_add_companion_matrix(var);
    }    

    for (k=0; k<var->neqns && !*err; k++) {
	int v = 0, lag = 0;

	pmod = var->models[k];

	/* store residuals in var->E */
	for (i=0; i<var->T; i++) {
	    gretl_matrix_set(var->E, i, k, pmod->uhat[pmod->t1 + i]);
	}

	/* store coefficients in var->A */
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

    if (!*err) {
	/* covariance matrix and related things */
	*err = VAR_add_stats(var);
    }

 bailout:

    if (*err) {
	gretl_VAR_free(var);
	var = NULL;
    } 

    return var;
}

static void johansen_serialize (JohansenInfo *jinfo, FILE *fp)
{
    fprintf(fp, "<gretl-johansen ID=\"%d\" code=\"%d\" rank=\"%d\" ", 
	    jinfo->ID, jinfo->code, jinfo->rank);
    fprintf(fp, "seasonals=\"%d\" nexo=\"%d\">\n", jinfo->seasonals,
	    jinfo->nexo);

    gretl_xml_put_tagged_list("list", jinfo->list, fp);
    gretl_xml_put_tagged_list("difflist", jinfo->difflist, fp);
    gretl_xml_put_tagged_list("biglist", jinfo->biglist, fp);
    gretl_xml_put_tagged_list("exolist", jinfo->exolist, fp);
    gretl_xml_put_tagged_list("levels_list", jinfo->levels_list, fp);
    gretl_xml_put_tagged_list("varlist", jinfo->varlist, fp);

    gretl_xml_put_matrix(jinfo->u, "u", fp);
    gretl_xml_put_matrix(jinfo->v, "v", fp);
    gretl_xml_put_matrix(jinfo->Suu, "Suu", fp);
    gretl_xml_put_matrix(jinfo->Suv, "Suv", fp);
    gretl_xml_put_matrix(jinfo->Beta, "Beta", fp);
    gretl_xml_put_matrix(jinfo->Alpha, "Alpha", fp);
    gretl_xml_put_matrix(jinfo->Bse, "Bse", fp);
    gretl_xml_put_matrix(jinfo->D, "D", fp);

    fputs("</gretl-johansen>\n", fp);
}

int gretl_VAR_serialize (const GRETL_VAR *var, SavedObjectFlags flags,
			 FILE *fp)
{
    int g = var->neqns;
    int m = g * g + g;
    int i, err = 0;

    fprintf(fp, "<gretl-VAR name=\"%s\" saveflags=\"%d\" ", 
	    (var->name == NULL)? "none" : var->name, (int) flags);

    fprintf(fp, "ci=\"%d\" neqns=\"%d\" order=\"%d\" ecm=\"%d\">\n",
	    var->ci, var->neqns, var->order, var->ecm);

    gretl_push_c_numeric_locale();

    if (var->Fvals != NULL) {
	gretl_xml_put_double_array("Fvals", var->Fvals, m, fp);
    }
    if (var->Ivals != NULL) {
	gretl_xml_put_double_array("Ivals", var->Ivals, N_IVALS, fp);
    }

    gretl_pop_c_numeric_locale();

    fputs("<equations>\n", fp);

    for (i=0; i<var->neqns; i++) {
	gretl_model_serialize(var->models[i], 0, fp);
    }

    fputs("</equations>\n", fp);

    if (var->jinfo != NULL) {
	johansen_serialize(var->jinfo, fp);
    }

    fputs("</gretl-VAR>\n", fp);

    return err;
}

#include "irfboot.c"
#include "varomit.c"
