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
#include "vartest.h"
#include "libset.h"
#include "transforms.h"
#include "gretl_xml.h"

#define VDEBUG 0

enum {
    VAR_ESTIMATE,
    VAR_LAGSEL,
    VECM_ESTIMATE,
    VECM_CTEST
};

static JohansenInfo *
johansen_info_new (GRETL_VAR *var, int rank, gretlopt opt);

static int make_VAR_global_lists (GRETL_VAR *var);

static gretlopt opt_from_jcode (JohansenCode jc);

static int VAR_add_models (GRETL_VAR *var)
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
		gretl_matrix_set(var->A, i, j, (j == i - var->neqns)? 1 : 0);
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

void gretl_VAR_clear (GRETL_VAR *var)
{
    var->ci = 0;
    var->refcount = 0;
    var->err = 0;
    var->neqns = var->order = 0;
    var->t1 = var->t2 = var->T = 0;
    var->ifc = var->ncoeff = 0;
    var->ecm = var->detflags = 0;
    var->robust = var->qr = 0;

    var->ylist = NULL;
    var->xlist = NULL;

    var->Y = NULL;
    var->X = NULL;
    var->B = NULL;
    var->XTX = NULL;
    var->A = NULL;
    var->L = NULL;
    var->E = NULL;
    var->C = NULL;
    var->S = NULL;
    var->F = NULL;

    var->models = NULL;
    var->Fvals = NULL;
    var->Ivals = NULL;

    var->ll = var->ldet = var->LR = NADBL;
    var->AIC = var->BIC = var->HQC = NADBL;

    var->jinfo = NULL;
    var->name = NULL;
}

/* Construct the common X matrix (composed of lags of the core
   variables plus other terms if applicable)
*/

void VAR_fill_X (GRETL_VAR *v, int p, const double **Z, 
		 const DATAINFO *pdinfo)
{
    int diff = (v->ci == VECM);
    int i, j, s, t, vi;
    int k = 0; /* X column */

    /* const first */

    if (v->detflags & DET_CONST) {
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    gretl_matrix_set(v->X, s++, k, 1.0);
	}
	k++;
    }    

    /* add lagged Ys */

    for (i=0; i<v->neqns; i++) {
	vi = v->ylist[i+1];
	for (j=1; j<=p; j++) {
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		if (diff) {
		    gretl_matrix_set(v->X, s++, k, Z[vi][t-j] - Z[vi][t-j-1]);
		} else {
		    gretl_matrix_set(v->X, s++, k, Z[vi][t-j]);
		}
	    }
	    k++;
	}
    }

    /* add any exogenous vars */

    if (v->xlist != NULL) {
	for (i=1; i<=v->xlist[0]; i++) {
	    vi = v->xlist[i];
	    s = 0;
	    for (t=v->t1; t<=v->t2; t++) {
		gretl_matrix_set(v->X, s++, k, Z[vi][t]);
	    }
	    k++;
	}
    }

    /* add other deterministics */

    if (v->detflags & DET_SEAS) {
	int per = get_subperiod(v->t1, pdinfo, NULL);
	int pd1 = pdinfo->pd - 1;
	double s0, s1;

	s1 = (v->ci == VECM)? 1.0 - 1.0 / pdinfo->pd : 1.0;
	s0 = (v->ci == VECM)? s1 - 1.0 : 0.0;
	
	for (t=0; t<v->T; t++) {
	    for (i=0; i<pd1; i++) {
		gretl_matrix_set(v->X, t, k+i, (per == i)? s1 : s0);
	    }
	    per = (per < pd1)? per + 1 : 0;
	}
	k += pd1;
    }

    if (v->detflags & DET_TREND) {
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    gretl_matrix_set(v->X, s++, k, (double) (t + 1));
	}
	k++;
    }

#if VDEBUG
    gretl_matrix_print(v->X, "X");
#endif
}

/* construct the matrix of dependent variables for VAR or
   VECM estimation */

static void VAR_fill_Y (GRETL_VAR *v, int mod, const double **Z)
{
    int i, vi, s, t;

    for (i=0; i<v->neqns; i++) {
	vi = v->ylist[i+1];
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    if (mod == DIFF) {
		gretl_matrix_set(v->Y, s++, i, Z[vi][t] - Z[vi][t-1]);
	    } else if (mod == LAGS) {
		gretl_matrix_set(v->Y, s++, i, Z[vi][t-1]);
	    } else {
		gretl_matrix_set(v->Y, s++, i, Z[vi][t]);
	    }
	}
    }

    if (mod == LAGS && restricted(v)) {
	int trend = (v->jinfo->code == J_REST_TREND);
	int col = v->neqns;

	for (t=0; t<v->T; t++) {
	    gretl_matrix_set(v->Y, t, col, (trend)? v->t1 + t : 1);
	}
    }

#if VDEBUG
    gretl_matrix_print(v->Y, "Y");
#endif
}

/* split the user-supplied list, if need be, and construct
   the lists of endogenous and (possibly) exogenous vars.
   Deterministic terms are handled separately, via option
   flags.
*/

static int VAR_make_lists (GRETL_VAR *v, const int *list,
			   const double **Z,
			   const DATAINFO *pdinfo)
{
    int err = 0;

    if (gretl_list_has_separator(list)) {
	err = gretl_list_split_on_separator(list, &v->ylist, &v->xlist);
    } else {
	v->ylist = gretl_list_copy(list);
	if (v->ylist == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	gretl_list_purge_const(v->ylist, Z, pdinfo);
	if (v->xlist != NULL) {
	    gretl_list_purge_const(v->xlist, Z, pdinfo);
	}
    }

#if VDEBUG
    printlist(v->ylist, "v->ylist");
    printlist(v->xlist, "v->xlist");
#endif

    return err;
}

/* Starting from the given sample range, construct the feasible
   estimation range for a VAR or VECM.  Flag an error if there
   are missing values within the sample period.
*/

static int VAR_set_sample (GRETL_VAR *v, const double **Z, 
			   const DATAINFO *pdinfo)
{
    int diff = (v->ci == VECM)? 1 : 0;
    int i, vi, t, err = 0;

    /* advance t1 if needed */

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	int miss = 0, s = t - (v->order + diff);

	for (i=1; i<=v->ylist[0] && !miss; i++) {
	    vi = v->ylist[i];
	    if (na(Z[vi][t]) || s < 0 || na(Z[vi][s])) {
		v->t1 += 1;
		miss = 1;
	    }
	}
	if (v->xlist != NULL && !miss) {
	    for (i=1; i<=v->xlist[0] && !miss; i++) {
		if (na(Z[v->xlist[i]][t])) {
		    v->t1 += 1;
		    miss = 1;
		}
	    }
	}
	if (!miss) {
	    break;
	}
    }

    /* retard t2 if needed */

    for (t=pdinfo->t2; t>=v->t1; t--) {
	int miss = 0;

	for (i=1; i<=v->ylist[0] && !miss; i++) {
	    if (na(Z[v->ylist[i]][t])) {
		v->t2 -= 1;
		miss = 1;
	    }
	}
	if (v->xlist != NULL && !miss) {
	    for (i=1; i<=v->xlist[0] && !miss; i++) {
		if (na(Z[v->xlist[i]][t])) {
		    v->t2 -= 1;
		    miss = 1;
		}
	    }
	}
	if (!miss) {
	    break;
	}
    }

    /* reject sample in case of internal missing values */

    for (t=v->t1; t<=v->t2 && !err; t++) {
	for (i=1; i<=v->ylist[0] && !err; i++) {
	    if (na(Z[v->ylist[i]][t])) {
		err = E_MISSDATA;
	    }
	}
	if (v->xlist != NULL && !err) {
	    for (i=1; i<=v->xlist[0] && !err; i++) {
		if (na(Z[v->xlist[i]][t])) {
		    err = E_MISSDATA;
		}
	    }
	}
    }

    return err;
}

/* Account for deterministic terms and check for
   non-negative degrees of freedom. 
*/

static int VAR_check_df_etc (GRETL_VAR *v, const DATAINFO *pdinfo,
			     gretlopt opt)
{
    int err = 0;

    v->T = v->t2 - v->t1 + 1;
    v->ncoeff = v->order * v->neqns;

    if (v->xlist != NULL) {
	/* user-specified exogenous vars */
	v->ncoeff += v->xlist[0];
    }

    if (v->ci == VECM) {
	if (!(opt & OPT_N) && !(opt & OPT_R)) {
	    v->detflags |= DET_CONST;
	    v->ncoeff += 1;
	    v->ifc = 1;
	}	    
    } else if (!(opt & OPT_N)) {
	v->detflags |= DET_CONST;
	v->ncoeff += 1;
	v->ifc = 1;
    }

    if ((opt & OPT_D) && pdinfo->pd != 1) {
	v->detflags |= DET_SEAS;
	v->ncoeff += pdinfo->pd - 1;
    }

    if (opt & OPT_T) {
	v->detflags |= DET_TREND;
	v->ncoeff += 1;
    }

    if (v->T < v->ncoeff) {
	err = E_DF;
    }

    return err;
}

/* Note: if "v" is a VECM, and it includes a restricted const
   or trend, we make the Y and B matrices big enough to
   accommodate the additional restricted term.
*/

static int VAR_add_basic_matrices (GRETL_VAR *v, gretlopt opt)
{
    int n = v->neqns;

    if (v->ci == VECM && (opt & (OPT_A | OPT_R))) {
	/* allow for restricted term */
	n++;
    }

    v->Y = gretl_matrix_alloc(v->T, n);
    v->E = gretl_matrix_alloc(v->T, v->neqns);
    if (v->Y == NULL || v->E == NULL) {
	return E_ALLOC;
    }  

    if (v->ncoeff > 0) {
	v->X = gretl_matrix_alloc(v->T, v->ncoeff);
	v->B = gretl_matrix_alloc(v->ncoeff, n);
	if (v->X == NULL || v->B == NULL) { 
	    return E_ALLOC;
	}  
    }

    return 0;
}

/* main function for constructing a new VAR struct, which
   may be used for estimating a VAR, a VECM, or various
   auxiliary tasks */

static GRETL_VAR *gretl_VAR_new (int code, int order, int rank,
				 const int *list,
				 const double **Z,
				 const DATAINFO *pdinfo,
				 gretlopt opt, int *errp)
{
    GRETL_VAR *var;
    int ci, err = 0;

    ci = (code >= VECM_ESTIMATE)? VECM : VAR;

    if ((ci == VAR && order < 1) || (ci == VECM && order < 0)) {
	sprintf(gretl_errmsg, "VAR: invalid lag order %d", order);
	*errp = E_DATA;
	return NULL;
    }

    var = malloc(sizeof *var);
    if (var == NULL) {
	*errp = E_ALLOC;
	return NULL;
    }

    gretl_VAR_clear(var);

    var->ci = ci;
    var->order = order;
    var->ecm = (ci == VECM);

    var->qr = get_use_qr();
    if (ci == VAR && (opt & OPT_R)) {
	var->robust = 1;
    }

    err = VAR_make_lists(var, list, Z, pdinfo);

    if (!err && var->ylist[0] < 2) {
	strcpy(gretl_errmsg, "VAR: needs at least two variables");
	err = E_DATA;
    }

    if (rank > var->ylist[0]) {
	sprintf(gretl_errmsg, _("vecm: rank %d is out of bounds"), rank);
	err = E_DATA;
    }

    if (!err) {
	var->neqns = var->ylist[0];
	var->t1 = pdinfo->t1;
	var->t2 = pdinfo->t2;
    }

    /* FIXME below: some of these allocations are (still) redundant
       for some uses of the VAR struct */

    if (!err) {
	err = VAR_set_sample(var, Z, pdinfo);
    }

    if (!err) {
	err = VAR_check_df_etc(var, pdinfo, opt);
    }

    if (!err) {
	err = VAR_add_basic_matrices(var, opt);
    }

    if (!err && var->ci == VAR) {
	err = VAR_add_companion_matrix(var);
    }

    if (!err && var->ci == VAR) {
	err = VAR_add_cholesky_matrix(var);
    }

    if (!err && code != VAR_LAGSEL) {
	err = VAR_add_models(var);
    }

    if (!err && code == VAR_ESTIMATE) {
	int m = var->neqns * var->neqns + var->neqns;
	
	var->Fvals = malloc(m * sizeof *var->Fvals);
	if (var->Fvals == NULL) {
	    err = 1;
	}
    }

    if (!err && code == VAR_ESTIMATE) {
	var->Ivals = malloc(N_IVALS * sizeof *var->Ivals);
	if (var->Ivals == NULL) {
	    err = 1;
	}
    }

    if (!err && var->ci == VECM) {
	var->jinfo = johansen_info_new(var, rank, opt);
	if (var->jinfo == NULL) {
	    err = E_ALLOC;
	} else if (var->detflags & DET_SEAS) {
	    var->jinfo->seasonals = pdinfo->pd - 1;
	}
    }	

    if (!err) {
	if (var->ci == VAR) {
	    /* this is deferred for a VECM */
	    VAR_fill_Y(var, 0, Z);
	}
	VAR_fill_X(var, order, Z, pdinfo);
    }

    if (err) {
	gretl_VAR_free(var);
	var = NULL;
    }

    *errp = err;

    return var;
}

static void johansen_info_free (JohansenInfo *jv)
{
    gretl_matrix_free(jv->R0);
    gretl_matrix_free(jv->R1);

    gretl_matrix_free(jv->S00);
    gretl_matrix_free(jv->S11);
    gretl_matrix_free(jv->S01);

    gretl_matrix_free(jv->Beta);
    gretl_matrix_free(jv->Alpha);
    gretl_matrix_free(jv->Bvar);
    gretl_matrix_free(jv->Bse);
    gretl_matrix_free(jv->R);
    gretl_matrix_free(jv->q);
    gretl_matrix_free(jv->Ra);
    gretl_matrix_free(jv->qa);

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

    free(var->ylist);
    free(var->xlist);

    gretl_matrix_free(var->Y);
    gretl_matrix_free(var->X);
    gretl_matrix_free(var->B);
    gretl_matrix_free(var->XTX);

    gretl_matrix_free(var->A);
    gretl_matrix_free(var->L);
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

static int VAR_add_fcast_variance (GRETL_VAR *var, gretl_matrix *F,
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
VAR_add_forecast (GRETL_VAR *var, int t0, int t1, int t2,
		  const double **Z, const DATAINFO *pdinfo, 
		  gretlopt opt)
{
    const MODEL *pmod;
    gretl_matrix *F;
    double fti, xti;
    int nf = t2 - t0 + 1;
    int staticfc = (opt & OPT_S);
    int i, j, k, s, t;
    int lag, vj, m = 0;
    int fcols;

    fcols = (staticfc)? var->neqns : 2 * var->neqns;

    /* rows = number of forecast periods; cols = 1 to hold forecast
       for each variable, plus 1 to hold variance for each variable
       if forecast is dynamic.
    */

    F = gretl_zero_matrix_new(nf, fcols);
    if (F == NULL) {
	return E_ALLOC;
    }

    if (var->detflags & DET_SEAS) {
	m = get_subperiod(t0, pdinfo, NULL);
    }

    for (t=t0, s=0; t<=t2; t++, s++) {
	int miss = 0;

	for (i=0; i<var->neqns; i++) {
	    pmod = var->models[i];
	    fti = 0.0;
	    k = 0;

	    /* constant, if present */
	    if (pmod->ifc) {
		fti += pmod->coeff[k++];
	    }

	    /* lags of stochastic vars */
	    for (j=0; j<var->neqns; j++) {
		vj = var->ylist[j+1];
		for (lag=1; lag<=var->order; lag++) {
		    if (t < t1 || staticfc || s - lag < 0) {
			/* pre-forecast value */
			if (t - lag < 0) {
			    xti = NADBL;
			} else {
			    xti = Z[vj][t-lag];
			}
		    } else {
			/* prior forecast value */
			xti = gretl_matrix_get(F, s - lag, j); /* ?? */
		    }
		    if (na(xti)) {
			miss = 1;
		    } else {
			fti += pmod->coeff[k] * xti;
		    }
		    k++;
		}
	    }

	    /* exogenous vars, if any */
	    if (!miss && var->xlist != NULL) {
		for (j=1; j<=var->xlist[0]; j++) {
		    vj = var->xlist[j];
		    xti = Z[vj][t];
		    if (na(xti)) {
			miss = 1;
		    } else {
			fti += pmod->coeff[k] * xti;
		    }
		    k++;
		}
	    }
	    
	    /* other deterministics, if present */
	    if (!miss) {
		if (var->detflags & DET_SEAS) {
		    if (m < pdinfo->pd - 1) {
			fti += pmod->coeff[k+m];
		    } 
		}
		if (var->detflags & DET_TREND) {
		    k = pmod->ncoeff - 1;
		    fti += pmod->coeff[k] * (t + 1);
		}		
	    }

	    if (miss) {
		fti = NADBL;
	    }

	    gretl_matrix_set(F, s, i, fti);
	}

	if (var->detflags & DET_SEAS) {
	    m = (m < pdinfo->pd - 1)? m + 1 : 0;
	}
    }

    /* now get variances, if not static */
    if (!staticfc) {
	VAR_add_fcast_variance(var, F, nf, t1 - t0);
    }

    gretl_matrix_set_t1(F, t0);
    gretl_matrix_set_t2(F, t2);

    var->F = F;

#if VDEBUG
    gretl_matrix_print(F, "var->F");
#endif

    return 0;
}

static int
VECM_add_forecast (GRETL_VAR *var, int t0, int t1, int t2,
		   const double **Z, DATAINFO *pdinfo, 
		   gretlopt opt)
{
    gretl_matrix *F = NULL;
    gretl_matrix *B = NULL;

    double s0 = 0, s1 = 1;
    int order = var->order + 1;
    int nexo = (var->xlist != NULL)? var->xlist[0] : 0;
    int nseas = var->jinfo->seasonals;
    int nf = t2 - t0 + 1;
    int staticfc = (opt & OPT_S);
    int i, j, k, vj, s, t;
    int fcols, m = 0;

    fcols = (staticfc)? var->neqns : 2 * var->neqns;

    F = gretl_zero_matrix_new(nf, fcols);
    if (F == NULL) {
	return E_ALLOC;
    }

    B = VAR_coeff_matrix_from_VECM(var);
    if (B == NULL) {
	gretl_matrix_free(F);
	return E_ALLOC;
    }

    if (nseas > 0) {
	m = get_subperiod(t0, pdinfo, NULL);
	s1 -= 1.0 / pdinfo->pd;
	s0 = s1 - 1;
    }

    for (t=t0, s=0; t<=t2; t++, s++) {

	for (i=0; i<var->neqns; i++) {
	    double bij, xtj, fti = 0.0;
	    int ft, col = 0;

	    /* unrestricted constant, if present */
	    if (var->ifc) {
		bij = gretl_matrix_get(B, i, col++);
		fti += bij;
	    }

	    /* lags of endogenous vars */
	    for (j=0; j<var->neqns; j++) {
		vj = var->ylist[j+1];
		for (k=1; k<=order; k++) {
		    if (t - k < 0) {
			fti = NADBL;
			break;
		    }			
		    bij = gretl_matrix_get(B, i, col++);
		    ft = s - k;
		    if (t >= t1 && ft >= 0 && !staticfc) {
			/* use prior forecast if available */
			xtj = gretl_matrix_get(F, ft, j);
		    } else {
			xtj = Z[vj][t-k];
		    }
		    if (na(xtj)) {
			fti = NADBL;
			break;
		    } else {
			fti += bij * xtj;
		    }
		}
		if (na(fti)) {
		    break;
		}
	    }

	    if (na(fti)) goto set_fcast;

	    /* exogenous vars, if present */
	    for (j=0; j<nexo; j++) {
		vj = var->xlist[j+1];
		xtj = Z[vj][t];
		if (na(xtj)) {
		    fti = NADBL;
		} else {
		    bij = gretl_matrix_get(B, i, col++);
		    fti += bij * xtj;
		}
	    }

	    if (na(fti)) goto set_fcast;

	    /* seasonals, if present */
	    for (j=0; j<nseas; j++) {
		xtj = (m == j)? s1 : s0;
		bij = gretl_matrix_get(B, i, col++);
		fti += bij * xtj;
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
		fti += bij * t; /* ?? */
	    }
	    
	set_fcast:

	    gretl_matrix_set(F, s, i, fti);
	}

	if (nseas > 0) {
	    m = (m < pdinfo->pd - 1)? m + 1 : 0;
	}
    }

    gretl_matrix_free(B);

    /* now get variances, if not static */
    if (!staticfc) {
	VAR_add_fcast_variance(var, F, nf, t1 - t0);
    }

    gretl_matrix_set_t1(F, t0);
    gretl_matrix_set_t2(F, t2);

    var->F = F;

#if VDEBUG
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
	    VECM_add_forecast(var, t0, t1, t2, Z, pdinfo, opt);
	} else {
	    VAR_add_forecast(var, t0, t1, t2, Z, pdinfo, opt);
	}
    }

    return var->F;
}

const gretl_matrix *
gretl_VAR_get_residual_matrix (const GRETL_VAR *var)
{
    return var->E;
}

static void VAR_dw_rho (MODEL *pmod)
{
    double ut, u1;
    double xx = 0;
    double ut1 = 0, u11 = 0;
    int t, s;

    for (t=pmod->t1; t<=pmod->t2; t++)  {
	s = t - 1;
	if (s >= 0 && !na(pmod->uhat[s])) {
	    ut = pmod->uhat[t];
	    u1 = pmod->uhat[s];
	    xx += (ut - u1) * (ut - u1);
	    ut1 += ut * u1;
	    u11 += u1 * u1;
	}
    }

    pmod->dw = xx / pmod->ess;
    pmod->rho = ut1 / u11;
}

/* set basic model statistics when estimation has been
   done by matrix methods */

int set_VAR_model_stats (GRETL_VAR *var, int i)
{
    MODEL *pmod = var->models[i];
    double *y = NULL;
    double u, x, SSR = 0, TSS = 0;
    int t;

    y = malloc(var->T * sizeof *y);
    if (y == NULL) {
	pmod->ybar = pmod->sdy = NADBL;
	pmod->rsq = NADBL;
	return E_ALLOC;
    }

    for (t=0; t<var->T; t++) {
	y[t] = gretl_matrix_get(var->Y, t, i);
    }

    pmod->ybar = gretl_mean(0, var->T - 1, y);
    pmod->sdy = gretl_stddev(0, var->T - 1, y);

    for (t=0; t<var->T; t++) {
	u = gretl_matrix_get(var->E, t, i);
	SSR += u * u;
	x = (var->ifc)? y[t] - pmod->ybar : y[t];
	TSS += x * x;
	pmod->uhat[t + pmod->t1] = u;
	pmod->yhat[t + pmod->t1] = y[t] - u;
    }

    pmod->ess = SSR;
    pmod->sigma = sqrt(SSR / pmod->dfd);
    pmod->tss = TSS;
    pmod->rsq = 1.0 - SSR / TSS;
    pmod->fstt = ((TSS - SSR) / pmod->dfn) / (SSR / pmod->dfd);

    pmod->adjrsq = pmod->lnL = NADBL;

    VAR_dw_rho(pmod);

    free(y);

    return 0;
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
    if (var->models != NULL && k >= 0 && k < var->neqns) {
	return var->models[k]->list[1];
    } else {
	return 0;
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

static int VAR_add_roots (GRETL_VAR *var)
{
    gretl_matrix *CompForm = NULL;
    int err = 0;

    if (var->A == NULL) {
	return 1;
    }

    var->L = NULL;

    CompForm = gretl_matrix_copy(var->A);
    if (CompForm == NULL) {
	err = E_ALLOC;
    }

    /* save eigenvalues of companion form matrix */
    if (!err) {
        var->L = gretl_general_matrix_eigenvals(CompForm, 0, &err);
#if 0
	gretl_matrix_print(var->A, "Companion form matrix");
	gretl_matrix_print(var->L, "Eigenvalues");
#endif
    }

    gretl_matrix_free(CompForm);

    if (err) {
	gretl_matrix_free(var->L);
	var->L = NULL;
    }

    return err;
}

const gretl_matrix *gretl_VAR_get_roots (GRETL_VAR *var)
{
    if (var->L == NULL) {
	/* roots not computed yet */
	VAR_add_roots(var);
    }

    return var->L;
}

double gretl_VAR_ldet (GRETL_VAR *var, int *err)
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

void gretl_VAR_param_names (GRETL_VAR *v, char **params, 
			    const DATAINFO *pdinfo)
{
    char lagstr[8];
    int i, j, n, k = 0;

    if (v->detflags & DET_CONST) {
	strcpy(params[k++], pdinfo->varname[0]);
    }     

    for (i=1; i<=v->ylist[0]; i++) {
	for (j=1; j<=v->order; j++) {
	    sprintf(lagstr, "_%d", j);
	    n = strlen(lagstr);
	    if (v->ci == VECM) {
		strcpy(params[k], "d_");
		n += 2;
	    }
	    strncat(params[k], pdinfo->varname[v->ylist[i]],
		    VNAMELEN - n - 1);
	    strncat(params[k], lagstr, n);
	    k++;
	}
    }

    if (v->xlist != NULL) {
	for (i=1; i<=v->xlist[0]; i++) {
	    strcpy(params[k++], pdinfo->varname[v->xlist[i]]);
	}
    }

    if (v->detflags & DET_SEAS) {
	for (i=1; i<pdinfo->pd; i++) {
	    sprintf(params[k++], "S%d", i);
	}	
    }

    if (v->detflags & DET_TREND) {
	strcpy(params[k++], "time");
    }

    if (v->ci == VECM) {
	int rank = jrank(v);

	for (i=0; i<rank; i++) {
	    sprintf(params[k++], "EC%d", i+1);
	}
    }
}

static int make_A_matrix (GRETL_VAR *var, int ifc)
{
    int i, j, v, lag;
    int dim = var->neqns * var->order;
    double bij;

    for (j=0; j<var->neqns; j++) {
	v = lag = 0;
	for (i=0; i<dim; i++) {
	    bij = gretl_matrix_get(var->B, i+ifc, j);
	    gretl_matrix_set(var->A, j, var->neqns * lag + v, bij);
	    if (lag < var->order - 1) {
		lag++;
	    } else {
		lag = 0;
		v++;
	    }
	}
    }

    return 0;
}

static int VAR_finalize (GRETL_VAR *var, int lagsel)
{
    int err = 0;

    if (var->detflags & DET_CONST) {
	var->ifc = 1;
    }

    if (!lagsel) {
	/* not needed if we're just doing lag selection */

	err = make_A_matrix(var, var->ifc);

	if (!err) {
	    err = VAR_wald_omit_tests(var, var->ifc);
	}

	if (!err && var->order > 1) {
	    err = last_lag_LR_prep(var, var->ifc);
	}
    }

    if (!err) {
	err = VAR_add_stats(var);
    }

    if (!err && !lagsel) {
	err = gretl_VAR_do_error_decomp(var->S, var->C);
    }

    return err;
}

/* transcribe the per-equation output from a VAR into
   MODEL structs, if wanted */

static int set_up_VAR_models (GRETL_VAR *var, 
			      const double **Z,
			      const DATAINFO *pdinfo)
{
    MODEL *pmod;
    char **params = NULL;
    int yno, N = pdinfo->n;
    const double *y;
    int i, j;
    int err = 0;

    params = strings_array_new_with_length(var->ncoeff, VNAMELEN);
    if (params == NULL) {
	return E_ALLOC;
    }

    gretl_VAR_param_names(var, params, pdinfo);

    for (i=0; i<var->neqns && !err; i++) {
	yno = var->ylist[i+1];
	y = Z[yno];

	pmod = var->models[i];
	pmod->ID = i + 1;
	pmod->ci = VAR;
	pmod->aux = AUX_VAR;

	pmod->full_n = N;
	pmod->nobs = var->T;
	pmod->t1 = var->t1;
	pmod->t2 = var->t2;
	pmod->ncoeff = var->ncoeff;
	pmod->dfd = pmod->nobs - pmod->ncoeff;
	pmod->ifc = (var->detflags & DET_CONST)? 1 : 0;
	pmod->dfn = var->ncoeff - pmod->ifc;

	err = gretl_model_allocate_storage(pmod);

	pmod->depvar = gretl_strdup(pdinfo->varname[yno]);

	if (i == 0) {
	    pmod->params = params;
	} else {
	    pmod->params = strings_array_dup(params, var->ncoeff);
	}
	pmod->nparams = var->ncoeff;

	pmod->list = gretl_list_new(1);
	pmod->list[1] = yno;

	set_VAR_model_stats(var, i);

	for (j=0; j<var->ncoeff; j++) {
	    pmod->coeff[j] = gretl_matrix_get(var->B, j, i);
	}
    }

    return err;
}

/**
 * gretl_VAR:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @Z: data array.
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
 * @errp: location to receive error code.
 *
 * Estimate a vector auto-regression (VAR), print and save
 * the results.
 *
 * Returns: pointer to VAR struct, which may be %NULL on error.
 */

GRETL_VAR *gretl_VAR (int order, int *list, 
		      const double **Z, const DATAINFO *pdinfo,
		      gretlopt opt, PRN *prn, 
		      int *errp)
{
    GRETL_VAR *var = NULL;
    int code = (opt & OPT_L)? VAR_LAGSEL : VAR_ESTIMATE;
    int err = 0;

    var = gretl_VAR_new(code, order, 0, list, Z, pdinfo, 
			opt, &err);
    if (var == NULL) {
	return NULL;
    }

    /* run the regressions: use QR or Cholesky */
    if (var->qr) {
	err = gretl_matrix_QR_ols(var->Y, var->X, 
				  var->B, var->E,
				  &var->XTX, NULL);
    } else {
	err = gretl_matrix_multi_ols(var->Y, var->X, 
				     var->B, var->E,
				     &var->XTX);
    }

    if (!err) {
	if (code == VAR_LAGSEL) {
	    err = VAR_finalize(var, 1);
	    if (!err) {
		err = VAR_do_lagsel(var, Z, pdinfo, prn);
	    }
	} else {
	    err = set_up_VAR_models(var, Z, pdinfo);
	    if (!err) {
		err = VAR_finalize(var, 0);
	    }
	    if (!err) {
		gretl_VAR_print(var, pdinfo, opt, prn);
	    }
	}
    }

    if (code == VAR_LAGSEL || (err && var != NULL)) {
	gretl_VAR_free(var);
	var = NULL;
    }

    *errp = err;

    return var;
}

static void 
print_johansen_sigmas (const JohansenInfo *jv, PRN *prn)
{
    int nr, nc;
    int i, j;

    pprintf(prn, "\n%s\n\n", _("Sample variance-covariance matrices for residuals"));

    nr = gretl_matrix_rows(jv->S00);
    pprintf(prn, " %s\n\n", _("VAR system in first differences"));
    for (i=0; i<nr; i++) {
	for (j=0; j<nr; j++) {
	    pprintf(prn, "%#12.5g", gretl_matrix_get(jv->S00, i, j));
	}
	pputc(prn, '\n');
    }

    nr = gretl_matrix_rows(jv->S11);
    pprintf(prn, "\n %s\n\n", _("System with levels as dependent variable"));
    for (i=0; i<nr; i++) {
	for (j=0; j<nr; j++) {
	    pprintf(prn, "%#12.5g", gretl_matrix_get(jv->S11, i, j));
	}
	pputc(prn, '\n');
    } 
    
    nr = gretl_matrix_rows(jv->S01);
    nc = gretl_matrix_cols(jv->S01);
    pprintf(prn, "\n %s\n\n", _("Cross-products"));
    for (i=0; i<nr; i++) {
	for (j=0; j<nc; j++) {
	    pprintf(prn, "%#12.5g", gretl_matrix_get(jv->S01, i, j));
	}
	pputc(prn, '\n');
    }     
}

int gretl_VECM_test (GRETL_VAR *vecm, 
		     const gretl_restriction_set *rset,
		     const DATAINFO *pdinfo, 
		     gretlopt opt,
		     PRN *prn)
{
    void *handle = NULL;
    int (*vecm_test_restriction) (GRETL_VAR *, 
				  const gretl_restriction_set *,
				  const DATAINFO *,
				  gretlopt opt,
				  PRN *);
    int err = 0;

    if (vecm->jinfo == NULL || rset == NULL) {
	return E_DATA;
    }    

    gretl_error_clear();

    vecm_test_restriction = 
	get_plugin_function("vecm_test_restriction", &handle);
    
    if (vecm_test_restriction == NULL) {
	err = 1;
    } else {
	err = (* vecm_test_restriction) (vecm, rset, pdinfo, 
					 opt, prn);
	close_plugin(handle);
    }

    return err;    
}

static int 
johansen_test_complete (GRETL_VAR *jvar, const DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn)
{
    void *handle = NULL;
    int (*johansen) (GRETL_VAR *, const DATAINFO *, gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();
    
    johansen = get_plugin_function("johansen_coint_test", &handle);

    if (johansen == NULL) {
	err = 1;
    } else {
	err = (* johansen) (jvar, pdinfo, opt, prn);
	close_plugin(handle);
    }
    
    return err;
}

static int 
johansen_estimate_complete (GRETL_VAR *jvar, 
			    const gretl_restriction_set *rset,
			    const double **Z, const DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn)
{
    void *handle = NULL;
    int (*johansen) (GRETL_VAR *, const gretl_restriction_set *,
		     const double **, const DATAINFO *, 
		     gretlopt, PRN *);
    int err = 0;

    gretl_error_clear();
    
    johansen = get_plugin_function("johansen_estimate", &handle);

    if (johansen == NULL) {
	err = 1;
    } else {
	err = (* johansen) (jvar, rset, Z, pdinfo, opt, prn);
	close_plugin(handle);
    }
    
    return err;
}

static int 
allocate_johansen_residual_matrices (GRETL_VAR *v)
{
    int p = v->neqns;
    int p1 = p + restricted(v);

    clear_gretl_matrix_err();

    /* N.B. allow for the possibility that this has already been done */

    if (v->jinfo->R0 == NULL) {
	v->jinfo->R0 = gretl_matrix_alloc(v->T, p);
	v->jinfo->R1 = gretl_matrix_alloc(v->T, p1);
    }

    if (v->jinfo->S00 == NULL) {
	v->jinfo->S00 = gretl_matrix_alloc(p, p);
	v->jinfo->S11 = gretl_matrix_alloc(p1, p1);
	v->jinfo->S01 = gretl_matrix_alloc(p, p1);
    }

    return get_gretl_matrix_err();
}

static void johansen_degenerate_stage_1 (GRETL_VAR *v, 
					 const double **Z)
{
    gretl_matrix *R0 = v->jinfo->R0;
    gretl_matrix *R1 = v->jinfo->R1;
    int i, vi, s, t;

    for (i=0; i<v->neqns; i++) {
	vi = v->ylist[i+1];
	s = 0;
	for (t=v->t1; t<=v->t2; t++) {
	    gretl_matrix_set(R0, s++, i, Z[vi][t] - Z[vi][t-1]);
	    gretl_matrix_set(R1, s++, i, Z[vi][t-1]);
	}
    }

    if (restricted(v)) {
	int trend = (jcode(v) == J_REST_TREND);

	for (t=0; t<v->T; t++) {
	    gretl_matrix_set(R1, t, v->neqns, (trend)? v->t1 + t : 1);
	}
    }	
}

/* For Johansen analysis: estimate VAR in differences along with the
   other auxiliary regressions required to compute the relevant
   matrices of residuals, for concentration of the log-likelihood.
   Then compute S00, S11, S01.
*/

int johansen_stage_1 (GRETL_VAR *jvar, 
		      const double **Z, const DATAINFO *pdinfo,
		      PRN *prn)
{
    int restr = restricted(jvar);
    int err;

    err = allocate_johansen_residual_matrices(jvar); 
    if (err) {
	return err;
    }

    if (jvar->ncoeff == 0) {
	/* nothing to concentrate out */
	johansen_degenerate_stage_1(jvar, Z);
	goto finish;
    }

    if (restr) {
	/* shrink by one column */
	gretl_matrix_reuse(jvar->Y, -1, jvar->neqns);
	gretl_matrix_reuse(jvar->B, -1, jvar->neqns);
    }

    VAR_fill_Y(jvar, DIFF, Z);

    /* (1) VAR in first differences */

    if (jvar->qr) {
	err = gretl_matrix_QR_ols(jvar->Y, jvar->X, jvar->B, 
				  jvar->jinfo->R0, NULL, NULL);
    } else {
	err = gretl_matrix_multi_ols(jvar->Y, jvar->X, jvar->B, 
				     jvar->jinfo->R0, NULL);
    }

    if (restr) {
	/* re-expand to full size */
	gretl_matrix_reuse(jvar->Y, -1, jvar->neqns + 1);
	gretl_matrix_reuse(jvar->B, -1, jvar->neqns + 1);
    }

    /* (2) System with lagged levels on LHS (may include
       restricted const or trend in last column)
    */

    VAR_fill_Y(jvar, LAGS, Z);

    if (jvar->qr) {
	err = gretl_matrix_QR_ols(jvar->Y, jvar->X, jvar->B, 
				  jvar->jinfo->R1, NULL, NULL);
    } else {
	err = gretl_matrix_multi_ols(jvar->Y, jvar->X, jvar->B, 
				     jvar->jinfo->R1, NULL);
    }  

    if (restr) {
	/* shrink to "final" size */
	gretl_matrix_reuse(jvar->Y, -1, jvar->neqns);
	gretl_matrix_reuse(jvar->B, -1, jvar->neqns);
    }

 finish:

    gretl_matrix_multiply_mod(jvar->jinfo->R0, GRETL_MOD_TRANSPOSE,
			      jvar->jinfo->R0, GRETL_MOD_NONE,
			      jvar->jinfo->S00, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(jvar->jinfo->R1, GRETL_MOD_TRANSPOSE,
			      jvar->jinfo->R1, GRETL_MOD_NONE,
			      jvar->jinfo->S11, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(jvar->jinfo->R0, GRETL_MOD_TRANSPOSE,
			      jvar->jinfo->R1, GRETL_MOD_NONE,
			      jvar->jinfo->S01, GRETL_MOD_NONE);

    gretl_matrix_divide_by_scalar(jvar->jinfo->S00, jvar->T);
    gretl_matrix_divide_by_scalar(jvar->jinfo->S11, jvar->T);
    gretl_matrix_divide_by_scalar(jvar->jinfo->S01, jvar->T);

#if VDEBUG
    fprintf(stderr, "johansen_stage_1: returning err = %d\n", err);
#endif    

    return err;
}

static gretlopt opt_from_jcode (JohansenCode jc)
{
    gretlopt opt = OPT_NONE;

    if (jc == J_NO_CONST) {
	opt = OPT_N;
    } else if (jc == J_UNREST_TREND) {
	opt = OPT_T;
    } else if (jc == J_REST_CONST) {
	opt =  OPT_R;
    } else if (jc == J_REST_TREND) {
	opt = OPT_A;
    }

    return opt;
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
johansen_info_new (GRETL_VAR *var, int rank, gretlopt opt)
{
    JohansenInfo *jv = malloc(sizeof *jv);

    if (jv == NULL) {
	return NULL;
    }

    jv->ID = 0;
    jv->code = jcode_from_opt(opt);
    jv->rank = rank;
    jv->seasonals = 0;

    jv->R0 = NULL;
    jv->R1 = NULL;
    jv->S00 = NULL;
    jv->S11 = NULL;
    jv->S01 = NULL;
    jv->Beta = NULL;
    jv->Alpha = NULL;
    jv->Bse = NULL;
    jv->Bvar = NULL;
    jv->R = NULL;
    jv->q = NULL;
    jv->Ra = NULL;
    jv->qa = NULL;

    jv->ll0 = NADBL;
    jv->bdf = 0;

    return jv;
}

/* Driver function for Johansen analysis.  An appropriately
   initialized "jvar" must have been set up already.
*/

static int
johansen_driver (GRETL_VAR *jvar, 
		 const gretl_restriction_set *rset,
		 const double **Z, const DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    PRN *varprn = (opt & OPT_V)? prn : NULL;
    
    jvar->err = johansen_stage_1(jvar, Z, pdinfo, varprn); 
    if (jvar->err) {
	return jvar->err;
    }

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
	jvar->err = VAR_add_models(jvar);
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

    /* Now get the johansen plugin to finish the job */
    if (!jvar->err) {
	if (jrank(jvar) > 0) {
	    jvar->err = johansen_estimate_complete(jvar, rset, Z, pdinfo, 
						   opt, prn);
	} else {
	    jvar->err = johansen_test_complete(jvar, pdinfo, opt, prn);
	}
    }

    return jvar->err;
}

static GRETL_VAR *
johansen_wrapper (int code, int order, int rank, 
		  const int *list, 
		  const gretl_restriction_set *rset, 
		  const double **Z, const DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar;
    int err = 0;

    jvar = gretl_VAR_new(code, order - 1, rank, list, Z, pdinfo,
			 opt, &err);
    if (jvar == NULL) {
	return NULL;
    }    

    if (jvar != NULL && !jvar->err) {
	jvar->err = johansen_driver(jvar, rset, Z, pdinfo, opt, prn);
    }

    return jvar;
}

/**
 * johansen_test:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @Z: data array.
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

GRETL_VAR *johansen_test (int order, const int *list, 
			  const double **Z, const DATAINFO *pdinfo,
			  gretlopt opt, PRN *prn)
{
    return johansen_wrapper(VECM_CTEST, order, 0, list, NULL, 
			    Z, pdinfo, opt, prn);
}

/**
 * johansen_test_simple:
 * @order: lag order for test.
 * @list: list of variables to test for cointegration.
 * @Z: data array.
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

int johansen_test_simple (int order, const int *list, 
			  const double **Z, const DATAINFO *pdinfo, 
			  gretlopt opt, PRN *prn)
{
    GRETL_VAR *jvar;
    int err;

    jvar = johansen_wrapper(VECM_CTEST, order, 0, list, NULL, 
			    Z, pdinfo, opt, prn);
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
 * gretl_VECM:
 * @order: lag order for test.
 * @rank: pre-specified cointegration rank.
 * @list: list of variables to test for cointegration.
 * @Z: pointer to data array.
 * @pdinfo: dataset information.
 * @opt:
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Returns: pointer to struct containing information on 
 * the VECM system or %NULL on failure.
 */

GRETL_VAR *gretl_VECM (int order, int rank, int *list, 
		       const double **Z, const DATAINFO *pdinfo,
		       gretlopt opt, PRN *prn, int *err)
{
    GRETL_VAR *jvar = NULL;

    if (rank <= 0) {
	sprintf(gretl_errmsg, _("vecm: rank %d is out of bounds"), rank);
	*err = E_DATA;
	return NULL;
    }

    jvar = johansen_wrapper(VECM_ESTIMATE, order, rank, list, NULL, 
			    Z, pdinfo, opt | OPT_S, prn);
    
    if (jvar != NULL) {
	if (!jvar->err) {
	    gretl_VAR_print(jvar, pdinfo, opt, prn);
	} else {
	    *err = jvar->err;
	}
    } else {
	*err = 1;
    }

    return jvar;
}

static int *rebuild_full_VAR_list (const GRETL_VAR *var, int *err)
{
    int *list = NULL;

    if (var->xlist == NULL) {
	list = gretl_list_copy(var->ylist);
    } else {
	int i, j, lsep = var->xlist[0] > 0;

	list = gretl_list_new(var->neqns + var->xlist[0] + lsep);
	if (list == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}

	j = 1;
	for (i=0; i<var->neqns; i++) {
	    list[j++] = var->ylist[i+1];
	}

	if (lsep) {
	    list[j++] = LISTSEP;
	}

	for (i=1; i<=var->xlist[0]; i++) {
	    list[j++] = var->xlist[i];
	}   
    } 

    return list;
}

/**
 * real_gretl_restricted_vecm:
 * @orig: orginal VECM model.
 * @rset: restriction information.
 * @Z: pointer to data array.
 * @pdinfo: dataset information.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Returns: pointer to struct containing information on 
 * the newly restricted VECM system or %NULL on failure.
 */

GRETL_VAR *
real_gretl_restricted_vecm (GRETL_VAR *orig, 
			    const gretl_restriction_set *rset,
			    const double **Z, const DATAINFO *pdinfo, 
			    PRN *prn, int *err)
{
    GRETL_VAR *jvar = NULL;
    gretlopt opt = OPT_S;
    int *list = NULL;

    if (orig == NULL || orig->jinfo == NULL || rset == NULL) {
	*err = E_DATA;
	return NULL;
    }   

    list = rebuild_full_VAR_list(orig, err);
    if (*err) {
	return NULL;
    }

    opt |= opt_from_jcode(orig->jinfo->code);

    if (orig->jinfo->seasonals > 0) {
	opt |= OPT_D;
    }

    jvar = johansen_wrapper(VECM_ESTIMATE, orig->order + 1, 
			    orig->jinfo->rank, list,
			    rset, Z, pdinfo, 
			    opt, prn);

    if (jvar != NULL) {
	jvar->jinfo->ll0 = orig->ll;
	if (!jvar->err) {
	    gretl_VAR_print(jvar, pdinfo, opt, prn);
	} else {
	    *err = jvar->err;
	}
    } else {
	*err = 1;
    }

    free(list);
    
    return jvar;    
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
    int i, vi, vmax = 0;

    if (var->ylist != NULL) {
	for (i=1; i<=var->ylist[0]; i++) {
	    vi = var->ylist[i];
	    if (vi > vmax) {
		vmax = vi;
	    }
	}
    }

    if (var->xlist != NULL) {
	for (i=1; i<=var->xlist[0]; i++) {
	    vi = var->xlist[i];
	    if (vi > vmax) {
		vmax = vi;
	    }
	}
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

int gretl_VECM_n_alpha (const GRETL_VAR *vecm)
{
    int na = 0;

    if (vecm->jinfo != NULL && vecm->jinfo->Alpha != NULL) {
	na = gretl_matrix_rows(vecm->jinfo->Alpha);
    }

    return na;
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
    return vecm->ylist;
}

int beta_restricted_VECM (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL && vecm->jinfo->R != NULL) {
	return 1;
    }

    return 0;
}

int alpha_restricted_VECM (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL && vecm->jinfo->Ra != NULL) {
	return 1;
    }

    return 0;
}

const gretl_matrix *gretl_VECM_R_matrix (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL && vecm->jinfo->R != NULL) {
	return vecm->jinfo->R;
    }

    return NULL;
}

const gretl_matrix *gretl_VECM_q_matrix (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL && vecm->jinfo->q != NULL) {
	return vecm->jinfo->q;
    }

    return NULL;
}   

const gretl_matrix *gretl_VECM_Ra_matrix (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL && vecm->jinfo->Ra != NULL) {
	return vecm->jinfo->Ra;
    }

    return NULL;
}

const gretl_matrix *gretl_VECM_qa_matrix (const GRETL_VAR *vecm)
{
    if (vecm->jinfo != NULL && vecm->jinfo->qa != NULL) {
	return vecm->jinfo->qa;
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
		src = var->jinfo->S00;
		break;
	    case M_JS11:
		src = var->jinfo->S11;
		break;
	    case M_JS01:
		src = var->jinfo->S01;
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

    gretl_VAR_clear(var);

    return var;
}

static int VAR_retrieve_jinfo (xmlNodePtr node, xmlDocPtr doc,
			       GRETL_VAR *var)
{
    xmlNodePtr cur = node->xmlChildrenNode;
    JohansenInfo *jinfo;
    int ID, code, rank;
    int seas;
    int got = 0;
    int err = 0;

    got += gretl_xml_get_prop_as_int(node, "ID", &ID);
    got += gretl_xml_get_prop_as_int(node, "code", &code);
    got += gretl_xml_get_prop_as_int(node, "rank", &rank);
    got += gretl_xml_get_prop_as_int(node, "seasonals", &seas);

    if (got != 4) {
	return E_DATA;
    }

    jinfo = johansen_info_new(var, rank, OPT_NONE);
    if (jinfo == NULL) {
	return E_ALLOC;
    }

    jinfo->ID = ID;
    jinfo->code = code;
    jinfo->rank = rank;
    jinfo->seasonals = seas;

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "u")) {
	    jinfo->R0 = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "v")) {
	    jinfo->R1 = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Suu")) {
	    jinfo->S00 = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Svv")) {
	    jinfo->S11 = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Suv")) {
	    jinfo->S01 = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Beta")) {
	    jinfo->Beta = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Alpha")) {
	    jinfo->Alpha = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Bvar")) {
	    jinfo->Bvar = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Bse")) {
	    jinfo->Bse = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "R")) {
	    jinfo->R = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "q")) {
	    jinfo->q = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Ra")) {
	    jinfo->Ra = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "qa")) {
	    jinfo->qa = gretl_xml_get_matrix(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "ll0")) {
	    gretl_xml_node_get_double(cur, doc, &jinfo->ll0);
	} else if (!xmlStrcmp(cur->name, (XUC) "bdf")) {
	    gretl_xml_node_get_int(cur, doc, &jinfo->bdf);
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

/* for backward compatibility */

static int make_VAR_global_lists (GRETL_VAR *var)
{
    MODEL *pmod = var->models[0];
    int *list = pmod->list;
    int p = var->order;
    int n = var->neqns;
    int nx, np = n * p;
    int i, j, ifc;

    if (list == NULL || list[0] < 3) {
	return E_DATA;
    }

    /* the _endogenous_ vars start in position 2 or 3 (3 if a constant
       is included), and there are (order * neqns) such terms */

    ifc = (list[2] == 0);
    nx = list[0] - 1 - ifc - np;

    var->ylist = gretl_list_new(n);
    if (var->ylist == NULL) {
	return E_ALLOC;
    }

    if (nx > 0) {
	var->xlist = gretl_list_new(nx);
	if (var->xlist == NULL) {
	    free(var->ylist);
	    var->ylist = NULL;
	    return E_ALLOC;
	}
    }

    for (i=0; i<n; i++) {
	pmod = var->models[i];
	var->ylist[i+1] = pmod->list[1];
    }

    j = 2 + ifc + np;
    for (i=0; i<nx; i++) {
	var->xlist[i+1] = list[j++];
    }

    return 0;    
}

GRETL_VAR *gretl_VAR_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    GRETL_VAR *var;
    MODEL *pmod;
    xmlNodePtr cur;
    int start = 0, rowmax = 0;
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

    gretl_xml_get_prop_as_int(node, "detflags", &var->detflags);

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
	if (!xmlStrcmp(cur->name, (XUC) "ylist")) {
	    var->ylist = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "xlist")) {
	    var->xlist = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "Fvals")) {
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

    /* Note: get "global" info from the first model (equation) */
    var->ncoeff = pmod->ncoeff;
    var->ifc = pmod->ifc;
    var->t1 = pmod->t1;
    var->t2 = pmod->t2;
    var->T = var->t2 - var->t1 + 1;

    if (var->ylist == NULL) {
	*err = make_VAR_global_lists(var);
    }

    if (!*err) {
	start = pmod->ifc;
	rowmax = var->neqns * var->order + start;

	/* set up storage for residuals */
	*err = VAR_add_residuals_matrix(var);
    }

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

static void johansen_serialize (JohansenInfo *j, FILE *fp)
{
    fprintf(fp, "<gretl-johansen ID=\"%d\" code=\"%d\" rank=\"%d\" ", 
	    j->ID, j->code, j->rank);
    fprintf(fp, "seasonals=\"%d\">\n", j->seasonals);

    gretl_xml_put_matrix(j->R0, "u", fp);
    gretl_xml_put_matrix(j->R1, "v", fp);
    gretl_xml_put_matrix(j->S00, "Suu", fp);
    gretl_xml_put_matrix(j->S11, "Svv", fp);
    gretl_xml_put_matrix(j->S01, "Suv", fp);
    gretl_xml_put_matrix(j->Beta, "Beta", fp);
    gretl_xml_put_matrix(j->Alpha, "Alpha", fp);
    gretl_xml_put_matrix(j->Bvar, "Bvar", fp);
    gretl_xml_put_matrix(j->Bse, "Bse", fp);
    gretl_xml_put_matrix(j->R, "R", fp);
    gretl_xml_put_matrix(j->q, "q", fp);
    gretl_xml_put_matrix(j->R, "Ra", fp);
    gretl_xml_put_matrix(j->q, "qa", fp);

    if (!na(j->ll0) && j->bdf > 0) {
	gretl_xml_put_double("ll0", j->ll0, fp);
	gretl_xml_put_int("bdf", j->bdf, fp);
    }

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

    fprintf(fp, "ci=\"%d\" neqns=\"%d\" order=\"%d\" detflags=\"%d\" ecm=\"%d\">\n",
	    var->ci, var->neqns, var->order, var->detflags, var->ecm);

    gretl_xml_put_tagged_list("ylist", var->ylist, fp);
    gretl_xml_put_tagged_list("xlist", var->xlist, fp);

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

