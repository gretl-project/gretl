/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/* estimate.c - basic gretl estimation procedures */

#include "libgretl.h"
#include "qr_estimate.h"
#include "gretl_panel.h"
#include "libset.h"
#include "compat.h"
#include "missing_private.h"
#include "estim_private.h"
#include "system.h"
#include "tsls.h"
#include "nls.h"

/* Comment on 'TINY': It's the minimum value for 'test' (see below)
   that libgretl's Cholesky decomposition routine will accept before
   rejecting a data matrix as too highly collinear.  If you set it too
   high, data sets for which Cholesky could produce reasonable
   estimates will be rejected.  If you set it too low (and 100 *
   DBL_EPSILON is definitely too low), gretl will produce more or less
   worthless coefficient estimates when given highly collinear data.
   Before making a permanent change to the value of TINY, check how
   gretl does on the NIST reference data sets for linear regression
   and ensure you're not getting any garbage results.  The current
   enables us to get decent results on the NIST nonlinear regression
   test suite; it might be a bit too low for some contexts.
*/

#define TINY      8.0e-09 /* was 2.1e-09 (last changed 2007/01/20) */
#define SMALL     1.0e-08 /* threshold for printing a warning for collinearity */
#define YBARZERO  0.5e-14 /* threshold for treating mean of dependent
			     variable as effectively zero */
#define ESSZERO   1.0e-22 /* threshold for considering a tiny error-sum-of-
			     squares value to be effectively zero */

#define XPX_DEBUG 0

static void regress (MODEL *pmod, double *xpy, 
		     double ysum, double ypy, 
		     const double **Z, 
		     double rho, gretlopt opt);
static int hatvar (MODEL *pmod, int n, const double **Z);
static void omitzero (MODEL *pmod, const double **Z, const DATAINFO *pdinfo,
		      gretlopt opt);
static int lagdepvar (const int *list, const double **Z, const DATAINFO *pdinfo); 
static int jackknife_vcv (MODEL *pmod, const double **Z);

/* compute statistics for the dependent variable in a model */

static void model_depvar_stats (MODEL *pmod, const double **Z)
{
    double xx, sum = 0.0;
    int yno = pmod->list[1];
    int t, dwt = 0;

    if (pmod->ci == WLS && gretl_model_get_int(pmod, "wt_dummy")) {
	dwt = pmod->nwt;
    }

    pmod->ybar = pmod->sdy = NADBL;

    if (pmod->nobs <= 0) {
	return;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}
	if (!model_missing(pmod, t)) {
	    sum += Z[yno][t];
	} 
    }

    pmod->ybar = sum / pmod->nobs;

    sum = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}	
	if (!model_missing(pmod, t)) {
	    sum += (Z[yno][t] - pmod->ybar); 
	}
    }

    pmod->ybar = pmod->ybar + sum / pmod->nobs;

    if (fabs(pmod->ybar) < YBARZERO) {
	pmod->ybar = 0.0;
    }

    sum = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (dwt && Z[pmod->nwt][t] == 0.0) {
	    continue;
	}
	if (!model_missing(pmod, t)) {
	    xx = Z[yno][t] - pmod->ybar;
	    sum += xx * xx;
	}
    }

    sum = (pmod->nobs > 1)? sum / (pmod->nobs - 1) : 0.0;

    pmod->sdy = (sum >= 0)? sqrt(sum) : NADBL;
}

/* determine the degrees of freedom for a model */

static int get_model_df (MODEL *pmod)
{
    int err = 0;

    pmod->ncoeff = pmod->list[0] - 1;
    pmod->dfd = pmod->nobs - pmod->ncoeff;

    if (pmod->dfd < 0) {
	pmod->errcode = E_DF;
        sprintf(gretl_errmsg, _("No. of obs (%d) is less than no. "
		"of parameters (%d)"), pmod->nobs, pmod->ncoeff);
	err = 1;
    } else {
	pmod->dfn = pmod->ncoeff - pmod->ifc;
    }

    return err;
}

#define LDDEBUG 0

static int
transcribe_ld_vcv (MODEL *targ, MODEL *src)
{
    int nv = targ->ncoeff;
    int nxpx = (nv * nv + nv) / 2;
    int i, j;

    if (makevcv(src, src->sigma)) {
	return 1;
    }

    if (targ->vcv == NULL) {
	targ->vcv = malloc(nxpx * sizeof *targ->vcv);
	if (targ->vcv == NULL) {
	    return 1;
	}
    }

    for (i=0; i<nv; i++) {
	for (j=i; j<nv; j++) {
	    targ->vcv[ijton(i, j, nv)] = 
		src->vcv[ijton(i, j, src->ncoeff)];
	}
    }

    return 0;
}

/* Calculate consistent standard errors (and VCV matrix) when doing
   AR(1) estimation of a model with lagged dependent variable.  See
   Ramanathan, Introductory Econometrics, 5e, p. 450.
*/

static int 
ldepvar_std_errors (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    MODEL emod;
    const double *x;
    int orig_t1 = pdinfo->t1;
    int orig_t2 = pdinfo->t2;
    double rho = gretl_model_get_double(pmod, "rho_in");
    int origv = pdinfo->v;
    int vnew = pmod->list[0] + 1 - pmod->ifc;
    int *list;
    int vi, vm;
    int i, t;
    int err = 0;

#if LDDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);
    printlist(pmod->list, "pmod->list");
    printf("vnew = %d\n", vnew);
    printf("rho = %g\n", rho);
#endif

    list = gretl_list_new(vnew + pmod->ifc);
    if (list == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    err = dataset_add_series(vnew, pZ, pdinfo);
    if (err) {
	free(list);
	pmod->errcode = E_ALLOC;
	return 1;
    }

    vi = origv;

    /* dependent var is residual from original model */
    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[vi][t] = pmod->uhat[t];
    }    
    strcpy(pdinfo->varname[vi], "eps");
    list[1] = vi++;

    /* indep vars are rho-differenced vars from original model */
    for (i=2; i<=pmod->list[0]; i++) {
	vm = pmod->list[i];
	if (vm == 0) {
	    list[i] = 0;
	    continue;
	}
	sprintf(pdinfo->varname[vi], "%.6s_r", pdinfo->varname[vm]);
	x = (*pZ)[vm];
	for (t=0; t<pdinfo->n; t++) {
	    if (t == 0 || na(x[t]) || na(x[t-1])) {
		(*pZ)[vi][t] = NADBL;
	    } else {
		(*pZ)[vi][t] = x[t] - rho * x[t-1];
	    }
	}
	list[i] = vi++;
    }

    /* last indep var is lagged u-hat */
    for (t=0; t<pdinfo->n; t++) {
	if (t == 0) {
	    (*pZ)[vi][t] = NADBL;
	} else { 
	    (*pZ)[vi][t] = (*pZ)[pmod->list[1]][t-1];
	}
	if (na((*pZ)[vi][t])) {
	    continue;
	}
	for (i=0; i<pmod->ncoeff; i++) {
	    x = (*pZ)[pmod->list[i+2]];
	    if (na(x[t-1])) {
		(*pZ)[vi][t] = NADBL;
		break;
	    } else {
		(*pZ)[vi][t] -= pmod->coeff[i] * x[t-1];
	    }
	}
    }

    list[list[0]] = vi;
    strcpy(pdinfo->varname[vi], "uhat_1");

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    emod = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (emod.errcode) {
	err = emod.errcode;
    } else {
#if LDDEBUG
	printmodel(&emod, pdinfo, OPT_NONE, prn);
	gretl_print_destroy(prn);
#endif
	for (i=0; i<pmod->ncoeff; i++) {
	    pmod->sderr[i] = emod.sderr[i];
	}

	err = transcribe_ld_vcv(pmod, &emod);
    }
    
    clear_model(&emod);
    
    free(list);
    dataset_drop_last_variables(vnew, pZ, pdinfo);

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    if (err) {
	pmod->errcode = err;
    }

    return err;
}

/* special computation of statistics for autoregressive models */

static int compute_ar_stats (MODEL *pmod, const double **Z, double rho)
{
    int i, t, yno = pmod->list[1];
    int pwe = (pmod->opt & OPT_P);
    double x, pw1 = 0.0;

    if (gretl_model_add_arinfo(pmod, 1)) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    if (pwe) {
	pw1 = sqrt(1.0 - rho * rho);
    }

    pmod->arinfo->arlist[0] = pmod->arinfo->arlist[1] = 1;
    pmod->arinfo->rho[0] = rho;
    gretl_model_set_double(pmod, "rho_in", rho);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (t == pmod->t1 && pwe) {
	    x = pw1 * Z[yno][t];
	    for (i=pmod->ifc; i<pmod->ncoeff; i++) {
		x -= pmod->coeff[i] * pw1 * Z[pmod->list[i+2]][t];
	    }
	    if (pmod->ifc) {
		x -= pw1 * pmod->coeff[0];
	    }
	} else {
	    x = Z[yno][t] - rho * Z[yno][t-1];
	    for (i=0; i<pmod->ncoeff; i++) {
		x -= pmod->coeff[i] * 
		    (Z[pmod->list[i+2]][t] - 
		     rho * Z[pmod->list[i+2]][t-1]);
	    }
	}
	pmod->uhat[t] = x;
	pmod->yhat[t] = Z[yno][t] - x;
    }

    pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, Z[yno], pmod->yhat);

    pmod->adjrsq = 
	1.0 - ((1.0 - pmod->rsq) * (pmod->t2 - pmod->t1) / 
	       (double) pmod->dfd);

    return 0;
}

/* calculation of WLS stats in agreement with GNU R */

static void get_wls_stats (MODEL *pmod, const double **Z)
{
    int dumwt = gretl_model_get_int(pmod, "wt_dummy");
    int t, wobs = pmod->nobs, yno = pmod->list[1];
    double x, dy, wmean = 0.0, wsum = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) { 
	if (model_missing(pmod, t)) {
	    continue;
	}
	if (Z[pmod->nwt][t] == 0.0 && !dumwt) {
	    wobs--;
	    pmod->dfd -= 1;
	} else {
	    wmean += Z[pmod->nwt][t] * Z[yno][t];
	    wsum += Z[pmod->nwt][t];
	}
    }

    wmean /= wsum;
    x = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t) || Z[pmod->nwt][t] == 0.0) {
	    continue;
	}	
	dy = Z[yno][t] - wmean;
	x += Z[pmod->nwt][t] * dy * dy;
    }

    pmod->fstt = ((x - pmod->ess) * pmod->dfd) / (pmod->dfn * pmod->ess);
    pmod->rsq = (1 - (pmod->ess / x));
    pmod->adjrsq = 1 - ((1 - pmod->rsq) * (pmod->nobs - 1)/pmod->dfd);
}

static void fix_wls_values (MODEL *pmod, double **Z)
{
    int t;

    if (gretl_model_get_int(pmod, "wt_dummy")) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    }
	}
    } else {
	double ess_orig = 0.0;
	double sw, sigma_orig;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    if (Z[pmod->nwt][t] == 0.0) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
		pmod->nobs -= 1;
	    } else {
		sw = sqrt(Z[pmod->nwt][t]);
		pmod->yhat[t] /= sw;
		pmod->uhat[t] /= sw;
		ess_orig += pmod->uhat[t] * pmod->uhat[t];
	    }
	}

	sigma_orig = sqrt(ess_orig / pmod->dfd);
	gretl_model_set_double(pmod, "ess_orig", ess_orig);
	gretl_model_set_double(pmod, "sigma_orig", sigma_orig);
    }
}

/* drop the weight var from the list of regressors (WLS) */

static void dropwt (int *list)
{
    int i;

    list[0] -= 1;
    for (i=1; i<=list[0]; i++) {
	list[i] = list[i+1];
    }
}

static void model_stats_init (MODEL *pmod)
{
    pmod->ess = NADBL;
    pmod->sigma = NADBL;
    pmod->fstt = pmod->lnL = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->chisq = NADBL;
    pmod->dw = NADBL;
}

#define SMPL_DEBUG 0

static int 
lsq_check_for_missing_obs (MODEL *pmod, gretlopt opts,
			   DATAINFO *pdinfo, const double **Z, 
			   int *misst)
{
    int ref_mask = reference_missmask_present();
    int missv = 0;
    int reject_missing = 0;

#if SMPL_DEBUG
    fprintf(stderr, "lsq_check_for_missing_obs: ref_mask = %d\n",
	    ref_mask);
#endif

    if (ref_mask) {
	int err = apply_reference_missmask(pmod);

	/* If there was a reference mask present, it was put there
	   as part of a hypothesis test on some original model, and
	   it should be respected in estimation of this model */

	if (err) {
	    pmod->errcode = E_ALLOC;
	    return 1;
	} else {
	    return 0;
	}
    } 

    /* can't do HAC VCV with missing obs in middle */
    if ((opts & OPT_R) && dataset_is_time_series(pdinfo) &&
	!libset_get_bool(FORCE_HC)) {
	reject_missing = 1;
    } 

    if (opts & OPT_M) {
	reject_missing = 1;
    }

    if (reject_missing) {
	/* reject missing obs within adjusted sample */
	missv = adjust_t1t2(pmod, pmod->list, &pmod->t1, &pmod->t2,
			    pdinfo->n, Z, misst);
    } else {
	/* we'll try to work around missing obs */
	missv = adjust_t1t2(pmod, pmod->list, &pmod->t1, &pmod->t2,
			    pdinfo->n, Z, NULL);
    }

#if SMPL_DEBUG
    if (1) {
	char t1s[OBSLEN], t2s[OBSLEN];
	int misscount = model_missval_count(pmod);

	ntodate(t1s, pmod->t1, pdinfo);
	ntodate(t2s, pmod->t2, pdinfo);
	fprintf(stderr, "*** after adjustment, t1=%d (%s), t2=%d (%s)\n", 
		pmod->t1, t1s, pmod->t2, t2s);
	fprintf(stderr, "Valid observations in range = %d\n", 
		pmod->t2 - pmod->t1 + 1 - misscount);
    }
#endif

    return missv;
}

static int
lagged_depvar_check (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    int ldv = lagdepvar(pmod->list, Z, pdinfo);

    if (ldv) {
	gretl_model_set_int(pmod, "ldepvar", ldv);
    } else if (gretl_model_get_int(pmod, "ldepvar")) {
	gretl_model_set_int(pmod, "ldepvar", 0);
    }

    return ldv;
}

static void 
log_depvar_ll (MODEL *pmod, const double **Z, const DATAINFO *pdinfo)
{
    char parent[VNAMELEN];

    if (is_log_variable(pmod->list[1], pdinfo, parent)) {
	double jll = pmod->lnL;
	int t;

	for (t=0; t<pdinfo->n; t++) {
	    if (!na(pmod->uhat[t])) {
		jll -= Z[pmod->list[1]][t];
	    }
	}
	gretl_model_set_double(pmod, "jll", jll);
	gretl_model_set_string_as_data(pmod, 
				       "log-parent", 
				       gretl_strdup(parent));
    }
}

static int check_weight_var (MODEL *pmod, const double *w, int *effobs)
{
    const char *wtzero = 
	N_("Weight variable is all zeros, aborting regression");
    const char *wtneg = 
	N_("Weight variable contains negative values");
    int t;

    if (gretl_iszero(pmod->t1, pmod->t2, w)) {
	strcpy(gretl_errmsg, _(wtzero));
	pmod->errcode = E_DATA;
	return 1;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (w[t] < 0.0) {
	    strcpy(gretl_errmsg, _(wtneg));
	    pmod->errcode = E_DATA;
	    return 1;
	}
    }

    *effobs = gretl_isdummy(pmod->t1, pmod->t2, w);

    if (*effobs) {
	/* the weight var is a dummy, with effobs 1s */
	gretl_model_set_int(pmod, "wt_dummy", 1);
    }

    return 0;
}

void 
maybe_shift_ldepvar (MODEL *pmod, const double **Z, DATAINFO *pdinfo)
{
    if (gretl_model_get_int(pmod, "ldepvar")) {
	lagged_depvar_check(pmod, Z, pdinfo);
    }
}

/**
 * XTX_XTy:
 * @list: list of variables in model.
 * @t1: starting observation.
 * @t2: ending observation.
 * @Z: data array.
 * @nwt: ID number of variable used as weight, or 0.
 * @rho: quasi-differencing coefficent, or 0.0;
 * @pwe: if non-zero, use Prais-Winsten for first observation.
 * @xpx: on output, X'X matrix as lower triangle, stacked by columns.
 * @xpy: on output, X'y vector (but see below).
 * @ysum: location to recieve \sum y_i, or %NULL.
 * @ypy: location to recieve (scalar) y'y, or %NULL.
 * @mask: missing observations mask, or %NULL.
 *
 * Calculates X'X and X'y, with various possible transformations
 * of the original data, depending on @nwt, @rho and @pwe.
 * If X'y is not required, @xpy can be given as %NULL.
 *
 * Note: the y-related pointer arguments @xpy, @ysum, and @ypy 
 * form a "package": either all should be given, or all should 
 * be %NULL.
 *
 * Returns: 0 on success, non-zero on error.
 */

static int XTX_XTy (const int *list, int t1, int t2, 
		    const double **Z, int nwt, double rho, int pwe,
		    double *xpx, double *xpy, 
		    double *ysum, double *ypy,
		    const char *mask)
{
    int yno = list[1];
    int lmin = (xpy != NULL)? 2 : 1;
    int lmax = list[0];
    int qdiff = (rho != 0.0);
    double x, pw1;
    int vi, vj, m;
    int i, j, t;

    /* Prais-Winsten term */
    if (qdiff && pwe) {
	pw1 = sqrt(1.0 - rho * rho);
    } else {
	pwe = 0;
	pw1 = 0.0;
    }

    if (xpy != NULL) {
	*ysum = *ypy = 0.0;

	for (t=t1; t<=t2; t++) {
	    if (masked(mask, t)) {
		continue;
	    }
	    x = Z[yno][t]; 
	    if (qdiff) {
		if (pwe && t == t1) {
		    x = pw1 * Z[yno][t];
		} else {
		    x -= rho * Z[yno][t-1];
		}
	    } else if (nwt) {
		x *= sqrt(Z[nwt][t]);
	    }
	    *ysum += x;
	    *ypy += x * x;
	}

	if (*ypy <= 0.0) {
	    /* error condition */
	    return yno; 
	} 
    }   

    m = 0;

    if (qdiff) {
	/* quasi-difference the data */
	for (i=lmin; i<=lmax; i++) {
	    vi = list[i];
	    for (j=i; j<=lmax; j++) {
		vj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (pwe && t == t1) {
			x += pw1 * Z[vi][t1] * pw1 * Z[vj][t];
		    } else {
			x += (Z[vi][t] - rho * Z[vi][t-1]) * 
			    (Z[vj][t] - rho * Z[vj][t-1]);
		    }
		}
		if (vi == vj && x < DBL_EPSILON)  {
		    return E_SINGULAR;
		}
		xpx[m++] = x;
	    }
	    if (xpy != NULL) {
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (pwe && t == t1) {
			x += pw1 * Z[yno][t] * pw1 * Z[vi][t];
		    } else {
			x += (Z[yno][t] - rho * Z[yno][t-1]) *
			    (Z[vi][t] - rho * Z[vi][t-1]);
		    }
		}
		xpy[i-2] = x;
	    }
	}
    } else if (nwt) {
	/* weight the data */
	for (i=lmin; i<=lmax; i++) {
	    vi = list[i];
	    for (j=i; j<=lmax; j++) {
		vj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!masked(mask, t)) {
			x += Z[nwt][t] * Z[vi][t] * Z[vj][t];
		    }
		}
		if (vi == vj && x < DBL_EPSILON) {
		    return E_SINGULAR;
		}   
		xpx[m++] = x;
	    }
	    if (xpy != NULL) {
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!masked(mask, t)) {
			x += Z[nwt][t] * Z[yno][t] * Z[vi][t];
		    }
		}
		xpy[i-2] = x;
	    }
	}
    } else {
	/* no quasi-differencing or weighting wanted */
	for (i=lmin; i<=lmax; i++) {
	    vi = list[i];
	    for (j=i; j<=lmax; j++) {
		vj = list[j];
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!masked(mask, t)) {
			x += Z[vi][t] * Z[vj][t];
		    }
		}
		if (vi == vj && x < DBL_EPSILON) {
		    return E_SINGULAR;
		}
		xpx[m++] = x;
	    }
	    if (xpy != NULL) {
		x = 0.0;
		for (t=t1; t<=t2; t++) {
		    if (!masked(mask, t)) {
			x += Z[yno][t] * Z[vi][t];
		    }
		}
		xpy[i-2] = x;
	    }
	}
    }

    return 0; 
}

/**
 * gretl_XTX:
 * @pmod: reference model.
 * @Z: data array.
 * @err: location to receive error code.
 *
 * (Re-)calculates X'X, with various possible transformations
 * of the original data depending on whether estimation of 
 * @pmod involved weighting or quasi-differencing.
 *
 * Returns: The vech of X'X or %NULL on error.
 */

double *gretl_XTX (const MODEL *pmod, const double **Z, int *err)
{
    int *xlist;
    double *xpx;
    double rho;
    int pwe = 0;
    int k, m;
    
    *err = 0;

    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	*err = E_DATA;
	return NULL;
    }

    k = xlist[0];
    m = k * (k + 1) / 2;

    xpx = malloc(m * sizeof *xpx);
    if (xpx == NULL) {
	*err = E_ALLOC;
	free(xlist);
	return NULL;
    }

    if (pmod->ci == AR1 && (pmod->opt & OPT_P)) {
	pwe = 1;
    }

    rho = gretl_model_get_double(pmod, "rho_in");
    if (na(rho)) {
	rho = 0.0;
    }

    *err = XTX_XTy(xlist, pmod->t1, pmod->t2, Z, pmod->nwt, 
		   rho, pwe, xpx, NULL, NULL, NULL, pmod->missmask);

    free(xlist);

    return xpx;
}

static int gretl_choleski_regress (MODEL *pmod, const double **Z, 
				   double rho, int pwe, gretlopt opt)
{
    double ysum = 0.0, ypy = 0.0;
    double *xpy;
    int k = pmod->ncoeff;
    int nxpx = k * (k + 1) / 2;
    int i;

    if (nxpx == 0) {
	fprintf(stderr, "problem: nxpx = 0\n");
	pmod->errcode = E_DATA;
	return pmod->errcode;
    }

    xpy = malloc(k * sizeof *xpy);
    if (xpy == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    pmod->xpx = malloc(nxpx * sizeof *pmod->xpx);
    pmod->coeff = malloc(k * sizeof *pmod->coeff);
    pmod->sderr = malloc(k * sizeof *pmod->sderr);

    if (pmod->yhat == NULL) {
	pmod->yhat = malloc(pmod->full_n * sizeof *pmod->yhat);
    } 

    if (pmod->uhat == NULL) {
	pmod->uhat = malloc(pmod->full_n * sizeof *pmod->uhat);
    }

    if (pmod->xpx == NULL || pmod->coeff == NULL ||
	pmod->sderr == NULL || pmod->yhat == NULL || 
	pmod->uhat == NULL) {
	free(xpy);
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    for (i=0; i<k; i++) {
	xpy[i] = 0.0;
    }
    for (i=0; i<nxpx; i++) {
	pmod->xpx[i] = 0.0;
    }

    /* calculate regression results, Cholesky style */
    pmod->errcode = XTX_XTy(pmod->list, pmod->t1, pmod->t2, Z, pmod->nwt, 
			    rho, pwe, pmod->xpx, xpy, 
			    &ysum, &ypy, pmod->missmask);

#if XPX_DEBUG
    for (i=0; i<k; i++) {
	fprintf(stderr, "xpy[%d] = %g\n", i, xpy[i]);
    }
    for (i=0; i<nxpx; i++) {
	fprintf(stderr, "xpx[%d] = %g\n", i, pmod->xpx[i]);
    }
    fputc('\n', stderr);
#endif

    if (!pmod->errcode) {
	regress(pmod, xpy, ysum, ypy, Z, rho, opt);
    }

    free(xpy);

    return pmod->errcode;
}

static int gretl_null_regress (MODEL *pmod, const double **Z)
{
    double yt;
    int t;

    if (pmod->yhat == NULL) {
	pmod->yhat = malloc(pmod->full_n * sizeof *pmod->yhat);
    } 

    if (pmod->uhat == NULL) {
	pmod->uhat = malloc(pmod->full_n * sizeof *pmod->uhat);
    }

    if (pmod->yhat == NULL || pmod->uhat == NULL) {
	pmod->errcode = E_ALLOC;
	return pmod->errcode;
    }

    pmod->nobs = 0;
    pmod->ifc = 0;
    pmod->ess = 0.0;
    pmod->rsq = pmod->adjrsq = 0.0;

    for (t=0; t<pmod->full_n; t++) {
	yt = Z[pmod->list[1]][t];
	if (t < pmod->t1 || t > pmod->t2 || na(yt)) {
	    pmod->uhat[t] = pmod->yhat[t] = NADBL;
	} else {
	    pmod->uhat[t] = yt;
	    pmod->yhat[t] = 0.0;
	    pmod->ess += yt * yt;
	    pmod->nobs += 1;
	}
    }

    if (pmod->ess == 0) {
	pmod->sigma = 0.0;
    } else if (pmod->nobs > 1) {
	pmod->sigma = sqrt(pmod->ess / (pmod->nobs - 1));
    } else {
	pmod->errcode = E_DATA;
    }

    return pmod->errcode;
}

/* check whether the variable with ID @yno is all zeros;
   return 1 if so, otherwise return 0 */

static int depvar_zero (const MODEL *pmod, int yno, const double **Z)
{
    double yt;
    int t, ret = 1;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}
	yt = Z[yno][t];
	if (pmod->nwt) {
	    yt *= Z[pmod->nwt][t];
	}
	if (yt != 0.0) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/* limited freeing of elements before passing a model
   on for QR estimation in the case of (near-)singularity 
*/

static void model_free_storage (MODEL *pmod)
{
    free(pmod->xpx);
    free(pmod->coeff);
    free(pmod->sderr);

    pmod->xpx = NULL;
    pmod->coeff = NULL;
    pmod->sderr = NULL;
}

/* as lsq() below, except that we allow for a non-zero value
   of the first-order quasi-differencing coefficient, rho,
   and there's no PRN.
*/

MODEL ar1_lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	       GretlCmdIndex ci, gretlopt opt, double rho)
{
    MODEL mdl;
    int effobs = 0;
    int missv = 0, misst = 0;
    int jackknife = 0;
    int pwe = (opt & OPT_P);
    int nullmod = 0, ldv = 0;
    int yno, i;

    gretl_error_clear();

    if (list == NULL || pZ == NULL || pdinfo == NULL) {
	fprintf(stderr, "E_DATA: lsq: list = %p, pZ = %p, pdinfo = %p\n",
		(void *) list, (void *) pZ, (void *) pdinfo);
	mdl.errcode = E_DATA;
        return mdl;
    }

    if (ci == HSK) {
	return hsk_func(list, pZ, pdinfo);
    } 

    gretl_model_init(&mdl);
    gretl_model_smpl_init(&mdl, pdinfo);
    model_stats_init(&mdl);

    if (ci == AR1) {
	if (opt & OPT_P) {
	    mdl.opt |= OPT_P;
	} else if (opt & OPT_H) {
	    mdl.opt |= OPT_H;
	}
    } 

    if (list[0] == 1 && ci == OLS && (opt & OPT_U)) {
	/* null model OK */
	nullmod = 1;
    } else if (list[0] == 1 || pdinfo->v == 1) {
	fprintf(stderr, "E_DATA: lsq: list[0] = %d, pdinfo->v = %d\n",
		list[0], pdinfo->v);
	mdl.errcode = E_DATA;
        return mdl;
    }

    /* preserve a copy of the list supplied, for future reference */
    mdl.list = gretl_list_copy(list);
    if (mdl.list == NULL) {
        mdl.errcode = E_ALLOC;
        return mdl;
    }

    mdl.t1 = pdinfo->t1;
    mdl.t2 = pdinfo->t2;
    mdl.full_n = pdinfo->n;
    mdl.ci = ci;

    /* Doing weighted least squares? */
    if (ci == WLS) { 
	check_weight_var(&mdl, (*pZ)[mdl.list[1]], &effobs);
	if (mdl.errcode) {
	    return mdl;
	}
	mdl.nwt = mdl.list[1];
    } else {
	mdl.nwt = 0;
    }

    /* sanity check */
    if (mdl.t1 < 0 || mdl.t2 > pdinfo->n - 1) {
        mdl.errcode = E_DATA;
        goto lsq_abort;
    }

    /* adjust sample range and check for missing obs: this
       may set the model errcode */
    missv = lsq_check_for_missing_obs(&mdl, opt, pdinfo,
				      (const double **) *pZ,
				      &misst);
    if (mdl.errcode) {
        goto lsq_abort;
    }

    /* react to presence of unhandled missing obs */
    if (missv) {
	if (dated_daily_data(pdinfo)) {
	    if (repack_missing_daily_obs(&mdl, *pZ, pdinfo)) {
		return mdl;
	    }
	} else {
	    sprintf(gretl_errmsg, _("Missing value encountered for "
		    "variable %d, obs %d"), missv, misst);
	    mdl.errcode = E_DATA;
	    return mdl;
	} 
    }

    if (ci == WLS) {
	dropwt(mdl.list);
    }

    yno = mdl.list[1];
    
    /* check for unknown vars in list */
    for (i=1; i<=mdl.list[0]; i++) {
        if (mdl.list[i] > pdinfo->v - 1) {
            mdl.errcode = E_UNKVAR;
            goto lsq_abort;
        }
    } 

    /* check for zero dependent var */
    if (depvar_zero(&mdl, yno, (const double **) *pZ)) {  
        mdl.errcode = E_ZERO;
        goto lsq_abort; 
    } 

    /* drop any vars that are all zero and repack the list */
    omitzero(&mdl, (const double **) *pZ, pdinfo, opt);

    /* if regressor list contains a constant, record this fact and 
       place it first among the regressors */
    mdl.ifc = reglist_check_for_const(mdl.list, (const double **) *pZ, 
				      pdinfo);

    /* Check for presence of lagged dependent variable? 
       (Don't bother if this is an auxiliary regression.) */
    if (!(opt & OPT_A)) {
	ldv = lagged_depvar_check(&mdl, (const double **) *pZ, pdinfo);
    }

    /* AR1: advance the starting observation by one? */
    if (rho != 0.0 && !pwe) {
	mdl.t1 += 1;
    }

    mdl.ncoeff = mdl.list[0] - 1; 
    if (effobs) {
	mdl.nobs = effobs; /* FIXME? */
    } else {
	mdl.nobs = mdl.t2 - mdl.t1 + 1;
	if (has_missing_obs(&mdl)) {
	    mdl.nobs -= model_missval_count(&mdl);
	}
    }

    /* check degrees of freedom */
    if (get_model_df(&mdl)) {
        goto lsq_abort; 
    }

    /* if df correction is not wanted, record this fact */
    if (opt & OPT_N) {
	mdl.opt |= OPT_N;
    }

    if (dataset_is_time_series(pdinfo)) {
	opt |= OPT_T;
    } 

    /* remove Durbin-Watson p-value flag if not appropriate */
    if (ldv || !(opt & OPT_T)) {
	opt &= ~OPT_I;
    }

    if ((opt & OPT_J) || ((opt & OPT_R) && libset_get_int(HC_VERSION) == 4)) {
	jackknife = 1;
    }

    if (nullmod) {
	gretl_null_regress(&mdl, (const double **) *pZ);
    } else if (!jackknife && (opt & (OPT_R | OPT_I | OPT_Q))) { 
	mdl.rho = rho;
	gretl_qr_regress(&mdl, (const double **) *pZ, pdinfo, opt);
    } else {
	gretl_choleski_regress(&mdl, (const double **) *pZ, rho, pwe, opt);
	if (mdl.errcode == E_SINGULAR && !jackknife) {
	    /* (near-) perfect collinearity is better handled by QR */
	    model_free_storage(&mdl);
	    mdl.rho = rho;
	    gretl_qr_regress(&mdl, (const double **) *pZ, pdinfo, opt);
	}
    }

    if (mdl.errcode) {
	goto lsq_abort;
    }

    /* get the mean and sd of dep. var. and make available */
    model_depvar_stats(&mdl, (const double **) *pZ);

    /* Doing an autoregressive procedure? */
    if (ci == AR1) {
	if (compute_ar_stats(&mdl, (const double **) *pZ, rho)) { 
	    goto lsq_abort;
	}
	if (ldv) {
	    if (ldepvar_std_errors(&mdl, pZ, pdinfo)) {
		goto lsq_abort;
	    }
	}
	if ((opt & OPT_H) && (opt & OPT_B)) {
	    gretl_model_set_int(&mdl, "no-corc", 1);
	}
    }

    /* weighted least squares: fix yhat and uhat; add calculation of
       ESS and sigma based on unweighted data
    */
    if (ci == WLS) {
	get_wls_stats(&mdl, (const double **) *pZ);
	fix_wls_values(&mdl, *pZ);
    }

    if (mdl.missmask == NULL && (opt & OPT_T) && !(opt & OPT_I)) {
	mdl.rho = rhohat(1, mdl.t1, mdl.t2, mdl.uhat);
	mdl.dw = dwstat(1, &mdl, (const double **) *pZ);
    } else if (!(opt & OPT_I)) {
	mdl.rho = mdl.dw = NADBL;
    }

    /* weird special case: degenerate model */
    if (mdl.ncoeff == 1 && mdl.ifc) {
	mdl.rsq = mdl.adjrsq = 0.0;
	mdl.fstt = NADBL;
    }

    /* Generate model selection statistics */
    ls_criteria(&mdl);
    if (!(opt & OPT_A) && !na(mdl.lnL)) {
	log_depvar_ll(&mdl, (const double **) *pZ, pdinfo);
    }

    /* hccm command or HC3a */
    if (jackknife) {
	mdl.errcode = jackknife_vcv(&mdl, (const double **) *pZ);
    }

 lsq_abort:

    /* If we reshuffled any missing observations, put them
       back in their right places now */
    if (gretl_model_get_int(&mdl, "daily_repack")) {
	undo_daily_repack(&mdl, *pZ, pdinfo);
    }

    if (!(opt & OPT_A)) {
	/* if it's not an auxiliary regression, set an ID number
	   on the model */
	set_model_id(&mdl);
    }

    return mdl;
}

/**
 * lsq:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @ci: one of the command indices in #LSQ_MODEL.
 * @opt: option flags: zero or more of the following --
 *   %OPT_R compute robust standard errors;
 *   %OPT_A treat as auxiliary regression (don't bother checking
 *     for presence of lagged dependent var, don't augment model count);
 *   %OPT_P use Prais-Winsten for first obs;
 *   %OPT_N don't use degrees of freedom correction for standard
 *      error of regression;
 *   %OPT_M reject missing observations within sample range;
 *   %OPT_Z (internal use) suppress the automatic elimination of 
 *      perfectly collinear variables.
 *   %OPT_X: compute "variance matrix" as just (X'X)^{-1}
 *   %OPT_B: don't compute R^2.
 *   %OPT_I: compute Durbin-Watson p-value.
 *   %OPT_Q: use QR decomposition (not necessarily robust VCV).
 *
 * Computes least squares estimates of the model specified by @list,
 * using an estimator determined by the value of @ci.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	   GretlCmdIndex ci, gretlopt opt)
{
    return ar1_lsq(list, pZ, pdinfo, ci, opt, 0.0);
}


static int make_ess (MODEL *pmod, const double **Z)
{
    int i, t, yno = pmod->list[1], l0 = pmod->list[0];
    int nwt = pmod->nwt;
    double yhat, resid;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (nwt && Z[nwt][t] == 0.0) {
	    continue;
	}
	if (model_missing(pmod, t)) {
	    continue;
	}
	yhat = 0.0;
	for (i=2; i<=l0; i++) {
	    yhat += pmod->coeff[i-2] * Z[pmod->list[i]][t];
	}
	resid = Z[yno][t] - yhat;
	if (nwt) {
	    resid *= sqrt(Z[nwt][t]);
	}
	pmod->ess += resid * resid;
    }

    return 0;
}

#define SMALLDIFF 9.0e-16

/* The heuristic used here is that the model effectively
   contains a constant or intercept if the means of y and
   yhat are the same, where "the same" means that the
   relative difference is less than SMALLDIFF.
*/

int check_for_effective_const (MODEL *pmod, const double *y)
{
    double x1 = 0.0, x2 = 0.0;
    double reldiff;
    int t, ret = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    x1 += pmod->yhat[t];
	    x2 += y[t];
	}
    }

    reldiff = fabs((x1 - x2) / x2);
#if 0
    fprintf(stderr, "check_for_effective_const: reldiff = %g\n", reldiff);
#endif

    if (reldiff < SMALLDIFF) {
	gretl_model_set_int(pmod, "effconst", 1);
	pmod->dfn -= 1;
	ret = 1;
    } else if (gretl_model_get_int(pmod, "effconst")) {
	gretl_model_set_int(pmod, "effconst", 0);
	pmod->dfn += 1;
    }

    return ret;
}

static void uncentered_r_squared (MODEL *pmod, const double *y)
{
    double y0 = y[pmod->t1];

    /* special computation for the case of TSS = 0, i.e.
       the dependent variable is a constant */

    if (y0 > 0) {
	double tss = pmod->nobs * y0 * y0;

	pmod->rsq = 1 - (pmod->ess / tss);
	gretl_model_set_int(pmod, "uncentered", 1);
    }
}

static void compute_r_squared (MODEL *pmod, const double *y, int *ifc)
{
    pmod->rsq = 1.0 - (pmod->ess / pmod->tss);

    if (*ifc == 0) {
	*ifc = check_for_effective_const(pmod, y);
    }

    if (pmod->dfd > 0) {
	double den = 0.0;

	if (*ifc) {
	    den = pmod->tss * pmod->dfd;
	    pmod->adjrsq = 1 - (pmod->ess * (pmod->nobs - 1) / den);
	} else {
	    /* model does not contain a constant */
	    int t;

	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!na(pmod->yhat[t])) {
		    den += y[t] * y[t];
		}
	    }

	    /* make the centered R^2 available for output */
	    gretl_model_set_double(pmod, "centered-R2", pmod->rsq);

	    /* but make the "official" figure the uncentered R^2,
	       as per NIST, R, Stata, SPSS, ... */
	    pmod->rsq = 1 - pmod->ess / den;
	    pmod->adjrsq = 
		1.0 - ((1.0 - pmod->rsq) * (pmod->nobs - 1.0) / pmod->dfd);
	} 
    }

    if (pmod->rsq < 0.0) {
	pmod->rsq = 0.0;
    }
}

/*
  diaginv: solves for the diagonal elements of X'X inverse.

  xpx = Cholesky-decomposed X'X
  tmp = work array, length >= nv
  diag = diagonal elements of {X'X}^{-1} (output)
  nv = number of regression coefficients = length of diag
*/

static void diaginv (const double *xpx, double *tmp, double *diag, 
		     int nv)
{
    int kk, l, m, k, i, j;
    double d, e;

    kk = 0;

    for (l=0; l<nv-1; l++) {
        d = xpx[kk];
        tmp[l] = d;
        e = d * d;
        m = 0;
        if (l > 0) {
	    for (j=1; j<=l; j++) {
		m += nv - j;
	    }
	}
        for (i=l+1; i<nv; i++) {
            d = 0.0;
            k = i + m;
            for (j=l; j<i; j++) {
                d += tmp[j] * xpx[k];
                k += nv - (j+1);
            }
            d = -d * xpx[k];
            tmp[i] = d;
            e += d * d;
        }
	kk += nv - l;
        diag[l] = e;
    }

    diag[nv-1] = xpx[kk] * xpx[kk];
}

/*
  cholbeta: in-place Cholesky decomposition of X'X (pmod->xpx) (lower
  triangular matrix stacked in columns) plus solution for the
  least-squares coefficient estimates.

  pmod->xpx = X'X on input and Cholesky decomposition on output
  xpy = the X'y vector on input and Cholesky-transformed t
        vector on output 
  rss = location to receive regression sum of squares

  The number of floating-point operations is basically 3.5 * nv^2
  plus (nv^3) / 3.
*/

static int cholbeta (MODEL *pmod, double *xpy, double *rss)
{
    int i, j, k, kk, l, jm1;
    double e, d, d1, d2, test, xx;
    double *xpx = pmod->xpx;
    double *coeff = pmod->coeff;
    int nc = pmod->ncoeff;

    if (xpx[0] <= 0.0) {
	fprintf(stderr, "%s %d: xpx <= 0.0\n", __FILE__, __LINE__);
	return E_NAN;
    }

    e = 1.0 / sqrt(xpx[0]);
    xpx[0] = e;
    xpy[0] *= e;
    for (i=1; i<nc; i++) {
	xpx[i] *= e;
    }

    kk = nc;

    for (j=1; j<nc; j++) {
	/* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j;

        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d1 += xx * xpy[l-1];
            d += xx * xx;
            k += nc - l;
        }

        d2 = xpx[kk] - d;
	test = d2 / xpx[kk];

	/* check for singularity */
        if (test < TINY) {
	    fprintf(stderr, "cholbeta: test[%d] = %g\n", j, test);
	    *rss = -1.0;
	    return E_SINGULAR;
        } else if (test < SMALL) {
	    gretl_model_set_int(pmod, "near-singular", 1);
	} 

        e = 1 / sqrt(d2);
        xpx[kk] = e;
        xpy[j] = (xpy[j] - d1) * e;

	/* off-diagonal elements */
        for (i=j+1; i<nc; i++) {
            kk++;
            d = 0.0;
            k = j;
            for (l=1; l<=jm1; l++) {
                d += xpx[k] * xpx[k-j+i];
                k += nc - l;
            }
            xpx[kk] = (xpx[kk] - d) * e;
        }
        kk++;
    }

    kk--;

    /* calculate regression sum of squares */
    d = 0.0;
    for (j=0; j<nc; j++) {
	d += xpy[j] * xpy[j];
    }
    *rss = d;

    /* back-solve for the coefficients */

    coeff[nc-1] = xpy[nc-1] * xpx[kk];

    for (j=nc-2; j>=0; j--) {
	d = xpy[j];
	for (i=nc-1; i>j; i--) {
	    d -= coeff[i] * xpx[--kk];
	}
	coeff[j] = d * xpx[--kk];
    }

    for (j=0; j<nc; j++) {
	if (isnan(coeff[j])) {
	    fprintf(stderr, "%s %d: coeff %d is NaN\n", __FILE__, __LINE__, j);
	    return E_NAN;
	}
    }

    return 0; 
}

/*
  regress: takes X'X (pmod->xpx) and X'y (@xpy) and
  computes OLS estimates and associated statistics.

  n = total number of observations per series in data set
  pmod->ifc = 1 if constant is present in model, else = 0

  ess = error sum of squares
  sigma = standard error of regression
  fstt = F-statistic
  pmod->coeff = array of regression coefficients
  pmod->sderr = corresponding array of standard errors
*/

static void regress (MODEL *pmod, double *xpy, 
		     double ysum, double ypy,
		     const double **Z, 
		     double rho, gretlopt opt)
{
    int yno = pmod->list[1];
    int ifc = pmod->ifc;
    int n = pmod->full_n;
    double zz, rss = 0.0;
    double s2 = 0.0;
    double *diag = NULL;
    int i, err = 0;

    for (i=0; i<n; i++) {
	pmod->yhat[i] = pmod->yhat[i] = NADBL;
    }    

    zz = ysum * ysum / pmod->nobs;
    pmod->tss = ypy - zz;

    /*  Cholesky-decompose X'X and find the coefficients */
    err = cholbeta(pmod, xpy, &rss);
    if (err) {
        pmod->errcode = err;
        return;
    }   
    
    if (rho != 0.0) {
	pmod->ess = ypy - rss;
    } else {
	make_ess(pmod, Z);
	rss = ypy - pmod->ess;
    }

    if (fabs(pmod->ess) < ESSZERO) {
	pmod->ess = 0.0;
    } else if (pmod->ess < 0.0) { 
	sprintf(gretl_errmsg, _("Error sum of squares (%g) is not > 0"),
		pmod->ess);
        return; 
    }

    if (pmod->dfd == 0) {
	pmod->sigma = 0.0;
	pmod->adjrsq = NADBL;
    } else {
	if (pmod->opt & OPT_N) {
	    /* no-df-corr */
	    s2 = pmod->ess / pmod->nobs;
	} else {
	    s2 = pmod->ess / pmod->dfd;
	}
	pmod->sigma = sqrt(s2);
    }

    if (pmod->tss < DBL_EPSILON) {
	pmod->rsq = pmod->adjrsq = NADBL;
    } 

    hatvar(pmod, n, Z); 
    if (pmod->errcode) return;

    if (!(opt & OPT_B)) {
	/* changed 2008-09-25 */
	if (pmod->tss > 0.0) {
	    compute_r_squared(pmod, Z[yno], &ifc);
	} else if (pmod->tss == 0.0) {
	    uncentered_r_squared(pmod, Z[yno]);
	}
    }

    if (s2 <= 0.0 || pmod->dfd == 0 || pmod->dfn == 0) {
	pmod->fstt = NADBL;
    } else if (pmod->rsq == 1.0) {
	pmod->fstt = NADBL;
    } else {
	pmod->fstt = (rss - zz * ifc) / (s2 * pmod->dfn);
	if (pmod->fstt < 0.0) {
	    pmod->fstt = 0.0;
	}
    }

    diag = malloc(pmod->ncoeff * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    /* Note: 'xpy' is used purely as workspace in diaginv() below, as
       a matter on economy/convenience; no values in this array are
       used before they are written. 
    */

    diaginv(pmod->xpx, xpy, diag, pmod->ncoeff);

    for (i=0; i<pmod->ncoeff; i++) {
	if (diag[i] >= 0.0) {
	    pmod->sderr[i] = pmod->sigma * sqrt(diag[i]);
	} else {
	    pmod->sderr[i] = 0.0;
	}
    }

    free(diag); 
}

/**
 * makevcv:
 * @pmod: pointer to model.
 * @sigma: square root of error variance, or 1.0 to
 * produce just X'X^{-1}.
 *
 * Inverts the Cholesky-decomposed X'X and computes the 
 * coefficient covariance matrix.
 * 
 * Returns: 0 on successful completion, non-zero code on error.
 */

int makevcv (MODEL *pmod, double sigma)
{
    int dec, mst, kk, i, j, kj, icnt, m, k, l = 0;
    int nv, nxpx;
    double d;

    if (pmod->vcv != NULL) {
	/* already done */
	return 0;
    }

    if (pmod->xpx == NULL) {
	/* raw material not available */
	fprintf(stderr, "makevcv: pmod->xpx = NULL\n");
	return E_BADSTAT;
    }

    nv = pmod->ncoeff;
    nxpx = (nv * nv + nv) / 2; 
    mst = nxpx;
    kk = nxpx - 1;

    pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nv; i++) {
	mst -= i;
	/* find diagonal element */
	d = pmod->xpx[kk];
	if (i > 0) {
	    for (j=kk+1; j<=kk+i; j++) {
		d -= pmod->xpx[j] * pmod->vcv[j];
	    }
	}
	pmod->vcv[kk] = d * pmod->xpx[kk];
	/* find off-diagonal elements indexed by kj */
	kj = kk;
	kk = kk - i - 2;
	if (i > nv - 2) {
	    continue;
	}
	for (j=i+1; j<nv; j++) {
	    icnt = i+1;
	    kj -= j;
	    d = 0.0;
	    m = mst + 1;
	    for (k=0; k<=j-1; k++) {
		if (icnt > 0) {
		    dec = 1;
		    icnt--;
		} else {
		    dec = k;
		}
		m -= dec;
		l = kj + i - k;
		d += pmod->vcv[m-1] * pmod->xpx[l];
	    }
	    pmod->vcv[kj] = -d * pmod->xpx[l-1];
	}
    }

    if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	sigma = 1.0;
    }

    if (sigma != 1.0) {
	double s2 = sigma * sigma;

	for (i=0; i<nxpx; i++) {
	    pmod->vcv[i] *= s2;
	}
    }

    return 0;
}

/**
 * dwstat:
 * @order: order of autoregression (usually 1).
 * @pmod: pointer to model.
 * @Z: data array.
 *
 * Computes the Durbin-Watson statistic for @pmod.
 * 
 * Returns: the D-W value, or #NADBL on error.
 */

double dwstat (int order, MODEL *pmod, const double **Z)
{
    double ut, u1;
    double num = 0.0;
    double den = 0.0;
    int t, t1;

    if (pmod->ess <= 0.0) {
	return NADBL;
    }

    t1 = pmod->t1 + order;

    if (pmod->nwt) {
	ut = pmod->uhat[t1 - 1];
	if (!na(ut)) {
	    den += ut * ut;
	}
    } else {
	den = pmod->ess;
    }

    for (t=t1; t<=pmod->t2; t++)  {
        ut = pmod->uhat[t];
        u1 = pmod->uhat[t-1];
        if (na(ut) || na(u1) ||
	    (pmod->nwt && (Z[pmod->nwt][t] == 0.0 || 
			   Z[pmod->nwt][t-1] == 0.0))) { 
	    continue;
	}
        num += (ut - u1) * (ut - u1);
	if (pmod->nwt) {
	    den += ut * ut;
	}
    }

    return num / den;
}

/* altrho: alternative calculation of rho */

static double altrho (int order, int t1, int t2, const double *uhat)
{
    double *ut, *u1;    
    int t, n, len = t2 - (t1 + order) + 1;
    double uht, uh1, rho;

    ut = malloc(len * sizeof *ut);
    if (ut == NULL) {
	return NADBL;
    }

    u1 = malloc(len * sizeof *u1);
    if (u1 == NULL) {
	free(ut);
	return NADBL;
    }

    n = 0;

    for (t=t1+order; t<=t2; t++) { 
        uht = uhat[t];
	uh1 = (t > 0)? uhat[t-1] : NADBL;
        if (!na(uht) && !na(uh1)) {
	    ut[n] = uht;
	    u1[n] = uh1;
	    n++;
	}
    }

    rho = gretl_corr(0, n - 1, ut, u1, NULL);

    free(ut);
    free(u1);

    return rho;
}

/**
 * rhohat:
 * @order: order of autoregression, usually 1.
 * @t1: start of sample range.
 * @t2: end of sample range.
 * @uhat: array of regression residuals.
 *
 * Computes the first order serial correlation coefficient
 * for @uhat, over the range @t1 to @t2.
 * 
 * Returns: the \hat{rho} value, or #NADBL on error.
 */

double rhohat (int order, int t1, int t2, const double *uhat)
{
    double ut, u1, uu = 0.0, xx = 0.0;
    double rho;
    int t;

    for (t=t1+order; t<=t2; t++) { 
        ut = uhat[t];
        u1 = uhat[t-1];
        if (na(ut) || na(u1)) {
	    continue;
	}
        uu += ut * u1;
        xx += u1 * u1;
    }

    if (xx < DBL_EPSILON) {
	return NADBL;
    }

    rho = uu / xx;

    if (rho > 1.0 || rho < -1.0) {
	rho = altrho(order, t1, t2, uhat);
    }

    return rho;
}

/* compute fitted values and residuals */

static int hatvar (MODEL *pmod, int n, const double **Z)
{
    int xno, i, t;
    int yno = pmod->list[1];
    double x;

    for (t=0; t<n; t++) {
	pmod->yhat[t] = pmod->uhat[t] = NADBL;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}
	pmod->yhat[t] = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            xno = pmod->list[i+2];
	    x = Z[xno][t];
	    if (pmod->nwt) {
		x *= sqrt(Z[pmod->nwt][t]);
	    }
            pmod->yhat[t] += pmod->coeff[i] * x;
        }
	x = Z[yno][t];
	if (pmod->nwt) {
	    x *= sqrt(Z[pmod->nwt][t]);
	}
        pmod->uhat[t] = x - pmod->yhat[t];                
    }

    return 0;
}

static int hilu_plot (double *ssr, double *rho, int n)
{
    FILE *fp;
    int i;

    if (gnuplot_init(PLOT_REGULAR, &fp)) {
	return E_FOPEN; 
    }

    fputs("# hildreth-lu\n", fp);
    fputs("set xlabel 'rho'\n", fp);

    fprintf(fp, "set ylabel '%s'\n", I_("ESS"));

    fputs("set nokey\n", fp);
    fputs("set xrange [-1.0:1.0]\n", fp);
    fputs("plot '-' using 1:2 w impulses\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	fprintf(fp, "%g %g\n", rho[i], ssr[i]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    gnuplot_make_graph();

    return 0;
}

static double autores (MODEL *pmod, const double **Z, gretlopt opt)
{
    int t, v, t1 = pmod->t1;
    double x, num = 0.0, den = 0.0;
    double rhohat;

    if (!(opt & OPT_P)) {
	/* not using Prais-Winsten */
	t1--;
    }

    for (t=t1; t<=pmod->t2; t++) {
	x = Z[pmod->list[1]][t];
	for (v=0; v<pmod->ncoeff; v++) {
	    x -= pmod->coeff[v] * Z[pmod->list[v+2]][t];
	}
	pmod->uhat[t] = x;
	if (t > t1) {
	    num += pmod->uhat[t] * pmod->uhat[t-1];
	    den += pmod->uhat[t-1] * pmod->uhat[t-1];
	}
    } 

    rhohat = num / den;

    return rhohat;
}

#define AR_DEBUG 0

/**
 * estimate_rho:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: option flags: may include %OPT_H to use Hildreth-Lu,
 *       %OPT_P to use Prais-Winsten, %OPT_B to suppress Cochrane-Orcutt
 *       fine-tuning of Hildreth-Lu results, %OPT_G to generate
 *       a gnuplot graph of the search in Hildreth-Lu case.
 * @prn: gretl printing struct.
 * @err: location to receeve error code.
 *
 * Estimate the quasi-differencing coefficient for use with the
 * Cochrane-Orcutt, Hildreth-Lu or Prais-Winsten procedures for
 * handling first-order serial correlation.  Print a trace of the
 * search for rho.
 * 
 * Returns: rho estimate on successful completion, #NADBL on error.
 */

double estimate_rho (const int *list, double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn, int *err)
{
    double rho = 0.0, rho0 = 0.0, diff;
    double finalrho = 0.0, essmin = 1.0e8;
    double ess, ssr[199], rh[199]; 
    int iter, nn = 0;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int missv = 0, misst = 0;
    gretlopt lsqopt = OPT_A;
    int quiet = (opt & OPT_Q);
    int ascii = !(opt & OPT_G);
    MODEL armod;

    gretl_error_clear();
    *err = 0;

    missv = adjust_t1t2(NULL, list, &pdinfo->t1, &pdinfo->t2, 
			pdinfo->n, (const double **) *pZ, &misst);
    if (missv) {
	sprintf(gretl_errmsg, _("Missing value encountered for "
				"variable %d, obs %d"), missv, misst);
	*err = E_DATA;
	goto bailout;
    }

    gretl_model_init(&armod);

    if (opt & OPT_P) {
	/* Prais-Winsten treatment of first observation */
	lsqopt |= OPT_P;
    } 

    if (opt & OPT_H) { 
	/* Do Hildreth-Lu first */
	for (rho = -0.990, iter = 0; rho < 1.0; rho += .01, iter++) {
	    clear_model(&armod);
	    armod = ar1_lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
	    if ((*err = armod.errcode)) {
		clear_model(&armod);
		goto bailout;
	    }
	    ess = armod.ess;
	    if (ascii && !quiet) {
		char num[16];
		int chk;
		
		if (iter == 0) {
		    pprintf(prn, "\n RHO       %s      RHO       %s      "
			    "RHO       %s      RHO       %s     \n",
			    _("ESS"), _("ESS"), _("ESS"), _("ESS"));
		}
		sprintf(num, "%f", 100 * fabs(rho));
		chk = atoi(num);
		if (chk == 99 || chk % 10 == 0) {
		    ssr[nn] = ess;
		    rh[nn++] = rho;
		    pprintf(prn, "%5.2f %10.4g", rho, ess);
		    if (nn % 4 == 0) {
			pputc(prn, '\n');
		    } else {
			bufspace(3, prn);
		    }
		} 
	    } else if (!quiet) {
		ssr[nn] = ess;
		rh[nn++] = rho;
	    }
	    if (iter == 0 || ess < essmin) {
		essmin = ess;
		finalrho = rho;
	    }
	} /* end of basic iteration */
	
	if (finalrho > 0.989) {
	    /* try exploring this funny region? */
	    for (rho = 0.99; rho <= 0.999; rho += .001) {
		clear_model(&armod);
		armod = ar1_lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
		if ((*err = armod.errcode)) {
		    clear_model(&armod);
		    goto bailout;
		}
		ess = armod.ess;
		if (ess < essmin) {
		    essmin = ess;
		    finalrho = rho;
		}
	    }
	}

	if (finalrho > 0.9989) {
	    /* this even funnier one? */
	    for (rho = 0.9991; rho <= 0.9999; rho += .0001) {
		clear_model(&armod);
		armod = ar1_lsq(list, pZ, pdinfo, OLS, OPT_A, rho);
		if ((*err = armod.errcode)) {
		    clear_model(&armod);
		    goto bailout;
		}
		ess = armod.ess;
		if (ess < essmin) {
		    essmin = ess;
		    finalrho = rho;
		}
	    }
	}

	rho0 = rho = finalrho;
	if (!quiet) {
	    pprintf(prn, _("\n\nESS is minimum for rho = %g\n\n"), rho);
	    if (ascii) {
		graphyx(ssr, rh, nn, "ESS", "RHO", prn); 
		pputs(prn, "\n");
	    } else {
		hilu_plot(ssr, rh, nn);
	    }
	}
    } else { 
	/* Go straight to Cochrane-Orcutt (or Prais-Winsten) */
	armod = lsq(list, pZ, pdinfo, OLS, OPT_A);
	if (!armod.errcode && armod.dfd == 0) {
	    armod.errcode = E_DF;
	}
	if ((*err = armod.errcode)) {
	    clear_model(&armod);
	    goto bailout;
	}
	rho0 = rho = armod.rho;
    }

    if (na(rho)) {
	*err = E_NOCONV;
	clear_model(&armod);
	goto bailout;
    }

    if (!(opt & OPT_H) || !(opt & OPT_B)) {

	if (!quiet) {
	    if (opt & OPT_H) {
		pputs(prn, _("\nFine-tune rho using the CORC procedure...\n\n"));
	    } else {
		pputs(prn, _("\nPerforming iterative calculation of rho...\n\n"));
	    }

	    pputs(prn, _("                 ITER       RHO        ESS"));
	    pputc(prn, '\n');
	}

	iter = 0;
	diff = 1.0;

	while (diff > 0.001) {
	    clear_model(&armod);
	    armod = ar1_lsq(list, pZ, pdinfo, OLS, lsqopt, rho);
#if AR_DEBUG
	    fprintf(stderr, "armod: t1=%d, first two uhats: %g, %g\n",
		    armod.t1, 
		    armod.uhat[armod.t1],
		    armod.uhat[armod.t1+1]);
#endif
	    if ((*err = armod.errcode)) {
		clear_model(&armod);
		goto bailout;
	    }
	    if (!quiet) {
		pprintf(prn, "          %10d %12.5f", ++iter, rho);
		pprintf(prn, "   %g\n", armod.ess);
	    }

	    rho = autores(&armod, (const double **) *pZ, opt);

#if AR_DEBUG
	    pputs(prn, "AR1 model (using rho-transformed data)\n");
	    printmodel(&armod, pdinfo, OPT_NONE, prn);
	    pprintf(prn, "autores gives rho = %g\n", rho);
#endif

	    if (rho > .99999 || rho < -.99999) {
		*err = E_NOCONV;
		clear_model(&armod);
		goto bailout;
	    }

	    diff = (rho > rho0)? rho - rho0 : rho0 - rho;
	    rho0 = rho;
	    if (iter == 30) break;
	}

	if (!quiet) {
	    pprintf(prn, "          %10d %12.5f", ++iter, rho);
	    pprintf(prn, "   %g\n", armod.ess);
	}
    }

    clear_model(&armod);

 bailout:

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    if (*err) {
	rho = NADBL;
    }

    return rho;
}

/**
 * augment_regression_list:
 * @orig: list giving original regression specification.
 * @aux: either %AUX_SQ, %AUX_LOG or %AUX_WHITE.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 *
 * Augment the regression list @orig with auxiliary terms.  If @aux 
 * is %AUX_SQ add the squares of the original regressors; if @aux
 * is %AUX_WHITE add squares and cross-products, or if @aux is
 * %AUX_LOG add the natural logs of the original regressors.
 * If they are not already present, these variables are added
 * to the data array.
 * 
 * Returns: the augmented list, or NULL on failure.
 */

int *augment_regression_list (const int *orig, int aux, 
			      double ***pZ, DATAINFO *pdinfo)
{
    int *list;
    int listlen;
    int cnum = 0;
    int i, k;

    if (aux == AUX_WHITE) {
	int cpos = gretl_list_const_pos(orig, 2, (const double **) *pZ, 
					pdinfo);
	int nt, trv = orig[0] - 1;

	if (cpos > 0) {
	    trv--;
	    cnum = orig[cpos];
	}
	nt = (trv * trv + trv) / 2;
	listlen = orig[0] + nt + 1;
    } else {
	listlen = 2 * orig[0];
    }

    list = malloc(listlen * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    /* transcribe original list */
    for (i=0; i<=orig[0]; i++) {
	list[i] = orig[i];
    }

    /* add squares, cross-products or logs of independent vars */
    k = list[0];
    for (i=2; i<=orig[0]; i++) {
	int vnew, vi = orig[i];

	if (vi == 0) {
	    continue;
	}

	if (aux == AUX_SQ || aux == AUX_WHITE) {
	    vnew = xpxgenr(vi, vi, pZ, pdinfo);
	    if (vnew > 0) {
		list[++k] = vnew;
	    }
	    if (aux == AUX_WHITE) {
		int j, vj;

		for (j=i+1; j<=orig[0]; j++) {
		    vj = orig[j];
		    if (vj == cnum) {
			continue;
		    }
		    vnew = xpxgenr(vi, vj, pZ, pdinfo);
		    if (vnew > 0) {
			/* ensure uniqueness of generated varnames */
			sprintf(pdinfo->varname[vnew], "X%d_X%d", i-1, j-1);
			list[++k] = vnew;
		    }
		}
	    }
	} else if (aux == AUX_LOG) {
	    vnew = loggenr(vi, pZ, pdinfo);
	    if (vnew > 0) {
		list[++k] = vnew;
	    }
	}
    }

    list[0] = k;

    return list;
}

/* For observation s, see if the regression list contains a variable
   that has a single non-zero value at that particular observation.
   We run this check on observations showing an OLS residual of zero.
*/

static int observation_is_dummied (const MODEL *pmod,
				   int *list, const double **Z,
				   int s)
{
    int i, t, v;
    int ret = 0;

    for (i=list[0]; i>=2; i--) {
	v = list[i];
	if (v == 0) {
	    continue;
	}
	ret = 1;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if ((t == s && Z[v][t] == 0.0) || (t != s && Z[v][t] != 0.0)) {
		ret = 0;
		break;
	    }
	}
	if (ret) {
	    gretl_list_delete_at_pos(list, i);
	    break;
	}
    }

    return ret;
}

/* get_hsk_weights: take the residuals from the model pmod, square them
   and take logs; find the fitted values for this series using an
   auxiliary regression including the original independent variables
   and their squares; exponentiate the fitted values; and add the
   resulting series to the data set.
*/

static int get_hsk_weights (MODEL *pmod, double ***pZ, DATAINFO *pdinfo)
{
    int oldv = pdinfo->v;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int *lcpy = NULL;
    int *list = NULL;
    int err = 0, shrink = 0;
    double xx;
    MODEL aux;

    lcpy = gretl_list_copy(pmod->list);
    if (lcpy == NULL) {
	return E_ALLOC;
    }

    /* allocate space for an additional variable */
    if (dataset_add_series(1, pZ, pdinfo)) {
	free(lcpy);
	return E_ALLOC;
    }

    /* add transformed pmod residuals to data set */
    for (t=0; t<pdinfo->n; t++) {
	xx = pmod->uhat[t];
	if (na(xx)) {
	    (*pZ)[oldv][t] = NADBL;
	} else if (xx == 0.0) {
	    if (observation_is_dummied(pmod, lcpy, (const double **) *pZ, t)) {
		(*pZ)[oldv][t] = NADBL;
	    } else {
		fprintf(stderr, "hsk: got a zero residual, could be a problem!\n");
		(*pZ)[oldv][t] = -1.0e16; /* ?? */
	    }
	} else {
	    (*pZ)[oldv][t] = log(xx * xx);
	}
    }

    /* build regression list, adding the squares of the original
       independent vars */
    list = augment_regression_list(lcpy, AUX_SQ, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    list[1] = oldv; /* the newly added log(uhat-squared) */

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    aux = lsq(list, pZ, pdinfo, OLS, OPT_A);
    err = aux.errcode;
    if (err) {
	shrink = pdinfo->v - oldv;
    } else {
	/* write into the data set the required weights */
	for (t=aux.t1; t<=aux.t2; t++) {
	    xx = aux.yhat[t];
	    if (na(xx)) {
		(*pZ)[oldv][t] = NADBL;
	    } else {
		(*pZ)[oldv][t] = 1.0 / exp(xx);
	    }
	}
	shrink = pdinfo->v - oldv - 1;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    clear_model(&aux);

    if (shrink > 0) {
	dataset_drop_last_variables(shrink, pZ, pdinfo);
    }

    free(list);
    free(lcpy);

    return err;
}

/**
 * hsk_func:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using a correction for
 * heteroskedasticity.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL hsk_func (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    int i, err;
    int orig_nvar = pdinfo->v;
    int *hsklist;
    MODEL hsk;

    gretl_error_clear();

    /* run initial OLS */
    hsk = lsq(list, pZ, pdinfo, OLS, OPT_A);
    if (hsk.errcode) {
	return hsk;
    }

    /* use the residuals from the initial OLS to form weights */
    err = get_hsk_weights(&hsk, pZ, pdinfo);
    if (err) {
	hsk.errcode = err;
	return hsk;
    }

    /* allocate regression list for weighted least squares */
    hsklist = gretl_list_new(list[0] + 1);
    if (hsklist == NULL) {
	hsk.errcode = E_ALLOC;
	return hsk;
    }

    /* the last variable in the dataset will be the weight var */
    hsklist[1] = pdinfo->v - 1;

    /* put the original dependent variable in at position 2 */
    hsklist[2] = list[1];

    /* add the original independent vars into the WLS list */
    for (i=3; i<=hsklist[0]; i++) {
	hsklist[i] = list[i-1];
    }

    clear_model(&hsk);
    hsk = lsq(hsklist, pZ, pdinfo, WLS, OPT_NONE);
    hsk.ci = HSK;

    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    free(hsklist);

    return hsk;
}

static double **allocate_hccm_p (int k, int n)
{
    double **p = malloc(k * sizeof *p);
    int i;

    if (p == NULL) return NULL;

    for (i=0; i<k; i++) {
	p[i] = malloc(n * sizeof **p);
	if (p[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) {
		free(p[j]);
	    }
	    free(p);
	    p = NULL;
	    break;
	}
    }

    return p;
}

static void free_hccm_p (double **p, int m)
{
    if (p != NULL) {
	int i;

	for (i=0; i<m; i++) {
	    free(p[i]);
	}
	free(p);
    }
}

static int jackknife_vcv (MODEL *pmod, const double **Z)
{
    double *st = NULL, *ustar = NULL;
    double **p = NULL;
    int nobs, tp, nc = 0;
    int i, j, k, t;
    int t1, t2;
    double xx;
    int err = 0;

    gretl_error_clear();

    t1 = pmod->t1;
    t2 = pmod->t2;
    nobs = pmod->nobs;
    nc = pmod->ncoeff;

    st = malloc(nc * sizeof *st);
    ustar = malloc(nobs * sizeof *ustar);
    p = allocate_hccm_p(nc, nobs);

    if (st == NULL || p == NULL || ustar == NULL) {
	err = E_ALLOC;
	goto bailout;
    }  

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
	pmod->vcv = NULL;
    }

    if (makevcv(pmod, 1.0)) {
	err = E_ALLOC;
	goto bailout;
    }

    /* form elements of (X'X)^{-1}X' */

    for (i=0; i<nc; i++) {
	tp = 0;
	for (t=t1; t<=t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    xx = 0.0;
	    for (j=0; j<nc; j++) {
		if (i <= j) {
		    k = ijton(i, j, nc);
		} else {
		    k = ijton(j, i, nc);
		}
		xx += pmod->vcv[k] * Z[pmod->list[j+2]][t];
	    }
	    p[i][tp++] = xx;
	}
    }

    tp = 0;
    for (t=t1; t<=t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}	
	xx = 0.0;
	for (i=0; i<nc; i++) {
	    xx += Z[pmod->list[i+2]][t] * p[i][tp];
	}
	if (floateq(xx, 1.0)) {
	    xx = 0.0;
	}
	ustar[tp++] = pmod->uhat[t] / (1.0 - xx);
    }

    for (i=0; i<nc; i++) {
	xx = 0.0;
	for (t=0; t<nobs; t++) {
	    xx += p[i][t] * ustar[t];
	}
	st[i] = xx;
    }

    for (t=0; t<nobs; t++) {
	for (i=0; i<nc; i++) {
	    p[i][t] *= ustar[t];
	}
    }

    /* MacKinnon and White, 1985, equation (13) */

    k = 0;
    for (i=0; i<nc; i++) {
	for (j=i; j<nc; j++) {
	    xx = 0.0;
	    for (t=0; t<nobs; t++) {
		xx += p[i][t] * p[j][t];
	    }
	    xx -= st[i] * st[j] / nobs;
	    /* MacKinnon and White: "It is tempting to omit the factor
	       (n - 1) / n from HC3" (1985, p. 309).  Here we leave it in
	       place, as in their simulations.
	    */
	    xx *= (nobs - 1.0) / nobs;
	    if (i == j) {
		pmod->sderr[i] = sqrt(xx);
	    }
	    pmod->vcv[k++] = xx;
	}
    }

    /* substitute robust F stat */
    if (pmod->dfd > 0 && pmod->dfn > 1) {
	pmod->fstt = wald_omit_F(NULL, pmod);
    }

    pmod->opt |= (OPT_R | OPT_J);
    gretl_model_set_vcv_info(pmod, VCV_HC, 4);

 bailout:

    free(st);
    free(ustar);
    free_hccm_p(p, nc);

    return err;
}

static void print_HET_1 (double z, double pval, PRN *prn)
{
    pprintf(prn, "\n%s\n", _("Pesaran-Taylor test for heteroskedasticity"));
    pprintf(prn, "\n%s: HET_1 = %f,\n", _("Test statistic"), z);
    pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n", 
	    _("with p-value"), z, pval);
}

static void print_whites_test (double LM, int df, double pval, 
			       gretlopt opt, PRN *prn)
{
    if (opt & OPT_B) {
	pprintf(prn, "\n%s\n", _("Breusch-Pagan test for heteroskedasticity"));
	if (opt & OPT_R) {
	    pprintf(prn, "%s\n", _("with Koenker robust variance estimator"));
	}
	pprintf(prn, "\n%s: LM = %f,\n", _("Test statistic"), LM);
    } else {
	if (opt & OPT_X) {
	    pprintf(prn, "\n%s\n", _("White's test for heteroskedasticity (no cross products)"));
	} else {
	    pprintf(prn, "\n%s\n", _("White's test for heteroskedasticity"));
	}
	pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic"), LM);
    }

    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
	    _("with p-value"), _("Chi-square"), 
	    df, LM, pval); 
}

/* For White's test, see if we have enough degrees of freedom
   to add squares -- and if so, if we have enough df to add
   cross-products also.
*/

static int get_whites_aux (const MODEL *pmod, const double **Z) 
{
    int aux = AUX_WHITE;
    int rem = pmod->ncoeff - pmod->ifc - 1;
    int nsq = 0, nxpx = 0;
    int i, v;

    for (i=2; i<=pmod->list[0]; i++) {
	v = pmod->list[i];
	if (v > 0) {
	    if (!gretl_isdummy(pmod->t1, pmod->t2, Z[v])) {
		nsq++;
	    }
	    nxpx += rem--;
	}
    }

    if (pmod->dfd - nsq < 1) {
	aux = AUX_NONE;
    } else if (pmod->dfd - nsq - nxpx < 1) {
	aux = AUX_SQ;
    }

    return aux;
}

#define PT_DEBUG 0

/**
 * tsls_hetero_test:
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save results to model; %OPT_Q
 * means don't print the auxiliary regression.
 * @prn: gretl printing struct.
 *
 * Runs Pesaran and Taylor's (1999) HET_1 test for heteroskedasticity
 * on the given tsls model. The statistic is just a t-statistic, so
 * under the null it is distributed as a standard normal.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

static int tsls_hetero_test (MODEL *pmod, double ***pZ, 
			     DATAINFO *pdinfo, gretlopt opt, 
			     PRN *prn)
{
    int pos, v = pmod->list[1], newv = pdinfo->v;
    int *auxlist = NULL, *testlist = NULL;
    int i, h, t;
    int savet1 = pdinfo->t1;
    int savet2 = pdinfo->t2;
    MODEL ptmod;
    double x;
    int err = 0;

    if (pmod->opt & (OPT_G | OPT_L)) {
	/* FIXME gmm, liml */
	return E_NOTIMP;
    }

    pos = gretl_list_separator_position(pmod->list);
    h = pmod->list[0] - pos;

#if PT_DEBUG
    pprintf(prn, "v = %d, h = %d\n", v, h);
#endif

    auxlist = gretl_list_new(h + 1);
    testlist = gretl_list_new(3);

    if (auxlist == NULL || testlist == NULL) {
	free(auxlist);
	free(testlist);
	return E_ALLOC;
    }

    auxlist[1] = v;
    for (i=2; i<=auxlist[0]; i++) {
	auxlist[i] = pmod->list[i + pos - 1];
    }	

    testlist[1] = newv;
    testlist[2] = 0;
    testlist[3] = newv + 1;

#if PT_DEBUG
    printlist(auxlist, "auxlist");
    printlist(testlist, "testlist");
#endif

    ptmod = lsq(auxlist, pZ, pdinfo, OLS, OPT_A);
    err = ptmod.errcode;
    if (err) {
	goto bailout;
    }

#if PT_DEBUG
    printmodel(&ptmod, pdinfo, OPT_S, prn);
#endif

    err = dataset_add_series(2, pZ, pdinfo);
    if (err) {
	clear_model(&ptmod);
	goto bailout;
    }

    strcpy(pdinfo->varname[newv+1], "yhat^2");

    for (t=pmod->t1; t<=pmod->t2; t++) {
	x = pmod->uhat[t];
	(*pZ)[newv][t] = x*x;
	x = ptmod.yhat[t];
	(*pZ)[newv+1][t] = x*x;
    }

    clear_model(&ptmod);

    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    ptmod = lsq(testlist, pZ, pdinfo, OLS, OPT_A);
    err = ptmod.errcode;

    if (!err) {
	double z = fabs(ptmod.coeff[1]) / ptmod.sderr[1];
	double pval = 2.0 * (1 - normal_cdf(z));

	if (opt & OPT_Q) {
	    print_HET_1(z, pval, prn);
	} else {
	    ptmod.aux = AUX_HET_1;
	    printmodel(&ptmod, pdinfo, OPT_NONE, prn);
	}

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_HET_1);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_Z);
		model_test_set_value(test, z);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	  
	}

	record_test_result(z, pval, _("HET_1"));
    }

    clear_model(&ptmod);

    dataset_drop_last_variables(2, pZ, pdinfo); 

 bailout:

    free(auxlist);
    free(testlist);

    pdinfo->t1 = savet1;
    pdinfo->t2 = savet2;

    return err;
}

/* Compute LM statistic as per Breusch and Pagan (Econometrica, 1979),
   with the option to use the robust variance estimator proposed by
   Koenker (Journal of Econometrics, 1981).
*/

static double get_BP_LM (MODEL *pmod, int *list, MODEL *aux,
			 double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, int *err)
{
    double s2, u2t, gt;
    double V = 0.0, LM = NADBL;
    int t, v = list[1];

    s2 = pmod->ess / pmod->nobs;

    if (opt & OPT_R) {
	/* calculate robust variance estimate a la Koenker */
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    u2t = pmod->uhat[t] * pmod->uhat[t];
	    V += (u2t - s2) * (u2t - s2);
	}
	V /= pmod->nobs;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	u2t = pmod->uhat[t] * pmod->uhat[t];
	if (opt & OPT_R) {
	    gt = u2t - s2;
	} else {
	    gt = u2t / s2;
	}
	(*pZ)[v][t] = gt;
    }

    *aux = lsq(list, pZ, pdinfo, OLS, OPT_A);
    *err = aux->errcode;

    if (!*err) {
	double RSS = aux->tss - aux->ess;

	if (RSS < 0) {
	    *err = E_DATA;
	} else {
	    if (opt & OPT_R) {
		aux->opt |= OPT_R;
		LM = RSS / V;
	    } else {
		LM = .5 * RSS;
	    }
	    gretl_model_set_double(aux, "BPLM", LM);
	    aux->aux = AUX_BP;
	}
    }

    return LM;
}

/**
 * whites_test:
 * @pmod: pointer to model.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save results to model; %OPT_Q
 * means don't print the auxiliary regression;  %OPT_B means
 * do the simpler Breusch-Pagan variant.
 * @prn: gretl printing struct.
 *
 * Runs White's test for heteroskedasticity on the given model.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int whites_test (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn)
{
    int lo, ncoeff, yno, t;
    int BP = (opt & OPT_B);
    int aux = AUX_NONE;
    int v = pdinfo->v;
    int *list = NULL;
    double zz, LM;
    MODEL white;
    int err = 0;

    if (pmod->ci == IVREG) {
	return tsls_hetero_test(pmod, pZ, pdinfo, opt, prn);
    }

    if (pmod->list == NULL || gretl_list_has_separator(pmod->list)) {
	return E_NOTIMP;
    }

    if (pmod->ci == NLS || pmod->ci == MLE ||
	pmod->ci == GMM || pmod->ci == ARMA || 
	pmod->ci == LOGISTIC) { 
	return E_NOTIMP;
    }

    if ((err = list_members_replaced(pmod->list, pdinfo, pmod->ID))) {
	return err;
    }

    /* what can we do, with the degrees of freedom available? */
    if (!BP) {
	if (opt & OPT_X) { 
	    aux = AUX_SQ;
	} else {
	    aux = get_whites_aux(pmod, (const double **) *pZ);
	    if (aux == AUX_NONE) {
		return E_DF;
	    }
	}
    }

    gretl_model_init(&white);

    lo = pmod->list[0];
    yno = pmod->list[1];
    ncoeff = pmod->ncoeff;    

    /* make space in data set */
    if (dataset_add_series(1, pZ, pdinfo)) {
	err = E_ALLOC;
    }

    if (!err) {
	/* get residuals, square and add to data set */
	for (t=0; t<pdinfo->n; t++) {
	    zz = pmod->uhat[t];
	    if (na(zz)) {
		(*pZ)[v][t] = NADBL;
	    } else {
		(*pZ)[v][t] = zz * zz;
	    }
	}
	strcpy(pdinfo->varname[v], "uhatsq");
    }

    if (!err) {
	if (BP) {
	    list = gretl_list_copy(pmod->list);
	} else {
	    /* build aux regression list, adding squares and
	       (possibly) cross-products of the original 
	       independent vars */
	    list = augment_regression_list(pmod->list, aux, pZ, pdinfo);
	}
	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    list[1] = v; /* the newly added uhat-squared */
	}
    }

    if (!err) {
	/* run auxiliary regression */
	if (BP) {
	    LM = get_BP_LM(pmod, list, &white, pZ, pdinfo, opt, &err);
	} else {
	    white = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_Q);
	    err = white.errcode;
	    if (!err) {
		LM = white.rsq * white.nobs;
		white.aux = AUX_WHITE;
	    }
	}
    }

    if (!err) {
	int df = white.ncoeff - 1;
	double pval = chisq_cdf_comp(df, LM);

	if (opt & OPT_Q) {
	    print_whites_test(LM, df, pval, opt, prn);
	} else {
	    printmodel(&white, pdinfo, OPT_NONE, prn);
	}

	if (opt & OPT_S) {
	    ModelTestType mt = BP? GRETL_TEST_BP : GRETL_TEST_WHITES;
	    ModelTest *test = model_test_new(mt);
	    
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LM);
		model_test_set_dfn(test, df);
		model_test_set_value(test, LM);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }	  
	}

	if (BP) {
	    record_test_result(LM, pval, "Breusch-Pagan");
	} else {
	    record_test_result(LM, pval, _("White's"));
	}
    }

    clear_model(&white);

    dataset_drop_last_variables(pdinfo->v - v, pZ, pdinfo);

    free(list);

    return err;
}

static int ar_list_max (const int *list) 
{
    int i, lmax = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] > lmax) {
	    lmax = list[i];
	}
    }

    return lmax;
}

/**
 * ar_func:
 * @list: list of lags plus dependent variable and list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list using the generalized 
 * Cochrane-Orcutt procedure for autoregressive errors.
 * 
 * Returns: #MODEL struct containing the results.
 */

MODEL ar_func (const int *list, double ***pZ, 
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    double diff, ess, tss, xx;
    int i, j, t, t1, t2, vc, yno, ryno, iter;
    int err, lag, maxlag, v = pdinfo->v;
    int *arlist = NULL, *rholist = NULL;
    int *reglist = NULL, *reglist2 = NULL;
    int pos, cpos;
    MODEL ar, rhomod;

    gretl_error_clear();

    gretl_model_init(&ar);
    gretl_model_init(&rhomod);

    pos = gretl_list_separator_position(list);

    arlist = malloc(pos * sizeof *arlist);
    reglist = malloc((list[0] - pos + 2) * sizeof *reglist);
    reglist2 = malloc((list[0] - pos + 2) * sizeof *reglist2);
    rholist = malloc((pos + 2) * sizeof *rholist);

    if (arlist == NULL || reglist == NULL || reglist2 == NULL ||
	rholist == NULL) {
	ar.errcode = E_ALLOC;
	goto bailout;
    }
    
    arlist[0] = pos - 1;
    for (i=1; i<pos; i++) {
	arlist[i] = list[i];
    }

    reglist2[0] = reglist[0] = list[0] - pos;
    for (i=1; i<=reglist[0]; i++) {
	reglist[i] = list[i + pos];
    }

    rholist[0] = arlist[0] + 1;
    maxlag = ar_list_max(arlist);

    cpos = reglist_check_for_const(reglist, (const double **) *pZ, pdinfo);

    /* special case: ar 1 ; ... => use AR1 */
    if (arlist[0] == 1 && arlist[1] == 1) {
	xx = estimate_rho(reglist, pZ, pdinfo, OPT_NONE, prn, &err);
	if (err) {
	    ar.errcode = err;
	} else {
	    ar = ar1_lsq(reglist, pZ, pdinfo, AR1, OPT_NONE, xx);
	}
	goto bailout;
    }

    /* first pass: estimate model via OLS: use OPT_M to generate an
       error in case of missing values within sample range 
    */
    ar = lsq(reglist, pZ, pdinfo, OLS, OPT_A | OPT_M);
    if (ar.errcode) {
	goto bailout;
    }

    /* allocate space for the uhat terms and transformed data */
    if (dataset_add_series(arlist[0] + 1 + reglist[0], pZ, pdinfo)) {
	ar.errcode = E_ALLOC;
	goto bailout;
    }

    yno = reglist[1];
    t1 = ar.t1; t2 = ar.t2;
    rholist[1] = v;

    pprintf(prn, "%s\n\n", _("Generalized Cochrane-Orcutt estimation"));
    bufspace(17, prn);
    /* xgettext:no-c-format */
    pputs(prn, _("ITER             ESS           % CHANGE"));
    pputs(prn, "\n\n");

    /* now loop while ess is changing */
    diff = 1.0e6;
    ess = 0.0;
    for (iter = 1; iter <= 20 && diff > 0.005; iter++) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t < t1 || t > t2) {
		(*pZ)[v][t] = NADBL;
	    } else {
		/* special computation of uhat */
		xx = (*pZ)[yno][t];
		for (j=0; j<reglist[0]-1; j++) {
		    xx -= ar.coeff[j] * (*pZ)[reglist[j+2]][t];
		}
		(*pZ)[v][t] = xx;
	    }
	}		
	for (i=1; i<=arlist[0]; i++) {
	    lag = arlist[i];
	    rholist[1+i] = v + i;
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 + lag || t > t2) {
		    (*pZ)[v+i][t] = NADBL;
		} else {
		    (*pZ)[v+i][t] = (*pZ)[v][t-lag];
		}
	    }
	}

	/* now estimate the rho terms */
	if (iter > 1) {
	    clear_model(&rhomod);
	}
	rhomod = lsq(rholist, pZ, pdinfo, OLS, OPT_A);

	/* and rho-transform the data */
	ryno = vc = v + i;
	for (i=1; i<=reglist[0]; i++) {
	    for (t=0; t<pdinfo->n; t++) {
		if (t < t1 + maxlag || t > t2) {
		    (*pZ)[vc][t] = NADBL;
		} else {
		    xx = (*pZ)[reglist[i]][t];
		    for (j=1; j<=arlist[0]; j++) {
			lag = arlist[j];
			xx -= rhomod.coeff[j-1] * (*pZ)[reglist[i]][t-lag];
		    }
		    (*pZ)[vc][t] = xx;
		}
	    }
	    reglist2[i] = vc++;
	}

	/* estimate the transformed model */
	clear_model(&ar);
	ar = lsq(reglist2, pZ, pdinfo, OLS, OPT_A);

        if (iter > 1) {
	    diff = 100 * (ar.ess - ess) / ess;
	}

        if (diff < 0.0) {
	    diff = -diff;
	}

	ess = ar.ess;
	pprintf(prn, "%16c%3d %20f ", ' ', iter, ess);

	if (iter > 1) {
	    pprintf(prn, "%13.3f\n", diff);
	} else {
	    pputc(prn, '\n');
	}
    } /* end "ess changing" loop */

    for (i=0; i<=reglist[0]; i++) {
	ar.list[i] = reglist[i];
    }
    if (cpos > 0) {
	ar.ifc = 1;
    }
    if (ar.ifc) {
	if (!gretl_model_get_int(&ar, "effconst")) {
	    ar.dfn -= 1;
	}
    }
    ar.ci = AR;

    /* special computation of fitted values */
    for (t=t1; t<=t2; t++) {
	xx = 0.0;
	for (j=2; j<=reglist[0]; j++) { 
	    xx += ar.coeff[j-2] * (*pZ)[reglist[j]][t];
	}
	ar.uhat[t] = (*pZ)[yno][t] - xx;
	for (j=1; j<=arlist[0]; j++) {
	    if (t - t1 >= arlist[j]) {
		xx += rhomod.coeff[j-1] * ar.uhat[t - arlist[j]];
	    }
	}
	ar.yhat[t] = xx;
    }

    for (t=t1; t<=t2; t++) { 
	ar.uhat[t] = (*pZ)[yno][t] - ar.yhat[t];
    }

    ar.rsq = gretl_corr_rsq(ar.t1, ar.t2, (*pZ)[reglist[1]], ar.yhat);
    ar.adjrsq = 1.0 - ((1.0 - ar.rsq) * (ar.nobs - 1.0) / ar.dfd);

    /* special computation of TSS */
    xx = gretl_mean(ar.t1, ar.t2, (*pZ)[ryno]);
    tss = 0.0;
    for (t=ar.t1; t<=ar.t2; t++) {
	tss += ((*pZ)[ryno][t] - xx) * ((*pZ)[ryno][t] - xx);
    }
    ar.fstt = ar.dfd * (tss - ar.ess) / (ar.dfn * ar.ess);
    ls_criteria(&ar);
    ar.dw = dwstat(maxlag, &ar, (const double **) *pZ);
    ar.rho = rhohat(maxlag, ar.t1, ar.t2, ar.uhat);

    dataset_drop_last_variables(arlist[0] + 1 + reglist[0], pZ, pdinfo);

    if (gretl_model_add_arinfo(&ar, maxlag)) {
	ar.errcode = E_ALLOC;
    } else {
	for (i=0; i<=arlist[0]; i++) { 
	    ar.arinfo->arlist[i] = arlist[i];
	    if (i >= 1) {
		ar.arinfo->rho[i-1] = rhomod.coeff[i-1];
		ar.arinfo->sderr[i-1] = rhomod.sderr[i-1];
	    }
	}
    }
    clear_model(&rhomod);

    if (!ar.errcode) {
	set_model_id(&ar);
    }  

 bailout:

    free(reglist);
    free(reglist2);
    free(rholist);
    free(arlist);

    return ar;
}

static int modelvar_iszero (const MODEL *pmod, const double *x)
{
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!model_missing(pmod, t) && floatneq(x[t], 0.0)) {
	    return 0;
	}
    }

    return 1;
}

/* From position 2 to end of list, omits variables with all zero
   observations and re-packs the rest of them */

static void omitzero (MODEL *pmod, const double **Z, const DATAINFO *pdinfo,
		      gretlopt opt)
{
    int *zlist = NULL;
    int i, v, offset;
    double x = 0.0;

    offset = (pmod->ci == WLS)? 3 : 2;

    if (!(opt & OPT_A)) {
	zlist = gretl_null_list();
    }

    for (i=pmod->list[0]; i>=offset; i--) {
        v = pmod->list[i];
        if (modelvar_iszero(pmod, Z[v])) {
	    if (zlist != NULL) {
		gretl_list_append_term(&zlist, v);
	    }
	    fprintf(stderr, "Deleting var %d at list pos %d: all zero\n", v, i);
	    gretl_list_delete_at_pos(pmod->list, i);
	}
    }

    if (pmod->nwt) {
	int t, wtzero;

	for (i=pmod->list[0]; i>=offset; i--) {
	    v = pmod->list[i];
	    wtzero = 1;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		x = Z[v][t] * Z[pmod->nwt][t];
		if (!model_missing(pmod, t) && floatneq(x, 0.0)) {
		    wtzero = 0;
		    break;
		}
	    }
	    if (wtzero) {
		if (zlist != NULL) {
		    gretl_list_append_term(&zlist, v);
		}
		gretl_list_delete_at_pos(pmod->list, i);
	    }
	}
    }

    if (zlist != NULL) {
	if (zlist[0] > 0) {
	    gretl_model_set_list_as_data(pmod, "zerolist", zlist);
	} else {
	    free(zlist);
	}
    }
}

/* lagdepvar: attempt to detect presence of a lagged dependent
   variable among the regressors -- if found, return the position of
   this lagged var in the list; otherwise return 0
*/

static int 
lagdepvar (const int *list, const double **Z, const DATAINFO *pdinfo) 
{
    char depvar[VNAMELEN], othervar[VNAMELEN];
    char *p;
    int i, t, ret = 0;

    strcpy(depvar, pdinfo->varname[list[1]]);

    for (i=2; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    break;
	}
	strcpy(othervar, pdinfo->varname[list[i]]);
	p = strrchr(othervar, '_');
	if (p != NULL && isdigit(*(p + 1))) {
	    /* looks like a lag */
	    size_t len = strlen(othervar) - strlen(p);

	    if (!strncmp(depvar, othervar, len)) {
		int gotlag = 1;

		/* strong candidate for lagged depvar, but make sure */
		for (t=pdinfo->t1+1; t<=pdinfo->t2; t++) {
		    if (Z[list[1]][t-1] != Z[list[i]][t]) {
			gotlag = 0;
			break;
		    }
		}
		if (gotlag) {
		    ret = i;
		    break;
		}
	    }
	}
    } 

    return ret;
}

static void arch_test_print_simple (int order, double LM,
				    double pval, PRN *prn)
{
    pprintf(prn, "\n%s %d\n", _("Test for ARCH of order"), order);
    pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic"), LM);
    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n", 
	    _("with p-value"), _("Chi-square"), 
	    order, LM, pval); 
}

static void 
arch_test_save_or_print (const gretl_matrix *b, const gretl_matrix *V,
			 int T, int order, double rsq, MODEL *pmod, 
			 gretlopt opt, PRN *prn)
{
    ModelTest *test = model_test_new(GRETL_TEST_ARCH);
    double LM = T * rsq;
    double pv = chisq_cdf_comp(order, LM);

    if (V != NULL) {
	int i, k = order + 1;
	double *se = malloc(k * sizeof *se);
	char **names;

	names = strings_array_new_with_length(k, 16);

	if (se != NULL && names != NULL) {
	    pputc(prn, '\n');
	    pprintf(prn, _("Test for ARCH of order %d"), order);
	    pputs(prn, "\n\n");

	    for (i=0; i<k; i++) {
		se[i] = sqrt(gretl_matrix_get(V, i, i));
		sprintf(names[i], "alpha(%d)", i);
	    }

	    print_coeffs(b->val, se, (const char **) names, 
			 k, T - k, ARCH, prn);
	}

	free(se);
	free_strings_array(names, k);
    }

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_LM);
	model_test_set_order(test, order);
	model_test_set_dfn(test, order);
	model_test_set_value(test, LM);
	model_test_set_pvalue(test, pv);

	if (opt & OPT_Q) {
	    arch_test_print_simple(order, LM, pv, prn);
	} else {
	    int heading = (V == NULL);

	    gretl_model_test_print_direct(test, heading, prn);
	}

	if (pmod != NULL && (opt & OPT_S)) {
	    maybe_add_test_to_model(pmod, test);
	} else {
	    model_test_free(test);
	}
    }	    

    record_test_result(LM, pv, "ARCH");
}

static int real_arch_test (const double *u, int T, int order, 
			   MODEL *pmod, const DATAINFO *pdinfo, 
			   gretlopt opt, PRN *prn)
{
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *V = NULL;
    int i, k, s, t;
    double x, s2, rsq;
    double *ps2 = NULL;
    int err = 0;

    gretl_error_clear();

    if (order < 1 || order > T - 1) {
	sprintf(gretl_errmsg, _("Invalid lag order for arch (%d)"), order);
	return E_DATA;
    }

    T -= order;
    k = order + 1;

    X = gretl_matrix_alloc(T, k);
    y = gretl_column_vector_alloc(T);
    b = gretl_column_vector_alloc(k);

    if (X == NULL || y == NULL || b == NULL) {
	gretl_matrix_free(X);
	gretl_matrix_free(y);
	gretl_matrix_free(b);
	return E_ALLOC;
    }

    if (!(opt & OPT_Q)) {
	V = gretl_matrix_alloc(k, k);
	if (V == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	ps2 = &s2;
    }

    /* fill out the matrices with squared residuals
       and lags of same */

    for (i=0; i<k; i++) {
	for (t=0; t<T; t++) {
	    s = t + order;
	    if (i == 0) {
		x = u[s];
		gretl_vector_set(y, t, x * x);
		gretl_matrix_set(X, t, i, 1.0);
	    } else {
		x = u[s - i];
		gretl_matrix_set(X, t, i, x * x);
	    }
	}
    }

    err = gretl_matrix_ols(y, X, b, V, NULL, ps2);

    if (!err) {
	rsq = gretl_matrix_r_squared(y, X, b, &err);
    }

    if (!err) {
	arch_test_save_or_print(b, V, T, order, rsq, pmod, 
				opt, prn);
    }

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

/**
 * arch_test:
 * @pmod: model to be tested.
 * @order: lag order for ARCH process.
 * @pdinfo: information on the data set.
 * @opt: if flags include %OPT_S, save test results to model;
 * if %OPT_Q, be less verbose.
 * @prn: gretl printing struct.
 *
 * Tests @pmod for AutoRegressive Conditional Heteroskedasticity.  
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int arch_test (MODEL *pmod, int order, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn)
{
    int err;

    if (pmod->missmask != NULL) {
	err = E_MISSDATA;
    } else {
	const double *u = pmod->uhat + pmod->t1;

	if (order == 0) {
	    /* use data frequency as default lag order */
	    order = pdinfo->pd;
	}

	err = real_arch_test(u, pmod->nobs, order, pmod, pdinfo, 
			     opt, prn);
    }

    return err;
}

int array_arch_test (const double *u, int n, int order, 
		     gretlopt opt, PRN *prn)
{
    return real_arch_test(u, n, order, NULL, NULL, opt, prn);
}

/**
 * arch_model:
 * @list: dependent variable plus list of regressors.
 * @order: lag order for ARCH process.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list via weighted least squares,
 * with the weights based on the predicted error variances from
 * an auxiliary regression of the squared residuals on their lagged
 * values.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arch_model (const int *list, int order, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    MODEL amod;
    int *wlist = NULL, *alist = NULL;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int oldv = pdinfo->v;
    int i, t, nwt, k, n = pdinfo->n;
    double *a = NULL;
    double *se = NULL;
    double xx;

    gretl_error_clear();
    gretl_model_init(&amod);

    if (order == 0) {
	/* use data frequency as default lag order */
	order = pdinfo->pd;
    }

    if (order < 1 || order > T - list[0]) {
	amod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, _("Invalid lag order for arch (%d)"), order);
	return amod;
    }

    if (dataset_add_series(order + 1, pZ, pdinfo)) {
	amod.errcode = E_ALLOC;
    } else {
	alist = gretl_list_new(order + 2);
	if (alist == NULL) {
	    amod.errcode = E_ALLOC;
	}
    }

    if (amod.errcode) {
	goto bailout;
    }

    /* start list for aux regression */
    alist[1] = pdinfo->v - order - 1;
    alist[2] = 0;

    /* run initial OLS and get squared residuals */
    amod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_M);
    if (amod.errcode) {
	goto bailout;
    }

    k = pdinfo->v - order - 1;
    strcpy(pdinfo->varname[k], "utsq");
    for (t=0; t<n; t++) {
	(*pZ)[k][t] = NADBL;
    }
    for (t=amod.t1; t<=amod.t2; t++) {
	xx = amod.uhat[t];
	(*pZ)[k][t] = xx * xx;
    }
    /* also lags of squared resids */
    for (i=1; i<=order; i++) {
	k =  pdinfo->v - order + i - 1;
	alist[i+2] = k;
	sprintf(pdinfo->varname[k], "utsq_%d", i);
	for (t=0; t<n; t++) {
	    (*pZ)[k][t] = NADBL;
	}
	for (t=amod.t1+i; t<=amod.t2; t++) {
	    (*pZ)[k][t] = (*pZ)[alist[1]][t-i];
	}
    }

    /* run auxiliary regression */
    clear_model(&amod);
    amod = lsq(alist, pZ, pdinfo, OLS, OPT_A);
    if (amod.errcode) {
	goto bailout;
    }

    /* steal the coefficients for reference */
    a = amod.coeff;
    amod.coeff = NULL;
    se = amod.sderr;
    amod.sderr = NULL;

    /* do weighted estimation */
    wlist = gretl_list_new(list[0] + 1);

    if (wlist == NULL) {
	amod.errcode = E_ALLOC;
    } else {
	/* construct the weight variable */
	nwt = wlist[1] = pdinfo->v - 1;
	strcpy(pdinfo->varname[nwt], "1/sigma");

	for (i=2; i<=wlist[0]; i++) {
	    wlist[i] = list[i-1];
	}

	k = pdinfo->v - order - 1;

	for (t=amod.t1; t<=amod.t2; t++) {
	    xx = amod.yhat[t];
	    if (xx <= 0.0) {
		xx = (*pZ)[k][t];
	    }
	    (*pZ)[nwt][t] = 1.0 / xx; /* FIXME is this right? */
	}

	clear_model(&amod);
	amod = lsq(wlist, pZ, pdinfo, WLS, OPT_NONE);
	amod.ci = ARCH;

	if (!amod.errcode) {
	    gretl_model_set_int(&amod, "arch_order", order);
	    gretl_model_set_data(&amod, "arch_coeff", a,
				 GRETL_TYPE_DOUBLE_ARRAY,
				 (order + 1) * sizeof *a);
	    gretl_model_set_data(&amod, "arch_sderr", se,
				 GRETL_TYPE_DOUBLE_ARRAY,
				 (order + 1) * sizeof *se);
	}
    }

 bailout:

    if (alist != NULL) free(alist);
    if (wlist != NULL) free(wlist);

    dataset_drop_last_variables(pdinfo->v - oldv, pZ, pdinfo); 

    return amod;
}

/**
 * lad:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 *
 * Estimate the model given in @list using the method of Least
 * Absolute Deviation (LAD).
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL lad (const int *list, double ***pZ, DATAINFO *pdinfo)
{
    MODEL lmod;
    void *handle;
    int (*lad_driver) (MODEL *, double **, DATAINFO *);

    /* run an initial OLS to "set the model up" and check for errors.
       the lad_driver function will overwrite the coefficients etc.
    */

    lmod = lsq(list, pZ, pdinfo, OLS, OPT_A);

    if (lmod.errcode) {
        return lmod;
    }

    lad_driver = get_plugin_function("lad_driver", &handle);

    if (lad_driver == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	lmod.errcode = E_FOPEN;
	return lmod;
    }

    (*lad_driver) (&lmod, *pZ, pdinfo);
    close_plugin(handle);

    if (lmod.errcode == 0) {
	set_model_id(&lmod);
    }

    return lmod;
}

/**
 * quantreg:
 * @parm: 
 * @list: model specification: dependent var and regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt:
 * @prn:
 *
 * Estimate the model given in @list using the method of 
 * quantile regression.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL quantreg (const char *parm, const int *list, 
		double ***pZ, DATAINFO *pdinfo,
		gretlopt opt, PRN *prn)
{
    MODEL qmod;
    void *handle;
    int (*rq_driver) (const char *, MODEL *, double ***, DATAINFO *,
		      gretlopt, PRN *);
    gretlopt olsopt = (OPT_A | OPT_M);

    /* Run an initial OLS to "set the model up" and check for errors.
       the driver function will overwrite the coefficients, etc.
       For now we make life easier by rejecting within-sample missing
       values (OPT_M).
    */

    if (opt & OPT_R) {
	olsopt |= OPT_R;
    }

    qmod = lsq(list, pZ, pdinfo, OLS, olsopt);

    if (qmod.errcode) {
        return qmod;
    }

    rq_driver = get_plugin_function("rq_driver", &handle);

    if (rq_driver == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	qmod.errcode = E_FOPEN;
	return qmod;
    }

    (*rq_driver) (parm, &qmod, pZ, pdinfo, opt, prn);
    close_plugin(handle);

    if (qmod.errcode == 0) {
	set_model_id(&qmod);
    }

    return qmod;
}

/**
 * arma:
 * @list: dependent variable, AR and MA orders, and any exogenous
 * regressors.
 * @pqspec: string giving specific non-seasonal AR/MA lags (or %NULL).
 * @Z: data array.
 * @pdinfo: information on the data set.
 * @opt: options: may include %OPT_S to suppress intercept, %OPT_V
 * for verbose results, %OPT_X to use X-12-ARIMA, %OPT_C to put
 * X-12-ARIMA into conditional maximum-likelihood mode.
 * @PRN: for printing details of iterations (or %NULL). 
 *
 * Calculate ARMA estimates, using either native gretl code or
 * by invoking X-12-ARIMA.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arma (const int *list, const char *pqspec, 
	    const double **Z, const DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn)
{
    MODEL armod;
    void *handle;
    MODEL (*arma_func) (const int *, const char *,
			const double **, const DATAINFO *, 
			gretlopt, PRN *);

    gretl_error_clear();

    if (opt & OPT_X && (pdinfo->t2 - pdinfo->t1) > 719) {
	strcpy(gretl_errmsg, _("X-12-ARIMA can't handle more than 720 observations.\n"
			       "Please select a smaller sample."));
	armod.errcode = E_DATA;
	return armod;
    }	

    if (opt & OPT_X) {
	arma_func = get_plugin_function("arma_x12_model", &handle);
    } else {
	arma_func = get_plugin_function("arma_model", &handle);
    }

    if (arma_func == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	gretl_model_init(&armod);
	armod.errcode = E_FOPEN;
	return armod;
    }

    armod = (*arma_func) (list, pqspec, Z, pdinfo, opt, prn);

    close_plugin(handle);
    set_model_id(&armod);

    return armod;
} 

/**
 * tobit_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @prn: printing struct for iteration info (or %NULL is this is not
 * wanted).
 *
 * Produce Tobit estimates of the model given in @list.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tobit_model (const int *list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL tmod;
    void *handle;
    MODEL (* tobit_estimate) (const int *, double ***, DATAINFO *, PRN *);

    gretl_error_clear();

    tobit_estimate = get_plugin_function("tobit_estimate", &handle);
    if (tobit_estimate == NULL) {
	gretl_model_init(&tmod);
	tmod.errcode = E_FOPEN;
	return tmod;
    }

    tmod = (*tobit_estimate) (list, pZ, pdinfo, prn);

    close_plugin(handle);
    set_model_id(&tmod);

    return tmod;
}

static int get_offset_var (int *list)
{
    int l0 = list[0];
    int ret = 0;

    if (list[l0 - 1] == LISTSEP) {
	ret = list[l0];
	list[0] -= 2;
    }

    return ret;
}

/**
 * poisson_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: printing struct for iteration info (or %NULL is this is not
 * wanted).
 *
 * Estimate the Poisson regression model given in @list using ML.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL poisson_model (const int *list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    MODEL pmodel;
    void *handle;
    int *listcpy;
    int offvar;
    int (* poisson_estimate) (MODEL *, int, double ***, DATAINFO *, PRN *);

    gretl_error_clear();

    gretl_model_init(&pmodel);

    listcpy = gretl_list_copy(list);
    if (listcpy == NULL) {
	pmodel.errcode = E_ALLOC;
        return pmodel;
    }

    offvar = get_offset_var(listcpy);

    /* run an initial OLS to "set the model up" and check for errors.
       the poisson_estimate_driver function will overwrite the
       coefficients etc.
    */

    pmodel = lsq(listcpy, pZ, pdinfo, OLS, OPT_A);
    free(listcpy);

    if (pmodel.errcode) {
        return pmodel;
    }

    poisson_estimate = get_plugin_function("poisson_estimate", &handle);

    if (poisson_estimate == NULL) {
	pmodel.errcode = E_FOPEN;
	return pmodel;
    }

    (*poisson_estimate) (&pmodel, offvar, pZ, pdinfo, prn);

    close_plugin(handle);

    set_model_id(&pmodel);

    return pmodel;
}

/**
 * heckit_model:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt: option flags (may include %OPT_V for verbose output).
 * @prn: printing struct for iteration info (or %NULL is this is not
 * wanted).
 *
 * Produce Heckit estimates of the model given in @list. The list must
 * include a separator to divide the main equation from the selection
 * equation.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL heckit_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn)
{
    MODEL hmod;
    void *handle;
    MODEL (* heckit_estimate) (const int *, double ***, DATAINFO *, 
			       gretlopt, PRN *);

    gretl_error_clear();

    heckit_estimate = get_plugin_function("heckit_estimate", &handle);
    if (heckit_estimate == NULL) {
	gretl_model_init(&hmod);
	hmod.errcode = E_FOPEN;
	return hmod;
    }

    hmod = (*heckit_estimate) (list, pZ, pdinfo, opt, prn);

    close_plugin(handle);

    set_model_id(&hmod);

    return hmod;
}

/**
 * garch:
 * @list: dependent variable plus arch and garch orders.
 * @pZ: pointer to data array.
 * @pdinfo: information on the data set.
 * @opt: can specify robust standard errors and VCV.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Calculate GARCH estimates.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL garch (const int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn)
{
    MODEL gmod;
    void *handle;
    PRN *myprn;
    MODEL (*garch_model) (const int *, double ***, DATAINFO *, PRN *,
			  gretlopt);

    gretl_error_clear();

    garch_model = get_plugin_function("garch_model", &handle);

    if (garch_model == NULL) {
	gretl_model_init(&gmod);
	gmod.errcode = E_FOPEN;
	return gmod;
    }

    if (opt & OPT_V) {
	myprn = prn;
    } else {
	myprn = NULL;
    }

    gmod = (*garch_model) (list, pZ, pdinfo, myprn, opt);

    close_plugin(handle);

    set_model_id(&gmod);

    return gmod;
} 

/**
 * mp_ols:
 * @list: specification of variables to use.
 * @Z: data array.
 * @pdinfo: information on the data set.
 *
 * Estimate an OLS model using multiple-precision arithmetic
 * via the GMP library.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL mp_ols (const int *list, const double **Z, DATAINFO *pdinfo)
{
    void *handle = NULL;
    int (*mplsq)(const int *, const int *, const int *, 
		 const double **, DATAINFO *, char *, 
		 MODEL *, gretlopt);
    MODEL mpmod;

    gretl_model_init(&mpmod);

    mplsq = get_plugin_function("mplsq", &handle);
    if (mplsq == NULL) {
	mpmod.errcode = 1;
	return mpmod;
    }

    if (gretl_list_has_separator(list)) {
	int *base = NULL;
	int *poly = NULL;

	mpmod.errcode = gretl_list_split_on_separator(list, &base, &poly);
	if (mpmod.errcode == 0 && (base == NULL || poly == NULL)) {
	    mpmod.errcode = E_ARGS;
	} else {
	    mpmod.errcode = (*mplsq)(base, poly, NULL, Z, pdinfo,  
				     gretl_errmsg, &mpmod, OPT_S);
	}
	free(base);
	free(poly);
    } else {
	mpmod.errcode = (*mplsq)(list, NULL, NULL, Z, pdinfo,  
				 gretl_errmsg, &mpmod, OPT_S); 
    }

    close_plugin(handle);

    set_model_id(&mpmod);

    return mpmod;
}

static int check_panel_options (gretlopt opt)
{
    if ((opt & OPT_U) && (opt & OPT_W)) {
	/* can't specify random effects + weighted least squares */
	return E_BADOPT;
    } else if ((opt & OPT_T) && !(opt & OPT_W)) {
	/* iterate option requires weighted least squares option */
	return E_BADOPT;
    } else if (incompatible_options(opt, OPT_B | OPT_U | OPT_P)) {
	/* mutually exclusive estimator requests */
	return E_BADOPT;
    }

    return 0;
}

/**
 * panel_model:
 * @list: regression list (dependent variable plus independent 
 * variables).
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the (panel) data set.
 * @opt: can include %OPT_Q (quiet estimation), %OPT_S
 * (silent estimation), %OPT_R (random effects model),
 * %OPT_W (weights based on the error variance for the
 * respective cross-sectional units), %OPT_T (iterate, only
 * available in conjunction with %OPT_W).
 * @prn: printing struct (or %NULL).
 *
 * Calculate estimates for a panel dataset, using fixed
 * effects (the default), random effects, or weighted
 * least squares based on the respective variances for the
 * cross-sectional units.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn)
{
    MODEL mod;

    gretl_error_clear();

    if (check_panel_options(opt)) {
	gretl_model_init(&mod);
	mod.errcode = E_BADOPT;
    } else if (opt & OPT_W) {
	mod = panel_wls_by_unit(list, pZ, pdinfo, opt, prn);
    } else {
	mod = real_panel_model(list, pZ, pdinfo, opt, prn);
    }

    return mod;
}

/**
 * ivreg:
 * @list: regression list, as per ivreg documentation.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the (panel) data set.
 * @opt: can include %OPT_Q (quiet estimation), %OPT_L
 * (LIML), ... (FIXME)
 *
 * Calculate IV estimates using either 2sls or LIML.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL ivreg (const int *list, double ***pZ, DATAINFO *pdinfo,
	     gretlopt opt)
{
    MODEL mod;

    gretl_error_clear();

    if ((opt & (OPT_T | OPT_I)) && !(opt & OPT_G)) {
	/* two-step and iterate options are GMM-only */
	gretl_model_init(&mod);
	mod.errcode = E_BADOPT;
	return mod;
    }

    if (opt & OPT_L) {
	mod = single_equation_liml(list, pZ, pdinfo, opt);
    } else if (opt & OPT_G) {
	mod = ivreg_via_gmm(list, pZ, pdinfo, opt);
    } else {
	mod = tsls(list, pZ, pdinfo, opt);
    }

    return mod;
}

/**
 * arbond_model:
 * @list: regression list.
 * @istr: may contain additional instrument specification.
 * @Z: data array.
 * @pdinfo: information on the (panel) data set.
 * @opt: to be hooked up.
 * @prn: printing struct (or %NULL).
 *
 * To be written.  This function is currently just for
 * testing.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL arbond_model (const int *list, const char *istr, const double **Z, 
		    const DATAINFO *pdinfo, gretlopt opt, 
		    PRN *prn)
{
    void *handle = NULL;
    MODEL (*arbond_estimate) (const int *, const char *, const double **, 
			      const DATAINFO *, gretlopt, PRN *);
    MODEL mod;

    gretl_model_init(&mod);

    arbond_estimate = get_plugin_function("arbond_estimate", &handle);
    if (arbond_estimate == NULL) {
	mod.errcode = 1;
	return mod;
    }

    mod = (*arbond_estimate)(list, istr, Z, pdinfo, opt, prn);

    close_plugin(handle);

    if (!mod.errcode) {
	set_model_id(&mod);
    }

    return mod;    
}

/**
 * groupwise_hetero_test:
 * @pmod: pooled OLS model to be tested.
 * @pZ: pointer to data array.
 * @pdinfo: information on the (panel) data set.
 * @prn: for printing details of iterations (or %NULL).
 *
 * Calculates iterated WLS estimates using weights based on the error
 * variance for the cross-sectional units and performs a Wald test
 * for the null hypothesis that the error variance is uniform
 * across the units.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int groupwise_hetero_test (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			   PRN *prn)
{
    MODEL wmod;
    int err;

    if (!POOLED_MODEL(pmod)) {
	return E_NOTIMP;
    }

    if (!dataset_is_panel(pdinfo)) {
	strcpy(gretl_errmsg, _("This test is only available for panel data"));
	return 1;
    }

    wmod = panel_wls_by_unit(pmod->list, pZ, pdinfo, OPT_T | OPT_A, prn);
    err = wmod.errcode;

    if (!err) {
	gretl_model_set_auxiliary(&wmod, AUX_GROUPWISE);
	printmodel(&wmod, pdinfo, OPT_NONE, prn);
    }

    clear_model(&wmod);

    return err;
}

