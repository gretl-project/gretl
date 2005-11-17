/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 by Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* GARCH plugin for gretl using the Fiorentini, Calzolari and 
   Panattoni mixed-gradient algorithm.
*/

#include "libgretl.h"
#include "libset.h"
#include "var.h"

#include "fcp.h"

#undef VPARM_DEBUG

#define VPARM_MAX 6            /* max number of variance parameters */
#define GARCH_PARAM_MAX 0.999

double vparm_init[VPARM_MAX];

static void add_garch_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				const int *list)
{
    int p = list[1];
    int q = list[2];
    int r = list[0] - 4;
    int i, j, np = 3 + p + q + r;

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    pmod->params = malloc(np * sizeof pmod->params);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    pmod->nparams = np;

    for (i=0; i<np; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    for (j=0; j<i; j++) free(pmod->params[j]);
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    pmod->errcode = E_ALLOC;
	    return;
	}
    }

    strcpy(pmod->params[0], pdinfo->varname[pmod->list[4]]);
    strcpy(pmod->params[1], pdinfo->varname[0]);

    j = 2;
    for (i=0; i<r; i++) {
	if (pmod->list[5+i] > 0) {
	    strcpy(pmod->params[j++], pdinfo->varname[pmod->list[5+i]]);
	}
    }

    strcpy(pmod->params[j++], "alpha(0)");

    for (i=0; i<q; i++) {
	sprintf(pmod->params[j++], "alpha(%d)", i + 1);
    }
    for (i=0; i<p; i++) {
	sprintf(pmod->params[j++], "beta(%d)", i + 1);
    }
}

static int make_packed_vcv (MODEL *pmod, double *vcv, int np,
			    int nc, double scale)
{
    const int nterms = np * (np + 1) / 2;
    double sfi, sfj;
    int i, j, k;

    free(pmod->vcv);
    pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return 1;  
    }

    for (i=0; i<np; i++) {
	if (i < nc) {
	    sfi = scale;
	} else if (i == nc) {
	    sfi = scale * scale;
	} else {
	    sfi = 1.0;
	}
	for (j=0; j<=i; j++) {
	    if (j < nc) {
		sfj = scale;
	    } else if (j == nc) {
		sfj = scale * scale;
	    } else {
		sfj = 1.0;
	    }
	    k = ijton(i, j, np);
	    pmod->vcv[k] = vcv[i + np * j] * sfi * sfj;
	}
    }

    return 0;
}

static int write_garch_stats (MODEL *pmod, const double **Z,
			      double scale, const DATAINFO *pdinfo,
			      const int *list, const double *theta, 
			      int nparam, int pad, const double *res,
			      const double *h)
{
    int err = 0;
    double *coeff, *sderr, *garch_h;
    double den;
    int ynum = list[4];
    int nvp = list[1] + list[2];
    int xvars = list[0] - 4;
    int i;

    coeff = realloc(pmod->coeff, nparam * sizeof *pmod->coeff);
    sderr = realloc(pmod->sderr, nparam * sizeof *pmod->sderr);

    if (coeff == NULL || sderr == NULL) return 1;

    for (i=0; i<nparam; i++) {
	coeff[i] = theta[i+1];
	sderr[i] = theta[i+nparam+1];
    }
    
    pmod->coeff = coeff;
    pmod->sderr = sderr;
   
    pmod->ncoeff = nparam;

    pmod->ess = 0.0;
    for (i=pmod->t1; i<=pmod->t2; i++) {
	pmod->uhat[i] = res[i + pad] * scale;
	pmod->ess += pmod->uhat[i] * pmod->uhat[i];
	pmod->yhat[i] =  Z[ynum][i] * scale - pmod->uhat[i];
    }

    /* set sigma to its unconditional or steady-state value */
    den = 1.0;
    for (i=0; i<nvp; i++) {
	den -= coeff[i+xvars+1];
    }
    pmod->sigma = sqrt(coeff[xvars] / den);

    pmod->adjrsq = NADBL; 
    pmod->fstt = NADBL;

    mle_aic_bic(pmod, 1);

    pmod->ci = GARCH;
    pmod->ifc = 1;
    
    add_garch_varnames(pmod, pdinfo, list);

    /* add predicted error variance to model */
    garch_h = malloc(pdinfo->n * sizeof *garch_h);
    if (garch_h != NULL) {
	for (i=0; i<pdinfo->n; i++) {
	    if (i < pmod->t1 || i > pmod->t2) {
		garch_h[i] = NADBL;
	    } else {
		garch_h[i] = h[i + pad] * scale * scale;
	    }
	}
	gretl_model_set_data(pmod, "garch_h", garch_h, 
			     pdinfo->n * sizeof *garch_h);
    }

    return err;
}

static int make_garch_dataset (const int *list, double **Z,
			       int bign, int pad, int nx,
			       double **py, double ***pX)
{
    double *y = NULL, **X = NULL;
    int i, k = 0, t;
    int xnum, ynum = list[4];

    /* If pad > 0 we have to create a newly allocated, padded
       dataset.  Otherwise we can use a virtual dataset, made
       up of pointers into the original dataset, Z. 
    */

    if (pad > 0) {
	y = malloc(bign * sizeof *y);
	if (y == NULL) return 1;
    } 

    if (nx > 0) {
	X = malloc(nx * sizeof *X);
	if (X == NULL) goto bailout;

	if (pad > 0) {
	    for (i=0; i<nx; i++) {
		X[i] = malloc(bign * sizeof **X);
		if (X[i] == NULL) {
		    for (t=0; t<i; t++) {
			free(X[t]);
		    }
		    free(X);
		    goto bailout;
		}
	    } 
	}  
    }

    if (pad > 0) {
	/* build padded dataset */
	for (t=0; t<bign; t++) {
	    if (t < pad) {
		y[t] = 0.0;
		for (i=0; i<nx; i++) {
		    X[i][t] = 0.0;
		}
	    } else {
		y[t] = Z[ynum][t-pad];
		if (nx > 0) k = 5;
		for (i=0; i<nx; i++) {
		    xnum = list[k++]; 
		    if (xnum == 0) xnum = list[k++];
		    X[i][t] = Z[xnum][t-pad];
		}
	    }
	}
	*py = y;
    } else {
	/* build virtual dataset */
	*py = Z[ynum];
	if (nx > 0) k = 5;
	for (i=0; i<nx; i++) {
	    xnum = list[k++]; 
	    if (xnum == 0) xnum = list[k++];
	    X[i] = Z[xnum];
	}
    }

    *pX = X;

    return 0;

 bailout:

    free(y);
    return E_ALLOC;
}

static int get_vopt (int robust)
{
    int vopt = get_garch_vcv_version();
    int ropt = get_garch_robust_vcv_version();

    /* The defaults: QML if "robust" option is in force,
       otherwise negative Hessian */
    if (vopt == VCV_UNSET) {
	if (robust) {
	    if (ropt == VCV_UNSET) {
		vopt = VCV_QML;
	    } else {
		vopt = ropt;
	    }
	} else {
	    vopt = VCV_HESSIAN;
	}
    }

    return vopt;
}

int do_fcp (const int *list, double **Z, double scale,
	    const DATAINFO *pdinfo, MODEL *pmod,
	    PRN *prn, gretlopt opt)
{
    int t1 = pmod->t1, t2 = pmod->t2;
    int ncoeff = pmod->ncoeff;
    int p = list[1];
    int q = list[2];
    double *y = NULL;
    double **X = NULL;
    double *h = NULL;
    double *amax = NULL; 
    double *res = NULL, *res2 = NULL;
    double *coeff = NULL, *b = NULL;
    double *vcv = NULL;
    int err = 0, iters = 0;
    int nobs, maxlag, bign, pad = 0;
    int i, nx, nparam, vopt;

    vopt = get_vopt(opt & OPT_R);

    nx = ncoeff - 1;
    maxlag = (p > q)? p : q; 
    nparam = ncoeff + p + q + 1;

    nobs = t2 + 1; /* number of obs in full dataset */

    if (maxlag > t1) {
	/* need to pad data series at start */
	pad = maxlag - t1;
    } 

    /* length of series to pass to garch_estimate */
    bign = nobs + pad;
	
    res2 = malloc(bign * sizeof *res2);
    res = malloc(bign * sizeof *res);
    h = malloc(bign * sizeof *h);
    amax = malloc(bign * sizeof *amax);
    if (res2 == NULL || res == NULL || 
	amax == NULL || h == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    for (i=0; i<bign; i++) {
	res2[i] = res[i] = amax[i] = 0.0;
    }   
 
    coeff = malloc(ncoeff * sizeof *coeff);
    b = malloc(ncoeff * sizeof *b);
    if (coeff == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;	
    }

    vcv = malloc((nparam * nparam) * sizeof *vcv);
    if (vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    for (i=0; i<nparam * nparam; i++) {
	vcv[i] = 0.0;
    } 

    /* create dataset for garch estimation */
    err = make_garch_dataset(list, Z, bign, pad, nx, &y, &X);
    if (err) {
	goto bailout;
    }

    /* initial coefficients from OLS */
    for (i=0; i<ncoeff; i++) {
	coeff[i] = pmod->coeff[i];
	b[i] = 0.0;
    }

    /* for compatibility with FCP... */
    amax[1] = q;
    amax[2] = p; 

    /* initialize variance parameters */
    amax[0] = vparm_init[0];
    for (i=0; i<p+q; i++) {
	amax[i+3] = vparm_init[i+1];
    }

    err = garch_estimate(t1 + pad, t2 + pad, bign, 
			 (const double **) X, nx, coeff, ncoeff, 
			 vcv, res2, res, h, y, amax, b, scale, &iters,
			 prn, vopt);

    if (err != 0) {
	pmod->errcode = err;
    } else {
	int nparam = ncoeff + p + q + 1;

	for (i=1; i<=nparam; i++) {
	    if (i <= ncoeff) {
		amax[i] *= scale;
		amax[i + nparam] *= scale;
	    } else if (i == ncoeff + 1) {
		amax[i] *= scale * scale;
		amax[i + nparam] *= scale * scale;
	    }
	    pprintf(prn, "theta[%d]: %#14.6g (%#.6g)\n", i-1, amax[i], 
		    amax[i + nparam]);
	}
	pputc(prn, '\n');

	pmod->lnL = amax[0];
	write_garch_stats(pmod, (const double **) Z, scale, pdinfo, 
			  list, amax, nparam, pad, res, h);
	make_packed_vcv(pmod, vcv, nparam, ncoeff, scale);
	gretl_model_set_int(pmod, "iters", iters);
	gretl_model_set_int(pmod, "garch_vcv", vopt);
    }

 bailout:

    free(res2);
    free(res);
    free(h);
    free(amax);    
    free(coeff);
    free(b);
    free(vcv); 

    if (pad > 0) {
	/* don't free y if it's just a pointer into Z */
	free(y);
    }
    
    if (X != NULL) {
	if (pad > 0) {
	    /* don't free the X[i] if they're just pointers into
	       the original data matrix */
	    for (i=0; i<nx; i++) {
		free(X[i]);
	    }
	}
	free(X);
    }

    return err;
}

static int 
add_uhat_squared (const MODEL *pmod, double scale,
		  double ***pZ, DATAINFO *pdinfo)
{
    int t, v = pdinfo->v;

    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    for (t=0; t<pdinfo->n; t++) {
	double u = pmod->uhat[t];

	if (na(u)) {
	    (*pZ)[v][t] = NADBL;
	} else {
	    u /= scale;
	    (*pZ)[v][t] = u * u;
	}
    }

    strcpy(pdinfo->varname[v], "uhat2");

    return 0;
}

/*
  p and q are the GARCH orders
  ao = max(q,p) is the ar order
  mo = q is the ma order

  it is assumed that armapar contains the arma parameters 
  in the following order:
  armapar[0] : intercept
  armapar[1..ao] : ar terms
  armapar[ao+1..ao+mo] : ma terms
*/

static void
garchpar_from_armapar (const double *armapar, int q, int p)
{
    double x, sum_ab = 0.0;
    int ao = (p > q)? p : q;
    int mo = q;
    int i;

#ifdef VPARM_DEBUG
    for (i=0; i<1+ao+mo; i++) {
	fprintf(stderr, "armapar[%d] = %#12.6g\n", i, armapar[i]);
    }
#endif

    for (i=1; i<=p; i++) {
	x = 0.0;
	if (i <= ao) {
	    x += armapar[i];
	} 
	if (i<=mo) {
	    x += armapar[p+i];
	} 
	vparm_init[i] = (x < 0.0) ? 0.01 : x;
	sum_ab += vparm_init[i];
    }

    for (i=1; i<=q; i++) {
	x = armapar[p+i];
	vparm_init[p+i] = (x > 0.0) ? 0 : -x;
	sum_ab += vparm_init[p+i];
    }

#ifdef VPARM_DEBUG
    fprintf(stderr, "sum_ab = %#12.6g\n", sum_ab);
#endif

    if (sum_ab > GARCH_PARAM_MAX) {
	for (i=1; i<=p+q; i++) {
	    vparm_init[i] *= GARCH_PARAM_MAX / sum_ab;
	}
	sum_ab = GARCH_PARAM_MAX;
    }

    vparm_init[0] = armapar[0];
}

static int 
garch_init_by_arma (const MODEL *pmod, const int *garchlist, 
		    double scale, double ***pZ, DATAINFO *pdinfo)
{
    MODEL amod;
    int q = garchlist[1], p = garchlist[2];
    int v = pdinfo->v;
    int *list = NULL;
    int err = 0;

    /* for now we'll try this only for GARCH up to (2,2) */
    if (q > 2 || p > 2) {
 	return 0;
    }

    /* add OLS uhat squared to dataset */
    if (add_uhat_squared(pmod, scale, pZ, pdinfo)) {
	return E_ALLOC;
    }

    list = gretl_list_copy(garchlist);
    if (list == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    list[1] = (q > p)? q : p;
    list[2] = q;
    /* dep var is squared OLS residual */
    list[4] = v;

    amod = arma(list, (const double **) *pZ, pdinfo, OPT_NONE, NULL);
    if (amod.errcode) {
	err = amod.errcode;
	goto bailout;
    } else {
	int i;

	model_count_minus();
	garchpar_from_armapar(amod.coeff, q, p);

	for (i=0; i<q+p+1; i++) {
	    fprintf(stderr, "from ARMA: vparm_init[%d] = %#12.6g\n", i, 
		    vparm_init[i]);
	}
    }

 bailout:

    dataset_drop_last_variables(pdinfo->v - v, pZ, pdinfo);
    free(list);

    return err;
}

/* sanity/dimension check */

static int *get_garch_list (const int *list, int *err)
{
    int *ret = NULL;
    int i, p = list[1], q = list[2];
    int add0 = 1;

    *err = 0;

    /* rule out pure AR in variance */
    if (p > 0 && q == 0) {
	gretl_errmsg_set(_("Error in garch command"));
	*err = E_DATA;
	return NULL;
    }

    /* rule out excessive total GARCH terms */
    else if (p + q > 5) {
	gretl_errmsg_set(_("Error in garch command"));
	*err = E_DATA;
	return NULL;
    }

    /* insert constant if not present */
    for (i=4; i<=list[0]; i++) {
	if (list[i] == 0) {
	    add0 = 0;
	    break;
	}
    }

    ret = malloc((list[0] + 1 + add0) * sizeof *ret);
    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	ret[0] = list[0] + add0;
	for (i=1; i<=list[0]; i++) {
	    ret[i] = list[i];
	}
	if (add0) {
	    ret[i] = 0;
	}
    }

    return ret;
}

/* make regresson list for initial OLS */

static int *make_ols_list (const int *list)
{
    int *olist;
    int i;

    olist = malloc((list[0] - 2) * sizeof *olist);
    if (olist == NULL) {
	return NULL;
    }

    olist[0] = list[0] - 3;
    for (i=4; i<=list[0]; i++) {
	olist[i-3] = list[i];
    }

    return olist;
}

#if 0
int garch_pretest (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		   PRN *prn)
{
    int err;

    err = autocorr_test(pmod, pdinfo->pd, pZ, pdinfo,
			OPT_S | OPT_Q, prn);

    if (!err) {
	char teststat[64];
	double LMF = get_last_test_statistic(teststat);
	double pval = get_last_pvalue(teststat);

	fprintf(stderr, "LMF = %g, pval = %g\n", LMF, pval);
    }

    return err;
}
#endif

#define GARCH_SCALE_SIGMA 1

/* the driver function for the plugin */

MODEL garch_model (const int *cmdlist, double ***pZ, DATAINFO *pdinfo,
		   PRN *prn, gretlopt opt) 
{
    MODEL model;
    int *list = NULL, *ols_list = NULL;
    double scale = 1.0;
    int t, err, init_err, yno = 0;

    gretl_model_init(&model);

    list = get_garch_list(cmdlist, &err);
    if (err) {
	model.errcode = err;
    }

    if (!err) {
	ols_list = make_ols_list(list);
	if (ols_list == NULL) {
	    err = model.errcode = E_ALLOC;
	}
    }

    /* run initial OLS */
    if (!err) {
	model = lsq(ols_list, pZ, pdinfo, OLS, OPT_A | OPT_M, 0.0);
	if (model.errcode) {
	    err = model.errcode;
	}
    }

#if 0
    /* pretest the residuals for autocorrelation */
    if (!err) {
	err = garch_pretest(&model, pZ, pdinfo, prn);
    }
#endif

#if GARCH_SCALE_SDY
    if (!err) {
	yno = ols_list[1];
	scale = model.sdy;
	for (t=0; t<pdinfo->n; t++) {
	    if (!na((*pZ)[yno][t])) {
		(*pZ)[yno][t] /= scale;
	    }
	}
	for (t=0; t<model.ncoeff; t++) {
	    model.coeff[t] *= scale;
	}
	model.ess /= scale * scale;
	model.sigma = sqrt(model.ess / model.dfd);
    } 
#elif GARCH_SCALE_SIGMA
    if (!err) {
	yno = ols_list[1];
	scale = model.sigma;
	for (t=0; t<pdinfo->n; t++) {
	    if (!na((*pZ)[yno][t])) {
		(*pZ)[yno][t] /= scale;
	    }
	}
	for (t=0; t<model.ncoeff; t++) {
	    model.coeff[t] *= scale;
	}
	model.ess /= scale * scale;
	model.sigma = 1.0;
    } 
#endif /* alternative scalings */

    /* default variance parameter initialization */
    vparm_init[0] = model.sigma * model.sigma;
    for (t=1; t<VPARM_MAX; t++) {
	vparm_init[t] = 0.1;
    }

    if (opt & OPT_A) {
	/* "--arma-init": try initializing params via ARMA */
	init_err = garch_init_by_arma(&model, list, scale, 
				      pZ, pdinfo);
    }

    if (!err) {
	do_fcp(list, *pZ, scale, pdinfo, &model, prn, opt); 
    }

    if (scale != 1.0) {
	/* undo scaling of dependent variable */
	for (t=0; t<pdinfo->n; t++) {
	    if (!na((*pZ)[yno][t])) {
		(*pZ)[yno][t] *= scale;
	    }
	}
    }

    free(ols_list);
    free(list);

    return model;
}

