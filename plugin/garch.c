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

/* GARCH plugin for gretl using the Fiorentini, Calzolari and 
   Panattoni mixed-gradient algorithm.
*/

#include "libgretl.h"
#include "libset.h"
#include "var.h"

#include "garch.h"

#define VPARM_DEBUG 0

#define VPARM_MAX 6            /* max number of variance parameters */
#define GARCH_PARAM_MAX 0.999

double vparm_init[VPARM_MAX];

static void add_garch_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				const int *list)
{
    int p = list[1];
    int q = list[2];
    int r = list[0] - 4;
    int i, j, np = 2 + p + q + r;

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    gretl_model_allocate_params(pmod, np);
    if (pmod->errcode) {
	return;
    }

    strcpy(pmod->params[0], pdinfo->varname[0]);

    j = 1;
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

static void rescale_results (double *theta, gretl_matrix *V,
			     double scale, int npar, int nc)
{
    double vij, sfi, sf, sc2 = scale * scale;
    int i, j;

    for (i=0; i<nc; i++) {
	theta[i] *= scale;
    }

    theta[nc] *= sc2;

    for (i=0; i<npar; i++) {
	sfi = (i < nc)? scale : (i == nc)? sc2 : 1.0;
	for (j=0; j<=i; j++) {
	    sf = (j < nc)? scale*sfi : (j == nc)? sc2*sfi : sfi;
	    vij = gretl_matrix_get(V, i, j) * sf;
	    gretl_matrix_set(V, i, j, vij);
	    gretl_matrix_set(V, j, i, vij);
	}
    }
}

static int 
write_garch_stats (MODEL *pmod, const int *list, 
		   const double **Z, const DATAINFO *pdinfo,
		   double *theta, gretl_matrix *V, double scale,
		   const double *e, const double *h, 
		   int npar, int nc, int pad, PRN *prn)
{
    int err = 0;
    double *coeff, *sderr, *vcv, *garch_h;
    double x, den;
    int ynum = list[4];
    int nvp = list[1] + list[2];
    int xvars = list[0] - 4;
    int nv = npar * (npar + 1) / 2;
    int i, j, k;

    coeff = realloc(pmod->coeff, npar * sizeof *pmod->coeff);
    sderr = realloc(pmod->sderr, npar * sizeof *pmod->sderr);
    vcv = realloc(pmod->vcv, nv * sizeof *pmod->vcv);

    if (coeff == NULL || sderr == NULL || vcv == NULL) {
	return E_ALLOC;
    }

    if (scale != 1.0) {
	rescale_results(theta, V, scale, npar, nc);
    }

    for (i=0; i<npar; i++) {
	coeff[i] = theta[i];
	x = gretl_matrix_get(V, i, i);
	sderr[i] = (x > 0.0)? sqrt(x) : 0.0;
	for (j=0; j<=i; j++) {
	    k = ijton(i, j, npar);
	    vcv[k] = gretl_matrix_get(V, i, j);
	}
    }

    /* verbose? */
    if (prn != NULL) {
	for (i=0; i<npar; i++) {
	    pprintf(prn, "theta[%d]: %#14.6g (%#.6g)\n", i, theta[i], 
		    sderr[i]);
	}
	pputc(prn, '\n'); 
    }   

    pmod->coeff = coeff;
    pmod->sderr = sderr;
    pmod->vcv = vcv;
    pmod->ncoeff = npar;

    pmod->ess = 0.0;
    for (i=pmod->t1; i<=pmod->t2; i++) {
	pmod->uhat[i] = e[i + pad] * scale;
	pmod->ess += pmod->uhat[i] * pmod->uhat[i];
	pmod->yhat[i] = Z[ynum][i] * scale - pmod->uhat[i];
    }

    /* set sigma to its unconditional or steady-state value */
    den = 1.0;
    for (i=0; i<nvp; i++) {
	den -= coeff[i+xvars+1];
    }
    pmod->sigma = sqrt(coeff[xvars] / den);

    pmod->adjrsq = NADBL; 
    pmod->fstt = NADBL;

    mle_criteria(pmod, 1);

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
			     GRETL_TYPE_DOUBLE_ARRAY,
			     pdinfo->n * sizeof *garch_h);
    }

    return err;
}

static int make_garch_dataset (const int *list, double **Z,
			       int bign, int pad, int nx,
			       double **py, double ***pX)
{
    double *y = NULL, **X = NULL;
    int vx, vy = list[4];
    int i, k, s, t;

    /* If pad > 0 we have to create a newly allocated, padded
       dataset.  Otherwise we can use a virtual dataset, made
       up of pointers into the original dataset, Z. 
    */

    if (pad > 0) {
	y = malloc(bign * sizeof *y);
	if (y == NULL) {
	    return E_ALLOC;
	}
	*py = y;
    } 

    if (pad) {
	X = doubles_array_new(nx, bign);
    } else {
	X = malloc(nx * sizeof *X);
    }

    if (X == NULL) {
	free(y);
	return E_ALLOC;
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
		s = t - pad;
		y[t] = Z[vy][s];
		k = 5;
		for (i=0; i<nx; i++) {
		    vx = list[k++]; 
		    X[i][t] = Z[vx][s];
		}
	    }
	}
    } else {
	/* build virtual dataset */
	*py = Z[vy];
	k = 5;
	for (i=0; i<nx; i++) {
	    vx = list[k++]; 
	    X[i] = Z[vx];
	}
    }

    *pX = X;

    return 0;
}

static int get_vopt (int robust)
{
    int vopt = libset_get_int(GARCH_VCV);
    int ropt = libset_get_int(GARCH_ROBUST_VCV);

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

static void garch_print_init (const double *theta, int k,
			      int p, int q, int manual, 
			      PRN *prn)
{
    int i, j = 0;

    pputc(prn, '\n');

    if (manual) {
	pputs(prn, "Manual initialization of parameters");
    } else {
	pputs(prn, "Automatic initialization of parameters");
    }

    pputs(prn, "\n\n Regression coefficients:\n");

    for (i=0; i<k; i++) {
	pprintf(prn, "  theta[%d] = %g\n", i, theta[j++]);
    }

    pputs(prn, "\n Variance parameters:\n");

    pprintf(prn, "  alpha[0] = %g\n", theta[j++]);
    for (i=0; i<p; i++) {
	pprintf(prn, "  alpha[%d] = %g\n", i+1, theta[j++]);
    }
    for (i=0; i<q; i++) {
	pprintf(prn, "   beta[%d] = %g\n", i, theta[j++]);
    }

    pputc(prn, '\n');
}

/* pick up any manually set initial values (if these
   have been set via "set initvals") */

static int garch_manual_init (double *theta, int k, int p, int q, 
			      PRN *prn)
{
    int mlen = n_init_vals();
    int n = k + p + q + 1;

    if (mlen != n) {
	if (mlen > 0) {
	    fprintf(stderr, "Number of initvals = %d, but we want %d "
		    "values for GARCH\n", mlen, n);
	}
	return 0;
    }

    /* if we're not using FCP, the following is handled 
       within the BFGS routine */

    if (libset_get_bool("fcp")) {
	const gretl_matrix *m = get_init_vals();
	int i;
	
	/* order: coeffs on regressors; variance params */

	for (i=0; i<n; i++) {
	    theta[i] = m->val[i];
	}

	garch_print_init(theta, k, p, q, 1, prn);

	free_init_vals();
    }

    return 1;
}

static int 
garch_driver (const int *list, double **Z, double scale,
	      const DATAINFO *pdinfo, MODEL *pmod,
	      PRN *prn, gretlopt opt)
{
    int t1 = pmod->t1, t2 = pmod->t2;
    int nc = pmod->ncoeff;
    int p = list[1];
    int q = list[2];
    double *y = NULL;
    double **X = NULL;
    double *h = NULL;
    double *e = NULL, *e2 = NULL;
    double *theta = NULL;
    double ll = NADBL;
    gretl_matrix *V = NULL;
    int fnc = 0, grc = 0, iters = 0;
    int nobs, maxlag, bign, pad = 0;
    int i, npar, vopt;
    int err = 0;

    vopt = get_vopt(opt & OPT_R);

    maxlag = (p > q)? p : q; 
    npar = nc + p + q + 1;

    nobs = t2 + 1; /* number of obs in full dataset */

    if (maxlag > t1) {
	/* need to pad data series at start */
	pad = maxlag - t1;
    } 

    /* length of series to pass to garch_estimate */
    bign = nobs + pad;
	
    e = malloc(bign * sizeof *e);
    e2 = malloc(bign * sizeof *e2);
    h = malloc(bign * sizeof *h);
    if (e == NULL || e2 == NULL || h == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<bign; i++) {
	e[i] = e2[i] = h[i] = 0.0;
    }   
 
    theta = malloc(npar * sizeof *theta);
    if (theta == NULL) {
	err = E_ALLOC;
	goto bailout;	
    }

    V = gretl_zero_matrix_new(npar, npar);
    if (V == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* create dataset for garch estimation */
    err = make_garch_dataset(list, Z, bign, pad, nc, &y, &X);
    if (err) {
	goto bailout;
    }

    if (!garch_manual_init(theta, nc, p, q, prn)) {

	/* initial coefficients from OLS */
	for (i=0; i<nc; i++) {
	    theta[i] = pmod->coeff[i];
	}

	/* initialize variance parameters */
	for (i=0; i<p+q+1; i++) {
	    theta[i+nc] = vparm_init[i];
	}

	if (opt & OPT_V) {
	    garch_print_init(theta, nc, p, q, 0, prn);
	}
    }

    if (libset_get_bool("fcp")) {
	err = garch_estimate(y, (const double **) X,
			     t1 + pad, t2 + pad, bign, nc, 
			     p, q, theta, V, e, e2, h,
			     scale, &ll, &iters, vopt, prn);
    } else {
	err = garch_estimate_mod(y, (const double **) X,
				 t1 + pad, t2 + pad, bign, nc, 
				 p, q, theta, V, e, e2, h, 
				 scale, &ll, &fnc, &grc, vopt, prn);
    }

    if (err != 0) {
	pmod->errcode = err;
    } else {
	pmod->lnL = ll;
	write_garch_stats(pmod, list, (const double **) Z, pdinfo, 
			  theta, V, scale, e, h, npar, nc, pad,
			  prn);
	if (iters > 0) {
	    gretl_model_set_int(pmod, "iters", iters);
	} else if (fnc > 0) {
	    gretl_model_set_int(pmod, "fncount", fnc);
	    gretl_model_set_int(pmod, "grcount", grc);
	}
	gretl_model_set_int(pmod, "ml_vcv", vopt);
    }

 bailout:

    free(e);
    free(e2);
    free(h);
    free(theta);
    gretl_matrix_free(V);

    if (pad > 0) {
	free(y);
	doubles_array_free(X, nc);
    } else {
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

#if VPARM_DEBUG
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

#if VPARM_DEBUG
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

    amod = arma(list, NULL, (const double **) *pZ, pdinfo, OPT_C, NULL);
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

static int *get_garch_list (const int *list, const double **Z,
			    const DATAINFO *pdinfo, int *err)
{
    int *glist = NULL;
    int i, p = list[1], q = list[2];
    int add0 = 1;

    *err = 0;

    /* rule out pure AR in variance: the model is unidentified */
    if (p > 0 && q == 0) {
	gretl_errmsg_set(_("Error in garch command"));
	*err = E_DATA;
	return NULL;
    }

    /* rule out excessive total GARCH terms */
    if (p + q > 5) {
	gretl_errmsg_set(_("Error in garch command"));
	*err = E_DATA;
	return NULL;
    }

    /* insert constant if not present */
    for (i=4; i<=list[0]; i++) {
	if (list[i] == 0 || true_const(list[i], Z, pdinfo)) {
	    add0 = 0;
	    break;
	}
    }

    glist = gretl_list_new(list[0] + add0);
    if (glist == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=1; i<=4; i++) {
	    glist[i] = list[i];
	}

	if (add0) {
	    glist[5] = 0;
	}

	for (i=5; i<=list[0]; i++) {
	    glist[i+add0] = list[i];
	}
    }

    return glist;
}

/* make regression list for initial OLS */

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

#define GARCH_AUTOCORR_TEST 1

#if GARCH_AUTOCORR_TEST

int garch_pretest (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		   double *LMF, double *pvF)
{
    int err;

    err = autocorr_test(pmod, pdinfo->pd, pZ, pdinfo,
			OPT_S | OPT_Q, NULL);

    if (!err) {
	*LMF = get_last_test_statistic(NULL);
	*pvF = get_last_pvalue(NULL);
    } 

    return err;
}

static void autocorr_message (double LMF, double pvF, int order, PRN *prn)
{
    if (!na(LMF) && pvF < 0.05) {
	pputs(prn, "\nConvergence was not reached.  One possible reason "
	      "for this is\nautocorrelation in the error term.\n");
	pprintf(prn, "After estimating the model by OLS, the following result "
		"was\nobtained for a test of autocorrelation of order %d:\n",
		order);
	pprintf(prn, "LMF = %g, with p-value %g\n", LMF, pvF);
    }
}

#endif /* GARCH_AUTOCORR_TEST */

#define GARCH_SCALE_SIGMA 1

/* the driver function for the plugin */

MODEL garch_model (const int *cmdlist, double ***pZ, DATAINFO *pdinfo,
		   PRN *prn, gretlopt opt) 
{
    MODEL model;
    int *list = NULL;
    int *ols_list = NULL;
    double LMF = NADBL;
    double pvF = NADBL;
    double scale = 1.0;
    int i, t, err, init_err, yno = 0;

    gretl_model_init(&model);

    list = get_garch_list(cmdlist, (const double **) *pZ,
			  pdinfo, &err);
    if (err) {
	model.errcode = err;
	return model;
    }

    ols_list = make_ols_list(list);
    if (ols_list == NULL) {
	model.errcode = E_ALLOC;
	return model;
    }

    /* run initial OLS */
    model = lsq(ols_list, pZ, pdinfo, OLS, OPT_A | OPT_M);
    if (model.errcode) {
	goto bailout;
    }

#if 0
    printmodel(&model, pdinfo, OPT_NONE, prn);
#endif

#if GARCH_AUTOCORR_TEST
    /* pretest the residuals for autocorrelation */
    if (prn != NULL) {
	garch_pretest(&model, pZ, pdinfo, &LMF, &pvF);
    }
#endif

#if GARCH_SCALE_SIGMA
    yno = ols_list[1];
    scale = model.sigma;
    for (t=0; t<pdinfo->n; t++) {
	if (!na((*pZ)[yno][t])) {
	    (*pZ)[yno][t] /= scale;
	}
    }
    for (t=0; t<model.ncoeff; t++) {
	model.coeff[t] /= scale;
    }
    model.ess /= scale * scale;
    model.sigma = 1.0;
#endif 

    /* default variance parameter initialization */
    int q = list[1], p = list[2];

    vparm_init[0] = model.sigma * model.sigma * 0.1;
    if (q>0) {
	vparm_init[p+1] = 0.7;
    }
    vparm_init[1] = (q==0 ? 0.9 : 0.2) / p;
    for (i=2; i<=p; i++) {
	vparm_init[i] = vparm_init[1];
    }

    if (opt & OPT_A) {
	/* "--arma-init": try initializing params via ARMA */
	init_err = garch_init_by_arma(&model, list, scale, 
				      pZ, pdinfo);
    }

    garch_driver(list, *pZ, scale, pdinfo, &model, prn, opt); 

    if (scale != 1.0) {
	/* undo scaling of dependent variable */
	for (t=0; t<pdinfo->n; t++) {
	    if (!na((*pZ)[yno][t])) {
		(*pZ)[yno][t] *= scale;
	    }
	}
    }

 bailout:

    free(ols_list);
    free(list);

#if GARCH_AUTOCORR_TEST
    if (!na(LMF)) {
	if (model.errcode == E_NOCONV) {
	    autocorr_message(LMF, pvF, pdinfo->pd, prn);
	} else {
	    gretl_model_destroy_tests(&model);
	}
    }
#endif

    return model;
}

