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
#include "version.h"
#include "libset.h"
#include "var.h"

#include "garch.h"

#define VPARM_DEBUG 0

#define PQ_MAX 7               /* max sum of GARCH p and q */
#define GARCH_PARAM_MAX 0.999

static void add_garch_varnames (MODEL *pmod, const DATASET *dset,
				const int *list)
{
    char tmp[16];
    int p = list[1];        /* GARCH beta terms */
    int q = list[2];        /* ARCH alpha terms > 0 */
    int r = list[0] - 4;    /* regressors */
    int np = 1 + p + q + r; /* the "1" is for alpha(0) */
    int i, j;

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    gretl_model_allocate_param_names(pmod, np);
    if (pmod->errcode) {
	return;
    }

    j = 0;

    for (i=0; i<r; i++) {
	gretl_model_set_param_name(pmod, j++, dset->varname[pmod->list[5+i]]);
    }

    gretl_model_set_param_name(pmod, j++, "alpha(0)");

    for (i=0; i<q; i++) {
	sprintf(tmp, "alpha(%d)", i + 1);
	gretl_model_set_param_name(pmod, j++, tmp);
    }

    for (i=0; i<p; i++) {
	sprintf(tmp, "beta(%d)", i + 1);
	gretl_model_set_param_name(pmod, j++, tmp);
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
	    sf = (j < nc)? scale * sfi : (j == nc)? sc2 * sfi : sfi;
	    vij = gretl_matrix_get(V, i, j) * sf;
	    gretl_matrix_set(V, i, j, vij);
	    gretl_matrix_set(V, j, i, vij);
	}
    }
}

static int 
write_garch_stats (MODEL *pmod, const int *list, const DATASET *dset,
		   double *theta, gretl_matrix *V, double scale,
		   const double *e, const double *h, 
		   int npar, int nc, int pad, int ifc, PRN *prn)
{
    double *garch_h;
    double den;
    const double *vcoef;
    int ynum = list[4];
    int nvp = list[1] + list[2];
    int xvars = list[0] - 4;
    int i, err;

    if (scale != 1.0) {
	rescale_results(theta, V, scale, npar, nc);
    }

    err = gretl_model_write_coeffs(pmod, theta, npar);

    if (!err) {
	gretl_model_write_vcv(pmod, V);
    }

    if (err) {
	return err;
    }

    /* verbose? */
    if (prn != NULL) {
	for (i=0; i<npar; i++) {
	    pprintf(prn, "theta[%d]: %#14.6g (%#.6g)\n", i, theta[i], 
		    pmod->sderr[i]);
	}
	pputc(prn, '\n'); 
    }   

    pmod->ess = 0.0;
    for (i=pmod->t1; i<=pmod->t2; i++) {
	pmod->uhat[i] = e[i + pad] * scale;
	pmod->ess += pmod->uhat[i] * pmod->uhat[i];
	pmod->yhat[i] = dset->Z[ynum][i] * scale - pmod->uhat[i];
    }

    vcoef = pmod->coeff + xvars;

    /* set sigma to its unconditional or steady-state value */
    den = 1.0;
    for (i=1; i<=nvp; i++) {
	den -= vcoef[i];
    }
    pmod->sigma = sqrt(vcoef[0] / den);

    pmod->adjrsq = NADBL; 
    pmod->fstt = NADBL;

    mle_criteria(pmod, 1);

    pmod->ci = GARCH;
    pmod->ifc = ifc;
    
    add_garch_varnames(pmod, dset, list);

    /* add predicted error variance to model */
    garch_h = malloc(dset->n * sizeof *garch_h);
    if (garch_h != NULL) {
	for (i=0; i<dset->n; i++) {
	    if (i < pmod->t1 || i > pmod->t2) {
		garch_h[i] = NADBL;
	    } else {
		garch_h[i] = h[i + pad] * scale * scale;
	    }
	}
	gretl_model_set_data(pmod, "garch_h", garch_h, 
			     GRETL_TYPE_DOUBLE_ARRAY,
			     dset->n * sizeof *garch_h);
    }

    return err;
}

static int make_garch_dataset (const int *list, const DATASET *dset,
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

    if (nx > 0) {
	if (pad) {
	    X = doubles_array_new(nx, bign);
	} else {
	    X = malloc(nx * sizeof *X);
	}

	if (X == NULL) {
	    free(y);
	    *py = NULL;
	    return E_ALLOC;
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
		s = t - pad;
		y[t] = dset->Z[vy][s];
		k = 5;
		for (i=0; i<nx; i++) {
		    vx = list[k++]; 
		    X[i][t] = dset->Z[vx][s];
		}
	    }
	}
    } else {
	/* build virtual dataset */
	*py = dset->Z[vy];
	k = 5;
	for (i=0; i<nx; i++) {
	    vx = list[k++]; 
	    X[i] = dset->Z[vx];
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
    if (vopt == ML_UNSET) {
	if (robust) {
	    if (ropt == ML_UNSET) {
		vopt = ML_QML;
	    } else {
		vopt = ropt;
	    }
	} else {
	    vopt = ML_HESSIAN;
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
    for (i=0; i<q; i++) {
	pprintf(prn, "  alpha[%d] = %g\n", i+1, theta[j++]);
    }
    for (i=0; i<p; i++) {
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
	/* initialization not done */
	return 0;
    }

    /* if we're _not_ using FCP, the following is handled 
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
garch_driver (const int *list, double scale,
	      const DATASET *dset, MODEL *pmod,
	      double *vparm, int ifc, gretlopt opt, 
	      PRN *prn)
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
    err = make_garch_dataset(list, dset, bign, pad, nc, &y, &X);
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
	    theta[i+nc] = vparm[i];
	}

	if (opt & OPT_V) {
	    garch_print_init(theta, nc, p, q, 0, prn);
	}
    }

    if ((opt & OPT_F) || libset_get_bool(USE_FCP)) {
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

    if (!err) {
	pmod->lnL = ll;
	write_garch_stats(pmod, list, dset, theta, V, scale, 
			  e, h, npar, nc, pad, ifc, prn);
	if (iters > 0) {
	    gretl_model_set_int(pmod, "iters", iters);
	} else if (grc > 0) {
	    gretl_model_set_int(pmod, "fncount", fnc);
	    gretl_model_set_int(pmod, "grcount", grc);
	} else {
	    gretl_model_set_int(pmod, "iters", fnc);
	}
	gretl_model_set_vcv_info(pmod, VCV_ML, vopt);
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

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }

    return err;
}

static int add_uhat_squared (const MODEL *pmod, double scale,
			     DATASET *dset)
{
    int t, v = dset->v;

    if (dataset_add_series(dset, 1)) {
	return E_ALLOC;
    }

    for (t=0; t<dset->n; t++) {
	double u = pmod->uhat[t];

	if (na(u)) {
	    dset->Z[v][t] = NADBL;
	} else {
	    u /= scale;
	    dset->Z[v][t] = u * u;
	}
    }

    strcpy(dset->varname[v], "uhat2");

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
garchpar_from_armapar (const double *armapar, int q, int p,
		       double *vparm)
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
	vparm[i] = (x < 0.0)? 0.01 : x;
	sum_ab += vparm[i];
    }

    for (i=1; i<=q; i++) {
	x = armapar[p+i];
	vparm[p+i] = (x > 0.0)? 0.0001 : -x;
	sum_ab += vparm[p+i];
    }

#if VPARM_DEBUG
    fprintf(stderr, "sum_ab = %#12.6g\n", sum_ab);
#endif

    if (sum_ab > GARCH_PARAM_MAX) {
	for (i=1; i<=p+q; i++) {
	    vparm[i] *= GARCH_PARAM_MAX / sum_ab;
	}
	sum_ab = GARCH_PARAM_MAX;
    }

    vparm[0] = armapar[0];
}

static int 
garch_init_by_arma (const MODEL *pmod, const int *glist, 
		    DATASET *dset, double scale, 
		    double *vparm)
{
    int p = glist[1], q = glist[2];
    int v = dset->v;
    int *list = NULL;
    int err = 0;

    /* for now we'll try this only for GARCH up to (2,2) */
    if (q > 2 || p > 2) {
 	return 0;
    }

    /* add OLS uhat squared to dataset */
    if (add_uhat_squared(pmod, scale, dset)) {
	return E_ALLOC;
    }

    list = gretl_list_copy(glist);

    if (list == NULL) {
	err = E_ALLOC;
    } else {
	MODEL amod;
	int i;

	list[1] = (q > p)? q : p;
	list[2] = q;
	/* dep var is squared OLS residual: last var added */
	list[4] = v;

	amod = arma(list, NULL, dset, OPT_C, NULL);
	err = amod.errcode;
	if (!err) {
	    model_count_minus();
	    garchpar_from_armapar(amod.coeff, p, q, vparm);
	    for (i=0; i<q+p+1; i++) {
		fprintf(stderr, "from ARMA: vparm_init[%d] = %#12.6g\n", i, 
			vparm[i]);
	    }
	}
	clear_model(&amod);
    }

    dataset_drop_last_variables(dset, dset->v - v);
    free(list);

    return err;
}

/* XPOS is the list position of the first regressor (if any).
   Note that the garch list structure is:

    0 1 2   3   4  5  6 ...
    # p q <sep> y x0 x1 ...

*/
#define XPOS 5

static int *get_garch_list (const int *list, const DATASET *dset, 
			    gretlopt opt, int *ifc, int *err)
{
    int *glist = NULL;
    int i, p, q;
    int cpos = 0;
    int add0 = 0;

    /* is the list well-formed? */
    if (list[0] < 4 || list[1] == LISTSEP ||
	list[2] == LISTSEP || list[3] != LISTSEP) {
	*err = E_PARSE;
	return NULL;
    }

    p = list[1];
    q = list[2];

    *err = 0;

    /* rule out pure AR in variance: the model is unidentified */
    if (p > 0 && q == 0) {
	gretl_errmsg_set(_("GARCH: p > 0 and q = 0: the model is unidentified"));
	*err = E_DATA;
	return NULL;
    }

    /* rule out excessive total GARCH-iness */
    if (p + q > PQ_MAX) {
	gretl_errmsg_sprintf(_("GARCH: p + q must not exceed %d"), PQ_MAX);
	*err = E_DATA;
	return NULL;
    }

    /* check for presence of constant among regressors */
    for (i=XPOS; i<=list[0]; i++) {
	if (list[i] == 0) {
	    /* got the constant: OK */
	    cpos = i;
	    break;
	}
    }

    /* OPT_N means don't auto-add a constant */
    if (cpos == 0 && !(opt & OPT_N)) {
	add0 = 1;
    } 

    *ifc = (cpos > 0 || add0);

    glist = gretl_list_new(list[0] + add0);

    if (glist == NULL) {
	*err = E_ALLOC;
    } else {
	int j = 1;

	/* transcribe first portion of original list */
	for (i=1; i<XPOS; i++) {
	    glist[j++] = list[i];
	}

	if (add0 || (cpos > 0 && cpos != XPOS)) {
	    /* insert constant here if not already present,
	       or if originally placed later */ 
	    glist[j++] = 0;
	} 

	/* transcribe the original regressors, if any */
	for (i=XPOS; i<=list[0]; i++) {
	    if (i == XPOS || list[i] != 0) {
		glist[j++] = list[i];
	    }
	}
    }

    return glist;
}

#define GARCH_AUTOCORR_TEST 1

#if GARCH_AUTOCORR_TEST

int garch_pretest (MODEL *pmod, DATASET *dset,
		   double *LMF, double *pvF)
{
    int err;

    err = autocorr_test(pmod, dset->pd, dset,
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

#if GARCH_SCALE_SIGMA

static double garch_scale_sigma (int yno, MODEL *pmod, DATASET *dset)
{
    double scale = pmod->sigma;
    int i;

    for (i=0; i<dset->n; i++) {
	if (!na(dset->Z[yno][i])) {
	    dset->Z[yno][i] /= scale;
	}
    }

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] /= scale;
    }

    pmod->ess /= scale * scale;
    pmod->sigma = 1.0;

    return scale;
}

static void garch_undo_scaling (int yno, double scale, DATASET *dset)
{
    int t;

    if (scale != 1.0) {
	for (t=0; t<dset->n; t++) {
	    if (!na(dset->Z[yno][t])) {
		dset->Z[yno][t] *= scale;
	    }
	}
    }
}

#endif

/* default variance parameter initialization */

static void garch_vparm_init (const int *list, double sigma,
			      double *vparm)
{
    int i, q = list[1], p = list[2];
    double den = 1.0;
    double tmp = (q>0) ? 0.2 : 0.8;

    if (p > 0) {
	for (i=1; i<=p; i++) {
	    vparm[i] = tmp / p;
	    den -= vparm[i];
	}
    }

    if (q > 0) {
	for (i=p+1; i<=p+q; i++) {
	    vparm[i] = 0.7 / q;
	    den -= vparm[i];
	}
    }

    vparm[0] = sigma * sigma * den;
}

/* make regression list for initial OLS: we skip three terms 
   from the GARCH list, namely p, q and the separator 
   that follows.
*/

static int *make_ols_list (const int *list, int *err)
{
    int *olist;
    int i;

    olist = gretl_list_new(list[0] - 3);

    if (olist == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=4; i<=list[0]; i++) {
	    olist[i-3] = list[i];
	}
    }

    return olist;
}

static MODEL garch_run_ols (const int *list, DATASET *dset, 
			    PRN *prn)
{
    int *ols_list;
    MODEL model;
    int err = 0;

    ols_list = make_ols_list(list, &err);
    if (err) {
	gretl_model_init(&model, NULL);
	model.errcode = err;
	return model;
    }

    model = lsq(ols_list, dset, OLS, OPT_A | OPT_M | OPT_U);

#if 0
    fprintf(stderr, "errcode=%d, ess=%g, sigma=%g\n",
	    model.errcode, model.ess, model.sigma);
    if (!model.errcode) {
	printmodel(&model, dset, OPT_NONE, prn);
    }
#endif

    free(ols_list);

    if (!model.errcode) {
	clear_model_xpx(&model);
    }

    return model;
}

static void clean_dropped_vars (MODEL mod, int *list)
{
    if (list[0] - mod.list[0] > 3) {
	int i;
        list[0] = mod.list[0] + 3;
	for (i=4; i<=list[0]; i++) {
	    list[i] = mod.list[i-3];
	}
    }
}

static void garch_add_lr_test (MODEL *pmod, double llr,
			       const int *list)
{
    if (!na(pmod->lnL) && llr <= pmod->lnL) {
	double LR = 2.0 * (pmod->lnL - llr);
	int LRdf = list[1] + list[2];

	gretl_model_set_double(pmod, "garch_LR", LR);
	gretl_model_set_int(pmod, "garch_LR_df", LRdf);
    }
}

static void garch_standardize_residuals (MODEL *pmod)
{
    double *h = gretl_model_get_data(pmod, "garch_h");

    if (h != NULL) {
	int t;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    pmod->uhat[t] /= sqrt(h[t]);
	}
	pmod->opt |= OPT_Z;
    }
}

/* the driver function for the plugin */

MODEL garch_model (const int *cmdlist, DATASET *dset,
		   PRN *prn, gretlopt opt) 
{
    MODEL model;
    int *list = NULL;
    double vparm[PQ_MAX+1] = {0};
    double LMF = NADBL;
    double pvF = NADBL;
    double llr = NADBL;
    double scale = 1.0;
    int ols_T, ifc, yno = 0;
    int err = 0;

    list = get_garch_list(cmdlist, dset, opt, &ifc, &err);
    if (err) {
	gretl_model_init(&model, NULL);
	model.errcode = err;
	return model;
    }

    /* run initial OLS */
    model = garch_run_ols(list, dset, prn);
    if (model.errcode) {
	free(list);
	return model;
    }

    clean_dropped_vars(model, list);
    llr = model.lnL;
    ols_T = model.nobs;

#if GARCH_AUTOCORR_TEST
    /* pretest the residuals for autocorrelation */
    if (prn != NULL) {
	garch_pretest(&model, dset, &LMF, &pvF);
    }
#endif

#if GARCH_SCALE_SIGMA
    yno = list[4];
    scale = garch_scale_sigma(yno, &model, dset);
#endif 

    /* variance parameter initialization */
    garch_vparm_init(list, model.sigma, vparm);

    if (opt & OPT_A) {
	/* "--arma-init": try initializing params via ARMA */
	garch_init_by_arma(&model, list, dset, scale, vparm);
    }

    garch_driver(list, scale, dset, &model, vparm,
		 ifc, opt, prn); 

#if GARCH_SCALE_SIGMA
    garch_undo_scaling(yno, scale, dset);
#endif 

#if GARCH_AUTOCORR_TEST
    if (!na(LMF)) {
	if (model.errcode == E_NOCONV) {
	    autocorr_message(LMF, pvF, dset->pd, prn);
	} else {
	    gretl_model_destroy_tests(&model);
	}
    }
#endif

    if (!model.errcode) {
	if (opt & OPT_Z) {
	    garch_standardize_residuals(&model);
	}
	if (!na(llr) && ols_T == model.nobs) {
	    garch_add_lr_test(&model, llr, cmdlist);
	}
    } 

    free(list);

    return model;
}
