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

#include "libgretl.h"
#include "libset.h"

#define INTDEBUG 0

#define INTERVAL_TOL 1.0e-12

enum {
    INT_LOW, /* no lower bound */
    INT_MID, /* both bounds */
    INT_HIGH, /* no upper bound */
    INT_POINT, /* point value */
};

typedef struct int_container_ int_container;

struct int_container_ {
    int hiv, lov;     /* variable numbers for limit series */
    double *hi, *lo;  /* limit series for dependent variable */
    int *obstype;     /* dependent variable classifier */
    double **X;       /* regressors */
    int *list;        /* dependent vars plus regular regressors */
    int t1;           /* beginning of sample */
    int t2;           /* end of sample */
    int nobs;         /* number of observations */
    int nx;           /* number of explanatory variables */
    int k;            /* total number of parameters */
    double *theta;    /* real parameter estimates */
    double *ndx;      /* index variable */
    double *dP;       /* probabilities */
    double *df;       /* densities */
    MODEL *pmod;      /* model struct, initially containing OLS */
    double **G;       /* score matrix by observation */
};

static void int_container_destroy (int_container *ICont)
{
    free(ICont->hi);
    free(ICont->lo);
    free(ICont->obstype);
    doubles_array_free(ICont->X, ICont->nx);
    free(ICont->list);
    free(ICont->theta);
    free(ICont->ndx);
    free(ICont->dP);
    free(ICont->df);
    doubles_array_free(ICont->G, ICont->k);

    free(ICont);
}

static int_container *int_container_new (int *list, double **Z, 
					 DATAINFO *pdinfo, MODEL *mod)
{
    int_container *ICont;
    int nobs, i, s, t, vlo = list[1], vhi = list[2];

    ICont = malloc(sizeof *ICont);
    if (ICont == NULL) {
	return NULL;
    }

    ICont->lov = vlo;
    ICont->hiv = vhi;
    ICont->pmod = mod;

    ICont->t1 = mod->t1;
    ICont->t2 = mod->t2;
    nobs = ICont->nobs = mod->nobs;
    ICont->nx = mod->ncoeff;
    ICont->k = mod->ncoeff + 1; /* beta + sigma */

    ICont->hi = NULL;
    ICont->lo = NULL;
    ICont->obstype = NULL;
    ICont->X = NULL;
    ICont->list = NULL;
    ICont->theta = NULL;
    ICont->ndx = NULL;
    ICont->dP = NULL;
    ICont->df = NULL;
    ICont->G = NULL;

    ICont->hi = malloc(nobs * sizeof *ICont->hi);
    ICont->lo = malloc(nobs * sizeof *ICont->lo);
    ICont->obstype = malloc(nobs * sizeof *ICont->obstype);
    ICont->X = doubles_array_new(ICont->nx, nobs);
    ICont->list = gretl_list_new(2 + ICont->nx);
    ICont->theta = malloc(ICont->k * sizeof *ICont->theta);

    ICont->ndx = malloc(nobs * sizeof *ICont->ndx);
    ICont->dP = malloc(nobs * sizeof *ICont->dP);
    ICont->df = malloc(nobs * sizeof *ICont->dP);
    ICont->G = doubles_array_new(ICont->k, nobs);

    if (ICont->lo == NULL || ICont->hi == NULL || 
	ICont->obstype == NULL || ICont->X == NULL || 
	ICont->list == NULL || ICont->theta == NULL || 
	ICont->ndx == NULL || ICont->dP == NULL || 
	ICont->df == NULL || ICont->G == NULL) {
	int_container_destroy(ICont);
	return NULL;
    }

    s = 0;
    double x0, x1;
    for (t=mod->t1; t<=mod->t2; t++) {
	if (!na(mod->uhat[t])) {

	    x0 = Z[vlo][t];
	    x1 = Z[vhi][t];
	    ICont->lo[s] = x0;
	    ICont->hi[s] = x1;

	    if (na(x0)) {
		ICont->obstype[s] = INT_LOW;
	    } else if (na(x1)) {
		ICont->obstype[s] = INT_HIGH;
	    } else if (x0==x1) {
		ICont->obstype[s] = INT_POINT;
	    } else {
		ICont->obstype[s] = INT_MID;
	    }

	    for (i=0; i<ICont->nx; i++) {
		ICont->X[i][s] = Z[list[i+3]][t];
	    }

	    s++;
	}
    }

    ICont->list[1] = vlo;
    ICont->list[2] = vhi;
    for (i=0; i<ICont->nx; i++) {
	ICont->list[i+3] = mod->list[i+3];
    }

    for (i=0; i<ICont->nx; i++) {
	ICont->theta[i] = mod->coeff[i];
    }
    ICont->theta[i] = log(mod->sigma);
    

#if INTDEBUG > 1
    fprintf(stderr, "nobs = %d\n", ICont->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", ICont->t1, ICont->t2);
    fprintf(stderr, "k = %d\n", ICont->k);
    fprintf(stderr, "nx = %d\n", ICont->nx);
#endif

    return ICont;
}


static int create_midpoint_y (int *list, double ***pZ, DATAINFO *pdinfo, 
			      int **initlist)
{
    int i, t;
    int n = list[0];
    int mpy = pdinfo->v;
    int err = dataset_add_series(1, pZ, pdinfo);

    if (err) {
	return err;
    }

    int lv, hv;
    double x0, x1;

    lv = list[1];
    hv = list[2];

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	x0 = (*pZ)[lv][t];
	x1 = (*pZ)[hv][t];

	if (na(x0)) {
	    (*pZ)[mpy][t] = x1;
	} else if (na(x1)) {
	    (*pZ)[mpy][t] = x0;
	} else {
	    (*pZ)[mpy][t] = 0.5*(x0 + x1);
	}
    }

    *initlist = gretl_list_new(n-1);
    (*initlist)[1] = mpy;
    for (i=3; i<=n; i++) {
	(*initlist)[i-1] = list[i];
    }

#if INTDEBUG 
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	fprintf(stderr, "%2d: %g\n", t, (*pZ)[mpy][t]);
    }
#endif

    return 0;
}

static double int_loglik (const double *theta, void *ptr)
{
    int_container *IC = (int_container *) ptr;
    double x0, x1, ndx, ll = 0.0;
    int err;
    int i, s, t, v;
    int k = IC->k;

    double sigma = exp(theta[k-1]);

    for (t=0; t<IC->nobs; t++) {
	x0 = IC->lo[t];
	x1 = IC->hi[t];

	ndx = 0;
	for (i=0; i<IC->nx; i++) {
	    ndx += theta[i] * IC->X[i][t];
	}

	IC->ndx[t] = ndx;

	switch (IC->obstype[t]) {
	case INT_LOW:
	    IC->dP[t] = normal_cdf((x1 - ndx)/sigma);
	    break;
	case INT_HIGH:
	    IC->dP[t] = 1.0 - normal_cdf((x0 - ndx)/sigma);
	    break;
	case INT_MID:
	    IC->dP[t] = normal_cdf((x1 - ndx)/sigma) - normal_cdf((x0 - ndx)/sigma);
	    break;
	case INT_POINT:
	    IC->dP[t] = normal_pdf((x0 - ndx)/sigma)/sigma;
	}

	ll += log(IC->dP[t]);
    }
    
#if INTDEBUG > 1
    fprintf(stderr, "ll = %16.10f\n", ll);
#endif

    return ll;
}

static void fill_model(int_container *IC, double *hess)
{
    int i, j, n;
    int k = IC->k;

    IC->pmod->ci = INTREG;

    for (i=0; i<IC->nx; i++) {
	IC->pmod->coeff[i] = IC->theta[i];
    }	

    double sigma = exp(IC->theta[k-1]);
    IC->pmod->sigma = sigma;

    IC->pmod->lnL = int_loglik(IC->theta, IC);

    n = 0;

    for (i=0; i<IC->nx; i++) {
	IC->pmod->sderr[i] = sqrt(hess[n]);
	n += (k-i);
    }

    IC->pmod->ybar = NADBL;
    IC->pmod->sdy = NADBL;

    gretl_model_set_int(IC->pmod, "lovar", IC->lov);
    gretl_model_set_int(IC->pmod, "hivar", IC->hiv);
}

static int do_interval (int *list, double **Z, DATAINFO *pdinfo, 
			MODEL *mod, gretlopt opt, PRN *prn) 
{
    int i, err;
    int_container *IC = NULL;
    double *hess; 

    IC = int_container_new(list, Z, pdinfo, mod);

    double Loglik;
    Loglik = int_loglik(IC->theta, IC);

    int fncount, grcount;

    err = BFGS_max(IC->theta, IC->k, 1000, INTERVAL_TOL, 
		   &fncount, &grcount, int_loglik, C_LOGLIK,
		   NULL, IC, (0 && prn != NULL)? OPT_V : OPT_NONE,
		   prn);

    if (!err) {
	hess = numerical_hessian(IC->theta, IC->k, int_loglik, IC, &err);
    }

    if (!err) {
	fill_model(IC, hess);
    }

    int_container_destroy(IC);
    return 0;
}

MODEL interval_estimate (int *list, double ***pZ, DATAINFO *pdinfo,
			 gretlopt opt, PRN *prn) 
{
    int err;
    MODEL model;
    int *initlist;

    gretl_model_init(&model);

    /* 
       create extra variable for model initialization and
       corresponding list
    */
    err = create_midpoint_y(list, pZ, pdinfo, &initlist);
    if (err) {
	return model;
    }

#if INTDEBUG
    printlist(initlist, "init list");
#endif

    /* run initial OLS */
    model = lsq(initlist, pZ, pdinfo, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "interval_estimate: initial OLS failed\n");
	free(initlist);
	return model;
    }

#if INTDEBUG
    pprintf(prn, "interval_estimate: initial OLS\n");
    printmodel(&model, pdinfo, OPT_S, prn);
#endif

    /* do the actual analysis */
    model.errcode = do_interval(list, (double **) *pZ, pdinfo, &model, opt, prn);

    /* clean up midpoint-y */
    dataset_drop_last_variables(1, pZ, pdinfo);

    free(initlist);

    return model;
}
