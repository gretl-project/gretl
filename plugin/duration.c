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
#include "matrix_extra.h"
#include "gretl_bfgs.h"

#include <errno.h>

typedef struct duration_info_ duration_info;

enum {
    SCORE_UPDATE_XB = 1
};

struct duration_info_ {
    int flags;             /* control info */
    double ll;             /* loglikelihood */
    int k;                 /* number of regressors */
    int T;                 /* number of observations */
    double *theta;         /* params array, length k + 1 */
    gretl_matrix_block *B; /* workspace */
    gretl_vector *y;       /* dependent variable */
    gretl_matrix *X;       /* regressors */
    gretl_matrix *beta;    /* coeffs on regressors */
    gretl_matrix *llt;     /* per-observation likelihood */
    gretl_matrix *Xb;      /* X \times \beta */
    gretl_matrix *G;       /* score */
    gretl_matrix *V;       /* covariance matrix */
    gretl_vector *cens;    /* censoring variable */
    PRN *prn;              /* verbose printer */
};

static void duration_free (duration_info *dinfo)
{
    gretl_matrix_block_destroy(dinfo->B);
    free(dinfo->theta);
    gretl_matrix_free(dinfo->cens);
}

static int duration_init (duration_info *dinfo, MODEL *pmod,
			  int censvar, const double **Z, 
			  gretlopt opt, PRN *prn)
{
    int n = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t, v;

    dinfo->flags = 0;

    dinfo->B = NULL;
    dinfo->cens = NULL;

    dinfo->theta = malloc((k + 1) * sizeof *dinfo->theta);
    if (dinfo->theta == NULL) {
	return E_ALLOC;
    }

    if (censvar > 0) {
	dinfo->cens = gretl_column_vector_alloc(n);
	if (dinfo->cens == NULL) {
	    return E_ALLOC;
	}
    }

    dinfo->B = gretl_matrix_block_new(&dinfo->y, n, 1,
				      &dinfo->X, n, k,
				      &dinfo->beta, k, 1,
				      &dinfo->Xb, n, 1,
				      &dinfo->llt, n, 1,
				      &dinfo->G, n, k+1,
				      &dinfo->V, k+1, k+1,
				      NULL);
    if (dinfo->B == NULL) {
	return E_ALLOC;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	v = pmod->list[1];
	dinfo->y->val[s] = Z[v][t];
	if (dinfo->cens != NULL) {
	    dinfo->cens->val[s] = Z[censvar][t];
	}
	for (i=0; i<k; i++) {
	    v = pmod->list[i+2];
	    gretl_matrix_set(dinfo->X, s, i, Z[v][t]);
	}
	s++;
    }

    dinfo->k = k;
    dinfo->T = n;

    for (i=0; i<k; i++) {
	/* initialize using OLS estimates */
	dinfo->theta[i] = pmod->coeff[i];
    }
    dinfo->theta[k] = 1.0; /* FIXME smart s initialization? */

    dinfo->ll = NADBL;
    dinfo->prn = (opt & OPT_V)? prn : NULL;

    return 0;
}

static void duration_update_Xb (duration_info *dinfo, const double *theta)
{
    int i;

    for (i=0; i<dinfo->k; i++) {
	dinfo->beta->val[i] = theta[i];
    }

    gretl_matrix_multiply(dinfo->X, dinfo->beta, dinfo->Xb);
}

static double duration_loglik (const double *theta, void *data)
{
    duration_info *dinfo = (duration_info *) data;
    double s = theta[dinfo->k];
    double *ll = dinfo->llt->val;
    double *Xb = dinfo->Xb->val;
    double *y = dinfo->y->val;
    double lns, zt;
    int t;

    if (s <= 0) {
	return NADBL;
    }

    duration_update_Xb(dinfo, theta);

    dinfo->ll = 0.0;
    errno = 0;
    lns = log(s);

    /* FIXME censoring */

    for (t=0; t<dinfo->T; t++) {
	zt = (y[t] - Xb[t]) / s;
	ll[t] = zt - lns - exp(zt);
	dinfo->ll += ll[t];
    }

    if (errno) {
	dinfo->ll = NADBL;
    }

    return dinfo->ll;
}

static int duration_score (double *theta, double *g, int nc, BFGS_CRIT_FUNC ll, 
			   void *data)
{
    duration_info *dinfo = (duration_info *) data;
    const double *y = dinfo->y->val;
    const double *Xb = dinfo->Xb->val;
    double s = theta[dinfo->k];
    double zt, wt, gti;
    int i, t, err = 0;

    if (dinfo->flags == SCORE_UPDATE_XB) {
	duration_update_Xb(dinfo, theta);
    }

    for (i=0; i<nc; i++) {
	g[i] = 0.0;
    }   

    /* FIXME censoring */

    for (t=0; t<dinfo->T; t++) {
	zt = (y[t] - Xb[t]) / s;
	wt = -(1 - exp(zt)) / s;
	for (i=0; i<nc; i++) {
	    if (i < dinfo->k) {
		gti = wt * gretl_matrix_get(dinfo->X, t, i);
	    } else {
		gti = -zt/s - 1/s + exp(zt) * zt/s;
	    }
	    gretl_matrix_set(dinfo->G, t, i, gti);
	    g[i] += gti;
	}
    }

    return err;
}

static gretl_matrix *duration_nhessian (double *theta, void *data, 
					BFGS_CRIT_FUNC ll, int *err)
{
    duration_info *dinfo = (duration_info *) data;
    gretl_matrix *H = NULL;
    double *g, *splus, *sminus;
    int nc = dinfo->k + 1;
    double x, eps = 1.0e-05;
    int i, j;
    
    splus  = malloc(nc * sizeof *splus);
    sminus = malloc(nc * sizeof *sminus);
    g      = malloc(nc * sizeof *g);
    if (splus == NULL || sminus == NULL || g == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    H = gretl_matrix_alloc(nc, nc);
    if (H == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    dinfo->flags = SCORE_UPDATE_XB;

    for (i=0; i<nc; i++) {
	double theta0 = theta[i];

	theta[i] = theta0 + eps;
	duration_score(theta, g, nc, ll, dinfo);
	for (j=0; j<nc; j++) {
	    splus[j] = g[j];
	}

	theta[i] = theta0 - eps;
	duration_score(theta, g, nc, ll, dinfo);
	for (j=0; j<nc; j++) {
	    sminus[j] = g[j];
	}

	theta[i] = theta0;
	for (j=0; j<nc; j++) {
	    x = -(splus[j] - sminus[j]) / (2*eps);
	    gretl_matrix_set(H, i, j, x);
	}
    }

    dinfo->flags = 0;

    gretl_matrix_xtr_symmetric(H);
    gretl_invert_symmetric_matrix(H);

 bailout:

    free(splus);
    free(sminus);
    free(g);

    return H;
}

/* OPG vcv matrix */

static int duration_OPG_vcv (MODEL *pmod, duration_info *dinfo)
{
    gretl_matrix *GG = NULL;
    int nc = dinfo->k + 1;
    int err = 0;

    GG = gretl_matrix_alloc(nc, nc);
    if (GG == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE,
				  dinfo->G, GRETL_MOD_NONE,
				  dinfo->V, GRETL_MOD_NONE);
	/* invert */
	err = gretl_invert_symmetric_matrix(dinfo->V);
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, dinfo->V);
	if (!err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_OP);
	    pmod->opt |= OPT_G;
	}
    }
    
    gretl_matrix_free(GG);

    return err;
}

/* QML sandwich VCV */

static int duration_robust_vcv (MODEL *pmod, duration_info *dinfo)
{
    gretl_matrix *H = NULL;
    gretl_matrix *GG = NULL;
    int nc = dinfo->k + 1;
    int err = 0;

    GG = gretl_matrix_alloc(nc, nc);
    if (GG == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE,
			      dinfo->G, GRETL_MOD_NONE,
			      GG, GRETL_MOD_NONE);

    H = duration_nhessian(dinfo->theta, dinfo, duration_loglik, &err);

    if (!err) {
	/* form sandwich */
	err = gretl_matrix_qform(H, GRETL_MOD_NONE,
				 GG, dinfo->V, GRETL_MOD_NONE);
    }	

    if (!err) {
	err = gretl_model_write_vcv(pmod, dinfo->V);
	if (!err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
	    pmod->opt |= OPT_R;
	}
    }
    
    gretl_matrix_free(H);
    gretl_matrix_free(GG);

    return err;
}

static int duration_hessian_vcv (MODEL *pmod, duration_info *dinfo)
{
    gretl_matrix *H;
    int err = 0;

    H = duration_nhessian(dinfo->theta, dinfo, duration_loglik, &err);

    if (!err) {
	err = gretl_model_write_vcv(pmod, H);
	if (!err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, VCV_HESSIAN);
	}
    }
    
    gretl_matrix_free(H);

    return err;
}

static int 
transcribe_duration_results (MODEL *pmod, duration_info *dinfo, 
			     const double **Z, const DATAINFO *pdinfo, 
			     int fncount, int grcount,
			     gretlopt opt)
{
    int nc = dinfo->k + 1;
    int i, s, t, err = 0;

    pmod->ci = NEGBIN; /* FIXME */

    gretl_model_set_int(pmod, "fncount", fncount);
    gretl_model_set_int(pmod, "grcount", grcount);

    pmod->ess = 0.0;

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    pmod->yhat[t] = dinfo->Xb->val[s]; /* FIXME */
	    pmod->uhat[t] = dinfo->y->val[s] - pmod->yhat[t];
	    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	    s++;
	}
    }

    if (!err) {
	err = gretl_model_allocate_params(pmod, nc);
	if (!err) {
	    int v;

	    for (i=0; i<dinfo->k; i++) {
		v = pmod->list[i+2];
		strcpy(pmod->params[i], pdinfo->varname[v]);
	    }
	    strcpy(pmod->params[nc-1], "lambda");
	}
    }

    pmod->dfd -= 1;
    pmod->dfn += 1;

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    err = gretl_model_write_coeffs(pmod, dinfo->theta, nc);
    
    if (!err) {
	if (opt & OPT_R) {
	    err = duration_robust_vcv(pmod, dinfo);
	} else if (opt & OPT_G) {
	    err = duration_OPG_vcv(pmod, dinfo);
	} else {
	    err = duration_hessian_vcv(pmod, dinfo);
	} 
    }

    if (!err) {
	pmod->lnL = dinfo->ll;
	mle_criteria(pmod, 0); 
	/* mask invalid statistics (FIXME) */
	pmod->fstt = pmod->chisq = NADBL;
	pmod->rsq = pmod->adjrsq = NADBL;
	pmod->ess = NADBL;
	pmod->sigma = NADBL;
    }

    return err;
}

int duration_estimate (MODEL *pmod, int censvar, const double **Z, 
		       const DATAINFO *pdinfo, gretlopt opt, 
		       PRN *prn)
{
    gretlopt max_opt = (opt & OPT_V)? OPT_V : OPT_NONE;
    duration_info dinfo;
    double toler;
    int maxit;
    int fncount = 0;
    int grcount = 0;
    int err = 0;

    err = duration_init(&dinfo, pmod, censvar, Z, opt, prn);

    if (!err) {
	BFGS_defaults(&maxit, &toler, NEGBIN); /* FIXME */
	err = BFGS_max(dinfo.theta, dinfo.k + 1, maxit, toler, 
		       &fncount, &grcount, duration_loglik, C_LOGLIK,
		       duration_score, &dinfo, NULL, max_opt, 
		       dinfo.prn);
    }

    if (!err) {
	err = transcribe_duration_results(pmod, &dinfo, Z, pdinfo, 
					  fncount, grcount, 
					  opt);
    }

    duration_free(&dinfo);

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }
    
    return pmod->errcode;
}


