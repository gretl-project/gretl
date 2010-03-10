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

#define DDEBUG 0

typedef struct duration_info_ duration_info;

enum {
    DUR_WEIBULL,
    DUR_EXPON,
    DUR_LOGLOG,
    DUR_LOGNORM
};

enum {
    DUR_UPDATE_XB  = 1 << 0,
    DUR_CONST_ONLY = 1 << 1
};

struct duration_info_ {
    int dist;              /* distribution type */
    int flags;             /* control info */
    int k;                 /* number of covariates (including constant) */
    int npar;              /* total number of parameters */
    int n;                 /* number of observations */
    double ll;             /* loglikelihood */    
    double *theta;         /* parameter array, length npar */
    gretl_matrix_block *B; /* workspace */
    gretl_vector *logt;    /* log of dependent variable (duration) */
    gretl_matrix *X;       /* covariates */
    gretl_vector *cens;    /* censoring variable (if needed) */
    gretl_matrix *beta;    /* coeffs on covariates */
    gretl_matrix *llt;     /* per-observation likelihood */
    gretl_matrix *Xb;      /* X \times \beta */
    gretl_matrix *G;       /* score */
    gretl_matrix *V;       /* covariance matrix */
    PRN *prn;              /* verbose printer */
};

static void duration_free (duration_info *dinfo)
{
    gretl_matrix_block_destroy(dinfo->B);
    free(dinfo->theta);
    gretl_matrix_free(dinfo->V);
}

/* initialize using OLS regression of the log of duration
   on the covariates (or a simpler variant if we're 
   estimating the constant-only model)
*/

static int duration_estimates_init (duration_info *dinfo)
{
    int err = 0;

    if (dinfo->flags & DUR_CONST_ONLY) {
	dinfo->theta[0] = gretl_vector_mean(dinfo->logt);
    } else {
	gretl_matrix *b = gretl_matrix_alloc(dinfo->k, 1);
	int j;

	if (b == NULL) {
	    return E_ALLOC;
	}

	err = gretl_matrix_ols(dinfo->logt, dinfo->X, b, 
			       NULL, NULL, NULL);

	if (!err) {
	    for (j=0; j<dinfo->k; j++) {
		dinfo->theta[j] = b->val[j];
	    }
	}

	gretl_matrix_free(b);
    }

    if (dinfo->dist != DUR_EXPON) {
	dinfo->theta[dinfo->k] = 1.0;
    }

    return err;
}

static int duration_init (duration_info *dinfo, MODEL *pmod,
			  int censvar, const double **Z, 
			  gretlopt opt, PRN *prn)
{
    int cn, n = pmod->nobs;
    int np, k = pmod->ncoeff;
    int i, j, t, v;
    int err = 0;

    dinfo->B = NULL;
    dinfo->theta = NULL;
    dinfo->V = NULL;

    if (opt & OPT_E) {
	/* exponential */
	dinfo->dist = DUR_EXPON;
	np = dinfo->npar = k;
    } else if (opt & OPT_L) {
	/* log-logistic */
	dinfo->dist = DUR_LOGLOG;
	np = dinfo->npar = k + 1;
    } else if (opt & OPT_Z) {
	/* log-normal */
	dinfo->dist = DUR_LOGNORM;
	np = dinfo->npar = k + 1;
    } else {
	/* default: Weibull */
	dinfo->dist = DUR_WEIBULL;
	np = dinfo->npar = k + 1;
    }

    dinfo->flags = 0;

    dinfo->theta = malloc(np * sizeof *dinfo->theta);
    dinfo->V = gretl_matrix_alloc(np, np);

    if (dinfo->theta == NULL || dinfo->V == NULL) {
	return E_ALLOC;
    }

    cn = (censvar > 0)? n : 0;

    dinfo->B = gretl_matrix_block_new(&dinfo->logt, n, 1,
				      &dinfo->X, n, k,
				      &dinfo->cens, cn, 1,
				      &dinfo->beta, k, 1,
				      &dinfo->Xb, n, 1,
				      &dinfo->llt, n, 1,
				      &dinfo->G, n, np,
				      NULL);
    if (dinfo->B == NULL) {
	return E_ALLOC;
    }

    if (cn == 0) {
	/* mask unused zero-size part of matrix block */
	dinfo->cens = NULL;
    }

    /* transcribe data into matrix form, taking the
       log of the duration measurements */

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	v = pmod->list[1];
	dinfo->logt->val[i] = log(Z[v][t]);
	if (dinfo->cens != NULL) {
	    dinfo->cens->val[i] = Z[censvar][t];
	}
	for (j=0; j<k; j++) {
	    v = pmod->list[j+2];
	    gretl_matrix_set(dinfo->X, i, j, Z[v][t]);
	}
	i++;
    }

    dinfo->k = k;
    dinfo->n = n;

    err = duration_estimates_init(dinfo);

    if (!err) {
	dinfo->ll = NADBL;
	dinfo->prn = (opt & OPT_V)? prn : NULL;
    }

    return err;
}

static void duration_update_Xb (duration_info *dinfo, const double *theta)
{
    int j;

    if (theta == NULL) {
	theta = dinfo->theta;
    }

    for (j=0; j<dinfo->k; j++) {
	dinfo->beta->val[j] = theta[j];
    }

    gretl_matrix_multiply(dinfo->X, dinfo->beta, dinfo->Xb);
}

#define uncensored(d,i) (d->cens == NULL || d->cens->val[i] == 0)

/* The approach taken here to the loglikelihood and score for duration
   models is that of Kalbfleisch and Prentice: see their Statistical
   Analysis of Failure Time Data, 2e (Wiley, 2002), pp. 68-70.
*/

static double duration_loglik (const double *theta, void *data)
{
    duration_info *dinfo = (duration_info *) data;
    double *ll = dinfo->llt->val;
    double *Xb = dinfo->Xb->val;
    double *logt = dinfo->logt->val;
    double wi, s = 1.0, lns = 0.0;
    double l1ew = 0.0;
    int i, di;

    if (dinfo->dist != DUR_EXPON) {
	s = theta[dinfo->k];
	if (s <= 0) {
	    return NADBL;
	}
	lns = log(s);
    } 

    duration_update_Xb(dinfo, theta);

    dinfo->ll = 0.0;
    errno = 0;

    for (i=0; i<dinfo->n; i++) {
	di = uncensored(dinfo, i);
	wi = (logt[i] - Xb[i]) / s;
	if (dinfo->dist == DUR_LOGLOG) {
	    l1ew = log(1 + exp(wi));
	    ll[i] = -l1ew;
	    if (di) {
		ll[i] += wi - l1ew - lns;
	    }
	} else if (dinfo->dist == DUR_LOGNORM) {
	    if (di) {
		/* density */
		ll[i] = -lns + log_normal_pdf(wi);
	    } else {
		/* survivor */
		ll[i] = log(normal_cdf(-wi));
	    }
	} else {
	    /* Weibull, exponential */
	    ll[i] = -exp(wi);
	    if (di) {
		ll[i] += wi - lns;
	    }
	}	
	dinfo->ll += ll[i];
    }

    if (errno) {
	dinfo->ll = NADBL;
    }

    return dinfo->ll;
}

/* normal hazard: ratio of density to survivor function */

static double normal_h (double w)
{
    return normal_pdf(w) / normal_cdf(-w);
}

static int duration_score (double *theta, double *g, int np, 
			   BFGS_CRIT_FUNC ll, void *data)
{
    duration_info *dinfo = (duration_info *) data;
    const double *logt = dinfo->logt->val;
    const double *Xb = dinfo->Xb->val;
    double wi, ewi, ai, xij, gij, s = 1.0;
    int i, j, di, err = 0;

    if (dinfo->flags == DUR_UPDATE_XB) {
	duration_update_Xb(dinfo, theta);
    }

    if (dinfo->dist != DUR_EXPON) {
	s = theta[dinfo->k];
    } 

    if (g != NULL) {
	for (j=0; j<np; j++) {
	    g[j] = 0.0;
	}
    }  

    for (i=0; i<dinfo->n; i++) {
	di = uncensored(dinfo, i);
	wi = (logt[i] - Xb[i]) / s;
	ewi = exp(wi);
	if (dinfo->dist == DUR_LOGLOG) {
	    ai = -di + (1 + di) * ewi / (1 + ewi);
	} else if (dinfo->dist == DUR_LOGNORM) {
	    ai = di ? wi : normal_h(wi);
	} else {
	    ai = ewi - di;
	}
	for (j=0; j<np; j++) {
	    if (j < dinfo->k) {
		/* covariates */
		xij = gretl_matrix_get(dinfo->X, i, j);
		gij = xij * ai;
	    } else {
		/* scale */
		gij = wi * ai - di;
	    }
	    gij /= s;
	    gretl_matrix_set(dinfo->G, i, j, gij);
	    if (g != NULL) {
		g[j] += gij;
	    }
	}
    }

    return err;
}

#define matrix_plus(m,i,j,x) (m->val[(j)*m->rows+(i)]+=x)

/* Analytical Hessian: see Kalbfleisch and Prentice, 2002, pp. 69-70.
   We're actually constructing the negative inverse of the Hessian here.
*/

static int duration_hessian (double *theta, 
			     gretl_matrix *H,
			     void *data)
{
    duration_info *dinfo = (duration_info *) data;
    const double *logt = dinfo->logt->val;
    const double *Xb = dinfo->Xb->val;
    double s, s2, Ai, wi, ewi, hwi;
    int np = dinfo->npar;
    double xij, xik, hjk;
    int i, j, k, di;
    int err = 0;

    gretl_matrix_zero(H);

    if (dinfo->dist == DUR_EXPON) {
	s2 = s = 1;
    } else {
	s = dinfo->theta[np - 1];
	s2 = s * s;
    }

    for (i=0; i<dinfo->n; i++) {
	di = uncensored(dinfo, i);
	wi = (logt[i] - Xb[i]) / s;
	ewi = exp(wi);

	if (dinfo->dist == DUR_LOGLOG) {
	    Ai = (1 + di) * ewi / ((1 + ewi) * (1 + ewi));
	} else if (dinfo->dist == DUR_LOGNORM) {
	    if (di) {
		Ai = 1;
	    } else {
		hwi = normal_h(wi);
		Ai = hwi * (hwi - wi);
	    }
	} else {
	    /* Weibull */
	    Ai = ewi;
	}

	for (j=0; j<np; j++) {
	    if (j < dinfo->k) {
		/* covariate coeffs cross-block */
		xij = gretl_matrix_get(dinfo->X, i, j);
		for (k=0; k<=j; k++) {
		    xik = gretl_matrix_get(dinfo->X, i, k);
		    hjk = xij * xik * Ai / s2;
		    matrix_plus(H, j, k, hjk);
		}
		if (dinfo->dist != DUR_EXPON) {
		    /* coeff j and scale */
		    hjk = xij * wi * Ai / s2;
		    hjk += gretl_matrix_get(dinfo->G, i, j) / s;
		    matrix_plus(H, np - 1, j, hjk);
		}
	    } else {
		/* scale */
		hjk = (wi * wi * Ai + di) / s2;
		hjk += (2/s) * gretl_matrix_get(dinfo->G, i, j) / s;
		matrix_plus(H, j, j, hjk);
	    }
	}
    }

    /* fill out upper triangle and invert */
    gretl_matrix_mirror(H, 'L');
    err = gretl_invert_symmetric_matrix(H);

    return err;
}

/* calculate the OPG matrix at the starting point and use 
   its inverse (if any) as initial curvature matrix for BFGS
*/

static gretl_matrix *duration_init_H (duration_info *dinfo)
{
    gretl_matrix *H = NULL;
    int err;

    dinfo->flags = DUR_UPDATE_XB;
    err = duration_score(dinfo->theta, NULL, dinfo->npar, 
			 NULL, dinfo);
    dinfo->flags = 0;

    if (!err) {
	H = gretl_matrix_GG_inverse(dinfo->G, &err);
    }

    return H;
}

/* OPG vcv matrix */

static int duration_OPG_vcv (duration_info *dinfo)
{
    gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE,
			      dinfo->G, GRETL_MOD_NONE,
			      dinfo->V, GRETL_MOD_NONE);

    return gretl_invert_symmetric_matrix(dinfo->V);
}

/* QML sandwich VCV */

static int duration_robust_vcv (duration_info *dinfo)
{
    gretl_matrix *H = NULL;
    gretl_matrix *GG = NULL;
    int err = 0;

    GG = gretl_matrix_alloc(dinfo->npar, dinfo->npar);
    if (GG == NULL) {
	return E_ALLOC;
    }

    H = gretl_matrix_alloc(dinfo->npar, dinfo->npar);
    if (H == NULL) {
	gretl_matrix_free(GG);
	return E_ALLOC;
    }   

    gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE,
			      dinfo->G, GRETL_MOD_NONE,
			      GG, GRETL_MOD_NONE);

    err = duration_hessian(NULL, H, dinfo);

    if (!err) {
	err = gretl_matrix_qform(H, GRETL_MOD_NONE,
				 GG, dinfo->V, GRETL_MOD_NONE);
    }	

    gretl_matrix_free(GG);
    gretl_matrix_free(H);

    return err;
}

/* This is the last thing we do, after transcribing the 
   MLE results to pmod, so we don't have to worry about 
   saving and then restoring all the original values
   attached to dinfo, which we will shortly destroy.
*/

static void
duration_overall_LR_test (MODEL *pmod, duration_info *dinfo)
{
    double llu = dinfo->ll;
    int err = 0;

    dinfo->k = 1;
    dinfo->npar = 1 + (dinfo->dist != DUR_EXPON);

    gretl_matrix_reuse(dinfo->X, -1, dinfo->k);
    gretl_matrix_reuse(dinfo->G, -1, dinfo->npar);
    gretl_matrix_reuse(dinfo->beta, dinfo->k, 1);

    dinfo->flags |= DUR_CONST_ONLY;
    err = duration_estimates_init(dinfo);

    if (!err) {
	int maxit, fncount = 0, grcount = 0;
	double toler;

	/* estimate constant-only model */

	BFGS_defaults(&maxit, &toler, DURATION); 
	err = BFGS_max(dinfo->theta, dinfo->npar, maxit, toler, 
		       &fncount, &grcount, duration_loglik, C_LOGLIK,
		       duration_score, dinfo, NULL, OPT_NONE, NULL);
    }

    if (!err && llu > dinfo->ll) {
	pmod->chisq = 2 * (llu - dinfo->ll);
    }
}

static void duration_set_predictions (MODEL *pmod, duration_info *dinfo,
				      const double **Z)
{
    const double *y = Z[pmod->list[1]];
    const double *logt = dinfo->logt->val;
    double St, G = 1.0;
    double s = 1.0, p = 1.0;
    double wi, Xbi, expXbi;
    int i, t;

    if (dinfo->dist != DUR_EXPON) {
	/* scale factor */
	s = dinfo->theta[dinfo->npar-1];
	p = 1/s;
    }

    if (dinfo->dist == DUR_WEIBULL) {
	/* agrees with Stata; R's "survreg" has this wrong? */
	G = gamma_function(1 + s);
    } else if (dinfo->dist == DUR_EXPON) {
	G = gamma_function(2.0);
    }   

    /* Below: we write into pmod->yhat E[t | X, theta] -- or
       in the case of the loglogistic and lognormal models,
       exp(E[log t | X, theta]).  And we write into pmod->uhat
       Cox-Snell generalized residuals, namely the integrated
       hazard function, which equals -log S(t).
    */

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->yhat[t])) {
	    continue;
	}

	Xbi = dinfo->Xb->val[i];
	expXbi = exp(Xbi);
	wi = (logt[i] - Xbi) / s;

	if (dinfo->dist == DUR_WEIBULL || dinfo->dist == DUR_EXPON) {
	    pmod->yhat[t] = expXbi * G;
	    St = exp(-exp(wi));
	} else if (dinfo->dist == DUR_LOGNORM) {
	    pmod->yhat[t] = expXbi;
	    St = normal_cdf(-wi);
	} else {
	    /* log-logistic */
	    pmod->yhat[t] = expXbi;
	    St = 1.0 / (1 + pow(y[t] / expXbi, p));
	}

	/* generalized (Cox-Snell) residual */
	pmod->uhat[t] = -log(St);

	i++;
    }
}

static int 
transcribe_duration_results (MODEL *pmod, duration_info *dinfo, 
			     const double **Z, const DATAINFO *pdinfo, 
			     int fncount, int grcount,
			     int censvar, gretlopt opt)
{
    int np = dinfo->npar;
    int j, v, err = 0;

    pmod->ci = DURATION;

    if (dinfo->dist == DUR_EXPON) {
	pmod->opt |= OPT_E;
    } else if (dinfo->dist == DUR_LOGLOG) {
	pmod->opt |= OPT_L;
    } else if (dinfo->dist == DUR_LOGNORM) {
	pmod->opt |= OPT_Z;
    }

    if (censvar > 0) {
	gretl_model_set_int(pmod, "cens_var", censvar);
    }

    gretl_model_set_int(pmod, "fncount", fncount);
    gretl_model_set_int(pmod, "grcount", grcount);

    if (!err) {
	err = gretl_model_allocate_params(pmod, np);
	if (!err) {
	    for (j=0; j<dinfo->k; j++) {
		v = pmod->list[j+2];
		strcpy(pmod->params[j], pdinfo->varname[v]);
	    }
	    if (dinfo->dist != DUR_EXPON) {
		strcpy(pmod->params[np-1], "sigma");
	    } 
	}
    }

    err = gretl_model_write_coeffs(pmod, dinfo->theta, np);
    
    if (!err) {
	if (opt & OPT_R) {
	    pmod->opt |= OPT_R;
	    err = duration_robust_vcv(dinfo);
	} else if (opt & OPT_G) {
	    pmod->opt |= OPT_G;
	    err = duration_OPG_vcv(dinfo);
	} else {
	    err = duration_hessian(NULL, dinfo->V, dinfo);
	} 
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, dinfo->V);
	if (!err) {
	    int vtype = (opt & OPT_G) ? VCV_OP :
		(opt & OPT_R)? VCV_QML : VCV_HESSIAN;
	    
	    gretl_model_set_vcv_info(pmod, VCV_ML, vtype);
	}
    }	

    if (!err) {
	duration_set_predictions(pmod, dinfo, Z);
	pmod->lnL = dinfo->ll;
	mle_criteria(pmod, 0); 
	/* mask invalid statistics */
	pmod->fstt = pmod->chisq = NADBL;
	pmod->rsq = pmod->adjrsq = NADBL;
	pmod->ess = NADBL;
	pmod->sigma = NADBL;
	/* but add overall LR test if possible */
	duration_overall_LR_test(pmod, dinfo);
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

#if DDEBUG
    if (!err && censvar > 0) {
	fprintf(stderr, "duration: using var %d (%s) for censoring info\n",
		censvar, pdinfo->varname[censvar]);
    }
#endif

    if (!err) {
	/* initialize BFGS curvature */
	gretl_matrix *H = duration_init_H(&dinfo);

	BFGS_defaults(&maxit, &toler, DURATION); 
	err = BFGS_max(dinfo.theta, dinfo.npar, maxit, toler, 
		       &fncount, &grcount, duration_loglik, C_LOGLIK,
		       duration_score, &dinfo, H, max_opt, 
		       dinfo.prn);
	gretl_matrix_free(H);
    }

    if (!err) {
	err = transcribe_duration_results(pmod, &dinfo, Z, pdinfo, 
					  fncount, grcount, 
					  censvar, opt);
    }

    duration_free(&dinfo);

    if (err && !pmod->errcode) {
	pmod->errcode = err;
    }
    
    return pmod->errcode;
}
