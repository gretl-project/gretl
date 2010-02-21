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
    DUR_EXPON
};

enum {
    DOING_HESSIAN = 1
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
    gretl_matrix *beta;    /* coeffs on covariates */
    gretl_matrix *llt;     /* per-observation likelihood */
    gretl_matrix *Xb;      /* X \times \beta */
    gretl_matrix *G;       /* score */
    gretl_matrix *V;       /* covariance matrix */
    gretl_vector *cens;    /* censoring variable (if needed) */
    PRN *prn;              /* verbose printer */
};

static void duration_free (duration_info *dinfo)
{
    gretl_matrix_block_destroy(dinfo->B);
    free(dinfo->theta);
    gretl_matrix_free(dinfo->V);
    gretl_matrix_free(dinfo->cens);
}

static int duration_nonpositive (const MODEL *pmod, const double **Z)
{
    const double *y = Z[pmod->list[1]];
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t]) && y[t] <= 0) {
	    gretl_errmsg_set("durations must be positive");
	    return 1;
	}
    }

    return 0;
}

/* initialize using OLS regression of the log of duration
   on the covariates */

static int duration_estimates_init (duration_info *dinfo)
{
    gretl_matrix *b;
    int j, err = 0;

    b = gretl_matrix_alloc(dinfo->k, 1);
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

    return err;
}

static int duration_init (duration_info *dinfo, MODEL *pmod,
			  int censvar, const double **Z, 
			  gretlopt opt, PRN *prn)
{
    int n = pmod->nobs;
    int np, k = pmod->ncoeff;
    int i, j, t, v;
    int err = 0;

    dinfo->B = NULL;
    dinfo->theta = NULL;
    dinfo->cens = NULL;
    dinfo->V = NULL;

    if (duration_nonpositive(pmod, Z)) {
	return E_DATA;
    }

    if (opt & OPT_E) {
	dinfo->dist = DUR_EXPON;
	np = dinfo->npar = k;
    } else {
	dinfo->dist = DUR_WEIBULL;
	np = dinfo->npar = k + 1;
    }

    dinfo->flags = 0;

    dinfo->theta = malloc(np * sizeof *dinfo->theta);
    dinfo->V = gretl_matrix_alloc(np, np);

    if (dinfo->theta == NULL || dinfo->V == NULL) {
	return E_ALLOC;
    }

    if (censvar > 0) {
	dinfo->cens = gretl_column_vector_alloc(n);
	if (dinfo->cens == NULL) {
	    return E_ALLOC;
	}
    }

    dinfo->B = gretl_matrix_block_new(&dinfo->logt, n, 1,
				      &dinfo->X, n, k,
				      &dinfo->beta, k, 1,
				      &dinfo->Xb, n, 1,
				      &dinfo->llt, n, 1,
				      &dinfo->G, n, np,
				      NULL);
    if (dinfo->B == NULL) {
	return E_ALLOC;
    }

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
	if (dinfo->dist == DUR_WEIBULL) {
	    dinfo->theta[k] = 1.0;
	}
	dinfo->ll = NADBL;
	dinfo->prn = (opt & OPT_V)? prn : NULL;
    }

    return err;
}

static void duration_update_Xb (duration_info *dinfo, const double *theta)
{
    int j;

    for (j=0; j<dinfo->k; j++) {
	dinfo->beta->val[j] = theta[j];
    }

    gretl_matrix_multiply(dinfo->X, dinfo->beta, dinfo->Xb);
}

#define uncensored(d,i) (d->cens == NULL || d->cens->val[i] == 0)

static double duration_loglik (const double *theta, void *data)
{
    duration_info *dinfo = (duration_info *) data;
    double *ll = dinfo->llt->val;
    double *Xb = dinfo->Xb->val;
    double *logt = dinfo->logt->val;
    double wi, s = 1.0, lns = 0.0;
    int i;

    if (dinfo->dist == DUR_WEIBULL) {
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
	wi = (logt[i] - Xb[i]) / s;
	ll[i] = -exp(wi);
	if (uncensored(dinfo, i)) {
	    ll[i] += wi - lns;
	} 
	dinfo->ll += ll[i];
    }

    if (errno) {
	dinfo->ll = NADBL;
    }

    return dinfo->ll;
}

static int duration_score (double *theta, double *g, int np, 
			   BFGS_CRIT_FUNC ll, void *data)
{
    duration_info *dinfo = (duration_info *) data;
    const double *logt = dinfo->logt->val;
    const double *Xb = dinfo->Xb->val;
    double wi, xij, gij, s = 1.0;
    int i, j, err = 0;

    if (dinfo->flags == DOING_HESSIAN) {
	duration_update_Xb(dinfo, theta);
    }

    if (dinfo->dist == DUR_WEIBULL) {
	s = theta[dinfo->k];
    } 

    if (g != NULL) {
	for (j=0; j<np; j++) {
	    g[j] = 0.0;
	}
    }  

    for (i=0; i<dinfo->n; i++) {
	wi = (logt[i] - Xb[i]) / s;
	for (j=0; j<np; j++) {
	    if (j < dinfo->k) {
		/* survival */
		xij = gretl_matrix_get(dinfo->X, i, j);
		gij = -exp(wi) * (-xij/s);
		if (uncensored(dinfo, i)) {
		    /* hazard */
		    gij += -xij/s;
		}
	    } else {
		/* survival */
		gij = -exp(wi) * (-wi/s);
		if (uncensored(dinfo, i)) {
		    /* hazard */
		    gij += -wi/s - 1/s;
		}
	    }
	    gretl_matrix_set(dinfo->G, i, j, gij);
	    if (g != NULL) {
		g[j] += gij;
	    }
	}
    }

    return err;
}

/* calculate the OPG matrix at the starting point and use 
   its inverse (if any) as initial curvature matrix for BFGS
*/

static gretl_matrix *duration_init_H (duration_info *dinfo)
{
    gretl_matrix *H = NULL;
    int np = dinfo->npar;
    int err;

    dinfo->flags = DOING_HESSIAN;
    err = duration_score(dinfo->theta, NULL, np, NULL, dinfo);
    dinfo->flags = 0;

    if (err) {
	return NULL;
    }

    H = gretl_matrix_alloc(np, np);
    if (H == NULL) {
	return NULL;
    }

    gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE, 
			      dinfo->G, GRETL_MOD_NONE, 
			      H, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(H); 
    if (err) {
	fprintf(stderr, "duration: init_H not pd\n");
	gretl_matrix_free(H);
	H = NULL;
    }

    return H;
}

/* Given the change of variables that we use for easy computation of
   the Weibull/exponential loglikelihood and score, we need to use the
   appropriate Jacobian to get the covariance matrix right.  
*/

static int duration_vcv_transform (duration_info *dinfo)
{
    gretl_matrix *J = NULL;
    gretl_matrix *JVJ = NULL;
    const double *theta = dinfo->theta;
    double a, s2, s = theta[dinfo->k];
    int np = dinfo->npar;
    int i, k = np - 1;
    int err = 0;

    J = gretl_zero_matrix_new(np, np);
    JVJ = gretl_matrix_alloc(np, np);

    if (J == NULL || JVJ == NULL) {
	return E_ALLOC;
    }

    a = 1/s;
    s2 = s * s;

    for (i=0; i<np; i++) {
	if (i < np - 1) {
	    gretl_matrix_set(J, i, i, a);
	    gretl_matrix_set(J, i, k, -theta[i] / s2);
	} else {
	    gretl_matrix_set(J, i, k, -1 / s2);
	}
    }

    gretl_matrix_qform(J, GRETL_MOD_NONE, dinfo->V, 
		       JVJ, GRETL_MOD_NONE);

    if (!err) {
	gretl_matrix_replace(&dinfo->V, JVJ);
	JVJ = NULL;
    }

    gretl_matrix_free(J);
    gretl_matrix_free(JVJ);

    return 0;
}

/* numerical Hessian using the analytical score */

static gretl_matrix *duration_nhessian (double *theta, void *data, 
					BFGS_CRIT_FUNC ll, int *err)
{
    duration_info *dinfo = (duration_info *) data;
    gretl_matrix *H = NULL;
    double *g, *splus, *sminus;
    double x, eps = 1.0e-05;
    int np = dinfo->npar;
    int i, j;
    
    splus  = malloc(np * sizeof *splus);
    sminus = malloc(np * sizeof *sminus);
    g      = malloc(np * sizeof *g);

    if (splus == NULL || sminus == NULL || g == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    H = gretl_matrix_alloc(np, np);
    if (H == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    dinfo->flags = DOING_HESSIAN;

    for (i=0; i<np && !*err; i++) {
	double theta0 = theta[i];

	theta[i] = theta0 + eps;
	*err += duration_score(theta, g, np, ll, dinfo);
	for (j=0; j<np; j++) {
	    splus[j] = g[j];
	}

	theta[i] = theta0 - eps;
	*err += duration_score(theta, g, np, ll, dinfo);
	for (j=0; j<np; j++) {
	    sminus[j] = g[j];
	}

	theta[i] = theta0;
	for (j=0; j<np; j++) {
	    x = -(splus[j] - sminus[j]) / (2*eps);
	    gretl_matrix_set(H, i, j, x);
	}
    }

    dinfo->flags = 0;

    if (*err) {
	*err = E_NAN;
    } else {
	gretl_matrix_xtr_symmetric(H);
	*err = gretl_invert_symmetric_matrix(H);
    }

 bailout:

    if (*err) {
	gretl_matrix_free(H);
	H = NULL;
    }

    free(splus);
    free(sminus);
    free(g);

    return H;
}

/* OPG vcv matrix */

static int duration_OPG_vcv (MODEL *pmod, duration_info *dinfo)
{
    gretl_matrix *GG = NULL;
    int err = 0;

    GG = gretl_matrix_alloc(dinfo->npar, dinfo->npar);
    if (GG == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE,
				  dinfo->G, GRETL_MOD_NONE,
				  dinfo->V, GRETL_MOD_NONE);
	err = gretl_invert_symmetric_matrix(dinfo->V);
    }

    gretl_matrix_free(GG);

    return err;
}

/* QML sandwich VCV */

static int duration_robust_vcv (MODEL *pmod, duration_info *dinfo)
{
    gretl_matrix *H = NULL;
    gretl_matrix *GG = NULL;
    int err = 0;

    GG = gretl_matrix_alloc(dinfo->npar, dinfo->npar);
    if (GG == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(dinfo->G, GRETL_MOD_TRANSPOSE,
			      dinfo->G, GRETL_MOD_NONE,
			      GG, GRETL_MOD_NONE);

    H = duration_nhessian(dinfo->theta, dinfo, duration_loglik, &err);

    if (!err) {
	err = gretl_matrix_qform(H, GRETL_MOD_NONE,
				 GG, dinfo->V, GRETL_MOD_NONE);
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
	gretl_matrix_replace(&dinfo->V, H);
    }

    return err;
}

static int 
transcribe_duration_results (MODEL *pmod, duration_info *dinfo, 
			     const double **Z, const DATAINFO *pdinfo, 
			     int fncount, int grcount,
			     int censvar, gretlopt opt)
{
    const double *y = Z[pmod->list[1]];
    double Et_mult = 1.0;
    int np = dinfo->npar;
    int i, j, t, err = 0;

    pmod->ci = DURATION;

    if (dinfo->dist == DUR_EXPON) {
	pmod->opt |= OPT_E;
    }

    if (censvar > 0) {
	gretl_model_set_int(pmod, "cens_var", censvar);
    }

    gretl_model_set_int(pmod, "fncount", fncount);
    gretl_model_set_int(pmod, "grcount", grcount);

    if (dinfo->dist == DUR_WEIBULL) {
	Et_mult = gamma_function(dinfo->theta[np-1] + 1);
    } else if (dinfo->dist == DUR_EXPON) {
	Et_mult = gamma_function(2.0);
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->yhat[t])) {
	    if (dinfo->dist == DUR_WEIBULL || dinfo->dist == DUR_WEIBULL) {
		pmod->yhat[t] = exp(dinfo->Xb->val[i]);
		pmod->yhat[t] *= Et_mult;
		/* FIXME generalized residual */
		pmod->uhat[t] = y[t] - pmod->yhat[t];
	    } else {
		/* FIXME other distributions */
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    }
	    i++;
	}
    }

    if (!err) {
	err = gretl_model_allocate_params(pmod, np);
	if (!err) {
	    int v;

	    for (j=0; j<dinfo->k; j++) {
		v = pmod->list[j+2];
		strcpy(pmod->params[j], pdinfo->varname[v]);
	    }
	    if (dinfo->dist == DUR_WEIBULL) {
		strcpy(pmod->params[np-1], "alpha");
	    }
	}
    }

    pmod->dfd -= (dinfo->npar - dinfo->k);
    pmod->dfn += (dinfo->npar - dinfo->k);

    err = gretl_model_write_coeffs(pmod, dinfo->theta, np);
    
    if (!err) {
	if (opt & OPT_R) {
	    pmod->opt |= OPT_R;
	    err = duration_robust_vcv(pmod, dinfo);
	} else if (opt & OPT_G) {
	    pmod->opt |= OPT_G;
	    err = duration_OPG_vcv(pmod, dinfo);
	} else {
	    err = duration_hessian_vcv(pmod, dinfo);
	} 
    }

    if (!err && dinfo->dist == DUR_WEIBULL) {
	err = duration_vcv_transform(dinfo);
	pmod->coeff[dinfo->k] = 1 / dinfo->theta[dinfo->k];
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, dinfo->V);
	if (!err) {
	    int vtype = (opt & OPT_G) ? VCV_OP :
		(opt & OPT_R)? VCV_QML : VCV_HESSIAN;
	    
	    gretl_model_set_vcv_info(pmod, VCV_ML, vtype);
	}
    }	

    for (j=0; j<dinfo->k; j++) {
	pmod->coeff[j] = -pmod->coeff[j];
	if (dinfo->dist == DUR_WEIBULL) {
	    pmod->coeff[j] /= dinfo->theta[dinfo->k];
	}
    }

    if (!err) {
	pmod->lnL = dinfo->ll;
	mle_criteria(pmod, 0); 
	/* mask invalid statistics (FIXME chisq?) */
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
    gretl_matrix *H = NULL;
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
	/* to initialize BFGS curvature */
	H = duration_init_H(&dinfo);
    }

    if (!err) {
	BFGS_defaults(&maxit, &toler, DURATION); 
	err = BFGS_max(dinfo.theta, dinfo.npar, maxit, toler, 
		       &fncount, &grcount, duration_loglik, C_LOGLIK,
		       duration_score, &dinfo, H, max_opt, 
		       dinfo.prn);
    }

    gretl_matrix_free(H);	

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


