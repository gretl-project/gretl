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
#include "gretl_bfgs.h"

#define ODEBUG 0

#define OPROBIT_TOL 1.0e-12

typedef struct op_container_ op_container;

struct op_container_ {
    int type;         /* model type (probit or logit) */
    gretlopt opt;     /* option flags */
    int *y;           /* dependent variable */
    double **Z;       /* data */
    int *list;        /* dependent var plus regular regressors */
    int M;            /* max of (possibly normalized) y */
    int t1;           /* beginning of sample */
    int t2;           /* end of sample */
    int nobs;         /* number of observations */
    int nx;           /* number of explanatory variables */
    int k;            /* total number of parameters */
    double *theta;    /* real parameter estimates */
    double *ndx;      /* index variable */
    double *dP;       /* probabilities */
    MODEL *pmod;      /* model struct, initially containing OLS */
    double **G;       /* score matrix by observation */
    double *g;        /* total score vector */
};

static double distfunc (double x, int type)
{
    switch (type) {
    case PROBIT:
	return normal_cdf(x);
    case LOGIT:
	return 1.0 / (1.0 + exp(-x));
    default:
	return NADBL;
    }
}

static double densfunc (double x, int type)
{
    double tmp;

    switch (type) {
    case PROBIT:
	return normal_pdf(x);
    case LOGIT:
	tmp = 1.0 + exp(-x);
	return (tmp - 1.0) / (tmp * tmp);
    default:
	return NADBL;
    }
}

static void op_container_destroy (op_container *OC)
{
    free(OC->y);
    free(OC->ndx);
    free(OC->dP);
    free(OC->list);
    doubles_array_free(OC->G, OC->k);
    free(OC->g);
    free(OC->theta);

    free(OC);
}

static op_container *op_container_new (int type, int ndum,
				       double **Z, MODEL *pmod,  
				       gretlopt opt)
{
    op_container *OC;
    int i, t, vy = pmod->list[1];
    int nobs = pmod->nobs;

    OC = malloc(sizeof *OC);
    if (OC == NULL) {
	return NULL;
    }

    OC->type = type;

    OC->Z = Z;
    OC->pmod = pmod;
    OC->t1 = pmod->t1;
    OC->t2 = pmod->t2;
    OC->nobs = nobs;
    OC->k = pmod->ncoeff;
    OC->M = ndum;
    OC->nx = OC->k - ndum;

    OC->opt = opt;

    OC->y = NULL;
    OC->ndx = NULL;
    OC->dP = NULL;
    OC->list = NULL;
    OC->G = NULL;
    OC->g = NULL;

    OC->y = malloc(nobs * sizeof *OC->y);
    OC->ndx = malloc(nobs * sizeof *OC->ndx);
    OC->dP = malloc(nobs * sizeof *OC->dP);

    OC->list = gretl_list_new(1 + OC->nx);
    OC->G = doubles_array_new(OC->k, nobs);
    OC->g = malloc(OC->k * sizeof *OC->g);
    OC->theta = malloc(OC->k * sizeof *OC->theta);

    if (OC->y == NULL || OC->ndx == NULL || 
	OC->dP == NULL || OC->list == NULL ||
	OC->G == NULL || OC->g == NULL ||
	OC->theta == NULL) {
	op_container_destroy(OC);
	return NULL;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    OC->y[i++] = (int) Z[vy][t];
	}
    }

    OC->list[1] = vy;
    for (i=0; i<OC->nx; i++) {
	OC->list[i+2] = pmod->list[i+2];
    }

#if ODEBUG
    fprintf(stderr, "nobs = %d\n", OC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", OC->t1, OC->t2);
    fprintf(stderr, "k = %d\n", OC->k);
    fprintf(stderr, "nx = %d\n", OC->nx);
    fprintf(stderr, "Max(y) = M = %d\n", OC->M);
    printlist(OC->list, "list, in op_container_new");
#endif

    return OC;
}

static int compute_score (op_container *OC, int yt, 
			  double ystar0, double ystar1,
			  double dP, int t, int s)
{
    double mills0;
    double mills1;
    double dm;
    int M = OC->M;
    int i, v;

    if (ystar1 < 6.0 || OC->type == LOGIT) {
	mills0 = (yt == 0)? 0.0 : densfunc(ystar0, OC->type) / dP;
	mills1 = (yt == M)? 0.0 : densfunc(ystar1, OC->type) / dP;
    } else { 
	/* L'Hopital-based approximation */
	mills0 = (yt == 0)? 0.0 : -ystar0;
	mills1 = (yt == M)? 0.0 : -ystar1;
    }

    dm = mills1 - mills0;

    for (i=0; i<OC->nx; i++) {
	v = OC->list[i+2];
	OC->G[i][s] = -dm * OC->Z[v][t];
	OC->g[i] += OC->G[i][s];
    }

    for (i=OC->nx; i<OC->k; i++) {
	OC->G[i][s] = 0.0;
	if (i == OC->nx + yt - 1) {
	    OC->G[i][s] = -mills0;
	    OC->g[i] += OC->G[i][s];
	}
	if (i == OC->nx + yt) {
	    OC->G[i][s] = mills1;
	    OC->g[i] += OC->G[i][s];
	}
    }

    return 0;
}

static int compute_probs (const double *theta, op_container *OC)
{
    double m0, m1, ystar0 = 0.0, ystar1 = 0.0;
    int M = OC->M;
    int nx = OC->nx;
    double P0, P1, h, adj, dP;
    int i, t, s, yt;
    int type = OC->type;

    /* initialize analytical score */
    for (i=0; i<OC->k; i++) {
	OC->g[i] = 0.0;
    }

    s = 0;

    for (t=OC->pmod->t1; t<=OC->pmod->t2; t++) {
	if (na(OC->pmod->uhat[t])) {
#if ODEBUG > 1
	    fprintf(stderr, "obs %4d excluded\n", t);
#endif
	    continue;
	}

	yt = OC->y[s];

	if (yt == 0) {
	    m0 = theta[nx];
	    ystar1 = OC->ndx[s] + m0;
	} else {
	    m0 = theta[nx + yt - 1];
	    ystar0 = OC->ndx[s] + m0;
	    if (yt < M) {
		m1 = theta[nx + yt];
		ystar1 = OC->ndx[s] + m1;
	    }
	} 

#if ODEBUG > 1
	fprintf(stderr, "t:%4d/%d s=%d y=%d, ndx = %10.6f, ystar0 = %9.7f, ystar1 = %9.7f\n", 
		t, OC->nobs, s, yt, OC->ndx[s], ystar0, ystar1);
#endif

	if (ystar0 < 6.0 || OC->type == LOGIT) {
	    P0 = (yt == 0)? 0.0 : distfunc(ystar0, type);
	    P1 = (yt == M)? 1.0 : distfunc(ystar1, type);
	    dP = P1 - P0;
	} else { 
	    /* Taylor-based 1st order approximation */
	    h = ystar1 - ystar0;
	    adj = densfunc(ystar1, type) + densfunc(ystar0, type);
	    dP =  0.5 * h * adj;
	}

	if (dP > 1.0e-15) {
	    OC->dP[s] = dP;
	} else {
#if ODEBUG > 1
	    fprintf(stderr, "very small dP at obs %d; y=%d, ndx = %10.6f, dP = %9.7f\n", 
 		    t, yt, OC->ndx[s], dP);
#endif
	    return 1;
	} 

	compute_score(OC, yt, ystar0, ystar1, dP, t, s);

	s++;
    }

    return 0;
}

/* Below: method for getting around the "non-increasing cut point"
   issue by construction: the 2nd and higher cut points are
   represented to BFGS in the form of the log-difference from the
   previous cut point.
   cases.  
*/

static void make_BFGS_theta (op_container *OC, double *theta)
{
    int i;

    for (i=0; i<=OC->nx; i++) {
	theta[i] = OC->theta[i];
    }

    for (i=OC->nx+1; i<OC->k; i++) {
	/* convert cut point 2 and higher to log-difference form */
	theta[i] = log(OC->theta[i] - OC->theta[i-1]);
    }
}

/* Inverse operation for the transformation done by
   make_BFGS_theta() */

static void get_real_theta (op_container *OC, const double *theta)
{
    int i;

    for (i=0; i<=OC->nx; i++) {
	OC->theta[i] = theta[i];
    }

    for (i=OC->nx+1; i<OC->k; i++) {
	/* retrieve cut point 2 and higher from log-difference form */
	OC->theta[i] = exp(theta[i]) + OC->theta[i-1];
    }
}

static double op_loglik (const double *theta, void *ptr)
{
    op_container *OC = (op_container *) ptr;
    double x, ll = 0.0;
    int i, s, t, v;
    int err;

    if (theta != OC->theta) {
	get_real_theta(OC, theta);
    }

    s = 0;
    for (t=OC->t1; t<=OC->t2; t++) {
	if (na(OC->pmod->uhat[t])) {
	    continue;
	}
	x = 0.0;
	for (i=0; i<OC->nx; i++) {
	    /* the independent variables */
	    v = OC->list[i+2];
	    x -= OC->theta[i] * OC->Z[v][t];
	}
	OC->ndx[s++] = x;
#if ODEBUG > 2
	fprintf(stderr, "t = %d, s = %d, x = %g\n", t, s, x);
#endif
    }
    
    err = compute_probs(OC->theta, OC);
    if (err) {
	ll = NADBL;
    } else {
	s = 0;
	for (t=OC->t1; t<=OC->t2; t++) {
	    if (!na(OC->pmod->uhat[t])) {
		ll += log(OC->dP[s++]);
	    }
	}
    }

#if ODEBUG > 1
    fprintf(stderr, "ll = %16.10f\n", ll);
#endif

    return ll;
}

static int op_score (double *theta, double *s, int npar, BFGS_CRIT_FUNC ll, 
		     void *ptr)
{
    op_container *OC = (op_container *) ptr;
    int i, j;

    for (i=0; i<npar; i++) {
	s[i] = OC->g[i];
    }

    for (i=OC->nx; i<npar; i++) {
	for (j=i+1; j<npar; j++) {
	    /* add effects of changes in subsequent cut points */
	    s[i] += OC->g[j];
	}
	if (i > OC->nx) {
	    s[i] *= exp(theta[i]); /* apply jacobian */
	}
    }

    return 1;
}

static int opg_from_ascore (op_container *OC, gretl_matrix *GG) 
{
    int s, t, i, j, k = OC->k;
    double x, ll, *g0;

    g0 = malloc(k * sizeof *g0);
    if (g0 == NULL) {
	return E_ALLOC;
    }

    ll = op_loglik(OC->theta, OC);

    for (j=0; j<k; j++) {
	g0[j] = OC->g[j];
    }

    for (i=0; i<k; i++) {
	for (j=i; j<k; j++) {
	    x = 0.0;
	    s = 0;
	    for (t=OC->t1; t<=OC->t2; t++) {
		if (na(OC->pmod->uhat[t])) {
		    continue;
		}
		x += OC->G[i][s] * OC->G[j][s];
		s++;
	    }
	    gretl_matrix_set(GG, j, i, x);
	    gretl_matrix_set(GG, i, j, x);
	}
    }

    free(g0);

    return 0;
}

static int ihess_from_ascore (op_container *OC, gretl_matrix *inH) 
{
    int i, j, err, k = OC->k;
    double smal = 1.0e-07;  /* "small" is some sort of macro on win32 */
    double smal2 = 2.0 * smal;
    double ti, x, ll, *g0;

    g0 = malloc(k * sizeof *g0);
    if (g0 == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	ti = OC->theta[i];
	OC->theta[i] -= smal;
	ll = op_loglik(OC->theta, OC);
	for (j=0; j<k; j++) {
	    g0[j] = OC->g[j];
	}
	OC->theta[i] += smal2;
	ll = op_loglik(OC->theta, OC);
	for (j=0; j<k; j++) {
	    x = (OC->g[j] - g0[j]) / smal2;
	    gretl_matrix_set(inH, i, j, -x);
	}
	/* restore original theta */
	OC->theta[i] = ti;
    }

    gretl_matrix_xtr_symmetric(inH);

    free(g0);

#if ODEBUG
    gretl_matrix_print(inH, "inverse of Hessian");
#endif

    err = gretl_invert_symmetric_matrix(inH);
    if (err) {
	fprintf(stderr, "ihess_from_ascore: failed to invert numerical Hessian\n");
    }

    return err;
}

/* compare X*b with the cut points to get a prediction for y */

static int get_prediction (op_container *OC, const MODEL *pmod,
			   double Xb)
{
    int k = OC->k - 1; /* position of greatest cut point */
    double cut;
    int i, pred = 0;

    for (i=OC->M; i>0; i--) {
	cut = pmod->coeff[k--];
	if (Xb > cut) {
	    pred = i;
	    break;
	}
    }

    return pred;
}

static double gen_resid (op_container *OC, const double *theta, int t) 
{
    double ndxt, m0, m1, ystar0, f0, f1;
    double ret, dP, ystar1 = 0.0;
    int M = OC->M;
    int nx = OC->nx;
    int yt;

    dP = OC->dP[t];
    yt = OC->y[t];
    ndxt = OC->ndx[t];

    if (yt == 0) {
	m0 = theta[nx];
	ystar1 = ndxt + m0;
    } else {
	m0 = theta[nx + yt - 1];
	ystar0 = ndxt + m0;
	if (yt < M) {
	    m1 = theta[nx + yt];
	    ystar1 = ndxt + m1;
	}
    } 

    if (ystar1 < 6.0 || OC->type == LOGIT || 1) {
	f0 = (yt == 0)? 0.0 : densfunc(ystar0, OC->type) / dP;
	f1 = (yt == M)? 0.0 : densfunc(ystar1, OC->type) / dP;
    } else { 
	/* L'Hopital-based approximation */
	f0 = (yt == 0)? 0.0 : -ystar0;
	f1 = (yt == M)? 0.0 : -ystar1;
    }

    ret = (f0 - f1);

    return ret;
} 

/* Initialize the cut-points by counting the occurrences of each value
   of the (normalized) dependent variable, finding the sample
   proportion (cumulating as we go), and taking the inverse of the
   normal CDF.
*/

static void cut_points_init (op_container *OC, const MODEL *pmod, 
			     const double **Z)
{
    const double *y = Z[pmod->list[1]];
    double p = 0.0;
    int i, j, t, nj;

    for (i=OC->nx, j=0; i<OC->k; i++, j++) {
	nj = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t]) && y[t] == j) {
		nj++;
	    }
	}	
	p += (double) nj / pmod->nobs;
	OC->theta[i] = normal_cdf_inverse(p);
    }
}

static void op_LR_test (MODEL *pmod, op_container *OC, 
			const double **Z)
{
    int nx = OC->nx;
    double L0;

    OC->k -= OC->nx;
    OC->nx = 0;

    cut_points_init(OC, pmod, Z);

    L0 = op_loglik(OC->theta, OC);

    if (!na(L0) && L0 <= pmod->lnL) {
	pmod->chisq = 2.0 * (pmod->lnL - L0);
	gretl_model_set_int(pmod, "lr_df", nx);
    }

    /* restore original data on OC */
    OC->nx = nx;
    OC->k += nx;
}

static void fill_op_model (MODEL *pmod, const int *list,
			   const double **Z, const DATAINFO *pdinfo, 
			   op_container *OC, gretl_matrix *V, 
			   int fnc, int grc)
{
    int npar = OC->k;
    int nx = OC->nx;
    int correct = 0;
    double x;
    int i, s, t, v;

    pmod->ci = OC->type;

    gretl_model_set_int(pmod, "ordered", 1);
    gretl_model_set_int(pmod, "fncount", fnc);
    gretl_model_set_int(pmod, "grcount", grc);

    pmod->ncoeff = npar;

    if (V != NULL) {
	pmod->errcode = gretl_model_write_vcv(pmod, V);
	if (pmod->errcode) {
	    return;
	}
    }

    for (i=0; i<npar; i++) {
	pmod->coeff[i] = OC->theta[i];
    }

    if (OC->opt & OPT_R) {
	gretl_model_set_vcv_info(pmod, VCV_ML, VCV_QML);
	pmod->opt |= OPT_R;
    }

    s = 0;
    for (t=OC->t1; t<=OC->t2; t++) {
	int pred;

	if (na(OC->pmod->uhat[t])) {
	    continue;
	}

	x = 0.0;
	for (i=0; i<OC->nx; i++) {
	    v = OC->list[i+2];
	    x += OC->theta[i] * OC->Z[v][t];
	}

	/* yhat = X\hat{beta} */
	pmod->yhat[t] = x;
	pred = get_prediction(OC, pmod, x);
	if (pred == OC->y[s]) {
	    correct++;
	}
	/* compute generalized residual */
	pmod->uhat[t] = gen_resid(OC, OC->theta, s);
	s++;
    }

    gretl_model_set_int(pmod, "correct", correct);

    pmod->lnL = op_loglik(OC->theta, OC);
    mle_criteria(pmod, 0);

    gretl_model_allocate_params(pmod, npar);

    if (pmod->errcode == 0) {
	for (i=0; i<nx; i++) {
	    v = OC->list[i+2];
	    strcpy(pmod->params[i], pdinfo->varname[v]);
	}
	s = 1;
	for (i=nx; i<npar; i++) {
	    sprintf(pmod->params[i], "cut%d", s++);
	}
    }

    /* trim the model list: remove references to the 'cut'
       dummy variables */
    for (i=pmod->list[0]; i>1; i--) {
	if (!in_gretl_list(list, pmod->list[i])) {
	    gretl_list_delete_at_pos(pmod->list, i);
	}
    }

    if (nx > 0) {
	op_LR_test(pmod, OC, Z);
	gretl_model_set_coeff_separator(pmod, NULL, nx);
    }
}

static gretl_matrix *oprobit_vcv (op_container *OC, int *err)
{
    gretl_matrix *V;
    int k = OC->k;

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* hessian from analytical score */
    *err = ihess_from_ascore(OC, V);

    if (!*err && (OC->opt & OPT_R)) {
	gretl_matrix *GG = NULL;
	gretl_matrix *Vr = NULL;

	/* sandwich of hessian and OPG */

	GG = gretl_matrix_alloc(k, k);
	Vr = gretl_matrix_alloc(k, k);

	if (GG == NULL || Vr == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = opg_from_ascore(OC, GG);
	    if (!*err) {
#if ODEBUG > 1
		gretl_matrix_print(GG, "OPG matrix");
#endif
		gretl_matrix_qform(V, GRETL_MOD_NONE,
				   GG, Vr, GRETL_MOD_NONE);
		gretl_matrix_copy_values(V, Vr);
	    }
	}

	gretl_matrix_free(GG);
	gretl_matrix_free(Vr);
    } 

    if (*err) {
	gretl_matrix_free(V);
	V = NULL;
    }

#if ODEBUG > 1
    gretl_matrix_print(V, "Covariance matrix");
#endif

    return V;
}

/* Main ordered estimation function */

static int do_ordered (int ci, int ndum, 
		       double **Z, DATAINFO *pdinfo, 
		       MODEL *pmod, const int *list, 
		       gretlopt opt, PRN *prn)
{
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;
    op_container *OC;
    int i, npar;
    gretl_matrix *V = NULL;
    double *theta = NULL;
    int err;

    OC = op_container_new(ci, ndum, Z, pmod, opt);
    if (OC == NULL) {
	return E_ALLOC;
    }

    npar = OC->k;

    /* transformed theta to pass to BFGS */
    theta = malloc(npar * sizeof *theta);
    if (theta == NULL) {
	op_container_destroy(OC);
	return E_ALLOC;
    }

    /* initialize slopes */
    for (i=0; i<OC->nx; i++) {
	OC->theta[i] = 0.0001;
    }

    /* initialize cut points */
    cut_points_init(OC, pmod, (const double **) Z);

    /* transform theta to log-diff form */
    make_BFGS_theta(OC, theta);

#if ODEBUG
    for (i=0; i<npar; i++) {
	fprintf(stderr, "theta[%d]: 'real' = %g, BFGS = %g\n", i, OC->theta[i], theta[i]);
    }
    fprintf(stderr, "\ninitial loglikelihood = %.12g\n", 
	    op_loglik(theta, OC));
#endif

    err = BFGS_max(theta, npar, maxit, OPROBIT_TOL, 
		   &fncount, &grcount, op_loglik, C_LOGLIK,
		   op_score, OC, (prn != NULL)? OPT_V : OPT_NONE,
		   prn);

    if (err) {
	goto bailout;
    } else {
	fprintf(stderr, "Number of iterations = %d (%d)\n", fncount, grcount);
    }

    /* transform back to 'real' theta */
    get_real_theta(OC, theta);

    V = oprobit_vcv(OC, &err);

    if (!err) {
	fill_op_model(pmod, list, (const double **) Z, pdinfo, 
		      OC, V, fncount, grcount);
    }

 bailout:

    free(theta);
    gretl_matrix_free(V);
    op_container_destroy(OC);

    return err;
}

struct sorter {
    double x;
    int t;
};

/* We want to ensure that the values of the dependent variable
   actually used in the analysis (after dropping any bad
   observations) form a zero-based series of consecutive 
   integers.
*/

static int maybe_fix_ordered_depvar (MODEL *pmod, double **Z,
				     const DATAINFO *pdinfo,
				     double **orig_y, int *ndum)
{
    struct sorter *s = NULL;
    double *sy = NULL;
    double nexty, bady;
    int dv = pmod->list[1];
    int i, t, n = 0;
    int fixit = 0;
    int err = 0;

    /* make a copy of the observations on the dep. var. that
       were actually used in the initial OLS, recording
       the observation numbers, then sort
    */

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    n++;
	}
    }

    s = malloc(n * sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    s[i].x = Z[dv][t];
	    s[i].t = t;
	    i++;
	}
    }

    qsort(s, n, sizeof *s, gretl_compare_doubles);

    /* normalize to minimum zero */
    if (s[0].x != 0) {
	double ymin = s[0].x;

	for (i=0; i<n && s[i].x == ymin; i++) {
	    s[i].x = 0.0;
	}
	fixit = 1;
    }

    /* ensure that the sorted values increase by steps of one */
    for (i=1; i<n; i++) {
	if (s[i].x != s[i-1].x) {
	    nexty = s[i-1].x + 1;
	    if (s[i].x != nexty) {
		bady = s[i].x;
		while (i < n && s[i].x == bady) {
		    s[i++].x = nexty;
		}
		i--; /* compensate for outer i++ */
		fixit = 1;
	    } 
	} 
    }

    /* the number of dummies actually used will equal the
       max of normalized y
    */
    *ndum = (int) s[n-1].x;

    if (fixit) {
	/* the dependent var needs transforming */
	sy = copyvec(Z[dv], pdinfo->n);
	if (sy == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		sy[s[i].t] = s[i].x;
	    }
	    /* back up the original dep. var and replace it for
	       the duration of the ordered analysis */
	    *orig_y = Z[dv];
	    Z[dv] = sy;
	}
    }

    if (fixit) {
	fputs("oprobit: using normalized y\n", stderr);
    } else {
	fputs("oprobit: using original y\n", stderr);
    }

    free(s);

    return err;
}

static void restore_depvar (double **Z, double *y, int v)
{
    free(Z[v]);
    Z[v] = y;
}

static int *make_dummies_list (const int *list, 
			       double ***pZ, DATAINFO *pdinfo,
			       int *err)
{
    int *dumlist = gretl_list_new(1);

    if (dumlist == NULL) {
	*err = E_ALLOC;
    } else {
	dumlist[1] = list[1];

	/* OPT_F -> drop first value */
	*err = list_dumgenr(&dumlist, pZ, pdinfo, OPT_F);
	if (*err) {
	    free(dumlist);
	    dumlist = NULL;
	}
    }

    return dumlist;
}

static int *make_big_list (const int *list, double ***pZ,
			   DATAINFO *pdinfo, int **pdumlist,
			   int *err)
{
    int *dumlist;
    int *biglist;
    int i, k, nv;

    dumlist = make_dummies_list(list, pZ, pdinfo, err);
    if (dumlist == NULL) {
	return NULL;
    }

    nv = list[0] + dumlist[0];

    biglist = gretl_list_new(nv);
    if (biglist == NULL) {
	free(dumlist);
	*err = E_ALLOC;
	return NULL;
    }

    k = 1;
    for (i=1; i<=list[0]; i++) {
	biglist[k++] = list[i];
    }
    for (i=1; i<=dumlist[0]; i++) {
	biglist[k++] = dumlist[i];
    }

    *pdumlist = dumlist;

    return biglist;
}

static void list_purge_const (int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
	if (list[i] == 0) {
	    gretl_list_delete_at_pos(list, i);
	    return;
	}
    }
}

static int check_for_missing_dummies (MODEL *pmod, int *dlist)
{
    int i, dv;

    for (i=1; i<=dlist[0]; i++) {
	dv = dlist[i];
	if (!in_gretl_list(pmod->list, dv)) {
	    fprintf(stderr, "check for missing dummies: var %d is gone!\n", dv);
	    return E_SINGULAR;
	}
    }

    return 0;
}

/* the driver function for the plugin: note, prn is non-NULL
   only if the verbose option has been selected
*/

MODEL ordered_estimate (int *list, int ci, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn) 
{
    MODEL model;
    int orig_v = pdinfo->v;
    double *orig_y = NULL;
    int *biglist = NULL;
    int *dumlist = NULL;
    int ndum = 0;

    gretl_model_init(&model);

    /* remove the constant from the incoming list, if present */
    list_purge_const(list);

    /* construct augmented regression list, including dummies 
       for the level of the dependent variable
    */
    biglist = make_big_list(list, pZ, pdinfo, &dumlist, &model.errcode);
    if (model.errcode) {
	return model;
    }

    /* run initial OLS, with dummies added */
    model = lsq(biglist, pZ, pdinfo, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "ordered_estimate: initial OLS failed\n");
	free(dumlist);
	free(biglist);
	return model;
    }

#if ODEBUG
    pprintf(prn, "ordered_estimate: initial OLS\n");
    printmodel(&model, pdinfo, OPT_S, prn);
#endif

    if (model.list[0] < biglist[0]) {
	/* One or more regressors were dropped in OLS, most likely due
	   to collinearity.  We can accept the dropping of
	   user-selected regressors, but all the added dummies must be
	   retained, so we'd better investigate.
	*/
	model.errcode = check_for_missing_dummies(&model, dumlist);
    }
	
    if (!model.errcode) {
	/* after accounting for any missing observations, check that
	   the dependent variable is acceptable 
	*/
	model.errcode = maybe_fix_ordered_depvar(&model, *pZ, pdinfo,
						 &orig_y, &ndum);
    }

    /* do the actual ordered probit analysis */
    if (!model.errcode) {
	model.errcode = do_ordered(ci, ndum, *pZ, pdinfo, &model, list, 
				   opt, prn);
    }

    free(dumlist);
    free(biglist);

    if (orig_y != NULL) {
	/* if we messed with the dependent var, put the original back */
	restore_depvar(*pZ, orig_y, list[1]);
    }

    if (pdinfo->v > orig_v) {
	/* clean up any automatically-added dummies */
	dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);
    }

    return model;
}
