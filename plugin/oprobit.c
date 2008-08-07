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

#define ODEBUG 0

#define OPROBIT_TOL 1.0e-12

typedef struct op_container_ op_container;

struct op_container_ {
    int type;         /* model type */
    gretlopt opt;     /* option flags */
    int *y;           /* dependent variable */
    double **Z;       /* data */
    int *list;
    int M;            /* max of (possibly normalized) y */
    int t1;           /* beginning of sample */
    int t2;           /* end of sample */
    int nobs;         /* number of observations */
    int nx;           /* number of explanatory variables */
    int k;            /* total number of parameters */
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
    OC->M = ndum + 1;
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

    if (OC->y == NULL || OC->ndx == NULL || 
	OC->dP == NULL || OC->list == NULL ||
	OC->G == NULL || OC->g == NULL) {
	op_container_destroy(OC);
	return NULL;
    }

    i = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t])) {
	    OC->y[i++] = (int) Z[vy][t];
	}
    }

#if ODEBUG
    fprintf(stderr, "nobs = %d\n", OC->nobs);
    fprintf(stderr, "t1-t2 = %d-%d\n", OC->t1, OC->t2);
    fprintf(stderr, "k = %d\n", OC->k);
    fprintf(stderr, "nx = %d\n", OC->nx);
    fprintf(stderr, "Max(y) = M = %d\n", OC->M);
#endif

    OC->list[1] = vy;
    for (i=0; i<OC->nx; i++) {
	OC->list[i+2] = pmod->list[i+2];
    }

#if ODEBUG
    printlist(OC->list, "list, in op_container_new");
#endif

    return OC;
}

static int compute_probs (const double *theta, op_container *OC)
{
    double m0, m1, ystar0 = 0.0, ystar1 = 0.0;
    int M = OC->M;
    int nx = OC->nx;
    double P0, P1, h, adj, dP;
    int i, t, s, yt, v;
    int type = OC->type;

    if (OC->opt & OPT_A) {
	/* analytical score */
	for (i=0; i<OC->k; i++) {
	    OC->g[i] = 0.0;
	}
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
	    ystar1 = OC->ndx[s];
	} else if (yt == 1) {
	    m1 = theta[nx];
	    ystar0 = OC->ndx[s];
	    ystar1 = OC->ndx[s] + m1;
	} else {
	    m0 = theta[nx + yt - 2];
	    ystar0 = OC->ndx[s] + m0;
	    if (yt < M) {
		m1 = theta[nx + yt - 1];
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

	if (OC->opt & OPT_A) {
	    double mills0;
	    double mills1;
	    double dm;

	    if (ystar1 < 6.0 || OC->type == LOGIT) {
		mills0 = (yt == 0)? 0.0 : densfunc(ystar0, type) / dP;
		mills1 = (yt == M)? 0.0 : densfunc(ystar1, type) / dP;
	    } else { 
		/* L'Hopital-based approximation */
		mills0 = (yt == 0)? 0.0 : -ystar0;
		mills1 = (yt == M)? 0.0 : -ystar1;
	    }

	    dm = mills1 - mills0;

	    for (i=0; i<nx; i++) {
		v = OC->list[i+2];
		OC->G[i][s] = -dm * OC->Z[v][t];
		OC->g[i] += OC->G[i][s];
	    }

	    for (i=nx; i<OC->k; i++) {
		OC->G[i][s] = 0.0;
		if (i == nx + yt - 2) {
		    OC->G[i][s] = -mills0;
		    OC->g[i] += OC->G[i][s];
		}
		if (i == nx + yt - 1) {
		    OC->G[i][s] = mills1;
		    OC->g[i] += OC->G[i][s];
		}
	    }
	}

	s++;

    }

    return 0;
}

static int bad_cutpoints (const double *theta, const op_container *OC)
{
    int i;
    
    if (theta[OC->nx] <= 0.0) {
#if ODEBUG
	fprintf(stderr, "theta[%d] non-positive!\n", OC->nx);
#endif
	return 1;
    }
    
    for (i=OC->nx+1; i<OC->k; i++) {
	if (theta[i] <= theta[i-1]) {
#if ODEBUG
	    fprintf(stderr, "theta[%d]: non-increasing cutpoint!\n", i);
#endif
	    return 1;
	}
    }
	
    return 0;
}

static double op_loglik (const double *theta, void *ptr)
{
    op_container *OC = (op_container *) ptr;
    double x, ll = 0.0;
    int err;
    int i, s, t, v;
    
    if (bad_cutpoints(theta, OC)) {
	return NADBL;
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
	    x -= theta[i] * OC->Z[v][t];
	}
	OC->ndx[s++] = x;
#if ODEBUG > 1
	fprintf(stderr, "t = %d, s = %d, x = %g\n", t, s, x);
#endif
    }
    
    err = compute_probs(theta, OC);
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
    int i;

    for (i=0; i<npar; i++) {
	s[i] = OC->g[i];
    }

    return 1;
}

static int opg_from_ascore (op_container *OC, double *theta, gretl_matrix *GG) 
{
    int s, t, i, j, k = OC->k;
    double x, ll;

    double *g0 = malloc(k * sizeof *g0);

    if (g0 == NULL) {
	return E_ALLOC;
    }

    ll = op_loglik(theta, OC);

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

static int ihess_from_ascore (op_container *OC, double *theta, gretl_matrix *inH) 
{
    int i, j, err, k = OC->k;
    double smal = 1.0e-07;  /* "small" is some sort of macro on win32 */
    double smal2 = 2.0 * smal;
    double x, ll;

    double *g0 = malloc(k * sizeof *g0);

    if (g0 == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	theta[i] -= smal;
	ll = op_loglik(theta, OC);
	for (j=0; j<k; j++) {
	    g0[j] = OC->g[j];
	}
	theta[i] += smal2;
	ll = op_loglik(theta, OC);
	for (j=0; j<k; j++) {
	    x = (OC->g[j] - g0[j]) / smal2;
	    gretl_matrix_set(inH, i, j, -x);
	}
	theta[i] -= smal;
    }

    gretl_matrix_xtr_symmetric(inH);

    free(g0);

#if ODEBUG
    gretl_matrix_print(inH, "inverse of Hessian");
#endif

    err = gretl_invert_symmetric_matrix(inH);

    return err;
}

static int 
numerical_ihess (op_container *OC, double *theta, gretl_matrix *invH)
{
    double vij, *V;
    int i, j, k, npar = OC->k;
    int err = 0;

    if (invH == NULL) {
	return E_ALLOC;
    }

    V = numerical_hessian(theta, npar, op_loglik, OC, &err);
    if (V == NULL) {
	return err;
    }

    k = 0;
    for (i=0; i<npar; i++) {
	for (j=i; j<npar; j++) {
	    vij = V[k++];
	    gretl_matrix_set(invH, i, j, vij);
	    if (i != j) {
		gretl_matrix_set(invH, j, i, vij);
	    }
	}
    }

    free(V);

    return 0;
}

static int get_pred (op_container *OC, const MODEL *pmod,
		     double Xb)
{
    int i, j, pred = 0;

    for (j=OC->M, i=OC->k-1; i>=OC->nx; i--, j++) {
	if (Xb >= pmod->coeff[i]) {
	    pred = j;
	    break;
	}
    }

    if (pred == 0 && Xb >= pmod->coeff[0]) {
	pred = 1;
    }

    return pred;
}

static double gen_resid (const double *theta, op_container *OC, int t) 
{
    double m0, m1, ystar0, f0, f1;
    double ret, ystar1 = 0.0;
    int ndxt = OC->ndx[t];
    int dP = OC->dP[t];
    int M = OC->M;
    int nx = OC->nx;
    int yt = OC->y[t];

    if (yt == 0) {
	ystar1 = ndxt;
    } else if (yt == 1) {
	m1 = theta[nx];
	ystar0 = ndxt;
	ystar1 = ndxt + m1;
    } else {
	m0 = theta[nx + yt - 2];
	ystar0 = ndxt + m0;
	if (yt < M) {
	    m1 = theta[nx + yt - 1];
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

static void fill_model (MODEL *pmod, const DATAINFO *pdinfo, 
			op_container *OC, double *theta, 
			gretl_matrix *V)
{
    int npar = OC->k;
    int nx = OC->nx;
    double x;
    int i, k, s, t, v;

    pmod->t1 = OC->t1;
    pmod->t2 = OC->t2;
    pmod->nobs = OC->nobs;
    pmod->ci = OC->type;
    gretl_model_set_int(pmod, "ordered", 1);

    pmod->ncoeff = npar;

    if (V != NULL) {
	pmod->errcode = gretl_model_write_vcv(pmod, V);
	if (pmod->errcode) {
	    return;
	}
    }

    for (i=0; i<npar; i++) {
	pmod->coeff[i] = theta[i];
#if ODEBUG > 1
	fprintf(stderr,"theta[%d] = %12.8f\n", i, theta[i]);
#endif
    }

    if (OC->opt & OPT_R) {
	gretl_model_set_int(pmod, "ml_vcv", VCV_QML);
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
	    x += theta[i] * OC->Z[v][t];
	}

	/* X\hat{beta} */
	pmod->yhat[t] = x;
	pred = get_pred(OC, pmod, x); /* should we do anything with this? */
	/* compute generalized residual */
	pmod->uhat[t] = gen_resid(theta, OC, s++);
    }

    pmod->lnL = op_loglik(theta, OC);

    gretl_model_allocate_params(pmod, npar);

    if (pmod->errcode == 0) {
	for (i=0; i<nx; i++) {
	    strcpy(pmod->params[i], pdinfo->varname[OC->list[i+2]]);
	}
	k = 1;
	for (i=nx; i<npar; i++) {
	    sprintf(pmod->params[i], "cut%d", k++);
	}
    }
}

static gretl_matrix *oprobit_vcv (op_container *OC, double *theta, int *err)
{
    gretl_matrix *V;
    int k = OC->k;

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (OC->opt & OPT_A) {
	/* hessian from analytical score */
	*err = ihess_from_ascore(OC, theta, V);
    } else {
	/* numerical hessian */
	*err = numerical_ihess(OC, theta, V);
    }

    if (!*err && (OC->opt & OPT_R)) {
	gretl_matrix *GG = NULL;
	gretl_matrix *Vr = NULL;

	/* sandwich of hessian and OPG */

	GG = gretl_matrix_alloc(k, k);
	Vr = gretl_matrix_alloc(k, k);

	if (GG == NULL || Vr == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = opg_from_ascore(OC, theta, GG);
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

    return V;
}

/* Main ordered estimation function */

static int do_ordered (int ci, int ndum, 
		       double **Z, DATAINFO *pdinfo, MODEL *pmod,
		       gretlopt opt, PRN *prn)
{
    op_container *OC;
    int i, npar;
    gretl_matrix *V = NULL;
    double *theta = NULL;
    int err;

    /* BFGS apparatus */
    int maxit = 1000;
    int fncount = 0;
    int grcount = 0;

    /* set analytical score option... */
    opt |= OPT_A;

    /* ...but if not analytical, remove robust option */
    if (!(opt & OPT_A)) {
	opt &= ~OPT_R;
    }

    OC = op_container_new(ci, ndum, Z, pmod, opt);
    if (OC == NULL) {
	return E_ALLOC;
    }

    npar = OC->k;

    theta = malloc(npar * sizeof *theta);
    if (theta == NULL) {
	op_container_destroy(OC);
	return E_ALLOC;
    }

    for (i=0; i<OC->nx; i++) {
	theta[i] = pmod->coeff[i];
#if ODEBUG
	fprintf(stderr, "x: init'd theta[%d] = ols[%d] = %g\n", i, i, theta[i]);
#endif
    }

    /* 
       the following is very heuristic, has no sound theoretical 
       justification but seems to work in practice. So much
       for the Scientific Method (TM) !!!
    */
 
    for (i=OC->nx; i<npar; i++) {
	theta[i] = 0.5 * pmod->coeff[i];
#if ODEBUG
	fprintf(stderr, "cuts: ols[%d] = %g -> theta[%d] = %g\n", 
		i, pmod->coeff[i], i, theta[i]);
#endif
    }

#if ODEBUG
    fprintf(stderr, "\ninitial loglikelihood = %.12g\n", 
	    op_loglik(theta, OC));
#endif

    err = BFGS_max(theta, npar, maxit, OPROBIT_TOL, 
		   &fncount, &grcount, op_loglik, C_LOGLIK,
		   (OC->opt & OPT_A)? op_score : NULL, 
		   OC, (prn != NULL)? OPT_V : OPT_NONE,
		   prn);

    if (err) {
	goto bailout;
    } else {
	fprintf(stderr, "Number of iterations = %d (%d)\n", fncount, grcount);
    }

    for (i=0; i<npar; i++) {
	pmod->coeff[i] = theta[i];
    }

    V = oprobit_vcv(OC, theta, &err);

#if ODEBUG > 1
    gretl_matrix_print(V, "Covariance matrix");
#endif

    if (!err) {
	fill_model(pmod, pdinfo, OC, theta, V);
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
				     double **orig_y)
{
    struct sorter *s = NULL;
    double *sy = NULL;
    int dv = pmod->list[1];
    int i, t, n = 0;
    int fixit = 0;
    int err = 0;

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

    /* ensure that the sorted values increase by one */
    for (i=1; i<n; i++) {
	if (s[i].x != s[i-1].x) {
	    double nexty = s[i-1].x + 1;

	    if (s[i].x != nexty) {
		double bady = s[i].x;

		while (i < n && s[i].x == bady) {
		    s[i++].x = nexty;
		}
		i--; /* compensate for outer i++ */
		fixit = 1;
	    } 
	} 
    }

    if (fixit) {
	sy = copyvec(Z[dv], pdinfo->n);
	if (sy == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		sy[s[i].t] = s[i].x;
	    }
	    /* back up the original dep. var and replace it */
	    *orig_y = Z[dv];
	    Z[dv] = sy;
	}
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
			   DATAINFO *pdinfo, int *ndum,
			   int *err)
{
    int *dumlist;
    int *biglist;
    int i, k, nv;

    dumlist = make_dummies_list(list, pZ, pdinfo, err);
    if (dumlist == NULL) {
	return NULL;
    }

    /* In dumlist, the first value will have been dropped
       automatically, but we want to drop the second one too.
    */
    *ndum = dumlist[0] - 1;
    nv = list[0] + *ndum;

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
    for (i=2; i<=dumlist[0]; i++) {
	biglist[k++] = dumlist[i];
    }

    free(dumlist);

    return biglist;
}

static int list_has_const (const int *list)
{
    int i;

    for (i=2; i<=list[0]; i++) {
	if (list[i] == 0) {
	    return 1;
	}
    }

    return 0;
}

/* the driver function for the plugin: note, prn is non-NULL
   only if the verbose option has been selected
*/

MODEL ordered_estimate (const int *list, int ci, double ***pZ, DATAINFO *pdinfo,
			gretlopt opt, PRN *prn) 
{
    MODEL model;
    double *orig_y = NULL;
    int *biglist = NULL;
    int ndum = 0;

    gretl_model_init(&model);

    /* check for obligatory constant in incoming list */
    if (!list_has_const(list)) {
	model.errcode = E_NOCONST;
	return model;
    }

    /* construct augmented regression list, including
       dummies for the level of the dependent variable
    */
    biglist = make_big_list(list, pZ, pdinfo, &ndum, &model.errcode);
    if (model.errcode) {
	return model;
    }

    /* run initial OLS, with dummies added */
    model = lsq(biglist, pZ, pdinfo, OLS, OPT_A);
    if (model.errcode) {
	fprintf(stderr, "ordered_estimate: initial OLS failed\n");
	goto bailout;
    }

#if ODEBUG > 1
    pprintf(prn, "oprobit_estimate: initial OLS\n");
    printmodel(&model, pdinfo, OPT_S, errprn);
#endif

    /* after accounting for any missing observations, check that
       the dependent variable is OK 
    */
    model.errcode = maybe_fix_ordered_depvar(&model, *pZ, pdinfo,
					     &orig_y);

    if (!model.errcode && orig_y != NULL) {
	/* we modified the dependent variable: re-initialize (?) */
	clear_model(&model);
	model = lsq(biglist, pZ, pdinfo, OLS, OPT_A);
	if (model.errcode) {
	    fprintf(stderr, "ordered_estimate: secondary OLS failed\n");
	    goto bailout;
	}
    }

    /* do the actual ordered probit analysis */
    if (!model.errcode) {
	model.errcode = do_ordered(ci, ndum, *pZ, pdinfo, &model, opt, prn);
    }

 bailout:

    if (orig_y != NULL) {
	restore_depvar(*pZ, orig_y, list[1]);
    }

    free(biglist);

    return model;
}
