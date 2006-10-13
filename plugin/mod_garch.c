/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2006 Allin Cottrell and Riccardo "Jack" Lucchetti
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

/*
  Modular GARCH routines by Jack Lucchetti, October 2006. For the 
  moment, meant to replace seamlessly fcp.c, but syntax should evolve in 
  the future.
*/

#include "libgretl.h"
#include "libset.h"
#include "gretl_matrix.h"
#include "mod_garch.h"
#include "nls.h"

#include "f2c.h"
#include "clapack_double.h"

#define GDEBUG 0

enum {
    INIT_VAR_THEO,
    INIT_VAR_OLS,
    INIT_VAR_RESID
};

enum {
    DIST_NORM,
    DIST_T
};

typedef struct garch_container_ garch_container;

struct garch_container_ {
    double *y;
    const double **X;
    int t1;
    int t2;
    int nobs;
    int ncm;
    int p;
    int q;
    int k;
    int init;
    int distrib;
    double *e;
    double *e2;
    double *h;
    int ascore;
    double **score_e;
    double **score_h;
    double **blockglue;
    double **G;
    double *tot_score;
};

#if 0
static void mark (int *n) 
{ 
    if (n == NULL) {
	fprintf(stderr,"Ha!\n"); 
    } else {
	int locn = *n;

	fprintf(stderr,"Ha! (%d)\n", locn); 
	*n = ++locn;
    }
}
#endif

#define MOD_DEBUG 0

static void free_eh_derivs (garch_container *DH)
{
    doubles_array_free(DH->score_e, DH->k);
    doubles_array_free(DH->score_h, DH->k);
    doubles_array_free(DH->G, DH->k);
    doubles_array_free(DH->blockglue, 2);
}

static int init_eh_derivs (garch_container *DH)
{
    DH->score_e = NULL;
    DH->score_h = NULL;
    DH->G = NULL;
    DH->blockglue = NULL;
    int err = 0;

    if (DH->ascore) {
	DH->score_e = doubles_array_new(DH->k, DH->nobs);
	DH->score_h = doubles_array_new(DH->k, DH->nobs);
	DH->G = doubles_array_new(DH->k, DH->nobs);
	DH->blockglue = doubles_array_new(2, DH->nobs);

	if (DH->score_e == NULL ||
	    DH->score_h == NULL ||
	    DH->G == NULL ||
	    DH->blockglue == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static garch_container *init_container (double *y, const double **X, 
					int t1, int t2, int nobs, int nx,
					int p, int q, int init_method, 
					double *res, double *h, 
					int analytical)
{
    garch_container *DH = malloc(sizeof *DH);

    if (DH == NULL) {
	return NULL;
    }

    DH->e2 = malloc(nobs * sizeof *DH->e2);
    if (DH->e2 == NULL) {
	free(DH);
	return NULL;
    }

    DH->y = y;
    DH->X = X;
    DH->t1 = t1;
    DH->t2 = t2;
    DH->nobs = nobs;
    DH->ncm = nx;
    DH->p = p;
    DH->q = q;
    DH->init = init_method;
    DH->e = res;
    DH->h = h;
    DH->ascore = analytical;
    DH->k = 1 + nx + 1 + p + q;

    init_eh_derivs(DH);

    return DH;
}

static void garch_container_destroy (garch_container *DH)
{
    if (DH->ascore) {
	free_eh_derivs(DH);
    }
    free(DH->e2);
    free(DH);
}

/* 
   Compute the *ARCH log-likelihood for Gaussian innovations. 
*/

static double normal_ll (const garch_container *DH, int *err)
{
    double e2t, ht, ll = 0.0;
    int t;

    for (t=DH->t1; t<=DH->t2; t++) {
	e2t = DH->e2[t];
	ht = DH->h[t];
	if (na(e2t) || na(ht)) {
	    *err = E_MISSDATA;
	    return NADBL;
	}
	ll -= log(ht) + e2t / ht;
    }

    ll *= 0.5;
    ll -= (DH->t2 - DH->t1 + 1) * LN_SQRT_2_PI;

    return ll;
} 

static void normal_score (const garch_container *DH)
{
    double ut;
    int t;

    for (t=DH->t1; t<=DH->t2; t++) {
	DH->blockglue[0][t] = ut = -DH->e[t] / DH->h[t];
	DH->blockglue[1][t] = 0.5 * (ut * ut - 1.0 / DH->h[t]);
    }
} 

/* Compute the GARCH quantities */

static int garch_etht (const double *par, void *ptr)
{
    garch_container *DH = (garch_container *) ptr;

    int t1 = DH->t1;
    int t2 = DH->t2;
    int p = DH->p;
    int q = DH->q;

    int maxlag = (p > q)? p : q;

    int i, j, k, ret = 0;
    int ncm = DH->ncm;

    double **dedq = DH->score_e;
    double **dhdq = DH->score_h;

    int t, T = t2 - t1 + 1;
    double et, ht, tmp, h0 = 0.0;
    double u_var = 0.0;

    /* compute residuals */

    tmp = 0.0;
    for (t = t1-maxlag; t <= t2; t++) {
	if (t<t1) {
	    et = 0.0;
	} else {
	    et = DH->y[t] - par[0];
	    if (DH->X != NULL) {
		for(i=0; i<ncm; i++) {
		    et -= DH->X[i][t]*par[i+1];
		}
	    }
	    DH->e[t] = et;
	    DH->e2[t] = et * et;
	    tmp += DH->e2[t]; 
	}
    }

    if (DH->ascore) {
	for (t=t-maxlag; t<t1; t++) {
	    for (i=0; i<DH->k; i++) {
		dedq[i][t] = 0.0;
	    }
	}
    }
	
    /* h0 and derivatives */

    switch (DH->init) {
    case INIT_VAR_OLS:
	h0 = 1.0;
	break;
    case INIT_VAR_RESID:
	h0 = tmp / T;
	break;
    case INIT_VAR_THEO:
	tmp = 1.0;
	for (i=ncm+2; i<DH->k; i++) {
	    tmp -= par[i];
	}
	u_var = par[ncm+1] / tmp;
	h0 = u_var;
	break;
    }

    for (t=t1-maxlag; t<t1; t++) {
	DH->h[t] = h0;
	DH->e2[t] = h0;
    }

    if (DH->ascore) {
	double dh0;

	switch (DH->init) {
	case INIT_VAR_OLS:
	    for (t=t1-maxlag; t<t1; t++) {
		for (i=0; i<DH->k; i++) {
		    dhdq[i][t] = 0.0;
		}
	    }
	    break;

	case INIT_VAR_RESID:
	    dh0 = 0.0;
	    for (t=t1; t<=t2; t++) {
		dh0 -= DH->e[t];
	    }
	    for (t=t1-maxlag; t<t1; t++) {
		dhdq[0][t] = dh0  * 2.0 / T;
	    }
	    
	    for (i=0; i<ncm; i++) {
		dh0 = 0.0;
		for (t=t1; t<=t2; t++) {
		    dh0 -= DH->e[t] * DH->X[i][t];
		}
		for (t=t1-maxlag; t<t1; t++) {
		    dhdq[i+1][t] = dh0  * 2.0 / T;
		}
	    }
	    
	    for (t=t1-maxlag; t<t1; t++) {
		for (i=ncm+1; i<DH->k; i++) {
		    dhdq[i][t] = 0.0;
		}
	    }
	    
	    break;
	    
	case INIT_VAR_THEO:
	    for (t=t1-maxlag; t<t1; t++) {
		for (i=0; i<=ncm; i++) {
		    dhdq[i][t] = 0.0;
		}
	    }
	    dh0 = u_var / par[ncm+1];
	    for (t=t1-maxlag; t<t1; t++) {
		dhdq[ncm+1][t] = dh0;
	    }
	    dh0 *= u_var;
	    for (t=t1-maxlag; t<t1; t++) {
		for (i=ncm+2; i<DH->k; i++) {
		    dhdq[i][t] = dh0;
		}
	    }
	    break;
	}
	
    }

    /* in-sample loop */

    for (t=t1; t<=t2; t++) {
	ht = par[ncm+1];

	for (i=1; i<=p; i++) {
	    ht += DH->e2[t-i] * par[ncm+i+1];
	}

	for (i=1; i<=q; i++) {
	    ht += DH->h[t-i] * par[ncm+i+p+1];
	}
	
	DH->h[t] = ht;
	    
	if (DH->ascore) {
	    
	    /* constant */
	    dedq[0][t] = -1.0;
	    k = ncm+1;
	    dhdq[0][t] = 0.0;
	    for (i=1; i<=p; i++) {
		if (t - p < t1 && DH->init == INIT_VAR_RESID) {
		    dhdq[0][t] += par[k+i] * dhdq[0][t1-1];
		} else {	
		    dhdq[0][t] += 2.0 * par[k+i] * (DH->e[t-i]) * dedq[0][t-i];
		}
	    }
	    
	    /* regressors */
	    for (i=1; i<=ncm; i++) {
		dedq[i][t] = -(DH->X[i-1][t]);
		k = ncm+1;
		dhdq[i][t] = 0.0;
		for (j=1; j<=p; j++) {
		    if (t - p < t1 && DH->init == INIT_VAR_RESID) { // add INIT_THEO here
			dhdq[i][t] += par[k+j] * dhdq[i][t1-1];
		    } else {	
			dhdq[i][t] += 2.0 * par[k+j] * (DH->e[t-j]) * dedq[i][t-j];
		    }
		}
	    }
	    
	    /* garch params: omega */
	    dedq[ncm+1][t] = 0.0;
	    dhdq[ncm+1][t] = 1.0;
	    if (t - p < t1 && DH->init == INIT_VAR_THEO) {
		for (i=1; i<=p; i++) {
		    dhdq[ncm+1][t] += par[ncm+1+i] * dhdq[ncm+1][t1-1];
		}
	    }
	    
	    /* garch params: alphas */
	    k = ncm + 2;
	    for (i=1; i<=p; i++) {
		dedq[k][t] = 0.0;
		dhdq[k][t] = DH->e2[t-i];
		if (t - p < t1 && DH->init == INIT_VAR_THEO) {
		    for (j=0; j<p; j++) {
			dhdq[k][t] += par[k+j] * dhdq[k][t1-1];
		    }
		}
		k++;
	    }
	    
	    /* garch params: betas */
	    k = ncm + p + 2;
	    for (i=1; i<=q; i++) {
		dedq[k][t] = 0.0;
		dhdq[k][t] = DH->h[t-i];
		if (t - p < t1 && DH->init == INIT_VAR_THEO) {
		    for (j=0; j<p; j++) {
			dhdq[k][t] += par[k+j-p] * dhdq[k][t1-1];
		    }
		}
		k++;
	    }
	    
	    /* "real" recursive part */
	    for (i=0; i<DH->k; i++) {
		k = ncm + p + 2;
		for (j=1; j<=q; j++) {
		    dhdq[i][t] += par[k++] * dhdq[i][t-j];
		}
	    }
	    
	}	
    }
    
#if MOD_DEBUG
    fputs("\n\n", stderr);
    for (i=0; i<DH->k; i++) {
	fprintf(stderr, "garch_etht: par[%d] = %9.6f ", i, par[i]);
    }
    fputc('\n', stderr);
    for (t=t1-maxlag; t<=20; t++) {
	if (t < t1) {
	    fputc('*', stderr); 
	} else {
	    fputc(' ', stderr);
	}
	fprintf(stderr, " t:%4d ", t);
	fprintf(stderr, " %8.4f", DH->e[t]);
	fprintf(stderr, " %8.4f", DH->e2[t]);
	fprintf(stderr, " %8.4f", DH->h[t]);
	fprintf(stderr, " %12.8f", dedq[ncm+2][t]);
	fprintf(stderr, " %12.8f", dhdq[ncm+2][t]);
	fputc('\n', stderr);
    }
#endif

    return ret;
} 

static void garch_infomat (garch_container *DH, double **info)
{
    double tmpi, tmpj, tmpx1, tmpx2, x;
    int i, j, t;

    for (i=0; i<DH->k; i++) {
	for (j=0; j<DH->k; j++) {
	    info[i][j] = 0.0;
	}
    }

    for (t=DH->t1; t<=DH->t2; t++) {
	for (i=0; i<=DH->ncm; i++) {
	    tmpi = DH->score_h[i][t] / DH->h[t];
	    tmpx1 = (i==0)? 1.0 : DH->X[i-1][t];
	    tmpx1 *= 2.0;
	    tmpx1 /= DH->h[t];
	    for (j=0; j<=i; j++) {
		tmpj = DH->score_h[j][t] / DH->h[t];
		tmpx2 = (j==0)? 1.0 : DH->X[j-1][t];
		x = tmpx1 * tmpx2 + tmpi * tmpj;
		info[i][j] += x;
		info[j][i] += x;
	    }
	}

	for (i=DH->ncm+1; i<DH->k; i++) {
	    tmpi = DH->score_h[i][t] / DH->h[t];
	    for (j=DH->ncm+1; j<=i; j++) {
		tmpj = DH->score_h[j][t] / DH->h[t];
		x = tmpi * tmpj;
		info[i][j] += x;
		info[j][i] += x;
	    }
	}
    }

    for (i=0; i<DH->k; i++) {
	for (j=0; j<DH->k; j++) {
	    info[i][j] *= 0.5;
	}
    }

#if GDEBUG
    fputs("Information matrix:\n", stderr);
    for (i=0; i<DH->k; i++) {
	for (j=0; j<DH->k; j++) {
	    fprintf(stderr, "%14.6f ", info[i][j]);
	}
	fputc('\n', stderr);
    }
#endif
}

static void garch_opg (garch_container *DH, double **GG)
{
    double tmpi, x;
    int t, i, j;

    for (i=0; i<DH->k; i++) {
	for (j=0; j<DH->k; j++) {
	    GG[i][j] = 0.0;
	}
    }

    for (t=DH->t1; t<=DH->t2; t++) {
	for (i=0; i<DH->k; i++) {
	    tmpi = DH->G[i][t];
	    for (j=0; j<=i; j++) {
		x = tmpi * DH->G[j][t]; 
		GG[i][j] += x;
		GG[j][i] += x;
	    }
	}
    }

#if GDEBUG
    fputs("OPG matrix:\n", stderr);
    for (i=0; i<DH->k; i++) {
	for (j=0; j<DH->k; j++) {
	    fprintf(stderr, "%14.6f ", GG[i][j]);
	}
	fputc('\n', stderr);
    }
#endif
}

static double loglik (const double *theta, void *ptr)
{
    garch_container *DH = (garch_container *) ptr;
    int err;

    err = garch_etht(theta, DH);
    return normal_ll(DH, &err);
}

static int score_fill_matrices (const double *theta, void *ptr)
{
    garch_container *DH = (garch_container *) ptr;
    int i, t, err;

    err = garch_etht(theta, DH);
    normal_score(DH);

    for (t=DH->t1; t<=DH->t2; t++) {
	for (i=0; i<DH->k; i++) {
	    DH->G[i][t] = DH->score_e[i][t] * DH->blockglue[0][t] + 
		DH->score_h[i][t] * DH->blockglue[1][t];
	}
    }

    return err;
}

static int anal_score (double *theta, double *s, int npar, BFGS_LL_FUNC ll, 
		       void *ptr)
{
    garch_container *DH = (garch_container *) ptr;
    double tmp;
    int t, i, ret = 0;

    score_fill_matrices(theta, DH);

    for (i=0; i<npar; i++) {
	tmp = 0.0;
	for (t=DH->t1;t<=DH->t2; t++) {
	    tmp += DH->G[i][t];
	}
	s[i] = tmp;
    }

    return ret;
}

/*
   Parameters to garch_estimate_mod()

   t1:    beginning of sample in auxiliary database
   t2:    end of sample in auxiliary database
   nobs:  total number of observations in auxiliary database
   X:     data matrix for auxiliary database (regressors, not needed on
          output)
   nx:    number of columns of X
   coeff: vector of coefficient for the conditional mean, normally
          initialised by OLS on input (not needed on output) -- does NOT
	  include the constant
   nc:    number of elements in coeff
   vcv:   n^2 vector (0 on input) to store covariance matrix of coeff
   res2:  vector of 0's on input, squared resids on output (not needed)
   e:     vector of 0's on input, resids on output
   h:     null pointer on input, conditional variances on output
   y:     on input, vector with dep. var., not needed on output
   amax:  vector; element 0 holds the garch intercept; 1 and 2 the
          arch & garch orders; from 3 onwards, the arch & garch 
          parameters
   b:     0 on input, holds vector of coefficient for the conditional
          mean on output
   scale: double used to scale dep. var.
   iters: int, 0 on input, holds number of iterations on output
   prn:   print handle for info on iterations and other diagnostic output

*/

int garch_estimate_mod (int t1, int t2, int nobs, 
			const double **X, int nx, double *coeff, int nc, 
			double *vcv, double *res2, double *res, double *h,
			double *y, double *amax, double *b, 
			double scale, int *iters, PRN *prn, int vopt)
{
    garch_container *DH;
    int p = amax[1];
    int q = amax[2];
    int npar = 1 + nx + 1 + p + q;
    int analytical = 1;
    double *theta = NULL;
    double **info = NULL;
    double *V, vij;
    int i, j, k, err = 0;

    /* BFGS apparatus */
    int maxit = 1000;
    double reltol = 1.0e-12;
    int fncount = 0;
    int grcount = 0;

    DH = init_container(y, X, t1, t2, nobs, nx, p, q, INIT_VAR_RESID,
			res, h, analytical);

    theta = malloc(npar * sizeof *theta);
    theta[0] = gretl_mean(t1, t2, y);
    for (i=1; i<=nx; i++) {
	theta[i] = coeff[i];
    }
    theta[nx+1] = amax[0];
    for (i=nx+2, j=3; i<npar; i++) {
	theta[i] = amax[j++];
    }

    if (DH->ascore) {
	fputs("\nUsing analytical score\n", stderr);
    } else {
	fputs("\nUsing numerical score\n", stderr);
    } 

    err = BFGS_max(theta, npar, maxit, reltol, 
		   &fncount, &grcount, loglik, 
		   (DH->ascore)? anal_score : NULL, 
		   DH, (prn != NULL)? OPT_V : OPT_NONE, 
		   prn);
    
    amax[0] = loglik(theta, DH) - (t2 - t1 + 1) * log(scale);
    
    fprintf(stderr, "iters = %d (%d)\n", fncount, grcount);

#if GDEBUG
    if (1) {
	double *testa, *testn;
	int testret;

	testa = malloc(npar * sizeof *testa);
	testn = malloc(npar * sizeof *testn);

	testret = anal_score(theta, testa, npar, loglik, DH);
	testret = BFGS_numeric_gradient(theta, testn, npar, loglik, DH);

	fprintf(stderr, "ret = %d:\n", testret);

	for (i=0; i<npar; i++) {
	    fprintf(stderr, "g[%d]: analytical = %14.8f, numerical: = %14.8f, \n", 
		    i, testa[i], testn[i]);
	}
	fputc('\n', stderr);

	free(testa);
	free(testn);
    }
#endif

    V = numerical_hessian(theta, npar, loglik, DH);

    k = 0;
    for (i=0; i<npar; i++) {
	for (j=i; j<npar; j++) {
	    vij = V[k++];
	    vcv[(i*npar) + j] = vij;
	    if (i != j) {
		vcv[(j*npar) + i] = vij;
	    }
	}
    }

#if GDEBUG
    k = 0;
    for (i=0; i<npar; i++) {
	for (j=i; j<npar; j++) {
	    fprintf(stderr, "V[%d,%d] = %g\t", i, j, V[k++]);
	}
	fputc('\n', stderr);
    }
#endif

    info = doubles_array_new(npar, npar);
    garch_infomat(DH, info);
    garch_opg(DH, info);

    doubles_array_free(info, npar);
    free(V);

    if (!err) {
	double sderr;

	/* transcribe coefficients and standard errors 
	   note: slightly different from fcp
	*/
	for (i=0, j=0; i<npar; i++) {
	    if (vcv[j] > 0.0) {
		sderr = sqrt(vcv[j]);
	    } else {
		sderr = 0.0;
	    }
	    j += npar + 1;
	    amax[i+1] = theta[i];
	    amax[i+1+npar] = sderr;
	}
    }

    free(theta);
    garch_container_destroy(DH);

    return err;
}
