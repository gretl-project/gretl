/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* discrete.c for gretl: logit and probit models */

#include "libgretl.h"
#include "gretl_private.h"

#define TINY 1.0e-13

static int neginv (const double *xpx, double *diag, int nv);
static int cholesky_decomp (double *xpx, int nv);

/* .......................................................... */

static double logit (double x)
{
    return 1.0 / (1.0 + exp(-x));
}

static double logit_pdf (double x)
{
    double z = exp(-x);

    return z / ((1.0 + z) * (1.0 + z));
}

/* .......................................................... */

static void Lr_chisq (MODEL *pmod, double **Z)
{
    int t, zeros, ones = 0, m = pmod->nobs;
    double Lr, chisq;
    
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (floateq(Z[pmod->list[1]][t], 1.0)) {
	    ones++;
	} 
    }
    zeros = m - ones;

    Lr = (double) ones * log((double) ones/ (double) m);
    Lr += (double) zeros * log((double) zeros/(double) m);

    chisq = 2.0 * (pmod->lnL - Lr);
    gretl_model_set_double(pmod, "chisq", chisq);
    
    /* McFadden pseudo-R^2 */
    pmod->rsq = 1.0 - pmod->lnL / Lr;
    pmod->adjrsq = NADBL;
}

/* .......................................................... */

static double 
logit_probit_llhood (const double *y, const MODEL *pmod, int opt)
{
    double q, lnL = 0.0;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->yhat[t])) {
	    continue;
	}
	q = 2.0 * y[t] - 1.0;
	if (opt == LOGIT) {
	    lnL += log(logit(q * pmod->yhat[t]));
	} else {
	    lnL += log(normal_cdf(q * pmod->yhat[t]));
	}
    }

    return lnL;
}

/* .......................................................... */

static int add_slopes_to_model (MODEL *pmod, double fbx)
{
    double *slopes;
    size_t ssize = pmod->ncoeff * sizeof *slopes;
    int i;

    slopes = malloc(ssize);

    if (slopes == NULL) {
	return 1;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	if (pmod->list[i+2] == 0) {
	    continue;
	}
	slopes[i] = pmod->coeff[i] * fbx;
    }

    if (gretl_model_set_data(pmod, "slopes", slopes, ssize)) {
	free(slopes);
	return 1;
    }

    return 0;
}

/* .......................................................... */

int dmod_isdummy (const double *x, int t1, int t2)
{
    int t, m = 0, goodobs = 0;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	if (x[t] != 0.0 && x[t] != 1.0) {
	    return 0;
	}
	if (x[t] == 1.0) {
	    m++;
	}
	goodobs++;
    }

    if (m < goodobs) return m;

    return 0;
} 

/* .......................................................... */

static double *hess_wts (MODEL *pmod, double **Z, int opt) 
{
    int i, t, tm, n = pmod->t2 - pmod->t1 + 1;
    double q, bx, xx, *w;

    w = malloc(n * sizeof *w);
    if (w == NULL) {
	return NULL;
    }

    for (t=pmod->t1; t<=pmod->t2; t++) {
	tm = t - pmod->t1;
	if (model_missing(pmod, t)) {
	    w[tm] = NADBL;
	    continue;
	}
	bx = 0.0;
	q = 2.0 * Z[pmod->list[1]][t] - 1.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    bx += pmod->coeff[i] * Z[pmod->list[i+2]][t];
	}
	if (opt == LOGIT) {
	    w[tm] = -1.0 * logit(bx) * (1.0 - logit(bx));
	} else {
	    xx = (q * normal_pdf(q * bx)) / normal_cdf(q * bx);
	    w[tm] = -xx * (xx + bx);
	}
    }

    return w;
}

/* .......................................................... */

static double *hessian (MODEL *pmod, double **Z, int opt) 
{
    int i, j, li, lj, m, t;
    const int l0 = pmod->list[0];
    double xx, *wt, *xpx;

    i = l0 - 1;
    m = i * (i + 1) / 2;

    xpx = malloc(m * sizeof *xpx);
    if (xpx == NULL) {
	return NULL;
    }

    wt = hess_wts(pmod, Z, opt);
    if (wt == NULL) {
	free(xpx);
	return NULL;
    }

    m = 0;
    for (i=2; i<=l0; i++) {
	li = pmod->list[i];
	for (j=i; j<=l0; j++) {
	    lj = pmod->list[j];
	    xx = 0.0;
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!model_missing(pmod, t)) {
		    xx += wt[t-pmod->t1] * Z[li][t] * Z[lj][t];
		}
	    }
	    if (floateq(xx, 0.0) && li == lj) {
		free(xpx);
		free(wt);
		return NULL;
	    }
	    xpx[m++] = -xx;
	}
    }

    free(wt);

    return xpx; 
}

/**
 * logit_probit:
 * @list: dependent variable plus list of regressors.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: option flag: If = LOGIT, perform logit regression, otherwise
 * perform probit regression.
 *
 * Computes estimates of the discrete model specified by @list,
 * using an estimator determined by the value of @opt.  Uses the
 * EM algorithm; see Ruud.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL logit_probit (int *list, double ***pZ, DATAINFO *pdinfo, int opt)
{
    int i, t, v, depvar = list[1];
    int itermax = 200;
    int *dmodlist = NULL;
    int dummy, n_correct;
    double xx, zz, fx, Fx, fbx, Lbak;
    double *xbar, *diag, *xpx = NULL;
    MODEL dmod;

    gretl_model_init(&dmod);
    
    /* check whether depvar is binary */
    dummy = dmod_isdummy((*pZ)[depvar], pdinfo->t1, pdinfo->t2);
    if (!dummy) {
	dmod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, _("The dependent variable '%s' is not a 0/1 "
				"variable.\n"), pdinfo->varname[depvar]);
	return dmod;
    }

    dmodlist = malloc((list[0] + 1) * sizeof *dmodlist);
    if (dmodlist == NULL) {
	dmod.errcode = E_ALLOC;
	return dmod;
    } 

    /* allocate space for means of indep vars */
    xbar = malloc(list[0] * sizeof *xbar);
    if (xbar == NULL) {
	free(dmodlist);
	dmod.errcode = E_ALLOC;
	return dmod;
    }

    /* make room for special expected value */
    if (dataset_add_vars(1, pZ, pdinfo)) {
	free(xbar);
	free(dmodlist);
	dmod.errcode = E_ALLOC;
	return dmod;
    }

    v = pdinfo->v - 1; /* the last created variable */

    dmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0);
    if (dmod.ifc == 0) {
	dmod.errcode = E_NOCONST;
    }
    if (dmod.errcode) {
	dataset_drop_vars(1, pZ, pdinfo);
	free(xbar);
	free(dmodlist);
	return dmod;
    }

    for (i=2; i<=list[0]; i++) {
	dmodlist[i] = list[i];
	xbar[i-2] = 0.0;
	for (t=dmod.t1; t<=dmod.t2; t++) {
	    if (!model_missing(&dmod, t)) {
		xbar[i-2] += (*pZ)[dmod.list[i]][t];
	    }
	}
	xbar[i-2] /= dmod.nobs;
    }

    dmodlist[0] = list[0];
    dmodlist[1] = v; /* dep var is the newly created one */

    /* Iterated least squares: see Ruud, "An Introduction to Classical 
       Econometric Theory", chapter 27 */
    Lbak = -1.0e8;
    for (i=0; i<itermax; i++) {
	for (t=dmod.t1; t<=dmod.t2; t++) {
	    xx = dmod.yhat[t];
	    if (!na(dmod.yhat[t])) {
		if (opt == LOGIT) {
		    fx = logit_pdf(xx);
		    Fx = logit(xx);
		} else {
		    fx = normal_pdf(xx);
		    Fx = normal_cdf(xx);
		}
		if (floateq((*pZ)[depvar][t], 0.0)) {
		    xx -= fx / (1.0 - Fx);
		} else {
		    xx += fx / Fx;
		}
	    }
	    (*pZ)[v][t] = xx;
	}

	dmod.lnL = logit_probit_llhood(&(*pZ)[v][0], &dmod, opt);
	if (fabs(dmod.lnL - Lbak) < .000005) {
	    break; 
	}

	/*  printf("Log likelihood = %g\n", dmod.lnL); */
	Lbak = dmod.lnL;
	clear_model(&dmod);
	dmod = lsq(dmodlist, pZ, pdinfo, OLS, OPT_A, 0);
	if (dmod.errcode) {
	    fprintf(stderr, "logit_probit: dmod errcode=%d\n", dmod.errcode);
	    dataset_drop_vars(1, pZ, pdinfo);
	    free(xbar);
	    free(dmodlist);
	    return dmod;
	}
    }

    dataset_drop_vars(1, pZ, pdinfo);
    dmod.list[1] = depvar;
    dmod.lnL = logit_probit_llhood((*pZ)[depvar], &dmod, opt);
    Lr_chisq(&dmod, *pZ);
    dmod.ci = opt;
    free(dmodlist);

    /* in case using QR, where vcv is precalculated */
    if (dmod.vcv != NULL) {
	free(dmod.vcv);
	dmod.vcv = NULL;
    }

    /* form the Hessian */
    if (dmod.xpx != NULL) {
	free(dmod.xpx);
    }
    dmod.xpx = hessian(&dmod, *pZ, opt);
    if (dmod.xpx == NULL) {
	free(xbar);
	dmod.errcode = E_ALLOC;
	strcpy(gretl_errmsg, _("Failed to construct Hessian matrix"));
	return dmod;
    } 

    /* obtain negative inverse of Hessian */
    cholesky_decomp(dmod.xpx, dmod.ncoeff); 
    diag = malloc(dmod.ncoeff * sizeof *diag); 
    if (diag == NULL) {
	free(xbar);
	dmod.errcode = E_ALLOC;
	return dmod;
    }

    xpx = copyvec(dmod.xpx, dmod.ncoeff * (dmod.ncoeff + 1) / 2);
    neginv(xpx, diag, dmod.ncoeff);

    /* std errors: square roots of diagonal elements */
    for (i=0; i<dmod.ncoeff; i++) {
	dmod.sderr[i] = sqrt(diag[i]);
    }
    free(diag);
    free(xpx);

    /* apparatus for calculating slopes at means */
    xx = 0.0;
    for (i=0; i<dmod.ncoeff; i++) {
	xx += dmod.coeff[i] * xbar[i];
    }
    free(xbar);

    if (opt == LOGIT) {
	fbx = logit_pdf(xx);
    } else {
	fbx = normal_pdf(xx);
    }

    if (add_slopes_to_model(&dmod, fbx)) {
	dmod.errcode = E_ALLOC;
	return dmod;
    }

    /* calculate additional statistics */
    xx = 0.0;
    n_correct = 0;
    for (t=dmod.t1; t<=dmod.t2; t++) {
	if (model_missing(&dmod, t)) {
	    continue;
	}
	zz = (*pZ)[depvar][t];
	xx += zz;
	n_correct += ((dmod.yhat[t] > 0.0 && floateq(zz, 1.0)) ||
		      (dmod.yhat[t] <= 0.0 && floateq(zz, 0.0)));
    }

    xx /= dmod.nobs;
    dmod.ybar = xx;
    dmod.sdy = fbx;
    gretl_model_set_int(&dmod, "correct", n_correct);

    mle_aic_bic(&dmod, 0);

    dmod.ID = model_count_plus();

    return dmod;
}

/* Solves for diagonal elements of X'X inverse matrix.
   X'X must be Cholesky-decomposed already.
*/

static int neginv (const double *xpx, double *diag, int nv)
{
    int kk, l, m, k, i, j;
    const int nxpx = nv * (nv + 1) / 2;
    double d, e, *tmp;

    tmp = malloc((nv + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return 1;
    }

    for (i=0; i<=nv; i++) {
	tmp[i] = 0.0;
    }

    kk = 0;

    for (l=1; l<=nv-1; l++) {
        d = xpx[kk];
        tmp[l] = d;
        e = d * d;
        m = 0;
        if (l > 1) 
	    for (j=1; j<=l-1; j++) m += nv - j;
        for (i=l+1; i<=nv; i++) {
            d = 0.0;
            k = i + m - 1;
            for (j=l; j<=i-1; j++) {
                d += tmp[j] * xpx[k];
                k += nv - j;
            }
            d = (-1.0) * d * xpx[k];
            tmp[i] = d;
            e += d * d;
        }
        kk += nv + 1 - l;
        diag[l-1] = e;
    }

    diag[nv-1] = xpx[nxpx-1] * xpx[nxpx-1];

    free(tmp);

    return 0;
}

/* Cholesky decomposition of X'X */

static int cholesky_decomp (double *xpx, int nv)
{
    int i, j, k, kk, l, jm1;
    double e, d, d1, test, xx;

    e = 1.0 / sqrt(xpx[0]);
    xpx[0] = e;

    for (i=1; i<nv; i++) {
	xpx[i] *= e;
    }

    kk = nv;

    for (j=2; j<=nv; j++) {
	/* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d += xx * xx;
            k += nv-l;
        }
        test = xpx[kk] - d;
        if (test / xpx[kk] < TINY) {
	    return 1;
	}
        e = 1 / sqrt(test);
        xpx[kk] = e;
        /* off-diagonal elements */
        for (i=j+1; i<=nv; i++) {
            kk++;
            d = 0.0;
            k = j - 1;
            for (l=1; l<=jm1; l++) {
                d += xpx[k] * xpx[k-j+i];
                k += nv - l;
            }
            xpx[kk] = (xpx[kk] - d) * e;
        }
        kk++;
    }

    return 0; 
}

