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
#include "internal.h"

#define TINY 1.0e-13

static double *hessian (MODEL *pmod, const double *Z, const int n, 
			const int opt); 
static double *hess_wts (MODEL *pmod, const double *Z, const int n, 
			 const int opt);
static int neginv (double *xpx, double *diag, int nv);
static int choleski (double *xpx, int nv);

/* .......................................................... */

static double _norm_pdf (const double xx)
{
    return (1.0/sqrt(2.0 * M_PI)) * exp(-0.5 * xx * xx);
}

/* .......................................................... */

static double _norm_cdf (const double xx)
{
    return 1.0 - normal(xx);
}

/* .......................................................... */

static double _logit (double xx)
{
    return 1.0/(1.0 + exp(-xx));
}

/* .......................................................... */

static double _logit_pdf (double xx)
{
    double zz = exp(-xx);

    return zz/((1.0 + zz) * (1.0 + zz));
}

/* .......................................................... */

static void Lr_chisq (MODEL *pmod, const double *Z, const int n)
{
    int t, zeros, ones = 0, m = pmod->nobs;
    double Lr;
    
    for (t=pmod->t1; t<=pmod->t2; t++) 
	if (floateq(Z[n*pmod->list[1] + t], 1.0)) ones++;
    zeros = m - ones;
    Lr = (double) ones * log((double) ones/ (double) m);
    Lr += (double) zeros * log((double) zeros/(double) m);
    pmod->chisq = 2.0 * (pmod->lnL - Lr);
}

/* .......................................................... */

static double _logit_probit_llhood (double *y, MODEL *pmod, int opt)
{
    int t;
    double q, lnL = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	q = 2.0 * y[t] - 1.0;
	if (opt == LOGIT) 
	    lnL += log(_logit(q * pmod->yhat[t]));
	else
	    lnL += log(_norm_cdf(q * pmod->yhat[t]));
    }
    return lnL;
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
 * using an estimator determined by the value of @opt.
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */


MODEL logit_probit (int *list, double **pZ, DATAINFO *pdinfo, int opt)
     /* EM algorithm, see Ruud */
{
    int i, t, v, depvar = list[1];
    int n = pdinfo->n, itermax = 200;
    double xx, zz, fx, Fx, fbx, Lbak;
    double *xbar, *diag, *xpx = NULL;
    MODEL dmod;

    init_model(&dmod);

    /* check that depvar is really a dummy */
    if (isdummy(depvar, pdinfo->t1, pdinfo->t2, *pZ, n) == 0) {
	dmod.errcode = E_UNSPEC;
	sprintf(gretl_errmsg, "The dependent variable '%s' is not a 0/1 "
		"variable.\n", pdinfo->varname[depvar]);
	return dmod;
    }

    /* allocate space for means of indep vars */
    xbar = malloc((list[0] - 1) * sizeof *xbar);
    if (xbar == NULL) {
	dmod.errcode = E_ALLOC;
	return dmod;
    }

    /* make room for special expected value */
    if (grow_Z(1, pZ, pdinfo)) {
	free(xbar);
	dmod.errcode = E_ALLOC;
	return dmod;
    }
    v = pdinfo->v - 1;

    dmod = lsq(list, pZ, pdinfo, OLS, 0, 0);
    if (dmod.ifc == 0) dmod.errcode = E_NOCONST;
    if (dmod.errcode) {
	(void) shrink_Z(1, pZ, pdinfo);
	free(xbar);
	return dmod;
    }
    for (i=2; i<=list[0]; i++) {
	xbar[i-2] = 0.0;
	for (t=dmod.t1; t<=dmod.t2; t++)
	    xbar[i-2] += (*pZ)[n*dmod.list[i] + t];
	xbar[i-2] /= dmod.nobs;
    }
    list[1] = v;

    /* Iterated least squares: see Ruud, "An Introduction to Classical 
       Econometric Theory", chapter 27 */
    Lbak = -9999.0;
    for (i=0; i<itermax; i++) {
	for (t=dmod.t1; t<=dmod.t2; t++) {
	    xx = dmod.yhat[t];
	    if (opt == LOGIT) {
		fx = _logit_pdf(xx);
		Fx = _logit(xx);
	    } else {
		fx = _norm_pdf(xx);
		Fx = _norm_cdf(xx);
	    }
	    if (floateq((*pZ)[n*depvar + t], 0.0)) 
		xx -= fx/(1.0 - Fx);
	    else 
		xx += fx/Fx;
	    (*pZ)[n*v + t] = xx;
	}
	dmod.lnL = _logit_probit_llhood(&(*pZ)[n*v], &dmod, opt);
	if (fabs(dmod.lnL - Lbak) < .000005) break; 
	/*  printf("Log likelihood = %f\n", dmod.lnL); */
	Lbak = dmod.lnL;
	clear_model(&dmod, NULL, NULL);
	dmod = lsq(list, pZ, pdinfo, OLS, 0, 0);
	if (dmod.errcode) {
	    (void) shrink_Z(1, pZ, pdinfo);
	    free(xbar);
	    return dmod;
	}
    }
    /* put back original dependent variable */
    dmod.list[1] = depvar;
    shrink_Z(1, pZ, pdinfo);
    dmod.lnL = _logit_probit_llhood(&(*pZ)[n*depvar], &dmod, opt);
    Lr_chisq(&dmod, *pZ, n);
    dmod.ci = opt;

    /* form the Hessian */
    if (dmod.xpx != NULL) free(dmod.xpx);
    dmod.xpx = hessian(&dmod, *pZ, n, opt);
    if (dmod.xpx == NULL) {
	free(xbar);
	dmod.errcode = E_ALLOC;
	strcpy(gretl_errmsg, "Failed to construct Hessian matrix");
	return dmod;
    } 
    /* obtain negative inverse of Hessian */
    choleski(dmod.xpx, dmod.ncoeff); 
    diag = malloc((dmod.ncoeff + 1) * sizeof(double)); 
    if (diag == NULL) {
	free(xbar);
	dmod.errcode = E_ALLOC;
	return dmod;
    }
    xpx = copyvec(dmod.xpx, 1 + dmod.ncoeff * (dmod.ncoeff + 1) / 2);
    neginv(xpx, diag, dmod.ncoeff);
    /* std errors: square roots of diagonal elements */
    for (i=1; i<=dmod.ncoeff; i++) 
	dmod.sderr[i] = sqrt(diag[i]);
    free(diag);
    free(xpx);

    /* apparatus for calculating slopes at means */
    xx = 0.0;
    for (i=0; i<dmod.ncoeff; i++) {
	xx += dmod.coeff[i+1] * xbar[i];
    }
    free(xbar);
    if (opt == LOGIT)
	fbx = _logit_pdf(xx);
    else
	fbx = _norm_pdf(xx);
    dmod.slope = malloc((dmod.ncoeff + 1) * sizeof(double));
    if (dmod.slope == NULL) {
	dmod.errcode = E_ALLOC;
	return dmod;
    }
    for (i=0; i<dmod.ncoeff; i++) {
	if (dmod.list[i+2] == 0) continue;
	dmod.slope[i+1] = dmod.coeff[i+1] * fbx;
    }

    /* calculate additional statistics */
    xx = 0.0;
    dmod.correct = 0;
    for (t=dmod.t1; t<=dmod.t2; t++) {
	zz = (*pZ)[n*depvar + t];
	xx += zz;
	dmod.correct += ((dmod.yhat[t] >= 0.5 && floateq(zz, 1.0)) ||
		    (dmod.yhat[t] < 0.5 && floateq(zz, 0.0)));
    }
    xx /= dmod.nobs;
    dmod.ybar = xx;
    dmod.sdy = fbx;

    return dmod;
}

/* .......................................................... */

static double *hessian (MODEL *pmod, const double *Z, const int n, 
			const int opt) 
{
    int i, j, li, lj, l0 = pmod->list[0], m, t;
    double xx, *wt, *xpx;

    i = l0 - 1;
    m = i * (i + 1) / 2;
    xpx = malloc((m+1) * sizeof *xpx);
    if (xpx == NULL) return NULL;
    if ((wt = hess_wts(pmod, Z, n, opt)) == NULL) {
	free(xpx);
	return NULL;
    }
    m = 0;
    for (i=2; i<=l0; i++) {
	li = pmod->list[i];
	for (j=i; j<=l0; j++) {
	    lj = pmod->list[j];
	    xx = 0.0;
	    for (t=pmod->t1; t<=pmod->t2; t++) 
		xx += wt[t] * Z(li, t) * Z(lj, t);
	    if (floateq(xx, 0.0) && li == lj) {
		free(xpx);
		free(wt);
		return NULL;
	    }
	    xpx[++m] = -xx;
	}
    }
/*      for (i=1; i<=m; i++)  */
/*    	printf("xpx[%d] = %10.6f\n", i, xpx[i]); */
    free(wt);
    return xpx; 
}

/* .......................................................... */

static double *hess_wts (MODEL *pmod, const double *Z, const int n, 
			 const int opt) 
{
    int i, t;
    double q, bx, xx, *wt;

    wt = malloc(n * sizeof *wt);
    if (wt == NULL) return NULL;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	bx = 0.0;
	q = 2.0 * Z(pmod->list[1], t) - 1.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    bx += pmod->coeff[i+1] * Z(pmod->list[i+2], t);
	}
	if (opt == LOGIT) {
	    wt[t] = -1.0 * _logit(bx) * (1.0 - _logit(bx));
	} else {
	    xx = (q * _norm_pdf(q * bx)) / _norm_cdf(q * bx);
	    wt[t] = -xx * (xx + bx);
	}
    }
    return wt;
}

/* ...............................................................    */

static int neginv (double *xpx, double *diag, int nv)
/*
  Solves for diagonal elements of X'X inverse matrix.
  X'X must be Choleski-decomposed already.
*/
{
    int kk = 1, l, m, nstop, k, i, j;
    double d, e, *tmp;

    tmp = malloc((nv + 1) * sizeof *tmp);
    if (tmp == NULL) return 1;
    for (i=0; i<=nv; i++) tmp[i] = 0.0;

    nstop = nv * (nv + 1)/2;
    for (l=1; l<=nv-1; l++) {
        d = xpx[kk];
        tmp[l] = d;
        e = d * d;
        m = 0;
        if (l > 1) 
	    for (j=1; j<=l-1; j++) m += nv - j;
        for (i=l+1; i<=nv; i++) {
            d = 0.0;
            k = i + m;
            for (j=l; j<=i-1; j++) {
                d += tmp[j] * xpx[k];
                k += nv - j;
            }
            d = (-1.0) * d * xpx[k];
            tmp[i] = d;
            e += d * d;
        }
        kk += nv + 1 - l;
        diag[l] = e;
    }
    diag[nv] = xpx[nstop] * xpx[nstop];
    free(tmp);
    return 0;
}

/* .......................................................... */

static int choleski (double *xpx, int nv)
     /* Choleski decomposition of X'X */
{
    int nm1, i, j, k, kk, l, jm1;
    double e, d, d1, test, xx;

    nm1 = nv - 1;
    e = 1/sqrt(xpx[1]);
    xpx[1] = e;
    for (i=2; i<=nv; i++) xpx[i] *= e;
    kk = nv + 1;

    for (j=2; j<=nv; j++) {
    /* diagonal elements */
        d = d1 = 0.0;
        k = j;
        jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d += xx * xx;
            k += nv-l;
        }
        test = xpx[kk] - d;
        if (test <= TINY) 
           return 1;
        e = 1/sqrt(test);
        xpx[kk] = e;
        /* off-diagonal elements */
        for (i=j+1; i<=nv; i++) {
            kk++;
            d = 0.0;
            k = j;
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

