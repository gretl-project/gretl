/* 
 * Copyright (C) 2005 Riccardo "Jack" Lucchetti
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


#include "libgretl.h"
#include "gretl_private.h"
#include "../cephes/libprob.h"

#undef PDEBUG 

#define POISSON_TOL 1.0e-10 
#define POISSON_MAX_ITER 100 

/* check whether a series contains nothing but non-negative
   integer values (some of which are > 1) */

static int is_count_variable (const double *x, int t1, int t2)
{
    int t, xi;
    int g1 = 0;
    int ret = 1;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    continue;
	}
	if (x[t] < 0.0) {
	    ret = 0;
	    break;
	}
	xi = x[t];
	if (x[t] != (double) xi) {
	    ret = 0;
	    break;
	}
	if (x[t] > 1.0) {
	    g1 = 1;
	}
    }

    if (g1 == 0) {
	ret = 0;
    }

    return ret;
}

static double poisson_ll (const double *y, const double *mu, 
			  int t1, int t2)
{
    int t;
    double loglik = 0.0;
    double mt, lmt, yt, lytfact, llt;

    for (t=t1; t<=t2; t++) {
	if (y[t] != NADBL) {
	    mt = mu[t];
	    lmt = log(mt);
	    yt = y[t];
	    lytfact = log(cephes_gamma(1.0 + yt));
	    llt = (-mt + yt*lmt - lytfact);
	    loglik += llt;
	}
    }  

    return loglik;
}

/* make covariance matrix based on the 'artificial' regression, and
   transcribe the relevant elements into the target model */

static int make_poisson_vcv (MODEL *targ, MODEL *src)
{
    int nc = targ->ncoeff;
    int nt = (nc * nc + nc) / 2;
    int tidx, sidx;
    int i, j;
    int err = 0;

    if (makevcv(src)) {
	err = 1;
    }

    if (!err && targ->vcv == NULL) {
	targ->vcv = malloc(nt * sizeof *targ->vcv);
	if (targ->vcv == NULL) {
	    err = 1;
	}
    }

    if (!err) {
	for (i=0; i<targ->ncoeff; i++) {
	    for (j=i; j<targ->ncoeff; j++) {
		tidx = ijton(i, j, targ->ncoeff);
		sidx = ijton(i, j, src->ncoeff);
		targ->vcv[tidx] = src->vcv[sidx];
	    }
	}
    }

    return err;
}

#undef POISSON_USE_GNR

#ifdef POISSON_USE_GNR

/* try the GNR-type approach in Davidson and MacKinnon, ETM, ch 11 */

static int 
transcribe_poisson_results (MODEL *targ, MODEL *src,
			    const double **Z, int yno, int munum,
			    int iter)
{
    int i, t;
    int err = 0;
    
    gretl_model_set_int(targ, "iters", iter);

    targ->ci = POISSON;

    targ->ess = 0.0;

    for (t=targ->t1; t<=targ->t2; t++) {
	targ->yhat[t] = Z[munum][t];
	targ->uhat[t] = Z[yno][t] - Z[munum][t];
	targ->ess += targ->uhat[t] * targ->uhat[t];
    }

    for (i=0; i<targ->ncoeff; i++) {
	targ->sderr[i] = src->sderr[i];
    }

    targ->lnL = poisson_ll(Z[yno], Z[munum], targ->t1, targ->t2);

#ifdef PDEBUG
    fprintf(stderr, "log-likelihood = %g\n", targ->lnL);
#endif

    mle_aic_bic(targ, 0); 

    /* mask invalid statistics */
    targ->rsq = NADBL;
    targ->adjrsq = NADBL;
    targ->sigma = NADBL;
    targ->fstt = NADBL;

    err = make_poisson_vcv(targ, src);

    return err;
}

static int 
do_poisson (MODEL *pmod, int offvar, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int i, t;
    int origv = pdinfo->v;
    int yno = pmod->list[1];
    int wtnum, depnum, munum;

    int iter = 0;
    double crit = 1.0;
    double xb;

    MODEL tmpmod;
    int *local_list = NULL;

    gretl_model_init(&tmpmod);

    local_list = gretl_list_new(pmod->list[0] + 1);
    if (local_list == NULL) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    if (dataset_add_vars(3, pZ, pdinfo)) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    wtnum = origv;
    depnum = origv + 1;
    munum = origv + 2;

    local_list[0] = pmod->list[0] + 1;
    local_list[1] = wtnum;
    local_list[2] = depnum;

    for (i=3; i<=local_list[0]; i++) { 
	/* original independent vars */
	local_list[i] = pmod->list[i-1];
    }    

    pmod->coeff[0] = log(pmod->ybar);
    for (i=1; i<pmod->ncoeff; i++) { 
	pmod->coeff[i] = 0.0;
    }

    xb = pmod->coeff[0];

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (offvar > 0) {
	    xb = pmod->coeff[0] + log((*pZ)[offvar][t]);
	}
	(*pZ)[depnum][t] = exp(-xb) * ((*pZ)[yno][t] - exp(xb));
	(*pZ)[wtnum][t] = exp(.5 * xb);
    }

    pputc(prn, '\n');

    while (iter < POISSON_MAX_ITER && crit > POISSON_TOL) {

	iter++;

	tmpmod = lsq(local_list, pZ, pdinfo, WLS, OPT_A | OPT_M, 0.0);

	if (tmpmod.errcode) {
	    fprintf(stderr, "poisson_estimate: lsq returned %d\n", 
		    tmpmod.errcode);
	    pmod->errcode = tmpmod.errcode;
	    break;
	}

	crit = tmpmod.nobs * tmpmod.rsq;

	pprintf(prn, "iter = %d\tcrit = %g\n", iter, crit);

	for (i=0; i<tmpmod.ncoeff; i++) { 
	    pmod->coeff[i] += tmpmod.coeff[i];
	    tmpmod.coeff[i] = pmod->coeff[i];
#ifdef PDEBUG
	    fprintf(stderr, "coeff[%d] = %g,\tgrad[%d] = %g\n", i, coeff[i], 
		    i, tmpmod.coeff[i]);
#endif
	}

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    xb = 0.0;
	    for (i=0; i<tmpmod.ncoeff; i++) { 
		int vi = pmod->list[i+2];

		xb += pmod->coeff[i] * (*pZ)[vi][t];
	    }
	    if (offvar > 0) {
		xb += log((*pZ)[offvar][t]);
	    }
	    (*pZ)[depnum][t] = exp(-xb) * ((*pZ)[yno][t] - exp(xb));
	    (*pZ)[wtnum][t] = exp(.5 * xb);
	    (*pZ)[munum][t] = exp(xb);
	}

	if (crit > POISSON_TOL) {
	    clear_model(&tmpmod);
	}
    }

    pputc(prn, '\n');

    if (crit > POISSON_TOL) {
	pmod->errcode = E_NOCONV;
    } 

    if (pmod->errcode == 0) {
	transcribe_poisson_results(pmod, &tmpmod, (const double **) *pZ, 
				   yno, munum, iter);
    }

 bailout:

    clear_model(&tmpmod);
    free(local_list);
    dataset_drop_vars(pdinfo->v - origv, pZ, pdinfo);

    return pmod->errcode;
}

#else

static int 
transcribe_poisson_results (MODEL *targ, MODEL *src,
			    const double **X, int nvars,
			    int n, int iter)
{
    double mt;
    int i, t, s;
    int err = 0;
    
    gretl_model_set_int(targ, "iters", iter);

    targ->ci = POISSON;

    targ->ess = 0.0;

    for (t=0; t<n; t++) {
	s = t + targ->t1;
	mt = X[nvars][t];
	targ->yhat[s] = mt;
	targ->uhat[s] = X[0][t] - mt;
	targ->ess += targ->uhat[s] * targ->uhat[s];
    }

    for (i=0; i<targ->ncoeff; i++) {
	targ->sderr[i] = src->sderr[i];
    }

    targ->lnL = poisson_ll(X[0], X[nvars], 0, n - 1);

#ifdef PDEBUG
    fprintf(stderr, "log-likelihood = %g\n", targ->lnL);
#endif

    mle_aic_bic(targ, 0); /* is the number of params right? */

    /* mask invalid statistics */
    targ->rsq = NADBL;
    targ->adjrsq = NADBL;
    targ->sigma = NADBL;
    targ->fstt = NADBL;

    err = make_poisson_vcv(targ, src);

    return err;
}

static int do_poisson (MODEL *pmod, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    int i, t, s;
    int iter = 0;
    double crit = 1.0;

    MODEL tmpmod;
    double **X = NULL;
    DATAINFO *tmpinfo = NULL;
    int *local_list = NULL;

    int nvars = pmod->list[0];
    int ncolX = nvars + 3;
    int ncoeff = nvars - 1;
    int n = pmod->nobs;
    
    double mt, et, wt;

#ifdef PDEBUG
    for (i=0; i<=pmod->list[0]; i++) { /* original vars */
	fprintf(stderr, "list[%d] = %d\n", i, pmod->list[i]);
    }
#endif

    gretl_model_init(&tmpmod);

    /* last 3 columns: mu_t, stdresid_t, wgt_t */
    
    local_list = gretl_list_new(nvars + 1);
    if (local_list == NULL) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

    local_list[0] = nvars + 1; /* original vars + weights */

    for (i=3; i<=local_list[0]; i++) { /* original vars */
	local_list[i] = i - 2;
    }

    local_list[1] = ncolX - 1; /* weights */
    local_list[2] = ncolX - 2; /* pseudo-residuals */

    /* note: the dataset X may be shorter than the original
       dataset, Z */

    tmpinfo = create_new_dataset(&X, ncolX, n, 0);
    if (tmpinfo == NULL) {
	pmod->errcode = E_ALLOC;
	goto bailout;
    }

#ifdef PDEBUG
    fprintf(stderr, "poisson_estimate: ybar = %g\n", pmod->ybar);
#endif

    mt = pmod->ybar;
    wt = sqrt(mt);

    pmod->coeff[0] = log(mt);
    for (i=1; i<ncoeff; i++) { 
	pmod->coeff[i] = 0.0;
    }

    for (t=0; t<n; t++) {
	s = t + pmod->t1;
	for (i=0; i<nvars; i++) { /* dep var + regressors */
	    X[i][t] = Z[pmod->list[i+1]][s];
	}
	et = X[0][t] / mt - 1.0;
	X[nvars][t] = mt;
	X[nvars + 1][t] = et;
	X[nvars + 2][t] = wt;
    }

    pputc(prn, '\n');

    while (iter < POISSON_MAX_ITER && crit > POISSON_TOL) {

	iter++;

	tmpmod = lsq(local_list, &X, tmpinfo, WLS, OPT_A | OPT_M, 0.0);

	if (tmpmod.errcode) {
	    fprintf(stderr, "poisson_estimate: lsq returned %d\n", 
		    tmpmod.errcode);
	    pmod->errcode = tmpmod.errcode;
	    break;
	}

	crit = n * tmpmod.rsq;

	pprintf(prn, "iter = %d\tcrit = %g\n", iter, crit);

	for (i=0; i<tmpmod.ncoeff; i++) { 
	    pmod->coeff[i] += tmpmod.coeff[i];
	    tmpmod.coeff[i] = pmod->coeff[i];
#ifdef PDEBUG
	    fprintf(stderr, "coeff[%d] = %g,\tgrad[%d] = %g\n", i, pmod->coeff[i], 
		    i, tmpmod.coeff[i]);
#endif
	}

	for (t=0; t<n; t++) {
	    X[nvars][t] *= exp(tmpmod.yhat[t]);
	    mt = X[nvars][t];
	    et = X[0][t] / mt - 1;
	    X[nvars + 1][t] = et;
	    X[nvars + 2][t] = sqrt(mt);
	}

	if (crit > POISSON_TOL) {
	    /* going round again: don't leak memory */
	    clear_model(&tmpmod);
	}
    }

    pputc(prn, '\n');

    if (crit > POISSON_TOL) {
	pmod->errcode = E_NOCONV;
    } 

    if (pmod->errcode == 0) {
	transcribe_poisson_results(pmod, &tmpmod, (const double **) X, 
				   nvars, n, iter);
    }

 bailout:

    clear_model(&tmpmod);

    free(local_list);		 
    free_Z(X, tmpinfo);

    clear_datainfo(tmpinfo, CLEAR_FULL);
    free(tmpinfo); 

    return pmod->errcode;
}

#endif

/* "offvar" = ID number of offset variable, i.e. a variable whose log
   should enter the equation with a coefficient of 1 (scaling).  This
   seems to be a fairly common thing.  I'm not sure how to implement it
   with the Amemiya approach, though I'm probably missing something
   simple.  AC
*/

int 
poisson_estimate (MODEL *pmod, int offvar, double ***pZ, DATAINFO *pdinfo,
		  PRN *prn) 
{
    int err = 0;

    if (!is_count_variable((*pZ)[pmod->list[1]], pmod->t1, pmod->t2)) {
	gretl_errmsg_set(_("poisson: the dependent variable must be count data"));
	err = pmod->errcode = E_DATA;
    } else {
#ifdef POISSON_USE_GNR
	err = do_poisson(pmod, offvar, pZ, pdinfo, prn);
#else
	err = do_poisson(pmod, (const double **) *pZ, pdinfo, prn);
#endif
    }

    return err;
}

