/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
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

/* Native gretl code for ARMA estimation.  Much of the code here
   was contributed by Riccardo "Jack" Lucchetti, the rest is due
   to Allin Cottrell; thanks also to Stephen Moshier for cephes.
*/

#include "libgretl.h"
#include "bhhh_max.h"

#include "../cephes/polrt.c"

#define ARMA_DEBUG 0

#define MAX_ARMA_ORDER 6

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

struct arma_info {
    int yno;      /* ID of dependent variable */
    int p;        /* non-seasonal AR order */
    int q;        /* non-seasonal MA order */
    int P;        /* seasonal AR order */
    int Q;        /* seasonal MA order */
    int maxlag;   /* longest lag in model */
    int r;        /* number of other regressors (ARMAX) */
    int ifc;      /* 1 for intercept included, otherwise 0 */
    int nc;       /* total number of coefficients */
    int t1;       /* starting observation */
    int t2;       /* ending observation */
    int seasonal; /* 1 if any seasonal terms, otherwise 0 */
    int pd;       /* periodicity of data */
};

/* check whether the MA estimates have gone out of bounds in the
   course of BHHH iterations */

static int 
ma_out_of_bounds (struct arma_info *ainfo, const double *ma_coeff,
		  const double *sma_coeff)
{
    double *temp = NULL, *tmp2 = NULL;
    double re, im, rt;
    cmplx *roots = NULL;
    int qmax = ainfo->q;
    int i, j, k;
    int err = 0, allzero = 1;

    for (i=0; i<ainfo->q && allzero; i++) {
	if (ma_coeff[i] != 0.0) {
	    allzero = 0;
	}    
    }  

    for (i=0; i<ainfo->Q && allzero; i++) {
	if (sma_coeff[i] != 0.0) {
	    allzero = 0;
	}    
    }  
    
    if (allzero) {
	return 0;
    }

    if (ainfo->seasonal) {
	qmax += ainfo->Q * ainfo->pd;
    }    

    /* we'll use a budget version of the "arma_roots" function here */

    temp  = malloc((qmax + 1) * sizeof *temp);
    tmp2  = malloc((qmax + 1) * sizeof *tmp2);
    roots = malloc(qmax * sizeof *roots);

    if (temp == NULL || tmp2 == NULL || roots == NULL) {
	free(temp);
	free(tmp2);
	free(roots);
	return 1;
    }

    temp[0] = 1.0;

    /* Is this right, with seasonal MA?  Should we be 
       calculating distinct seasonal and non-seasonal
       roots?? 
    */
    

    /* initialize to non-seasonal MA or zero */
    for (i=0; i<qmax; i++) {
	if (i < ainfo->q) {
	    temp[i+1] = ma_coeff[i];
	} else {
	    temp[i+1] = 0.0;
	}
    }

    /* add seasonal MA and interaction */
    for (i=0; i<ainfo->Q; i++) {
	k = (i + 1) * ainfo->pd;
	temp[k] += sma_coeff[i];
	for (j=0; j<ainfo->q; j++) {
	    int m = k + j + 1;

	    temp[m] += sma_coeff[i] * ma_coeff[j];
	}
    }

    polrt(temp, tmp2, qmax, roots);

    for (i=0; i<qmax; i++) {
	re = roots[i].r;
	im = roots[i].i;
	rt = re * re + im * im;
	if (rt > DBL_EPSILON && rt <= 1.0) {
	    fprintf(stderr, "MA root %d = %g\n", i, rt);
	    err = 1;
	    break;
	}
    }

    free(temp);
    free(tmp2);
    free(roots);

    return err;
}

static void do_MA_partials (double *drv,
			    struct arma_info *ainfo,
			    const double *ma_coeff,
			    const double *sma_coeff,
			    int t)
{
    int i, j, s, p;

    for (i=0; i<ainfo->q; i++) {
	s = t - (i + 1);
	drv[t] -= ma_coeff[i] * drv[s];
    }

    for (i=0; i<ainfo->Q; i++) {
	s = t - ainfo->pd * (i + 1);
	drv[t] -= sma_coeff[i] * drv[s];
	for (j=0; j<ainfo->q; j++) {
	    p = s - (j + 1);
	    drv[t] -= sma_coeff[i] * ma_coeff[j] * drv[p];
	}
    }
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a "callback". */

static int arma_ll (double *coeff, 
		    const double **bhX, double **Z, 
		    model_info *arma,
		    int do_score)
{
    int i, j, s, t;
    int t1 = model_info_get_t1(arma);
    int t2 = model_info_get_t2(arma);
    int n = t2 - t1 + 1;

    const double *y = bhX[0];
    const double **X = bhX + 1;
    double **series = model_info_get_series(arma);
    double *e = series[0];
    double **de = series + 1;
    double **de_a, **de_sa, **de_m, **de_sm, **de_r;
    const double *ar_coeff, *sar_coeff;
    const double *ma_coeff, *sma_coeff;
    const double *reg_coeff;

    struct arma_info *ainfo;

    double ll, s2 = 0.0;
    int err = 0;

    /* retrieve ARMA-specific information */
    ainfo = model_info_get_extra_info(arma);

    /* pointers to blocks of coefficients */
    ar_coeff = coeff + ainfo->ifc;
    sar_coeff = ar_coeff + ainfo->p;
    ma_coeff = sar_coeff + ainfo->P;
    sma_coeff = ma_coeff + ainfo->q;
    reg_coeff = sma_coeff + ainfo->Q;

    /* pointers to blocks of derivatives */
    de_a = de + ainfo->ifc;
    de_sa = de_a + ainfo->p;
    de_m = de_sa + ainfo->P;
    de_sm = de_m + ainfo->q;
    de_r = de_sm + ainfo->Q;

#if ARMA_DEBUG
    fprintf(stderr, "arma_ll: p=%d, q=%d, P=%d, Q=%d, pd=%d\n",
	    ainfo->p, ainfo->q, ainfo->P, ainfo->Q, ainfo->pd);
#endif

    if (ma_out_of_bounds(ainfo, ma_coeff, sma_coeff)) {
	fputs("arma: MA estimate(s) out of bounds\n", stderr);
	return 1;
    }

    /* update forecast errors */

    for (t=t1; t<=t2; t++) {
	int p;

	e[t] = y[t];

	/* intercept */
	if (ainfo->ifc) {
	    e[t] -= coeff[0];
	} 

	/* non-seasonal AR component */
	for (i=0; i<ainfo->p; i++) {
	    s = t - (i + 1);
	    e[t] -= ar_coeff[i] * y[s];
	}

	/* seasonal AR component plus interactions */
	for (i=0; i<ainfo->P; i++) {
	    s = t - ainfo->pd * (i + 1);
	    e[t] -= sar_coeff[i] * y[s];
	    for (j=0; j<ainfo->p; j++) {
		p = s - (j + 1);
		e[t] -= sar_coeff[i] * ar_coeff[j] * y[p];
	    }
	}

	/* non-seasonal MA component */
	for (i=0; i<ainfo->q; i++) {
	    s = t - (i + 1);
	    if (s >= t1) {
		e[t] -= ma_coeff[i] * e[s];
	    }
	}

	/* seasonal MA component plus interactions */
	for (i=0; i<ainfo->Q; i++) {
	    s = t - ainfo->pd * (i + 1);
	    if (s >= t1) {
		e[t] -= sma_coeff[i] * e[s];
		for (j=0; j<ainfo->q; j++) {
		    p = s - (j + 1);
		    if (p >= t1) {
			e[t] -= sma_coeff[i] * ma_coeff[j] * e[p];
		    }
		}
	    }
	}

	/* exogenous regressors */
	for (i=0; i<ainfo->r; i++) {
	    e[t] -= reg_coeff[i] * X[i][t];
	}

	s2 += e[t] * e[t];
    }

    /* get error variance and log-likelihood */

    s2 /= (double) n;

    ll = -n * (0.5 * log(s2) + LN_SQRT_2_PI_P5);
    model_info_set_ll(arma, ll, do_score);

    if (do_score) {
	int lag, xlag;
	double x;

	for (t=t1; t<=t2; t++) {

	    /* the constant term (de_0) */
	    if (ainfo->ifc) {
		de[0][t] = -1.0;
		do_MA_partials(de[0], ainfo, ma_coeff, sma_coeff, t);
	    }

	    /* non-seasonal AR terms (de_a) */
	    for (j=0; j<ainfo->p; j++) {
		lag = j + 1;
		if (t >= lag) {
		    de_a[j][t] = -y[t-lag];
		    /* cross-partial with seasonal AR */
		    for (i=0; i<ainfo->P; i++) {
			xlag = lag + ainfo->pd * (i + 1);
			if (t >= xlag) {
			    de_a[j][t] -= sar_coeff[i] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_a[j], ainfo, ma_coeff, sma_coeff, t);
		}
	    }

	    /* seasonal AR terms (de_sa) */
	    for (j=0; j<ainfo->P; j++) {
		lag = ainfo->pd * (j + 1);
		if (t >= lag) {
		    de_sa[j][t] = -y[t-lag];
		    /* cross-partial with non-seasonal AR */
		    for (i=0; i<ainfo->p; i++) {
			xlag = lag + (i + 1);
			if (t >= xlag) {
			    de_sa[j][t] -= ar_coeff[i] * y[t-xlag];
			}
		    }
		    do_MA_partials(de_sa[j], ainfo, ma_coeff, sma_coeff, t);
		}
	    }

	    /* non-seasonal MA terms (de_m) */
	    for (j=0; j<ainfo->q; j++) {
		lag = j + 1;
		if (t >= lag) {
		    de_m[j][t] = -e[t-lag];
		    /* cross-partial with seasonal MA */
		    for (i=0; i<ainfo->Q; i++) {
			xlag = lag + ainfo->pd * (i + 1);
			if (t >= xlag) {
			    de_m[j][t] -= sma_coeff[i] * e[t-xlag];
			}
		    }
		    do_MA_partials(de_m[j], ainfo, ma_coeff, sma_coeff, t);
		}
	    }

	    /* seasonal MA terms (de_sm) */
	    for (j=0; j<ainfo->Q; j++) {
		lag = ainfo->pd * (j + 1);
		if (t >= lag) {
		    de_sm[j][t] = -e[t-lag];
		    /* cross-partial with non-seasonal MA */
		    for (i=0; i<ainfo->q; i++) {
			xlag = lag + (i + 1);
			if (t >= xlag) {
			    de_sm[j][t] -= ma_coeff[i] * e[t-xlag];
			}
		    }
		    do_MA_partials(de_sm[j], ainfo, ma_coeff, sma_coeff, t);
		}
	    }

	    /* exogenous regressors (de_r) */
	    for (j=0; j<ainfo->r; j++) {
		de_r[j][t] = -X[j][t]; 
		do_MA_partials(de_r[j], ainfo, ma_coeff, sma_coeff, t);
	    }

	    /* update OPG data set */
	    x = e[t] / s2; /* sqrt(s2)? does it matter? */
	    for (i=0; i<ainfo->nc; i++) {
		Z[i+1][t] = -de[i][t] * x;
	    }
	}
    }

    if (isnan(ll)) {
	err = 1;
    }

    return err;
}

/*
  Given an ARMA process $A(L)B(L) y_t = C(L)D(L) \epsilon_t$, returns the 
  roots of the four polynomials -- or just two polynomials if seasonal
  AR and MA effects, B(L) and D(L) are not present.

  ainfo: gives various pieces of information on the ARMA model,
  including seasonal and non-seasonal AR and MA orders.

  coeff: p+q+P+Q+ifc vector of coefficients (if an intercept is present
  it is element 0 and is ignored)

  returns: the p + P + q + Q roots (AR part first)
*/

static cmplx *arma_roots (struct arma_info *ainfo, const double *coeff) 
{
    const double *ar_coeff = coeff + ainfo->ifc;
    const double *sar_coeff = ar_coeff + ainfo->p;
    const double *ma_coeff = sar_coeff + ainfo->P;
    const double *sma_coeff = ma_coeff + ainfo->q;

    int nr = ainfo->p + ainfo->P + ainfo->q + ainfo->Q;
    int pmax, qmax, lmax;
    double *temp = NULL, *temp2 = NULL;
    cmplx *rptr, *roots = NULL;
    int i;

    pmax = (ainfo->p > ainfo->P)? ainfo->p : ainfo->P;
    qmax = (ainfo->q > ainfo->Q)? ainfo->q : ainfo->Q;
    lmax = (pmax > qmax)? pmax : qmax;

    temp  = malloc((lmax + 1) * sizeof *temp);
    temp2 = malloc((lmax + 1) * sizeof *temp2);
    roots = malloc(nr * sizeof *roots);

    if (temp == NULL || temp2 == NULL || roots == NULL) {
	free(temp);
	free(temp2);
	free(roots);
	return NULL;
    }

    temp[0] = 1.0;
    rptr = roots;

    if (ainfo->p > 0) {
	/* A(L), non-seasonal */
	for (i=0; i<ainfo->p; i++) {
	    temp[i+1] = -ar_coeff[i];
	}
	polrt(temp, temp2, ainfo->p, rptr);
	rptr += ainfo->p;
    }

    if (ainfo->P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo->P; i++) {
	    temp[i+1] = -sar_coeff[i];
	}    
	polrt(temp, temp2, ainfo->P, rptr);
	rptr += ainfo->P;
    }

    if (ainfo->q > 0) {
	/* C(L), non-seasonal */
	for (i=0; i<ainfo->q; i++) {
	    temp[i+1] = ma_coeff[i];
	}  
	polrt(temp, temp2, ainfo->q, rptr);
	rptr += ainfo->q;
    }

    if (ainfo->Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo->Q; i++) {
	    temp[i+1] = sma_coeff[i];
	}  
	polrt(temp, temp2, ainfo->Q, rptr);
    }
    
    free(temp);
    free(temp2);

    return roots;
}

/* package the various statistics from ARMA estimation in the
   form of a gretl MODEL struct */

static void rewrite_arma_model_stats (MODEL *pmod, model_info *arma,
				      const int *list, const double *y, 
				      const double *theta, 
				      struct arma_info *ainfo)
{
    int i, t;
    double **series = model_info_get_series(arma);
    const double *e = series[0];
    double mean_error;

    pmod->ci = ARMA;
    pmod->ifc = ainfo->ifc;

    pmod->lnL = model_info_get_ll(arma);

    pmod->dfn = ainfo->nc - ainfo->ifc;
    pmod->dfd = pmod->nobs - pmod->dfn;
    pmod->ncoeff = ainfo->nc;

    for (i=0; i<pmod->ncoeff; i++) {
	pmod->coeff[i] = theta[i];
    }

    free(pmod->list);
    pmod->list = gretl_list_copy(list);

    pmod->ybar = gretl_mean(pmod->t1, pmod->t2, y);
    pmod->sdy = gretl_stddev(pmod->t1, pmod->t2, y);

    mean_error = pmod->ess = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->uhat[t] = e[t];
	pmod->yhat[t] = y[t] - pmod->uhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
	mean_error += pmod->uhat[t];
    }

    mean_error /= pmod->nobs;
    gretl_model_set_double(pmod, "mean_error", mean_error);

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->tss += (y[t] - pmod->ybar) * (y[t] - pmod->ybar);
    }

    if (pmod->tss > pmod->ess) {
	pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);
    } else {
	pmod->fstt = NADBL;
    }

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    mle_aic_bic(pmod, 1);

    if (ainfo->seasonal) {
	gretl_model_set_int(pmod, "arma_P", ainfo->P);
	gretl_model_set_int(pmod, "arma_Q", ainfo->Q);
	gretl_model_set_int(pmod, "arma_pd", ainfo->pd);	
    }

    if (ainfo->r > 0) {
	gretl_model_set_int(pmod, "armax", 1);
    }
}

/* remove the intercept from list of regressors */

static int remove_const (int *list, int seasonal)
{
    int xstart, ret = 0;
    int i, j;

    xstart = (seasonal)? 8 : 5;

    for (i=xstart; i<=list[0]; i++) {
	if (list[i] == 0) {
	    for (j=i; j<list[0]; j++) {
		list[j] = list[j+1];
	    }
	    list[0] -= 1;
	    ret = 1;
	    break;
	}
    }

    return ret;
}

#define has_seasonals(l) (l[0] > 5 && l[3] == LISTSEP && l[6] == LISTSEP)

static int 
check_arma_list (int *list, gretlopt opt, struct arma_info *ainfo)
{
    int armax = 0;
    int hadconst = 0;
    int err = 0;

    /* FIXME need more checks for seasonal case, e.g. non-
       seasonal order must be less than pd */

    ainfo->seasonal = has_seasonals(list);
    ainfo->p = ainfo->q = 0;
    ainfo->P = ainfo->Q = 0;
    ainfo->r = ainfo->nc = 0;

    if (ainfo->seasonal) {
	armax = (list[0] > 7);
    } else {
	armax = (list[0] > 4);
    }

    if (list[0] < 4) {
	err = 1;
    } else if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[1] + list[2] == 0 && !ainfo->seasonal) {
	err = 1;
    }

    if (!err) {
	ainfo->p = list[1];
	ainfo->q = list[2];
    }

    if (!err && ainfo->seasonal) {
	if (list[0] < 7) {
	    err = 1;
	} else if (list[4] < 0 || list[4] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[4] + list[5] == 0) {
	    err = 1;
	}
    }

    if (!err && ainfo->seasonal) {
	ainfo->P = list[4];
	ainfo->Q = list[5];
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_S (meaning: suppress the intercept) we'll flag this by
       setting ifc = 0.  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.
    */

    if (!err) {
	if (armax) {
	    hadconst = remove_const(list, ainfo->seasonal);
	}
	if ((opt & OPT_S) || (armax && !hadconst)) {
	    ainfo->ifc = 0;
	} else {
	    ainfo->ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command"));
    } else {
	ainfo->r = list[0] - ((ainfo->seasonal)? 7 : 4);
	ainfo->nc = ainfo->p + ainfo->q + ainfo->P + ainfo->Q
	    + ainfo->r + ainfo->ifc;
    }

    return err;
}

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

static const double **
make_armax_X (int *list, struct arma_info *ainfo, const double **Z)
{
    const double **X;
    int offset, nv;
    int v, i;

    offset = (ainfo->seasonal)? 7 : 4;
    nv = list[0] - offset;

    X = malloc((nv + 1) * sizeof *X);
    if (X == NULL) {
	return NULL;
    }

    /* the dependent variable */
    X[0] = Z[list[offset]];

    /* the independent variables */
    for (i=1; i<=nv; i++) {
	v = list[i + offset];
	X[i] = Z[v];
    }

    return X;
}

/* Run an initial OLS to get initial values for the AR coefficients */

static int ar_init_by_ols (const int *list, double *coeff,
			   const double **Z, const DATAINFO *pdinfo,
			   struct arma_info *ainfo)
{
    int an = pdinfo->t2 - ainfo->t1 + 1;
    int av = ainfo->p + ainfo->P + ainfo->r + 2;
    int np = ainfo->p, nq = ainfo->q;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int offset, xstart;
    int i, j, t, err = 0;

    gretl_model_init(&armod);  

    alist = gretl_list_new(av);
    if (alist == NULL) {
	return 1;
    }

    alist[1] = 1;

    if (ainfo->ifc) {
	alist[2] = 0;
	offset = 3;
    } else {
	alist[0] -= 1;
	offset = 2;
    }

    for (i=0; i<ainfo->p; i++) {
	alist[i + offset] = i + 2;
    }

    for (i=0; i<ainfo->P; i++) {
	alist[i + offset + ainfo->p] = i + ainfo->p + 2;
   }

    for (i=0; i<ainfo->r; i++) {
	alist[i + offset + ainfo->p + ainfo->P] = i + ainfo->p + ainfo->P + 2;
    }

    adinfo = create_new_dataset(&aZ, av, an, 0);
    if (adinfo == NULL) {
	free(alist);
	return 1;
    }

    xstart = (ainfo->seasonal)? 8 : 5;

    /* build temporary dataset including lagged vars */
    for (t=0; t<an; t++) {
	int k, s;

	s = t + ainfo->t1;
	aZ[1][t] = Z[ainfo->yno][s];

	for (i=0; i<ainfo->p; i++) {
	    s = t + ainfo->t1 - (i + 1);
	    k = i + 2;
	    aZ[k][t] = Z[ainfo->yno][s];
	}

	for (i=0; i<ainfo->P; i++) {
	    s = t + ainfo->t1 - ainfo->pd * (i + 1);
	    k = i + ainfo->p + 2;
	    aZ[k][t] = Z[ainfo->yno][s];
	}

	s = t + ainfo->t1;

	for (i=0; i<ainfo->r; i++) {
	    k = i + ainfo->p + ainfo->P + 2;
	    j = list[xstart + i];
	    aZ[k][t] = Z[j][s];
	}
    }

    if (ainfo->seasonal) {
	np += ainfo->P;
	nq += ainfo->Q;
    }

    /* run the OLS */
    armod = lsq(alist, &aZ, adinfo, OLS, OPT_A | OPT_Z, 0.0);
    err = armod.errcode;
    if (!err) {
	j = 0;
	for (i=0; i<armod.ncoeff; i++) {
	    if (i == np + ainfo->ifc) {
		j += nq; /* reserve space for MA coeffs */
	    }
	    coeff[j++] = armod.coeff[i];
	}
	for (i=0; i<nq; i++) {
	    /* insert zeros for MA coeffs */
	    coeff[i + np + ainfo->ifc] = 0.0;
	} 
    }

#if ARMA_DEBUG
    fprintf(stderr, "OLS init: ncoeff = %d, nobs = %d\n", 
	    armod.ncoeff, armod.nobs);
    for (i=0; i<armod.ncoeff; i++) {
	fprintf(stderr, " coeff[%d] = %g\n", i, armod.coeff[i]);
    }
#endif

    /* clear everything up */
    free(alist);
    free_Z(aZ, adinfo);
    clear_datainfo(adinfo, CLEAR_FULL);
    free(adinfo);

    clear_model(&armod);

    return err;
}

static int 
adjust_sample (DATAINFO *pdinfo, const double **Z, const int *list,
	       struct arma_info *ainfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int an, i, v, t, t1min = 0;
    int vstart, pmax, anymiss;

    vstart = (ainfo->seasonal)? 7 : 4;

    for (t=0; t<=pdinfo->t2; t++) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) {
	    t1min++;
        } else {
	    break;
	}
    }

    t1min += ainfo->maxlag;

    if (t1 < t1min) {
	t1 = t1min;
    }

    for (t=pdinfo->t2; t>=t1; t--) {
	anymiss = 0;
	for (i=vstart; i<=list[0]; i++) {
	    v = list[i];
	    if (na(Z[v][t])) {
		anymiss = 1;
		break;
	    }
	}
	if (anymiss) {
	    t2--;
        } else {
	    break;
	}
    }

    pmax = ainfo->p;
    if (ainfo->P > 0) {
	pmax += ainfo->P * ainfo->pd;
    }

    for (t=t1-pmax; t<t2; t++) {
	for (i=vstart; i<=list[0]; i++) {
	    if (t < t1 && i > vstart) {
		continue;
	    }
	    v = list[i];
	    if (na(Z[v][t])) {
		char msg[64];

		sprintf(msg, _("Missing value encountered for "
			       "variable %d, obs %d"), v, t + 1);
		gretl_errmsg_set(msg);
		return 1;
	    }
	}
    }

    an = t2 - t1 + 1;
    if (an <= ainfo->nc) {
	return 1; 
    }

    ainfo->t1 = t1;
    ainfo->t2 = t2;

    return 0;
}

/* set up a model_info struct for passing to bhhh_max */

static model_info *
set_up_arma_model_info (struct arma_info *ainfo)
{
    model_info *arma;

    arma = model_info_new(ainfo->nc, ainfo->t1, ainfo->t2, 1.0e-6);

    if (arma == NULL) return NULL;

    model_info_set_opts(arma, PRESERVE_OPG_MODEL);
    model_info_set_n_series(arma, ainfo->nc + 1);

    /* add pointer to ARMA-specific details */
    model_info_set_extra_info(arma, ainfo);

    return arma;
}

static void calc_max_lag (struct arma_info *ainfo)
{
    int pmax = ainfo->p;
    int qmax = ainfo->q;

    if (ainfo->seasonal) {
	pmax += ainfo->P * ainfo->pd;
	qmax += ainfo->Q * ainfo->pd;
    }

    ainfo->maxlag = (pmax > qmax)? pmax : qmax;
}

MODEL arma_model (const int *list, const double **Z, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn)
{
    double *coeff = NULL;
    const double **X = NULL;
    int *alist = NULL;
    PRN *aprn = NULL;
    model_info *arma = NULL;
    MODEL armod;
    struct arma_info ainfo;
    int err = 0;

    if (opt & OPT_V) {
	aprn = prn;
    } 

    gretl_model_init(&armod); 
    gretl_model_smpl_init(&armod, pdinfo);

    alist = gretl_list_copy(list);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    if (check_arma_list(alist, opt, &ainfo)) {
	armod.errcode = E_UNSPEC;
	goto bailout;
    }

    /* dependent variable */
    ainfo.yno = (ainfo.seasonal)? alist[7] : alist[4];

    /* periodicity of data */
    ainfo.pd = pdinfo->pd;

    /* calculate maximum lag */
    calc_max_lag(&ainfo);

    /* adjust sample range if need be */
    if (adjust_sample(pdinfo, Z, alist, &ainfo)) {
        armod.errcode = E_DATA;
	goto bailout;
    }

    /* allocate initial coefficient vector */
    coeff = malloc(ainfo.nc * sizeof *coeff);
    if (coeff == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    /* create model_info struct to feed to bhhh_max() */
    arma = set_up_arma_model_info(&ainfo);
    if (arma == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    /* initialize the coefficients: AR and regression part by OLS, 
       MA at 0 */
    err = ar_init_by_ols(alist, coeff, Z, pdinfo, &ainfo);
    if (err) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }	

    /* construct virtual dataset for dep var, real regressors */
    X = make_armax_X(alist, &ainfo, Z);
    if (X == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    /* call BHHH conditional ML function (OPG regression) */
    err = bhhh_max(arma_ll, X, coeff, arma, aprn);

    if (err) {
	fprintf(stderr, "arma: bhhh_max returned %d\n", err);
	armod.errcode = E_NOCONV;
    } else {
	MODEL *pmod = model_info_capture_OPG_model(arma);
	double *theta = model_info_get_theta(arma);
	cmplx *roots;

	rewrite_arma_model_stats(pmod, arma, alist, Z[ainfo.yno], 
				 theta, &ainfo);

	/* compute and save polynomial roots */
	roots = arma_roots(&ainfo, theta);
	if (roots != NULL) {
	    gretl_model_set_data(pmod, "roots", roots,
				 (ainfo.p + ainfo.q) * sizeof *roots);
	}

	gretl_model_add_arma_varnames(pmod, pdinfo, ainfo.yno, ainfo.p,
				      ainfo.q, ainfo.P, ainfo.Q, ainfo.r);

	armod = *pmod;
	free(pmod);
    }

 bailout:

    free(alist);
    free(coeff);
    free(X);
    model_info_free(arma);

    return armod;
}
