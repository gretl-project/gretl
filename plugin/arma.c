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
#include "libset.h"

#include "../cephes/polrt.c"

#define ARMA_DEBUG 0

/* ln(sqrt(2*pi)) + 0.5 */
#define LN_SQRT_2_PI_P5 1.41893853320467274178

#include "arma_common.c"

static PRN *errprn;

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

    if (arma_has_seasonal(ainfo)) {
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
	    pprintf(errprn, "MA root %d = %g\n", i, rt);
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
    ar_coeff = coeff + arma_has_const(ainfo);
    sar_coeff = ar_coeff + ainfo->p;
    ma_coeff = sar_coeff + ainfo->P;
    sma_coeff = ma_coeff + ainfo->q;
    reg_coeff = sma_coeff + ainfo->Q;

    /* pointers to blocks of derivatives */
    de_a = de + arma_has_const(ainfo);
    de_sa = de_a + ainfo->p;
    de_m = de_sa + ainfo->P;
    de_sm = de_m + ainfo->q;
    de_r = de_sm + ainfo->Q;

#if ARMA_DEBUG
    fprintf(stderr, "arma_ll: p=%d, q=%d, P=%d, Q=%d, pd=%d\n",
	    ainfo->p, ainfo->q, ainfo->P, ainfo->Q, ainfo->pd);
#endif

    if (ma_out_of_bounds(ainfo, ma_coeff, sma_coeff)) {
	pputs(errprn, "arma: MA estimate(s) out of bounds\n");
	fputs("arma: MA estimate(s) out of bounds\n", stderr);
	return 1;
    }

    /* update forecast errors */

    for (t=t1; t<=t2; t++) {
	int p;

	e[t] = y[t];

	/* intercept */
	if (arma_has_const(ainfo)) {
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
		e[t] += sar_coeff[i] * ar_coeff[j] * y[p];
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
	    if (arma_has_const(ainfo)) {
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
			    de_a[j][t] += sar_coeff[i] * y[t-xlag];
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
			    de_sa[j][t] += ar_coeff[i] * y[t-xlag];
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
    const double *ar_coeff = coeff + arma_has_const(ainfo);
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

    if (pmax == 0 && qmax == 0) {
	return NULL;
    }

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

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

static const double **
make_armax_X (int *list, struct arma_info *ainfo, const double **Z)
{
    const double **X;
    int ypos, nx;
    int v, i;

    ypos = arma_list_y_position(ainfo);
    nx = list[0] - ypos;

#if ARMA_DEBUG
    fprintf(stderr, "make_armax_X: allocating %d series pointers\n",
	    nx + 1);
#endif    

    X = malloc((nx + 1) * sizeof *X);
    if (X == NULL) {
	return NULL;
    }

    /* the dependent variable */
    if (ainfo->dy != NULL) {
	X[0] = ainfo->dy;
    } else {
	X[0] = Z[list[ypos]];
    }

    /* the independent variables */
    for (i=1; i<=nx; i++) {
	v = list[i + ypos];
	X[i] = Z[v];
    }

    return X;
}

static int arma_get_nls_model (MODEL *amod, struct arma_info *ainfo,
			       const int *alist, int axstart,
			       double ***pZ, DATAINFO *pdinfo) 
{
#if ARMA_DEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
    gretlopt opt = OPT_V;
#else
    PRN *prn = NULL;
    gretlopt opt = OPT_NONE;
#endif
    char fnstr[MAXLINE]; /* FIXME? */
    char term[32];
    nls_spec *spec;
    int *plist = NULL;
    int v = pdinfo->v;
    int nparam;
    int i, j, k, err = 0;

    spec = nls_spec_new(NLS, pdinfo);
    if (spec == NULL) {
	return E_ALLOC;
    }

    nparam = arma_has_const(ainfo) + ainfo->p + ainfo->P + ainfo->r;

    plist = gretl_list_new(nparam);
    if (plist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = dataset_add_scalars(nparam, pZ, pdinfo); 
    if (err) {
	goto bailout;
    }

    strcpy(fnstr, "y=");
    k = 1;

    if (arma_has_const(ainfo)) {
	double ybar = gretl_mean(0, pdinfo->n - 1, (*pZ)[1]);

	strcat(fnstr, "b0+");
	strcpy(pdinfo->varname[v], "b0");
	(*pZ)[v][0] = ybar;
	plist[k++] = v++;
    }

    if (!arma_has_const(ainfo)) {
	(*pZ)[v][0] = 1.0; /* FIXME arbitrary */
    }    

    for (i=1; i<=ainfo->p; i++) {
	if (i > 1) {
	    strcat(fnstr, "+");
	}
	sprintf(term, "phi_%d*y_%d", i, i);
	strcat(fnstr, term);
	sprintf(pdinfo->varname[v], "phi_%d", i);
	plist[k++] = v++;
    }

    if (ainfo->p > 0) {
	strcat(fnstr, "+");
    } else if (!arma_has_const(ainfo)) {
	(*pZ)[v][0] = 1.0; /* FIXME arbitrary */
    }      

    for (i=1; i<=ainfo->P; i++) {
	if (i > 1) {
	    strcat(fnstr, "+");
	}
	sprintf(term, "Phi_%d*y_%d", i, i * ainfo->pd);
	strcat(fnstr, term);
	sprintf(pdinfo->varname[v], "Phi_%d", i);
	plist[k++] = v++;
    }

    for (i=1; i<=ainfo->P; i++) {
	for (j=1; j<=ainfo->p; j++) {
	    sprintf(term, "-phi_%d*Phi_%d*y_%d", j, i, i * ainfo->pd + j);
	    strcat(fnstr, term);
	}
    }  

    for (i=0; i<ainfo->r; i++) {
	j = alist[axstart + i];
	sprintf(term, "+b_%d*%s", i + 1, pdinfo->varname[j]);
	strcat(fnstr, term);
	sprintf(pdinfo->varname[v], "b%d", i + 1);
	plist[k++] = v++;
    }	

    err = nls_spec_set_regression_function(spec, fnstr, pdinfo);

    if (!err) {
	err = nls_spec_add_param_list(spec, plist, (const double **) *pZ,
				      pdinfo);
    }

    if (!err) {
	*amod = model_from_nls_spec(spec, pZ, pdinfo, opt, prn);
	err = amod->errcode;
#if ARMA_DEBUG
	if (!err) {
	    printmodel(amod, pdinfo, OPT_NONE, prn);
	}
	gretl_print_destroy(prn);
#endif
    }

    bailout:

    nls_spec_destroy(spec);
    free(plist);

    return err;
}

/* Run a least squares model to get initial values for the AR
   coefficients */

static int ar_init_by_ls (const int *list, double *coeff,
			  const double **Z, const DATAINFO *pdinfo,
			  struct arma_info *ainfo)
{
    int an = pdinfo->t2 - ainfo->t1 + 1;
    int np = ainfo->p, nq = ainfo->q;
    int nmixed = ainfo->p * ainfo->P;
    int ptotal = ainfo->p + ainfo->P + nmixed;
    int av = ptotal + ainfo->r + 2;
    int ifc = arma_has_const(ainfo);
    const double *y;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *alist = NULL;
    MODEL armod;
    int xstart, axstart;
    int lag, offset;
    int i, j, t, err = 0;

    gretl_model_init(&armod); 

    alist = gretl_list_new(av);
    if (alist == NULL) {
	return 1;
    }

    alist[1] = 1;

    if (ifc) {
	alist[2] = 0;
	offset = 3;
    } else {
	alist[0] -= 1;
	offset = 2;
    }

    for (i=0; i<ainfo->p; i++) {
	alist[i + offset] = 2 + i;
    }

    for (i=0; i<ainfo->P; i++) {
	alist[i + offset + ainfo->p] = ainfo->p + 2 + i;
    }

    for (i=0; i<ainfo->P; i++) {
	for (j=0; j<ainfo->p; j++) {
	    alist[i + offset + ainfo->p + ainfo->P] = 
		ainfo->p + ainfo->P + 2 + i;
	}
    }

    axstart = offset + ptotal;
    for (i=0; i<ainfo->r; i++) {
	alist[axstart + i] = ptotal + 2 + i;
    }

    adinfo = create_new_dataset(&aZ, av, an, 0);
    if (adinfo == NULL) {
	free(alist);
	return 1;
    }

    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    if (ainfo->dy != NULL) {
	y = ainfo->dy;
    } else {
	y = Z[ainfo->yno];
    }

    /* build temporary dataset including lagged vars */
    for (t=0; t<an; t++) {
	int k, s;

	s = t + ainfo->t1;
	aZ[1][t] = y[s];
	strcpy(adinfo->varname[1], "y");

	for (i=0; i<ainfo->p; i++) {
	    lag = i + 1;
	    s = t + ainfo->t1 - lag;
	    k = 2 + i;
	    aZ[k][t] = y[s];
	    sprintf(adinfo->varname[k], "y_%d", lag);
	}

	for (i=0; i<ainfo->P; i++) {
	    lag = ainfo->pd * (i + 1);
	    s = t + ainfo->t1 - lag;
	    k = ainfo->p + 2 + i;
	    aZ[k][t] = y[s];
	    sprintf(adinfo->varname[k], "y_%d", lag);
	    for (j=0; j<ainfo->p; j++) {
		lag = ainfo->pd * (i + 1) + (j + 1);
		s = t + ainfo->t1 - lag;
		k = ainfo->p + ainfo->P + 2 + j;
		aZ[k][t] = y[s];
		sprintf(adinfo->varname[k], "y_%d", lag);
	    }
	}

	s = t + ainfo->t1;

	for (i=0; i<ainfo->r; i++) {
	    k = ptotal + 2 + i;
	    j = list[xstart + i];
	    aZ[k][t] = Z[j][s];
	    strcpy(adinfo->varname[k], pdinfo->varname[j]);
	}
    }

    if (arma_has_seasonal(ainfo)) {
	np += ainfo->P;
	nq += ainfo->Q;
    }

#if ARMA_DEBUG
    printlist(alist, "'alist' in ar_init_by_ls");
#endif

    if (alist[0] == 1) {
	/* arma 0, q model */
	for (i=0; i<nq; i++) {
	    /* insert zeros for MA coeffs */
	    coeff[i + np + ifc] = 0.0;
	} 
	goto exit_init;
    }

    if (nmixed > 0) {
	/* mixed: need to use nonlinear least squares */
	err = arma_get_nls_model(&armod, ainfo, alist, axstart, &aZ, adinfo);
    } else {
	/* just use OLS */
	armod = lsq(alist, &aZ, adinfo, OLS, OPT_A | OPT_Z, 0.0);
	err = armod.errcode;
    }

    if (!err) {
	j = 0;
	for (i=0; i<armod.ncoeff; i++) {
	    if (i == np + ifc) {
		j += nq; /* reserve space for MA coeffs */
	    }
	    coeff[j++] = armod.coeff[i];
	}
	for (i=0; i<nq; i++) {
	    /* insert zeros for MA coeffs */
	    coeff[i + np + ifc] = 0.0;
	} 
    }

#if ARMA_DEBUG
    if (!err) {
	fprintf(stderr, "LS init: ncoeff = %d, nobs = %d\n", 
		armod.ncoeff, armod.nobs);
	for (i=0; i<armod.ncoeff; i++) {
	    fprintf(stderr, " coeff[%d] = %g\n", i, armod.coeff[i]);
	}
    } else {
	fprintf(stderr, "LS init: armod.errcode = %d\n", err);
    }
#endif

 exit_init:

    /* clear everything up */
    free(alist);
    destroy_dataset(aZ, adinfo);
    clear_model(&armod);

    return err;
}

/* set up a model_info struct for passing to bhhh_max */

static model_info *
set_up_arma_model_info (struct arma_info *ainfo)
{
    double tol = get_bhhh_toler();
    model_info *arma;

    if (na(tol)) {
	tol = 1.0e-6;
    }

    arma = model_info_new(ainfo->nc, ainfo->t1, ainfo->t2, ainfo->T, tol);

    if (arma == NULL) return NULL;

    model_info_set_opts(arma, PRESERVE_OPG_MODEL);
    model_info_set_n_series(arma, ainfo->nc + 1);

    /* add pointer to ARMA-specific details */
    model_info_set_extra_info(arma, ainfo);

    return arma;
}

MODEL arma_model (const int *list, const double **Z, const DATAINFO *pdinfo, 
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
	errprn = prn;
    } else {
	errprn = NULL;
    }

    arma_info_init(&ainfo, 0, pdinfo);
    gretl_model_init(&armod); 
    gretl_model_smpl_init(&armod, pdinfo);

    alist = gretl_list_copy(list);
    if (alist == NULL) {
	armod.errcode = E_ALLOC;
	goto bailout;
    }

    err = arma_check_list(alist, opt, Z, pdinfo, &ainfo);
    if (err) {
	armod.errcode = err;
	goto bailout;
    } 

    /* calculate maximum lag */
    calc_max_lag(&ainfo);

    /* adjust sample range if need be */
    if (arma_adjust_sample(pdinfo, Z, alist, &ainfo)) {
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

    /* create differenced series if needed */
    if (ainfo.d > 0 || ainfo.D > 0) {
	err = arima_difference(Z[ainfo.yno], &ainfo);
    }

    /* initialize the coefficients: AR and regression part by least
       squares, MA at 0 */
    err = ar_init_by_ls(alist, coeff, Z, pdinfo, &ainfo);
    if (err) {
	armod.errcode = err;
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

	write_arma_model_stats(pmod, arma, alist, Z, theta, &ainfo);

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

    errprn = NULL;

    return armod;
}
