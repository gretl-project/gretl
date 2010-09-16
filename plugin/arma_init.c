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
#include "gretl_scalar.h"
#include "arma_priv.h"

#define AINIT_DEBUG 0

/* Given an estimate of the ARMA constant via OLS, convert to the form
   wanted for initializing the Kalman filter.  Note: the 'b' array
   goes: const, phi, Phi, theta, Theta, beta.
*/

static void transform_arma_const (double *b, arma_info *ainfo)
{
    const double *phi = b + 1;
    const double *Phi = phi + ainfo->np;
    double narfac = 1.0;
    double sarfac = 1.0;
    int i, k = 0;

#if AINIT_DEBUG
    fprintf(stderr, "transform_arma_const: initially = %g\n", b[0]);
#endif

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    narfac -= phi[k++];
	}
    }

    for (i=0; i<ainfo->P; i++) {
	sarfac -= Phi[i];
    }

    b[0] /= (narfac * sarfac);
}

static int use_preprocessed_y (arma_info *ainfo)
{
    if (arma_xdiff(ainfo)) {
	/* for initialization, use the level of y */
	return 0;
    } else {
	/* use preprocessed y if available */
	return (ainfo->y != NULL);
    }
}

#define HR_MINLAGS 16

static int hr_transcribe_coeffs (arma_info *ainfo,
				 MODEL *pmod, double *b)
{
    const double *theta = NULL;
    const double *Theta = NULL;
    int j = ainfo->nexo + ainfo->ifc;
    int i, k = 0;
    int err = 0;

    if (ainfo->ifc) {
	b[0] = pmod->coeff[0];
	if (arma_xdiff(ainfo)) {
	    b[0] /= ainfo->T;
	}
	k = 1;
    } 

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    b[k++] = pmod->coeff[j++];
	}
    }

    for (i=0; i<ainfo->P; i++) { 
	b[k++] = pmod->coeff[j];
	j += ainfo->np + 1; /* assumes ainfo->p < pd */
    }

    theta = pmod->coeff + j;

    for (i=0; i<ainfo->q; i++) {
	if (MA_included(ainfo, i)) {
	    b[k++] = pmod->coeff[j++];
	}
    }

    Theta = pmod->coeff + j;

    for (i=0; i<ainfo->Q; i++) {
	b[k++] = pmod->coeff[j];
	j += ainfo->nq + 1; /* assumes ainfo->q < pd */
    }

    j = ainfo->ifc;

    for (i=0; i<ainfo->nexo; i++) {
	b[k++] = pmod->coeff[j++];
    }

    /* check MA values? */
    if (ainfo->q > 0 || ainfo->Q > 0) {
	err = ma_out_of_bounds(ainfo, theta, Theta);
	bounds_checker_cleanup();
    }

    return err;
}

/* Hannan-Rissanen ARMA initialization via two OLS passes. In the
   first pass we run an OLS regression of y on the exogenous vars plus
   a certain (biggish) number of lags. In the second we estimate the
   ARMA model by OLS, substituting innovations and corresponding lags
   with the first-pass residuals.
*/

static int real_hr_arma_init (double *coeff, const double **Z, 
			      const DATAINFO *pdinfo,
			      arma_info *ainfo, PRN *prn)
{
    const int *list = ainfo->alist;
    int np = ainfo->p, nq = ainfo->q;
    int nP = ainfo->P, nQ = ainfo->Q;
    int ptotal = np + nP + np * nP;
    int qtotal = nq + nQ + nq * nQ;
    int nexo = ainfo->nexo;
    int pass1lags, pass1v;
    const double *y;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *pass1list = NULL;
    int *pass2list = NULL;
    int *arlags = NULL;
    int *malags = NULL;
    MODEL armod;
    int xstart;
    int m, pos, s;
    int i, j, t;
    int err = 0;

    pass1lags = (ainfo->Q + ainfo->P) * pdinfo->pd;
    if (pass1lags < HR_MINLAGS) {
	pass1lags = HR_MINLAGS;
    }
    pass1v = pass1lags + nexo + 2;

    /* dependent variable */
    if (use_preprocessed_y(ainfo)) {
	y = ainfo->y;
    } else {
	y = Z[ainfo->yno];
    }

    adinfo = create_auxiliary_dataset(&aZ, pass1v + qtotal, ainfo->T);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

#if AINIT_DEBUG
    fprintf(stderr, "hr_arma_init: dataset allocated: %d vars, %d obs\n", 
	    pass1v + qtotal, ainfo->T);
#endif

    /* in case we bomb before estimating a model */
    gretl_model_init(&armod);

    /* Start building stuff for pass 1 */

    pass1list = gretl_list_new(pass1v);
    if (pass1list == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
	
    pass1list[1] = 1;
    pass1list[2] = 0;
    for (i=2; i<pass1v; i++) {
	pass1list[i+1] = i;
    }

    /* variable names */

    strcpy(adinfo->varname[1], "y");
    for (i=0; i<nexo; i++) { 
	/* exogenous vars */
	sprintf(adinfo->varname[i+1], "x%d", i);
    }
    for (i=1; i<=pass1lags; i++) { 
	/* lags */
	sprintf(adinfo->varname[i+1+nexo], "y_%d", i);
    }

     /* Fill the dataset with the data for pass 1 */

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    for (t=0; t<ainfo->T; t++) {
	s = t + ainfo->t1;
	aZ[1][t] = y[s];
	for (i=0, pos=2; i<nexo; i++) {
	    m = list[xstart + i];
	    aZ[pos++][t] = Z[m][s];
	}
	for (i=1; i<=pass1lags; i++) {
	    s = t + ainfo->t1 - i;
	    aZ[pos++][t] = (s >= 0)? y[s] : NADBL;
	}
    }

    /* pass 1 proper */

    armod = lsq(pass1list, aZ, adinfo, OLS, OPT_A);
    if (armod.errcode) {
	err = armod.errcode;
	goto bailout;
    } 

#if AINIT_DEBUG
    fprintf(stderr, "pass1 model: t1=%d, t2=%d, nobs=%d, ncoeff=%d, dfd = %d\n", 
	    armod.t1, armod.t2, armod.nobs, armod.ncoeff, armod.dfd);
#endif

    /* allocations for pass 2 */

    if (qtotal > 0) {
	malags = malloc(qtotal * sizeof *malags);
	if (malags == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0, pos=0; i<nq; i++) {
		malags[pos++] = i+1;
	    }
	    for (i=0; i<ainfo->Q; i++) {
		for (j=0; j<=nq; j++) {
		    malags[pos++] = (i+1) * pdinfo->pd + j;
		}
	    }
	}
    }

    if (ptotal > 0 && !err) {
	arlags = malloc(ptotal * sizeof *arlags);
	if (arlags == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0, pos=0; i<np; i++) {
		arlags[pos++] = i+1;
	    }
	    for (i=0; i<ainfo->P; i++) {
		for (j=0; j<=np; j++) {
		    arlags[pos++] = (i+1) * pdinfo->pd + j;
		}
	    }
	}
    }

    if (!err) {
	pass2list = gretl_list_new(2 + nexo + ptotal + qtotal);
	if (pass2list == NULL) {
	    err = E_ALLOC;
	}
    }

    /* handle error in pass2 allocations */
    if (err) {
	goto bailout;
    }

    /* stick lagged residuals into temp dataset */
    pos = pass1v;
    for (i=0; i<qtotal; i++) {
	sprintf(adinfo->varname[pos], "e_%d", malags[i]);
	for (t=0; t<ainfo->T; t++) {
	    s = t - malags[i];
	    aZ[pos][t] = (s >= 0)? armod.uhat[s] : NADBL;
	}
	pos++;
    }

    /* compose pass 2 regression list */
    for (i=1, pos=1; i<=nexo+2; i++) {
	pass2list[pos++] = pass1list[i];
    }
    for (i=0; i<ptotal; i++) {
	/* FIXME? */
	if (AR_included(ainfo,i)) {
	    pass2list[pos++] = arlags[i] + nexo + 1;
	}
    }
    for (i=0; i<qtotal; i++) {
	/* FIXME? */
	if (MA_included(ainfo,i)) {
	    pass2list[pos++] = pass1v + i;
	}
    }
    
    /* now do pass2 */
    clear_model(&armod);
    armod = lsq(pass2list, aZ, adinfo, OLS, OPT_A);

    if (armod.errcode) {
	err = armod.errcode;
    } else {
#if AINIT_DEBUG
	PRN *modprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	printmodel(&armod, adinfo, OPT_S, modprn);
	gretl_print_destroy(modprn);
#endif
	err = hr_transcribe_coeffs(ainfo, &armod, coeff);

	if (!err && arma_exact_ml(ainfo) && 
	    ainfo->ifc && ainfo->nexo == 0) {
	    transform_arma_const(coeff, ainfo);
	}
    }

#if AINIT_DEBUG
    if (!err) {
	fprintf(stderr, "HR init:\n");
	for (i=0; i<ainfo->nc; i++) {
	    fprintf(stderr, "coeff[%d] = %g\n", i, coeff[i]);
	}
    }
#endif

 bailout:

    free(pass1list);
    free(pass2list);
    free(arlags);
    free(malags);
    destroy_dataset(aZ, adinfo);
    clear_model(&armod);

    if (!err && prn != NULL) {
	pputs(prn, "\narma initialization: using Hannan-Rissanen method\n\n");
    }

    return err;
}

/* Do we have enough observations to do Hannan-Rissanen? */

static int hr_df_check (arma_info *ainfo, const DATAINFO *pdinfo)
{
    int nobs = ainfo->T;
    int nlags = (ainfo->P + ainfo->Q) * pdinfo->pd;
    int ncoeff, df;
    int ok = 1;

    if (nlags < HR_MINLAGS) {
	nlags = HR_MINLAGS;
    }

    ncoeff = nlags + ainfo->nexo + ainfo->ifc;
    nobs -= nlags;
    df = nobs - ncoeff;

    if (df < 1) {
	ok = 0;
    }

#if AINIT_DEBUG
    fprintf(stderr, "hr_init_check: ncoeff=%d, nobs=%d, 'df'=%d\n", 
	    ncoeff, nobs, df);
#endif

    return ok;
}

int hr_arma_init (double *coeff, const double **Z, 
		  const DATAINFO *pdinfo,
		  arma_info *ainfo, int *done)
{
    int ok = hr_df_check(ainfo, pdinfo);
    int err = 0;

    if (ok) {
	err = real_hr_arma_init(coeff, Z, pdinfo, ainfo, ainfo->prn);
	if (!err) {
	    *done = 1;
	}
    }

#if AINIT_DEBUG
    if (*done) {
	fputs("*** hr_arma_init OK\n", stderr);
    } else {
	fputs("*** hr_arma_init failed, will try ar_arma_init\n", stderr);
    } 
#endif

    return err;
}

/* try to avoid numerical problems when doing exact ML: 
   scale the dependent variable if it's "too big" 
*/

static void maybe_rescale_y (arma_info *ainfo, const double **Z,
			     const DATAINFO *pdinfo)
{
    double ybar;
    int t, doit = 0;

    if (ainfo->y != NULL) {
	ybar = gretl_mean(ainfo->t1, ainfo->t2, ainfo->y);
	doit = (fabs(ybar) > 250);
    } else {
	const double *y = Z[ainfo->yno];

	ybar = gretl_mean(ainfo->t1, ainfo->t2, y);
	if (fabs(ybar) > 250) {
	    ainfo->y = malloc(pdinfo->n * sizeof *ainfo->y);
	    if (ainfo->y != NULL) {
		for (t=0; t<pdinfo->n; t++) {
		    ainfo->y[t] = y[t];
		}
		doit = 1;
	    }
	}
    }

    if (doit) {
	fprintf(stderr, "arma: ybar = %g, rescaling y\n", ybar);
	for (t=0; t<=ainfo->t2; t++) {
	    if (!na(ainfo->y[t])) {
		ainfo->y[t] /= ybar;
	    }
	}
	ainfo->yscale = ybar;
    }
}

/* transcribe coeffs from the OLS or NLS model used for initializing,
   into the array @b that will be passed to the maximizer.
*/

static void arma_init_transcribe_coeffs (arma_info *ainfo,
					 MODEL *pmod, double *b)
{
    int q0 = ainfo->ifc + ainfo->np + ainfo->P;
    int Q0 = q0 + ainfo->nq;
    int i, j = 0;

    for (i=0; i<pmod->ncoeff; i++) {
	if (i == q0) {
	    /* reserve space for nonseasonal MA */
	    j += ainfo->nq;
	} 
	if (i == Q0) {
	    /* and for seasonal MA */
	    j += ainfo->Q;
	}
	b[j++] = pmod->coeff[i];
    }

    if (arma_xdiff(ainfo) && ainfo->ifc) {
	/* is this a good idea? */
	b[0] /= ainfo->T;
    }

    /* insert near-zeros for nonseasonal MA */
    for (i=0; i<ainfo->nq; i++) {
	b[q0 + i] = 0.0001;
    } 

    /* and also seasonal MA */
    for (i=0; i<ainfo->Q; i++) {
	b[Q0 + i] = 0.0001;
    }	
}

/* compose variable names for temporary dataset */

static void arma_init_add_varnames (arma_info *ainfo, 
				    int ptotal, int narmax, 
				    DATAINFO *adinfo)
{
    int i, j, k, kx, ky;
    int lag, k0 = 2;

    strcpy(adinfo->varname[1], "y");

    k = k0;
    kx = ptotal + ainfo->nexo + k0;

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(adinfo->varname[k++], "y_%d", lag);
	    for (j=0; j<narmax; j++) {
		sprintf(adinfo->varname[kx++], "x%d_%d", j+1, lag);
	    }
	}
    }

    ky = ainfo->np + ainfo->P + k0;

    for (j=0; j<ainfo->P; j++) {
	lag = (j + 1) * ainfo->pd;
	k = k0 + ainfo->np + j;
	sprintf(adinfo->varname[k], "y_%d", lag);
	for (i=0; i<narmax; i++) {
	    sprintf(adinfo->varname[kx++], "x%d_%d", i+1, lag);
	}
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		lag = (j + 1) * ainfo->pd + (i + 1);
		sprintf(adinfo->varname[ky++], "y_%d", lag);
		for (k=0; k<narmax; k++) {
		    sprintf(adinfo->varname[kx++], "x%d_%d", k+1, lag);
		}
	    }
	}
    }

    kx = ptotal + k0;

    for (i=0; i<ainfo->nexo; i++) {
	sprintf(adinfo->varname[kx++], "x%d", i+1);
    }
}

/* Build temporary dataset including lagged vars: if we're doing exact
   ML on an ARMAX model we need lags of the exogenous variables as
   well as lags of y_t.  Note that the auxiliary dataset has "t = 0"
   at an offset of ainfo->t1 into the "real", external dataset.
*/

static void arma_init_build_dataset (arma_info *ainfo, 
				     int ptotal, int narmax, 
				     const int *list,
				     const double **Z,
				     double **aZ, 
				     DATAINFO *adinfo)
{
    const double *y;
    int i, j, k, kx, ky;
    int t, s, m, k0 = 2;
    int lag, xstart;

    /* add variable names to auxiliary dataset */
    arma_init_add_varnames(ainfo, ptotal, narmax, adinfo);

    /* dependent variable */
    if (use_preprocessed_y(ainfo)) {
	y = ainfo->y;
    } else {
	y = Z[ainfo->yno];
    }

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    for (t=0; t<adinfo->n; t++) {
	int realt = t + ainfo->t1;
	int miss = 0;

	aZ[1][t] = y[realt];

	k = k0;
	kx = ptotal + ainfo->nexo + k0;

	for (i=0; i<ainfo->p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    lag = i + 1;
	    s = realt - lag;
	    if (s < 0) {
		miss = 1;
		aZ[k++][t] = NADBL;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k++][t] = y[s];
		for (j=0; j<narmax; j++) {
		    m = list[xstart + j];
		    aZ[kx++][t] = Z[m][s];
		}
	    }
	}

	ky = ainfo->np + ainfo->P + k0;

	for (j=0; j<ainfo->P; j++) {
	    lag = (j + 1) * ainfo->pd;
	    s = realt - lag;
	    k = ainfo->np + k0 + j;
	    if (s < 0) {
		miss = 1;
		aZ[k][t] = NADBL;
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k][t] = y[s];
		for (k=0; k<narmax; k++) {
		    m = list[xstart + k];
		    aZ[kx++][t] = Z[m][s];
		}
	    }
	    for (i=0; i<ainfo->p; i++) {
		if (!AR_included(ainfo, i)) {
		    continue;
		}
		lag = (j + 1) * ainfo->pd + (i + 1);
		s = realt - lag;
		if (s < 0) {
		    miss = 1;
		    aZ[ky++][t] = NADBL;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = NADBL;
		    }
		} else {
		    aZ[ky++][t] = y[s];
		    for (k=0; k<narmax; k++) {
			m = list[xstart + k];
			aZ[kx++][t] = Z[m][s];
		    }
		}
	    }
	}

	kx = ptotal + k0;

	for (i=0; i<ainfo->nexo; i++) {
	    m = list[xstart + i];
	    aZ[kx++][t] = Z[m][realt];
	}

	if (miss) {
	    adinfo->t1 = t + 1;
	}	
    }

#if AINIT_DEBUG
    fprintf(stderr, "arma init dataset:\n");
    for (i=0; i<adinfo->v; i++) {
	fprintf(stderr, "var %d '%s', obs[0] = %g\n", i, adinfo->varname[i], 
		aZ[i][0]);
    }
#endif
}

static void nls_kickstart (MODEL *pmod, double **Z, 
			   DATAINFO *pdinfo,
			   double *b0, double *by1)
{
    int list[4];

    if (b0 != 0) {
	list[0] = 3;
	list[1] = 1;
	list[2] = 0;
	list[3] = 2;
    } else {
	list[0] = 2;
	list[1] = 1;
	list[2] = 2;
    }

    *pmod = lsq(list, Z, pdinfo, OLS, OPT_A | OPT_Z);

    if (!pmod->errcode) {
	if (b0 != 0) {
	    *b0 = pmod->coeff[0];
	    *by1 = pmod->coeff[1];
	} else {
	    *by1 = pmod->coeff[0];
	}
    }

    clear_model(pmod);
}

static int add_to_spec (char *targ, const char *src)
{
    if (strlen(src) + strlen(targ) > MAXLINE - 1) {
	return 1;
    } else {
	strcat(targ, src);
	return 0;
    }
}

/* for ARMAX: write the component of the NLS specification
   that takes the form (y_{t-i} - X_{t-i} \beta)
*/

static int y_Xb_at_lag (char *spec, arma_info *ainfo, 
			int narmax, int lag)
{
    char chunk[32];
    int i, nt;
    int err = 0;

    if (narmax == 0) {
	sprintf(chunk, "y_%d", lag);
	return add_to_spec(spec, chunk);
    }

    nt = ainfo->ifc + narmax;

    sprintf(chunk, "(y_%d-", lag);

    if (nt > 1) {
	strcat(chunk, "(");
    }

    if (ainfo->ifc) {
	strcat(chunk, "b0");
    }

    err = add_to_spec(spec, chunk);

    for (i=0; i<narmax && !err; i++) {
	if (ainfo->ifc || i > 0) {
	    err += add_to_spec(spec, "+");
	} 
	sprintf(chunk, "b%d*x%d_%d", i+1, i+1, lag);
	err += add_to_spec(spec, chunk); 
    }

    if (nt > 1) {
	err += add_to_spec(spec, "))");
    } else {
	err += add_to_spec(spec, ")");
    }

    return err;
}

static int arma_get_nls_model (MODEL *amod, arma_info *ainfo,
			       int narmax, const double *coeff,
			       double ***pZ, DATAINFO *pdinfo,
			       PRN *prn) 
{
    gretlopt nlsopt = OPT_A;
    char fnstr[MAXLINE];
    char term[32];
    nlspec *spec;
    double *parms = NULL;
    char **pnames = NULL;
    double *b0 = NULL, *by1 = NULL;
    int nparam, lag;
    int i, j, k, err = 0;

    spec = nlspec_new(NLS, pdinfo);
    if (spec == NULL) {
	return E_ALLOC;
    }

    if (arma_least_squares(ainfo)) {
	/* respect verbose option */
	if (prn != NULL) {
	    nlsopt |= OPT_V;
	}
    } else {
#if AINIT_DEBUG
	nlsopt |= OPT_V;
#else
	/* don't bother with standard errors */
	nlsopt |= OPT_C;
#endif
    }

    nlspec_set_t1_t2(spec, 0, ainfo->t2 - ainfo->t1); /* ?? */

    nparam = ainfo->ifc + ainfo->np + ainfo->P + ainfo->nexo;

    parms = malloc(nparam * sizeof *parms);
    if (parms == NULL) {
	err = E_ALLOC;
	goto bailout;
    }	

    pnames = strings_array_new_with_length(nparam, VNAMELEN);
    if (pnames == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* make names for the parameters; construct the param list;
       and do some rudimentary fall-back initialization */

    for (i=0; i<nparam; i++) {
	parms[i] = 0.0;
    }

    k = 0;

    if (ainfo->ifc) {
	if (coeff != NULL) {
	    parms[k] = coeff[k];
	} else {
	    parms[k] = gretl_mean(0, pdinfo->n - 1, (*pZ)[1]);
	}
	b0 = &parms[k];
	strcpy(pnames[k++], "b0");
    }

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    if (by1 == NULL) {
		by1 = &parms[k];
		if (coeff == NULL) {
		    parms[k] = 0.1;
		}
	    }
	    if (coeff != NULL) {
		parms[k] = coeff[k];
	    }
	    sprintf(pnames[k++], "phi%d", i+1);
	}
    }

    for (i=0; i<ainfo->P; i++) {
	if (by1 == NULL) {
	    by1 = &parms[k];
	    if (coeff == NULL) {
		parms[k] = 0.1;
	    }
	}
	if (coeff != NULL) {
	    parms[k] = coeff[k];
	}
	sprintf(pnames[k++], "Phi%d", i+1);
    }

    for (i=0; i<ainfo->nexo; i++) {
	if (coeff != NULL) {
	    parms[k] = coeff[k];
	}
	sprintf(pnames[k++], "b%d", i+1);
    }

    /* construct NLS specification */

    strcpy(fnstr, "y=");

    if (ainfo->ifc) {
	strcat(fnstr, "b0");
    } else {
	strcat(fnstr, "0");
    } 

    for (i=0; i<ainfo->p && !err; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(term, "+phi%d*", lag);
	    err = add_to_spec(fnstr, term);
	    if (!err) {
		err = y_Xb_at_lag(fnstr, ainfo, narmax, lag);
	    }
	}
    }

    for (j=0; j<ainfo->P && !err; j++) {
	sprintf(term, "+Phi%d*", j+1);
	strcat(fnstr, term);
	lag = (j + 1) * ainfo->pd;
	y_Xb_at_lag(fnstr, ainfo, narmax, lag);
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		sprintf(term, "-phi%d*Phi%d*", i+1, j+1);
		err = add_to_spec(fnstr, term);
		if (!err) {
		    lag = (j+1) * ainfo->pd + (i+1);
		    y_Xb_at_lag(fnstr, ainfo, narmax, lag);
		}
	    }
	}
    }

    for (i=0; i<ainfo->nexo && !err; i++) {
	sprintf(term, "+b%d*x%d", i+1, i+1);
	err = add_to_spec(fnstr, term);
    }

    if (!err) {
	if (coeff == NULL) {
	    nls_kickstart(amod, *pZ, pdinfo, b0, by1);
	}

#if AINIT_DEBUG
	fprintf(stderr, "initting using NLS spec:\n %s\n", fnstr);
	for (i=0; i<nparam; i++) {
	    fprintf(stderr, "initial NLS b[%d] = %g (%s)\n",
		    i, parms[i], pnames[i]);
	}
#endif

	err = nlspec_set_regression_function(spec, fnstr, pdinfo);
    }

    if (!err) {
	set_auxiliary_scalars();
	err = nlspec_add_param_list(spec, nparam, parms, pnames,
				    pZ, pdinfo);

	if (!err) {
	    *amod = model_from_nlspec(spec, pZ, pdinfo, nlsopt, prn);
	    err = amod->errcode;
#if AINIT_DEBUG
	    if (!err) {
		printmodel(amod, pdinfo, OPT_NONE, prn);
	    }
#endif
	}
	unset_auxiliary_scalars();
    }

 bailout:

    nlspec_destroy(spec);
    free(parms);
    free_strings_array(pnames, nparam);

    return err;
}

/* compose the regression list for the case where we're initializing
   ARMA via plain OLS (not NLS)
*/

static int *make_ar_ols_list (arma_info *ainfo, int av)
{
    int *list = gretl_list_new(av);
    int i, k, vi;

    if (list == NULL) {
	return NULL;
    }

    list[1] = 1;

    if (ainfo->ifc) {
	list[2] = 0;
	k = 3;
    } else {
	list[0] -= 1;
	k = 2;
    }

    /* allow for const and y */
    vi = 2;

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    list[k++] = vi++;
	}
    }

    for (i=0; i<ainfo->P; i++) {
	list[k++] = vi++;
    }

    for (i=0; i<ainfo->nexo; i++) {
	list[k++] = vi++;
    }

    return list;
}

/* Run a least squares model to get initial values for the AR
   coefficients, either OLS or NLS.  We use NLS if there is
   nonlinearity due to either (a) the presence of both a seasonal and
   a non-seasonal AR component or (b) the presence of exogenous
   variables in the context of a non-zero AR order, where estimation
   will be via exact ML.  

   In this initialization any MA coefficients are simply set to
   near-zero.
*/

int ar_arma_init (double *coeff, const double **Z, 
		  const DATAINFO *pdinfo,
		  arma_info *ainfo, MODEL *pmod)
{
    PRN *prn = ainfo->prn;
    int *list = ainfo->alist;
    int nmixed = ainfo->np * ainfo->P;
    int ptotal = ainfo->np + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *arlist = NULL;
    MODEL armod;
    int narmax, nonlin = 0;
    int i, err = 0;

#if AINIT_DEBUG
    fprintf(stderr, "ar_arma_init: pdinfo->t1=%d, pdinfo->t2=%d (pdinfo->n=%d);\n"
	    " ainfo->t1=%d, ainfo->t2=%d, ",
	    pdinfo->t1, pdinfo->t2, pdinfo->n, ainfo->t1, ainfo->t2);
    fprintf(stderr, "nmixed = %d, ptotal = %d\n", nmixed, ptotal);
#endif

    if (ptotal == 0 && ainfo->nexo == 0 && !ainfo->ifc) {
	/* special case of pure MA model */
	for (i=0; i<ainfo->nq + ainfo->Q; i++) {
	    coeff[i] = 0.0001; 
	} 
#if AINIT_DEBUG
	fprintf(stderr, " pure MA: just setting small coeff value(s)\n");
#endif
	return 0;
    }

    gretl_model_init(&armod); 

    narmax = (arma_exact_ml(ainfo))? ainfo->nexo : 0;
    if (narmax > 0) {
	/* ARMAX-induced lags of exog vars */
	av += ainfo->nexo * ptotal;
    } 

    if (arma_exact_ml(ainfo) && ainfo->ifc) {
	maybe_rescale_y(ainfo, Z, pdinfo);
    }

    adinfo = create_auxiliary_dataset(&aZ, av, ainfo->T);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

    if (ptotal > 0 && (narmax > 0 || nmixed > 0)) {
	/* we'll have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	arlist = make_ar_ols_list(ainfo, av);
    }

    /* build temporary dataset */
    arma_init_build_dataset(ainfo, ptotal, narmax, list,
			    Z, aZ, adinfo);

    if (nonlin) {
	PRN *dprn = NULL;

#if AINIT_DEBUG
	fprintf(stderr, "arma:_init_by_ls: doing NLS\n");
	dprn = prn;
#endif
	err = arma_get_nls_model(&armod, ainfo, narmax, NULL, &aZ, adinfo,
				 dprn);
    } else {
#if AINIT_DEBUG
	printlist(arlist, "'arlist' in ar_arma_init (OLS)");
#endif
	armod = lsq(arlist, aZ, adinfo, OLS, OPT_A | OPT_Z);
	err = armod.errcode;
    }

#if AINIT_DEBUG
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

    if (!err) {
	arma_init_transcribe_coeffs(ainfo, &armod, coeff);
    }

    /* handle the case where we need to translate from an
       estimate of the regression constant to the
       unconditional mean of y_t
    */
    if (!err && arma_exact_ml(ainfo) && ainfo->ifc && 
	(!nonlin || ainfo->nexo == 0)) {
	transform_arma_const(coeff, ainfo);
    }

    if (!err && prn != NULL) {
	if (nonlin) {
	    pprintf(prn, "\n%s: %s\n\n", _("ARMA initialization"),
		    _("using nonlinear AR model"));
	} else {
	    pprintf(prn, "\n%s: %s\n\n", _("ARMA initialization"),
		    _("using linear AR model"));
	}
    }

    /* clean up */
    clear_model(&armod);
    free(arlist);
    destroy_dataset(aZ, adinfo);

    return err;
}

int arma_by_ls (const double *coeff, 
		const double **Z, const DATAINFO *pdinfo,
		arma_info *ainfo, MODEL *pmod)
{
    PRN *prn = ainfo->prn;
    int *list = ainfo->alist;
    int nmixed = ainfo->np * ainfo->P;
    int ptotal = ainfo->np + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    double **aZ = NULL;
    DATAINFO *adinfo = NULL;
    int *arlist = NULL;
    int nonlin = 0;

    adinfo = create_auxiliary_dataset(&aZ, av, ainfo->T);
    if (adinfo == NULL) {
	return E_ALLOC;
    }

    if (ptotal > 0 && nmixed > 0) {
	/* we'll have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	arlist = make_ar_ols_list(ainfo, av);
    }

    /* build temporary dataset */
    arma_init_build_dataset(ainfo, ptotal, 0, list,
			    Z, aZ, adinfo);

    if (nonlin) {
	pmod->errcode = arma_get_nls_model(pmod, ainfo, 0, coeff, &aZ, adinfo,
					   prn);
    } else {
	*pmod = lsq(arlist, aZ, adinfo, OLS, OPT_A | OPT_Z);
    }

    /* clean up */
    free(arlist);
    destroy_dataset(aZ, adinfo);

    if (!pmod->errcode && pmod->full_n < pdinfo->n) {
	/* the model series are short */
	double *uhat = malloc(pdinfo->n * sizeof *uhat);
	double *yhat = malloc(pdinfo->n * sizeof *yhat);
	int s, t;

	if (uhat == NULL || yhat == NULL) {
	    free(uhat);
	    free(yhat);
	    pmod->errcode = E_ALLOC;
	} else {
	    for (t=0; t<pdinfo->n; t++) {
		uhat[t] = yhat[t] = NADBL;
	    }
	    t = ainfo->t1;
	    for (s=0; s<pmod->full_n; s++, t++) {
		uhat[t] = pmod->uhat[s];
		yhat[t] = pmod->yhat[s];
	    }
	    free(pmod->uhat);
	    pmod->uhat = uhat;
	    free(pmod->yhat);
	    pmod->yhat = yhat;
	}
    }

    return pmod->errcode;
}



