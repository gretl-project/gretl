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
#include "uservar.h"
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

    if (ainfo->np == 0 && ainfo->P == 0) {
	return;
    }

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

static int real_hr_arma_init (double *coeff, const DATASET *dset,
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
    DATASET *aset = NULL;
    int *pass1list = NULL;
    int *pass2list = NULL;
    int *arlags = NULL;
    int *malags = NULL;
    MODEL armod;
    int xstart;
    int m, pos, s;
    int i, j, t;
    int err = 0;

    pass1lags = (ainfo->Q + ainfo->P) * dset->pd;
    if (pass1lags < HR_MINLAGS) {
	pass1lags = HR_MINLAGS;
    }
    pass1v = pass1lags + nexo + 2;

    /* dependent variable */
    if (arma_xdiff(ainfo)) {
	/* for initialization, use the level of y */
	y = dset->Z[ainfo->yno];
    } else { 
	y = ainfo->y;
    } 

    aset = create_auxiliary_dataset(pass1v + qtotal, ainfo->T, 0);
    if (aset == NULL) {
	return E_ALLOC;
    }

#if AINIT_DEBUG
    fprintf(stderr, "hr_arma_init: dataset allocated: %d vars, %d obs\n", 
	    pass1v + qtotal, ainfo->T);
#endif

    /* in case we bomb before estimating a model */
    gretl_model_init(&armod, dset);

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

    strcpy(aset->varname[1], "y");
    for (i=0; i<nexo; i++) { 
	/* exogenous vars */
	sprintf(aset->varname[i+1], "x%d", i);
    }
    for (i=1; i<=pass1lags; i++) { 
	/* lags */
	sprintf(aset->varname[i+1+nexo], "y_%d", i);
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
	aset->Z[1][t] = y[s];
	for (i=0, pos=2; i<nexo; i++) {
	    m = list[xstart + i];
	    aset->Z[pos++][t] = dset->Z[m][s];
	}
	for (i=1; i<=pass1lags; i++) {
	    s = t + ainfo->t1 - i;
	    aset->Z[pos++][t] = (s >= 0)? y[s] : NADBL;
	}
    }

    /* pass 1 proper */

    armod = lsq(pass1list, aset, OLS, OPT_A);
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
		    malags[pos++] = (i+1) * dset->pd + j;
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
		    arlags[pos++] = (i+1) * dset->pd + j;
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
	sprintf(aset->varname[pos], "e_%d", malags[i]);
	for (t=0; t<ainfo->T; t++) {
	    s = t - malags[i];
	    aset->Z[pos][t] = (s >= 0)? armod.uhat[s] : NADBL;
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
    armod = lsq(pass2list, aset, OLS, OPT_A);

    if (armod.errcode) {
	err = armod.errcode;
    } else {
#if AINIT_DEBUG
	PRN *modprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	printmodel(&armod, aset, OPT_S, modprn);
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
    destroy_dataset(aset);
    clear_model(&armod);

    if (!err && prn != NULL) {
	pprintf(prn, "\n%s: %s\n\n", _("ARMA initialization"), 
		_("Hannan-Rissanen method"));
    }

    return err;
}

/* Do we have enough observations to do Hannan-Rissanen? */

static int hr_df_check (arma_info *ainfo, const DATASET *dset)
{
    int nobs = ainfo->T;
    int nlags = (ainfo->P + ainfo->Q) * dset->pd;
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

int hr_arma_init (double *coeff, const DATASET *dset,
		  arma_info *ainfo, int *done)
{
    int ok = hr_df_check(ainfo, dset);
    int err = 0;

    if (ok) {
	err = real_hr_arma_init(coeff, dset, ainfo, ainfo->prn);
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
   arrange for scaling of the dependent variable if it's 
   "too big" 
*/

static void maybe_set_yscale (arma_info *ainfo)
{
    double ybar = gretl_mean(ainfo->t1, ainfo->t2, ainfo->y);

    if (fabs(ybar) > 250) {
	if (arima_levels(ainfo)) {
	    set_arma_avg_ll(ainfo); /* is this a good idea? */
	} else {
	    ainfo->yscale = 10 / ybar;
	}
    }
}

/* transcribe coeffs from the OLS or NLS model used for initializing,
   into the array @b that will be passed to the maximizer.
*/

static void arma_init_transcribe_coeffs (arma_info *ainfo,
					 MODEL *pmod, double *b)
{
    int q0 = ainfo->ifc + ainfo->np + ainfo->P;
    int totq = ainfo->nq + ainfo->Q;
    int i, j = 0;

    for (i=0; i<pmod->ncoeff; i++) {
	if (i == q0 && totq > 0) {
	    /* reserve space for MA terms */
	    j += totq;
	} 
	if (j < ainfo->nc) {
	    b[j++] = pmod->coeff[i];
	}
    }

    if (arma_xdiff(ainfo) && ainfo->ifc) {
	/* is this a good idea? */
	b[0] /= ainfo->T;
    }

    /* insert near-zeros for MA terms */
    for (i=0; i<totq; i++) {
	b[q0 + i] = 0.0001;
    } 
}

/* compose variable names for temporary dataset */

static void arma_init_add_varnames (arma_info *ainfo, 
				    int ptotal, int narmax, 
				    DATASET *aset)
{
    int i, j, k, kx, ky;
    int lag, k0 = 2;

    strcpy(aset->varname[1], "y");

    k = k0;
    kx = ptotal + ainfo->nexo + k0;

    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(aset->varname[k++], "y_%d", lag);
	    for (j=0; j<narmax; j++) {
		sprintf(aset->varname[kx++], "x%d_%d", j+1, lag);
	    }
	}
    }

    ky = ainfo->np + ainfo->P + k0;

    for (j=0; j<ainfo->P; j++) {
	lag = (j + 1) * ainfo->pd;
	k = k0 + ainfo->np + j;
	sprintf(aset->varname[k], "y_%d", lag);
	for (i=0; i<narmax; i++) {
	    sprintf(aset->varname[kx++], "x%d_%d", i+1, lag);
	}
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		lag = (j + 1) * ainfo->pd + (i + 1);
		sprintf(aset->varname[ky++], "y_%d", lag);
		for (k=0; k<narmax; k++) {
		    sprintf(aset->varname[kx++], "x%d_%d", k+1, lag);
		}
	    }
	}
    }

    kx = ptotal + k0;

    for (i=0; i<ainfo->nexo; i++) {
	sprintf(aset->varname[kx++], "x%d", i+1);
    }
}

/* experimental: when initializing an AR(I)MA model via
   NLS, work around interior NAs by adding observation-
   specific dummies to the dataset
*/

static int arma_init_add_dummies (arma_info *ainfo,
				  DATASET *dset)
{
    int *misslist = NULL;
    int t1 = dset->t1;
    int i, t, err = 0;

    /* if we have a block of leading NAs, skip it */

    for (t=t1; t<=dset->t2 && !err; t++) {
	int miss = 0;

	for (i=1; i<dset->v; i++) {
	    if (na(dset->Z[i][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) {
	    t1++;
	} else {
	    break;
	}
    }

    /* form list of observation indices of interior NAs */

    for (t=t1; t<=dset->t2 && !err; t++) {
	for (i=1; i<dset->v; i++) {
	    if (na(dset->Z[i][t])) {
		misslist = gretl_list_append_term(&misslist, t);
		if (misslist == NULL) {
		    err = E_ALLOC;
		}
		break;
	    }
	}
    }

#if AINIT_DEBUG
    printlist(misslist, "arma_init_add_dummies: misslist");
#endif

    if (misslist != NULL) {
	/* For each observation with any missing values, add
	   a specific dummy and zero out the missing data.
	*/
	int origv = dset->v;
	int j, v, nd = misslist[0];

	err = dataset_add_series(dset, nd);
	if (!err) {
	    for (i=1; i<=misslist[0]; i++) {
		v = origv + i - 1;
		t = misslist[i];
		sprintf(dset->varname[v], "d%d", i);
		dset->Z[v][t] = 1.0;
		for (j=1; j<origv; j++) {
		    if (na(dset->Z[j][t])) {
			dset->Z[j][t] = 0.0;
		    }
		}
	    }
	}
    }

    ainfo->misslist = misslist;

    return err;
}

static const int *xlist;

/* X, if non-NULL, holds the differenced regressors */

static double get_xti (const DATASET *dset, int i, int t, 
		       const gretl_matrix *X)
{
    if (X != NULL) {
	return gretl_matrix_get(X, t, i);
    } else {
	return dset->Z[xlist[i]][t];
    }
}

/* Build temporary dataset including lagged vars: if we're doing exact
   ML on an ARMAX model we need lags of the exogenous variables as
   well as lags of y_t.  Note that the auxiliary dataset has "t = 0"
   at an offset of ainfo->t1 into the "real", external dataset.
*/

static int arma_init_build_dataset (arma_info *ainfo, 
				    int ptotal, int narmax, 
				    const int *list,
				    const DATASET *dset,
				    DATASET *aset,
				    int nonlin)
{
    double **aZ = aset->Z;
    const double *y;
    const gretl_matrix *X = NULL;
    int i, j, k, kx, ky;
    int t, s, k0 = 2;
    int undo_diff = 0;
    int xstart;
    int err = 0;

    if (arima_levels(ainfo)) {
	/* we'll need differences for initialization */
	err = arima_difference(ainfo, dset, 1);
	if (err) {
	    return err;
	}
	undo_diff = 1;
	y = ainfo->y;
	X = ainfo->dX;
    } else if (arma_xdiff(ainfo)) {
	/* run init in levels (FIXME?) */
	y = dset->Z[ainfo->yno];
    } else {
	y = ainfo->y;
    }

    /* add variable names to auxiliary dataset */
    arma_init_add_varnames(ainfo, ptotal, narmax, aset);

    /* starting position for reading exogeneous vars */
    if (ainfo->d > 0 || ainfo->D > 0) {
	xstart = (arma_has_seasonal(ainfo))? 10 : 6;
    } else {
	xstart = (arma_has_seasonal(ainfo))? 8 : 5;
    }

    /* set "local" globals */
    xlist = list + xstart;

    for (t=0; t<aset->n; t++) {
	int realt = t + ainfo->t1;
	int miss = 0;

	if (ainfo->yscale != 1.0 && !na(y[realt])) {
	    aZ[1][t] = y[realt] * ainfo->yscale;
	} else {
	    aZ[1][t] = y[realt];
	}

	k = k0;
	kx = ptotal + ainfo->nexo + k0;

	for (i=0; i<ainfo->p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    s = realt - (i + 1);
	    if (s < 0) {
		miss = 1;
		aZ[k++][t] = NADBL;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k][t] = y[s];
		if (ainfo->yscale != 1.0 && !na(y[s])) {
		    aZ[k][t] *= ainfo->yscale;
		}
		k++;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = get_xti(dset, j, s, X);
		}
	    }
	}

	ky = ainfo->np + ainfo->P + k0;

	for (j=0; j<ainfo->P; j++) {
	    s = realt - (j + 1) * ainfo->pd;
	    k = ainfo->np + k0 + j;
	    if (s < 0) {
		miss = 1;
		aZ[k][t] = NADBL;
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k][t] = y[s];
		if (ainfo->yscale != 1.0 && !na(y[s])) {
		    aZ[k][t] *= ainfo->yscale;
		}		
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = get_xti(dset, k, s, X);
		}
	    }
	    for (i=0; i<ainfo->p; i++) {
		if (!AR_included(ainfo, i)) {
		    continue;
		}
		s = realt - ((j + 1) * ainfo->pd + (i + 1));
		if (s < 0) {
		    miss = 1;
		    aZ[ky++][t] = NADBL;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = NADBL;
		    }
		} else {
		    aZ[ky][t] = y[s];
		    if (ainfo->yscale != 1.0 && !na(y[s])) {
			aZ[ky][t] *= ainfo->yscale;
		    }
		    ky++;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = get_xti(dset, k, s, X);
		    }
		}
	    }
	}

	kx = ptotal + k0;

	for (i=0; i<ainfo->nexo; i++) {
	    aZ[kx++][t] = get_xti(dset, i, realt, X);
	}

	if (miss) {
	    aset->t1 = t + 1;
	}	
    }

    if (nonlin && arma_missvals(ainfo)) {
	err = arma_init_add_dummies(ainfo, aset);
    }

    if (undo_diff) {
	arima_difference_undo(ainfo, dset);
    }

#if AINIT_DEBUG
    fprintf(stderr, "arma init dataset:\n");
    for (i=0; i<aset->v; i++) {
	fprintf(stderr, "var %d '%s', obs[0] = %g\n", i, aset->varname[i], 
		aset->Z[i][0]);
    }
#endif

    return err;
}

static void nls_kickstart (MODEL *pmod, DATASET *dset,
			   double *b0, double *by1)
{
    int list[4];

    if (b0 != NULL) {
	list[0] = 3;
	list[1] = 1;
	list[2] = 0;
	list[3] = 2;
    } else {
	list[0] = 2;
	list[1] = 1;
	list[2] = 2;
    }

    *pmod = lsq(list, dset, OLS, OPT_A | OPT_Z);

    if (!pmod->errcode) {
	if (b0 != NULL) {
	    *b0 = pmod->coeff[0];
	    *by1 = pmod->coeff[1];
	} else {
	    *by1 = pmod->coeff[0];
	}
	if (*by1 >= 1.0) {
	   *by1 = 0.95;
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
			       DATASET *dset, PRN *prn) 
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

    spec = nlspec_new(NLS, dset);
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

    nlspec_set_t1_t2(spec, 0, ainfo->T - 1);

    nparam = ainfo->ifc + ainfo->np + ainfo->P + ainfo->nexo;

    if (ainfo->misslist != NULL) {
	nparam += ainfo->misslist[0];
    }

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
	    parms[k] = gretl_mean(0, dset->n - 1, dset->Z[1]);
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

    if (ainfo->misslist != NULL) {
	for (i=1; i<=ainfo->misslist[0]; i++) {
	    j = ainfo->misslist[i];
	    parms[k] = dset->Z[1][j];
	    sprintf(pnames[k++], "c%d", i);
	}
    }

    /* construct NLS specification string */

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

    if (!err && ainfo->misslist != NULL) {
	for (i=1; i<=ainfo->misslist[0]; i++) {
	    sprintf(term, "+c%d*d%d", i, i);
	    err = add_to_spec(fnstr, term);
	}
    }

    if (!err) {
	if (coeff == NULL) {
	    nls_kickstart(amod, dset, b0, by1);
	}

#if AINIT_DEBUG
	fprintf(stderr, "initting using NLS spec:\n %s\n", fnstr);
	for (i=0; i<nparam; i++) {
	    fprintf(stderr, "initial NLS b[%d] = %g (%s)\n",
		    i, parms[i], pnames[i]);
	}
#endif

	err = nlspec_set_regression_function(spec, fnstr, dset);
    }

    if (!err) {
	set_auxiliary_scalars();
	err = aux_nlspec_add_param_list(spec, nparam, parms, pnames);
	if (!err) {
	    *amod = model_from_nlspec(spec, dset, nlsopt, prn);
	    err = amod->errcode;
#if AINIT_DEBUG
	    if (!err) {
		printmodel(amod, dset, OPT_NONE, prn);
	    }
#endif
	}
	unset_auxiliary_scalars();
    }

 bailout:

    nlspec_destroy(spec);
    free(parms);
    strings_array_free(pnames, nparam);

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

int ar_arma_init (double *coeff, const DATASET *dset,
		  arma_info *ainfo, MODEL *pmod)
{
    PRN *prn = ainfo->prn;
    int *list = ainfo->alist;
    int nmixed = ainfo->np * ainfo->P;
    int ptotal = ainfo->np + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    DATASET *aset = NULL;
    int *arlist = NULL;
    MODEL armod;
    int narmax, nonlin = 0;
    int i, err = 0;

#if AINIT_DEBUG
    fprintf(stderr, "ar_arma_init: dset->t1=%d, dset->t2=%d (dset->n=%d);\n"
	    " ainfo->t1=%d, ainfo->t2=%d, ",
	    dset->t1, dset->t2, dset->n, ainfo->t1, ainfo->t2);
    fprintf(stderr, "nmixed = %d, ptotal = %d, ifc = %d, nexo = %d\n", 
	    nmixed, ptotal, ainfo->ifc, ainfo->nexo);
#endif

    if (ptotal == 0 && ainfo->nexo == 0 && !ainfo->ifc) {
	/* special case of pure MA model */
	for (i=0; i<ainfo->nq + ainfo->Q; i++) {
	    coeff[i] = 0.0001; 
	} 
	pprintf(ainfo->prn, "\n%s: %s\n\n", _("ARMA initialization"), 
		_("small MA values"));
	return 0;
    }

    gretl_model_init(&armod, dset); 

    narmax = arma_exact_ml(ainfo) ? ainfo->nexo : 0;
    if (narmax > 0 && ptotal > 0) {
	/* ARMAX-induced lags of exog vars */
	av += ainfo->nexo * ptotal;
    } 

    if (arma_exact_ml(ainfo) && ainfo->ifc) {
	maybe_set_yscale(ainfo);
    }

    aset = create_auxiliary_dataset(av, ainfo->fullT, 0);
    if (aset == NULL) {
	return E_ALLOC;
    }

    if (ptotal > 0 && (narmax > 0 || nmixed > 0)) {
	/* we'll have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	arlist = make_ar_ols_list(ainfo, av);
    }

    /* build temporary dataset, dset -> aset */
    arma_init_build_dataset(ainfo, ptotal, narmax, list,
			    dset, aset, nonlin);

    if (nonlin) {
	PRN *dprn = NULL;

#if AINIT_DEBUG
	fprintf(stderr, "arma:_init_by_ls: doing NLS\n");
	dprn = prn;
#endif
	err = arma_get_nls_model(&armod, ainfo, narmax, NULL, aset,
				 dprn);
    } else {
#if AINIT_DEBUG
	printlist(arlist, "'arlist' in ar_arma_init (OLS)");
#endif
	armod = lsq(arlist, aset, OLS, OPT_A | OPT_Z);
	err = armod.errcode;
    }

#if AINIT_DEBUG
    if (!err) {
	pputs(prn, "\n*** armod, in ar_arma_init\n");
	printmodel(&armod, aset, OPT_NONE, prn);
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
    destroy_dataset(aset);

    return err;
}

int arma_by_ls (const double *coeff, const DATASET *dset,
		arma_info *ainfo, MODEL *pmod)
{
    PRN *prn = ainfo->prn;
    int *list = ainfo->alist;
    int nmixed = ainfo->np * ainfo->P;
    int ptotal = ainfo->np + ainfo->P + nmixed;
    int av = ptotal + ainfo->nexo + 2;
    DATASET *aset = NULL;
    int *arlist = NULL;
    int nonlin = 0;

    aset = create_auxiliary_dataset(av, ainfo->T, 0);
    if (aset == NULL) {
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
			    dset, aset, nonlin);

    if (nonlin) {
	pmod->errcode = arma_get_nls_model(pmod, ainfo, 0, coeff, aset,
					   prn);
    } else {
	*pmod = lsq(arlist, aset, OLS, OPT_A | OPT_Z);
    }

    /* clean up */
    free(arlist);
    destroy_dataset(aset);

    if (!pmod->errcode && pmod->full_n < dset->n) {
	/* the model series are short */
	double *uhat = malloc(dset->n * sizeof *uhat);
	double *yhat = malloc(dset->n * sizeof *yhat);
	int s, t;

	if (uhat == NULL || yhat == NULL) {
	    free(uhat);
	    free(yhat);
	    pmod->errcode = E_ALLOC;
	} else {
	    for (t=0; t<dset->n; t++) {
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
