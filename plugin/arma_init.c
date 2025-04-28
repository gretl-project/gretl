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
#include "libset.h"
#include "arma_priv.h"

#define AINIT_DEBUG 0

/* Given an estimate of the ARMA constant via OLS, convert to the form
   wanted for initializing the Kalman filter.  Note: the @b array
   goes: const, phi, Phi, theta, Theta, beta.
*/

void transform_arma_const (double *b, arma_info *ainfo)
{
    const double *phi = b + 1;
    const double *Phi = phi + ainfo->np;
    double narfac = 1.0; /* nonseasonal AR factor */
    double sarfac = 1.0; /* seasonal AR factor */
    int i, k = 0;

    if (ainfo->np == 0 && ainfo->P == 0) {
	/* nothing to be done */
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

#if AINIT_DEBUG
    fprintf(stderr, "transform_arma_const: revised = %g (ns=%g, s=%g)\n",
	    b[0], narfac, sarfac);
#endif
}

static int init_transform_const (arma_info *ainfo)
{
    return ainfo->ifc && arma_exact_ml(ainfo);
}

void maybe_set_yscale (arma_info *ainfo)
{
    double ybar, sdy;
    int err;

    if (arima_levels(ainfo)) {
	/* not sure about this clause! */
	ybar = gretl_mean(ainfo->t1, ainfo->t2, ainfo->y);
	if (fabs(ybar) > 250) {
	    set_arma_avg_ll(ainfo);
	}
	return;
    }

    err = gretl_moments(ainfo->t1, ainfo->t2, ainfo->y,
			NULL, &ybar, &sdy, NULL, NULL, 1);

    if (!err && sdy > 0) {
	/* try to catch cases where (a) the overall scale is
	   too big or (b) the coefficient of variation is
	   too small: in such cases set up conversion to
	   (y - ybar) / sdy + 1 = (y - (ybar - sdy)) * 1/sdy.
	*/
	double abs_ybar = fabs(ybar);
	double hi = 200, lo = 0.01;
	// double hi = 1.25, lo = 0.75;

	if (abs_ybar > hi || abs_ybar < lo || sdy/abs_ybar < lo) {
	    ainfo->yshift = ybar - sdy; /* subtract */
	    ainfo->yscale = 1 / sdy;    /* multiply */
#if 0
	    fprintf(stderr, "scale: ybar %g, sdy %g: subtract %g, mul by %g\n",
		    ybar, sdy, ainfo->yshift, ainfo->yscale);
#endif
	}
    }

    if (!err && ainfo->prn != NULL && ainfo->yscale != 1.0) {
	pputc(ainfo->prn, '\n');
	pprintf(ainfo->prn, _("Shifting y by %g, scaling by %g\n"),
		ainfo->yshift, ainfo->yscale);
    }
}

#define apply_yscaling(a,x) (arma_exact_ml(a) && !na(x))

#define HR_MINLAGS 16

/* complex inversion of @z */

static gretl_matrix *cinv (gretl_matrix *z)
{
    gretl_matrix *tmp, *ret = NULL;
    int n = z->rows;
    int i, err = 0;

    tmp = gretl_zero_matrix_new(n, 2);
    for (i=0; i<n; i++) {
	tmp->val[i] = 1.0;
    }
    ret = gretl_matrix_complex_divide(tmp, z, 1, &err);
    gretl_matrix_free(tmp);

    return ret;
}

static void copy_row (gretl_matrix *targ, int it,
		      const gretl_matrix *src, int is,
		      int neg)
{
    double d;
    int j;

    for (j=0; j<src->cols; j++) {
	d = gretl_matrix_get(src, is, j);
	gretl_matrix_set(targ, it, j, neg ? -d : d);
    }
}

/* computes a polynomial from its roots: @r is assumed
   to be n x 2 (complex)
*/

static gretl_matrix *pol_from_roots (const gretl_matrix *r)
{
    gretl_matrix *tmp, *ret = NULL;
    int n = r->rows;
    int err = 0;

    tmp = gretl_matrix_alloc(1, 2);

    if (n == 0) {
	tmp->val[0] = 1;
	tmp->val[1] = 0;
	ret = tmp;
    } else {
	copy_row(tmp, 0, r, n-1, 0);
	if (tmp->val[0] == 0 && tmp->val[1] == 0) {
	    tmp->val[0] = tmp->val[1] = NADBL;
	    ret = tmp;
        } else {
	    gretl_matrix *ix = cinv(tmp);
	    int i;

            if (n == 1) {
		/* hansl: ret = {1,0} | -ix */
		ret = gretl_zero_matrix_new(ix->rows + 1, 2);
		ret->val[0] = 1;
		for (i=0; i<ix->rows; i++) {
		    copy_row(ret, i+1, ix, i, 1);
		}
            } else {
		gretl_matrix *rslice; /* hansl: = r[1:n-1,] */
		gretl_matrix *tmp1, *tmp2;
		double d0, d1;

		rslice = gretl_matrix_alloc(n-1, 2);
		for (i=0; i<rslice->rows; i++) {
		    copy_row(rslice, i, r, i, 0);
		}
		gretl_matrix_free(tmp);
                tmp = pol_from_roots(rslice);
		/* hansl: ret = tmp | {0,0} */
		ret = gretl_zero_matrix_new(tmp->rows + 1, 2);
		for (i=0; i<tmp->rows; i++) {
		    copy_row(ret, i, tmp, i, 0);
		}
		/* hansl: ix = transp(mshape(ix, 2, n)) */
		tmp1 = gretl_matrix_shape(ix, 2, n, &err);
		gretl_matrix_transpose_in_place(tmp1);
		/* hansl: tmp = cmult(tmp1, -ix) */
		gretl_matrix_multiply_by_scalar(tmp1, -1.0);
		tmp2 = gretl_matrix_complex_multiply(tmp, tmp1, 1, &err);
		/* hansl: ret[2:,] += tmp */
		for (i=1; i<ret->rows; i++) {
		    d0 = gretl_matrix_get(ret, i, 0);
		    d0 += gretl_matrix_get(tmp2, i-1, 0);
		    d1 = gretl_matrix_get(ret, i, 1);
		    d1 += gretl_matrix_get(tmp2, i-1, 1);
		    gretl_matrix_set(ret, i, 0, d0);
		    gretl_matrix_set(ret, i, 1, d1);
		}
		gretl_matrix_free(tmp1);
		gretl_matrix_free(tmp2);
		gretl_matrix_free(rslice);
	    }
	    gretl_matrix_free(ix);
	}
    }

    if (ret == tmp) {
	tmp = NULL;
    }
    gretl_matrix_free(tmp);

    return ret;
}

/* for non-seasonal "gappy" coeff vector: expand
   by inserting zeros as needed */

static gretl_matrix *poly_from_coeff (const double *coeff,
				      const char *mask,
				      int n, int ar)
{
    gretl_matrix *ret;
    int i, k = 0;

    ret = gretl_zero_matrix_new(n + 1, 1);
    ret->val[0] = 1.0;

    for (i=0; i<n; i++) {
        if (mask[i] == '1') {
	    ret->val[i+1] = ar ? -coeff[k] : coeff[k];
	    k++;
        }
    }

    return ret;
}

/* checks if the polynomial given by @coeff is fundamental,
   and modifies it if that is not the case.
*/

int flip_poly (double *coeff, arma_info *ainfo,
	       int ar, int seasonal)
{
    gretl_matrix *tmp = NULL;
    gretl_matrix *r = NULL;
    const char *mask;
    double re, im;
    int n, n_inside = 0;
    int i, err = 0;

    if (ar) {
	n = seasonal ? ainfo->P : ainfo->p;
	mask = seasonal ? NULL : ainfo->pmask;
    } else {
	n = seasonal ? ainfo->Q : ainfo->q;
	mask = seasonal ? NULL : ainfo->qmask;
    }

    if (mask == NULL) {
	/* no expansion needed */
	tmp = gretl_matrix_alloc(n + 1, 1);
	tmp->val[0] = 1.0;
	for (i=0; i<n; i++) {
	    tmp->val[i+1] = ar ? -coeff[i] : coeff[i];
	}
    } else {
	/* expand to handle gappiness */
	tmp = poly_from_coeff(coeff, mask, n, ar);
    }

    /* force legacy form of complex output, for now */
    r = gretl_matrix_polroots(tmp, 1, 1, &err);

    if (err) {
	goto bailout;
    }

    gretl_matrix_zero(tmp);
    for (i=0; i<r->rows; i++) {
	re = gretl_matrix_get(r, i, 0);
	im = gretl_matrix_get(r, i, 1);
	if (re*re + im*im < 1.0) {
	    /* record row reference */
	    tmp->val[i] = 1;
	    n_inside++;
	}
    }

    if (n_inside > 0) {
	gretl_matrix *rfix, *ifix;
	double ci1;
	int k = 0;

	/* compose sub-matrix */
	rfix = gretl_matrix_alloc(n_inside, 2);
	for (i=0; i<r->rows; i++) {
	    if (tmp->val[i] == 1) {
		copy_row(rfix, k++, r, i, 0);
	    }
	}
	/* complex inversion */
	ifix = cinv(rfix);
	/* replace the inverted portion of r */
	k = 0;
	for (i=0; i<r->rows; i++) {
	    if (tmp->val[i] == 1) {
		copy_row(r, i, ifix, k++, 0);
	    }
	}
	gretl_matrix_free(tmp);
        tmp = pol_from_roots(r);
	if (mask != NULL) {
	    /* shrink to coeff */
	    k = 0;
	    for (i=0; i<n; i++) {
		if (mask[i] == '1') {
		    ci1 = tmp->val[i+1];
		    coeff[k++] = ar ? -ci1 : ci1;
		}
	    }
	} else {
	    /* just copy to coeff */
	    for (i=0; i<n; i++) {
		ci1 = tmp->val[i+1];
		coeff[i] = ar ? -ci1 : ci1;
	    }
	}
	gretl_matrix_free(rfix);
	gretl_matrix_free(ifix);
    }

 bailout:

    gretl_matrix_free(r);
    gretl_matrix_free(tmp);

    return err;
}

/* @pmod->coeff contains coefficients from step 2 of
   the H-R procedure, in the order: intercept, exogenous
   vars, phi, Phi, theta, Theta (in each case, if present).
   The array @b has to be filled in the order: intercept,
   phi, Phi, theta, Theta, exogenous vars.
*/

static int hr_transcribe_coeffs (arma_info *ainfo,
				 MODEL *pmod, double *b)
{
    double *phi = NULL;
    double *Phi = NULL;
    double *theta = NULL;
    double *Theta = NULL;
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

    phi = b + k;
    for (i=0; i<ainfo->p; i++) {
	if (AR_included(ainfo, i)) {
	    b[k++] = pmod->coeff[j++];
	}
    }

    Phi = b + k;
    for (i=0; i<ainfo->P; i++) {
	b[k++] = pmod->coeff[j];
	j += ainfo->np + 1; /* assumes ainfo->p < pd */
    }

    theta = b + k;
    for (i=0; i<ainfo->q; i++) {
	if (MA_included(ainfo, i)) {
	    b[k++] = pmod->coeff[j++];
	}
    }

    Theta = b + k;
    for (i=0; i<ainfo->Q; i++) {
	b[k++] = pmod->coeff[j];
	j += ainfo->nq + 1; /* assumes ainfo->q < pd */
    }

    j = ainfo->ifc;
    for (i=0; i<ainfo->nexo; i++) {
	b[k++] = pmod->coeff[j++];
    }

    if (ainfo->p > 0) {
	flip_poly(phi, ainfo, 1, 0);
    }
    if (ainfo->P > 0) {
	flip_poly(Phi, ainfo, 1, 1);
    }
    if (ainfo->q > 0) {
	flip_poly(theta, ainfo, 0, 0);
    }
    if (ainfo->Q > 0) {
	flip_poly(Theta, ainfo, 0, 1);
    }

    return err;
}

static int pre_sample_count (arma_info *ainfo,
			     const double *y,
			     int maxlag)
{
    int t, n = 0;

    for (t=ainfo->t1-1; t>=0; t--) {
	if (!na(y[t])) {
	    n++;
	    if (n == maxlag) {
		break;
	    }
	} else {
	    break;
	}
    }

    return n;
}

static double *prescale_y (double *y, arma_info *ainfo, int n)
{
    double *ys = copyvec(y, n);
    int t;

    for (t=0; t<n; t++) {
	if (!na(ys[t])) {
	    ys[t] -= ainfo->yshift;
	    ys[t] *= ainfo->yscale;
	}
    }

    return ys;
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
    double *y, *dest, *src;
    DATASET *aset = NULL;
    MODEL armod;
    int *hrlist = NULL;
    char *done = NULL;
    int maxp2p, maxp2q, p1lags;
    int nv, nv1, nv2;
    int np2p, np2q;
    int m, xpos, pos, mapos;
    size_t datalen;
    int free_y = 0;
    int i, j, T, t1;
    int err = 0;

    /* the dependent variable (of length dset->n) */
    if (arma_xdiff(ainfo)) {
	/* for initialization, use the level of y */
	y = (double *) dset->Z[ainfo->yno];
    } else {
	y = ainfo->y;
    }

    /* the greatest AR lag-length in pass 2 */
    maxp2p = ainfo->p + ainfo->pd * ainfo->P;
    /* the greatest MA lag-length in pass 2 */
    maxp2q = ainfo->q + ainfo->pd * ainfo->Q;

    /* the actual number of lags in pass 2 */
    np2p = ainfo->np + ainfo->P + ainfo->p * ainfo->P;
    np2q = ainfo->nq + ainfo->Q + ainfo->q * ainfo->Q;
    nv2 = 2 + ainfo->nexo + np2p + np2q;

    /* do we need more AR lags for H-R? */
    p1lags = maxp2p < HR_MINLAGS ? HR_MINLAGS : maxp2p;
    nv1 = 2 + ainfo->nexo + p1lags;

    /* how many variables do we need to allocate for? */
    nv = nv1 > nv2 ? nv1 : nv2;

    /* and how many observations? */
    T = ainfo->T - p1lags;
    if (ainfo->t1 > 0) {
	/* use non-missing pre-sample obs? */
	T += pre_sample_count(ainfo, y, p1lags);
    }
    /* benchmark position for reading from dset */
    t1 = ainfo->t1 + ainfo->T - T;

    /* allocate sufficient storage for both passes */
    aset = create_auxiliary_dataset(nv, T, 0);
    if (aset == NULL) {
	return E_ALLOC;
    }

#if AINIT_DEBUG
    fprintf(stderr, "hr_arma_init: dataset allocated: %d vars, %d obs\n",
	    nv, T);
#endif

    if (ainfo->yscale != 1.0) {
	y = prescale_y(y, ainfo, dset->n);
	if (y == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
	    free_y = 1;
	}
    }

    /* in case we fail before estimating a model */
    gretl_model_init(&armod, dset);

    /* regression list */
    hrlist = gretl_list_new(nv);
    if (hrlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    } else {
	hrlist[1] = 1;
	hrlist[2] = 0;
	for (i=2; i<nv; i++) {
	    hrlist[i+1] = i;
	}
    }

    /* adjust the list for pass 1 */
    hrlist[0] = nv1;

    /* recorder array, etc. */
    done = calloc(p1lags, 1);
    datalen = T * sizeof(double);

    /* dependent var */
    strcpy(aset->varname[1], "y");
    memcpy(aset->Z[1], y + t1, datalen);
    /* exogenous vars */
    pos = 2;
    for (i=0; i<ainfo->nexo; i++) {
	sprintf(aset->varname[pos], "x%d", i);
	xpos = ainfo->xlist[i+1];
	memcpy(aset->Z[pos], dset->Z[xpos] + t1, datalen);
	pos++;
    }
    /* pass2 non-seasonal AR lags first */
    for (i=1; i<=ainfo->p; i++) {
	if (AR_included(ainfo, i-1)) {
	    sprintf(aset->varname[pos], "y_%d", i);
	    memcpy(aset->Z[pos], y + t1 - i, datalen);
	    done[i-1] = 1;
	    pos++;
	}
    }
    /* then pass2 seasonal AR lags */
    for (j=1; j<=ainfo->P; j++) {
	m = j * ainfo->pd;
	sprintf(aset->varname[pos], "y_%d", m);
	memcpy(aset->Z[pos], y + t1 - m, datalen);
	done[m-1] = 1;
	pos++;
    }
    /* then pass2 AR interactions */
    for (j=1; j<=ainfo->P; j++) {
	for (i=1; i<=ainfo->p; i++) {
	    if (AR_included(ainfo, i-1)) {
		m = j * ainfo->pd + i;
		sprintf(aset->varname[pos], "y_%d", m);
		memcpy(aset->Z[pos], y + t1 - m, datalen);
		done[m-1] = 1;
		pos++;
	    }
	}
    }
    mapos = pos; /* insertion point for pass2 MA lags */
    /* then any "extra" AR lags for pass 1 only */
    for (i=1; i<=p1lags; i++) {
	if (!done[i-1]) {
	    sprintf(aset->varname[pos], "y_%d", i);
	    memcpy(aset->Z[pos], y + t1 - i, datalen);
	    pos++;
	}
    }

    /* pass 1 estimation (FIXME constant?) */
    armod = lsq(hrlist, aset, OLS, OPT_A);
    if (armod.errcode) {
	err = armod.errcode;
	goto bailout;
    }

#if AINIT_DEBUG
    fprintf(stderr, "pass1 model: t1=%d, t2=%d, nobs=%d, ncoeff=%d, dfd = %d\n",
	    armod.t1, armod.t2, armod.nobs, armod.ncoeff, armod.dfd);
#endif

    /* revise hrlist if needed (FIXME constant?) */
    hrlist[0] = nv2;
    /* position for insertion of MA terms: at least some of
       these will likely overwrite AR terms that were wanted
       only for pass 1
    */
    pos = mapos;

    /* pass 2 sample: skip the leading observations for
       which we don't have lagged residuals from pass 1
    */
    aset->t1 = maxp2q;
    datalen = (T - maxp2q) * sizeof(double);

    for (i=1; i<=ainfo->q; i++) {
	if (MA_included(ainfo, i-1)) {
	    sprintf(aset->varname[pos], "e_%d", i);
	    dest = aset->Z[pos] + aset->t1;
	    src = armod.uhat + aset->t1 - i;
	    memcpy(dest, src, datalen);
	    pos++;
	}
    }
    for (j=1; j<=ainfo->Q; j++) {
	m = j * ainfo->pd;
	sprintf(aset->varname[pos], "e_%d", m);
	dest = aset->Z[pos] + aset->t1;
	src = armod.uhat + aset->t1 - m;
	memcpy(dest, src, datalen);
	pos++;
    }
    for (j=1; j<=ainfo->Q; j++) {
	for (i=1; i<=ainfo->q; i++) {
	    if (MA_included(ainfo, i-1)) {
		m = j * ainfo->pd + i;
		sprintf(aset->varname[pos], "e_%d", m);
		dest = aset->Z[pos] + aset->t1;
		src = armod.uhat + aset->t1 - m;
		memcpy(dest, src, datalen);
		pos++;
	    }
	}
    }

    /* pass 2 estimation */
    clear_model(&armod);
    armod = lsq(hrlist, aset, OLS, OPT_A);

    if (armod.errcode) {
	err = armod.errcode;
    } else {
#if AINIT_DEBUG
	PRN *modprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	printmodel(&armod, aset, OPT_S, modprn);
	gretl_print_destroy(modprn);
#endif
	err = hr_transcribe_coeffs(ainfo, &armod, coeff);
	if (!err && ainfo->nexo == 0 && init_transform_const(ainfo)) {
	    transform_arma_const(coeff, ainfo);
	}
    }

#if AINIT_DEBUG
    if (!err) {
	fprintf(stderr, "HR init (nobs=%d):\n", armod.nobs);
	for (i=0; i<ainfo->nc; i++) {
	    fprintf(stderr, "coeff[%d] = %g\n", i, coeff[i]);
	}
    }
#endif

 bailout:

    free(hrlist);
    free(done);
    destroy_dataset(aset);
    clear_model(&armod);
    if (free_y) {
	free(y);
    }

    if (!err) {
	ainfo->init = INI_HR;
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
		  arma_info *ainfo)
{
    int ok = hr_df_check(ainfo, dset);
    int err = 0;

    if (ok) {
	err = real_hr_arma_init(coeff, dset, ainfo, ainfo->prn);
    }

#if AINIT_DEBUG
    if (ainfo->init) {
	fputs("*** hr_arma_init OK\n", stderr);
    } else {
	fputs("*** hr_arma_init failed, will try ar_arma_init\n", stderr);
    }
#endif

    return err;
}

static double get_y_mean (arma_info *ainfo)
{
    double ysum = 0.0;
    int t, T = 0;

    for (t=ainfo->t1; t<=ainfo->t2; t++) {
	if (!na(ainfo->y[t])) {
	    if (ainfo->yscale != 1.0) {
		ysum += (ainfo->y[t] - ainfo->yshift) * ainfo->yscale;
	    } else {
		ysum += ainfo->y[t];
	    }
	    T++;
	}
    }

    return ysum / T;
}

#define MA_SMALL  0.0001

/* Transcribe coeffs from the OLS or NLS model used for initializing,
   into the array @b that will be passed to the maximizer. While
   we're at it, check that the AR polynomial(s) are in the
   stationary zone.
*/

static void ar_init_transcribe_coeffs (arma_info *ainfo,
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
	b[q0 + i] = MA_SMALL;
    }

    /* stationarity check/fix */
    if (ainfo->p > 0) {
	flip_poly(b + ainfo->ifc, ainfo, 1, 0);
    }
    if (ainfo->P > 0) {
	flip_poly(b + ainfo->ifc + ainfo->np, ainfo, 1, 1);
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

/* X, if non-NULL, holds the differenced regressors */

static double get_xti (const DATASET *dset,
		       int i, int t,
		       const int *xlist,
		       const gretl_matrix *X)
{
    if (X != NULL) {
	return gretl_matrix_get(X, t, i);
    } else {
        return dset->Z[xlist[i+1]][t];
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
    const int *xlist = ainfo->xlist;
    int i, j, k, kx, ky;
    int t, s, k0 = 2;
    int undo_diff = 0;
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

    for (t=0; t<aset->n; t++) {
	int realt = t + ainfo->t1;
	int miss = 0;

	if (apply_yscaling(ainfo, y[realt])) {
	    aZ[1][t] = (y[realt] - ainfo->yshift) * ainfo->yscale;
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
		if (apply_yscaling(ainfo, y[s])) {
		    aZ[k][t] -= ainfo->yshift;
		    aZ[k][t] *= ainfo->yscale;
		}
		k++;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = get_xti(dset, j, s, xlist, X);
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
		if (apply_yscaling(ainfo, y[s])) {
		    aZ[k][t] -= ainfo->yshift;
		    aZ[k][t] *= ainfo->yscale;
		}
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = get_xti(dset, k, s, xlist, X);
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
		    if (apply_yscaling(ainfo, y[s])) {
			aZ[ky][t] -= ainfo->yshift;
			aZ[ky][t] *= ainfo->yscale;
		    }
		    ky++;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = get_xti(dset, k, s, xlist, X);
		    }
		}
	    }
	}

	kx = ptotal + k0;

	for (i=0; i<ainfo->nexo; i++) {
	    aZ[kx++][t] = get_xti(dset, i, realt, xlist, X);
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

#if AINIT_DEBUG > 1
    PRN *eprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

    if (eprn != NULL) {
	printdata(NULL, NULL, aset, OPT_O, eprn);
	gretl_print_destroy(eprn);
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
	strcat(chunk, "_b0");
    }

    err = add_to_spec(spec, chunk);

    for (i=0; i<narmax && !err; i++) {
	if (ainfo->ifc || i > 0) {
	    err += add_to_spec(spec, "+");
	}
	sprintf(chunk, "_b%d*x%d_%d", i+1, i+1, lag);
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
	strcpy(pnames[k++], "_b0");
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
	    sprintf(pnames[k++], "_phi%d", i+1);
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
	sprintf(pnames[k++], "_Phi%d", i+1);
    }

    for (i=0; i<ainfo->nexo; i++) {
	if (coeff != NULL) {
	    parms[k] = coeff[k];
	}
	sprintf(pnames[k++], "_b%d", i+1);
    }

    if (ainfo->misslist != NULL) {
	for (i=1; i<=ainfo->misslist[0]; i++) {
	    j = ainfo->misslist[i];
	    parms[k] = dset->Z[1][j];
	    sprintf(pnames[k++], "_c%d", i);
	}
    }

    /* construct NLS specification string */

    strcpy(fnstr, "y=");

    if (ainfo->ifc) {
	strcat(fnstr, "_b0");
    } else {
	strcat(fnstr, "0");
    }

    for (i=0; i<ainfo->p && !err; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(term, "+_phi%d*", lag);
	    err = add_to_spec(fnstr, term);
	    if (!err) {
		err = y_Xb_at_lag(fnstr, ainfo, narmax, lag);
	    }
	}
    }

    for (j=0; j<ainfo->P && !err; j++) {
	sprintf(term, "+_Phi%d*", j+1);
	strcat(fnstr, term);
	lag = (j + 1) * ainfo->pd;
	y_Xb_at_lag(fnstr, ainfo, narmax, lag);
	for (i=0; i<ainfo->p; i++) {
	    if (AR_included(ainfo, i)) {
		sprintf(term, "-_phi%d*_Phi%d*", i+1, j+1);
		err = add_to_spec(fnstr, term);
		if (!err) {
		    lag = (j+1) * ainfo->pd + (i+1);
		    y_Xb_at_lag(fnstr, ainfo, narmax, lag);
		}
	    }
	}
    }

    for (i=0; i<ainfo->nexo && !err; i++) {
	sprintf(term, "+_b%d*x%d", i+1, i+1);
	err = add_to_spec(fnstr, term);
    }

    if (!err && ainfo->misslist != NULL) {
	for (i=1; i<=ainfo->misslist[0]; i++) {
	    sprintf(term, "+_c%d*d%d", i, i);
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
	double save_tol = libset_get_double(NLS_TOLER);

	libset_set_double(NLS_TOLER, 1.0e-5);
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
	libset_set_double(NLS_TOLER, save_tol);
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

/* Apply least squares to get initial values for the AR coefficients,
   either OLS or NLS.  We use NLS if there is nonlinearity due to
   either (a) the presence of both a seasonal and a non-seasonal AR
   component or (b) the presence of exogenous variables in the context
   of a non-zero AR order, where estimation will be via exact ML.

   In this initialization any MA coefficients are simply set to
   "near-zero" (MA_SMALL).
*/

int ar_arma_init (double *coeff, const DATASET *dset,
		  arma_info *ainfo, MODEL *pmod,
		  gretlopt opt)
{
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
	    coeff[i] = MA_SMALL;
	}
	ainfo->init = INI_SMALL;
	return 0;
    }

    gretl_model_init(&armod, dset);

    narmax = arma_exact_ml(ainfo) ? ainfo->nexo : 0;
    if (narmax > 0 && ptotal > 0) {
	/* ARMAX-induced lags of exog vars */
	av += ainfo->nexo * ptotal;
    }

    if (ptotal == 0 && ainfo->nexo == 0 && ainfo->ifc) {
	/* straight MA model with constant */
	coeff[0] = get_y_mean(ainfo);
	for (i=1; i<=ainfo->nq + ainfo->Q; i++) {
	    coeff[i] = MA_SMALL;
	}
	ainfo->init = INI_SMALL;
	return 0;
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
	dprn = ainfo->prn;
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
    if (err) {
	fprintf(stderr, "LS init: armod.errcode = %d\n", err);
    }
#endif

    if (!err) {
	ar_init_transcribe_coeffs(ainfo, &armod, coeff);
    }

    /* handle the case where we need to translate from an
       estimate of the regression constant to the
       unconditional mean of y_t
    */
    if (!err && (!nonlin || ainfo->nexo == 0) &&
	init_transform_const(ainfo)) {
	transform_arma_const(coeff, ainfo);
    }

    if (!err) {
	ainfo->init = nonlin ? INI_NLS : INI_OLS;
    }

    /* clean up */
    clear_model(&armod);
    destroy_dataset(aset);
    free(arlist);

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

    aset = create_auxiliary_dataset(av, ainfo->fullT, 0);
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
	gretlopt opt = OPT_A | OPT_Z;

	if (ainfo->nc == 0) {
	    opt |= OPT_U;
	}
	*pmod = lsq(arlist, aset, OLS, opt);
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
