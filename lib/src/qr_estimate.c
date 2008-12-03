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
#include "qr_estimate.h"
#include "gretl_matrix.h"
#include "matrix_extra.h"
#include "libset.h"
#include "gretl_panel.h"
#include "estim_private.h"

#include "gretl_f2c.h"
#include "clapack_double.h"

#define QR_RCOND_MIN 1e-15 /* experiment with this? */
#define ESSZERO      1e-22 /* SSR less than this counts as zero */

enum {
    VCV_SIMPLE,
    VCV_ROBUST,
    VCV_XPX
};

/* General note: in fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are contiguous.
*/

static double qr_get_tss (MODEL *pmod, const double *y, int *ifc,
			  int *yconst)
{
    int pwe = gretl_model_get_int(pmod, "pwe");
    double y0 = 0.0, ymean = 0.0;
    double x, tss = 0.0;
    double ctss = 0.0;
    int t;

    if (*ifc == 0) {
	*ifc = check_for_effective_const(pmod, y);
    }

    *yconst = 1;

    if (pmod->rho != 0.0) {
	double ry, d;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    ry = y[t];
	    if (t == pmod->t1 && pwe) {
		ry *= sqrt(1.0 - pmod->rho * pmod->rho);
	    } else {
		ry -= pmod->rho * y[t-1];
	    }
	    ymean += ry;
	}
	ymean /= pmod->nobs;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    ry = y[t];
	    if (t == pmod->t1 && pwe) {
		ry *= sqrt(1.0 - pmod->rho * pmod->rho);
	    } else {
		ry -= pmod->rho * y[t-1];
	    }
	    if (t == pmod->t1) {
		y0 = ry;
	    } else if (ry != y0) {
		*yconst = 0;
	    }
	    d = ry - ymean;
	    if (*ifc) {
		tss += d * d;
	    } else {
		tss += ry * ry;
		ctss += d * d;
	    }
	}
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->yhat[t])) {
		ymean += y[t];
	    }
	}
	ymean /= pmod->nobs;

	y0 = y[pmod->t1];

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->yhat[t])) {
		if (y[t] != y0) {
		    *yconst = 0;
		}
		x = y[t] - ymean;
		if (*ifc) {
		    tss += x * x;
		} else {
		    tss += y[t] * y[t];
		    ctss += x * x;
		}
	    }
	}
    } 

    if (!*ifc && ctss > 0) {
	double cR2 = 1 - (pmod->ess / ctss);

	gretl_model_set_double(pmod, "centered-R2", cR2);
    }

    return tss;
}

static void qr_compute_stats (MODEL *pmod, const double *y, int n,
			      gretlopt opt)
{
    int yconst, ifc = pmod->ifc;

    pmod->tss = qr_get_tss(pmod, y, &ifc, &yconst);

    if (yconst && pmod->dfd > 0) {
	double y0 = y[pmod->t1];
    
	if (y0 > 0) {
	    double tss = pmod->nobs * y0 * y0;

	    pmod->rsq = 1 - (pmod->ess / tss);
	    gretl_model_set_int(pmod, "uncentered", 1);
	}
    } else if (pmod->dfd > 0) {
	double den = pmod->tss * pmod->dfd;

	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	pmod->adjrsq = 1.0 - (pmod->ess * (n - 1) / den);
    } else {
	pmod->rsq = 1.0;
    }

    if (pmod->ncoeff == 1 && pmod->ifc) {
	pmod->fstt = NADBL;
	return;
    }

    if (pmod->dfd > 0 && pmod->dfn > 0) {
	if (opt & OPT_R) {
	    pmod->fstt = robust_omit_F(NULL, pmod);
	} else if (pmod->rsq == 1.0) {
	    pmod->fstt = NADBL;
	} else {
	    pmod->fstt = (pmod->tss - pmod->ess) * pmod->dfd / 
		(pmod->ess * pmod->dfn);
	}
    } else {
	pmod->fstt = NADBL;
    }
}

static int qr_make_vcv (MODEL *pmod, gretl_matrix *v, int flag)
{
    int k = pmod->ncoeff;
    int m = k * (k + 1) / 2;
    double x;
    int i, j, idx;

    pmod->vcv = malloc(m * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) {
	return E_ALLOC;
    }

    if (flag == VCV_XPX) {
	gretl_model_set_int(pmod, "vcv_xpx", 1);
    }

    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    idx = ijton(i, j, k);
	    x = gretl_matrix_get(v, i, j);
	    if (flag == VCV_SIMPLE) {
		x *= pmod->sigma * pmod->sigma;
	    }
	    pmod->vcv[idx] = x;
	    if (i == j) {
		pmod->sderr[i] = sqrt(x);
		if (flag == VCV_XPX) {
		    pmod->sderr[i] *= pmod->sigma;
		}
	    }
	}
    }

    return 0;
}

static void get_resids_and_SSR (MODEL *pmod, const double **Z,
				gretl_matrix *yhat, int fulln)
{
    int t, i = 0;
    int dwt = gretl_model_get_int(pmod, "wt_dummy");
    int qdiff = (pmod->rho != 0.0);
    int pwe = gretl_model_get_int(pmod, "pwe");
    int yvar = pmod->list[1];
    double y;

    if (dwt) {
	dwt = pmod->nwt;
    }

    pmod->ess = 0.0;

    if (qdiff) {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    } else {
		y = Z[yvar][t];
		if (t == pmod->t1 && pwe) {
		    y *= sqrt(1.0 - pmod->rho * pmod->rho);
		} else {
		    y -= pmod->rho * Z[yvar][t-1];
		}
		pmod->yhat[t] = yhat->val[i];
		pmod->uhat[t] = y - yhat->val[i];
		pmod->ess += pmod->uhat[t] * pmod->uhat[t];
		i++;
	    }
	}
    } else if (pmod->nwt) {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2 || model_missing(pmod, t)) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    } else {
		y = Z[yvar][t];
		if (dwt && Z[dwt][t] == 0.0) {
		    pmod->yhat[t] = NADBL;
		} else {
		    if (!dwt) {
			y *= sqrt(Z[pmod->nwt][t]);
		    }
		    pmod->yhat[t] = yhat->val[i];
		    pmod->uhat[t] = y - yhat->val[i];
		    pmod->ess += pmod->uhat[t] * pmod->uhat[t];
		    i++;
		}
	    }
	}
    } else {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2 || model_missing(pmod, t)) {
		pmod->yhat[t] = pmod->uhat[t] = NADBL;
	    } else {
		pmod->yhat[t] = yhat->val[i];
		pmod->uhat[t] = Z[yvar][t] - yhat->val[i];
		pmod->ess += pmod->uhat[t] * pmod->uhat[t];
		i++;
	    }
	}
    }

    /* if SSR is small enough, treat it as zero */
    if (fabs(pmod->ess) < ESSZERO) {
	pmod->ess = 0.0;
    } 
}

static void 
get_data_X (gretl_matrix *X, const MODEL *pmod, const double **Z)
{
    int wt = pmod->nwt;
    int i, j, t, vi;

    /* copy independent vars into matrix X */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	vi = pmod->list[i];
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!model_missing(pmod, t)) {
		if (wt) {
		    X->val[j++] = sqrt(Z[wt][t]) * Z[vi][t];
		} else {
		    X->val[j++] = Z[vi][t];
		}
	    }
	}
    }
}

static gretl_matrix *make_data_X (const MODEL *pmod, const double **Z)
{
    gretl_matrix *X;

    X = gretl_matrix_alloc(pmod->nobs, pmod->ncoeff);
    if (X != NULL) {
	get_data_X(X, pmod, Z);
    }

    return X;
}

/* Calculate W(t)-transpose * W(t-lag) */

static void wtw (gretl_matrix *wt, gretl_matrix *X, 
		 int n, int t, int lag)
{
    int i, j;
    double xi, xj;

    for (i=0; i<n; i++) {
	xi = gretl_matrix_get(X, t, i);
	for (j=0; j<n; j++) {
	    xj = gretl_matrix_get(X, t - lag, j);
	    gretl_matrix_set(wt, i, j, xi * xj);
	}
    }
}

double qs_hac_weight (double bt, int i)
{
    double di = i / bt;
    double mi = 6 * M_PI * di / 5;
    double w;

    w = 25 / (12 * M_PI * M_PI * di * di);
    w *= sin(mi) / mi - cos(mi);

    return w;
}

double hac_weight (int kern, int h, int i)
{
    double ai = fabs((double) i) / (h + 1.0);
    double w;

    if (kern == KERNEL_PARZEN) {
	if (ai <= 0.5) {
	    w = 1.0 - 6*ai*ai + 6*pow(ai, 3.0);
	} else {
	    w = 2.0 * pow(1.0 - ai, 3.0);
	}
    } else {
	/* Bartlett kernel */
	w = 1.0 - ai;
    }

    return w;
}

/* Newey and West's data-based bandwidth selection, based on
   the exposition in A. Hall, "Generalized Method of Moments"
   (Oxford, 2005), p. 82.
*/

int newey_west_bandwidth (const gretl_matrix *f, int kern, int *h, double *bt)
{
    const double cg[] = { 1.4117, 2.6614, 1.3221 };
    const double v[] = { 1, 2, 2 };
    double g, p, s0, sv;
    double *s = NULL, *c = NULL;
    int n, T, q;
    int i, j, t;
    int err = 0;

    if (f == NULL) {
	return E_ALLOC;
    }

    T = f->rows;
    q = f->cols;

    if (kern == KERNEL_BARTLETT) {
	n = (int) pow((double) T, 2.0 / 9);
    } else if (kern == KERNEL_PARZEN) {
	n = (int) pow((double) T, 4.0 / 25);
    } else {
	n = (int) pow((double) T, 2.0 / 25);
    }

    s = malloc((n + 1) * sizeof *s);
    c = malloc(T * sizeof *c);

    if (s == NULL || c == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (t=0; t<T; t++) {
	c[t] = 0.0;
	for (i=0; i<q; i++) {
	    /* equal weighting here? */
	    c[t] += gretl_matrix_get(f, t, i);
	}
    }

    for (j=0; j<=n; j++) {
	s[j] = 0.0;
	for (t=j; t<T; t++) {
	    s[j] += c[t] * c[t-j];
	}
	s[j] /= T;
    }

    s0 = s[0];
    sv = 0.0;
    
    for (j=1; j<=n; j++) {
	if (kern == KERNEL_BARTLETT) {
	    sv += 2.0 * j * s[j];
	} else {
	    sv += 2.0 * j * j * s[j];
	}
	s0 += 2 * s[j];
    }

    p = 1.0 / (2.0 * v[kern] + 1);
    g = cg[kern] * pow((sv / s0) * (sv / s0), p);

    *bt = g * pow((double) T, p);
    *h = (int) floor(*bt);

#if 0
    fprintf(stderr, "bt = %g, h = %d\n", *bt, *h);
#endif
    
 bailout:

    free(s);
    free(c);

    return err;
}

static int prewhiten_uhat (double **pu, int T, double *pa)
{
    double a, num = 0.0, den = 0.0;
    double *uw, *u = *pu;
    int sgn, t;

    uw = malloc(T * sizeof *uw);
    if (uw == NULL) {
	return E_ALLOC;
    }

    for (t=1; t<T; t++) {
	num += u[t-1] * u[t];
	den += u[t-1] * u[t-1];
    }

    a = num / den;
    sgn = (a < 0.0)? -1 : 1;
    if (fabs(a) > 0.97) {
	a = sgn * 0.97;
    }

    for (t=1; t<T; t++) {
	uw[t] = u[t] - a * u[t-1];
    } 

    *pa = a;
    *pu = uw;

    return 0;
}

/* Calculate HAC covariance matrix.  Algorithm and (basically)
   notation taken from Davidson and MacKinnon (DM), Econometric Theory
   and Methods, chapter 9.
*/

static int qr_make_hac (MODEL *pmod, const double **Z, gretl_matrix *xpxinv)
{
    gretl_matrix *vcv = NULL, *wtj = NULL, *gammaj = NULL;
    gretl_matrix *X;
    int prewhiten = libset_get_bool(PREWHITEN);
    int kern = libset_get_int(HAC_KERNEL);
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int free_uhat = 0;
    int p, j, t;
    double wj, uu;
    double a = 0, bt = 0;
    double *uhat;
    int err = 0;

    X = make_data_X(pmod, Z);
    if (X == NULL) return 1;

    /* pmod->uhat is a full-length series: we must take an offset
       into it, equal to the offset of the data on which the model
       is actually estimated.
    */
    uhat = pmod->uhat + pmod->t1;

    if (prewhiten) {
	err = prewhiten_uhat(&uhat, T, &a);
	if (err) {
	    return err;
	}
	free_uhat = 1;
    }

    vcv = gretl_matrix_alloc(k, k);
    wtj = gretl_matrix_alloc(k, k);
    gammaj = gretl_matrix_alloc(k, k);

    if (vcv == NULL || wtj == NULL || gammaj == NULL) {
	err = 1;
	goto bailout;
    }

    /* determine the bandwidth setting */

    if (data_based_hac_bandwidth()) {
	gretl_matrix u;

	u.rows = T;
	u.cols = 1;
	u.val = uhat;
	err = newey_west_bandwidth(&u, kern, &p, &bt);
	if (err) {
	    goto bailout;
	}
    } else if (kern == KERNEL_QS) {
	bt = libset_get_double(QS_BANDWIDTH);
	p = pmod->nobs - 1;
    } else {
	p = get_hac_lag(T);
    }

    gretl_model_set_int(pmod, "hac_kernel", kern);
    gretl_model_set_int(pmod, "hac_prewhiten", prewhiten);
    if (kern == KERNEL_QS) {
	gretl_model_set_double(pmod, "qs_bandwidth", bt);
    } else {
	gretl_model_set_int(pmod, "hac_lag", p);
    }

    gretl_matrix_zero(vcv);

    for (j=0; j<=p; j++) {
	/* cumulate running sum of Gamma-hat terms */
	gretl_matrix_zero(gammaj);
	for (t=j; t<T; t++) {
	    /* W(t)-transpose * W(t-j) */
	    wtw(wtj, X, k, t, j);
	    uu = uhat[t] * uhat[t-j];
	    gretl_matrix_multiply_by_scalar(wtj, uu);
	    /* DM equation (9.36), p. 363 */
	    gretl_matrix_add_to(gammaj, wtj);
	}

	if (j > 0) {
	    /* Gamma(j) = Gamma(j) + Gamma(j)-transpose */
	    gretl_matrix_add_self_transpose(gammaj);
	    if (kern == KERNEL_QS) {
		wj = qs_hac_weight(bt, j);
	    } else {
		wj = hac_weight(kern, p, j);
	    }
	    gretl_matrix_multiply_by_scalar(gammaj, wj);
	}

	/* DM equation (9.38), p. 364 */
	gretl_matrix_add_to(vcv, gammaj);
    }

    if (prewhiten) {
	/* re-color */
	gretl_matrix_divide_by_scalar(vcv, (1-a) * (1-a));
    }

    gretl_matrix_copy_values(wtj, vcv);
    gretl_matrix_qform(xpxinv, GRETL_MOD_TRANSPOSE, wtj,
		       vcv, GRETL_MOD_NONE);

    /* Transcribe vcv into triangular representation */
    err = qr_make_vcv(pmod, vcv, VCV_ROBUST);

 bailout:

    gretl_matrix_free(wtj);
    gretl_matrix_free(gammaj);
    gretl_matrix_free(vcv);
    gretl_matrix_free(X);

    if (free_uhat) {
	free(uhat);
    }

    return err;
}

/* Multiply X transpose into D, treating D as if it were
   a diagonal matrix -- although in fact it is just a vector,
   for economy of storage.  Result into R.
*/

static void do_X_prime_diag (const gretl_matrix *X,
			     const gretl_vector *D,
			     gretl_matrix *R)
{
    double x;
    int i, j;

    for (i=0; i<R->rows; i++) {
	for (j=0; j<R->cols; j++) {
	    x = gretl_matrix_get(X, j, i);
	    gretl_matrix_set(R, i, j, x * D->val[j]);
	}
    }
}

/* Heteroskedasticity-Consistent Covariance Matrices: See Davidson
   and MacKinnon, Econometric Theory and Methods, chapter 5, esp.
   page 200.  Implements HC0, HC1, HC2 and HC3.
*/

static int qr_make_hccme (MODEL *pmod, const double **Z, 
			  gretl_matrix *Q, gretl_matrix *xpxinv)
{
    gretl_matrix *X;
    gretl_matrix *diag = NULL;
    gretl_matrix *tmp1 = NULL, *tmp2 = NULL, *vcv = NULL;
    int T = pmod->nobs; 
    int k = pmod->list[0] - 1;
    int hc_version;
    int i, t;
    int err = 0;

    X = make_data_X(pmod, Z);
    if (X == NULL) return 1;

    diag = gretl_vector_from_array(pmod->uhat + pmod->t1, T,
				   GRETL_MOD_SQUARE);
    if (diag == NULL) {
	err = 1;
	goto bailout;
    }  

    tmp1 = gretl_matrix_alloc(k, T);
    tmp2 = gretl_matrix_alloc(k, k);
    vcv = gretl_matrix_alloc(k, k);
    if (tmp1 == NULL || tmp2 == NULL || vcv == NULL) {
	err = 1;
	goto bailout;
    }  

    hc_version = libset_get_int(HC_VERSION);
    gretl_model_set_int(pmod, "hc", 1);
    if (hc_version > 0) {
	gretl_model_set_int(pmod, "hc_version", hc_version);
    }

    if (hc_version == 1) {
	for (t=0; t<T; t++) {
	    diag->val[t] *= (double) T / (T - k);
	}
    } else if (hc_version > 1) {
	/* do the h_t calculations */
	for (t=0; t<T; t++) {
	    double q, ht = 0.0;

	    for (i=0; i<k; i++) {
		q = gretl_matrix_get(Q, t, i);
		ht += q * q;
	    }
	    if (hc_version == 2) {
		diag->val[t] /= (1.0 - ht);
	    } else { /* HC3 */
		diag->val[t] /= (1.0 - ht) * (1.0 - ht);
	    }
	}
    }

    do_X_prime_diag(X, diag, tmp1);
    gretl_matrix_multiply(tmp1, X, tmp2);
    gretl_matrix_qform(xpxinv, GRETL_MOD_NONE, tmp2,
		       vcv, GRETL_MOD_NONE);

    /* Transcribe vcv into triangular representation */
    err = qr_make_vcv(pmod, vcv, VCV_ROBUST);

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(diag);
    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);
    gretl_matrix_free(vcv);

    return err;
}

static int qr_dw_stats (MODEL *pmod, const double **Z,
			gretl_matrix *X, gretl_matrix *u)
{
    double DW, pv;
    int t, s, err = 0;

    get_data_X(X, pmod, Z);

    for (s=0, t=pmod->t1; t<=pmod->t2; s++, t++) {
	gretl_vector_set(u, s, pmod->uhat[t]);
    }
    
    pv = dw_pval(u, X, &DW, &err);

    if (!err) {
	pmod->dw = DW;
	gretl_model_set_double(pmod, "dw_pval", pv);
    }
    
    return err;
}

static int qr_make_regular_vcv (MODEL *pmod, gretl_matrix *v,
				gretlopt opt)
{
    int flag = (opt & OPT_X)? VCV_XPX : VCV_SIMPLE;

    return qr_make_vcv(pmod, v, flag);
}

static void get_model_data (MODEL *pmod, const double **Z, 
			    gretl_matrix *Q, gretl_matrix *y)
{
    int i, j, t;
    double x;
    int dwt = gretl_model_get_int(pmod, "wt_dummy");
    int qdiff = (pmod->rho != 0.0);
    int pwe = gretl_model_get_int(pmod, "pwe");
    double pw1 = 0.0;

    if (pwe) {
	pw1 = sqrt(1.0 - pmod->rho * pmod->rho);
    }

    if (dwt) {
	dwt = pmod->nwt;
    }

    /* copy independent vars into matrix Q */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    x = Z[pmod->list[i]][t];
	    if (dwt) {
		if (Z[dwt][t] == 0.0) continue;
	    } else if (pmod->nwt) {
		x *= sqrt(Z[pmod->nwt][t]);
	    } else if (qdiff) {
		if (pwe && t == pmod->t1) {
		    x *= pw1;
		} else {
		    x -= pmod->rho * Z[pmod->list[i]][t-1];
		}
	    }
	    Q->val[j++] = x;
	}
    }

    if (y != NULL) {
	/* copy dependent variable into y vector */
	j = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }		
	    x = Z[pmod->list[1]][t];
	    if (dwt) {
		if (Z[dwt][t] == 0.0) continue;
	    } else if (pmod->nwt) {
		x *= sqrt(Z[pmod->nwt][t]);
	    } else if (qdiff) {
		if (pwe && t == pmod->t1) {
		    x *= pw1;
		} else {
		    x -= pmod->rho * Z[pmod->list[1]][t-1];
		}
	    }
	    y->val[j++] = x;
	}
    }
}

static int 
allocate_model_arrays (MODEL *pmod, int k, int T)
{
    if (pmod->sderr == NULL) {
	pmod->sderr = malloc(k * sizeof *pmod->sderr);
    }
    if (pmod->yhat == NULL) {
	pmod->yhat = malloc(T * sizeof *pmod->yhat);
    }
    if (pmod->uhat == NULL) {
	pmod->uhat = malloc(T * sizeof *pmod->uhat);
    }

    if (pmod->sderr == NULL || pmod->yhat == NULL || pmod->uhat == NULL) {
	return 1;
    }

    return 0;
}

/* perform QR decomposition plus some additional tasks */

static int QR_decomp_plus (gretl_matrix *Q, gretl_matrix *R, int *rank)
{
    integer k = gretl_matrix_rows(R);
    int r, err;

    /* basic decomposition */
    err = gretl_matrix_QR_decomp(Q, R);
    if (err) {
	return err;
    }

    /* check rank of QR */
    r = gretl_check_QR_rank(R, &err);
    if (err) {
	return err;
    }

    if (r < k) {
	err = E_SINGULAR;
    } else {
	/* invert R */
	char uplo = 'U';
	char diag = 'N';
	integer info = 0;

	dtrtri_(&uplo, &diag, &k, R->val, &k, &info);
	if (info != 0) {
	    fprintf(stderr, "dtrtri: info = %d\n", (int) info);
	    err = 1;
	}
    }

    if (rank != NULL) {
	*rank = r;
    }

    return err;
}

#define REDEBUG 0

static void
drop_redundant_vars (MODEL *pmod, gretl_matrix *R, int rank, gretlopt opt)
{
    int *dlist = NULL;
    double d;
    int nd = R->rows - rank;
    int i, j;

    if (!(opt & OPT_A)) {
	dlist = gretl_list_new(nd);
	if (dlist != NULL) {
	    dlist[0] = 0;
	}
    }

#if REDEBUG
    printlist(pmod->list, "pmod->list, into drop_redundant_vars");
    fprintf(stderr, "rank = %d\n", rank);
#endif

    j = 2;
    nd = 0;
    for (i=0; i<R->rows; i++) {
	d = gretl_matrix_get(R, i, i);
	if (fabs(d) < R_DIAG_MIN) {
	    if (dlist != NULL) {
		dlist[0] += 1;
		dlist[dlist[0]] = pmod->list[j];
	    }
	    gretl_list_delete_at_pos(pmod->list, j--);
	    nd++;
	}
	j++;
    }

    pmod->ncoeff -= nd;
    pmod->dfd = pmod->nobs - pmod->ncoeff;
    pmod->dfn = pmod->ncoeff - pmod->ifc;

    if (dlist != NULL) {
	gretl_model_set_list_as_data(pmod, "droplist", dlist);
    }   
}

int gretl_qr_regress (MODEL *pmod, const double **Z, DATAINFO *pdinfo,
		      gretlopt opt)
{
    integer T, k;
    gretl_matrix *Q = NULL, *y = NULL;
    gretl_matrix *R = NULL, *g = NULL, *b = NULL;
    gretl_matrix *V = NULL;
    int rank, err = 0;

    T = pmod->nobs;               /* # of rows (observations) */
    k = pmod->list[0] - 1;        /* # of cols (variables) */

    Q = gretl_matrix_alloc(T, k);
    R = gretl_matrix_alloc(k, k);
    V = gretl_matrix_alloc(k, k);

    if (y == NULL) {
	y = gretl_matrix_alloc(T, 1);
    }

    if (Q == NULL || R == NULL || V == NULL || y == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    get_model_data(pmod, Z, Q, y);
    err = QR_decomp_plus(Q, R, &rank);

    /* handling of (near-)perfect collinearity */
    if (err == E_SINGULAR && !(opt & OPT_Z)) {
	drop_redundant_vars(pmod, R, rank, opt);
	k = pmod->list[0] - 1;
	gretl_matrix_reuse(Q, T, k);
	gretl_matrix_reuse(R, k, k);
	gretl_matrix_reuse(V, k, k);
	get_model_data(pmod, Z, Q, y);
	err = QR_decomp_plus(Q, R, NULL);
	if (!err) {
	    maybe_shift_ldepvar(pmod, Z, pdinfo);
	}
    }

    if (err) {
	goto qr_cleanup;
    }

    /* allocate temporary arrays */
    g = gretl_matrix_alloc(k, 1);
    b = gretl_matrix_alloc(k, 1);
    if (g == NULL || b == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    if (allocate_model_arrays(pmod, k, pdinfo->n)) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* make "g" into gamma-hat */    
    gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
			      y, GRETL_MOD_NONE, 
			      g, GRETL_MOD_NONE);

    /* OLS coefficients */
    gretl_matrix_multiply(R, g, b);
    pmod->coeff = gretl_matrix_steal_data(b);

    /* write vector of fitted values into y */
    gretl_matrix_multiply(Q, g, y);    

    /* get vector of residuals and SSR */
    get_resids_and_SSR(pmod, Z, y, pdinfo->n);

    /* standard error of regression */
    if (T - k > 0) {
	if (gretl_model_get_int(pmod, "no-df-corr")) {
	    pmod->sigma = sqrt(pmod->ess / T);
	} else {
	    pmod->sigma = sqrt(pmod->ess / (T - k));
	}
    } else {
	pmod->sigma = 0.0;
    }

    /* create (X'X)^{-1} */
    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      V, GRETL_MOD_NONE);

    /* VCV and standard errors */
    if (opt & OPT_R) { 
	gretl_model_set_int(pmod, "robust", 1);
	if ((opt & OPT_T) && !libset_get_bool(FORCE_HC)) {
	    qr_make_hac(pmod, Z, V);
	} else {
	    qr_make_hccme(pmod, Z, Q, V);
	}
    } else {
	qr_make_regular_vcv(pmod, V, opt);
    }

    /* get R^2, F */
    qr_compute_stats(pmod, Z[pmod->list[1]], T, opt);

    /* D-W stat and p-value */
    if ((opt & OPT_I) && pmod->missmask == NULL) {
	qr_dw_stats(pmod, Z, Q, y);
    }	

 qr_cleanup:

    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(y);

    gretl_matrix_free(g);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    pmod->errcode = err;

    return err;    
}

int qr_tsls_vcv (MODEL *pmod, const double **Z, const DATAINFO *pdinfo,
		 gretlopt opt)
{
    integer T, k;
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *V = NULL;
    int err = 0;

    T = pmod->nobs;               /* # of rows (observations) */
    k = pmod->list[0] - 1;        /* # of cols (variables) */

    Q = make_data_X(pmod, Z);
    R = gretl_matrix_alloc(k, k);
    V = gretl_matrix_alloc(k, k);

    if (Q == NULL || R == NULL || V == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    err = QR_decomp_plus(Q, R, NULL);
    if (err) {
	goto qr_cleanup;
    }

    /* create (X'X)^{-1} */
    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      V, GRETL_MOD_NONE);

    /* VCV and standard errors */
    if (opt & OPT_R) {
	if (dataset_is_panel(pdinfo)) {
	    err = qr_make_regular_vcv(pmod, V, OPT_X);
	    if (!err) {
		err = panel_tsls_robust_vcv(pmod, Z, pdinfo);
	    }
	} else if (dataset_is_time_series(pdinfo) && 
		   !libset_get_bool(FORCE_HC)) {
	    gretl_model_set_int(pmod, "robust", 1);
	    err = qr_make_hac(pmod, Z, V);
	} else {
	    gretl_model_set_int(pmod, "robust", 1);
	    err = qr_make_hccme(pmod, Z, Q, V);
	}
    } else {
	qr_make_regular_vcv(pmod, V, OPT_NONE);
    }
    
 qr_cleanup:

    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(V);

    pmod->errcode = err;

    return err;    
}

