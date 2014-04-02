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

#define ESSZERO 1e-22 /* SSR less than this counts as zero */

enum {
    VCV_SIMPLE,
    VCV_ROBUST,
    VCV_XPX
};

static int qr_make_cluster_vcv (MODEL *pmod, int ci,
				const DATASET *dset,
				gretl_matrix *XX,
				gretlopt opt);

/* General note: in fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are contiguous.
*/

static double qr_get_tss (MODEL *pmod, const DATASET *dset,
			  int *ifc, int *yconst)
{
    int yno = pmod->list[1];
    int pwe = (pmod->opt & OPT_P);
    double yt, y0 = 0.0, ymean = 0.0;
    double x, tss = 0.0;
    double ctss = 0.0;
    int t;

    if (*ifc == 0) {
	*ifc = check_for_effective_const(pmod, dset);
    }

    *yconst = 1;

    if (pmod->rho != 0.0) {
	double ry, d;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    ry = dset->Z[yno][t];
	    if (t == pmod->t1 && pwe) {
		ry *= sqrt(1.0 - pmod->rho * pmod->rho);
	    } else {
		ry -= pmod->rho * dset->Z[yno][t-1];
	    }
	    ymean += ry;
	}
	ymean /= pmod->nobs;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    ry = dset->Z[yno][t];
	    if (t == pmod->t1 && pwe) {
		ry *= sqrt(1.0 - pmod->rho * pmod->rho);
	    } else {
		ry -= pmod->rho * dset->Z[yno][t-1];
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
		ymean += dset->Z[yno][t];
	    }
	}
	ymean /= pmod->nobs;

	y0 = dset->Z[yno][pmod->t1];

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->yhat[t])) {
		if (dset->Z[yno][t] != y0) {
		    *yconst = 0;
		}
		x = dset->Z[yno][t] - ymean;
		if (*ifc) {
		    tss += x * x;
		} else {
		    yt = dset->Z[yno][t];
		    tss += yt * yt;
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

static void qr_compute_stats (MODEL *pmod, const DATASET *dset,
			      int n, gretlopt opt)
{
    int yconst, ifc = pmod->ifc;
    int yno = pmod->list[1];

    pmod->tss = qr_get_tss(pmod, dset, &ifc, &yconst);
    pmod->chisq = NADBL;

    if (yconst && pmod->dfd > 0) {
	double y0 = dset->Z[yno][pmod->t1];
    
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
	    pmod->fstt = wald_omit_F(NULL, pmod);
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

static void get_resids_and_SSR (MODEL *pmod, const DATASET *dset,
				gretl_matrix *yhat, int fulln)
{
    int dwt = gretl_model_get_int(pmod, "wt_dummy");
    int qdiff = (pmod->rho != 0.0);
    int pwe = (pmod->opt & OPT_P);
    int yvar = pmod->list[1];
    double *u = pmod->uhat;
    double y;
    int t, i = 0;

    if (dwt) {
	dwt = pmod->nwt;
    }

    pmod->ess = 0.0;

    if (qdiff) {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		pmod->yhat[t] = u[t] = NADBL;
	    } else {
		y = dset->Z[yvar][t];
		if (t == pmod->t1 && pwe) {
		    y *= sqrt(1.0 - pmod->rho * pmod->rho);
		} else {
		    y -= pmod->rho * dset->Z[yvar][t-1];
		}
		pmod->yhat[t] = yhat->val[i];
		u[t] = y - yhat->val[i];
		pmod->ess += u[t] * u[t];
		i++;
	    }
	}
    } else if (pmod->nwt) {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2 || model_missing(pmod, t)) {
		pmod->yhat[t] = u[t] = NADBL;
	    } else {
		y = dset->Z[yvar][t];
		if (dwt && dset->Z[dwt][t] == 0.0) {
		    pmod->yhat[t] = NADBL;
		} else {
		    if (!dwt) {
			y *= sqrt(dset->Z[pmod->nwt][t]);
		    }
		    pmod->yhat[t] = yhat->val[i];
		    u[t] = y - yhat->val[i];
		    pmod->ess += u[t] * u[t];
		    i++;
		}
	    }
	}
    } else {
	for (t=0; t<fulln; t++) {
	    if (t < pmod->t1 || t > pmod->t2 || model_missing(pmod, t)) {
		pmod->yhat[t] = u[t] = NADBL;
	    } else {
		pmod->yhat[t] = yhat->val[i];
		u[t] = dset->Z[yvar][t] - yhat->val[i];
		pmod->ess += u[t] * u[t];
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
get_data_X (gretl_matrix *X, const MODEL *pmod, const DATASET *dset)
{
    int wt = pmod->nwt;
    int i, t, vi, s = 0;

    /* copy independent vars into matrix X */

    if (!wt && pmod->missmask == NULL) {
	/* use faster method */
	int T = pmod->t2 - pmod->t1 + 1;
	size_t sz = T * sizeof(double);

	for (i=2; i<=pmod->list[0]; i++) {
	    vi = pmod->list[i];
	    memcpy(X->val + s, dset->Z[vi] + pmod->t1, sz);
	    s += T;
	}
    } else {
	for (i=2; i<=pmod->list[0]; i++) {
	    vi = pmod->list[i];
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!model_missing(pmod, t)) {
		    if (wt) {
			X->val[s++] = sqrt(dset->Z[wt][t]) * dset->Z[vi][t];
		    } else {
			X->val[s++] = dset->Z[vi][t];
		    }
		}
	    }
	}
    }
}

static gretl_matrix *make_data_X (const MODEL *pmod, 
				  const DATASET *dset)
{
    gretl_matrix *X;

    X = gretl_matrix_alloc(pmod->nobs, pmod->ncoeff);
    if (X != NULL) {
	get_data_X(X, pmod, dset);
    }

    return X;
}

/* Calculate W(t)-transpose * W(t-lag) */

static void wtw (gretl_matrix *wt, const gretl_matrix *X, 
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

static double *prewhiten_uhat (const double *u, int T, double *pa)
{
    double a, num = 0.0, den = 0.0;
    double *uw;
    int sgn, t;

    uw = malloc(T * sizeof *uw);
    if (uw == NULL) {
	return NULL;
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

    return uw;
}

gretl_matrix *HAC_XOX (const gretl_matrix *uhat, const gretl_matrix *X,
		       VCVInfo *vi, int *err)
{
    gretl_matrix *XOX = NULL, *Wtj = NULL, *Gj = NULL;
    int prewhiten = libset_get_bool(PREWHITEN);
    int kern = libset_get_int(HAC_KERNEL);
    int T = X->rows;
    int k = X->cols;
    int p, j, t;
    double wj, uu;
    double a = 0, bt = 0;
    double *u = NULL;

    if (prewhiten) {
	u = prewhiten_uhat(uhat->val, T, &a);
	if (u == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}
    } else {
	u = uhat->val;
    }

    XOX = gretl_zero_matrix_new(k, k);
    Wtj = gretl_matrix_alloc(k, k);
    Gj = gretl_matrix_alloc(k, k);

    if (XOX == NULL || Wtj == NULL || Gj == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* determine the bandwidth setting */

    if (data_based_hac_bandwidth()) {
	gretl_matrix umat;

	gretl_matrix_init(&umat);
	umat.rows = T;
	umat.cols = 1;
	umat.val = u;
	*err = newey_west_bandwidth(&umat, kern, &p, &bt);
	if (*err) {
	    goto bailout;
	}
    } else if (kern == KERNEL_QS) {
	bt = libset_get_double(QS_BANDWIDTH);
	p = T - 1;
    } else {
	p = get_hac_lag(T);
    }

    for (j=0; j<=p; j++) {
	/* cumulate running sum of Gamma-hat terms */
	gretl_matrix_zero(Gj);
	for (t=j; t<T; t++) {
	    /* W(t)-transpose * W(t-j) */
	    wtw(Wtj, X, k, t, j);
	    uu = u[t] * u[t-j];
	    gretl_matrix_multiply_by_scalar(Wtj, uu);
	    /* DM equation (9.36), p. 363 */
	    gretl_matrix_add_to(Gj, Wtj);
	}

	if (j > 0) {
	    /* Gamma(j) = Gamma(j) + Gamma(j)-transpose */
	    gretl_matrix_add_self_transpose(Gj);
	    if (kern == KERNEL_QS) {
		wj = qs_hac_weight(bt, j);
	    } else {
		wj = hac_weight(kern, p, j);
	    }
	    gretl_matrix_multiply_by_scalar(Gj, wj);
	}

	/* DM equation (9.38), p. 364 */
	gretl_matrix_add_to(XOX, Gj);
    }

    if (prewhiten) {
	/* re-color */
	gretl_matrix_divide_by_scalar(XOX, (1-a) * (1-a));
    }

    vi->vmaj = VCV_HAC;
    vi->vmin = kern;
    vi->flags = prewhiten;

    if (kern == KERNEL_QS) {
	vi->order = 0;
	vi->bw = bt;
    } else {
	vi->order = p;
	vi->bw = NADBL;
    }

 bailout:

    gretl_matrix_free(Wtj);
    gretl_matrix_free(Gj);

    if (u != uhat->val) {
	free(u);
    }

    if (*err && XOX != NULL) {
	gretl_matrix_free(XOX);
	XOX = NULL;
    }

    return XOX;
}

/* Calculate HAC covariance matrix.  Algorithm and (basically)
   notation taken from Davidson and MacKinnon (DM), Econometric Theory
   and Methods, chapter 9.
*/

static int qr_make_hac (MODEL *pmod, const DATASET *dset, 
			gretl_matrix *XTXi)
{
    gretl_matrix *X, *XOX, *V = NULL;
    gretl_matrix umat;
    VCVInfo vi;
    int T = pmod->nobs;
    int err = 0;

    X = make_data_X(pmod, dset);
    if (X == NULL) {
	return E_ALLOC;
    }

    /* pmod->uhat is a full-length series: we must take an offset
       into it, equal to the offset of the data on which the model
       is actually estimated.
    */
    gretl_matrix_init(&umat);
    umat.rows = T;
    umat.cols = 1;
    umat.val = pmod->uhat + pmod->t1;

    XOX = HAC_XOX(&umat, X, &vi, &err);

    if (!err) {
	V = gretl_matrix_alloc(XOX->rows, XOX->rows);
	if (V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	gretl_matrix_qform(XTXi, GRETL_MOD_TRANSPOSE, XOX,
			   V, GRETL_MOD_NONE);
	/* Transcribe vcv into triangular representation */
	err = qr_make_vcv(pmod, V, VCV_ROBUST);
    }

    if (!err) {
	gretl_model_set_full_vcv_info(pmod, vi.vmaj, vi.vmin,
				      vi.order, vi.flags,
				      vi.bw);
    }	

    gretl_matrix_free(X);
    gretl_matrix_free(XOX);
    gretl_matrix_free(V);

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

static int qr_make_hccme (MODEL *pmod, const DATASET *dset,
			  gretl_matrix *Q, gretl_matrix *XTXi)
{
    gretl_matrix *X;
    gretl_matrix *diag = NULL;
    gretl_matrix *tmp1 = NULL, *tmp2 = NULL, *vcv = NULL;
    int T = pmod->nobs; 
    int k = pmod->list[0] - 1;
    int hc_version;
    int i, t;
    int err = 0;

    X = make_data_X(pmod, dset);
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
    gretl_model_set_vcv_info(pmod, VCV_HC, hc_version);

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
    gretl_matrix_qform(XTXi, GRETL_MOD_NONE, tmp2,
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

static int qr_dw_stats (MODEL *pmod, const DATASET *dset,
			gretl_matrix *X, gretl_matrix *u)
{
    double DW, pv;
    int t, s, err = 0;

    get_data_X(X, pmod, dset);

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

static void get_model_data (MODEL *pmod, const DATASET *dset,
			    gretl_matrix *Q, gretl_matrix *y)
{
    int dwt = gretl_model_get_int(pmod, "wt_dummy");
    int qdiff = (pmod->rho != 0.0);
    int pwe = (pmod->opt & OPT_P);
    double x, pw1 = 0.0;
    int i, s, t;

    if (pmod->missmask == NULL && !pwe && !qdiff && !dwt && !pmod->nwt) {
	/* simple case: no missing values and no data transformation
	   called for, so use faster procedure 
	*/
	int T = pmod->t2 - pmod->t1 + 1;
	size_t sz = T * sizeof(double);

	s = 0;
	for (i=2; i<=pmod->list[0]; i++) {
	    memcpy(Q->val + s, dset->Z[pmod->list[i]] + pmod->t1, sz);
	    s += T;
	}
	if (y != NULL) {
	    memcpy(y->val, dset->Z[pmod->list[1]] + pmod->t1, sz);
	}
	return;
    }

    if (pwe) {
	pw1 = sqrt(1.0 - pmod->rho * pmod->rho);
    }

    if (dwt) {
	dwt = pmod->nwt;
    }

    /* copy independent vars into matrix Q */
    s = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	int vi = pmod->list[i];

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }
	    x = dset->Z[vi][t];
	    if (dwt) {
		if (dset->Z[dwt][t] == 0.0) continue;
	    } else if (pmod->nwt) {
		x *= sqrt(dset->Z[pmod->nwt][t]);
	    } else if (qdiff) {
		if (pwe && t == pmod->t1) {
		    x *= pw1;
		} else {
		    x -= pmod->rho * dset->Z[vi][t-1];
		}
	    }
	    Q->val[s++] = x;
	}
    }

    if (y != NULL) {
	/* copy dependent variable into y vector */
	int vy = pmod->list[1];

	s = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }		
	    x = dset->Z[vy][t];
	    if (dwt) {
		if (dset->Z[dwt][t] == 0.0) continue;
	    } else if (pmod->nwt) {
		x *= sqrt(dset->Z[pmod->nwt][t]);
	    } else if (qdiff) {
		if (pwe && t == pmod->t1) {
		    x *= pw1;
		} else {
		    x -= pmod->rho * dset->Z[vy][t-1];
		}
	    }
	    y->val[s++] = x;
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

    if (pmod->sderr == NULL || pmod->yhat == NULL || 
	pmod->uhat == NULL) {
	return 1;
    }

    return 0;
}

#define RCOND_WARN 1.0e-07

/* perform QR decomposition plus some additional tasks */

static int QR_decomp_plus (gretl_matrix *Q, gretl_matrix *R, int *rank,
			   int *warn)
{
    integer k = gretl_matrix_rows(R);
    double rcond = 0;
    int r, err;

    if (warn != NULL) {
	*warn = 0;
    }

    /* basic decomposition */
    err = gretl_matrix_QR_decomp(Q, R);
    if (err) {
	return err;
    }

    /* check rank of QR */
    r = gretl_check_QR_rank(R, &err, &rcond);
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
	} else if (rcond < RCOND_WARN && warn != NULL) {
	    *warn = 1;
	}
    }

    if (rank != NULL) {
	*rank = r;
    }

    return err;
}

#define REDEBUG 0

static void
drop_redundant_vars (MODEL *pmod, DATASET *dset, gretl_matrix *R, 
		     int rank, gretlopt opt)
{
    int *droplist = NULL;
    int i, vi, pos, nd;
    double d;

#if REDEBUG
    printlist(pmod->list, "pmod->list, into drop_redundant_vars");
    fprintf(stderr, "rank = %d\n", rank);
#endif

    pos = 2;
    nd = 0;
    for (i=0; i<R->rows; i++) {
	d = gretl_matrix_get(R, i, i);
	if (fabs(d) < R_DIAG_MIN) {
	    vi = pmod->list[pos];
	    gretl_list_append_term(&droplist, vi);
	    fprintf(stderr, "dropping redundant variable %d (%s)\n",
		    vi, dset->varname[vi]);
	    gretl_list_delete_at_pos(pmod->list, pos--);
	    nd++;
	}
	pos++;
    }

    pmod->ncoeff -= nd;
    pmod->dfd = pmod->nobs - pmod->ncoeff;
    pmod->dfn = pmod->ncoeff - pmod->ifc;

    if (droplist != NULL) {
	gretl_model_set_list_as_data(pmod, "droplist", droplist);
    }
}

int gretl_qr_regress (MODEL *pmod, DATASET *dset, gretlopt opt)
{
    integer T, k;
    gretl_matrix *Q = NULL, *y = NULL;
    gretl_matrix *R = NULL, *g = NULL, *b = NULL;
    gretl_matrix *V = NULL;
    int rank, warn = 0, err = 0;

    T = pmod->nobs;        /* # of rows (observations) */
    k = pmod->list[0] - 1; /* # of cols (variables) */

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

    get_model_data(pmod, dset, Q, y);
    err = QR_decomp_plus(Q, R, &rank, &warn);

    /* handling of near-perfect collinearity */
    if (err == E_SINGULAR && !(opt & OPT_Z)) {
	drop_redundant_vars(pmod, dset, R, rank, opt);
	k = pmod->list[0] - 1;
	gretl_matrix_reuse(Q, T, k);
	gretl_matrix_reuse(R, k, k);
	gretl_matrix_reuse(V, k, k);
	get_model_data(pmod, dset, Q, y);
	err = QR_decomp_plus(Q, R, NULL, &warn);
	if (!err) {
	    maybe_shift_ldepvar(pmod, dset);
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

    if (allocate_model_arrays(pmod, k, dset->n)) {
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
    get_resids_and_SSR(pmod, dset, y, dset->n);

    /* standard error of regression */
    if (T - k > 0) {
	if (pmod->opt & OPT_N) {
	    /* no-df-corr */
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
	pmod->opt |= OPT_R;
	if (opt & OPT_C) {
	    err = qr_make_cluster_vcv(pmod, OLS, dset, V, opt);
	} else if ((opt & OPT_T) && !libset_get_bool(FORCE_HC)) {
	    err = qr_make_hac(pmod, dset, V);
	} else {
	    err = qr_make_hccme(pmod, dset, Q, V);
	}
    } else {
	err = qr_make_regular_vcv(pmod, V, opt);
    }

    if (!err) {
	/* get R^2, F-stat */
	qr_compute_stats(pmod, dset, T, opt);

	/* D-W stat and p-value */
	if ((opt & OPT_I) && pmod->missmask == NULL) {
	    qr_dw_stats(pmod, dset, Q, y);
	}

	/* near-singularity? */
	if (warn) {
	    gretl_model_set_int(pmod, "near-singular", 1);
	}
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

int qr_tsls_vcv (MODEL *pmod, const DATASET *dset, gretlopt opt)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *V = NULL;
    int k, err = 0;

    k = pmod->list[0] - 1;

    Q = make_data_X(pmod, dset);
    R = gretl_matrix_alloc(k, k);
    V = gretl_matrix_alloc(k, k);

    if (Q == NULL || R == NULL || V == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    err = QR_decomp_plus(Q, R, NULL, NULL);
    if (err) {
	goto qr_cleanup;
    }

    /* create (X'X)^{-1} */
    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      V, GRETL_MOD_NONE);

    /* VCV and standard errors */
    if (opt & OPT_R) {
	if (opt & OPT_C) {
	    pmod->opt |= OPT_R;
	    err = qr_make_cluster_vcv(pmod, IVREG, dset, V, opt);
	} else if (dataset_is_panel(dset)) {
	    err = qr_make_regular_vcv(pmod, V, OPT_X);
	    if (!err) {
		err = panel_tsls_robust_vcv(pmod, dset);
	    }
	} else if (dataset_is_time_series(dset) && 
		   !libset_get_bool(FORCE_HC)) {
	    pmod->opt |= OPT_R;
	    err = qr_make_hac(pmod, dset, V);
	} else {
	    pmod->opt |= OPT_R;
	    err = qr_make_hccme(pmod, dset, Q, V);
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

static int cval_count (MODEL *pmod, double cvi, const double *cZ)
{
    int t, cc = 0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(pmod->uhat[t]) && cZ[t] == cvi) {
	    cc++;
	}
    }

    return cc;
}

static int cval_count_max (MODEL *pmod, const gretl_matrix *cvals, 
			   const double *cZ)
{
    int n = gretl_vector_get_length(cvals);
    int i, cc, cmax = 0;

    for (i=0; i<n; i++) {
	cc = cval_count(pmod, cvals->val[i], cZ);
	if (cc > cmax) {
	    cmax = cc;
	}
    }

    return cmax;
}

#define CDEBUG 0

static gretl_matrix *cluster_vcv_calc (MODEL *pmod,
				       int cvar,
				       gretl_matrix *cvals, 
				       gretl_matrix *XX,
				       const DATASET *dset,
				       int *err)

{
    gretl_matrix *V = NULL;
    gretl_matrix *W = NULL;
    gretl_matrix *XXW = NULL;
    gretl_vector *ei = NULL;
    gretl_matrix *Xi = NULL;
    gretl_vector *eXi = NULL;
    const double *cZ;
    int n_c, M, N, k = pmod->ncoeff;
    int total_obs = 0;
    int i, j, v, t;

    cZ = dset->Z[cvar];    
    N = cval_count_max(pmod, cvals, cZ);
#if CDEBUG
    fprintf(stderr, "max cval count = %d\n", N);
#endif

    V   = gretl_matrix_alloc(k, k);
    W   = gretl_zero_matrix_new(k, k);
    XXW = gretl_zero_matrix_new(k, k);
    ei  = gretl_column_vector_alloc(N);
    Xi  = gretl_matrix_alloc(N, k);
    eXi = gretl_vector_alloc(k);

    if (V == NULL || W == NULL || XXW == NULL || 
	ei == NULL || Xi == NULL || eXi == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    M = gretl_vector_get_length(cvals);
    n_c = 0;

    for (i=0; i<M; i++) {
	double cvi = cvals->val[i];
	int Ni = cval_count(pmod, cvi, cZ);
	int s = 0;

	if (Ni == 0) {
	    continue;
	}

#if CDEBUG
	fprintf(stderr, "i=%d, cvi=%g, Ni=%d\n", i, cvi, Ni);
#endif
	ei = gretl_matrix_reuse(ei, Ni, -1);
	Xi = gretl_matrix_reuse(Xi, Ni, -1);

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t]) && cZ[t] == cvi) {
		gretl_vector_set(ei, s, pmod->uhat[t]);
		for (j=0; j<k; j++) {
		    v = pmod->list[j+2];
		    gretl_matrix_set(Xi, s, j, dset->Z[v][t]);
		}
		s++;
	    }
	    if (s == Ni) {
		/* we've filled this matrix */
		break;
	    }
	}

	gretl_matrix_multiply_mod(ei, GRETL_MOD_TRANSPOSE,
				  Xi, GRETL_MOD_NONE,
				  eXi, GRETL_MOD_NONE);
	gretl_matrix_multiply_mod(eXi, GRETL_MOD_TRANSPOSE,
				  eXi, GRETL_MOD_NONE,
				  W, GRETL_MOD_CUMULATE);
#if CDEBUG > 1
	gretl_matrix_print(ei, "e(i)");
	gretl_matrix_print(Xi, "X(i)");
	gretl_matrix_print(W, "W");
#endif
	n_c++;
	total_obs += s;
    }

    if (n_c < 2) {
	gretl_errmsg_set("Invalid clustering variable");
	*err = E_DATA;
	goto bailout;
    } else if (total_obs < pmod->nobs) {
	*err = E_MISSDATA;
	goto bailout;
    }

    /* form V(W) = (X'X)^{-1} W (X'X)^{-1} */
    gretl_matrix_multiply(XX, W, XXW);
    gretl_matrix_multiply(XXW, XX, V);
    gretl_matrix_xtr_symmetric(V);

#if CDEBUG
    gretl_matrix_print(XX, "X'X^{-1}");
    gretl_matrix_print(W, "W");
    gretl_matrix_print(V, "V");
#endif

    if (!(pmod->opt & OPT_N)) {
	/* apply df adjustment a la Stata */
	double dfadj;

	N = pmod->nobs;
	dfadj = (M/(M-1.0)) * (N-1.0)/(N-k);
	gretl_matrix_multiply_by_scalar(V, dfadj);
#if CDEBUG > 1
	gretl_matrix_print(V, "V(adjusted)");
#endif

    }

 bailout:

    gretl_matrix_free(W);
    gretl_matrix_free(XXW);
    gretl_matrix_free(ei);
    gretl_matrix_free(Xi);
    gretl_matrix_free(eXi);

    if (*err) {
	gretl_matrix_free(V);
	V = NULL;
    }

    return V;
}

/* Get the sorted values of the clustering series, @cvar, checking
   for missing values as we go. This is a little more complicated
   if there are interior missing values for the regressand or
   regressors: in that case we need to construct a temporary
   array to hold the relevant values of the clustering series.
*/

static gretl_matrix *cluster_var_values (const double *cvar,
					 MODEL *pmod,
					 int *err)
{
    gretl_matrix *cvals = NULL;
    int t;

    if (pmod->missmask != NULL) {
	double *ctmp = malloc(pmod->nobs * sizeof *ctmp);

	if (ctmp == NULL) {
	    *err = E_ALLOC;
	} else {
	    int i = 0;

	    for (t=pmod->t1; t<=pmod->t2 && !*err; t++) {
		if (pmod->missmask[t] != '1') {
		    if (na(cvar[t])) {
			*err = E_MISSDATA;
		    } else {
			ctmp[i++] = cvar[t];
		    }
		}
	    }

	    if (!*err) {
		cvals = gretl_matrix_values(ctmp, pmod->nobs, OPT_S, err);
	    }
	    free(ctmp);
	}
    } else {
	/* no interior missing values for y or X */
	for (t=pmod->t1; t<=pmod->t2 && !*err; t++) {
	    if (na(cvar[t])) {
		*err = E_MISSDATA;
	    }
	}
	if (!*err) {
	    cvals = gretl_matrix_values(cvar + pmod->t1, pmod->nobs, 
					OPT_S, err);
	}
    }

    return cvals;
}

/**
 * qr_make_cluster_vcv:
 * @pmod: pointer to model.
 * @ci: command index (right now, OLS or IVREG).
 * @dset: pointer to dataset.
 * @XX: X'X matrix.
 * 
 * Compute and set on @pmod a variance matrix that is "clustered"
 * by the value of a selected variable via the --cluster=foo 
 * command-line option.
 *
 * Returns: 0 on success, non-zero code on error.
 */

static int qr_make_cluster_vcv (MODEL *pmod, int ci,
				const DATASET *dset,
				gretl_matrix *XX,
				gretlopt opt)
{
    gretl_matrix *cvals = NULL;
    gretl_matrix *V = NULL;
    const char *cname;
    int cvar, n_c = 0;
    int err = 0;

    if (pmod->ci != OLS && pmod->ci != IVREG) {
	/* relax this? */
	return E_NOTIMP;
    }

    cname = get_optval_string(ci, OPT_C); 
    if (cname == NULL) {
	return E_PARSE;
    }

    cvar = current_series_index(dset, cname);
    if (cvar < 1 || cvar >= dset->v) {
	err = E_UNKVAR;
    }

    if (!err) {
	cvals = cluster_var_values(dset->Z[cvar], pmod, &err);
	if (!err) {
	    n_c = gretl_vector_get_length(cvals);
	    if (n_c < 2) {
		err = E_DATA;
	    }
	}
    }

#if CDEBUG
    fprintf(stderr, "qr_make_cluster_vcv: err = %d\n", err);
    fprintf(stderr, "cluster var = %s (%d)\n", cname, cvar);
    gretl_matrix_print(cvals, "cvals");
#endif

    if (!err) {
	V = cluster_vcv_calc(pmod, cvar, cvals, XX, dset, &err);
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, V);
    }

    if (!err) {
	gretl_model_set_vcv_info(pmod, VCV_CLUSTER, cvar);
	gretl_model_set_int(pmod, "n_clusters", n_c);
    }

    gretl_matrix_free(V);
    gretl_matrix_free(cvals);

    return err;
}
