/*
 *   Copyright (c) 2003-2004 by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "qr_estimate.h"
#include "gretl_matrix.h"
#include "libset.h"
#include "estim_private.h"

#include "f2c.h"
#include "clapack_double.h"

#define QR_RCOND_MIN 1e-15 /* experiment with this? */
#define ESSZERO      1e-22 /* SSR less than this counts as zero */

/* General note: in fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are contiguous.
*/

static double qr_get_tss (MODEL *pmod, const double *y, int *ifc)
{
    int pwe = gretl_model_get_int(pmod, "pwe");
    double ymean = 0.0;
    double x, tss = 0.0;
    int t;

    if (*ifc == 0) {
	*ifc = check_for_effective_const(pmod, y);
    }

    if (pmod->rho != 0.0) {
	double ry, d;

	if (*ifc) {
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
	}

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    ry = y[t];
	    if (t == pmod->t1 && pwe) {
		ry *= sqrt(1.0 - pmod->rho * pmod->rho);
	    } else {
		ry -= pmod->rho * y[t-1];
	    }
	    d = ry - ymean;
	    tss += d * d;
	}
    } else {
	if (*ifc) {
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		if (!na(pmod->yhat[t])) {
		    ymean += y[t];
		}
	    }
	    ymean /= pmod->nobs;
	}

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->yhat[t])) {
		x = y[t] - ymean;
		tss += x * x;
	    }
	}
    } 

    return tss;
}

static void qr_compute_stats (MODEL *pmod, const double *y, int n,
			      gretlopt opt)
{
    int ifc = pmod->ifc;

    pmod->tss = qr_get_tss(pmod, y, &ifc);

    if (pmod->dfd > 0) {
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

static int get_vcv_index (MODEL *pmod, int i, int j, int n)
{
    int k, vi, vj;

    if (pmod->ifc) {
	vi = (i + 1) % n;
	vj = (j + 1) % n;
    } else {
	vi = i;
	vj = j;
    }

    k = ijton(vi, vj, n);

    return k;
}

static int qr_make_vcv (MODEL *pmod, gretl_matrix *v, int robust)
{
    int k = pmod->ncoeff;
    int nterms = k * (k + 1) / 2;
    double x;
    int i, j, idx;

    pmod->vcv = malloc(nterms * sizeof *pmod->vcv);
    if (pmod->vcv == NULL) return 1;  

    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    idx = get_vcv_index(pmod, i, j, k);
	    x = gretl_matrix_get(v, i, j);
	    if (!robust) {
		x *= pmod->sigma * pmod->sigma;
	    }
	    pmod->vcv[idx] = x;
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
    int i, j, t;
    int start = (pmod->ifc)? 3 : 2;

    /* copy independent vars into matrix X */
    j = 0;
    for (i=start; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!model_missing(pmod, t)) {
		if (pmod->nwt) {
		    X->val[j++] = sqrt(Z[pmod->nwt][t]) * Z[pmod->list[i]][t];
		} else {
		    X->val[j++] = Z[pmod->list[i]][t];
		}
	    }
	}
    }

    /* insert constant */
    if (pmod->ifc) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!model_missing(pmod, t)) {
		if (pmod->nwt) {
		    X->val[j++] = sqrt(Z[pmod->nwt][t]);
		} else {
		    X->val[j++] = 1.0;
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
	xi = X->val[mdx(X, t, i)];
	for (j=0; j<n; j++) {
	    xj = X->val[mdx(X, t - lag, j)];
	    wt->val[mdx(wt, i, j)] = xi * xj;
	}
    }
}

/* Calculate the Newey-West HAC covariance matrix.  Algorithm and
   (basically) notation taken from Davidson and MacKinnon (DM), 
   Econometric Theory and Methods, chapter 9.
*/

static int qr_make_hac (MODEL *pmod, const double **Z, gretl_matrix *xpxinv)
{
    gretl_matrix *vcv = NULL, *wtj = NULL, *gammaj = NULL;
    gretl_matrix *X;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int p, i, j, t;
    double weight, uu;
    double *uhat;
    int err = 0;

    X = make_data_X(pmod, Z);
    if (X == NULL) return 1;

    /* pmod->uhat is a full-length series: we must take an offset
       into it, equal to the offset of the data on which the model
       is actually estimated.
    */
    uhat = pmod->uhat + pmod->t1;

    /* get the user's preferred maximum lag setting */
    p = get_hac_lag(T);
    gretl_model_set_int(pmod, "hac_lag", p);

    vcv = gretl_matrix_alloc(k, k);
    wtj = gretl_matrix_alloc(k, k);
    gammaj = gretl_matrix_alloc(k, k);

    if (vcv == NULL || wtj == NULL || gammaj == NULL) {
	err = 1;
	goto bailout;
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
	    weight = 1.0 - (double) j / (p + 1.0);
	    /* multiply by Newey-West weight */
	    gretl_matrix_multiply_by_scalar(gammaj, weight);
	}

	/* DM equation (9.38), p. 364 */
	gretl_matrix_add_to(vcv, gammaj);
    }

    gretl_matrix_multiply_mod(xpxinv, GRETL_MOD_TRANSPOSE,
			      vcv, GRETL_MOD_NONE,
			      wtj, GRETL_MOD_NONE);
    gretl_matrix_multiply(wtj, xpxinv, vcv);

    /* vcv now holds HAC */
    for (i=0; i<k; i++) {
	double x = gretl_matrix_get(vcv, i, i);

	j = (pmod->ifc)? (i + 1) % k : i;
	pmod->sderr[j] = sqrt(x);
    }

    /* Transcribe vcv into triangular representation */
    err = qr_make_vcv(pmod, vcv, 1);

 bailout:

    gretl_matrix_free(wtj);
    gretl_matrix_free(gammaj);
    gretl_matrix_free(vcv);
    gretl_matrix_free(X);

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
    const double *d;
    double x;
    int i, j;

    for (i=0; i<R->rows; i++) {
	d = D->val;
	for (j=0; j<R->cols; j++) {
	    x = X->val[mdx(X, j, i)];
	    R->val[mdx(R, i, j)] = x * (*d++);
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
    gretl_matrix *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL;
    int T = pmod->nobs; 
    int k = pmod->list[0] - 1;
    int hc_version;
    int i, j, t;
    int err = 0;

    X = make_data_X(pmod, Z);
    if (X == NULL) return 1;

    diag = gretl_column_vector_from_array(pmod->uhat + pmod->t1, T,
					  GRETL_MOD_SQUARE);
    if (diag == NULL) {
	err = 1;
	goto bailout;
    }  

    tmp1 = gretl_matrix_alloc(k, T);
    tmp2 = gretl_matrix_alloc(k, k);
    tmp3 = gretl_matrix_alloc(k, k);
    if (tmp1 == NULL || tmp2 == NULL || tmp3 == NULL) {
	err = 1;
	goto bailout;
    }  

    hc_version = get_hc_version();
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
		q = Q->val[mdx(Q, t, i)];
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
    gretl_matrix_multiply(xpxinv, tmp2, tmp3); 
    gretl_matrix_multiply(tmp3, xpxinv, tmp2);

    /* tmp2 now holds HCCM */
    for (i=0; i<k; i++) {
	double x = tmp2->val[mdx(tmp2, i, i)];

	j = (pmod->ifc)? (i + 1) % k : i;
	pmod->sderr[j] = sqrt(x);
    }

    err = qr_make_vcv(pmod, tmp2, 1);

    bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(diag);
    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);
    gretl_matrix_free(tmp3);

    return err;
}

static int qr_make_regular_vcv (MODEL *pmod, gretl_matrix *v)
{
    int i, j, k = pmod->ncoeff;

    for (i=0; i<k; i++) {
	double x = v->val[mdx(v, i, i)];

	j = (pmod->ifc)? (i + 1) % k : i;
	pmod->sderr[j] = pmod->sigma * sqrt(x);
    }

    return qr_make_vcv(pmod, v, 0);
}

static void get_model_data (MODEL *pmod, const double **Z, 
			    gretl_matrix *Q, gretl_matrix *y)
{
    int i, j, t, start;
    double x;
    int dwt = gretl_model_get_int(pmod, "wt_dummy");
    int qdiff = (pmod->rho != 0.0);
    int pwe = gretl_model_get_int(pmod, "pwe");
    double pw1 = 0.0;

    if (pwe) {
	pw1 = sqrt(1.0 - pmod->rho * pmod->rho);
    } 

    start = (pmod->ifc)? 3 : 2;

    if (dwt) {
	dwt = pmod->nwt;
    }

    /* copy independent vars into matrix Q */
    j = 0;
    for (i=start; i<=pmod->list[0]; i++) {
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

    /* insert constant last (numerical issues) */
    if (pmod->ifc) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (model_missing(pmod, t)) {
		continue;
	    }

	    /* in some special cases the constant is pre-transformed */
	    x = Z[0][t];

	    if (dwt) {
		if (Z[dwt][t] == 0.0) continue;
	    } else if (pmod->nwt) {
		x = sqrt(Z[pmod->nwt][t]);
	    } else if (qdiff) {
		if (pwe && t == pmod->t1) {
		    x = pw1;
		} else {
		    x = 1.0 - pmod->rho;
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

static void save_coefficients (MODEL *pmod, gretl_matrix *b,
			       int n)
{
    pmod->coeff = gretl_matrix_steal_data(b);

    if (pmod->ifc) {
	int i;
	double tmp = pmod->coeff[n-1];

	for (i=n-1; i>0; i--) {
	    pmod->coeff[i] = pmod->coeff[i-1];
	}
	pmod->coeff[0] = tmp;
    }
}

static int 
allocate_model_arrays (MODEL *pmod, int k, int T)
{
    pmod->sderr = malloc(k * sizeof *pmod->sderr);
    pmod->yhat = malloc(T * sizeof *pmod->yhat);
    pmod->uhat = malloc(T * sizeof *pmod->uhat);

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

static int *make_droplist (MODEL *pmod, gretl_matrix *R, 
			   int rank, int *err)
{
    int *dlist = NULL;
    double d;
    int nd = R->rows - rank;
    int i, j;

    dlist = gretl_list_new(nd);
    if (dlist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    dlist[0] = 0;

    j = 2;
    for (i=0; i<R->rows; i++) {
	d = gretl_matrix_get(R, i, i);
	if (fabs(d) < R_DIAG_MIN) {
	    dlist[0] += 1;
	    dlist[dlist[0]] = pmod->list[j];
	    gretl_list_delete_at_pos(pmod->list, j--);
	}
	j++;
    }

    pmod->ncoeff -= dlist[0];
    pmod->dfd = pmod->nobs - pmod->ncoeff;
    pmod->dfn = pmod->ncoeff - pmod->ifc;

    return dlist;
}

int gretl_qr_regress (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		      gretlopt opts)
{
    integer T, k;
    gretl_matrix *Q = NULL, *y = NULL;
    gretl_matrix *R = NULL, *g = NULL, *b = NULL;
    gretl_matrix *xpxinv = NULL;
    int *droplist = NULL;
    int rank, err = 0;

    T = pmod->nobs;               /* # of rows (observations) */
    k = pmod->list[0] - 1;        /* # of cols (variables) */

    Q = gretl_matrix_alloc(T, k);
    R = gretl_matrix_alloc(k, k);
    xpxinv = gretl_matrix_alloc(k, k);

    if (y == NULL) {
	y = gretl_matrix_alloc(T, 1);
    }

    if (Q == NULL || R == NULL || xpxinv == NULL || y == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

 trim_var:    

    get_model_data(pmod, (const double **) *pZ, Q, y);
    err = QR_decomp_plus(Q, R, &rank);

    /* alternate ways of handling (near-)perfect collinearity */

#if 1
    if (err == E_SINGULAR && !(opts & OPT_Z)) {
	err = 0;
	droplist = make_droplist(pmod, R, rank, &err);
	if (!err) {
	    k = pmod->list[0] - 1;
	    gretl_matrix_reuse(Q, T, k);
	    gretl_matrix_reuse(R, k, k);
	    gretl_matrix_reuse(xpxinv, k, k);
	    goto trim_var;
	}
    }
#else
    if (err == E_SINGULAR && !(opts & OPT_Z) &&
	redundant_var(pmod, pZ, pdinfo, &droplist)) {
	err = 0;
	k = pmod->list[0] - 1;
	gretl_matrix_reuse(Q, T, k);
	gretl_matrix_reuse(R, k, k);
	gretl_matrix_reuse(xpxinv, k, k);
	goto trim_var;
    } 
#endif

    if (err) {
	goto qr_cleanup;
    }

    if (droplist != NULL) {
	maybe_shift_ldepvar(pmod, (const double **) *pZ, pdinfo);
	gretl_model_set_list_as_data(pmod, "droplist", droplist);
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
    save_coefficients(pmod, b, k);

    /* write vector of fitted values into y */
    gretl_matrix_multiply(Q, g, y);    

    /* get vector of residuals and SSR */
    get_resids_and_SSR(pmod, (const double **) *pZ, y, pdinfo->n);

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
			      xpxinv, GRETL_MOD_NONE);

    /* VCV and standard errors */
    if (opts & OPT_R) { 
	gretl_model_set_int(pmod, "robust", 1);
	if ((opts & OPT_T) && !get_force_hc()) {
	    qr_make_hac(pmod, (const double **) *pZ, xpxinv);
	} else {
	    qr_make_hccme(pmod, (const double **) *pZ, Q, xpxinv);
	}
    } else {
	qr_make_regular_vcv(pmod, xpxinv);
    }

    /* get R^2, F */
    qr_compute_stats(pmod, (*pZ)[pmod->list[1]], T, opts);

 qr_cleanup:

    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(y);

    gretl_matrix_free(g);
    gretl_matrix_free(b);
    gretl_matrix_free(xpxinv);

    pmod->errcode = err;

    return err;    
}

int qr_tsls_vcv (MODEL *pmod, const double **Z, gretlopt opts)
{
    integer T, k;
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *xpxinv = NULL;
    int err = 0;

    T = pmod->nobs;               /* # of rows (observations) */
    k = pmod->list[0] - 1;        /* # of cols (variables) */

    Q = make_data_X(pmod, Z);
    R = gretl_matrix_alloc(k, k);
    xpxinv = gretl_matrix_alloc(k, k);

    if (Q == NULL || R == NULL || xpxinv == NULL) {
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
			      xpxinv, GRETL_MOD_NONE);

    /* VCV and standard errors */
    if (opts & OPT_R) { 
	gretl_model_set_int(pmod, "robust", 1);
	if (opts & OPT_T) {
	    qr_make_hac(pmod, Z, xpxinv);
	} else {
	    qr_make_hccme(pmod, Z, Q, xpxinv);
	}
    } else {
	qr_make_regular_vcv(pmod, xpxinv);
    }

 qr_cleanup:

    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(xpxinv);

    pmod->errcode = err;

    return err;    
}

