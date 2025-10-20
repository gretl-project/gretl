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

/*  Apparatus for factorized variant of OLS */

#include "libgretl.h"
#include "matrix_extra.h"
#include "version.h"

static void add_ols_stuff (MODEL *pmod,
                           gretl_matrix *y,
                           gretl_matrix *b,
                           gretl_matrix *u,
                           gretl_matrix *V,
                           int T, int k,
                           int nfac)
{
    double cybar = 0.0;
    double ctss = 0.0;
    double d, s2;
    int t, s;

    pmod->ci = FOLS;
    pmod->dfn = k + nfac - 1;
    pmod->dfd = T - (k + nfac);

    pmod->uhat = malloc(pmod->full_n * sizeof *pmod->uhat);
    pmod->yhat = malloc(pmod->full_n * sizeof *pmod->yhat);
    pmod->coeff = b->val;
    b->val = NULL; /* donated */

    pmod->ess = 0.0;
    s = pmod->t1;
    for (t=0; t<T; t++) {
        pmod->ess += u->val[t] * u->val[t];
        pmod->uhat[s] = u->val[t];
        pmod->yhat[s] = y->val[t] - u->val[t];
        cybar += y->val[t];
        s++;
    }
    cybar /= T;
    for (t=0; t<T; t++) {
        d = y->val[t] - cybar;
        ctss += d * d;
    }

    /* rsq, fstt */
    pmod->rsq = (1 - (pmod->ess / pmod->tss));
    pmod->adjrsq = (1 - (pmod->ess / ctss));
    pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);

    s2 = pmod->ess / pmod->dfd;
    pmod->sigma = sqrt(s2);
    gretl_matrix_multiply_by_scalar(V, s2);
    gretl_model_write_vcv(pmod, V);
}

static void add_depvar_stats (MODEL *pmod, int yv, DATASET *dset)
{
    const double *y = dset->Z[yv];
    double d;
    int t;

    pmod->ybar = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        pmod->ybar += y[t];
    }
    pmod->ybar /= pmod->nobs;

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
        d = y[t] - pmod->ybar;
        pmod->tss += d * d;
    }
    pmod->sdy = sqrt(pmod->tss / (pmod->nobs - 1));
}

static void fols_basics (MODEL *pmod, int yv, DATASET *dset)
{
    pmod->t1 = dset->t1;
    pmod->t2 = dset->t2;
    pmod->nobs = sample_size(dset);
    pmod->full_n = dset->n;
    add_depvar_stats(pmod, yv, dset);
}

int fols_estimate (MODEL *pmod, int yv, int *xlist, int facv,
                   DATASET *dset, gretlopt opt, PRN *prn)
{
    gretl_matrix *fvals = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *u = NULL;
    gretl_matrix *V = NULL;
    const double *fac;
    double *means = NULL;
    int *njs = NULL;
    double x, fvi;
    int i, j, t, k, T;
    int nfvals;
    int err = 0;

    T = sample_size(dset);
    pmod->ncoeff = k = xlist[0];
    fols_basics(pmod, yv, dset);

    means = malloc((k+1) * sizeof *means);
    njs = malloc((k+1) * sizeof *njs);
    if (means == NULL || njs == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    fac = dset->Z[facv] + dset->t1;
    fvals = gretl_matrix_values(fac, T, OPT_S, &err);
    if (!err) {
        y = gretl_vector_from_series(dset->Z[yv], dset->t1, dset->t2);
        X = gretl_matrix_data_subset(xlist, dset, dset->t1, dset->t2,
                                     M_MISSING_ERROR, &err);
    }
    if (err) {
        goto bailout;
    }

    nfvals = gretl_vector_get_length(fvals);
    gretl_model_set_int(pmod, "nfvals", nfvals);

    for (i=0; i<nfvals; i++) {
        fvi = fvals->val[i];
        for (j=0; j<=k; j++) {
            means[j] = 0.0;
            njs[j] = 0;
        }
        for (t=0; t<T; t++) {
            if (fac[t] == fvi) {
                for (j=0; j<=k; j++) {
                    if (j == 0) {
                        means[j] += y->val[t];
                    } else {
                        means[j] += gretl_matrix_get(X, t, j-1);
                    }
                    njs[j] += 1;
                }
            }
        }
        for (j=0; j<=k; j++) {
            means[j] /= njs[j];
        }
        for (t=0; t<T; t++) {
            if (fac[t] == fvi) {
                for (j=0; j<=k; j++) {
                    if (j == 0) {
                        y->val[t] -= means[j];
                    } else {
                        x = gretl_matrix_get(X, t, j-1);
                        gretl_matrix_set(X, t, j-1, x - means[j]);
                    }
                }
            }
        }
    }

    b = gretl_matrix_alloc(k, 1);
    u = gretl_matrix_alloc(T, 1);
    V = gretl_matrix_alloc(k, k);
    err = gretl_matrix_ols(y, X, b, V, u, NULL);
    if (!err) {
        add_ols_stuff(pmod, y, b, u, V, T, k, nfvals);
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(fvals);
    gretl_matrix_free(u);
    gretl_matrix_free(V);
    free(means);
    free(njs);
    free(xlist); /* or set on pmod? */

    pmod->errcode = err;

    return err;
}
