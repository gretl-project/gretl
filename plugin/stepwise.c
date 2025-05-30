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

/* Forward stepwise plugin for gretl, using progressive QR decomposition */

#include "libgretl.h"
#include "version.h"
#include "matrix_extra.h"

/* Drop/cut a single column from matrix @m, provided that the given
   column index is valid. Return non-zero if the specification is not
   valid.
*/

static int matrix_drop_column (gretl_matrix *m, int drop)
{
    int i, j, k = drop * m->rows;
    double x;

    if (drop < 0 || drop >= m->cols) {
        fprintf(stderr, "matrix_drop_column: invalid index\n");
        return E_INVARG;
    }

    for (j=drop+1; j<m->cols; j++) {
        for (i=0; i<m->rows; i++) {
            x = gretl_matrix_get(m, i, j);
            m->val[k++] = x;
        }
    }

    /* don't leak column names */
    gretl_matrix_set_colnames(m, NULL);
    m->cols -= 1;

    return 0;
}

static void ssr2crit (int *best, double *xbest,
                      const gretl_matrix *ssr,
                      int T, int k, int crit)
{
    int n = ssr->cols;
    int j;

    *best = 0;
    *xbest = 1.0e200;

    if (crit == 1) {
        /* just using SSR */
        for (j=0; j<n; j++) {
            if (ssr->val[j] < *xbest) {
                *best = j;
                *xbest = ssr->val[j];
            }
        }
    } else {
        /* using AIC, BIC or HQC */
        double llj, icj;
        double l0 = -T/2.0;
        double l1 = log(2*M_PI);
        double C0 = 0;

        if (crit == 3) {
            C0 = log((double) T) * k; /* BIC */
        } else if (crit == 4) {
            C0 = 2 * log(log((double) T)) * k; /* HQC */
        }
        for (j=0; j<n; j++) {
            llj = l0 * (l1 + log(ssr->val[j]/T) + 1);
            if (crit == 2) {
                /* AIC */
                icj = 2.0 * (k - llj);
            } else {
                /* BIC or HQC */
                icj = C0 - 2 * llj;
            }
            if (icj < *xbest) {
                *xbest = icj;
                *best = j;
            }
        }
    }
}

static int qr_update (gretl_matrix **pQ,
                      gretl_matrix **pR,
                      const gretl_matrix *Z,
                      const gretl_matrix *e,
                      int crit, int *best,
                      double *xbest)
{
    gretl_matrix_block *Blk;
    gretl_matrix *Q = *pQ;
    gretl_matrix *R = *pR;
    gretl_matrix *B;
    gretl_matrix *E;
    gretl_matrix *den;
    gretl_matrix *std;
    gretl_matrix *stdres;
    gretl_matrix *ssr;
    gretl_matrix *tmp;
    double xij, ee;
    int n = Q->rows;
    int k = Q->cols;
    int zc = Z->cols;
    int i, j;
    int err = 0;

    Blk = gretl_matrix_block_new(&B, k, zc,
                                 &E, n, zc,
                                 &den, 1, zc,
                                 &std, 1, zc,
                                 &stdres, n, zc,
                                 &ssr, 1, zc,
                                 NULL);
    if (Blk == NULL) {
        return E_ALLOC;
    }

    /* B = Q'Z */
    gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
                              Z, GRETL_MOD_NONE,
                              B, GRETL_MOD_NONE);
    /* E = Z - Q*B */
    gretl_matrix_copy_values(E, Z);
    gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                              B, GRETL_MOD_NONE,
                              E, GRETL_MOD_DECREMENT);

    for (j=0; j<zc; j++) {
        /* den = sumc(E.^2) */
        den->val[j] = 0;
        for (i=0; i<n; i++) {
            xij = gretl_matrix_get(E, i, j);
            den->val[j] += xij * xij;
        }
        /* std = sqrt(den) */
        std->val[j] = sqrt(den->val[j]);
        /* stdres = E ./ std */
        for (i=0; i<n; i++) {
            xij = gretl_matrix_get(E, i, j);
            gretl_matrix_set(stdres, i, j, xij / std->val[j]);
        }
    }

    /* ee = e'e */
    ee = 0;
    for (i=0; i<n; i++) {
        ee += e->val[i] * e->val[i];
    }

    /* num2 = (e'Z).^2 ; gain = num2./den ; ssr = e'e - gain */
    gretl_matrix_multiply_mod(e, GRETL_MOD_TRANSPOSE,
                              Z, GRETL_MOD_NONE,
                              ssr, GRETL_MOD_NONE);
    for (j=0; j<zc; j++) {
        xij = ssr->val[j];
        ssr->val[j] = ee - xij * xij / den->val[j];
    }

    ssr2crit(best, xbest, ssr, n, k + 1, crit);

    /* update Q -> Q ~ -stdres[,best] */
    tmp = gretl_matrix_alloc(n, k + 1);
    memcpy(tmp->val, Q->val, n * k * sizeof *tmp->val);
    for (i=0; i<stdres->rows; i++) {
        xij = gretl_matrix_get(stdres, i, *best);
        gretl_matrix_set(tmp, i, k, -xij);
    }
    gretl_matrix_free(Q);
    *pQ = tmp;

    /* update R = (R|0)  ~ (B[,best] | -std[best]) */
    tmp = gretl_matrix_alloc(k + 1, k + 1);
    for (j=0; j<k; j++) {
        for (i=0; i<k; i++) {
            xij = gretl_matrix_get(R, i, j);
            gretl_matrix_set(tmp, i, j, xij);
        }
        gretl_matrix_set(tmp, k, j, 0.0);
    }
    for (i=0; i<k; i++) {
        xij = gretl_matrix_get(B, i, *best);
        gretl_matrix_set(tmp, i, k, xij);
    }
    xij = std->val[*best];
    gretl_matrix_set(tmp, k, k, -xij);
    gretl_matrix_free(R);
    *pR = tmp;

    gretl_matrix_block_destroy(Blk);

    return err;
}

static const char *cstrs[] = {
    "SSR", "AIC", "BIC", "HQC"
};

static const char *crit_string (int crit)
{
    return cstrs[crit-1];
}

int *forward_stepwise (MODEL *pmod,
                       const int *zlist,
                       DATASET *dset,
                       int crit,
                       double alpha,
                       int verbose,
                       PRN *prn,
                       int *err)
{
    const char *cstr;
    gretl_matrix *e;
    gretl_matrix *Q;
    gretl_matrix *R;
    gretl_matrix *mZ;
    gretl_matrix *my;
    gretl_matrix *tmp;
    gretl_matrix *ssr;
    int conv = 0;
    int added = 0;
    int best = 0;
    int bpos;
    double cur, prev;
    int *xlist = NULL;
    int *aux = NULL;
    int *ret = NULL;
    int yvar;
    int t1, t2;
    int T, k, i, nz;

    T = pmod->nobs;
    k = pmod->ncoeff;
    t1 = pmod->t1;
    t2 = pmod->t2;
    e = gretl_matrix_alloc(T, 1);
    for (i=0; i<pmod->nobs; i++) {
        e->val[i] = pmod->uhat[i];
    }
    ssr = gretl_matrix_from_scalar(pmod->ess);
    ssr2crit(&best, &prev, ssr, T, k, crit);
    gretl_matrix_free(ssr);

    yvar = pmod->list[1];
    xlist = gretl_list_new(pmod->list[0] - 1);
    for (i=2; i<=pmod->list[0]; i++) {
        xlist[i-1] = pmod->list[i];
    }

    Q  = gretl_matrix_data_subset(xlist, dset, t1, t2, M_MISSING_ERROR, err);
    mZ = gretl_matrix_data_subset(zlist, dset, t1, t2, M_MISSING_ERROR, err);
    my = gretl_vector_from_series(dset->Z[yvar], t1, t2);
    R = gretl_zero_matrix_new(k, k);

    gretl_matrix_QR_decomp(Q, R);
    aux = gretl_list_copy(zlist);
    nz = zlist[0];
    cstr = crit_string(crit);

    if (verbose) {
        printf("\nBaseline:            %s = %#g\n", cstr, prev);
    }

    while (!conv && added < nz) {
        *err = qr_update(&Q, &R, mZ, e, crit, &best, &cur);
        if (*err) {
            printf("got error %d from qr_update\n", *err);
            break;
        }

        /* e = my - Q*(Q'my) */
        tmp = gretl_matrix_alloc(Q->cols, 1);
        gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
                                  my, GRETL_MOD_NONE,
                                  tmp, GRETL_MOD_NONE);
        gretl_matrix_copy_values(e, my);
        gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                                  tmp, GRETL_MOD_NONE,
                                  e, GRETL_MOD_DECREMENT);
        gretl_matrix_free(tmp);

        if (crit == 1) {
            /* SSR */
            double parm[1] = {1.0};
            double Xcrit = gretl_get_cdf_inverse(D_CHISQ, parm, 1.0 - alpha);
            double W = T * (prev/cur - 1.0);

            conv = W < Xcrit;
        } else {
            /* AIC, etc. */
            conv = cur > prev;
        }

        if (verbose) {
            if (conv && added < nz) {
                printf("[%-19s %s = %#g]\n", dset->varname[aux[best+1]],
                       cstr, cur);
            } else {
                printf("Add %-16s %s = %#g\n", dset->varname[aux[best+1]],
                       cstr, cur);
            }
        }

        if (!conv) {
            bpos = best + 1;
            matrix_drop_column(mZ, best);
            gretl_list_append_term(&ret, aux[bpos]);
            gretl_list_delete_at_pos(aux, bpos);
            added++;
            prev = cur;
        }
    }

    free(aux);
    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(mZ);
    gretl_matrix_free(my);
    gretl_matrix_free(e);

    return ret;
}

/* In case of forward stepwise regression, check @zlist (the list of
   candidate regressors to be added) for the possibility that it may
   contain one or more of the baseline regressors.
*/

static int stepwise_check_zlist (MODEL *pmod,
                                 const int *zlist,
                                 DATASET *dset)
{
    int i;

    for (i=1; i<=zlist[0]; i++) {
        if (in_gretl_list(pmod->list, zlist[i])) {
            return E_ADDDUP;
        }
    }

    return 0;
}

/* Process the --silent, --quiet and --auto options to the "add"
   command. In the --auto case we should have either the standard
   abbreviation for one of the Information Criteria, or an alpha
   value for use with the plain SSR criterion.
*/

static int process_stepwise_options (gretlopt opt,
                                     int *crit,
                                     double *alpha,
                                     int *verbosity)
{
    const char *s = get_optval_string(ADD, OPT_A);
    int i, err = 0;

    if (opt & OPT_I) {
        *verbosity = 0;
    } else if (opt & OPT_Q) {
        *verbosity = 1;
    } else {
        *verbosity = 2;
    }

    for (i=1; i<4; i++) {
        /* AIC, BIC or HQC */
        if (!strcmp(s, cstrs[i])) {
            *crit = i + 1;
        }
    }

    if (*crit == 0) {
        /* SSR */
        *alpha = gretl_double_from_string(s, &err);
        if (!err && (*alpha < 0.001 || *alpha > 0.99)) {
            err = E_INVARG;
        } else {
            *crit = 1;
        }
    }

    return err;
}

/* compose the final OLS list based on the original model plus the
   'best' list of added regressors
*/

static int *compose_list (MODEL *pmod, int *best)
{
    int n = pmod->list[0] + best[0];
    int i, j = 1;
    int *list;

    list = gretl_list_new(n);
    for (i=1; i<=pmod->list[0]; i++) {
        list[j++] = pmod->list[i];
    }
    for (i=1; i<=best[0]; i++) {
        list[j++] = best[i];
    }    
   
    return list;
}

/* implement "add --auto=<param>" using forward stepwise procedure */

MODEL stepwise_add (MODEL *pmod,
                    const int *zlist,
                    DATASET *dset,
                    gretlopt opt,
                    PRN *prn)
{
    MODEL model;
    int *best = NULL;
    int crit = 0;
    int verbosity = 2;
    int ols_done = 0;
    double alpha = 0;
    int err;

    err = stepwise_check_zlist(pmod, zlist, dset);

    if (!err) {
        err = process_stepwise_options(opt, &crit, &alpha, &verbosity);
    }

    if (!err) {
        best = forward_stepwise(pmod, zlist, dset, crit, alpha,
                                verbosity, prn, &err);
    }

    if (!err) {
        gretlopt ols_opt = OPT_NONE;
        int *list = compose_list(pmod, best);

        if (verbosity == 0) {
            ols_opt = OPT_Q;
        }
        model = lsq(list, dset, OLS, ols_opt);
        free(list);
        ols_done = 1;
    }

    free(best);

    if (!ols_done) {
        gretl_model_init(&model, NULL);
        model.errcode = err;
    }

    return model;
}
