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

/* qr_wspace: workspace matrices whose row dimension will remain
   unchanged but whose column dimension @zc will shrink on each call
   to qr_update().
*/

typedef struct qr_wspace_ {
    gretl_matrix_block *B; /* holder */
    gretl_matrix *E;       /* n x zc */
    gretl_matrix *den;     /* 1 x zc */
    gretl_matrix *std;     /* 1 x zc */
    gretl_matrix *stdres;  /* n x zc */
    gretl_matrix *ssr;     /* 1 x zc */
} qr_wspace;

static int qr_wspace_alloc (qr_wspace *mm, int n, int zc)
{
    mm->B = gretl_matrix_block_new(&mm->E,   n, zc,
                                   &mm->den, 1, zc,
                                   &mm->std, 1, zc,
                                   &mm->stdres, n, zc,
                                   &mm->ssr, 1, zc,
                                   NULL);
    if (mm->B == NULL) {
        return E_ALLOC;
    } else {
        return 0;
    }
}

static void qr_wspace_shrink (qr_wspace *mm)
{
    mm->E->cols   -= 1;
    mm->den->cols -= 1;
    mm->std->cols -= 1;
    mm->stdres->cols -= 1;
    mm->ssr->cols -= 1;
}

static void qr_wspace_free (qr_wspace *mm)
{
    gretl_matrix_block_destroy(mm->B);
}

/* Drop/cut a single column from matrix @m */

static void matrix_drop_column (gretl_matrix *m, int drop)
{
    int i, j, k = drop * m->rows;
    double x;

    for (j=drop+1; j<m->cols; j++) {
        for (i=0; i<m->rows; i++) {
            x = gretl_matrix_get(m, i, j);
            m->val[k++] = x;
        }
    }

    /* don't leak column names */
    if (m->info != NULL) {
        gretl_matrix_destroy_info(m);
    }

    m->cols -= 1;
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

static int qr_update (gretl_matrix *Q,
                      gretl_matrix *R,
                      const gretl_matrix *Z,
                      const gretl_matrix *e,
                      qr_wspace *mm,
                      int crit, int *best,
                      double *xbest)
{
    gretl_matrix *B;
    double *dest;
    double *src;
    double xij, ee;
    size_t sz;
    int n = Q->rows;
    int k = Q->cols;
    int zc = Z->cols;
    int i, j, p;
    int err = 0;

    B = gretl_matrix_alloc(k, zc);
    if (B == NULL) {
        return E_ALLOC;
    }

    /* B = Q'Z */
    gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
                              Z, GRETL_MOD_NONE,
                              B, GRETL_MOD_NONE);
    /* E = Z - Q*B */
    gretl_matrix_copy_values(mm->E, Z);
    gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                              B, GRETL_MOD_NONE,
                              mm->E, GRETL_MOD_DECREMENT);

    for (j=0; j<zc; j++) {
        /* den = sumc(E.^2) */
        mm->den->val[j] = 0;
        for (i=0; i<n; i++) {
            xij = gretl_matrix_get(mm->E, i, j);
            mm->den->val[j] += xij * xij;
        }
        /* std = sqrt(den) */
        mm->std->val[j] = sqrt(mm->den->val[j]);
        /* stdres = E ./ std */
        for (i=0; i<n; i++) {
            xij = gretl_matrix_get(mm->E, i, j);
            gretl_matrix_set(mm->stdres, i, j, xij / mm->std->val[j]);
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
                              mm->ssr, GRETL_MOD_NONE);
    for (j=0; j<zc; j++) {
        xij = mm->ssr->val[j];
        mm->ssr->val[j] = ee - xij * xij / mm->den->val[j];
    }

    ssr2crit(&p, xbest, mm->ssr, n, k + 1, crit);
    *best = p;

    /* update Q -> Q ~ -stdres[,best] */
    gretl_matrix_realloc(Q, n, k + 1);
    for (i=0; i<mm->stdres->rows; i++) {
        xij = gretl_matrix_get(mm->stdres, i, p);
        gretl_matrix_set(Q, i, k, -xij);
    }

    /* update R = (R|0)  ~ (B[,best] | -std[best]) */
    gretl_matrix_realloc(R, k+1, k+1);
    sz = k * sizeof *dest;
    for (j=k-1; j>0; j--) {
        src = R->val + k * j;
        dest = src + j;
        memmove(dest, src, sz);
        dest[k] = 0.0;
    }
    /* fill the last column */
    src = B->val + k * p;
    dest = R->val + k*(k-1);
    memcpy(dest, src, sz);
    xij = mm->std->val[p];
    gretl_matrix_set(R, k, k, -xij);

    gretl_matrix_free(B);

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
                       int addlen,
                       int namelen,
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
    qr_wspace mm = {0};
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

    if (!*err && (my == NULL || R == NULL)) {
        *err = E_ALLOC;
    }
    if (*err) {
        goto bailout;
    }

    gretl_matrix_QR_decomp(Q, R);
    aux = gretl_list_copy(zlist);
    nz = zlist[0];
    cstr = crit_string(crit);
    qr_wspace_alloc(&mm, Q->rows, nz);

    if (verbose) {
        pprintf(prn, "\n%-*s %s = %#g\n", addlen + namelen + 1,
                _("Baseline"), cstr, prev);
    }

    while (!conv && added < nz) {
        *err = qr_update(Q, R, mZ, e, &mm, crit, &best, &cur);
        if (*err) {
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
                pprintf(prn, "[%-*s %s = %#g]\n", namelen + addlen,
                        dset->varname[aux[best+1]], cstr, cur);
            } else {
                pprintf(prn, "%s %-*s %s = %#g\n", _("Add"), namelen,
                        dset->varname[aux[best+1]], cstr, cur);
            }
        }

        if (!conv) {
            bpos = best + 1;
            matrix_drop_column(mZ, best);
            qr_wspace_shrink(&mm);
            gretl_list_append_term(&ret, aux[bpos]);
            gretl_list_delete_at_pos(aux, bpos);
            added++;
            prev = cur;
        }
    }

 bailout:

    free(aux);
    free(xlist);
    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(mZ);
    gretl_matrix_free(my);
    gretl_matrix_free(e);
    qr_wspace_free(&mm);

    return ret;
}

/* In case of forward stepwise regression, check @zlist (list of
   candidate regressors) for the possibility that it may contain one
   or more of the baseline regressors. While we're at it, get the
   maximum length of the names of the @zlist members.
*/

static int stepwise_check_zlist (MODEL *pmod,
                                 const int *zlist,
                                 DATASET *dset,
                                 int *len)
{
    int i, n;

    for (i=1; i<=zlist[0]; i++) {
        if (in_gretl_list(pmod->list, zlist[i])) {
            return E_ADDDUP;
        }
        n = strlen(dset->varname[zlist[i]]);
        if (n > *len) {
            *len = n;
        }
    }

    return 0;
}

/* Process the --auto option to the "add" command. We should have
   either the standard abbreviation for one of the Information
   Criteria, or an alpha value for use with the SSR criterion.
*/

static int process_stepwise_option (gretlopt opt,
                                    int *crit,
                                    double *alpha)
{
    const char *s = get_optval_string(ADD, OPT_A);
    int i, err = 0;

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

/* Compose the final OLS list based on the original model plus the
   input @zlist and the @best list of added regressors.
*/

static int *compose_list (MODEL *pmod, const int *best,
                          const int *zlist)
{
    int n = pmod->list[0] + best[0];
    int i, j = 1;
    int *list;

    list = gretl_list_new(n);
    for (i=1; i<=pmod->list[0]; i++) {
        list[j++] = pmod->list[i];
    }
    for (i=1; i<=zlist[0]; i++) {
        if (in_gretl_list(best, zlist[i])) {
            list[j++] = zlist[i];
        }
    }

    return list;
}

static void do_overall_test (MODEL *orig, MODEL *revised)
{
    double W = orig->nobs * (revised->rsq - orig->rsq);
    int dk = revised->ncoeff - orig->ncoeff;

    if (dk > 0) {
        /* should we print this? */
        double parm[] = {dk};
        double pval = gretl_get_pvalue(D_CHISQ, parm, W);
        record_test_result(W, pval);
    }
}

/* Implement "add --auto=..." using forward stepwise procedure */

MODEL stepwise_add (MODEL *pmod,
                    const int *zlist,
                    DATASET *dset,
                    gretlopt opt,
                    PRN *prn)
{
    MODEL model;
    int *best = NULL;
    int crit = 0;
    int ols_done = 0;
    int namelen = 0;
    double alpha = 0;
    int err;

    err = stepwise_check_zlist(pmod, zlist, dset, &namelen);

    if (!err) {
        err = process_stepwise_option(opt, &crit, &alpha);
    }

    if (!err) {
        int verbose = (opt & OPT_Q)? 0 : 1;
        int addlen = verbose ? g_utf8_strlen(_("Add"), -1) : 0;

        best = forward_stepwise(pmod, zlist, dset, crit, alpha,
                                verbose, addlen, namelen + 2,
                                prn, &err);
    }

    if (!err) {
        gretlopt ols_opt = OPT_NONE;
        int *list = compose_list(pmod, best, zlist);

        if (opt & OPT_I) {
            /* --silent */
            ols_opt |= OPT_Q;
        }
        if (opt & OPT_O) {
            /* --vcv */
            ols_opt |= OPT_O;
        }
        model = lsq(list, dset, OLS, ols_opt);
        free(list);
        ols_done = 1;
        if (!model.errcode) {
            do_overall_test(pmod, &model);
        }
    }

    free(best);

    if (!ols_done) {
        gretl_model_init(&model, NULL);
        model.errcode = err;
    }

    return model;
}
