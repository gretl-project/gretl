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

enum {SSR = 1, AIC, BIC, HQC};

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

static double get_info_crit (double ssr, int T, int k, int crit)
{
    /* using AIC, BIC or HQC */
    double l1 = log(2*M_PI);
    double ll, C0 = 0;

    if (crit == BIC) {
        C0 = k * log((double) T);
    } else if (crit == HQC) {
        C0 = 2 * k * log(log((double) T));
    }
    ll = -0.5 * T * (l1 + log(ssr/T) + 1);
    if (crit == AIC) {
        return 2.0 * (k - ll);
    } else {
        return C0 - 2 * ll;
    }
}

static double best_ssr2crit (int *best,
                             const gretl_matrix *ssr,
                             int T, int k, int crit)
{
    int j, n = ssr->cols;
    double ssr_min = 1.0e200;

    *best = 0;

    /* check SSR first */
    for (j=0; j<n; j++) {
        if (ssr->val[j] < ssr_min) {
            *best = j;
            ssr_min = ssr->val[j];
        }
    }

    if (crit == SSR) {
        return ssr_min;
    } else {
        return get_info_crit(ssr_min, T, k, crit);
    }
}

static double ssr2crit (double ssr, int T, int k, int crit)
{
    if (crit == SSR) {
        return ssr;
    } else {
        return get_info_crit(ssr, T, k, crit);
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

    *xbest = best_ssr2crit(&p, mm->ssr, n, k + 1, crit);
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

static const char *crit_string (int ci, int crit)
{
    if (ci == OMIT && crit == 1) {
        return "P-value";
    } else {
        return cstrs[crit-1];
    }
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
    prev = ssr2crit(pmod->ess, T, k, crit);

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
    cstr = crit_string(ADD, crit);
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

        if (crit == SSR) {
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

/* Returns the index of the column of the X matrix which seems to be the
   prime candidate for dropping next: the one associated with the
   smallest absolute t-ratio. We're ignoring the constant, if present,
   and if @zlist is non-NULL we're also ignoring any regressors that
   are not specified therein.
*/

static int tval_min_pos (const double *b,
                         const double *se,
                         const int *xlist,
                         const int *zlist,
                         int ifc, int T, int k,
                         double *pval)
{
    double tval, tmin = 1.0e200;
    int i, vi, ret = 0;

    for (i=ifc; i<k; i++) {
        if (zlist != NULL) {
            vi = xlist[i+1];
            if (!in_gretl_list(zlist, vi)) {
                /* vi is not a candidate for dropping */
                continue;
            }
        }
        tval = fabs(b[i]) / se[i];
        if (tval < tmin) {
            tmin = tval;
            ret = i;
        }
    }

    if (pval != NULL) {
        *pval = student_pvalue_2(T - k, tmin);
    }

    return ret;
}

int *backward_stepwise (MODEL *pmod,
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
    gretl_matrix *y;
    gretl_matrix *X;
    gretl_matrix *b;
    gretl_matrix *V;
    int *xlist = NULL;
    double *se = NULL;
    double cur, prev;
    double ssr, pval;
    double s2;
    int ifc = pmod->ifc;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int t1 = pmod->t1;
    int t2 = pmod->t2;
    int yvar = pmod->list[1];
    int conv = 0;
    int trycol, delvar;
    int nz, dropped;
    int i;

    prev = ssr2crit(pmod->ess, T, k, crit);
    xlist = gretl_model_get_x_list(pmod);
    trycol = tval_min_pos(pmod->coeff, pmod->sderr,
                          xlist, zlist, ifc, T, k, &pval);
    if (crit == SSR && pval < alpha) {
        printf("Nothing to be dropped\n");
        return NULL;
    }

    k--;
    y = gretl_vector_from_series(dset->Z[yvar], t1, t2);
    X = gretl_matrix_data_subset(xlist, dset, t1, t2, M_MISSING_ERROR, err);
    b = gretl_matrix_alloc(k, 1);
    V = gretl_matrix_alloc(k, k);
    se = malloc(k * sizeof *se);

    dropped = 0;
    nz = zlist != NULL ? zlist[0] : xlist[0] - ifc;
    cstr = crit_string(OMIT, crit);

    while (!conv && dropped < nz) {
        delvar = xlist[trycol+1];
        matrix_drop_column(X, trycol);
        k = X->cols;
        gretl_matrix_reuse(b, k, 1);
        gretl_matrix_reuse(V, k, k);
        *err = gretl_matrix_ols(y, X, b, V, NULL, &s2);
        if (*err) {
            break;
        }
        ssr = (T - k) * s2;
        if (crit > SSR) {
            /* AIC, etc. */
            cur = ssr2crit(ssr, T, k, crit);
            conv = cur > prev;
        } else {
            cur = pval;
        }
        if (verbose) {
            if (conv && dropped < nz) {
                pprintf(prn, "[%-*s %s = %#g]\n", namelen + addlen,
                        dset->varname[delvar], cstr, cur);
            } else {
                pprintf(prn, "%s %-*s %s = %#g\n", _("Drop"), namelen,
                        dset->varname[delvar], cstr, cur);
            }
        }
        if (!conv) {
            for (i=ifc; i<V->cols; i++) {
                se[i] = sqrt(gretl_matrix_get(V, i, i));
            }
            gretl_list_delete_at_pos(xlist, trycol+1);
            trycol = tval_min_pos(b->val, se, xlist, zlist, ifc, T, k, &pval);
            if (crit == SSR && pval < alpha) {
                conv = 1;
            }
            dropped++;
            prev = cur;
        }
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);
    free(se);

    if (*err) {
        free(xlist);
        xlist = NULL;
    }

    return xlist;
}

/* In case of forward stepwise regression, check @zlist (list of
   candidate regressors) for the possibility that it may contain one
   or more of the baseline regressors. While we're at it, get the
   maximum length of the names of the @zlist members.
*/

static int forward_stepwise_check_zlist (MODEL *pmod,
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

/* In case of backward stepwise regression, check @zlist (list of
   candidate for dropping) for the possibility that it may contain one
   or more terms that are not among the original. While we're at it, get
   the maximum length of the names of the relevant regressors members.
*/

static int backward_stepwise_check_zlist (MODEL *pmod,
                                          const int *zlist,
                                          DATASET *dset,
                                          int *len)
{
    int i, n;

    if (zlist != NULL) {
        for (i=1; i<=zlist[0]; i++) {
            if (zlist[i] == pmod->list[1] ||
                !in_gretl_list(pmod->list, zlist[i])) {
                return E_INVARG;
            }
            n = strlen(dset->varname[zlist[i]]);
            if (n > *len) {
                *len = n;
            }
        }
    } else {
        for (i=2; i<=pmod->list[0]; i++) {
            n = strlen(dset->varname[pmod->list[i]]);
            if (n > *len) {
                *len = n;
            }
        }
    }

    return 0;
}

/* Process the --auto option to the "add" command. We should have
   either the standard abbreviation for one of the Information
   Criteria, or an alpha value for use with the SSR criterion.
*/

static int process_stepwise_option (int ci,
                                    gretlopt opt,
                                    int *crit,
                                    double *alpha)
{
    const char *s = get_optval_string(ci, OPT_A);
    int i, err = 0;

    if (ci == OMIT && s == NULL) {
        /* omit: the --auto parameter is optional */
        *alpha = 0.10;
        return 0;
    }

    for (i=1; i<4; i++) {
        /* AIC, BIC or HQC */
        if (!strcmp(s, cstrs[i])) {
            *crit = i + 1;
        }
    }

    if (*crit == 0) {
        /* not yet determined */
        *alpha = gretl_double_from_string(s, &err);
        if (!err && (*alpha < 0.001 || *alpha > 0.99)) {
            err = E_INVARG;
        } else {
            *crit = SSR;
        }
    }

    return err;
}

/* Compose the final OLS list based on the original model plus the
   input @zlist and the @best list of added regressors.
*/

static int *compose_forward_list (const MODEL *pmod,
                                  const int *best,
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

static int *compose_backward_list (const MODEL *pmod,
                                   const int *xlist)
{
    int *list = gretl_list_new(1 + xlist[0]);
    int i;

    list[1] = pmod->list[1];
    for (i=1; i<=xlist[0]; i++) {
        list[i+1] = xlist[i];
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

    err = forward_stepwise_check_zlist(pmod, zlist, dset, &namelen);

    if (!err) {
        err = process_stepwise_option(ADD, opt, &crit, &alpha);
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
        int *list = compose_forward_list(pmod, best, zlist);

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

/* Implement "omit --auto=..." using backward stepwise procedure */

MODEL stepwise_omit (MODEL *pmod,
                     const int *zlist,
                     DATASET *dset,
                     gretlopt opt,
                     PRN *prn)
{
    MODEL model;
    int *xlist = NULL;
    int crit = 0;
    int ols_done = 0;
    int namelen = 0;
    double alpha = 0;
    int err;

    err = backward_stepwise_check_zlist(pmod, zlist, dset, &namelen);

    if (!err) {
        err = process_stepwise_option(OMIT, opt, &crit, &alpha);
    }

    if (!err) {
        int verbose = (opt & OPT_Q)? 0 : 1;
        int addlen = verbose ? g_utf8_strlen(_("Drop"), -1) : 0;

        xlist = backward_stepwise(pmod, zlist, dset, crit, alpha,
                                  verbose, addlen, namelen + 2,
                                  prn, &err);
    }

    if (!err) {
        gretlopt ols_opt = OPT_NONE;
        int *list = compose_backward_list(pmod, xlist);

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

    free(xlist);

    if (!ols_done) {
        gretl_model_init(&model, NULL);
        model.errcode = err;
    }

    return model;
}
