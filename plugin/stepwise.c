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

#define C_SSR C_MAX

/* fwd_wspace: workspace matrices whose row dimension will remain
   unchanged but whose column dimension @zc will shrink on each call
   to qr_augment().
*/

typedef struct fwd_wspace_ {
    gretl_matrix_block *B; /* holder */
    gretl_matrix *E;       /* n x zc */
    gretl_matrix *den;     /* 1 x zc */
    gretl_matrix *std;     /* 1 x zc */
    gretl_matrix *stdres;  /* n x zc */
    gretl_matrix *ssr;     /* 1 x zc */
} fwd_wspace;

/* bwd_wspace: workspace matrices whose @k dimension will shrink on
   each call to qr_reduce().
*/

typedef struct bwd_wspace_ {
    gretl_matrix_block *B; /* holder */
    gretl_matrix *iR;      /* R-inverse: k x k */
    gretl_matrix *g;       /* k x 1 */
    gretl_matrix *b;       /* k x 1 */
    gretl_matrix *V;       /* k x k */
    gretl_matrix *e;       /* T x 1 */
} bwd_wspace;

static int fwd_wspace_alloc (fwd_wspace *mm, int n, int zc)
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

static void fwd_wspace_shrink (fwd_wspace *mm)
{
    mm->E->cols   -= 1;
    mm->den->cols -= 1;
    mm->std->cols -= 1;
    mm->stdres->cols -= 1;
    mm->ssr->cols -= 1;
}

static void fwd_wspace_free (fwd_wspace *mm)
{
    gretl_matrix_block_destroy(mm->B);
}

static int bwd_wspace_alloc (bwd_wspace *mm, int T, int k)
{
    mm->B = gretl_matrix_block_new(&mm->iR, k, k,
                                   &mm->g, k, 1,
                                   &mm->b, k, 1,
                                   &mm->V, k, k,
                                   &mm->e, T, 1,
                                   NULL);
    if (mm->B == NULL) {
        return E_ALLOC;
    } else {
        return 0;
    }
}

static void bwd_wspace_shrink (bwd_wspace *mm)
{
    mm->g->rows -= 1;
    mm->b->rows -= 1;
    mm->V->rows -= 1;
    mm->V->cols -= 1;
    mm->iR->rows -= 1;
    mm->iR->cols -= 1;
}

static void bwd_wspace_free (bwd_wspace *mm)
{
    gretl_matrix_block_destroy(mm->B);
}

/* Drop/cut a single column from matrix @m */

static void matrix_drop_column (gretl_matrix *m, int j)
{
    size_t sz = m->rows * sizeof(double);
    int ntrail = m->cols - j - 1;
    double *dest = m->val + j * m->rows;
    double *src = dest + m->rows;

    if (m->info != NULL) {
        /* don't leak column names */
        gretl_matrix_destroy_info(m);
    }

    memmove(dest, src, ntrail * sz);
    m->cols -= 1;
}

static void matrix_drop_index (gretl_matrix *m, int idx)
{
    size_t sz = m->rows * sizeof(double);
    int ntrail = m->cols - idx - 1;
    double *dest = m->val + idx * m->rows;
    double *src = dest + m->rows;
    int n, j;

    if (m->info != NULL) {
        /* don't leak column names */
        gretl_matrix_destroy_info(m);
    }

    /* drop the relevant column */
    memmove(dest, src, ntrail * sz);
    m->cols -= 1;

    n = m->rows * m->cols - idx - 1;
    sz = n * sizeof(double);

    /* drop corresponding row */
    dest = m->val + idx;
    for (j=0; j<m->cols; j++) {
        src = dest + 1;
        memmove(dest, src, sz);
        dest += m->rows - 1;
        n -= m->rows;
        sz = n * sizeof(double);
    }
    m->rows -= 1;
}

static double get_info_crit (double ssr, int T, int k, int crit)
{
    /* using AIC, BIC or HQC */
    double l0 = log(2*M_PI);
    double ll, C0 = 0;

    if (crit == C_BIC) {
        C0 = k * log((double) T);
    } else if (crit == C_HQC) {
        C0 = 2 * k * log(log((double) T));
    }

    ll = -0.5 * T * (l0 + log(ssr/T) + 1);

    return crit == C_AIC ? 2.0 * (k - ll) : C0 - 2 * ll;
}

static double best_ssr2crit (int *best,
                             const gretl_matrix *ssr,
                             int T, int k, int crit)
{
    int j, n = ssr->cols;
    double ssr_min = 1.0e200;

    *best = 0;

    /* find minimum SSR first */
    for (j=0; j<n; j++) {
        if (ssr->val[j] < ssr_min) {
            *best = j;
            ssr_min = ssr->val[j];
        }
    }

    if (crit == C_SSR) {
        return ssr_min;
    } else {
        return get_info_crit(ssr_min, T, k, crit);
    }
}

static double ssr2crit (double ssr, int T, int k, int crit)
{
    if (crit == C_SSR) {
        return ssr;
    } else {
        return get_info_crit(ssr, T, k, crit);
    }
}

static int qr_augment (gretl_matrix *Q,
                       gretl_matrix *R,
                       const gretl_matrix *Z,
                       const gretl_matrix *e,
                       fwd_wspace *mm,
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

/* The following function fixes up the trailing columns of the QR
   decomposition of the matrix @A, following deletion of column @jmin of
   @A when @jmin was not the last column. We assume that @Q contains
   the correct decomposition of the submatrix of @A to the left of the
   @jmin.
*/

void qr_fixup (gretl_matrix *Q,
               gretl_matrix *R,
               gretl_matrix *A,
               int jmin)
{
    double *aj;
    size_t sz;
    double tmp;
    int n = A->rows;
    int m = R->rows;
    int i, j, k;

    aj = malloc(n * sizeof *aj);
    sz = n * sizeof(double);

    /* complete Q */
    for (j=jmin; j<m; j++) {
        memcpy(aj, A->val + j*n, sz);
        for (i=0; i<j; i++) {
            tmp = 0.0;
            for (k=0; k<n; k++) {
                tmp += aj[k] * gretl_matrix_get(Q, k, i);
            }
            for (k=0; k<n; k++) {
                aj[k] -= tmp * gretl_matrix_get(Q, k, i);
            }
        }
        tmp = 0.0;
        for (i=0; i<n; i++) {
            tmp += aj[i] * aj[i];
        }
        tmp = sqrt(tmp);
        for (i=0; i<n; i++) {
            gretl_matrix_set(Q, i, j, aj[i]/tmp);
        }
    }

    /* complete R */
    for (i=jmin; i<m; i++) {
        for (j=i; j<m; j++) {
            tmp = 0.0;
            for (k=0; k<n; k++) {
                tmp += gretl_matrix_get(A, k, j) * gretl_matrix_get(Q, k, i);
            }
            gretl_matrix_set(R, i, j, tmp);
        }
    }

    free(aj);
}

static int qr_reduce (gretl_matrix *Q,
                      gretl_matrix *R,
                      gretl_matrix *A,
                      int delcol,
                      const gretl_matrix *y,
                      bwd_wspace *mm,
                      double *ssr)
{
    int i, k = R->rows;
    int T = Q->rows;
    double s2;
    int err = 0;

    matrix_drop_column(Q, delcol);
    matrix_drop_column(A, delcol);
    matrix_drop_index(R, delcol);

    if (delcol < Q->cols) {
        qr_fixup(Q, R, A, delcol);
    }

    bwd_wspace_shrink(mm);

    gretl_matrix_copy_values(mm->iR, R);
    err = gretl_invert_triangular_matrix(mm->iR, 'U');
    if (err) {
        return err;
    }

    /* g = Q'y */
    gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
                              y, GRETL_MOD_NONE,
                              mm->g, GRETL_MOD_NONE);
    /* b = inv(R) * g */
    gretl_matrix_multiply(mm->iR, mm->g, mm->b);

    /* e = residuals */
    gretl_matrix_copy_values(mm->e, y);
    gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                              mm->g, GRETL_MOD_NONE,
                              mm->e, GRETL_MOD_DECREMENT);

    /* V [= (X'X)^{-1}] = inv(R)*inv(R)' */
    gretl_matrix_multiply_mod(mm->iR, GRETL_MOD_NONE,
                              mm->iR, GRETL_MOD_TRANSPOSE,
                              mm->V, GRETL_MOD_NONE);

    *ssr = 0.0;
    for (i=0; i<T; i++) {
        *ssr += mm->e->val[i] * mm->e->val[i];
    }
    s2 = *ssr / (T - k);

    /* g <- std errors */
    for (i=0; i<k; i++) {
        mm->g->val[i] = sqrt(s2 * gretl_matrix_get(mm->V, i, i));
    }

    return err;
}

static const char *cstrs[] = {
    "AIC", "BIC", "HQC", "SSR"
};

static const char *crit_string (int crit)
{
    return cstrs[crit];
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
    fwd_wspace mm = {0};
    int conv = 0;
    int added = 0;
    int best = 0;
    int bpos;
    double cur, prev;
    int *xlist = NULL;
    int *aux = NULL;
    int *ret = NULL;
    int yvar;
    int t1 = pmod->t1;
    int t2 = pmod->t2;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int i, nz;

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
    cstr = crit_string(crit);
    fwd_wspace_alloc(&mm, Q->rows, nz);

    if (verbose) {
        pprintf(prn, "\n %-*s %s = %#g\n", addlen + namelen + 1,
                _("Baseline"), cstr, prev);
    }

    while (!conv && added < nz) {
        *err = qr_augment(Q, R, mZ, e, &mm, crit, &best, &cur);
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

        if (crit == C_SSR) {
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
                pprintf(prn, " [%-*s %s = %#g]\n", namelen + addlen,
                        dset->varname[aux[best+1]], cstr, cur);
            } else {
                pprintf(prn, " %s %-*s %s = %#g\n", _("Add"), namelen,
                        dset->varname[aux[best+1]], cstr, cur);
            }
        }
        if (!conv) {
            bpos = best + 1;
            matrix_drop_column(mZ, best);
            fwd_wspace_shrink(&mm);
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
    fwd_wspace_free(&mm);

    if (!*err && added == 0) {
        *err = E_NOADD;
    }
    if (*err) {
        if (verbose) {
            pputc(prn, '\n');
        }
        free(ret);
        ret = NULL;
    }

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
                         int ifc, int k)
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

    return ret;
}

int *backward_stepwise (MODEL *pmod,
                        const int *zlist,
                        DATASET *dset,
                        int crit,
                        int verbose,
                        int droplen,
                        int namelen,
                        PRN *prn,
                        int *err)
{
    const char *cstr;
    gretl_matrix *Q;
    gretl_matrix *R;
    gretl_matrix *y;
    gretl_matrix *A;
    int *xlist = NULL;
    bwd_wspace mm = {0};
    double cur, prev;
    int ifc = pmod->ifc;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int t1 = pmod->t1;
    int t2 = pmod->t2;
    int yvar = pmod->list[1];
    int conv = 0;
    int trycol, delvar;
    int nz, dropped = 0;

    prev = ssr2crit(pmod->ess, T, k, crit);
    xlist = gretl_model_get_x_list(pmod);
    trycol = tval_min_pos(pmod->coeff, pmod->sderr,
                          xlist, zlist, ifc, k);

    y = gretl_vector_from_series(dset->Z[yvar], t1, t2);
    Q = gretl_matrix_data_subset(xlist, dset, t1, t2, M_MISSING_ERROR, err);
    R = gretl_zero_matrix_new(k, k);
    A = gretl_matrix_copy(Q);

    if (!*err && (y == NULL || Q == NULL || R == NULL || A == NULL)) {
        *err = E_ALLOC;
    }
    if (*err) {
        goto bailout;
    }

    gretl_matrix_QR_decomp(Q, R);

    nz = zlist != NULL ? zlist[0] : xlist[0] - ifc;
    cstr = crit_string(crit);
    bwd_wspace_alloc(&mm, Q->rows, R->rows);

    if (verbose) {
        pprintf(prn, " %-*s %s = %#g\n", droplen + namelen + 1,
                _("Baseline"), cstr, prev);
    }

    while (!conv && dropped < nz) {
        double ssr;

        delvar = xlist[trycol+1];
        *err = qr_reduce(Q, R, A, trycol, y, &mm, &ssr);
        if (*err) {
            fprintf(stderr, "error %d in qr_reduce\n", *err);
            break;
        }
        k = R->rows;
        cur = ssr2crit(ssr, T, k, crit);
        conv = cur > prev;
        if (verbose) {
            if (conv && dropped < nz) {
                pprintf(prn, " [%-*s %s = %#g]\n", namelen + droplen,
                        dset->varname[delvar], cstr, cur);
            } else {
                pprintf(prn, " %s %-*s %s = %#g\n", _("Drop"), namelen,
                        dset->varname[delvar], cstr, cur);
            }
        }
        if (!conv) {
            gretl_list_delete_at_pos(xlist, trycol+1);
            trycol = tval_min_pos(mm.b->val, mm.g->val, xlist, zlist, ifc, k);
            dropped++;
            prev = cur;
        }
    }

 bailout:

    gretl_matrix_free(Q);
    gretl_matrix_free(R);
    gretl_matrix_free(y);
    gretl_matrix_free(A);
    bwd_wspace_free(&mm);

    if (!*err && dropped == 0) {
        *err = E_NOOMIT;
    }
    if (*err) {
        if (verbose) {
            pputc(prn, '\n');
        }
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

/* Compose the final OLS list based on the original model plus @xlist,
   which contains the regressors that survived automatic elimination.
*/

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
                    int crit,
                    double alpha,
                    DATASET *dset,
                    gretlopt opt,
                    PRN *prn)
{
    MODEL model;
    int *best = NULL;
    int ols_done = 0;
    int namelen = 0;
    int err;

    err = forward_stepwise_check_zlist(pmod, zlist, dset, &namelen);

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

/* Implement "omit --auto=..." using backward stepwise procedure.
   Note: @alpha is unused at present, since backward_stepwise() above
   doesn't handle the p-value criterion for omission. (That's handled
   in lib/src/compare.c.)
*/

MODEL stepwise_omit (MODEL *pmod,
                     const int *zlist,
                     int crit,
                     double alpha,
                     DATASET *dset,
                     gretlopt opt,
                     PRN *prn)
{
    MODEL model;
    int *xlist = NULL;
    int ols_done = 0;
    int namelen = 0;
    int err;

    err = backward_stepwise_check_zlist(pmod, zlist, dset, &namelen);

    if (!err) {
        int verbose = (opt & OPT_Q)? 0 : 1;
        int droplen = verbose ? g_utf8_strlen(_("Drop"), -1) : 0;

	if (verbose) {
	    // pputc(prn, '\n');
	    pprintf(prn, _("Sequential elimination using %s"), crit_string(crit));
	    pputs(prn, "\n\n");
	}
        xlist = backward_stepwise(pmod, zlist, dset, crit,
                                  verbose, droplen, namelen + 2,
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
