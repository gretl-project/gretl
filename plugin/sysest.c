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
#include "version.h"
#include "gretl_matrix.h"
#include "libset.h"
#include "system.h"
#include "tsls.h"
#include "sysml.h"

#define SDEBUG 0

#define sys_ols_ok(s) (s->method == SYS_METHOD_SUR || \
                       s->method == SYS_METHOD_OLS || \
                       s->method == SYS_METHOD_WLS)

/* Insert the elements of sub-matrix @M, multiplied by @scale, in the
   appropriate position within the big matrix @X. We're building
   a symmetric matrix here so we insert off-diagonal elements both
   above and below the diagonal.
*/

static void
insert_sys_X_block (gretl_matrix *X, const gretl_matrix *M,
                    int startrow, int startcol, double scale)
{
    int r, c, i, j;
    double mij;

    for (i=0; i<M->rows; i++) {
        r = startrow + i;
        for (j=0; j<M->cols; j++) {
            c = startcol + j;
            mij = gretl_matrix_get(M, i, j);
            gretl_matrix_set(X, r, c, mij * scale);
        }
    }

    if (startrow != startcol) {
        for (i=0; i<M->rows; i++) {
            c = startrow + i;
            for (j=0; j<M->cols; j++) {
                r = startcol + j;
                mij = gretl_matrix_get(M, i, j);
                gretl_matrix_set(X, r, c, mij * scale);
            }
        }
    }
}

/* Retrieve the data placed on @pmod in the course of 2sls estimation:
   for endogenous regressors these values are the first-stage fitted
   values, otherwise they're just the original data.  Note that @endog
   is a mask with nonzero values identifying the endogenous
   regressors, and the special array @X contains a column for each
   endogenous variable.
*/

double *model_get_Xi (const MODEL *pmod, DATASET *dset, int i)
{
    gretl_matrix *endog = gretl_model_get_data(pmod, "endog");
    double *xi = NULL;

    if (endog == NULL || endog->val[i] == 0) {
        /* use the original data */
        xi = dset->Z[pmod->list[i+2]];
    } else {
        double **X = gretl_model_get_data(pmod, "tslsX");

        if (X != NULL) {
            /* find and return the correct column of X */
            int j, k = 0;

            for (j=0; j<i; j++) {
                if (endog->val[j] != 0) k++;
            }
            xi = X[k];
        }
    }

    return xi;
}

/* retrieve the special k-class data wanted for LIML estimation: these
   represent a further transformation of the data placed on the model
   via 2sls estimation.
*/

static int make_liml_X_block (gretl_matrix *X, const MODEL *pmod,
                              DATASET *dset, int t1)
{
    const double *Xi;
    int i, t, err = 0;

    X->cols = pmod->ncoeff;

    for (i=0; i<X->cols && !err; i++) {
        Xi = model_get_Xi(pmod, dset, i);
        if (Xi == NULL) {
            err = E_DATA;
        } else {
            for (t=0; t<X->rows; t++) {
                gretl_matrix_set(X, t, i, Xi[t+t1]);
            }
        }
    }

    return err;
}

/* construct the X data block pertaining to a specific equation, using
   either the original data or fitted values from regression on a set
   of instruments
*/

static int
make_sys_X_block (gretl_matrix *X, const MODEL *pmod,
                  DATASET *dset, int t1, int method)
{
    const double *Xi;
    int i, t, err = 0;

    X->cols = pmod->ncoeff;

    for (i=0; i<X->cols && !err; i++) {
        if (method == SYS_METHOD_3SLS ||
            method == SYS_METHOD_FIML ||
            method == SYS_METHOD_TSLS) {
            Xi = model_get_Xi(pmod, dset, i);
        } else {
            Xi = dset->Z[pmod->list[i+2]];
        }
        if (Xi == NULL) {
            err = E_DATA;
        } else {
            for (t=0; t<X->rows; t++) {
                gretl_matrix_set(X, t, i, Xi[t+t1]);
            }
        }
    }

    return err;
}

/* Populate the cross-equation covariance matrix, @S, based on the
   per-equation residuals: if @do_diag is non-zero we also want to compute
   the Breusch-Pagan test for diagonality of the matrix. For the robust
   variant see Halunga, Orme and Yamagata, "A heteroskasticity robust
   Breusch-Pagan test for contemporaneous correlation...", Journal of
   Econometrics, 198 (2017) pp. 209-230.
*/

static int
gls_sigma_from_uhat (equation_system *sys, gretl_matrix *S,
                     int do_diag)
{
    int geomean = system_vcv_geomean(sys);
    double *rdenom = NULL;
    double eti, etj, sij, dij;
    int m = sys->neqns;
    int robust = 0;
    int i, j, t, k;

    if (do_diag && (sys->flags & SYSTEM_ROBUST)) {
        /* denominator of robust B-P test */
        rdenom = malloc((m * m - m) / 2 * sizeof *rdenom);
        if (rdenom != NULL) {
            robust = 1;
        }
    }

    k = 0;
    for (i=0; i<m; i++) {
        for (j=i; j<m; j++) {
            sij = dij = 0.0;
            for (t=0; t<sys->T; t++) {
                eti = gretl_matrix_get(sys->E, t, i);
                etj = gretl_matrix_get(sys->E, t, j);
                /* increment sum of cross-products */
                sij += eti * etj;
                if (i != j && robust) {
                    /* also sum of cross-products of squares */
                    dij += eti * etj * eti * etj;
                }
            }
            if (i != j && robust) {
                rdenom[k++] = dij;
            }
            gretl_matrix_set(S, i, j, sij);
            if (j != i) {
                gretl_matrix_set(S, j, i, sij);
            }
        }
    }
    if (do_diag) {
        /* compute B-P test statistic */
        double sii, sjj;

        sys->diag_test = 0.0;

        k = 0;
        for (i=0; i<m-1; i++) {
            sii = gretl_matrix_get(S, i, i);
            for (j=i+1; j<m; j++) {
                sij = gretl_matrix_get(S, i, j);
                sjj = gretl_matrix_get(S, j, j);
                if (robust) {
                    /* cumulate \hat{gamma}_{ij}^2 */
                    sys->diag_test += (sij * sij) / rdenom[k++];
                } else {
                    /* cumulate \hat{rho}_{ij}^2 */
                    sys->diag_test += (sij * sij) / (sii * sjj);
                }
            }
        }
        if (robust) {
            free(rdenom);
        } else {
            /* supply the missing factor of T */
            sys->diag_test *= sys->T;
        }
    }

    /* scale the S matrix */
    if (geomean) {
        for (j=0; j<S->cols; j++) {
            for (i=j; i<S->rows; i++) {
                sij = gretl_matrix_get(S, i, j);
                sij /= system_vcv_denom(sys, i, j);
                gretl_matrix_set(S, i, j, sij);
                if (i != j) {
                    gretl_matrix_set(S, j, i, sij);
                }
            }
        }
    } else {
        gretl_matrix_divide_by_scalar(S, sys->T);
    }

    return 0;
}

static void maybe_correct_model_dfd (equation_system *sys,
                                     MODEL *pmod)
{
    int nr = system_n_restrictions(sys);

    if (nr > 0) {
        /* adjust the df correction for the number of
           restrictions, averaged across the equations
        */
        double avr = nr / (double) sys->neqns;

        pmod->dfd = sys->T - pmod->ncoeff + floor(avr);
    }
}

/* Compute residuals and related quantities, for all cases
   other than FIML
*/

static void
sys_resids (equation_system *sys, int eq, const DATASET *dset)
{
    MODEL *pmod = sys->models[eq];
    int yno = pmod->list[1];
    double yh;
    int i, t;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
        yh = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            yh += pmod->coeff[i] * dset->Z[pmod->list[i+2]][t];
        }
        pmod->yhat[t] = yh;
        pmod->uhat[t] = dset->Z[yno][t] - yh;
        /* for cross-equation vcv */
        gretl_matrix_set(sys->E, t - pmod->t1, pmod->ID, pmod->uhat[t]);
        pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    /* df correction? */
    if (system_want_df_corr(sys)) {
        maybe_correct_model_dfd(sys, pmod);
        pmod->sigma = sqrt(pmod->ess / pmod->dfd);
    } else {
        pmod->sigma = sqrt(pmod->ess / pmod->nobs);
    }

    if (pmod->ifc && pmod->tss > 0) {
        /* R-squared based on actual-fitted correlation */
        double r;

        pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, dset->Z[yno],
                                   pmod->yhat);
        r = 1 - pmod->rsq;
        pmod->adjrsq = 1.0 - (r * (pmod->nobs - 1) / pmod->dfd);
    } else {
        pmod->rsq = pmod->adjrsq = NADBL;
    }
}

/* use the per-equation residual standard deviations for scaling
   the covariance matrix, block by diagonal block
*/

static void
single_eqn_scale_vcv (MODEL *pmod, gretl_matrix *V, int offset)
{
    double vij, s2 = pmod->sigma * pmod->sigma;
    int i, j, ii, jj;

    for (i=0; i<pmod->ncoeff; i++) {
        ii = i + offset;
        for (j=i; j<pmod->ncoeff; j++) {
            jj = j + offset;
            vij = gretl_matrix_get(V, ii, jj);
            vij *= s2;
            gretl_matrix_set(V, ii, jj, vij);
            gretl_matrix_set(V, jj, ii, vij);
        }
    }
}

/* write results from system estimation into the individual
   model structs
*/

static void transcribe_sys_results (equation_system *sys,
                                    const DATASET *dset,
                                    int do_iters)
{
    int do_bdiff = (sys->method == SYS_METHOD_3SLS && do_iters);
    double bij, oldb, bnum = 0, bden = 0;
    int offset = 0;
    int i, j, k;

    for (i=0; i<sys->neqns; i++) {
        MODEL *pmod = sys->models[i];

        /* update the coefficients */
        for (j=0; j<pmod->ncoeff; j++) {
            k = j + offset;
            bij = gretl_vector_get(sys->b, k);
            if (do_bdiff) {
                oldb = pmod->coeff[j];
                bnum += (bij - oldb) * (bij - oldb);
                bden += oldb * oldb;
            }
            pmod->coeff[j] = bij;
        }
        /* update residuals (and sigma-hat) */
        sys_resids(sys, i, dset);
        /* for single-equation methods, vcv needs scaling by
           an estimate of the residual variance
        */
        if (sys->method == SYS_METHOD_OLS ||
            sys->method == SYS_METHOD_TSLS ||
            sys->method == SYS_METHOD_LIML) {
            single_eqn_scale_vcv(pmod, sys->vcv, offset);
        }
        /* update standard errors */
        for (j=0; j<pmod->ncoeff; j++) {
            k = j + offset;
            pmod->sderr[j] = sqrt(gretl_matrix_get(sys->vcv, k, k));
        }
        offset += pmod->ncoeff;
    }

    if (do_bdiff) {
        sys->bdiff = sqrt(bnum / bden);
    }
}

/* compute SUR, 3SLS or LIML parameter estimates (or restricted OLS,
   TSLS, WLS)
*/

static int
calculate_sys_coeffs (equation_system *sys, const DATASET *dset,
                      gretl_matrix *X, gretl_matrix *y,
                      int mk, int nr, int do_iters)
{
    gretl_matrix *V;
    int posdef = 0;
    int err = 0;

    V = gretl_matrix_copy(X);
    if (V == NULL) {
        return E_ALLOC;
    }

#if SDEBUG
    gretl_matrix_print(X, "sys X");
    gretl_matrix_print(y, "sys y");
#endif

    /* calculate the coefficients */

    if (sys->R == NULL) {
        /* hopefully X may be positive definite */
        err = gretl_cholesky_decomp_solve(X, y);
        if (!err) {
            posdef = 1;
            err = gretl_cholesky_invert(X);
        } else {
            /* nope: fall back to the LU solver */
            gretl_matrix_copy_values(X, V);
            err = 0;
        }
    }

    if (!posdef) {
        err = gretl_LU_solve_invert(X, y);
    }

    if (!err) {
        /* save the coefficient vector and covariance matrix */
        gretl_matrix *b = gretl_matrix_copy(y);

        if (b == NULL) {
            err = E_ALLOC;
        } else {
            system_attach_coeffs(sys, b);
            gretl_matrix_copy_values(V, X);
            system_attach_vcv(sys, V);
            /* transcribe stuff to the included models */
            transcribe_sys_results(sys, dset, do_iters);
        }
    }

    if (err) {
        gretl_matrix_free(V);
    }

    return err;
}

static int *
system_model_list (equation_system *sys, int i, int *freeit)
{
    int *list = NULL;

    *freeit = 0;

    if (sys_ols_ok(sys)) {
        return system_get_list(sys, i);
    }

    if (sys->method != SYS_METHOD_FIML) {
        list = system_get_list(sys, i);
    }

    if (sys->method == SYS_METHOD_3SLS ||
        sys->method == SYS_METHOD_TSLS ||
        sys->method == SYS_METHOD_LIML) {
        /* Is the list already in ivreg form?
           If not, ignore it. */
        if (list != NULL && !gretl_list_has_separator(list)) {
            list = NULL;
        }
    }

    if ((sys->method == SYS_METHOD_FIML ||
         sys->method == SYS_METHOD_LIML ||
         sys->method == SYS_METHOD_3SLS ||
         sys->method == SYS_METHOD_TSLS)
        && list == NULL) {
        list = compose_ivreg_list(sys, i);
        *freeit = 1;
    }

    return list;
}

/* Hansen-Sargan overidentification test for the system as a whole, as
   in Davidson and MacKinnon, ETM: p. 511 and equation (12.25) for the
   case of SUR; p. 532 and equation (12.61) for the case of 3SLS.  See
   also D & M, Estimation and Inference in Econometrics, equation
   (18.60), for a more computation-friendly statement of the criterion
   function.
*/

static int hansen_sargan_test (equation_system *sys,
                               const DATASET *dset)
{
    const int *exlist = system_get_instr_vars(sys);
    const double *Wi, *Wj;
    int nx = exlist[0];
    int m = sys->neqns;
    int T = sys->T;
    int df = system_get_overid_df(sys);
    gretl_matrix_block *B;
    gretl_matrix *WTW, *eW, *tmp;
    double x, X2;
    int i, j, t;
    int err = 0;

    if (df <= 0) {
        return 1;
    }

    B = gretl_matrix_block_new(&WTW, nx, nx,
                               &eW, m, nx,
                               &tmp, m, nx,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    /* construct W-transpose W in WTW */
    for (i=0; i<nx; i++) {
        Wi = dset->Z[exlist[i+1]] + sys->t1;
        for (j=i; j<nx; j++) {
            Wj = dset->Z[exlist[j+1]] + sys->t1;
            x = 0.0;
            for (t=0; t<sys->T; t++) {
                x += Wi[t] * Wj[t];
            }
            gretl_matrix_set(WTW, i, j, x);
            if (i != j) {
                gretl_matrix_set(WTW, j, i, x);
            }
        }
    }

    err = gretl_invert_symmetric_matrix(WTW);
#if SDEBUG
    fprintf(stderr, "hansen_sargan: on invert, err=%d\n", err);
#endif
    if (err) {
        sys->X2 = NADBL;
        goto bailout;
    }

    /* set up vectors of SUR or 3SLS residuals, transposed,
       times W; these are stacked in the m * nx matrix eW
    */
    for (i=0; i<m; i++) {
        for (j=0; j<nx; j++) {
            Wj = dset->Z[exlist[j+1]] + sys->t1;
            x = 0.0;
            for (t=0; t<T; t++) {
                x += gretl_matrix_get(sys->E, t, i) * Wj[t];
            }
            gretl_matrix_set(eW, i, j, x);
        }
    }

    /* multiply these vectors into (WTW)^{-1} */
    for (i=0; i<m; i++) {
        for (j=0; j<nx; j++) {
            x = 0.0;
            for (t=0; t<nx; t++) {
                x += gretl_matrix_get(eW, i, t) *
                    gretl_matrix_get(WTW, t, j);
            }
            gretl_matrix_set(tmp, i, j, x);
        }
    }

    /* cumulate the Chi-square value */
    X2 = 0.0;
    for (i=0; i<m; i++) {
        for (j=0; j<m; j++) {
            x = 0.0;
            for (t=0; t<nx; t++) {
                x += gretl_matrix_get(tmp, i, t) *
                    gretl_matrix_get(eW, j, t); /* transposed */
            }
            X2 += gretl_matrix_get(sys->S, i, j) * x;
        }
    }

#if SDEBUG
    fprintf(stderr, "Hansen-Sargan: Chi-square(%d) = %g (p-value %g)\n",
            df, X2, chisq_cdf_comp(df, X2));
#endif
    sys->X2 = X2;

 bailout:

    gretl_matrix_block_destroy(B);

    return err;
}

static int basic_system_allocate (equation_system *sys,
                                  int mk, int nr,
                                  gretl_matrix **X,
                                  gretl_matrix **y)
{
    int m = sys->neqns;
    int T = sys->T;
    int ldx = mk + nr;

    /* allocate a model for each stochastic equation */
    sys->models = gretl_model_array_new(m);
    if (sys->models == NULL) {
        return E_ALLOC;
    }

    sys->E = gretl_matrix_alloc(T, m);
    if (sys->E == NULL) {
        return E_ALLOC;
    }

    sys->S = gretl_matrix_alloc(m, m);
    if (sys->S == NULL) {
        return E_ALLOC;
    }

    if (X != NULL) {
        *X = gretl_matrix_alloc(ldx, ldx);
        if (*X == NULL) {
            return E_ALLOC;
        }
    }

    if (y != NULL) {
        *y = gretl_column_vector_alloc(ldx);
        if (*y == NULL) {
            return E_ALLOC;
        }
    }

    return 0;
}

/* compute log-likelihood for iterated SUR estimator */

double sur_loglik (equation_system *sys)
{
    int m = sys->neqns;
    int T = sys->T;
    gretl_matrix *tmp;
    double ldet;
    int err = 0;

    tmp = gretl_matrix_alloc(m, m);
    if (tmp == NULL) {
        return NADBL;
    }

    gls_sigma_from_uhat(sys, tmp, 0);
    ldet = gretl_vcv_log_determinant(tmp, &err);

    if (na(ldet)) {
        sys->ll = NADBL;
    } else {
        sys->ll = -(m * T / 2.0) * (LN_2_PI + 1.0);
        sys->ll -= (T / 2.0) * ldet;
    }

    gretl_matrix_free(tmp);

    return sys->ll;
}

/* When estimating with a specified set of linear restrictions,
   Rb = q, augment the X matrix with R and R-transpose.
*/

static void
augment_X_with_restrictions (gretl_matrix *X, int mk,
                             equation_system *sys)
{
    int i, j, nr = sys->R->rows;

    /* place the R matrix */
    insert_sys_X_block(X, sys->R, mk, 0, 1.0);

    /* zero the bottom right-hand block */
    for (i=mk; i<mk+nr; i++) {
        for (j=mk; j<mk+nr; j++) {
            gretl_matrix_set(X, i, j, 0.0);
        }
    }
}

/* When estimating with a specified set of linear restrictions,
   Rb = q, augment the y matrix with q.
*/

static int
augment_y_with_restrictions (gretl_matrix *y, int mk, int nr,
                             equation_system *sys)
{
    int i;

    if (sys->q == NULL) {
        return 1;
    }

    for (i=0; i<nr; i++) {
        gretl_vector_set(y, mk + i, gretl_vector_get(sys->q, i));
    }

    return 0;
}

#define SYS_MAX_ITER 250
#define SYS_LL_TOL 1.0e-12
#define SYS_BDIFF_TOL 1.0e-9

/* check for convergence of iteration: we use the change in the
   log-likelihood when iterating SUR to the ML solution, or
   a measure of the change in the coefficients when iterating
   3SLS
*/

static int sys_converged (equation_system *sys, double *llbak,
                          gretlopt opt, PRN *prn, int *err)
{
    double crit, tol = 0.0;
    int met = 0;

    if (sys->method == SYS_METHOD_SUR ||
        sys->method == SYS_METHOD_WLS) {
        double ll = sur_loglik(sys);

        tol = SYS_LL_TOL;
        crit = ll - *llbak;

        if (opt & OPT_V) {
            pprintf(prn, _("Iteration %3d, ll = %#.8g\n"), sys->iters, ll);
        }
        if (crit <= tol) {
            met = 1;
        } else if (sys->iters < SYS_MAX_ITER) {
            *llbak = ll;
        }
    } else if (sys->method == SYS_METHOD_3SLS) {
        tol = SYS_BDIFF_TOL;
        crit = sys->bdiff;
        if (opt & OPT_V) {
            pprintf(prn, _("Iteration %3d, criterion = %g\n"), sys->iters, crit);
        }
        if (crit <= tol) {
            met = 1;
        }
    }

    if (met && tol > 0 && (opt & OPT_V)) {
        pprintf(prn, _("Tolerance of %g is met\n"), tol);
    }

    if (!met && sys->iters >= SYS_MAX_ITER) {
        pprintf(prn, _("Reached %d iterations without meeting "
                "tolerance of %g\n"), sys->iters, tol);
        *err = E_NOCONV;
    }

    return met;
}

static void clean_up_models (equation_system *sys)
{
    int i;

    if (sys->models == NULL) {
        /* "can't happen" */
        return;
    }

    sys->ess = 0.0;

    for (i=0; i<sys->neqns; i++) {
        sys->ess += sys->models[i]->ess;
        if (sys->method == SYS_METHOD_3SLS ||
            sys->method == SYS_METHOD_FIML ||
            sys->method == SYS_METHOD_TSLS ||
            sys->method == SYS_METHOD_LIML) {
            tsls_free_data(sys->models[i]);
        }
        if (sys->neqns > 1) {
            gretl_model_free(sys->models[i]);
        }
    }

    if (sys->flags & SYSTEM_LIML1) {
        /* ivreg (single-equation LIML) */
        sys->models[0]->rho = sys->models[0]->dw = NADBL;
    } else {
        free(sys->models);
        sys->models = NULL;
    }
}

#define REJECT_COLLINEAR 1

#if REJECT_COLLINEAR

static int perfect_collinearity_check (MODEL *pmod,
                                       DATASET *dset,
                                       int i)
{
    const int *d1 = gretl_model_get_list(pmod, "droplist");
    const int *d2 = gretl_model_get_list(pmod, "inst_droplist");

    if (d1 != NULL) {
        gretl_errmsg_sprintf(_("Equation %d exhibits perfect collinearity.\n"
                             "The regressor %s cannot be included. Please "
                             "respecify the system."), i, dset->varname[d1[1]]);
        return E_SINGULAR;
    } else if (d2 != NULL) {
        gretl_errmsg_sprintf(_("Equation %d exhibits perfect collinearity.\n"
                             "The instrument %s cannot be included. Please "
                             "respecify the system."), i, dset->varname[d2[1]]);
        return E_SINGULAR;
    } else {
        return 0;
    }
}

#else

/* Given a record of redundant regressors or instruments dropped
   at the initial OLS or TSLS stage, namely @droplist, purge these
   series IDs from both the relevant system lists. The @mk pointer
   argument tells us whether we're looking at regular regressors
   (mk non-NULL) or instruments.
*/

static int drop_redundant_variables (equation_system *sys,
                                     const int *droplist,
                                     int eqn, int *mk)
{
    int i, j, di, pmax, pmin = 0;
    int *eqnlist = sys->lists[eqn];
    int insts = (mk == NULL);
    int err = 0;

    if (insts) {
        for (i=1; i<=droplist[0]; i++) {
            di = droplist[i];
            j = in_gretl_list(sys->ilist, di);
            if (j > 0) {
                gretl_list_delete_at_pos(sys->ilist, j);
            } else {
                err = 1;
            }
        }
    }

    j = gretl_list_separator_position(eqnlist);

    if (insts) {
        /* instruments */
        if (j > 0) {
            pmin = j + 1;
            pmax = eqnlist[0];
        }
    } else {
        pmin = 2;
        pmax = (j > 0) ? j - 1 : eqnlist[0];
    }

    if (pmin > 0) {
        for (i=1; i<=droplist[0]; i++) {
            di = droplist[i];
            for (j=pmax; j>=pmin; j--) {
                if (eqnlist[j] == di) {
                    gretl_list_delete_at_pos(eqnlist, j);
                    if (mk != NULL) {
                        mk -= 1;
                    }
                    break;
                }
            }
        }
    }

    return err;
}

#endif /* REJECT_COLLINEAR or not */

/* options to be passed in running initial 2SLS */

static gretlopt sys_tsls_opt (const equation_system *sys,
                              gretlopt opt)
{
    gretlopt tsls_opt = (OPT_E | OPT_A);

    if (!(sys->flags & SYSTEM_DFCORR)) {
        tsls_opt |= OPT_N; /* suppress df correction */
    }

    if (sys->flags & SYSTEM_LIML1) {
        tsls_opt |= OPT_H; /* add "hatlist" of instrumented vars */
    }

    return tsls_opt;
}

/* options to be passed in running initial OLS; @nr is
   the number of restrictions imposed
*/

static gretlopt sys_ols_opt (equation_system *sys,
                             int nr, gretlopt opt)
{
    gretlopt ols_opt = OPT_S; /* flag as part of system */

    if (sys->method == SYS_METHOD_OLS && !(sys->flags & SYSTEM_DFCORR)) {
        /* suppress df correction */
        ols_opt |= OPT_N;
    } else if (sys->method == SYS_METHOD_WLS) {
        if (sys->R == NULL && !(opt & OPT_N)) {
            /* equivalent to OLS */
            sys->flags |= SYSTEM_DFCORR;
        } else {
            ols_opt |= OPT_N;
        }
    }

    if (sys->method == SYS_METHOD_OLS && nr == 0) {
        /* initial OLS will supply the estimates */
        if (sys->flags & SYSTEM_ROBUST) {
             ols_opt |= OPT_R;
        }
    } else {
        /* treat initial OLS as auxiliary */
         ols_opt |= OPT_A;
    }

    return  ols_opt;
}

static int allocate_Xi_etc (gretl_matrix **Xi,
                            gretl_matrix **Xj,
                            gretl_matrix **M,
                            int T, int k)
{
    *Xi = gretl_matrix_alloc(T, k);
    *Xj = gretl_matrix_alloc(T, k);
    *M = gretl_matrix_alloc(k, k);

    if (*Xi == NULL || *Xj == NULL || *M == NULL) {
        return E_ALLOC;
    } else {
        return 0;
    }
}

static int ols_data_to_sys (equation_system *sys, int mk)
{
    gretl_matrix *B, *V;
    MODEL *pmod;
    double vij;
    int i, j, k;
    int mi, mj;
    int vi, vj;
    int err = 0;

    B = gretl_matrix_alloc(mk, 1);
    V = gretl_zero_matrix_new(mk, mk);
    if (B == NULL || V == NULL) {
        return E_ALLOC;
    }

    j = vi = vj = 0;
    for (i=0; i<sys->neqns; i++) {
        pmod = sys->models[i];
        if (pmod->vcv == NULL) {
            err = makevcv(pmod, pmod->sigma);
            if (err) {
                break;
            }
        }
        k = pmod->ncoeff;
        for (mi=0; mi<k; mi++) {
            B->val[j++] = pmod->coeff[mi];
            for (mj=0; mj<k; mj++) {
                vij = gretl_model_get_vcv_element(pmod, mi, mj, k);
                gretl_matrix_set(V, vi+mi, vj+mj, vij);
            }
        }
        vi += k;
        vj += k;
    }

    if (err) {
        gretl_matrix_free(B);
        gretl_matrix_free(V);
    } else {
        system_attach_coeffs(sys, B);
        system_attach_vcv(sys, V);
    }

    return err;
}

/* general function that forms the basis for all specific system
   estimators */

int system_estimate (equation_system *sys, DATASET *dset,
                     gretlopt opt, PRN *prn)
{
    int i, j, k, T, t;
    int v, l, mk, krow, nr;
    int orig_t1 = dset->t1;
    int orig_t2 = dset->t2;
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *Xi = NULL;
    gretl_matrix *Xj = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix **pX = NULL;
    gretl_matrix **py = NULL;
    MODEL **models = NULL;
    int method = sys->method;
    double llbak = -1.0e9;
    int single_equation = 0;
    int do_iteration = 0;
    int plain_ols = 0;
    int rsingle = 0;
    int do_diag = 0;
    int err = 0;

    sys->iters = 0;

    if (sys->flags & SYSTEM_ITERATE) {
        do_iteration = 1;
    }

    nr = system_n_restrictions(sys);

    if (method == SYS_METHOD_OLS || method == SYS_METHOD_TSLS ||
        method == SYS_METHOD_LIML || method == SYS_METHOD_WLS) {
        single_equation = 1;
    }

    if ((method == SYS_METHOD_OLS || method == SYS_METHOD_WLS) &&
        sys->R == NULL) {
        plain_ols = 1;
    } else {
        pX = &X;
        py = &y;
    }

    if (nr > 0) {
        if (method == SYS_METHOD_3SLS) {
            /* doing 3SLS with restrictions: we need to obtain
               restricted TSLS estimates as a starting point
            */
            rsingle = 1;
        } else if (method == SYS_METHOD_SUR ||
                   method == SYS_METHOD_WLS) {
            /* doing SUR or WLS with restrictions: we need to obtain
               restricted OLS estimates as a starting point
            */
            rsingle = 1;
        }
    }

    /* get uniform sample starting and ending points and check for
       missing data */
    err = system_adjust_t1t2(sys, dset);
    if (err) {
        return err;
    }

    /* set sample for auxiliary regressions */
    dset->t1 = sys->t1;
    dset->t2 = sys->t2;

    /* max indep vars per equation */
    k = system_max_indep_vars(sys);

    /* total indep vars, all equations */
    mk = system_n_indep_vars(sys);

    /* set sample for auxiliary regressions */
    dset->t1 = sys->t1;
    dset->t2 = sys->t2;

    /* number of observations per series */
    T = sys->T;

    /* allocate models etc */
    err = basic_system_allocate(sys, mk, nr, pX, py);
    if (err) {
        goto cleanup;
    }

    /* convenience pointers */
    models = sys->models;

    if ((method == SYS_METHOD_FIML ||
         method == SYS_METHOD_LIML) && !(opt & OPT_Q)) {
        print_equation_system_info(sys, dset, OPT_H, prn);
    }

    /* First estimate the equations separately (either by OLS or
       TSLS), and put the single-equation residuals into the uhat
       matrix.  Note that at this stage we are not in a position to
       impose any cross-equation restrictions, since we're doing
       straight equation-by-equation estimation.
    */

    for (i=0; i<sys->neqns; i++) {
#if !REJECT_COLLINEAR
        const int *droplist = NULL;
#endif
        int freeit = 0;
        int *list = system_model_list(sys, i, &freeit);
        MODEL *pmod = models[i];
        gretlopt eq_opt;

        if (list == NULL) {
            err = 1;
            break;
        }

        if (sys_ols_ok(sys)) {
            eq_opt = sys_ols_opt(sys, nr, opt);
            *pmod = lsq(list, dset, OLS, eq_opt);
        } else {
            eq_opt = sys_tsls_opt(sys, opt);
            *pmod = tsls(list, dset, eq_opt);
        }

        if (freeit) {
            free(list);
        }

        if ((err = pmod->errcode)) {
            fprintf(stderr, "system_estimate: failed to estimate equation %d: "
                    "err = %d\n", i+1, err);
            break;
        }

#if REJECT_COLLINEAR
        err = perfect_collinearity_check(pmod, dset, i+1);
        if (err) break;
#else
        droplist = gretl_model_get_list(pmod, "droplist");
        if (droplist != NULL) {
            drop_redundant_variables(sys, droplist, i, &mk);
        }
        droplist = gretl_model_get_list(pmod, "inst_droplist");
        if (droplist != NULL) {
            drop_redundant_variables(sys, droplist, i, NULL);
        }
#endif

        pmod->ID = i;
        pmod->aux = AUX_SYS;
        gretl_model_set_int(pmod, "method", method);
        if (eq_opt & OPT_N) {
            gretl_model_set_int(pmod, "asy", 1);
        }

        /* save sigma-squared for an LR test for diagonal
           covariance matrix */
        if (method == SYS_METHOD_SUR && do_iteration && nr == 0) {
            gretl_model_set_double(pmod, "ols_sigma_squared",
                                   pmod->ess / pmod->nobs);
        }

        /* do we want this with @rsingle? */
        for (t=0; t<T; t++) {
            gretl_matrix_set(sys->E, t, i, pmod->uhat[t + sys->t1]);
        }
    }

    if (err) {
        fprintf(stderr, "system_estimate: after single-equation "
                "estimation, err = %d\n", err);
        goto cleanup;
    }

    if (method == SYS_METHOD_LIML) {
        /* compute the minimum eigenvalues and generate the
           k-class data matrices */
        err = liml_driver(sys, dset, prn);
        if (err) goto cleanup;
    }

    if (plain_ols) {
        gls_sigma_from_uhat(sys, sys->S, 1);
        err = ols_data_to_sys(sys, mk);
        goto save_etc;
    }

    /* marker for iterated versions of SUR, WLS, or 3SLS; also for
       loopback in case of restricted 3SLS, where we want to compute
       restricted TSLS estimates first
    */
 iteration_start:

    gls_sigma_from_uhat(sys, sys->S, 0);

    if (method == SYS_METHOD_WLS) {
        gretl_matrix_zero(X);
        err = gretl_invert_diagonal_matrix(sys->S);
    } else if (single_equation || rsingle) {
        gretl_matrix_zero(X);
    } else {
        err = gretl_invert_symmetric_matrix(sys->S);
    }

#if SDEBUG
    fprintf(stderr, "system_estimate: on invert, err=%d\n", err);
    gretl_matrix_print(sys->S, "sys->S");
#endif

    if (!err && Xi == NULL) {
        /* the test against NULL here allows for the possibility
           that we're iterating
        */
        err = allocate_Xi_etc(&Xi, &Xj, &M, T, k);
    }

    if (err) goto cleanup;

    /* form the big stacked X matrix: Xi = data matrix for equation i,
       specified in lists[i]
    */

    krow = 0;
    for (i=0; i<sys->neqns && !err; i++) {
        int kcol = 0;

        err = make_sys_X_block(Xi, models[i], dset, sys->t1, method);

        for (j=0; j<=i && !err; j++) {
            const gretl_matrix *Xk;
            double sij;

            if (i != j) {
                if (single_equation || rsingle) {
                    kcol += models[j]->ncoeff;
                    continue;
                }
                err = make_sys_X_block(Xj, models[j], dset, sys->t1, method);
                Xk = Xj;
            } else if (method == SYS_METHOD_LIML) {
                err = make_liml_X_block(Xj, models[i], dset, sys->t1);
                Xk = Xj;
            } else {
                Xk = Xi;
            }

            M->rows = Xi->cols;
            M->cols = Xk->cols;

            err = gretl_matrix_multiply_mod(Xi, GRETL_MOD_TRANSPOSE,
                                            Xk, GRETL_MOD_NONE,
                                            M, GRETL_MOD_NONE);

            if (rsingle || (single_equation && method != SYS_METHOD_WLS)) {
                sij = 1.0;
            } else {
                sij = gretl_matrix_get(sys->S, i, j);
            }

            insert_sys_X_block(X, M, krow, kcol, sij);
            kcol += models[j]->ncoeff;
        }

        krow += models[i]->ncoeff;
    }

    if (err) {
        fprintf(stderr, "after trying to make X matrix: err = %d\n", err);
        goto cleanup;
    }

    if (nr > 0) {
        /* there are restrictions to be imposed */
        augment_X_with_restrictions(X, mk, sys);
    }

    if (!do_iteration && !rsingle) {
        /* we're not coming back this way, so free some storage */
        gretl_matrix_free(Xj);
        Xj = NULL;
        gretl_matrix_free(M);
        M = NULL;
    }

    /* form stacked Y column vector (m x k) */

    v = 0;
    for (i=0; i<sys->neqns; i++) {

        /* loop over the m vertically-arranged blocks */
        make_sys_X_block(Xi, models[i], dset, sys->t1, method);

        for (j=0; j<models[i]->ncoeff; j++) {
            /* loop over the rows within each of the m blocks */
            double yv = 0.0;
            int lmin = 0, lmax = sys->neqns;

            if (single_equation || rsingle) {
                /* no cross terms wanted */
                lmin = i;
                lmax = i + 1;
            }

            for (l=lmin; l<lmax; l++) {
                /* loop over the components that must be
                   added to form each element */
                const double *yl = NULL;
                double xx = 0.0;

                if (method == SYS_METHOD_LIML) {
                    yl = gretl_model_get_data(models[l], "liml_y");
                } else {
                    yl = dset->Z[system_get_depvar(sys, l)];
                }

                /* multiply X'[l] into y */
                for (t=0; t<T; t++) {
                    xx += gretl_matrix_get(Xi, t, j) * yl[t + sys->t1];
                }

                if (rsingle || (single_equation && method != SYS_METHOD_WLS)) {
                    ; /* leave xx unmodified */
                } else if (method == SYS_METHOD_WLS && i == l) {
                    xx *= gretl_matrix_get(sys->S, i, i);
                } else {
                    xx *= gretl_matrix_get(sys->S, i, l);
                }

                yv += xx;
            }
            gretl_vector_set(y, v++, yv);
        }
    }

    if (nr > 0) {
        /* there are restrictions */
        augment_y_with_restrictions(y, mk, nr, sys);
    }

    /* The estimates calculated below will be SUR, 3SLS or LIML,
       depending on how the data matrices above were constructed --
       unless, that is, we're just doing restricted OLS, WLS or TSLS
       estimates.
    */
    err = calculate_sys_coeffs(sys, dset, X, y, mk, nr,
                               do_iteration);

    if (!err && rsingle) {
        /* take one more pass */
        rsingle = 0;
        goto iteration_start;
    }

    if (!err && do_iteration) {
        /* check for convergence */
        if (!sys_converged(sys, &llbak, opt, prn, &err)) {
            if (!err) {
                sys->iters += 1;
                goto iteration_start;
            }
        }
    }

    if (err) goto cleanup;

    if (nr == 0 && (method == SYS_METHOD_3SLS || method == SYS_METHOD_SUR)) {
        /* compute this test while we have sigma-inverse available */
        hansen_sargan_test(sys, dset);
    }

    do_diag = !(method == SYS_METHOD_SUR && do_iteration);

    /* refresh sigma (non-inverted) */
    gls_sigma_from_uhat(sys, sys->S, do_diag);

    if (method == SYS_METHOD_FIML) {
        /* compute FIML estimates */
        err = fiml_driver(sys, dset, opt, prn);
    }

 save_etc:

    if (!err && !(sys->flags & SYSTEM_LIML1)) {
        err = system_save_and_print_results(sys, dset, opt, prn);
    }

 cleanup:

    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);
    gretl_matrix_free(M);
    gretl_matrix_free(X);
    gretl_matrix_free(y);

    clean_up_models(sys);

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    return err;
}
