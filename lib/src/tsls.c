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

/* tsls.c - Two-Stage Least Squares for gretl */

#include "libgretl.h"
#include "qr_estimate.h"
#include "libset.h"
#include "estim_private.h"
#include "missing_private.h"
#include "matrix_extra.h"
#include "tsls.h"

#define TDEBUG 0

struct iv_info_ {
    MODEL *pmod;      /* pointer to the model to be returned */
    gretl_matrix *Q;  /* matrix for creating "hat" variables */
    gretl_matrix *g;  /* auxiliary matrix */
    char *missmask;   /* missing observations mask */
    int *reglist;     /* full list of regressors */
    int *instlist;    /* list of instruments */
    int *hatlist;     /* list of fitted values of endogenous regressors */
    int *s2list;      /* regression list for second stage */
    int *exolist;     /* list of exogenous regressors */
    int *endolist;    /* list of endogenous regressors */
    int *idroplist;   /* list of dropped instruments */
    int overid;       /* over-identification rank */
};

typedef struct iv_info_ iv_info;

static iv_info *iv_info_new (MODEL *pmod)
{
    iv_info *ivi = malloc(sizeof *ivi);

    if (ivi != NULL) {
	ivi->pmod = pmod;
	ivi->Q = ivi->g = NULL;
	ivi->missmask = NULL;
	ivi->reglist = NULL;
	ivi->instlist = NULL;
	ivi->hatlist = NULL;
	ivi->s2list = NULL;
	ivi->exolist = NULL;
	ivi->endolist = NULL;
	ivi->idroplist = NULL;
	ivi->overid = 0;
    }

    return ivi;
}

static void iv_info_destroy (iv_info *ivi)
{
    if (ivi == NULL) {
	return;
    }

    gretl_matrix_free(ivi->Q);
    gretl_matrix_free(ivi->g);

    free(ivi->missmask);
    free(ivi->reglist);
    free(ivi->instlist);
    free(ivi->hatlist);
    free(ivi->s2list);
    free(ivi->exolist);
    free(ivi->endolist);
    free(ivi->idroplist);

    free(ivi);
}

static void tsls_omitzero (int *list, const char *mask,
			   const DATASET *dset)
{
    int i, v, t, allzero;

    for (i=list[0]; i>1; i--) {
        v = list[i];
        allzero = 1;
        for (t=dset->t1; t<=dset->t2; t++) {
            if (mask != NULL && mask[t - dset->t1]) {
                continue;
            }
            if (dset->Z[v][t] != 0.0) {
                allzero = 0;
                break;
            }
        }
        if (allzero) {
            gretl_list_delete_at_pos(list, i);
        }
    }
}

/**
 * tsls_free_data:
 * @pmod: model to operate on.
 *
 * Frees the first-stage data attached to a two-stage least squares
 * model that was estimated using %OPT_E (part of a system).
 */

void tsls_free_data (const MODEL *pmod)
{
    gretl_matrix *endog = gretl_model_get_data(pmod, "endog");
    double **X = gretl_model_get_data(pmod, "tslsX");
    int i, m = 0;

    if (endog != NULL && X != NULL) {
        for (i=0; i<pmod->ncoeff; i++) {
            if (endog->val[i] != 0) m++;
        }
        for (i=0; i<m; i++) {
            free(X[i]);
        }
    }
}

/* when tsls is called for initial estimation of an equation
   within a system of equations, save the first-stage fitted
   values on the model under the key "tslsX".
*/

static int tsls_save_data (iv_info *ivi, DATASET *dset)
{
    MODEL *pmod = ivi->pmod;
    gretl_matrix *endog = NULL;
    double **X = NULL;
    size_t Xsize, xs_old;
    int recycle_X = 0;
    int recycle_e = 0;
    int i, v, err = 0;

    Xsize = ivi->hatlist[0] * sizeof *X;

    /* re-use old pointers if applicable */

    X = gretl_model_get_data_full(pmod, "tslsX", NULL, NULL, &xs_old);
    if (X != NULL) {
        if (Xsize == xs_old) {
            recycle_X = 1;
        } else {
            tsls_free_data(pmod);
            gretl_model_detach_data_item(pmod, "tslsX");
            free(X);
            X = NULL;
        }
    }

    if (!recycle_X && Xsize > 0) {
        X = malloc(Xsize);
        if (X == NULL) {
            return E_ALLOC;
        }
    }

    endog = gretl_model_get_data(pmod, "endog");
    if (endog != NULL) {
        if (gretl_vector_get_length(endog) == pmod->ncoeff) {
            recycle_e = 1;
        } else {
            gretl_model_detach_data_item(pmod, "endog");
            gretl_matrix_free(endog);
            endog = NULL;
        }
    }

    if (!recycle_e) {
        endog = gretl_matrix_alloc(pmod->ncoeff, 1);
        if (endog == NULL) {
            free(X);
            return E_ALLOC;
        }
    }

    /* steal the appropriate columns from Z */
    for (i=1; i<=ivi->hatlist[0]; i++) {
        v = ivi->hatlist[i];
        X[i-1] = dset->Z[v];
        dset->Z[v] = NULL;
    }

    for (i=0; i<pmod->ncoeff; i++) {
        v = pmod->list[i+2];
        endog->val[i] = !in_gretl_list(ivi->exolist, v);
    }

     /* now attach X and endog to the model */
    if (!recycle_X && X != NULL) {
        gretl_model_set_data(pmod, "tslsX", X, GRETL_TYPE_DOUBLE_ARRAY,
                             Xsize);
    }
    if (!recycle_e && endog != NULL) {
        gretl_model_set_matrix_as_data(pmod, "endog", endog);
    }

    return err;
}

static void add_tsls_var (int *list, int v, gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i;

    if (opt & OPT_T) {
        /* add as instrument */
        list[0] += 1;
        list[list[0]] = v;
    } else if (opt & OPT_B) {
        /* add as exogenous regressor */
        list[0] += 2;
        for (i=list[0]-1; i>pos; i--) {
            list[i] = list[i-1];
        }
        list[pos] = v;
        list[list[0]] = v;
    } else {
        /* add as endogenous regressor */
        list[0] += 1;
        for (i=list[0]; i>pos; i--) {
            list[i] = list[i-1];
        }
        list[pos] = v;
    }
}

static void delete_tsls_var (int *list, int v, gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i, j;

    if (opt & OPT_T) {
        /* delete as instrument */
        for (i=pos+1; i<=list[0]; i++) {
            if (list[i] == v) {
                for (j=i; j<list[0]; j++) {
                    list[j] = list[j + 1];
                }
                list[0] -= 1;
                break;
            }
        }
    } else if (opt & OPT_B) {
        /* delete from both sub-lists */
        for (i=2; i<=list[0]; i++) {
            if (list[i] == v) {
                for (j=i; j<list[0]; j++) {
                    list[j] = list[j + 1];
                }
                list[0] -= 1;
            }
        }
    } else {
        /* delete from regressor list */
        for (i=2; i<pos; i++) {
            if (list[i] == v) {
                for (j=i; j<list[0]; j++) {
                    list[j] = list[j + 1];
                }
                list[0] -= 1;
                break;
            }
        }
    }
}

static int in_ivreg_list (const int *list, int v, gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i, imin, imax;

    if (pos == 0) {
        return 0;
    }

    if (opt & OPT_T) {
        imin = pos + 1;
        imax = list[0];
    } else if (opt & OPT_B) {
        imin = 2;
        imax = list[0];
    } else {
        imin = 2;
        imax = pos - 1;
    }

    for (i=imin; i<=imax; i++) {
        if (list[i] == v) {
            return i;
        }
    }

    return 0;
}

/**
 * ivreg_list_omit:
 * @orig: list specifying original IV model.
 * @drop: list of variables to be omitted.
 * @opt: may include %OPT_T (omit variable only as instrument)
 * or %OPT_B (omit both as regressor and instrument).  The default
 * is to omit (only) as regressor.
 * @err: pointer to receive error code.
 *
 * Creates a new IVREG specification list, by first copying @orig
 * then deleting from the copy the variables found in @drop.
 *
 * Returns: the new, reduced list or %NULL on error.
 */

int *
ivreg_list_omit (const int *orig, const int *drop, gretlopt opt, int *err)
{
    int *newlist;
    int i;

    *err = incompatible_options(opt, OPT_T | OPT_B);
    if (*err) {
        return NULL;
    }

    newlist = gretl_list_copy(orig);

    for (i=1; i<=drop[0]; i++) {
        if (in_ivreg_list(orig, drop[i], opt)) {
            delete_tsls_var(newlist, drop[i], opt);
        } else {
            *err = E_UNSPEC;
        }
    }

    if (*err) {
        free(newlist);
        newlist = NULL;
    }

    return newlist;
}

/**
 * ivreg_list_add:
 * @orig: list specifying original IV model.
 * @add: list of variables to be added.
 * @opt: may include %OPT_T (add variable only as instrument)
 * or %OPT_B (add both as regressor and instrument).  The default
 * is to add as an endogenous regressor.
 * @err: pointer to receive error code.
 *
 * Creates a new IVREG specification list, by copying @orig
 * and adding the variables found in @add.
 *
 * Returns: the new, augmented list or %NULL on error.
 */

int *
ivreg_list_add (const int *orig, const int *add, gretlopt opt, int *err)
{
    int norig = orig[0];
    int nadd = add[0];
    int *newlist;
    int i;

    *err = incompatible_options(opt, OPT_T | OPT_B);
    if (*err) {
        return NULL;
    }

    if (opt & OPT_B) {
        nadd *= 2;
    }

    newlist = gretl_list_new(norig + nadd);

    for (i=0; i<=norig; i++) {
        newlist[i] = orig[i];
    }

    for (i=1; i<=add[0]; i++) {
        if (in_ivreg_list(orig, add[i], opt)) {
            *err = E_ADDDUP;
        } else {
            add_tsls_var(newlist, add[i], opt);
        }
    }

    if (*err) {
        free(newlist);
        newlist = NULL;
    }

    return newlist;
}

/* prepend constant to list of instruments */

static int zlist_prepend_const (int **pzlist)
{
    int clist[2] = {1, 0};

    return gretl_list_insert_list(pzlist, clist, 1);
}

/**
 * tsls_make_endolist:
 * @reglist: regression specification.
 * @instlist: predetermined and exogenous variables.
 * @err: location to receive error code.
 *
 * Determines which variables in @reglist, when checked against the
 * predetermined and exogenous vars in @instlist, need to be
 * instrumented, and populates the returned list accordingly.
 *
 * Returns: allocated list of variables to be instrumented or
 * NULL if there are no such variables.
 */

int *tsls_make_endolist (const int *reglist,
			 const int *instlist,
                         int *err)
{
    int *ret = NULL;
    int i, nendo = 0;

    for (i=2; i<=reglist[0]; i++) {
        if (!in_gretl_list(instlist, reglist[i])) {
	    nendo++;
	}
    }

    if (nendo == 0) {
	return ret;
    }

    ret = gretl_list_new(nendo);
    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	int vi, j = 1;

	for (i=2; i<=reglist[0]; i++) {
	    vi = reglist[i];
	    if (!in_gretl_list(instlist, vi)) {
		ret[j++] = vi;
	    }
	}
    }

    return ret;
}

/* Fill the residuals matrix for tsls likelihood calculation.
   We run a multivariate regression of the dependent variable
   and the endogenous regressors on all the instruments.
*/

static int fill_E_matrix (gretl_matrix *E, iv_info *ivi,
                          const DATASET *dset)
{
    MODEL *pmod = ivi->pmod;
    gretl_matrix *Y = NULL;
    gretl_matrix *X = NULL;
    int ny, *ylist;
    int i, err = 0;

    /* form list of LHS series */
    ny = ivi->endolist[0] + 1;
    ylist = gretl_list_new(ny);
    ylist[1] = ivi->reglist[1];
    for (i=2; i<=ny; i++) {
        ylist[i] = ivi->endolist[i-1];
    }

    /* LHS: depvar and endogenous regressors */
    Y = gretl_matrix_data_subset_masked(ylist, dset, pmod->t1, pmod->t2,
                                        ivi->missmask, &err);
    if (!err) {
        /* RHS: all instruments */
        X = gretl_matrix_data_subset_masked(ivi->instlist, dset,
					    pmod->t1, pmod->t2,
                                            ivi->missmask, &err);
    }

    if (!err) {
        /* populate the residuals matrix, E */
        err = gretl_matrix_multi_ols(Y, X, NULL, E, NULL);
    }

    gretl_matrix_free(Y);
    gretl_matrix_free(X);
    free(ylist);

    return err;
}

/* calculate the log determinant */

static double tsls_get_ldet (gretl_matrix *S, int *err)
{
    double ret = NADBL;
    char *mask;

    /* allow for the possibility that S is rank-deficient? */
    mask = gretl_matrix_rank_mask(S, err);
    if (mask != NULL) {
        *err = gretl_matrix_cut_rows_cols(S, mask);
    }
    if (!*err) {
        ret = gretl_matrix_log_determinant(S, err);
    }

    return ret;
}

/* calculate the maximized log-likelihood for a just-identified
   model */

static int tsls_loglik (iv_info *ivi, DATASET *dset)
{
    gretl_matrix *E, *S;
    int T = ivi->pmod->nobs;
    int k = 1 + ivi->endolist[0];
    int err = 0;

    ivi->pmod->lnL = NADBL;

    E = gretl_matrix_alloc(T, k);
    S = gretl_matrix_alloc(k, k);

    if (E == NULL || S == NULL) {
        err = E_ALLOC;
    } else {
        err = fill_E_matrix(E, ivi, dset);
    }

    if (!err) {
        err = gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
                                        E, GRETL_MOD_NONE,
                                        S, GRETL_MOD_NONE);
    }

    if (!err) {
        double ldet = tsls_get_ldet(S, &err);

        if (na(ldet)) {
            ivi->pmod->lnL = NADBL;
        } else {
            /* Davidson and MacKinnon, ETM, p. 538, taking
               kappa = 1 for the just-identified case
            */
            ivi->pmod->lnL = -(T / 2.0) * (LN_2_PI + ldet);
        }
    }

    mle_criteria(ivi->pmod, 0);

    gretl_matrix_free(E);
    gretl_matrix_free(S);

    return err;
}

/* perform the Sargan overidentification test for an IV model */

static int ivreg_sargan_test (iv_info *ivi, DATASET *dset)
{
    MODEL *pmod = ivi->pmod;
    MODEL smod;
    int i, t, err = 0;
    int *OT_list = NULL;
    int nv = dset->v;

    if (pmod->list[0] == 2 && pmod->list[2] == 0) {
        /* degenerate model with const only */
        return 0;
    }

    err = dataset_add_series(dset, 1);
    if (err) {
        return err;
    }

    /* add the @pmod residual to the dataset */
    for (t=pmod->t1; t<=pmod->t2; t++) {
        dset->Z[nv][t] = pmod->uhat[t];
    }

    OT_list = gretl_list_new(ivi->instlist[0] + 1);
    if (OT_list == NULL) {
        dataset_drop_last_variables(dset, 1);
        return E_ALLOC;
    }

    OT_list[1] = nv;
    for (i=2; i<=OT_list[0]; i++) {
        OT_list[i] = ivi->instlist[i-1];
    }

    /* OLS: the dependent variable is the residual from @pmod;
       the regressors are the (full set of) instruments
    */
    smod = lsq(OT_list, dset, OLS, OPT_A);

    if (smod.errcode) {
        fprintf(stderr, "ivreg_sargan_test: smod.errcode = %d\n", smod.errcode);
        err = smod.errcode;
    } else {
        ModelTest *test = model_test_new(GRETL_TEST_SARGAN);
        double OTest = smod.rsq * smod.nobs;

        if (test != NULL) {
            model_test_set_teststat(test, GRETL_STAT_LM);
            model_test_set_dfn(test, ivi->overid);
            model_test_set_value(test, OTest);
            model_test_set_pvalue(test, chisq_cdf_comp(ivi->overid, OTest));
            maybe_add_test_to_model(pmod, test);
        }
    }

    clear_model(&smod);
    dataset_drop_last_variables(dset, 1);
    free(OT_list);

    return err;
}

static int hausman_drop_check (MODEL *hmod,
                               int *hatlist,
                               int *endolist,
                               DATASET *dset)
{
    int *dlist = gretl_model_get_list(hmod, "droplist");
    int n_drop = 0;

    if (dlist != NULL) {
        int i, pos;

        for (i=1; i<=dlist[0]; i++) {
            pos = in_gretl_list(hatlist, dlist[i]);
            if (pos > 0 && pos <= endolist[0]) {
                n_drop++;
                gretl_warnmsg_sprintf("%s does not seem to be endogenous",
                                      dset->varname[endolist[pos]]);
            }
        }
    }

    return n_drop;
}

static double tsls_robust_hausman_fixup (const int *hatlist,
                                         MODEL *pmod,
                                         int *df)
{
    int *wlist = gretl_list_copy(hatlist);
    double H = NADBL;

    if (wlist != NULL) {
        int i;

        for (i=wlist[0]; i>0; i--) {
            if (!in_gretl_list(pmod->list, wlist[i])) {
                gretl_list_delete_at_pos(wlist, i);
            }
        }
        H = wald_omit_chisq(wlist, pmod);
        *df = wlist[0];
        free(wlist);
    }

    return H;
}

static void add_hausman_to_model (MODEL *pmod, double HTest, int df)
{
    ModelTest *test = model_test_new(GRETL_TEST_IV_HAUSMAN);

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
        model_test_set_dfn(test, df);
        model_test_set_value(test, HTest);
        model_test_set_pvalue(test, chisq_cdf_comp(df, HTest));
        maybe_add_test_to_model(pmod, test);
    }
}

/* Hausman test via the matrix method (both plain and robust).  See
   Alecos Papadopoulos, 'A New Matrix Statistic for the Hausman
   Endogeneity Test under Heteroskedasticity', Econometrics 11:23,
   2023. https://doi.org/10.3390/econometrics11040023
*/

static int tsls_matrix_hausman_test (iv_info *ivi,
				     gretlopt opt,
				     DATASET *dset)
{
    const char *mask = ivi->missmask;
    MODEL *pmod = ivi->pmod;
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *X1 = NULL;
    gretl_matrix *Z = NULL;
    gretl_matrix *u = NULL;
    gretl_matrix *E = NULL;
    gretl_matrix *X1u = NULL;
    gretl_matrix *A = NULL;
    int *xlist = NULL;
    int i, ylist[2];
    double q;
    int err = 0;

    /* Step 1: obtain residuals from estimation of the specified
       model via plain OLS.
    */

    xlist = gretl_list_new(ivi->reglist[0] - 1);
    for (i=1; i<=xlist[0]; i++) {
	xlist[i] = ivi->reglist[i+1];
    }
    ylist[0] = 1;
    ylist[1] = ivi->reglist[1];

    /* all regressors */
    X = gretl_matrix_data_subset_masked(xlist, dset, pmod->t1, pmod->t2,
					mask, &err);
    if (!err) {
	y = gretl_matrix_data_subset_masked(ylist, dset, pmod->t1, pmod->t2,
					    mask, &err);
    }
    if (!err) {
	/* put OLS residuals into @u */
	u = gretl_matrix_alloc(X->rows, 1);
	gretl_matrix_ols(y, X, NULL, NULL, u, NULL);
    }

    /* Step 2: fill @X1 with \hat{X1} */
    X1 = gretl_matrix_data_subset_masked(ivi->hatlist, dset, pmod->t1, pmod->t2,
					 mask, &err);

    if (!err) {
	/* some more required allocation */
	X1u = gretl_matrix_alloc(X1->cols, 1);
	E = gretl_matrix_alloc(X1->rows, X1->cols);
	A = gretl_matrix_alloc(X1->cols, X1->cols);
	if (X1u == NULL || E == NULL || A == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	goto bailout;
    }

    /* Step 3: form X1u = \hat{X}_1'\hat{u} */
    gretl_matrix_multiply_mod(X1, GRETL_MOD_TRANSPOSE,
			      u, GRETL_MOD_NONE,
			      X1u, GRETL_MOD_NONE);

    /* Step 4: form E = M_x \hat{X}_1 */
    err = gretl_matrix_multi_ols(X1, X, NULL, E, NULL);

    if (!err) {
	if (opt & OPT_R) {
	    /* Step 5 (robust): A = qform(M_x \hat{X_1}', diag(u.^2)) */
	    double ui;

	    for (i=0; i<u->rows; i++) {
		ui = u->val[i];
		u->val[i] = ui * ui;
	    }
	    gretl_matrix_diag_qform(E, GRETL_MOD_TRANSPOSE, u,
				    A, GRETL_MOD_NONE);
	} else {
	    /* Step 5 (plain): A = \hat{X_1}' M_x \hat{X_1} */
	    double s2 = 0;

	    for (i=0; i<u->rows; i++) {
		s2 += u->val[i] * u->val[i];
	    }
	    s2 /= u->rows;
	    gretl_matrix_multiply_mod(X1, GRETL_MOD_TRANSPOSE,
				      E, GRETL_MOD_NONE,
				      A, GRETL_MOD_NONE);
	    gretl_matrix_multiply_by_scalar(A, s2);
	}
	err = gretl_invert_symmetric_matrix(A);
	if (err) {
	    gretl_errmsg_set("Hausman test matrix could not be inverted: "
			     "this indicates invalid instruments");
	}
    }

    if (err) {
	q = NADBL;
    } else {
	q = gretl_scalar_qform(X1u, A, &err);
	if (!err) {
	    add_hausman_to_model(pmod, q, X1->cols);
	}
    }

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(X1);
    gretl_matrix_free(Z);
    gretl_matrix_free(u);
    gretl_matrix_free(E);
    gretl_matrix_free(X1u);
    gretl_matrix_free(A);
    free(xlist);

    return err;
}

/* Hausman test via the regression method */

static int tsls_hausman_test (iv_info *ivi, gretlopt opt, DATASET *dset)
{
    MODEL hmod;
    double URSS;
    gretlopt uopt;
    int *HT_list = NULL;
    int t, df;
    int nv = dset->v;
    int n_drop = 0;
    int ku = 0;
    int err = 0;

    /* create regression list for unrestricted model by appending
       @hatlist to @reglist */
    HT_list = gretl_list_add(ivi->reglist, ivi->hatlist, &err);
    if (err == E_NOADD) {
        err = 0; /* more of a no-op than an error */
    } else if (err) {
        return err;
    }

    /* option for unrestricted regression */
    uopt = (opt & OPT_R)? (OPT_A | OPT_R) : OPT_A;

    /* estimate the unrestricted model */
    hmod = lsq(HT_list, dset, OLS, uopt);
    if (hmod.errcode) {
        fprintf(stderr, "tsls_hausman_test: hmod.errcode (U) = %d\n", hmod.errcode);
        err = hmod.errcode;
        goto bailout;
    }

    n_drop = hausman_drop_check(&hmod, ivi->hatlist, ivi->endolist, dset);

    if (hmod.dfd == 0) {
        /* perfect fit, can't do test */
        goto bailout;
    }

    if (opt & OPT_R) {
        /* 2023-08-31: if the --robust option was in force for
           tsls estimation, compute the Hausman test as a Wald
           test on the robust variance matrix for the U model.
        */
        double HTest = NADBL;

        if (n_drop > 0) {
            HTest = tsls_robust_hausman_fixup(ivi->hatlist, &hmod, &df);
        } else {
            HTest = wald_omit_chisq(ivi->hatlist, &hmod);
            df = ivi->hatlist[0];
        }
        if (!na(HTest)) {
            add_hausman_to_model(ivi->pmod, HTest, df);
        }
        goto bailout;
    }

    /* record U-model info */
    URSS = hmod.ess;
    ku = hmod.ncoeff;

    /* add fitted values from unrestricted model to dataset */
    err = dataset_add_series(dset, 1);
    if (err) {
        goto bailout;
    } else {
        for (t=hmod.t1; t<=hmod.t2; t++) {
            dset->Z[nv][t] = hmod.yhat[t];
        }
    }

    clear_model(&hmod);

    /* redefine HT_list: dep var is the fitted value from the
       U-model, regressors are from @reglist */
    free(HT_list);
    HT_list = gretl_list_copy(ivi->reglist);
    HT_list[1] = nv;

#if TDEBUG
    fprintf(stderr, "running regression for Hausman test, 2\n");
#endif

    /* regress the U fitted values on the regressors of the restricted
       model (so the sum of squares equals RRSS-URSS)
    */
    hmod = lsq(HT_list, dset, OLS, OPT_A);

    if (hmod.errcode) {
        fprintf(stderr, "tsls_hausman_test: hmod.errcode (D) = %d\n", hmod.errcode);
        err = hmod.errcode;
    } else if ((df = ku - hmod.ncoeff) > 0) {
        /* Degrees of freedom for the Hausman test need not be equal to the
           number of added regressors, since some of them may not really be
           endogenous.
        */
        double DRSS = hmod.ess;
        double HTest = (DRSS/URSS) * hmod.nobs;

        add_hausman_to_model(ivi->pmod, HTest, df);
    }

 bailout:

    dataset_drop_last_variables(dset, dset->v - nv);
    clear_model(&hmod);
    free(HT_list);

    return err;
}

/* form matrix of instruments and perform QR decomposition
   of this matrix; attach return Q to @ivi; return 0 on
   success, error code on error.
*/

static int tsls_add_Q (iv_info *ivi, const DATASET *dset)
{
    gretl_matrix *R = NULL;
    int rank, ndrop = 0;
    int *dlist = NULL;
    double test;
    int i, j, k;
    int err = 0;

    ivi->Q = gretl_matrix_data_subset_masked(ivi->instlist, dset,
					     dset->t1, dset->t2,
					     ivi->missmask, &err);
    if (err) {
        return err;
    }

    k = gretl_matrix_cols(ivi->Q);

    if (k > ivi->Q->rows) {
        /* can't do QR decomp! */
        gretl_errmsg_sprintf(_("Number of instruments (%d) exceeds the "
                               "number of observations (%d)"),
                             k, ivi->Q->rows);
        err = E_DF;
        goto bailout;
    }

    R = gretl_matrix_alloc(k, k);
    if (R == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    err = gretl_matrix_QR_decomp(ivi->Q, R);
    if (err) {
        goto bailout;
    }

    rank = gretl_check_QR_rank(R, &err, NULL);

#if TDEBUG
    fprintf(stderr, "tsls_Q: k = Q->cols = %d, rank = %d\n", k, rank);
    if (rank == 0) {
        gretl_matrix_print(R, "R");
    }
#endif

    if (err) {
        goto bailout;
    } else if (rank < k) {
        fprintf(stderr, "tsls_Q: k = %d, rank = %d\n", k, rank);
        ndrop = k - rank;
    }

    if (ndrop > 0) {
        dlist = gretl_list_new(ndrop);
        if (dlist != NULL) {
            dlist[0] = 0;
        }
        j = 1;
        for (i=0; i<k; i++) {
            test = gretl_matrix_get(R, i, i);
            if (fabs(test) < R_DIAG_MIN) {
                if (dlist != NULL) {
                    dlist[0] += 1;
                    dlist[dlist[0]] = ivi->instlist[j];
                }
                fprintf(stderr, "tsls_Q: dropping redundant instrument %d (%s)\n",
                        ivi->instlist[j], dset->varname[ivi->instlist[j]]);
                gretl_list_delete_at_pos(ivi->instlist, j--);
            }
            j++;
        }

        k = ivi->instlist[0];
        gretl_matrix_free(ivi->Q);
        ivi->Q = gretl_matrix_data_subset_masked(ivi->instlist, dset,
						 dset->t1, dset->t2,
						 ivi->missmask, &err);
        if (!err) {
            R = gretl_matrix_reuse(R, k, k);
            err = gretl_matrix_QR_decomp(ivi->Q, R);
        }
    }

 bailout:

    gretl_matrix_free(R);

    if (err) {
        free(dlist);
        gretl_matrix_free(ivi->Q);
        ivi->Q = NULL;
    } else {
        ivi->idroplist = dlist;
    }

    return err;
}

static int tsls_form_xhat (iv_info *ivi, DATASET *dset,
			   int v0, int v1)
{
    const char *mask = ivi->missmask;
    const double *x = dset->Z[v0];
    double *xhat = dset->Z[v1];
    int k = gretl_matrix_cols(ivi->Q);
    int allzero = 1;
    int i, t, s = 0;
    int err = 0;

#if TDEBUG > 1
    fprintf(stderr, "tsls_form_xhat: v0=%d, v1=%d, t1=%d, t2=%d, mask=%p\n",
            v0, v1, dset->t1, dset->t2, (void *) mask);
#endif

    /* form g = Q'x */
    for (i=0; i<k; i++) {
        ivi->g->val[i] = 0.0;
        s = 0;
        for (t=dset->t1; t<=dset->t2; t++) {
            if (mask != NULL && mask[t - dset->t1]) {
                continue;
            }
            ivi->g->val[i] += gretl_matrix_get(ivi->Q, s++, i) * x[t];
        }
    }

    /* form \hat{x} = Qg = QQ'x */
    s = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
        if (mask != NULL && mask[t - dset->t1]) {
            xhat[t] = NADBL;
            continue;
        }
        xhat[t] = 0.0;
        for (i=0; i<k; i++) {
            xhat[t] += gretl_matrix_get(ivi->Q, s, i) * ivi->g->val[i];
        }
        if (xhat[t] != 0) {
            allzero = 0;
        }
        s++;
    }

    if (allzero) {
        gretl_errmsg_sprintf(_("The first-stage fitted values for %s are all zero"),
                             dset->varname[v0]);
        err = E_DATA;
    } else {
        /* name the fitted series according to the original */
        strcpy(dset->varname[v1], "h_");
        strncat(dset->varname[v1], dset->varname[v0], VNAMELEN - 3);
    }

    return err;
}

static void tsls_residuals (iv_info *ivi, gretlopt opt,
                            const DATASET *dset)
{
    MODEL *pmod = ivi->pmod;
    double yh, sigma0 = pmod->sigma;
    int yno = ivi->reglist[1];
    int i, t;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
        if (model_missing(pmod, t)) {
            continue;
        }
        yh = 0.0;
        for (i=0; i<pmod->ncoeff; i++) {
            yh += pmod->coeff[i] * dset->Z[ivi->reglist[i+2]][t];
        }
        pmod->yhat[t] = yh;
        pmod->uhat[t] = dset->Z[yno][t] - yh;
        pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    if (pmod->ess <= 0.0) {
        pmod->sigma = 0.0;
    } else {
        int den = (opt & OPT_N)? pmod->nobs : pmod->dfd;

        pmod->sigma = sqrt(pmod->ess / den);
    }

    if (sigma0 > 0.0) {
        double corr = pmod->sigma / sigma0;

        for (i=0; i<pmod->ncoeff; i++) {
            pmod->sderr[i] *= corr;
        }
    }
}

static void tsls_extra_stats (iv_info *ivi, gretlopt opt,
			      const DATASET *dset)
{
    MODEL *pmod = ivi->pmod;
    int yno = ivi->reglist[1];
    double r;

    /* squared correlation between y and yhat */
    pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, dset->Z[yno], pmod->yhat);

    r = 1 - pmod->rsq;
    pmod->adjrsq = 1.0 - (r * (pmod->nobs - 1) / pmod->dfd);

    pmod->fstt = pmod->chisq = NADBL;

    if (pmod->ncoeff > 1) {
        if ((opt & OPT_N) || (ivi->overid && pmod->ncoeff == 2)) {
            pmod->chisq = wald_omit_chisq(NULL, pmod);
        } else {
            pmod->fstt = wald_omit_F(NULL, pmod);
        }
    }

    /* tsls is not an ML estimator unless it's exactly identified,
       and in that case we require a special calculation */
    pmod->lnL = NADBL;
    pmod->criterion[C_AIC] = NADBL;
    pmod->criterion[C_BIC] = NADBL;
    pmod->criterion[C_HQC] = NADBL;

    if (dataset_is_time_series(dset) && !model_has_missing_obs(pmod)) {
        /* time series, no missing obs within sample range */
        pmod->rho = rhohat(1, pmod->t1, pmod->t2, pmod->uhat);
        pmod->dw = dwstat(1, pmod, dset);
    } else {
        pmod->rho = pmod->dw = NADBL;
    }
}

static void tsls_recreate_full_list (iv_info *ivi)
{
    int *full_list;

    if (ivi->instlist != NULL && ivi->instlist[0] > 0) {
        full_list = gretl_lists_join_with_separator(ivi->reglist, ivi->instlist);
    } else {
        full_list = gretl_list_copy(ivi->reglist);
    }

    if (full_list == NULL) {
        ivi->pmod->errcode = E_ALLOC;
    } else {
        free(ivi->pmod->list);
        ivi->pmod->list = full_list;
    }
}

static void replace_list_element (int *list, int targ, int repl)
{
    int i;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == targ) {
            list[i] = repl;
            break;
        }
    }
}

static int reglist_remove_redundant_vars (iv_info *ivi)
{
    int *dlist = gretl_model_get_list(ivi->pmod, "droplist");
    int i, pos;

    if (dlist == NULL) {
        return 1;
    }

    for (i=1; i<=dlist[0]; i++) {
        pos = in_gretl_list(ivi->s2list, dlist[i]);
        if (pos > 1) {
            gretl_list_delete_at_pos(ivi->reglist, pos);
        }
        if (ivi->endolist != NULL) {
            /* note: if endolist is NULL, so is hatlist */
            pos = in_gretl_list(ivi->hatlist, dlist[i]);
            if (pos > 1) {
                /* First replace dlist[i] with the ID of the original
                   regressor from @endolist, in place of the "hatlist"
                   variable (which will be deleted when the model is
                   returned), so that the printout of regressors that are
                   dropped due to exact collinearity will show the
                   appropriate names.

                   Second, delete the redundant term from both endolist
                   and hatlist, so that subsequent calculations will not
                   get messed up.
                */
                dlist[i] = ivi->endolist[pos];
                gretl_list_delete_at_pos(ivi->endolist, pos);
                gretl_list_delete_at_pos(ivi->hatlist, pos);
            }
        }
    }

    return 0;
}

/* Stock-Yogo: compute test statistic for weakness of instruments.
   This is the Cragg-Donald statistic: the minimum eigenvalue of the
   matrix defined below (which reduces to the first-stage F-statistic
   when there is just one endogenous regressor).  See Stock and Yogo,
   "Testing for Weak Instruments in Linear IV Regressions" (originally
   NBER Technical Working Paper 284; revised February 2003), at
   http://ksghome.harvard.edu/~JStock/pdf/rfa_6.pdf
*/

static int compute_stock_yogo (iv_info *ivi, const DATASET *dset)
{
    const double **dZ = (const double **) dset->Z;
    MODEL *pmod = ivi->pmod;
    gretl_matrix_block *B = NULL;
    gretl_matrix_block *B2 = NULL;
    gretl_matrix *G, *S, *Y;
    gretl_matrix *Ya, *YPY;
    gretl_matrix *X, *Z, *E = NULL;
    int T = pmod->nobs;
    int n = ivi->endolist[0];
    int K1 = pmod->ncoeff - n;
    int K2 = ivi->instlist[0] - K1;
    int i, s, vi, vj, t;
    int err = 0;

    B = gretl_matrix_block_new(&G, n, n,
                               &S, n, n,
                               &Y, T, n,
                               &Ya, T, n,
                               &YPY, n, n,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

#if TDEBUG
    fprintf(stderr, "stock_yogo: pmod->ncoeff=%d, n=%d, K1=%d, K2=%d\n",
            pmod->ncoeff, n, K1, K2);
    printlist(endolist, "endolist, in stock-yogo");
#endif

    if (K1 > 0) {
        /* we have some exogenous regressors */
        gretl_matrix *M;
        int j, ix = 0, iz = 0;

        B2 = gretl_matrix_block_new(&X, T, K1,
                                    &Z, T, K2,
                                    &E, T, K2,
                                    NULL);
        if (B2 == NULL) {
            gretl_matrix_block_destroy(B);
            return E_ALLOC;
        }

        /* form the X and Z matrices */
        for (i=1; i<=ivi->instlist[0]; i++) {
            vi = ivi->instlist[i];
            if (in_gretl_list(ivi->pmod->list, vi)) {
                /* included exogenous var -> X */
                M = X;
                j = ix++;
            } else {
                /* excluded exogenous var -> Z */
                M = Z;
                j = iz++;
            }
            s = 0;
            for (t=pmod->t1; t<=pmod->t2; t++) {
                if (!na(pmod->uhat[t])) {
                    gretl_matrix_set(M, s++, j, dZ[vi][t]);
                }
            }
        }
    }

    /* form Y matrix */
    for (i=0; i<n; i++) {
        vi = ivi->endolist[i+1];
        s = 0;
        for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t])) {
                gretl_matrix_set(Y, s++, i, dZ[vi][t]);
            }
        }
    }

    if (K1 == 0) {
        /* form Y-hat matrix in Ya */
        for (i=0; i<n; i++) {
            vi = ivi->hatlist[i+1];
            s = 0;
            for (t=pmod->t1; t<=pmod->t2; t++) {
                if (!na(pmod->uhat[t])) {
                    gretl_matrix_set(Ya, s++, i, dZ[vi][t]);
                }
            }
        }
    } else {
        int bdim = (K1 > K2)? K1 : K2;
        gretl_matrix *b;

        /* the leading dimension of b is the number of regressors */
        b = gretl_matrix_alloc(bdim, K2);
        if (b == NULL) {
            err = E_ALLOC;
        }

        if (!err) {
            /* partial X out of Y: Y <- M_x Y */
            gretl_matrix_reuse(b, K1, n);
            gretl_matrix_reuse(E, -1, n);
            err = gretl_matrix_multi_ols(Y, X, b, E, NULL);
            gretl_matrix_copy_values(Y, E);
        }

        if (!err) {
            /* partial X out of Z: Z <- M_x Z */
            gretl_matrix_reuse(b, K1, K2);
            gretl_matrix_reuse(E, -1, K2);
            err = gretl_matrix_multi_ols(Z, X, b, E, NULL);
            gretl_matrix_copy_values(Z, E);
        }

        if (!err) {
            /* form projection P_z Y, in Ya */
            gretl_matrix_reuse(b, K2, n);
            gretl_matrix_reuse(E, -1, n);
            err = gretl_matrix_multi_ols(Y, Z, b, E, NULL);
            gretl_matrix_copy_values(Ya, Y);
            gretl_matrix_subtract_from(Ya, E);
        }

        gretl_matrix_free(b);
    }

    if (!err) {
        /* form Y' P_z Y */
        err = gretl_matrix_multiply_mod(Y, GRETL_MOD_TRANSPOSE,
                                        Ya, GRETL_MOD_NONE,
                                        YPY, GRETL_MOD_NONE);
    }

    /* now write first-stage residuals into Ya */
    for (i=0; i<n && !err; i++) {
        vi = ivi->endolist[i+1];
        vj = ivi->hatlist[i+1];
        s = 0;
        for (t=pmod->t1; t<=pmod->t2; t++) {
            if (!na(pmod->uhat[t])) {
                gretl_matrix_set(Ya, s++, i, dZ[vi][t] - dZ[vj][t]);
            }
        }
    }

    if (!err) {
        /* form S = \hat{\Sigma}_{vv} */
        err = gretl_matrix_multiply_mod(Y, GRETL_MOD_TRANSPOSE,
                                        Ya, GRETL_MOD_NONE,
                                        S, GRETL_MOD_NONE);
        if (!err) {
            gretl_matrix_xtr_symmetric(S);
            gretl_matrix_divide_by_scalar(S, T - K1 - K2);
            for (i=0; i<S->rows; i++) {
                if (gretl_matrix_get(S, i, i) <= 0) {
                    err = E_SINGULAR;
                    break;
                }
            }
        }
    }

    if (!err) {
        double rc = gretl_symmetric_matrix_rcond(S, &err);

        if (!err && (na(rc) || rc < 1.0e-7)) {
            err = E_SINGULAR;
        }
    }

    if (!err) {
        /* S^{-1/2}: invert S and Cholesky-decompose */
        err = gretl_invert_symmetric_matrix(S);
        if (!err) {
            err = gretl_matrix_cholesky_decomp(S);
        }
    }

    if (!err) {
        /* finally, form big sandwich */
        err = gretl_matrix_qform(S, GRETL_MOD_TRANSPOSE,
                                 YPY, G, GRETL_MOD_NONE);
        if (!err) {
            gretl_matrix_divide_by_scalar(G, K2);
        }
    }

    if (!err) {
        double gmin = gretl_symm_matrix_lambda_min(G, &err);

        if (!err) {
            gretl_model_set_double(pmod, "gmin", gmin);
            gretl_model_set_int(pmod, "n", n);
            gretl_model_set_int(pmod, "K2", K2);
        }
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(B2);

#if 0
    if (err) {
        fprintf(stderr, "compute_stock_yogo: err = %d\n", err);
    }
#endif

    return err;
}

/* Calculate first-stage F-test as diagnostic for weak instruments.
   This variant can only handle one endogenous regressor, but it
   supports robust estimation if the model in question uses a robust
   covariance estimator.
*/

static int compute_first_stage_F (iv_info *ivi, gretlopt opt,
				  DATASET *dset)
{
    MODEL mod1;
    int *list1 = NULL;
    int *flist = NULL;
    double F;
    int i, ev, vi;
    int err = 0;

    gretl_model_init(&mod1, dset);
    ev = ivi->endolist[1];

    /* The regression list, list1, has endogenous regressor ev as the
       dependent variable and includes all instruments; the list for
       omission, flist, includes all the excluded instruments.
    */

    list1 = gretl_list_append_term(&list1, ev);
    if (list1 == NULL) {
        return E_ALLOC;
    }

    for (i=1; i<=ivi->instlist[0]; i++) {
        vi = ivi->instlist[i];
        list1 = gretl_list_append_term(&list1, vi);
        if (list1 == NULL) {
            err = E_ALLOC;
            break;
        }
        if (!in_gretl_list(ivi->reglist, vi)) {
            flist = gretl_list_append_term(&flist, vi);
            if (flist == NULL) {
                err = E_ALLOC;
                break;
            }
        }
    }

    if (!err) {
        gretlopt myopt = (opt & OPT_R)? OPT_R : OPT_NONE;

        mod1 = lsq(list1, dset, OLS, myopt | OPT_A);
        err = mod1.errcode;
        if (err) {
            fprintf(stderr, "compute_first_stage F: lsq failed\n");
        }
    }

    if (!err) {
        F = wald_omit_F(flist, &mod1);
        if (na(F)) {
            err = 1;
        }
    }

    if (!err) {
        gretl_model_set_double(ivi->pmod, "stage1-F", F);
        gretl_model_set_int(ivi->pmod, "stage1-dfn", flist[0]);
        gretl_model_set_int(ivi->pmod, "stage1-dfd", mod1.dfd);
        if (!(opt & OPT_R)) {
            /* flag lookup of Stock-Yogo critical values */
            gretl_model_set_double(ivi->pmod, "gmin", F);
        }
    }

    clear_model(&mod1);
    free(list1);
    free(flist);

    return err;
}

static int tsls_adjust_sample (const int *list, DATASET *dset,
                               char **pmask)
{
    int i, t, t1min = dset->t1, t2max = dset->t2;
    char *mask = NULL;
    int T, vi, missobs;
    int err = 0;

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
        missobs = 0;
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || vi == LISTSEP) {
                continue;
            }
            if (na(dset->Z[vi][t])) {
                missobs = 1;
                break;
            }
        }
        if (missobs) {
            t1min++;
        } else {
            break;
        }
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
        missobs = 0;
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || vi == LISTSEP) {
                continue;
            }
            if (na(dset->Z[vi][t])) {
                missobs = 1;
                break;
            }
        }
        if (missobs) {
            t2max--;
        } else {
            break;
        }
    }

    /* count missing values in mid-range of data */
    missobs = 0;
    for (t=t1min; t<=t2max; t++) {
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || vi == LISTSEP) {
                continue;
            }
            if (na(dset->Z[vi][t])) {
                missobs++;
                break;
            }
        }
    }

    T = t2max - t1min + 1;

    if (missobs == T) {
        err = E_MISSDATA;
    } else if (missobs > 0) {
        mask = calloc(T, 1); /* all NUL bytes */
        if (mask == NULL) {
            err = E_ALLOC;
        } else {
            for (t=t1min; t<=t2max; t++) {
                for (i=1; i<=list[0]; i++) {
                    vi = list[i];
                    if (vi == 0 || vi == LISTSEP) {
                        continue;
                    }
                    if (na(dset->Z[vi][t])) {
                        mask[t - t1min] = 1;
                        break;
                    }
                }
            }
        }
    }

#if TDEBUG
    fprintf(stderr, "tsls_adjust_sample: t1=%d, t2=%d, missobs=%d, ok obs=%d\n",
            t1min, t2max, missobs, t2max - t1min + 1 - missobs);
#endif

    if (!err) {
        dset->t1 = t1min;
        dset->t2 = t2max;
    }

    *pmask = mask;

    return err;
}

/**
 * ivreg_process_lists:
 * @list: original specification list.
 * @reglist: location to receive regression list.
 * @instlist: location to receive list of instruments.
 *
 * Split the incoming list into its two components and perform
 * some basic checks; if the checks fail the two created
 * lists are destroyed.
 *
 * Returns: 0, on success, non-zero error code on failure.
 */

int ivreg_process_lists (const int *list, int **reglist, int **instlist)
{
    int *rlist = NULL;
    int *zlist = NULL;
    int i, oc, err;

    err = gretl_list_split_on_separator(list, &rlist, &zlist);
    if (err) {
        return err;
    }

    if (rlist[0] < 2 || zlist == NULL || zlist[0] < 1) {
	/* The regression list must have at least two members, and the
	   instrument list must not be empty
	*/
        err = E_ARGS;
    } else {
        for (i=1; i<=zlist[0] && !err; i++) {
            if (zlist[i] == rlist[1]) {
                gretl_errmsg_set(_("You can't use the dependent variable as an instrument"));
                err = E_DATA;
            }
        }
    }

    if (!err && in_gretl_list(rlist, 0) > 1 && !in_gretl_list(zlist, 0)) {
	/* Don't treat the constant as if it were endogenous: if
	   it's in the regressor list ensure it's also in the
	   instrument list.
	*/
	err = zlist_prepend_const(&zlist);
    }

    if (!err && (oc = zlist[0] - rlist[0] + 1) < 0) {
	gretl_errmsg_sprintf(_("The order condition for identification is not satisfied.\n"
			       "At least %d more instruments are needed."), -oc);
	err = E_DATA;
    }

if (err) {
    free(rlist);
    free(zlist);
 } else {
    *reglist = rlist;
    *instlist = zlist;
 }

    return err;
}

/**
 * tsls:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: may contain %OPT_R for robust VCV, %OPT_N for no df correction,
 * %OPT_A if this is an auxiliary reression, %OPT_E if we're
 * estimating one equation within a system, %OPT_H to
 * add "hatlist" to model even under %OPT_E, %OPT_X to skip
 * the set of tests that are usually calculated. Also %OPT_C
 * for clustered standard errors.
 *
 * Estimate the model given in @list by means of Two-Stage Least
 * Squares.  If %OPT_E is given, fitted values from the first-stage
 * regressions are saved with the model.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tsls (const int *list, DATASET *dset, gretlopt opt)
{
    MODEL tsls;
    iv_info *ivi = NULL;
    int orig_t1 = dset->t1, orig_t2 = dset->t2;
    int orig_nvar = dset->v;
    int sysest = (opt & OPT_E);
    int no_tests = 0;
    int no_hausman = 0;
    int nendo = 0;
    int nreg = 0;
    int i, err = 0;

    /* initialize model in case we bail out early on error */
    gretl_model_init(&tsls, dset);

    /* allocate container struct */
    ivi = iv_info_new(&tsls);
    if (ivi == NULL) {
	tsls.errcode = E_ALLOC;
	return tsls;
    }

    if (opt & OPT_X) {
        /* supplementary tests wanted? */
        const char *s = get_optval_string(IVREG, OPT_X);

        if (s != NULL && !strcmp(s, "hausman")) {
            /* skip the hausman test only */
            no_hausman = 1;
        } else {
            /* skip all extra tests */
            no_tests = 1;
        }
    }

    /* reglist: dependent var plus list of regressors
       instlist: list of instruments
    */
    tsls.errcode = ivreg_process_lists(list, &ivi->reglist, &ivi->instlist);
    if (tsls.errcode) {
	iv_info_destroy(ivi);
        return tsls;
    }

    /* adjust sample range for missing observations */
    err = tsls_adjust_sample(list, dset, &ivi->missmask);

    if (!err) {
        /* allocate second stage regression list */
        ivi->s2list = gretl_list_copy(ivi->reglist);
        if (ivi->s2list == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        goto bailout;
    }

    /* in case we drop any instruments as redundant */
    if (sysest) {
        ivi->exolist = gretl_list_copy(ivi->instlist);
    }

    /* drop any vars that are all zero, and reshuffle the constant
       into first position among the independent vars
    */
    tsls_omitzero(ivi->reglist, ivi->missmask, dset);
    reglist_check_for_const(ivi->reglist, dset);
    tsls_omitzero(ivi->instlist, ivi->missmask, dset);

#if TDEBUG
    printlist(ivi->reglist, "reglist");
    printlist(ivi->instlist, "instlist");
#endif

    /* Determine the list of variables (endolist) for which we need to
       obtain fitted values in the first stage.  Note that we accept
       the condition where there are no such variables (and hence
       TSLS is in fact redundant); this may sometimes be required
       when tsls is called in the context of system estimation.
    */
    ivi->endolist = tsls_make_endolist(ivi->reglist,
				       ivi->instlist,
				       &err);
    if (err) {
        goto bailout;
    }

#if TDEBUG
    printlist(ivi->endolist, "endolist, in tsls");
#endif

    if (ivi->endolist != NULL) {
        nendo = ivi->endolist[0];
        ivi->hatlist = gretl_list_copy(ivi->endolist);
        if (ivi->hatlist == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        err = tsls_add_Q(ivi, dset);
    }

    if (err) {
        goto bailout;
    }

    /* check (again) for order condition: we do this after tsls_Q in case
       any instruments get dropped at that stage
    */
    nreg = ivi->reglist[0] - 1;
    ivi->overid = ivi->instlist[0] - nreg;
    if (ivi->overid < 0) {
        gretl_errmsg_sprintf(_("The order condition for identification is not satisfied.\n"
                               "At least %d more instruments are needed."), -ivi->overid);
        err = E_DATA;
        goto bailout;
    }

#if TDEBUG
    fprintf(stderr, "nreg = %d, OverIdRank = %d\n", nreg, ivi->overid);
#endif

    /* allocate storage for fitted vars (etc.), if needed */
    if (nendo > 0) {
        err = dataset_add_series(dset, nendo);
        if (!err) {
            ivi->g = gretl_vector_alloc(ivi->Q->cols);
            if (ivi->g == NULL) {
                err = E_ALLOC;
            }
        }
        if (err) {
            goto bailout;
        }
    }

    /* Deal with the variables for which instruments are needed: loop
       across the list of variables to be instrumented (endolist),
       form the fitted values as QQ'x_i, and add these fitted values
       into the data array, dset->Z.
    */

    for (i=0; i<nendo; i++) {
        int v0 = ivi->endolist[i+1];
        int v1 = orig_nvar + i;

        /* form xhat_i = QQ'x_i */
        err = tsls_form_xhat(ivi, dset, v0, v1);
        if (err) {
            goto bailout;
        }

        /* substitute v1 into the right place in the second-stage
           regression list */
        replace_list_element(ivi->s2list, v0, v1);

        /* and update hatlist */
        ivi->hatlist[i+1] = v1;
    }

    /* second-stage regression */
    tsls = lsq(ivi->s2list, dset, OLS, (sysest)? OPT_Z : OPT_NONE);
    if (tsls.errcode) {
        fprintf(stderr, "tsls, stage 2: tsls.errcode = %d\n", tsls.errcode);
        goto bailout;
    }

    if (tsls.list[0] < ivi->s2list[0]) {
        /* Were collinear regressors dropped? If so, adjustments are needed */
        ivi->overid += ivi->s2list[0] - tsls.list[0];
        reglist_remove_redundant_vars(ivi);
        if (ivi->endolist != NULL) {
            nendo = ivi->endolist[0];
        }
    }

    if (ivi->instlist != NULL && ivi->instlist[0] > 0) {
        /* record the number of instruments used */
        gretl_model_set_int(&tsls, "ninst", ivi->instlist[0]);
    }

    if (!no_tests) {
        if (!sysest || (opt & OPT_H)) {
            if (nendo == 1) {
                /* handles robust estimation, for single endogenous regressor */
                compute_first_stage_F(ivi, opt, dset);
            } else if (nendo > 0 && !(opt & OPT_R)) {
                /* at present, only handles case of i.i.d. errors */
                compute_stock_yogo(ivi, dset);
            }
        }
    }

    if (nendo > 0) {
        /* special: we need to use the original RHS vars to compute
           residuals and associated statistics */
        tsls_residuals(ivi, opt, dset);
    }

    if (opt & OPT_C) {
        /* clustered implies robust */
        opt |= OPT_R;
    }

    if (opt & OPT_R) {
        /* robust standard errors called for */
        qr_tsls_vcv(&tsls, dset, opt);
        if (tsls.errcode) {
            fprintf(stderr, "qr_tsls_vcv: err = %d\n", tsls.errcode);
            goto bailout;
        }
    }

    if (nendo > 0) {
        /* compute additional statistics (R^2, F, etc.) */
        tsls_extra_stats(ivi, opt, dset);
    }

    if (!sysest && !no_tests) {
        if (!no_hausman && nendo > 0 && ivi->hatlist != NULL) {
	    if (opt & OPT_M) {
		tsls_matrix_hausman_test(ivi, opt, dset);
	    } else {
		tsls_hausman_test(ivi, opt, dset);
	    }
        }
        if (ivi->overid > 0) {
            ivreg_sargan_test(ivi, dset);
        } else if (nendo > 0) {
            tsls_loglik(ivi, dset);
        }
    }

    /* set command code on the model */
    tsls.ci = IVREG;

    /* write the full 2SLS list (dep. var. and regressors, followed by
       instruments, possibly purged of redundant elements) into the model
       for future reference
    */
    tsls_recreate_full_list(ivi);

    if (ivi->idroplist != NULL) {
        gretl_model_set_list_as_data(&tsls, "inst_droplist", ivi->idroplist);
	ivi->idroplist = NULL; /* don't free */
    }

 bailout:

    if (err && !tsls.errcode) {
        tsls.errcode = err;
    }

    if (!tsls.errcode && nendo > 0) {
        if (sysest) {
            /* save first-stage fitted values */
            tsls_save_data(ivi, dset);
        }
        if (!sysest || (opt & OPT_H)) {
            /* save list of endogenous regressors on model */
            gretl_model_set_list_as_data(&tsls, "endolist", ivi->endolist);
            ivi->endolist = NULL; /* model takes ownership */
            gretl_model_set_list_as_data(&tsls, "instlist", ivi->instlist);
            ivi->instlist = NULL; /* ditto */
        }
    }

    /* delete any first-stage fitted values from dataset */
    dataset_drop_last_variables(dset, dset->v - orig_nvar);

    iv_info_destroy(ivi);

    if ((opt & OPT_A) || tsls.errcode) {
        model_count_minus(&tsls); /* OK? */
    }

    if (!tsls.errcode && !(opt & OPT_N)) {
        gretl_model_set_int(&tsls, "dfcorr", 1);
    }

    /* restore original sample range */
    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    return tsls;
}
