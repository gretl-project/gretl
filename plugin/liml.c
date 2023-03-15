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
#include "gretl_matrix.h"
#include "matrix_extra.h"
#include "system.h"
#include "sysml.h"

#define LDEBUG 0

/* Note: it sort-of seems that to produce proper restricted LIML
   estimates one would somehow have to impose appropriate restrictions
   at the stage of the computations below, i.e. where the matrices of
   residuals E0 and E1 are generated.  These are residuals from the
   regression of the dependent and endogenous RHS variables on the
   exogenous vars (first just the included exogenous vars, then all of
   the instruments).  But since we're calculating E0 and E1
   equation-by-equation, we can't impose any cross-equation
   restrictions (and it's not clear to me what they would look like,
   anyway).  This means that the E0 and E1 estimates will be
   invariant, for a given equation, regardless of whether or not we're
   imposing restrictions at the level of the subsequent solution for
   the k-class estimator.  So the minimum eigenvalue and
   log-likelihood, for each equation, will also be invariant with
   respect to any restrictions.  This seems inconsistent.

   But maybe I'm wrong.  Clearly, the invariance mentioned above would
   produce nonsense if we were trying to conduct an LR test of the
   restrictions, but in fact we do an F-test, based on the covariance
   matrix of the unrestricted LIML estimates.  So perhaps it's OK...
*/

/* compose E0 or E1 as in Greene, 4e, p. 686, looping across the
   endogenous vars in the model list
*/

static int resids_to_E (gretl_matrix *E, MODEL *lmod, int *reglist,
			const int *exlist, const int *list,
			DATASET *dset)
{
    int i, t, j = 0;
    int T = E->rows;
    int t1 = dset->t1;
    int err = 0;

    for (i=1; i<=list[0] && !err; i++) {
	if (in_gretl_list(exlist, list[i])) {
	    continue;
	}
        /* set the dependent variable */
	reglist[1] = list[i];

	if (reglist[0] == 1) {
	    /* null model! */
	    int v = reglist[1];

            for (t=0; t<T; t++) {
                gretl_matrix_set(E, t, j, dset->Z[v][t + t1]);
            }
            j++;
	    continue;
	}

	/* regress on the specified set of instruments */
	*lmod = lsq(reglist, dset, OLS, OPT_A);
        err = lmod->errcode;
	if (!err) {
            /* put residuals into appropriate column of E and
               increment the column */
            for (t=0; t<T; t++) {
                gretl_matrix_set(E, t, j, lmod->uhat[t + t1]);
            }
            j++;
        }
	clear_model(lmod);
    }

    return err;
}

/* construct the regression list for the auxiliary regressions
   needed as a basis for LIML */

static int *
liml_make_reglist (const equation_system *sys,
		   DATASET *dset, const int *list,
		   const int *exlist, int *k,
		   int *err)
{
    int nexo = exlist[0];
    int *reglist;
    int i, j, vi;

    reglist = gretl_list_new(nexo + 1);
    if (reglist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

#if LDEBUG
    fprintf(stderr, "liml_make_reglist: found %d exog vars\n", nexo);
    printlist(exlist, "exog list");
#endif

    /* at first, put all _included_ exog vars in reglist */
    *k = 1;
    reglist[0] = 1;
    reglist[1] = 0;
    j = 2;
    for (i=2; i<=list[0]; i++) {
	vi = list[i];
	if (in_gretl_list(exlist, vi)) {
	    reglist[0] += 1;
	    reglist[j++] = vi;
	} else {
	    /* an endogenous var */
	    *k += 1;
	}
    }

#if LDEBUG
    printlist(reglist, "liml_make_reglist, reglist");
#endif

    return reglist;
}

/* set the special LIML k-class data on the model: these data will be
   retrieved when calculating the LIML coefficients and their
   covariance matrix (in sysest.c)
*/

static int
liml_set_model_data (MODEL *pmod, const gretl_matrix *E,
		     const int *exlist, const int *list,
		     int T, double lmin, DATASET *dset)
{
    double *Xi = NULL;
    double *ymod = NULL;
    double yt, xit, eit;
    int pos, m;
    int i, vi, j, s, t;
    int err = 0;

    pos = gretl_list_separator_position(list);
    m = (pos > 0)? (pos - 2) : list[0] - 1;

    ymod = malloc(dset->n * sizeof *ymod);
    if (ymod == NULL) {
	return 1;
    }

    for (t=0; t<dset->n; t++) {
	ymod[t] = NADBL;
    }

    for (t=0; t<T; t++) {
	s = t + dset->t1;
	yt = dset->Z[list[1]][s];
	eit = gretl_matrix_get(E, t, 0);
	ymod[t + dset->t1] = yt - lmin * eit;
	j = 1;
	for (i=0; i<m; i++) {
	    vi = list[i+2];
	    if (in_gretl_list(exlist, vi)) {
		continue;
	    }
	    Xi = model_get_Xi(pmod, dset, i);
	    if (Xi == NULL) {
		err = 1;
		break;
	    }
	    xit = dset->Z[vi][s];
	    eit = gretl_matrix_get(E, t, j++);
	    Xi[s] = xit - lmin * eit;
	}
	if (err) break;
    }

    if (!err) {
	err = gretl_model_set_data(pmod, "liml_y", ymod,
				   GRETL_TYPE_DOUBLE_ARRAY,
				   dset->n * sizeof *ymod);
    }

    if (err) {
	free(ymod);
    }

    return err;
}

static double liml_get_ldet (gretl_matrix *W1, int *err)
{
    double ret = NADBL;
    char *mask;

    /* allow for the possibility that W1 is rank-deficient? */
    mask = gretl_matrix_rank_mask(W1, err);
    if (mask != NULL) {
	fprintf(stderr, "note: LIML W1 is rank deficient\n");
        *err = gretl_matrix_cut_rows_cols(W1, mask);
    }
    if (!*err) {
        ret = gretl_matrix_log_determinant(W1, err);
    }

    return ret;
}

static int liml_eqn_get_lists (equation_system *sys, int eq,
                               int **plist, int **pexlist,
                               int *freelists)
{
    int *list = system_get_list(sys, eq);
    int err = 0;

    if (gretl_list_has_separator(list)) {
        /* got a TSLS-style list */
        err = gretl_list_split_on_separator(list, plist, pexlist);
        *freelists = 1;
    } else {
        *plist = list;
        *pexlist = system_get_instr_vars(sys);
    }

    return err;
}

static int liml_do_equation (equation_system *sys, int eq,
			     DATASET *dset, PRN *prn)
{
    int *list = NULL;
    int *exlist = NULL;
    int *reglist = NULL;
    gretl_matrix *E = NULL;
    gretl_matrix *W0 = NULL;
    gretl_matrix *W1 = NULL;
    double lmin = 1.0;
    MODEL *pmod;
    MODEL lmod;
    int freelists = 0;
    int idf, i, k;
    int T = sys->T;
    int err = 0;

    err = liml_eqn_get_lists(sys, eq, &list, &exlist, &freelists);
    if (err) {
        return err;
    }

#if LDEBUG
    fprintf(stderr, "\nWorking on equation for %s\n", dset->varname[list[1]]);
#endif

    /* get pointer to model (initialized via TSLS) */
    pmod = system_get_model(sys, eq);

    /* degrees of freedom for over-identification test: total
       exogenous vars minus the number of parameters in the equation
       (unless we're estimating subject to specified restrictions, in
       which case we skip the usual over-id test)
    */
    if (system_n_restrictions(sys) == 0) {
	idf = exlist[0] - pmod->ncoeff;
    } else {
	idf = -1;
	gretl_model_set_int(pmod, "restricted", 1);
    }

    /* first make regression list using only included instruments */
    reglist = liml_make_reglist(sys, dset, list, exlist, &k, &err);
    if (err) {
	if (freelists) {
	    free(list);
	    free(exlist);
	}
	return err;
    }

#if LDEBUG
    fprintf(stderr, "number of endogenous vars in equation: k = %d\n", k);
#endif

    E = gretl_matrix_alloc(T, k);
    W0 = gretl_matrix_alloc(k, k);
    W1 = gretl_matrix_alloc(k, k);

    if (E == NULL || W0 == NULL || W1 == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = resids_to_E(E, &lmod, reglist, exlist, list, dset);
    if (!err) {
	err = gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
					E, GRETL_MOD_NONE,
					W0, GRETL_MOD_NONE);
    }
#if LDEBUG
    if (!err) gretl_matrix_print(W0, "W0");
#endif

    if (!err) {
	/* re-make the regression list using all instruments */
	reglist[0] = 1 + exlist[0];
	for (i=2; i<=reglist[0]; i++) {
	    reglist[i] = exlist[i-1];
	}
	err = resids_to_E(E, &lmod, reglist, exlist, list, dset);
    }
    if (!err) {
	err = gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
					E, GRETL_MOD_NONE,
					W1, GRETL_MOD_NONE);
    }
#if LDEBUG
    if (!err) gretl_matrix_print(W1, "W1");
#endif

    if (!err) {
        /* determine the minimum eigenvalue of W1^{-1} * W0 */
        gretl_matrix *L = gretl_gensymm_eigenvals(W1, W0, NULL, &err);

        if (!err) {
            lmin = 1.0 / L->val[k-1];
        }
        gretl_matrix_free(L);
    }

    if (!err) {
	gretl_model_set_double(pmod, "lmin", lmin);
	gretl_model_set_int(pmod, "idf", idf);
#if LDEBUG
	fprintf(stderr, "lmin = %g, idf = %d\n", lmin, idf);
#endif
	err = liml_set_model_data(pmod, E, exlist, list, T,
				  lmin, dset);
	if (err) {
	    fprintf(stderr, "error in liml_set_model_data()\n");
	}
    }

    if (!err) {
	/* compute and set log-likelihood, etc. */
	double ldet = liml_get_ldet(W1, &err);

	if (na(ldet)) {
	    pmod->lnL = NADBL;
	} else {
	    /* Davidson and MacKinnon, ETM, p. 538 */
	    pmod->lnL = -(T/2.0) * (sys->neqns * LN_2_PI + log(lmin) + ldet);
	}
	mle_criteria(pmod, 0);
    }

 bailout:

    gretl_matrix_free(E);
    gretl_matrix_free(W0);
    gretl_matrix_free(W1);

    free(reglist);
    if (freelists) {
	free(list);
	free(exlist);
    }

    return err;
}

/* Driver function for LIML: calculate the minimum eigenvalue per
   equation, and set the suitably transformed data on the respective
   models
*/

int liml_driver (equation_system *sys, DATASET *dset, PRN *prn)
{
    int i, err = 0;

#if LDEBUG
    fprintf(stderr, "\n *** liml driver called: sys = %p\n", (void *) sys);
#endif

    for (i=0; i<sys->neqns && !err; i++) {
#if LDEBUG > 1
	if (prn != NULL) {
	    printmodel(system_get_model(sys, i), dset, OPT_NONE, prn);
	}
#endif
	err = liml_do_equation(sys, i, dset, prn);
    }

    return err;
}
