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

    for (i=1; i<=list[0]; i++) {
	int vi = list[i];

	if (in_gretl_list(exlist, vi)) {
	    continue;
	}

	reglist[1] = vi;

	/* regress the given endogenous var on the specified
	   set of instruments */
	*lmod = lsq(reglist, dset, OLS, OPT_A);
	if ((err = lmod->errcode)) {
	    clear_model(lmod);
	    break;
	}

	/* put residuals into appropriate column of E and
	   increment the column */
	for (t=0; t<T; t++) {
	    gretl_matrix_set(E, t, j, lmod->uhat[t + t1]);
	}
	j++;

	clear_model(lmod);
    }

    return err;
}

/* construct the regression list for the auxiliary regressions
   needed as a basis for LIML */

static int *
liml_make_reglist (const equation_system *sys, const int *list, 
		   const int *exlist, int *k, int *err)
{
    int nexo = exlist[0];
    int *reglist;
    int i, j;

    reglist = gretl_list_new(nexo + 1);
    if (reglist == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

#if LDEBUG
    fprintf(stderr, "liml_make_reglist: found %d exog vars\n", nexo);
#endif

    /* at first, put all _included_ exog vars in reglist */
    *k = 1;
    reglist[0] = 1;
    reglist[1] = 0;
    j = 2;
    for (i=2; i<=list[0]; i++) {
	if (in_gretl_list(exlist, list[i])) {
	    reglist[0] += 1;
	    reglist[j++] = list[i];
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

static int liml_do_equation (equation_system *sys, int eq, 
			     DATASET *dset, PRN *prn)
{
    int *list = system_get_list(sys, eq);
    int *exlist = NULL;
    gretl_matrix_block *B = NULL;
    gretl_matrix *E, *W0, *W1, *W2, *Inv;
    gretl_matrix *lambda = NULL;
    double lmin = 1.0;
    MODEL *pmod;
    MODEL lmod;
    int *reglist;
    int freelists = 0;
    int idf, i, k;
    int T = sys->T;
    int err = 0;

#if LDEBUG
    fprintf(stderr, "\nWorking on equation for %s\n", dset->varname[list[1]]);
    printlist(list, "original equation list");
#endif

    if (gretl_list_has_separator(list)) {
	/* we got a TSLS-style list */
	const int *elist = list;

	list = NULL;
	err = gretl_list_split_on_separator(elist, &list, &exlist);
	if (err) {
	    return err;
	}
	freelists = 1;
    } else {
	exlist = system_get_instr_vars(sys);
    }

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
    reglist = liml_make_reglist(sys, list, exlist, &k, &err);
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

    B = gretl_matrix_block_new(&E, T, k,
			       &W0, k, k,
			       &W1, k, k,
			       &W2, k, k,
			       &Inv, k, k,
			       NULL);

    if (B == NULL) {
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
    gretl_matrix_print(W0, "W0");
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
    gretl_matrix_print(W1, "W1");
#endif

    /* form a symmetric matrix that has the same eigenvalues
       as W1^{-1} * W0, and find its minimum eigenvalue
    */

    if (!err) {
	gretl_matrix_copy_values(Inv, W1);
	err = gretl_matrix_cholesky_decomp(Inv) ||
	    gretl_invert_triangular_matrix(Inv, 'L');
    }

    if (!err) {
	err = gretl_matrix_qform(Inv, GRETL_MOD_NONE, W0,
				 W2, GRETL_MOD_NONE);
	if (!err) {
	    lmin = gretl_symm_matrix_lambda_min(W2, &err);
	}
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
	double ldet = gretl_matrix_log_determinant(W1, &err);
	int g = sys->neqns;

	if (err) {
	    pmod->lnL = NADBL;
	} else {
	    /* Davidson and MacKinnon, ETM, p. 538 */
	    pmod->lnL = -(T / 2.0) * (g * LN_2_PI + log(lmin) + ldet);
	}
	mle_criteria(pmod, 0); /* check the "0" (additional params) here */
    }

 bailout:

    free(reglist);
    gretl_matrix_block_destroy(B);
    gretl_matrix_free(lambda);

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
    fprintf(stderr, "\n *** liml driver called\n");
#endif

    for (i=0; i<sys->neqns && !err; i++) {
#if LDEBUG > 1
	printmodel(system_get_model(sys, i), dset, OPT_NONE, prn);
#endif
	err = liml_do_equation(sys, i, dset, prn);
    }

    return err;
}



    
	

    
    
    

    
    
