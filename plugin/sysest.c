/*
 *  Copyright (c) 2002-2004 by Allin Cottrell
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
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"
#include "system.h"

/* fiml.c */
extern int fiml_driver (gretl_equation_system *sys, double ***pZ, 
			gretl_matrix *sigma, DATAINFO *pdinfo, 
			PRN *prn);

/* liml.c */
extern int liml_driver (gretl_equation_system *sys, double ***pZ, 
			gretl_matrix *sigma, DATAINFO *pdinfo, 
			PRN *prn);

static void 
print_system_vcv (const gretl_matrix *m, int triangle, PRN *prn)
{
    int jmax = (triangle)? 1 : m->cols;
    char numstr[16];
    double x;
    int i, j;

    pprintf(prn, "%s\n(%s)\n\n",
	    _("Cross-equation VCV for residuals"),
	    _("correlations above the diagonal"));

    for (i=0; i<m->rows; i++) {
	for (j=0; j<jmax; j++) {
	    pprintf(prn, "%#10.5g ", gretl_matrix_get(m, i, j));
	}
	for (j=jmax; j<m->cols; j++) {
	    x = gretl_matrix_get(m, i, i) * gretl_matrix_get(m, j, j);
	    x = sqrt(x);
	    x = gretl_matrix_get(m, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    pprintf(prn, "%11s", numstr);
	}
	pputc(prn, '\n');
	if (triangle && jmax < m->cols) {
	    jmax++;
	}
    }

    pputc(prn, '\n');
}

static void kronecker_place (gretl_matrix *X, 
			     const gretl_matrix *M,
			     int startrow, int startcol,
			     double scale)
{
    int i, j;
    int row, col;
    double x;
    
    for (i=0; i<M->rows; i++) {
	row = startrow + i;
	for (j=0; j<M->cols; j++) {
	    col = startcol + j;
	    x = gretl_matrix_get(M, i, j);
	    gretl_matrix_set(X, row, col, x * scale);
	}
    }
}

static int make_sys_X_block (gretl_matrix *X,
			     const MODEL *pmod,
			     const double **Z, 
			     int t1, int systype)
{
    int i, t;
    const double *Xi;

    X->cols = pmod->ncoeff;

    for (i=0; i<X->cols; i++) {
	if (systype == THREESLS || systype == FIML) {
	    Xi = tsls_get_Xi(pmod, Z, i);
	} else {
	    Xi = Z[pmod->list[i+2]];
	}
	if (Xi == NULL) {
	    return 1;
	}
	for (t=0; t<X->rows; t++) {
	    gretl_matrix_set(X, t, i, Xi[t+t1]);
	}
    }

    return 0;
}

/* populate the cross-equation covariance matrix based on
   the per-equation residuals
*/

static int
gls_sigma_from_uhat (gretl_matrix *sigma, const gretl_matrix *e, 
		     int m, int T)
{
    int i, j, t;
    double xx;

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=i; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gretl_matrix_get(e, i, t) * gretl_matrix_get(e, j, t);
	    }
	    xx /= T;
	    gretl_matrix_set(sigma, i, j, xx);
	    if (j != i) {
		gretl_matrix_set(sigma, j, i, xx);
	    }
	}
    }

    return 0;
}

static void sys_resids (MODEL *pmod, const double **Z, gretl_matrix *uhat)
{
    double yh;
    int i, t;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	yh = 0.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    yh += pmod->coeff[i] * Z[pmod->list[i+2]][t];
	}
	pmod->yhat[t] = yh;
	pmod->uhat[t] = Z[pmod->list[1]][t] - yh;
	/* for cross-equation vcv */
	gretl_matrix_set(uhat, pmod->ID, t - pmod->t1, pmod->uhat[t]);
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    /* df correction could be applied here */
    pmod->sigma = sqrt(pmod->ess / pmod->nobs);
}

static int 
calculate_sys_coefficients (MODEL **models, const double **Z,
			    gretl_matrix *X, gretl_matrix *uhat,
			    double *tmp_y, int m, int mk)
{
    gretl_vector *coeff;
    gretl_matrix *vcv;
    int i, j, start;

    coeff = gretl_vector_alloc(mk);
    if (coeff == NULL) return 1;

    for (i=0; i<mk; i++) {
	gretl_vector_set(coeff, i, tmp_y[i]);
    }

    vcv = gretl_matrix_copy(X);
    if (vcv == NULL) {
	gretl_vector_free(coeff);
	return 1;
    }

    gretl_LU_solve(X, coeff);
    gretl_invert_general_matrix(vcv); 

    start = 0;
    for (i=0; i<m; i++) {
	int nci = (models[i])->ncoeff;

	for (j=0; j<nci; j++) {
	    (models[i])->coeff[j] = gretl_vector_get(coeff, start + j);
	    (models[i])->sderr[j] = 
		sqrt(gretl_matrix_get(vcv, start + j, start + j));
	}
	sys_resids(models[i], Z, uhat);
	start += nci;
    }

    gretl_vector_free(coeff);
    gretl_matrix_free(vcv);

    return 0;
}

/* FIXME below: naming and labeling of added variables */

static void add_results_to_dataset (gretl_equation_system *sys, 
				    MODEL *pmod, int i, int *pj,
				    double **Z, DATAINFO *pdinfo,
				    int systype)
{
    int t;

    if (system_save_uhat(sys)) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		Z[*pj][t] = NADBL;
	    } else {
		Z[*pj][t] = pmod->uhat[t];
	    }
	}
	sprintf(pdinfo->varname[*pj], "uhat_s%02d", i + 1);
	if (systype == SUR) {
	    sprintf(VARLABEL(pdinfo, *pj), _("SUR residual, equation %d"), 
		    i + 1);
	} else if (systype == THREESLS) {
	    sprintf(VARLABEL(pdinfo, *pj), _("3SLS residual, equation %d"), 
		    i + 1);
	} else {
	    sprintf(VARLABEL(pdinfo, *pj), "system residual, equation %d", 
		    i + 1);
	}
	*pj += 1;
    }

    if (system_save_yhat(sys)) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		Z[*pj][t] = NADBL;
	    } else {
		Z[*pj][t] = pmod->yhat[t];
	    }
	}	
	sprintf(pdinfo->varname[*pj], "yhat_s%02d", i + 1);
	if (systype == SUR) {
	    sprintf(VARLABEL(pdinfo, *pj), _("SUR fitted value, equation %d"), 
		    i + 1);
	} else if (systype == THREESLS) {
	    sprintf(VARLABEL(pdinfo, *pj), _("3SLS fitted value, equation %d"), 
		    i + 1);
	} else {
	    sprintf(VARLABEL(pdinfo, *pj), "system fitted value, equation %d", 
		    i + 1);
	}	    
	*pj += 1;
    }
}

static int in_list (const int *list, int k)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (k == list[i]) return 1;
    }

    return 0;
}

static int *
system_model_list (gretl_equation_system *sys, int i, int *freeit)
{
    int systype = system_get_type(sys);
    int *list = NULL;

    *freeit = 0;

    if (systype == SUR || systype == THREESLS) {
	list = system_get_list(sys, i);
    }

    if (systype == THREESLS) {
	/* is list already in tsls form? */
	if (list != NULL && !in_list(list, LISTSEP)) {
	    list = NULL;
	}
    }

    if (systype == FIML || systype == LIML ||
	(systype == THREESLS && list == NULL)) {
	list = compose_tsls_list(sys, i);
	*freeit = 1;
    }

    return list;
}

static void
print_fiml_overidentification_test (const gretl_equation_system *sys, 
				    PRN *prn)
{
    int df = system_get_df(sys);

    if (df > 0) {
	double ll = system_get_ll(sys);
	double llu = system_get_llu(sys);
	double X2;

	/* let's not print rubbish */
	if (na(ll) || na(llu) || ll == 0.0 || llu == 0.0) {
	    return;
	}

	X2 = 2.0 * (llu - ll);

	pprintf(prn, "\n%s:\n", _("LR over-identification test"));
	pprintf(prn, "  %s = %g\n", _("Restricted log-likelihood"), ll);
	pprintf(prn, "  %s = %g\n", _("Unrestricted log-likelihood"), llu);
	pprintf(prn, "  %s(%d) = %g %s %g\n", _("Chi-square"),
		df, X2, _("with p-value"), chisq(X2, df));
    }
}

static int basic_system_allocate (int systype, int m, int T, int mk,
				  MODEL ***models,
				  gretl_matrix **uhat, 
				  gretl_matrix **sigma,
				  gretl_matrix **X)
{
    MODEL **pmods;
    int i, j;

    /* allocate a model for each stochastic equation */
    pmods = malloc(m * sizeof *pmods);
    if (pmods == NULL) {
	return E_ALLOC;
    }
    for (i=0; i<m; i++) {
	pmods[i] = gretl_model_new();
	if (pmods[i] == NULL) {
	    for (j=0; j<i; j++) {
		free_model(pmods[j]);
	    }
	    free(pmods);
	    return E_ALLOC;
	}
    }

    *models = pmods;

    *uhat = gretl_matrix_alloc(m, T);
    if (*uhat == NULL) {
	return E_ALLOC;
    }

    *sigma = gretl_matrix_alloc(m, m);
    if (*sigma == NULL) {
	return E_ALLOC;
    }    
    
    if (systype != LIML) {
	*X = gretl_matrix_alloc(mk, mk);
	if (*X == NULL) {
	    return E_ALLOC;
	}
    }

    return 0;
}

static int 
save_and_print_results (gretl_equation_system *sys, const gretl_matrix *sigma,
			MODEL **models, double ***pZ, DATAINFO *pdinfo,
			PRN *prn)
{
    int systype = system_get_type(sys);
    int m = system_n_equations(sys);
    int i, j = 0;
    int err = 0;

    if (system_save_uhat(sys)) {
	j = pdinfo->v;
	err = dataset_add_vars(m, pZ, pdinfo);
    }
    if (system_save_yhat(sys)) {
	if (j == 0) {
	    j = pdinfo->v;
	}
	err = dataset_add_vars(m, pZ, pdinfo);
    }

    for (i=0; i<m; i++) {
	printmodel(models[i], pdinfo, OPT_NONE, prn);
	add_results_to_dataset(sys, models[i], i, &j, *pZ, pdinfo, systype);
	if (systype == THREESLS || systype == FIML) {
	    tsls_free_data(models[i]);
	}
    }

    if (!err) {
	print_system_vcv(sigma, 1, prn);
	if (systype == FIML) {
	    print_fiml_overidentification_test(sys, prn);
	}
    }

    return err;
}

/* general function that forms the basis for SUR, 3SLS, FIML and LIML
   estimates
*/

int system_estimate (gretl_equation_system *sys, double ***pZ, DATAINFO *pdinfo, 
		     PRN *prn)
{
    int i, j, k, m, T, t;
    int v, l, mk, krow;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int orig_t1 = t1, orig_t2 = t2;
    gretl_matrix *X = NULL;
    gretl_matrix *uhat = NULL;
    gretl_matrix *sigma = NULL;
    gretl_matrix *Xi = NULL, *Xj = NULL, *M = NULL;
    double *tmp_y = NULL;
    MODEL **models = NULL;
    int systype = system_get_type(sys);
    int err = 0;

    /* get uniform sample starting and ending points */
    if (system_adjust_t1t2(sys, &t1, &t2, (const double **) *pZ)) {
	return E_DATA;
    }   

    /* reset sample temporarily */
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    /* number of equations */
    m = system_n_equations(sys);

    /* max indep vars per equation */
    k = system_max_indep_vars(sys);

    /* total indep vars, all equations */
    mk = system_n_indep_vars(sys);

    /* number of observations per series */
    T = t2 - t1 + 1;
    system_set_n_obs(sys, T);

    /* allocate models etc */
    err = basic_system_allocate(systype, m, T, mk, &models,
				&uhat, &sigma, &X);
    if (err) goto cleanup;
	    
    if (systype == FIML || systype == LIML) {
	print_equation_system_info(sys, pdinfo, prn);
    }

    /* first estimate the equations separately, and put the 
       single-equation residuals into the uhat matrix
    */
    for (i=0; i<m; i++) {
	int freeit = 0;
	int *list = system_model_list(sys, i, &freeit);
	MODEL *pmod = models[i];

	if (list == NULL) {
	    err = 1;
	    break;
	}

	if (systype == SUR) {
	    *pmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	} else if (systype == THREESLS || systype == FIML) {
	    *pmod = tsls_func(list, 0, pZ, pdinfo, OPT_S);
	} else if (systype == LIML) {
	    *pmod = tsls_func(list, 0, pZ, pdinfo, OPT_N);
	}

	if (freeit) {
	    free(list);
	}

	if ((err = pmod->errcode)) {
	    fprintf(stderr, "model failed on lists[%d], code=%d\n",
		    i, err);
	    break;
	}

	pmod->ID = i;
	pmod->aux = AUX_SYS;
	gretl_model_set_int(pmod, "systype", systype);

	for (t=0; t<T; t++) {
	    gretl_matrix_set(uhat, i, t, pmod->uhat[t + pdinfo->t1]);
	}
    }

    if (err) goto cleanup;

    if (systype == LIML) {
	system_attach_models(sys, models);
	err = liml_driver(sys, pZ, sigma, pdinfo, prn);
	system_unattach_models(sys);
	if (!err) {
	    err = save_and_print_results(sys, sigma, models,
					 pZ, pdinfo, prn);
	}
	goto cleanup;
    }

    gls_sigma_from_uhat(sigma, uhat, m, T);
    err = gretl_invert_symmetric_matrix(sigma);    
    if (err) {
	goto cleanup;
    }

    Xi = gretl_matrix_alloc(T, k);
    Xj = gretl_matrix_alloc(T, k);
    M = gretl_matrix_alloc(k, k);
    if (Xi == NULL || Xj == NULL || M == NULL) {
	err = E_ALLOC;
	goto cleanup;
    }	

    /* Xi = data matrix for equation i, specified in lists[i] */
    krow = 0;
    for (i=0; i<m && !err; i++) {
	const gretl_matrix *Xk;
	double sij;
	int kcol = 0;

	err = make_sys_X_block(Xi, models[i], (const double **) *pZ, 
			       pdinfo->t1, systype);

	for (j=0; j<m && !err; j++) { 
	    if (i != j) {
		err = make_sys_X_block(Xj, models[j], (const double **) *pZ, 
				       pdinfo->t1, systype);
		Xk = Xj;
	    } else {
		Xk = Xi;
	    }
	    M->rows = Xi->cols;
	    M->cols = Xk->cols;
	    err = gretl_matrix_multiply_mod (Xi, GRETL_MOD_TRANSPOSE,
					     Xk, GRETL_MOD_NONE, 
					     M);
	    sij = gretl_matrix_get(sigma, i, j);
	    kronecker_place (X, M, krow, kcol, sij);
	    kcol += models[j]->ncoeff;
	}

	krow += models[i]->ncoeff;
    }

    if (err) goto cleanup;

    gretl_matrix_free(Xj);
    Xj = NULL;
    gretl_matrix_free(M);
    M = NULL;

    tmp_y = malloc(mk * sizeof *tmp_y);
    if (tmp_y == NULL) {
	err = E_ALLOC;
	goto cleanup;
    }

    /* form Y column vector (m x k) */
    v = 0;
    for (i=0; i<m; i++) { /* loop over the m vertically arranged
			     blocks in the final column vector */
	const double *y;
	double xx;

	make_sys_X_block(Xi, models[i], (const double **) *pZ, 
			 pdinfo->t1, systype);

	for (j=0; j<models[i]->ncoeff; j++) { /* loop over the rows within each of 
						 the m blocks */
	    tmp_y[v] = 0.0;
	    for (l=0; l<m; l++) { /* loop over the m components that
				     must be added to form each element */
		y = (*pZ)[system_get_depvar(sys, l)];
		/* multiply X'[l] into y */
		xx = 0.0;
		for (t=0; t<T; t++) {
		    xx += gretl_matrix_get(Xi, t, j) * y[t + pdinfo->t1];
		}
		xx *= gretl_matrix_get(sigma, i, l);
		tmp_y[v] += xx;
	    }
	    v++;
	}
    }

    /* depending on whether we used OLS or TSLS at the first stage above,
       these estimates will be either SUR or 3SLS
    */
    calculate_sys_coefficients(models, (const double **) *pZ, X, uhat, 
			       tmp_y, m, mk);

    gls_sigma_from_uhat(sigma, uhat, m, T);

    if (systype == FIML) {
	/* set up convenience pointers */
	system_attach_uhat(sys, uhat);
	system_attach_models(sys, models);

	err = fiml_driver(sys, pZ, sigma, pdinfo, prn);

	/* discard convenience pointers */
	system_unattach_uhat(sys);
	system_unattach_models(sys);
    } 

    if (!err) {
	err = save_and_print_results(sys, sigma, models, 
				     pZ, pdinfo, prn);
    }

 cleanup:

    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);
    gretl_matrix_free(M);
    gretl_matrix_free(X);
    gretl_matrix_free(sigma);
    gretl_matrix_free(uhat);

    free(tmp_y);

    if (models != NULL) {
	for (i=0; i<m; i++) {
	    free_model(models[i]);
	}
	free(models);
    }

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return err;
}
