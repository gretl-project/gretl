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

static int in_list (const int *list, int k)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (k == list[i]) return 1;
    }

    return 0;
}

static void 
print_system_vcv (const gretl_matrix *m, int triangle, PRN *prn)
{
    int i, j, jmax;
    double x;
    char numstr[16];

    jmax = (triangle)? 1 : m->cols;

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
	if (triangle && jmax < m->cols) jmax++;
    }
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

static int
gls_sigma_from_uhat (gretl_matrix *sigma, const gretl_matrix *e, 
		     int m, int T)
{
    int i, j, t;
    double xx;

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gretl_matrix_get(e, i, t) * gretl_matrix_get(e, j, t);
	    }
	    gretl_matrix_set(sigma, i, j, xx / T);
	}
    }

    return 0;
}

static gretl_matrix *
gls_sigma_inverse_from_uhat (const gretl_matrix *e, int m, int T)
{
    gretl_matrix *sigma;

    sigma = gretl_matrix_alloc(m, m);
    if (sigma == NULL) return NULL;

    gls_sigma_from_uhat(sigma, e, m, T);
    gretl_invert_general_matrix(sigma);

    return sigma;
}

static void sys_resids (MODEL *pmod, const double **Z, gretl_matrix *uhat)
{
    int i, t;
    double fit;

    pmod->ess = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	fit = 0.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    fit += pmod->coeff[i] * Z[pmod->list[i+2]][t];
	}
	pmod->yhat[t] = fit;
	pmod->uhat[t] = Z[pmod->list[1]][t] - fit;
	/* for cross-equation vcv */
	gretl_matrix_set(uhat, pmod->ID, t - pmod->t1, pmod->uhat[t]);
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);
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

static void restore_sample (DATAINFO *pdinfo, int t1, int t2)
{
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;
}

static int *
system_model_list (gretl_equation_system *sys, int i, int *freeit)
{
    int systype = system_get_type(sys);

    if (systype == SUR) {
	*freeit = 0;
	return system_get_list(sys, i);
    } 

    if (systype == THREESLS) {
	int *list = system_get_list(sys, i);
	
	if (in_list(list, LISTSEP)) {
	    *freeit = 0;
	    return list;
	}
    }

    if (systype == THREESLS || systype == FIML) {
	*freeit = 1;
	return compose_tsls_list(sys, i);
    }

    return NULL;
}

int system_estimate (gretl_equation_system *sys, double ***pZ, DATAINFO *pdinfo, 
		     PRN *prn)
{
    int i, j, k, m, T, t, v, l, mk, krow;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int orig_t1 = t1, orig_t2 = t2;
    gretl_matrix *X = NULL;
    gretl_matrix *uhat = NULL, *sigma = NULL;
    gretl_matrix *Xi = NULL, *Xj = NULL, *M = NULL;
    double *tmp_y = NULL;
    int systype = system_get_type(sys);
    MODEL **models = NULL;
    int err = 0;

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

    models = malloc(m * sizeof *models);
    if (models == NULL) {
	restore_sample(pdinfo, orig_t1, orig_t2);
	return E_ALLOC;
    }

    for (i=0; i<m; i++) {
	models[i] = gretl_model_new();
	if (models[i] == NULL) {
	    for (j=0; j<i; j++) {
		free_model(models[j]);
	    }
	    restore_sample(pdinfo, orig_t1, orig_t2);
	    return E_ALLOC;
	}
    }

    X = gretl_matrix_alloc(mk, mk);
    uhat = gretl_matrix_alloc(m, T);

    if (X == NULL || uhat == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* just testing for now */
    if (systype == FIML) {
	print_fiml_sys_info(sys, pdinfo, prn);
	pprintf(prn, "Sorry, FIML is not really implmented yet.\n");
    }

    /* first grab the single-equation residuals */
    for (i=0; i<m; i++) {
	int freeit = 0;
	int *list = system_model_list(sys, i, &freeit);

	if (list == NULL) {
	    err = 1;
	    goto bailout;
	}

	if (systype == SUR) {
	    *models[i] = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	} else if (systype == THREESLS || systype == FIML) {
	    *models[i] = tsls_func(list, 0, pZ, pdinfo, OPT_S);
	}

	if (freeit) free(list);

	if ((models[i])->errcode) {
	    fprintf(stderr, "model failed on lists[%d], code=%d\n",
		    i, (models[i])->errcode);
	    err = (models[i])->errcode;
	    goto bailout;
	}

	(models[i])->ID = i;
	(models[i])->aux = AUX_SYS;
	gretl_model_set_int(models[i], "systype", systype);

	for (t=0; t<T; t++) {
	    gretl_matrix_set(uhat, i, t, (models[i])->uhat[t + pdinfo->t1]);
	}
    }

    sigma = gls_sigma_inverse_from_uhat(uhat, m, T);
    if (sigma == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    Xi = gretl_matrix_alloc(T, k);
    Xj = gretl_matrix_alloc(T, k);
    M = gretl_matrix_alloc(k, k);
    if (Xi == NULL || Xj == NULL || M == NULL) {
	err = E_ALLOC;
	goto bailout;
    }	

    /* Xi = data matrix for equation i, specified in lists[i] */
    krow = 0;
    for (i=0; i<m && !err; i++) {
	int kcol = 0;
	gretl_matrix *Xk;

	err = make_sys_X_block(Xi, models[i], (const double **) *pZ, 
			       pdinfo->t1, systype);
	kcol = 0;
	for (j=0; j<m && !err; j++) { 
	    if (i != j) {
		err = make_sys_X_block(Xj, models[j], (const double **) *pZ, 
				       pdinfo->t1, systype);
		Xk = Xj;
	    } else {
		Xk = Xi;
	    }
	    M->rows = Xi->cols;
	    M->cols = Xj->cols;
	    err = gretl_matrix_multiply_mod ((const gretl_matrix *) Xi, 
					     GRETL_MOD_TRANSPOSE,
					     (const gretl_matrix *) Xk, 
					     GRETL_MOD_NONE, 
					     M);
	    kronecker_place (X, (const gretl_matrix *) M,
			     krow, kcol, 
			     gretl_matrix_get(sigma, i, j)); 
	    kcol += (models[j])->ncoeff;
	}
	krow += (models[i])->ncoeff;
    }

    if (err) goto bailout;

    gretl_matrix_free(Xj);
    Xj = NULL;
    gretl_matrix_free(M);
    M = NULL;

    tmp_y = malloc(mk * sizeof *tmp_y);
    if (tmp_y == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* form Y column vector (m x k) */
    v = 0;
    for (i=0; i<m; i++) { /* loop over the m vertically arranged
			     blocks in the final column vector */
	double xx;
	const double *y;

	make_sys_X_block(Xi, models[i], (const double **) *pZ, 
			 pdinfo->t1, systype);
	for (j=0; j<(models[i])->ncoeff; j++) { /* loop over the rows within each of 
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

    calculate_sys_coefficients(models, (const double **) *pZ, X, uhat, 
			       tmp_y, m, mk);

    gls_sigma_from_uhat(sigma, uhat, m, T);

    /* FIXME: implement FIML here, after getting 3SLS starting values */
    if (systype == FIML) {
	pprintf(prn, "Sorry, FIML is not implmented yet.\n"
		"The following are just 3SLS estimates.\n");
    }

    j = 0;
    if (system_save_uhat(sys)) {
	j = pdinfo->v;
	err = dataset_add_vars(m, pZ, pdinfo);
    }
    if (system_save_yhat(sys)) {
	if (j == 0) j = pdinfo->v;
	err = dataset_add_vars(m, pZ, pdinfo);
    }

    for (i=0; i<m; i++) {
	printmodel(models[i], pdinfo, prn);
	add_results_to_dataset(sys, models[i], i, &j, *pZ, pdinfo, systype);
	if (systype == THREESLS || systype == FIML) {
	    tsls_free_data(models[i]);
	}
    }

    if (!err) {
	pprintf(prn, "%s\n(%s)\n\n",
		_("Cross-equation VCV for residuals"),
		_("correlations above the diagonal"));
	print_system_vcv(sigma, 1, prn);
    }

 bailout:

    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);
    gretl_matrix_free(M);

    gretl_matrix_free(X);
    gretl_matrix_free(sigma);
    gretl_matrix_free(uhat);

    free(tmp_y);

    for (i=0; i<m; i++) {
	free_model(models[i]);
    }
    free(models);

    restore_sample(pdinfo, orig_t1, orig_t2);

    return err;
}
