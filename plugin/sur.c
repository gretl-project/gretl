/*
 *  Copyright (c) by Allin Cottrell
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

static void 
print_sur_vcv (const gretl_matrix *m, int triangle, PRN *prn)
{
    int i, j, jmax;
    int mcols = gretl_matrix_cols(m);
    int mrows = gretl_matrix_rows(m);
    double x;
    char numstr[16];

    jmax = (triangle)? 1 : mcols;

    for (i=0; i<mrows; i++) {
	for (j=0; j<jmax; j++) {
	    pprintf(prn, "%#10.5g ", gretl_matrix_get(m, i, j));
	}
	for (j=jmax; j<mcols; j++) {
	    x = gretl_matrix_get(m, i, i) * gretl_matrix_get(m, j, j);
	    x = sqrt(x);
	    x = gretl_matrix_get(m, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    pprintf(prn, "%11s", numstr);
	}
	pputs(prn, "\n");
	if (triangle && jmax < mcols) jmax++;
    }
}

static void kronecker_place (gretl_matrix *X, 
			     const gretl_matrix *M,
			     int startrow, int startcol,
			     int k, double scale)
{
    int i, j;
    int row, col;
    double x;
    
    for (i=0; i<k; i++) {
	row = startrow * k + i;
	for (j=0; j<k; j++) {
	    col = startcol * k + j;
	    x = gretl_matrix_get(M, i, j);
	    gretl_matrix_set(X, row, col, x * scale);
	}
    }
}

static void make_Xi_from_Z (gretl_matrix *X, double **Z, int *list, int T)
{
    int i, t;

    for (i=2; i<=list[0]; i++) {
	for (t=0; t<T; t++) {
	    gretl_matrix_set(X, t, i-2, Z[list[i]][t]);
	}
    }
}

static int
gls_sigma_from_uhat (gretl_matrix *sigma, const gretl_matrix *e, int m, int T)
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
	    gretl_matrix_set (sigma, i, j, xx / T);
	}
    }

    return 0;
}

static gretl_matrix *
gls_sigma_inverse_from_uhat (const gretl_matrix *e, int m, int T)
{
    int i, j, t;
    double xx;
    gretl_matrix *sigma;

    sigma = gretl_matrix_alloc (m, m);

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gretl_matrix_get(e, i, t) * gretl_matrix_get(e, j, t);
	    }
	    gretl_matrix_set (sigma, i, j, xx / T);
	}
    }

    gretl_invert_general_matrix(sigma);

    return sigma;
}

/* m = number of equations 
   k = number of indep vars per equation 
*/

static void sur_resids (MODEL *pmod, double **Z, gretl_matrix *uhat)
{
    int i, t;
    int k = pmod->ncoeff, T = pmod->nobs;
    double fit;

    for (t=0; t<T; t++) {
	fit = 0.0;
	for (i=0; i<k; i++) {
	    fit += pmod->coeff[i] * Z[pmod->list[i+2]][t];
	}
	pmod->yhat[t] = fit;
	pmod->uhat[t] = Z[pmod->list[1]][t] - fit;
	/* for cross-equation vcv */
	gretl_matrix_set(uhat, pmod->ID, t, pmod->uhat[t]);
    }

    pmod->ess = 0.0;
    for (t=0; t<T; t++) {
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }
    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    /* pmod->rsq = 1.0 - (pmod->ess / pmod->tss); hmm... */
}

static int 
calculate_sur_coefficients (MODEL **models, double **Z,
			    gretl_matrix *X, gretl_matrix *uhat,
			    double *tmp_y, int m, int k)
{
    gretl_vector *coeff;
    gretl_matrix *vcv;
    int i, j;
    int ncoeff = m * k;

    coeff = gretl_vector_alloc(ncoeff);
    if (coeff == NULL) return 1;

    for (i=0; i<ncoeff; i++) {
	gretl_vector_set(coeff, i, tmp_y[i]);
    }

    vcv = gretl_matrix_copy(X);
    gretl_LU_solve (X, coeff);
    gretl_invert_general_matrix(vcv); 

    for (i=0; i<m; i++) {
	for (j=0; j<k; j++) {
	    (models[i])->coeff[j] = gretl_vector_get(coeff, i * k + j);
	    (models[i])->sderr[j] = 
		sqrt(gretl_matrix_get(vcv, i * k + j, i * k + j));
	}
	sur_resids(models[i], Z, uhat);
    }

    gretl_vector_free (coeff);
    gretl_matrix_free (vcv);

    return 0;
}

static void add_results_to_dataset (gretl_equation_system *sys, 
				    MODEL *pmod, int i, int *pj,
				    double **Z, DATAINFO *pdinfo)
{
    int t;

    if (sys->flags & GRETL_SYSTEM_SAVE_UHAT) {
	for (t=0; t<pdinfo->n; t++) {
	    Z[*pj][t] = pmod->uhat[t];
	}
	sprintf(pdinfo->varname[*pj], "uhat_s%02d", i + 1);
	sprintf(VARLABEL(pdinfo, *pj), _("SUR residual, equation %d"), 
		i + 1);
	*pj += 1;
    }
    if (sys->flags & GRETL_SYSTEM_SAVE_YHAT) {
	for (t=0; t<pdinfo->n; t++) {
	    Z[*pj][t] = pmod->yhat[t];
	}	
	sprintf(pdinfo->varname[*pj], "yhat_s%02d", i + 1);
	sprintf(VARLABEL(pdinfo, *pj), _("SUR fitted value, equation %d"), 
		i + 1);
	*pj += 1;
    }
}

int sur (gretl_equation_system *sys, 
	 double ***pZ, DATAINFO *pdinfo, 
	 PRN *prn)
{
    int i, j, k, m, T, t, l;
    gretl_matrix *X, *Xi, *Xj, *M;
    gretl_matrix *uhat, *sigma;
    double *tmp_y, *y;
    int v, bigrows;
    MODEL **models;
    int err = 0;

    /* number of equations */
    m = sys->n_equations;

    /* number of indep vars per equation */
    k = sys->lists[0][0] - 1;

    /* number of observations per series */
    T = pdinfo->t2 - pdinfo->t1 + 1;

    bigrows = m * k;

    models = malloc(m * sizeof *models);
    if (models == NULL) return E_ALLOC;

    for (i=0; i<m; i++) {
	models[i] = gretl_model_new(pdinfo);
	if (models[i] == NULL) return E_ALLOC;
    }

    X = gretl_matrix_alloc(bigrows, bigrows);
    Xi = gretl_matrix_alloc(T, k);
    Xj = gretl_matrix_alloc(T, k);
    M = gretl_matrix_alloc(k, k);
    uhat = gretl_matrix_alloc(m, T);

    /* first grab the OLS residuals */
    for (i=0; i<m; i++) {
	*models[i] = lsq(sys->lists[i], pZ, pdinfo, OLS, 1, 0.0);
	if ((models[i])->errcode) {
	    fprintf(stderr, "model failed on lists[%d], code=%d\n",
		    i, (models[i])->errcode);
	    return 1;
	}
	(models[i])->ID = i;
	(models[i])->aux = AUX_SUR;
	for (t=0; t<T; t++) {
	    gretl_matrix_set(uhat, i, t, (models[i])->uhat[t]);
	}
    }

    sigma = gls_sigma_inverse_from_uhat (uhat, m, T);

#ifdef LDEBUG 
    pprintf(prn, "gls sigma inverse matrix\n");
    simple_matrix_print(sigma, m, m, prn);
    /* OK so far, it seems */
#endif

    /* Xi = data matrix for equation i, specified in lists[i] */
    for (i=0; i<m; i++) {
	const gretl_matrix *Y;

#ifdef SUR_DEBUG
	fprintf(stderr, "doing make_Xi_from_Z(), i=%d\n", i);
#endif
	make_Xi_from_Z(Xi, *pZ, sys->lists[i], T);
	for (j=0; j<m; j++) { 
	    if (i != j) {
		make_Xi_from_Z(Xj, *pZ, sys->lists[j], T);
	    }
#ifdef LDEBUG
	    pprintf(prn, "Xi:\n");
	    simple_matrix_print(Xi, k, k, prn);	    
#endif
	    Y = (i == j)? Xi : Xj;
	    gretl_matrix_multiply_mod ((const gretl_matrix *) Xi, 
				       GRETL_MOD_TRANSPOSE,
				       Y, GRETL_MOD_NONE, M);
#ifdef LDEBUG
	    pprintf(prn, "M:\n");
	    simple_matrix_print(M, k, k, prn);	    
#endif
	    kronecker_place (X, (const gretl_matrix *) M,
			     i, j, k, 
			     gretl_matrix_get(sigma, i, j)); 
	}
    }

#ifdef LDEBUG 
    pprintf(prn, "big X matrix\n");
    simple_matrix_print(X, bigrows, bigrows, prn);
#endif

    tmp_y = malloc((m * k) * sizeof *tmp_y);

    /* form Y column vector (m x k) */
    v = 0;
    for (i=0; i<m; i++) { /* loop over the m vertically arranged
			     blocks in the final column vector */
	double xx;

#ifdef SUR_DEBUG
	fprintf(stderr, "working on block %d\n", i);
#endif
	make_Xi_from_Z(Xi, *pZ, sys->lists[i], T);
	for (j=0; j<k; j++) { /* loop over the k rows within each of 
				 the m blocks */
#ifdef SUR_DEBUG
	    fprintf(stderr, " working on row %d\n", i * k + j);
#endif
	    tmp_y[v] = 0.0;
	    for (l=0; l<m; l++) { /* loop over the m components that
				     must be added to form each element */
#ifdef SUR_DEBUG
		fprintf(stderr, "  component %d of row %d\n", 
		       l+1, i * k + j + 1);
		fprintf(stderr, "    sigma(%d, %d) * ", i, l);
		fprintf(stderr, "X'_%d[%d] * ", i, j);
		fprintf(stderr, "y_%d\n", l);
#endif
		y = (*pZ)[sys->lists[l][1]];
		/* multiply X'[l] into y */
		xx = 0.0;
		for (t=0; t<T; t++) {
		    xx += gretl_matrix_get(Xi, t, j) * y[t];
		}
		xx *= gretl_matrix_get(sigma, i, l);
		tmp_y[v] += xx;
	    }
#ifdef SUR_DEBUG
	    fprintf(stderr, " finished row %d\n", i * k + j);
#endif
	    v++;
	}
#ifdef SUR_DEBUG
	fprintf(stderr, "finished block %d\n", i);
#endif	
    }

    calculate_sur_coefficients (models, *pZ, X, uhat, tmp_y, m, k);
    gls_sigma_from_uhat (sigma, uhat, m, T);

    j = 0;
    if (sys->flags & GRETL_SYSTEM_SAVE_UHAT) {
	j = pdinfo->v;
	err = dataset_add_vars(m, pZ, pdinfo);
    }
    if (sys->flags & GRETL_SYSTEM_SAVE_YHAT) {
	if (j == 0) j = pdinfo->v;
	err = dataset_add_vars(m, pZ, pdinfo);
    }

    for (i=0; i<m; i++) {
	printmodel(models[i], pdinfo, prn);
	add_results_to_dataset(sys, models[i], i, &j, *pZ, pdinfo);
	free_model(models[i]);
    }

    pprintf(prn, "%s\n(%s)\n\n",
	    _("Cross-equation VCV for residuals"),
	    _("correlations above the diagonal"));

    print_sur_vcv(sigma, 1, prn);

    gretl_matrix_free(X);
    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);
    gretl_matrix_free(M);
    gretl_matrix_free(sigma);
    gretl_matrix_free(uhat);

    free(tmp_y);
    free(models);

    return err;
}
