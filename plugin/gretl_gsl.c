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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgretl.h"
#include <gsl/gsl_linalg.h>

static void gretl_gsl_matrix_print (gsl_matrix *X, int rows, int cols,
				    int triangle, PRN *prn)
{
    int i, j, jmax;
    double x;
    char numstr[16];

    jmax = (triangle)? 1 : cols;

    for (i=0; i<rows; i++) {
	for (j=0; j<jmax; j++) {
	    pprintf(prn, "%#10.5g ", gsl_matrix_get(X, i, j));
	}
	for (j=jmax; j<cols; j++) {
	    x = gsl_matrix_get(X, i, i) * gsl_matrix_get(X, j, j);
	    x = sqrt(x);
	    x = gsl_matrix_get(X, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    pprintf(prn, "%11s", numstr);
	}
	pputs(prn, "\n");
	if (triangle && jmax < cols) jmax++;
    }
}

static void kronecker_place (gsl_matrix *X, 
			     const gsl_matrix *M,
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
	    x = gsl_matrix_get(M, i, j);
	    gsl_matrix_set(X, row, col, x * scale);
	}
    }
}

static void make_Xi_from_Z (gsl_matrix *X, double **Z, int *list, int T)
{
    int i, t;

    for (i=2; i<=list[0]; i++) {
	for (t=0; t<T; t++) {
	    gsl_matrix_set(X, t, i-2, Z[list[i]][t]);
	}
    }
}

static int
gls_sigma_from_uhat (gsl_matrix *sigma, const gsl_matrix *e, int m, int T)
{
    int i, j, t;
    double xx;

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gsl_matrix_get(e, i, t) * gsl_matrix_get(e, j, t);
	    }
	    gsl_matrix_set (sigma, i, j, xx / T);
	}
    }

    return 0;
}

static gsl_matrix *
gls_sigma_inverse_from_uhat (const gsl_matrix *e, int m, int T)
{
    int i, j, t;
    double xx;
    gsl_matrix *sigma, *inverse;
    gsl_permutation *p;
    int sign;

    sigma = gsl_matrix_alloc (m, m);
    inverse = gsl_matrix_alloc (m, m);
    p = gsl_permutation_alloc (m);  

    /* construct sigma: s_{ij} = e'_i * e_j / T  */
    for (i=0; i<m; i++) {
	for (j=0; j<m; j++) {
	    xx = 0.0;
	    for (t=0; t<T; t++) {
		xx += gsl_matrix_get(e, i, t) * gsl_matrix_get(e, j, t);
	    }
	    gsl_matrix_set (sigma, i, j, xx / T);
	}
    }

    gsl_linalg_LU_decomp (sigma, p, &sign);
    gsl_linalg_LU_invert (sigma, p, inverse);

    gsl_matrix_free (sigma);
    gsl_permutation_free (p); 

    return inverse;
}

/* m = number of equations 
   k = number of indep vars per equation 
*/

static int ApB (const gsl_matrix * A, const gsl_matrix * B,
		gsl_matrix *C) 
{
    int err;

    err = gsl_linalg_matmult_mod (A, GSL_LINALG_MOD_TRANSPOSE,
				  B, GSL_LINALG_MOD_NONE,
				  C);
    return err;
}

static void sur_resids (MODEL *pmod, double **Z, gsl_matrix *uhat)
{
    int i, t;
    int k = pmod->ncoeff, T = pmod->nobs;
    double fit;

    for (t=0; t<T; t++) {
	fit = 0.0;
	for (i=0; i<k; i++) {
	    fit += pmod->coeff[i+1] * Z[pmod->list[i+2]][t];
	}
	pmod->yhat[t] = fit;
	pmod->uhat[t] = Z[pmod->list[1]][t] - fit;
	/* for cross-equation vcv */
	gsl_matrix_set(uhat, pmod->ID, t, pmod->uhat[t]);
    }

    pmod->ess = 0.0;
    for (t=0; t<T; t++) {
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }
    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    /* pmod->rsq = 1.0 - (pmod->ess / pmod->tss); hmm... */
}

static int calculate_coefficients (MODEL **models, double **Z,
				   gsl_matrix *X, gsl_matrix *uhat,
				   double *tmp_y, int m, int k)
{
    gsl_vector *y;
    gsl_vector *coeff;
    gsl_permutation *p; 
    gsl_matrix *vcv;
    int i, j, sign;
    int ncoeff = m * k;

    y = gsl_vector_alloc (ncoeff);
    coeff = gsl_vector_calloc (ncoeff);
    p = gsl_permutation_alloc (ncoeff);
    vcv = gsl_matrix_alloc (ncoeff, ncoeff);

    for (i=0; i<ncoeff; i++) {
	gsl_vector_set(y, i, tmp_y[i]);
    }

    gsl_linalg_LU_decomp (X, p, &sign); 
    gsl_linalg_LU_solve (X, p, y, coeff);
    gsl_linalg_LU_invert (X, p, vcv);

    for (i=0; i<m; i++) {
	for (j=0; j<k; j++) {
	    (models[i])->coeff[j+1] = gsl_vector_get(coeff, i * k + j);
	    (models[i])->sderr[j+1] = 
		sqrt(gsl_matrix_get(vcv, i * k + j, i * k + j));
	}
	sur_resids(models[i], Z, uhat);
    }

    gsl_vector_free (y);
    gsl_vector_free (coeff);
    gsl_permutation_free (p);
    gsl_matrix_free (vcv);

    return 0;
}

int sur (gretl_equation_system *sys, double ***pZ,
	 DATAINFO *pdinfo, PRN *prn)
{
    int i, j, k, m, T, t, l;
    gsl_matrix *X, *Xi, *Xj, *M;
    gsl_matrix *uhat, *sigma;
    double *tmp_y, *y;
    int v, bigrows;
    MODEL **models;

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

    X = gsl_matrix_alloc (bigrows, bigrows);
    Xi = gsl_matrix_alloc (T, k);
    Xj = gsl_matrix_alloc (T, k);
    M = gsl_matrix_alloc(k, k);
    uhat = gsl_matrix_alloc(m, T);

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
	    gsl_matrix_set(uhat, i, t, (models[i])->uhat[t]);
	}
    }

    sigma = gls_sigma_inverse_from_uhat (uhat, m, T);

    /* Xi = data matrix for equation i, specified in lists[i] */
    for (i=0; i<m; i++) {
#ifdef SUR_DEBUG
	fprintf(stderr, "doing make_Xi_from_Z(), i=%d\n", i);
#endif
	make_Xi_from_Z(Xi, *pZ, sys->lists[i], T);
	for (j=0; j<m; j++) { 
	    if (i != j) {
		make_Xi_from_Z(Xj, *pZ, sys->lists[j], T);
	    }
	    ApB ((const gsl_matrix *) Xi, 
		 (i == j)? (const gsl_matrix *) Xi : 
		 (const gsl_matrix *) Xj, M);
	    kronecker_place (X, (const gsl_matrix *) M,
			     i, j, k, 
			     gsl_matrix_get(sigma, i, j)); 
	}
    }

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
		    xx += gsl_matrix_get(Xi, t, j) * y[t];
		}
		xx *= gsl_matrix_get(sigma, i, l);
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

    calculate_coefficients (models, *pZ, X, uhat, tmp_y, m, k);
    gls_sigma_from_uhat (sigma, uhat, m, T);

    for (i=0; i<m; i++) {
	printmodel(models[i], pdinfo, prn);
	free_model(models[i]);
    }

    pputs(prn, "Cross-equation VCV for residuals\n"
	  "(correlations above the diagonal)\n\n");

    gretl_gsl_matrix_print(sigma, m, m, 1, prn);

    gsl_matrix_free(X);
    gsl_matrix_free(Xi);
    gsl_matrix_free(Xj);
    gsl_matrix_free(M);
    gsl_matrix_free(sigma);
    gsl_matrix_free(uhat);

    free(tmp_y);
    free(models);

    return 0;
}


