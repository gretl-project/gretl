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

#include <f2c.h> 
#include <cblas.h>
#include <clapack_double.h>

enum {
    GRETL_MATRIX_OK = 0,
    GRETL_MATRIX_NOMEM,
    GRETL_MATRIX_NON_CONFORM,
    GRETL_MATRIX_RANGE
} gretl_matrix_errors;

enum {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE = 1
} gretl_matrix_mods;

typedef struct _gretl_matrix gretl_matrix;
typedef struct _gretl_matrix gretl_vector;

struct _gretl_matrix {
    int rows;
    int cols;
    double *val;
};

/* Note: gretl matrices are stored in column major order, for 
   compatibility with the Fortran-style matrix storage used by
   Lapack.
*/

#define mdx(m,i,j) (j)*(m)->rows+(i)
#define mdxtr(m,i,j) (i)*(m)->cols+(j)

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))
#define gretl_vector_free(v) gretl_matrix_free(v)
#define gretl_vector_get(v,i) gretl_matrix_get((v),0,(i))
#define gretl_vector_set(v,i,x) gretl_matrix_set((v),0,(i),(x))


static const char *wspace_fail = "Workspace query failed\n";

static gretl_matrix *gretl_matrix_alloc (int rows, int cols)
{
    gretl_matrix *m;

    m = malloc(sizeof *m);
    if (m == NULL) return m;

    m->val = malloc(rows * cols * sizeof *m->val);

    if (m->val == NULL) {
	free(m);
	return NULL;
    }

    m->rows = rows;
    m->cols = cols;

    return m;
}

static void gretl_matrix_free (gretl_matrix *m)
{
    if (m == NULL) return;

    free(m->val);
    free(m);
}

static double gretl_matrix_get (const gretl_matrix *m, int i, int j)
{
    if (m == NULL || m->val == NULL) return -999.0;

    if (i >= m->rows || j >= m->cols) return -999.0;

    return m->val[mdx(m, i, j)];
}

static int gretl_matrix_set (gretl_matrix *m, int i, int j, double x)
{
    if (m == NULL || m->val == NULL) return 1;

    if (i >= m->rows || j >= m->cols) return 1;

    m->val[mdx(m, i, j)] = x;

    return 0;
}

static gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
						 int rows, int cols)
{
    int i, j, p;
    gretl_matrix *m;

    m = gretl_matrix_alloc(rows, cols);
    if (m == NULL) return m;

    p = 0;
    for (j=0; j<rows; j++) {
	for (i=0; i<cols; i++) {
	    m->val[p++] = X[i][j];
	}
    } 

    return m;
}

static int gretl_matmult_mod (const gretl_matrix *a, int aflag,
			      const gretl_matrix *b, int bflag,
			      gretl_matrix *c)
{
    int i, j, k;
    int lrows, lcols;
    int rrows, rcols;
    double x, y;

    lrows = (aflag == GRETL_MOD_NONE)? a->rows : a->cols;
    lcols = (aflag == GRETL_MOD_NONE)? a->cols : a->rows;
    rrows = (bflag == GRETL_MOD_NONE)? b->rows : b->cols;
    rcols = (bflag == GRETL_MOD_NONE)? b->cols : b->rows;

    if (lcols != rrows) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (c->cols != lcols || c->rows != rrows) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<rrows; i++) {
	for (j=0; j<lcols; j++) {
	    c->val[mdx(c, i, j)] = 0.0;
	    for (k=0; k<lcols; k++) {
		x = (aflag == GRETL_MOD_NONE)? 
		    a->val[mdx(a, i, k)] : a->val[mdxtr(a, i, k)];
		y = (bflag == GRETL_MOD_NONE)? 
		    b->val[mdx(b, k, j)] : b->val[mdxtr(b, k, j)];
		c->val[mdx(c, i, j)] += x * y;
	    }
	}
    }

    return GRETL_MATRIX_OK;
}

static int gretl_matmult (const gretl_matrix *a, const gretl_matrix *b,
			  gretl_matrix *c)
{
    return gretl_matmult_mod(a, GRETL_MOD_NONE,
			     b, GRETL_MOD_NONE,
			     c);
}

static void gretl_matrix_print (gretl_matrix *X, int rows, int cols,
				int triangle, PRN *prn)
{
    int i, j, jmax;
    double x;
    char numstr[16];

    jmax = (triangle)? 1 : cols;

    for (i=0; i<rows; i++) {
	for (j=0; j<jmax; j++) {
	    pprintf(prn, "%#10.5g ", gretl_matrix_get(X, i, j));
	}
	for (j=jmax; j<cols; j++) {
	    x = gretl_matrix_get(X, i, i) * gretl_matrix_get(X, j, j);
	    x = sqrt(x);
	    x = gretl_matrix_get(X, i, j) / x;
	    sprintf(numstr,"(%.3f)", x); 
	    pprintf(prn, "%11s", numstr);
	}
	pputs(prn, "\n");
	if (triangle && jmax < cols) jmax++;
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
    gretl_matrix *sigma, *inverse;

    sigma = gretl_matrix_alloc (m, m);
    inverse = gretl_matrix_alloc (m, m);

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

#if 0 /* FIXME */
    gretl_LU_decomp (sigma, p, &sign);
    gretl_LU_invert (sigma, p, inverse);
#endif

    gretl_matrix_free (sigma);

    return inverse;
}

/* m = number of equations 
   k = number of indep vars per equation 
*/

static int ApB (const gretl_matrix * A, const gretl_matrix * B,
		gretl_matrix *C) 
{
    int err;

    err = gretl_matmult_mod (A, GRETL_MOD_TRANSPOSE,
			     B, GRETL_MOD_NONE,
			     C);
    return err;
}

static void sur_resids (MODEL *pmod, double **Z, gretl_matrix *uhat)
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
	gretl_matrix_set(uhat, pmod->ID, t, pmod->uhat[t]);
    }

    pmod->ess = 0.0;
    for (t=0; t<T; t++) {
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }
    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    /* pmod->rsq = 1.0 - (pmod->ess / pmod->tss); hmm... */
}

static int calculate_coefficients (MODEL **models, double **Z,
				   gretl_matrix *X, gretl_matrix *uhat,
				   double *tmp_y, int m, int k)
{
    gretl_vector *y;
    gretl_vector *coeff;
    gretl_matrix *vcv;
    int i, j;
    int ncoeff = m * k;

    y = gretl_vector_alloc (ncoeff);
    coeff = gretl_vector_alloc (ncoeff);
    vcv = gretl_matrix_alloc (ncoeff, ncoeff);

    for (i=0; i<ncoeff; i++) {
	gretl_vector_set(y, i, tmp_y[i]);
	gretl_vector_set(coeff, i, 0.0);
    }

#if 0 /* FIXME */
    gretl_LU_decomp (X, p, &sign); 
    gretl_LU_solve (X, p, y, coeff);
    gretl_LU_invert (X, p, vcv);
#endif

    for (i=0; i<m; i++) {
	for (j=0; j<k; j++) {
	    (models[i])->coeff[j+1] = gretl_vector_get(coeff, i * k + j);
	    (models[i])->sderr[j+1] = 
		sqrt(gretl_matrix_get(vcv, i * k + j, i * k + j));
	}
	sur_resids(models[i], Z, uhat);
    }

    gretl_vector_free (y);
    gretl_vector_free (coeff);
    gretl_matrix_free (vcv);

    return 0;
}

int sur (gretl_equation_system *sys, double ***pZ,
	 DATAINFO *pdinfo, PRN *prn)
{
    int i, j, k, m, T, t, l;
    gretl_matrix *X, *Xi, *Xj, *M;
    gretl_matrix *uhat, *sigma;
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

    X = gretl_matrix_alloc (bigrows, bigrows);
    Xi = gretl_matrix_alloc (T, k);
    Xj = gretl_matrix_alloc (T, k);
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
	    ApB ((const gretl_matrix *) Xi, 
		 (i == j)? (const gretl_matrix *) Xi : 
		 (const gretl_matrix *) Xj, M);
	    kronecker_place (X, (const gretl_matrix *) M,
			     i, j, k, 
			     gretl_matrix_get(sigma, i, j)); 
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

    calculate_coefficients (models, *pZ, X, uhat, tmp_y, m, k);
    gls_sigma_from_uhat (sigma, uhat, m, T);

    for (i=0; i<m; i++) {
	printmodel(models[i], pdinfo, prn);
	free_model(models[i]);
    }

    pputs(prn, "Cross-equation VCV for residuals\n"
	  "(correlations above the diagonal)\n\n");

    gretl_matrix_print(sigma, m, m, 1, prn);

    gretl_matrix_free(X);
    gretl_matrix_free(Xi);
    gretl_matrix_free(Xj);
    gretl_matrix_free(M);
    gretl_matrix_free(sigma);
    gretl_matrix_free(uhat);

    free(tmp_y);
    free(models);

    return 0;
}

static int lapack_invert_gretl_matrix (gretl_matrix *m)
{
    integer n = m->rows;    
    integer info;
    integer lwork;
    integer *ipiv;

    double *work;

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) {
	return 1;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }    

    dgetrf_(&n, &n, m->val, &n, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	return info;
    }

    lwork = -1;
    dgetri_(&n, m->val, &n, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	free(ipiv);
	return 1;
    }

    lwork = (integer) work[0];

#ifdef LAPACK_DEBUG
    printf("dgetri: workspace = %d\n", (int) lwork);
#endif

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }  

    dgetri_(&n, m->val, &n, ipiv, work, &lwork, &info);

#ifdef LAPACK_DEBUG
    printf("dgetri: info = %d\n", (int) info);
#endif

    return info;
}

static double *gretl_lapack_eigenvals (gretl_matrix *m) 
{
    integer n = m->rows;
    integer info, sdim;
    integer lwork;
    integer one = 1;

    double *work;
    double *wr, *wi;

    char job = 'N', sort = 'N';

    work = malloc(sizeof *work);
    if (work == NULL) {
	return NULL;
    }

    wr = malloc(n * sizeof *wr);
    if (wr == NULL) {
	free(work);
	return NULL;
    }

    wi = malloc(n * sizeof *wi);
    if (wi == NULL) {
	free(work);
	free(wr);
	return NULL;
    }

    lwork = -1; /* find optimal workspace size */
    dgees_(&job, &sort, NULL, &n, 
	   m->val, &n, &sdim, wr, wi, NULL, &one, 
	   work, &lwork, NULL, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	free(work);
	free(wr);
	free(wi);
	return NULL;
    }	

    lwork = (integer) work[0];

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(wr);
	free(wi);
	return NULL;
    }    

    dgees_(&job, &sort, NULL, &n, 
	   m->val, &n, &sdim, wr, wi, NULL, &one, 
	   work, &lwork, NULL, &info);

    if (info != 0) {
	free(wr);
	wr = NULL;
    }

    free(wi);
    free(work);

    return wr;
}

#if 0
static void simple_print_gretl_matrix (gretl_matrix *m)
{
    int i, j;

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	   printf("%#13.7g ", gretl_matrix_get(m, i, j));
	}
	putchar('\n');
    } 
    putchar('\n');
}
#endif

int johansen_eigenvals (const double **X, const double **Y, const double **Z, 
			int k, double *evals)
{
    gretl_matrix *Suu, *Svv, *Suv;
    gretl_matrix *Inv, *TmpL, *TmpR, *M;
    double *eigvals;
    int err = 0;

    Suu = gretl_matrix_from_2d_array(X, k, k);
    Svv = gretl_matrix_from_2d_array(Y, k, k);
    Suv = gretl_matrix_from_2d_array(Z, k, k);

    Inv = gretl_matrix_alloc(k, k);
    TmpL = gretl_matrix_alloc(k, k);
    TmpR = gretl_matrix_alloc(k, k);
    M = gretl_matrix_alloc(k, k);

    /* calculate Suu^{-1} Suv */
    lapack_invert_gretl_matrix(Suu);
    gretl_matmult(Suu, Suv, TmpR);

    /* calculate Svv^{-1} Suv' */
    lapack_invert_gretl_matrix(Svv);
    gretl_matmult_mod(Svv, GRETL_MOD_NONE,
		      Suv, GRETL_MOD_TRANSPOSE, 
		      TmpL);

    gretl_matmult(TmpL, TmpR, M);
    eigvals = gretl_lapack_eigenvals(M);

    if (eigvals != NULL) {
	int i;

	for (i=0; i<k; i++) {
	    evals[i] = eigvals[i];
	}
	free(eigvals);
    }

    /* free stuff */
    gretl_matrix_free(Svv);
    gretl_matrix_free(Suu);
    gretl_matrix_free(Suv);

    gretl_matrix_free(Inv);
    gretl_matrix_free(TmpL);
    gretl_matrix_free(TmpR);
    gretl_matrix_free(M);

    return err;
}
