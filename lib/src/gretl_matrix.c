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

static const char *wspace_fail = "gretl_matrix: workspace query failed\n";

/* ....................................................... */

gretl_matrix *gretl_matrix_alloc (int rows, int cols)
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

/* ....................................................... */

gretl_matrix *gretl_matrix_copy (gretl_matrix *m)
{
    gretl_matrix *c;
    int i, n = m->rows * m->cols;

    c = malloc(sizeof *c);
    if (c == NULL) return c;

    c->val = malloc(n * sizeof *c->val);

    if (c->val == NULL) {
	free(c);
	return NULL;
    }

    c->rows = m->rows;
    c->cols = m->cols;

    for (i=0; i<n; i++) {
	c->val[i] = m->val[i];
    }

    return c;
}

/* ....................................................... */

void gretl_matrix_free (gretl_matrix *m)
{
    if (m == NULL) return;

    if (m->val != NULL) free(m->val);
    free(m);
}

/* ....................................................... */

double *gretl_matrix_steal_data (gretl_matrix *m)
{
    double *vals = NULL;

    if (m != NULL) {
	vals = m->val;
	m->val = NULL;
    }
    return vals;
}

/* ....................................................... */

double gretl_matrix_get (const gretl_matrix *m, int i, int j)
{
    if (m == NULL || m->val == NULL) return -999.0;

    if (i >= m->rows || j >= m->cols) return -999.0;

    return m->val[mdx(m, i, j)];
}

/* ....................................................... */

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x)
{
    if (m == NULL || m->val == NULL) return 1;

    if (i >= m->rows || j >= m->cols) return 1;

    m->val[mdx(m, i, j)] = x;

    return 0;
}

#ifdef LDEBUG
static void simple_matrix_print (gretl_matrix *X, int rows, int cols,
				 PRN *prn)
{
    int i, j;

    pprintf(prn, "printing %d x %d matrix...\n", rows, cols);

    for (i=0; i<rows; i++) {
	for (j=0; j<cols; j++) {
	    pprintf(prn, "%#10.5g ", gretl_matrix_get(X, i, j));
	}
	pputs(prn, "\n");
    }
    pputs(prn, "\n");
}
#endif

int gretl_LU_solve (gretl_matrix *a, gretl_vector *b)
{
    /* Solves ax = b.  On exit, b is replaced by the solution vector */
    char trans = 'N';
    integer info;
    integer m = a->rows;
    integer n = a->cols;
    integer nrhs = 1;
    integer ldb = gretl_vector_get_length(b);
    integer *ipiv;

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) return 1;

    dgetrf_(&m, &n, a->val, &n, ipiv, &info);

    if (info != 0) {
	free(ipiv);
	return info;
    }

    dgetrs_(&trans, &n, &nrhs, a->val, &n, ipiv, b->val, &ldb, &info);

    free(ipiv);

    return info;
}

gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
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

int gretl_matmult_mod (const gretl_matrix *a, int aflag,
		       const gretl_matrix *b, int bflag,
		       gretl_matrix *c)
{
    int i, j, k;
    double x, y;
    int lrows, lcols;
    int rrows, rcols;
    int atr = (aflag == GRETL_MOD_TRANSPOSE);
    int btr = (bflag == GRETL_MOD_TRANSPOSE);
    int bmax = b->rows * b->cols;

    lrows = (atr)? a->cols : a->rows;
    lcols = (atr)? a->rows : a->cols;
    rrows = (btr)? b->cols : b->rows;
    rcols = (btr)? b->rows : b->cols;

    if (lcols != rrows) {
	fprintf(stderr, "gretl_matmult_mod: matrices not conformable\n");
	fprintf(stderr, "left-hand cols = %d, right-hand rows = %d\n",
		lcols, rrows);	
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (c->rows != lrows || c->cols != rcols) {
	fprintf(stderr, "gretl_matmult_mod: matrices not conformable\n");
	fprintf(stderr, "Product cols = %d, left-hand cols = %d;\n"
		"Product rows = %d, right-hand rows = %d\n",
		c->cols, lcols, c->rows, rrows);
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<lrows; i++) {
	for (j=0; j<rcols; j++) {
	    c->val[mdx(c, i, j)] = 0.0;
	    for (k=0; k<lcols; k++) {
		x = (atr)? a->val[mdxtr(a,i,k)] : a->val[mdx(a,i,k)];
		y = (btr)? b->val[mdxtr(b,k,j)] : b->val[mdx(b,k,j)];
		if (mdx(b,k,j) >= bmax) {
		    fprintf(stderr, "gretl_matmult_mod: Bmax = %d exceeded\n", 
			    bmax);
		    return 1;
		}
		c->val[mdx(c, i, j)] += x * y;
	    }
	}
    }

    return GRETL_MATRIX_OK;
}

int gretl_matmult (const gretl_matrix *a, const gretl_matrix *b,
		   gretl_matrix *c)
{
    return gretl_matmult_mod(a, GRETL_MOD_NONE,
			     b, GRETL_MOD_NONE,
			     c);
}

int gretl_invert_general_matrix (gretl_matrix *a)
{
    integer m = a->rows;
    integer n = a->cols;
    integer info;
    integer lwork;
    integer *ipiv;
    int lipiv;

    double *work;

    if (m <= n) lipiv = m;
    else lipiv = n;

    ipiv = malloc(lipiv * sizeof *ipiv);
    if (ipiv == NULL) {
	return 1;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return 1;
    }    

    dgetrf_(&m, &n, a->val, &m, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	return info;
    }

    lwork = -1;
    dgetri_(&n, a->val, &n, ipiv, work, &lwork, &info);

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

    dgetri_(&n, a->val, &n, ipiv, work, &lwork, &info);

#ifdef LAPACK_DEBUG
    printf("dgetri: info = %d\n", (int) info);
#endif

    free(work);
    free(ipiv);

    return info;
}

double *gretl_general_matrix_eigenvals (gretl_matrix *m) 
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

