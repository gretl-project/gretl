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
#include "gretl_matrix_private.h"

static const char *wspace_fail = "gretl_matrix: workspace query failed\n";

#define mdx(a,i,j)   ((j)*(a)->rows+(i))
#define mdxtr(a,i,j) ((i)*(a)->rows+(j))

static int packed_idx (int nrows, int i, int j);

/* ....................................................... */

static gretl_matrix *real_gretl_matrix_alloc (int rows, int cols,
					      int packed)
{
    gretl_matrix *m;

    m = malloc(sizeof *m);
    if (m == NULL) return m;

    if (packed) { /* symmetric, only triangle stored */
	int n = (rows * rows + rows) / 2;

	m->val = malloc(n * sizeof *m->val);
    } else {
	m->val = malloc(rows * cols * sizeof *m->val);
    }

    if (m->val == NULL) {
	free(m);
	return NULL;
    }

    m->rows = rows;
    m->cols = cols;
    m->packed = packed;
    m->t = 0;

    return m;
}

/* ....................................................... */

gretl_matrix *gretl_matrix_alloc (int rows, int cols)
{
    return real_gretl_matrix_alloc(rows, cols, 0);
}

/* ....................................................... */

gretl_matrix *gretl_packed_matrix_alloc (int rows)
{
    return real_gretl_matrix_alloc(rows, rows, 1);
}

/* ....................................................... */

gretl_matrix *gretl_diagonal_matrix (const double *d, int n, int mod)
{
    gretl_matrix *m;
    double x;
    int i, j;

    m = real_gretl_matrix_alloc(n, n, 0);
    if (m == NULL) return NULL;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    if (i == j) {
		x = *d++;
		if (mod == GRETL_MOD_SQUARE) {
		    gretl_matrix_set(m, i, j, x * x);
		} else {
		    gretl_matrix_set(m, i, j, x);
		}
	    }
	    else gretl_matrix_set(m, i, j, 0.0);
	}
    }

    return m;
}

/* ....................................................... */

static gretl_matrix *gretl_matrix_copy_mod (gretl_matrix *m, int mod)
{
    gretl_matrix *c;
    int i, j, n;

    if (m->packed) {
	n = (m->rows * m->rows + m->rows) / 2;
    } else {
	n = m->rows * m->cols;
    }

    c = malloc(sizeof *c);
    if (c == NULL) return c;

    c->val = malloc(n * sizeof *c->val);

    if (c->val == NULL) {
	free(c);
	return NULL;
    }

    c->rows = m->rows;
    c->cols = m->cols;

    c->packed = m->packed;

    if (mod == GRETL_MOD_TRANSPOSE) {
	for (i=0; i<m->rows; i++) {
	    for (j=0; j<m->cols; j++) {
		if (m->packed) { 
		    c->val[packed_idx(m->rows, i, j)] = 
			m->val[packed_idx(m->rows, j, i)];
		} else {
		    c->val[mdx(m, i, j)] = 
			m->val[mdx(m, j, i)];
		}
	    }
	}
    } else { /* not transposing */
	for (i=0; i<n; i++) {
	    c->val[i] = m->val[i];
	}
    }

    return c;
}

gretl_matrix *gretl_matrix_copy (gretl_matrix *m)
{
    return gretl_matrix_copy_mod(m, GRETL_MOD_NONE);
}

/* ....................................................... */

void gretl_matrix_free (gretl_matrix *m)
{
    if (m == NULL) return;

    if (m->val != NULL) free(m->val);
    free(m);
}

/* ....................................................... */

void gretl_matrix_zero (gretl_matrix *m)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    if (m->packed) {
	n = (m->rows * m->rows + m->rows) / 2;
    } else {
	n = m->rows * m->cols;
    }
    
    for (i=0; i<n; i++) m->val[i] = 0.0;
}

/* ....................................................... */

void gretl_matrix_multiply_by_scalar (gretl_matrix *m, double x)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    if (m->packed) {
	n = (m->rows * m->rows + m->rows) / 2;
    } else {
	n = m->rows * m->cols;
    }
    
    for (i=0; i<n; i++) m->val[i] *= x;
}

/* ....................................................... */

void gretl_matrix_divide_by_scalar (gretl_matrix *m, double x)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    if (m->packed) {
	n = (m->rows * m->rows + m->rows) / 2;
    } else {
	n = m->rows * m->cols;
    }
    
    for (i=0; i<n; i++) m->val[i] /= x;
}

/* ....................................................... */

int gretl_matrix_copy_values (gretl_matrix *targ, 
			      const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (targ->packed != src->packed) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (src->packed) {
	n = (src->rows * src->rows + src->rows) / 2;
    } else {
	n = src->rows * src->cols;
    }
    
    for (i=0; i<n; i++) targ->val[i] = src->val[i];

    return GRETL_MATRIX_OK;
}

/* ....................................................... */

int 
gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (targ->packed != src->packed) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (src->packed) {
	n = (src->rows * src->rows + src->rows) / 2;
    } else {
	n = src->rows * src->cols;
    }
    
    for (i=0; i<n; i++) targ->val[i] += src->val[i];

    return GRETL_MATRIX_OK;
}

/* On input, general matrix M; on output, the symmetric matrix
   S = M + M'
*/

int gretl_matrix_add_self_transpose (gretl_matrix *m)
{
    int i, j;
    double x1, x2;

    if (m->rows != m->cols) {
	fputs("gretl_matrix_add_self_transpose: matrix must be square\n", 
	      stderr);
	return GRETL_MATRIX_ERR;
    }

    for (i=0; i<m->rows; i++) {
	for (j=i; j<m->rows; j++) {
	    if (j == i) {
		m->val[mdx(m, i, j)] *= 2.0;
	    } else {
		x1 = m->val[mdx(m, i, j)];
		x2 = m->val[mdx(m, j, i)];
		m->val[mdx(m, i, j)] = 
		    m->val[mdx(m, j, i)] = x1 + x2;
	    }
	}
    }

    return GRETL_MATRIX_OK;
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

static int packed_idx (int nrows, int i, int j)
{
    int idx;

    if (i > j) {
	int tmp = i;

	i = j;
	j = tmp;
    }

    idx = nrows * i + j - i - ((i - 1) * i/2);
    return idx;
}

/* ....................................................... */

double gretl_matrix_get (const gretl_matrix *m, int i, int j)
{
    double x;

    if (m == NULL || m->val == NULL) return NADBL;

    if (i >= m->rows || j >= m->cols) return NADBL;

    if (m->packed) {
	x = m->val[packed_idx(m->rows, i, j)];
    } else {
	x = m->val[mdx(m, i, j)];
    }

    return x;
}

/* ....................................................... */

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x)
{
    if (m == NULL || m->val == NULL) return 1;

    if (i >= m->rows || j >= m->cols) return 1;

    if (m->packed) {
	m->val[packed_idx(m->rows, i, j)] = x;
    } else {
	m->val[mdx(m, i, j)] = x;
    }

    return 0;
}

void gretl_matrix_print (gretl_matrix *m, const char *msg, PRN *prn)
{
    int i, j;

    if (msg != NULL && *msg != '\0') {
	pprintf(prn, "%s\n\n", msg);
    }

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(m, i, j));
	}
	pputs(prn, "\n");
    }
    pputs(prn, "\n");
}

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
#if 0
    for (j=0; j<rows; j++) {
	for (i=0; i<cols; i++) {
	    m->val[p++] = X[i][j];
	}
    } 
#else
    for (j=0; j<cols; j++) {
	for (i=0; i<rows; i++) {
	    m->val[p++] = X[j][i];
	}
    }
#endif

    return m;
}

int gretl_matrix_multiply_mod (const gretl_matrix *a, int aflag,
			       const gretl_matrix *b, int bflag,
			       gretl_matrix *c)
{
    int i, j, k;
    int lrows, lcols;
    int rrows, rcols;
    int atr = (aflag == GRETL_MOD_TRANSPOSE);
    int btr = (bflag == GRETL_MOD_TRANSPOSE);
    int aidx, amax = a->rows * a->cols;
    int bidx, bmax = b->rows * b->cols;

    if (a == c || b == c) {
	fputs("gretl_matrix_multiply:\n product matrix must be "
	      "distinct from both input matrices\n", stderr);
	return GRETL_MATRIX_ERR;
    }

    lrows = (atr)? a->cols : a->rows;
    lcols = (atr)? a->rows : a->cols;
    rrows = (btr)? b->cols : b->rows;
    rcols = (btr)? b->rows : b->cols;

    if (lcols != rrows) {
	fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
	fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
		lrows, lcols, rrows, rcols, c->rows, c->cols);
	return GRETL_MATRIX_NON_CONFORM;
    }

    if (c->rows != lrows || c->cols != rcols) {
	fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
	fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
		lrows, lcols, rrows, rcols, c->rows, c->cols);
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<lrows; i++) {
	for (j=0; j<rcols; j++) {
	    c->val[mdx(c, i, j)] = 0.0;
	    for (k=0; k<lcols; k++) {
		aidx = (atr)? mdxtr(a,i,k) : mdx(a,i,k);
		bidx = (btr)? mdxtr(b,k,j) : mdx(b,k,j);
		if (aidx >= amax || bidx >= bmax) {
		    fputs("gretl_matrix_multiply_mod: index out of bounds\n", 
			  stderr);
		    return 1;
		}
		c->val[mdx(c,i,j)] += a->val[aidx] * b->val[bidx];
	    }
	}
    }

    return GRETL_MATRIX_OK;
}

static double 
gretl_matrix_column_mean (const gretl_matrix *m, int col)
{
    int i;
    double sum = 0.0;

    for (i=0; i<m->rows; i++) {
	sum += m->val[mdx(m, i, col)];
    }

    return sum / (double) m->rows;
}

/* Form a VCV matrix from matrix m (which is expected to have rows >= cols).
   Returns NULL on failure, allocated VCV on success.  Note that m is
   overwritten, the column means being subtracted.  It is up to the
   caller to free both m and VCV.
*/

gretl_matrix *gretl_matrix_vcv (gretl_matrix *m)
{
    int i, j, err = 0;
    double colmean;
    gretl_matrix *v;

    if (m->cols > m->rows) {
	fputs("gretl_matrix_vcv: expected rows >= cols\n", stderr);
	return NULL;
    }

    v = gretl_matrix_alloc(m->cols, m->cols);
    if (v == NULL) return NULL;

    /* subtract the column means from the column elements */
    for (j=0; j<m->cols; j++) {
	colmean = gretl_matrix_column_mean(m, j);
	for (i=0; i<m->rows; i++) {
	    m->val[mdx(m, i, j)] -= colmean;
	}
    }

    err = gretl_matrix_multiply_mod(m, GRETL_MOD_TRANSPOSE,
				    m, GRETL_MOD_NONE,
				    v);

    gretl_matrix_divide_by_scalar(v, (double) m->rows);

    if (err) {
	gretl_matrix_free(v);
	return NULL;
    }

    return v;
}

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c)
{
    return gretl_matrix_multiply_mod(a, GRETL_MOD_NONE,
				     b, GRETL_MOD_NONE,
				     c);
}

int gretl_matrix_cholesky_decomp (gretl_matrix *a)
{
    char uplo = 'L';
    integer n = a->rows;
    integer lda = a->rows;
    integer info;

    dpotrf_(&uplo, &n, a->val, &lda, &info);

#ifdef LAPACK_DEBUG
    printf("dpotrf: info = %d\n", (int) info);
#endif

    return (info != 0);
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

double *gretl_symmetric_matrix_eigenvals (gretl_matrix *m,
					  int eigenvecs) 
{
    integer n = m->rows;
    integer info;
    integer lwork;

    double *work;
    double *w;

    char jobz, uplo = 'U'; 

    if (eigenvecs) jobz = 'V';
    else jobz = 'N';

    work = malloc(sizeof *work);
    if (work == NULL) {
	return NULL;
    }

    w = malloc(n * sizeof *w);
    if (w == NULL) {
	free(work);
	return NULL;
    }

    lwork = -1; /* find optimal workspace size */
    dsyev_(&jobz, &uplo, &n, m->val, &n, 
	   w, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	free(work);
	free(w);
	return NULL;
    }	

    lwork = (integer) work[0];

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(w);
	return NULL;
    }    

    dsyev_(&jobz, &uplo, &n, m->val, &n, 
	   w, work, &lwork, &info);

    if (info != 0) {
	free(w);
	w = NULL;
    }

    free(work);

    return w;
}

void gretl_matrix_set_int (gretl_matrix *m, int t)
{
    m->t = t;
}

int gretl_matrix_get_int (const gretl_matrix *m)
{
    return m->t;
}

int gretl_vector_get_length (const gretl_vector *v) 
{
    return v->cols;
}

int gretl_matrix_cols (const gretl_matrix *m)
{
    return m->cols;
}

int gretl_matrix_rows (const gretl_matrix *m)
{
    return m->rows;
}
