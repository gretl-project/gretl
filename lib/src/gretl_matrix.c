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

/* LAPACK-based matrix routines for gretl */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "gretl_matrix_private.h"

static const char *wspace_fail = "gretl_matrix: workspace query failed\n";

static int packed_idx (int nrows, int i, int j);

#define gretl_is_vector(v) (v->rows == 1 || v->cols == 1)

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

/**
 * gretl_matrix_alloc:
 * @rows: desired number of rows in matrix.
 * @cols: desired number of columns.
 *
 * Returns: pointer to a newly allocated gretl_matrix, or %NULL
 * on failure.  Note that the actual data storage is not
 * initialized.
 * 
 */

gretl_matrix *gretl_matrix_alloc (int rows, int cols)
{
    return real_gretl_matrix_alloc(rows, cols, 0);
}

/**
 * gretl_packed_matrix_alloc:
 * @rows: desired number of rows and columns in (symmetric) matrix.
 *
 * Returns: pointer to a newly allocated gretl_matrix, or %NULL
 * on failure.  The matrix is in packed (triangular) form.
 * 
 */

gretl_matrix *gretl_packed_matrix_alloc (int rows)
{
    return real_gretl_matrix_alloc(rows, rows, 1);
}

/**
 * gretl_column_vector_from_array:
 * @x: pointer to array of elements.
 * @n: number of elements.
 * @mod: modifier flag: either %GRETL_MOD_NONE, or %GRETL_MOD_SQUARE
 * to use the squares of the elements of @x.
 *
 * Returns: pointer to a newly allocated gretl_vector containing
 * the elements of x (or their squares), or %NULL on failure.  
 * Missing elements of x are skipped.
 * 
 */

gretl_vector *
gretl_column_vector_from_array (const double *x, int n, int mod)
{
    gretl_matrix *v;
    double xi;
    int i = 0;
    
    v = gretl_column_vector_alloc(n);
    if (v == NULL) return NULL;

    while (i < n) {
	xi = *x++;
	if (!na(xi)) {
	    if (mod == GRETL_MOD_SQUARE) {
		v->val[i] = xi * xi;
	    } else {
		v->val[i] = xi;
	    }
	    i++;
	}
    }

    return v;
}

/* ....................................................... */

static gretl_matrix *
gretl_matrix_copy_mod (const gretl_matrix *m, int mod)
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

/**
 * gretl_matrix_copy:
 * @m: source matrix to be copied.
 *
 * Returns: an allocated copy of matrix @m, or %NULL on failure.  
 * 
 */

gretl_matrix *gretl_matrix_copy (const gretl_matrix *m)
{
    return gretl_matrix_copy_mod(m, GRETL_MOD_NONE);
}

/**
 * gretl_matrix_free:
 * @m: matrix to be freed.
 *
 * Frees the allocated storage in @m, then @m itself.
 * 
 */

void gretl_matrix_free (gretl_matrix *m)
{
    if (m == NULL) return;

    if (m->val != NULL) free(m->val);
    free(m);
}

/**
 * gretl_matrix_zero:
 * @m: matrix to be set to zero.
 *
 * Sets all elements of @m to zero.
 * 
 */

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

static int gretl_matrix_zero_triangle (gretl_matrix *m, char t)
{
    int i, j;

    if (m == NULL || m->val == NULL || m->packed) 
	return GRETL_MATRIX_ERR;

    if (m->rows != m->cols)
	return GRETL_MATRIX_NON_CONFORM;

    if (t == 'U') {
	for (i=0; i<m->rows-1; i++) {
	    for (j=i+1; j<m->cols; j++) {
		m->val[mdx(m, i, j)] = 0.0;
	    }
	}
    } else {
	for (i=1; i<m->rows; i++) {
	    for (j=0; j<i; j++) {
		m->val[mdx(m, i, j)] = 0.0;
	    }
	}
    }

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_zero_upper:
 * @m: square matrix to operate on.
 *
 * Sets the elements of @m outside of the lower triangle to zero.
 * 
 * Returns: %GRETL_MATRIX_OK on success, non-zero error code otherwise.
 * 
 */

int gretl_matrix_zero_upper (gretl_matrix *m)
{
    return gretl_matrix_zero_triangle(m, 'U');
}

/**
 * gretl_matrix_zero_lower:
 * @m: square matrix to operate on.
 *
 * Sets the elements of @m outside of the upper triangle to zero.
 * 
 * Returns: %GRETL_MATRIX_OK on success, non-zero error code otherwise.
 * 
 */

int gretl_matrix_zero_lower (gretl_matrix *m)
{
    return gretl_matrix_zero_triangle(m, 'L');
}

/**
 * gretl_matrix_multiply_by_scalar:
 * @m: matrix to operate on.
 * @x: scalar by which to multiply.
 *
 * Multiplies all elements of @m by @x.
 * 
 */

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

/**
 * gretl_matrix_divide_by_scalar:
 * @m: matrix to operate on.
 * @x: scalar by which to divide.
 *
 * Divides all elements of @m by @x.
 * 
 */

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

/**
 * gretl_matrix_copy_values:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Copies the elements of @src into the corresponding elements
 * of @targ.
 * 
 * Returns: %GRETL_MATRIX_OK on successful completion, or
 * %GRETL_MATRIX_NON_CONFORM if the two matrices are not
 * conformable for the operation.
 * 
 */

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

/**
 * gretl_matrix_add_to:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Adds the elements of @src to the corresponding elements
 * of @targ.
 * 
 * Returns: %GRETL_MATRIX_OK on successful completion, or
 * %GRETL_MATRIX_NON_CONFORM if the two matrices are not
 * conformable for the operation.
 * 
 */

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

/**
 * gretl_matrix_subtract_from:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Subtracts the elements of @src from the corresponding elements
 * of @targ.
 * 
 * Returns: %GRETL_MATRIX_OK on successful completion, or
 * %GRETL_MATRIX_NON_CONFORM if the two matrices are not
 * conformable for the operation.
 * 
 */

int 
gretl_matrix_subtract_from (gretl_matrix *targ, const gretl_matrix *src)
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
    
    for (i=0; i<n; i++) targ->val[i] -= src->val[i];

    return GRETL_MATRIX_OK;
}

/**
 * gretl_square_matrix_transpose:
 * @m: square matrix to operate on.
 *
 * Transposes the matrix @m.
 *
 * Returns: %GRETL_MATRIX_OK on success, non-zero error code otherwise.
 */

int gretl_square_matrix_transpose (gretl_matrix *m)
{
    int i, j;
    double x;

    if (m->rows != m->cols) {
	fputs("gretl_square_matrix_transpose: matrix must be square\n", 
	      stderr);
	return GRETL_MATRIX_ERR;
    }

    for (i=0; i<m->rows-1; i++) {
	for (j=i+1; j<m->rows; j++) {
	    x = m->val[mdx(m, i, j)];
	    m->val[mdx(m, i, j)] = m->val[mdx(m, j, i)];
	    m->val[mdx(m, j, i)] = x;
	}
    }

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_add_self_transpose:
 * @m: (square) matrix to operate on.
 *
 * Adds the transpose of @m to @m itself, yielding a symmetric
 * matrix.
 * 
 * Returns: %GRETL_MATRIX_OK on successful completion, or
 * %GRETL_MATRIX_ERR if the source matrix is not square.
 * 
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

/**
 * gretl_matrix_steal_data:
 * @m: matrix to operate on.
 *
 * "Steals" the allocated data from @m, which is left with a
 * %NULL data pointer.
 * 
 * Returns: a pointer to the "stolen" data.
 * 
 */

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

/**
 * gretl_matrix_get:
 * @m: matrix to get data from.
 * @i: row index (zero-based).
 * @j: column index (zero-based).
 *
 * Retrieves the data value from the specified row and column
 * of @m.
 * 
 * Returns: the data value, or #NADBL if the indices are out
 * of range for @m.
 * 
 */

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

/**
 * gretl_vector_get:
 * @v: vector to get data from.
 * @i: index (zero-based).
 *
 * Retrieves the ith value from the specified vector.
 * 
 * Returns: the data value, or #NADBL if the index is out
 * of range for @v.
 * 
 */

double gretl_vector_get (const gretl_vector *v, int i)
{
    if (v == NULL || v->val == NULL) return NADBL;

    if (i >= v->rows && i >= v->cols) return NADBL;

    return v->val[i];
}

/**
 * gretl_matrix_set:
 * @m: matrix to operate on.
 * @i: row index (zero-based).
 * @j: column index (zero-based).
 * @x: value to set.
 *
 * Sets element i, j of @m to the value @x.
 * 
 * Returns: %GRETL_MATRIX_OK, or %GRETL_MATRIX_ERR if the
 * indices are out of range for @m.
 * 
 */

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x)
{
    if (m == NULL || m->val == NULL) return GRETL_MATRIX_ERR;

    if (i >= m->rows || j >= m->cols) return GRETL_MATRIX_ERR;

    if (m->packed) {
	m->val[packed_idx(m->rows, i, j)] = x;
    } else {
	m->val[mdx(m, i, j)] = x;
    }

    return GRETL_MATRIX_OK;
}

/**
 * gretl_vector_set:
 * @v: vector to operate on.
 * @i: index (zero-based).
 * @x: value to set.
 *
 * Sets element i of @v to the value @x.
 * 
 * Returns: %GRETL_MATRIX_OK, or %GRETL_MATRIX_ERR if the
 * given index is out of range for @v.
 * 
 */

int gretl_vector_set (gretl_vector *v, int i, double x)
{
    if (v == NULL || v->val == NULL) return GRETL_MATRIX_ERR;

    if (i >= v->rows && i >= v->cols) return GRETL_MATRIX_ERR;

    v->val[i] = x;

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_print:
 * @m: matrix to print.
 * @msg: accompanying message text (or %NULL if no message is wanted).
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn.
 * 
 */

void gretl_matrix_print (const gretl_matrix *m, const char *msg, PRN *prn)
{
    int i, j;
    PRN myprn;

    if (prn == NULL) {
	myprn.fp = stdout;
	myprn.fpaux = NULL;
	myprn.buf = NULL;
	prn = &myprn;
    }

    if (msg != NULL && *msg != '\0') {
	pprintf(prn, "%s\n\n", msg);
    }

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%#12.5g ", gretl_matrix_get(m, i, j));
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/* calculate determinant using LU factorization.
   if logdet != 0, return the log of the determinant.
   if logdet != 0 and absval != 0, return the log of the
   absolute value of the determinant.  
*/   

static double gretl_LU_determinant (gretl_matrix *a, int logdet, int absval)
{
    integer info;
    integer m = a->rows;
    integer n = a->cols;
    integer *ipiv;
    double det;
    int i;

    if (a->rows != a->cols) {
	fputs("gretl_LU_determinant: matrix must be square\n", stderr);
	return NADBL;
    }

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) return NADBL;

    dgetrf_(&m, &n, a->val, &n, ipiv, &info);

    if (info != 0) {
	fprintf(stderr, "gretl_LU_determinant: dgetrf gave info = %d\n", 
		(int) info);
	free(ipiv);
	return NADBL;
    }

    if (logdet) {
	int negcount = 0;

	/* not sure if this is worth the bother: but do we get a bit
	   more numerical accuracy in some cases by adding logs,
	   rather than by multiplying terms then taking the log of the
	   product? 
	*/

	det = 0.0;
	for (i=0; i<n; i++) {
	    double aii = a->val[mdx(a, i, i)];

	    if (aii == 0.0) {
		fputs("gretl_matrix_log_determinant: determinant = 0\n", stderr);
		det = NADBL;
		break;
	    }

	    if (ipiv[i] != i + 1) {
		aii = -aii;
	    }
	    if (aii < 0) {
		aii = -aii;
		negcount++;
	    }
	    det += log(aii);
	}
	if (!absval && negcount % 2) {
	    fputs("gretl_matrix_log_determinant: determinant is < 0\n", stderr);
	    det = NADBL;
	}
    } else {
	det = 1.0;
	for (i=0; i<n; i++) {
	    if (ipiv[i] != i + 1) {
		det = -det;
	    }
	    det *= a->val[mdx(a, i, i)];
	}
    }

    free(ipiv);

    return det;
}

/**
 * gretl_matrix_determinant:
 * @a: gretl_matrix.
 *
 * Compute the determinant of the square matrix @a using the LU
 * factorization.  Matrix @a is not preserved: it is overwritten
 * by the factorization.
 * 
 * Returns: the determinant, or #NABDL on failure.
 * 
 */

double gretl_matrix_determinant (gretl_matrix *a)
{
    return gretl_LU_determinant(a, 0, 0);
}

/**
 * gretl_matrix_log_determinant:
 * @a: gretl_matrix.
 *
 * Compute the log of the determinant of the square matrix @a using LU
 * factorization.  Matrix @a is not preserved: it is overwritten
 * by the factorization.
 * 
 * Returns: the determinant, or #NABDL on failure.
 * 
 */

double gretl_matrix_log_determinant (gretl_matrix *a)
{
    return gretl_LU_determinant(a, 1, 0);
}

/**
 * gretl_matrix_log_abs_determinant:
 * @a: gretl_matrix.
 *
 * Compute the log of the absolute value of the determinant of the 
 * square matrix @a using LU factorization.  Matrix @a is not 
 * preserved: it is overwritten by the factorization.
 * 
 * Returns: the determinant, or #NABDL on failure.
 * 
 */

double gretl_matrix_log_abs_determinant (gretl_matrix *a)
{
    return gretl_LU_determinant(a, 1, 1);
}

/**
 * gretl_LU_solve:
 * @a: gretl_matrix.
 * @b: gretl_vector.
 *
 * Solves ax = b for the unknown vector x, using LU decomposition.
 * On exit, b is replaced by the solution and a is replaced by its 
 * LU decomposition.
 * 
 * Returns: 0 on successful completion, or a non-zero error code
 * (from the lapack function %dgetrs) on error.
 * 
 */

int gretl_LU_solve (gretl_matrix *a, gretl_vector *b)
{
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
	fprintf(stderr, "gretl_LU_solve: dgetrf gave info = %d\n", 
		(int) info);
	free(ipiv);
	return info;
    }

    dgetrs_(&trans, &n, &nrhs, a->val, &n, ipiv, b->val, &ldb, &info);

    free(ipiv);

    return info;
}

/**
 * gretl_matrix_from_2d_array:
 * @X: two-dimensional array of doubles.
 * @rows: number of rows in target matrix.
 * @cols: number of columns in target matrix.
 *
 * Returns: allocated gretl_matrix, the elements of which are set to
 * the values in @X, or %NULL on allocation failure.
 * 
 */

gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
					  int rows, int cols)
{
    int i, j, p;
    gretl_matrix *m;

    m = gretl_matrix_alloc(rows, cols);
    if (m == NULL) return m;

    p = 0;
    for (j=0; j<cols; j++) {
	for (i=0; i<rows; i++) {
	    m->val[p++] = X[j][i];
	}
    }

    return m;
}

/**
 * gretl_matrix_multiply_mod:
 * @a: left-hand matrix.
 * @aflag: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @b: right-hand matrix.
 * @bflag: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @c: matrix to hold the product.
 * 
 * Multiplies @a (or a-transpose) into @b (or b transpose),
 * with the result written into @c.
 *
 * Returns: %GRETL_MATRIX_OK on success; non-zero error code on
 * failure.
 * 
 */

int gretl_matrix_multiply_mod (const gretl_matrix *a, int aflag,
			       const gretl_matrix *b, int bflag,
			       gretl_matrix *c)
{
    register int i, j, k;
    int lrows, lcols;
    int rrows, rcols;
    const int atr = (aflag == GRETL_MOD_TRANSPOSE);
    const int btr = (bflag == GRETL_MOD_TRANSPOSE);
    int aidx, bidx;
    double targ;

    if (a == c || b == c) {
	fputs("gretl_matrix_multiply:\n product matrix must be "
	      "distinct from both input matrices\n", stderr);
	fprintf(stderr, "a = %p, b = %p, c = %p\n", 
		(void *) a, (void *) b, (void *) c);
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
	    targ = 0.0;
	    for (k=0; k<lcols; k++) {
		aidx = (atr)? mdxtr(a,i,k) : mdx(a,i,k);
		bidx = (btr)? mdxtr(b,k,j) : mdx(b,k,j);
		targ += a->val[aidx] * b->val[bidx];
	    }
	    c->val[mdx(c,i,j)] = targ;
	}
    }

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_dot_product:
 * @a: left-hand matrix.
 * @aflag: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @b: right-hand matrix.
 * @bflag: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @err: pointer to integer error-code variable.
 * 
 * Returns: The dot (scalar) product of @a (or a-transpose) and
 * @b (or b-transpose).  @err points to %GRETL_MATRIX_ERR on
 * failure.
 * 
 */

double gretl_matrix_dot_product (const gretl_matrix *a, int aflag,
				 const gretl_matrix *b, int bflag,
				 int *err)
{
    gretl_matrix *c;
    double ret = NADBL;

    c = gretl_matrix_alloc(1, 1);
    if (c == NULL) {
	*err = GRETL_MATRIX_ERR;
	return ret;
    }

    *err = gretl_matrix_multiply_mod(a, aflag, b, bflag, c);
    if (*err == GRETL_MATRIX_OK) {
	ret = c->val[0];
    }
	
    gretl_matrix_free(c);

    return ret;
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

/**
 * gretl_matrix_vcv:
 * @m: source matrix (expected to have rows >= cols).
 *
 * Forms a variance-covariance matrix based on @m, thus:
 * - subtract the column means from the column elements of @m.
 * - multiply @m-transpose into @m.
 * - divide the elements of the product by the number of rows
 *   in @m.
 * 
 * Returns: the allocated variance-covariance matrix, or %NULL
 * on failure.
 * 
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

    /* v = m'm */
    err = gretl_matrix_multiply_mod(m, GRETL_MOD_TRANSPOSE,
				    m, GRETL_MOD_NONE,
				    v);

    if (err) {
	gretl_matrix_free(v);
	return NULL;
    } else {
	gretl_matrix_divide_by_scalar(v, (double) m->rows);
    }

    return v;
}

/**
 * gretl_matrix_multiply:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @c: matrix to hold the product.
 * 
 * Multiplies @a into @b, with the result written into @c.
 *
 * Returns: %GRETL_MATRIX_OK on success; non-zero error code on
 * failure.
 * 
 */

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c)
{
    return gretl_matrix_multiply_mod(a, GRETL_MOD_NONE,
				     b, GRETL_MOD_NONE,
				     c);
}

/**
 * gretl_matrix_cholesky_decomp:
 * @a: matrix to operate on.
 * 
 * Computes the Cholesky factorization of the symmetric,
 * positive definite matrix @a.  On exit the lower triangle of 
 * @a is replaced by the factor L, as in a = LL'.
 * Uses the lapack function %dpotrf.
 *
 * Returns: %GRETL_MATRIX_OK on success; %GRETL_MATRIX_ERR on 
 * failure.
 * 
 */

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

    return (info == 0)? GRETL_MATRIX_OK : GRETL_MATRIX_ERR;
}

/**
 * gretl_invert_general_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a general matrix using LU
 * factorization.  On exit @a is overwritten with the inverse.
 * Uses the lapack functions %dgetrf and %dgetri.
 *
 * Returns: 0 on success; non-zero error code on failure.
 * 
 */

int gretl_invert_general_matrix (gretl_matrix *a)
{
    integer m = a->rows;
    integer n = a->cols;
    integer info;
    integer lwork;
    integer *ipiv;
    int lipiv;
    int err = 0;

    double *work;

    if (m <= n) lipiv = m;
    else lipiv = n;

    ipiv = malloc(lipiv * sizeof *ipiv);
    if (ipiv == NULL) {
	return GRETL_MATRIX_NOMEM;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return GRETL_MATRIX_NOMEM;
    }    

    dgetrf_(&m, &n, a->val, &m, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	return GRETL_MATRIX_SINGULAR;
    }

    lwork = -1;
    dgetri_(&n, a->val, &n, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	free(ipiv);
	return GRETL_MATRIX_ERR;
    }

    lwork = (integer) work[0];

#ifdef LAPACK_DEBUG
    printf("dgetri: workspace = %d\n", (int) lwork);
#endif

    work = realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return GRETL_MATRIX_NOMEM;
    }  

    dgetri_(&n, a->val, &n, ipiv, work, &lwork, &info);

#ifdef LAPACK_DEBUG
    printf("dgetri: info = %d\n", (int) info);
#endif

    free(work);
    free(ipiv);

    if (info) {
	err = GRETL_MATRIX_SINGULAR;
    }

    return err;
}

/* In the case of symmetric matrices, the lapack functions tend
   to process only either the upper or lower triangle.  This
   function "expands" the solution, reconstituting the matrix
   as symmetric. 
*/

static int 
gretl_symmetric_matrix_expand (gretl_matrix *m, char uplo)
{
    int i, j, n;
    double x;

    if (m->cols != m->rows) {
	fputs("gretl_symmetric_matrix_expand: input is not square\n",
	      stderr);
	return 1;
    }

    n = m->rows;

    for (i=0; i<n; i++) {
	for (j=i+1; j<n; j++) {
	    if (uplo == 'U') {
		x = m->val[mdx(m, i, j)];
		m->val[mdx(m, j, i)] = x;
	    } else {
		x = m->val[mdx(m, j, i)];
		m->val[mdx(m, i, j)] = x;
	    }
	}
    }

    return 0;
}

/**
 * gretl_invert_symmetric_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization.  On exit @a is overwritten with 
 * the inverse. Uses the lapack functions %dpotrf and %dpotri.
 *
 * Returns: 0 on success; non-zero error code on failure.
 * 
 */

int gretl_invert_symmetric_matrix (gretl_matrix *a)
{
    integer n, info;
    char uplo = 'U';
    int err = 0;

    if (a->cols != a->rows) {
	fputs("gretl_invert_symmetric_matrix: input is not square\n",
	      stderr);
	return GRETL_MATRIX_NON_CONFORM;
    }

    n = a->cols;

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fputs("gretl_invert_symmetric_matrix: dpotrf failed\n", stderr);
	return GRETL_MATRIX_SINGULAR;
    }

    dpotri_(&uplo, &n, a->val, &n, &info);

#ifdef LAPACK_DEBUG
    printf("dpotri: info = %d\n", (int) info);
#endif
    
    if (info != 0) {
	err = GRETL_MATRIX_SINGULAR;
	fputs("gretl_invert_symmetric_matrix: dpotrf failed\n", stderr);
    } else {
	gretl_symmetric_matrix_expand(a, uplo);
    }

    return err;
}

/**
 * gretl_general_matrix_eigenvals:
 * @m: matrix to operate on.
 * @ev: matrix to store eigenvectors, or %NULL if the eigenvectors
 * are not required.
 * 
 * Computes the eigenvalues of the general matrix @m.  If @ev is
 * non-%NULL, write the right eigenvectors of @m into @ev.
 * Uses the lapack function %dgeev.
 *
 * Returns: allocated storage containing the eigenvalues, or %NULL
 * on failure.
 * 
 */

double *gretl_general_matrix_eigenvals (gretl_matrix *m, gretl_matrix *ev) 
{
    integer n = m->rows;
    integer info;
    integer lwork;
    integer nvr, nvl = 2;
    char jvr, jvl = 'N';
    double *work, *work2;
    double *wr = NULL, *wi = NULL, *vr = NULL;
    double nullvl[2] = {0.0};
    double nullvr[2] = {0.0};

    work = malloc(sizeof *work);
    if (work == NULL) {
	return NULL;
    }

    wr = malloc(n * sizeof *wr);
    wi = malloc(n * sizeof *wi);
    if (wr == NULL || wi == NULL) {
	goto bailout;
    }

    if (ev != NULL) {
	/* eigenvectors wanted */
	vr = ev->val;
	nvr = n;
	jvr = 'V';
    } else {
	vr = nullvr;
	nvr = 2;
	jvr = 'N';
    }	

    lwork = -1; /* find optimal workspace size */
    dgeev_(&jvl, &jvr, &n, m->val, &n, wr, wi, nullvl, 
	   &nvl, vr, &nvr, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	goto bailout;
    } else {
	work = work2;
    }

    dgeev_(&jvl, &jvr, &n, m->val, &n, wr, wi, nullvl, 
	   &nvl, vr, &nvr, work, &lwork, &info);

    if (info != 0) {
	goto bailout;
    } 

    free(wi);
    free(work);

    return wr;

 bailout:
    free(work);
    free(wr);
    free(wi);

    return NULL;    
}

/**
 * gretl_symmetric_matrix_eigenvals:
 * @m: matrix to operate on.
 * @eigenvecs: non-zero to calculate eigenvectors, 0 to omit.
 * 
 * Computes the eigenvalues of the real symmetric matrix @m.  
 * If @eigenvecs is non-zero, also compute the orthonormal
 * eigenvectors of @m, which are stored in @m. Uses the lapack 
 * function %dsyev.
 *
 * Returns: allocated storage containing the eigenvalues, or %NULL
 * on failure.
 * 
 */

double *
gretl_symmetric_matrix_eigenvals (gretl_matrix *m, int eigenvecs) 
{
    integer n = m->rows;
    integer info;
    integer lwork;

    double *work, *work2;
    double *w;

    char uplo = 'U', jobz = (eigenvecs)? 'V' : 'N';

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

    work2 = realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	free(work);
	free(w);
	return NULL;
    } else {
	work = work2;
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

/**
 * gretl_matrix_set_int:
 * @m: matrix to operate on.
 * @t: value to set
 * 
 * Sets an integer value on the gretl_matrix (used for internal
 * information).  
 * 
 */

void gretl_matrix_set_int (gretl_matrix *m, int t)
{
    m->t = t;
}

/**
 * gretl_matrix_get_int:
 * @m: matrix to read from.
 * 
 * Returns the integer that has been set on the matrix, or zero.
 * 
 */

int gretl_matrix_get_int (const gretl_matrix *m)
{
    return m->t;
}

int gretl_vector_get_length (const gretl_vector *v) 
{
    return (v->cols > v->rows)? v->cols : v->rows;
}

/**
 * gretl_matrix_cols:
 * @m: matrix to query.
 * 
 * Returns: the number of columns in @m. 
 * 
 */

int gretl_matrix_cols (const gretl_matrix *m)
{
    return m->cols;
}

/**
 * gretl_matrix_rows:
 * @m: matrix to query.
 * 
 * Returns: the number of rows in @m. 
 * 
 */

int gretl_matrix_rows (const gretl_matrix *m)
{
    return m->rows;
}

/* ....................................................... */

static int
get_svd_ols_vcv (const gretl_matrix *A, const gretl_matrix *B,
		 const double *s, gretl_matrix *vcv, double *s2)
{
    double sigma2;
    double aik, ajk, vij;
    int m = A->cols;
    int T = A->rows;
    int i, j, k;

    sigma2 = 0.0;
    for (i=m; i<T; i++) {
	sigma2 += B->val[i] * B->val[i];
    }
    sigma2 /= T - m;

    /* Get X'X, based on the work done by the SV decomp:
       squares of singular values, premultiplied by V and
       postmultiplied by V-transpose 
    */

    for (i=0; i<m; i++) {
	for (j=i; j<m; j++) {
	    vij = 0.0;
	    for (k=0; k<m; k++) {
		aik = A->val[mdx(A, k, i)];
		ajk = A->val[mdx(A, k, j)];
		vij += aik * s[k] * s[k] * ajk;
	    }
	    vcv->val[mdx(vcv, i, j)] = vij;
	    if (j != i) {
		vcv->val[mdx(vcv, j, i)] = vij;
	    }
	}
    }

    if (gretl_invert_symmetric_matrix(vcv)) return 1;

    gretl_matrix_multiply_by_scalar(vcv, sigma2);  

    if (s2 != NULL) {
	*s2 = sigma2;
    }

    return 0;
}

static int
get_ols_vcv (const gretl_vector *y, const gretl_matrix *X,
	     const gretl_vector *b, gretl_matrix *vcv,
	     double *s2)
{
    double u, sigma2 = 0.0;
    int k = X->cols;
    int n = X->rows;
    int i, j;

    if (gretl_invert_symmetric_matrix(vcv)) return 1;

    for (i=0; i<n; i++) {
	u = y->val[i];
	for (j=0; j<k; j++) {
	    u -= X->val[mdx(X, i, j)] * b->val[j];
	}
	sigma2 += u * u;
    }

    sigma2 /= (n - k);

    gretl_matrix_multiply_by_scalar(vcv, sigma2);  

    if (s2 != NULL) {
	*s2 = sigma2;
    }

    return 0;
}

static void
get_ols_uhat (const gretl_vector *y, const gretl_matrix *X,
	      const gretl_vector *b, gretl_vector *uhat)
{
    int ncoeff = gretl_vector_get_length(b);
    int n = gretl_vector_get_length(uhat);
    int i, j;
    double uh;

    for (i=0; i<n; i++) {
	uh = y->val[i];
	for (j=0; j<ncoeff; j++) {
	    uh -= b->val[j] * X->val[mdx(X, i, j)];
	}
	uhat->val[i] = uh;
    }
}

int gretl_matrix_svd_ols (const gretl_vector *y, const gretl_matrix *X,
			  gretl_vector *b, gretl_matrix *vcv,
			  gretl_vector *uhat, double *s2)
{
    gretl_vector *A = NULL;
    gretl_matrix *B = NULL;

    int T = X->rows;
    int k = X->cols;

    integer m = T;
    integer n = k;
    integer nrhs = 1;
    integer lda = T;
    integer ldb = T;
    integer lwork = -1;
    integer rank;
    integer info;

    double rcond = -1.0;
    double *work = NULL;
    double *work2 = NULL;
    double *s = NULL;

    int err = GRETL_MATRIX_OK;

    if (gretl_vector_get_length(b) != k) {
	err = GRETL_MATRIX_NON_CONFORM;
	goto bailout;
    }

    A = gretl_matrix_copy(X);
    if (A == NULL) {
	err = GRETL_MATRIX_NOMEM;
	goto bailout;
    }

    B = gretl_matrix_copy(y);
    if (B == NULL) {
	err = GRETL_MATRIX_NOMEM;
	goto bailout;
    }

    /* for singular values of A */
    s = malloc(k * sizeof *s);
    if (s == NULL) {
	err = GRETL_MATRIX_NOMEM; 
	goto bailout;
    }

    work = malloc(sizeof *work);
    if (work == NULL) {
	err = GRETL_MATRIX_NOMEM; 
	goto bailout;
    } 

    /* workspace query */
    dgelss_(&m, &n, &nrhs, A->val, &lda, B->val, &ldb, s, &rcond,
	    &rank, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	fputs(wspace_fail, stderr);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	err = GRETL_MATRIX_NOMEM; 
	goto bailout;
    } 

    work = work2;

    /* get actual solution */
    dgelss_(&m, &n, &nrhs, A->val, &lda, B->val, &ldb, s, &rcond,
	    &rank, work, &lwork, &info);

    if (info != 0) {
	err = GRETL_MATRIX_ERR;
    }

    if (rank < k) {
	fprintf(stderr, "dgelss: rank of X = %d\n", (int) rank);
    }

    if (!err) {
	int i;

	for (i=0; i<k; i++) {
	    b->val[i] = B->val[i];
	}
	if (vcv != NULL) {
	    err = get_svd_ols_vcv(A, B, s, vcv, s2);
	}
	if (uhat != NULL) {
	    get_ols_uhat(y, X, b, uhat);
	}
    }

    bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(B);
    free(work);
    free(s);

    return err;
}

/**
 * gretl_matrix_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or %NULL if this is not needed.
 * @uhat: vector to hold the regression residuals, or %NULL if 
 * these are not needed.
 * @s2: pointer to receive residual variance, or %NULL.
 *
 * Computes OLS estimates using LU factorization, and puts the
 * coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv and the residuals in @uhat.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 * 
 */

int gretl_matrix_ols (const gretl_vector *y, const gretl_matrix *X,
		      gretl_vector *b, gretl_matrix *vcv,
		      gretl_vector *uhat, double *s2)
{
    gretl_vector *XTy = NULL;
    gretl_matrix *XTX = NULL;
    int k = X->cols;
    int err = GRETL_MATRIX_OK;

    if (gretl_vector_get_length(b) != k) {
	err = GRETL_MATRIX_NON_CONFORM;
    }

    if (!err) {
	XTy = gretl_column_vector_alloc(k);
	if (XTy == NULL) err = GRETL_MATRIX_NOMEM;
    }

    if (!err) {
	XTX = gretl_matrix_alloc(k, k);
	if (XTX == NULL) err = GRETL_MATRIX_NOMEM;
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					y, GRETL_MOD_NONE,
					XTy);
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					X, GRETL_MOD_NONE,
					XTX);
    }

    if (!err && vcv != NULL) {
	err = gretl_matrix_copy_values(vcv, XTX);
    }

    if (!err) {
	err = gretl_LU_solve(XTX, XTy);
    }

    if (!err) {
	int i;
	
	for (i=0; i<k; i++) {
	    b->val[i] = XTy->val[i];
	}
	if (vcv != NULL) {
	    err = get_ols_vcv(y, X, b, vcv, s2);
	}
	if (uhat != NULL) {
	    get_ols_uhat(y, X, b, uhat);
	}
    }

    if (XTy != NULL) gretl_vector_free(XTy);
    if (XTX != NULL) gretl_matrix_free(XTX);

    return err;
}

/**
 * gretl_scalar_b_prime_X_b:
 * @b: column k-vector.
 * @X: k x k matrix.
 * @err: pointer to error code variable.
 *
 * Computes the scalar product, @b transpose times @X times @b.
 * On success, *err = 0, otherwise it is non-zero.
 * 
 * Returns: scalar product.
 * 
 */

double gretl_scalar_b_prime_X_b (const gretl_vector *b, const gretl_matrix *X,
				 int *err)
{
    gretl_matrix *tmp;
    double ret = NADBL;

    *err = 0;

    if (b->rows != X->rows ||
	X->rows != X->cols ||
	b->cols != 1) {
	*err = GRETL_MATRIX_NON_CONFORM;
	return ret;
    }

    tmp = gretl_matrix_alloc(1, b->rows);
    if (tmp == NULL) {
	*err = GRETL_MATRIX_NOMEM;
	return ret;
    }

    *err = gretl_matrix_multiply_mod(b, GRETL_MOD_TRANSPOSE,
				     X, GRETL_MOD_NONE,
				     tmp);

    if (!*err) {
	ret = gretl_matrix_dot_product(tmp, GRETL_MOD_NONE,
				       b, GRETL_MOD_NONE,
				       err);
    }

    if (*err) ret = NADBL;

    gretl_matrix_free(tmp);

    return ret;
}

/**
 * gretl_matrix_A_X_A_prime:
 * @A: m * k matrix.
 * @X: k * k matrix.
 * @err: pointer to error code variable.
 *
 * Computes A * X * A'.
 * 
 * Returns: m * m matrix product, or %NULL on error (in which case
 * @err will point to a non-zero error code).
 * 
 */

gretl_matrix *
gretl_matrix_A_X_A_prime (const gretl_matrix *A, const gretl_matrix *X,
			  int *err)
{
    gretl_matrix *tmp = NULL;
    gretl_matrix *ret = NULL;
    int m = A->rows;
    int k = A->cols;

    *err = 0;

    if (X->rows != k || X->cols != k) {
	*err = GRETL_MATRIX_NON_CONFORM;
	return NULL;
    }

    tmp = gretl_matrix_alloc(m, k);
    ret = gretl_matrix_alloc(m, m);

    if (tmp == NULL || ret == NULL) {
	gretl_matrix_free(tmp);
	gretl_matrix_free(ret);
	*err = GRETL_MATRIX_NOMEM;
	return NULL;
    }

    *err = gretl_matrix_multiply_mod(A, GRETL_MOD_NONE,
				     X, GRETL_MOD_NONE,
				     tmp);

    if (!*err) {
	*err = gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
					 A, GRETL_MOD_TRANSPOSE,
					 ret);
    }

    gretl_matrix_free(tmp);

    if (*err) {
	gretl_matrix_free(ret);
	ret = NULL;
    }

    return ret;
}

static int count_selection (const char *s, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
	if (s[i] != 0) c++;
    }

    return c;
}

/**
 * gretl_vcv_matrix_from_model:
 * @pmod: pointer to model
 * @select: char array indicating which rows and colums to select
 * (or %NULL for the full matrix).
 *
 * Produces all or part of the covariance matrix for @pmod 
 * in the form of a gretl_matrix.  Storage is allocated, to be freed
 * by the caller.  If @select is not %NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows (and columns), and zero elements otherwise.
 * 
 * Returns: the covariance matrix, or %NULL on error.
 * 
 */

gretl_matrix *
gretl_vcv_matrix_from_model (MODEL *pmod, const char *select)
{
    gretl_matrix *vcv;
    int i, j, idx, nc;
    int ii, jj;
    int k = pmod->ncoeff;

    /* first ensure the model _has_ a vcv */
    if (makevcv(pmod)) return NULL;

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }
    
    if (nc == 0) return NULL;

    vcv = gretl_matrix_alloc(nc, nc);
    if (vcv == NULL) return NULL;

    ii = 0;
    for (i=0; i<k; i++) {
	if (select != NULL && !select[i]) continue;
	jj = 0;
	for (j=0; j<=i; j++) {
	    if (select != NULL && !select[j]) continue;
	    idx = ijton(i, j, pmod->ncoeff);
	    gretl_matrix_set(vcv, ii, jj, pmod->vcv[idx]);
	    if (jj != ii) {
		gretl_matrix_set(vcv, jj, ii, pmod->vcv[idx]);
	    }
	    jj++;
	}
	ii++;
    }

    return vcv;
}

/**
 * gretl_coeff_vector_from_model:
 * @pmod: pointer to model
 * @select: char array indicating which rows to select
 * (or %NULL for the full vector).
 *
 * Produces all or part of the coefficient vector for @pmod  
 * in the form of a gretl_vector.  Storage is allocated, to be freed
 * by the caller.  If @select is non-%NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows and zero elements otherwise.
 * 
 * Returns: the coefficient vector, or %NULL on error.
 * 
 */

gretl_vector *
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select)
{
    gretl_vector *b;
    int i, ii, nc;
    int k = pmod->ncoeff;

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }
    
    if (nc == 0) return NULL;

    b = gretl_column_vector_alloc(nc);
    if (b == NULL) return NULL;

    ii = 0;
    for (i=0; i<k; i++) {
	if (select != NULL && !select[i]) continue;
	gretl_vector_set(b, ii++, pmod->coeff[i]);
    }

    return b;
}

/**
 * gretl_is_identity_matrix:
 * @m: matrix to examine.
 *
 * Returns: 1 if @m is an identity matrix, 0 otherwise.
 */

int gretl_is_identity_matrix (const gretl_matrix *m)
{
    int i, j, idx;

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    if (m->packed) {
		if (i > j) continue;
		idx = packed_idx(m->rows, i, j);
	    } else {
		idx = mdx(m, i, j);
	    }
	    if (i == j && m->val[idx] != 1.0) return 0;
	    if (i != j && m->val[idx] != 0.0) return 0;
	}
    }

    return 1;
}

/**
 * gretl_is_zero_vector:
 * @v: vector to examine.
 *
 * Returns: 1 if @v is a zero vector, 0 otherwise.
 */

int gretl_is_zero_vector (const gretl_vector *v)
{
    int i, n = gretl_vector_get_length(v);

    for (i=0; i<n; i++) {
	if (v->val[i] != 0.0) return 0;
    }

    return 1;
}

#if 0

/* Some of the following may be useful later. */

int gretl_vector_normal_fill (gretl_vector *v, double *sigma)
{
    int i, n = gretl_vector_get_length(v);

    if (opt & OPT_N) {
	gretl_normal_dist(v->val, 0, n - 1);
    } else {
	gretl_uniform_dist(v->val, 0, n - 1);
    }

    if (sigma != 1.0) {
	for (i=0; i<n; i++) {
	    v->val *= sigma;
	}
    }

    return 0;
}

int gretl_matrix_set_column_value (gretl_matrix *m,
				   int column,
				   double val)
{
    for (i=0; i<m->rows; i++) {
	gretl_matrix_set(m, i, column, val);
    }

    return 0;
}

int gretl_matrix_set_column_values (gretl_matrix *m,
				    int column,
				    const gretl_vector *v)
{
    int i, n = gretl_vector_get_length(v);

    if (n != m->rows) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<m->rows; i++) {
	gretl_matrix_set(m, i, column, v->val[i]);
    }

    return 0;
}

int gretl_matrix_drop_initial_rows (gretl_matrix *targ,
				    const gretl_matrix *src,
				    int rows)
{
    int i, j;
    double srcv;

    if (rows >= src->rows) {
	return 1;
    }

    if (targ->rows != src->rows - rows ||
	targ->cols != src->cols) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<targ->rows; i++) {
	for (j=0; j<targ->cols; j++) {
	    srcv = gretl_matrix_get(src, i + rows, j);
	    gretk_matrix_set(targ, i, j, srcv);
	}
    }

    return 0;
}

int gretl_vector_drop_initial_vals (gretl_vector *targ,
				    const gretl_vector *src,
				    int drop)
{
    int i;
    int lsrc = gretl_vector_get_length(src);
    int ltarg = gretl_vector_get_length(targ);

    if (drop >= lsrc) {
	return 1;
    }

    if (ltarg != lsrc - drop) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<ltarg; i++) {
	targ->val[i] = src->val[i + drop];
    }    

    return 0;
}

#endif

