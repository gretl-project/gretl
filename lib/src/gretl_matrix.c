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

/* LAPACK-related matrix routines for gretl */

#include "libgretl.h"
#include "gretl_matrix.h"

#include "f2c.h"
#include "clapack_double.h"

static const char *wspace_fail = "gretl_matrix: workspace query failed\n";

#define gretl_is_vector(v) (v->rows == 1 || v->cols == 1)

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
 */

double gretl_matrix_get (const gretl_matrix *m, int i, int j)
{
    if (m == NULL || m->val == NULL ||
	i >= m->rows || j >= m->cols) {
	return NADBL;
    }

    return m->val[mdx(m, i, j)];
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
 */

double gretl_vector_get (const gretl_vector *v, int i)
{
    if (v == NULL || v->val == NULL ||
	(i >= v->rows && i >= v->cols)) {
	return NADBL;
    }

    return v->val[i];
}

/**
 * gretl_matrix_set:
 * @m: matrix to operate on.
 * @i: row index (zero-based).
 * @j: column index (zero-based).
 * @x: value to set.
 *
 * Sets element @i, @j of @m to the value @x, if the 
 * row and column are within bounds.
 *
 * Returns: %GRETL_MATRIX_OK, or %GRETL_MATRIX_ERR if the
 * indices are out of range for @m.
 */

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x)
{
    if (m == NULL || m->val == NULL ||
	i >= m->rows || j >= m->cols) {
	return GRETL_MATRIX_ERR;
    }

    m->val[mdx(m, i, j)] = x;

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
 */

int gretl_vector_set (gretl_vector *v, int i, double x)
{
    if (v == NULL || v->val == NULL ||
	(i >= v->rows && i >= v->cols)) {
	return GRETL_MATRIX_ERR;
    }

    v->val[i] = x;

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_alloc:
 * @rows: desired number of rows in matrix.
 * @cols: desired number of columns.
 *
 * Returns: pointer to a newly allocated gretl_matrix, or %NULL
 * on failure.  Note that the actual data storage is not
 * initialized.
 */

gretl_matrix *gretl_matrix_alloc (int rows, int cols)
{
    gretl_matrix *m = malloc(sizeof *m);

    if (m == NULL) {
	return m;
    }

    m->val = malloc(rows * cols * sizeof *m->val);
    if (m->val == NULL) {
	free(m);
	return NULL;
    }

    m->rows = rows;
    m->cols = cols;
    m->t = 0;

    return m;
}

/**
 * gretl_matrix_reuse:
 * @m: matrix to reuse.
 * @rows: desired number of rows in "new" matrix.
 * @cols: desired number of columns in "new" matrix.
 *
 * An "experts only" memory-conservation trick. If @m is an 
 * already-allocated gretl matrix, you can "resize" it by 
 * specifying a new number of rows and columns.  This works 
 * only if the product of @rows and @cols is less than or equal 
 * to the product of the number of rows and columns in the 
 * matrix as originally allocated; no actual reallocation of memory 
 * is performed.  If you "reuse" with an excessive number of rows
 * or columns you will surely crash your program or smash the
 * stack. Note also that the matrix-pointer returned is not really 
 * new, and the when the matrix is to be freed, gretl_matrix_free()
 * should be applied only once. 
 *
 * Returns: pointer to the "resized" gretl_matrix, or %NULL
 * if the product of @rows and @cols is out of bounds.
 */

gretl_matrix *gretl_matrix_reuse (gretl_matrix *m, int rows, int cols)
{
    m->rows = rows;
    m->cols = cols;

    return m;
}

/**
 * gretl_identity_matrix_new:
 * @n: desired number of rows and columns in the matrix.
 *
 * Returns: pointer to a newly allocated identity matrix, or %NULL
 * on failure.
 */

gretl_matrix *gretl_identity_matrix_new (int n)
{
    gretl_matrix *m;
    int i, j;

    m = gretl_matrix_alloc(n, n);

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		if (i == j) {
		    m->val[mdx(m, i, j)] = 1.0;
		} else {
		    m->val[mdx(m, i, j)] = 0.0;
		}
	    }
	}
    }

    return m;
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
 */

gretl_vector *
gretl_column_vector_from_array (const double *x, int n, GretlMatrixMod mod)
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

/**
 * gretl_data_series_to_vector:
 * @Z: data array.
 * @varno: ID number of variable
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Returns: a newly allocated gretl_vector containing the values
 * of the given variable (data series) for the given range,
 * or %NULL on failure.  
 */

gretl_vector *gretl_data_series_to_vector (const double **Z, int varno, 
					   int t1, int t2)
{
    gretl_matrix *v;
    int t, n = t2 - t1 + 1;

    if (n <= 0) {
	return NULL;
    }

    v = gretl_column_vector_alloc(n);
    if (v == NULL) {
	return NULL;
    }

    for (t=0; t<n; t++) {
	v->val[t] = Z[varno][t + t1];
    }

    return v;
}

static gretl_matrix *
gretl_matrix_copy_mod (const gretl_matrix *m, int mod)
{
    gretl_matrix *c;
    int rows, cols;
    int i, j, n;

    if (mod == GRETL_MOD_TRANSPOSE) {
	rows = m->cols;
	cols = m->rows;
    } else {
	rows = m->rows;
	cols = m->cols;
    }

    c = gretl_matrix_alloc(rows, cols);
    if (c == NULL) {
	return NULL;
    }

    if (mod == GRETL_MOD_TRANSPOSE) {
	for (i=0; i<c->rows; i++) {
	    for (j=0; j<c->cols; j++) {
		c->val[mdx(c, i, j)] = m->val[mdx(m, j, i)];
	    }
	}
    } else { 
	/* not transposing */
	n = rows * cols;
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
 */

gretl_matrix *gretl_matrix_copy (const gretl_matrix *m)
{
    return gretl_matrix_copy_mod(m, GRETL_MOD_NONE);
}

/**
 * gretl_matrix_copy_transpose:
 * @m: source matrix to be copied.
 *
 * Returns: an allocated copy of the tranpose of @m, or %NULL on failure.  
 */

gretl_matrix *gretl_matrix_copy_transpose (const gretl_matrix *m)
{
    return gretl_matrix_copy_mod(m, GRETL_MOD_TRANSPOSE);
}

/**
 * gretl_matrix_free:
 * @m: matrix to be freed.
 *
 * Frees the allocated storage in @m, then @m itself.
 */

void gretl_matrix_free (gretl_matrix *m)
{
    if (m == NULL) return;

    if (m->val != NULL) {
	free(m->val);
    }
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

    n = m->rows * m->cols;
    
    for (i=0; i<n; i++) {
	m->val[i] = 0.0;
    }
}

/**
 * gretl_matrix_log:
 * @m: input matrix.
 *
 * Sets all elements of @m to logs of original values.
 */

void gretl_matrix_log (gretl_matrix *m)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    n = m->rows * m->cols;
    
    for (i=0; i<n; i++) {
	m->val[i] = log(m->val[i]);
    }
}

/**
 * gretl_vector_mean:
 * @v: input vector.
 *
 * Returns: the arithmetic mean of the elements of @v, or
 * #NADBL on failure.
 */

double gretl_vector_mean (const gretl_vector *v)
{
    double ret = 0.0;
    int i, n, den;

    if (v == NULL || v->val == NULL) {
	return NADBL;
    }

    if (v->rows > 1 && v->cols > 1) {
	return NADBL;
    }
    
    if (v->rows > 1) {
	den = n = v->rows;
    } else {
	den = n = v->cols;
    }

    for (i=0; i<n; i++) {
	if (!na(v->val[i])) {
	    ret += v->val[i];
	} else {
	    den--;
	}
    }

    if (den > 0) {
	ret /= den;
    } else {
	ret = NADBL;
    }

    return ret;
}

/**
 * gretl_vector_variance:
 * @v: input vector.
 *
 * Returns: the variance of the elements of @v, or
 * #NADBL on failure.
 */

double gretl_vector_variance (const gretl_vector *v)
{
    double s2 = 0.0;
    double x, xbar;
    int i, n, den;

    if (v == NULL || v->val == NULL) {
	return NADBL;
    }

    den = n = gretl_vector_get_length(v);
    xbar = gretl_vector_mean(v);

    for (i=0; i<n; i++) {
	x = v->val[i];
	if (!na(x)) {
	    x -= xbar;
	    s2 += x * x;
	} else {
	    den--;
	}
    }

    if (den > 0) {
	s2 /= den;
    } else {
	s2 = NADBL;
    }

    return s2;
}

static int gretl_matrix_zero_triangle (gretl_matrix *m, char t)
{
    int i, j;

    if (m == NULL || m->val == NULL) 
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
 */

void gretl_matrix_multiply_by_scalar (gretl_matrix *m, double x)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    n = m->rows * m->cols;
    
    for (i=0; i<n; i++) {
	m->val[i] *= x;
    }
}

/**
 * gretl_matrix_divide_by_scalar:
 * @m: matrix to operate on.
 * @x: scalar by which to divide.
 *
 * Divides all elements of @m by @x.
 */

void gretl_matrix_divide_by_scalar (gretl_matrix *m, double x)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    n = m->rows * m->cols;
    
    for (i=0; i<n; i++) {
	m->val[i] /= x;
    }
}

/**
 * gretl_matrix_dot_pow:
 * @m: matrix to operate on.
 * @x: scalar to use for exponentiation.
 *
 * Raises all elements of @m to the power @x.
 */

void gretl_matrix_dot_pow (gretl_matrix *m, double x)
{
    int i, n;

    if (m == NULL || m->val == NULL) return;

    n = m->rows * m->cols;
    
    for (i=0; i<n; i++) {
	m->val[i] = pow(m->val[i], x);
    }
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
 */

int gretl_matrix_copy_values (gretl_matrix *targ, 
			      const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    n = src->rows * src->cols;
    
    for (i=0; i<n; i++) {
	targ->val[i] = src->val[i];
    }

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
 */

int 
gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    n = src->rows * src->cols;
    
    for (i=0; i<n; i++) {
	targ->val[i] += src->val[i];
    }

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
 */

int 
gretl_matrix_subtract_from (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    n = src->rows * src->cols;
    
    for (i=0; i<n; i++) {
	targ->val[i] -= src->val[i];
    }

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_transpose:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Fills out @targ (which must be pre-allocated and of the right
 * dimensions) with the transpose of @src.
 *
 * Returns: %GRETL_MATRIX_OK on success, non-zero error code otherwise.
 */

int gretl_matrix_transpose (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, j;
    double x;

    if (targ->rows != src->cols || targ->cols != src->rows) {
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<src->rows; i++) {
	for (j=0; j<src->cols; j++) {
	    x = src->val[mdx(src, i, j)];
	    targ->val[mdx(targ, j, i)] = x;
	}
    }

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

/**
 * gretl_matrix_print_to_prn:
 * @m: matrix to print.
 * @msg: accompanying message text (or %NULL if no message is wanted).
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn.
 */

void 
gretl_matrix_print_to_prn (const gretl_matrix *m, const char *msg, PRN *prn)
{
    int i, j;

    if (prn == NULL) {
	return;
    }

    if (msg != NULL && *msg != '\0') {
	pprintf(prn, "%s\n\n", msg);
    }

    if (m == NULL) {
	pputs(prn, " matrix is NULL\n");
    } else {
	for (i=0; i<m->rows; i++) {
	    for (j=0; j<m->cols; j++) {
		pprintf(prn, "%#12.5g ", gretl_matrix_get(m, i, j));
	    }
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    }
}

/**
 * gretl_matrix_print:
 * @m: matrix to print.
 * @msg: accompanying message text (or %NULL if no message is wanted).
 *
 * Prints the matrix @m to stderr.
 */

void gretl_matrix_print (const gretl_matrix *m, const char *msg)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);

    gretl_matrix_print_to_prn(m, msg, prn);
    gretl_print_destroy(prn);
}

#define EQTOL 1.0e-12

static int sneq (double x, double y)
{
    double reldiff;

    if (x == 0.0) {
	reldiff = fabs(y);
    } else if (y == 0.0) {
	reldiff = fabs(x);
    } else if (x > y) {
	reldiff = fabs((x - y) / y);
    } else {
	reldiff = fabs((y - x) / x);
    }

    if (reldiff > EQTOL) {
	fprintf(stderr, "relative difference = %g\n", reldiff);
    }

    return reldiff > EQTOL;
}

static int matrix_is_symmetric (const gretl_matrix *m)
{
    int i, j;

    for (i=1; i<m->rows; i++) {
	for (j=0; j<i; j++) {
	    if (sneq(m->val[mdx(m, i, j)], m->val[mdx(m, j, i)])) {
		double x = m->val[mdx(m, i, j)];
		double y = m->val[mdx(m, j, i)];

		fprintf(stderr, "M(%d,%d) = %.16g but M(%d,%d) = %.16g\n",
			i, j, x, j, i, y);
		gretl_matrix_print(m, "matrix_is_symmetric()");
		return 0;
	    }
	}
    }

    return 1;
}

static double gretl_matrix_one_norm (const gretl_matrix *m)
{
    double colsum, colmax = 0.0;
    int i, j;

    for (j=0; j<m->cols; j++) {
	colsum = 0.0;
	for (i=0; i<m->rows; i++) {
	    colsum += fabs(m->val[mdx(m, i, j)]);
	}
	if (colsum > colmax) {
	    colmax = colsum;
	}
    }

    return colmax;
}

/**
 * gretl_vcv_log_determinant:
 * @m: gretl_matrix.
 *
 * Compute the log determinant of the symmetric positive-definite
 * matrix @m using Cholesky decomposition.  
 * 
 * Returns: the log determinant, or #NABDL on failure.
 */

double gretl_vcv_log_determinant (const gretl_matrix *m)
{
    gretl_matrix *a = NULL;
    char uplo = 'U';
    integer info;
    integer n = m->rows;
    double det = NADBL;
    int i;

    if (m->rows != m->cols) {
	fputs("gretl_vcv_log_determinant: matrix must be square\n", stderr);
	return det;
    }

    if (!matrix_is_symmetric(m)) {
	fputs("gretl_vcv_log_determinant: matrix is not symmetric\n", stderr);
	return det;
    }

    a = gretl_matrix_copy(m);
    if (a == NULL) {
	fputs("gretl_vcv_log_determinant: out of memory\n", stderr);
	return det;
    }

    dpotrf_(&uplo, &n, a->val, &n, &info);

    if (info != 0) {
	if (info > 0) {
	    fputs("gretl_vcv_log_determinant: matrix not positive definite\n", 
		  stderr);
	} else {
	    fputs("gretl_vcv_log_determinant: illegal argument to dpotrf\n", 
		  stderr);
	}
    } else {
	det = 1.0;
	for (i=0; i<n; i++) {
	    det *= a->val[mdx(a, i, i)] * a->val[mdx(a, i, i)];
	}
	det = log(det);
    }

    gretl_matrix_free(a);

    return det;
}

/* calculate determinant using LU factorization.
   if logdet != 0, return the log of the determinant.
   if logdet != 0 and absval != 0, return the log of the
   absolute value of the determinant.  
*/   

static double gretl_LU_determinant (gretl_matrix *a, int logdet, int absval)
{
    integer info;
    integer n = a->rows;
    integer *ipiv;
    double det;
    int i;

    if (a->rows != a->cols) {
	fputs("gretl_LU_determinant: matrix must be square\n", stderr);
	return NADBL;
    }

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) return NADBL;

    dgetrf_(&n, &n, a->val, &n, ipiv, &info);

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
 * On exit, @b is replaced by the solution and @a is replaced by its 
 * LU decomposition.
 * 
 * Returns: 0 on successful completion, or a non-zero error code
 * (from the lapack function %dgetrs) on error.
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

static int 
matrix_multiply_self_transpose (const gretl_matrix *a,
				gretl_matrix *c)
{
    register int i, j, k;
    int nc = a->cols;
    int nr = a->rows;
    double targ;

    if (gretl_is_vector(a)) {
	fprintf(stderr, "matrix_multiply_self_transpose: got vector\n");
    }

    for (i=0; i<nc; i++) {
	for (j=i; j<nc; j++) {
	    targ = 0.0;
	    for (k=0; k<nr; k++) {
		targ += a->val[mdxtr(a,i,k)] * a->val[mdx(a,k,j)];
	    }
	    c->val[mdx(c,i,j)] = targ;
	    c->val[mdx(c,j,i)] = targ;
	}
    } 

    return GRETL_MATRIX_OK;
}

/**
 * gretl_matrix_multiply_mod:
 * @a: left-hand matrix.
 * @amod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @b: right-hand matrix.
 * @bmod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @c: matrix to hold the product.
 * 
 * Multiplies @a (or a-transpose) into @b (or b transpose),
 * with the result written into @c.
 *
 * Returns: %GRETL_MATRIX_OK on success; non-zero error code on
 * failure.
 * 
 */

int gretl_matrix_multiply_mod (const gretl_matrix *a, GretlMatrixMod amod,
			       const gretl_matrix *b, GretlMatrixMod bmod,
			       gretl_matrix *c)
{
    register int i, j, k;
    int lrows, lcols;
    int rrows, rcols;
    const int atr = (amod == GRETL_MOD_TRANSPOSE);
    const int btr = (bmod == GRETL_MOD_TRANSPOSE);
    int aidx, bidx;
    double targ;

    if (a == c || b == c) {
	fputs("gretl_matrix_multiply:\n product matrix must be "
	      "distinct from both input matrices\n", stderr);
	fprintf(stderr, "a = %p, b = %p, c = %p\n", 
		(void *) a, (void *) b, (void *) c);
	return GRETL_MATRIX_ERR;
    }

    if (a == b && atr && !btr && c->rows == a->cols && c->cols == a->cols) {
	return matrix_multiply_self_transpose(a, c);
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
 * gretl_matrix_kronecker_product:
 * @A: left-hand matrix.
 * @B: right-hand matrix.
 * 
 * Returns: A newly allocated matrix which is the Kronecker 
 * product of matrices @a and @b, or %NULL on failure.
 */

gretl_matrix *
gretl_matrix_kronecker_product (const gretl_matrix *A, const gretl_matrix *B)
{
    gretl_matrix *K;
    double aij, bkl;
    int p = A->rows;
    int q = A->cols;
    int r = B->rows;
    int s = B->cols;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;
    
    K = gretl_matrix_alloc(p * r, q * s);

    if (K != NULL) {
	for (i=0; i<p; i++) {
	    ioff = i * r;
	    for (j=0; j<q; j++) {
		/* block ij is an r * s matrix, a_{ij} * B */
		aij = A->val[mdx(A, i, j)];
		joff = j * s;
		for (k=0; k<r; k++) {
		    Ki = ioff + k;
		    for (l=0; l<s; l++) {
			bkl = B->val[mdx(B, k, l)];
			Kj = joff + l;
			K->val[mdx(K, Ki, Kj)] = aij * bkl;
		    }
		}
	    }
	}
    }

    return K;
}

/**
 * gretl_vector_dot_product:
 * @a: first vector.
 * @b: second vector.
 * @errp: pointer to receive error code (zero on success,
 * non-zero on failure), or %NULL.
 * 
 * Returns: The dot (scalar) product of @a and @b, or #NADBL on
 * failure.
 */

double gretl_vector_dot_product (const gretl_vector *a, const gretl_vector *b,
				 int *errp)
{
    int dima = (a->rows > 1)? a->rows : a->cols;
    int dimb = (b->rows > 1)? b->rows : b->cols;
    double dp = 0.0;
    int i;

    if (!gretl_is_vector(a) || !gretl_is_vector(b) || dima != dimb) {
	if (errp != NULL) {
	    *errp = GRETL_MATRIX_NON_CONFORM;
	}
	dp = NADBL;
    } else {
	for (i=0; i<dima; i++) {
	    dp += a->val[i] * b->val[i];
	}
    }

    return dp;
}

/**
 * gretl_matrix_dot_product:
 * @a: left-hand matrix.
 * @amod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @b: right-hand matrix.
 * @bmod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @errp: pointer to receive error code (zero on success,
 * non-zero on failure), or %NULL.
 * 
 * Returns: The dot (scalar) product of @a (or @a-transpose) and
 * @b (or @b-transpose), or #NADBL on failure.
 */

double gretl_matrix_dot_product (const gretl_matrix *a, GretlMatrixMod amod,
				 const gretl_matrix *b, GretlMatrixMod bmod,
				 int *errp)
{
    gretl_matrix *c = NULL;
    double ret = NADBL;
    int err = 0;

    if (gretl_is_vector(a) && gretl_is_vector(b)) {
	return gretl_vector_dot_product(a, b, errp);
    }

    c = gretl_matrix_alloc(1, 1);
    if (c == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(a, amod, b, bmod, c);
    }

    if (!err) {
	ret = c->val[0];
    }
	
    gretl_matrix_free(c);

    if (errp != NULL) {
	*errp = err;
    }

    return ret;
}

/**
 * gretl_matrix_dot_multiply:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * 
 * Returns: a new matrix, each of whose elements is the product of the
 * corresponding elements of the matrices @a and @b (or %NULL on
 * failure).
 */

gretl_matrix *gretl_matrix_dot_multiply (const gretl_matrix *a, 
					 const gretl_matrix *b)
{
    gretl_matrix *c;
    int i, n;

    if (a->rows != b->rows || a->cols != b->cols) {
	fputs("gretl_matrix_dot_multiply: matrices not conformable\n", stderr);
	return NULL;
    }

    c = gretl_matrix_alloc(a->rows, a->cols);
    if (c == NULL) {
	return NULL;
    }

    n = a->rows * a->cols;

    for (i=0; i<n; i++) {
	c->val[i] = a->val[i] * b->val[i];
    }

    return c;
}

/**
 * gretl_matrix_row_mean:
 * @m: source matrix.
 * @col: zero-based index of row.
 *
 * Returns: the mean of the elements in row @row of @m,
 * or #NADBL if the row is out of bounds.
 */

double gretl_matrix_row_mean (const gretl_matrix *m, int row)
{
    double sum = 0.0;
    int j;

    if (row >= m->rows) {
	return NADBL;
    }

    for (j=0; j<m->cols; j++) {
	sum += m->val[mdx(m, row, j)];
    }

    return sum / (double) m->cols;
}

/**
 * gretl_matrix_column_mean:
 * @m: source matrix.
 * @col: zero-based index of column.
 *
 * Returns: the mean of the elements in column @col of @m,
 * or #NADBL if the column is out of bounds.
 */

double gretl_matrix_column_mean (const gretl_matrix *m, int col)
{
    double sum = 0.0;
    int i;

    if (col >= m->cols) {
	return NADBL;
    }

    for (i=0; i<m->rows; i++) {
	sum += m->val[mdx(m, i, col)];
    }

    return sum / (double) m->rows;
}

/**
 * gretl_matrix_demean_by_row:
 * @m: matrix on which to operate.
 * 
 * For each row of @m, subtracts the row mean from each 
 * element on the row.
 */

void gretl_matrix_demean_by_row (gretl_matrix *m)
{
    double rowmean;  
    int i, j;

    for (i=0; i<m->rows; i++) {
	rowmean = gretl_matrix_row_mean(m, i);
	for (j=0; j<m->cols; j++) {
	    m->val[mdx(m, i, j)] -= rowmean;
	}
    }    
}

/**
 * gretl_matrix_demean_by_column:
 * @m: matrix on which to operate.
 * 
 * For each column of @m, subtracts the column mean from each 
 * element in the column.
 */

void gretl_matrix_demean_by_column (gretl_matrix *m)
{
    double colmean; 
    int i, j;

    for (j=0; j<m->cols; j++) {
	colmean = gretl_matrix_column_mean(m, j);
	for (i=0; i<m->rows; i++) {
	    m->val[mdx(m, i, j)] -= colmean;
	}
    }    
}

/**
 * gretl_matrix_vcv:
 * @m: source matrix (expected to have rows >= cols).
 *
 * Forms a variance-covariance matrix based on @m, thus:
 * (1) subtract the column means from the column elements of @m;
 * (2) multiply @m-transpose into @m; and
 * (3) divide the elements of the product by the number of rows
 *   in @m.
 * 
 * Returns: the allocated variance-covariance matrix, or %NULL
 * on failure.  Note that on return the column means have
 * been subtracted from @m.
 */

gretl_matrix *gretl_matrix_vcv (gretl_matrix *m)
{
    gretl_matrix *v;
    int err = 0;

    if (m->cols > m->rows) {
	fputs("gretl_matrix_vcv: expected rows >= cols\n", stderr);
	return NULL;
    }

    v = gretl_matrix_alloc(m->cols, m->cols);
    if (v == NULL) return NULL;

    gretl_matrix_demean_by_column(m);

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
 */

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c)
{
    return gretl_matrix_multiply_mod(a, GRETL_MOD_NONE,
				     b, GRETL_MOD_NONE,
				     c);
}

/**
 * gretl_symmetric_matrix_rcond:
 * @m: matrix to examine.
 * 
 * Estimates the reciprocal condition number of the real symmetric
 * positive definite matrix @m (in the 1-norm), using the lapack 
 * functions %dpotrf and %dpocon.
 *
 * Returns: the estimate, or #NADBL on failure to allocate memory.
 */

double gretl_symmetric_matrix_rcond (const gretl_matrix *m)
{
    gretl_matrix *a = NULL;
    char uplo = 'L';
    integer n = a->rows;
    integer lda = a->rows;
    integer info, *iwork = NULL;
    double *work = NULL;
    double anorm, rcond = NADBL;

    a = gretl_matrix_copy(m);
    work = malloc((3 * n) * sizeof *work);
    iwork = malloc(n * sizeof *iwork);

    if (a == NULL || work == NULL || iwork == NULL) {
	goto bailout;
    }

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_symmetric_matrix_rcond:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	rcond = 0.0;
    } else {
	anorm = gretl_matrix_one_norm(a);
	dpocon_(&uplo, &n, a->val, &lda, &anorm, &rcond, work, iwork, &info);
	if (info != 0) {
	    rcond = NADBL;
	}
    }

 bailout:

    free(work);
    free(iwork);
    gretl_matrix_free(a);

    return rcond;
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
 */

int gretl_matrix_cholesky_decomp (gretl_matrix *a)
{
    char uplo = 'L';
    integer n = a->rows;
    integer lda = a->rows;
    integer info;

    dpotrf_(&uplo, &n, a->val, &lda, &info);

    if (info != 0) {
	if (info > 0) {
	    fprintf(stderr, "n = %d, info = %d\n", (int) n, (int) info);
	    fputs("gretl_matrix_cholesky_decomp: matrix not positive definite\n", 
		  stderr);
	} else {
	    fputs("gretl_matrix_cholesky_decomp: illegal argument to dpotrf\n", 
		  stderr);
	}
    } 

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
 * gretl_invert_diagonal_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a diagonal matrix.
 * On exit @a is overwritten with the inverse. 
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_diagonal_matrix (gretl_matrix *a)
{
    double x;
    int i;

    if (a->cols != a->rows) {
	fputs("gretl_invert_diagonal_matrix: input is not square\n",
	      stderr);
	return GRETL_MATRIX_NON_CONFORM;
    }

    for (i=0; i<a->rows; i++) {
	x = a->val[mdx(a,i,i)];
	a->val[mdx(a,i,i)] = 1.0 / x;
    }

    return 0;
}

#if 0
static int spd_matrix_check_scaling (const gretl_matrix *a)
{
    integer n = a->rows;
    integer info;
    double *s;
    double scond, amax;

    s = malloc(n * sizeof *s);
    if (s == NULL) {
	return 1;
    }

    dpoequ_(&n, a->val, &n, s, &scond, &amax, &info);

    gretl_matrix_print(a, "a: check for scaling");
    fprintf(stderr, "info = %d, scond = %g, amax = %g\n",
	    (int) info, scond, amax);

    free(s);

    return 0;
}
#endif

/**
 * gretl_invert_symmetric_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization.  On exit @a is overwritten with 
 * the inverse. Uses the lapack functions %dpotrf and %dpotri.
 *
 * Returns: 0 on success; non-zero error code on failure.
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

    if (!matrix_is_symmetric(a)) {
	fputs("gretl_invert_symmetric_matrix: matrix is not symmetric\n", stderr);
	return GRETL_MATRIX_ERR;
    }

    n = a->cols;

    if (n == 1) {
	a->val[0] = 1.0 / a->val[0];
	return 0;
    }

#if 0
    spd_matrix_check_scaling(a);
#endif

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is not positive definite\n", stderr);
	}
	return GRETL_MATRIX_SINGULAR;
    } 

    dpotri_(&uplo, &n, a->val, &n, &info);

    if (info != 0) {
	err = GRETL_MATRIX_SINGULAR;
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotri failed with info = %d\n", (int) info);
    } else {
	gretl_symmetric_matrix_expand(a, uplo);
    }

    return err;
}

static int inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da < *db) - (*da > *db);
}

/**
 * gretl_eigen_sort:
 * @evals: array of eigenvalues.
 * @evecs: matrix of eigenvectors.
 * @rank: desired number of columns in output.
 * 
 * Sorts the real components of the eigenvalues in @evals from 
 * largest to smallest, and rearranges the columns in @evecs 
 * correspondingly.  If @rank is greater than zero and less than 
 * the number of columns in @evecs, then on output @evecs is shrunk 
 * so that it contains only the columns associated with the
 * largest @rank eigenvalues.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_eigen_sort (double *evals, gretl_matrix *evecs, int rank)
{
    struct esort {
	double vr;
	double vi;
	int idx;
    }; 
    struct esort *es;
    gretl_matrix *tmp;
    double x;
    int i, j, h, n;

    h = n = evecs->rows;
    if (rank > 0 && rank < h) {
	h = rank;
    }

    es = malloc(n * sizeof *es);
    if (es == NULL) {
	return E_ALLOC;
    }

    tmp = gretl_matrix_alloc(n, h);
    if (tmp == NULL) {
	free(es);
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	es[i].vr = evals[i];
	es[i].vi = evals[i + n];
	es[i].idx = i;
    }

    qsort(es, n, sizeof *es, inverse_compare_doubles);

    for (i=0; i<n; i++) {
	evals[i] = es[i].vr;
	evals[i + n] = es[i].vi;
    }

    for (j=0; j<h; j++) {
	for (i=0; i<n; i++) {
	    x = evecs->val[mdx(evecs, i, es[j].idx)];
	    tmp->val[mdx(tmp, i, j)] = x;
	}
    }

    free(evecs->val);
    evecs->val = tmp->val;
    evecs->cols = tmp->cols;
    tmp->val = NULL;
    free(tmp);

    free(es);

    return 0;
}

/**
 * gretl_general_matrix_eigenvals:
 * @m: square matrix on which to operate.
 * @ev: matrix in which to store the eigenvectors, or %NULL if 
 * the eigenvectors are not required.
 * 
 * Computes the eigenvalues of the general matrix @m.  If @ev is
 * non-%NULL, write the right eigenvectors of @m into @ev.
 * Uses the lapack function %dgeev.
 *
 * Returns: allocated storage containing the eigenvalues, or %NULL
 * on failure.  The returned array, on successful completion,
 * has 2n elements (where n = the number of rows and columns in the
 * matrix @m), the first n of which contain the real components of 
 * the eigenvalues of @m, and the remainder of which hold the 
 * imaginary components.
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

    if (m->rows != m->cols) {
	fprintf(stderr, "gretl_general_matrix_eigenvals:\n"
		" matrix must be square, is %d x %d\n", m->rows, m->cols);
	return NULL;
    }

    if (ev != NULL) {
	if (ev->rows != n || ev->cols != n) {
	    fprintf(stderr, "gretl_general_matrix_eigenvals:\n"
		    " matrix to hold eigenvalues should be %d x %d, is %d x %d\n",
		    m->rows, m->rows, ev->rows, ev->cols);
	    return NULL;
	}  
    }  

    work = malloc(sizeof *work);
    if (work == NULL) {
	return NULL;
    }

    wr = malloc(2 * n * sizeof *wr);
    if (wr == NULL) {
	goto bailout;
    } else {
	wi = wr + n;
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

    free(work);

    return wr;

 bailout:
    free(work);
    free(wr);

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

    if (!matrix_is_symmetric(m)) {
	fputs("gretl_symmetric_matrix_eigenvals: matrix is not symmetric\n", stderr);
	return NULL;
    }

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

/* return the row-index of the element in column col of
   matrix X that has the greatest absolute magnitude
*/

static int max_abs_index (const gretl_matrix *X, int col) 
{
    double aij, tmp = 0.0;
    int i, idx = 0;

    for (i=0; i<X->rows; i++) {
	aij = fabs(X->val[mdx(X, i, col)]);
	if (aij > tmp) {
	    tmp = aij;
	    idx = i;
	}
    }

    return idx;
}

/* given M = I - X'(XX')^{-1}X, normalize M and reduce to rank:
   return result in newly allocate matrix
*/

static gretl_matrix *find_base (gretl_matrix *M)
{
    gretl_matrix *C = NULL;
    char *keep = NULL;
    double TOL = 1.0e-15;
    double x, tmp;
    int idx, rank = 0;
    int i, j, k;

    keep = calloc(M->cols, 1);
    if (keep == NULL) {
	return NULL;
    }

    for (j=0; j<M->cols; j++) {
	idx = max_abs_index(M, j);
	tmp = M->val[mdx(M, idx, j)];

	if (fabs(tmp) > TOL) {
	    rank++;
	    keep[j] = 1;
	    /* normalize column j */
	    for (i=0; i<M->rows; i++) {
		M->val[mdx(M, i, j)] /= tmp;
	    }
	    /* normalize columns to the right of j */
	    for (k=j+1; k<M->cols; k++) {
		tmp = M->val[mdx(M, idx, k)];
		if (fabs(tmp) > TOL) {
		    for (i=0; i<M->rows; i++) {
			x = M->val[mdx(M, i, j)] * tmp;
			M->val[mdx(M, i, k)] -= x;
		    }
		}
	    }
	} 
    }

    if (rank > 0) {
	C = gretl_matrix_alloc(M->rows, rank);
    }

    if (C != NULL) {
	k = 0;
	for (j=0; j<M->cols && k<rank; j++) {
	    if (keep[j]) {
		for (i=0; i<M->rows; i++) {
		    C->val[mdx(C, i, k)] = M->val[mdx(M, i, j)];
		}
		k++;
	    }
	}
    }

    free(keep);

    return C;
}

/**
 * gretl_matrix_right_nullspace:
 * @M: matrix to operate on.
 * 
 * Given an m x n matrix @M, construct a conformable matrix
 * R such that MR = 0 (that is, all the columns of R are
 * orthogonal to the space spanned by the rows of @M).
 *
 * Returns: the allocated matrix R, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_right_nullspace (const gretl_matrix *M)
{
    gretl_matrix *A = NULL;
    gretl_matrix *B = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *R = NULL;
    int m = gretl_matrix_rows(M);
    int n = gretl_matrix_cols(M);
    double x;
    int i, j, err = 0;

    A = gretl_matrix_alloc(n, m);
    B = gretl_matrix_alloc(m, m);
    C = gretl_matrix_alloc(n, n);

    if (A == NULL || B == NULL || C == NULL) {
	err = 1;
    }

    if (!err) {
	/* B = MM' */
	err = gretl_matrix_multiply_mod(M, GRETL_MOD_NONE,
					M, GRETL_MOD_TRANSPOSE,
					B);
    }

    if (!err) {
	/* B = (MM')^{-1} */
	err = gretl_invert_symmetric_matrix(B);
    }

    if (!err) {
	/* A = M'(MM')^{-1} */
	err = gretl_matrix_multiply_mod(M, GRETL_MOD_TRANSPOSE,
					B, GRETL_MOD_NONE,
					A);
    }

    if (!err) {
	/* C = M'(MM')^{-1}M */
	err = gretl_matrix_multiply(A, M, C);
    } 

    if (!err) {
	/* make C = I - M'(MM')^{-1}M */
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		if (i == j) {
		    x = 1.0 - C->val[mdx(C, i, j)];
		} else {
		    x = - C->val[mdx(C, i, j)];
		}
		C->val[mdx(C, i, j)] = x;
	    }
	}
    }

    if (!err) {
	/* now make R into the normalized matrix */
	R = find_base(C);
    }

    gretl_matrix_free(A);
    gretl_matrix_free(B);
    gretl_matrix_free(C);
    
    if (err) {
	gretl_matrix_free(R);
	R = NULL;
    }

    return R;
}

/**
 * gretl_matrix_set_int:
 * @m: matrix to operate on.
 * @t: value to set
 * 
 * Sets an integer value on the gretl_matrix (used for internal
 * information).  
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
 */

int gretl_matrix_get_int (const gretl_matrix *m)
{
    return m->t;
}

static int
get_svd_ols_vcv (const gretl_matrix *A, const gretl_matrix *B,
		 const double *s, gretl_matrix *vcv, double *s2)
{
    double aik, ajk, vij;
    int m = A->cols;
    int i, j, k;

    /* Get X'X{-1}, based on the work done by the SV decomp:
       reciprocals of the squares of the (positive) singular values,
       premultiplied by V and postmultiplied by V-transpose
    */

    for (i=0; i<m; i++) {
	for (j=i; j<m; j++) {
	    vij = 0.0;
	    for (k=0; k<m; k++) {
		if (s[k] > 0.0) {
		    aik = A->val[mdx(A, k, i)];
		    ajk = A->val[mdx(A, k, j)];
		    vij += aik * ajk / (s[k] * s[k]);
		}
	    }
	    vcv->val[mdx(vcv, i, j)] = vij;
	    if (j != i) {
		vcv->val[mdx(vcv, j, i)] = vij;
	    }
	}
    }

    if (s2 != NULL) {
	double sigma2 = 0.0;
	int T = A->rows;

	for (i=m; i<T; i++) {
	    sigma2 += B->val[i] * B->val[i];
	}
	sigma2 /= T - m;
	gretl_matrix_multiply_by_scalar(vcv, sigma2);  
	*s2 = sigma2;
    }

    return 0;
}

static int
get_ols_vcv (const gretl_vector *y, const gretl_matrix *X,
	     const gretl_vector *b, gretl_matrix *vcv,
	     double *s2)
{
    if (gretl_invert_general_matrix(vcv)) {
	gretl_matrix_print(vcv, "vcv: inversion failed");
	return 1;
    }

    if (s2 != NULL) {
	double u, sigma2 = 0.0;
	int k = X->cols;
	int n = X->rows;
	int i, j;

	for (i=0; i<n; i++) {
	    u = y->val[i];
	    for (j=0; j<k; j++) {
		u -= X->val[mdx(X, i, j)] * b->val[j];
	    }
	    sigma2 += u * u;
	}

	sigma2 /= (n - k);

	gretl_matrix_multiply_by_scalar(vcv, sigma2);  

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

/**
 * gretl_matrix_svd_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or %NULL if this is not needed.
 * @uhat: vector to hold the regression residuals, or %NULL if 
 * these are not needed.
 * @s2: pointer to receive residual variance, or %NULL.  Note:
 * if @s2 is %NULL, the vcv estimate will be plain (X'X)^{-1}.
 *
 * Computes OLS estimates using SVD decomposition, and puts the
 * coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv and the residuals in @uhat.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

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
	fprintf(stderr, "gretl_matrix_svd_ols:\n"
		" dgelss: rank of data matrix X = %d (rows = %d, cols = %d)\n", 
		(int) rank, X->rows, X->cols);
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
 * gretl_SVD_invert_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a general matrix using SVD
 * factorization.  Uses the lapack function %dgesvd.
 * Prints to stderr an estimate of the reciprocal
 * condition number of the matrix if it is less then 1e-9.
 * On exit the original matrix is overwritten by the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_SVD_invert_matrix (gretl_matrix *a)
{
    gretl_matrix *u = NULL;
    gretl_matrix *vt = NULL;
    integer n = a->rows;
    char jobu = 'A';
    char jobvt = 'A';
    integer lwork = -1;
    integer info;

    double *work2, *work = NULL;
    double *s = NULL;

    int i, j, err = GRETL_MATRIX_OK;

    if (a->rows != a->cols) {
	err = GRETL_MATRIX_NON_CONFORM;
	goto bailout;
    }	

    /* a = USV' ; a^{-1} = VWU' where W holds inverse of diag elements of S */

    s = malloc(n * sizeof *s);
    u = gretl_matrix_alloc(n, n);
    vt = gretl_matrix_alloc(n, n);
    work = malloc(sizeof *work);

    if (s == NULL || u == NULL || vt == NULL || work == NULL) {
	err = GRETL_MATRIX_NOMEM; 
	goto bailout;
    }	

    /* workspace query */
    dgesvd_(&jobu, &jobvt, &n, &n, a->val, &n, s, u->val, &n, vt->val, &n, 
	    work, &lwork, &info);

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

    /* actual computation */
    dgesvd_(&jobu, &jobvt, &n, &n, a->val, &n, s, u->val, &n, vt->val, &n, 
	    work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "gretl_matrix_SVD_inverse: info = %d\n", (int) info);
	err = GRETL_MATRIX_ERR;
    }

    if (!err) {
	/* invert singular values */
	for (i=0; i<n; i++) {
	    double wi = 0.0;

	    if (s[i] != 0.0) {
		wi = 1.0 / s[i];
	    }
	    for (j=0; j<n; j++) {
		u->val[mdx(u, j, i)] *= wi;
	    }
	}
	if (s[n-1] == 0.0) {
	    fputs("gretl_SVD_invert_matrix: matrix is not of full rank\n", stderr);
	} else if (s[n-1] / s[0] < 1e-9) {
	    fprintf(stderr, "rcond = %g\n", s[n-1] / s[0]);
	}	
	gretl_matrix_multiply_mod(vt, GRETL_MOD_TRANSPOSE,
				  u, GRETL_MOD_TRANSPOSE,
				  a);
    }

 bailout:
    
    free(s);
    free(work);
    gretl_matrix_free(u);
    gretl_matrix_free(vt);

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
 * @s2: pointer to receive residual variance, or %NULL.  Note:
 * if @s2 is %NULL, the vcv estimate will be plain (X'X)^{-1}.
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
 * gretl_matrix_restricted_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @R: left-hand restriction matrix, as in Rb = q.
 * @q: right-hand restriction vector.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or %NULL if this is not needed.
 * @s2: pointer ro receive residual variance, or NULL.  If vcv is non-NULL
 * and s2 is NULL, the vcv estimate is just W^{-1}.
 *
 * Computes OLS estimates restricted by R and q, using LU factorization, 
 * and puts the coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 * 
 */

int 
gretl_matrix_restricted_ols (const gretl_vector *y, const gretl_matrix *X,
			     const gretl_matrix *R, const gretl_vector *q,
			     gretl_vector *b, gretl_matrix *vcv,
			     double *s2)
{
    gretl_matrix *XTX = NULL;
    gretl_vector *V = NULL;
    gretl_matrix *W = NULL;
    int k = X->cols;
    int nr = R->rows;
    int ldW = k + nr;
    int err = GRETL_MATRIX_OK;
    int i, j;

    if (gretl_vector_get_length(b) != k) {
	err = GRETL_MATRIX_NON_CONFORM;
    }

    if (!err) {
	XTX = gretl_matrix_alloc(k, k);
	V = gretl_column_vector_alloc(ldW);
	W = gretl_matrix_alloc(ldW, ldW);
	if (XTX == NULL || V == NULL || W == NULL) {
	    err = GRETL_MATRIX_NOMEM;
	}
    }

    /* construct V matrix: X'y augmented by q (or by a 0
       matrix if q is NULL)
    */

    if (!err) {
	V->rows = k;
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					y, GRETL_MOD_NONE,
					V);
	V->rows = ldW;
    }

    if (!err) {
	for (i=k; i<ldW; i++) {
	    if (q != NULL) {
		V->val[i] = q->val[i-k];
	    } else {
		V->val[i] = 0.0;
	    }
	}
    }

    /* construct W matrix: X'X augmented by R and R' */

    if (!err) {
	gretl_matrix_zero(W);
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					X, GRETL_MOD_NONE,
					XTX);
    }

    if (!err) {
	for (i=0; i<XTX->rows; i++) {
	    for (j=0; j<XTX->cols; j++) {
		W->val[mdx(W,i,j)] = XTX->val[mdx(XTX,i,j)];
	    }
	}
	for (i=0; i<R->rows; i++) {
	    for (j=0; j<R->cols; j++) {
		W->val[mdx(W,i+k,j)] = R->val[mdx(R,i,j)];
	    }
	}
	for (i=0; i<R->cols; i++) {
	    for (j=0; j<R->rows; j++) {
		W->val[mdx(W,i,j+k)] = R->val[mdx(R,j,i)];
	    }
	}
    }

    if (!err && vcv != NULL) {
	err = gretl_matrix_copy_values(vcv, W);
    }

    if (!err) {
	err = gretl_LU_solve(W, V);
    }

    if (!err) {
	int i;
	
	for (i=0; i<k; i++) {
	    b->val[i] = V->val[i];
	}
	if (vcv != NULL) {
	    err = get_ols_vcv(y, X, b, vcv, s2);
	}
    }

    if (XTX != NULL) gretl_matrix_free(XTX);
    if (V != NULL) gretl_vector_free(V);
    if (W != NULL) gretl_matrix_free(W);

    return err;
}

/**
 * gretl_scalar_b_X_b:
 * @b: k-vector.
 * @bmod: either %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @X: k x k matrix.
 * @errp: pointer to receive error code, or %NULL.
 *
 * If @bmod = %GRETL_MOD_NONE, computes the scalar product
 * b * X * b'; if @bmod = %GRETL_MOD_TRANSPOSE, computes
 * the scalar product b' * X * b.
 * If @errp is not %NULL its content is set to 0 on success, non-zero
 * on failure.
 * 
 * Returns: the scalar product, or #NADBL on failure.
 */

double gretl_scalar_b_X_b (const gretl_vector *b, 
			   GretlMatrixMod bmod,
			   const gretl_matrix *X,
			   int *errp)
{
    gretl_matrix *tmp = NULL;
    double ret = NADBL;
    int tmpdim = (bmod == GRETL_MOD_TRANSPOSE)?
	b->rows : b->cols;
    int chk = (bmod == GRETL_MOD_TRANSPOSE)? 
	b->cols : b->rows;
    int err = 0;

    if (X->rows != X->cols || tmpdim != X->rows || chk != 1) {
	err = GRETL_MATRIX_NON_CONFORM;
    }

    if (!err) {
	tmp = gretl_matrix_alloc(1, tmpdim);
	if (tmp == NULL) {
	    err = GRETL_MATRIX_NOMEM;
	}
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(b, bmod,
					X, GRETL_MOD_NONE,
					tmp);
    }

    if (!err) {
	GretlMatrixMod rmod = (bmod == GRETL_MOD_TRANSPOSE)?
	    GRETL_MOD_NONE : GRETL_MOD_TRANSPOSE;

	ret = gretl_matrix_dot_product(tmp, GRETL_MOD_NONE,
				       b, rmod, &err);
    }

    gretl_matrix_free(tmp);

    if (err) {
	ret = NADBL;
    }

    if (errp != NULL) {
	*errp = err;
    }

    return ret;
}

#define mod_complement(m) ((m)? GRETL_MOD_NONE : GRETL_MOD_TRANSPOSE)

/**
 * gretl_matrix_A_X_A:
 * @A: m * k matrix or k * m matrix, depending on @amod.
 * @amod: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE: in the first
 * case @A should be m * k; in the second, k * m;
 * @X: k * k matrix.
 * @errp: pointer to receive error code, or %NULL.
 *
 * Computes either A * X * A' (if amod = %GRETL_MOD_NONE) or
 * A' * X * A (if amod = %GRETL_MOD_TRANSPOSE).
 * If @errp is not %NULL its content is set to 0 on success, non-zero
 * on failure.
 * 
 * Returns: allocated m * m matrix product, or %NULL on error.
 */

gretl_matrix *
gretl_matrix_A_X_A (const gretl_matrix *A, GretlMatrixMod amod,
		    const gretl_matrix *X, int *errp)
{
    gretl_matrix *tmp = NULL;
    gretl_matrix *ret = NULL;
    int m = (amod)? A->cols : A->rows;
    int k = (amod)? A->rows : A->cols;
    int err = 0;

    if (errp != NULL) {
	*errp = 0;
    }

    if (X->rows != k || X->cols != k) {
	if (errp != NULL) {
	    *errp = GRETL_MATRIX_NON_CONFORM;
	}
	return NULL;
    }

    tmp = gretl_matrix_alloc(m, k);
    ret = gretl_matrix_alloc(m, m);

    if (tmp == NULL || ret == NULL) {
	gretl_matrix_free(tmp);
	gretl_matrix_free(ret);
	if (errp != NULL) {
	    *errp = GRETL_MATRIX_NOMEM;
	}
	return NULL;
    }

    err = gretl_matrix_multiply_mod(A, amod,
				    X, GRETL_MOD_NONE,
				    tmp);

    if (!err) {
	err = gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
					A, mod_complement(amod),
					ret);
    }

    gretl_matrix_free(tmp);

    if (err) {
	gretl_matrix_free(ret);
	ret = NULL;
	if (errp != NULL) {
	    *errp = err;
	}
    }

    return ret;
}

/**
 * gretl_matrix_diagonal_sandwich:
 * @d: k-vector.
 * @X: k * k matrix.
 * @DXD: target k * k matrix.
 *
 * Computes in @DXD (which must be pre-allocated), the matrix
 * product D * X * D, where D is a diagonal matrix with elements
 * given by the vector @d.
 * 
 * Returns: 0 on success, non-zero code on error.
 */

int
gretl_matrix_diagonal_sandwich (const gretl_vector *d, const gretl_matrix *X,
				gretl_matrix *DXD)
{
    int dim = (d->rows == 1)? d->cols : d->rows;
    double xij;
    int i, j, err = 0;

    if (dim != X->rows || dim != X->cols ||
	dim != DXD->rows || dim != DXD->cols) {
	err = GRETL_MATRIX_NON_CONFORM;
    } else {
	for (i=0; i<dim; i++) {
	    for (j=0; j<dim; j++) {
		xij = X->val[mdx(X, i, j)];
		DXD->val[mdx(DXD, i, j)] = xij * d->val[i] * d->val[j];
	    }
	}
    }

    return err;
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
    if (makevcv(pmod)) {
	return NULL;
    }

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }
    
    if (nc == 0) {
	return NULL;
    }

    vcv = gretl_matrix_alloc(nc, nc);
    if (vcv == NULL) {
	return NULL;
    }

    ii = 0;
    for (i=0; i<k; i++) {
	if (select != NULL && !select[i]) {
	    continue;
	}
	jj = 0;
	for (j=0; j<=i; j++) {
	    if (select != NULL && !select[j]) {
		continue;
	    }
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
 * in the form of a gretl column vector.  Storage is allocated, to be freed
 * by the caller.  If @select is non-%NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows and zero elements otherwise.
 * 
 * Returns: the coefficient vector, or %NULL on error.
 */

gretl_vector *
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select)
{
    gretl_vector *b;
    int i, j, nc;
    int k = pmod->ncoeff;

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }
    
    if (nc == 0) {
	return NULL;
    }

    b = gretl_column_vector_alloc(nc);
    if (b == NULL) {
	return NULL;
    }

    j = 0;
    for (i=0; i<k; i++) {
	if (select != NULL && !select[i]) {
	    continue;
	}
	b->val[j++] = pmod->coeff[i];
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
	    idx = mdx(m, i, j);
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

/**
 * gretl_covariance_matrix_from_varlist:
 * @list: list of variables by ID number.
 * @Z: data array.
 * @pdinfo: pointer to data information struct.
 * @means: pointer to pick up vector of means, or %NULL to discard.
 * @errp: pointer to receive non-zero error code in case of
 * failure, or %NULL.
 *
 * Returns: the variance-covariance matrix of the listed variables
 * (over the currently defined data sample), or %NULL in case of
 * failure.
 */

gretl_matrix *
gretl_covariance_matrix_from_varlist (const int *list, const double **Z, 
				      const DATAINFO *pdinfo, 
				      gretl_matrix **means,
				      int *errp)
{
    gretl_matrix *vcv;
    gretl_vector *xbar;

    int k = list[0];
    int t, i, j;

    double vv, x, y;
    double xbi, xbj;
    int nv, err = 0;
    
    if (errp != NULL) {
	*errp = 0;
    }

    vcv = gretl_matrix_alloc(k, k);
    if (vcv == NULL) {
	if (errp != NULL) {
	    *errp = E_ALLOC;
	}
	return NULL;
    }

    xbar = gretl_vector_alloc(k);
    if (xbar == NULL) {
	if (errp != NULL) {
	    *errp = E_ALLOC;
	}
	gretl_matrix_free(vcv);
	return NULL;
    }
    
    for (i=0; i<k && !err; i++) {
	xbi = gretl_mean(pdinfo->t1, pdinfo->t2, Z[list[i+1]]);
	if (na(xbi)) {
	    err = E_DATA;
	} else {
	    gretl_vector_set(xbar, i, xbi);
	}
    }

    for (i=0; i<k && !err; i++) {
	xbi = gretl_vector_get(xbar, i);
	for (j=i; j<k; j++) {
	    xbj = gretl_vector_get(xbar, j);
	    vv = 0.0;
	    nv = 0;
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		x = Z[list[i+1]][t];
		y = Z[list[j+1]][t];
		if (na(x) || na(y)) {
		    continue;
		}
		vv += (x - xbi) * (y - xbj);
		nv++;
	    }
	    if (nv < 2) {
		err = E_DATA;
		vv = NADBL;
	    } else {
		vv /= (nv - 1); /* plain nv? */
	    }
	    gretl_matrix_set(vcv, i, j, vv);
	    gretl_matrix_set(vcv, j, i, vv);
	}
    }

    if (means != NULL && !err) {
	*means = xbar;
    } else {
	gretl_vector_free(xbar);
    }

    if (err) {
	gretl_matrix_free(vcv);
	vcv = NULL;
	if (errp != NULL) {
	    *errp = err;
	}
    }

    return vcv;
}

/**
 * gretl_matrix_row_to_array:
 * @m: source matrix.
 * @i: the row from which values should be copied.
 * @x: array of doubles large enough to hold a row from @m.
 *
 * Copies the values from row @i of matrix @m into the array
 * @x, which should already be allocated to the correct size.
 *
 * Returns: 0 on sucess, 1 if the row is out of bounds.
 */

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x)
{
    int j, err = 0;

    if (i >= gretl_matrix_rows(m)) {
	err = 1;
    } else {
	for (j=0; j<m->cols; j++) {
	    x[j] = m->val[mdx(m, i, j)];
	}
    }

    return err;
}

/**
 * gretl_matrix_array_alloc:
 * @n: number of matrices.
 *
 * Allocates an array of @n gretl matrix pointers. On successful
 * allocation of the array, each element is initialized to %NULL.
 *
 * Returns: pointer on sucess, %NULL on failure.
 */

gretl_matrix **gretl_matrix_array_alloc (int n)
{
    gretl_matrix **A = malloc(n * sizeof *A);
    int i;

    if (A != NULL) {
	for (i=0; i<n; i++) {
	    A[i] = NULL;
	}
    }

    return A;
}

/**
 * gretl_matrix_array_alloc_with_size:
 * @n: number of matrices.
 * @rows: number of rows in each matrix.
 * @cols: number of columns in each matrix.
 *
 * Allocates an array of @n gretl matrix pointers, each one
 * with size @rows * @cols.
 *
 * Returns: pointer on sucess, %NULL on failure.
 */

gretl_matrix **
gretl_matrix_array_alloc_with_size (int n, int rows, int cols)
{
    gretl_matrix **A = malloc(n * sizeof *A);
    int i, j;

    if (A != NULL) {
	for (i=0; i<n; i++) {
	    A[i] = gretl_matrix_alloc(rows, cols);
	    if (A[i] == NULL) {
		for (j=0; j<i; j++) {
		    gretl_matrix_free(A[i]);
		}
		free(A);
		A = NULL;
		break;
	    }
	}
    }

    return A;
}


/**
 * gretl_matrix_array_free:
 * @A: dyamically allocated array of gretl matrices.
 * @n: number of matrices in array.
 *
 * Frees each of the @n gretl matrices in the array @A, and
 * the array itself.  See also gretl_matrix_array_alloc().
 */

void gretl_matrix_array_free (gretl_matrix **A, int n)
{
    int i;

    if (A != NULL) {
	for (i=0; i<n; i++) {
	    gretl_matrix_free(A[i]);
	}
	free(A);
    }
}

/**
 * gretl_matrix_data_subset:
 * @list: list of variable to process.
 * @Z: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @Z specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns.
 *
 * Returns: allocated matrix or %NULL on failure. 
 */

gretl_matrix *gretl_matrix_data_subset (const int *list, const double **Z,
					int t1, int t2)
{
    gretl_matrix *M;
    int T = t2 - t1 + 1;
    int n = list[0];
    int j, t, v;

    if (T <= 0 || n <= 0) {
	return NULL;
    }

    M = gretl_matrix_alloc(T, n);
    if (M == NULL) {
	return NULL;
    }

    for (j=0; j<n; j++) {
	v = list[j + 1];
	for (t=0; t<T; t++) {
	    gretl_matrix_set(M, t, j, Z[v][t + t1]);
	}
    }

    return M;
}
