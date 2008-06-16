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

/* Matrix routines for gretl, several of which are LAPACK-related */

#include "libgretl.h"
#include "gretl_matrix.h"

#include <errno.h>
#include <assert.h>

#include "gretl_f2c.h"
#include "clapack_double.h"
#include "../../cephes/libprob.h"

#define gretl_is_vector(v) (v->rows == 1 || v->cols == 1)
#define matrix_is_scalar(m) (m->rows == 1 && m->cols == 1)

#define mdx(a,i,j)   ((j)*a->rows+(i))

#define matrix_transp_get(m,i,j) (m->val[(i)*m->rows+(j)])
#define matrix_transp_set(m,i,j,x) (m->val[(i)*m->rows+(j)]=x)

#define SVD_SMIN 1.0e-9

static int add_scalar_to_matrix (gretl_matrix *targ, double x);

/* Central accounting for error in matrix allocation */

static int gretl_matrix_err;

/* get, and clear, the matrix error code */

int get_gretl_matrix_err (void)
{
    int ret = gretl_matrix_err;
    gretl_matrix_err = 0;
    return ret;
}

void clear_gretl_matrix_err (void)
{
    gretl_matrix_err = 0;
}

static void set_gretl_matrix_err (int err)
{
    if (gretl_matrix_err == 0) {
	gretl_matrix_err = err;
    }
}

static int wspace_fail (integer info, double w0)
{
    int iinfo = (int) info;

    fprintf(stderr, "gretl_matrix: workspace query failed: info = %d, "
	    "work[0] = %g\n", iinfo, w0);
    return E_DATA;
}

/* An efficient means of allocating temporary storage for lapack
   operations: this should be used _only_ for temporary allocations
   that would ordinarily be freed before returning from the function
   in question.  In this mode we keep the chunk around for future use,
   expanding it as needed.
*/

static void *lapack_mem_chunk;
static size_t lapack_mem_sz;

static void *lapack_malloc (size_t sz)
{
    void *mem = NULL;

    if (sz > lapack_mem_sz) {
	void *chunk = realloc(lapack_mem_chunk, sz);

	if (chunk != NULL) {
	    lapack_mem_chunk = mem = chunk;
	    lapack_mem_sz = sz;
	}
    } else {
	mem = lapack_mem_chunk;
    }

    return mem;
}

static void *lapack_realloc (void *p, size_t sz)
{
    return lapack_malloc(sz);
}

void lapack_mem_free (void)
{
    free(lapack_mem_chunk);
    lapack_mem_chunk = NULL;
    lapack_mem_sz = 0;
}

void lapack_free (void *p)
{
    return;
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
    gretl_matrix *m;

    if (rows <= 0 || cols <= 0) {
	return NULL;
    }

    m = malloc(sizeof *m);
    if (m == NULL) {
	set_gretl_matrix_err(E_ALLOC);
	return m;
    }

    m->val = malloc(rows * cols * sizeof *m->val);
    if (m->val == NULL) {
	free(m);
	return NULL;
    }

    m->rows = rows;
    m->cols = cols;
    m->t1 = m->t2 = 0;

    return m;
}

int gretl_matrix_xna_check (const gretl_matrix *m)
{
    if (m == NULL) {
	return E_DATA;
    } else {
	int i, n = m->rows * m->cols;

	for (i=0; i<n; i++) {
	    if (xna(m->val[i])) {
		gretl_errmsg_set(_("Matrix is not finite"));
		return E_NAN;
	    }
	}
    }

    return 0;
}

int gretl_matrix_get_structure (const gretl_matrix *m)
{
    int ret = 0;

    if (gretl_is_null_matrix(m)) {
	return 0;
    }

    if (m != NULL) {
	if (m->rows == m->cols) {
	    ret = GRETL_MATRIX_SQUARE;
	    if (m->rows == 1) {
		ret = GRETL_MATRIX_SCALAR;
	    }
	} 
    }

    if (ret == GRETL_MATRIX_SQUARE) {
	double x;
	int uzero = 1;
	int lzero = 1;
	int symm = 1;
	int udiag = 1;
	int i, j;
	

	for (i=0; i<m->rows; i++) {
	    for (j=0; j<m->cols; j++) {
		x = gretl_matrix_get(m,i,j);
		if (j > i) {
		    if (x != 0.0) {
			uzero = 0;
		    }
		} else if (i > j) {
		    if (x != 0.0) {
			lzero = 0;
		    }
		} else if (i == j) {
		    if (x != 1.0) {
			udiag = 0;
		    }
		}
		if (j != i && x != gretl_matrix_get(m,j,i)) {
		    symm = 0;
		}
		if (!uzero && !lzero && !symm) {
		    break;
		}
	    }
	    if (!uzero && !lzero && !symm) {
		break;
	    }
	}

	if (udiag && uzero && lzero) {
	    ret = GRETL_MATRIX_IDENTITY;
	} else if (uzero && lzero) {
	    ret = GRETL_MATRIX_DIAGONAL;
	} else if (uzero) {
	    ret = GRETL_MATRIX_LOWER_TRIANGULAR;
	} else if (lzero) {
	    ret = GRETL_MATRIX_UPPER_TRIANGULAR;
	} else if (symm) {
	    ret = GRETL_MATRIX_SYMMETRIC;
	}
    }

    return ret;
}

/**
 * gretl_matrix_reuse:
 * @m: matrix to reuse.
 * @rows: desired number of rows in "new" matrix, or -1
 * to leave the current value unchanged.
 * @cols: desired number of columns in "new" matrix, or -1
 * to leave the current value unchanged.
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
 * Returns: pointer to the "resized" gretl_matrix.
 */

gretl_matrix *gretl_matrix_reuse (gretl_matrix *m, int rows, int cols)
{
    if (rows > 0) {
	m->rows = rows;
    }

    if (cols > 0) {
	m->cols = cols;
    }

    return m;
}

/**
 * gretl_matrix_realloc:
 * @m: matrix to reallocate.
 * @rows: desired number of rows in "new" matrix.
 * @cols: desired number of columns in "new" matrix.
 *
 * Reallocates the storage in @m to the specified dimensions.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int gretl_matrix_realloc (gretl_matrix *m, int rows, int cols)
{
    int n = rows * cols;
    double *x;

    if (rows == m->rows && cols == m->cols) {
	/* no-op */
	return 0;
    }

    x = realloc(m->val, n * sizeof *m->val);
    if (x == NULL) {
	return E_ALLOC;
    }

    m->val = x;
    m->rows = rows;
    m->cols = cols;

    return 0;
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
    int i, k;

    if (n < 0) {
	return NULL;
    } else if (n == 0) {
	return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(n, n);

    if (m != NULL) {
	k = n * n;
	n++;
	for (i=0; i<k; i++) {
	    m->val[i] = (i % n)? 0.0 : 1.0;
	}
    }


    return m;
}

/**
 * gretl_zero_matrix_new:
 * @r: desired number of rows in the matrix.
 * @c: desired number of columns in the matrix.
 *
 * Returns: pointer to a newly allocated zero matrix, or %NULL
 * on failure.
 */

gretl_matrix *gretl_zero_matrix_new (int r, int c)
{
    gretl_matrix *m;
    int i, n = r * c;

    if (r < 0 || c < 0) {
	return NULL;
    } else if (r == 0 || c == 0) {
	return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(r, c);

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    m->val[i] = 0.0;
	}
    }

    return m;
}

/**
 * gretl_unit_matrix_new:
 * @r: desired number of rows in the matrix.
 * @c: desired number of columns in the matrix.
 *
 * Returns: pointer to a newly allocated matrix, all
 * of whose elements equal 1, or %NULL on failure.
 */

gretl_matrix *gretl_unit_matrix_new (int r, int c)
{
    gretl_matrix *m;
    int i, n = r * c;

    if (r < 0 || c < 0) {
	return NULL;
    } else if (r == 0 || c == 0) {
	return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(r, c);

    if (m != NULL) {
	for (i=0; i<n; i++) {
	    m->val[i] = 1.0;
	}
    }

    return m;
}

/**
 * gretl_null_matrix_new:
 *
 * Returns: pointer to a newly allocated null matrix 
 * (for use in declaration of a variable as a matrix),
 * %NULL on failure.
 */

gretl_matrix *gretl_null_matrix_new (void)
{
    gretl_matrix *m = malloc(sizeof *m);

    if (m == NULL) {
	set_gretl_matrix_err(E_ALLOC);
	return m;
    }

    m->val = NULL;
    m->rows = m->cols = 0;
    m->t1 = m->t2 = 0;

    return m;
}

/**
 * gretl_matrix_seq:
 *
 * @start: first element.
 * @end: last element.
 *
 * Returns: pointer to a row vector, containing the numbers from
 * @start to @end, in decreasing order if @start > @end --
 * or %NULL on failure.
 */

gretl_matrix *gretl_matrix_seq (int start, int end)
{
    int reverse = (start > end);
    int i, k, n = 1 + (reverse ? (start-end) : (end-start));

    if (n == 0) {
	return NULL;
    }

    gretl_matrix *v = gretl_vector_alloc(n);
  
    if (v == NULL) {
	return v;
    }

    k = start;
    for (i=0; i<n; i++) {
	v->val[i] = k;
	if (reverse) {
	    k--;
	} else {
	    k++;
	}
    }

    return v;
}

/**
 * gretl_matrix_fill:
 * @m: matrix to fill.
 * @x: value with which to fill.
 *
 * Sets all entries in @m to the value @x.
 */

void gretl_matrix_fill (gretl_matrix *m, double x)
{
    if (m != NULL) {
	int i, n = m->rows * m->cols;

	for (i=0; i<n; i++) {
	    m->val[i] = x;
	}
    }
}

static gretl_matrix *
gretl_matrix_copy_mod (const gretl_matrix *m, int mod)
{
    gretl_matrix *c;
    int rows, cols;
    int i, j;

    if (m == NULL) {
	return NULL;
    }

    if (gretl_is_null_matrix(m)) {
	return gretl_null_matrix_new();
    }

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
	double x;
	int k = 0;

	for (j=0; j<m->cols; j++) {
	    for (i=0; i<m->rows; i++) {
		x = m->val[k++];
		gretl_matrix_set(c, j, i, x);
	    }
	}
    } else { 
	/* not transposing */
	int n = rows * cols;

	memcpy(c->val, m->val, n * sizeof *m->val);
	c->t1 = m->t1;
	c->t2 = m->t2;
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
 * gretl_matrix_copy_row:
 * @dest: destination matrix.
 * @di: row to copy into.
 * @src: source matrix.
 * @si: row to copy from.
 *
 * Copies the values from row @si of @src into row @di
 * of @dest, provided @src and @dest have the same number
 * of columns.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_matrix_copy_row (gretl_matrix *dest, int di,
			   const gretl_matrix *src, int si)
{
    int err = 0;

    if (dest == NULL || src == NULL ||
	gretl_is_null_matrix(dest) ||
	gretl_is_null_matrix(src)) {
	err = E_DATA;
    } else if (dest->cols != src->cols) {
	err = E_NONCONF;
    } else {
	double x;
	int j;
	
	for (j=0; j<src->cols; j++) {
	    x = gretl_matrix_get(src, si, j);
	    gretl_matrix_set(dest, di, j, x);
	}
    }

    return err;
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
 */

void gretl_matrix_zero (gretl_matrix *m)
{
    int i, n = m->rows * m->cols;

    for (i=0; i<n; i++) {
	m->val[i] = 0.0;
    }
}

/**
 * gretl_matrix_get_diagonal:
 * @m: input matrix.
 * @err: location to receive error code.
 *
 * Returns: a column vector containing the diagonal elements of
 * @m, otherwise %NULL.  A non-zero value is assigned via @err 
 * on failure.
 */

gretl_matrix *gretl_matrix_get_diagonal (const gretl_matrix *m, int *err)
{
    gretl_matrix *d = NULL;
    int i, n = 0;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
	d = gretl_null_matrix_new();
    } else {
	n = (m->rows < m->cols)? m->rows : m->cols;
	d = gretl_column_vector_alloc(n);
    }

    if (d == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<n; i++) {
	    d->val[i] = gretl_matrix_get(m, i, i);
	}
    } 

    return d;
}

/**
 * gretl_matrix_trace:
 * @m: square input matrix.
 * @err: location to receive error code.
 *
 * Returns: the trace (sum of diagonal elements) of @m, if 
 * @m is square, otherwise #NADBL.
 */

double gretl_matrix_trace (const gretl_matrix *m, int *err)
{
    double tr = 0.0;
    int i;

    *err = 0;
    
    if (gretl_is_null_matrix(m)) {
	return NADBL;
    }
    
    if (m->rows != m->cols) {
	*err = E_NONCONF;
	return NADBL;
    }

    for (i=0; i<m->rows; i++) {
	tr += gretl_matrix_get(m, i, i);
    }

    return tr;
}

/**
 * gretl_matrix_random_fill:
 * @m: input matrix.
 * @dist: either %D_UNIFORM or %D_NORMAL.
 *
 * Fills @m with pseudo-random values from either the uniform
 * or the standard normal distribution.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_matrix_random_fill (gretl_matrix *m, int dist)
{
    int n;

    if (m == NULL || (dist != D_UNIFORM && dist != D_NORMAL)) {
	return 1;
    }

    n = m->rows * m->cols;

    if (n == 0) {
	return 0;
    }

    if (dist == D_NORMAL) {
	gretl_rand_normal(m->val, 0, n - 1);
    } else if (dist == D_UNIFORM) {
	gretl_rand_uniform(m->val, 0, n - 1);
    }

    return 0;
}

/**
 * gretl_random_matrix_new:
 * @r: number of rows.
 * @c: number of columns.
 * @dist: either %D_UNIFORM or %D_NORMAL.
 *
 * Creates a new $r x @c matrix and filles it with pseudo-random 
 * values from either the uniform or the standard normal 
 * distribution.
 *
 * Returns: allocated matrix or %NULL on failure.
 */

gretl_matrix *gretl_random_matrix_new (int r, int c, int dist)
{
    gretl_matrix *m;

    if (dist != D_UNIFORM && dist != D_NORMAL) {
	return NULL;
    }

    if (r == 0 || c == 0) {
	return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(r, c);
    if (m == NULL) {
	return NULL;
    }

    if (dist == D_NORMAL) {
	gretl_rand_normal(m->val, 0, r * c - 1);
    } else if (dist == D_UNIFORM) {
	gretl_rand_uniform(m->val, 0, r * c - 1);
    }

    return m;
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
    int i, n, den = 0;

    if (gretl_is_null_matrix(v)) {
	return NADBL;
    }

    n = gretl_vector_get_length(v);
    if (n == 0) {
	return NADBL;
    }

    for (i=0; i<n; i++) {
	if (!na(v->val[i])) {
	    ret += v->val[i];
	    den++;
	} 
    }

    return (den > 0)? (ret / den) : NADBL;
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
    double x, xbar = 0.0;
    int i, n, den = 0;

    if (gretl_is_null_matrix(v)) {
	return NADBL;
    }

    n = gretl_vector_get_length(v);
    if (n == 0) {
	return NADBL;
    }    

    for (i=0; i<n; i++) {
	if (!na(v->val[i])) {
	    xbar += v->val[i];
	    den++;
	} 
    }

    if (den == 0) {
	return NADBL;
    }

    xbar /= den;

    for (i=0; i<n; i++) {
	x = v->val[i];
	if (!na(x)) {
	    x -= xbar;
	    s2 += x * x;
	} 
    }

    return s2 / den;
}

/**
 * gretl_matrix_resample:
 * @m: input matrix.
 * @r: number of rows in output.
 * @err: location to receive error code.
 *
 * Returns: a new matrix consisting of a random re-sampling
 * (with replacement) of the rows of @m, or %NULL on
 * failure.
 */

gretl_matrix *gretl_matrix_resample (const gretl_matrix *m, int r,
				     int *err)
{
    gretl_matrix *R = NULL;
    double x, *z = NULL;
    int i, j, k, r0;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (r <= 0) {
	r = m->rows;
    }

    R = gretl_matrix_alloc(r, m->cols);
    z = malloc(r * sizeof *z);

    if (R == NULL || z == NULL) {
	gretl_matrix_free(R);
	free(z);
	*err = E_ALLOC;
	return NULL;
    }

    r0 = m->rows;

    /* generate uniform random series */
    gretl_rand_uniform(z, 0, r - 1);

    /* sample from source matrix based on row indices */
    for (i=0; i<r; i++) {
	k = (int) (r0 * z[i]);
	if (k > r0 - 1) {
	    k = r0 - 1;
	}
	for (j=0; j<m->cols; j++) {
	    x = gretl_matrix_get(m, k, j);
	    gretl_matrix_set(R, i, j, x);
	}
    }

    free(z);

    return R;
}

static int gretl_matrix_zero_triangle (gretl_matrix *m, char t)
{
    int i, j;

    if (gretl_is_null_matrix(m)) 
	return E_DATA;

    if (m->rows != m->cols) {
	return E_NONCONF;
    }

    if (t == 'U') {
	for (i=0; i<m->rows-1; i++) {
	    for (j=i+1; j<m->cols; j++) {
		gretl_matrix_set(m, i, j, 0.0);
	    }
	}
    } else {
	for (i=1; i<m->rows; i++) {
	    for (j=0; j<i; j++) {
		gretl_matrix_set(m, i, j, 0.0);
	    }
	}
    }

    return 0;
}

/**
 * gretl_matrix_zero_upper:
 * @m: square matrix to operate on.
 *
 * Sets the elements of @m outside of the lower triangle to zero.
 * 
 * Returns: 0 on success, non-zero error code otherwise.
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
 * Returns: 0 on success, non-zero error code otherwise.
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
    int i, n = m->rows * m->cols;

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
 *
 * Returns: 0 on success, 1 if x = 0.
 */

int gretl_matrix_divide_by_scalar (gretl_matrix *m, double x)
{
    if (x == 0.0) {
	return 1;
    } else {
	gretl_matrix_multiply_by_scalar(m, 1.0 / x);
	return 0;
    }
}

/**
 * gretl_matrix_switch_sign:
 * @m: matrix to operate on.
 *
 * Changes the sign of each element of @m.
 */

void gretl_matrix_switch_sign (gretl_matrix *m)
{
    if (!gretl_is_null_matrix(m)) {
	int i, n = m->rows * m->cols;

	for (i=0; i<n; i++) {
	    m->val[i] = -m->val[i];
	}
    }
}

/**
 * gretl_matrix_raise:
 * @m: matrix to operate on.
 * @x: exponent.
 *
 * Raises each element of @m to the power @x.
 */

void gretl_matrix_raise (gretl_matrix *m, double x)
{
    if (!gretl_is_null_matrix(m)) {
	int i, n = m->rows * m->cols;

	for (i=0; i<n; i++) {
	    m->val[i] = pow(m->val[i], x);
	}
    }
}

/* if x can be represented as an integer and is an exact power of 2,
   return the exact log_2 of x, otherwise use log() to compute an
   approximation.
*/

static double log_2 (double x)
{
    int i, s;

    if (x <= 0) {
	return log(x);
    }

    if (floor(x) != x || x < 2 || x > (double) INT_MAX) {
	return log(x) / log(2.0);
    }

    s = floor(x);

    for (i=1; ; i++) {
	if (s % 2) {
	    break;
	}
	s /= 2;
	if (s == 1) {
	    return (double) i;
	}
    }

    return log(x) / log(2.0);
}

static double mexp_error_eps (int q)
{
    double x1, x2, x3, x4;
    double qf = x_factorial(q);

    x1 = pow(2.0, 3.0 - (q + q));
    x2 = qf * qf;
    x3 = x_factorial(2 * q);
    x4 = (2 * q + 1) * x3;
    x3 *= x4;

    return x1 * (x2 / x3);
}

/**
 * gretl_matrix_exp:
 * @m: square matrix to operate on.
 * @err: location to receive error code.
 *
 * Calculates the matrix exponential of @m, using algorithm
 * 11.3.1 from Golub and Van Loan, "Matrix Computations", 3e.
 *
 * Returns: the exponential, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_exp (const gretl_matrix *m, int *err)
{
    gretl_matrix *A = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *N = NULL;
    gretl_matrix *D = NULL;
    gretl_matrix *W = NULL;
    
    double xa, c, j, delta = 1.0e-13;
    int q, k, n;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    n = m->rows;

    A = gretl_matrix_copy(m);
    X = gretl_identity_matrix_new(n);
    N = gretl_identity_matrix_new(n);
    D = gretl_identity_matrix_new(n);
    W = gretl_matrix_alloc(n, n);

    if (A == NULL || X == NULL || N == NULL || 
	D == NULL || W == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    xa = gretl_matrix_infinity_norm(A);

    j = floor(log_2(xa));
    if (j < 0) {
	j = 0;
    }

    gretl_matrix_divide_by_scalar(A, pow(2.0, j));

    for (q=1; q<16; q++) {
	c = mexp_error_eps(q);
	if (c * xa <= delta) {
	    break;
	}
    }

    c = 1.0;

    for (k=1; k<=q; k++) {
	c *= (q - k + 1.0) / ((2.0*q - k + 1) * k);
	/* X = AX */
	gretl_matrix_multiply(A, X, W);
	gretl_matrix_copy_values(X, W);
	/* N = N + cX */
	gretl_matrix_multiply_by_scalar(W, c);
	gretl_matrix_add_to(N, W);
	/* D = D + (-1)^k cX */
	if (k % 2) {
	    gretl_matrix_subtract_from(D, W);
	} else {
	    gretl_matrix_add_to(D, W);
	}
    }

    /* "solve DF = N for F" */
    *err = gretl_LU_solve(D, N);

    if (!*err) {
	for (k=0; k<j; k++) {
	    gretl_matrix_multiply(N, N, W);
	    gretl_matrix_copy_values(N, W);
	}
    }

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(X);
    gretl_matrix_free(D);
    gretl_matrix_free(W);

    if (*err) {
	gretl_matrix_free(N);
	N = NULL;
    }

    return N;
}

/**
 * gretl_matrix_polroots:
 * @a: vector of coefficients.
 * @err: location to receive error code.
 *
 * Calculates the roots of the polynomial with coefficients
 * given by @a.  If the degree of the polynomial is p, then
 * @a should contain p + 1 coefficients in ascending order,
 * i.e. starting with the constant and ending with the
 * coefficient on x^p.
 *
 * Returns: a p-vector if all the roots are real, otherwise a
 * p x 2 matrix with the real parts in the first column and
 * the imaginary parts in the second.  Or %NULL on failure.
 */

gretl_matrix *gretl_matrix_polroots (const gretl_matrix *a,
				     int *err)
{
    gretl_matrix *r = NULL;
    double *xcof = NULL, *cof = NULL;
    cmplx *roots = NULL;
    int i, m, order, polerr;

    *err = 0;

    m = gretl_vector_get_length(a);

    if (m < 2) {
	*err = E_DATA;
	return NULL;
    }

    order = m - 1;

    xcof = malloc(m * sizeof *xcof);
    cof = malloc(m * sizeof *cof);
    roots = malloc(order * sizeof *roots);

    if (xcof == NULL || cof == NULL || roots == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<m; i++) {
	xcof[i] = a->val[i];
    }

    polerr = polrt(xcof, cof, order, roots);

    if (polerr) {
	*err = E_DATA;
    } else {
	int allreal = 1;

	for (i=0; i<order; i++) {
	    if (roots[i].i != 0) {
		allreal = 0;
		break;
	    }
	}

	if (allreal) {
	    r = gretl_matrix_alloc(order, 1);
	} else {
	    r = gretl_matrix_alloc(order, 2);
	}

	if (r == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}

	for (i=0; i<order; i++) {
	    gretl_matrix_set(r, i, 0, roots[i].r);
	    if (!allreal) {
		gretl_matrix_set(r, i, 1, roots[i].i);
	    }
	}
    }

 bailout:

    free(xcof);
    free(cof);
    free(roots);

    return r;
}

/**
 * gretl_matrix_copy_values:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Copies the elements of @src into the corresponding elements
 * of @targ.
 * 
 * Returns: 0 on successful completion, or
 * %E_NONCONF if the two matrices are not
 * conformable for the operation.
 */

int gretl_matrix_copy_values (gretl_matrix *targ, 
			      const gretl_matrix *src)
{
    int n;

    if (src == NULL) {
	fprintf(stderr, "gretl_matrix_copy_values: src is NULL\n");
	return E_DATA;
    }

    if (targ == src) {
	/* no-op */
	return 0;
    }

    if (targ->rows != src->rows || targ->cols != src->cols) {
	fprintf(stderr, "gretl_matrix_copy_values: targ is %d x %d but src is %d x %d\n",
		targ->rows, targ->cols, src->rows, src->cols);
	return E_NONCONF;
    }

    n = src->rows * src->cols;

    if (n > 0) {
	memcpy(targ->val, src->val, n * sizeof *targ->val);
    }

    return 0;
}

/**
 * gretl_matrix_copy_values_shaped:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Copies the elements of @src into @targ, column
 * by column.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if 
 * the two matrices do not contain the same number of
 * elements.
 */

int gretl_matrix_copy_values_shaped (gretl_matrix *targ, 
				     const gretl_matrix *src)
{
    int n = targ->rows * targ->cols;

    if (src->rows * src->cols != n) {
	fprintf(stderr, "gretl_matrix_copy_values_shaped: targ is %d x %d but src is %d x %d\n",
		targ->rows, targ->cols, src->rows, src->cols);
	return E_NONCONF;
    }

    if (n > 0) {
	memcpy(targ->val, src->val, n * sizeof *targ->val);
    }

    return 0;
}

static int add_scalar_to_matrix (gretl_matrix *targ, double x)
{
    int i, n = targ->rows * targ->cols;

    for (i=0; i<n; i++) {
	targ->val[i] += x;
    }

    return 0;
}

static int subtract_scalar_from_matrix (gretl_matrix *targ, double x)
{
    int i, n = targ->rows * targ->cols;
    
    for (i=0; i<n; i++) {
	targ->val[i] -= x;
    }

    return 0;
}

/**
 * gretl_matrix_add_to:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Adds the elements of @src to the corresponding elements
 * of @targ.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if the 
 * two matrices are not conformable for the operation.
 * In the special case where @src is in fact a scalar, the 
 * operation always goes through OK, with the scalar being 
 * added to each element of @targ.
 */

int 
gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	if (matrix_is_scalar(src)) {
	    return add_scalar_to_matrix(targ, src->val[0]);
	} else {
	    fprintf(stderr, "gretl_matrix_add_to: adding %d x %d to %d x %d\n",
		    src->rows, src->cols, targ->rows, targ->cols);
	    return E_NONCONF;
	}
    }

    n = src->rows * src->cols;
    
    for (i=0; i<n; i++) {
	targ->val[i] += src->val[i];
    }

    return 0;
}

/**
 * gretl_matrix_add_transpose_to:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Adds the elements of @src, transposed, to the corresponding 
 * elements of @targ.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if the 
 * two matrices are not conformable for the operation.
 */

int gretl_matrix_add_transpose_to (gretl_matrix *targ, 
				   const gretl_matrix *src)
{
    int i, j, k = 0;

    if (targ->rows != src->cols || targ->cols != src->rows) {
	fprintf(stderr, "gretl_matrix_add_transpose_to: "
		"adding %d x %d to %d x %d\n",
		src->cols, src->rows, targ->rows, targ->cols);
	return E_NONCONF;
    }

    /* note: the k index follows column-major order */

    for (i=0; i<src->rows; i++) {
	for (j=0; j<src->cols; j++) {
	    targ->val[k++] += gretl_matrix_get(src, i, j);
	}
    } 

    return 0;
}

/**
 * gretl_matrix_subtract_from:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Subtracts the elements of @src from the corresponding elements
 * of @targ.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if the 
 * two matrices are not conformable for the operation.
 * In the special case where @src is in fact a scalar, the 
 * operation always goes through OK, with the scalar being 
 * subtracted from each element of @targ.
 */

int 
gretl_matrix_subtract_from (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n;

    if (targ->rows != src->rows || targ->cols != src->cols) {
	if (matrix_is_scalar(src)) {
	    return subtract_scalar_from_matrix(targ, src->val[0]);
	} else {
	    return E_NONCONF;
	}
    }

    n = src->rows * src->cols;
    
    for (i=0; i<n; i++) {
	targ->val[i] -= src->val[i];
    }

    return 0;
}

/**
 * gretl_matrix_subtract_reversed:
 * @a: m x n matrix.
 * @b: m x n matrix.
 *
 * Operates on @b such that b_{ij} = a_{ij} - b_{ij}.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if the 
 * two matrices are not conformable for the operation.
 */

int 
gretl_matrix_subtract_reversed (const gretl_matrix *a, gretl_matrix *b)
{
    int i, n;

    if (a->rows != b->rows || a->cols != b->cols) {
	return E_NONCONF;
    }

    n = a->rows * b->cols;
    
    for (i=0; i<n; i++) {
	b->val[i] = a->val[i] - b->val[i];
    }

    return 0;
}

/**
 * gretl_matrix_I_minus:
 * @m: original square matrix, n x n.
 *
 * Rewrites @m as (I - m), where I is the n x n identity
 * matrix.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if @m is
 * not square.
 */

int gretl_matrix_I_minus (gretl_matrix *m)
{
    double x;
    int i, j;

    if (m->rows != m->cols) {
	return E_NONCONF;
    }

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    x = gretl_matrix_get(m, i, j);
	    if (i == j) {
		gretl_matrix_set(m, i, j, 1.0 - x);
	    } else if (x != 0.0) {
		gretl_matrix_set(m, i, j, -x);
	    }
	}
    }
    
    return 0;
}

/**
 * gretl_matrix_inscribe_I:
 * @m: original matrix.
 * @row: top row for insertion.
 * @col: leftmost row for insertion.
 * @n: dimension (rows and columns) of identity matrix.
 * 
 * Writes an n x n identity matrix into matrix @m, the top left-hand
 * corner of the insertion being given by @row and @col (which are
 * 0-based).
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if an identity
 * matrix of the specified size cannot be fitted into @m at the
 * specified location.
 */

int gretl_matrix_inscribe_I (gretl_matrix *m, int row, int col, int n)
{
    int i, j, mi, mj;

    if (n <= 0) {
	return E_NONCONF;
    }

    if (row < 0 || row + n > m->rows) {
	return E_NONCONF;
    }

    if (col < 0 || col + n > m->cols) {
	return E_NONCONF;
    }

    for (i=0; i<n; i++) {
	mi = row + i;
	for (j=0; j<n; j++) {
	    mj = col + j;
	    gretl_matrix_set(m, mi, mj, (i == j)? 1.0 : 0.0);
	}
    }
    
    return 0;
}

/**
 * gretl_matrix_transpose:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Fills out @targ (which must be pre-allocated and of the right
 * dimensions) with the transpose of @src.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_matrix_transpose (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, j, k = 0;
    double x;

    if (targ->rows != src->cols || targ->cols != src->rows) {
	return E_NONCONF;
    }

    for (j=0; j<src->cols; j++) {
	for (i=0; i<src->rows; i++) {
	    x = src->val[k++];
	    gretl_matrix_set(targ, j, i, x);
	}
    }

    return 0;
}

/**
 * gretl_square_matrix_transpose:
 * @m: square matrix to operate on.
 *
 * Transposes the matrix @m.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_square_matrix_transpose (gretl_matrix *m)
{
    double x, y;
    int mij, mji;
    int i, j;

    if (m->rows != m->cols) {
	fputs("gretl_square_matrix_transpose: matrix must be square\n", 
	      stderr);
	return 1;
    }

    for (i=0; i<m->rows-1; i++) {
	for (j=i+1; j<m->rows; j++) {
	    mij = mdx(m,i,j);
	    mji = mdx(m,j,i);
	    x = m->val[mij];
	    y = m->val[mji];
	    m->val[mij] = y;
	    m->val[mji] = x;
	}
    }

    return 0;
}

/**
 * gretl_matrix_xtr_symmetric:
 * @m: gretl_matrix.
 *
 * Computes the symmetric part of @m by averaging its off-diagonal 
 * elements.
 */

void gretl_matrix_xtr_symmetric (gretl_matrix *m)
{
    double x;
    int mij, mji;
    int i, j;

    for (i=0; i<m->rows; i++) {
	for (j=0; j<i; j++) {
	    mij = mdx(m,i,j);
	    mji = mdx(m,j,i);
	    x = m->val[mij];
	    x += m->val[mji];
	    m->val[mij] = m->val[mji] = 0.5 * x;
	}
    }
}

/**
 * gretl_matrix_add_self_transpose:
 * @m: (square) matrix to operate on.
 *
 * Adds the transpose of @m to @m itself, yielding a symmetric
 * matrix.
 * 
 * Returns: 0 on successful completion, or
 * 1 if the source matrix is not square.
 */

int gretl_matrix_add_self_transpose (gretl_matrix *m)
{
    double x;
    int mij, mji;
    int i, j;

    if (m->rows != m->cols) {
	fputs("gretl_matrix_add_self_transpose: matrix must be square\n", 
	      stderr);
	return E_NONCONF;
    }

    for (i=0; i<m->rows; i++) {
	for (j=i; j<m->rows; j++) {
	    mij = mdx(m,i,j);
	    mji = mdx(m,j,i);
	    x = m->val[mij];
	    x += m->val[mji];
	    m->val[mij] = m->val[mji] = x;
	}
    }

    return 0;
}

/**
 * gretl_matrix_vectorize:
 * @targ: target vector, (m * n) x 1.
 * @src: source matrix, m x n.
 *
 * Writes into @targ vec(@src), that is, a column vector
 * formed by stacking the columns of @src.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if
 * @targ is not correctly dimensioned.
 */

int 
gretl_matrix_vectorize (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n = src->rows * src->cols;

    if (targ->cols != 1 || targ->rows != n) {
	return E_NONCONF;
    }

    for (i=0; i<n; i++) {
	targ->val[i] = src->val[i];
    }

    return 0;
}

/**
 * gretl_matrix_unvectorize:
 * @targ: target matrix, m x n.
 * @src: source vector, (m * n) x 1.
 *
 * Writes successive blocks of length m from @src into 
 * the successive columns of @targ (that is, performs the
 * inverse of the vec() operation).
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if
 * @targ is not correctly dimensioned.
 */

int 
gretl_matrix_unvectorize (gretl_matrix *targ, const gretl_matrix *src)
{
    int i, n = targ->rows * targ->cols;

    if (src->cols != 1 || src->rows != n) {
	return E_NONCONF;
    }

    for (i=0; i<n; i++) {
	targ->val[i] = src->val[i];
    }

    return 0;
}

/**
 * gretl_matrix_vectorize_h:
 * @targ: target vector, (m * (m+1)/2) x 1.
 * @src: source square matrix, m x m.
 *
 * Writes into @targ vech(@src), that is, a column vector
 * containing the lower-triangular elements of @src.
 * This is only useful for symmetric matrices, but for the 
 * sake of performance we don't check for that.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if
 * @targ is not correctly dimensioned.
 */

int 
gretl_matrix_vectorize_h (gretl_matrix *targ, const gretl_matrix *src)
{
    int n = src->rows;
    int m = n * (n+1) / 2;
    int i, j;

    if (targ->cols != 1 || targ->rows != m) {
	return E_NONCONF;
    }

    m = 0;
    for (i=0; i<n; i++) {
	for (j=i; j<n; j++) {
	    targ->val[m++] = gretl_matrix_get(src, i, j);
	}
    }

    return 0;
}

/**
 * gretl_matrix_unvectorize_h:
 * @targ: target matrix, n x n.
 * @src: source vector, m x 1.
 *
 * Rearranges successive blocks of decreasing length from @src into 
 * the successive columns of @targ (that is, performs the
 * inverse of the vech() operation): @targ comes out symmetric.
 * 
 * Returns: 0 on successful completion, or %E_NONCONF if
 * @targ is not correctly dimensioned.
 */

int 
gretl_matrix_unvectorize_h (gretl_matrix *targ, const gretl_matrix *src)
{
    int n = targ->rows, m = src->rows;
    int i, j;
    double x;

    if (src->cols != 1 || (n * (n + 1) != 2 * m)) {
	return E_NONCONF;
    }

    m = 0;
    for (j=0; j<n; j++) {
	for (i=j; i<n; i++) {
	    x = src->val[m++];
	    gretl_matrix_set(targ, i, j, x);
	    gretl_matrix_set(targ, j, i, x);
	}
    }

    return 0;
}

/**
 * gretl_matrix_inscribe_matrix:
 * @targ: target matrix.
 * @src: source matrix.
 * @row: row offset for insertion (0-based).
 * @col: column offset for insertion.
 * @mod: either %GRETL_MOD_TRANSPOSE or %GRETL_MOD_NONE.
 *
 * Writes @src into @targ, starting at offset @row, @col.
 * The @targ matrix must be large enough to sustain the
 * inscription of @src at the specified point.  If @mod
 * is %GRETL_MOD_TRANSPOSE it is in fact the transpose of
 * @src that is written into @targ.
 * 
 * Returns: 0 on success, %E_NONCONF if the matrices are
 * not conformable for the operation.
 */

int gretl_matrix_inscribe_matrix (gretl_matrix *targ,
				  const gretl_matrix *src,
				  int row, int col,
				  GretlMatrixMod mod)
{
    int m = (mod == GRETL_MOD_TRANSPOSE)? src->cols : src->rows;
    int n = (mod == GRETL_MOD_TRANSPOSE)? src->rows : src->cols;
    double x;
    int i, j;

    if (row < 0 || col < 0) {
	return E_NONCONF;
    }

    if (row + m > targ->rows ||
	col + n > targ->cols) {
	fprintf(stderr, "gretl_matrix_inscribe_matrix: out of bounds\n");
	return E_NONCONF;
    }

    for (i=0; i<m; i++) {
	for (j=0; j<n; j++) {
	    if (mod == GRETL_MOD_TRANSPOSE) {
		x = matrix_transp_get(src, i, j);
	    } else {
		x = gretl_matrix_get(src, i, j);
	    }
	    gretl_matrix_set(targ, row + i, col + j, x);
	}
    }

    return 0;
}

/**
 * gretl_matrix_extract_matrix:
 * @targ: target matrix.
 * @src: source matrix.
 * @row: row offset for extraction (0-based).
 * @col: column offset for extraction.
 * @mod: either %GRETL_MOD_TRANSPOSE or %GRETL_MOD_NONE.
 *
 * Writes into @targ a sub-matrix of @src, taken from the
 * offset @row, @col.  The @targ matrix must be large enough
 * to provide a sub-matrix of the dimensions of @src.
 * If @mod is %GRETL_MOD_TRANSPOSE it is in fact the transpose 
 * of the sub-matrix that that is written into @targ.
 * 
 * Returns: 0 on success, %E_NONCONF if the matrices are
 * not conformable for the operation.
 */

int gretl_matrix_extract_matrix (gretl_matrix *targ,
				 const gretl_matrix *src,
				 int row, int col,
				 GretlMatrixMod mod)
{
    int m = (mod == GRETL_MOD_TRANSPOSE)? targ->cols : targ->rows;
    int n = (mod == GRETL_MOD_TRANSPOSE)? targ->rows : targ->cols;
    double x;
    int i, j, si, sj;

    if (row < 0 || col < 0) {
	return E_NONCONF;
    }

    if (row + m > src->rows ||
	col + n > src->cols) {
	return E_NONCONF;
    }

    si = row;
    for (i=0; i<m; i++) {
	sj = col;
	for (j=0; j<n; j++) {
	    x = gretl_matrix_get(src, si, sj++);
	    if (mod == GRETL_MOD_TRANSPOSE) {
		matrix_transp_set(targ, i, j, x);
	    } else {
		gretl_matrix_set(targ, i, j, x);
	    }
	}
	si++;
    }

    return 0;
}

/**
 * gretl_matrix_steal_data:
 * @m: matrix to operate on.
 *
 * "Steals" the allocated data from @m, which is left with a
 * %NULL data pointer.
 * 
 * Returns: a pointer to the "stolen" data.
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
 * gretl_matrix_print:
 * @m: matrix.
 * @msg: message to print with matrix, or %NULL.
 *
 * Prints the given matrix to stderr.  
 */

void gretl_matrix_print (const gretl_matrix *m, const char *msg) 
{
    int i, j;

    if (m == NULL || m->val == NULL) {
	if (msg != NULL && *msg != '\0') {
	    fprintf(stderr, "%s: matrix is NULL\n", msg);
	} else {
	    fputs("matrix is NULL\n", stderr);
	}
	return;
    }

    if (msg != NULL && *msg != '\0') {
	fprintf(stderr, "%s (%d x %d)", msg, m->rows, m->cols);
	if (m->t1 != 0 || m->t2 != 0) {
	    fprintf(stderr, " [t1 = %d, t2 = %d]\n\n", m->t1 + 1, m->t2 + 1);
	} else {
	    fputs("\n\n", stderr);
	}
    }

    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    fprintf(stderr, "%#12.5g ", gretl_matrix_get(m, i, j));
	}
	fputc('\n', stderr);
    }

    fputc('\n', stderr);
}

#define DEFAULT_EQTOL 1.5e-12

static double eq_tol = DEFAULT_EQTOL;

/**
 * gretl_matrix_set_equals_tolerance:
 * @tol: tolerance value.
 *
 * Sets the tolerance for judging whether or not a matrix is symmetric
 * (see gretl_matrix_is_symmetric() and also 
 * gretl_invert_symmetric_matrix()).  The tolerance is the maximum
 * relative difference between corresponding off-diagonal elements that
 * is acceptable in a supposedly "symmetric" matrix.  The default
 * value is 1.5e-12.
 */

void gretl_matrix_set_equals_tolerance (double tol)
{
    eq_tol = tol;
}

/**
 * gretl_matrix_unset_equals_tolerance:
 *
 * Sets the tolerance for judging whether or not a matrix is symmetric
 * to its default value.  See also gretl_matrix_set_equals_tolerance().
 */

void gretl_matrix_unset_equals_tolerance (void)
{
    eq_tol = DEFAULT_EQTOL;
}

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

    if (reldiff > eq_tol) {
	fprintf(stderr, "relative difference = %g\n", reldiff);
    }

    return reldiff > eq_tol;
}

static int real_gretl_matrix_is_symmetric (const gretl_matrix *m,
					   int verbose)
{
    double x, y;
    int i, j;

    if (gretl_is_null_matrix(m)) {
	return 0;
    }

    for (i=1; i<m->rows; i++) {
	for (j=0; j<i; j++) {
	    x = gretl_matrix_get(m, i, j);
	    y = gretl_matrix_get(m, j, i);
	    if (sneq(x, y)) {
		if (verbose) {
		    fprintf(stderr, "M(%d,%d) = %.16g but M(%d,%d) = %.16g\n",
			    i, j, x, j, i, y);
		    if (m->rows < 100) {
			gretl_matrix_print(m, "gretl_matrix_is_symmetric()");
		    }
		}
		return 0;
	    }
	}
    }

    return 1;
}

/**
 * gretl_matrix_is_symmetric:
 * @m: gretl_matrix.
 *
 * Returns: 1 if @m is symmetric (with a small relative tolerance
 * for asymmetry), otherwise 0.
 */

int gretl_matrix_is_symmetric (const gretl_matrix *m)
{
    return real_gretl_matrix_is_symmetric(m, 0);
}

/**
 * gretl_matrix_is_idempotent:
 * @m: gretl_matrix.
 *
 * Returns: 1 if @m is idempotent, otherwise 0.
 */

int gretl_matrix_is_idempotent (const gretl_matrix *m)
{
    gretl_matrix *b;
    int k, ret, err;

    if (gretl_is_null_matrix(m)) {
	return 0;
    }

    k = m->rows;

    if (m->cols != k) {
	return 0;
    }

    b = gretl_matrix_alloc(k, k);
    if (b == NULL) {
	return 0;
    }

    gretl_matrix_multiply(m, m, b);
    ret = gretl_matrices_are_equal(m, b, &err);
    gretl_matrix_free(b);

    return ret;
}

/**
 * gretl_matrix_infinity_norm:
 * @m: gretl_matrix.
 *
 * Returns: the infinity-norm of @m (the maximum value across
 * the rows of @m of the sum of the absolute values of
 * the elements in the given row).
 */

double gretl_matrix_infinity_norm (const gretl_matrix *m)
{
    double rsum, rmax = 0.0;
    int i, j;

    if (gretl_is_null_matrix(m)) {
	return NADBL;
    }

    for (i=0; i<m->rows; i++) {
	rsum = 0.0;
	for (j=0; j<m->cols; j++) {
	    rsum += fabs(gretl_matrix_get(m, i, j));
	}
	if (rsum > rmax) {
	    rmax = rsum;
	}
    }

    return rmax;
}

/**
 * gretl_matrix_one_norm:
 * @m: gretl_matrix.
 *
 * Returns: the 1-norm of @m (the maximum value across
 * the columns of @m of the sum of the absolute values of
 * the elements in the given column).
 */

double gretl_matrix_one_norm (const gretl_matrix *m)
{
    double csum, cmax = 0.0;
    int i, j;

    if (gretl_is_null_matrix(m)) {
	return NADBL;
    }

    for (j=0; j<m->cols; j++) {
	csum = 0.0;
	for (i=0; i<m->rows; i++) {
	    csum += fabs(gretl_matrix_get(m, i, j));
	}
	if (csum > cmax) {
	    cmax = csum;
	}
    }

    return cmax;
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
    char uplo = 'L';
    integer n, info;
    double det = NADBL;
    int i;

    if (gretl_is_null_matrix(m)) {
	return NADBL;
    }    

    n = m->rows;

    if (m->rows != m->cols) {
	fputs("gretl_vcv_log_determinant: matrix must be square\n", stderr);
	return det;
    }

    if (!real_gretl_matrix_is_symmetric(m, 1)) {
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
	double x;

	det = 1.0;
	for (i=0; i<n; i++) {
	    x = gretl_matrix_get(a, i, i);
	    det *= x * x;
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

static double gretl_LU_determinant (gretl_matrix *a, int logdet, int absval,
				    int *err)
{
    integer n, info;
    integer *ipiv;
    double det;
    int i;

    if (gretl_is_null_matrix(a)) {
	*err = E_DATA;
	return NADBL;
    }  

    *err = 0;
    n = a->rows;

    if (a->rows != a->cols) {
	fputs("gretl_LU_determinant: matrix must be square\n", stderr);
	*err = E_NONCONF;
	return NADBL;
    }

    if (a->rows == 1) {
	if (!logdet) {
	    return a->val[0];
	} else if (a->val[0] > 0) {
	    return log(a->val[0]);
	} else {
#if 0
	    fputs("gretl_matrix_log_determinant: determinant is <= 0\n", stderr);
#endif
	    *err = 1;
	    return NADBL;
	}
    }

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) {
	*err = E_ALLOC;
	return NADBL;
    }

    dgetrf_(&n, &n, a->val, &n, ipiv, &info);

    if (info > 0) {
	return logdet ? NADBL : 0;
    } else if (info < 0) {
	fprintf(stderr, "gretl_LU_determinant: dgetrf gave info = %d\n", 
		(int) info);
	free(ipiv);
	*err = 1;
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
	    double aii = gretl_matrix_get(a, i, i);

	    if (aii == 0.0) {
		fputs("gretl_matrix_log_determinant: determinant = 0\n", stderr);
		*err = 1;
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
	    *err = 1;
	    det = NADBL;
	}
    } else {
	det = 1.0;
	for (i=0; i<n; i++) {
	    if (ipiv[i] != i + 1) {
		det = -det;
	    }
	    det *= gretl_matrix_get(a, i, i);
	}
    }

    free(ipiv);

    return det;
}

/**
 * gretl_matrix_determinant:
 * @a: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the determinant of the square matrix @a using the LU
 * factorization.  Matrix @a is not preserved: it is overwritten
 * by the factorization.
 * 
 * Returns: the determinant, or #NABDL on failure.
 */

double gretl_matrix_determinant (gretl_matrix *a, int *err)
{
    return gretl_LU_determinant(a, 0, 0, err);
}

/**
 * gretl_matrix_log_determinant:
 * @a: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the log of the determinant of the square matrix @a using LU
 * factorization.  Matrix @a is not preserved: it is overwritten
 * by the factorization.
 * 
 * Returns: the determinant, or #NABDL on failure.
 */

double gretl_matrix_log_determinant (gretl_matrix *a, int *err)
{
    return gretl_LU_determinant(a, 1, 0, err);
}

/**
 * gretl_matrix_log_abs_determinant:
 * @a: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the log of the absolute value of the determinant of the 
 * square matrix @a using LU factorization.  Matrix @a is not 
 * preserved: it is overwritten by the factorization.
 * 
 * Returns: the determinant, or #NABDL on failure.
 */

double gretl_matrix_log_abs_determinant (gretl_matrix *a, int *err)
{
    return gretl_LU_determinant(a, 1, 1, err);
}

static void matrix_grab_content (gretl_matrix *targ, gretl_matrix *src)
{
    free(targ->val);
    targ->val = src->val;
    src->val = NULL;

    targ->rows = src->rows;
    targ->cols = src->cols;
    targ->t1 = src->t1;
    targ->t2 = src->t2;

    src->rows = src->cols = 0;
    src->t1 = src->t2 = 0;
}

/* least squares solution using QR */

static int overdet_solve (gretl_matrix *a, gretl_matrix *b,
			  integer m, integer n, integer nrhs,
			  integer ldb)
{
    char trans = 'N';
    integer info;
    integer lda = m;
    integer lwork = -1;
    double *work;
    int err = 0;

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	return E_ALLOC;
    }  
    
    dgels_(&trans, &m, &n, &nrhs, a->val, &lda, b->val, &ldb, 
	   work, &lwork, &info);
    if (info != 0) {
	err = wspace_fail(info, work[0]);
    }

    lwork = (integer) work[0];

    work = lapack_realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	return E_ALLOC;
    }      

    dgels_(&trans, &m, &n, &nrhs, a->val, &lda, b->val, &ldb, 
	   work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "overdet_solve: dgels gave info = %d\n", 
		(int) info);
	err = E_DATA;
    } 

    /* Note: we could retrieve the sum(s) of squared errors from
       b on output if we wanted to.  But for now we'll trim
       b to contain just the solution. */

    if (!err) {
	gretl_matrix *c;

	c = gretl_matrix_trim_rows(b, 0, m - n, &err);
	if (!err) {
	    matrix_grab_content(b, c);
	    gretl_matrix_free(c);
	}
    }

    lapack_free(work);

    return err;
}

/**
 * gretl_LU_solve:
 * @a: gretl_matrix.
 * @b: gretl_matrix.
 *
 * Solves ax = b for the unknown x, using LU decomposition.
 * On exit, @b is replaced by the solution and @a is replaced 
 * by its LU decomposition.
 * 
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_LU_solve (gretl_matrix *a, gretl_matrix *b)
{
    char trans = 'N';
    integer info;
    integer m, n;
    integer ldb, nrhs = 1;
    integer *ipiv;
    int err = 0;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	return E_DATA;
    }

    m = a->rows;
    n = a->cols;

    if (b->cols == 1) {
	ldb = b->rows;
    } else if (b->rows == 1) {
	ldb = b->cols;
    } else {
	nrhs = b->cols;
	ldb = b->rows;
    }

    if (m > n) {
	return overdet_solve(a, b, m, n, nrhs, ldb);
    } else if (m < n) {
	return E_DATA;
    }

    ipiv = malloc(n * sizeof *ipiv);
    if (ipiv == NULL) {
	return E_ALLOC;
    }

    dgetrf_(&m, &n, a->val, &n, ipiv, &info);

    if (info != 0) {
	fprintf(stderr, "gretl_LU_solve: dgetrf gave info = %d\n", 
		(int) info);
	err = (info < 0)? E_DATA : E_SINGULAR;
    }

    if (!err) {
	dgetrs_(&trans, &n, &nrhs, a->val, &n, ipiv, b->val, &ldb, &info);
	if (info != 0) {
	    fprintf(stderr, "gretl_LU_solve: dgetrs gave info = %d\n", 
		    (int) info);
	    err = E_DATA;
	}
    }    

    free(ipiv);

    return err;
}

/**
 * gretl_cholesky_decomp_solve:
 * @a: symmetric positive-definite matrix.
 * @b: vector 'x'.
 *
 * Solves ax = b for the unknown vector x, using Cholesky decomposition.
 * On exit, @b is replaced by the solution and @a is replaced by its 
 * Cholesky decomposition.
 * 
 * Returns: 0 on successful completion, or non-zero code on error.
 */

int gretl_cholesky_decomp_solve (gretl_matrix *a, gretl_matrix *b)
{
    integer n, m, info = 0;
    char uplo = 'L';

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	return E_DATA;
    }    

    n = a->cols;
    m = b->cols;

    dpotrf_(&uplo, &n, a->val, &n, &info);   
    if (info != 0) {
	fprintf(stderr, "gretl_cholesky_solve:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is not positive definite\n", stderr);
	}
	return E_SINGULAR;
    } 

    dpotrs_(&uplo, &n, &m, a->val, &n, b->val, &n, &info);
    if (info != 0) {
	fprintf(stderr, "gretl_cholesky_solve:\n"
		" dpotrs failed with info = %d (n = %d)\n", (int) info, (int) n);
	return E_SINGULAR;
    }     

    return 0;
}

/**
 * gretl_cholesky_solve:
 * @a: Cholesky-decomposed symmetric positive-definite matrix.
 * @b: vector 'x'.
 *
 * Solves ax = b for the unknown vector x, using the pre-computed
 * Cholesky decomposition of @a. On exit, @b is replaced by the 
 * solution.
 * 
 * Returns: 0 on successful completion, or non-zero code on error.
 */

int gretl_cholesky_solve (const gretl_matrix *a, gretl_vector *b)
{
    integer n, info, one = 1;
    char uplo = 'L';

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	return E_DATA;
    }    

    n = a->cols;

    dpotrs_(&uplo, &n, &one, a->val, &n, b->val, &n, &info);
    if (info != 0) {
	fprintf(stderr, "gretl_cholesky_solve:\n"
		" dpotrs failed with info = %d (n = %d)\n", (int) info, (int) n);
	return E_SINGULAR;
    }     

    return 0;
} 

#define gretl_matrix_cum(m,i,j,x) (m->val[(j)*m->rows+(i)]+=x)

#define BLAS_DEBUG 0

static int blas_nmk_min = -1; 

void set_blas_nmk_min (int n)
{
    blas_nmk_min = n;
}

int get_blas_nmk_min (void)
{
    return blas_nmk_min;
}

static int USE_BLAS (int n, int m, int k)
{
    if (blas_nmk_min >= 0) {
	double nmk = (double) n * m * k;

	if (nmk > (double) blas_nmk_min) {
#if BLAS_DEBUG
	    fprintf(stderr, "nmk = %g, using blas\n", nmk);
#endif
	    return 1;
	}
    }

    return 0;
}

static int 
gretl_blas_dsyrk (const gretl_matrix *a, int atr,
		  gretl_matrix *c, GretlMatrixMod cmod)
{
    char uplo = 'U';
    char tr = (atr)? 'T' : 'N';
    integer n = c->rows;
    integer k = (atr)? a->rows : a->cols;
    integer lda = a->rows;
    double x, alpha = 1.0, beta = 0.0;
    int i, j;

    if (cmod == GRETL_MOD_CUMULATE) {
	beta = 1.0;
    } else if (cmod == GRETL_MOD_DECUMULATE) {
	alpha = -1.0;
	beta = 1.0;
    }

    dsyrk_(&uplo, &tr, &n, &k, &alpha, a->val, &lda,
	   &beta, c->val, &n);

    for (i=0; i<n; i++) {
	for (j=i+1; j<n; j++) {
	    x = gretl_matrix_get(c, i, j);
	    gretl_matrix_set(c, j, i, x);
	}
    }

    return 0;
}

#define gretl_st_result(c,i,j,x,m)			\
    do {						\
	if (m==GRETL_MOD_CUMULATE) {			\
	    c->val[(j)*c->rows+(i)]+=x;			\
	    if (i!=j) c->val[(i)*c->rows+(j)]+=x;	\
	} else if (m==GRETL_MOD_DECUMULATE) {		\
	    c->val[(j)*c->rows+(i)]-=x;			\
	    if (i!=j) c->val[(i)*c->rows+(j)]-=x;	\
	} else {					\
	    gretl_matrix_set(c,i,j,x);			\
	    gretl_matrix_set(c,j,i,x);			\
	}						\
    } while (0);      


static int 
matrix_multiply_self_transpose (const gretl_matrix *a, int atr,
				gretl_matrix *c, GretlMatrixMod cmod)
{
    register int i, j, k;
    int nc = (atr)? a->cols : a->rows;
    int nr = (atr)? a->rows : a->cols;
    int idx1, idx2;
    double x;

    if (c->rows != nc) {
	return E_NONCONF;
    }

    if (USE_BLAS(nc, nc, nr)) {
	return gretl_blas_dsyrk(a, atr, c, cmod);
    }

    if (c->rows == 1) {
	k = a->cols * a->rows;
	if (cmod != GRETL_MOD_CUMULATE) {
	    c->val[0] = 0.0;
	}
	for (i=0; i<k; i++) {
	    c->val[0] += a->val[i] * a->val[i];
	}
	return 0;
    }

    if (atr) {
	for (i=0; i<nc; i++) {
	    for (j=i; j<nc; j++) {
		idx1 = i * a->rows;
		idx2 = j * a->rows;
		x = 0.0;
		for (k=0; k<nr; k++) {
		    x += a->val[idx1++] * a->val[idx2++];
		}
		gretl_st_result(c,i,j,x,cmod);
	    }
	} 
    } else {
	for (i=0; i<nc; i++) {
	    for (j=i; j<nc; j++) {
		idx1 = i;
		idx2 = j;
		x = 0.0;
		for (k=0; k<nr; k++) {
		    x += a->val[idx1] * a->val[idx2];
		    idx1 += a->rows;
		    idx2 += a->rows;
		}
		gretl_st_result(c,i,j,x,cmod);
	    }
	} 
    }

    return 0;
}

/**
 * gretl_matrix_XTX_new:
 * @X: matrix to process.
 *
 * Returns: a newly allocated matrix containing X'X, or
 * %NULL on error.
*/

gretl_matrix *gretl_matrix_XTX_new (const gretl_matrix *X)
{
    gretl_matrix *XTX = NULL;

    if (!gretl_is_null_matrix(X)) {
	XTX = gretl_matrix_alloc(X->cols, X->cols);
    }

    if (XTX != NULL) {
	matrix_multiply_self_transpose(X, 1, XTX, GRETL_MOD_NONE);
    }

    return XTX;
}

static void gretl_blas_dgemm (const gretl_matrix *a, int atr,
			      const gretl_matrix *b, int btr,
			      gretl_matrix *c, GretlMatrixMod cmod,
			      int m, int n, int k)
{
    char TransA = (atr)? 'T' : 'N';
    char TransB = (btr)? 'T' : 'N';
    double alpha = 1.0, beta = 0.0;

    if (cmod == GRETL_MOD_CUMULATE) {
	beta = 1.0;
    } else if (cmod == GRETL_MOD_DECUMULATE) {
	alpha = -1.0;
	beta = 1.0;
    }

    dgemm_(&TransA, &TransB, &m, &n, &k, 
	   &alpha, a->val, &a->rows, b->val, &b->rows, &beta, 
	   c->val, &c->rows);
}

/* below: a native C re-write of netlib BLAS dgemm.f: note that
   for gretl's purposes we do not support values of 'beta'
   other than 0 or 1 */

static void gretl_dgemm (const gretl_matrix *a, int atr,
			 const gretl_matrix *b, int btr,
			 gretl_matrix *c, GretlMatrixMod cmod,
			 int m, int n, int k)
{
    const double *A = a->val;
    const double *B = b->val;
    double *C = c->val;
    double x, alpha = 1.0;
    int beta = 0;
    int ar = a->rows;
    int br = b->rows;
    int cr = c->rows;
    int i, j, l;

    if (cmod == GRETL_MOD_CUMULATE) {
	beta = 1;
    } else if (cmod == GRETL_MOD_DECUMULATE) {
	alpha = -1.0;
	beta = 1;
    }

    if (!btr) {
	if (!atr) {
	    /* C := alpha*A*B + beta*C */
	    for (j=0; j<n; j++) {
		if (beta == 0) {
		    for (i=0; i<m; i++) {
			C[j*cr+i] = 0.0;
		    }
		} 
		for (l=0; l<k; l++) {
		    if (B[j*br+l] != 0.0) {
			x = alpha * B[j*br+l];
			for (i=0; i<m; i++) {
			    C[j*cr+i] += x * A[l*ar+i];
			}
		    }
		}
	    }
	} else {
	    /* C := alpha*A'*B + beta*C */
	    for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
		    x = 0.0;
		    for (l=0; l<k; l++) {
			x += A[i*ar+l] * B[j*br+l];
		    }
		    if (beta == 0) {
			C[j*cr+i] = alpha * x;
		    } else {
			C[j*cr+i] += alpha * x;
		    }
		}
	    }
	} 
    } else {
	if (!atr) {
	    /* C := alpha*A*B' + beta*C */
	    for (j=0; j<n; j++) {
		if (beta == 0) {
		    for (i=0; i<m; i++) {
			C[j*cr+i] = 0.0;
		    }
		} 
		for (l=0; l<k; l++) {
		    if (B[l*br+j] != 0.0) {
			x = alpha * B[l*br+j];
			for (i=0; i<m; i++) {
			    C[j*cr+i] += x * A[l*ar+i];
			}
		    }
		}
	    }
	} else {
	    /* C := alpha*A'*B' + beta*C */
	    for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
		    x = 0.0;
		    for (l=0; l<k; l++) {
			x += A[i*ar+l] * B[l*br+j];
		    }
		    if (beta == 0) {
			C[j*cr+i] = alpha * x;
		    } else {
			C[j*cr+i] += alpha * x;
		    }
		}
	    }
	}
    }
}  

/**
 * gretl_matrix_multiply_mod:
 * @a: left-hand matrix.
 * @amod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @b: right-hand matrix.
 * @bmod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE.
 * @c: matrix to hold the product.
 * @cmod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_CUMULATE to
 * add the result to the existing value of @c.
 *
 * Multiplies @a (or a-transpose) into @b (or b transpose),
 * with the result written into @c.
 *
 * Returns: 0 on success; non-zero error code on
 * failure.
*/

int gretl_matrix_multiply_mod (const gretl_matrix *a, GretlMatrixMod amod,
			       const gretl_matrix *b, GretlMatrixMod bmod,
			       gretl_matrix *c, GretlMatrixMod cmod)
{
    const int atr = (amod == GRETL_MOD_TRANSPOSE);
    const int btr = (bmod == GRETL_MOD_TRANSPOSE);
    int lrows, lcols;
    int rrows, rcols;

    if (gretl_is_null_matrix(a) ||
	gretl_is_null_matrix(b) ||
	gretl_is_null_matrix(c)) {
	return E_DATA;
    }

    if (a == c || b == c) {
	fputs("gretl_matrix_multiply:\n product matrix must be "
	      "distinct from both input matrices\n", stderr);
	fprintf(stderr, "a = %p, b = %p, c = %p\n", 
		(void *) a, (void *) b, (void *) c);
	return 1;
    }

    if (a == b && atr != btr && c->rows == c->cols) {
	return matrix_multiply_self_transpose(a, atr, c, cmod);
    }

    lrows = (atr)? a->cols : a->rows;
    lcols = (atr)? a->rows : a->cols;
    rrows = (btr)? b->cols : b->rows;
    rcols = (btr)? b->rows : b->cols;

    if (lcols != rrows) {
	fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
	fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
		lrows, lcols, rrows, rcols, c->rows, c->cols);
	return E_NONCONF;
    }

    if (c->rows != lrows || c->cols != rcols) {
	fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
	fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
		lrows, lcols, rrows, rcols, c->rows, c->cols);
	return E_NONCONF;
    }

    if (USE_BLAS(lrows, rcols, lcols)) { 
	gretl_blas_dgemm(a, atr, b, btr, c, cmod, lrows, rcols, lcols);
    } else {
	gretl_dgemm(a, atr, b, btr, c, cmod, lrows, rcols, lcols);
    }

    return 0;
}

/**
 * gretl_matrix_I_kronecker:
 * @p: dimension of left-hand identity matrix.
 * @B: right-hand matrix, r x s.
 * @K: target matrix, (p * r) x (p * s).
 *
 * Writes the Kronecker product of the identity matrix 
 * of order @r and @B into @K.
 *
 * Returns: 0 on success, %E_NONCONF if matrix @K is
 * not correctly dimensioned for the operation.
 */

int
gretl_matrix_I_kronecker (int p, const gretl_matrix *B,
			  gretl_matrix *K)
{
    double x, aij, bkl;
    int r, s;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;

    if (gretl_is_null_matrix(B)) {
	return E_DATA;
    }

    r = B->rows;
    s = B->cols;

    if (K->rows != p * r || K->cols != p * s) {
	return E_NONCONF;
    }
    
    for (i=0; i<p; i++) {
	ioff = i * r;
	for (j=0; j<p; j++) {
	    /* block ij is an r * s matrix, I_{ij} * B */
	    aij = (i == j)? 1 : 0;
	    joff = j * s;
	    for (k=0; k<r; k++) {
		Ki = ioff + k;
		for (l=0; l<s; l++) {
		    bkl = gretl_matrix_get(B, k, l);
		    Kj = joff + l;
		    x = aij * bkl;
		    if (x == -0.0) {
			x = 0.0;
		    }
		    gretl_matrix_set(K, Ki, Kj, x);
		}
	    }
	}
    }

    return 0;
}

/**
 * gretl_matrix_I_kronecker_new:
 * @p: dimension of left-hand identity matrix.
 * @B: right-hand matrix, r x s.
 * @err: location to receive error code.
 *
 * Writes the Kronecker product of the identity matrix 
 * of order @r and @B into a newly allocated matrix.
 *
 * Returns: the new matrix, or %NULL on failure.
 */

gretl_matrix *
gretl_matrix_I_kronecker_new (int p, const gretl_matrix *B, int *err)
{
    gretl_matrix *K;

    if (gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    K = gretl_matrix_alloc(p * B->rows, p * B->cols);

    if (K == NULL) {
	*err = E_ALLOC;
    } else {
	gretl_matrix_I_kronecker(p, B, K);
    }

    return K;
}

/**
 * gretl_matrix_kronecker_I:
 * @A: left-hand matrix, p x q.
 * @r: dimension of right-hand identity matrix.
 * @K: target matrix, (p * r) x (q * r).
 *
 * Writes the Kronecker product of @A and the identity
 * matrix of order @r into @K.
 *
 * Returns: 0 on success, %E_NONCONF if matrix @K is
 * not correctly dimensioned for the operation.
 */

int
gretl_matrix_kronecker_I (const gretl_matrix *A, int r,
			  gretl_matrix *K)
{
    double x, aij, bkl;
    int p, q;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;

    if (gretl_is_null_matrix(A)) {
	return E_DATA;
    }

    p = A->rows;
    q = A->cols;

    if (K->rows != p * r || K->cols != q * r) {
	return E_NONCONF;
    }
    
    for (i=0; i<p; i++) {
	ioff = i * r;
	for (j=0; j<q; j++) {
	    /* block ij is an r * r matrix, a_{ij} * I_r */
	    aij = gretl_matrix_get(A, i, j);
	    joff = j * r;
	    for (k=0; k<r; k++) {
		Ki = ioff + k;
		for (l=0; l<r; l++) {
		    bkl = (k == l)? 1 : 0;
		    Kj = joff + l;
		    x = aij * bkl;
		    if (x == -0.0) {
			x = 0.0;
		    }
		    gretl_matrix_set(K, Ki, Kj, x);
		}
	    }
	}
    }

    return 0;
}

/**
 * gretl_matrix_kronecker_I_new:
 * @A: left-hand matrix, p x q.
 * @r: dimension of right-hand identity matrix.
 * @err: location to receive error code.
 *
 * Writes into a newl allocated matrix the Kronecker 
 * product of @A and the identity matrix of order @r.
 *
 * Returns: the new matrix, or %NULL on failure.
 */

gretl_matrix *
gretl_matrix_kronecker_I_new (const gretl_matrix *A, int r, int *err)
{
    gretl_matrix *K;

    if (gretl_is_null_matrix(A)) {
	*err = E_DATA;
	return NULL;
    }

    K = gretl_matrix_alloc(A->rows * r, A->cols * r);

    if (K == NULL) {
	*err = E_ALLOC;
    } else {
	gretl_matrix_kronecker_I(A, r, K);
    }

    return K;
}

/**
 * gretl_matrix_kronecker_product:
 * @A: left-hand matrix, p x q.
 * @B: right-hand matrix, r x s.
 * @K: target matrix, (p * r) x (q * s).
 *
 * Writes the Kronecker product of @A and @B into @K.
 *
 * Returns: 0 on success, %E_NONCONF if matrix @K is
 * not correctly dimensioned for the operation.
 */

int
gretl_matrix_kronecker_product (const gretl_matrix *A, const gretl_matrix *B,
				gretl_matrix *K)
{
    double x, aij, bkl;
    int p, q, r, s;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;

    if (gretl_is_null_matrix(A) || 
	gretl_is_null_matrix(B) || 
	gretl_is_null_matrix(K)) {
	return E_DATA;
    }

    p = A->rows;
    q = A->cols;
    r = B->rows;
    s = B->cols;

    if (K->rows != p * r || K->cols != q * s) {
	return E_NONCONF;
    }
    
    for (i=0; i<p; i++) {
	ioff = i * r;
	for (j=0; j<q; j++) {
	    /* block ij is an r * s matrix, a_{ij} * B */
	    aij = gretl_matrix_get(A, i, j);
	    joff = j * s;
	    for (k=0; k<r; k++) {
		Ki = ioff + k;
		for (l=0; l<s; l++) {
		    bkl = gretl_matrix_get(B, k, l);
		    Kj = joff + l;
		    x = aij * bkl;
		    if (x == -0.0) {
			x = 0.0;
		    }
		    gretl_matrix_set(K, Ki, Kj, x);
		}
	    }
	}
    }

    return 0;
}

/**
 * gretl_matrix_kronecker_product_new:
 * @A: left-hand matrix, p x q.
 * @B: right-hand matrix, r x s.
 * @err: location to receive error code.
 * 
 * Returns: A newly allocated (p * r) x (q * s) matrix which 
 * is the Kronecker product of matrices @A and @B, or %NULL 
 * on failure.
 */

gretl_matrix *
gretl_matrix_kronecker_product_new (const gretl_matrix *A, 
				    const gretl_matrix *B,
				    int *err)
{
    gretl_matrix *K;
    int p, q, r, s;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    p = A->rows;
    q = A->cols;
    r = B->rows;
    s = B->cols;    
    
    K = gretl_matrix_alloc(p * r, q * s);

    if (K == NULL) {
	*err = E_ALLOC;
    } else {
	gretl_matrix_kronecker_product(A, B, K);
    }

    return K;
}

/* 
   returns the sequence of bits in the binary expansion of s
*/

static char *binary_expansion (int s, int *t, int *pow2)
{
    char *bits = NULL;
    double l2 = log_2(s);
    int k = (int) floor(l2);

    if (l2 == k) {
	*pow2 = 1;
    }

    *t = k;
    bits = calloc(k + 1, 1);
    if (bits == NULL) {
	return NULL;
    }

    while (1) {
	bits[k] = 1;
	s -= pow(2.0, k);
	if (s == 0) {
	    break;
	}
	l2 = log_2(s);
	k = (int) floor(l2);
    }

    return bits;
}

/**
 * gretl_matrix_pow:
 * @A: square source matrix.
 * @s: exponent >= 0.
 * @err: location to receive error code.
 * 
 * Calculates the matrix A^k using Golub and Van Loan's Algorithm
 * 11.2.2 ("Binary Powering").
 *
 * Returns: allocated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_pow (const gretl_matrix *A, 
				int s, int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *W = NULL;
    char *bits = NULL;
    int t, pow2 = 0;

    if (gretl_is_null_matrix(A) || s < 0) {
	*err = E_DATA;
	return NULL;
    }

    if (A->rows != A->cols) {
	*err = E_NONCONF;
	return NULL;
    }

    if (s == 0) {
	B = gretl_identity_matrix_new(A->rows);
    } else if (s == 1) {
	B = gretl_matrix_copy(A);
    }

    if (s < 2) {
	if (B == NULL) {
	    *err = E_ALLOC;
	}
	return B;
    }	

    bits = binary_expansion(s, &t, &pow2);
    if (bits == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    B = gretl_matrix_copy(A);
    C = gretl_matrix_alloc(A->rows, A->cols);

    if (!pow2) {
	W = gretl_matrix_alloc(A->rows, A->cols);
    }

    if (B == NULL || C == NULL || (W == NULL && !pow2)) {
	gretl_matrix_free(C);
	C = NULL;
	*err = E_ALLOC;
    }

    if (!*err) {
	int k, q = 0;

	while (bits[q] == 0) {
	    /* B = B^2 */
	    gretl_matrix_multiply(B, B, C);
	    gretl_matrix_copy_values(B, C);
	    q++;
	}

	if (pow2) {
	    goto done;
	}

	gretl_matrix_copy_values(C, B);

	for (k=q+1; k<=t; k++) {
	    /* B = B^2 */
	    gretl_matrix_multiply(B, B, W);
	    gretl_matrix_copy_values(B, W);
	    if (bits[k]) {
		/* C = CB */
		gretl_matrix_multiply(C, B, W);
		gretl_matrix_copy_values(C, W);
	    }
	}
    }

 done:

    gretl_matrix_free(B);
    gretl_matrix_free(W);
    free(bits);

    return C;
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
    int i, dima, dimb;
    double dp = 0.0;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	return NADBL;
    }

    dima = (a->rows > 1)? a->rows : a->cols;
    dimb = (b->rows > 1)? b->rows : b->cols;

    if (!gretl_is_vector(a) || !gretl_is_vector(b) || dima != dimb) {
	if (errp != NULL) {
	    *errp = E_NONCONF;
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

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	return NADBL;
    }    

    if (gretl_is_vector(a) && gretl_is_vector(b)) {
	return gretl_vector_dot_product(a, b, errp);
    }

    c = gretl_matrix_alloc(1, 1);
    if (c == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(a, amod, b, bmod, 
					c, GRETL_MOD_NONE);
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

enum {
    CONF_NONE = 0,
    CONF_ELEMENTS,
    CONF_A_COLVEC,
    CONF_B_COLVEC,
    CONF_A_ROWVEC,
    CONF_B_ROWVEC,
    CONF_A_SCALAR,
    CONF_B_SCALAR,
    CONF_AC_BR,
    CONF_AR_BC
};

/**
 * dot_op_conf:
 * @ra: rows of A
 * @ca: columns of A
 * @rb: rows of B
 * @cb: columns of B
 * @r: pointer to rows of result
 * @c: pointer to columns of result
 * 
 * Used to establish the dimensions of the result of a "dot" operation.
 *
 * Returns: a numeric code identifying the convention to be used;
 * %CONF_NONE indicates non-conformability.
 */

static int dot_op_conf (int ra, int ca, int rb, int cb, int *r, int *c)
{
    int confr = (ra == rb);
    int confc = (ca == cb);
    int rowva = (ra == 1);
    int colva = (ca == 1);
    int rowvb = (rb == 1);
    int colvb = (cb == 1);
    int ret = CONF_NONE;

    if (confr && confc) {
	/* element-by-element operation */
	ret = CONF_ELEMENTS;
	*r = ra;
	*c = ca;
    } else if (confr && colva) {
	/* rows match; A is a column vector */
	ret = CONF_A_COLVEC;
	*r = ra;
	*c = (colva)? cb : ca;
    } else if (confr && colvb) {
	/* rows match; B is a column vector */
	ret = CONF_B_COLVEC;
	*r = ra;
	*c = (colva)? cb : ca;
    } else if (confc && rowva) {
	/* columns match; A is a row vector */
	ret = CONF_A_ROWVEC;
	*r = (rowva)? rb : ra;
	*c = ca;
    } else if (confc && rowvb) {
	/* columns match; B is a row vector */
	ret = CONF_B_ROWVEC;
	*r = (rowva)? rb : ra;
	*c = ca;
    } else if (rowva && colva) {
	/* A is a scalar in disguise */
	ret = CONF_A_SCALAR;
	*r = rb;
	*c = cb;
    } else if (rowvb && colvb) {
	/* B is a scalar in disguise */
	ret = CONF_B_SCALAR;
	*r = ra;
	*c = ca;
    } else if (colva && rowvb) {
	/* A is a column and B is a row */
	ret = CONF_AC_BR;
	*r = ra;
	*c = cb;
    } else if (rowva && colvb) {
	/* A is a row and B is a column */
	ret = CONF_AR_BC;
	*r = rb;
	*c = ca;
    }

    return ret;
}

static double x_op_y (double x, double y, int op)
{
    switch (op) {
    case '*': 
	return x * y;
    case '/':
	return x / y;
    case '+':
	return x + y;
    case '-':
	return x - y;
    case '^':
	return pow(x, y);
    case '=':
	return x == y;
    case '>':
	return x > y;
    case '<':
	return x < y;
    default:
	return 0;
    }
}

/**
 * gretl_matrix_dot_op:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @op: operator.
 * @err: location to receive error code.
 * 
 * Returns: a new matrix, each of whose elements is the result
 * of (x op y), where x and y are the  corresponding elements of 
 * the matrices @a and @b (or %NULL on failure).
 */

gretl_matrix *gretl_matrix_dot_op (const gretl_matrix *a, 
				   const gretl_matrix *b,
				   int op, int *err)
{
    gretl_matrix *c = NULL;
    double x, y;
    int m, n, p, q, nr, nc;
    int i, j, nv;
    int conftype;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	*err = E_DATA;
	return NULL;
    }

    m = a->rows;
    n = a->cols;
    p = b->rows;
    q = b->cols;    

    conftype = dot_op_conf(m, n, p, q, &nr, &nc);

    if (conftype == CONF_NONE) {
	fputs("gretl_matrix_dot_op: matrices not conformable\n", stderr);
	*err = E_NONCONF;
	return NULL;
    }	

    c = gretl_matrix_alloc(nr, nc);
    if (c == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    errno = 0;

    switch (conftype) {
    case CONF_ELEMENTS:
	nv = m * n;
	for (i=0; i<nv; i++) {
	    c->val[i] = x_op_y(a->val[i], b->val[i], op);
	}
	break;
    case CONF_A_COLVEC:
	for (i=0; i<nr; i++) {
	    x = gretl_vector_get(a, i);
	    for (j=0; j<nc; j++) {
		y = gretl_matrix_get(b, i, j);
		y = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, y);
	    }
	}
	break;
    case CONF_B_COLVEC:
	for (i=0; i<nr; i++) {
	    y = gretl_vector_get(b, i);
	    for (j=0; j<nc; j++) {
		x = gretl_matrix_get(a, i, j);
		x = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, x);
	    }
	}
	break;
    case CONF_A_ROWVEC:
	for (j=0; j<nc; j++) {
	    x = gretl_vector_get(a, j);
	    for (i=0; i<nr; i++) {
		y = gretl_matrix_get(b, i, j);
		y = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, y);
	    }
	}
	break;
    case CONF_B_ROWVEC:
	for (j=0; j<nc; j++) {
	    y = gretl_vector_get(b, j);
	    for (i=0; i<nr; i++) {
		x = gretl_matrix_get(a, i, j);
		x = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, x);
	    }
	}
	break;
    case CONF_A_SCALAR:
	x = a->val[0];
	for (i=0; i<nr; i++) {
	    for (j=0; j<nc; j++) {
		y = gretl_matrix_get(b, i, j);
		y = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, y);
	    }
	}
	break;
    case CONF_B_SCALAR:
	y = b->val[0];
	for (i=0; i<nr; i++) {
	    for (j=0; j<nc; j++) {
		x = gretl_matrix_get(a, i, j);
		x = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, x);
	    }
	}
	break;
    case CONF_AC_BR:
	for (i=0; i<nr; i++) {
	    x = a->val[i];
	    for (j=0; j<nc; j++) {
		y = b->val[j];
		y = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, y);
	    }
	}
	break;
    case CONF_AR_BC:
	for (j=0; j<nc; j++) {
	    x = a->val[j];
	    for (i=0; i<nr; i++) {
		y = b->val[i];
		y = x_op_y(x, y, op);
		gretl_matrix_set(c, i, j, y);
	    }
	}
	break;
    }

    if (errno) {
	gretl_matrix_free(c);
	c = NULL;
	*err = E_DATA;
	gretl_errmsg_set_from_errno("gretl_matrix_dot_op");
    }

    return c;
}

static gretl_matrix *
gretl_matrix_complex_multdiv (const gretl_matrix *a, 
			      const gretl_matrix *b,
			      int multiply, int *err)
{
    gretl_matrix *c = NULL;
    double *ar, *ai;
    double *br, *bi;
    double *cr, *ci;
    int m, n, p, q;
    int i, izero = 1;
    double r2;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	*err = E_DATA;
	return NULL;
    }

    m = a->rows;
    n = a->cols;
    p = b->rows;
    q = b->cols;    

    if (m != p) {
	*err = E_NONCONF;
	return NULL;
    }

    if ((n != 1 && n != 2) || (q != 1 && q != 2)) {
	*err = E_NONCONF;
	return NULL;
    }

    p = (n == 1 && q == 1)? 1 : 2;

    c = gretl_matrix_alloc(m, p);
    if (c == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    errno = 0;

    ar = a->val; 
    ai = (a->cols == 2)? ar + m : NULL;

    br = b->val; 
    bi = (b->cols == 2)? br + m : NULL;

    cr = c->val; 
    ci = (c->cols == 2)? cr + m : NULL;

    for (i=0; i<m; i++) {
	cr[i] = ar[i] * br[i];
	if (multiply) {
	    if (ai != NULL && bi != NULL) {
		cr[i] -= ai[i] * bi[i];
	    }
	    if (ci != NULL) {
		ci[i] = 0.0;
	    }
	    if (bi != NULL) {
		ci[i] += ar[i] * bi[i];
	    }
	    if (ai != NULL) {
		ci[i] += br[i] * ai[i];
	    }
	} else {
	    r2 = br[i]*br[i] + bi[i]*bi[i];

	    if (ai != NULL && bi != NULL) {
		cr[i] += ai[i] * bi[i];
	    }
	    if (ci != NULL) {
		ci[i] = 0.0;
	    }
	    if (bi != NULL) {
		ci[i] -= ar[i] * bi[i];
	    }
	    if (ai != NULL) {
		ci[i] += br[i] * ai[i];
	    }
	    cr[i] /= r2;
	    ci[i] /= r2;
	}
	if (ci != NULL && ci[i] != 0.0) {
	    izero = 0;
	}
    }

    if (errno) {
	gretl_matrix_free(c);
	c = NULL;
	*err = E_DATA;
	gretl_errmsg_set_from_errno("gretl_matrix_complex_multdiv");
    } else if (c->cols == 2 && izero) {
	*err = gretl_matrix_realloc(c, c->rows, 1);
	if (*err) {
	    gretl_matrix_free(c);
	    c = NULL;
	}	    
    }

    return c;
}

/**
 * gretl_matrix_complex_multiply:
 * @a: m x (1 or 2) matrix.
 * @b: m x (1 or 2) matrix.
 * @err: location to receive error code.
 *
 * Computes the complex product of @a and @b.  The first
 * column in these matrices is assumed to contain real
 * values, and the second column (if present) imaginary 
 * coefficients.
 * 
 * Returns: an m x 2 matrix with the result of the multiplication 
 * of the two vectors of complex numbers. If both @a and @b have no 
 * imaginary part, the return value will be m x 1.  Or %NULL on 
 * failure.
 */

gretl_matrix *gretl_matrix_complex_multiply (const gretl_matrix *a, 
					     const gretl_matrix *b,
					     int *err)
{
    return gretl_matrix_complex_multdiv(a, b, 1, err);
}

/**
 * gretl_matrix_complex_divide:
 * @a: m x (1 or 2) matrix.
 * @b: m x (1 or 2) matrix.
 * @err: location to receive error code.
 *
 * Computes the complex division of @a over @b.  The first
 * column in these matrices is assumed to contain real
 * values, and the second column (if present) imaginary 
 * coefficients.
 * 
 * Returns: an m x 2 matrix with the result of the division 
 * of the two vectors of complex numbers. If both @a and @b have no 
 * imaginary part, the return value will be m x 1.  Or %NULL on 
 * failure.
 */

gretl_matrix *gretl_matrix_complex_divide (const gretl_matrix *a, 
					   const gretl_matrix *b,
					   int *err)
{
    return gretl_matrix_complex_multdiv(a, b, 0, err);
}

/* return sum of elements in row i */

static double row_sum (const gretl_matrix *m, int i)
{
    double x = 0.0;
    int j;

    if (i < 0 || i >= m->rows) {
	return NADBL;
    }

    for (j=0; j<m->cols; j++) {
	x += gretl_matrix_get(m, i, j);
    }

    return x;
}

/* return sum of elements in column j */

static double col_sum (const gretl_matrix *m, int j)
{
    double x = 0.0;
    int i;

    if (j < 0 || j >= m->cols) {
	return NADBL;
    }

    for (i=0; i<m->rows; i++) {
	x += gretl_matrix_get(m, i, j);
    }

    return x;
}

/* return col vector containing row sums, or row vector containing
   column sums */

static gretl_matrix *gretl_matrix_sum (const gretl_matrix *m, int bycol,
				       int *err)
{

    gretl_matrix *s = NULL;
    int dim, i;
    double x;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (bycol) {
	dim = m->cols;
	s = gretl_matrix_alloc(1, dim);
    } else {
	dim = m->rows;
	s = gretl_matrix_alloc(dim, 1);
    }
    
    if (s == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<dim; i++) {
	    x = (bycol)? col_sum(m, i) : row_sum(m, i);
	    gretl_vector_set(s, i, x);
	}
    }

    return s;
}

/**
 * gretl_matrix_row_sum:
 * @m: source matrix.
 *
 * Returns: a column vector containing the sums of
 * the rows of @m, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_row_sum (const gretl_matrix *m, int *err)
{
    return gretl_matrix_sum(m, 0, err);
}

/**
 * gretl_matrix_column_sum:
 * @m: source matrix.
 *
 * Returns: a row vector containing the sums of
 * the columns of @m, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_column_sum (const gretl_matrix *m, int *err)
{
    return gretl_matrix_sum(m, 1, err);
}

/**
 * gretl_matrix_row_mean:
 * @m: source matrix.
 * @err: location to receive error code.
 *
 * Returns: a column vector containing the means of
 * the rows of @m, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_row_mean (const gretl_matrix *m, int *err)
{
    gretl_matrix *s = gretl_matrix_sum(m, 0, err);

    if (s != NULL) {
	int i;

	for (i=0; i<m->rows; i++) {
	    s->val[i] /= m->cols;
	}
    }

    return s;
}

/**
 * gretl_matrix_column_mean:
 * @m: source matrix.
 * @err: location to receive error code.
 *
 * Returns: a row vector containing the means of
 * the columns of @m, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_column_mean (const gretl_matrix *m, int *err)
{
    gretl_matrix *s = gretl_matrix_sum(m, 1, err);

    if (s != NULL) {
	int j;

	for (j=0; j<m->cols; j++) {
	    s->val[j] /= m->rows;
	}
    }

    return s;
}

/**
 * gretl_matrix_column_sd:
 * @m: source matrix.
 * @err: location to receive error code.
 *
 * Returns: a row vector containing the standard deviations of
 * the columns of @m (without a degrees of freedom correction), 
 * or %NULL on failure.
 */

gretl_matrix *gretl_matrix_column_sd (const gretl_matrix *m, int *err)
{
    gretl_matrix *s;
    int i, j;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    s = gretl_matrix_alloc(1, m->cols);

    if (s == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (j=0; j<m->cols; j++) {
	double dev, v = 0.0, xbar = 0.0;

	for (i=0; i<m->rows; i++) {
	    xbar += gretl_matrix_get(m, i, j);
	}

	xbar /= m->rows;

	for (i=0; i<m->rows; i++) {
	    dev = gretl_matrix_get(m, i, j) - xbar;
	    v += dev * dev;
	}

	v /= m->rows;

	s->val[j] = sqrt(v);
    }

    return s;
}

/**
 * gretl_matrix_row_i_mean:
 * @m: source matrix.
 * @row: zero-based index of row.
 *
 * Returns: the mean of the elements in row @row of @m,
 * or #NADBL if the row is out of bounds.
 */

double gretl_matrix_row_i_mean (const gretl_matrix *m, int i)
{
    double x = row_sum(m, i);

    if (!na(x)) {
	x /= (double) m->cols;
    }

    return x;
}

/**
 * gretl_matrix_column_j_mean:
 * @m: source matrix.
 * @col: zero-based index of column.
 *
 * Returns: the mean of the elements in column @col of @m,
 * or #NADBL if the column is out of bounds.
 */

double gretl_matrix_column_j_mean (const gretl_matrix *m, int j)
{
    double x = col_sum(m, j);

    if (!na(x)) {
	x /= (double) m->rows;
    }    

    return x;
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
	rowmean = gretl_matrix_row_i_mean(m, i);
	for (j=0; j<m->cols; j++) {
	    gretl_matrix_cum(m, i, j, -rowmean);
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
	colmean = gretl_matrix_column_j_mean(m, j);
	for (i=0; i<m->rows; i++) {
	    gretl_matrix_cum(m, i, j, -colmean);
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

    if (gretl_is_null_matrix(m)) {
	return NULL;
    }

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
				    v, GRETL_MOD_NONE);

    if (err) {
	gretl_matrix_free(v);
	return NULL;
    } else {
	gretl_matrix_divide_by_scalar(v, (double) m->rows);
    }

    return v;
}

/**
 * gretl_matrix_quantiles:
 * @m: matrix on which to operate.
 * @p: probability.
 * @err: location to receive error code.
 * 
 * Returns: a row vector containing the @p quantiles
 * of the columns of @m, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_quantiles (const gretl_matrix *m,
				      double p, int *err)
{
    gretl_matrix *q;
    double *a = NULL;
    int i, j;

    if (gretl_is_null_matrix(m) || p <= 0 || p >= 1) {
	*err = E_DATA;
	return NULL;
    }

    q = gretl_matrix_alloc(1, m->cols);
    if (q == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    a = malloc(m->rows * sizeof *a);
    if (a == NULL) {
	*err = E_ALLOC;
	gretl_matrix_free(q);
	return NULL;
    }

    for (j=0; j<m->cols; j++) {
	for (i=0; i<m->rows; i++) {
	    a[i] = gretl_matrix_get(m, i, j);
	}
	q->val[j] = gretl_array_quantile(a, m->rows, p);
	if (na(q->val[j])) {
	    *err = E_DATA;
	    break;
	}
    }

    if (*err) {
	gretl_matrix_free(q);
	q = NULL;
    }

    free(a);

    return q;
}

/**
 * gretl_matrix_multiply:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @c: matrix to hold the product.
 * 
 * Multiplies @a into @b, with the result written into @c.
 *
 * Returns: 0 on success; non-zero error code on
 * failure.
 */

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c)
{
    int err = 0;

    if (gretl_is_null_matrix(a) || 
	gretl_is_null_matrix(b) ||
	gretl_is_null_matrix(c)) {
	return E_DATA;
    }

    if (matrix_is_scalar(a)) {
	err = gretl_matrix_copy_values(c, b);
	if (!err) {
	    gretl_matrix_multiply_by_scalar(c, a->val[0]);
	}
    } else if (matrix_is_scalar(b)) {
	err = gretl_matrix_copy_values(c, a);
	if (!err) {
	    gretl_matrix_multiply_by_scalar(c, b->val[0]);
	}
    } else {
	err = gretl_matrix_multiply_mod(a, GRETL_MOD_NONE,
					b, GRETL_MOD_NONE,
					c, GRETL_MOD_NONE);
    }

    return err;
}

/**
 * gretl_matrix_multiply_new:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @err: location for error code.
 * 
 * Multiplies @a into @b, with the result written into a newly
 * allocated matrix.
 *
 * Returns: matrix product on success, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_multiply_new (const gretl_matrix *a, 
					 const gretl_matrix *b,
					 int *err)
{
    gretl_matrix *c;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	*err = E_DATA;
	return NULL;
    }

    if (a->cols != b->rows) {
	fprintf(stderr, "gretl_matrix_multiply_new: requested (%d x %d) * (%d x %d)\n",
		a->rows, a->cols, b->rows, b->cols);
	*err = E_NONCONF;
	return NULL;
    }

    c = gretl_matrix_alloc(a->rows, b->cols);
    if (c == NULL) {
	*err = E_ALLOC;
	return NULL;
    }
    
    *err = gretl_matrix_multiply_mod(a, GRETL_MOD_NONE,
				     b, GRETL_MOD_NONE,
				     c, GRETL_MOD_NONE);

    if (*err) {
	gretl_matrix_free(c);
	c = NULL;
    }

    return c;
}

/**
 * gretl_general_matrix_rcond:
 * @m: matrix to examine.
 * @err: location to receive error code.
 * 
 * Estimates the reciprocal condition number of the general
 * real matrix @m (in the 1-norm), using the lapack 
 * functions %dgetrf and %dgecon.
 *
 * Returns: the estimate, or #NADBL on failure to allocate memory.
 */

static double gretl_general_matrix_rcond (const gretl_matrix *A, 
					  int *err)
{
    gretl_matrix *a = NULL;
    char norm = '1';
    integer m, n, lda, info;
    integer *iwork = NULL;
    integer *ipiv = NULL;
    double *work = NULL;
    double anorm, rcond = NADBL;

    *err = 0;

    if (gretl_is_null_matrix(A)) {
	return NADBL;
    }

    m = A->rows;
    n = A->cols;
    lda = A->rows;

    a = gretl_matrix_copy(A);
    work = malloc((4 * n) * sizeof *work);
    iwork = malloc(n * sizeof *iwork);
    ipiv = malloc(min(m, n) * sizeof *ipiv);

    if (a == NULL || work == NULL || iwork == NULL || ipiv == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    anorm = gretl_matrix_one_norm(a);

    dgetrf_(&m, &n, a->val, &lda, ipiv, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_general_matrix_rcond:\n"
		" dgetrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	rcond = 0.0;
    } else {
	dgecon_(&norm, &n, a->val, &lda, &anorm, &rcond, work, iwork, &info);
	if (info != 0) {
	    *err = 1;
	    rcond = NADBL;
	}
    }

 bailout:

    free(work);
    free(iwork);
    free(ipiv);
    gretl_matrix_free(a);

    return rcond;
}

/**
 * gretl_symmetric_matrix_rcond:
 * @m: matrix to examine.
 * @err: location to receive error code.
 * 
 * Estimates the reciprocal condition number of the real symmetric
 * positive definite matrix @m (in the 1-norm), using the lapack 
 * functions %dpotrf and %dpocon.
 *
 * Returns: the estimate, or #NADBL on failure to allocate memory.
 */

double gretl_symmetric_matrix_rcond (const gretl_matrix *m, int *err)
{
    gretl_matrix *a = NULL;
    char uplo = 'L';
    integer n, lda;
    integer info, *iwork = NULL;
    double *work = NULL;
    double anorm, rcond = NADBL;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
	return NADBL;
    }

    n = m->rows;
    lda = m->rows;

    a = gretl_matrix_copy(m);
    work = malloc((3 * n) * sizeof *work);
    iwork = malloc(n * sizeof *iwork);

    if (a == NULL || work == NULL || iwork == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    anorm = gretl_matrix_one_norm(a);

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_symmetric_matrix_rcond:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	rcond = 0.0;
    } else {
	dpocon_(&uplo, &n, a->val, &lda, &anorm, &rcond, work, iwork, &info);
	if (info != 0) {
	    *err = 1;
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
 * gretl_matrix_rcond:
 * @m: matrix to examine.
 * @err: location to receive error code.
 * 
 * Estimates the reciprocal condition number of the real 
 * matrix @m (in the 1-norm).
 *
 * Returns: the estimate, or #NADBL on failure to allocate memory.
 */

double gretl_matrix_rcond (const gretl_matrix *m, int *err)
{
    if (gretl_matrix_is_symmetric(m)) {
	return gretl_symmetric_matrix_rcond(m, err);
    } else {
	return gretl_general_matrix_rcond(m, err);
    }
}

/**
 * gretl_matrix_cholesky_decomp:
 * @a: matrix to operate on.
 * 
 * Computes the Cholesky factorization of the symmetric,
 * positive definite matrix @a.  On exit the lower triangle of 
 * @a is replaced by the factor L, as in a = LL', and the
 * upper triangle is set to zero.  Uses the lapack function 
 * %dpotrf.
 *
 * Returns: 0 on success; 1 on failure.
 */

int gretl_matrix_cholesky_decomp (gretl_matrix *a)
{
    char uplo = 'L';
    integer n, lda;
    integer info;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    n = a->rows;
    lda = a->rows;

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
    } else {
	gretl_matrix_zero_upper(a);
    }

    return (info == 0)? 0 : 1;
}

#if 0 /* experimental */

int gretl_matrix_QR_pivot_decomp (gretl_matrix *M, gretl_matrix *R,
				  int **order)
{
    integer m = M->rows;
    integer n = M->cols;

    integer info = 0;
    integer lwork = -1;
    integer lda = m;
    integer *iwork = NULL;
    doublereal *tau = NULL;
    doublereal *work = NULL;
    doublereal *work2;
    integer *jpvt = NULL;

    int i, j;
    int moved = 0;
    int err = 0;

    if (R == NULL || R->rows != n || R->cols != n) {
	return E_NONCONF;
    }

    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = malloc(sizeof *work);
    iwork = malloc(n * sizeof *iwork);

    if (tau == NULL || work == NULL || iwork == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* workspace size query */
    jpvt = malloc(n * sizeof *jpvt);
    if (jpvt == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    dgeqp3_(&m, &n, M->val, &lda, jpvt, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqp3: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    }

    /* optimally sized work array */
    lwork = (integer) work[0];
    work2 = realloc(work, (size_t) lwork * sizeof *work);
    if (work2 == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    work = work2;

    /* run actual QR factorization */
    dgeqp3_(&m, &n, M->val, &lda, jpvt, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqp3: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    }

    /* copy the upper triangular R out of M */
    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    if (i <= j) {
		gretl_matrix_set(R, i, j, 
				 gretl_matrix_get(M, i, j));
	    } else {
		gretl_matrix_set(R, i, j, 0.0);
	    }
	}
    }

    /* obtain the real "Q" matrix (in M) */
    dorgqr_(&m, &n, &n, M->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dorgqr: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    } 

 bailout:

    free(tau);
    free(work);
    free(iwork);

    for (i=0; i<n; i++) {
	if (jpvt[i] != i + 1) {
	    fprintf(stderr, "column was moved\n");
	    moved = 1;
	}
    }

    if (moved && order != NULL) {
	*order = malloc(n * sizeof **order);
	if (*order == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		(*order)[i] = jpvt[i] - 1;
	    }
	}
    }

    free(jpvt);

    return err;
}

#endif

/**
 * gretl_matrix_QR_decomp:
 * @M: m x n matrix to be decomposed.
 * @R: n x n matrix into which to write R, as in M = Q * R,
 * or %NULL if this is not wanted.
 * 
 * Computes the QR factorization of @M.  On successful exit
 * the matrix @M holds Q, and, if @R is not %NULL, the upper 
 * triangle of @R holds R.  Uses the lapack functions 
 * %dgeqrf and %dorgqr.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_matrix_QR_decomp (gretl_matrix *M, gretl_matrix *R)
{
    integer m, n, lda;
    integer info = 0;
    integer lwork = -1;
    doublereal *tau = NULL;
    doublereal *work = NULL;
    doublereal *work2;
    int i, j;
    int err = 0;

    if (gretl_is_null_matrix(M)) {
	return E_DATA;
    }

    lda = m = M->rows;
    n = M->cols;

    if (R != NULL && (R->rows != n || R->cols != n)) {
	return E_NONCONF;
    }

    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = lapack_malloc(sizeof *work);

    if (tau == NULL || work == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* workspace size query */
    dgeqrf_(&m, &n, M->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    }

    /* optimally sized work array */
    lwork = (integer) work[0];
    work2 = lapack_realloc(work, (size_t) lwork * sizeof *work);
    if (work2 == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    work = work2;

    /* run actual QR factorization */
    dgeqrf_(&m, &n, M->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    }

    if (R != NULL) {
	/* copy the upper triangular R out of M */
	double x;

	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		if (i <= j) {
		    x = gretl_matrix_get(M, i, j);
		    gretl_matrix_set(R, i, j, x);
		} else {
		    gretl_matrix_set(R, i, j, 0.0);
		}
	    }
	}
    }

    /* obtain the real "Q" matrix (in M) */
    dorgqr_(&m, &n, &n, M->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dorgqr: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    } 

 bailout:

    free(tau);
    lapack_free(work);

    return err;
}

#define QR_RCOND_MIN 1e-15 /* experiment with this? */

static int get_R_rank (const gretl_matrix *R)
{
    double d;
    int i, rank = R->rows;

#if 0
    gretl_matrix_print(R, "R, in get_R_rank");
#endif

    for (i=0; i<R->rows; i++) {
	d = gretl_matrix_get(R, i, i);
	if (isnan(d) || isinf(d) || fabs(d) < R_DIAG_MIN) {
	    rank--;
	}
    }

    return rank;
}

/**
 * gretl_check_QR_rank:
 * @R: matrix R from QR decomposition.
 * @err: pointer to receive error code.
 * 
 * Checks the reciprocal condition number of R and calculates 
 * the rank of the matrix QR.
 *
 * Returns: on success, the rank of QR.
 */

int gretl_check_QR_rank (const gretl_matrix *R, int *err)
{
    integer *iwork = NULL;
    doublereal *work = NULL;
    integer n, info = 0;
    doublereal rcond;
    char uplo = 'U';
    char diag = 'N';
    char norm = '1';
    int rank;

    if (gretl_is_null_matrix(R)) {
	*err = E_DATA;
	return 0;
    }

    *err = 0;

    rank = n = R->rows;
    work = lapack_malloc(3 * n * sizeof *work);
    iwork = malloc(n * sizeof *iwork);

    if (work == NULL || iwork == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    dtrcon_(&norm, &uplo, &diag, &n, R->val, &n, &rcond, work, 
	    iwork, &info);

    if (info != 0) {
	fprintf(stderr, "dtrcon: info = %d\n", (int) info);
	*err = 1;
	goto bailout;
    }

    if (rcond < QR_RCOND_MIN) {
	fprintf(stderr, "gretl_matrix_QR_rank: rcond = %g\n", rcond);
	rank = get_R_rank(R);
    }

 bailout:

    lapack_free(work);
    free(iwork);

    return rank;
}

/**
 * gretl_matrix_rank:
 * @a: matrix to examine.
 * @err: location to receive error code on failure.
 * 
 * Computes the rank of @a via its SV decomposition.
 *
 * Returns: the rank of @a, or 0 on failure.
 */

int gretl_matrix_rank (const gretl_matrix *a, int *err)
{
    gretl_matrix *S = NULL;
    double smin = SVD_SMIN;
    double macheps = 2.0e-16;
    int i, k, rank = 0;

    if (gretl_is_null_matrix(a)) {
	return 0;
    }

    k = (a->rows < a->cols)? a->rows : a->cols;

    *err = gretl_matrix_SVD(a, NULL, &S, NULL);

    if (!*err) {
	smin = 1.0e4 * macheps * gretl_matrix_infinity_norm(a);
	for (i=0; i<k; i++) {
	    if (S->val[i] > smin) {
		rank++;
	    }
	}
    } 

    gretl_matrix_free(S);

    return rank;
}

/**
 * gretl_invert_triangular_matrix:
 * @a: triangular matrix to invert.
 * @uplo: 'L' for lower triangular @a, 'U' for upper.
 * 
 * Computes the inverse of a triangular matrix.  On exit
 * @a is overwritten with the inverse.  Uses the lapack 
 * function %dtrtri.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_triangular_matrix (gretl_matrix *a, char uplo)
{
    char diag = 'N';
    integer n, info = 0;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    n = a->rows;

    if (a->rows != a->cols) {
	return E_NONCONF;
    }

    dtrtri_(&uplo, &diag, &n, a->val, &n, &info);

    if (info < 0) {
	err = E_DATA;
    } else if (info > 0) {
	err = E_SINGULAR;
    }

    return err;
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
    integer m, n;
    integer info;
    integer lwork;
    integer *ipiv;
    int lipiv;
    int err = 0;
    double *work;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }
    
    m = a->rows;
    n = a->cols;

    lipiv = (m <= n)? m : n;

    ipiv = malloc(lipiv * sizeof *ipiv);
    if (ipiv == NULL) {
	return E_ALLOC;
    }

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return E_ALLOC;
    }    

    dgetrf_(&m, &n, a->val, &m, ipiv, &info);   

    if (info != 0) {
	free(ipiv);
	fprintf(stderr, "dgetrf: matrix is singular\n");
	return E_SINGULAR;
    }

    lwork = -1;
    dgetri_(&n, a->val, &n, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	free(ipiv);
	return wspace_fail(info, work[0]);
    }

    lwork = (integer) work[0];

#ifdef LAPACK_DEBUG
    printf("dgetri: workspace = %d\n", (int) lwork);
#endif

    work = lapack_realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	free(ipiv);
	return E_ALLOC;
    }  

    dgetri_(&n, a->val, &n, ipiv, work, &lwork, &info);

#ifdef LAPACK_DEBUG
    printf("dgetri: info = %d\n", (int) info);
#endif

    lapack_free(work);
    free(ipiv);

    if (info != 0) {
	fprintf(stderr, "dgetri: matrix is singular\n");
	err = E_SINGULAR;
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
		x = gretl_matrix_get(m, i, j);
		gretl_matrix_set(m, j, i, x);
	    } else {
		x = gretl_matrix_get(m, j, i);
		gretl_matrix_set(m, i, j, x);
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

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    if (a->cols != a->rows) {
	fputs("gretl_invert_diagonal_matrix: input is not square\n",
	      stderr);
	return E_NONCONF;
    }

    for (i=0; i<a->rows; i++) {
	if (gretl_matrix_get(a, i, i) == 0.0) {
	    return E_SINGULAR;
	}
    }

    for (i=0; i<a->rows; i++) {
	x = gretl_matrix_get(a, i, i);
	gretl_matrix_set(a, i, i, 1.0 / x);
    }

    return 0;
}

/**
 * gretl_invert_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of matrix @a: on exit @a is 
 * overwritten with the inverse.  If @a is diagonal
 * or symmetric, appropriate simple inversion routines 
 * are called.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_matrix (gretl_matrix *a)
{
    int s, err = 0;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    s = gretl_matrix_get_structure(a);
    
    if (s == GRETL_MATRIX_IDENTITY) {
	return 0;
    } else if (s == GRETL_MATRIX_DIAGONAL) {
	err = gretl_invert_diagonal_matrix(a);
    } else if (s == GRETL_MATRIX_SYMMETRIC) {
	err = gretl_invert_symmetric_matrix(a);
	if (err) {
	    err = gretl_invert_symmetric_indef_matrix(a);
	}
    } else if (s == GRETL_MATRIX_LOWER_TRIANGULAR) {
	err = gretl_invert_triangular_matrix(a, 'L');
    } else if (s == GRETL_MATRIX_UPPER_TRIANGULAR) {
	err = gretl_invert_triangular_matrix(a, 'U');
    } else if (s >= GRETL_MATRIX_SQUARE) {
	err = gretl_invert_general_matrix(a);
    } else {
	err = E_NONCONF;
    }

    return err;
}

#define RS_RCOND_MIN 1.0e-15

/**
 * gretl_invert_symmetric_indef_matrix:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a real symmetric matrix via the 
 * Bunch-Kaufman diagonal pivoting method.  Uses the lapack 
 * functions %dsytrf and %dsytri.  On exit @a is overwritten
 * with the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_symmetric_indef_matrix (gretl_matrix *a)
{
    char uplo = 'U';
    integer n, info;
    integer *ipiv;
    integer *iwork;
    integer lwork = -1;
    double anorm, rcond;
    double *work;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    if (a->cols != a->rows) {
	fputs("gretl_invert_symmetric_indef_matrix: input is not square\n",
	      stderr);
	return E_NONCONF;
    }

    n = a->rows;
    ipiv = malloc(n * sizeof *ipiv);
    iwork = malloc(n * sizeof *iwork);
    work = lapack_malloc(sizeof *work);

    if (ipiv == NULL || iwork == NULL || work == NULL) {
	err = E_ALLOC;
	goto bailout;
    }  

    anorm = gretl_matrix_one_norm(a);

    /* workspace query */
    dsytrf_(&uplo, &n, a->val, &n, ipiv, work, &lwork, &info);   
    if (info != 0 || work[0] <= 0.0) {
	err = wspace_fail(info, work[0]);
	goto bailout;
    }

    lwork = (integer) work[0];
#ifdef LAPACK_DEBUG
    printf("dsytrf: workspace = %d\n", (int) lwork);
#endif

    if (lwork < 2 * n) {
	lwork = 2 * n;
    }

    work = lapack_realloc(work, lwork * sizeof *work);
    if (work == NULL) {
	err = E_ALLOC;
	goto bailout;
    }  

    /* decompose */
    dsytrf_(&uplo, &n, a->val, &n, ipiv, work, &lwork, &info); 
    if (info != 0) {
	fprintf(stderr, "dsytrf: matrix is singular\n");
	err = E_SINGULAR;
	goto bailout;
    }

    /* check condition number */
    dsycon_(&uplo, &n, a->val, &n, ipiv, &anorm, &rcond,
	    work, iwork, &info);
    if (info != 0) {
	fprintf(stderr, "dsycon: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    } else if (rcond < RS_RCOND_MIN) {
	fprintf(stderr, "dsycon: rcond = %g\n", rcond);
	err = E_SINGULAR;
	goto bailout;
    }

    /* invert */
    dsytri_(&uplo, &n, a->val, &n, ipiv, work, &info);

#ifdef LAPACK_DEBUG
    printf("dsytri: info = %d\n", (int) info);
#endif

 bailout:

    lapack_free(work);
    free(ipiv);
    free(iwork);

    if (!err) {
	if (info != 0) {
	    fputs("dsytri: matrix is singular\n", stderr);
	    err = E_SINGULAR;
	} else {
	    gretl_symmetric_matrix_expand(a, uplo);
	}
    }

    return err;
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
 */

int gretl_invert_symmetric_matrix (gretl_matrix *a)
{
    integer n, info;
    double *aval = NULL;
    size_t bytes;
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    if (a->cols != a->rows) {
	fputs("gretl_invert_symmetric_matrix: input is not square\n",
	      stderr);
	return E_NONCONF;
    }

    n = a->cols;

    if (n == 1) {
	a->val[0] = 1.0 / a->val[0];
	return 0;
    }

    if (!real_gretl_matrix_is_symmetric(a, 1)) {
	fputs("gretl_invert_symmetric_matrix: matrix is not symmetric\n", stderr);
	return 1;
    }

    /* back-up, just in case */
    bytes = n * n * sizeof *aval;
    aval = lapack_malloc(bytes);
    if (aval == NULL) {
	return E_ALLOC;
    }

    memcpy(aval, a->val, bytes);

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is not positive definite\n", stderr);
	}
	err = E_SINGULAR;
    } 

    if (!err) {
	dpotri_(&uplo, &n, a->val, &n, &info);
	if (info != 0) {
	    err = E_SINGULAR;
	    fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		    " dpotri failed with info = %d\n", (int) info);
	} else {
	    gretl_symmetric_matrix_expand(a, uplo);
	}
    }

    if (err) {
	memcpy(a->val, aval, bytes);
    }
    
    lapack_free(aval);

    return err;
}

/**
 * gretl_invpd:
 * @a: matrix to invert.
 * 
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization.  On exit @a is overwritten with 
 * the inverse. Uses the lapack functions %dpotrf and %dpotri.
 * Little checking is done, for speed: we assume the caller
 * knows what he's doing.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invpd (gretl_matrix *a)
{
    integer n, info;
    char uplo = 'L';
    int err = 0;

    if (a->cols != a->rows) {
	fputs("gretl_invert_symmetric_matrix: input is not square\n",
	      stderr);
	return E_NONCONF;
    }

    n = a->cols;

    if (n == 1) {
	a->val[0] = 1.0 / a->val[0];
	return 0;
    }

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is not positive definite\n", stderr);
	}
	err = E_SINGULAR;
    } 

    if (!err) {
	dpotri_(&uplo, &n, a->val, &n, &info);
	if (info != 0) {
	    err = E_SINGULAR;
	    fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		    " dpotri failed with info = %d\n", (int) info);
	} else {
	    gretl_symmetric_matrix_expand(a, uplo);
	}
    }

    return err;
}

/**
 * gretl_inverse_from_cholesky_decomp:
 * @targ: matrix to hold inverse.
 * @src: Cholesky-decomposed matrix.
 * 
 * Computes in @targ the inverse of a symmetric positive definite 
 * matrix, on the assumption that the original matrix this has already 
 * been Cholesky-decomposed in @src.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_inverse_from_cholesky_decomp (gretl_matrix *targ, 
					const gretl_matrix *src)
{
    integer info, n;
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(targ) || gretl_is_null_matrix(src)) {
	return E_DATA;
    }

    n = src->cols;

    if (n != src->rows || targ->cols != targ->rows || targ->cols != n) {
	return E_NONCONF;
    }

    memcpy(targ->val, src->val, n * n * sizeof *src->val);

    dpotri_(&uplo, &n, targ->val, &n, &info);

    if (info != 0) {
	err = E_SINGULAR;
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotri failed with info = %d\n", (int) info);
    } else {
	gretl_symmetric_matrix_expand(targ, uplo);
    }

    return err;
}

/**
 * gretl_invert_symmetric_matrix2:
 * @a: matrix to invert.
 * @ldet: location to recieve log determinant, or %NULL.
 * 
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization, computing the log-determinant 
 * in the process.  On exit @a is overwritten with the inverse 
 * and if @ldet is not %NULL the log-determinant is written to 
 * that location.  Uses the lapack functions %dpotrf and %dpotri.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_symmetric_matrix2 (gretl_matrix *a, double *ldet)
{
    integer n, info;
    char uplo = 'L';
    int i, err = 0;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    if (a->cols != a->rows) {
	fputs("gretl_invert_symmetric_matrix: input is not square\n",
	      stderr);
	return E_NONCONF;
    }

    n = a->cols;

    if (n == 1) {
	if (ldet != NULL) {
	    *ldet = log(a->val[0]);
	}
	a->val[0] = 1.0 / a->val[0];
	return 0;
    }

    if (!real_gretl_matrix_is_symmetric(a, 1)) {
	fputs("gretl_invert_symmetric_matrix: matrix is not symmetric\n", stderr);
	return 1;
    }

    dpotrf_(&uplo, &n, a->val, &n, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is not positive definite\n", stderr);
	}
	return E_SINGULAR;
    } 

    if (ldet != NULL) {
	double x = 0.0;

	for (i=0; i<n; i++) {
	    x += log(gretl_matrix_get(a,i,i));
	}
	*ldet = 2.0 * x;
    }

    dpotri_(&uplo, &n, a->val, &n, &info);

    if (info != 0) {
	err = E_SINGULAR;
	fprintf(stderr, "gretl_invert_symmetric_matrix:\n"
		" dpotri failed with info = %d\n", (int) info);
    } else {
	gretl_symmetric_matrix_expand(a, uplo);
    }

    return err;
}

#if 0
static int invert_packed_symm_indef_matrix (gretl_matrix *v,
					    integer n)
{
    integer info;
    integer *ipiv = NULL;
    double *work = NULL;
    char uplo = 'L';
    int err = 0;

    ipiv = malloc(n * sizeof *ipiv);
    work = malloc(n * sizeof *work);

    if (ipiv == NULL || work == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    dsptrf_(&uplo, &n, v->val, ipiv, &info);   

    if (info != 0) {
	fprintf(stderr, "invert_packed_symm_indef_matrix:\n"
		" dsptrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is singular\n", stderr);
	}
	err = E_SINGULAR;
	goto bailout;
    } 

    dsptri_(&uplo, &n, v->val, ipiv, work, &info);
    if (info != 0) {
	err = E_SINGULAR;
	fprintf(stderr, "invert_packed_symm_indef_matrix:\n"
		" dsptri failed with info = %d\n", (int) info);
    } 

 bailout:

    free(ipiv);
    free(work);

    return err;
}
#endif

/**
 * gretl_invert_packed_symmetric_matrix:
 * @v: symmetric matrix in vech form (lower triangle packed
 * as a column vector).
 * 
 * Computes the inverse of a symmetric positive definite matrix,
 * stored in vech form, using Cholesky factorization.  On exit 
 * @v is overwritten with the lower triangle of the inverse.
 * Uses the lapack functions %dpptrf and %dpptri.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_packed_symmetric_matrix (gretl_matrix *v)
{
    integer info, n;
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(v)) {
	return E_DATA;
    }

    if (v->cols != 1) {
	fprintf(stderr, "gretl_invert_packed_symmetric_matrix:\n"
		" matrix is not in vech form\n");
	return E_DATA;
    }

    if (v->rows == 1) {
	v->val[0] = 1.0 / v->val[0];
	return 0;
    }

    n = (integer) ((sqrt(1.0 + 8.0 * v->rows) - 1.0) / 2.0);

    dpptrf_(&uplo, &n, v->val, &info);   

    if (info != 0) {
	fprintf(stderr, "gretl_invert_packed_symmetric_matrix:\n"
		" dpptrf failed with info = %d (n = %d)\n", (int) info, (int) n);
	if (info > 0) {
	    fputs(" matrix is not positive definite\n", stderr);
	}
	return E_SINGULAR;
    } 

    dpptri_(&uplo, &n, v->val, &info);

    if (info != 0) {
	err = E_SINGULAR;
	fprintf(stderr, "gretl_invert_packed_symmetric_matrix:\n"
		" dpptri failed with info = %d\n", (int) info);
    } 

    return err;
}

static int real_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs, int rank,
			    int flag)
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
	es[i].vr = evals->val[i];
	if (flag == GRETL_MATRIX_SYMMETRIC) {
	    es[i].vi = 0.0;
	} else {
	    es[i].vi = evals->val[i + n];
	}
	es[i].idx = i;
    }

    qsort(es, n, sizeof *es, gretl_inverse_compare_doubles);

    for (i=0; i<n; i++) {
	evals->val[i] = es[i].vr;
	if (flag != GRETL_MATRIX_SYMMETRIC) {
	    evals->val[i + n] = es[i].vi;
	}
    }

    for (j=0; j<h; j++) {
	for (i=0; i<n; i++) {
	    x = gretl_matrix_get(evecs, i, es[j].idx);
	    gretl_matrix_set(tmp, i, j, x);
	}
    }

    free(evecs->val);
    evecs->val = tmp->val;
    tmp->val = NULL;

    evecs->cols = tmp->cols;

    free(tmp);
    free(es);

    return 0;
}

/**
 * gretl_general_eigen_sort:
 * @evals: array of eigenvalues, from general (not necessarily
 * symmetric) matrix.
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

int gretl_general_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs, 
			      int rank)
{
    return real_eigen_sort(evals, evecs, rank, 0);
}

/**
 * gretl_symmetric_eigen_sort:
 * @evals: array of real eigenvalues, from symmetric matrix.
 * @evecs: matrix of eigenvectors.
 * @rank: desired number of columns in output.
 * 
 * Sorts the eigenvalues in @evals from largest to smallest, and 
 * rearranges the columns in @evecs correspondingly.  If @rank is 
 * greater than zero and less than the number of columns in @evecs, 
 * then on output @evecs is shrunk so that it contains only the 
 * columns associated with the largest @rank eigenvalues.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_symmetric_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs, 
				int rank)
{
    return real_eigen_sort(evals, evecs, rank, GRETL_MATRIX_SYMMETRIC);    
}

static void gretl_matrix_replace_data (gretl_matrix *m, double *x)
{
    free(m->val);
    m->val = x;
}

/**
 * gretl_general_matrix_eigenvals:
 * @m: square matrix on which to operate.
 * @eigenvecs: non-zero to calculate eigenvectors, 0 to omit.
 * @err: location to receive error code.
 * 
 * Computes the eigenvalues of the general matrix @m.  
 * If @eigenvecs is non-zero, also compute the right
 * eigenvectors of @m, which are stored in @m. Uses the lapack 
 * function %dgeev.
 * 
 * Returns: allocated matrix containing the eigenvalues, or %NULL
 * on failure.  The returned matrix, on successful completion,
 * is n x 2 (where n = the number of rows and columns in the
 * matrix @m); the first column contains the real parts of 
 * the eigenvalues of @m, and the second holds the 
 * imaginary parts.
 */

gretl_matrix *
gretl_general_matrix_eigenvals (gretl_matrix *m, int eigenvecs, int *err) 
{
    gretl_matrix *evals = NULL;
    integer n, info, lwork;
    integer nvr, nvl = 2;
    char jvr, jvl = 'N';
    double *work, *work2;
    double *wr = NULL, *wi = NULL, *vr = NULL;
    double nullvl[2] = {0.0};
    double nullvr[2] = {0.0};

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (m->rows != m->cols) {
	fprintf(stderr, "gretl_general_matrix_eigenvals:\n"
		" matrix must be square, is %d x %d\n", m->rows, m->cols);
	*err = E_NONCONF;
	return NULL;
    }

    n = m->rows;

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    evals = gretl_zero_matrix_new(n, 2);
    if (evals == NULL) {
	*err = E_ALLOC;
	goto bailout;
    } else {
	wr = evals->val;
	wi = wr + n;
    }

    if (eigenvecs) {
	/* eigenvectors wanted */
	vr = malloc(n * n * sizeof *vr);
	if (vr == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
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
	*err = wspace_fail(info, work[0]);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = lapack_realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	*err = E_ALLOC;
	goto bailout;
    } else {
	work = work2;
    }

    dgeev_(&jvl, &jvr, &n, m->val, &n, wr, wi, nullvl, 
	   &nvl, vr, &nvr, work, &lwork, &info);

    if (info != 0) {
	*err = 1;
    } 

 bailout:

    lapack_free(work);

    if (*err) {
	gretl_matrix_free(evals);
	evals = NULL;
	if (vr != NULL) {
	    free(vr);
	}
    } else if (eigenvecs) {
	gretl_matrix_replace_data(m, vr);
    }	

    return evals;
}

/**
 * gretl_symmetric_matrix_eigenvals:
 * @m: n x n matrix to operate on.
 * @eigenvecs: non-zero to calculate eigenvectors, 0 to omit.
 * @err: location to receive error code.
 * 
 * Computes the eigenvalues of the real symmetric matrix @m.  
 * If @eigenvecs is non-zero, also compute the orthonormal
 * eigenvectors of @m, which are stored in @m. Uses the lapack 
 * function %dsyev.
 *
 * Returns: n x 1 matrix containing the eigenvalues, or %NULL
 * on failure.
 */

gretl_matrix *
gretl_symmetric_matrix_eigenvals (gretl_matrix *m, int eigenvecs, int *err) 
{
    integer n, info, lwork;
    gretl_matrix *evals = NULL;
    double *work = NULL;
    double *work2 = NULL;
    double *w = NULL;
    char uplo = 'U', jobz = (eigenvecs)? 'V' : 'N';

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (!real_gretl_matrix_is_symmetric(m, 1)) {
	fputs("gretl_symmetric_matrix_eigenvals: matrix is not symmetric\n", stderr);
	*err = E_NONCONF;
	return NULL;
    }

    n = m->rows;

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    evals = gretl_column_vector_alloc(n);
    if (evals == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    w = evals->val;

    lwork = -1L; /* find optimal workspace size */
    dsyev_(&jobz, &uplo, &n, m->val, &n, 
	   w, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	*err = wspace_fail(info, work[0]);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = lapack_realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	*err = E_ALLOC;
	goto bailout;
    } 

    if (!*err) {
	work = work2;
	dsyev_(&jobz, &uplo, &n, m->val, &n, 
	       w, work, &lwork, &info);
	if (info != 0) {
	    *err = 1;
	}
    }

 bailout:

    lapack_free(work);

    if (*err && evals != NULL) {
	gretl_matrix_free(evals);
	evals = NULL;
    }

    return evals;
}

static int gensymm_conformable (const gretl_matrix *A,
				const gretl_matrix *B)
{
    if (!real_gretl_matrix_is_symmetric(A, 1)) {
	fputs("gretl_gensymm_eigenvals: matrix A is not symmetric\n", 
	      stderr);
	return 0;
    }

    if (!real_gretl_matrix_is_symmetric(B, 1)) {
	fputs("gretl_gensymm_eigenvals: matrix B is not symmetric\n", 
	      stderr);
	return 0;
    }

    if (B->rows != A->rows) {
	fputs("gretl_gensymm_eigenvals: matrices A and B have different size\n", 
	      stderr);
	return 0;
    }

    return 1;
}

#define GSDEBUG 0

/**
 * gretl_gensymm_eigenvals:
 * @A: symmetric matrix.
 * @B: symmetric positive definite matrix.
 * @V: matrix to hold the generalized eigenvectors, or %NULL if
 * these are not required.
 * @err: location to receive error code.
 * 
 * Solves the generalized eigenvalue problem
 * | A - \lambda B | = 0 , where both A and B are symmetric
 * and B is positive definite.
 *
 * Returns: allocated storage containing the eigenvalues, or %NULL
 * on failure.
 */

gretl_matrix *gretl_gensymm_eigenvals (const gretl_matrix *A, 
				       const gretl_matrix *B, 
				       gretl_matrix *V, 
				       int *err)
{
    gretl_matrix *K = NULL;
    gretl_matrix *tmp = NULL;
    gretl_matrix *evals = NULL;
    int n;

#if GSDEBUG
    gretl_matrix_print(A, "A");
    gretl_matrix_print(B, "B");
#endif

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    if (!gensymm_conformable(A, B)) {
	*err = E_NONCONF;
	return NULL;
    }
    
    n = A->rows;
    K = gretl_matrix_copy(B);
    tmp = gretl_matrix_alloc(n, n);

    if (K == NULL || tmp == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    *err = gretl_matrix_cholesky_decomp(K);
    if (*err) {
	fputs("gretl_gensymm_eigenvals: matrix B not p.d.\n", 
	      stderr);
	*err = E_NONCONF;
	goto bailout;
    }

    *err = gretl_invert_triangular_matrix(K, 'L');
    if (*err) {
	fputs("gretl_gensymm_eigenvals: matrix B only p.s.d.\n", 
	      stderr);
	*err = E_NONCONF;
	goto bailout;
    }

    gretl_matrix_qform(K, GRETL_MOD_NONE, A, tmp, GRETL_MOD_NONE);

#if GSDEBUG
    gretl_matrix_print(tmp, "tmp");
#endif

    evals = gretl_symmetric_matrix_eigenvals(tmp, 1, err);
    if (*err) {
	goto bailout;
    }

    if (V != NULL) {
	*err = gretl_matrix_multiply_mod(K, GRETL_MOD_TRANSPOSE, 
					 tmp, GRETL_MOD_NONE, 
					 V, GRETL_MOD_NONE);
#if GSDEBUG
	gretl_matrix_print(V, "V");
#endif
    }

 bailout:

    gretl_matrix_free(K);
    gretl_matrix_free(tmp);

    if (*err && evals != NULL) {
	gretl_matrix_free(evals);
	evals = NULL;
    }

    return evals;
}

enum {
    SVD_THIN,
    SVD_FULL
};

/**
 * gretl_matrix_SVD:
 * @a: matrix to decompose.
 * @pu: location for matrix U, or %NULL if not wanted.
 * @ps: location for vector of singular values, or %NULL if not wanted.
 * @pvt: location for matrix V (transposed), or %NULL if not wanted.
 * 
 * Computes SVD factorization of a general matrix using the lapack
 * function %dgesvd. A = u * diag(s) * vt.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

static int 
real_gretl_matrix_SVD (const gretl_matrix *a, gretl_matrix **pu, 
		       gretl_vector **ps, gretl_matrix **pvt,
		       int smod)
{
    integer m, n, lda;
    integer ldu = 1, ldvt = 1;
    integer lwork = -1L;
    integer info;
    gretl_matrix *b = NULL;
    gretl_matrix *s = NULL;
    gretl_matrix *u = NULL;
    gretl_matrix *vt = NULL;
    double xu, xvt;
    char jobu = 'N', jobvt = 'N';
    double *uval = &xu, *vtval = &xvt;
    double *work = NULL, *work2 = NULL;
    int k, err = 0;

    if (pu == NULL && ps == NULL && pvt == NULL) {
	/* no-op */
	return 0;
    }

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    lda = m = a->rows;
    n = a->cols;

    if (smod == SVD_THIN && m < n) {
	fprintf(stderr, "real_gretl_matrix_SVD: a is %d x %d, should be 'thin'\n",
		a->rows, a->cols);
	return E_NONCONF;
    }

    b = gretl_matrix_copy(a);
    if (b == NULL) {
	return E_ALLOC;
    }

    k = (m < n)? m : n;

    s = gretl_vector_alloc(k);
    if (s == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (pu != NULL) {
	if (smod == SVD_FULL) {
	    u = gretl_matrix_alloc(m, m);
	} else {
	    u = gretl_matrix_alloc(m, n);
	}
	if (u == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
	    ldu = m;
	    uval = u->val;
	    jobu = (smod == SVD_FULL)? 'A' : 'S';
	}
    } 

    if (pvt != NULL) {
	vt = gretl_matrix_alloc(n, n);
	if (vt == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	} else {
	    ldvt = n;
	    vtval = vt->val;
	    jobvt = 'A';
	}
    }

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	err = E_ALLOC; 
	goto bailout;
    }	

    /* workspace query */
    dgesvd_(&jobu, &jobvt, &m, &n, b->val, &lda, s->val, uval, &ldu, 
	    vtval, &ldvt, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	err = wspace_fail(info, work[0]);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = lapack_realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	err = E_ALLOC; 
	goto bailout;
    } 
    work = work2;

    /* actual computation */
    dgesvd_(&jobu, &jobvt, &m, &n, b->val, &lda, s->val, uval, &ldu, 
	    vtval, &ldvt, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "gretl_matrix_SVD: info = %d\n", (int) info);
	err = 1;
	goto bailout;
    }

    if (ps != NULL) {
	*ps = s;
	s = NULL;
    }

    if (pu != NULL) {
	*pu = u;
	u = NULL;
    }

    if (pvt != NULL) {
	*pvt = vt;
	vt = NULL;
    }

 bailout:
    
    lapack_free(work);
    gretl_matrix_free(b);
    gretl_matrix_free(s);
    gretl_matrix_free(u);
    gretl_matrix_free(vt);

    return err;
}

/**
 * gretl_matrix_SVD:
 * @a: matrix to decompose.
 * @pu: location for matrix U, or %NULL if not wanted.
 * @ps: location for vector of singular values, or %NULL if not wanted.
 * @pvt: location for matrix V (transposed), or %NULL if not wanted.
 * 
 * Computes SVD factorization of a general matrix using the lapack
 * function %dgesvd. A = u * diag(s) * vt.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_SVD (const gretl_matrix *a, gretl_matrix **pu, 
		      gretl_vector **ps, gretl_matrix **pvt)
{
    return real_gretl_matrix_SVD(a, pu, ps, pvt, SVD_FULL);
}

/**
 * gretl_matrix_SVD_johansen_solve:
 * @R0: T x p matrix of residuals.
 * @R1: T x p1 matrix of residuals.
 * @evals: vector to receive eigenvals, or %NULL if not wanted.
 * @B: matrix to hold \beta, or %NULL if not wanted.
 * @A: matrix to hold \alpha, or %NULL if not wanted.
 * @jrank: cointegration rank, <= p.
 * 
 * Solves the Johnsen generalized eigenvalue problem via
 * SVD decomposition.  See J. A. Doornik and R. J. O'Brien, 
 * "Numerically stable cointegration analysis", Computational 
 * Statistics and Data Analysis, 41 (2002), pp. 185-193, 
 * Algorithm 4.
 *
 * If @B is non-null it should be p1 x p on input; it will
 * be trimmed to p1 x @jrank on output if @jrank < p.
 * If @A is non-null it should be p x p on input; it will
 * be trimmed to p x @jrank on output if @jrank < p.
 * @evals should be a vector of length @jrank.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_SVD_johansen_solve (const gretl_matrix *R0,
				     const gretl_matrix *R1,
				     gretl_matrix *evals,
				     gretl_matrix *B,
				     gretl_matrix *A,
				     int jrank)
{
    gretl_matrix *U0 = NULL;
    gretl_matrix *U1 = NULL;
    gretl_matrix *Uz = NULL;
    gretl_matrix *S1 = NULL;
    gretl_matrix *Sz = NULL;
    gretl_matrix *V1 = NULL;
    gretl_matrix *Z = NULL;
    int T = R0->rows;
    int p = R0->cols;
    int p1 = R1->cols;
    int r, err;

    if (evals == NULL && B == NULL && A == NULL) {
	/* no-op */
	return 0;
    }

    r = (jrank == 0)? p : jrank;

    if (r < 1 || r > p) {
	fprintf(stderr, "Johansen SVD: r is wrong (%d)\n", r);
	return E_NONCONF;
    }

    if (evals != NULL && gretl_vector_get_length(evals) < r) {
	fprintf(stderr, "Johansen SVD: evals is too short\n");
	return E_NONCONF;
    } 

    if (B != NULL && (B->rows != p1 || B->cols != p)) {
	fprintf(stderr, "Johansen SVD: B is wrong size\n");
	return E_NONCONF;
    }

    if (A != NULL && (A->rows != p || A->cols != p)) {
	fprintf(stderr, "Johansen SVD: A is wrong size\n");
	return E_NONCONF;
    }    

    err = real_gretl_matrix_SVD(R0, &U0, NULL, NULL, SVD_THIN);

    if (!err) {
	err = real_gretl_matrix_SVD(R1, &U1, &S1, &V1, SVD_THIN);
    }

    if (!err) {
	Z = gretl_matrix_alloc(p1, p);
	if (Z == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_multiply_mod(U1, GRETL_MOD_TRANSPOSE,
					    U0, GRETL_MOD_NONE,
					    Z, GRETL_MOD_NONE);
	}
    }

    if (!err) {
	err = real_gretl_matrix_SVD(Z, &Uz, &Sz, NULL, SVD_THIN);
    }

    if (!err) {
	double x, si;
	int i, j;
	
	if (evals != NULL) {
	    for (i=0; i<r; i++) {
		evals->val[i] = Sz->val[i] * Sz->val[i];
	    }
	}

	if (B != NULL) {
	    /* \hat{\beta} = T^{1/2} V_1 {\Sigma_1}^{-1} U_z */
	    
	    for (i=0; i<p1; i++) {
		si = S1->val[i];
		for (j=0; j<p1; j++) {
		    if (si > SVD_SMIN) {
			x = gretl_matrix_get(V1, i, j);
			gretl_matrix_set(V1, i, j, x / si);
		    } else {
			gretl_matrix_set(V1, i, j, 0);
		    }
		}
	    }

	    gretl_matrix_multiply_mod(V1, GRETL_MOD_TRANSPOSE,
				      Uz, GRETL_MOD_NONE,
				      B, GRETL_MOD_NONE);
	    gretl_matrix_multiply_by_scalar(B, sqrt((double) T));
	    if (r < p) {
		gretl_matrix_reuse(B, -1, r);
	    }
	}

	if (A != NULL) {
	    /* \hat{\alpha} = T^{-1/2} R_0' U_1 U_z */

	    gretl_matrix_reuse(Z, p, p1);
	    gretl_matrix_multiply_mod(R0, GRETL_MOD_TRANSPOSE,
				      U1, GRETL_MOD_NONE,
				      Z, GRETL_MOD_NONE);
	    gretl_matrix_multiply(Z, Uz, A);
	    gretl_matrix_divide_by_scalar(A, sqrt((double) T));
	    if (r < p) {
		gretl_matrix_reuse(A, -1, r);
	    }
	}
    }

    gretl_matrix_free(U0);
    gretl_matrix_free(U1);
    gretl_matrix_free(Uz);
    gretl_matrix_free(S1);
    gretl_matrix_free(Sz); 
    gretl_matrix_free(V1);
    gretl_matrix_free(Z); 

    return err;
}

/* return the row-index of the element in column col of
   matrix X that has the greatest absolute magnitude
*/

static int max_abs_index (const gretl_matrix *X, int col) 
{
    double aij, tmp = 0.0;
    int i, idx = 0;

    for (i=0; i<X->rows; i++) {
	aij = fabs(gretl_matrix_get(X, i, col));
	if (aij > tmp) {
	    tmp = aij;
	    idx = i;
	}
    }

    return idx;
}

#define NSMIN 1.0e-16

static void normalize_nullspace (gretl_matrix *M)
{
    int i, j, k, idx;
    double x, y; 

    /* FIXME? */

    if (M->cols == 1) {
	j = 0;
	idx = max_abs_index(M, j);
	x = gretl_matrix_get(M, idx, j);
	for (i=0; i<M->rows; i++) {
	    y = gretl_matrix_get(M, i, j);
	    y /= x;
	    if (fabs(y) < NSMIN) y = 0.0;
	    gretl_matrix_set(M, i, j, y);
	}
    }

    /* remove ugliness for printing */
    k = M->rows * M->cols;
    for (i=0; i<k; i++) {
	if (M->val[i] == -0) {
	    M->val[i] = 0;
	}
    }
}

/**
 * gretl_matrix_right_nullspace:
 * @M: matrix to operate on.
 * @err: location to receive error code.
 * 
 * Given an m x n matrix @M, construct a conformable matrix
 * R such that MR = 0 (that is, all the columns of R are
 * orthogonal to the space spanned by the rows of @M).
 *
 * Returns: the allocated matrix R, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_right_nullspace (const gretl_matrix *M, int *err)
{
    gretl_matrix *R = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *S = NULL;
    double x;
    int m, n, r;
    int i, j, k;

    if (gretl_is_null_matrix(M)) {
	
	*err = E_DATA;
	return NULL;
    }

    m = M->rows;
    n = M->cols;
    r = (m < n)? m : n;

    *err = gretl_matrix_SVD(M, NULL, &S, &V);

    if (!*err) {
	/* rank plus nullity = n */
	k = n;
	for (i=0; i<r; i++) {
	    if (S->val[i] > SVD_SMIN) {
		k--;
	    }
	}

	if (k == 0) {
	    R = gretl_null_matrix_new();
	} else {
	    R = gretl_matrix_alloc(n, k);
	}

	if (R == NULL) {
	    *err = E_ALLOC;
	} else if (k > 0) {
	    for (i=0; i<n; i++) {
		for (j=0; j<k; j++) {
		    x = gretl_matrix_get(V, j + n - k, i);
		    gretl_matrix_set(R, i, j, x);
		}
	    }
	    normalize_nullspace(R);
	}
    }

#if 0
    gretl_matrix_print(V, "V'");
    gretl_matrix_print(R, "R");
#endif

    gretl_matrix_free(S);
    gretl_matrix_free(V);

    return R;
}

/**
 * gretl_matrix_left_nullspace:
 * @M: matrix to operate on.
 * @mod: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE
 * @err: location to receive error code.
 * 
 * Given an m x n matrix @M, construct a conformable matrix
 * L such that LM = 0 (that is, all the columns of @M are
 * orthogonal to the space spanned by the rows of L).
 *
 * Returns: the allocated matrix L, or if @mod is
 * %GRETL_MOD_TRANSPOSE, L', or %NULL on failure.
 */

gretl_matrix *gretl_matrix_left_nullspace (const gretl_matrix *M, 
					   GretlMatrixMod mod,
					   int *err)
{
    gretl_matrix *Tmp = NULL;
    gretl_matrix *L = NULL;

    if (gretl_is_null_matrix(M)) {
	*err = E_DATA;
	return NULL;
    }

    Tmp = gretl_matrix_copy_transpose(M);
    if (Tmp == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    L = gretl_matrix_right_nullspace(Tmp, err);
    gretl_matrix_free(Tmp);

    if (!*err && mod == GRETL_MOD_TRANSPOSE) {
	Tmp = gretl_matrix_copy_transpose(L);
	if (Tmp == NULL) {
	    *err = E_ALLOC;
	} else {
	    gretl_matrix_free(L);
	    L = Tmp;
	}
    }

    return L;
}

/**
 * gretl_matrix_row_concat:
 * @a: upper source matrix (m x n).
 * @b: lower source matrix (p x n).
 * @err: location to receive error code.
 * 
 * Returns: newly allocated matrix ((m+p) x n) that results from 
 * the row-wise concatenation of @a and @b, or %NULL on failure.
 */

gretl_matrix *
gretl_matrix_row_concat (const gretl_matrix *a, const gretl_matrix *b,
			 int *err)
{
    gretl_matrix *c = NULL;

    if (a == NULL || b == NULL) {
	*err = 1;
	return NULL;
    }

    if (gretl_is_null_matrix(a) && gretl_is_null_matrix(b)) {
	c = gretl_null_matrix_new();
    } else if (gretl_is_null_matrix(a)) {
	c = gretl_matrix_copy(b);
    } else if (gretl_is_null_matrix(b)) {
	c = gretl_matrix_copy(a);
    } else {
	double x;
	int i, j, k;
	
	if (a->cols != b->cols) {
	    *err = E_NONCONF;
	    return NULL;
	}

	c = gretl_matrix_alloc(a->rows + b->rows, a->cols);

	if (c != NULL) {
	    for (i=0; i<a->rows; i++) {
		for (j=0; j<a->cols; j++) {
		    x = gretl_matrix_get(a, i, j);
		    gretl_matrix_set(c, i, j, x);
		}
	    }  

	    k = a->rows;
	    for (i=0; i<b->rows; i++) {
		for (j=0; j<b->cols; j++) {
		    x = gretl_matrix_get(b, i, j);
		    gretl_matrix_set(c, k, j, x);
		}
		k++;
	    }  
	}
    }

    if (c == NULL) {
	*err = E_ALLOC;
    }    

    return c;
}

/**
 * gretl_matrix_col_concat:
 * @a: left-hand source matrix (m x n).
 * @b: right-hand source matrix (m x p).
 * @err: location to receive error code.
 * 
 * Returns: newly allocated matrix (m x (n+p)) that results from 
 * the column-wise concatenation of @a and @b, or %NULL on failure.
 */

gretl_matrix *
gretl_matrix_col_concat (const gretl_matrix *a, const gretl_matrix *b,
			 int *err)
{
    gretl_matrix *c = NULL;

    if (a == NULL || b == NULL) {
	*err = 1;
	return NULL;
    }

    if (gretl_is_null_matrix(a) && gretl_is_null_matrix(b)) {
	c = gretl_null_matrix_new();
    } else if (gretl_is_null_matrix(a)) {
	c = gretl_matrix_copy(b);
    } else if (gretl_is_null_matrix(b)) {
	c = gretl_matrix_copy(a);
    } else {
	size_t asize, bsize;
	int anelem;

	if (a->rows != b->rows) {
	    *err = E_NONCONF;
	    return NULL;
	}

	c = gretl_matrix_alloc(a->rows, a->cols + b->cols);
	if (c != NULL) {
	    anelem = a->rows * a->cols;
	    asize = anelem * sizeof *a->val;
	    bsize = b->rows * b->cols * sizeof *b->val;

	    memcpy(c->val, a->val, asize);
	    memcpy(c->val + anelem, b->val, bsize);
	}
    }

    if (c == NULL) {
	*err = E_ALLOC;
    }

    return c;
}

/**
 * gretl_matrix_inplace_colcat:
 * @a: matrix to be enlarged (m x n).
 * @b: matrix from which columns should be added (m x p).
 * @mask: char array, of length p, with 1s in positions
 * corresponding to columns of @b that are to be added
 * to @a, 0s elsewhere; or %NULL to add all columns
 * of @b.
 *
 * Concatenates onto @a the selected columns of @b, if the
 * two matrices are conformable.
 * 
 * Returns: 0 on success, non-zero code on error. 
 */

int gretl_matrix_inplace_colcat (gretl_matrix *a, 
				 const gretl_matrix *b,
				 const char *mask)
{
    double x;
    int addc;
    int i, j, k;

    if (a == NULL || b == NULL) {
	return E_DATA;
    }

    if (a->rows != b->rows) {
	return E_NONCONF;
    }

    if (mask == NULL) {
	addc = b->cols;
    } else {
	addc = 0;
	for (j=0; j<b->cols; j++) {
	    if (mask[j]) addc++;
	}
	if (addc == 0) {
	    return 0;
	}
    }

    k = a->cols;

    if (gretl_matrix_realloc(a, a->rows, k + addc)) {
	return E_ALLOC;
    }

    if (mask == NULL) {
	size_t bsize = b->rows * b->cols * sizeof *b->val;

	memcpy(a->val + a->rows * k, b->val, bsize);
    } else {
	for (j=0; j<b->cols; j++) {
	    if (mask[j]) {
		for (i=0; i<b->rows; i++) {
		    x = gretl_matrix_get(b, i, j);
		    gretl_matrix_set(a, i, k, x);
		}
		k++;
	    }
	}
    }

    return 0;
}

/**
 * gretl_matrix_cumcol:
 * @m: source matrix.
 * @err: error code.
 * 
 * Returns: a matrix of the same dimensions as @m, containing 
 * the cumulated columns of @m.
 */

gretl_matrix *gretl_matrix_cumcol (const gretl_matrix *m, int *err)
{
    gretl_matrix *a;
    double x;
    int t, i;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
	return NULL;
    }

    a = gretl_matrix_alloc(m->rows, m->cols);
    if (a == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<m->cols; i++) {
	x = 0;
	for (t=0; t<m->rows; t++) {
	    x += gretl_matrix_get(m, t, i);
	    gretl_matrix_set(a, t, i, x);
	}
    }

    return a;
}

/**
 * gretl_matrix_diffcol:
 * @m: source matrix.
 * @missval: value to represent missing observations.
 * @err: error code.
 * 
 * Returns: a matrix of the same dimensions as @m, containing 
 * @missval in the first row and the difference between consecutive 
 * rows of @m afterwards.
 */

gretl_matrix *gretl_matrix_diffcol (const gretl_matrix *m, 
				    double missval, int *err)
{
    gretl_matrix *a;
    double x, xlag;
    int t, i;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
	return NULL;
    }

    a = gretl_matrix_alloc(m->rows, m->cols);
    if (a == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<m->cols; i++) {
	gretl_matrix_set(a, 0, i, missval);
    }

    for (i=0; i<m->cols; i++) {
	xlag = gretl_matrix_get(m, 0, i);
	for (t=1; t<m->rows; t++) {
	    x = gretl_matrix_get(m, t, i);
	    gretl_matrix_set(a, t, i, x - xlag);
	    xlag = x;
	}
    }

    return a;
}

/**
 * gretl_matrix_lag:
 * @m: source matrix.
 * @k: lag order (> 0 for lags, < 0 for leads).
 * @missval: value to represent missing observations.
 * 
 * Returns: A matrix of the same dimensions as @m, containing lags
 * of the variables in the columns of @m, with missing values set
 * to @missval.
 */

gretl_matrix *gretl_matrix_lag (const gretl_matrix *m, int k, 
				double missval)
{
    gretl_matrix *a;
    double x;
    int s, t, i;

    if (gretl_is_null_matrix(m)) {
	return NULL;
    }

    a = gretl_matrix_alloc(m->rows, m->cols);
    if (a == NULL) {
	return NULL;
    }

    for (t=0; t<m->rows; t++) {
	s = t - k;
	if (s < 0 || s >= m->rows) {
	    for (i=0; i<m->cols; i++) {
		gretl_matrix_set(a, t, i, missval);
	    }
	} else {
	    for (i=0; i<m->cols; i++) {
		x = gretl_matrix_get(m, s, i);
		gretl_matrix_set(a, t, i, x);
	    }
	}
    }

    return a;
}

/**
 * gretl_matrix_inplace_lag:
 * @targ: target matrix.
 * @src: source matrix.
 * @k: lag order (> 0 for lags, < 0 for leads).
 *
 * Fills out @targ (if it is of the correct dimensions),
 * with (columnwise) lags of @src, using 0 for missing
 * values.
 * 
 * Returns: 0 on success, non-zero code otherwise.
 */

int gretl_matrix_inplace_lag (gretl_matrix *targ,
			      const gretl_matrix *src,
			      int k)
{
    int m, n;
    double x;
    int s, t, i;

    if (gretl_is_null_matrix(targ) || gretl_is_null_matrix(src)) {
	return E_DATA;
    }

    m = src->rows;
    n = src->cols;

    if (targ->rows != m || targ->cols != n) {
	return E_NONCONF;
    }

    for (t=0; t<m; t++) {
	s = t - k;
	if (s < 0 || s >= m) {
	    for (i=0; i<n; i++) {
		gretl_matrix_set(targ, t, i, 0.0);
	    }
	} else {
	    for (i=0; i<n; i++) {
		x = gretl_matrix_get(src, s, i);
		gretl_matrix_set(targ, t, i, x);
	    }
	}
    }

    return 0;
}

/**
 * gretl_matrix_set_t1:
 * @m: matrix to operate on.
 * @t: integer value to set.
 * 
 * Sets an integer value on the %t1 member of the gretl_matrix 
 * (used for internal information).  
 */

void gretl_matrix_set_t1 (gretl_matrix *m, int t)
{
    m->t1 = t;
}

/**
 * gretl_matrix_set_t2:
 * @m: matrix to operate on.
 * @t: integer value to set.
 * 
 * Sets an integer value on the %t2 member of the gretl_matrix 
 * (used for internal information).  
 */

void gretl_matrix_set_t2 (gretl_matrix *m, int t)
{
    m->t2 = t;
}

/**
 * gretl_matrix_get_t1:
 * @m: matrix to read from.
 * 
 * Returns: the integer that has been set on the %t1 member
 * of the matrix struct, or zero.
 */

int gretl_matrix_get_t1 (const gretl_matrix *m)
{
    if (m != NULL) {
	return m->t1;
    } else {
	return 0;
    }
}

/**
 * gretl_matrix_get_t2:
 * @m: matrix to read from.
 * 
 * Returns: the integer that has been set on the %t2 member
 * of the matrix struct, or zero.
 */

int gretl_matrix_get_t2 (const gretl_matrix *m)
{
    if (m != NULL) {
	return m->t2;
    } else {
	return 0;
    }
}

static int
get_svd_ols_vcv (const gretl_matrix *A, const gretl_matrix *B,
		 const double *s, gretl_matrix *V, double *s2)
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
		    aik = gretl_matrix_get(A, k, i);
		    ajk = gretl_matrix_get(A, k, j);
		    vij += aik * ajk / (s[k] * s[k]);
		}
	    }
	    gretl_matrix_set(V, i, j, vij);
	    if (j != i) {
		gretl_matrix_set(V, j, i, vij);
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
	gretl_matrix_multiply_by_scalar(V, sigma2);  
	*s2 = sigma2;
    }

    return 0;
}

static double 
get_ols_error_variance (const gretl_vector *y, const gretl_matrix *X,
			const gretl_vector *b, int nr)
{
    double u, s2 = 0.0;
    int k = X->cols;       /* number of regressors */
    int n = X->rows;       /* number of observations */
    int i, j;

    for (i=0; i<n; i++) {
	u = y->val[i];
	for (j=0; j<k; j++) {
	    u -= gretl_matrix_get(X, i, j) * b->val[j];
	}
	s2 += u * u;
    }

    s2 /= (n - k + nr); /* nr = number of restrictions */

    return s2;
}

static int
get_ols_vcv (const gretl_vector *y, const gretl_matrix *X,
	     const gretl_vector *b, gretl_matrix *V,
	     double *s2)
{
    if (gretl_invert_general_matrix(V)) {
	gretl_matrix_print(V, "get_ols_vcv: inversion failed");
	return 1;
    }

    if (s2 != NULL) {
	gretl_matrix_multiply_by_scalar(V, *s2);
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
	    uh -= b->val[j] * gretl_matrix_get(X, i, j);
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
    int T, k;
    integer m, n;
    integer nrhs = 1;
    integer lda, ldb;
    integer lwork = -1L;
    integer rank;
    integer info;
    double rcond = -1.0;
    double *work = NULL;
    double *work2 = NULL;
    double *s = NULL;
    int err = 0;

    if (gretl_is_null_matrix(y) ||
	gretl_is_null_matrix(X) ||
	gretl_is_null_matrix(b)) {
	return E_DATA;
    }

    lda = ldb = m = T = X->rows;
    n = k = X->cols;

    if (gretl_vector_get_length(b) != k) {
	err = E_NONCONF;
	goto bailout;
    }

    A = gretl_matrix_copy(X);
    if (A == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    B = gretl_matrix_copy(y);
    if (B == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* for singular values of A */
    s = malloc(k * sizeof *s);
    if (s == NULL) {
	err = E_ALLOC; 
	goto bailout;
    }

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	err = E_ALLOC; 
	goto bailout;
    } 

    /* workspace query */
    dgelss_(&m, &n, &nrhs, A->val, &lda, B->val, &ldb, s, &rcond,
	    &rank, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	err = wspace_fail(info, work[0]);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = lapack_realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	err = E_ALLOC; 
	goto bailout;
    } 

    work = work2;

    /* get actual solution */
    dgelss_(&m, &n, &nrhs, A->val, &lda, B->val, &ldb, s, &rcond,
	    &rank, work, &lwork, &info);

    if (info != 0) {
	err = 1;
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
    lapack_free(work);
    free(s);

    return err;
}

/**
 * gretl_matrix_multi_SVD_ols:
 * @y: T x g matrix of dependent variables.
 * @X: T x k matrix of independent variables.
 * @B: k x g matrix to hold coefficient estimates.
 * @E: T x g matrix to hold the regression residuals, or %NULL if these are 
 * not needed.
 * @XTXi: location to receive (X'X)^{-1}, or %NULL if this is not needed.
 *
 * Computes OLS estimates using SVD decomposition, and puts the
 * coefficient estimates in @B.  Optionally, calculates the
 * residuals in @E, (X'X)^{-1} in @XTXi.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_multi_SVD_ols (const gretl_matrix *Y, 
				const gretl_matrix *X,
				gretl_matrix *B, 
				gretl_matrix *E,
				gretl_matrix **XTXi)
{
    int g, k, T;
    gretl_matrix *A = NULL;
    gretl_matrix *C = NULL;
    integer m, n, nrhs;
    integer lda, ldb;
    integer lwork = -1L;
    integer rank;
    integer info;
    double rcond = -1.0;
    double *work = NULL;
    double *work2 = NULL;
    double *s = NULL;
    int err = 0;

    /* null_matrix */
    if (gretl_is_null_matrix(Y) ||
	gretl_is_null_matrix(X) ||
	gretl_is_null_matrix(B)) {
	return E_DATA;
    }

    nrhs = g = Y->cols;
    n = k = X->cols;
    lda = ldb = m = T = X->rows;

    if (B->rows != k || B->cols != g) {
	err = E_NONCONF;
    } else if (Y->rows != T) {
	err = E_NONCONF;
    } else if (E != NULL && (E->cols != g || E->rows != T)) {
	err = E_NONCONF;
    } else if (k > T) {
	err = E_DF;
    }

    A = gretl_matrix_copy(X);
    if (A == NULL) {
	return E_ALLOC;
    }

    C = gretl_matrix_copy(Y);
    if (C == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* for singular values of A */
    s = malloc(k * sizeof *s);
    if (s == NULL) {
	err = E_ALLOC; 
	goto bailout;
    }

    work = lapack_malloc(sizeof *work);
    if (work == NULL) {
	err = E_ALLOC; 
	goto bailout;
    } 

    /* workspace query */
    dgelss_(&m, &n, &nrhs, A->val, &lda, C->val, &ldb, s, &rcond,
	    &rank, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
	err = wspace_fail(info, work[0]);
	goto bailout;
    }	

    lwork = (integer) work[0];

    work2 = lapack_realloc(work, lwork * sizeof *work);
    if (work2 == NULL) {
	err = E_ALLOC; 
	goto bailout;
    } 

    work = work2;

    /* get actual solution */
    dgelss_(&m, &n, &nrhs, A->val, &lda, C->val, &ldb, s, &rcond,
	    &rank, work, &lwork, &info);

    if (info != 0) {
	err = 1;
    }

    if (rank < k) {
	fprintf(stderr, "gretl_matrix_svd_ols:\n"
		" dgelss: rank of data matrix X = %d (rows = %d, cols = %d)\n", 
		(int) rank, T, k);
    }

    if (!err) {
	/* coeffs: extract the first k rows from C */
	double bij;
	int i, j;

	for (i=0; i<k; i++) {
	    for (j=0; j<g; j++) {
		bij = gretl_matrix_get(C, i, j);
		gretl_matrix_set(B, i, j, bij);
	    }
	}
    }

    if (!err && E != NULL) {
	/* compute residuals, if wanted */
	int i, imax = E->rows * E->cols;

	gretl_matrix_multiply(X, B, E);
	for (i=0; i<imax; i++) {
	    E->val[i] = Y->val[i] - E->val[i];
	}
    }

    if (!err && XTXi != NULL) {
	/* build (X'X)^{-1}, if wanted */
	*XTXi = gretl_matrix_alloc(k, k);
	if (*XTXi == NULL) {
	    err = E_ALLOC;
	} else {
	    err = get_svd_ols_vcv(A, C, s, *XTXi, NULL);
	}
    }

    bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(C);
    lapack_free(work);
    free(s);

    return err;
}

/**
 * gretl_matrix_moore_penrose:
 * @a: m x n matrix.
 * 
 * Computes the generalized inverse of matrix @a via its SVD
 * factorization, with the help of the lapack function 
 * %dgesvd.  On exit the original matrix is overwritten by 
 * the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_moore_penrose (gretl_matrix *A)
{
    gretl_matrix *U = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *Vt = NULL;
    gretl_matrix *SUt = NULL;
    double x;
    int m, n;
    int i, j;
    int err = 0;

    if (gretl_is_null_matrix(A)) {
	return E_DATA;
    }

    m = A->rows;
    n = A->cols;

    err = gretl_matrix_SVD(A, &U, &S, &Vt);

    if (!err) {
	int nsv = (m < n)? m : n;

	SUt = gretl_zero_matrix_new(n, m);
	if (SUt == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}

	/* invert singular values and multiply into U' */
	for (i=0; i<nsv; i++) {
	    if (S->val[i] > SVD_SMIN) {
		for (j=0; j<m; j++) {
		    x = gretl_matrix_get(U, j, i);
		    gretl_matrix_set(SUt, i, j, x / S->val[i]);
		}
	    }
	}

	/* A^{+} = VS^{-1}U' */
	A->rows = n;
	A->cols = m;
	err = gretl_matrix_multiply_mod(Vt, GRETL_MOD_TRANSPOSE,
					SUt, GRETL_MOD_NONE,
					A, GRETL_MOD_NONE);
    }

 bailout:
    
    gretl_matrix_free(U);
    gretl_matrix_free(S);
    gretl_matrix_free(Vt);
    gretl_matrix_free(SUt);

    return err;
}

/**
 * gretl_SVD_invert_matrix:
 * @a: n x n matrix to invert.
 * 
 * Computes the inverse (or generalized inverse) of a general square 
 * matrix using SVD factorization, with the help of the lapack function 
 * %dgesvd.  If any of the singular values of @a are less than 1.0e-9
 * the Moore-Penrose generalized inverse is computed instead of the
 * standard inverse.  On exit the original matrix is overwritten by 
 * the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_SVD_invert_matrix (gretl_matrix *a)
{
    gretl_matrix *u = NULL;
    gretl_matrix *s = NULL;
    gretl_matrix *vt = NULL;
    double x;
    int i, j, k, n;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
	return E_DATA;
    }

    if (a->rows != a->cols) {
	err = E_NONCONF;
	goto bailout;
    }

    n = a->rows;

    /* a = USV' ; a^{-1} = VWU' where W holds inverse of diag elements of S */

    err = gretl_matrix_SVD(a, &u, &s, &vt);

    if (!err) {
	k = 0;
	for (i=0; i<n; i++) {
	    if (s->val[i] < SVD_SMIN) {
		break;
	    }
	    k++;
	}
	if (k < n) {
	    gretl_matrix *vt2;

	    fprintf(stderr, "gretl_SVD_invert_matrix: rank = %d (dim = %d)\n", 
		    k, (int) n);
	    fputs("Warning: computing Moore-Penrose generalized inverse\n", stderr);

	    vt2 = gretl_matrix_alloc(k, n);
	    if (vt2 == NULL) {
		err = E_ALLOC;
		goto bailout;
	    }
	    for (i=0; i<k; i++) {
		for (j=0; j<n; j++) {
		    x = gretl_matrix_get(vt, i, j);
		    gretl_matrix_set(vt2, i, j, x);
		}
	    }
	    gretl_matrix_free(vt);
	    vt = vt2;
	    gretl_matrix_reuse(u, n, k);
	}	    
    }

    if (!err) {
	/* invert singular values */
	for (j=0; j<k; j++) {
	    for (i=0; i<n; i++) {
		x = gretl_matrix_get(u, i, j);
		gretl_matrix_set(u, i, j, x / s->val[j]);
	    }
	}
	err = gretl_matrix_multiply_mod(vt, GRETL_MOD_TRANSPOSE,
					u, GRETL_MOD_TRANSPOSE,
					a, GRETL_MOD_NONE);
    }

 bailout:
    
    gretl_matrix_free(u);
    gretl_matrix_free(s);
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
 * if @s2 is %NULL, the "vcv" estimate will be plain (X'X)^{-1}.
 *
 * Computes OLS estimates using LU factorization, and puts the
 * coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv and the residuals in @uhat.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_ols (const gretl_vector *y, const gretl_matrix *X,
		      gretl_vector *b, gretl_matrix *vcv, 
		      gretl_vector *uhat, double *s2)
{
    gretl_matrix *XTX = NULL;
    int k, err = 0;

    if (gretl_is_null_matrix(y) ||
	gretl_is_null_matrix(X) ||
	gretl_is_null_matrix(b)) {
	return E_DATA;
    }

    k = X->cols;

    if (gretl_vector_get_length(b) != k) {
	err = E_NONCONF;
    }

    if (vcv != NULL && (vcv->rows != k || vcv->cols != k)) {
	err = E_NONCONF;
    }    

    if (!err) {
	XTX = gretl_matrix_XTX_new(X);
	if (XTX == NULL) err = E_ALLOC;
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					y, GRETL_MOD_NONE,
					b, GRETL_MOD_NONE);
    }

    if (!err && vcv != NULL) {
	err = gretl_matrix_copy_values(vcv, XTX);
    }

    if (!err) {
#if 0
	err = gretl_cholesky_decomp_solve(XTX, b);
#else
	err = gretl_LU_solve(XTX, b);
#endif
    }

    if (!err) {
	if (s2 != NULL) {
	    *s2 = get_ols_error_variance(y, X, b, 0);
	}
	if (vcv != NULL) {
	    err = get_ols_vcv(y, X, b, vcv, s2);
	}
	if (uhat != NULL) {
	    get_ols_uhat(y, X, b, uhat);
	}
    }

    if (XTX != NULL) gretl_matrix_free(XTX);

    return err;
}

/**
 * gretl_matrix_multi_ols:
 * @y: T x g matrix of dependent variables.
 * @X: T x k matrix of independent variables.
 * @B: k x g matrix to hold coefficient estimates.
 * @E: T x g matrix to hold the regression residuals, or %NULL if these are 
 * not needed.
 * @XTXi: location to receive (X'X)^{-1} on output, or %NULL if this is not
 * needed.
 *
 * Computes OLS estimates using Cholesky decomposition, and puts the
 * coefficient estimates in @B.  Optionally, calculates the
 * residuals in @E.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_multi_ols (const gretl_matrix *Y, const gretl_matrix *X,
			    gretl_matrix *B, gretl_matrix *E,
			    gretl_matrix **XTXi)
{
    gretl_matrix *XTX = NULL;
    int g, T, k;
    int err = 0;

    if (gretl_is_null_matrix(Y) ||
	gretl_is_null_matrix(X) ||
	gretl_is_null_matrix(B)) {
	return E_DATA;
    }

    g = Y->cols;
    T = X->rows;
    k = X->cols;

    if (B->rows != k || B->cols != g) {
	err = E_NONCONF;
    } else if (Y->rows != T) {
	err = E_NONCONF;
    } else if (E != NULL && (E->cols != g || E->rows != T)) {
	err = E_NONCONF;
    } else if (k > T) {
	err = E_DF;
    } 

    if (!err) {
	XTX = gretl_matrix_XTX_new(X);
	if (XTX == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					Y, GRETL_MOD_NONE,
					B, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_cholesky_decomp_solve(XTX, B);
    }

    if (!err && E != NULL) {
	gretl_matrix_multiply(X, B, E);
	gretl_matrix_switch_sign(E);
	gretl_matrix_add_to(E, Y);
    }

    if (!err && XTXi != NULL) {
	integer info = 0, ik = k;
	char uplo = 'L';

	dpotri_(&uplo, &ik, XTX->val, &ik, &info);
	gretl_symmetric_matrix_expand(XTX, uplo);
	*XTXi = XTX;
    } else {
	gretl_matrix_free(XTX);
    }

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
 * @uhat: vector to hold residuals, if wanted.
 * @s2: pointer ro receive residual variance, or NULL.  If vcv is non-NULL
 * and s2 is NULL, the vcv estimate is just W^{-1}.
 *
 * Computes OLS estimates restricted by R and q, using LU factorization, 
 * and puts the coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int 
gretl_matrix_restricted_ols (const gretl_vector *y, const gretl_matrix *X,
			     const gretl_matrix *R, const gretl_vector *q,
			     gretl_vector *b, gretl_matrix *vcv,
			     gretl_vector *uhat, double *s2)
{
    gretl_matrix *XTX = NULL;
    gretl_vector *V = NULL;
    gretl_matrix *W = NULL;
    gretl_matrix *S = NULL;
    double x;
    int k = X->cols;
    int nr = R->rows;
    int ldW = k + nr;
    int err = 0;
    int i, j;

    if (gretl_vector_get_length(b) != k) {
	fprintf(stderr, "gretl_matrix_restricted_ols: "
		"b should be a %d-vector\n", k);
	err = E_NONCONF;
    }

    if (!err) {
	XTX = gretl_matrix_XTX_new(X);
	V = gretl_column_vector_alloc(ldW);
	W = gretl_zero_matrix_new(ldW, ldW);
	if (XTX == NULL || V == NULL || W == NULL) {
	    err = E_ALLOC;
	}
    }

    /* construct V matrix: X'y augmented by q (or by a 0
       matrix if q is NULL)
    */

    if (!err) {
	V->rows = k;
	err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
					y, GRETL_MOD_NONE,
					V, GRETL_MOD_NONE);
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
	for (i=0; i<XTX->rows; i++) {
	    for (j=0; j<XTX->cols; j++) {
		x = gretl_matrix_get(XTX, i, j);
		gretl_matrix_set(W, i, j, x);
	    }
	}
	for (i=0; i<R->rows; i++) {
	    for (j=0; j<R->cols; j++) {
		x = gretl_matrix_get(R, i, j);
		gretl_matrix_set(W, i+k, j, x);
	    }
	}
	for (i=0; i<R->cols; i++) {
	    for (j=0; j<R->rows; j++) {
		x = gretl_matrix_get(R, j, i);
		gretl_matrix_set(W, i, j+k, x);
	    }
	}
    } 

    if (!err && vcv != NULL) {
	S = gretl_matrix_copy(W);
	if (S == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_LU_solve(W, V);
    }

    if (!err) {
	int i;
	
	for (i=0; i<k; i++) {
	    b->val[i] = V->val[i];
	}
	if (s2 != NULL) {
	    *s2 = get_ols_error_variance(y, X, b, nr);
	}
	if (S != NULL) {
	    err = get_ols_vcv(y, X, b, S, s2);
	    if (!err) {
		for (i=0; i<k; i++) {
		    for (j=0; j<k; j++) {
			x = gretl_matrix_get(S, i, j);
			gretl_matrix_set(vcv, i, j, x);
		    }
		}		
	    }
	    gretl_matrix_free(S);
	}
	if (uhat != NULL) {
	    get_ols_uhat(y, X, b, uhat);
	}
    }

    if (XTX != NULL) gretl_matrix_free(XTX);
    if (V != NULL) gretl_vector_free(V);
    if (W != NULL) gretl_matrix_free(W);

    return err;
}

static int QR_OLS_work (gretl_matrix *Q, gretl_matrix *R)
{
    integer k = gretl_matrix_rows(R);
    int r, err;

    /* basic decomposition */
    err = gretl_matrix_QR_decomp(Q, R);
    if (err) {
	return err;
    }

    /* check rank of QR */
    r = gretl_check_QR_rank(R, &err);
    if (err) {
	return err;
    }

    if (r < k) {
	err = E_SINGULAR;
    } else {
	/* invert R */
	char uplo = 'U';
	char diag = 'N';
	integer info = 0;

	dtrtri_(&uplo, &diag, &k, R->val, &k, &info);
	if (info != 0) {
	    fprintf(stderr, "dtrtri: info = %d\n", (int) info);
	    err = 1;
	}
    }

    return err;
}

/**
 * gretl_matrix_QR_ols:
 * @y: T x g matrix of dependent variables.
 * @X: T x k matrix of independent variables.
 * @B: k x g matrix to hold coefficient estimates.
 * @E: T x g matrix to hold the regression residuals, or %NULL if these are 
 * not needed.
 * @XTXi: location to receive (X'X)^{-1}, or %NULL if this is not needed.
 * @Qout: location to receive Q on output, or %NULL.
 *
 * Computes OLS estimates using QR decomposition, and puts the
 * coefficient estimates in @B.  Optionally, calculates the
 * residuals in @E, (X'X)^{-1} in @XTXi, and/or the matrix Q
 * in @Qout.
 * 
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_QR_ols (const gretl_matrix *Y,
			 const gretl_matrix *X,
			 gretl_matrix *B,
			 gretl_matrix *E,
			 gretl_matrix **XTXi,
			 gretl_matrix **Qout)
{
    int g = Y->cols;
    int k = X->cols;
    int T = X->rows;
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *G = NULL;
    int err = 0;

    if (B->rows != k || B->cols != g) {
	err = E_NONCONF;
    } else if (Y->rows != T) {
	err = E_NONCONF;
    } else if (E != NULL && (E->cols != g || E->rows != T)) {
	err = E_NONCONF;
    } else if (k > T) {
	err = E_DF;
    }

    if (!err) {
	Q = gretl_matrix_copy(X);
	R = gretl_matrix_alloc(k, k);
	G = gretl_matrix_alloc(k, g);
	if (Q == NULL || R == NULL || G == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = QR_OLS_work(Q, R);
    }

    if (!err) {
	/* make "G" into gamma-hat */    
	gretl_matrix_multiply_mod(Q, GRETL_MOD_TRANSPOSE,
				  Y, GRETL_MOD_NONE, 
				  G, GRETL_MOD_NONE);
    }

    if (!err) {
	/* OLS coefficients */
	gretl_matrix_multiply(R, G, B);
    }

    if (!err && E != NULL) {
	/* compute residuals */
	int i, imax = E->rows * E->cols;

	gretl_matrix_multiply(X, B, E);
	for (i=0; i<imax; i++) {
	    E->val[i] = Y->val[i] - E->val[i];
	}
    }

    /* create (X'X)^{-1} = RR' */
    if (!err && XTXi != NULL) {
	*XTXi = gretl_matrix_alloc(k, k);

	if (*XTXi == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
				      R, GRETL_MOD_TRANSPOSE,
				      *XTXi, GRETL_MOD_NONE);
	}
    }

    if (!err && Qout != NULL) {
	*Qout = Q;
    } else {
	gretl_matrix_free(Q);
    }

    gretl_matrix_free(R);
    gretl_matrix_free(G);

    return err;    
}

/**
 * gretl_matrix_r_squared:
 * @y: dependent variable, T-vector.
 * @X: independent variables matrix, T x k.
 * @b: coefficients, k-vector.
 * @err: location to receive error code.
 *
 * Returns: the unadjusted R-squared, based one the regression 
 * represented by @y, @X and @b, or #NADBL on failure.
 */

double gretl_matrix_r_squared (const gretl_matrix *y,
			       const gretl_matrix *X,
			       const gretl_matrix *b,
			       int *err)
{
    double ess = 0.0, tss = 0.0;
    double xx, ybar;
    int i, j;

    if (gretl_vector_get_length(y) != X->rows ||
	gretl_vector_get_length(b) != X->cols) {
	*err = E_NONCONF;
	return NADBL;
    }

    ybar = gretl_vector_mean(y);

    for (i=0; i<X->rows; i++) {
	xx = y->val[i];
	for (j=0; j<X->cols; j++) {
	    xx -= b->val[j] * gretl_matrix_get(X, i, j);
	}
	ess += xx * xx;
	xx = y->val[i] - ybar;
	tss += xx * xx;
    }

    return 1.0 - ess / tss;
}

/**
 * gretl_matrix_columwise_product:
 * @A: T x k matrix.
 * @B: T x n matrix.
 * @C: T x (k*n) matrix to hold the product.
 *
 * Computes a columnwise product in k blocks, each of n columns.
 * The first block consists of the Hadamard product of the first
 * column of @A and the matrix @B, the second block holds the
 * Hadamard product of the second column of @A and matrix @B, 
 * and so on.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_columnwise_product (const gretl_matrix *A,
				     const gretl_matrix *B,
				     gretl_matrix *C)
{
    int k, n, T;
    double x, y;
    int i, j, t, p;

    if (gretl_is_null_matrix(A) ||
	gretl_is_null_matrix(B) ||
	gretl_is_null_matrix(C)) {
	return E_DATA;
    }

    k = A->cols;
    n = B->cols;
    T = A->rows;

    if (B->rows != T || C->rows != T) {
	return E_NONCONF;
    }

    if (C->cols != k * n) {
	return E_NONCONF;
    }

    p = 0;
    for (i=0; i<k; i++) {
	for (j=0; j<n; j++) {
	    for (t=0; t<T; t++) {
		x = gretl_matrix_get(A, t, i);
		y = gretl_matrix_get(B, t, j);
		gretl_matrix_set(C, t, p, x * y);
	    }
	    p++;
	}
    }
	    
    return 0;
}

/**
 * gretl_matrix_qform:
 * @A: m * k matrix or k * m matrix, depending on @amod.
 * @amod: %GRETL_MOD_NONE or %GRETL_MOD_TRANSPOSE: in the first
 * case @A should be m * k; in the second, k * m;
 * @X: symmetric k * k matrix.
 * @C: matrix to hold the product.
 * @cmod: modifier: %GRETL_MOD_NONE or %GRETL_MOD_CUMULATE to
 * add the result to the existing value of @C.
 *
 * Computes either A * X * A' (if amod = %GRETL_MOD_NONE) or
 * A' * X * A (if amod = %GRETL_MOD_TRANSPOSE), with the result 
 * written into @C.  The matrix @X must be symmetric, but this
 * is not checked, to save time.  If you are in doubt on this
 * point you can call gretl_matrix_is_symmetric() first.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

#define QFORM_SMALL 1.0e-20 

int gretl_matrix_qform (const gretl_matrix *A, GretlMatrixMod amod,
			const gretl_matrix *X, gretl_matrix *C, 
			GretlMatrixMod cmod)
{
    register int i, j, ii, jj;
    double xi, xj, xij, xx;
    int m, k;

    if (gretl_is_null_matrix(A) ||
	gretl_is_null_matrix(X) ||
	gretl_is_null_matrix(C)) {
	return E_DATA;
    }

    m = (amod)? A->cols : A->rows;
    k = (amod)? A->rows : A->cols;

    if (k != X->rows) {
	fprintf(stderr, "gretl_matrix_qform: %s is (%d x %d) but X is (%d x %d)\n", 
		(amod)? "A'" : "A", m, k, X->rows, X->cols);
	return E_NONCONF;
    }

    if (C->rows != m || C->cols != m) {
	fputs("gretl_matrix_qform: destination matrix not conformable\n", stderr);
	return E_NONCONF;
    }

    if (amod) {
	for (i=0; i<m; i++) {
	    for (j=i; j<m; j++) {
		xx = 0.0;
		for (ii=0; ii<k; ii++) {
		    xi = gretl_matrix_get(A,ii,i);
		    if (fabs(xi) > QFORM_SMALL) {
			for (jj=0; jj<k; jj++) {
			    xj = gretl_matrix_get(A,jj,j);
			    xij = gretl_matrix_get(X,ii,jj);
			    xx += xij * xi * xj;
			}
		    }
		}
		if (cmod == GRETL_MOD_CUMULATE) {
		    xx += gretl_matrix_get(C, i, j);
		} else if (cmod == GRETL_MOD_DECUMULATE) {
		    xx = gretl_matrix_get(C, i, j) - xx;
		}
		gretl_matrix_set(C, i, j, xx);
		gretl_matrix_set(C, j, i, xx);
	    }
	}
    } else {
	for (i=0; i<m; i++) {
	    for (j=i; j<m; j++) {
		xx = 0.0;
		for (ii=0; ii<k; ii++) {
		    xi = gretl_matrix_get(A,i,ii);
		    if (fabs(xi) > QFORM_SMALL) {
			for (jj=0; jj<k; jj++) {
			    xj = gretl_matrix_get(A,j,jj);
			    xij = gretl_matrix_get(X,ii,jj);
			    xx += xij * xi * xj;
			}
		    }
		}
		if (cmod == GRETL_MOD_CUMULATE) {
		    xx += gretl_matrix_get(C, i, j);
		} else if (cmod == GRETL_MOD_DECUMULATE) {
		    xx = gretl_matrix_get(C, i, j) - xx;
		}
		gretl_matrix_set(C, i, j, xx);
		gretl_matrix_set(C, j, i, xx);
	    }
	}
    }

    return 0;
}

/**
 * gretl_scalar_qform:
 * @b: k-vector.
 * @X: symmetric k x k matrix.
 * @errp: pointer to receive error code.
 *
 * Computes the scalar product bXb', or b'Xb if @b is a column
 * vector.  The content of @errp is set to 0 on success,
 * or a non-zero code on failure.
 * 
 * Returns: the scalar product, or #NADBL on failure.
 */

double gretl_scalar_qform (const gretl_vector *b, 
			   const gretl_matrix *X,
			   int *errp)
{
    gretl_matrix *tmp = NULL;
    double ret = NADBL;
    int mod, k;

    if (gretl_is_null_matrix(b) || gretl_is_null_matrix(X)) {
	*errp = E_DATA;
	return NADBL;
    } 

    k = gretl_vector_get_length(b);

    if (k == 0 || X->rows != k || X->cols != k) {
	*errp = E_NONCONF;
	return NADBL;
    } 

    mod = (b->rows > 1)? GRETL_MOD_TRANSPOSE : GRETL_MOD_NONE;

    tmp = gretl_matrix_alloc(1, 1);
    if (tmp == NULL) {
	*errp = E_ALLOC;
    } else {
	tmp->val[0] = 0.0;
	*errp = gretl_matrix_qform(b, mod, X, tmp, GRETL_MOD_NONE);
	if (!*errp) {
	    ret = tmp->val[0];
	}
	gretl_matrix_free(tmp);
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
    double x, xij;
    int i, j, err = 0;

    if (dim != X->rows || dim != X->cols ||
	dim != DXD->rows || dim != DXD->cols) {
	err = E_NONCONF;
    } else {
	for (i=0; i<dim; i++) {
	    for (j=0; j<dim; j++) {
		xij = gretl_matrix_get(X, i, j);
		x = xij * d->val[i] * d->val[j];
		gretl_matrix_set(DXD, i, j, x);
	    }
	}
    }

    return err;
}

/**
 * gretl_is_identity_matrix:
 * @m: matrix to examine.
 *
 * Returns: 1 if @m is an identity matrix, 0 otherwise.
 */

int gretl_is_identity_matrix (const gretl_matrix *m)
{
    double x;
    int i, j;

    if (gretl_is_null_matrix(m)) {
	return 0;
    }

    for (j=0; j<m->cols; j++) {
	for (i=0; i<m->rows; i++) {
	    x = gretl_matrix_get(m, i, j);
	    if (i == j && x != 1.0) return 0;
	    if (i != j && x != 0.0) return 0;
	}
    }

    return 1;
}

/**
 * gretl_is_zero_matrix:
 * @m: matrix to examine.
 *
 * Returns: 1 if @m is a zero matrix, 0 otherwise.
 */

int gretl_is_zero_matrix (const gretl_matrix *m)
{
    int i, n;

    if (gretl_is_null_matrix(m)) {
	return 0;
    }

    n = m->rows * m->cols;

    for (i=0; i<n; i++) {
	if (m->val[i] != 0.0) return 0;
    }

    return 1;
}

/**
 * gretl_matrices_are_equal:
 * @a: first matrix in comparison.
 * @b: second matrix in comparison.
 * @err: location to receive error code.
 *
 * Returns: 1 if the matrices @a and @b compare equal, 0 if they
 * differ, and -1 if the comparison is invalid, in which case
 * %E_NONCONF is written to @err.
 */

int gretl_matrices_are_equal (const gretl_matrix *a, const gretl_matrix *b,
			      int *err)
{
    double ax, bx;
    int i, j;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	*err = E_DATA;
	return -1;
    }

    if (a->rows != b->rows || a->cols != b->cols) {
	*err = E_NONCONF;
	return -1;
    }

    for (i=0; i<a->rows; i++) {
	for (j=0; j<a->cols; j++) {
	    ax = gretl_matrix_get(a, i, j);
	    bx = gretl_matrix_get(b, i, j);
	    if (ax != bx) {
		fprintf(stderr, "gretl_matrices_are_equal:\n "
			"a(%d,%d) = %.15g but b(%d,%d) = %.15g\n",
			i, j, ax, i, j, bx);
		return 0;
	    }
	}
    }

    return 1;
}

/* if pxbar, pssx are non-NULL, they get the vectors of column
   means and sums of squares */

static gretl_matrix *
real_gretl_covariance_matrix (const gretl_matrix *m, int corr,
			      gretl_matrix **pxbar, gretl_matrix **pssx,
			      int *errp)
{
    gretl_matrix *V = NULL;
    gretl_vector *xbar = NULL;
    gretl_vector *ssx = NULL;
    int k, n;
    int t, i, j;
    double vv, x, y;
    int err = 0;

    if (gretl_is_null_matrix(m)) {
	err = E_DATA;
	goto bailout;
    }	
    
    if (errp != NULL) {
	*errp = 0;
    }

    k = m->cols;
    n = m->rows;

    if (n < 2) {
	err = E_DATA;
	goto bailout;
    }

    V = gretl_matrix_alloc(k, k);
    xbar = gretl_vector_alloc(k);

    if (V == NULL || xbar == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (corr) {
	ssx = gretl_vector_alloc(k);
	if (ssx == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    for (i=0; i<k; i++) {
	xbar->val[i] = 0.0;
	for (t=0; t<n; t++) {
	    xbar->val[i] += gretl_matrix_get(m, t, i);
	}
	xbar->val[i] /= n;
    }

    if (ssx != NULL) {
	for (i=0; i<k; i++) {
	    ssx->val[i] = 0.0;
	    for (t=0; t<n; t++) {
		x = gretl_matrix_get(m, t, i) - xbar->val[i];
		ssx->val[i] += x * x;
	    }
	}
    }	

    for (i=0; i<k; i++) {
	for (j=i; j<k; j++) {
	    vv = 0.0;
	    for (t=0; t<n; t++) {
		x = gretl_matrix_get(m, t, i);
		y = gretl_matrix_get(m, t, j);
		vv += (x - xbar->val[i]) * (y - xbar->val[j]);
	    }
	    if (ssx != NULL) {
		if (vv != 0.0) {
		    x = ssx->val[i] * ssx->val[j];
		    vv /= sqrt(x);
		}
	    } else {
		vv /= n - 1;
	    }
	    gretl_matrix_set(V, i, j, vv);
	    gretl_matrix_set(V, j, i, vv);
	}
    }

 bailout:

    if (!err && pxbar != NULL) {
	*pxbar = xbar;
    } else {
	gretl_vector_free(xbar);
    }

    if (!err && pssx != NULL) {
	*pssx = ssx;
    } else {
	gretl_vector_free(ssx);
    }

    if (err) {
	if (errp != NULL) {
	    *errp = err;
	}
	gretl_matrix_free(V);
	V = NULL;
    }

    return V;
}

/**
 * gretl_covariance_matrix:
 * @m: (x x n) matrix containing n observations on each of k 
 * variables.
 * @corr: flag for computing correlations.
 * @errp: pointer to receive non-zero error code in case of
 * failure, or %NULL.
 *
 * Returns: the covariance matrix of variables in the columns of
 * @m, or the correlation matrix if @corr is non-zero.
 */

gretl_matrix *gretl_covariance_matrix (const gretl_matrix *m, int corr,
				       int *errp)
{
    return real_gretl_covariance_matrix(m, corr, NULL, NULL, errp);
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
 * gretl_matrix_values:
 * @x: array to process.
 * @n: length of array.
 * @err: location to receive error code.
 *
 * Returns: an allocated matrix containing the distinct
 * values in array @x, or %NULL on failure. 
 */

gretl_matrix *gretl_matrix_values (const double *x, int n,
				   int *err)
{
    gretl_matrix *v = NULL;
    int *sorted = NULL;
    int i, k, m, last;

    sorted = malloc(n * sizeof *sorted);
    if (sorted == NULL) {
	*err = E_ALLOC;
	return NULL;
    }
	
    k = 0;
    for (i=0; i<n; i++) {
	if (!na(x[i])) {
	    sorted[k++] = (int) x[i];
	}
    }

    if (k == 0) {
	*err = E_DATA;
	goto bailout;
    }
	
    qsort(sorted, k, sizeof *sorted, gretl_compare_ints); 
    m = count_distinct_int_values(sorted, k);

    v = gretl_column_vector_alloc(m);
    if (v == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    v->val[0] = last = sorted[0];
    m = 1;

    for (i=1; i<k; i++) {
	if (sorted[i] != last) {
	    last = sorted[i];
	    v->val[m++] = sorted[i];
	}
    }

 bailout:

    free(sorted);

    return v;
}

/**
 * gretl_matrix_shape:
 * @A: array to process.
 * @r: rows of target matrix.
 * @c: columns of target matrix.
 *
 * Creates an (r x c) matrix containing the re-arranged
 * values of A.  Elements are read from A by column and
 * written into the target, also by column.  If A contains 
 * less elements than n = r*c, they are repeated cyclically; 
 * if A has more elements, only the first n are used.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_shape (const gretl_matrix *A, 
				  int r, int c)
{
    gretl_matrix *B;
    int i, k, nA, nB;

    if (r <= 0 || c <= 0 || gretl_is_null_matrix(A)) {
	return NULL;
    }
    
    B = gretl_matrix_alloc(r, c);
    if (B == NULL) {
	return NULL;
    }

    nA = A->rows * A->cols;
    nB = r * c;

    k = 0;
    for (i=0; i<nB; i++) {
	B->val[i] = A->val[k++];
	if (k == nA) {
	    k = 0;
	}
    }

    return B;
}

/**
 * gretl_matrix_trim_rows:
 * @A: array to process.
 * @ttop: rows to trim at top.
 * @tbot: rows to trim at bottom.
 * @err: location to receive error code.
 *
 * Creates a new matrix which is a copy of @A with @ttop rows 
 * trimmed from the top and @tbot rows trimmed from the
 * bottom.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_trim_rows (const gretl_matrix *A, 
				      int ttop, int tbot,
				      int *err)
{
    gretl_matrix *B;
    double x;
    int i, j, m;

    if (gretl_is_null_matrix(A)) {
	*err = E_DATA;
	return NULL;
    }

    m = A->rows - (ttop + tbot);

    if (ttop < 0 || tbot < 0 || m <= 0) {
	*err = E_DATA;
	return NULL;
    }
    
    B = gretl_matrix_alloc(m, A->cols);
    if (B == NULL) {
	return NULL;
    }

    for (j=0; j<A->cols; j++) {
	for (i=0; i<m; i++) {
	    x = gretl_matrix_get(A, i + ttop, j);
	    gretl_matrix_set(B, i, j, x);
	}
    }

    return B;
}

/**
 * gretl_matrix_minmax:
 * @A: m x n matrix to examine.
 * @mm: 0 for minimum, 1 for maximum.
 * @rc: 0 for row, 1 for column.
 * @idx: 0 for values, 1 for indices.
 *
 * Creates a matrix holding the row or column mimima or
 * maxima from @A, either as values or as location indices.
 * For example, if @mm = 0, @rc = 0, and @idx = 0, the
 * created matrix is m x 1 and holds the values of the row
 * minima.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_minmax (const gretl_matrix *A, 
				   int mm, int rc, int idx,
				   int *err)
{
    gretl_matrix *B;
    double d, x;
    int i, j, k;
    
    if (gretl_is_null_matrix(A)) {
	*err = E_DATA;
	return NULL;
    }

    if (rc == 0) {
	B = gretl_matrix_alloc(A->rows, 1);
    } else {
	B = gretl_matrix_alloc(1, A->cols);
    }

    if (B == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (rc == 0) {
	for (i=0; i<A->rows; i++) {
	    d = gretl_matrix_get(A, i, 0);
	    k = 0;
	    for (j=1; j<A->cols; j++) {
		x = gretl_matrix_get(A, i, j);
		if ((mm > 0 && x > d) || (mm == 0 && x < d)) {
		    d = x;
		    k = j;
		} 
	    }
	    B->val[i] = (idx > 0)? (double) k + 1 : d;
	}
    } else {
	for (j=0; j<A->cols; j++) {
	    d = gretl_matrix_get(A, 0, j);
	    k = 0;
	    for (i=1; i<A->rows; i++) {
		x = gretl_matrix_get(A, i, j);
		if ((mm > 0 && x > d) || (mm == 0 && x < d)) {
		    d = x;
		    k = i;
		} 
	    }
	    B->val[j] = (idx > 0)? (double) k + 1 : d;
	}
    }	

    return B;
}

/**
 * gretl_matrix_pca:
 * @X: T x m data matrix.
 * @p: number of principal components to return: 0 < p <= m,
 * or p = -1 to return components for eigenvalues > 1.0.
 * @err: location to receive error code.
 *
 * Carries out a Principal Components analysis of @X and
 * returns the first @p components: the component corresponding
 * to the largest eigenvalue of the correlation matrix of @X
 * is placed in column 1, and so on.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_pca (const gretl_matrix *X, int p, int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *P = NULL;
    gretl_matrix *xbar = NULL;
    gretl_matrix *ssx = NULL;
    gretl_matrix *evals = NULL;
    int T, m;
    double x, load, val;
    int i, j, k;

    if (gretl_is_null_matrix(X)) {
	*err = E_DATA;
	return NULL;
    }

    T = X->rows;
    m = X->cols;

    if (m == 1) {
	/* match wit to wit */
	P = gretl_matrix_copy(X);
	if (P == NULL) {
	    *err = E_ALLOC;
	}
	return P;
    }

    if (p <= 0) {
	p = 1;
    } else if (p > m) {
	p = m;
    }

    C = real_gretl_covariance_matrix(X, 1, &xbar, &ssx, err);
    if (*err) {
	return NULL;
    }

    evals = gretl_symmetric_matrix_eigenvals(C, 1, err);
    if (*err) {
	goto bailout;
    }

    gretl_symmetric_eigen_sort(evals, C, p);

    /* make matrix to contain the first p components */
    P = gretl_matrix_alloc(T, p);
    if (P == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* convert ssx to std deviations */
    for (i=0; i<m; i++) {
	x = ssx->val[i];
	ssx->val[i] = sqrt(x / (T - 1));
    }

    /* compute the PCs */
    for (j=0; j<p; j++) {
	for (i=0; i<T; i++) {
	    x = 0.0;
	    for (k=0; k<m; k++) {
		load = gretl_matrix_get(C, k, j);
		val = gretl_matrix_get(X, i, k);
		x += load * (val - xbar->val[k]) / ssx->val[k];
	    }
	    gretl_matrix_set(P, i, j, x);
	}
    }

 bailout:
    
    gretl_matrix_free(xbar);
    gretl_matrix_free(ssx);
    gretl_matrix_free(C);
    gretl_matrix_free(evals);

    return P;
}

#define complete_obs(x,y,t) (!na(x[t]) && !na(y[t]))

static int ok_xy_count (int t1, int t2, const double *x, const double *y)
{
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (complete_obs(x, y, t)) {
	    n++;
	}
    }
    
    return n;
}

static void make_matrix_xtab (double **X, int n, 
			      const gretl_matrix *vx,
			      const gretl_matrix *vy,
			      gretl_matrix *tab)
{
    int xr, xc, rndx, cndx;
    int counter, i;

    qsort(X, n, sizeof *X, compare_xtab_rows);

    /* compute frequencies by going through sorted X */

    counter = rndx = cndx = 0;
    xr = (int) gretl_vector_get(vx, 0);
    xc = (int) gretl_vector_get(vy, 0);

    for (i=0; i<n; i++) {
	while (X[i][0] > xr) { 
	    /* skip row */
	    gretl_matrix_set(tab, rndx, cndx, counter);
	    counter = 0;
	    xr = gretl_vector_get(vx, ++rndx);
	    cndx = 0;
	    xc = gretl_vector_get(vy, 0);
	}
	while (X[i][1] > xc) { 
	    /* skip column */
	    gretl_matrix_set(tab, rndx, cndx, counter);
	    counter = 0;
	    xc = gretl_vector_get(vy, ++cndx);
	}
	counter++;
    }
    gretl_matrix_set(tab, rndx, cndx, counter);
}

/**
 * matrix_matrix_xtab:
 * @x: data vector
 * @y: data vector
 * @err: error code
 *
 * Computes the cross tabulation of the values contained in the
 * vectors x (by row) and y (by column). These must be integer values.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *matrix_matrix_xtab (const gretl_matrix *x,
				  const gretl_matrix *y,
				  int *err)
{
    gretl_matrix *tab = NULL;
    gretl_matrix *vx = NULL;
    gretl_matrix *vy = NULL;
    double **X = NULL;
    int i, nx, ny;

    *err = 0;

    nx = gretl_vector_get_length(x);
    ny = gretl_vector_get_length(y);

    if (nx < 2 || ny != nx) {
	*err = E_NONCONF;
	return NULL;
    }

    vx = gretl_matrix_values(x->val, nx, err);
    if (*err) {
	return NULL;
    }

    vy = gretl_matrix_values(y->val, ny, err);
    if (*err) {
	goto bailout;
    }

    tab = gretl_zero_matrix_new(vx->rows, vy->rows);
    if (tab == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    X = doubles_array_new(nx, 2);
    if (X == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<nx; i++) {
	X[i][0] = (int) x->val[i];
	X[i][1] = (int) y->val[i];
    }

    make_matrix_xtab(X, nx, vx, vy, tab);

 bailout:

    gretl_matrix_free(vx);
    gretl_matrix_free(vy);
    doubles_array_free(X, nx);

    return tab;
}

/**
 * gretl_matrix_xtab:
 * @x: data vector
 * @y: data vector
 * @t1: start
 * @t2: end
 * @err: error code
 *
 * Computes the cross tabulation of the values contained in the
 * vectors x (by row) and y (by column). These must be integer values.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_xtab (int t1, int t2, const double *x, 
				 const double *y, int *err)
{
    gretl_matrix *tab = NULL;
    gretl_matrix *vx = NULL;
    gretl_matrix *vy = NULL;
    double *tmp = NULL;
    double **X = NULL;
    int i, t, nmax = t2 - t1 + 1;

    *err = 0;

    nmax = ok_xy_count(t1, t2, x, y);
    if (nmax < 2) {
	*err = E_MISSDATA;
	return NULL;
    }

    tmp = malloc(nmax * sizeof *tmp);
    if (tmp == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    i = 0;
    for (t=t1; t<=t2; t++) {
	if (complete_obs(x, y, t)) {
	    tmp[i++] = x[t];
	}
    }

    vx = gretl_matrix_values(tmp, nmax, err);
    if (*err) {
	free(tmp);
	return NULL;
    }

    i = 0;
    for (t=t1; t<=t2; t++) {
	if (complete_obs(x, y, t)) {
	    tmp[i++] = y[t];
	}
    }

    vy = gretl_matrix_values(tmp, nmax, err);
    if (*err) {
	goto bailout;
    }

    tab = gretl_zero_matrix_new(gretl_matrix_rows(vx), 
				gretl_matrix_rows(vy));
    if (tab == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    X = doubles_array_new(nmax, 2);
    if (X == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    i = 0;
    for (t=t1; t<=t2; t++) {
	if (complete_obs(x, y, t)) {
	    X[i][0] = (int) x[t];
	    X[i][1] = (int) y[t];
	    i++;
	}
    }

    make_matrix_xtab(X, nmax, vx, vy, tab);

 bailout:

    free(tmp);
    gretl_matrix_free(vx);
    gretl_matrix_free(vy);
    doubles_array_free(X, nmax);

    return tab;
}

/**
 * gretl_matrix_bool_sel:
 * @A: matrix.
 * @sel: selection vector.
 * @rowsel: row/column mode selector.
 * @err: location to receive error code.
 *
 * If @rowsel = 1, constructs a matrix which contains the rows 
 * of A corresponding to non-zero values in the vector @sel; 
 * if @rowsel = 0, does the same thing but column-wise.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_bool_sel (const gretl_matrix *A, 
				     const gretl_matrix *sel, 
				     int rowsel, int *err)
{
    gretl_matrix *ret = NULL;
    int ra, ca, rs, cs;
    int nonzero = 0; 
    int i, j, k, n;
    double x;

    *err = 0;

    if (gretl_is_null_matrix(A)) {
	return gretl_null_matrix_new();
    }

    ra = A->rows;
    ca = A->cols;
    rs = sel->rows;
    cs = sel->cols;

    /* check dimensions */

    if (rowsel) {
	if ((ra != rs) || (cs > 1)) {
	    *err = E_NONCONF;
	    return NULL;
	}
    } else {
	if ((ca != cs) || (rs > 1)) {
	    *err = E_NONCONF;
	    return NULL;
	}
    }

    /* count nonzeros */

    n = (rowsel)? rs : cs ;
    for (i=0; i<n; i++) {
	if (gretl_vector_get(sel, i) != 0) {
	    nonzero++;
	}
    }

    /* check for extreme cases */

    if (nonzero == n) {
	ret = gretl_matrix_copy(A);
	goto bailout;
    } else if (nonzero == 0) {
	ret = gretl_null_matrix_new();
	goto bailout;
    } 

    /* copy selected row/columns */
	
    if (rowsel) {
	ret = gretl_matrix_alloc(nonzero, ca);
	if (ret == NULL) {
	    goto bailout;
	}
	k = 0;
	for (i=0; i<ra; i++) {
	    if (gretl_vector_get(sel, i) != 0) {
		for (j=0; j<ca; j++) {
		    x = gretl_matrix_get(A, i, j);
		    gretl_matrix_set(ret, k, j, x);
		}
		k++;
	    }
	}
    } else {
	ret = gretl_matrix_alloc(ra, nonzero);
	if (ret == NULL) {
	    goto bailout;
	}
	for (i=0; i<ra; i++) {
	    k = 0;
	    for (j=0; j<ca; j++) {
		if (gretl_vector_get(sel, j) != 0) {
		    x = gretl_matrix_get(A, i, j);
		    gretl_matrix_set(ret, i, k++, x);
		}
	    }
	}
    }

 bailout:

    if (ret == NULL) {
	*err = E_ALLOC;
    }

    return ret;
}

/**
 * gretl_matrix_sort_by_column:
 * @m: matrix.
 * @k: column by which to sort.
 * @err: location to receive error code.
 *
 * Produces a matrix which contains the rows of @m, re-
 * ordered by increasing value of the elements in column
 * @k.  
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_sort_by_column (const gretl_matrix *m, 
					   int k, int *err)
{
    struct rsort {
	double x;
	int row;
    } *rs; 
    gretl_matrix *a;
    double x;
    int i, j;

    if (gretl_is_null_matrix(m) || k < 0 || k >= m->cols) {
	*err = E_DATA;
	return NULL;
    }

    rs = malloc(m->rows * sizeof *rs);
    if (rs == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    a = gretl_matrix_copy(m);
    if (a == NULL) {
	free(rs);
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<m->rows; i++) {
	rs[i].x = gretl_matrix_get(m, i, k);
	rs[i].row = i;
    }

    qsort(rs, m->rows, sizeof *rs, gretl_compare_doubles);

    for (j=0; j<m->cols; j++) {
	for (i=0; i<m->rows; i++) {
	    x = gretl_matrix_get(m, rs[i].row, j);
	    gretl_matrix_set(a, i, j, x);
	}
    }

    free(rs);

    return a;
}
