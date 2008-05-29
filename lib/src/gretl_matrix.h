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

#ifndef GRETL_MATRIX_H
#define GRETL_MATRIX_H

/* #define LDEBUG 1 */

/* minimum value of diagonal element of R (as in X = QR) that counts
   as non-zero for the purpose of determining the rank of X */

#define R_DIAG_MIN 1.0e-8

typedef enum {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE,
    GRETL_MOD_SQUARE,
    GRETL_MOD_CUMULATE,
    GRETL_MOD_DECUMULATE
} GretlMatrixMod;

typedef enum {
    GRETL_MATRIX_SQUARE = 1,
    GRETL_MATRIX_LOWER_TRIANGULAR,
    GRETL_MATRIX_UPPER_TRIANGULAR,
    GRETL_MATRIX_SYMMETRIC,
    GRETL_MATRIX_DIAGONAL,
    GRETL_MATRIX_IDENTITY,
    GRETL_MATRIX_SCALAR,
} GretlMatrixStructure;

typedef struct _gretl_matrix gretl_matrix;
typedef struct _gretl_matrix gretl_vector;

struct _gretl_matrix {
    int rows;
    int cols;
    int t1, t2;
    double *val;
};

/**
 * gretl_matrix_get:
 * @m: matrix.
 * @i: row.
 * @j: column.
 * 
 * Retrieves the @i, @j element of @m.
 */

#define gretl_matrix_get(m,i,j) (m->val[(j)*m->rows+(i)])

/**
 * gretl_vector_get:
 * @v: vector.
 * @i: index.
 * 
 * Gives element @i of @v.
 */

#define gretl_vector_get(v,i) (v->val[i])

/**
 * gretl_matrix_set:
 * @m: matrix.
 * @i: row.
 * @j: column.
 * @x: value to set.
 * 
 * Sets the @i, @j element of @m to @x.
 */

#define gretl_matrix_set(m,i,j,x) ((m)->val[(j)*(m)->rows+(i)]=x)

/**
 * gretl_vector_set:
 * @v: vector.
 * @i: index.
 * @x: value to set.
 * 
 * Sets element @i of @v to @x.
 */

#define gretl_vector_set(v,i,x) ((v)->val[i]=x)

/**
 * gretl_matrix_cols:
 * @m: matrix to query.
 * 
 * Gives the number of columns in @m. 
 */

#define gretl_matrix_cols(m) ((m == NULL)? 0 : m->cols)

/**
 * gretl_matrix_rows:
 * @m: matrix to query.
 * 
 * Gives the number of rows in @m. 
 */

#define gretl_matrix_rows(m) ((m == NULL)? 0 : m->rows)

/**
 * gretl_vector_get_length:
 * @v: vector to examine.
 * 
 * Gives the length of vector @v (without regard to whether
 * it is a row or column vector).
 */

#define gretl_vector_get_length(v) ((v == NULL)? 0 : \
                                    (v->cols == 1)? v->rows : \
                                    (v->rows == 1)? v->cols : 0)

/**
 * gretl_vector_alloc:
 * @i: number of columns.
 *
 * Allocates a new #gretl_vector with @i columns.
 */

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))

/**
 * gretl_column_vector_alloc:
 * @i: number of rows.
 *
 * Allocates a new column gretl_vector with @i rows.
 */

#define gretl_column_vector_alloc(i) gretl_matrix_alloc((i),1)

/**
 * gretl_vector_free:
 * @v: %gretl_vector to free.
 *
 * Frees the vector @v and its associated storage.
 */

#define gretl_vector_free(v) gretl_matrix_free(v)

/**
 * gretl_matrix_is_scalar:
 * @m: matrix to test.
 *
 * Gives 1 if @m is 1 x 1, else 0.
 */

#define gretl_matrix_is_scalar(m) ((m) != NULL && \
                                   (m)->rows == 1 && \
                                   (m)->cols == 1)

#define gretl_is_null_matrix(m) (m == NULL || m->rows == 0 || m->cols == 0)

int get_gretl_matrix_err (void);

void clear_gretl_matrix_err (void);

void gretl_matrix_print (const gretl_matrix *m, const char *msg);

int gretl_matrix_is_finite (const gretl_matrix *m);

int gretl_matrix_is_symmetric (const gretl_matrix *m);

int gretl_matrix_is_idempotent (const gretl_matrix *m);

void gretl_matrix_xtr_symmetric (gretl_matrix *m);

void gretl_matrix_set_equals_tolerance (double tol);

void gretl_matrix_unset_equals_tolerance (void);

gretl_matrix *gretl_matrix_alloc (int rows, int cols);

gretl_matrix *gretl_matrix_reuse (gretl_matrix *m, int rows, int cols);

int gretl_matrix_realloc (gretl_matrix *m, int rows, int cols);

gretl_matrix *gretl_identity_matrix_new (int n);

gretl_matrix *gretl_zero_matrix_new (int r, int c);

gretl_matrix *gretl_unit_matrix_new (int r, int c);

gretl_matrix *gretl_null_matrix_new (void);

gretl_matrix *gretl_matrix_seq (int start, int end);

gretl_matrix *gretl_matrix_copy (const gretl_matrix *m);

int gretl_matrix_copy_row (gretl_matrix *dest, int di,
			   const gretl_matrix *src, int si);

int gretl_matrix_inscribe_I (gretl_matrix *m, int row, int col, int n);

gretl_matrix *gretl_matrix_copy_transpose (const gretl_matrix *m);

gretl_matrix *gretl_matrix_get_diagonal (const gretl_matrix *m, int *err);

double gretl_matrix_trace (const gretl_matrix *m, int *err);

int gretl_matrix_random_fill (gretl_matrix *m, int dist);

gretl_matrix *gretl_random_matrix_new (int r, int c, int dist);

gretl_matrix *gretl_matrix_resample (const gretl_matrix *m, int r,
				     int *err);

double gretl_vector_mean (const gretl_vector *v);

double gretl_vector_variance (const gretl_vector *v);

void gretl_matrix_zero (gretl_matrix *m);

int gretl_matrix_zero_upper (gretl_matrix *m);

int gretl_matrix_zero_lower (gretl_matrix *m);

void gretl_matrix_fill (gretl_matrix *m, double x);

void gretl_matrix_multiply_by_scalar (gretl_matrix *m, double x);

int gretl_matrix_divide_by_scalar (gretl_matrix *m, double x);

void gretl_matrix_switch_sign (gretl_matrix *m);

gretl_matrix *
gretl_matrix_dot_op (const gretl_matrix *a, const gretl_matrix *b,
		     int op, int *err);

gretl_matrix *gretl_matrix_complex_multiply (const gretl_matrix *a, 
					     const gretl_matrix *b,
					     int *err);

gretl_matrix *gretl_matrix_complex_divide (const gretl_matrix *a, 
					   const gretl_matrix *b,
					   int *err);

gretl_matrix *gretl_matrix_exp (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_polroots (const gretl_matrix *a,
				     int *err);

void gretl_matrix_raise (gretl_matrix *m, double x);

void gretl_matrix_free (gretl_matrix *m);

double *gretl_matrix_steal_data (gretl_matrix *m);

int gretl_matrix_copy_values (gretl_matrix *targ, 
			      const gretl_matrix *src);

int gretl_matrix_copy_values_shaped (gretl_matrix *targ, 
				     const gretl_matrix *src);

int gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src);

int gretl_matrix_add_transpose_to (gretl_matrix *targ, 
				   const gretl_matrix *src);

int 
gretl_matrix_subtract_from (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_subtract_reversed (const gretl_matrix *a, gretl_matrix *b);

int gretl_matrix_I_minus (gretl_matrix *m);

int gretl_matrix_transpose (gretl_matrix *targ, const gretl_matrix *src);

int gretl_square_matrix_transpose (gretl_matrix *m);

int gretl_matrix_add_self_transpose (gretl_matrix *m);

int 
gretl_matrix_vectorize (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_unvectorize (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_vectorize_h (gretl_matrix *targ, const gretl_matrix *src);

int 
gretl_matrix_unvectorize_h (gretl_matrix *targ, const gretl_matrix *src);

int gretl_matrix_inscribe_matrix (gretl_matrix *targ,
				  const gretl_matrix *src,
				  int row, int col,
				  GretlMatrixMod mod);

int gretl_matrix_extract_matrix (gretl_matrix *targ,
				 const gretl_matrix *src,
				 int row, int col,
				 GretlMatrixMod mod);

int gretl_matrix_multiply_mod (const gretl_matrix *a, GretlMatrixMod amod,
			       const gretl_matrix *b, GretlMatrixMod bmod,
			       gretl_matrix *c, GretlMatrixMod cmod);

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c);

gretl_matrix *gretl_matrix_multiply_new (const gretl_matrix *a, 
					 const gretl_matrix *b,
					 int *err);

int
gretl_matrix_kronecker_product (const gretl_matrix *A, const gretl_matrix *B,
				gretl_matrix *K);

gretl_matrix *
gretl_matrix_kronecker_product_new (const gretl_matrix *A, 
				    const gretl_matrix *B,
				    int *err);

int
gretl_matrix_I_kronecker (int p, const gretl_matrix *B,
			  gretl_matrix *K);

gretl_matrix *
gretl_matrix_I_kronecker_new (int p, const gretl_matrix *B, int *err);

int
gretl_matrix_kronecker_I (const gretl_matrix *A, int r,
			  gretl_matrix *K);

gretl_matrix *
gretl_matrix_kronecker_I_new (const gretl_matrix *A, int r, int *err);

gretl_matrix *gretl_matrix_pow (const gretl_matrix *A, 
				int s, int *err);

double gretl_matrix_dot_product (const gretl_matrix *a, GretlMatrixMod amod,
				 const gretl_matrix *b, GretlMatrixMod bmod,
				 int *errp);

double gretl_vector_dot_product (const gretl_vector *a, const gretl_vector *b,
				 int *errp);

gretl_matrix *gretl_matrix_row_sum (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_column_sum (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_row_mean (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_column_mean (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_column_sd (const gretl_matrix *m, int *err);

double gretl_matrix_row_i_mean (const gretl_matrix *m, int row);

double gretl_matrix_column_j_mean (const gretl_matrix *m, int col);

void gretl_matrix_demean_by_row (gretl_matrix *m);

void gretl_matrix_demean_by_column (gretl_matrix *m);

gretl_matrix *gretl_matrix_vcv (gretl_matrix *m);

gretl_matrix *gretl_matrix_quantiles (const gretl_matrix *m,
				      double p, int *err);

double gretl_matrix_determinant (gretl_matrix *a, int *err);

double gretl_matrix_log_determinant (gretl_matrix *a, int *err);

double gretl_matrix_log_abs_determinant (gretl_matrix *a, int *err);

double gretl_vcv_log_determinant (const gretl_matrix *m);

double gretl_matrix_one_norm (const gretl_matrix *m);

double gretl_matrix_infinity_norm (const gretl_matrix *m);

int gretl_LU_solve (gretl_matrix *a, gretl_matrix *b);

int gretl_cholesky_decomp_solve (gretl_matrix *a, gretl_matrix *b);

int gretl_cholesky_solve (const gretl_matrix *a, gretl_vector *b);

int gretl_inverse_from_cholesky_decomp (gretl_matrix *targ, 
					const gretl_matrix *src);

int gretl_invert_general_matrix (gretl_matrix *a);

int gretl_invert_symmetric_indef_matrix (gretl_matrix *a);

int gretl_invert_symmetric_matrix (gretl_matrix *a);

int gretl_invert_symmetric_matrix2 (gretl_matrix *a, double *ldet);

int gretl_invert_packed_symmetric_matrix (gretl_matrix *v);

int gretl_invert_triangular_matrix (gretl_matrix *a, char uplo);

int gretl_invert_diagonal_matrix (gretl_matrix *a);

int gretl_invert_matrix (gretl_matrix *a);

int gretl_matrix_moore_penrose (gretl_matrix *A);

int gretl_SVD_invert_matrix (gretl_matrix *a);

int gretl_invpd (gretl_matrix *a);

int gretl_matrix_SVD (const gretl_matrix *a, gretl_matrix **pu, 
		      gretl_vector **ps, gretl_matrix **pvt);

double gretl_symmetric_matrix_rcond (const gretl_matrix *m, int *err);

double gretl_matrix_rcond (const gretl_matrix *m, int *err);

int gretl_general_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs, 
			      int rank);

int gretl_symmetric_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs, 
				int rank);

gretl_matrix *
gretl_general_matrix_eigenvals (gretl_matrix *m,
				int eigenvecs, 
				int *err);

gretl_matrix *
gretl_symmetric_matrix_eigenvals (gretl_matrix *m,
				  int eigenvecs, 
				  int *err);

gretl_matrix *
gretl_gensymm_eigenvals (const gretl_matrix *A, 
			 const gretl_matrix *B, 
			 gretl_matrix *V, 
			 int *err);

gretl_matrix *gretl_matrix_right_nullspace (const gretl_matrix *M, 
					    int *err);

gretl_matrix *gretl_matrix_left_nullspace (const gretl_matrix *M,
					   GretlMatrixMod mod,
					   int *err);

gretl_matrix *
gretl_matrix_row_concat (const gretl_matrix *a, const gretl_matrix *b,
			 int *err);

gretl_matrix *
gretl_matrix_col_concat (const gretl_matrix *a, const gretl_matrix *b,
			 int *err);

int
gretl_matrix_inplace_colcat (gretl_matrix *a, const gretl_matrix *b,
			     const char *mask);

gretl_matrix *gretl_matrix_cumcol (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_diffcol (const gretl_matrix *m, 
				    double missval, int *err);

gretl_matrix *gretl_matrix_lag (const gretl_matrix *m, int k, 
				double missval);

int gretl_matrix_inplace_lag (gretl_matrix *targ,
			      const gretl_matrix *src,
			      int k);

int gretl_matrix_cholesky_decomp (gretl_matrix *a);

int gretl_matrix_QR_decomp (gretl_matrix *M, gretl_matrix *R);

int gretl_check_QR_rank (const gretl_matrix *R, int *err);

int gretl_matrix_rank (const gretl_matrix *a, int *err);

int gretl_matrix_ols (const gretl_vector *y, const gretl_matrix *X,
		      gretl_vector *b, gretl_matrix *vcv,
		      gretl_vector *uhat, double *s2);

int gretl_matrix_multi_ols (const gretl_matrix *Y, 
			    const gretl_matrix *X,
			    gretl_matrix *B, 
			    gretl_matrix *E,
			    gretl_matrix **XTXi);

int gretl_matrix_multi_SVD_ols (const gretl_matrix *Y, 
				const gretl_matrix *X,
				gretl_matrix *B, 
				gretl_matrix *E,
				gretl_matrix **XTXi);

int gretl_matrix_QR_ols (const gretl_matrix *Y,
			 const gretl_matrix *X,
			 gretl_matrix *B,
			 gretl_matrix *E,
			 gretl_matrix **XTXi,
			 gretl_matrix **Qout);

double gretl_matrix_r_squared (const gretl_matrix *y,
			       const gretl_matrix *X,
			       const gretl_matrix *b,
			       int *err);

int gretl_matrix_SVD_johansen_solve (const gretl_matrix *R0,
				     const gretl_matrix *R1,
				     gretl_matrix *evals,
				     gretl_matrix *B,
				     gretl_matrix *A,
				     int jrank);

int 
gretl_matrix_restricted_ols (const gretl_vector *y, const gretl_matrix *X,
			     const gretl_matrix *R, const gretl_vector *q,
			     gretl_vector *b, gretl_matrix *vcv,
			     gretl_vector *uhat, double *s2);

int gretl_matrix_svd_ols (const gretl_vector *y, const gretl_matrix *X,
			  gretl_vector *b, gretl_matrix *vcv,
			  gretl_vector *uhat, double *s2);

int gretl_matrix_qform (const gretl_matrix *A, GretlMatrixMod amod,
			const gretl_matrix *X, gretl_matrix *C, 
			GretlMatrixMod cmod);

double gretl_scalar_qform (const gretl_vector *b, 
			   const gretl_matrix *X,
			   int *errp);

int gretl_matrix_columnwise_product (const gretl_matrix *A,
				     const gretl_matrix *B,
				     gretl_matrix *C);

int
gretl_matrix_diagonal_sandwich (const gretl_vector *d, const gretl_matrix *X,
				gretl_matrix *DXD);

void gretl_matrix_set_t1 (gretl_matrix *m, int t);

void gretl_matrix_set_t2 (gretl_matrix *m, int t);

int gretl_matrix_get_t1 (const gretl_matrix *m);

int gretl_matrix_get_t2 (const gretl_matrix *m);

int gretl_is_identity_matrix (const gretl_matrix *m);

int gretl_is_zero_matrix (const gretl_matrix *m);

int gretl_matrix_get_structure (const gretl_matrix *m);

int gretl_matrices_are_equal (const gretl_matrix *a, const gretl_matrix *b,
			      int *err);

gretl_matrix *gretl_covariance_matrix (const gretl_matrix *m, int corr,
				       int *errp);

gretl_matrix **gretl_matrix_array_alloc (int n);

gretl_matrix **
gretl_matrix_array_alloc_with_size (int n, int rows, int cols);

void gretl_matrix_array_free (gretl_matrix **A, int n);

gretl_matrix *gretl_matrix_values (const double *x, int n,
				   int *err);

gretl_matrix *gretl_matrix_shape (const gretl_matrix *A, 
				  int r, int c);

gretl_matrix *gretl_matrix_trim_rows (const gretl_matrix *A, 
				      int ttop, int tbot,
				      int *err);

gretl_matrix *gretl_matrix_minmax (const gretl_matrix *A, 
				   int mm, int rc, int idx,
				   int *err);

gretl_matrix *gretl_matrix_pca (const gretl_matrix *X, int p, int *err);

gretl_matrix *gretl_matrix_xtab (int t1, int t2, const double *x, 
				 const double *y, int *err);

gretl_matrix *matrix_matrix_xtab (const gretl_matrix *x,
				  const gretl_matrix *y,
				  int *err);

gretl_matrix *gretl_matrix_bool_sel(const gretl_matrix *A, 
				    const gretl_matrix *sel, 
				    int rowsel, int *err);

gretl_matrix *gretl_matrix_sort_by_column (const gretl_matrix *m, 
					   int k, int *err);

void lapack_mem_free (void);

#endif /* GRETL_MATRIX_H */
