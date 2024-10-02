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

#ifdef  __cplusplus
extern "C" {
#endif

#include <stdarg.h>
#include <complex.h>

/* minimum value of diagonal element of R (as in X = QR) that counts
   as non-zero for the purpose of determining the rank of X */

#define R_DIAG_MIN 1.0e-8

typedef enum {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE,
    GRETL_MOD_SQUARE,
    GRETL_MOD_CUMULATE,
    GRETL_MOD_DECREMENT,
    GRETL_MOD_CTRANSP
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

typedef enum {
    V_SUM,
    V_PROD,
    V_MEAN
} GretlVecStat;

typedef enum {
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
} ConfType;

typedef struct gretl_matrix_ gretl_vector;

typedef struct matrix_info_ matrix_info;

/**
 * gretl_matrix:
 * @rows: number of rows in matrix
 * @cols: number of columns
 * @val: flat array of double-precision values
 *
 * The basic libgretl matrix type; #gretl_vector is an alias
 * that can be used for matrices with @rows or @cols = 1.
 */

typedef struct gretl_matrix_ {
    int rows;
    int cols;
    double *val;
    double _Complex *z; /* was "complex" */
    int is_complex;
    /*< private >*/
    matrix_info *info;
} gretl_matrix;

typedef struct gretl_matrix_block_ gretl_matrix_block;

/**
 * gretl_matrix_get:
 * @m: matrix.
 * @i: row.
 * @j: column.
 *
 * Returns: the @i, @j element of @m.
 */

#define gretl_matrix_get(m,i,j) (m->val[(j)*m->rows+(i)])
#define gretl_cmatrix_get(m,i,j) (m->z[(j)*m->rows+(i)])

/**
 * gretl_vector_get:
 * @v: vector.
 * @i: index.
 *
 * Returns: element @i of @v.
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
#define gretl_cmatrix_set(m,i,j,x) ((m)->z[(j)*(m)->rows+(i)]=x)

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
 * Returns: the number of columns in @m.
 */

#define gretl_matrix_cols(m) ((m == NULL)? 0 : m->cols)

/**
 * gretl_matrix_rows:
 * @m: matrix to query.
 *
 * Returns: the number of rows in @m.
 */

#define gretl_matrix_rows(m) ((m == NULL)? 0 : m->rows)

/**
 * gretl_vector_get_length:
 * @v: vector to examine.
 *
 * Returns: the length of vector @v (without regard to whether
 * it is a row or column vector).
 */

#define gretl_vector_get_length(v) ((v == NULL)? 0 : \
				    ((v)->cols == 1)? (v)->rows :	\
				    ((v)->rows == 1)? (v)->cols : 0)

/**
 * gretl_vector_alloc:
 * @i: number of columns.
 *
 * Returns: a new #gretl_vector with @i columns.
 */

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))

/**
 * gretl_column_vector_alloc:
 * @i: number of rows.
 *
 * Returns: a new column gretl_vector with @i rows.
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
 * Returns: 1 if @m is 1 x 1, else 0.
 */

#define gretl_matrix_is_scalar(m) ((m) != NULL && \
				   (m)->is_complex == 0 && \
                                   (m)->rows == 1 && \
                                   (m)->cols == 1)

#define gretl_matrix_is_cscalar(m) ((m) != NULL && \
				    (m)->is_complex && \
                                    (m)->rows == 1 && \
                                    (m)->cols == 1)

#define gretl_is_null_matrix(m) (m == NULL || m->rows == 0 || m->cols == 0)

#define gretl_is_complex(m) (m != NULL && m->is_complex == 1)

int gretl_matrix_set_complex (gretl_matrix *m, int c);

int gretl_matrix_set_complex_full (gretl_matrix *m, int c);

int get_gretl_matrix_err (void);

void clear_gretl_matrix_err (void);

void gretl_matrix_print (const gretl_matrix *m, const char *msg);

void gretl_matrix_print2 (const gretl_matrix *m, const char *msg);

int gretl_matrix_na_check (const gretl_matrix *m);

int gretl_matrix_is_symmetric (const gretl_matrix *m);

int gretl_matrix_is_idempotent (const gretl_matrix *m, double tol);

void gretl_matrix_xtr_symmetric (gretl_matrix *m);

void gretl_matrix_set_equals_tolerance (double tol);

void gretl_matrix_unset_equals_tolerance (void);

gretl_matrix *gretl_matrix_alloc (int rows, int cols);

gretl_matrix *gretl_cmatrix_new (int r, int c);

gretl_matrix *gretl_cmatrix_new0 (int r, int c);

gretl_matrix *gretl_cmatrix_new1 (int r, int c);

gretl_matrix *gretl_matching_matrix_new (int r, int c,
					 const gretl_matrix *m);

gretl_matrix *gretl_matrix_reuse (gretl_matrix *m, int rows, int cols);

int gretl_matrix_realloc (gretl_matrix *m, int rows, int cols);

gretl_matrix *gretl_matrix_init (gretl_matrix *m);

gretl_matrix *gretl_matrix_init_full (gretl_matrix *m,
				      int rows, int cols,
				      double *val);

gretl_matrix *gretl_matrix_replace (gretl_matrix **pa,
				    gretl_matrix *b);

int gretl_matrix_replace_content (gretl_matrix *targ,
				  gretl_matrix *donor);

void gretl_matrix_block_destroy (gretl_matrix_block *B);

void gretl_matrix_block_zero (gretl_matrix_block *B);

gretl_matrix_block *gretl_matrix_block_new (gretl_matrix **pm, ...);

int gretl_matrix_block_n_matrices (gretl_matrix_block *B);

gretl_matrix *gretl_matrix_block_get_matrix (gretl_matrix_block *B,
					     int i);

gretl_matrix *gretl_identity_matrix_new (int n);

gretl_matrix *gretl_DW_matrix_new (int n);

gretl_matrix *gretl_zero_matrix_new (int r, int c);

gretl_matrix *gretl_unit_matrix_new (int r, int c);

gretl_matrix *gretl_null_matrix_new (void);

gretl_matrix *gretl_matrix_seq (double start, double end,
				double step, int *err);

gretl_matrix *gretl_matrix_copy (const gretl_matrix *m);

int gretl_matrix_copy_row (gretl_matrix *dest, int di,
			   const gretl_matrix *src, int si);

int gretl_matrix_inscribe_I (gretl_matrix *m, int row, int col, int n);

gretl_matrix *gretl_matrix_copy_transpose (const gretl_matrix *m);

gretl_matrix *gretl_matrix_reverse_rows (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_reverse_cols (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_get_diagonal (const gretl_matrix *m, int *err);

int gretl_matrix_set_diagonal (gretl_matrix *targ,
			       const gretl_matrix *src,
			       double x);

gretl_matrix *gretl_matrix_get_triangle (const gretl_matrix *m,
					 int upper, int *err);

int gretl_matrix_set_triangle (gretl_matrix *targ,
			       const gretl_matrix *src,
			       double x, int upper);

int gretl_matrix_get_row (const gretl_matrix *m, int i, gretl_vector *v);

double gretl_matrix_trace (const gretl_matrix *m);

int gretl_matrix_random_fill (gretl_matrix *m, int dist);

int correlated_normal_vec (gretl_vector *X, const gretl_matrix *L,
                           int prec);

gretl_matrix *gretl_random_matrix_new (int r, int c, int dist);

gretl_matrix *gretl_matrix_resample (const gretl_matrix *m,
				     int draws, int *err);

int gretl_matrix_resample2 (gretl_matrix *targ,
			    const gretl_matrix *src);

gretl_matrix *gretl_matrix_block_resample (const gretl_matrix *m,
					   int blocklen, int draws,
					   int *err);

int gretl_matrix_block_resample2 (gretl_matrix *targ,
				  const gretl_matrix *src,
				  int blocklen, int *z);

double gretl_vector_mean (const gretl_vector *v);

double gretl_vector_variance (const gretl_vector *v);

void gretl_matrix_zero (gretl_matrix *m);

int gretl_matrix_zero_upper (gretl_matrix *m);

int gretl_matrix_zero_lower (gretl_matrix *m);

int gretl_matrix_mirror (gretl_matrix *m, char uplo);

void gretl_matrix_fill (gretl_matrix *m, double x);

void gretl_matrix_multiply_by_scalar (gretl_matrix *m, double x);

int gretl_matrix_divide_by_scalar (gretl_matrix *m, double x);

void gretl_matrix_switch_sign (gretl_matrix *m);

gretl_matrix *gretl_matrix_dot_op (const gretl_matrix *a,
				   const gretl_matrix *b,
				   int op, int *err);

ConfType dot_operator_conf (const gretl_matrix *A,
			    const gretl_matrix *B,
			    int *r, int *c);

gretl_matrix *gretl_matrix_complex_multiply (const gretl_matrix *a,
					     const gretl_matrix *b,
					     int force_complex,
					     int *err);

gretl_matrix *gretl_matrix_divide (const gretl_matrix *a,
				   const gretl_matrix *b,
				   GretlMatrixMod mod,
				   int *err);

gretl_matrix *gretl_matrix_complex_divide (const gretl_matrix *a,
					   const gretl_matrix *b,
					   int force_complex,
					   int *err);

gretl_matrix *gretl_matrix_exp (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_polroots (const gretl_matrix *a,
				     int force_complex,
				     int legacy,
				     int *err);

void gretl_matrix_raise (gretl_matrix *m, double x);

void gretl_matrix_free (gretl_matrix *m);

double *gretl_matrix_steal_data (gretl_matrix *m);

int gretl_vector_copy_values (gretl_vector *targ,
			      const gretl_vector *src);

int gretl_matrix_copy_values (gretl_matrix *targ,
			      const gretl_matrix *src);

int gretl_matrix_copy_data (gretl_matrix *targ,
			    const gretl_matrix *src);

int gretl_matrix_copy_values_shaped (gretl_matrix *targ,
				     const gretl_matrix *src);

int gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src);

int gretl_matrix_add (const gretl_matrix *a, const gretl_matrix *b,
		      gretl_matrix *c);

int gretl_matrix_add_transpose_to (gretl_matrix *targ,
				   const gretl_matrix *src);

int
gretl_matrix_subtract_from (gretl_matrix *targ, const gretl_matrix *src);

int
gretl_matrix_subtract (const gretl_matrix *a, const gretl_matrix *b,
		       gretl_matrix *c);

int
gretl_matrix_subtract_reversed (const gretl_matrix *a, gretl_matrix *b);

int gretl_matrix_I_minus (gretl_matrix *m);

int gretl_matrix_transpose_in_place (gretl_matrix *m);

int gretl_matrix_transpose (gretl_matrix *targ, const gretl_matrix *src);

int gretl_square_matrix_transpose (gretl_matrix *m);

int gretl_matrix_add_self_transpose (gretl_matrix *m);

int
gretl_matrix_vectorize (gretl_matrix *targ, const gretl_matrix *src);

gretl_matrix *gretl_matrix_vectorize_new (const gretl_matrix *m);

int gretl_matrix_unvectorize (gretl_matrix *targ, const gretl_matrix *src);

int
gretl_matrix_vectorize_h (gretl_matrix *targ, const gretl_matrix *src);

int
gretl_matrix_vectorize_h_skip (gretl_matrix *targ, const gretl_matrix *src);

int
gretl_matrix_unvectorize_h (gretl_matrix *targ, const gretl_matrix *src);

int
gretl_matrix_unvectorize_h_diag (gretl_matrix *targ,
				 const gretl_matrix *src,
				 double diag);

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

int gretl_matrix_multiply_mod_single (const gretl_matrix *a,
				      GretlMatrixMod amod,
				      const gretl_matrix *b,
				      GretlMatrixMod bmod,
				      gretl_matrix *c,
				      GretlMatrixMod cmod);

int gretl_matrix_multiply (const gretl_matrix *a,
			   const gretl_matrix *b,
			   gretl_matrix *c);

int gretl_matrix_multiply_single (const gretl_matrix *a,
				  const gretl_matrix *b,
				  gretl_matrix *c);

gretl_matrix *gretl_matrix_multiply_new (const gretl_matrix *a,
					 const gretl_matrix *b,
					 int *err);

void gretl_blas_dgemm (const gretl_matrix *a, int atr,
		       const gretl_matrix *b, int btr,
		       gretl_matrix *c, GretlMatrixMod cmod,
		       int m, int n, int k);

void gretl_blas_dsymm (const gretl_matrix *a, int asecond,
		       const gretl_matrix *b, int upper,
		       gretl_matrix *c, GretlMatrixMod cmod,
		       int m, int n);

void gretl_blas_dtrmm (const gretl_matrix *a,
                       gretl_matrix *b,
                       const char *flags);

int
gretl_matrix_kronecker_product (const gretl_matrix *A,
				const gretl_matrix *B,
				gretl_matrix *K);

gretl_matrix *
gretl_matrix_kronecker_product_new (const gretl_matrix *A,
				    const gretl_matrix *B,
				    int *err);

int gretl_matrix_hdproduct (const gretl_matrix *A,
			    const gretl_matrix *B,
			    gretl_matrix *C);

gretl_matrix *
gretl_matrix_hdproduct_new (const gretl_matrix *A,
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
				double s, int *err);

double gretl_matrix_dot_product (const gretl_matrix *a, GretlMatrixMod amod,
				 const gretl_matrix *b, GretlMatrixMod bmod,
				 int *err);

double gretl_vector_dot_product (const gretl_vector *a, const gretl_vector *b,
				 int *err);

gretl_matrix *gretl_rmatrix_vector_stat (const gretl_matrix *m,
					 GretlVecStat vs, int rowwise,
					 int skip_na, int *err);

gretl_matrix *gretl_matrix_column_sd (const gretl_matrix *m,
				      int df, int skip_na,
				      int *err);

void gretl_matrix_demean_by_row (gretl_matrix *m);

int gretl_matrix_standardize (gretl_matrix *m, int dfcorr, int skip_na);

int gretl_matrix_center (gretl_matrix *m, int skip_na);

gretl_matrix *gretl_matrix_quantiles (const gretl_matrix *m,
				      const gretl_matrix *p,
				      int *err);

double gretl_matrix_determinant (gretl_matrix *a, int *err);

double gretl_matrix_log_determinant (gretl_matrix *a, int *err);

double gretl_matrix_log_abs_determinant (gretl_matrix *a, int *err);

double gretl_vcv_log_determinant (const gretl_matrix *m, int *err);

double gretl_matrix_one_norm (const gretl_matrix *m);

double gretl_matrix_infinity_norm (const gretl_matrix *m);

int gretl_LU_solve (gretl_matrix *a, gretl_matrix *b);

int gretl_LU_solve_invert (gretl_matrix *a, gretl_matrix *b);

int gretl_cholesky_decomp_solve (gretl_matrix *a, gretl_matrix *b);

int gretl_cholesky_solve (const gretl_matrix *a, gretl_vector *b);

int gretl_cholesky_invert (gretl_matrix *a);

int cholesky_factor_of_inverse (gretl_matrix *a);

gretl_vector *gretl_toeplitz_solve (const gretl_vector *c,
				    const gretl_vector *r,
				    const gretl_vector *b,
				    double *det,
				    int *err);

gretl_matrix *gretl_matrix_XTX_new (const gretl_matrix *X);

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

int gretl_matrix_moore_penrose (gretl_matrix *A, double tol);

int gretl_SVD_invert_matrix (gretl_matrix *a);

int gretl_invpd (gretl_matrix *a);

int gretl_matrix_SVD (const gretl_matrix *a, gretl_matrix **pu,
		      gretl_vector **ps, gretl_matrix **pvt,
		      int full);

double gretl_symmetric_matrix_rcond (const gretl_matrix *m, int *err);

double gretl_matrix_rcond (const gretl_matrix *m, int *err);

double gretl_matrix_cond_index (const gretl_matrix *m, int *err);

int gretl_symmetric_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs,
				int rank);

gretl_matrix *
gretl_general_matrix_eigenvals (const gretl_matrix *m, int *err);

gretl_matrix *
gretl_symmetric_matrix_eigenvals (gretl_matrix *m,
				  int eigenvecs,
				  int *err);

double gretl_symmetric_matrix_min_eigenvalue (const gretl_matrix *m);

gretl_matrix *
gretl_symm_matrix_eigenvals_descending (gretl_matrix *m,
					int eigenvecs,
					int *err);

gretl_matrix *
gretl_gensymm_eigenvals (const gretl_matrix *A,
			 const gretl_matrix *B,
			 gretl_matrix *V,
			 int *err);

gretl_matrix *gretl_dgeev (const gretl_matrix *A,
			   gretl_matrix *VR,
			   gretl_matrix *VL,
			   int *err);

gretl_matrix *old_eigengen (const gretl_matrix *m,
                            gretl_matrix *VR,
                            gretl_matrix *VL,
                            int *err);

double gretl_symm_matrix_lambda_min (const gretl_matrix *m, int *err);

double gretl_symm_matrix_lambda_max (const gretl_matrix *m, int *err);

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

gretl_matrix *gretl_matrix_direct_sum (const gretl_matrix *a,
				       const gretl_matrix *b,
				       int *err);

int
gretl_matrix_inplace_colcat (gretl_matrix *a, const gretl_matrix *b,
			     const char *mask);

gretl_matrix *gretl_matrix_cumcol (const gretl_matrix *m, int *err);

gretl_matrix *gretl_matrix_diffcol (const gretl_matrix *m,
				    double missval, int *err);

gretl_matrix *gretl_matrix_lag (const gretl_matrix *m,
				const gretl_vector *k,
				gretlopt opt,
				double missval);

int gretl_matrix_inplace_lag (gretl_matrix *targ,
			      const gretl_matrix *src,
			      int k);

int gretl_matrix_cholesky_decomp (gretl_matrix *a);

int gretl_matrix_psd_root (gretl_matrix *a, int check);

int gretl_matrix_QR_decomp (gretl_matrix *M,
			    gretl_matrix *R);

int gretl_matrix_QR_pivot_decomp (gretl_matrix *M,
				  gretl_matrix *R,
				  gretl_matrix *P);

double gretl_triangular_matrix_rcond (const gretl_matrix *A,
				      char uplo, char diag);

int gretl_check_QR_rank (const gretl_matrix *R,
			 int *err,
			 double *rcnd);

int gretl_matrix_rank (const gretl_matrix *a, double eps,
		       int *err);

int gretl_matrix_ols (const gretl_vector *y,
		      const gretl_matrix *X,
		      gretl_vector *b,
		      gretl_matrix *vcv,
		      gretl_vector *uhat,
		      double *s2);

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
gretl_matrix_restricted_ols (const gretl_vector *y,
			     const gretl_matrix *X,
			     const gretl_matrix *R,
			     const gretl_vector *q,
			     gretl_vector *b,
			     gretl_matrix *vcv,
			     gretl_vector *uhat,
			     double *s2);

int
gretl_matrix_restricted_multi_ols (const gretl_matrix *Y,
				   const gretl_matrix *X,
				   const gretl_matrix *R,
				   const gretl_matrix *q,
				   gretl_matrix *B,
				   gretl_matrix *U,
				   gretl_matrix **W);

int gretl_matrix_SVD_ols (const gretl_vector *y,
			  const gretl_matrix *X,
			  gretl_vector *b,
			  gretl_matrix *vcv,
			  gretl_vector *uhat,
			  double *s2);

int gretl_matrix_qform (const gretl_matrix *A,
			GretlMatrixMod amod,
			const gretl_matrix *X,
			gretl_matrix *C,
			GretlMatrixMod cmod);

int gretl_matrix_diag_qform (const gretl_matrix *A,
			     GretlMatrixMod amod,
			     const gretl_vector *X,
			     gretl_matrix *C,
			     GretlMatrixMod cmod);

double gretl_scalar_qform (const gretl_vector *b,
			   const gretl_matrix *X,
			   int *err);

int gretl_matrix_columnwise_product (const gretl_matrix *A,
				     const gretl_matrix *B,
				     const gretl_matrix *S,
				     gretl_matrix *C);

int
gretl_matrix_diagonal_sandwich (const gretl_vector *d,
				const gretl_matrix *X,
				gretl_matrix *DXD);

int gretl_matrix_set_t1 (gretl_matrix *m, int t);

int gretl_matrix_set_t2 (gretl_matrix *m, int t);

int gretl_matrix_get_t1 (const gretl_matrix *m);

int gretl_matrix_get_t2 (const gretl_matrix *m);

int gretl_matrix_is_dated (const gretl_matrix *m);

int gretl_is_identity_matrix (const gretl_matrix *m);

int gretl_is_zero_matrix (const gretl_matrix *m);

gretl_matrix *gretl_matrix_isfinite (const gretl_matrix *m, int *err);

int gretl_matrix_get_structure (const gretl_matrix *m);

int gretl_matrices_are_equal (const gretl_matrix *a,
			      const gretl_matrix *b,
			      double tol, int *err);

gretl_matrix *gretl_covariance_matrix (const gretl_matrix *m,
				       int corr, int dfc,
				       int *err);

gretl_matrix *gretl_matrix_GG_inverse (const gretl_matrix *G,
				       int *err);

gretl_matrix *gretl_matrix_varsimul (const gretl_matrix *A,
				     const gretl_matrix *U,
				     const gretl_matrix *x0,
				     int *err);

gretl_matrix **gretl_matrix_array_new (int n);

gretl_matrix **
gretl_matrix_array_new_with_size (int n, int rows, int cols);

void gretl_matrix_array_free (gretl_matrix **A, int n);

gretl_matrix *gretl_matrix_values (const double *x, int n,
				   gretlopt opt, int *err);

gretl_matrix *gretl_matrix_values_full (const double *x, int n,
					gretlopt opt, int *missvals,
					int *err);

int gretl_matrix_n_values (const double *x, int n, int *err);

gretl_matrix *gretl_matrix_shape (const gretl_matrix *A,
				  int r, int c, int *err);

gretl_matrix *gretl_matrix_trim_rows (const gretl_matrix *A,
				      int ttop, int tbot,
				      int *err);

gretl_matrix *gretl_matrix_minmax (const gretl_matrix *A,
				   int mm, int rc, int idx,
				   int skip_na, int *err);

double gretl_matrix_global_minmax (const gretl_matrix *A,
				   int mm, int *err);

double gretl_matrix_global_sum (const gretl_matrix *A,
				int *err);

gretl_matrix *gretl_matrix_pca (const gretl_matrix *X, int p,
				gretlopt opt, int *err);

gretl_matrix *gretl_matrix_xtab (const double *x,
				 const double *y,
				 int n, int *err);

gretl_matrix *gretl_matrix_bool_sel(const gretl_matrix *A,
				    const gretl_matrix *sel,
				    int rowsel, int *err);

gretl_matrix *gretl_matrix_sort_by_column (const gretl_matrix *m,
					   int k, int *err);

gretl_matrix *gretl_vector_sort (const gretl_matrix *v,
				 int descending,
				 int *err);

gretl_matrix *gretl_matrix_covariogram (const gretl_matrix *X,
					const gretl_matrix *u,
					const gretl_matrix *w,
					int p, int *err);

gretl_matrix *gretl_matrix_commute(gretl_matrix *A, int r, int c,
				   int pre, int add_id, int *err);

void gretl_matrix_transcribe_obs_info (gretl_matrix *targ,
				       const gretl_matrix *src);

int gretl_matrix_set_colnames (gretl_matrix *m, char **S);

int gretl_matrix_set_rownames (gretl_matrix *m, char **S);

const char **gretl_matrix_get_colnames (const gretl_matrix *m);

const char **gretl_matrix_get_rownames (const gretl_matrix *m);

void gretl_matrix_destroy_info (gretl_matrix *m);

void lapack_mem_free (void);

void set_blas_mnk_min (int mnk);

int get_blas_mnk_min (void);

void set_simd_k_max (int k);

int get_simd_k_max (void);

void set_simd_mn_min (int mn);

int get_simd_mn_min (void);

#ifdef  __cplusplus
}
#endif

#endif /* GRETL_MATRIX_H */
