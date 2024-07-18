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

#ifndef MATRIX_EXTRA_H
#define MATRIX_EXTRA_H

typedef enum {
    M_MISSING_OK,
    M_MISSING_ERROR,
    M_MISSING_SKIP,
    M_MISSING_TRIM
} MMissingCode;

typedef enum {
    QUAD_GHERMITE = 1,
    QUAD_LEGENDRE,
    QUAD_LAGUERRE,
    QUAD_INVALID
} QuadMethod;

gretl_vector *
gretl_vector_from_array (const double *x, int n, GretlMatrixMod mod);

gretl_vector *gretl_vector_from_series (const double *x, 
					int t1, int t2);

gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
					  int rows, int cols);

gretl_matrix *gretl_matrix_from_scalar (double x);

gretl_matrix *
gretl_vcv_matrix_from_model (MODEL *pmod, const char *select, int *err);

gretl_vector *
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select, int *err);

gretl_matrix *
gretl_covariance_matrix_from_varlist (const int *list, 
				      const DATASET *dset, 
				      gretl_matrix **means,
				      int *errp);

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x);

double **gretl_matrix_get_columns (const gretl_matrix *m, int *err);

gretl_matrix *
gretl_matrix_data_subset_masked (const int *list, 
				 const DATASET *dset,
				 int t1, int t2, const char *mask, 
				 int *err);

gretl_matrix *gretl_matrix_data_subset (const int *list, 
					const DATASET *dset,
					int t1, int t2, int missop, 
					int *err);

gretl_matrix *
gretl_matrix_data_subset_special (const int *list, 
				  const DATASET *dset,
				  const gretl_matrix *mmask,
				  int *err);

DATASET *gretl_dataset_from_matrix (const gretl_matrix *m, 
				    const int *list,
				    gretlopt opt,
				    int *err);

DATASET *gretl_dataset_from_matrices (const gretl_matrix *m1,
				      const gretl_matrix *m2,
				      int *err);

int write_matrix_as_dataset (const char *fname,
			     gretlopt opt,
			     PRN *prn);

int gretl_plotfit_matrices (const double *yvar, const double *xvar,
			    FitType fit, int t1, int t2, 
			    gretl_matrix **py, gretl_matrix **pX);

gretl_matrix *gretl_matrix_read_from_file (const char *fname, 
					   int import, int *err);

int gretl_matrix_write_to_file (gretl_matrix *A, const char *fname,
				int use_dotdir);

void gretl_matrix_print_to_prn (const gretl_matrix *m,
				const char *msg,
				PRN *prn);

void gretl_matrix_print_range (const gretl_matrix *m,
			       const char *msg,
			       int rmin, int rmax,
			       PRN *prn);

void gretl_matrix_print_with_col_heads (const gretl_matrix *m, 
					const char *title,
					const char **heads,
					const DATASET *dset,
					PRN *prn);

void gretl_matrix_print_with_format (const gretl_matrix *m, 
				     const char *fmt,
				     int wid, int prec,
				     PRN *prn);

gchar *gretl_matrix_write_constructor (const gretl_matrix *m);

void debug_print_matrix (const gretl_matrix *m, const char *msg);

int gretl_matrix_cut_cols (gretl_matrix *m, const char *mask);

int gretl_matrix_cut_rows (gretl_matrix *m, const char *mask);

int gretl_matrix_cut_rows_cols (gretl_matrix *m, const char *mask);

char *gretl_matrix_zero_row_mask (const gretl_matrix *m, int *err);

char *gretl_matrix_zero_col_mask (const gretl_matrix *m, int *err);

char *gretl_matrix_zero_diag_mask (const gretl_matrix *m, int *err);

char *gretl_matrix_rank_mask (const gretl_matrix *m, int *err);

int gretl_matrix_mp_ols (const gretl_vector *y, const gretl_matrix *X,
			 gretl_vector *b, gretl_matrix *vcv, 
			 gretl_vector *uhat, double *s2);

gretl_matrix *gretl_quadrule_matrix_new (int n, int method, 
					 double a, double b,
					 int *err);

gretl_matrix *gretl_gauss_hermite_matrix_new (int n, int *err);

gretl_matrix *gretl_matrix_2d_convolution (const gretl_matrix *A,
					   const gretl_matrix *B,
					   int *err);

gretl_matrix *vector_from_strings (char **S, int ns,
				   const char *fmt,
				   int *nvals,
				   int *err);

#endif /* MATRIX_EXTRA_H */
