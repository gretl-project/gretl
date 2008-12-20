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
gretl_covariance_matrix_from_varlist (const int *list, const double **Z, 
				      const DATAINFO *pdinfo, 
				      gretl_matrix **means,
				      int *errp);

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x);


gretl_matrix *gretl_matrix_data_subset (const int *list, const double **Z,
					int t1, int t2, const char *mask,
					int *err);

gretl_matrix *
gretl_matrix_data_subset_no_missing (const int *list, const double **Z,
				     int t1, int t2, int *err);

gretl_matrix *
gretl_matrix_data_subset_skip_missing (const int *list, const double **Z,
				       int t1, int t2, int *err);

int gretl_plotfit_matrices (int yno, int xno, FitType fit,
			    const double **Z, int t1, int t2, 
			    gretl_matrix **py, gretl_matrix **pX);

int gretl_matrix_delete_columns (gretl_matrix *X, int *list);

gretl_matrix *gretl_matrix_read_from_text (const char *fname, int *err);

int gretl_matrix_write_as_text (gretl_matrix *A, const char *fname);

void 
gretl_matrix_print_to_prn (const gretl_matrix *m, const char *msg, PRN *prn);

void gretl_packed_matrix_print (const gretl_matrix *m, const char *msg);

void debug_print_matrix (const gretl_matrix *m, const char *msg);

void gretl_matrix_print_with_col_heads (const gretl_matrix *m, 
					const char *title,
					const char **heads,
					PRN *prn);

void gretl_matrix_print_with_format (const gretl_matrix *m, 
				     const char *fmt,
				     int wid, int prec,
				     PRN *prn);

int gretl_matrix_cut_rows (gretl_matrix *m, const char *mask);

int gretl_matrix_cut_rows_cols (gretl_matrix *m, const char *mask);

char *gretl_matrix_zero_row_mask (const gretl_matrix *m, int *err);

char *gretl_matrix_rank_mask (const gretl_matrix *m, int *err);

int gretl_matrix_mp_ols (const gretl_vector *y, const gretl_matrix *X,
			 gretl_vector *b, gretl_matrix *vcv, 
			 gretl_vector *uhat, double *s2);

#endif /* MATRIX_EXTRA_H */
