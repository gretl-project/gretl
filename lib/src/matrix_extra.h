/*
 *  Copyright (c) by Allin Cottrell and Riccardo "Jack" Lucchetti
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
gretl_vcv_matrix_from_model (MODEL *pmod, const char *select);

gretl_vector *
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select);

gretl_matrix *
gretl_covariance_matrix_from_varlist (const int *list, const double **Z, 
				      const DATAINFO *pdinfo, 
				      gretl_matrix **means,
				      int *errp);

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x);


gretl_matrix *gretl_matrix_data_subset (const int *list, const double **Z,
					int t1, int t2, const char *mask);

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

#endif /* MATRIX_EXTRA_H */
