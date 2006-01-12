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

#ifndef USERMAT_H_
#define USERMAT_H_

gretl_matrix *get_matrix_by_name (const char *name);

gretl_matrix *get_matrix_by_name_at_level (const char *name, int level);

int named_matrix_get_scalar (const char *name, double *x);

int add_or_replace_user_matrix (gretl_matrix *M, const char *name,
				const char *mask, gretl_matrix **R);

int copy_named_matrix_as (const char *orig, const char *new);

int user_matrix_reconfigure (gretl_matrix *M, char *newname, int level);

void destroy_user_matrices (void);

int destroy_user_matrices_at_level (int level);

int is_user_matrix (gretl_matrix *m);

gretl_matrix *user_matrix_get_slice (const char *s, int *err);

int matrix_command (const char *line, double ***pZ, DATAINFO *pdinfo, PRN *prn);

gretl_matrix *matrix_calc_AB (gretl_matrix *A, gretl_matrix *B, 
			      char op, int *err);

double user_matrix_get_determinant (gretl_matrix *m);

double user_matrix_get_log_determinant (gretl_matrix *m);

gretl_matrix *user_matrix_get_determinant_as_matrix (gretl_matrix *m);

gretl_matrix *user_matrix_get_log_determinant_as_matrix (gretl_matrix *m);

gretl_matrix *user_matrix_get_inverse (gretl_matrix *m);

gretl_matrix *
user_matrix_get_transformation (gretl_matrix *m, GretlMathFunc fn);

int reposition_transpose_symbol (char *s);

#endif /* USERMAT_H_ */
