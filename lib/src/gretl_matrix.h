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

#ifndef GRETL_MATRIX_H
#define GRETL_MATRIX_H

#include "f2c.h"
#include "clapack_double.h"

/* #define LDEBUG 1 */

enum gretl_matrix_errors {
    GRETL_MATRIX_OK = 0,
    GRETL_MATRIX_NOMEM,
    GRETL_MATRIX_NON_CONFORM,
    GRETL_MATRIX_RANGE,
    GRETL_MATRIX_ERR 
};

enum gretl_matrix_mods {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE,
    GRETL_MOD_SQUARE
};

typedef struct _gretl_matrix gretl_matrix;
typedef struct _gretl_matrix gretl_vector;

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))
#define gretl_column_vector_alloc(i) gretl_matrix_alloc((i),1)
#define gretl_vector_free(v) gretl_matrix_free(v)

gretl_matrix *gretl_matrix_alloc (int rows, int cols);

gretl_matrix *gretl_packed_matrix_alloc (int rows);

gretl_matrix *gretl_matrix_copy (const gretl_matrix *m);

gretl_matrix *gretl_diagonal_matrix (const double *d, int n, int mod);

gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
					  int rows, int cols);

void gretl_matrix_zero (gretl_matrix *m);

int gretl_matrix_zero_upper (gretl_matrix *m);

int gretl_matrix_zero_lower (gretl_matrix *m);

void gretl_matrix_multiply_by_scalar (gretl_matrix *m, double x);

void gretl_matrix_divide_by_scalar (gretl_matrix *m, double x);

void gretl_matrix_free (gretl_matrix *m);

double *gretl_matrix_steal_data (gretl_matrix *m);

int gretl_matrix_copy_values (gretl_matrix *targ, 
			      const gretl_matrix *src);

double gretl_matrix_get (const gretl_matrix *m, int i, int j);

double gretl_vector_get (const gretl_vector *v, int i);

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x);

int gretl_vector_set (gretl_vector *v, int i, double x);

int gretl_matrix_add_to (gretl_matrix *targ, const gretl_matrix *src);

int gretl_square_matrix_transpose (gretl_matrix *m);

int gretl_matrix_add_self_transpose (gretl_matrix *m);

int gretl_matrix_multiply_mod (const gretl_matrix *a, int aflag,
			       const gretl_matrix *b, int bflag,
			       gretl_matrix *c);

int gretl_matrix_multiply (const gretl_matrix *a, const gretl_matrix *b,
			   gretl_matrix *c);

double gretl_matrix_dot_product (const gretl_matrix *a, int aflag,
				 const gretl_matrix *b, int bflag,
				 int *err);

gretl_matrix *gretl_matrix_vcv (gretl_matrix *m);

double gretl_LU_determinant (gretl_matrix *a);

int gretl_LU_solve (gretl_matrix *a, gretl_vector *b);

int gretl_invert_general_matrix (gretl_matrix *a);

int gretl_invert_symmetric_matrix (gretl_matrix *a);

double *gretl_general_matrix_eigenvals (gretl_matrix *m, gretl_matrix *ev);

double *gretl_symmetric_matrix_eigenvals (gretl_matrix *m,
					  int eigenvecs);

int gretl_matrix_cholesky_decomp (gretl_matrix *a);

int gretl_matrix_ols (const gretl_vector *y, const gretl_matrix *X,
		      gretl_vector *b, gretl_matrix *vcv);

double gretl_scalar_b_prime_X_b (const gretl_vector *b, const gretl_matrix *X,
				 int *err);

gretl_matrix *
gretl_vcv_matrix_from_model (const MODEL *pmod, const char *select);

void gretl_matrix_print (gretl_matrix *m, const char *msg, PRN *prn);

void gretl_matrix_set_int (gretl_matrix *m, int t);

int gretl_matrix_get_int (const gretl_matrix *m);

int gretl_vector_get_length (const gretl_vector *v);

int gretl_matrix_cols (const gretl_matrix *m);

int gretl_matrix_rows (const gretl_matrix *m);

#endif /* GRETL_MATRIX_H */
