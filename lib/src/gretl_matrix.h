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
    GRETL_MATRIX_RANGE
};

enum gretl_matrix_mods {
    GRETL_MOD_NONE = 0,
    GRETL_MOD_TRANSPOSE = 1
};

typedef struct _gretl_matrix gretl_matrix;
typedef struct _gretl_matrix gretl_vector;

struct _gretl_matrix {
    int packed;
    int rows;
    int cols;
    double *val;
};

#define mdx(a,i,j)   (j)*(a)->rows+i
#define mdxtr(a,i,j) (i)*(a)->rows+j

#define gretl_vector_alloc(i) gretl_matrix_alloc(1,(i))
#define gretl_vector_free(v) gretl_matrix_free(v)
#define gretl_vector_get(v,i) gretl_matrix_get((v),0,(i))
#define gretl_vector_set(v,i,x) gretl_matrix_set((v),0,(i),(x))
#define gretl_vector_get_length(v) (v)->cols


gretl_matrix *gretl_matrix_alloc (int rows, int cols);

gretl_matrix *gretl_packed_matrix_alloc (int rows);

gretl_matrix *gretl_matrix_copy (gretl_matrix *m);

gretl_matrix *gretl_matrix_from_2d_array (const double **X, 
					  int rows, int cols);

void gretl_matrix_free (gretl_matrix *m);

double *gretl_matrix_steal_data (gretl_matrix *m);

double gretl_matrix_get (const gretl_matrix *m, int i, int j);

int gretl_matrix_set (gretl_matrix *m, int i, int j, double x);

int gretl_matmult_mod (const gretl_matrix *a, int aflag,
		       const gretl_matrix *b, int bflag,
		       gretl_matrix *c);

int gretl_matmult (const gretl_matrix *a, const gretl_matrix *b,
		   gretl_matrix *c);

int gretl_LU_solve (gretl_matrix *a, gretl_vector *b);

int gretl_invert_general_matrix (gretl_matrix *m);

double *gretl_general_matrix_eigenvals (gretl_matrix *m);

double *gretl_symmetric_matrix_eigenvals (gretl_matrix *m,
					  int eigenvecs);

#endif /* GRETL_MATRIX_H */
