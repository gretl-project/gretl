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

#ifndef USERMAT_H_
#define USERMAT_H_

#define MSEL_MAX -999

enum {
    SEL_RANGE,
    SEL_ELEMENT,
    SEL_MATRIX,
    SEL_DIAG,
    SEL_ALL,
    SEL_NULL
};

typedef struct matrix_subspec_ matrix_subspec;

union msel {
    int range[2];
    gretl_matrix *m;
};

struct matrix_subspec_ {
    int type[2];
    union msel sel[2];
    int *rslice;
    int *cslice;
};

#define mspec_get_row_index(m) (m->sel[0].range[0])
#define mspec_get_col_index(m) (m->sel[1].range[0])

#define mspec_set_row_index(m,i) (m->sel[0].range[0] = m->sel[0].range[1] = (i))
#define mspec_set_col_index(m,j) (m->sel[1].range[0] = m->sel[1].range[1] = (j))

#define gretl_is_matrix(s) (get_matrix_by_name(s) != NULL)

GList *get_named_matrix_list (void);

gretl_matrix *get_matrix_by_name (const char *name);

gretl_matrix *get_matrix_copy_by_name (const char *name, int *err);

gretl_matrix *steal_matrix_by_name (const char *name);

int assign_scalar_to_submatrix (gretl_matrix *M, double x,
				matrix_subspec *spec);

int user_matrix_replace_submatrix (const char *mname, 
				   const gretl_matrix *S,
				   matrix_subspec *spec);

int umatrix_set_names_from_string (gretl_matrix *M, 
				   const char *s,
				   int byrow);

int umatrix_set_names_from_list (gretl_matrix *M, 
				 const int *list,
				 const DATASET *dset,
				 int byrow);

char *user_matrix_get_column_name (const gretl_matrix *M, int col,
				   int *err);

double user_matrix_get_determinant (gretl_matrix *m, int tmpmat,
				    int f, int *err);

gretl_matrix *user_matrix_matrix_func (gretl_matrix *m, int tmpmat,
				       int f, int *err);

gretl_matrix *user_matrix_vec (const gretl_matrix *m, int *err);

gretl_matrix *user_matrix_vech (const gretl_matrix *m, int *err);

gretl_matrix *user_matrix_unvech (const gretl_matrix *m, int *err);

gretl_matrix *
user_matrix_QR_decomp (const gretl_matrix *m, const char *rname, 
		       int *err);

gretl_matrix *user_matrix_SVD (const gretl_matrix *m, 
			       const char *uname, 
			       const char *vname, 
			       int *err);

gretl_matrix *user_matrix_ols (const gretl_matrix *Y, 
			       const gretl_matrix *X, 
			       const char *Uname, 
			       const char *Vname,
			       gretlopt opt,
			       int *err);

gretl_matrix *user_matrix_rls (const gretl_matrix *Y, 
			       const gretl_matrix *X,
			       const gretl_matrix *R,
			       const gretl_matrix *Q,
			       const char *Uname, 
			       const char *Vname, 
			       int *err);

gretl_matrix *
user_matrix_eigen_analysis (const gretl_matrix *m, const char *rname, int symm,
			    int *err);

gretl_matrix *user_gensymm_eigenvals (const gretl_matrix *A, 
				      const gretl_matrix *B,
				      const char *rname,
				      int *err);

double matrix_get_element (const gretl_matrix *M, int i, int j,
			   int *err);

int check_matrix_subspec (matrix_subspec *spec, const gretl_matrix *m);

gretl_matrix *matrix_get_submatrix (const gretl_matrix *M, 
				    matrix_subspec *spec,
				    int prechecked,
				    int *err);

gretl_matrix *user_matrix_get_submatrix (const char *name, 
					 matrix_subspec *spec,
					 int *err);

int matrix_invert_in_place (gretl_matrix *m);

int matrix_cholesky_in_place (gretl_matrix *m);

int matrix_transpose_in_place (gretl_matrix *m);

int matrix_XTX_in_place (gretl_matrix *m);

#endif /* USERMAT_H_ */
