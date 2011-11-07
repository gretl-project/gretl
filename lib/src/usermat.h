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

typedef struct user_matrix_ user_matrix;

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

int n_user_matrices (void);

const char *get_matrix_name_by_index (int idx);

gretl_matrix *get_matrix_by_name (const char *name);

gretl_matrix *get_matrix_copy_by_name (const char *name, int *err);

gretl_matrix *steal_matrix_by_name (const char *name);

gretl_matrix *get_matrix_by_name_at_level (const char *name, int level);

gretl_matrix *user_matrix_get_matrix (user_matrix *u);

user_matrix *get_user_matrix_by_name (const char *name);

user_matrix *get_user_matrix_by_index (int idx);

user_matrix *get_user_matrix_by_data (const gretl_matrix *M);

int user_matrix_add (gretl_matrix *M, const char *name);

int private_matrix_add (gretl_matrix *M, const char *name);

int matrix_add_as_shell (gretl_matrix *M, const char *name);

int user_matrix_destroy_by_name (const char *name, PRN *prn);

int user_matrix_destroy (user_matrix *u);

int user_matrix_get_level (user_matrix *u);

int user_matrix_adjust_level (user_matrix *u, int adj);

const char *user_matrix_get_name (user_matrix *u);

int user_matrix_set_name (user_matrix *u, const char *name);

int user_matrix_replace_matrix (user_matrix *u, gretl_matrix *M);

int user_matrix_replace_matrix_by_name (const char *name, 
					gretl_matrix *M);

int user_matrix_replace_submatrix (const char *mname, 
				   const gretl_matrix *S,
				   matrix_subspec *spec);

int add_or_replace_user_matrix (gretl_matrix *M, const char *name);

int copy_named_matrix_as (const char *orig, const char *newname);

int copy_matrix_as (const gretl_matrix *m, const char *newname);

int umatrix_set_names_from_string (gretl_matrix *M, 
				   const char *s,
				   int byrow);

int umatrix_set_names_from_list (gretl_matrix *M, 
				 const int *list,
				 const DATASET *dset,
				 int byrow);

const char **user_matrix_get_names (const gretl_matrix *M,
				    int byrow);

char *user_matrix_get_column_name (const gretl_matrix *M, int col,
				   int *err);

void destroy_user_matrices (void);

int destroy_user_matrices_at_level (int level);

int destroy_private_matrices (void);

double user_matrix_get_determinant (gretl_matrix *m, int f, int *err);

gretl_matrix *user_matrix_matrix_func (gretl_matrix *m, int f, int *err);

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

int user_eigen_sort (gretl_matrix *evals, gretl_matrix *evecs,
		     int rank);

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

void write_matrices_to_file (FILE *fp);

void set_matrix_add_callback (void (*callback));

#endif /* USERMAT_H_ */
