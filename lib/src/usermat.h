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

#define MSEL_MAX INT_MIN

typedef enum {
    SEL_NULL,    /* nothing supplied */
    SEL_RANGE,   /* integer range p:q provided */
    SEL_MATRIX,  /* selection matrix provided */
    SEL_ALL,     /* comma-separated blank */
    SEL_DIAG,    /* the "diag" dummy constant */
    SEL_UPPER,   /* the "upper" dummy constant */
    SEL_LOWER,   /* the "lower" dummy constant */
    SEL_REAL,    /* the "real" dummy constant */
    SEL_IMAG,    /* the "imag" dummy constant */
    SEL_ELEMENT, /* derived: selection is a single element */
    SEL_CONTIG,  /* derived: selection is contiguous */
    SEL_EXCL,    /* single exclusion (negative index) */
    SEL_SINGLE,  /* derived: degenerate range + null */
    SEL_STR      /* for use with bundles only */
} SelType;

#define is_sel_dummy(s) (s >= SEL_DIAG && s <= SEL_IMAG)

/* Note SEL_EXCL is flagged only in the case of a single negative
   index. SEL_MATRIX can also do exclusion, if all the elements
   of the vector are negative.
*/

typedef struct matrix_subspec_ matrix_subspec;

union msel {
    int range[2];
    gretl_matrix *m;
    char *str;
};

struct matrix_subspec_ {
    int checked;
    SelType ltype, rtype;
    union msel lsel, rsel;
    int *rslice;
    int *cslice;
};

#define mspec_get_row_index(m) (m->lsel.range[0])
#define mspec_get_col_index(m) (m->rsel.range[0])

#define mspec_set_row_index(m,i) (m->lsel.range[0] = m->lsel.range[1] = (i))
#define mspec_set_col_index(m,j) (m->rsel.range[0] = m->rsel.range[1] = (j))

#define mspec_get_element(m) (m->lsel.range[0])

matrix_subspec *matrix_subspec_new (void);

GList *get_named_matrix_list (void);

gretl_matrix *get_matrix_by_name (const char *name);

gretl_matrix *get_matrix_copy_by_name (const char *name, int *err);

gretl_matrix *steal_matrix_by_name (const char *name);

int assign_scalar_to_submatrix (gretl_matrix *M,
				const gretl_matrix *S,
				double x,
				matrix_subspec *spec);

int matrix_replace_submatrix (gretl_matrix *M,
			      const gretl_matrix *S,
			      matrix_subspec *spec);

int umatrix_set_names_from_string (gretl_matrix *M,
				   const char *s,
				   int byrow,
                                   PRN *prn);

int umatrix_set_names_from_array (gretl_matrix *M,
				  void *data,
				  int byrow,
                                  PRN *prn);

int umatrix_set_names_from_list (gretl_matrix *M,
				 const int *list,
				 const DATASET *dset,
				 int byrow,
                                 PRN *prn);

char *user_matrix_get_column_name (const gretl_matrix *M, int col,
				   int *err);

char *user_matrix_get_row_name (const gretl_matrix *M, int row,
				int *err);

double user_matrix_get_determinant (gretl_matrix *m, int tmpmat,
				    int ldet, int *err);

gretl_matrix *user_matrix_vec (const gretl_matrix *m, int *err);

gretl_matrix *user_matrix_vech (const gretl_matrix *m,
				int omit_diag, int *err);

gretl_matrix *user_matrix_unvech (const gretl_matrix *m,
				  double d, int *err);

gretl_matrix *user_matrix_QR_decomp (const gretl_matrix *m,
				     gretl_matrix *R,
				     gretl_matrix *P,
				     int *err);

gretl_matrix *user_matrix_SVD (const gretl_matrix *m,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       int *err);

gretl_matrix *user_matrix_ols (const gretl_matrix *Y,
			       const gretl_matrix *X,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       gretlopt opt,
			       int *err);

gretl_matrix *user_matrix_rls (const gretl_matrix *Y,
			       const gretl_matrix *X,
			       const gretl_matrix *R,
			       const gretl_matrix *Q,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       int *err);

gretl_matrix *user_matrix_GHK (const gretl_matrix *C,
			       const gretl_matrix *A,
			       const gretl_matrix *B,
			       const gretl_matrix *U,
			       gretl_matrix *dP,
			       int *err);

gretl_matrix *user_matrix_eigensym (const gretl_matrix *m,
				    gretl_matrix *R,
				    int *err);

gretl_matrix *user_gensymm_eigenvals (const gretl_matrix *A,
				      const gretl_matrix *B,
				      gretl_matrix *V,
				      int *err);

double matrix_get_element (const gretl_matrix *M, int i, int *err);

gretl_matrix *matrix_get_chunk (const gretl_matrix *M,
				matrix_subspec *spec,
				int *err);

int *mspec_make_list (int type, union msel *sel, int n,
		      int *err);

int check_matrix_subspec (matrix_subspec *spec, const gretl_matrix *m);

const char *mspec_get_string (matrix_subspec *spec, int i);

gretl_matrix *matrix_get_submatrix (const gretl_matrix *M,
				    matrix_subspec *spec,
				    int prechecked,
				    int *err);

int matrix_invert_in_place (gretl_matrix *m);

int matrix_cholesky_in_place (gretl_matrix *m);

int matrix_transpose_in_place (gretl_matrix *m);

int matrix_XTX_in_place (gretl_matrix *m);

int gretl_matrix_set_part (gretl_matrix *targ,
			   const gretl_matrix *src,
			   double x, SelType sel);


#endif /* USERMAT_H_ */
