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

#include "libgretl.h"
#include "gretl_matrix.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "genparse.h"
#include "uservar.h"

#define MDEBUG 0

/**
 * get_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_by_name (const char *name)
{
    gretl_matrix *ret = NULL;

    if (name != NULL && *name != '\0') {
	user_var *u;

	u = get_user_var_of_type_by_name(name, GRETL_TYPE_MATRIX);
	if (u != NULL) {
	    ret = user_var_get_value(u);
	}
    }

    return ret;
}

/**
 * steal_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name and if found,
 * grabs the matrix, leaving the matrix pointer on the
 * named matrix as %NULL.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *steal_matrix_by_name (const char *name)
{
    gretl_matrix *ret = NULL;

    if (name != NULL && *name != '\0') {
	user_var *u;

	u = get_user_var_of_type_by_name(name, GRETL_TYPE_MATRIX);

	if (u != NULL) {
	    ret = user_var_steal_value(u);
	}
    }

    return ret;
}

/**
 * get_matrix_copy_by_name:
 * @name: name of the matrix.
 * @err: location to receive error code;
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: a copy of the named matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_copy_by_name (const char *name, int *err)
{
    gretl_matrix *m = get_matrix_by_name(name);
    gretl_matrix *ret = NULL;

    if (m == NULL) {
	*err = E_UNKVAR;
    } else {
	ret = gretl_matrix_copy(m);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

static int msel_out_of_bounds (int *range, int n)
{
    int *bad = NULL;

    if (range[0] < 1 || range[0] > n) {
	bad = &range[0];
    } else if (range[1] < 1 || range[1] > n) {
	bad = &range[1];
    }

    if (bad != NULL) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 
			     *bad);
	return 1;
    } else {
	return 0;
    }
}

static int vec_is_const (const gretl_matrix *m, int n)
{
    int i;

    for (i=1; i<n; i++) {
	if (m->val[i] != m->val[i-1]) {
	    return 0;
	}
    }

    return 1;
}

/* convert a matrix subspec component into list of rows 
   or columns */

static int *mspec_make_list (int type, union msel *sel, int n,
			     int *err)
{
    int *slice = NULL;
    int i, ns = 0;

    if (type == SEL_ALL || type == SEL_NULL) {
	return NULL;
    }

    if (type == SEL_MATRIX) {
	if (sel->m == NULL) {
	    gretl_errmsg_set(_("Range is non-positive!"));
	    *err = E_DATA;
	} else {
	    ns = gretl_vector_get_length(sel->m);
	    if (ns == 0) {
		fprintf(stderr, "selection matrix is %d x %d\n", 
			sel->m->rows, sel->m->cols);
		*err = E_NONCONF;
	    } else if (vec_is_const(sel->m, ns)) {
		ns = 1;
	    }
	}
    } else {
	/* range or element */
	if (sel->range[1] == MSEL_MAX) {
	    sel->range[1] = n;
	}
	if (msel_out_of_bounds(sel->range, n)) {
	    *err = E_DATA;
	} else {
	    ns = sel->range[1] - sel->range[0] + 1;
	    if (ns <= 0) {
		gretl_errmsg_sprintf(_("Range %d to %d is non-positive!"),
				     sel->range[0], sel->range[1]); 
		*err = E_DATA;
	    }
	}
    } 

    if (*err) {
	return NULL;
    }

    slice = gretl_list_new(ns);
    if (slice == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<slice[0]; i++) {
	if (type == SEL_MATRIX) {
	    slice[i+1] = sel->m->val[i];
	} else {
	    slice[i+1] = sel->range[0] + i;
	} 
    }

    for (i=1; i<=slice[0] && !*err; i++) {
	if (slice[i] < 1 || slice[i] > n) {
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 
				 slice[i]);
	    *err = 1;
	}
    }

    if (*err) {
	free(slice);
	slice = NULL;
    }

    return slice;
}

/* catch the case of an implicit column or row specification for
   a sub-matrix of an (n x 1) or (1 x m) matrix; also catch the
   error of giving just one row/col spec for a matrix that has
   more than one row and more than one column
*/

int check_matrix_subspec (matrix_subspec *spec, const gretl_matrix *m)
{
    int err = 0;

    if (spec->type[1] == SEL_NULL) {
	/* we got only one row/col spec */
	if (m->cols == 1) {
	    /* OK: implicitly col = 1 */
	    spec->type[1] = SEL_RANGE;
	    mspec_set_col_index(spec, 1);
	} else if (m->rows == 1) {
	    /* OK: implicitly row = 1, and transfer the single 
	       given spec to the column dimension */
	    spec->type[1] = spec->type[0];
	    if (spec->type[1] == SEL_MATRIX) {
		spec->sel[1].m = spec->sel[0].m;
	    } else {
		spec->sel[1].range[0] = spec->sel[0].range[0];
		spec->sel[1].range[1] = spec->sel[0].range[1];
	    }		
	    spec->type[0] = SEL_RANGE;
	    mspec_set_row_index(spec, 1);
	} else {
	    gretl_errmsg_set(_("Ambiguous matrix index"));
	    err = E_DATA;
	}	    
    }

    if (spec->type[0] == SEL_RANGE && spec->type[1] == SEL_RANGE) {
	if (spec->sel[0].range[0] == spec->sel[0].range[1] &&
	    spec->sel[1].range[0] == spec->sel[1].range[1]) {
	    spec->type[0] = spec->type[1] = SEL_ELEMENT;
	}
    }

    return err;
}

static int get_slices (matrix_subspec *spec, 
		       const gretl_matrix *M)
{
    int err = 0;

    spec->rslice = mspec_make_list(spec->type[0], &spec->sel[0], 
				   M->rows, &err);

    if (!err) {
	spec->cslice = mspec_make_list(spec->type[1], &spec->sel[1], 
				       M->cols, &err);
    }

    return err;
}

int assign_scalar_to_submatrix (gretl_matrix *M, double x,
				matrix_subspec *spec)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int i, err = 0;

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (spec->type[0] == SEL_DIAG) {
	int n = (mr < mc)? mr : mc;

	for (i=0; i<n; i++) {
	    gretl_matrix_set(M, i, i, x);
	}
	return 0;
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	/* parse mspec into lists of affected rows and columns */
	err = get_slices(spec, M);
    }

    if (!err) {
	int sr = (spec->rslice == NULL)? mr : spec->rslice[0];
	int sc = (spec->cslice == NULL)? mc : spec->cslice[0];
	int j, l, k = 0;
	int mi, mj;

	for (i=0; i<sr; i++) {
	    mi = (spec->rslice == NULL)? k++ : spec->rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<sc; j++) {
		mj = (spec->cslice == NULL)? l++ : spec->cslice[j+1] - 1;
		gretl_matrix_set(M, mi, mj, x);
	    }
	}
    }

    return err;
}

static int matrix_insert_diagonal (gretl_matrix *M, 
				   const gretl_matrix *S,
				   int mr, int mc)
{
    int i, n = gretl_vector_get_length(S);
    int k = (mr < mc)? mr : mc;

    if (n != k) {
	return E_NONCONF;
    }

    for (i=0; i<n; i++) {
	gretl_matrix_set(M, i, i, S->val[i]);
    }
    
    return 0;
}

/* @M is the target for partial replacement, @S is the source to
   substitute, and @spec tells how/where to make the
   substitution.
*/

static int matrix_replace_submatrix (gretl_matrix *M,
				     const gretl_matrix *S,
				     matrix_subspec *spec)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int sr = gretl_matrix_rows(S);
    int sc = gretl_matrix_cols(S);
    int sscalar = 0;
    int err = 0;

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (sr > mr || sc > mc) {
	/* the replacement matrix won't fit into M */
	fprintf(stderr, "matrix_replace_submatrix: target is %d x %d but "
		"replacement part is %d x %d\n", mr, mc, sr, sc);
	return E_NONCONF;
    }

    if (spec->type[0] == SEL_DIAG) {
	return matrix_insert_diagonal(M, S, mr, mc);
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	/* parse mspec into lists of affected rows and columns */
	err = get_slices(spec, M);
	if (err) {
	    return err;
	}
    }

#if MDEBUG
    printlist(spec->rslice, "rslice (rows list)");
    printlist(spec->cslice, "cslice (cols list)");
    fprintf(stderr, "orig M = %d x %d, S = %d x %d\n", mr, mc, sr, sc);
#endif

    if (sr == 1 && sc == 1) {
	/* selection matrix is a scalar */
	sscalar = 1;
	sr = (spec->rslice == NULL)? mr : spec->rslice[0];
	sc = (spec->cslice == NULL)? mc : spec->cslice[0];
    } else if (spec->rslice != NULL && spec->rslice[0] != sr) {
	fprintf(stderr, "mspec has %d rows but substitute matrix has %d\n", 
		spec->rslice[0], sr);
	err = E_NONCONF;
    } else if (spec->cslice != NULL && spec->cslice[0] != sc) {
	fprintf(stderr, "mspec has %d cols but substitute matrix has %d\n", 
		spec->cslice[0], sc);
	err = E_NONCONF;
    }

    if (!err) {
	int i, j, l, k = 0;
	int mi, mj;
	double x;

	x = (sscalar)? S->val[0] : 0.0;

	for (i=0; i<sr; i++) {
	    mi = (spec->rslice == NULL)? k++ : spec->rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<sc; j++) {
		mj = (spec->cslice == NULL)? l++ : spec->cslice[j+1] - 1;
		if (!sscalar) {
		    x = gretl_matrix_get(S, i, j);
		}
		gretl_matrix_set(M, mi, mj, x);
	    }
	}
    }

    return err;
}

gretl_matrix *matrix_get_submatrix (const gretl_matrix *M, 
				    matrix_subspec *spec,
				    int prechecked, 
				    int *err)
{
    gretl_matrix *S;
    int r, c;

    if (gretl_is_null_matrix(M)) {
	*err = E_DATA;
	return NULL;
    }

    if (!prechecked) {
	*err = check_matrix_subspec(spec, M);
	if (*err) {
	    return NULL;
	}
    }

    if (spec->type[0] == SEL_DIAG) {
	return gretl_matrix_get_diagonal(M, err);
    } else if (spec->type[0] == SEL_ELEMENT) {
	int i = mspec_get_row_index(spec);
	int j = mspec_get_col_index(spec);
	double x = matrix_get_element(M, i, j, err);

	return (*err)? NULL : gretl_matrix_from_scalar(x);
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	*err = get_slices(spec, M);
	if (*err) {
	    return NULL;
	}
    }

#if MDEBUG
    printlist(spec->rslice, "rslice");
    printlist(spec->cslice, "cslice");
    fprintf(stderr, "M = %d x %d\n", M->rows, M->cols);
#endif

    r = (spec->rslice == NULL)? M->rows : spec->rslice[0];
    c = (spec->cslice == NULL)? M->cols : spec->cslice[0];

    S = gretl_matrix_alloc(r, c);
    if (S == NULL) {
	*err = E_ALLOC;	
    }

    if (S != NULL) {
	int i, j, k, l;
	int mi, mj;
	double x;

	k = 0;
	for (i=0; i<r && !*err; i++) {
	    mi = (spec->rslice == NULL)? k++ : spec->rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<c && !*err; j++) {
		mj = (spec->cslice == NULL)? l++ : spec->cslice[j+1] - 1;
		x = gretl_matrix_get(M, mi, mj);
		gretl_matrix_set(S, i, j, x);
	    }
	}
    }

    if (S != NULL && S->rows == M->rows && gretl_matrix_is_dated(M)) {
	int mt1 = gretl_matrix_get_t1(M);
	int mt2 = gretl_matrix_get_t2(M);

	gretl_matrix_set_t1(S, mt1);
	gretl_matrix_set_t2(S, mt2);
    }

    return S;
}

double matrix_get_element (const gretl_matrix *M, int i, int j,
			   int *err)
{
    double x = NADBL;

    /* The incoming i and j are from userspace, and will
       be 1-based.
    */
    i--;
    j--;

    if (M == NULL) {
	*err = E_DATA;
    } else if (i < 0 || i >= M->rows || j < 0 || j >= M->cols) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 
			     (i < 0 || i >= M->rows)? (i+1) : (j+1));
	*err = 1;
    } else {
	x = gretl_matrix_get(M, i, j);
    }

    return x;
}

gretl_matrix *user_matrix_get_submatrix (const char *name, 
					 matrix_subspec *spec,
					 int *err)
{
    gretl_matrix *S = NULL;
    gretl_matrix *M;

    M = user_var_get_value_by_name(name);
    if (M == NULL) {
	*err = E_UNKVAR;
    } else {
	S = matrix_get_submatrix(M, spec, 0, err);
    }

    return S;
}

/* Look up the existing matrix called @name, and substitute
   the matrix @S for part of the original, as specified by
   @spec.
*/

int user_matrix_replace_submatrix (const char *mname, 
				   const gretl_matrix *S,
				   matrix_subspec *spec)
{
    gretl_matrix *M = user_var_get_value_by_name(mname);

    if (M == NULL) {
	return E_UNKVAR;
    } else {
	return matrix_replace_submatrix(M, S, spec);
    }
}

int umatrix_set_names_from_string (gretl_matrix *M, 
				   const char *s,
				   int byrow)
{
    int n, err = 0;

    n = (byrow)? M->rows : M->cols;

    if (s == NULL || *s == '\0') {
	if (byrow) {
	    gretl_matrix_set_rownames(M, NULL);
	} else {
	    gretl_matrix_set_colnames(M, NULL);
	}
    } else {
	char **S;
	int ns;

	S = gretl_string_split(s, &ns, " \n\t");
	if (S == NULL) {
	    err = E_ALLOC;
	} else if (ns != n) {
	    err = E_NONCONF;
	    strings_array_free(S, ns);
	} else if (byrow) {
	    gretl_matrix_set_rownames(M, S);
	} else {
	    gretl_matrix_set_colnames(M, S);
	} 
    }

    return err;
}

int umatrix_set_names_from_list (gretl_matrix *M, 
				 const int *list,
				 const DATASET *dset,
				 int byrow)
{
    int i, n, err = 0;

    n = (byrow)? M->rows : M->cols;

    if (list == NULL || list[0] == 0) {
	if (byrow) {
	    gretl_matrix_set_rownames(M, NULL);
	} else {
	    gretl_matrix_set_colnames(M, NULL);
	}
    } else if (list[0] != n) {
	err = E_NONCONF;
    } else {
	char **S = strings_array_new(n);

	if (S == NULL) {
	    err = E_ALLOC;
	}

	for (i=0; i<n && !err; i++) {
	    S[i] = gretl_strndup(dset->varname[list[i+1]], 12);
	    if (S[i] == NULL) {
		err = E_ALLOC;
	    }
	}

	if (err) {
	    strings_array_free(S, n);
	} else if (byrow) {
	    gretl_matrix_set_rownames(M, S);
	} else {
	    gretl_matrix_set_colnames(M, S);
	}
    }

    return err;
}

char *user_matrix_get_column_name (const gretl_matrix *M, int col,
				   int *err)
{
    char *ret = NULL;

    if (M == NULL || col < 1 || col > M->cols) {
	*err = E_DATA;
    } else {
	const char **S = gretl_matrix_get_colnames(M);

	if (S == NULL) {
	    ret = gretl_strdup("");
	} else {
	    ret = gretl_strdup(S[col-1]);
	}
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

double 
user_matrix_get_determinant (gretl_matrix *m, int tmpmat, 
			     int f, int *err)
{
    gretl_matrix *R = NULL;
    double d = NADBL;

    if (gretl_is_null_matrix(m)) {
	return d;
    } else if (tmpmat) {
	/* it's OK to overwrite @m */
	R = m;
    } else {
	/* @m should not be over-written! */
	R = gretl_matrix_copy(m);
    }

    if (R != NULL) {
	if (f == F_LDET) {
	    d = gretl_matrix_log_determinant(R, err);
	} else {
	    d = gretl_matrix_determinant(R, err);
	}
	if (R != m) {
	    gretl_matrix_free(R);
	}
    }

    return d;
}

gretl_matrix *user_matrix_matrix_func (gretl_matrix *m, int tmpmat, 
				       int f, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
    } else if (tmpmat) {
	/* it's OK to overwrite @m */
	R = m;
    } else {
	/* @m should not be over-written! */
	R = gretl_matrix_copy(m);
	if (R == NULL) {
	    *err = E_ALLOC;
	}	
    }

    if (R != NULL) {
	if (f == F_CDEMEAN) {
	    gretl_matrix_demean_by_column(R);
	} else if (f == F_CHOL) {
	    *err = gretl_matrix_cholesky_decomp(R);
	} else if (f == F_PSDROOT) {
	    *err = gretl_matrix_psd_root(R);
	} else if (f == F_INVPD) {
	    *err = gretl_invpd(R);
	} else if (f == F_GINV) {
	    *err = gretl_matrix_moore_penrose(R);
	} else if (f == F_INV) {
	    *err = gretl_invert_matrix(R);
	} else if (f == F_UPPER) {
	    *err = gretl_matrix_zero_lower(R);
	} else if (f == F_LOWER) {
	    *err = gretl_matrix_zero_upper(R);
	} else {
	    *err = E_DATA;
	}
	if (*err && R != m) {
	    gretl_matrix_free(R);
	    R = NULL;
	}
    } 
   
    return R;
}

static void matrix_cannibalize (gretl_matrix *targ, gretl_matrix *src)
{
    targ->rows = src->rows;
    targ->cols = src->cols;

    free(targ->val);
    targ->val = src->val;
    src->val = NULL;
}

int matrix_invert_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_copy(m);
    int err = 0;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_invert_matrix(R);
	if (!err) {
	    matrix_cannibalize(m, R);
	}
	gretl_matrix_free(R);
    } 

    return err;
}

int matrix_cholesky_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_copy(m);
    int err = 0;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_cholesky_decomp(R);
	if (!err) {
	    matrix_cannibalize(m, R);
	}
	gretl_matrix_free(R);
    } 

    return err;
}

int matrix_transpose_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_copy_transpose(m);
    int err = 0;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	matrix_cannibalize(m, R);
	gretl_matrix_free(R);
    }

    return err;
}

int matrix_XTX_in_place (gretl_matrix *m)
{
    gretl_matrix *R = gretl_matrix_alloc(m->cols, m->cols);
    int err;

    if (R == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_multiply_mod(m, GRETL_MOD_TRANSPOSE,
					m, GRETL_MOD_NONE,
					R, GRETL_MOD_NONE);
    }

    if (!err) {
	matrix_cannibalize(m, R);
    }

    gretl_matrix_free(R);

    return err;
}

gretl_matrix *user_matrix_vec (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else {
	R = gretl_matrix_alloc(m->rows * m->cols, 1);
	if (R != NULL) {
	    gretl_matrix_vectorize(R, m);
	} 
    }

    if (R == NULL) {
	*err = E_ALLOC;
    }

    return R;
}

gretl_matrix *user_matrix_vech (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else if (m->rows != m->cols) {
	*err = E_NONCONF;
    } else {
	int n = m->rows;
	int k = n * (n + 1) / 2;

	R = gretl_matrix_alloc(k, 1);
	if (R != NULL) {
	    *err = gretl_matrix_vectorize_h(R, m);
	}
    } 

    if (R == NULL && !*err) {
	*err = E_ALLOC;
    }

    return R;
}

gretl_matrix *user_matrix_unvech (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else if (m->cols != 1) {
	*err = E_NONCONF;
    } else {
	int n = (int) ((sqrt(1.0 + 8.0 * m->rows) - 1.0) / 2.0);

	R = gretl_matrix_alloc(n, n);
	if (R != NULL) {
	    *err = gretl_matrix_unvectorize_h(R, m);
	} 
    }

    if (R == NULL && !*err) {
	*err = E_ALLOC;
    }

    return R;
}

static int 
real_user_matrix_QR_decomp (const gretl_matrix *m, gretl_matrix **Q, 
			    gretl_matrix **R)
{
    int mc = gretl_matrix_cols(m);
    int err = 0;

    *Q = gretl_matrix_copy(m);

    if (*Q == NULL) {
	err = E_ALLOC;
    } else if (R != NULL) {
	*R = gretl_matrix_alloc(mc, mc);
	if (*R == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_QR_decomp(*Q, (R == NULL)? NULL : *R);
    }

    if (err) {
	gretl_errmsg_set(_("Matrix decomposition failed"));
	gretl_matrix_free(*Q);
	*Q = NULL;
	if (R != NULL) {
	    gretl_matrix_free(*R);
	    *R = NULL;
	}
    }

    return err;
}

#define nullarg(s) (s == NULL || !strcmp(s, "null"))

gretl_matrix *
user_matrix_QR_decomp (const gretl_matrix *m, const char *rname, int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix **pR = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (!nullarg(rname)) {
	if (get_matrix_by_name(rname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), rname);
	    *err = E_UNKVAR;
	} else {
	    pR = &R;
	}
    }

    if (!*err) {
	*err = real_user_matrix_QR_decomp(m, &Q, pR);
    }

    if (!*err && R != NULL) {
	user_matrix_replace_matrix_by_name(rname, R);
    }

    return Q;
}

static int revise_SVD_V (gretl_matrix **pV, int r, int c)
{
    gretl_matrix *V;
    double x;
    int i, j;

    V = gretl_matrix_alloc(r, c);
    if (V == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    x = gretl_matrix_get((*pV), i, j);
	    gretl_matrix_set(V, i, j, x);
	}
    }

    gretl_matrix_free(*pV);
    *pV = V;

    return 0;
}

gretl_matrix *user_matrix_SVD (const gretl_matrix *m, 
			       const char *uname, 
			       const char *vname, 
			       int *err)
{
    gretl_matrix *U = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix **pU = NULL;
    gretl_matrix **pV = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (!nullarg(uname)) {
	if (get_matrix_by_name(uname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), uname);
	    *err = E_UNKVAR;
	} else {
	    pU = &U;
	}
    }

    if (!*err && !nullarg(vname)) {
	if (get_matrix_by_name(vname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), vname);
	    *err = E_UNKVAR;
	} else {
	    pV = &V;
	}
    }

    if (!*err) {
	*err = gretl_matrix_SVD(m, pU, &S, pV);
    }

    if (!*err && (U != NULL || V != NULL)) {
	int tall = m->rows - m->cols;
	int minrc = (m->rows > m->cols)? m->cols : m->rows;

	if (U != NULL) {
	    if (tall > 0) {
		*err = gretl_matrix_realloc(U, m->rows, minrc);
	    }
	    if (!*err) {
		user_matrix_replace_matrix_by_name(uname, U);
	    }
	}
	if (V != NULL) {
	    if (tall < 0) {
		*err = revise_SVD_V(&V, minrc, m->cols);
	    } 
	    if (!*err) {
		user_matrix_replace_matrix_by_name(vname, V);
	    }
	}
    }

    return S;
}

/* Here we're looking for a matrix that was passed by address to
   mols() or similar. If we can't find it, that's an error.
   Otherwise, set @newmat depending on whether the existing matrix has
   the appropriate dimensions to accept the result (@newmat = 0) or
   not (in which case we need to allocate a new gretl_matrix to
   replace the old one, and we set @newmat = 1).
*/

static gretl_matrix *get_ols_matrix (const char *mname, 
				     int r, int c,
				     int *newmat, 
				     int *err)
{
    gretl_matrix *m = get_matrix_by_name(mname);

    if (m == NULL) {
	gretl_errmsg_sprintf(_("'%s': no such matrix"), mname);
	*err = E_UNKVAR;
    } else if (m->rows != r || m->cols != c) {
	m = gretl_matrix_alloc(r, c);
	if (m == NULL) {
	    *err = E_ALLOC;
	} else {
	    *newmat = 1;
	}
    }

    return m;
}

static int check_vcv_arg (const char *mname)
{
    gretl_matrix *m = get_matrix_by_name(mname);

    if (m == NULL) {
	gretl_errmsg_sprintf(_("'%s': no such matrix"), mname);
	return E_UNKVAR;
    } else {
	return 0;
    }
}

gretl_matrix *user_matrix_ols (const gretl_matrix *Y, 
			       const gretl_matrix *X, 
			       const char *Uname, 
			       const char *Vname, 
			       gretlopt opt,
			       int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *U = NULL;
    gretl_matrix *V = NULL;
    int newU = 0, newV = 0;
    double s2, *ps2 = NULL;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || X == NULL) {
	*err = E_DATA;
	return NULL;
    }

    T = Y->rows;
    k = X->cols;
    g = Y->cols;

    if (X->rows != T) {
	*err = E_NONCONF;
	return NULL;
    }

    if (g > 1 && (opt & OPT_M)) {
	/* multiple precision: we accept only one y var */
	*err = E_DATA;
	return NULL;
    }

    if (!nullarg(Uname)) {
	U = get_ols_matrix(Uname, T, g, &newU, err);
	if (*err) {
	    return NULL;
	}
    } 

    if (!nullarg(Vname)) {
	if (g > 1) {
	    /* multiple dependent variables */
	    *err = check_vcv_arg(Vname);
	    if (!*err) {
		newV = 1;
	    }
	} else {
	    /* a single dependent variable */
	    int nv = g * k;

	    V = get_ols_matrix(Vname, nv, nv, &newV, err);
	    if (!*err) {
		ps2 = &s2;
	    }
	}
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	if (gretl_is_null_matrix(X)) {
	    if (U != NULL) {
		gretl_matrix_copy_values(U, Y);
	    }
	    if (!nullarg(Vname)) {
		V = gretl_null_matrix_new();
		if (V == NULL) {
		    *err = E_ALLOC;
		}
	    }
	} else if (g == 1) {
	    if (opt & OPT_M) {
		/* use multiple precision */
		*err = gretl_matrix_mp_ols(Y, X, B, V, U, ps2);
	    } else {
		*err = gretl_matrix_ols(Y, X, B, V, U, ps2);
	    }
	} else {
	    if (newV) {
		/* note: "V" will actually be (X'X)^{-1} */
		*err = gretl_matrix_multi_ols(Y, X, B, U, &V);
	    } else {
		*err = gretl_matrix_multi_ols(Y, X, B, U, NULL);
	    }
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
	if (newU) gretl_matrix_free(U);
	if (newV) gretl_matrix_free(V);
    } else {
	if (newU) {
	    user_matrix_replace_matrix_by_name(Uname, U);
	}
	if (newV) {
	    user_matrix_replace_matrix_by_name(Vname, V);
	}	
    }

    return B;
}

gretl_matrix *user_matrix_rls (const gretl_matrix *Y, 
			       const gretl_matrix *X,
			       const gretl_matrix *R,
			       const gretl_matrix *Q,
			       const char *Uname, 
			       const char *Vname, 
			       int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *U = NULL;
    gretl_matrix *V = NULL;
    int newU = 0, newV = 0;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || gretl_is_null_matrix(X)) {
	*err = E_DATA;
	return NULL;
    }

    T = Y->rows;
    k = X->cols;
    g = Y->cols;

    if (X->rows != T) {
	*err = E_NONCONF;
	return NULL;
    }

    if (!nullarg(Uname)) {
	U = get_ols_matrix(Uname, T, g, &newU, err);
	if (*err) {
	    return NULL;
	}
    } 

    if (!nullarg(Vname)) {
	*err = check_vcv_arg(Vname);
	if (!*err) {
	    newV = 1;
	}
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	if (newV) {
	    /* note: "V" will actually be M (X'X)^{-1} M' */
	    *err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B, U, &V);
	} else {
	    *err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B, U, NULL);
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
	if (newU) gretl_matrix_free(U);
	if (newV) gretl_matrix_free(V);
    } else {
	if (newU) {
	    user_matrix_replace_matrix_by_name(Uname, U);
	}
	if (newV) {
	    user_matrix_replace_matrix_by_name(Vname, V);
	}	
    }

    return B;
}

static void maybe_eigen_trim (gretl_matrix *E)
{
    double x;
    int i, allreal = 1;

    for (i=0; i<E->rows; i++) {
	x = gretl_matrix_get(E, i, 1);
	if (x != 0.0) {
	    allreal = 0;
	    break;
	}
    }

    if (allreal) {
	gretl_matrix_reuse(E, -1, 1);
    }
}

gretl_matrix *
user_matrix_eigen_analysis (const gretl_matrix *m, const char *rname, int symm,
			    int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *E = NULL;
    int vecs = 0;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_xna_check(m)) {
	*err = E_NAN;
	return NULL;
    }

    if (!nullarg(rname)) {
	vecs = 1;
	if (get_matrix_by_name(rname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), rname);
	    *err = E_UNKVAR;
	    return NULL;
	}
    }

    C = gretl_matrix_copy(m);
    if (C == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	if (symm) {
	    E = gretl_symmetric_matrix_eigenvals(C, vecs, err);
	} else {
	    E = gretl_general_matrix_eigenvals(C, vecs, err);
	    if (E != NULL && E->cols == 2) {
		maybe_eigen_trim(E);
	    }
	}
    }

    if (!*err && vecs) {
	user_matrix_replace_matrix_by_name(rname, C);
    }

    if (!vecs) {
	gretl_matrix_free(C);
    }

    return E;
}

gretl_matrix *user_gensymm_eigenvals (const gretl_matrix *A, 
				      const gretl_matrix *B,
				      const char *rname,
				      int *err)
{
    gretl_matrix *E = NULL, *V = NULL;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_xna_check(A) || gretl_matrix_xna_check(B)) {
	*err = E_NAN;
	return NULL;
    }

    if (!nullarg(rname)) {
	if (get_matrix_by_name(rname) == NULL) {
	    gretl_errmsg_sprintf(_("'%s': no such matrix"), rname);
	    *err = E_UNKVAR;
	    return NULL;
	} else {
	    V = gretl_matrix_alloc(B->cols, A->rows);
	    if (V == NULL) {
		*err = E_ALLOC;
		return NULL;
	    }
	}
    }

    E = gretl_gensymm_eigenvals(A, B, V, err);

    if (V != NULL) {
	if (*err) {
	    gretl_matrix_free(V);
	} else {
	    user_matrix_replace_matrix_by_name(rname, V);
	}
    }

    return E;
}
