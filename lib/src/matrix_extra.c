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
#include "usermat.h"
#include "libset.h"
#include "gretl_cmatrix.h"
#include "matrix_extra.h"
#include "swap_bytes.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <errno.h>

/**
 * SECTION:matrix_extra
 * @short_description: more matrix operations
 * @title: Matrix extra
 * @include: gretl/libgretl.h, gretl/matrix_extra.h
 *
 * Functions that convert between gretl_matrix and other
 * datatypes; also printing of gretl_matrices and some
 * functions relating to masking and cutting rows and/or
 * columns of matrices.
 */

/**
 * gretl_vector_from_array:
 * @x: pointer to array of elements.
 * @n: number of elements.
 * @mod: modifier flag: either %GRETL_MOD_NONE, or %GRETL_MOD_SQUARE
 * to use the squares of the elements of @x.
 *
 * Returns: pointer to a newly allocated gretl_vector containing
 * the elements of x (or their squares), or NULL on failure.
 * Missing values in @x are skipped.
 */

gretl_vector *
gretl_vector_from_array (const double *x, int n, GretlMatrixMod mod)
{
    gretl_matrix *v;
    double xi;

    v = gretl_column_vector_alloc(n);

    if (v != NULL) {
	int i = 0, j = 0;

	while (j < n) {
	    xi = x[i++];
	    if (!na(xi)) {
		if (mod == GRETL_MOD_SQUARE) {
		    v->val[j] = xi * xi;
		} else {
		    v->val[j] = xi;
		}
		j++;
	    }
	}
    }

    return v;
}

/**
 * gretl_vector_from_series:
 * @x: series from data array.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Returns: a newly allocated gretl_vector containing the values
 * of the given data series for the given range, or NULL on failure.
 * Any missing values in the input array are preserved as NaNs in
 * the output.
 */

gretl_vector *gretl_vector_from_series (const double *x,
					int t1, int t2)
{
    gretl_matrix *v = NULL;
    int n = t2 - t1 + 1;

    if (n > 0) {
	size_t sz = n * sizeof *x;

	v = gretl_column_vector_alloc(n);
	if (v != NULL) {
	    memcpy(v->val, x + t1, sz);
	    gretl_matrix_set_t1(v, t1);
	    gretl_matrix_set_t2(v, t2);
	}
    }

    return v;
}

/**
 * gretl_matrix_from_2d_array:
 * @X: two-dimensional array of doubles.
 * @rows: number of rows in target matrix.
 * @cols: number of columns in target matrix.
 *
 * Returns: allocated gretl_matrix, the elements of which are set to
 * the values in @X, or NULL on allocation failure.
 */

gretl_matrix *gretl_matrix_from_2d_array (const double **X,
					  int rows, int cols)
{
    gretl_matrix *m;
    int i, j, k = 0;

    m = gretl_matrix_alloc(rows, cols);

    if (m != NULL) {
	for (j=0; j<cols; j++) {
	    for (i=0; i<rows; i++) {
		m->val[k++] = X[j][i];
	    }
	}
    }

    return m;
}

/**
 * gretl_matrix_from_scalar:
 * @x: scalar to be "promoted".
 *
 * Returns: allocated 1 x 1 gretl_matrix, the single element
 * of which is set to @x, or NULL on allocation failure.
 */

gretl_matrix *gretl_matrix_from_scalar (double x)
{
    gretl_matrix *m = gretl_matrix_alloc(1, 1);

    m->val[0] = x;

    return m;
}

static int count_selection (const char *s, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
	if (s[i] != 0) c++;
    }

    return c;
}

/**
 * gretl_vcv_matrix_from_model:
 * @pmod: pointer to model
 * @select: char array indicating which rows and colums to select
 * (or NULL for the full matrix).
 * @err: location to receive error code.
 *
 * Produces all or part of the covariance matrix for @pmod
 * in the form of a gretl_matrix.  Storage is allocated, to be freed
 * by the caller.  If @select is not NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows (and columns), and zero elements otherwise.
 *
 * Returns: the covariance matrix, or NULL on error.
 */

gretl_matrix *
gretl_vcv_matrix_from_model (MODEL *pmod, const char *select, int *err)
{
    gretl_matrix *vcv;
    int i, j, idx, nc;
    int ii, jj;
    int k = pmod->ncoeff;

    /* first ensure the model _has_ a vcv */
    *err = makevcv(pmod, pmod->sigma);
    if (*err) {
	fprintf(stderr, "gretl_vcv_matrix_from_model: no $vcv!\n");
	return NULL;
    }

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }

    if (nc == 0) {
	*err = E_DATA;
	return NULL;
    }

    vcv = gretl_matrix_alloc(nc, nc);
    if (vcv == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    ii = 0;
    for (i=0; i<k; i++) {
	if (select != NULL && !select[i]) {
	    continue;
	}
	jj = 0;
	for (j=0; j<=i; j++) {
	    if (select != NULL && !select[j]) {
		continue;
	    }
	    idx = ijton(i, j, pmod->ncoeff);
	    gretl_matrix_set(vcv, ii, jj, pmod->vcv[idx]);
	    if (jj != ii) {
		gretl_matrix_set(vcv, jj, ii, pmod->vcv[idx]);
	    }
	    jj++;
	}
	ii++;
    }

    return vcv;
}

/**
 * gretl_coeff_vector_from_model:
 * @pmod: pointer to model
 * @select: char array indicating which rows to select
 * (or NULL for the full vector).
 * @err: location to receive error code.
 *
 * Produces all or part of the coefficient vector for @pmod
 * in the form of a gretl column vector.  Storage is allocated, to be freed
 * by the caller.  If @select is non-NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows and zero elements otherwise.
 *
 * Returns: the coefficient vector, or NULL on error.
 */

gretl_vector *
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select, int *err)
{
    gretl_vector *b;
    int i, j, nc;
    int k = pmod->ncoeff;

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }

    if (nc == 0) {
	*err = E_DATA;
	return NULL;
    }

    b = gretl_column_vector_alloc(nc);
    if (b == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    j = 0;
    for (i=0; i<k; i++) {
	if (select != NULL && !select[i]) {
	    continue;
	}
	b->val[j++] = pmod->coeff[i];
    }

    return b;
}

/**
 * gretl_covariance_matrix_from_varlist:
 * @list: list of variables by ID number.
 * @dset: dataset struct.
 * @means: pointer to pick up vector of means, or NULL to discard.
 * @errp: pointer to receive non-zero error code in case of
 * failure, or NULL.
 *
 * Returns: the variance-covariance matrix of the listed variables
 * (over the currently defined data sample), or NULL in case of
 * failure.
 */

gretl_matrix *
gretl_covariance_matrix_from_varlist (const int *list,
				      const DATASET *dset,
				      gretl_matrix **means,
				      int *errp)
{
    gretl_matrix *V;
    gretl_vector *xbar;
    int k = list[0];
    int t, i, j, vi, vj;
    double vv, x, y;
    int n, err = 0;

    V = gretl_matrix_alloc(k, k);
    xbar = gretl_vector_alloc(k);

    if (V == NULL || xbar == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<k && !err; i++) {
	xbar->val[i] = gretl_mean(dset->t1, dset->t2, dset->Z[list[i+1]]);
	if (na(xbar->val[i])) {
	    err = E_DATA;
	}
    }

    for (i=0; i<k && !err; i++) {
	vi = list[i+1];
	for (j=i; j<k; j++) {
	    vj = list[j+1];
	    vv = 0.0;
	    n = 0;
	    for (t=dset->t1; t<=dset->t2; t++) {
		x = dset->Z[vi][t];
		y = dset->Z[vj][t];
		if (na(x) || na(y)) {
		    continue;
		}
		vv += (x - xbar->val[i]) * (y - xbar->val[j]);
		n++;
	    }
	    if (n < 2) {
		err = E_DATA;
		vv = NADBL;
	    } else {
		vv /= (n - 1); /* or plain n? */
	    }
	    gretl_matrix_set(V, i, j, vv);
	    gretl_matrix_set(V, j, i, vv);
	}
    }

 bailout:

    if (means != NULL && !err) {
	*means = xbar;
    } else {
	gretl_vector_free(xbar);
    }

    if (err) {
	gretl_matrix_free(V);
	V = NULL;
	if (errp != NULL) {
	    *errp = err;
	}
    }

    return V;
}

/**
 * gretl_matrix_row_to_array:
 * @m: source matrix.
 * @i: the row from which values should be copied.
 * @x: array of doubles large enough to hold a row from @m.
 *
 * Copies the values from row @i of matrix @m into the array
 * @x, which should already be allocated to the correct size.
 *
 * Returns: 0 on sucess, non-zero if the row is out of bounds.
 */

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x)
{
    int j, err = 0;

    if (i < 0 || i >= gretl_matrix_rows(m)) {
	err = E_DATA;
    } else {
	for (j=0; j<m->cols; j++) {
	    x[j] = gretl_matrix_get(m, i, j);
	}
    }

    return err;
}

/**
 * gretl_matrix_get_columns:
 * @m: source matrix.
 * @err: location to receive error code.
 *
 * Constructs a two-dimensional array in which each sub-array
 * is a pointer to a column of @m. The content of these arrays
 * belongs to @m and must not be freed; only the returned
 * pointer itself should be freed.
 *
 * Returns: the constructed array on success or NULL on failure.
 */

double **gretl_matrix_get_columns (const gretl_matrix *m, int *err)
{
    double **X = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
    } else {
	double *val = m->val;
	int j;

	X = doubles_array_new(m->cols, 0);

	if (X == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (j=0; j<m->cols; j++) {
		X[j] = val;
		val += m->rows;
	    }
	}
    }

    return X;
}

static int get_mask_count (const char *mask, int n)
{
    int i, k = 0;

    for (i=0; i<n; i++) {
	if (mask[i]) k++;
    }

    return k;
}

static void add_dataset_colnames (gretl_matrix *M,
				  const int *list,
				  const DATASET *dset)
{
    int i, vi, nv;
    char **S = NULL;
    int err = 0;

    if (list != NULL) {
	nv = list[0];
    } else {
	/* all series except "const" */
	nv = dset->v - 1;
    }

    /* we won't treat errors here as fatal */

    if (nv != M->cols) {
	/* "can't happen" */
	return;
    }

    S = strings_array_new(nv);
    if (S == NULL) {
	return;
    }

    for (i=1; i<=nv; i++) {
	vi = list == NULL ? i : list[i];
	S[i-1] = gretl_strdup(dset->varname[vi]);
	if (S[i-1] == NULL) {
	    strings_array_free(S, nv);
	    return;
	}
    }

    err = gretl_matrix_set_colnames(M, S);
    if (err) {
	strings_array_free(S, nv);
    }
}

static gretl_matrix *
real_gretl_matrix_data_subset (const int *list,
			       const DATASET *dset,
			       int t1, int t2,
			       const char *mask,
			       int op, int *err)
{
    gretl_matrix *M = NULL;
    int T, Tmax = t2 - t1 + 1;
    int do_obs_info = 1;
    int k, mt1 = 0, mt2 = 0;
    int skip;
    int j, vj, s, t;

    if (list != NULL) {
	k = list[0];
    } else {
	k = dset->v - 1;
    }

    if (t1 < 0 && t2 < 0) {
	/* use full dataset */
	t1 = 0;
	t2 = dset->n - 1;
	T = Tmax = dset->n;
	do_obs_info = 0;
    }

    if (k <= 0 || Tmax <= 0 || t2 >= dset->n) {
	*err = E_DATA;
	return NULL;
    }

    *err = 0;
    T = Tmax;

    if (mask != NULL) {
	T -= get_mask_count(mask, Tmax);
    } else if (op == M_MISSING_TRIM) {
	*err = list_adjust_sample(list, &t1, &t2, dset, NULL);
	if (*err) {
	    return NULL;
	} else {
	    Tmax = T = t2 - t1 + 1;
	}
    } else if (op == M_MISSING_SKIP || op == M_MISSING_ERROR) {
	for (t=t1; t<=t2 && !*err; t++) {
	    for (j=1; j<=k; j++) {
		vj = list == NULL ? j : list[j];
		if (na(dset->Z[vj][t])) {
		    if (op == M_MISSING_SKIP) {
			T--;
		    } else {
			*err = E_MISSDATA;
			return NULL;
		    }
		    break;
		}
	    }
	}
    }

    if (T < 0) {
	/* "can't happen" */
	*err = E_DATA;
    } else {
	M = gretl_matrix_alloc(T, k);
	if (M == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err || T == 0) {
	return M;
    }

    if (do_obs_info) {
	mt1 = t1; mt2 = t2;
    }

    s = 0;
    for (t=t1; t<=t2; t++) {
	skip = 0;
	if (mask != NULL) {
	    skip = mask[t - t1];
	} else if (op == M_MISSING_SKIP) {
	    for (j=1; j<=k; j++) {
		vj = list == NULL ? j : list[j];
		if (na(dset->Z[vj][t])) {
		    skip = 1;
		    break;
		}
	    }
	}
	if (!skip) {
	    for (j=1; j<=k; j++) {
		vj = list == NULL ? j : list[j];
		gretl_matrix_set(M, s, j-1, dset->Z[vj][t]);
	    }
	    if (s == 0) {
		mt1 = t;
	    } else if (s == T - 1) {
		mt2 = t;
		break;
	    }
	    s++;
	}
    }

    if (!*err && (mask != NULL || op == M_MISSING_OK)) {
	int n = k * T;

	for (j=0; j<n && !*err; j++) {
	    if (na(M->val[j])) {
		if (op != M_MISSING_OK) {
		    *err = E_MISSDATA;
		}
	    }
	}
    }

    if (*err) {
	gretl_matrix_free(M);
	M = NULL;
    } else {
	if (do_obs_info && T == mt2 - mt1 + 1) {
	    gretl_matrix_set_t1(M, mt1);
	    gretl_matrix_set_t2(M, mt2);
	}
	if (dset->varname != NULL) {
	    add_dataset_colnames(M, list, dset);
	}
    }

    return M;
}

/**
 * gretl_matrix_data_subset:
 * @list: list of series to process, or %NULL to include
 * all series except "const".
 * @dset: pointer to dataset struct.
 * @t1: starting observation.
 * @t2: ending observation.
 * @missop: how to handle missing observations.
 * @err: location to receive error code.
 *
 * Creates a matrix holding the subset of variables from
 * @dset specified by @list, over the sample range @t1 to @t2,
 * inclusive. Series are in columns. The @missop flag can
 * be %M_MISSING_OK to indicate that it's OK to include
 * missing values in the matrix, %M_MISSING_ERROR (it's an
 * error if any missing values are found), or %M_MISSING_SKIP
 * (observations with any missing values are omitted from the
 * matrix).
 *
 * Returns: allocated matrix or NULL on failure.
 */

gretl_matrix *gretl_matrix_data_subset (const int *list,
					const DATASET *dset,
					int t1, int t2,
					int missop,
					int *err)
{
    return real_gretl_matrix_data_subset(list, dset, t1, t2, NULL,
					 missop, err);
}

/**
 * gretl_matrix_data_subset_masked:
 * @list: list of variable to process.
 * @dset: dataset struct.
 * @t1: starting observation.
 * @t2: ending observation.
 * @mask: missing observations mask, or NULL.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @dset specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns.  @mask should be an
 * array of char of length (@t2 - @t1 + 1) with 1s in the positions
 * of observations to exclude from the subset and zeros elsewhere.
 * This apparatus can be used to exclude missing observations.
 *
 * Returns: allocated matrix or NULL on failure.
 */

gretl_matrix *
gretl_matrix_data_subset_masked (const int *list,
				 const DATASET *dset,
				 int t1, int t2,
				 const char *mask,
				 int *err)
{
    return real_gretl_matrix_data_subset(list, dset, t1, t2, mask,
					 M_MISSING_ERROR, err);
}

/* count the number of selected rows in the current
   sample range */

static int mmask_row_count (const gretl_matrix *mask,
			    const DATASET *dset)
{
    int t, m = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (mask->val[t] != 0) {
	    m++;
	}
    }

    return m;
}

/**
 * gretl_matrix_data_subset_special:
 * @list: list of variables to process.
 * @dset: dataset struct.
 * @mmask: matrix holding desired observations mask.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @dset specified by @list, using the observations (data rows)
 * selected by non-zero elements in @mmask. This is designed
 * to support the gretl "libset" variable %matrix_mask.
 *
 * The length of the vector @mmask must equal the number of
 * observations in the dataset, the member %n of @dset.
 *
 * Returns: allocated matrix or NULL on failure.
 */

gretl_matrix *
gretl_matrix_data_subset_special (const int *list,
				  const DATASET *dset,
				  const gretl_matrix *mmask,
				  int *err)
{
    int n = gretl_vector_get_length(mmask);
    gretl_matrix *X = NULL;

    if (list == NULL || n != dset->n) {
	*err = E_DATA;
    } else if (list[0] == 0) {
	X = gretl_null_matrix_new();
    } else {
	int T = mmask_row_count(mmask, dset);
	int k = list[0];

	if (T == 0) {
	    X = gretl_null_matrix_new();
	} else {
	    X = gretl_matrix_alloc(T, k);
	}

	if (X != NULL && T > 0) {
	    const double *xi;
	    int i, s, t;

	    for (i=0; i<k; i++) {
		xi = dset->Z[list[i+1]];
		s = 0;
		for (t=dset->t1; t<=dset->t2; t++) {
		    if (mmask->val[t] != 0) {
			if (s == 0) {
			    gretl_matrix_set_t1(X, t);
			} else if (s == T - 1) {
			    gretl_matrix_set_t2(X, t);
			}
			gretl_matrix_set(X, s++, i, xi[t]);
		    }
		}
	    }
	}
    }

    if (X == NULL && *err == 0) {
	*err = E_ALLOC;
    }

    return X;
}

static void check_matrix_varnames (DATASET *dset)
{
    int i, j, err = 0;

    for (i=1; i<dset->v; i++) {
	for (j=1; j<i; j++) {
	    if (!strcmp(dset->varname[i], dset->varname[j])) {
		err = 1;
		break;
	    }
	}
    }

    if (err) {
	for (i=1; i<dset->v; i++) {
	    sprintf(dset->varname[i], "v%d", i);
	}
    }
}

/**
 * gretl_dataset_from_matrix:
 * @m: source matrix.
 * @list: list of columns (1-based) to include, or NULL to
 * include all columns.
 * @opt: may include OPT_B to attempt "borrowing" of data;
 * OPT_N to use plain numbers as variable names; OPT_R to
 * use row-names as observation markers, if present; OPT_S
 * when called from the "store" command.
 * @err: location to receive error code.
 *
 * Creates a gretl dataset from matrix @m, either using the
 * columns specified in @list or using all columns if @list
 * is NULL.
 *
 * Returns: pointer to new dataset information struct on success,
 * or NULL on failure.
 */

DATASET *gretl_dataset_from_matrix (const gretl_matrix *m,
				    const int *list,
				    gretlopt opt,
				    int *err)
{
    DATASET *dset = NULL;
    const char **cnames = NULL;
    const char **rnames = NULL;
    int use_rownames;
    int i, t, col, nv, T;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    use_rownames = opt & (OPT_R | OPT_S);

    T = m->rows;
    nv = m->cols;

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    col = list[i];
	    if (col < 1 || col > nv) {
		gretl_errmsg_sprintf("Variable number %d is out of bounds", col);
		*err = E_DATA;
		break;
	    }
	}
	nv = list[0];
    }

    if (!*err) {
	gretlopt dsopt = OPT_NONE;

	if (opt & OPT_B) {
	    dsopt |= OPT_B;
	}
	if (use_rownames) {
	    rnames = gretl_matrix_get_rownames(m);
	    if (rnames != NULL) {
		dsopt |= OPT_M;
	    }
	}
	dset = create_auxiliary_dataset(nv + 1, T, dsopt);
	if (dset == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	return NULL;
    }

    cnames = gretl_matrix_get_colnames(m);

    for (i=1; i<=nv; i++) {
	double *src;

	col = (list != NULL)? list[i] - 1 : i - 1;
	src = m->val + T * col;
	if (opt & OPT_B) {
	    /* "borrowing" */
	    dset->Z[i] = src;
	} else {
	    /* copying */
	    memcpy(dset->Z[i], src, T * sizeof *src);
	}
	if (cnames != NULL) {
	    if (gretl_normalize_varname(dset->varname[i], cnames[col], 0, i)) {
		series_record_display_name(dset, i, cnames[col]);
	    }
	} else if (opt & OPT_N) {
	    sprintf(dset->varname[i], "%d", col + 1);
	} else if (opt & OPT_R) {
	    sprintf(dset->varname[i], "v%d", i);
	} else {
	    sprintf(dset->varname[i], "col%d", col + 1);
	}
    }

    if (rnames != NULL) {
	for (t=0; t<T; t++) {
	    dset->S[t][0] = '\0';
	    strncat(dset->S[t], rnames[t], OBSLEN-1);
	}
    }

    if (cnames != NULL && (opt & OPT_S)) {
	/* we must have valid varnames for "store" */
	check_matrix_varnames(dset);
    }

    return dset;
}

/**
 * gretl_dataset_from_matrices:
 * @m1: first source matrix.
 * @m2: second source matrix.
 * @err: location to receive error code.
 *
 * Creates a gretl dataset from the side-by-side concatenation
 * of @m1 and @m2. Storage for the values is "borrowed" from the
 * columns of the matrices. Column names are used to label the
 * series, if present; otherwise the labels are "v1", "v2" and so on.
 *
 * The two matrices must have the same number of rows. Rownames on
 * either matrix are disregarded.
 *
 * Returns: pointer to new dataset on success, or NULL on failure.
 */

DATASET *gretl_dataset_from_matrices (const gretl_matrix *m1,
				      const gretl_matrix *m2,
				      int *err)
{
    DATASET *dset = NULL;
    const char **cnames1 = NULL;
    const char **cnames2 = NULL;
    const char **cnames = NULL;
    const gretl_matrix *m;
    int i, col, nv, T;

    if (gretl_is_null_matrix(m1) ||
	gretl_is_null_matrix(m2)) {
	*err = E_DATA;
	return NULL;
    } else if (m1->rows != m2->rows) {
	*err = E_NONCONF;
	return NULL;
    }

    T = m1->rows;
    nv = m1->cols + m2->cols;

    dset = create_auxiliary_dataset(nv + 1, T, OPT_B);
    if (dset == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    cnames1 = gretl_matrix_get_colnames(m1);
    cnames2 = gretl_matrix_get_colnames(m2);
    m = m1;
    cnames = cnames1;
    col = 0;

    for (i=1; i<=nv; i++) {
	dset->Z[i] = m->val + T * col;
	if (cnames != NULL) {
	    if (gretl_normalize_varname(dset->varname[i], cnames[col], 0, i)) {
		series_record_display_name(dset, i, cnames[col]);
	    }
	} else {
	    sprintf(dset->varname[i], "v%d", i);
	}
	if (m == m1 && i == m1->cols) {
	    m = m2;
	    cnames = cnames2;
	    col = 0;
	} else {
	    col++;
	}
    }

    return dset;
}

int write_matrix_as_dataset (const char *fname,
			     gretlopt opt,
			     PRN *prn)
{
    const char *mname;
    gretl_matrix *m;
    DATASET *mdset;
    int err = 0;

    mname = get_optval_string(STORE, OPT_A);
    m = get_matrix_by_name(mname);
    if (m == NULL) {
	return E_DATA;
    }

    mdset = gretl_dataset_from_matrix(m, NULL, OPT_S, &err);
    if (!err) {
	opt &= ~OPT_A;
	err = write_data(fname, NULL, mdset, opt, prn);
    }

    destroy_dataset(mdset);

    return err;
}

/**
 * gretl_plotfit_matrices:
 * @yvar: the y variable.
 * @xvar: the x variable.
 * @fit: type of fit sought.
 * @t1: starting observation.
 * @t2: ending observation.
 * @py: location to receive y vector.
 * @pX: location to receive X matrix.
 *
 * Creates a vector y and matrix X based on the input @yvar,
 * @xvar and @fit, using the given sample range.  An observation
 * is skipped if either @yvar or @xvar is missing at that
 * observation.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plotfit_matrices (const double *yvar, const double *xvar,
			    FitType fit, int t1, int t2,
			    gretl_matrix **py, gretl_matrix **pX)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    char *mask = NULL;
    double xt;
    int T = t2 - t1 + 1;
    int n = 0;
    int i, j, k, s, t;
    int err = 0;

    if (T <= 0) {
	return E_DATA;
    }

    if (fit == PLOT_FIT_LOGLIN && !gretl_ispositive(t1, t2, yvar, 1)) {
	gretl_errmsg_set(_("Non-positive values encountered"));
	return E_DATA;
    } else if (fit == PLOT_FIT_LINLOG && !gretl_ispositive(t1, t2, xvar, 1)) {
	gretl_errmsg_set(_("Non-positive values encountered"));
	return E_DATA;
    }

    mask = calloc(T, 1);
    if (mask == NULL) {
	return E_ALLOC;
    }

    for (s=0; s<T; s++) {
	t = s + t1;
	if (na(yvar[t]) || (xvar != NULL && na(xvar[t]))) {
	    mask[s] = 1;
	} else {
	    n++;
	}
    }

    if (n == 0) {
	free(mask);
	return E_MISSDATA;
    }

    if (fit == PLOT_FIT_CUBIC) {
	k = 4;
    } else if (fit == PLOT_FIT_QUADRATIC) {
	k = 3;
    } else if (fit == PLOT_FIT_LOESS) {
	k = 1;
    } else {
	k = 2;
    }

    y = gretl_column_vector_alloc(n);
    X = gretl_matrix_alloc(n, k);
    if (y == NULL || X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    i = 0;
    for (s=0; s<T; s++) {
	t = s + t1;
	if (!mask[s]) {
	    j = 0;
	    if (fit == PLOT_FIT_LOGLIN) {
		y->val[i] = log(yvar[t]);
	    } else {
		y->val[i] = yvar[t];
	    }
	    if (fit != PLOT_FIT_LOESS) {
		gretl_matrix_set(X, i, j++, 1.0);
	    }
	    xt = (xvar != NULL)? xvar[t] : s;
	    if (fit == PLOT_FIT_INVERSE) {
		gretl_matrix_set(X, i, j++, 1.0 / xt);
	    } else if (fit == PLOT_FIT_LINLOG) {
		gretl_matrix_set(X, i, j++, log(xt));
	    } else {
		gretl_matrix_set(X, i, j++, xt);
	    }
	    if (fit == PLOT_FIT_QUADRATIC || fit == PLOT_FIT_CUBIC) {
		gretl_matrix_set(X, i, j++, xt * xt);
	    }
	    if (fit == PLOT_FIT_CUBIC) {
		gretl_matrix_set(X, i, j, xt * xt * xt);
	    }
	    i++;
	}
    }

 bailout:

    free(mask);

    if (err) {
	gretl_matrix_free(y);
	gretl_matrix_free(X);
    } else {
	*py = y;
	*pX = X;
    }

    return err;
}

static int skip_matrix_comment (FILE *fp, int *rows, int *cols,
				int *is_complex, int *err)
{
    int c, ret = 0;

    c = fgetc(fp);

    if (c == '#') {
	char *p, test[16] = {0};
	int n, i = 0;

	ret = 1;
	while (c != '\n' && c != EOF && c != -1) {
	    c = fgetc(fp);
	    if (i < 15) {
		test[i] = c;
	    }
	    i++;
	}
	if ((p = strstr(test, "rows:")) != NULL) {
	    if (sscanf(p + 5, "%d", &n) == 1) {
		*rows = n;
	    }
	} else if ((p = strstr(test, "columns:")) != NULL) {
	    if (sscanf(p + 8, "%d", &n) == 1) {
		*cols = n;
	    }
	} else if ((p = strstr(test, "complex:")) != NULL) {
	    if (sscanf(p + 8, "%d", &n) == 1) {
		*is_complex = (n == 1);
	    }
	}
    } else {
	ungetc(c, fp);
    }

    if (c == EOF || c == -1) {
	fprintf(stderr, "reached premature end of file\n");
	*err = E_DATA;
    }

    return ret;
}

#if G_BYTE_ORDER == G_BIG_ENDIAN

static void matrix_swap_endianness (gretl_matrix *A)
{
    int i, n = A->rows * A->cols;

    for (i=0; i<n; i++) {
	reverse_double(A->val[i]);
    }
}

#endif

static gretl_matrix *read_binary_matrix_file (FILE *fp, int *err)
{
    gretl_matrix *A = NULL;
    char header[20] = {0};
    int is_complex = 0;
    gint32 dim[2];

    if (fread(header, 1, 19, fp) < 19) {
	*err = E_DATA;
    } else if (fread(dim, sizeof *dim, 2, fp) < 2) {
	*err = E_DATA;
    } else if (dim[0] <= 0 || dim[1] <= 0) {
	*err = E_DATA;
    }

    if (!*err) {
	if (!strcmp(header, "gretl_binary_matrix")) {
	    ; /* OK */
	} else if (!strcmp(header, "gretl_binar_cmatrix")) {
	    is_complex = 1;
	} else {
	    *err = E_DATA;
	}
    }

    if (!*err) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
	reverse_int(dim[0]);
	reverse_int(dim[1]);
#endif
	A = gretl_matrix_alloc(dim[0], dim[1]);
	if (A == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	size_t n = dim[0] * dim[1];

	if (fread(A->val, sizeof *A->val, n, fp) < n) {
	    *err = E_DATA;
	    gretl_matrix_free(A);
	    A = NULL;
	} else if (is_complex) {
	    gretl_matrix_set_complex_full(A, 1);
	}
    }

#if G_BYTE_ORDER == G_BIG_ENDIAN
    if (A != NULL) {
	matrix_swap_endianness(A);
    }
#endif

    fclose(fp);

    return A;
}

#ifndef WIN32

/* In reading matrices, accept "NA" and "." (Stata) for NaN? */

static double unix_scan_NA (FILE *fp, int *err)
{
    char test[3] = {0};

    if (fscanf(fp, "%2s", test) < 1) {
	*err = E_DATA;
	return 0;
    } else if (!strcmp(test, "NA") || !strcmp(test, ".")) {
	return NADBL;
    } else {
	gretl_errmsg_sprintf(_("got invalid field '%s'"), test);
	*err = E_DATA;
	return 0;
    }
}

#endif

/**
 * gretl_matrix_read_from_file:
 * @fname: name of file.
 * @import: non-zero means we're importing via dotdir.
 * @err: location to receive error code.
 *
 * Reads a matrix from a file by the name @fname; if it's
 * a text file the column separator must be space or tab. It is
 * assumed that the dimensions of the matrix (number of rows and
 * columns) are found on the first line of the file, so no
 * heuristics are necessary. In case of error @err is filled
 * appropriately.
 *
 * Returns: The matrix as read from file, or NULL.
 */

gretl_matrix *gretl_matrix_read_from_file (const char *fname,
					   int import, int *err)
{
    char fullname[FILENAME_MAX];
    gchar *unzname = NULL;
    int r = -1, c = -1, n = 0;
    int gz, bin = 0;
    gretl_matrix *A = NULL;
    int is_complex = 0;
    FILE *fp = NULL;

    gz = has_suffix(fname, ".gz");

    if (!gz) {
	bin = has_suffix(fname, ".bin");
    }

    if (import) {
	gretl_build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
	strcpy(fullname, fname);
	gretl_maybe_prepend_dir(fullname);
    }

    if (gz) {
	char tmp[FILENAME_MAX];

	gretl_build_path(tmp, gretl_dotdir(), "_unztmp_.mat", NULL);
	*err = gretl_gunzip(fullname, tmp);
	if (!*err) {
	    fp = gretl_fopen(tmp, "rb");
	    unzname = g_strdup(tmp);
	}
    } else {
	fp = gretl_fopen(fullname, "rb");
    }

    if (fp == NULL) {
	*err = E_FOPEN;
	g_free(unzname);
	return NULL;
    }

    if (bin) {
	return read_binary_matrix_file(fp, err);
    }

    /* 'Eat' any leading comment lines starting with '#',
       but while we go, scan for Matlab-style rows/columns
       specification.
    */
    while (!*err && skip_matrix_comment(fp, &r, &c, &is_complex, err)) {
	;
    }

    if (!*err) {
	if (r >= 0 && c >= 0) {
	    /* we got dimensions from the preamble */
	    n = 2;
	} else {
	    n = fscanf(fp, "%d %ds\n", &r, &c);
	}
	if (n < 2 || r < 0 || c < 0) {
	    fprintf(stderr, "error reading rows, cols (r=%d, c=%d)\n",
		    r, c);
	    *err = E_DATA;
	} else {
	    A = gretl_matrix_alloc(r, c);
	    if (A == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    if (!*err) {
	long pos = 0;
	double x;
	int i, j;

	gretl_push_c_numeric_locale();

	for (i=0; i<r && !*err; i++) {
	    for (j=0; j<c && !*err; j++) {
		pos = ftell(fp);
		n = fscanf(fp, "%lf", &x);
		if (n != 1) {
		    fseek(fp, pos, SEEK_SET);
#ifdef WIN32
		    x = win32_fscan_nonfinite(fp, err);
#else
		    x = unix_scan_NA(fp, err);
#endif
		}
		if (*err) {
		    fprintf(stderr, "error reading row %d, column %d\n", i+1, j+1);
		} else {
		    gretl_matrix_set(A, i, j, x);
		}
	    }
	}

	gretl_pop_c_numeric_locale();
    }

    fclose(fp);

    if (unzname != NULL) {
	gretl_remove(unzname);
	g_free(unzname);
    }

    if (!*err && is_complex) {
	gretl_matrix_set_complex_full(A, 1);
    }

    if (*err && A != NULL) {
	gretl_matrix_free(A);
	A = NULL;
    }

    return A;
}

#ifdef WIN32

static void win32_xna_out (double x, char pad, PRN *prn)
{
    if (isnan(x) || na(x)) {
	pputs(prn, "nan");
    } else if (x < 0) {
	pputs(prn, "-inf");
    } else {
	pputs(prn, "inf");
    }
    pputc(prn, pad);
}

#endif

static int matrix_to_csv (const gretl_matrix *A, const char *fname)
{
    const char *na_str;
    char delim;
    FILE *fp;
    double aij;
    int i, j;

    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    na_str = get_csv_na_write_string();
    delim = get_data_export_delimiter();
    if (delim == '\t') {
	/* we can't allow this! */
	delim = ',';
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<A->rows; i++) {
	for (j=0; j<A->cols; j++) {
	    aij = gretl_matrix_get(A, i, j);
	    if (na(aij)) {
		fputs(na_str, fp);
	    } else {
		fprintf(fp, "%.17g", aij);
	    }
	    if (j < A->cols - 1) {
		fputc(delim, fp);
	    } else {
		fputc('\n', fp);
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    fclose(fp);

    return 0;
}

/**
 * gretl_matrix_write_to_file:
 * @A: matrix to write.
 * @fname: name of file.
 * @use_dotdir: non-zero means we're exporting via dotdir.
 *
 * Writes the matrix @A to a file by the name @fname; if it's a text
 * file, the column separator is tab. The number of rows and columns
 * are written on the first line of the file.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_matrix_write_to_file (gretl_matrix *A, const char *fname,
				int use_dotdir)
{
    char targ[FILENAME_MAX];
    int r, c, i, j;
    int is_complex = 0;
    PRN *prn = NULL;
    FILE *fp = NULL;
    int csv = 0;
    int gz, bin = 0;
    int err = 0;

    gz = has_suffix(fname, ".gz");

    if (!gz) {
	bin = has_suffix(fname, ".bin");
    }

    if (!gz && !bin) {
	csv = has_suffix(fname, ".csv");
	if (csv && A->is_complex) {
	    return E_CMPLX;
	}
    }

    if (use_dotdir) {
	gretl_build_path(targ, gretl_dotdir(), fname, NULL);
    } else {
	fname = gretl_maybe_switch_dir(fname);
	strcpy(targ, fname);
    }

    if (csv) {
	return matrix_to_csv(A, targ);
    }

    if (bin) {
	fp = gretl_fopen(targ, "wb");
    } else if (gz) {
	prn = gretl_gzip_print_new(targ, -1, &err);
    } else {
	prn = gretl_print_new_with_filename(targ, &err);
    }

    if (fp == NULL && prn == NULL) {
	return E_FOPEN;
    }

    if (A->is_complex) {
	/* reset to "pretend real" for writing */
	is_complex = 1;
	gretl_matrix_set_complex_full(A, 0);
    }

    r = A->rows;
    c = A->cols;

    if (bin) {
	const char *header = is_complex ? "gretl_binar_cmatrix" :
	    "gretl_binary_matrix";
	gint32 dim[2] = {r, c};
	size_t n = r * c;

#if G_BYTE_ORDER == G_BIG_ENDIAN
	double x;
	int k;

	fwrite(header, 1, strlen(header), fp);
	for (i=0; i<2; i++) {
	    k = dim[i];
	    reverse_int(k);
	    fwrite(&k, sizeof k, 1, fp);
	}
	for (i=0; i<n; i++) {
	    x = A->val[i];
	    reverse_double(x);
	    fwrite(&x, sizeof x, 1, fp);
	}
#else
	fwrite(header, 1, strlen(header), fp);
	fwrite(dim, sizeof *dim, 2, fp);
	fwrite(A->val, sizeof *A->val, n, fp);
#endif
	fclose(fp);
    } else {
	/* !bin: textual representation */
	int format_g = libset_get_bool(MWRITE_G);
	char pad, d = '\t';
	double x;

	if (is_complex) {
	    pprintf(prn, "# rows: %d\n", r);
	    pprintf(prn, "# columns: %d\n", c);
	    pprintf(prn, "# complex: %d\n", 1);
	} else {
	    pprintf(prn, "%d%c%d\n", r, d, c);
	}

	gretl_push_c_numeric_locale();

	for (i=0; i<r; i++) {
	    for (j=0; j<c; j++) {
		pad = (j == c-1)? '\n' : d;
		x = gretl_matrix_get(A, i, j);
#ifdef WIN32
		if (na(x)) {
		    win32_xna_out(x, pad, prn);
		    continue;
		}
#endif
		if (format_g) {
		    pprintf(prn, "%g", x);
		} else {
		    pprintf(prn, "%26.18E", x);
		}
		pputc(prn, pad);
	    }
	}

	gretl_pop_c_numeric_locale();
	gretl_print_destroy(prn);
    }

    if (is_complex) {
	/* reset A's original status */
	gretl_matrix_set_complex_full(A, 1);
    }

    return err;
}

static void make_numstr (char *s, double x)
{
    if (x == -0.0) {
	x = 0.0;
    }

    if (isnan(x)) {
	strcpy(s, "nan");
	return;
    }

    sprintf(s, "%#.5g", x);

    /* remove surplus, or add deficient, zeros */

    if (strstr(s, ".00000")) {
	s[strlen(s) - 1] = 0;
    } else {
	char *p = s;
	int n = 0;

	while (*p) {
	    if (isdigit(*p)) {
		n++;
	    } else if (isalpha(*p)) {
		n = 0;
		break;
	    }
	    p++;
	}
	if (n > 0 && n < 5) {
	    strcat(s, "0");
	}
    }
}

static int max_numchars (const gretl_matrix *m)
{
    char s[24];
    int i, n = m->rows * m->cols;
    int c, cmax = 0;

    for (i=0; i<n && cmax<6; i++) {
	sprintf(s, "%g", m->val[i]);
	c = strlen(s);
	if (c > cmax) {
	    cmax = c;
	}
    }

    return cmax;
}

static int max_label_length (const char **names, int n)
{
    int i, len, maxlen = 0;

    for (i=0; i<n; i++) {
	if (names[i] != NULL) {
	    len = g_utf8_strlen(names[i], -1);
	    if (len > maxlen) {
		maxlen = len;
	    }
	}
    }

    return maxlen;
}

static void
real_matrix_print_to_prn (const gretl_matrix *m,
			  const char *msg,
			  int plain,
			  const char **colnames,
			  const DATASET *dset,
			  int imin, int imax,
			  PRN *prn)
{
    const char **rownames = NULL;
    char numstr[32];
    char obs[OBSLEN];
    double x;
    int strwidth = 12;
    int rnamelen = 0;
    int dated = 0;
    int cmax = 0, cpad = 2;
    int mt1 = 0, mt2 = 0;
    int i, j;

    if (prn == NULL) {
	return;
    }

    if (m == NULL) {
	if (msg != NULL && *msg != '\0') {
	    pprintf(prn, _("%s: matrix is NULL"), msg);
	} else {
		pputs(prn, _("matrix is NULL"));
	}
	pputc(prn, '\n');
	return;
    } else if (m->rows == 0 || m->cols == 0) {
	if (msg != NULL && *msg != '\0') {
	    pprintf(prn, _("%s: matrix is empty (%d x %d)"),
		    msg, m->rows, m->cols);
	} else {
	    pprintf(prn, _("matrix is empty (%d x %d)"),
		    m->rows, m->cols);
	}
	pputc(prn, '\n');
	return;
    }

    dated = gretl_matrix_is_dated(m);
    if (dated) {
	mt1 = gretl_matrix_get_t1(m);
	mt2 = gretl_matrix_get_t2(m);
    }

    /* @plain means skip the header stuff */

    if (msg != NULL && *msg != '\0' && !plain) {
	pprintf(prn, "%s (%d x %d)", msg, m->rows, m->cols);
	if (imin > 0 && imax > 0) {
	    pprintf(prn, " rows %d to %d\n\n", imin+1, imax);
	} else if (dated && dset == NULL) {
	    pprintf(prn, " [t1 = %d, t2 = %d]\n\n", mt1 + 1, mt2 + 1);
	} else {
	    pputs(prn, "\n\n");
	}
    }

    cmax = max_numchars(m);
    if (cmax > 5) {
	cmax = 0;
    }

    if (colnames == NULL) {
	colnames = gretl_matrix_get_colnames(m);
    }

    rownames = gretl_matrix_get_rownames(m);

    if (rownames != NULL) {
	rnamelen = max_label_length(rownames, m->rows);
    } else if (dated && dset != NULL) {
	rnamelen = max_obs_marker_length(dset);
    }

    if (colnames != NULL) {
	int numeric_names = 1;

	if (rnamelen > 0) {
	    bufspace(rnamelen + 1, prn);
	}
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%*.*s ", strwidth, strwidth, colnames[j]);
	    if (!numeric_string(colnames[j])) {
		numeric_names = 0;
	    }
	}
	pputc(prn, '\n');
	if (numeric_names) {
	    pputc(prn, '\n');
	}
    }

    if (imin < 0) imin = 0;
    if (imax < 0) imax = m->rows;

    for (i=imin; i<imax; i++) {
	if (rnamelen > 0) {
	    if (rownames != NULL) {
		pprintf(prn, "%*s ", rnamelen, rownames[i]);
	    } else {
		ntolabel(obs, i + mt1, dset);
		pprintf(prn, "%*s ", rnamelen, obs);
	    }
	}
	for (j=0; j<m->cols; j++) {
	    x = gretl_matrix_get(m, i, j);
	    if (cmax > 0) {
		if (colnames != NULL) {
		    bufspace(strwidth - cmax - cpad, prn);
		}
		if (isnan(x)) {
		    strcpy(numstr, "nan");
		} else {
		    sprintf(numstr, "%g", x);
		}
		pprintf(prn, "%*s ", cmax + cpad, numstr);
	    } else {
		make_numstr(numstr, x);
		pprintf(prn, "%*s ", strwidth, numstr);
	    }
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n');
}

/**
 * gretl_matrix_print_to_prn:
 * @m: matrix to print.
 * @msg: accompanying message text (or NULL if no message is wanted).
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn.
 */

void
gretl_matrix_print_to_prn (const gretl_matrix *m,
			   const char *msg, PRN *prn)
{
    real_matrix_print_to_prn(m, msg, 0, NULL, NULL,
			     -1, -1, prn);
}

/**
 * gretl_matrix_print_range:
 * @m: matrix to print.
 * @rmin: first row to print.
 * @prn: last row to print.
 *
 * Prints the specified rows of matrix @m to @prn.
 */

void gretl_matrix_print_range (const gretl_matrix *m,
			       const char *msg,
			       int rmin, int rmax,
			       PRN *prn)
{
    if (m->is_complex) {
	gretl_cmatrix_print_range(m, msg, rmin, rmax, prn);
    } else {
	real_matrix_print_to_prn(m, msg, 0, NULL, NULL,
				 rmin, rmax, prn);
    }
}

/**
 * gretl_matrix_print_with_col_heads:
 * @m: T x k matrix to print.
 * @title: accompanying title (or NULL if no title is wanted).
 * @heads: array of k strings to identify the columns.
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn, with column headings given
 * by @heads.
 */

void gretl_matrix_print_with_col_heads (const gretl_matrix *m,
					const char *title,
					const char **heads,
					const DATASET *dset,
					PRN *prn)
{
    real_matrix_print_to_prn(m, title, 0, heads, dset,
			     -1, -1, prn);
}

static void maybe_print_col_heads (const gretl_matrix *m,
				   const char *fmt,
				   int wid, int prec,
				   int icast, int llen,
				   PRN *prn)
{
    const char **heads;

    heads = gretl_matrix_get_colnames(m);

    if (heads != NULL) {
	int numeric = 1;
	const char *s;
	char wtest[32];
	double x;
	int j, n, d;

	x = gretl_matrix_get(m, 0, 0);

	if (icast) {
	    if (wid >= 0 && prec >= 0) {
		snprintf(wtest, 32, fmt, wid, prec, (int) x);
	    } else if (wid >= 0 || prec >= 0) {
		n = (wid >= 0)? wid : prec;
		snprintf(wtest, 32, fmt, n, (int) x);
	    } else {
		snprintf(wtest, 32, fmt, (int) x);
	    }
	} else {
	    if (wid >= 0 && prec >= 0) {
		snprintf(wtest, 32, fmt, wid, prec, x);
	    } else if (wid >= 0 || prec >= 0) {
		n = (wid >= 0)? wid : prec;
		snprintf(wtest, 32, fmt, n, x);
	    } else {
		snprintf(wtest, 32, fmt, x);
	    }
	}

	if (llen > 0) {
	    bufspace(llen + 1, prn);
	}

	n = strlen(wtest);
	for (j=0; j<m->cols; j++) {
	    s = heads[j];
	    d = strlen(s) - g_utf8_strlen(s, -1);
	    pprintf(prn, "%*s", n + d, s);
	    if (!numeric_string(s)) {
		numeric = 0;
	    }
	}
	pputc(prn, '\n');
	if (numeric) {
	    pputc(prn, '\n');
	}
    }
}

/* Return an array of char holding 1s for columns in
   @m that are all integers, 0s otherwise. Except
   that if we find that either none or all the columns
   if @m are composed of integers we scratch the
   array and signal the result via the @intcast pointer.
*/

static char *get_intcols (const gretl_matrix *m,
			  int *intcast)
{
    char *ret = malloc(m->cols);

    if (ret != NULL) {
	double x;
	int i, j, n = 0;

	for (j=0; j<m->cols; j++) {
	    ret[j] = 1;
	    for (i=0; i<m->rows; i++) {
		x = gretl_matrix_get(m, i, j);
		if (x != floor(x)) {
		    ret[j] = 0;
		    break;
		}
	    }
	    n += ret[j];
	}
	if (n == 0 || n == m->cols) {
	    /* no, or all, integer columns */
	    free(ret);
	    ret = NULL;
	    *intcast = n > 0;
	}
    }

    return ret;
}

/* Take a format which is designed for floating-point
   values and convert it for use with integers: preserve
   the width but discard the precision, and change the
   trailing 'f', 'g' or whatever to 'd'.
*/

static void revise_integer_format (char *ifmt)
{
    char *p = strchr(ifmt, '.');

    if (p != NULL) {
	*p = '\0';
	strcat(p, "d");
    } else {
	ifmt[strlen(ifmt)-1] = 'd';
    }
}

/**
 * gretl_matrix_print_with_format:
 * @m: matrix to print.
 * @fmt: a printf-type format string.
 * @wid: an integer width, or 0.
 * @prec: an integer precision, or 0.
 * @prn: pointer to gretl printing struct.
 *
 * Prints the matrix @m to @prn, with the elements of @m
 * formatted as specified by @fmt, @wid and @prec.
 * The arguments @wid and/or @prec are required only if
 * @fmt contains one or more placeholders ("*") to be filled out.
 * For example, if @fmt is "\%14.6f" then neither @wid nor
 * @prec is needed and both should be set to 0; if
 * @fmt is "\%*.6f" then a @wid value is needed but @prec
 * should be 0.
 */

void gretl_matrix_print_with_format (const gretl_matrix *m,
				     const char *fmt,
				     int wid, int prec,
				     PRN *prn)
{
    if (prn == NULL) {
	return;
    }

    if (gretl_is_null_matrix(m) || fmt == NULL || *fmt == '\0') {
	real_matrix_print_to_prn(m, NULL, 1, NULL, NULL, -1, -1, prn);
    } else if (m->is_complex) {
	gretl_cmatrix_printf(m, fmt, prn);
    } else {
	const char **rownames = NULL;
	char *ifmt = NULL, *xfmt = NULL;
	char *intcols = NULL;
	int llen = 0, intcast = 0;
	int cpos = strlen(fmt) - 1;
	int variant = 0;
	char c = fmt[cpos];
	int i, j, k;
	double x;

	xfmt = gretl_strdup(fmt);

	if (c == 'd' || c == 'u' || c == 'x' || c == 'l') {
	    intcast = 1;
	} else if (fmt[cpos-1] == 'v') {
	    /* new-style "variant" column format */
	    intcols = get_intcols(m, &intcast);
	    gretl_delchar('v', xfmt);
	    variant = 1;
	} else if (c == 'v') {
	    /* old-style "variant" column format */
	    intcols = get_intcols(m, &intcast);
	    xfmt[cpos] = 'g';
	    variant = 1;
	}

	if (intcast || intcols) {
	    ifmt = gretl_strdup(xfmt);
	    if (variant) {
		revise_integer_format(ifmt);
	    }
	}

	rownames = gretl_matrix_get_rownames(m);
	if (rownames != NULL) {
	    llen = max_label_length(rownames, m->rows);
	}

	maybe_print_col_heads(m, xfmt, wid, prec, intcast, llen, prn);

	for (i=0; i<m->rows; i++) {
	    if (rownames != NULL) {
		pprintf(prn, "%*s ", llen, rownames[i]);
	    }
	    for (j=0; j<m->cols; j++) {
		x = gretl_matrix_get(m, i, j);
		if (intcols != NULL) {
		    intcast = intcols[j];
		}
		if (intcast) {
		    if (wid > 0 && prec > 0) {
			pprintf(prn, ifmt, wid, prec, (int) x);
		    } else if (wid > 0 || prec > 0) {
			k = (wid > 0)? wid : prec;
			pprintf(prn, ifmt, k, (int) x);
		    } else {
			pprintf(prn, ifmt, (int) x);
		    }
		} else {
		    if (wid > 0 && prec > 0) {
			pprintf(prn, xfmt, wid, prec, x);
		    } else if (wid > 0 || prec > 0) {
			k = (wid > 0)? wid : prec;
			pprintf(prn, xfmt, k, x);
		    } else {
			pprintf(prn, xfmt, x);
		    }
		}
	    }
	    pputc(prn, '\n');
	}

	free(intcols);
	free(ifmt);
	free(xfmt);
    }
}

/* gretl_matrix_write_constructor() is intended for internal use in
   gretl_bundle.c, for creating a string representation of a matrix
   contained in a bundle when the bundle is passed as argument to
   gretl_bundle_write_constructor(). In that context it will not be
   appropriate to show the constructor for a large matrix in extenso,
   hence the size limitation here. However it's not clear what exactly
   we should do with an oversize matrix.
*/

gchar *gretl_matrix_write_constructor (const gretl_matrix *m)
{
    gchar *ret = NULL;
    int n;

    if (gretl_is_null_matrix(m)) {
	return g_strdup("{}");
    }

    n = m->rows * m->cols;

    if (n > 16) {
	/* should we do something else here? */
	return g_strdup("{}");
    } else {
	GString *gs;
	gsize sz = 4 * n;
	double mij;
	int i, j;

	gs = g_string_sized_new(sz);
	g_string_append_c(gs, '{');
	for (i=0; i<m->rows; i++) {
	    for (j=0; j<m->cols; j++) {
		mij = gretl_matrix_get(m, i, j);
		if (na(mij)) {
		    g_string_append(gs, "NA");
		} else {
		    g_string_append_printf(gs, "%g", mij);
		}
		if (j < m->cols - 1) {
		    g_string_append_c(gs, ',');
		} else if (i < m->rows - 1) {
		    g_string_append_c(gs, ';');
		}
	    }
	}
	g_string_append_c(gs, '}');
	ret = g_string_free(gs, FALSE);
    }

    return ret;
}

/**
 * debug_print_matrix:
 * @m: matrix to print.
 * @msg: accompanying message text (or NULL if no message is wanted).
 *
 * Prints the matrix @m to stderr with high precision, along
 * with the address of the matrix struct.
 */

void debug_print_matrix (const gretl_matrix *m, const char *msg)
{
    char full[64] = {0};

    if (msg != NULL) {
	strncpy(full, msg, 32);
	sprintf(full + strlen(full), " (%p)", (void *) m);
    } else {
	sprintf(full, " (%p)", (void *) m);
    }

    if (m != NULL) {
	int i, n = m->rows * m->cols;
	int d = (int) ceil(log10((double) n));

	fprintf(stderr, "%s\n", full);
	for (i=0; i<n; i++) {
	    fprintf(stderr, "val[%0*d] = % .10E\n", d, i, m->val[i]);
	}
    } else {
	PRN *prn;
	int err = 0;

	prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
	if (!err) {
	    gretl_matrix_print_to_prn(m, full, prn);
	    gretl_print_destroy(prn);
	}
    }
}

static int count_unmasked_elements (const char *mask, int n)
{
    int i, c = 0;

    for (i=0; i<n; i++) {
	c += (mask[i] == 0);
    }

    return c;
}

/**
 * gretl_matrix_cut_rows:
 * @m: matrix to process.
 * @mask: character array of length equal to the rows of @m,
 * with 1s indicating rows to be cut, 0s for rows to be
 * retained.
 *
 * In-place reduction of @m based on @mask: the masked rows
 * are cut out of @m.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_rows (gretl_matrix *m, const char *mask)
{
    int i, j, k, n;
    double x;

    if (m == NULL || mask == NULL) {
	return E_DATA;
    }

    n = count_unmasked_elements(mask, m->rows);

    for (j=0; j<m->cols; j++) {
	k = 0;
	for (i=0; i<m->rows; i++) {
	    if (!mask[i]) {
		x = gretl_matrix_get(m, i, j);
		m->val[j * n + k] = x;
		k++;
	    }
	}
    }

    m->rows = n;

    return 0;
}

/**
 * gretl_matrix_cut_cols:
 * @m: matrix to process.
 * @mask: character array of length equal to the cols of @m,
 * with 1s indicating cols to be cut, 0s for cols to be
 * retained.
 *
 * In-place reduction of @m based on @mask: the masked cols
 * are cut out of @m.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_cols (gretl_matrix *m, const char *mask)
{
    int i, j, k, n;
    double x;

    if (m == NULL || mask == NULL) {
	return E_DATA;
    }

    n = count_unmasked_elements(mask, m->cols);

    k = 0;
    for (j=0; j<m->cols; j++) {
	if (!mask[j]) {
	    for (i=0; i<m->rows; i++) {
		x = gretl_matrix_get(m, i, j);
		m->val[k++] = x;
	    }
	}
    }

    m->cols = n;

    return 0;
}

/**
 * gretl_matrix_cut_rows_cols:
 * @m: square matrix to process.
 * @mask: character array of length equal to the dimension
 * of @m, with 1s indicating rows and columns to be cut, 0s
 * for rows/columns to be retained.
 *
 * In-place reduction of @m based on @mask: the masked rows
 * and columns are cut out of @m. (The data array within @m
 * is not "physically" resized.)
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_rows_cols (gretl_matrix *m, const char *mask)
{
    gretl_matrix *tmp;
    double x;
    int i, j, k, l, n;

    if (m == NULL || mask == NULL) {
	return E_DATA;
    } else if (m->rows != m->cols) {
	return E_NONCONF;
    }

    n = count_unmasked_elements(mask, m->rows);
    if (n == 0) {
	gretl_matrix_reuse(m, 0, 0);
	return 0;
    }

    /* create smaller temporary matrix */
    tmp = gretl_matrix_alloc(n, n);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    /* copy unmasked values into temp matrix */
    k = 0;
    for (i=0; i<m->rows; i++) {
	if (!mask[i]) {
	    l = 0;
	    for (j=0; j<m->cols; j++) {
		if (!mask[j]) {
		    x = gretl_matrix_get(m, i, j);
		    gretl_matrix_set(tmp, k, l++, x);
		}
	    }
	    k++;
	}
    }

    /* redimension the original matrix, copy the values back,
       and free the temp matrix */
    gretl_matrix_reuse(m, n, n);
    gretl_matrix_copy_values(m, tmp);
    gretl_matrix_free(tmp);

    return 0;
}

/**
 * gretl_matrix_zero_row_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Checks matrix @m for rows that are all zero.  If there are
 * any such rows, constructs a mask of length equal to the
 * number of rows in @m, with 1s indicating zero rows, 0s
 * elsewhere.  If there are no such rows, returns NULL.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char *gretl_matrix_zero_row_mask (const gretl_matrix *m, int *err)
{
    char *mask = NULL;
    int row0, any0 = 0;
    int i, j;

    mask = calloc(m->rows, 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<m->rows; i++) {
	row0 = 1;
	for (j=0; j<m->cols; j++) {
	    if (gretl_matrix_get(m, i, j) != 0.0) {
		row0 = 0;
		break;
	    }
	}
	if (row0) {
	    mask[i] = 1;
	    any0 = 1;
	}
    }

    if (!any0) {
	free(mask);
	mask = NULL;
    }

    return mask;
}

/**
 * gretl_matrix_zero_col_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Checks matrix @m for columns that are all zero.  If there are
 * any such columns, constructs a mask of length equal to the
 * number of columns in @m, with 1s indicating zero columns, 0s
 * elsewhere.  If there are no such columns, returns NULL.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char *gretl_matrix_zero_col_mask (const gretl_matrix *m, int *err)
{
    char *mask = NULL;
    int col0, any0 = 0;
    int i, j;

    mask = calloc(m->cols, 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (j=0; j<m->cols; j++) {
	col0 = 1;
	for (i=0; i<m->rows; i++) {
	    if (gretl_matrix_get(m, i, j) != 0.0) {
		col0 = 0;
		break;
	    }
	}
	if (col0) {
	    mask[j] = 1;
	    any0 = 1;
	}
    }

    if (!any0) {
	free(mask);
	mask = NULL;
    }

    return mask;
}

/**
 * gretl_matrix_zero_diag_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Checks square matrix @m for diagonal elements that are zero.
 * If there are any such, constructs a mask of length equal to
 * the number of rows (and columns) of @m, with 1s indicating
 * the zero diagonal entries. If there are no zero diagonal
 * elements, returns NULL.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char *gretl_matrix_zero_diag_mask (const gretl_matrix *m, int *err)
{
    char *mask = NULL;
    int i, trim = 0;

    if (gretl_is_null_matrix(m)) {
	return NULL;
    }

    if (m->rows != m->cols) {
	*err = E_NONCONF;
	return NULL;
    }

    for (i=0; i<m->rows; i++) {
	if (gretl_matrix_get(m, i, i) == 0.0) {
	    trim = 1;
	    break;
	}
    }

    if (trim) {
	mask = calloc(m->rows, 1);
	if (mask == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=0; i<m->rows; i++) {
		if (gretl_matrix_get(m, i, i) == 0.0) {
		    mask[i] = 1;
		}
	    }
	}
    }

    return mask;
}

/**
 * gretl_matrix_rank_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Performs a QR decomposition of matrix @m and uses this
 * to assess the rank of @m.  If @m is not of full rank,
 * constructs a mask of length equal to the numbers of
 * columns in @m, with 1s in positions corresponding
 * to diagonal elements of R that are effectively 0, and
 * 0s elsewhere.  If @m is of full column rank, NULL is
 * returned.
 *
 * E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or NULL.
 */

char *gretl_matrix_rank_mask (const gretl_matrix *m, int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    char *mask = NULL;
    double Ri;
    int fullrank = 1;
    int i, n = m->cols;

    Q = gretl_matrix_copy(m);
    if (Q == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    R = gretl_matrix_alloc(n, n);
    if (R == NULL) {
	gretl_matrix_free(Q);
	*err = E_ALLOC;
	return NULL;
    }

    *err = gretl_matrix_QR_decomp(Q, R);

    if (!*err) {
	mask = calloc(n, 1);
	if (mask == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	for (i=0; i<n; i++) {
	    Ri = gretl_matrix_get(R, i, i);
	    if (fabs(Ri) < R_DIAG_MIN) {
		mask[i] = 1;
		fullrank = 0;
	    }
	}
    }

    if (*err || fullrank) {
	free(mask);
	mask = NULL;
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(R);

    return mask;
}

/**
 * gretl_matrix_mp_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or NULL if this is not needed.
 * @uhat: vector to hold the regression residuals, or NULL if
 * these are not needed.
 * @s2: pointer to receive residual variance, or NULL.  Note:
 * if @s2 is NULL, the "vcv" estimate will be plain (X'X)^{-1}.
 *
 * Computes OLS estimates using Cholesky factorization, via
 * the GMP multiple-precision library, and puts the
 * coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv and the residuals in @uhat.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_mp_ols (const gretl_vector *y, const gretl_matrix *X,
			 gretl_vector *b, gretl_matrix *vcv,
			 gretl_vector *uhat, double *s2)
{
    int (*matrix_mp_ols) (const gretl_vector *,
			  const gretl_matrix *,
			  gretl_vector *,
			  gretl_matrix *,
			  gretl_vector *,
			  double *);

    matrix_mp_ols = get_plugin_function("matrix_mp_ols");
    if (matrix_mp_ols == NULL) {
	return 1;
    }

    return (*matrix_mp_ols)(y, X, b, vcv, uhat, s2);
}

/* following: quadrature sources based on John Burkhart's GPL'd
   code at http://people.sc.fsu.edu/~jburkardt
*/

#define sign(x) (x < 0.0 ? -1.0 : 1.0)

/* Diagonalize a Jacobi (symmetric tridiagonal) matrix. On entry, @d
   contains the diagonal and @e the sub-diagonal. On exit, @d contains
   quadrature nodes or abscissae and @z the square roots of the
   corresponding weights; @e is used as workspace.

   This is a C adaptation of subroutine IMTQLX, by Sylvan Elhay,
   Jaroslav Kautsky and John Burkardt. See
   http://people.sc.fsu.edu/~jburkardt/f_src/quadrature_weights/
   specifically, qw_golub_welsch.f90
*/

static int diag_jacobi (int n, double *d, double *e, double *z)
{
    double b, c, f, g;
    double p, r, s;
    int j, k, l, m = 0;
    int maxiter = 30;
    int i, ii, mml;
    int err = 0;

    if (n == 1) {
	return 0;
    }

    e[n-1] = 0.0;

    errno = 0;

    for (l=1; l<=n; l++) {
	j = 0;
	for ( ; ; ) {
	    for (m=l; m<=n; m++) {
		if (m == n) {
		    break;
		}
		p = fabs(d[m-1]) + fabs(d[m]);
		if (fabs(e[m-1]) <= DBL_EPSILON * p) {
		    break;
		}
	    }
	    p = d[l-1];
	    if (m == l) {
		break;
	    }
	    if (j >= maxiter) {
		fprintf(stderr, "diag_jacobi: iteration limit exceeded\n");
		return E_NOCONV;
	    }
	    j++;
	    g = (d[l] - p) / (2.0 * e[l-1]);
	    r =  sqrt(g * g + 1.0);
	    g = d[m-1] - p + e[l-1] / (g + fabs(r) * sign(g));
	    s = 1.0;
	    c = 1.0;
	    p = 0.0;
	    mml = m - l;

	    for (ii=1; ii<=mml; ii++) {
		i = m - ii;
		f = s * e[i-1];
		b = c * e[i-1];
		if (fabs(g) <= fabs(f)) {
		    c = g / f;
		    r =  sqrt(c * c + 1.0);
		    e[i] = f * r;
		    s = 1.0 / r;
		    c = c * s;
		} else {
		    s = f / g;
		    r =  sqrt(s * s + 1.0);
		    e[i] = g * r;
		    c = 1.0 / r;
		    s = s * c;
		}
		g = d[i] - p;
		r = (d[i-1] - g) * s + 2.0 * c * b;
		p = s * r;
		d[i] = g + p;
		g = c * r - b;
		f = z[i];
		z[i] = s * z[i-1] + c * f;
		z[i-1] = c * z[i-1] - s * f;
	    }

	    d[l-1] = d[l-1] - p;
	    e[l-1] = g;
	    e[m-1] = 0.0;
	}
    }

    /* sorting */

    for (ii=2; ii<=m; ii++) {
	i = ii - 1;
	k = i;
	p = d[i-1];
	for (j=ii; j<=n; j++) {
	    if (d[j-1] < p) {
		k = j;
		p = d[j-1];
	    }
	}
	if (k != i) {
	    d[k-1] = d[i-1];
	    d[i-1] = p;
	    p = z[i-1];
	    z[i-1] = z[k-1];
	    z[k-1] = p;
	}
    }

    if (errno) {
	/* some calculation went badly */
	err = E_NAN;
	errno = 0;
    }

    return err;
}

/* Construct the Jacobi matrix for a quadrature rule: @d
   holds the diagonal, @e the subdiagonal.
*/

static double make_jacobi (int n, int kind, double *d, double *e)
{
    double z0 = 0.0;
    int i;

    if (kind == QUAD_LEGENDRE) {
	double xi, xj;

	z0 = 2.0;
	for (i=1; i<=n; i++) {
	    xi = i;
	    xj = 2.0 * i;
	    e[i-1] = sqrt(xi * xi / (xj * xj - 1.0));
	}
    } else if (kind == QUAD_LAGUERRE) {
	z0 = 1.0; /* tgamma(1.0) */
	for (i=1; i<=n; i++) {
	    d[i-1] = 2.0 * i - 1.0;
	    e[i-1] = i;
	}
    } else if (kind == QUAD_GHERMITE) {
	z0 = tgamma(0.5);
	for (i=1; i<=n; i++) {
	    e[i-1] = sqrt(i / 2.0);
	}
    }

    return z0;
}

/* compute basic Gauss quadrature formula */

static int gauss_quad_basic (int n, int kind, double *x, double *w)
{
    double z0, *tmp;
    int i, err;

    tmp = malloc(n * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    z0 = make_jacobi(n, kind, x, tmp);
    w[0] = sqrt(z0);
    err = diag_jacobi(n, x, tmp, w);

    if (!err) {
	for (i=0; i<n; i++) {
	    w[i] = w[i] * w[i];
	}
    }

    free(tmp);

    return err;
}

/* Scale Legendre quadrature rule to the interval (a, b) */

static int legendre_scale (int n, double *x, double *w,
			   double a, double b)
{
    double shft, slp;
    int i;

    if (na(a) || na(b)) {
	return E_INVARG;
    }

    if (fabs(b - a) <= DBL_EPSILON) {
	fprintf(stderr, "legendre: |b - a| too small\n");
	return E_DATA;
    }

    shft = (a + b) / 2.0;
    slp = (b - a) / 2.0;

    for (i=0; i<n; i++) {
	x[i] = shft + slp * x[i];
	w[i] = w[i] * slp;
    }

    return 0;
}

static int hermite_scale (int n, double *x, double *w,
			  double mu, double sigma)
{
    double rtpi = sqrt(M_PI);
    int i;

    if (na(mu) || na(sigma) || sigma <= 0.0) {
	return E_INVARG;
    }

    for (i=0; i<n; i++) {
	x[i] = mu + M_SQRT2 * sigma * x[i];
	w[i] /= rtpi;
    }

    return 0;
}

/**
 * gretl_quadrule_matrix_new:
 * @n: the order (i.e. the number of abscissae and weights).
 * @method: should be one of the values in #QuadMethod.
 * @a: optional parameter (see below).
 * @b: optional parameter (see below).
 * @err: location to receive error code.
 *
 * Calculates a quadrature "rule" (i.e. a set of abscissae or
 * nodes and associated weights) for use in numerical integration.
 * The three supported methods are Gauss-Hermite, Gauss-Legendre
 * and Gauss-Laguerre.
 *
 * If the optional parameters are not to be used, both should
 * be given the value NADBL.
 *
 * If the method is QUAD_GHERMITE, the default weights are W(x)
 * = exp(-x^2), but @a and @b can be used to apply a normal
 * scaling with mean @a and standard deviation @b. The limits
 * of integration are in principle minus infinity and plus
 * infinity.
 *
 * If the method is QUAD_LEGENDRE the weights are W(x) = 1 and
 * the default limits of integration are -1 and 1, but @a and @b
 * can be used to adjust the lower and upper limit respectively.
 *
 * In the QUAD_LAGUERRE case W(x) = exp(-x) and the limits
 * of integration are 0 and plus infinity; @a and @b are
 * ignored.
 *
 * Returns: an @n x 2 matrix containing the abscissae in the first
 * column and the weights in the second, or NULL on failure.
 */

gretl_matrix *gretl_quadrule_matrix_new (int n, int method,
					 double a, double b,
					 int *err)
{
    gretl_matrix *m;

    if (method < QUAD_GHERMITE || method >= QUAD_INVALID) {
	*err = E_DATA;
	return NULL;
    }

    if (n < 0) {
	*err = E_DATA;
	return NULL;
    } else if (n == 0) {
	return gretl_null_matrix_new();
    }

    m = gretl_zero_matrix_new(n, 2);

    if (m == NULL) {
	*err = E_ALLOC;
    } else {
	double *x = m->val;
	double *w = x + n;

	*err = gauss_quad_basic(n, method, x, w);

	if (!*err) {
	    if (method == QUAD_LEGENDRE) {
		if (na(a) && na(b)) {
		    ; /* default */
		} else if (a == -1.0 && b == 1.0) {
		    ; /* default */
		} else {
		    *err = legendre_scale(n, x, w, a, b);
		}
	    } else if (method == QUAD_GHERMITE) {
		if (!na(a) || !na(b)) {
		    *err = hermite_scale(n, x, w, a, b);
		}
	    } else if (method == QUAD_LAGUERRE) {
		; /* a and b are ignored */
	    }
	}
    }

    if (*err && m != NULL) {
	gretl_matrix_free(m);
	m = NULL;
    }

    return m;
}

gretl_matrix *gretl_gauss_hermite_matrix_new (int n, int *err)
{
    return gretl_quadrule_matrix_new(n, QUAD_GHERMITE, 0, 1, err);
}

/**
 * gretl_matrix_2d_convolution:
 * @A: first matrix.
 * @B: second matrix.
 * @err: location to receive error code.
 *
 * Computes the 2-dimensional convolution of the matrices
 * @A and @B. The code borrows from Octave's conv2.m.
 *
 * Returns: newly allocated convolution matrix, or %NULL on
 * failure.
 */

gretl_matrix *gretl_matrix_2d_convolution (const gretl_matrix *A,
					   const gretl_matrix *B,
					   int *err)
{
    gretl_matrix *C;
    double sum;
    int ci, cj, bj, aj, bi, ai;
#if 0
    /* Octave: "Convolution works fastest if we choose the A matrix to be
     *  the largest" -- but we're not really seeing any difference.
     */
    int Am, An, B, Bn, Cm, Cn;
    int edgM, edgN;

    if (A->rows * A->cols < B->rows * B->cols) {
	gretl_matrix *tmp = A;
	B = A;
	A = tmp;
    }

    Am = A->rows;
    An = A->cols;
    Bm = B->rows;
    Bn = B->cols;
    Cm = Am + Bm - 1;
    Cn = An + Bn - 1;
    edgM = Bm - 1;
    edgN = Bn - 1;
#else
    /* leave the input matrices alone */
    int Am = A->rows;
    int An = A->cols;
    int Bm = B->rows;
    int Bn = B->cols;
    int Cm = Am + Bm - 1;
    int Cn = An + Bn - 1;
    int edgM = Bm - 1;
    int edgN = Bn - 1;
#endif

    C = gretl_matrix_alloc(Cm, Cn);

    for (ci=0; ci < Cm; ci++) {
	for (cj=0; cj < Cn; cj++) {
            sum = 0;
            for (bj = Bn - 1 - MAX(0, edgN-cj), aj = MAX(0, cj-edgN);
		 bj >= 0 && aj < An; bj--, aj++) {
		bi = Bm - 1 - MAX(0, edgM-ci);
		ai = MAX(0, ci-edgM);
		for ( ; bi >= 0 && ai < Am; bi--, ai++) {
		    sum += gretl_matrix_get(A, ai, aj) * gretl_matrix_get(B, bi, bj);
		}
            }
            gretl_matrix_set(C, ci, cj, sum);
	}
    }

    return C;
}

/* scan format string: should be "%<num>m" or just "%m" */

static int check_matrix_scan_format (const char *fmt, int ns, int *n)
{
    int err = 0;

    if (*fmt == '%') {
	char *p;
	long len = strtol(++fmt, &p, 10);

	if (p == fmt) {
	    /* no "max rows" given */
	    *n = ns;
	} else if (len < 0) {
	    err = E_INVARG;
	} else if (len > ns) {
	    *n = ns;
	} else {
	    *n = len;
	}
	if (!err && strcmp(p, "m")) {
	    err = E_INVARG;
	}
    } else {
	err = E_INVARG;
    }

    return err;
}

gretl_matrix *vector_from_strings (char **S, int ns,
				   const char *fmt,
				   int *nvals,
				   int *err)
{
    gretl_matrix *m = NULL;
    int n = 0;

    *nvals = 0;
    *err = check_matrix_scan_format(fmt, ns, &n);

    if (!*err) {
	if (n == 0) {
	    m = gretl_null_matrix_new();
	} else {
	    m = gretl_column_vector_alloc(n);
	}
	if (m == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err && n > 0) {
	char *s, *p;
	double x;
	int i;

	gretl_push_c_numeric_locale();
	for (i=0; i<n; i++) {
	    s = S[i];
	    if (s == NULL || *s == '\0') {
		m->val[i] = NADBL;
	    } else {
		x = strtod(s, &p);
		if (*p == '\0') {
		    m->val[i] = x;
		    *nvals += 1;
		} else {
		    m->val[i] = NADBL;
		}
	    }
	}
	gretl_pop_c_numeric_locale();
    }

    return m;
}
