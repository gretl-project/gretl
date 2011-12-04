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
#include "matrix_extra.h"

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
 */

gretl_vector *gretl_vector_from_series (const double *x, 
					int t1, int t2)
{
    gretl_matrix *v = NULL;
    int T = t2 - t1 + 1;

    if (T > 0) {
	v = gretl_column_vector_alloc(T);
	if (v != NULL) {
	    memcpy(v->val, x + t1, T * sizeof *x);
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
 * If @x = #NADBL the matrix value is #M_NA.
 */

gretl_matrix *gretl_matrix_from_scalar (double x) 
{
    gretl_matrix *m = gretl_matrix_alloc(1, 1);

    if (na(x)) {
	m->val[0] = M_NA;
    } else {
	m->val[0] = x;
    }

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
    int t, i, j;
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
	for (j=i; j<k; j++) {
	    vv = 0.0;
	    n = 0;
	    for (t=dset->t1; t<=dset->t2; t++) {
		x = dset->Z[list[i+1]][t];
		y = dset->Z[list[j+1]][t];
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
		vv /= (n - 1); /* plain n? */
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

static gretl_matrix *
real_gretl_matrix_data_subset (const int *list, 
			       const DATASET *dset,
			       int t1, int t2, const char *mask,
			       int op, int *err)
{
    gretl_matrix *M;
    double x;
    int T, Tmax = t2 - t1 + 1;
    int k = list[0];
    int skip;
    int j, s, t;

    if (k <= 0 || Tmax <= 0) {
	*err = E_DATA;
	return NULL;
    }	

    *err = 0;

    T = Tmax;

    if (mask != NULL) {
	T -= get_mask_count(mask, Tmax);
    } else if (op == M_MISSING_TRIM) {
	*err = list_adjust_sample(list, &t1, &t2, dset);
	if (*err) {
	    return NULL;
	} else {
	    Tmax = T = t2 - t1 + 1;
	}
    } else if (op == M_MISSING_SKIP || op == M_MISSING_ERROR) {
	for (t=t1; t<=t2 && !*err; t++) {
	    for (j=0; j<k; j++) {
		x = dset->Z[list[j+1]][t];
		if (na(x)) {
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

    if (T <= 0) {
	*err = E_DATA;
	return NULL;
    }

    M = gretl_matrix_alloc(T, k);
    if (M == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    s = 0;
    for (t=t1; t<=t2; t++) {
	skip = 0;
	if (mask != NULL) {
	    skip = mask[t - t1];
	} else if (op == M_MISSING_SKIP) {
	    for (j=0; j<k; j++) {
		x = dset->Z[list[j+1]][t];
		if (na(x)) {
		    skip = 1;
		    break;
		}
	    }
	}
	if (!skip) {
	    for (j=0; j<k; j++) {
		x = dset->Z[list[j+1]][t];
		gretl_matrix_set(M, s, j, x);
	    }
	    if (s == 0) {
		gretl_matrix_set_t1(M, t);
	    } else if (s == T - 1) {
		gretl_matrix_set_t2(M, t);
		break;
	    }
	    s++;
	}
    }

    if (!*err && (mask != NULL || op == M_MISSING_OK)) {
	for (j=0; j<k && !*err; j++) {
	    for (t=0; t<T && !*err; t++) {
		x = gretl_matrix_get(M, t, j);
		if (na(x)) {
		    if (op == M_MISSING_OK) {
			/* convert NA to NAN */
			gretl_matrix_set(M, t, j, M_NA);
		    } else {
			*err = E_MISSDATA;
		    }
		}
	    }
	}
    }	

    if (*err) {
	gretl_matrix_free(M);
	M = NULL;
    } else if (T == Tmax) {
	gretl_matrix_set_t1(M, t1);
	gretl_matrix_set_t2(M, t2);
    }

    return M;
}

/**
 * gretl_matrix_data_subset:
 * @list: list of variable to process.
 * @dset: dataset struct.
 * @t1: starting observation.
 * @t2: ending observation.
 * @missop: how to handle missing observations.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @Z specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns. The @missop flag
 * can be %M_MISSING_OK to indicate that it's OK to include
 * missing values in the matrix (these become NaNs),
 * %M_MISSING_ERROR (it's an error of any missing values are
 * found), or %M_MISSING_SKIP (observations with any missing
 * values are omitted from the matrix).
 *
 * Returns: allocated matrix or NULL on failure. 
 */

gretl_matrix *gretl_matrix_data_subset (const int *list, 
					const DATASET *dset,
					int t1, int t2, int missop, 
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
 * @Z specified by @list, over the sample range @t1 to @t2,
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
				 int t1, int t2, const char *mask,
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
 * @Z specified by @list, using the observations (data rows)
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
	    double xti;
	    int i, s, t;

	    for (i=0; i<k; i++) {
		xi = dset->Z[list[i+1]];
		s = 0;
		for (t=dset->t1; t<=dset->t2; t++) {
		    if (mmask->val[t] != 0) {
			xti = xi[t];
			if (na(xti)) {
			    xti = M_NA;
			}
			if (s == 0) {
			    gretl_matrix_set_t1(X, t);
			} else if (s == T - 1) {
			    gretl_matrix_set_t2(X, t);
			}
			gretl_matrix_set(X, s++, i, xti);
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

/**
 * gretl_dataset_from_matrix:
 * @m: source matrix.
 * @list: list of columns (1-based) to include, or NULL.
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
				    int *err)
{
    DATASET *dset = NULL;
    const char **names;
    double x;
    int i, t, col, nv, T;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    T = gretl_matrix_rows(m);
    nv = gretl_matrix_cols(m);

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    col = list[i];
	    if (col < 1 || col > nv) {
		fprintf(stderr, "column %d: out of bounds\n", col);
		*err = E_DATA;
		break;
	    } 
	}
	nv = list[0];
    } 

    if (!*err) {
	dset = create_auxiliary_dataset(nv + 1, T);
	if (dset == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	return NULL;
    }

    names = gretl_matrix_get_colnames(m);

    for (i=1; i<=nv; i++) {
	col = (list != NULL)? list[i] - 1 : i - 1;
	for (t=0; t<T; t++) {
	    x = gretl_matrix_get(m, t, col);
	    dset->Z[i][t] = x;
	}
	if (names != NULL) {
	    strcpy(dset->varname[i], names[col]);
	} else {
	    sprintf(dset->varname[i], "col%d", col + 1);
	}
    }

    return dset;
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

    if (fit == PLOT_FIT_LOGLIN && !gretl_ispositive(t1, t2, yvar)) {
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

/**
 * gretl_matrix_read_from_text:
 * @fname: name of text file.
 * @err: location to receive error code.
 *
 * Reads a matrix from a text file by the name @fname; the column
 * separator must be space or tab. It is assumed that the dimensions of
 * the matrix (number of rows and columns) are found on the first line
 * of the csv file, so no heuristics are necessary. In case of error,
 * an empty matrix is returned and @err is filled appropriately.
 *
 * Returns: The matrix read from file, or NULL.
 */

gretl_matrix *gretl_matrix_read_from_text (const char *fname, int *err)
{
    int r, c, i, j, n;
    double x;
    gretl_matrix *A = NULL;
    FILE *fp;

    fp = gretl_read_user_file(fname);

    if (fp == NULL) {
	*err = E_FOPEN;
	return NULL;
    }

    n = fscanf(fp, "%d %d\n", &r, &c);

    if (n < 2 || r <= 0 || c <= 0) {
	*err = E_DATA;
	fclose(fp);
	return NULL;
    }

    A = gretl_matrix_alloc(r, c);
    if (A == NULL) {
	*err = E_ALLOC;
	fclose(fp);
	return NULL;
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<r && !*err; i++) {
	for (j=0; j<c && !*err; j++) {
	    if (fscanf(fp, "%lf", &x) != 1) {
		*err = E_DATA;
		fprintf(stderr, "error reading row %d, column %d\n", i+1, j+1);
	    } else {
		gretl_matrix_set(A, i, j, x);
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    
    fclose(fp);

    if (*err) {
	gretl_matrix_free(A);
	A = NULL;
    }

    return A;
}

/**
 * gretl_matrix_write_as_text:
 * @A: matrix.
 * @fname: name of file to write.
 *
 * Writes the matrix @A to a plain text file by the name @fname; the
 * column separator is the tab. The number of rows and columns are
 * written on the first line of the file (which comes in handy for
 * reading the matrix).
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_matrix_write_as_text (gretl_matrix *A, const char *fname)
{
    int r = A->rows;
    int c = A->cols;
    int i, j, err = 0;
    FILE *fp;
    char d = '\t';

    fname = gretl_maybe_switch_dir(fname);

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	return E_FOPEN;
    }

    fprintf(fp, "%d%c%d\n", r, d, c);
    
    gretl_push_c_numeric_locale();

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    fprintf(fp, "%26.18E", gretl_matrix_get(A, i, j));
	    if (j == c-1) {
		fputc('\n', fp); 
	    } else {
		fputc(d, fp);
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    
    fclose(fp);

    return err;
}

static void make_numstr (char *s, double x)
{
    if (x == -0.0) {
	x = 0.0;
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
	    strncat(s, "0", 1);
	}
    }
}

static int max_numchars (const gretl_matrix *m)
{
    char s[24];
    int i, n = m->rows * m->cols;
    int c, cmax = 0;

    for (i = 0; i < n && cmax < 6; i++) {
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
	    len = strlen(names[i]);
	    if (len > maxlen) {
		maxlen = len;
	    }
	}
    }

    return maxlen;
}

static void 
real_matrix_print_to_prn (const gretl_matrix *m, const char *msg, 
			  int packed, int errout, int plain,
			  const char **heads, PRN *prn)
{
    const char **rownames = NULL;
    char numstr[32];
    double x;
    int i, j, llen = 0;

    if (prn == NULL) {
	return;
    }

    if (m == NULL || m->val == NULL) {
	if (msg != NULL && *msg != '\0') {
	    pprintf(prn, "%s: matrix is NULL\n", msg);
	} else {
	    pputs(prn, "matrix is NULL\n");
	}
	return;
    }

    if (msg != NULL && *msg != '\0' && !plain) {
	pprintf(prn, "%s (%d x %d)", msg, m->rows, m->cols);
	if (gretl_matrix_is_dated(m)) {
	    int mt1 = gretl_matrix_get_t1(m);
	    int mt2 = gretl_matrix_get_t2(m);

	    pprintf(prn, " [t1 = %d, t2 = %d]\n\n", mt1 + 1, mt2 + 1);
	} else {
	    pputs(prn, "\n\n");
	}
    }

    if (heads == NULL) {
	heads = gretl_matrix_get_colnames(m);
    }

    if (!packed) {
	rownames = gretl_matrix_get_rownames(m);
	if (rownames != NULL) {
	    llen = max_label_length(rownames, m->rows);
	}
    }

    if (heads != NULL) {
	if (rownames != NULL) {
	    bufspace(llen + 1, prn);
	}
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%12.12s ", heads[j]);
	}
	pputc(prn, '\n');
    }

    if (packed) {
	int v, n;

	v = gretl_vector_get_length(m);
	if (v == 0) {
	    pputs(prn, " not a packed matrix\n");
	    return;
	}

	n = (sqrt(1.0 + 8.0 * v) - 1.0) / 2.0;

	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		x = m->val[ijton(i, j, n)];
		make_numstr(numstr, x);
		pprintf(prn, "%12s ", numstr);
	    }
	    pputc(prn, '\n');
	}
    } else {
	int cmax = 0;

	if (heads == NULL) {
	    cmax = max_numchars(m);
	    if (cmax > 5) {
		cmax = 0;
	    }
	} 

	for (i=0; i<m->rows; i++) {
	    if (rownames != NULL) {
		pprintf(prn, "%*s ", llen, rownames[i]);
	    }
	    for (j=0; j<m->cols; j++) {
		x = gretl_matrix_get(m, i, j);
		if (cmax) {
		    sprintf(numstr, "%g", x);
		    pprintf(prn, "%*s ", cmax + 2, numstr);
		} else {
		    make_numstr(numstr, x);
		    pprintf(prn, "%12s ", numstr);
		}
	    }
	    pputc(prn, '\n');
	}
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
gretl_matrix_print_to_prn (const gretl_matrix *m, const char *msg, PRN *prn)
{
    real_matrix_print_to_prn(m, msg, 0, 0, 0, NULL, prn);
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
					PRN *prn)
{
    real_matrix_print_to_prn(m, title, 0, 0, 0, heads, prn);
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
	char wtest[32];
	double x;
	int j, n;

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
	    pprintf(prn, "%*s", n, heads[j]);
	}
	pputc(prn, '\n');
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
	real_matrix_print_to_prn(m, NULL, 0, 0, 1, NULL, prn);
    } else {
	const char **rownames = NULL;
	int llen = 0, intcast = 0;
	double x;
	int i, j, c;

	c = fmt[strlen(fmt)-1];
	if (c == 'd' || c == 'u' || c == 'x' || c == 'l') {
	    intcast = 1;
	}

	rownames = gretl_matrix_get_rownames(m);
	if (rownames != NULL) {
	    llen = max_label_length(rownames, m->rows);
	}

	maybe_print_col_heads(m, fmt, wid, prec, intcast, llen, prn);

	for (i=0; i<m->rows; i++) {
	    if (rownames != NULL) {
		pprintf(prn, "%*s ", llen, rownames[i]);
	    }
	    for (j=0; j<m->cols; j++) {
		x = gretl_matrix_get(m, i, j);
		if (intcast) {
		    if (wid >= 0 && prec >= 0) {
			pprintf(prn, fmt, wid, prec, (int) x);
		    } else if (wid >= 0 || prec >= 0) {
			c = (wid >= 0)? wid : prec;
			pprintf(prn, fmt, c, (int) x);
		    } else {
			pprintf(prn, fmt, (int) x);
		    }
		} else {
		    if (wid >= 0 && prec >= 0) {
			pprintf(prn, fmt, wid, prec, x);
		    } else if (wid >= 0 || prec >= 0) {
			c = (wid >= 0)? wid : prec;
			pprintf(prn, fmt, c, x);
		    } else {
			pprintf(prn, fmt, x);
		    }
		}
	    }
	    pputc(prn, '\n');
	}
    }
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
    for (i=0; i<m->cols; i++) {
	if (!mask[i]) {
	    for (j=0; j<m->rows; j++) {
		x = gretl_matrix_get(m, j, i);
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
    void *handle = NULL;
    int (*matrix_mp_ols) (const gretl_vector *,
			  const gretl_matrix *,
			  gretl_vector *,
			  gretl_matrix *,
			  gretl_vector *,
			  double *);
    int err;

    matrix_mp_ols = get_plugin_function("matrix_mp_ols", &handle);
    if (matrix_mp_ols == NULL) {
	return 1;
    }

    err = (*matrix_mp_ols)(y, X, b, vcv, uhat, s2);

    close_plugin(handle);

    return err;
}


