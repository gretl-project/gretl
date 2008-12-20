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

/* functions that convert between gretl_matrix and other
   datatypes; also printing of gretl_matrices and some
   functions relating to masking and cutting rows and/or
   columns of matrices.
*/

/**
 * gretl_vector_from_array:
 * @x: pointer to array of elements.
 * @n: number of elements.
 * @mod: modifier flag: either %GRETL_MOD_NONE, or %GRETL_MOD_SQUARE
 * to use the squares of the elements of @x.
 *
 * Returns: pointer to a newly allocated gretl_vector containing
 * the elements of x (or their squares), or %NULL on failure.  
 * Missing valies in @x are skipped.
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
 * of the given data series for the given range, or %NULL on failure.  
 */

gretl_vector *gretl_vector_from_series (const double *x, 
					int t1, int t2)
{
    gretl_matrix *v;
    int t, T = t2 - t1 + 1;

    if (T <= 0) {
	return NULL;
    }

    v = gretl_column_vector_alloc(T);

    if (v != NULL) {
	for (t=0; t<T; t++) {
	    v->val[t] = x[t + t1];
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
 * the values in @X, or %NULL on allocation failure.
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
 * Returns: allocated 1x1 gretl_matrix, the single element
 * of which is set to @x, or %NULL on allocation failure
 * or if @x = #NADBL.
 */

gretl_matrix *gretl_matrix_from_scalar (double x) 
{
    gretl_matrix *m = NULL;

    if (!na(x)) {
	m = gretl_matrix_alloc(1, 1);
	if (m != NULL) {
	    m->val[0] = x;
	}
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
 * (or %NULL for the full matrix).
 * @err: location to receive error code.
 *
 * Produces all or part of the covariance matrix for @pmod 
 * in the form of a gretl_matrix.  Storage is allocated, to be freed
 * by the caller.  If @select is not %NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows (and columns), and zero elements otherwise.
 * 
 * Returns: the covariance matrix, or %NULL on error.
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
 * (or %NULL for the full vector).
 * @err: location to receive error code.
 *
 * Produces all or part of the coefficient vector for @pmod  
 * in the form of a gretl column vector.  Storage is allocated, to be freed
 * by the caller.  If @select is non-%NULL, it should be an array
 * with non-zero elements in positions corresponding to the
 * desired rows and zero elements otherwise.
 * 
 * Returns: the coefficient vector, or %NULL on error.
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
 * @Z: data array.
 * @pdinfo: pointer to data information struct.
 * @means: pointer to pick up vector of means, or %NULL to discard.
 * @errp: pointer to receive non-zero error code in case of
 * failure, or %NULL.
 *
 * Returns: the variance-covariance matrix of the listed variables
 * (over the currently defined data sample), or %NULL in case of
 * failure.
 */

gretl_matrix *
gretl_covariance_matrix_from_varlist (const int *list, const double **Z, 
				      const DATAINFO *pdinfo, 
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
	xbar->val[i] = gretl_mean(pdinfo->t1, pdinfo->t2, Z[list[i+1]]);
	if (na(xbar->val[i])) {
	    err = E_DATA;
	} 
    }

    for (i=0; i<k && !err; i++) {
	for (j=i; j<k; j++) {
	    vv = 0.0;
	    n = 0;
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		x = Z[list[i+1]][t];
		y = Z[list[j+1]][t];
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
 * Returns: 0 on sucess, 1 if the row is out of bounds.
 */

int gretl_matrix_row_to_array (const gretl_matrix *m, int i, double *x)
{
    int j, err = 0;

    if (i < 0 || i >= gretl_matrix_rows(m)) {
	err = 1;
    } else {
	for (j=0; j<m->cols; j++) {
	    x[j] = gretl_matrix_get(m, i, j);
	}
    }

    return err;
}

static int get_mask_count (const char *mask, int n)
{
    int i, k = 0;

    for (i=0; i<n; i++) {
	if (mask[i]) k++;
    }

    return k;
}

enum {
    OBS_MISS_MASK,
    OBS_MISS_ERR,
    OBS_MISS_SKIP
};

static gretl_matrix *
real_gretl_matrix_data_subset (const int *list, const double **Z,
			       int t1, int t2, const char *mask,
			       int op, int *err)
{
    gretl_matrix *M;
    double x;
    int T = t2 - t1 + 1;
    int k = list[0];
    int skip;
    int j, s, t;

    if (k <= 0 || T <= 0) {
	*err = E_DATA;
	return NULL;
    }	

    *err = 0;

    if (mask != NULL) {
	T -= get_mask_count(mask, T);
    } else if (op == OBS_MISS_SKIP || op == OBS_MISS_ERR) {
	for (t=t1; t<=t2 && !*err; t++) {
	    for (j=0; j<k; j++) {
		x = Z[list[j+1]][t];
		if (na(x)) {
		    if (op == OBS_MISS_SKIP) {
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
	} else if (op == OBS_MISS_SKIP) {
	    for (j=0; j<k; j++) {
		x = Z[list[j+1]][t];
		if (na(x)) {
		    skip = 1;
		    break;
		}
	    }
	}
	if (!skip) {
	    for (j=0; j<k; j++) {
		x = Z[list[j+1]][t];
		gretl_matrix_set(M, s, j, x);
	    }
	    if (s == 0) {
		M->t1 = t;
	    } else if (s == T - 1) {
		M->t2 = t;
		break;
	    }
	    s++;
	}
    }

    if (!*err) {
	for (j=0; j<k && !*err; j++) {
	    for (t=0; t<T && !*err; t++) {
		x = gretl_matrix_get(M, t, j);
		if (na(x)) {
		    *err = E_MISSDATA;
		}
	    }
	}
    }	

    if (*err) {
	gretl_matrix_free(M);
	M = NULL;
    } 

    return M;
}

/**
 * gretl_matrix_data_subset:
 * @list: list of variable to process.
 * @Z: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 * @mask: missing observations mask, or %NULL.
 * @err: location to recieve error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @Z specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns.  If @mask is not
 * %NULL then it should be an array of char of length (@t2 - @t1
 * + 1) with 1s in the positions of observations to exclude
 * from the subset and zeros elsewhere. This apparatus can be
 * used to exclude missing observations.
 *
 * Returns: allocated matrix or %NULL on failure. 
 */

gretl_matrix *gretl_matrix_data_subset (const int *list, const double **Z,
					int t1, int t2, const char *mask,
					int *err)
{
    return real_gretl_matrix_data_subset(list, Z, t1, t2, mask, 
					 OBS_MISS_MASK, err);
}

/**
 * gretl_matrix_data_subset_no_missing:
 * @list: list of variable to process.
 * @Z: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @Z specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns.  If any missing 
 * values are encountered this constitutes an error.
 *
 * Returns: allocated matrix or %NULL on failure. 
 */

gretl_matrix *
gretl_matrix_data_subset_no_missing (const int *list, const double **Z,
				     int t1, int t2, int *err)
{
    return real_gretl_matrix_data_subset(list, Z, t1, t2, NULL, 
					 OBS_MISS_ERR, err);
}

/**
 * gretl_matrix_data_subset_skip_missing:
 * @list: list of variable to process.
 * @Z: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 * @err: location to receive error code.
 *
 * Creates a gretl matrix holding the subset of variables from
 * @Z specified by @list, over the sample range @t1 to @t2,
 * inclusive.  Variables are in columns.  If there is a missing
 * value for any variable on a given row, that row is skipped.
 *
 * Returns: allocated matrix or %NULL on failure. 
 */

gretl_matrix *
gretl_matrix_data_subset_skip_missing (const int *list, const double **Z,
				       int t1, int t2, int *err)
{
    return real_gretl_matrix_data_subset(list, Z, t1, t2, NULL, 
					 OBS_MISS_SKIP, err);
}

/**
 * gretl_plotfit_matrices:
 * @yno: ID number of the y variable.
 * @xno: ID number of the y variable.
 * @fit: type of fit sought.
 * @Z: data array.
 * @t1: starting observation.
 * @t2: ending observation.
 * @py: location to receive y vector.
 * @pX: location to receive X matrix.
 *
 * Creates a vector y and matrix X based on the input @yno, @xno
 * and @fit, using the given sample range.  An observation is
 * skipped if any of the variables in @list are missing at that
 * observation.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plotfit_matrices (int yno, int xno, FitType fit,
			    const double **Z, int t1, int t2, 
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

    mask = calloc(T, 1);
    if (mask == NULL) {
	return E_ALLOC;
    }

    for (s=0; s<T; s++) {
	t = s + t1;
	if (na(Z[yno][t]) || na(Z[xno][t])) {
	    mask[s] = 1;
	} else {
	    n++;
	}
    }

    if (n == 0) {
	free(mask);
	return E_MISSDATA;
    }

    if (fit == PLOT_FIT_QUADRATIC) {
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
	    y->val[i] = Z[yno][t];
	    if (fit != PLOT_FIT_LOESS) {
		gretl_matrix_set(X, i, j++, 1.0);
	    }
	    xt = Z[xno][t];
	    if (fit == PLOT_FIT_INVERSE) {
		gretl_matrix_set(X, i, j++, 1.0 / xt);
	    } else {
		gretl_matrix_set(X, i, j++, xt);
	    }
	    if (fit == PLOT_FIT_QUADRATIC) {
		gretl_matrix_set(X, i, j, xt * xt);
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

/* delete the columns of X specified in @list */

int gretl_matrix_delete_columns (gretl_matrix *X, int *list)
{
    size_t csz = X->rows * sizeof *X->val;
    void *dest, *src;
    int i, j, n, col;

    for (i=1; i<=list[0]; i++) {
	col = list[i];
	if (col < 0 || col >= X->cols) {
	    return E_NONCONF;
	}
    }

    for (i=1; i<=list[0]; i++) {
	col = list[i];
	dest = X->val + col * X->rows;
	src = X->val + (col + 1) * X->rows;
	n = X->cols - col - 1;
	if (n > 0) {
	    memmove(dest, src, n * csz);
	}
	for (j=i+1; j<=list[0]; j++) {
	    list[j] -= 1;
	}
    }

    X->cols -= list[0];

    return 0;
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
 * Returns: The matrix read from file, or %NULL.
 */

gretl_matrix *gretl_matrix_read_from_text (const char *fname, int *err)
{
    int r, c, i, j, n;
    double x;
    gretl_matrix *A = NULL;
    FILE *fp;

    fp = gretl_fopen(fname, "r");

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
		gretl_matrix_free(A);
		A = NULL;
	    } else {
		gretl_matrix_set(A, i, j, x);
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    
    fclose(fp);

    return A;
}

/**
 * gretl_matrix_write_as text:
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

    fname = gretl_maybe_switch_dir(fname);

    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	return E_FOPEN;
    }

    fprintf(fp, "%d\t%d\n", r, c);
    
    gretl_push_c_numeric_locale();

    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    fprintf(fp, "%26.18E\t", gretl_matrix_get(A, i, j));
	}
	fputc('\n', fp);
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

static void 
real_matrix_print_to_prn (const gretl_matrix *m, const char *msg, 
			  int packed, int errout, int plain,
			  const char **heads, PRN *prn)
{
    char numstr[32];
    double x;
    int i, j;

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
	if (!(m->t1 == 0 && m->t2 == 0)) {
	    pprintf(prn, " [t1 = %d, t2 = %d]\n\n", m->t1 + 1, m->t2 + 1);
	} else {
	    pputs(prn, "\n\n");
	}
    }

    if (heads == NULL) {
	heads = user_matrix_get_column_names(m);
    }

    if (heads != NULL) {
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
 * @msg: accompanying message text (or %NULL if no message is wanted).
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
 * @m: matrix to print.
 * @title: accompanying title (or %NULL if no title is wanted).
 * @heads: array of strings to identify the columns.
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
				   int icast, PRN *prn)
{
    const char **heads;

    heads = user_matrix_get_column_names(m);
    
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

	n = strlen(wtest);
	for (j=0; j<m->cols; j++) {
	    pprintf(prn, "%*s", n, heads[j]);
	}
	pputc(prn, '\n');
    }
}

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
	int intcast = 0;
	double x;
	int i, j, c;

	c = fmt[strlen(fmt)-1];
	if (c == 'd' || c == 'u' || c == 'x' || c == 'l') {
	    intcast = 1;
	}

	maybe_print_col_heads(m, fmt, wid, prec, intcast, prn);

	for (i=0; i<m->rows; i++) {
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
 * gretl_packed_matrix_print:
 * @m: packed matrix to print.
 * @msg: accompanying message text (or %NULL if no message is wanted).
 *
 * Prints the symmetric matrix @m (packed as lower triangle)
 * to stderr.
 */

void gretl_packed_matrix_print (const gretl_matrix *m, const char *msg)
{
    PRN *prn;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_STDERR, &err);

    if (!err) {
	real_matrix_print_to_prn(m, msg, 1, 0, 0, NULL, prn);
	gretl_print_destroy(prn);
    }
}

/**
 * debug_print_matrix:
 * @m: matrix to print.
 * @msg: accompanying message text (or %NULL if no message is wanted).
 *
 * Prints the matrix @m to stderr, as with gretl_matrix_print(), but
 * appends the address of the matrix struct.
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

static int count_unmasked_rows (const char *mask, int n)
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

    n = count_unmasked_rows(mask, m->rows);

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
 * gretl_matrix_cut_rows_cols:
 * @m: square matrix to process.
 * @mask: character array of length equal to the dimension
 * of @m, with 1s indicating rows and columns to be cut, 0s
 * for rows/columns to be retained.
 *
 * In-place reduction of @m based on @mask: the masked rows
 * and columns are cut out of @m.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_cut_rows_cols (gretl_matrix *m, const char *mask)
{
    gretl_matrix *tmp;
    double x;
    int i, j, k, l, n;

    if (m == NULL || mask == NULL || m->rows != m->cols) {
	return E_DATA;
    }

    n = count_unmasked_rows(mask, m->rows);

    /* allocate smaller temp matrix */
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

    /* resize the original matrix, copy the values back,
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
 * %E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or %NULL.
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
 * gretl_matrix_rank_mask:
 * @m: matrix to process.
 * @err: location to receive error code.
 *
 * Performs a QR decomposition of matrix @m and uses this
 * to assess the rank of @m.  If @m is not of full rank,
 * constructs a mask of length equal to the numbers of
 * columns in @m, with 1s in positions corresponding 
 * to diagonal elements of R that are effectively 0, and
 * 0s elsewhere.  If @m is of full column rank, %NULL is
 * returned.
 * 
 * %E_ALLOC is written to @err in case a mask should have
 * been constructed but allocation failed.
 *
 * Returns: allocated mask or %NULL.
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
 * or %NULL if this is not needed.
 * @uhat: vector to hold the regression residuals, or %NULL if 
 * these are not needed.
 * @s2: pointer to receive residual variance, or %NULL.  Note:
 * if @s2 is %NULL, the "vcv" estimate will be plain (X'X)^{-1}.
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


