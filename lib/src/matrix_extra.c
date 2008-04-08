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
   datatypes; also printing of gretl_matrices */

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
gretl_vcv_matrix_from_model (MODEL *pmod, const char *select)
{
    gretl_matrix *vcv;
    int i, j, idx, nc;
    int ii, jj;
    int k = pmod->ncoeff;

    /* first ensure the model _has_ a vcv */
    if (makevcv(pmod, pmod->sigma)) {
	return NULL;
    }

    if (select == NULL) {
	nc = k;
    } else {
	nc = count_selection(select, k);
    }
    
    if (nc == 0) {
	return NULL;
    }

    vcv = gretl_matrix_alloc(nc, nc);
    if (vcv == NULL) {
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
gretl_coeff_vector_from_model (const MODEL *pmod, const char *select)
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
	return NULL;
    }

    b = gretl_column_vector_alloc(nc);
    if (b == NULL) {
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
    int r, c, i, j;
    int ret;
    double x;
    gretl_matrix *A = NULL;
    FILE *f;

    f = fopen(fname, "r");

    if (f == NULL) {
	*err = E_FOPEN;
	return gretl_null_matrix_new();
    }

    ret = fscanf(f, "%d %d\n", &r, &c);

    if (ret < 2 || r <= 0 || c <= 0) {
	*err = E_DATA;
	A = gretl_null_matrix_new();
    } else {
	A = gretl_matrix_alloc(r, c);
    }

    if (A == NULL) {
	*err = E_ALLOC;
    }

    if (*err) {
	fclose(f);
	return A;
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<r && !*err; i++) {
	for (j=0; j<c && !*err; j++) {
	    if (fscanf(f, "%lf", &x) != 1) {
		*err = E_DATA;
		fprintf(stderr, "error reading row %d, column %d\n", i+1, j+1);
	    } else {
		gretl_matrix_set(A, i, j, x);
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    
    fclose(f);

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

    fp = fopen(fname, "w");

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

