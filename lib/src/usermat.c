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
#include "gretl_normal.h"
#include "usermat.h"
#include "genparse.h"
#include "uservar.h"

#define MDEBUG 0
#define CONTIG_DEBUG 0

#define mspec_get_offset(m) (m->sel[0].range[0])
#define mspec_get_n_elem(m) (m->sel[0].range[1])

#define mspec_set_offset(m,o) (m->sel[0].range[0] = o)
#define mspec_set_n_elem(m,n) (m->sel[0].range[1] = n)

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

static int vec_is_exclusion (const gretl_matrix *m, int n,
			     gretl_matrix **pv, int *err)
{
    int len = gretl_vector_get_length(m);
    int i, k, neg = 0;

    for (i=0; i<len; i++) {
	if (m->val[i] < 0) {
	    neg++;
	    k = (int) fabs(m->val[i]);
	    if (k > n) {
		gretl_errmsg_sprintf(_("Index value %d is out of bounds"), k);
		*err = E_DATA;
		return 0;
	    }
	}
    }

    if (neg == len) {
	int excluded, j, nsel = n;

	for (i=1; i<=n; i++) {
	    for (j=0; j<len; j++) {
		if (m->val[j] == -i) {
		    nsel--;
		    break;
		}
	    }
	}

	if (nsel == 0) {
	    return 1;
	}

	*pv = gretl_vector_alloc(nsel);
	if (*pv == NULL) {
	    *err = E_ALLOC;
	    return 0;
	}

	k = 0;
	for (i=1; i<=n; i++) {
	    excluded = 0;
	    for (j=0; j<len; j++) {
		if (m->val[j] == -i) {
		    excluded = 1;
		    break;
		}
	    }
	    if (!excluded) {
		(*pv)->val[k++] = i;
	    }
	}
	return 1;
    }

    return 0;
}

static int bad_sel_vector (const gretl_matrix *v, int n)
{
    int i, len = gretl_vector_get_length(v);
    int exclude = 0;
    int vvi, neg = 0;

    for (i=0; i<len; i++) {
	if (v->val[i] < 0) {
	    neg++;
	}
    }

    if (neg == len) {
	exclude = 1;
    }

    for (i=0; i<len; i++) {
	vvi = exclude ? -v->val[i] : v->val[i];
	if (vvi < 1 || vvi > n) {
	    gretl_errmsg_sprintf(_("Index value %g is out of bounds"),
				 v->val[i]);
	    return E_INVARG;
	}
    }

    return 0;
}

static int bad_sel_range (int *range, int n)
{
    int i, k;

    for (i=0; i<2; i++) {
	k = range[i];
	if (k != MSEL_MAX && (k < 1 || k > n)) {
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), k);
	    return E_INVARG;
	}
    }

    return 0;
}

static int bad_sel_single (int k, int n)
{
    if (k != MSEL_MAX && (k < 1 || k > n)) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), k);
	return E_INVARG;
    }

    return 0;
}

/* convert a matrix subspec component into list of rows
   or columns */

static int *mspec_make_list (int type, union msel *sel, int n,
			     int *err)
{
    gretl_vector *ivec = NULL;
    int *slice = NULL;
    int exclude = 0;
    int i, ns = 0;

    if (type == SEL_ALL || type == SEL_NULL) {
	return NULL;
    }

    if (type == SEL_MATRIX) {
	if (sel->m == NULL) {
	    gretl_errmsg_set(_("Range is non-positive!"));
	    *err = E_DATA;
	} else {
	    if (vec_is_exclusion(sel->m, n, &ivec, err)) {
		ns = (ivec == NULL)? 0 : gretl_vector_get_length(ivec);
	    } else {
		ns = gretl_vector_get_length(sel->m);
	    }
	}
    } else {
	/* range or exclusion */
	int sr0 = sel->range[0];
	int sr1 = sel->range[1];

	if (sr1 == MSEL_MAX) {
	    sr1 = sel->range[1] = n;
	}
	if (sr0 < 0 && sr1 == sr0) {
	    /* excluding a single row or column? */
	    sr0 = -sr0;
	    if (sr0 > n) {
		gretl_errmsg_sprintf(_("Index value %d is out of bounds"),
				     sr0);
		*err = E_DATA;
	    } else {
		ns = n - 1;
		exclude = sr0;
	    }
	} else if (bad_sel_range(sel->range, n)) {
	    *err = E_DATA;
	} else {
	    ns = sr1 - sr0 + 1;
	    if (ns <= 0) {
		gretl_errmsg_sprintf(_("Range %d to %d is non-positive!"),
				     sr0, sr1);
		*err = E_DATA;
	    }
	}
    }

    if (!*err) {
	if (ns == 0) {
	    slice = gretl_null_list();
	} else {
	    slice = gretl_list_new(ns);
	}
	if (slice == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	goto bailout;
    }

    if (exclude) {
	int k = 1;

	for (i=1; i<=slice[0]; i++) {
	    if (i == exclude) {
		k++;
	    }
	    slice[i] = k++;
	}
    } else {
	for (i=0; i<slice[0]; i++) {
	    if (ivec != NULL) {
		slice[i+1] = ivec->val[i];
	    } else if (type == SEL_MATRIX) {
		slice[i+1] = sel->m->val[i];
	    } else {
		slice[i+1] = sel->range[0] + i;
	    }
	}
    }

    for (i=1; i<=slice[0] && !*err; i++) {
	if (slice[i] < 1 || slice[i] > n) {
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"),
				 slice[i]);
	    *err = 1;
	}
    }

 bailout:

    if (ivec != NULL) {
	gretl_matrix_free(ivec);
    }

    if (*err) {
	free(slice);
	slice = NULL;
    }

    return slice;
}

#define all_or_null(t) (t == SEL_ALL || t == SEL_NULL)
#define singleton_range(sel) (sel.range[0] == sel.range[1])

#define lhs_is_scalar(s,m) (s->type[0] == SEL_ELEMENT ||		\
			    (m->rows == 1 && all_or_null(s->type[0])) || \
			    (s->type[0] == SEL_RANGE && singleton_range(s->sel[0])))

#define rhs_is_scalar(s,m) (s->type[1] == SEL_ELEMENT ||		\
			    (m->cols == 1 && all_or_null(s->type[1])) || \
			    (s->type[1] == SEL_RANGE && singleton_range(s->sel[1])))

static int set_element_index (matrix_subspec *spec,
			      const gretl_matrix *m)

{
    int i = spec->sel[0].range[0];
    int j = spec->sel[1].range[0];
    int k;

    if (spec->type[0] == SEL_ELEMENT ||
	(spec->type[0] == SEL_RANGE && spec->type[1] == SEL_RANGE)) {
	k = (j-1) * m->rows + (i-1);
    } else {
	k = i > j ? (i-1) : (j-1);
    }

#if CONTIG_DEBUG
    fprintf(stderr, "element_get_index: i=%d, j=%d, k=%d\n",
	    i, j, k);
#endif

    /* write single (zero-based) index into sel[0] */
    spec->sel[0].range[0] = k;

    return 0;
}

static int spec_check_dimensions (matrix_subspec *spec,
				  const gretl_matrix *m)
{
    int colvec = m->cols == 1 && m->rows > 1;
    int st0 = spec->type[0];
    int st1 = spec->type[1];
    int err = 0;

    if (colvec || st1 != SEL_NULL) {
	/* spec->sel[0] must pertain to rows */
	if (st0 == SEL_MATRIX) {
	    err = bad_sel_vector(spec->sel[0].m, m->rows);
	} else if (st0 == SEL_RANGE || st0 == SEL_ELEMENT) {
	    err = bad_sel_range(spec->sel[0].range, m->rows);
	}
	if (!err) {
	    /* check spec->sel[1] against cols */
	    if (st1 == SEL_MATRIX) {
		err = bad_sel_vector(spec->sel[1].m, m->cols);
	    } else if (st1 == SEL_RANGE || st1 == SEL_ELEMENT) {
		err = bad_sel_range(spec->sel[1].range, m->cols);
	    }
	}
    } else {
	/* spec->sel[0] must pertain to cols */
	if (st0 == SEL_MATRIX) {
	    err = bad_sel_vector(spec->sel[0].m, m->cols);
	} else if (st0 == SEL_RANGE || st0 == SEL_ELEMENT) {
	    err = bad_sel_range(spec->sel[0].range, m->cols);
	}
    }

    return err;
}

/* Catch the case of an implicit column or row specification for
   a sub-matrix of an (n x 1) or (1 x m) matrix; also catch the
   error of giving just one row/col spec for a matrix that has
   more than one row and more than one column.
*/

int check_matrix_subspec (matrix_subspec *spec, const gretl_matrix *m)
{
    int isvec = gretl_vector_get_length(m);
    int rh_scalar, get_contig = 0;
    int err = 0;

#if CONTIG_DEBUG
    fprintf(stderr, "check_matrix_subspec: types = (%d,%d), ",
	    spec->type[0], spec->type[1]);
    if (spec->type[0] == SEL_MATRIX) {
	fputs("vector sel, ", stderr);
    } else if (spec->type[0] != SEL_NULL) {
	fprintf(stderr, "S0 = (%d,%d), ", spec->sel[0].range[0],
		spec->sel[0].range[1]);
    }
    if (spec->type[1] == SEL_MATRIX) {
	fputs("vector sel, ", stderr);
    } else if (spec->type[1] != SEL_NULL) {
	fprintf(stderr, "S1 = (%d,%d), ", spec->sel[1].range[0],
		spec->sel[1].range[1]);
    }
    fprintf(stderr, "m is %d x %d\n", m->rows, m->cols);
    fprintf(stderr, "lh scalar %d, rh scalar %d\n",
	    lhs_is_scalar(spec, m), rhs_is_scalar(spec, m));
#endif

    if (spec->type[0] == SEL_DIAG) {
	/* there's nothing to check, since this is OK even
	   for an empty matrix argument, and for a non-empty
	   matrix the data will not be contiguous
	*/
	return 0;
    }

    if (spec->type[0] == SEL_SINGLE) {
	if (!isvec) {
	    gretl_errmsg_set(_("Ambiguous matrix index"));
	    err = E_DATA;
	} else {
	    err = bad_sel_single(spec->sel[0].range[0], isvec);
	    if (!err) {
		spec->type[0] = spec->type[1] = SEL_ELEMENT;
		spec->sel[0].range[0] -= 1; /* 0-based offset */
		spec->sel[1].range[0] = 1;
	    }
	}
	/* dims checked, nothing more to do */
	return err;
    }

    err = spec_check_dimensions(spec, m);
    if (err) {
	return err;
    }

    if (m->rows == 0 || m->cols == 0) {
	/* dims checked, nothing more to do */
	return 0;
    }

    rh_scalar = rhs_is_scalar(spec, m);

    if (lhs_is_scalar(spec, m)) {
	if (rh_scalar || (m->rows == 1 && spec->type[1] == SEL_NULL)) {
	    /* we're looking at just one element */
	    set_element_index(spec, m);
	    spec->type[0] = spec->type[1] = SEL_ELEMENT;
	    spec->sel[1].range[0] = 1;
	    return 0;
	}
    }

    if ((isvec || rh_scalar) &&
	(spec->type[0] == SEL_RANGE || spec->type[0] == SEL_ALL)) {
	/* flag as contiguous values provided the rh spec doesn't
	   take the form of a matrix or column exclusion
	*/
	get_contig = spec->type[1] != SEL_MATRIX && spec->type[1] != SEL_EXCL;
#if CONTIG_DEBUG
	if (get_contig) {
	    fprintf(stderr, "got_contig: isvec %d, rh_scalar %d, rh range %d:%d\n",
		    isvec, rh_scalar, spec->sel[1].range[0], spec->sel[1].range[0]);
	}
#endif
    }

    if (spec->type[1] == SEL_NULL) {
	/* we got only one row/col spec */
	if (!isvec) {
	    gretl_errmsg_set(_("Ambiguous matrix index"));
	    return E_DATA;
	} else if (m->cols == 1) {
	    /* OK: implicitly col = 1 */
	    spec->type[1] = SEL_RANGE;
	    mspec_set_col_index(spec, 1);
	} else {
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
	}
    } else if (spec->type[1] == SEL_ALL) {
	/* select all columns */
	spec->type[1] = SEL_RANGE;
	spec->sel[1].range[0] = 1;
	spec->sel[1].range[1] = m->cols;
    } else if (spec->type[0] == SEL_ALL) {
	/* select all rows */
	spec->type[0] = SEL_RANGE;
	spec->sel[0].range[0] = 1;
	spec->sel[0].range[1] = m->rows;
    }

    if (get_contig) {
	int i, j, n;

	if (spec->sel[0].range[1] == MSEL_MAX) {
	    spec->sel[0].range[1] = m->rows;
	}
	if (spec->sel[1].range[1] == MSEL_MAX) {
	    spec->sel[1].range[1] = m->cols;
	}

	if (!isvec) {
	    if (spec->type[0] == SEL_RANGE) {
		i = spec->sel[0].range[0];
		j = spec->sel[1].range[0];
		n = spec->sel[0].range[1] - i + 1;
	    } else {
		/* must be SEL_ALL */
		i = 1;
		j = spec->sel[1].range[0];
		n = m->rows;
	    }
	} else {
	    /* we're looking at a vector */
            if (m->cols == 1) {
                if (spec->type[0] == SEL_ALL) {
                    spec->sel[0].range[0] = 1;
                    spec->sel[0].range[1] = m->rows;
                }
                if (spec->type[1] == SEL_ALL) {
		    mspec_set_col_index(spec, 1);
                }
            }
	    i = spec->sel[0].range[0];
	    if (m->rows == 1) {
		j = spec->sel[1].range[0];
		n = spec->sel[1].range[1] - j + 1;
	    } else {
		j = spec->sel[1].range[1];
		n = spec->sel[0].range[1] - i + 1;
	    }
	}
	spec->type[0] = SEL_CONTIG;
	mspec_set_offset(spec, (j-1) * m->rows + (i-1));
	mspec_set_n_elem(spec, n);
	if (spec->sel[0].range[0] < 0 || n <= 0) {
	    fprintf(stderr, "*** offset = %d, n = %d (i=%d, j=%d, m: %dx%d) ***\n",
		    mspec_get_offset(spec), n, i, j, m->rows, m->cols);
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

#if MDEBUG
    if (err) {
	fprintf(stderr, "matrix: get_slices, err=%d\n", err);
    }
#endif

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

    if (spec->type[0] == SEL_CONTIG) {
	int ini = mspec_get_offset(spec);
	int fin = ini + mspec_get_n_elem(spec);

	for (i=ini; i<fin; i++) {
	    M->val[i] = x;
	}
	return 0;
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

matrix_subspec *matrix_subspec_new (void)
{
    matrix_subspec *spec = calloc(1, sizeof *spec);

    if (spec != NULL) {
	spec->rslice = spec->cslice = NULL;
    }

    return spec;
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

int matrix_replace_submatrix (gretl_matrix *M,
			      const gretl_matrix *S,
			      matrix_subspec *spec)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int sr = gretl_matrix_rows(S);
    int sc = gretl_matrix_cols(S);
    int sscalar = 0;
    int err = 0;

#if MDEBUG
    fprintf(stderr, "\nmatrix_replace_submatrix\n");
#endif

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (spec->type[0] == SEL_CONTIG) {
	int ini = mspec_get_offset(spec);
	int n = mspec_get_n_elem(spec);

	if (gretl_vector_get_length(S) != n) {
	    return E_NONCONF;
	} else {
	    memcpy(M->val + ini, S->val, n * sizeof(double));
	    return 0;
	}
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
#if MDEBUG
	fprintf(stderr, "calling get_slices\n");
#endif
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

    if (!err && spec->rslice == NULL && spec->cslice != NULL) {
	/* the target is just specified by column(s) */
	int j, mcol, nr = M->rows;

	if (sscalar) {
	    /* write scalar into all rows of each selected col */
	    double x = S->val[0];
	    int i;

	    for (j=1; j<=spec->cslice[0]; j++) {
		mcol = spec->cslice[j] - 1;
		for (i=0; i<nr; i++) {
		    gretl_matrix_set(M, i, mcol, x);
		}
	    }
	} else {
	    /* zap from cols of S into selected cols of M */
	    double *src = S->val;
	    size_t colsize = nr * sizeof *src;

	    for (j=1; j<=spec->cslice[0]; j++) {
		mcol = spec->cslice[j] - 1;
		memcpy(M->val + mcol * nr, src, colsize);
		src += nr;
	    }
	}
    } else if (!err) {
	double x = sscalar ? S->val[0] : 0.0;
	int i, j, l, k = 0;
	int mi, mj;

	for (j=0; j<sc; j++) {
	    mj = (spec->cslice == NULL)? k++ : spec->cslice[j+1] - 1;
	    l = 0;
	    for (i=0; i<sr; i++) {
		mi = (spec->rslice == NULL)? l++ : spec->rslice[i+1] - 1;
		if (!sscalar) {
		    x = gretl_matrix_get(S, i, j);
		}
		gretl_matrix_set(M, mi, mj, x);
	    }
	}
    }

    return err;
}

static void matrix_transcribe_dates (gretl_matrix *targ,
				     const gretl_matrix *src)
{
    int mt1 = gretl_matrix_get_t1(src);
    int mt2 = gretl_matrix_get_t2(src);

    gretl_matrix_set_t1(targ, mt1);
    gretl_matrix_set_t2(targ, mt2);
}

gretl_matrix *matrix_get_submatrix (const gretl_matrix *M,
				    matrix_subspec *spec,
				    int prechecked,
				    int *err)
{
    gretl_matrix *S = NULL;
    int r, c;

    if (M == NULL || spec == NULL) {
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
    } else if (spec->type[0] == SEL_CONTIG) {
	return matrix_get_chunk(M, spec, err);
    } else if (spec->type[0] == SEL_ELEMENT) {
	int i = mspec_get_element(spec);
	double x = matrix_get_element(M, i, err);

	if (!*err) {
	    S = gretl_matrix_from_scalar(x);
	}
	return S;
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
	int j, mj;

	if (spec->cslice != NULL && spec->rslice == NULL) {
	    /* copying entire columns */
	    double *dest = S->val;
	    size_t csize = r * sizeof *dest;

	    for (j=0; j<c; j++) {
		mj = spec->cslice[j+1] - 1;
		memcpy(dest, M->val + mj * r, csize);
		dest += r;
	    }
	} else {
	    int i, mi, l, k = 0;
	    double x;

	    for (j=0; j<c; j++) {
		mj = (spec->cslice == NULL)? k++ : spec->cslice[j+1] - 1;
		l = 0;
		for (i=0; i<r; i++) {
		    mi = (spec->rslice == NULL)? l++ : spec->rslice[i+1] - 1;
		    x = gretl_matrix_get(M, mi, mj);
		    gretl_matrix_set(S, i, j, x);
		}
	    }
	}
    }

    if (S != NULL) {
	/* try transcribing metadata on @M if applicable */
	if (S->rows == M->rows && gretl_matrix_is_dated(M)) {
	    matrix_transcribe_dates(S, M);
	}
	if (S->cols == M->cols) {
	    const char **cnames = gretl_matrix_get_colnames(M);

	    if (cnames != NULL) {
		char **cpy = strings_array_dup((char **) cnames, M->cols);

		if (cpy != NULL) {
		    gretl_matrix_set_colnames(S, cpy);
		}
	    }
	}
    }

    return S;
}

/* Return the element of @M that occupies 0-based position @i in
   its internal vec representation.
*/

double matrix_get_element (const gretl_matrix *M, int i, int *err)
{
    double x = NADBL;

    if (M == NULL) {
	*err = E_DATA;
    } else if (i < 0 || i >= M->rows * M->cols) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	*err = E_INVARG;
    } else {
	x = M->val[i];
    }

    return x;
}

/* Copy a contiguous chunk of data out of @M, the offset and
   extent of which are given by @spec. At present we come
   here only if it's known that @M is a vector.
*/

gretl_matrix *matrix_get_chunk (const gretl_matrix *M,
				matrix_subspec *spec,
				int *err)
{
    int offset = spec->sel[0].range[0];
    int n = spec->sel[0].range[1];
    size_t sz = n * sizeof(double);
    gretl_matrix *ret;

    if (offset < 0) {
	fprintf(stderr, "matrix_get_chunk: offset = %d\n", offset);
	*err = E_DATA;
	return NULL;
    }

    if (M->rows == 1) {
	ret = gretl_matrix_alloc(1, n);
    } else {
	ret = gretl_matrix_alloc(n, 1);
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	memcpy(ret->val, M->val + offset, sz);
	if (M->rows > 1 && n == M->rows && offset == 0 &&
	    gretl_matrix_is_dated(M)) {
	    matrix_transcribe_dates(ret, M);
	}
    }

    return ret;
}

/* Handle the case where we got a single string as argument
   to colnames() or rownames(), for a matrix with more than
   one column or row: construct specific names by appending
   a column or row index.
*/

static char **expand_names (const char *s, int n, int *err)
{
    char **S = NULL;

    if (!gretl_is_ascii(s)) {
	*err = E_INVARG;
	return NULL;
    }

    S = strings_array_new(n);

    if (S == NULL) {
	*err = E_ALLOC;
    } else {
	char tmp[12];
	int i, m;

	for (i=0; i<n && !*err; i++) {
	    sprintf(tmp, "%d", i+1);
	    m = strlen(tmp);
	    sprintf(tmp, "%.*s%d", 9-m, s, i+1);
	    S[i] = gretl_strdup(tmp);
	    if (S[i] == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (*err) {
	   strings_array_free(S, n);
	   S = NULL;
	}
    }

    return S;
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
	} else if (ns == 1 && n > 1) {
	    strings_array_free(S, ns);
	    S = expand_names(s, n, &err);
	} else if (ns != n) {
	    err = E_NONCONF;
	    strings_array_free(S, ns);
	}

	if (!err) {
	    if (byrow) {
		gretl_matrix_set_rownames(M, S);
	    } else {
		gretl_matrix_set_colnames(M, S);
	    }
	}
    }

    return err;
}

int umatrix_set_names_from_array (gretl_matrix *M,
				  void *data,
				  int byrow)
{
    gretl_array *A = data;
    int n, err = 0;

    n = (byrow)? M->rows : M->cols;

    if (A == NULL || gretl_array_get_length(A) == 0) {
	if (byrow) {
	    gretl_matrix_set_rownames(M, NULL);
	} else {
	    gretl_matrix_set_colnames(M, NULL);
	}
    } else {
	char **AS;
	int i, ns;

	AS = gretl_array_get_strings(A, &ns);

	if (ns != n) {
	    err = E_NONCONF;
	} else {
	    for (i=0; i<ns && !err; i++) {
		if (AS[i] == NULL || AS[i][0] == '\0') {
		    fprintf(stderr, "Missing string in colnames/rownames\n");
		    err = E_INVARG;
		}
	    }
	}

	if (!err) {
	    char **S = strings_array_dup(AS, ns);

	    if (S == NULL) {
		err = E_ALLOC;
	    } else if (byrow) {
		gretl_matrix_set_rownames(M, S);
	    } else {
		gretl_matrix_set_colnames(M, S);
	    }
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

char *user_matrix_get_row_name (const gretl_matrix *M, int row,
				int *err)
{
    char *ret = NULL;

    if (M == NULL || row < 1 || row > M->rows) {
	*err = E_DATA;
    } else {
	const char **S = gretl_matrix_get_rownames(M);

	if (S == NULL) {
	    ret = gretl_strdup("");
	} else {
	    ret = gretl_strdup(S[row-1]);
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

static void matrix_cannibalize (gretl_matrix *targ, gretl_matrix *src)
{
    gretl_matrix_destroy_info(targ);

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

static gretl_matrix *get_sized_matrix (gretl_matrix *m,
				       int r, int c,
				       int *err)
{
    if (m->rows == r && m->cols == c) {
	return m;
    } else {
	gretl_matrix *m1 = gretl_matrix_alloc(r, c);

	if (m1 == NULL) {
	    *err = E_ALLOC;
	}
	return m1;
    }
}

gretl_matrix *user_matrix_QR_decomp (const gretl_matrix *m,
				     gretl_matrix *R,
				     int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *Rarg = NULL;
    gretl_matrix **pR = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (R != NULL) {
	pR = &Rarg;
    }

    if (!*err) {
	*err = real_user_matrix_QR_decomp(m, &Q, pR);
    }

    if (Rarg != NULL) {
	if (!*err) {
	    gretl_matrix_replace_content(R, Rarg);
	}
	gretl_matrix_free(Rarg);
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
			       gretl_matrix *U,
			       gretl_matrix *V,
			       int *err)
{
    gretl_matrix *S = NULL;
    gretl_matrix *Uarg = NULL;
    gretl_matrix *Varg = NULL;
    gretl_matrix **pU = NULL;
    gretl_matrix **pV = NULL;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (U != NULL) {
	pU = &Uarg;
    }

    if (V != NULL) {
	pV = &Varg;
    }

    if (!*err) {
	*err = gretl_matrix_SVD(m, pU, &S, pV);
    }

    if (!*err && (Uarg != NULL || Varg != NULL)) {
	int tall = m->rows - m->cols;
	int minrc = (m->rows > m->cols)? m->cols : m->rows;

	if (Uarg != NULL) {
	    if (tall > 0) {
		*err = gretl_matrix_realloc(Uarg, m->rows, minrc);
	    }
	    if (!*err) {
		gretl_matrix_replace_content(U, Uarg);
		gretl_matrix_free(Uarg);
	    }
	}
	if (Varg != NULL) {
	    if (tall < 0) {
		*err = revise_SVD_V(&Varg, minrc, m->cols);
	    }
	    if (!*err) {
		gretl_matrix_replace_content(V, Varg);
		gretl_matrix_free(Varg);
	    }
	}
    }

    return S;
}

gretl_matrix *user_matrix_ols (const gretl_matrix *Y,
			       const gretl_matrix *X,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       gretlopt opt,
			       int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *Uarg = NULL;
    gretl_matrix *Varg = NULL;
    gretl_matrix **Vp = NULL;
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

    if (U != NULL) {
	Uarg = get_sized_matrix(U, T, g, err);
	if (*err) {
	    return NULL;
	}
    }

    if (V != NULL) {
	if (g > 1) {
	    /* multiple dependent variables */
	    Vp = &Varg;
	} else {
	    /* a single dependent variable */
	    int nv = g * k;

	    Varg = get_sized_matrix(V, nv, nv, err);
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
	    /* null model: FIXME */
	    if (Uarg != NULL) {
		gretl_matrix_copy_values(U, Y);
	    }
	    if (V != NULL) {
		Varg = gretl_null_matrix_new();
		if (Varg == NULL) {
		    *err = E_ALLOC;
		}
	    }
	} else if (g == 1) {
	    /* single regressand */
	    if (opt & OPT_M) {
		/* use multiple precision */
		*err = gretl_matrix_mp_ols(Y, X, B, Varg, Uarg, ps2);
	    } else {
		*err = gretl_matrix_ols(Y, X, B, Varg, Uarg, ps2);
	    }
	} else {
	    /* multiple regressands: @V will actually be (X'X)^{-1} */
	    *err = gretl_matrix_multi_ols(Y, X, B, Uarg, Vp);
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
    }

    if (Uarg != NULL && Uarg != U) {
	if (!*err) {
	    gretl_matrix_replace_content(U, Uarg);
	}
	gretl_matrix_free(Uarg);
    }
    if (Varg != NULL && Varg != V) {
	if (!*err) {
	    gretl_matrix_replace_content(V, Varg);
	}
	gretl_matrix_free(Varg);
    }

    return B;
}

gretl_matrix *user_matrix_rls (const gretl_matrix *Y,
			       const gretl_matrix *X,
			       const gretl_matrix *R,
			       const gretl_matrix *Q,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       int *err)
{
    gretl_matrix *B = NULL;
    gretl_matrix *Uarg = NULL;
    gretl_matrix *Varg = NULL;
    gretl_matrix **Vp = NULL;
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

    if (U != NULL) {
	Uarg = get_sized_matrix(U, T, g, err);
	if (*err) {
	    return NULL;
	}
    }

    if (V != NULL) {
	Vp = &Varg;
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* note: &V will actually be M (X'X)^{-1} M' */
	*err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B,
						 Uarg, Vp);
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
    }

    if (Uarg != NULL && Uarg != U) {
	if (!*err) {
	    gretl_matrix_replace_content(U, Uarg);
	}
	gretl_matrix_free(Uarg);
    }
    if (Varg != NULL && Varg != V) {
	if (!*err) {
	    gretl_matrix_replace_content(V, Varg);
	}
	gretl_matrix_free(Varg);
    }

    return B;
}

gretl_matrix *user_matrix_GHK (const gretl_matrix *C,
			       const gretl_matrix *A,
			       const gretl_matrix *B,
			       const gretl_matrix *U,
			       gretl_matrix *dP,
			       int *err)
{
    gretl_matrix *P = NULL;
    gretl_matrix *dParg = NULL;
    int n, m, npar;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(C)) {
	*err = E_DATA;
	return NULL;
    }

    n = A->rows;
    m = C->rows;
    npar = m + m + m*(m+1)/2;

    if (dP != NULL) {
	dParg = get_sized_matrix(dP, n, npar, err);
    }

    if (!*err) {
	P = gretl_GHK2(C, A, B, U, dParg, err);
    }

    if (dParg != NULL && dParg != dP) {
	if (!*err) {
	    gretl_matrix_replace_content(dP, dParg);
	}
	gretl_matrix_free(dParg);
    }

    return P;
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
user_matrix_eigen_analysis (const gretl_matrix *m,
			    gretl_matrix *R,
			    int symm, int *err)
{
    gretl_matrix *C = NULL;
    gretl_matrix *E = NULL;
    int vecs = 0;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_na_check(m)) {
	*err = E_NAN;
	return NULL;
    }

    if (R != NULL) {
	vecs = 1;
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
	gretl_matrix_replace_content(R, C);
    }

    if (!vecs) {
	gretl_matrix_free(C);
    }

    return E;
}

gretl_matrix *user_gensymm_eigenvals (const gretl_matrix *A,
				      const gretl_matrix *B,
				      gretl_matrix *V,
				      int *err)
{
    gretl_matrix *E = NULL;
    gretl_matrix *Varg = NULL;
    int newV = 0;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_na_check(A) || gretl_matrix_na_check(B)) {
	*err = E_NAN;
	return NULL;
    }

    if (V != NULL) {
	Varg = get_sized_matrix(V, B->cols, A->rows, err);
    }

    if (!*err) {
	E = gretl_gensymm_eigenvals(A, B, Varg, err);
    }

    if (Varg != NULL && Varg != V) {
	if (!*err) {
	    gretl_matrix_replace_content(V, Varg);
	}
	gretl_matrix_free(Varg);
    }

    return E;
}
