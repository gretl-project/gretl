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
#include "gretl_cmatrix.h"
#include "gretl_normal.h"
#include "usermat.h"
#include "uservar.h"

#define MDEBUG 0
#define CONTIG_DEBUG 0

/* mspec convenience macros */

#define mspec_get_offset(m) (m->lsel.range[0])
#define mspec_get_n_elem(m) (m->lsel.range[1])

#define mspec_set_offset(m,o) (m->lsel.range[0] = o)
#define mspec_set_n_elem(m,n) (m->lsel.range[1] = n)

#define mspec_set_element(m,i) (m->lsel.range[0] = i)

#define spec_is_single(s) (s->ltype==SEL_SINGLE || s->rtype==SEL_SINGLE)
#define spec_single_val(s) (s->ltype==SEL_SINGLE ? s->lsel.range[0] : \
	                    s->rsel.range[0])

#define singleton_left_range(s) (s->ltype == SEL_RANGE && \
				 s->lsel.range[0] == s->lsel.range[1])
#define singleton_right_range(s) (s->rtype == SEL_RANGE && \
				  s->rsel.range[0] == s->rsel.range[1])

#define all_or_null(t) (t == SEL_ALL || t == SEL_NULL)

#define lhs_is_scalar(s,m) (s->ltype == SEL_ELEMENT || \
			    (m->rows == 1 && s->ltype == SEL_ALL) || \
			    singleton_left_range(s))
#define rhs_is_scalar(s,m) (s->rtype == SEL_ELEMENT || \
			    (m->cols == 1 && all_or_null(s->rtype)) || \
			    singleton_right_range(s))

#define rowmax(s,m) (s->range[1] == MSEL_MAX ? m->rows : s->range[1])
#define colmax(s,m) (s->range[1] == MSEL_MAX ? m->cols : s->range[1])

/* end mspec convenience macros */

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

/* @s is the user-supplied selection vector.
   @n is the max possible value of row or col.
   @pslice is location to receive positive selection
   list.
*/

static int handle_vector_exclusion (const gretl_vector *s,
				    int n, int **pslice)
{
    int slen = gretl_vector_get_length(s);
    int i, j, k, nsel = n;
    int err = 0;

    /* do we need this check here? */
    for (i=0; i<slen; i++) {
	k = (int) fabs(s->val[i]);
	if (k > n) {
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), k);
	    return E_BOUNDS;
	}
    }

    /* how many of the @n values are selected? */
    for (i=1; i<=n; i++) {
	for (j=0; j<slen; j++) {
	    if (s->val[j] == -i) {
		/* element i is dropped */
		nsel--;
		break;
	    }
	}
    }

    if (nsel == 0) {
	*pslice = gretl_null_list();
    } else {
	/* compose positive selection list */
	int *slice = gretl_list_new(nsel);
	int excluded, k = 1;

	if (slice != NULL) {
	    for (i=1; i<=n; i++) {
		excluded = 0;
		for (j=0; j<slen; j++) {
		    if (s->val[j] == -i) {
			excluded = 1;
			break;
		    }
		}
		if (!excluded) {
		    slice[k++] = i;
		}
	    }
	    *pslice = slice;
	}
    }

    if (*pslice == NULL) {
	err = E_ALLOC;
    }

    return err;
}

static int bad_sel_vector (const gretl_vector *v, int n)
{
    int i, k, len = gretl_vector_get_length(v);
    int exclude = v->val[0] < 0;

    for (i=0; i<len; i++) {
	k = exclude ? -v->val[i] : v->val[i];
	if (k < 1 || k > n) {
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), k);
	    return E_BOUNDS;
	}
    }

    return 0;
}

static int bad_sel_range (int *range, int n)
{
    int i, k, err = 0;

    for (i=0; i<2; i++) {
	k = range[i];
	if (k != MSEL_MAX && (k < 1 || k > n)) {
	    err = E_BOUNDS;
	    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), k);
	    break;
	}
    }

    return err;
}

static int bad_sel_single (int *pk, int n)
{
    int err = 0;

    if (*pk != MSEL_MAX && (*pk < 1 || *pk > n)) {
	err = E_BOUNDS;
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), *pk);
    }

    return err;
}

/* convert a matrix subspec component into list of rows
   or columns: n is the maximal dimension of selectable items
*/

int *mspec_make_list (int type, union msel *sel, int n,
		      int *err)
{
    int *slice = NULL;
    int single_exclude = 0;
    int i, ns = 0;

    if (type == SEL_ALL || type == SEL_NULL) {
	return NULL;
    }

    if (type == SEL_MATRIX) {
	if (sel->m == NULL) {
	    gretl_errmsg_set(_("Range is non-positive!"));
	    *err = E_DATA;
	} else if (sel->m->val[0] < 0) {
	    *err = handle_vector_exclusion(sel->m, n, &slice);
	    return slice; /* we're done */
	} else {
	    ns = gretl_vector_get_length(sel->m);
	}
    } else if (type == SEL_STR) {
        *err = E_TYPES;
        return NULL;
    } else {
	/* range or single exclusion */
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
		*err = E_BOUNDS;
	    } else {
		ns = n - 1;
		single_exclude = sr0;
	    }
	} else {
	    *err = bad_sel_range(sel->range, n);
	    if (!*err) {
		ns = sel->range[1] - sel->range[0] + 1;
		if (ns <= 0) {
		    gretl_errmsg_sprintf(_("Range %d to %d is non-positive!"),
					 sel->range[0], sel->range[1]);
		    *err = E_INVARG;
		}
	    }
	}
    }

    if (!*err) {
	slice = gretl_list_new(ns);
	if (slice == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* compose and check slice */
	if (single_exclude > 0) {
	    int k = 1;

	    for (i=1; i<=slice[0]; i++) {
		if (i == single_exclude) {
		    k++;
		}
		slice[i] = k++;
	    }
	} else {
	    for (i=0; i<slice[0]; i++) {
		if (type == SEL_MATRIX) {
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
		*err = E_BOUNDS;
	    }
	}
    }

    if (*err) {
	free(slice);
	slice = NULL;
    }

    return slice;
}

static int set_element_index (matrix_subspec *spec,
			      const gretl_matrix *m)

{
    int i = spec->lsel.range[0];
    int j = spec->rsel.range[0];
    int k;

    if (spec->ltype == SEL_ALL && spec->rtype == SEL_ALL) {
	k = 0;
    } else if (spec->ltype == SEL_ELEMENT ||
	(spec->ltype == SEL_RANGE && spec->rtype == SEL_RANGE)) {
	k = (j-1) * m->rows + (i-1);
    } else {
	if (spec->ltype == SEL_CONTIG) {
	    fprintf(stderr, "CONTIG on left\n");
	}
	if (spec->rtype == SEL_CONTIG) {
	    fprintf(stderr, "CONTIG on right\n");
	}
	k = i > j ? (i-1) : (j-1);
    }

#if CONTIG_DEBUG
    fprintf(stderr, "set_element_index: i=%d, j=%d -> k=%d\n",
	    i, j, k);
#endif

    /* record single (zero-based) index */
    mspec_set_offset(spec, k);
    mspec_set_n_elem(spec, 1);

    return 0;
}

static int spec_check_dimensions (matrix_subspec *spec,
				  const gretl_matrix *m)
{
    int lt = spec->ltype;
    int rt = spec->rtype;
    int err = 0;

    /* check spec->lsel against rows */
    if (lt == SEL_MATRIX) {
	err = bad_sel_vector(spec->lsel.m, m->rows);
    } else if (lt == SEL_RANGE || lt == SEL_ELEMENT) {
	err = bad_sel_range(spec->lsel.range, m->rows);
    }
    if (!err) {
	/* check spec->rsel against cols */
	if (rt == SEL_MATRIX) {
	    err = bad_sel_vector(spec->rsel.m, m->cols);
	} else if (rt == SEL_RANGE || rt == SEL_ELEMENT) {
	    err = bad_sel_range(spec->rsel.range, m->cols);
	}
    }

    return err;
}

#if CONTIG_DEBUG
# define IN_USERMAT
# include "mspec_debug.c"
#endif

/* When we come here we've already screened out DIAG
   selection and single-element selection.
*/

static int submatrix_contig (matrix_subspec *spec, const gretl_matrix *m)
{
    int contig = 0;

    if (spec->ltype == SEL_MATRIX || spec->rtype == SEL_MATRIX ||
	spec->ltype == SEL_EXCL || spec->rtype == SEL_EXCL) {
	; /* cannot treat as contig */
    } else if (m->rows == 1 || m->cols == 1) {
	/* must be contig */
	contig = 1;
    } else if (singleton_right_range(spec)) {
	/* single column selected */
	contig = 1;
    } else if ((spec->rtype == SEL_ALL || spec->rtype == SEL_RANGE) &&
	       spec->ltype == SEL_ALL) {
	/* range of columns, all rows */
	contig = 1;
    }

    return contig;
}

/* Determine the offset from the start of m->val, and the
   number of elements, in a matrix slice that has been found
   to consist of contiguous data.
*/

static int get_offset_nelem (matrix_subspec *spec,
			     const gretl_matrix *m,
			     int *poff, int *pn)
{
    SelType stype = spec->ltype;
    union msel *sel = &spec->lsel;
    int offset = 0, nelem = 0;

    if (m->cols == 1) {
	/* column vector */
	if (stype == SEL_ALL) {
	    nelem = m->rows;
	} else {
	    int rmax = rowmax(sel, m);

	    offset = sel->range[0] - 1;
	    nelem = rmax - sel->range[0] + 1;
	}
    } else if (m->rows == 1) {
	/* row vector */
	if (spec->rtype != SEL_NULL) {
	    stype = spec->rtype;
	    sel = &spec->rsel;
	}
	if (stype == SEL_ALL) {
	    nelem = m->cols;
	} else {
	    int cmax = colmax(sel, m);

	    offset = sel->range[0] - 1;
	    nelem = cmax - sel->range[0] + 1;
	}
    } else if (singleton_right_range(spec)) {
	/* a single matrix column */
	if (stype == SEL_ALL) {
	    nelem = m->rows;
	} else {
	    int rmax = rowmax(sel, m);

	    offset = sel->range[0] - 1;
	    nelem = rmax - sel->range[0] + 1;
	}
	offset += m->rows * (spec->rsel.range[0] - 1);
    } else {
	/* multiple columns, all rows */
	if (spec->rtype == SEL_ALL) {
	    nelem = m->rows * m->cols;
	} else {
	    int cmax, nc;

	    sel = &spec->rsel;
	    cmax = colmax(sel, m);
	    nc = cmax - sel->range[0] + 1;
	    offset = m->rows * (sel->range[0] - 1);
	    nelem = m->rows * nc;
	}
    }

    *poff = offset;
    *pn = nelem;

    return 0;
}

/* to make clear that spec->lsel pertains to columns, move
   it to the right position */

static void commute_selectors (matrix_subspec *spec)
{
    spec->rtype = spec->ltype;
    spec->rsel = spec->lsel;
    memset(&spec->lsel, 0, sizeof spec->lsel);
    spec->ltype = SEL_ALL;
}

/* Catch the case of an implicit column or row specification for
   a sub-matrix of an (n x 1) or (1 x m) matrix; also catch the
   error of giving just one row/col spec for a matrix that has
   more than one row and more than one column.
*/

int check_matrix_subspec (matrix_subspec *spec, const gretl_matrix *m)
{
    int veclen, err = 0;

    if (spec->checked) {
	return 0;
    } else if (spec->ltype == SEL_REAL && !m->is_complex) {
	spec->ltype = spec->rtype = SEL_ALL;
	goto check_done;
    } else if (is_sel_dummy(spec->ltype)) {
	/* SEL_DIAG, SEL_UPPER, SEL_LOWER, SEL_IMAG:
	   nothing to be done here?
	*/
	goto check_done;
    }

#if CONTIG_DEBUG
    subspec_debug_print(spec, m);
#endif

    if (spec->rtype == SEL_NULL && m->rows == 1 && m->cols > 1) {
	/* row vector: transfer spec to column dimension */
	commute_selectors(spec);
    }

    veclen = gretl_vector_get_length(m);

    if (spec->rtype == SEL_NULL && m->rows > 1 && m->cols > 1) {
	gretl_errmsg_set(_("Ambiguous matrix index"));
	return E_DATA;
    }

    if (spec_is_single(spec)) {
	int k = spec_single_val(spec);

	err = bad_sel_single(&k, veclen);
	if (!err) {
	    spec->ltype = spec->rtype = SEL_ELEMENT;
	    mspec_set_offset(spec, k - 1);
	    mspec_set_n_elem(spec, 1);
	}
	/* nothing more to do */
	goto check_done;
    }

    err = spec_check_dimensions(spec, m);
    if (err) {
	return err;
    }

    if (m->rows == 0 || m->cols == 0) {
	/* nothing more to do */
	goto check_done;
    }

    if (lhs_is_scalar(spec, m) && rhs_is_scalar(spec, m)) {
	/* we're looking at just one element */
	set_element_index(spec, m);
	spec->ltype = spec->rtype = SEL_ELEMENT;
	goto check_done;
    }

    if (submatrix_contig(spec, m)) {
	int offset, nelem;

	get_offset_nelem(spec, m, &offset, &nelem);
	spec->ltype = SEL_CONTIG;
	mspec_set_offset(spec, offset);
	mspec_set_n_elem(spec, nelem);
	if (offset < 0 || nelem <= 0) {
	    fprintf(stderr, "*** offset = %d, nelem = %d (m: %dx%d) ***\n",
		    offset, nelem, m->rows, m->cols);
            gretl_errmsg_set(_("Invalid submatrix specification"));
	    return E_DATA;
	}
    }

 check_done:

    if (!err) {
	spec->checked = 1;
    }

    return err;
}

const char *mspec_get_string (matrix_subspec *spec, int i)
{
    if (spec == NULL || i < 0 || i > 1) {
	return NULL;
    } else if (i == 0) {
	return spec->ltype != SEL_STR ? NULL : spec->lsel.str;
    } else {
	return spec->rtype != SEL_STR ? NULL : spec->rsel.str;
    }
}

static int get_slices (matrix_subspec *spec,
		       const gretl_matrix *M)
{
    int err = 0;

    spec->rslice = mspec_make_list(spec->ltype, &spec->lsel,
				   M->rows, &err);
    if (!err) {
	spec->cslice = mspec_make_list(spec->rtype, &spec->rsel,
				       M->cols, &err);
    }

#if MDEBUG
    if (err) {
	fprintf(stderr, "usermat: get_slices, err=%d\n", err);
    }
#endif

    return err;
}

int assign_scalar_to_submatrix (gretl_matrix *M,
				const gretl_matrix *S,
				double x,
				matrix_subspec *spec)
{
    double complex z = NADBL;
    int mr = M->rows;
    int mc = M->cols;
    int i, err = 0;

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (M->is_complex) {
	z = (S != NULL)? S->z[0] : x;
    }

    if (spec->ltype == SEL_CONTIG) {
	int ini = mspec_get_offset(spec);
	int fin = ini + mspec_get_n_elem(spec);

	for (i=ini; i<fin; i++) {
	    if (M->is_complex) {
		M->z[i] = z;
	    } else {
		M->val[i] = x;
	    }
	}
	return 0;
    }

    if (is_sel_dummy(spec->ltype)) {
	if (M->is_complex) {
	    err = gretl_matrix_set_part(M, S, x, spec->ltype);
	} else {
	    err = gretl_matrix_set_part(M, NULL, x, spec->ltype);
	}
	return err; /* we're done */
    }

    if (spec->rslice == NULL && spec->cslice == NULL) {
	/* parse @spec into lists of affected rows and columns */
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
		if (M->is_complex) {
		    gretl_cmatrix_set(M, mi, mj, z);
		} else {
		    gretl_matrix_set(M, mi, mj, x);
		}
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

static int contig_cols (matrix_subspec *spec,
			const gretl_matrix *m)
{
    if (spec->rtype == SEL_ALL) {
	return m->cols;
    } else if (spec->rtype == SEL_RANGE) {
	int r1 = spec->rsel.range[1];
	int cmax = r1 == MSEL_MAX ? m->cols : r1;

	return cmax - spec->rsel.range[0] + 1;
    } else {
	return 1;
    }
}

static void transcribe_cols_8 (gretl_matrix *M,
			       const gretl_matrix *S,
			       const int *cslice,
			       int sscalar)
{
    int i, j, mcol, nr = M->rows;

    if (sscalar) {
	/* write RHS scalar into all rows of each selected col */
	for (j=1; j<=cslice[0]; j++) {
	    mcol = cslice[j] - 1;
	    for (i=0; i<nr; i++) {
		gretl_matrix_set(M, i, mcol, S->val[0]);
	    }
	}
    } else {
	/* zap from (cols of) S into selected cols of M */
	double *xtarg, *xsrc = S->val;
	size_t colsize = nr * sizeof *M->val;

	for (j=1; j<=cslice[0]; j++) {
	    mcol = cslice[j] - 1;
	    xtarg = M->val + mcol * nr;
	    memcpy(xtarg, xsrc, colsize);
	    xsrc += nr;
	}
    }
}

static void transcribe_cols_16 (gretl_matrix *M,
			       const gretl_matrix *S,
			       const int *cslice,
			       int sscalar)
{
    int i, j, mcol, nr = M->rows;

    if (sscalar) {
	/* write RHS scalar into all rows of each selected col */
	double complex z = S->is_complex ? S->z[0] : S->val[0];

	for (j=1; j<=cslice[0]; j++) {
	    mcol = cslice[j] - 1;
	    for (i=0; i<nr; i++) {
		gretl_cmatrix_set(M, i, mcol, z);
	    }
	}
    } else {
	/* zap from (cols of) S into selected cols of M */
	double complex *ztarg, *zsrc = S->z;
	double *xsrc = S->val;
	size_t colsize = nr * sizeof *M->z;

	for (j=1; j<=cslice[0]; j++) {
	    mcol = cslice[j] - 1;
	    ztarg = M->z + mcol * nr;
	    if (!S->is_complex) {
		for (i=0; i<nr; i++) {
		    ztarg[i] = xsrc[i];
		}
		xsrc += nr;
	    } else {
		memcpy(ztarg, zsrc, colsize);
		zsrc += nr;
	    }
	}
    }
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
    int i, j;
    int err = 0;

    if (!M->is_complex && S->is_complex) {
	fputs("matrix_replace_submatrix: M is real but S is complex\n", stderr);
	return E_MIXED;
    }

#if MDEBUG
    fprintf(stderr, "\nmatrix_replace_submatrix\n");
#endif

    if (spec == NULL) {
	fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n");
	return E_DATA;
    }

    if (is_sel_dummy(spec->ltype)) {
	return gretl_matrix_set_part(M, S, 0, spec->ltype);
    }

    if (spec->ltype == SEL_CONTIG) {
	int ini = mspec_get_offset(spec);
	int n = mspec_get_n_elem(spec);
	int ccols = contig_cols(spec, M);

	if (M->is_complex && S->rows == 1 && S->cols == 1) {
	    for (i=0; i<n; i++) {
		M->z[ini + i] = S->z[0];
	    }
	} else if (S->rows * S->cols != n) {
	    err = E_NONCONF;
	} else if (ccols > 1 && M->rows != S->rows) {
	    err = E_NONCONF;
	} else if (M->is_complex && !S->is_complex) {
	    /* can't use memcpy here! */
	    for (i=0; i<n; i++) {
		M->z[ini + i] = S->val[i];
	    }
	} else if (M->is_complex) {
	    memcpy(M->z + ini, S->z, n * sizeof *M->z);
	} else {
	    memcpy(M->val + ini, S->val, n * sizeof *M->val);
	}
	return err;
    }

    if (sr > mr || sc > mc) {
	/* the replacement matrix won't fit into M */
	fprintf(stderr, "matrix_replace_submatrix: target is %d x %d but "
		"replacement part is %d x %d\n", mr, mc, sr, sc);
	return E_NONCONF;
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
	/* the replacement is a scalar */
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
	if (M->is_complex) {
	    transcribe_cols_16(M, S, spec->cslice, sscalar);
	} else {
	    transcribe_cols_8(M, S, spec->cslice, sscalar);
	}
    } else if (!err) {
	/* the general case, no special shortcuts */
	double complex z = 0;
	double x = 0;
	int l, k = 0;
	int mi, mj;

	if (sscalar && S->is_complex) {
	    z = S->z[0];
	} else if (sscalar) {
	    x = S->val[0];
	    z = x; /* in case M is complex */
	}

	for (j=0; j<sc; j++) {
	    mj = (spec->cslice == NULL)? k++ : spec->cslice[j+1] - 1;
	    l = 0;
	    for (i=0; i<sr; i++) {
		mi = (spec->rslice == NULL)? l++ : spec->rslice[i+1] - 1;
		if (!sscalar) {
		    if (S->is_complex) {
			z = gretl_cmatrix_get(S, i, j);
		    } else {
			x = gretl_matrix_get(S, i, j);
		    }
		}
		if (M->is_complex) {
		    gretl_cmatrix_set(M, mi, mj, z);
		} else {
		    gretl_matrix_set(M, mi, mj, x);
		}
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
    gretl_matrix *S = NULL;
    int r = -1;
    int c = -1;

    if (M == NULL || spec == NULL) {
	*err = E_DATA;
	return NULL;
    }

#if MDEBUG
    fprintf(stderr, "matrix_get_submatrix, M = %d x %d\n", M->rows, M->cols);
    fprintf(stderr, " ltype: %d, rtype %d, prechecked %d\n", spec->ltype,
            spec->rtype, prechecked);
#endif

    if (!prechecked) {
	*err = check_matrix_subspec(spec, M);
	if (*err) {
	    return NULL;
	}
    }

    if (spec->ltype == SEL_DIAG) {
	return gretl_matrix_get_diagonal(M, err);
    } else if (spec->ltype == SEL_UPPER || spec->ltype == SEL_LOWER) {
	int upper = (spec->ltype == SEL_UPPER);

	return gretl_matrix_get_triangle(M, upper, err);
    } else if (spec->ltype == SEL_REAL || spec->ltype == SEL_IMAG) {
	int im = (spec->ltype == SEL_IMAG);

	return gretl_cmatrix_extract(M, im, err);
    } else if (spec->ltype == SEL_CONTIG) {
	return matrix_get_chunk(M, spec, err);
    } else if (spec->ltype == SEL_ELEMENT) {
	int i = mspec_get_element(spec);

	if (M->is_complex) {
	    S = gretl_cmatrix_from_scalar(M->z[i], err);
	} else {
	    S = gretl_matrix_from_scalar(M->val[i]);
	}
	if (S == NULL) {
	    *err = E_ALLOC;
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
#endif

    if (M->rows == 1 && M->cols == 1) {
        if (spec->ltype == SEL_EXCL && spec->rtype == SEL_NULL) {
            r = c = 0;
        } else if (spec->ltype == SEL_EXCL && spec->rtype == SEL_EXCL) {
            r = c = 0;
        } else if (spec->ltype == SEL_EXCL && spec->rtype == SEL_ALL) {
            r = 0; c = 1;
        } else if (spec->ltype == SEL_ALL && spec->rtype == SEL_EXCL) {
            r = 1; c = 0;
        }
    }
    if (r < 0 && c < 0) {
        /* dimensions not yet determined */
        r = (spec->rslice == NULL)? M->rows : spec->rslice[0];
        c = (spec->cslice == NULL)? M->cols : spec->cslice[0];
    }

    S = gretl_matching_matrix_new(r, c, M);

    if (S == NULL) {
	*err = E_ALLOC;
    } else if (r * c == 0) {
        /* empty matrix: nothing more to be done */
        return S;
    } else {
	int j, mj;

	if (spec->cslice != NULL && spec->rslice == NULL) {
	    /* copying entire columns */
	    if (M->is_complex) {
		double complex *dest = S->z;
		size_t csize = r * sizeof *dest;

		for (j=0; j<c; j++) {
		    mj = spec->cslice[j+1] - 1;
		    memcpy(dest, M->z + mj * r, csize);
		    dest += r;
		}
	    } else {
		double *dest = S->val;
		size_t csize = r * sizeof *dest;

		for (j=0; j<c; j++) {
		    mj = spec->cslice[j+1] - 1;
		    memcpy(dest, M->val + mj * r, csize);
		    dest += r;
		}
	    }
	} else {
	    int i, mi, l, k = 0;
	    double complex z;
	    double x;

	    for (j=0; j<c; j++) {
		mj = (spec->cslice == NULL)? k++ : spec->cslice[j+1] - 1;
		l = 0;
		for (i=0; i<r; i++) {
		    mi = (spec->rslice == NULL)? l++ : spec->rslice[i+1] - 1;
		    if (M->is_complex) {
			z = gretl_cmatrix_get(M, mi, mj);
			gretl_cmatrix_set(S, i, j, z);
		    } else {
			x = gretl_matrix_get(M, mi, mj);
			gretl_matrix_set(S, i, j, x);
		    }
		}
	    }
	}
    }

    if (S != NULL) {
	/* try transcribing metadata on @M if applicable */
	if (S->rows == M->rows && gretl_matrix_is_dated(M)) {
	    gretl_matrix_transcribe_obs_info(S, M);
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
	*err = E_BOUNDS;
    } else {
	x = M->val[i];
    }

    return x;
}

/* Copy a contiguous chunk of data out of @M, the offset and
   extent of which are given by @spec.
*/

gretl_matrix *matrix_get_chunk (const gretl_matrix *M,
				matrix_subspec *spec,
				int *err)
{
    int offset = spec->lsel.range[0];
    int nelem = spec->lsel.range[1];
    gretl_matrix *ret;
    int rows, cols;

    if (offset < 0) {
	fprintf(stderr, "matrix_get_chunk: offset = %d\n", offset);
	*err = E_DATA;
	return NULL;
    }

    if (M->cols > 1 && M->rows > 1) {
	cols = contig_cols(spec, M);
	rows = nelem / cols;
    } else if (M->rows == 1) {
	rows = 1;
	cols = nelem;
    } else {
	rows = nelem;
	cols = 1;
    }

    ret = gretl_matching_matrix_new(rows, cols, M);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	size_t sz;

	if (M->is_complex) {
	    sz = nelem * sizeof *M->z;
	    memcpy(ret->z, M->z + offset, sz);
	} else {
	    sz = nelem * sizeof *M->val;
	    memcpy(ret->val, M->val + offset, sz);
	}
	if (M->rows > 1 && rows == M->rows && offset % rows == 0 &&
	    gretl_matrix_is_dated(M)) {
	    gretl_matrix_transcribe_obs_info(ret, M);
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

static int real_umatrix_set_names (gretl_matrix *M,
				   char **S,
				   int byrow)
{
    int n = byrow ? M->rows : M->cols;
    int err = 0;

    if (S == NULL) {
	if (byrow) {
	    gretl_matrix_set_rownames(M, NULL);
	} else {
	    gretl_matrix_set_colnames(M, NULL);
	}
    } else {
	int i;

	for (i=0; i<n && !err; i++) {
	    if (S[i] == NULL || S[i][0] == '\0') {
		gretl_errmsg_sprintf("Missing string in %s",
				     byrow? "rnameset" : "cnameset");
		err = E_INVARG;
	    }
	}
	if (!err) {
	    if (byrow) {
		gretl_matrix_set_rownames(M, S);
	    } else {
		gretl_matrix_set_colnames(M, S);
	    }
	}
    }

    if (err) {
	strings_array_free(S, n);
    }

    return err;
}

static int n_names_check (int n, int ns, int byrow)
{
    if (ns != n) {
	gretl_errmsg_sprintf("%s: got %d names but matrix has %d %s",
			     byrow ? "rnameset" : "cnameset",
			     ns, n, byrow ? "row(s)" : "column(s)");
	return E_INVARG;
    } else {
	return 0;
    }
}

int umatrix_set_names_from_string (gretl_matrix *M,
				   const char *s,
				   int byrow)
{
    char **S = NULL;
    int n, ns = 0;
    int err = 0;

    n = byrow ? M->rows : M->cols;

    if (s != NULL && *s != '\0') {
	S = gretl_string_split(s, &ns, " \n\t");

	if (S == NULL) {
	    err = E_ALLOC;
	} else if (ns == 1 && n > 1) {
	    strings_array_free(S, ns);
	    S = expand_names(s, n, &err);
	} else {
	    err = n_names_check(n, ns, byrow);
	    if (err) {
		strings_array_free(S, ns);
	    }
	}
    }

    if (!err) {
	err = real_umatrix_set_names(M, S, byrow);
    }

    return err;
}

int umatrix_set_names_from_array (gretl_matrix *M,
				  void *data,
				  int byrow)
{
    gretl_array *A = data;
    char **S = NULL;
    int n, ns = 0;
    int err = 0;

    n = byrow ? M->rows : M->cols;

    if (A != NULL && gretl_array_get_length(A) > 0) {
	char **AS = gretl_array_get_strings(A, &ns);

	err = n_names_check(n, ns, byrow);
	if (!err) {
	    S = strings_array_dup(AS, n);
	    if (S == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	err = real_umatrix_set_names(M, S, byrow);
    }

    return err;
}

int umatrix_set_names_from_list (gretl_matrix *M,
				 const int *list,
				 const DATASET *dset,
				 int byrow)
{
    char **S = NULL;
    int n, err = 0;

    n = byrow ? M->rows : M->cols;

    if (list != NULL && list[0] > 0) {
	int i;

	err = n_names_check(n, list[0], byrow);
	if (!err) {
	    S = strings_array_new(n);
	    if (S == NULL) {
		err = E_ALLOC;
	    }
	}
	for (i=0; i<n && !err; i++) {
	    S[i] = gretl_strndup(dset->varname[list[i+1]], 12);
	    if (S[i] == NULL) {
		err = E_ALLOC;
	    }
	}
	if (err) {
	    strings_array_free(S, n);
	}
    }

    if (!err) {
	err = real_umatrix_set_names(M, S, byrow);
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
			     int ldet, int *err)
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
	if (ldet) {
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

static void maybe_replace_content (gretl_matrix *targ,
				   gretl_matrix *tmp,
				   int err)
{
    if (!err) {
	gretl_matrix_replace_content(targ, tmp);
    }
    gretl_matrix_free(tmp);
}

int matrix_invert_in_place (gretl_matrix *m)
{
    gretl_matrix *tmp = gretl_matrix_copy(m);
    int err = 0;

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_invert_matrix(tmp);
	maybe_replace_content(m, tmp, err);
    }

    return err;
}

int matrix_cholesky_in_place (gretl_matrix *m)
{
    gretl_matrix *tmp = gretl_matrix_copy(m);
    int err = 0;

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_cholesky_decomp(tmp);
	maybe_replace_content(m, tmp, err);
    }

    return err;
}

int matrix_transpose_in_place (gretl_matrix *m)
{
    gretl_matrix *tmp = gretl_matrix_copy_transpose(m);
    int err = 0;

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	maybe_replace_content(m, tmp, 0);
    }

    return err;
}

int matrix_XTX_in_place (gretl_matrix *m)
{
    gretl_matrix *tmp = gretl_matrix_alloc(m->cols, m->cols);
    int err;

    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	err = gretl_matrix_multiply_mod(m, GRETL_MOD_TRANSPOSE,
					m, GRETL_MOD_NONE,
					tmp, GRETL_MOD_NONE);
	maybe_replace_content(m, tmp, err);
    }

    return err;
}

gretl_matrix *user_matrix_vec (const gretl_matrix *m, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else {
	R = gretl_matching_matrix_new(m->rows * m->cols, 1, m);
	if (R != NULL) {
	    gretl_matrix_vectorize(R, m);
	}
    }

    if (R == NULL) {
	*err = E_ALLOC;
    }

    return R;
}

gretl_matrix *user_matrix_vech (const gretl_matrix *m,
				int omit_diag, int *err)
{
    gretl_matrix *R = NULL;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else if (m->rows != m->cols) {
	*err = E_NONCONF;
    } else {
	int k, n = m->rows;

	if (omit_diag) {
	    k = n * (n - 1) / 2;
	} else {
	    k = n * (n + 1) / 2;
	}
	R = gretl_matching_matrix_new(k, 1, m);
	if (R != NULL) {
	    if (omit_diag) {
		*err = gretl_matrix_vectorize_h_skip(R, m);
	    } else {
		*err = gretl_matrix_vectorize_h(R, m);
	    }
	}
    }

    if (R == NULL && !*err) {
	*err = E_ALLOC;
    }

    return R;
}

gretl_matrix *user_matrix_unvech (const gretl_matrix *m,
				  double diag, int *err)
{
    gretl_matrix *R = NULL;
    int k;

    if (gretl_is_null_matrix(m)) {
	R = gretl_null_matrix_new();
    } else if ((k = gretl_vector_get_length(m)) == 0) {
	*err = E_NONCONF;
    } else {
	int n;

	if (na(diag)) {
	    n = (int) ((sqrt(1.0 + 8.0 * m->rows) - 1.0) / 2.0);
	} else {
	    n = (int) ((sqrt(1.0 + 8.0 * m->rows) + 1.0) / 2.0);
	}
	R = gretl_matching_matrix_new(n, n, m);
	if (R != NULL) {
	    if (na(diag)) {
		*err = gretl_matrix_unvectorize_h(R, m);
	    } else {
		*err = gretl_matrix_unvectorize_h_diag(R, m, diag);
	    }
	}
    }

    if (R == NULL && !*err) {
	*err = E_ALLOC;
    }

    return R;
}

static int real_user_matrix_QR_decomp (const gretl_matrix *m,
				       gretl_matrix **Q,
				       gretl_matrix **R,
				       gretl_matrix **P)
{
    int mc = gretl_matrix_cols(m);
    int err = 0;

    *Q = gretl_matrix_copy(m);

    if (*Q == NULL) {
	err = E_ALLOC;
    } else {
	if (R != NULL) {
	    *R = gretl_matrix_alloc(mc, mc);
	    if (*R == NULL) {
		err = E_ALLOC;
	    }
	}
	if (!err && P != NULL) {
	    *P = gretl_matrix_alloc(1, mc);
	    if (*P == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	if (P != NULL) {
	    err = gretl_matrix_QR_pivot_decomp(*Q, (R == NULL)? NULL : *R, *P);
	} else {
	    err = gretl_matrix_QR_decomp(*Q, (R == NULL)? NULL : *R);
	}
    }

    if (err) {
	gretl_errmsg_set(_("Matrix decomposition failed"));
	gretl_matrix_free(*Q);
	*Q = NULL;
	if (R != NULL) {
	    gretl_matrix_free(*R);
	    *R = NULL;
	}
	if (P != NULL) {
	    gretl_matrix_free(*P);
	    *P = NULL;
	}
    }

    return err;
}

gretl_matrix *user_matrix_QR_decomp (const gretl_matrix *m,
				     gretl_matrix *R,
				     gretl_matrix *P,
				     int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *Rtmp = NULL;
    gretl_matrix *Ptmp = NULL;
    gretl_matrix **pR;
    gretl_matrix **pP;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    if (m->cols > m->rows) {
	gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns"));
	*err = E_NONCONF;
	return NULL;
    }

    pR = (R != NULL)? &Rtmp : NULL;
    pP = (P != NULL)? &Ptmp : NULL;

    *err = real_user_matrix_QR_decomp(m, &Q, pR, pP);

    if (Rtmp != NULL) {
	maybe_replace_content(R, Rtmp, *err);
    }
    if (Ptmp != NULL) {
	maybe_replace_content(P, Ptmp, *err);
    }

    return Q;
}

gretl_matrix *user_matrix_SVD (const gretl_matrix *m,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       int *err)
{
    gretl_matrix *S = NULL;
    gretl_matrix *Utmp = NULL;
    gretl_matrix *Vtmp = NULL;
    gretl_matrix **pU, **pV;

    if (gretl_is_null_matrix(m)) {
	*err = E_DATA;
	return NULL;
    }

    pU = (U != NULL)? &Utmp : NULL;
    pV = (V != NULL)? &Vtmp : NULL;

    if (m->is_complex) {
	*err = gretl_cmatrix_SVD(m, pU, &S, pV, 0);
    } else {
	*err = gretl_matrix_SVD(m, pU, &S, pV, 0);
    }
    if (!*err) {
	if (Utmp != NULL) {
	    maybe_replace_content(U, Utmp, *err);
	}
	if (Vtmp != NULL) {
	    maybe_replace_content(V, Vtmp, *err);
	}
    }

    return S;
}

static gretl_matrix *null_OLS (const gretl_matrix *Y,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       int *err)
{
    gretl_matrix *B = NULL;

    if (U != NULL) {
	/* residuals = Y */
	gretl_matrix *Ycpy;

	Ycpy = gretl_matrix_copy(Y);
	if (Ycpy == NULL) {
	    *err = E_ALLOC;
	} else {
	    gretl_matrix_replace_content(U, Ycpy);
	    gretl_matrix_free(Ycpy);
	}
    }

    if (!*err && V != NULL) {
	/* variance is empty */
	gretl_matrix *nullV;

	nullV = gretl_null_matrix_new();
	gretl_matrix_replace_content(V, nullV);
	gretl_matrix_free(nullV);
    }

    if (!*err) {
	B = gretl_matrix_alloc(0, Y->cols);
    }

    return B;
}

gretl_matrix *user_matrix_ols (const gretl_matrix *Y,
			       const gretl_matrix *X,
			       gretl_matrix *U,
			       gretl_matrix *V,
			       gretlopt opt,
			       int *err)
{
    gretl_matrix *B = NULL;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || X == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_is_complex(Y) || gretl_is_complex(X) ||
	gretl_is_complex(U) || gretl_is_complex(V)) {
	fprintf(stderr, "E_CMPLX in user_matrix_ols\n");
	*err = E_CMPLX;
	return NULL;
    }

    T = Y->rows;
    k = X->cols;
    g = Y->cols;

    if (X->rows != T) {
	*err = E_NONCONF;
	return NULL;
    }

    if (X->cols == 0) {
	/* handle the null model case */
	return null_OLS(Y, U, V, err);
    }

    if (g > 1 && (opt & OPT_M)) {
	/* multiple precision: we accept only one y var */
	*err = E_DATA;
	return NULL;
    }

    if (U != NULL) {
	*err = gretl_matrix_realloc(U, T, g);
	if (*err) {
	    return NULL;
	}
    }

    if (V != NULL && g == 1) {
	/* single regressand case */
	*err = gretl_matrix_realloc(V, k, k);
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	if (g == 1) {
	    /* single regressand */
	    double s2 = 0;
	    double *ps2 = (V != NULL)? &s2 : NULL;

	    if (opt & OPT_M) {
		/* use multiple precision */
		*err = gretl_matrix_mp_ols(Y, X, B, V, U, ps2);
	    } else {
		*err = gretl_matrix_ols(Y, X, B, V, U, ps2);
	    }
	} else {
	    /* multiple regressands: @V will actually be (X'X)^{-1} */
	    gretl_matrix *Vtmp = NULL;
	    gretl_matrix **Vp = (V != NULL)? &Vtmp : NULL;

	    *err = gretl_matrix_multi_ols(Y, X, B, U, Vp);
	    if (Vtmp != NULL) {
		maybe_replace_content(V, Vtmp, *err);
	    }
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
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
    int g, k, T;

    if (gretl_is_null_matrix(Y) || gretl_is_null_matrix(X)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_is_complex(Y) || gretl_is_complex(X) ||
	gretl_is_complex(R) || gretl_is_complex(Q)) {
	fprintf(stderr, "E_CMPLX in user_matrix_rls\n");
	*err = E_CMPLX;
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
	*err = gretl_matrix_realloc(U, T, g);
	if (*err) {
	    return NULL;
	}
    }

    if (!*err) {
	B = gretl_matrix_alloc(k, g);
	if (B == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* note: &V will actually be M (X'X)^{-1} M' */
	gretl_matrix *Vtmp = NULL;
	gretl_matrix **Vp = (V != NULL)? &Vtmp : NULL;

	*err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B,
						 U, Vp);
	if (Vtmp != NULL) {
	    maybe_replace_content(V, Vtmp, *err);
	}
    }

    if (*err) {
	gretl_matrix_free(B);
	B = NULL;
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
    int n, m, npar;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(C)) {
	*err = E_DATA;
	return NULL;
    }

    n = A->rows;
    m = C->rows;
    npar = m + m + m*(m+1)/2;

    if (dP != NULL) {
	*err = gretl_matrix_realloc(dP, n, npar);
    }

    if (!*err) {
	P = gretl_GHK2(C, A, B, U, dP, err);
    }

    return P;
}

gretl_matrix *user_matrix_eigensym (const gretl_matrix *m,
				    gretl_matrix *R,
				    int *err)
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
	/* computing (right) eigenvectors */
	vecs = 1;
    }

    /* @m would be destroyed by this operation */
    C = gretl_matrix_copy(m);
    if (C == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	if (m->is_complex) {
	    E = gretl_zheev(C, vecs, err);
	} else {
	    E = gretl_symmetric_matrix_eigenvals(C, vecs, err);
	}
    }

    if (!*err && vecs) {
	maybe_replace_content(R, C, 0);
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

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
	*err = E_DATA;
	return NULL;
    }

    if (gretl_matrix_na_check(A) || gretl_matrix_na_check(B)) {
	*err = E_NAN;
	return NULL;
    }

    if (V != NULL) {
	*err = gretl_matrix_realloc(V, B->cols, A->rows);
    }

    if (!*err) {
	E = gretl_gensymm_eigenvals(A, B, V, err);
    }

    return E;
}

int gretl_matrix_set_part (gretl_matrix *targ,
			   const gretl_matrix *src,
			   double x, SelType sel)
{
    int err = 0;

    if (sel == SEL_DIAG) {
	if (targ->is_complex) {
	    err = gretl_cmatrix_set_diagonal(targ, src, x);
	} else {
	    err = gretl_matrix_set_diagonal(targ, src, x);
	}
    } else if (sel == SEL_LOWER || sel == SEL_UPPER) {
	int upper = (sel == SEL_UPPER);

	if (targ->is_complex) {
	    err = gretl_cmatrix_set_triangle(targ, src, x, upper);
	} else {
	    err = gretl_matrix_set_triangle(targ, src, x, upper);
	}
    } else if (sel == SEL_REAL || sel == SEL_IMAG) {
	if (targ->is_complex) {
	    int im = (sel == SEL_IMAG);

	    err = gretl_cmatrix_set_part(targ, src, x, im);
	} else {
	    err = E_TYPES;
	}
    }

    return err;
}
