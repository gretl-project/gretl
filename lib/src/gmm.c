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
#include "nlspec.h"
#include "libset.h"
#include "qr_estimate.h"
#include "gretl_bfgs.h"
#include "uservar.h"

#include "../../minpack/minpack.h"

#include <errno.h>

#define GMM_DEBUG 0

typedef struct colsrc_ colsrc;
typedef struct hac_info_ hac_info;

struct colsrc_ {
    int v;                 /* ID number of variable in column */
    char mname[VNAMELEN];  /* or name of source matrix */
    int j;                 /* and column within matrix */
};

struct hac_info_ {
    int kern;              /* kernel type */
    int h;                 /* integer bandwidth (Bartlett, Parzen) */
    double bt;             /* QS bandwidth */
    int whiten;            /* pre-whiten (1) or not (0) */
};

struct ocset_ {
    gretl_matrix *e;      /* GMM residual, or LHS term in O.C. */
    gretl_matrix *Z;      /* instruments, or RHS terms in O.C. */
    gretl_matrix *W;      /* matrix of weights */
    gretl_matrix *tmp;    /* holds columnwise product of e and Z */
    gretl_matrix *sum;    /* holds column sums of tmp */
    gretl_matrix *S;      /* selector matrix for computing tmp */
    colsrc *ecols;        /* info on provenance of columns of 'e' */
    int noc;              /* total number of orthogonality conds. */
    int step;             /* number of estimation steps */
    int userwts;          /* got weights matrix from user (1/0) */
    double scale_wt;      /* auto-scaling for weights matrix  */
    hac_info hinfo;       /* HAC characteristics */
    char Wname[VNAMELEN]; /* name of weighting matrix */
    char **lnames;        /* names of LHS terms in O.C.s */
    char **rnames;        /* names of RHS terms in O.C.s */
    int n_names;          /* number of right/left name pairs */
};

#define using_HAC(s) (s->oc->hinfo.kern >= KERNEL_BARTLETT)

/* destructor apparatus for set of O.C. info */

void oc_set_destroy (ocset *oc)
{
    if (oc == NULL) {
	return;
    }

    gretl_matrix_free(oc->e);
    gretl_matrix_free(oc->Z);
    gretl_matrix_free(oc->tmp);
    gretl_matrix_free(oc->sum);
    gretl_matrix_free(oc->S);

    free(oc->ecols);

    if (oc->lnames != NULL) {
        strings_array_free(oc->lnames, oc->n_names);
    }
    if (oc->rnames != NULL) {
	strings_array_free(oc->rnames, oc->n_names);
    }
    if (!oc->userwts) {
	/* we used auto-generated weights */
	gretl_matrix_free(oc->W);
    }

    free(oc);
}

/* constructor for same */

static ocset *oc_set_new (void)
{
    ocset *oc = malloc(sizeof *oc);

    if (oc != NULL) {
	oc->e = NULL;
	oc->Z = NULL;
	oc->W = NULL;
	oc->tmp = NULL;
	oc->sum = NULL;
	oc->S = NULL;
	oc->ecols = NULL;
	oc->noc = 0;
	oc->step = 0;
	oc->userwts = 0;
	oc->scale_wt = 1;
	*oc->Wname = '\0';
	oc->lnames = NULL;
	oc->rnames = NULL;
	oc->n_names = 0;
    }

    return oc;
}

static int gmm_unkvar (const char *s)
{
    gretl_errmsg_sprintf(_("Unknown variable '%s'"), s);
    return E_UNKVAR;
}

/* Determine the type of an element in an "orthog" specification:
   series, matrix or list of series
*/

static int
oc_get_type (const char *name, const DATASET *dset, int *err)
{
    int *list = NULL;
    int j, v = current_series_index(dset, name);
    int ret = GRETL_TYPE_NONE;

#if GMM_DEBUG
    fprintf(stderr, "oc_get_type: looking at '%s' (v = %d, dset->v = %d)\n",
	    name, v, dset->v);
#endif

    if (v >= 0) {
	/* try for a series first */
	ret = GRETL_TYPE_SERIES;
    } else if (get_matrix_by_name(name) != NULL) {
	/* then a matrix */
	ret = GRETL_TYPE_MATRIX;
    } else if ((list = get_list_by_name(name)) != NULL) {
	/* failing that, a list */
	for (j=1; j<=list[0]; j++) {
	    if (list[j] < 0 || list[j] >= dset->v) {
		*err = E_DATA;
		return GRETL_TYPE_NONE;
	    }
	}
	ret = GRETL_TYPE_LIST;
    } else {
	*err = gmm_unkvar(name);
    }

    return ret;
}

/* See if column j in matrix b is present anywhere in matrix a.  If
   so, set *colno to the matching (0-based) column number in 'a' and
   return 1; if not, return 0.
*/

static int
col_present (const gretl_matrix *a, const gretl_matrix *b, int j,
	     int *colno)
{
    double xa, xb;
    int i, k, match;
    int ret = 0;

    for (k=0; k<a->cols; k++) {
	match = 1;
	for (i=0; i<a->rows; i++) {
	    xa = gretl_matrix_get(a, i, k);
	    xb = gretl_matrix_get(b, i, j);
	    if (xa != xb) {
		match = 0;
		break;
	    }
	}
	if (match) {
	    *colno = k;
	    ret = 1;
	    break;
	}
    }

#if GMM_DEBUG
    if (ret) {
	fprintf(stderr, "col_present: col %d of new matrix: "
		"matched old col %d\n", j, *colno);
    } else {
	fprintf(stderr, "col_present: col %d of new matrix: no match\n", j);
    }
#endif

    return ret;
}

#define Sidx(i,j,m)   ((j)*m+(i))

/* Expand the matrix that keeps track of which columns of e should
   be paired with which columns of Z for the computation of the
   GMM criterion function.  The Z matrix has already been expanded
   at this point.
*/

static int
expand_selector_matrix (nlspec *s, int sm, int sn, const char *mask)
{
    gretl_matrix *Z = s->oc->Z;
    gretl_matrix *S = s->oc->S;
    int ecols = s->oc->e->cols;
    int Zcols = Z->cols;
    double xij;
    int i, j, err = 0;

#if GMM_DEBUG
    fprintf(stderr, "expand_selector_matrix\n");
    gretl_matrix_print(S, "S, on entry");
#endif

    if (S == NULL) {
	/* starting from scratch */
	S = s->oc->S = gretl_unit_matrix_new(sm, sn);
	if (S == NULL) {
	    return E_ALLOC;
	}
    }

    err = gretl_matrix_realloc(S, ecols, Zcols);
    if (err) {
	return err;
    }

    /* transcribe original entries */
    for (j=sn-1; j>=0; j--) {
	for (i=sm-1; i>=0; i--) {
	    xij = S->val[Sidx(i,j,sm)];
	    gretl_matrix_set(S, i, j, xij);
	}
    }

    /* 0s in upper-right block */
    for (i=0; i<sm; i++) {
	for (j=sn; j<S->cols; j++) {
	    gretl_matrix_set(S, i, j, 0);
	}
    }

    /* 1s in lower-right block */
    for (i=sm; i<S->rows; i++) {
	for (j=sn; j<S->cols; j++) {
	    gretl_matrix_set(S, i, j, 1);
	}
    }

    /* mixed 0/1 in lower-left region */
    for (i=sm; i<S->rows; i++) {
	for (j=0; j<sn; j++) {
	    gretl_matrix_set(S, i, j, mask[j] ? 1 : 0);
	}
    }

#if GMM_DEBUG
    gretl_matrix_print(S, "S, on exit");
#endif

    return err;
}

/* Add zero or more columns to the 'Z' (instrument) matrix in the
   nlspec (GMM).  We add columns only if they are not already present.
*/

static int
add_new_cols_to_Z (nlspec *s, const gretl_matrix *M, int oldecols)
{
    gretl_matrix *Z = s->oc->Z;
    int oldZcols = Z->cols;
    char *mask;
    int j, k, n;
    int err = 0;

    if (Z->rows != M->rows) {
	return E_NONCONF;
    }

    n = Z->cols + M->cols;

#if GMM_DEBUG
    fprintf(stderr, "add_new_cols_to_Z: old Z cols = %d, M->cols = %d\n",
	    oldZcols, M->cols);
#endif

#if GMM_DEBUG > 1
    gretl_matrix_print(M, "Add matrix, M");
#endif

    mask = malloc(n);
    if (mask == NULL) {
	return E_ALLOC;
    }

    for (j=M->cols; j<n; j++) {
	mask[j] = 0;
    }

    n = M->cols;

    for (j=0; j<M->cols; j++) {
	if (col_present(Z, M, j, &k)) {
	    mask[j] = 0;
	    mask[M->cols + k] = 1;
	    n--;
	} else {
	    mask[j] = 1;
	}
    }

    if (n > 0) {
#if GMM_DEBUG
	fprintf(stderr, "add_new_cols_to_Z: adding %d new columns\n", n);
#endif
	err = gretl_matrix_inplace_colcat(Z, M, mask);
    }

    if (!err) {
	err = expand_selector_matrix(s, oldecols, oldZcols, mask + M->cols);
    }

    free(mask);

    return err;
}

/* expand the array containing info on the columns on the left
   of a set of GMM orthogonality conditions.
*/

static int add_oc_cols (nlspec *s, int current_cols, int add_cols)
{
    colsrc *cols;
    int err = 0;

    cols = realloc(s->oc->ecols, (current_cols + add_cols) * sizeof *cols);

    if (cols == NULL) {
	err = E_ALLOC;
    } else {
	int j, k;

	for (j=0; j<add_cols; j++) {
	    k = current_cols + j;
	    cols[k].v = 0;
	    cols[k].j = j;
	    cols[k].mname[0] = '\0';
	}

	s->oc->ecols = cols;
    }

    return err;
}

/* Record the source (ID number of series) for a column on
   the LHS of a set of GMM orthogonality conditions.
*/

static int push_col_source_series (nlspec *s, int v)
{
    int n = (s->oc->e != NULL)? s->oc->e->cols : 0;
    int err = add_oc_cols(s, n, 1);

    if (!err) {
	s->oc->ecols[n].v = v;
    }

    return err;
}

/* Record the source (named matrix) for a set of columns on
   the LHS of a set of GMM orthogonality conditions.
*/

static int push_col_source_matrix (nlspec *s, const char *mname)
{
    int j, n = (s->oc->e != NULL)? s->oc->e->cols : 0;
    gretl_matrix *m = get_matrix_by_name(mname);
    int err = add_oc_cols(s, n, m->cols);

    if (!err) {
	for (j=0; j<m->cols; j++) {
	    strcpy(s->oc->ecols[n+j].mname, mname);
	}
    }

    return err;
}

/* Record the source (list of series) for a set of columns on
   the LHS of a set of GMM orthogonality conditions.
*/

static int push_col_source_list (nlspec *s, const int *list)
{
    int j, n = (s->oc->e != NULL)? s->oc->e->cols : 0;
    int err = add_oc_cols(s, n, list[0]);

    if (!err) {
	for (j=0; j<list[0]; j++) {
	    s->oc->ecols[n+j].v = list[j+1];
	}
    }

    return err;
}

int nlspec_add_ivreg_oc (nlspec *s, int lhv, const int *rlist,
			 const double **Z)
{
    gretl_matrix *e = NULL;
    gretl_matrix *M = NULL;
    int i, k, t, v;
    double x;
    int err = 0;

    s->oc = oc_set_new();
    if (s->oc == NULL) {
	return E_ALLOC;
    }

    /* the left-hand side: series -> vector */

    e = gretl_column_vector_alloc(s->nobs);
    if (e == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_set_t1(e, s->t1);
	gretl_matrix_set_t2(e, s->t2);
	for (t=0; t<s->nobs; t++) {
	    x = Z[lhv][t + s->t1];
	    gretl_vector_set(e, t, x);
	}
	err = push_col_source_series(s, lhv);
    }

    if (err) {
	return err;
    }

    /* the right-hand side: list -> matrix */

    k = rlist[0];
    M = gretl_matrix_alloc(s->nobs, k);
    if (M == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_set_t1(M, s->t1);
	gretl_matrix_set_t2(M, s->t2);
	for (i=0; i<k && !err; i++) {
	    v = rlist[i + 1];
	    for (t=0; t<s->nobs; t++) {
		x = Z[v][t + s->t1];
		gretl_matrix_set(M, t, i, x);
	    }
	}
    }

    if (!err) {
	s->oc->e = e;
	s->oc->Z = M;
	s->oc->noc = k;
    } else {
	oc_set_destroy(s->oc);
	s->oc = NULL;
    }

    return err;
}

static int oc_add_matrices (nlspec *s, int ltype, const char *lname,
			    int rtype, const char *rname,
			    const DATASET *dset)
{
    gretl_matrix *e = NULL;
    gretl_matrix *M = NULL;
    int oldecols = 0;
    int i, k, t, v;
    double x;
    int err = 0;

    if (s->oc->e != NULL) {
	oldecols = s->oc->e->cols;
    }

    /* First work on the left-hand side: if we've been given a series,
       convert to a vector; if a list, convert to a matrix.
    */

#if GMM_DEBUG
    fprintf(stderr, "oc_add_matrices: lname = '%s'\n", lname);
#endif

    if (ltype == GRETL_TYPE_MATRIX) {
	e = get_matrix_copy_by_name(lname, &err);
	if (!err) {
	    err = push_col_source_matrix(s, lname);
	}
    } else if (ltype == GRETL_TYPE_LIST) {
	int *list = get_list_by_name(lname);

	k = list[0];
	e = gretl_matrix_alloc(s->nobs, k);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_set_t1(e, s->t1);
	    gretl_matrix_set_t2(e, s->t2);
	    for (i=0; i<k; i++) {
		v = list[i+1];
		for (t=0; t<s->nobs; t++) {
		    x = dset->Z[v][t + s->t1];
		    gretl_matrix_set(e, t, i, x);
		}
	    }
	    err = push_col_source_list(s, list);
	}
    } else if (ltype == GRETL_TYPE_SERIES) {
	v = series_index(dset, lname);
	e = gretl_column_vector_alloc(s->nobs);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_set_t1(e, s->t1);
	    gretl_matrix_set_t2(e, s->t2);
	    for (t=0; t<s->nobs; t++) {
		x = dset->Z[v][t + s->t1];
		gretl_vector_set(e, t, x);
	    }
	    err = push_col_source_series(s, v);
	}
    } else {
	err = gmm_unkvar(lname);
    }

    if (err) {
	return err;
    }

    k = 0;

    /* Now process the right-hand side: should be a matrix, a
       list, or a single series.
    */

#if GMM_DEBUG
    fprintf(stderr, "oc_add_matrices: rname = '%s'\n", rname);
#endif

    if (rtype == GRETL_TYPE_MATRIX) {
	M = get_matrix_copy_by_name(rname, &err);
	if (!err) {
	    k = gretl_matrix_cols(M);
	}
    } else if (rtype == GRETL_TYPE_LIST) {
	int *list = get_list_by_name(rname);

	k = list[0];
	M = gretl_matrix_alloc(s->nobs, k);
	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_set_t1(M, s->t1);
	    gretl_matrix_set_t2(M, s->t2);
	    for (i=0; i<k && !err; i++) {
		v = list[i + 1];
		for (t=0; t<s->nobs; t++) {
		    x = dset->Z[v][t + s->t1];
		    gretl_matrix_set(M, t, i, x);
		}
	    }
	}
    } else {
	/* name of single series */
	v = series_index(dset, rname);
	k = 1;
	M = gretl_column_vector_alloc(s->nobs);
	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_matrix_set_t1(M, s->t1);
	    gretl_matrix_set_t2(M, s->t2);
	    for (t=0; t<s->nobs; t++) {
		x = dset->Z[v][t + s->t1];
		gretl_vector_set(M, t, x);
	    }
	}
    }

    if (!err) {
	if (s->oc->e == NULL) {
	    /* this is the first set of O.C.s */
	    s->oc->e = e;
	    s->oc->Z = M;
#if GMM_DEBUG > 1
	    fprintf(stderr, "oc_add_matrices: added first entries\n");
	    gretl_matrix_print(s->oc->e, "e");
	    gretl_matrix_print(s->oc->Z, "Z");
#endif
	} else {
	    /* it's an additional set of O.C.s */
	    err = gretl_matrix_inplace_colcat(s->oc->e, e, NULL);
#if GMM_DEBUG > 1
	    fprintf(stderr, "oc_add_matrices: expanded e\n");
	    gretl_matrix_print(s->oc->e, "e");
#endif
	    /* don't leak temporary matrix */
	    gretl_matrix_free(e);

	    if (!err) {
		err = add_new_cols_to_Z(s, M, oldecols);
#if GMM_DEBUG > 1
		fprintf(stderr, "oc_add_matrices: expanded Z\n");
		gretl_matrix_print(s->oc->Z, "Z");
#endif
	    }
	    gretl_matrix_free(M);
	}
    }

    if (!err) {
	s->oc->noc += k * (s->oc->e->cols - oldecols); /* FIXME? */
    }

    return err;
}

/* construct a list from two or more series names */

static int get_oc_list (gchar **ps, GretlType *pt,
                        const DATASET *dset, int rhs)
{
    int *list = NULL;
    gchar *lname = NULL;
    int err = 0;

    list = gretl_list_from_varnames(*ps, dset, &err);

    if (!err) {
        if (rhs) {
            lname = g_strdup("oc_rhs_tmp___");
        } else {
            lname = g_strdup("oc_lhs_tmp___");
        }
        err = user_var_add(lname, GRETL_TYPE_LIST, list);
        user_var_privatize_by_name(lname);
    }

    if (!err) {
        g_free(*ps);
        *ps = lname;
        *pt = GRETL_TYPE_LIST;
    }

    return err;
}

/* Add a set of orthogonality conditions to a GMM specification.
   The source command takes the form:

   orthog { series | list | matrix } ; { series | list | matrix }
*/

int
nlspec_add_orthcond (nlspec *s, const char *str,
		     const DATASET *dset)
{
    const char *p;
    gchar *lhs = NULL;
    gchar *rhs = NULL;
    GretlType ltype = GRETL_TYPE_NONE;
    GretlType rtype = GRETL_TYPE_NONE;
    int nl, nr;
    int n, err = 0;

    if (s->ci != GMM) {
	return E_TYPES;
    }

    if (s->oc != NULL && s->oc->W != NULL) {
	gretl_errmsg_set(_("Orthogonality conditions must come before the "
			   "weights matrix"));
	return E_DATA;
    }

#if GMM_DEBUG
    fprintf(stderr, "nlspec_add_orthcond: line = '%s'\n", str);
#endif

    str += strspn(str, " ");

    p = strchr(str, ';');
    if (p == NULL) {
	return E_PARSE;
    }

    nl = p - str;
    nr = strlen(str) - nl - 1;
    if (nl == 0 || nr == 0) {
	return E_PARSE;
    }

    lhs = g_strstrip(g_strndup(str, nl));
    rhs = g_strstrip(g_strndup(p+1, nr));

    if (strchr(lhs, ' ') != NULL) {
        err = get_oc_list(&lhs, &ltype, dset, 0);
        if (err) {
            goto bailout;
        }
    }
    if (strchr(rhs, ' ') != NULL) {
        err = get_oc_list(&rhs, &rtype, dset, 1);
        if (err) {
            goto bailout;
        }
    }

    if (s->oc == NULL) {
	s->oc = oc_set_new();
	if (s->oc == NULL) {
	    return E_ALLOC;
	}
    }

    if (ltype == GRETL_TYPE_NONE) {
        ltype = oc_get_type(lhs, dset, &err);
    }

    if (!err && rtype == GRETL_TYPE_NONE) {
	rtype = oc_get_type(rhs, dset, &err);
    }

    if (!err) {
	err = oc_add_matrices(s, ltype, lhs, rtype, rhs, dset);
    }

    if (!err && gretl_in_gui_mode()) {
	n = s->oc->n_names;
	strings_array_add(&s->oc->lnames, &s->oc->n_names, lhs);
	strings_array_add(&s->oc->rnames, &n, rhs);
    }

    if (err) {
	nlspec_destroy_arrays(s);
	oc_set_destroy(s->oc);
	s->oc = NULL;
    }

 bailout:

    g_free(lhs);
    g_free(rhs);

    return err;
}

#define matrix_t1(m,t) ((m->t1 == 0 && m->t2 == 0)? t : m->t1)

static int gmm_matrix_resize (gretl_matrix **pA, nlspec *s, int oldt1)
{
    gretl_matrix *A = *pA;
    gretl_matrix *B = NULL;
    int T = s->t2 - s->t1 + 1;
    int At1 = gretl_matrix_get_t1(A);
    int At2 = gretl_matrix_get_t2(A);
    int m = A->cols;
    double x;
    int j, t, offt;

    B = gretl_matrix_alloc(T, m);
    if (B == NULL) {
	return E_ALLOC;
    }

    if (At1 == 0 && At2 == 0) {
	offt = s->t1 - oldt1;
    } else {
	offt = s->t1 - At1;
    }

    for (j=0; j<m; j++) {
	for (t=0; t<T; t++) {
	    x = gretl_matrix_get(A, t + offt, j);
	    gretl_matrix_set(B, t, j, x);
	}
    }

    gretl_matrix_set_t1(B, s->t1);
    gretl_matrix_set_t2(B, s->t2);

    gretl_matrix_replace(pA, B);

    return 0;
}

#define intmax3(a,b,c) (((a)>(b))? (((c)>(a))? (c) : (a)) : (((c)>(b))? (c) : (b)))
#define intmin3(a,b,c) (((a)<(b))? (((c)<(a))? (c) : (a)) : (((c)<(b))? (c) : (b)))
#define intdiff3(a,b,c) (a != b || a != c)

static int gmm_fix_datarows (nlspec *s)
{
    int et1 = gretl_matrix_get_t1(s->oc->e);
    int et2 = gretl_matrix_get_t2(s->oc->e);
    int Zt1 = gretl_matrix_get_t1(s->oc->Z);
    int Zt2 = gretl_matrix_get_t2(s->oc->Z);
    int oldt1 = s->t1;
    int err = 0;

    if ((et1 == 0 && et2 == 0) || (Zt1 == 0 && Zt2 == 0)) {
	/* we're jiggered: can't figure out alignment */
	return E_DATA;
    }

    s->t1 = intmax3(s->t1, et1, Zt1);
    s->t2 = intmin3(s->t2, et2, Zt2);
    s->nobs = s->t2 - s->t1 + 1;

    if (s->oc->e->rows > s->nobs) {
#if GMM_DEBUG
	fprintf(stderr, "gmm_fix_datarows: resizing e to %d rows\n", s->nobs);
#endif
	err = gmm_matrix_resize(&s->oc->e, s, oldt1);
    }

    if (s->oc->Z->rows > s->nobs && !err) {
#if GMM_DEBUG
	fprintf(stderr, "gmm_fix_datarows: resizing Z to %d rows\n", s->nobs);
#endif
	err = gmm_matrix_resize(&s->oc->Z, s, oldt1);
    }

    return err;
}

static int gmm_add_workspace (nlspec *s)
{
    int k = s->oc->noc;
    int err = 0;

    s->oc->tmp = gretl_matrix_alloc(s->nobs, k);
    s->oc->sum = gretl_column_vector_alloc(k);

    if (s->oc->tmp == NULL || s->oc->sum == NULL) {
	err = E_ALLOC;
    }

    return err;
}

/* Add the matrix of weights to the GMM specification: this must
   come after all the O.C.s are given, so that we can check the
   dimensions of the supplied matrix.
*/

int nlspec_add_weights (nlspec *s, const char *str)
{
    int k, err = 0;

    if (s->ci != GMM) {
	return E_TYPES;
    }

    if (s->oc == NULL) {
	gretl_errmsg_set(_("Weights must come after orthogonality conditions"));
	return E_DATA;
    }

    if (s->oc->W != NULL) {
	gretl_errmsg_set(_("Weights are already defined"));
	return E_DATA;
    }

#if GMM_DEBUG
    fprintf(stderr, "nlspec_add_weights: line = '%s'\n", str);
#endif

    if (gretl_scan_varname(str, s->oc->Wname) != 1) {
	return E_PARSE;
    }

    s->oc->W = get_matrix_by_name(s->oc->Wname);
    if (s->oc->W == NULL) {
	return gmm_unkvar(s->oc->Wname);
    }
    s->oc->userwts = 1;

    k = s->oc->noc;

    /* is the weight matrix of the correct dimensions? */
    if (s->oc->W->rows != k || s->oc->W->cols != k) {
	gretl_errmsg_sprintf(_("Weight matrix is of wrong size: should be "
			       "%d x %d"), k, k);
	err = E_DATA;
    }

#if GMM_DEBUG > 1
    fprintf(stderr, "weights matrix is %d x %d: %s\n",
	    s->oc->W->rows, s->oc->W->cols, err? "BAD" : "OK");
    if (!err) {
	gretl_matrix_print(s->oc->W, "Weights");
    }
#endif

    if (!err) {
	if (s->oc->e->rows != s->oc->Z->rows) {
	    err = gmm_fix_datarows(s);
	}
    }

    /* now we're ready to add the workspace matrices */

    if (!err) {
	err = gmm_add_workspace(s);
    }

    return err;
}

void maybe_add_gmm_residual (MODEL *pmod, const nlspec *spec,
			     const DATASET *dset)
{
    if (spec->oc != NULL && spec->oc->e != NULL && spec->oc->e->cols == 1) {
	int t, s = 0;

	if (pmod->uhat != NULL) {
	    free(pmod->uhat);
	}

	pmod->uhat = malloc(dset->n * sizeof *pmod->uhat);

	if (pmod->uhat != NULL) {
	    for (t=0; t<dset->n; t++) {
		if (t >= spec->t1 && t <= spec->t2) {
		    pmod->uhat[t] = spec->oc->e->val[s++];
		} else {
		    pmod->uhat[t] = NADBL;
		}
	    }
	    pmod->full_n = dset->n;
	}
    }
}

void nlspec_print_gmm_info (const nlspec *spec, PRN *prn)
{
    int i;

    if (spec->oc == NULL || spec->oc->lnames == NULL ||
	spec->oc->rnames == NULL) {
	return;
    }

    for (i=0; i<spec->oc->n_names; i++) {
	pprintf(prn, "orthog %s ; %s\n", spec->oc->lnames[i],
		spec->oc->rnames[i]);
    }

    pprintf(prn, "weights %s\n", spec->oc->Wname);
}

/* sanity check before proceeding with GMM estimation */

int check_gmm_requirements (nlspec *spec)
{
    int err = 0;

    if (spec->oc == NULL) {
	gretl_errmsg_set(_("No orthogonality conditions have been specified"));
	err = E_DATA;
    } else if (spec->oc->W == NULL) {
	/* add automatic weights matrix */
	spec->oc->W = gretl_identity_matrix_new(spec->oc->noc);
	if (spec->oc->W == NULL) {
	    err = E_ALLOC;
	} else {
	    strcpy(spec->oc->Wname, "auto");
	    err = gmm_add_workspace(spec);
	}
    } else {
	spec->oc->userwts = 1;
    }

    return err;
}

/* Update the column(s) of the LHS matrix in the O.C. set,
   after adusting the parameter values */

static int gmm_update_e (nlspec *s)
{
    gretl_matrix *e;
    double etj;
    int j, t, v;
    int err = 0;

    for (j=0; j<s->oc->e->cols && !err; j++) {
	v = s->oc->ecols[j].v;
	if (v > 0) {
#if GMM_DEBUG > 1
	    fprintf(stderr, "gmm_update_e: updating col %d from var %d\n",
		    j, v);
#endif
	    /* transcribe from series */
	    for (t=0; t<s->nobs; t++) {
		etj = s->dset->Z[v][t + s->t1];
		gretl_matrix_set(s->oc->e, t, j, etj);
	    }
	} else {
#if GMM_DEBUG > 1
	    fprintf(stderr, "gmm_update_e: updating col %d from matrix %s\n",
		    j, s->oc->ecols[j].mname);
#endif
	    /* transcribe from matrix */
	    e = get_matrix_by_name(s->oc->ecols[j].mname);
	    if (e == NULL) {
		err = 1;
	    } else if (e != s->oc->e) {
		for (t=0; t<s->nobs; t++) {
		    etj = gretl_matrix_get(e, t, s->oc->ecols[j].j);
		    gretl_matrix_set(s->oc->e, t, j, etj);
		}
	    }
	}
    }

    return err;
}

/* Carry out the required O.C. columnwise multiplications:
   e_i * Z_j */

static int gmm_multiply_ocs (nlspec *s)
{
    int err;

    err = gretl_matrix_columnwise_product(s->oc->e,
					  s->oc->Z,
					  s->oc->S,
					  s->oc->tmp);
    if (err) {
	fprintf(stderr, "gmm_multiply_ocs: err = %d from "
		"gretl_matrix_columnwise_product\n", err);
    }

#if GMM_DEBUG > 1
    gretl_matrix_print(s->oc->e, "gmm_multiply_ocs: s->oc->e");
    gretl_matrix_print(s->oc->Z, "gmm_multiply_ocs: s->oc->Z");
    gretl_matrix_print(s->oc->S, "gmm_multiply_ocs: s->oc->S");
    gretl_matrix_print(s->oc->tmp, "gmm_multiply_ocs: s->oc->tmp");
    fprintf(stderr, "err = %d\n", err);
#endif

    return err;
}

static double gmm_criterion (nlspec *s)
{
    double crit = 0.0;
    gretl_matrix *sum = s->oc->sum;
    int k = s->oc->noc;
    int i, j, err;

    err = gmm_multiply_ocs(s);
    if (err) {
	return NADBL;
    }

    /* compute column sums */
    for (j=0; j<k; j++) {
	sum->val[j] = 0.0;
	for (i=0; i<s->nobs; i++) {
	    sum->val[j] += gretl_matrix_get(s->oc->tmp, i, j);
	}
    }

    crit = gretl_scalar_qform(sum, s->oc->W, &err);
    if (!err) {
	crit = -crit;
    }

    return crit;
}

/* calculate the value of the GMM criterion given the current
   parameter values */

static double get_gmm_crit (const double *b, void *data)
{
    nlspec *s = (nlspec *) data;
    int err;

    update_coeff_values(b, s);

    err = nl_calculate_fvec(s);
    if (err) {
	fprintf(stderr, "get_gmm_crit: err on nl_calculate_fvec = %d\n", err);
	return NADBL;
    }

    err = gmm_update_e(s);
    if (err) {
	fprintf(stderr, "get_gmm_crit: err on gmm_update_e = %d\n", err);
	return NADBL;
    }

    s->crit = gmm_criterion(s);

#if GMM_DEBUG > 2
    gretl_matrix_print(s->oc->sum, "GMM: s->oc->sum");
    fprintf(stderr, "GMM: s->crit = %g\n", s->crit);
#endif

    return s->crit;
}

static int gmm_jacobian_calc (int m, int n, double *x, double *f,
			      int *iflag, void *p)
{
    nlspec *s = (nlspec *) p;
    double fac;
    int i, t, T;

    update_coeff_values(x, p);

    if (nl_calculate_fvec(s)) {
	*iflag = -1;
	return 1;
    }

    if (gmm_update_e(s)) {
	*iflag = -1;
	return 1;
    }

    if (gmm_multiply_ocs(s)) {
	*iflag = -1;
	return 1;
    }

    T = s->nobs;
    fac = sqrt((double) T) / T;

    for (i=0; i<m; i++) {
	f[i] = 0.0;
	for (t=0; t<T; t++) {
	    f[i] += gretl_matrix_get(s->oc->tmp, t, i);
	}
	f[i] *= fac;
    }

    return 0;
}

static int HAC_prewhiten (gretl_matrix *E, gretl_matrix *A)
{
    gretl_matrix_block *B;
    gretl_matrix *y, *X, *XTX;
    gretl_matrix *b, *e;
    int T = E->rows;
    int k = E->cols;
    double xtj, eti;
    int i, j, t;
    int err = 0;

    B = gretl_matrix_block_new(&y, T-1, 1,
			       &X, T-1, k,
			       &XTX, k, k,
			       &b, k, 1,
			       &e, k, 1,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    /* make matrix of RHS vars */
    for (j=0; j<k; j++) {
	for (t=1; t<T; t++) {
	    xtj = gretl_matrix_get(E, t-1, j);
	    gretl_matrix_set(X, t-1, j, xtj);
	}
    }

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
			      X, GRETL_MOD_NONE,
			      XTX, GRETL_MOD_NONE);

    err = gretl_matrix_cholesky_decomp(XTX);

    /* loop across LHS vars and compute coeffs */
    for (i=0; i<k && !err; i++) {
	for (t=1; t<T; t++) {
	    y->val[t-1] = gretl_matrix_get(E, t, i);
	}
	gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE,
				  y, GRETL_MOD_NONE,
				  b, GRETL_MOD_NONE);
	err = gretl_cholesky_solve(XTX, b);
	if (!err) {
	    for (j=0; j<k; j++) {
		gretl_matrix_set(A, i, j, b->val[j]);
	    }
	}
    }

#if GMM_DEBUG
    /* original VAR coefficients */
    gretl_matrix_print(A, "A");
#endif

    if (!err) {
	err = maybe_limit_VAR_coeffs(A, NULL, NULL, NULL);
	if (err) {
	    goto bailout;
	}

#if GMM_DEBUG
	/* (possibly) modified VAR coefficients */
	gretl_matrix_print(A, "A~");
#endif

	for (i=0; i<k; i++) {
	    /* get starting E values */
	    e->val[i] = gretl_matrix_get(E, 0, i);
	}

	/* Now "whiten" E using A~ */
	for (t=1; t<T; t++) {
	    /* re-use @b for fitted values */
	    gretl_matrix_multiply(A, e, b);
	    for (i=0; i<k; i++) {
		/* retrieve current value */
		eti = gretl_matrix_get(E, t, i);
		/* save lagged value for next round */
		e->val[i] = eti;
		/* substitute prediction error */
		gretl_matrix_set(E, t, i, eti - b->val[i]);
	    }
	}
    }

 bailout:

    gretl_matrix_block_destroy(B);

    return err;
}

static int gmm_HAC (gretl_matrix *E, gretl_matrix *V, hac_info *hinfo)
{
    static gretl_matrix *W;
    static gretl_matrix *Tmp;
    static gretl_matrix *A;
    static gretl_matrix *E2;
    int T, k;
    double w;
    int i, err = 0;

    if (E == NULL) {
	/* cleanup signal */
	gretl_matrix_free(W);
	gretl_matrix_free(Tmp);
	gretl_matrix_free(A);
	gretl_matrix_free(E2);
	W = Tmp = A = E2 = NULL;
	return 0;
    }

    T = E->rows;
    k = E->cols;

    if (W == NULL) {
	W = gretl_matrix_alloc(T, k);
	Tmp = gretl_matrix_alloc(k, k);
	if (W == NULL || Tmp == NULL) {
	    return E_ALLOC;
	}
	if (hinfo->whiten) {
	    A = gretl_matrix_alloc(k, k);
	    E2 = gretl_matrix_alloc(T, k);
	    if (A == NULL || E2 == NULL) {
		return E_ALLOC;
	    }
	}
    }

    if (hinfo->whiten) {
	gretl_matrix_copy_values(E2, E); /* keep a back-up */
	err = HAC_prewhiten(E, A);
	if (err) {
	    return err;
	}
    }

    if (data_based_hac_bandwidth()) {
	err = newey_west_bandwidth(E, NULL, hinfo->kern, hinfo->whiten,
				   &hinfo->h, &hinfo->bt);
	if (err) {
	    return err;
	}
    }

    gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
			      E, GRETL_MOD_NONE,
			      V, GRETL_MOD_NONE);

    for (i=1; i<=hinfo->h; i++) {
	if (hinfo->kern == KERNEL_QS) {
	    w = qs_hac_weight(hinfo->bt, i);
	} else {
	    w = hac_weight(hinfo->kern, hinfo->h, i);
	}
	gretl_matrix_inplace_lag(W, E, i);
	gretl_matrix_multiply_by_scalar(W, 2*w);
	gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
				  W, GRETL_MOD_NONE,
				  Tmp, GRETL_MOD_NONE);
	gretl_matrix_xtr_symmetric(Tmp);
	gretl_matrix_add_to(V, Tmp);
    }

    if (!gretl_matrix_is_symmetric(V)) {
	/* should we do this? */
	fprintf(stderr, "gmm_HAC: V is not symmetric\n");
	gretl_matrix_xtr_symmetric(V);
    }

    if (hinfo->whiten) {
	/* now we have to recolor */
	double aij;
	int j;

	/* make A into (I - A) */
	for (i=0; i<k; i++) {
	    for (j=0; j<k; j++) {
		aij = gretl_matrix_get(A, i, j);
		aij = (i == j)? 1 - aij: -aij;
		gretl_matrix_set(A, i, j, aij);
	    }
	}

	err = gretl_invert_general_matrix(A);

	if (!err) {
	    /* form V = (I-A)^{-1} * V * (I-A)^{-1}' */
	    gretl_matrix_copy_values(Tmp, V);
	    gretl_matrix_qform(A, GRETL_MOD_NONE,
			       Tmp, V, GRETL_MOD_NONE);
	}

	/* put back the original E values */
	gretl_matrix_copy_values(E, E2);
    }

    return err;
}

static void gmm_HAC_cleanup (void)
{
    gmm_HAC(NULL, NULL, NULL);
}

/* Calculate the covariance matrix for the GMM estimator.  We do an HC
   variant for non-time series, and a HAC variant for time-series,
   subject to (some) control via the libset.c apparatus.
*/

int gmm_add_vcv (MODEL *pmod, nlspec *s)
{
    gretl_matrix_block *B;
    gretl_matrix *V, *J, *S;
    gretl_matrix *m1, *m2, *m3;
    int i, j, k = s->ncoeff;
    int T = s->nobs;
    double *wa4;
    int m, n;
    int iflag = 0;
    double *f;
    int err = 0;

    m = gretl_matrix_cols(s->oc->tmp);
    n = s->ncoeff;

    wa4 = malloc(m * sizeof *wa4);
    if (wa4 == NULL) {
	return E_ALLOC;
    }

    B = gretl_matrix_block_new(&V, k, k,
			       &J, m, k,
			       &S, m, m,
			       &m1, k, m,
			       &m2, k, k,
			       &m3, k, k,
			       NULL);
    if (B == NULL) {
	free(wa4);
	return E_ALLOC;
    }

    if (using_HAC(s)) {
	err = gmm_HAC(s->oc->tmp, S, &s->oc->hinfo);
	gmm_HAC_cleanup();
    } else {
	err = gretl_matrix_multiply_mod(s->oc->tmp, GRETL_MOD_TRANSPOSE,
					s->oc->tmp, GRETL_MOD_NONE,
					S, GRETL_MOD_NONE);
    }

    if (!err) {
	double Tfac = sqrt((double) T) / T;

	gretl_matrix_divide_by_scalar(S, T);
	f = s->oc->sum->val;

	for (i=0; i<m; i++) {
	    f[i] = 0.0;
	    for (j=0; j<T; j++) {
		f[i] += gretl_matrix_get(s->oc->tmp, j, i);
	    }
	    f[i] *= Tfac;
	}

	fdjac2_(gmm_jacobian_calc, m, n, 0, s->coeff, f,
		J->val, m, &iflag, 0.0, wa4, s);

	if (iflag != 0) {
	    fprintf(stderr, "fdjac2_: iflag = %d\n", (int) iflag);
	    err = 1;
	}
    }

    if (!err) {
	err = gretl_matrix_multiply_mod(J, GRETL_MOD_TRANSPOSE,
					s->oc->W, GRETL_MOD_NONE,
					m1, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_qform(J, GRETL_MOD_TRANSPOSE,
				 s->oc->W, m2, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(m2);
    }

    if (!err) {
	err = gretl_matrix_qform(m1, GRETL_MOD_NONE,
				 S, m3, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_matrix_qform(m2, GRETL_MOD_NONE,
				 m3, V, GRETL_MOD_NONE);
    }

    if (!err) {
	err = gretl_model_write_vcv(pmod, V);
    }

    if (!err && using_HAC(s)) {
	if (s->oc->hinfo.kern == KERNEL_QS) {
	    gretl_model_set_hac_vcv_info(pmod, s->oc->hinfo.kern,
					 0, s->oc->hinfo.whiten,
					 s->oc->hinfo.bt);
	} else {
	    gretl_model_set_hac_vcv_info(pmod, s->oc->hinfo.kern,
					 s->oc->hinfo.h, s->oc->hinfo.whiten,
					 NADBL);
	}
    }

    if (!err) {
	/* set additional GMM info */
	double TGcrit = - s->crit / s->nobs;
	int l = s->oc->noc;

	pmod->ess = TGcrit / s->nobs; /* note the borrowing! */

	if (l > k && ((s->opt & OPT_V) || s->oc->step > 1)) {
	    gretl_model_set_int(pmod, "J_df", l - k);
	    gretl_model_set_double(pmod, "J_test", TGcrit);
	}

	if (s->oc->step > 1) {
	    gretl_model_set_int(pmod, "step", s->oc->step);
	}

	if (s->opt & OPT_T) {
	    pmod->opt |= OPT_T;
	} else if (s->opt & OPT_I) {
	    pmod->opt |= OPT_I;
	}
    }

    gretl_matrix_block_destroy(B);
    free(wa4);

    return err;
}

/* Compute summary criterion for change in parameter estimates.
   (Context: whether or not to continue GMM iteration, if the
   iterate option is selected.)
*/

static double gmm_get_icrit (nlspec *s, double *oldcoeff)
{
    double db, x = 0.0;
    int i;

    for (i=0; i<s->ncoeff; i++) {
	db = s->coeff[i] - oldcoeff[i];
	x += db * db;
	oldcoeff[i] = s->coeff[i];
    }

    return x;
}

/* In case we're proceeding to another round of GMM iteration,
   recompute the weights matrix */

static int gmm_recompute_weights (nlspec *s)
{
    gretl_matrix *W = s->oc->W;
    int err = 0;

    if (using_HAC(s)) {
	err = gmm_HAC(s->oc->tmp, W, &s->oc->hinfo);
    } else {
	err = gretl_matrix_multiply_mod(s->oc->tmp, GRETL_MOD_TRANSPOSE,
					s->oc->tmp, GRETL_MOD_NONE,
					W, GRETL_MOD_NONE);
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(W, s->nobs);
	err = gretl_invert_symmetric_matrix(W);
    }

    return err;
}

static void gmm_print_oc (nlspec *s, PRN *prn)
{
    gretl_matrix *V;
    int i, k = s->oc->noc;
    int T = s->nobs;
    int err = 0;

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	pprintf(prn, "gmm_print_oc: allocation failed!\n");
	return;
    }

    if (using_HAC(s)) {
	err = gmm_HAC(s->oc->tmp, V, &s->oc->hinfo);
    } else {
	err = gretl_matrix_multiply_mod(s->oc->tmp, GRETL_MOD_TRANSPOSE,
					s->oc->tmp, GRETL_MOD_NONE,
					V, GRETL_MOD_NONE);
    }

    pprintf(prn, "\n%s\n",
	    _("Orthogonality conditions - descriptive statistics"));
    pprintf(prn, "\n%10s  %10s %10s\n\n", _("OC"),
            _("mean"), _("std. dev"));

    for (i=0; i<k; i++) {
	pprintf(prn, "%10d    %10.6f", i, gretl_vector_get(s->oc->sum, i)/T);
	if (!err) {
	    pprintf(prn, " %10.6f\n", sqrt(gretl_matrix_get(V, i, i))/T);
	} else {
	    pprintf(prn, " %10s\n", "NA");
	}
    }
    pputc(prn, '\n');

    gretl_matrix_free(V);
}

#define NORM_DEBUG 0

static double gmm_log_10 (double x)
{
    double y;

    errno = 0;

    y = log(fabs(x));

    if (errno) {
	y = NADBL;
    } else {
	y /= log(10.0);
    }

    return y;
}

/* The user didn't specify a matrix of initial weights in
   GMM, and gretl assigned an identity matrix; here we try
   to see if we should scale the weights for the sake of
   keeping computer arithmetic in bounds.
*/

static int gmm_normalize_wts_1 (nlspec *s)
{
    double *coeff;
    double crit;

    coeff = copyvec(s->coeff, s->ncoeff);
    if (coeff == NULL) {
	return E_ALLOC;
    }

    crit = -1 * get_gmm_crit(coeff, s);

    if (crit > 0 && !na(crit)) {
	double m, lc = gmm_log_10(crit);

#if NORM_DEBUG
	fprintf(stderr, "maybe_preadjust_weights: crit=%g, lc=%g\n",
		crit, lc);
	gretl_matrix_print(s->oc->sum, "s->oc->sum");
#endif

	if (!na(lc) && (lc > 4 || lc < -4)) { /* was 5 */
	    if (lc > 0) {
		m = floor(lc/2);
	    } else {
		m = ceil(lc/3);
	    }
	    m = pow(10, -m);
	    fprintf(stderr, "GMM weights matrix: scaling by %g\n", m);
	    gretl_matrix_multiply_by_scalar(s->oc->W, m);
	    s->oc->scale_wt = m;
	}
    }

    free(coeff);

    return 0;
}

/* the idea here is to rescale rows and columns of the initial
   weights matrix so that 1-step becomes independent of the
   order of magnitude of the individual orthogonality conditions
*/

static int gmm_normalize_wts_2 (nlspec *s)
{
    double *coeff;
    double inicrit;
    int err = 0;

    /* don't mess with what might be a user-specified
       initialization of s->coeff */
    coeff = copyvec(s->coeff, s->ncoeff);
    if (coeff == NULL) {
	return E_ALLOC;
    }

    inicrit = -1 * get_gmm_crit(coeff, s);

    if (inicrit > 0 && !na(inicrit)) {
	int k = s->oc->noc;
	int n = s->oc->tmp->rows;
	gretl_vector *qvec;

	qvec = gretl_column_vector_alloc(k);

	if (qvec == NULL) {
	    err = E_ALLOC;
	} else {
	    double *q = qvec->val;
	    double x, sq = 0;
	    double scale = 0;
	    int i, j;
#if NORM_DEBUG
	    gretl_matrix_print(s->oc->W, "Weights (before)");
#endif
	    for (j=0; j<k; j++) {
		q[j] = 0;
		for (i=0; i<n; i++) {
		    x = gretl_matrix_get(s->oc->tmp, i, j);
		    q[j] += x*x;
		}
		sq += q[j];
		scale += gretl_matrix_get(s->oc->W, j, j);
	    }

	    scale /= n * sq;
	    fprintf(stderr, "gmm_normalize: scale = %g\n", scale);

	    for (i=0; i<k; i++) {
		for (j=0; j<k; j++) {
		    x = gretl_matrix_get(s->oc->W, i, j);
		    x /= scale * sqrt(q[i]*q[j]) * n;
		    gretl_matrix_set(s->oc->W, i, j, x);
		}
	    }

#if NORM_DEBUG
	    gretl_matrix_print(qvec, "q");
	    gretl_matrix_print(s->oc->W, "Weights (after)");
#endif
	    gretl_matrix_free(qvec);
	}
    }

    free(coeff);

    return err;
}

/* BFGS driver, for the case of GMM estimation.  We may have to handle
   both "inner" (BFGS) and "outer" iterations here: the "outer"
   iterations (if applicable) are re-runs of BFGS using an updated
   weights matrix.
*/

int gmm_calculate (nlspec *s, PRN *prn)
{
    static int normalize = 0;
    int full_fncount = 0;
    int full_grcount = 0;
    double itol = 1.0e-12, icrit = 1;
    double *oldcoeff = NULL;
    int maxit = 500;
    gretlopt iopt = s->opt;
    int outer_iters = 0;
    int outer_max = 1;
    int converged = 0;
    int err = 0;

    if (normalize == 0) {
	if (getenv("GMM_NORMALIZE_2") != NULL) {
	    normalize = 2;
	} else {
	    normalize = 1;
	}
    }

    if (s->opt & OPT_I) {
	/* iterate (but note: don't pass this option to BFGS,
	   since that would switch maximize to minimize)
	*/
	oldcoeff = copyvec(s->coeff, s->ncoeff);
	if (oldcoeff == NULL) {
	    err = E_ALLOC;
	} else {
	    outer_max = libset_get_int(GMM_MAXITER);
	}
	iopt ^= OPT_I;
    } else if (s->opt & OPT_T) {
	/* two-step */
	outer_max = 2;
    }

    /* experimental, 2013-12-30 */
    if (normalize == 1) {
	if (!s->oc->userwts) {
	    gmm_normalize_wts_1(s);
	}
    } else if (normalize == 2) {
	err = gmm_normalize_wts_2(s);
    }

    while (!err && outer_iters < outer_max && !converged) {

#if GMM_DEBUG
	fprintf(stderr, "GMM calling BFGS: outer_iters = %d\n",
		outer_iters);
#endif
	s->crit = 0.0;

	err = BFGS_max(s->coeff, s->ncoeff, maxit, s->tol,
		       &s->fncount, &s->grcount,
		       get_gmm_crit, C_GMM, NULL, s,
		       NULL, iopt, s->prn);

#if GMM_DEBUG
	fprintf(stderr, "GMM BFGS: err = %d\n", err);
#endif

	/* don't keep displaying certain things */
	iopt |= OPT_Q;

	full_fncount += s->fncount;
	full_grcount += s->grcount;

	if (!err && outer_max > 1) {
	    if (outer_max > 2) {
		icrit = gmm_get_icrit(s, oldcoeff);
		if (icrit < itol && outer_iters > 0) {
		    fprintf(stderr, "Breaking on icrit = %g at iteration %d\n",
			    icrit, outer_iters);
		    converged = 1;
		}
	    }
	    if (!converged && outer_iters < outer_max - 1) {
		err = gmm_recompute_weights(s);
	    }
	}

	if (err) {
	    fprintf(stderr, "Breaking on err = %d\n", err);
	} else if (!converged) {
	    if (++outer_iters == outer_max) {
		if (outer_max > 2) {
		    fprintf(stderr, "Breaking on max outer iter\n");
		    err = E_NOCONV;
		}
	    }
	}
    }

    if (!err) {
	s->oc->step = outer_iters;
#if GMM_DEBUG
	fprintf(stderr, "Total function evaluations = %d\n", full_fncount);
	fprintf(stderr, "Total gradient evaluations = %d\n", full_grcount);
#endif
    }

    if (!err && normalize == 1 && s->oc->scale_wt != 1) {
	/* we auto-scaled the identity matrix for weighting
	   purposes: undo the scaling here
	*/
	s->crit /= s->oc->scale_wt;
    }

    if (oldcoeff != NULL) {
	free(oldcoeff);
    }

    if (s->opt & OPT_V) {
	gmm_print_oc(s, prn);
    }

    gmm_HAC_cleanup();

    return err;
}

static int needs_rejigging (gretl_matrix *m, int t1, int t2)
{
    int T = t2 - t1 + 1;

    if (m->rows != T) {
	/* size is wrong */
	return 1;
    }

    if (gretl_matrix_is_dated(m)) {
	int mt1 = gretl_matrix_get_t1(m);
	int mt2 = gretl_matrix_get_t2(m);

	if (mt1 != t1 || mt2 != t2) {
	    /* offset is wrong */
	    return 1;
	}
    }

    return 0;
}

/* if we discover missing values after setting up the OC
   matrices, we'll have to shrink them down
*/

static int resize_oc_matrices (nlspec *s, int oldt1)
{
    int T = s->t2 - s->t1 + 1;
    int err = 0;

#if GMM_DEBUG
    fprintf(stderr, "resize_oc_matrices: t1=%d, t2=%d\n", s->t1, s->t2);
#endif

    if (needs_rejigging(s->oc->e, s->t1, s->t2)) {
#if GMM_DEBUG
	fprintf(stderr, "resize_oc_matrices: resizing e\n");
#endif
	err = gmm_matrix_resize(&s->oc->e, s, oldt1);
    }

    if (needs_rejigging(s->oc->Z, s->t1, s->t2) && !err) {
#if GMM_DEBUG
	fprintf(stderr, "resize_oc_matrices: resizing Z\n");
#endif
	err = gmm_matrix_resize(&s->oc->Z, s, oldt1);
    }

    if (!err) {
	gretl_matrix_reuse(s->oc->tmp, T, 0);
    }

    return 0;
}

static void gmm_set_HAC_info (nlspec *s)
{
    hac_info *hinfo = &s->oc->hinfo;

    if (dataset_is_time_series(s->dset) && !libset_get_bool(FORCE_HC)) {
	hinfo->whiten = libset_get_bool(PREWHITEN);
	hinfo->kern = libset_get_int(HAC_KERNEL);
	if (hinfo->kern == KERNEL_QS) {
	    hinfo->bt = libset_get_double(QS_BANDWIDTH);
	    hinfo->h = s->nobs - 1;
	} else {
	    hinfo->h = get_hac_lag(s->nobs);
	    hinfo->bt = 0.0;
	}
    } else {
	hinfo->kern = -1;
	hinfo->h = 0;
	hinfo->bt = 0.0;
	hinfo->whiten = 0;
    }
}

/* Check the e and Z matrices for missing values; while we're
   at it, if there's no error condition then set the HAC
   info for the GMM run.
*/

int gmm_missval_check_etc (nlspec *s)
{
    int orig_t1 = s->t1;
    int orig_t2 = s->t2;
    int t, t1, t2;
    int i, m, k, p, all_ok;
    double x;
    int err = 0;

#if GMM_DEBUG
    fprintf(stderr, "gmm_missval_check: initial s->t1=%d, s->t2=%d\n",
	    s->t1, s->t2);
#endif

#if 0
    /* thought: allow matrix version of GMM for vectors longer
       then the dataset -- but not ready yet */
    if (!gretl_matrix_is_dated(s->oc->e) &&
	!gretl_matrix_is_dated(s->oc->Z)) {
	int nr1 = s->oc->e->rows;
	int nr2 = s->oc->Z->rows;

	if (nr1 == nr2) {
	    s->t1 = 0;
	    s->t2 = nr1 - 1;
	}
    }
#endif

    if (gretl_matrix_is_dated(s->oc->e)) {
	int et1 = gretl_matrix_get_t1(s->oc->e);
	int et2 = gretl_matrix_get_t2(s->oc->e);

	if (et1 > s->t1) {
	    s->t1 = et1;
	}
	if (et2 < s->t2) {
	    s->t2 = et2;
	}
    }

    if (gretl_matrix_is_dated(s->oc->Z)) {
	int Zt1 = gretl_matrix_get_t1(s->oc->Z);
	int Zt2 = gretl_matrix_get_t2(s->oc->Z);

	if (Zt1 > s->t1) {
	    s->t1 = Zt1;
	}
	if (Zt2 < s->t2) {
	    s->t2 = Zt2;
	}
    }

#if GMM_DEBUG
    fprintf(stderr, " after checking matrix limits: t1=%d, t2=%d\n",
	    s->t1, s->t2);
#endif

    s->nobs = s->t2 - s->t1 + 1;

    err = nl_calculate_fvec(s);
    if (!err) {
	err = gmm_update_e(s);
    }

    if (err) {
	return err;
    }

    m = s->oc->e->cols;
    p = s->oc->Z->cols;

    for (t1=s->t1; t1<=s->t2; t1++) {
	k = t1 - s->t1;
	all_ok = 1;
	for (i=0; i<m; i++) {
	    x = gretl_matrix_get(s->oc->e, k, i);
	    if (na(x)) {
		all_ok = 0;
		break;
	    }
	}
	if (all_ok) {
	    for (i=0; i<p; i++) {
		x = gretl_matrix_get(s->oc->Z, k, i);
		if (na(x)) {
		    all_ok = 0;
		    break;
		}
	    }
	}
	if (all_ok) {
	    break;
	}
    }

    for (t2=s->t2; t2>=s->t1; t2--) {
	k = t2 - s->t1;
	all_ok = 1;
	for (i=0; i<m; i++) {
	    x = gretl_matrix_get(s->oc->e, k, i);
	    if (na(x)) {
		all_ok = 0;
		break;
	    }
	}
	if (all_ok) {
	    for (i=0; i<p; i++) {
		x = gretl_matrix_get(s->oc->Z, k, i);
		if (na(x)) {
		    all_ok = 0;
		    break;
		}
	    }
	}
	if (all_ok) {
	    break;
	}
    }

    if (t2 - t1 + 1 < s->ncoeff) {
	err = E_DF;
    }

    for (t=t1; t<=t2 && !err; t++) {
	k = t - s->t1;
	for (i=0; i<m && !err; i++) {
	    x = gretl_matrix_get(s->oc->e, k, i);
	    if (na(x)) {
		fprintf(stderr, "  after setting t1=%d, t2=%d, "
			"got NA for e(%d) at obs %d\n", t1, t2, i, t);
		err = E_MISSDATA;
	    }
	}
	if (!err) {
	    for (i=0; i<p && !err; i++) {
		x = gretl_matrix_get(s->oc->Z, k, i);
		if (na(x)) {
		    fprintf(stderr, "  after setting t1=%d, t2=%d, "
			    "got NA for Z(%d) at obs %d\n", t1, t2, i, t);
		    err = E_MISSDATA;
		}
	    }
	}
    }

    s->t1 = t1;
    s->t2 = t2;
    s->nobs = t2 - t1 + 1;

    if (!err && (s->t1 > orig_t1 || s->t2 < orig_t2)) {
	err = resize_oc_matrices(s, orig_t1);
    }

    if (!err) {
	gmm_set_HAC_info(s);
    }

#if GMM_DEBUG
    fprintf(stderr, "gmm_missval_check: on exit t1=%d, t2=%d\n",
	    s->t1, s->t2);
#endif

    return err;
}
