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

#include "../../minpack/minpack.h"

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
    gretl_matrix *e;     /* GMM residual, or LHS term in O.C. */
    gretl_matrix *Z;     /* instruments, or RHS terms in O.C. */
    gretl_matrix *W;     /* matrix of weights */
    gretl_matrix *tmp;   /* holds columnwise product of e and Z */
    gretl_matrix *sum;   /* holds column sums of tmp */
    gretl_matrix *S;     /* selector matrix for computing tmp */
    colsrc *ecols;       /* info on provenance of columns of 'e' */
    int noc;             /* total number of orthogonality conds. */
    int step;            /* number of estimation steps */
    int free_e;          /* indicator: should 'e' be freed? */
    int free_Z;          /* indicator: should 'Z' be freed? */
    hac_info hinfo;      /* HAC characteristics */
};

/* destructor for set of O.C. info */

void oc_set_destroy (ocset *oc)
{
    if (oc == NULL) {
	return;
    }

    if (oc->free_e) {
	gretl_matrix_free(oc->e);
    }

    if (oc->free_Z) {
	gretl_matrix_free(oc->Z);
    }

    gretl_matrix_free(oc->tmp);
    gretl_matrix_free(oc->sum);
    gretl_matrix_free(oc->S);

    free(oc->ecols);
	
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

	oc->free_e = 0;
	oc->free_Z = 0;
    }

    return oc;
}

static int gmm_unkvar (const char *s)
{
    sprintf(gretl_errmsg, _("Unknown variable '%s'"), s);
    return E_UNKVAR;
}

static int matrix_is_dated (const gretl_matrix *m)
{
    return !(m->t1 == 0 && m->t2 == 0);
}

/* Determine the type of an element in an "orthog" specification */

static int oc_get_type (const char *name, const DATAINFO *pdinfo,
			int rhs, int *err)
{
    int *list = NULL;
    int v = varindex(pdinfo, name);

#if GMM_DEBUG
    fprintf(stderr, "oc_get_type: looking at '%s' (v = %d, pdinfo->v = %d)\n", 
	    name, v, pdinfo->v);
#endif

    /* try for a series first */
    if (v >= 0 && v < pdinfo->v) {
	if (var_is_scalar(pdinfo, v)) {
	    *err = E_TYPES;
	    return ARG_NONE;
	} else {
	    return ARG_SERIES;
	}
    }

    /* then a matrix */
    if (get_matrix_by_name(name) != NULL) {
	return ARG_MATRIX;
    }

    /* failing that, a list */
    list = get_list_by_name(name);
    if (list != NULL) {
	int j;

	for (j=1; j<=list[0]; j++) {
	    if (list[j] < 0 || list[j] >= pdinfo->v) {
		*err = E_DATA;
		return ARG_NONE;
	    } else if (var_is_scalar(pdinfo, list[j])) {
		*err = E_TYPES;
		return ARG_NONE;
	    }
	}
	return ARG_LIST;
    }

    *err = gmm_unkvar(name);

    return ARG_NONE;
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

/* Record the source (ID number of variable, or name of matrix and
   column within matrix) for a column on the left-hand side of the set
   of GMM orthogonality conditions.  We need this later for updating
   the values in the given column, after adjusting the parameters.
*/

static int 
push_column_source (nlspec *s, int v, const char *mname)
{
    gretl_matrix *m;
    colsrc *cols;
    int j, k, n, p;

    m = (mname != NULL)? get_matrix_by_name(mname) : NULL;
    n = (s->oc->e != NULL)? s->oc->e->cols : 0;
    k = (m != NULL)? m->cols : 1;

    cols = realloc(s->oc->ecols, (n + k) * sizeof *cols);
    if (cols == NULL) {
	return E_ALLOC;
    }

    s->oc->ecols = cols;

    for (j=0; j<k; j++) {
	p = n + j;
	cols[p].v = (v < 0)? 0 : v;
	cols[p].j = j;

	if (mname == NULL) {
	    cols[p].mname[0] = '\0';
	} else {
	    strcpy(cols[p].mname, mname);
	}

#if GMM_DEBUG
	fprintf(stderr, "push_column_source: added source %d: v=%d, "
		"mname='%s', j=%d\n", p, cols[p].v, cols[p].mname, 
		cols[p].j);
#endif
    }

    return 0;
}

static int oc_add_matrices (nlspec *s, int ltype, const char *lname,
			    int rtype, const char *rname,
			    const double **Z, const DATAINFO *pdinfo)
{
    gretl_matrix *e = NULL;
    gretl_matrix *M = NULL;
    int oldecols = 0;
    int i, k, t, v;
    int free_e = 0;
    int free_M = 0;
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

    if (ltype == ARG_MATRIX) {
	e = get_matrix_by_name(lname);
	err = push_column_source(s, 0, lname);
    } else if (ltype == ARG_LIST) {
	int *list = get_list_by_name(lname);

	k = list[0];
	e = gretl_matrix_alloc(s->nobs, k);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    e->t1 = s->t1;
	    e->t2 = s->t2;
	    for (i=0; i<k && !err; i++) {
		v = list[i + 1];
		for (t=0; t<s->nobs; t++) {
		    x = Z[v][t + s->t1];
		    gretl_matrix_set(e, t, i, x);
		}
		err = push_column_source(s, v, NULL);
	    }
	    free_e = s->oc->free_e = 1;
	}	
    } else if (ltype == ARG_SERIES) {
	v = varindex(pdinfo, lname);
	e = gretl_column_vector_alloc(s->nobs);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    e->t1 = s->t1;
	    e->t2 = s->t2;
	    for (t=0; t<s->nobs; t++) {
		x = Z[v][t + s->t1];
		gretl_vector_set(e, t, x);
	    }
	    free_e = s->oc->free_e = 1;
	    err = push_column_source(s, v, NULL);
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

    if (rtype == ARG_MATRIX) {
	M = get_matrix_by_name(rname);
	k = gretl_matrix_cols(M);
    } else if (rtype == ARG_LIST) {
	int *list = get_list_by_name(rname);

	k = list[0];
	M = gretl_matrix_alloc(s->nobs, k);
	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    M->t1 = s->t1;
	    M->t2 = s->t2;
	    for (i=0; i<k && !err; i++) {
		v = list[i + 1];
		for (t=0; t<s->nobs; t++) {
		    x = Z[v][t + s->t1];
		    gretl_matrix_set(M, t, i, x);
		}
	    }
	    free_M = s->oc->free_Z = 1;
	}
    } else {
	/* name of single series */
	v = varindex(pdinfo, rname);
	k = 1;
	M = gretl_column_vector_alloc(s->nobs);
	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    M->t1 = s->t1;
	    M->t2 = s->t2;
	    for (t=0; t<s->nobs; t++) {
		x = Z[v][t + s->t1];
		gretl_vector_set(M, t, x);
	    }
	    free_M = s->oc->free_Z = 1;
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
	    if (!s->oc->free_e) {
		/* current 'e' is pointer to external matrix: this
		   must now be changed */
		s->oc->e = gretl_matrix_copy(s->oc->e);
		if (s->oc->e == NULL) {
		    err = E_ALLOC;
		}
	    }

	    if (!err) {
		err = gretl_matrix_inplace_colcat(s->oc->e, e, NULL);
#if GMM_DEBUG > 1
		fprintf(stderr, "oc_add_matrices: expanded e\n");
		gretl_matrix_print(s->oc->e, "e");
#endif
	    }

	    if (free_e) {
		/* don't leak temporary matrix */
		gretl_matrix_free(e);
	    }

	    if (!err) {
		if (!s->oc->free_Z) {
		    /* current 'Z' is pointer to external matrix */
		    s->oc->Z = gretl_matrix_copy(s->oc->Z);
		    if (s->oc->Z == NULL) {
			err = E_ALLOC;
		    }
		}
		if (!err) {
		    err = add_new_cols_to_Z(s, M, oldecols);
#if GMM_DEBUG > 1
		    fprintf(stderr, "oc_add_matrices: expanded Z\n");
		    gretl_matrix_print(s->oc->Z, "Z");
#endif
		}
	    }

	    if (free_M) {
		gretl_matrix_free(M);
	    }

	    s->oc->free_e = s->oc->free_Z = 1;
	}
    }

    if (!err) {
	s->oc->noc += k * (s->oc->e->cols - oldecols); /* FIXME? */
    }

    return err;
}

/* Add a set of orthogonality conditions to a GMM specification.
   The source command takes the form:

   orthog { series | list | matrix } ; { series | list | matrix }
*/

int 
nlspec_add_orthcond (nlspec *s, const char *str,
		     const double **Z, const DATAINFO *pdinfo)
{
    char lname[VNAMELEN];
    char rname[VNAMELEN];
    int ltype = ARG_NONE, rtype = ARG_NONE;
    int err = 0;

    if (s->ci != GMM) {
	return E_TYPES;
    }

    if (s->oc != NULL && s->oc->W != NULL) {
	strcpy(gretl_errmsg, _("Orthogonality conditions must come before the "
	       "weights matrix"));
	return E_DATA;
    }

#if GMM_DEBUG
    fprintf(stderr, "nlspec_add_orthcond: line = '%s'\n", str);
#endif

    str += strspn(str, " ");

    if (sscanf(str, "%15[^; ] ; %15s", lname, rname) != 2) {
	return E_PARSE;
    }

    if (s->oc == NULL) {
	s->oc = oc_set_new();
	if (s->oc == NULL) {
	    return E_ALLOC;
	}
    }

    ltype = oc_get_type(lname, pdinfo, 0, &err);

    if (!err) {
	rtype = oc_get_type(rname, pdinfo, 1, &err);
    }

    if (!err) {
	err = oc_add_matrices(s, ltype, lname, rtype, rname, 
			      Z, pdinfo);
    }

    if (err) {
	nlspec_destroy_arrays(s);
	oc_set_destroy(s->oc);
	s->oc = NULL;
    } 

    return err;
}

static int matrix_t1 (const gretl_matrix *m, int t1)
{
    if (m->t1 == 0 && m->t2 == 0) {
	return t1;
    } else {
	return m->t1;
    }
}

static int gmm_matrix_resize (gretl_matrix **pA, nlspec *s, int oldt1)
{
    gretl_matrix *A = *pA;
    gretl_matrix *B = NULL;
    int T = s->t2 - s->t1 + 1;
    int m = A->cols;
    double x;
    int j, t, offt;

    B = gretl_matrix_alloc(T, m);
    if (B == NULL) {
	return E_ALLOC;
    }

    offt = s->t1 - matrix_t1(A, oldt1); /* ?? */

    for (j=0; j<m; j++) {
	for (t=0; t<T; t++) {
	    x = gretl_matrix_get(A, t + offt, j);
	    gretl_matrix_set(B, t, j, x);
	}
    }

    B->t1 = s->t1;
    B->t2 = s->t2;

    gretl_matrix_free(A);
    *pA = B;

    return 0;
}


#define intmax3(a,b,c) (((a)>(b))? (((c)>(a))? (c) : (a)) : (((c)>(b))? (c) : (b)))
#define intmin3(a,b,c) (((a)<(b))? (((c)<(a))? (c) : (a)) : (((c)<(b))? (c) : (b)))
#define intdiff3(a,b,c) (a != b || a != c)

static int gmm_fix_datarows (nlspec *s)
{
    int et1 = s->oc->e->t1;
    int et2 = s->oc->e->t2;
    int Zt1 = s->oc->Z->t1;
    int Zt2 = s->oc->Z->t2;
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

/* Add the matrix of weights to the GMM specification: this must
   come after all the O.C.s are given, so that we can check the
   dimensions of the given matrix.
*/

int nlspec_add_weights (nlspec *s, const char *str)
{
    char name[VNAMELEN];
    int k, err = 0;

    if (s->ci != GMM) {
	return E_TYPES;
    }

    if (s->oc == NULL) {
	strcpy(gretl_errmsg, _("Weights must come after orthogonality conditions"));
	return E_DATA;
    }

    if (s->oc->W != NULL) {
	strcpy(gretl_errmsg, _("Weights are already defined"));
	return E_DATA;
    }

#if GMM_DEBUG
    fprintf(stderr, "nlspec_add_weights: line = '%s'\n", str);
#endif

    if (sscanf(str, "%15s", name) != 1) {
	return E_PARSE;
    } 

    s->oc->W = get_matrix_by_name(name);
    if (s->oc->W == NULL) {
	return gmm_unkvar(name);
    }

    k = s->oc->noc;

    /* is the weight matrix of the correct dimensions? */
    if (s->oc->W->rows != k || s->oc->W->cols != k) {
	sprintf(gretl_errmsg, _("Weight matrix is of wrong size: should be "
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
	s->oc->tmp = gretl_matrix_alloc(s->nobs, k);
	s->oc->sum = gretl_column_vector_alloc(k);
	if (s->oc->tmp == NULL || s->oc->sum == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

/* sanity check before proceeding with GMM estimation */

int check_gmm_requirements (nlspec *spec)
{
    if (spec->oc == NULL) {
	strcpy(gretl_errmsg, _("No orthogonality conditions have been specified"));
	return 1;
    }

    if (spec->oc->W == NULL) {
	strcpy(gretl_errmsg, _("No weights have been specified"));
	return 1;
    }

    return 0;
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
#if GMM_DEBUG
	    fprintf(stderr, "gmm_update_e: updating col %d from var %d\n",
		    j, v);
#endif
	    /* transcribe from series */
	    for (t=0; t<s->nobs; t++) {
		etj = (*s->Z)[v][t + s->t1];
		gretl_matrix_set(s->oc->e, t, j, etj);
	    }
	} else {
#if GMM_DEBUG
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
    int err = 0;

    if (s->oc->S == NULL) {
	/* simple case: no selection mechanism needed */
	err = gretl_matrix_columnwise_product(s->oc->e,
					      s->oc->Z,
					      s->oc->tmp);
    } else {
	/* select the wanted i, j pairs */
	int i, j, k, m, t, p = 0;
	double x, y;

	m = s->oc->e->cols;
	k = s->oc->Z->cols;

	for (i=0; i<m; i++) {
	    for (j=0; j<k; j++) {
		if (gretl_matrix_get(s->oc->S, i, j) != 0) {
#if GMM_DEBUG
		    fprintf(stderr, "O.C. %d: multiplying col %d of e "
			    "into col %d of Z\n", p, i, j);
#endif
		    for (t=0; t<s->nobs; t++) {
			x = gretl_matrix_get(s->oc->e, t, i);
			y = gretl_matrix_get(s->oc->Z, t, j);
			gretl_matrix_set(s->oc->tmp, t, p, x * y);
		    }
		    p++;
		}
	    }
	}
    }

#if GMM_DEBUG > 1
    gretl_matrix_print(s->oc->e, "gmm_multiply_ocs: s->oc->e");
    gretl_matrix_print(s->oc->Z, "gmm_multiply_ocs: s->oc->Z");
    gretl_matrix_print(s->oc->tmp, "gmm_multiply_ocs: s->oc->tmp");
    fprintf(stderr, "err = %d\n", err);
#endif

    return err;
}

/* calculate the value of the GMM criterion given the current
   parameter values */

static double get_gmm_crit (const double *b, void *p)
{
    nlspec *s = (nlspec *) p;
    int i, k, t;
    int err = 0;

    update_coeff_values(b, s);

    if (nl_calculate_fvec(s)) {
	return NADBL;
    }

    if (gmm_update_e(s)) {
	return NADBL;
    }

    if (gmm_multiply_ocs(s)) {
	return NADBL;
    }
    
    k = s->oc->noc;

    for (i=0; i<k; i++) {
	s->oc->sum->val[i] = 0.0;
	for (t=0; t<s->nobs; t++) {
	    s->oc->sum->val[i] += gretl_matrix_get(s->oc->tmp, t, i);
	}
    }

    s->crit = gretl_scalar_qform(s->oc->sum, s->oc->W, &err);
    if (!err) {
	s->crit = -s->crit;
    }

#if GMM_DEBUG > 2
    gretl_matrix_print(s->oc->sum, "s->oc->sum");
    fprintf(stderr, "GMM: s->crit = %g\n", s->crit);
#endif

    if (err) {
	s->crit = NADBL;
    }
	
    return s->crit;
}

/* GMM callback for Minpack's fdjac2 function */

static int 
gmm_jacobian_calc (integer *m, integer *n, double *x, double *f,
		   integer *iflag, void *p)
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

    for (i=0; i<*m; i++) {
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
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *XTX = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *e = NULL;
    int T = E->rows;
    int k = E->cols;
    double xtj, eti;
    int i, j, t;
    int err = 0;

    y = gretl_column_vector_alloc(T - 1);
    X = gretl_matrix_alloc(T - 1, k);
    XTX = gretl_matrix_alloc(k, k);
    b = gretl_column_vector_alloc(k);
    e = gretl_column_vector_alloc(k);

    if (y == NULL || X == NULL || XTX == NULL ||
	b == NULL || e == NULL) {
	err = E_ALLOC;
	goto bailout;
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
	gretl_matrix *U = NULL;
	gretl_matrix *S = NULL;
	gretl_matrix *V = NULL;
	gretl_matrix *tmp = NULL;
	int amod = 0;

	err = gretl_matrix_SVD(A, &U, &S, &V);
	if (err) {
	    goto bailout;
	}

	for (i=0; i<k; i++) {
	    if (S->val[i] > 0.97) {
		S->val[i] = 0.97;
		amod = 1;
	    }
	}

	if (amod) {
	    tmp = gretl_matrix_dot_op(U, S, '*', &err);
	    gretl_matrix_multiply(tmp, V, A);
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
	    /* re-use 'b' for fitted values */
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

	gretl_matrix_free(U);
	gretl_matrix_free(S);
	gretl_matrix_free(V);
	gretl_matrix_free(tmp);
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(XTX);
    gretl_matrix_free(b);
    gretl_matrix_free(e);

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
	err = newey_west_bandwidth(E, hinfo->kern, &hinfo->h,
				   &hinfo->bt);
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
	/* now we "re-color" */
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

	/* form V = (I-A)^{-1} * V * (I-A)^{-1}' */
	if (!err) {
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
    gretl_matrix *V = NULL;
    gretl_matrix *J = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *m1 = NULL;
    gretl_matrix *m2 = NULL;
    gretl_matrix *m3 = NULL;
    int i, j, k = s->ncoeff;
    int T = s->nobs;

    doublereal *wa4 = NULL;
    doublereal epsfcn = 0.0; /* could be > 0.0 */
    integer m, n, ldjac;
    integer iflag = 0;

    double x, *f;
    int err = 0;

    m = gretl_matrix_cols(s->oc->tmp);
    n = s->ncoeff;
    ldjac = m; /* leading dimension of jac array */

    wa4 = malloc(m * sizeof *wa4);
    V = gretl_matrix_alloc(k, k);
    J = gretl_matrix_alloc(m, k);
    S = gretl_matrix_alloc(m, m);
    m1 = gretl_matrix_alloc(k, m);
    m2 = gretl_matrix_alloc(k, k);
    m3 = gretl_matrix_alloc(k, k);

    if (V == NULL || J == NULL || S == NULL ||
	m1 == NULL || m2 == NULL || m3 == NULL ||
	wa4 == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (s->oc->hinfo.kern != 0) {
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

	fdjac2_(gmm_jacobian_calc, &m, &n, s->coeff, f, 
		J->val, &ldjac, &iflag, &epsfcn, wa4, s);

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
	for (i=0; i<k; i++) {
	    for (j=0; j<=i; j++) {
		x = gretl_matrix_get(V, i, j);
		pmod->vcv[ijton(i, j, k)] = x;
		if (j == i) {
		   pmod->sderr[i] = sqrt(x);
		} 
	    }
	}
	if (s->oc->hinfo.kern) {
	    gretl_model_set_int(pmod, "hac_kernel", s->oc->hinfo.kern);
	    if (s->oc->hinfo.kern == KERNEL_QS) {
		gretl_model_set_double(pmod, "qs_bandwidth", s->oc->hinfo.bt);
	    } else {
		gretl_model_set_int(pmod, "hac_lag", s->oc->hinfo.h);
	    }
	    if (s->oc->hinfo.whiten) {
		gretl_model_set_int(pmod, "hac_prewhiten", 1);
	    }
	}
    }

    if (!err) {
	/* set additional GMM info */
	int l = s->oc->noc;
	
	pmod->ess = - s->crit; /* note the borrowing! */

	if (l > k && ((s->opt & OPT_V) || s->oc->step > 1)) {
	    gretl_model_set_int(pmod, "J_df", l - k);
	    gretl_model_set_double(pmod, "J_test", pmod->ess / s->nobs);
	}

	if (s->oc->step > 1) {
	    gretl_model_set_int(pmod, "step", s->oc->step);
	}
    }

 bailout:

    gretl_matrix_free(V);
    gretl_matrix_free(J);
    gretl_matrix_free(S);
    gretl_matrix_free(m1);
    gretl_matrix_free(m2);
    gretl_matrix_free(m3);
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

    if (s->oc->hinfo.kern) {
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

    if (s->oc->hinfo.kern) {
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

/* Driver for BFGS, case of GMM estimation.  We may have to handle
   both "inner" (BFGS) and "outer" iterations here: the "outer"
   iterations (if applicable) are re-runs of BFGS using an updated
   weights matrix.
*/

int gmm_calculate (nlspec *s, double *fvec, double *jac, PRN *prn)
{
    int full_fncount = 0;
    int full_grcount = 0;
    double itol = 1.0e-12, icrit = 1;
    double *oldcoeff = NULL;
    int maxit = get_bfgs_maxiter();
    int outer_iters = 0;
    int outer_max = 1;
    int converged = 0;
    int err = 0;

    if (s->opt & OPT_I) {
	/* iterate */
	oldcoeff = copyvec(s->coeff, s->ncoeff);
	if (oldcoeff == NULL) {
	    err = E_ALLOC;
	} else {
	    outer_max = 200; /* arbitrary! */
	}
    } else if (s->opt & OPT_T) {
	/* two-step */
	outer_max = 2;
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
		       s->opt, s->prn);

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
	    outer_iters++;
	    if (outer_iters == outer_max) {
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

    if (matrix_is_dated(m)) {
	if (m->t1 != t1 || m->t2 != t2) {
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

    if (dataset_is_time_series(s->dinfo) && !get_force_hc()) {
	hinfo->whiten = get_hac_prewhiten();
	hinfo->kern = get_hac_kernel();
	if (hinfo->kern == KERNEL_QS) {
	    hinfo->bt = get_qs_bandwidth();
	    hinfo->h = s->nobs - 1;
	} else {
	    hinfo->h = get_hac_lag(s->nobs);
	    hinfo->bt = 0.0;
	}
    } else {
	hinfo->kern = 0;
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

    if (matrix_is_dated(s->oc->e)) {
	if (s->oc->e->t1 > s->t1) {
	    s->t1 = s->oc->e->t1;
	}
	if (s->oc->e->t2 < s->t2) {
	    s->t2 = s->oc->e->t2;
	}
    }

    if (matrix_is_dated(s->oc->Z)) {
	if (s->oc->Z->t1 > s->t1) {
	    s->t1 = s->oc->Z->t1;
	}
	if (s->oc->Z->t2 < s->t2) {
	    s->t2 = s->oc->Z->t2;
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
