/*
 *  Copyright (c) 2007 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h" 
#include "usermat.h"
#include "nlspec.h"
#include "libset.h"

#include "../../minpack/minpack.h"

#define GMM_DEBUG 0

typedef struct colsrc_ colsrc;

struct colsrc_ {
    int v;                 /* ID number of variable n column */
    char mname[VNAMELEN];  /* or name of source matrix */
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

/* Determine the type of an element in an "orthog" specification */

static int oc_get_type (const char *name, const DATAINFO *pdinfo,
			int rhs, int *err)
{
    int *list = NULL;
    int v = varindex(pdinfo, name);

#if GMM_DEBUG
    fprintf(stderr, "oc_get_type: looking at '%s'\n", name);
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
    if (list == NULL) {
	*err = E_UNKVAR;
	return ARG_NONE;
    }
    
    if (rhs) {
	/* lists are OK on the right-hand side */
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
    } else {
	/* but not on the left */
	*err = E_TYPES;
    }

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
expand_selector_matrix (nlspec *s, int oldZcols, const char *mask)
{
    gretl_matrix *Z = s->oc->Z;
    gretl_matrix *S = s->oc->S;
    int Zcols = Z->cols;
    int j, err = 0;

#if GMM_DEBUG
    fprintf(stderr, "expand_selector_matrix\n");
    gretl_matrix_print(S, "S, on entry");
#endif    

    if (S == NULL) {
	/* starting from scratch: one O.C. will already be present */
	S = s->oc->S = gretl_matrix_alloc(2, Zcols);
	if (s->oc->S == NULL) {
	    err = E_ALLOC;
	} else {
	    for (j=0; j<Zcols; j++) {
		gretl_matrix_set(S, 0, j, (j < oldZcols)? 1 : 0);
		gretl_matrix_set(S, 1, j, (j < oldZcols)? 
				 (mask[j] ? 1 : 0) : 1);
	    }
	}
    } else {
	/* enlarging an existing S matrix */
	int i, k = s->oc->e->cols;
	double xij, *x;

	x = realloc(S->val, k * Zcols * sizeof *x);
	if (x == NULL) {
	    err = E_ALLOC;
	} else {
	    /* record original dimensions */
	    int sm = S->rows;
	    int sn = S->cols;

	    S->val = x;
	    S->rows = k;
	    S->cols = Zcols;

	    /* transcribe original entries */
	    for (j=sn-1; j>=0; j--) {
		for (i=sm-1; i>=0; i--) {
		    xij = x[Sidx(i,j,sm)];
		    gretl_matrix_set(S, i, j, xij);
		}
	    }

	    /* 0s in upper-right block */
	    for (i=0; i<sm; i++) {
		for (j=sn; j<Zcols; j++) {
		    gretl_matrix_set(S, i, j, 0);
		}
	    }

	    /* 1s in lower-right block */
	    for (j=sn; j<Zcols; j++) {
		gretl_matrix_set(S, k-1, j, 1);
	    }

	    /* mixed 0/1 in lower-left region */
	    for (j=0; j<sn; j++) {
		gretl_matrix_set(S, k-1, j, mask[j] ? 1 : 0);
	    }
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
add_new_cols_to_Z (nlspec *s, const gretl_matrix *M)
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
	err = expand_selector_matrix(s, oldZcols, mask + M->cols);
    }
    
    free(mask);

    return err;
}

/* Record the source (ID number of variable, or name of matrix) for a
   column on the left-hand side of the set of GMM orthogonality
   conditions.  We need this later for updating the values in the
   given column, after adjusting the parameters.
*/

static int 
push_column_source (nlspec *s, int v, const char *mname)
{
    colsrc *cols;
    int n = 0;

    if (s->oc->e != NULL) {
	n = s->oc->e->cols;
    }

    cols = realloc(s->oc->ecols, (n+1) * sizeof *cols);
    if (cols == NULL) {
	return E_ALLOC;
    }

    s->oc->ecols = cols;

    if (v <= 0) {
	cols[n].v = 0;
    } else {
	cols[n].v = v;
    }

    if (mname == NULL) {
	cols[n].mname[0] = '\0';
    } else {
	strcpy(cols[n].mname, mname);
    }

#if GMM_DEBUG
    fprintf(stderr, "push_column_source: added source %d: v=%d, mname='%s'\n",
	    n, cols[n].v, cols[n].mname);
#endif

    return 0;
}

static int oc_add_matrices (nlspec *s, int ltype, const char *lname,
			    int rtype, const char *rname,
			    const double **Z, const DATAINFO *pdinfo)
{
    gretl_matrix *e = NULL;
    gretl_matrix *M = NULL;
    int i, k, t, v;
    int free_e = 0;
    int free_M = 0;
    double x;
    int err = 0;

    /* First work on the left-hand side: have we been given a vector
       or a series?
    */

#if GMM_DEBUG
    fprintf(stderr, "oc_add_matrices: lname = '%s'\n", lname);
#endif

    if (ltype == ARG_MATRIX) {
	e = get_matrix_by_name(lname);
	err = push_column_source(s, 0, lname);
    } else if (ltype == ARG_SERIES) {
	v = varindex(pdinfo, lname);
	e = gretl_column_vector_alloc(s->nobs);
	if (e == NULL) {
	    err = E_ALLOC;
	} else {
	    for (t=0; t<s->nobs; t++) {
		x = Z[v][t + s->t1];
		gretl_vector_set(e, t, x);
	    }
	    free_e = s->oc->free_e = 1;
	    err = push_column_source(s, v, NULL);
	}
    } else {
	err = E_UNKVAR;
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
		    err = add_new_cols_to_Z(s, M);
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
	s->oc->noc += k;
    }

    return err;
}

/* Add a set of orthogonality conditions to a GMM specification.
   The source command takes the form:

   orthog { series | vector } ; { series | list | matrix }
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

    if (sscanf(str, "%15s ; %15s", lname, rname) != 2) {
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
	return E_UNKVAR;
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
    gretl_matrix *col;
    double etj;
    int j, m, t, v;
    int err = 0;

    m = s->oc->e->cols;

    for (j=0; j<m && !err; j++) {
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
	    fprintf(stderr, "gmm_update_e: updating col %d from vector %s\n",
		    j, s->oc->ecols[j].mname);
#endif
	    /* transcribe from vector */
	    col = get_matrix_by_name(s->oc->ecols[j].mname);
	    if (col == NULL) {
		err = 1;
	    } else if (col != s->oc->e) {
		for (t=0; t<s->nobs; t++) {
		    etj = col->val[t];
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
	s->crit = - s->crit;
    }

#if GMM_DEBUG > 2
    gretl_matrix_print(s->oc->e, "s->oc->e");
    gretl_matrix_print(s->oc->tmp, "s->oc->tmp");
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

static int newey_west (const gretl_matrix *E, int h,
		       gretl_matrix *V)
{
    gretl_matrix *W;
    double w;
    int i;

    W = gretl_matrix_alloc(E->rows, E->cols);
    if (W == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_zero(V);

    for (i=-h; i<=h; i++) {
	w = 1.0 - fabs((double) i) / (h + 1.0);
	gretl_matrix_inplace_lag(W, E, i);
	gretl_matrix_multiply_by_scalar(W, w);
	gretl_matrix_multiply_mod(E, GRETL_MOD_TRANSPOSE,
				  W, GRETL_MOD_NONE,
				  V, GRETL_MOD_CUMULATE);
    }

    gretl_matrix_free(W);

    return 0;
}

/* Calculate the covariance matrix for the GMM estimator.  Right now
   we do an HC variant for non-time series, and a HAC variant for
   time-series, subject to (some) control for the libset.c apparatus.
   Perhaps we should allow for selection of a "straight" asymptotic
   version?
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
    doublereal epsfcn = 0.0;
    integer m, n, ldjac;
    integer iflag = 0;
    double x, *f;
    int hac_lag = 0;
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

    if (dataset_is_time_series(s->dinfo) && !get_force_hc()) {
	hac_lag = get_hac_lag(s->nobs);
    }

    if (hac_lag == 0) {
	err = gretl_matrix_multiply_mod(s->oc->tmp, GRETL_MOD_TRANSPOSE,
					s->oc->tmp, GRETL_MOD_NONE,
					S, GRETL_MOD_NONE);
    } else {
	err = newey_west(s->oc->tmp, hac_lag, S);
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(S, s->nobs);

	f = s->oc->sum->val;

	for (i=0; i<m; i++) {
	    f[i] = 0.0;
	    for (j=0; j<T; j++) {
		f[i] += gretl_matrix_get(s->oc->tmp, j, i);
	    }
	    f[i] *= sqrt((double) T) / T;
	}

	fdjac2_(gmm_jacobian_calc, &m, &n, s->coeff, f, 
		J->val, &ldjac, &iflag, &epsfcn, wa4, s);

	if (iflag != 0) {
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
	if (hac_lag > 0) {
	    gretl_model_set_int(pmod, "hac_lag", hac_lag);
	}
    }

    if (!err) {
	/* set additional GMM info */
	int l = s->oc->noc;
	
	pmod->ess = - s->crit; /* note the borrowing! */
	if (l > k) {
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
    int hac_lag = 0;
    int err = 0;

    if (dataset_is_time_series(s->dinfo) && !get_force_hc()) {
	hac_lag = get_hac_lag(s->nobs);
    }

    if (hac_lag == 0) {
	err = gretl_matrix_multiply_mod(s->oc->tmp, GRETL_MOD_TRANSPOSE,
					s->oc->tmp, GRETL_MOD_NONE,
					W, GRETL_MOD_NONE);
    } else {
	err = newey_west(s->oc->tmp, hac_lag, W);
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
    int hac_lag = 0;
    int err = 0;

    if (dataset_is_time_series(s->dinfo) && !get_force_hc()) {
	hac_lag = get_hac_lag(s->nobs);
    }

    V = gretl_matrix_alloc(k, k);
    if (V == NULL) {
	pprintf(prn, "gmm_print_oc: allocation failed!\n");
	return;
    }

    if (hac_lag == 0) {
	err = gretl_matrix_multiply_mod(s->oc->tmp, GRETL_MOD_TRANSPOSE,
					s->oc->tmp, GRETL_MOD_NONE,
					V, GRETL_MOD_NONE);
    } else {
	err = newey_west(s->oc->tmp, hac_lag, V);
    }

    pprintf(prn, "  %s\n", _("Orthgonality"));
    pprintf(prn, "   %10s %10s %10s\n\n", _("condition"), 
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
	    outer_max = 50; /* arbitrary! */
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

	err = BFGS_max(s->coeff, s->ncoeff, maxit, s->tol, 
		       &s->fncount, &s->grcount, 
		       get_gmm_crit, C_GMM, NULL, s,
		       s->opt, s->prn);

	if (!err && outer_max > 1) {
	    if (outer_max > 2) {
		icrit = gmm_get_icrit(s, oldcoeff);
		if (icrit < itol && outer_iters > 0) {
		    fprintf(stderr, "Breaking on icrit = %g\n", icrit);
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
		}
	    }
	}
    }

    if (!err) {
	s->oc->step = outer_iters;
    }

    if (oldcoeff != NULL) {
	free(oldcoeff);
    }

    if (s->opt & OPT_V) {
	gmm_print_oc(s, prn);
    }

    return err;    
}

/* if we discover missing values after setting up the OC
   matrices, we'll have to shrink them down
*/

static int resize_oc_matrices (nlspec *s, int t1, int t2)
{
    gretl_matrix *e = NULL;
    gretl_matrix *Z = NULL;
    gretl_matrix *M = NULL;
    int T = t2 - t1 + 1;
    int m = s->oc->e->cols;
    int k = s->oc->Z->cols;
    double x;
    int j, t, offt;

    e = gretl_matrix_alloc(T, m);
    Z = gretl_matrix_alloc(T, k);
    M = gretl_matrix_alloc(T, s->oc->noc);

    if (e == NULL || Z == NULL || M == NULL) {
	gretl_matrix_free(e);
	gretl_matrix_free(Z);
	gretl_matrix_free(M);
	return E_ALLOC;
    }

    offt = t1 - s->t1;

    for (j=0; j<m; j++) {
	for (t=0; t<T; t++) {
	    x = gretl_matrix_get(s->oc->e, t + offt, j);
	    gretl_matrix_set(e, t, j, x);
	}
    }

    for (j=0; j<k; j++) {
	for (t=0; t<T; t++) {
	    x = gretl_matrix_get(s->oc->Z, t + offt, j);
	    gretl_matrix_set(Z, t, j, x);
	}
    }

    gretl_matrix_free(s->oc->e);
    gretl_matrix_free(s->oc->Z);
    gretl_matrix_free(s->oc->tmp);

    s->oc->e = e;
    s->oc->Z = Z;
    s->oc->tmp = M;

    return 0;
}

/* Check the e and Z matrices for missing values */

int gmm_missval_check (nlspec *s)
{
    int i, t, t1 = s->t1, t2 = s->t2;
    int m, k, p, all_ok;
    double x;
    int err = 0;

#if GMM_DEBUG
    fprintf(stderr, "gmm_missval_check: initial t1=%d, t2=%d\n",
	    s->t1, s->t2);
#endif

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

    if (!err && (t1 > s->t1 || t2 < s->t2)) {
	err = resize_oc_matrices(s, t1, t2);
    }

    if (!err) {
	s->t1 = t1;
	s->t2 = t2;
	s->nobs = t2 - t1 + 1;
    }

    return err;
}
