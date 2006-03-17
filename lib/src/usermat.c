/*
 *  Copyright (c) by Allin Cottrell
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
#include "gretl_matrix.h"
#include "gretl_func.h"
#include "libset.h"
#include "usermat.h"

#include <errno.h>

#define MDEBUG 0

#define LEVEL_AUTO -1

enum {
    TRANSPOSE_NOT_OK,
    TRANSPOSE_OK
};

typedef struct user_matrix_ user_matrix;

struct user_matrix_ {
    gretl_matrix *M;
    int level;
    char name[VNAMELEN];
};

static user_matrix **matrices;
static int n_matrices;

static int 
make_slices (const char *s, int m, int n, int **rslice, int **cslice);
static int delete_matrix_by_name (const char *name);

static double ***gZ;
static DATAINFO *gdinfo;

static void 
usermat_publish_dataset (double ***pZ, DATAINFO *pdinfo)
{
    gZ = pZ;
    gdinfo = pdinfo;
}

static void usermat_unpublish_dataset (double ***pZ)
{
    *pZ = *gZ;
    gZ = NULL;
    gdinfo = NULL;
}

static user_matrix *user_matrix_new (gretl_matrix *M, const char *name)
{
    user_matrix *u;

    u = malloc(sizeof *u);
    if (u == NULL) {
	return NULL;
    }

    u->M = M;
    u->level = gretl_function_stack_depth();
    *u->name = '\0';
    strncat(u->name, name, VNAMELEN - 1);

    return u;
}

static int add_user_matrix (gretl_matrix *M, const char *name)
{
    user_matrix **tmp;

    if (M == NULL) {
	return 0;
    }

    if (check_varname(name)) {
	return E_DATA;
    }

    tmp = realloc(matrices, (n_matrices + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    } else {
	matrices = tmp;
    }

    if (is_user_matrix(M)) {
	/* ensure uniqueness of matrix pointers */
	gretl_matrix *Mcpy = gretl_matrix_copy(M);

	if (Mcpy == NULL) {
	    return E_ALLOC;
	}
	M = Mcpy;
    }

    matrices[n_matrices] = user_matrix_new(M, name);

    if (matrices[n_matrices] == NULL) {
	return E_ALLOC;
    }

    n_matrices++;

#if MDEBUG
    fprintf(stderr, "add_user_matrix: allocated '%s' at %p (M at %p, %dx%d)\n",
	    name, (void *) matrices[n_matrices-1], (void *) M,
	    gretl_matrix_rows(M), gretl_matrix_cols(M));
#endif

    return 0;
}

static int 
matrix_insert_diagonal (gretl_matrix *M, const gretl_matrix *S,
			int mr, int mc)
{
    double x;
    int i, len = gretl_vector_get_length(S);

    if (mr != mc || len != mr) {
	return E_NONCONF;
    }

    for (i=0; i<len; i++) {
	x = gretl_vector_get(S, i);
	gretl_matrix_set(M, i, i, x);
    }
    
    return 0;
}

static int matrix_insert_submatrix (gretl_matrix *M, const gretl_matrix *S,
				    const char *mask)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int sr = gretl_matrix_rows(S);
    int sc = gretl_matrix_cols(S);
    int *rslice = NULL;
    int *cslice = NULL;
    int err = 0;

    if (sr > mr || sc > mc) {
	err = E_NONCONF;
    }

    if (!err && !strcmp(mask, "[diag]")) {
	return matrix_insert_diagonal(M, S, mr, mc);
    }

    if (!err) {
	err = make_slices(mask, mr, mc, &rslice, &cslice);
#if MDEBUG
	printlist(rslice, "rslice");
	printlist(cslice, "cslice");
	fprintf(stderr, "M = %d x %d, S = %d x %d\n", mr, mc, sr, sc);
#endif
    }

    if (!err) {
	if (rslice != NULL && rslice[0] != sr) {
	    err = E_NONCONF;
	} else if (cslice != NULL && cslice[0] != sc) {
	    err = E_NONCONF;
	}
    }

    if (!err) {
	int i, j, k, l;
	int mi, mj;
	double x;

	k = 0;
	for (i=0; i<sr; i++) {
	    mi = (rslice == NULL)? k++ : rslice[i+1] - 1;
	    l = 0;
	    for (j=0; j<sc; j++) {
		mj = (cslice == NULL)? l++ : cslice[j+1] - 1;
		x = gretl_matrix_get(S, i, j);
		gretl_matrix_set(M, mi, mj, x);
	    }
	}
    }

    free(rslice);
    free(cslice);

    return err;
}

static int replace_user_matrix (user_matrix *u, gretl_matrix *M,
				gretl_matrix **R, const char *mask)
{
    int err = 0;

    if (R != NULL) {
#if MDEBUG
	fprintf(stderr, "replace_user_matrix: u->M = %p (matrix '%s')\n",
		u->M, u->name);
#endif
	*R = u->M;
    }  

    if (M == NULL) {
	return 0;
    }

    if (mask != NULL && *mask != '\0') {
	/* the new matrix M is actally a submatrix, to be written
	   into the original matrix, u->M */
	err = matrix_insert_submatrix(u->M, M, mask);
	if (!is_user_matrix(M)) {
	    /* is this always right? */
	    gretl_matrix_free(M);
	}
    } else if (M != u->M) {
#if MDEBUG
	fprintf(stderr, " freeing u->M at %p, replacing with matrix at %p\n", u->M, M);
#endif
	gretl_matrix_free(u->M);
	u->M = M;
    } else {
	fprintf(stderr, " no-op: 'replacing' matrix with itself!\n");
    }

    return err;
}

/**
 * is_user_matrix:
 * @m: gretl_matrix to test.
 *
 * Returns: 1 if the matrix @m is saved on the stack of matrices,
 * else 0.
 */

int is_user_matrix (const gretl_matrix *m)
{
    int i;

    for (i=0; i<n_matrices; i++) {
	if (m == matrices[i]->M) {
	    return 1;
	}
    }

    return 0;
}

static int name_is_variable (const char *name, const DATAINFO *pdinfo)
{
    int v = varindex(pdinfo, name);
    int ret = 0;

    if (v < pdinfo->v) {
	ret = 1;
    }

    return ret;
}

/* If slevel is LEVEL_AUTO search at the current function execution
   depth, otherwise search at the function execution depth given by
   slevel.
*/

static gretl_matrix *
real_get_matrix_by_name (const char *name, int slevel, const DATAINFO *pdinfo)
{
    int level, i;

    if (slevel == LEVEL_AUTO) {
	level = gretl_function_stack_depth();
    } else {
	level = slevel;
    }

    /* if a series or scalar has been created with the same
       name as the requested matrix, delete the matrix
    */
    if (name_is_variable(name, pdinfo)) {
	delete_matrix_by_name(name);
	return NULL;
    }

    for (i=0; i<n_matrices; i++) {
	if (!strcmp(name, matrices[i]->name) &&
	    matrices[i]->level == level) {
	    return matrices[i]->M;
	}
    }

    return NULL;
}

static user_matrix *get_user_matrix_by_name (const char *name)
{
    int level = gretl_function_stack_depth();
    int i;

    for (i=0; i<n_matrices; i++) {
	if (!strcmp(name, matrices[i]->name) &&
	    matrices[i]->level == level) {
	    return matrices[i];
	}
    }

    return NULL;
}

static user_matrix *get_user_matrix_by_data (const gretl_matrix *M)
{
    int level = gretl_function_stack_depth();
    int i;

    for (i=0; i<n_matrices; i++) {
	if (matrices[i]->M == M && matrices[i]->level == level) {
	    return matrices[i];
	}
    }

    return NULL;
}

/**
 * get_matrix_by_name:
 * @name: name of the matrix.
 * @pdinfo: dataset information.
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_by_name (const char *name, const DATAINFO *pdinfo)
{
    return real_get_matrix_by_name(name, LEVEL_AUTO, pdinfo);
}

/**
 * get_matrix_transpose_by_name:
 * @name: name of the matrix.
 * @pdinfo: dataset information.
 *
 * If @name ends with the transpose symbol ('), and if
 * a stacked matrix M exists with a name equal to @name
 * without the trailing symbol, return an allocated
 * copy of the transpose of M.
 *
 * Returns: newly allocated matrix, or %NULL if not found.
 */

gretl_matrix *
get_matrix_transpose_by_name (const char *name, const DATAINFO *pdinfo)
{
    gretl_matrix *M = NULL;
    char test[VNAMELEN];

    *test = '\0';
    strncat(test, name, VNAMELEN - 1);

    if (test[strlen(test) - 1] == '\'') {
	test[strlen(test) - 1] = '\0';
	M = get_matrix_by_name(test, pdinfo);
	if (M != NULL) {
	    return gretl_matrix_copy_transpose(M);
	}
    }

    return NULL;
}

/**
 * get_matrix_by_name_at_level:
 * @name: name of the matrix.
 * @level: level of function execution at which to search.
 * @pdinfo: dataset information.
 *
 * Looks up a user-defined matrix by name, at the given
 * level fo function execution.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_by_name_at_level (const char *name, int level,
					   const DATAINFO *pdinfo)
{
    return real_get_matrix_by_name(name, level, pdinfo);
}

/**
 * copy_named_matrix_as:
 * @orig: the name of the original matrix.
 * @new: the name to be given to the copy.
 *
 * If a saved matrix is found by the name @orig, a copy of
 * this matrix is added to the stack of saved matrices under the
 * name @new.  This is intended for use when a matrix is given
 * as the argument to a user-defined function: it is copied
 * under the name assigned by the function's parameter list.
 *
 * Returns: 0 on success, non-zero on error.
 */

int copy_named_matrix_as (const char *orig, const char *new)
{
    user_matrix *u;
    int err = 0;

    u = get_user_matrix_by_name(orig);
    if (u == NULL) {
	err = 1;
    } else {
	gretl_matrix *M = gretl_matrix_copy(u->M);

	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    err = add_user_matrix(M, new);
	}
	if (!err) {
	    /* for use in functions: increment level of last-added matrix */
	    u = matrices[n_matrices - 1];
	    u->level += 1;
	}
    }

    return err;
}

/**
 * user_matrix_set_name_and_level:
 * @m: matrix to be reconfigured.
 * @name: new name to be given to matrix.
 * @level: function execution level to be assigned to matrix.
 *
 * If matrix @m is found on the stack of saved matrices, change
 * its name to @newname and change its level to @level.
 *
 * Returns: 0 on success, 1 if the matrix is not found.
 */

int user_matrix_set_name_and_level (const gretl_matrix *M, char *name, 
				    int level)
{
    user_matrix *u = get_user_matrix_by_data(M);
    int err = 0;

    if (u == NULL) {
	err = 1;
    } else {
	*u->name = '\0';
	strncat(u->name, name, VNAMELEN - 1);
	u->level = level;
    }

    return err;
}

/**
 * named_matrix_get_variable:
 * @mspec: name of matrix, possibly followed by selection string
 * in square backets (e.g. "m[,1]" to get column 1 of m).
 * @Z: data array.
 * @pdinfo: dataset information.
 * @px: location to receive allocated series (array).
 * @plen: on input, if this contains 1 then only a 1 x 1 matrix
 * is acceptable; on output, the length of the array in @px.
 *
 * If there exists a matrix of the given name at the current
 * level of function execution, try to assign to @x the
 * selection specified in @mspec.  This works only if (a)
 * the matrix is in fact a vector of length equal to either the 
 * full length of the dataset or the length of the current sample 
 * range, or (b) the selection string "extracts" such a vector,
 * or (c) the specified matrix is 1 x 1, in effect yielding a
 * scalar result (array of length 1).
 *
 * If the sample is currently restricted while the selection has 
 * a number of elements equal to the full length of the dataset, 
 * the "out of sample" observations are set to #NADBL.
 *
 * Returns: 0 on success or non-zero if no matrix is found, or
 * if the matrix selection does not have the required number of
 * elements.
 */

int named_matrix_get_variable (const char *mspec, 
			       double ***pZ, DATAINFO *pdinfo,
			       double **px, int *plen)
{
    double *x = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *S = NULL;
    int sn = pdinfo->t2 - pdinfo->t1 + 1;
    int i, len = 0;
    int err = 0;

    if (strchr(mspec, '[')) {
	S = user_matrix_get_slice(mspec, pZ, pdinfo, &err);
	if (!err) {
	    M = S;
	}
    } else if (strchr(mspec, '\'')) {
	/* we won't handle transposes in this context */
	return 1;
    } else {
	M = get_matrix_by_name(mspec, pdinfo);
    }

    if (M != NULL) {
	len = gretl_vector_get_length(M);
    } else {
	err = E_UNKVAR;
    }

    if (!err) {
	if (*plen == 1 && len != 1) {
	    err = 1;
	} else if (len != 1 && len != pdinfo->n && len != sn) {
	    err = E_NONCONF;
	}
    }

    if (!err) {
	if (len == 1) {
	    x = malloc(sizeof *x);
	    if (x == NULL) {
		err = E_ALLOC;
	    } else {
		*x = gretl_vector_get(M, 0);
		*px = x;
	    }
	} else {
	    int mi;

	    x = malloc(pdinfo->n * sizeof *x);
	    if (x == NULL) {
		err = E_ALLOC;
	    } else {
		for (i=0; i<pdinfo->n; i++) {
		    x[i] = NADBL;
		}
		for (i=0; i<sn; i++) {
		    mi = (len == sn)? i : i + pdinfo->t1;
		    x[i + pdinfo->t1] = gretl_vector_get(M, mi);
		}
		*px = x;
	    }
	}
    }

    *plen = len;

    if (S != NULL) {
	gretl_matrix_free(S);
    }   

    return err;
}

static int *slice_from_index_vector (const gretl_matrix *v, int *err)
{
    int sn = gretl_vector_get_length(v);
    int *slice = NULL;
    int i;

    if (sn > 0) {
	slice = gretl_list_new(sn);
	if (slice == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=0; i<sn; i++) {
		slice[i+1] = gretl_vector_get(v, i);
	    }
	}
    } else {
	*err = E_DATA;
    }

    return slice;
}

static int *slice_from_scalar (const char *s, int *err)
{
    int *slice = NULL;
    double x = NADBL;

    if (gZ == NULL || gdinfo == NULL) {
	*err = E_DATA;
	return NULL;
    }

    *err = get_generated_value(s, &x, gZ, gdinfo, 0);

#if MDEBUG
    fprintf(stderr, "slice_from_scalar: tried s='%s', got err = %d\n", 
	    s, *err);
#endif

    if (!*err && (na(x) || x < 0)) {
	*err = E_DATA;
    }

    if (!*err) {
	slice = gretl_list_new(1);
	slice[1] = x;
    }

    return slice;
}

/* A "slice" is a gretl list: the first element holds the count of the
   following elements, and the following elements indicate from which
   row (or column) of the source matrix to draw the successive rows
   (columns) of the submatrix.  A NULL value is also OK, indicating
   that all rows (columns) should be used.
*/

static int *parse_slice_spec (const char *s, int n, int *err)
{
    int *slice = NULL;
    int i, sn = 0;

    if (*s == '\0') {
	/* null spec: use all rows or columns */
	return NULL;
    }

    if (isalpha(*s) && strchr(s, ':') == NULL) {
	/* could be the name of an index matrix? */
	gretl_matrix *v = get_matrix_by_name(s, gdinfo);

#if MDEBUG
	fprintf(stderr, "index matrix? s='%s', v=%p\n", s, (void *) v);
#endif
	if (v != NULL) {
	    slice = slice_from_index_vector(v, err);
	} else {
	    /* scalar var or formula? */
	    slice = slice_from_scalar(s, err);
	}
    }

    if (slice == NULL && !*err) {
	/* not found yet: keep trying, for p:q or plain p */
	char f1[32], f2[32]; /* arbitrary? */
	double x[2];

	if (sscanf(s, "%31[^:]:%31s", f1, f2) == 2) {
	    *err = get_generated_value(f1, &x[0], gZ, gdinfo, 0);
	    if (!*err) {
		*err = get_generated_value(f2, &x[1], gZ, gdinfo, 0);
	    }
	    if (*err && (na(x[0]) || na(x[1]))) {
		*err = 1;
	    }
	    if (!*err && (x[0] < 1 || x[1] < 1 || x[1] < x[0])) {
		*err = 1;
	    }
	    if (!*err) {
		sn = x[1] - x[0] + 1;
	    }
	} else if (sscanf(s, "%31s", f1) == 1) {
	    *err = get_generated_value(f1, &x[0], gZ, gdinfo, 0);
	    if (!*err && (na(x[0]) || x[0] < 1)) {
		*err = 1;
	    }
	    if (!*err) {
		sn = 1;
	    }
	} else {
	    *err = 1;
	}

	if (!*err) {
	    slice = gretl_list_new(sn);
	    if (slice == NULL) {
		*err = E_ALLOC;
	    } else {
		for (i=0; i<sn; i++) {
		    slice[i+1] = x[0] + i;
		}
	    }
	}
    }

    if (slice != NULL) {
#if MDEBUG
	printlist(slice, "slice");
#endif
	for (i=1; i<=slice[0] && !*err; i++) {
	    if (slice[i] < 0 || slice[i] > n) {
		fprintf(stderr, "index value %d is out of bounds\n", 
			slice[i]);
		*err = 1;
	    }
	}
    }

    if (*err && slice != NULL) {
	free(slice);
	slice = NULL;
    }

    return slice;    
}

enum {
    RSLICE,
    CSLICE,
    VSLICE
};

static int get_slice_string (const char *s, char *spec, int i)
{
    int err = 0;

    *spec = '\0';

    if (i == RSLICE) {
	s = strchr(s, '[');
	if (s == NULL || strchr(s, ',') == NULL) {
	    err = 1;
	} else {
	    sscanf(s + 1, "%15[^,]]", spec);
	}
    } else if (i == CSLICE) {
	s = strchr(s, ',');
	if (s == NULL || strchr(s, ']') == NULL) {
	    err = 1;
	} else {	
	    sscanf(s + 1, "%15[^]]", spec);
	}	    
    } else if (i == VSLICE) {
	s = strchr(s, '[');
	if (s == NULL || strchr(s, ']') == NULL) {
	    err = 1;
	} else {
	    sscanf(s + 1, "%15[^]]", spec);
	}
    }

#if MDEBUG
    fprintf(stderr, "get_slice_string: spec = '%s', err = %d\n",
	    spec, err);
#endif

    return err;
}

static int 
make_slices (const char *s, int m, int n, int **rslice, int **cslice)
{
    char spec[VNAMELEN];
    int err = 0;

    *rslice = *cslice = NULL;

    if ((m == 1 || n == 1) && strchr(s, ',') == NULL) {
	/* vector: a single slice is acceptable */
	err = get_slice_string(s, spec, VSLICE);
	if (!err) {
	    if (m == 1) {
		*cslice = parse_slice_spec(spec, n, &err);
	    } else {
		*rslice = parse_slice_spec(spec, m, &err);
	    }
	}
    } else {
	/* not a vector: must have both slice strings "in principle",
	   even if one or the other is empty */
	err = get_slice_string(s, spec, RSLICE);
	if (!err) {
	    *rslice = parse_slice_spec(spec, m, &err);
	}

	if (!err) {
	    err = get_slice_string(s, spec, CSLICE);
	    if (!err) {
		*cslice = parse_slice_spec(spec, n, &err);
	    } 
	}
    }   

    if (err) {
	free(*rslice);
	*rslice = NULL;
	free(*cslice);
	*cslice = NULL;
    }

    return err;
}

static gretl_matrix *
try_diagonal (const gretl_matrix *M, const char *s, int *err)
{
    gretl_matrix *d = NULL;

    s = strchr(s, '[');

    if (s != NULL) {
	char test[5];

	if (sscanf(s + 1, "%4[^]]]", test)) {
	    if (!strcmp(test, "diag")) {
		d = gretl_matrix_get_diagonal(M, err);
	    }
	}
    }

    return d;
}

/* This supports the extraction of scalars, vectors or sub-matrices.
   E.g. B[1,2] extracts a scalar, B[1,] extracts a row vector, and
   B[,2] extracts a column vector.  B[1:2,2:3] extracts a sub-matrix
   composed of the intersection of rows 1 and 2 with columns 2 and 3.
   As a special case B[diag] extracts the diagonal from a square
   matrix B.
*/

gretl_matrix *
matrix_get_submatrix (const gretl_matrix *M, const char *s, 
		      double ***pZ, DATAINFO *pdinfo, 
		      int *err)
{
    gretl_matrix *S;
    const char *p;
    int *rslice = NULL;
    int *cslice = NULL;
    int m = gretl_matrix_rows(M);
    int n = gretl_matrix_cols(M);
    int nr, nc;

    /* the selection string should end with ']' */
    p = strrchr(s, ']');
    if (p == NULL || *(p+1) != '\0') {
	*err = E_SYNTAX;
	return NULL;
    }    
    
    /* special case of "A[diag]": get the diagonal as vector */
    if (strstr(s, "diag")) {
	S = try_diagonal(M, s, err);
	if (S != NULL || *err != 0) {
	    return S;
	}
    }

    usermat_publish_dataset(pZ, pdinfo);

    *err = make_slices(s, m, n, &rslice, &cslice);
    if (*err) {
	return NULL;
    }

    nr = (rslice == NULL)? m : rslice[0];
    nc = (cslice == NULL)? n : cslice[0];

    S = gretl_matrix_alloc(nr, nc);
    if (S == NULL) {
	*err = E_ALLOC;	
    }

    if (S != NULL) {
	int i, j, k, l;
	int mi, mj;
	double x;

	k = 0;
	for (i=0; i<nr && !*err; i++) {
	    mi = (rslice == NULL)? k++ : rslice[i+1] - 1;
	    if (mi < 0) {
		*err = 1;
	    }
	    l = 0;
	    for (j=0; j<nc && !*err; j++) {
		mj = (cslice == NULL)? l++ : cslice[j+1] - 1;
		if (mj < 0) {
		    *err = 1;
		} else {
		    x = gretl_matrix_get(M, mi, mj);
		    gretl_matrix_set(S, i, j, x);
		}
	    }
	}
    }

    usermat_unpublish_dataset(pZ);

    free(rslice);
    free(cslice);

    if (*err) {
	gretl_matrix_free(S);
	S = NULL;
    }
	
    return S;
}

/**
 * user_matrix_get_slice:
 * @s: string specifying a sub-matrix.
 * @Z: data array.
 * @pdinfo: dataset information.
 * @err: location to receive error code.
 *
 * If @s specifies a valid "slice" of an existing named
 * matrix at the current level of function execution, 
 * constructs a newly allocated sub-matrix.
 *
 * Returns: allocated sub-matrix on success, %NULL on failure.
 */

gretl_matrix *user_matrix_get_slice (const char *s, 
				     double ***pZ, 
				     DATAINFO *pdinfo,
				     int *err)
{
    gretl_matrix *M = NULL;
    gretl_matrix *S = NULL;
    char test[VNAMELEN];
    int len = strcspn(s, "[");

    if (len < VNAMELEN) {
	*test = '\0';
	strncat(test, s, len);
	M = get_matrix_by_name(test, pdinfo);
	if (M != NULL) {
	    S = matrix_get_submatrix(M, s, pZ, pdinfo, err);
	} else {
	    *err = E_UNKVAR;
	}
    }

    return S;
}

/**
 * add_or_replace_user_matrix:
 * @M: gretl matrix.
 * @name: name for the matrix.
 * @mask: submatrix specification (or empty string).
 * @R: location to receive address of matrix that was
 * replaced, if any (or %NULL).
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 *
 * Checks whether a matrix of the given @name already exists.
 * If so, the original matrix is replaced by @M; if not, the
 * the matrix @M is added to the stack of user-defined
 * matrices.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int add_or_replace_user_matrix (gretl_matrix *M, const char *name,
				const char *mask, gretl_matrix **R,
				double ***pZ, DATAINFO *pdinfo)
{
    user_matrix *u;
    int err = 0;

    u = get_user_matrix_by_name(name);

    if (u != NULL) {
	if (pZ != NULL && pdinfo != NULL) {
	    usermat_publish_dataset(pZ, pdinfo);
	}
	err = replace_user_matrix(u, M, R, mask);
	usermat_unpublish_dataset(pZ);
    } else {
	int v;

	err = add_user_matrix(M, name);
	/* if we created a matrix with the same name as an existing
	   scalar, delete the scalar */
	v = varindex(pdinfo, name);
	if (v < pdinfo->v && !pdinfo->vector[v]) {
	    dataset_drop_variable(v, pZ, pdinfo);
	}
    }

    return err;
}

static void destroy_user_matrix (user_matrix *u)
{
    if (u == NULL) {
	return;
    }

#if MDEBUG
    fprintf(stderr, "destroy_user_matrix: freeing matrix at %p...", 
	    (void *) u->M);
#endif
    gretl_matrix_free(u->M);
    free(u);
#if MDEBUG
    fprintf(stderr, " done\n");
#endif
}

/**
 * destroy_user_matrices_at_level:
 * @level: stack level of function execution.
 *
 * Destroys and removes from the stack of user matrices all
 * matrices that were created at the given @level.  This is 
 * part of the cleanup that is performed when a user-defined
 * function terminates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_matrices_at_level (int level)
{
    user_matrix **tmp;
    int i, j, nm = 0;
    int err = 0;

#if MDEBUG
    fprintf(stderr, "destroy_user_matrices_at_level: level = %d, "
	    "n_matrices = %d\n", level, n_matrices);
#endif

    for (i=0; i<n_matrices; i++) {
	if (matrices[i] == NULL) {
	    break;
	}
	if (matrices[i]->level == level) {
#if MDEBUG
	    fprintf(stderr, "destroying matrix[%d] (M at %p)\n",
		    i, (void *) matrices[i]->M);
#endif
	    destroy_user_matrix(matrices[i]);
	    for (j=i; j<n_matrices - 1; j++) {
		matrices[j] = matrices[j+1];
	    }
	    matrices[n_matrices - 1] = NULL;
	} else {
	    nm++;
	}
    }

    if (nm < n_matrices) {
	n_matrices = nm;
	if (nm == 0) {
	    free(matrices);
	    matrices = NULL;
	} else {
	    tmp = realloc(matrices, nm * sizeof *tmp);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		matrices = tmp;
	    }
	}
    }

    return err;
}

/**
 * destroy_user_matrices:
 *
 * Frees all resources associated with the stack of user-
 * defined matrices.
 */

void destroy_user_matrices (void)
{
    int i;

#if MDEBUG
    fprintf(stderr, "destroy_user_matrices called, n_matrices = %d\n",
	    n_matrices);
#endif

    if (matrices == NULL) {
	return;
    }

    for (i=0; i<n_matrices; i++) {
#if MDEBUG
	fprintf(stderr, "destroying user_matrix %d (%s) at %p\n", i,
		matrices[i]->name, (void *) matrices[i]);
#endif
	destroy_user_matrix(matrices[i]);
    }

    free(matrices);
    matrices = NULL;
    n_matrices = 0;
}

static int real_delete_user_matrix (user_matrix *u)
{
    int err = 0;

    if (u == NULL) {
	err = E_DATA;
    } else {
	int i, j, nm = n_matrices - 1;

	for (i=0; i<n_matrices; i++) {
	    if (matrices[i] == u) {
		destroy_user_matrix(matrices[i]);
		for (j=i; j<n_matrices - 1; j++) {
		    matrices[j] = matrices[j+1];
		}
		matrices[nm] = NULL;
		break;
	    }
	} 

	if (nm == 0) {
	    free(matrices);
	    matrices = NULL;
	} else {
	    user_matrix **tmp = realloc(matrices, nm * sizeof *tmp);

	    n_matrices = nm;
	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		matrices = tmp;
	    }
	}
    }

    return err;
}

static int delete_matrix_by_name (const char *name)
{
    user_matrix *u = get_user_matrix_by_name(name);

    return real_delete_user_matrix(u);
}

static int first_field_is_series (const char *s, const DATAINFO *pdinfo)
{
    char word[16];
    int v, ret = 0;

    s += strspn(s, " {");

    *word = '\0';
    sscanf(s, "%15[^,;} ]", word);
    v = varindex(pdinfo, word);
    if (v < pdinfo->v && pdinfo->vector[v]) {
	ret = 1;
    }

    return ret;
}

static int 
get_rows_cols (const char *str, const DATAINFO *pdinfo,
	       int *r, int *c, int *series)
{
    const char *p = str;
    char sepstr[2] = ";";
    char sep = ';';
    char *s;
    int nf0 = 0, nf1;
    int i, len;
    int err = 0;

    if (first_field_is_series(p, pdinfo)) {
	*series = 1;
	*sepstr = ',';
	sep = ',';
    }

    *r = 1;
    while (*p) {
	if (*p == sep) {
	    *r += 1;
	}
	p++;
    }

    for (i=0; i<*r && !err; i++) {
	len = strcspn(str, sepstr);

	s = gretl_strndup(str, len);
	if (s == NULL) {
	    err = E_ALLOC;
	    break;
	}

	/* clean up */
	charsub(s, ',', ' ');
	charsub(s, '}', ' ');
	charsub(s, '\'', ' ');

	nf1 = count_fields(s);
	if (i > 0 && nf1 != nf0) {
	    strcpy(gretl_errmsg, "Inconsistent matrix specification");
	    err = 1;
	    break;
	}

	nf0 = nf1;
	str += len + 1;
	free(s);
    }

    if (!err) {
	if (*series) {
	    *c = *r;
	    *r = nf0;
	} else {
	    *c = nf0;
	}
    }

    /* FIXME series and matrix orientation */

    return err;
}

static int 
get_varnum (const char **s, const DATAINFO *pdinfo, int *err)
{
    char vname[12];
    int v = 0;

    if (sscanf(*s, "%11s", vname)) {
	int len = strlen(vname);
	int vlen = gretl_varchar_spn(vname);

	if (vlen < len) {
	    vname[vlen] = '\0';
	}

	v = varindex(pdinfo, vname);
	if (v < pdinfo->v) {
	    *s += len;
	} else {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	    *err = 1;
	} 
    } else {
	*err = 1;
    }

    return v;
}

static double get_var_double (const char **s, const double **Z,
			      const DATAINFO *pdinfo, int *err)
{
    char vname[VNAMELEN];
    double x = NADBL;
    int v, len = gretl_varchar_spn(*s);

    if (len > VNAMELEN - 1) {
	*err = 1;
    } else {
	*vname = '\0';
	strncat(vname, *s, len);
	v = varindex(pdinfo, vname);
	if (v == pdinfo->v) {
	    *err = E_UNKVAR;
	} else if (pdinfo->vector[v]) {
	    *err = E_DATA;
	} else {
	    x = Z[v][0];
	    if (na(x)) {
		*err = E_MISSDATA;
	    } else {
		*s += len;
	    }
	}
    }

    return x;
}

static double get_numeric_double (const char **s, int *err)
{
    double x = NADBL;
    char *p;

    x = strtod(*s, &p);

    if (!strcmp(*s, p)) {
	sprintf(gretl_errmsg, _("'%s' -- no numeric conversion performed!"), *s);
	*err = 1;
    } else if (*p != '\0' && *p != ',' && *p != ';' && *p != ' ' && *p != '}') {
	if (isprint(*p)) {
	    sprintf(gretl_errmsg, _("Extraneous character '%c' in data"), *p);
	} else {
	    sprintf(gretl_errmsg, _("Extraneous character (0x%x) in data"), *p);
	}
	*err = 1;
    } else if (errno == ERANGE) {
	sprintf(gretl_errmsg, _("'%s' -- number out of range!"), *s);
	*err = 1;
    }

    *s = p;

    return x;
}

static int 
fill_matrix_from_scalars (gretl_matrix *M, const char *s, 
			  int r, int c, int transp,
			  const double **Z, const DATAINFO *pdinfo)
{
    double x;
    int i, j;
    int err = 0;

    gretl_push_c_numeric_locale();

    if (transp) {
	for (j=0; j<c && !err; j++) {
	    for (i=0; i<r && !err; i++) {
		if (isalpha(*s)) {
		    x = get_var_double(&s, Z, pdinfo, &err);
		} else {
		    x = get_numeric_double(&s, &err);
		}
		if (!err) {
		    gretl_matrix_set(M, i, j, x);
		    s += strspn(s, " ,;");
		}
	    }
	}
    } else {
	for (i=0; i<r && !err; i++) {
	    for (j=0; j<c && !err; j++) {
		if (isalpha(*s)) {
		    x = get_var_double(&s, Z, pdinfo, &err);
		} else {
		    x = get_numeric_double(&s, &err);
		}
		if (!err) {
		    gretl_matrix_set(M, i, j, x);	
		    s += strspn(s, " ,;");
		}
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    return err;
}

/* Fill a matrix with values from one or more series: since
   gretl's matrix methods cannot handle missing values, it is
   an error if any missing values are encountered in the
   given series.  Each series occupies a row by default
   (unless transp is non-zero).
*/

static int fill_matrix_from_series (gretl_matrix *M, const char *s,
				    int r, int c, int transp,
				    const double **Z, const DATAINFO *pdinfo)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int nvr = r, nvc = c;
    int i, j, k, t, v;
    double x;
    int err = 0;

    s += strspn(s, " ");

    if (!transp) {
	nvr /= T;
	for (j=0; j<c && !err; j++) {
	    for (i=0; i<nvr && !err; i++) {
		v = get_varnum(&s, pdinfo, &err);
		if (!err) {
		    for (t=0; t<T; t++) {
			if (pdinfo->vector[v]) {
			    x = Z[v][t + pdinfo->t1];
			} else {
			    x = Z[v][0];
			}
			if (na(x)) {
			    err = E_MISSDATA;
			} else {
			    k = i * T + t;
			    gretl_matrix_set(M, k, j, x);
			}
		    }
		    s += strspn(s, " ,;");
		}
	    }
	}
    } else {
	nvc /= T;
	for (i=0; i<r && !err; i++) {
	    for (j=0; j<nvc && !err; j++) {
		v = get_varnum(&s, pdinfo, &err);
		if (!err) {
		    for (t=0; t<T; t++) {
			if (pdinfo->vector[v]) {
			    x = Z[v][t + pdinfo->t1]; 
			} else {
			    x = Z[v][0];
			}
			if (na(x)) {
			    err = E_MISSDATA;
			} else {			
			    k = j * T + t;
			    gretl_matrix_set(M, i, k, x);
			}
		    }
		    s += strspn(s, " ,;");
		}
	    }
	}
    } 

    return err;
}

static int matrix_genr (const char *name, const char *mask, 
			const char *s, double ***pZ, 
			DATAINFO *pdinfo, PRN *prn)
{
    char genline[MAXLINE];
    int err;

    if (strlen(name) + strlen(mask) + strlen(s) + 7 > MAXLINE - 1) {
	err = 1;
    } else {
	sprintf(genline, "genr %s%s %s", name, mask, s);
	err = generate(genline, pZ, pdinfo, OPT_M, prn);
    }

    return err;
}

gretl_matrix *fill_matrix_from_list (const char *s, const double **Z,
				     const DATAINFO *pdinfo, int transp,
				     int *err)
{
    gretl_matrix *M = NULL;
    char word[VNAMELEN];
    char *mask = NULL;
    const int *list;
    int len;

    while (isspace(*s)) s++;

    len = gretl_varchar_spn(s);
    if (len == 0 || len > VNAMELEN - 1) {
	return NULL;
    }

    *word = '\0';
    strncat(word, s, len);
    list = get_list_by_name(word);
    if (list == NULL) {
	return NULL;
    }

    M = gretl_matrix_data_subset(list, Z, pdinfo->t1, pdinfo->t2, &mask);

    if (M == NULL) {
	*err = 1;
    }

    if (mask != NULL) {
	*err = E_MISSDATA;
	free(mask);
	gretl_matrix_free(M);
	M = NULL;
    }

    if (M != NULL && transp) {
	gretl_matrix *R = gretl_matrix_copy_transpose(M);

	if (R == NULL) {
	    *err = E_ALLOC;
	    gretl_matrix_free(M);
	    M = NULL;
	} else {
	    gretl_matrix_free(M);
	    M = R;
	}
    }

    return M;
}

static int name_is_series (const char *name, const DATAINFO *pdinfo)
{
    int v = varindex(pdinfo, name);
    int ret = 0;

    if (v < pdinfo->v && pdinfo->vector[v]) {
	sprintf(gretl_errmsg, _("'%s' is the name of a data series"), name);
	ret = 1;
    }

    return ret;
}

/* Currently we can create a user matrix in any one of four ways (but
   we can't mix these in a single matrix specification).

   1. Full specification of scalar elements, either numerical values or
      by reference to scalar variables.

   2. Specification of individual data series to place in the matrix.

   3. Specification of one named list of variables to place in matrix.

   4. Use of a "genr"-type expression referring to existing matrices
      and/or the matrix-from-scratch functions such as I(n), which 
      generates an n x n identity matrix.

   If the name supplied for a matrix is already taken by an "ordinary"
   series variable, the attempt to create a matrix of the same name
   fails with an error message.
*/

static int create_matrix (const char *name, const char *mask,
			  const char *s, double ***pZ, DATAINFO *pdinfo,
			  PRN *prn)
{
    gretl_matrix *M = NULL;
    char *p;
    int nm = n_matrices;
    int series = 0;
    int transp = 0;
    int r = 0, c = 0;
    int err = 0;

    if (name_is_series(name, pdinfo)) {
	/* can't overwrite data series with matrix */
	return E_TYPES;
    }

    p = strchr(s, '{');
    if (p == NULL) {
	err = matrix_genr(name, mask, s, pZ, pdinfo, prn);
	goto finalize;
    }
    s = p + 1;

    if (!err) {
	p = strchr(s, '}');
	if (p == NULL) {
	    err = 1;
	}
	if (*(p+1) == '\'') {
	    transp = 1;
	}
    }

    if (!err) {
	M = fill_matrix_from_list(s, (const double **) *pZ, pdinfo, 
				  transp, &err);
	if (err) {
	    goto finalize;
	} else if (M != NULL) {
	    err = add_or_replace_user_matrix(M, name, mask, NULL,
					     pZ, pdinfo);
	    goto finalize;
	}
    }

    if (!err) {
	err = get_rows_cols(s, pdinfo, &r, &c, &series);
	if (!err && c == 0) {
	    err = 1;
	}
    }

    if (!err && series) {
	r *= pdinfo->t2 - pdinfo->t1 + 1;
    }

    if (!err && transp) {
	int tmp = r;

	r = c;
	c = tmp;
    }

#if MDEBUG
    fprintf(stderr, "r=%d, c=%d, transp=%d, series=%d, s='%s'\n",
	    r, c, transp, series, s);
#endif

    if (!err) {
	M = gretl_matrix_alloc(r, c);
	if (M == NULL) {
	   err = E_ALLOC;
	} 
    }

    if (!err) {
	if (series) {
	    err = fill_matrix_from_series(M, s, r, c, transp, 
					  (const double **) *pZ, 
					  pdinfo);
	} else {
	    err = fill_matrix_from_scalars(M, s, r, c, transp,
					   (const double **) *pZ,
					   pdinfo);
	}
    }
    
    if (!err) {
	err = add_or_replace_user_matrix(M, name, mask, NULL,
					 pZ, pdinfo);
    }

 finalize:

    if (!err) {
	if (gretl_messages_on()) {
	    if (n_matrices > nm) {
		pprintf(prn, "Added matrix '%s'\n", name);
	    } else {
		pprintf(prn, "Replaced matrix '%s'\n", name);
	    }
	}
    } else {
	pprintf(prn, "Error adding matrix '%s'\n", name);
    }

    return err;
}

static int 
print_matrix_by_name (const char *name, const char *mask, 
		      double ***pZ, DATAINFO *pdinfo,
		      PRN *prn)
{
    gretl_matrix *M;
    gretl_matrix *S;
    int transp = 0;
    int err = 0;

    M = get_matrix_by_name(name, pdinfo);

    if (M == NULL) {
	transp = 1;
	M = get_matrix_transpose_by_name(name, pdinfo);
    }

    if (M == NULL) {
	pprintf(prn, _("'%s': no such matrix\n"), name);
	err = 1;
    } else {
	if (*mask != '\0') {
	    S = matrix_get_submatrix(M, mask, pZ, pdinfo, &err);
	    if (!err) {
		char mspec[96];

		sprintf(mspec, "%s%s", name, mask);
		gretl_matrix_print_to_prn(S, mspec, prn);
	    }
	    gretl_matrix_free(S);
	} else {
	    gretl_matrix_print_to_prn(M, name, prn);
	}
	if (transp) {
	    /* we got a transpose, created on the fly */
	    gretl_matrix_free(M);
	}
    }

    return err;
}

static int 
print_matrix_address_by_name (const char *name, const DATAINFO *pdinfo, PRN *prn)
{
    const gretl_matrix *M = get_matrix_by_name(name, pdinfo);
    int err = 0;

    if (M == NULL) {
	pprintf(prn, _("'%s': no such matrix\n"), name);
    } else {
	pprintf(prn, "matrix '%s' is at %p\n", name, (void *) M);
    }

    return err;
}

/**
 * matrix_command:
 * @s: string that specifies matrix command.
 * @prn: pointer to printing struct.
 *
 * To be written.
 * 
 * Returns: 0 on success, non-zero code on failure.
 */

int matrix_command (const char *line, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    char mask[48] = {0};
    char name[48];
    char word[9];
    char *p;
    int err = 0;

    if (!strncmp(line, "matrix ", 7)) line += 7;
    while (isspace(*line)) line++;

    if (!sscanf(line, "%47[^ =]", name)) {
	return 1;
    } 

    line += strlen(name);
    while (isspace(*line)) line++;

    /* extract submatrix "mask", if any */
    p = strchr(name, '[');
    if (p != NULL) {
	strncat(mask, p, 47);
	*p = '\0';
	if (get_matrix_by_name(name, pdinfo) == NULL) {
	    return E_UNKVAR;
	}
    } 

    if (*line == '=') {
	/* defining a matrix */
	err = create_matrix(name, mask, line, pZ, pdinfo, prn);
    } else {
	*word = '\0';
	sscanf(line, "%8s", word);
	if (*word == '\0' || !strcmp(word, "print")) {
	    err = print_matrix_by_name(name, mask, pZ, pdinfo, prn);
	} else if (!strcmp(word, "delete")) {
	    err = delete_matrix_by_name(name);
	} else if (!strcmp(word, "address")) {
	    err = print_matrix_address_by_name(name, pdinfo, prn);
	} else {
	    /* no other commands available yet */
	    err = 1;
	}
    } 

    return err;
}

static gretl_matrix *
matrix_test_equality (const gretl_matrix *A, const gretl_matrix *B, int *err)
{
    gretl_matrix *C = NULL;
    int eq;

    eq = gretl_matrices_are_equal(A, B, err);
    if (!*err) {
	C = gretl_matrix_from_scalar((double) eq);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
    }

    return C;
}

#if MDEBUG
static void print_calc_input_info (gretl_matrix *A, gretl_matrix *B, char op)
{
    fprintf(stderr, "matrix_calc_AB: A = %p", (void *) A);
    if (A != NULL) {
	fprintf(stderr, " (%dx%d)", gretl_matrix_rows(A), gretl_matrix_cols(A));
    } 
    fprintf(stderr, ", B = %p", (void *) B);
    if (B != NULL) {
	fprintf(stderr, " (%dx%d), ", gretl_matrix_rows(B), gretl_matrix_cols(B));
    }
    if (isprint(op)) {
	fprintf(stderr, "op='%c'\n", op);
    } else {
	fprintf(stderr, "op=%d\n", (int) op);
    }
# if MDEBUG > 1
    if (A != NULL) gretl_matrix_print(A, "input A");
    if (B != NULL) gretl_matrix_print(B, "input B");
# endif
}
#endif

/* for use in genr, for matrices */

gretl_matrix *matrix_calc_AB (gretl_matrix *A, gretl_matrix *B, 
			      char op, int *err) 
{
    gretl_matrix *C = NULL;
    gretl_matrix *D = NULL;
    double x;
    int ra, ca;
    int rb, cb;

    *err = 0;

#if MDEBUG
    print_calc_input_info(A, B, op);
#endif

    switch (op) {
    case '\0':
	C = B;
	break;
    case '=':
	C = matrix_test_equality(A, B, err);
	break;
    case '+':
    case '-':
	C = gretl_matrix_copy(A);
	if (C == NULL) {
	    *err = E_ALLOC;
	} else if (op == '+') {
	    *err = gretl_matrix_add_to(C, B);
	} else {
	    *err = gretl_matrix_subtract_from(C, B);
	}
	break;
    case '~':
	/* column-wise concatenation */
	C = gretl_matrix_col_concat(A, B, err);
	break;
    case '*':
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	if (ra == 1 && ca == 1) {
	    C = gretl_matrix_copy(B);
#if MDEBUG
	    fprintf(stderr, " allocated copy 'C' of B (%dx%d) at %p\n", 
		    rb, cb, (void *) C);
#endif
	    if (C == NULL) {
		*err = E_ALLOC;	
	    } else {
		x = gretl_matrix_get(A, 0, 0);
		gretl_matrix_multiply_by_scalar(C, x);
	    }
	} else if (rb == 1 && cb == 1) {
	    C = gretl_matrix_copy(A);
#if MDEBUG
	    fprintf(stderr, " allocated copy 'C' of A (%dx%d) at %p\n", 
		    ra, ca, (void *) C);
#endif
	    if (C == NULL) {
		*err = E_ALLOC;	
	    } else {
		x = gretl_matrix_get(B, 0, 0);
		gretl_matrix_multiply_by_scalar(C, x);
	    }
	} else {
	    C = gretl_matrix_alloc(ra, cb);
#if MDEBUG
	    fprintf(stderr, " allocated new 'C' (%dx%d) at %p\n", 
		    ra, cb, (void *) C);
#endif
	    if (C == NULL) {
		*err = E_ALLOC;
	    } else {
		*err = gretl_matrix_multiply(A, B, C);
	    }
	}
	break;
    case '/':
	/* matrix "division" */
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	C = gretl_matrix_copy(A);
	if (C == NULL) {
	    *err = E_ALLOC;
	} else {
	    if (rb == 1 && cb == 1) {
		x = gretl_vector_get(B, 0);
		*err = gretl_matrix_divide_by_scalar(C, x);
	    } else {
		D = gretl_matrix_copy(B);
		if (D == NULL) {
		    gretl_matrix_free(C);
		    C = NULL;
		    *err = E_ALLOC;
		} else {	
		    *err = gretl_LU_solve(D, C);
		    gretl_matrix_free(D);
		}
	    }
	}
	break;
    case OP_DOTMULT:
	/* element-wise multiplication */
	C = gretl_matrix_dot_multiply(A, B, err);
	break;
    case OP_DOTDIV:
	/* element-wise division */
	C = gretl_matrix_dot_divide(A, B, err);
	break;
    case OP_DOTPOW:
	/* element-wise exponentiation */
	if (!gretl_matrix_is_scalar(B)) {
	    *err = E_NONCONF;
	} else {
	    C = gretl_matrix_copy(A);
	    if (C == NULL) {
		*err = E_ALLOC;
	    } else {
		x = gretl_matrix_get(B, 0, 0);
		gretl_matrix_dot_pow(C, x);
	    }
	}
	break;
    case OP_KRON:
	/* Kronecker product */
	C = gretl_matrix_kronecker_product(A, B);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
	break;
    default:
	*err = E_TYPES;
	break;
    } 

    if (*err && C != NULL) {
	gretl_matrix_free(C);
	C = NULL;
    }

    return C;
}

static double 
real_user_matrix_get_determinant (const gretl_matrix *m, int log, int *err)
{
    double d = NADBL;

    if (m != NULL) {
	gretl_matrix *tmp = gretl_matrix_copy(m);

	if (tmp != NULL) {
	    if (log) {
		d = gretl_matrix_log_determinant(tmp, err);
	    } else {
		d = gretl_matrix_determinant(tmp, err);
	    }
	    gretl_matrix_free(tmp);
	}
    }

    return d;
}

double user_matrix_get_determinant (const gretl_matrix *m, int *err)
{
    return real_user_matrix_get_determinant(m, 0, err);
}

double user_matrix_get_log_determinant (const gretl_matrix *m, int *err)
{
    return real_user_matrix_get_determinant(m, 1, err);
}

gretl_matrix *user_matrix_get_inverse (const gretl_matrix *m)
{
    gretl_matrix *R = NULL;

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R != NULL) {
	    if (gretl_invert_matrix(R)) {
		gretl_matrix_free(R);
		R = NULL;
	    }
	} 
    }

    if (R == NULL) {
	strcpy(gretl_errmsg, _("Matrix inversion failed"));
    }

    return R;
}

gretl_matrix *user_matrix_cholesky_decomp (const gretl_matrix *m)
{
    gretl_matrix *R = NULL;

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R != NULL) {
	    if (gretl_matrix_cholesky_decomp(R)) {
		gretl_matrix_free(R);
		R = NULL;
	    }
	} 
    }

    if (R == NULL) {
	strcpy(gretl_errmsg, _("Matrix decomposition failed"));
    }

    return R;
}

gretl_matrix *user_matrix_get_row_sum (const gretl_matrix *m)
{
    gretl_matrix *s = NULL;

    if (m != NULL) {
	s = gretl_matrix_row_sum(m);
    }

    return s;
}

gretl_matrix *user_matrix_get_column_sum (const gretl_matrix *m)
{
    gretl_matrix *s = NULL;

    if (m != NULL) {
	s = gretl_matrix_column_sum(m);
    }

    return s;
}

gretl_matrix *user_matrix_column_demean (const gretl_matrix *m)
{
    gretl_matrix *R = NULL;

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R != NULL) {
	    gretl_matrix_demean_by_column(R);
	} 
    }

    return R;
}

static int 
real_user_matrix_QR_decomp (const gretl_matrix *m, gretl_matrix **Q, 
			    gretl_matrix **R)
{
    int mc, err = 0;

    *Q = NULL;
    *R = NULL;

    if (m != NULL) {
	mc = gretl_matrix_cols(m);
	*Q = gretl_matrix_copy(m);
	*R = gretl_matrix_alloc(mc, mc);

	if (*Q == NULL || *R == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_QR_decomp(*Q, *R);
	} 
    }

    if (err) {
	strcpy(gretl_errmsg, _("Matrix decomposition failed"));
	gretl_matrix_free(*Q);
	gretl_matrix_free(*R);
	*Q = NULL;
	*R = NULL;
    }

    return err;
}

static int get_two_matrix_names (const char *s, char *lstr, char *rstr,
				 const DATAINFO *pdinfo)
{
    int err = 0;

    if (sscanf(s, "%15[^,],%15s", lstr, rstr) != 2) {
	err = 1;
    } else {
#if MDEBUG
	fprintf(stderr, "left-hand matrix = '%s'\n", lstr);
	fprintf(stderr, "right-hand matrix = '%s'\n", rstr);
#endif
	if (!strcmp(rstr, "null")) {
	    *rstr = 0;
	} 

	if (!err && name_is_series(lstr, pdinfo)) {
	    err = 1;
	}

	if (!err && *rstr && name_is_series(rstr, pdinfo)) {
	    err = 1;
	}
    }

    return err;
}

static int add_or_replace_aux_matrix (gretl_matrix *A,
				      const char *aname,
				      int old_nm,
				      double ***pZ,
				      DATAINFO *pdinfo,
				      PRN *prn)
{
    int err;

    err = add_or_replace_user_matrix(A, aname, NULL, NULL, pZ, pdinfo);

    if (!err && gretl_messages_on()) {
	if (n_matrices > old_nm) {
	    pprintf(prn, "Added matrix '%s'\n", aname);
	} else {
	    pprintf(prn, "Replaced matrix '%s'\n", aname);
	}
    }

    return err;
}

gretl_matrix *
user_matrix_QR_decomp (const char *str, double ***pZ, DATAINFO *pdinfo,
		       PRN *prn, int *err)
{
    int nm = n_matrices;
    gretl_matrix *M = NULL;
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    char qstr[VNAMELEN] = {0};
    char rstr[VNAMELEN] = {0};

    *err = get_two_matrix_names(str, qstr, rstr, pdinfo);

    if (!*err) {
	M = get_matrix_by_name(qstr, pdinfo);
	if (M == NULL) {
	    *err = E_UNKVAR;
	} else {
	    *err = real_user_matrix_QR_decomp(M, &Q, &R);
	}
	if (R != NULL) {
	    if (*rstr) {
		*err = add_or_replace_aux_matrix(R, rstr, nm, pZ, 
						 pdinfo, prn);
	    } else {
		gretl_matrix_free(R);
	    }
	}
    }

    return Q;
}

gretl_matrix *
user_matrix_eigen_analysis (const char *str, double ***pZ, DATAINFO *pdinfo,
			    PRN *prn, int *err, int symm)
{
    int nm = n_matrices;
    gretl_matrix *M = NULL;
    gretl_matrix *C = NULL;
    gretl_matrix *E = NULL;
    char lstr[VNAMELEN] = {0};
    char rstr[VNAMELEN] = {0};
    double *ev = NULL;
    int vecs = 0;

    *err = get_two_matrix_names(str, lstr, rstr, pdinfo);

    if (!*err) {
	int en = 0;

	vecs = (*rstr != 0);

	M = get_matrix_by_name(lstr, pdinfo);
	if (M == NULL) {
	    *err = E_UNKVAR;
	} else {
	    en = gretl_matrix_rows(M);
	    if (!symm) {
		en *= 2; /* allow for imaginary components */
	    }
	    C = gretl_matrix_copy(M);
	    if (C == NULL) {
		*err = E_ALLOC;
	    }
	}

	if (!*err) {
	    if (symm) {
		ev = gretl_symmetric_matrix_eigenvals(C, vecs, err);
	    } else {
		ev = gretl_general_matrix_eigenvals(C, vecs, err);
	    }
	}

	if (ev != NULL) {
	    E = gretl_vector_from_array(ev, en, GRETL_MOD_NONE);
	    free(ev);
	}

	if (!*err && vecs) {
	    *err = add_or_replace_aux_matrix(C, rstr, nm, pZ, 
					     pdinfo, prn);
	}

	if (!vecs) {
	    gretl_matrix_free(C);
	}
    }

    return E;
}

gretl_matrix *
user_matrix_get_transformation (const gretl_matrix *m, GretlMathFunc fn)
{
    gretl_matrix *R = NULL;

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R != NULL) {
	    if (gretl_matrix_transform_elements(R, fn)) {
		gretl_matrix_free(R);
		R = NULL;
	    }
	}
    }

    return R;
}

gretl_matrix *
user_matrix_get_sorted_vector (const gretl_matrix *m, int s, int *err)
{
    gretl_matrix *R = NULL;
    int len = gretl_vector_get_length(m);

    if (len == 0) {
	*err = E_NONCONF;
	return NULL;
    }

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R == NULL) {
	    *err = E_ALLOC;
	} else {
	    if (len > 1) {
		if (s == SORT_DESCENDING) {
		    qsort(R->val, len, sizeof *R->val, 
			  gretl_inverse_compare_doubles);
		} else {
		    qsort(R->val, len, sizeof *R->val, 
			  gretl_compare_doubles);
		}
	    }
	}
    }

    return R;
}  

/* move transpose symbol ' in front of parenthesized matrix expression
   (possibly with a leading function word), so genr can handle it as
   if it were a function
*/

int reposition_transpose_symbol (char *s)
{
    int pc, len = strlen(s);
    int offset;
    int i, j, k, sz;
    int err = 0;

    for (i=3; i<len; i++) {
	if (s[i] == '\'' && s[i-1] == ')') {
	    pc = sz = 1;
	    /* back up to matching left paren */
	    for (j=i-2; j>=0; j--) {
		if (s[j] == ')') {
		    pc++;
		} else if (s[j] == '(') {
		    pc--;
		}
		sz++;
		if (pc == 0) {
		    /* back up to start of function word, if any */
		    for (k=j-1; k>=0; k--) {
			if (!isalpha(s[k]) && !isdigit(s[k]) && s[k] != '_') {
			    break;
			}
			sz++;
		    }
		    offset = i - sz;
		    memmove(s + offset + 1, s + offset, sz);
		    s[offset] = '\'';
		    i++;
		    break;
		}
	    }
	    if (j <= 0 && pc != 0) {
		err = E_UNBAL;
		break;
	    }
	}
    }

    return err;
}
