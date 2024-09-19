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

#define FULL_XML_HEADERS 1

#include "libgretl.h"
#include "uservar.h"
#include "gretl_func.h"
#include "gretl_xml.h"
#include "gretl_typemap.h"
#include "matrix_extra.h"
#include "gretl_cmatrix.h"
#include "gretl_array.h"

/**
 * gretl_array:
 *
 * An opaque type; use the relevant accessor functions.
 */

struct gretl_array_ {
    GretlType type;  /* type of data */
    int n;           /* number of elements */
    void **data;     /* actual data array */
    double *mdata;   /* for matrix block */
};

static void gretl_array_destroy_data (gretl_array *A)
{
    int i;

    if (A->data != NULL) {
	if (A->type == GRETL_TYPE_STRINGS ||
	    A->type == GRETL_TYPE_LISTS) {
	    /* a simple "free" will do */
	    for (i=0; i<A->n; i++) {
		free(A->data[i]);
	    }
	} else if (A->type == GRETL_TYPE_MATRICES) {
	    if (A->mdata != NULL) {
		free(A->mdata);
		for (i=0; i<A->n; i++) {
		    free(A->data[i]);
		}
	    } else {
		for (i=0; i<A->n; i++) {
		    gretl_matrix_free(A->data[i]);
		}
	    }
	} else if (A->type == GRETL_TYPE_BUNDLES) {
	    for (i=0; i<A->n; i++) {
		gretl_bundle_destroy(A->data[i]);
	    }
	} else if (A->type == GRETL_TYPE_ARRAYS) {
	    for (i=0; i<A->n; i++) {
		gretl_array_destroy(A->data[i]);
	    }
	}
	free(A->data);
	A->data = NULL;
    }
}

/* Destroy the whole array, including freeing its
   contents.
*/

void gretl_array_destroy (gretl_array *A)
{
    if (A != NULL) {
	gretl_array_destroy_data(A);
	free(A);
    }
}

/* Reduce the array to empty status, freeing the
   contents but not the array structure itself.
*/

void gretl_array_void_content (gretl_array *A)
{
    if (A != NULL) {
	gretl_array_destroy_data(A);
	A->n = 0;
    }
}

/* Reduce the array to empty status, setting the
   entire data array to NULL without freeing anything.
   Makes sense only if the data array is actually
   "owned" elsewhere.
*/

void gretl_array_nullify_content (gretl_array *A)
{
    if (A != NULL) {
	A->data = NULL;
	A->n = 0;
    }
}

/* Set the array's element pointers to NULL, without
   freeing them. Makes sense only if the element
   pointers are actually "owned" by something else
   (have not been copied into the array).
*/

void gretl_array_nullify_elements (gretl_array *A)
{
    if (A != NULL && A->data != NULL) {
	int i;

	for (i=0; i<A->n; i++) {
	    A->data[i] = NULL;
	}
    }
}

static int array_allocate_storage (gretl_array *A)
{
    int i, err = 0;

    A->data = malloc(A->n * sizeof *A->data);

    if (A->data == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<A->n; i++) {
	    A->data[i] = NULL;
	}
    }

    return err;
}

static int array_extend_content (gretl_array *A, int plus)
{
    if (plus == 0) {
	return 0; /* no-op */
    } else if (plus < 0) {
	return E_DATA;
    } else {
	void **data;
	int n = A->n + plus;
	int i, err = 0;

	data = realloc(A->data, n * sizeof *A->data);

	if (data == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=A->n; i<n; i++) {
		data[i] = NULL;
	    }
	    A->data = data;
	    A->n = n;
	}

	return err;
    }
}

/* Create a new array of type @type. The idea is that it's OK
   to have n = 0: in that case the array starts empty but can
   be extended by the += operator. Or it can be sized in
   advance.
*/

gretl_array *gretl_array_new (GretlType type, int n, int *err)
{
    gretl_array *A;

    if (type != GRETL_TYPE_STRINGS &&
	type != GRETL_TYPE_MATRICES &&
	type != GRETL_TYPE_BUNDLES &&
	type != GRETL_TYPE_LISTS &&
	type != GRETL_TYPE_ARRAYS &&
	type != GRETL_TYPE_ANY) {
	*err = E_TYPES;
	return NULL;
    } else if (n < 0) {
	*err = E_DATA;
	return NULL;
    }

    A = malloc(sizeof *A);

    if (A == NULL) {
	*err = E_ALLOC;
    } else {
	A->type = type;
	A->n = n;
	A->data = NULL;
	A->mdata = NULL;
	if (n > 0) {
	    *err = array_allocate_storage(A);
	    if (*err) {
		gretl_array_destroy(A);
		A = NULL;
	    }
	}
    }

    return A;
}

gretl_array *gretl_singleton_array (void *ptr, GretlType atype,
				    int copy, int *err)
{
    gretl_array *A = gretl_array_new(atype, 1, err);

    if (A != NULL) {
	GretlType t = gretl_type_get_singular(atype);

	*err = gretl_array_set_element(A, 0, ptr, t, copy);
	if (*err) {
	    free(A);
	    A = NULL;
	}
    }

    return A;
}

gretl_array *gretl_array_from_strings (char **S, int n,
				       int copy, int *err)
{
    gretl_array *A;

    A = gretl_array_new(GRETL_TYPE_STRINGS, 0, err);

    if (A != NULL && n > 0) {
	if (copy) {
	    A->data = (void **) strings_array_dup(S, n);
	    if (A->data == NULL) {
		*err = E_ALLOC;
	    }
	} else {
	    A->data = (void **) S;
	}
	if (!*err) {
	    A->n = n;
	}
    }

    return A;
}

#define COMMON_BLOCK 1

gretl_array *gretl_matrix_array_sized (int n, int r, int c,
				       int *err)
{
    gretl_array *A;

    A = gretl_array_new(GRETL_TYPE_MATRICES, n, err);

#if COMMON_BLOCK
    if (A != NULL && n > 0) {
	size_t msize = n * r * c * sizeof(double);
	double *ai_val;
	gretl_matrix *ai;
	int i, rc = r * c;

	/* common block for matrix data */
	ai_val = A->mdata = malloc(msize);
	if (A->mdata == NULL) {
	    *err = E_ALLOC;
	}

	for (i=0; i<n && !*err; i++) {
	    ai = gretl_null_matrix_new();
	    if (ai == NULL) {
		*err = E_ALLOC;
		break;
	    }
	    ai->val = ai_val;
	    ai->rows = r;
	    ai->cols = c;
	    A->data[i] = ai;
	    ai_val += rc;
	}
    }
#else
    if (A != NULL && n > 0) {
	int i;

	for (i=0; i<n && !*err; i++) {
	    A->data[i] = gretl_matrix_alloc(r, c);
	    if (A->data[i] == NULL) {
		*err = E_ALLOC;
	    }
	}
    }
#endif

    if (*err && A != NULL) {
	gretl_array_destroy(A);
	A = NULL;
    }

    return A;
}

/* When we're returning an array of strings, ensure
   that any NULL elements are converted to empty
   strings.
*/

static int strings_array_null_check (gretl_array *A)
{
    int i;

    for (i=0; i<A->n; i++) {
	if (A->data[i] == NULL) {
	    A->data[i] = gretl_strdup("");
	    if (A->data[i] == NULL) {
		return E_ALLOC;
	    }
	}
    }

    return 0;
}

/* note: don't modify the returned value */

char **gretl_array_get_strings (gretl_array *A, int *ns)
{
    char **AS = NULL;

    *ns = 0;

    if (A != NULL && A->type == GRETL_TYPE_STRINGS) {
	int err = strings_array_null_check(A);

	if (!err) {
	    *ns = A->n;
	    AS = (char **) A->data;
	}
    }

    return AS;
}

/* note: take ownership of the returned value */

char **gretl_array_steal_strings (gretl_array *A, int *ns)
{
    char **AS = NULL;

    *ns = 0;

    if (A != NULL && A->type == GRETL_TYPE_STRINGS) {
	int err = strings_array_null_check(A);

	if (!err) {
	    *ns = A->n;
	    AS = (char **) A->data;
	    A->n = 0;
	    A->data = NULL;
	}
    }

    return AS;
}

char **gretl_array_get_stringify_strings (gretl_array *A,
					  int nreq, int *pns,
					  int *err)
{
    char **S = NULL;

    *pns = 0;

    if (A == NULL) {
	*err = E_DATA;
    } else if (A->type != GRETL_TYPE_STRINGS) {
	*err = E_TYPES;
    } else if (A->n < nreq) {
	gretl_errmsg_sprintf("Too few strings: %d given but %d needed",
			     A->n, nreq);
	*err = E_DATA;
    }

    if (!*err) {
	char **AS = (char **) A->data;

	S = strings_array_new(A->n);
	if (S == NULL) {
	    *err = E_ALLOC;
	} else {
	    int myerr = 0;
	    int ndone = 0;
	    int i, j;

	    for (i=0; i<A->n && !myerr; i++) {
		if (AS[i] == NULL || AS[i][0] == '\0') {
		    myerr = E_DATA;
		} else {
		    S[i] = gretl_strdup(AS[i]);
		    if (i > 0) {
			for (j=0; j<i; j++) {
			    if (!strcmp(AS[j], S[i])) {
				gretl_errmsg_sprintf("Duplicated string '%s'", S[i]);
				myerr = E_DATA;
			    }
			}
		    }
		}
		if (!myerr) {
		    ndone++;
		}
	    }
	    if (myerr && ndone < nreq) {
		*err = myerr;
		strings_array_free(S, A->n);
		S = NULL;
	    } else {
		*pns = ndone;
	    }
	}
    }

    return S;
}

/* note: the return value is newly allocated, and owned by the caller */

char *gretl_strings_array_flatten (gretl_array *A,
                                   const char *sep,
                                   int *err)
{
    char *s = NULL;

    if (A == NULL) {
	*err = E_DATA;
    } else if (A->type != GRETL_TYPE_STRINGS) {
	*err = E_TYPES;
    } else {
        int ns = strlen(sep);
	int i, len = 1; /* for terminating NUL */

	for (i=0; i<A->n; i++) {
	    if (A->data[i] == NULL) {
		*err = E_MISSDATA;
		break;
	    } else {
		len += strlen(A->data[i]);
                if (i < A->n - 1) {
                    len += ns;
                }
	    }
	}

	if (!*err) {
	    s = calloc(len, 1);
	    if (s == NULL) {
		*err = E_ALLOC;
	    } else {
		for (i=0; i<A->n; i++) {
		    strcat(s, A->data[i]);
		    if (i < A->n - 1) {
			strcat(s, sep);
		    }
		}
	    }
	}
    }

    return s;
}

/* Return a column vector holding the position(s) in
   the strings array @A at which the string @s is
   matched -- or an empty matrix in case of no matches.
*/

gretl_matrix *gretl_strings_array_pos (gretl_array *A,
				       const char *s,
				       int *err)
{
    gretl_matrix *ret = NULL;
    const char *si;
    int i, np = 0;

    for (i=0; i<A->n; i++) {
	si = A->data[i] == NULL ? "" : A->data[i];
	if (strcmp(si, s) == 0) {
	    np++;
	}
    }

    if (np == 0) {
	ret = gretl_null_matrix_new();
    } else {
	ret = gretl_matrix_alloc(np, 1);
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    } else if (np > 0) {
	int j = 0;

	for (i=0; i<A->n; i++) {
	    si = A->data[i] == NULL ? "" : A->data[i];
	    if (strcmp(si, s) == 0) {
		ret->val[j++] = i+1;
	    }
	}
    }

    return ret;
}

/* Return 1 if @A includes @s, otherwise 0 */

int gretl_strings_array_includes (gretl_array *A, const char *s)
{
    int ret = 0;
    const char *si;
    int i;

    for (i=0; i<A->n; i++) {
	si = A->data[i] == NULL ? "" : A->data[i];
	if (strcmp(si, s) == 0) {
	    ret = 1;
	    break;
	}
    }

    return ret;
}

void *gretl_array_get_element (gretl_array *A, int i,
			       GretlType *type,
			       int *err)
{
    void *ret = NULL;

    /* Note that the data returned here are not "deep copied",
       we just pass the pointer. It's up to the caller to
       decide if a copy has to be made, given that the
       pointer from here should not be modified.
    */

    if (A == NULL) {
	*err = E_UNKVAR;
    } else if (i < 0 || i >= A->n) {
	*err = E_INVARG;
    } else if (A->type == GRETL_TYPE_ANY) {
	*err = E_TYPES;
    } else {
	if (type != NULL) {
	    *type = gretl_type_get_singular(A->type);
	}
	if (A->type == GRETL_TYPE_STRINGS) {
	    if (A->data[i] == NULL) {
		A->data[i] = gretl_strdup("");
	    }
	} else if (A->type == GRETL_TYPE_MATRICES) {
	    if (A->data[i] == NULL) {
		A->data[i] = gretl_null_matrix_new();
	    }
	} else if (A->type == GRETL_TYPE_BUNDLES) {
	    if (A->data[i] == NULL) {
		A->data[i] = gretl_bundle_new();
	    }
	} else if (A->type == GRETL_TYPE_ARRAYS) {
	    if (A->data[i] == NULL) {
		A->data[i] = gretl_array_new(GRETL_TYPE_ANY, 0, err);
	    }
	} else if (A->type == GRETL_TYPE_LISTS) {
	    if (A->data[i] == NULL) {
		A->data[i] = gretl_list_new(0);
	    }
	}
	ret = A->data[i];
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

/* Note: no type-checking, the caller is supposed to
   be on the ball. Use with caution.
*/

void *gretl_array_get_data (gretl_array *A, int i)
{
    if (A == NULL || i < 0 || i >= A->n) {
	return NULL;
    } else {
	return A->data[i];
    }
}

int gretl_array_set_data (gretl_array *A, int i, void *ptr)
{
    if (A == NULL || i < 0 || i >= A->n) {
	return E_DATA;
    } else {
	A->data[i] = ptr;
	return 0;
    }
}

void *gretl_array_get_all_data (gretl_array *A)
{
    if (A == NULL) {
	return NULL;
    } else {
	return A->data;
    }
}

static int check_list_bounds (int *list, int arrdim)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] <= 0 || list[i] > arrdim) {
	    return 0;
	}
    }

    return 1;
}

gretl_array *gretl_array_copy_subspec (gretl_array *A,
				       int *list,
				       int *err)
{
    gretl_array *C = NULL;

    if (A == NULL) {
	*err = E_DATA;
    } else if (!check_list_bounds(list, A->n)) {
	*err = E_INVARG;
    } else {
	int i, k, m = list[0];

	C = gretl_array_new(A->type, m, err);

	for (i=0; i<m && !*err; i++) {
	    k = list[i+1] - 1;
	    if (A->data[k] != NULL) {
		if (A->type == GRETL_TYPE_STRINGS) {
		    C->data[i] = gretl_strdup(A->data[k]);
		} else if (A->type == GRETL_TYPE_MATRICES) {
		    C->data[i] = gretl_matrix_copy(A->data[k]);
		} else if (A->type == GRETL_TYPE_BUNDLES) {
		    C->data[i] = gretl_bundle_copy(A->data[k], err);
		} else if (A->type == GRETL_TYPE_ARRAYS) {
		    C->data[i] = gretl_array_copy(A->data[k], err);
		} else if (A->type == GRETL_TYPE_LISTS) {
		    C->data[i] = gretl_list_copy(A->data[k]);
		}
		if (!*err && C->data[i] == NULL) {
		    *err = E_ALLOC;
		}
	    }
	}
    }

    if (*err && C != NULL) {
	gretl_array_destroy(C);
	C = NULL;
    }

    return C;
}

gretl_bundle *gretl_array_get_bundle (gretl_array *A, int i)
{
    gretl_bundle *b = NULL;

    /* The bundle returned here is not "deep copied",
       we just pass the pointer. It's up to the caller to
       decide if a copy has to be made, given that the
       pointer from here should not be modified.
    */

    if (A != NULL && i >= 0 && i < A->n &&
	A->type == GRETL_TYPE_BUNDLES) {
	b = A->data[i];
    }

    return b;
}

/* "flatten" an array of matrices, yielding a single matrix.
   Concatenation may be horizontal (@vcat = 0, default), vertical
   (@vcat = 1) or (for lack of a better term) "vec-wise" (@vcat =
   2). An error is flagged if the matrices are not conformable for the
   operation.
*/

gretl_matrix *gretl_matrix_array_flatten (gretl_array *A,
					  int vcat,
					  int *err)
{
    gretl_matrix *ret = NULL;
    gretl_matrix *m;
    int common_r = 1;
    int common_c = 1;
    int sum_r = 0;
    int sum_c = 0;
    int r0 = 0, c0 = 0;
    int cmplx = 0;
    int real = 0;
    int i, tot_elem;

    if (A->type != GRETL_TYPE_MATRICES) {
	*err = E_TYPES;
	return NULL;
    }

    for (i=0; i<A->n; i++) {
	m = A->data[i];
	if (!gretl_is_null_matrix(m)) {
	    if (m->is_complex) {
		cmplx = 1;
	    } else {
		real = 1;
	    }
	}
    }

    for (i=0; i<A->n; i++) {
	m = A->data[i];
	if (!gretl_is_null_matrix(m)) {
	    if (c0 == 0) {
		r0 = m->rows;
		c0 = m->cols;
	    } else {
		if (m->rows != r0) {
		    common_r = 0;
		}
		if (m->cols != c0) {
		    common_c = 0;
		}
		if (!common_r && !common_c) {
		    *err = E_NONCONF;
		    break;
		}
	    }
	    sum_r += m->rows;
	    sum_c += m->cols;
	}
    }

    tot_elem = m->rows * m->cols;

    if (!*err) {
	if ((vcat && !common_c) || (!vcat && !common_r)) {
	    *err = E_NONCONF;
	}
    }

    if (!*err && sum_r == 0) {
	ret = gretl_null_matrix_new();
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err || ret != NULL) {
	return ret;
    }

    if (vcat == 0) {
	if (cmplx) {
	    ret = gretl_cmatrix_new(r0, sum_c);
	} else {
	    ret = gretl_matrix_alloc(r0, sum_c);
	}
    } else if (vcat == 1) {
	if (cmplx) {
	    ret = gretl_cmatrix_new(sum_r, c0);
	} else {
	    ret = gretl_matrix_alloc(sum_r, c0);
	}
    } else if (vcat == 2) {
	if (cmplx) {
	    ret = gretl_cmatrix_new(tot_elem, A->n);
	} else {
	    ret = gretl_matrix_alloc(tot_elem, A->n);
	}
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    } else if (vcat == 0) {
	if (cmplx && real) {
	    /* horizontal concatenation: mixed matrices */
	    double complex *dest = ret->z;
	    int n, k = 0;

	    for (i=0; i<A->n; i++) {
		m = A->data[i];
		if (!gretl_is_null_matrix(m)) {
		    if (m->is_complex) {
			memcpy(dest, m->z, tot_elem * sizeof *dest);
		    } else {
			real_to_complex_fill(ret, m, 0, k);
		    }
		    dest += tot_elem;
		    k += m->cols;
		}
	    }
	} else {
	    /* horizontal concatenation: the easy case */
	    double *dest = ret->val;
	    int n, p = cmplx ? 2 : 1;

	    for (i=0; i<A->n; i++) {
		m = A->data[i];
		if (!gretl_is_null_matrix(m)) {
		    n = p * tot_elem;
		    memcpy(dest, m->val, n * sizeof *dest);
		    dest += n;
		}
	    }
	}
    } else if (vcat == 1) {
	/* vertical concatenation */
	int ii, j, k, rpos = 0;
	double complex z;
	double x;

	for (k=0; k<A->n; k++) {
	    m = A->data[k];
	    if (!gretl_is_null_matrix(m)) {
		for (j=0; j<c0; j++) {
		    for (i=0, ii=rpos; i<m->rows; i++, ii++) {
			if (cmplx && !m->is_complex) {
			    z = gretl_matrix_get(m, i, j);
			    gretl_cmatrix_set(ret, ii, j, z);
			} else if (cmplx) {
			    z = gretl_cmatrix_get(m, i, j);
			    gretl_cmatrix_set(ret, ii, j, z);
			} else {
			    x = gretl_matrix_get(m, i, j);
			    gretl_matrix_set(ret, ii, j, x);
			}
		    }
		}
		rpos += m->rows;
	    }
	}
    } else if (vcat == 2) {
	/* vec-wise */
	if (cmplx && real) {
	    /* vec-wise concatenation: mixed matrices */
	    double complex *dest = ret->z;
	    int k = 0;
	    double complex z;

	    for (i=0; i<A->n; i++) {
		m = A->data[i];
		if (!gretl_is_null_matrix(m)) {
		    if (m->is_complex) {
			memcpy(dest, m->z, tot_elem * sizeof *dest);
		    } else {
			for (k=0; k<tot_elem; k++) {
			    z = m->val[k];
			    gretl_cmatrix_set(ret, k, i, z);
			}
		    }
		    dest += tot_elem;
		}
	    }
	} else {
	    /* the easy case */
	    double *dest = ret->val;
	    int n, p = cmplx ? 2 : 1;

	    for (i=0; i<A->n; i++) {
		m = A->data[i];
		if (!gretl_is_null_matrix(m)) {
		    n = p * tot_elem;
		    memcpy(dest, m->val, n * sizeof *dest);
		    dest += n;
		}
	    }
	}
    }

    return ret;
}

static int split_matrix_by_chunks (gretl_array *A, int n,
				   const gretl_matrix *X,
				   int chunk, int colwise)
{
    int rows = colwise ? X->rows : chunk;
    int cols = colwise ? chunk : X->cols;
    gretl_matrix *tmp;
    int i, err = 0;

    /* allocate the matrices */
    for (i=0; i<n && !err; i++) {
	tmp = gretl_matrix_alloc(rows, cols);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_array_set_element(A, i, tmp, GRETL_TYPE_MATRIX, 0);
	}
    }

    /* fill the matrices */
    if (!err && colwise) {
	int nelem = rows * cols;
	size_t csize = nelem * sizeof *X->val;
	const double *src = X->val;

	for (i=0; i<n; i++) {
	    tmp = A->data[i];
	    memcpy(tmp->val, src, csize);
	    src += nelem;
	}
    } else if (!err) {
	int k, j, ii = 0;
	double x;

	for (k=0; k<n; k++) {
	    tmp = A->data[k];
	    for (j=0; j<cols; j++) {
		for (i=0; i<rows; i++) {
		    x = gretl_matrix_get(X, ii+i, j);
		    gretl_matrix_set(tmp, i, j, x);
		}
	    }
	    ii += rows;
	}
    }

    return err;
}

static int split_matrix_by_vector (gretl_array *A, int n,
				   const gretl_matrix *X,
				   int dim, int colwise,
				   gretl_matrix *vals,
				   const double *sel)
{
    int *nj = NULL;
    gretl_matrix *tmp;
    int rows, cols;
    int i, j, k;
    int err = 0;

    /* How many rows or columns will each matrix need? */
    nj = calloc(n, sizeof *nj);
    if (nj == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<dim; i++) {
	    j = sel[i] - 1;
	    nj[j] += 1;
	}
    }
    for (i=0; i<n && !err; i++) {
	if (nj[i] > 0) {
	    rows = colwise ? X->rows : nj[i];
	    cols = colwise ? nj[i] : X->cols;
	    tmp = gretl_matrix_alloc(rows, cols);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		gretl_array_set_element(A, i, tmp, GRETL_TYPE_MATRIX, 0);
		/* @nj will serve as an array of counters below */
		nj[i] = 0;
	    }
	}
    }

    /* fill the individual matrices */
    if (!err && colwise) {
	size_t csize = X->rows * sizeof *X->val;
	const double *src = X->val;
	double *targ;

	for (i=0; i<X->cols; i++) {
	    k = sel[i] - 1;
	    tmp = A->data[k];
	    targ = tmp->val + nj[k];
	    memcpy(targ, src, csize);
	    /* advance the read and write positions */
	    src += X->rows;
	    nj[k] += X->rows;
	}
    } else if (!err) {
	double x;
	int jj;

	for (i=0; i<X->rows; i++) {
	    k = sel[i] - 1;
	    tmp = A->data[k];
	    jj = 0;
	    for (j=0; j<X->cols; j++) {
		x = gretl_matrix_get(X, i, j);
		gretl_matrix_set(tmp, nj[k], jj++, x);
	    }
	    /* advance the write position */
	    nj[k] += 1;
	}
    }

    free(nj);

    return err;
}

gretl_array *gretl_matrix_split_by (const gretl_matrix *X,
				    const gretl_matrix *v,
				    int colwise, int chunks,
                                    int *err)
{
    gretl_array *ret = NULL;
    int i, dim, nm = 0;
    int chunksize = 0;
    gretl_matrix *vals = NULL;
    const double *sel;
    double x;

    dim = colwise ? X->cols : X->rows;

    if (chunks) {
	/* interpret single element of @v as chunk size */
	chunksize = v->val[0];
	if (chunksize <= 0 || chunksize > dim) {
	    *err = E_INVARG;
	} else if (dim % chunksize != 0) {
	    *err = E_NONCONF;
	}
    } else if (gretl_vector_get_length(v) == dim) {
	/* vector of indices: only positive integers allowed */
	for (i=0; i<dim; i++) {
	    x = v->val[i];
	    if (x != floor(x) || x <= 0 || x >= INT_MAX) {
		*err = E_INVARG;
		break;
	    }
	}
    } else {
	*err = E_INVARG;
    }

    if (*err) {
	return NULL;
    }

    /* How many matrices do we need ? */
    if (chunksize > 0) {
	nm = dim / chunksize;
    } else {
	sel = v->val;
	vals = gretl_matrix_values(sel, dim, OPT_NONE, err);
	if (!*err) {
	    /* get the maximum index value */
	    for (i=0; i<vals->rows; i++) {
		nm = vals->val[i] > nm ? vals->val[i] : nm;
	    }
	}
    }

    if (!*err) {
	ret = gretl_array_new(GRETL_TYPE_MATRICES, nm, err);
    }

    if (!*err && chunksize > 0) {
	/* the easier case */
	*err = split_matrix_by_chunks(ret, nm, X, chunksize, colwise);
    } else if (!*err) {
	/* the fiddly but general one */
	*err = split_matrix_by_vector(ret, nm, X, dim, colwise, vals, sel);
    }

    gretl_matrix_free(vals);

    if (*err && ret != NULL) {
	gretl_array_destroy(ret);
	ret = NULL;
    }

    return ret;
}

int gretl_array_set_type (gretl_array *A, GretlType type)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (type != GRETL_TYPE_STRINGS &&
	       type != GRETL_TYPE_MATRICES &&
	       type != GRETL_TYPE_BUNDLES &&
	       type != GRETL_TYPE_LISTS &&
	       type != GRETL_TYPE_ARRAYS) {
	err = E_TYPES;
    } else if (type == A->type) {
	/* no-op */
	return 0;
    } else if (A->n > 0) {
	/* we can (re-)set the type only if no data have
	   been entered already
	*/
	int i;

	for (i=0; i<A->n; i++) {
	    if (A->data[i] != NULL) {
		err = E_DATA;
		break;
	    }
	}
    }

    if (!err) {
	A->type = type;
    }

    return err;
}

GretlType gretl_array_get_type (gretl_array *A)
{
    return (A != NULL)? A->type : GRETL_TYPE_NONE;
}

GretlType gretl_array_get_content_type (gretl_array *A)
{
    return (A != NULL)? gretl_type_get_singular(A->type) : GRETL_TYPE_NONE;
}

int gretl_array_get_length (gretl_array *A)
{
    return (A != NULL)? A->n : 0;
}

/* Return the 0-based index of the first empty slot
   in @A, or -1 on failure (if @A is NULL or all its
   slots are already filled).
*/

int gretl_array_get_next_index (gretl_array *A)
{
    int ret = -1;

    if (A != NULL) {
	int i;

	for (i=0; i<A->n; i++) {
	    if (A->data[i] == NULL) {
		ret = i;
		break;
	    }
	}
    }

    return ret;
}

static int set_string (gretl_array *A, int i,
		       char *s, int copy)
{
    int err = 0;

    if (copy) {
	A->data[i] = gretl_strdup(s);
	if (A->data[i] == NULL) {
	    err = E_ALLOC;
	}
    } else {
	A->data[i] = s;
    }

    return err;
}

/* The static set_*() functions assume that error-checking
   has already been done.
*/

static int set_matrix (gretl_array *A, int i,
		       gretl_matrix *m, int copy)
{
    int err = 0;

    if (copy) {
	A->data[i] = gretl_matrix_copy(m);
	if (A->data[i] == NULL) {
	    err = E_ALLOC;
	}
    } else {
	A->data[i] = m;
    }

    return err;
}

static int set_bundle (gretl_array *A, int i,
		       gretl_bundle *b, int copy)
{
    int err = 0;

    if (copy) {
	A->data[i] = gretl_bundle_copy(b, &err);
    } else {
	A->data[i] = b;
    }

    return err;
}

static int set_array (gretl_array *A, int i,
		      gretl_array *a, int copy)
{
    int err = 0;

    if (copy) {
	A->data[i] = gretl_array_copy(a, &err);
    } else {
	A->data[i] = a;
    }

    return err;
}

static int set_list (gretl_array *A, int i,
		     int *L, int copy)
{
    int err = 0;

    if (copy) {
	A->data[i] = gretl_list_copy(L);
	if (A->data[i] == NULL) {
	    err = E_ALLOC;
	}
    } else {
	A->data[i] = L;
    }

    return err;
}

static int set_type_error (gretl_array *A,
			   GretlType required)
{
    if (A->type == GRETL_TYPE_ANY) {
	A->type = required;
	return 0;
    } else if (A->type != required) {
	GretlType reqt = gretl_type_get_singular(required);

	gretl_errmsg_sprintf("Cannot add %s to array of %s",
			     gretl_type_get_name(reqt),
			     gretl_type_get_name(A->type));
	return 1;
    } else {
	return 0;
    }
}

/* In the functions below we assume the @copy parameter
   will be set appropriately by "genr": if the incoming
   object is a named variable in its own right it should
   be copied, but if it's on on-the-fly thing there's
   no need to copy it, just "donate" it to the array.
*/

/* respond to A[i] = s */

int gretl_array_set_string (gretl_array *A, int i,
			    char *s, int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_STRINGS)) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else if (s != A->data[i]) {
	free(A->data[i]);
	err = set_string(A, i, s, copy);
    }

    return err;
}

/* respond to A += s */

int gretl_array_append_string (gretl_array *A,
			       char *s,
			       int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_STRINGS)) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A, 1);
	if (!err) {
	    err = set_string(A, A->n - 1, s, copy);
	}
    }

    return err;
}

/* respond to A[i] = m */

int gretl_array_set_matrix (gretl_array *A, int i,
			    gretl_matrix *m,
			    int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_MATRICES)) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else if (m != A->data[i]) {
	gretl_matrix_free(A->data[i]);
	err = set_matrix(A, i, m, copy);
    }

    return err;
}

/* respond to A += m */

int gretl_array_append_matrix (gretl_array *A,
			       gretl_matrix *m,
			       int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_MATRICES)) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A, 1);
	if (!err) {
	    err = set_matrix(A, A->n - 1, m, copy);
	}
    }

    return err;
}

/* respond to A[i] = b */

int gretl_array_set_bundle (gretl_array *A, int i,
			    gretl_bundle *b,
			    int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_BUNDLES)) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else if (b != A->data[i]) {
	gretl_bundle_destroy(A->data[i]);
	err = set_bundle(A, i, b, copy);
    }

    return err;
}

/* respond to A[i] = a */

int gretl_array_set_array (gretl_array *A, int i,
			   gretl_array *a,
			   int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_ARRAYS)) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else if (a != A->data[i]) {
	gretl_array_destroy(A->data[i]);
	err = set_array(A, i, a, copy);
    }

    return err;
}

/* respond to A += b */

int gretl_array_append_bundle (gretl_array *A,
			       gretl_bundle *b,
			       int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_BUNDLES)) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A, 1);
	if (!err) {
	    err = set_bundle(A, A->n - 1, b, copy);
	}
    }

    return err;
}

/* respond to A += a */

int gretl_array_append_array (gretl_array *A,
			      gretl_array *a,
			      int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_ARRAYS)) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A, 1);
	if (!err) {
	    err = set_array(A, A->n - 1, a, copy);
	}
    }

    return err;
}

/* respond to A[i] = L */

int gretl_array_set_list (gretl_array *A, int i,
			  int *L, int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_LISTS)) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else if (L != A->data[i]) {
	free(A->data[i]);
	err = set_list(A, i, L, copy);
    }

    return err;
}

/* respond to A += L */

int gretl_array_append_list (gretl_array *A,
			     int *L, int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (set_type_error(A, GRETL_TYPE_LISTS)) {
	err = E_TYPES;

    } else {
	err = array_extend_content(A, 1);
	if (!err) {
	    err = set_list(A, A->n - 1, L, copy);
	}
    }

    return err;
}

int gretl_array_set_element (gretl_array *A, int i,
			     void *ptr, GretlType type,
			     int copy)
{
    int err = 0;

    if (type == GRETL_TYPE_MATRIX) {
	err = gretl_array_set_matrix(A, i, ptr, copy);
    } else if (type == GRETL_TYPE_STRING) {
	err = gretl_array_set_string(A, i, ptr, copy);
    } else if (type == GRETL_TYPE_BUNDLE) {
	err = gretl_array_set_bundle(A, i, ptr, copy);
    } else if (type == GRETL_TYPE_LIST) {
	err = gretl_array_set_list(A, i, ptr, copy);
    } else if (type == GRETL_TYPE_ARRAY) {
	err = gretl_array_set_array(A, i, ptr, copy);
    }

    return err;
}

static void free_array_element (gretl_array *A, int i)
{
    if (A->type == GRETL_TYPE_STRINGS ||
	A->type == GRETL_TYPE_LISTS) {
	free(A->data[i]);
    } else if (A->type == GRETL_TYPE_MATRICES) {
	gretl_matrix_free(A->data[i]);
    } else if (A->type == GRETL_TYPE_BUNDLES) {
	gretl_bundle_destroy(A->data[i]);
    } else if (A->type == GRETL_TYPE_ARRAYS) {
	gretl_array_destroy(A->data[i]);
    }
}

int gretl_array_delete_element (gretl_array *A, int i)
{
    if (A == NULL) {
	return E_DATA;
    } else if (i < 0 || i >= A->n) {
	return E_INVARG;
    } else {
	int j;

	if (A->data[i] != NULL) {
	    free_array_element(A, i);
	}
	/* shift the higher-numbered elements down */
	for (j=i; j<A->n-1; j++) {
	    A->data[j] = A->data[j+1];
	}
	/* decrement the element count */
	A->n -= 1;

	return 0;
    }
}

/* Drop any/all occurrences of string @s from array @A */

int gretl_array_drop_string (gretl_array *A, const char *s)
{
    if (A->type != GRETL_TYPE_STRINGS) {
	return E_TYPES;
    } else {
	int i, j, rem = A->n;
	int n_orig = A->n;
	size_t sz;

	for (i=0; ; ) {
	    if (A->data[i] != NULL && !strcmp(A->data[i], s)) {
		free(A->data[i]);
		j = i + 1;
		rem = A->n - j;
		sz = rem * sizeof *A->data;
		memmove(A->data + i, A->data + j, sz);
		A->n -= 1;
	    } else {
		i++;
	    }
	    if (i == A->n || rem == 0) {
		break;
	    }
	}
	if (A->n == 0) {
	    free(A->data);
	    A->data = NULL;
	} else if (A->n < n_orig) {
	    A->data = realloc(A->data, A->n * sizeof *A->data);
	}
    }

    return 0;
}

int gretl_array_drop_null (gretl_array *A)
{
    int i, j, rem = A->n;
    int n_orig = A->n;
    size_t sz;

    for (i=0; ; ) {
	if (A->data[i] == NULL) {
	    j = i + 1;
	    rem = A->n - j;
	    sz = rem * sizeof *A->data;
	    memmove(A->data + i, A->data + j, sz);
	    A->n -= 1;
	} else {
	    i++;
	}
	if (i == A->n || rem == 0) {
	    break;
	}
    }
    if (A->n == 0) {
	free(A->data);
	A->data = NULL;
    } else if (A->n < n_orig) {
	A->data = realloc(A->data, A->n * sizeof *A->data);
    }

    return 0;
}

/* @ptr must be pre-checked as matching the array type */

int gretl_array_append_object (gretl_array *A,
			       void *ptr,
			       int copy)
{
    int err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (A->type == GRETL_TYPE_MATRICES) {
	gretl_array_append_matrix(A, ptr, copy);
    } else if (A->type == GRETL_TYPE_STRINGS) {
	gretl_array_append_string(A, ptr, copy);
    } else if (A->type == GRETL_TYPE_BUNDLES) {
	gretl_array_append_bundle(A, ptr, copy);
    } else if (A->type == GRETL_TYPE_LISTS) {
	gretl_array_append_list(A, ptr, copy);
    } else if (A->type == GRETL_TYPE_ARRAYS) {
	gretl_array_append_array(A, ptr, copy);
    }

    return err;
}

static int
gretl_array_copy_content (gretl_array *Acpy, const gretl_array *A,
			  int write_offset)
{
    int i, j, err = 0;

    for (i=0; i<A->n && !err; i++) {
	if (A->data[i] != NULL) {
	    j = i + write_offset;
	    if (A->type == GRETL_TYPE_STRINGS) {
		Acpy->data[j] = gretl_strdup(A->data[i]);
	    } else if (A->type == GRETL_TYPE_MATRICES) {
		Acpy->data[j] = gretl_matrix_copy(A->data[i]);
	    } else if (A->type == GRETL_TYPE_BUNDLES) {
		Acpy->data[j] = gretl_bundle_copy(A->data[i], &err);
	    } else if (A->type == GRETL_TYPE_ARRAYS) {
		Acpy->data[j] = gretl_array_copy(A->data[i], &err);
	    } else if (A->type == GRETL_TYPE_LISTS) {
		Acpy->data[j] = gretl_list_copy(A->data[i]);
	    } else {
		err = E_TYPES;
	    }
	    if (!err && Acpy->data[j] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

gretl_array *gretl_array_copy (const gretl_array *A,
			       int *err)
{
    gretl_array *Acpy = NULL;

    if (A != NULL) {
	Acpy = gretl_array_new(A->type, A->n, err);
	if (!*err) {
	    *err = gretl_array_copy_content(Acpy, A, 0);
	}
    }

    return Acpy;
}

static int inverse_compare_strings (const void *a, const void *b)
{
    const char **sa = (const char **) a;
    const char **sb = (const char **) b;

    return -g_strcmp0(*sa, *sb);
}

gretl_array *gretl_strings_sort (const gretl_array *A,
				 int descending,
				 int *err)
{
    gretl_array *Asrt = NULL;

    if (A != NULL) {
	if (A->type != GRETL_TYPE_STRINGS) {
	    *err = E_TYPES;
	} else {
	    Asrt = gretl_array_new(A->type, A->n, err);
	}
	if (!*err) {
	    *err = gretl_array_copy_content(Asrt, A, 0);
	}
	if (!*err) {
	    qsort(Asrt->data, Asrt->n, sizeof *Asrt->data,
		  descending ? inverse_compare_strings : gretl_compare_strings);
	}
    }

    return Asrt;
}

/* respond to A1 += A2 */

int gretl_array_copy_into (gretl_array *A1,
			   const gretl_array *A2)
{
    gretl_array *Tmp = NULL;
    int old_n = 0, err = 0;

    if (A1 == NULL || A2 == NULL) {
	err = E_DATA;
    } else if (A1->type != A2->type) {
	err = E_TYPES;
    } else if (A1 == A2) {
        Tmp = gretl_array_copy(A1, &err);
        if (!err) {
            A2 = Tmp;
        }
    }

    if (!err) {
	old_n = A1->n;
	err = array_extend_content(A1, A2->n);
    }
    if (!err) {
	err = gretl_array_copy_content(A1, A2, old_n);
    }
    if (Tmp != NULL) {
        gretl_array_destroy(Tmp);
    }

    return err;
}

/* respond to C = A + B */

gretl_array *gretl_arrays_join (gretl_array *A,
				gretl_array *B,
				int *err)
{
    gretl_array *C = NULL;

    if (A == NULL || B == NULL) {
	*err = E_DATA;
    } else if (A->type != B->type) {
	*err = E_TYPES;
    } else {
	int n = A->n + B->n;

	C = gretl_array_new(A->type, n, err);
    }

    if (!*err) {
	*err = gretl_array_copy_content(C, A, 0);
    }

    if (!*err) {
	*err = gretl_array_copy_content(C, B, A->n);
    }

    if (*err && C != NULL) {
	gretl_array_destroy(C);
	C = NULL;
    }

    return C;
}

/* respond to C = A || B */

gretl_array *gretl_arrays_union (gretl_array *A,
				 gretl_array *B,
				 int *err)
{
    gretl_array *C = NULL;
    const char *sa, *sb;
    char *copy = NULL;
    int i, j, n = 0;

    if (A == NULL || B == NULL) {
	*err = E_DATA;
    } else if (A->type != GRETL_TYPE_STRINGS ||
	       B->type != GRETL_TYPE_STRINGS) {
	*err = E_TYPES;
    } else {
	if (B->n > 0) {
	    copy = calloc(1, B->n);
	}
	n = A->n;
	for (j=0; j<B->n; j++) {
	    sb = B->data[j];
	    if (sb == NULL || *sb == '\0') {
		continue;
	    }
	    copy[j] = 1;
	    for (i=0; i<A->n; i++) {
		sa = A->data[i];
		if (sa == NULL || *sa == '\0') {
		    continue;
		}
		if (strcmp(sa, sb) == 0) {
		    copy[j] = 0;
		    break;
		}
	    }
	    n += copy[j];
	}
	C = gretl_array_new(A->type, n, err);
    }

    if (!*err) {
	*err = gretl_array_copy_content(C, A, 0);
    }

    if (!*err && n > A->n) {
	i = A->n;
	for (j=0; j<B->n; j++) {
	    if (copy[j]) {
		C->data[i++] = gretl_strdup(B->data[j]);
	    }
	}
    }
    free(copy);

    if (*err && C != NULL) {
	gretl_array_destroy(C);
	C = NULL;
    }

    return C;
}

/* respond to C = A && B */

gretl_array *gretl_arrays_intersection (gretl_array *A,
					gretl_array *B,
					int *err)
{
    gretl_array *C = NULL;
    const char *sa, *sb;
    char *copy = NULL;
    int i, j, n = 0;

    if (A == NULL || B == NULL) {
	*err = E_DATA;
    } else if (A->type != GRETL_TYPE_STRINGS ||
	       B->type != GRETL_TYPE_STRINGS) {
	*err = E_TYPES;
    } else {
	if (A->n > 0) {
	    copy = calloc(1, A->n);
	}
	for (i=0; i<A->n; i++) {
	    sa = A->data[i];
	    if (sa == NULL || *sa == '\0') {
		continue;
	    }
	    for (j=0; j<B->n; j++) {
		sb = B->data[j];
		if (sb == NULL || *sb == '\0') {
		    continue;
		}
		if (strcmp(sa, sb) == 0) {
		    copy[i] = 1;
		    break;
		}
	    }
	    n += copy[i];
	}
	C = gretl_array_new(A->type, n, err);
    }

    if (!*err && n > 0) {
	j = 0;
	for (i=0; i<A->n; i++) {
	    if (copy[i]) {
		C->data[j++] = gretl_strdup(A->data[i]);
	    }
	}
    }
    free(copy);

    if (*err && C != NULL) {
	gretl_array_destroy(C);
	C = NULL;
    }

    return C;
}

/* respond to C = A ~ B, for strings only so far */

gretl_array *gretl_arrays_concat (const gretl_array *A,
                                  const gretl_array *B,
                                  int *err)
{
    gretl_array *C = NULL;

    if (A == NULL || A->type != GRETL_TYPE_STRINGS ||
        B == NULL || B->type != GRETL_TYPE_STRINGS ||
        A->n != B->n) {
        *err = E_TYPES;
        return NULL;
    }

    C = gretl_array_new(GRETL_TYPE_STRINGS, A->n, err);

    if (C != NULL) {
        char *si, *sa, *sb;
        int i;

        for (i=0; i<A->n; i++) {
            sa = A->data[i] == NULL ? "" : (char *) A->data[i];
            sb = B->data[i] == NULL ? "" : (char *) B->data[i];
            si = calloc(strlen(sa) + strlen(sb) + 1, 1);
            sprintf(si, "%s%s", sa, sb);
            C->data[i] = si;
        }
    }

    return C;
}

/**
 * gretl_array_copy_as:
 * @name: name of source array.
 * @copyname: name for copy.
 * @copytpe: the type specified for the copied array, or 0.
 *
 * Look for a saved array specified by @name, and if found,
 * make a full copy and save it under @copyname. This is
 * called from geneval.c on completion of assignment to a
 * array named @copyname, where the returned value on the
 * right-hand side is a pre-existing array.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_array_copy_as (const char *name, const char *copyname,
			 GretlType copytype)
{
    gretl_array *A0, *A1 = NULL;
    user_var *u;
    int err = 0;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_ARRAY);
    if (u == NULL) {
	return E_DATA;
    } else {
	A0 = user_var_get_value(u);
    }

    if (copytype > 0 && A0->type != copytype) {
	return E_TYPES;
    }

    /* is there a pre-existing array named @copyname? */
    u = get_user_var_of_type_by_name(copyname, GRETL_TYPE_ARRAY);
    if (u != NULL) {
	A1 = user_var_get_value(u);
    }

    if (A1 != NULL) {
	if (A1->type != A0->type) {
	    err = E_TYPES;
	} else {
	    gretl_array_void_content(A1);
	    A1->n = A0->n;
	    err = array_allocate_storage(A1);
	    if (!err) {
		err = gretl_array_copy_content(A1, A0, 0);
	    }
	}
    } else {
	A1 = gretl_array_copy(A0, &err);
	if (!err) {
	    err = user_var_add(copyname, A1->type, A1);
	}
    }

    return err;
}

/**
 * get_array_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved array, if found, else NULL.
 */

gretl_array *get_array_by_name (const char *name)
{
    gretl_array *a = NULL;

    if (name != NULL && *name != '\0') {
	user_var *u =
	    get_user_var_of_type_by_name(name, GRETL_TYPE_ARRAY);

	if (u != NULL) {
	    a = user_var_get_value(u);
	}
    }

    return a;
}

/**
 * get_strings_array_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved array of strings, if found, else NULL.
 */

gretl_array *get_strings_array_by_name (const char *name)
{
    gretl_array *ret = NULL;

    if (name != NULL && *name != '\0') {
	user_var *u =
	    get_user_var_of_type_by_name(name, GRETL_TYPE_ARRAY);

	if (u != NULL) {
	    ret = user_var_get_value(u);
	    if (ret->type != GRETL_TYPE_STRINGS) {
		ret = NULL;
	    }
	}
    }

    return ret;
}

/**
 * get_strings_array_from_series:
 * @dset: gretl dataset.
 * @v: series ID number.
 * @err: location to receive error code.
 *
 * Returns: a newly allocated array of strings holding the string
 * values of series @v over the current sample range, or NULL on
 * failure.
 */

gretl_array *get_strings_array_from_series (DATASET *dset,
					    int v, int *err)
{
    gretl_array *ret = NULL;
    const char *st;
    int i, t, n;

    if (!is_string_valued(dset, v)) {
	*err = E_TYPES;
	return NULL;
    }

    n = dset->t2 - dset->t1 + 1;
    ret = gretl_array_new(GRETL_TYPE_STRINGS, n, err);
    if (ret == NULL) {
	return NULL;
    }

    for (t=dset->t1, i=0; t<=dset->t2; t++, i++) {
	st = series_get_string_for_obs(dset, v, t);
	ret->data[i] = gretl_strdup(st);
    }

    return ret;
}

gretl_array *gretl_array_pull_from_stack (const char *name,
					  int *err)
{
    gretl_array *a = NULL;
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_ARRAY);

    if (u != NULL) {
	a = user_var_unstack_value(u);
    }

    if (a == NULL) {
	*err = E_DATA;
    }

    return a;
}

static void print_array_string (const char *s, PRN *prn)
{
    int n = strcspn(s, "\r\n");
    int m = strlen(s);

    if (n > 72) {
	pprintf(prn, "\"%.69s...\"\n", s);
    } else if (n < m) {
	pprintf(prn, "\"%.*s...\"\n", n, s);
    } else {
	pprintf(prn, "\"%s\"\n", s);
    }
}

static void print_array_elements (gretl_array *A,
				  int imin, int imax,
				  int range, PRN *prn)
{
    int i, lim = MIN(A->n, imax);

    for (i=imin; i<lim; i++) {
	pprintf(prn, "[%d] ", i+1);
	if (A->data[i] == NULL) {
	    pputs(prn, "null\n");
	} else if (A->type == GRETL_TYPE_STRINGS) {
	    const char *s = A->data[i];

	    print_array_string(s, prn);
	} else if (A->type == GRETL_TYPE_MATRICES) {
	    const gretl_matrix *m = A->data[i];

	    pprintf(prn, "%d x %d\n", m->rows, m->cols);
	} else if (A->type == GRETL_TYPE_LISTS) {
	    const int *list = A->data[i];

	    gretl_list_print(list, NULL, prn);
	}
    }

    if (!range && A->n > lim) {
	pputs(prn, "...\n\n");
    } else {
	pputc(prn, '\n');
    }
}

int gretl_array_print (gretl_array *A, PRN *prn)
{
    if (A != NULL) {
	const char *s = gretl_type_get_name(A->type);
	int nmax = 10;

	pprintf(prn, _("Array of %s, length %d\n"), s, A->n);

	if (A->n > 0 &&
	    A->type != GRETL_TYPE_BUNDLES && A->type != GRETL_TYPE_ARRAYS) {
	    print_array_elements(A, 0, nmax, 0, prn);
	}
    }

    return 0;
}

int gretl_array_print_range (gretl_array *A, int imin, int imax, PRN *prn)
{
    if (A != NULL) {
	const char *s = gretl_type_get_name(A->type);

	pprintf(prn, _("Array of %s, length %d\n"), s, A->n);

	if (A->type != GRETL_TYPE_BUNDLES && A->type != GRETL_TYPE_ARRAYS) {
	    print_array_elements(A, imin, imax, 1, prn);
	}
    }

    return 0;
}

/* Called from gretl_bundle.c when serializing a bundle
   which contains one or more arrays.
*/

void gretl_array_serialize (gretl_array *A, PRN *prn)
{
    GretlType type;
    const char *subname;
    void *ptr;
    int i;

    type = gretl_type_get_singular(A->type);
    subname = gretl_type_get_name(type);

    pprintf(prn, "<gretl-array type=\"%s\" length=\"%d\">\n",
	    gretl_type_get_name(A->type), A->n);

    for (i=0; i<A->n; i++) {
	ptr = A->data[i];
	if (ptr == NULL) {
	    pprintf(prn, "<%s placeholder=\"true\"/>\n", subname);
	} else if (type == GRETL_TYPE_STRING) {
	    gretl_xml_put_tagged_string("string", ptr, prn);
	} else if (type == GRETL_TYPE_MATRIX) {
	    gretl_matrix_serialize(ptr, NULL, prn);
	} else if (type == GRETL_TYPE_BUNDLE) {
	    gretl_bundle_serialize(ptr, NULL, prn);
	} else if (type == GRETL_TYPE_ARRAY) {
	    gretl_array_serialize(ptr, prn);
	} else if (type == GRETL_TYPE_LIST) {
	    gretl_xml_put_tagged_list("list", ptr, prn);
	}
    }

    pputs(prn, "</gretl-array>\n");
}

static int name_matches_array_type (char *s, GretlType type)
{
    if (!strncmp(s, "gretl-", 6)) {
	/* we're expecting "short-form" type strings */
	s += 6;
    }

    return (gretl_type_from_string(s) == type);
}

static int deserialize_array_elements (gretl_array *A,
				       xmlNodePtr cur,
				       xmlDocPtr doc)
{
    GretlType type = gretl_type_get_singular(A->type);
    int i = 0, err = 0;

    while (cur != NULL && !err && i < A->n) {
	if (!name_matches_array_type((char *) cur->name, type)) {
	    fprintf(stderr, "deserialize array: mismatched element '%s'\n",
		    (char *) cur->name);
	    err = E_DATA;
	} else if (gretl_xml_get_prop_as_bool(cur, "placeholder")) {
	    ; /* null array element: no-op */
	} else if (A->type == GRETL_TYPE_STRINGS) {
	    A->data[i] = gretl_xml_get_string(cur, doc);
	} else if (A->type == GRETL_TYPE_MATRICES) {
	    A->data[i] = gretl_xml_get_matrix(cur, doc, &err);
	} else if (A->type == GRETL_TYPE_BUNDLES) {
	    A->data[i] = gretl_bundle_deserialize(cur, doc, &err);
	} else if (A->type == GRETL_TYPE_ARRAYS) {
	    A->data[i] = gretl_array_deserialize(cur, doc, &err);
	} else if (A->type == GRETL_TYPE_LISTS) {
	    A->data[i] = gretl_xml_get_list(cur, doc, &err);
	}
	/* note: arrays of scalars not handled */
	i++;
	cur = cur->next;
    }

    if (!err && i != A->n) {
	fprintf(stderr, "deserialize array: array is corrupted\n");
	err = E_DATA;
    }

    return err;
}

/* For internal use only: @p1 should be of type xmlNodePtr and @p2
   should be an xmlDocPtr. We suppress the actual pointer types in the
   prototype so that it's possible for a module to include
   gretl_array.h without including the full libxml headers.
*/

gretl_array *gretl_array_deserialize (void *p1, void *p2,
				      int *err)
{
    xmlNodePtr node = p1;
    xmlDocPtr doc = p2;
    GretlType type = 0;
    int n = 0;
    gretl_array *A = NULL;

    if (xmlStrcmp(node->name, (XUC) "gretl-array")) {
	fprintf(stderr, "deserialize array: node is not gretl-array!\n");
	*err = E_DATA;
    } else {
	type = gretl_xml_get_type_property(node);
	if (type == 0) {
	    fprintf(stderr, "deserialize array: couldn't get array type\n");
	    *err = E_DATA;
	}
    }

    if (!*err && !gretl_xml_get_prop_as_int(node, "length", &n)) {
	fprintf(stderr, "deserialize array: couldn't get length\n");
	*err = E_DATA;
    }

    if (!*err) {
	A = gretl_array_new(type, n, err);
    }

    if (A != NULL && n > 0) {
	*err = deserialize_array_elements(A, node->xmlChildrenNode, doc);
	if (*err) {
	    gretl_array_destroy(A);
	    A = NULL;
	}
    }

    return A;
}

/* In case a matrix is too wide for comfortable printing, split it by
   column. If there are one or more "leading" columns that should be
   displayed with each chunk of the matrix this should be signalled
   via a non-zero value for @leadcols. The maximum number of columns
   per chunk (including the leading columns, if any) is set by the
   @maxcols argument.

   The return value is an array of suitably sized matrices. If @m
   has column names attached these are distributed to the matrices
   in the array.
*/

gretl_array *gretl_matrix_col_split (const gretl_matrix *m,
				     int leadcols, int maxcols,
				     int *err)
{
    gretl_array *a = NULL;
    int maincols, nm, rem;

    if (gretl_is_null_matrix(m)) {
	*err = E_INVARG;
    } else {
	maincols = m->cols - leadcols;
	if (maincols <= 0) {
	    *err = E_INVARG;
	}
    }

    if (!*err) {
	/* how many matrices will we need? */
	nm = maincols / (maxcols - leadcols);
	rem = maincols % (maxcols - leadcols);
	nm += rem > 0;
	if (nm == 1) {
	    /* nothing to be done */
	    *err = E_INVARG;
	}
    }

    if (!*err) {
	a = gretl_array_new(GRETL_TYPE_MATRICES, nm, err);
    }

    if (!*err) {
	const char **S0 = gretl_matrix_get_colnames(m);
	char **Si = NULL;
	size_t rsize = m->rows * sizeof(double);
	const double *src;
	double *targ;
	gretl_matrix *ai;
	int i, j, cols, spos;

	/* initial read position for "main" columns */
	src = m->val + leadcols * m->rows;
	spos = leadcols;

	for (i=0; i<nm && !*err; i++) {
	    cols = (i == nm-1 && rem > 0)? rem + leadcols : maxcols;
	    ai = gretl_zero_matrix_new(m->rows, cols);
	    if (ai == NULL) {
		*err = E_ALLOC;
	    } else {
		Si = S0 == NULL ? NULL : strings_array_new(cols);
		/* initial write position */
		targ = ai->val;
		if (leadcols > 0) {
		    memcpy(targ, m->val, leadcols * rsize);
		    targ += leadcols * m->rows;
		    if (Si != NULL) {
			/* transcribe column names */
			for (j=0; j<leadcols; j++) {
			    Si[j] = gretl_strdup(S0[j]);
			}
		    }
		}
		/* now handle the "main" columns */
		cols -= leadcols;
		memcpy(targ, src, cols * rsize);
		/* advance the read position */
		src += cols * m->rows;
		if (Si != NULL) {
		    for (j=0; j<cols; j++) {
			Si[leadcols+j] = gretl_strdup(S0[spos++]);
		    }
		    gretl_matrix_set_colnames(ai, Si);
		}
		/* stick matrix @i into the array */
		gretl_array_set_matrix(a, i, ai, 0);
	    }
	}
    }

    if (*err && a != NULL) {
	gretl_array_destroy(a);
	a = NULL;
    }

    return a;
}

/* is_strings_array_element() tests @str for representation of an
   element of an array of strings. This should work if @str takes the
   form

   <arrayname>[<index>]

   where @arrayname identifies a strings array and <index> represents
   a valid (1-based) index into this array. In that case @arrayname is
   written into @aname (which should be of length 32 bytes or more),
   @index is written into @idx (should be 8 bytes or more), and 1 is
   returned. On failure, 0 is returned.
*/

int is_strings_array_element (const char *str,
			      char *aname,
			      char *idx)
{
    int ret = 0;

    if (strchr(str, '[') != NULL) {
	gretl_array *A = NULL;
	int i, err = 0;

	if (sscanf(str, "%31[^[][%7[^]]", aname, idx) == 2) {
	    A = get_strings_array_by_name(aname);
	}
	if (A != NULL) {
	    i = generate_int(idx, NULL, &err);
	    if (!err && i > 0 && i <= A->n) {
		ret = 1;
	    }
	}
    }

    return ret;
}

/* start apparatus for sorting a gretl array via qsort */

#define AS_DEBUG 0

struct asorter {
    fncall *call;
    GretlType type;
    DATASET *dset;
    PRN *prn;
    int err;
};

static struct asorter asort;

static void asort_cleanup (void)
{
    fncall_destroy(asort.call);
    asort.call = NULL;
    asort.type = 0;
    asort.dset = NULL;
    asort.prn = NULL;
    asort.err = 0;
}

static int asort_setup (const char *fname, GretlType type,
			DATASET *dset, PRN *prn)
{
    fncall *call = NULL;
    GretlType rtype;
    ufunc *func;
    int i, err = 0;

    func = get_user_function_by_name(fname);
    if (func == NULL) {
	return E_DATA;
    }

    call = fncall_new(func, 1);
    if (call == NULL) {
	return E_ALLOC;
    }

    rtype = fncall_get_return_type(call);
    if (rtype != GRETL_TYPE_DOUBLE) {
	fncall_destroy(call);
	return E_TYPES;
    }

    for (i=0; i<2 && !err; i++) {
	err = push_anon_function_arg(call, type, NULL);
    }

    if (err) {
	fncall_destroy(call);
    } else {
	asort.call = call;
	asort.type = type;
	asort.dset = dset;
	asort.prn = prn;
	asort.err = 0;
    }

    return err;
}

static int compare_array_elements (const void *a, const void *b)
{
    void *va = *(void **) a;
    void *vb = *(void **) b;
    double d = 0;
    double *pd = &d;

    if (asort.err) {
	/* get out asap */
	return 0;
    }

#if AS_DEBUG
    fprintf(stderr, "compare: va = %p, vb = %p\n", va, vb);
#endif
    set_anon_function_arg(asort.call, 0, asort.type, va);
    set_anon_function_arg(asort.call, 1, asort.type, vb);
    asort.err = gretl_function_exec(asort.call, GRETL_TYPE_DOUBLE,
				    asort.dset, &pd, asort.prn);
    if (asort.err) {
	fprintf(stderr, " compare_bundles: err=%d, d=%g\n", asort.err, *pd);
    }

    return (int) *pd;
}

int gretl_array_qsort (gretl_array *a, const char *fname,
		       DATASET *dset, PRN *prn)
{
    int n = gretl_array_get_length(a);
    GretlType atype;
    void *aptr;
    int err;

    if (n < 2) {
	/* nothing to be done */
	return 0;
    }

    atype = gretl_array_get_content_type(a);
    err = asort_setup(fname, atype, dset, prn);
    if (err) {
	return err;
    }

    aptr = gretl_array_get_all_data(a);

    set_user_qsorting(1);
    qsort(aptr, n, sizeof(void *), compare_array_elements);
    set_user_qsorting(0);

    err = asort.err;
    asort_cleanup();

    return err;
}

/* end apparatus for sorting a gretl array via qsort */

/* Check for equality between two arrays. Returns 1 if the two arrays
   are of the same type and length, and all their elements compare
   equal, else returns 0.
*/

int gretl_arrays_are_equal (const gretl_array *a,
			    const gretl_array *b)
{
    if (a == b) {
	return 1;
    } else if (a->type != b->type || a->n != b->n) {
	return 0;
    } else {
	int i, eq, nulls;
	int err = 0;

	for (i=0; i<a->n; i++) {
	    if (a->type == GRETL_TYPE_STRINGS) {
		nulls = (a->data[i] == NULL) + (b->data[i] == NULL);
		if (nulls == 1) {
		    return 0;
		} else if (nulls == 0 && strcmp(a->data[i], b->data[i])) {
		    return 0;
		}
	    } else if (a->type == GRETL_TYPE_MATRICES) {
		eq = gretl_matrices_are_equal(a->data[i], b->data[i], 0, &err);
		if (eq != 1) {
		    /* could be 0 or non-conformable (-1) */
		    return eq; /* or return 0? */
		}
	    } else if (a->type == GRETL_TYPE_BUNDLES) {
		if (!gretl_bundles_are_equal(a->data[i], b->data[i])) {
		    return 0;
		}
	    } else if (a->type == GRETL_TYPE_ARRAYS) {
		if (!gretl_arrays_are_equal(a->data[i], b->data[i])) {
		    return 0;
		}
	    } else if (a->type == GRETL_TYPE_LISTS) {
		if (gretl_list_cmp(a->data[i], b->data[i])) {
		    return 0;
		}
	    }
	}
	return 1;
    }
}

/* arglist_validate() is used only in geneval.c but it's defined here
   because (a) geneval.c is already overburdened and (b) the job can
   be done more efficiently with access to the internals of
   gretl_array. We're checking that @args is an array of strings
   holding exactly the strings in @keys minus "arglist" itself.
*/

int arglist_validate (gretl_array *keys, gretl_array *args)
{
    const char *sk, *sa;
    int n_keys = keys->n - 1;
    int i, j, found = 0;

    if (args->type != GRETL_TYPE_STRINGS || args->n != n_keys) {
	return 0;
    }

    for (i=0; i<keys->n; i++) {
	sk = keys->data[i];
	for (j=0; j<args->n; j++) {
	    sa = args->data[j];
	    if (sa == NULL || !strcmp(sa, "arglist")) {
		return 0;
	    } else if (!strcmp(sa, sk)) {
		found++;
		break;
	    }
	}
    }

    return found == n_keys;
}
