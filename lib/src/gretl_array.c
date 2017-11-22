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
	    for (i=0; i<A->n; i++) {
		gretl_matrix_free(A->data[i]);
	    }
	} else {
	    for (i=0; i<A->n; i++) {
		gretl_bundle_destroy(A->data[i]);
	    }
	}   
	free(A->data);
	A->data = NULL;
    }
}

void gretl_array_destroy (gretl_array *A)
{
    if (A != NULL) {
	gretl_array_destroy_data(A);
	free(A);
    }
}

void gretl_array_void_content (gretl_array *A)
{
    if (A != NULL) {
	gretl_array_destroy_data(A);
	A->n = 0;
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
	    for (i=A->n; i<n-1; i++) {
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
	type != GRETL_TYPE_ANY) {
	*err = E_TYPES;
	return NULL;
    } else if (n < 0) {
	*err = E_DATA;
	return NULL;
    } else if (type == GRETL_TYPE_ANY && n > 0) {
	/* an array of unspecified type must be
	   empty in the first instance
	*/
	*err = E_TYPES;
	return NULL;
    }

    A = malloc(sizeof *A);
    
    if (A == NULL) {
	*err = E_ALLOC;
    } else {
	A->type = type;
	A->n = n;
	A->data = NULL;
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

gretl_array *gretl_array_from_strings (char **S, int n,
				       int copy, int *err)
{
    gretl_array *A;

    A = gretl_array_new(GRETL_TYPE_STRINGS, 0, err);

    if (A != NULL) {
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

    if (A->type == GRETL_TYPE_STRINGS) {
	int err = strings_array_null_check(A);

	if (!err) {
	    *ns = A->n;
	    AS = (char **) A->data;
	}
    }

    return AS;
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
	} else {
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

gretl_array *gretl_array_copy_range (gretl_array *A,
				     int r1, int r2,
				     int *err)
{
    gretl_array *C = NULL;

    if (A == NULL) {
	*err = E_DATA;
    } else if (r1 < 0 || r1 >= A->n || r2 < r1 || r2 >= A->n) {
	*err = E_INVARG;
    } else {
	C = gretl_array_new(A->type, r2 - r1 + 1, err);
	if (!*err) {
	    int i, j = 0;

	    for (i=r1; i<=r2 && !*err; i++) {
		if (A->data[i] != NULL) {
		    if (A->type == GRETL_TYPE_STRINGS) {
			C->data[j] = gretl_strdup(A->data[i]);
		    } else if (A->type == GRETL_TYPE_MATRICES) {
			C->data[j] = gretl_matrix_copy(A->data[i]);
		    } else if (A->type == GRETL_TYPE_BUNDLES) {
			C->data[j] = gretl_bundle_copy(A->data[i], err);
		    } else {
			C->data[j] = gretl_list_copy(A->data[i]);
		    }
		    if (!*err && C->data[j] == NULL) {
			*err = E_ALLOC;
		    }
		}
		j++;
	    }
	}
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
    int i;
    
    if (A->type != GRETL_TYPE_MATRICES) {
	*err = E_TYPES;
	return NULL;
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

    if (vcat) {
	ret = gretl_matrix_alloc(sum_r, c0);
    } else {
	ret = gretl_matrix_alloc(r0, sum_c);
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    } else if (vcat) {
	/* vertical concatenation */
	int ii, j, k, rpos = 0;
	double x;

	for (k=0; k<A->n; k++) {
	    m = A->data[k];
	    if (!gretl_is_null_matrix(m)) {
		for (j=0; j<c0; j++) {
		    for (i=0, ii=rpos; i<m->rows; i++, ii++) {
			x = gretl_matrix_get(m, i, j);
			gretl_matrix_set(ret, ii, j, x);
		    }
		}
		rpos += m->rows;
	    }
	}
    } else {
	/* horizontal concatenation */
	double *val = ret->val;
	int n;

	for (i=0; i<A->n; i++) {
	    m = A->data[i];
	    if (!gretl_is_null_matrix(m)) {
		n = m->rows * m->cols;
		memcpy(val, m->val, n * sizeof *val);
		val += n;
	    }
	}
    }
    
    return ret;
}

int gretl_array_set_type (gretl_array *A, GretlType type)
{
    if (A == NULL || A->n > 0) {
	return E_DATA;
    } else if (type != GRETL_TYPE_STRINGS &&
	       type != GRETL_TYPE_MATRICES &&
	       type != GRETL_TYPE_BUNDLES &&
	       type != GRETL_TYPE_LISTS) {
	return E_TYPES;
    } else {
	A->type = type;
	return 0;
    }
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

/* In the functions below I assume the @copy parameter
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
    } else if (A->type != GRETL_TYPE_STRINGS) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else {
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
    } else if (A->type != GRETL_TYPE_STRINGS) {
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
    } else if (A->type != GRETL_TYPE_MATRICES) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else {
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
    } else if (A->type != GRETL_TYPE_MATRICES) {
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
    } else if (A->type != GRETL_TYPE_BUNDLES) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else {
	gretl_bundle_destroy(A->data[i]);
	err = set_bundle(A, i, b, copy);
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
    } else if (A->type != GRETL_TYPE_BUNDLES) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A, 1);
	if (!err) {
	    err = set_bundle(A, A->n - 1, b, copy);
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
    } else if (A->type != GRETL_TYPE_LISTS) {
	err = E_TYPES;
    } else if (i < 0 || i >= A->n) {
	gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i+1);
	err = E_DATA;
    } else {
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
    } else if (A->type != GRETL_TYPE_LISTS) {
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
    }

    return err;
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
	    } else {
		Acpy->data[j] = gretl_list_copy(A->data[i]);
	    }
	    if (Acpy->data[j] == NULL) {
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

/* respond to A1 += A2 */

int gretl_array_append_array (gretl_array *A1,
			      const gretl_array *A2)
{
    int old_n = 0, err = 0;

    if (A1 == NULL || A2 == NULL) {
	err = E_DATA;
    } else if (A1->type != A2->type) {
	err = E_TYPES;
    } else {
	old_n = A1->n;
	err = array_extend_content(A1, A2->n);	
    }

    if (!err) {
	err = gretl_array_copy_content(A1, A2, old_n);
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

static void print_array_elements (gretl_array *A, PRN *prn)
{
    int i;

    for (i=0; i<A->n; i++) {
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
    
    pputc(prn, '\n');
}

int gretl_array_print (gretl_array *A, PRN *prn)
{
    if (A != NULL) {
	const char *s = gretl_type_get_name(A->type);

	pprintf(prn, _("Array of %s, length %d\n"), s, A->n);

	if (A->n > 0 && A->n < 10 && A->type != GRETL_TYPE_BUNDLES) {
	    print_array_elements(A, prn);
	}
    }

    return 0;
}

/* Called from gretl_bundle.c when serializing a bundle
   which contains one or more arrays.
*/

void gretl_array_serialize (gretl_array *A, FILE *fp)
{
    GretlType type;
    const char *subname;
    void *ptr;
    int i;

    type = gretl_type_get_singular(A->type);
    subname = gretl_type_get_name(type);

    fprintf(fp, "<gretl-array type=\"%s\" length=\"%d\">\n",
	    gretl_type_get_name(A->type), A->n); 

    for (i=0; i<A->n; i++) {
	ptr = A->data[i];
	if (ptr == NULL) {
	    fprintf(fp, "<%s placeholder=\"true\"/>\n", subname);
	} else if (type == GRETL_TYPE_STRING) {
	    gretl_xml_put_tagged_string("string", ptr, fp);
	} else if (type == GRETL_TYPE_MATRIX) {
	    gretl_matrix_serialize(ptr, NULL, fp);
	} else if (type == GRETL_TYPE_BUNDLE) {
	    gretl_bundle_serialize(ptr, NULL, fp);
	} else if (type == GRETL_TYPE_LIST) {
	    gretl_xml_put_tagged_list("list", ptr, fp);
	}
    }	    

    fputs("</gretl-array>\n", fp); 
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
	} else if (A->type == GRETL_TYPE_LISTS) {
	    A->data[i] = gretl_xml_get_list(cur, doc, &err);
	}
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
