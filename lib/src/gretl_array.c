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
#include "uservar.h"
#include "gretl_func.h"
#include "gretl_array.h"

union array_data {
    char **S;
    gretl_matrix **M;
    gretl_bundle **B;
};

/**
 * gretl_array:
 *
 * An opaque type; use the relevant accessor functions.
 */

struct gretl_array_ {
    GretlType type;        /* type of data */
    int n;                 /* number of elements */
    union array_data data; /* data array variants */
};

static void gretl_array_destroy_data (gretl_array *A)
{
    int i;

    if (A->type == GRETL_TYPE_STRING_ARRAY) {
	if (A->data.S != NULL) {
	    for (i=0; i<A->n; i++) {
		free(A->data.S[i]);
	    }
	    free(A->data.S);
	    A->data.S = NULL;
	}
    } else if (A->type == GRETL_TYPE_MATRIX_ARRAY) {
	if (A->data.M != NULL) {
	    for (i=0; i<A->n; i++) {
		gretl_matrix_free(A->data.M[i]);
	    }
	    free(A->data.M);
	    A->data.M = NULL;
	}
    } else {
	if (A->data.B != NULL) {
	    for (i=0; i<A->n; i++) {
		gretl_bundle_destroy(A->data.B[i]);
	    }
	    free(A->data.B);
	    A->data.B = NULL;
	}
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

static int array_allocate_content (gretl_array *A)
{
    int i;

    if (A->type == GRETL_TYPE_STRING_ARRAY) {
	A->data.S = malloc(A->n * sizeof *A->data.S);
	if (A->data.S == NULL) {
	    return E_ALLOC;
	} else {
	    for (i=0; i<A->n; i++) {
		A->data.S[i] = NULL;
	    }
	}
    } else if (A->type == GRETL_TYPE_MATRIX_ARRAY) {
	A->data.M = malloc(A->n * sizeof *A->data.M);
	if (A->data.M == NULL) {
	    return E_ALLOC;
	} else {
	    for (i=0; i<A->n; i++) {
		A->data.M[i] = NULL;
	    }
	}
    } else {
	A->data.B = malloc(A->n * sizeof *A->data.B);
	if (A->data.B == NULL) {
	    return E_ALLOC;
	} else {
	    for (i=0; i<A->n; i++) {
		A->data.B[i] = NULL;
	    }
	}	
    }

    return 0;
}

static int array_extend_content (gretl_array *A)
{
    int n = A->n + 1;
    int err = 0;

    if (A->type == GRETL_TYPE_STRING_ARRAY) {
	char **S = 
	    realloc(A->data.S, n * sizeof *A->data.S);

	if (S == NULL) {
	    err = E_ALLOC;
	} else {
	    S[n-1] = NULL;
	    A->data.S = S;
	}
    } else if (A->type == GRETL_TYPE_MATRIX_ARRAY) {
	gretl_matrix **M = 
	    realloc(A->data.M, n * sizeof *A->data.M);

	if (M == NULL) {
	    err = E_ALLOC;
	} else {
	    M[n-1] = NULL;
	    A->data.M = M;
	}
    } else {
	gretl_bundle **B = 
	    realloc(A->data.B, n * sizeof *A->data.B);

	if (B == NULL) {
	    err = E_ALLOC;
	} else {
	    B[n-1] = NULL;
	    A->data.B = B;
	}	
    }

    if (!err) {
	A->n = n;
    }

    return err;
}

/* Create a new array of type @type. The idea is that it's OK
   to have n = 0: in that case the array starts empty but can
   be extended by the += operator. Or it can be sized in
   advance.
*/

gretl_array *gretl_array_new (GretlType type, int n, int *err)
{
    gretl_array *A;

    if (type != GRETL_TYPE_STRING_ARRAY &&
	type != GRETL_TYPE_MATRIX_ARRAY &&
	type != GRETL_TYPE_BUNDLE_ARRAY) {
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
	if (type == GRETL_TYPE_STRING_ARRAY) {
	    A->data.S = NULL;
	} else if (type == GRETL_TYPE_MATRIX_ARRAY) {
	    A->data.M = NULL;
	} else {
	    A->data.B = NULL;
	}
	if (n > 0) {
	    *err = array_allocate_content(A);
	    if (*err) {
		gretl_array_destroy(A);
		A = NULL;
	    }
	}
    }

    return A;
}

void *gretl_array_get_element (gretl_array *A, int i, 
			       GretlType *type,
			       int *err)
{
    void *ret = NULL;

    /* FIXME: between here and geneval.c, decide what exactly
       we want to do about copying or not when we grab an
       array element, to avoid both memory corruption and
       wasted cycles/leakage. Right now we copy here,
       unconditionally, but that's probably a bad idea.
    */

    if (A == NULL || i < 0 || i >= A->n) {
	*err = E_DATA;
    } else {
	if (A->type == GRETL_TYPE_STRING_ARRAY) {
	    *type = GRETL_TYPE_STRING;
	    if (A->data.S[i] == NULL) {
		A->data.S[i] = gretl_strdup("");
	    }
	    ret = gretl_strdup(A->data.S[i]);
	} else if (A->type == GRETL_TYPE_MATRIX_ARRAY) {
	    *type = GRETL_TYPE_MATRIX;
	    if (A->data.M[i] == NULL) {
		A->data.M[i] = gretl_null_matrix_new();
	    }	    
	    ret = gretl_matrix_copy(A->data.M[i]);
	} else {
	    *type = GRETL_TYPE_BUNDLE;
	    if (A->data.B[i] == NULL) {
		A->data.B[i] = gretl_bundle_new();
	    }	    
	    ret = gretl_bundle_copy(A->data.B[i], err);
	}
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

GretlType gretl_array_get_type (gretl_array *A)
{
    if (A != NULL) {
	return A->type;
    } else {
	return GRETL_TYPE_NONE;
    }
}

static int set_string (gretl_array *A, int i,
		       char *s, int copy)
{
    int err = 0;

    if (copy) {
	A->data.S[i] = gretl_strdup(s);
	if (A->data.S[i] == NULL) {
	    err = E_ALLOC;
	}
    } else {
	A->data.S[i] = s;
    }

    return err;
}

static int set_matrix (gretl_array *A, int i,
		       gretl_matrix *m, int copy)
{
    int err = 0;

    if (copy) {
	A->data.M[i] = gretl_matrix_copy(m);
	if (A->data.M[i] == NULL) {
	    err = E_ALLOC;
	}
    } else {
	A->data.M[i] = m;
    }

    return err;
}

static int set_bundle (gretl_array *A, int i,
		       gretl_bundle *b, int copy)
{
    int err = 0;

    if (copy) {
	A->data.B[i] = gretl_bundle_copy(b, &err);
    } else {
	A->data.B[i] = b;
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

    if (A == NULL || i < 0 || i >= A->n) {
	err = E_DATA;
    } else if (A->type != GRETL_TYPE_STRING_ARRAY) {
	err = E_TYPES;
    } else {
	free(A->data.S[i]);
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
    } else if (A->type != GRETL_TYPE_STRING_ARRAY) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A);
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

    if (A == NULL || i < 0 || i >= A->n) {
	err = E_DATA;
    } else if (A->type != GRETL_TYPE_MATRIX_ARRAY) {
	err = E_TYPES;
    } else {
	gretl_matrix_free(A->data.M[i]);
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
    } else if (A->type != GRETL_TYPE_MATRIX_ARRAY) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A);
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

    if (A == NULL || i < 0 || i >= A->n) {
	err = E_DATA;
    } else if (A->type != GRETL_TYPE_BUNDLE_ARRAY) {
	err = E_TYPES;
    } else {
	gretl_bundle_destroy(A->data.B[i]);
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
    } else if (A->type != GRETL_TYPE_BUNDLE_ARRAY) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A);
	if (!err) {
	    err = set_bundle(A, A->n - 1, b, copy);
	}
    }

    return err;
}

static int 
gretl_array_copy_content (gretl_array *Acpy, const gretl_array *A)
{
    int i, err = 0;

    if (A->type == GRETL_TYPE_STRING_ARRAY) {
	for (i=0; i<A->n && !err; i++) {
	    if (A->data.S[i] != NULL) {
		Acpy->data.S[i] = gretl_strdup(A->data.S[i]);
		if (Acpy->data.S[i] == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    } else if (A->type == GRETL_TYPE_MATRIX_ARRAY) {
	for (i=0; i<A->n && !err; i++) {
	    if (A->data.M[i] != NULL) {
		Acpy->data.M[i] = gretl_matrix_copy(A->data.M[i]);
		if (Acpy->data.M[i] == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    } else {
	for (i=0; i<A->n && !err; i++) {
	    if (A->data.B[i] != NULL) {
		Acpy->data.B[i] = 
		    gretl_bundle_copy(A->data.B[i], &err);
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
	    *err = gretl_array_copy_content(Acpy, A);
	}
    }

    return Acpy;
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
	    err = array_allocate_content(A1);
	    if (!err) {
		err = gretl_array_copy_content(A1, A0);
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
	user_var *u = get_user_var_of_type_by_name(name, GRETL_TYPE_ARRAY);

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

int gretl_array_print (gretl_array *A, PRN *prn)
{
    if (A != NULL) {
	const char *s = gretl_arg_type_name(A->type);

	pprintf(prn, _("Array of %s, length %d\n"), s, A->n);
    }

    return 0;
}
