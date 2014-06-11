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

void gretl_array_destroy (gretl_array *A)
{
    if (A != NULL) {
	int i;

	if (A->type == GRETL_TYPE_STRING) {
	    if (A->data.S != NULL) {
		for (i=0; i<A->n; i++) {
		    free(A->data.S[i]);
		}
		free(A->data.S);
	    }
	} else if (A->type == GRETL_TYPE_MATRIX) {
	    if (A->data.M != NULL) {
		for (i=0; i<A->n; i++) {
		    gretl_matrix_free(A->data.M[i]);
		}
		free(A->data.M);
	    }
	} else {
	    if (A->data.B != NULL) {
		for (i=0; i<A->n; i++) {
		    gretl_bundle_destroy(A->data.B[i]);
		}
		free(A->data.B);
	    }
	}	    

	free(A);
    }
}

static int array_allocate_content (gretl_array *A)
{
    int i;

    if (A->type == GRETL_TYPE_STRING) {
	A->data.S = malloc(A->n * sizeof *A->data.S);
	if (A->data.S == NULL) {
	    return E_ALLOC;
	} else {
	    for (i=0; i<A->n; i++) {
		A->data.S[i] = NULL;
	    }
	}
    } else if (A->type == GRETL_TYPE_MATRIX) {
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

    if (A->type == GRETL_TYPE_STRING) {
	char **S = 
	    realloc(A->data.S, n * sizeof *A->data.S);

	if (S == NULL) {
	    err = E_ALLOC;
	} else {
	    S[n-1] = NULL;
	    A->data.S = S;
	}
    } else if (A->type == GRETL_TYPE_MATRIX) {
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

    if (type != GRETL_TYPE_STRING &&
	type != GRETL_TYPE_MATRIX &&
	type != GRETL_TYPE_BUNDLE) {
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
	if (type == GRETL_TYPE_STRING) {
	    A->data.S = NULL;
	} else if (type == GRETL_TYPE_MATRIX) {
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

void *gretl_array_get_member (gretl_array *A, int i, 
			      GretlType *type,
			      int *err)
{
    void *ret = NULL;

    if (A == NULL || i < 0 || i >= A->n) {
	*err = E_DATA;
    } else {
	*type = A->type;
	if (A->type == GRETL_TYPE_STRING) {
	    ret = A->data.S[i];
	} else if (A->type == GRETL_TYPE_MATRIX) {
	    ret = A->data.M[i];
	} else {
	    ret = A->data.B[i];
	}
    }

    return ret;
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
    } else if (A->type != GRETL_TYPE_STRING) {
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
    } else if (A->type != GRETL_TYPE_STRING) {
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
    } else if (A->type != GRETL_TYPE_MATRIX) {
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
    } else if (A->type != GRETL_TYPE_MATRIX) {
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
    } else if (A->type != GRETL_TYPE_BUNDLE) {
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
    } else if (A->type != GRETL_TYPE_BUNDLE) {
	err = E_TYPES;
    } else {
	err = array_extend_content(A);
	if (!err) {
	    err = set_bundle(A, A->n - 1, b, copy);
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
	    int i;

	    if (A->type == GRETL_TYPE_STRING) {
		for (i=0; i<A->n && !*err; i++) {
		    Acpy->data.S[i] = gretl_strdup(A->data.S[i]);
		    if (Acpy->data.S[i] == NULL) {
			*err = E_ALLOC;
		    }
		}
	    } else if (A->type == GRETL_TYPE_MATRIX) {
		for (i=0; i<A->n && !*err; i++) {
		    Acpy->data.M[i] = gretl_matrix_copy(A->data.M[i]);
		    if (Acpy->data.M[i] == NULL) {
			*err = E_ALLOC;
		    }
		}
	    } else {
		for (i=0; i<A->n && !*err; i++) {
		    Acpy->data.B[i] = 
			gretl_bundle_copy(A->data.B[i], err);
		}
	    }
	}
    }

    return Acpy;
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
    pputs(prn, "array\n"); /* FIXME */
    return 0;
}
