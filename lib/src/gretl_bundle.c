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
#include "gretl_func.h"
#include "libset.h"
#include "gretl_scalar.h"
#include "gretl_string_table.h"
#include "gretl_bundle.h"

#include <glib.h>

#define BDEBUG 0
#define BUNDLE_RETVAL 999

typedef struct gretl_bundle_ gretl_bundle;
typedef struct bundle_value_ bundle_value;

struct gretl_bundle_ {
    char name[VNAMELEN]; 
    GHashTable *ht;      /* holds key/value pairs */
    int level;           /* level of function execution */
};

struct bundle_value_ {
    GretlType type;
    gpointer data;
};

static gretl_bundle **bundles;
static int n_bundles;

/* allocate and fill out a 'value' (type plus data pointer) that will
   be inserted into a bundle's hash table */

static bundle_value *bundle_value_new (GretlType type, void *ptr, int *err)
{
    bundle_value *val = malloc(sizeof *val);

    if (val == NULL) {
	*err = E_ALLOC;
    } else {
	val->type = type;
	switch (val->type) {
	case GRETL_TYPE_DOUBLE:
	    val->data = malloc(sizeof(double));
	    if (val->data != NULL) {
		double *dp = val->data;

		*dp = *(double *) ptr;
	    }
	    break;
	case GRETL_TYPE_STRING:	
	    val->data = gretl_strdup((char *) ptr);
	    break;
	case GRETL_TYPE_MATRIX:
	    val->data = gretl_matrix_copy((gretl_matrix *) ptr);
	    break;
	default:
	    *err = E_TYPES;
	    break;
	}

	if (!*err && val->data == NULL) {
	    free(val);
	    val = NULL;
	    *err = E_ALLOC;
	}
    }

    return val;
}

/* callback invoked when a bundle's hash table is destroyed */

static void bundle_value_destroy (gpointer data)
{
    bundle_value *val = data;

    switch (val->type) {
    case GRETL_TYPE_DOUBLE:
    case GRETL_TYPE_STRING:	
	free(val->data);
	break;
    case GRETL_TYPE_MATRIX:
	gretl_matrix_free((gretl_matrix *) val->data);
	break;
    default:
	break;
    }

    free(val);
}

static void gretl_bundle_free (gretl_bundle *b)
{
    if (b != NULL) {
	if (b->ht != NULL) {
	    g_hash_table_destroy(b->ht);
	}
	free(b);
    }
}

static void set_n_bundles (int n)
{
    n_bundles = n;

    if (n_bundles == 0) {
	free(bundles);
	bundles = NULL;
    }
}

static int reallocate_bundles (int n)
{
    gretl_bundle **tmp;

    tmp = realloc(bundles, n * sizeof *tmp);

    if (tmp == NULL) {
	return E_ALLOC;
    } else {
	bundles = tmp;
	return 0;
    }
}

static gretl_bundle *get_bundle_pointer (const char *name, int level)
{
    gretl_bundle *b;
    int i;

    for (i=0; i<n_bundles; i++) {
	b = bundles[i];
	if (b->level == level && !strcmp(name, b->name)) {
	    return b;
	}
    }

    return NULL;
}

static int gretl_bundle_push (gretl_bundle *b)
{
    int n = n_bundles + 1;

    if (reallocate_bundles(n)) {
	free(b);
	return E_ALLOC;
    }

    bundles[n-1] = b;
    set_n_bundles(n);

    return 0;
}

static int real_delete_bundle (int i)
{
    int n = n_bundles - 1;
    int err = 0;

    if (i < 0 || i > n) {
	return E_DATA;
    }

    gretl_bundle_free(bundles[i]);
    bundles[i] = NULL;

    if (n == 0) {
	set_n_bundles(0);
    } else {
	int j;

	for (j=i; j<n; j++) {
	    bundles[j] = bundles[j+1];
	}
	if (reallocate_bundles(n)) {
	    err = E_ALLOC;
	} else {
	    set_n_bundles(n);
	}
    }

    return err;
}

/* Determine whether @name is the name of a saved bundle. */

int gretl_is_bundle (const char *name)
{
    if (name == NULL || *name == '\0') {
	return 0;
    } else {
	return (get_bundle_pointer(name, gretl_function_depth()) != NULL);
    }
}

/**
 * gretl_bundle_get_data:
 * @name: name of bundle.
 * @key: name of key to access.
 * @type: location to receive data type.
 *
 * Returns: the data pointer associated with @key in the
 * bundle given by @name, or NULL on failure.
 */

void *gretl_bundle_get_data (const char *name, const char *key,
			     GretlType *type)
{
    void *ret = NULL;
    gretl_bundle *b;

    b = get_bundle_pointer(name, gretl_function_depth());

    if (b != NULL) {
	gpointer p = g_hash_table_lookup(b->ht, key);

	if (p != NULL) {
	    bundle_value *val = p;
	    
	    *type = val->type;
	    ret = val->data;
	}
    }

    return ret;
}

static int real_gretl_bundle_set_data (gretl_bundle *b, const char *key,
				       void *ptr, GretlType type)
{
    int err = 0;
    bundle_value *val = bundle_value_new(type, ptr, &err);

    if (!err) {
	gchar *k = g_strdup(key);

	if (g_hash_table_lookup(b->ht, key) != NULL) {
	    g_hash_table_replace(b->ht, k, val);
	} else {
	    g_hash_table_insert(b->ht, k, val);
	}
    }

    return err;
}

/**
 * gretl_bundle_set_data:
 * @name: name of bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: type of data.
 * 
 * Sets the data type and pointer to be associated with @key in 
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_data (const char *name, const char *key,
			   void *ptr, GretlType type)
{
    gretl_bundle *b;
    int err;

    b = get_bundle_pointer(name, gretl_function_depth());
    
    if (b == NULL) {
	err = E_UNKVAR;
    } else {
	err = real_gretl_bundle_set_data(b, key, ptr, type);
    }

    return err;
}

int gretl_bundle_add (const char *name)
{
    int level = gretl_function_depth();
    gretl_bundle *b;

    b = get_bundle_pointer(name, level);
    if (b != NULL) {
	fprintf(stderr, "*** gretl_bundle_add: there's already a '%s' at level %d\n", 
		name, level);
	return E_DATA;
    }

    b = malloc(sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    strcpy(b->name, name);
    b->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
				  g_free, bundle_value_destroy);
    b->level = level;

    return gretl_bundle_push(b);
}

void copy_bundle_value (gpointer key, gpointer value, gpointer p)
{
    bundle_value *bval = (bundle_value *) value;
    gretl_bundle *targ = (gretl_bundle *) p;

    real_gretl_bundle_set_data(targ, (const char *) key,
			       bval->data, bval->type);
}

/* Called from geneval.c on completion of assignment to a
   bundle named @cpyname, where the returned value on the
   right-hand side has not been produced by a user-function
*/

int gretl_bundle_copy_as (const char *name, const char *cpyname)
{
    int level = gretl_function_depth();
    gretl_bundle *b0, *b1;
    int err = 0;

    b0 = get_bundle_pointer(name, level);
    if (b0 == NULL) {
	return E_UNKVAR;
    }    

    b1 = get_bundle_pointer(cpyname, level);

    if (b1 != NULL) {
	g_hash_table_destroy(b1->ht);
	b1->ht = NULL;
    } else {
	b1 = malloc(sizeof *b1);
	if (b1 == NULL) {
	    err = E_ALLOC;
	} else {
	    strcpy(b1->name, cpyname);
	    b1->ht = NULL;
	    b1->level = level;
	    err = gretl_bundle_push(b1);
	}
    }

    if (!err) {
	b1->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
				       g_free, bundle_value_destroy);
	g_hash_table_foreach(b0->ht, copy_bundle_value, b1);
    }

    return err;
}

/* Called from gretl_func.c on return, to mark a given
   bundle as a function return value; we do this by
   setting the bundle's level to a special value.
*/

int gretl_bundle_mark_as_return (const char *name)
{
    gretl_bundle *b;

#if BDEBUG
    fprintf(stderr, "gretl_bundle_mark_as_return: '%s' (depth %d)\n",
	    name, gretl_function_depth());
#endif

    b = get_bundle_pointer(name, gretl_function_depth());
    
    if (b == NULL) {
	return E_DATA;
    } else {
	b->level = BUNDLE_RETVAL;
	return 0;
    }
}

static int get_bundle_index (gretl_bundle *b)
{
    int i;

    for (i=0; i<n_bundles; i++) {
	if (b == bundles[i]) {
	    return i;
	}
    }

    return -1;
}

/* Called from geneval.c to operate on a bundle that has been marked
   as a user-function return value (indentfied by its special "level"
   value). We see if there's already a bundle of the given name at
   caller level, and either overwrite an existing bundle or adjust the
   bundle's name and level.
*/

int gretl_bundle_name_return (const char *name)
{
    gretl_bundle *b0, *b1 = NULL;
    int level = gretl_function_depth();
    int i;

#if BDEBUG
    fprintf(stderr, "gretl_bundle_name_return: '%s'\n", name);
#endif

    for (i=0; i<n_bundles; i++) {
	if (bundles[i]->level == BUNDLE_RETVAL) {
	    b1 = bundles[i];
	    break;
	}
    }

    if (b1 == NULL) {
	fprintf(stderr, "bundle_name_return: no returned bundle found\n");
	return E_DATA;
    }

    b0 = get_bundle_pointer(name, level);
    
    if (b0 != NULL) {
	/* replace */
	g_hash_table_destroy(b0->ht);
	b0->ht = b1->ht;
	b1->ht = NULL;
	real_delete_bundle(get_bundle_index(b1));
    } else {
	/* rename and set level */
	strcpy(b1->name, name);
	b1->level = level;
    }

    return 0;
}

/**
 * gretl_bundle_localize:
 * @origname: name of bundle at caller level.
 * @localname: name to be used within function.
 *
 * On entry to a function, renames the named bundle (provided 
 * as an argument) and sets its level so that is is accessible
 * within the function.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int gretl_bundle_localize (const char *origname,
			   const char *localname)
{
    gretl_bundle *b;
    int err = 0;

    b = get_bundle_pointer(origname, gretl_function_depth());

    if (b == NULL) {
	err = E_DATA;
    } else {
	strcpy(b->name, localname);
	b->level += 1;
    }

    return err;
}

/**
 * gretl_bundle_unlocalize:
 * @localname: name of bundle within function.
 * @origname: name used at caller level.
 *
 * On exit from a function, restores the original name and
 * level of a bundle which has been made available as an argument. 
 * 
 * Returns: 0 on success, non-zero on error.
 */

int gretl_bundle_unlocalize (const char *localname,
			     const char *origname)
{
    gretl_bundle *b;
    int err = 0;

    b = get_bundle_pointer(localname, gretl_function_depth());

    if (b == NULL) {
	err = E_DATA;
    } else {
	strcpy(b->name, origname);
	b->level -= 1;
    }

    return err;
}

/**
 * gretl_bundle_delete:
 * @name: name of bundle to delete.
 * @prn: gretl printing struct.
 *
 * Deletes the named bundle, and prints a message to @prn if
 * messages are turned on.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int gretl_bundle_delete (const char *name, PRN *prn)
{
    int level = gretl_function_depth();
    int i, err = E_UNKVAR;

    for (i=0; i<n_bundles; i++) {
	if (bundles[i]->level == level && 
	    !strcmp(name, bundles[i]->name)) {
	    err = real_delete_bundle(i);
	    break;
	}
    }

    if (!err && prn != NULL && gretl_messages_on()) {
	pprintf(prn, _("Deleted bundle %s"), name);
	pputc(prn, '\n');
    }

    return err;
}

/**
 * destroy_user_bundles_at_level.
 * @level: depth of function execution.
 *
 * Deletes any saved bundles at a particular level of function
 * execution; used on exiting a user-defined function.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_bundles_at_level (int level)
{
    int nb = n_bundles;
    int bmax = nb - 1;
    int i, j;
    int err = 0;

    for (i=bmax; i>=0; i--) {
	if (bundles[i]->level == level) {
	    gretl_bundle_free(bundles[i]);
	    bundles[i] = NULL;
	    for (j=i; j<bmax; j++) {
		bundles[j] = bundles[j+1];
	    }
	    bundles[bmax] = NULL;
	    nb--;
	} 
    }

    if (nb < n_bundles) {
	set_n_bundles(nb);
	if (nb > 0) {
	    err = reallocate_bundles(nb);
	}
    }

    return err;
}

/**
 * destroy_user_bundles:
 *
 * Frees all resources associated with the stack of user-
 * defined bundles.
 */

void destroy_user_bundles (void)
{
    int i;

    if (bundles == NULL) {
	return;
    }

    for (i=0; i<n_bundles; i++) {
	gretl_bundle_free(bundles[i]);
    }

    set_n_bundles(0);
}

