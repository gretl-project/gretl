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

/**
 * gretl_bundle:
 *
 * An opaque type; use the relevant accessor functions.
 */

struct gretl_bundle_ {
    char name[VNAMELEN]; 
    GHashTable *ht;      /* holds key/value pairs */
    int level;           /* level of function execution */
};

/**
 * bundled_item:
 *
 * An item of data within a gretl_bundle. This is an
 * opaque type; use the relevant accessor functions.
 */

struct bundled_item_ {
    GretlType type;
    int size;
    gpointer data;
};

static gretl_bundle **bundles;
static int n_saved_bundles;

/* given a bundle pointer, return its 0-based stack position
   or -1 if it's not on the stack */

static int get_bundle_index (gretl_bundle *b)
{
    if (b != NULL) {
	int i;

	for (i=0; i<n_saved_bundles; i++) {
	    if (b == bundles[i]) {
		return i;
	    }
	}
    }

    return -1;
}

/* given a name and function exec level, return the bundle
   pointer corresponding to @name or NULL if there's no 
   such bundle */

static gretl_bundle *get_bundle_pointer (const char *name, int level)
{
    gretl_bundle *b;
    int i;

    for (i=0; i<n_saved_bundles; i++) {
	b = bundles[i];
	if (b->level == level && !strcmp(name, b->name)) {
	    return b;
	}
    }

    return NULL;
}

static int bundle_index_by_name (const char *name, int level)
{
    return get_bundle_index(get_bundle_pointer(name, level));
}

/* allocate and fill out a 'value' (type plus data pointer) that will
   be inserted into a bundle's hash table */

static bundled_item *bundled_item_new (GretlType type, void *ptr, 
				       int size, int *err)
{
    bundled_item *item = malloc(sizeof *item);

    if (item == NULL) {
	*err = E_ALLOC;
    } else {
	item->type = type;
	item->size = 0;
	switch (item->type) {
	case GRETL_TYPE_DOUBLE:
	    item->data = malloc(sizeof(double));
	    if (item->data != NULL) {
		double *dp = item->data;

		*dp = *(double *) ptr;
	    }
	    break;
	case GRETL_TYPE_STRING:	
	    item->data = gretl_strdup((char *) ptr);
	    break;
	case GRETL_TYPE_MATRIX:
	    item->data = gretl_matrix_copy((gretl_matrix *) ptr);
	    break;
	case GRETL_TYPE_SERIES:
	    item->data = copyvec((const double *) ptr, size);
	    item->size = size;
	    break;
	default:
	    *err = E_TYPES;
	    break;
	}

	if (!*err && item->data == NULL) {
	    free(item);
	    item = NULL;
	    *err = E_ALLOC;
	}
    }

    return item;
}

/* callback invoked when a bundle's hash table is destroyed */

static void bundled_item_destroy (gpointer data)
{
    bundled_item *item = data;

    switch (item->type) {
    case GRETL_TYPE_DOUBLE:
    case GRETL_TYPE_STRING:
    case GRETL_TYPE_SERIES:
	free(item->data);
	break;
    case GRETL_TYPE_MATRIX:
	gretl_matrix_free((gretl_matrix *) item->data);
	break;
    default:
	break;
    }

    free(item);
}

/**
 * gretl_bundle_destroy:
 * @bundle: bundle to destroy.
 *
 * Frees all contents of @bundle as well as the pointer itself.
 */

void gretl_bundle_destroy (gretl_bundle *bundle)
{
    if (bundle != NULL) {
	if (bundle->ht != NULL) {
	    g_hash_table_destroy(bundle->ht);
	}
	free(bundle);
    }
}

/*
 * gretl_bundle_new:
 *
 * Returns: a newly allocated, empty bundle.
 */

static gretl_bundle *gretl_bundle_new (void)
{
    gretl_bundle *b = malloc(sizeof *b);

    if (b != NULL) {
	b->level = 0;
	*b->name = '\0';
	b->ht = NULL;
    }

    return b;
}

static void set_n_saved_bundles (int n)
{
    n_saved_bundles = n;

    if (n_saved_bundles == 0) {
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

/* push a bundle onto the saved stack: note that if this
   fails we destroy the bundle and return non-zero */

static int gretl_bundle_push (gretl_bundle *b)
{
    int n = n_saved_bundles + 1;

    if (reallocate_bundles(n)) {
	gretl_bundle_destroy(b);
	return E_ALLOC;
    }

    bundles[n-1] = b;
    set_n_saved_bundles(n);

    return 0;
}

/* destroy a stacked bundle indentified by its stack position,
   and resize the stack */

static int real_delete_bundle (int i)
{
    int n = n_saved_bundles - 1;
    int err = 0;

    if (i < 0 || i > n) {
	return E_DATA;
    }

    gretl_bundle_destroy(bundles[i]);
    bundles[i] = NULL;

    if (n == 0) {
	set_n_saved_bundles(0);
    } else {
	int j;

	for (j=i; j<n; j++) {
	    bundles[j] = bundles[j+1];
	}
	if (reallocate_bundles(n)) {
	    err = E_ALLOC;
	} else {
	    set_n_saved_bundles(n);
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
 * get_gretl_bundle_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved gretl bundle, if found, else NULL.
 */

gretl_bundle *get_gretl_bundle_by_name (const char *name)
{
    if (name == NULL || *name == '\0') {
	return NULL;
    } else {
	return get_bundle_pointer(name, gretl_function_depth());
    }
}

static int gretl_bundle_has_data (gretl_bundle *b, const char *key)
{
    gpointer p = g_hash_table_lookup(b->ht, key);

    return (p != NULL);
}

/**
 * gretl_bundle_get_data:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0).
 *
 * Returns: the data pointer associated with @key in the
 * specified @bundle, or NULL on failure.
 */

void *gretl_bundle_get_data (gretl_bundle *bundle, const char *key,
			     GretlType *type, int *size)
{
    void *ret = NULL;

    if (bundle != NULL) {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;
	    
	    *type = item->type;
	    ret = item->data;
	    *size = item->size;
	}
    }

    return ret;
}

/**
 * bundle_item_get_data:
 * @item: bundled item to access.
 * @type: location to receive data type.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0).
 *
 * Returns: the data pointer associated with @item, or
 * NULL on failure.
 */

void *bundled_item_get_data (bundled_item *item, GretlType *type,
			     int *size)
{
    *type = item->type;
    *size = item->size;
    return item->data;
}

/**
 * gretl_bundle_get_content:
 * @bundle: bundle to access.
 *
 * Returns: the content of @bundle, which is in fact
 * a GHashTable object.
 */

void *gretl_bundle_get_content (gretl_bundle *bundle)
{
    return (bundle == NULL)? NULL : (void *) bundle->ht;
}

static int real_gretl_bundle_set_data (gretl_bundle *b, const char *key,
				       void *ptr, GretlType type,
				       int size)
{
    int err = 0;
    bundled_item *item = bundled_item_new(type, ptr, size, &err);

    if (!err) {
	gchar *k = g_strdup(key);

	if (g_hash_table_lookup(b->ht, key) != NULL) {
	    g_hash_table_replace(b->ht, k, item);
	} else {
	    g_hash_table_insert(b->ht, k, item);
	}
    }

    return err;
}

/**
 * gretl_bundle_set_data:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @ptr: data pointer.
 * @type: type of data.
 * @size: if @type == GRETL_TYPE_SERIES, the length of
 * the series, otherwise 0.
 * 
 * Sets the data type and pointer to be associated with @key in 
 * the bundle given by @name. If @key is already present in
 * the bundle's hash table the original value is replaced
 * and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_data (gretl_bundle *bundle, const char *key,
			   void *ptr, GretlType type, int size)
{
    int err;

    if (bundle == NULL) {
	err = E_UNKVAR;
    } else {
	err = real_gretl_bundle_set_data(bundle, key, ptr, type, size);
    }

    return err;
}

/**
 * gretl_bundle_add:
 * @name: name to give the bundle.
 *
 * Creates an empty gretl bundle with the given name and pushes it onto
 * the stack of saved bundles.
 */

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
				  g_free, bundled_item_destroy);
    b->level = level;

    return gretl_bundle_push(b);
}

/* For use by geneval: take @bundle, created on the fly, and stack
   it under @name.  If there's already a bundle with the given name
   destroy and replace it.
*/

int gretl_bundle_add_or_replace (gretl_bundle *bundle, const char *name)
{
    int level = gretl_function_depth();
    int b0idx = bundle_index_by_name(name, level);
    int err = 0;

    strcpy(bundle->name, name);
    bundle->level = level;

    if (b0idx >= 0) {
	/* just replace on stack */
	gretl_bundle_destroy(bundles[b0idx]);
	bundles[b0idx] = bundle;
    } else {
	err = gretl_bundle_push(bundle);
    }

    return err;
}

/* replicate on a target bundle a bundled_item from some other
   other bundle, provided the target bundle does not already
   have a bundled_item under the same key
*/

static void copy_new_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = (bundled_item *) value;
    gretl_bundle *targ = (gretl_bundle *) p;

    if (!gretl_bundle_has_data(targ, (const char *) key)) {
	real_gretl_bundle_set_data(targ, (const char *) key,
				   item->data, item->type,
				   item->size);
    }
}

/* replicate on a target bundle a bundled_item from some other
   bundle */

static void copy_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = (bundled_item *) value;
    gretl_bundle *targ = (gretl_bundle *) p;

    real_gretl_bundle_set_data(targ, (const char *) key,
			       item->data, item->type,
			       item->size);
}

/* Create a new bundle as the union of two existing bundles:
   we first copy bundle1 in its entirety, then append any elements
   of bundle2 whose keys that are not already present in the
   copy-target. In case bundle1 and bundle2 share any keys, the
   value copied to the target is therefore that from bundle1.
*/

gretl_bundle *gretl_bundle_union (const gretl_bundle *bundle1,
				  const gretl_bundle *bundle2,
				  int *err)
{
    gretl_bundle *b = gretl_bundle_new();

    if (b == NULL) {
	*err = E_ALLOC;
    } else {
	b->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
				      g_free, bundled_item_destroy);
	g_hash_table_foreach(bundle1->ht, copy_bundled_item, b);
	g_hash_table_foreach(bundle2->ht, copy_new_bundled_item, b);
    }
    
    return b;
}

/* Called from geneval.c on completion of assignment to a
   bundle named @cpyname, where the returned value on the
   right-hand side is a pre-existing named bundle that has 
   not been produced by a user-function
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
				       g_free, bundled_item_destroy);
	g_hash_table_foreach(b0->ht, copy_bundled_item, b1);
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

    for (i=0; i<n_saved_bundles; i++) {
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
 * gretl_bundle_delete_by_name:
 * @name: name of bundle to delete.
 * @prn: gretl printing struct.
 *
 * Deletes the named bundle, and prints a message to @prn if
 * messages are turned on.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int gretl_bundle_delete_by_name (const char *name, PRN *prn)
{
    int level = gretl_function_depth();
    int i, err = E_UNKVAR;

    for (i=0; i<n_saved_bundles; i++) {
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
 * destroy_saved_bundles_at_level.
 * @level: depth of function execution.
 *
 * Deletes any saved bundles at a particular level of function
 * execution; used on exiting a user-defined function.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int destroy_saved_bundles_at_level (int level)
{
    int nb = n_saved_bundles;
    int bmax = nb - 1;
    int i, j;
    int err = 0;

    for (i=bmax; i>=0; i--) {
	if (bundles[i]->level == level) {
	    gretl_bundle_destroy(bundles[i]);
	    bundles[i] = NULL;
	    for (j=i; j<bmax; j++) {
		bundles[j] = bundles[j+1];
	    }
	    bundles[bmax] = NULL;
	    nb--;
	} 
    }

    if (nb < n_saved_bundles) {
	set_n_saved_bundles(nb);
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

    for (i=0; i<n_saved_bundles; i++) {
	gretl_bundle_destroy(bundles[i]);
    }

    set_n_saved_bundles(0);
}

