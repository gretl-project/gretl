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
#include "gretl_func.h"
#include "libset.h"
#include "uservar.h"
#include "gretl_xml.h"
#include "gretl_foreign.h"

#include <glib.h>

#define BDEBUG 0

/**
 * gretl_bundle:
 *
 * An opaque type; use the relevant accessor functions.
 */

struct gretl_bundle_ {
    GHashTable *ht;  /* holds key/value pairs */
    char *creator;   /* name of function that built the bundle */   
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
    char *note;
};

static gretl_bundle *sysinfo_bundle;

static int real_bundle_set_data (gretl_bundle *b, const char *key,
				 void *ptr, GretlType type,
				 int size, const char *note);

int gretl_bundle_set_name (gretl_bundle *b, const char *name)
{
    user_var *u = get_user_var_by_data(b);
    
    if (u == NULL) {
	return E_UNKVAR;
    } else {
	user_var_set_name(u, name);
	return 0;
    }
}

int gretl_bundle_is_stacked (gretl_bundle *b)
{
    user_var *u = get_user_var_by_data(b);

    return u != NULL;
}

int gretl_bundle_get_n_keys (gretl_bundle *b)
{
    int n_items = 0;

    if (b != NULL && b->ht != NULL) {
	n_items = g_hash_table_size(b->ht);
    }

    return n_items;
}

int type_can_be_bundled (GretlType type)
{
    if (type == GRETL_TYPE_INT ||
	type == GRETL_TYPE_BOOL) {
	type = GRETL_TYPE_DOUBLE;
    }

    return (type == GRETL_TYPE_DOUBLE ||
	    type == GRETL_TYPE_STRING ||
	    type == GRETL_TYPE_MATRIX ||
	    type == GRETL_TYPE_MATRIX_REF ||
	    type == GRETL_TYPE_SERIES ||
	    type == GRETL_TYPE_BUNDLE);
}

/* allocate and fill out a 'value' (type plus data pointer) that will
   be inserted into a bundle's hash table */

static bundled_item *bundled_item_new (GretlType type, void *ptr, 
				       int size, const char *note,
				       int *err)
{
    bundled_item *item = malloc(sizeof *item);

    if (item == NULL) {
	*err = E_ALLOC;
    } else {
	item->type = type;
	item->size = 0;
	item->note = NULL;

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
	case GRETL_TYPE_MATRIX_REF:
	    item->data = ptr;
	    break;
	case GRETL_TYPE_SERIES:
	    item->data = copyvec((const double *) ptr, size);
	    item->size = size;
	    break;
	case GRETL_TYPE_BUNDLE:
	    item->data = gretl_bundle_copy((gretl_bundle *) ptr, err);
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

	if (item != NULL && note != NULL) {
	    item->note = gretl_strdup(note);
	}	
    }

    return item;
}

/* callback invoked when a bundle's hash table is destroyed */

static void bundled_item_destroy (gpointer data)
{
    bundled_item *item = data;

#if BDEBUG
    fprintf(stderr, "bundled_item_destroy: type %d\n", item->type);
    if (item->type == GRETL_TYPE_STRING) {
	fprintf(stderr, " string: '%s'\n", (char *) item->data);
    }
#endif

    switch (item->type) {
    case GRETL_TYPE_DOUBLE:
    case GRETL_TYPE_STRING:
    case GRETL_TYPE_SERIES:
	free(item->data);
	break;
    case GRETL_TYPE_MATRIX:
	gretl_matrix_free((gretl_matrix *) item->data);
	break;
    case GRETL_TYPE_MATRIX_REF:
	item->data = NULL;
	break;
    case GRETL_TYPE_BUNDLE:
	gretl_bundle_destroy((gretl_bundle *) item->data);
	break;
    default:
	break;
    }

    free(item->note);
    free(item);
}

static void bundle_key_destroy (gpointer data)
{
#if BDEBUG
    fprintf(stderr, "freeing key '%s'\n", (char *) data);
#endif
    free(data);
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
	free(bundle->creator);
	free(bundle);
    }
}

/**
 * gretl_bundle_void_content:
 * @bundle: target bundle.
 *
 * Frees all contents of @bundle.
 */

void gretl_bundle_void_content (gretl_bundle *bundle)
{
    if (bundle != NULL) {
	if (bundle->ht != NULL) {
	    g_hash_table_destroy(bundle->ht);
	}
	free(bundle->creator);
	bundle->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
					   bundle_key_destroy, 
					   bundled_item_destroy);
	bundle->creator = NULL;
    }
}

/**
 * gretl_bundle_new:
 *
 * Returns: a newly allocated, empty gretl bundle.
 */

gretl_bundle *gretl_bundle_new (void)
{
    gretl_bundle *b = malloc(sizeof *b);

    if (b != NULL) {
	b->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
				      bundle_key_destroy, 
				      bundled_item_destroy);
	b->creator = NULL;
    }

    return b;
}

/* Determine whether @name is the name of a saved bundle. */

int gretl_is_bundle (const char *name)
{
    return get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE) != NULL;
}

/**
 * get_bundle_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved bundle, if found, else NULL.
 */

gretl_bundle *get_bundle_by_name (const char *name)
{
    gretl_bundle *b = NULL;

    if (name != NULL && *name != '\0') {
	user_var *u = get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE);

	if (u != NULL) {
	    b = user_var_get_value(u);
	}
    }

    return b;
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
 * @type: location to receive data type, or NULL.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0), or NULL.
 * @err: location to receive error code, or NULL.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 * If @err is non-NULL, its content is set to a non-zero
 * value if @bundle contains no item with key @key. If
 * the intent is simply to determine whether @bundle contains 
 * an item under the specified @key, @err should generally be 
 * left NULL.
 *
 * Note that the value returned is the actual data pointer from
 * within the bundle, not a copy of the data; so the pointer
 * must not be freed, and in general its content should not
 * be modified.
 */

void *gretl_bundle_get_data (gretl_bundle *bundle, const char *key,
			     GretlType *type, int *size, int *err)
{
    void *ret = NULL;

    if (bundle == NULL) {
	if (err != NULL) {
	    *err = E_DATA;
	}
    } else {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;

	    ret = item->data;
	    if (type != NULL) {
		*type = item->type;
	    }
	    if (size != NULL) {
		*size = item->size;
	    }
	} else if (err != NULL) {
	    gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	    *err = E_DATA;
	}
    }

    return ret;
}

/**
 * gretl_bundle_steal_data:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0), or NULL.
 * @err: location to receive error code, or NULL.
 *
 * Works like gretl_bundle_get_data() except that the data
 * pointer in question is removed from @bundle before it is
 * returned to the caller; so in effect the caller assumes
 * ownership of the item.
 *
 * Returns: the item pointer associated with @key in the
 * specified @bundle, or NULL if there is no such item.
 */

void *gretl_bundle_steal_data (gretl_bundle *bundle, const char *key,
			       GretlType *type, int *size, int *err)
{
    void *ret = NULL;

    if (bundle == NULL) {
	if (err != NULL) {
	    *err = E_DATA;
	}
    } else {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    GList *keys = g_hash_table_get_keys(bundle->ht);
	    bundled_item *item = p;
	    gchar *keycpy = NULL;

	    ret = item->data;
	    if (type != NULL) {
		*type = item->type;
	    }
	    if (size != NULL) {
		*size = item->size;
	    }
	    while (keys) {
		if (!strcmp(keys->data, key)) {
		    keycpy = keys->data;
		    break;
		} else if (keys->next) {
		    keys = keys->next;
		} else {
		    break;
		}
	    }
	    g_hash_table_steal(bundle->ht, key);
	    g_free(keycpy);
	    g_list_free(keys);
	    free(item);
	} else if (err != NULL) {
	    gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	    *err = E_DATA;
	}	    
    }

    return ret;
}

/**
 * gretl_bundle_get_type:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err:location to receive error code.
 *
 * Returns: the data type associated with @key in the
 * specified @bundle, or 0 on failure.
 */

GretlType gretl_bundle_get_type (gretl_bundle *bundle, const char *key,
				 int *err)
{
    GretlType ret = GRETL_TYPE_NONE;

    if (bundle == NULL) {
	*err = E_DATA;
    } else {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;
	    
	    ret = item->type;
	} else {
	    gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
	    *err = E_DATA;
	}
				 
    }

    return ret;
}

/**
 * gretl_bundle_get_matrix:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the matrix associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

gretl_matrix *gretl_bundle_get_matrix (gretl_bundle *bundle,
				       const char *key,
				       int *err)
{
    gretl_matrix *m = NULL;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (!*err && type != GRETL_TYPE_MATRIX && 
	type != GRETL_TYPE_MATRIX_REF) {
	*err = E_TYPES;
    }

    if (!*err) {
	m = (gretl_matrix *) ptr;
    }

    return m;
}

/**
 * gretl_bundle_get_series:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @n: location to receive length of series.
 * @err: location to receive error code.
 *
 * Returns: the series associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

double *gretl_bundle_get_series (gretl_bundle *bundle,
				 const char *key,
				 int *n, int *err)
{
    double *x = NULL;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, n, err);
    if (!*err && type != GRETL_TYPE_SERIES) {
	*err = E_TYPES;
    }

    if (!*err) {
	x = (double *) ptr;
    }

    return x;
}

/**
 * gretl_bundle_get_payload_matrix:
 * @bundle: bundle to access.
 * @err: location to receive error code.
 *
 * Returns: a copy of the matrix associated with the key 
 * "payload" in the specified @bundle, if any; otherwise NULL.
 */

gretl_matrix *gretl_bundle_get_payload_matrix (gretl_bundle *bundle,
					       int *err)
{
    gretl_matrix *m = NULL;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, "payload", &type, 
				NULL, err);
    if (!*err && type != GRETL_TYPE_MATRIX) {
	*err = E_TYPES;
    }

    if (!*err) {
	m = gretl_matrix_copy((gretl_matrix *) ptr);
	if (m == NULL) {
	    *err = E_ALLOC;
	}
    }

    return m;
}

/**
 * gretl_bundle_set_payload_matrix:
 * @bundle: bundle to access.
 * @m: matrix to set as payload.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_bundle_set_payload_matrix (gretl_bundle *bundle,
				     const gretl_matrix *m)
{
    return real_bundle_set_data(bundle, "payload", (void *) m,
				GRETL_TYPE_MATRIX, 0, NULL);
}

/**
 * gretl_bundle_get_scalar:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the scalar value associated with @key in the
 * specified @bundle, if any; otherwise #NADBL.
 */

double gretl_bundle_get_scalar (gretl_bundle *bundle,
				const char *key,
				int *err)
{
    double x = NADBL;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (!*err && type != GRETL_TYPE_DOUBLE) {
	*err = E_TYPES;
    }

    if (!*err) {
	double *px = (double *) ptr;

	x = *px;
    }

    return x;
}

/**
 * gretl_bundle_get_string:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the string value associated with @key in the
 * specified @bundle, if any; otherwise #NADBL.
 */

const char *gretl_bundle_get_string (gretl_bundle *bundle,
				     const char *key,
				     int *err)
{
    const char *ret = NULL;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (!*err && type != GRETL_TYPE_STRING) {
	*err = E_TYPES;
    }

    if (!*err) {
	ret = (const char *) ptr;
    }

    return ret;
}

/**
 * gretl_bundle_get_note:
 * @bundle: bundle to access.
 * @key: name of key to access.
 *
 * Returns: the note associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const char *gretl_bundle_get_note (gretl_bundle *bundle, 
				   const char *key)
{
    const char *ret = NULL;

    if (bundle != NULL) {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p != NULL) {
	    bundled_item *item = p;

	    ret = item->note;
	}
    }

    return ret;
}

/**
 * gretl_bundle_get_creator:
 * @bundle: bundle to access.
 *
 * Returns: the name of the package that created @bundle, if any, 
 * otherwise NULL.
 */

const char *gretl_bundle_get_creator (gretl_bundle *bundle)
{
    return (bundle != NULL)? bundle->creator : NULL;
}

/**
 * bundled_item_get_data:
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
 * bundled_item_get_note:
 * @item: bundled item.
 *
 * Returns: the note associated with @item (may be NULL).
 */

const char *bundled_item_get_note (bundled_item *item)
{
    return item->note;
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

static int real_bundle_set_data (gretl_bundle *b, const char *key,
				 void *ptr, GretlType type,
				 int size, const char *note)
{
    bundled_item *item = NULL;
    int err;

    err = strlen(key) >= VNAMELEN ? E_DATA : 0;

    if (err) {
	gretl_errmsg_sprintf("'%s': invalid key string", key);
    } else {
	item = bundled_item_new(type, ptr, size, note, &err);
    }

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
	err = real_bundle_set_data(bundle, key, ptr, type, size, NULL);
    }

    return err;
}

/**
 * gretl_bundle_set_string:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @str: the string to set.
 * 
 * Sets @str as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_string (gretl_bundle *bundle, const char *key,
			     const char *str)
{
    return gretl_bundle_set_data(bundle, key, (void *) str, 
				 GRETL_TYPE_STRING, 0);
}

/**
 * gretl_bundle_set_scalar:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the value to set.
 * 
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_scalar (gretl_bundle *bundle, const char *key,
			     double val)
{
    return gretl_bundle_set_data(bundle, key, &val, 
				 GRETL_TYPE_DOUBLE, 0);
}

/**
 * gretl_bundle_set_series:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @x: array of doubles.
 * @n: the length of @x.
 * 
 * Sets @x as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_series (gretl_bundle *bundle, const char *key,
			     const double *x, int n)
{
    return gretl_bundle_set_data(bundle, key, (void *) x, 
				 GRETL_TYPE_SERIES, n);
}

/**
 * gretl_bundle_set_matrix:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @m: gretl matrix.
 * 
 * Sets @m as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original 
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_matrix (gretl_bundle *bundle, const char *key,
			     const gretl_matrix *m)
{
    return gretl_bundle_set_data(bundle, key, (void *) m, 
				 GRETL_TYPE_MATRIX, 0);
}

/**
 * gretl_bundle_delete_data:
 * @bundle: target bundle.
 * @key: name of key to delete.
 * 
 * Deletes the data item under @key from @bundle, if
 * such an item is present.
 */

int gretl_bundle_delete_data (gretl_bundle *bundle, const char *key)
{
    int err = 0;

    if (bundle == NULL) {
	err = E_DATA;
    } else {
	gboolean ok = g_hash_table_remove(bundle->ht, key);
	
	if (!ok) {
	    err = E_DATA;
	}
    }
    
    return err;
}

/**
 * gretl_bundle_set_note:
 * @bundle: target bundle.
 * @key: name of key to access.
 * @note: note to add.
 * 
 * Adds a descriptive note to the item under @key in @bundle.
 * If a note is already present it is replaced by the new one.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_note (gretl_bundle *bundle, const char *key,
			   const char *note)
{
    int err = 0;

    if (bundle == NULL) {
	err = E_UNKVAR;
    } else {
	gpointer p = g_hash_table_lookup(bundle->ht, key);

	if (p == NULL) {
	    err = E_DATA;
	} else {
	    bundled_item *item = p;

	    free(item->note);
	    item->note = gretl_strdup(note);
	}
    }

    return err;
}

/**
 * gretl_bundle_set_creator:
 * @bundle: target bundle.
 * @name: name of function package that built the bundle.
 * 
 * Sets the "creator" attribute of @bundle. This is called 
 * automatically when a bundle is returned to top-level 
 * userspace by a packaged function.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_creator (gretl_bundle *bundle, const char *name)
{
    int err = 0;

    if (bundle == NULL) {
	err = E_DATA;
    } else {
	free(bundle->creator);
	if (name == NULL) {
	    bundle->creator = NULL;
	} else {
	    bundle->creator = gretl_strdup(name);
	    if (bundle->creator == NULL) {
		err = E_ALLOC;
	    }
	}
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
	real_bundle_set_data(targ, (const char *) key,
			     item->data, item->type,
			     item->size, item->note);
    }
}

/* replicate on a target bundle a bundled_item from some other
   bundle */

static void copy_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = (bundled_item *) value;
    gretl_bundle *targ = (gretl_bundle *) p;

    real_bundle_set_data(targ, (const char *) key,
			 item->data, item->type,
			 item->size, item->note);
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
	g_hash_table_foreach(bundle1->ht, copy_bundled_item, b);
	g_hash_table_foreach(bundle2->ht, copy_new_bundled_item, b);
    }
    
    return b;
}

/**
 * gretl_bundle_copy_as:
 * @name: name of source bundle.
 * @copyname: name for copy.
 *
 * Look for a saved bundle specified by @name, and if found,
 * make a full copy and save it under @copyname. This is
 * called from geneval.c on completion of assignment to a
 * bundle named @copyname, where the returned value on the
 * right-hand side is a pre-existing saved bundle.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_bundle_copy_as (const char *name, const char *copyname)
{
    gretl_bundle *b0, *b1 = NULL;
    user_var *u;
    int err = 0;

    if (!strcmp(name, "$sysinfo")) {
	b0 = sysinfo_bundle;
    } else {
	u = get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE);
	if (u == NULL) {
	    return E_DATA;
	} else {
	    b0 = user_var_get_value(u);
	}
    }

    u = get_user_var_of_type_by_name(copyname, GRETL_TYPE_BUNDLE);
    if (u != NULL) {
	b1 = user_var_get_value(u);
    }

    if (b1 != NULL) {
	g_hash_table_destroy(b1->ht);
	b1->ht = g_hash_table_new_full(g_str_hash, 
				       g_str_equal, 
				       bundle_key_destroy, 
				       bundled_item_destroy);
    } else {
	b1 = gretl_bundle_new();
	if (b1 == NULL) {
	    err = E_ALLOC;
	} else {
	    err = user_var_add(copyname, GRETL_TYPE_BUNDLE, b1);
	}
    }

    if (!err) {
	g_hash_table_foreach(b0->ht, copy_bundled_item, b1);
    }

    return err;
}

/**
 * gretl_bundle_copy:
 * @bundle: gretl bundle to be copied.
 * @err: location to receive error code.
 *
 * Returns: a "deep copy" of @bundle (all the items in @bundle
 * are themselves copied), or NULL on failure.
 */

gretl_bundle *gretl_bundle_copy (const gretl_bundle *bundle, int *err)
{
    gretl_bundle *bcpy = NULL;

    if (bundle == NULL) {
	*err = E_DATA;
    } else {
	bcpy = gretl_bundle_new();
	if (bcpy == NULL) {
	    *err = E_ALLOC;
	} else {
	    g_hash_table_foreach(bundle->ht, copy_bundled_item, bcpy);
	}
    }

    return bcpy;
}

static void print_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = value;
    const gchar *kstr = key;
    gretl_matrix *m;
    double x;
    char *s;
    PRN *prn = p;

    switch (item->type) {
    case GRETL_TYPE_DOUBLE:
	x = *(double *) item->data;
	if (na(x)) {
	    pprintf(prn, " %s = NA", kstr);
	} else {
	    pprintf(prn, " %s = %g", kstr, x);
	}	
	break;
    case GRETL_TYPE_STRING:
	s = (char *) item->data;
	if (strlen(s) < 64) {
	    pprintf(prn, " %s = %s", kstr, s);
	} else {
	    pprintf(prn, " %s (%s)", kstr, 
		    gretl_arg_type_name(item->type));
	}	
	break;
    case GRETL_TYPE_BUNDLE:
	pprintf(prn, " %s (%s)", kstr, 
		gretl_arg_type_name(item->type));
	break;
    case GRETL_TYPE_MATRIX:
    case GRETL_TYPE_MATRIX_REF:
	m = item->data;
	if (m->rows == 1 && m->cols == 1) {
	    pprintf(prn, " %s = %g", kstr, m->val[0]);
	} else {
	    pprintf(prn, " %s (%s: %d x %d)", kstr, 
		    gretl_arg_type_name(item->type), 
		    m->rows, m->cols);
	}
	break;
    case GRETL_TYPE_SERIES:
	pprintf(prn, " %s (%s: length %d)", kstr, 
		gretl_arg_type_name(item->type), item->size);
	break;
    default:
	break;
    }

    if (item->note != NULL) {
	pprintf(prn, " %s\n", item->note);
    } else {
	pputc(prn, '\n');
    }
}

/**
 * gretl_bundle_print:
 * @bundle: gretl bundle.
 * @prn: gretl printer.
 *
 * Prints to @prn a list of the keys defined in @bundle, along
 * with descriptive notes, if any.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int gretl_bundle_print (gretl_bundle *bundle, PRN *prn)
{
    if (bundle == NULL) {
	return E_DATA;
    } else {
	int n_items = g_hash_table_size(bundle->ht);
	user_var *u = get_user_var_by_data(bundle);
	const char *name = NULL;

	if (u != NULL) {
	    name = user_var_get_name(u);
	} else {
	    name = "anonymous";
	}

	if (n_items == 0) {
	    pprintf(prn, "bundle %s: empty\n", name);
	} else {
	    if (bundle->creator != NULL) {
		pprintf(prn, "bundle %s, created by %s:\n", 
			name, bundle->creator);
	    } else {
		pprintf(prn, "bundle %s:\n", name);
	    }
	    g_hash_table_foreach(bundle->ht, print_bundled_item, prn);
	    pputc(prn, '\n');
	}
	return 0;
    }
}

static gboolean match_by_data (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = value;

    return item->data == p;
}

int bundle_contains_data (gretl_bundle *b, void *data)
{
    return g_hash_table_find(b->ht, match_by_data, data) != NULL;
}

/* Called from gretl_func.c on return, to remove
   a given bundle from the stack of named bundles in
   preparation for handing it over to the caller,
   who will take ownership of it.
*/

gretl_bundle *gretl_bundle_pull_from_stack (const char *name,
					    int *err)
{
    gretl_bundle *b = NULL;
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_BUNDLE);

    if (u != NULL) {
	b = user_var_unstack_value(u);
    }

    if (b == NULL) {
	*err = E_DATA;
    } 

    return b;
}

/* serialize a gretl bundled item as XML */

static void xml_put_bundled_item (gpointer keyp, gpointer value, gpointer p)
{
    const char *key = keyp;
    bundled_item *item = value;
    PRN *prn = p;
    int j;

    if (item->type == GRETL_TYPE_STRING) {
	char *s = item->data;

	if (s == NULL || *s == '\0') {
	    fprintf(stderr, "bundle -> XML: skipping empty string %s\n", key);
	    return;
	} 
    }

    pprintf(prn, "<bundled-item key=\"%s\" type=\"%s\"", key,
	    gretl_arg_type_name(item->type));

    if (item->note != NULL) {
	pprintf(prn, " note=\"%s\"", item->note);
    } 

    if (item->size > 0) {
	pprintf(prn, " size=\"%d\">\n", item->size);
    } else if (item->type == GRETL_TYPE_STRING) {
	pputc(prn, '>');
    } else {
	pputs(prn, ">\n");
    }

    if (item->type == GRETL_TYPE_DOUBLE) {
	double x = *(double *) item->data;
	
	if (na(x)) {
	    pputs(prn, "NA");
	} else {
	    pprintf(prn, "%.15g", x);
	}
    } else if (item->type == GRETL_TYPE_SERIES) {
	double *vals = (double *) item->data;

	for (j=0; j<item->size; j++) {
	    if (na(vals[j])) {
		pputs(prn, "NA ");
	    } else {
		pprintf(prn, "%.15g ", vals[j]);
	    }
	}	    
    } else if (item->type == GRETL_TYPE_STRING) {
	pputs(prn, (char *) item->data);
    } else if (item->type == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = (gretl_matrix *) item->data;

	gretl_xml_put_matrix_to_prn(m, NULL, prn);
    } else {
	fprintf(stderr, "bundle -> XML: skipping %s\n", key);
    }

    pputs(prn, "</bundled-item>\n");    
}

void xml_put_bundle (gretl_bundle *b, const char *name, PRN *prn)
{
    if (b->creator != NULL && *b->creator != '\0') {
	pprintf(prn, "<gretl-bundle name=\"%s\" creator=\"%s\">\n", 
		name, b->creator);
    } else {
	pprintf(prn, "<gretl-bundle name=\"%s\">\n", name);
    }
    g_hash_table_foreach(b->ht, xml_put_bundled_item, prn);
    pputs(prn, "</gretl-bundle>\n"); 
}

static int load_bundled_items (gretl_bundle *b, xmlNodePtr cur, xmlDocPtr doc)
{
    char *key, *typename;
    int err = 0;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "bundled-item")) {
	    key = (char *) xmlGetProp(cur, (XUC) "key");
	    typename = (char *) xmlGetProp(cur, (XUC) "type");
	    if (key == NULL || typename == NULL) {
		err = E_DATA;
	    } else {
		int type = gretl_type_from_string(typename);
		int size = 0;

		if (type == GRETL_TYPE_DOUBLE) {
		    double x;

		    if (!gretl_xml_node_get_double(cur, doc, &x)) {
			err = E_DATA;
		    } else {
			err = gretl_bundle_set_data(b, key, &x, type, size);
		    }
		} else if (type == GRETL_TYPE_STRING) {
		    char *s;

		    if (!gretl_xml_node_get_trimmed_string(cur, doc, &s)) {
			err = E_DATA;
		    } else {
			err = gretl_bundle_set_data(b, key, s, type, size);
			free(s);
		    }
		} else if (type == GRETL_TYPE_SERIES) {
		    double *xvec = gretl_xml_get_double_array(cur, doc, &size, &err);

		    if (!err) {
			err = gretl_bundle_set_data(b, key, xvec, type, size);
			free(xvec);
		    }
		} else if (type == GRETL_TYPE_MATRIX) {
		    xmlNodePtr child = cur->xmlChildrenNode;
		    gretl_matrix *m;

		    if (child == NULL) {
			err = E_DATA;
		    } else {
			m = gretl_xml_get_matrix(child, doc, &err);
			if (!err) {
			    err = gretl_bundle_set_data(b, key, m, type, size);
			    gretl_matrix_free(m);
			}
		    }
		} else {
		    fprintf(stderr, "bundle: ignoring unhandled type %d\n", type);
		}

		if (!err) {
		    char *note = (char *) xmlGetProp(cur, (XUC) "note");

		    if (note != NULL) {
			gretl_bundle_set_note(b, key, note);
			xmlFree(note);
		    }
		}

		xmlFree(key);
		xmlFree(typename);
	    }
	}
	cur = cur->next;
    }

    return err;
}

/* For internal use only, called from uservar.c: @p1 should be of
   type xmlNodePtr and @p2 should be an xmlDocPtr. We suppress the
   actual pointer types in the prototype so that it's possible for a
   module to include gretl_bundle.h without including the full libxml
   headers.  
*/

int load_bundle_from_xml (void *p1, void *p2, const char *name,
			  const char *creator)
{
    xmlNodePtr node = (xmlNodePtr) p1;
    xmlDocPtr doc = (xmlDocPtr) p2;
    xmlNodePtr cur = node->xmlChildrenNode;
    gretl_bundle *b;
    int err = 0;

    b = gretl_bundle_new();
    if (b == NULL) {
	return E_ALLOC;
    }

    fprintf(stderr, "load_bundle_from_xml: '%s'\n", name);

    if (creator != NULL && *creator != '\0') {
	b->creator = gretl_strdup(creator);
    }

    err = load_bundled_items(b, cur, doc);

    if (!err) {
	err = user_var_add(name, GRETL_TYPE_BUNDLE, b);
	fprintf(stderr, "gretl_bundle_push: err = %d\n", err);
    } else {
	gretl_bundle_destroy(b);
	fprintf(stderr, "bundle is broken (err = %d)\n", err);
    }

    return err;
}

int gretl_bundle_write_as_xml (gretl_bundle *b, const char *fname,
			       int to_dotdir)
{
    char fullname[FILENAME_MAX];
    PRN *prn;
    int err = 0;

    if (to_dotdir) {
	build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
	strcpy(fullname, fname);
	gretl_maybe_switch_dir(fullname);
    }

    prn = gretl_print_new_with_filename(fullname, &err);

    if (!err) {
	pputs(prn, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	xml_put_bundle(b, "temp", prn);
	gretl_print_destroy(prn);
    }

    return err;
}

gretl_bundle *gretl_bundle_read_from_xml (const char *fname, int import,
					  int *err)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    char fullname[FILENAME_MAX];
    gretl_bundle *b;

    b = gretl_bundle_new();
    if (b == NULL) {
	*err = E_ALLOC;
	return NULL;
    }    

    if (import) {
	build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
	strcpy(fullname, fname);
    }

    *err = gretl_xml_open_doc_root(fullname, "gretl-bundle", &doc, &cur);

    if (!*err) {
	gretl_push_c_numeric_locale();
	cur = cur->xmlChildrenNode;
	*err = load_bundled_items(b, cur, doc);
	gretl_pop_c_numeric_locale();
	xmlFreeDoc(doc);
    }

    if (*err) {
	gretl_bundle_destroy(b);
	b = NULL;
    }

    return b;
}

gretl_bundle *gretl_bundle_read_from_buffer (const char *buf, int len,
					     int *err)
{
    xmlDocPtr doc = NULL;
    gretl_bundle *b;

    b = gretl_bundle_new();
    if (b == NULL) {
	*err = E_ALLOC;
	return NULL;
    }    

    xmlKeepBlanksDefault(0);
    doc = xmlParseMemory(buf, len);

    if (doc == NULL) {
	gretl_errmsg_set(_("xmlParseMemory failed"));
	*err = 1;
    } else {
	xmlNodePtr cur = xmlDocGetRootElement(doc);

	if (cur == NULL) {
	    gretl_errmsg_set(_("xmlDocGetRootElement failed"));
	    *err = 1;
	} else {
	    gretl_push_c_numeric_locale();
	    cur = cur->xmlChildrenNode;
	    *err = load_bundled_items(b, cur, doc);
	    gretl_pop_c_numeric_locale();
	}
	xmlFreeDoc(doc);
    }

    if (*err) {
	gretl_bundle_destroy(b);
	b = NULL;
    }

    return b;
}

gretl_bundle *get_sysinfo_bundle (int *err)
{
    if (sysinfo_bundle == NULL) {
	gretl_bundle *b = gretl_bundle_new();

	if (b == NULL) {
	    *err = E_ALLOC;
	} else {
	    int ival = 0;

#if HAVE_MPI
	    ival = check_for_mpiexec();
#endif	    
	    gretl_bundle_set_scalar(b, "mpi", (double) ival);
	    ival = gretl_max_mpi_processes();
	    gretl_bundle_set_scalar(b, "mpimax", (double) ival);
	    ival = gretl_n_processors();
	    gretl_bundle_set_scalar(b, "nproc", (double) ival);
	    ival = 0;
#ifdef _OPENMP
	    ival = 1;
#endif
	    gretl_bundle_set_scalar(b, "omp", (double) ival);
	    ival = sizeof(void*) == 8 ? 64 : 32;
	    gretl_bundle_set_scalar(b, "wordlen", (double) ival);
	    gretl_bundle_set_scalar(b, "omp_num_threads", get_omp_n_threads());
#if defined(G_OS_WIN32)
	    gretl_bundle_set_string(b, "os", "windows");
#elif defined(OS_OSX) 
	    gretl_bundle_set_string(b, "os", "osx");
#elif defined(linux)
	    gretl_bundle_set_string(b, "os", "linux");
#else
	    gretl_bundle_set_string(b, "os", "other");
#endif
	    gretl_bundle_set_string(b, "hostname", g_get_host_name());
	}
	sysinfo_bundle = b;
    }

    return sysinfo_bundle;
}
