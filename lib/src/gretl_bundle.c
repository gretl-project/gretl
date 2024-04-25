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
#include "uservar.h"
#include "gretl_mt.h"
#include "gretl_xml.h"
#include "gretl_foreign.h"
#include "gretl_typemap.h"
#include "gen_public.h"
#include "kalman.h"
#include "var.h"
#include "system.h"
#include "matrix_extra.h"
#include "gretl_array.h"
#include "libset.h"
#include "build.h"
#include "gretl_bundle.h"

#define BDEBUG 0

/**
 * gretl_bundle:
 *
 * An opaque type; use the relevant accessor functions.
 */

struct gretl_bundle_ {
    BundleType type; /* see enum in gretl_bundle.h */
    GHashTable *ht;  /* holds key/value pairs */
    char *creator;   /* name of function that built the bundle */
    void *data;      /* holds pointer to struct for some uses */
};

static gretl_bundle *sysinfo_bundle;

static int real_bundle_set_data (gretl_bundle *b, const char *key,
                                 void *ptr, GretlType type,
                                 int size, int copy,
                                 const char *note);

int gretl_bundle_is_stacked (gretl_bundle *b)
{
    user_var *u = get_user_var_by_data(b);

    return u != NULL;
}

/* gets the number of keys in the bundle's hash table */

int gretl_bundle_get_n_keys (gretl_bundle *b)
{
    int n_items = 0;

    if (b != NULL && b->ht != NULL) {
        n_items = g_hash_table_size(b->ht);
    }

    return n_items;
}

/* gets total number of members including any "special"
   contents outside of the hash table */

int gretl_bundle_get_n_members (gretl_bundle *b)
{
    int nmemb = 0;

    if (b != NULL) {
        if (b->type == BUNDLE_KALMAN) {
            nmemb += kalman_bundle_n_members(b);
        }
        if (b->ht != NULL) {
            nmemb += g_hash_table_size(b->ht);
        }
    }

    return nmemb;
}

static void maybe_append_list (gpointer key, gpointer value, gpointer p)
{
    bundled_item *item = (bundled_item *) value;
    GList **plist = (GList **) p;

    if (item->type == GRETL_TYPE_LIST) {
        *plist = g_list_append(*plist, item->data);
    }
}

/* If bundle @b contains any members of type list, compose
   and return a GList that contains all of them.
*/

GList *gretl_bundle_get_lists (gretl_bundle *b)
{
    GList *list = NULL;

    g_hash_table_foreach(b->ht, maybe_append_list, &list);

    return list;
}

int gretl_bundle_has_content (gretl_bundle *b)
{
    int ret = 0;

    if (b != NULL && b->ht != NULL &&
        (b->type == BUNDLE_KALMAN ||
         g_hash_table_size(b->ht) > 0)) {
        ret = 1;
    }

    return ret;
}

int type_can_be_bundled (GretlType type)
{
    if (type == GRETL_TYPE_INT ||
        type == GRETL_TYPE_UNSIGNED ||
        type == GRETL_TYPE_BOOL) {
        type = GRETL_TYPE_DOUBLE;
    }

    return (type == GRETL_TYPE_DOUBLE ||
            type == GRETL_TYPE_STRING ||
            type == GRETL_TYPE_MATRIX ||
            type == GRETL_TYPE_SERIES ||
            type == GRETL_TYPE_BUNDLE ||
            type == GRETL_TYPE_ARRAY ||
            type == GRETL_TYPE_LIST);
}

#define bundled_scalar(t) (t == GRETL_TYPE_DOUBLE ||    \
                           t == GRETL_TYPE_INT ||       \
                           t == GRETL_TYPE_UNSIGNED)

static int bundled_item_copy_in_data (bundled_item *item,
                                      void *ptr, int size)
{
    int err = 0;

    switch (item->type) {
    case GRETL_TYPE_DOUBLE:
        item->data = malloc(sizeof(double));
        if (item->data != NULL) {
            double *dp = item->data;

            *dp = *(double *) ptr;
        }
        break;
    case GRETL_TYPE_INT:
        item->data = malloc(sizeof(int));
        if (item->data != NULL) {
            int *ip = item->data;

            *ip = *(int *) ptr;
        }
        break;
    case GRETL_TYPE_UNSIGNED:
        item->data = malloc(sizeof(guint32));
        if (item->data != NULL) {
            guint32 *up = item->data;

            *up = *(guint32 *) ptr;
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
    case GRETL_TYPE_LIST:
        item->data = gretl_list_copy((const int *) ptr);
        break;
    case GRETL_TYPE_BUNDLE:
        item->data = gretl_bundle_copy((gretl_bundle *) ptr, &err);
        break;
    case GRETL_TYPE_ARRAY:
        item->data = gretl_array_copy((gretl_array *) ptr, &err);
        break;
    default:
        err = E_TYPES;
        break;
    }

    if (!err && item->data == NULL) {
        err = E_ALLOC;
    }

    return err;
}

/* allocate and fill out a 'value' (type plus data pointer) that will
   be inserted into a bundle's hash table */

static bundled_item *bundled_item_new (GretlType type, void *ptr,
                                       int size, int copy,
                                       const char *note,
                                       int *err)
{
    bundled_item *item = malloc(sizeof *item);

    if (item == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    item->type = type;
    item->size = 0;
    item->key = NULL;
    item->note = NULL;

    if (!copy && !bundled_scalar(item->type)) {
        item->data = ptr;
        if (item->type == GRETL_TYPE_SERIES) {
            item->size = size;
        }
    } else {
        *err = bundled_item_copy_in_data(item, ptr, size);
        if (*err) {
            free(item);
            item = NULL;
        }
    }

    if (item != NULL && note != NULL) {
        item->note = gretl_strdup(note);
    }

    return item;
}

static void bundled_item_free_data (GretlType type, void *data)
{
    if (type == GRETL_TYPE_MATRIX) {
        gretl_matrix_free((gretl_matrix *) data);
    } else if (type == GRETL_TYPE_BUNDLE) {
        gretl_bundle_destroy((gretl_bundle *) data);
    } else if (type == GRETL_TYPE_ARRAY) {
        gretl_array_destroy((gretl_array *) data);
    } else {
        free(data);
    }
}

/* Note: we come here only if the replacement type is the
   same as the original type, apart from the case of
   inter-conversion of the various scalar types.
*/

static int bundled_item_replace_data (bundled_item *item,
                                      void *ptr, GretlType src_t,
                                      int size, int copy)
{
    int err = 0;

    if (ptr == item->data) {
        return 0;
    }

    if (bundled_scalar(item->type)) {
        /* storage is already available */
        if (item->type == GRETL_TYPE_DOUBLE) {
            double *dp = item->data;

            if (src_t == GRETL_TYPE_DOUBLE) {
                *dp = *(double *) ptr;
            } else if (src_t == GRETL_TYPE_INT) {
                *dp = (double) *(int *) ptr;
            } else {
                *dp = (double) *(guint32 *) ptr;
            }
        } else if (item->type == GRETL_TYPE_INT) {
            int *ip = item->data;

            if (src_t == GRETL_TYPE_INT) {
                *ip = *(int *) ptr;
            } else if (src_t == GRETL_TYPE_DOUBLE) {
                *ip = gretl_int_from_double(*(double *) ptr, &err);
            } else {
                *ip = (int) *(guint32 *) ptr;
            }
        } else if (item->type == GRETL_TYPE_UNSIGNED) {
            guint32 *up = item->data;

            if (src_t == GRETL_TYPE_UNSIGNED) {
                *up = *(guint32 *) ptr;
            } else if (src_t == GRETL_TYPE_DOUBLE) {
                *up = gretl_unsigned_from_double(*(double *) ptr, &err);
            } else {
                *up = (guint32) *(int *) ptr;
            }
        }
    } else {
        /* free then copy or donate */
        bundled_item_free_data(item->type, item->data);
        if (copy) {
            err = bundled_item_copy_in_data(item, ptr, size);
        } else {
            item->data = ptr;
        }
    }

    if (!err && item->data == NULL) {
        err = E_ALLOC;
    }

    if (item->note != NULL) {
        free(item->note);
        item->note = NULL;
    }

    return err;
}

/* callback invoked when a bundle's hash table is destroyed */

static void bundle_item_destroy (gpointer data)
{
    bundled_item *item = data;

#if BDEBUG
    fprintf(stderr, "bundled_item_destroy: %s\t'%s'\tdata %p\n",
            gretl_type_get_name(item->type), item->key, (void *) item->data);
    if (item->type == GRETL_TYPE_STRING) {
        fprintf(stderr, " string: '%s'\n", (char *) item->data);
    }
#endif

    bundled_item_free_data(item->type, item->data);
    g_free(item->key);
    free(item->note);
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
        free(bundle->creator);
        if (bundle->type == BUNDLE_KALMAN) {
            kalman_free(bundle->data);
        }
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
    if (bundle == NULL) {
        return;
    }

    if (bundle->creator != NULL) {
        free(bundle->creator);
        bundle->creator = NULL;
    }

    if (bundle->ht != NULL && g_hash_table_size(bundle->ht) > 0) {
        g_hash_table_destroy(bundle->ht);
        bundle->ht = g_hash_table_new_full(g_str_hash, g_str_equal,
                                           NULL, bundle_item_destroy);
    }

    if (bundle->type == BUNDLE_KALMAN) {
        kalman_free(bundle->data);
        bundle->data = NULL;
        bundle->type = BUNDLE_PLAIN;
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
        b->type = BUNDLE_PLAIN;
        b->ht = g_hash_table_new_full(g_str_hash, g_str_equal,
                                      NULL, bundle_item_destroy);
        b->creator = NULL;
        b->data = NULL;
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

/* Supports gretl_bundle_get_data(), gretl_bundle_get_element()
   and gretl_bundle_get_target.
*/

static void *real_gretl_bundle_get_data (gretl_bundle *bundle,
                                         const char *key,
                                         GretlType *type,
                                         int *size,
                                         int *ownit,
                                         int get_target,
                                         int *err)
{
    void *ret = NULL;
    int reserved = 0;
    int myerr = 0;

    if (bundle == NULL) {
        myerr = E_DATA;
        goto finish;
    }

    if (!get_target && bundle->type == BUNDLE_KALMAN) {
        ret = maybe_retrieve_kalman_element(bundle->data, key,
                                            type, &reserved,
                                            ownit, &myerr);
    }

    if (!myerr && ret == NULL && !reserved) {
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
        } else {
            if (err != NULL) {
                gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
            }
            myerr = E_DATA;
        }
    }

 finish:

    if (err != NULL) {
        *err = myerr;
    }

    return ret;
}

/**
 * gretl_bundle_get_element:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @type: location to receive data type, or NULL.
 * @size: location to receive size of data (= series
 * length for GRETL_TYPE_SERIES, otherwise 0), or NULL.
 * @ownit: location to receive indication that the
 * element in question is newly allocated, not just a
 * pointer to something in the bundle.
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
 * If the content of @ownit is non-zero on return, the caller
 * has ownership of the retuned data; otherwise the return
 * value is the actual data pointer from within the bundle.
 */

void *gretl_bundle_get_element (gretl_bundle *bundle, const char *key,
                                GretlType *type, int *size,
                                int *ownit, int *err)
{
    return real_gretl_bundle_get_data(bundle, key, type,
                                      size, ownit, 0, err);
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
    return real_gretl_bundle_get_data(bundle, key, type,
                                      size, NULL, 0, err);
}

/* As gretl_bundle_get_data() except that the caller is looking
   for a "target" object within @bundle with a view to modifying
   it. Called only from geneval.c, in response to, for example,
   b.X = Y.
*/

void *gretl_bundle_get_target (gretl_bundle *bundle, const char *key,
                               GretlType *type, int *size, int *err)
{
    return real_gretl_bundle_get_data(bundle, key, type,
                                      size, NULL, 1, err);
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

void *gretl_bundle_get_private_data (gretl_bundle *bundle)
{
    return bundle->data;
}

BundleType gretl_bundle_get_type (gretl_bundle *bundle)
{
    return bundle->type;
}

/**
 * gretl_bundle_get_member_type:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err:location to receive error code, or NULL.
 *
 * Returns: the data type associated with @key in the
 * specified @bundle, or 0 on failure. Set the @err
 * argument to NULL if you do not want an error flagged
 * in case there's no such key in @bundle.
 */

GretlType gretl_bundle_get_member_type (gretl_bundle *bundle,
                                        const char *key,
                                        int *err)
{
    GretlType ret = GRETL_TYPE_NONE;
    int reserved = 0;
    int myerr = 0;

    if (bundle == NULL) {
        myerr = E_DATA;
    } else if (bundle->type == BUNDLE_KALMAN) {
        maybe_retrieve_kalman_element(bundle->data, key,
                                      &ret, &reserved,
                                      NULL, &myerr);
    }

    if (!myerr && ret == GRETL_TYPE_NONE && !reserved) {
        gpointer p = g_hash_table_lookup(bundle->ht, key);

        if (p != NULL) {
            bundled_item *item = p;

            ret = item->type;
        } else if (err != NULL) {
            gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
            myerr = E_DATA;
        }
    }

    if (err != NULL) {
        *err = myerr;
    }

    return ret;
}

/**
 * gretl_bundle_has_key:
 * @bundle: bundle to access.
 * @key: name of key to test.
 *
 * Returns: 1 if there is an item under the given @key in the
 * specified @bundle, 0 otherwise.
 */

int gretl_bundle_has_key (gretl_bundle *bundle,
                          const char *key)
{
    int ret = 0;

    if (bundle != NULL && key != NULL) {
        gpointer p = g_hash_table_lookup(bundle->ht, key);

        ret = (p != NULL);
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
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);

    if (ptr != NULL && type != GRETL_TYPE_MATRIX) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        m = (gretl_matrix *) ptr;
    }

    if (err != NULL) {
        *err = myerr;
    }

    return m;
}

/**
 * gretl_bundle_get_array:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the array associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

gretl_array *gretl_bundle_get_array (gretl_bundle *bundle,
                                     const char *key,
                                     int *err)
{
    gretl_array *a = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_ARRAY) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        a = (gretl_array *) ptr;
    }

    if (err != NULL) {
        *err = myerr;
    }

    return a;
}

/**
 * gretl_bundle_get_bundle:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the bundle associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

gretl_bundle *gretl_bundle_get_bundle (gretl_bundle *bundle,
                                       const char *key,
                                       int *err)
{
    gretl_bundle *ret = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_BUNDLE) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        ret = (gretl_bundle *) ptr;
    }

    if (err != NULL) {
        *err = myerr;
    }

    return ret;
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
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, n, err);
    if (ptr != NULL && type != GRETL_TYPE_SERIES) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        x = (double *) ptr;
    }

    if (err != NULL) {
        *err = myerr;
    }

    return x;
}

/**
 * gretl_bundle_get_list:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the list associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

int *gretl_bundle_get_list (gretl_bundle *bundle,
                            const char *key,
                            int *err)
{
    int *list = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_LIST) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        list = (int *) ptr;
    }

    if (err != NULL) {
        *err = myerr;
    }

    return list;
}

/**
 * gretl_bundle_get_scalar:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code, or NULL.
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
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);

    if (ptr == NULL) {
        myerr = E_DATA;
    } else if (type != GRETL_TYPE_DOUBLE &&
               type != GRETL_TYPE_INT &&
               type != GRETL_TYPE_UNSIGNED &&
               type != GRETL_TYPE_MATRIX) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        if (type == GRETL_TYPE_DOUBLE) {
            double *px = (double *) ptr;

            x = *px;
        } else if (type == GRETL_TYPE_INT) {
            int *pi = (int *) ptr;

            x = (double) *pi;
        } else if (type == GRETL_TYPE_UNSIGNED) {
            guint32 *pu = (guint32 *) ptr;

            x = (double) *pu;
        } else {
            /* must be a matrix */
            gretl_matrix *m = ptr;

            if (gretl_matrix_is_scalar(m)) {
                x = m->val[0];
            } else {
                myerr = E_TYPES;
            }
        }
    }

    if (err != NULL) {
        *err = myerr;
    }

    return x;
}

/**
 * gretl_bundle_get_int:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the integer value associated with @key in the
 * specified @bundle, if any; otherwise 0.
 */

int gretl_bundle_get_int (gretl_bundle *bundle,
                          const char *key,
                          int *err)
{
    int i = 0;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);

    if (ptr != NULL) {
        if (type == GRETL_TYPE_INT) {
            int *pi = (int *) ptr;

            i = *pi;
        } else if (type == GRETL_TYPE_UNSIGNED) {
            guint32 u, *pu = (guint32 *) ptr;

            u = *pu;
            if (u <= INT_MAX) {
                i = (int) u;
            } else {
                myerr = E_TYPES;
            }
        } else if (type == GRETL_TYPE_DOUBLE) {
            double x, *px = (double *) ptr;

            x = *px;
            if (x >= INT_MIN && x <= INT_MAX && x == floor(x)) {
                i = (int) x;
            } else {
                myerr = E_TYPES;
            }
        } else {
            myerr = E_TYPES;
        }
    }

    if (err != NULL) {
        *err = myerr;
    }

    return i;
}

/**
 * gretl_bundle_get_int_deflt:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @deflt: default integer value.
 *
 * Returns: the integer value associated with @key in the
 * specified @bundle, if any; otherwise @deflt.
 */

int gretl_bundle_get_int_deflt (gretl_bundle *bundle,
                                const char *key,
                                int deflt)
{
    int val = deflt;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, NULL);

    if (ptr != NULL) {
        if (type == GRETL_TYPE_INT) {
            int *pi = (int *) ptr;

            val = *pi;
        } else if (type == GRETL_TYPE_UNSIGNED) {
            guint32 *pu = (guint32 *) ptr;

            val = (int) *pu;
        } else if (type == GRETL_TYPE_DOUBLE) {
            double *px = (double *) ptr;

            val = (int) *px;
        }
    }

    return val;
}

/**
 * gretl_bundle_get_bool:
 * @bundle: bundle to access.
 * @key: name of key to access, if present.
 * @deflt: default boolean value.
 *
 * Returns: the boolean value associated with @key in the
 * specified @bundle, if any; otherwise @deflt.
 */

int gretl_bundle_get_bool (gretl_bundle *bundle,
                           const char *key,
                           int deflt)
{
    int val = deflt;
    GretlType type;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, NULL);

    if (ptr != NULL) {
        if (type == GRETL_TYPE_INT) {
            int *pi = (int *) ptr;

            val = (*pi != 0);
        } else if (type == GRETL_TYPE_UNSIGNED) {
            guint32 *pu = (guint32 *) ptr;

            val = (*pu != 0);
        } else if (type == GRETL_TYPE_DOUBLE) {
            double *px = (double *) ptr;

            val = (*px != 0);
        }
    }

    return val;
}

/**
 * gretl_bundle_get_unsigned:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the unsigned integer value associated with @key in the
 * specified @bundle, if any; otherwise 0.
 */

guint32 gretl_bundle_get_unsigned (gretl_bundle *bundle,
                                   const char *key,
                                   int *err)
{
    guint32 u = 0;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);

    if (ptr != NULL) {
        if (type == GRETL_TYPE_UNSIGNED) {
            guint32 *pu = (guint32 *) ptr;

            u = *pu;
        } else if (type == GRETL_TYPE_INT) {
            int i, *pi = (int *) ptr;

            i = *pi;
            if (i >= 0) {
                u = (guint32) i;
            } else {
                myerr = E_TYPES;
            }
        } else if (type == GRETL_TYPE_DOUBLE) {
            double x, *px = (double *) ptr;

            x = *px;
            if (x >= 0 && x <= UINT_MAX && x == floor(x)) {
                u = (guint32) x;
            } else {
                myerr = E_TYPES;
            }
        } else {
            myerr = E_TYPES;
        }
    }

    if (err != NULL) {
        *err = myerr;
    }

    return u;
}

/**
 * gretl_bundle_get_string:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @err: location to receive error code.
 *
 * Returns: the string value associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const char *gretl_bundle_get_string (gretl_bundle *bundle,
                                     const char *key,
                                     int *err)
{
    const char *ret = NULL;
    GretlType type;
    void *ptr;
    int myerr = 0;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, err);
    if (ptr != NULL && type != GRETL_TYPE_STRING) {
        myerr = E_TYPES;
    }

    if (ptr != NULL && !myerr) {
        ret = (const char *) ptr;
    }

    if (err != NULL) {
        *err = myerr;
    }

    return ret;
}

/**
 * gretl_bundle_get_strings:
 * @bundle: bundle to access.
 * @key: name of key to access.
 * @ns: location to receive the number of strings, if any.
 *
 * Returns: the array of strings associated with @key in the
 * specified @bundle, if any; otherwise NULL.
 */

const char **gretl_bundle_get_strings (gretl_bundle *bundle,
                                       const char *key,
                                       int *ns)
{
    const char **ret = NULL;
    GretlType type = 0;
    void *ptr;

    ptr = gretl_bundle_get_data(bundle, key, &type, NULL, NULL);
    if (ptr != NULL && type == GRETL_TYPE_ARRAY) {
        gretl_array *a = ptr;

        ret = (const char **) gretl_array_get_strings(a, ns);
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

/**
 * gretl_bundles_swap_content:
 * @b1: first bundle.
 * @b2: second bundle.
 *
 * Exchanges the entire contents of the two bundles.
 *
 * Returns: 0 on success, error code on failure.
 */

int gretl_bundles_swap_content (gretl_bundle *b1, gretl_bundle *b2)
{
    if (b1 == NULL || b1->type != BUNDLE_PLAIN ||
        b2 == NULL || b2->type != BUNDLE_PLAIN) {
        return E_DATA;
    } else {
        void *tmp = b1->ht;

        b1->ht = b2->ht;
        b2->ht = tmp;
        return 0;
    }
}

static int real_bundle_set_data (gretl_bundle *b, const char *key,
                                 void *ptr, GretlType type,
                                 int size, int copy,
                                 const char *note)
{
    int err = 0, done = 0;

    if (key == NULL || key[0] == '\0') {
        gretl_errmsg_sprintf("real_bundle_set_data: missing key string");
        return E_DATA;
    }

    if (b->type == BUNDLE_KALMAN) {
        done = maybe_set_kalman_element(b->data, key,
                                        ptr, type, copy,
                                        &err);
    }

    if (!done && !err) {
        bundled_item *item = g_hash_table_lookup(b->ht, key);
        int replace = 0;

        if (item != NULL) {
            replace = 1;
            if (item->type == type) {
                /* we can take a shortcut */
                return bundled_item_replace_data(item, ptr, type, size, copy);
            }
        }

        item = bundled_item_new(type, ptr, size, copy, note, &err);

        if (!err) {
            item->key = g_strdup(key);
            if (replace) {
                g_hash_table_replace(b->ht, item->key, item);
            } else {
                g_hash_table_insert(b->ht, item->key, item);
            }
        }
    }

    return err;
}

/**
 * gretl_bundle_donate_data:
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
 * and destroyed. The value of @ptr is transcribed into the
 * bundle, which therefore "takes ownership" of the data;
 * compare gretl_bundle_set_data().
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_donate_data (gretl_bundle *bundle, const char *key,
                              void *ptr, GretlType type, int size)
{
    if (key != NULL && ptr == NULL) {
        gretl_errmsg_sprintf("'%s': got NULL data value", key);
        return E_DATA;
    }

    return real_bundle_set_data(bundle, key, ptr, type, size, 0, NULL);
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
 * and destroyed. The content of @ptr is copied into the
 * bundle; compare gretl_bundle_donate_data().
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
        err = real_bundle_set_data(bundle, key, ptr, type,
                                   size, 1, NULL);
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
 * gretl_bundle_set_int:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the integer value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_int (gretl_bundle *bundle, const char *key,
                          int val)
{
    return gretl_bundle_set_data(bundle, key, &val,
                                 GRETL_TYPE_INT, 0);
}

/**
 * gretl_bundle_set_unsigned:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @val: the unsigned integer value to set.
 *
 * Sets @val as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_unsigned (gretl_bundle *bundle, const char *key,
                               guint32 val)
{
    return gretl_bundle_set_data(bundle, key, &val,
                                 GRETL_TYPE_UNSIGNED, 0);
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
 * gretl_bundle_set_list:
 * @bundle: target bundle.
 * @key: name of key to create or replace.
 * @list: gretl list.
 *
 * Sets @list as a member of @bundle under the name @key.
 * If @key is already present in the bundle the original
 * value is replaced and destroyed.
 *
 * Returns: 0 on success, error code on error.
 */

int gretl_bundle_set_list (gretl_bundle *bundle, const char *key,
                           const int *list)
{
    return gretl_bundle_set_data(bundle, key, (void *) list,
                                 GRETL_TYPE_LIST, 0);
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
    int done = 0;
    int err = 0;

    if (bundle == NULL) {
        return E_DATA;
    }

    if (bundle->type == BUNDLE_KALMAN) {
        done = maybe_delete_kalman_element(bundle->data, key, &err);
    }

    if (!done && !err) {
        done = g_hash_table_remove(bundle->ht, key);
        if (!done) {
            err = E_DATA;
        }
    }

    return err;
}

/**
 * gretl_bundle_rekey_data:
 * @bundle: target bundle.
 * @oldkey: name of key to access.
 * @newkey: revised key.
 *
 * If @bundle contains an item under the key @oldkey, changes
 * the string by which the item is identified to @newkey.
 *
 * Returns: 0 on success, error code on failure.
 */

int gretl_bundle_rekey_data (gretl_bundle *bundle,
                             const char *oldkey,
                             const char *newkey)
{
    bundled_item *item = NULL;
    int err = 0;

    if (strcmp(oldkey, newkey) == 0) {
        return 0;
    }

    if (bundle != NULL) {
        item = g_hash_table_lookup(bundle->ht, oldkey);
    }

    if (item == NULL) {
        err = E_DATA;
    } else {
        g_hash_table_steal(bundle->ht, oldkey);
        g_free(item->key);
        item->key = g_strdup(newkey);
        g_hash_table_insert(bundle->ht, item->key, item);
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
                             item->size, 1, item->note);
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
                         item->size, 1, item->note);
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
    gretl_bundle *b = NULL;

    if (bundle2->type == BUNDLE_KALMAN) {
        gretl_errmsg_set(_("bundle union: the right-hand operand cannot "
                           "be a kalman bundle"));
        *err = E_DATA;
    } else {
        b = gretl_bundle_copy(bundle1, err);
    }

    if (!*err) {
        g_hash_table_foreach(bundle2->ht, copy_new_bundled_item, b);
    }

    return b;
}

/* argument-bundle checking apparatus */

struct bchecker_ {
    gretl_bundle *b;
    char **bounds_keys;
    gretl_matrix *bounds_mat;
    int n_bounds;
    char **ignore;
    int n_ignores;
    int *ret;
    int *err;
    PRN *prn;
};

typedef struct bchecker_ bchecker;

static void bchecker_free (bchecker *bc)
{
    if (bc->n_bounds > 0) {
        gretl_matrix_free(bc->bounds_mat);
        strings_array_free(bc->bounds_keys, bc->n_bounds);
    }
    free(bc);
}

static int unpack_bounds_info (bchecker *bc, gretl_bundle *bounds);

static bchecker *bchecker_new (gretl_bundle *b, int *ret, int *err,
                               gretl_array *ignore, PRN *prn)
{
    bchecker *bc = calloc(1, sizeof *bc);

    if (bc == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    bc->b = b;
    bc->ret = ret;
    bc->err = err;
    bc->prn = prn;

    if (gretl_bundle_has_key(b, "bounds")) {
        gretl_bundle *bb;

        bb = gretl_bundle_get_bundle(b, "bounds", err);
        if (!*err) {
            *err = unpack_bounds_info(bc, bb);
        }
    }

    if (!*err && ignore != NULL) {
        char **S = NULL;
        int ns = 0;

        /* note: borrowing here */
        S = gretl_array_get_strings(ignore, &ns);
        if (S != NULL) {
            bc->ignore = S;
            bc->n_ignores = ns;
        }
    }

    return bc;
}

/* Check for mismatched types in @targ vs @src, allowing for
   convertibility between scalar types.
*/

static int bundled_types_mismatch (bundled_item *targ,
                                   bundled_item *src)
{
    int err = 0;

    if (src->type == targ->type) {
        return 0; /* OK */
    } else if (gretl_is_scalar_type(targ->type) &&
               gretl_is_scalar_type(src->type)) {
        /* may be inter-convertible: the main potential
           problem, diagnosed here, lies in conversion
           from double to integer types
        */
        if (src->type == GRETL_TYPE_DOUBLE) {
            double *px = (double *) src->data;

            if (targ->type == GRETL_TYPE_BOOL) {
                err = na(*px);
            } else if (targ->type == GRETL_TYPE_INT) {
                gretl_int_from_double(*px, &err);
            } else if (targ->type == GRETL_TYPE_UNSIGNED) {
                gretl_unsigned_from_double(*px, &err);
            }
        }
    } else {
        /* other mismatches will not work */
        err = 1;
    }

    return err;
}

static int should_ignore (const char *key,
                          bchecker *bc)
{
    int i;

    for (i=0; i<bc->n_ignores; i++) {
        if (!strcmp(key, bc->ignore[i])) {
            return 1;
        }
    }

    return 0;
}

static double bundled_item_get_scalar (bundled_item *item)
{
    if (item->type == GRETL_TYPE_DOUBLE) {
        return *(double *) item->data;
    } else if (item->type == GRETL_TYPE_INT) {
        return (double) *(int *) item->data;
    } else {
        return NADBL;
    }
}

/* check @src, a presumably scalar bundled input value,
   against any relevant bounds specification in @bc
*/

static int bcheck_bounds (const char *key,
                          bundled_item *src,
                          bchecker *bc)
{
    double lo, hi, val;
    const char *s;
    int i;

    for (i=0; i<bc->n_bounds; i++) {
        s = bc->bounds_keys[i];
        if (!strcmp(s, src->key)) {
            /* found an associated limits item */
            if (!gretl_is_scalar_type(src->type)) {
                pprintf(bc->prn, _("%s: should be scalar\n"), src->key);
                return 1;
            }
            val = bundled_item_get_scalar(src);
            if (na(val)) {
                pprintf(bc->prn, _("%s: invalid (missing) value\n"), src->key);
                return 1;
            }
            lo = gretl_matrix_get(bc->bounds_mat, i, 0);
            hi = gretl_matrix_get(bc->bounds_mat, i, 1);
            if (!na(lo) && val < lo) {
                pprintf(bc->prn, _("%s: value %g is out of bounds (lower limit %g)\n"),
                        src->key, val, lo);
                return 1;
            } else if (!na(hi) && val > hi) {
                pprintf(bc->prn, _("%s: value %g is out of bounds (upper limit %g)\n"),
                        src->key, val, hi);
                return 1;
            }
            break;
        }
    }

    return 0;
}

/* Given a sub-bundle named 'bounds', unpack it to form an
   array of n >= 1 key-strings and an associated n x 2 matrix
   holding lower and upper bounds.
*/

static int unpack_bounds_info (bchecker *bc, gretl_bundle *bounds)
{
    gretl_matrix *bi, *bb = NULL;
    char **S = NULL;
    int i, ns = 0;
    int err = 0;

    S = gretl_bundle_get_keys_raw(bounds, &ns);
    if (S == NULL) {
        err = 1;
    } else {
        bb = gretl_matrix_alloc(ns, 2);
        if (bb == NULL) {
            err = 1;
        }
    }

    for (i=0; i<ns && !err; i++) {
        bi = gretl_bundle_get_matrix(bounds, S[i], &err);
        if (!err && gretl_vector_get_length(bi) != 2) {
            err = E_INVARG;
        }
        if (!err) {
            gretl_matrix_set(bb, i, 0, bi->val[0]);
            gretl_matrix_set(bb, i, 1, bi->val[1]);
        }
    }

    if (!err) {
        bc->bounds_keys = S;
        bc->bounds_mat = bb;
        bc->n_bounds = ns;
    }

    return err;
}

static void check_bundled_item (gpointer key, gpointer value, gpointer p)
{
    bundled_item *targ, *src = (bundled_item *) value;
    bchecker *bc = (bchecker *) p;

    if (*bc->ret || *bc->err) {
        /* don't waste time if we already hit an error */
        return;
    }

    if (bc->ignore != NULL && should_ignore(key, bc)) {
        return;
    }

    if (bc->n_bounds > 0) {
        *bc->ret = bcheck_bounds(key, src, bc);
        if (*bc->ret) {
            return;
        }
    }

    /* look up @key (from input) in the template bundle */
    targ = g_hash_table_lookup(bc->b->ht, (const char *) key);

    if (targ == NULL) {
        /* extraneous key in input */
        pprintf(bc->prn, "bcheck: unrecognized key '%s'\n", key);
        *bc->ret = 2;
    } else if (targ->type == GRETL_TYPE_MATRIX &&
               src->type == GRETL_TYPE_DOUBLE) {
        double x = *(double *) src->data;
        gretl_matrix *m = gretl_matrix_from_scalar(x);

        *bc->err = bundled_item_replace_data(targ, m, GRETL_TYPE_MATRIX,
                                             0, 0);
    } else if (targ->type == GRETL_TYPE_DOUBLE &&
               src->type == GRETL_TYPE_MATRIX) {
        gretl_matrix *m = src->data;

        if (gretl_matrix_is_scalar(m)) {
            *bc->err = bundled_item_replace_data(targ, &m->val[0],
                                                 GRETL_TYPE_DOUBLE,
                                                 0, 0);
        } else {
            *bc->ret = 3;
        }
    } else if (bundled_types_mismatch(targ, src)) {
        *bc->ret = 3;
    } else {
        /* transcribe input value -> template */
        *bc->err = bundled_item_replace_data(targ, src->data, src->type, 0, 1);
        if (*bc->err) {
            pprintf(bc->prn, "bcheck: failed to copy '%s'\n", key);
        }
    }

    if (*bc->ret == 3) {
        pprintf(bc->prn, "bcheck: '%s' should be %s, is %s\n", key,
                gretl_type_get_name(targ->type),
                gretl_type_get_name(src->type));
    }
}

/**
 * gretl_bundle_extract_args:
 * @defaults: bundle containing keys for all supported
 * inputs, both optional and required (if any).
 * @input: bundle supplied by caller.
 * @reqd: array of strings identifying required keys, if any
 * (or NULL).
 * @ignore: array of strings identifying keys that are not
 * associated with options but whose presence should not
 * provoke an error.
 * @prn: pointer to printing struct for display of error
 * messages.
 * @err:location to receive error code.
 *
 * This function checks @input against @defaults. It flags an
 * error (a) if a required element is missing, (b) if @input
 * contains an unrecognized key, or (c) if the type of any element
 * in @input fails to match the type of the corresponding element
 * in @defaults. If there's no error the content of @defaults
 * is updated from @input; that is, default values are replaced by
 * values selected by the caller.
 *
 * Note the the @err pointer receives a non-zero value only
 * if this function fails, which is distinct from the case
 * where it successfully diagnoses an error in @input.
 *
 * Returns: 0 if @input is deemed valid in light of @defaults,
 * non-zero otherwise.
 */

int gretl_bundle_extract_args (gretl_bundle *defaults,
                               gretl_bundle *input,
                               gretl_array *reqd,
                               gretl_array *ignore,
                               PRN *prn, int *err)
{
    int ret = 0;

    if (reqd != NULL && gretl_array_get_type(reqd) != GRETL_TYPE_STRINGS) {
        *err = E_TYPES;
        return *err;
    }

    if (input == NULL) {
        if (reqd == NULL) {
            /* OK, it's a no-op */
            return 0;
        } else {
            /* null @input is not acceptable given @reqd */
            pprintf(prn, "bcheck: required input bundle is missing\n");
            return 1;
        }
    }

    if (reqd != NULL) {
        int i, n = gretl_array_get_length(reqd);
        const char *key;

        for (i=0; i<n; i++) {
            key = gretl_array_get_data(reqd, i);
            if (!gretl_bundle_has_key(input, key)) {
                pprintf(prn, "bcheck: required key '%s' is missing\n", key);
                ret = 1;
                break;
            }
        }
    }

    if (ret == 0) {
        bchecker *bc = bchecker_new(defaults, &ret, err, ignore, prn);

        if (bc != NULL) {
            g_hash_table_foreach(input->ht, check_bundled_item, bc);
            bchecker_free(bc);
        }
    }

    return ret;
}

int gretl_bundle_append (gretl_bundle *bundle1,
                         const gretl_bundle *bundle2)
{
    if (bundle2->type == BUNDLE_KALMAN) {
        gretl_errmsg_set(_("bundle union: the right-hand operand cannot "
                           "be a kalman bundle"));
        return E_DATA;
    } else {
        g_hash_table_foreach(bundle2->ht, copy_new_bundled_item, bundle1);
        return 0;
    }
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
        if (bundle->type == BUNDLE_KALMAN) {
            bcpy = kalman_bundle_copy(bundle, err);
        } else {
            bcpy = gretl_bundle_new();
            if (bcpy == NULL) {
                *err = E_ALLOC;
            }
        }
        if (!*err) {
            g_hash_table_foreach(bundle->ht, copy_bundled_item, bcpy);
        }
    }

    return bcpy;
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
    int prev = 0;
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
        b1 = user_var_steal_value(u);
        if (b1 != NULL) {
            gretl_bundle_destroy(b1);
        }
        prev = 1;
    }

    b1 = gretl_bundle_copy(b0, &err);

    if (!err) {
        if (prev) {
            err = user_var_replace_value(u, b1, GRETL_TYPE_BUNDLE);
        } else {
            err = user_var_add(copyname, GRETL_TYPE_BUNDLE, b1);
        }
    }

    return err;
}

struct b_item_printer {
    PRN *prn;
    int indent;
    int tree;
};

static void print_bundled_item (gpointer value, gpointer p)
{
    bundled_item *item = value;
    GretlType t = item->type;
    const gchar *kstr = item->key;
    struct b_item_printer *bip = p;
    PRN *prn = bip->prn;
    int indent;

    if (bip->tree && t == GRETL_TYPE_BUNDLE) {
        return;
    }

    indent = 2 + 2 * bip->indent;
    bufspace(indent, prn);

    if (t == GRETL_TYPE_DOUBLE) {
        double x = *(double *) item->data;

        if (na(x)) {
            pprintf(prn, "%s = NA", kstr);
        } else {
            pprintf(prn, "%s = %g", kstr, x);
        }
    } else if (t == GRETL_TYPE_INT) {
        int i = *(int *) item->data;

        pprintf(prn, "%s = %d", kstr, i);
    } else if (t == GRETL_TYPE_UNSIGNED) {
        guint32 u = *(guint32 *) item->data;

        pprintf(prn, "%s = %u", kstr, u);
    } else if (t == GRETL_TYPE_STRING) {
        char *s = (char *) item->data;
        int n = strlen(s);

        if (n + indent < 70) {
            if (strchr(s, '"') != NULL) {
                pprintf(prn, "%s = '%s'", kstr, s);
            } else {
                pprintf(prn, "%s = \"%s\"", kstr, s);
            }
        } else {
            pprintf(prn, "%s (%s, %d bytes)", kstr,
                    gretl_type_get_name(t), n);
        }
    } else if (t == GRETL_TYPE_BUNDLE) {
        pprintf(prn, "%s (%s)", kstr, gretl_type_get_name(t));
    } else if (t == GRETL_TYPE_MATRIX) {
        gretl_matrix *m = item->data;

        if (m->rows == 1 && m->cols == 1) {
            pprintf(prn, "%s = %g", kstr, m->val[0]);
        } else {
            pprintf(prn, "%s (%s: %d x %d)", kstr,
                    gretl_type_get_name(t),
                    m->rows, m->cols);
        }
    } else if (t == GRETL_TYPE_SERIES) {
        pprintf(prn, "%s (%s: length %d)", kstr,
                gretl_type_get_name(t), item->size);
    } else if (t == GRETL_TYPE_LIST) {
        pprintf(prn, "%s (%s)", kstr, gretl_type_get_name(t));
    } else if (t == GRETL_TYPE_ARRAY) {
        gretl_array *a = item->data;
        int n = gretl_array_get_length(a);

        t = gretl_array_get_type(a);
        pprintf(prn, "%s = array of %s, length %d", kstr,
                gretl_type_get_name(t), n);
    }

    if (item->note != NULL) {
        pprintf(prn, " %s\n", item->note);
    } else {
        pputc(prn, '\n');
    }
}

static void bundle_header (const char *name,
                           const char *creator,
                           int empty, PRN *prn)
{
    if (name != NULL) {
        if (empty) {
            pprintf(prn, "bundle %s: empty\n", name);
        } else if (creator != NULL) {
            pprintf(prn, "bundle %s, created by %s:\n", name, creator);
        } else {
            pprintf(prn, "bundle %s:\n", name);
        }
    } else {
        if (empty) {
            pputs(prn, "bundle: empty\n");
        } else if (creator != NULL) {
            pprintf(prn, "bundle created by %s:\n", creator);
        } else {
            pputs(prn, "bundle:\n");
        }
    }
}

static int real_bundle_print (gretl_bundle *bundle, int indent,
                              int tree, PRN *prn)
{
    struct b_item_printer bip = {prn, indent, tree};
    GList *L;

    if (bundle == NULL) {
        return E_DATA;
    }

    L = gretl_bundle_get_sorted_items(bundle);

    if (indent > 0) {
        /* child, when printing tree */
        int n_items = g_hash_table_size(bundle->ht);

        if (bundle->type == BUNDLE_PLAIN && n_items == 0) {
            pputs(prn, "empty\n");
        } else if (bundle->type == BUNDLE_KALMAN) {
            print_kalman_bundle_info(bundle->data, prn);
            if (n_items > 0) {
                pputs(prn, "\nOther content\n");
                g_list_foreach(L, print_bundled_item, &bip);
            }
        } else if (n_items > 0) {
            g_list_foreach(L, print_bundled_item, &bip);
        }
    } else {
        int n_items = g_hash_table_size(bundle->ht);
        user_var *u = get_user_var_by_data(bundle);
        const char *name = NULL;

        if (u != NULL) {
            name = user_var_get_name(u);
        }

        if (bundle->type == BUNDLE_PLAIN && n_items == 0) {
            bundle_header(name, NULL, 1, prn);
        } else {
            bundle_header(name, bundle->creator, 0, prn);
            if (bundle->type == BUNDLE_KALMAN) {
                print_kalman_bundle_info(bundle->data, prn);
                if (n_items > 0) {
                    pputs(prn, "\nOther content\n");
                    g_list_foreach(L, print_bundled_item, &bip);
                }
            } else if (n_items > 0) {
                g_list_foreach(L, print_bundled_item, &bip);
            }
        }
    }

    g_list_free(L);

    return 0;
}

static int print_bundle_tree (gretl_bundle *b, int level, PRN *prn)
{
    gretl_array *K = gretl_bundle_get_keys(b, NULL);
    gretl_array *A;
    gretl_bundle *bj;
    char **keys;
    int i, j, nb, n = 0;

    real_bundle_print(b, level, 1, prn);

    keys = gretl_array_get_strings(K, &n);
    if (n == 0) {
        bufspace(2*(level+1), prn);
        pputs(prn, "empty\n");
    }

    for (i=0; i<n; i++) {
        const char *key = keys[i];
        bundled_item *item = g_hash_table_lookup(b->ht, key);

        if (item->type == GRETL_TYPE_BUNDLE) {
            int nmemb = gretl_bundle_get_n_members(item->data);

            bufspace(2*(level+1), prn);
            if (nmemb == 0) {
                pprintf(prn, "%s = bundle (empty)\n", key);
            } else {
                pprintf(prn, "%s = bundle (%d members)\n", key, nmemb);
                print_bundle_tree((gretl_bundle *) item->data, level + 1, prn);
            }
            continue;
        } else if (item->type == GRETL_TYPE_ARRAY) {
            GretlType atype;

            A = (gretl_array *) item->data;
            atype = gretl_array_get_content_type(A);
            if (atype == GRETL_TYPE_BUNDLE) {
                nb = gretl_array_get_length(A);
                for (j=0; j<nb; j++) {
                    bj = gretl_array_get_bundle(A, j);
                    bufspace(2*(level+2), prn);
                    pprintf(prn, "%s[%d]:\n", key, j + 1);
                    print_bundle_tree(bj, level + 2, prn);
                }
                continue;
            }
        }
    }

    gretl_array_destroy(K);

    return 0;
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
    int err = real_bundle_print(bundle, 0, 0, prn);

    if (!err) {
        pputc(prn, '\n');
    }

    return err;
}

/* helper function for gretl_bundle_write_constructor() */

static void write_item_constructor (gpointer value, gpointer p)
{
    bundled_item *item = value;
    const gchar *kstr = item->key;
    GretlType t = item->type;
    GString *gs = p;

    if (t == GRETL_TYPE_DOUBLE) {
        double x = *(double *) item->data;

        if (na(x)) {
            g_string_append_printf(gs, "%s=NA,", kstr);
        } else {
            g_string_append_printf(gs, "%s=%g,", kstr, x);
        }
    } else if (t == GRETL_TYPE_INT) {
        int i = *(int *) item->data;

        g_string_append_printf(gs, "%s=%d,", kstr, i);
    } else if (t == GRETL_TYPE_UNSIGNED) {
        guint32 u = *(guint32 *) item->data;

        g_string_append_printf(gs, "%s=%u,", kstr, u);
    } else if (t == GRETL_TYPE_STRING) {
        char *s = (char *) item->data;

	g_string_append_printf(gs, "%s=\"%s\",", kstr, s);
    } else if (t == GRETL_TYPE_MATRIX) {
	gretl_matrix *m = (gretl_matrix *) item->data;
	gchar *ms = gretl_matrix_write_constructor(m);

	g_string_append_printf(gs, "%s=%s,", kstr, ms);
	g_free(ms);
    } else {
	fprintf(stderr, "write bundle constructor: skipping "
		"unhandled member '%s'\n", kstr);
    }
}

/**
 * gretl_bundle_write_constructor:
 * @bundle: gretl bundle.
 *
 * For internal use. This function is designed to create a
 * string that can be used in a hansl command or function
 * call to create a copy of the @bundle argument. In that
 * context, only rather simple bundle elements are supported:
 * scalars, strings and small matrices. Large matrices, as well
 * as members of type bundle, array, series or list, are
 * ignored.
 *
 * Returns: a GLib-allocated string; use g_free() to destroy it.
 */

gchar *gretl_bundle_write_constructor (gretl_bundle *bundle)
{
    int n_items = g_hash_table_size(bundle->ht);
    gchar *ret = NULL;

    if (n_items == 0) {
	return g_strdup("_()");
    } else {
	GList *L = g_hash_table_get_values(bundle->ht);
	gsize sz = 4 + 8 * n_items;
	GString *gs;

	gs = g_string_sized_new(sz);
	g_string_append(gs, "_(");
	gretl_push_c_numeric_locale();
	g_list_foreach(L, write_item_constructor, gs);
	gretl_pop_c_numeric_locale();
	sz = strlen(gs->str);
	if (gs->str[sz-1] == ',') {
	    gs->str[sz-1] = ')';
	} else {
	    g_string_append_c(gs, ')');
	}
	g_list_free(L);
	ret = g_string_free(gs, FALSE);
    }

    return ret;
}

/**
 * gretl_bundle_debug_print:
 * @bundle: gretl bundle.
 * @msg: extra string to print, or NULL.
 *
 * Prints to stderr a list of the keys defined in @bundle, along
 * with descriptive notes, if any. If @msg is non-NULL it is
 * printed first.
 */

void gretl_bundle_debug_print (gretl_bundle *bundle, const char *msg)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

    if (msg != NULL) {
        pputs(prn, msg);
        pputc(prn, '\n');
    }
    real_bundle_print(bundle, 0, 0, prn);
    pputc(prn, '\n');
    gretl_print_destroy(prn);
}

int gretl_bundle_print_tree (gretl_bundle *bundle, PRN *prn)
{
    int err;

    if (bundle == NULL) {
        err = E_DATA;
    } else {
        err = print_bundle_tree(bundle, 0, prn);
        if (!err) {
            pputc(prn, '\n');
        }
    }

    return err;
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
        if (item->data == NULL) {
            fprintf(stderr, "bundle -> XML: skipping NULL string %s\n", key);
            return;
        }
    }

    pprintf(prn, "<bundled-item key=\"%s\" type=\"%s\"", key,
            gretl_type_get_name(item->type));

    if (item->note != NULL) {
        pprintf(prn, " note=\"%s\"", item->note);
    }

    if (item->size > 0) {
        pprintf(prn, " size=\"%d\">\n", item->size);
    } else if (item->type == GRETL_TYPE_STRING ||
               gretl_is_scalar_type(item->type)) {
        pputc(prn, '>');
    } else {
        pputs(prn, ">\n");
    }

    if (item->type == GRETL_TYPE_DOUBLE) {
        double x = *(double *) item->data;

        if (na(x)) {
            pputs(prn, "NA");
        } else {
            pprintf(prn, "%.16g", x);
        }
    } else if (item->type == GRETL_TYPE_INT) {
        int i = *(int *) item->data;

        pprintf(prn, "%d", i);
    } else if (item->type == GRETL_TYPE_UNSIGNED) {
        guint32 u = *(guint32 *) item->data;

        pprintf(prn, "%u", u);
    } else if (item->type == GRETL_TYPE_SERIES) {
        double *vals = (double *) item->data;

        for (j=0; j<item->size; j++) {
            if (na(vals[j])) {
                pputs(prn, "NA ");
            } else {
                pprintf(prn, "%.16g ", vals[j]);
            }
        }
    } else if (item->type == GRETL_TYPE_STRING) {
        gretl_xml_put_string((char *) item->data, prn);
    } else if (item->type == GRETL_TYPE_MATRIX) {
        gretl_matrix *m = (gretl_matrix *) item->data;

        gretl_matrix_serialize(m, NULL, prn);
    } else if (item->type == GRETL_TYPE_BUNDLE) {
        gretl_bundle *b = (gretl_bundle *) item->data;

        gretl_bundle_serialize(b, NULL, prn);
    } else if (item->type == GRETL_TYPE_ARRAY) {
        gretl_array *a = (gretl_array *) item->data;

        gretl_array_serialize(a, prn);
    } else if (item->type == GRETL_TYPE_LIST) {
        int *list = (int *) item->data;

        gretl_list_serialize(list, NULL, prn);
    } else {
        fprintf(stderr, "bundle -> XML: skipping %s\n", key);
    }

    pputs(prn, "</bundled-item>\n");
}

void gretl_bundle_serialize (gretl_bundle *b, const char *name,
                             PRN *prn)
{
    pputs(prn, "<gretl-bundle");
    if (name != NULL) {
        pprintf(prn, " name=\"%s\"", name);
    }
    if (b->creator != NULL && *b->creator != '\0') {
        pprintf(prn, " creator=\"%s\"", b->creator);
    }
    if (b->type == BUNDLE_KALMAN) {
        pputs(prn, " type=\"kalman\"");
    }
    pputs(prn, ">\n");

    if (b->type == BUNDLE_KALMAN) {
        kalman_serialize(b->data, prn);
    }

    if (b->ht != NULL) {
        g_hash_table_foreach(b->ht, xml_put_bundled_item, prn);
    }

    pputs(prn, "</gretl-bundle>\n");
}

static int load_bundled_items (gretl_bundle *b, xmlNodePtr cur, xmlDocPtr doc)
{
    GretlType type;
    char *key;
    int err = 0;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "bundled-item")) {
            key = (char *) xmlGetProp(cur, (XUC) "key");
            type = gretl_xml_get_type_property(cur);
            if (key == NULL || type == 0) {
                err = E_DATA;
            } else {
                int size = 0;

                if (type == GRETL_TYPE_DOUBLE) {
                    double x;

                    if (!gretl_xml_node_get_double(cur, doc, &x)) {
                        err = E_DATA;
                    } else {
                        err = gretl_bundle_set_data(b, key, &x, type, size);
                    }
                } else if (type == GRETL_TYPE_INT) {
                    int i;

                    if (!gretl_xml_node_get_int(cur, doc, &i)) {
                        err = E_DATA;
                    } else {
                        err = gretl_bundle_set_data(b, key, &i, type, size);
                    }
                } else if (type == GRETL_TYPE_UNSIGNED) {
                    guint32 u;

                    if (!gretl_xml_node_get_unsigned(cur, doc, &u)) {
                        err = E_DATA;
                    } else {
                        err = gretl_bundle_set_data(b, key, &u, type, size);
                    }
                } else if (type == GRETL_TYPE_STRING) {
                    char *s;

                    if (!gretl_xml_node_get_trimmed_string(cur, doc, &s)) {
                        err = E_DATA;
                    } else {
                        err = gretl_bundle_donate_data(b, key, s, type, size);
                    }
                } else if (type == GRETL_TYPE_SERIES) {
                    double *xvec = gretl_xml_get_double_array(cur, doc, &size, &err);

                    if (!err) {
                        err = gretl_bundle_donate_data(b, key, xvec, type, size);
                    }
                } else if (type == GRETL_TYPE_MATRIX) {
                    xmlNodePtr child = cur->xmlChildrenNode;
                    gretl_matrix *m;

                    if (child == NULL) {
                        err = E_DATA;
                    } else {
                        m = gretl_xml_get_matrix(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, m, type, size);
                        }
                    }
                } else if (type == GRETL_TYPE_BUNDLE) {
                    xmlNodePtr child = cur->xmlChildrenNode;
                    gretl_bundle *baby;

                    if (child == NULL) {
                        err = E_DATA;
                    } else {
                        baby = gretl_bundle_deserialize(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, baby, type, size);
                        }
                    }
                } else if (type == GRETL_TYPE_ARRAY) {
                    xmlNodePtr child = cur->xmlChildrenNode;
                    gretl_array *a;

                    if (child == NULL) {
                        err = E_DATA;
                    } else {
                        a = gretl_array_deserialize(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, a, type, size);
                        }
                    }
                } else if (type == GRETL_TYPE_LIST) {
                    xmlNodePtr child = cur->xmlChildrenNode;
                    int *list;

                    if (child == NULL) {
                        err = E_DATA;
                    } else {
                        list = gretl_xml_get_list(child, doc, &err);
                        if (!err) {
                            err = gretl_bundle_donate_data(b, key, list, type, size);
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
            }
        }
        cur = cur->next;
    }

    return err;
}

/* For internal use only: @p1 should be of type xmlNodePtr and @p2
   should be an xmlDocPtr. We suppress the actual pointer types in the
   signature so that it's possible for a module to include
   gretl_bundle.h without including the full libxml headers.
*/

gretl_bundle *gretl_bundle_deserialize (void *p1, void *p2,
                                        int *err)
{
    xmlNodePtr cur, node = p1;
    xmlDocPtr doc = p2;
    gretl_bundle *b = NULL;
    char *btype = NULL;
    char *creator = NULL;

    btype = (char *) xmlGetProp(node, (XUC) "type");
    creator = (char *) xmlGetProp(node, (XUC) "creator");

    cur = node->xmlChildrenNode;

    if (btype != NULL && !strcmp(btype, "kalman")) {
        b = kalman_deserialize(cur, doc, err);
    } else {
        b = gretl_bundle_new();
        if (b == NULL) {
            *err = E_ALLOC;
        } else if (creator != NULL) {
            b->creator = creator;
            creator = NULL;
        }
    }

    free(btype);
    free(creator);

    if (b != NULL) {
        *err = load_bundled_items(b, cur, doc);
        if (*err) {
            fprintf(stderr, "deserialize bundle: "
                    "bundle is broken (err = %d)\n", *err);
            gretl_bundle_destroy(b);
            b = NULL;
        }
    }

    return b;
}

static int call_bundle_to_json (gretl_bundle *b,
                                const char *fname,
                                int control)
{
    int (*jfunc) (gretl_bundle *, const char *, gretlopt);

    jfunc = get_plugin_function("bundle_to_json");
    if (jfunc == NULL) {
        return E_FOPEN;
    } else {
        gretlopt opt = OPT_NONE;

        if (control & 2) {
            opt |= OPT_A;
        }
        if (control & 4) {
            opt |= OPT_P;
        }
        return jfunc(b, fname, opt);
    }
}

int gretl_bundle_write_to_file (gretl_bundle *b,
                                const char *fname,
                                int control)
{
    char fullname[FILENAME_MAX];
    PRN *prn = NULL;
    int err = 0;

    if (control & 1) {
        gretl_build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
        strcpy(fullname, fname);
        gretl_maybe_switch_dir(fullname);
    }

    if (has_suffix(fname, ".json") || has_suffix(fname, ".geojson")) {
        return call_bundle_to_json(b, fullname, control);
    }

    if (has_suffix(fname, ".gz")) {
        prn = gretl_gzip_print_new(fullname, -1, &err);
    } else {
        prn = gretl_print_new_with_filename(fullname, &err);
    }

    if (prn != NULL) {
        gretl_push_c_numeric_locale();
        gretl_xml_header(prn);
        gretl_bundle_serialize(b, NULL, prn);
        gretl_print_destroy(prn);
        gretl_pop_c_numeric_locale();
    }

    return err;
}

char *gretl_bundle_write_to_buffer (gretl_bundle *b,
                                    int rank,
                                    int *bytes,
                                    int *err)
{
    char *buf = NULL;
    PRN *prn;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, err);

    if (!*err) {
        gretl_push_c_numeric_locale();
        gretl_xml_header(prn);
        gretl_bundle_serialize(b, NULL, prn);
        buf = gretl_print_steal_buffer(prn);
        *bytes = strlen(buf) + 1;
        gretl_pop_c_numeric_locale();
        gretl_print_destroy(prn);
    }

    return buf;
}

static gretl_bundle *read_json_bundle (const char *fname,
                                       int *err)
{
    gretl_bundle *(*jfunc) (const char *, const char *, int *);
    gretl_bundle *b = NULL;
    GError *gerr = NULL;
    gchar *JSON = NULL;
    gsize len = 0;
    gboolean ok;

    ok = g_file_get_contents(fname, &JSON, &len, &gerr);

    if (ok) {
        jfunc = get_plugin_function("json_get_bundle");
        if (jfunc == NULL) {
            *err = E_FOPEN;
        } else {
            b = jfunc(JSON, NULL, err);
        }
        g_free(JSON);
    } else if (gerr != NULL) {
        *err = E_FOPEN;
        gretl_errmsg_set(gerr->message);
        g_error_free(gerr);
    }

    return b;
}

static gretl_bundle *read_shapefile_bundle (const char *fname,
                                            int *err)
{
    gretl_bundle *(*bfunc) (const char *, int *);
    gretl_bundle *b = NULL;

    bfunc = get_plugin_function("shp_get_bundle");
    if (bfunc == NULL) {
        *err = E_FOPEN;
    } else {
        b = bfunc(fname, err);
    }

    return b;
}

gretl_bundle *gretl_bundle_read_from_file (const char *fname,
                                           int from_dotdir,
                                           int *err)
{
    char fullname[FILENAME_MAX];
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    gretl_bundle *b = NULL;

    if (2 * strlen(fname) >= sizeof fullname) {
	gretl_errmsg_sprintf("Filename too long, %d bytes",
			     (int) strlen(fname));
	*err = E_DATA;
	return NULL;
    }

    if (from_dotdir) {
        gretl_build_path(fullname, gretl_dotdir(), fname, NULL);
    } else {
        strcpy(fullname, fname);
        gretl_maybe_prepend_dir(fullname);
    }

    if (has_suffix(fname, ".json") || has_suffix(fname, ".geojson")) {
        b = read_json_bundle(fullname, err);
    } else if (has_suffix(fname, ".shp")) {
        b = read_shapefile_bundle(fullname, err);
    } else {
        *err = gretl_xml_open_doc_root(fullname, "gretl-bundle", &doc, &cur);
        if (!*err) {
            gretl_push_c_numeric_locale();
            b = gretl_bundle_deserialize(cur, doc, err);
            gretl_pop_c_numeric_locale();
            xmlFreeDoc(doc);
        }
    }

    if (*err && b != NULL) {
        gretl_bundle_destroy(b);
        b = NULL;
    }

    return b;
}

gretl_bundle *gretl_bundle_read_from_buffer (const char *buf,
                                             int len,
                                             int *err)
{
    xmlDocPtr doc = NULL;
    gretl_bundle *b = NULL;

    doc = xmlReadMemory(buf, len, NULL, NULL,
                        XML_PARSE_HUGE | XML_PARSE_NOBLANKS);

    if (doc == NULL) {
        gretl_errmsg_set(_("xmlReadMemory failed"));
        *err = 1;
    } else {
        xmlNodePtr cur = xmlDocGetRootElement(doc);

        if (cur == NULL) {
            gretl_errmsg_set(_("xmlDocGetRootElement failed"));
            *err = 1;
        } else {
            gretl_push_c_numeric_locale();
            b = gretl_bundle_deserialize(cur, doc, err);
            gretl_pop_c_numeric_locale();
        }
        xmlFreeDoc(doc);
    }

    if (*err && b != NULL) {
        gretl_bundle_destroy(b);
        b = NULL;
    }

    return b;
}

/* get the key strings from @b in the form of a gretl_array */

gretl_array *gretl_bundle_get_keys (gretl_bundle *b, int *err)
{
    gretl_array *A = NULL;
    int myerr = 0;

    if (b == NULL || b->ht == NULL) {
        myerr = E_DATA;
    } else {
        GList *keys = g_hash_table_get_keys(b->ht);
        guint n;

        if (keys != NULL && (n = g_list_length(keys)) > 0) {
            A = gretl_array_new(GRETL_TYPE_STRINGS, n, &myerr);
            if (!myerr) {
                GList *L = g_list_first(keys);
                int i = 0;

                while (L != NULL) {
                    gretl_array_set_string(A, i++, L->data, 1);
                    L = L->next;
                }
            }
        } else {
            A = gretl_array_new(GRETL_TYPE_STRINGS, 0, &myerr);
        }
        if (keys != NULL) {
            g_list_free(keys);
        }
    }

    if (err != NULL) {
        *err = myerr;
    }

    return A;
}

/* Get the key strings from @b in the form of a "raw" array
   of C type char *, terminated by a NULL pointer.
*/

char **gretl_bundle_get_keys_raw (gretl_bundle *b, int *ns)
{
    char **S = NULL;

    *ns = 0;

    if (b != NULL && b->ht != NULL) {
        GList *keys = g_hash_table_get_keys(b->ht);
        guint n;

        if (keys != NULL && (n = g_list_length(keys)) > 0) {
            S = strings_array_new(n + 1);
            if (S != NULL) {
                GList *L = g_list_first(keys);
                int i = 0;

                while (L != NULL) {
                    S[i++] = gretl_strdup(L->data);
                    L = g_list_next(L);
                }
                *ns = n;
                S[i] = NULL; /* terminator */
            }
        }
        if (keys != NULL) {
            g_list_free(keys);
        }
    }

    return S;
}

gretl_bundle *get_sysinfo_bundle (int *err)
{
    gretl_matrix *memvals = NULL;

    if (sysinfo_bundle == NULL) {
        gretl_bundle *b = gretl_bundle_new();

        if (b == NULL) {
            *err = E_ALLOC;
        } else {
            gretl_bundle *fb;
            char *s1, *s2, *s3 = NULL;
            int ival = 0;

#if HAVE_MPI
            ival = check_for_mpiexec();
#endif
            gretl_bundle_set_scalar(b, "mpi", (double) ival);
            ival = gretl_max_mpi_processes();
            gretl_bundle_set_scalar(b, "mpimax", (double) ival);
            ival = gretl_n_processors();
            gretl_bundle_set_scalar(b, "nproc", (double) ival);
            ival = gretl_n_physical_cores();
            gretl_bundle_set_scalar(b, "ncores", (double) ival);
            ival = gretl_in_gui_mode();
            gretl_bundle_set_scalar(b, "gui_mode", (double) ival);
            ival = 0;
#ifdef _OPENMP
            ival = 1;
#endif
            gretl_bundle_set_scalar(b, "omp", (double) ival);
            ival = sizeof(void*) == 8 ? 64 : 32;
            gretl_bundle_set_scalar(b, "wordlen", (double) ival);
#if defined(G_OS_WIN32)
            gretl_bundle_set_string(b, "os", "windows");
#elif defined(OS_OSX)
            gretl_bundle_set_string(b, "os", "macos");
#elif defined(linux)
            gretl_bundle_set_string(b, "os", "linux");
#else
            gretl_bundle_set_string(b, "os", "other");
#endif
            gretl_bundle_set_string(b, "hostname", g_get_host_name());
            gretl_bundle_set_string(b, "blas", blas_variant_string());
	    gretl_bundle_set_string(b, "gnuplot", gnuplot_version_string());
            if (get_blas_details(&s1, &s2, &s3)) {
                gretl_bundle_set_string(b, "blascore", s1);
                gretl_bundle_set_string(b, "blas_parallel", s2);
                if (s3 != NULL) {
                    gretl_bundle_set_string(b, "blas_version", s3);
                }
            }
#if defined(COMPILER_IDENT)
            gretl_bundle_set_string(b, "compiler", COMPILER_IDENT);
#endif
	    /* allocated string */
	    gretl_bundle_donate_data(b, "cpuid", get_cpu_details(),
				     GRETL_TYPE_STRING, 0);
	    /* information pertaining to 'foreign' programs */
            fb = foreign_info();
            if (fb != NULL) {
                gretl_bundle_donate_data(b, "foreign", fb,
                                         GRETL_TYPE_BUNDLE, 0);
            }
        }
        sysinfo_bundle = b;
    }

    memvals = gretl_matrix_alloc(1, 2);
    if (memvals != NULL) {
        char **S = malloc(2 * sizeof *S);

        memory_stats(memvals->val);
        S[0] = gretl_strdup("MBtotal");
        S[1] = gretl_strdup("MBfree");
        gretl_matrix_set_colnames(memvals, S);
        gretl_bundle_donate_data(sysinfo_bundle, "mem", memvals,
                                 GRETL_TYPE_MATRIX, 0);
    }

    return sysinfo_bundle;
}

void *sysinfo_bundle_get_data (const char *key, GretlType *type,
                               int *err)
{
    gretl_bundle *b = get_sysinfo_bundle(err);
    void *ret = NULL;

    if (b != NULL) {
        ret = gretl_bundle_get_data(b, key, type, NULL, err);
    }

    return ret;
}

/* For a single-equation model, create a bundle containing
   all the data available via $-accessors.
*/

gretl_bundle *bundle_from_model (MODEL *pmod,
                                 DATASET *dset,
                                 int *err)
{
    gretl_bundle *b = NULL;
    gretl_matrix *m;
    gretl_array *a;
    double *x;
    double val;
    int *list;
    const char *s;
    const char *key;
    int i, t, berr;

    if (pmod == NULL) {
        GretlObjType type = 0;
        void *p = get_genr_model(&type);

        if (p == NULL || type != GRETL_OBJ_EQN) {
            p = get_last_model(&type);
            if (p == NULL || type != GRETL_OBJ_EQN) {
                gretl_errmsg_sprintf(_("%s: no data available"), "$model");
                *err = E_DATA;
                return NULL;
            }
        }
        pmod = p;
    }

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    b = gretl_bundle_new();
    if (b == NULL) {
        free(x);
        *err = E_ALLOC;
        return NULL;
    }

    for (i=M_ESS; i<M_SCALAR_MAX && !*err; i++) {
        /* model scalars */
        berr = 0;
        if (i == M_DWPVAL) {
            /* Durbin-Watson p-value: don't include this unless
               it has already been computed and attached to the
               model, since it may involve heavy computation.
            */
            val = gretl_model_get_double(pmod, "dwpval");
            if (!na(val)) {
                *err = gretl_bundle_set_scalar(b, "dwpval", val);
            }
        } else {
            val = gretl_model_get_scalar(pmod, i, dset, &berr);
            if (!berr) {
                key = mvarname(i) + 1;
                *err = gretl_bundle_set_scalar(b, key, val);
            }
        }
    }

    for (i=M_SCALAR_MAX+1; i<M_SERIES_MAX && !*err; i++) {
        /* model series */
        for (t=0; t<dset->n; t++) {
            x[t] = NADBL;
        }
        berr = gretl_model_get_series(x, pmod, dset, i);
        if (!berr) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_set_series(b, key, x, dset->n);
        }
    }

    for (i=M_SERIES_MAX+1; i<M_MATRIX_MAX && !*err; i++) {
        /* model matrices */
        berr = 0;
        m = gretl_model_get_matrix(pmod, i, &berr);
        if (!berr) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_donate_data(b, key, m,
                                            GRETL_TYPE_MATRIX,
                                            0);
        }
    }

    for (i=M_MBUILD_MAX+1; i<M_LIST_MAX && !*err; i++) {
        /* model lists */
        list = NULL;
        if (i == M_XLIST) {
            list = gretl_model_get_x_list(pmod);
        } else if (i == M_YLIST) {
            list = gretl_model_get_y_list(pmod);
        }
        if (list != NULL) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_donate_data(b, key, list,
                                            GRETL_TYPE_LIST,
                                            0);
        }
    }

    for (i=M_LIST_MAX+1; i<M_PARNAMES && !*err; i++) {
        /* model strings */
        s = NULL;
        if (i == M_DEPVAR) {
            s = gretl_model_get_depvar_name(pmod, dset);
        } else if (i == M_COMMAND) {
            s = gretl_command_word(pmod->ci);
        }
        if (s != NULL && *s != '\0') {
            key = mvarname(i) + 1;
            *err = gretl_bundle_set_string(b, key, s);
        }
    }

    for (i=M_PARNAMES; i<M_MAX && !*err; i++) {
        /* model arrays */
        a = NULL;
        if (i == M_PARNAMES) {
            a = gretl_model_get_param_names(pmod, dset, err);
        }
        if (a != NULL) {
            key = mvarname(i) + 1;
            *err = gretl_bundle_donate_data(b, key, a,
                                            GRETL_TYPE_ARRAY,
                                            0);
        }
    }

    if (!*err) {
        /* miscellaneous data items */
        *err = bundlize_model_data_items(pmod, b);
    }

    free(x);

    /* don't return a broken bundle */
    if (*err && b != NULL) {
        gretl_bundle_destroy(b);
        b = NULL;
    }

    return b;
}

/* For an estimated system of some sort, create a bundle containing
   relevant data.
*/

gretl_bundle *bundle_from_system (void *ptr,
                                  int type,
                                  DATASET *dset,
                                  int *err)
{
    GRETL_VAR *var = NULL;
    GretlObjType otype = type;
    equation_system *sys = NULL;
    gretl_bundle *b = NULL;

    if (ptr == NULL) {
        ptr = get_genr_model(&otype);
        if (ptr == NULL || otype == GRETL_OBJ_EQN) {
            ptr = get_last_model(&otype);
        }
    }

    if (ptr == NULL) {
        gretl_errmsg_sprintf(_("%s: no data available"), "$system");
        *err = E_DATA;
    } else if (otype == GRETL_OBJ_VAR) {
        var = (GRETL_VAR *) ptr;
    } else if (otype == GRETL_OBJ_SYS) {
        sys = (equation_system *) ptr;
    } else {
        gretl_errmsg_sprintf(_("%s: no data available"), "$system");
        *err = E_DATA;
    }

    if (!*err) {
        b = gretl_bundle_new();
        if (b == NULL) {
            *err = E_ALLOC;
            return NULL;
        }
    }

    if (var != NULL) {
        *err = gretl_VAR_bundlize(var, dset, b);
    } else if (sys != NULL) {
        *err = equation_system_bundlize(sys, b);
    }

    /* don't return a broken bundle */
    if (*err && b != NULL) {
        gretl_bundle_destroy(b);
        b = NULL;
    }

    return b;
}

gretl_bundle *kalman_bundle_new (gretl_matrix *M[],
                                 int copy[], int nmat,
                                 int dkvar, int *err)
{
    gretl_bundle *b = gretl_bundle_new();

    if (b == NULL) {
        *err = E_ALLOC;
    } else {
        b->type = BUNDLE_KALMAN;
        b->data = kalman_new_minimal(M, copy, nmat, dkvar, err);
    }

    /* don't return a broken bundle */
    if (*err && b != NULL) {
        gretl_bundle_destroy(b);
        b = NULL;
    }

    return b;
}

static gint sort_bundled_items (const void *a, const void *b)
{
    const bundled_item *ia = a;
    const bundled_item *ib = b;
    int ta = gretl_type_get_order(ia->type);
    int tb = gretl_type_get_order(ib->type);
    int ret = ta - tb;

    if (ret == 0) {
        ret = g_ascii_strcasecmp(ia->key, ib->key);
    }

    return ret;
}

GList *gretl_bundle_get_sorted_items (gretl_bundle *b)
{
    GList *blist;

    blist = g_hash_table_get_values(b->ht);
    blist = g_list_sort(blist, sort_bundled_items);

    return blist;
}

int gretl_bundles_are_equal (gretl_bundle *b1,
			     gretl_bundle *b2)
{
    GList *k1 = g_hash_table_get_keys(b1->ht);
    void *v1, *v2;
    const char *s1;
    GretlType t1, t2;
    int sz1, sz2;
    int err;
    int ret = 1;

    while (k1 != NULL && ret) {
	s1 = k1->data;
	v1 = gretl_bundle_get_data(b1, s1, &t1, &sz1, &err);
	if (!gretl_bundle_has_key(b2, s1)) {
	    ret = 0;
	} else {
	    v2 = gretl_bundle_get_data(b2, s1, &t2, &sz2, &err);
	    if (t2 != t1 || sz1 != sz2) {
		ret = 0;
	    } else if (gretl_is_scalar_type(t1)) {
		double d1 = gretl_bundle_get_scalar(b1, s1, NULL);
		double d2 = gretl_bundle_get_scalar(b2, s1, NULL);
		int nas = na(d1) + na(d2);

		if (d1 != d2 && nas != 2) {
		    ret = 0;
		}
	    } else if (t1 == GRETL_TYPE_SERIES) {
		double *x1 = v1;
		double *x2 = v2;
		int i, nas;

		for (i=0; i<sz1; i++) {
		    nas = na(x1[i]) + na(x2[i]);
		    if (x1[i] != x2[i] && nas != 2) {
			ret = 0;
			break;
		    }
		}
	    } else if (t1 == GRETL_TYPE_MATRIX) {
		int eq = gretl_matrices_are_equal(v1, v2, 0, &err);

		if (eq != 1) {
		    ret = 0;
		}
	    } else if (t1 == GRETL_TYPE_STRING) {
		if (strcmp(v1, v2)) {
		    ret = 0;
		}
	    } else if (t1 == GRETL_TYPE_BUNDLE) {
		if (!gretl_bundles_are_equal(v1, v2)) {
		    ret = 0;
		}
	    } else if (t1 == GRETL_TYPE_ARRAY) {
		if (!gretl_arrays_are_equal(v1, v2)) {
		    ret = 0;
		}
	    } else {
		/* must be GRETL_TYPE_LIST */
		if (gretl_list_cmp(v1, v2)) {
		    ret = 0;
		}
	    }
	}
	k1 = k1->next;
    }

    if (ret) {
	GList *k2 = g_hash_table_get_keys(b2->ht);

	while (k2 != NULL && ret) {
	    if (!gretl_bundle_has_key(b1, k2->data)) {
		ret = 0;
	    }
	    k2 = k2->next;
	}
	g_list_free(k2);
    }

    g_list_free(k1);

    return ret;
}

void gretl_bundle_cleanup (void)
{
    if (sysinfo_bundle != NULL) {
        gretl_bundle_destroy(sysinfo_bundle);
        sysinfo_bundle = NULL;
    }
}
