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

typedef struct gretl_bundle_ gretl_bundle;
typedef struct bundle_value_ bundle_value;

struct gretl_bundle_ {
    char name[VNAMELEN];
    GHashTable *ht;
    int level;
};

struct bundle_value_ {
    GretlType type;
    gpointer data;
};

static gretl_bundle **bundles;
static int n_bundles;
static int bundle_imin;

static bundle_value *bundle_value_new (GretlType type, void *ptr, int *err)
{
    bundle_value *val = malloc(sizeof *val);

    if (val == NULL) {
	*err = E_ALLOC;
    } else {
	double *dp;

	val->type = type;
	switch (val->type) {
	case GRETL_TYPE_DOUBLE:
	    val->data = malloc(sizeof(double));
	    if (val->data != NULL) {
		dp = val->data;
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

static void free_hash_key (gpointer data)
{
    return;
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

static gretl_bundle *get_bundle_pointer (const char *s, int level)
{
    int i;

    for (i=bundle_imin; i<n_bundles; i++) {
	if (bundles[i]->level == level &&
	    !strcmp(s, bundles[i]->name)) {
	    return bundles[i];
	}
    }

    return NULL;
}

static int gretl_bundle_push (gretl_bundle *s)
{
    int n = n_bundles + 1;

    if (reallocate_bundles(n)) {
	free(s);
	return E_ALLOC;
    }

    bundles[n-1] = s;
    set_n_bundles(n);

    return 0;
}

static int real_delete_bundle (int i)
{
    int n = n_bundles - 1;
    int err = 0;

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

int n_saved_bundles (void)
{
    return n_bundles;
}

int gretl_is_bundle (const char *name)
{
    int ret = 0;

    if (name == NULL || *name == '\0') {
	return 0;
    }

    ret = (get_bundle_pointer(name, gretl_function_depth()) != NULL);

    if (!ret) {
	ret = const_lookup(name);
    }

    return ret;
}

int gretl_bundle_get_index (const char *name, int *err)
{
    int level = gretl_function_depth();
    int i;

    for (i=0; i<n_bundles; i++) {
	if (level == bundles[i]->level &&
	    !strcmp(name, bundles[i]->name)) {
	    return i;
	}
    }

    *err = E_UNKVAR;

    return -1;
}

const char *gretl_bundle_get_name (int i)
{
    if (i >= 0 && i < n_bundles) {
	return bundles[i]->name;
    } else {
	return NULL;
    }
}

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

int gretl_bundle_set_data (const char *name, const char *key,
			   void *ptr, GretlType type)
{
    gretl_bundle *b;
    int err = 0;

    b = get_bundle_pointer(name, gretl_function_depth());
    
    if (b == NULL) {
	err = E_UNKVAR;
    } else {
	bundle_value *val = bundle_value_new(type, ptr, &err);

	if (!err) {
	    if (g_hash_table_lookup(b->ht, key) != NULL) {
		g_hash_table_replace(b->ht, (gpointer) key, val);
	    } else {
		g_hash_table_insert(b->ht, (gpointer) key, val);
	    }
	}
    }

    return err;
}

/* check that we're not colliding with a pre-existing object of
   different type */

static int gretl_bundle_check_name (const char *s, const DATAINFO *pdinfo)
{
    int v, err = 0;

    if ((v = series_index(pdinfo, s)) < pdinfo->v) {
	err = E_TYPES;
    } else if (get_matrix_by_name(s)) {
	err = E_TYPES;
    } else if (get_list_by_name(s)) {
	err = E_TYPES;
    } else if (get_string_by_name(s)) {
	err = E_TYPES;
    } else if (gretl_is_scalar(s)) {
	err = E_TYPES;
    } else {
	err = check_varname(s);
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
		name, b->level);
	return E_DATA;
    }

    b = malloc(sizeof *b);
    if (b == NULL) {
	return E_ALLOC;
    }

    strcpy(b->name, name);
    b->ht = g_hash_table_new_full(g_str_hash, g_str_equal, 
				  free_hash_key, bundle_value_destroy);
    b->level = gretl_function_depth();

    return gretl_bundle_push(b);
}

int gretl_bundle_add_with_check (const char *name, const DATAINFO *pdinfo)
{
    int err = gretl_bundle_check_name(name, pdinfo);

    if (!err) {
	err = gretl_bundle_add(name);
    }

    return err;
}

int gretl_bundle_set_local_name (int i, const char *name)
{
    if (i >= 0 && i < n_bundles) {
	bundles[i]->name[0] = '\0';
	strncat(bundles[i]->name, name, VNAMELEN - 1);
	bundles[i]->level += 1;
	return 0;
    } else {
	return E_DATA;
    }
}

int gretl_bundle_restore_name (int i, const char *name)
{
    if (i >= 0 && i < n_bundles) {
	bundles[i]->name[0] = '\0';
	strncat(bundles[i]->name, name, VNAMELEN - 1);
	bundles[i]->level -= 1;
	return 0;
    } else {
	return E_DATA;
    }
}

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

