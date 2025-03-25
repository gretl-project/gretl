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
#include "gretl_data_store.h"

static GHashTable *data_store;

static void destroy_stored_dset (gpointer data)
{
    destroy_dataset(data);
}

void gretl_data_store_add (DATASET *dset, const char *key)
{
    if (data_store == NULL) {
        data_store = g_hash_table_new_full(g_str_hash, g_str_equal,
                                           g_free, destroy_stored_dset);
    }
    if (data_store != NULL) {
        gchar *id = g_strdup(key);

        g_hash_table_insert(data_store, id, dset);
    }
}

DATASET *gretl_data_store_get (const char *key, int *err)
{
    DATASET *dset = NULL;

    if (data_store == NULL) {
        *err = E_DATA;
    } else {
        dset = g_hash_table_lookup(data_store, key);
        if (dset == NULL) {
            *err = E_DATA;
        }
    }

    return dset;
}

int gretl_data_store_contains (const char *key)
{
    if (data_store != NULL) {
        return g_hash_table_lookup(data_store, key) != NULL;
    } else {
        return 0;
    }
}

int gretl_data_store_get_size (void)
{
    if (data_store != NULL) {
        return g_hash_table_size(data_store);
    } else {
        return 0;
    }
}

void gretl_data_store_remove (const char *key)
{
    if (data_store != NULL) {
        g_hash_table_remove(data_store, key);
    }
}

void gretl_data_store_destroy (void)
{
    if (data_store != NULL) {
        g_hash_table_destroy(data_store);
        data_store = NULL;
    }
}

gchar *gretl_data_store_new_id (void)
{
    gint64 t = g_get_monotonic_time();

    return g_strdup_printf("ds_%" G_GINT64_FORMAT, t);
}
