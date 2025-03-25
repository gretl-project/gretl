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

#ifndef GRETL_DATA_STORE_H
#define GRETL_DATA_STORE_H

#ifdef  __cplusplus
extern "C" {
#endif

void gretl_data_store_add (DATASET *dset, const char *key);

DATASET *gretl_data_store_get (const char *key, int *err);

int gretl_data_store_contains (const char *key);

int gretl_data_store_get_size (void);

void gretl_data_store_remove (const char *key);

void gretl_data_store_destroy (void);

gchar *gretl_data_store_new_id (void);

#ifdef  __cplusplus
}
#endif

#endif /* GRETL_DATA_STORE_H */
