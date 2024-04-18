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

#ifndef GRETL_TYPEMAP_H_
#define GRETL_TYPEMAP_H_

#define NUMERIC_TYPE(t) (t == GRETL_TYPE_INT || \
			 t == GRETL_TYPE_DOUBLE || \
			 t == GRETL_TYPE_SERIES || \
			 t == GRETL_TYPE_USERIES || \
			 t == GRETL_TYPE_LIST || \
			 t == GRETL_TYPE_MATRIX)

GretlType gretl_type_get_plural (GretlType type);

GretlType gretl_type_get_singular (GretlType type);

GretlType gretl_type_get_ref_type (GretlType type);

GretlType gretl_type_get_plain_type (GretlType type);

const char *gretl_type_get_name (GretlType type);

GretlType gretl_type_from_string (const char *s);

GretlType gretl_get_gen_type (const char *s);

int gretl_type_get_order (GretlType type);

int gretl_is_array_type (GretlType type);

int gretl_is_array_ref_type (GretlType type);

int gretl_is_arrayable_type (GretlType type);

int gretl_is_scalar_type (GretlType type);

int gretl_is_series_type (GretlType type);

int gretl_type_mismatch (GretlType t1, GretlType t2);

void gretl_typemap_cleanup (void);

gchar *name_conflict_message (const char *name, GretlType type);

#endif /* GRETL_TYPEMAP_H_ */
