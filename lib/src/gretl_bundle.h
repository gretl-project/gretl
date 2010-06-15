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

int gretl_is_bundle (const char *name);

void *gretl_bundle_get_data (const char *name, const char *key,
			     GretlType *type);

int gretl_bundle_set_data (const char *name, const char *key,
			   void *ptr, GretlType type);

int gretl_bundle_add (const char *name);

int gretl_bundle_delete (const char *name, PRN *prn);

int gretl_bundle_mark_as_return (const char *name);

int gretl_bundle_name_return (const char *name);

int destroy_user_bundles_at_level (int level);

void destroy_user_bundles (void);



