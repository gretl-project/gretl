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

#ifndef FNSAVE_H
#define FNSAVE_H

void edit_function_package (const char *fname);

void edit_new_function_package (void);

int save_function_package (const char *fname, gpointer p);

int save_function_package_as_script (const char *fname, gpointer p);

int save_function_package_spec (const char *fname, gpointer p);

int no_user_functions_check (void);

void get_default_package_name (char *fname, gpointer p, int mode);

gchar *package_sample_get_script (windata_t *vwin);

void update_sample_script (windata_t *vwin);

int update_func_code (windata_t *vwin);

void edit_package_at_startup (const char *fname);

void revise_function_package (void *p);

#endif /* FNSAVE_H */
