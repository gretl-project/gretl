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

/* Information from 'genr' which is needed for documentation
   and certain other purposes
*/

#ifndef GEN_PUBLIC_H
#define GEN_PUBLIC_H

/* name lookup functions */
const char *constname (int c);
const char *dvarname (int t);
const char *mvarname (int t);
const char *bvarname (int t);
const char *dumname (int t);
int is_gretl_accessor (const char *s);
int mvar_lookup (const char *s);

/* helper functions for manual, gretl.lang file */
int gen_func_count (void);
const char *gen_func_name (int i);
int model_var_count (void);
const char *model_var_name (int i);
int data_var_count (void);
const char *data_var_name (int i);
int bundle_var_count (void);
const char *bundle_var_name (int i);
int gretl_const_count (void);
const char *gretl_const_name (int i);

/* helpers for gretl_func.c */
int install_function_override (const char *funname,
			       const char *pkgname,
			       gpointer data);
int delete_function_override (const char *funname,
			      const char *pkgname);

#endif /* GEN_PUBLIC_H */

