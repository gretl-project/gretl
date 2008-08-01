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

#ifndef GRETL_SCALAR_H
#define GRETL_SCALAR_H

int gretl_is_scalar (const char *name);

double gretl_scalar_get_value (const char *name);

void gretl_scalar_set_value (const char *name, double val);

int gretl_scalar_add (const char *name, double val);

int gretl_scalar_delete (const char *name);

int destroy_user_scalars_at_level (int level);

void destroy_user_scalars (void);

void set_auxiliary_scalars (void);

void unset_auxiliary_scalars (void);

#endif /* GRETL_SCALAR_H */
