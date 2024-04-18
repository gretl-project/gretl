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

#ifndef GFN_ARGLISTS_H
#define GFN_ARGLISTS_H

/* name for ad hoc list based on main-window selection */
#define AUTOLIST "LTmp___"

typedef struct arglist_ arglist;

arglist *arglist_new (const char *pkgname, const void *func, int argc);

arglist *arglist_lookup (const char *pkgname, const void *func);

int arglist_record_arg (arglist *a, int i, const char *val);

const char *arglist_lookup_val (arglist *a, int i);

void arglist_cleanup (void);

#endif /* GFN_ARGLISTS_H */
