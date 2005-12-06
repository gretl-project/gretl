/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef TSLS_H
#define TSLS_H

double *tsls_get_Xi (const MODEL *pmod, double **Z, int i);

void tsls_free_data (const MODEL *pmod);

int *tsls_list_omit (const int *orig, const int *drop, gretlopt opt, int *err);

int *tsls_list_add (const int *orig, const int *add, gretlopt opt, int *err);

MODEL tsls_func (const int *list, int ci, 
		 double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt);

#endif /* TSLS_H */
