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

#ifndef TSLS_H
#define TSLS_H

void tsls_free_data (const MODEL *pmod);

int ivreg_process_lists (const int *list, int **reglist, int **instlist);

int *ivreg_list_omit (const int *orig, const int *drop, gretlopt opt, int *err);

int *ivreg_list_add (const int *orig, const int *add, gretlopt opt, int *err);

int *tsls_make_endolist (const int *reglist, const int *instlist,
			 int *err);

MODEL tsls (const int *list, DATASET *dset, gretlopt opt);

#endif /* TSLS_H */
