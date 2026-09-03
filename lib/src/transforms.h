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

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

int default_lag_order (const DATASET *dset);

int standard_lag_of (int v, int parent, const DATASET *dset);

int is_standard_diff (int v, const DATASET *dset, int *parent);

int diffgenr (int v, int ci, DATASET *dset);

int laggenr (int v, int lag, DATASET *dset);

int loggenr (int v, DATASET *dset);

int invgenr (int v, DATASET *dset);

int xpxgenr (int vi, int vj, DATASET *dset);

int list_diffgenr (int *list, int ci, DATASET *dset);

int list_orthdev (int *list, DATASET *dset);

int list_resample (int *list, DATASET *dset);

int list_dropcoll (int *list, double eps, DATASET *dset);

int list_laggenr (int **plist, int lmin, int lmax,
		  const gretl_matrix *lvec,
		  DATASET *dset, int compfac,
		  gretlopt opt);

int *laggenr_from_to (int v, int fromlag, int tolag, 
		      DATASET *dset, int *err);

int list_loggenr (int *list, DATASET *dset);

int list_stdgenr (int *list, DATASET *dset, gretlopt opt);

int list_xpxgenr (int **plist, DATASET *dset, gretlopt opt);

int list_dumgenr (int **plist, DATASET *dset, gretlopt opt);

int dumgenr_with_oddval (int **plist, DATASET *dset, double oddval);

int auto_dummify_list (int **plist, DATASET *dset);

int list_makediscrete (const int *list, DATASET *dset, gretlopt opt);

int hf_list_diffgenr (int *list, int ci, double parm, DATASET *dset);

void gretl_transforms_cleanup (void);

#endif /* TRANSFORMS_H */
