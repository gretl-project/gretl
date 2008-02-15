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

int default_lag_order (const DATAINFO *pdinfo);

int is_standard_lag (int v, const DATAINFO *pdinfo, int *parent);

int is_dummy_child (int v, const DATAINFO *pdinfo, int *parent);

int diffgenr (int v, int ci, double ***pZ, DATAINFO *pdinfo);

int laggenr (int v, int lag, double ***pZ, DATAINFO *pdinfo);

int loggenr (int v, double ***pZ, DATAINFO *pdinfo);

int invgenr (int v, double ***pZ, DATAINFO *pdinfo);

int xpxgenr (int vi, int vj, double ***pZ, DATAINFO *pdinfo);

int list_diffgenr (int *list, int ci, double ***pZ, DATAINFO *pdinfo);

int list_laggenr (int **plist, int order, double ***pZ, DATAINFO *pdinfo);

int *laggenr_from_to (int v, int minlag, int maxlag, double ***pZ, 
		      DATAINFO *pdinfo, int *err);

int list_loggenr (int *list, double ***pZ, DATAINFO *pdinfo);

int list_xpxgenr (int **plist, double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt);

int list_dumgenr (int **plist, double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt);

int list_makediscrete (const int *list, DATAINFO *pdinfo, gretlopt opt);

int gettrend (double ***pZ, DATAINFO *pdinfo, int square);

void gretl_transforms_cleanup (void);

#endif /* TRANSFORMS_H */
