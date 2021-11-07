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

#ifndef GRETL_TDISAGG_H
#define GRETL_TDISAGG_H

gretl_matrix *
get_tdisagg_matrix (int ynum, const int *ylist, const double *yval,
		    int xnum, const int *xlist, const double *xval,
		    int xmidas, int fac, GretlType targ,
		    gretl_matrix *Y, gretl_matrix *X,
		    DATASET *dset, gretl_bundle *b,
		    gretl_bundle *r, PRN *prn, int *err);

gretl_matrix *matrix_chowlin (const gretl_matrix *Y,
			      const gretl_matrix *X,
			      int s, int *err);

#endif /* GRETL_TDISAGG_H */ 
