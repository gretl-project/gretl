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

#ifndef ADF_KPSS_H
#define ADF_KPSS_H

int adf_test (int order, const int *list, DATASET *dset, 
	      gretlopt opt, PRN *prn);

int kpss_test (int order, const int *list, DATASET *dset, 
	       gretlopt opt, PRN *prn);

int levin_lin_test (int vnum, const int *plist, DATASET *dset, 
		    gretlopt opt, PRN *prn);

int engle_granger_test (int order, const int *list, DATASET *dset, 
			gretlopt opt, PRN *prn);

double get_urc_pvalue (double tau, int n, int niv, int itv);

gretl_matrix *kpss_critvals (int T, int trend, int *err);

#endif /* ADF_KPSS_H */
