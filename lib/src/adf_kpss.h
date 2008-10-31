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

typedef enum {
    UR_NO_CONST = 1,
    UR_CONST,
    UR_TREND,
    UR_QUAD_TREND,
    UR_MAX
} AdfCode;

int adf_test (int order, const int *list, double ***pZ,
	      DATAINFO *pdinfo, gretlopt opt, PRN *prn);

int kpss_test (int order, const int *list, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn);

int coint (int order, const int *list, double ***pZ, 
	   DATAINFO *pdinfo, gretlopt opt, PRN *prn);

double df_pvalue_from_plugin (double tau, int n, int niv, 
			      int itv, gretlopt opt);

#endif /* ADF_KPSS_H */
