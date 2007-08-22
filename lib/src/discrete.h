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

#ifndef DISCRETE_H
#define DISCRETE_H

MODEL logit_probit (const int *list, double ***pZ, DATAINFO *pdinfo, 
		    int ci, gretlopt opt, PRN *prn);

MODEL logistic_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		      const char *param); 

int logistic_ymax_lmax (const double *y, const DATAINFO *pdinfo,
			double *ymax, double *lmax);

int fishers_exact_test (const Xtab *tab, PRN *prn);

#endif /* DISCRETE_H */



