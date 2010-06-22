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

MODEL binary_logit (int *list, double ***pZ, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn);

MODEL binary_probit (int *list, double ***pZ, DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn);

MODEL ordered_logit (int *list, double ***pZ, DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn);

MODEL ordered_probit (int *list, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn);

MODEL multinomial_logit (int *list, double ***pZ, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn);

MODEL logistic_model (const int *list, double lmax,
		      double ***pZ, DATAINFO *pdinfo);

MODEL interval_model (int *list, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt, PRN *prn);

MODEL tobit_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn);

MODEL duration_model (const int *list, double ***pZ, 
		      DATAINFO *pdinfo, gretlopt opt, 
		      PRN *prn);

MODEL count_model (const int *list, int ci,
		   double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, PRN *prn);

MODEL heckit_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn);

int fishers_exact_test (const Xtab *tab, PRN *prn);

double ordered_model_prediction (const MODEL *pmod, double Xb);

int logistic_ymax_lmax (const double *y, const DATAINFO *pdinfo,
			double *ymax, double *lmax);

#endif /* DISCRETE_H */



