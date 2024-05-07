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

MODEL binary_logit (const int *list, DATASET *dset, 
		    gretlopt opt, PRN *prn);

MODEL binary_probit (const int *list, DATASET *dset, 
		     gretlopt opt, PRN *prn);

MODEL ordered_logit (int *list, DATASET *dset, 
		     gretlopt opt, PRN *prn);

MODEL ordered_probit (int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn);

MODEL multinomial_logit (int *list, DATASET *dset, 
			 gretlopt opt, PRN *prn);

MODEL biprobit_model (int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn);

MODEL reprobit_model (const int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn);

MODEL logistic_model (const int *list, double lmax,
		      DATASET *dset, gretlopt opt);

MODEL interval_model (int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn);

MODEL tobit_model (const int *list, double llim, double rlim,
		   DATASET *dset, gretlopt opt, PRN *prn);

MODEL duration_model (const int *list, DATASET *dset, 
		      gretlopt opt, PRN *prn);

MODEL count_model (const int *list, int ci, DATASET *dset, 
		   gretlopt opt, PRN *prn);

MODEL heckit_model (const int *list, DATASET *dset, 
		    gretlopt opt, PRN *prn);

int fishers_exact_test (const Xtab *tab, PRN *prn);

double ordered_model_prediction (const MODEL *pmod, double Xb,
				 int ymin);

int logistic_ymax_lmax (const double *y, const DATASET *dset,
			double *ymax, double *lmax);

gretl_matrix *mn_logit_probabilities (const MODEL *pmod,
				      int t1, int t2,
				      const DATASET *dset,
				      int *err);

gretl_matrix *ordered_probabilities (const MODEL *pmod,
				     const double *zhat,
				     int t1, int t2,
				     const DATASET *dset,
				     int *err);

double mn_logit_prediction (const gretl_matrix *Xt,
			    const double *b,
			    const gretl_matrix *yvals);

void binary_model_hatvars (MODEL *pmod, 
			   const gretl_matrix *ndx,
			   const int *y,
			   gretlopt opt);

gretl_matrix *bit_permutations (int n, int k, int *err);

#endif /* DISCRETE_H */



