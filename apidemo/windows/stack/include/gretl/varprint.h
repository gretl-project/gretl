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

#ifndef VARPRINT_H_
#define VARPRINT_H_

int gretl_VAR_print_sigma (const GRETL_VAR *var, PRN *prn);

int gretl_VAR_print (GRETL_VAR *var, const DATASET *dset, gretlopt opt,
		     PRN *prn);

int gretl_VAR_print_all_fcast_decomps (GRETL_VAR *var, const DATASET *dset, 
				       int horizon, PRN *prn);

int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATASET *dset, 
			      PRN *prn);

int gretl_VAR_print_all_impulse_responses (GRETL_VAR *var, const DATASET *dset, 
					   int horizon, PRN *prn);

char *vecm_beta_varname (char *vname,
			 const GRETL_VAR *v,
			 const DATASET *dset,
			 int i);

#endif /* VARPRINT_H_ */
