/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef VAR_H_
#define VAR_H_

#include "gretl_matrix.h"

int list_diffgenr (const LIST list, 
		   double ***pZ, DATAINFO *pdinfo);

int list_ldiffgenr (const LIST list, 
		    double ***pZ, DATAINFO *pdinfo);

int simple_var (int order, const LIST list, 
		double ***pZ, DATAINFO *pdinfo,
		int pause, PRN *prn);

GRETL_VAR *full_var (int order, const LIST list, 
		     double ***pZ, DATAINFO *pdinfo,
		     PRN *prn);

void gretl_var_free (GRETL_VAR *var);

void gretl_var_free_unnamed (GRETL_VAR *var);

int coint (int order, const LIST list, 
	   double ***pZ, DATAINFO *pdinfo, 
	   PRN *prn);

int adf_test (int order, const int varno, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn);

int ma_model (LIST list, double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn);

int johansen_test (int order, const LIST list, 
		   double ***pZ, DATAINFO *pdinfo,
		   int verbose, PRN *prn);

double *
gretl_var_get_impulse_responses (GRETL_VAR *var, int targ, int shock,
				 int periods);

int 
gretl_var_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  PRN *prn);

int gretl_var_print (GRETL_VAR *var, const DATAINFO *pdinfo, PRN *prn);

void gretl_var_assign_name (GRETL_VAR *var);

void gretl_var_assign_specific_name (GRETL_VAR *var, const char *name);

const char *gretl_var_get_name (const GRETL_VAR *var);

int gretl_var_get_variable_number (const GRETL_VAR *var, int k);

int gretl_var_get_n_equations (const GRETL_VAR *var);

#endif /* VAR_H_ */

