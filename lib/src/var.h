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

typedef enum {
    UR_NO_CONST = 1,
    UR_CONST,
    UR_TREND,
    UR_TREND_SQUARED
} AdfCode;

typedef enum {
    J_NO_CONST = 0,
    J_REST_CONST,
    J_UNREST_CONST,
    J_REST_TREND,
    J_UNREST_TREND
} JohansenCode;

int var_max_order (const int *list, const DATAINFO *pdinfo);

int simple_VAR (int order, const int *list, 
		double ***pZ, DATAINFO *pdinfo,
		gretlopt opt, PRN *prn);

GRETL_VAR *full_VAR (int order, const int *list, 
		     double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn);

const gretl_matrix *
gretl_VAR_get_forecast_matrix (GRETL_VAR *var, int t1, int t2,
			       const double **Z, 
			       const DATAINFO *pdinfo,
			       gretlopt opt);

const gretl_matrix *
gretl_VAR_get_residual_matrix (const GRETL_VAR *var);

int gretl_VAR_print_VCV (const GRETL_VAR *var, PRN *prn);

int gretl_VAR_autocorrelation_test (GRETL_VAR *var, int order, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn);

int gretl_VAR_arch_test (GRETL_VAR *var, int order, 
			 double ***pZ, DATAINFO *pdinfo,
			 PRN *prn);

int gretl_VAR_normality_test (const GRETL_VAR *var, PRN *prn);

void gretl_VAR_free (GRETL_VAR *var);

void gretl_VAR_free_unnamed (GRETL_VAR *var);

int coint (int order, const int *list, 
	   double ***pZ, DATAINFO *pdinfo, 
	   gretlopt opt, PRN *prn);

int adf_test (int order, int varno, 
	      double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn);

int kpss_test (int order, int varno, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn);

int johansen_test (int order, const int *list, 
		   double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn);

int default_VAR_horizon (const DATAINFO *pdinfo);

gretl_matrix *
gretl_VAR_get_impulse_response (GRETL_VAR *var, 
				int targ, int shock,
				int periods,
				const double **Z,
				const DATAINFO *pdinfo);

int gretl_VAR_print (GRETL_VAR *var, const DATAINFO *pdinfo, gretlopt opt,
		     PRN *prn);

int 
gretl_VAR_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  int pause, PRN *prn);

int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      int pause, PRN *prn);

void gretl_VAR_assign_name (GRETL_VAR *var);

void gretl_VAR_assign_specific_name (GRETL_VAR *var, const char *name);

const char *gretl_VAR_get_name (const GRETL_VAR *var);

int gretl_VAR_get_variable_number (const GRETL_VAR *var, int k);

int gretl_VAR_get_n_equations (const GRETL_VAR *var);

int gretl_VAR_get_t1 (const GRETL_VAR *var);

int gretl_VAR_get_t2 (const GRETL_VAR *var);

const MODEL *gretl_VAR_get_model (const GRETL_VAR *var, int i);

int gretl_VAR_add_resids_to_dataset (GRETL_VAR *var, int eqnum,
				     double ***pZ, DATAINFO *pdinfo);

int gretl_VAR_get_highest_variable (const GRETL_VAR *var,
				    const DATAINFO *pdinfo);

#endif /* VAR_H_ */

