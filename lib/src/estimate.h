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

/* estimate.h for gretl:  set out structs and functions employed */

#include <stdio.h>

/* functions follow */

MODEL lsq (int *list, 
	   double ***pZ, DATAINFO *pdinfo, 
	   int ci, gretlopt opts, double rho);

double estimate_rho (int *list, double ***pZ, DATAINFO *pdinfo,
		     int batch, int opt, int *err, PRN *prn);

int hilu_corc (double *toprho, int *list, 
	       double ***pZ, DATAINFO *pdinfo,
	       PATHS *ppaths, int batch,
	       int opt, PRN *prn);

MODEL lad (int *list, double ***pZ, DATAINFO *pdinfo); 

MODEL arma (int *list, const double **Z, DATAINFO *pdinfo, 
	    PRN *prn);

MODEL arma_x12 (int *list, const double **Z, DATAINFO *pdinfo, 
		PRN *prn, const PATHS *ppaths);

MODEL logistic_model (int *list, double ***pZ, DATAINFO *pdinfo,
		      const char *param);

MODEL tobit_model (int *list, double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn);

MODEL garch (int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn);

MODEL pooled (int *list, double ***pZ, DATAINFO *pdinfo,
	      gretlopt opt, PRN *prn);

int groupwise_hetero_test (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			   PRN *prn);

const double *tsls_get_Xi (const MODEL *pmod, const double **Z, int i);

void tsls_free_data (const MODEL *pmod);

MODEL tsls_func (int *list, int pos_in, 
		 double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt);

MODEL hsk_func (int *list, double ***pZ, DATAINFO *pdinfo);

int whites_test (MODEL *pmod, 
		 double ***pZ, DATAINFO *pdinfo, 
		 PRN *prn, GRETLTEST *test);

MODEL ar_func (int *list, int pos, 
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

MODEL arch (int order, int *list, 
	    double ***pZ, DATAINFO *pdinfo, 
	    GRETLTEST *test, gretlopt opt, PRN *prn);

int makevcv (MODEL *pmod);

VCV *get_vcv (MODEL *pmod);

void free_vcv (VCV *vcv);

