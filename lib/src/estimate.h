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

/* estimate.h for gretl:  set out structs and functions employed */

#ifndef ESTIMATE_H
#define ESTIMATE_H

#include "gretl_matrix.h"

MODEL lsq (const int *list, double **Z, DATAINFO *pdinfo, 
	   GretlCmdIndex ci, gretlopt opt);

MODEL ar_model (const int *list,  
		double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn);

MODEL ar1_model (const int *list,  
		 double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn);

MODEL lad (const int *list, double **Z, DATAINFO *pdinfo); 

MODEL quantreg (const gretl_matrix *tau, const int *list, 
		double **Z, DATAINFO *pdinfo,
		gretlopt opt, PRN *prn);

MODEL arma (const int *list, const int *pqlags,
	    const double **Z, const DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn);

MODEL garch (const int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn);

MODEL mp_ols (const int *list, const double **Z, DATAINFO *pdinfo);

MODEL panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn);

MODEL ivreg (const int *list, double ***pZ, DATAINFO *pdinfo,
	     gretlopt opt);

MODEL arbond_model (const int *list, const char *ispec, const double **Z, 
		    const DATAINFO *pdinfo, gretlopt opt, PRN *prn);

MODEL dpd_model (const int *list, const int *laglist,
		 const char *ispec, const double **Z, 
		 const DATAINFO *pdinfo, gretlopt opt, 
		 PRN *prn);

MODEL hsk_model (const int *list, double ***pZ, DATAINFO *pdinfo);

MODEL arch_model (const int *list, int order, 
		  double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn);

int whites_test (MODEL *pmod, 
		 double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn);

int arch_test (MODEL *pmod, int order, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int array_arch_test (const double *u, int n, int order, 
		     gretlopt opt, PRN *prn);

int makevcv (MODEL *pmod, double sigma);

int *augment_regression_list (const int *orig, int aux, 
			      double ***pZ, DATAINFO *pdinfo);

double *gretl_XTX (const MODEL *pmod, const double **Z, int *err);

int anova (const int *list, const double **Z, const DATAINFO *pdinfo, 
	   gretlopt opt, PRN *prn);

#endif /* ESTIMATE_H */


