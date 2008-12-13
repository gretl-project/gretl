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

MODEL lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	   GretlCmdIndex ci, gretlopt opt);

MODEL ar1_lsq (const int *list, double ***pZ, DATAINFO *pdinfo, 
	    GretlCmdIndex ci, gretlopt opt, double rho);

double estimate_rho (const int *list, double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn, int *err);

MODEL lad (const int *list, double ***pZ, DATAINFO *pdinfo); 

MODEL quantreg (const char *parm, const int *list, 
		double ***pZ, DATAINFO *pdinfo,
		gretlopt opt, PRN *prn);

MODEL arma (const int *list, const char *pqspec,
	    const double **Z, const DATAINFO *pdinfo, 
	    gretlopt opt, PRN *prn);

MODEL tobit_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
		   PRN *prn);

MODEL poisson_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
		     PRN *prn);

MODEL heckit_model (const int *list, double ***pZ, DATAINFO *pdinfo, 
		    gretlopt opt, PRN *prn);

MODEL garch (const int *list, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	     PRN *prn);

MODEL mp_ols (const int *list, const double **Z, DATAINFO *pdinfo);

MODEL panel_model (const int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn);

MODEL ivreg (const int *list, double ***pZ, DATAINFO *pdinfo,
	     gretlopt opt);

MODEL arbond_model (const int *list, const char *istr, const double **Z, 
		    const DATAINFO *pdinfo, gretlopt opt, 
		    PRN *prn);

int groupwise_hetero_test (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
			   PRN *prn);

MODEL hsk_func (const int *list, double ***pZ, DATAINFO *pdinfo);

int whites_test (MODEL *pmod, 
		 double ***pZ, DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn);

MODEL ar_func (const int *list,  
	       double ***pZ, DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int arch_test (MODEL *pmod, int order, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int array_arch_test (const double *u, int n, int order, 
		     gretlopt opt, PRN *prn);

MODEL arch_model (const int *list, int order, 
		  double ***pZ, DATAINFO *pdinfo, 
		  gretlopt opt, PRN *prn);

int makevcv (MODEL *pmod, double sigma);

int *augment_regression_list (const int *orig, int aux, 
			      double ***pZ, DATAINFO *pdinfo);

int gretl_XTX_XTy (const int *list, int t1, int t2, 
		   const double **Z, int nwt, double rho, int pwe,
		   double *xpx, double *xpy, const char *mask);

#endif /* ESTIMATE_H */


