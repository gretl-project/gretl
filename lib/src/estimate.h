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
 
MODEL lsq (LIST list, 
	   double ***pZ, DATAINFO *pdinfo, 
	   int ci, int opt, double rho);

int hilu_corc (double *toprho, LIST list, 
	       double ***pZ, DATAINFO *pdinfo, 
	       int opt, PRN *prn);

MODEL lad (LIST list, double ***pZ, DATAINFO *pdinfo); 

MODEL tsls_func (LIST list, const int pos, 
		 double ***pZ, DATAINFO *pdinfo);

MODEL hsk_func (LIST list, 
		double ***pZ, DATAINFO *pdinfo);

MODEL hccm_func (LIST list, 
		 double ***pZ, DATAINFO *pdinfo);

int whites_test (MODEL *pmod, 
		 double ***pZ, DATAINFO *pdinfo, 
		 PRN *prn, GRETLTEST *test);

MODEL ar_func (LIST list, int pos, 
	       double ***pZ, DATAINFO *pdinfo, 
	       int *model_count, PRN *prn);

MODEL arch (int order, LIST list, 
	    double ***pZ, DATAINFO *pdinfo, 
	    int *model_count, PRN *prn, 
	    GRETLTEST *test);

int makevcv (MODEL *pmod);
