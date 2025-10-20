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

MODEL lsq (const int *list, DATASET *dset, 
	   GretlCmdIndex ci, gretlopt opt);

MODEL ar_model (const int *list, DATASET *dset, 
		gretlopt opt, PRN *prn);

MODEL ar1_model (const int *list, DATASET *dset, 
		 gretlopt opt, PRN *prn);

MODEL lad_model (const int *list, DATASET *dset, gretlopt opt);

MODEL fols_model (const int *list, DATASET *dset,
                  gretlopt opt, PRN *prn);

MODEL quantreg (const gretl_matrix *tau, const int *list, 
		DATASET *dset, gretlopt opt, PRN *prn);

MODEL arma (const int *list, const int *pqlags,
	    DATASET *dset, gretlopt opt, PRN *prn);

MODEL garch (const int *list, DATASET *dset, gretlopt opt,
	     PRN *prn);

MODEL mp_ols (const int *list, DATASET *dset, gretlopt opt);

MODEL panel_model (const int *list, DATASET *dset,
		   gretlopt opt, PRN *prn);

MODEL ivreg (const int *list, DATASET *dset, gretlopt opt);

MODEL dpd_model (const int *list, const int *laglist,
		 const char *ispec, const DATASET *dset, 
		 gretlopt opt, PRN *prn);

MODEL hsk_model (const int *list, DATASET *dset, gretlopt opt);

MODEL arch_model (const int *list, int order, DATASET *dset, 
		  gretlopt opt);

int whites_test (MODEL *pmod, DATASET *dset, 
		 gretlopt opt, PRN *prn);

int arch_test (MODEL *pmod, int order, const DATASET *dset,
	       gretlopt opt, PRN *prn);

int array_arch_test (const double *u, int n, int order, 
		     gretlopt opt, PRN *prn);

int makevcv (MODEL *pmod, double sigma);

int *augment_regression_list (const int *orig, int aux, 
			      DATASET *dset, int *err);

int anova (const int *list, const DATASET *dset, 
	   gretlopt opt, PRN *prn);

#endif /* ESTIMATE_H */


