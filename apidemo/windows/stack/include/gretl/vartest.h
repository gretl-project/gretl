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

#ifndef VARTEST_H_
#define VARTEST_H_

#define N_IVALS 3
#define VAR_S_DFCORR 0

void gretl_VAR_clear (GRETL_VAR *var);

void VAR_fill_X (GRETL_VAR *v, int p, const DATASET *dset);

void VAR_write_A_matrix (GRETL_VAR *v, GretlMatrixMod mod);

gretl_matrix *decompanionize (const gretl_matrix *A, int neqns,
			      GretlMatrixMod mod);

int johansen_stage_1 (GRETL_VAR *jvar, const DATASET *dset,
		      gretlopt opt, PRN *prn);

double gretl_VAR_ldet (GRETL_VAR *var, const gretl_matrix *E,
		       int *err);

int VAR_LR_lag_test (GRETL_VAR *var, const gretl_matrix *E);

int VAR_portmanteau_test (GRETL_VAR *var);

int VAR_do_lagsel (GRETL_VAR *var, const DATASET *dset, 
		   gretlopt opt, PRN *prn,
		   gretl_matrix **pm);

int VAR_wald_omit_tests (GRETL_VAR *var);

gretl_matrix *VAR_coeff_matrix_from_VECM (GRETL_VAR *var,
					  int Atrans);

gretl_matrix *reorder_responses (const GRETL_VAR *var, int *err);

void copy_north_west (gretl_matrix *targ,
		      const gretl_matrix *src,
		      int add);

gretl_matrix *irf_bootstrap (GRETL_VAR *var, 
			     int targ, int shock,
			     int periods, double alpha,
			     gretl_matrix *point,
			     const DATASET *dset,
			     int *err);

#endif /* VARTEST_H_ */
