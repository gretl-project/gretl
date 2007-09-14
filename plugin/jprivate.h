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

#ifndef JPRIVATE_H
#define JPRIVATE_H

enum {
    V_ALPHA,
    V_BETA
};

int johansen_ll_calc (GRETL_VAR *jvar, const gretl_matrix *evals);

int
johansen_LR_calc (const GRETL_VAR *jvar, const gretl_matrix *evals, 
		  const gretl_matrix *H, int job, PRN *prn);

void print_beta_alpha_Pi (const GRETL_VAR *jvar,
			  const DATAINFO *pdinfo,
			  PRN *prn);

int 
general_vecm_analysis (GRETL_VAR *jvar, 
		       const gretl_restriction *rset,
		       const DATAINFO *pdinfo,
		       PRN *prn);

int vecm_alpha_test (GRETL_VAR *jvar, 
		     const gretl_restriction *rset,
		     const DATAINFO *pdinfo, 
		     gretlopt opt,
		     PRN *prn);

const char *beta_vname (const GRETL_VAR *v,
			const DATAINFO *pdinfo,
			int i);

#endif /* JPRIVATE_H */
