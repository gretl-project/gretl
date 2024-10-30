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

#ifndef GRETL_RESTRICT_H
#define GRETL_RESTRICT_H

typedef struct gretl_restriction_ gretl_restriction;

gretl_restriction *
restriction_set_start (const char *target, gretlopt opt, int *err);

gretl_restriction *
cross_restriction_set_start (const char *line, equation_system *sys);

gretl_restriction *
var_restriction_set_start (const char *line, GRETL_VAR *var);

gretl_restriction *
eqn_restriction_set_start (const char *line, MODEL *pmod, 
			   const DATASET *dset,
			   gretlopt opt);

gretl_restriction *rset_from_VECM (GRETL_VAR *var, int *err);

int 
restriction_set_parse_line (gretl_restriction *rset, const char *line,
			    const DATASET *dset);

int
gretl_restriction_finalize (gretl_restriction *rset, 
			    const DATASET *dset,
			    gretlopt opt, PRN *prn);

int
gretl_restriction_finalize_full (ExecState *state,
				 gretl_restriction *rset, 
				 const DATASET *dset,
				 gretlopt opt,
				 PRN *prn);

GRETL_VAR *
gretl_restricted_vecm (gretl_restriction *rset, 
		       const DATASET *dset,
		       gretlopt opt, 
		       PRN *prn,
		       int *err);

void print_restriction_from_matrices (const gretl_matrix *R,
				      const gretl_matrix *q,
				      char letter, int npar, 
				      PRN *prn);

void destroy_restriction_set (gretl_restriction *rset);

int gretl_sum_test (const int *list, MODEL *pmod, DATASET *dset,
		    gretlopt opt, PRN *prn);

const gretl_matrix *rset_get_R_matrix (const gretl_restriction *rset);

const gretl_matrix *rset_get_q_matrix (const gretl_restriction *rset);

const gretl_matrix *rset_get_Ra_matrix (const gretl_restriction *rset);

const gretl_matrix *rset_get_qa_matrix (const gretl_restriction *rset);

int gretl_restriction_set_boot_params (int B, gretlopt opt);

void gretl_restriction_get_boot_params (int *pB, gretlopt *popt);

gretlopt gretl_restriction_get_options (const gretl_restriction *rset);

GretlObjType gretl_restriction_get_type (const gretl_restriction *rset);

int rset_VECM_bcols (const gretl_restriction *rset);

int rset_VECM_acols (const gretl_restriction *rset);

void rset_add_results (gretl_restriction *rset,
		       double test, double pval,
		       double lnl);

void rset_record_LR_result (gretl_restriction *rset);

#endif /* GRETL_RESTRICT_H */


