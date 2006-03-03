/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_EQUATION_SYSTEM_H
#define GRETL_EQUATION_SYSTEM_H

typedef struct _gretl_equation_system gretl_equation_system;

typedef enum {
    SYS_SUR = 0,
    SYS_3SLS,
    SYS_FIML,
    SYS_LIML,
    SYS_OLS,
    SYS_TSLS,
    SYS_WLS,
    SYS_MAX
} gretl_system_methods;

gretl_equation_system *system_start (const char *line, gretlopt opt);

char *get_system_name_from_line (const char *s);

int gretl_equation_system_append (gretl_equation_system *sys, 
				  const int *list);

int system_parse_line (gretl_equation_system *sys,
		       const char *line, 
		       const DATAINFO *pdinfo);

int gretl_equation_system_finalize (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn);

int 
gretl_equation_system_estimate (gretl_equation_system *sys, 
				double ***pZ, DATAINFO *pdinfo, 
				gretlopt opt, PRN *prn);

int estimate_named_system (const char *line, double ***pZ, DATAINFO *pdinfo, 
			   gretlopt opt, PRN *prn);

int estimate_saved_equation_system (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn);

void gretl_equation_system_destroy (gretl_equation_system *sys);

const char *system_get_full_string (const gretl_equation_system *sys);

int system_save_uhat (const gretl_equation_system *sys);
int system_save_yhat (const gretl_equation_system *sys);

int system_save_vcv (const gretl_equation_system *sys);

int system_want_df_corr (const gretl_equation_system *sys);

int system_n_equations (const gretl_equation_system *sys);
int system_n_identities (const gretl_equation_system *sys);
int system_n_restrictions (const gretl_equation_system *sys);

int system_n_obs (const gretl_equation_system *sys);
void system_set_n_obs (gretl_equation_system *sys, int n);

int system_iters (const gretl_equation_system *sys);
void system_set_iters (gretl_equation_system *sys, int n);

int system_max_indep_vars (const gretl_equation_system *sys);
int system_n_indep_vars (const gretl_equation_system *sys);

int system_adjust_t1t2 (gretl_equation_system *sys,
			int *t1, int *t2, const double **Z);

int *system_get_list (const gretl_equation_system *sys, int i);

int *compose_tsls_list (gretl_equation_system *sys, int i);

int system_get_depvar (const gretl_equation_system *sys, int i);

const char *gretl_system_short_string (const MODEL *pmod);

void gretl_system_set_name (gretl_equation_system *sys, const char *name);
const char *gretl_system_get_name (const gretl_equation_system *sys);

const gretl_matrix *system_get_R_matrix (const gretl_equation_system *sys);
const gretl_matrix *system_get_q_matrix (const gretl_equation_system *sys);

int system_get_method (const gretl_equation_system *sys);

int gretl_system_method_from_string (const char *s);
const char *system_method_full_string (int method);
const char *system_method_short_string (int method);

int *system_get_endog_vars (const gretl_equation_system *sys);
int *system_get_instr_vars (const gretl_equation_system *sys);

void system_attach_uhat (gretl_equation_system *sys, gretl_matrix *uhat);
gretl_matrix *system_get_uhat (const gretl_equation_system *sys);

void system_attach_sigma (gretl_equation_system *sys, gretl_matrix *sigma);
gretl_matrix *system_get_sigma (const gretl_equation_system *sys);

void system_attach_coeffs (gretl_equation_system *sys, gretl_matrix *b);
void system_attach_vcv (gretl_equation_system *sys, gretl_matrix *vcv);

MODEL *system_get_model (const gretl_equation_system *sys, int i);

double system_get_ll (const gretl_equation_system *sys);
double system_get_llu (const gretl_equation_system *sys);
double system_get_X2 (const gretl_equation_system *sys);
double system_get_ess (const gretl_equation_system *sys);
double system_get_diag_stat (const gretl_equation_system *sys);

void system_set_ll (gretl_equation_system *sys, double ll);
void system_set_llu (gretl_equation_system *sys, double llu);
void system_set_X2 (gretl_equation_system *sys, double X2);
void system_set_ess (gretl_equation_system *sys, double ess);
void system_set_diag_stat (gretl_equation_system *sys, double s);

int system_get_overid_df (const gretl_equation_system *sys);

int system_vcv_geomean (const gretl_equation_system *sys);
double system_vcv_denom (const gretl_equation_system *sys, 
			 int i, int j);

int rhs_var_in_identity (const gretl_equation_system *sys, int lhsvar,
			 int rhsvar);

void 
print_equation_system_info (const gretl_equation_system *sys, 
			    const DATAINFO *pdinfo, PRN *prn);

void 
system_set_restriction_matrices (gretl_equation_system *sys,
				 gretl_matrix *R, gretl_matrix *q);

int 
system_normality_test (const gretl_equation_system *sys, PRN *prn);

int gretl_system_add_resids_to_dataset (const char *sysname, int eqnum,
					double ***pZ, DATAINFO *pdinfo);

double *
gretl_equation_system_get_series (const gretl_equation_system *sys, 
				  const DATAINFO *pdinfo,
				  int idx, const char *key, int *err);

gretl_matrix *
gretl_equation_system_get_matrix (const gretl_equation_system *sys, int idx, 
				  int *err);

int highest_numbered_var_in_system (const gretl_equation_system *sys, 
				    const DATAINFO *pdinfo);

#endif /* GRETL_EQUATION_SYSTEM_H */
