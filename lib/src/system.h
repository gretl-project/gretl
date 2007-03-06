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

enum {
    GRETL_SYSTEM_SAVE_UHAT = 1 << 0,
    GRETL_SYSTEM_SAVE_YHAT = 1 << 1,
    GRETL_SYSTEM_DFCORR    = 1 << 2,
    GRETL_SYS_VCV_GEOMEAN  = 1 << 3,
    GRETL_SYS_SAVE_VCV     = 1 << 4,
    GRETL_SYS_RESTRICT     = 1 << 5,
    GRETL_SYS_ITERATE      = 1 << 6,
    GRETL_SYS_SAVEIT       = 1 << 7
};

typedef struct id_atom_ id_atom;
typedef struct identity_ identity;

struct gretl_equation_system_ {
    char *name;                 /* user-specified name for system, or NULL */
    int refcount;               /* for saving/deleting */
    int t1;                     /* starting observation number */
    int t2;                     /* ending observation number */
    int method;                 /* estimation method */
    int n_equations;            /* number of stochastic equations */
    int n_identities;           /* number of identities */
    int n_obs;                  /* number of observations per equation */
    int iters;                  /* number of iterations taken */
    char flags;                 /* to record options (e.g. save residuals) */
    double ll;                  /* log-likelihood (restricted) */
    double llu;                 /* unrestricted log-likelihood */
    double X2;                  /* chi-square test value */
    double ess;                 /* total error sum of squares */
    double diag;                /* test stat for diagonal covariance matrix */
    double bdiff;               /* summary stat for change in coefficients */
    int **lists;                /* regression lists for stochastic equations */
    int *endog_vars;            /* list of endogenous variables */
    int *instr_vars;            /* list of instruments (exogenous vars) */
    identity **idents;          /* set of identities */
    gretl_matrix *b;            /* coefficient estimates */
    gretl_matrix *vcv;          /* covariance matrix of coefficients */
    gretl_matrix *sigma;        /* cross-equation covariance matrix */
    gretl_matrix *R;            /* LHS of any linear restrictions */
    gretl_matrix *q;            /* RHS of any linear restrictions */  
    gretl_matrix *uhat;         /* residuals, all equations */
    MODEL **models;             /* set of pointers to per-equation models: just
				   convenience pointers -- these should NOT be
				   freed as part of sys cleanup
				*/
};

void make_system_data_info (gretl_equation_system *sys, int eqn, 
			    DATAINFO *pdinfo, int v, int code);

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

void gretl_equation_system_destroy (gretl_equation_system *sys);

const char *system_get_full_string (const gretl_equation_system *sys);

int system_save_uhat (const gretl_equation_system *sys);
int system_save_yhat (const gretl_equation_system *sys);

int system_save_vcv (const gretl_equation_system *sys);

int system_want_df_corr (const gretl_equation_system *sys);

int system_n_restrictions (const gretl_equation_system *sys);

int system_max_indep_vars (const gretl_equation_system *sys);
int system_n_indep_vars (const gretl_equation_system *sys);

int system_adjust_t1t2 (gretl_equation_system *sys,
			int *t1, int *t2, const double **Z);

int *system_get_list (const gretl_equation_system *sys, int i);

int *compose_tsls_list (gretl_equation_system *sys, int i);

int system_get_depvar (const gretl_equation_system *sys, int i);

const char *gretl_system_short_string (const MODEL *pmod);

void gretl_system_set_name (gretl_equation_system *sys, const char *name);

int gretl_system_method_from_string (const char *s);
const char *system_method_full_string (int method);
const char *system_method_short_string (int method);

int *system_get_endog_vars (const gretl_equation_system *sys);
int *system_get_instr_vars (const gretl_equation_system *sys);

void system_attach_uhat (gretl_equation_system *sys, gretl_matrix *uhat);

void system_attach_sigma (gretl_equation_system *sys, gretl_matrix *sigma);

void system_attach_coeffs (gretl_equation_system *sys, gretl_matrix *b);
void system_attach_vcv (gretl_equation_system *sys, gretl_matrix *vcv);

MODEL *system_get_model (const gretl_equation_system *sys, int i);

int system_get_overid_df (const gretl_equation_system *sys);

int system_vcv_geomean (const gretl_equation_system *sys);
double system_vcv_denom (const gretl_equation_system *sys, 
			 int i, int j);

int rhs_var_in_identity (const gretl_equation_system *sys, int lhsvar,
			 int rhsvar);

void 
print_equation_system_info (const gretl_equation_system *sys, 
			    const DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn);

void 
system_set_restriction_matrices (gretl_equation_system *sys,
				 gretl_matrix *R, gretl_matrix *q);

int 
system_normality_test (const gretl_equation_system *sys, PRN *prn);

int gretl_system_add_resids_to_dataset (gretl_equation_system *sys, int eqnum,
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

int gretl_system_serialize (gretl_equation_system *sys, 
			    SavedObjectFlags flags,
			    FILE *fp);

void gretl_system_set_save_flag (gretl_equation_system *sys);

void gretl_system_unset_save_flag (gretl_equation_system *sys);

int gretl_system_save_flag_set (gretl_equation_system *sys);

#ifndef GRETLCLI

gretl_equation_system *
gretl_system_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err);

#endif


#endif /* GRETL_EQUATION_SYSTEM_H */
