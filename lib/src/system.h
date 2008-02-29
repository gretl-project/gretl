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

#ifndef EQUATION_SYSTEM_H
#define EQUATION_SYSTEM_H

typedef enum {
    SYS_METHOD_SUR = 0,
    SYS_METHOD_3SLS,
    SYS_METHOD_FIML,
    SYS_METHOD_LIML,
    SYS_METHOD_OLS,
    SYS_METHOD_TSLS,
    SYS_METHOD_WLS,
    SYS_METHOD_MAX
} GretlSystemMethods;

enum {
    SYSTEM_SAVE_UHAT   = 1 << 0,
    SYSTEM_SAVE_YHAT   = 1 << 1,
    SYSTEM_DFCORR      = 1 << 2,
    SYSTEM_VCV_GEOMEAN = 1 << 3,
    SYSTEM_RESTRICT    = 1 << 4,
    SYSTEM_ITERATE     = 1 << 5,
    SYSTEM_SAVEIT      = 1 << 6
};

typedef struct id_atom_ id_atom;
typedef struct identity_ identity;
typedef struct predet_ predet;

struct equation_system_ {
    char *name;                 /* user-specified name for system, or NULL */
    int refcount;               /* for saving/deleting */
    int t1;                     /* starting observation number */
    int t2;                     /* ending observation number */
    int method;                 /* estimation method */
    int neqns;                  /* number of stochastic equations */
    int nidents;                /* number of identities */
    int n_obs;                  /* number of observations per equation */
    int n_predet;               /* number of predetermined regressors */
    int maxlag;                 /* max lag of endogenous variable */
    int iters;                  /* number of iterations taken */
    char flags;                 /* to record options (e.g. save residuals) */
    double ll;                  /* log-likelihood (restricted) */
    double llu;                 /* unrestricted log-likelihood */
    double X2;                  /* chi-square test value */
    double ess;                 /* total error sum of squares */
    double diag;                /* test stat for diagonal covariance matrix */
    double bdiff;               /* summary stat for change in coefficients */
    int **lists;                /* regression lists for stochastic equations */
    int *elist;                 /* list of endogenous variables */
    int *ilist;                 /* list of instruments */
    int *xlist;                 /* list of truly exogenous variables */
    predet *pre_vars;           /* array of info on predetermined regressors */
    identity **idents;          /* set of identities */
    gretl_matrix *b;            /* coefficient estimates */
    gretl_matrix *vcv;          /* covariance matrix of coefficients */
    gretl_matrix *sigma;        /* cross-equation covariance matrix */
    gretl_matrix *R;            /* LHS of any linear restrictions */
    gretl_matrix *q;            /* RHS of any linear restrictions */  
    gretl_matrix *uhat;         /* residuals, all equations */
    gretl_matrix *yhat;         /* fitted values, all equations */
    gretl_matrix *Gamma;        /* structural form Gamma matrix (endog + identities)*/
    gretl_matrix *B;            /* structural form B matrix (exogenous) */
    gretl_matrix *A;            /* structural form A matrix (lagged endogenous) */
    gretl_matrix *F;            /* forecast matrix */
    MODEL **models;             /* set of pointers to per-equation models: just
				   convenience pointers -- these should NOT be
				   freed as part of sys cleanup
				*/
};

equation_system *equation_system_start (const char *line, 
					gretlopt opt);

char *get_system_name_from_line (const char *s);

int equation_system_append (equation_system *sys, 
			    const int *list);

int system_parse_line (equation_system *sys,
		       const char *line, 
		       const DATAINFO *pdinfo);

int equation_system_finalize (equation_system *sys, 
			      double ***pZ, DATAINFO *pdinfo,
			      PRN *prn);

int 
equation_system_estimate (equation_system *sys, 
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, PRN *prn);

int estimate_named_system (const char *line, double ***pZ, DATAINFO *pdinfo, 
			   gretlopt opt, PRN *prn);

void equation_system_destroy (equation_system *sys);

const char *system_get_full_string (const equation_system *sys);

int system_want_df_corr (const equation_system *sys);

int system_n_restrictions (const equation_system *sys);

int system_max_indep_vars (const equation_system *sys);
int system_n_indep_vars (const equation_system *sys);

int system_adjust_t1t2 (equation_system *sys,
			int *t1, int *t2, const double **Z);

int *system_get_list (const equation_system *sys, int i);

int system_get_list_length (const equation_system *sys, int i);

int *compose_tsls_list (equation_system *sys, int i);

int system_get_depvar (const equation_system *sys, int i);

const char *system_short_string (const MODEL *pmod);

void equation_system_set_name (equation_system *sys, const char *name);

int system_method_from_string (const char *s);
const char *system_method_full_string (int method);
const char *system_method_short_string (int method);

int *system_get_endog_vars (const equation_system *sys);
int *system_get_instr_vars (const equation_system *sys);

void system_attach_uhat (equation_system *sys, gretl_matrix *uhat);

void system_attach_sigma (equation_system *sys, gretl_matrix *sigma);

void system_attach_coeffs (equation_system *sys, gretl_matrix *b);
void system_attach_vcv (equation_system *sys, gretl_matrix *vcv);

MODEL *system_get_model (const equation_system *sys, int i);

int system_get_overid_df (const equation_system *sys);

int system_vcv_geomean (const equation_system *sys);
double system_vcv_denom (const equation_system *sys, 
			 int i, int j);

int rhs_var_in_identity (const equation_system *sys, int lhsvar,
			 int rhsvar);

void 
print_equation_system_info (const equation_system *sys, 
			    const DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn);

void 
system_set_restriction_matrices (equation_system *sys,
				 gretl_matrix *R, gretl_matrix *q);

int 
system_normality_test (const equation_system *sys, PRN *prn);

int system_add_resids_to_dataset (equation_system *sys, int eqnum,
				  double ***pZ, DATAINFO *pdinfo);

double *
equation_system_get_series (const equation_system *sys, 
			    const DATAINFO *pdinfo,
			    int idx, const char *key, int *err);

gretl_matrix *
equation_system_get_matrix (const equation_system *sys, int idx, 
			    int *err);

int highest_numbered_var_in_system (const equation_system *sys, 
				    const DATAINFO *pdinfo);

int equation_system_serialize (equation_system *sys, 
			       SavedObjectFlags flags,
			       FILE *fp);

void system_set_save_flag (equation_system *sys);

void system_unset_save_flag (equation_system *sys);

int system_save_flag_is_set (equation_system *sys);

const gretl_matrix *
system_get_forecast_matrix (equation_system *sys, int t0, int t1, int t2,
			    const double **Z, DATAINFO *pdinfo, 
			    gretlopt opt, int *err);

#ifndef GRETLCLI

equation_system *
equation_system_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err);

#endif

int 
system_save_and_print_results (equation_system *sys,
			       double ***pZ, DATAINFO *pdinfo,
			       gretlopt opt, PRN *prn);

#endif /* EQUATION_SYSTEM_H */
