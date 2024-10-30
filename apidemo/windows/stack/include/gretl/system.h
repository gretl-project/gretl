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
    /* note: allow for obsoleted flags in saved sessions --
       treat this enumeration as append-only 
    */
    SYSTEM_DFCORR      = 1 << 2,
    SYSTEM_VCV_GEOMEAN = 1 << 3,
    SYSTEM_RESTRICT    = 1 << 4,
    SYSTEM_ITERATE     = 1 << 5,
    SYSTEM_SAVEIT      = 1 << 6,
    SYSTEM_LIML1       = 1 << 7,
    SYSTEM_QUIET       = 1 << 8,
    SYSTEM_ROBUST      = 1 << 9
};

typedef struct id_atom_ id_atom;
typedef struct identity_ identity;
typedef struct predet_ predet;
typedef struct liml_data_ liml_data;

struct equation_system_ {
    char *name;                 /* user-specified name for system, or NULL */
    int refcount;               /* for saving/deleting */
    int fd;                     /* function execution depth */
    int t1;                     /* starting observation number */
    int t2;                     /* ending observation number */
    int smpl_t1;                /* first obs in sample range */
    int smpl_t2;                /* last obs in sample range */
    int T;                      /* number of observations per equation */
    int df;                     /* T - average coeffs per equation */
    int method;                 /* estimation method */
    int neqns;                  /* number of stochastic equations */
    int nidents;                /* number of identities */
    int order;                  /* max lag of endogenous variable */
    int iters;                  /* number of iterations taken */
    int flags;                  /* to record options */
    double ll;                  /* log-likelihood (restricted) */
    double llu;                 /* unrestricted log-likelihood */
    double X2;                  /* chi-square test value */
    double ess;                 /* total error sum of squares */
    double diag_test;           /* test stat for diagonal covariance matrix */
    double bdiff;               /* summary stat for change in coefficients */
    double ldet;                /* log-determinant of covariance matrix */
    int **lists;                /* regression lists for stochastic equations */
    int **tslists;              /* back-up of TSLS-style lists */
    int *ylist;                 /* list of endogenous variables */
    int *ilist;                 /* list of instruments */
    int *xlist;                 /* list of truly exogenous variables */
    int *plist;                 /* list of predetermined variables */
    int *biglist;               /* list of all variables, for data checking */
    predet *pre_vars;           /* array of info on predetermined regressors */
    identity **idents;          /* set of identities */
    gretl_matrix *b;            /* coefficient estimates */
    gretl_matrix *vcv;          /* covariance matrix of coefficients */
    gretl_matrix *S;            /* cross-equation covariance matrix */
    gretl_matrix *R;            /* LHS of any linear restrictions */
    gretl_matrix *q;            /* RHS of any linear restrictions */  
    gretl_matrix *E;            /* residuals, all equations */
    gretl_matrix *yhat;         /* fitted values, all equations */
    gretl_matrix *Gamma;        /* structural form Gamma matrix (endog + identities)*/
    gretl_matrix *B;            /* structural form B matrix (exogenous) */
    gretl_matrix *A;            /* structural form A matrix (lagged endogenous) */
    gretl_matrix *F;            /* forecast matrix */
    gretl_matrix *Sr;           /* reduced-form error covariance matrix */
    MODEL **models;             /* set of pointers to per-equation models */
    liml_data *ldata;           /* extra info from LIML estimation */
};

equation_system *equation_system_start (const char *param, 
					const char *name,
					gretlopt opt,
					int *err);

char *get_system_name_from_line (const char *s);

equation_system *get_anonymous_equation_system (void);

int equation_system_append (equation_system *sys, const int *list);

int equation_system_append_multi (equation_system *sys, 
				  const char *parm1,
				  const char *parm2,
				  const DATASET *dset);

int system_parse_line (equation_system *sys,
		       const char *line,
		       DATASET *dset);

int equation_system_finalize (equation_system *sys, 
			      DATASET *dset,
			      gretlopt opt, PRN *prn);

int 
equation_system_estimate (equation_system *sys, 
			  DATASET *dset, 
			  gretlopt opt, PRN *prn);

int estimate_named_system (const char *sysname, 
			   const char *param,
			   DATASET *dset, 
			   gretlopt opt, PRN *prn);

void equation_system_destroy (equation_system *sys);

void delete_anonymous_equation_system (int level);

int system_want_df_corr (const equation_system *sys);

int system_n_restrictions (const equation_system *sys);

int system_max_indep_vars (const equation_system *sys);

int system_n_indep_vars (const equation_system *sys);

int *system_get_list (const equation_system *sys, int i);

int system_get_list_length (const equation_system *sys, int i);

int *compose_ivreg_list (const equation_system *sys, int i);

int system_get_depvar (const equation_system *sys, int i);

const char *system_short_string (const MODEL *pmod);

void equation_system_set_name (equation_system *sys, const char *name);

int system_method_from_string (const char *s);
const char *system_method_full_string (int method);
const char *system_method_short_string (int method);

int *system_get_endog_vars (const equation_system *sys);
int *system_get_instr_vars (const equation_system *sys);

void system_attach_uhat (equation_system *sys, gretl_matrix *E);

void system_attach_sigma (equation_system *sys, gretl_matrix *S);

void system_attach_coeffs (equation_system *sys, gretl_matrix *b);
void system_attach_vcv (equation_system *sys, gretl_matrix *vcv);

MODEL *system_get_model (const equation_system *sys, int i);

int system_get_overid_df (const equation_system *sys);

int system_vcv_geomean (const equation_system *sys);
double system_vcv_denom (const equation_system *sys, 
			 int i, int j);

int rhs_var_in_identity (const equation_system *sys,
			 int lhsvar, int rhsvar);

void 
print_equation_system_info (const equation_system *sys, 
			    const DATASET *dset, 
			    gretlopt opt, PRN *prn);

void
system_set_restriction_matrices (equation_system *sys,
				 gretl_matrix *R,
				 gretl_matrix *q);

int 
system_normality_test (const equation_system *sys,
		       gretlopt opt, PRN *prn);

int multi_eqn_wald_test (const gretl_matrix *b,
			 const gretl_matrix *V,
			 const gretl_matrix *R,
			 const gretl_matrix *q,
			 int dfu, gretlopt opt,
			 PRN *prn);

int system_wald_test (const equation_system *sys, 
		      const gretl_matrix *R,
		      const gretl_matrix *q,
		      gretlopt opt,
		      PRN *prn);

int system_diag_test (const equation_system *sys,
		      double *test, double *pval);

double *system_get_resid_series (equation_system *sys, int eqnum,
				 DATASET *dset, int *err);

double *
equation_system_get_series (const equation_system *sys, 
			    const DATASET *dset,
			    int idx, const char *key,
			    int *err);

gretl_matrix *
equation_system_get_matrix (const equation_system *sys, int idx, 
			    int *err);

int highest_numbered_var_in_system (const equation_system *sys, 
				    const DATASET *dset);

int equation_system_serialize (equation_system *sys, 
			       SavedObjectFlags flags,
			       PRN *prn);

int equation_system_bundlize (equation_system *sys,
			      gretl_bundle *b);

int gretl_system_print (equation_system *sys, 
			const DATASET *dset, 
			gretlopt opt, PRN *prn);

int system_print_sigma (const equation_system *sys, PRN *prn);

const gretl_matrix *
system_get_forecast_matrix (equation_system *sys, int t1, int t2,
			    DATASET *dset, gretlopt opt, 
			    int *err);

gretl_matrix *sys_get_fitted_values (equation_system *sys,
				     int v, int t1, int t2,
				     const DATASET *dset,
				     int *err);

int system_adjust_t1t2 (equation_system *sys, 
			const DATASET *dset);

int system_supports_method (equation_system *sys, int method);

#ifdef FULL_XML_HEADERS

equation_system *equation_system_from_XML (xmlNodePtr node, 
					   xmlDocPtr doc, 
					   int *err);

#endif

int 
system_save_and_print_results (equation_system *sys,
			       DATASET *dset, gretlopt opt, 
			       PRN *prn);

int system_autocorrelation_test (equation_system *sys, int order, 
				 gretlopt opt, PRN *prn);

int system_arch_test (equation_system *sys, int order, 
		      gretlopt opt, PRN *prn);

MODEL single_equation_liml (const int *list, DATASET *dset, 
			    gretlopt opt);

int gretl_system_get_sample (const equation_system *sys,
			     int *t1, int *t2);

#endif /* EQUATION_SYSTEM_H */
