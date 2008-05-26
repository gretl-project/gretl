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

#ifndef VAR_H_
#define VAR_H_

#include "gretl_matrix.h"
#include "gretl_restrict.h"

typedef struct JohansenInfo_ JohansenInfo;

struct GRETL_VAR_ {
    int ci;              /* command index (VAR or VECM) */
    int refcount;        /* for saving/deleting */
    int err;             /* error code */
    int neqns;           /* number of equations in system */
    int order;           /* lag order */
    int t1;              /* starting observation */
    int t2;              /* ending observation */
    int T;               /* number of observations */
    int df;              /* T - average coeffs per equation */
    int ifc;             /* equations include a constant (1) or not (0) */
    int ncoeff;          /* total coefficients per equation */
    int *ylist;          /* list of stochastic vars */
    int *xlist;          /* list of exogenous variables */
    int *rlist;          /* restricted exogenous variables (VECM only) */
    int detflags;        /* record of automatic deterministic vars added */
    int robust;          /* computing robust std errors? */
    int qr;              /* using QR decomposition? */
    gretl_matrix *Y;     /* matrix of dependent variables */
    gretl_matrix *X;     /* matrix of independent variables */
    gretl_matrix *B;     /* basic coefficient matrix */
    gretl_matrix *XTX;   /* X'X */
    gretl_matrix *A;     /* augmented coefficient matrix (companion form) */
    gretl_matrix *L;     /* lambda: inverse roots of A(L) polynomial */
    gretl_matrix *E;     /* residuals matrix */
    gretl_matrix *C;     /* augmented Cholesky-decomposed error matrix */
    gretl_matrix *S;     /* cross-equation variance matrix */
    gretl_matrix *vcv;   /* parameter covariance matrix */
    gretl_matrix *F;     /* optional forecast matrix */
    MODEL **models;      /* pointers to individual equation estimates */
    double *Fvals;       /* hold results of F-tests */
    double *Ivals;       /* hold results of info criteria comparisons */
    double ldet;         /* log-determinant of S */
    double ll;           /* log-likelihood */
    double AIC;          /* Akaike criterion */
    double BIC;          /* Bayesian criterion */
    double HQC;          /* Hannan-Quinn criterion */
    double LR;           /* for likelihood-ratio testing */
    double LB;           /* Ljung-Box (Portmanteau) test statistic */
    int LBs;             /* order for for Portmanteau test */
    JohansenInfo *jinfo; /* extra information for VECMs */
    char *name;          /* for use in session management */
};

int var_max_order (const int *list, const DATAINFO *pdinfo);

GRETL_VAR *gretl_VAR (int order, int *list, 
		      const double **Z, const DATAINFO *pdinfo,
		      gretlopt opt, PRN *prn, int *errp);

GRETL_VAR *gretl_VECM (int order, int rank, int *list, 
		       const double **Z, const DATAINFO *pdinfo,
		       gretlopt opt, PRN *prn, int *err);

const gretl_matrix *
gretl_VAR_get_forecast_matrix (GRETL_VAR *var, int t1, int t2, 
			       const double **Z, DATAINFO *pdinfo,
			       gretlopt opt, int *err);

const gretl_matrix *
gretl_VAR_get_residual_matrix (const GRETL_VAR *var);

gretl_matrix *
gretl_VAR_get_fcast_decomp (GRETL_VAR *var, int targ, int periods,
			    int *errp);

int 
gretl_VAR_do_error_decomp (const gretl_matrix *S, gretl_matrix *C);

const gretl_matrix *gretl_VAR_get_roots (GRETL_VAR *var, int *err);

int gretl_VAR_autocorrelation_test (GRETL_VAR *var, int order, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn);

int gretl_VAR_arch_test (GRETL_VAR *var, int order, 
			 DATAINFO *pdinfo, PRN *prn);

int gretl_VAR_normality_test (const GRETL_VAR *var, PRN *prn);

int set_VAR_model_stats (GRETL_VAR *var, int i);

const int *gretl_VAR_get_exo_list (const GRETL_VAR *var);

GRETL_VAR *gretl_VAR_omit_test (const int *omitvars, const GRETL_VAR *orig, 
				const double **Z, DATAINFO *pdinfo, 
				PRN *prn, int *err);

double *gretl_VECM_get_EC (GRETL_VAR *vecm, int j, const double **Z, 
			   const DATAINFO *pdinfo, int *err);

void gretl_VAR_free (GRETL_VAR *var);

int default_VAR_horizon (const DATAINFO *pdinfo);

gretl_matrix *
gretl_VAR_get_impulse_response (GRETL_VAR *var, 
				int targ, int shock,
				int periods,
				const double **Z,
				const DATAINFO *pdinfo);

void gretl_VAR_set_name (GRETL_VAR *var, const char *name);

const char *gretl_VAR_get_name (const GRETL_VAR *var);

int gretl_VAR_get_variable_number (const GRETL_VAR *var, int k);

int gretl_VAR_get_n_equations (const GRETL_VAR *var);

int gretl_VAR_get_t1 (const GRETL_VAR *var);

int gretl_VAR_get_t2 (const GRETL_VAR *var);

const MODEL *gretl_VAR_get_model (const GRETL_VAR *var, int i);

int gretl_VAR_add_resids_to_dataset (GRETL_VAR *var, int eqnum,
				     double ***pZ, DATAINFO *pdinfo);

int gretl_VAR_do_irf (GRETL_VAR *var, const char *line,
		      const double **Z, const DATAINFO *pdinfo);

int gretl_VAR_get_highest_variable (const GRETL_VAR *var);

int gretl_VECM_n_beta (const GRETL_VAR *vecm);

int gretl_VECM_n_alpha (const GRETL_VAR *vecm);

int gretl_VECM_test (GRETL_VAR *vecm, 
		     gretl_restriction *rset,
		     const DATAINFO *pdinfo, 
		     gretlopt opt,
		     PRN *prn);

GRETL_VAR *
real_gretl_restricted_vecm (GRETL_VAR *orig, 
			    gretl_restriction *rset,
			    const double **Z, const DATAINFO *pdinfo, 
			    PRN *prn, int *err);

int gretl_VECM_rank (const GRETL_VAR *vecm);

const int *gretl_VECM_list (const GRETL_VAR *vecm);

int beta_restricted_VECM (const GRETL_VAR *vecm);

int alpha_restricted_VECM (const GRETL_VAR *vecm);

int restricted_VECM (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_R_matrix (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_q_matrix (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_Ra_matrix (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_qa_matrix (const GRETL_VAR *vecm);

double *gretl_VAR_get_series (const GRETL_VAR *var, const DATAINFO *pdinfo, 
			      int idx, const char *key, int *err);

gretl_matrix *gretl_VAR_get_matrix (const GRETL_VAR *var, int idx, 
				    int *err);

void gretl_VAR_param_names (GRETL_VAR *v, char **params, 
			    const DATAINFO *pdinfo);

int gretl_VAR_serialize (const GRETL_VAR *var, SavedObjectFlags flags,
			 FILE *fp);

int transcribe_VAR_models (GRETL_VAR *var, 
			   const double **Z,
			   const DATAINFO *pdinfo,
			   const gretl_matrix *XTX);

#ifndef GRETLCLI

GRETL_VAR *gretl_VAR_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err);

#endif

#endif /* VAR_H_ */

