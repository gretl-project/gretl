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

enum Detflags {
    DET_CONST = 1 << 0,
    DET_TREND = 1 << 1,
    DET_SEAS  = 1 << 2
};

enum VAR_robust {
    VAR_HC = 1,
    VAR_HAC
};

typedef struct JohansenInfo_ JohansenInfo;

struct GRETL_VAR_ {
    int ci;              /* command index (VAR or VECM) */
    int refcount;        /* for saving/deleting */
    int err;             /* error code */
    int neqns;           /* number of equations in system */
    int order;           /* maximum lag order */
    int t1;              /* starting observation */
    int t2;              /* ending observation */
    int T;               /* number of observations */
    int df;              /* T - average coeffs per equation */
    int ifc;             /* equations include a constant (1) or not (0) */
    int ncoeff;          /* total coefficients per equation */
    int *lags;           /* list of specific lags */
    int *ylist;          /* list of stochastic vars */
    int *xlist;          /* list of exogenous variables */
    int *rlist;          /* restricted exogenous variables (VECM only) */
    int detflags;        /* record of automatic deterministic vars added */
    int robust;          /* computing robust std errors? */
    int xcols;           /* full column size of X matrix (VECM special) */
    gretl_matrix *Y;     /* matrix of dependent variables */
    gretl_matrix *X;     /* matrix of independent variables */
    gretl_matrix *B;     /* basic coefficient matrix */
    gretl_matrix *XTX;   /* X'X inverse */
    gretl_matrix *A;     /* augmented coefficient matrix (companion form) */
    gretl_matrix *L;     /* lambda: inverse roots of A(L) polynomial */
    gretl_matrix *E;     /* residuals matrix */
    gretl_matrix *C;     /* Cholesky-decomposed covariance matrix */
    gretl_matrix *S;     /* cross-equation variance matrix */
    gretl_matrix *F;     /* optional forecast matrix */
    gretl_matrix *V;     /* full parameter covariance matrix */
    gretl_matrix *ord;   /* optional Cholesky-ordering vector */
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
    int LBs;             /* order for Portmanteau test */
    JohansenInfo *jinfo; /* extra information for VECMs */
    char *name;          /* for use in session management */
};

#define var_n_lags(v) ((v->lags != NULL)? v->lags[0] : v->order)
#define var_max_lag(v) ((v->lags != NULL)? v->lags[v->lags[0]] : v->order)
#define levels_order(v) (v->order+(v->ci==VECM))

int var_max_order (const int *list, const DATASET *dset);

GRETL_VAR *gretl_VAR (int order, int *laglist, int *list,
		      const DATASET *dset, gretlopt opt,
		      PRN *prn, int *err);

GRETL_VAR *gretl_VECM (int order, int rank, int *list,
		       const DATASET *dset, gretlopt opt,
		       PRN *prn, int *err);

int gui_VAR_lagsel (int order, int *list,
		    const DATASET *dset,
		    gretlopt opt, PRN *prn,
		    gretl_matrix **pm);

const gretl_matrix *
gretl_VAR_get_forecast_matrix (GRETL_VAR *var, int t1, int t2,
			       DATASET *dset, gretlopt opt,
			       int *err);

const gretl_matrix *
gretl_VAR_get_residual_matrix (const GRETL_VAR *var);

gretl_matrix *
gretl_VAR_get_fcast_decomp (const GRETL_VAR *var,
			    int targ, int periods,
			    int *errp);

gretl_matrix *
gretl_VAR_get_vma_matrix (const GRETL_VAR *var, const DATASET *dset,
			  int *err);

gretl_matrix *gretl_VAR_get_FEVD_matrix (const GRETL_VAR *var,
					 int targ, int shock,
					 int horizon,
					 const DATASET *dset,
					 int *err);

gretl_matrix *
gretl_VAR_get_full_FEVD_matrix (const GRETL_VAR *var,
				const DATASET *dset,
				int *err);

gretl_matrix *gretl_FEVD_from_bundle (gretl_bundle *b,
				      int targ, int shock,
				      const DATASET *dset,
				      int *err);

gretl_matrix *gretl_IRF_from_bundle (gretl_bundle *b,
				     int targ, int shock,
				     double alpha,
				     const DATASET *dset,
				     int *err);

gretl_matrix *VECM_get_EC_matrix (const GRETL_VAR *v,
				  const DATASET *dset,
				  int *err);

int
gretl_VAR_do_error_decomp (const gretl_matrix *S, gretl_matrix *C,
			   const gretl_matrix *ord);

const gretl_matrix *gretl_VAR_get_roots (GRETL_VAR *var, int *err);

int gretl_VAR_autocorrelation_test (GRETL_VAR *var, int order,
				    DATASET *dset, gretlopt opt,
				    PRN *prn);

int gretl_VAR_arch_test (GRETL_VAR *var, int order,
			 DATASET *dset, gretlopt opt,
			 PRN *prn);

int gretl_VAR_normality_test (const GRETL_VAR *var,
			      gretlopt opt, PRN *prn);

int set_VAR_model_stats (GRETL_VAR *var, int i);

const int *gretl_VAR_get_exo_list (const GRETL_VAR *var);

const int *gretl_VAR_get_endo_list (const GRETL_VAR *var);

GRETL_VAR *gretl_VAR_omit_test (GRETL_VAR *var, const int *omitlist,
				DATASET *dset, gretlopt opt,
				PRN *prn, int *err);

int gretl_VAR_wald_omit_test (GRETL_VAR *var, const int *omitlist,
			      DATASET *dset, gretlopt opt,
			      PRN *prn);

double *gretl_VECM_get_EC (GRETL_VAR *vecm, int j, const DATASET *dset,
			   int *err);

void gretl_VAR_free (GRETL_VAR *var);

int default_VAR_horizon (const DATASET *dset);

gretl_matrix *
gretl_VAR_get_impulse_response (GRETL_VAR *var,
				int targ, int shock,
				int periods, double alpha,
				const DATASET *dset,
				int *err);

void gretl_VAR_set_name (GRETL_VAR *var, const char *name);

const char *gretl_VAR_get_name (const GRETL_VAR *var);

int gretl_VAR_get_variable_number (const GRETL_VAR *var, int k);

int gretl_VAR_get_n_equations (const GRETL_VAR *var);

int gretl_VAR_get_t1 (const GRETL_VAR *var);

int gretl_VAR_get_t2 (const GRETL_VAR *var);

int gretl_var_get_sample (const GRETL_VAR *var, int *t1, int *t2);

const MODEL *gretl_VAR_get_model (const GRETL_VAR *var, int i);

double *gretl_VAR_get_resid_series (GRETL_VAR *var, int eqnum,
				    int *err);

int gretl_VAR_set_ordering (GRETL_VAR *var, gretl_matrix *ord);

int gretl_VAR_do_irf (GRETL_VAR *var, const char *line,
		      const DATASET *dset);

int gretl_VAR_get_highest_variable (const GRETL_VAR *var);

GRETL_VAR *johansen_test (int order, const int *list,
			  const DATASET *dset, gretlopt opt,
			  PRN *prn);

int johansen_test_simple (int order, const int *list,
			  const DATASET *dset,
			  gretlopt opt, PRN *prn);

int gretl_VECM_n_beta (const GRETL_VAR *vecm);

int gretl_VECM_n_alpha (const GRETL_VAR *vecm);

int gretl_VECM_test (GRETL_VAR *vecm,
		     gretl_restriction *rset,
		     const DATASET *dset,
		     gretlopt opt,
		     PRN *prn);

int gretl_VAR_test (GRETL_VAR *var,
		    gretl_restriction *rset,
		    const DATASET *dset,
		    gretlopt opt,
		    PRN *prn);

GRETL_VAR *
real_gretl_restricted_vecm (GRETL_VAR *orig,
			    gretl_restriction *rset,
			    const DATASET *dset,
			    PRN *prn, int *err);

int gretl_VECM_rank (const GRETL_VAR *vecm);

int beta_restricted_VECM (const GRETL_VAR *vecm);

int alpha_restricted_VECM (const GRETL_VAR *vecm);

int restricted_VECM (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_R_matrix (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_q_matrix (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_Ra_matrix (const GRETL_VAR *vecm);

const gretl_matrix *gretl_VECM_qa_matrix (const GRETL_VAR *vecm);

double *gretl_VAR_get_series (const GRETL_VAR *var, const DATASET *dset,
			      int idx, const char *key, int *err);

gretl_matrix *gretl_VAR_get_matrix (const GRETL_VAR *var, int idx,
				    int *err);

void gretl_VAR_param_names (GRETL_VAR *v, char **params,
			    const DATASET *dset);

int gretl_VAR_serialize (const GRETL_VAR *var, SavedObjectFlags flags,
			 PRN *prn);

int gretl_VAR_bundlize (const GRETL_VAR *var, DATASET *dset,
			gretl_bundle *b);

int transcribe_VAR_models (GRETL_VAR *var,
			   const DATASET *dset,
			   const gretl_matrix *XTX);

#ifdef FULL_XML_HEADERS

GRETL_VAR *gretl_VAR_from_XML (xmlNodePtr node, xmlDocPtr doc,
			       const DATASET *dset,
			       int *err);

#endif

gretl_matrix *vma_rep (gretl_matrix *A, gretl_matrix *C,
		       int horizon, int *err);

#endif /* VAR_H_ */
