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

#ifndef GENMAIN_H
#define GENMAIN_H

#include "gretl_matrix.h"

typedef enum {
    R_NOBS = 1,   /* number of observations in current sample range */
    R_NVARS,      /* number of variables in dataset (including the constant) */
    R_PD,         /* periodicity of dataset */
    R_T1,         /* start of current sample range */
    R_T2,         /* end of current sample range */
    R_DATATYPE,   /* dataset structure (x-section, time-series, panel) */
    R_WINDOWS,    /* running on MS Windows (1) or not (0) */
    R_VERSION,    /* gretl version number */
    R_ERRNO,      /* internal gretl error code */
    R_SEED,       /* RNG seed */
    R_HUGE,       /* conventional "huge" number, like eg 1.0e100 */
    R_DSET_MAX,   /* separator */
    R_TEST_LNL,   /* log-likelihood from last test (if applicable) */
    R_KLNL,       /* log-likelihood from Kalman filter (if applicable) */
    R_KS2,        /* variance estimate from Kalman filter (if applicable) */
    R_KSTEP,      /* current Kalman time-step (if applicable) */
    R_STOPWATCH,  /* stopwatch */ 
    R_SCALAR_MAX, /* separator: scalars vs series */
    R_INDEX,      /* consecutive observations index */
    R_PUNIT,      /* 1-based panel unit index */
    R_OBSMAJ,     /* major component of observation (e.g. year) */
    R_OBSMIN,     /* minor component of observation (e.g. quarter, month) */
    R_OBSMIC,     /* micro component of observation (e.g. day) */
    R_DATES,      /* ISO 8601 "basic" dates series */
    R_SERIES_MAX, /* separator: series vs matrices */
    R_TEST_STAT,  /* last test statistic(s) (scalar or matrix) */
    R_TEST_PVAL,  /* last test p-value(s) (scalar or matrix) */
    R_MAX
} RetrievalIndex;

/**
 * ModelDataIndex:
 * 
 * Symbolic names for various scalars, series and matrices
 * that can be retrieved from gretl models.
 */

typedef enum {
    M_ESS = R_MAX + 1, /* error sum of squares */
    M_T,          /* observations used */
    M_RSQ,        /* R-squared */
    M_SIGMA,      /* standard error of residuals */
    M_DF,         /* degrees of freedom */
    M_NCOEFF,     /* total number of estimated coefficients */
    M_LNL,        /* log-likelihood */
    M_GMMCRIT,    /* GMM criterion */
    M_AIC,        /* Akaike info criterion */
    M_BIC,        /* Bayesian info criterion */
    M_HQC,        /* Hannan-Quinn criterion */
    M_TRSQ,       /* T * R-squared, last model */
    M_DWPVAL,     /* Durbin-Watson p-value, last model */
    M_FSTT,       /* overall F-statistic, last model */
    M_CHISQ,      /* overall chi-square stat, last model */
    M_DIAGTEST,   /* system test for diagonal covariance matrix */
    M_DIAGPVAL,   /* p-value for the above */
    M_SCALAR_MAX, /* -- SEPARATOR, scalars/series -- */
    M_UHAT,       /* residuals */
    M_YHAT,       /* fitted values */
    M_LLT,        /* per-observation loglikelihood */
    M_AHAT,       /* per-unit intercepts in panel model */
    M_H,          /* GARCH predicted variances */
    M_SAMPLE,     /* observations used in estimation */
    M_UHAT2,      /* squared residuals */
    M_SERIES_MAX, /* -- SEPARATOR, series/matrices -- */
    M_COEFF,      /* parameter estimates */
    M_SE,         /* parameter standard errors */
    M_VCV,        /* parameter covariance matrix */
    M_RHO,        /* autoregressive coefficients */
    M_COMPAN,     /* VAR companion matrix */
    M_XTXINV,     /* VARs, VECMs: X'X^{-1} */
    M_VECG,       /* VECMs: the Gamma matrices */
    M_EVALS,      /* VECMs: eigenvalues */
    M_JALPHA,     /* Johansen's alpha */
    M_JBETA,      /* Johansen's beta */
    M_JVBETA,     /* Covariance matrix for Johansen's normalized beta */
    M_JS00,       /* VECM residual covariance matrix (1st differences) */
    M_JS11,       /* VECM residual covariance matrix (levels) */
    M_JS01,       /* VECM residual cross-product matrix */
    M_HAUSMAN,    /* Hausman test after tsls or fixed effects */
    M_SARGAN,     /* Sargan over-identification test after tsls */
    M_SYSGAM,     /* Parameter matrix Gamma (simultaneous systems) */
    M_SYSA,       /* Parameter matrix A (simultaneous systems) */
    M_SYSB,       /* Parameter matrix B (simultaneous systems) */
    M_FCAST,      /* last forecast generated via fcast command */
    M_FCERR,      /* standard errors associated with M_FCAST */
    M_COEFF_CI,   /* (asymmetric) confidence intervals for coeffs */
    M_KLLT,       /* Kalman log-likelihood, per time-step */
    M_KUHAT,      /* Kalman: current prediction error */
    M_EHAT,       /* ARMA: vector of estimated innovations */
    M_MATRIX_MAX, /* -- SEPARATOR, end of matrices -- */
    M_EC,         /* VECM error-correction terms */
    M_VMA,        /* VARs, VECMs: vector moving average representation */
    M_FEVD,       /* VAR variance decomposition */
    M_MNLPROBS,   /* case probabilities for multinomial logit */
    M_MBUILD_MAX, /* -- SEPARATOR, end of matrix-builders -- */
    M_XLIST,      /* list of regressors */
    M_YLIST,      /* list of endogenous variables */
    M_LIST_MAX,   /* -- SEPARATOR, end of lists -- */
    M_COMMAND,    /* model command word */
    M_DEPVAR,     /* name of dependent variable */
    M_MAX         /* sentinel */
} ModelDataIndex;

#define model_data_scalar(i) (i > R_MAX && i < M_SCALAR_MAX)
#define model_data_series(i) (i > M_SCALAR_MAX && i < M_SERIES_MAX)
#define model_data_matrix(i) (i > M_SERIES_MAX && i < M_MATRIX_MAX)
#define model_data_matrix_builder(i) (i > M_MATRIX_MAX && i < M_MBUILD_MAX)
#define model_data_list(i)   (i > M_MBUILD_MAX && i < M_LIST_MAX)
#define model_data_string(i) (i > M_LIST_MAX && i < M_MAX)

typedef struct parser_ GENERATOR;

int generate (const char *line, DATASET *dset, 
	      gretlopt opt, PRN *prn);

GENERATOR *genr_compile (const char *s, DATASET *dset, 
			 gretlopt opt, int *err);

int execute_genr (GENERATOR *genr, DATASET *dset, PRN *prn);

void destroy_genr (GENERATOR *genr);

int genr_get_output_type (const GENERATOR *genr);

int genr_get_output_varnum (const GENERATOR *genr);

double genr_get_output_scalar (const GENERATOR *genr);

int genr_get_last_output_type (void);

gretl_matrix *genr_get_output_matrix (const GENERATOR *genr);

int series_index (const DATASET *dset, const char *varname);

int series_greatest_index (const DATASET *dset, const char *varname);

int current_series_index (const DATASET *dset, const char *vname);

int extract_varname (char *targ, const char *src, int *len);

int genr_fit_resid (const MODEL *pmod, DATASET *dset,
		    ModelDataIndex idx);

double evaluate_if_cond (GENERATOR *genr, DATASET *dset, int *err);

double generate_scalar (const char *s, DATASET *dset, int *err);

double *generate_series (const char *s, DATASET *dset, PRN *prn,
			 int *err);

gretl_matrix *generate_matrix (const char *s, DATASET *dset, 
			       int *err);

char *generate_string (const char *s, DATASET *dset, int *err);

int *generate_list (const char *s, DATASET *dset, int *err);

int gretl_is_series (const char *name, const DATASET *dset);

int gretl_reserved_word (const char *str);

int genr_special_word (const char *s);

int genr_function_word (const char *s);

int genr_is_print (const GENERATOR *p);

int genr_is_autoregressive (const GENERATOR *p);

void genr_set_na_check (GENERATOR *genr);

void genr_unset_na_check (GENERATOR *genr);

int genr_get_series_max (GENERATOR *genr);

int function_from_string (const char *s);

int function_lookup (const char *s);

int const_lookup (const char *s);

double get_const_by_name (const char *name, int *err);

const char *gretl_function_complete (const char *s);

void gretl_function_hash_cleanup (void);

void set_mpi_rank_and_size (int rank, int size);

#endif /* GENMAIN_H */

