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
    R_DSET_MAX,   /* separator */
    R_TEST_STAT,  /* test statistic from last explicit test performed */
    R_TEST_PVAL,  /* p-value from last explicit test performed */
    R_TEST_LNL,   /* log-likelihood from last test (if applicable) */
    R_STOPWATCH,  /* stopwatch */ 
    R_NSCAN,      /* number of items scanned via sscanf */
    R_SCALAR_MAX, /* separator: scalars vs series */
    R_INDEX,      /* consecutive observations index */
    R_PUNIT,      /* 1-based panel unit index */
    R_MAX
} RetrievalIndex;

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
    M_SCALAR_MAX, /* -- separator, scalars/series -- */
    M_UHAT,       /* residuals */
    M_YHAT,       /* fitted values */
    M_AHAT,       /* per-unit intercepts in panel model */
    M_H,          /* GARCH predicted variances */
    M_SAMPLE,     /* observations used in estimation */
    M_SERIES_MAX, /* -- separator, series/matrices -- */
    M_COEFF,      /* parameter estimates */
    M_SE,         /* parameter standard errors */
    M_VCV,        /* parameter covariance matrix */
    M_RHO,        /* autoregressive coefficients */
    M_COMPAN,     /* VAR companion matrix */
    M_JALPHA,     /* Johansen's alpha */
    M_JBETA,      /* Johansen's beta */
    M_JVBETA,     /* Covariance matrix for Johansen's normalized beta */
    M_JS00,       /* VECM residual covariance matrix (1st differences) */
    M_JS11,       /* VECM residual covariance matrix (levels) */
    M_JS01,       /* VECM residual cross-product matrix */
    M_HAUSMAN,    /* Hausman test after tsls */
    M_SARGAN,     /* Sargan over-identification test after tsls */
    M_SYSGAM,     /* Parameter matrix Gamma (simultaneous systems) */
    M_SYSA,       /* Parameter matrix A (simultaneous systems) */
    M_SYSB,       /* Parameter matrix B (simultaneous systems) */
    M_FCAST,      /* last forecast generated via fcasterr */
    M_FCERR,      /* standard errors associated with M_FCAST */
    M_COEFF_CI,   /* (asymmetric) confidence intervals for coeffs */
    M_MATRIX_MAX, /* -- separator, matrices/lists -- */
    M_XLIST,      /* list of regressors */
    M_MAX         /* sentinel */
} ModelDataIndex;

#define model_data_scalar(i) (i > R_MAX && i < M_SCALAR_MAX)
#define model_data_series(i) (i > M_SCALAR_MAX && i < M_SERIES_MAX)
#define model_data_matrix(i) (i > M_SERIES_MAX && i < M_MATRIX_MAX)
#define model_data_list(i)   (i > M_MATRIX_MAX && i < M_MAX)

typedef struct parser_ GENERATOR;

int generate (const char *line, double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn);

GENERATOR *genr_compile (const char *s, double ***pZ, DATAINFO *pdinfo, 
			 gretlopt opt, int *err);

int execute_genr (GENERATOR *genr, double ***pZ, DATAINFO *pdinfo,
		  gretlopt opt, PRN *prn);

void destroy_genr (GENERATOR *genr);

void genr_set_loopline (GENERATOR *genr, int i);

int genr_get_loopline (GENERATOR *genr);

int genr_get_output_varnum (const GENERATOR *genr);

const char *genr_get_formula (const GENERATOR *genr);

gretl_matrix *genr_get_output_matrix (const GENERATOR *genr);

int varindex (const DATAINFO *pdinfo, const char *varname);

int extract_varname (char *targ, const char *src, int *len);

int genr_fit_resid (const MODEL *pmod, 
		    double ***pZ, DATAINFO *pdinfo,
		    int code, int undo);

double generate_scalar (const char *s, double ***pZ, 
			DATAINFO *pdinfo, int *err);

double *generate_series (const char *s, double ***pZ, 
			 DATAINFO *pdinfo, int *err);

gretl_matrix *generate_matrix (const char *s, double ***pZ, 
			       DATAINFO *pdinfo, int *err);

char *generate_string (const char *s, double ***pZ, 
		       DATAINFO *pdinfo, int *err);

int *generate_list (const char *s, double ***pZ, 
		    DATAINFO *pdinfo, int *err);

int print_object_var (const char *oname, const char *param,
		      double ***pZ, DATAINFO *pdinfo,
		      PRN *prn);

int gretl_reserved_word (const char *str);

int genr_special_word (const char *s);

int genr_function_word (const char *s);

int genr_is_print (const GENERATOR *p);

void genr_set_na_check (GENERATOR *genr);

void genr_unset_na_check (GENERATOR *genr);

int function_from_string (const char *s);

int is_gretl_function_call (const char *s);

int function_lookup (const char *s);

int const_lookup (const char *s);

void gretl_function_hash_cleanup (void);

#endif /* GENMAIN_H */

