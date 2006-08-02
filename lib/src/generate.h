/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef GENERATE_H
#define GENERATE_H

enum {
    SORT_ASCENDING,
    SORT_DESCENDING
};

enum {
    HIGHNUM = 5000,
    TNUM,
    INDEXNUM
};

typedef enum {
    R_NOBS = 1,  /* number of observations in current sample range */
    R_NVARS,     /* number of variables in dataset (including the constant) */
    R_PD,        /* periodicity of dataset */
    R_DSET_MAX,  /* separator */
    R_TEST_STAT, /* test statistic from last explicit test performed */
    R_TEST_PVAL, /* p-value from last explicit test performed */
    R_MAX
} RetrievalIndex;

typedef struct _GENERATOR GENERATOR;

int generate (const char *line, double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn); 

GENERATOR *
genr_compile (const char *line, double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn);

int execute_genr (GENERATOR *genr, int oldv);

void destroy_genr (GENERATOR *genr);

int genr_get_varnum (const GENERATOR *genr);

int genr_get_err (const GENERATOR *genr);

int varindex (const DATAINFO *pdinfo, const char *varname);

int genr_fit_resid (const MODEL *pmod, 
		    double ***pZ, DATAINFO *pdinfo,
		    int code, int undo);

int get_generated_value (const char *argv, double *val,
			 double ***pZ, DATAINFO *pdinfo,
			 int t);

int gretl_reserved_word (const char *str);

int genr_function_from_string (const char *s);

int get_t_from_obs_string (char *s, const double **Z, 
			   const DATAINFO *pdinfo);

const char *get_retriever_word (int idx);

/* genrfuncs.c, public functions */

int hp_filter (const double *x, double *hp, const DATAINFO *pdinfo, 
	       gretlopt opt);

int bkbp_filter (const double *y, double *bk, const DATAINFO *pdinfo);

int dummy (double ***pZ, DATAINFO *pdinfo, int center);

int panel_unit_first_obs (int t, const DATAINFO *pdinfo);

int panel_dummies (double ***pZ, DATAINFO *pdinfo, gretlopt opt);

int genrtime (double ***pZ, DATAINFO *pdinfo, int tm);

const double *gretl_plotx (const DATAINFO *pdinfo);

const char *get_model_stat_word (int idx);

#endif /* GENERATE_H */

