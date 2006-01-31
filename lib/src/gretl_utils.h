/*
 * Copyright (C) 1999-2004 Ramu Ramanathan and Allin Cottrell
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

#ifndef GRETL_UTILS_H
#define GRETL_UTILS_H

#include "libgretl.h"

typedef enum {
    C_AIC,
    C_BIC,
    C_HQC,
    C_MAX
} model_selection_criteria;

#include <float.h>

#define floateq(x, y)  (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y)  ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y)  ((y) - (x) > DBL_EPSILON)

/* functions follow */

void libgretl_init (void);

void libgretl_cleanup (void);

double date (int nt, int pd, const double sd0);

/* checks on variables */

int gretl_isdummy (int t1, int t2, const double *x);

int gretl_iszero (int t1, int t2, const double *x);

int gretl_isconst (int t1, int t2, const double *x);

int gretl_isunits (int t1, int t2, const double *x);

int true_const (int v, const double **Z, const DATAINFO *pdinfo);

/* setting observations */

char *format_obs (char *obs, int maj, int min, int pd);

int set_obs (const char *line, DATAINFO *pdinfo, gretlopt opt);

/* other */

void printlist (const int *list, const char *msg);

int gretl_int_from_string (const char *s, const double **Z, 
			   const DATAINFO *pdinfo, int *err);

int positive_int_from_string (const char *s);

int varnum_from_string (const char *str, DATAINFO *pdinfo);

int rename_var_by_id (const char *str, const char *vname, 
		      DATAINFO *pdinfo);

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo);

double *copyvec (const double *src, int n);

int ijton (int i, int j, int nrows);

int ztox (int i, double *px, const double **Z, const DATAINFO *pdinfo);

double get_xvalue (int i, const double **Z, const DATAINFO *pdinfo);

int gretl_compare_doubles (const void *a, const void *b);

int gretl_inverse_compare_doubles (const void *a, const void *b);

int gretl_copy_file (const char *src, const char *dest);

#ifndef WIN32
int gretl_spawn (char *cmdline);
#endif

/* model selection criteria */

int gretl_calculate_criteria (double ess, int nobs, int ncoeff,
			      double *ll, double *aic, double *bic, 
			      double *hqc);

int gretl_print_criteria (double ess, int nobs, int ncoeff, PRN *prn);

int ls_criteria (MODEL *pmod);

/* panel data utilities */

int get_panel_structure (const DATAINFO *pdinfo, int *nunits, int *T);

int set_panel_structure (gretlopt opt, DATAINFO *pdinfo, PRN *prn); 

int balanced_panel (const DATAINFO *pdinfo);

/* multiple-precision utilities */

void free_gretl_mp_results (mp_results *mpvals);

mp_results *gretl_mp_results_new (int totvar);

int allocate_mp_varnames (mp_results *mpvals);

/* hypothesis tests mechanism */

void record_test_result (double teststat, double pval, char *blurb);

double get_last_test_statistic (char *blurb);

double get_last_pvalue (char *blurb);

#endif /* GRETL_UTILS_H */
