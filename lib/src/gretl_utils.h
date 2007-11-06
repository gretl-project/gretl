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

#ifndef GRETL_UTILS_H
#define GRETL_UTILS_H

#include "libgretl.h"

#include <float.h>
#include <limits.h>

#define floateq(x, y)  (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y)  ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y)  ((y) - (x) > DBL_EPSILON)

#define ok_int(x) (x <= (double) INT_MAX && x >= (double) INT_MIN)

void libgretl_init (void);

void libgretl_session_cleanup (void);

void libgretl_cleanup (void);

double date (int nt, int pd, const double sd0);

/* checks on variables */

int gretl_isdummy (int t1, int t2, const double *x);

int gretl_iszero (int t1, int t2, const double *x);

int gretl_isconst (int t1, int t2, const double *x);

int gretl_isunits (int t1, int t2, const double *x);

int gretl_isint (int t1, int t2, const double *x);

int gretl_isdiscrete (int t1, int t2, const double *x);

int gretl_is_oprobit_ok (int t1, int t2, const double *x);

int true_const (int v, const double **Z, const DATAINFO *pdinfo);

/* setting observations */

char *format_obs (char *obs, int maj, int min, int pd);

int set_obs (const char *line, double **Z, DATAINFO *pdinfo, 
	     gretlopt opt);

/* sorting and comparison */

int gretl_compare_doubles (const void *a, const void *b);

int gretl_inverse_compare_doubles (const void *a, const void *b);

int count_distinct_values (const double *x, int n);

int count_distinct_int_values (const int *x, int n);

int rearrange_id_array (double *x, int m, int n);

int gretl_compare_ints (const void *a, const void *b);

/* miscellaneous */

void printlist (const int *list, const char *msg);

int gretl_int_from_string (const char *s, const double **Z, 
			   const DATAINFO *pdinfo, int *err);

int positive_int_from_string (const char *s);

int varnum_from_string (const char *str, DATAINFO *pdinfo);

int rename_var_by_id (const char *idstr, const char *vname, 
		      DATAINFO *pdinfo);

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo);

double *copyvec (const double *src, int n);

void doubles_array_free (double **X, int m);

double **doubles_array_new (int m, int n);

double **data_array_from_model (const MODEL *pmod, double **Z, 
				int missv);

int ijton (int i, int j, int nrows);

int ztox (int i, double *px, const double **Z, const DATAINFO *pdinfo);

double get_xvalue (int i, const double **Z, const DATAINFO *pdinfo);

int gretl_copy_file (const char *src, const char *dest);

#ifndef WIN32
int gretl_spawn (char *cmdline);
#endif

/* model selection criteria */

int gretl_calculate_criteria (double ess, int n, int k,
			      double *ll, double *aic, double *bic, 
			      double *hqc);

int gretl_print_criteria (double ess, int nobs, int ncoeff, PRN *prn);

int ls_criteria (MODEL *pmod);

/* hypothesis tests mechanism */

void record_test_result (double teststat, double pval, char *blurb);

void record_LR_test_result (double teststat, double pval, double lnl,
			    char *blurb);

double get_last_test_statistic (char *blurb);

double get_last_pvalue (char *blurb);

double get_last_lnl (char *blurb);

/* timer */

double gretl_stopwatch (void);

#endif /* GRETL_UTILS_H */
