/* gretl - The Gnu Regression, Econometrics and Time-series Library
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

#include <stdio.h>

enum model_selection_criteria {
    C_SGMASQ = 0,
    C_AIC,
    C_FPE,
    C_HQ,
    C_BIC,
    C_SHIBATA,
    C_GCV,
    C_RICE,
    C_MAX
};

/**
 * dataset_is_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data (1) or not (0).
 */
#define dataset_is_time_series(p) (p->time_series == TIME_SERIES)

/**
 * dataset_is_daily:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains daily time series
 * data (1) or not (0).
 */
#define dataset_is_daily(p) (p->time_series == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7))

/**
 * dataset_is_weekly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains weekly time series
 * data (1) or not (0).
 */
#define dataset_is_weekly(p) (p->time_series == TIME_SERIES \
                              && p->pd == 52)

/**
 * dataset_is_hourly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains hourly time series
 * data (1) or not (0).
 */
#define dataset_is_hourly(p) (p->time_series == TIME_SERIES \
                              && p->pd == 24)

/**
 * dated_daily_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p->time_series == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                             && p->sd0 > 10000.0)

/**
 * dated_seven_day_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily 
 * (seven-day) time series data (1) or not (0).
 */
#define dated_seven_day_data(p) (p->time_series == TIME_SERIES \
                                 && p->pd == 7 && \
                                 p->sd0 > 10000.0)

/**
 * dataset_is_panel:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) (p->time_series == STACKED_TIME_SERIES \
                            || p->time_series == STACKED_CROSS_SECTION)


#include <float.h>
#define floateq(x, y) (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y) ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y) ((y) - (x) > DBL_EPSILON)
#define na(x) (fabs((x) + 999.0) < DBL_EPSILON)
#define NADBL -999.0

/* functions follow */
 
double date (int nt, int pd, const double sd0);

int isdummy (const double *x, int t1, int t2);

void printlist (const int *list, const char *msg);

int print_list_to_buffer (const int *list, char *buf, size_t len);

void list_exclude (int n, int *list);

char *format_obs (char *obs, int maj, int min, int pd);

int set_obs (const char *line, DATAINFO *pdinfo, gretlopt opt);

void set_miss (LIST list, const char *param, double **Z,
	       DATAINFO *pdinfo, PRN *prn);

int *copylist (const int *src);

int grow_nobs (int newobs, double ***pZ, DATAINFO *pdinfo);

int dataset_add_vars (int newvars, double ***pZ, DATAINFO *pdinfo);

int dataset_add_allocated_var (double *x, double ***pZ, DATAINFO *pdinfo);

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo);

int varnum_from_string (const char *str, DATAINFO *pdinfo);

int dataset_drop_listed_vars (const int *list, double ***pZ, 
			      DATAINFO *pdinfo, int *renumber);

int dataset_drop_vars (int delvars, double ***pZ, DATAINFO *pdinfo);

int rename_var_by_id (const char *str, const char *vname, 
		      DATAINFO *pdinfo);

int hidden_var (int i, const DATAINFO *pdinfo);

FITRESID *get_fit_resid (const MODEL *pmod, double ***pZ, 
			 DATAINFO *pdinfo);

FITRESID *get_fcast_with_errs (const char *str, const MODEL *pmod, 
			       double ***pZ, DATAINFO *pdinfo, PRN *prn);

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn,
		     PATHS *ppaths, int plot);

int is_model_cmd (const char *line);

int is_model_ref_cmd (int ci);

int is_quiet_model_test (int ci, gretlopt opt);

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo);

double *copyvec (const double *src, int n);

int ijton (int i, int j, int nrows);

int ztox (int i, double *px, 
	  double **Z, const DATAINFO *pdinfo);

int get_panel_structure (DATAINFO *pdinfo, int *nunits, int *T);

int set_panel_structure (gretlopt opt, DATAINFO *pdinfo, PRN *prn); 

int balanced_panel (const DATAINFO *pdinfo);

double get_xvalue (int i, const double **Z, const DATAINFO *pdinfo);

void free_gretl_mp_results (mp_results *mpvals);

mp_results *gretl_mp_results_new (int totvar);

int allocate_mp_varnames (mp_results *mpvals);

FITRESID *fit_resid_new (int n, int errs);

void free_fit_resid (FITRESID *fr);

CONFINT *get_model_confints (const MODEL *pmod);

void free_confint (CONFINT *cf);

#ifndef WIN32
int gretl_spawn (const char *cmdline);
int gretl_spawn_quiet (const char *cmdline);
#endif

#endif /* GRETL_UTILS_H */
