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

enum model_selection_criteria {
    C_AIC,
    C_BIC,
    C_MAX
};

/**
 * dataset_is_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data (1) or not (0).
 */
#define dataset_is_time_series(p) ((p)->structure == TIME_SERIES || \
				   (p)->structure == SPECIAL_TIME_SERIES)

/**
 * custom_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data with custom (non-standard) frequency (1) or not (0).
 */
#define custom_time_series(p) ((p)->structure == SPECIAL_TIME_SERIES)

/**
 * dataset_is_daily:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains daily time series
 * data (1) or not (0).
 */
#define dataset_is_daily(p) (p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7))

/**
 * dataset_is_weekly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains weekly time series
 * data (1) or not (0).
 */
#define dataset_is_weekly(p) (p->structure == TIME_SERIES \
                              && p->pd == 52)

/**
 * dataset_is_hourly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains hourly time series
 * data (1) or not (0).
 */
#define dataset_is_hourly(p) (p->structure == TIME_SERIES \
                              && p->pd == 24)

/**
 * dataset_is_decennial:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains decennial time series
 * data (1) or not (0).
 */
#define dataset_is_decennial(p) (p->structure == TIME_SERIES \
                                 && p->pd == 10)

/**
 * dated_daily_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                             && p->sd0 > 10000.0)

/**
 * dated_seven_day_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily 
 * (seven-day) time series data (1) or not (0).
 */
#define dated_seven_day_data(p) (p->structure == TIME_SERIES \
                                 && p->pd == 7 && \
                                 p->sd0 > 10000.0)

/**
 * dated_weekly_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated weekly 
 * time series data (1) or not (0).
 */
#define dated_weekly_data(p) (p->structure == TIME_SERIES \
                              && p->pd == 52 && \
                              p->sd0 > 10000.0)

/**
 * calendar_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set uses calendar
 * dates for observation strings (1) or not (0).
 */
#define calendar_data(p) (p->structure == TIME_SERIES && \
                          (p->pd == 5 || p->pd == 6 || p->pd == 7 \
                           || p->pd == 52) && p->sd0 > 10000.0) 
                          
/**
 * dataset_is_panel:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) ((p)->structure == STACKED_TIME_SERIES || \
                             (p)->structure == STACKED_CROSS_SECTION)

#include <float.h>

#define floateq(x, y) (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y) ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y) ((y) - (x) > DBL_EPSILON)

/* functions follow */

void libgretl_init (CMD *cmd);

void libgretl_cleanup (CMD *cmd);
 
double date (int nt, int pd, const double sd0);

int isdummy (const double *x, int t1, int t2);

void printlist (const int *list, const char *msg);

int print_list_to_buffer (const int *list, char *buf, size_t len);

void list_exclude (int n, int *list);

char *format_obs (char *obs, int maj, int min, int pd);

int set_obs (const char *line, DATAINFO *pdinfo, gretlopt opt);

int *copylist (const int *src);

int grow_nobs (int newobs, double ***pZ, DATAINFO *pdinfo);

int dataset_add_vars (int newvars, double ***pZ, DATAINFO *pdinfo);

int dataset_add_allocated_var (double *x, double ***pZ, DATAINFO *pdinfo);

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo);

int positive_int_from_string (const char *s);

int varnum_from_string (const char *str, DATAINFO *pdinfo);

int dataset_drop_listed_vars (const int *list, double ***pZ, 
			      DATAINFO *pdinfo, int *renumber);

int dataset_destroy_hidden_vars (double ***pZ, DATAINFO *pdinfo);

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
		     int plot);

int is_model_cmd (const char *line);

int is_model_ref_cmd (int ci);

int is_quiet_model_test (int ci, gretlopt opt);

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo);

double *copyvec (const double *src, int n);

int ijton (int i, int j, int nrows);

int ztox (int i, double *px, const double **Z, const DATAINFO *pdinfo);

int get_panel_structure (const DATAINFO *pdinfo, int *nunits, int *T);

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

int ls_aic_bic (MODEL *pmod);

double gretl_mean (int t1, int t2, const double *x);

double gretl_stddev (int t1, int t2, const double *x);

double gretl_variance (int t1, int t2, const double *x);

double gretl_sst (int t1, int t2, const double *x);

#ifndef WIN32
int gretl_spawn (const char *cmdline);
#endif

#endif /* GRETL_UTILS_H */
