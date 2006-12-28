/*
 *  Copyright (c) 2005 by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef DATASET_H
#define DATASET_H

typedef enum {
    NO_MARKERS = 0,
    REGULAR_MARKERS,
    DAILY_DATE_STRINGS
} DatasetMarkerType;

typedef enum {
    VAR_DISCRETE   = 1 << 0,
    VAR_SCALAR     = 1 << 1,
    VAR_HIDDEN     = 1 << 2,
    VAR_GENERATED  = 1 << 3,
    VAR_SETCONST   = 1 << 4
} VarinfoFlags;

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
 * dataset_is_seasonal:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains seasonal time series
 * data (1) or not (0).
 */
#define dataset_is_seasonal(p) (((p)->structure == TIME_SERIES || \
                                (p)->structure == SPECIAL_TIME_SERIES) && \
                                (p)->pd > 1)

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
 * quarterly_or_monthly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set is a quarterly
 * or monthly time series (1), or something else (0).
 */
#define quarterly_or_monthly(p) (p->structure == TIME_SERIES && \
                                 (p->pd == 4 || p->pd == 12))

/**
 * dataset_is_panel:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) ((p)->structure == STACKED_TIME_SERIES)

/**
 * var_is_discrete:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether a variable should be treated as discrete
 * or not.
 */
#define var_is_discrete(p, i) ((p)->varinfo[i]->flags & VAR_DISCRETE)

/**
 * var_is_scalar:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable is a scalar.
 */
#define var_is_scalar(p, i) ((p)->varinfo[i]->flags & VAR_SCALAR)

/**
 * var_is_series:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable is a series (as opposed
 * to a scalar).
 */
#define var_is_series(p, i) (!((p)->varinfo[i]->flags & VAR_SCALAR))

/**
 * var_is_hidden:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable is hidden.
 */
#define var_is_hidden(p, i) ((p)->varinfo[i]->flags & VAR_HIDDEN)

/**
 * var_is_generated:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable was generated using
 * a formula or transformation function.
 */
#define var_is_generated(p, i) ((p)->varinfo[i]->flags & VAR_GENERATED)

/**
 * var_is_const:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable has been marked as
 * "const".
 */
#define var_is_const(p, i) (i == 0 || ((p)->varinfo[i]->flags & VAR_SETCONST))

/**
 * set_var_const:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Set the "const" flag on the given variable.
 */
#define set_var_const(p, i) ((p)->varinfo[i]->flags |= VAR_SETCONST)

/**
 * unset_var_const:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Remove the "const" flag from the given variable.
 */
#define unset_var_const(p, i) ((p)->varinfo[i]->flags &= ~VAR_SETCONST)


void free_Z (double **Z, DATAINFO *pdinfo);

DATAINFO *datainfo_new (void);

DATAINFO *create_new_dataset (double ***pZ, /* data matrix */
			      int nvar,     /* number of variables */
			      int nobs,     /* observations per variable */
			      int markers   /* case markers or not? */
			      );

void destroy_dataset (double **Z, DATAINFO *pdinfo);

void clear_datainfo (DATAINFO *pdinfo, int code);

int allocate_Z (double ***pZ, const DATAINFO *pdinfo);

int dataset_allocate_varnames (DATAINFO *pdinfo);

int dataset_allocate_obs_markers (DATAINFO *pdinfo);

void dataset_destroy_obs_markers (DATAINFO *pdinfo);

int dataset_allocate_panel_info (DATAINFO *pdinfo);

void dataset_destroy_panel_info (DATAINFO *pdinfo);

int dataset_add_default_panel_indices (DATAINFO *pdinfo);

int dataset_finalize_panel_indices (DATAINFO *pdinfo);

void dataset_obs_info_default (DATAINFO *pdinfo);

void copy_dataset_obs_info (DATAINFO *targ, const DATAINFO *src);

void copy_varinfo (VARINFO *targ, const VARINFO *src);

void set_sorted_markers (DATAINFO *pdinfo, int v, char **S);

void dataset_set_regular_markers (DATAINFO *pdinfo);

int start_new_Z (double ***pZ, DATAINFO *pdinfo, int resample);

int is_trend_variable (const double *x, int n);

int is_periodic_dummy (const double *x, const DATAINFO *pdinfo);

int dataset_add_observations (int newobs, double ***pZ, DATAINFO *pdinfo,
			      gretlopt opt);

int dataset_drop_observations (int n, double ***pZ, DATAINFO *pdinfo);

int dataset_shrink_obs_range (double ***pZ, DATAINFO *pdinfo);

int dataset_add_series (int newvars, double ***pZ, DATAINFO *pdinfo);

int dataset_add_allocated_series (double *x, double ***pZ, 
				  DATAINFO *pdinfo);

int dataset_add_scalars (int n, double ***pZ, DATAINFO *pdinfo);

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo);

int dataset_add_scalar_as (double x, const char *newname,
			   double ***pZ, DATAINFO *pdinfo);

int dataset_add_series_as (double *x, const char *newname,
			   double ***pZ, DATAINFO *pdinfo);

int dataset_copy_variable_as (int v, const char *newname,
			      double ***pZ, DATAINFO *pdinfo);

int overwrite_err (const DATAINFO *pdinfo, int v);

int dataset_drop_listed_variables (int *list, double ***pZ, 
				   DATAINFO *pdinfo, int *renumber);

int dataset_drop_variable (int v, double ***pZ, DATAINFO *pdinfo); 

int dataset_destroy_hidden_variables (double ***pZ, DATAINFO *pdinfo,
				      int vmin);

int dataset_drop_last_variables (int delvars, double ***pZ, DATAINFO *pdinfo);

int dataset_stack_variables (const char *vname, const char *line,
			     double ***pZ, DATAINFO *pdinfo, 
			     PRN *prn);

int is_log_variable (int i, const DATAINFO *pdinfo, char *parent);

void set_var_discrete (DATAINFO *pdinfo, int i, int s);

void set_var_scalar (DATAINFO *pdinfo, int i, int s);

void set_var_hidden (DATAINFO *pdinfo, int i);

void var_set_linewidth (DATAINFO *pdinfo, int i, int w);

int var_get_linewidth (const DATAINFO *pdinfo, int i);

#endif /* DATASET_H */
