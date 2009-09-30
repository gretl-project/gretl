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

#ifndef DATASET_H
#define DATASET_H

typedef enum {
    DATA_NONE,
    DATA_XSECT,
    DATA_TS,
    DATA_PANEL
} DatasetStructure;

typedef enum {
    NO_MARKERS = 0,
    REGULAR_MARKERS,
    DAILY_DATE_STRINGS
} DatasetMarkerType;

typedef enum {
    VAR_DISCRETE   = 1 << 0,
    VAR_HIDDEN     = 1 << 1,
    VAR_GENERATED  = 1 << 2,
    VAR_LISTARG    = 1 << 3
} VarinfoFlags;

typedef enum {
    DS_NONE,
    DS_ADDOBS,
    DS_COMPACT,
    DS_EXPAND,
    DS_TRANSPOSE,
    DS_DELETE,
    DS_KEEP,
    DS_SORTBY,
    DS_DSORTBY,
    DS_RESAMPLE,
    DS_RESTORE,
    DS_CLEAR
} DatasetOp;

/**
 * dataset_is_cross_section:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains cross-sectional
 * data (1) or not (0).
 */
#define dataset_is_cross_section(p) (p != NULL && p->structure == CROSS_SECTION)

/**
 * dataset_is_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data (1) or not (0).
 */
#define dataset_is_time_series(p) (p != NULL && (p->structure == TIME_SERIES || \
						 p->structure == SPECIAL_TIME_SERIES))

/**
 * dataset_is_seasonal:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains seasonal time series
 * data (1) or not (0).
 */
#define dataset_is_seasonal(p) (p != NULL && (p->structure == TIME_SERIES || \
                                p->structure == SPECIAL_TIME_SERIES) && \
                                p->pd > 1)

/**
 * custom_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data with custom (non-standard) frequency (1) or not (0).
 */
#define custom_time_series(p) (p != NULL && p->structure == SPECIAL_TIME_SERIES)

/**
 * dataset_is_daily:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains daily time series
 * data (1) or not (0).
 */
#define dataset_is_daily(p) (p != NULL && p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7))

/**
 * dataset_is_weekly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains weekly time series
 * data (1) or not (0).
 */
#define dataset_is_weekly(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 52)

/**
 * dataset_is_hourly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains hourly time series
 * data (1) or not (0).
 */
#define dataset_is_hourly(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 24)

/**
 * dataset_is_decennial:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains decennial time series
 * data (1) or not (0).
 */
#define dataset_is_decennial(p) (p != NULL && p->structure == TIME_SERIES \
                                 && p->pd == 10)

/**
 * dated_daily_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p != NULL && p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                             && p->sd0 > 10000.0)

/**
 * dated_seven_day_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily 
 * (seven-day) time series data (1) or not (0).
 */
#define dated_seven_day_data(p) (p != NULL && p->structure == TIME_SERIES \
                                 && p->pd == 7 && \
                                 p->sd0 > 10000.0)

/**
 * dated_weekly_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated weekly 
 * time series data (1) or not (0).
 */
#define dated_weekly_data(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 52 && \
                              p->sd0 > 10000.0)

/**
 * calendar_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set uses calendar
 * dates for observation strings (1) or not (0).
 */
#define calendar_data(p) (p != NULL && p->structure == TIME_SERIES && \
                          (p->pd == 5 || p->pd == 6 || p->pd == 7 \
                           || p->pd == 52) && p->sd0 > 10000.0) 

/**
 * quarterly_or_monthly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set is a quarterly
 * or monthly time series (1), or something else (0).
 */
#define quarterly_or_monthly(p) (p != NULL && p->structure == TIME_SERIES && \
                                 (p->pd == 4 || p->pd == 12))

/**
 * dataset_is_panel:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) (p != NULL && p->structure == STACKED_TIME_SERIES)

/**
 * var_is_discrete:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether a variable should be treated as discrete
 * or not.
 */
#define var_is_discrete(p, i) (p->varinfo[i]->flags & VAR_DISCRETE)

/**
 * var_is_hidden:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable is hidden.
 */
#define var_is_hidden(p, i) (p->varinfo[i]->flags & VAR_HIDDEN)

/**
 * var_is_generated:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable was generated using
 * a formula or transformation function.
 */
#define var_is_generated(p, i) (p->varinfo[i]->flags & VAR_GENERATED)

/**
 * var_is_listarg:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Determine whether or not a variable has been marked as
 * belonging to a list argument to a function.
 */
#define var_is_listarg(p, i) (p->varinfo[i]->flags & VAR_LISTARG)

/**
 * set_var_listarg:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Set the "listarg" flag on the given variable.
 */
#define set_var_listarg(p, i) (p->varinfo[i]->flags |= VAR_LISTARG)

/**
 * unset_var_listarg:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 *
 * Remove the "listarg" flag from the given variable.
 */
#define unset_var_listarg(p, i) (p->varinfo[i]->flags &= ~VAR_LISTARG)

/**
 * series_set_flag:
 * @p: pointer to data information struct.
 * @i: index number of variable.
 * @f: flag to set.
 *
 * Set the given flag on the given (series) variable.
 */
#define series_set_flag(p, i, f) (p->varinfo[i]->flags |= f)

#define sample_size(p) ((p == NULL)? 0 : (p->t2 - p->t1 + 1))

void free_Z (double **Z, DATAINFO *pdinfo);

DATAINFO *datainfo_new (void);

DATAINFO *create_new_dataset (double ***pZ, /* data matrix */
			      int nvar,     /* number of variables */
			      int nobs,     /* observations per variable */
			      int markers   /* case markers or not? */
			      );

DATAINFO *
create_auxiliary_dataset (double ***pZ, int nvar, int nobs);

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

int dataset_add_series_as (double *x, const char *newname,
			   double ***pZ, DATAINFO *pdinfo);

int dataset_copy_variable_as (int v, const char *newname,
			      double ***pZ, DATAINFO *pdinfo);

int overwrite_err (const char *name);

int series_is_parent (const DATAINFO *pdinfo, int v);

int dataset_rename_series (DATAINFO *pdinfo, int v, const char *name);

int dataset_drop_listed_variables (int *list, double ***pZ, 
				   DATAINFO *pdinfo, int *renumber,
				   PRN *prn);

int dataset_drop_variable (int v, double ***pZ, DATAINFO *pdinfo); 

int dataset_destroy_hidden_variables (double ***pZ, DATAINFO *pdinfo,
				      int vmin);

int dataset_drop_last_variables (int delvars, double ***pZ, DATAINFO *pdinfo);

int maybe_prune_dataset (double ***pZ, DATAINFO **ppdinfo, void *p);

int dataset_stack_variables (const char *vname, const char *line,
			     double ***pZ, DATAINFO *pdinfo, 
			     PRN *prn);

int dataset_sort_by (int v, double **Z, DATAINFO *pdinfo, gretlopt opt);

int is_log_variable (int i, const DATAINFO *pdinfo, char *parent);

void set_var_discrete (DATAINFO *pdinfo, int i, int s);

void set_var_scalar (DATAINFO *pdinfo, int i, int s);

void set_var_hidden (DATAINFO *pdinfo, int i);

void var_set_linewidth (DATAINFO *pdinfo, int i, int w);

int var_get_linewidth (const DATAINFO *pdinfo, int i);

int var_set_display_name (DATAINFO *pdinfo, int i,
			  const char *s); 

int var_set_description (DATAINFO *pdinfo, int i,
			 const char *s); 

int var_set_compact_method (DATAINFO *pdinfo, int i,
			    int method);

const char *var_get_graph_name (const DATAINFO *pdinfo, int i);

unsigned int get_resampling_seed (void);

int dataset_resample (int n, unsigned int seed,
		      double ***pZ, DATAINFO *pdinfo);

int dataset_op_from_string (const char *s);

int modify_dataset (int op, const int *list, const char *s, 
		    double ***pZ, DATAINFO *pdinfo, 
		    PRN *prn);

int dataset_get_structure (const DATAINFO *pdinfo);

int dataset_purge_missing_rows (double **Z, DATAINFO *pdinfo);

int check_dataset_is_changed (void);

void set_dataset_is_changed (void);

#endif /* DATASET_H */
