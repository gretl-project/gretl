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
    VAR_LISTARG    = 1 << 3,
    VAR_TIMECOL    = 1 << 4
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
    DS_CLEAR,
    DS_RENUMBER,
    DS_INSOBS
} DatasetOp;

typedef enum {
    DS_COPY_VALUES,
    DS_GRAB_VALUES
} DataCopyFlag;

/**
 * CompactMethod:
 * @COMPACT_NONE: no data compaction
 * @COMPACT_SUM: take sum of higher frequency data
 * @COMPACT_AVG: take mean of higher frequency data
 * @COMPACT_SOP: use start-of-period value
 * @COMPACT_EOP: use end-of-period value
 * @COMPACT_WDAY: use a specified day of the week
 * @COMPACT_MAX: sentinel value
 *
 * Symbolic codes for various methods of compacting data
 * series (i.e. converting from a higher to a lower
 * frequency). %COMPACT_WDAY is applicable only when
 * converting from daily to weekly frequency.
 */

typedef enum {
    COMPACT_NONE,
    COMPACT_SUM,
    COMPACT_AVG,
    COMPACT_SOP,
    COMPACT_EOP,
    COMPACT_WDAY,
    COMPACT_MAX
} CompactMethod;

typedef struct _series_table series_table;

/**
 * dataset_is_cross_section:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains cross-sectional
 * data (1) or not (0).
 */
#define dataset_is_cross_section(p) (p != NULL && p->structure == CROSS_SECTION)

/**
 * dataset_is_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains time series
 * data (1) or not (0).
 */
#define dataset_is_time_series(p) (p != NULL && (p->structure == TIME_SERIES || \
						 p->structure == SPECIAL_TIME_SERIES))

/**
 * dataset_is_seasonal:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains seasonal time series
 * data (1) or not (0).
 */
#define dataset_is_seasonal(p) (p != NULL && (p->structure == TIME_SERIES || \
                                p->structure == SPECIAL_TIME_SERIES) && \
                                p->pd > 1)

/**
 * custom_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains time series
 * data with custom (non-standard) frequency (1) or not (0).
 */
#define custom_time_series(p) (p != NULL && p->structure == SPECIAL_TIME_SERIES)

/**
 * dataset_is_daily:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains daily time series
 * data (1) or not (0).
 */
#define dataset_is_daily(p) (p != NULL && p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7))

/**
 * dataset_is_weekly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains weekly time series
 * data (1) or not (0).
 */
#define dataset_is_weekly(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 52)

/**
 * dataset_is_hourly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains hourly time series
 * data (1) or not (0).
 */
#define dataset_is_hourly(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 24)

/**
 * dataset_is_decennial:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains decennial time series
 * data (1) or not (0).
 */
#define dataset_is_decennial(p) (p != NULL && p->structure == TIME_SERIES \
                                 && p->pd == 10)

/**
 * dated_daily_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p != NULL && p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                             && p->sd0 > 10000.0)

/**
 * dated_seven_day_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains dated daily 
 * (seven-day) time series data (1) or not (0).
 */
#define dated_seven_day_data(p) (p != NULL && p->structure == TIME_SERIES \
                                 && p->pd == 7 && \
                                 p->sd0 > 10000.0)

/**
 * dated_weekly_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains dated weekly 
 * time series data (1) or not (0).
 */
#define dated_weekly_data(p) (p != NULL && p->structure == TIME_SERIES \
                              && p->pd == 52 && \
                              p->sd0 > 10000.0)

/**
 * calendar_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset uses calendar
 * dates for observation strings (1) or not (0).
 */
#define calendar_data(p) (p != NULL && p->structure == TIME_SERIES && \
                          (p->pd == 5 || p->pd == 6 || p->pd == 7 \
                           || p->pd == 52) && p->sd0 > 10000.0) 

/**
 * quarterly_or_monthly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset is a quarterly
 * or monthly time series (1), or something else (0).
 */
#define quarterly_or_monthly(p) (p != NULL && p->structure == TIME_SERIES && \
                                 (p->pd == 4 || p->pd == 12))

/**
 * annual_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset is an annual
 * time series (1), or something else (0).
 */
#define annual_data(p) (p != NULL && p->structure == TIME_SERIES && \
			p->pd == 1)

/**
 * decennial_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset is a decemmial
 * time series (1), or something else (0).
 */
#define decennial_data(p) (p != NULL && p->structure == TIME_SERIES && \
			   p->pd == 10 && p->sd0 > 1000)

/**
 * dataset_is_panel:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a dataset contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) (p != NULL && p->structure == STACKED_TIME_SERIES)

/**
 * dataset_has_markers:
 * @p: pointer to data information struct.
 *
 * Determine whether a dataset has observation marker strings (1)
 * or not (0).
 */
#define dataset_has_markers(p) (p != NULL && p->markers && p->S != NULL)

/**
 * dataset_has_panel_time:
 * @p: pointer to data information struct.
 *
 * Determine whether a panel dataset has information on its time
 * dimension recorded (1) or not (0).
 */
#define dataset_has_panel_time(p) (p != NULL && \
				   p->structure == STACKED_TIME_SERIES && \
				   p->panel_pd > 0 && p->panel_sd0 > 0.0)

/**
 * sample_size:
 * @p: pointer to data information struct.
 *
 * Retrieves the length of the current sample range.
 */
#define sample_size(p) ((p == NULL)? 0 : (p->t2 - p->t1 + 1))

/**
 * dset_get_data:
 * @d: pointer to dataset.
 * @i: index number of variable.
 * @t: observation number.
 *
 * Gets the value of series @i at observation @t.
 */
#define dset_get_data(d,i,t) (d->Z[i][t])

/**
 * dset_set_data:
 * @d: pointer to dataset.
 * @i: index number of variable.
 * @t: observation number.
 * @x: value to set.
 *
 * Sets the value of series @i at observation @t.
 */
#define dset_set_data(d,i,t,x) (d->Z[i][t]=x)

void free_Z (DATASET *dset);

DATASET *datainfo_new (void);

void datainfo_init (DATASET *dset);

DATASET *create_new_dataset (int nvar,     /* number of variables */
			     int nobs,     /* observations per variable */
			     int markers   /* case markers or not? */
			     );

DATASET *create_auxiliary_dataset (int nvar, int nobs, gretlopt opt);

void destroy_dataset (DATASET *dset);

void clear_datainfo (DATASET *dset, int code);

int allocate_Z (DATASET *dset, gretlopt opt);

int dataset_allocate_varnames (DATASET *dset);

int dataset_allocate_obs_markers (DATASET *dset);

void dataset_destroy_obs_markers (DATASET *dset);

void dataset_obs_info_default (DATASET *dset);

void copy_dataset_obs_info (DATASET *targ, const DATASET *src);

void copy_varinfo (VARINFO *targ, const VARINFO *src);

void set_sorted_markers (DATASET *dset, int v, char **S);

void dataset_set_regular_markers (DATASET *dset);

int start_new_Z (DATASET *dset, gretlopt opt);

int is_trend_variable (const double *x, int n);

int is_periodic_dummy (const double *x, const DATASET *dset);

int dataset_add_observations (DATASET *dset, int n, gretlopt opt);

int dataset_drop_observations (DATASET *dset, int n);

int dataset_shrink_obs_range (DATASET *dset);

int dataset_add_series (DATASET *dset, int newvars);

int dataset_add_NA_series (DATASET *dset);

int dataset_add_allocated_series (DATASET *dset, double *x);

int dataset_add_series_as (DATASET *dset, double *x, const char *name);

int dataset_copy_series_as (DATASET *dset, int v, const char *name);

int overwrite_err (const char *name);

int series_is_parent (const DATASET *dset, int v);

int dataset_replace_series (DATASET *dset, int v,
			    double *x, const char *descrip,
			    DataCopyFlag flag);

int dataset_rename_series (DATASET *dset, int v, const char *name);

int dataset_drop_listed_variables (int *list, DATASET *dset, 
				   int *renumber, PRN *prn);

int dataset_drop_variable (int v, DATASET *dset); 

int dataset_destroy_hidden_variables (DATASET *dset, int vmin);

int dataset_drop_last_variables (DATASET *dset, int delvars);

int dataset_renumber_variable (int v_old, int v_new, 
			       DATASET *dset);

int renumber_series_with_checks (const char *s, int fixmax,
				 DATASET *dset, PRN *prn);

int maybe_prune_dataset (DATASET **pdset, void *p);

int dataset_stack_variables (const char *vname, const char *line,
			     DATASET *dset, PRN *prn);

int dataset_sort_by (DATASET *dset, const int *list, gretlopt opt);

int dataset_set_matrix_name (DATASET *dset, const char *name);

const char *dataset_get_matrix_name (const DATASET *dset);

const char *dataset_period_label (const DATASET *dset);

int series_is_log (const DATASET *dset, int i, char *parent);

void series_set_discrete (DATASET *dset, int i, int s);

int series_record_display_name (DATASET *dset, int i,
				const char *s); 

int series_record_label (DATASET *dset, int i,
			 const char *s); 

const char *series_get_graph_name (const DATASET *dset, int i);

unsigned int get_resampling_seed (void);

int dataset_resample (DATASET *dset, int n, unsigned int seed);

int dataset_op_from_string (const char *s);

int modify_dataset (DATASET *dset, int op, const int *list, 
		    const char *s, PRN *prn);

int dataset_get_structure (const DATASET *dset);

int panel_sample_size (const DATASET *dset);

int multi_unit_panel_sample (const DATASET *dset);

int dataset_purge_missing_rows (DATASET *dset);

int check_dataset_is_changed (void);

void set_dataset_is_changed (void);

int dataset_set_time_series (DATASET *dset, int pd, 
			     int yr0, int minor0);

int series_is_discrete (const DATASET *dset, int i);

int series_is_hidden (const DATASET *dset, int i);

int series_is_generated (const DATASET *dset, int i);

int series_is_listarg (const DATASET *dset, int i);

int series_get_flags (DATASET *dset, int i);

void series_set_flag (DATASET *dset, int i, int flag);

void series_unset_flag (DATASET *dset, int i, int flag);

void series_zero_flags (DATASET *dset, int i);

const char *series_get_label (const DATASET *dset, int i);

const char *series_get_display_name (const DATASET *dset, int i);

const char *series_get_parent_name (const DATASET *dset, int i);

int series_get_parent_id (const DATASET *dset, int i);

int series_get_compact_method (const DATASET *dset, int i);

int series_get_stack_level (const DATASET *dset, int i);

int series_get_transform (const DATASET *dset, int i);

int series_get_lag (const DATASET *dset, int i);

void series_set_label (DATASET *dset, int i, 
		       const char *s);

void series_set_display_name (DATASET *dset, int i, 
			      const char *s);

void series_set_compact_method (DATASET *dset, int i, 
				int method);

void series_set_parent (DATASET *dset, int i, 
			const char *parent);

void series_set_transform (DATASET *dset, int i, 
			   int transform);

void series_set_lag (DATASET *dset, int i, int lag);

void series_set_stack_level (DATASET *dset, int i, int level);

void series_increment_stack_level (DATASET *dset, int i);

void series_decrement_stack_level (DATASET *dset, int i);

void series_attach_string_table (DATASET *dset, int i, void *ptr);

void series_destroy_string_table (DATASET *dset, int i);

int is_string_valued (const DATASET *dset, int i);

series_table *series_get_string_table (const DATASET *dset, int i);

const char *series_get_string_for_obs (const DATASET *dset, int i, 
				       int t);

const char *series_get_string_for_value (const DATASET *dset, int i,
					 double val);

int series_set_string_val (DATASET *dset, int i, int t, const char *s);

double series_decode_string (const DATASET *dset, int i, const char *s);

const char **series_get_string_vals (const DATASET *dset, int i,
				     int *n_strs);

int steal_string_table (DATASET *l_dset, int lvar,
			DATASET *r_dset, int rvar);

int set_panel_groups_name (DATASET *dset, const char *vname);

char const **get_panel_group_names (DATASET *dset);

int panel_group_names_ok (const DATASET *dset);

const char *panel_group_names_varname (const DATASET *dset);

int is_panel_group_names_series (const DATASET *dset, int v);

#endif /* DATASET_H */
