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

/* utils.h for gretl */

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

/**
 * free_model:
 * @p: pointer to #MODEL.
 *
 * Free allocated content of MODEL then the pointer itself.
 */
#define free_model(p) if (p != NULL) { \
                             clear_model(p, NULL); \
                             free(p); \
                          }

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
                             && (p->pd == 5 || p->pd == 7))

/**
 * dated_daily_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p->time_series == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 7) && \
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

int isdummy (double *x, int t1, int t2);

void printlist (const int *list, const char *msg);

int print_list_to_buffer (const int *list, char *buf, size_t len);

void list_exclude (int n, int *list);

int set_obs (char *line, DATAINFO *pdinfo, unsigned char opt);

void set_miss (LIST list, const char *param, double **Z,
	       DATAINFO *pdinfo, PRN *prn);

char *addpath (char *fname, PATHS *ppaths, int script);

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 int setpath, int script);

int catchflag (char *line, unsigned char *oflag);

MODEL *gretl_model_new (const DATAINFO *pdinfo);

void exchange_smpl (MODEL *pmod, DATAINFO *pdinfo);

void clear_model (MODEL *pmod, const DATAINFO *pdinfo);

void debug_print_model_info (const MODEL *pmod, const char *msg);

void show_paths (PATHS *ppaths);

int set_paths (PATHS *ppaths, int defaults, int gui);

const char *fetch_gretl_lib_path (void);

int copylist (int **target, const int *src);

int grow_nobs (int newobs, double ***pZ, DATAINFO *pdinfo);

int dataset_add_vars (int newvars, double ***pZ, DATAINFO *pdinfo);

int dataset_add_allocated_var (double *x, double ***pZ, DATAINFO *pdinfo);

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo);

int varnum_from_string (const char *str, DATAINFO *pdinfo);

int dataset_drop_listed_vars (const int *list, double ***pZ, 
			      DATAINFO *pdinfo, int *renumber);

int dataset_drop_var (int varno, double ***pZ, DATAINFO *pdinfo);

int dataset_drop_vars (int delvars, double ***pZ, DATAINFO *pdinfo);

int rename_var_by_id (const char *str, const char *vname, 
		      DATAINFO *pdinfo);

int hidden_var (int i, const DATAINFO *pdinfo);

int copy_model (MODEL *targ, const MODEL *src, const DATAINFO *pdinfo);

int swap_models (MODEL **targ, MODEL **src);

FITRESID *get_fit_resid (const MODEL *pmod, double ***pZ, 
			 DATAINFO *pdinfo);

FITRESID *get_fcast_with_errs (const char *str, const MODEL *pmod, 
			       double ***pZ, DATAINFO *pdinfo, PRN *prn);

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, PRN *prn,
		     PATHS *ppaths, int plot);

int is_model_cmd (const char *line);

int is_model_ref_cmd (int ci);

int save_model_spec (MODEL *pmod, MODELSPEC *spec, DATAINFO *fullinfo);

int re_estimate (char *model_spec, MODEL *tmpmod, 
		 double ***pZ, DATAINFO *pdinfo);

double *copyvec (const double *src, int n);

int ijton (int i, int j, int nrows);

int ztox (int i, double *px, 
	  double **Z, const DATAINFO *pdinfo);

int get_panel_structure (DATAINFO *pdinfo, int *nunits, int *T);

int set_panel_structure (unsigned char flag, DATAINFO *pdinfo, PRN *prn); 

int balanced_panel (const DATAINFO *pdinfo);

double get_xvalue (int i, double **Z, const DATAINFO *pdinfo);

void free_gretl_mp_results (mp_results *mpvals);

mp_results *gretl_mp_results_new (int totvar);

int allocate_mp_varnames (mp_results *mpvals);

FITRESID *fit_resid_new (int n, int errs);

void free_fit_resid (FITRESID *fr);

CONFINT *get_model_confints (const MODEL *pmod);

void free_confint (CONFINT *cf);

#ifndef WIN32
int gretl_spawn (const char *cmdline);
#endif

#endif /* UTILS_H */
