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

enum flagvals {
    OPT_O = 1,
    OPT_M,
    OPT_C,
    OPT_R,
    OPT_S,
    OPT_L,
    OPT_Z,
    OPT_FA,
    OPT_RESID,
    OPT_R_ALT,
    OPT_FAZ,
    OPT_RESIDZ
};

#define free_model(p) if (p != NULL) { \
                             clear_model(p, NULL, NULL); \
                             free(p); \
                          }
#define dataset_add_vars(n, p1, p2) grow_Z((n), p1, p2)
#define dataset_drop_vars(n, p1, p2) shrink_Z((n), p1, p2)
#define dataset_is_time_series(p) (p->time_series == 1)
#define dataset_is_panel(p) (p->time_series == 2 || p->time_series == 3)

#include <float.h>
#define floateq(x, y) (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y) ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y) ((y) - (x) > DBL_EPSILON)
#define na(x) (fabs((x) + 999.0) < DBL_EPSILON)
#define NADBL -999.0

/* functions follow */
 
double date (const int nt, const int pd, const double sd0);

int isdummy (const int varnum, const int t1, const int t2, 
	     const double *Z, const int n);

void printlist (const int *list, const char *msg);

void list_exclude (const int n, int *list);

int set_obs (char *line, DATAINFO *pdinfo, int opt, char *msg);

char *addpath (char *filename, PATHS *ppaths, int script);

int getopenfile (const char *line, char *fname, PATHS *ppaths,
		 int setpath, int script);

char getflag (int opt);

int catchflag (char *line, int *oflag);

MODEL *gretl_model_new (void);

int clear_model (void *ptr, SESSION *psession, session_t *rebuild);

void show_paths (PATHS *ppaths);

int set_paths (PATHS *ppaths, const int defaults, const int gui);

int copylist (int **target, const int *src);

int grow_nobs (const int newobs, double **pZ, DATAINFO *pdinfo);

int grow_Z (const int newvars, double **pZ, DATAINFO *pdinfo);

int shrink_Z (const int delvars, double **pZ, DATAINFO *pdinfo);

int hidden_var (const int i, const DATAINFO *pdinfo);

int copy_model (MODEL *targ, const MODEL *src, const DATAINFO *pdinfo);

int swap_models (MODEL **targ, MODEL **src);

int fcast_with_errs (const char *str, const MODEL *pmod, 
		     DATAINFO *pdinfo, double **pZ, print_t *prn,
		     const PATHS *ppaths, const int plot, char *msg);

int is_model_cmd (const char *line);

int is_model_ref_cmd (const int ci);

int save_model_spec (MODEL *pmod, MODELSPEC *spec, DATAINFO *fullinfo);

int re_estimate (char *model_spec, MODEL *tmpmod, DATAINFO *pdinfo, 
		 double **pZ);

double *copyvec (const double *src, const int n);

int ijton (const int i, const int j, const int lo);

int ztox (const int i, double *px, const DATAINFO *pdinfo, 
	  const double *Z);

#endif /* UTILS_H */
