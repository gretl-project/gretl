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

#ifndef LAGSELECT_H
#define LAGSELECT_H

enum {
    LAG_X = LISTSEP + 1, /* lags set for regular variable context */
    LAG_Y_X,             /* lags for dependent variable */
    LAG_W,               /* lags set for variable as instrument */
    LAG_Y_W,             /* lags for dependent var as instrument */
    LAG_Y_V              /* lags of endoenous vars in VAR */
} LagContext;

#define VDEFLT -1

void enable_lags_for_context (int context, gboolean s);

gboolean lags_enabled_for_context (int context);

void destroy_lag_preferences (void);

int set_lag_prefs_from_list (int v, int *llist, char context,
			     int *changed);

int set_lag_prefs_from_minmax (int v, int lmin, int lmax,
			       char context, int *changed);

void set_null_lagpref (int v, char context, int *changed);

int set_lag_prefs_from_model (int dv, int *xlist, int *zlist);

int set_lag_prefs_from_VAR (const int *lags, int *xlist);

int *get_lag_pref_as_list (int v, char context);

int remove_specific_lag (int v, int lag, char context);

int is_lag_dummy (int v, int lag, char context);

const int *get_VAR_lags_list (void);

void set_VAR_max_lag (int lmax);

void get_lag_preference (int v, int *lmin, int *lmax, const int **laglist,
			 char context, selector *sr);

#endif
