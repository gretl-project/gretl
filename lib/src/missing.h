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

#ifndef MISSING_H
#define MISSING_H

#include <float.h>

#define NEW_NA 0

#define NADBL DBL_MAX
#define na(x) ((x) == NADBL)

/* xna = "extended NA", including regular NA for missing data
   as well as NaNs and infinities */

#ifndef isfinite
# define isfinite(x) (!isnan(x) && !isinf(x))
# define xna(x) ((x) == NADBL || isnan(x) || isinf(x))
#else
# define xna(x) ((x) == NADBL || !isfinite(x))
#endif

#define model_missing(m,t) ((m)->missmask != NULL && (m)->missmask[t] == '1')

int model_has_missing_obs (const MODEL *pmod);

int first_missing_index (const double *x, int t1, int t2);

int series_adjust_sample (const double *x, int *t1, int *t2);

int list_adjust_sample (const int *list, int *t1, int *t2, 
			const DATASET *dset, int *nmiss);

int set_miss (const int *list, const char *param,
	      DATASET *dset, PRN *prn);

double missing_obs_fraction (const DATASET *dset);

int any_missing_user_values (const DATASET *dset);

#endif /* MISSING_H */
