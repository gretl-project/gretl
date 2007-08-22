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

int model_missval_count (const MODEL *pmod);

int list_adjust_t1t2 (const int *list, const double **Z, 
		      DATAINFO *pdinfo);

int array_adjust_t1t2 (const double *x, int *t1, int *t2);

int varlist_adjust_sample (const int *list, int *t1, int *t2, 
			   const double **Z);

int set_reference_missmask_from_list (const int *list, const double **Z,
				      const DATAINFO *pdinfo);

int check_for_missing_obs (const int *list, int *t1, int *t2,
			   const double **Z, int *misst);

int set_miss (const int *list, const char *param, double **Z,
	      DATAINFO *pdinfo, PRN *prn);

#endif /* MISSING_H */
