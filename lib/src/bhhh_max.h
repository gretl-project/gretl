/* 
 * Copyright (C) 2004 Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
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

typedef struct _model_info model_info;

struct _model_info {
    int k;
    int n;
    int n_series;
    double ll;
    double ll2;
    int *list;
    double *theta;
    double *delta;
    double *deltmp;
    double **series;
};

int bhhh_max (int (*loglik) (double *, const double **, double **,
			     model_info *, int), 
	      const double **X, const double *init_coeff,
	      int n_init_coeff, int n_series, int nobs,
	      double *ll, double *params, gretl_matrix **pVCV,
	      PRN *prn);
