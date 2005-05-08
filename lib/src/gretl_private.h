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

/* functions shared internally by library translation units */

#ifndef GRETL_PRIVATE_H
#define GRETL_PRIVATE_H

extern int newlag; /* transforms.c */

int dataset_stack_vars (double ***pZ, DATAINFO *pdinfo, 
			char *newvar, char *s);

void gretl_printxs (double xx, int n, int ci, PRN *prn);

void bufspace (int n, PRN *prn);

void gretl_print_ar (MODEL *pmod, PRN *prn);

int gretl_criteria (double ess, int nobs, int ncoeff, PRN *prn);

int calculate_criteria (double *x, double ess, int nobs, int ncoeff);

int get_t_from_obs_string (char *s, const double **Z, 
			   const DATAINFO *pdinfo);

int model_mask_leaves_balanced_panel (const MODEL *pmod,
				      const DATAINFO *pdinfo);

int gretl_forecast (int t1, int t2, int nv, 
		    const MODEL *pmod, double ***pZ);

int gretl_is_reserved (const char *str);

int real_list_laggenr (const int *list, double ***pZ, 
		       DATAINFO *pdinfo, int maxlag, 
		       int **lagnums);

int lagvarnum (int v, int l, const DATAINFO *pdinfo);

int get_genr_function (const char *s);

char *copy_submask (const char *src, int n);

int get_vcv_index (MODEL *pmod, int i, int j, int n);

int takenotes (int quit_opt);

char *get_month_name (char *mname, int m);

double *testvec (int n);

#endif /* GRETL_PRIVATE_H */
