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

int laggenr (int parent, int lag, int opt, double ***pZ, 
	     DATAINFO *pdinfo);

int gretl_multiply (char *s, int *list, char *sfx, double ***pZ,
		    DATAINFO *pdinfo);

void gretl_print_add (const COMPARE *add, const int *addvars, 
		      const DATAINFO *pdinfo, int aux_code, PRN *prn,
		      gretlopt opt);

void gretl_print_omit (const COMPARE *omit, const int *omitvars, 
		       const DATAINFO *pdinfo, PRN *prn,
		       gretlopt opt);

void graphyzx (const int *list, const double *zy1, const double *zy2, 
	       const double *zx, int n, const char *yname, 
	       const char *xname, const DATAINFO *pdinfo, 
	       gretlopt oflag, PRN *prn);

void gretl_printxs (double xx, int n, int ci, PRN *prn);

void bufspace (int n, PRN *prn);

void gretl_print_ar (MODEL *pmod, PRN *prn);

void gretl_delete (char *str, int indx, int count);

int numeric_string (const char *str);

void gretl_trunc (char *str, size_t n);

int count_fields (const char *str);

void shift_left (char *str, size_t move);

double gretl_corr (int n, const double *zx, const double *zy);

double gretl_covar (int n, const double *zx, const double *zy);

int gretl_iszero (int t1, int t2, const double *x);

int gretl_isconst (int t1, int t2, const double *x);

double gretl_mean (int t1, int t2, const double *x);

void gretl_minmax (int t1, int t2, const double zx[], 
	      double *min, double *max);

int gretl_hasconst (const int *list);

int gretl_compare_doubles (const void *a, const void *b);

double gretl_stddev (int t1, int t2, const double *x);

double gretl_variance (int t1, int t2, const double *x);

double gretl_sst (int t1, int t2, const double *x);

void gretl_aic_etc (MODEL *pmod);

void gretl_criteria (double ess, int nobs, int ncoeff, PRN *prn);

int adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		 const double **Z, int *misst);

int list_dups (const int *list, int ci);

int gretl_forecast (int t1, int t2, int nv, 
		    const MODEL *pmod, double ***pZ);

int z_to_xy (int v1, int v2, double *px, double *py, 
	     const DATAINFO *pdinfo, double **Z);

int gretl_is_reserved (const char *str);

void gretl_test_init (GRETLTEST *test);

void rearrange_list (int *list);

int vars_identical (const double *x, const double *y, int n);

int get_function (const char *s);

void gretl_varinfo_init (VARINFO *vinfo);

double corrrsq (int nobs, const double *y, const double *yhat);

int *big_list (const int *orig, const int *add, const DATAINFO *pdinfo, 
	       int model_count, int *err);

int get_hac_lag (int m);

int get_hc_version (void);

int get_use_qr (void);

char *copy_subdum (const char *src, int n);

int get_vcv_index (MODEL *pmod, int i, int j, int n);

int path_append (char *file, const char *path);

#endif /* GRETL_PRIVATE_H */
