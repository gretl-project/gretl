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

#ifndef GRETL_INTERNAL_H
#define GRETL_INTERNAL_H

int laggenr (const int iv, const int lag, const int opt, double **pZ, 
	     DATAINFO *pdinfo);

int multiply (char *s, int *list, char *sfx, double **pZ,
	      DATAINFO *pdinfo);

void print_add (const COMPARE *add, const int *addvars, 
		const DATAINFO *pdinfo, const int aux_code, print_t *prn);

void print_omit (const COMPARE *omit, const int *omitvars, 
		 const DATAINFO *pdinfo, print_t *prn);

void print_aicetc (const MODEL *pmod, print_t *prn);

void printcorr (int *list, const CORRMAT corrmat, 
		const DATAINFO *pdinfo, print_t *prn);

void graphyzx (const int *list, const double *zy1, const double *zy2, 
	       const double *zx, int n, const char *yname, 
	       const char *xname, const DATAINFO *pdinfo, 
	       int oflag, print_t *prn);

void printxs (double xx, int n, int ci, print_t *prn);

void print_ar (MODEL *pmod, print_t *prn);

void delete (char *str, const int indx, const int count);

int isnumber (const char *str);

void esl_trunc (char *str, const int n);

int count_fields (const char *str);

void shiftleft (char *str, size_t move);

double corr (const int n, const double *zx, const double *zy);

double covar (const int n, const double *zx, const double *zy);

int ijton (const int i, const int j, const int lo);

int iszero (const int t1, const int t2, const double *x);

void list_exclude (const int n, int *list);

int isconst (const int t1, const int t2, const double *x);

double esl_mean (const int t1, const int t2, const double *x);

void minmax (const int t1, const int t2, const double zx[], 
	     double *min, double *max);

int hasconst (const int *list);

int compare_doubles (const void *a, const void *b);

double esl_stddev (const int t1, const int t2, const double *x);

double esl_variance (const int t1, const int t2, const double *x);

void aicetc (MODEL *pmod);

void criteria (const double ess, const int nobs, const int ncoeff, 
	       print_t *prn);

CORRMAT corrlist (int *list, const double *Z, const DATAINFO *pdinfo);

int adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		 const double *Z, const int n, int *misst);

int list_dups (const int *list, int ci);

double getvalue (const char *s, const double *Z, const DATAINFO *pdinfo);

int lagvarnum (const int iv, const int lag, const DATAINFO *pdinfo);

int forecast (int t1, const int t2, const int nv, 
	      const MODEL *pmod, DATAINFO *pdinfo, double **pZ);

int full_model_list (MODEL *pmod, int **plist);

void init_model (MODEL *pmod);

double tcrit95 (const int df);

int allocate_case_markers (char ***S, int n);

int _ztox (const int nlv, double *px, const DATAINFO *pdinfo, 
	   const double *Z);

int _ztoxy (const int v1, const int v2, double *px, double *py, 
	    const DATAINFO *pdinfo, const double *Z);

void _pgbreak (const int n, int *lineno, const int batch);

int get_precision (double *x, int n);

#endif /* GRETL_INTERNAL_H */
