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

int _laggenr (int iv, int lag, int opt, double ***pZ, 
	      DATAINFO *pdinfo);

int _multiply (char *s, int *list, char *sfx, double ***pZ,
	       DATAINFO *pdinfo);

void gretl_print_add (const COMPARE *add, const int *addvars, 
		      const DATAINFO *pdinfo, int aux_code, PRN *prn);

void gretl_print_omit (const COMPARE *omit, const int *omitvars, 
		       const DATAINFO *pdinfo, PRN *prn);

void _graphyzx (const int *list, const double *zy1, const double *zy2, 
		const double *zx, int n, const char *yname, 
		const char *xname, const DATAINFO *pdinfo, 
		int oflag, PRN *prn);

void _printxs (double xx, int n, int ci, PRN *prn);

void _bufspace (int n, PRN *prn);

void _print_ar (MODEL *pmod, PRN *prn);

void _delete (char *str, int indx, int count);

int _isnumber (const char *str);

void _esl_trunc (char *str, int n);

int _count_fields (const char *str);

void _shiftleft (char *str, size_t move);

double _corr (int n, const double *zx, const double *zy);

double _covar (int n, const double *zx, const double *zy);

int _iszero (int t1, int t2, const double *x);

int _isconst (int t1, int t2, const double *x);

double _esl_mean (int t1, int t2, const double *x);

void _minmax (int t1, int t2, const double zx[], 
	      double *min, double *max);

int _hasconst (const int *list);

int _compare_doubles (const void *a, const void *b);

double _esl_stddev (int t1, int t2, const double *x);

double _esl_variance (int t1, int t2, const double *x);

void _aicetc (MODEL *pmod);

void _criteria (const double ess, int nobs, int ncoeff, 
		PRN *prn);

int _adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		  double **Z, int *misst);

int _list_dups (const int *list, int ci);

int _lagvarnum (int iv, int lag, const DATAINFO *pdinfo);

int _forecast (int t1, int t2, int nv, 
	       const MODEL *pmod, double ***pZ);

int _full_model_list (MODEL *pmod, int **plist);

void _init_model (MODEL *pmod, const DATAINFO *pdinfo);

double _tcrit95 (int df);

int _ztoxy (int v1, int v2, double *px, double *py, 
	    const DATAINFO *pdinfo, double **Z);

int _reserved (const char *str);

void gretl_test_init (GRETLTEST *test);

#endif /* GRETL_INTERNAL_H */
