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

#define _grow_Z(n, p1, p2) dataset_add_vars((n), p1, p2)
#define _shrink_Z(n, p1, p2) dataset_drop_vars((n), p1, p2)

int _laggenr (const int iv, const int lag, const int opt, double **pZ, 
	      DATAINFO *pdinfo);

int _multiply (char *s, int *list, char *sfx, double **pZ,
	       DATAINFO *pdinfo);

void _print_add (const COMPARE *add, const int *addvars, 
		 const DATAINFO *pdinfo, const int aux_code, PRN *prn);

void _print_omit (const COMPARE *omit, const int *omitvars, 
		  const DATAINFO *pdinfo, PRN *prn);

void _graphyzx (const int *list, const double *zy1, const double *zy2, 
		const double *zx, int n, const char *yname, 
		const char *xname, const DATAINFO *pdinfo, 
		int oflag, PRN *prn);

void _printxs (double xx, int n, int ci, PRN *prn);

void _bufspace (int n, PRN *prn);

void _print_ar (MODEL *pmod, PRN *prn);

void _delete (char *str, const int indx, const int count);

int _isnumber (const char *str);

void _esl_trunc (char *str, const int n);

int _count_fields (const char *str);

void _shiftleft (char *str, size_t move);

double _corr (const int n, const double *zx, const double *zy);

double _covar (const int n, const double *zx, const double *zy);

int _iszero (const int t1, const int t2, const double *x);

int _isconst (const int t1, const int t2, const double *x);

double _esl_mean (const int t1, const int t2, const double *x);

void _minmax (const int t1, const int t2, const double zx[], 
	      double *min, double *max);

int _hasconst (const int *list);

int _compare_doubles (const void *a, const void *b);

double _esl_stddev (const int t1, const int t2, const double *x);

double _esl_variance (const int t1, const int t2, const double *x);

void _aicetc (MODEL *pmod);

void _criteria (const double ess, const int nobs, const int ncoeff, 
		PRN *prn);

int _adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		  const double *Z, const int n, int *misst);

int _list_dups (const int *list, int ci);

int _lagvarnum (const int iv, const int lag, const DATAINFO *pdinfo);

int _forecast (int t1, const int t2, const int nv, 
	       const MODEL *pmod, DATAINFO *pdinfo, double **pZ);

int _full_model_list (MODEL *pmod, int **plist);

void _init_model (MODEL *pmod);

double _tcrit95 (const int df);

int _allocate_case_markers (char ***S, int n);

int _ztoxy (const int v1, const int v2, double *px, double *py, 
	    const DATAINFO *pdinfo, const double *Z);

int _get_precision (double *x, int n);

#endif /* GRETL_INTERNAL_H */
