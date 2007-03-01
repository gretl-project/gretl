/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2004 Allin Cottrell
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

#ifndef NONPARAM_H
#define NONPARAM_H
 
int spearman (const int *list, 
	      const double **Z, const DATAINFO *pdinfo,
	      gretlopt opt, PRN *prn);

double lockes_test (const double *x, int t1, int t2);

int runs_test (int varno, const double **Z, const DATAINFO *pdinfo, 
	       PRN *prn);

int diff_test (const int *list, const double **Z, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int sort_pairs_by_x (gretl_matrix *x, gretl_matrix *y, int **order,
		     char **labels);

gretl_matrix *loess_fit (const gretl_matrix *x, const gretl_matrix *y,
			 int d, double q, gretlopt opt, int *err);

#endif /* NONPARAM_H */

