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

#ifndef GRETL_BFGS_H
#define GRETL_BFGS_H

#define RSTEPS 4

typedef double (*BFGS_CRIT_FUNC) (const double *, void *);
typedef int (*BFGS_GRAD_FUNC) (double *, double *, int, 
			       BFGS_CRIT_FUNC, void *);
typedef double *(*BFGS_SCORE_FUNC) (const double *, int, void *);

int BFGS_max (double *b, int n, int maxit, double reltol,
	      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	      gretlopt opt, PRN *prn);

int LBFGS_max (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn);

int BFGS_numeric_gradient (double *b, double *g, int n,
			   BFGS_CRIT_FUNC func, void *data);

gretl_matrix *build_OPG_matrix (double *b, int k, int T,
				BFGS_SCORE_FUNC scorefun,
				void *data, int *err);

double *numerical_hessian (double *b, int n, BFGS_CRIT_FUNC func, 
			   void *data, int *err);

double user_BFGS (gretl_matrix *b, const char *fncall,
		  double ***pZ, DATAINFO *pdinfo,
		  PRN *prn, int *err);

gretl_matrix *fdjac (gretl_matrix *theta, const char *fncall,
		     double ***pZ, DATAINFO *pdinfo,
		     int *err);

#endif /* GRETL_BFGS_H */
