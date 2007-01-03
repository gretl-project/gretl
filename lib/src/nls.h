/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_NLS_H
#define GRETL_NLS_H

typedef struct _nlspec nlspec;

typedef double (*BFGS_CRIT_FUNC) (const double *, void *);
typedef int (*BFGS_GRAD_FUNC) (double *, double *, int, 
			       BFGS_CRIT_FUNC, void *);
typedef double *(*BFGS_SCORE_FUNC) (const double *, int, void *);


nlspec *nlspec_new (int ci, const DATAINFO *pdinfo);

void nlspec_destroy (nlspec *spec);

int 
nlspec_add_param_with_deriv (nlspec *spec, 
			     const char *dstr,
			     const double **Z, 
			     const DATAINFO *pdinfo);

int nlspec_add_param_list (nlspec *spec, const int *list,
			   const double **Z, const DATAINFO *pdinfo);

int 
nlspec_set_regression_function (nlspec *spec, 
				const char *fnstr, 
				const DATAINFO *pdinfo);

void nlspec_set_t1_t2 (nlspec *spec, int t1, int t2);

int nls_parse_line (int ci, const char *line, const double **Z,
		    const DATAINFO *pdinfo, PRN *prn);

MODEL nls (double ***pZ, DATAINFO *pdinfo, gretlopt opt, PRN *prn);

MODEL model_from_nlspec (nlspec *spec, double ***pZ, 
			 DATAINFO *pdinfo, gretlopt opt, 
			 PRN *prn);

double get_default_nls_toler (void);

int BFGS_max (double *b, int n, int maxit, double reltol,
	      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	      BFGS_GRAD_FUNC gradfunc, void *data, 
	      gretlopt opt, PRN *prn);

int BFGS_numeric_gradient (double *b, double *g, int n,
			   BFGS_CRIT_FUNC func, void *data);

gretl_matrix *build_OPG_matrix (double *b, int k, int T,
				BFGS_SCORE_FUNC scorefun,
				void *data, int *err);

double *numerical_hessian (double *b, int n, BFGS_CRIT_FUNC func, 
			   void *data);

double user_BFGS (gretl_matrix *b, const char *fncall,
		  double ***pZ, DATAINFO *pdinfo,
		  PRN *prn, int *err);

gretl_matrix *fdjac (gretl_matrix *theta, const char *fncall,
		     double ***pZ, DATAINFO *pdinfo,
		     int *err);

#endif /* GRETL_NLS_H */
