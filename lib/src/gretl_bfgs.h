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

typedef enum {
    BHHH_MAX,
    BFGS_MAX,
    LBFGS_MAX,
    LM_MAX
} OptimizerCode;

typedef double (*BFGS_CRIT_FUNC) (const double *, void *);
typedef int (*BFGS_GRAD_FUNC) (double *, double *, int, 
			       BFGS_CRIT_FUNC, void *);
typedef const double *(*BFGS_LLT_FUNC) (const double *, int, void *);
typedef int (*HESS_FUNC) (double *, gretl_matrix *, void *);

int BFGS_max (double *b, int n, int maxit, double reltol,
	      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	      gretl_matrix *A0, gretlopt opt, PRN *prn);

int LBFGS_max (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn);

int newton_raphson_max (double *b, int n, int maxit, 
			double crittol, double gradtol, 
			int *itercount, int crittype, 
			BFGS_CRIT_FUNC cfunc,
			BFGS_GRAD_FUNC gradfunc, 
			HESS_FUNC hessfunc,
			void *data, gretlopt opt, PRN *prn);

int BFGS_numeric_gradient (double *b, double *g, int n,
			   BFGS_CRIT_FUNC func, void *data);

gretl_matrix *numerical_score_matrix (double *b, int T, int k,
				      BFGS_LLT_FUNC lltfun,
				      void *data, int *err);

int hessian_from_score (double *b, gretl_matrix *H,
			BFGS_GRAD_FUNC gradfunc, 
			BFGS_CRIT_FUNC cfunc,
			void *data);

gretl_matrix *hessian_inverse_from_score (double *b, int n,
					  BFGS_GRAD_FUNC gradfunc,
					  BFGS_CRIT_FUNC cfunc,
					  void *data, int *err);

gretl_matrix *numerical_hessian_inverse (const double *b, int n, 
					 BFGS_CRIT_FUNC func, 
					 void *data, int *err);

double user_BFGS (gretl_matrix *b, 
		  const char *fncall,
		  const char *gradcall,
		  DATASET *dset, 
		  PRN *prn, 
		  int *err);

double user_NR (gretl_matrix *b, 
		const char *fncall,
		const char *gradcall, 
		const char *hesscall,
		DATASET *dset,
		PRN *prn, 
		int *err);

double user_simann (gretl_matrix *b, 
		    const char *fncall,
		    int maxit, 
		    DATASET *dset,
		    PRN *prn, 
		    int *err);

gretl_matrix *fdjac (gretl_matrix *theta, const char *fncall,
		     DATASET *dset, int *err);

int gretl_simann (double *theta, int n, int maxit,
		  BFGS_CRIT_FUNC cfunc, void *data,
		  gretlopt opt, PRN *prn);

void BFGS_defaults (int *maxit, double *tol, int ci);

int optimizer_get_matrix_name (const char *fncall, char *name);

#endif /* GRETL_BFGS_H */
