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
    SIMANN_MAX,
    NM_MAX,
    GSS_MAX,
    ROOT_FIND
} MaxMethod;

typedef double (*BFGS_CRIT_FUNC) (const double *, void *);
typedef int (*BFGS_GRAD_FUNC) (double *, double *, int, 
			       BFGS_CRIT_FUNC, void *);
typedef double (*BFGS_COMBO_FUNC) (double *, double *, int, void *);
typedef const double *(*BFGS_LLT_FUNC) (const double *, int, void *);
typedef int (*HESS_FUNC) (double *, gretl_matrix *, void *);
typedef double (*ZFUNC) (double, void *);

int BFGS_max (double *b, int n, int maxit, double reltol,
	      int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	      int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	      const gretl_matrix *A0, gretlopt opt, PRN *prn);

int LBFGS_max (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc,
	       BFGS_COMBO_FUNC combo_func, void *data,
	       const gretl_matrix *bounds, gretlopt opt, PRN *prn);

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

int numerical_hessian (double *b, gretl_matrix *H,
		       BFGS_CRIT_FUNC func, void *data,
		       int neg, double d);

gretl_matrix *numerical_hessian_inverse (const double *b, int n, 
					 BFGS_CRIT_FUNC func, 
					 void *data, double d,
					 int *err);

double user_BFGS (gretl_matrix *b, 
		  const char *fncall,
		  const char *gradcall,
		  DATASET *dset,
		  const gretl_matrix *bounds,
		  int minimize, PRN *prn, 
		  int *err);

double user_NR (gretl_matrix *b, 
		const char *fncall,
		const char *gradcall, 
		const char *hesscall,
		DATASET *dset,
		int minimize, PRN *prn, 
		int *err);

double deriv_free_optimize (MaxMethod method,
			    gretl_matrix *b,
			    const char *fncall,
			    int maxit,
			    double tol,
			    int minimize,
			    DATASET *dset,
			    PRN *prn,
			    int *err);

gretl_matrix *user_fdjac (gretl_matrix *theta, const char *fncall,
			  double eps, DATASET *dset, int *err);

gretl_matrix *user_numhess (gretl_matrix *b, const char *fncall,
			    double d, DATASET *dset, int *err);

int gretl_simann (double *theta, int n, int maxit,
		  BFGS_CRIT_FUNC cfunc, void *data,
		  gretlopt opt, PRN *prn);

int gretl_amoeba (double *theta, int n, int maxit,
		  BFGS_CRIT_FUNC cfunc, void *data,
		  gretlopt opt, PRN *prn);

int gretl_fzero (double *bracket, double tol,
		 ZFUNC zfunc, void *data, double *px,
		 gretlopt opt, PRN *prn);

int gretl_gss (double *theta, double tol, int *ic,
	       BFGS_CRIT_FUNC cfunc, void *data,
	       gretlopt opt, PRN *prn);

void BFGS_defaults (int *maxit, double *tol, int ci);

int optimizer_get_matrix_name (const char *fncall, char *name);

int numgrad_in_progress (void);

#endif /* GRETL_BFGS_H */
