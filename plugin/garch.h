#ifndef GARCH_H
#define GARCH_H

int garch_estimate (const double *y, const double **X, 
		    int t1, int t2, int nobs, int nc,
		    int p, int q, double *theta, gretl_matrix *V, 
		    double *e, double *e2, double *h,
		    double scale, double *pll, int *iters, 
		    int vopt, PRN *prn);

gretl_matrix *
garch_analytical_hessian (const double *y, const double **X, 
			  int t1, int t2, int nobs, int nc,
			  int p, int q, double *theta, 
			  double *e, double *e2, double *h,
			  double scale, int *err);

int garch_estimate_mod (const double *y, const double **X,
			int t1, int t2, int nobs, int nc, 
			int p, int q, double *theta,  gretl_matrix *V,
			double *e, double *e2, double *h,
			double scale, double *pll, 
			int *fncount, int *grcount,
			int vopt, PRN *prn);

#endif /* GARCH_H */




