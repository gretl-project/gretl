#ifndef FCP_H
#define FCP_H

int garch_estimate (const double *y, const double **X, 
		    int t1, int t2, int nobs, int nc,
		    int p, int q, double *theta, gretl_matrix *V, 
		    double *e, double *e2, double *h,
		    double scale, double *pll, int *iters, 
		    int vopt, PRN *prn);

#endif /* FCP_H */




