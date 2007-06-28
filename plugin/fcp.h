#ifndef FCP_H
#define FCP_H

int garch_estimate (int t1, int t2, int nobs, 
		    const double **X, int nx, double *coeff, int nc, 
		    gretl_matrix *V, double *res2, double *res, double *h,
		    const double *y, double *amax, double *b, 
		    double scale, int *iters, PRN *prn, int vopt);

#endif /* FCP_H */




