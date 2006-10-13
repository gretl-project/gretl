#ifndef MOD_GARCH_H
#define MOD_GARCH_H

int garch_estimate_mod (int t1, int t2, int nobs, 
			const double **X, int nx, double *coeff, int nc, 
			double *vcv, double *res2, double *res, double *h,
			double *y, double *amax, double *b, 
			double scale, int *iters, PRN *prn, int vopt);

#endif /* MOD_GARCH_H */
