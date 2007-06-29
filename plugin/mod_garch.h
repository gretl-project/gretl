#ifndef MOD_GARCH_H
#define MOD_GARCH_H

int garch_estimate_mod (const double *y, const double **X,
			int t1, int t2, int nobs, int nc,
			int p, int q, double *theta,  gretl_matrix *V,
			double *e, double *e2, double *h,
			double scale, double *pll, 
			int *fncount, int *grcount,
			int vopt, PRN *prn);

#endif /* MOD_GARCH_H */
