#ifndef FCP_H
#define FCP_H

int vsanal_(int ninit, int nfinsm, double *yobs,
	    int iread, const double **X, int nexo, 
	    double *ydet, double *yy, double *coeff, int ncoeff,
	    double *oldc, double *vc, double *res2,
	    double *res, double *ystoc,
	    double *amax, double *b, int *iters, 
	    PRN *prn);

#endif /* FCP_H */




