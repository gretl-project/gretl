#ifndef FCP_H
#define FCP_H

int vsanal_(int ninit, int nfinsm, double *yobs,
	    int iread, double *xobs, int nexo, double *umc,
	    double *ydet, double *yy, double *coeff, int ncoeff,
	    double *d__, double *oldc, double *vc, double *res2,
	    double *res, double *sigma, double *ystoc,
	    double *amax, double *b, int *ncoefb,
	    int *iters, int *info, PRN *prn);

#endif /* FCP_H */




