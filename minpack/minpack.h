#ifndef LMDER1_H
#define LMDER1_H

extern int lmder1_(U_fp fcn, integer *m, integer *n, 
		   doublereal *x, doublereal *fvec, doublereal *fjac, 
		   integer *ldfjac, doublereal *tol, integer *info, 
		   integer *ipvt, doublereal *wa, integer *lwa);

#endif

