#ifndef GRETL_MINPACK_H
#define GRETL_MINPACK_H

extern int chkder_(integer *m, integer *n, doublereal *x,
		   doublereal *fvec, doublereal *fjac, integer *ldfjac,
		   doublereal *xp, doublereal *fvecp, integer *mode,
		   doublereal *err);

extern int lmder1_(U_fp fcn, integer *m, integer *n, 
		   doublereal *x, doublereal *fvec, doublereal *fjac, 
		   integer *ldfjac, doublereal *tol, integer *info, 
		   integer *ipvt, doublereal *wa, integer *lwa);

#endif

