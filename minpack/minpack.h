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

extern int lmdif_(S_fp fcn, integer *m, integer *n, doublereal *x, 
		  doublereal *fvec, doublereal *ftol, doublereal *xtol, 
		  doublereal *gtol, integer *maxfev, doublereal *epsfcn, 
		  doublereal *diag, integer *mode, doublereal *factor, 
		  integer *nprint, integer *info, integer *nfev, 
		  doublereal *fjac, integer *ldfjac, integer *ipvt, 
		  doublereal *qtf, doublereal *wa1, doublereal *wa2, 
		  doublereal *wa3, doublereal *wa4);

extern int fdjac2_(S_fp fcn, integer *m, integer *n, doublereal *x, 
		   doublereal *fvec, doublereal *fjac, integer *ldfjac, 
		   integer *iflag, doublereal *epsfcn, doublereal *wa);


extern doublereal dpmpar_(integer *i__);


#endif

