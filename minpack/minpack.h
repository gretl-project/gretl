#ifndef GRETL_MINPACK_H
#define GRETL_MINPACK_H

typedef int (*S_fp6)(int, int, double *, double *, int *, void *);
typedef int (*S_fp8)(int, int, double *, double *, double *,
                     int, int *, void *);

#ifndef min
# define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
# define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

double enorm_(int n, double *x);

int fdjac2_(S_fp6 fcn, int m, int n, int quality, double *x,
	    double *fvec, double *fjac, int ldfjac, int *iflag,
	    double epsfcn, double *wa, void *p);

int lmder_(S_fp8 fcn, int m, int n, double *x, double *fvec,
	   double *fjac, int ldfjac, double ftol, double xtol,
	   double gtol, int maxfev, double *diag, int mode,
	   double factor, int nprint, int *info, int *nfev,
	   int *njev, int *ipvt, double *qtf, double *wa1,
	   double *wa2, double *wa3, double *wa4, void *p);

int lmdif_(S_fp6 fcn, int m, int n, double *x, double *fvec,
	   double ftol, double xtol, double gtol, int maxfev,
	   double epsfcn, double *diag, int mode, double factor,
	   int nprint, int *info, int *nfev, double *fjac,
	   int ldfjac, int *ipvt, double *qtf, double *wa1,
	   double *wa2, double *wa3, double *wa4, void *p);

int lmpar_(int n, double *r, int ldr, int *ipvt, double *diag,
	   double *qtb, double delta, double *par, double *x,
	   double *sdiag, double *wa1, double *wa2);

int qrfac_(int m, int n, double *a, int lda, int *ipvt,
	   double *rdiag, double *acnorm, double *wa);

int qrsolv_(int n, double *r, int ldr, int *ipvt,
	    double *diag, double *qtb, double *x,
	    double *sdiag, double *wa);

int setulb_(int *n, int *m, double *x, double *l, double *u,
	    int *nbd, double *f, double *g, double *factr,
	    double *pgtol, double *wa, int *iwa,
	    char *task, char *csave, int *lsave,
	    int *isave, double *dsave);

#endif /* GRETL_MINPACK_H */
