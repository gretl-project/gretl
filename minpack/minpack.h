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

/* L-BFGS-B (lbfgsb.c) */

/* Callback to compute the function value *f and gradient g at x.
   The callback may alter x temporarily (e.g. for numerical
   differentiation) but must restore it exactly. Return non-zero to
   stop the minimization.
*/
typedef int (*LBFGSB_FG) (double *x, double *f, double *g, void *data);

/* Callback invoked at the starting point (iter = 0) and after each
   iteration, with the number of iterations completed, the number of
   function evaluations so far, the current function value, iterate
   and gradient, and the step length in the last line search. Return
   non-zero to stop the minimization.
*/
typedef int (*LBFGSB_ITER) (int iter, int nfg, double f,
			    const double *x, const double *g,
			    double stp, void *data);

enum {
    LBFGSB_CONVERGED,
    LBFGSB_ABNORMAL,
    LBFGSB_STOPPED,
    LBFGSB_BADINPUT,
    LBFGSB_NOMEM,
    LBFGSB_NONFINITE
};

int lbfgsb (int n, int m, double *x, const double *l, const double *u,
	    const int *nbd, double factr, double pgtol, double *f,
	    double *g, LBFGSB_FG fgfunc, LBFGSB_ITER iterfunc,
	    void *data, int *iter, int *nfg, char *msg);

#endif /* GRETL_MINPACK_H */
