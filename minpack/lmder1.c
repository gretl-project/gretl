#include "minpack.h"

/*
c     lmder1:
c
c     simplified driver for lmder, which sets various control
c     parameters to defalt values; see lmder for details
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int lmder1_(S_fp fcn, int m, int n, double *x, 
	    double *fvec, double *fjac, int ldfjac, double *tol, 
	    int *info, int *ipvt, double *wa, int lwa, void *p)
{
    /* default control values */
    const double factor = 100;
    int maxfev = (n + 1) * 100;
    double ftol = *tol;
    double xtol = *tol;
    double gtol = 0.0;
    int nprint = 0;
    int mode = 1;
    /* iteration counts */
    int nfev, njev;

    *info = 0;

    /* check the input parameters for errors */
    if (n <= 0 || m < n || ldfjac < m || 
	*tol < 0.0 || lwa < n * 5 + m) {
	return 0;
    }

    lmder_((S_fp) fcn, m, n, x, fvec, fjac, ldfjac, 
	   ftol, xtol, gtol, maxfev, wa, mode, factor, nprint, 
	   info, &nfev, &njev, ipvt, wa + n, wa + 2*n, 
	   wa + 3*n, wa + 4*n, wa + 5*n, p);

    if (*info == 8) {
	*info = 4;
    }

    return 0;
}
