/* 
   This source file based on Minpack: initially converted from 
   fortran using f2c, then rendered into relatively idiomatic
   C with zero-based indexing throughout and pass-by-value for
   parameters that do not function as pointers. We also rely
   on <float.h> for the machine precision rather than Minpack's
   dpmpar().

   See README in this directory for the Minpack Copyright.

   Allin Cottrell, Wake Forest University, April 2012
*/

#include "minpack.h"
#include <math.h>
#include <float.h>

/*
c     lmder:
c
c     the purpose of lmder is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions and the jacobian.
c
c     the subroutine statement is
c
c       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
c                        maxfev,diag,mode,factor,nprint,info,nfev,
c                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions and the jacobian. fcn must
c         be declared in an external statement in the user
c         calling program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
c         integer m,n,ldfjac,iflag
c         double precision x(n),fvec(m),fjac(ldfjac,n)
c         ----------
c         if iflag = 1 calculate the functions at x and
c         return this vector in fvec. do not alter fjac.
c         if iflag = 2 calculate the jacobian at x and
c         return this matrix in fjac. do not alter fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmder.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn with iflag = 1
c         has reached maxfev.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.).100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x, fvec, and fjac
c         available for printing. fvec and fjac should not be
c         altered. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn with iflag = 1 has
c                   reached maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn with iflag = 1.
c
c       njev is an integer output variable set to the number of
c         calls to fcn with iflag = 2.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... enorm,lmpar,qrfac
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
*/

int lmder_(S_fp fcn, int m, int n, double *x, double *fvec, 
	   double *fjac, int ldfjac, double ftol, double xtol, 
	   double gtol, int maxfev, double *diag, int mode, 
	   double factor, int nprint, int *info, int *nfev, 
	   int *njev, int *ipvt, double *qtf, double *wa1, 
	   double *wa2, double *wa3, double *wa4, void *p)
{
    const double p1 = .1;
    const double p5 = .5;
    const double p25 = .25;
    const double p75 = .75;
    const double p0001 = 1e-4;
    const double epsmch = DBL_EPSILON;
    double d, par, sum;
    double temp1, temp2;
    double ratio, delta = 0;
    double fnorm, gnorm, pnorm, fnorm1;
    double actred, dirder, prered;
    double xnorm = 0.0, temp = 0.0;
    int iter, iflag = 0;
    int i, j, l;

    *info = 0;
    *nfev = 0;
    *njev = 0;

    /* check the input parameters for errors */

    if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0 || 
	    gtol < 0.0 || maxfev <= 0 || factor <= 0.0) {
	goto bailout;
    }

    if (mode == 2) {
	for (j = 0; j < n; ++j) {
	    if (diag[j] <= 0.0) {
		goto bailout;
	    }
	}
    }

    /* evaluate the function at the starting point
       and calculate its norm */

    iflag = 1;
    (*fcn)(m, n, x, fvec, fjac, ldfjac, &iflag, p);
    *nfev = 1;
    if (iflag < 0) {
	goto bailout;
    }
    fnorm = enorm_(m, fvec);

    /* initialize Levenberg-Marquardt parameter and iter counter */
    par = 0.0;
    iter = 1;

    /* beginning of the outer loop */
 outer_start:

    /* calculate the jacobian matrix */

    iflag = 2;
    (*fcn)(m, n, x, fvec, fjac, ldfjac, &iflag, p);
    *njev += 1;
    if (iflag < 0) {
	goto bailout;
    }

    /* if requested, call fcn to enable printing of iterates */
    if (nprint > 0) {
	iflag = 0;
	if ((iter - 1) % nprint == 0) {
	    (*fcn)(m, n, x, fvec, fjac, ldfjac, &iflag, p);
	}
	if (iflag < 0) {
	    goto bailout;
	}
    }

    /* compute the QR factorization of the jacobian */
    qrfac_(m, n, fjac, ldfjac, ipvt, wa1, wa2, wa3);

    /* on the first iteration and if mode is 1, scale according
       to the norms of the columns of the initial jacobian 
    */

    if (iter == 1) {
	if (mode != 2) {
	    for (j = 0; j < n; ++j) {
		diag[j] = wa2[j];
		if (wa2[j] == 0.0) {
		    diag[j] = 1.0;
		}
	    }
	}
	/* on the first iteration, calculate the norm of the scaled x
	   and initialize the step bound delta 
        */
	for (j = 0; j < n; ++j) {
	    wa3[j] = diag[j] * x[j];
	}
	xnorm = enorm_(n, wa3);
	delta = factor * xnorm;
	if (delta == 0.0) {
	    delta = factor;
	}
    }

    /* form Q'*fvec and store the first n components in qtf */

    for (i = 0; i < m; ++i) {
	wa4[i] = fvec[i];
    }
    for (j = 0; j < n; ++j) {
	if (fjac[j + j * ldfjac] != 0.0) {
	    sum = 0.0;
	    for (i = j; i < m; ++i) {
		sum += fjac[i + j * ldfjac] * wa4[i];
	    }
	    temp = -sum / fjac[j + j * ldfjac];
	    for (i = j; i < m; ++i) {
		wa4[i] += fjac[i + j * ldfjac] * temp;
	    }
	}
	fjac[j + j * ldfjac] = wa1[j];
	qtf[j] = wa4[j];
    }

    /* compute the norm of the scaled gradient */

    gnorm = 0.0;
    if (fnorm != 0.0) {
	for (j = 0; j < n; ++j) {
	    l = ipvt[j];
	    if (wa2[l] != 0.0) {
		sum = 0.0;
		for (i = 0; i <= j; ++i) {
		    sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
		}
		d = fabs(sum / wa2[l]);
		gnorm = max(gnorm, d);
	    }
	}
    }

    /* test for convergence of the gradient norm */
    if (gnorm <= gtol) {
	*info = 4;
    }
    if (*info != 0) {
	goto bailout;
    }

    /* rescale if necessary */
    if (mode != 2) {
	for (j = 0; j < n; ++j) {
	    diag[j] = max(diag[j], wa2[j]);
	}
    }

    /* beginning of the inner loop */

    while (1) {

	/* determine the Levenberg-Marquardt parameter */
	lmpar_(n, fjac, ldfjac, ipvt, diag, qtf, delta, 
	       &par, wa1, wa2, wa3, wa4);

	/* store the direction p and x + p; calculate the norm of p */
	for (j = 0; j < n; ++j) {
	    wa1[j] = -wa1[j];
	    wa2[j] = x[j] + wa1[j];
	    wa3[j] = diag[j] * wa1[j];
	}
	pnorm = enorm_(n, wa3);

	/* on the first iteration, adjust the initial step bound */
	if (iter == 1) {
	    delta = min(delta, pnorm);
	}

	/* evaluate the function at x + p and calculate its norm */
	iflag = 1;
	(*fcn)(m, n, wa2, wa4, fjac, ldfjac, &iflag, p);
	*nfev += 1;
	if (iflag < 0) {
	    goto bailout;
	}
	fnorm1 = enorm_(m, wa4);

	/* compute the scaled actual reduction */
	actred = -1.0;
	if (p1 * fnorm1 < fnorm) {
	    d = fnorm1 / fnorm;
	    actred = 1.0 - d * d;
	}

	/* compute the scaled predicted reduction and
	   the scaled directional derivative 
	*/
	for (j = 0; j < n; ++j) {
	    wa3[j] = 0.0;
	    l = ipvt[j];
	    temp = wa1[l];
	    for (i = 0; i <= j; ++i) {
		wa3[i] += fjac[i + j * ldfjac] * temp;
	    }
	}
	temp1 = enorm_(n, wa3) / fnorm;
	temp2 = sqrt(par) * pnorm / fnorm;
	prered = temp1 * temp1 + temp2 * temp2 / p5;
	dirder = -(temp1 * temp1 + temp2 * temp2);

	/* compute the ratio of the actual to the predicted
	   reduction */

	ratio = 0.0;
	if (prered != 0.0) {
	    ratio = actred / prered;
	}

	/* update the step bound */
	if (ratio <= p25) {
	    if (actred >= 0.0) {
		temp = p5;
	    }
	    if (actred < 0.0) {
		temp = p5 * dirder / (dirder + p5 * actred);
	    }
	    if (p1 * fnorm1 >= fnorm || temp < p1) {
		temp = p1;
	    }
	    d = pnorm / p1;
	    delta = temp * min(delta, d);
	    par /= temp;
	} else if (par == 0.0 || ratio >= p75) {
	    delta = pnorm / p5;
	    par = p5 * par;
	}

	/* test for successful iteration */

	if (ratio >= p0001) {
	    /* successful iteration: update x, fvec, and their norms */
	    for (j = 0; j < n; ++j) {
		x[j] = wa2[j];
		wa2[j] = diag[j] * x[j];
	    }
	    for (i = 0; i < m; ++i) {
		fvec[i] = wa4[i];
	    }
	    xnorm = enorm_(n, wa2);
	    fnorm = fnorm1;
	    ++iter;
	}

	/* tests for convergence */

	if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0) {
	    *info = 1;
	}
	if (delta <= xtol * xnorm) {
	    *info = 2;
	}
	if (fabs(actred) <= ftol && prered <= ftol && 
	    p5 * ratio <= 1.0 && *info == 2) {
	    *info = 3;
	}
	if (*info != 0) {
	    goto bailout;
	}

	/* tests for termination and stringent tolerances */

	if (*nfev >= maxfev) {
	    *info = 5;
	}
	if (fabs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1.0) {
	    *info = 6;
	}
	if (delta <= epsmch * xnorm) {
	    *info = 7;
	}
	if (gnorm <= epsmch) {
	    *info = 8;
	}
	if (*info != 0) {
	    goto bailout;
	}

	/* end of inner loop: repeat (only) if iteration unsuccessful */
	if (ratio >= p0001) {
	    break;
	}
    }

    /* end of the outer loop */
    goto outer_start;

 bailout:

    /* termination, either normal or user-imposed */

    if (iflag < 0) {
	*info = iflag;
    }
    iflag = 0;
    if (nprint > 0) {
	(*fcn)(m, n, x, fvec, fjac, ldfjac, &iflag, p);
    }

    return 0;
}

