/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* The beginings of a GARCH plugin using the Fiorentini, Calzolari and
   Panattoni (fcp) fortran code.
*/

#include "libgretl.h"
#include "internal.h"

#include "fcp.h"

/* vsanal error codes (in "info"):
   0 OK
   1 Insufficient dimensions in FORTRAN arrays
   2 No convergence due to bad gradient
   3 ML iteration did not converge
*/

int do_fcp (const int *list, const double **Z, 
	    const DATAINFO *pdinfo, PRN *prn)
{
    int t1, t2;
    double *yobs, *xobs;
    int nexo, nobs;
    double *ydet;
    int ncoeff, ncoefb;
    double *coeff, *d;
    double *vc, *res, *res2;
    double *ystoc;
    double *amax;
    double *b;
    double oldc, yy, umc, sigma;
    int info, iters = 0;
    
    int i, p, q, ynum;

    /* FIXME: how exactly should these be set? */
    nexo = 0;
    ncoeff = 1;
    ncoefb = 1;

    t2 = nobs = pdinfo->t2 - pdinfo->t1 + 1;
    t1 = 1;

    p = list[1];
    q = list[2];
    ynum = list[4];

    yobs = malloc(nobs * sizeof *yobs);
    ydet = malloc(nobs * sizeof *ydet);
    ystoc = malloc(nobs * sizeof *ystoc);

    res2 = malloc(nobs * sizeof *res2);
    for (i=0; i<nobs; i++) {
	res2[i] = 0.0;
    }

    res = malloc(nobs * sizeof *res);
    for (i=0; i<nobs; i++) {
	res[i] = 0.0;
    }    

    sigma = umc = oldc = yy = 0.0;

    amax = malloc(nobs * sizeof *amax);
    for (i=0; i<nobs; i++) {
	amax[i] = 0.0;
    }

    coeff = malloc(ncoeff * sizeof *coeff);
    b = malloc(ncoeff * sizeof *b);
    for (i=0; i<ncoeff; i++) {
	coeff[i] = b[i] = 0.0;
    }    

    d = malloc((ncoeff + 1) * sizeof *d);
    vc = malloc((ncoeff * ncoeff + 1) * sizeof *vc);
    for (i=1; i<=ncoeff; i++) {
	d[i] = vc[i] = 0.0;
    } 

    if (nexo > 0) {
	xobs = malloc((nobs * nexo + 1) * sizeof *xobs);
	/* now fill in exog var values */
    } else {
	xobs = malloc((nobs + 1) * sizeof *xobs);
    }

    for (i=0; i<nobs; i++) {
	ystoc[i] = ydet[i] = yobs[i] = Z[ynum][i];
    }

    /* initialize at unconditional mean of y */
    coeff[0] = _esl_mean(0, nobs - 1, Z[ynum]);

    /* initialize elements of alpha, beta such that 
       alpha_0/(1 - alpha_1 - beta_1) = unconditional
       variance of y (????)
    */
    amax[0] = _esl_variance(0, nobs - 1, Z[ynum]);
    amax[1] = p;
    amax[2] = q; 
    for (i=0; i<p+q; i++) {
	/* initial alpha, beta values */
	amax[3+i] = 0.1;
    }

    vsanal_(t1, t2, yobs, nobs,
	    &xobs[1], nexo, &umc, ydet, 
	    &yy, coeff, ncoeff, &d[1], &oldc, 
	    &vc[1], res2, res, &sigma, ystoc, 
	    amax, b, &ncoefb, &iters, &info, prn);

    if (info != 0) {
	fprintf(stderr, "vsanal returned with info = %d\n", info);
    }

    pprintf(prn, "Number of iterations = %d\n", iters);

    if (info == 0) {
	int k, nparam = ncoeff + p + q + 1;

	pprintf(prn, "Convergence reached, with tolerance = %g\n", 
	       amax[0]);
	pputs(prn, "\nRegression coefficient estimates:\n");
	for (i=1; i<=1+nexo; i++) {
	    pprintf(prn, "    A[%d]: %#14.6g (%#.6g)\n", i, amax[i],
		   amax[i+ncoeff]);
	}
	k = i;
	pputs(prn, "\nGARCH coefficient estimates:\n");
	for (i=k; i<=nparam; i++) {
	    if (i==k) {
		pputs(prn, "alpha[0]: ");
	    } else if (i-k <= p) {
		pprintf(prn, "alpha[%d]: ", i-k);
	    } else {
		pprintf(prn, " beta[%d]: ", i-k-p);
	    }
	    pprintf(prn, "%#14.6g (%#.6g)\n", amax[i], amax[i+nparam]);
	}
	pputs(prn, "\n(Standard errors in parentheses)\n");
    }

    free(yobs);
    free(xobs);
    free(ydet);
    free(coeff);
    free(d);
    free(vc);
    free(res2);
    free(res);
    free(ystoc);
    free(amax);
    free(b);

    return 0;
}

/* the driver function for the plugin */

MODEL garch_model (int *list, const double **Z, DATAINFO *pdinfo,
		   PRN *prn) 
{
    MODEL model;
    PRN *myprn;

#if 0
    /* run initial OLS */
    model = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
    if (model.errcode) {
        return model;
    } 
#endif 

    myprn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    do_fcp(list, Z, pdinfo, myprn); 

    gretl_print_destroy(myprn);

    model.errcode = 1; /* bodge for now */
    return model;
}


