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

#include <f2c.h>
#include "fcp.P"

/* vsanal error codes (in "info"):
   0 OK
   1 Insufficient dimensions in FORTRAN arrays
   2 No convergence due to bad gradient
   3 ML iteration did not converge
*/

int do_fcp (const int *list, const double **Z, 
	    const DATAINFO *pdinfo, PRN *prn)
{
    integer t1, t2;
    doublereal *yobs, *xobs;
    integer nexo, nobs;
    doublereal *umc, *ydet, *yy;
    integer ncoeff, ncoefb;
    doublereal *coeff, *d, *oldc;
    doublereal *vc, *res2;
    doublereal *res, *sigma;
    doublereal *ystoc;
    doublereal *a, *amax, *amin;
    doublereal *b;
    integer iters, info;
    
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

    yobs = malloc((nobs + 1) * sizeof *yobs);
    ydet = malloc((nobs + 1) * sizeof *ydet);
    ystoc = malloc((nobs + 1) * sizeof *ystoc);

    res2 = malloc((nobs + 1) * sizeof *res2);
    for (i=1; i<=nobs; i++) {
	res2[i] = 0.0;
    }

    umc = malloc(2 * sizeof *umc);
    umc[1] = 0.0;

    res = malloc((nobs + 1) * sizeof *res);
    for (i=1; i<=nobs; i++) {
	res[i] = 0.0;
    }    

    sigma = malloc(2 * sizeof *sigma);
    sigma[1] = 0.0;

    oldc = malloc(2 * sizeof *oldc);
    yy = malloc(2 * sizeof *yy);
    oldc[1] = yy[1] = 0.0;

    a = malloc(2 * sizeof *a);
    a[1] = 0.0;

    amax = malloc((nobs + 1) * sizeof *amax);
    amin = malloc((nobs + 1) * sizeof *amin);
    for (i=1; i<=nobs; i++) {
	amax[i] = amin[i] = 0.0;
    }

    coeff = malloc((ncoeff + 1) * sizeof *coeff);
    b = malloc((ncoeff + 1) * sizeof *b);
    for (i=1; i<=ncoeff; i++) {
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
	xobs = NULL;
    }

    for (i=1; i<=nobs; i++) {
	ystoc[i] = ydet[i] = yobs[i] = Z[ynum][i-1];
    }

    /* initialize at unconditional mean of y */
    coeff[1] = _esl_mean(0, nobs - 1, Z[ynum]);

    /* initialize elements of alpha, beta such that 
       alpha_0/(1 - alpha_1 - beta_1) = unconditional
       variance of y (????)
    */
    amax[1] = _esl_variance(0, nobs - 1, Z[ynum]);
    amax[2] = p;
    amax[3] = q; 
    amax[4] = 0.1; /* initial alpha_1 */
    amax[5] = 0.1; /* initial beta_1 */

    printf("t1=%d\n"
	   "t2=%d\n"
	   "nobs=%d\n"
	   "nexo=%d\n"
	   "ncoeff=%d\n",
	   (int) t1, (int) t2, (int) nobs, 
	   (int) nexo, (int) ncoeff);

    vsanal_(&t1, &t2, &yobs[1], &nobs,
	    xobs, &nexo, &umc[1], &ydet[1], 
	    &yy[1], &coeff[1], &ncoeff, &d[1], &oldc[1], 
	    &vc[1], &res2[1], &res[1], &sigma[1], &a[1], &ystoc[1], 
	    &amax[1], &amin[1], &b[1], &ncoefb, &iters, &info, prn);

    if (info != 0) {
	fprintf(stderr, "vsanal returned with info = %d\n", (int) info);
    }

    pprintf(prn, "Number of iterations = %d\n", (int) iters);

    if (info == 0) {
	int k;

	pprintf(prn, "Convergence reached, with tolerance = %g\n", 
	       amax[1]);
	pputs(prn, "\nRegression coefficient estimates:\n");
	for (i=1; i<=1+nexo; i++) {
	    pprintf(prn, "     [%d]: %#14.6g (%#.6g)\n", i, amax[i+1],
		   amax[i+1+ncoeff]);
	}
	k = i;
	pputs(prn, "\nGARCH coefficient estimates:\n");
	for (i=k; i<=ncoeff; i++) {
	    if (i==k) {
		pputs(prn, "alpha[0]: ");
	    } else if (i-k <= p) {
		pprintf(prn, "alpha[%d]: ", i-k);
	    } else {
		pprintf(prn, " beta[%d]: ", i-k-p);
	    }
	    pprintf(prn, "%#14.6g (%#.6g)\n", amax[i+1], amax[i+1+ncoeff]);
	}
	pputs(prn, "\n(Standard errors in parentheses)\n");
    }

    free(yobs);
    free(xobs);
    free(umc);
    free(ydet);
    free(yy);
    free(coeff);
    free(d);
    free(oldc);
    free(vc);
    free(res2);
    free(res);
    free(sigma);
    free(a);
    free(ystoc);
    free(amax);
    free(amin);
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


