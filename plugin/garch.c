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
    integer ninit, nfinsm, ifrom, ito;
    integer inyear, nfinyr;
    doublereal *yobs, *xobs;
    integer nend, nexo;
    integer iread;
    integer nstoch;
    doublereal *umc, *ydet, *yy;
    integer ncoeff;
    doublereal *coeff, *d, *oldc;
    doublereal *vc, *res2;
    doublereal *res, *sigma;
    doublereal *a;
    doublereal *ystoc;
    doublereal *amax, *amin;
    doublereal *b;
    integer ncoefb;
    integer iters, info;
    
    int i, nobs;
    int p, q, ynum;

    inyear = 1;
    nfinyr = 99;
    ifrom = 1;
    nend = 1;
    nexo = 0;
    nstoch = 1;
    ncoeff = 1;
    ncoefb = 1;

    /* FIXME */
    ito = nobs = pdinfo->t2 - pdinfo->t1 + 1;

    p = list[1];
    q = list[2];
    ynum = list[4];

    if (ifrom < inyear) inyear = ifrom;
    if (ito > nfinyr) nfinyr = ito;
    iread = nfinyr - inyear + 1;

    ninit = ifrom - inyear + 1;
    nfinsm = ito - inyear + 1;
    
    yobs = malloc((nend * iread + 1) * sizeof *yobs);
    ydet = malloc((nend * iread + 1) * sizeof *ydet);
    ystoc = malloc((nend * iread + 1) * sizeof *ystoc);

    res2 = malloc((nend * iread + 1) * sizeof *res2);
    for (i=1; i<=nend*iread; i++) {
	res2[i] = 0.0;
    }

    umc = malloc((nstoch + 1) * sizeof *umc);
    for (i=1; i<=nstoch; i++) {
	umc[i] = 0.0;
    }

    res = malloc((nstoch * iread + 1) * sizeof *res);
    for (i=1; i<=nstoch*iread; i++) {
	res[i] = 0.0;
    }    

    sigma = malloc((nstoch * nstoch + 1) * sizeof *sigma);
    for (i=1; i<=nstoch*nstoch; i++) {
	sigma[i] = 0.0;
    }

    oldc = malloc((nend + 1) * sizeof *oldc);
    yy = malloc((nend + 1) * sizeof *yy);
    for (i=1; i<=nend; i++) {
	oldc[i] = yy[i] = 0.0;
    } 

    a = malloc((nend * nend + 1) * sizeof *a);
    for (i=1; i<=nend*nend; i++) {
	a[i] = 0.0;
    }

    amax = malloc((nend * iread + 1) * sizeof *amax);
    amin = malloc((nend * iread + 1) * sizeof *amin);
    for (i=1; i<=nend*iread; i++) {
	amax[i] = amin[i] = 0.0;
    }

    coeff = malloc((ncoeff + 1) * sizeof *coeff);
    b = malloc((ncoeff + 1) * sizeof *b);
    for (i=1; i<=ncoeff; i++) {
	coeff[i] = b[i] = 0.0;
    }    

    d = malloc((nend * ncoeff + 1) * sizeof *d);
    vc = malloc((ncoeff * ncoeff + 1) * sizeof *vc);
    for (i=1; i<=nend*ncoeff; i++) {
	d[i] = vc[i] = 0.0;
    } 

    xobs = NULL;

    for (i=1; i<=iread; i++) {
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

    vsanal_(&ninit, &nfinsm, &yobs[1], &nend, &iread,
	    xobs, &nexo, &umc[1], &nstoch, &ydet[1], 
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
	for (i=1; i<=nend+nexo; i++) {
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


