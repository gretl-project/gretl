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

int do_fcp (const int *list, const double **Z, 
	    const DATAINFO *pdinfo, PRN *prn)
{
    int t1, t2;
    double *yobs, **X;
    int nx, nobs;
    double *ydet;
    int ncoeff;
    double *coeff;
    double *vc, *res, *res2;
    double *ystoc;
    double *amax;
    double *b;
    double oldc, yy;
    int err = 0, iters = 0;
    int maxlag;
    int nobsmod;
    
    int i, p, q, ynum;

    /* FIXME: how exactly should these be set? */
    nx = 0;
    ncoeff = 1;

    p = list[1];
    q = list[2];
    ynum = list[4];
    maxlag = (p > q)? p : q; 

    t1 = 0; 
    t2 = pdinfo->t2;
    nobs = t2 + 1;
    nobsmod = nobs + maxlag;

    yobs = malloc(nobsmod * sizeof *yobs);
    ydet = malloc(nobsmod * sizeof *ydet);
    ystoc = malloc(nobsmod * sizeof *ystoc);

    res2 = malloc(nobsmod * sizeof *res2);
    for (i=0; i<nobsmod; i++) {
	res2[i] = 0.0;
    }

    res = malloc(nobsmod * sizeof *res);
    for (i=0; i<nobsmod; i++) {
	res[i] = 0.0;
    }    

    oldc = yy = 0.0;

    amax = malloc(nobsmod * sizeof *amax);
    for (i=0; i<nobsmod; i++) {
	amax[i] = 0.0;
    }

    coeff = malloc(ncoeff * sizeof *coeff);
    b = malloc(ncoeff * sizeof *b);
    for (i=0; i<ncoeff; i++) {
	coeff[i] = b[i] = 0.0;
    }    

    vc = malloc((ncoeff * ncoeff) * sizeof *vc);
    for (i=0; i<ncoeff; i++) {
	vc[i] = 0.0;
    } 

    if (nx > 0) {
	X = malloc(nx * sizeof *X);
	/* FIXME */
    } else {
	X = NULL;
    }

    for (i=0; i<maxlag; i++) {
	ystoc[i] = ydet[i] = yobs[i] = 0.0;
    }    
    for (i=maxlag; i<nobsmod; i++) {
	ystoc[i] = ydet[i] = yobs[i] = Z[ynum][i-maxlag];
    }

    /* initialize at unconditional mean of y */
    coeff[0] = _esl_mean(t1, t2, Z[ynum]);

    /* initialize elements of alpha, beta such that 
       alpha_0/(1 - alpha_1 - beta_1) = unconditional
       variance of y (?)
    */
    amax[0] = _esl_variance(t1, t2, Z[ynum]);
    amax[1] = p;
    amax[2] = q; 
    for (i=0; i<p+q; i++) {
	/* initial alpha, beta values */
	amax[3+i] = 0.1;
    }

    /* Need to set t1 high enough to allow for lags? */

    err = vsanal_(t1 + maxlag, t2 + maxlag, 
	    yobs, nobsmod, 
	    (const double **) X, nx, 
	    ydet, &yy, 
	    coeff, ncoeff, 
	    &oldc, vc, 
	    res2, 
	    res, 
	    ystoc, 
	    amax, b, &iters, prn);

    if (err != 0) {
	fprintf(stderr, "vsanal returned %d\n", err);
    }

    pprintf(prn, "Number of iterations = %d\n", iters);

    if (err == 0) {
	int nparam = ncoeff + p + q + 1;

	pprintf(prn, "Convergence reached, with tolerance = %g\n", 
	       amax[0]);

	for (i=1; i<=nparam; i++) {
	    pprintf(prn, "param[%d]: %#14.6g (%#.6g)\n", i, amax[i], amax[i+nparam]);
	}

	pputs(prn, "\n(Standard errors in parentheses)\n");
    }

    free(yobs);
    free(X); /* FIXME */
    free(ydet);
    free(coeff);
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
