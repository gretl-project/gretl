/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"

/* To do: 
   - restrict this test to OLS
   - fix x-axis of graph (obs number when sub-sampled; time-series)
   - do something about menu item when lapack not available?
   - fix printout of h and influence
*/
  
static int leverage_plot (int n, const double *uhat, const double *h, 
			  const DATAINFO *pdinfo, PATHS *ppaths)
{
    FILE *fp = NULL;
    int t;

    if (gnuplot_init(ppaths, &fp)) return 1;

    fputs("# leverage/influence plot\n", fp);
    fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set xrange [%g:%g]\n", 0.5, n + 0.5);
    fputs("set nokey\n", fp); 

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    /* upper plot: leverage factor */
    fputs("set origin 0.0,0.50\n", fp);
    fputs("set yrange [0:1]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("leverage"));
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);
    for (t=0; t<n; t++) {
	fprintf(fp, "%d %g\n", t+1, h[t]);
    }
    fputs("e\n", fp);

    /* lower plot: influence factor */
    fputs("set origin 0.0,0.0\n", fp);
    fputs("set yrange [*:*]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("influence")); 
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);
    for (t=0; t<n; t++) {
	fprintf(fp, "%d %g\n", t+1, uhat[t] * h[t] / (1.0 - h[t]));
    }
    fputs("e\n", fp);
    fputs("set nomultiplot\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fprintf(fp, "pause -1\n");
#endif
    fclose(fp);
    return 0;
}

/* In fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are
   contiguous.
*/

int model_leverage (const MODEL *pmod, const double **Z, 
		    const DATAINFO *pdinfo, PRN *prn,
		    PATHS *ppaths)
{
    integer info, lwork;
    integer m, n, lda;
    gretl_matrix *Q;
    doublereal *tau, *work;
    double *h = NULL;
    int i, j, t;
    int err = 0;

    m = pmod->t2 - pmod->t1 + 1; /* # of rows = # of observations */
    lda = m;                     /* leading dimension of Q */
    n = pmod->list[0] - 1;       /* # of cols = # of variables */

    Q = gretl_matrix_alloc(m, n);
    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = malloc(sizeof *work);

    if (Q == NULL || tau == NULL || work == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* copy independent var values into Q */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    Q->val[j++] = Z[pmod->list[i]][t];
	}
    }

    /* do a workspace size query */
    lwork = -1;
    info = 0;
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* set up optimally sized work array */
    lwork = (integer) work[0];
    work = realloc(work, (size_t) lwork * sizeof *work);
    if (work == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* run actual QR factorization */
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* obtain the real "Q" matrix */
    dorgqr_(&m, &n, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	err = 1;
	goto qr_cleanup;
    }

    h = malloc(m * sizeof *h);
    if (h == NULL) {
	err = 1;
	goto qr_cleanup;
    }

    for (t=0; t<m; t++) {
	double infl;
	int tmod = t + pmod->t1;

	h[t] = 0.0;
	for (i=0; i<n; i++) {
	    double q = gretl_matrix_get(Q, t, i);

	    h[t] += q * q;
	}
	infl = pmod->uhat[tmod] * h[t] / (1.0 - h[t]);
	pprintf(prn, "%d: uhat = %6.3f h = %.3f, ratio = %6.3f\n", 
		t, pmod->uhat[tmod], h[t], infl);
    }

    leverage_plot(m, &pmod->uhat[pmod->t1], h, pdinfo, ppaths);
    free(h);

 qr_cleanup:
    gretl_matrix_free(Q);
    free(tau); free(work);

    return err;    
}

