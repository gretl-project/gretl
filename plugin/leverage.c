/*
 *  Copyright (c) 2003 by Allin Cottrell
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
   - do something about menu item when lapack not available?
*/
  
static int leverage_plot (int n, int tstart, const double *uhat, 
			  const double *h, double ***pZ, DATAINFO *pdinfo, 
			  PATHS *ppaths)
{
    FILE *fp = NULL;
    int t, tmod;
    int timeplot = 0;

    if (gnuplot_init(ppaths, &fp)) return 1;

    if (dataset_is_time_series(pdinfo) && 
	(pdinfo->pd == 1 || pdinfo->pd == 4 || pdinfo->pd == 12)) {
	char per[8];

	if (pdinfo->pd == 1) strcpy(per, "annual");
	else if (pdinfo->pd == 4) strcpy(per, "qtrs");
	else if (pdinfo->pd == 12) strcpy(per, "months");
	plotvar(pZ, pdinfo, per);
	timeplot = varindex(pdinfo, per);
    }

    fputs("# leverage/influence plot\n", fp);
    fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (!timeplot) { 
	fprintf(fp, "set xrange [%g:%g]\n", 
		tstart + 0.5, tstart + n + 0.5);
    }
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
	if (timeplot) {
	    tmod = t + tstart;
	    fprintf(fp, "%g %g\n", (*pZ)[timeplot][tmod], h[t]);
	} else { 
	    fprintf(fp, "%d %g\n", t+tstart+1, h[t]);
	}
    }
    fputs("e\n", fp);

    /* lower plot: influence factor */
    fputs("set origin 0.0,0.0\n", fp);
    fputs("set yrange [*:*]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("influence")); 
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);
    for (t=0; t<n; t++) {
	double f;

	tmod = t + tstart;
	f = uhat[tmod] * h[t] / (1.0 - h[t]);
	if (timeplot) {
	    fprintf(fp, "%g %g\n", (*pZ)[timeplot][tmod], f);
	} else {
	    fprintf(fp, "%d %g\n", tmod+1, f);
	}
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

int model_leverage (const MODEL *pmod, double ***pZ, 
		    DATAINFO *pdinfo, PRN *prn,
		    PATHS *ppaths)
{
    integer info, lwork;
    integer m, n, lda;
    gretl_matrix *Q;
    doublereal *tau, *work;
    double *hvec = NULL;
    double lp;
    int i, j, t;
    int err = 0, gotlp = 0;

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
	    Q->val[j++] = (*pZ)[pmod->list[i]][t];
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

    if (ppaths != NULL) {
	hvec = malloc(m * sizeof *hvec);
	if (hvec == NULL) {
	    err = 1;
	    goto qr_cleanup;
	}
    }

    pputs(prn, "        ");
    pprintf(prn, "%*s", UTF_WIDTH(_("residual"), 13), _("residual"));
    pprintf(prn, "%*s", UTF_WIDTH(_("leverage"), 13), _("leverage"));
    pprintf(prn, "%*s", UTF_WIDTH(_("influence"), 13), _("influence"));
    pputs(prn, "\n        ");
    pputs(prn, "         u         0<=h<=1    u*h/(1-h)\n\n");

    lp = 2.0 * n / m;

    for (t=0; t<m; t++) {
	double f, h = 0.0;
	int tmod = t + pmod->t1;

	for (i=0; i<n; i++) {
	    double q = gretl_matrix_get(Q, t, i);

	    h += q * q;
	}
	if (h > lp) gotlp = 1;
	f = pmod->uhat[tmod] * h / (1.0 - h);
	print_obs_marker(tmod, pdinfo, prn);
	pprintf(prn, "%12.5g %11.3f%s %12.5g\n", pmod->uhat[tmod], h, 
		(h > lp)? "*" : " ", f);
	if (hvec != NULL) hvec[t] = h;
    }

    if (gotlp) {
	pprintf(prn, "\n%s\n", _("('*' indicates a leverage point)"));
    } else {
	pprintf(prn, "\n%s\n", _("No leverage points were found"));
    }

    if (ppaths != NULL) {
	leverage_plot(m, pmod->t1, &pmod->uhat[pmod->t1], hvec, 
		      pZ, pdinfo, ppaths);
	free(hvec);
    }

 qr_cleanup:
    gretl_matrix_free(Q);
    free(tau); free(work);

    return err;    
}

