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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"

static void get_range_and_mean (int t1, int t2, double *x,
				double *range, double *mean)
{
    double me = 0.0, mi = x[t1], ma = x[t1];
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) continue;
	ma = (x[t] > ma)? x[t] : ma;
	mi = (x[t] < mi)? x[t] : mi;
	me += x[t];
	n++;
    }

    if (n > 0) {
	*mean = me / n;
	*range = ma - mi;
    } else {
	*mean = NADBL;
	*range = NADBL;
    }
}

static int do_range_mean_plot (int n, double **Z, double *yhat,
			       const char *varname, PATHS *ppaths)
{
    FILE *fp = NULL;
    int t;

    if (gnuplot_init(ppaths, PLOT_RANGE_MEAN, &fp)) return E_FOPEN;

    fprintf(fp, "# range-mean plot for %s\n", varname);
    fputs("set nokey\n", fp);
    fprintf(fp, "set title '%s %s %s'\n", 
	    I_("range-mean plot for"), varname, 
	    (yhat == NULL)? "" : I_("with least squares fit"));
    fprintf(fp, "set xlabel '%s'\nset ylabel '%s'\n",
	    I_("mean"), I_("range"));
    fputs("plot \\\n'-' using 1:2 w points", fp);
    if (yhat != NULL) {
	fputs(" ,\\\n'-' using 1:2 w lines\n", fp);
    } else {
	fputc('\n', fp);
    }

    /* send data inline */
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    for (t=0; t<n; t++) {
	fprintf(fp, "%g %g\n", Z[2][t], Z[1][t]);
    }
    fputs("e\n", fp);
    if (yhat != NULL) {
	for (t=0; t<n; t++) {
	    fprintf(fp, "%g %g\n", Z[2][t], yhat[t]);
	}
	fputs("e\n", fp);
    }
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return 0;
}

static int adjust_t1t2 (int varnum, double **Z, int *t1, int *t2)
     /* drop first/last observations from sample if missing obs 
        encountered */
{
    int t, t1min = *t1, t2max = *t2;
    double xx;

    for (t=t1min; t<t2max; t++) {
	xx = Z[varnum][t];
	if (na(xx)) t1min += 1;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	xx = Z[varnum][t];
	if (na(xx)) t2max -= 1;
	else break;
    }

    *t1 = t1min; *t2 = t2max;
    return 0;
}

int range_mean_graph (int varnum, double **Z, DATAINFO *pdinfo, 
		      PRN *prn, PATHS *ppaths)
{
    double **rmZ;
    DATAINFO *rminfo;
    MODEL rmmod;
    int rmlist[4] = { 3, 1, 0, 2 };
    int k, t, m, nsamp, err = 0;
    int start, end, extra, len;
    double mean, range, tpval, *yhat = NULL;
    char startdate[OBSLEN], enddate[OBSLEN];
    int t1, t2;

    t1 = pdinfo->t1;
    t2 = pdinfo->t2;
    adjust_t1t2(varnum, Z, &t1, &t2);

    nsamp = t2 - t1 + 1;

    if (nsamp < 16) {
	pputs(prn, _("Sample is too small for range-mean graph\n"));
	errmsg(err, prn);
	return 1;
    } 	

    if (pdinfo->pd > 1 && nsamp >= 3 * pdinfo->pd) {
	k = pdinfo->pd;
    } else {
	k = (int) sqrt((double) nsamp);
    }

    extra = nsamp % k;
    m = (nsamp / k) + ((extra >= 3)? 1 : 0);

    rminfo = create_new_dataset(&rmZ, 3, m, 0);
    if (rminfo == NULL) return E_ALLOC;

    pprintf(prn, _("Range-mean statistics for %s\n"), 
	    pdinfo->varname[varnum]);
    pprintf(prn, _("using %d sub-samples of size %d\n\n"),
	    m, k);

    ntodate(startdate, t1, pdinfo);
    len = strlen(startdate) * 2 + 3 + 11;

    pprintf(prn, "%*s%16s\n", len, _("range"), _("mean"));

    /* find group means and ranges */
    for (t=0; t<m; t++) {
	start = t1 + t * k;
	end = start + k - 1;
	if (end > t2) {
	    end = t2;
	} else if (t2 - end <= extra && extra < 3) {
	    end += extra;
	}

	get_range_and_mean(start, end, Z[varnum], &range, &mean);
	rmZ[1][t] = range;
	rmZ[2][t] = mean;

	ntodate(startdate, start, pdinfo);
	ntodate(enddate, end, pdinfo);
	pprintf(prn, "%s - %s  ", startdate, enddate);
	gretl_print_fullwidth_double(rmZ[1][t], GRETL_DIGITS, prn);
	gretl_print_fullwidth_double(rmZ[2][t], GRETL_DIGITS, prn);
	pputs(prn, "\n");
    }

    strcpy(rminfo->varname[1], "range");
    strcpy(rminfo->varname[2], "mean");

    rmmod = lsq(rmlist, &rmZ, rminfo, OLS, OPT_A, 0.0);

    if ((err = rmmod.errcode)) {
	pputs(prn, _("Error estimating range-mean model\n"));
	errmsg(err, prn);
    } else {
	pputs(prn, "\n");
	pprintf(prn, _("slope of range against mean = %g\n"),
		rmmod.coeff[1]);
	if (rmmod.sderr[1] > 0) {
	    tpval = tprob(rmmod.coeff[1] / rmmod.sderr[1], rmmod.dfd);
	    pprintf(prn, _("p-value for H0: slope = 0 is %g\n"), tpval);
	} else {
	    tpval = 1.0;
	}

	if (tpval < .10) {
	    yhat = rmmod.yhat;
	} 
    }

    err = do_range_mean_plot(m, rmZ, yhat, pdinfo->varname[varnum], 
			     ppaths);

    clear_model(&rmmod);
    free_Z(rmZ, rminfo);
    clear_datainfo(rminfo, CLEAR_FULL);
    free(rminfo);

    return err;
}
