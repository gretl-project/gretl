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

    if (gnuplot_init(ppaths, &fp)) return E_FOPEN;

    /* FIXME title and axis labels */

    fprintf(fp, "# range-mean plot for %s\n", varname);
    fprintf(fp, "set nokey\n");
    fprintf(fp, "set title '%s %s'\n", I_("range-mean plot for"),
	    varname);
    fprintf(fp, "set xlabel '%s'\nset ylabel '%s\n",
	    I_("mean"), I_("range"));
    fprintf(fp, "plot \\\n'-' using 1:2 w points");
    if (yhat != NULL) {
	fprintf(fp, " ,\\\n'-' using 1:2 w lines\n");
    } else {
	fprintf(fp, "\n");
    }

    /* send data inline */
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    for (t=0; t<n; t++) {
	fprintf(fp, "%g %g\n", Z[2][t], Z[1][t]);
    }
    fprintf(fp, "e\n");
    if (yhat != NULL) {
	for (t=0; t<n; t++) {
	    fprintf(fp, "%g %g\n", Z[2][t], yhat[t]);
	}
    }
    fprintf(fp, "e\n");
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fprintf(fp, "pause -1\n");
#endif
    fclose(fp);
    return 0;
}

int range_mean_graph (int varnum, double **Z, DATAINFO *pdinfo, 
		      PRN *prn, PATHS *ppaths)
{
    double **rmZ;
    DATAINFO *rminfo;
    MODEL rmmod;
    int rmlist[4] = { 3, 1, 2, 0};
    int k, t, m, nsamp, err = 0;
    int start, end, extra, len;
    double mean, range, tpval, *yhat;
    char startdate[9], enddate[9];

    nsamp = pdinfo->t2 - pdinfo->t1 + 1;

    if (nsamp < 16) {
	pprintf(prn, _("Sample is too small for range-mean graph\n"));
	errmsg(err, prn);
	return 1;
    } 	

    if (pdinfo->pd > 1 && nsamp >= 3 * pdinfo->pd) {
	k = pdinfo->pd;
    } else if (nsamp >= 30) {
	k = 10;
    } else {
	k = 5;
    }

    extra = nsamp % k;
    m = (nsamp / k) + ((extra >= 3)? 1 : 0);

    rminfo = create_new_dataset(&rmZ, 3, m, 0);
    if (rminfo == NULL) return E_ALLOC;
    rminfo->extra = 1;

    pprintf(prn, _("Range-mean statistics for %s\n"), 
	    pdinfo->varname[varnum]);
    pprintf(prn, _("using %d sub-samples of size %d\n\n"),
	    m, k);

    ntodate(startdate, pdinfo->t1, pdinfo);
    len = strlen(startdate) * 2 + 3 + 11;

    pprintf(prn, "%*s%16s\n", len, _("range"), _("mean"));

    /* find group means and ranges */
    for (t=0; t<m; t++) {
	start = pdinfo->t1 + t * k;
	end = start + k;
	if (end > pdinfo->t2) {
	    end = pdinfo->t2;
	} else if (pdinfo->t2 - end <= extra && extra < 3) {
	    end += extra;
	}
	end = (start + k > pdinfo->t2)? pdinfo->t2 : start + k;

	get_range_and_mean(start, end, Z[varnum], &range, &mean);
	rmZ[1][t] = range;
	rmZ[2][t] = mean;

	ntodate(startdate, start, pdinfo);
	ntodate(enddate, end, pdinfo);
	pprintf(prn, "%s - %s  ", startdate, enddate);
	gretl_print_fullwidth_double(rmZ[1][t], GRETL_DIGITS, prn);
	gretl_print_fullwidth_double(rmZ[2][t], GRETL_DIGITS, prn);
	pprintf(prn, "\n");
    }

    strcpy(rminfo->varname[1], "range");
    strcpy(rminfo->varname[2], "mean");

    rmmod = lsq(rmlist, &rmZ, rminfo, OLS, 0, 0.0);
    if ((err = rmmod.errcode)) {
	pprintf(prn, _("Error estimating range-mean model\n"));
	errmsg(err, prn);
    } else { /* just testing */
	/* printmodel(&rmmod, rminfo, prn); */
    }

    pprintf(prn, "\n");
    pprintf(prn, _("slope of range against mean = %g\n"),
	    rmmod.coeff[1]);
    tpval = tprob(rmmod.coeff[1] / rmmod.sderr[1], rmmod.dfd);
    pprintf(prn, _("p-value for H0: slope = 0 is %g\n"), tpval);

    if (tpval < .10) {
	yhat = rmmod.yhat;
    } else {
	yhat = NULL;
    }

    err = do_range_mean_plot(m, rmZ, yhat, pdinfo->varname[varnum], 
			     ppaths);

    clear_model(&rmmod, NULL);
    free_Z(rmZ, rminfo);
    clear_datainfo(rminfo, CLEAR_FULL);
    free(rminfo);

    return err;
}
