/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "libgretl.h"
#include "libset.h"

static void get_range_and_mean (int t1, int t2, const double *x,
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

static int 
do_range_mean_plot (int n, const double **Z, double a, double b,
		    const char *vname)
{
    FILE *fp = NULL;
    int fitline = 0;
    int t, err;

    if ((err = gnuplot_init(PLOT_RANGE_MEAN, &fp))) {
	return err;
    }

    if (!na(a) && !na(b)) {
	fitline = 1;
    }

    fprintf(fp, "# for %s\n", vname);
    fputs("set nokey\n", fp);
    fprintf(fp, "set title '%s %s %s'\n", 
	    G_("range-mean plot for"), vname, 
	    (fitline)? G_("with least squares fit") : "");
    fprintf(fp, "set xlabel '%s'\nset ylabel '%s'\n",
	    G_("mean"), G_("range"));
    fputs("plot \\\n", fp);

    gretl_push_c_numeric_locale();

    if (fitline) {
	fprintf(fp, "%g+%g*x notitle w lines lt 2 ,\\\n", a, b);
    }
    fputs("'-' using 1:2 w points lt 1\n", fp);

    for (t=0; t<n; t++) {
	fprintf(fp, "%g %g\n", Z[2][t], Z[1][t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

/* drop first/last observations from sample if missing obs 
   encountered */

static int 
rm_adjust_sample (int v, const double **Z, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;

    for (t=t1min; t<t2max; t++) {
	if (na(Z[v][t])) t1min++;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	if (na(Z[v][t])) t2max--;
	else break;
    }

    *t1 = t1min; *t2 = t2max;

    return 0;
}

int range_mean_graph (int vnum, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    double **rmZ;
    DATAINFO *rminfo;
    MODEL rmmod;
    int rmlist[4] = { 3, 1, 0, 2 };
    int k, t, m, nsamp, err = 0;
    int start, end, extra;
    double mean, range, tpval;
    char startdate[OBSLEN], enddate[OBSLEN];
    double a, b;
    int t1, t2;

    t1 = pdinfo->t1;
    t2 = pdinfo->t2;
    rm_adjust_sample(vnum, Z, &t1, &t2);

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

    rminfo = create_auxiliary_dataset(&rmZ, 3, m);
    if (rminfo == NULL) {
	return E_ALLOC;
    }

    pprintf(prn, _("Range-mean statistics for %s\n"), 
	    pdinfo->varname[vnum]);
    pprintf(prn, _("using %d sub-samples of size %d\n\n"),
	    m, k);

    pprintf(prn, "%30s%16s\n", _("range"), _("mean"));

    /* find group means and ranges */
    for (t=0; t<m; t++) {
	char obsstr[32];

	start = t1 + t * k;
	end = start + k - 1;
	if (end > t2) {
	    end = t2;
	} else if (t2 - end <= extra && extra < 3) {
	    end += extra;
	}

	get_range_and_mean(start, end, Z[vnum], &range, &mean);
	rmZ[1][t] = range;
	rmZ[2][t] = mean;

	ntodate(startdate, start, pdinfo);
	ntodate(enddate, end, pdinfo);
	sprintf(obsstr, "%s - %s", startdate, enddate);
	pputs(prn, obsstr);
	bufspace(20 - strlen(obsstr), prn);
	gretl_print_fullwidth_double(rmZ[1][t], GRETL_DIGITS, prn);
	gretl_print_fullwidth_double(rmZ[2][t], GRETL_DIGITS, prn);
	pputc(prn, '\n');
    }

    strcpy(rminfo->varname[1], "range");
    strcpy(rminfo->varname[2], "mean");

    rmmod = lsq(rmlist, &rmZ, rminfo, OLS, OPT_A);

    a = b = NADBL;

    if ((err = rmmod.errcode)) {
	pputs(prn, _("Error estimating range-mean model\n"));
	errmsg(err, prn);
    } else {
	pputc(prn, '\n');
	pprintf(prn, _("slope of range against mean = %g\n"),
		rmmod.coeff[1]);
	if (rmmod.sderr[1] > 0) {
	    tpval = student_pvalue_2(rmmod.coeff[1] / rmmod.sderr[1], rmmod.dfd);
	    pprintf(prn, _("p-value for H0: slope = 0 is %g\n"), tpval);
	} else {
	    tpval = 1.0;
	}

	if (tpval < .10) {
	    a = rmmod.coeff[0];
	    b = rmmod.coeff[1];
	} 
    }

    if (!gretl_in_batch_mode() && !gretl_looping()) {
	err = do_range_mean_plot(m, (const double **) rmZ, a, b, 
				 pdinfo->varname[vnum]);
    }

    clear_model(&rmmod);

    destroy_auxiliary_dataset(rmZ, rminfo);

    return err;
}
