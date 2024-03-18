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
#include "version.h"
#include "libset.h"

static int get_range_and_mean (int t1, int t2, const double *x,
			       double *range, double *mean,
			       gretlopt opt)
{
    double xsum = 0.0, xmin = x[t1], xmax = x[t1];
    int t, n = 0;
    int err = 0;

    if (opt & OPT_T) {
	/* trim the extrema */
	int T = t2 - t1 + 1;
	double *xs = malloc(T * sizeof *xs);

	if (xs == NULL) {
	    err = E_ALLOC;
	} else {
	    int s = 0;

	    for (t=t1; t<=t2; t++) {
		if (!na(x[t])) {
		    xs[s++] = x[t];
		    n++;
		}
	    }

	    if (n < 4) {
		err = E_TOOFEW;
		n = 0;
	    } else {
		qsort(xs, n, sizeof *x, gretl_compare_doubles);
		xmin = xs[1];
		xmax = xs[n-2];
		for (s=1; s<n-1; s++) {
		    xsum += xs[s];
		}
		n -= 2;
	    }

	    free(xs);
	}
    } else {
	/* no trimming: use all available data */
	for (t=t1; t<=t2; t++) {
	    if (!na(x[t])) {
		xmax = (x[t] > xmax) ? x[t] : xmax;
		xmin = (x[t] < xmin) ? x[t] : xmin;
		xsum += x[t];
		n++;
	    }
	}
    }

    if (n > 0) {
	*mean = xsum / n;
	*range = xmax - xmin;
    } else {
	*mean = NADBL;
	*range = NADBL;
	if (!err) {
	    err = E_DATA;
	}
    }

    return err;
}

static int do_range_mean_plot (const gretl_matrix *y,
			       const gretl_matrix *X,
			       double b0, double b1,
			       const char *vname)
{
    FILE *fp;
    int fitline = 0;
    int t, err = 0;

    fp = open_plot_input_file(PLOT_RANGE_MEAN, 0, &err);
    if (err) {
	return err;
    }

    if (!na(b0) && !na(b1)) {
	fitline = 1;
    }

    fprintf(fp, "# for %s\n", vname);
    fputs("set nokey\n", fp);
    fprintf(fp, "set title '%s %s %s'\n", 
	    _("range-mean plot for"), vname, 
	    (fitline)? _("with least squares fit") : "");
    fprintf(fp, "set xlabel '%s'\nset ylabel '%s'\n",
	    _("mean"), _("range"));
    fputs("plot \\\n", fp);

    gretl_push_c_numeric_locale();

    if (fitline) {
	fprintf(fp, "%g+%g*x notitle w lines lt 2 ,\\\n", b0, b1);
    }
    fputs("'-' using 1:2 w points lt 1\n", fp);

    for (t=0; t<X->rows; t++) {
	fprintf(fp, "%g %g\n", gretl_matrix_get(X, t, 1),
		gretl_vector_get(y, t));
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* drop first/last observations from sample if missing obs 
   encountered */

static int 
rm_adjust_sample (int v, const DATASET *dset, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;

    for (t=t1min; t<t2max; t++) {
	if (na(dset->Z[v][t])) t1min++;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	if (na(dset->Z[v][t])) t2max--;
	else break;
    }

    *t1 = t1min; 
    *t2 = t2max;

    return 0;
}

static int get_n_samples (int T, int pd, int mmin)
{
    int k = (int) sqrt((double) T);

    if (k < mmin) {
	k = mmin;
    } else {
	/* If pd >= mmin and is "not too far" from k, it may make sense to
	   use pd instead of sqrt(T) as the subsample length. This will,
	   e.g., make each subsample a year, with quarterly or monthly
	   data.
	*/
	if (pd >= mmin && T/pd >= 4 && 
	    pd >= 2.0*k/3.0 && pd <= 3.0*k/2.0) {
	    k = pd;
	}
    }

    return k;
}

int range_mean_graph (int vnum, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    gretl_matrix_block *B;
    gretl_matrix *y, *X, *b, *V;
    int k, t, m, T, mmin, rem;
    int quiet = (opt & OPT_Q);
    char startdate[OBSLEN], enddate[OBSLEN];
    int t1 = dset->t1;
    int t2 = dset->t2;
    int digits;
    int err = 0;

    rm_adjust_sample(vnum, dset, &t1, &t2);
    T = t2 - t1 + 1;

    /* We need at least 4 sub-samples to do this, and each
       must have at least @mmin observations, the specific
       value depending on whether or not we're trimming the
       extrema in each sub-sample.
    */
    mmin = (opt & OPT_T)? 6 : 4;
    if (T < 4 * mmin) {
	pputs(prn, _("Sample is too small for range-mean graph\n"));
	return E_TOOFEW;
    } 

    k = get_n_samples(T, dset->pd, mmin);
    rem = T % k;
    m = (T / k) + ((rem >= mmin)? 1 : 0);

    B = gretl_matrix_block_new(&y, m, 1,
			       &X, m, 2,
			       &b, 2, 1,
			       &V, 2, 2,
			       NULL);
    if (B == NULL) {
	return E_ALLOC;
    }

    if (!quiet) {
	pprintf(prn, _("Range-mean statistics for %s\n"), 
		dset->varname[vnum]);
	pprintf(prn, _("using %d sub-samples of size %d\n\n"),
		m, k);
	pprintf(prn, "%30s%16s\n", _("range"), _("mean"));
    }

    digits = get_gretl_digits();

    /* find sub-sample means and ranges */

    for (t=0; t<m; t++) {
	int start = t1 + t * k;
	int end = start + k - 1; 
	double range, mean;

	if (end > t2) {
	    end = t2;
	} else if (t2 - end <= rem && rem < mmin) {
	    end += rem;
	}

	err = get_range_and_mean(start, end, dset->Z[vnum], 
				 &range, &mean, opt);
	if (err) {
	    break;
	}

	gretl_vector_set(y, t, range);
	gretl_matrix_set(X, t, 0, 1.0);
	gretl_matrix_set(X, t, 1, mean);

	if (!quiet) {
	    int len;

	    ntolabel(startdate, start, dset);
	    ntolabel(enddate, end, dset);
	    len = pprintf(prn, "%s - %s", startdate, enddate);
	    bufspace(20 - len, prn);
	    gretl_print_fullwidth_double(range, digits, prn);
	    gretl_print_fullwidth_double(mean, digits, prn);
	    pputc(prn, '\n');
	}
    }

    if (!err) {
	double b0 = NADBL, b1 = NADBL; 
	double s2;

	err = gretl_matrix_ols(y, X, b, V, NULL, &s2);

	if (err) {
	    pputs(prn, _("Error estimating range-mean model\n"));
	    errmsg(err, prn);
	} else {
	    double bse = sqrt(gretl_matrix_get(V, 1, 1));
	    double tstat, pv;

	    if (!quiet) {
		pputc(prn, '\n');
		pprintf(prn, _("slope of range against mean = %g\n"),
			b->val[1]);
	    }

	    if (bse > 0) {
		tstat = b->val[1] / bse;
		pv = student_pvalue_2(m - 2, tstat);
		record_test_result(tstat, pv);
		if (!quiet) {
		    pprintf(prn, _("p-value for H0: slope = 0 is %g\n"), pv);
		}
		if (pv < .10) {
		    /* arrange for fitted line to be drawn */
		    b0 = b->val[0];
		    b1 = b->val[1];
		} 
	    }
	}

	if (!err && !quiet) {
	    if (gnuplot_graph_wanted(PLOT_RANGE_MEAN, opt, &err)) {
		err = do_range_mean_plot(y, X, b0, b1, dset->varname[vnum]);
	    }
	}
    }

    gretl_matrix_block_destroy(B);

    return err;
}
