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

#define MINSAMP 8
#define HDEBUG 0

static int
do_hurst_plot (int n, DATASET *dset, const MODEL *pmod,
	       const char *vname)
{
    FILE *fp;
    int t, err = 0;

    fp = open_plot_input_file(PLOT_HURST, 0, &err);
    if (err) {
	return err;
    }

    fprintf(fp, "# for %s\n", vname);
    fputs("set nokey\n", fp);
    fprintf(fp, "set title '%s %s'\n", _("Rescaled-range plot for"), vname);
    fprintf(fp, "set xlabel '%s'\n", _("log(sample size)"));
    fprintf(fp, "set ylabel '%s'\n", _("log(RS)"));

    gretl_push_c_numeric_locale();

    fputs("plot \\\n", fp);
    fprintf(fp, "%.10g+%.10g*x notitle w lines lt 2 ,\\\n",
	    pmod->coeff[0], pmod->coeff[1]);
    fputs("'-' using 1:2 w points lt 1\n", fp);

    for (t=0; t<n; t++) {
	fprintf(fp, "%.10g %.10g\n", dset->Z[2][t], dset->Z[1][t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static double cum_range (const double *x, int n, double xbar)
{
    double w, wmin, wmax;
    int i;

    w = wmin = wmax = 0.0;

    for (i=1; i<n; i++) {
	w += x[i-1] - xbar;
	if (w > wmax) {
	    wmax = w;
	} else if (w < wmin) {
	    wmin = w;
	}
    }

    return wmax - wmin;
}

static double stdev (const double *x, int n, double xbar)
{
    double dev, ssx = 0.0;
    int i;

    for (i=0; i<n; i++) {
	dev = x[i] - xbar;
	ssx += dev * dev;
    }

    return ssx > 0.0 ? sqrt(ssx / n) : 0.0;
}

static int hurst_calc (const double *x, int n, int depth,
		       DATASET *dset, PRN *prn)
{
    const char *heads[] = {
	N_("Size"),
	N_("RS(avg)"),
	N_("log(Size)"),
	N_("log(RS)")
    };
    int cw[] = { 5, 11, 11, 11 };
    int m, i, j;

    get_column_widths(heads, cw, 4);

    pprintf(prn, "%*s %*s %*s %*s\n",
	    UTF_WIDTH(_(heads[0]), cw[0]), _(heads[0]),
	    UTF_WIDTH(_(heads[1]), cw[1]), _(heads[1]),
	    UTF_WIDTH(_(heads[2]), cw[2]), _(heads[2]),
	    UTF_WIDTH(_(heads[3]), cw[3]), _(heads[3]));

    for (i=0, m=n; i<depth; i++, m/=2) {
	double RS = 0.0;
	int nsub = n / m;

#if HDEBUG
	fprintf(stderr, "nsub = %d\n", nsub);
	fprintf(stderr, "calculating at m = %d...\n", m);
#endif

	for (j=0; j<nsub; j++) {
	    double xbar, r, s;

	    xbar = gretl_mean(0, m-1, x + j*m);
	    r = cum_range(x + j*m, m, xbar);
	    s = stdev(x + j*m, m, xbar);
#if HDEBUG
	    fprintf(stderr, "range x + %d (%d) = %g\n", j*m, m, r);
	    fprintf(stderr, "stdev x + %d (%d) = %g\n", j*m, m, s);
#endif
	    RS += r / s;
	}

	RS /= nsub;

	dset->Z[1][i] = log2(RS);
	dset->Z[2][i] = log2(m);

	pprintf(prn, "%*d %#*.5g %#*.5g %#*.5g\n", cw[0], m, cw[1], RS,
		cw[2], dset->Z[2][i], cw[3], dset->Z[1][i]);
    }

    return 0;
}

static int get_depth (int T)
{
    int m = T;
    int depth = 0;

    while (m >= MINSAMP) {
	m /= 2;
	depth++;
    }

    return depth;
}

/* drop first/last observations from sample if missing obs
   encountered */

static int h_adjust_t1t2 (int v, const DATASET *dset,
			  int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;
    int miss = 0;

    for (t=t1min; t<t2max; t++) {
	if (na(dset->Z[v][t])) t1min++;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	if (na(dset->Z[v][t])) t2max--;
	else break;
    }

    *t1 = t1min; *t2 = t2max;

    for (t=t1min; t<t2max; t++) {
	if (na(dset->Z[v][t])) {
	    miss = 1;
	    break;
	}
    }

    return miss;
}

int hurst_exponent (int vnum, const DATASET *dset, gretlopt opt,
		    PRN *prn)
{
    DATASET *hset;
    MODEL hmod;
    int hlist[4] = {3, 1, 0, 2};
    int k, T;
    int t1, t2;
    int missing;
    int err = 0;

    t1 = dset->t1;
    t2 = dset->t2;

    missing = h_adjust_t1t2(vnum, dset, &t1, &t2);

    if (missing) {
	pputs(prn, _("There were missing data values"));
	pputc(prn, '\n');
	return E_MISSDATA;
    }

    T = t2 - t1 + 1;

    if (T < 96) {
	pputs(prn, _("Sample is too small for Hurst exponent"));
	pputc(prn, '\n');
	return E_TOOFEW;
    }

    k = get_depth(T);

    hset = create_auxiliary_dataset(3, k, 0);
    if (hset == NULL) {
	return E_ALLOC;
    }

    pprintf(prn, _("Rescaled range figures for %s"),
	    dset->varname[vnum]);
    pputc(prn, '\n');
    pputs(prn, _("(logs are to base 2)"));
    pputs(prn, "\n\n");

    /* do the rescaled range calculations */
    hurst_calc(dset->Z[vnum] + t1, T, k, hset, prn);

    strcpy(hset->varname[1], "RSavg");
    strcpy(hset->varname[2], "size");

    hmod = lsq(hlist, hset, OLS, OPT_A);

    if ((err = hmod.errcode)) {
	pputs(prn, _("Error estimating Hurst exponent model\n"));
	errmsg(err, prn);
    } else {
	gretl_matrix *result = NULL;
	char **colnames = NULL;
	const char *heads[] = {
	    N_("Intercept"),
	    N_("coeff"),
	    N_("std. error")
	};
	int cw[] = {12, 12, 12};

	get_column_widths(heads, cw, 3);

	pprintf(prn, "\n%s (n = %d)\n\n", _("Regression results"), k);
	pprintf(prn, "%*s %*s %*s\n", cw[0], "",
		UTF_WIDTH(_(heads[1]), cw[1]), _(heads[1]),
		UTF_WIDTH(_(heads[2]), cw[2]), _(heads[2]));
	pprintf(prn, "%*s %#*.5g %#*.5g\n",
		UTF_WIDTH(_(heads[0]), cw[0]), _(heads[0]),
		cw[1], hmod.coeff[0],
		cw[2], hmod.sderr[0]);
	pprintf(prn, "%*s %#*.5g %#*.5g\n",
		UTF_WIDTH(_("Slope"), cw[0]), _("Slope"),
		cw[1], hmod.coeff[1],
		cw[2], hmod.sderr[1]);
	pputc(prn, '\n');
	pprintf(prn, "%s = %g\n", _("Estimated Hurst exponent"), hmod.coeff[1]);

	/* store estimate for $result accessor */
	result = gretl_matrix_alloc(1, 2);
	gretl_matrix_set(result, 0, 0, hmod.coeff[1]);
	gretl_matrix_set(result, 0, 1, hmod.sderr[1]);
	colnames = strings_array_new(2);
	colnames[0] = gretl_strdup("hurst");
	colnames[1] = gretl_strdup("se");
	gretl_matrix_set_colnames(result, colnames);
	set_last_result_data(result, GRETL_TYPE_MATRIX);
    }

    if (!err && gnuplot_graph_wanted(PLOT_HURST, opt, &err)) {
	err = do_hurst_plot(k, hset, &hmod, dset->varname[vnum]);
    }

    clear_model(&hmod);
    destroy_dataset(hset);

    return err;
}
