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

/* boxplots for gretl via gnuplot */

#include "libgretl.h"
#include "boxplots.h"
#include "libset.h"

typedef struct {
    double mean;
    double median;
    double conf[2];    /* bounds of confidence interval */
    double uq, lq;     /* upper and lower quartiles */
    double max, min;   /* data max and min */
    double wmax, wmin; /* whisker max and min */
    int n;
    char varname[VNAMELEN];
    char *bool;
    gretl_matrix *outliers;
} BOXPLOT;

enum {
    BOX_SHOW_MEAN = 1 << 0,
    BOX_INTERVALS = 1 << 1,
    BOX_WHISKBARS = 1 << 2,
    BOX_XTICS     = 1 << 3,
    BOX_FACTOR    = 1 << 4
};

typedef struct {
    int flags;
    int nplots;
    const int *list;
    char *title;
    int n_bools;
    BOXPLOT *plots;
    double gmin, gmax;
    double limit;
    double *x;
    gretl_matrix *dvals;
    char xlabel[VNAMELEN];
    char ylabel[VNAMELEN];
    const char *literal;
} PLOTGROUP;

#define show_mean(g)    (g->flags & BOX_SHOW_MEAN)
#define do_intervals(g) (g->flags & BOX_INTERVALS)
#define whiskerbars(g)  (g->flags & BOX_WHISKBARS)
#define use_xtics(g)    (g->flags & BOX_XTICS)
#define factorized(g)   (g->flags & BOX_FACTOR)

#define set_show_mean(g)    (g->flags |= BOX_SHOW_MEAN)
#define set_do_intervals(g) (g->flags |= BOX_INTERVALS)
#define set_use_xtics(g)    (g->flags |= BOX_XTICS)
#define set_box_factor(g)   (g->flags |= BOX_FACTOR)

#define unset_show_mean(g)    (g->flags &= ~BOX_SHOW_MEAN)
#define unset_do_intervals(g) (g->flags &= ~BOX_INTERVALS)

#define BPSTRLEN 128

static char boxplots_string[BPSTRLEN];

const char *get_last_boxplots_string (void)
{
    return boxplots_string;
}

static void quartiles_etc (double *x, int n, BOXPLOT *plt,
			   double limit)
{
    double p[3] = {0.25, 0.5, 0.75};

    plt->wmin = plt->min = x[0];
    plt->wmax = plt->max = x[n-1];

    gretl_array_quantiles(x, n, p, 3);

    plt->lq = p[0];
    plt->median = p[1];
    plt->uq = p[2];

    if (limit > 0) {
	double d = limit * (plt->uq - plt->lq);
	double xlo = plt->lq - d;
	double xhi = plt->uq + d;
	int i, j, nout = 0;

	for (i=0; i<n; i++) {
	    if (x[i] < xlo) {
		nout++;
	    } else if (x[i] > xhi) {
		nout++;
		if (i > 0 && x[i-1] <= xhi) {
		    /* the top of the upper whisker */
		    plt->wmax = x[i-1];
		}
	    } else if (i > 0 && x[i-1] < xlo) {
		/* the bottom of the lower whisker */
		plt->wmin = x[i-1];
	    }
	}

	if (nout > 0) {
	    plt->outliers = gretl_vector_alloc(nout);
	}

	if (plt->outliers != NULL) {
	    j = 0;
	    for (i=0; i<n; i++) {
		if (x[i] < xlo || x[i] > xhi) {
		    gretl_vector_set(plt->outliers, j++, x[i]);
		}
	    }
	}
    }
}

static void six_numbers (BOXPLOT *plt, int do_mean, PRN *prn)
{
    if (do_mean) {
	pprintf(prn, "%#9.5g", plt->mean);
    }

    pprintf(prn, "%#9.5g%#9.5g%#9.5g%#9.5g%#9.5g", plt->min, 
	    plt->lq, plt->median, plt->uq, plt->max);

    if (plt->n > 0) {
	pprintf(prn, "  (n=%d)\n", plt->n);
    } else {
	pputc(prn, '\n');
    }
}

static void five_numbers_with_interval (BOXPLOT *plt, PRN *prn)
{
    char tmp[32];

    sprintf(tmp, "%#.5g - %#.5g", plt->conf[0], plt->conf[1]);
    pprintf(prn, "%#8.5g%#10.5g%#10.5g%17s%#10.5g%#10.5g\n",
	    plt->min, plt->lq, plt->median,
	    tmp, plt->uq, plt->max);
}

static void box_stats_leader (PLOTGROUP *grp, int i, int offset, 
			      PRN *prn)
{
    BOXPLOT *plt = &grp->plots[i];

    if (factorized(grp)) {
	char numstr[32];

	sprintf(numstr, "%g", gretl_vector_get(grp->dvals, i));
	pprintf(prn, "%-*s", offset, numstr);
    } else if (plt->bool != NULL) {
	pprintf(prn, "%s\n %-*s", plt->varname, offset - 1, plt->bool);
    } else {
	pprintf(prn, "%-*s", offset, plt->varname);
    }
}

static int has_mean (PLOTGROUP *grp)
{
    int i;

    for (i=0; i<grp->nplots; i++) {
	if (na(grp->plots[i].mean)) {
	    return 0;
	}
    }

    return 1;
}

static int get_format_offset (PLOTGROUP *grp)
{
    int L = 6;
    int i, n;

    for (i=0; i<grp->nplots; i++) {
	if (grp->plots[i].bool != NULL) {
	    n = strlen(grp->plots[i].bool);
	} else {
	    n = strlen(grp->plots[i].varname);
	}
	if (n > L) {
	    L = n;
	}	    
    }

    return L + 2;
}

static int boxplot_print_stats (PLOTGROUP *grp, PRN *prn) 
{
    int do_mean = 0;
    int offset, pad, i;

    if (factorized(grp)) {
	offset = strlen(grp->xlabel) + 2;
    } else {
	offset = get_format_offset(grp);
    }

    if (do_intervals(grp)) {
	pprintf(prn, "%s\n\n", _("Numerical summary with bootstrapped confidence "
				 "interval for median"));
	pad = offset + 8;
    } else {
	if (*grp->ylabel != '\0') {
	    pprintf(prn, _("Numerical summary for %s"), grp->ylabel);
	    pputs(prn, "\n\n");
	} else {
	    pprintf(prn, "%s\n\n", _("Numerical summary"));
	}
	do_mean = has_mean(grp);
	pad = (do_mean)? (offset + 9) : (offset + 10);
    }

    if (factorized(grp)) {
	pputs(prn, grp->xlabel);
	pad -= strlen(grp->xlabel);
    }    

    if (do_intervals(grp)) {
	pprintf(prn, "%*s%10s%10s%17s%10s%10s\n", pad,
		"min", "Q1", _("median"), 
		/* xgettext:no-c-format */
		_("(90% interval)"), 
		"Q3", "max");
	for (i=0; i<grp->nplots; i++) {
	    box_stats_leader(grp, i, offset, prn);
	    five_numbers_with_interval(&grp->plots[i], prn);
	}
    } else {	
	if (do_mean) {
	    pprintf(prn, "%*s%9s%9s%9s%9s%9s\n", pad, _("mean"),
		    "min", "Q1", _("median"), "Q3", "max");
	} else {
	    pprintf(prn, "%*s%10s%10s%10s%10s\n", pad,
		    "min", "Q1", _("median"), "Q3", "max");
	}	    
	for (i=0; i<grp->nplots; i++) {
	    box_stats_leader(grp, i, offset, prn);
	    six_numbers(&grp->plots[i], do_mean, prn);
	}
    } 

    return 0;
}

#define ITERS 560
#define CONFIDENCE 0.90

/* obtain bootstrap estimate of 90% confidence interval
   for the sample median of data series x; return low and
   high values in 'low' and 'high' */

static int 
median_interval (double *x, int n, double *low, double *high)
{
    double *medians, *samp;
    double p[2];
    int i, j, t;
    int err;

    medians = malloc(ITERS * sizeof *medians);
    if (medians == NULL) return 1;

    samp = malloc(n * sizeof *samp);
    if (samp == NULL) {
	free(medians);
	return 1;
    }

    for (i=0; i<ITERS; i++) {
	/* sample with replacement from x */
	for (j=0; j<n; j++) {
	    t = gretl_rand_int_max(n);
	    samp[j] = x[t];
	}
	/* find the median of the sample */
	medians[i] = gretl_array_quantile(samp, n, 0.5);
    }

    p[0] = (1.0 - CONFIDENCE) / 2.0;
    p[1] = 1.0 - p[0];

    err = gretl_array_quantiles(medians, ITERS, p, 2);

    *low = p[0];
    *high = p[1];

    free(samp);
    free(medians);

    return err;
}

static int special_varcount (const char *s)
{
    char test[36];
    int n = 0;

    while (sscanf(s, "%35s", test) == 1) {
	if (*test != '(') n++;
	s += strspn(s, " ");
	s += strlen(test);
    }

    return n;
}

static void plotgroup_destroy (PLOTGROUP *grp)
{
    int i;

    if (grp == NULL) {
	return;
    }

    for (i=0; i<grp->nplots; i++) { 
	gretl_matrix_free(grp->plots[i].outliers);
	free(grp->plots[i].bool);
    }

    free(grp->title);
    free(grp->x);
    free(grp->plots);
    gretl_matrix_free(grp->dvals);

    free(grp);
}

static void boxplot_init (BOXPLOT *bp)
{
    bp->mean = NADBL;
    bp->median = NADBL;
    bp->conf[0] = bp->conf[1] = NADBL;
    bp->uq = bp->lq = NADBL;
    bp->max = bp->min = NADBL;
    bp->wmax = bp->wmin = NADBL;
    bp->n = 0;
    bp->varname[0] = '\0';
    bp->bool = NULL;
    bp->outliers = NULL;
}

static int plotgroup_factor_setup (PLOTGROUP *grp, 
				   const DATASET *dset)
{
    int T = sample_size(dset);
    int dv = grp->list[2];
    int err = 0;

    grp->dvals = gretl_matrix_values(dset->Z[dv] + dset->t1, T, 
				     OPT_S, &err);
    if (!err) {
	grp->nplots = gretl_vector_get_length(grp->dvals);
    }

    return err;
}

static PLOTGROUP *plotgroup_new (const int *list, 
				 const char *literal,
				 const DATASET *dset,
				 double limit,
				 gretlopt opt)
{
    PLOTGROUP *grp;
    int i, n, err = 0;

    grp = malloc(sizeof *grp);
    if (grp == NULL) {
	return NULL;
    }

    grp->nplots = 0;
    grp->list = NULL;
    grp->plots = NULL;
    grp->x = NULL;
    grp->dvals = NULL;
    grp->title = NULL;

    grp->literal = literal;

    if (na(limit)) {
	grp->limit = 1.5;
    } else {
	grp->limit = limit;
    }

    *grp->xlabel = '\0';
    *grp->ylabel = '\0';
    
    if (dset != NULL) {
	grp->list = list;
	n = sample_size(dset);
	grp->x = malloc(n * sizeof *grp->x);
	if (grp->x == NULL) {
	    err = E_ALLOC;
	}	
    }

    if (!err) {
	if (opt & OPT_Z) {
	    if (dset != NULL) {
		strcpy(grp->ylabel, dset->varname[list[1]]);
		strcpy(grp->xlabel, dset->varname[list[2]]);
		err = plotgroup_factor_setup(grp, dset);
	    }
	} else {
	    grp->nplots = list[0];
	}
	if ((opt & OPT_P) && dset != NULL) {
	    strcpy(grp->ylabel, series_get_label(dset, 1));
	    strcpy(grp->xlabel, _("group"));
	}
    }

    if (!err) {
	grp->plots = malloc(grp->nplots * sizeof *grp->plots);
	if (grp->plots == NULL) {
	    err = E_ALLOC;
	    grp->nplots = 0;
	}
    }

    if (!err) {
	for (i=0; i<grp->nplots; i++) {
	    boxplot_init(&grp->plots[i]);
	}
	grp->flags = 0;
	grp->n_bools = 0;
	grp->gmax = grp->gmin = 0.0;
    }

    if (err) {
	plotgroup_destroy(grp);
	grp = NULL;
    }

    return grp;
}

static int test_for_outliers (PLOTGROUP *grp)
{
    int i, n = 0;

    for (i=0; i<grp->nplots; i++) {
	if (grp->plots[i].outliers != NULL) {
	    n += gretl_vector_get_length(grp->plots[i].outliers);
	}
    }

    return n;
}

static int test_for_bool (PLOTGROUP *grp)
{
    int i;

    for (i=0; i<grp->nplots; i++) {
	if (grp->plots[i].bool != NULL) {
	    return 1;
	}
    }

    return 0;
}

static void get_box_y_range (PLOTGROUP *grp, double gyrange,
			     double *ymin, double *ymax)
{
    if (grp->gmin >= 0 && grp->gmin < 0.05 &&
	grp->gmax <= 1 && grp->gmax > 0.95) {
	/* proportions? */
	*ymin = 0.0;
	*ymax = 1.0;
    } else {
	*ymin = grp->gmin - 0.04 * gyrange;
	*ymax = grp->gmax + 0.04 * gyrange;
    }
}

static int write_gnuplot_boxplot (PLOTGROUP *grp, 
				  gretlopt opt, 
				  const char **labels)
{
    FILE *fp = NULL;
    BOXPLOT *bp;
    int n = grp->nplots;
    double lpos, h, gyrange;
    double ymin, ymax;
    int anybool, n_outliers;
    int lwidth = 2;
    int fmt, qtype = 2;
    int i, err = 0;

    fp = open_plot_input_file(PLOT_BOXPLOTS, &err);

    if (err) {
	return err;
    }

    fmt = specified_gp_output_format();
    if (fmt == GP_TERM_EPS) {
	/* line type for main outline, monochrome */
	qtype = 1;
    }    

    anybool = test_for_bool(grp);
    n_outliers = test_for_outliers(grp);

    gyrange = grp->gmax - grp->gmin;
    h = gyrange / 20;
    lpos = grp->gmin - h;
    get_box_y_range(grp, gyrange, &ymin, &ymax);

    if (opt & (OPT_Z | OPT_P)) {
	if (n > 1 && n <= 30) {
	    set_use_xtics(grp);
	}
    } 

    if (grp->title != NULL) {
	fprintf(fp, "set title \"%s\"\n", grp->title);
    } 

    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xrange [0:%d]\n", n + 1);
    fprintf(fp, "set yrange [%g:%g]\n", ymin, ymax);

    if (*grp->ylabel != '\0') {
	fprintf(fp, "set ylabel \"%s\"\n", grp->ylabel);
    }

    if (*grp->xlabel != '\0') {
	fprintf(fp, "set xlabel \"%s\"\n", grp->xlabel);
    }    

    if (n > 30) {
	lwidth = 1;
    } else {
	fputs("set ytics nomirror\n", fp);
	fputs("set border 2\n", fp);

	if (use_xtics(grp)) {
	    fputs("set xtics (", fp);
	    for (i=0; i<n; i++) {
		if (opt & OPT_Z) {
		    if (labels==NULL) {
			fprintf(fp, "\"%g\" %d", grp->dvals->val[i], i+1);
		    } else {
			fprintf(fp, "\"%s\" %d", labels[i], i+1);
		    }
		} else {
		    fprintf(fp, "\"%s\" %d", grp->plots[i].varname, i+1);
		}
		if (i < n - 1) {
		    fputs(", ", fp);
		} else {
		    fputs(") scale 0.0\n", fp);
		} 
	    }
	    fputs("set xtics nomirror\n", fp);
	} else {
	    fputs("set noxtics\n", fp);
	    if (n > 1 || anybool) {
		fprintf(fp, "set bmargin %d\n", (anybool)? 4 : 3);
		for (i=0; i<n; i++) {
		    fprintf(fp, "set label \"%s\" at %d,%g center\n", 
			    grp->plots[i].varname, i+1, lpos);
		    if (grp->plots[i].bool != NULL) {
			/* FIXME positioning? */
			fprintf(fp, "set label \"%s\" at %d,%g center\n", 
				grp->plots[i].bool, i+1, lpos - .8 * h);
		    }
		}
	    }
	}
    }

    fputs("set boxwidth 0.3 relative\n", fp);

    if (grp->literal != NULL) {
	print_gnuplot_literal_lines(grp->literal, fp);
    }

    fputs("plot \\\n", fp);
    /* the quartiles and extrema */
    fprintf(fp, "'-' using 1:3:2:5:4 w candlesticks lt %d lw %d "
	    "notitle%s, \\\n", qtype, lwidth, (whiskerbars(grp))? 
	    " whiskerbars 0.5" : "");
    /* the median */
    fputs("'-' using 1:2:2:2:2 w candlesticks lt -1 notitle", fp);

    if (show_mean(grp) || do_intervals(grp) || n_outliers > 0) {
	/* continue plot lines array */
	fputs(", \\\n", fp);
    } else {
	fputc('\n', fp);
    }

    if (show_mean(grp)) {
	/* plot the mean as point */
	fputs("'-' using 1:2 w points pt 1 notitle", fp);
    } else if (do_intervals(grp)) {
	/* upper and lower bounds of median interval */
	fputs("'-' using 1:2:2:2:2 w candlesticks lt 0 notitle, \\\n", fp);
	fputs("'-' using 1:2:2:2:2 w candlesticks lt 0 notitle", fp);
    }

    if (n_outliers > 0) {
	fputs(", \\\n", fp);
	fprintf(fp, "'-' using 1:2 w points lt 0 pt 7 ps 0.25 notitle\n");
    } else {
	fputc('\n', fp);
    }

    /* data for quartiles and whiskers */
    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	if (na(bp->wmin)) bp->wmin = bp->min;
	if (na(bp->wmax)) bp->wmax = bp->max;
	fprintf(fp, "%d %g %g %g %g\n", i+1, bp->wmin, bp->lq, bp->uq, bp->wmax);
    }
    fputs("e\n", fp);

    /* data for median */
    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	fprintf(fp, "%d %g\n", i+1, bp->median);
    }
    fputs("e\n", fp);

    if (show_mean(grp)) {
	/* data for mean */
	for (i=0; i<n; i++) {
	    bp = &grp->plots[i];
	    fprintf(fp, "%d %g\n", i+1, bp->mean);
	}
	fputs("e\n", fp);
    } else if (do_intervals(grp)) {
	/* data for lower median bound */
	for (i=0; i<n; i++) {
	    bp = &grp->plots[i];
	    fprintf(fp, "%d %g\n", i+1, bp->conf[0]);
	}
	fputs("e\n", fp);
	/* data for upper median bound */
	for (i=0; i<n; i++) {
	    bp = &grp->plots[i];
	    fprintf(fp, "%d %g\n", i+1, bp->conf[1]);
	}
	fputs("e\n", fp);
    }

    if (n_outliers > 0) {
	int j, nout;

	fprintf(fp, "# auxdata %d 2\n", n_outliers);
	for (i=0; i<n; i++) {
	    bp = &grp->plots[i];
	    if (bp->outliers != NULL) {
		nout = gretl_vector_get_length(bp->outliers);
		for (j=0; j<nout; j++) {
		    fprintf(fp, "%d %.10g\n", i+1, 
			    gretl_vector_get(bp->outliers, j));
		}
	    }
	}
	fputs("e\n", fp);    
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static int transcribe_array_factorized (PLOTGROUP *grp, int i,
					const DATASET *dset) 
{
    const double *y = dset->Z[grp->list[1]];
    const double *d = dset->Z[grp->list[2]];
    int t, n = 0;

    for (t=dset->t1; t<=dset->t2; t++) {
	if (d[t] == grp->dvals->val[i] && !na(y[t])) {
	    grp->x[n++] = y[t];
	}
    }

    return n;
}

static int factorized_boxplot_check (const int *list, 
				     char **bools, 
				     int *haslabels,
				     const DATASET *dset)
{
    int err = 0;

    if (list[0] != 2 || bools != NULL) {
	err = E_DATA;
    } else {
	int v2 = list[2];
	
	if (!series_is_discrete(dset, v2) &&
	    !gretl_isdiscrete(dset->t1, dset->t2, dset->Z[v2])) {
	    gretl_errmsg_set(_("You must supply two variables, the second of "
			       "which is discrete"));
	    err = E_DATA;
	}

	*haslabels = is_string_valued(dset, v2);
    }

    return err;
}

static int all_missing (int t1, int t2, const double *x)
{
    int t;

    for (t=t1; t<=t2; t++) {
	if (!na(x[t])) {
	    return 0;
	}
    }

    return 1;
}

/* build the boxplot group data structure, then pass it to 
   gnuplot_do_boxplot to write the command file */

static int real_boxplots (const int *list, char **bools, 
			  const char *literal,
			  const DATASET *dset,
			  gretlopt opt)
{
    PLOTGROUP *grp;
    double lim = NADBL;
    int i, k, np;
    int err = 0;
    int haslabels;
    const char **labels = NULL;

    if (opt & OPT_L) {
	/* we should have an explicit outlier limit */
	lim = get_optval_double(BXPLOT, OPT_L);
	if (na(lim) || lim < 0) {
	    return E_BADOPT;
	} 
    }

    for (i=1; i<=list[0]; i++) {
	int v = list[i];

	if (gretl_isconst(dset->t1, dset->t2, dset->Z[v])) {
	    gretl_errmsg_sprintf(_("%s is a constant"), dset->varname[v]);
	    return E_DATA;
	} else if (all_missing(dset->t1, dset->t2, dset->Z[v])) {
	    return E_MISSDATA;
	}
    }

    if (opt & OPT_Z) {
	err = factorized_boxplot_check(list, bools, &haslabels, dset);
	if (err) {
	    return err;
	} 

	if (haslabels) {
	    int n_labels;
	    labels = series_get_string_vals(dset, list[2], &n_labels);
	}
    } 

    grp = plotgroup_new(list, literal, dset, lim, opt);
    if (grp == NULL) {
	return E_ALLOC;
    }

    if (*grp->ylabel != '\0' && *grp->xlabel != '\0') {
	grp->title = gretl_strdup_printf(_("Distribution of %s by %s"),
					 grp->ylabel, grp->xlabel);
    } else if (list[0] == 1 && bools == NULL) {
	/* doing a single, standard plot */
	const char *s = series_get_label(dset, list[1]);
	char yname[VNAMELEN];

	if (s != NULL && strncmp(s, "residual for ", 13) == 0 &&
	    gretl_scan_varname(s + 13, yname) == 1) {
	    grp->title = gretl_strdup_printf(_("Regression residuals (= observed - fitted %s)"), 
					     yname);
	} else {
	    grp->title = gretl_strdup(dset->varname[list[1]]);
	}
    }

    if (opt & OPT_O) {
	set_do_intervals(grp);
    }

    np = grp->nplots;
    k = 0;

    for (i=0; i<np && !err; i++) {
	BOXPLOT *plot;
	const char *vname;
	int n;

	if (opt & OPT_Z) {
	    n = transcribe_array_factorized(grp, i, dset);
	    vname = dset->varname[list[1]];
	} else {
	    n = transcribe_array(grp->x, dset->Z[list[i+1]], dset);
	    vname = dset->varname[list[i+1]];
	}

	if (n < 2) {
	    if (grp->nplots == 1) {
		/* we'll be out of data */
		err = E_DATA;
		break;
	    } else {
		/* discard this plot */
		grp->nplots -= 1;
		continue;
	    }
	}

	plot = &grp->plots[k]; /* note: k rather than i */
	plot->mean = gretl_mean(0, n - 1, grp->x);
	qsort(grp->x, n, sizeof *grp->x, gretl_compare_doubles);
	quartiles_etc(grp->x, n, plot, grp->limit);
	
	plot->n = n;

	if (k == 0 || plot->min < grp->gmin) {
	    grp->gmin = plot->min;
	}

	if (k == 0 || plot->max > grp->gmax) {
	    grp->gmax = plot->max;
	}

	if (do_intervals(grp)) {
	    if (median_interval(grp->x, n, &plot->conf[0], &plot->conf[1])) {
		plot->conf[0] = plot->conf[1] = NADBL;
		unset_do_intervals(grp);
	    }
	} 

	strcpy(plot->varname, vname);

	if (opt & OPT_Z) {
	    plot->bool = gretl_strdup_printf("%g", grp->dvals->val[i]);
	    grp->n_bools += 1;
	} else if (bools != NULL && bools[i] != NULL) {
	    plot->bool = bools[i];
	    grp->n_bools += 1;
	} 

	k++; /* increment the actually used plot index */
    }

    if (!err) {
	if (!do_intervals(grp)) {
	    set_show_mean(grp);
	}
	err = write_gnuplot_boxplot(grp, opt, labels);
    }

    plotgroup_destroy(grp);   
    
    return err;
}

/* boxplots using a single selected series, by panel group */

static int panel_group_boxplots (int vnum, const DATASET *dset,
				 gretlopt opt)
{
    DATASET *gdset;
    int nunits, T = dset->pd;
    int *list = NULL;
    int nv, u0, i, t, s, s0;
    int err = 0;

    nunits = panel_sample_size(dset);
    nv = nunits + 1;
    u0 = dset->t1 / T;

    gdset = create_auxiliary_dataset(nv, T, 0);
    if (gdset == NULL) {
	return E_ALLOC;
    }

    list = gretl_consecutive_list_new(1, nv - 1);
    if (list == NULL) {
	destroy_dataset(gdset);
	return E_ALLOC;
    }

    s0 = dset->t1 * dset->pd;

    /* record the name of the original variable */
    series_set_label(gdset, 1, dset->varname[vnum]);

    for (i=0; i<nunits; i++) {
	sprintf(gdset->varname[i+1], "%d", u0+i+1);
	s = s0 + i * T;
	for (t=0; t<T; t++) {
	    gdset->Z[i+1][t] = dset->Z[vnum][s++];
	}
    }

    err = real_boxplots(list, NULL, NULL, gdset, opt);

    destroy_dataset(gdset);
    free(list);

    return err;
}

/*
 * boxplots:
 * @list: list of variables to plot.
 * @dset: dataset struct.
 * @opt: may include OPT_P for a plot by panel group (in which
 * case @list should have just one member), or OPT_Z for a
 * "factorized" plot (@list should have exactly three members),
 * and/or OPT_L to set the limit for outliers.
 *
 * Constructs a boxplot command file for use with gnuplot.
 *
 * Returns: 0 on success, error code on error.
 */

int boxplots (const int *list, const char *literal, 
	      const DATASET *dset, gretlopt opt)
{
    int err;

    if (opt & OPT_P) {
	if (!multi_unit_panel_sample(dset) ||
	    list[0] > 1 || (opt & OPT_Z)) {
	    err = E_BADOPT;
	} else {
	    err = panel_group_boxplots(list[1], dset, opt);
	}
    } else {
	err = real_boxplots(list, NULL, literal, dset, opt);
    }

    return err;
}

/* remove extra spaces around operators in boxplots line */

static char *boxplots_fix_parentheses (const char *line, int *err)
{
    char *s, *p;
    int inparen = 0;
    int i, flen, len = strlen(line);

    /* make room to insert space before parens, if needed */
    flen = len;
    for (i=0; i<len; i++) {
	if (i > 0 && line[i] == '(' && line[i-1] != ' ') {
	    flen++;
	}
    }

    s = malloc(flen + 1);
    if (s == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    p = s;
    for (i=0; i<len; i++) {
	if (line[i] == '(') {
	    if (i > 0 && line[i-1] != ' ') {
		*p++ = ' ';
	    }
	    inparen = 1;
	}
	if (line[i] == ')') {
	    if (inparen == 1) {
		inparen = 0;
	    } else {
		*err = E_PARSE;
		free(s);
		return NULL;
	    }
	}
	if (inparen && line[i] == ' ') ;
	else *p++ = line[i];
    }

    *p = 0;

    return s;
}

int boolean_boxplots (const char *line, const char *literal,
		      DATASET *dset, gretlopt opt)
{
    int i, k, v, nvars, nbool;
    int n = dset->n;
    int origv = dset->v;
    char *tok, *s = NULL;
    char **bools = NULL;
    int *list = NULL;
    int err = 0;

    if (!strncmp(line, "boxplots ", 9)) {
	line += 9;
    } else if (!strncmp(line, "boxplot ", 8)) {
	line += 8;
    }

    s = boxplots_fix_parentheses(line, &err);

    if (!err) {
	nvars = special_varcount(s);
	if (nvars == 0) {
	    err = E_ARGS;
	}
    }

    if (!err) {
	list = gretl_list_new(nvars);
	bools = strings_array_new(nvars);
	if (list == NULL || bools == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	free(list);
	free(bools);
	free(s);
	return err;
    }

    i = nbool = 0;

    /* record the command string */
    *boxplots_string = '\0';
    strncat(boxplots_string, s, BPSTRLEN - 1);

    while (!err && (tok = strtok((i)? NULL : s, " "))) {
	if (tok[0] == '(') {
	    if (i) {
		bools[i-1] = gretl_strdup(tok);
		nbool++;
	    } else {
		err = E_DATA;
	    }
	} else {
	    if (isdigit(tok[0])) { 
		v = atoi(tok);
		if (v < origv) {
		    list[++i] = v;
		} else {
		    gretl_errmsg_sprintf(_("got invalid variable number %d"), v);
		    err = E_DATA;
		}
	    } else if (isalpha(tok[0])) {
		v = current_series_index(dset, tok);
		if (v >= 0) {
		    list[++i] = v;
		} else {
		    gretl_errmsg_sprintf(_("got invalid varname '%s'"), tok);
		    err = E_DATA;
		}
	    } else {
		gretl_errmsg_sprintf(_("got invalid field '%s'"), tok);
		err = E_DATA; 
	    }
	}
    }

    /* Now we add nbool new variables, with ID numbers origv,
       origv + 1, and so on.  These are the original variables
       that have boolean conditions attached, masked by those
       conditions. 
    */

    k = origv;
    nbool = 0;

    for (i=0; i<nvars && !err; i++) {
	char formula[80];
	
	if (bools[i] == NULL) {
	    continue;
	}

	sprintf(formula, "bool_%d = %s", i, bools[i]);
	err = generate(formula, dset, OPT_P, NULL);

	if (!err) {
	    int t, vi = list[i+1];

	    for (t=0; t<n; t++) {
		if (dset->Z[k][t] == 1.0) {
		    dset->Z[k][t] = dset->Z[vi][t];
		} else { 
		    dset->Z[k][t] = NADBL;
		}
	    }
	    strcpy(dset->varname[k], dset->varname[vi]);
	    list[i+1] = k++;
	    nbool++;
	}
    }

    if (!err) {
	err = real_boxplots(list, bools, literal, dset, opt);
    } 
    
    free(list);
    free(bools); /* the bool[i]s are now attached to the plots */
    free(s);

    if (nbool) {
	dataset_drop_last_variables(dset, nbool);
    }
    
    return err;
}

/* parse "factor" values out of the xtics line in a gretl
   boxplot file, e.g. ("0" 1, "1" 2, "4", 3) */

static int get_vals_from_tics (PLOTGROUP *grp, const char *s)
{
    double x;
    int i, err = 0;

    grp->dvals = gretl_column_vector_alloc(grp->nplots);
    if (grp->dvals == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<grp->nplots && !err; i++) {
	s = strchr(s, '"');
	if (s == NULL) {
	    err = E_DATA;
	} else if (sscanf(s+1, "%lf", &x) != 1) {
	    err = E_DATA;
	} else {
	    gretl_vector_set(grp->dvals, i, x);
	    s++; /* skip the current '"' */
	    s += strcspn(s, "\""); /* skip to closing */
	    s++; 
	}
    }

    if (err) {
	gretl_matrix_free(grp->dvals);
	grp->dvals = NULL;
    }

    return err;
}

/* read outliers data out of a boxplots file and revise the
   min and max values for the individual plots */

static int read_plot_outliers (PLOTGROUP *grp, char *line,
			       size_t lsize, FILE *fp)
{
    BOXPLOT *plot;
    double x;
    int i, k, n;
    int err = 0;

    if (sscanf(line + 9, "%d", &n) != 1) {
	return E_DATA;
    }

    for (i=0; i<n && !err; i++) {
	if (fgets(line, lsize, fp) == NULL) {
	    err = E_DATA;
	} else if (sscanf(line, "%d %lf", &k, &x) != 2) {
	    err = E_DATA;
	} else if (k < 1 || k > grp->nplots) {
	    err = E_DATA;
	} else {
	    plot = &grp->plots[k-1];
	    if (x < plot->min) {
		plot->min = x;
	    } else if (x > plot->max) {
		plot->max = x;
	    }
	}
    }

    return err;
}

/* read data out of a boxplots file in order to reconstruct
   numerical summary */

static int read_plot_data_block (PLOTGROUP *grp, char *line, 
				 size_t lsize, int *block, 
				 FILE *fp)
{
    BOXPLOT *plot;
    int i, j, nf, targ;
    int err = 0;

    for (i=0; i<grp->nplots && !err; i++) {
	plot = &grp->plots[i];
	targ = 2; /* the number of values we should read */
	nf = 0;   /* the number of values actually read */
	if (*block == 0) {
	    /* min, Q1, Q3, max */
	    targ = 5;
	    nf = sscanf(line, "%d %lf %lf %lf %lf", 
			&j, &plot->min, &plot->lq, &plot->uq, &plot->max);
	} else if (*block == 1) {
	    /* median */
	    nf = sscanf(line, "%d %lf", &j, &plot->median);
	} else if (*block == 2) {
	    /* mean, or lower bound of c.i. */
	    if (do_intervals(grp)) {
		nf = sscanf(line, "%d %lf", &j, &plot->conf[0]);
	    } else {
		nf = sscanf(line, "%d %lf", &j, &plot->mean);
	    }
	} else if (*block == 3) {
	    /* upper bound of c.i. */
	    nf = sscanf(line, "%d %lf", &j, &plot->conf[1]);
	}
	if (nf != targ || j != i + 1) {
	    /* data read screwed up */
	    err = E_DATA;
	} else if (fgets(line, lsize, fp) == NULL) {
	    /* can't get data line for next plot */
	    err = E_DATA;
	}
    }

    *block += 1;

    return err;
}

static void title_to_varname (PLOTGROUP *grp, const char *title)
{
    int len = strlen(title);

    if (len < 16 && len == gretl_namechar_spn(title)) {
	strcpy(grp->plots[0].varname, title);
    } 
}

/* Given a gnuplot boxplot file created by gretl, parse out the
   "numerical summary" information.  Note that this requires a
   strictly regimented plot file, so if you make any changes to the
   way boxplot files are printed you should check this function for 
   breakage.
*/

int boxplot_numerical_summary (const char *fname, PRN *prn)
{
    PLOTGROUP *grp = NULL;
    int factorized = 0;
    int plotlines = 0;
    char yvar[VNAMELEN];
    char xvar[VNAMELEN];
    char title[128];
    char fmt[16], tfmt[16], line[512];
    int parsing = 1;
    FILE *fp;
    int n = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    *yvar = *xvar = *title = '\0';
    sprintf(fmt, "%%%d[^\\\"]", VNAMELEN-1);
    sprintf(tfmt, "%%%d[^\\\"]", 127);

    /* first pass: count the plots, using labels as heuristic */

    while (fgets(line, sizeof line, fp)) {
	if (parsing && !strncmp(line, "# start literal", 15)) {
	    parsing = 0;
	} else if (!parsing && !strncmp(line, "# end literal", 13)) {
	    parsing = 1;
	}
	if (!parsing) {
	    continue;
	}
	if (!strncmp(line, "1 ", 2)) {
	    /* reached first data block */
	    break;
	} else if (!strncmp(line, "set ylabel ", 11)) {
	    sscanf(line + 12, fmt, yvar);
	} else if (!strncmp(line, "set xlabel ", 11)) {
	    sscanf(line + 12, fmt, xvar);
	} else if (!strncmp(line, "set title ", 10)) {
	    sscanf(line + 11, tfmt, title);
	} else if (!strncmp(line, "set label ", 10)) {
	    if (!strncmp(line + 10, "\"(", 2)) {
		; /* boolean specifier */
	    } else {
		/* should correspond to a plot */
		n++;
	    }
	} else if (!strncmp(line, "'-'", 3)) {
	    plotlines++;
	}
    }

    if (n == 0 && *yvar != '\0' && *xvar != '\0') {
	/* try getting a count of plots via the first data block */
	while (fgets(line, sizeof line, fp)) {
	    sscanf(line, "%d", &n);
	    if (*line == 'e') {
		/* end of block */
		break;
	    }
	}
	if (n > 0) {
	    factorized = 1;
	}
    } else if (n == 0 && *title != '\0') {
	/* plotting a single variable, no labels */
	n = 1;
    }

    if (n == 0) {
	err = E_DATA;
    } else {
	int list[2] = {n, 0};

	grp = plotgroup_new(list, NULL, NULL, NADBL, OPT_NONE);
	if (grp == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && plotlines > 3) {
	/* note: this is provisional */
	set_do_intervals(grp);
    }

    if (!err) {
	/* second pass: read the data */
	int i = 0, block = 0, outliers = 0;
	char *label;

	rewind(fp);
	gretl_push_c_numeric_locale();
	parsing = 1;

	while (fgets(line, sizeof line, fp) && !err) {
	    if (parsing && !strncmp(line, "# start literal", 15)) {
		parsing = 0;
	    } else if (!parsing && !strncmp(line, "# end literal", 13)) {
		parsing = 1;
	    }
	    if (!parsing) {
		continue;
	    }
	    if (factorized && !strncmp(line, "set xtics (", 11)) {
		err = get_vals_from_tics(grp, line + 10);
	    } else if (!factorized && !strncmp(line, "set label ", 10)) {
		label = gretl_quoted_string_strdup(line + 10, NULL);
		if (label == NULL) {
		    err = E_DATA;
		} else if (*label == '(') {
		    /* should be boolean specifier, following varname */
		    if (i > 0) {
			grp->plots[i-1].bool = label;
		    } else {
			free(label);
			err = E_DATA;
		    }
		} else {
		    strncat(grp->plots[i++].varname, label, VNAMELEN-1);
		    free(label);
		}
	    } else if (!strncmp(line, "# auxdata", 9)) {
		outliers = 1;
		err = read_plot_outliers(grp, line, sizeof line, fp);
	    } else if (!strncmp(line, "1 ", 2)) {
		/* start of a data block */
		err = read_plot_data_block(grp, line, sizeof line, &block, fp);
	    } 
	}

	gretl_pop_c_numeric_locale();

	if (!err && outliers && plotlines == 4) {
	    /* we didn't really get median intervals after all */
	    for (i=0; i<grp->nplots; i++) {
		grp->plots[i].mean = grp->plots[i].conf[0];
		grp->plots[i].conf[0] = NADBL;
	    }
	    unset_do_intervals(grp);
	}
    }

    fclose(fp);

    if (!err) {
	if (n == 1 && grp->plots[0].varname[0] == '\0' &&
	    *title != '\0') {
	    title_to_varname(grp, title);
	}
	strcpy(grp->ylabel, yvar);
	strcpy(grp->xlabel, xvar);
	if (factorized) {
	    set_box_factor(grp);
	}
	if (!na(grp->plots[0].conf[0])) {
	    set_do_intervals(grp);
	}
	boxplot_print_stats(grp, prn);
    }

    plotgroup_destroy(grp);

    return err;
}
