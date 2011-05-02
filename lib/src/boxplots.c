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
    double conf[2];
    double uq, lq;
    double max, min;
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
    char *numbers;
    BOXPLOT *plots;
    double gmin, gmax;
    double limit;
    double *x;
    gretl_matrix *dvals;
    char xlabel[VNAMELEN];
    char ylabel[VNAMELEN];
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

    gretl_array_quantiles(x, n, p, 3);

    plt->lq = p[0];
    plt->median = p[1];
    plt->uq = p[2];

    if (limit > 0) {
	double d = limit * (plt->uq - plt->lq);
	double xlo = plt->lq - d;
	double xhi = plt->uq + d;
	int i, nout = 0;

	for (i=0; i<n; i++) {
	    if (x[i] < xlo || x[i] > xhi) {
		nout++;
	    }
	}

	if (nout > 0) {
	    plt->outliers = gretl_vector_alloc(nout);
	}

	if (plt->outliers != NULL) {
	    int j = 0;

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
    free(grp->numbers);
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
    bp->n = 0;
    bp->varname[0] = '\0';
    bp->bool = NULL;
    bp->outliers = NULL;
}

static int plotgroup_factor_setup (PLOTGROUP *grp, 
				   const double **Z,
				   const DATAINFO *pdinfo)
{
    int T = sample_size(pdinfo);
    int dv = grp->list[2];
    int err = 0;

    grp->dvals = gretl_matrix_values(Z[dv] + pdinfo->t1, T, &err);
    if (!err) {
	grp->nplots = gretl_vector_get_length(grp->dvals);
    }

    return err;
}

static PLOTGROUP *plotgroup_new (const int *list, 
				 const double **Z,
				 const DATAINFO *pdinfo,
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
    grp->numbers = NULL;

    if (na(limit)) {
	grp->limit = 1.5;
    } else {
	grp->limit = limit;
    }

    *grp->xlabel = '\0';
    *grp->ylabel = '\0';
    
    if (pdinfo != NULL) {
	grp->list = list;
	n = sample_size(pdinfo);
	grp->x = malloc(n * sizeof *grp->x);
	if (grp->x == NULL) {
	    err = E_ALLOC;
	}	
    }

    if (!err) {
	if (opt & OPT_Z) {
	    if (pdinfo != NULL) {
		strcpy(grp->ylabel, pdinfo->varname[list[1]]);
		strcpy(grp->xlabel, pdinfo->varname[list[2]]);
		err = plotgroup_factor_setup(grp, Z, pdinfo);
	    }
	} else {
	    grp->nplots = list[0];
	}
	if ((opt & OPT_P) && pdinfo != NULL) {
	    strcpy(grp->ylabel, VARLABEL(pdinfo, 1));
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

/* note: @fname is non-NULL only when we're doing the
   backward-compatibility thing of reading an old-style
   gretl boxplot file and converting it to gnuplot.
*/

static int write_gnuplot_boxplot (PLOTGROUP *grp, const char *fname,
				  gretlopt opt)
{
    FILE *fp = NULL;
    BOXPLOT *bp;
    int n = grp->nplots;
    double lpos, h, gyrange;
    double ymin, ymax;
    int anybool, n_outliers;
    int lwidth = 2;
    int fmt = 0;
    int qtype = 2;
    int i, err = 0;

    if (fname != NULL) {
	/* converting old-style file */
	fp = gretl_fopen(fname, "w");
	if (fp == NULL) {
	    err = E_FOPEN;
	}
    } else if (gretl_in_batch_mode()) {
	/* batch mode: auto-named file */
	fp = get_gnuplot_batch_stream(PLOT_BOXPLOTS, &err);
	if (!err) {
	    fmt = specified_gp_output_format();
	    if (fmt == GP_TERM_EPS || fmt == GP_TERM_PLT) {
		/* line type for main outline, monochrome */
		qtype = 1;
	    }
	}	    
    } else {
	/* displaying graph: auto-named temp file */
	fp = get_plot_input_stream(PLOT_BOXPLOTS, &err);
    }

    if (err) {
	return err;
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
		    fprintf(fp, "\"%g\" %d", grp->dvals->val[i], i+1);
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
	double ymin, ymax;

	bp = &grp->plots[i];
	if (bp->outliers != NULL) {
	    ymin = bp->lq - grp->limit * (bp->uq - bp->lq);
	    ymax = bp->uq + grp->limit * (bp->uq - bp->lq);
	} else {
	    ymin = bp->min;
	    ymax = bp->max;
	}
	fprintf(fp, "%d %g %g %g %g\n", i+1, ymin, bp->lq, bp->uq, ymax);
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

    fclose(fp);

    return err;
}

static int gnuplot_do_boxplot (PLOTGROUP *grp, gretlopt opt)
{
    int err = write_gnuplot_boxplot(grp, NULL, opt);

    if (!err) {
	err = gnuplot_make_graph();
    } 

    return err;
}

int transcribe_array_factorized (PLOTGROUP *grp, int i,
				 const double **Z, 
				 const DATAINFO *pdinfo) 
{
    const double *y = Z[grp->list[1]];
    const double *d = Z[grp->list[2]];
    int t, n = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (d[t] == grp->dvals->val[i] && !na(y[t])) {
	    grp->x[n++] = y[t];
	}
    }

    return n;
}

static int factorized_boxplot_check (const int *list, 
				     char **bools, 
				     const double **Z,
				     const DATAINFO *pdinfo)
{
    int err = 0;

    if (list[0] != 2 || bools != NULL) {
	err = E_DATA;
    } else {
	int v2 = list[2];
	
	if (!var_is_discrete(pdinfo, v2) &&
	    !gretl_isdiscrete(pdinfo->t1, pdinfo->t2, Z[v2])) {
	    gretl_errmsg_set(_("You must supply two variables, the second of "
			       "which is discrete"));
	    err = E_DATA;
	}
    }

    return err;
}

static int real_boxplots (const int *list, char **bools, 
			  const double **Z, 
			  const DATAINFO *pdinfo,
			  gretlopt opt)
{
    PLOTGROUP *grp;
    double lim = NADBL;
    int i, k, np;
    int err = 0;

    if (opt & OPT_L) {
	/* we should have an explicit outlier limit */
	lim = get_optval_double(BXPLOT, OPT_L);
	if (na(lim) || lim < 0) {
	    return E_BADOPT;
	} 
    }

    if (opt & OPT_Z) {
	err = factorized_boxplot_check(list, bools, Z, pdinfo);
	if (err) {
	    return err;
	} 
    } 

    grp = plotgroup_new(list, Z, pdinfo, lim, opt);
    if (grp == NULL) {
	return E_ALLOC;
    }

    if (*grp->ylabel != '\0' && *grp->xlabel != '\0') {
	grp->title = gretl_strdup_printf(_("Distribution of %s by %s"),
					 grp->ylabel, grp->xlabel);
    } else if (list[0] == 1 && bools == NULL) {
	/* doing a single straight plot */
	grp->title = gretl_strdup(pdinfo->varname[list[1]]);
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
	    n = transcribe_array_factorized(grp, i, Z, pdinfo);
	    vname = pdinfo->varname[list[1]];
	} else {
	    n = transcribe_array(grp->x, Z[list[i+1]], pdinfo);
	    vname = pdinfo->varname[list[i+1]];
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
	plot->min = grp->x[0];
	plot->max = grp->x[n-1];
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
	err = gnuplot_do_boxplot(grp, opt);
    }

    plotgroup_destroy(grp);   
    
    return err;
}

/* boxplots using a single selected series, by panel group */

static int panel_group_boxplots (int vnum, const double **Z, 
				 const DATAINFO *pdinfo,
				 gretlopt opt)
{
    double **gZ = NULL;
    DATAINFO *gdinfo;
    int nunits, T = pdinfo->pd;
    int *list = NULL;
    int nv, u0, i, t, s, s0;
    int err = 0;

    nunits = panel_sample_size(pdinfo);
    nv = nunits + 1;
    u0 = pdinfo->t1 / T;

    gdinfo = create_auxiliary_dataset(&gZ, nv, T);
    if (gdinfo == NULL) {
	return E_ALLOC;
    }

    list = gretl_consecutive_list_new(1, nv - 1);
    if (list == NULL) {
	destroy_dataset(gZ, gdinfo);
	return E_ALLOC;
    }

    s0 = pdinfo->t1 * pdinfo->pd;

    /* record the name of the original variable */
    strcpy(VARLABEL(gdinfo, 1), pdinfo->varname[vnum]);

    for (i=0; i<nunits; i++) {
	sprintf(gdinfo->varname[i+1], "%d", u0+i+1);
	s = s0 + i * T;
	for (t=0; t<T; t++) {
	    gZ[i+1][t] = Z[vnum][s++];
	}
    }

    err = real_boxplots(list, NULL, (const double **) gZ, gdinfo, 
			opt);

    destroy_dataset(gZ, gdinfo);
    free(list);

    return err;
}

int boxplots (const int *list, const double **Z, const DATAINFO *pdinfo, 
	      gretlopt opt)
{
    int err;

    if (opt & OPT_P) {
	if (!multi_unit_panel_sample(pdinfo) ||
	    list[0] > 1 || (opt & OPT_Z)) {
	    err = E_BADOPT;
	} else {
	    err = panel_group_boxplots(list[1], Z, pdinfo, opt);
	}
    } else {
	err = real_boxplots(list, NULL, Z, pdinfo, opt);
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

int boolean_boxplots (const char *str, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt)
{
    int i, k, v, nvars, nbool;
    int n = pdinfo->n;
    int origv = pdinfo->v;
    char *tok, *s = NULL;
    char **bools = NULL;
    int *list = NULL;
    int err = 0;

    if (!strncmp(str, "boxplots ", 9)) {
	str += 9;
    } else if (!strncmp(str, "boxplot ", 8)) {
	str += 8;
    }

    s = boxplots_fix_parentheses(str, &err);

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
		v = current_series_index(pdinfo, tok);
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
	err = generate(formula, pZ, pdinfo, OPT_P, NULL);

	if (!err) {
	    int t, vi = list[i+1];

	    for (t=0; t<n; t++) {
		if ((*pZ)[k][t] == 1.0) {
		    (*pZ)[k][t] = (*pZ)[vi][t];
		} else { 
		    (*pZ)[k][t] = NADBL;
		}
	    }
	    strcpy(pdinfo->varname[k], pdinfo->varname[vi]);
	    list[i+1] = k++;
	    nbool++;
	}
    }

    if (!err) {
	err = real_boxplots(list, bools, (const double **) *pZ, 
			    pdinfo, opt);
    } 
    
    free(list);
    free(bools); /* the bool[i]s are now attached to the plots */
    free(s);

    if (nbool) {
	dataset_drop_last_variables(nbool, pZ, pdinfo);
    }
    
    return err;
}

static int check_confidence_interval (BOXPLOT *plt) 
{
    if (plt->conf[0] == plt->conf[1]) {
	plt->conf[0] = NADBL;
	plt->conf[1] = NADBL;
    }

    return !na(plt->conf[0]) && !na(plt->conf[1]);
}

static int plot_retrieve_outliers (BOXPLOT *plt, int n, FILE *fp)
{
    char line[80];
    double x;
    int i;

    plt->outliers = gretl_vector_alloc(n);
    if (plt->outliers == NULL) {
	return E_ALLOC;
    }

    if (fgets(line, sizeof line, fp) == NULL) {
	/* old rmax, rmin line */
	return E_DATA;
    }

    for (i=0; i<n; i++) {
	if (fgets(line, sizeof line, fp) == NULL ||
	    sscanf(line, "%*d vals %lf", &x) != 1) {
	    return E_DATA;
	} else {
	    gretl_vector_set(plt->outliers, i, x);
	}
    }

    return 0;
}

/* return non-zero if we get a mean */

static int maybe_get_plot_mean (BOXPLOT *plt, FILE *fp)
{
    char line[80];
    long pos = ftell(fp);
    double x = NADBL;
    int n, ok = 0;
    
    if (fgets(line, sizeof line, fp) == NULL) {
	return 0;
    } else if (sscanf(line, "%*d mean = %lf", &x) == 1) {
	plt->mean = x;
	ok = 1;
    } else {
	/* back up */
	fseek(fp, pos, SEEK_SET);
	return 0;
    }

    pos = ftell(fp);

    if (fgets(line, sizeof line, fp) == NULL) {
	return ok;
    } else if (sscanf(line, "%*d nobs = %d", &n) == 1) {
	plt->n = n;
    } else {
	/* back up */
	fseek(fp, pos, SEEK_SET);
    }

    return ok;
}

static void get_bool_from_line (const char *line, BOXPLOT *plt)
{
    char boolstr[128];

    if (sscanf(line, "%*d varname = %*s %127s", boolstr) == 1) {
	plt->bool = gretl_strdup(boolstr);
    }
}

/* Convert an old-style gretl boxplot file to a gnuplot command file:
   we use this when opening an old gretl session file that contains
   boxplots.
*/

int gnuplot_from_boxplot (const char *fname)
{
    FILE *fp;
    int i, j, nout;
    int gotmean, gotnotch;
    PLOTGROUP *grp = NULL;
    BOXPLOT *plt = NULL;
    char line[80], numstr[24];
    double gmax, gmin;
    int nplots = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", fname);
	return E_FOPEN;
    }

    for (i=0; i<6 && fgets(line, 79, fp) && !err; i++) {
	if (i == 1 && sscanf(line, "nplots = %d", &nplots) != 1) {
	    err = E_DATA;
	} else if (i == 2 && sscanf(line, "numbers = %7s", numstr) != 1) {
	    err = E_DATA;
	} else if (i == 5 && sscanf(line, "gmax gmin = %lf %lf", 
				    &gmax, &gmin) != 2) {
	    err = E_DATA;
	}
    }

    if (!err) {
	int list[2] = {nplots, 0};

	grp = plotgroup_new(list, NULL, NULL, NADBL, OPT_NONE);
	if (grp == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	goto bailout;
    }

    if (strcmp(numstr, "NULL")) {
	grp->numbers = gretl_strdup(numstr);
    }

    grp->gmax = gmax;
    grp->gmin = gmin;
    set_show_mean(grp);
    set_do_intervals(grp);

    for (i=0; i<grp->nplots && !err; i++) {
	plt = &grp->plots[i];
	nout = 0;

	for (j=0; j<7 && fgets(line, 79, fp) && !err; j++) {
	    if (j == 0 && 
		sscanf(line, "%*d median = %lf", &plt->median) != 1) {
		err = E_DATA;
	    } else if (j == 1 && 
		     sscanf(line, "%*d conf = %lf %lf", 
			    &plt->conf[0], & plt->conf[1]) != 2) {
		err = E_DATA;
	    } else if (j == 2 && 
		     sscanf(line, "%*d quartiles = %lf %lf", 
			    &plt->uq, &plt->lq) != 2) {
		err = E_DATA;
	    } else if (j == 3 && 
		     sscanf(line, "%*d maxmin = %lf %lf", 
			    &plt->max, &plt->min) != 2) {
		err = E_DATA;
	    } else if (j == 5) {
		if (sscanf(line, "%*d varname = %15s", plt->varname) != 1) {
		    err = E_DATA;
		} else {
		    get_bool_from_line(line, plt);
		}
	    } else if (j == 6 && 
		       sscanf(line, "%*d n_outliers = %d", &nout) != 1) {
		err = E_DATA;
	    }
	}

	if (err) {
	    break;
	}

	if (!err && do_intervals(grp)) {
	    gotnotch = check_confidence_interval(plt);
	    if (!gotnotch) {
		unset_do_intervals(grp);
	    }
	}	    

	/* any outliers? */
	if (nout > 0) {
	    err = plot_retrieve_outliers(plt, nout, fp);
	} 

	if (!err && show_mean(grp)) {
	    gotmean = maybe_get_plot_mean(plt, fp);
	    if (!gotmean) {
		unset_show_mean(grp);
	    }
	}
    }

 bailout:

    fclose(fp);

    if (!err) {
	if (show_mean(grp) && do_intervals(grp)) {
	    unset_show_mean(grp);
	}
	err = write_gnuplot_boxplot(grp, fname, OPT_NONE);
    } else if (err == E_DATA) {
	gretl_errmsg_set(_("boxplot file is corrupt"));
    }

    plotgroup_destroy(grp);

    return err;
}

/* e.g. ("0" 1, "1" 2, "4", 3) */

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

/* Given a gnuplot boxplot file created by gretl, parse out the
   "numerical summary" information.  Note that this requires a
   strictly regimented plot file, so if you make any changes to the
   way boxplot files are printed (see gnuplot_from_boxplot above) you
   should check this function for breakage.
*/

int boxplot_numerical_summary (const char *fname, PRN *prn)
{
    PLOTGROUP *grp = NULL;
    int factorized = 0;
    int plotlines = 0;
    char yvar[VNAMELEN];
    char xvar[VNAMELEN];
    char line[512];
    FILE *fp;
    int n = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    *yvar = *xvar = '\0';

    /* first pass: count the plots, using labels as heuristic */

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "1 ", 2)) {
	    /* reached first data block */
	    break;
	} else if (!strncmp(line, "set ylabel ", 11)) {
	    sscanf(line + 12, "%15[^\"]", yvar);
	} else if (!strncmp(line, "set xlabel ", 11)) {
	    sscanf(line + 12, "%15[^\"]", xvar);
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

	while (fgets(line, sizeof line, fp) && !err) {
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

    if (!err) {
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

    fclose(fp);

    plotgroup_destroy(grp);

    return err;
}
