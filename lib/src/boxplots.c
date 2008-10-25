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

typedef struct {
    int n;
    double *vals;
    double rmax;
    double rmin;
} OUTLIERS;

typedef struct {
    double mean;
    double median;
    double conf[2];
    double uq, lq;
    double max, min;
    int n;
    char varname[VNAMELEN];
    char *bool;
    OUTLIERS *outliers;
} BOXPLOT;

typedef struct {
    int nplots;
    int show_mean;
    int do_notches;
    int show_outliers;
    int n_bools;
    char *numbers;
    BOXPLOT *plots;
    double gmin, gmax;
} PLOTGROUP;

#define BPSTRLEN 128

static char boxplots_string[BPSTRLEN];

const char *get_last_boxplots_string (void)
{
    return boxplots_string;
}

static void quartiles (double *x, const int n, BOXPLOT *box)
{
    double p[3] = {0.25, 0.5, 0.75};

    gretl_array_quantiles(x, n, p, 3);

    box->lq = p[0];
    box->median = p[1];
    box->uq = p[2];
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
    int n = 0;
    char test[36];

    while (sscanf(s, "%35s", test) == 1) {
	if (*test != '(') n++;
	s += strspn(s, " ");
	s += strlen(test);
    }
    return n;
}

static void destroy_boxplots (PLOTGROUP *grp)
{
    int i;

    if (grp == NULL) {
	return;
    }

    for (i=0; i<grp->nplots; i++) { 
	free(grp->plots[i].outliers);
	free(grp->plots[i].bool);
    }

    free(grp->plots);
    free(grp->numbers);
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

static PLOTGROUP *plotgroup_new (int nplots)
{
    PLOTGROUP *grp;
    int i;

    grp = malloc(sizeof *grp);
    if (grp == NULL) {
	return NULL;
    }

    grp->plots = malloc(nplots * sizeof *grp->plots);
    if (grp->plots == NULL) {
	free(grp);
	grp = NULL;
    } else {
	for (i=0; i<nplots; i++) {
	    boxplot_init(&grp->plots[i]);
	}
	grp->nplots = nplots;
	grp->numbers = NULL;
	grp->show_mean = 0;
	grp->do_notches = 0;
	grp->show_outliers = 0;
	grp->n_bools = 0;
	grp->gmax = grp->gmin = 0.0;
    }

    return grp;
}

/* FIXME outliers */

static int write_gnuplot_boxplot (PLOTGROUP *grp)
{
    FILE *fp = NULL;
    BOXPLOT *bp;
    int n = grp->nplots;
    double loff;
    int i, err = 0;

    err = gnuplot_init(PLOT_BOXPLOTS, &fp);
    if (err) {
	return err;
    }

    fprintf(fp, "set xrange [0:%d]\n", n + 1);
    fputs("set ytics nomirror\n"
	  "set noxtics\n"
	  "set border 2\n"
	  "set bmargin 3\n", fp);

    loff = -(grp->gmax - grp->gmin) / 25;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	fprintf(fp, "set label \"%s\" at %d,%g center\n", 
		grp->plots[i].varname, i+1, loff);
	if (grp->plots[i].bool != NULL) {
	    /* FIXME positioning? */
	    fprintf(fp, "set label \"%s\" at %d,%g center\n", 
		    grp->plots[i].bool, i+1, loff * 2);
	}
    }

    fprintf(fp, "set boxwidth %g absolute\n", 1.0 / (n + 2));

    fputs("plot \\\n", fp);
    /* the quartiles and extrema */
    fputs("'-' using 1:3:2:5:4 w candlesticks lt 2 lw 2 "
	  "notitle whiskerbars 0.5, \\\n", fp);
    /* the median */
    fputs("'-' using 1:2:2:2:2 w candlesticks lt -1 notitle", fp);

    if (grp->show_mean || grp->do_notches) {
	/* continue plot lines array */
	fputs(", \\\n", fp);
    } else {
	fputc('\n', fp);
    }

    if (grp->show_mean) {
	/* plot the mean as point */
	fputs("'-' using 1:2 w points pt 1 notitle\n", fp);
    } else if (grp->do_notches) {
	/* upper and lower bounds of median interval */
	fputs("'-' using 1:2:2:2:2 w candlesticks lt 0 notitle, \\\n", fp);
	fputs("'-' using 1:2:2:2:2 w candlesticks lt 0 notitle\n", fp);
    }

    /* data for quartiles and extrema */
    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	fprintf(fp, "%d %g %g %g %g\n", i+1, bp->min, bp->lq, bp->uq, bp->max);
    }
    fputs("e\n", fp);

    /* data for median */
    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	fprintf(fp, "%d %g\n", i+1, bp->median);
    }
    fputs("e\n", fp);

    if (grp->show_mean) {
	/* data for mean */
	for (i=0; i<n; i++) {
	    bp = &grp->plots[i];
	    fprintf(fp, "%d %g\n", i+1, bp->mean);
	}
	fputs("e\n", fp);
    } else if (grp->do_notches) {
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

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

static int gnuplot_do_boxplot (PLOTGROUP *grp)
{
    int err;

    err = write_gnuplot_boxplot(grp);

    if (!err) {
	err = gnuplot_make_graph();
    }

    return err;
}

static int real_boxplots (int *list, char **bools, 
			  double ***pZ, const DATAINFO *pdinfo, 
			  gretlopt opt)
{
    int i, j, n = pdinfo->t2 - pdinfo->t1 + 1;
    double *x = NULL;
    PLOTGROUP *grp = NULL;
    int err = 0;

    x = malloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    grp = plotgroup_new(list[0]);
    if (grp == NULL) {
	free(x);
	return E_ALLOC;
    }

    if (opt & OPT_O) {
	grp->do_notches = 1;
    }

    for (i=0, j=0; i<grp->nplots; i++, j++) {
	BOXPLOT *plot = &grp->plots[i];
	int vi = list[i+1];

	n = ztox(vi, x, (const double **) *pZ, pdinfo);

	if (n < 2) {
	    gretl_list_delete_at_pos(list, i+1);
	    if (list[0] == 0) {
		err = E_DATA;
	    } else {
		grp->nplots -= 1;
		i--;
		continue;
	    }
	}

	if (err) {
	    break;
	}

	plot->mean = gretl_mean(0, n - 1, x);
	qsort(x, n, sizeof *x, gretl_compare_doubles);
	plot->min = x[0];
	plot->max = x[n-1];
	quartiles(x, n, &grp->plots[i]);
	plot->n = n;

	if (i == 0 || plot->min < grp->gmin) {
	    grp->gmin = plot->min;
	}

	if (i == 0 || plot->max > grp->gmax) {
	    grp->gmax = plot->max;
	}

	if (grp->do_notches) {
	    if (median_interval(x, n, &plot->conf[0], &plot->conf[1])) {
		plot->conf[0] = plot->conf[1] = NADBL;
		grp->do_notches = 0;
	    }
	} 

	strcpy(plot->varname, pdinfo->varname[vi]);

	if (bools != NULL) { 
	    plot->bool = bools[j];
	    if (plot->bool != NULL) {
		grp->n_bools += 1;
	    }
	} 
    }

    if (!err) {
	if (!grp->do_notches) {
	    grp->show_mean = 1;
	}
	err = gnuplot_do_boxplot(grp);
    }

    destroy_boxplots(grp);   
    free(x);
    
    return err;
}

int boxplots (int *list, double ***pZ, const DATAINFO *pdinfo, 
	      gretlopt opt)
{
    return real_boxplots(list, NULL, pZ, pdinfo, opt);
}

/* remove extra spaces around operators in boxplots line */

static char *boxplots_fix_parentheses (const char *line)
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

    s = boxplots_fix_parentheses(str);
    if (s == NULL) {
	return 1;
    }

    nvars = special_varcount(s);
    if (nvars == 0) {
	free(s);
	return 1;
    }

    list = gretl_list_new(nvars);
    bools = malloc(nvars * sizeof *bools);

    if (list == NULL || bools == NULL) {
	free(list);
	free(bools);
	free(s);
	return E_ALLOC;
    }

    for (i=0; i<nvars; i++) {
	bools[i] = NULL;
    }

    list[0] = nvars;
    i = 0;
    nbool = 0;

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
		v = series_index(pdinfo, tok);
		if (v < origv) {
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

    /* now we add nbool new variables, with ID numbers origv,
       origv + 1, and so on.  These are the original variables
       that have boolean conditions attached, masked by those
       conditions 
    */

    k = origv;
    nbool = 0;

    for (i=1; i<=list[0] && !err; i++) {
	char formula[80];
	int t;
	
	if (bools[i-1] == NULL) {
	    continue;
	}
	sprintf(formula, "bool_%d = %s", i-1, bools[i-1]);
	err = generate(formula, pZ, pdinfo, OPT_P, NULL);
	if (err) {
	    gretl_errmsg_sprintf(_("boxplots: generation of dummy variable failed\n%s"),
				 gretl_errmsg_get());
	    err = 1;
	} else {
	    for (t=0; t<n; t++) {
		if ((*pZ)[k][t] == 1.0) {
		    (*pZ)[k][t] = (*pZ)[list[i]][t];
		} else { 
		    (*pZ)[k][t] = NADBL;
		}
	    }
	    strcpy(pdinfo->varname[k], pdinfo->varname[list[i]]);
	    list[i] = k++;
	    nbool++;
	}
    }

    if (!err) {
	err = real_boxplots(list, bools, pZ, pdinfo, opt);
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
    int i;

    plt->outliers = malloc(sizeof *plt->outliers);
    if (plt->outliers == NULL) {
	return E_ALLOC;
    }

    plt->outliers->vals = malloc(n * sizeof *plt->outliers->vals);
    if (plt->outliers->vals == NULL) {
	return E_ALLOC;
    }

    plt->outliers->n = n;

    if (fgets(line, sizeof line, fp) == NULL) {
	return E_DATA;
    }

    if (sscanf(line, "%*d rmax rmin = %lf %lf", &plt->outliers->rmax, 
	       &plt->outliers->rmin) != 2) {
	return E_DATA;
    }

    for (i=0; i<n; i++) {
	if (fgets(line, sizeof line, fp) == NULL ||
	    sscanf(line, "%*d vals %lf", &plt->outliers->vals[i]) != 1) {
	    return E_DATA;
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
	grp = plotgroup_new(nplots);
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
    grp->show_mean = 1;
    grp->do_notches = 1;

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

	if (!err && grp->do_notches) {
	    gotnotch = check_confidence_interval(plt);
	    if (!gotnotch) {
		grp->do_notches = 0;
	    }
	}	    

	/* any outliers? */
	if (nout > 0) {
	    err = plot_retrieve_outliers(plt, nout, fp);
	} 

	if (!err && grp->show_mean) {
	    gotmean = maybe_get_plot_mean(plt, fp);
	    if (!gotmean) {
		grp->show_mean = 0;
	    }
	}
    }

 bailout:

    fclose(fp);

    if (!err) {
	if (grp->show_mean && grp->do_notches) {
	    grp->show_mean = 0;
	}
	err = write_gnuplot_boxplot(grp);
	if (!err) {
	    const char *pname = gretl_plotfile();

	    remove(fname);
	    gretl_copy_file(pname, fname);
	    remove(pname);
	    set_gretl_plotfile("");
	}
    } else if (err == E_DATA) {
	gretl_errmsg_set(_("boxplot file is corrupt"));
    }

    destroy_boxplots(grp);

    return err;
}
