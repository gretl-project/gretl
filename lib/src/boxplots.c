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
    double xbase;
    char varname[VNAMELEN];
    char *bool;
    OUTLIERS *outliers;
} BOXPLOT;

typedef struct {
    int nplots;
    int show_outliers;
    char *numbers;
    BOXPLOT *plots;
} PLOTGROUP;

static void destroy_boxplots (PLOTGROUP *grp)
{
    int i;

    for (i=0; i<grp->nplots; i++) { 
	free(grp->plots[i].outliers);
	free(grp->plots[i].bool);
    }

    free(grp->plots);
    free(grp->numbers);
    free(grp);
}

static int gnuplot_do_boxplot (PLOTGROUP *grp, const int *list,
			       const DATAINFO *pdinfo)
{
    FILE *fp;
    BOXPLOT *bp;
    int n = grp->nplots;
    double loff, ymax;
    int i, err = 0;

    err = gnuplot_init(PLOT_BOXPLOTS, &fp);
    if (err) {
	return err;
    }

    fprintf(fp, "set xrange[0:%d]\n", n + 1);
    fputs("set ytics nomirror\n"
	  "set noxtics\n"
	  "set border 2\n"
	  "set bmargin 3\n", fp);

    ymax = grp->plots[0].max;

    for (i=1; i<n; i++) {
	if (grp->plots[i].max > ymax) {
	    ymax = grp->plots[i].max;
	}
    }

    loff = -ymax / 25;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	fprintf(fp, "set label \"%s\" at %d,%g center\n", 
		pdinfo->varname[list[i+1]], i+1, loff);
    }

    fprintf(fp, "set boxwidth %g absolute\n", 1.0 / (n + 2));

    fputs("plot \\\n", fp);
    fputs("'-' using 1:3:2:5:4 with candlesticks lt 2 lw 2 "
	  "notitle whiskerbars 0.5, \\\n", fp);
    fputs("'-' using 1:2:2:2:2 with candlesticks lt -1 notitle, \\\n", fp);
    fputs("'-' using 1:2 with points pt 1 notitle\n", fp);

    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	fprintf(fp, "%d %g %g %g %g\n", i+1, bp->min, bp->lq, bp->uq, bp->max);
    }
    fputs("e\n", fp);

    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	fprintf(fp, "%d %g\n", i+1, bp->median);
    }
    fputs("e\n", fp);

    for (i=0; i<n; i++) {
	bp = &grp->plots[i];
	fprintf(fp, "%d %g\n", i+1, bp->mean);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();
    
    fclose(fp);

    return gnuplot_make_graph();
}

int boxplots (int *list, char **bools, 
	      double ***pZ, const DATAINFO *pdinfo, 
	      gretlopt opt)
{
    int i, j, n = pdinfo->t2 - pdinfo->t1 + 1;
    double *x;
    PLOTGROUP *plotgrp;
    int width = 576, height = 448;
    int err = 0;

    x = mymalloc(n * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    plotgrp = mymalloc(sizeof *plotgrp);
    if (plotgrp == NULL) {
	free(x);
	return E_ALLOC;
    }

    plotgrp->saved = 0;
    plotgrp->show_outliers = 0;

    plotgrp->nplots = list[0];
    plotgrp->plots = mymalloc(plotgrp->nplots * sizeof *plotgrp->plots);
    if (plotgrp->plots == NULL) {
	free(plotgrp);
	free(x);
	return E_ALLOC;
    }

    for (i=0, j=0; i<plotgrp->nplots; i++, j++) {
	n = ztox(list[i+1], x, (const double **) *pZ, pdinfo);
	if (n < 2) {
	    errbox(_("Dropping %s: insufficient observations"),
		   pdinfo->varname[list[i+1]]);
	    gretl_list_delete_at_pos(list, i+1);
	    if (list[0] == 0) {
		free(plotgrp->plots);
		free(plotgrp);
		free(x);
		return E_DATA;
	    } else {
		plotgrp->nplots -= 1;
		i--;
		continue;
	    }
	}

	plotgrp->plots[i].outliers = NULL;
	plotgrp->plots[i].mean = gretl_mean(0, n-1, x);
	qsort(x, n, sizeof *x, gretl_compare_doubles);
	plotgrp->plots[i].min = x[0];
	plotgrp->plots[i].max = x[n-1];
	quartiles(x, n, &plotgrp->plots[i]);
	plotgrp->plots[i].n = n;

	if (opt & OPT_O) {
	    /* notched boxplots wanted */
	    if (median_interval(x, n, &plotgrp->plots[i].conf[0],
				&plotgrp->plots[i].conf[1])) {
		errbox(_("Couldn't obtain confidence interval"));
		plotgrp->plots[i].conf[0] = 
		    plotgrp->plots[i].conf[1] = NADBL;
	    }
	} else {
	    plotgrp->plots[i].conf[0] = plotgrp->plots[i].conf[1] = NADBL;
	}

	strcpy(plotgrp->plots[i].varname, pdinfo->varname[list[i+1]]);

	if (bools) { 
	    plotgrp->plots[i].bool = bools[j];
	} else {
	    plotgrp->plots[i].bool = NULL;
	}
    }

    err = gnuplot_do_boxplot(plotgrp, list, pdinfo);

    destroy_boxplots(plotgrp);   
    
    return err;
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
    int i, k, v, nvars, nbool, err = 0;
    int n = pdinfo->n, origv = pdinfo->v;
    char *tok, *s = NULL, **bools = NULL;
    int *list = NULL;

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

    list = malloc((nvars + 1) * sizeof *list);
    bools = malloc(nvars * sizeof *bools);

    if (list == NULL || bools == NULL) {
	free(list);
	free(bools);
	free(s);
	return 1;
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
		bools[i-1] = malloc(strlen(tok) + 1);
		strcpy(bools[i-1], tok);
		nbool++;
	    } else {
		err = 1;
	    }
	} else {
	    if (isdigit(tok[0])) { 
		v = atoi(tok);
		if (v < origv) {
		    list[++i] = v;
		} else {
		    errbox(_("got invalid variable number %d"), v);
		    err = 1;
		}
	    } else if (isalpha(tok[0])) {
		v = series_index(pdinfo, tok);
		if (v < origv) {
		    list[++i] = v;
		} else {
		    errbox(_("got invalid varname '%s'"), tok);
		    err = 1;
		}
	    } else {
		errbox(_("got invalid field '%s'"), tok);
		err = 1; 
	    }
	}
    }

    /* now we add nbool new variables, with ID numbers origv,
       origv + 1, and so on.  These are the original variables
       that have boolean conditions attached, masked by those
       conditions */

    k = origv;
    nbool = 0;
    for (i=1; i<=list[0] && !err; i++) {
	if (bools[i-1] != NULL) {
	    char formula[80];
	    int t;
	    
	    sprintf(formula, "bool_%d = %s", i-1, bools[i-1]);
	    err = generate(formula, pZ, pdinfo, OPT_P, NULL);
	    if (err) {
		char *msg = 
		    g_strdup_printf(_("boxplots: generation of dummy variable failed\n%s"), 
				    gretl_errmsg_get());

		errbox(msg);
		g_free(msg);
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
    }

    if (!err) {
	err = boxplots(list, bools, pZ, pdinfo, opt);
    } 
    
    free(list);
    free(bools); /* the bool[i]s are now attached to the plots */
    free(s);

    if (nbool) {
	dataset_drop_last_variables(nbool, pZ, pdinfo);
    }
    
    return err;
}
