/*
 *  Copyright (c) by Allin Cottrell and Riccardo Lucchetti
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
#include "plotspec.h"

/* FIXME cf. graphing.c */
static const char *auto_fit_string = "# plot includes automatic fitted line\n";

void plotspec_label_init (GPT_LABEL *lbl)
{
    lbl->text[0] = '\0';
    lbl->just = GP_JUST_LEFT;
    lbl->pos[0] = NADBL;
    lbl->pos[1] = NADBL;
}

GPT_SPEC *plotspec_new (void)
{
    GPT_SPEC *spec;
    int i;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->lines = malloc(MAX_PLOT_LINES * sizeof *spec->lines);

    if (spec->lines == NULL) {
	free(spec);
	return NULL;
    }

    for (i=0; i<MAX_PLOT_LINES; i++) {
	spec->lines[i].varnum = 0;
	spec->lines[i].title[0] = 0;
	spec->lines[i].formula[0] = 0;
	spec->lines[i].style[0] = 0;
	spec->lines[i].scale[0] = 0;
	spec->lines[i].yaxis = 1;
	spec->lines[i].type = 0;
	spec->lines[i].width = 1;
	spec->lines[i].ncols = 0;
    }

    for (i=0; i<4; i++) {
	spec->titles[i][0] = 0;
    }

    spec->literal = NULL;
    spec->n_literal = 0;

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	plotspec_label_init(&(spec->labels[i]));
    }

    spec->xtics[0] = 0;
    spec->mxtics[0] = 0;
    spec->fname[0] = 0;
    strcpy(spec->keyspec, "left top");
    spec->xzeroaxis = 0;

    for (i=0; i<3; i++) {
	spec->range[i][0] = NADBL;
	spec->range[i][1] = NADBL;
    }

    spec->b_ols = NULL;
    spec->b_quad = NULL;

    spec->code = PLOT_REGULAR;
    spec->flags = 0;
    spec->fit = PLOT_FIT_NONE;
    spec->fp = NULL;
    spec->data = NULL;
    spec->markers = NULL;
    spec->n_markers = 0;
    spec->labeled = NULL;
    spec->ptr = NULL;
    spec->reglist = NULL;
    spec->n_lines = 0;
    spec->nobs = 0;
    spec->okobs = 0;
    spec->boxwidth = 0;

    spec->termtype[0] = 0;

    return spec;
}

void plotspec_destroy (GPT_SPEC *spec)
{
    if (spec == NULL) {
	return;
    }

    if (spec->lines != NULL) {
	free(spec->lines);
    }

    if (spec->data != NULL) {
	free(spec->data);
    }

    if (spec->reglist != NULL) {
	free(spec->reglist);
    }

    if (spec->literal != NULL) {
	free_strings_array(spec->literal, spec->n_literal);
    }

    if (spec->markers != NULL) {
	free_strings_array(spec->markers, spec->n_markers);
    }

    if (spec->labeled != NULL) {
	free(spec->labeled);
    }

    gretl_matrix_free(spec->b_ols);
    gretl_matrix_free(spec->b_quad);

    free(spec);
}

int plotspec_add_line (GPT_SPEC *spec)
{
    GPT_LINE *lines;
    int n = spec->n_lines;

    lines = realloc(spec->lines, (n + 1) * sizeof *lines);
    if (lines == NULL) {
	return E_ALLOC;
    }

    spec->lines = lines;
    spec->n_lines += 1;

    lines[n].varnum = 0;
    lines[n].title[0] = 0;
    lines[n].formula[0] = 0;
    lines[n].style[0] = 0;
    lines[n].scale[0] = 0;
    lines[n].yaxis = 1;
    lines[n].type = 0;
    lines[n].width = 1;
    lines[n].ncols = 0;

    return 0;
}

static char *escape_quotes (const char *s)
{
    if (strchr(s, '"') == NULL) {
	return NULL;
    } else {
	int qcount = 0;
	char *ret, *r;
	const char *p = s;

	while (*p) {
	    if (*p == '"') qcount++;
	    p++;
	}

	ret = malloc(strlen(s) + 1 + qcount);
	if (ret == NULL) {
	    return NULL;
	}

	r = ret;
	while (*s) {
	    if (*s == '"') {
		*r++ = '\\';
		*r++ = '"';
	    } else {
		*r++ = *s;
	    }
	    s++;
	}
	*r = 0;

	return ret;
    }
}

static void 
gp_string (FILE *fp, const char *fmt, const char *s, int png)
{
#ifdef ENABLE_NLS  
    if (png && use_latin_2()) {
	char htmlstr[128];

	sprint_l2_to_html(htmlstr, s, sizeof htmlstr);
	fprintf(fp, fmt, htmlstr);
    } else {
	fprintf(fp, fmt, s); 
    }
#else
    fprintf(fp, fmt, s);
#endif
}

const char *gp_justification_string (int j)
{
    if (j == GP_JUST_LEFT) {
	return "left";
    } else if (j == GP_JUST_CENTER) {
	return "center";
    } else if (j == GP_JUST_RIGHT) {
	return "right";
    } else {
	return "left";
    }
}

static void 
print_plot_labelspec (const GPT_LABEL *lbl, int png, FILE *fp)
{
    char *label = escape_quotes(lbl->text);

    gp_string(fp, "set label \"%s\" ", (label != NULL)? 
	      label : lbl->text, png);

    gretl_push_c_numeric_locale();

    fprintf(fp, "at %g,%g %s%s\n", 
	    lbl->pos[0], lbl->pos[1],
	    gp_justification_string(lbl->just),
	    gnuplot_label_front_string());

    gretl_pop_c_numeric_locale();

    if (label != NULL) {
	free(label);
    }
}

static int print_data_labels (const GPT_SPEC *spec, FILE *fp)
{
    const double *x, *y;
    double xrange, yrange;
    double yoff;
    int t;

    if (spec->n_lines > 2 || 
	spec->lines[0].ncols != 2 || 
	spec->markers == NULL) {
	return 1;
    }

    x = spec->data;
    y = x + spec->nobs;

    xrange = spec->range[0][1] - spec->range[0][0];
    yrange = spec->range[1][1] - spec->range[1][0];

    if (xrange == 0.0 || yrange == 0.0) {
	double ymin = 1.0e+16, ymax = -1.0e+16;
	double xmin = 1.0e+16, xmax = -1.0e+16;

	for (t=0; t<spec->nobs; t++) {
	    if (!na(x[t]) && !na(y[t])) {
		if (yrange == 0.0) {
		    if (y[t] < ymin) {
			ymin = y[t];
		    } else if (y[t] > ymax) {
			ymax = y[t];
		    }
		}
		if (xrange == 0.0) {
		    if (x[t] < xmin) {
			xmin = x[t];
		    } else if (x[t] > xmax) {
			xmax = x[t];
		    }
		}		
	    }
	}
	if (yrange == 0.0) {
	    yrange = ymax - ymin;
	}
	if (xrange == 0.0) {
	    xrange = xmax - xmin;
	}
    }   

    yoff = 0.03 * yrange;

    fputs("# printing data labels\n", fp);

    gretl_push_c_numeric_locale();

    for (t=0; t<spec->nobs; t++) {
	double xoff = 0.0;

	if (spec->labeled != NULL && !spec->labeled[t]) {
	    /* printing only specified labels */
	    continue;
	}

	if (!na(x[t]) && !na(y[t])) {
	    if (x[t] > .90 * xrange) {
		xoff = -.02 * xrange;
	    }
	    fprintf(fp, "set label '%s' at %.8g,%.8g\n", spec->markers[t],
		    x[t] + xoff, y[t] + yoff);
	}
    }

    gretl_pop_c_numeric_locale();

    return 0;
}



static void print_fitted_title (const GPT_SPEC *spec, FILE *fp, int png)
{
    const char *s1 = spec->titles[1];
    const char *s2 = spec->titles[2];
    char title[128];

    if (*s1 == '\0' || *s2 == '\0') {
	return;
    }

    *title = '\0';

    if (spec->fit == PLOT_FIT_OLS) {
	sprintf(title, I_("%s versus %s (with least squares fit)"),
		s2, s1);
    } else if (spec->fit == PLOT_FIT_QUADRATIC) {
	sprintf(title, I_("%s versus %s (with quadratic fit)"),
		s2, s1);
    } else if (spec->fit == PLOT_FIT_LOESS) {
	sprintf(title, I_("%s versus %s (with lowess fit)"),
		s2, s1);
    }
	
    if (*title != '\0') {
	gp_string(fp, "set title '%s'\n", title, png);
    }
}

#define show_fit(s) (s->fit == PLOT_FIT_OLS || \
                     s->fit == PLOT_FIT_QUADRATIC || \
                     s->fit == PLOT_FIT_LOESS)

int plotspec_print (const GPT_SPEC *spec, FILE *fp)
{
    int i, k, t;
    int png = get_png_output(spec);
    int started_data_lines = 0;
    int n_lines = spec->n_lines;
    double *x[4];
    int any_y2 = 0;
    int miss = 0;

    if (spec->flags & GPT_FIT_HIDDEN) {
	n_lines--;
    }

    if (show_fit(spec)) {
	print_fitted_title(spec, fp, png);
    } else if (!string_is_blank(spec->titles[0])) {
	gp_string(fp, "set title '%s'\n", spec->titles[0], png);
    }

    if (!string_is_blank(spec->titles[1])) {
	gp_string(fp, "set xlabel '%s'\n", spec->titles[1], png);
    }

    if (!string_is_blank(spec->titles[2])) {
	gp_string(fp, "set ylabel '%s'\n", spec->titles[2], png);
    }

    if ((spec->flags & GPT_Y2AXIS) && !string_is_blank(spec->titles[3])) {
	gp_string(fp, "set y2label '%s'\n", spec->titles[3], png);
    }

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	if (!string_is_blank(spec->labels[i].text)) {
	    print_plot_labelspec(&(spec->labels[i]), png, fp);
	}
    }

    if (spec->xzeroaxis) {
	fputs("set xzeroaxis\n", fp);
    }

    fputs("set missing \"?\"\n", fp);

    if (strcmp(spec->keyspec, "none") == 0) {
	fputs("set nokey\n", fp);
    } else {
	fprintf(fp, "set key %s\n", spec->keyspec);
    }

    k = (spec->flags & GPT_Y2AXIS)? 3 : 2;

    gretl_push_c_numeric_locale();

    for (i=0; i<k; i++) {
	if (na(spec->range[i][0]) || na(spec->range[i][1]) ||
	    spec->range[i][0] == spec->range[i][1]) {
	    continue;
	}
	fprintf(fp, "set %srange [%.7g:%.7g]\n",
		(i == 0)? "x" : (i == 1)? "y" : "y2",
		spec->range[i][0], spec->range[i][1]);
    }

    if (spec->boxwidth > 0 && spec->boxwidth < 1) {
	fprintf(fp, "set boxwidth %.3f\n", (double) spec->boxwidth);
    }

    gretl_pop_c_numeric_locale();

    /* customized xtics? */
    if (!string_is_blank(spec->xtics)) {
	fprintf(fp, "set xtics %s\n", spec->xtics);
    }
    if (!string_is_blank(spec->mxtics)) {
	fprintf(fp, "set mxtics %s\n", spec->mxtics);
    }

    if (spec->flags & GPT_Y2AXIS) {
	/* using two y axes */
	fputs("set ytics nomirror\n", fp);
	fputs("set y2tics\n", fp);
    } else if (spec->flags & GPT_MINIMAL_BORDER) {
	/* suppressing part of border */
	fputs("set border 3\n", fp);
	if (string_is_blank(spec->xtics)) {
	    fputs("set xtics nomirror\n", fp);
	}
	fputs("set ytics nomirror\n", fp);
    }

    /* in case of plots that are editable (see gui client), it is
       important to write out the comment string that identifies the
       sort of graph, so that it will be recognized by type when
       it is redisplayed */

    write_plot_type_string(spec->code, fp);

    if (spec->n_literal > 0) {
	fprintf(fp, "# literal lines = %d\n", spec->n_literal);
	for (i=0; i<spec->n_literal; i++) {
	    if (spec->literal[i] != NULL && *spec->literal[i] != '\0') {
		fprintf(fp, "%s\n", spec->literal[i]);
	    } else {
		fputs("# empty line!\n", fp);
	    }
	}
    }

    if (spec->flags & GPT_AUTO_FIT) {
	fputs(auto_fit_string, fp);
    }

    if ((spec->code == PLOT_FREQ_SIMPLE ||
	 spec->code == PLOT_FREQ_NORMAL ||
	 spec->code == PLOT_FREQ_GAMMA) && gnuplot_has_style_fill()) {
	fputs("set style fill solid 0.8\n", fp);
    }  

    if (spec->flags & GPT_ALL_MARKERS) {
	print_data_labels(spec, fp);
    }

    fputs("plot \\\n", fp);

    for (i=0; i<n_lines; i++) {
	if ((spec->flags & GPT_Y2AXIS) && spec->lines[i].yaxis != 1) {
	    any_y2 = 1;
	    break;
	}
    }

    for (i=0; i<n_lines; i++) {

	if (strcmp(spec->lines[i].scale, "NA")) {
	    if (!strcmp(spec->lines[i].scale, "1.0")) {
		fputs("'-' using 1", fp);
		for (k=2; k<=spec->lines[i].ncols; k++) {
		    fprintf(fp, ":%d", k);
		}
		fputc(' ', fp);
	    } else {
		fprintf(fp, "'-' using 1:($2*%s) ", 
			spec->lines[i].scale);
	    }
	} else {
	    fprintf(fp, "%s ", spec->lines[i].formula); 
	}

	if ((spec->flags & GPT_Y2AXIS) && spec->lines[i].yaxis != 1) {
	    fprintf(fp, "axes x1y%d ", spec->lines[i].yaxis);
	}

	gp_string(fp, "title '%s", spec->lines[i].title, png);

	if (any_y2) {
	    if (spec->lines[i].yaxis == 1) {
		fprintf(fp, " (%s)' ", I_("left"));
	    } else {
		fprintf(fp, " (%s)' ", I_("right"));
	    }
	} else {
	    fputs("' ", fp);
	}

	fprintf(fp, "w %s", spec->lines[i].style);
	if (spec->lines[i].type != 0) {
	    fprintf(fp, " lt %d", spec->lines[i].type);
	}
	if (spec->lines[i].width != 1) {
	    fprintf(fp, " lw %d", spec->lines[i].width);
	}

	if (i == n_lines - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp);
	}
    } 

    miss = 0;

    /* supply the data to gnuplot inline */

    gretl_push_c_numeric_locale();

    for (i=0; i<spec->n_lines; i++) { 
	int j;

	if (spec->lines[i].ncols == 0) {
	    continue;
	}

	if (!started_data_lines) {
	    x[0] = spec->data;
	    x[1] = x[0] + spec->nobs;
	    started_data_lines = 1;
	} 

	x[2] = x[1] + spec->nobs;
	x[3] = x[2] + spec->nobs;

	for (t=0; t<spec->nobs; t++) {
	    if (na(x[0][t])) {
		fputs("? ", fp);
		miss = 1;
	    } else {
		fprintf(fp, "%.8g ", x[0][t]);
	    }

	    for (j=1; j<spec->lines[i].ncols; j++) {
		if (na(x[j][t])) {
		    fputs("? ", fp);
		    miss = 1;
		} else {
		    fprintf(fp, "%.8g ",  x[j][t]);
		}
	    }

	    if (spec->markers != NULL && i == 0) {
		fprintf(fp, " # %s", spec->markers[t]);
	    }

	    fputc('\n', fp);
	}

	fputs("e\n", fp);

	x[1] += (spec->lines[i].ncols - 1) * spec->nobs;
    }

    gretl_pop_c_numeric_locale();

    return miss;
}

int plotspec_get_term_string (const GPT_SPEC *spec, char *termstr)
{
    int cmds = 0;

    if (!strcmp(spec->termtype, "postscript color")) {
	strcpy(termstr, "postscript eps color"); 
    } else if (!strcmp(spec->termtype, "postscript")) {
	strcpy(termstr, "postscript eps"); 
    } else if (!strcmp(spec->termtype, "PDF")) {
	strcpy(termstr, "pdf");
    } else if (!strcmp(spec->termtype, "fig")) {
	strcpy(termstr, "fig");
    } else if (!strcmp(spec->termtype, "latex")) {
	strcpy(termstr, "latex");
    } else if (!strcmp(spec->termtype, "png")) { 
	const char *png_str = 
	    get_gretl_png_term_line(spec->code, 0);

	strcpy(termstr, png_str + 9);
    } else if (!strcmp(spec->termtype, "emf color")) {
	const char *emf_str = 
	    get_gretl_emf_term_line(spec->code, 1);

	strcpy(termstr, emf_str + 9);
    } else if (!strcmp(spec->termtype, "plot commands")) { 
	cmds = 1;
    } else {
	strcpy(termstr, spec->termtype);
    }

    return cmds;
}

int plotspec_ship_out (GPT_SPEC *spec, char *fname)
{
    FILE *fp = NULL;
    int dump = 0;
    int err = 0, miss;
    char termstr[72];

    dump = plotspec_get_term_string(spec, termstr);

    if (dump) {  
	/* dump of gnuplot commands to named file */
	if (fname == NULL) {
	    return 1;  /* impossible */
	}
	fp = gretl_fopen(fname, "w");
	if (fp == NULL) {
	    return 1;
	}
    } else {     
	/* output to gnuplot, for screen or other "term" */
	if (spec->fp == NULL) {
	    fp = gretl_fopen(gretl_plotfile(), "w");
	}
	if (fp == NULL) {
	    return 1;
	}

	if (fname != NULL) { 
#ifdef ENABLE_NLS
	    fprint_gnuplot_encoding(termstr, fp);
#endif 
	    /* file, not screen display */
	    fprintf(fp, "set term %s\n", termstr);
	    fprintf(fp, "set output '%s'\n", fname);
	}
    }

    if (strstr(termstr, "png")) {
	set_png_output(spec);
    }

    miss = plotspec_print(spec, fp);
    fflush(fp);

    if (dump) {
	/* we're finished */
	fclose(fp);
    }
    
    if (!dump) {
	char plotcmd[MAXLEN];
#ifdef WIN32
	int winshow = 0;

	if (fname == NULL) { 
	    /* sending plot to screen */
	    fputs("pause -1\n", fp);
	    winshow = 1;
	} 
#endif
	fclose(fp);
	spec->fp = NULL;
	sprintf(plotcmd, "\"%s\" \"%s\"", gretl_gnuplot_path(), gretl_plotfile());
#ifdef WIN32
	if (winshow) {
	    err = (WinExec(plotcmd, SW_SHOWNORMAL) < 32);
	} else {
	    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
	}
#else
	if (gretl_spawn(plotcmd)) err = 1;
#endif 
    }

    if (miss) {
	err = 2;
    }

    return err;
}

static void set_fitted_line (GPT_SPEC *spec, FitType f)
{
    char *formula = spec->lines[1].formula;
    char *title = spec->lines[1].title;
    const double *b;

    if (f == PLOT_FIT_OLS) {
	b = spec->b_ols->val;
	sprintf(title, "Y = %#.3g %c %#.3gX", b[0],
		(b[1] > 0)? '+' : '-', fabs(b[1]));
	gretl_push_c_numeric_locale();
	sprintf(formula, "%g + %g*x", b[0], b[1]);
	gretl_pop_c_numeric_locale();
    } else if (f == PLOT_FIT_QUADRATIC) {
	b = spec->b_quad->val;
	sprintf(title, "Y = %#.3g %c %#.3gX %c %#.3gX^2", b[0],
		(b[1] > 0)? '+' : '-', fabs(b[1]),
		(b[2] > 0)? '+' : '-', fabs(b[2]));
	gretl_push_c_numeric_locale();
	sprintf(formula, "%g + %g*x + %g*x**2", b[0], b[1], b[2]);
	gretl_pop_c_numeric_locale();
    }

    strcpy(spec->lines[1].scale, "NA");
    strcpy(spec->lines[1].style, "lines");
    *spec->titles[0] = '\0';

    spec->fit = f;
}

int plotspec_add_fit (GPT_SPEC *spec, FitType f)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    const double *py;
    const double *px;
    int T = spec->okobs;
    int i, t, k;
    int err = 0;

    if ((f == PLOT_FIT_OLS && spec->b_ols != NULL) ||
	(f == PLOT_FIT_QUADRATIC && spec->b_quad != NULL)) {
	/* just activate existing setup */
	set_fitted_line(spec, f);
	return 0;
    }
	
    if (f == PLOT_FIT_OLS) {
	k = 2;
    } else if (f == PLOT_FIT_QUADRATIC) {
	k = 3;
    } else if (f == PLOT_FIT_LOESS) {
	k = 1;
	fprintf(stderr, "LOESS: not ready\n");
	return 1;
    } else {
	return E_DATA;
    }

    y = gretl_column_vector_alloc(T);
    X = gretl_matrix_alloc(T, k);
    b = gretl_column_vector_alloc(k);

    if (y == NULL || X == NULL || b == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    px = spec->data;
    py = px + spec->nobs;

    i = 0;
    for (t=0; t<spec->nobs; t++) {
	if (!na(py[t]) && !na(px[t])) {
	    y->val[i] = py[t];
	    if (f == PLOT_FIT_LOESS) {
		gretl_matrix_set(X, i, 0, px[t]);
	    } else {
		gretl_matrix_set(X, i, 0, 1.0);
		gretl_matrix_set(X, i, 1, px[t]);
		if (f == PLOT_FIT_QUADRATIC) {
		    gretl_matrix_set(X, i, 2, px[t] * px[t]);
		}
	    }
	    i++;
	}
    }

#if 0
    gretl_matrix_print(y, "y");
    gretl_matrix_print(X, "X");
#endif

    err = gretl_matrix_ols(y, X, b, NULL, NULL, NULL);

    if (!err && spec->n_lines == 1) {
	err = plotspec_add_line(spec);
    }

    if (!err) {
	if (f == PLOT_FIT_OLS) {
	    spec->b_ols = b;
	} else if (f == PLOT_FIT_QUADRATIC) {
	    spec->b_quad = b;
	}
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);

    if (err) {
	gretl_matrix_free(b);
    } else {
	set_fitted_line(spec, f);
    }

    return err;
}
