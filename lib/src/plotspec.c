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
#include "plotspec.h"

GPT_SPEC *plotspec_new (void)
{
    GPT_SPEC *spec;
    int i;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->lines = NULL;
    spec->n_lines = 0;

    spec->labels = NULL;
    spec->n_labels = 0;

    spec->literal = NULL;
    spec->n_literal = 0;

    for (i=0; i<4; i++) {
	spec->titles[i][0] = 0;
    }

    *spec->xvarname = '\0';
    *spec->yvarname = '\0';

    spec->xtics[0] = 0;
    spec->mxtics[0] = 0;
    spec->fname[0] = 0;
    strcpy(spec->keyspec, "left top");

    for (i=0; i<4; i++) {
	spec->range[i][0] = NADBL;
	spec->range[i][1] = NADBL;
	if (i < 3) {
	    spec->logbase[i] = 0.0;
	}
    }

    spec->b_ols = NULL;
    spec->b_quad = NULL;
    spec->b_inv = NULL;

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
    spec->nobs = 0;
    spec->okobs = 0;
    spec->pd = 0;
    spec->boxwidth = 0;
    spec->samples = 0;
    spec->border = GP_BORDER_DEFAULT;
    spec->bmargin = 0;

    spec->termtype = GP_TERM_NONE;

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

    if (spec->labels != NULL) {
	free(spec->labels);
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
    gretl_matrix_free(spec->b_inv);

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
    lines[n].title[0] = '\0';
    lines[n].formula[0] = '\0';
    lines[n].style[0] = '\0';
    lines[n].scale[0] = '\0';
    lines[n].rgb[0] = '\0';
    lines[n].yaxis = 1;
    lines[n].type = LT_NONE;
    lines[n].ptype = 0;
    lines[n].width = 1;
    lines[n].ncols = 0;
    lines[n].whiskwidth = 0;
    lines[n].flags = 0;

    return 0;
}

static void copy_line_content (GPT_LINE *targ, GPT_LINE *src)
{
    targ->varnum = src->varnum;
    strcpy(targ->title, src->title);
    strcpy(targ->formula, src->formula);
    strcpy(targ->style, src->style);
    strcpy(targ->scale, src->scale);
    strcpy(targ->rgb, src->rgb);
    targ->yaxis = src->yaxis;
    targ->type = src->type;
    targ->ptype = src->ptype;
    targ->width = src->width;
    targ->ncols = src->ncols;
    targ->whiskwidth = src->whiskwidth;
    targ->flags = src->flags;
}

int plotspec_delete_line (GPT_SPEC *spec, int i)
{
    GPT_LINE *lines = spec->lines;
    int j, n = spec->n_lines;

    if (i < 0 || i >= n) {
	return E_DATA;
    }

    for (j=i; j<n-1; j++) {
	copy_line_content(&lines[j], &lines[j+1]);
    }

    spec->n_lines -= 1;

    lines = realloc(spec->lines, (n - 1) * sizeof *lines);
    if (lines == NULL) {
	return E_ALLOC;
    }

    spec->lines = lines;

    return 0;
}

GPT_LINE *plotspec_clone_lines (GPT_SPEC *spec, int *err)
{
    GPT_LINE *lines = NULL;
    int i;

    if (spec->n_lines == 0) {
	/* no-op */
	return NULL;
    }

    lines = malloc(spec->n_lines * sizeof *lines);
    if (lines == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<spec->n_lines; i++) {
	    copy_line_content(&lines[i], &spec->lines[i]);
	}
    }

    return lines;
}

int plotspec_add_label (GPT_SPEC *spec)
{
    GPT_LABEL *labels;
    int n = spec->n_labels;

    labels = realloc(spec->labels, (n + 1) * sizeof *labels);
    if (labels == NULL) {
	return E_ALLOC;
    }

    spec->labels = labels;
    spec->n_labels += 1;

    labels[n].text[0] = '\0';
    labels[n].pos[0] = NADBL;
    labels[n].pos[1] = NADBL;
    labels[n].just = GP_JUST_LEFT;

    return 0;
}

static void copy_label_content (GPT_LABEL *targ, GPT_LABEL *src)
{
    strcpy(targ->text, src->text);
    targ->pos[0] = src->pos[0];
    targ->pos[1] = src->pos[1];
    targ->just = src->just;
}

int plotspec_delete_label (GPT_SPEC *spec, int i)
{
    GPT_LABEL *labels = spec->labels;
    int j, n = spec->n_labels;
    int err = 0;

    if (i < 0 || i >= n) {
	return E_DATA;
    }

    for (j=i; j<n-1; j++) {
	copy_label_content(&labels[j], &labels[j+1]);
    }

    spec->n_labels -= 1;

    if (spec->n_labels == 0) {
	free(spec->labels);
	spec->labels = NULL;
    } else {
	labels = realloc(spec->labels, (n - 1) * sizeof *labels);
	if (labels == NULL) {
	    err = E_ALLOC;
	} else {
	    spec->labels = labels;
	}
    }

    return err;
}

GPT_LABEL *plotspec_clone_labels (GPT_SPEC *spec, int *err)
{
    GPT_LABEL *labels = NULL;
    int i;

    if (spec->n_labels == 0) {
	/* no-op */
	return NULL;
    }

    labels = malloc(spec->n_labels * sizeof *labels);
    if (labels == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<spec->n_labels; i++) {
	    copy_label_content(&labels[i], &spec->labels[i]);
	}
    }

    return labels;
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
print_plot_labelspec (const GPT_LABEL *lbl, FILE *fp)
{
    char *label = escape_quotes(lbl->text);

    fprintf(fp, "set label \"%s\" ", (label != NULL)? 
	    label : lbl->text);

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
	    fprintf(fp, "set label \"%s\" at %.8g,%.8g\n", spec->markers[t],
		    x[t] + xoff, y[t] + yoff);
	}
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

void print_auto_fit_string (FitType fit, FILE *fp)
{
    if (fit == PLOT_FIT_OLS) {
	fputs("# plot includes automatic fit: OLS\n", fp);
    } else if (fit == PLOT_FIT_QUADRATIC) {
	fputs("# plot includes automatic fit: quadratic\n", fp);
    } else if (fit == PLOT_FIT_INVERSE) {
	fputs("# plot includes automatic fit: inverse\n", fp);
    } else if (fit == PLOT_FIT_LOESS) {
	fputs("# plot includes automatic fit: loess\n", fp);
    }
}

void print_plot_ranges_etc (const GPT_SPEC *spec, FILE *fp)
{
    const char *rstrs[] = {
	"x", "y", "y2", "t"
    };
    int i;

    gretl_push_c_numeric_locale();

    for (i=0; i<4; i++) {
	if (i < 3 && spec->logbase[i] > 0.0) {
	    fprintf(fp, "set logscale %s %g\n", rstrs[i], spec->logbase[i]);
	}

	if (na(spec->range[i][0]) || na(spec->range[i][1]) ||
	    spec->range[i][0] == spec->range[i][1]) {
	    continue;
	}

	if ((i == 2 && !(spec->flags & GPT_Y2AXIS)) ||
	    (i == 3 && !(spec->flags & GPT_PARAMETRIC))) {
	    continue;
	}

	fprintf(fp, "set %srange [%.7g:%.7g]\n", rstrs[i], 
		spec->range[i][0], spec->range[i][1]);

	if (i == 4 && spec->code == PLOT_PROB_DIST && spec->samples == 0) {
	    int ns = spec->range[i][1] - spec->range[i][0] + 1;

	    fprintf(fp, "set samples %d\n", ns);
	}
    }

    if (spec->boxwidth > 0 && spec->boxwidth < 1) {
	fprintf(fp, "set boxwidth %.3f\n", (double) spec->boxwidth);
    } else if (spec->boxwidth < 0 && spec->boxwidth > -1) {
	fprintf(fp, "set boxwidth %g absolute\n", (double) -spec->boxwidth);
    }

    gretl_pop_c_numeric_locale();
}

static int blank_user_line (const GPT_SPEC *spec, int i)
{
    return ((spec->lines[i].flags & GP_LINE_USER) && 
	    spec->lines[i].formula[0] == '\0');

}

static void print_user_lines_info (const GPT_SPEC *spec, FILE *fp)
{
    int i, n = 0;

    for (i=0; i<spec->n_lines; i++) {
	if (spec->lines[i].flags & GP_LINE_USER) {
	    n++;
	}
    }

    if (n > 0) {
	fprintf(fp, "# %d user-defined lines: ", n);
	for (i=0; i<spec->n_lines; i++) {
	    if (spec->lines[i].flags & GP_LINE_USER) {
		fprintf(fp, "%d ", i);
	    }
	}
	fputc('\n', fp);
    }
}

static int more_lines (const GPT_SPEC *spec, int i, int skipline)
{
    if (i == spec->n_lines - 1) {
	/* at last line */
	return 0;
    } else {
	int j;

	for (j=i+1; j<spec->n_lines; j++) {
	    if (j != skipline && !blank_user_line(spec, j)) {
		return 1;
	    }
	}
	return 0;
    }
}

static void print_linestyle (const GPT_SPEC *spec, int i, FILE *fp)
{
    GPT_LINE *line;
    int done = 0;

    fprintf(fp, "set style line %d ", i+1);

    if (i < spec->n_lines) {
	line = &spec->lines[i];
	if (*line->rgb != '\0' && line->type != LT_NONE) {
	    fprintf(fp, "lc rgb \"%s\" lt %d\n", line->rgb, line->type);
	    done = 1;
	} else if (*line->rgb != '\0') {
	    fprintf(fp, "lc rgb \"%s\"\n", line->rgb);
	    done = 1;
	}
    }

    if (!done) {
	const gretlRGB *color = get_graph_color(i);
	char cstr[8];

	print_rgb_hash(cstr, color);
	fprintf(fp, "lc rgb \"%s\"\n", cstr);
    }
}

static void write_styles_from_plotspec (const GPT_SPEC *spec, FILE *fp)
{
    char cstr[8];
    int i;
    
    if (frequency_plot_code(spec->code)) {
	const gretlRGB *color = get_graph_color(BOXCOLOR);

	print_rgb_hash(cstr, color);
	fprintf(fp, "set style line 1 lc rgb \"%s\"\n", cstr);
	fputs("set style line 2 lc rgb \"#000000\"\n", fp);
    } else if (spec->code == PLOT_RQ_TAU) {
	fputs("set style line 1 lc rgb \"#000000\"\n", fp);
	for (i=1; i<BOXCOLOR; i++) {
	    print_linestyle(spec, i, fp);
	}
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    print_linestyle(spec, i, fp);
	}
    }

    fputs("set style increment user\n", fp);
}

static int print_point_type (GPT_LINE *line)
{
    return line->ptype != 0 && 
	(!strcmp(line->style, "points") ||
	 !strcmp(line->style, "linespoints"));
}

#define show_fit(s) (s->fit == PLOT_FIT_OLS || \
                     s->fit == PLOT_FIT_QUADRATIC || \
                     s->fit == PLOT_FIT_INVERSE || \
                     s->fit == PLOT_FIT_LOESS)

int plotspec_print (const GPT_SPEC *spec, FILE *fp)
{
    int i, k, t;
    int png = get_png_output(spec);
    int mono = (spec->flags & GPT_MONO);
    int started_data_lines = 0;
    double et, yt;
    double *x[5];
    int skipline = -1;
    int any_y2 = 0;
    int miss = 0;

    if (spec->pd > 0) {
	fprintf(fp, "# timeseries %d", spec->pd);
	if (spec->flags & GPT_LETTERBOX) {
	    fputs(" (letterbox)\n", fp);
	} else {
	    fputc('\n', fp);
	}
    }

    if (!mono && gnuplot_has_rgb()) {
	write_styles_from_plotspec(spec, fp);
    }

    if (!string_is_blank(spec->titles[0])) {
	fprintf(fp, "set title \"%s\"\n", spec->titles[0]);
    }

    if (!string_is_blank(spec->titles[1])) {
	fprintf(fp, "set xlabel \"%s\"\n", spec->titles[1]);
    }

    if (!string_is_blank(spec->titles[2])) {
	fprintf(fp, "set ylabel \"%s\"\n", spec->titles[2]);
    }

    if ((spec->flags & GPT_Y2AXIS) && !string_is_blank(spec->titles[3])) {
	fprintf(fp, "set y2label \"%s\"\n", spec->titles[3]);
    }

    for (i=0; i<spec->n_labels; i++) {
	if (!string_is_blank(spec->labels[i].text)) {
	    print_plot_labelspec(&spec->labels[i], fp);
	}
    }

    if (spec->flags & GPT_XZEROAXIS) {
	fputs("set xzeroaxis\n", fp);
    }

    if (spec->flags & GPT_YZEROAXIS) {
	fputs("set yzeroaxis\n", fp);
    }    

    if (spec->flags & GPT_PARAMETRIC) {
	fputs("set parametric\n", fp);
	if (spec->samples > 0) {
	    fprintf(fp, "set samples %d\n", spec->samples);
	}
    }

    gnuplot_missval_string(fp);

    if (strcmp(spec->keyspec, "none") == 0) {
	fputs("set nokey\n", fp);
    } else {
	fprintf(fp, "set key %s\n", spec->keyspec);
    }

    print_plot_ranges_etc(spec, fp);

    /* customized xtics? */
    if (!strcmp(spec->xtics, "none")) {
	fputs("set noxtics\n", fp);
    } else {
	if (!string_is_blank(spec->xtics)) {
	    fprintf(fp, "set xtics %s\n", spec->xtics);
	}
	if (!string_is_blank(spec->mxtics)) {
	    fprintf(fp, "set mxtics %s\n", spec->mxtics);
	}
    }

    if (spec->flags & GPT_Y2AXIS) {
	/* using two y axes */
	fputs("set ytics nomirror\n", fp);
	fputs("set y2tics\n", fp);
    } else if (spec->border != GP_BORDER_DEFAULT) {
	/* suppressing all or part of border */
	if (spec->border == 0) {
	    fputs("unset border\n", fp);
	} else {
	    fprintf(fp, "set border %d\n", spec->border);
	}
	if (string_is_blank(spec->xtics)) {
	    fputs("set xtics nomirror\n", fp);
	}
	fputs("set ytics nomirror\n", fp);
    }

    if (spec->bmargin > 0) {
	fprintf(fp, "set bmargin %d\n", spec->bmargin);
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

    if (show_fit(spec)) {
	print_auto_fit_string(spec->fit, fp);
    }

    print_user_lines_info(spec, fp);

    if ((spec->code == PLOT_FREQ_SIMPLE ||
	 spec->code == PLOT_FREQ_NORMAL ||
	 spec->code == PLOT_FREQ_GAMMA) && gnuplot_has_style_fill()) {
	if (mono) {
	    fputs("set style fill solid 0.3\n", fp);
	} else {
	    fputs("set style fill solid 0.6\n", fp);
	}
    } else if (spec->code == PLOT_FORECAST) {
	fputs("set style fill solid 0.4\n", fp);
    }

    if (spec->flags & GPT_ALL_MARKERS) {
	print_data_labels(spec, fp);
    }

    if (spec->flags & GPT_FIT_HIDDEN) {
	skipline = 1;
    }

    fputs("plot \\\n", fp);

    for (i=0; i<spec->n_lines; i++) {
	if (i == skipline || blank_user_line(spec, i)) {
	    continue;
	}
	if ((spec->flags & GPT_Y2AXIS) && spec->lines[i].yaxis != 1) {
	    any_y2 = 1;
	    break;
	}
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<spec->n_lines; i++) {
	GPT_LINE *line = &spec->lines[i];

	if (i == skipline || blank_user_line(spec, i)) {
	    if (i < spec->n_lines - 1) {
		continue;
	    } else {
		break;
	    }
	}

	if (!strcmp(line->scale, "NA")) {
	    fprintf(fp, "%s ", line->formula); 
	} else {
	    if (!strcmp(line->scale, "1.0")) {
		fputs("'-' using 1", fp);
		if (line->ncols == 5) {
		    /* Note: boxplot candlesticks, hard-wired! */
		    fputs(":3:2:5:4", fp);
		} else if (line->ncols == 2) {
		    if (line->flags & GP_LINE_BOXDATA) {
			/* boxplot median, hard-wired */
			fputs(":2:2:2:2", fp);
		    } else {
			fputs(":($2)", fp);
		    }
		} else {
		    for (k=2; k<=line->ncols; k++) {
			fprintf(fp, ":%d", k);
		    }
		}
		fputc(' ', fp);
	    } else {
		fprintf(fp, "'-' using 1:($2*%s) ", line->scale);
	    }
	} 

	if ((spec->flags & GPT_Y2AXIS) && line->yaxis != 1) {
	    fprintf(fp, "axes x1y%d ", line->yaxis);
	}

	if (*line->title == '\0') {
	    fputs("notitle ", fp);
	} else {
	    fprintf(fp, "title \"%s", line->title);
	    if (any_y2) {
		if (line->yaxis == 1) {
		    fprintf(fp, " (%s)\" ", G_("left"));
		} else {
		    fprintf(fp, " (%s)\" ", G_("right"));
		}
	    } else {
		fputs("\" ", fp);
	    }
	}

	fprintf(fp, "w %s", line->style);

	if (line->type != LT_NONE) {
	    fprintf(fp, " lt %d", line->type);
	}

	if (print_point_type(line)) {
	    fprintf(fp, " pt %d", line->ptype);
	}

	if (line->width != 1) {
	    fprintf(fp, " lw %d", line->width);
	}

	if (line->whiskwidth > 0) {
	    fprintf(fp, " whiskerbars %g", line->whiskwidth);
	}

	if (more_lines(spec, i, skipline)) {
	    fputs(", \\", fp);
	} 

	fputc('\n', fp);
    } 

    miss = 0;

    /* supply the data to gnuplot inline */

    for (i=0; i<spec->n_lines; i++) { 
	int j, ncols = spec->lines[i].ncols;

	if (ncols == 0) {
	    /* no data to print */
	    continue;
	}

	if (!started_data_lines) {
	    x[0] = spec->data;
	    /* see below for subsequent adjustment of x[1] */
	    x[1] = x[0] + spec->nobs;
	    started_data_lines = 1;
	} 

	x[2] = x[1] + spec->nobs;
	x[3] = x[2] + spec->nobs;
	x[4] = x[3] + spec->nobs;

	for (t=0; t<spec->nobs; t++) {
	    /* print x-axis value */
	    if (na(x[0][t])) {
		fputs("? ", fp);
		miss = 1;
	    } else {
		fprintf(fp, "%.8g ", x[0][t]);
	    }

	    /* conversion, if needed, between (y, ydelta) and
	       (ylow, yhigh)
	     */
	    if (spec->lines[i].ncols == 3) {
		if (spec->flags & GPT_FILL_SWITCH) { 
		    yt = x[1][t];
		    et = x[2][t];
		    if (!na(yt) && !na(et)) {
			x[1][t] = yt - et;
			x[2][t] = yt + et;
		    }
		} else if (spec->flags & GPT_ERR_SWITCH) {
		    if (!na(x[1][t]) && !na(x[2][t])) {
			et = (x[2][t] - x[1][t]) / 2.0;
			yt = x[1][t] + et;
			x[1][t] = yt;
			x[2][t] = et;
		    }
		}
	    } 

	    /* print y-axis value(s) */
	    for (j=1; j<ncols; j++) {
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

	x[1] += (ncols - 1) * spec->nobs;
    }

    gretl_pop_c_numeric_locale();

    if (png && gnuplot_has_bbox()) {
	print_plot_bounding_box_request(fp);
    }

    return miss;
}

static int set_loess_fit (GPT_SPEC *spec, int d, double q, gretl_matrix *x,
			  gretl_matrix *y, gretl_matrix *yh)
{
    int t, T = gretl_vector_get_length(y);
    double *data;

    data = realloc(spec->data, 3 * T * sizeof *data);
    if (data == NULL) {
	return E_ALLOC;
    }

    for (t=0; t<T; t++) {
	data[t] = x->val[t];
	data[t+T] = y->val[t];
	data[t+2*T] = yh->val[t];
    }

    spec->data = data;
    spec->nobs = spec->okobs = T;

    sprintf(spec->lines[1].title, G_("loess fit, d = %d, q = %g"), d, q);
    strcpy(spec->lines[1].scale, "1.0");
    strcpy(spec->lines[1].style, "lines");
    spec->lines[1].ncols = 2;

    spec->fit = PLOT_FIT_LOESS;

    return 0;
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
    } else if (f == PLOT_FIT_INVERSE) {
	b = spec->b_inv->val;
	sprintf(title, "Y = %#.3g %c %#.3g(1/X)", b[0],
		(b[1] > 0)? '+' : '-', fabs(b[1]));
	gretl_push_c_numeric_locale();
	sprintf(formula, "%g + %g/x", b[0], b[1]);
	gretl_pop_c_numeric_locale();
    }

    strcpy(spec->lines[1].scale, "NA");
    strcpy(spec->lines[1].style, "lines");
    spec->lines[1].ncols = 0;

    spec->fit = f;
}

int plotspec_add_fit (GPT_SPEC *spec, FitType f)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *yh = NULL;
    const double *py;
    const double *px;
    int T = spec->okobs;
    double q = 0.5;
    int d = 1;
    int i, t, k;
    int err = 0;

    if ((f == PLOT_FIT_OLS && spec->b_ols != NULL) ||
	(f == PLOT_FIT_QUADRATIC && spec->b_quad != NULL) ||
	(f == PLOT_FIT_INVERSE && spec->b_inv != NULL)) {
	/* just activate existing setup */
	set_fitted_line(spec, f);
	return 0;
    }
	
    if (f == PLOT_FIT_OLS || f == PLOT_FIT_INVERSE) {
	k = 2;
    } else if (f == PLOT_FIT_QUADRATIC) {
	k = 3;
    } else if (f == PLOT_FIT_LOESS) {
	k = 1;
    } else {
	return E_DATA;
    }

    y = gretl_column_vector_alloc(T);
    X = gretl_matrix_alloc(T, k);
    if (y == NULL || X == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (f != PLOT_FIT_LOESS) {
	b = gretl_column_vector_alloc(k);
	if (b == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
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
		if (f == PLOT_FIT_INVERSE) {
		    gretl_matrix_set(X, i, 1, 1.0 / px[t]);
		} else {
		    gretl_matrix_set(X, i, 1, px[t]);
		}
		if (f == PLOT_FIT_QUADRATIC) {
		    gretl_matrix_set(X, i, 2, px[t] * px[t]);
		}
	    }
	    i++;
	}
    }

    if (f == PLOT_FIT_LOESS) {
	err = sort_pairs_by_x(X, y, NULL, spec->markers);
	if (!err) {
	    yh = loess_fit(X, y, d, q, OPT_R, &err);
	}
    } else {
	err = gretl_matrix_ols(y, X, b, NULL, NULL, NULL);
    }

    if (!err && spec->n_lines == 1) {
	err = plotspec_add_line(spec);
    }

    if (!err) {
	if (f == PLOT_FIT_OLS) {
	    spec->b_ols = b;
	    b = NULL;
	} else if (f == PLOT_FIT_QUADRATIC) {
	    spec->b_quad = b;
	    b = NULL;
	} else if (f == PLOT_FIT_INVERSE) {
	    spec->b_inv = b;
	    b = NULL;
	}	    
    }

    if (!err) {
	if (f == PLOT_FIT_LOESS) {
	    set_loess_fit(spec, d, q, X, y, yh);
	} else {
	    set_fitted_line(spec, f);
	}
    }

 bailout:
    
    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(yh);

    return err;
}
