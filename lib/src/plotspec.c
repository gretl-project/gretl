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

struct plotbars_ {
    int n;         /* number of available start-stop "bar" pairs */
    double t1;     /* current min. value on time axis */
    double t2;     /* current max. value on time axis */
    double ymin;   /* current min. value on y axis */
    double ymax;   /* current max. value on y axis */
    double **dx;   /* start-stop time pairs */
};

/* for handling "recession bars" or similar */
static int n_bars_shown (double xmin, double xmax, plotbars *bars);
static void print_bars_header (GPT_SPEC *spec, FILE *fp);
static void print_plotbars (GPT_SPEC *spec, FILE *fp);

static void print_filledcurve_color (FILE *fp);

GPT_SPEC *plotspec_new (void)
{
    GPT_SPEC *spec;
    int i;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->datacols = 0;
    spec->heredata = 0;

    spec->lines = NULL;
    spec->n_lines = 0;

    spec->labels = NULL;
    spec->n_labels = 0;

    spec->arrows = NULL;
    spec->n_arrows = 0;

    spec->literal = NULL;
    spec->n_literal = 0;

    spec->multi_xtics = NULL;
    spec->n_xtics = 0;

    for (i=0; i<5; i++) {
	spec->titles[i] = NULL;
    }

    *spec->xvarname = '\0';
    *spec->yvarname = '\0';

    spec->xticstr = NULL;
    spec->x2ticstr = NULL;

    spec->xfmt[0] = 0;
    spec->xtics[0] = 0;
    spec->mxtics[0] = 0;
    spec->yfmt[0] = 0;
    spec->ytics[0] = 0;
    spec->fname[0] = 0;
    spec->keyspec = GP_KEY_LEFT_TOP;

    for (i=0; i<5; i++) {
	spec->range[i][0] = NADBL;
	spec->range[i][1] = NADBL;
	if (i < 3) {
	    spec->logbase[i] = 0.0;
	}
    }

    spec->b_ols = NULL;
    spec->b_quad = NULL;
    spec->b_cub = NULL;
    spec->b_inv = NULL;
    spec->b_log = NULL;
    spec->b_linlog = NULL;

    spec->code = PLOT_REGULAR;
    spec->flags = 0;
    spec->fit = PLOT_FIT_NONE;
    spec->fp = NULL;
    spec->data = NULL;
    spec->auxdata = NULL;
    spec->markers = NULL;
    spec->n_markers = 0;
    spec->scale = 1.0;
    spec->labeled = NULL;
    spec->ptr = NULL;
    spec->bars = NULL;
    spec->fontstr = NULL;
    spec->reglist = NULL;
    spec->nobs = 0;
    spec->okobs = 0;
    spec->pd = 0;
    spec->nbars = 0;
    spec->boxwidth = 0;
    spec->fillfrac = 0;
    spec->samples = 0;
    spec->border = GP_BORDER_DEFAULT;
    spec->border_lc[0] = '\0';
    spec->bmargin = 0;

    spec->termtype = GP_TERM_NONE;

    return spec;
}

static void plotbars_free (plotbars *bars)
{
    if (bars != NULL) {
	doubles_array_free(bars->dx, bars->n);
	free(bars);
    }
}

static void free_plotspec_lines (GPT_SPEC *spec)
{
    int i;

    for (i=0; i<spec->n_lines; i++) {
	free(spec->lines[i].ustr);
	free(spec->lines[i].mcols);
    }

    free(spec->lines);
}

void plotspec_destroy (GPT_SPEC *spec)
{
    if (spec == NULL) {
	return;
    }

    if (spec->lines != NULL) {
	free_plotspec_lines(spec);
    }

    if (spec->labels != NULL) {
	free(spec->labels);
    }

    if (spec->arrows != NULL) {
	free(spec->arrows);
    }

    if (spec->data != NULL) {
	gretl_matrix_free(spec->data);
    }

    if (spec->auxdata != NULL) {
	gretl_matrix_free(spec->auxdata);
    }

    if (spec->reglist != NULL) {
	free(spec->reglist);
    }

    if (spec->literal != NULL) {
	strings_array_free(spec->literal, spec->n_literal);
    }

    if (spec->multi_xtics != NULL) {
	strings_array_free(spec->multi_xtics, spec->n_xtics);
    }

    if (spec->markers != NULL) {
	strings_array_free(spec->markers, spec->n_markers);
    }

    if (spec->labeled != NULL) {
	free(spec->labeled);
    }

    if (spec->bars != NULL) {
	plotbars_free(spec->bars);
    }

    if (spec->fontstr != NULL) {
	free(spec->fontstr);
    }

    if (spec->xticstr != NULL) {
	free(spec->xticstr);
    }

    if (spec->x2ticstr != NULL) {
	free(spec->x2ticstr);
    }

    gretl_matrix_free(spec->b_ols);
    gretl_matrix_free(spec->b_quad);
    gretl_matrix_free(spec->b_cub);
    gretl_matrix_free(spec->b_inv);
    gretl_matrix_free(spec->b_log);
    gretl_matrix_free(spec->b_linlog);

    free(spec);
}

static void plotspec_get_xt1_xt2 (const GPT_SPEC *spec,
				  double *xt1,
				  double *xt2)
{
    if (spec->data != NULL) {
	*xt1 = spec->data->val[0];
	*xt2 = spec->data->val[spec->nobs - 1];
    }
}

static gp_style_spec style_specs[] = {
    { GP_STYLE_LINES,        "lines",        N_("lines") },
    { GP_STYLE_POINTS,       "points",       N_("points") },
    { GP_STYLE_LINESPOINTS,  "linespoints",  N_("lines/points") },
    { GP_STYLE_IMPULSES,     "impulses",     N_("impulses") },
    { GP_STYLE_DOTS,         "dots",         N_("dots") },
    { GP_STYLE_STEPS,        "steps",        N_("steps") },
    { GP_STYLE_BOXES,        "boxes",        N_("boxes") },
    { GP_STYLE_ERRORBARS,    "errorbars",    N_("error bars") },
    { GP_STYLE_FILLEDCURVE,  "filledcurve",  N_("filled curve") },
    { GP_STYLE_CANDLESTICKS, "candlesticks", N_("candlesticks") },
    { 0, NULL, NULL }
};

const char *gp_line_style_display_name (int t)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
	if (t == style_specs[i].id) {
	    return style_specs[i].trname;
	}
    }

    return N_("lines");
}

static const char *gp_line_style_name (int t)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
	if (t == style_specs[i].id) {
	    return style_specs[i].name;
	}
    }

    return "lines";
}

int gp_style_index_from_name (const char *s)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
	if (!strcmp(s, style_specs[i].name)) {
	    return style_specs[i].id;
	}
    }

    /* allowed abbreviations */

    if (!strcmp(s, "l")) {
	return GP_STYLE_LINES;
    } else if (!strcmp(s, "p")) {
	return GP_STYLE_POINTS;
    } else if (!strcmp(s, "lp")) {
	return GP_STYLE_LINESPOINTS;
    } else if (!strcmp(s, "i")) {
	return GP_STYLE_IMPULSES;
    } else if (!strcmp(s, "b")) {
	return GP_STYLE_BOXES;
    } else if (!strncmp(s, "can", 3)) {
	return GP_STYLE_CANDLESTICKS;
    }

    /* fallback */
    return GP_STYLE_AUTO;
}

int gp_style_index_from_display_name (const char *s)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
	if (!strcmp(s, _(style_specs[i].trname))) {
	    return style_specs[i].id;
	}
    }

    return GP_STYLE_LINES;
}

gp_style_spec *get_style_spec (int t)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
	if (t == style_specs[i].id) {
	    return &style_specs[i];
	}
    }

    return NULL;
}

static gp_key_spec key_specs[] = {
    { GP_KEY_LEFT_TOP,     N_("left top") },
    { GP_KEY_RIGHT_TOP,    N_("right top") },
    { GP_KEY_LEFT_BOTTOM,  N_("left bottom") },
    { GP_KEY_RIGHT_BOTTOM, N_("right bottom") },
    { GP_KEY_OUTSIDE,      N_("outside") },
    { GP_KEY_NONE,         N_("none") },
    /* in reverse order */
    { GP_KEY_LEFT_TOP,     N_("top left") },
    { GP_KEY_RIGHT_TOP,    N_("top right") },
    { GP_KEY_LEFT_BOTTOM,  N_("bottom left") },
    { GP_KEY_RIGHT_BOTTOM, N_("bottom right") },
    { -1,                  NULL }
};

static const char *gp_keypos_string (int t)
{
    int i;

    for (i=0; key_specs[i].str != NULL; i++) {
	if (t == key_specs[i].id) {
	    return key_specs[i].str;
	}
    }

    return N_("none");
}

void print_keypos_string (int t, FILE *fp)
{
    const char *s = gp_keypos_string(t);

    if (!strcmp(s, "none")) {
	fputs("set nokey\n", fp);
    } else {
	fprintf(fp, "set key %s\n", s);
    }
}

int gp_keypos_from_name (const char *s)
{
    int i;

    for (i=0; key_specs[i].id >= 0; i++) {
	if (!strcmp(s, key_specs[i].str)) {
	    return key_specs[i].id;
	}
    }

    return GP_KEY_NONE;
}

int gp_keypos_from_display_name (const char *s)
{
    int i;

    for (i=0; key_specs[i].id >= 0; i++) {
	if (!strcmp(s, _(key_specs[i].str))) {
	    return key_specs[i].id;
	}
    }

    return GP_KEY_NONE;
}

gp_key_spec *get_keypos_spec (int t)
{
    int i;

    for (i=0; key_specs[i].id >= 0; i++) {
	if (t == key_specs[i].id) {
	    return &key_specs[i];
	}
    }

    return NULL;
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
    lines[n].style = 0;
    lines[n].pscale = 1.0;
    lines[n].rgb[0] = '\0';
    lines[n].yaxis = 1;
    lines[n].type = LT_AUTO;
    lines[n].ptype = 0;
    lines[n].dtype = 0;
    lines[n].width = 1.0;
    lines[n].ncols = 0;
    lines[n].whiskwidth = 0;
    lines[n].flags = 0;

    lines[n].title = NULL;
    lines[n].formula = NULL;
    lines[n].ustr = NULL;
    lines[n].mcols = NULL;

    return 0;
}

static void copy_line_content (GPT_LINE *targ, GPT_LINE *src)
{
    targ->varnum = src->varnum;
    targ->style = src->style;
    targ->pscale = src->pscale;
    strcpy(targ->rgb, src->rgb);
    targ->yaxis = src->yaxis;
    targ->type = src->type;
    targ->ptype = src->ptype;
    targ->dtype = src->dtype;
    targ->width = src->width;
    targ->ncols = src->ncols;
    targ->whiskwidth = src->whiskwidth;
    targ->flags = src->flags;

    if (src->title != NULL) {
	g_free(targ->title);
	targ->title = g_strdup(src->title);
    }
    if (src->formula != NULL) {
	g_free(targ->formula);
	targ->formula = g_strdup(src->formula);
    }
    if (src->ustr != NULL) {
	free(targ->ustr);
	targ->ustr = gretl_strdup(src->ustr);
    }
    if (src->mcols != NULL) {
	free(targ->mcols);
	targ->mcols = gretl_list_copy(src->mcols);
    }
}

int plotspec_delete_line (GPT_SPEC *spec, int i)
{
    GPT_LINE *lines = spec->lines;
    int j, n = spec->n_lines;

    if (i < 0 || i >= n) {
	return E_DATA;
    }

    g_free(lines[i].title);
    g_free(lines[i].formula);
    free(lines[i].ustr);
    free(lines[i].mcols);

    /* move the remaining line contents down */
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

    lines = calloc(spec->n_lines, sizeof *lines);

    if (lines == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<spec->n_lines; i++) {
	    lines[i].ustr = NULL;
	    lines[i].mcols = NULL;
	    copy_line_content(&lines[i], &spec->lines[i]);
	}
    }

    return lines;
}

int plotspec_max_line_width (GPT_SPEC *spec)
{
    int i, lw = 0;

    for (i=0; i<spec->n_lines; i++) {
	if (spec->lines[i].width > lw) {
	    lw = spec->lines[i].width;
	}
    }

    return lw;
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

int plotspec_add_arrow (GPT_SPEC *spec)
{
    GPT_ARROW *arrows;
    int n = spec->n_arrows;

    arrows = realloc(spec->arrows, (n + 1) * sizeof *arrows);
    if (arrows == NULL) {
	return E_ALLOC;
    }

    spec->arrows = arrows;
    spec->n_arrows += 1;
    arrows[n].x0 = 0;
    arrows[n].y0 = 0;
    arrows[n].x1 = 0;
    arrows[n].y1 = 0;
    arrows[n].flags = 0;

    return 0;
}

static void copy_arrow_content (GPT_ARROW *targ, GPT_ARROW *src)
{
    targ->x0 = src->x0;
    targ->y0 = src->y0;
    targ->x1 = src->x1;
    targ->y1 = src->y1;
    targ->flags = src->flags;
}

int plotspec_delete_arrow (GPT_SPEC *spec, int i)
{
    GPT_ARROW *arrows = spec->arrows;
    int j, n = spec->n_arrows;
    int err = 0;

    if (i < 0 || i >= n) {
	return E_DATA;
    }

    for (j=i; j<n-1; j++) {
	copy_arrow_content(&arrows[j], &arrows[j+1]);
    }

    spec->n_arrows -= 1;

    if (spec->n_arrows == 0) {
	free(spec->arrows);
	spec->arrows = NULL;
    } else {
	arrows = realloc(spec->arrows, (n - 1) * sizeof *arrows);
	if (arrows == NULL) {
	    err = E_ALLOC;
	} else {
	    spec->arrows = arrows;
	}
    }

    return err;
}

GPT_ARROW *plotspec_clone_arrows (GPT_SPEC *spec, int *err)
{
    GPT_ARROW *arrows = NULL;
    int i;

    if (spec->n_arrows == 0) {
	/* no-op */
	return NULL;
    }

    arrows = malloc(spec->n_arrows * sizeof *arrows);
    if (arrows == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<spec->n_arrows; i++) {
	    copy_arrow_content(&arrows[i], &spec->arrows[i]);
	}
    }

    return arrows;
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
	return N_("left");
    } else if (j == GP_JUST_CENTER) {
	return N_("center");
    } else if (j == GP_JUST_RIGHT) {
	return N_("right");
    } else {
	return N_("left");
    }
}

static void print_plot_labelspec (const GPT_LABEL *lbl, FILE *fp)
{
    char *label = escape_quotes(lbl->text);

    fprintf(fp, "set label \"%s\" ", (label != NULL)?
	    label : lbl->text);

    gretl_push_c_numeric_locale();

    fprintf(fp, "at %g,%g %s front\n",
	    lbl->pos[0], lbl->pos[1],
	    gp_justification_string(lbl->just));

    gretl_pop_c_numeric_locale();

    if (label != NULL) {
	free(label);
    }
}

static void print_plot_arrow (const GPT_ARROW *arrow, FILE *fp)
{
    gretl_push_c_numeric_locale();
    fprintf(fp, "set arrow from %g,%g to %g,%g %s",
	    arrow->x0, arrow->y0, arrow->x1, arrow->y1,
	    (arrow->flags & GP_ARROW_HEAD)? "head" : "nohead");
    if (arrow->flags & GP_ARROW_DOTS) {
	fputs(" lt 0", fp);
    }
    fputc('\n', fp);
    gretl_pop_c_numeric_locale();
}

static int print_polar_labels (const GPT_SPEC *spec, FILE *fp)
{
    double x, y;
    int t;

    fputs("# printing data labels\n", fp);

    gretl_push_c_numeric_locale();

    for (t=0; t<spec->nobs; t++) {
	if (spec->labeled != NULL && !spec->labeled[t]) {
	    /* printing only specified labels */
	    continue;
	}
	if (sscanf(spec->markers[t], "%lf,%lf", &x, &y) == 2) {
	    fprintf(fp, "set label \"%.3f,%.3f\" at %.3f,%.3f\n",
		    x, y, x, y + .04);
	}
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

static int usable_obs (const double *x, const double *y0,
		       const double *y1, int t,
		       const double **py)
{
    *py = y0;

    if (!na(x[t]) && !na(y0[t])) {
	return 1;
    } else if (!na(x[t]) && y1 != NULL && !na(y1[t])) {
	*py = y1;
	return 1;
    } else {
	return 0;
    }
}

int plotspec_line_is_formula (const GPT_SPEC *spec, int i)
{
    int ret = 0;

    if (i >= 0 && i < spec->n_lines) {
	GPT_LINE *line = &spec->lines[i];

	if (!string_is_blank(line->formula)) {
	    ret = 1;
	} else if (line->flags & GP_LINE_USER) {
	    ret = 1;
	} else if (i == 1 && (spec->flags & GPT_AUTO_FIT)) {
	    ret = 1;
	}
    }

    return ret;
}

static int print_data_labels (const GPT_SPEC *spec, FILE *fp)
{
    const double *x, *y, *y0;
    const double *y1 = NULL;
    double xrange, yrange;
    double xmin, xmax;
    double ymin, ymax;
    double yoff;
    int nlines;
    int i, t;

    if (spec->markers == NULL) {
	/* "can't happen" */
	return 1;
    }

    if (spec->flags & GPT_POLAR) {
	return print_polar_labels(spec, fp);
    }

    nlines = 0;
    for (i=0; i<spec->n_lines; i++) {
	if (!plotspec_line_is_formula(spec, i)) {
	    nlines++;
	}
    }
    if (nlines > 2 || spec->lines[0].ncols != 2) {
	return 1;
    }

    x = spec->data->val;
    y0 = x + spec->nobs;

    if (spec->code == PLOT_FACTORIZED) {
	y1 = y0 + spec->nobs;
    }

    xmin = spec->range[GP_X_RANGE][0];
    xmax = spec->range[GP_X_RANGE][1];
    xrange = (na(xmin) || na(xmax)) ? 0.0 : xmax - xmin;

    ymin = spec->range[GP_Y_RANGE][0];
    ymax = spec->range[GP_Y_RANGE][1];
    yrange = (na(ymin) || na(ymax)) ? 0.0 : ymax - ymin;

    if (xrange == 0.0 || yrange == 0.0) {
	ymin = 1.0e+50; ymax = -1.0e+50;
	xmin = 1.0e+50; xmax = -1.0e+50;

	for (t=0; t<spec->nobs; t++) {
	    if (usable_obs(x, y0, y1, t, &y)) {
		if (yrange == 0.0) {
		    if (y[t] < ymin) {
			ymin = y[t];
		    }
		    if (y[t] > ymax) {
			ymax = y[t];
		    }
		}
		if (xrange == 0.0) {
		    if (x[t] < xmin) {
			xmin = x[t];
		    }
		    if (x[t] > xmax) {
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

	if (usable_obs(x, y0, y1, t, &y)) {
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
    } else if (fit == PLOT_FIT_CUBIC) {
	fputs("# plot includes automatic fit: cubic\n", fp);
    } else if (fit == PLOT_FIT_INVERSE) {
	fputs("# plot includes automatic fit: inverse\n", fp);
    } else if (fit == PLOT_FIT_LOESS) {
	fputs("# plot includes automatic fit: loess\n", fp);
    } else if (fit == PLOT_FIT_LOGLIN) {
	fputs("# plot includes automatic fit: semilog\n", fp);
    } else if (fit == PLOT_FIT_LINLOG) {
	fputs("# plot includes automatic fit: linlog\n", fp);
    }
}

void print_plot_ranges_etc (const GPT_SPEC *spec, FILE *fp)
{
    const char *rstrs[] = {
	"x", "y", "y2", "t", "x2"
    };
    int i;

    gretl_push_c_numeric_locale();

    for (i=0; i<5; i++) {
	if (i < 3 && spec->logbase[i] > 0.0) {
	    fprintf(fp, "set logscale %s %g\n", rstrs[i], spec->logbase[i]);
	}

	if (na(spec->range[i][0]) || na(spec->range[i][1]) ||
	    spec->range[i][0] == spec->range[i][1]) {
	    continue;
	}

	if ((i == GP_Y2_RANGE && !(spec->flags & GPT_Y2AXIS)) ||
	    (i == GP_T_RANGE && !(spec->flags & GPT_PARAMETRIC))) {
	    continue;
	}

	fprintf(fp, "set %srange [%.10g:%.10g]\n", rstrs[i],
		spec->range[i][0], spec->range[i][1]);

	if (i == 4 && spec->code == PLOT_PROB_DIST && spec->samples == 0) {
	    int ns = spec->range[i][1] - spec->range[i][0] + 1;

	    fprintf(fp, "set samples %d\n", ns);
	}
    }

    if (spec->boxwidth > 0 && spec->boxwidth < 1) {
	fprintf(fp, "set boxwidth %g relative\n", (double) spec->boxwidth);
    } else if (spec->boxwidth < 0 && spec->boxwidth > -1) {
	fprintf(fp, "set boxwidth %g absolute\n", (double) -spec->boxwidth);
    }

    gretl_pop_c_numeric_locale();
}

static int blank_user_line (GPT_LINE *line)
{
    return ((line->flags & GP_LINE_USER) &&
	    string_is_blank(line->formula));
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
	    if (j != skipline && !blank_user_line(&spec->lines[j])) {
		return 1;
	    }
	}
	return 0;
    }
}

static void print_linestyle (const GPT_SPEC *spec,
			     int targ, int src,
			     FILE *fp)
{
    GPT_LINE *line;
    int done = 0;

    fprintf(fp, "set linetype %d ", targ + 1);

    if (src < spec->n_lines) {
	line = &spec->lines[src];
	if (*line->rgb != '\0' && line->style != GP_STYLE_FILLEDCURVE) {
	    fprintf(fp, "lc rgb \"%s\"\n", line->rgb);
	    done = 1;
	}
    }

    if (!done) {
	/* hmm, should we be doing this? */
	gretlRGB color = get_graph_color(targ);
	char cstr[10];

	print_rgb_hash(cstr, color);
	fprintf(fp, "lc rgb \"%s\"\n", cstr);
    }
}

static void write_styles_from_plotspec (const GPT_SPEC *spec, FILE *fp)
{
    char cstr[8];
    int i;

    if (frequency_plot_code(spec->code)) {
	gretlRGB color = get_boxcolor();

	print_rgb_hash(cstr, color);
	fprintf(fp, "set linetype 1 lc rgb \"%s\"\n", cstr);
	fputs("set linetype 2 lc rgb \"#000000\"\n", fp);
    } else if (spec->code == PLOT_RQ_TAU) {
	fputs("set linetype 1 lc rgb \"#000000\"\n", fp);
	for (i=1; i<N_GP_LINETYPES; i++) {
	    print_linestyle(spec, i, i, fp);
	}
    } else if (spec->code == PLOT_BOXPLOTS) {
	/* the first two elements use the same color */
	print_linestyle(spec, 0, 0, fp);
	print_linestyle(spec, 1, 2, fp);
    } else {
	/* note: handle the case where "line 0" is filledcurve,
	   the color of which is handled separately
	*/
	int offset = 0;

	if (spec->lines[0].style == GP_STYLE_FILLEDCURVE) {
	    offset = 1;
	}
	for (i=0; i<N_GP_LINETYPES; i++) {
	    print_linestyle(spec, i, i + offset, fp);
	}
    }
}

static void maybe_print_point_info (const GPT_LINE *line, FILE *fp)
{
    if (line->style == GP_STYLE_POINTS ||
	line->style == GP_STYLE_LINESPOINTS) {
	if (line->ptype != 0) {
	    fprintf(fp, " pt %d", line->ptype);
	}
	if (line->pscale != 1.0) {
	    fprintf(fp, " ps %g", (double) line->pscale);
	}
    }

    if (line->style == GP_STYLE_LINES ||
	line->style == GP_STYLE_LINESPOINTS) {
	if (line->dtype > 1) {
	    fprintf(fp, " dt %d", line->dtype);
	}
    }
}

#define show_fit(s) (s->fit == PLOT_FIT_OLS ||		\
                     s->fit == PLOT_FIT_QUADRATIC ||	\
		     s->fit == PLOT_FIT_CUBIC ||	\
                     s->fit == PLOT_FIT_INVERSE ||	\
                     s->fit == PLOT_FIT_LOESS ||	\
		     s->fit == PLOT_FIT_LOGLIN ||	\
		     s->fit == PLOT_FIT_LINLOG)

int gp_line_data_columns (GPT_SPEC *spec, int i)
{
    if (spec->lines[i].flags & GP_LINE_AUXDATA) {
	return 0;
    } else {
	return spec->lines[i].ncols;
    }
}

static void plotspec_print_data (GPT_SPEC *spec,
				 int skipline,
				 int *miss,
				 FILE *fp)
{
    double *x[5];
    int started_data_lines = 0;
    int i, j, t;

    for (i=0; i<spec->n_lines; i++) {
	int ncols = gp_line_data_columns(spec, i);

	if (ncols == 0 || i == skipline) {
	    /* no (regular) data to print */
	    continue;
	}

	if (!started_data_lines) {
	    x[0] = spec->data->val;
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
		write_gp_dataval(NADBL, fp, 0);
		*miss = 1;
	    } else if (spec->flags & GPT_TIMEFMT) {
		fprintf(fp, "%.0f ", x[0][t]);
	    } else {
		fprintf(fp, "%.10g ", x[0][t]);
	    }
	    /* print y-axis value(s) */
	    for (j=1; j<ncols; j++) {
		if (na(x[j][t])) {
		    *miss = 1;
		}
		write_gp_dataval(x[j][t], fp, 0);
	    }
	    /* append markers if applicable */
	    if (spec->markers != NULL && i == 0) {
		fprintf(fp, " # %s", spec->markers[t]);
	    }
	    fputc('\n', fp);
	}

	fputs("e\n", fp);

	x[1] += (ncols - 1) * spec->nobs;
    }

    if (spec->auxdata != NULL) {
	int rows = gretl_matrix_rows(spec->auxdata);
	int cols = gretl_matrix_cols(spec->auxdata);

	fprintf(fp, "# auxdata %d %d\n", rows, cols);

	for (i=0; i<rows; i++) {
	    for (j=0; j<cols; j++) {
		fprintf(fp, "%.10g ", gretl_matrix_get(spec->auxdata, i, j));
	    }
	    fputc('\n', fp);
	}
	fputs("e\n", fp);
    }
}

static void plotspec_print_heredata (GPT_SPEC *spec,
				     int *miss,
				     FILE *fp)
{
    gretl_matrix *m = NULL;
    double mti;
    int i, t;

    if (spec->auxdata != NULL) {
	int r = spec->auxdata->rows;
	int c = spec->auxdata->cols;
	int j;

	fprintf(fp, "# auxdata %d %d\n", r, c);
	fputs("$aux << EOA\n", fp);

	for (i=0; i<r; i++) {
	    for (j=0; j<c; j++) {
		fprintf(fp, "%.10g ", gretl_matrix_get(spec->auxdata, i, j));
	    }
	    fputc('\n', fp);
	}
	fputs("EOA\n", fp);
    }

    if (spec->data == NULL || spec->datacols == 0) {
	return;
    }

    m = spec->data;

    fputs("$data << EOD\n", fp);

    for (t=0; t<spec->nobs; t++) {
	for (i=0; i<m->cols; i++) {
	    mti = gretl_matrix_get(m, t, i);
	    if (na(mti)) {
		write_gp_dataval(NADBL, fp, 0);
		*miss = 1;
	    } else if (i == 0 && (spec->flags & GPT_TIMEFMT)) {
		fprintf(fp, "%.0f ", mti);
	    } else {
		fprintf(fp, "%.10g ", mti);
	    }
	}
	if (spec->markers != NULL) {
	    fprintf(fp, " # %s", spec->markers[t]);
	}
	fputc('\n', fp);
    }

    fputs("EOD\n", fp);
}

int plotspec_print (GPT_SPEC *spec, FILE *fp)
{
    int png = get_png_output(spec);
    int mono = (spec->flags & GPT_MONO);
    double xt1 = 0, xt2 = 0;
    char src[16] = "'-'";
    int skipline = -1;
    int ycol, any_y2 = 0;
    int anydata = 0;
    int i, k, miss = 0;

    if (spec->pd > 0) {
	fprintf(fp, "# timeseries %d", spec->pd);
	if (spec->flags & GPT_LETTERBOX) {
	    fputs(" (letterbox)\n", fp);
	} else {
	    fputc('\n', fp);
	}
    }

    if (spec->scale != 1.0) {
	gretl_push_c_numeric_locale();
	fprintf(fp, "# scale = %.1f\n", spec->scale);
	gretl_pop_c_numeric_locale();
    }

    if (mono) {
	fputs("set mono\n", fp);
    } else {
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

    if (!string_is_blank(spec->titles[4])) {
	fprintf(fp, "set x2label \"%s\"\n", spec->titles[4]);
    }

    for (i=0; i<spec->n_labels; i++) {
	if (!string_is_blank(spec->labels[i].text)) {
	    print_plot_labelspec(&spec->labels[i], fp);
	}
    }

    for (i=0; i<spec->n_arrows; i++) {
	print_plot_arrow(&spec->arrows[i], fp);
    }

    if (spec->flags & GPT_XZEROAXIS) {
	fputs("set xzeroaxis\n", fp);
    }

    if (spec->flags & GPT_YZEROAXIS) {
	fputs("set yzeroaxis\n", fp);
    }

    if (spec->flags & (GPT_GRID_Y | GPT_GRID_X)) {
	if (!(spec->flags & GPT_GRID_X)) {
	    fputs("set grid ytics\n", fp);
	} else if (!(spec->flags & GPT_GRID_Y)) {
	    fputs("set grid xtics\n", fp);
	} else {
	    /* default, both */
	    fputs("set grid\n", fp);
	}
    }

    if (spec->flags & GPT_PARAMETRIC) {
	fputs("set parametric\n", fp);
	if (spec->samples > 0) {
	    fprintf(fp, "set samples %d\n", spec->samples);
	}
    }

    if (spec->flags & GPT_POLAR) {
	fputs("unset ytics\n", fp);
	fputs("set size square\n", fp);
	fputs("set polar\n", fp);
    }

    gnuplot_missval_string(fp);

    if (spec->keyspec == GP_KEY_NONE) {
	fputs("set nokey\n", fp);
    } else {
	fprintf(fp, "set key %s\n", gp_keypos_string(spec->keyspec));
    }

    print_plot_ranges_etc(spec, fp);

    /* using time format for x-axis? */
    if (spec->flags & GPT_TIMEFMT) {
	fputs("set xdata time\n", fp);
	fputs("set timefmt \"%s\"\n", fp);
    }

    /* special x and/or y format? */
    if (*spec->xfmt != '\0') {
	fprintf(fp, "set format x \"%s\"\n", spec->xfmt);
    }
    if (*spec->yfmt != '\0') {
	fprintf(fp, "set format y \"%s\"\n", spec->yfmt);
    }

    /* customized xtics? */
    if (!strcmp(spec->xtics, "none")) {
	fputs("set noxtics\n", fp);
    } else {
	if (!string_is_blank(spec->xtics)) {
	    /* may contain "nomirror" */
	    fprintf(fp, "set xtics %s\n", spec->xtics);
	}
	if (spec->xticstr != NULL) {
	    fprintf(fp, "set xtics %s\n", spec->xticstr);
	}
	if (!string_is_blank(spec->mxtics)) {
	    fprintf(fp, "set mxtics %s\n", spec->mxtics);
	}
    }

    /* using x2tics? */
    if (spec->x2ticstr != NULL) {
	fprintf(fp, "set x2tics %s\n", spec->x2ticstr);
    }

    /* customized ytics? */
    if (!strcmp(spec->ytics, "none")) {
	fputs("set noytics\n", fp);
    } else {
	if (!string_is_blank(spec->ytics)) {
	    fprintf(fp, "set ytics %s\n", spec->ytics);
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
	} else if (spec->border_lc[0] != '\0') {
	    fprintf(fp, "set border %d lc rgb \"%s\"\n",
		    spec->border, spec->border_lc);
	} else {
	    fprintf(fp, "set border %d\n", spec->border);
	}
	if (string_is_blank(spec->xtics)) {
	    fputs("set xtics nomirror\n", fp);
	}
	if (string_is_blank(spec->ytics)) {
	    fputs("set ytics nomirror\n", fp);
	}
    }

    if (spec->bmargin > 0) {
	fprintf(fp, "set bmargin %d\n", spec->bmargin);
    }

    /* In case of plots that are editable (see gui client), it is
       important to write out the comment string that identifies the
       sort of graph, so that it will be recognized by type when
       it is redisplayed.
    */

    write_plot_type_string(spec->code, spec->flags, fp);

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

    if (spec->n_xtics > 0) {
	fprintf(fp, "# xtics lines = %d\n", spec->n_xtics);
	for (i=0; i<spec->n_xtics; i++) {
	    if (spec->multi_xtics[i] != NULL && *spec->multi_xtics[i] != '\0') {
		fprintf(fp, "%s\n", spec->multi_xtics[i]);
	    }
	}
    }

    if (show_fit(spec)) {
	print_auto_fit_string(spec->fit, fp);
    }

    if (spec->bars != NULL) {
	plotspec_get_xt1_xt2(spec, &xt1, &xt2);
	spec->nbars = n_bars_shown(xt1, xt2, spec->bars);
	if (spec->nbars > 0) {
	    fprintf(fp, "# n_bars = %d\n", spec->nbars);
	}
    }

    print_user_lines_info(spec, fp);

    if ((spec->code == PLOT_FREQ_SIMPLE ||
	 spec->code == PLOT_FREQ_NORMAL ||
	 spec->code == PLOT_FREQ_GAMMA)) {
	if (mono) {
	    fputs("set style fill solid 0.3\n", fp);
	} else {
	    fputs("set style fill solid 0.6\n", fp);
	}
    } else if (spec->fillfrac > 0) {
	gretl_push_c_numeric_locale();
	fprintf(fp, "set style fill solid %g\n", (double) spec->fillfrac);
	gretl_pop_c_numeric_locale();
    }

    if (spec->flags & GPT_PRINT_MARKERS) {
	print_data_labels(spec, fp);
    }

    if (spec->flags & GPT_FIT_HIDDEN) {
	skipline = 1;
    }

    for (i=0; i<spec->n_lines; i++) {
	if (i == skipline || blank_user_line(&spec->lines[i])) {
	    continue;
	}
	if ((spec->flags & GPT_Y2AXIS) && spec->lines[i].yaxis != 1) {
	    any_y2 = 1;
	    break;
	}
    }

    if ((spec->bars != NULL && spec->nbars > 0) ||
	spec->data != NULL || spec->auxdata != NULL) {
	anydata = 1;
    }

    if (spec->heredata && anydata) {
	fputs("# start inline data\n", fp);
	if (spec->bars != NULL && spec->nbars > 0) {
	    print_plotbars(spec, fp);
	}
	if (spec->data != NULL || spec->auxdata != NULL) {
	    plotspec_print_heredata(spec, &miss, fp);
	}
	fputs("# end inline data\n", fp);
	strcpy(src, "$data");
    }

    fputs("plot \\\n", fp);

    gretl_push_c_numeric_locale();

    if (spec->nbars > 0) {
	print_bars_header(spec, fp);
    }

    ycol = 2;

    for (i=0; i<spec->n_lines; i++) {
	GPT_LINE *line = &spec->lines[i];

	if (i == skipline || blank_user_line(line)) {
	    if (i < spec->n_lines - 1) {
		continue;
	    } else {
		break;
	    }
	}

	if (!string_is_blank(line->formula)) {
	    fprintf(fp, "%s ", line->formula);
	} else if (line->ustr != NULL) {
	    fprintf(fp, "%s using %s ", src, line->ustr);
	} else if (line->mcols != NULL) {
	    int j;

	    fprintf(fp, "%s using ", src);
	    for (j=1; j<=line->mcols[0]; j++) {
		fprintf(fp, "%d", line->mcols[j]);
		fputc(j < line->mcols[0] ? ':' : ' ', fp);
	    }
	} else {
	    fprintf(fp, "%s using 1", src);
	    if (line->ncols == 5) {
		/* Note: boxplot candlesticks, hard-wired! */
		fprintf(fp, ":%d:%d:%d:%d", ycol+1, ycol, ycol+3, ycol+2);
	    } else if (line->ncols == 2) {
		if (line->flags & GP_LINE_BOXDATA) {
		    /* boxplot median, hard-wired */
		    fprintf(fp, ":%d:%d:%d:%d", ycol, ycol, ycol, ycol);
		} else {
		    fprintf(fp, ":%d", ycol);
		}
	    } else {
		for (k=2; k<=line->ncols; k++) {
		    fprintf(fp, ":%d", ycol + k - 2);
		}
	    }
	    fputc(' ', fp);
	}

	if ((spec->flags & GPT_Y2AXIS) && line->yaxis != 1) {
	    fprintf(fp, "axes x1y%d ", line->yaxis);
	}

	if (string_is_blank(line->title)) {
	    fputs("notitle ", fp);
	} else {
	    fprintf(fp, "title \"%s", line->title);
	    if (any_y2 && strchr(line->title, '(') == NULL) {
		if (line->yaxis == 1) {
		    fprintf(fp, " (%s)\" ", _("left"));
		} else {
		    fprintf(fp, " (%s)\" ", _("right"));
		}
	    } else {
		fputs("\" ", fp);
	    }
	}

	if (line->style == GP_STYLE_FILLEDCURVE && line->type == LT_AUTO) {
	    if (*line->rgb != '\0') {
		fprintf(fp, "lc rgb \"%s\" ", line->rgb);
	    } else {
		print_filledcurve_color(fp);
	    }
	}

	if (line->style != GP_STYLE_AUTO) {
	    fprintf(fp, "w %s", gp_line_style_name(line->style));
	}

	/* FIXME handling of linetype for PLOT_BAND? */
	if (spec->code != PLOT_BAND && line->type != LT_AUTO) {
	    fprintf(fp, " lt %d", line->type);
	} else if (spec->nbars > 0 && line->style != GP_STYLE_FILLEDCURVE) {
	    fprintf(fp, " lt %d", i + 1);
	}

	maybe_print_point_info(line, fp);

	if (line->style != GP_STYLE_FILLEDCURVE) {
	    if (line->width == 1.0 && spec->scale > 1.0) {
		fprintf(fp, " lw 2");
	    } else if (line->width != 1) {
		fprintf(fp, " lw %g", (double) line->width);
	    }
	}

	if (line->whiskwidth > 0) {
	    fprintf(fp, " whiskerbars %g", (double) line->whiskwidth);
	}

	if (spec->heredata) {
	    ycol += line->ncols - 1;
	}

	if (more_lines(spec, i, skipline)) {
	    fputs(", \\", fp);
	}

	fputc('\n', fp);
    }

    /* supply the data to gnuplot inline, old style */

    if (!spec->heredata && anydata) {
	if (spec->bars != NULL && spec->nbars > 0) {
	    print_plotbars(spec, fp);
	}
	if (spec->data != NULL || spec->auxdata != NULL) {
	    plotspec_print_data(spec, skipline, &miss, fp);
	}
    }

    gretl_pop_c_numeric_locale();

    if (png) {
	write_plot_bounding_box_request(fp);
    }

    return miss;
}

static int set_loess_fit (GPT_SPEC *spec, int d, double q, gretl_matrix *x,
			  gretl_matrix *y, gretl_matrix *yh)
{
    int i, t, T = gretl_vector_get_length(y);
    const double *src;
    gretl_matrix *m;

    m = gretl_matrix_alloc(T, 3);
    if (m == NULL) {
	return E_DATA;
    }

    gretl_matrix_free(spec->data);
    spec->data = m;

    for (i=0; i<3; i++) {
	src = i == 0 ? x->val : i == 1 ? y->val : yh->val;
	for (t=0; t<T; t++) {
	    gretl_matrix_set(m, t, i, src[t]);
	}
    }

    spec->data = m;
    spec->datacols = 3;
    spec->nobs = spec->okobs = T;

    g_free(spec->lines[1].title);
    spec->lines[1].title = g_strdup_printf(_("loess fit, d = %d, q = %g"), d, q);
    spec->lines[1].pscale = 1.0;
    spec->lines[1].style = GP_STYLE_LINES;
    spec->lines[1].ncols = 2;

    spec->fit = PLOT_FIT_LOESS;

    return 0;
}

static void set_plotfit_formula (GPT_LINE *line, FitType f, const double *b,
				 double x0, double pd)
{
    gchar *s = NULL;
    char x[32];

    g_free(line->formula);
    line->formula = NULL;
    gretl_push_c_numeric_locale();

    /* If x0 is not NA that indicates that we've calculated a fit
       against time (t = 1,2,3,...,T), and we have to transform the
       equation so that it produces correct fitted values given
       x-data in the form year[.fraction].
    */
    if (!na(x0)) {
	if (pd > 1) {
	    sprintf(x, "(%.1f*(x-%.8g)+1)", pd, x0);
	} else {
	    sprintf(x, "(x-%.8g+1)", x0);
	}
    } else {
	strcpy(x, "x");
    }

    if (f == PLOT_FIT_OLS) {
	s = g_strdup_printf("%.8g + %.8g*%s", b[0], b[1], x);
    } else if (f == PLOT_FIT_QUADRATIC) {
	s = g_strdup_printf("%.8g + %.8g*%s + %.8g*%s**2", b[0], b[1], x, b[2], x);
    } else if (f == PLOT_FIT_CUBIC) {
	s = g_strdup_printf("%.8g + %.8g*%s + %.8g*%s**2 + %.8g*%s**3",
			    b[0], b[1], x, b[2], x, b[3], x);
    } else if (f == PLOT_FIT_INVERSE) {
	s = g_strdup_printf("%.8g + %.8g/%s", b[0], b[1], x);
    } else if (f == PLOT_FIT_LOGLIN) {
	s = g_strdup_printf("exp(%.8g + %.8g*%s)", b[0], b[1], x);
    } else if (f == PLOT_FIT_LINLOG) {
	if (!na(x0)) {
	    s = g_strdup_printf("%.8g + %.8g*log%s", b[0], b[1], x);
	} else {
	    s = g_strdup_printf("%.8g + %.8g*log(x)", b[0], b[1]);
	}
    }

    gretl_pop_c_numeric_locale();
    line->formula = s;
}

void set_plotfit_line (GPT_LINE *line, FitType f, const double *b,
		       double x0, double pd)
{
    char xc = (na(x0)) ? 'X' : 't';
    gchar *s = NULL;

    g_free(line->title);
    line->title = NULL;

    /* first compose the key string for the fitted line */
    if (f == PLOT_FIT_OLS) {
	s = g_strdup_printf("Y = %#.3g %c %#.3g%c", b[0],
			    (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    } else if (f == PLOT_FIT_QUADRATIC) {
	s = g_strdup_printf("Y = %#.3g %c %#.3g%c %c %#.3g%c^2", b[0],
			    (b[1] > 0)? '+' : '-', fabs(b[1]), xc,
			    (b[2] > 0)? '+' : '-', fabs(b[2]), xc);
    } else if (f == PLOT_FIT_CUBIC) {
	s = g_strdup_printf("Y = %#.3g %c %#.3g%c %c %#.3g%c^2 %c %#.3g%c^3",
			    b[0], (b[1] > 0)? '+' : '-', fabs(b[1]), xc,
			    (b[2] > 0)? '+' : '-', fabs(b[2]), xc,
			    (b[3] > 0)? '+' : '-', fabs(b[3]), xc);
    } else if (f == PLOT_FIT_INVERSE) {
	s = g_strdup_printf("Y = %#.3g %c %#.3g(1/%c)", b[0],
			    (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    } else if (f == PLOT_FIT_LOGLIN) {
	s = g_strdup_printf("logY = %#.3g %c %#.3g%c", b[0],
			    (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
	if (xc == 't' && (pd == 1 || pd == 4 || pd == 12)) {
	    /* display annual growth rate in title */
	    double g = 100 * (pow(exp(b[1]), pd) - 1);
	    gchar *gstr, *tmp;

	    gstr = g_strdup_printf("\\n(%s %.2f%%)", _("annual growth"), g);
	    tmp = g_strdup_printf("%s%s", s, gstr);
	    g_free(s);
	    s = tmp;
	    g_free(gstr);
	}
    } else if (f == PLOT_FIT_LINLOG) {
	s = g_strdup_printf("Y = %#.3g %c %#.3glog(%c)", b[0],
			    (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    }

    line->title = s;

    /* then set the formula itself */
    set_plotfit_formula(line, f, b, x0, pd);
}

static void plotspec_set_fitted_line (GPT_SPEC *spec, FitType f,
				      double x0)
{
    double pd = spec->pd;
    const double *b;

    if (f == PLOT_FIT_OLS) {
	b = spec->b_ols->val;
    } else if (f == PLOT_FIT_QUADRATIC) {
	b = spec->b_quad->val;
    } else if (f == PLOT_FIT_CUBIC) {
	b = spec->b_cub->val;
    } else if (f == PLOT_FIT_INVERSE) {
	b = spec->b_inv->val;
    } else if (f == PLOT_FIT_LOGLIN) {
	b = spec->b_log->val;
    } else if (f == PLOT_FIT_LINLOG) {
	b = spec->b_linlog->val;
    } else {
	return;
    }

    set_plotfit_line(&spec->lines[1], f, b, x0, pd);

    spec->fit = f;
    spec->lines[1].pscale = 1.0;
    spec->lines[1].style = GP_STYLE_LINES;
    spec->lines[1].ncols = 0;
}

int plotspec_add_fit (GPT_SPEC *spec, FitType f)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *yh = NULL;
    const double *px = spec->data->val;
    const double *py;
    int T = spec->okobs;
    double x0 = NADBL;
    double xt, q = 0.5;
    int d = 1;
    int i, t, k;
    int err = 0;

    if ((spec->flags & GPT_TS) && f != PLOT_FIT_LOESS) {
	if (spec->pd == 1 || spec->pd == 4 || spec->pd == 12) {
	    /* Annual, quarterly or monthly data: we'll estimate
	       the fit against t = 1,2,3,...,T to give a more
	       comprehensible displayed equation. So we have
	       to record the initial x-value (e.g. 2008.25) to
	       get the fitted line right on the plot.
	    */
	    x0 = px[0];
	}
    }

    if ((f == PLOT_FIT_OLS && spec->b_ols != NULL) ||
	(f == PLOT_FIT_QUADRATIC && spec->b_quad != NULL) ||
	(f == PLOT_FIT_CUBIC && spec->b_cub != NULL) ||
	(f == PLOT_FIT_INVERSE && spec->b_inv != NULL) ||
	(f == PLOT_FIT_LOGLIN && spec->b_log != NULL) ||
	(f == PLOT_FIT_LINLOG && spec->b_linlog != NULL)) {
	/* just activate existing setup */
	plotspec_set_fitted_line(spec, f, x0);
	return 0;
    }

    if (f == PLOT_FIT_OLS || f == PLOT_FIT_INVERSE ||
	f == PLOT_FIT_LOGLIN || f == PLOT_FIT_LINLOG) {
	k = 2;
    } else if (f == PLOT_FIT_QUADRATIC) {
	k = 3;
    } else if (f == PLOT_FIT_CUBIC) {
	k = 4;
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

    px = spec->data->val;
    py = px + spec->nobs;

    i = 0;
    for (t=0; t<spec->nobs; t++) {
	if (!na(x0)) {
	    /* use a 1-based time index */
	    xt = t + 1;
	} else {
	    xt = px[t];
	}
	if (!na(py[t]) && !na(xt)) {
	    if (f == PLOT_FIT_LOGLIN) {
		/* y must be positive */
		if (py[t] > 0.0) {
		    y->val[i] = log(py[t]);
		} else {
		    gretl_errmsg_set(_("Non-positive values encountered"));
		    err = E_DATA;
		    goto bailout;
		}
	    } else {
		y->val[i] = py[t];
	    }
	    if (f == PLOT_FIT_LOESS) {
		gretl_matrix_set(X, i, 0, xt);
	    } else {
		gretl_matrix_set(X, i, 0, 1.0);
		if (f == PLOT_FIT_LINLOG) {
		    /* x must be positive */
		    if (xt > 0.0) {
			gretl_matrix_set(X, i, 1, log(xt));
		    } else {
			gretl_errmsg_set(_("Non-positive values encountered"));
			err = E_DATA;
			goto bailout;
		    }
		} else if (f == PLOT_FIT_INVERSE) {
		    gretl_matrix_set(X, i, 1, 1.0 / xt);
		} else {
		    gretl_matrix_set(X, i, 1, xt);
		}
		if (f == PLOT_FIT_QUADRATIC || f == PLOT_FIT_CUBIC) {
		    gretl_matrix_set(X, i, 2, xt * xt);
		}
		if (f == PLOT_FIT_CUBIC) {
		    gretl_matrix_set(X, i, 3, xt * xt * xt);
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
    } else if (f == PLOT_FIT_LOGLIN) {
	double s2;

	err = gretl_matrix_ols(y, X, b, NULL, NULL, &s2);
	if (!err) {
	    b->val[0] += s2 / 2;
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
	} else if (f == PLOT_FIT_CUBIC) {
	    spec->b_cub = b;
	    b = NULL;
	} else if (f == PLOT_FIT_INVERSE) {
	    spec->b_inv = b;
	    b = NULL;
	} else if (f == PLOT_FIT_LOGLIN) {
	    spec->b_log = b;
	    b = NULL;
	} else if (f == PLOT_FIT_LINLOG) {
	    spec->b_linlog = b;
	    b = NULL;
	}
    }

    if (!err) {
	if (f == PLOT_FIT_LOESS) {
	    set_loess_fit(spec, d, q, X, y, yh);
	    g_free(spec->lines[1].formula);
	    spec->lines[1].formula = NULL;
	} else {
	    plotspec_set_fitted_line(spec, f, x0);
	}
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(yh);

#if 0
    fprintf(stderr, "plotspec_add_fit: returning %d\n", err);
#endif

    return err;
}

/* next: apparatus for producing shaded bars in time-series plots */

#define BDEBUG 0

static plotbars *plotbars_new (int n)
{
    plotbars *bars = malloc(sizeof *bars);

    if (bars != NULL) {
	bars->dx = doubles_array_new(n, 2);
	if (bars->dx == NULL) {
	    free(bars);
	    bars = NULL;
	} else {
	    bars->n = n;
	    bars->t1 = bars->t2 = 0;
	    bars->ymin = bars->ymax = 0;
	}
    }

    return bars;
}

/* Read and check a plain text "plotbars" file.  Such a file may
   contain comment lines starting with '#'; other than comments, each
   line should contain a pair of space-separated date strings in gretl
   monthly date format (YYYY:MM).  These strings represent,
   respectively, the start and end of an episode of some kind
   (paradigm case: peak and trough of business cycle).  Each pair may
   be used to draw a vertical bar on a time-series plot. The file can
   contain any number of such pairs.  An example is provided in
   <prefix>/share/data/plotbars/nber.txt.

   The dates are converted internally into the format, year plus
   decimal fraction of year, and may be used when plotting
   annual, quarterly or monthly time series.
*/

static plotbars *parse_bars_file (const char *fname,
				  int *err)
{
    FILE *fp;
    plotbars *bars = NULL;
    char line[128];
    int ncolon = 0;
    int ndash = 0;
    int nyrs = 0;
    int nother = 0;
    int ntypes = 0;
    int d1, d2, d3, d4;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	*err = E_FOPEN;
	return NULL;
    }

#if BDEBUG > 1
    fprintf(stderr, "parse_bars_file: '%s'\n", fname);
#endif

    /* first count the data lines */

    while (fgets(line, sizeof line, fp)) {
	if (*line == '#' || string_is_blank(line)) {
	    continue;
	} else if (sscanf(line, "%d:%d %d:%d", &d1, &d2, &d3, &d4) == 4) {
	    ncolon++;
	} else if (sscanf(line, "%d-%d %d-%d", &d1, &d2, &d3, &d4) == 4) {
	    ndash++;
	} else if (sscanf(line, "%d %d", &d1, &d3) == 2) {
	    nyrs++;
	} else {
	    nother++;
	}
    }

    ntypes = (ncolon > 0) + (ndash > 0) + (nyrs > 0);

    /* initial check and allocation */

    if (nother > 0) {
	/* non-comment, non-blank, non-parseable cruft found */
	*err = E_DATA;
    } else if (ntypes == 0 || ntypes > 1) {
	/* no parseable data found, or a mash-up */
	*err = E_DATA;
    } else {
	bars = plotbars_new(ncolon + ndash + nyrs);
	if (bars == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err == 0) {
	double x0, x1;
	int i = 0;

	rewind(fp);

	/* now read, check, convert and record the date pairs */

	while (fgets(line, sizeof line, fp) && !*err) {
	    if (*line == '#' || string_is_blank(line)) {
		continue;
	    }
	    if (nyrs) {
		sscanf(line, "%d %d", &d1, &d3);
		d2 = 1;
		d4 = 12;
	    } else if (ncolon) {
		sscanf(line, "%d:%d %d:%d", &d1, &d2, &d3, &d4);
	    } else {
		sscanf(line, "%d-%d %d-%d", &d1, &d2, &d3, &d4);
	    }
	    x0 = d1 + (d2 - 1.0) / 12;
	    x1 = d3 + (nyrs ? d4 / 12.0 : (d4 - 1.0) / 12);
#if BDEBUG > 1
	    fprintf(stderr, "%.4f %.4f\n", x0, x1);
#endif
	    if (x1 < x0) {
		*err = E_DATA;
	    } else {
		bars->dx[i][0] = x0;
		bars->dx[i][1] = x1;
		i++;
	    }
	}
    }

    fclose(fp);

    if (*err == E_DATA) {
	gretl_errmsg_set(_("Dates file does not conform to the specification"));
    }

    if (*err && bars != NULL) {
	plotbars_free(bars);
	bars = NULL;
    }

    return bars;
}

/* output data representing vertical shaded bars: each
   bar is represented by a pair of rows */

static void print_plotbars (GPT_SPEC *spec, FILE *fp)
{
    plotbars *bars = spec->bars;
    double y0 = bars->ymin;
    double y1 = bars->ymax;
    int i, started, stopped;

    if (spec->heredata) {
	fputs("$bars << EOB\n", fp);
    }

    for (i=0; i<bars->n; i++) {
	if (bars->dx[i][1] < bars->t1) {
	    continue;
	}
	if (bars->dx[i][0] > bars->t2) {
	    break;
	}
	started = stopped = 0;
	if (bars->dx[i][0] >= bars->t1) {
	    /* start is in range */
	    fprintf(fp, "%.3f %.10g %.10g\n", bars->dx[i][0], y0, y1);
	    started = 1;
	}
	if (bars->dx[i][1] <= bars->t2) {
	    /* stop is in range */
	    if (!started) {
		/* but the start was not */
		fprintf(fp, "%.3f %.10g %.10g\n", bars->t1, y0, y1);
	    }
	    fprintf(fp, "%.3f %.10g %.10g\n", bars->dx[i][1], y0, y1);
	    if (!spec->heredata) {
		fputs("e\n", fp);
	    }
	    stopped = 1;
	}
	if (started && !stopped) {
	    /* truncated */
	    fprintf(fp, "%.3f %.10g %.10g\n", bars->t2, y0, y1);
	    if (!spec->heredata) {
		fputs("e\n", fp);
	    }
	}
    }

    if (spec->heredata) {
	fputs("EOB\n", fp);
    }
}

/* output special plot lines for creating vertical shaded bars
   in a time-series graph: the corresponding data will take
   the form x:ymin:ymax */

static void print_bars_header (GPT_SPEC *spec, FILE *fp)
{
    char cstr[8];
    int i;

    print_rgb_hash(cstr, get_shadecolor());

    if (spec->heredata) {
	int j = 0;

	for (i=0; i<spec->nbars; i++) {
	    fprintf(fp, "'$bars' every ::%d::%d using 1:2:3 notitle lc rgb \"%s\" "
		    "w filledcurve, \\\n", j, j+1, cstr);
	    j += 2;
	}
    } else {
	for (i=0; i<spec->nbars; i++) {
	    fprintf(fp, "'-' using 1:2:3 notitle lc rgb \"%s\" w filledcurve, \\\n",
		    cstr);
	}
    }
}

static void print_filledcurve_color (FILE *fp)
{
    char cstr[8];

    print_rgb_hash(cstr, get_shadecolor());
    fprintf(fp, "lc rgb \"%s\" ", cstr);
}

/* given the info in @bars, calculate how many of its start-stop
   pairs fall (at least partially) within the current x-axis
   range */

static int n_bars_shown (double xmin, double xmax, plotbars *bars)
{
    double **dx = bars->dx;
    int i, n = 0;

    for (i=0; i<bars->n; i++) {
	if (dx[i][1] < xmin) {
	    continue;
	} else if (dx[i][0] > xmax) {
	    break;
	} else if (dx[i][0] >= xmin || dx[i][1] <= xmax) {
	    n++;
	}
    }

    bars->t1 = xmin;
    bars->t2 = xmax;

#if BDEBUG
    fprintf(stderr, "n_bars_shown: n = %d\n", n);
#endif

    return n;
}

/* Given the limits of the data area of a plot, @xmin et al,
   and the name of a plain text file containing "bars"
   information in the form of a set of start-stop pairs,
   try attaching this information to @spec. We fail if
   there's something wrong with the info in the file, or
   if none of the bars defined therein would be visible
   in the current range, @xmin to @xmax.
*/

int plotspec_add_bars_info (GPT_SPEC *spec,
			    double xmin, double xmax,
			    double ymin, double ymax,
			    const char *fname)
{
    plotbars *bars = NULL;
    int err = 0;

    if (spec->bars != NULL) {
	/* clear out any stale info */
	plotbars_free(bars);
	spec->bars = NULL;
	spec->nbars = 0;
    }

    bars = parse_bars_file(fname, &err);

    if (!err) {
	int n = n_bars_shown(xmin, xmax, bars);

#if BDEBUG
	fprintf(stderr, "plotspec_add_bars_info:\n"
		" xmin=%g, xmax=%g, n=%d, ymin=%g, ymax=%g\n",
		xmin, xmax, n, ymin, ymax);
#endif
	if (n > 0) {
	    spec->bars = bars;
	    bars->ymin = ymin;
	    bars->ymax = ymax;
	    spec->nbars = n;
	} else {
	    plotbars_free(bars);
	}
    }

    return err;
}

/* the following is a public interface because it's also
   used when reconstituting @spec from a gnuplot command
   file.
*/

int plotspec_allocate_bars (GPT_SPEC *spec)
{
    int err = 0;

    spec->bars = plotbars_new(spec->nbars);
    if (spec->bars == NULL) {
	err = E_ALLOC;
    }

    return err;
}

/* for use when reconstituting @spec from a gnuplot
   command file: set the start/stop values for a
   particular bar @i.
*/

int plotspec_set_bar_info (GPT_SPEC *spec, int i,
			   double t1, double t2)
{
    if (i < spec->bars->n) {
	spec->bars->dx[i][0] = t1;
	spec->bars->dx[i][1] = t2;
	return 0;
    } else {
	return E_DATA;
    }
}

void plotspec_set_bars_limits (GPT_SPEC *spec,
			       double t1, double t2,
			       double ymin, double ymax)
{
    if (spec->bars != NULL) {
	spec->bars->t1 = t1;
	spec->bars->t2 = t2;
	spec->bars->ymin = ymin;
	spec->bars->ymax = ymax;
    }
}

void plotspec_remove_bars (GPT_SPEC *spec)
{
    plotbars_free(spec->bars);
    spec->bars = NULL;
    spec->nbars = 0;
}
