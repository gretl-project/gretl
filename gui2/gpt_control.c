/*
 *  Copyright (c) by Allin Cottrell
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

/* gpt_control.c for gretl -- gnuplot controller */

#include "gretl.h"
#include "gpt_control.h"
#include "session.h"
#include "gpt_dialog.h"

#undef POINTS_DEBUG

#ifdef G_OS_WIN32
# include <windows.h>
# include <io.h>
#endif

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk/gdkkeysyms.h>

#ifdef PNG_COMMENTS
# include <png.h>
#endif

#if !GLIB_CHECK_VERSION(2,0,0)
# define OLD_GTK
#endif

enum {
    PLOT_SAVED          = 1 << 0,
    PLOT_HAS_CONTROLLER = 1 << 1,
    PLOT_ZOOMED         = 1 << 2,
    PLOT_ZOOMING        = 1 << 3,
    PLOT_NO_LABELS      = 1 << 4,
    PLOT_PNG_COORDS     = 1 << 5,
    PLOT_DONT_ZOOM      = 1 << 6,
    PLOT_DONT_EDIT      = 1 << 7,
    PLOT_DONT_MOUSE     = 1 << 8,
    PLOT_POSITIONING    = 1 << 9
} plot_status_flags;

enum {
    PLOT_TITLE          = 1 << 0,
    PLOT_XLABEL         = 1 << 1,
    PLOT_YLABEL         = 1 << 2,
    PLOT_Y2AXIS         = 1 << 3,
    PLOT_Y2LABEL        = 1 << 4,
    PLOT_LABELS_UP      = 1 << 5,
} plot_format_flags;

#define plot_is_saved(p)        (p->status & PLOT_SAVED)
#define plot_has_controller(p)  (p->status & PLOT_HAS_CONTROLLER)
#define plot_is_zoomed(p)       (p->status & PLOT_ZOOMED)
#define plot_is_zooming(p)      (p->status & PLOT_ZOOMING)
#define plot_has_no_labels(p)   (p->status & PLOT_NO_LABELS)
#define plot_has_png_coords(p)  (p->status & PLOT_PNG_COORDS)
#define plot_not_zoomable(p)    (p->status & PLOT_DONT_ZOOM)
#define plot_not_editable(p)    (p->status & PLOT_DONT_EDIT)
#define plot_doing_position(p)  (p->status & PLOT_POSITIONING)

#define plot_has_title(p)       (p->format & PLOT_TITLE)
#define plot_has_xlabel(p)      (p->format & PLOT_XLABEL)
#define plot_has_ylabel(p)      (p->format & PLOT_YLABEL)
#define plot_has_y2axis(p)      (p->format & PLOT_Y2AXIS)
#define plot_has_y2label(p)     (p->format & PLOT_Y2LABEL)
#define plot_has_data_labels(p) (p->format & PLOT_LABELS_UP)

#define plot_is_range_mean(p)   (p->spec->code == PLOT_RANGE_MEAN)

enum {
    PNG_START,
    PNG_ZOOM,
    PNG_UNZOOM,
    PNG_REDISPLAY
} png_zoom_codes;

struct png_zoom_t {
    double xmin, xmax;
    double ymin, ymax;
    int screen_xmin, screen_ymin;
};

struct png_plot_t {
    GtkWidget *shell;
    GtkWidget *canvas;
    GtkWidget *popup;
    GtkWidget *statusarea;    
    GtkWidget *statusbar;
    GtkWidget *cursor_label;
    GtkWidget *labelpos_entry;
    GdkPixmap *pixmap;
    GdkGC *invert_gc;
    GPT_SPEC *spec;
    double xmin, xmax;
    double ymin, ymax;
    int pixel_xmin, pixel_xmax;
    int pixel_ymin, pixel_ymax;
    int xint, yint;
    int pd;
    guint cid;
    struct png_zoom_t *zoom;
    unsigned long status; 
    unsigned char format;
#ifndef OLD_GTK
    char *labeled;
#endif
};

static void render_pngfile (png_plot *plot, int view);
static int zoom_unzoom_png (png_plot *plot, int view);
static void create_selection_gc (png_plot *plot);
static int get_plot_ranges (png_plot *plot);

#ifdef G_OS_WIN32
enum {
    WIN32_TO_CLIPBOARD,
    WIN32_TO_PRINTER
};
#endif /* G_OS_WIN32 */

#ifdef PNG_COMMENTS

enum {
    GRETL_PNG_OK,
    GRETL_PNG_NO_OPEN,
    GRETL_PNG_NOT_PNG,
    GRETL_PNG_NO_COMMENTS,
    GRETL_PNG_BAD_COMMENTS,
    GRETL_PNG_NO_COORDS
};

typedef struct png_bounds_t png_bounds;

struct png_bounds_t {
    int xleft;
    int xright;
    int ybot;
    int ytop;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

static int get_png_bounds_info (png_bounds *bounds);

#endif /* PNG_COMMENTS */

#define PLOTSPEC_DETAILS_IN_MEMORY(s)  (s->data != NULL)

static void terminate_plot_positioning (png_plot *plot)
{
    plot->status ^= PLOT_POSITIONING;
    plot->labelpos_entry = NULL;
    gdk_window_set_cursor(plot->canvas->window, NULL);
    gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
    raise_gpt_control_window();
}

void plot_remove_controller (png_plot *plot) 
{
    if (plot_has_controller(plot)) {
	plot->status ^= PLOT_HAS_CONTROLLER;
	if (plot_doing_position(plot)) {
	    terminate_plot_positioning(plot);
	}
    }
}

int plot_is_mouseable (const png_plot *plot)
{
    return !(plot->status & PLOT_DONT_MOUSE);
}

void set_plot_has_y2_axis (png_plot *plot, gboolean s)
{
    if (s == TRUE) {
	plot->format |= PLOT_Y2AXIS;
    } else {
	plot->format &= ~PLOT_Y2AXIS;
    }
}

void plot_label_position_click (GtkWidget *w, GPT_SPEC *spec)
{
    png_plot *plot;

    plot = (png_plot *) spec->ptr;
    if (plot != NULL) {
	GtkWidget *entry;
	GdkCursor* cursor;

	cursor = gdk_cursor_new(GDK_CROSSHAIR);
	gdk_window_set_cursor(plot->canvas->window, cursor);
	gdk_cursor_destroy(cursor);
	entry = g_object_get_data(G_OBJECT(w), "labelpos_entry");
	plot->labelpos_entry = entry;
	plot->status |= PLOT_POSITIONING;
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid, 
			   _(" Click to set label position"));
    }
}

static void line_to_file (const char *s, FILE *fp, int l2)
{
#ifdef ENABLE_NLS
    if (l2 == -1) {
	print_as_html(s, fp);
    } else if (l2 == 1) {
	print_as_locale(s, fp);
    } else {
	fputs(s, fp);
    }
#else
    fputs(s, fp);
#endif
}

static FILE *open_gp_file (const char *fname, const char *mode)
{
    FILE *fp = fopen(fname, mode);

    if (fp == NULL) {
	if (*mode == 'w') {
	    sprintf(errtext, _("Couldn't write to %s"), fname);
	} else {
	    sprintf(errtext, _("Couldn't open %s"), fname);
	}
        errbox(errtext);
    }

    return fp;
}

static int commented_term_line (const char *s)
{
    return !strncmp(s, "# set term png", 14);
}

static int set_output_line (const char *s)
{
    return !strncmp(s, "set output", 10);
}

static int add_or_remove_png_term (const char *fname, int add, GPT_SPEC *spec)
{
    FILE *fsrc, *ftmp;
    char temp[MAXLEN], fline[MAXLEN];
    char restore_line[MAXLEN] = {0};
    int png_line_saved = 0;
#ifdef ENABLE_NLS
    int l2 = use_latin_2();
#else
    int l2 = 0;
#endif

    sprintf(temp, "%sgpttmp.XXXXXX", paths.userdir);
    if (mktemp(temp) == NULL) return 1;

    ftmp = open_gp_file(temp, "w");
    if (ftmp == NULL) return 1;

    fsrc = open_gp_file(fname, "r");
    if (fsrc == NULL) {
	fclose(ftmp);
	return 1;
    }

    if (add && spec == NULL) {
	/* see if there's a commented out term setting to restore */
	while (fgets(fline, sizeof fline, fsrc)) {
	    if (commented_term_line(fline)) {
		strcat(restore_line, fline + 2);
		break;
	    }
	}
	rewind(fsrc);
    }

    if (add) {
	int need_term_line = 1;

	if (spec == NULL) {
	    /* we're reconstituting a session graph from a 
	       saved gnuplot command file */
	    if (*restore_line) {
		/* found a saved png term specification */
		fputs(restore_line, ftmp);
		need_term_line = 0;
	    }
	}
	if (need_term_line) {
	    fprintf(ftmp, "%s\n",
		    get_gretl_png_term_line(PLOT_REGULAR));
	}	    
	fprintf(ftmp, "set output '%sgretltmp.png'\n", 
		paths.userdir);
    }

    /* now for the body of the plot file */
    while (fgets(fline, sizeof fline, fsrc)) {
	if (add) {
	    if (!commented_term_line(fline) && !set_output_line(fline)) {
		line_to_file(fline, ftmp, -l2);
	    }
	} else {
	    /* we're removing the png term line */
	    int printit = 1;

	    if (!strncmp(fline, "set term png", 12)) {
		if (!png_line_saved) {
		    /* comment it out, for future reference */
		    fprintf(ftmp, "# %s", fline);
		    png_line_saved = 1;
		} 
		printit = 0;
	    }
	    else if (commented_term_line(fline)) {
		printit = 0;
	    } else if (set_output_line(fline)) {
		printit = 0;
	    } else if (spec != NULL && (spec->flags & GPTSPEC_OLS_HIDDEN)
		       && is_auto_ols_string(fline)) {
		printit = 0;
	    }
	    if (printit) {
		line_to_file(fline, ftmp, l2);
	    }
	}
    }

    fclose(fsrc);
    fclose(ftmp);

    remove(fname);
    return rename(temp, fname);
}

static int add_png_term_to_plotfile (const char *fname)
{
    return add_or_remove_png_term(fname, 1, NULL);
}

int remove_png_term_from_plotfile (const char *fname, GPT_SPEC *spec)
{
    /* called from session.c when saving a graph file */
    return add_or_remove_png_term(fname, 0, spec);
}

void mark_plot_as_saved (GPT_SPEC *spec)
{
    png_plot *plot = (png_plot *) spec->ptr;

    plot->status |= PLOT_SAVED;
}

static int gnuplot_png_init (GPT_SPEC *spec, FILE **fpp)
{
    *fpp = fopen(spec->fname, "w");
    if (*fpp == NULL) {
	sprintf(errtext, _("Couldn't write to %s"), spec->fname);
	errbox(errtext);
	return 1;
    }
    fprintf(*fpp, "%s\n", get_gretl_png_term_line(spec->code));
    fprintf(*fpp, "set output '%sgretltmp.png'\n", paths.userdir);
    return 0;
}

void display_session_graph_png (char *fname) 
{
    char *myfname = fname;
    gchar *plotcmd;
    int err = 0;

    /* take saved plot source file and make PNG from it, then display
       the PNG */
    if (add_png_term_to_plotfile(myfname)) {
	return;
    }

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, myfname);
#ifdef G_OS_WIN32
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    err = gretl_spawn(plotcmd);
#endif
    g_free(plotcmd);

    if (err) {
	errbox(_("Gnuplot error creating graph"));
    } else {
	gnuplot_show_png(myfname, NULL, 1);
    }
}

#ifdef G_OS_WIN32
static void win32_process_graph (GPT_SPEC *spec, int color, int dest);
#endif

void save_this_graph (GPT_SPEC *plot, const char *fname)
{
    FILE *fq;
    PRN *prn;
    char plottmp[MAXLEN], plotline[MAXLEN], termstr[MAXLEN];
    gchar *plotcmd = NULL;
    int cmds, err;

    if (!user_fopen("gptout.tmp", plottmp, &prn)) {
	return;
    }

    fq = fopen(plot->fname, "r");
    if (fq == NULL) {
	errbox(_("Couldn't access graph info"));
	gretl_print_destroy(prn);
	return;
    }
 
    cmds = get_termstr(plot, termstr);
  
    if (cmds) {
	if (copyfile(plot->fname, fname)) { 
	    errbox(_("Failed to copy graph file"));
	}
	return;
    } else {
#ifdef ENABLE_NLS
	pprint_gnuplot_encoding(termstr, prn);
#endif /* ENABLE_NLS */
	pprintf(prn, "set term %s\n", termstr);
	pprintf(prn, "set output '%s'\n", fname);
	while (fgets(plotline, MAXLEN-1, fq)) {
	    if (strncmp(plotline, "set term", 8) && 
		strncmp(plotline, "set output", 10))
		pputs(prn, plotline);
	}
    }

    gretl_print_destroy(prn);
    fclose(fq);

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
			      plottmp);

#ifdef G_OS_WIN32
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    err = gretl_spawn(plotcmd);
#endif /* G_OS_WIN32 */

    remove(plottmp);
    g_free(plotcmd);

    if (err) {
	errbox(_("Gnuplot error creating graph"));
    } else {
	infobox(_("Graph saved"));
    }
}

/* chop trailing comma, if present; return 1 if comma chopped,
   zero otherwise */

static int chop_comma (char *str)
{
    size_t i, n = strlen(str);

    for (i=n-1; i>0; i--) {
	if (isspace((unsigned char) str[i])) {
	    continue;
	}
	if (str[i] == ',') {
	    str[i] = 0;
	    return 1;
	} else {
	    break;
	}
    }
		
    return 0;
}

/* ........................................................... */

static int get_gpt_label (const char *line, char *label)
{
    const char *p = strchr(line, '#');
    char format[6];

    if (p != NULL) {
	sprintf(format, "%%%ds", OBSLEN - 1);
	sscanf(p + 1, format, label);
	return 0;
    }

    return 1;
}

/* ........................................................... */

static void get_gpt_data (char *line, double *x, double *y)
{
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    if (x != NULL) {
	if (sscanf(line, "%lf %lf", x, y) == 2) {
	    ; /* we're fine */
	} else {
	    if (sscanf(line, "%lf", x) == 1) {
		*y = NADBL;
	    } else if (sscanf(line, "%*s %lf", y) == 1) {
		*x = NADBL;
	    } else {
		*x = *y = NADBL;
	    }
	}
    } else {
	if (sscanf(line, "%*s %lf", y) != 1) {
	    *y = NADBL;
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif
}

/* ........................................................... */

static int cant_zoom (const char *line)
{
    if (strncmp(line, "# multi", 7) == 0) {
	/* multiple scatterplots */
	return 1;
    }
    if (strncmp(line, "# CUSUM", 7) == 0) {
	/* graph from CUSUM test */
	return 1;
    }
    if (strncmp(line, "# frequ", 7) == 0) {
	/* frequency distribution plot */
	return 1;
    }
    if (strncmp(line, "# sampl", 7) == 0) {
	/* sampling distribution graph (stats calculator) */
	return 1;
    }
    if (strncmp(line, "# corre", 7) == 0) {
	/* correlogram */
	return 1;
    }
    if (strncmp(line, "# perio", 7) == 0) {
	/* periodogram */
	return 1;
    }
    if (strncmp(line, "# range", 7) == 0) {
	/* range-mean plot */
	return 1;
    }
    if (strncmp(line, "# lever", 7) == 0) {
	/* leverage plot */
	return 1;
    }
    if (strstr(line, "no auto-parse")) {
	/* general prohibition on gui plot editing */
	return 1;
    }
    return 0;
}

/* ........................................................... */

static int cant_edit (const char *line)
{
    if (strncmp(line, "# mult", 6) == 0) {
	/* multiple scatterplots */
	return 1;
    }
    if (strncmp(line, "# CUSUM", 7) == 0) {
	/* graph from CUSUM test */
	return 1;
    }
    if (strncmp(line, "# corre", 7) == 0) {
	/* correlogram */
	return 1;
    }
    if (strncmp(line, "# sampl", 7) == 0) {
	/* sampling distribution graph (stats calculator) */
	return 1;
    }
    if (strncmp(line, "# lever", 7) == 0) {
	/* leverage plot */
	return 1;
    }
    if (strstr(line, "no auto-parse")) {
	/* general prohibition on gui plot editing */
	return 1;
    }
    return 0;
}

/* ........................................................... */

static GPT_SPEC *plotspec_new (void)
{
    GPT_SPEC *spec;
    int i;

    spec = mymalloc(sizeof *spec);
    if (spec == NULL) {
	return NULL;
    }

    spec->lines = mymalloc(MAX_PLOT_LINES * sizeof *spec->lines);

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
    }

    for (i=0; i<4; i++) {
	spec->titles[i][0] = 0;
	spec->literal[i] = NULL;
    }

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	strcpy(spec->text_labels[i].text, "");
	strcpy(spec->text_labels[i].just, "left");
	strcpy(spec->text_labels[i].pos, "0,0");
    }

    spec->xtics[0] = 0;
    spec->mxtics[0] = 0;
    spec->fname[0] = 0;
    strcpy(spec->keyspec, "left top");

    for (i=0; i<3; i++) {
	strcpy(spec->range[i][0], "*");
	strcpy(spec->range[i][1], "*");
    }

    spec->code = PLOT_REGULAR;
    spec->flags = 0;
    spec->fp = NULL;
    spec->data = NULL;
    spec->labels = NULL;
    spec->nlabels = 0;
    spec->ptr = NULL;
    spec->nlines = 0;
    spec->n_y_series = 0;

    spec->termtype[0] = 0;
    spec->t1 = spec->t2 = 0;

    return spec;
}

/* ........................................................... */

static int parse_label_line (GPT_SPEC *spec, const char *line, int i)
{
    const char *p, *s;
    int n, x, y;
    char coord[8];

    /* e.g. set label 'foobar' at [ screen | graph ] 1500,350 left */

    if (i >= MAX_PLOT_LABELS) return 1;

    strcpy(spec->text_labels[i].text, "");
    strcpy(spec->text_labels[i].pos, "");
    strcpy(spec->text_labels[i].just, "");

    p = strstr(line, "'");
    if (p == NULL) return 1;
    p++;
    s = p;

    /* get the label text */
    while (*s) {
	if (*s == '\'') {
	    int len = s - p;

	    if (len > PLOT_LABEL_TEXT_LEN) {
		len = PLOT_LABEL_TEXT_LEN;
	    }
	    strncat(spec->text_labels[i].text, p, len);
	    break;
	}
	s++;
    }

    /* get the position */
    p = strstr(s, "at");
    if (p == NULL) {
	strcpy(spec->text_labels[i].text, "");
	return 1;
    }
    p += 2;

    /* coordinate system? */
    *coord = 0;
    s = strstr(p, "graph");
    if (s != NULL) {
	strcpy(coord, "graph ");
	p = s + 5;
    } else {
	s = strstr(p, "screen");
	if (s != NULL) {
	    strcpy(coord, "screen ");
	    p = s + 6;
	}
    }

    /* actual coordinates */
    n = sscanf(p, "%d,%d", &x, &y);
    if (n != 2) {
	strcpy(spec->text_labels[i].text, "");
	return 1;
    }
    sprintf(spec->text_labels[i].pos, "%s%d,%d", coord, x, y);

    /* justification */
    if (strstr(p, "left")) {
	strcpy(spec->text_labels[i].just, "left");
    } else if (strstr(p, "right")) {
	strcpy(spec->text_labels[i].just, "right");
    } else if (strstr(p, "center")) {
	strcpy(spec->text_labels[i].just, "center");
    } else {
	strcpy(spec->text_labels[i].just, "left");
    }	

    return 0;
}

/* ........................................................... */

static int parse_gp_set_line (GPT_SPEC *spec, const char *line,
			      int *i, int *labelno)
{
    char set_thing[12], setting[MAXLEN], range[32];
    size_t n, j;

    if (strstr(line, "encoding")) return 0;

    if (sscanf(line + 4, "%11s", set_thing) != 1) {
	errbox(_("Failed to parse gnuplot file"));
	fprintf(stderr, "plotfile line: '%s'\n", line);
	return 1;
    }

    if (strcmp(set_thing, "y2tics") == 0) {
	spec->flags |= GPTSPEC_Y2AXIS;
	return 0;
    }
    else if (strcmp(set_thing, "border 3") == 0) {
	spec->flags |= GPTSPEC_BORDER_HIDDEN;
	return 0;
    }    

    n = strlen(set_thing);
    strcpy(setting, line + 4 + n);
    top_n_tail(setting);

    if (strstr(set_thing, "range")) {
	strcpy(range, setting + 1);
	n = strlen(range);
	for (j=0; j<n; j++) {
	    if (range[j] == ':') {
		range[j] = 0;
		strcpy(spec->range[*i][0], range);
		break;
	    }
	}
	strcpy(range, strchr(setting, ':') + 1);
	delchar(']', range);
	strcpy(spec->range[*i][1], range);
	*i += 1;
    }
    else if (strcmp(set_thing, "title") == 0) 
	strcpy(spec->titles[0], setting);
    else if (strcmp(set_thing, "xlabel") == 0)
	strcpy(spec->titles[1], setting);
    else if (strcmp(set_thing, "ylabel") == 0)
	strcpy(spec->titles[2], setting);
    else if (strcmp(set_thing, "y2label") == 0)
	strcpy(spec->titles[3], setting);
    else if (strcmp(set_thing, "key") == 0)
	strcpy(spec->keyspec, setting);
    else if (strcmp(set_thing, "nokey") == 0)
	strcpy(spec->keyspec, "none");
    else if (strcmp(set_thing, "xtics") == 0) 
	safecpy(spec->xtics, setting, 15);
    else if (strcmp(set_thing, "mxtics") == 0) 
	safecpy(spec->mxtics, setting, 3);
    else if (strcmp(set_thing, "label") == 0) {
	parse_label_line(spec, line, *labelno);
	*labelno += 1;
    }

    return 0;
}

/* ........................................................... */

static int allocate_plotspec_labels (GPT_SPEC *spec, int plot_n)
{
    int i;

    spec->labels = mymalloc(plot_n * sizeof *spec->labels);
    if (spec->labels == NULL) {
	return 1;
    }

    for (i=0; i<plot_n; i++) {
	spec->labels[i] = malloc(OBSLEN);
	if (spec->labels[i] == NULL) {
	    free(spec->labels);
	    spec->nlabels = 0;
	    return 1;
	}
	spec->labels[i][0] = 0;
    }

    spec->nlabels = plot_n;

    return 0;
}

/* ........................................................... */

static int get_plot_n (FILE *fp, int *got_labels)
{
    int n = 0, started = -1;
    char line[MAXLEN], label[OBSLEN];
    char *p;

    *got_labels = 0;

    while (fgets(line, MAXLEN - 1, fp)) {
	if (!strncmp(line, "plot", 4)) {
	    started = 0;
	}
	if (started == 0 && strchr(line, '\\') == NULL) {
	    started = 1;
	    continue;
	}
	if (started == 1) {
	    if (*got_labels == 0 && (p = strchr(line, '#')) != NULL) {
		if (sscanf(p + 1, "%8s", label) == 1) {
		    *got_labels = 1;
		}
	    }
	    if (*line == 'e') break;
	    n++;
	}
    }

    return n;
}

/* Read in plotspec struct from gnuplot command file.  This is _not_ a
   general parser for gnuplot files; it is designed specifically for
   files auto-generated by gretl.
*/

static int read_plotspec_from_file (GPT_SPEC *spec)
{
    int i, j, t, n, plot_n, done, labelno;
    int got_labels = 0;
    char line[MAXLEN], *got = NULL, *p = NULL;
    double *tmpy = NULL;
    size_t diff;
    FILE *fp;

    /* check: are we already done? */
    if (PLOTSPEC_DETAILS_IN_MEMORY(spec)) {
	return 0;
    }

    /* open the plot file */
    fp = fopen(spec->fname, "r");
    if (fp == NULL) {
	errbox(_("Couldn't open graph file"));
	return 1;
    }

    /* get the number of data-points, and check for labels */
    plot_n = get_plot_n(fp, &got_labels);
    if (plot_n == 0) {
	fclose(fp);
	return 1;
    }
    rewind(fp);

    /* get the preamble and "set" lines */
    i = 0;
    labelno = 0;
    while ((got = fgets(line, MAXLEN - 1, fp))) {
	if (cant_edit(line)) {
	    goto plot_bailout;
	}
	if (strncmp(line, "# timeseries", 12) == 0) {
	    spec->flags |= GPTSPEC_TS;
	    continue;
	}
	if (strncmp(line, "# freq", 6) == 0 ||
	    strncmp(line, "# peri", 6) == 0) {
	    /* special cases */
	    if (line[2] == 'f') {
		if (strstr(line, "normal")) {
		    spec->code = PLOT_FREQ_NORMAL;
		} else if (strstr(line, "gamma")) {
		    spec->code = PLOT_FREQ_GAMMA;
		} else {
		    spec->code = PLOT_FREQ_SIMPLE;
		}
	    } else {
		spec->code = PLOT_PERIODOGRAM;
	    }
	    if (spec->code != PLOT_FREQ_SIMPLE) {
		/* grab special plot lines */
		for (j=0; j<4; j++) {
		    spec->literal[j] = mymalloc(MAXLEN);
		    if (spec->literal[j] == NULL) {
			return 1;
		    }
		    if (!fgets(spec->literal[j], MAXLEN - 1, fp)) {
			errbox(_("Plot file is corrupted"));
			free(spec->literal[j]);
			spec->literal[j] = NULL;
			goto plot_bailout;
		    }
		    top_n_tail(spec->literal[j]);
		}
	    }
	    continue;
	}
	if (strncmp(line, "# forecast", 10) == 0) {
	    spec->code = PLOT_FORECAST;
	    continue;
	}
	if (strstr(line, "automatic OLS")) {
	    spec->flags |= GPTSPEC_AUTO_OLS;
	    continue;
	}

	/* ignore an unknown comment line */
	if (strncmp(line, "# ", 2) == 0) {
	    continue;
	}

	if (strncmp(line, "set ", 4)) {
	    break;
	}
	if (parse_gp_set_line(spec, line, &i, &labelno)) {
	    goto plot_bailout;
	}
    }

    if (got == NULL) {
	goto plot_bailout;
    }

    for (i=0; i<4; i++) {
	if (spec->titles[i][0] != 0) {
	    delchar('\'', spec->titles[i]);
	}
    }

    if (spec->keyspec[0] == 0) {
	strcpy(spec->keyspec, "none");
    }

    /* then get the "plot" lines */
    if (strncmp(line, "plot ", 4) ||
	(strlen(line) < 10 && fgets(line, MAXLEN - 1, fp) == NULL)) {	
	errbox(_("Failed to parse gnuplot file"));
	fprintf(stderr, "plotfile line: '%s'\n", line);
	goto plot_bailout;
    }

    i = 0;
    done = 0;
    while (1) {
	top_n_tail(line);

	if (!chop_comma(line)) {
	    /* line did not end with comma -> no contination of
	       the plot command */
	    done = 1;
	} 

	/* scale, [yaxis,] style */
	p = strstr(line, "using");
	if (p && p[11] == '*') {
            safecpy(spec->lines[i].scale, p + 12, 7);
	    charsub(spec->lines[i].scale, ')', '\0');
	} else {
	    if (p) { 
		strcpy(spec->lines[i].scale, "1.0");
	    } else {
		strcpy(spec->lines[i].scale, "NA");
		p = strstr(line, "axes");
		if (p == NULL) {
		    p = strstr(line, "title");
		}
		if (p != NULL) {
		    diff = p - line;
		    strncpy(spec->lines[i].formula, line, diff);
		    spec->lines[i].formula[diff - 1] = 0;
		} 
	    }
	}

	spec->lines[i].yaxis = 1;
	if (strstr(line, "axes x1y2")) {
	    spec->lines[i].yaxis = 2;
	}

	p = strstr(line, "title");
	if (p != NULL) {
	    sscanf(p + 7, "%79[^']'", spec->lines[i].title);
	}

	p = strstr(line, " w ");
	if (p != NULL) {
	    sscanf(p + 3, "%15[^, ]", spec->lines[i].style);
	} else {
	    strcpy(spec->lines[i].style, "points");
	}

	if (done) {
	    break;
	}

	i++;

	if ((got = fgets(line, MAXLEN - 1, fp)) == NULL) {
	    break;
	}
    }

    if (got == NULL) {
	goto plot_bailout;
    }

    spec->nlines = i + 1; /* i is a zero-based index */

    /* free any unused lines */
    if (spec->nlines < MAX_PLOT_LINES) {
	spec->lines = myrealloc(spec->lines, 
				spec->nlines * sizeof *spec->lines);
    }

    /* finally, get the plot data, and labels if any */
    spec->data = mymalloc(plot_n * (i + 2) * sizeof *spec->data);
    tmpy = mymalloc(plot_n * sizeof *tmpy);
    if (spec->data == NULL || tmpy == NULL) {
	goto plot_bailout;
    }

    if (got_labels) {
	if (allocate_plotspec_labels(spec, plot_n)) {
	    free(spec->data);
	    spec->data = NULL;
	    free(tmpy);
	    goto plot_bailout;
	}
    }

    /* Below: read the data from the plot.  There may be more
       than one y series. */
    j = 1;
    t = 0;
    n = 0;
    while (fgets(line, MAXLEN - 1, fp)) {
	if (line[0] == 'e') {
	    n = t;
	    t = 0;
	    j++;
	    continue;
	}
	if (strncmp(line, "pause", 5) == 0) {
	    break;
	}
	if (j == 1) { 
	    /* first set: read both x and y (and label?) */ 
	    get_gpt_data(line, &(spec->data[t]), &(tmpy[t]));
	    if (got_labels) {
		get_gpt_label(line, spec->labels[t]);
	    }
	} else {      
	    /* any subsequent sets: read y only */ 
	    get_gpt_data(line, NULL, &(spec->data[j*n + t]));
	}
	t++;
    }

    spec->n_y_series = j - 1;

    /* put "tmpy" in as last data column */
    for (t=0; t<n; t++) {
	spec->data[n + t] = tmpy[t];
    }
    free(tmpy);

    spec->t1 = 0;
    spec->t2 = n - 1;

    /* see if we really have any labels */
    if (spec->labels != NULL && !got_labels) {
	for (i=0; i<plot_n; i++) {
	    free(spec->labels[i]);
	}
	free(spec->labels);
	spec->labels = NULL;
	spec->nlabels = 0;
    }

#if 0
    fprintf(stderr, "spec->labels = %p, spec->nlabels = %d\n",
	    (void *) spec->labels, spec->nlabels);
#endif

    fclose(fp);
    return 0;

 plot_bailout:

    fclose(fp);
    return 1;
}

/* Size of drawing area */
#define PLOT_PIXEL_WIDTH  640   /* try 576? 608? */
#define PLOT_PIXEL_HEIGHT 480   /* try 432? 456? */

#ifdef USE_GNOME
extern void gnome_print_graph (const char *fname);
#endif

static int get_data_xy (png_plot *plot, int x, int y, 
			double *data_x, double *data_y)
{
    double xmin, xmax;
    double ymin, ymax;

    if (plot_is_zoomed(plot)) {
	xmin = plot->zoom->xmin;
	xmax = plot->zoom->xmax;
	ymin = plot->zoom->ymin;
	ymax = plot->zoom->ymax;
    } else {
	xmin = plot->xmin;
	xmax = plot->xmax;
	ymin = plot->ymin;
	ymax = plot->ymax;
    }

#ifdef POINTS_DEBUG
    if (plot_doing_position(plot)) {
	fprintf(stderr, "get_data_xy:\n"
		" plot->xmin=%g, plot->xmax=%g, plot->ymin=%g, plot->ymax=%g\n",
		plot->xmin, plot->xmax, plot->ymin, plot->ymax);
    }
#endif

    if (xmin == 0.0 && xmax == 0.0) { /* unknown x range */
	fprintf(stderr, "get_data_xy: unknown x range\n");
	*data_x = NADBL;
    } else {
	*data_x = xmin + ((double) x - plot->pixel_xmin) / 
	    (plot->pixel_xmax - plot->pixel_xmin) * (xmax - xmin);
    }

    if (!na(*data_x)) {
	if (ymin == 0.0 && ymax == 0.0) { /* unknown y range */
	    fprintf(stderr, "get_data_xy: unknown y range\n");
	    *data_y = NADBL;
	} else {
	    *data_y = ymax - ((double) y - plot->pixel_ymin) / 
		(plot->pixel_ymax - plot->pixel_ymin) * (ymax - ymin);
	}
    }

    return (!na(*data_x) && !na(*data_y));
}

static void x_to_date (double x, int pd, char *str)
{
    int yr = (int) x;
    double t, frac = 1.0 / pd;
    int subper = (int) ((x - yr + frac) * pd);
    static int decpoint;

    if (decpoint == 0) decpoint = get_local_decpoint();

    t = yr + subper / ((pd < 10)? 10.0 : 100.0);
    sprintf(str, "%.*f", (pd < 10)? 1 : 2, t);
    charsub(str, decpoint, ':');
}

static void create_selection_gc (png_plot *plot)
{
    if (plot->invert_gc == NULL) {
	plot->invert_gc = gdk_gc_new(plot->canvas->window);
	gdk_gc_set_function(plot->invert_gc, GDK_INVERT);
    }
}

static void draw_selection_rectangle (png_plot *plot,
				      int x, int y)
{
    int rx, ry, rw, rh;

    rx = (plot->zoom->screen_xmin < x)? plot->zoom->screen_xmin : x;
    ry = (plot->zoom->screen_ymin < y)? plot->zoom->screen_ymin : y;
    rw = x - plot->zoom->screen_xmin;
    rh = y - plot->zoom->screen_ymin;
    if (rw < 0) rw = -rw;
    if (rh < 0) rh = -rh;    

    /* draw one time to make the rectangle appear */
    gdk_draw_rectangle(plot->pixmap,
		       plot->invert_gc,
		       FALSE,
		       rx, ry, rw, rh);
    /* show the modified pixmap */
    gdk_window_copy_area(plot->canvas->window,
			 plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
			 0, 0,
			 plot->pixmap,
			 0, 0,
			 PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT);
    /* draw (invert) again to erase the rectangle */
    gdk_draw_rectangle(plot->pixmap,
		       plot->invert_gc,
		       FALSE,
		       rx, ry, rw, rh);
}

#ifdef OLD_GTK

static void
write_label_to_plot (png_plot *plot, const gchar *label,
		     gint x, gint y)
{
    static GdkFont *label_font;

    if (plot->invert_gc == NULL) {
	create_selection_gc(plot);
    }

    if (label_font == NULL) {
	label_font = gdk_font_load("fixed");
    }

    /* draw the label */
    gdk_draw_text (plot->pixmap,
		   label_font,
		   plot->invert_gc,
		   x, y,
		   label,
		   strlen(label));

    /* show the modified pixmap */
    gdk_window_copy_area(plot->canvas->window,
			 plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
			 0, 0,
			 plot->pixmap,
			 0, 0,
			 PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT);

    /* draw (invert) again to erase the text */
    gdk_draw_text (plot->pixmap,
		   label_font,
		   plot->invert_gc,
		   x, y,
		   label,
		   strlen(label));
}

#else

static void
write_label_to_plot (png_plot *plot, const gchar *label,
		     gint x, gint y)
{
    PangoContext *context;
    PangoLayout *pl;

    if (plot->invert_gc == NULL) {
	create_selection_gc(plot);
    }

    context = gtk_widget_get_pango_context(plot->shell);
    pl = pango_layout_new(context);
    pango_layout_set_text(pl, label, -1);

    /* draw the label */
    gdk_draw_layout (plot->pixmap, plot->invert_gc, x, y, pl);

    /* show the modified pixmap */
    gdk_window_copy_area(plot->canvas->window,
			 plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
			 0, 0,
			 plot->pixmap,
			 0, 0,
			 PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT);

    /* trash the pango layout */
    g_object_unref (G_OBJECT(pl));

    /* record that a label is shown */
    plot->format |= PLOT_LABELS_UP;
}

#endif /* GTK versions */

#define TOLDIST 0.01

static void
identify_point (png_plot *plot, int pixel_x, int pixel_y,
		double x, double y) 
{
    double xrange, yrange;
    double xdiff, ydiff;
    double min_xdist, min_ydist;
    int best_match = -1;
    int t, plot_n;
    const double *data_x, *data_y = NULL;

    if (!PLOTSPEC_DETAILS_IN_MEMORY(plot->spec)) {
	if (read_plotspec_from_file(plot->spec)) return;
    }

    /* no labels to show */
    if (plot->spec->labels == NULL) {
	plot->status |= PLOT_NO_LABELS;	
	return;
    }

    plot_n = plot->spec->t2 - plot->spec->t1 + 1;

#ifndef OLD_GTK
    /* need array to keep track of which points are labeled */
    if (plot->labeled == NULL) {
	plot->labeled = mymalloc(plot_n);
	if (plot->labeled == NULL) return;
	memset(plot->labeled, 0, plot_n);
    }
#endif

    if (plot_is_zoomed(plot)) {
	min_xdist = xrange = plot->zoom->xmax - plot->zoom->xmin;
	min_ydist = yrange = plot->zoom->ymax - plot->zoom->ymin;
    } else {
	min_xdist = xrange = plot->xmax - plot->xmin;
	min_ydist = yrange = plot->ymax - plot->ymin;
    }

    data_x = &plot->spec->data[0];
    /* there's an ambiguity here: in case of more than one y series,
       what do we want to use as "data_y"? */
    if (plot_has_y2axis(plot)) {
	/* use first y-var that's on y1 axis */
	int i;

	for (i=0; i<plot->spec->nlines; i++) {
	    if (plot->spec->lines[i].yaxis == 1) {
		data_y = &plot->spec->data[(i + 1) * plot_n];
		break;
	    }
	}
    } 

    if (data_y == NULL) {
	data_y = &plot->spec->data[plot->spec->n_y_series * plot_n];
    }

    /* try to find the best-matching data point */
    for (t=0; t<plot_n; t++) {
	if (na(data_x[t]) || na(data_y[t])) continue;
	xdiff = fabs(data_x[t] - x);
	ydiff = fabs(data_y[t] - y);
	if (xdiff <= min_xdist && ydiff <= min_ydist) {
	    min_xdist = xdiff;
	    min_ydist = ydiff;
	    best_match = t;
	}
    }

#ifndef OLD_GTK
    /* if the point is already labeled, skip */
    if (plot->labeled[best_match]) return;
#endif

    /* if the match is good enough, show the label */
    if (best_match >= 0 && min_xdist < TOLDIST * xrange &&
	min_ydist < TOLDIST * yrange) {
	write_label_to_plot(plot, plot->spec->labels[best_match],
			    pixel_x, pixel_y);
#ifndef OLD_GTK
	/* flag the point as labeled already */
	plot->labeled[best_match] = 1;
#endif
    }
}

#define MAX_LABELS_SHOWN 250

static gint
motion_notify_event (GtkWidget *widget, GdkEventMotion *event,
		     png_plot *plot)
{
    int x, y;
    GdkModifierType state;
    gchar label[32], label_y[16];

    if (event->is_hint) {
        gdk_window_get_pointer(event->window, &x, &y, &state);
    } else {
        x = event->x;
        y = event->y;
        state = event->state;
    }

    *label = 0;

    if (x > plot->pixel_xmin && x < plot->pixel_xmax && 
	y > plot->pixel_ymin && y < plot->pixel_ymax) {
	double data_x, data_y;

	get_data_xy(plot, x, y, &data_x, &data_y);
	if (na(data_x)) return TRUE;

	if (!(plot_has_no_labels(plot)) && 
	    datainfo->t2 - datainfo->t1 < MAX_LABELS_SHOWN &&
	    !plot_is_zooming(plot) &&
	    !na(data_y)) {
	    identify_point(plot, x, y, data_x, data_y);
	}

	if (plot->pd == 4 || plot->pd == 12) {
	    x_to_date(data_x, plot->pd, label);
	} else {
	    sprintf(label, (plot->xint)? "%7.0f" : "%7.4g", data_x);
	}

	if (!na(data_y)) {
	    if (plot_has_png_coords(plot)) {
		sprintf(label_y, (plot->yint)? " %-7.0f" : " %-7.4g", data_y);
	    } else {
		/* pretty much guessing at y coordinate here */
		sprintf(label_y, (plot->yint)? " %-7.0f" : " %-6.3g", data_y);
	    }
	    strcat(label, label_y);
	}

	if (plot_is_zooming(plot) && (state & GDK_BUTTON1_MASK)) {
	    draw_selection_rectangle(plot, x, y);
	}
    } 

    gtk_label_set_text(GTK_LABEL(plot->cursor_label), label);
  
    return TRUE;
}

static void set_plot_format_flags (png_plot *plot)
{
    plot->format = 0;

    if (!string_is_blank(plot->spec->titles[0])) {
	plot->format |= PLOT_TITLE;
    }
    if (!string_is_blank(plot->spec->titles[1])) {
	plot->format |= PLOT_XLABEL;
    }
    if (!string_is_blank(plot->spec->titles[2])) {
	plot->format |= PLOT_YLABEL;
    }
    if (!string_is_blank(plot->spec->titles[3])) {
	plot->format |= PLOT_Y2LABEL;
    }
    if (plot->spec->flags & GPTSPEC_Y2AXIS) {
	plot->format |= PLOT_Y2AXIS;
    }
}

static void start_editing_png_plot (png_plot *plot)
     /* called from png plot popup menu */
{
    /* the spec struct is not yet filled out by reference
       to the gnuplot source file 
    */
    if (read_plotspec_from_file(plot->spec)) {
	return;
    }

    if (show_gnuplot_dialog(plot->spec) == 0) { /* OK */
	plot->status |= PLOT_HAS_CONTROLLER;
    }
}

#ifdef HAVE_AUDIO
static void audio_render_plot (png_plot *plot)
{
    void *handle;
    int (*midi_play_graph) (const char *, const char *, const char *);

    if (plot_not_editable(plot)) {
	return;
    }

    midi_play_graph = gui_get_plugin_function("midi_play_graph", 
					      &handle);
    if (midi_play_graph == NULL) {
        return;
    }

# ifdef G_OS_WIN32
    (*midi_play_graph) (plot->spec->fname, paths.userdir, NULL);
# else
    (*midi_play_graph) (plot->spec->fname, paths.userdir, midiplayer);
# endif

    close_plugin(handle);
}
#endif

static gint color_popup_activated (GtkWidget *w, gpointer data)
{
    gchar *item = (gchar *) data;
    gpointer ptr = g_object_get_data(G_OBJECT(w), "plot");
    png_plot *plot = (png_plot *) ptr;
    gint color = strcmp(item, _("monochrome"));
    GtkWidget *parent = (GTK_MENU(w->parent))->parent_menu_item;
    gchar *parent_item = g_object_get_data(G_OBJECT(parent), "string");

    if (!strcmp(parent_item, _("Save as postscript (EPS)..."))) {
	strcpy(plot->spec->termtype, "postscript");
	if (color) {
	    strcat(plot->spec->termtype, " color");
	}
	file_selector(_("Save gnuplot graph"), SAVE_THIS_GRAPH, 
		      plot->spec);
    } 
    else if (!strcmp(parent_item, _("Save as Windows metafile (EMF)..."))) {
	strcpy(plot->spec->termtype, "emf");
	if (color) {
	    strcat(plot->spec->termtype, " color");
	}
	file_selector(_("Save gnuplot graph"), SAVE_THIS_GRAPH, 
		      plot->spec);
    } 
#ifdef G_OS_WIN32
    else if (!strcmp(parent_item, _("Copy to clipboard"))) {
	win32_process_graph(plot->spec, color, WIN32_TO_CLIPBOARD);
    } 
    else if (!strcmp(parent_item, _("Print"))) {
	win32_process_graph(plot->spec, color, WIN32_TO_PRINTER);
    }    
#endif   

    return TRUE;
}

static gint plot_popup_activated (GtkWidget *w, gpointer data)
{
    gchar *item = (gchar *) data;
    gpointer ptr = g_object_get_data(G_OBJECT(w), "plot");
    png_plot *plot = (png_plot *) ptr;
    int killplot = 0;

    gtk_widget_destroy(plot->popup);
    plot->popup = NULL;

    if (!strcmp(item, _("Save as PNG..."))) {
	strcpy(plot->spec->termtype, "png");
        file_selector(_("Save gnuplot graph"), SAVE_THIS_GRAPH, plot->spec);
    }
    else if (!strcmp(item, _("Save to session as icon"))) { 
	add_graph_to_session(plot->spec, GRETL_GNUPLOT_GRAPH, NULL);
    }
    else if (plot_is_range_mean(plot) && !strcmp(item, _("Help"))) { 
	context_help(NULL, GINT_TO_POINTER(RMPLOT));
    }
#ifndef OLD_GTK
    else if (!strcmp(item, _("Clear data labels"))) { 
	zoom_unzoom_png(plot, PNG_START);
    }
#endif
    else if (!strcmp(item, _("Zoom..."))) { 
	GdkCursor* cursor;

	cursor = gdk_cursor_new(GDK_CROSSHAIR);
	gdk_window_set_cursor(plot->canvas->window, cursor);
	gdk_cursor_destroy(cursor);
	plot->status |= PLOT_ZOOMING;
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid, 
			   _(" Drag to define zoom rectangle"));
	create_selection_gc(plot);
    }
    else if (!strcmp(item, _("Restore full view"))) { 
	zoom_unzoom_png(plot, PNG_UNZOOM);
    }
#ifdef USE_GNOME 
    else if (!strcmp(item, _("Print..."))) { 
	gnome_print_graph(plot->spec->fname);
    }
#endif 
    else if (!strcmp(item, _("Edit"))) { 
	start_editing_png_plot(plot);
    }
    else if (!strcmp(item, _("Close"))) { 
        killplot = 1;
    } 

    if (killplot) gtk_widget_destroy(plot->shell);

    return TRUE;
}

static void attach_color_popup (GtkWidget *w, png_plot *plot)
{
    GtkWidget *item, *cpopup;
    const char *color_items[] = {
	N_("color"),
	N_("monochrome")
    };
    int i;

    cpopup = gtk_menu_new();

    for (i=0; i<2; i++) {
	item = gtk_menu_item_new_with_label(_(color_items[i]));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(color_popup_activated),
			 _(color_items[i]));
	g_object_set_data(G_OBJECT(item), "plot", plot);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(cpopup), item);
    } 

    gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), cpopup);
}

static void build_plot_menu (png_plot *plot)
{
    GtkWidget *item;    
    const char *regular_items[] = {
#ifdef G_OS_WIN32
	N_("Save as Windows metafile (EMF)..."),
#endif
	N_("Save as PNG..."),
        N_("Save as postscript (EPS)..."),
#ifndef G_OS_WIN32
	N_("Save as Windows metafile (EMF)..."),
#endif
#ifdef G_OS_WIN32
	N_("Copy to clipboard"),
#endif
	N_("Save to session as icon"),
#ifndef OLD_GTK
	N_("Clear data labels"),
#endif
	N_("Zoom..."),
#ifdef USE_GNOME
	N_("Print..."),
#endif
#ifdef G_OS_WIN32
	N_("Print"),
#endif
	N_("Edit"),
	N_("Help"),
        N_("Close"),
        NULL
    };
    const char *zoomed_items[] = {
	N_("Restore full view"),
	N_("Close"),
	NULL
    };
    const char **plot_items;
    int i;

    plot->popup = gtk_menu_new();

    if (plot_is_zoomed(plot)) {
	plot_items = zoomed_items;
    } else {
	plot_items = regular_items;
    }

    i = 0;
    while (plot_items[i]) {
	if (plot_not_zoomable(plot) &&
	    !strcmp(plot_items[i], "Zoom...")) {
	    i++;
	    continue;
	}
	if (!plot_is_range_mean(plot) &&
	    !strcmp(plot_items[i], "Help")) {
	    i++;
	    continue;
	}
	if (plot_is_saved(plot) &&
	    !strcmp(plot_items[i], "Save to session as icon")) {
	    i++;
	    continue;
	}
	if ((plot_has_controller(plot) || plot_not_editable(plot)) &&
	    !strcmp(plot_items[i], "Edit")) {
	    i++;
	    continue;
	}
#ifndef OLD_GTK
	if (!plot_has_data_labels(plot) &&
	    !strcmp(plot_items[i], "Clear data labels")) {
	    i++;
	    continue;
	}
#endif

        item = gtk_menu_item_new_with_label(_(plot_items[i]));
        g_object_set_data(G_OBJECT(item), "plot", plot);
        gtk_widget_show(item);
        gtk_menu_shell_append(GTK_MENU_SHELL(plot->popup), item);

	/* items with color sub-menu */
	if (!strcmp(plot_items[i], "Save as Windows metafile (EMF)...") ||
	    !strcmp(plot_items[i], "Save as postscript (EPS)...") ||
	    !strcmp(plot_items[i], "Copy to clipboard") ||
	    !strcmp(plot_items[i], "Print")) {
	    attach_color_popup(item, plot);
	    g_object_set_data(G_OBJECT(item), "string", _(plot_items[i]));
	} else {
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(plot_popup_activated),
			     _(plot_items[i]));
	}
        i++;
    }

    g_signal_connect(G_OBJECT(plot->popup), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), 
		     &plot->popup);
}

int redisplay_edited_png (png_plot *plot)
{
    gchar *plotcmd;
    FILE *fp;
    int err = 0;

    /* open file in which to dump plot specification */
    gnuplot_png_init(plot->spec, &fp);
    if (fp == NULL) return 1;

    /* dump the edited plot details to file */
    set_png_output(plot->spec);
    print_plotspec_details(plot->spec, fp);
    fclose(fp);

    /* get gnuplot to create a new PNG graph */
    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
			      plot->spec->fname);
#ifdef G_OS_WIN32
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else    
    err = gretl_spawn(plotcmd);
#endif
    g_free(plotcmd);

    if (err) {
	errbox(_("Failed to generate PNG file"));
	return 1;
    }

    /* grab (possibly modified) data ranges */
    get_plot_ranges(plot);

    /* put the newly created PNG onto the plot canvas */
    render_pngfile(plot, PNG_REDISPLAY);

    return 0;
}

static int zoom_unzoom_png (png_plot *plot, int view)
{
    int err = 0;
    char fullname[MAXLEN];
    gchar *plotcmd = NULL;

    if (view == PNG_ZOOM) {
	FILE *fpin, *fpout;
	char line[MAXLEN];

	fpin = fopen(plot->spec->fname, "r");
	if (fpin == NULL) return 1;

	build_path(paths.userdir, "zoomplot.gp", fullname, NULL);
	fpout = fopen(fullname, "w");
	if (fpout == NULL) {
	    fclose(fpin);
	    return 1;
	}

	/* switch to zoomed data range */
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
	fprintf(fpout, "set xrange [%g:%g]\n", plot->zoom->xmin,
		plot->zoom->xmax);
	fprintf(fpout, "set yrange [%g:%g]\n", plot->zoom->ymin,
		plot->zoom->ymax);
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif

	while (fgets(line, MAXLEN-1, fpin)) {
	    if (strncmp(line, "set xrange", 10) &&
		strncmp(line, "set yrange", 10))
		fputs(line, fpout);
	}

	fclose(fpout);
	fclose(fpin);

	plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
				  fullname);
    } else { /* PNG_UNZOOM or PNG_START */
	plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
				  plot->spec->fname);
    }

#ifdef G_OS_WIN32
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    err = gretl_spawn(plotcmd);
#endif
    g_free(plotcmd);  

    if (view == PNG_ZOOM) {
	remove(fullname);
    }

    if (err) {
	errbox(_("Failed to generate PNG file"));
	return 1;
    }

    render_pngfile(plot, view);

    return 0;
}

static gint plot_button_release (GtkWidget *widget, GdkEventButton *event, 
				 png_plot *plot)
{
    if (plot_is_zooming(plot)) {
	double z;

	if (!get_data_xy(plot, event->x, event->y, 
			 &plot->zoom->xmax, &plot->zoom->ymax)) {
	    return TRUE;
	}

	/* flip the selected rectangle if required */
	if (plot->zoom->xmin > plot->zoom->xmax) {
	    z = plot->zoom->xmax;
	    plot->zoom->xmax = plot->zoom->xmin;
	    plot->zoom->xmin = z;
	}

	if (plot->zoom->ymin > plot->zoom->ymax) {
	    z = plot->zoom->ymax;
	    plot->zoom->ymax = plot->zoom->ymin;
	    plot->zoom->ymin = z;
	}

	if (plot->zoom->xmin != plot->zoom->xmax &&
	    plot->zoom->ymin != plot->zoom->ymax) {
	    zoom_unzoom_png(plot, PNG_ZOOM);
	}

	plot->status ^= PLOT_ZOOMING;
	gdk_window_set_cursor(plot->canvas->window, NULL);
	gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
    }

    return TRUE;
}

static gint plot_button_press (GtkWidget *widget, GdkEventButton *event, 
			       png_plot *plot)
{
    if (plot_is_zooming(plot)) {
	/* think about this */
	if (get_data_xy(plot, event->x, event->y, 
			&plot->zoom->xmin, &plot->zoom->ymin)) {
	    plot->zoom->screen_xmin = event->x;
	    plot->zoom->screen_ymin = event->y;
	}
	return TRUE;
    }

    if (plot_doing_position(plot)) {
	if (plot->labelpos_entry != NULL) {
	    double dx, dy;
	    
	    if (get_data_xy(plot, event->x, event->y, &dx, &dy)) {
		gchar *posstr;

#ifdef ENABLE_NLS
		setlocale(LC_NUMERIC, "C");
#endif
		posstr = g_strdup_printf("%g,%g", dx, dy);
#ifdef ENABLE_NLS
		setlocale(LC_NUMERIC, "");
#endif
		gtk_entry_set_text(GTK_ENTRY(plot->labelpos_entry), posstr);
		g_free(posstr);
	    }
	} 
	terminate_plot_positioning(plot);
	return TRUE;
    }

    if (plot->popup != NULL) {
	gtk_widget_destroy(plot->popup);
	plot->popup = NULL;
    }

    build_plot_menu(plot);
    gtk_menu_popup(GTK_MENU(plot->popup), NULL, NULL, NULL, NULL,
                   event->button, event->time);

    return TRUE;
}

static gboolean 
plot_key_handler (GtkWidget *w, GdkEventKey *key, png_plot *plot)
{
    switch (key->keyval) {
    case GDK_q:
    case GDK_Q:
	gtk_widget_destroy(w);
	break;
    case GDK_s:
    case GDK_S:
	add_graph_to_session(plot->spec, GRETL_GNUPLOT_GRAPH, NULL);
	break;
#ifdef G_OS_WIN32
    case GDK_c:
	win32_process_graph(plot->spec, 1, WIN32_TO_CLIPBOARD);
	break;
#endif
#ifdef HAVE_AUDIO
    case GDK_a:
    case GDK_A:
	audio_render_plot(plot);
	break;
#endif
    default:
	break;
    }

    return TRUE;
}

static 
void plot_expose (GtkWidget *widget, GdkEventExpose *event,
		  GdkPixmap *dbuf_pixmap)
{
    /* Don't repaint entire window on each exposure */
    gdk_window_set_back_pixmap(widget->window, NULL, FALSE);

    /* Refresh double buffer, then copy the "dirtied" area to
       the on-screen GdkWindow */
    gdk_window_copy_area(widget->window,
			 widget->style->fg_gc[GTK_STATE_NORMAL],
			 event->area.x, event->area.y,
			 dbuf_pixmap,
			 event->area.x, event->area.y,
			 event->area.width, event->area.height);
}

static void render_pngfile (png_plot *plot, int view)
{
    gint width;
    gint height;
    GdkPixbuf *pbuf;
    char pngname[MAXLEN];
#ifndef OLD_GTK
    GError *error = NULL;
#endif

    build_path(paths.userdir, "gretltmp.png", pngname, NULL);

#ifdef OLD_GTK
    pbuf = gdk_pixbuf_new_from_file(pngname);
    if (pbuf == NULL) {
	errbox(_("Failed to create pixbuf from file"));
	remove(pngname);
	return;
    }
#else
    pbuf = gdk_pixbuf_new_from_file(pngname, &error);
    if (pbuf == NULL) {
        errbox(error->message);
        g_error_free(error);
	remove(pngname);
	return;
    }
#endif

    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);

    if (width == 0 || height == 0) {
	errbox(_("Malformed PNG file for graph"));
#ifdef OLD_GTK
	gdk_pixbuf_unref(pbuf);
#else
	g_object_unref(pbuf);
#endif
	remove(pngname);
	return;
    }

#ifndef OLD_GTK
    /* scrap any old record of which points are labeled */
    if (plot->labeled != NULL) {
	free(plot->labeled);
	plot->labeled = NULL;
	plot->format ^= PLOT_LABELS_UP;
    }
#endif

    gdk_pixbuf_render_to_drawable(pbuf, plot->pixmap, 
				  plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
				  0, 0, 0, 0, width, height,
				  GDK_RGB_DITHER_NONE, 0, 0);

#ifdef OLD_GTK
    gdk_pixbuf_unref(pbuf);
#else
    g_object_unref(pbuf);
#endif
    remove(pngname);
    
    if (view != PNG_START) { 
	/* we're changing the view, so refresh the whole canvas */
	gdk_window_copy_area(plot->canvas->window,
			     plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
			     0, 0,
			     plot->pixmap,
			     0, 0,
			     PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT);
	if (view == PNG_ZOOM) {
	    plot->status |= PLOT_ZOOMED;
	} else if (view == PNG_UNZOOM) {
	    plot->status ^= PLOT_ZOOMED;
	}
    }
}

static void destroy_png_plot (GtkWidget *w, png_plot *plot)
{
    /* delete temporary plot source file? */
    if (!plot_is_saved(plot)) {
	remove(plot->spec->fname);
    }

    if (plot_has_controller(plot)) {
	/* if the png plot has a controller, destroy it too */
	plot->spec->ptr = NULL;
	destroy_gpt_control_window();
    } else {
	/* no controller: take responsibility for freeing the
	   plot specification */
	free_plotspec(plot->spec);
    }

    /* free allocated elements of png_plot struct */
    if (plot->zoom != NULL) {
	free(plot->zoom);
    }
#ifndef OLD_GTK
    if (plot->labeled != NULL) {
	free(plot->labeled);
    }
#endif
    if (plot->invert_gc != NULL) {
	gdk_gc_destroy(plot->invert_gc);
    }

    free(plot);
}

static void set_approx_pixel_bounds (png_plot *plot, 
				     int max_num_width,
				     int max_num2_width)
{
    if (PLOTSPEC_DETAILS_IN_MEMORY(plot->spec)) {
	set_plot_format_flags(plot);
    }

    if (plot_has_xlabel(plot)) {
	plot->pixel_ymax = PLOT_PIXEL_HEIGHT - 36;
    } else {
	plot->pixel_ymax = PLOT_PIXEL_HEIGHT - 24;
    }

    if (plot_has_title(plot)) {
	plot->pixel_ymin = 36;
    } else {
	plot->pixel_ymin = 14;
    }

    plot->pixel_xmin = 27 + 7 * max_num_width;
    if (plot_has_ylabel(plot)) {
	plot->pixel_xmin += 12;
    }

    plot->pixel_xmax = PLOT_PIXEL_WIDTH - 20; 
    if (plot_has_y2axis(plot)) {
	plot->pixel_xmax -= 7 * (max_num2_width + 1);
    }
    if (plot_has_y2label(plot)) {
	plot->pixel_xmax -= 11;
    }

#ifdef POINTS_DEBUG
    fprintf(stderr, "set_approx_pixel_bounds():\n"
	    " xmin=%d xmax=%d ymin=%d ymax=%d\n", 
	    plot->pixel_xmin, plot->pixel_xmax,
	    plot->pixel_ymin, plot->pixel_ymax);
    fprintf(stderr, "set_approx_pixel_bounds():\n"
	    " max_num_width=%d max_num2_width=%d\n", 
	    max_num_width, max_num2_width);
#endif
}

/* Attempt to read y-range info from the ascii representation
   of a gnuplot graph (the "dumb" terminal): return 0 on
   success, non-zero on failure.
*/

static int get_dumb_plot_yrange (png_plot *plot)
{
    FILE *fpin, *fpout;
    char line[MAXLEN], dumbgp[MAXLEN], dumbtxt[MAXLEN];
    gchar *plotcmd = NULL;
    int err = 0, x2axis = 0;
    int max_ywidth = 0;
    int max_y2width = 0;

    fpin = fopen(plot->spec->fname, "r");
    if (fpin == NULL) {
	return 1;
    }

    build_path(paths.userdir, "dumbplot.gp", dumbgp, NULL);
    build_path(paths.userdir, "gptdumb.txt", dumbtxt, NULL);
    fpout = fopen(dumbgp, "w");
    if (fpout == NULL) {
	fclose(fpin);
	return 1;
    }

    /* switch to the "dumb" (ascii) terminal in gnuplot */
    while (fgets(line, MAXLEN-1, fpin)) {
	if (strstr(line, "set term")) {
	    fputs("set term dumb\n", fpout);
	} else if (strstr(line, "set output")) { 
	    fprintf(fpout, "set output '%s'\n", dumbtxt);
	} else {
	    fputs(line, fpout);
	}
	if (strstr(line, "x2range")) {
	    x2axis = 1;
	}
    }

    fclose(fpin);
    fclose(fpout);

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot,
			      dumbgp);

#ifdef G_OS_WIN32
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    err = gretl_spawn(plotcmd);
#endif
    
    g_free(plotcmd);
    remove(dumbgp);

    if (err) {
#ifdef POINTS_DEBUG
	fputs("get_dumb_plot_yrange(): plot command failed\n", stderr);
#endif
	return 1;
    } else {
	double y[16] = {0};
	int y_numwidth[16] = {0};
	int y2_numwidth[16] = {0};
	char numstr[32];
	int i, j, k, imin;

	fpin = fopen(dumbtxt, "r");
	if (fpin == NULL) {
	    return 1;
	}

	/* read the y-axis min and max from the ascii graph */
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "C");
#endif
	i = j = 0;
	while (i < 16 && fgets(line, MAXLEN-1, fpin)) {
	    const char *s = line;
	    int nsp = 0;

	    while (isspace((unsigned char) *s)) {
	        nsp++;
	        s++;
            }
	    if (nsp > 5) {
		/* not a y-axis number */
		continue; 
	    }
	    if (sscanf(s, "%lf", &y[i]) == 1) {
#ifdef POINTS_DEBUG
		fprintf(stderr, "from text plot: read y[%d]=%g\n",
			i, y[i]);
#endif
		sscanf(s, "%31s", numstr);
		y_numwidth[i++] = strlen(numstr);
	    }
	    if (plot_has_y2axis(plot) && j < 16) {
		double y2;

		s = strrchr(s, ' ');
		if (s != NULL && sscanf(s, "%lf", &y2) == 1) {
		    sscanf(s, "%31s", numstr);
		    y2_numwidth[j++] = strlen(numstr);
		}
	    }
	}
#ifdef ENABLE_NLS
	setlocale(LC_NUMERIC, "");
#endif

	fclose(fpin);
#ifndef POINTS_DEBUG
	remove(dumbtxt);
#endif

	imin = (x2axis)? 1 : 0;

	if (i > (imin + 2) && y[imin] > y[i-1]) {
	    plot->ymin = y[i-1];
	    plot->ymax = y[imin];
	    for (k=imin; k<i-1; k++) {
		if (y_numwidth[k] > max_ywidth) {
		    max_ywidth = y_numwidth[k];
		}
	    }
	}	    

#ifdef POINTS_DEBUG
	fprintf(stderr, "Reading y range from text plot: plot->ymin=%g, "
		"plot->ymax=%g\n", plot->ymin, plot->ymax);
#endif

	if (plot_has_y2axis(plot)) {
	    for (k=imin; k<j-2; k++) {
		if (y2_numwidth[k] > max_y2width) {
		    max_y2width = y2_numwidth[k];
		}
	    }
	}
    }

    if (plot->ymax <= plot->ymin) {
	err = 1;
    }

    if (!err) {
	set_approx_pixel_bounds(plot, max_ywidth, max_y2width);
    }
    
    return err;
}

/* note: return 0 on failure */

static int get_plot_ranges (png_plot *plot)
{
    FILE *fp;
    char line[MAXLEN];
    int noedit = 0;
    int got_x = 0, got_y = 0, got_pd = 0;
#ifdef PNG_COMMENTS
    png_bounds b;
#endif

    plot->xmin = plot->xmax = 0.0;
    plot->ymin = plot->ymax = 0.0;   
    plot->xint = plot->yint = 0;
    plot->pd = 0;

    fp = fopen(plot->spec->fname, "r");
    if (fp == NULL) return 0;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    while (fgets(line, MAXLEN-1, fp)) {
	if (cant_edit(line)) {
	    plot->status |= PLOT_DONT_EDIT;
	    plot->status |= PLOT_DONT_ZOOM;
	    noedit = 1;
	    break;
	}
	if (strstr(line, "# forecasts with 95")) {
	    /* auto-parse can't handle the error bars */
	    plot->status |= PLOT_DONT_EDIT;
	}	
	else if (strstr(line, "# range-mean")) {
	    plot->spec->code = PLOT_RANGE_MEAN;
	}
	else if (sscanf(line, "set xrange [%lf:%lf]", 
			&plot->xmin, &plot->xmax) == 2) { 
	    got_x = 1;
	} 
	else if (sscanf(line, "# timeseries %d", &plot->pd) == 1) {
	    got_pd = 1;
	}
	if (!PLOTSPEC_DETAILS_IN_MEMORY(plot->spec)) {
	    if (!strncmp(line, "set tit", 7)) {
		plot->format |= PLOT_TITLE;
	    } else if (!strncmp(line, "set xla", 7)) {
		plot->format |= PLOT_XLABEL;
	    } else if (!strncmp(line, "set yla", 7)) {
		plot->format |= PLOT_YLABEL;
	    } else if (!strncmp(line, "set y2la", 8)) {
		plot->format |= PLOT_Y2LABEL;
	    } else if (!strncmp(line, "set y2ti", 8)) {
		plot->format |= PLOT_Y2AXIS;
	    } else if (cant_zoom(line)) {
		plot->status |= PLOT_DONT_ZOOM;
	    } 
	}
	if (!strncmp(line, "plot ", 5)) break;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    if (noedit) {
	return 0;
    }

#ifdef PNG_COMMENTS
    /* now try getting accurate coordinate info from PNG file */
    if (get_png_bounds_info(&b) == GRETL_PNG_OK) {
	plot->status |= PLOT_PNG_COORDS;
	got_x = got_y = 1;
	plot->pixel_xmin = b.xleft;
	plot->pixel_xmax = b.xright;
	plot->pixel_ymin = PLOT_PIXEL_HEIGHT - b.ytop;
	plot->pixel_ymax = PLOT_PIXEL_HEIGHT - b.ybot;
	plot->xmin = b.xmin;
	plot->xmax = b.xmax;
	plot->ymin = b.ymin;
	plot->ymax = b.ymax;
# ifdef POINTS_DEBUG
	fprintf(stderr, "get_png_bounds_info():\n"
		" xmin=%d xmax=%d ymin=%d ymax=%d\n", 
		plot->pixel_xmin, plot->pixel_xmax,
		plot->pixel_ymin, plot->pixel_ymax);
# endif
    } else {
	fprintf(stderr, "get_png_bounds_info(): failed\n");
    }
#endif /* PNG_COMMENTS */

    /* If got_x = 0 at this point, we didn't an x-range out of
       the gnuplot source file OR the PNG file, so we might as
       well give up.
    */

    if (got_x) {
	int err = 0;

	/* get the "dumb" coordinates only if we haven't got
	   more accurate ones from the PNG file */
	if (!plot_has_png_coords(plot)) { 
	    err = get_dumb_plot_yrange(plot);
	}

	if (!err) {
	    got_y = 1;
	    if ((plot->xmax - plot->xmin) / 
		(plot->pixel_xmax - plot->pixel_xmin) >= 1.0) {
		plot->xint = 1;
	    }
	    if ((plot->ymax - plot->ymin) / 
		(plot->pixel_ymax - plot->pixel_ymin) >= 1.0) {
		plot->yint = 1;
	    }
	}
    }

    if (!got_x || !got_y) {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
#ifdef POINTS_DEBUG 
	fputs("get_plot_ranges: setting PLOT_DONT_ZOOM, PLOT_DONT_MOUSE\n", 
	      stderr);
#endif
    }

    return (got_x && got_y);
}

static png_plot *png_plot_new (void)
{
    png_plot *plot;

    plot = mymalloc(sizeof *plot);
    if (plot == NULL) {
	return NULL;
    }

    plot->shell = NULL;
    plot->canvas = NULL;
    plot->popup = NULL;
    plot->statusarea = NULL;    
    plot->statusbar = NULL;
    plot->cursor_label = NULL;
    plot->pixmap = NULL;
    plot->invert_gc = NULL;
    plot->spec = NULL;
    plot->xmin = plot->xmax = 0.0;
    plot->ymin = plot->ymax = 0.0;
    plot->xint = plot->yint = 0;
    plot->pd = 0;
    plot->cid = 0;
    plot->zoom = NULL;
    plot->status = 0;
    plot->format = 0;
#ifndef OLD_GTK
    plot->labeled = NULL;
#endif

    return plot;
}

int gnuplot_show_png (const char *plotfile, GPT_SPEC *spec, int saved)
{
    png_plot *plot;
    int plot_has_xrange;

    GtkWidget *vbox;
    GtkWidget *canvas_hbox;
    GtkWidget *label_frame = NULL;
    GtkWidget *status_hbox = NULL;

    plot = png_plot_new();
    if (plot == NULL) {
	return 1;
    }

    if (spec != NULL) {
	plot->spec = spec;
    } else {
	plot->spec = plotspec_new();
	if (plot->spec == NULL) {
	    free(plot);
	    return 1;
	}
	strcpy(plot->spec->fname, plotfile);
    }

    /* make png plot struct accessible via spec */
    plot->spec->ptr = plot;

    plot->zoom = mymalloc(sizeof *plot->zoom);
    if (plot->zoom == NULL) {
	if (spec == NULL) {
	    free_plotspec(plot->spec);
	}
	free(plot);
	return 1;
    }

    if (saved) {
	plot->status |= PLOT_SAVED;
    }

    /* parse this file for x range */
    plot_has_xrange = get_plot_ranges(plot);

#ifdef OLD_GTK
    gtk_widget_push_visual(gdk_rgb_get_visual());
    gtk_widget_push_colormap(gdk_rgb_get_cmap());
    plot->shell = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_widget_pop_visual();
    gtk_widget_pop_colormap();
#else
    plot->shell = gtk_window_new(GTK_WINDOW_TOPLEVEL);
#endif

    gtk_widget_ref(plot->shell);

    gtk_window_set_title(GTK_WINDOW(plot->shell), _("gretl: gnuplot graph")); 
#ifdef OLD_GTK
    gtk_window_set_policy(GTK_WINDOW(plot->shell), FALSE, FALSE, FALSE);
#else
    gtk_window_set_resizable(GTK_WINDOW(plot->shell), FALSE);
#endif

    vbox = gtk_vbox_new(FALSE, 2);
    gtk_container_add(GTK_CONTAINER(plot->shell), vbox);

    g_signal_connect(G_OBJECT(plot->shell), "destroy",
		     G_CALLBACK(destroy_png_plot), plot);
    g_signal_connect(G_OBJECT(plot->shell), "key_press_event", 
		     G_CALLBACK(plot_key_handler), plot);

    /* box to hold canvas */
    canvas_hbox = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), canvas_hbox, TRUE, TRUE, 0);
    gtk_widget_show(canvas_hbox);

    /*  eventbox and hbox for status area  */
    plot->statusarea = gtk_event_box_new();
    gtk_box_pack_start(GTK_BOX(vbox), plot->statusarea, FALSE, FALSE, 0);

    status_hbox = gtk_hbox_new (FALSE, 2);
    gtk_container_add(GTK_CONTAINER(plot->statusarea), status_hbox);
    gtk_widget_show (status_hbox);
    gtk_container_set_resize_mode (GTK_CONTAINER (status_hbox),
				   GTK_RESIZE_QUEUE);

    /* Create drawing-area widget */
    plot->canvas = gtk_drawing_area_new();
#ifdef OLD_GTK
    gtk_drawing_area_size(GTK_DRAWING_AREA(plot->canvas), 
			  PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT);
#else
    gtk_widget_set_size_request(GTK_WIDGET(plot->canvas), 
				PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT);
#endif
    gtk_widget_set_events (plot->canvas, GDK_EXPOSURE_MASK
                           | GDK_LEAVE_NOTIFY_MASK
                           | GDK_BUTTON_PRESS_MASK
                           | GDK_BUTTON_RELEASE_MASK
                           | GDK_POINTER_MOTION_MASK
                           | GDK_POINTER_MOTION_HINT_MASK);

    GTK_WIDGET_SET_FLAGS (plot->canvas, GTK_CAN_FOCUS);

    g_signal_connect(G_OBJECT(plot->canvas), "button_press_event", 
		     G_CALLBACK(plot_button_press), plot);
    g_signal_connect(G_OBJECT(plot->canvas), "button_release_event", 
		     G_CALLBACK(plot_button_release), plot);

    /* create the contents of the status area */
    if (plot_has_xrange) {
	/* cursor label (graph position indicator) */
	label_frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(label_frame), GTK_SHADOW_IN);

	plot->cursor_label = gtk_label_new(" ");
	gtk_container_add(GTK_CONTAINER(label_frame), plot->cursor_label);
	gtk_widget_show(plot->cursor_label);
    }

    /* the statusbar */
    plot->statusbar = gtk_statusbar_new();
#ifdef OLD_GTK
    gtk_widget_set_usize(plot->statusbar, 1, -1);
#else
    gtk_widget_set_size_request(plot->statusbar, 1, -1);
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(plot->statusbar), FALSE);
#endif
    gtk_container_set_resize_mode(GTK_CONTAINER (plot->statusbar),
				  GTK_RESIZE_QUEUE);
    plot->cid = gtk_statusbar_get_context_id (GTK_STATUSBAR (plot->statusbar),
					      "plot_message");
    gtk_statusbar_push (GTK_STATUSBAR (plot->statusbar),
			plot->cid, _(" Click on graph for pop-up menu"));
    
    if (plot_has_xrange) {
	g_signal_connect(G_OBJECT(plot->canvas), "motion_notify_event",
			 G_CALLBACK(motion_notify_event), plot);
    }

    /* pack the widgets */
    gtk_box_pack_start(GTK_BOX(canvas_hbox), plot->canvas, FALSE, FALSE, 0);

    /* fill the status area */
    if (plot_has_xrange) {
	gtk_box_pack_start(GTK_BOX(status_hbox), label_frame, FALSE, FALSE, 0);
    }

    gtk_box_pack_start(GTK_BOX(status_hbox), plot->statusbar, TRUE, TRUE, 0); 

    /* show stuff */
    gtk_widget_show(plot->canvas);

    if (plot_has_xrange) {
	gtk_widget_show(label_frame);
    }

    gtk_widget_show(plot->statusbar);
    gtk_widget_show(plot->statusarea);

    gtk_widget_realize (plot->canvas);
    gdk_window_set_back_pixmap (plot->canvas->window, NULL, FALSE);

    if (plot_has_xrange) {
	gtk_widget_realize (plot->cursor_label);
#ifdef OLD_GTK
	gtk_widget_set_usize (plot->cursor_label, 140, -1);
#else
	gtk_widget_set_size_request (plot->cursor_label, 140, -1);
#endif
    }

    gtk_widget_show(vbox);
    gtk_widget_show(plot->shell);       

    /* set the focus to the canvas area */
    gtk_widget_grab_focus(plot->canvas);  

    plot->pixmap = gdk_pixmap_new(plot->shell->window, 
				  PLOT_PIXEL_WIDTH, PLOT_PIXEL_HEIGHT, 
				  -1);
    g_signal_connect(G_OBJECT(plot->canvas), "expose_event",
		     G_CALLBACK(plot_expose), plot->pixmap);

    render_pngfile(plot, PNG_START);

    return 0;
}

/* apparatus for getting coordinate info out of PNG files created using
   Allin Cottrell's modified version of gnuplot, which writes such info
   into the PNG comment fields
*/

#ifdef PNG_COMMENTS

static int get_png_plot_bounds (const char *str, png_bounds *bounds)
{
    int ret = GRETL_PNG_OK;

    if (sscanf(str, "xleft=%d xright=%d ybot=%d ytop=%d", 
	       &bounds->xleft, &bounds->xright,
	       &bounds->ybot, &bounds->ytop) != 4) {
	ret = GRETL_PNG_BAD_COMMENTS;
    } else if (bounds->xleft == 0 && bounds->xright == 0 &&
	       bounds->ybot == 0 && bounds->ytop == 0) {
	ret = GRETL_PNG_NO_COORDS;
    } 

    return ret;
}

static int get_png_data_bounds (char *str, png_bounds *bounds)
{
    char *p = str;
    int ret = GRETL_PNG_OK;

    while (*p) {
	if (*p == ',') *p = '.';
	p++;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    if (sscanf(str, "xmin=%lf xmax=%lf ymin=%lf ymax=%lf", 
	       &bounds->xmin, &bounds->xmax,
	       &bounds->ymin, &bounds->ymax) != 4) {
	ret = GRETL_PNG_BAD_COMMENTS;
    } else if (bounds->xmin == 0.0 && bounds->xmax == 0.0 &&
	       bounds->ymin == 0.0 && bounds->ymax == 0.0) {
	ret = GRETL_PNG_NO_COORDS;
    } 

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    return ret;
}

#define PNG_CHECK_BYTES 4

static int get_png_bounds_info (png_bounds *bounds)
{
    FILE *fp;
    char header[PNG_CHECK_BYTES];
    char pngname[MAXLEN];
    png_structp png_ptr;
    png_infop info_ptr;
    png_text *text_ptr = NULL;
    int i, num_text;
    volatile int ret = GRETL_PNG_OK;

    build_path(paths.userdir, "gretltmp.png", pngname, NULL); 

    fp = fopen(pngname, "rb");
    if (fp == NULL) {
	return GRETL_PNG_NO_OPEN;
    }

    fread(header, 1, PNG_CHECK_BYTES, fp);

    if (png_sig_cmp(header, 0, PNG_CHECK_BYTES)) {
	fclose(fp);
	sprintf(errtext, "Bad PNG header: Got bytes %x %x %x %x", 
		header[0],header[1],header[2],header[3]);
	errbox(errtext);
	return GRETL_PNG_NOT_PNG;
    }

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 
				     NULL, NULL, NULL);
    if (png_ptr == NULL) {
	fclose(fp);
	return GRETL_PNG_NO_OPEN;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        png_destroy_read_struct(&png_ptr, (png_infopp) NULL, 
				(png_infopp) NULL);
	fclose(fp);
        return GRETL_PNG_NO_OPEN;
    }

    if (setjmp(png_ptr->jmpbuf)) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fclose(fp);
        return GRETL_PNG_NO_OPEN;
    }

    png_init_io(png_ptr, fp);

    png_set_sig_bytes(png_ptr, PNG_CHECK_BYTES);
    png_read_info(png_ptr, info_ptr);

    num_text = png_get_text(png_ptr, info_ptr, &text_ptr, &num_text);

    if (num_text > 1) {
	int plot_ret = -1, data_ret = -1;

	for (i=1; i<num_text; i++) {
	    if (!strcmp(text_ptr[i].key, "plot bounds")) {
		plot_ret = get_png_plot_bounds(text_ptr[i].text, bounds);
	    }
	    if (!strcmp(text_ptr[i].key, "data bounds")) {
		data_ret = get_png_data_bounds(text_ptr[i].text, bounds);
	    }
	}
	if (plot_ret == GRETL_PNG_NO_COORDS && data_ret == GRETL_PNG_NO_COORDS) {
	    /* comments were present and correct, but all zero */
	    ret = GRETL_PNG_NO_COORDS;
	}
	else if (plot_ret != GRETL_PNG_OK || data_ret != GRETL_PNG_OK) {
	    /* one or both set of coordinates bad or missing */
	    if (plot_ret >= 0 || data_ret >= 0) {
		ret = GRETL_PNG_BAD_COMMENTS;
	    } else {
		ret = GRETL_PNG_NO_COMMENTS;
	    }
	}
    } else {
	/* no coordinates comments present */
	ret = GRETL_PNG_NO_COMMENTS;
    }

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    fclose(fp);
    
    return ret;
}

#endif /* PNG_COMMENTS */

#ifdef G_OS_WIN32

/* win32: copy plot to clipboard by generating an EMF file (enhanced
   metafile), reading it into a buffer, and putting it on the
   clipboard.

   Weirdness: when an emf is put on the clipboard as below, Word 2000
   behaves thus: a straight "Paste" puts in a version of the graph
   with squashed up numbers on the axes and no legend text; but a
   "Paste special" (where one accepts the default of pasting it as an
   enhanced metafile) puts in an accurate version with correct text.
   Go figure.  (This is on win98)
*/

#include <gdk/gdkwin32.h>
#include "guiprint.h"

static int emf_to_clip (char *emfname)
{
    HWND mainw;
    HENHMETAFILE hemf, hemfclip;
    HANDLE htest;

    mainw = GDK_WINDOW_HWND(mdata->w->window);
    if (mainw == NULL) {
	errbox("Got NULL HWND");
	return 1;
    }	

    if (!OpenClipboard(mainw)) {
	errbox(_("Cannot open the clipboard"));
	return 1;
    }

    EmptyClipboard();

    hemf = GetEnhMetaFile(emfname);
    if (hemf == NULL) {
	errbox("Couldn't get handle to graphic metafile");
	return 1;
    }

    hemfclip = CopyEnhMetaFile(hemf, NULL);
    if (hemfclip == NULL) {
	errbox("Couldn't copy graphic metafile");
	return 1;
    }    

    htest = SetClipboardData(CF_ENHMETAFILE, hemfclip);
    if (htest == NULL) {
	errbox("Failed to put data on clipboard");
	return 1;
    }  	

    CloseClipboard();

    DeleteEnhMetaFile(hemf);

    return 0;
}

static void win32_process_graph (GPT_SPEC *spec, int color, int dest)
{
    FILE *fq;
    PRN *prn;
    char plottmp[MAXLEN], plotline[MAXLEN];
    gchar *plotcmd = NULL;
    gchar *emfname = NULL;
    int err;

    /* create temporary file to hold the special gnuplot commands */
    if (!user_fopen("gptout.tmp", plottmp, &prn)) return;

    /* open the gnuplot source file for the graph */
    fq = fopen(spec->fname, "r");
    if (fq == NULL) {
	errbox(_("Couldn't access graph info"));
	gretl_print_destroy(prn);
	return;
    }

    /* generate gnuplot source file to make emf */
    pprintf(prn, "%s\n", get_gretl_emf_term_line(spec->code, color));
    emfname = g_strdup_printf("%sgpttmp.emf", paths.userdir);
    pprintf(prn, "set output '%s'\n", emfname);
    pprintf(prn, "set size 0.8,0.8\n");
    while (fgets(plotline, MAXLEN-1, fq)) {
	if (strncmp(plotline, "set term", 8) && 
	    strncmp(plotline, "set output", 10))
	    pprintf(prn, "%s", plotline);
    }

    gretl_print_destroy(prn);
    fclose(fq);

    /* get gnuplot to create the emf file */
    plotcmd = g_strdup_printf("\"%s\" \"%s\"", paths.gnuplot, 
			      plottmp);
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
    g_free(plotcmd);
    remove(plottmp);
    
    if (err) {
        errbox(_("Gnuplot error creating graph"));
    } else if (dest == WIN32_TO_CLIPBOARD) {
	err = emf_to_clip(emfname);
	if (!err) {
	    infobox(_("To paste, use Edit/Paste special.../Enhanced metafile"));
	}
    } else if (dest == WIN32_TO_PRINTER) {
	err = winprint_graph(emfname);
    }

    remove(emfname);
    g_free(emfname);
}

#endif /* G_OS_WIN32 */

