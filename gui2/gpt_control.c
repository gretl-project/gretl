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

/* gpt_control.c for gretl -- gnuplot controller */

#include "gretl.h"
#include "plotspec.h"
#include "gpt_control.h"
#include "session.h"
#include "gpt_dialog.h"
#include "fileselect.h"
#include "calculator.h"
#include "guiprint.h"
#include "textbuf.h"
#include "graphics.h"
#include "boxplots.h"
#include "dlgutils.h"
#include "winstack.h"

#define GPDEBUG 0
#define POINTS_DEBUG 0

#ifdef G_OS_WIN32
# include <io.h>
# include "gretlwin32.h"
#endif

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk/gdkkeysyms.h>

enum {
    PLOT_SAVED          = 1 << 0,
    PLOT_ZOOMED         = 1 << 1,
    PLOT_ZOOMING        = 1 << 2,
    PLOT_NO_MARKERS     = 1 << 3,
    PLOT_PNG_COORDS     = 1 << 4,
    PLOT_HAS_XRANGE     = 1 << 5,
    PLOT_HAS_YRANGE     = 1 << 6,
    PLOT_DONT_ZOOM      = 1 << 7,
    PLOT_DONT_EDIT      = 1 << 8,
    PLOT_DONT_MOUSE     = 1 << 9,
    PLOT_POSITIONING    = 1 << 10,
    PLOT_CURSOR_LABEL   = 1 << 11
} plot_status_flags;

enum {
    PLOT_TITLE          = 1 << 0,
    PLOT_XLABEL         = 1 << 1,
    PLOT_YLABEL         = 1 << 2,
    PLOT_Y2AXIS         = 1 << 3,
    PLOT_Y2LABEL        = 1 << 4,
    PLOT_MARKERS_UP     = 1 << 5,
    PLOT_POLAR          = 1 << 6
} plot_format_flags;

#define MAX_MARKERS 250

#define plot_is_zoomed(p)       (p->status & PLOT_ZOOMED)
#define plot_is_zooming(p)      (p->status & PLOT_ZOOMING)
#define plot_has_png_coords(p)  (p->status & PLOT_PNG_COORDS)
#define plot_has_xrange(p)      (p->status & PLOT_HAS_XRANGE)
#define plot_has_yrange(p)      (p->status & PLOT_HAS_YRANGE)
#define plot_not_editable(p)    (p->status & PLOT_DONT_EDIT)
#define plot_is_editable(p)     (!(p->status & PLOT_DONT_EDIT))
#define plot_doing_position(p)  (p->status & PLOT_POSITIONING)

#define plot_has_title(p)        (p->format & PLOT_TITLE)
#define plot_has_xlabel(p)       (p->format & PLOT_XLABEL)
#define plot_has_ylabel(p)       (p->format & PLOT_YLABEL)
#define plot_has_y2axis(p)       (p->format & PLOT_Y2AXIS)
#define plot_has_y2label(p)      (p->format & PLOT_Y2LABEL)
#define plot_labels_shown(p)     (p->format & PLOT_MARKERS_UP)
#define plot_is_polar(p)         (p->format & PLOT_POLAR)

#define plot_is_range_mean(p)   (p->spec->code == PLOT_RANGE_MEAN)
#define plot_is_hurst(p)        (p->spec->code == PLOT_HURST)
#define plot_is_roots(p)        (p->spec->code == PLOT_ROOTS)

#define plot_has_regression_list(p) (p->spec->reglist != NULL)

#define labels_frozen(p)        (p->spec->flags & GPT_PRINT_MARKERS)
#define cant_do_labels(p)       (p->status & PLOT_NO_MARKERS)

#define plot_show_cursor_label(p) (p->status & PLOT_CURSOR_LABEL)

#define plot_has_controller(p) (p->editor != NULL)

enum {
    PNG_START,
    PNG_ZOOM,
    PNG_UNZOOM,
    PNG_REDISPLAY
} png_zoom_codes;

struct png_plot_t {
    GtkWidget *shell;
    GtkWidget *canvas;
    GtkWidget *popup;
    GtkWidget *statusarea;    
    GtkWidget *statusbar;
    GtkWidget *cursor_label;
    GtkWidget *pos_entry;
    GtkWidget *editor;
    GtkWidget *up_icon;
    GtkWidget *down_icon;
    GdkWindow *window;
    cairo_t *cr;
#if GTK_MAJOR_VERSION >= 3
    cairo_surface_t *cs;
#else
    GdkPixmap *pixmap;
    GdkPixbuf *savebuf;
#endif
    GPT_SPEC *spec;
    double xmin, xmax;
    double ymin, ymax;
    int pixel_width, pixel_height;
    int pixel_xmin, pixel_xmax;
    int pixel_ymin, pixel_ymax;
    int xint, yint;
    int pd;
    int err;
    guint cid;
    double zoom_xmin, zoom_xmax;
    double zoom_ymin, zoom_ymax;
    int screen_x0, screen_y0; /* to define selection */
    unsigned long status; 
    unsigned char format;
};

static int render_pngfile (png_plot *plot, int view);
static int repaint_png (png_plot *plot, int view);
static int zoom_replaces_plot (png_plot *plot);
static void prepare_for_zoom (png_plot *plot);
static int get_plot_ranges (png_plot *plot, PlotType ptype);
static void graph_display_pdf (GPT_SPEC *spec);
#ifdef G_OS_WIN32
static void win32_process_graph (GPT_SPEC *spec, int dest);
#endif
static void plot_do_rescale (png_plot *plot, int mod);

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

typedef struct linestyle_ linestyle;

struct linestyle_ {
    char rgb[8];
    int type;
};

#define MAX_STYLES N_GP_COLORS

static int get_png_bounds_info (png_bounds *bounds);

#define PLOTSPEC_DETAILS_IN_MEMORY(s) (s->data != NULL)

static void terminate_plot_positioning (png_plot *plot)
{
    if (plot->status & PLOT_POSITIONING) {
	plot->status ^= PLOT_POSITIONING;
	plot->pos_entry = NULL;
	gdk_window_set_cursor(plot->window, NULL);
	gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
	if (plot->editor != NULL) {
	    gtk_window_present(GTK_WINDOW(plot->editor));
	}
    }
}

GtkWidget *plot_get_shell (png_plot *plot) 
{
    return plot->shell;
}

GPT_SPEC *plot_get_spec (png_plot *plot) 
{
    return plot->spec;
}

int plot_is_saved (const png_plot *plot)
{
    return (plot->status & PLOT_SAVED);
}

int plot_is_mouseable (const png_plot *plot)
{
    return !(plot->status & PLOT_DONT_MOUSE);
}

double plot_get_xmin (png_plot *plot)
{
    return (plot != NULL)? plot->xmin : -1;
}

double plot_get_ymin (png_plot *plot)
{
    return (plot != NULL)? plot->ymin : -1;
}

/* apparatus for graph toolbar */

static void gp_button_sizing (GtkWidget *w)
{
    static int style_done;

    gtk_widget_set_name(w, "gp_button");

    if (!style_done) {
	gtk_rc_parse_string("style \"gp-style\"\n{\n"
			    "  GtkWidget::focus-padding = 1\n"
			    "  GtkWidget::focus-line-width = 0\n"
			    "  xthickness = 2\n"
			    "  ythickness = 1\n"
			    "}\n"
			    "widget \"*.gp_button\" style \"gp-style\"");
	style_done = 1;
    }
}

static GtkWidget *small_tool_button (GretlToolItem *item,
				     png_plot *plot)
{
    GtkWidget *img, *button = gtk_button_new();

    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    gp_button_sizing(button);
    img = gtk_image_new_from_stock(item->icon, GTK_ICON_SIZE_MENU);
    gtk_container_add(GTK_CONTAINER(button), img);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button), _(item->tip));
    if (!strcmp(item->icon, GRETL_STOCK_WINLIST)) {
	g_signal_connect(G_OBJECT(button), "button-press-event", 
			 item->func, plot->shell);
    } else {
	g_signal_connect(G_OBJECT(button), "clicked", 
			 item->func, plot);
    }
    gtk_widget_show_all(button);

    return button;
}

static void graph_enlarge_callback (GtkWidget *w, png_plot *plot)
{
    plot_do_rescale(plot, 1);
}

static void graph_shrink_callback (GtkWidget *w, png_plot *plot)
{
    plot_do_rescale(plot, -1);
}

static void graph_edit_callback (GtkWidget *w, png_plot *plot)
{
    start_editing_png_plot(plot);
}

static void graph_zoom_callback (GtkWidget *w, png_plot *plot)
{
    prepare_for_zoom(plot);
}

static void show_pdf_callback (GtkWidget *w, png_plot *plot)
{
    graph_display_pdf(plot->spec);
}

static void close_plot_callback (GtkWidget *w, png_plot *plot)
{
    gtk_widget_destroy(plot->shell);
}

static GretlToolItem plotbar_items[] = {
    { N_("Close"),       GTK_STOCK_CLOSE,     G_CALLBACK(close_plot_callback), 0 },
    { N_("Windows"),     GRETL_STOCK_WINLIST, G_CALLBACK(window_list_popup), 0 },
    { N_("Edit"),        GTK_STOCK_EDIT,      G_CALLBACK(graph_edit_callback), 0 },
    { N_("Zoom..."),     GTK_STOCK_ZOOM_IN,   G_CALLBACK(graph_zoom_callback), 0 },
    { N_("Display PDF"), GRETL_STOCK_PDF,     G_CALLBACK(show_pdf_callback), 0 },
    { N_("Smaller"),     GRETL_STOCK_SMALLER, G_CALLBACK(graph_shrink_callback), 0 },
    { N_("Bigger"),      GRETL_STOCK_BIGGER,  G_CALLBACK(graph_enlarge_callback), 0 }
};

static void add_graph_toolbar (GtkWidget *hbox, png_plot *plot)
{
    GtkWidget *button;
    GretlToolItem *item;
    int i, n = G_N_ELEMENTS(plotbar_items);

    for (i=0; i<n; i++) {
	item = &plotbar_items[i];
	if (item->func == G_CALLBACK(graph_edit_callback) &&
	    plot_not_editable(plot)) {
	    continue;
	}
	if (item->func == G_CALLBACK(graph_zoom_callback)) {
	    /* not ready yet */
	    continue;
	}
	if ((item->func == G_CALLBACK(graph_enlarge_callback) ||
	     item->func == G_CALLBACK(graph_shrink_callback)) &&
	    (plot_not_editable(plot) || 
	     gnuplot_png_terminal() != GP_PNG_CAIRO)) {
	    continue;
	}
	button = small_tool_button(item, plot);
	gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	if (item->func == G_CALLBACK(graph_enlarge_callback)) {
	    plot->up_icon = button;
	} else if (item->func == G_CALLBACK(graph_shrink_callback)) {
	    plot->down_icon = button;
	}
    }	
}

/* end apparatus for graph toolbar */

/* Provide the data coordinates for a gretl/gnuplot
   graph, if they are all positive, otherwise return 
   non-zero.
*/

int plot_get_coordinates (png_plot *plot,
			  double *xmin,
			  double *xmax,
			  double *ymin,
			  double *ymax)
{
    int err = 0;

    if (plot != NULL && plot->xmin > 0 && 
	plot->xmax > plot->xmin &&
	plot->ymax > plot->ymin) {
	*xmin = plot->xmin;
	*xmax = plot->xmax;
	*ymin = plot->ymin;
	*ymax = plot->ymax;
    } else {
	err = E_DATA;
	gretl_errmsg_set("Couldn't get plot coordinates");
    }

    return err;
}

void set_plot_has_y2_axis (png_plot *plot, gboolean s)
{
    if (s == TRUE) {
	plot->format |= PLOT_Y2AXIS;
    } else {
	plot->format &= ~PLOT_Y2AXIS;
    }
}

void plot_position_click (GtkWidget *w, png_plot *plot)
{
    if (plot != NULL) {
	GtkWidget *entry;
	GdkCursor* cursor;

	cursor = gdk_cursor_new(GDK_CROSSHAIR);
	gdk_window_set_cursor(plot->window, cursor);
	gdk_cursor_unref(cursor);
	entry = g_object_get_data(G_OBJECT(w), "pos_entry");
	plot->pos_entry = entry;
	plot->status |= PLOT_POSITIONING;
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid, 
			   _(" Click to set position"));
    }
}

static FILE *open_gp_file (const char *fname, const char *mode)
{
    FILE *fp = gretl_fopen(fname, mode);

    if (fp == NULL) {
	if (*mode == 'w') {
	    file_write_errbox(fname);
	} else {
	    file_read_errbox(fname);
	}
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

static int set_print_line (const char *s)
{
    return (!strncmp(s, "set print ", 10) ||
	    !strncmp(s, "print \"pixe", 11) ||
	    !strncmp(s, "print \"data", 11));
}

enum {
    REMOVE_PNG,
    ADD_PNG
};

static int is_png_term_line (const char *s)
{
    return !strncmp(s, "set term png", 12);
}

#ifdef G_OS_WIN32

static void check_win32_png_spec (char *s)
{
    if (!strncmp("set term png ", s, 13)) {
	/* should now be pngcairo */
	*s = '\0';
    }
}

#endif

static int 
add_or_remove_png_term (const char *fname, int action, GPT_SPEC *spec)
{
    FILE *fsrc, *ftmp;
    char temp[MAXLEN], fline[MAXLEN];
    GptFlags flags = 0;

    sprintf(temp, "%sgpttmp", gretl_dotdir());
    ftmp = gretl_tempfile_open(temp);
    if (ftmp == NULL) {
	return 1;
    }

    fsrc = open_gp_file(fname, "r");
    if (fsrc == NULL) {
	fclose(ftmp);
	return 1;
    }
    
    if (action == ADD_PNG) {
	/* see if there's already a png term setting, possibly commented
	   out, that can be reused */
	char restore_line[MAXLEN];
	int add_line_styles = 1;

	*restore_line = '\0';

	while (fgets(fline, sizeof fline, fsrc)) {
	    if (is_png_term_line(fline) && *restore_line == '\0') {
		strcat(restore_line, fline);
	    } else if (commented_term_line(fline) && *restore_line == '\0') {
		strcat(restore_line, fline + 2);
	    } else if (strstr(fline, "letterbox")) {
		flags = GPT_LETTERBOX;
	    } else if (strstr(fline, "large")) {
		flags = GPT_XL;
	    } else if (strstr(fline, "extra-large")) {
		flags = GPT_XXL;
	    } else if (!strncmp(fline, "set style line", 14)) {
		add_line_styles = 0;
	    } else if (!strncmp(fline, "plot", 4)) {
		break;
	    }
	}

	rewind(fsrc);

#ifdef G_OS_WIN32
	/* check for obsolete png term specification (as may be found
	   in an old session file) */
	check_win32_png_spec(restore_line);
#endif

	if (*restore_line != '\0') {
	    fputs(restore_line, ftmp);
	} else if (spec != NULL) {
	    fprintf(ftmp, "%s\n", get_png_line_for_plotspec(spec));
	} else {
	    fprintf(ftmp, "%s\n", get_gretl_png_term_line(PLOT_REGULAR, flags));
	}
	write_plot_output_line(NULL, ftmp);

	if (add_line_styles) {
	    write_plot_line_styles(PLOT_REGULAR, ftmp);
	}

	/* now for the body of the plot file */

	while (fgets(fline, sizeof fline, fsrc)) {
	    if (set_print_line(fline)) {
		; /* skip it (portability) */
	    } else if (is_png_term_line(fline)) {
		; /* handled above */
	    } else if (commented_term_line(fline)) {
		; /* handled above */
	    } else if (set_output_line(fline)) {
		; /* handled above */
	    } else {
		fputs(fline, ftmp);
	    }
	}
	write_plot_bounding_box_request(ftmp);
    } else {
	/* not ADD_PNG: we're removing the png term line */
	int printit, png_line_saved = 0;
	
	while (fgets(fline, sizeof fline, fsrc)) {
	    printit = 1;
	    if (is_png_term_line(fline)) {
		if (!png_line_saved) {
		    /* comment it out, for future reference */
		    fprintf(ftmp, "# %s", fline);
		    png_line_saved = 1;
		} 
		printit = 0;
	    } else if (commented_term_line(fline)) {
		if (png_line_saved) {
		    printit = 0;
		}
	    } else if (set_output_line(fline)) {
		printit = 0;
	    } else if (spec != NULL && (spec->flags & GPT_FIT_HIDDEN)
		       && is_auto_fit_string(fline)) {
		printit = 0;
	    } else if (set_print_line(fline)) {
		printit = 0;
	    }
	    if (printit) {
		fputs(fline, ftmp);
	    }
	}
    }

    fclose(fsrc);
    fclose(ftmp);
    gretl_remove(fname);

    return gretl_rename(temp, fname);
}

static int add_png_term_to_plot (const char *fname)
{
    return add_or_remove_png_term(fname, ADD_PNG, NULL);
}

static int remove_png_term_from_plot (const char *fname, GPT_SPEC *spec)
{
    return add_or_remove_png_term(fname, REMOVE_PNG, spec);
}

/* public because called from session.c when editing plot commands */

int remove_png_term_from_plot_by_name (const char *fname)
{
    return add_or_remove_png_term(fname, REMOVE_PNG, NULL);
}

static void mark_plot_as_saved (GPT_SPEC *spec)
{
    png_plot *plot = (png_plot *) spec->ptr;

    plot->status |= PLOT_SAVED;
}

static int gnuplot_png_init (png_plot *plot, FILE **fpp)
{
    GPT_SPEC *spec = plot->spec;
    char fname[FILENAME_MAX];

    if (plot_is_saved(plot)) {
	/* session graph: ensure we're writing to the correct
	   directory */
	session_graph_make_path(fname, spec->fname);
    } else {
	strcpy(fname, spec->fname);
    }

    *fpp = gretl_fopen(fname, "w");

    if (*fpp == NULL) {
	file_write_errbox(fname);
	return 1;
    }

    fprintf(*fpp, "%s\n", get_png_line_for_plotspec(spec));
    write_plot_output_line(NULL, *fpp);

    return 0;
}

int gp_term_code (gpointer p, int action)
{
    GPT_SPEC *spec;
    
    if (action == SAVE_GNUPLOT) {
	spec = (GPT_SPEC *) p;
    } else {
	/* EPS/PDF saver */
	spec = graph_saver_get_plotspec(p);
    }

    return spec->termtype;
}

static void get_full_term_string (const GPT_SPEC *spec, char *termstr) 
{
    if (spec->termtype == GP_TERM_EPS) {
	const char *s = get_gretl_eps_term_line(spec->code, spec->flags);

	strcpy(termstr, s);
    } else if (spec->termtype == GP_TERM_PDF) {
	const char *s = get_gretl_pdf_term_line(spec->code, spec->flags);

	strcpy(termstr, s);
    } else if (spec->termtype == GP_TERM_FIG) {
	strcpy(termstr, "set term fig");
    } else if (spec->termtype == GP_TERM_PNG) { 
	strcpy(termstr, get_png_line_for_plotspec(spec)); 
    } else if (spec->termtype == GP_TERM_EMF) {
	int mono = (spec->flags & GPT_MONO);

	strcpy(termstr, get_gretl_emf_term_line(spec->code, !mono));
    } 
}

/* Sometimes we want to put \pi into tic-marks or a graph
   title. We just use the UTF-8 character, but if we're writing
   a plot in EPS format this has to be changed to {/Symbol p}.
*/

static char *eps_replace_pi (unsigned char *s)
{
    char *repl = NULL;
    int i, picount = 0;

    for (i=0; s[i] != '\0'; i++) {
	if (s[i] == 0xcf && s[i+1] == 0x80) {
	    picount++;
	}
    }

    if (picount > 0) {
	i = picount * 9 + strlen((const char *) s) + 1;
	repl = calloc(i, 1);
    }

    if (repl != NULL) {
	i = 0;
	while (*s) {
	    if (*s == 0xcf && *(s+1) == 0x80) {
		strcat(repl, "{/Symbol p}");
		i += 11;
		s++;
	    } else {
		repl[i++] = *s;
	    }
	    s++;
	}
    }

    return repl;
}

/* for postscript output, e.g. in Latin-2, or EMF output in CPXXXX */

static int maybe_recode_gp_line (char *s, int ttype, FILE *fp)
{
    int err = 0;

    if (!gretl_is_ascii(s) && g_utf8_validate(s, -1, NULL)) {
	char *tmp;

	if (ttype == GP_TERM_EMF) {
	    tmp = utf8_to_cp(s);
	} else {
	    char *repl = eps_replace_pi((unsigned char *) s);

	    if (repl != NULL) {
		tmp = utf8_to_latin(repl);
		free(repl);
	    } else {
		tmp = utf8_to_latin(s);
	    }
	}

	if (tmp == NULL) {
	    err = 1;
	} else {
	    fputs(tmp, fp);
	    free(tmp);
	}
    } else {
	fputs(s, fp);
    }

    return err;
}

/* Check for non-ASCII strings in plot file: these may
   require special treatment. */

static int non_ascii_gp_file (FILE *fp)
{
    char pline[512];
    int ret = 0;

    while (fgets(pline, sizeof pline, fp)) {
	if (set_print_line(pline)) {
	    break;
	}
	if (*pline == '#') {
	    continue;
	}
	if (!gretl_is_ascii(pline)) {
	    ret = 1;
	    break;
	}
    }

    rewind(fp);

    return ret;
}

static int term_uses_utf8 (int ttype)
{
    if (ttype == GP_TERM_PNG || 
	ttype == GP_TERM_SVG ||
	ttype == GP_TERM_PLT) {
	return 1;
    } else if (ttype == GP_TERM_PDF && 
	       gnuplot_pdf_terminal() == GP_PDF_CAIRO) {
	return 1;
    } else {
	return 0;
    }
}

#define is_color_line(s) (strstr(s, "set style line") && strstr(s, "rgb"))

void filter_gnuplot_file (int ttype, int latin, int mono,
			  FILE *fpin, FILE *fpout)
{
    char pline[512];
    int err, err_shown = 0;

    while (fgets(pline, sizeof pline, fpin)) {
	if (set_print_line(pline)) {
	    break;
	}

	if (!strncmp(pline, "set term", 8) ||
	    !strncmp(pline, "set enco", 8) ||
	    !strncmp(pline, "set outp", 8)) {
	    continue;
	}

	if (mono) {
	    if (is_color_line(pline)) {
		continue;
	    } else if (strstr(pline, "set style fill solid")) {
		fputs("set style fill solid 0.3\n", fpout);
		continue;
	    }
	}

	if (latin && *pline != '#') {
	    err = maybe_recode_gp_line(pline, ttype, fpout);
	    if (err && !err_shown) {
		gui_errmsg(err);
		err_shown = 1;
	    }
	} else {
	    fputs(pline, fpout);
	} 
    }
}

/* for non-UTF-8 plot formats: print a "set encoding" string
   if appropriate, but only if gnuplot won't choke on it.
*/

static void maybe_print_gp_encoding (int ttype, int latin, FILE *fp)
{
    if (ttype == GP_TERM_EMF) {
	if (latin == 2) {
	    fputs("set encoding cp1250\n", fp);
	} else if (latin == 9) {
	    fputs("set encoding cp1254\n", fp);
	} else if (chinese_locale() && gnuplot_has_cp950()) {
	    fputs("set encoding cp950\n", fp);
	}
    } else if (latin == 1 || latin == 2 || latin == 15 || latin == 9) {
	/* supported by gnuplot >= 4.4.0 */
	fprintf(fp, "set encoding iso_8859_%d\n", latin);
    }
} 

static int revise_plot_file (GPT_SPEC *spec, 
			     const char *pltname,
			     const char *outtarg,
			     const char *setterm)
{
    FILE *fpin = NULL;
    FILE *fpout = NULL;
    int mono = (spec->flags & GPT_MONO);
    int ttype = spec->termtype;
    int latin = 0;
    int err = 0;

    fpin = gretl_fopen(spec->fname, "r");
    if (fpin == NULL) {
	file_read_errbox(spec->fname);
	return 1;
    }

    fpout = gretl_fopen(pltname, "w");
    if (fpout == NULL) {
	fclose(fpin);
	file_write_errbox(pltname);
	return 1;
    }

    if (non_ascii_gp_file(fpin)) {
	/* plot contains UTF-8 strings */
	if (!term_uses_utf8(ttype)) {
	    latin = iso_latin_version();
	    maybe_print_gp_encoding(ttype, latin, fpout);
	} else {
	    fputs("set encoding utf8\n", fpout);
	}
    }

    if (outtarg != NULL && *outtarg != '\0') {
	fprintf(fpout, "%s\n", setterm);
	write_plot_output_line(outtarg, fpout);
    }	

    filter_gnuplot_file(ttype, latin, mono, fpin, fpout);

    fclose(fpin);
    fclose(fpout);

    return err;
}

void save_graph_to_file (gpointer data, const char *fname)
{
    GPT_SPEC *spec = (GPT_SPEC *) data;
    char pltname[FILENAME_MAX];
    char setterm[MAXLEN];
    int err = 0;

    get_full_term_string(spec, setterm);

    sprintf(pltname, "%sgptout.tmp", gretl_dotdir());
    err = revise_plot_file(spec, pltname, fname, setterm);

    if (!err) {
	gchar *plotcmd;

	plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
				  gretl_gnuplot_path(), 
				  pltname);
	err = gretl_spawn(plotcmd);
	gretl_remove(pltname);
	g_free(plotcmd);
	if (err) {
	    gui_errmsg(err);
	} 
    }
}

#define GRETL_PDF_TMP "gretltmp.pdf"
#define GRETL_EPS_TMP "gretltmp.eps"

static void graph_display_pdf (GPT_SPEC *spec)
{
    char pdfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    char setterm[128];
    gchar *plotcmd;
    int err = 0;

    spec->termtype = GP_TERM_PDF;

    strcpy(setterm, get_gretl_pdf_term_line(spec->code, spec->flags));

    build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);
    build_path(pdfname, gretl_dotdir(), GRETL_PDF_TMP, NULL);

    err = revise_plot_file(spec, plttmp, pdfname, setterm);
    if (err) {
	return;
    }

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
			      gretl_gnuplot_path(), 
			      plttmp);
    err = gretl_spawn(plotcmd);
    gretl_remove(plttmp);
    g_free(plotcmd);

    if (err) {
	gui_errmsg(err);
	return;
    } 

#if defined(G_OS_WIN32)
    win32_open_file(pdfname);
#elif defined(OS_OSX)
    osx_open_file(pdfname);
#else
    gretl_fork("viewpdf", pdfname);
#endif
}

void saver_preview_graph (GPT_SPEC *spec, char *termstr)
{
    char grfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    gchar *plotcmd;
    int err = 0;

    build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);

    if (spec->termtype == GP_TERM_EPS) {
	build_path(grfname, gretl_dotdir(), GRETL_EPS_TMP, NULL);
    } else {
	build_path(grfname, gretl_dotdir(), GRETL_PDF_TMP, NULL);
    }

    err = revise_plot_file(spec, plttmp, grfname, termstr);
    if (err) {
	return;
    }

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", gretl_gnuplot_path(), 
			      plttmp);
    err = gretl_spawn(plotcmd);
    gretl_remove(plttmp);
    g_free(plotcmd);

    if (err) {
	gui_errmsg(err);
	return;
    } 

#if defined(G_OS_WIN32)
    win32_open_file(grfname);
#elif defined(OS_OSX)
    osx_open_file(grfname);
#else
    if (spec->termtype == GP_TERM_EPS) {
	gretl_fork("viewps", grfname);
    } else {
	gretl_fork("viewpdf", grfname);
    }
#endif
}

int saver_save_graph (GPT_SPEC *spec, char *termstr, const char *fname)
{
    char plttmp[FILENAME_MAX];
    int err;

    build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);

    err = revise_plot_file(spec, plttmp, fname, termstr);

    if (!err) {
	gchar *plotcmd;

	plotcmd = g_strdup_printf("\"%s\" \"%s\"", gretl_gnuplot_path(), 
			      plttmp);
	err = gretl_spawn(plotcmd);
	gretl_remove(plttmp);
	g_free(plotcmd);

	if (err) {
	    gui_errmsg(err);
	}
    }

    return err;
}

/* we're looking for an uncommented "set term png" or some such */

static int is_batch_term_line (const char *s)
{
    int ret = 0;

    while (isspace(*s)) s++;

    if (*s != '#') {
	s = strstr(s, "set term");

	if (s != NULL) {
	    s += 8;
	    s += strspn(s, "inal");
	    while (isspace(*s)) s++;
	    if (strncmp(s, "win", 3))
		ret = 1;
	}
    }

    return ret;
}

/* dump_plot_buffer: this is used when we're taking the material from
   an editor window containing gnuplot commands, and either (a)
   sending it to gnuplot for execution, or (b) saving it to a "user
   file".  

   This function handles the addition of "pause -1" on MS Windows,
   if @addpause is non-zero.
*/

int dump_plot_buffer (const char *buf, const char *fname,
		      int addpause)
{
    FILE *fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	file_write_errbox(fname);
	return E_FOPEN;
    }

    if (!addpause) {
	/* nice and simple! */
	fputs(buf, fp);
    } else {
#ifdef G_OS_WIN32
	int gotpause = 0;
#endif
	char bufline[512];

	bufgets_init(buf);

	while (bufgets(bufline, sizeof bufline, buf)) {
	    if (addpause && is_batch_term_line(bufline)) {
		addpause = 0;
	    }
	    fputs(bufline, fp);
#ifdef G_OS_WIN32
	    if (addpause && strstr(bufline, "pause -1")) {
		gotpause = 1;
	    }
#endif
	}

	bufgets_finalize(buf);

#ifdef G_OS_WIN32
	/* sending directly to gnuplot on MS Windows */
	if (addpause && !gotpause) {
	    fprintf(stderr, "adding 'pause -1'\n");
	    fputs("pause -1\n", fp);
	}
#endif
    }

    fclose(fp);

    return 0;
}

#ifdef G_OS_WIN32

static void real_send_to_gp (const char *tmpfile)
{
    gchar *cmd;
    int err;

    cmd = g_strdup_printf("\"%s\" \"%s\"", 
			  gretl_gnuplot_path(), 
			  tmpfile);
    err = (WinExec(cmd, SW_SHOWNORMAL) < 32);
    g_free(cmd);

    if (err) {
	win_show_last_error();
    }
}

#else

#include <sys/types.h>
#include <sys/wait.h>

static void real_send_to_gp (const char *tmpfile)
{
    GError *error = NULL;
    gchar *argv[4];
    GPid pid = 0;
    gint fd = -1;
    gboolean run;

    argv[0] = g_strdup(gretl_gnuplot_path());
    argv[1] = g_strdup("-persist");
    argv[2] = g_strdup(tmpfile);
    argv[3] = NULL;

    run = g_spawn_async_with_pipes(NULL, argv, NULL, 
				   G_SPAWN_SEARCH_PATH | 
				   G_SPAWN_DO_NOT_REAP_CHILD,
				   NULL, NULL, &pid, NULL, NULL,
				   &fd, &error);

    if (error != NULL) {
	errbox(error->message);
	g_error_free(error);
    } else if (!run) {
	errbox(_("gnuplot command failed"));
    } else if (pid > 0) {
	int status = 0, err = 0;

	/* bodge below: try to give gnuplot time to bomb
	   out, if it's going to -- but we don't want to
	   hold things up if we're doing OK */
	
	sleep(1);
	waitpid(pid, &status, WNOHANG);
	if (WIFEXITED(status)) {
	    err = WEXITSTATUS(status);
	} 

	if (err && fd > 0) {
	    char buf[128] = {0};

	    if (read(fd, buf, 127) > 0) {
		errbox(buf);
	    }
	}
    }

    if (fd > 0) {
	close(fd);
    }

    g_spawn_close_pid(pid);

    g_free(argv[0]);
    g_free(argv[1]);
    g_free(argv[2]);
}

#endif

/* callback for execute icon in window editing gnuplot
   commands */

void run_gp_script (gchar *buf)
{
#ifdef G_OS_WIN32
    int addpause = 1;
#else
    int addpause = 0;
#endif
    gchar *tmpfile;
    int err = 0;

    tmpfile = g_strdup_printf("%showtmp.gp", gretl_dotdir());
    err = dump_plot_buffer(buf, tmpfile, addpause);

    if (!err) {
	real_send_to_gp(tmpfile);
    }   

    gretl_remove(tmpfile);
    g_free(tmpfile);
}

#ifdef G_OS_WIN32

/* common code for sending an EMF file to the clipboard,
   or printing an EMF, on MS Windows */

static void win32_process_graph (GPT_SPEC *spec, int dest)
{
    char emfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    gchar *plotcmd;
    const char *setterm;
    int color, err = 0;

    build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);
    build_path(emfname, gretl_dotdir(), "gpttmp.emf", NULL);

    color = !(spec->flags & GPT_MONO);

    setterm = get_gretl_emf_term_line(spec->code, color);
    
    err = revise_plot_file(spec, plttmp, emfname, setterm);
    if (err) {
	return;
    }

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
			      gretl_gnuplot_path(), 
			      plttmp);
    err = gretl_spawn(plotcmd);
    g_free(plotcmd);
    gretl_remove(plttmp);
    
    if (err) {
        errbox(_("Gnuplot error creating graph"));
    } else if (dest == WIN32_TO_CLIPBOARD) {
	err = emf_to_clipboard(emfname);
    } else if (dest == WIN32_TO_PRINTER) {
	err = win32_print_graph(emfname);
    }

    gretl_remove(emfname);
}

#endif

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

static int get_gpt_marker (const char *line, char *label, 
			   const char *format)
{
    const char *p = strchr(line, '#');

#if GPDEBUG > 1
    fprintf(stderr, "get_gpt_marker, p='%s'\n", p);
#endif

    if (p != NULL) {
	sscanf(p + 2, format, label);
#if GPDEBUG > 1
	fprintf(stderr, "read marker: '%s'\n", label);
#endif
	return 0;
    }

    return 1;
}

/* special graphs for which editing via GUI is not supported */

#define cant_edit(p) (p == PLOT_CORRELOGRAM || \
                      p == PLOT_LEVERAGE || \
                      p == PLOT_MULTI_IRF || \
                      p == PLOT_MULTI_SCATTER || \
                      p == PLOT_PANEL || \
                      p == PLOT_TRI_GRAPH || \
                      p == PLOT_BI_GRAPH || \
		      p == PLOT_ELLIPSE || \
		      p == PLOT_RQ_TAU)

/* graphs where we don't attempt to find data coordinates */

#define no_readback(p) (p == PLOT_CORRELOGRAM || \
                        p == PLOT_LEVERAGE || \
                        p == PLOT_MULTI_IRF || \
                        p == PLOT_MULTI_SCATTER || \
                        p == PLOT_PANEL || \
                        p == PLOT_TRI_GRAPH || \
                        p == PLOT_BI_GRAPH)

static int get_gpt_data (GPT_SPEC *spec, int datacols, int do_markers, 
			 const char *buf)
{
    char s[MAXLEN];
    char *got = NULL;
    double *x[5] = { NULL };
    char test[5][32];
    char obsfmt[12] = {0};
    int started_data_lines = 0;
    int i, j, t, imin = 0;
    int err = 0;

    spec->okobs = spec->nobs;

    gretl_push_c_numeric_locale();

    /* first handle "shaded bars" info, if present */

    for (i=0; i<spec->nbars && !err; i++) {
	double y1, y2, dx[2];

	/* start date */
	got = bufgets(s, sizeof s, buf);
	if (got == NULL ||
	    sscanf(s, "%lf %lf %lf", &dx[0], &y1, &y2) != 3) {
	    err = 1;
	    break;
	}
	/* stop date */
	got = bufgets(s, sizeof s, buf);
	if (got == NULL ||
	    sscanf(s, "%lf %lf %lf", &dx[1], &y1, &y2) != 3) {
	    err = 1;
	    break;
	}
	/* trailing 'e' */
	got = bufgets(s, sizeof s, buf);
	if (got) {
	    plotspec_set_bar_info(spec, i, dx[0], dx[1]);
	} else {
	    err = 1;
	}
    }

    if (do_markers) {
	sprintf(obsfmt, "%%%d[^\r\n]", OBSLEN - 1);
    }

    /* then get the regular plot data */

    for (i=0; i<spec->n_lines && !err; i++) {
	int ncols = gp_line_data_columns(spec, i);
	int okobs = spec->nobs;
	int offset = 1;

	if (ncols == 0) {
	    if (i == 0) {
		imin = 1;
	    }
	    continue;
	}

#if GPDEBUG
	fprintf(stderr, "reading data, line %d\n", i);
#endif

	if (!started_data_lines) {
	    offset = 0;
	    x[0] = spec->data;
	    x[1] = x[0] + spec->nobs;
	    started_data_lines = 1;
	} 

	x[2] = x[1] + spec->nobs;
	x[3] = x[2] + spec->nobs;
	x[4] = x[3] + spec->nobs;

	for (t=0; t<spec->nobs && !err; t++) {
	    int missing = 0;
	    int nf = 0;

	    got = bufgets(s, sizeof s, buf);
	    if (got == NULL) {
		err = 1;
		break;
	    }

	    nf = 0;
	    if (ncols == 5) {
		nf = sscanf(s, "%31s %31s %31s %31s %31s", test[0], test[1], test[2], 
			    test[3], test[4]);
	    } else if (ncols == 4) {
		nf = sscanf(s, "%31s %31s %31s %31s", test[0], test[1], test[2], test[3]);
	    } else if (ncols == 3) {
		nf = sscanf(s, "%31s %31s %31s", test[0], test[1], test[2]);
	    } else if (ncols == 2) {
		nf = sscanf(s, "%31s %31s", test[0], test[1]);
	    }

	    if (nf != ncols) {
		err = 1;
	    }

	    for (j=offset; j<nf && !err; j++) {
		if (test[j][0] == '?') {
		    x[j][t] = NADBL;
		    missing++;
		} else if (j == 0 && (spec->flags & GPT_TIMEFMT)) {
		    x[j][t] = gnuplot_time_from_date(test[j], spec->timefmt);
		    if (na(x[j][t])) {
			err = E_DATA;
		    }
		} else {
		    x[j][t] = atof(test[j]);
		}
	    }

	    if (missing) {
		okobs--;
	    }

	    if (i <= imin && do_markers) {
		get_gpt_marker(s, spec->markers[t], obsfmt);
		if (spec->code == PLOT_FACTORIZED && imin == 0) {
		    imin = 1;
		}
	    }
	}

	if (okobs < spec->okobs) {
	    spec->okobs = okobs;
	} 

	/* trailer line for data block */
	bufgets(s, sizeof s, buf);

	/* shift 'y' writing location */
	x[1] += (ncols - 1) * spec->nobs;
    }

#if GPDEBUG
    for (i=0; i<spec->nobs*datacols; i++) {
	fprintf(stderr, "spec->data[%d] = %g\n", i, spec->data[i]);
    }
#endif

    gretl_pop_c_numeric_locale();

    return err;
}

/* read a gnuplot source line specifying a text label */

static int parse_label_line (GPT_SPEC *spec, const char *line)
{
    const char *p, *s;
    char *text = NULL;
    double x, y;
    int nc, just = GP_JUST_LEFT;
    int err = 0;

    /* Examples:
       set label "this is a label" at 1998.26,937.557 left front
       set label 'foobar' at 1500,350 left 
    */

    /* find first double or single quote */
    p = strchr(line, '"');
    if (p == NULL) {
	p = strchr(line, '\'');
    }

    if (p == NULL) {
	/* no label text found */
	return 1;
    }

    text = gretl_quoted_string_strdup(p, &s);
    if (text == NULL) {
	return 1;
    }

    /* get the position */
    p = strstr(s, "at");
    if (p == NULL) {
	err = E_DATA;
    } else {
	p += 2;
	gretl_push_c_numeric_locale();
	nc = sscanf(p, "%lf,%lf", &x, &y);
	gretl_pop_c_numeric_locale();
	if (nc != 2) {
	    err = E_DATA;
	} 
    }

    if (!err) {
	/* justification */
	if (strstr(p, "right")) {
	    just = GP_JUST_RIGHT;
	} else if (strstr(p, "center")) {
	    just = GP_JUST_CENTER;
	}
    }
    
    if (!err) {
	err = plotspec_add_label(spec);
	if (!err) {
	    int i = spec->n_labels - 1;

	    strncat(spec->labels[i].text, text, PLOT_LABEL_TEXT_LEN);
	    spec->labels[i].pos[0] = x;
	    spec->labels[i].pos[1] = y;
	    spec->labels[i].just = just;
	}
    }

    free(text);

    return err;
}

/* read a gnuplot source line specifying an arrow */

static int parse_arrow_line (GPT_SPEC *spec, const char *line)
{
    double x0, y0, x1, y1;
    char head[12] = {0};
    int n, err = 0;

    /* Example:
       set arrow from 2000,150 to 2000,500 nohead [ lt 0 ]
    */

    gretl_push_c_numeric_locale();
    n = sscanf(line, "set arrow from %lf,%lf to %lf,%lf %s", 
	       &x0, &y0, &x1, &y1, head);
    gretl_pop_c_numeric_locale();

    if (n < 5) {
	err = E_DATA;
    } else {
	err = plotspec_add_arrow(spec);
    }

    if (!err) {
	int flags = 0;

	if (strcmp(head, "nohead")) {
	    flags |= GP_ARROW_HEAD;
	}
	if (strstr(line, "lt 0")) {
	    flags |= GP_ARROW_DOTS;
	}	
	n = spec->n_arrows - 1;
	spec->arrows[n].x0 = x0;
	spec->arrows[n].y0 = y0;
	spec->arrows[n].x1 = x1;
	spec->arrows[n].y1 = y1;
	spec->arrows[n].flags = flags;
    }

    return err;
}

static int 
read_plotspec_range (const char *obj, const char *s, GPT_SPEC *spec)
{
    double r0, r1;
    int i = 0, err = 0;

    if (!strcmp(obj, "xrange")) {
	i = 0;
    } else if (!strcmp(obj, "yrange")) {
	i = 1;
    } else if (!strcmp(obj, "y2range")) {
	i = 2;
    } else if (!strcmp(obj, "trange")) {
	i = 3;
    } else if (!strcmp(obj, "x2range")) {
	i = 4;
    } else {
	err = 1;
    }

    if (!strcmp(s, "[*:*]")) {
	r0 = r1 = NADBL;
    } else {
	gretl_push_c_numeric_locale();
	if (!err && sscanf(s, "[%lf:%lf]", &r0, &r1) != 2) {
	    err = 1;
	}
	gretl_pop_c_numeric_locale();
    }

    if (!err) {
	spec->range[i][0] = r0;
	spec->range[i][1] = r1;
    }

    return err;
}

static int read_plot_logscale (const char *s, GPT_SPEC *spec)
{
    char axis[3] = {0};
    double base = 0;
    int i, n, err = 0;

    n = sscanf(s, "%2s %lf", axis, &base);

    if (n < 1 || (n == 2 && base < 1.1)) {
	err = 1;
    } else {
	if (n == 1) {
	    base = 10.0;
	}
	if (!strcmp(axis, "x")) {
	    i = 0;
	} else if (!strcmp(axis, "y")) {
	    i = 1;
	} else if (!strcmp(axis, "y2")) {
	    i = 2;
	} else {
	    err = 1;
	}
    }
    
    if (!err) {
	spec->logbase[i] = base;
    }

    return err;
}

static int read_plot_format (const char *s, GPT_SPEC *spec,
			     int timefmt)
{
    char axis, fmt[16];
    int n, err = 0;

    n = sscanf(s, "%c \"%15[^\"]", &axis, fmt);

    if (n < 2) {
	err = 1;
    } else if (timefmt) {
	if (axis == 'x') {
	    *spec->timefmt = '\0';
	    strncat(spec->timefmt, fmt, 15);
	}
    } else if (axis == 'x') {
	*spec->xfmt = '\0';
	strncat(spec->xfmt, fmt, 15);
    } else if (axis == 'y') {
	*spec->yfmt = '\0';
	strncat(spec->yfmt, fmt, 15);
    }

    return err;
}

static int catch_value (char *targ, const char *src, int maxlen)
{
    int i, n;

    src += strspn(src, " \t\r\n");
    if (*src == '\'' || *src == '"') {
	src++;
    }

    *targ = '\0';

    if (*src != '\0') {
	strncat(targ, src, maxlen - 1);
	n = strlen(targ);

	for (i=n-1; i>=0; i--) {
	    if (isspace((unsigned char) targ[i])) {
		targ[i] = '\0';
	    } else {
		break;
	    }
	}  
	if (targ[i] == '\'' || targ[i] == '"') {
	    targ[i] = '\0';
	}
    }

    return (*targ != '\0');
}

static void fix_old_roots_plot (GPT_SPEC *spec)
{
    strings_array_free(spec->literal, spec->n_literal);
    spec->literal = NULL;
    spec->n_literal = 0;

    spec->flags |= (GPT_POLAR | GPT_XZEROAXIS | GPT_YZEROAXIS);
    spec->border = 0;
    spec->keyspec = GP_KEY_NONE;
    strcpy(spec->xtics, "none");
    strcpy(spec->ytics, "none");
}

static int parse_gp_unset_line (GPT_SPEC *spec, const char *s)
{
    char key[16] = {0};
    int err = 0;

    if (sscanf(s + 6, "%15s", key) != 1) {
	err = 1;
    } else if (!strcmp(key, "border")) {
	spec->border = 0;
    } else if (!strcmp(key, "key")) {
	spec->keyspec = GP_KEY_NONE;
    } else if (!strcmp(key, "xtics")) {
	strcpy(spec->xtics, "none");
    } else if (!strcmp(key, "ytics")) {
	strcpy(spec->ytics, "none");
    } else {
	/* unrecognized "unset" parameter */
	err = 1;
    }

    if (err) {
	errbox(_("Failed to parse gnuplot file"));
	fprintf(stderr, "parse_gp_unset_line: bad line '%s'\n", s);
    }

    return err;
}

static void read_xtics_setting (GPT_SPEC *spec, 
				const char *key,
				const char *val)
{
    if (!strcmp(key, "x2tics")) {
	spec->x2ticstr = gretl_strdup(val);
    } else if (*val == '(') {
	spec->xticstr = gretl_strdup(val);
    } else {
	*spec->xtics = '\0';
	strncat(spec->xtics, val, sizeof(spec->xtics) - 1);
    }
}

static int parse_gp_set_line (GPT_SPEC *spec, const char *s, 
			      linestyle *styles)
{
    char key[16] = {0};
    char val[MAXLEN] = {0};

    if (!strncmp(s, "set style line", 14)) {
	/* e.g. set style line 1 lc rgb "#ff0000 lt 6" */
	int n, idx = 0, lt = LT_AUTO;
	char rgb[8] = {0};

	n = sscanf(s + 14, " %d lc rgb \"%7s\" lt %d", &idx, rgb, &lt);
	if (n >= 2 && idx > 0 && idx <= MAX_STYLES) {
	    strcpy(styles[idx-1].rgb, rgb);
	    styles[idx-1].type = (n == 3)? lt : LT_AUTO;
	}
	return 0;
    }

    if (sscanf(s + 4, "%11s", key) != 1) {
	errbox(_("Failed to parse gnuplot file"));
	fprintf(stderr, "parse_gp_set_line: bad line '%s'\n", s);
	return 1;
    }

#if GPDEBUG
    fprintf(stderr, "parse_gp_set_line: key = '%s'\n", key);
#endif

    if (!strcmp(key, "term") || 
	!strcmp(key, "output") ||
	!strcmp(key, "encoding") ||
	!strcmp(key, "datafile")) {
	/* we ignore these */
	return 0;
    }

    /* first, settings that don't require a parameter */

    if (!strcmp(key, "y2tics")) {
	spec->flags |= GPT_Y2AXIS;
	return 0;
    } else if (!strcmp(key, "parametric")) {
	spec->flags |= GPT_PARAMETRIC;
	return 0;
    } else if (!strcmp(key, "polar")) {
	spec->flags |= GPT_POLAR;
	return 0;
    } else if (!strcmp(key, "xzeroaxis")) {
	spec->flags |= GPT_XZEROAXIS;
	return 0;
    } else if (!strcmp(key, "yzeroaxis")) {
	spec->flags |= GPT_YZEROAXIS;
	return 0;
    } else if (!strcmp(key, "noxtics")) {
	strcpy(spec->xtics, "none");
	return 0;
    } else if (!strcmp(key, "noytics")) {
	strcpy(spec->ytics, "none");
	return 0;
    } else if (!strcmp(key, "nokey")) {
	spec->keyspec = GP_KEY_NONE;
	return 0;
    } else if (!strcmp(key, "label")) {
	parse_label_line(spec, s);
	return 0;
    } else if (!strcmp(key, "arrow")) {
	parse_arrow_line(spec, s);
	return 0;
    }

    /* grid lines: parameter is optional */
    if (!strcmp(key, "grid")) {
	if (catch_value(val, s + 4 + strlen(key), MAXLEN)) {
	    if (!strcmp(val, "ytics")) {
		spec->flags |= GPT_GRID_Y;
	    } else if (!strcmp(val, "xtics")) {
		spec->flags |= GPT_GRID_X;
	    }
	} else {
	    spec->flags |= GPT_GRID_Y;
	    spec->flags |= GPT_GRID_X;
	}
    }

    /* now catch value for settings that need a parameter */

    if (!catch_value(val, s + 4 + strlen(key), MAXLEN)) {
	return 0;
    }

#if GPDEBUG
    fprintf(stderr, " value = '%s'\n", val);
#endif

    if (strstr(key, "range")) {
	if (read_plotspec_range(key, val, spec)) {
	    errbox(_("Failed to parse gnuplot file"));
	    fprintf(stderr, "parse_gp_set_line: bad line '%s'\n", s);
	    return 1;
	}
    } else if (!strcmp(key, "logscale")) {
	if (read_plot_logscale(val, spec)) {
	    errbox(_("Failed to parse gnuplot file"));
	    fprintf(stderr, "parse_gp_set_line: bad line '%s'\n", s);
	    return 1;
	}
    } else if (!strcmp(key, "format")) {
	if (read_plot_format(val, spec, 0)) {
	    errbox(_("Failed to parse gnuplot file"));
	    fprintf(stderr, "parse_gp_set_line: bad line '%s'\n", s);
	    return 1;
	}
    } else if (!strcmp(key, "timefmt")) {
	if (read_plot_format(val, spec, 1)) {
	    errbox(_("Failed to parse gnuplot file"));
	    fprintf(stderr, "parse_gp_set_line: bad line '%s'\n", s);
	    return 1;
	}	
    } else if (!strcmp(key, "title")) {
	strcpy(spec->titles[0], val);
    } else if (!strcmp(key, "xlabel")) {
	strcpy(spec->titles[1], val);
	*spec->xvarname = '\0';
	strncat(spec->xvarname, val, MAXDISP-1);
    } else if (!strcmp(key, "ylabel")) {
	strcpy(spec->titles[2], val);
	*spec->yvarname = '\0';
	strncat(spec->yvarname, val, MAXDISP-1);
    } else if (!strcmp(key, "y2label")) {
	strcpy(spec->titles[3], val);
    } else if (!strcmp(key, "x2label")) {
	strcpy(spec->titles[4], val);
    } else if (!strcmp(key, "key")) {
	spec->keyspec = gp_keypos_from_name(val);
    } else if (!strcmp(key, "xtics") || !strcmp(key, "x2tics")) { 
	read_xtics_setting(spec, key, val);
    } else if (!strcmp(key, "mxtics")) {
	*spec->mxtics = '\0';
	strncat(spec->mxtics, val, sizeof(spec->mxtics) - 1);
    } else if (!strcmp(key, "ytics")) {
	*spec->ytics = '\0';
	strncat(spec->ytics, val, sizeof(spec->ytics) - 1);
    } else if (!strcmp(key, "border")) {
	spec->border = atoi(val);
    } else if (!strcmp(key, "bmargin")) {
	spec->bmargin = atoi(val);
    } else if (!strcmp(key, "boxwidth")) {
	spec->boxwidth = (float) atof(val);
	if (strstr(s, "absolute")) {
	    spec->boxwidth = -spec->boxwidth;
	}
    } else if (!strcmp(key, "samples")) {
	spec->samples = atoi(val);
    } else if (!strcmp(key, "xdata")) {
	if (!strncmp(val, "time", 4)) {
	    spec->flags |= GPT_TIMEFMT;
	}
    }

    return 0;
}

/* allocate markers for identifying particular data points */

static int plotspec_allocate_markers (GPT_SPEC *spec)
{
    spec->markers = strings_array_new_with_length(spec->nobs,
						  OBSLEN);
    if (spec->markers == NULL) {
	spec->n_markers = 0;
	return E_ALLOC;
    } else {
	spec->n_markers = spec->nobs;
	return 0;
    }
}

/* Determine the number of data points in a plot.  While we're at it,
   determine the type of plot, check whether there are any
   data-point markers along with the data, and see if there are
   are definitions of time-series vertical bars.
*/

static void get_plot_nobs (GPT_SPEC *spec, const char *buf, int *do_markers)
{
    int n = 0, started = -1;
    int startmin = 1;
    char line[MAXLEN], test[9];
    char *p;

    spec->code = PLOT_REGULAR;
    *do_markers = 0;

    while (bufgets(line, MAXLEN - 1, buf)) {

	if (*line == '#' && spec->code == PLOT_REGULAR) {
	    tailstrip(line);
	    spec->code = plot_type_from_string(line);
	}

	if (sscanf(line, "# n_bars = %d", &spec->nbars) == 1) {
	    startmin += spec->nbars;
	    continue;
	}

	if (!strncmp(line, "plot", 4)) {
	    started = 0;
	}

	if (started == 0 && strchr(line, '\\') == NULL) {
	    started = 1;
	    continue;
	}

	if (started > 0 && started < startmin) {
	    if (*line == 'e') {
		started++;
		continue;
	    }
	}

	if (started == startmin) {
	    if (*do_markers == 0 && (p = strchr(line, '#')) != NULL) {
		if (sscanf(p + 1, "%8s", test) == 1) {
		    *do_markers = 1;
		}
	    }
	    if (*line == 'e' || !strncmp(line, "set ", 4)) {
		/* end of data, or onto "set print" for bounds */
		break;
	    }
	    n++;
	}
    }

    if (spec->code == PLOT_BOXPLOTS) {
	/* we stopped reading at the first "e" line, but the
	   file may contain auxiliary outlier data, in which
	   case we need to read this into an auxiliary matrix
	*/
	double x, y;
	int i, r = 0, c = 0;
	int err = 0;

	while (bufgets(line, MAXLEN - 1, buf)) {
	    if (!strncmp(line, "# auxdata", 9)) {
		if (sscanf(line + 9, "%d %d", &r, &c) == 2 &&
		    r > 0 && c > 0) {
		    spec->auxdata = gretl_matrix_alloc(r, c);
		    if (spec->auxdata != NULL) {
			gretl_push_c_numeric_locale();
			for (i=0; i<r && !err; i++) {
			    if (bufgets(line, MAXLEN - 1, buf) == NULL) {
				err = E_DATA;
			    } else if (sscanf(line, "%lf %lf", &x, &y) != 2) {
				err = E_DATA;
			    } else {
				gretl_matrix_set(spec->auxdata, i, 0, x);
				gretl_matrix_set(spec->auxdata, i, 1, y);
			    }
			}
			gretl_pop_c_numeric_locale();
			if (err) {
			    gretl_matrix_free(spec->auxdata);
			    spec->auxdata = NULL;
			}
		    }
		}
		break;
	    }
	}
    }

    spec->nobs = n;
}

static int grab_fit_coeffs (GPT_SPEC *spec, const char *s)
{
    gretl_matrix *b = NULL;
    int k, f = spec->fit;
    int n, err = 0;

    k = (f == PLOT_FIT_CUBIC)? 4 : (f == PLOT_FIT_QUADRATIC)? 3 : 2;

    b = gretl_column_vector_alloc(k);
    if (b == NULL) {
	return E_ALLOC;
    }

    if (!strncmp(s, "exp(", 4)) {
	s += 4;
    }

    gretl_push_c_numeric_locale();
    if (k == 2) {
	n = sscanf(s, "%lf + %lf", &b->val[0], &b->val[1]);
    } else if (k == 3) {
	n = sscanf(s, "%lf + %lf*x + %lf", &b->val[0],
		   &b->val[1], &b->val[2]);
    } else if (k == 4) {
	n = sscanf(s, "%lf + %lf*x + %lf*x**2 + %lf", &b->val[0], 
		   &b->val[1], &b->val[2], &b->val[3]);
    }
    gretl_pop_c_numeric_locale();

    if (n != k) {
	err = E_DATA;
	gretl_matrix_free(b);
	spec->flags &= ~GPT_AUTO_FIT;
    } else if (f == PLOT_FIT_OLS) {
	spec->b_ols = b;
    } else if (f == PLOT_FIT_INVERSE) {
	spec->b_inv = b;
    } else if (f == PLOT_FIT_LOGLIN) {
	spec->b_log = b;
    } else if (f == PLOT_FIT_QUADRATIC) {
	spec->b_quad = b;
    } else if (f == PLOT_FIT_CUBIC) {
	spec->b_cub = b;
    }

    return err;
}

/* scan the stuff after "title '" or 'title "' */

static void grab_line_title (char *targ, const char *src)
{
    char *fmt;

    if (*src == '\'') {
	fmt = "%79[^']'";
    } else {
	fmt = "%79[^\"]\"";
    }

    sscanf(src + 1, fmt, targ);
}

/* parse the "using..." portion of plot specification for a
   given plot line: full form is like:
  
     using XX axes XX title XX w XX lt XX lw XX
*/

static int parse_gp_line_line (const char *s, GPT_SPEC *spec,
			       int auto_linewidth)
{
    GPT_LINE *line;
    char tmp[16];
    const char *p;
    int i, err;

    err = plotspec_add_line(spec);
    if (err) {
	return err;
    }

    i = spec->n_lines - 1;
    line = &spec->lines[i];

    if (!strncmp(s, "plot ", 5)) {
	s += 5;
    }

    s += strspn(s, " ");

    if ((p = strstr(s, " using "))) {
	/* data column spec */
	p += 7;
	if (strstr(p, "1:3:2:5:4")) {
	    line->ncols = 5;
	} else if (strstr(p, "1:2:2:2:2")) {
	    line->ncols = 2;
	    line->flags |= GP_LINE_BOXDATA;
	} else if (strstr(p, "1:2:3:4")) {
	    line->ncols = 4;
	} else if (strstr(p, "1:2:3")) {
	    line->ncols = 3;
	} else if ((p = strstr(s, "($2*"))) {
	    sscanf(p + 4, "%15[^)]", tmp);
	    line->scale = dot_atof(tmp);
	    line->ncols = 2;
	} else {
	    line->ncols = 2;
	}
    } else if (*s == '\'' || *s == '"') {
	/* name of data file, without 'using' */
	if (*(s+1) != '-') {
	    fprintf(stderr, "plotting datafile, not supported\n");
	} else {
	    line->ncols = 2;
	}
    } else {
	/* absence of "using" should mean that the line plots 
	   a formula, not a set of data columns 
	*/
	line->scale = NADBL;
	/* get the formula: it runs up to "title" or "notitle" */
	p = strstr(s, " title");
	if (p == NULL) {
	    p = strstr(s, " notitle");
	}
	if (p != NULL) {
	    strncat(line->formula, s, p - s);
	    if (i == 1 && spec->flags & GPT_AUTO_FIT) {
		grab_fit_coeffs(spec, line->formula);
	    }
	} else {
	    /* title must be implicit? */
	    char fmt[8];

	    sprintf(fmt, "%%%ds", GP_MAXFORMULA - 1);
	    sscanf(s, fmt, line->formula);
	    strcpy(fmt, "%79s");
	    sscanf(s, fmt, line->title);
	}
    }

    if (strstr(s, "axes x1y2")) {
	line->yaxis = 2;
    } 

    if ((p = strstr(s, " title "))) {
	grab_line_title(line->title, p + 7);
    }

    if ((p = strstr(s, " w "))) {
	sscanf(p + 3, "%15[^, ]", tmp);
	line->style = gp_style_index_from_name(tmp);
    } 

    if ((p = strstr(s, " lt "))) {
	sscanf(p + 4, "%d", &line->type);
    }

    if ((p = strstr(s, " ps "))) {
	sscanf(p + 4, "%lf", &line->pscale);
    }

    if ((p = strstr(s, " pt "))) {
	sscanf(p + 4, "%d", &line->ptype);
    }

    if (!auto_linewidth && (p = strstr(s, " lw "))) {
	sscanf(p + 4, "%d", &line->width);
    } 

    if ((p = strstr(s, " whiskerbars "))) {
	double ww;

	sscanf(p + 13, "%lf", &ww);
	line->whiskwidth = (float) ww;
    } 

    if (line->ncols == 0 && line->formula[0] == '\0') {
	/* got neither data column spec nor formula */
	err = 1;
    }

#if GPDEBUG
    fprintf(stderr, "parse_gp_line_line, returning %d\n", err);
#endif

    return err;
}

/* We got a special comment supposedly indicating the name and ID
   number of the independent variable, X, in an OLS fitted line.  Here
   we check this info for validity: @vname should be the name of a
   bona fide series variable, and its ID number should match the given
   @v.  We return 1 if this checks out OK, otherwise 0.
*/

static int plot_ols_var_ok (const char *vname)
{
    int v = current_series_index(dataset, vname);

    if (v >= 0 && !strcmp(dataset->varname[v], vname)) {
	return 1;
    }

    return 0;
}

static void maybe_set_add_fit_ok (GPT_SPEC *spec)
{
    if (spec->n_lines == 2 && spec->fit != PLOT_FIT_NONE) {
	; /* already OK */
    } else if (spec->data != NULL &&
	       spec->code == PLOT_REGULAR &&
	       spec->n_lines == 1 &&
	       spec->lines[0].ncols == 2) {
	if (spec->flags & GPT_TS) {
	    spec->fit = dataset_is_time_series(dataset) ?
		PLOT_FIT_NONE : PLOT_FIT_NA;
	} else {
	    spec->fit = PLOT_FIT_NONE;
	}
    } else {
	spec->fit = PLOT_FIT_NA;
    }
}

static int 
plot_get_data_and_markers (GPT_SPEC *spec, const char *buf, 
			   int datacols, int do_markers)
{
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "allocating: nobs=%d, datacols=%d, size=%d\n", 
	    spec->nobs, datacols, spec->nobs * datacols * sizeof *spec->data);
#endif  

    /* allocate for the plot data... */
    spec->data = mymalloc(spec->nobs * datacols * sizeof *spec->data);
    if (spec->data == NULL) {
	err = E_ALLOC;
    }

    /* and markers if any */
    if (!err && do_markers) {
	err = plotspec_allocate_markers(spec);
    }

    /* and time-series bars, if any */
    if (!err && spec->nbars > 0) {
	err = plotspec_allocate_bars(spec);
    }

    /* read the data (and perhaps markers) from the plot file */
    if (!err) {
	err = get_gpt_data(spec, datacols, do_markers, buf);
    } else {
	free(spec->data);
	spec->data = NULL;
    }

#if GPDEBUG
    fprintf(stderr, "plot_get_data_and_markers:\n"
	    " spec->data = %p, spec->markers = %p, spec->n_markers = %d, err = %d\n",
	    spec->data, (void *) spec->markers, spec->n_markers, err);
#endif

    return err;
}

/* Parse special comment to get the 0-based index numbers of
   any user-defined lines that have been added to the plot.
   The comment looks like:

   # %d user-defined lines: %d %d ...
*/

static int *get_user_lines_list (const char *s)
{
    int *list = NULL;
    char id[6];
    int i, n = 0;

    sscanf(s, "# %d", &n);

    if (n > 0) {
	s = strchr(s, ':');
	if (s != NULL) {
	    list = gretl_list_new(n);
	    if (list != NULL) {
		s++;
		for (i=0; i < n && *s != '\0'; i++) {
		    s++;
		    if (sscanf(s, "%5[^ ]", id)) {
			list[i+1] = atoi(id);
			s += strlen(id);
		    }
		}
	    }
	}
    }

    return list;
}

static FitType recognize_fit_string (const char *s)
{
    if (strstr(s, "OLS")) {
	return PLOT_FIT_OLS;
    } else if (strstr(s, "quadratic")) {
	return PLOT_FIT_QUADRATIC;
    } else if (strstr(s, "cubic")) {
	return PLOT_FIT_CUBIC;
    } else if (strstr(s, "inverse")) {
	return PLOT_FIT_INVERSE;
    } else if (strstr(s, "loess")) {
	return PLOT_FIT_LOESS;
    } else if (strstr(s, "semilog")) {
	return PLOT_FIT_LOGLIN;
    } else {
	return PLOT_FIT_NONE;
    }
}

static void check_for_plot_size (GPT_SPEC *spec, gchar *buf)
{
    char line[128];
    int i = 0;

    while (bufgets(line, sizeof line, buf) && i < 6) {
	if (!strncmp(line, "# multiple ", 11)) {
	    if (strstr(line, "extra-large")) {
		spec->flags |= GPT_XXL;
		break;
	    } else if (strstr(line, "large")) {
		spec->flags |= GPT_XL;
		break;
	    }		
	}
	i++;
    }
}

static void linestyle_init (linestyle *ls)
{
    ls->rgb[0] = '\0';
    ls->type = LT_AUTO;
}

#define plot_needs_obs(c) (c != PLOT_ELLIPSE && \
                           c != PLOT_PROB_DIST && \
                           c != PLOT_CURVE && \
			   c != PLOT_USER)

/* Read plotspec struct from gnuplot command file.  This is not a
   general parser for gnuplot files; it is designed specifically for
   files auto-generated by gretl and even then there are some sorts
   of plot that are not fully handled.
*/

static int read_plotspec_from_file (GPT_SPEC *spec, int *plot_pd)
{
    linestyle styles[MAX_STYLES];
    int do_markers = 0;
    int datacols = 0;
    int auto_linewidth = 0;
    int reglist[4] = {0};
    int *uservec = NULL;
    char gpline[MAXLEN];
    gchar *buf = NULL;
    char *got = NULL;
    int i, done;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "read_plotspec_from_file: spec=%p\n", 
	    (void *) spec);
#endif

    /* check: are we already done? */
    if (PLOTSPEC_DETAILS_IN_MEMORY(spec)) {
#if GPDEBUG
	fprintf(stderr, " info already in memory, returning 0\n");
#endif
	return 0;
    }

    /* get the content of the plot file */
    err = gretl_file_get_contents(spec->fname, &buf, NULL);
    if (err) {
	return err;
    }

    bufgets_init(buf);

    /* get the number of data-points, plot type, and check for 
       observation markers and time-series bars */
    get_plot_nobs(spec, buf, &do_markers);

    if (spec->nobs == 0 && plot_needs_obs(spec->code)) {
	/* failed reading plot data */
#if GPDEBUG
	fprintf(stderr, " got spec->nobs = 0\n");
#endif
	err = 1;
	goto bailout;
    }

    if (spec->nobs > MAX_MARKERS && do_markers) {
	do_markers = 0;
    }

    if (cant_edit(spec->code)) {
	fprintf(stderr, "read_plotspec_from_file: plot is not editable\n");
	if (maybe_big_multiplot(spec->code)) {
	    buf_rewind(buf);
	    check_for_plot_size(spec, buf);
	}
	goto bailout;
    }

    for (i=0; i<MAX_STYLES; i++) {
	linestyle_init(&styles[i]);
    }

    buf_rewind(buf);

    /* get the preamble and "set" lines */

    while ((got = bufgets(gpline, sizeof gpline, buf))) {
	char vname[VNAMELEN];
	int v;

#if GPDEBUG
	tailstrip(gpline);
	fprintf(stderr, "gpline: '%s'\n", gpline);
#endif

	if (!strncmp(gpline, "# timeseries", 12)) {
	    if (sscanf(gpline, "# timeseries %d", &spec->pd)) {
		*plot_pd = spec->pd;
	    }
	    spec->flags |= GPT_TS;
	    if (strstr(gpline, "letterbox")) {
		spec->flags |= GPT_LETTERBOX;
	    }
	    continue;
	} else if (!strncmp(gpline, "# multiple timeseries", 21)) {
	    if (sscanf(gpline, "# multiple timeseries %d", &spec->pd)) {
		*plot_pd = spec->pd;
	    }
	    spec->flags |= GPT_TS;
	    continue;
	} else if (!strncmp(gpline, "# scale = ", 10)) {
	    double x = dot_atof(gpline + 10);
	    
	    if (x >= 0.5 && x <= 2.0) {
		spec->scale = x;
	    }
	    continue;
	} else if (!strncmp(gpline, "# auto linewidth", 16)) {
	    auto_linewidth = 1;
	    continue;
	} else if (!strncmp(gpline, "# boxplots", 10)) {
	    continue;
	}

	if (!strncmp(gpline, "# X = ", 6) || !strncmp(gpline, "# Y = ", 6)) {
	    char fmt[16];

	    sprintf(fmt, "'%%%d[^\\']' (%%d)", VNAMELEN - 1);

	    if (sscanf(gpline + 6, fmt, vname, &v) == 2) {
		if (gpline[2] == 'X') {
		    if (plot_ols_var_ok(vname)) {
			reglist[2] = v;
		    }
		} else { 
		    /* 'Y' */
		    if (reglist[2] > 0 && plot_ols_var_ok(vname)) {
			reglist[0] = 3;
			reglist[1] = v;
		    }
		}
	    }
	    continue;
	}
	
	if (sscanf(gpline, "# literal lines = %d", &spec->n_literal)) {
	    spec->literal = strings_array_new(spec->n_literal);
	    if (spec->literal == NULL) {
		err = E_ALLOC;
		goto bailout;
	    }
	    for (i=0; i<spec->n_literal; i++) {
		if (!bufgets(gpline, MAXLEN - 1, buf)) {
		    errbox(_("Plot file is corrupted"));
		} else {
		    top_n_tail(gpline, 0, NULL);
		    spec->literal[i] = g_strdup(gpline);
		}
	    }
	    continue;
	}

	if (strstr(gpline, "automatic fit")) {
	    spec->flags |= GPT_AUTO_FIT;
	    spec->fit = recognize_fit_string(gpline);
	    continue;
	}

	if (strstr(gpline, "user-defined")) {
	    uservec = get_user_lines_list(gpline);
	    continue;
	}

	if (strstr(gpline, "printing data labels")) {
	    spec->flags |= GPT_PRINT_MARKERS;
	    continue;
	}

	if (!strncmp(gpline, "# ", 2)) {
	    /* ignore unknown comment lines */
	    continue;
	}

	if (strncmp(gpline, "set ", 4) && strncmp(gpline, "unset ", 6)) {
	    /* done reading "set" lines */
	    break;
	}

	if (!strncmp(gpline, "set ", 4)) {
	    err = parse_gp_set_line(spec, gpline, styles);
	} else if (!strncmp(gpline, "unset ", 6)) {
	    err = parse_gp_unset_line(spec, gpline);
	} else {
	    /* done reading "set" lines */
	    break;
	}	    

	if (err) {
	    goto bailout;
	}
    }

    if (got == NULL) {
	err = 1;
	goto bailout;
    }

    if ((spec->flags & GPT_TIMEFMT) && *spec->timefmt == '\0') {
	fprintf(stderr, "got 'xdata time' but no timefmt\n");
	err = E_DATA;
	goto bailout;
    }

    for (i=0; i<4; i++) {
	if (spec->titles[i][0] != '\0') {
	    gretl_delchar('"', spec->titles[i]);
	}
    }

    /* then get the "plot" lines */
    if (strncmp(gpline, "plot ", 5) ||
	(strlen(gpline) < 10 && bufgets(gpline, MAXLEN - 1, buf) == NULL)) {	
	errbox(_("Failed to parse gnuplot file"));
	fprintf(stderr, "bad gnuplot line: '%s'\n", gpline);
	err = 1;
	goto bailout;
    }

    /* if there are any time-series bars, skip past them */
    for (i=0; i<spec->nbars; i++) {
	if ((got = bufgets(gpline, MAXLEN - 1, buf)) == NULL) {
	    err = E_DATA;
	}
    }

    done = 0;

    while (!err) {
	top_n_tail(gpline, 0, NULL);
	if (!chop_comma(gpline)) {
	    /* line did not end with comma -> no continuation of
	       the plot command */
	    done = 1;
	} 

#if GPDEBUG
	fprintf(stderr, "calling parse_gp_line_line\n");
#endif

	err = parse_gp_line_line(gpline, spec, auto_linewidth);

	if (err || done || (got = bufgets(gpline, MAXLEN - 1, buf)) == NULL) {
	    break;
	}
    }

    if (err || got == NULL) {
	err = 1;
	goto bailout;
    }

    /* determine total number of required data columns, etc.,
       and transcribe styles info into lines for use in the
       GUI editor */

    for (i=0; i<spec->n_lines; i++) {
	if (uservec != NULL && in_gretl_list(uservec, i)) {
	    spec->lines[i].flags |= GP_LINE_USER;
	}
	if (i < MAX_STYLES) {
	    strcpy(spec->lines[i].rgb, styles[i].rgb);
	}
	if (spec->auxdata != NULL && i == spec->n_lines - 1) {
	    /* the last "line" doesn't use the regular
	       data mechanism */
	    spec->lines[i].flags = GP_LINE_AUXDATA;
	    continue;
	}
 	if (spec->lines[i].ncols == 0) {
	    continue;
	}
	if (datacols == 0) {
	    datacols = spec->lines[i].ncols;
	} else {
	    datacols += spec->lines[i].ncols - 1;
	}
    }

    err = plot_get_data_and_markers(spec, buf, datacols, do_markers);

    if (!err && reglist != NULL && reglist[0] > 0) {
	spec->reglist = gretl_list_copy(reglist);
    }

    if (!err && spec->fit == PLOT_FIT_NONE) {
	maybe_set_add_fit_ok(spec);
    }

    if (spec->code == PLOT_ROOTS && spec->n_literal == 8) {
	/* backward compatibility */
	fprintf(stderr, "old roots plot: fixing\n");
	fix_old_roots_plot(spec);
    }   

 bailout:

    bufgets_finalize(buf);
    g_free(buf);
    free(uservec);

    return err;
}

#define has_log_axis(s) (s->logbase[0] != 0 || s->logbase[1] != 0)

static int get_data_xy (png_plot *plot, int x, int y, 
			double *data_x, double *data_y)
{
    double xmin, xmax;
    double ymin, ymax;
    double dx = NADBL;
    double dy = NADBL;
    int ok = 1;

    if (plot_is_zoomed(plot)) {
	xmin = plot->zoom_xmin;
	xmax = plot->zoom_xmax;
	ymin = plot->zoom_ymin;
	ymax = plot->zoom_ymax;
    } else {
	xmin = plot->xmin;
	xmax = plot->xmax;
	ymin = plot->ymin;
	ymax = plot->ymax;
    }

#if POINTS_DEBUG
    if (plot_doing_position(plot)) {
	fprintf(stderr, "get_data_xy:\n"
		" plot->xmin=%g, plot->xmax=%g, plot->ymin=%g, plot->ymax=%g\n",
		plot->xmin, plot->xmax, plot->ymin, plot->ymax);
    }
#endif

    if (xmin == 0.0 && xmax == 0.0) { 
	fprintf(stderr, "get_data_xy: unknown x range\n");
    } else {
	dx = xmin + ((double) x - plot->pixel_xmin) / 
	    (plot->pixel_xmax - plot->pixel_xmin) * (xmax - xmin);
    }

    if (!na(dx)) {
	if (ymin == 0.0 && ymax == 0.0) { 
	    fprintf(stderr, "get_data_xy: unknown y range\n");
	} else {
	    dy = ymax - ((double) y - plot->pixel_ymin) / 
		(plot->pixel_ymax - plot->pixel_ymin) * (ymax - ymin);
	}
    }

    if (na(dx) || na(dx)) {
	ok = 0;
    } else if (has_log_axis(plot->spec)) {
	double base, dprop, lr;

	base = plot->spec->logbase[0];
	if (base != 0) {
	    if (xmin > 0) {
		dprop = (dx - xmin) / (xmax - xmin);
		lr = log(xmax / xmin) / log(base);
		dx = pow(base, dprop * lr);
	    } else {
		dx = NADBL;
		ok = 0;
	    }
	}
	base = plot->spec->logbase[1];
	if (base != 0) {
	    if (ymin > 0) {
		dprop = (dy - ymin) / (ymax - ymin);
		lr = log(ymax / ymin) / log(base);
		dy = pow(base, dprop * lr);
	    } else {
		dy = NADBL;
		ok = 0;
	    }
	}
    } else if (plot_is_polar(plot)) {
	double px = atan2(dy, dx);
	double py = sqrt(dx * dx + dy * dy);

	dx = px;
	dy = py;
    }

    *data_x = dx;
    *data_y = dy;

    return ok;
}

static void x_to_date (double x, int pd, char *str)
{
    int yr = (int) x;
    double t, frac = 1.0 / pd;
    int subper = (int) ((x - yr + frac) * pd);
    static int decpoint;

    if (decpoint == 0) {
	decpoint = get_local_decpoint();
    }

    t = yr + subper / ((pd < 10)? 10.0 : 100.0);
    sprintf(str, "%.*f", (pd < 10)? 1 : 2, t);
    gretl_charsub(str, decpoint, ':');
}

static void redraw_plot_rectangle (png_plot *plot, GdkRectangle *r)
{
    if (r != NULL) {
	gdk_window_invalidate_rect(plot->window, r, FALSE);
    } else {
	GdkRectangle rt = {
	    0, 0, plot->pixel_width, plot->pixel_height
	};

	gdk_window_invalidate_rect(plot->window, &rt, FALSE);
    }  
}

static void make_cairo_rectangle (png_plot *plot,
				  GdkRectangle *r)
{
    cairo_move_to(plot->cr, r->x, r->y);
    cairo_line_to(plot->cr, r->x, r->y + r->height);
    cairo_line_to(plot->cr, r->x + r->width, r->y + r->height);
    cairo_line_to(plot->cr, r->x + r->width, r->y);
    cairo_close_path(plot->cr);
    cairo_stroke(plot->cr);
}

# if GTK_MAJOR_VERSION == 2

static void copy_state_to_pixbuf (png_plot *plot)
{
    if (plot->savebuf != NULL) {
	g_object_unref(plot->savebuf);
    } 

    plot->savebuf = 
	gdk_pixbuf_get_from_drawable(NULL,
				     plot->pixmap,
				     NULL,
				     0, 0, 0, 0,
				     plot->pixel_width,
				     plot->pixel_height);
}

# endif

/* Note that the screen coordinates as of the last mouse
   button press are recorded in plot->screen_x0 and
   plot->screen_y0. These represent the constant corner of
   the box that should be drawn.
*/

static void draw_selection_box (png_plot *plot, int x, int y)
{
    GdkRectangle r;

    r.x = MIN(plot->screen_x0, x);
    r.y = MIN(plot->screen_y0, y);
    r.width = abs(x - plot->screen_x0);
    r.height = abs(y - plot->screen_y0);

#if GTK_MAJOR_VERSION >= 3
    plot->cr = gdk_cairo_create(plot->window);
    cairo_set_source_rgba(plot->cr, 0.3, 0.3, 0.3, 0.3);
    make_cairo_rectangle(plot, &r);
    redraw_plot_rectangle(plot, &r);
    cairo_destroy(plot->cr);
#else
    plot->cr = gdk_cairo_create(plot->pixmap);
    if (plot->savebuf != NULL) {
	/* restore state prior to zoom start */
	gdk_cairo_set_source_pixbuf(plot->cr, plot->savebuf, 0, 0);
	cairo_paint(plot->cr);
    } else {
	copy_state_to_pixbuf(plot);
    }
    cairo_set_source_rgba(plot->cr, 0.3, 0.3, 0.3, 0.3);
    make_cairo_rectangle(plot, &r);
    redraw_plot_rectangle(plot, NULL);
    cairo_destroy(plot->cr);
#endif
}

static int make_alt_label (gchar *alt, const gchar *label)
{
    double x, y;
    int err = 0;

    gretl_push_c_numeric_locale();

    if (sscanf(label, "%lf,%lf", &x, &y) != 2) {
	err = 1;
    }

    gretl_pop_c_numeric_locale();

    if (!err) {
	if (get_local_decpoint() != '.') {
	    sprintf(alt, "%.2f %.2f", x, y);
	} else {
	    sprintf(alt, "%.2f,%.2f", x, y);
	}
    }

    return err;
}

#if GTK_MAJOR_VERSION >= 3

/* given a GdkPixbuf read from file, create a corresponding cairo
   surface, attached to @plot as plot->cs. This will be used as
   the "backing store" for re-draws of the plot.
*/

static int copy_pixbuf_to_surface (png_plot *plot,
				   GdkPixbuf *pixbuf)
{
    cairo_format_t  format;
    guchar *p_data, *s_data;
    int nc, ps, ss;
    int width, height;
    int i, j;

    width = gdk_pixbuf_get_width(pixbuf);
    height = gdk_pixbuf_get_height(pixbuf);
    nc = gdk_pixbuf_get_n_channels(pixbuf);
    ps = gdk_pixbuf_get_rowstride(pixbuf);
    format = (nc == 3)? CAIRO_FORMAT_RGB24 : CAIRO_FORMAT_ARGB32;

    if (plot->cs != NULL) {
	cairo_surface_destroy(plot->cs);
    }

    plot->cs = cairo_image_surface_create(format, width, height);

    if (plot->cs == NULL || 
	cairo_surface_status(plot->cs) != CAIRO_STATUS_SUCCESS) {
	fprintf(stderr, "copy_pixbuf_to_surface: failed\n");
	return 1;
    }

    ss = cairo_image_surface_get_stride(plot->cs);
    p_data = gdk_pixbuf_get_pixels(pixbuf);
    s_data = cairo_image_surface_get_data(plot->cs);

    for (j=0; j<height; j++) {
        guchar *p_iter = p_data + j * ps,
               *s_iter = s_data + j * ss;

        for (i=0; i<width; i++) {
#if G_BYTE_ORDER == G_LITTLE_ENDIAN
            /* BGR(A) -> RGB(A) */
            s_iter[0] = p_iter[2];
            s_iter[1] = p_iter[1];
            s_iter[2] = p_iter[0];
            if (nc == 4) {
                s_iter[3] = p_iter[3];
	    }
#else
            /* (A)RGB -> RGB(A) */
            if (nc == 4) {
                s_iter[0] = p_iter[3];
	    }
            s_iter[1] = p_iter[0];
            s_iter[2] = p_iter[1];
            s_iter[3] = p_iter[2];
#endif
            p_iter += nc;
            s_iter += 4;
        }
    }

    return 0;
}

#endif

/* implements "brushing-in" of data-point labels with the mouse */

static void
write_label_to_plot (png_plot *plot, int i, gint x, gint y)
{
    GdkRectangle r = {x, y, 0, 0};
    const gchar *label = plot->spec->markers[i];
    PangoContext *context;
    PangoLayout *pl;

    if (plot_is_roots(plot)) {
	gchar alt_label[12];
	
	if (make_alt_label(alt_label, label)) {
	    return;
	}

	label = alt_label;
    }

    context = gtk_widget_get_pango_context(plot->shell);
    pl = pango_layout_new(context);
    pango_layout_set_text(pl, label, -1);

    /* add the label, then show the modified image */

#if GTK_MAJOR_VERSION >= 3
    plot->cr = cairo_create(plot->cs);
#else
    plot->cr = gdk_cairo_create(plot->pixmap);
#endif
    cairo_move_to(plot->cr, x, y);
    pango_cairo_show_layout(plot->cr, pl);
    cairo_fill(plot->cr);
    pango_layout_get_pixel_size(pl, &r.width, &r.height);
    redraw_plot_rectangle(plot, &r);
    cairo_destroy(plot->cr);

    /* trash the pango layout */
    g_object_unref(G_OBJECT(pl));

    /* record that a label is shown */
    plot->format |= PLOT_MARKERS_UP;
}

#define TOLDIST 0.01

static gint identify_point (png_plot *plot, int pixel_x, int pixel_y,
			    double x, double y) 
{
    const double *data_x = NULL;
    const double *data_y = NULL;
    double xrange, yrange;
    double xdiff, ydiff;
    double dist, mindist = NADBL;
    int best_match = -1;
    int done = 0, bank = 1;
    int t;

#if GPDEBUG > 2
    fprintf(stderr, "identify_point: pixel_x = %d (x=%g), pixel_y = %d (y=%g)\n",
	    pixel_x, x, pixel_y, y);
#endif

    if (plot->err) {
	return TRUE;
    }

    /* no markers to show */
    if (plot->spec->markers == NULL) {
	plot->status |= PLOT_NO_MARKERS;	
	return TRUE;
    }

    /* need array to keep track of which points are labeled */
    if (plot->spec->labeled == NULL) {
	plot->spec->labeled = calloc(plot->spec->nobs, 1);
	if (plot->spec->labeled == NULL) {
	    return TRUE;
	}
    }

    if (plot_is_zoomed(plot)) {
	xrange = plot->zoom_xmax - plot->zoom_xmin;
	yrange = plot->zoom_ymax - plot->zoom_ymin;
    } else {
	xrange = plot->xmax - plot->xmin;
	yrange = plot->ymax - plot->ymin;
    }

    data_x = plot->spec->data;
    data_y = data_x + plot->spec->nobs;

    if (plot_has_y2axis(plot)) {
	/* use first y-var that's on y1 axis, if any */
	int i, got_y = 0;

	for (i=0; i<plot->spec->n_lines; i++) {
	    if (plot->spec->lines[i].yaxis == 1) {
		got_y = 1;
		break;
	    }
	    if (plot->spec->lines[i].ncols > 0) {
		data_y += (plot->spec->lines[i].ncols - 1) * plot->spec->nobs;
	    }
	}
	if (!got_y) {
	    data_y = NULL;
	    plot->status |= PLOT_NO_MARKERS;	
	    return TRUE;
	}
    }

 repeat:

    /* find the best-matching data point */
    for (t=0; t<plot->spec->nobs; t++) {
	if (na(data_x[t]) || na(data_y[t])) {
	    continue;
	}
#if GPDEBUG > 4
	fprintf(stderr, "considering t=%d: x=%g, y=%g\n", t, data_x[t], data_y[t]);
#endif
	xdiff = fabs(data_x[t] - x);
	ydiff = fabs(data_y[t] - y);
	dist = sqrt(xdiff * xdiff + ydiff * ydiff);
	if (dist < mindist) {
	    mindist = dist;
	    best_match = t;
	    bank = done ? 2 : 1;
	}
    }

    if (plot->spec->code == PLOT_FACTORIZED && !done) {
	/* check the second "bank" of data also */
	data_y += plot->spec->nobs;
	done = 1;
	goto repeat;
    }

    t = best_match;
    data_y = data_x + bank * plot->spec->nobs;

    /* if the closest point is already labeled, get out */
    if (t >= 0 && plot->spec->labeled[t]) {
	return TRUE;
    }

#if GPDEBUG > 2
    fprintf(stderr, " best_match=%d, with data_x[%d]=%g, data_y[%d]=%g\n", 
	    t, t, data_x[t], t, data_y[t]);
#endif

    /* if the match is good enough, show the label */
    if (best_match >= 0) {
	xdiff = fabs(data_x[t] - x);
	ydiff = fabs(data_y[t] - y);
	if (xdiff < TOLDIST * xrange && ydiff < TOLDIST * yrange) {
	    write_label_to_plot(plot, t, pixel_x, pixel_y);
	    /* flag the point as labeled already */
	    plot->spec->labeled[t] = 1;
	}
    } 

    return TRUE;
}

#define float_fmt(i,x) ((i) && fabs(x) < 1.0e7)

static gint
plot_motion_callback (GtkWidget *widget, GdkEventMotion *event, png_plot *plot)
{
    GdkModifierType state;
    gchar label[48], label_y[24];
    const char *xfmt = NULL;
    const char *yfmt = NULL;
    int do_label;
    int x, y;

    if (plot->err) {
	return TRUE;
    }

    if (event->is_hint) {
        gdk_window_get_pointer(event->window, &x, &y, &state);
    } else {
        x = event->x;
        y = event->y;
        state = event->state;
    }

    if (plot->spec->xfmt[0] != '\0') {
	xfmt = plot->spec->xfmt;
    }

    if (plot->spec->yfmt[0] != '\0') {
	yfmt = plot->spec->yfmt;
    }

    *label = '\0';
    do_label = plot_show_cursor_label(plot);

    if (x > plot->pixel_xmin && x < plot->pixel_xmax && 
	y > plot->pixel_ymin && y < plot->pixel_ymax) {
	double data_x, data_y;

	get_data_xy(plot, x, y, &data_x, &data_y);
	if (na(data_x)) {
	    return TRUE;
	}

	if (!cant_do_labels(plot) && !labels_frozen(plot) &&
	    !plot_is_zooming(plot) &&
	    !na(data_y)) {
	    identify_point(plot, x, y, data_x, data_y);
	}

	if (do_label) {
	    if (plot->pd == 4 || plot->pd == 12) {
		x_to_date(data_x, plot->pd, label);
	    } else if (xfmt != NULL) {
		if (plot->spec->flags & GPT_TIMEFMT) {
		    date_from_gnuplot_time(label, sizeof label, 
					   xfmt, data_x);
		} else {
		    sprintf(label, xfmt, data_x);
		}
	    } else {
		sprintf(label, (float_fmt(plot->xint, data_x))? "%7.0f" : 
			"%#7.4g", data_x);
	    }
	}

	if (do_label && !na(data_y)) {
	    if (plot_has_png_coords(plot)) {
		if (yfmt != NULL) {
		    sprintf(label_y, yfmt, data_y);
		} else {
		    sprintf(label_y, (float_fmt(plot->yint, data_y))? " %-7.0f" : 
			    " %#-7.4g", data_y);
		}
	    } else {
		/* pretty much guessing at y coordinate here */
		sprintf(label_y, (float_fmt(plot->yint, data_y))? " %-7.0f" : 
				  " %#-6.3g", data_y);
	    }
	    if (strlen(label) + strlen(label_y) < sizeof label) {
		strcat(label, label_y);
	    } else {
		fprintf(stderr, "label='%s', label_y='%s'\n", label, label_y);
		strcpy(label, label_y);
	    }
	}

	if (plot_is_zooming(plot) && (state & GDK_BUTTON1_MASK)) {
	    draw_selection_box(plot, x, y);
	}
    }

    if (do_label) {
	gtk_label_set_text(GTK_LABEL(plot->cursor_label), label);
    }
  
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
    if (plot->spec->flags & GPT_Y2AXIS) {
	plot->format |= PLOT_Y2AXIS;
    }
    if (plot->spec->flags & GPT_PRINT_MARKERS) {
	plot->format |= PLOT_MARKERS_UP;
    }
    if (plot->spec->flags & GPT_POLAR) {
	plot->format |= PLOT_POLAR;
    }    
}

/* called from png plot popup menu, also toolbar item */

void start_editing_png_plot (png_plot *plot)
{
#if GPDEBUG
    fprintf(stderr, "start_editing_png_plot: plot = %p\n", (void *) plot);
#endif

    if (!PLOTSPEC_DETAILS_IN_MEMORY(plot->spec)) {
	errbox(_("Couldn't access graph info"));
	plot->err = 1;
	return;
    }

    if (plot->editor != NULL) {
	gtk_window_present(GTK_WINDOW(plot->editor));
    } else {
	plot->editor = plot_add_editor(plot);
	if (plot->editor != NULL) {
	    g_signal_connect(G_OBJECT(plot->editor), "destroy",
			     G_CALLBACK(gtk_widget_destroyed),
			     &plot->editor);
	    g_signal_connect_swapped(G_OBJECT(plot->editor), "destroy",
				     G_CALLBACK(terminate_plot_positioning),
				     plot);
	}
    }
}

static int substitute_graph_font (char *line, const gchar *fstr)
{
    char *p = strstr(line, " font ");

    if (p != NULL) {
	char tmp[256];

	*tmp = '\0';
	strncat(tmp, line, p - line + 7);
	strcat(tmp, fstr);      /* replacement font string */
	p = strchr(p + 7, '"'); /* closing quote of original font string */
	if (p != NULL) {
	    strcat(tmp, p);
	    strcpy(line, tmp);
	} 
    } 

    return (p == NULL)? E_DATA : 0;
}

/* Filter the original font spec out of the plot file, substituting
   the user's choice from a font selection dialog, then recreate
   the on-screen PNG. We do this for plots for which we can't offer
   the full-blown graph editor.
*/

void activate_plot_font_choice (png_plot *plot, const char *grfont)
{
    FILE *fp = NULL;
    FILE *ftmp = NULL;
    char tmpname[FILENAME_MAX];
    char line[256], fontname[128];
    gchar *fstr = NULL;
    int nf, fsize = 0;
    int gotterm = 0;
    int err = 0;

    nf = split_graph_fontspec(grfont, fontname, &fsize);

    if (nf == 2) {
	fstr = g_strdup_printf("%s,%d", fontname, fsize);
    } else if (nf == 1) {
	fstr = g_strdup(fontname);
    }

    fp = gretl_fopen(plot->spec->fname, "r");
    if (fp == NULL) {
	file_read_errbox(plot->spec->fname);
	return;
    }

    sprintf(tmpname, "%sgpttmp", gretl_dotdir());
    ftmp = gretl_tempfile_open(tmpname);
    if (ftmp == NULL) {
	fclose(fp);
	gui_errmsg(E_FOPEN);
	return;	
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (!gotterm && strncmp(line, "set term", 8) == 0) {
	    err = substitute_graph_font(line, fstr);
	    fputs(line, ftmp);
	    gotterm = 1;
	} else {
	    fputs(line, ftmp);
	}
    }

    fclose(fp);
    fclose(ftmp);

    if (err) {
	gretl_remove(tmpname);
    } else {
	gchar *plotcmd;

	gretl_rename(tmpname, plot->spec->fname);
	plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
				  gretl_gnuplot_path(), 
				  plot->spec->fname);
	err = gretl_spawn(plotcmd);
	g_free(plotcmd);
    }

    if (err) {
	gui_errmsg(err);
    } else {
	free(plot->spec->fontstr);
	plot->spec->fontstr = gretl_strdup(fstr);
	render_pngfile(plot, PNG_REDISPLAY);
    } 

    g_free(fstr);
}

#ifdef HAVE_AUDIO

static void audio_render_plot (png_plot *plot)
{
# ifdef G_OS_WIN32
    const char *player = NULL;
# else
    const char *player = midiplayer;
# endif
    int (*midi_play_graph) (const char *, const char *, const char *);
    void *handle;

    if (plot_not_editable(plot)) {
	return;
    }

    midi_play_graph = gui_get_plugin_function("midi_play_graph", 
					      &handle);
    if (midi_play_graph == NULL) {
        return;
    }

    (*midi_play_graph) (plot->spec->fname, gretl_dotdir(), player);

    close_plugin(handle);
}

#endif

static const gchar *menu_item_get_text (GtkMenuItem *item)
{
    GtkWidget *label = gtk_bin_get_child(GTK_BIN(item));

    return gtk_label_get_text(GTK_LABEL(label));
}

static gint color_popup_activated (GtkMenuItem *item, gpointer data)
{
    png_plot *plot = g_object_get_data(G_OBJECT(item), "plot");
    GtkWidget *parent = data;
    const gchar *item_string = NULL;
    const gchar *menu_string = NULL;

    item_string = menu_item_get_text(item);
    menu_string = menu_item_get_text(GTK_MENU_ITEM(parent));

    if (!strcmp(item_string, _("monochrome"))) {
	plot->spec->flags |= GPT_MONO;
    }

    if (!strcmp(menu_string, _("Save as Windows metafile (EMF)..."))) {
	plot->spec->termtype = GP_TERM_EMF;
	file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC, 
				  plot->spec, plot->shell);
    } 

#ifdef G_OS_WIN32
    else if (!strcmp(menu_string, _("Copy to clipboard"))) {
	win32_process_graph(plot->spec, WIN32_TO_CLIPBOARD);
    } else if (!strcmp(menu_string, _("Print"))) {
	win32_process_graph(plot->spec, WIN32_TO_PRINTER);
    }    
#endif 

    plot->spec->flags &= ~GPT_MONO;

    return TRUE;
}

static void show_numbers_from_markers (GPT_SPEC *spec)
{
    PRN *prn;
    double x, y;
    double mod, freq;
    int i, err = 0;

    if (bufopen(&prn)) {
	return;
    } 

    pputs(prn, _("roots (real, imaginary, modulus, frequency)"));
    pputs(prn, "\n\n");

    if (get_local_decpoint() != '.') {
	gretl_push_c_numeric_locale();
	for (i=0; i<spec->n_markers; i++) {
	    if (sscanf(spec->markers[i], "%lf,%lf", &x, &y) == 2) {
		freq = spec->data[i] / (2.0 * M_PI);
		mod = spec->data[spec->nobs + i];
		gretl_pop_c_numeric_locale();
		pprintf(prn, "%2d: (%7.4f  %7.4f  %7.4f  %7.4f)\n", i+1, 
			x, y, mod, freq);
		gretl_push_c_numeric_locale();
	    } else {
		err = E_DATA;
		break;
	    }
	}
	gretl_pop_c_numeric_locale();
    } else {
	for (i=0; i<spec->n_markers; i++) {
	    if (sscanf(spec->markers[i], "%lf,%lf", &x, &y) == 2) {
		freq = spec->data[i] / (2.0 * M_PI);
		mod = spec->data[spec->nobs + i];
		pprintf(prn, "%2d: (%7.4f, %7.4f, %7.4f, %7.4f)\n", i+1, 
			x, y, mod, freq);
	    } else {
		err = E_DATA;
		break;
	    }
	}
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	gchar *title = g_strdup_printf("gretl: %s", _("roots"));

	view_buffer(prn, 72, 340, title, PRINT, NULL);
	g_free(title);	
    }
}

static void boxplot_show_summary (GPT_SPEC *spec)
{
    PRN *prn = NULL;
    int err;

    if (bufopen(&prn)) {
	return;
    } 

    err = boxplot_numerical_summary(spec->fname, prn);
    
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 240, _("gretl: boxplot data"), PRINT, NULL);
    }	
}

static void add_to_session_callback (GPT_SPEC *spec)
{
    char fullname[MAXLEN] = {0};
    int err, type;

    type = (spec->code == PLOT_BOXPLOTS)? GRETL_OBJ_PLOT :
	GRETL_OBJ_GRAPH;

    err = gui_add_graph_to_session(spec->fname, fullname, type);

    if (!err) {
	remove_png_term_from_plot(fullname, spec);
	mark_plot_as_saved(spec);
    }
}

static double graph_scales[] = {
    0.8, 1.0, 1.1, 1.2, 1.4
};

int get_graph_scale (int i, double *s)
{
    int n = G_N_ELEMENTS(graph_scales);

    if (i >= 0 && i < n) {
	*s = graph_scales[i];
	return 1;
    } else {
	return 0;
    }
}

static void plot_do_rescale (png_plot *plot, int mod)
{
    int n = G_N_ELEMENTS(graph_scales);
    FILE *fp = NULL;

    if (mod == 0) {
	/* reset to default */
	if (plot->spec->scale == 1.0) {
	    /* no-op */
	    return;
	} else {
	    plot->spec->scale = 1.0;
	}
    } else {
	/* enlarge or shrink */
	int i;

	for (i=0; i<n; i++) {
	    if (plot->spec->scale == graph_scales[i]) {
		break;
	    }
	}

	if (mod == 1 && i < n - 1) {
	    plot->spec->scale = graph_scales[i+1];
	} else if (mod == -1 && i > 0) {
	    plot->spec->scale = graph_scales[i-1];
	} else {
	    gdk_window_beep(plot->window);
	    return;
	}
    }

    gnuplot_png_init(plot, &fp);

    if (fp == NULL) {
	gui_errmsg(E_FOPEN);
	return;
    }

    gtk_widget_set_sensitive(plot->up_icon, plot->spec->scale != graph_scales[n-1]);
    gtk_widget_set_sensitive(plot->down_icon, plot->spec->scale != graph_scales[0]);

    set_png_output(plot->spec);
    plotspec_print(plot->spec, fp);
    fclose(fp);
    unset_png_output(plot->spec);

    repaint_png(plot, PNG_REDISPLAY); 
}

static void show_all_labels (png_plot *plot)
{
    FILE *fp;

    if (plot->spec->labeled != NULL) {
	free(plot->spec->labeled);
	plot->spec->labeled = NULL;
    }

    plot->spec->flags |= GPT_PRINT_MARKERS;

    gnuplot_png_init(plot, &fp);

    if (fp == NULL) {
	gui_errmsg(E_FOPEN);
	return;
    }

    plotspec_print(plot->spec, fp);
    fclose(fp);

    repaint_png(plot, PNG_REDISPLAY); 
    plot->format |= PLOT_MARKERS_UP;
}

static void clear_labels (png_plot *plot)
{
    if (plot->spec->flags & GPT_PRINT_MARKERS) {
	FILE *fp;

	plot->spec->flags &= ~GPT_PRINT_MARKERS;
	gnuplot_png_init(plot, &fp);
	if (fp == NULL) {
	    gui_errmsg(E_FOPEN);
	    return;
	}
	plotspec_print(plot->spec, fp);
	fclose(fp);
    }

    repaint_png(plot, PNG_REDISPLAY); 
    plot->format &= ~PLOT_MARKERS_UP;
}

static void prepare_for_zoom (png_plot *plot)
{
    GdkCursor* cursor = gdk_cursor_new(GDK_CROSSHAIR);

    gdk_window_set_cursor(plot->window, cursor);
    gdk_cursor_unref(cursor);
    plot->status |= PLOT_ZOOMING;
    gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid, 
		       _(" Drag to define zoom rectangle"));
}

static gint plot_popup_activated (GtkMenuItem *item, gpointer data)
{
    png_plot *plot = (png_plot *) data;
    const gchar *item_string = NULL;
    int killplot = 0;

    item_string = menu_item_get_text(item);

    if (!strcmp(item_string, _("Add another curve..."))) {
	dist_graph_add(plot);
    } else if (!strcmp(item_string, _("Save as PNG..."))) {
	plot->spec->termtype = GP_TERM_PNG;
        file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC, 
				  plot->spec, plot->shell);
    } else if (!strcmp(item_string, _("Save as PDF..."))) {
	plot->spec->termtype = GP_TERM_PDF;
	pdf_ps_dialog(plot->spec, plot->shell);
    } else if (!strcmp(item_string, _("Save as postscript (EPS)..."))) {
	plot->spec->termtype = GP_TERM_EPS;
	pdf_ps_dialog(plot->spec, plot->shell);
    } else if (!strcmp(item_string, _("Save to session as icon"))) { 
	add_to_session_callback(plot->spec);
    } else if (plot_is_range_mean(plot) && !strcmp(item_string, _("Help"))) { 
	show_gui_help(RMPLOT);
    } else if (plot_is_hurst(plot) && !strcmp(item_string, _("Help"))) { 
	show_gui_help(HURST);
    } else if (!strcmp(item_string, _("Freeze data labels"))) {
	plot->spec->flags |= GPT_PRINT_MARKERS;
	redisplay_edited_plot(plot);
    } else if (!strcmp(item_string, _("Clear data labels"))) { 
	clear_labels(plot);
    } else if (!strcmp(item_string, _("All data labels"))) { 
	show_all_labels(plot);
    } else if (!strcmp(item_string, _("Zoom..."))) { 
	prepare_for_zoom(plot);
    } else if (!strcmp(item_string, _("Restore full view"))) { 
	repaint_png(plot, PNG_UNZOOM);
    } else if (!strcmp(item_string, _("Replace full view"))) {
	zoom_replaces_plot(plot);
    }
#ifndef G_OS_WIN32
    else if (!strcmp(item_string, _("Print..."))) { 
	gtk_print_graph(plot->spec->fname, plot->shell);
    }
#endif 
    else if (!strcmp(item_string, _("Display PDF"))) { 
	graph_display_pdf(plot->spec);
    } else if (!strcmp(item_string, _("OLS estimates"))) { 
	if (plot->spec != NULL) {
	    do_graph_model(plot->spec->reglist, plot->spec->fit);
	}
    } else if (!strcmp(item_string, _("Numerical values"))) {
	show_numbers_from_markers(plot->spec);
    } else if (!strcmp(item_string, _("Numerical summary"))) {
	boxplot_show_summary(plot->spec);
    } else if (!strcmp(item_string, _("Edit"))) { 
	start_editing_png_plot(plot);
    } else if (!strcmp(item_string, _("Font"))) { 
	plot_add_font_selector(plot, plot->spec->fontstr);
    } else if (!strcmp(item_string, _("Close"))) { 
        killplot = 1;
    } 

    if (killplot) {
	gtk_widget_destroy(plot->shell);
    }

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
			 G_CALLBACK(color_popup_activated), w);
	g_object_set_data(G_OBJECT(item), "plot", plot);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(cpopup), item);
    } 

    gtk_menu_item_set_submenu(GTK_MENU_ITEM(w), cpopup);
}

#define showing_all_labels(p) (p->spec != NULL && \
			       (p->spec->flags & GPT_PRINT_MARKERS) &&	\
			       p->spec->labeled == NULL)

#define graph_model_ok(f) (f == PLOT_FIT_OLS || \
                           f == PLOT_FIT_QUADRATIC || \
                           f == PLOT_FIT_INVERSE)

#define plot_not_zoomable(p) ((p->status & PLOT_DONT_ZOOM) || \
			      (p->spec != NULL && \
			       p->spec->code == PLOT_ROOTS))

static void build_plot_menu (png_plot *plot)
{
    GtkWidget *item;    
    const char *regular_items[] = {
	N_("Add another curve..."),
#ifdef G_OS_WIN32
	N_("Save as Windows metafile (EMF)..."),
#endif
	N_("Save as PNG..."),
        N_("Save as postscript (EPS)..."),
	N_("Save as PDF..."),
#ifndef G_OS_WIN32
	N_("Save as Windows metafile (EMF)..."),
#endif
#ifdef G_OS_WIN32
	N_("Copy to clipboard"),
#endif
	N_("Save to session as icon"),
	N_("Freeze data labels"),
	N_("All data labels"),
	N_("Clear data labels"),
	N_("Zoom..."),
#ifdef G_OS_WIN32
	N_("Print"),
#else
	N_("Print..."),
#endif
	N_("Display PDF"),
	N_("OLS estimates"),
	N_("Numerical values"),
	N_("Numerical summary"),
	N_("Edit"),
	N_("Font"),
	N_("Help"),
        N_("Close"),
        NULL
    };
    const char *zoomed_items[] = {
	N_("Restore full view"),
	N_("Replace full view"),
	N_("Close"),
	NULL
    };
    const char **plot_items;
    static int pdf_ok = -1;
    static int pngcairo = -1;
    int i;

    if (pdf_ok < 0) {
	pdf_ok = gnuplot_pdf_terminal();
    }

    if (pngcairo < 0) {
	pngcairo = (gnuplot_png_terminal() == GP_PNG_CAIRO);
    }    

    plot->popup = gtk_menu_new();

    if (plot_is_zoomed(plot)) {
	plot_items = zoomed_items;
    } else {
	plot_items = regular_items;
    }

    i = 0;
    while (plot_items[i]) {
	if (plot->spec->code != PLOT_PROB_DIST &&
	    !strcmp(plot_items[i], "Add another curve...")) {
	    i++;
	    continue;
	}
	if (plot_not_zoomable(plot) &&
	    !strcmp(plot_items[i], "Zoom...")) {
	    i++;
	    continue;
	}
	if (!(plot_is_range_mean(plot) || plot_is_hurst(plot)) &&
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
	if ((!pngcairo || plot_has_controller(plot) || 
	     plot_is_editable(plot)) && !strcmp(plot_items[i], "Font")) {
	    i++;
	    continue;
	}
	if (!pdf_ok && (!strcmp(plot_items[i], "Save as PDF...") ||
			!strcmp(plot_items[i], "Display PDF"))) {
	    i++;
	    continue;
	}
	if (pdf_ok && !strcmp(plot_items[i], "Print...")) {
	    /* Print... is currently very funky for graphs.  If
	       we're able to display PDF, bypass this option */
	    i++;
	    continue;
	}	    
	if (!plot_labels_shown(plot) &&
	    (!strcmp(plot_items[i], "Freeze data labels") ||
	     !strcmp(plot_items[i], "Clear data labels"))) {
	    /* no labels displayed, so these items are not relevant */
	    i++;
	    continue;
	}
	if (labels_frozen(plot) && 
	    !strcmp(plot_items[i], "Freeze data labels")) {
	    /* labels are frozen so this item inapplicable */
	    i++;
	    continue;
	}
	if ((cant_edit(plot->spec->code) || 
	     cant_do_labels(plot) || 
	     showing_all_labels(plot)) && 
	    !strcmp(plot_items[i], "All data labels")) {
	    i++;
	    continue;
	}
	if ((!plot_has_regression_list(plot) || 
	     !graph_model_ok(plot->spec->fit)) && 
	    !strcmp(plot_items[i], "OLS estimates")) {
	    i++;
	    continue;
	}
	if (!plot_is_roots(plot) && 
	    !strcmp(plot_items[i], "Numerical values")) {
	    i++;
	    continue;
	}
	if (plot->spec->code != PLOT_BOXPLOTS && 
	    !strcmp(plot_items[i], "Numerical summary")) {
	    i++;
	    continue;
	}

	item = gtk_menu_item_new_with_label(_(plot_items[i]));

	if (!strcmp(plot_items[i], "Save as Windows metafile (EMF)...") ||
	    !strcmp(plot_items[i], "Copy to clipboard") ||
	    !strcmp(plot_items[i], "Print")) {
	    /* special: items with color sub-menu */
	    attach_color_popup(item, plot);
	} else {
	    /* all other menu items */
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(plot_popup_activated),
			     plot);
	}

        gtk_widget_show(item);
        gtk_menu_shell_append(GTK_MENU_SHELL(plot->popup), item);

        i++;
    }

    g_signal_connect(G_OBJECT(plot->popup), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), 
		     &plot->popup);
}

int redisplay_edited_plot (png_plot *plot)
{
    gchar *plotcmd;
    FILE *fp;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "redisplay_edited_plot: plot = %p\n", (void *) plot);
#endif

    /* open file in which to dump plot specification */
    gnuplot_png_init(plot, &fp);
    if (fp == NULL) {
	return 1;
    }

    /* dump the edited plot details to file */
    set_png_output(plot->spec);
    plotspec_print(plot->spec, fp);
    fclose(fp);

    /* get gnuplot to create a new PNG graph */
    plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
			      gretl_gnuplot_path(), 
			      plot->spec->fname);
    err = gretl_spawn(plotcmd);
    g_free(plotcmd);

    if (err) {
	gui_errmsg(err);
	return err;
    }

    /* reset format flags */
    set_plot_format_flags(plot);

    /* grab (possibly modified) data ranges */
    get_plot_ranges(plot, plot->spec->code);

    /* put the newly created PNG onto the plot canvas */
    return render_pngfile(plot, PNG_REDISPLAY);
}

/* preparation for redisplaying graph: here we handle the case where
   we're switching to a zoomed view (by use of a temporary gnuplot
   source file); then we get gnuplot to create a new PNG.
 */

static int repaint_png (png_plot *plot, int view)
{
    char zoomname[MAXLEN];
    gchar *plotcmd = NULL;
    int err = 0;

    if (view == PNG_ZOOM) {
	FILE *fpin, *fpout;
	char line[MAXLEN];

	fpin = gretl_fopen(plot->spec->fname, "r");
	if (fpin == NULL) {
	    return 1;
	}

	build_path(zoomname, gretl_dotdir(), "zoomplot.gp", NULL);
	fpout = gretl_fopen(zoomname, "w");
	if (fpout == NULL) {
	    fclose(fpin);
	    return 1;
	}

	/* write zoomed range into auxiliary gnuplot source file */

	gretl_push_c_numeric_locale();
	fprintf(fpout, "set xrange [%g:%g]\n", plot->zoom_xmin,
		plot->zoom_xmax);
	fprintf(fpout, "set yrange [%g:%g]\n", plot->zoom_ymin,
		plot->zoom_ymax);
	gretl_pop_c_numeric_locale();

	while (fgets(line, MAXLEN-1, fpin)) {
	    if (strncmp(line, "set xrange", 10) &&
		strncmp(line, "set yrange", 10))
		fputs(line, fpout);
	}

	fclose(fpout);
	fclose(fpin);

	plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
				  gretl_gnuplot_path(), 
				  zoomname);
    } else { 
	/* PNG_UNZOOM, PNG_START or PNG_REDISPLAY */
	plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
				  gretl_gnuplot_path(),
				  plot->spec->fname);
    }

    err = gretl_spawn(plotcmd);
    g_free(plotcmd);  

    if (view == PNG_ZOOM) {
	gretl_remove(zoomname);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    return render_pngfile(plot, view);
}

/* with a zoomed version of the current plot in place,
   replace the full version wiuth the zoom
*/

static int zoom_replaces_plot (png_plot *plot)
{
    FILE *fpin, *fpout;
    char temp[MAXLEN], line[MAXLEN];
    int err = 0;

    fpin = gretl_fopen(plot->spec->fname, "r");
    if (fpin == NULL) {
	return 1;
    }

    sprintf(temp, "%szoomtmp", gretl_dotdir());
    fpout = gretl_tempfile_open(temp);
    if (fpout == NULL) {
	fclose(fpin);
	return 1;
    }

    /* write zoomed range into temporary file */

    gretl_push_c_numeric_locale();
    fprintf(fpout, "set xrange [%g:%g]\n", plot->zoom_xmin,
	    plot->zoom_xmax);
    fprintf(fpout, "set yrange [%g:%g]\n", plot->zoom_ymin,
	    plot->zoom_ymax);
    gretl_pop_c_numeric_locale();

    while (fgets(line, MAXLEN-1, fpin)) {
	if (strncmp(line, "set xrange", 10) &&
	    strncmp(line, "set yrange", 10)) {
	    fputs(line, fpout);
	}
    }

    fclose(fpout);
    fclose(fpin);

    /* and copy over original graph source file */

    err = gretl_copy_file(temp, plot->spec->fname);
    if (err) {
	gui_errmsg(err);
    } else {
	plot->xmin = plot->zoom_xmin;
	plot->xmax = plot->zoom_xmax;
	plot->ymin = plot->zoom_ymin;
	plot->ymax = plot->zoom_ymax;
	plot->zoom_xmin = plot->zoom_xmax = 0.0;
	plot->zoom_ymin = plot->zoom_ymax = 0.0;
	plot->status ^= PLOT_ZOOMED;
    }

    gretl_remove(temp);

    return err;
}

static gint plot_button_release (GtkWidget *widget, GdkEventButton *event, 
				 png_plot *plot)
{
    if (plot_is_zooming(plot)) {
	double z;

	if (!get_data_xy(plot, event->x, event->y, 
			 &plot->zoom_xmax, &plot->zoom_ymax)) {
	    return TRUE;
	}

	/* flip the selected rectangle if required */
	if (plot->zoom_xmin > plot->zoom_xmax) {
	    z = plot->zoom_xmax;
	    plot->zoom_xmax = plot->zoom_xmin;
	    plot->zoom_xmin = z;
	}

	if (plot->zoom_ymin > plot->zoom_ymax) {
	    z = plot->zoom_ymax;
	    plot->zoom_ymax = plot->zoom_ymin;
	    plot->zoom_ymin = z;
	}

	if (plot->zoom_xmin != plot->zoom_xmax &&
	    plot->zoom_ymin != plot->zoom_ymax) {
	    repaint_png(plot, PNG_ZOOM);
	}

	plot->status ^= PLOT_ZOOMING;
	gdk_window_set_cursor(plot->window, NULL);
	gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
    }

    return TRUE;
}

static gint plot_button_press (GtkWidget *widget, GdkEventButton *event, 
			       png_plot *plot)
{
    if (plot_is_zooming(plot)) {
	if (get_data_xy(plot, event->x, event->y, 
			&plot->zoom_xmin, &plot->zoom_ymin)) {
	    plot->screen_x0 = event->x;
	    plot->screen_y0 = event->y;
	}
	return TRUE;
    }

    if (plot_doing_position(plot)) {
	if (plot->pos_entry != NULL) {
	    double dx, dy;
	    
	    if (get_data_xy(plot, event->x, event->y, &dx, &dy)) {
		gchar *posstr;

		posstr = g_strdup_printf("%g %g", dx, dy);
		gtk_entry_set_text(GTK_ENTRY(plot->pos_entry), posstr);
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

    if (!plot->err) {
	if (right_click(event)) {
	    build_plot_menu(plot);
	    gtk_menu_popup(GTK_MENU(plot->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	}
    }

    return TRUE;
}

static gboolean 
plot_key_handler (GtkWidget *w, GdkEventKey *key, png_plot *plot)
{
    guint k = key->keyval;

    if (gnuplot_png_terminal() == GP_PNG_CAIRO &&
	plot_is_editable(plot) && 
	(k == GDK_plus || k == GDK_greater ||
	 k == GDK_minus || k == GDK_less ||
	 k == GDK_equal || k == GDK_0)) {
	int rk = 1;

	if (k == GDK_minus || k == GDK_less) {
	    rk = -1;
	} else if (k == GDK_0) {
	    rk = 0;
	}
	plot_do_rescale(plot, rk);
	return TRUE;
    }

    switch (k) {
    case GDK_q:
    case GDK_Q:
#ifdef MAC_NATIVE
    case GDK_w:
    case GDK_W:
#endif
	gtk_widget_destroy(w);
	break;
    case GDK_s:
    case GDK_S:
	add_to_session_callback(plot->spec);
	break;
    case GDK_z:
    case GDK_Z:	
	prepare_for_zoom(plot);
	break;
#ifdef G_OS_WIN32
    case GDK_c:
	win32_process_graph(plot->spec, WIN32_TO_CLIPBOARD);
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

#if GTK_MAJOR_VERSION >= 3

static 
void plot_draw (GtkWidget *canvas, cairo_t *cr, gpointer data)
{
    png_plot *plot = data;

    cairo_set_source_surface(cr, plot->cs, 0, 0);
    cairo_paint(cr);
}

#else /* transitional use of cairo */

static 
void plot_expose (GtkWidget *canvas, GdkEventExpose *event, 
		  gpointer data)
{
    png_plot *plot = data;

    plot->cr = gdk_cairo_create(plot->window);
    gdk_cairo_set_source_pixmap(plot->cr, plot->pixmap, 0, 0);
    gdk_cairo_rectangle(plot->cr, &event->area);
    cairo_fill(plot->cr);
    cairo_destroy(plot->cr);
}

#endif

static GdkPixbuf *gretl_pixbuf_new_from_file (const gchar *fname)
{
    GdkPixbuf *pbuf = NULL;
    gchar *trfname = NULL;
    int err;

    err = validate_filename_for_glib(fname, &trfname);

    if (!err) {
	GError *gerr = NULL;

	if (trfname != NULL) {
	    pbuf = gdk_pixbuf_new_from_file(trfname, &gerr);
	    g_free(trfname);
	} else {
	    pbuf = gdk_pixbuf_new_from_file(fname, &gerr);
	}

	if (gerr != NULL) {
	    errbox(gerr->message);
	    g_error_free(gerr);
	}
    }    

    return pbuf;
}

static int resize_png_plot (png_plot *plot, int width, int height)
{
    png_bounds b;

    plot->pixel_width = width;
    plot->pixel_height = height;

    gtk_widget_set_size_request(GTK_WIDGET(plot->canvas), 
				plot->pixel_width, plot->pixel_height);

#if GTK_MAJOR_VERSION == 2
    g_object_unref(plot->pixmap);
    plot->pixmap = gdk_pixmap_new(plot->window, 
				  plot->pixel_width, 
				  plot->pixel_height, 
				  -1);
#endif

    if (plot->status & (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE)) {
	return 0;
    }

    /* try revising the gnuplot bounds info? */

    if (plot_has_png_coords(plot) && 
	get_png_bounds_info(&b) == GRETL_PNG_OK) {
	plot->status |= PLOT_PNG_COORDS;
	plot->pixel_xmin = b.xleft;
	plot->pixel_xmax = b.xright;
	plot->pixel_ymin = plot->pixel_height - b.ytop;
	plot->pixel_ymax = plot->pixel_height - b.ybot;
	plot->xmin = b.xmin;
	plot->xmax = b.xmax;
	plot->ymin = b.ymin;
	plot->ymax = b.ymax;
    } else {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
    }

    return 0;
}

/* The last step in displaying a graph (or redisplaying after some
   change has been made): grab the gnuplot-generated PNG file, make a
   pixbuf out of it, and draw the pixbuf onto the canvas of the plot
   window.
*/

static int render_pngfile (png_plot *plot, int view)
{
    gint width, height;
    GdkPixbuf *pbuf;
    char pngname[MAXLEN];

    build_path(pngname, gretl_dotdir(), "gretltmp.png", NULL);
    pbuf = gretl_pixbuf_new_from_file(pngname);

    if (pbuf == NULL) {
	gretl_remove(pngname);
	return 1;
    }

    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);

    if (width == 0 || height == 0) {
	errbox(_("Malformed PNG file for graph"));
	g_object_unref(pbuf);
	gretl_remove(pngname);
	return 1;
    }

    if (width != plot->pixel_width || height != plot->pixel_height) {
	resize_png_plot(plot, width, height);
    }

    /* scrap any old record of which points are labeled */
    if (plot->spec->labeled != NULL) {
	free(plot->spec->labeled);
	plot->spec->labeled = NULL;
	if (!(plot->spec->flags & GPT_PRINT_MARKERS)) {
	    /* any markers will have disappeared on reprinting */
	    plot->format &= ~PLOT_MARKERS_UP;
	} 
    }

#if GTK_MAJOR_VERSION >= 3
    copy_pixbuf_to_surface(plot, pbuf);
#else
    plot->cr = gdk_cairo_create(plot->pixmap);
    gdk_cairo_set_source_pixbuf(plot->cr, pbuf, 0, 0);
    cairo_paint(plot->cr);
    cairo_destroy(plot->cr);
    if (plot->savebuf != NULL) {
	g_object_unref(plot->savebuf);
	plot->savebuf = NULL;
    }
#endif

    g_object_unref(pbuf);
    gretl_remove(pngname);
   
    if (view != PNG_START) { 
	/* we're changing the view, so refresh the whole canvas */
	redraw_plot_rectangle(plot, NULL);	
	if (view == PNG_ZOOM) {
	    plot->status |= PLOT_ZOOMED;
	} else if (view == PNG_UNZOOM) {
	    plot->status ^= PLOT_ZOOMED;
	}
    }

#ifdef G_OS_WIN32
    /* somehow the plot can end up underneath */
    gtk_window_present(GTK_WINDOW(plot->shell));
#endif

    return 0;
}

static void destroy_png_plot (GtkWidget *w, png_plot *plot)
{
    /* delete temporary plot source file? */
    if (!plot_is_saved(plot)) {
	gretl_remove(plot->spec->fname);
    }

#if GPDEBUG
    fprintf(stderr, "destroy_png_plot: plot = %p, spec = %p\n",
	    (void *) plot, (void *) plot->spec);
#endif

    plotspec_destroy(plot->spec);

#if GTK_MAJOR_VERSION >= 3
    if (plot->cs != NULL) {
	cairo_surface_destroy(plot->cs);
    }
#else
    if (plot->savebuf != NULL) {
	g_object_unref(plot->savebuf);
    }
#endif

    g_object_unref(plot->shell);

    free(plot);
}

/* Do a partial parse of the gnuplot source file: enough to determine
   the data ranges so we can read back the mouse pointer coordinates
   when the user moves the pointer over the graph.
*/

static int get_plot_ranges (png_plot *plot, PlotType ptype)
{
    FILE *fp;
    char line[MAXLEN];
    int got_x = 0;
    int got_y = 0;
    png_bounds b;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "get_plot_ranges: plot=%p, plot->spec=%p\n", 
	    (void *) plot, (void *) plot->spec);
#endif    

    plot->xmin = plot->xmax = 0.0;
    plot->ymin = plot->ymax = 0.0;   
    plot->xint = plot->yint = 0;

    if (no_readback(plot->spec->code)) {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
	return 1;
    }

    fp = gretl_fopen(plot->spec->fname, "r");
    if (fp == NULL) {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
	return 1;
    }

    gretl_push_c_numeric_locale();

    while (fgets(line, MAXLEN-1, fp) && strncmp(line, "plot ", 5)) {
	if (sscanf(line, "set xrange [%lf:%lf]", 
		   &plot->xmin, &plot->xmax) == 2) { 
	    got_x = 1;
	} 
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    /* now try getting accurate coordinate info from 
       auxiliary file (or maybe PNG file)
    */
    if (get_png_bounds_info(&b) == GRETL_PNG_OK) {
	plot->status |= PLOT_PNG_COORDS;
	got_x = got_y = 1;
	plot->pixel_xmin = b.xleft;
	plot->pixel_xmax = b.xright;
	plot->pixel_ymin = plot->pixel_height - b.ytop;
	plot->pixel_ymax = plot->pixel_height - b.ybot;
	plot->xmin = b.xmin;
	plot->xmax = b.xmax;
	plot->ymin = b.ymin;
	plot->ymax = b.ymax;
# if POINTS_DEBUG
	fprintf(stderr, "get_png_bounds_info():\n"
		" xmin=%d xmax=%d ymin=%d ymax=%d\n", 
		plot->pixel_xmin, plot->pixel_xmax,
		plot->pixel_ymin, plot->pixel_ymax);
	fprintf(stderr, "using px_height %d, px_width %d\n",
		plot->pixel_height, plot->pixel_width);
# endif
	fprintf(stderr, "get_png_bounds_info(): OK\n");
    } else {
	fprintf(stderr, "get_png_bounds_info(): failed\n");
    }

    /* If got_x = 0 at this point, we didn't get an x-range out of 
       gnuplot, so we might as well give up.
    */
    if (!got_x) {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
	return 1;
    }    

    if (plot_has_png_coords(plot)) {
	plot->status |= PLOT_HAS_XRANGE;
	plot->status |= PLOT_HAS_YRANGE;
	if (ptype != PLOT_ROOTS && ptype != PLOT_QQ) {
	    plot->status |= PLOT_CURSOR_LABEL;
	}
	if ((plot->xmax - plot->xmin) / 
	    (plot->pixel_xmax - plot->pixel_xmin) >= 1.0) {
	    plot->xint = 1;
	}
	if ((plot->ymax - plot->ymin) / 
	    (plot->pixel_ymax - plot->pixel_ymin) >= 1.0) {
	    plot->yint = 1;
	}
    } else {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
#if POINTS_DEBUG 
	fputs("get_plot_ranges: setting PLOT_DONT_ZOOM, PLOT_DONT_MOUSE\n", 
	      stderr);
#endif
    }

    return err;
}

static png_plot *png_plot_new (void)
{
    png_plot *plot = mymalloc(sizeof *plot);

    if (plot == NULL) {
	return NULL;
    }

    plot->shell = NULL;
    plot->canvas = NULL;
    plot->popup = NULL;
    plot->statusarea = NULL;    
    plot->statusbar = NULL;
    plot->cursor_label = NULL;
    plot->cr = NULL;
#if GTK_MAJOR_VERSION >= 3
    plot->cs = NULL;
#else
    plot->pixmap = NULL;
    plot->savebuf = NULL;
#endif
    plot->spec = NULL;
    plot->editor = NULL;
    plot->window = NULL;
    plot->up_icon = NULL;
    plot->down_icon = NULL;

    plot->pixel_width = GP_WIDTH;
    plot->pixel_height = GP_HEIGHT;

    plot->xmin = plot->xmax = 0.0;
    plot->ymin = plot->ymax = 0.0;
    plot->xint = plot->yint = 0;

    plot->zoom_xmin = plot->zoom_xmax = 0.0;
    plot->zoom_ymin = plot->zoom_ymax = 0.0;
    plot->screen_x0 = plot->screen_y0 = 0;
    plot->pd = 0;
    plot->err = 0;
    plot->cid = 0;
    plot->status = 0;
    plot->format = 0;

    return plot;
}

/* note: @fname is the name of the file containing the 
   plot commands.
*/

static int gnuplot_show_png (const char *fname, const char *name,
			     GPT_SPEC *spec, int saved)
{
    GtkWidget *vbox;
    GtkWidget *canvas_hbox;
    GtkWidget *status_hbox = NULL;
    gchar *title = NULL;
    png_plot *plot;
    int err = 0;

    if (*fname == '\0') {
	return 0;
    }

    gretl_error_clear(); 

#if GPDEBUG
    fprintf(stderr, "gnuplot_show_png:\n fname='%s', spec=%p, saved=%d\n",
	    fname, (void *) spec, saved);
#endif

    plot = png_plot_new();
    if (plot == NULL) {
	return E_ALLOC;
    }

    if (spec != NULL) {
	plot->spec = spec;
    } else {
	plot->spec = plotspec_new();
	if (plot->spec == NULL) {
	    free(plot);
	    return E_ALLOC;
	}
	strcpy(plot->spec->fname, fname);
    }

    if (saved) {
	plot->status |= PLOT_SAVED;
    }

    /* make png plot struct accessible via spec */
    plot->spec->ptr = plot;

    /* Parse the gnuplot source file.  If we hit errors here,
       flag this, but it's not necessarily a show-stopper in
       terms of simply displaying the graph -- unless we get
       E_FOPEN.
    */
    plot->err = read_plotspec_from_file(plot->spec, &plot->pd);

    if (plot->err == E_FOPEN) {
	if (plot->spec != spec) {
	    plotspec_destroy(plot->spec);
	}
	free(plot);
	return E_FOPEN;
    }

#if GPDEBUG 
    fprintf(stderr, "gnuplot_show_png: read_plotspec_from_file returned %d\n",
	    plot->err);
#endif

    if (plot->err || cant_edit(plot->spec->code)) {
	plot->status |= (PLOT_DONT_EDIT | PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
    } else {
	set_plot_format_flags(plot);
    } 

    if (plot->spec->code == PLOT_ROOTS ||
	plot->spec->code == PLOT_QQ) {
	plot->pixel_width = plot->pixel_height = GP_SQ_SIZE;
    }

    if (plot->spec->flags & GPT_LETTERBOX) {
	plot->pixel_width = GP_LB_WIDTH;
	plot->pixel_height = GP_LB_HEIGHT;
    } else if (plot->spec->flags & GPT_XXL) {
	plot->pixel_width = GP_XXL_WIDTH;
	plot->pixel_height = GP_XXL_HEIGHT;
    } else if (plot->spec->flags & GPT_XL) {
	plot->pixel_width = GP_XL_WIDTH;
	plot->pixel_height = GP_XL_HEIGHT;
    }

    if (plot->spec->scale != 1.0) {
	plot_get_scaled_dimensions(&plot->pixel_width,
				   &plot->pixel_height,
				   plot->spec->scale);
    }

    if (!plot->err) {
	int range_err = get_plot_ranges(plot, plot->spec->code);

#if GPDEBUG
	fprintf(stderr, "range_err = %d\n", range_err);
#endif
	if (plot->spec->nbars > 0) {
	    if (range_err) {
		plot->spec->nbars = 0;
	    } else {
		plotspec_set_bars_limits(plot->spec, 
					 plot->xmin, plot->xmax,
					 plot->ymin, plot->ymax);
	    }
	}
    }

    plot->shell = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    /* note need for corresponding unref */
    g_object_ref(plot->shell);

    if (name != NULL) {
	title = g_strdup_printf("gretl: %s", name);
    } else {
	title = g_strdup(_("gretl: graph"));
    }

    gtk_window_set_title(GTK_WINDOW(plot->shell), title);
    gtk_window_set_resizable(GTK_WINDOW(plot->shell), FALSE);
    g_signal_connect(G_OBJECT(plot->shell), "destroy",
		     G_CALLBACK(destroy_png_plot), plot);
    g_free(title);

    vbox = gtk_vbox_new(FALSE, 2);
    gtk_container_add(GTK_CONTAINER(plot->shell), vbox);

    /* box to hold canvas */
    canvas_hbox = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), canvas_hbox, TRUE, TRUE, 0);

    /* eventbox and hbox for status area  */
    plot->statusarea = gtk_event_box_new();
    gtk_box_pack_start(GTK_BOX(vbox), plot->statusarea, FALSE, FALSE, 0);
    status_hbox = gtk_hbox_new(FALSE, 2);
    gtk_container_add(GTK_CONTAINER(plot->statusarea), status_hbox);

    /* Create drawing-area widget */
    plot->canvas = gtk_drawing_area_new();
    gtk_widget_set_size_request(GTK_WIDGET(plot->canvas), 
				plot->pixel_width, plot->pixel_height);
    gtk_widget_set_events (plot->canvas, GDK_EXPOSURE_MASK
                           | GDK_LEAVE_NOTIFY_MASK
                           | GDK_BUTTON_PRESS_MASK
                           | GDK_BUTTON_RELEASE_MASK
                           | GDK_POINTER_MOTION_MASK
                           | GDK_POINTER_MOTION_HINT_MASK);
    gtk_widget_set_can_focus(plot->canvas, TRUE);
    g_signal_connect(G_OBJECT(plot->canvas), "button-press-event", 
		     G_CALLBACK(plot_button_press), plot);
    g_signal_connect(G_OBJECT(plot->canvas), "button-release-event", 
		     G_CALLBACK(plot_button_release), plot);
    gtk_box_pack_start(GTK_BOX(canvas_hbox), plot->canvas, FALSE, FALSE, 0);

    if (plot_show_cursor_label(plot)) {
	/* cursor label (graph position indicator) */
	GtkWidget *frame = gtk_frame_new(NULL);

	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	plot->cursor_label = gtk_label_new(" ");
	gtk_widget_set_size_request(plot->cursor_label, 160, -1);
	gtk_container_add(GTK_CONTAINER(frame), plot->cursor_label);
	gtk_box_pack_start(GTK_BOX(status_hbox), frame, FALSE, FALSE, 0);
    }

    /* the statusbar */
    plot->statusbar = gtk_statusbar_new();
    gtk_box_pack_start(GTK_BOX(status_hbox), plot->statusbar, TRUE, TRUE, 0);
    add_graph_toolbar(status_hbox, plot);

#if GTK_MAJOR_VERSION < 3
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(plot->statusbar), FALSE);
#endif
    plot->cid = gtk_statusbar_get_context_id(GTK_STATUSBAR(plot->statusbar),
					     "plot_message");
    if (!plot->err) {
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar),
			   plot->cid, _(" Right-click on graph for menu"));
    }
    
    if (plot_has_xrange(plot)) {
	g_signal_connect(G_OBJECT(plot->canvas), "motion-notify-event",
			 G_CALLBACK(plot_motion_callback), plot);
    }

    gtk_widget_realize(plot->canvas);
    plot->window = gtk_widget_get_window(plot->canvas);

#if GTK_MAJOR_VERSION < 3
    gdk_window_set_back_pixmap(plot->window, NULL, FALSE);
#endif

    /* finish setup of plot->shell */
    g_object_set_data(G_OBJECT(plot->shell), "plot-filename", 
		      plot->spec->fname);
    window_list_add(plot->shell, GNUPLOT);
    g_signal_connect(G_OBJECT(plot->shell), "key-press-event", 
		     G_CALLBACK(plot_key_handler), plot);
    gtk_widget_show_all(plot->shell);

#if 0
    /* set the focus to the canvas area */
    gtk_widget_grab_focus(plot->canvas);
#endif

#if GTK_MAJOR_VERSION >= 3
    g_signal_connect(G_OBJECT(plot->canvas), "draw",
		     G_CALLBACK(plot_draw), plot);
#else
    plot->pixmap = gdk_pixmap_new(plot->window, 
				  plot->pixel_width, 
				  plot->pixel_height, 
				  -1);
    g_signal_connect(G_OBJECT(plot->canvas), "expose-event",
		     G_CALLBACK(plot_expose), plot);
#endif

    err = render_pngfile(plot, PNG_START);

    if (err) {
	gtk_widget_destroy(plot->shell);
    } else {
	g_object_set_data(G_OBJECT(plot->shell), "object", plot);
    }

    return err;
}

/* Called on a newly created PNG graph; note that
   the filename returned by gretl_plotfile() is the
   input file (set of gnuplot commands).
*/

void register_graph (void)
{
    gnuplot_show_png(gretl_plotfile(), NULL, NULL, 0);
}

/* @fname is the name of a plot command file from the
   current session, and @title is its display name */

void display_session_graph (const char *fname, const char *title) 
{
    char fullname[MAXLEN];
    gchar *plotcmd;
    int err = 0;

    if (g_path_is_absolute(fname)) {
	strcpy(fullname, fname);
    } else {
	sprintf(fullname, "%s%s", gretl_dotdir(), fname);
    }

    if (add_png_term_to_plot(fullname)) {
	return;
    }

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
			      gretl_gnuplot_path(), 
			      fullname);
    err = gretl_spawn(plotcmd);
    g_free(plotcmd);

    if (err) {
	/* display the bad plot file */
	view_file(fullname, 0, 0, 78, 350, VIEW_FILE);
    } else {
	err = gnuplot_show_png(fullname, title, NULL, 1);
    }

    if (err) {
	gui_errmsg(err);
    }
}

static int get_png_plot_bounds (const char *str, png_bounds *bounds)
{
    int ret = GRETL_PNG_OK;

    bounds->xleft = bounds->xright = 0;
    bounds->ybot = bounds->ytop = 0;

    if (sscanf(str, "pixel_bounds: %d %d %d %d",
	       &bounds->xleft, &bounds->xright,
	       &bounds->ybot, &bounds->ytop) != 4) {
	ret = GRETL_PNG_BAD_COMMENTS;
    } 

    if (ret == GRETL_PNG_OK && bounds->xleft == 0 && 
	bounds->xright == 0 && bounds->ybot == 0 && 
	bounds->ytop == 0) {
	ret = GRETL_PNG_NO_COORDS;
    }

#if POINTS_DEBUG
    fprintf(stderr, "Got: xleft=%d, xright=%d, ybot=%d, ytop=%d\n",
	    bounds->xleft, bounds->xright, bounds->ybot, bounds->ytop);
#endif

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

    bounds->xmin = bounds->xmax = 0.0;
    bounds->ymin = bounds->ymax = 0.0;

    gretl_push_c_numeric_locale();

    if (sscanf(str, "data_bounds: %lf %lf %lf %lf",
	       &bounds->xmin, &bounds->xmax,
	       &bounds->ymin, &bounds->ymax) != 4) {
	ret = GRETL_PNG_BAD_COMMENTS;
    } 

    if (ret == GRETL_PNG_OK && bounds->xmin == 0.0 && 
	bounds->xmax == 0.0 && bounds->ymin == 0.0 && 
	bounds->ymax == 0.0) {
	ret = GRETL_PNG_NO_COORDS;
    } 

#if POINTS_DEBUG
    fprintf(stderr, "Got: xmin=%g, xmax=%g, ymin=%g, ymax=%g\n",
	    bounds->xmin, bounds->xmax, bounds->ymin, bounds->ymax);
#endif

    gretl_pop_c_numeric_locale();

    return ret;
}

static int get_png_bounds_info (png_bounds *bounds)
{
    FILE *fp;
    char bbname[MAXLEN], line[128];
    int plot_ret = -1, data_ret = -1;
    int ret = GRETL_PNG_OK;

    build_path(bbname, gretl_dotdir(), "gretltmp.png", ".bounds"); 
    fp = gretl_fopen(bbname, "r");

    if (fp == NULL) {
	fprintf(stderr, "couldn't open %s\n", bbname);
	return GRETL_PNG_NO_COMMENTS;
    }

    if (fgets(line, sizeof line, fp) == NULL) {
	plot_ret = GRETL_PNG_NO_COMMENTS;
    } else {
	plot_ret = get_png_plot_bounds(line, bounds);
    }

    if (fgets(line, sizeof line, fp) == NULL) {
	data_ret = GRETL_PNG_NO_COMMENTS;
    } else {
	data_ret = get_png_data_bounds(line, bounds);
    }

    if (plot_ret == GRETL_PNG_NO_COORDS && data_ret == GRETL_PNG_NO_COORDS) {
	/* comments were present and correct, but all zero */
	ret = GRETL_PNG_NO_COORDS;
    } else if (plot_ret != GRETL_PNG_OK || data_ret != GRETL_PNG_OK) {
	/* one or both set of coordinates bad or missing */
	if (plot_ret >= 0 || data_ret >= 0) {
	    ret = GRETL_PNG_BAD_COMMENTS;
	} else {
	    ret = GRETL_PNG_NO_COMMENTS;
	}
    }

    fclose(fp);
    gretl_remove(bbname);

    return ret;
}

#ifndef G_OS_WIN32

#include <errno.h>

static int get_terminal (char *s)
{
    const gchar *terms[] = {
	"xterm",
	"rxvt",
	"gnome-terminal",
	"kterm",
	"urxvt",
	NULL
    };
    gchar *test;
    int i;

    for (i=0; terms[i] != NULL; i++) {
	test = g_find_program_in_path(terms[i]);
	if (test != NULL) {
	    g_free(test);
	    strcpy(s, terms[i]);
	    return 0;
	}
    }

    errbox(_("Couldn't find a usable terminal program"));

    return 1;
}

#endif /* !G_OS_WIN32 */

void launch_gnuplot_interactive (const char *plotfile)
{
#if defined(G_OS_WIN32)
    gchar *gpline;

    if (plotfile == NULL) {
	gpline = g_strdup_printf("\"%s\"", gretl_gnuplot_path());
    } else {
	gpline = g_strdup_printf("\"%s\" \"%s\" -", 
				 gretl_gnuplot_path(),
				 plotfile);
    }	
    create_child_process(gpline);
    g_free(gpline);
#elif defined(MAC_NATIVE)
    gchar *gpline;

    if (plotfile == NULL) {
	gpline = g_strdup_printf("open -a Terminal.app \"%s.sh\"",
				 gretl_gnuplot_path());
    } else {
	gpline = g_strdup_printf("open -a Terminal.app \"%s.sh\" \"%s\"",
				 gretl_gnuplot_path(), plotfile);
    }	
    system(gpline);
    g_free(gpline);    
#else 
    char term[16];
    char fname[MAXLEN];
    int err = 0;

    if (plotfile != NULL) {
	strcpy(fname, plotfile);
    } else {
	*fname = '\0';
    }

    if (*fname != '\0' && gnuplot_has_wxt()) {
	*term = '\0';
    } else {
	err = get_terminal(term);
    }

    if (!err) {
	const char *gp = gretl_gnuplot_path();
	GError *error = NULL;
	gchar *argv[12];
	int ok;

	if (*term == '\0') {
	    /* no controller is needed */
	    argv[0] = (char *) gp;
	    argv[1] = fname;
	    argv[2] = "-persist";
	    argv[3] = NULL;
	} else if (strstr(term, "gnome")) {
	    /* gnome-terminal */
	    argv[0] = term;
	    if (*fname != '\0') {
		argv[1] = "--geometry=40x4";
		argv[2] = "--title=\"gnuplot: type q to quit\"";
		argv[3] = "-x";
		argv[4] = (char *) gp;
		argv[5] = fname;
		argv[6] = "-";
		argv[7] = NULL;
	    } else {
		argv[1] = "--title=\"gnuplot: type q to quit\"";
		argv[2] = "-x";
		argv[3] = (char *) gp;
		argv[4] = NULL;
	    }		
	} else {	    
	    /* xterm, rxvt, kterm */
	    argv[0] = term;
	    if (*fname != '\0') {
		argv[1] = "+sb";
		argv[2] = "+ls";
		argv[3] = "-geometry";
		argv[4] = "40x4";
		argv[5] = "-title";
		argv[6] = "gnuplot: type q to quit";
		argv[7] = "-e";
		argv[8] = (char *) gp;
		argv[9] = fname;
		argv[10] = "-";
		argv[11] = NULL;
	    } else {
		argv[1] = "-title";
		argv[2] = "gnuplot: type q to quit";
		argv[3] = "-e";
		argv[4] = (char *) gp;
		argv[5] = NULL;
	    }
	} 

	ok = g_spawn_async(NULL, /* working dir */
			   argv,
			   NULL, /* env */
			   G_SPAWN_SEARCH_PATH,
			   NULL, /* child_setup */
			   NULL, /* user_data */
			   NULL, /* child_pid ptr */
			   &error);
	if (!ok) {
	    errbox(error->message);
	    g_error_free(error);
	} 
    }
#endif /* !(G_OS_WIN32 or MAC_NATIVE) */
}
