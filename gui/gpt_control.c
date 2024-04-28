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
#include "libset.h"
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
#include "toolbar.h"
#include "clipboard.h"

#ifdef G_OS_WIN32
# include <io.h>
# include "gretlwin32.h"
#endif

#ifdef OS_OSX
# include "osx_open.h"
#endif

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk/gdkkeysyms.h>
#include <errno.h>

#define GPDEBUG 0
#define WINDEBUG 0
#define POINTS_DEBUG 0
#define COLLDEBUG 0
#define LT_DEBUG 0

/* the following needs more testing */
#define HANDLE_HEREDATA 1

/* pager for plot collection: use spin button? */
#define SPIN_PAGER 0

enum {
    PLOT_SAVED          = 1 << 0,
    PLOT_ZOOMED         = 1 << 1,
    PLOT_ZOOMING        = 1 << 2,
    PLOT_PNG_COORDS     = 1 << 3,
    PLOT_HAS_XRANGE     = 1 << 4,
    PLOT_HAS_YRANGE     = 1 << 5,
    PLOT_DONT_ZOOM      = 1 << 6,
    PLOT_DONT_EDIT      = 1 << 7,
    PLOT_DONT_MOUSE     = 1 << 8,
    PLOT_POSITIONING    = 1 << 9,
    PLOT_CURSOR_LABEL   = 1 << 10,
    PLOT_TERM_HIDDEN    = 1 << 11
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
#define plot_is_geomap(p)       (p->spec->code == PLOT_GEOMAP)

#define plot_has_regression_list(p) (p->spec->reglist != NULL)

#define labels_frozen(p)        (p->spec->flags & GPT_PRINT_MARKERS)
#define cant_do_labels(p)       (p->err || p->spec->markers == NULL)

#define plot_show_cursor_label(p) (p->status & PLOT_CURSOR_LABEL)

#define plot_has_controller(p) (p->editor != NULL)

typedef enum {
    PNG_START,
    PNG_ZOOM,
    PNG_UNZOOM,
    PNG_REDISPLAY,
    PNG_REPLACE
} viewcode;

#if COLLDEBUG
static const char *viewstrs[] = {
    "START", "ZOOM", "UNZOOM", "REDISPLAY", "REPLACE"
};

static const char *viewstr (viewcode i)
{
    return viewstrs[i];
}
#endif

struct multiplot_ {
    gint64 mtime;
    GList *list;
    int current;
    int id;
#if SPIN_PAGER
    GtkWidget *sb;
#endif
};

typedef struct multiplot_ multiplot;

struct png_plot_t {
    GtkWidget *shell;     /* top-level GTK window */
    multiplot *mp;        /* collection-specific info */
    GtkWidget *canvas;    /* area in which plot is drawn */
    GtkWidget *popup;     /* transient popup menu */
    GtkWidget *statusbar;
    GtkWidget *cursor_label;
    GtkWidget *editor;    /* state-dependent pointer */
    GtkWidget *pos_entry; /* state-dependent pointer */
    GtkWidget *toolbar;   /* toolbar at bottom right of window */
    GtkWidget *up_icon;   /* resizing button */
    GtkWidget *down_icon; /* resizing button */
    GdkWindow *window;
    GdkPixbuf *pbuf;
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
    guint32 status;
    guint8 format;
};

static int render_png (png_plot *plot, viewcode view);
static int repaint_png (png_plot *plot, int view);
static int zoom_replaces_plot (png_plot *plot);
static void prepare_for_zoom (png_plot *plot);
static int get_plot_ranges (png_plot *plot, PlotType ptype);
static void graph_display_pdf (png_plot *plot);
static void plot_do_rescale (png_plot *plot, int mod);
#ifdef G_OS_WIN32
static void win32_process_graph (png_plot *plot, int dest);
#else
static void set_plot_for_copy (png_plot *plot);
#endif
static void build_plot_menu (png_plot *plot);
static int plot_add_shell (png_plot *plot, const char *name);
static GdkPixbuf *pixbuf_from_file (png_plot *plot);
static int resize_png_plot (png_plot *plot, int width, int height,
			    int follower);
static png_plot *widget_get_plot (gpointer p);
static void destroy_png_plot (GtkWidget *w, png_plot *plot);
#if GTK_MAJOR_VERSION == 2
static void plot_nullify_surface (png_plot *plot);
#endif

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
    int width;
    int height;
};

typedef struct linestyle_ linestyle;

struct linestyle_ {
    char lc[20]; /* color */
    float lw;    /* line width */
    int dt;      /* dash type */
    int pt;      /* point type */
};

static int get_png_bounds_info (png_bounds *bounds);

#define PLOTSPEC_DETAILS_IN_MEMORY(s) (s->lines != NULL)

/* file-scope globals */
static png_plot *plot_collection;
static int collection_count;

#include "gpt_collect.c"

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

static void plot_invalidate_pixbuf (png_plot *plot)
{
    if (plot->pbuf != NULL) {
	g_object_unref(G_OBJECT(plot->pbuf));
	plot->pbuf = NULL;
    }
}

static png_plot *plot_get_current (png_plot *plot)
{
    if (in_collection(plot)) {
	return g_list_nth_data(plot->mp->list, plot->mp->current);
    } else {
	return plot;
    }
}

static png_plot *widget_get_plot (gpointer p)
{
    return g_object_get_data(G_OBJECT(p), "plot");
}

gboolean is_shell_for_plotfile (GtkWidget *w,
				const char *fname)
{
    png_plot *plot = widget_get_plot(w);

    return plot != NULL && plot->spec != NULL &&
	!strcmp(fname, plot->spec->fname);
}

GPT_SPEC *plot_get_spec (png_plot *plot)
{
    return plot->spec;
}

/* called from calculator.c */

GtkWidget *plot_get_shell (png_plot *plot)
{
    return plot->shell;
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

void plot_get_pixel_dims (png_plot *plot,
			  double *pw,
			  double *ph)
{
    *pw = plot->pixel_width;
    *ph = plot->pixel_height;
}

/* apparatus for graph toolbar */

static void graph_rescale_callback (GtkWidget *w, png_plot *plot)
{
    int mod = (w == plot->up_icon)? 1 : -1;

    plot_do_rescale(plot, mod);
}

static void graph_popup_callback (GtkWidget *w, png_plot *plot)
{
    build_plot_menu(plot);
    gtk_menu_popup(GTK_MENU(plot->popup), NULL, NULL, NULL, NULL,
		   1, gtk_get_current_event_time());
}

static void plot_winlist_popup (GtkWidget *w, png_plot *plot)
{
    window_list_popup(w, NULL, plot->shell);
}

static GretlToolItem plotbar_items[] = {
    { N_("Menu"),        GRETL_STOCK_MENU,    G_CALLBACK(graph_popup_callback), 0 },
    { N_("Bigger"),      GRETL_STOCK_BIGGER,  G_CALLBACK(graph_rescale_callback), 32 },
    { N_("Smaller"),     GRETL_STOCK_SMALLER, G_CALLBACK(graph_rescale_callback), 0 },
    { N_("Windows"),     GRETL_STOCK_WINLIST, G_CALLBACK(plot_winlist_popup), 0 }
};

static void add_graph_toolbar (GtkWidget *hbox, png_plot *plot)
{
    GtkWidget *tbar, *button, **pb;
    GretlToolItem *item;
    int i, n = G_N_ELEMENTS(plotbar_items);

    plot->toolbar = tbar = gretl_toolbar_new(NULL);

    for (i=0; i<n; i++) {
	pb = NULL;
	item = &plotbar_items[i];
	if (item->func == G_CALLBACK(graph_rescale_callback)) {
	    if (plot_not_editable(plot)) {
		continue;
	    } else {
		pb = (item->flag == 32)? &plot->up_icon :
		    &plot->down_icon;
	    }
	}
	button = gretl_toolbar_insert(tbar, item, item->func, plot, -1);
	if (pb != NULL) {
	    *pb = button;
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), tbar, FALSE, FALSE, 5);
}

/* end apparatus for graph toolbar */

/* Provide the data coordinates for a gretl/gnuplot
   graph, if they are all positive, otherwise return
   non-zero. Called by the plot editor.
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
	GdkCursor *cursor;

	cursor = gdk_cursor_new(GDK_CROSSHAIR);
	if (cursor != NULL) {
	    gdk_window_set_cursor(plot->window, cursor);
	    gdk_cursor_unref(cursor);
	}
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

static int set_encoding_line (const char *s)
{
    return !strncmp(s, "set encod", 9);
}

static int set_print_line (const char *s)
{
    return (!strncmp(s, "set print ", 10) ||
	    !strncmp(s, "print \"pixe", 11) ||
	    !strncmp(s, "print \"data", 11) ||
	    !strncmp(s, "print \"term", 11));
}

enum {
    REMOVE_PNG,
    ADD_PNG
};

static int is_png_term_line (const char *s)
{
    return !strncmp(s, "set term png", 12);
}

static void real_remove_png_term_line (FILE *ftarg, FILE *fsrc,
				       GPT_SPEC *spec)
{
    char line[MAXLEN];
    int png_line_saved = 0;
    int printit;

    while (fgets(line, sizeof line, fsrc)) {
	printit = 1;
	if (is_png_term_line(line)) {
	    if (!png_line_saved) {
		/* comment it out, for future reference */
		fprintf(ftarg, "# %s", line);
		png_line_saved = 1;
	    }
	    printit = 0;
	} else if (commented_term_line(line)) {
	    if (png_line_saved) {
		printit = 0;
	    }
	} else if (set_output_line(line)) {
	    printit = 0;
	} else if (spec != NULL && (spec->flags & GPT_FIT_HIDDEN)
		   && is_auto_fit_string(line)) {
	    printit = 0;
	} else if (set_print_line(line)) {
	    printit = 0;
	}
	if (printit) {
	    fputs(line, ftarg);
	}
    }
}

static void add_new_term_line (FILE *fp, GPT_SPEC *spec,
			       int ix, int iy, GptFlags flags)
{
    const char *term_line;

    if (ix > 0 && iy > 0) {
	set_special_plot_size(ix, iy);
    }

    if (spec != NULL) {
	term_line = get_png_line_for_plotspec(spec);
    } else {
	term_line = gretl_gnuplot_term_line(GP_TERM_PNG,
					    PLOT_REGULAR,
					    flags, NULL);
    }

    fprintf(fp, "%s\n", term_line);
    if (strstr(term_line, "encoding") == NULL) {
	fputs("set encoding utf8\n", fp);
    }
}

static void real_add_png_term_line (FILE *ftarg, FILE *fsrc,
				    GPT_SPEC *spec)
{
    char line[MAXLEN];
    gchar *restore = NULL;
    GptFlags flags = 0;
    int add_line_styles = 1;
    int ix = 0, iy = 0;

    /* First pass: see if we can find an existing PNG term
       specification (possible commented out) to restore.
       Also see if we need to add line styles.
    */

    while (fgets(line, sizeof line, fsrc)) {
	if (restore == NULL && is_png_term_line(line)) {
	    restore = g_strdup(line);
	} else if (restore == NULL && commented_term_line(line)) {
	    restore = g_strdup(line + 2);
	} else if (strstr(line, "letterbox")) {
	    flags = GPT_LETTERBOX;
	} else if (strstr(line, "large")) {
	    flags = GPT_XL;
	} else if (strstr(line, "extra-large")) {
	    flags = GPT_XXL;
	} else if (strstr(line, "extra-wide")) {
	    flags = GPT_XW;
	} else if (sscanf(line, "# geoplot %d %d", &ix, &iy) == 2) {
	    ; /* OK */
	} else if (!strncmp(line, "set linetype", 12)) {
	    /* line styles already present */
	    add_line_styles = 0;
	} else if (!strncmp(line, "plot", 4)) {
	    break;
	}
    }

    /* prepare for second pass */
    rewind(fsrc);

    if (restore != NULL) {
	/* we got an existing term line: restore it */
	fputs(restore, ftarg);
	fputs("set encoding utf8\n", ftarg);
	g_free(restore);
    } else {
	/* we didn't get a previous term line */
	add_new_term_line(ftarg, spec, ix, iy, flags);
    }

    /* Note: write_plot_output_line() must come after "set encoding"
       since the output filename may be UTF-8.
    */
    write_plot_output_line(NULL, ftarg);

    if (spec != NULL && add_line_styles) {
	write_plot_line_styles(spec->code, ftarg);
    }

    /* Now for the rest of the plot file: skip elements that are
       already handled, otherwise transcribe from @fsrc.
    */
    while (fgets(line, sizeof line, fsrc)) {
	if (set_print_line(line)) {
	    ; /* skip it (portability) */
	} else if (is_png_term_line(line)) {
	    ; /* handled above */
	} else if (commented_term_line(line)) {
	    ; /* handled above */
	} else if (set_encoding_line(line)) {
	    ; /* handled above */
	} else if (set_output_line(line)) {
	    ; /* handled above */
	} else {
	    fputs(line, ftarg);
	}
    }

    write_plot_bounding_box_request(ftarg);
}

static int
add_or_remove_png_term (const char *fname, int action, GPT_SPEC *spec)
{
    FILE *ftarg = NULL;
    FILE *fsrc = NULL;
    gchar *temp = NULL;
    int err = 0;

    /* open stream to receive the revised content */
    temp = gretl_make_dotpath("gpttmp.XXXXXX");
    ftarg = gretl_tempfile_open(temp);
    if (ftarg == NULL) {
	g_free(temp);
	return 1;
    }

    /* open original file for reading */
    fsrc = open_gp_file(fname, "r");
    if (fsrc == NULL) {
	fclose(ftarg);
	g_free(temp);
	return 1;
    }

#if WINDEBUG
    fprintf(stderr, "add_or_remove_png_term\n src='%s'\n temp='%s'\n",
	    fname, temp);
#endif

    if (action == ADD_PNG) {
	real_add_png_term_line(ftarg, fsrc, spec);
    } else {
	real_remove_png_term_line(ftarg, fsrc, spec);
    }

    fclose(fsrc);
    fclose(ftarg);

    /* delete the original file */
    gretl_remove(fname);

    /* rename the revised version to the original name */
    err = gretl_rename(temp, fname);
    if (err) {
	fprintf(stderr, "warning: rename failed\n");
    }

    g_free(temp);

    return err;
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
	png_plot *plot = (png_plot *) p;

	spec = plot->spec;
    } else {
	/* EPS/PDF saver */
	spec = graph_saver_get_plotspec(p);
    }

    return spec->termtype;
}

static void emit_alt_datafile_line (const char *inpname,
				    FILE *fp)
{
    gchar *tmp = g_strdup_printf("%s.dat", inpname);

#ifdef G_OS_WIN32
    gchar *s = tmp;

    while (*s) {
	if (*s == '\\') *s = '/';
	s++;
    }
#endif

    fprintf(fp, "datafile = \"%s\"\n", tmp);
    g_free(tmp);
}

void filter_gnuplot_file (int mono, const char *inpname,
			  FILE *fpin, FILE *fpout)
{
    char pline[512];

    while (fgets(pline, sizeof pline, fpin)) {
	if (set_print_line(pline)) {
	    break;
	}

	if (!strncmp(pline, "set term", 8) ||
	    !strncmp(pline, "set enco", 8) ||
	    !strncmp(pline, "set outp", 8)) {
	    continue;
	}

	if (!strncmp(pline, "datafile = sprintf", 17)) {
	    /* geomap special */
	    emit_alt_datafile_line(inpname, fpout);
	    continue;
	}

	if (mono) {
	    if (strstr(pline, "set linetype") && strstr(pline, "rgb")) {
		continue;
	    } else if (strstr(pline, "set style fill solid")) {
		fputs("set style fill solid 0.3\n", fpout);
		continue;
	    }
	}

	fputs(pline, fpout);
    }
}

/* note @termstr will be non-NULL only if we're being
   called from the interactive PDF/EPS preview or save
   context
*/

static int revise_plot_file (GPT_SPEC *spec,
			     const char *inpname,
			     const char *outname,
			     const char *termstr)
{
    FILE *fpin = NULL;
    FILE *fpout = NULL;
    int mono = (spec->flags & GPT_MONO);
    int err = 0;

    fpin = gretl_fopen(spec->fname, "rb");
    if (fpin == NULL) {
	file_read_errbox(spec->fname);
	return 1;
    }

    fpout = gretl_fopen(inpname, "wb");
    if (fpout == NULL) {
	fclose(fpin);
	file_write_errbox(inpname);
	return 1;
    }

    if (outname != NULL && *outname != '\0') {
	if (termstr == NULL) {
	    termstr = gretl_gnuplot_term_line(spec->termtype,
					      spec->code,
					      spec->flags,
					      spec->fontstr);
	}
	fprintf(fpout, "%s\n", termstr);
	if (mono) {
	    fputs("set mono\n", fpout);
	}
	if (strstr(termstr, "set encoding") == NULL) {
	    fputs("set encoding utf8\n", fpout);
	}
	write_plot_output_line(outname, fpout);
    }

    filter_gnuplot_file(mono, spec->fname, fpin, fpout);

    fclose(fpin);
    fclose(fpout);

    return err;
}

void save_graph_to_file (gpointer data, const char *fname)
{
    png_plot *plot = (png_plot *) data;
    GPT_SPEC *spec = plot->spec;
    char pltname[FILENAME_MAX];
    int err = 0;

    sprintf(pltname, "%sgptout.tmp", gretl_dotdir());

    if (plot_is_geomap(plot)) {
	set_special_plot_size(plot->pixel_width, plot->pixel_height);
    }
    err = revise_plot_file(spec, pltname, fname, NULL);

    if (!err) {
	err = gnuplot_make_image(pltname);
	if (err) {
	    gui_errmsg(err);
	}
	gretl_remove(pltname);
    }
}

static gchar *map_pdf_termstr (png_plot *plot)
{
    double w = plot->pixel_width;
    double h = plot->pixel_height;
    double hr = h / w;

    if (hr > 8.5 / 5.5) {
	/* tall and skinny */
	h = 8.5;
	w = h / hr;
    } else {
	w = 5.5;
	h = w * hr;
    }

    return g_strdup_printf("set term pdfcairo font \"sans,12\" size %g,%g", w, h);
}

#define GRETL_PDF_TMP "gretltmp.pdf"
#define GRETL_EPS_TMP "gretltmp.eps"

static void graph_display_pdf (png_plot *plot)
{
    GPT_SPEC *spec = plot->spec;
    char pdfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    int err = 0;

    spec->termtype = GP_TERM_PDF;
    gretl_build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);
    gretl_build_path(pdfname, gretl_dotdir(), GRETL_PDF_TMP, NULL);

    if (plot_is_geomap(plot)) {
	gchar *termstr = map_pdf_termstr(plot);

	err = revise_plot_file(spec, plttmp, pdfname, termstr);
	g_free(termstr);
    } else {
	err = revise_plot_file(spec, plttmp, pdfname, NULL);
    }
    if (err) {
	return;
    }

    err = gnuplot_make_image(plttmp);
    gretl_remove(plttmp);

    if (err) {
	gui_errmsg(err);
	return;
    }

#if defined(G_OS_WIN32)
    win32_open_file(pdfname);
#elif defined(OS_OSX)
    osx_open_file(pdfname);
#else
    gretl_fork("viewpdf", pdfname, NULL);
#endif
}

void saver_preview_graph (GPT_SPEC *spec, char *termstr)
{
    char grfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    int err = 0;

    gretl_build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);

    if (spec->termtype == GP_TERM_EPS) {
	gretl_build_path(grfname, gretl_dotdir(), GRETL_EPS_TMP, NULL);
    } else {
	gretl_build_path(grfname, gretl_dotdir(), GRETL_PDF_TMP, NULL);
    }

    err = revise_plot_file(spec, plttmp, grfname, termstr);
    if (err) {
	return;
    }

    err = gnuplot_make_image(plttmp);
    gretl_remove(plttmp);

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
	gretl_fork("viewps", grfname, NULL);
    } else {
	gretl_fork("viewpdf", grfname, NULL);
    }
#endif
}

int saver_save_graph (GPT_SPEC *spec, char *termstr, const char *fname)
{
    char plttmp[FILENAME_MAX];
    int err;

    gretl_build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);
    err = revise_plot_file(spec, plttmp, fname, termstr);

    if (!err) {
	err = gnuplot_make_image(plttmp);
	gretl_remove(plttmp);
	if (err) {
	    gui_errmsg(err);
	}
    }

    return err;
}

/* we're looking for an uncommented "set term ..." */

static int is_term_line (const char *s, int *batch)
{
    int ret = 0;

    while (isspace(*s)) s++;

    if (*s != '#') {
	s = strstr(s, "set term");
	if (s != NULL) {
	    *batch = ret = 1;
	    s += 8;
	    s += strcspn(s, " "); /* skip "inal"? */
	    s += strspn(s, " ");  /* skip space */
	    if (!strncmp(s, "win", 3) ||
		!strncmp(s, "x11", 3) ||
		!strncmp(s, "wxt", 3) ||
		!strncmp(s, "qt", 2)) {
		/* these are all interactive */
		*batch = 0;
	    }
	}
    }

    return ret;
}

/* Check whether (a) we might want to prepend a "set
   term" statement (not if the user has already given
   one) and (b) whether we might want to append
   "pause mouse close" (not if a batch-type terminal
   has been specified).
*/

static void pre_test_plot_buffer (const char *buf,
				  int *addpause,
				  int *putterm,
				  int *dim)
{
    char bufline[512];
    int batch = 0;

    bufgets_init(buf);

    while (bufgets(bufline, sizeof bufline, buf)) {
	if (is_term_line(bufline, &batch)) {
	    *putterm = 0;
	    if (*addpause && batch) {
		*addpause = 0;
	    }
	} else if (!strncmp(bufline, "# geoplot", 9)) {
	    int w, h;

	    if (sscanf(bufline + 10, "%d %d", &w, &h) == 2) {
		dim[0] = w; dim[1] = h;
	    }
	} else if (*addpause && strstr(bufline, "pause ")) {
	    *addpause = 0;
	}
    }

    bufgets_finalize(buf);
}

static void geoplot_dump_revise (const char *buf,
				 const char *src,
				 FILE *fp)
{
    char bufline[512];

    bufgets_init(buf);
    while (bufgets(bufline, sizeof bufline, buf)) {
	if (!strncmp(bufline, "datafile = sp", 13)) {
	    emit_alt_datafile_line(src, fp);
	} else {
	    fputs(bufline, fp);
	}
    }
    bufgets_finalize(buf);
}

/* dump_plot_buffer: this is used when we're taking the material from
   an editor window containing gnuplot commands, and either (a)
   sending it to gnuplot for execution, or (b) saving it to a "user
   file".  In the saving-to-file case @addpause will be 0.

   This function handles the addition of "pause mouse close",
   if @addpause is non-zero.
*/

int dump_plot_buffer (const char *buf, const char *fname,
		      int addpause, const char *src)
{
    FILE *fp = gretl_fopen(fname, "wb");
    int dim[2] = {0};
    int putterm = addpause;
    int wxt_ok = 0;

    if (fp == NULL) {
	file_write_errbox(fname);
	return E_FOPEN;
    }

    if (addpause) {
	wxt_ok = gnuplot_has_wxt();
	pre_test_plot_buffer(buf, &addpause, &putterm, dim);
	if (putterm && wxt_ok) {
	    if (dim[0] > 0 && dim[1] > 0) {
		fprintf(fp, "set term wxt size %d,%d noenhanced\n",
			dim[0], dim[1]);
	    } else {
		fputs("set term wxt size 640,420 noenhanced\n", fp);
	    }
	}
    }

    if (src != NULL && strstr(buf, "datafile = sp")) {
	geoplot_dump_revise(buf, src, fp);
    } else {
	fputs(buf, fp);
    }

    if (addpause) {
	fputs("pause mouse close\n", fp);
    }

    fclose(fp);

    return 0;
}

#ifdef G_OS_WIN32

static int real_send_to_gp (const char *fname, int persist)
{
    return win32_run_async(gretl_gnuplot_path(), fname);
}

#else

static void gnuplot_done (GPid pid, gint status, gpointer p)
{
    if (p != NULL) {
	gint err_fd = GPOINTER_TO_INT(p);

	if (err_fd > 0) {
	    if (status != 0) {
		char buf[128] = {0};

		if (read(err_fd, buf, 127) > 0) {
		    errbox(buf);
		}
	    }
	    close(err_fd);
	}
    }

    g_spawn_close_pid(pid);
}

static int real_send_to_gp (const char *fname, int persist)
{
    const char *gp = gretl_gnuplot_path();
    GError *error = NULL;
    gchar *argv[4];
    GPid pid = 0;
    gint fd = -1;
    gboolean run;
    int err = 0;

    argv[0] = g_strdup(gp);
    argv[1] = g_strdup(fname);
    argv[2] = persist ? g_strdup("-persist") : NULL;
    argv[3] = NULL;

    run = g_spawn_async_with_pipes(NULL, argv, NULL,
				   G_SPAWN_SEARCH_PATH |
				   G_SPAWN_DO_NOT_REAP_CHILD,
				   NULL, NULL, &pid,
				   NULL, NULL,
				   &fd, &error);

    if (error != NULL) {
	errbox(error->message);
	g_error_free(error);
	err = 1;
    } else if (!run) {
	errbox(_("gnuplot command failed"));
	err = 1;
    } else if (pid > 0) {
	gpointer p = fd > 0 ? GINT_TO_POINTER(fd) : NULL;

	g_child_watch_add(pid, gnuplot_done, p);
    }

    g_free(argv[0]);
    g_free(argv[1]);
    g_free(argv[2]);

    return err;
}

#endif /* Windows vs other */

/* Callback for execute icon in window editing gnuplot
   commands: send script in @buf to gnuplot itself
*/

void run_gnuplot_script (gchar *buf, windata_t *vwin)
{
    gchar *src = NULL;
    gchar *tmpfile;
    int err;

    if (vwin != NULL && vwin->data != NULL) {
	src = session_graph_get_filename(vwin->data);
    }

    tmpfile = gretl_make_dotpath("showtmp.gp");
    err = dump_plot_buffer(buf, tmpfile, 1, src);

    if (!err) {
	err = real_send_to_gp(tmpfile, 1);
    }

    g_free(tmpfile);
    g_free(src);
}

#ifdef G_OS_WIN32

/* common code for sending an EMF file to the clipboard,
   or printing an EMF, on MS Windows
*/

static void win32_process_graph (png_plot *plot, int dest)
{
    GPT_SPEC *spec = plot->spec;
    char emfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    int err = 0;

    spec->termtype = GP_TERM_EMF;
    gretl_build_path(plttmp, gretl_dotdir(), "gptout.tmp", NULL);
    gretl_build_path(emfname, gretl_dotdir(), "gpttmp.emf", NULL);

    if (plot_is_geomap(plot)) {
	set_special_plot_size(plot->pixel_width, plot->pixel_height);
    }

    err = revise_plot_file(spec, plttmp, emfname, NULL);
    if (err) {
	return;
    }

    err = gnuplot_make_image(plttmp);
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

#else /* ! MS Windows */

static png_plot *copyplot;
static gboolean cb_image_mono;

/* Here we're just posting the information that an image
   should be available for pasting (and also recording
   if the user wanted it to be monochrome).
*/

static void set_plot_for_copy (png_plot *plot)
{
    copyplot = plot;
    cb_image_mono = plot->spec->flags & GPT_MONO ? 1 : 0;
    flag_image_available();
}

/* Here we're responding to a request to paste the plot
   advertised above; @target tells us which of the posted
   formats the application has selected.
*/

int write_plot_for_copy (int target)
{
    GPT_SPEC *spec = copyplot->spec;
    char outname[FILENAME_MAX];
    char inpname[FILENAME_MAX];
    GptFlags saveflags;
    double savescale;
    int saveterm;
    int err = 0;

    if (spec == NULL) {
	fprintf(stderr, "retrieve_plot_for_copy: no data\n");
	return 1;
    }

    saveflags = spec->flags;
    saveterm = spec->termtype;
    savescale = spec->scale;

    /* FIXME dimensions for PLOT_GEOMAP */

    if (target == TARGET_SVG) {
	spec->termtype = GP_TERM_SVG;
    } else if (target == TARGET_HTM) {
	spec->termtype = GP_TERM_HTM;
    } else if (target == TARGET_EMF) {
	spec->termtype = GP_TERM_EMF;
    } else if (target == TARGET_PNG) {
	spec->termtype = GP_TERM_PNG;
	spec->scale = 0.8;
    } else if (target == TARGET_EPS) {
	spec->termtype = GP_TERM_EPS;
    } else if (target == TARGET_PDF) {
	spec->termtype = GP_TERM_PDF;
    } else {
	fprintf(stderr, "write_plot_for_copy: unsupported type\n");
	return 1;
    }

    if (cb_image_mono) {
	spec->flags |= GPT_MONO;
    }

    gretl_build_path(inpname, gretl_dotdir(), "gptinp.tmp", NULL);
    gretl_build_path(outname, gretl_dotdir(), "gptout.tmp", NULL);
    err = revise_plot_file(spec, inpname, outname, NULL);

    spec->flags = saveflags;
    spec->termtype = saveterm;
    spec->scale = savescale;

    if (!err) {
	err = gnuplot_make_image(inpname);
	gretl_remove(inpname);
    }

    if (err) {
        errbox(_("Gnuplot error creating graph"));
    } else {
	err = image_file_to_clipboard(outname);
    }

    gretl_remove(outname);

    return err;
}

#endif /* Windows vs not */

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

/* returns non-zero on obtaining a marker */

static int get_gpt_marker (const char *line, char *label,
			   const char *format)
{
    const char *p = strchr(line, '#');

#if GPDEBUG > 1
    fprintf(stderr, "get_gpt_marker, p='%s'\n", p);
#endif

    *label = '\0';

    if (p != NULL) {
	sscanf(p + 2, format, label);
#if GPDEBUG > 1
	fprintf(stderr, "read marker: '%s'\n", label);
#endif
    }

    return *label != '\0';
}

/* special graphs for which editing via GUI is not supported */

#define cant_edit(p) (p == PLOT_CORRELOGRAM || \
                      p == PLOT_LEVERAGE || \
                      p == PLOT_MULTI_IRF || \
                      p == PLOT_MULTI_BASIC || \
                      p == PLOT_PANEL || \
                      p == PLOT_TRI_GRAPH || \
                      p == PLOT_BI_GRAPH || \
		      p == PLOT_ELLIPSE || \
		      p == PLOT_3D || \
		      p == PLOT_HEATMAP || \
		      p == PLOT_GEOMAP || \
		      p == PLOT_GRIDPLOT)

/* graphs where we don't attempt to find data coordinates */

#define no_readback(p) (p == PLOT_CORRELOGRAM || \
                        p == PLOT_LEVERAGE || \
                        p == PLOT_MULTI_IRF || \
                        p == PLOT_MULTI_BASIC || \
                        p == PLOT_PANEL || \
                        p == PLOT_TRI_GRAPH || \
                        p == PLOT_BI_GRAPH || \
			p == PLOT_STACKED_BAR || \
			p == PLOT_BAR || \
			p == PLOT_3D)

#define gp_missing(s) (s[0] == '?' || !strcmp(s, "NaN"))

static int get_gpt_data (GPT_SPEC *spec,
			 int *do_markers,
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

    if (*do_markers) {
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
	    x[0] = spec->data->val;
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
		if (gp_missing(test[j])) {
		    x[j][t] = NADBL;
		    missing++;
		} else {
		    x[j][t] = atof(test[j]);
		}
	    }

	    if (missing) {
		okobs--;
	    }

	    if (i <= imin && *do_markers) {
		*do_markers = get_gpt_marker(s, spec->markers[t], obsfmt);
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
    gretl_matrix_print(spec->data, "spec->data");
#endif

    gretl_pop_c_numeric_locale();

    return err;
}

static int get_gpt_heredata (GPT_SPEC *spec,
			     int *do_markers,
			     long barpos,
			     long datapos,
			     const char *buf)
{
    gretl_matrix *m = spec->data;
    char line[MAXLEN], test[32];
    char obsfmt[12] = {0};
    char *s, *got = NULL;
    double xij;
    int i, j;
    int err = 0;

    spec->okobs = spec->nobs;

    gretl_push_c_numeric_locale();

    /* first handle "shaded bars" info, if present */

    if (barpos > 0) {
	bufseek(buf, barpos);
	for (i=0; i<spec->nbars && !err; i++) {
	    double y1, y2, dx[2];

	    /* start date */
	    got = bufgets(line, sizeof line, buf);
	    if (got == NULL ||
		sscanf(line, "%lf %lf %lf", &dx[0], &y1, &y2) != 3) {
		err = 1;
		break;
	    }
	    /* stop date */
	    got = bufgets(line, sizeof line, buf);
	    if (got == NULL ||
		sscanf(line, "%lf %lf %lf", &dx[1], &y1, &y2) != 3) {
		err = 1;
		break;
	    }
	    plotspec_set_bar_info(spec, i, dx[0], dx[1]);
	}
    }

    if (*do_markers) {
	sprintf(obsfmt, "%%%d[^\r\n]", OBSLEN - 1);
    }

    /* then get the regular plot data */

    bufseek(buf, datapos);

    for (i=0; i<m->rows && !err; i++) {
	got = bufgets(line, sizeof line, buf);
	if (got == NULL) {
	    err = 1;
	}
	s = line;
	for (j=0; j<m->cols && !err; j++) {
	    s += strspn(s, " ");
	    sscanf(s, "%31s", test);
	    if (gp_missing(test)) {
		xij = NADBL;
	    } else {
		xij = atof(test);
	    }
	    gretl_matrix_set(m, i, j, xij);
	    s += strlen(test);
	    if (*do_markers) {
		*do_markers = get_gpt_marker(s, spec->markers[i], obsfmt);
	    }
	}
    }

    gretl_pop_c_numeric_locale();

#if GPDEBUG
    gretl_matrix_print(m, "gp heredata");
#endif

    gretl_pop_c_numeric_locale();

    return err;
}

static int line_starts_heredata (const char *line,
				 char **eod)
{
    int ret = 0;

    if (*line == '$') {
	const char *p = strstr(line, "<<");

	if (p != NULL) {
	    int len;

	    p += 2;
	    p += strspn(p, " ");
	    len = gretl_namechar_spn(p);
	    if (len > 0) {
		*eod = gretl_strndup(p, len);
		ret = 1;
	    }
	}
    }

    return ret;
}

/* This recognizes the limitations of the GPT_LABEL struct (see
   lib/src/plotspec.h). In principle that struct could be rejigged
   to hold information on "level" (front/back), point type and point
   size.
*/

static int unhandled_label_spec (const char *s)
{
    if (strstr(s, " front") || strstr(s, " back")) {
	return 1;
    } else if (strstr(s, " point") || strstr(s, " pt") ||
	       strstr(s, "ps")) {
	return 1;
    } else {
	return 0;
    }
}

/* read a gnuplot source line specifying a text label */

static int parse_label_line (GPT_SPEC *spec, const char *line,
			     int *unhandled)
{
    const char *p, *s;
    char *text = NULL;
    double x, y;
    int nc, just = GP_JUST_LEFT;
    int err = 0;

    /* Examples:
       (a) set label "this is a label" at 1998.26,937.557 left front
       (b) set label 'foobar' at 1500,350 left
       (c) set label "" at 0.0127,0.6308 front center point pt 8 ps 2
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

    /* note: text shouldn't be NULL but it may be empty */
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

    /* We need to watch out for case like (c) above, where the
       specification exceeds what we're currently able to pack into
       a GPT_LABEL struct in @spec.
    */
    if (!err && unhandled != NULL && unhandled_label_spec(p)) {
	*unhandled = 1;
	return 1;
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

    if (!err) {
        gretl_push_c_numeric_locale();
        if (!strcmp(s, "[*:*]")) {
            r0 = r1 = NADBL;
        } else if (strstr(s, "[*:")) {
            r0 = NADBL;
            r1 = atof(s + 3);
        } else if (strstr(s, ":*]")) {
            r0 = atof(s + 1);
            r1 = NADBL;
        } else if (sscanf(s, "[%lf:%lf]", &r0, &r1) != 2) {
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
    char fmt[16];
    char axis = '\0';
    int n, err = 0;

    if (timefmt) {
        if (strcmp(s, "%s") != 0) {
            gretl_errmsg_set("Invalid 'timefmt' setting");
            err = 1;
        } else {
            spec->flags |= GPT_TIMEFMT;
        }
    } else {
	n = sscanf(s, "%c \"%15[^\"]", &axis, fmt);
	err = n < 2;
    }

    if (!err) {
	if (axis == 'x') {
	    *spec->xfmt = '\0';
	    strncat(spec->xfmt, fmt, 15);
	} else if (axis == 'y') {
	    *spec->yfmt = '\0';
	    strncat(spec->yfmt, fmt, 15);
	}
    }

    return err;
}

static int catch_value (char *targ, const char *src, int maxlen)
{
    int i, n, q = 0;
    int ret = 0;

    src += strspn(src, " \t\r\n");
    if (*src == '\'' || *src == '"') {
	q = *src;
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
	    if (i == 0 && targ[i] == q) {
		ret = 1; /* accept empty string as "value" */
	    }
	    targ[i] = '\0';
	}
    }

    return ret || (*targ != '\0');
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

static int unhandled_gp_line_error (const char *s)
{
    int n = strlen(s);

    if (string_is_blank(s)) {
	return 0;
    }

    if (s[n-1] == '\n') {
	char *tmp = gretl_strndup(s, n-1);

	fprintf(stderr, "unhandled gnuplot command: '%s'\n", tmp);
	free(tmp);
    } else {
	fprintf(stderr, "unhandled gnuplot command: '%s'\n", s);
    }

    return 1;
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
	unhandled_gp_line_error(s);
    }

    return err;
}

/* here we're appending a second specification (e.g. rotation)
   to a prior tics spec (e.g. "nomirror")
*/

static void xtics_concat (char *targ, const char *src, size_t sz)
{
    int n = strlen(targ);

    if (targ[n-1] != ' ') {
	gchar *tmp = g_strdup_printf(" %s", src);

	strncat(targ, tmp, sz - n);
	g_free(tmp);
    } else {
	strncat(targ, src, sz - n);
    }
}

static void read_xtics_setting (GPT_SPEC *spec,
				const char *key,
				const char *val)
{
    if (!strcmp(key, "x2tics")) {
	spec->x2ticstr = gretl_strdup(val);
    } else if (*val == '(') {
	/* the 'special' form of xtics */
	int len = strlen(val);

	if (val[len-1] == '\\') {
	    fprintf(stderr, "got continuing xtics\n");
	}
	spec->xticstr = gretl_strdup(val);
    } else {
	/* regular xtics specification: this might come in
	   two parts if rotation is specified
	*/
	size_t sz = sizeof(spec->xtics) - 1;

	if (*spec->xtics == '\0') {
	    strncat(spec->xtics, val, sz);
	} else if (!strcmp(spec->xtics, "none")) {
	    *spec->xtics = '\0';
	    strncat(spec->xtics, val, sz);
	} else {
	    /* concatenate */
	    xtics_concat(spec->xtics, val, sz);
	}
    }
}

/* Try to accept variant RGB specifications besides the
   one that's standard in gretl plot files, namely
   "#RRGGBB" in hex. The candidate string -- which might
   be a color name of up to 17 characters -- is in @rgb,
   and on success it's ovewritten by a standard gnuplot
   hex string.
*/

static int verify_rgb (char *rgb)
{
    char *s = rgb;
    char delim = s[0];
    int err = 0;

    if (delim == '"' || delim == '\'') {
	/* valid input will be quoted */
	const char *p;

	s++;
	p = strchr(s, delim);
	if (p != NULL) {
	    gretlRGB color;
	    char test[18] = {0};
	    int len = p - s;

	    if (len >= 3 && len <= 17) {
		strncat(test, s, len);
		color = numeric_color_from_string(test, &err);
		if (!err) {
		    sprintf(rgb, "#%x", color);
		}
	    } else {
		err = E_DATA;
	    }
	} else {
	    err = E_DATA;
	}
    }

    return err;
}

/* e.g. set linetype 1 pt 6  lc rgb "#1B9E77" # dark teal
        set linetype 2 pt 7  lc rgb "#D95F02" # dark orange
*/

static int parse_linetype (const char *s, linestyle *styles)
{
    const char *p;
    int n, i, idx = 0;
    int err = 0;

    /* read the 1-based index of the linetype */
    n = sscanf(s, " %d", &idx);

    if (n != 1 || idx <= 0 || idx > N_GP_LINETYPES) {
	/* get out on error */
	return 1;
    } else {
	/* convert index to zero-based */
	i = idx - 1;
    }

    if (!err && (p = strstr(s, " lc ")) != NULL) {
	/* check for a color specification (20 bytes
	   allows for quoted named colors)
	*/
	char lc[20];

	p += 4;
	if (!strncmp(p, "rgb ", 4)) {
	    p += 4;
	}
	p += strspn(p, " ");
	if (sscanf(p, "%19s", lc)) {
	    /* validate and convert to hex if needed */
	    err = verify_rgb(lc);
	    if (!err) {
		strcpy(styles[i].lc, lc);
	    }
	} else {
	    err = 1;
	}
    }
    if (!err && (p = strstr(s, " lw ")) != NULL) {
	/* check line-width specification */
	float lw = 0.0;

	if (sscanf(p + 4, "%f", &lw)) {
	    styles[i].lw = lw;
	} else {
	    err = 1;
	}
    }
    if (!err && (p = strstr(s, " dt ")) != NULL) {
	/* check dash-type specification */
	int dt = 0;

	if (sscanf(p + 4, "%d", &dt)) {
	    styles[i].dt = dt;
	} else {
	    err = 1;
	}
    }
    if (!err && (p = strstr(s, " pt ")) != NULL) {
	/* check point-type specification */
	int pt = 0;

	if (sscanf(p + 4, "%d", &pt)) {
	    styles[i].pt = pt;
	} else {
	    err = 1;
	}
    }

#if LT_DEBUG
    fprintf(stderr, "styles[%d]: idx %d lc '%s' pt %d lw %g dt %d\n",
	    i, idx, styles[i].lc, styles[i].pt, (double) styles[i].lw,
	    styles[i].dt);
#endif

    return err;
}

static int parse_gp_set_line (GPT_SPEC *spec,
			      const char *s,
			      linestyle *styles,
			      int *unhandled)
{
    char *p, key[16] = {0};
    char val[MAXLEN] = {0};
    char *extra = NULL;
    int lt_pos = 0;
    int err = 0;

    if (!strncmp(s, "set linetype", 12)) {
	lt_pos = 12;
    } else if (!strncmp(s, "set style fill solid", 20)) {
	/* special with multi-word key */
	strncat(val, s + 20, 8);
	spec->fillfrac = (float) atof(val);
	val[0] = '\0';
    }

    if (lt_pos > 0) {
	err = parse_linetype(s + lt_pos, styles);
	if (err && unhandled != NULL) {
	    *unhandled = 1;
	}
	return err;
    }

    if (unhandled != NULL) {
	/* we're peeking at "literal" lines */
	if (strstr(s, "set style data histogram")) {
	    spec->code = PLOT_BAR;
	    *unhandled = 1;
	    return 0;
	} else if (strstr(s, "histogram rowstacked")) {
	    spec->code = PLOT_STACKED_BAR;
	    *unhandled = 1;
	    return 0;
	}
    }

    if (sscanf(s + 4, "%11s", key) != 1) {
	return unhandled_gp_line_error(s);
    }

    if ((p = strchr(key, '[')) != NULL) {
	/* try fixing this */
	*p = '\0';
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
    } else if (!strcmp(key, "mono") ||
	       !strcmp(key, "monochrome")) {
	spec->flags |= GPT_MONO;
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
	parse_label_line(spec, s, unhandled);
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
	    return unhandled_gp_line_error(s);
	}
    } else if (!strcmp(key, "logscale")) {
	if (read_plot_logscale(val, spec)) {
	    return unhandled_gp_line_error(s);
	}
    } else if (!strcmp(key, "format")) {
	if (read_plot_format(val, spec, 0)) {
	    return unhandled_gp_line_error(s);
	}
    } else if (!strcmp(key, "timefmt")) {
	if (read_plot_format(val, spec, 1)) {
	    return unhandled_gp_line_error(s);
	}
    } else if (!strcmp(key, "title")) {
	spec->titles[0] = g_strdup(val);
    } else if (!strcmp(key, "xlabel")) {
	spec->titles[1] = g_strdup(val);
	*spec->xvarname = '\0';
	strncat(spec->xvarname, val, MAXDISP-1);
    } else if (!strcmp(key, "ylabel")) {
	spec->titles[2] = g_strdup(val);
	*spec->yvarname = '\0';
	strncat(spec->yvarname, val, MAXDISP-1);
    } else if (!strcmp(key, "y2label")) {
	spec->titles[3] = g_strdup(val);
    } else if (!strcmp(key, "x2label")) {
	spec->titles[4] = g_strdup(val);
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
	if ((extra = strstr(val, "lc rgb \"")) != NULL) {
	    spec->border_lc[0] = '\0';
	    strncat(spec->border_lc, extra + 8, 7);
	}
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
    } else if (unhandled != NULL) {
	*unhandled = 1;
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

static void plotspec_destroy_markers (GPT_SPEC *spec)
{
    strings_array_free(spec->markers, spec->n_markers);
    spec->markers = NULL;
    spec->n_markers = 0;
}

static int get_boxplot_auxdata (GPT_SPEC *spec, char *line,
				const char *buf)
{
    double x, y;
    int i, n, r = 0, c = 0;
    int err = 0;

    n = sscanf(line + 9, "%d %d", &r, &c);

    if (n == 2 && r > 0 && c == 2) {
	spec->auxdata = gretl_matrix_alloc(r, c);
	if (spec->auxdata != NULL) {
	    gretl_push_c_numeric_locale();
	    if (spec->heredata) {
		/* skip a line: "$aux << EOA" */
		if (bufgets(line, MAXLEN - 1, buf) == NULL) {
		    err = E_DATA;
		}
	    }
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

    return err;
}

/* Determine the number of data points in a plot. While we're at it,
   determine the type of plot, check whether there are any
   data-point markers along with the data, and see if there are
   are definitions of time-series vertical bars.
*/

static void get_plot_nobs (png_plot *plot,
			   const char *buf,
			   int *do_markers,
			   long *barpos,
			   long *datapos)
{
    GPT_SPEC *spec = plot->spec;
    int nobs = 0, started = -1;
    int startmin = 1;
    long auxpos = 0;
    char line[MAXLEN], test[12];
    char *eod = NULL;
    char *p = NULL;

    spec->code = PLOT_REGULAR;
    *do_markers = 0;

    while (bufgets(line, MAXLEN - 1, buf)) {
	if (*line == '#' && spec->code == PLOT_REGULAR) {
	    tailstrip(line);
	    spec->code = plot_type_from_string(line);
	    if (spec->code == PLOT_GEOMAP) {
		int w, h;

		if (sscanf(line + 9, "%d %d", &w, &h) == 2) {
		    plot->pixel_width = w;
		    plot->pixel_height = h;
		}
		break;
	    }
	}
	if (sscanf(line, "# n_bars = %d", &spec->nbars) == 1) {
	    /* for old-style data representation */
	    startmin += spec->nbars;
	    continue;
	}
	if (spec->nbars > 0 && !strncmp(line, "$bars <", 7)) {
	    *barpos = buftell(buf);
	}
	if (!strncmp(line, "# auxdata", 9)) {
	    auxpos = buftell(buf);
	}
	if (started < 0 && line_starts_heredata(line, &eod)) {
	    /* "heredoc" method of handling plot data */
	    spec->heredata = 1;
	    *datapos = buftell(buf);
	    continue;
	}
	if (eod != NULL) {
	    if (!strncmp(line, eod, strlen(eod))) {
		free(eod);
		eod = NULL;
		break;
	    } else if (*do_markers == 0 && (p = strchr(line, '#')) != NULL) {
		if (sscanf(p + 1, "%8s", test) == 1) {
		    *do_markers = 1;
		}
	    }
	    nobs++;
	} else {
	    if (!strncmp(line, "plot", 4) || !strncmp(line, "splot", 5)) {
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
		if (*line == 'e' || !strncmp(line, "set ", 4)) {
		    /* end of data, or onto "set print" for bounds */
		    break;
		} else if (*do_markers == 0 && (p = strchr(line, '#')) != NULL) {
		    if (sscanf(p + 1, "%11s", test) == 1) {
			*do_markers = 1;
		    }
		}
		nobs++;
	    }
	}
    }

    spec->nobs = nobs;

    if (spec->code == PLOT_BOXPLOTS) {
	/* In the case of boxplots the gnuplot file may contain extra
	   (outlier) data, which have to be read into an auxiliary
	   matrix.
	*/
	if (auxpos > 0) {
	    /* we already went past it */
	    bufseek(buf, auxpos);
	    buf_back_lines(buf, 1);
	    bufgets(line, MAXLEN - 1, buf);
	} else if (!spec->heredata) {
	    /* try reading further */
	    while (bufgets(line, MAXLEN - 1, buf)) {
		if (!strncmp(line, "# auxdata", 9)) {
		    auxpos = buftell(buf);
		    break;
		}
	    }
	}
	if (auxpos > 0) {
	    get_boxplot_auxdata(spec, line, buf);
	}
    }
}

static int grab_fit_coeffs (GPT_SPEC *spec, const char *s)
{
    gretl_matrix *b = NULL;
    int k, f = spec->fit;
    int n = 0, err = 0;

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
    } else if (f == PLOT_FIT_LINLOG) {
	spec->b_linlog = b;
    }

    return err;
}

/* scan the stuff after "title '" or 'title "' */

static void grab_line_title (GPT_LINE *line, const char *src)
{
    char tmp[128];

    *tmp = '\0';

    if (*src == '\'') {
	sscanf(src + 1, "%127[^']'", tmp);
    } else {
	sscanf(src + 1, "%127[^\"]\"", tmp);
    }

    if (*tmp != '\0') {
	line->title = g_strdup(tmp);
    }
}

static void grab_line_rgb (char *targ, const char *src)
{
    if (*src == '"') {
	sscanf(src + 1, "%7s", targ);
    }
}

/* Examples:

   using 1:2
   using 1:n:m:p
   using 1:xtic(...)
   using 1:2:xtic(""):ytic("")
   using 1:($2+x)

   we need to extract the data column numbers, and if we
   find special stuff inlinerecord the 'using' string as
   a whole
*/

static int process_using_spec (const char **ps,
			       GPT_SPEC *spec,
			       int i)
{
    GPT_LINE *line = &spec->lines[i];
    const char *f, *p, *s = *ps;
    int *cols = NULL;
    int inparen = 0;
    int anyparen = 0;
    int ticspecs = 0;
    int n_uniq = 0;
    int k, n = 0;
    int err = 0;

    if (isdigit(*s)) {
	k = atoi(s);
	gretl_list_append_term(&cols, k);
	n_uniq++;
    }

    p = s;

    while (*s) {
	if (*s == '(') {
	    anyparen = 1;
	    inparen++;
	} else if (*s == ')') {
	    inparen--;
	} else if (!inparen) {
	    if (*s == ':') {
		f = s + 1;
		if (*f == 'x' || *f == 'y') {
		    ticspecs++;
		} else {
		    k = (*f == '$')? atoi(f + 1) : atoi(f);
		    if (k > 0) {
			if (!in_gretl_list(cols, k)) {
			    n_uniq++;
			}
			gretl_list_append_term(&cols, k);
		    }
		}
	    } else if (isspace(*s) || *s == ',') {
		break;
	    }
	} else if (*s == '$' && isdigit(*(s+1))) {
	    /* inparen: a tic specification may make reference
	       to a data column */
	    k = atoi(s + 1);
	    if (k > 0) {
		if (!in_gretl_list(cols, k)) {
		    n_uniq++;
		}
		gretl_list_append_term(&cols, k);
	    }
	}
	s++;
	n++;
    }

    /* pass back pointer to remainder of line, if any */
    *ps = s;

    line->style = GP_STYLE_AUTO;

    if (n > 0 && (anyparen || ticspecs > 0)) {
	/* record the 'using' string as is */
	char *ustr = gretl_strndup(p, n);
#if GPDEBUG
	fprintf(stderr, "saving ustr = '%s'\n", ustr);
#endif
	line->ustr = ustr;
    }

    if (cols != NULL) {
	line->ncols = n_uniq;
	if (cols[0] == 5 && n_uniq == 2) {
	    /* boxplot special */
	    line->flags |= GP_LINE_BOXDATA;
	}
	/* cumulate total data columns */
	if (spec->datacols > 0 && in_gretl_list(cols, 1)) {
	    /* col 1 should already be counted? */
	    spec->datacols += line->ncols - 1;
	} else {
	    spec->datacols += line->ncols;
	}
#if GPDEBUG
	fprintf(stderr, "number of unique columns %d\n", n_uniq);
	printlist(cols, "cols list");
	fprintf(stderr, "spec->datacols now = %d\n", spec->datacols);
#endif
	if (line->ustr == NULL) {
	    line->mcols = cols;
	} else {
	    free(cols);
	}
    }

    return err;
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

#if GPDEBUG
    fprintf(stderr, "parse_gp_line_line, starting\n");
#endif

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
	err = process_using_spec(&p, spec, i);
	s = p; /* remainder of line */
    } else if (*s == '\'' || *s == '"') {
	/* name of data file, without 'using': implicitly
	   using cols 1 and 2
	*/
	if (*(s+1) != '-') {
	    fprintf(stderr, "plotting datafile, not supported\n");
	} else {
	    line->ncols = 2;
	    spec->datacols += 2;
	}
    } else {
	/* absence of "using" should mean that the line plots
	   a formula, not a set of data columns
	*/
	/* get the formula: it runs up to "title" or "notitle" */
	p = strstr(s, " title");
	if (p == NULL) {
	    p = strstr(s, " notitle");
	}
	if (p != NULL) {
	    line->formula = g_strndup(s, p - s);
	    if (i == 1 && spec->flags & GPT_AUTO_FIT) {
		grab_fit_coeffs(spec, line->formula);
	    }
	} else {
	    /* title must be implicit? */
	    gchar *tmp = g_strstrip(g_strdup(s));

	    line->formula = tmp;
	    line->title = g_strdup(tmp);
	}
    }

    if (strstr(s, "axes x1y2")) {
	line->yaxis = 2;
    }
    if ((p = strstr(s, " title "))) {
	grab_line_title(line, p + 7);
    }
    if ((p = strstr(s, " w "))) {
	sscanf(p + 3, "%15[^, ]", tmp);
	line->style = gp_style_index_from_name(tmp);
    }
    if ((p = strstr(s, " lt "))) {
	sscanf(p + 4, "%d", &line->type);
    } else if ((p = strstr(s, " lc rgb "))) {
	grab_line_rgb(line->rgb, p + 8);
    }
    if ((p = strstr(s, " ps "))) {
	sscanf(p + 4, "%f", &line->pscale);
    }
    if ((p = strstr(s, " pt "))) {
	sscanf(p + 4, "%d", &line->ptype);
    }
    if ((p = strstr(s, " dt "))) {
	sscanf(p + 4, "%d", &line->dtype);
    }
    if (!auto_linewidth && (p = strstr(s, " lw "))) {
	sscanf(p + 4, "%f", &line->width);
    }
    if ((p = strstr(s, " whiskerbars "))) {
	double ww;

	sscanf(p + 13, "%lf", &ww);
	line->whiskwidth = (float) ww;
    }

    if (line->ncols == 0 && line->formula == NULL) {
	/* got neither data column spec nor formula */
	err = 1;
    }

#if GPDEBUG
    fprintf(stderr, "parse_gp_line_line, returning %d\n", err);
#endif

    return err;
}

/* We got a special comment supposedly indicating the name and ID
   number of the X or Y variable in an OLS fitted line. Here we check
   this info for validity: @vname should be the name of a bona fide
   series variable, and its ID number should match the given @v.  We
   return 1 if this works out OK, otherwise 0.
*/

static int plot_ols_var_ok (const char *vname, int v)
{
    int vcheck = current_series_index(dataset, vname);

    return vcheck == v;
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

static int plot_get_data_and_markers (GPT_SPEC *spec,
				      const char *buf,
				      int *do_markers,
				      long barpos,
				      long datapos)
{
    int err = 0;

    if (spec->nobs == 0 || spec->datacols == 0) {
	/* nothing to be done */
	return 0;
    }

#if GPDEBUG
    fprintf(stderr, "plot_get_data, allocating: nobs=%d, datacols=%d\n",
	    spec->nobs, spec->datacols);
#endif

    /* allocate for the plot data... */
    spec->data = gretl_matrix_alloc(spec->nobs, spec->datacols);
    if (spec->data == NULL) {
	err = E_ALLOC;
    }

    /* and markers if any */
    if (!err && *do_markers) {
	err = plotspec_allocate_markers(spec);
    }

    /* and time-series bars, if any */
    if (!err && spec->nbars > 0) {
	err = plotspec_allocate_bars(spec);
    }

    /* read the data (and perhaps markers) from the plot file */
    if (!err) {
	if (spec->heredata) {
	    err = get_gpt_heredata(spec, do_markers,
				   barpos, datapos, buf);
	} else {
	    err = get_gpt_data(spec, do_markers, buf);
	}
    } else {
	gretl_matrix_free(spec->data);
	spec->data = NULL;
    }

    if (spec->markers != NULL && *do_markers == 0) {
	/* something must have gone wrong: clean up */
	plotspec_destroy_markers(spec);
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
    } else if (strstr(s, "linlog")) {
	return PLOT_FIT_LINLOG;
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
	} else if (strstr(line, "extra-wide")) {
	    spec->flags |= GPT_XW;
	    break;
	}
	i++;
    }
}

static void check_for_plot_time_data (GPT_SPEC *spec, gchar *buf)
{
    char line[128];

    while (bufgets(line, sizeof line, buf)) {
	if (!strncmp(line, "set timefmt ", 12)) {
	    spec->flags |= GPT_TIMEFMT;
	} else if (*line == '#' && strstr(line, "letterbox")) {
	    spec->flags |= GPT_LETTERBOX;
	} else if (!strncmp(line, "plot", 4)) {
	    break;
	}
    }
}

static void linestyle_init (linestyle *ls)
{
    ls->lc[0] = '\0';
    ls->lw = 0;
    ls->dt = 0;
    ls->pt = 0;
}

static int push_z_row (gretl_matrix *z, int i, int n, char *line)
{
    char *p = line;
    double x;
    int j, err = 0;

    errno = 0;

    for (j=0; j<n && *p; j++) {
	if (*p == ' ') {
	    p++;
	}
	if (*p == '?') {
	    x = NADBL;
	    p += 2;
	} else if (!strncmp(p, "NaN", 3)) {
	    x = NADBL;
	    p += 4;
	} else {
	    x = strtod(p, &p);
	    if (errno) {
		err = 1;
		break;
	    }
	}
	gretl_matrix_set(z, i, j, x);
    }

    return err;
}

static void get_heatmap_matrix (GPT_SPEC *spec, gchar *buf,
				char *line, size_t len,
				long datapos)
{
    gretl_matrix *z;
    int n = spec->nobs;
    int i, err = 0;

    z = gretl_matrix_alloc(n, n);
    if (z == NULL) {
	return;
    }

    bufseek(buf, datapos);

    gretl_push_c_numeric_locale();

    for (i=0; i<n && !err; i++) {
	bufgets(line, len, buf);
	err = push_z_row(z, i, n, line);
    }

    gretl_pop_c_numeric_locale();

    if (err) {
	gretl_matrix_free(z);
    } else {
	spec->auxdata = z;
    }
}

/* Here we're seeing if we should "promote" literal plot lines to
   supplement or override "set" variables that are recognized by gretl
   and form part of @spec. If we don't do this the GUI plot editor may
   be broken.
*/

static void maybe_promote_literal_lines (GPT_SPEC *spec,
					 linestyle *styles)
{
    const char *line;
    int *rmlines = NULL;
    int unhandled;
    int i, err = 0;

    for (i=0; i<spec->n_literal; i++) {
	line = spec->literal[i];
	if (!strncmp(line, "set ", 4)) {
	    unhandled = 0;
	    gretl_push_c_numeric_locale();
	    err = parse_gp_set_line(spec, line, styles, &unhandled);
	    gretl_pop_c_numeric_locale();
	    if (!err && !unhandled) {
		gretl_list_append_term(&rmlines, i);
	    }
	} else if (!strncmp(line, "unset ", 6)) {
	    err = parse_gp_unset_line(spec, line);
	    if (!err) {
		gretl_list_append_term(&rmlines, i);
	    }
	}
    }

    if (rmlines != NULL) {
	/* at least one "literal" line has been promoted */
	if (rmlines[0] == spec->n_literal) {
	    strings_array_free(spec->literal, spec->n_literal);
	    spec->literal = NULL;
	    spec->n_literal = 0;
	} else {
	    char **S = spec->literal;
	    int ns, orig_ns = spec->n_literal;

	    ns = orig_ns - rmlines[0];
	    spec->literal = strings_array_new(ns);

	    if (spec->literal == NULL) {
		/* unlikely, but... */
		spec->literal = S;
	    } else {
		/* reduced version of the "literal" array */
		int j = 0;

		for (i=0; i<orig_ns; i++) {
		    if (in_gretl_list(rmlines, i)) {
			free(S[i]);
		    } else {
			spec->literal[j++] = S[i];
		    }
		    S[i] = NULL;
		}
		strings_array_free(S, orig_ns);
		spec->n_literal = ns;
	    }
	}
    }
}

static void plot_get_ols_info (const char *line,
			       int *reglist)
{
    char fmt[16], vname[VNAMELEN];
    int v;

    sprintf(fmt, "'%%%d[^\\']' (%%d)", VNAMELEN - 1);

    if (sscanf(line + 6, fmt, vname, &v) == 2) {
	if (line[2] == 'X') {
	    if (plot_ols_var_ok(vname, v)) {
		reglist[3] = v;
	    }
	} else {
	    /* 'Y' */
	    if (reglist[3] > 0 && plot_ols_var_ok(vname, v)) {
		reglist[0] = 3;
		reglist[1] = v;
	    }
	}
    }
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

static int read_plotspec_from_file (png_plot *plot)
{
    GPT_SPEC *spec = plot->spec;
    linestyle styles[N_GP_LINETYPES];
    int do_markers = 0;
    int auto_linewidth = 0;
    int reglist[4] = {0};
    int *uservec = NULL;
    long barpos = 0;
    long datapos = 0;
    char gpline[MAXLEN];
    gchar *buf = NULL;
    char *got = NULL;
#if HANDLE_HEREDATA
    char *eod = NULL;
#endif
    int ignore = 0;
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

    spec->datacols = 0;
    bufgets_init(buf);

    /* get the number of data-points, plot type, and check for
       observation markers and time-series bars */
    get_plot_nobs(plot, buf, &do_markers, &barpos, &datapos);

    if (spec->code == PLOT_GEOMAP) {
	/* FIXME */
	goto bailout;
    }

    if (spec->nobs == 0 && plot_needs_obs(spec->code)) {
	/* failed reading plot data */
#if GPDEBUG
	fprintf(stderr, " got spec->nobs = 0\n");
#endif
	err = 1;
	goto bailout;
    }

    if (do_markers && spec->nobs > MAX_MARKERS) {
	do_markers = 0;
    }

    if (cant_edit(spec->code)) {
	fprintf(stderr, "read_plotspec_from_file: plot is not editable\n");
	if (maybe_big_multiplot(spec->code)) {
	    buf_rewind(buf);
	    check_for_plot_size(spec, buf);
	} else if (spec->code == PLOT_BAND) {
	    buf_rewind(buf);
	    check_for_plot_time_data(spec, buf);
	} else if (spec->code == PLOT_HEATMAP) {
	    buf_rewind(buf);
	    get_heatmap_matrix(spec, buf, gpline, sizeof gpline, datapos);
	}
	goto bailout;
    } else {
	buf_rewind(buf);
	check_for_plot_size(spec, buf);
    }

    for (i=0; i<N_GP_LINETYPES; i++) {
	linestyle_init(&styles[i]);
    }

    buf_rewind(buf);

    /* get the preamble and "set" lines */

    while ((got = bufgets(gpline, sizeof gpline, buf)) && !err) {
#if GPDEBUG
	tailstrip(gpline);
	fprintf(stderr, "gpline: '%s'\n", gpline);
#endif
	if (!strncmp(gpline, "plot ", 5)) {
	    /* we're done with 'set' */
	    break;
	}
	if (ignore) {
	    if (!strncmp(gpline, "# end inline", 12)) {
		ignore = 0;
	    }
	} else if (!strncmp(gpline, "# start inline", 14)) {
	    ignore = 1;
	} else if (!strncmp(gpline, "# timeseries", 12)) {
	    if (sscanf(gpline, "# timeseries %d", &spec->pd)) {
		plot->pd = spec->pd;
	    }
	    spec->flags |= GPT_TS;
	    if (strstr(gpline, "letterbox")) {
		spec->flags |= GPT_LETTERBOX;
	    }
	} else if (!strncmp(gpline, "# multiple timeseries", 21)) {
	    if (sscanf(gpline, "# multiple timeseries %d", &spec->pd)) {
		plot->pd = spec->pd;
	    }
	    spec->flags |= GPT_TS;
	} else if (!strncmp(gpline, "# scale = ", 10)) {
	    double x = dot_atof(gpline + 10);

	    if (x >= 0.5 && x <= 2.0) {
		spec->scale = x;
	    }
	} else if (!strncmp(gpline, "# auto linewidth", 16)) {
	    auto_linewidth = 1;
	} else if (!strncmp(gpline, "# fontspec: ", 12)) {
	    free(spec->fontstr);
	    spec->fontstr = gretl_strdup(gpline + 12);
	    tailstrip(spec->fontstr);
	} else if (!strncmp(gpline, "# boxplots", 10)) {
	    ; /* OK */
	} else if (!strncmp(gpline, "# X = ", 6) ||
		   !strncmp(gpline, "# Y = ", 6)) {
	    plot_get_ols_info(gpline, reglist);
	} else if (sscanf(gpline, "# literal lines = %d", &spec->n_literal)) {
	    spec->literal = strings_array_new(spec->n_literal);
	    if (spec->literal == NULL) {
		err = E_ALLOC;
	    } else {
		for (i=0; i<spec->n_literal; i++) {
		    if (!bufgets(gpline, MAXLEN - 1, buf)) {
			errbox(_("Plot file is corrupted"));
		    } else {
			g_strstrip(gpline);
			spec->literal[i] = gretl_strdup(gpline);
		    }
		}
	    }
	} else if (sscanf(gpline, "# xtics lines = %d", &spec->n_xtics)) {
	    spec->multi_xtics = strings_array_new(spec->n_xtics);
	    if (spec->multi_xtics == NULL) {
		err = E_ALLOC;
	    } else {
		for (i=0; i<spec->n_xtics; i++) {
		    if (!bufgets(gpline, MAXLEN - 1, buf)) {
			errbox(_("Plot file is corrupted"));
		    } else {
			g_strstrip(gpline);
			spec->multi_xtics[i] = gretl_strdup(gpline);
		    }
		}
	    }
	} else if (!strncmp(gpline, "# start literal lines", 21) &&
		   spec->literal == NULL) {
	    for (i=0; ; i++) {
		if (!bufgets(gpline, MAXLEN - 1, buf)) {
		    errbox(_("Plot file is corrupted"));
		} else if (!strncmp(gpline, "# end literal lines", 19)) {
		    break;
		} else {
		    top_n_tail(gpline, 0, NULL);
		    strings_array_add(&spec->literal, &spec->n_literal,
				      gpline);
		}
	    }
	} else if (strstr(gpline, "automatic fit")) {
	    spec->flags |= GPT_AUTO_FIT;
	    spec->fit = recognize_fit_string(gpline);
	} else 	if (strstr(gpline, "user-defined")) {
	    uservec = get_user_lines_list(gpline);
	} else if (strstr(gpline, "printing data labels")) {
	    spec->flags |= GPT_PRINT_MARKERS;
	} else if (!strncmp(gpline, "# ", 2)) {
	    ; /* ignore unknown comment lines */
	} else if (!strncmp(gpline, "set ", 4)) {
	    gretl_push_c_numeric_locale();
	    err = parse_gp_set_line(spec, gpline, styles, NULL);
	    gretl_pop_c_numeric_locale();
	} else if (!strncmp(gpline, "unset ", 6)) {
	    err = parse_gp_unset_line(spec, gpline);
	} else {
	    /* done reading "set" lines? */
	    break;
	}
    }

    if (!err && got == NULL) {
	err = E_DATA;
    }

    if (err) {
	goto bailout;
    }

    for (i=0; i<4; i++) {
	if (!string_is_blank(spec->titles[i])) {
	    gretl_delchar('"', spec->titles[i]);
	}
    }

#if HANDLE_HEREDATA
    /* Hmm, what's this doing here? (2020-07-13) */
    if (line_starts_heredata(gpline, &eod)) {
	while (bufgets(gpline, MAXLEN - 1, buf) != NULL) {
	    if (!strncmp(gpline, eod, strlen(eod))) {
		bufgets(gpline, MAXLEN - 1, buf);
		break;
	    }
	}
	free(eod);
    }
#endif

    /* then get the "plot" lines */
    if (strncmp(gpline, "plot ", 5) ||
	(strlen(gpline) < 10 && bufgets(gpline, MAXLEN - 1, buf) == NULL)) {
	err = unhandled_gp_line_error(gpline);
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
	err = parse_gp_line_line(gpline, spec, auto_linewidth);
	if (err || done || (got = bufgets(gpline, MAXLEN - 1, buf)) == NULL) {
	    break;
	}
    }

    if (!err && got == NULL) {
	err = E_DATA;
    }

    if (!err && spec->literal != NULL) {
	maybe_promote_literal_lines(spec, styles);
    }

    /* transcribe styles info into lines for use in the
       GUI plot-editor */

    for (i=0; i<spec->n_lines && !err; i++) {
	GPT_LINE *line = &spec->lines[i];
	int idx = line->type; /* this will be 1-based */

	if (idx == LT_AUTO) {
	    /* automatic sequential line-types */
	    idx = i + 1;
	}
	if (uservec != NULL && in_gretl_list(uservec, i)) {
	    line->flags |= GP_LINE_USER;
	}
	if (idx > 0 && idx < N_GP_LINETYPES) {
	    int j = idx - 1;

	    if (line->rgb[0] == '\0') {
		strcpy(line->rgb, styles[j].lc);
	    }
	    if (line->width == 1 && styles[j].lw > 0) {
		line->width = styles[j].lw;
	    }
	    if (line->dtype == 0 && styles[j].dt > 0) {
		line->dtype = styles[j].dt;
	    }
	    if (line->ptype == 0 && styles[j].pt > 0) {
		line->ptype = styles[j].pt;
	    }
	}
#if LT_DEBUG
	fprintf(stderr, "spec->lines[%d]: type %d idx %d rgb '%s' pt %d lw %g\n",
		i, line->type, idx, line->rgb, line->ptype, (double) line->width);
#endif
	if (spec->auxdata != NULL && i == spec->n_lines - 1) {
	    /* the last "line" doesn't use the regular
	       data mechanism */
	    line->flags = GP_LINE_AUXDATA;
	    continue;
	}
    }

#if GPDEBUG
    fprintf(stderr, "plotspec, after line transcription, err=%d\n", err);
#endif

    if (!err) {
	err = plot_get_data_and_markers(spec, buf, &do_markers,
					barpos, datapos);
    }

    if (!err && reglist[0] > 0) {
	spec->reglist = gretl_list_copy(reglist);
    }

    if (!err && spec->fit == PLOT_FIT_NONE) {
	maybe_set_add_fit_ok(spec);
    }

    if (!err && spec->code == PLOT_ROOTS && spec->n_literal == 8) {
	/* backward compatibility */
	fprintf(stderr, "old roots plot: fixing\n");
	fix_old_roots_plot(spec);
    }

 bailout:

    if (err) {
	fprintf(stderr, "read_plotspec_from_file: err = %d\n", err);
    }

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

    if (na(dx) || na(dy)) {
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

    if (plot->cs == NULL) {
	fprintf(stderr, "copy_pixbuf_to_surface: failed\n");
	return 1;
    } else if (cairo_surface_status(plot->cs) != CAIRO_STATUS_SUCCESS) {
	fprintf(stderr, "copy_pixbuf_to_surface: failed\n");
	cairo_surface_destroy(plot->cs);
	plot->cs = NULL;
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
    gchar alt_label[16] = {0};
    PangoContext *context;
    PangoLayout *pl;

    if (plot_is_roots(plot)) {
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
    double dist, mindist = DBL_MAX;
    int best_match = -1;
    int done = 0, bank = 1;
    int t;

#if GPDEBUG > 2
    fprintf(stderr, "identify_point: pixel_x = %d (x=%g), pixel_y = %d (y=%g)\n",
	    pixel_x, x, pixel_y, y);
#endif

    if (plot->err || plot->spec->markers == NULL) {
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

    data_x = plot->spec->data->val;
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
	    return TRUE;
	}
    }

 repeat:

    /* find the best-matching data point */
    for (t=0; t<plot->spec->nobs; t++) {
	if (na(data_x[t]) || na(data_y[t])) {
	    continue;
	}
	xdiff = fabs(data_x[t] - x);
	ydiff = fabs(data_y[t] - y);
	dist = sqrt(xdiff * xdiff + ydiff * ydiff);
#if GPDEBUG > 4
	fprintf(stderr, " obs %d: x=%g, y=%g, dist=%g\n",
		t, data_x[t], data_y[t], dist);
#endif
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
    if (t >= 0) {
	fprintf(stderr, " best_match=%d, with data_x[%d]=%g, data_y[%d]=%g\n",
		t, t, data_x[t], t, data_y[t]);
    } else {
	fprintf(stderr, " no 'best match' for %g, %g\n", x, y);
    }
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

static inline int use_integer_format (int xint, double xdata)
{
    return (xint && fabs(xdata) < 1.0e7);
}

static void heatmap_show_z (png_plot *plot, double x, double y,
			    gchar *label)
{
    gretl_matrix *z = plot->spec->auxdata;
    double zij;
    int i = nearbyint(y);
    int j = nearbyint(x);
    int n = z->rows;

    if (i >= 0 && i < n && j >= 0 && j < n) {
	zij = gretl_matrix_get(z, i, j);
	if (isnan(zij)) {
	    gtk_label_set_text(GTK_LABEL(plot->cursor_label), "");
	} else {
	    sprintf(label, "r = %.3f", zij);
	    gtk_label_set_text(GTK_LABEL(plot->cursor_label), label);
	}
    }
}

static gint
plot_motion_callback (GtkWidget *widget, GdkEventMotion *event,
		      png_plot *plot)
{
    GdkModifierType state;
    gchar label[48], label_y[24];
    const char *xfmt = NULL;
    const char *yfmt = NULL;
    int do_label;
    int x, y;

    plot = plot_get_current(plot);
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

	if (plot->spec->code == PLOT_HEATMAP) {
	    if (plot->spec->auxdata != NULL) {
		heatmap_show_z(plot, data_x, data_y, label);
	    }
	    return TRUE;
	}

	if (!cant_do_labels(plot) && !labels_frozen(plot) &&
	    !plot_is_zooming(plot) && !na(data_y)) {
	    identify_point(plot, x, y, data_x, data_y);
	}

	if (do_label) {
	    if (plot->pd == 4 || plot->pd == 12) {
		x_to_date(data_x, plot->pd, label);
	    } else if (plot->spec->flags & GPT_TIMEFMT) {
                gretl_strftime(label, sizeof label, "%Y-%m-%d",
                               (gint64) data_x, NADBL);
	    } else if (xfmt != NULL) {
		sprintf(label, xfmt, data_x);
	    } else if (use_integer_format(plot->xint, data_x)) {
		sprintf(label, "%7d", (int) data_x);
	    } else {
		sprintf(label, "%#7.4g", data_x);
	    }
	}

	if (do_label && !na(data_y)) {
	    if (plot_has_png_coords(plot)) {
		if (yfmt != NULL) {
		    sprintf(label_y, yfmt, data_y);
		} else if (use_integer_format(plot->yint, data_y)) {
		    sprintf(label_y, " %-7d", (int) data_y);
		} else {
		    sprintf(label_y, " %#-7.4g", data_y);
		}
	    } else if (use_integer_format(plot->yint, data_y)) {
		sprintf(label_y, " %-7d", (int) data_y);
	    } else {
		sprintf(label_y, " %#-6.3g", data_y);
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
    gchar *tmpname = NULL;
    char line[256], fontspec[128];
    int gotterm = 0;
    int err = 0;

    adjust_fontspec_string(fontspec, grfont, ADD_COMMA);

#if GPDEBUG
    fprintf(stderr, "font choice: grfont='%s', fontspec='%s'\n",
	    grfont, fontspec);
#endif

    fp = gretl_fopen(plot->spec->fname, "r");
    if (fp == NULL) {
	file_read_errbox(plot->spec->fname);
	return;
    }

    tmpname = gretl_make_dotpath("gpttmp.XXXXXX");
    ftmp = gretl_tempfile_open(tmpname);
    if (ftmp == NULL) {
	g_free(tmpname);
	fclose(fp);
	gui_errmsg(E_FOPEN);
	return;
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (!gotterm && strncmp(line, "set term", 8) == 0) {
	    err = substitute_graph_font(line, fontspec);
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
	err = gretl_rename(tmpname, plot->spec->fname);
	if (!err) {
	    err = gnuplot_make_image(plot->spec->fname);
	}
    }

    g_free(tmpname);

    if (err) {
	gui_errmsg(err);
    } else {
	free(plot->spec->fontstr);
	plot->spec->fontstr = gretl_strdup(fontspec);
	render_png(plot, PNG_REDISPLAY);
    }
}

static const gchar *menu_item_get_text (GtkMenuItem *item)
{
    GtkWidget *label = gtk_bin_get_child(GTK_BIN(item));

    return gtk_label_get_text(GTK_LABEL(label));
}

static gint color_popup_activated (GtkMenuItem *item, gpointer data)
{
    png_plot *plot = widget_get_plot(item);
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
				  plot, plot->shell);
    } else if (!strcmp(menu_string, _("Copy to clipboard"))) {
#ifdef G_OS_WIN32
	win32_process_graph(plot, WIN32_TO_CLIPBOARD);
#else
	set_plot_for_copy(plot);
#endif
    }

#ifdef G_OS_WIN32
    else if (!strcmp(menu_string, _("Print"))) {
	win32_process_graph(plot, WIN32_TO_PRINTER);
    }
#endif

    plot->spec->flags &= ~GPT_MONO;

    return TRUE;
}

static void show_numbers_from_markers (GPT_SPEC *spec)
{
    PRN *prn;
    double x, y, pi2;
    double mod, freq;
    int dcomma = 0;
    int i, err = 0;

    if (bufopen(&prn)) {
	return;
    }

    pi2 = 2.0 * M_PI;
    dcomma = get_local_decpoint() != '.';

    pputs(prn, _("roots (real, imaginary, modulus, frequency)"));
    pputs(prn, "\n\n");

    if (dcomma) {
	gretl_push_c_numeric_locale();
    }

    for (i=0; i<spec->n_markers; i++) {
	if (sscanf(spec->markers[i], "%lf,%lf", &x, &y) == 2) {
	    freq = gretl_matrix_get(spec->data, i, 0) / pi2;
	    mod = gretl_matrix_get(spec->data, i, 1);
	    if (dcomma) {
		gretl_pop_c_numeric_locale();
		pprintf(prn, "%2d: (%7.4f  %7.4f  %7.4f  %7.4f)\n", i+1,
			x, y, mod, freq);
		gretl_push_c_numeric_locale();
	    } else {
		pprintf(prn, "%2d: (%7.4f, %7.4f, %7.4f, %7.4f)\n", i+1,
			x, y, mod, freq);
	    }
	} else {
	    err = E_DATA;
	    break;
	}
    }

    if (dcomma) {
	gretl_pop_c_numeric_locale();
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

static void add_to_session_callback (png_plot *plot)
{
    char fullname[MAXLEN] = {0};
    int err, type;

    type = (plot->spec->code == PLOT_BOXPLOTS)? GRETL_OBJ_PLOT :
	GRETL_OBJ_GRAPH;

    err = gui_add_graph_to_session(plot->spec->fname, fullname, type);

    if (!err) {
	remove_png_term_from_plot(fullname, plot->spec);
	plot->status |= PLOT_SAVED;
	plot->status |= PLOT_TERM_HIDDEN;
    }
}

static GList *plot_get_siblings (png_plot *plot)
{
    if (in_collection(plot)) {
	return plot->mp->list;
    } else {
	return NULL;
    }
}

static int real_plot_rescale (png_plot *plot, double scale)
{
    FILE *fp = NULL;
    int err = 0;

    plot->spec->scale = scale;
    gnuplot_png_init(plot, &fp);

    if (fp == NULL) {
	err = E_FOPEN;
	gui_errmsg(err);
    } else {
	set_png_output(plot->spec);
	plotspec_print(plot->spec, fp);
	fclose(fp);
	unset_png_output(plot->spec);
    }

    return err;
}

static int regenerate_pixbuf (png_plot *plot)
{
    int err = gnuplot_make_image(plot->spec->fname);

    if (!err) {
	plot_invalidate_pixbuf(plot);
	plot->pbuf = pixbuf_from_file(plot);
	if (plot->pbuf == NULL) {
	    err = 1;
	}
    }

    return err;
}

static int rescale_siblings (png_plot *plot, double scale)
{
    GList *L = plot_get_siblings(plot);
    png_plot *sib;
    int w, h, err = 0;

    while (L != NULL && !err) {
	sib = L->data;
	if (sib != plot) {
	    err = real_plot_rescale(sib, scale);
	    if (!err) {
		err = regenerate_pixbuf(sib);
	    }
	    if (!err) {
		w = gdk_pixbuf_get_width(sib->pbuf);
		h = gdk_pixbuf_get_height(sib->pbuf);
		if (w != sib->pixel_width || h != sib->pixel_height) {
		    resize_png_plot(sib, w, h, 1);
		}
	    }
	}
	L = L->next;
    }

    return err;
}

static void sensitize_scale_buttons (png_plot *plot)
{
    if (plot_is_zoomed(plot)) {
	gtk_widget_set_sensitive(plot->up_icon, FALSE);
	gtk_widget_set_sensitive(plot->down_icon, FALSE);
    } else {
	double scale = plot->spec->scale;

	gtk_widget_set_sensitive(plot->up_icon, scale < max_graph_scale());
	gtk_widget_set_sensitive(plot->down_icon, scale > min_graph_scale());
    }
}

static void plot_do_rescale (png_plot *plot, int mod)
{
    double scale = 1.0;
    int err = 0;

    plot = plot_get_current(plot);

    if (mod == 0) {
	/* reset to default */
	if (plot->spec->scale == 1.0) {
	    /* no-op */
	    return;
	}
    } else {
	/* enlarge (mod = 1) or shrink (mod = -1) */
	scale = next_graph_scale(plot->spec->scale, mod);
	if (na(scale)) {
	    gdk_window_beep(plot->window);
	    return;
	}
    }

    err = real_plot_rescale(plot, scale);

    if (!err) {
	repaint_png(plot, PNG_REDISPLAY);
	if (in_collection(plot)) {
	    rescale_siblings(plot, scale);
	}
    }

    sensitize_scale_buttons(plot);
}

static void show_all_labels (png_plot *plot)
{
    FILE *fp = NULL;

    if (plot->spec->labeled != NULL) {
	/* Destroy record of labeling of specific points via
	   "brushing", if present.
	*/
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
    FILE *fp = NULL;

    if (plot->spec->flags & GPT_PRINT_MARKERS) {
	plot->spec->flags &= ~GPT_PRINT_MARKERS;
    }

    gnuplot_png_init(plot, &fp);
    if (fp == NULL) {
	gui_errmsg(E_FOPEN);
	return;
    }
    plotspec_print(plot->spec, fp);
    fclose(fp);

    repaint_png(plot, PNG_REDISPLAY);
    plot->format &= ~PLOT_MARKERS_UP;
}

static void prepare_for_zoom (png_plot *plot)
{
    GdkCursor* cursor = gdk_cursor_new(GDK_CROSSHAIR);

    if (cursor != NULL) {
	gdk_window_set_cursor(plot->window, cursor);
	gdk_cursor_unref(cursor);
    }
    plot->status |= PLOT_ZOOMING;
    gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid,
		       _(" Drag to define zoom rectangle"));
}

static gint plot_popup_activated (GtkMenuItem *item, gpointer data)
{
    png_plot *plot = (png_plot *) data;
    GtkWidget *shell = plot->shell;
    const gchar *item_string;
    int killplot = 0;

    plot = plot_get_current(plot);
    item_string = menu_item_get_text(item);

    if (!strcmp(item_string, _("Add another curve..."))) {
	dist_graph_add(plot);
    } else if (!strcmp(item_string, _("Save as PNG..."))) {
	plot->spec->termtype = GP_TERM_PNG;
        file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC,
				  plot, shell);
    } else if (!strcmp(item_string, _("Save as PDF..."))) {
	plot->spec->termtype = GP_TERM_PDF;
	pdf_ps_dialog(plot->spec, shell);
    } else if (!strcmp(item_string, _("Save as postscript (EPS)..."))) {
	plot->spec->termtype = GP_TERM_EPS;
	pdf_ps_dialog(plot->spec, shell);
    } else if (!strcmp(item_string, _("Save to session as icon"))) {
	add_to_session_callback(plot);
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
    } else if (!strcmp(item_string, _("Display PDF"))) {
	graph_display_pdf(plot);
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
	plot_show_font_selector(plot, plot->spec->fontstr);
    } else if (!strcmp(item_string, _("Close"))) {
	killplot = 1;
    } else if (!strcmp(item_string, _("Delete plot"))) {
        killplot = 2;
    } else if (!strcmp(item_string, _("Close collection"))) {
	killplot = 1;
    } else if (!strcmp(item_string, _("Extract plot"))) {
	killplot = 3;
    }

    if (killplot == 1) {
	/* trash a singleton plot or entire collection */
	gtk_widget_destroy(shell);
    } else if (killplot == 2) {
	/* trash a specific plot in a collection */
	plot_collection_remove_plot(plot, 1);
    } else if (killplot == 3) {
	/* pull a plot out of a collection */
	plot_collection_remove_plot(plot, 0);
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
	N_("Copy to clipboard"),
	N_("Save to session as icon"),
	N_("Freeze data labels"),
	N_("All data labels"),
	N_("Clear data labels"),
	N_("Zoom..."),
#ifdef G_OS_WIN32
	N_("Print"),
#endif
	N_("Display PDF"),
	N_("OLS estimates"),
	N_("Numerical values"),
	N_("Numerical summary"),
	N_("Edit"),
	N_("Font"),
	N_("Help"),
        N_("Close"),
	N_("Delete plot"),
	N_("Close collection"),
	N_("Extract plot"),
        NULL
    };
    const char *zoomed_items[] = {
	N_("Restore full view"),
	N_("Replace full view"),
	N_("Close"),
	NULL
    };
    const char **plot_items;
    int i, in_coll;

    plot->popup = gtk_menu_new();

    if (plot_is_zoomed(plot)) {
	plot_items = zoomed_items;
    } else {
	plot_items = regular_items;
    }

    /* geoplot FIXME: for several menu items below: either
       support them for geoplot or don't show them
    */

    in_coll = in_collection(plot);

    i = 0;
    while (plot_items[i]) {
	int colorpop = 0;

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
	if ((plot_has_controller(plot) || plot_is_editable(plot)) &&
	    !strcmp(plot_items[i], "Font")) {
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
	if (in_coll) {
	    if (!strcmp(plot_items[i], "Close") ||
		!strncmp(plot_items[i], "Save to", 7)) {
		i++;
		continue;
	    }
	} else if (!strncmp(plot_items[i], "Close ", 6) ||
		   !strncmp(plot_items[i], "Delete ", 7) ||
		   !strncmp(plot_items[i], "Extract", 7)) {
	    i++;
	    continue;
	}

	item = gtk_menu_item_new_with_label(_(plot_items[i]));

#ifdef G_OS_WIN32
	if (!strcmp(plot_items[i], "Copy to clipboard") ||
	    !strcmp(plot_items[i], "Save as Windows metafile (EMF)...") ||
	    !strcmp(plot_items[i], "Print")) {
	    colorpop = 1;
	}
#else
	if (!strcmp(plot_items[i], "Copy to clipboard") ||
	    !strcmp(plot_items[i], "Save as Windows metafile (EMF)...")) {
	    colorpop = 1;
	}
#endif
	if (colorpop) {
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
    FILE *fp;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "redisplay_edited_plot: plot = %p\n", (void *) plot);
#endif

    /* get the actual target of the edit */
    plot = plot_get_current(plot);

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
    err = gnuplot_make_image(plot->spec->fname);
    if (err) {
	gui_errmsg(err);
	return err;
    }

    /* reset format flags */
    set_plot_format_flags(plot);

    /* grab (possibly modified) data ranges */
    get_plot_ranges(plot, plot->spec->code);

    /* put the newly created PNG onto the plot canvas */
    return render_png(plot, PNG_REDISPLAY);
}

static gchar *recover_plot_header (char *line, FILE *fp)
{
    gchar *setterm = NULL;

    while (fgets(line, MAXLEN-1, fp)) {
	if (setterm == NULL && commented_term_line(line)) {
	    gchar *dpath = gretl_make_dotpath("gretltmp.png");
	    GString *gs = g_string_new(line + 2);

	    g_string_append(gs, "set encoding utf8\n");
	    g_string_append_printf(gs, "set output \"%s\"\n", dpath);
	    g_free(dpath);
	    setterm = g_string_free(gs, FALSE);
	} else if (!strncmp(line, "plot", 4)) {
	    break;
	}
    }

    return setterm;
}

/* preparation for redisplaying graph: here we handle the case where
   we're switching to a zoomed view (by use of a temporary gnuplot
   source file); then we get gnuplot to create a new PNG.
 */

static int repaint_png (png_plot *plot, int view)
{
    gchar *altname = NULL;
    int do_zoom = view == PNG_ZOOM;
    int err = 0;

    if (do_zoom || (plot->status & PLOT_TERM_HIDDEN)) {
	/* we'll have to rewrite the plot file */
	gchar *setterm = NULL;
	char line[MAXLEN];
	FILE *fpin, *fpout;
	int gotterm = 0;

	fpin = gretl_fopen(plot->spec->fname, "r");
	if (fpin == NULL) {
	    return E_FOPEN;
	}
	if (plot->status & PLOT_TERM_HIDDEN) {
	    setterm = recover_plot_header(line, fpin);
	    rewind(fpin);
	}

	altname = gretl_make_dotpath("altplot.gp");
	fpout = gretl_fopen(altname, "w");
	if (fpout == NULL) {
	    fclose(fpin);
	    g_free(altname);
	    return E_FOPEN;
	}

	if (do_zoom) {
	    /* write revised range into auxiliary gnuplot source file */
	    gretl_push_c_numeric_locale();
	    fprintf(fpout, "set xrange [%g:%g]\n", plot->zoom_xmin,
		    plot->zoom_xmax);
	    fprintf(fpout, "set yrange [%g:%g]\n", plot->zoom_ymin,
		    plot->zoom_ymax);
	    gretl_pop_c_numeric_locale();
	}

	while (fgets(line, MAXLEN-1, fpin)) {
	    if (setterm != NULL && !gotterm && commented_term_line(line)) {
		fputs(setterm, fpout);
		gotterm = 1;
	    } else if (!do_zoom) {
		fputs(line, fpout);
	    } else if (strncmp(line, "set xrange", 10) &&
		       strncmp(line, "set yrange", 10)) {
		fputs(line, fpout);
	    }
	}

	fclose(fpin);
	fclose(fpout);

	err = gnuplot_make_image(altname);
    } else {
	err = gnuplot_make_image(plot->spec->fname);
    }

    if (altname != NULL) {
	gretl_remove(altname);
	g_free(altname);
    }

    if (err) {
	gui_errmsg(err);
	return err;
    }

    return render_png(plot, view);
}

/* On replacing "full" view of a plot with zoomed view:
   sync the zoomed x and y ranges to both the plot and
   its attached GPT_SPEC pointer.
*/

static void plot_sync_xy_ranges (png_plot *plot)
{
    double *rx = plot->spec->range[GP_X_RANGE];
    double *ry = plot->spec->range[GP_Y_RANGE];

    rx[0] = plot->xmin = plot->zoom_xmin;
    rx[1] = plot->xmax = plot->zoom_xmax;
    ry[0] = plot->ymin = plot->zoom_ymin;
    ry[1] = plot->ymax = plot->zoom_ymax;

    plot->zoom_xmin = plot->zoom_xmax = 0.0;
    plot->zoom_ymin = plot->zoom_ymax = 0.0;
}

/* with a zoomed version of the current plot in place,
   replace the full version with the zoom
*/

static int zoom_replaces_plot (png_plot *plot)
{
    FILE *fpin, *fpout;
    gchar *temp = NULL;
    char line[MAXLEN];
    int err = 0;

    fpin = gretl_fopen(plot->spec->fname, "r");
    if (fpin == NULL) {
	return 1;
    }

    temp = gretl_make_dotpath("zoomtmp.XXXXXX");
    fpout = gretl_tempfile_open(temp);
    if (fpout == NULL) {
	fclose(fpin);
	g_free(temp);
	return 1;
    }

    /* first copy into temp file, replacing data ranges */

    gretl_push_c_numeric_locale();
    while (fgets(line, MAXLEN-1, fpin)) {
	if (!strncmp(line, "set xrange", 10)) {
	    fprintf(fpout, "set xrange [%g:%g]\n", plot->zoom_xmin,
		    plot->zoom_xmax);
	} else if (!strncmp(line, "set yrange", 10)) {
	    fprintf(fpout, "set yrange [%g:%g]\n", plot->zoom_ymin,
		    plot->zoom_ymax);
	} else {
	    fputs(line, fpout);
	}
    }
    gretl_pop_c_numeric_locale();

    fclose(fpout);
    fclose(fpin);

    /* then replace the original graph source file */

    err = gretl_copy_file(temp, plot->spec->fname);
    if (err) {
	gui_errmsg(err);
    } else {
	plot_sync_xy_ranges(plot);
	plot->status ^= PLOT_ZOOMED;
    }

    gretl_remove(temp);
    g_free(temp);

    sensitize_scale_buttons(plot);

    return err;
}

static gint plot_button_release (GtkWidget *widget, GdkEventButton *event,
				 png_plot *plot)
{
    plot = plot_get_current(plot);

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
    plot = plot_get_current(plot);

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

    if (right_click(event)) {
	build_plot_menu(plot);
	gtk_menu_popup(GTK_MENU(plot->popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
    }

    return TRUE;
}

static gboolean
plot_key_handler (GtkWidget *w, GdkEventKey *event, png_plot *plot)
{
    int Ctrl = (event->state & GDK_CONTROL_MASK);
    guint k = event->keyval;

#ifdef OS_OSX
    if (!Ctrl && cmd_key(event)) {
	/* treat Command as Ctrl */
	Ctrl = 1;
    }
#endif

    if (plot_is_editable(plot) &&
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

    if (Ctrl && k == GDK_c) {
#ifdef G_OS_WIN32
	win32_process_graph(plot, WIN32_TO_CLIPBOARD);
#else
	set_plot_for_copy(plot);
#endif
	return TRUE;
    }

    if (in_collection(plot)) {
	return TRUE;
    }

    switch (k) {
    case GDK_q:
    case GDK_Q:
#ifdef OS_OSX
    case GDK_w:
    case GDK_W:
#endif
	gtk_widget_destroy(plot->shell);
	break;
    case GDK_s:
    case GDK_S:
	add_to_session_callback(plot);
	break;
    case GDK_z:
    case GDK_Z:
	if (!plot_not_zoomable(plot)) {
	    prepare_for_zoom(plot);
	}
	break;
    default:
	break;
    }

    return TRUE;
}

#if GTK_MAJOR_VERSION >= 3

static
void plot_draw (GtkWidget *canvas, cairo_t *cr, gpointer data)
{
    png_plot *plot = plot_get_current(data);

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

static void record_coordinate_info (png_plot *plot, png_bounds *b)
{
    if (b->height > 0) {
	plot->pixel_height = b->height;
    }
    if (b->width > 0) {
	plot->pixel_width = b->width;
    }
    plot->status |= PLOT_PNG_COORDS;
    plot->pixel_xmin = b->xleft;
    plot->pixel_xmax = b->xright;
    plot->pixel_ymin = plot->pixel_height - b->ytop;
    plot->pixel_ymax = plot->pixel_height - b->ybot;
    plot->xmin = b->xmin;
    plot->xmax = b->xmax;
    plot->ymin = b->ymin;
    plot->ymax = b->ymax;

#if POINTS_DEBUG
    fprintf(stderr, "get_png_bounds_info():\n"
	    " xmin=%d xmax=%d ymin=%d ymax=%d\n",
	    plot->pixel_xmin, plot->pixel_xmax,
	    plot->pixel_ymin, plot->pixel_ymax);
    fprintf(stderr, "using px_height %d, px_width %d\n",
	    plot->pixel_height, plot->pixel_width);
#endif
}

#if GTK_MAJOR_VERSION == 2

static void pixmap_sync (png_plot *plot)
{
    GList *L = plot->mp->list;
    png_plot *sibling;

    while (L != NULL) {
	sibling = L->data;
	if (sibling != plot) {
	    sibling->pixmap = plot->pixmap;
	}
	L = L->next;
    }
}

#endif

static int resize_png_plot (png_plot *plot, int width, int height,
			    int follower)
{
    png_bounds b;

    plot->pixel_width = width;
    plot->pixel_height = height;

    if (!follower) {
	gtk_widget_set_size_request(GTK_WIDGET(plot->canvas),
				    plot->pixel_width, plot->pixel_height);
#if GTK_MAJOR_VERSION == 2
	g_object_unref(plot->pixmap);
	plot->pixmap = gdk_pixmap_new(plot->window,
				      plot->pixel_width,
				      plot->pixel_height,
				      -1);
	if (in_collection(plot)) {
	    pixmap_sync(plot);
	}
#endif
    }

    if (plot->status & (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE)) {
	return 0;
    }

    b.width = width;
    b.height = height;

    /* try revising the gnuplot bounds info? */
    if (plot_has_png_coords(plot) &&
	get_png_bounds_info(&b) == GRETL_PNG_OK) {
	record_coordinate_info(plot, &b);
    } else {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
    }

    return 0;
}

static GdkPixbuf *pixbuf_from_file (png_plot *plot)
{
    char fname[MAXLEN];
    GError *gerr = NULL;
    GdkPixbuf *pbuf;

    gretl_build_path(fname, gretl_dotdir(), "gretltmp.png", NULL);
    pbuf = gdk_pixbuf_new_from_file(fname, &gerr);

    if (gerr != NULL) {
	errbox(gerr->message);
	g_error_free(gerr);
    } else if (pbuf == NULL) {
	file_read_errbox(fname);
	gretl_remove(fname);
    } else if (gdk_pixbuf_get_width(pbuf) == 0 ||
	       gdk_pixbuf_get_height(pbuf) == 0) {
	errbox(_("Malformed PNG file for graph"));
	g_object_unref(pbuf);
	pbuf = NULL;
    }

    return pbuf;
}

/* The last step in displaying a graph (or redisplaying after some
   change has been made): grab the gnuplot-generated PNG file, make a
   pixbuf out of it, and draw the pixbuf onto the canvas of the plot
   window.
*/

static int render_png (png_plot *plot, viewcode view)
{
    char pngname[MAXLEN];
    gint width, height;
    GdkPixbuf *pbuf;

    if (view == PNG_REDISPLAY || view == PNG_ZOOM || view == PNG_UNZOOM) {
	/* we need to read a revised PNG file */
	plot_invalidate_pixbuf(plot);
    }

#if COLLDEBUG
    fprintf(stderr, "\nrender_png: plot %p, pixbuf %p, view = %s\n",
	    (void *) plot, (void *) plot->pbuf, viewstr(view));
#endif

    if (plot->pbuf != NULL) {
	pbuf = plot->pbuf;
    } else {
	pbuf = pixbuf_from_file(plot);
	if (pbuf == NULL) {
	    return 1;
	}
    }

    width = gdk_pixbuf_get_width(pbuf);
    height = gdk_pixbuf_get_height(pbuf);

    if (width != plot->pixel_width || height != plot->pixel_height) {
	resize_png_plot(plot, width, height, 0);
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

    /* FIXME case of singleton plot? */
    plot->pbuf = pbuf;

    if (view != PNG_REPLACE) {
	gretl_build_path(pngname, gretl_dotdir(), "gretltmp.png", NULL);
	gretl_remove(pngname);
    }

    if (view != PNG_START) {
	/* we're changing the view, so refresh the whole canvas */
	redraw_plot_rectangle(plot, NULL);
	if (view == PNG_ZOOM) {
	    plot->status |= PLOT_ZOOMED;
	    sensitize_scale_buttons(plot);
	} else if (view == PNG_UNZOOM) {
	    plot->status ^= PLOT_ZOOMED;
	    sensitize_scale_buttons(plot);
	}
    }

#ifdef G_OS_WIN32
    /* somehow the plot can end up underneath */
    gtk_window_present(GTK_WINDOW(plot->shell));
#endif

    return 0;
}

#if GTK_MAJOR_VERSION == 2

static void plot_nullify_surface (png_plot *plot)
{
    plot->pixmap = NULL;
    plot->savebuf = NULL;
}

#endif

static void plot_destroy_surface (png_plot *plot)
{
#if GTK_MAJOR_VERSION >= 3
    if (plot->cs != NULL) {
	cairo_surface_destroy(plot->cs);
    }
#else
    if (plot->pixmap != NULL) {
	g_object_unref(plot->pixmap);
    }
    if (plot->savebuf != NULL) {
	g_object_unref(plot->savebuf);
    }
#endif
}

static void destroy_png_plot (GtkWidget *w, png_plot *plot)
{
    if (in_collection(plot)) {
	GList *L = g_list_first(plot->mp->list);
	png_plot *sibling;

	while (L != NULL) {
	    sibling = L->data;
	    if (sibling != plot) {
#if GTK_MAJOR_VERSION == 2
		plot_nullify_surface(sibling);
#endif
		sibling->shell = NULL;
		sibling->mp = NULL;
		destroy_png_plot(NULL, sibling);
	    }
	    L = L->next;
	}
    }

    if (!plot_is_saved(plot) && plot->spec != NULL) {
	/* delete temporary plot source file? */
	if (strstr(plot->spec->fname, gretl_dotdir())) {
	    gretl_remove(plot->spec->fname);
	}
    }

#if GPDEBUG
    fprintf(stderr, "destroy_png_plot: plot = %p, spec = %p\n",
	    (void *) plot, (void *) plot->spec);
#endif

#ifndef G_OS_WIN32
    if (copyplot == plot) {
	copyplot = NULL;
    }
#endif

    if (plot->spec != NULL) {
	plotspec_destroy(plot->spec);
    }

    plot_destroy_surface(plot);

    if (plot->pbuf != NULL) {
	g_object_unref(plot->pbuf);
    }

    if (plot->shell != NULL) {
	g_object_unref(plot->shell);
    }

    if (plot == plot_collection) {
	plot_collection = NULL;
    }

    if (plot->mp != NULL) {
	g_list_free(plot->mp->list);
	free(plot->mp);
    }

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
    int annual = 0;
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
	} else if (!strncmp(line, "# timeseries 1", 13)) {
	    annual = 1;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    /* now try getting accurate coordinate info from
       auxiliary file (or maybe PNG file)
    */
    if (get_png_bounds_info(&b) == GRETL_PNG_OK) {
	record_coordinate_info(plot, &b);
	got_x = got_y = 1;
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
	if (annual) {
	    /* years on x-axis: show as integer */
	    plot->xint = 1;
	} else if ((plot->xmax - plot->xmin) /
	    (plot->pixel_xmax - plot->pixel_xmin) >= 1.0) {
	    /* show x-axis variable as integer */
	    plot->xint = 1;
	}
	if ((plot->ymax - plot->ymin) /
	    (plot->pixel_ymax - plot->pixel_ymin) >= 1.0) {
	    /* show y-axis variable as integer */
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
    plot->mp = NULL;
    plot->canvas = NULL;
    plot->popup = NULL;
    plot->statusbar = NULL;
    plot->cursor_label = NULL;
    plot->pbuf = NULL;
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
    plot->toolbar = NULL;
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

static int plot_add_shell (png_plot *plot, const char *name)
{
    GtkWidget *vbox;
    GtkWidget *canvas_hbox;
    GtkWidget *status_hbox;
    GtkWidget *status_area;
    gchar *title;

    plot->shell = gretl_gtk_window();
    g_object_ref(plot->shell);
    g_object_set_data(G_OBJECT(plot->shell), "plot", plot);
#if 0 /* do we want this? */
    gtk_window_set_position(GTK_WINDOW(plot->shell), GTK_WIN_POS_MOUSE);
#endif

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
    status_area = gtk_event_box_new();
    gtk_box_pack_start(GTK_BOX(vbox), status_area, FALSE, FALSE, 0);
    status_hbox = gtk_hbox_new(FALSE, 2);
    gtk_container_add(GTK_CONTAINER(status_area), status_hbox);

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

    /* more on canvas and window */
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
    g_signal_connect(G_OBJECT(plot->canvas), "key-press-event",
		     G_CALLBACK(plot_key_handler), plot);
    gtk_widget_show_all(plot->shell);

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

    return 0;
}

static void plot_handle_specials (png_plot *plot)
{
    if (plot->spec->code == PLOT_PERIODOGRAM) {
	/* the x2 axis gets broken, and also the x axis if it's
	   in degrees or radians, on zooming
	*/
	plot->status |= PLOT_DONT_ZOOM;
    } else if (plot->spec->code == PLOT_ROOTS ||
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
    } else if (plot->spec->flags & GPT_XW) {
	plot->pixel_width = GP_XW_WIDTH;
	plot->pixel_height = GP_HEIGHT;
    }

    if (plot->spec->scale != 1.0) {
	plot_get_scaled_dimensions(&plot->pixel_width,
				   &plot->pixel_height,
				   plot->spec->scale);
    }
}

/* note: @fname is the name of the file containing the
   plot commands.
*/

static int gnuplot_show_png (const char *fname,
			     const char *name,
			     void *session_ptr,
			     png_plot *coll)
{
    png_plot *plot;
    int err = 0;

    if (*fname == '\0') {
	return 0;
    }

    gretl_error_clear();

#if GPDEBUG || WINDEBUG
    fprintf(stderr, "gnuplot_show_png:\n fname='%s', saved=%d\n",
	    fname, (session_ptr != NULL));
#endif

    plot = png_plot_new();
    if (plot == NULL) {
	return E_ALLOC;
    }

    plot->spec = plotspec_new();
    if (plot->spec == NULL) {
	free(plot);
	return E_ALLOC;
    }

    strcpy(plot->spec->fname, fname);

    if (session_ptr != NULL) {
	plot->status |= PLOT_SAVED;
    }

    /* make png plot struct accessible via spec */
    plot->spec->ptr = plot;

    /* Parse the gnuplot source file.  If we hit errors here,
       flag this, but it's not necessarily a show-stopper in
       terms of simply displaying the graph -- unless we get
       E_FOPEN.
    */
    plot->err = read_plotspec_from_file(plot);

    if (plot->err == E_FOPEN) {
	plotspec_destroy(plot->spec);
	free(plot);
	return E_FOPEN;
    }

#if GPDEBUG || WINDEBUG
    fprintf(stderr, "gnuplot_show_png: read_plotspec_from_file returned %d\n",
	    plot->err);
#endif

    if (plot->err || cant_edit(plot->spec->code)) {
	plot->status |= (PLOT_DONT_EDIT | PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
    } else {
	set_plot_format_flags(plot);
    }

    plot_handle_specials(plot);

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

    if (coll != NULL) {
	/* add this plot to current collection? */
	int attached = 0;

	if (!err && plot_can_be_collected(coll, plot)) {
	    err = plot_collection_add_plot(coll, plot);
	    attached = (err == 0);
	}
	if (err || attached) {
	    return err;
	}
    }

    plot_add_shell(plot, name);
    err = render_png(plot, PNG_START);
    if (err) {
	gtk_widget_destroy(plot->shell);
    } else {
	g_object_set_data(G_OBJECT(plot->shell), "object", plot);
	if (session_ptr != NULL) {
	    g_object_set_data(G_OBJECT(plot->shell),
			      "session-ptr", session_ptr);
	} else if (coll == NULL && do_collect_plots()) {
	    set_plot_collection(plot);
	}
    }

    return err;
}

int gnuplot_show_map (gretl_bundle *mb)
{
    const char *fname = NULL;
    const char *mapname = NULL;
    const gretl_matrix *dims = NULL;
    png_plot *plot;
    png_bounds b;
    int err = 0;

    gretl_error_clear();

    fname = gretl_bundle_get_string(mb, "plotfile", &err);
    if (!err) {
	dims = gretl_bundle_get_matrix(mb, "dims", &err);
    }
    if (err) {
	gui_errmsg(err);
	return err;
    }

    gretl_error_clear();

    plot = png_plot_new();
    if (plot == NULL) {
	return E_ALLOC;
    }

    plot->spec = plotspec_new();
    if (plot->spec == NULL) {
	free(plot);
	return E_ALLOC;
    }

    strcpy(plot->spec->fname, fname);
    plot->spec->ptr = plot;

    err = gretl_test_fopen(plot->spec->fname, "rb");
    if (err) {
	gretl_errmsg_sprintf(_("Couldn't read '%s'"), plot->spec->fname);
	err = E_FOPEN;
    } else {
	err = gnuplot_make_image(plot->spec->fname);
    }

    if (err) {
	plotspec_destroy(plot->spec);
	free(plot);
	return err;
    }

    plot->spec->code = PLOT_GEOMAP;
    plot->status |= (PLOT_DONT_EDIT | PLOT_DONT_ZOOM);
    plot->pixel_width = dims->val[0];
    plot->pixel_height = dims->val[1];
#if GPDEBUG
    fprintf(stderr, "map, via bundle: pixwidth %d, pixheight %d\n",
	    plot->pixel_width, plot->pixel_height);
#endif

    if (get_png_bounds_info(&b) == GRETL_PNG_OK) {
	record_coordinate_info(plot, &b);
	plot->status |= PLOT_CURSOR_LABEL;
	plot->status |= PLOT_PNG_COORDS;
	plot->status |= PLOT_HAS_XRANGE;
	plot->status |= PLOT_HAS_YRANGE;
#if GPDEBUG
	fprintf(stderr, "map, via png bounds: %d, %d\n",
		plot->pixel_width, plot->pixel_height);
#endif
    } else {
	plot->status |= PLOT_DONT_MOUSE;
    }

    mapname = gretl_bundle_get_string(mb, "mapname", NULL);
    plot_add_shell(plot, mapname != NULL ? mapname : "map");
    err = render_png(plot, PNG_START);

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
    int dcp = do_collect_plots();
    png_plot *pp = NULL;

    if (dcp && plot_collection != NULL) {
	if (dcp > 1) {
	    /* explicitly "on" */
	    pp = plot_collection;
	} else {
	    /* "auto": time-limited */
	    gint64 now = gretl_monotonic_time();
	    gint64 ptm = plot_collection_get_mtime();

	    if (now - ptm < 1.25e6) {
		/* time limit 1.25 seconds */
		pp = plot_collection;
	    }
	}
    }

    gnuplot_show_png(gretl_plotfile(), NULL, NULL, pp);
}

/* @fname is the name of a plot command file from the
   current session and @title is its display name;
   @session_ptr is a pointer to the session object.
*/

void display_session_graph (const char *fname,
			    const char *title,
			    void *session_ptr)
{
    gchar *fullname = NULL;
    int err = 0;

    if (g_path_is_absolute(fname)) {
	fullname = g_strdup(fname);
    } else {
	fullname = gretl_make_dotpath(fname);
    }

#if WINDEBUG
    fprintf(stderr, "display_session_graph:\n fullname = '%s'\n",
	    fullname);
    fprintf(stderr, " gnuplot path = '%s'\n", gretl_gnuplot_path());
    fprintf(stderr, " call add_png_term_to_plot\n");
#endif

    err = add_png_term_to_plot(fullname);
    if (err) {
	g_free(fullname);
	return;
    }

    err = gnuplot_make_image(fullname);
#if WINDEBUG
    fprintf(stderr, " gnuplot_make_image returned %d\n", err);
#endif

    if (err) {
	/* display the bad plot file */
	view_file(fullname, 0, 0, 78, 350, VIEW_FILE);
    } else {
	err = gnuplot_show_png(fullname, title, session_ptr, NULL);
    }

    g_free(fullname);

    if (err) {
	gui_errmsg(err);
    }
}

static int get_png_plot_bounds (const char *str, png_bounds *bounds)
{
    int ret = GRETL_PNG_OK;
    double bb[4];

    bounds->xleft = bounds->xright = 0;
    bounds->ybot = bounds->ytop = 0;

    gretl_push_c_numeric_locale();

    if (sscanf(str, "pixel_bounds: %lf %lf %lf %lf",
	       &bb[0], &bb[1], &bb[2], &bb[3]) != 4) {
	ret = GRETL_PNG_BAD_COMMENTS;
        fprintf(stderr, "bad bounds string = '%s'\n", str);
    } else {
        bounds->xleft  = nearbyint(bb[0]);
        bounds->xright = nearbyint(bb[1]);
        bounds->ybot   = nearbyint(bb[2]);
        bounds->ytop   = nearbyint(bb[3]);
    }

    gretl_pop_c_numeric_locale();

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
    int ret = GRETL_PNG_OK;
    char *p = str;

    /* ensure decimal dot */
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

static int get_png_size (char *str, png_bounds *bounds)
{
    int ret = GRETL_PNG_OK;
    int pw, ph, sc;

    bounds->width = bounds->height = 0;

    if (sscanf(str, "term_size: %d %d %d", &pw, &ph, &sc) != 3 ||
	pw <= 0 || ph <= 0 || sc <= 0) {
	ret = GRETL_PNG_BAD_COMMENTS;
    } else {
	pw /= sc; ph /= sc;
	if (pw % 2) pw++;
	if (ph % 2) ph++;
#if 0
	fprintf(stderr, "Got size: %d x %d\n", pw, ph);
#endif
	bounds->width = pw;
	bounds->height = ph;
    }

    return ret;
}

static int get_png_bounds_info (png_bounds *bounds)
{
    FILE *fp;
    char bbname[MAXLEN], line[128];
    int plot_ret = -1, data_ret = -1;
    int ret = GRETL_PNG_OK;

    gretl_build_path(bbname, gretl_dotdir(), "gretltmp.png.bounds", NULL);
    fp = gretl_fopen(bbname, "rb");

    if (fp == NULL) {
	fprintf(stderr, "couldn't open %s\n", bbname);
	return GRETL_PNG_NO_COMMENTS;
    }

    /* bounding box of plot */
    if (fgets(line, sizeof line, fp) == NULL) {
	fprintf(stderr, "bounds file: couldn't get plot dims\n");
	plot_ret = GRETL_PNG_NO_COMMENTS;
    } else {
	plot_ret = get_png_plot_bounds(line, bounds);
    }

    /* data ranges */
    if (fgets(line, sizeof line, fp) == NULL) {
	fprintf(stderr, "bounds file: couldn't get data dims\n");
	data_ret = GRETL_PNG_NO_COMMENTS;
    } else {
	data_ret = get_png_data_bounds(line, bounds);
    }

    /* overall size of plot */
    if (fgets(line, sizeof line, fp) == NULL) {
	fprintf(stderr, "bounds file: couldn't get size info\n");
    } else {
	get_png_size(line, bounds);
    }

    if (plot_ret == GRETL_PNG_NO_COORDS && data_ret == GRETL_PNG_NO_COORDS) {
	/* comments were present and correct, but all zero */
	fprintf(stderr, "bounds file: no actual coordinates\n");
	ret = GRETL_PNG_NO_COORDS;
    } else if (plot_ret != GRETL_PNG_OK || data_ret != GRETL_PNG_OK) {
	/* one or both set of coordinates bad or missing */
	if (plot_ret >= 0 || data_ret >= 0) {
	    fprintf(stderr, "bounds file: bad data: plot_ret %d, data_ret %d\n",
                    plot_ret, data_ret);
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

#ifdef OS_OSX
    /* fallback for XQuartz: may not be in PATH */
    strcpy(s, "/opt/X11/bin/xterm");
    if (gretl_file_exists(s)) {
	return 0;
    } else {
	*s = '\0';
    }
#endif

    errbox(_("Couldn't find a usable terminal program"));

    return 1;
}

#endif /* !G_OS_WIN32 */

#ifdef OS_OSX

static void mac_do_gp_script (const char *plotfile)
{
    gchar *buf = NULL;
    gsize sz = 0;

    if (g_file_get_contents(plotfile, &buf, &sz, NULL)) {
	run_gnuplot_script(buf, NULL);
	g_free(buf);
    }
}

#endif

/* Callback for "Gnuplot" item in Tools menu: open a
   gnuplot session and let the user do whatever. In
   this case we need a controlling terminal window if
   we're not on MS Windows.
*/

void launch_gnuplot_interactive (void)
{
#if defined(G_OS_WIN32)
    win32_run_async(gretl_gnuplot_path(), NULL);
#elif defined(OS_OSX)
    const char *gppath = gretl_gnuplot_path();
    gchar *gpline;

# ifdef PKGBUILD
    /* call driver script to set environment correctly -- and
       in addition prepend a full path spec if necessary
    */
    if (g_path_is_absolute(gppath)) {
	gpline = g_strdup_printf("open -a Terminal.app \"%s.sh\"",
				 gppath);
    } else {
	gpline = g_strdup_printf("open -a Terminal.app \"%s%s.sh\"",
				 gretl_bindir(), gppath);
    }
# else
    gpline = g_strdup_printf("open -a Terminal.app \"%s\"", gppath);
# endif
    system(gpline);
    g_free(gpline);
#else /* neither WIN32 nor MAC */
    char term[32];
    int err;

    err = get_terminal(term);

    if (!err) {
	const char *gp = gretl_gnuplot_path();
	GError *error = NULL;
	gchar *argv[6];

# ifdef OS_OSX
	char *altgp = g_strdup_printf("%s.sh", gp);

	if (gretl_file_exists(altgp)) {
	    gp = altgp;
	} else {
	    g_free(altgp);
	    altgp = NULL;
	}
# endif

	if (strstr(term, "gnome")) {
	    /* gnome-terminal */
	    argv[0] = term;
	    argv[1] = "--title=\"gnuplot: type q to quit\"";
	    argv[2] = "-x";
	    argv[3] = (char *) gp;
	    argv[4] = NULL;
	} else {
	    /* xterm, rxvt, kterm */
	    argv[0] = term;
	    argv[1] = "-title";
	    argv[2] = "gnuplot: type q to quit";
	    argv[3] = "-e";
	    argv[4] = (char *) gp;
	    argv[5] = NULL;
	}

	g_spawn_async(NULL, /* working dir */
		      argv,
		      NULL, /* env */
		      G_SPAWN_SEARCH_PATH,
		      NULL, /* child_setup */
		      NULL, /* user_data */
		      NULL, /* child_pid ptr */
		      &error);

	if (error != NULL) {
	    errbox(error->message);
	    g_error_free(error);
	}

# ifdef OS_OSX
	g_free(altgp);
# endif
    }
#endif /* !(G_OS_WIN32 or MAC) */
}

void gnuplot_view_3d (const char *plotfile)
{
#if defined(G_OS_WIN32)
    win32_run_async(gretl_gnuplot_path(), plotfile);
#elif defined(OS_OSX) && !defined(GNUPLOT3D)
    mac_do_gp_script(plotfile);
#else
    real_send_to_gp(plotfile, 0);
#endif /* !(G_OS_WIN32 or MAC) */
}
