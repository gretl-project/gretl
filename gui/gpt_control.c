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

#ifdef __APPLE__
# include "macos_open.h"
#endif

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk/gdkkeysyms.h>
#include <errno.h>

#define GPDEBUG 0
#define WINDEBUG 0
#define POINTS_DEBUG 0
#define COLLDEBUG 0

/* pager for plot collection: use spin button? */
#define SPIN_PAGER 0

enum {
    PLOT_ZOOMED         = 1 << 0,
    PLOT_ZOOMING        = 1 << 1,
    PLOT_PNG_COORDS     = 1 << 2,
    PLOT_HAS_XRANGE     = 1 << 3,
    PLOT_HAS_YRANGE     = 1 << 4,
    PLOT_DONT_ZOOM      = 1 << 5,
    PLOT_DONT_EDIT      = 1 << 6,
    PLOT_DONT_MOUSE     = 1 << 7,
    PLOT_POSITIONING    = 1 << 8,
    PLOT_CURSOR_LABEL   = 1 << 9,
    PLOT_TERM_HIDDEN    = 1 << 10
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
    GtkToolItem *scaler;  /* scale selector combo */
    GdkWindow *window;
    GdkPixbuf *pbuf;
    cairo_t *cr;
#if GTK_MAJOR_VERSION >= 3
    cairo_surface_t *cs;
#else
    GdkPixmap *pixmap;
    GdkPixbuf *savebuf;
#endif
    SESSION_GRAPH *grf;
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
static void plot_do_rescale (png_plot *plot, double scale);
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
    return plot->grf != NULL;
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

static void sensitize_scale_selector (png_plot *plot)
{
    gboolean blocked = plot_is_zoomed(plot) ||
        plot_not_editable(plot);

    gtk_widget_set_sensitive(GTK_WIDGET(plot->scaler), !blocked);
}

static void graph_rescale_callback (GtkComboBox *cb, png_plot *plot)
{
    gchar *s = combo_box_get_active_text(cb);
    double scale = atof(s);

    plot_do_rescale(plot, scale);
    g_free(s);
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

static GtkToolItem *plot_scale_combo (png_plot *plot)
{
    GtkToolItem *item = gtk_tool_item_new();
    GtkWidget *combo = gtk_combo_box_text_new();
    double ps = min_graph_scale();
    char ps_str[16];
    int i = 0;

    for (i=0; !na(ps); i++) {
        sprintf(ps_str, "%.1f", ps);
        combo_box_append_text(combo, ps_str);
        if (ps == plot->spec->scale) {
            gtk_combo_box_set_active(GTK_COMBO_BOX(combo), i);
        }
        ps = next_graph_scale(ps, 1);
    }

    g_signal_connect(G_OBJECT(combo), "changed",
                     G_CALLBACK(graph_rescale_callback),
                     plot);
    gtk_container_add(GTK_CONTAINER(item), combo);
    gtk_widget_show(combo);

    return item;
}

static GretlToolItem plotbar_items[] = {
    { N_("Menu"),        GRETL_STOCK_MENU,    G_CALLBACK(graph_popup_callback), 0 },
    { N_("Windows"),     GRETL_STOCK_WINLIST, G_CALLBACK(plot_winlist_popup), 0 }
};

static void add_graph_toolbar (GtkWidget *hbox, png_plot *plot)
{
    GtkWidget *tbar;
    GtkToolItem *item;

    plot->toolbar = tbar = gretl_toolbar_new(NULL);

    /* hamburger menu */
    gretl_toolbar_insert(tbar, &plotbar_items[0], plotbar_items[0].func, plot, -1);
    /* plot scale selector */
    plot->scaler = item = plot_scale_combo(plot);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_tool_item_set_tooltip_text(item, _("scale"));
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, -1);
    sensitize_scale_selector(plot);
    /* window list button */
    gretl_toolbar_insert(tbar, &plotbar_items[1], plotbar_items[1].func, plot, -1);

    gtk_box_pack_start(GTK_BOX(hbox), tbar, FALSE, FALSE, 5);
}

/* end apparatus for graph toolbar */

/* Provide the data coordinates for a gretl/gnuplot graph, if they are
   all positive, otherwise return non-zero. Called by the plot editor.
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

static int is_plot_line (const char *s)
{
    return !strncmp(s, "plot", 4) || !strncmp(s, "splot", 5);
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
        } else if (is_plot_line(line)) {
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

    snprintf(pltname, sizeof pltname, "%sgptout.tmp", gretl_dotdir());

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
#elif defined(__APPLE__)
    macos_open_file(pdfname);
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
#elif defined(__APPLE__)
    macos_open_file(grfname);
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

/* Check whether (a) we might want to prepend a "set term" statement
   (not if the user has already given one) and (b) whether we might want
   to append "pause mouse close" (not if a batch-type terminal has been
   specified).
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

/* dump_plot_buffer: this is used when we're taking the material from an
   editor window containing gnuplot commands, and either (a) sending it
   to gnuplot for execution, or (b) saving it to a "user file".  In the
   saving-to-file case @addpause will be 0.

   This function handles the addition of "pause mouse close", if
   @addpause is non-zero.
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

void gnuplot_view_session_graph (const char *fname)
{
    gchar *buf = NULL;

    if (g_file_get_contents(fname, &buf, NULL, NULL)) {
        gchar *tmpfile = gretl_make_dotpath("showtmp.gp");
        int err;

        err = dump_plot_buffer(buf, tmpfile, 1, fname);
        if (!err) {
            err = real_send_to_gp(tmpfile, 1);
        }
        g_free(tmpfile);
        g_free(buf);
    }
}

#ifdef G_OS_WIN32

/* Common code for sending an EMF file to the clipboard, or printing an
   EMF, on MS Windows.
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

/* Here we're just posting the information that an image should be
   available for pasting (and also recording if the user wanted it to be
   monochrome).
*/

static void set_plot_for_copy (png_plot *plot)
{
    copyplot = plot;
    cb_image_mono = plot->spec->flags & GPT_MONO ? 1 : 0;
    flag_image_available();
}

/* Here we're responding to a request to paste the plot advertised
   above; @target tells us which of the posted formats the application
   has selected.
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

/* Read plotspec struct from gnuplot command file: the parsing is
   done by plotspec_read_file() in libgretl; here we just pass in the
   current dataset and transfer the extra information to @plot.
*/

static int read_plotspec_from_file (png_plot *plot)
{
    GPT_SPEC *spec = plot->spec;
    PlotReadInfo info;
    int err;

    /* check: are we already done? */
    if (PLOTSPEC_DETAILS_IN_MEMORY(spec)) {
#if GPDEBUG
        fprintf(stderr, "read_plotspec_from_file: info already in memory\n");
#endif
        return 0;
    }

    info.pd = plot->pd;
    info.pixel_width = plot->pixel_width;
    info.pixel_height = plot->pixel_height;

    gretl_error_clear();
    err = plotspec_read_file(spec, dataset, &info);

    plot->pd = info.pd;
    plot->pixel_width = info.pixel_width;
    plot->pixel_height = info.pixel_height;

    if ((err == E_DATA || err == E_FOPEN) && *gretl_errmsg_get() != '\0') {
        errbox(gretl_errmsg_get());
    }

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

/* Note that the screen coordinates as of the last mouse button press
   are recorded in plot->screen_x0 and plot->screen_y0. These represent
   the constant corner of the box that should be drawn.
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
   surface, attached to @plot as plot->cs. This will be used as the
   "backing store" for re-draws of the plot.
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
	size_t prefix_len = (p - line) + 7;
	size_t used;

	if (prefix_len >= sizeof tmp) {
	    return E_DATA;
	}
	memcpy(tmp, line, prefix_len);
	tmp[prefix_len] = '\0';
	used = prefix_len;

	if (used + strlen(fstr) >= sizeof tmp) {
	    return E_DATA;
	}
	strcpy(tmp + used, fstr); /* replacement font string */
	used += strlen(fstr);

	p = strchr(p + 7, '"'); /* closing quote of original font string */
	if (p != NULL) {
	    if (used + strlen(p) >= sizeof tmp) {
		return E_DATA;
	    }
	    strcpy(tmp + used, p);
	    /* caller's @line buffer is the same size as @tmp,
	       declared char[256] in activate_plot_font_choice
	    */
	    strcpy(line, tmp);
	}
    }

    return (p == NULL)? E_DATA : 0;
}

/* Filter the original font spec out of the plot file, substituting the
   user's choice from a font selection dialog, then recreate the
   on-screen PNG. We do this for plots for which we can't offer the
   full-blown graph editor.
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
        view_buffer(prn, 84, 240, _("gretl: boxplot data"), PRINT, NULL);
    }
}

static int add_to_session_callback (png_plot *plot)
{
    char fullname[FILENAME_MAX] = {0};
    SESSION_GRAPH *grf = NULL;
    int err, type;

    type = (plot->spec->code == PLOT_BOXPLOTS)? GRETL_OBJ_PLOT :
        GRETL_OBJ_GRAPH;

    err = gui_add_graph_to_session(plot->spec->fname, fullname,
                                   type, &grf);

    if (!err) {
        remove_png_term_from_plot(fullname, plot->spec);
        plot->status |= PLOT_TERM_HIDDEN;
        plot->grf = grf;
    }

    return err;
}

static void call_view_plot_commands (png_plot *plot)
{
    int err = 0;

    if (plot->grf == NULL) {
        err = add_to_session_callback(plot);
    }
    if (!err) {
        view_plot_commands(plot->grf);
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

static void plot_do_rescale (png_plot *plot, double scale)
{
    int err = 0;

    plot = plot_get_current(plot);

    if (plot->spec->scale == scale) {
        /* no-op */
        return;
    }

    err = real_plot_rescale(plot, scale);

    if (!err) {
        repaint_png(plot, PNG_REDISPLAY);
        if (in_collection(plot)) {
            rescale_siblings(plot, scale);
        }
    }
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
    } else if (!strcmp(item_string, _("View gnuplot commands"))) {
        call_view_plot_commands(plot);
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
                               (p->spec->flags & GPT_PRINT_MARKERS) &&  \
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
        N_("View gnuplot commands"),
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

/* Preparation for redisplaying graph: here we handle the case where
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

/* On replacing "full" view of a plot with zoomed view: sync the zoomed
   x and y ranges to both the plot and its attached GPT_SPEC pointer.
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

/* With a zoomed version of the current plot in place, replace the full
   version with the zoom.
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
    sensitize_scale_selector(plot);

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

#ifdef __APPLE__
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
#ifdef __APPLE__
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
    char fname[FILENAME_MAX];
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
    char pngname[FILENAME_MAX];
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
            sensitize_scale_selector(plot);
        } else if (view == PNG_UNZOOM) {
            plot->status ^= PLOT_ZOOMED;
            sensitize_scale_selector(plot);
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
    plot->grf = NULL;
    plot->spec = NULL;
    plot->editor = NULL;
    plot->window = NULL;
    plot->toolbar = NULL;
    plot->scaler = NULL;

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

/* @fname is the name of the file containing the plot commands */

static int gnuplot_show_png (const char *fname,
                             const char *name,
                             SESSION_GRAPH *grf,
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
            fname, (grf != NULL));
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

    if (strlen(fname) >= sizeof plot->spec->fname) {
	plotspec_destroy(plot->spec);
	free(plot);
	gretl_errmsg_set(_("Filename is too long"));
	return E_DATA;
    }
    strcpy(plot->spec->fname, fname);

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
        if (grf != NULL) {
            plot->grf = grf;
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

    if (strlen(fname) >= sizeof plot->spec->fname) {
	plotspec_destroy(plot->spec);
	free(plot);
	gretl_errmsg_set(_("Filename is too long"));
	return E_DATA;
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

/* Called on a newly created PNG graph; note that the filename returned
   by gretl_plotfile() is the input file (set of gnuplot commands).
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

/* @fname is the name of a plot command file from the current session
   and @title is its display name; @grf is a pointer to the session
   object.
*/

void display_session_graph (const char *fname,
                            const char *title,
                            SESSION_GRAPH *grf)
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
        err = gnuplot_show_png(fullname, title, grf, NULL);
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
    char bbname[FILENAME_MAX], line[128];
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

#ifdef __APPLE__
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

#ifdef __APPLE__

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

/* Callback for "Gnuplot" item in Tools menu: open a gnuplot session and
   let the user do whatever. In this case we need a controlling terminal
   window if we're not on MS Windows.
*/

void launch_gnuplot_interactive (void)
{
#if defined(G_OS_WIN32)
    win32_run_async(gretl_gnuplot_path(), NULL);
#elif defined(__APPLE__)
    const char *gppath = gretl_gnuplot_path();
    gchar *script;
    gchar *argv[5];
    GError *error = NULL;

# ifdef PKGBUILD
    /* call driver script to set environment correctly -- and
       in addition prepend a full path spec if necessary
    */
    script = g_path_is_absolute(gppath) ?
	g_strdup_printf("%s.sh", gppath) :
	g_strdup_printf("%s%s.sh", gretl_bindir(), gppath);
# else
    script = g_strdup(gppath);
# endif
    argv[0] = "open";
    argv[1] = "-a";
    argv[2] = "Terminal.app";
    argv[3] = script;
    argv[4] = NULL;

    g_spawn_async(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
		  NULL, NULL, NULL, &error);
    if (error != NULL) {
	errbox(error->message);
	g_error_free(error);
    }
    g_free(script);
#else /* neither WIN32 nor MAC */
    char term[32];
    int err;

    err = get_terminal(term);

    if (!err) {
        const char *gp = gretl_gnuplot_path();
        GError *error = NULL;
        gchar *argv[6];

# ifdef __APPLE__
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

# ifdef __APPLE__
        g_free(altgp);
# endif
    }
#endif /* !(G_OS_WIN32 or MAC) */
}

void gnuplot_view_3d (const char *plotfile)
{
#if defined(G_OS_WIN32)
    win32_run_async(gretl_gnuplot_path(), plotfile);
#elif defined(__APPLE__) && !defined(GNUPLOT3D)
    mac_do_gp_script(plotfile);
#else
    real_send_to_gp(plotfile, 0);
#endif /* !(G_OS_WIN32 or MAC) */
}
