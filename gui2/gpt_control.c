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

#include "boxplots.h"

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
    PLOT_POSITIONING    = 1 << 10
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

#define plot_is_saved(p)        (p->status & PLOT_SAVED)
#define plot_is_zoomed(p)       (p->status & PLOT_ZOOMED)
#define plot_is_zooming(p)      (p->status & PLOT_ZOOMING)
#define plot_has_png_coords(p)  (p->status & PLOT_PNG_COORDS)
#define plot_has_xrange(p)      (p->status & PLOT_HAS_XRANGE)
#define plot_has_yrange(p)      (p->status & PLOT_HAS_YRANGE)
#define plot_not_zoomable(p)    (p->status & PLOT_DONT_ZOOM)
#define plot_not_editable(p)    (p->status & PLOT_DONT_EDIT)
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
#define plot_is_roots(p)        (p->spec->code == PLOT_VAR_ROOTS)

#define plot_has_regression_list(p) (p->spec->reglist != NULL)

#define labels_frozen(p)        (p->spec->flags & GPT_PRINT_MARKERS)
#define cant_do_labels(p)       (p->status & PLOT_NO_MARKERS)

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
    GtkWidget *labelpos_entry;
    GtkWidget *editor;
    GdkPixmap *pixmap;
    GdkGC *invert_gc;
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
    int screen_xmin, screen_ymin;
    unsigned long status; 
    unsigned char format;
};

static int render_pngfile (png_plot *plot, int view);
static int repaint_png (png_plot *plot, int view);
static void create_selection_gc (png_plot *plot);
static int get_plot_ranges (png_plot *plot);
static void graph_display_pdf (GPT_SPEC *spec);
#ifdef G_OS_WIN32
static void win32_process_graph (GPT_SPEC *spec, int dest);
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
};

typedef struct linestyle_ linestyle;

struct linestyle_ {
    char rgb[8];
    int type;
};

#define MAX_STYLES 6

static int get_png_bounds_info (png_bounds *bounds);

#define PLOTSPEC_DETAILS_IN_MEMORY(s) (s->data != NULL)

static void terminate_plot_positioning (png_plot *plot)
{
    if (!(plot->status & PLOT_POSITIONING)) {
	return;
    }
	    
    plot->status ^= PLOT_POSITIONING;
    plot->labelpos_entry = NULL;
    gdk_window_set_cursor(plot->canvas->window, NULL);
    gtk_statusbar_pop(GTK_STATUSBAR(plot->statusbar), plot->cid);
    if (plot->editor != NULL) {
	gtk_window_present(GTK_WINDOW(plot->editor));
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

void set_plot_has_y2_axis (png_plot *plot, gboolean s)
{
    if (s == TRUE) {
	plot->format |= PLOT_Y2AXIS;
    } else {
	plot->format &= ~PLOT_Y2AXIS;
    }
}

void plot_label_position_click (GtkWidget *w, png_plot *plot)
{
    if (plot != NULL) {
	GtkWidget *entry;
	GdkCursor* cursor;

	cursor = gdk_cursor_new(GDK_CROSSHAIR);
	gdk_window_set_cursor(plot->canvas->window, cursor);
	gdk_cursor_unref(cursor);
	entry = g_object_get_data(G_OBJECT(w), "labelpos_entry");
	plot->labelpos_entry = entry;
	plot->status |= PLOT_POSITIONING;
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid, 
			   _(" Click to set label position"));
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
    char restore_line[MAXLEN];
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
    
    *restore_line = '\0';

    if (action == ADD_PNG) {
	/* see if there's already a png term setting, possibly commented
	   out, that can be reused */
	while (fgets(fline, sizeof fline, fsrc)) {
	    if (is_png_term_line(fline) && *restore_line == '\0') {
		strcat(restore_line, fline);
	    } else if (commented_term_line(fline) && *restore_line == '\0') {
		strcat(restore_line, fline + 2);
	    } else if (strstr(fline, "letterbox")) {
		flags = GPT_LETTERBOX;
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

	if (*restore_line) {
	    fputs(restore_line, ftmp);
	} else {
	    int ptype = (spec != NULL)? spec->code : PLOT_REGULAR;
	    const char *pline;

	    if (spec != NULL) {
		flags = spec->flags;
	    }
	    pline = get_gretl_png_term_line(ptype, flags);
	    fprintf(ftmp, "%s\n", pline);
	}	    
	fprintf(ftmp, "set output '%sgretltmp.png'\n", gretl_dotdir());
    }

    /* now for the body of the plot file */

    if (action == ADD_PNG) {
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
	if (gnuplot_has_bbox()) {
	    print_plot_bounding_box_request(ftmp);
	}
    } else {
	/* we're removing the png term line */
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

static int add_png_term_to_plotfile (const char *fname)
{
    return add_or_remove_png_term(fname, ADD_PNG, NULL);
}

static int remove_png_term_from_plotfile (const char *fname, GPT_SPEC *spec)
{
    return add_or_remove_png_term(fname, REMOVE_PNG, spec);
}

/* public because called from session.c when editing plot commands */

int remove_png_term_from_plotfile_by_name (const char *fname)
{
    return add_or_remove_png_term(fname, REMOVE_PNG, NULL);
}

static void mark_plot_as_saved (GPT_SPEC *spec)
{
    png_plot *plot = (png_plot *) spec->ptr;

    plot->status |= PLOT_SAVED;
}

static int gnuplot_png_init (GPT_SPEC *spec, FILE **fpp)
{
    *fpp = gretl_fopen(spec->fname, "w");

    if (*fpp == NULL) {
	file_write_errbox(spec->fname);
	return 1;
    }

    fprintf(*fpp, "%s\n", get_gretl_png_term_line(spec->code, spec->flags));
    fprintf(*fpp, "set output '%sgretltmp.png'\n", gretl_dotdir());

    return 0;
}

int gp_term_code (gpointer p)
{
    GPT_SPEC *spec = (GPT_SPEC *) p;

    return spec->termtype;
}

#define PDF_CAIRO_STRING "set term pdfcairo font \"sans,5\""

static void 
get_full_term_string (const GPT_SPEC *spec, char *termstr) 
{
    int mono = (spec->flags & GPT_MONO);
    char fontname[64] = {0};
    int psz = 0;

    if (spec->termtype == GP_TERM_EPS || spec->termtype == GP_TERM_PDF) {
	const char *pngfont = gretl_png_font();

	if (*pngfont != '\0') {
	    split_graph_fontspec(pngfont, fontname, &psz);
	}
    }

    if (spec->termtype == GP_TERM_EPS) {
	if (mono) {
	    strcpy(termstr, "set term postscript eps mono"); 
	} else {
	    strcpy(termstr, "set term postscript eps color solid");
	} 
    } else if (spec->termtype == GP_TERM_PDF) {
	if (gnuplot_pdf_terminal() == GP_PDF_CAIRO) {
	    if (*fontname != '\0') {
		sprintf(termstr, "set term pdfcairo font \"%s\"\n",
			fontname);
	    } else {
		strcpy(termstr, "set term pdfcairo");
	    }
	} else {
	    strcpy(termstr, "set term pdf");
	}
    } else if (spec->termtype == GP_TERM_FIG) {
	strcpy(termstr, "set term fig");
    } else if (spec->termtype == GP_TERM_PNG) { 
	const char *png_str = 
	    get_gretl_png_term_line(spec->code, spec->flags);

	strcpy(termstr, png_str); 
    } else if (spec->termtype == GP_TERM_EMF) {
	const char *emf_str = 
	    get_gretl_emf_term_line(spec->code, !mono);

	strcpy(termstr, emf_str);
    } 
}

static char *gp_contd_string (char *s)
{
    char *p = strstr(s, ", \\");
    int n = 0;

    if (p != NULL) {
	n = 3;
    } else {
	p = strstr(s, ",\\");
	if (p != NULL) {
	    n = 2;
	}
    }

    if (p != NULL) {
	/* ensure we've really got '\' at end of line */
	char c = *(p+n);

	if (c != '\0' && c != '\n' && c != '\r') {
	    p = NULL;
	}
    }

    return p;
}

static char *get_insert_point (char *s, char *p)
{
    if (p == NULL) {
	p = s + strlen(s) - 1;
    }

    if (p - s > 0 && *(p-1) == ' ') {
	p--;
    }

    return p;
}

/* old gnuplot: can't do "rgb" line spec; we'll do what we can,
   namely switch line type 2 for 3
*/

static void maybe_recolor_line (char *s, int lnum)
{
    const gretlRGB *color = get_graph_color(lnum - 1);

    if (color != NULL) {
	char *contd = gp_contd_string(s);
	char cstr[8];

	print_rgb_hash(cstr, color);

#if GPDEBUG
	fprintf(stderr, "lnum=%d, cstr='%s'\n", lnum, cstr);
	fprintf(stderr, "s='%s'\n", s);
#endif
    
	if (lnum == 2 && strcmp(cstr, "#00ff00") && !strstr(s, " lt ")) {
	    char *p = get_insert_point(s, contd);

	    *p = '\0';
	    strcpy(p, " lt 3");
	    if (contd != NULL) {
		strcat(s, ", \\\n");
	    } else {
		strcat(s, "\n");
	    }
	} 
    } 
}

static void dataline_check (char *s, int *d)
{
    if (!strncmp(s, "plot \\", 6)) {
	*d = 0;
    } else {
	if (!strncmp(s, "plot ", 5)) {
	    *d = 0;
	}
	if (*d == 0 && !gp_contd_string(s)) {
	    *d = 1;
	}
    }
}

/* for postscript output, e.g. in Latin-2, or EMF output in CP125X */

static int maybe_recode_gp_line (char *s, int ttype, FILE *fp)
{
    int err = 0;

    if (!gretl_is_ascii(s) && g_utf8_validate(s, -1, NULL)) {
	char *tmp;
	
	if (ttype == GP_TERM_EMF) {
	    tmp = utf8_to_cp(s);
	} else {
	    tmp = utf8_to_latin(s);
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

/* check for non-ASCII strings in plot file: these may
   require special treatment */

static int non_ascii_gp_file (FILE *fp)
{
    char pline[512];
    int dataline = -1;
    int ret = 0;

    while (fgets(pline, sizeof pline, fp) && dataline <= 0) {
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
	dataline_check(pline, &dataline);
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

void filter_gnuplot_file (int ttype, int latin, int mono, int recolor, 
			  FILE *fpin, FILE *fpout)
{
    char pline[512];
    int dataline = -1;
    int lnum = -1;
    int err = 0;

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

	if (recolor) {
	    if (!strncmp(pline, "plot ", 5)) {
		lnum = 0;
	    } else if (lnum >= 0) {
		lnum++;
	    }
	    if (lnum > 0 && dataline <= 0) {
		maybe_recolor_line(pline, lnum);
	    }
	}

	if (latin && dataline <= 0 && *pline != '#') {
	    err += maybe_recode_gp_line(pline, ttype, fpout);
	    if (err == 1) {
		gui_errmsg(err);
	    }
	} else {
	    fputs(pline, fpout);
	} 

	dataline_check(pline, &dataline);
    }
}

/* for non-UTF-8 plot formats: print a "set encoding" string
   if appropriate, but only if gnuplot won't choke on it.
*/

static void maybe_print_gp_encoding (int ttype, int latin, FILE *fp)
{
    if (ttype == GP_TERM_EMF) {
	if (latin == 2 && gnuplot_has_cp1250()) {
	    fputs("set encoding cp1250\n", fp);
	} else if (latin == 9 && gnuplot_has_cp1254()) {
	    fputs("set encoding cp1254\n", fp);
	}
    } else {
	if (latin != 1 && latin != 2 && latin != 15 && latin != 9) {
	    /* unsupported by gnuplot */
	    latin = 0;
	}
	if (latin == 9 && !gnuplot_has_latin5()) {
	    /* Turkish not supported */
	    latin = 0;
	}
	if (latin) {
	    fprintf(fp, "set encoding iso_8859_%d\n", latin);
	}
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
    int recolor = 0;
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
	} else if (gnuplot_has_utf8()) {
	    fputs("set encoding utf8\n", fpout);
	}
    }

    /* FIXME font name */

    if (outtarg != NULL && *outtarg != '\0') {
	fprintf(fpout, "%s\n", setterm);
	fprintf(fpout, "set output '%s'\n", outtarg);
    }	

    if (!mono && (ttype == GP_TERM_EPS || ttype == GP_TERM_PDF)
	&& !gnuplot_has_rgb()) {
	recolor = 1;
    }

    filter_gnuplot_file(ttype, latin, mono, recolor, fpin, fpout);

    fclose(fpin);
    fclose(fpout);

    return err;
}

void save_graph_to_file (gpointer data, const char *fname)
{
    GPT_SPEC *spec = (GPT_SPEC *) data;
    char setterm[MAXLEN];
    char pltname[MAXLEN];
    int err = 0;

    get_full_term_string(spec, setterm);

    build_path(pltname, gretl_dotdir(), "gptout.tmp", NULL);

    err = revise_plot_file(spec, pltname, fname, setterm);

    if (!err) {
	gchar *plotcmd;

	plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
				  gretl_gnuplot_path(), 
				  pltname);
	err = gretl_spawn(plotcmd);
	/* gretl_remove(pltname); */
	g_free(plotcmd);
	if (err) {
	    gui_errmsg(err);
	} 
    }
}

#define GRETL_PDF_TMP "gretltmp.pdf"

static void graph_display_pdf (GPT_SPEC *spec)
{
    char pdfname[FILENAME_MAX];
    char plttmp[FILENAME_MAX];
    static char setterm[64];
    gchar *plotcmd;
    int err = 0;

    spec->termtype = GP_TERM_PDF;

    if (*setterm == '\0') {
	if (gnuplot_pdf_terminal() == GP_PDF_CAIRO) {
	    fprintf(stderr, "gnuplot: using pdfcairo driver\n");
	    strcpy(setterm, PDF_CAIRO_STRING);
	} else {
	    strcpy(setterm, "set term pdf");
	}
    }

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
#elif defined(OSX_BUILD)
    osx_open_file(pdfname);
#else
    gretl_fork("viewpdf", pdfname);
#endif
}

/* dump_plot_buffer: this is used when we're taking the material from
   an editor window containing gnuplot commands, and either (a)
   sending it to gnuplot for execution, or (b) saving it to a "user
   file".  

   There's a question over what we should do with non-ascii strings in
   the plot file.  These will be in UTF-8 in the GTK editor window.
   That's OK if gnuplot supports UTF-8, but otherwise it seems that
   the best thing is to determine the character set for the current
   locale (using g_get_charset) and if it is not UTF-8, recode to the
   locale. 

   This function also handles the addition of "pause -1" on MS
   Windows, if @addpause is non-zero.
*/

int dump_plot_buffer (const char *buf, const char *fname,
		      int addpause)
{
    const gchar *cset;
    FILE *fp;
    int recode = 0;

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	file_write_errbox(fname);
	return E_FOPEN;
    }

    if (!g_get_charset(&cset) && !gnuplot_has_utf8()) {
	/* we're on a non-UTF-8 platform and we need to convert */
	recode = 1;
    } else if (!gnuplot_has_utf8()) {
	/* we're screwed -- what do we recode to? */
	fprintf(stderr, "Warning: gnuplot does not support UTF-8\n");
    }

    if (!recode && !addpause) {
	/* nice and simple! */
	fputs(buf, fp);
    } else {
	gchar *trbuf;
	char bufline[512];
	int gotpause = 0;
	int handled;

	bufgets_init(buf);

	while (bufgets(bufline, sizeof bufline, buf)) {
	    handled = 0;
	    if (recode) {
		trbuf = gp_locale_from_utf8(bufline);
		if (trbuf != NULL) {
		    fputs(trbuf, fp);
		    g_free(trbuf);
		    handled = 1;
		}
	    }
	    if (!handled) {
		fputs(bufline, fp);
	    }
	    if (addpause && strstr(bufline, "pause -1")) {
		gotpause = 1;
	    }
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
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
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

static int get_gpt_marker (const char *line, char *label)
{
    const char *p = strchr(line, '#');
    char format[6];

    if (p != NULL) {
	sprintf(format, "%%%ds", OBSLEN - 1);
	sscanf(p + 1, format, label);
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
                      p == PLOT_VAR_ROOTS || \
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

static int get_gpt_data (GPT_SPEC *spec, int do_markers, 
			 const char *buf)
{
    char s[MAXLEN];
    char *got;
    double *x[5] = { NULL };
    char test[5][32];
    int started_data_lines = 0;
    int i, j, t;
    int err = 0;

    spec->okobs = spec->nobs;

    gretl_push_c_numeric_locale();

    for (i=0; i<spec->n_lines && !err; i++) {
	int ncols = spec->lines[i].ncols;
	int okobs = spec->nobs;
	int offset = 1;

	if (ncols == 0) {
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

	for (t=0; t<spec->nobs; t++) {
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

	    for (j=offset; j<nf; j++) {
		if (test[j][0] == '?') {
		    x[j][t] = NADBL;
		    missing++;
		} else {
		    x[j][t] = atof(test[j]);
		}
	    }

	    if (missing) {
		okobs--;
	    }

	    if (i == 0 && do_markers) {
		get_gpt_marker(s, spec->markers[t]);
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

static int read_plot_format (const char *s, GPT_SPEC *spec)
{
    char axis, fmt[16];
    int n, err = 0;

    n = sscanf(s, "%c \"%15[^\"]", &axis, fmt);

    if (n < 2) {
	err = 1;
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

static int parse_gp_set_line (GPT_SPEC *spec, const char *s, 
			      linestyle *styles)
{
    char key[16] = {0};
    char val[MAXLEN] = {0};

    if (!strncmp(s, "set style line", 14)) {
	/* e.g. set style line 1 lc rgb "#ff0000 lt 6" */
	int n, idx = 0, lt = LT_NONE;
	char rgb[8];

	n = sscanf(s + 14, " %d lc rgb \"%7s\" lt %d", &idx, rgb, &lt);
	if (n >= 2 && idx > 0 && idx <= MAX_STYLES) {
	    strcpy(styles[idx-1].rgb, rgb);
	    styles[idx-1].type = (n == 3)? lt : LT_NONE;
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
    } else if (!strcmp(key, "xzeroaxis")) {
	spec->flags |= GPT_XZEROAXIS;
	return 0;
    } else if (!strcmp(key, "yzeroaxis")) {
	spec->flags |= GPT_YZEROAXIS;
	return 0;
    } else if (!strcmp(key, "noxtics")) {
	strcpy(spec->xtics, "none");
	return 0;
    } else if (!strcmp(key, "nokey")) {
	spec->keyspec = GP_KEY_NONE;
	return 0;
    } else if (!strcmp(key, "label")) {
	parse_label_line(spec, s);
	return 0;
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
	if (read_plot_format(val, spec)) {
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
    } else if (!strcmp(key, "key")) {
	spec->keyspec = gp_keypos_from_string(val);
    } else if (!strcmp(key, "xtics")) { 
	safecpy(spec->xtics, val, sizeof(spec->xtics) - 1);
    } else if (!strcmp(key, "mxtics")) { 
	safecpy(spec->mxtics, val, sizeof(spec->mxtics) - 1);
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
    } 

    return 0;
}

/* allocate markers for identifying particular data points */

static int allocate_plotspec_markers (GPT_SPEC *spec)
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
   determine the type of plot, and check whether there are any
   data-point markers along with the data.
*/

static int get_plot_nobs (const char *buf, PlotType *ptype, int *do_markers)
{
    int n = 0, started = -1;
    char line[MAXLEN], test[9];
    char *p;

    *ptype = PLOT_REGULAR;
    *do_markers = 0;

    while (bufgets(line, MAXLEN - 1, buf)) {

	if (*line == '#' && *ptype == PLOT_REGULAR) {
	    tailstrip(line);
	    *ptype = plot_type_from_string(line);
	}

	if (!strncmp(line, "plot", 4)) {
	    started = 0;
	}

	if (started == 0 && strchr(line, '\\') == NULL) {
	    started = 1;
	    continue;
	}

	if (started == 1) {
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

    return n;
}

static int grab_fit_coeffs (GPT_SPEC *spec, const char *s)
{
    int n, err = 0;

    if (spec->fit == PLOT_FIT_OLS) {
	spec->b_ols = gretl_column_vector_alloc(2);
	if (spec->b_ols == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_push_c_numeric_locale();
	    n = sscanf(s, "%lf + %lf", &spec->b_ols->val[0],
		       &spec->b_ols->val[1]);
	    gretl_pop_c_numeric_locale();
	    if (n != 2) {
		gretl_matrix_free(spec->b_ols);
		spec->b_ols = NULL;
		err = E_DATA;
	    }
	}
    } else if (spec->fit == PLOT_FIT_INVERSE) {
	spec->b_inv = gretl_column_vector_alloc(2);
	if (spec->b_inv == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_push_c_numeric_locale();
	    n = sscanf(s, "%lf + %lf", &spec->b_inv->val[0],
		       &spec->b_inv->val[1]);
	    gretl_pop_c_numeric_locale();
	    if (n != 2) {
		gretl_matrix_free(spec->b_inv);
		spec->b_inv = NULL;
		err = E_DATA;
	    }
	}
    } else if (spec->fit == PLOT_FIT_QUADRATIC) {
	spec->b_quad = gretl_column_vector_alloc(3);
	if (spec->b_quad == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_push_c_numeric_locale();
	    n = sscanf(s, "%lf + %lf*X + %lf", &spec->b_quad->val[0],
		       &spec->b_quad->val[1], &spec->b_quad->val[2]);
	    gretl_pop_c_numeric_locale();
	    if (n != 3) {
		gretl_matrix_free(spec->b_quad);
		spec->b_quad = NULL;
		err = E_DATA;
	    }
	}
    }	

    if (err) {
	spec->flags &= ~GPT_AUTO_FIT;
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

static int parse_gp_line_line (const char *s, GPT_SPEC *spec)
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
    } else {
	/* absence of "using" means the line plots a formula, not a
	   set of data columns 
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
	line->style = gp_style_from_string(tmp);
    } 

    if ((p = strstr(s, " lt "))) {
	sscanf(p + 4, "%d", &line->type);
    }

    if ((p = strstr(s, " pt "))) {
	sscanf(p + 4, "%d", &line->ptype);
    }

    if ((p = strstr(s, " lw "))) {
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
    int v = current_series_index(datainfo, vname);

    if (v >= 0 && !strcmp(datainfo->varname[v], vname)) {
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
	spec->lines[0].ncols == 2 &&
	!(spec->flags & (GPT_IMPULSES|GPT_LINES|GPT_RESIDS|GPT_TS))) {
	spec->fit = PLOT_FIT_NONE;
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
	err = 1;
    }

    /* and markers if any */
    if (!err && do_markers) {
	if (allocate_plotspec_markers(spec)) {
	    free(spec->data);
	    spec->data = NULL;
	    err = 1;
	}
    }

    /* Read the data (and perhaps markers) from the plot file */
    if (!err) {
	err = get_gpt_data(spec, do_markers, buf);
    }

#if GPDEBUG
    fprintf(stderr, "plot_get_data_and_markers:\n"
	    " spec->data = %p, spec->markers = %p, spec->n_markers = %d, err = %d\n",
	    spec->data, (void *) spec->markers, spec->n_markers, err);
#endif

    return err;
}

/* Get data markers from a "non-editable" plot and set the polar flag,
   if relevant.  The only case where we're able to use brush-on
   markers with non-editable plots is where we're showing numerical
   values for the roots in a VAR roots plot.
*/

static int uneditable_get_markers (GPT_SPEC *spec, const char *buf, int *polar)
{
    char line[256];
    long offset = 0;
    int gotit = 0;
    int err = 0;

    buf_rewind(buf);

    /* advance to the right line (with data plus markers) */
    while (bufgets(line, sizeof line, buf)) {
	if ((isdigit(*line) || *line == '-') && strchr(line, '#')) {
	    gotit = 1;
	    break;
	} else if (strstr(line, "set polar")) {
	    *polar = 1;
	}
	offset = buftell(buf);
    }

    if (!gotit) {
	return 1;
    } else {
	/* not found: back up */
	bufseek(buf, offset);
    }

    err = plotspec_add_line(spec);
    
    if (!err) {
	spec->lines[0].ncols = 2;
	err = plot_get_data_and_markers(spec, buf, 2, 1);
    }

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
    } else if (strstr(s, "inverse")) {
	return PLOT_FIT_INVERSE;
    } else if (strstr(s, "loess")) {
	return PLOT_FIT_LOESS;
    } else {
	return PLOT_FIT_NONE;
    }
}

static void linestyle_init (linestyle *ls)
{
    ls->rgb[0] = '\0';
    ls->type = LT_NONE;
}

#define plot_needs_obs(c) (c != PLOT_ELLIPSE && \
                           c != PLOT_PROB_DIST && \
                           c != PLOT_CURVE)

/* Read plotspec struct from gnuplot command file.  This is _not_ a
   general parser for gnuplot files; it is designed specifically for
   files auto-generated by gretl.
*/

static int read_plotspec_from_file (GPT_SPEC *spec, int *plot_pd, int *polar)
{
    linestyle styles[MAX_STYLES];
    int i, done;
    int do_markers = 0;
    int datacols = 0;
    int reglist[4] = {0};
    int *uservec = NULL;
    char gpline[MAXLEN];
    gchar *buf = NULL;
    char *got = NULL;
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
    err = gretl_file_get_contents(spec->fname, &buf);
    if (err) {
	gui_errmsg(err);
	return err;
    }

    bufgets_init(buf);

    /* get the number of data-points, plot type, and check for markers */
    spec->nobs = get_plot_nobs(buf, &spec->code, &do_markers);
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
	if (do_markers) {
	    uneditable_get_markers(spec, buf, polar);
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
	} else if (!strncmp(gpline, "# boxplots", 10)) {
	    continue;
	}

	if (sscanf(gpline, "# X = '%15[^\']' (%d)", vname, &v) == 2) {
	    if (plot_ols_var_ok(vname)) {
		reglist[2] = v;
	    }
	    continue;
	} else if (sscanf(gpline, "# Y = '%15[^\']' (%d)", vname, &v) == 2) {
	    if (reglist[2] > 0 && plot_ols_var_ok(vname)) {
		reglist[0] = 3;
		reglist[1] = v;
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

	if (strncmp(gpline, "set ", 4)) {
	    /* done reading "set" lines */
	    break;
	}

	if (parse_gp_set_line(spec, gpline, styles)) {
	    err = 1;
	    goto bailout;
	}
    }

    if (got == NULL) {
	err = 1;
	goto bailout;
    }

    for (i=0; i<4; i++) {
	if (spec->titles[i][0] != '\0') {
	    delchar('"', spec->titles[i]);
	}
    }

    /* then get the "plot" lines */
    if (strncmp(gpline, "plot ", 5) ||
	(strlen(gpline) < 10 && bufgets(gpline, MAXLEN - 1, buf) == NULL)) {	
	errbox(_("Failed to parse gnuplot file"));
	fprintf(stderr, "bad plotfile line: '%s'\n", gpline);
	err = 1;
	goto bailout;
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

	err = parse_gp_line_line(gpline, spec);

	if (err || done || (got = bufgets(gpline, MAXLEN - 1, buf)) == NULL) {
	    break;
	}
    }

    if (err || got == NULL) {
	err = 1;
	goto bailout;
    }

    /* determine total number of required data columns, etc. */
    for (i=0; i<spec->n_lines; i++) {
	if (uservec != NULL && in_gretl_list(uservec, i)) {
	    spec->lines[i].flags |= GP_LINE_USER;
	}
	if (i < MAX_STYLES) {
	    strcpy(spec->lines[i].rgb, styles[i].rgb);
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

    rx = (plot->screen_xmin < x)? plot->screen_xmin : x;
    ry = (plot->screen_ymin < y)? plot->screen_ymin : y;
    rw = x - plot->screen_xmin;
    rh = y - plot->screen_ymin;
    if (rw < 0) rw = -rw;
    if (rh < 0) rh = -rh;    

    /* draw one time to make the rectangle appear */
    gdk_draw_rectangle(plot->pixmap,
		       plot->invert_gc,
		       FALSE,
		       rx, ry, rw, rh);
    /* show the modified pixmap */
    gdk_draw_drawable(plot->canvas->window,
		      plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
		      plot->pixmap,
		      0, 0,
		      0, 0,
		      plot->pixel_width, plot->pixel_height);
    /* draw (invert) again to erase the rectangle */
    gdk_draw_rectangle(plot->pixmap,
		       plot->invert_gc,
		       FALSE,
		       rx, ry, rw, rh);
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

static void
write_label_to_plot (png_plot *plot, int i, gint x, gint y)
{
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

    if (plot->invert_gc == NULL) {
	create_selection_gc(plot);
    }

    context = gtk_widget_get_pango_context(plot->shell);
    pl = pango_layout_new(context);
    pango_layout_set_text(pl, label, -1);

    /* draw the label */
    gdk_draw_layout(plot->pixmap, plot->invert_gc, x, y, pl);

    /* show the modified pixmap */
    gdk_draw_drawable(plot->canvas->window,
		      plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
		      plot->pixmap,
		      0, 0,
		      0, 0,
		      plot->pixel_width, plot->pixel_height);

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
    double min_xdist, min_ydist;
    int best_match = -1;
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
	min_xdist = xrange = plot->zoom_xmax - plot->zoom_xmin;
	min_ydist = yrange = plot->zoom_ymax - plot->zoom_ymin;
    } else {
	min_xdist = xrange = plot->xmax - plot->xmin;
	min_ydist = yrange = plot->ymax - plot->ymin;
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

    /* try to find the best-matching data point */
    for (t=0; t<plot->spec->nobs; t++) {
	if (na(data_x[t]) || na(data_y[t])) {
	    continue;
	}
#if GPDEBUG > 2
	fprintf(stderr, "considering t=%d: x=%g, y=%g\n", t, data_x[t], data_y[t]);
#endif
	xdiff = fabs(data_x[t] - x);
	ydiff = fabs(data_y[t] - y);
	if (xdiff <= min_xdist && ydiff <= min_ydist) {
	    min_xdist = xdiff;
	    min_ydist = ydiff;
	    best_match = t;
	}
    }

    /* if the point is already labeled, skip */
    if (plot->spec->labeled[best_match]) {
	return TRUE;
    }

#if GPDEBUG > 2
    fprintf(stderr, " best_match=%d, with data_x[%d]=%g, data_y[%d]=%g\n", 
	    best_match, best_match, data_x[best_match], 
	    best_match, data_y[best_match]);
#endif

    /* if the match is good enough, show the label */
    if (best_match >= 0 && 
	min_xdist < TOLDIST * xrange &&
	min_ydist < TOLDIST * yrange) {
	write_label_to_plot(plot, best_match, pixel_x, pixel_y);
	/* flag the point as labeled already */
	plot->spec->labeled[best_match] = 1;
    }

    return TRUE;
}

#define float_fmt(i,x) ((i) && fabs(x) < 1.0e7)

static gint
motion_notify_event (GtkWidget *widget, GdkEventMotion *event, png_plot *plot)
{
    GdkModifierType state;
    gchar label[48], label_y[24];
    const char *xfmt = NULL;
    const char *yfmt = NULL;
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

    *label = 0;

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

	if (plot->pd == 4 || plot->pd == 12) {
	    x_to_date(data_x, plot->pd, label);
	} else if (xfmt != NULL) {
	    sprintf(label, xfmt, data_x);
	} else {
	    sprintf(label, (float_fmt(plot->xint, data_x))? "%7.0f" : 
		    "%#7.4g", data_x);
	}

	if (!na(data_y)) {
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
    if (plot->spec->flags & GPT_Y2AXIS) {
	plot->format |= PLOT_Y2AXIS;
    }
    if (plot->spec->flags & GPT_PRINT_MARKERS) {
	plot->format |= PLOT_MARKERS_UP;
    }
}

/* called from png plot popup menu */

static void start_editing_png_plot (png_plot *plot)
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

static gint color_popup_activated (GtkWidget *w, gpointer data)
{
    gchar *item = (gchar *) data;
    gpointer ptr = g_object_get_data(G_OBJECT(w), "plot");
    png_plot *plot = (png_plot *) ptr;
    GtkWidget *parent = (GTK_MENU(w->parent))->parent_menu_item;
    gchar *up_item = g_object_get_data(G_OBJECT(parent), "string");

    if (!strcmp(item, _("monochrome"))) {
	plot->spec->flags |= GPT_MONO;
    }

    if (!strcmp(up_item, _("Save as postscript (EPS)..."))) {
	plot->spec->termtype = GP_TERM_EPS;
	file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC, 
				  plot->spec, plot->shell);
    } else if (!strcmp(up_item, _("Save as Windows metafile (EMF)..."))) {
	plot->spec->termtype = GP_TERM_EMF;
	file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC, 
				  plot->spec, plot->shell);
    } 
#ifdef G_OS_WIN32
    else if (!strcmp(up_item, _("Copy to clipboard"))) {
	win32_process_graph(plot->spec, WIN32_TO_CLIPBOARD);
    } else if (!strcmp(up_item, _("Print"))) {
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

    pputs(prn, _("VAR roots (real, imaginary, modulus, frequency)"));
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
	gchar *title = g_strdup_printf("gretl: %s", _("VAR roots"));

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

    err = add_graph_to_session(spec->fname, fullname, type);

    if (!err) {
	remove_png_term_from_plotfile(fullname, spec);
	mark_plot_as_saved(spec);
    }
}

static void show_all_labels (png_plot *plot)
{
    FILE *fp;

    if (plot->spec->labeled != NULL) {
	free(plot->spec->labeled);
	plot->spec->labeled = NULL;
    }

    plot->spec->flags |= GPT_PRINT_MARKERS;

    gnuplot_png_init(plot->spec, &fp);

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
	gnuplot_png_init(plot->spec, &fp);
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

static gint plot_popup_activated (GtkWidget *w, gpointer data)
{
    gchar *item = (gchar *) data;
    gpointer ptr = g_object_get_data(G_OBJECT(w), "plot");
    png_plot *plot = (png_plot *) ptr;
    int killplot = 0;

    gtk_widget_destroy(plot->popup);
    plot->popup = NULL;

    if (!strcmp(item, _("Add another curve..."))) {
	dist_graph_add(plot);
    } else if (!strcmp(item, _("Save as PNG..."))) {
	plot->spec->termtype = GP_TERM_PNG;
        file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC, 
				  plot->spec, plot->shell);
    } else if (!strcmp(item, _("Save as PDF..."))) {
	plot->spec->termtype = GP_TERM_PDF;
        file_selector_with_parent(SAVE_GNUPLOT, FSEL_DATA_MISC, 
				  plot->spec, plot->shell);
    } else if (!strcmp(item, _("Save to session as icon"))) { 
	add_to_session_callback(plot->spec);
    } else if (plot_is_range_mean(plot) && !strcmp(item, _("Help"))) { 
	context_help(NULL, GINT_TO_POINTER(RMPLOT));
    } else if (plot_is_hurst(plot) && !strcmp(item, _("Help"))) { 
	context_help(NULL, GINT_TO_POINTER(HURST));
    } else if (!strcmp(item, _("Freeze data labels"))) {
	plot->spec->flags |= GPT_PRINT_MARKERS;
	redisplay_edited_plot(plot);
    } else if (!strcmp(item, _("Clear data labels"))) { 
	clear_labels(plot);
    } else if (!strcmp(item, _("All data labels"))) { 
	show_all_labels(plot);
    } else if (!strcmp(item, _("Zoom..."))) { 
	GdkCursor* cursor;

	cursor = gdk_cursor_new(GDK_CROSSHAIR);
	gdk_window_set_cursor(plot->canvas->window, cursor);
	gdk_cursor_unref(cursor);
	plot->status |= PLOT_ZOOMING;
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar), plot->cid, 
			   _(" Drag to define zoom rectangle"));
	create_selection_gc(plot);
    } else if (!strcmp(item, _("Restore full view"))) { 
	repaint_png(plot, PNG_UNZOOM);
    }
#ifdef GTK_PRINTING
    else if (!strcmp(item, _("Print..."))) { 
	gtk_print_graph(plot->spec->fname);
    }
#endif 
    else if (!strcmp(item, _("Display PDF"))) { 
	graph_display_pdf(plot->spec);
    } else if (!strcmp(item, _("OLS estimates"))) { 
	if (plot->spec != NULL) {
	    do_graph_model(plot->spec->reglist, plot->spec->fit);
	}
    } else if (!strcmp(item, _("Numerical values"))) {
	show_numbers_from_markers(plot->spec);
    } else if (!strcmp(item, _("Numerical summary"))) {
	boxplot_show_summary(plot->spec);
    } else if (!strcmp(item, _("Edit"))) { 
	start_editing_png_plot(plot);
    } else if (!strcmp(item, _("Close"))) { 
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
			 G_CALLBACK(color_popup_activated),
			 _(color_items[i]));
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
#if defined GTK_PRINTING
	N_("Print..."),
#elif defined(G_OS_WIN32)
	N_("Print"),
#endif
	N_("Display PDF"),
	N_("OLS estimates"),
	N_("Numerical values"),
	N_("Numerical summary"),
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
    static int pdf_ok = -1;
    int i;

    if (pdf_ok == -1) {
	pdf_ok = gnuplot_pdf_terminal();
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

int redisplay_edited_plot (png_plot *plot)
{
    gchar *plotcmd;
    FILE *fp;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "redisplay_edited_plot: plot = %p\n", (void *) plot);
#endif

    /* open file in which to dump plot specification */
    gnuplot_png_init(plot->spec, &fp);
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
    get_plot_ranges(plot);

    /* put the newly created PNG onto the plot canvas */
    return render_pngfile(plot, PNG_REDISPLAY);
}

/* preparation for redisplaying graph: here we handle the case where
   we're switching to a zoomed view (by use of a temporary gnuplot
   source file); then we get gnuplot to create a new PNG.
 */

static int repaint_png (png_plot *plot, int view)
{
    int err = 0;
    char zoomname[MAXLEN];
    gchar *plotcmd = NULL;

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
	gdk_window_set_cursor(plot->canvas->window, NULL);
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
	    plot->screen_xmin = event->x;
	    plot->screen_ymin = event->y;
	}
	return TRUE;
    }

    if (plot_doing_position(plot)) {
	if (plot->labelpos_entry != NULL) {
	    double dx, dy;
	    
	    if (get_data_xy(plot, event->x, event->y, &dx, &dy)) {
		gchar *posstr;

		posstr = g_strdup_printf("%g %g", dx, dy);
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

    if (!plot->err) {
	build_plot_menu(plot);
	gtk_menu_popup(GTK_MENU(plot->popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
    }

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
	add_to_session_callback(plot->spec);
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

static 
void plot_expose (GtkWidget *widget, GdkEventExpose *event,
		  GdkPixmap *dbuf_pixmap)
{
    /* Don't repaint entire window on each exposure */
    gdk_window_set_back_pixmap(widget->window, NULL, FALSE);

    /* Refresh double buffer, then copy the "dirtied" area to
       the on-screen GdkWindow */
    gdk_draw_drawable(widget->window,
		      widget->style->fg_gc[GTK_STATE_NORMAL],
		      dbuf_pixmap,
		      event->area.x, event->area.y,
		      event->area.x, event->area.y,
		      event->area.width, event->area.height);
}

static GdkPixbuf *gretl_pixbuf_new_from_file (const gchar *fname)
{
    GdkPixbuf *pbuf;
    GError *gerr = NULL;

    pbuf = gdk_pixbuf_new_from_file(fname, &gerr);

    if (pbuf == NULL) {
	verbose_gerror_report(gerr, "gdk_pixbuf_new_from_file");
	if (g_error_matches(gerr, G_FILE_ERROR, G_FILE_ERROR_INVAL)) {
	    gchar *trfname = NULL;
	    gsize bytes;

	    g_error_free(gerr);
	    gerr = NULL;

	    if (!g_utf8_validate(fname, -1, NULL)) {
		fprintf(stderr, "Trying g_locale_to_utf8 on filename\n");
		trfname = g_locale_to_utf8(fname, -1, NULL, &bytes, &gerr);
		if (trfname == NULL) {
		    verbose_gerror_report(gerr, "g_locale_to_utf8");
		}
	    } else {
		fprintf(stderr, "Trying g_locale_from_utf8 on filename\n");
		trfname = g_locale_from_utf8(fname, -1, NULL, &bytes, &gerr);
		if (trfname == NULL) {
		    verbose_gerror_report(gerr, "g_locale_from_utf8");
		}
	    }

	    if (trfname != NULL) {
		pbuf = gdk_pixbuf_new_from_file(trfname, &gerr);
		g_free(trfname);
		if (pbuf == NULL) {
		    verbose_gerror_report(gerr, "gdk_pixbuf_new_from_file");
		}
	    }
	}
	if (gerr != NULL) {
	    errbox(gerr->message);
	    g_error_free(gerr);
	}
    } 

    return pbuf;
}

/* The last step of redisplaying a graph after some change has been
   made: grab the gulot-generated PNG file, make a pixbuf out of it,
   and draw the pixbuf onto the canvas of the plot window.
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

    /* scrap any old record of which points are labeled */
    if (plot->spec->labeled != NULL) {
	free(plot->spec->labeled);
	plot->spec->labeled = NULL;
	if (!(plot->spec->flags & GPT_PRINT_MARKERS)) {
	    /* any markers will have disappeared on reprinting */
	    plot->format &= ~PLOT_MARKERS_UP;
	} 
    }

    gdk_draw_pixbuf(plot->pixmap, 
		    plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
		    pbuf, 0, 0, 0, 0, width, height,
		    GDK_RGB_DITHER_NONE, 0, 0);

    g_object_unref(pbuf);
    gretl_remove(pngname);
   
    if (view != PNG_START) { 
	/* we're changing the view, so refresh the whole canvas */
	gdk_draw_drawable(plot->canvas->window,
			  plot->canvas->style->fg_gc[GTK_STATE_NORMAL],
			  plot->pixmap,
			  0, 0,
			  0, 0,
			  plot->pixel_width, plot->pixel_height);

	if (view == PNG_ZOOM) {
	    plot->status |= PLOT_ZOOMED;
	} else if (view == PNG_UNZOOM) {
	    plot->status ^= PLOT_ZOOMED;
	}
    }

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

    if (plot->invert_gc != NULL) {
	g_object_unref(plot->invert_gc);
    }

    g_object_unref(plot->shell);

    free(plot);
}

static void set_approx_pixel_bounds (png_plot *plot, 
				     int max_num_width,
				     int max_num2_width)
{
    if (plot_has_xlabel(plot)) {
	plot->pixel_ymax = plot->pixel_height - 36;
    } else {
	plot->pixel_ymax = plot->pixel_height - 24;
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

    plot->pixel_xmax = plot->pixel_width - 20; 
    if (plot_has_y2axis(plot)) {
	plot->pixel_xmax -= 7 * (max_num2_width + 1);
    }
    if (plot_has_y2label(plot)) {
	plot->pixel_xmax -= 11;
    }

#if POINTS_DEBUG
    fprintf(stderr, "set_approx_pixel_bounds():\n"
	    " xmin=%d xmax=%d ymin=%d ymax=%d\n", 
	    plot->pixel_xmin, plot->pixel_xmax,
	    plot->pixel_ymin, plot->pixel_ymax);
    fprintf(stderr, "set_approx_pixel_bounds():\n"
	    " max_num_width=%d max_num2_width=%d\n", 
	    max_num_width, max_num2_width);
#endif
}

int ok_dumb_line (const char *s)
{
    if (strstr(s, "x2tics")) return 0;
    if (strstr(s, "set style line")) return 0;
    if (strstr(s, "set style inc")) return 0;
    return 1;
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

    fpin = gretl_fopen(plot->spec->fname, "r");
    if (fpin == NULL) {
	return 1;
    }

    build_path(dumbgp, gretl_dotdir(), "dumbplot.gp", NULL);
    build_path(dumbtxt, gretl_dotdir(), "gptdumb.txt", NULL);
    fpout = gretl_fopen(dumbgp, "w");
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
	} else if (ok_dumb_line(line)) {
	    fputs(line, fpout);
	}
	if (strstr(line, "x2range")) {
	    x2axis = 1;
	}
    }

    fclose(fpin);
    fclose(fpout);

    plotcmd = g_strdup_printf("\"%s\" \"%s\"", 
			      gretl_gnuplot_path(), 
			      dumbgp);
    err = gretl_spawn(plotcmd);
    g_free(plotcmd);

    gretl_remove(dumbgp);

    if (err) {
#if POINTS_DEBUG
	fputs("get_dumb_plot_yrange(): plot command failed\n", stderr);
#endif
	gretl_error_clear();
	return 1;
    } else {
	double y[16] = {0};
	int y_numwidth[16] = {0};
	int y2_numwidth[16] = {0};
	char numstr[32];
	int i, j, k, imin;

	fpin = gretl_fopen(dumbtxt, "r");
	if (fpin == NULL) {
	    return 1;
	}

	/* read the y-axis min and max from the ascii graph */

	gretl_push_c_numeric_locale();

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
#if POINTS_DEBUG
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

	gretl_pop_c_numeric_locale();

	fclose(fpin);
#if (POINTS_DEBUG == 0)
	gretl_remove(dumbtxt);
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

#if POINTS_DEBUG
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

/* Do a partial parse of the gnuplot source file: enough to determine
   the data ranges so we can read back the mouse pointer coordinates
   when the user moves the pointer over the graph.
*/

static int get_plot_ranges (png_plot *plot)
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

    if (got_x) {
	plot->status |= PLOT_HAS_XRANGE;
    } else {
	plot->status |= (PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
	return 1;
    }    

    /* get the "dumb" y coordinates only if we haven't got
       more accurate ones already */
    if (!plot_has_png_coords(plot)) { 
	err = get_dumb_plot_yrange(plot);
    }

    if (!err) {
	plot->status |= PLOT_HAS_YRANGE;
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
    plot->pixmap = NULL;
    plot->invert_gc = NULL;
    plot->spec = NULL;
    plot->editor = NULL;

    plot->pixel_width = 640;
    plot->pixel_height = 480;

    plot->xmin = plot->xmax = 0.0;
    plot->ymin = plot->ymax = 0.0;
    plot->xint = plot->yint = 0;

    plot->zoom_xmin = plot->zoom_xmax = 0.0;
    plot->zoom_ymin = plot->zoom_ymax = 0.0;
    plot->screen_xmin = plot->screen_ymin = 0;

    plot->pd = 0;
    plot->err = 0;
    plot->cid = 0;
    plot->status = 0;
    plot->format = 0;

    return plot;
}

static int gnuplot_show_png (const char *plotfile, const char *name,
			     GPT_SPEC *spec, int saved)
{
    GtkWidget *vbox;
    GtkWidget *canvas_hbox;
    GtkWidget *label_frame = NULL;
    GtkWidget *status_hbox = NULL;
    gchar *title = NULL;
    png_plot *plot;
    int polar = 0;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "gnuplot_show_png:\n plotfile='%s', spec=%p, saved=%d\n",
	    plotfile, (void *) spec, saved);
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
	strcpy(plot->spec->fname, plotfile);
    }

    if (saved) {
	plot->status |= PLOT_SAVED;
    }

    /* make png plot struct accessible via spec */
    plot->spec->ptr = plot;

    /* Parse the gnuplot source file.  If we hit errors here,
       flag this, but it's not necessarily a show-stopper in
       terms of simply displaying the graph. 
    */
    plot->err = read_plotspec_from_file(plot->spec, &plot->pd, &polar);

#if GPDEBUG
    fprintf(stderr, "gnuplot_show_png: read_plotspec_from_file returned %d\n",
	    plot->err);
#endif

    if (plot->err) {
	plot->status |= (PLOT_DONT_EDIT | PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
    } else if (cant_edit(plot->spec->code)) {
	if (plot->spec->n_markers > 0) {
	    plot->status |= (PLOT_DONT_EDIT | PLOT_DONT_ZOOM);
	    if (polar) {
		plot->format |= PLOT_POLAR;
	    }
	} else {
	    plot->status |= (PLOT_DONT_EDIT | PLOT_DONT_ZOOM | PLOT_DONT_MOUSE);
	}
    } else {
	set_plot_format_flags(plot);
    } 

    if (plot->spec->code == PLOT_VAR_ROOTS) {
	plot->pixel_width = plot->pixel_height;
    }

    if (plot->spec->flags & GPT_LETTERBOX) {
	plot->pixel_width = 680;
	plot->pixel_height = 400;
    }

    if (!plot->err) {
	get_plot_ranges(plot);
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
    g_free(title);

    gtk_window_set_resizable(GTK_WINDOW(plot->shell), FALSE);

    vbox = gtk_vbox_new(FALSE, 2);
    gtk_container_add(GTK_CONTAINER(plot->shell), vbox);

    g_signal_connect(G_OBJECT(plot->shell), "destroy",
		     G_CALLBACK(destroy_png_plot), plot);
    g_signal_connect(G_OBJECT(plot->shell), "key-press-event", 
		     G_CALLBACK(plot_key_handler), plot);

    /* box to hold canvas */
    canvas_hbox = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), canvas_hbox, TRUE, TRUE, 0);
    gtk_widget_show(canvas_hbox);

    /* eventbox and hbox for status area  */
    plot->statusarea = gtk_event_box_new();
    gtk_box_pack_start(GTK_BOX(vbox), plot->statusarea, FALSE, FALSE, 0);

    status_hbox = gtk_hbox_new (FALSE, 2);
    gtk_container_add(GTK_CONTAINER(plot->statusarea), status_hbox);
    gtk_widget_show (status_hbox);
    gtk_container_set_resize_mode (GTK_CONTAINER (status_hbox),
				   GTK_RESIZE_QUEUE);

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

    GTK_WIDGET_SET_FLAGS(plot->canvas, GTK_CAN_FOCUS);

    g_signal_connect(G_OBJECT(plot->canvas), "button-press-event", 
		     G_CALLBACK(plot_button_press), plot);
    g_signal_connect(G_OBJECT(plot->canvas), "button-release-event", 
		     G_CALLBACK(plot_button_release), plot);

    /* create the contents of the status area */
    if (plot_has_xrange(plot)) {
	/* cursor label (graph position indicator) */
	label_frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(label_frame), GTK_SHADOW_IN);

	plot->cursor_label = gtk_label_new(" ");
	gtk_container_add(GTK_CONTAINER(label_frame), plot->cursor_label);
	gtk_widget_show(plot->cursor_label);
    }

    /* the statusbar */
    plot->statusbar = gtk_statusbar_new();

    gtk_widget_set_size_request(plot->statusbar, 1, -1);
    gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(plot->statusbar), FALSE);

    gtk_container_set_resize_mode(GTK_CONTAINER (plot->statusbar),
				  GTK_RESIZE_QUEUE);
    plot->cid = gtk_statusbar_get_context_id(GTK_STATUSBAR(plot->statusbar),
					     "plot_message");

    if (!plot->err) {
	gtk_statusbar_push(GTK_STATUSBAR(plot->statusbar),
			   plot->cid, _(" Click on graph for pop-up menu"));
    }
    
    if (plot_has_xrange(plot)) {
	g_signal_connect(G_OBJECT(plot->canvas), "motion-notify-event",
			 G_CALLBACK(motion_notify_event), plot);
    }

    /* pack the widgets */
    gtk_box_pack_start(GTK_BOX(canvas_hbox), plot->canvas, FALSE, FALSE, 0);

    /* fill the status area */
    if (plot_has_xrange(plot)) {
	gtk_box_pack_start(GTK_BOX(status_hbox), label_frame, FALSE, FALSE, 0);
    }

    gtk_box_pack_start(GTK_BOX(status_hbox), plot->statusbar, TRUE, TRUE, 0); 

    /* show stuff */
    gtk_widget_show(plot->canvas);

    if (plot_has_xrange(plot)) {
	gtk_widget_show(label_frame);
    }

    gtk_widget_show(plot->statusbar);
    gtk_widget_show(plot->statusarea);

    gtk_widget_realize(plot->canvas);
    gdk_window_set_back_pixmap(plot->canvas->window, NULL, FALSE);

    if (plot_has_xrange(plot)) {
	gtk_widget_realize(plot->cursor_label);
	gtk_widget_set_size_request(plot->cursor_label, 160, -1);
    }

    gtk_widget_show(vbox);
    gtk_widget_show(plot->shell);       

    /* set the focus to the canvas area */
    gtk_widget_grab_focus(plot->canvas);  

    plot->pixmap = gdk_pixmap_new(plot->shell->window, 
				  plot->pixel_width, plot->pixel_height, 
				  -1);
    g_signal_connect(G_OBJECT(plot->canvas), "expose-event",
		     G_CALLBACK(plot_expose), plot->pixmap);

    err = render_pngfile(plot, PNG_START);
    if (err) {
	gtk_widget_destroy(plot->shell);
	plot = NULL;
    }

    return err;
}

/* @fname is the name of a pre-made PNG file */

int display_graph_file (const char *fname)
{
    return gnuplot_show_png(fname, NULL, NULL, 0);
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

    if (add_png_term_to_plotfile(fullname)) {
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
    }

    if (!err) {
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
    char bbname[MAXLEN];
    FILE *fp;
    char line[128];
    int plot_ret = -1, data_ret = -1;
    int ret = GRETL_PNG_OK;

    build_path(bbname, gretl_dotdir(), "gretltmp.png", ".bounds"); 
    fp = gretl_fopen(bbname, "r");

    if (fp == NULL) {
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


