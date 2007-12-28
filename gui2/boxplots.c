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

/* boxplots for gretl */

#include "gretl.h"
#include "session.h"
#include "boxplots.h"
#include "fileselect.h"

#ifdef G_OS_WIN32
# include <windows.h>
#endif

#include "gtkplot-lite.h"

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
    int width, height;
    double boxwidth;
    double gmax, gmin;
    GtkWidget *window, *area, *popup;
    GdkPixmap *pixmap;
    int saved;
} PLOTGROUP;

static char boxplottmp[MAXLEN];

static double headroom = 0.24;
static double scalepos = 75.0; /* was 60 */
static char boxfont[64] = "Helvetica";
static int boxfontsize = 12;

static int six_numbers (gpointer data);
static int dump_boxplot (PLOTGROUP *grp);

#ifdef G_OS_WIN32
static int cb_copy_image (gpointer data);
#else
int plot_to_xpm (const char *fname, gpointer data);
#endif


const char *get_boxdump_name (void)
{
    return boxplottmp;
}

/* Create a new backing pixmap of the appropriate size */

static gint
configure_event (GtkWidget *widget, GdkEventConfigure *event, gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    if (grp->pixmap) {
	g_object_unref(G_OBJECT(grp->pixmap));
    }

    grp->pixmap = gdk_pixmap_new(widget->window,
				 widget->allocation.width,
				 widget->allocation.height,
				 -1);

    gdk_draw_rectangle (grp->pixmap,
			widget->style->white_gc,
			TRUE,
			0, 0,
			widget->allocation.width,
			widget->allocation.height);

    return TRUE;
}

/* Redraw the screen from the backing pixmap */

static gint
expose_event (GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    gdk_draw_pixmap(widget->window,
		    widget->style->fg_gc[GTK_STATE_NORMAL],
		    grp->pixmap,
		    event->area.x, event->area.y,
		    event->area.x, event->area.y,
		    event->area.width, event->area.height);

    return FALSE;
}

static gboolean 
box_key_handler (GtkWidget *w, GdkEventKey *key, gpointer data)
{
    if (key->keyval == GDK_q) {
	gtk_widget_destroy(w);
    } else if (key->keyval == GDK_s) {
        file_selector(_("Save boxplot file"), SAVE_BOXPLOT_EPS, 
		      FSEL_DATA_MISC, data);
    } else if (key->keyval == GDK_p) {  
	six_numbers(data);
    }
#ifdef G_OS_WIN32
    else if (key->keyval == GDK_c) {  
	cb_copy_image(data);
    }
#else
    else if (key->keyval == GDK_c) { 
        file_selector(_("Save boxplot file"), SAVE_BOXPLOT_XPM, 
		      FSEL_DATA_MISC, data);	
    }
#endif
    return TRUE;
}

static gint box_popup_activated (GtkWidget *w, gpointer data)
{
    gchar *item = (gchar *) data;
    gpointer ptr = g_object_get_data(G_OBJECT(w), "group");
    PLOTGROUP *grp = (PLOTGROUP *) ptr;

    if (!strcmp(item, _("Numerical summary"))) {
        six_numbers(grp);
    } else if (!strcmp(item, _("Save to session as icon"))) {
        if (dump_boxplot(grp) == 0) {
	    add_boxplot_to_session(boxplottmp);
	    remove(boxplottmp);
	    grp->saved = 1;
	}
    } else if (!strcmp(item, _("Save as EPS..."))) {
        file_selector(_("Save boxplot file"), SAVE_BOXPLOT_EPS, 
		      FSEL_DATA_MISC, ptr);
    } else if (!strcmp(item, _("Save as PS..."))) {
        file_selector(_("Save boxplot file"), SAVE_BOXPLOT_PS, 
		      FSEL_DATA_MISC, ptr);
    } else if (!strcmp(item, _("Help"))) {
	context_help(NULL, GINT_TO_POINTER(BXPLOT));
    } else if (!strcmp(item, _("Close"))) { 
	gtk_widget_destroy(grp->popup);
	grp->popup = NULL;
        gtk_widget_destroy(grp->window);
	return TRUE;
    }

#ifdef G_OS_WIN32
    else if (!strcmp(item, _("Copy to clipboard"))) {
	cb_copy_image(ptr);
    }
#else
    else if (!strcmp(item, _("Save as XPM..."))) {
        file_selector(_("Save boxplot file"), SAVE_BOXPLOT_XPM, 
		      FSEL_DATA_MISC, ptr);
    }
#endif

    gtk_widget_destroy(grp->popup);
    grp->popup = NULL;

    return TRUE;
}

static GtkWidget *build_menu (PLOTGROUP *grp)
{
    GtkWidget *menu, *item;    
    static char *items[] = {
        N_("Numerical summary"),
	N_("Save to session as icon"),
        N_("Save as EPS..."),
        N_("Save as PS..."),
#ifdef G_OS_WIN32
	N_("Copy to clipboard"),
#else
	N_("Save as XPM..."),
#endif
	N_("Help"),
        N_("Close"),
	NULL
    };
    int i = 0;

    menu = gtk_menu_new();

    while (items[i]) {
	if (grp->saved && !strcmp(items[i], "Save to session as icon")) {
	    i++;
	    continue;
	}
	item = gtk_menu_item_new_with_label(_(items[i]));
        g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(box_popup_activated),
			 _(items[i]));
	g_object_set_data(G_OBJECT(item), "group", grp);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	gtk_widget_show(item);
	i++;
    }
    return menu;
}

static gint box_popup (GtkWidget *widget, GdkEventButton *event, 
		       gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    if (grp->popup) g_free(grp->popup);
    grp->popup = build_menu(grp);
    gtk_menu_popup(GTK_MENU(grp->popup), NULL, NULL, NULL, NULL,
		   event->button, event->time);
    return TRUE;
}

static void
setup_text (GtkWidget *area, GdkPixmap *pixmap,
	    GdkGC *gc, GtkPlotPC *pc, const char *text, 
	    double x, double y, unsigned just)
{
    if (pc != NULL) {
	GdkColor black, white;
	
	black.red = black.green = black.blue = 0;
	white.red = white.green = white.blue = 65535.0;

	gtk_plot_pc_draw_string (pc,
				 x, y + 4,        /* alignment bodge */
				 0,               /* angle */
				 &black, &white,  /* fore, back */
				 FALSE,           /* transparent? */
				 GTK_PLOT_BORDER_NONE, 0, 1, 0,
				 boxfont, boxfontsize,
				 just,
				 text);
    } else {
	PangoContext *context;
	PangoLayout *pl;
	int width, height;

	context = gtk_widget_get_pango_context(area);
	pl = pango_layout_new(context);
	pango_layout_set_text(pl, text, -1);
	pango_layout_get_pixel_size(pl, &width, &height);

	if (just == GTK_JUSTIFY_CENTER) {
	    x -= width / 2.0;
	}
	else if (just == GTK_JUSTIFY_RIGHT) {
	    y -= (double) height / 2.0; 
	    x -= width;
	}

	gdk_draw_layout(pixmap, gc, x, y, pl);
	g_object_unref(G_OBJECT(pl));
    }
}

static void 
draw_line (double *points, GtkWidget *area, GdkPixmap *pixmap,
	   GdkGC *gc, GtkPlotPC *pc)
{
    if (pc != NULL) {
	gtk_plot_pc_draw_line(pc, points[0], points[1], points[2], points[3]); 
    } else {
	gdk_draw_line(pixmap, gc, points[0], points[1], points[2], points[3]);
    }
}

static void
draw_outlier (double x, double y,
	      GtkWidget *area, GdkPixmap *pixmap, 
	      GdkGC *gc, GtkPlotPC *pc)
{
    int fullcircle = 360 * 64;

    if (pc != NULL) {
	gtk_plot_pc_draw_circle(pc, FALSE, x, y, 8);
    } else {
	gdk_draw_arc(pixmap, gc, FALSE, x - 4, y - 4, 8, 8, 0, fullcircle);
    }
}

static void 
draw_mean (double x, double y, 
	   GtkWidget *area, GdkPixmap *pixmap, 
	   GdkGC *gc, GtkPlotPC *pc)
{
    if (pc != NULL) {
	gtk_plot_pc_draw_line(pc, x - 4, y, x + 4, y);
	gtk_plot_pc_draw_line(pc, x, y + 4, x, y - 4);
    } else {
	gdk_draw_line(pixmap, gc, x - 4, y, x + 4, y);
	gdk_draw_line(pixmap, gc, x, y + 4, x, y - 4);
    }
}

static void 
place_plots (PLOTGROUP *plotgrp)
     /* calculate placement of plots based on the number of them
	and their respective max and min values */
{
    int i;
    double start = scalepos + plotgrp->width * (headroom / 2.0);
    double xrange = (1.0 - headroom) * plotgrp->width - scalepos;
    double boxwidth = xrange / (2.0 * plotgrp->nplots - 1.0);

    plotgrp->boxwidth = boxwidth;

    /* divide horizontal space between plots; also get global
       maximum and minimum y-values */

    for (i=0; i<plotgrp->nplots; i++) {
	plotgrp->plots[i].xbase = start + (2.0 * (i+1.0) - 2.0) * boxwidth;
	if (plotgrp->plots[i].max > plotgrp->gmax) 
	    plotgrp->gmax = plotgrp->plots[i].max;
	if (plotgrp->plots[i].min < plotgrp->gmin) 
	    plotgrp->gmin = plotgrp->plots[i].min;
    }
}

static void 
gtk_boxplot_yscale (PLOTGROUP *grp, GtkPlotPC *pc)
{
    double points[4];
    double top = (headroom / 2.0) * grp->height;
    double bottom = (1.0 - headroom / 2.0) * grp->height;
    gchar *tmp = NULL;
    GdkGC *gc = NULL;

    if (pc == NULL) {
	gc = grp->window->style->fg_gc[GTK_STATE_NORMAL];
    }

    /* draw vertical line */
    points[0] = points[2] = scalepos;
    points[1] = top;
    points[3] = bottom;
    draw_line(points, grp->area, grp->pixmap, gc, pc);

    /* draw backticks top and bottom */
    points[2] = points[0] - 5;
    points[1] = points[3] = top;
    draw_line (points, grp->area, grp->pixmap, gc, pc);
    points[1] = points[3] = bottom;
    draw_line(points, grp->area, grp->pixmap, gc, pc);

    /* draw backtick at middle */
    points[1] = points[3] = top + (bottom - top) / 2.0;
    draw_line(points, grp->area, grp->pixmap, gc, pc);
    
    /* mark three values on scale */
    tmp = g_strdup_printf("%.4g", grp->gmax);
    gretl_fix_exponent(tmp);
    setup_text(grp->area, grp->pixmap, gc, pc, tmp, scalepos - 8, top, 
	       GTK_JUSTIFY_RIGHT);
    g_free(tmp);

    tmp = g_strdup_printf("%.4g", grp->gmin);
    gretl_fix_exponent(tmp);
    setup_text(grp->area, grp->pixmap, gc, pc, tmp, scalepos - 8, bottom, 
	       GTK_JUSTIFY_RIGHT);
    g_free(tmp);

    tmp = g_strdup_printf("%.4g", (grp->gmax + grp->gmin) / 2.0);
    gretl_fix_exponent(tmp);
    setup_text(grp->area, grp->pixmap, gc, pc, tmp, scalepos - 8, 
	       top + (bottom - top) / 2.0, GTK_JUSTIFY_RIGHT);
    g_free(tmp);
}

static void 
gtk_area_boxplot (BOXPLOT *plot, GtkWidget *area, GdkPixmap *pixmap,
		  GtkStyle *style, GtkPlotPC *pc,
		  int height, double boxwidth, 
		  double gmax, double gmin, char *numbers)
{
    double points[4];
    double ybase = height * headroom / 2.0;
    double xcenter = plot->xbase + boxwidth / 2.0;
    double scale = (1.0 - headroom) * height / (gmax - gmin);
    double mean, median, uq, lq, maxval, minval;
    double conflo = 0., confhi = 0.;
    double nameoff = headroom / 4.0;
    GdkGC *gc = NULL;

    if (pc == NULL) {
	gc = style->fg_gc[GTK_STATE_NORMAL];
    }

    if (na(plot->mean)) {
	mean = NADBL;
    } else {
	mean = ybase + (gmax - plot->mean) * scale;
    }

    median = ybase + (gmax - plot->median) * scale;
    uq = ybase + (gmax - plot->uq) * scale;
    lq = ybase + (gmax - plot->lq) * scale;

    if (plot->outliers == NULL) {
	maxval = ybase + (gmax - plot->max) * scale;
	minval = ybase + (gmax - plot->min) * scale;
    } else {
	maxval = ybase + (gmax - plot->outliers->rmax) * scale;
	minval = ybase + (gmax - plot->outliers->rmin) * scale;
    }
	
    if (!na(plot->conf[0])) { /* confidence intervals defined */
	if (plot->conf[1] > plot->uq) {
	    confhi = uq;
	} else {
	    confhi = ybase + (gmax - plot->conf[1]) * scale;
	}
	if (plot->conf[0] < plot->lq) { 
	    conflo = lq;
	} else {
	    conflo = ybase + (gmax - plot->conf[0]) * scale;
	}
    }

    if (confhi == 0.) {
	/* no notches: draw simple inter-quartile box */
	if (pc != NULL) { 
	    gtk_plot_pc_draw_rectangle (pc,
					FALSE,
					plot->xbase, uq,
					boxwidth,
					lq - uq);
	} else {
	    gdk_draw_rectangle (pixmap, 
				gc, 
				FALSE, /* filled ? */
				plot->xbase, 
				uq, 
				boxwidth, 
				lq - uq);
	}
    } else { 
	/* draw notched boxes */
	if (pc != NULL) {
	    GtkPlotPoint points[10];

	    points[0].x = points[6].x = points[7].x = points[9].x = 
		plot->xbase;
	    points[0].y = points[1].y = uq;
	    points[1].x = points[2].x = points[4].x = points[5].x =
		plot->xbase + boxwidth;
	    points[5].y = points[6].y = lq;
	    points[3].y = points[8].y = median;
	    points[2].y = points[9].y = confhi;
	    points[4].y = points[7].y = conflo;
	    points[8].x = plot->xbase + .10 * boxwidth;
	    points[3].x = plot->xbase + .90 * boxwidth;

	    gtk_plot_pc_draw_polygon (pc,
				      FALSE, /* filled ? */
				      points,
				      10);
	} else {
	    GdkPoint points[10];

	    points[0].x = points[6].x = points[7].x = points[9].x = 
		plot->xbase;
	    points[0].y = points[1].y = uq;
	    points[1].x = points[2].x = points[4].x = points[5].x =
		plot->xbase + boxwidth;
	    points[5].y = points[6].y = lq;
	    points[3].y = points[8].y = median;
	    points[2].y = points[9].y = confhi;
	    points[4].y = points[7].y = conflo;
	    points[8].x = plot->xbase + .10 * boxwidth;
	    points[3].x = plot->xbase + .90 * boxwidth;

	    gdk_draw_polygon (pixmap,
			      gc,
			      FALSE, /* filled ? */
			      points,
			      10);
	}
    }

    /* draw line at median */
    points[0] = plot->xbase + ((confhi > 0.)? (0.1 * boxwidth) : 0.0);
    points[1] = points[3] = median;
    points[2] = plot->xbase + ((confhi > 0.)? (0.9 * boxwidth) : boxwidth);
    draw_line(points, area, pixmap, gc, pc); /* was whitegc */

    /* draw '+' at mean */
    if (!na(mean)) {
	draw_mean(xcenter, mean, area, pixmap, gc, pc);
    }

    /* insert numerical values for median and quartiles? */
    if (numbers != NULL) {
	char numstr[12];
	double x = xcenter + .55 * boxwidth;

	sprintf(numstr, numbers, plot->uq);
	setup_text (area, pixmap, gc, pc, numstr, 
		    x, uq, GTK_JUSTIFY_LEFT);
	sprintf(numstr, numbers, plot->median);
	setup_text (area, pixmap, gc, pc, numstr, 
		    x, median, GTK_JUSTIFY_LEFT);
	sprintf(numstr, numbers, plot->lq);
	setup_text (area, pixmap, gc, pc, numstr, 
		    x, lq, GTK_JUSTIFY_LEFT);	
    }

    /* draw line to maximum value */
    points[0] = points[2] = xcenter;
    points[1] = maxval;
    points[3] = uq;
    draw_line(points, area, pixmap, gc, pc);

    /* plus a little crossbar */
    points[0] = xcenter - 5.0;
    points[2] = xcenter + 5.0;
    points[1] = points[3] = maxval;
    draw_line(points, area, pixmap, gc, pc);

    /* draw line to minimum value */
    points[0] = points[2] = xcenter;
    points[1] = lq;
    points[3] = minval;
    draw_line(points, area, pixmap, gc, pc);

    /* plus a little crossbar */
    points[0] = xcenter - 5.0;
    points[2] = xcenter + 5.0;
    points[1] = points[3] = minval;
    draw_line(points, area, pixmap, gc, pc);

    /* draw outliers, if any */
    if (plot->outliers != NULL) {
	double y, ybak = NADBL;
	int i;

	for (i=0; i<plot->outliers->n; i++) {
	    /* fprintf(stderr, "outlier: %g\n", plot->outliers->vals[i]); */
	    y = ybase + (gmax - plot->outliers->vals[i]) * scale;
	    if (y == ybak) y += 1.0;
	    draw_outlier(xcenter, y, area, pixmap, gc, pc);
	    ybak = y;
	}
    }

    /* write name of variable beneath */
    if (plot->bool) {
	nameoff = headroom / 3.5;
    }
    setup_text(area, pixmap, gc, pc, plot->varname, xcenter, 
	       height * (1.0 - nameoff), GTK_JUSTIFY_CENTER);
    if (plot->bool) {
	setup_text(area, pixmap, gc, pc, plot->bool, xcenter, 
		   height * (1.0 - headroom/6.0), GTK_JUSTIFY_CENTER);
    }
}

static void 
destroy_boxplots (GtkWidget *w, gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    int i;

    for (i=0; i<grp->nplots; i++) { 
	free(grp->plots[i].outliers);
	free(grp->plots[i].bool);
    }

    free(grp->plots);
    free(grp->numbers);
    g_object_unref(G_OBJECT(grp->pixmap));
    free(grp);

    remove(boxplottmp);
}

static GtkWidget *
make_area (PLOTGROUP *grp)
{
    grp->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(grp->window), _("gretl: boxplots"));
    gtk_window_set_resizable(GTK_WINDOW(grp->window), FALSE);

    /* Create the drawing area */
    grp->area = gtk_drawing_area_new ();

    gtk_widget_set_events (grp->area, GDK_EXPOSURE_MASK
			   | GDK_LEAVE_NOTIFY_MASK
			   | GDK_BUTTON_PRESS_MASK);
    
    gtk_widget_set_sensitive(grp->area, TRUE);

    g_signal_connect(G_OBJECT(grp->area), "configure_event",
		     G_CALLBACK(configure_event), grp);

    g_signal_connect(G_OBJECT(grp->area), "expose_event",
		     G_CALLBACK(expose_event), grp);

    g_signal_connect(G_OBJECT(grp->area), "button_press_event", 
		     G_CALLBACK(box_popup), grp);

    g_signal_connect(G_OBJECT(grp->window), "key_press_event", 
		     G_CALLBACK(box_key_handler), grp);

    g_signal_connect(G_OBJECT(grp->window), "destroy",
		     G_CALLBACK(destroy_boxplots), grp);

    gtk_widget_set_size_request (GTK_WIDGET(grp->area),
				 grp->width, grp->height); 

    gtk_widget_show (grp->area);

    gtk_container_add(GTK_CONTAINER(grp->window), grp->area);

    return grp->window;
}

static double median (double *x, const int n)
{
    int n2;
    double xx;

    qsort(x, n, sizeof *x, gretl_compare_doubles);

    n2 = n/2;
    xx = (n % 2)? x[n2] : 0.5 * (x[n2 - 1] + x[n2]);
    return xx;
}

static double 
quartiles (const double *x, const int n, BOXPLOT *box)
{
    int n2;
    double xx;

    n2 = n/2;
    xx = (n % 2)? x[n2] : 0.5 * (x[n2 - 1] + x[n2]);

    if (box != NULL) {
	box->median = xx;
	if (n % 2) {
	    box->lq = quartiles(x, n2 + 1, NULL);
	    box->uq = quartiles(x + n2, n2 + 1, NULL);
	} else {
	    box->lq = quartiles(x, n2, NULL);
	    box->uq = quartiles(x + n2, n2, NULL);
	}
    }
    return xx;
}

static int 
add_outliers (const double *x, const int n, BOXPLOT *box)
{
    int i, nout = 0;
    double iqr = 1.5 * (box->uq - box->lq);
    
    box->outliers = NULL;

    for (i=0; i<n; i++) 
	if (x[i] < (box->lq - iqr) || x[i] > (box->uq + iqr)) 
	    nout++;
 
    if (nout > 0) {
	int nlo = 0, nhi = 0;

	box->outliers = malloc(sizeof *box->outliers);
	if (box->outliers == NULL) return 1;
	box->outliers->n = nout;
	box->outliers->vals = malloc (nout * sizeof(double));
	if (box->outliers->vals == NULL) {
	    free(box->outliers);
	    box->outliers = NULL;
	    return 1;
	}
	nout = 0;
	for (i=0; i<n; i++) {
	    if (x[i] < (box->lq - iqr)) {
		box->outliers->vals[nout++] = x[i];
		nlo++;
	    }
	    else if (x[i] > (box->uq + iqr)) {
		box->outliers->vals[nout++] = x[i];
		nhi++;
	    }
	}
	box->outliers->rmin = (nlo)? x[nlo] : box->min;
	box->outliers->rmax = (nhi)? x[n-nhi-1] : box->max;
    }
    return 0;
}

#define ITERS 560
#define CONFIDENCE 90

static int 
median_interval (double *x, int n, double *low, double *high)
     /* obtain bootstrap estimate of 90% confidence interval
	for the sample median of data series x; return low and
	high values in 'low' and 'high' */
{
    double *medians, *samp;
    int i, j, t;

    medians = malloc (ITERS * sizeof *medians);
    if (medians == NULL) return 1;

    samp = malloc (n * sizeof *samp);
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
	medians[i] = median(samp, n);
    }

    /* sort the sample medians */
    qsort(medians, ITERS, sizeof *medians, gretl_compare_doubles);
    
    /* return the right values */
    j = 100 / ((100 - CONFIDENCE) / 2);
    *low = medians[ITERS / j];
    *high = medians[ITERS - ITERS/j];

    free(samp);
    free(medians);

    return 0;
}

int ps_print_plots (const char *fname, int flag, gpointer data) 
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GtkPlotPS *ps;
    int i, eps = 1, orient = GTK_PLOT_PORTRAIT;
    gdouble pscale = 0.8;

    if (flag == SAVE_BOXPLOT_PS) {
	pscale = 1.2;
	eps = 0;
	orient = GTK_PLOT_LANDSCAPE;
    }

    ps = GTK_PLOT_PS(gtk_plot_ps_new (fname, 
				      orient, 
				      eps, 
				      GTK_PLOT_LETTER, 
				      pscale, pscale));

    if (ps == NULL) return 1;

    if (eps) {
	ps->page_width = pscale * grp->width;
	ps->page_height = pscale * grp->height;
    }

    if (!psinit(GTK_PLOT_PC(ps))) {
	return 1;
    }

    gtk_psfont_init();

    for (i=0; i<grp->nplots; i++)
	gtk_area_boxplot (&grp->plots[i], 
			  NULL, NULL, NULL, &ps->pc, 
			  grp->height, grp->boxwidth, 
			  grp->gmax, grp->gmin, grp->numbers);
    
    gtk_boxplot_yscale(grp, &ps->pc);

    psleave(GTK_PLOT_PC(ps));

    gtk_object_destroy(GTK_OBJECT(ps));
    gtk_psfont_unref();

    return 0;
}

static void real_six_numbers (BOXPLOT *plt, int offset, int do_mean,
			      PRN *prn)
{
    if (plt->bool != NULL) {
	pprintf(prn, "%s\n %-*s", plt->varname, offset - 1, plt->bool);
    } else {
	pprintf(prn, "%-*s", offset, plt->varname);
    }

    if (do_mean) {
	pprintf(prn, "%9.5g", plt->mean);
    }

    pprintf(prn, "%9.5g%9.5g%9.5g%9.5g%9.5g", plt->min, 
	    plt->lq, plt->median, plt->uq, plt->max);

    if (plt->n > 0) {
	pprintf(prn, "  (n=%d)\n", plt->n);
    } else {
	pputc(prn, '\n');
    }
}

static void five_numbers_with_interval (BOXPLOT *plt, int offset, PRN *prn)
{
    char tmp[32];

    sprintf(tmp, "%.5g - %.5g", plt->conf[0], plt->conf[1]);

    if (plt->bool != NULL) {
	pprintf(prn, "%s\n %-*s", plt->varname, offset - 1, plt->bool);
    } else {
	pprintf(prn, "%-*s", offset, plt->varname);
    }
	    
    pprintf(prn, "%8.5g%10.5g%10.5g%17s%10.5g%10.5g\n",
	    plt->min, plt->lq, plt->median,
	    tmp, plt->uq, plt->max);
}

static int has_mean (PLOTGROUP *grp)
{
    int i;

    for (i=0; i<grp->nplots; i++) {
	if (na(grp->plots[i].mean)) {
	    return 0;
	}
    }

    return 1;
}

static int get_format_offset (PLOTGROUP *grp)
{
    int L = 6;
    int i, n;

    for (i=0; i<grp->nplots; i++) {
	if (grp->plots[i].bool != NULL) {
	    n = strlen(grp->plots[i].bool);
	} else {
	    n = strlen(grp->plots[i].varname);
	}
	if (n > L) {
	    L = n;
	}	    
    }

    return L + 2;
}

static int six_numbers (gpointer data) 
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    PRN *prn;
    int offset;
    int i;

    if (bufopen(&prn)) {
	return 1;
    }

    offset = get_format_offset(grp);

    if (na(grp->plots[0].conf[0])) { 
	/* no confidence intervals */
	int do_mean = has_mean(grp);

	pprintf(prn, "%s\n\n", _("Numerical summary"));

	if (do_mean) {
	    pprintf(prn, "%*s%9s%9s%9s%9s%9s\n", offset + 9, _("mean"),
		    "min", "Q1", _("median"), "Q3", "max");
	} else {
	    pprintf(prn, "%*s%10s%10s%10s%10s\n", offset + 10,
		    "min", "Q1", _("median"), "Q3", "max");
	}	    

	for (i=0; i<grp->nplots; i++) {
	    real_six_numbers(&grp->plots[i], offset, do_mean, prn);
	}
    } else { 
	/* with confidence intervals */
	pprintf(prn, "%s\n\n", _("Numerical summary with bootstrapped confidence "
				 "interval for median"));	 

	pprintf(prn, "%*s%10s%10s%17s%10s%10s\n",
		offset + 8, "min", "Q1", _("median"), 
		/* xgettext:no-c-format */
		_("(90% interval)"), 
		"Q3", "max");

	for (i=0; i<grp->nplots; i++) {
	    five_numbers_with_interval(&grp->plots[i], offset, prn);
	}
    }

    view_buffer(prn, 78, 240, _("gretl: boxplot data"), BXPLOT, NULL);

    return 0;
}

static void read_boxrc (PLOTGROUP *grp);

int boxplots (int *list, char **bools, double ***pZ, const DATAINFO *pdinfo, 
	      gretlopt opt)
{
    int i, j, n = pdinfo->t2 - pdinfo->t1 + 1;
    double *x;
    PLOTGROUP *plotgrp;
    int width = 576, height = 448;

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
	plotgrp->plots[i].mean = gretl_mean(0, n, x);
	qsort(x, n, sizeof *x, gretl_compare_doubles);
	plotgrp->plots[i].min = x[0];
	plotgrp->plots[i].max = x[n-1];
	quartiles(x, n, &plotgrp->plots[i]);
	plotgrp->plots[i].n = n;

	if (opt & OPT_O) {
	    /* notched boxplots wanted */
	    if (median_interval(x, n, &plotgrp->plots[i].conf[0],
				&plotgrp->plots[i].conf[1])) {
		errbox (_("Couldn't obtain confidence interval"));
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

    plotgrp->height = height;
    plotgrp->width = width;
    plotgrp->numbers = NULL;

    plotgrp->show_outliers = 0;
    read_boxrc(plotgrp);

    if (plotgrp->show_outliers) {
	for (i=0; i<plotgrp->nplots; i++) {
	    n = ztox(list[i+1], x, (const double **) *pZ, pdinfo);
	    qsort(x, n, sizeof *x, gretl_compare_doubles);
	    add_outliers(x, n, &plotgrp->plots[i]);
	}
    }

    free(x);

    plotgrp->popup = NULL;
    plotgrp->pixmap = NULL;
    place_plots(plotgrp);

    if (make_area(plotgrp) == NULL) {
	free(plotgrp->plots);
	free(plotgrp);
	return 1;
    }

    gtk_widget_show(plotgrp->window);

    for (i=0; i<plotgrp->nplots; i++)
	gtk_area_boxplot (&plotgrp->plots[i], 
			  plotgrp->area, plotgrp->pixmap,
			  plotgrp->window->style, NULL, 
			  plotgrp->height, plotgrp->boxwidth, 
			  plotgrp->gmax, plotgrp->gmin,
			  plotgrp->numbers);
    
    gtk_boxplot_yscale(plotgrp, NULL);

    if (opt & OPT_B) {
	dump_boxplot(plotgrp);
    }
    
    return 0;
}

/* copy functions */

#ifdef G_OS_WIN32

static void image_copy_err (const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    errbox(_("Failed to set clipboard data"));
}

static int cb_copy_image (gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GdkImage *image;
    int i;
    guint32 pixel, white_pixel;
    size_t image_bytes, dibsize, linelen;
    size_t palsize = sizeof(RGBQUAD) * 2;
    HANDLE hDIB;
    BOOL ret;
    BITMAPINFOHEADER *hdr;

    image = gdk_image_get(grp->pixmap, 0, 0, 
			  grp->width, grp->height);

    white_pixel = pow(2, image->depth) - 1;
    linelen = ((grp->width/8 - 1)/4 + 1) * 4;
    image_bytes = grp->height * linelen;

    /* allocate room for DIB */
    dibsize = sizeof(BITMAPINFOHEADER) + palsize + image_bytes;

    hDIB = GlobalAlloc(GMEM_MOVEABLE | GMEM_DDESHARE, dibsize);
    if (hDIB == NULL) {
	image_copy_err("Failed to allocate DIB");
	return FALSE;
    }

    /* fill header info */
    ret = FALSE;
    hdr = GlobalLock(hDIB);
    if (hdr) {
	hdr->biSize = sizeof(BITMAPINFOHEADER);
	hdr->biWidth = grp->width;
	hdr->biHeight = grp->height; 
	hdr->biPlanes = 1;
	hdr->biBitCount = 1;
	hdr->biCompression = BI_RGB; /* none */
	hdr->biSizeImage = image_bytes; /* not needed? */
	hdr->biXPelsPerMeter = 0;
	hdr->biYPelsPerMeter = 0;
	hdr->biClrUsed = 0; 
	hdr->biClrImportant = 0; /* all */
          
	GlobalUnlock(hDIB);
	ret = TRUE;
    } else {
	image_copy_err("Failed to lock DIB Header");
    }

    /* fill color map */
    if (ret) {
	char *bmp = GlobalLock(hDIB);
      
	ret = FALSE;
	if (bmp) {
	    RGBQUAD *pal = (RGBQUAD *)(bmp + sizeof(BITMAPINFOHEADER));

	    /* white, for cleared bits */
	    pal[0].rgbBlue = pal[0].rgbGreen = pal[0].rgbRed = 255;
	    pal[0].rgbReserved = 0;
	    /* black, for set bits */
	    pal[1].rgbBlue = pal[1].rgbGreen = pal[1].rgbRed = 0;
	    pal[1].rgbReserved = 0;
	  
	    ret = TRUE;
	    GlobalUnlock(hDIB);
	} else {
	    image_copy_err("Failed to lock DIB Palette");
	}
    } 
  
    /* copy data to DIB */
    if (ret) {
	unsigned char *data = GlobalLock(hDIB);
      
	ret = FALSE;
	if (data) {
	    unsigned char c;
	    int x, y;
	    int pad = linelen - grp->width/8 - (grp->width % 8 > 0);
	    
	    data += (sizeof(BITMAPINFOHEADER) + palsize);

	    for (y=grp->height-1; y>=0; y--) {
		i = 0; c = 0;
		for (x=0; x<grp->width; x++) {
		    pixel = gdk_image_get_pixel(image, x, y);
		    if (pixel != white_pixel) c |= (1 << (7-i)); 
		    i++;
		    if (i == 8) { /* done 8 bits -> ship out char */
			*data++ = c;
			i = 0;
			c = 0;
		    }
		}
		for (i=0; i<pad; i++) *data++ = 0;
	    }
	    ret = TRUE;
	    GlobalUnlock (hDIB);
	} else {
	    image_copy_err("Failed to lock DIB Data");
	}
    } /* copy data to DIB */
  
    /* copy DIB to ClipBoard */
    if (ret) {      
	if (!OpenClipboard (NULL)) {
	    image_copy_err("Cannot open the Clipboard!");
	    ret = FALSE;
	} else {
	    if (ret && !EmptyClipboard ()) {
		image_copy_err("Cannot empty the Clipboard");
		ret = FALSE;
	    }
	    if (ret) {
		if (NULL != SetClipboardData(CF_DIB, hDIB)) {
		    hDIB = NULL; /* data now owned by clipboard */
		} else {
		    errbox(_("Failed to set clipboard data"));
		}
	    }
	    if (!CloseClipboard()) {
		errbox(_("Failed to close Clipboard"));
	    }
	}
    }

    if (hDIB) GlobalFree(hDIB);
  
    g_object_unref(G_OBJECT(image));
  
    return ret;
} 

#else /* end G_OS_WIN32 */

int plot_to_xpm (const char *fname, gpointer data)
{    
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GdkImage *image;
    int i, j;
    guint32 pixel, white_pixel;
    FILE *fp;
    
    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	file_write_errbox(fname);
	return 1;
    }

    fprintf(fp, "/* XPM */\n"
	    "static char *boxplot[] = {\n"
	    "/* width height ncolors chars_per_pixel */\n"
	    "\"%d %d 2 1\",\n"
	    "/* colors */\n"
	    "\"  c white\",\n"
	    "\". c black\",\n"
	    "/* pixels */\n", grp->width, grp->height);

    image = gdk_image_get (grp->pixmap, 0, 0, 
			   grp->width, grp->height);

    white_pixel = pow(2, image->depth) - 1;

    for (i=0; i<grp->height; i++) {
	fprintf(fp, "\"");
	for (j=0; j<grp->width; j++) {
	    pixel = gdk_image_get_pixel(image, j, i);
	    if (pixel == white_pixel) fprintf(fp, " ");
	    else fprintf(fp, ".");
	}
	fprintf(fp, "\"%s\n", (i<grp->height-1)? "," : "};");
    }

    g_object_unref(G_OBJECT(image));

    fclose(fp);

    return 0;
}

#endif 

static FILE *get_custom_file (void)
{
    char boxrc[FILENAME_MAX];
    FILE *fp;
    
    fp = fopen(".boxplotrc", "r");

    if (fp == NULL) {
	sprintf(boxrc, "%s.boxplotrc", paths.workdir);
	fp = gretl_fopen(boxrc, "r");
    }

    if (fp == NULL) {
	fp = fopen("plotconfig.txt", "r");
    }

    if (fp == NULL) {
	sprintf(boxrc, "%splotconfig.txt", paths.workdir);
	fp = gretl_fopen(boxrc, "r");
    }    

    return fp;
}

static void read_boxrc (PLOTGROUP *grp)
{
    FILE *fp;

    grp->gmax = grp->gmin = NADBL;

    fp = get_custom_file();

    if (fp != NULL) {
	char line[80], key[18], val[32];

	while (fgets(line, 79, fp)) {
	    if (line[0] == '#') continue;
	    if (sscanf(line, "%17s = %31s", key, val) == 2) {
		key[17] = '\0';
		val[31] = '\0';
		if (!strcmp(key, "max")) { 
		    grp->gmax = atof(val);
		}
		else if (!strcmp(key, "min")) { 
		    grp->gmin = atof(val);
		}
		else if (!strcmp(key, "font")) { 
		    strncpy(boxfont, val, 63);
		    boxfont[63] = '\0';
		}
		else if (!strcmp(key, "fontsize")) {
		    boxfontsize = atoi(val);
		}
		else if (!strcmp(key, "width") && atoi(val) > 0) {
		    grp->width = atoi(val);
		}
		else if (!strcmp(key, "height") && atoi(val) > 0) { 
		    grp->height = atoi(val);
		}
		else if (!strcmp(key, "numbers") && 
			 (grp->numbers = malloc(8))) {
		    grp->numbers[0] = 0;
		    strncat(grp->numbers, val, 7);
		}
		else if (!strcmp(key, "outliers") && !strcmp(val, "true")) {
		    grp->show_outliers = 1;
		}
	    }
	}
	fclose (fp);
    }

    if (na(grp->gmax)) {
	grp->gmax = grp->plots[0].max;
    }

    if (na(grp->gmin)) {
	grp->gmin = grp->plots[0].min;
    }
}

static void 
prep_confint_for_printing (BOXPLOT *plt, double *c0, double *c1)
{
    if (na(plt->conf[0]) || na(plt->conf[1]) || plt->conf[0] == plt->conf[1]) {
	*c0 = 0.0;
	*c1 = 0.0;
    } else {
	*c0 = plt->conf[0];
	*c1 = plt->conf[1];
    }
}

static int dump_boxplot (PLOTGROUP *grp)
{
    FILE *fp;
    int i;
    BOXPLOT *plt;

    build_path(boxplottmp, paths.dotdir, "boxdump.tmp", NULL);

    fp = gretl_fopen(boxplottmp, "w");
    if (fp == NULL) {
	file_write_errbox(boxplottmp);
	return 1;
    }

    fprintf(fp, "# boxplot generated by gretl\n");
    fprintf(fp, "nplots = %d\n", grp->nplots);
    fprintf(fp, "numbers = %s\n", (grp->numbers == NULL)? 
	    "NULL" : grp->numbers);
    fprintf(fp, "width height = %d %d\n", grp->width, grp->height);
    fprintf(fp, "boxwidth = %g\n", grp->boxwidth);
    fprintf(fp, "gmax gmin = %g %g\n", grp->gmax, grp->gmin);

    for (i=0; i<grp->nplots; i++) {
	double c0, c1;

	plt = &grp->plots[i];
	prep_confint_for_printing(plt, &c0, &c1);
	fprintf(fp, "%d median = %g\n", i, plt->median);
	fprintf(fp, "%d conf = %g %g\n", i, c0, c1);
	fprintf(fp, "%d quartiles = %g %g\n", i, plt->uq, plt->lq);
	fprintf(fp, "%d maxmin = %g %g\n", i, plt->max, plt->min);
	fprintf(fp, "%d xbase = %g\n", i, plt->xbase);
	fprintf(fp, "%d varname = %s", i, plt->varname);
	if (plt->bool != NULL) {
	    fprintf(fp, " %s\n", plt->bool);
	} else {
	    fprintf(fp, "\n");
	}
	if (plt->outliers != NULL) {
	    int j;

	    fprintf(fp, "%d n_outliers = %d\n", i, plt->outliers->n);
	    fprintf(fp, "%d rmax rmin = %g %g\n", i, 
		    plt->outliers->rmax, plt->outliers->rmin);
	    for (j=0; j<plt->outliers->n; j++) 
		fprintf(fp, "%d vals %g\n", i, plt->outliers->vals[j]); 
	} else {
	    fprintf(fp, "%d n_outliers = 0\n", i);
	}
	fprintf(fp, "%d mean = %g\n", i, plt->mean);
	fprintf(fp, "%d nobs = %d\n", i, plt->n);
    }

    fclose(fp);

    return 0;
}

static void get_bool_from_line (const char *line, BOXPLOT *plt)
{
    char boolstr[128];

    if (sscanf(line, "%*d varname = %*s %127s", boolstr) == 1) {
	plt->bool = g_strdup(boolstr);
    }
}

static void screen_confidence_interval (BOXPLOT *plt) 
{
    if (plt->conf[0] == plt->conf[1]) {
	plt->conf[0] = NADBL;
	plt->conf[1] = NADBL;
    }
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
	free(plt->outliers);
	plt->outliers = NULL;
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

static void maybe_get_plot_mean (BOXPLOT *plt, FILE *fp)
{
    char line[80];
    long pos = ftell(fp);
    double x = NADBL;
    int n = 0;

    plt->mean = NADBL;
    plt->n = 0;

    if (fgets(line, sizeof line, fp) == NULL) {
	return;
    } else if (sscanf(line, "%*d mean = %lf", &x) == 1) {
	plt->mean = x;
    } else {
	fseek(fp, pos, SEEK_SET);
	return;
    }

    pos = ftell(fp);

    if (fgets(line, sizeof line, fp) == NULL) {
	return;
    } else if (sscanf(line, "%*d nobs = %d", &n) == 1) {
	plt->mean = x;
    } else {
	fseek(fp, pos, SEEK_SET);
    }
}

int retrieve_boxplot (const char *fname)
{
    FILE *fp;
    int i, j, nout;
    PLOTGROUP *grp = NULL;
    BOXPLOT *plt = NULL;
    char line[80], numstr[24];
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	file_read_errbox(fname);
	return 1;
    }

    grp = mymalloc(sizeof *grp);
    if (grp == NULL) {
	fclose(fp);
	return 1;
    }

    grp->saved = 1;
    grp->numbers = NULL;

    for (i=0; i<6 && fgets(line, 79, fp); i++) {
	if (i == 1 && sscanf(line, "nplots = %d", &grp->nplots) != 1) {
	    goto corrupt;
	} else if (i == 2 && sscanf(line, "numbers = %7s", numstr) != 1) {
	    goto corrupt;
	} else if (i == 3 && sscanf(line, "width height = %d %d", 
				    &grp->width, &grp->height) != 2) {
	    goto corrupt;
	} else if (i == 4 && sscanf(line, "boxwidth = %lf", 
				    &grp->boxwidth) != 1) {
	    goto corrupt;
	} else if (i == 5 && sscanf(line, "gmax gmin = %lf %lf", 
				    &grp->gmax, &grp->gmin) != 2) {
	    goto corrupt;
	}
    }

    if (strcmp(numstr, "NULL")) {
	grp->numbers = malloc(strlen(numstr) + 1);
	if (grp->numbers != NULL) {
	    strcpy(grp->numbers, numstr);
	}
    }

    grp->plots = malloc(grp->nplots * sizeof *grp->plots);
    if (grp->plots == NULL) {
	free(grp);
	fclose(fp);
	return 1;
    }

    for (i=0; i<grp->nplots; i++) {
	plt = &grp->plots[i];
	plt->outliers = NULL;
	plt->bool = NULL;
	nout = 0;

	for (j=0; j<7 && fgets(line, 79, fp); j++) {
	    if (j == 0 && 
		sscanf(line, "%*d median = %lf", &plt->median) != 1) {
		goto corrupt;
	    } else if (j == 1 && 
		     sscanf(line, "%*d conf = %lf %lf", 
			    &plt->conf[0], & plt->conf[1]) != 2) {
		goto corrupt;
	    } else if (j == 2 && 
		     sscanf(line, "%*d quartiles = %lf %lf", 
			    &plt->uq, &plt->lq) != 2) {
		goto corrupt;
	    } else if (j == 3 && 
		     sscanf(line, "%*d maxmin = %lf %lf", 
			    &plt->max, &plt->min) != 2) {
		goto corrupt;
	    } else if (j == 4 && 
		       sscanf(line, "%*d xbase = %lf", &plt->xbase) != 1) {
		goto corrupt;
	    } else if (j == 5) {
		if (sscanf(line, "%*d varname = %15s", plt->varname) != 1) {
		    goto corrupt;
		} else {
		    get_bool_from_line(line, plt);
		}
	    } else if (j == 6 && 
		       sscanf(line, "%*d n_outliers = %d", &nout) != 1) {
		goto corrupt;
	    }
	}

	screen_confidence_interval(plt);

	/* any outliers? */
	if (nout > 0) {
	    err = plot_retrieve_outliers(plt, nout, fp);
	    if (err == E_DATA) {
		goto corrupt;
	    }
	} 

	maybe_get_plot_mean(plt, fp);
    }
    
    fclose(fp);

    grp->popup = NULL;
    grp->pixmap = NULL;
    place_plots(grp);

    if (make_area(grp) == NULL) {
	free(grp->plots);
	free(grp);
	return 1;
    }

    gtk_widget_show(grp->window);

    for (i=0; i<grp->nplots; i++)
	gtk_area_boxplot (&grp->plots[i], 
			  grp->area, grp->pixmap,
			  grp->window->style, NULL, 
			  grp->height, grp->boxwidth, 
			  grp->gmax, grp->gmin,
			  grp->numbers);
    
    gtk_boxplot_yscale(grp, NULL);
    
    return 0;

 corrupt:

    errbox(_("boxplot file is corrupt"));
    fclose(fp);

    return 1;
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
    if (s == NULL) return NULL;

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

#define BPSTRLEN 128

static char boxplots_string[BPSTRLEN];

const char *get_boxplots_string (void)
{
    return boxplots_string;
}

int boolean_boxplots (const char *str, double ***pZ, DATAINFO *pdinfo, 
		      gretlopt opt)
{
    int i, k, v, nvars, nbool, err = 0;
    int n = pdinfo->n, origv = pdinfo->v;
    char *tok, *s = NULL, **bools = NULL;
    int *list = NULL;

    if (!strncmp(str, "boxplots ", 9)) str += 9;
    else if (!strncmp(str, "boxplot ", 8)) str += 8;

    s = boxplots_fix_parentheses(str);
    if (s == NULL) return 1;

    nvars = special_varcount(s);
    if (nvars == 0) {
	free(s);
	return 1;
    }

    list = malloc((nvars + 1) * sizeof *list);
    bools = malloc(nvars * sizeof *bools);
    if (list == NULL || bools == NULL) {
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
		if (v < origv) list[++i] = v;
		else {
		    errbox(_("got invalid variable number %d"), v);
		    err = 1;
		}
	    } else if (isalpha(tok[0])) {
		v = varindex(pdinfo, tok);
		if (v < origv) list[++i] = v;
		else {
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
		    if ((*pZ)[k][t] == 1.0) 
			(*pZ)[k][t] = (*pZ)[list[i]][t];
		    else 
			(*pZ)[k][t] = NADBL;
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


