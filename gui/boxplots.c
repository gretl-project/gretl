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

/* boxplots for gretl */

#include "gretl.h"

#ifndef G_OS_WIN32
# include <gtkextra/gtkextra.h>
#else
# include "gtkextra.h"
# include <windows.h>
#endif

typedef struct {
    double median;
    double conf[2];
    double uq, lq;
    double max, min;
    double xbase;
    char varname[9];
} BOXPLOT;

typedef struct {
    int nplots;
    BOXPLOT *plots;
    int width, height;
    double boxwidth;
    double gmax, gmin;
    GtkWidget *window, *area, *popup;
    GdkPixmap *pixmap;
} PLOTGROUP;

double headroom = 0.24;
double scalepos = 60.0;
char boxfont[64] = "Helvetica";
int boxfontsize = 12;

extern void file_selector (char *msg, char *startdir, int action, 
                           gpointer data);
extern char remember_dir[MAXLEN];

int ps_print_plots (const char *fname, int flag, gpointer data);
static int five_numbers (gpointer data);

#ifdef G_OS_WIN32
static int cb_copy_image (gpointer data);
#else
int plot_to_xpm (const char *fname, gpointer data);
#endif

/* ............................................................. */

/* Create a new backing pixmap of the appropriate size */
static gint
configure_event (GtkWidget *widget, GdkEventConfigure *event, gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    if (grp->pixmap)
	gdk_pixmap_unref(grp->pixmap);

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

/* ............................................................. */

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

/* ............................................................. */

static gboolean 
box_key_handler (GtkWidget *w, GdkEventKey *key, gpointer data)
{
    if (key->keyval == GDK_q) {
	gtk_widget_destroy(w);
    }
    else if (key->keyval == GDK_s) {
        file_selector("Save boxplot file", 
		      (remember_dir[0])? remember_dir : paths.userdir, 
                      SAVE_BOXPLOT_EPS, data);
    }
    else if (key->keyval == GDK_p) {  
	five_numbers(data);
    }
#ifdef G_OS_WIN32
    else if (key->keyval == GDK_c) {  
	cb_copy_image(data);
    }
#else
    else if (key->keyval == GDK_c) { 
        file_selector("Save boxplot file", 
		      (remember_dir[0])? remember_dir : paths.userdir, 
                      SAVE_BOXPLOT_XPM, data);	
    }
#endif
    return TRUE;
}

/* ........................................................... */

static gint popup_activated (GtkWidget *w, gpointer data)
{
    gchar *item = (gchar *) data;
    gpointer ptr = gtk_object_get_data(GTK_OBJECT(w), "group");
    PLOTGROUP *grp = (PLOTGROUP *) ptr;

    if (!strcmp(item, "Five-number summary")) 
        five_numbers(grp);
    else if (!strcmp(item, "Save as EPS...")) 
        file_selector("Save boxplot file", 
		      (remember_dir[0])? remember_dir : paths.userdir, 
                      SAVE_BOXPLOT_EPS, ptr);
    else if (!strcmp(item, "Save as PS...")) 
        file_selector("Save boxplot file", 
		      (remember_dir[0])? remember_dir : paths.userdir, 
                      SAVE_BOXPLOT_PS, ptr);
#ifdef G_OS_WIN32
    else if (!strcmp(item, "Copy to clipboard"))
	cb_copy_image(ptr);
#else
    else if (!strcmp(item, "Save as XPM..."))
        file_selector("Save boxplot file", 
		      (remember_dir[0])? remember_dir : paths.userdir, 
                      SAVE_BOXPLOT_XPM, ptr);
#endif
    else if (!strcmp(item, "Help"))
	do_help (NULL, GR_BOX, NULL);
    else if (!strcmp(item, "Close")) { 
	gtk_widget_destroy(grp->popup);
	grp->popup = NULL;
        gtk_widget_destroy(grp->window);
	return TRUE;
    }

    gtk_widget_destroy(grp->popup);
    grp->popup = NULL;
    return TRUE;
}

/* ........................................................... */

static GtkWidget *build_menu (gpointer data)
{
    GtkWidget *menu, *item;    
    static char *items[] = {
        "Five-number summary",
        "Save as EPS...",
        "Save as PS...",
#ifdef G_OS_WIN32
	"Copy to clipboard",
#else
	"Save as XPM...",
#endif
	"Help",
        "Close",
	NULL
    };
    int i = 0;

    menu = gtk_menu_new();

    while (items[i]) {
        item = gtk_menu_item_new_with_label(items[i]);
        gtk_signal_connect(GTK_OBJECT(item), "activate",
                           (GtkSignalFunc) popup_activated,
                           items[i]);
	gtk_object_set_data(GTK_OBJECT(item), "group", data);
        GTK_WIDGET_SET_FLAGS (item, GTK_SENSITIVE | GTK_CAN_FOCUS);
        gtk_widget_show(item);
        gtk_menu_append(GTK_MENU(menu), item);
	i++;
    }
    return menu;
}

/* ........................................................... */

static gint box_popup (GtkWidget *widget, GdkEventButton *event, 
		       gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    if (grp->popup) g_free(grp->popup);
    grp->popup = build_menu(data);
    gtk_menu_popup(GTK_MENU(grp->popup), NULL, NULL, NULL, NULL,
		   event->button, event->time);
    return TRUE;
}

/* ............................................................. */

static void
setup_text (GtkWidget *area, GdkPixmap *pixmap,
	    GdkGC *gc, GtkPlotPC *pc, char *text, 
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
	GdkRectangle rect; 
	size_t len = strlen(text);
	int cw = gdk_char_width(fixed_font, 'X');
	
	if (just == GTK_JUSTIFY_CENTER) {
	    x -= len * cw / 2.0;
	}
	else if (just == GTK_JUSTIFY_RIGHT) {
	    int ch = gdk_char_height(fixed_font, '1');

	    y += (double) ch / 2.0;
	    x -= len * cw;
	}
	rect.x = x;
	rect.y = y;
	rect.width = 80;
	rect.height = 10;	
	gdk_draw_string (pixmap, fixed_font, gc, x, y, text);
	gtk_widget_draw (area, &rect);
    }
}

/* ............................................................. */

static void 
draw_line (double *points, GtkWidget *area, GdkPixmap *pixmap,
	   GdkGC *gc, GtkPlotPC *pc)
{
    if (pc != NULL) {
	gtk_plot_pc_draw_line (pc, points[0], points[1], points[2], points[3]); 
    } else {
	GdkRectangle rect;

	rect.x = points[0];
	rect.y = points[1];
	rect.width = points[2] - points[1] + 1;
	rect.height = points[3] - points[1] + 1;
	gdk_draw_line (pixmap, gc, points[0], points[1], points[2], points[3]);
	gtk_widget_draw (area, &rect);
    }
}

/* ............................................................. */

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

/* ............................................................. */

static void 
gtk_boxplot_yscale (PLOTGROUP *grp, GtkPlotPC *pc)
{
    double points[4];
    double top = (headroom / 2.0) * grp->height;
    double bottom = (1.0 - headroom / 2.0) * grp->height;
    char numstr[16];
    GdkGC *gc = NULL;

    if (pc == NULL) gc = grp->window->style->fg_gc[GTK_STATE_NORMAL];
    /* gc = grp->window->style->bg_gc[GTK_STATE_SELECTED]; */

    /* draw vertical line */
    points[0] = points[2] = scalepos;
    points[1] = top;
    points[3] = bottom;
    draw_line (points, grp->area, grp->pixmap, gc, pc);

    /* draw backticks top and bottom */
    points[2] = points[0] - 5;
    points[1] = points[3] = top;
    draw_line (points, grp->area, grp->pixmap, gc, pc);
    points[1] = points[3] = bottom;
    draw_line (points, grp->area, grp->pixmap, gc, pc);

    /* draw backtick at middle */
    points[1] = points[3] = top + (bottom - top) / 2.0;
    draw_line (points, grp->area, grp->pixmap, gc, pc);
    
    /* mark max and min values on scale */
    sprintf(numstr, "%.4g", grp->gmax);
    setup_text (grp->area, grp->pixmap, gc, pc, numstr, scalepos - 10, top, 
		GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", grp->gmin);
    setup_text (grp->area, grp->pixmap, gc, pc, numstr, scalepos - 10, bottom, 
		GTK_JUSTIFY_RIGHT);
    sprintf(numstr, "%.4g", (grp->gmax + grp->gmin) / 2.0);
    setup_text (grp->area, grp->pixmap, gc, pc, numstr, scalepos - 10, 
		top + (bottom - top) / 2.0, GTK_JUSTIFY_RIGHT);

    /* special on-screen string for notched plots */
    if (pc == NULL && grp->plots[0].conf[0] != -999.0 && grp->width >=460) {
	setup_text (grp->area, grp->pixmap, gc, pc, 
		    "notches show bootstrapped 90% confidence intervals "
		    "for medians", 
		    grp->width / 2.0,
		    grp->height * headroom / 3.0,
		    GTK_JUSTIFY_CENTER);
    }
}

/* ............................................................. */

static void 
gtk_area_boxplot (BOXPLOT *plot, GtkWidget *area, GdkPixmap *pixmap,
		  GtkStyle *style, GtkPlotPC *pc,
		  int height, double boxwidth, 
		  double gmax, double gmin)
{
    double points[4];
    double ybase = height * headroom / 2.0;
    double xcenter = plot->xbase + boxwidth / 2.0;
    double scale = (1.0 - headroom) * height / (gmax - gmin);
    double median, uq, lq, maxval, minval;
    double conflo = 0., confhi = 0.;
    GdkRectangle rect;
    GdkGC *gc = NULL;
    /* GdkGC *whitegc = NULL; */

    if (pc == NULL) {
	gc = style->fg_gc[GTK_STATE_NORMAL];
	/* gc = style->bg_gc[GTK_STATE_SELECTED]; */
	/* whitegc = style->fg_gc[GTK_STATE_SELECTED]; */
    }

    median = ybase + (gmax - plot->median) * scale;
    uq = ybase + (gmax - plot->uq) * scale;
    lq = ybase + (gmax - plot->lq) * scale;
    maxval = ybase + (gmax - plot->max) * scale;
    minval = ybase + (gmax - plot->min) * scale;

    if (plot->conf[0] != -999.0) { /* confidence intervals defined */
	if (plot->conf[1] > plot->uq) 
	    confhi = uq;
	else
	    confhi = ybase + (gmax - plot->conf[1]) * scale;
	if (plot->conf[0] < plot->lq) 
	    conflo = lq;
	else
	    conflo = ybase + (gmax - plot->conf[0]) * scale;
    }

    /* no notches: draw simple inter-quartile box */
    if (confhi == 0.) {
	if (pc != NULL) 
	    gtk_plot_pc_draw_rectangle (pc,
					FALSE,
					plot->xbase, uq,
					boxwidth,
					lq - uq);
	else {
	    rect.x = plot->xbase;
	    rect.y = uq;
	    rect.width = boxwidth;
	    rect.height = uq - lq; 	
	    gdk_draw_rectangle (pixmap, 
				gc, 
				FALSE, /* filled ? */
				plot->xbase, 
				uq, 
				boxwidth, 
				lq - uq);
	    gtk_widget_draw (area, &rect);
	}
    } else { /* draw notched boxes */
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

	    rect.x = plot->xbase;
	    rect.y = uq;
	    rect.width = boxwidth;
	    rect.height = uq - lq; 

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
	    gtk_widget_draw (area, &rect);
	}
    }

    /* draw line at median */
    points[0] = plot->xbase + ((confhi > 0.)? (0.1 * boxwidth) : 0.0);
    points[1] = points[3] = median;
    points[2] = plot->xbase + ((confhi > 0.)? (0.9 * boxwidth) : boxwidth);
    draw_line(points, area, pixmap, gc, pc); /* was whitegc */

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

    /* write name of variable beneath */
    setup_text (area, pixmap, gc, pc, plot->varname, xcenter, 
		height * (1.0 - headroom/4.0), GTK_JUSTIFY_CENTER);
}

/* ............................................................. */

static void 
destroy_boxplots (GtkWidget *w, gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;

    free(grp->plots);
    gdk_pixmap_unref(grp->pixmap);
    free(grp);
}

/* ............................................................. */

static GtkWidget *
make_area (PLOTGROUP *grp)
{
    grp->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(grp->window), "gretl: boxplots");

    /* Create the drawing area */
    grp->area = gtk_drawing_area_new ();

    gtk_widget_set_events (grp->area, GDK_EXPOSURE_MASK
			   | GDK_LEAVE_NOTIFY_MASK
			   | GDK_BUTTON_PRESS_MASK);
    
    gtk_widget_set_sensitive(grp->area, TRUE);

    gtk_signal_connect(GTK_OBJECT(grp->area), "configure_event",
		       GTK_SIGNAL_FUNC(configure_event), grp);

    gtk_signal_connect(GTK_OBJECT(grp->area), "expose_event",
		       GTK_SIGNAL_FUNC(expose_event), grp);

    gtk_signal_connect(GTK_OBJECT(grp->area), "button_press_event", 
		       GTK_SIGNAL_FUNC(box_popup), grp);

    gtk_signal_connect(GTK_OBJECT(grp->window), "key_press_event", 
		       GTK_SIGNAL_FUNC(box_key_handler), grp);

    gtk_signal_connect(GTK_OBJECT(grp->window), "destroy",
		       GTK_SIGNAL_FUNC(destroy_boxplots), grp);

    gtk_drawing_area_size (GTK_DRAWING_AREA(grp->area), 
			   grp->width, grp->height); 
    gtk_widget_show (grp->area);

    gtk_container_add(GTK_CONTAINER(grp->window), grp->area);

    return grp->window;
}

/* ............................................................. */

static int 
compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;
     
    return (*da > *db) - (*da < *db);
}

/* ............................................................. */

static double 
median (double *x, const int n)
{
    int n2;
    double xx;

    qsort(x, n, sizeof *x, compare_doubles);

    n2 = n/2;
    xx = (n % 2)? x[n2] : 0.5 * (x[n2 - 1] + x[n2]);
    return xx;
}

/* ............................................................. */

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

/* ............................................................. */

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
	    t = rand() / (RAND_MAX / n + 1);
	    samp[j] = x[t];
	}
	/* find the median of the sample */
	medians[i] = median(samp, n);
    }

    /* sort the sample medians */
    qsort(medians, ITERS, sizeof *medians, compare_doubles);
    
    /* return the right values */
    j = 100 / ((100 - CONFIDENCE) / 2);
    *low = medians[ITERS / j];
    *high = medians[ITERS - ITERS/j];

    free(samp);
    free(medians);

    return 0;
}

/* At this point we need to borrow a few functions from gtkplotps.c,
   which is part of gtkextra.  These functions are declared as
   static in that context, but we want direct access to them */

#include "plotps.c"

/* ............................................................. */

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
				      pscale, pscale)
		     );

    if (ps == NULL) return 1;

    if (eps) {
	ps->page_width = pscale * grp->width;
	ps->page_height = pscale * grp->height;
#ifdef notdef
	ps->page_width = ps->pc.width = pscale * grp->width;
	ps->page_height = ps->pc.height = pscale * grp->height;
#endif
    }

    if (!psinit (ps)) return 1;
    gtk_psfont_init();

    for (i=0; i<grp->nplots; i++)
	gtk_area_boxplot (&grp->plots[i], 
			  NULL, NULL, NULL, &ps->pc, 
			  grp->height, grp->boxwidth, 
			  grp->gmax, grp->gmin);
    
    gtk_boxplot_yscale (grp, &ps->pc);

    psleave (ps);
    gtk_object_destroy(GTK_OBJECT(ps));
    gtk_psfont_unref();

    return 0;
}

/* ............................................................. */

static int
five_numbers (gpointer data) 
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    int i;
    print_t prn;
    extern GtkItemFactoryEntry view_items[];

    if (bufopen(&prn)) return 1;

    if (grp->plots[0].conf[0] == -999.0) { /* no confidence intervals */
	pprintf(&prn, "Five-number summar%s\n\n"
		"%20s%10s%10s%10s%10s\n",
		(grp->nplots > 1)? "ies" : "y",
		"min", "Q1", "median", "Q3", "max");

	for (i=0; i<grp->nplots; i++) {
	    pprintf(&prn, "%-10s%10g%10g%10g%10g%10g\n",
		    grp->plots[i].varname, grp->plots[i].min, 
		    grp->plots[i].lq, grp->plots[i].median,
		    grp->plots[i].uq, grp->plots[i].max);
	}
    } else { /* confidence intervals */
	char intstr[24];

	pprintf(&prn, "Five-number summar%s with bootstrapped confidence "
		"interval for median\n\n"
		"%18s%10s%10s%17s%10s%10s\n",
		(grp->nplots > 1)? "ies" : "y",
		"min", "Q1", "median", "(90% interval)", "Q3", "max");

	for (i=0; i<grp->nplots; i++) {
	    sprintf(intstr, "%g - %g", grp->plots[i].conf[0], 
		    grp->plots[i].conf[1]);
	    pprintf(&prn, "%-10s%8g%10g%10g%17s%10g%10g\n",
		    grp->plots[i].varname, grp->plots[i].min, 
		    grp->plots[i].lq, grp->plots[i].median,
		    intstr,
		    grp->plots[i].uq, grp->plots[i].max);
	}
    }

    (void) view_buffer(&prn, 78, 240, "gretl: 5 numbers", BXPLOT,
                       view_items);

    return 0;
}

static void read_boxrc (PLOTGROUP *grp);

/* ............................................................. */

int boxplots (int *list, double **pZ, const DATAINFO *pdinfo, int notches)
{
    int i, n = pdinfo->t2 - pdinfo->t1 + 1;
    double *x;
    PLOTGROUP *plotgrp;
    int width = 576, height = 448;

    x = mymalloc(n * sizeof *x);
    if (x == NULL) return 1;

    plotgrp = mymalloc(sizeof *plotgrp);
    if (plotgrp == NULL) {
	free(x);
	return 1;
    }

    plotgrp->nplots = list[0];
    plotgrp->plots = mymalloc (plotgrp->nplots * sizeof(BOXPLOT));
    if (plotgrp->plots == NULL) {
	free(plotgrp);
	free(x);
	return 1;
    }

    for (i=0; i<plotgrp->nplots; i++) {
	n = ztox(list[i+1], x, pdinfo, *pZ);
	if (n < 2) {
	    char s[63];

	    sprintf(s, "Dropping %s: inadequate observations for boxplot",
		    pdinfo->varname[list[i+1]]);
	    errbox(s);
	    list_exclude(i+1, list);
	    if (list[0] == 0) {
		free(plotgrp->plots);
		free(plotgrp);
		free(x);
		return 1;
	    } else {
		plotgrp->nplots -= 1;
		i--;
		continue;
	    }
	}
	qsort(x, n, sizeof *x, compare_doubles);
	plotgrp->plots[i].min = x[0];
	plotgrp->plots[i].max = x[n-1];
	quartiles(x, n, &plotgrp->plots[i]);
	/* notched boxplots wanted? */
	if (notches) {
	    if (median_interval(x, n, &plotgrp->plots[i].conf[0],
				&plotgrp->plots[i].conf[1])) {
		errbox ("Couldn't obtain confidence interval");
		plotgrp->plots[i].conf[0] = 
		    plotgrp->plots[i].conf[1] = -999.0;
	    }
	} else 
	    plotgrp->plots[i].conf[0] = plotgrp->plots[i].conf[1] = -999.0;
	strcpy(plotgrp->plots[i].varname, pdinfo->varname[list[i+1]]);
    }
    free(x);

    plotgrp->height = height;
    plotgrp->width = width;

    read_boxrc(plotgrp);

    plotgrp->popup = NULL;
    plotgrp->pixmap = NULL;
    place_plots (plotgrp);

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
			  height, plotgrp->boxwidth, 
			  plotgrp->gmax, plotgrp->gmin);
    
    gtk_boxplot_yscale(plotgrp, NULL);
    
    return 0;
}

/* copy functions */

#ifdef G_OS_WIN32

static int cb_copy_image (gpointer data)
{
    PLOTGROUP *grp = (PLOTGROUP *) data;
    GdkImage *image;
    int i, j;
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
	    errbox ("Failed to allocate DIB");
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
    } else 
	errbox("Failed to lock DIB Header");

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
	} else
	    errbox ("Failed to lock DIB Palette");
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
	} else 
	    errbox("Failed to lock DIB Data");
    } /* copy data to DIB */
  
    /* copy DIB to ClipBoard */
    if (ret) {      
	if (!OpenClipboard (NULL)) {
	    errbox ("Cannot open the Clipboard!");
	    ret = FALSE;
	} else {
	    if (ret && !EmptyClipboard ()) {
		errbox ("Cannot empty the Clipboard");
		ret = FALSE;
	    }
	    if (ret) {
		if (NULL != SetClipboardData (CF_DIB, hDIB))
		    hDIB = NULL; /* data now owned by clipboard */
		else
		    errbox ("Failed to set clipboard data");
	    }
	    if (!CloseClipboard ())
		errbox ("Failed to close Clipboard");
	}
    }
    /* done */
    if (hDIB) GlobalFree(hDIB);
  
    gdk_image_destroy (image);
  
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
    
    fp = fopen(fname, "w");
    if (fp == NULL) {
	errbox ("Couldn't open XPM file for writing");
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

    gdk_image_destroy(image);

    fclose(fp);

    return 0;
}

#endif 

static void read_boxrc (PLOTGROUP *grp)
{
    FILE *fp;

    grp->gmax = grp->gmin = -999.0;

    fp = fopen(".boxplotrc", "r");
    if (fp == NULL) {
	char *homedir, boxrc[MAXLEN];

	homedir = getenv("HOME");
	if (homedir == NULL) homedir = paths.userdir;
	sprintf(boxrc, "%s%s.boxplotrc", homedir,
		(homedir[strlen(homedir)-1] != SLASH)? SLASHSTR : "");
	fp = fopen(boxrc, "r");
    }

    if (fp != NULL) {
	char line[80], key[18], val[32];

	while (fgets(line, 79, fp)) {
	    if (line[0] == '#') continue;
	    if (sscanf(line, "%17s = %31s", key, val) == 2) {
		key[17] = '\0';
		val[31] = '\0';
		if (!strcmp(key, "max")) 
		    grp->gmax = atof(val);
		else if (!strcmp(key, "min")) 
		    grp->gmin = atof(val);
		else if (!strcmp(key, "font")) { 
		    strncpy(boxfont, val, 63);
		    boxfont[63] = '\0';
		}
		else if (!strcmp(key, "fontsize")) 
		    boxfontsize = atoi(val);
		else if (!strcmp(key, "width") && atoi(val) > 0) 
		    grp->width = atoi(val);
		else if (!strcmp(key, "height") && atoi(val) > 0) 
		    grp->height = atoi(val);
	    }
	}
	fclose (fp);
    }

    if (grp->gmax == -999.0) grp->gmax = grp->plots[0].max;
    if (grp->gmin == -999.0) grp->gmin = grp->plots[0].min;
}
